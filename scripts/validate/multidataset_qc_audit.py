#!/usr/bin/env python3
"""QC audit for multidataset RRBS region matrices."""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler

ROOT = Path("/home/zdq-as/mouse_methyl_work")
META_FILE = ROOT / "metadata" / "unified_sample_metadata.csv"
GSM_RE = re.compile(r"(GSM\d+)")
GSE60012_TILE_RE = re.compile(r"^(GSE60012_tile_\d{3})")


def resolve_matrix_sample_id(value: object) -> str:
    text = str(value)
    tile = GSE60012_TILE_RE.match(text)
    if tile:
        return tile.group(1)
    gsm = GSM_RE.search(text)
    if gsm:
        return gsm.group(1)
    return text


def eta_squared(values: np.ndarray, labels: pd.Series) -> float:
    valid = np.isfinite(values) & labels.notna().to_numpy()
    if valid.sum() < 3:
        return np.nan
    y = values[valid]
    lab = labels.to_numpy()[valid]
    overall = float(np.mean(y))
    ss_total = float(np.sum((y - overall) ** 2))
    if ss_total == 0:
        return 0.0
    ss_between = 0.0
    for group in pd.unique(lab):
        mask = lab == group
        ss_between += float(mask.sum() * (np.mean(y[mask]) - overall) ** 2)
    return ss_between / ss_total


def age_bins(age_weeks: pd.Series) -> pd.Series:
    return pd.cut(
        age_weeks,
        bins=[0, 4, 13, 26, 52, 104, np.inf],
        labels=["0-4w", "4-13w", "13-26w", "26-52w", "52-104w", "104w+"],
        include_lowest=True,
    ).astype(str)


def align_meta(matrix: pd.DataFrame, metadata_path: Path) -> pd.DataFrame:
    meta = pd.read_csv(metadata_path)
    sample_ids = [resolve_matrix_sample_id(col) for col in matrix.columns]
    meta = meta.set_index("sample_id").reindex(sample_ids).reset_index()
    meta["age_bin"] = age_bins(meta["age_weeks"])
    return meta


def write_qc_summary(matrix: pd.DataFrame, meta: pd.DataFrame, out_dir: Path) -> pd.DataFrame:
    X = matrix.T.to_numpy(dtype=np.float32)
    sample_qc = meta[["sample_id", "dataset_batch", "tissue", "age_weeks", "age_bin", "intervention"]].copy()
    sample_qc["sample_missing_fraction"] = np.isnan(X).mean(axis=1)
    sample_qc["sample_mean_beta"] = np.nanmean(X, axis=1)
    sample_qc["sample_sd_beta"] = np.nanstd(X, axis=1)
    sample_qc.to_csv(out_dir / "sample_qc_metrics.csv", index=False)

    summary = (
        sample_qc.groupby(["dataset_batch", "tissue", "age_bin"], dropna=False)
        .agg(
            n_samples=("sample_id", "count"),
            age_weeks_min=("age_weeks", "min"),
            age_weeks_median=("age_weeks", "median"),
            age_weeks_max=("age_weeks", "max"),
            missingness_mean=("sample_missing_fraction", "mean"),
            mean_beta=("sample_mean_beta", "mean"),
            sd_beta=("sample_sd_beta", "mean"),
        )
        .reset_index()
    )
    summary.to_csv(out_dir / "multidataset_qc_summary.csv", index=False)
    return sample_qc


def write_feature_shift(matrix: pd.DataFrame, meta: pd.DataFrame, out_dir: Path) -> None:
    X = matrix.T.to_numpy(dtype=np.float32)
    global_mean = np.nanmean(X, axis=0)
    rows = []
    for dataset in sorted(meta["dataset_batch"].dropna().unique()):
        mask = meta["dataset_batch"].eq(dataset).to_numpy()
        subset = X[mask]
        group_mean = np.nanmean(subset, axis=0)
        abs_shift = np.abs(group_mean - global_mean)
        top_idx = int(np.nanargmax(abs_shift))
        rows.append(
            {
                "dataset_batch": dataset,
                "n_samples": int(mask.sum()),
                "mean_sample_missing_fraction": float(np.isnan(subset).mean()),
                "mean_abs_region_shift": float(np.nanmean(abs_shift)),
                "median_abs_region_shift": float(np.nanmedian(abs_shift)),
                "top_shift_region": str(matrix.index[top_idx]),
                "top_shift_abs_beta": float(abs_shift[top_idx]),
                "mean_region_variance": float(np.nanmean(np.nanvar(subset, axis=0))),
            }
        )
    pd.DataFrame(rows).to_csv(out_dir / "dataset_feature_shift.csv", index=False)


def write_pca(matrix: pd.DataFrame, meta: pd.DataFrame, out_dir: Path) -> None:
    X = matrix.T.to_numpy(dtype=np.float32)
    X = SimpleImputer(strategy="median").fit_transform(X)
    X = StandardScaler().fit_transform(X)
    n_components = min(10, X.shape[0] - 1)
    pca = PCA(n_components=n_components, svd_solver="randomized", random_state=42)
    scores = pca.fit_transform(X)

    score_cols = [f"PC{i + 1}" for i in range(n_components)]
    score_df = meta[["sample_id", "dataset_batch", "tissue", "age_weeks", "age_bin", "intervention"]].copy()
    for idx, col in enumerate(score_cols):
        score_df[col] = scores[:, idx]
    score_df.to_csv(out_dir / "pca_batch_separation.csv", index=False)

    eta_rows = []
    for idx, col in enumerate(score_cols):
        eta_rows.append(
            {
                "component": col,
                "explained_variance_ratio": float(pca.explained_variance_ratio_[idx]),
                "dataset_eta_squared": eta_squared(scores[:, idx], meta["dataset_batch"]),
                "tissue_eta_squared": eta_squared(scores[:, idx], meta["tissue"]),
                "age_bin_eta_squared": eta_squared(scores[:, idx], meta["age_bin"]),
            }
        )
    pd.DataFrame(eta_rows).to_csv(out_dir / "pca_eta_squared.csv", index=False)


def write_residual_summary(predictions: Path, out_dir: Path) -> None:
    if not predictions.exists():
        return
    pred = pd.read_csv(predictions)
    if "age_bin" not in pred.columns:
        pred["age_bin"] = age_bins(pred["age_weeks_true"])
    if "abs_error_weeks" not in pred.columns:
        pred["abs_error_weeks"] = (pred["age_weeks_true"] - pred["age_weeks_pred"]).abs()
    summary = (
        pred.groupby(["dataset_batch", "tissue", "age_bin"], dropna=False)
        .agg(
            n_samples=("sample_id", "count"),
            mae_weeks=("abs_error_weeks", "mean"),
            median_abs_error_weeks=("abs_error_weeks", "median"),
            residual_weeks_mean=("residual_weeks", "mean"),
            residual_weeks_median=("residual_weeks", "median"),
        )
        .reset_index()
    )
    summary.to_csv(out_dir / "residual_by_dataset_tissue_agebin.csv", index=False)


def write_gse60012_mapping_audit(out_dir: Path, metadata_path: Path) -> None:
    parse_stats = ROOT / "results" / "multidataset" / "GSE60012_sample_parse_stats.csv"
    manifest_path = ROOT / "results" / "multidataset" / "GSE60012_matrix_manifest.json"
    v71_metadata = ROOT / "metadata" / "gse60012_header_sample_metadata.csv"
    v71_exclusions = ROOT / "results" / "multidataset_v7_1" / "GSE60012_header_mapping_exclusions.csv"
    rows = []
    if metadata_path.name == "model_sample_metadata_v7_1.csv" and v71_metadata.exists():
        parsed = pd.read_csv(v71_metadata)
        exclusions = pd.read_csv(v71_exclusions) if v71_exclusions.exists() else pd.DataFrame()
        rows.append(
            {
                "dataset": "GSE60012",
                "status": "mapped_header_synthetic_non_gsm",
                "reason": "official_tile_matrix_header_parsed_without_gsm_claim",
                "n_tile_columns": int(len(parsed) + len(exclusions)),
                "duplicate_source_columns": None,
                "metadata_overlap": int(len(parsed)),
                "excluded_columns": int(len(exclusions)),
            }
        )
    elif parse_stats.exists():
        stats = pd.read_csv(parse_stats)
        duplicate_source_columns = int(stats["source_column"].duplicated().sum()) if "source_column" in stats.columns else None
        rows.append(
            {
                "dataset": "GSE60012",
                "status": "blocked",
                "reason": "official_tile_matrix_columns_are_not_unique_gsm_accessions",
                "n_tile_columns": int(len(stats)),
                "duplicate_source_columns": duplicate_source_columns,
                "metadata_overlap": 0,
            }
        )
    elif manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        rows.append(
            {
                "dataset": "GSE60012",
                "status": manifest.get("status"),
                "reason": manifest.get("metadata_mapping_status", "unknown"),
                "n_tile_columns": manifest.get("raw_sample_columns"),
                "duplicate_source_columns": None,
                "metadata_overlap": manifest.get("n_samples_metadata_overlap"),
            }
        )
    pd.DataFrame(rows).to_csv(out_dir / "gse60012_mapping_audit.csv", index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--matrix",
        default=str(ROOT / "results" / "multidataset" / "all_rrbs_region_matrix_5kb.parquet"),
    )
    parser.add_argument("--metadata_path", default=str(META_FILE))
    parser.add_argument(
        "--predictions",
        default=str(ROOT / "results" / "benchmark_v6_groupkfold" / "lgbm_500_median" / "predictions.csv"),
    )
    parser.add_argument("--out_dir", default=str(ROOT / "results" / "benchmark_v7_harmonization" / "qc"))
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    matrix = pd.read_parquet(args.matrix)
    matrix.columns = [resolve_matrix_sample_id(col) for col in matrix.columns]
    meta = align_meta(matrix, Path(args.metadata_path))

    write_qc_summary(matrix, meta, out_dir)
    write_feature_shift(matrix, meta, out_dir)
    write_pca(matrix, meta, out_dir)
    write_residual_summary(Path(args.predictions), out_dir)
    write_gse60012_mapping_audit(out_dir, Path(args.metadata_path))
    print(f"[QC] Wrote audit outputs -> {out_dir}")


if __name__ == "__main__":
    main()

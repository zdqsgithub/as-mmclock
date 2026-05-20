#!/usr/bin/env python3
"""GSE121141-focused coverage and shift diagnostics for v7.5.

This script is diagnostic only. Its all-dataset region tables must not be used
as a benchmark feature mask, because that would expose held-out dataset coverage
to GroupKFold training.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("/home/zdq-as/mouse_methyl_work")
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


def age_bins(age_weeks: pd.Series) -> pd.Series:
    return pd.cut(
        age_weeks,
        bins=[0, 4, 13, 26, 52, 104, np.inf],
        labels=["0-4w", "4-13w", "13-26w", "26-52w", "52-104w", "104w+"],
        include_lowest=True,
    ).astype(str)


def safe_pearson(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    valid = np.isfinite(x) & np.isfinite(y)
    if valid.sum() < 3 or np.std(x[valid]) == 0 or np.std(y[valid]) == 0:
        return 0.0, 1.0
    r, pval = stats.pearsonr(x[valid], y[valid])
    if np.isnan(r):
        return 0.0, 1.0
    return float(r), float(pval)


def load_aligned(matrix_path: Path, metadata_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    matrix = pd.read_parquet(matrix_path)
    matrix.columns = [resolve_matrix_sample_id(col) for col in matrix.columns]
    meta = pd.read_csv(metadata_path)
    meta = meta[meta["age_days"].notna()].copy()
    common = [sample for sample in meta["sample_id"] if sample in matrix.columns]
    meta = meta[meta["sample_id"].isin(common)].reset_index(drop=True)
    meta["age_bin"] = age_bins(meta["age_weeks"])
    return matrix[common], meta


def write_dataset_sample_summary(matrix: pd.DataFrame, meta: pd.DataFrame, out_dir: Path) -> pd.DataFrame:
    X = matrix.T.to_numpy(dtype=np.float32)
    sample_qc = meta[
        ["sample_id", "dataset_batch", "tissue", "sex", "age_weeks", "age_bin", "intervention"]
    ].copy()
    sample_qc["sample_missing_fraction"] = np.isnan(X).mean(axis=1)
    sample_qc["sample_mean_beta"] = np.nanmean(X, axis=1)
    sample_qc["sample_sd_beta"] = np.nanstd(X, axis=1)
    sample_qc.to_csv(out_dir / "sample_level_coverage_qc.csv", index=False)

    summary = (
        sample_qc.groupby(["dataset_batch", "tissue", "age_bin"], dropna=False)
        .agg(
            n_samples=("sample_id", "count"),
            age_min=("age_weeks", "min"),
            age_median=("age_weeks", "median"),
            age_max=("age_weeks", "max"),
            missing_mean=("sample_missing_fraction", "mean"),
            missing_median=("sample_missing_fraction", "median"),
            mean_beta=("sample_mean_beta", "mean"),
            sd_beta=("sample_sd_beta", "mean"),
        )
        .reset_index()
    )
    summary.to_csv(out_dir / "dataset_tissue_agebin_sample_summary.csv", index=False)
    return sample_qc


def write_region_dataset_diagnostics(matrix: pd.DataFrame, meta: pd.DataFrame, out_dir: Path) -> pd.DataFrame:
    X = matrix.T.to_numpy(dtype=np.float32)
    datasets = sorted(meta["dataset_batch"].dropna().unique())
    rows = pd.DataFrame({"region_id": matrix.index.astype(str)})
    mean_cols = []
    presence_cols = []
    for dataset in datasets:
        mask = meta["dataset_batch"].eq(dataset).to_numpy()
        subset = X[mask]
        presence_col = f"{dataset}_presence"
        mean_col = f"{dataset}_mean_beta"
        rows[presence_col] = np.isfinite(subset).mean(axis=0)
        with np.errstate(invalid="ignore"):
            rows[mean_col] = np.nanmean(subset, axis=0)
        presence_cols.append(presence_col)
        mean_cols.append(mean_col)

    mean_values = rows[mean_cols].to_numpy(dtype=float)
    presence_values = rows[presence_cols].to_numpy(dtype=float)
    rows["all_dataset_min_presence"] = np.nanmin(presence_values, axis=1)
    rows["all_dataset_mean_presence"] = np.nanmean(presence_values, axis=1)
    rows["all_dataset_max_mean_shift"] = np.nanmax(mean_values, axis=1) - np.nanmin(mean_values, axis=1)

    if "GSE121141_mean_beta" in rows.columns:
        other_mean_cols = [col for col in mean_cols if col != "GSE121141_mean_beta"]
        other_means = rows[other_mean_cols].to_numpy(dtype=float)
        with np.errstate(invalid="ignore"):
            rows["non_gse121141_mean_beta"] = np.nanmean(other_means, axis=1)
        rows["gse121141_abs_shift_vs_others"] = (
            rows["GSE121141_mean_beta"] - rows["non_gse121141_mean_beta"]
        ).abs()
    else:
        rows["non_gse121141_mean_beta"] = np.nan
        rows["gse121141_abs_shift_vs_others"] = np.nan

    rows["stable_all_data_p095_shift015"] = (
        (rows["all_dataset_min_presence"] >= 0.95) & (rows["all_dataset_max_mean_shift"] <= 0.15)
    )
    rows["stable_all_data_p08_shift015"] = (
        (rows["all_dataset_min_presence"] >= 0.8) & (rows["all_dataset_max_mean_shift"] <= 0.15)
    )
    rows.to_csv(out_dir / "region_dataset_presence_shift.csv", index=False)

    summary = {
        "n_regions": int(len(rows)),
        "stable_all_data_p095_shift015": int(rows["stable_all_data_p095_shift015"].sum()),
        "stable_all_data_p08_shift015": int(rows["stable_all_data_p08_shift015"].sum()),
        "gse121141_presence_mean": float(rows.get("GSE121141_presence", pd.Series(dtype=float)).mean()),
        "gse121141_presence_median": float(rows.get("GSE121141_presence", pd.Series(dtype=float)).median()),
        "gse121141_shift_gt_015": int((rows["gse121141_abs_shift_vs_others"] > 0.15).sum()),
        "gse121141_shift_gt_010": int((rows["gse121141_abs_shift_vs_others"] > 0.10).sum()),
        "diagnostic_only": True,
    }
    (out_dir / "region_dataset_presence_shift_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    return rows


def write_pairwise_region_overlap(datasets: list[str], out_dir: Path) -> pd.DataFrame:
    region_sets: dict[str, set[str]] = {}
    for dataset in datasets:
        if dataset == "GSE120137":
            path = ROOT / "results" / "phase0" / "region_matrix_5kb.parquet"
        else:
            path = ROOT / "results" / "multidataset" / f"{dataset}_region_matrix_5kb.parquet"
        if not path.exists():
            continue
        matrix = pd.read_parquet(path, columns=[])
        region_sets[dataset] = set(matrix.index.astype(str))

    rows = []
    for left in sorted(region_sets):
        for right in sorted(region_sets):
            if left >= right:
                continue
            inter = region_sets[left] & region_sets[right]
            union = region_sets[left] | region_sets[right]
            rows.append(
                {
                    "dataset_left": left,
                    "dataset_right": right,
                    "left_regions": len(region_sets[left]),
                    "right_regions": len(region_sets[right]),
                    "intersection_regions": len(inter),
                    "union_regions": len(union),
                    "jaccard": round(len(inter) / len(union), 6) if union else 0.0,
                }
            )
    overlap = pd.DataFrame(rows)
    overlap.to_csv(out_dir / "dataset_pair_region_overlap.csv", index=False)
    return overlap


def write_prediction_failure_summary(predictions: Path, sample_qc: pd.DataFrame, out_dir: Path) -> None:
    if not predictions.exists():
        return
    pred = pd.read_csv(predictions)
    pred["abs_error_weeks"] = (pred["age_weeks_true"] - pred["age_weeks_pred"]).abs()
    pred["age_bin"] = age_bins(pred["age_weeks_true"])
    pred = pred.merge(
        sample_qc[["sample_id", "sample_missing_fraction", "sample_mean_beta", "sample_sd_beta"]],
        on="sample_id",
        how="left",
    )
    pred.to_csv(out_dir / "predictions_with_sample_qc.csv", index=False)
    summary = (
        pred.groupby(["dataset_batch", "tissue", "age_bin"], dropna=False)
        .agg(
            n_samples=("sample_id", "count"),
            mae_weeks=("abs_error_weeks", "mean"),
            median_abs_error_weeks=("abs_error_weeks", "median"),
            residual_mean=("residual_weeks", "mean"),
            sample_missing_mean=("sample_missing_fraction", "mean"),
        )
        .reset_index()
    )
    summary.to_csv(out_dir / "prediction_error_by_dataset_tissue_agebin.csv", index=False)

    rows = []
    for dataset, sub in pred.groupby("dataset_batch"):
        r_missing, p_missing = safe_pearson(
            sub["sample_missing_fraction"].to_numpy(float), sub["abs_error_weeks"].to_numpy(float)
        )
        r_age, p_age = safe_pearson(sub["age_weeks_true"].to_numpy(float), sub["abs_error_weeks"].to_numpy(float))
        rows.append(
            {
                "dataset_batch": dataset,
                "n_samples": int(len(sub)),
                "missing_vs_abs_error_r": round(r_missing, 4),
                "missing_vs_abs_error_p": p_missing,
                "age_vs_abs_error_r": round(r_age, 4),
                "age_vs_abs_error_p": p_age,
                "mae_weeks": round(float(sub["abs_error_weeks"].mean()), 3),
            }
        )
    pd.DataFrame(rows).to_csv(out_dir / "sample_qc_error_correlations.csv", index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--matrix",
        default=str(ROOT / "results" / "multidataset" / "all_rrbs_region_matrix_5kb.parquet"),
    )
    parser.add_argument(
        "--metadata_path",
        default=str(ROOT / "metadata" / "model_sample_metadata_v7_1.csv"),
    )
    parser.add_argument(
        "--predictions",
        default=str(
            ROOT
            / "results"
            / "benchmark_v7_4_stability"
            / "05_lgbm_robust_p095_top1000_shift015"
            / "predictions.csv"
        ),
    )
    parser.add_argument(
        "--datasets",
        default="GSE120137,GSE80672,GSE93957,GSE121141",
        help="Comma-separated datasets for pairwise raw region overlap diagnostics.",
    )
    parser.add_argument(
        "--out_dir",
        default=str(ROOT / "results" / "benchmark_v7_5_gse121141_harmonization" / "qc"),
    )
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    matrix, meta = load_aligned(Path(args.matrix), Path(args.metadata_path))
    sample_qc = write_dataset_sample_summary(matrix, meta, out_dir)
    region_diag = write_region_dataset_diagnostics(matrix, meta, out_dir)
    overlap = write_pairwise_region_overlap([item.strip() for item in args.datasets.split(",") if item.strip()], out_dir)
    write_prediction_failure_summary(Path(args.predictions), sample_qc, out_dir)

    summary = {
        "matrix": str(Path(args.matrix)),
        "metadata_path": str(Path(args.metadata_path)),
        "predictions": str(Path(args.predictions)),
        "n_samples": int(matrix.shape[1]),
        "n_regions": int(matrix.shape[0]),
        "datasets": sorted(meta["dataset_batch"].dropna().unique().tolist()),
        "gse121141_samples": int(meta["dataset_batch"].eq("GSE121141").sum()),
        "stable_all_data_p095_shift015": int(region_diag["stable_all_data_p095_shift015"].sum()),
        "stable_all_data_p08_shift015": int(region_diag["stable_all_data_p08_shift015"].sum()),
        "pairwise_overlap_rows": int(len(overlap)),
        "diagnostic_only_no_training_mask": True,
    }
    (out_dir / "v7_5_gse121141_audit_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"[v7.5 audit] Wrote outputs -> {out_dir}")


if __name__ == "__main__":
    main()


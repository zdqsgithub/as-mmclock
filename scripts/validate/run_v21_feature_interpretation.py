#!/usr/bin/env python3
"""Generate v21 top-region interpretation and confounding audit tables."""
from __future__ import annotations

import argparse
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats


LOCAL_ROOT = Path("/home/zdq-as/mouse_methyl_work")
REMOTE_ROOT = Path("/root/autodl-tmp/mouse_methyl_work")
ROOT = REMOTE_ROOT if REMOTE_ROOT.exists() else LOCAL_ROOT
DEFAULT_MATRIX = ROOT / "results" / "multidataset_v8_3_ablation" / "all6" / "all_rrbs_region_matrix_5kb.parquet"
DEFAULT_METADATA = ROOT / "metadata" / "model_sample_metadata_v8.csv"
DEFAULT_ML_FEATURES = ROOT / "results" / "autoresearch_v19_ml_baseline" / "final_all_data_model" / "selected_features.csv"
DEFAULT_DL_FEATURES = ROOT / "results" / "autoresearch_v20_deep_learning" / "final_all_data_model" / "selected_features.csv"
DEFAULT_OUT = ROOT / "results" / "ralph_v21_raid_raw_clock" / "feature_interpretation"
GSM_RE = re.compile(r"(GSM\d+)")
GSE60012_TILE_RE = re.compile(r"^(GSE60012_tile_\d{3})")


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str), encoding="utf-8")


def resolve_matrix_sample_id(value: object) -> str:
    text = str(value)
    tile = GSE60012_TILE_RE.match(text)
    if tile:
        return tile.group(1)
    gsm = GSM_RE.search(text)
    if gsm:
        return gsm.group(1)
    return text


def parse_region(region_id: str) -> dict[str, Any]:
    match = re.match(r"^(chr[^:]+):(\d+)-(\d+)$", str(region_id))
    if not match:
        return {"chrom": "", "start": None, "end": None, "cluster_5mb": ""}
    chrom = match.group(1)
    start = int(match.group(2))
    end = int(match.group(3))
    cluster_start = (start // 5_000_000) * 5_000_000
    return {"chrom": chrom, "start": start, "end": end, "cluster_5mb": f"{chrom}:{cluster_start}-{cluster_start + 4_999_999}"}


def safe_corr(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    valid = np.isfinite(x) & np.isfinite(y)
    if valid.sum() < 5 or np.nanstd(x[valid]) == 0 or np.nanstd(y[valid]) == 0:
        return 0.0, 1.0
    r, p = stats.pearsonr(x[valid], y[valid])
    if not np.isfinite(r):
        return 0.0, 1.0
    return float(r), float(p)


def load_selected(path: Path, source: str) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=["feature_id", f"{source}_rank", f"{source}_selected"])
    df = pd.read_csv(path)
    if "feature_id" not in df.columns:
        raise SystemExit(f"selected feature file lacks feature_id: {path}")
    df = df[["feature_id"]].drop_duplicates().reset_index(drop=True)
    df[f"{source}_rank"] = np.arange(1, len(df) + 1)
    df[f"{source}_selected"] = True
    return df


def build_feature_universe(ml_features: Path, dl_features: Path, top_n: int) -> pd.DataFrame:
    ml = load_selected(ml_features, "ml")
    dl = load_selected(dl_features, "dl")
    if ml.empty and dl.empty:
        raise SystemExit("No selected feature files found for interpretation.")
    merged = ml.merge(dl, on="feature_id", how="outer")
    merged["ml_selected"] = merged.get("ml_selected", False).fillna(False).astype(bool)
    merged["dl_selected"] = merged.get("dl_selected", False).fillna(False).astype(bool)
    merged["ml_rank"] = pd.to_numeric(merged.get("ml_rank"), errors="coerce")
    merged["dl_rank"] = pd.to_numeric(merged.get("dl_rank"), errors="coerce")
    merged["min_rank"] = merged[["ml_rank", "dl_rank"]].min(axis=1, skipna=True)
    merged = merged.sort_values(["min_rank", "feature_id"]).head(top_n).reset_index(drop=True)
    merged["interpretation_rank"] = np.arange(1, len(merged) + 1)
    return merged


def classify_region(
    *,
    dataset_support_count: int,
    tissue_support_count: int,
    n_datasets: int,
    n_tissues: int,
    missing_fraction: float,
    dataset_mean_shift: float,
    tissue_mean_shift: float,
) -> str:
    if missing_fraction > 0.25:
        return "coverage_driven"
    if dataset_mean_shift > 0.25 and dataset_support_count < max(2, n_datasets // 2):
        return "dataset_confounded"
    if tissue_mean_shift > 0.25 and tissue_support_count < max(2, n_tissues // 2):
        return "tissue_specific"
    if dataset_support_count >= max(2, n_datasets // 2) and tissue_support_count >= max(2, n_tissues // 2):
        return "cross_tissue_stable"
    return "mixed_or_low_support"


def interpret_features(matrix: pd.DataFrame, meta: pd.DataFrame, features: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    age_weeks = meta["age_days"].astype(float).to_numpy() / 7.0
    datasets = meta["dataset_batch"].fillna("unknown").astype(str).to_numpy()
    tissues = meta.get("tissue", pd.Series(["unknown"] * len(meta))).fillna("unknown").astype(str).to_numpy()
    unique_datasets = sorted(set(datasets))
    unique_tissues = sorted(set(tissues))

    for _, feature in features.iterrows():
        fid = str(feature["feature_id"])
        if fid not in matrix.index:
            continue
        x = matrix.loc[fid].astype(float).to_numpy()
        r, p = safe_corr(x, age_weeks)
        sign = 1 if r >= 0 else -1
        dataset_corrs = {}
        tissue_corrs = {}
        for dataset in unique_datasets:
            mask = datasets == dataset
            dataset_corrs[dataset] = safe_corr(x[mask], age_weeks[mask])[0] if mask.sum() >= 5 else 0.0
        for tissue in unique_tissues:
            mask = tissues == tissue
            tissue_corrs[tissue] = safe_corr(x[mask], age_weeks[mask])[0] if mask.sum() >= 5 else 0.0
        dataset_support = sum(1 for value in dataset_corrs.values() if np.sign(value) == sign and abs(value) >= 0.1)
        tissue_support = sum(1 for value in tissue_corrs.values() if np.sign(value) == sign and abs(value) >= 0.1)
        missing_fraction = float(np.mean(~np.isfinite(x)))
        dataset_means = [float(np.nanmean(x[datasets == dataset])) for dataset in unique_datasets if np.isfinite(x[datasets == dataset]).any()]
        tissue_means = [float(np.nanmean(x[tissues == tissue])) for tissue in unique_tissues if np.isfinite(x[tissues == tissue]).any()]
        dataset_shift = float(np.nanmax(dataset_means) - np.nanmin(dataset_means)) if dataset_means else 0.0
        tissue_shift = float(np.nanmax(tissue_means) - np.nanmin(tissue_means)) if tissue_means else 0.0
        source_count = int(bool(feature.get("ml_selected"))) + int(bool(feature.get("dl_selected")))
        stability = (dataset_support / max(1, len(unique_datasets)) + tissue_support / max(1, len(unique_tissues)) + source_count / 2) / 3
        ml_rank = feature.get("ml_rank")
        dl_rank = feature.get("dl_rank")
        ml_importance_proxy = 0.0 if pd.isna(ml_rank) else 1.0 / float(ml_rank)
        dl_importance_proxy = 0.0 if pd.isna(dl_rank) else 1.0 / float(dl_rank)
        row = {
            **parse_region(fid),
            "region_id": fid,
            "interpretation_rank": int(feature["interpretation_rank"]),
            "selected_by_ml": bool(feature.get("ml_selected")),
            "selected_by_dl": bool(feature.get("dl_selected")),
            "ml_rank": None if pd.isna(ml_rank) else int(ml_rank),
            "dl_rank": None if pd.isna(dl_rank) else int(dl_rank),
            "age_pearson_r": round(r, 6),
            "age_pearson_p": p,
            "age_direction": "hypermethylating_with_age" if r >= 0 else "hypomethylating_with_age",
            "dataset_support_count": int(dataset_support),
            "tissue_support_count": int(tissue_support),
            "n_datasets_evaluated": int(len(unique_datasets)),
            "n_tissues_evaluated": int(len(unique_tissues)),
            "missing_fraction": round(missing_fraction, 6),
            "dataset_mean_shift": round(dataset_shift, 6),
            "tissue_mean_shift": round(tissue_shift, 6),
            "fold_stability_proxy": round(float(stability), 6),
            "permutation_importance_proxy": round(float(abs(r) * (1 - missing_fraction)), 8),
            "ml_rank_importance_proxy": round(float(ml_importance_proxy), 8),
            "dl_gradient_importance_proxy": round(float(dl_importance_proxy * abs(r)), 8),
            "dl_occlusion_importance_proxy": round(float(dl_importance_proxy * np.nanstd(x)), 8),
        }
        row["confounding_label"] = classify_region(
            dataset_support_count=dataset_support,
            tissue_support_count=tissue_support,
            n_datasets=len(unique_datasets),
            n_tissues=len(unique_tissues),
            missing_fraction=missing_fraction,
            dataset_mean_shift=dataset_shift,
            tissue_mean_shift=tissue_shift,
        )
        row["interpretation_note"] = "age-associated methylation signal; mechanism not inferred"
        rows.append(row)
    return pd.DataFrame(rows)


def write_report(out_dir: Path, top: pd.DataFrame, cluster: pd.DataFrame, audit: pd.DataFrame, summary: dict[str, Any]) -> None:
    lines = [
        "# v21 Interpretable Clock Report",
        "",
        f"Date: {utc_now()}",
        "",
        "## Scope",
        "",
        "Region-level interpretation of selected v19 ML and v20 DL clock features. The labels describe age-associated methylation signals and confounding risk; they do not establish aging mechanisms.",
        "",
        "## Summary",
        "",
        f"- top regions interpreted: `{len(top)}`",
        f"- cross-tissue stable: `{int((top['confounding_label'] == 'cross_tissue_stable').sum()) if not top.empty else 0}`",
        f"- tissue-specific: `{int((top['confounding_label'] == 'tissue_specific').sum()) if not top.empty else 0}`",
        f"- dataset-confounded: `{int((top['confounding_label'] == 'dataset_confounded').sum()) if not top.empty else 0}`",
        f"- coverage-driven: `{int((top['confounding_label'] == 'coverage_driven').sum()) if not top.empty else 0}`",
        "",
        "## Important Caveat",
        "",
        "Permutation, DL gradient, and DL occlusion columns are deterministic interpretation proxies unless a downstream run supplies model-native importance tensors. They are adequate for triage and reporting guardrails, not mechanistic claims.",
        "",
        "## Outputs",
        "",
        f"- `{out_dir / 'top_region_interpretation.csv'}`",
        f"- `{out_dir / 'region_cluster_summary.csv'}`",
        f"- `{out_dir / 'feature_confounding_audit.csv'}`",
        f"- `{out_dir / 'v21_feature_interpretation_summary.json'}`",
    ]
    (out_dir / "v21_interpretable_clock_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    doc_report = ROOT / "doc" / "20_analysis" / "54_20260521_v21_interpretable_clock_report.md"
    doc_report.parent.mkdir(parents=True, exist_ok=True)
    doc_report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    write_json(out_dir / "v21_feature_interpretation_summary.json", summary)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--ml-selected-features", type=Path, default=DEFAULT_ML_FEATURES)
    parser.add_argument("--dl-selected-features", type=Path, default=DEFAULT_DL_FEATURES)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--top-n", type=int, default=100)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    features = build_feature_universe(args.ml_selected_features, args.dl_selected_features, args.top_n)
    matrix = pd.read_parquet(args.matrix)
    matrix.columns = [resolve_matrix_sample_id(col) for col in matrix.columns]
    metadata = pd.read_csv(args.metadata)
    metadata = metadata[metadata["age_days"].notna()].copy()
    common = [sample for sample in metadata["sample_id"].astype(str) if sample in matrix.columns]
    metadata = metadata[metadata["sample_id"].astype(str).isin(common)].reset_index(drop=True)
    matrix = matrix[common]
    top = interpret_features(matrix, metadata, features)
    top.to_csv(args.out_dir / "top_region_interpretation.csv", index=False)

    if top.empty:
        cluster = pd.DataFrame()
        audit = pd.DataFrame()
    else:
        cluster = (
            top.groupby(["chrom", "cluster_5mb", "confounding_label"], dropna=False)
            .agg(
                n_regions=("region_id", "count"),
                mean_abs_age_r=("age_pearson_r", lambda s: float(np.mean(np.abs(s)))),
                mean_stability=("fold_stability_proxy", "mean"),
                mean_missing_fraction=("missing_fraction", "mean"),
            )
            .reset_index()
            .sort_values(["n_regions", "mean_abs_age_r"], ascending=[False, False])
        )
        audit = top[
            top["confounding_label"].isin(["dataset_confounded", "coverage_driven", "tissue_specific"])
        ].copy()
    cluster.to_csv(args.out_dir / "region_cluster_summary.csv", index=False)
    audit.to_csv(args.out_dir / "feature_confounding_audit.csv", index=False)
    summary = {
        "timestamp": utc_now(),
        "matrix": str(args.matrix),
        "metadata": str(args.metadata),
        "n_samples": int(len(metadata)),
        "n_regions_interpreted": int(len(top)),
        "n_clusters": int(len(cluster)),
        "n_confounding_audit_rows": int(len(audit)),
        "labels": top["confounding_label"].value_counts().to_dict() if not top.empty else {},
    }
    write_report(args.out_dir, top, cluster, audit, summary)
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

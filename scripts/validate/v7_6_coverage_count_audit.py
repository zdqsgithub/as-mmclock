#!/usr/bin/env python3
"""Audit coverage/count structure for v7.6 coverage-weighted matrices."""
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


def resolve_matrix_sample_id(value: object) -> str:
    text = str(value)
    match = GSM_RE.search(text)
    if match:
        return match.group(1)
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


def load_sample_aligned_matrix(path: Path, meta: pd.DataFrame) -> pd.DataFrame:
    matrix = pd.read_parquet(path)
    matrix.columns = [resolve_matrix_sample_id(col) for col in matrix.columns]
    common = [sample for sample in meta["sample_id"] if sample in matrix.columns]
    return matrix[common]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--coverage_matrix",
        default=str(
            ROOT
            / "results"
            / "multidataset_v7_6_weighted"
            / "GSE121141_coverage_weighted_region_total_coverage_5kb.parquet"
        ),
    )
    parser.add_argument(
        "--cpg_count_matrix",
        default=str(
            ROOT
            / "results"
            / "multidataset_v7_6_weighted"
            / "GSE121141_coverage_weighted_region_cpg_counts_5kb.parquet"
        ),
    )
    parser.add_argument("--metadata_path", default=str(ROOT / "metadata" / "model_sample_metadata_v7_1.csv"))
    parser.add_argument("--dataset", default="GSE121141")
    parser.add_argument(
        "--predictions",
        default=str(
            ROOT
            / "results"
            / "validation_v7_5_gse121141_all_except"
            / "05_lgbm_quantile_p08_top1000_agebin08"
            / "predictions.csv"
        ),
    )
    parser.add_argument("--out_dir", default=str(ROOT / "results" / "benchmark_v7_6_coverage_weighted" / "coverage_qc"))
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(args.metadata_path)
    meta = meta[(meta["dataset_batch"] == args.dataset) & meta["age_days"].notna()].copy().reset_index(drop=True)
    meta["age_bin"] = age_bins(meta["age_weeks"])

    coverage = load_sample_aligned_matrix(Path(args.coverage_matrix), meta)
    cpg_counts = load_sample_aligned_matrix(Path(args.cpg_count_matrix), meta)
    meta = meta[meta["sample_id"].isin(coverage.columns)].copy().reset_index(drop=True)

    coverage_sample = coverage.T
    cpg_sample = cpg_counts.T
    sample_qc = meta[
        ["sample_id", "dataset_batch", "tissue", "sex", "age_weeks", "age_bin", "intervention"]
    ].copy()
    sample_qc["mean_region_total_coverage"] = coverage_sample.mean(axis=1, skipna=True).to_numpy(float)
    sample_qc["median_region_total_coverage"] = coverage_sample.median(axis=1, skipna=True).to_numpy(float)
    sample_qc["mean_region_cpg_count"] = cpg_sample.mean(axis=1, skipna=True).to_numpy(float)
    sample_qc["median_region_cpg_count"] = cpg_sample.median(axis=1, skipna=True).to_numpy(float)
    sample_qc["region_presence_fraction"] = coverage_sample.notna().mean(axis=1).to_numpy(float)
    sample_qc.to_csv(out_dir / f"{args.dataset}_sample_coverage_count_qc.csv", index=False)

    grouped = (
        sample_qc.groupby(["tissue", "age_bin"], dropna=False)
        .agg(
            n_samples=("sample_id", "count"),
            age_min=("age_weeks", "min"),
            age_median=("age_weeks", "median"),
            age_max=("age_weeks", "max"),
            mean_total_coverage=("mean_region_total_coverage", "mean"),
            median_total_coverage=("median_region_total_coverage", "mean"),
            mean_cpg_count=("mean_region_cpg_count", "mean"),
            region_presence=("region_presence_fraction", "mean"),
        )
        .reset_index()
    )
    grouped.to_csv(out_dir / f"{args.dataset}_coverage_by_tissue_agebin.csv", index=False)

    old_mask = sample_qc["age_weeks"] >= 104
    non_old_mask = sample_qc["age_weeks"] < 104
    region_shift = pd.DataFrame({"region_id": coverage.index.astype(str)})
    region_shift["old_mean_coverage"] = coverage.loc[:, old_mask.to_numpy()].mean(axis=1, skipna=True).to_numpy(float)
    region_shift["non_old_mean_coverage"] = coverage.loc[:, non_old_mask.to_numpy()].mean(axis=1, skipna=True).to_numpy(float)
    region_shift["old_minus_non_old_coverage"] = (
        region_shift["old_mean_coverage"] - region_shift["non_old_mean_coverage"]
    )
    region_shift["old_abs_log2_coverage_ratio"] = np.abs(
        np.log2((region_shift["old_mean_coverage"] + 1.0) / (region_shift["non_old_mean_coverage"] + 1.0))
    )
    region_shift.to_csv(out_dir / f"{args.dataset}_old_age_region_coverage_shift.csv", index=False)

    error_corr_rows = []
    pred_path = Path(args.predictions)
    if pred_path.exists():
        pred = pd.read_csv(pred_path)
        pred["abs_error_weeks"] = (pred["age_weeks_true"] - pred["age_weeks_pred"]).abs()
        merged = sample_qc.merge(pred[["sample_id", "age_weeks_true", "age_weeks_pred", "abs_error_weeks", "residual_weeks"]], on="sample_id", how="inner")
        merged.to_csv(out_dir / f"{args.dataset}_coverage_count_with_predictions.csv", index=False)
        for col in [
            "mean_region_total_coverage",
            "median_region_total_coverage",
            "mean_region_cpg_count",
            "region_presence_fraction",
            "age_weeks",
        ]:
            r, pval = safe_pearson(merged[col].to_numpy(float), merged["abs_error_weeks"].to_numpy(float))
            error_corr_rows.append(
                {
                    "variable": col,
                    "abs_error_pearson_r": round(r, 4),
                    "abs_error_pearson_p": pval,
                }
            )
    pd.DataFrame(error_corr_rows).to_csv(out_dir / f"{args.dataset}_coverage_error_correlations.csv", index=False)

    summary = {
        "dataset": args.dataset,
        "n_samples": int(len(sample_qc)),
        "n_regions": int(coverage.shape[0]),
        "old_age_samples_ge_104w": int(old_mask.sum()),
        "mean_region_total_coverage": float(sample_qc["mean_region_total_coverage"].mean()),
        "mean_region_cpg_count": float(sample_qc["mean_region_cpg_count"].mean()),
        "regions_abs_log2_old_nonold_coverage_ratio_gt_1": int(
            (region_shift["old_abs_log2_coverage_ratio"] > 1.0).sum()
        ),
    }
    (out_dir / f"{args.dataset}_coverage_count_audit_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(f"[v7.6 coverage audit] Wrote outputs -> {out_dir}")


if __name__ == "__main__":
    main()


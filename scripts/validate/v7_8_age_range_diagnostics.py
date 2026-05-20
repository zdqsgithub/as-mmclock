#!/usr/bin/env python3
"""Age-range support diagnostics and post-hoc calibration probes for v7.8.

This script does not train clocks. Calibration outputs are diagnostic-only and
must not be reported as valid held-out benchmark performance.
"""
from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import ConstantInputWarning
from sklearn.linear_model import LinearRegression

ROOT = Path("/home/zdq-as/mouse_methyl_work")

DEFAULT_RUNS = [
    {
        "run_id": "v7_4_groupkfold_best",
        "predictions": ROOT
        / "results"
        / "benchmark_v7_4_stability"
        / "05_lgbm_robust_p095_top1000_shift015"
        / "predictions.csv",
        "mode": "groupkfold_dataset",
        "result": ROOT
        / "results"
        / "benchmark_v7_4_stability"
        / "05_lgbm_robust_p095_top1000_shift015"
        / "result.json",
    },
    {
        "run_id": "v7_5_gse121141_heldout",
        "predictions": ROOT
        / "results"
        / "validation_v7_5_gse121141_all_except"
        / "05_lgbm_quantile_p08_top1000_agebin08"
        / "predictions.csv",
        "mode": "heldout",
        "result": ROOT
        / "results"
        / "validation_v7_5_gse121141_all_except"
        / "05_lgbm_quantile_p08_top1000_agebin08"
        / "result.json",
    },
    {
        "run_id": "v7_5_gse80672_cr_heldout",
        "predictions": ROOT
        / "results"
        / "validation_v7_5_gse80672_cr_all_except"
        / "05_lgbm_quantile_p08_top1000_agebin08"
        / "predictions.csv",
        "mode": "heldout",
        "result": ROOT
        / "results"
        / "validation_v7_5_gse80672_cr_all_except"
        / "05_lgbm_quantile_p08_top1000_agebin08"
        / "result.json",
    },
    {
        "run_id": "v7_6_gse121141_weighted_heldout",
        "predictions": ROOT
        / "results"
        / "validation_v7_6_gse121141_all_except"
        / "02_weighted_v75_quantile_p08_top1000_agebin08"
        / "predictions.csv",
        "mode": "heldout",
        "result": ROOT
        / "results"
        / "validation_v7_6_gse121141_all_except"
        / "02_weighted_v75_quantile_p08_top1000_agebin08"
        / "result.json",
    },
]


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
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ConstantInputWarning)
        r, pval = stats.pearsonr(x[valid], y[valid])
    if np.isnan(r):
        return 0.0, 1.0
    return float(r), float(pval)


def regression_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict:
    r, pval = safe_pearson(y_true, y_pred)
    mae = float(np.mean(np.abs(y_true - y_pred)))
    medae = float(np.median(np.abs(y_true - y_pred)))
    rmse = float(np.sqrt(np.mean((y_true - y_pred) ** 2)))
    ss_res = float(np.sum((y_true - y_pred) ** 2))
    ss_tot = float(np.sum((y_true - y_true.mean()) ** 2))
    r2 = 1 - ss_res / ss_tot if ss_tot > 1e-12 else 0.0
    return {
        "pearson_r": round(r, 4),
        "pearson_pval": pval,
        "mae_weeks": round(mae, 3),
        "medae_weeks": round(medae, 3),
        "rmse_weeks": round(rmse, 3),
        "r2": round(float(r2), 4),
    }


def load_json(path: Path) -> dict:
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    return {}


def selector_mask(meta: pd.DataFrame, selector: str) -> pd.Series:
    selector = str(selector or "all").strip()
    if selector.startswith("all_except:"):
        excluded = {item.strip() for item in selector.split(":", 1)[1].split(",") if item.strip()}
        return ~meta["dataset_batch"].isin(excluded)
    if selector in {"all", "*"}:
        return pd.Series(True, index=meta.index)
    allowed = {item.strip() for item in selector.split(",") if item.strip()}
    return meta["dataset_batch"].isin(allowed)


def train_meta_for_prediction_row(
    row: pd.Series,
    meta: pd.DataFrame,
    mode: str,
    result: dict,
) -> pd.DataFrame:
    if mode == "groupkfold_dataset":
        mask = meta["dataset_batch"].ne(row["dataset_batch"])
        return meta[mask & meta["age_weeks"].notna()].copy()
    selector = result.get("train_dataset") or result.get("train_datasets") or "all"
    return meta[selector_mask(meta, selector) & meta["age_weeks"].notna()].copy()


def add_age_support_flags(pred: pd.DataFrame, meta: pd.DataFrame, mode: str, result: dict) -> pd.DataFrame:
    records = []
    for _, row in pred.iterrows():
        train = train_meta_for_prediction_row(row, meta, mode, result)
        same_tissue = train[train["tissue"].eq(row.get("tissue"))]
        same_dataset_family = train[train["dataset_batch"].eq(row.get("dataset_batch"))]
        true_age = float(row["age_weeks_true"])
        global_min = float(train["age_weeks"].min()) if not train.empty else np.nan
        global_max = float(train["age_weeks"].max()) if not train.empty else np.nan
        tissue_min = float(same_tissue["age_weeks"].min()) if not same_tissue.empty else np.nan
        tissue_max = float(same_tissue["age_weeks"].max()) if not same_tissue.empty else np.nan
        records.append(
            {
                "train_n": int(len(train)),
                "train_global_age_min": global_min,
                "train_global_age_max": global_max,
                "train_same_tissue_n": int(len(same_tissue)),
                "train_same_tissue_age_min": tissue_min,
                "train_same_tissue_age_max": tissue_max,
                "test_dataset_seen_in_train": bool(len(same_dataset_family) > 0),
                "outside_global_train_age_range": bool(
                    np.isfinite(global_min) and np.isfinite(global_max) and (true_age < global_min or true_age > global_max)
                ),
                "above_global_train_max": bool(np.isfinite(global_max) and true_age > global_max),
                "below_global_train_min": bool(np.isfinite(global_min) and true_age < global_min),
                "same_tissue_missing_in_train": bool(same_tissue.empty),
                "outside_same_tissue_train_age_range": bool(
                    same_tissue.empty
                    or (np.isfinite(tissue_min) and np.isfinite(tissue_max) and (true_age < tissue_min or true_age > tissue_max))
                ),
                "above_same_tissue_train_max": bool(same_tissue.empty or (np.isfinite(tissue_max) and true_age > tissue_max)),
                "same_tissue_age_margin_weeks": (
                    np.nan if same_tissue.empty or not np.isfinite(tissue_max) else round(float(tissue_max - true_age), 3)
                ),
            }
        )
    flags = pd.DataFrame(records)
    return pd.concat([pred.reset_index(drop=True), flags], axis=1)


def prediction_compression_summary(df: pd.DataFrame, run_id: str) -> pd.DataFrame:
    rows = []
    for dataset, sub in df.groupby("dataset_batch", dropna=False):
        if len(sub) < 3:
            continue
        true = sub["age_weeks_true"].to_numpy(float)
        pred = sub["age_weeks_pred"].to_numpy(float)
        lr = LinearRegression().fit(true.reshape(-1, 1), pred)
        metrics = regression_metrics(true, pred)
        rows.append(
            {
                "run_id": run_id,
                "dataset_batch": dataset,
                "n_samples": int(len(sub)),
                "true_age_min": round(float(true.min()), 3),
                "true_age_max": round(float(true.max()), 3),
                "pred_age_min": round(float(pred.min()), 3),
                "pred_age_max": round(float(pred.max()), 3),
                "pred_vs_true_slope": round(float(lr.coef_[0]), 4),
                "pred_vs_true_intercept": round(float(lr.intercept_), 3),
                "predicted_age_range_weeks": round(float(pred.max() - pred.min()), 3),
                "true_age_range_weeks": round(float(true.max() - true.min()), 3),
                **metrics,
            }
        )
    return pd.DataFrame(rows)


def age_support_summary(df: pd.DataFrame, run_id: str) -> pd.DataFrame:
    rows = []
    group_cols = ["dataset_batch", "tissue", "age_bin"]
    for keys, sub in df.groupby(group_cols, dropna=False):
        true = sub["age_weeks_true"].to_numpy(float)
        pred = sub["age_weeks_pred"].to_numpy(float)
        rows.append(
            {
                "run_id": run_id,
                "dataset_batch": keys[0],
                "tissue": keys[1],
                "age_bin": keys[2],
                "n_samples": int(len(sub)),
                "true_age_min": round(float(true.min()), 3),
                "true_age_median": round(float(np.median(true)), 3),
                "true_age_max": round(float(true.max()), 3),
                "pred_age_median": round(float(np.median(pred)), 3),
                "mae_weeks": round(float(np.mean(np.abs(true - pred))), 3),
                "mean_residual_weeks": round(float(np.mean(true - pred)), 3),
                "pct_above_global_train_max": round(float(sub["above_global_train_max"].mean()), 4),
                "pct_above_same_tissue_train_max": round(float(sub["above_same_tissue_train_max"].mean()), 4),
                "same_tissue_train_age_max": (
                    None if sub["train_same_tissue_age_max"].isna().all() else round(float(sub["train_same_tissue_age_max"].max()), 3)
                ),
                "global_train_age_max": (
                    None if sub["train_global_age_max"].isna().all() else round(float(sub["train_global_age_max"].max()), 3)
                ),
            }
        )
    return pd.DataFrame(rows)


def oracle_linear_calibration(df: pd.DataFrame, run_id: str, dataset: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    sub = df[df["dataset_batch"].eq(dataset)].copy()
    if len(sub) < 5:
        return pd.DataFrame(), pd.DataFrame()
    true = sub["age_weeks_true"].to_numpy(float)
    pred = sub["age_weeks_pred"].to_numpy(float)
    model = LinearRegression().fit(pred.reshape(-1, 1), true)
    calibrated = model.predict(pred.reshape(-1, 1))
    overall = pd.DataFrame(
        [
            {
                "run_id": run_id,
                "dataset_batch": dataset,
                "calibration_type": "oracle_in_sample_true_age_from_pred_age",
                "invalid_as_benchmark": True,
                "calibration_intercept": round(float(model.intercept_), 4),
                "calibration_slope": round(float(model.coef_[0]), 4),
                **regression_metrics(true, calibrated),
            }
        ]
    )

    sub["age_months"] = (sub["age_weeks_true"].astype(float) / 4.34524).round(1)
    rows = []
    for age_months, test in sub.groupby("age_months"):
        train = sub[~sub["age_months"].eq(age_months)]
        if len(train) < 5 or len(test) < 1:
            continue
        loo_model = LinearRegression().fit(
            train["age_weeks_pred"].to_numpy(float).reshape(-1, 1),
            train["age_weeks_true"].to_numpy(float),
        )
        y_true = test["age_weeks_true"].to_numpy(float)
        y_cal = loo_model.predict(test["age_weeks_pred"].to_numpy(float).reshape(-1, 1))
        rows.append(
            {
                "run_id": run_id,
                "dataset_batch": dataset,
                "heldout_age_months": age_months,
                "n_test": int(len(test)),
                "invalid_as_benchmark": True,
                "calibration_intercept": round(float(loo_model.intercept_), 4),
                "calibration_slope": round(float(loo_model.coef_[0]), 4),
                **regression_metrics(y_true, y_cal),
            }
        )
    return overall, pd.DataFrame(rows)


def warnings_for_run(df: pd.DataFrame, compression: pd.DataFrame, run_id: str) -> list[dict]:
    warnings = []
    for _, row in compression.iterrows():
        if row["pred_vs_true_slope"] < 0.6 and row["true_age_range_weeks"] >= 50:
            warnings.append(
                {
                    "run_id": run_id,
                    "dataset_batch": row["dataset_batch"],
                    "warning": "prediction_age_range_compression",
                    "detail": f"slope={row['pred_vs_true_slope']}, true_range={row['true_age_range_weeks']}, pred_range={row['predicted_age_range_weeks']}",
                }
            )
    grouped = df.groupby(["dataset_batch", "age_bin"], dropna=False)
    for (dataset, age_bin), sub in grouped:
        mae = float(np.mean(np.abs(sub["age_weeks_true"] - sub["age_weeks_pred"])))
        if age_bin == "104w+" and mae >= 50:
            warnings.append(
                {
                    "run_id": run_id,
                    "dataset_batch": dataset,
                    "warning": "high_old_age_error",
                    "detail": f"age_bin={age_bin}, mae_weeks={mae:.3f}, n={len(sub)}",
                }
            )
        if sub["above_same_tissue_train_max"].mean() >= 0.5:
            warnings.append(
                {
                    "run_id": run_id,
                    "dataset_batch": dataset,
                    "warning": "test_age_exceeds_same_tissue_train_support",
                    "detail": f"age_bin={age_bin}, pct={sub['above_same_tissue_train_max'].mean():.3f}, n={len(sub)}",
                }
            )
    return warnings


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata_path", default=str(ROOT / "metadata" / "model_sample_metadata_v7_1.csv"))
    parser.add_argument("--out_dir", default=str(ROOT / "results" / "benchmark_v7_8_age_range_diagnostics"))
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(args.metadata_path)
    meta = meta[meta["age_days"].notna()].copy()

    all_predictions = []
    all_age_support = []
    all_compression = []
    all_calibration = []
    all_calibration_age_loo = []
    all_warnings = []
    manifest_runs = []

    for run in DEFAULT_RUNS:
        pred_path = Path(run["predictions"])
        if not pred_path.exists():
            continue
        result = load_json(Path(run["result"]))
        pred = pd.read_csv(pred_path)
        pred["run_id"] = run["run_id"]
        pred["age_bin"] = age_bins(pred["age_weeks_true"])
        pred["abs_error_weeks"] = (pred["age_weeks_true"] - pred["age_weeks_pred"]).abs()
        pred = add_age_support_flags(pred, meta, run["mode"], result)
        all_predictions.append(pred)
        support = age_support_summary(pred, run["run_id"])
        compression = prediction_compression_summary(pred, run["run_id"])
        all_age_support.append(support)
        all_compression.append(compression)
        cal, cal_loo = oracle_linear_calibration(pred, run["run_id"], "GSE121141")
        if not cal.empty:
            all_calibration.append(cal)
        if not cal_loo.empty:
            all_calibration_age_loo.append(cal_loo)
        all_warnings.extend(warnings_for_run(pred, compression, run["run_id"]))
        manifest_runs.append(
            {
                "run_id": run["run_id"],
                "predictions": str(pred_path),
                "mode": run["mode"],
                "result": str(run["result"]),
                "n_predictions": int(len(pred)),
            }
        )

    predictions_out = pd.concat(all_predictions, ignore_index=True)
    predictions_out.to_csv(out_dir / "predictions_with_age_range_flags.csv", index=False)
    age_support_out = pd.concat(all_age_support, ignore_index=True)
    age_support_out.to_csv(out_dir / "age_range_support_by_dataset_tissue_agebin.csv", index=False)
    compression_out = pd.concat(all_compression, ignore_index=True)
    compression_out.to_csv(out_dir / "prediction_age_compression_summary.csv", index=False)
    if all_calibration:
        pd.concat(all_calibration, ignore_index=True).to_csv(
            out_dir / "posthoc_oracle_linear_calibration_metrics.csv", index=False
        )
    if all_calibration_age_loo:
        pd.concat(all_calibration_age_loo, ignore_index=True).to_csv(
            out_dir / "posthoc_leave_age_level_out_calibration_metrics.csv", index=False
        )
    warnings_df = pd.DataFrame(all_warnings)
    warnings_df.to_csv(out_dir / "age_range_warning_flags.csv", index=False)

    summary = {
        "status": "completed",
        "metadata_path": str(Path(args.metadata_path)),
        "runs": manifest_runs,
        "n_prediction_rows": int(len(predictions_out)),
        "n_warning_rows": int(len(warnings_df)),
        "diagnostic_only_calibration_outputs": True,
        "warning_counts": warnings_df["warning"].value_counts().to_dict() if not warnings_df.empty else {},
    }
    (out_dir / "v7_8_age_range_diagnostic_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(f"[v7.8 age-range diagnostics] Wrote outputs -> {out_dir}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Benchmark metrics for mouse methylation age clocks.

The evaluator accepts either the v2 standard prediction schema
(`age_days_true`, `age_days_pred`) or legacy week/day columns and normalizes
units before computing metrics.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.metrics import f1_score, roc_auc_score

TARGETS = {
    "pearson_r": (">", 0.90),
    "mae_weeks": ("<", 3.5),
    "medae_weeks": ("<", 3.0),
    "r2": (">", 0.80),
    "cr_detection_auc": (">", 0.80),
    "cr_detection_f1": (">", 0.75),
    "cr_cohens_d": (">", 0.80),
    "cross_dataset_mae": ("<", 5.0),
}


def as_jsonable(value):
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    return value


def safe_pearson(y_true: np.ndarray, y_pred: np.ndarray) -> tuple[float, float]:
    if len(y_true) < 3 or np.std(y_true) == 0 or np.std(y_pred) == 0:
        return 0.0, 1.0
    r, pval = stats.pearsonr(y_true, y_pred)
    if np.isnan(r):
        return 0.0, 1.0
    return float(r), float(pval)


def compute_regression_metrics(y_true_days: np.ndarray, y_pred_days: np.ndarray) -> dict:
    y_true_wk = y_true_days / 7
    y_pred_wk = y_pred_days / 7
    r, pval = safe_pearson(y_true_wk, y_pred_wk)
    mae = float(np.mean(np.abs(y_true_wk - y_pred_wk)))
    medae = float(np.median(np.abs(y_true_wk - y_pred_wk)))
    rmse = float(np.sqrt(np.mean((y_true_wk - y_pred_wk) ** 2)))
    ss_res = float(np.sum((y_true_wk - y_pred_wk) ** 2))
    ss_tot = float(np.sum((y_true_wk - y_true_wk.mean()) ** 2))
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return {
        "pearson_r": round(r, 4),
        "pearson_pval": pval,
        "mae_weeks": round(mae, 3),
        "medae_weeks": round(medae, 3),
        "rmse_weeks": round(rmse, 3),
        "r2": round(float(r2), 4),
        "mae_pct_lifespan": round(mae / 100.0, 4),
    }


def compute_age_acceleration(
    y_pred_days: np.ndarray, y_true_days: np.ndarray, batch_labels: np.ndarray
) -> np.ndarray:
    """Return biological-age acceleration residuals.

    Positive values mean predicted methylation age is older than expected for
    chronological age within the same dataset batch. Negative values mean
    younger-than-expected predicted methylation age, which is the direction
    expected for CR validation.
    """
    acceleration = np.full(len(y_true_days), np.nan, dtype=float)
    for batch in np.unique(batch_labels):
        mask = batch_labels == batch
        if mask.sum() < 3 or np.std(y_true_days[mask]) == 0:
            acceleration[mask] = y_pred_days[mask] - y_true_days[mask]
            continue
        lr = LinearRegression()
        lr.fit(y_true_days[mask].reshape(-1, 1), y_pred_days[mask])
        acceleration[mask] = y_pred_days[mask] - lr.predict(y_true_days[mask].reshape(-1, 1))
    return acceleration


def compute_intervention_metrics(
    acceleration: np.ndarray, intervention_labels: np.ndarray, intervention: str
) -> dict:
    key = "cr" if intervention == "CR" else intervention.lower()
    is_intervention = (intervention_labels == intervention).astype(int)
    is_control = (intervention_labels == "control").astype(int)
    mask = (is_intervention == 1) | (is_control == 1)
    if mask.sum() < 4 or is_intervention[mask].sum() < 2 or is_control[mask].sum() < 2:
        return {
            f"{key}_detection_auc": None,
            f"{key}_detection_f1": None,
            f"{key}_cohens_d": None,
            f"{key}_mannwhitney_p": None,
        }

    accel_sub = acceleration[mask]
    label_sub = is_intervention[mask]
    valid = np.isfinite(accel_sub)
    accel_sub = accel_sub[valid]
    label_sub = label_sub[valid]
    if len(np.unique(label_sub)) < 2:
        return {
            f"{key}_detection_auc": None,
            f"{key}_detection_f1": None,
            f"{key}_cohens_d": None,
            f"{key}_mannwhitney_p": None,
        }

    score = -accel_sub
    auc = float(roc_auc_score(label_sub, score))
    pred_binary = (accel_sub < 0).astype(int)
    f1 = float(f1_score(label_sub, pred_binary, zero_division=0))

    intervention_accel = accel_sub[label_sub == 1]
    control_accel = accel_sub[label_sub == 0]
    pooled_std = np.sqrt((intervention_accel.std(ddof=1) ** 2 + control_accel.std(ddof=1) ** 2) / 2)
    if not np.isfinite(pooled_std) or pooled_std == 0:
        cohens_d = None
    else:
        cohens_d = float((control_accel.mean() - intervention_accel.mean()) / pooled_std)

    try:
        _, mw_p = stats.mannwhitneyu(control_accel, intervention_accel, alternative="greater")
        mw_p_value = float(mw_p)
    except Exception:
        mw_p_value = None

    return {
        f"{key}_detection_auc": round(auc, 4),
        f"{key}_detection_f1": round(f1, 4),
        f"{key}_cohens_d": round(cohens_d, 4) if cohens_d is not None else None,
        f"{key}_mannwhitney_p": round(mw_p_value, 6) if mw_p_value is not None else None,
    }


def normalize_prediction_columns(pred_df: pd.DataFrame) -> pd.DataFrame:
    df = pred_df.copy()
    if {"age_days_true", "age_days_pred"}.issubset(df.columns):
        pass
    elif {"age_weeks_true", "age_weeks_pred"}.issubset(df.columns):
        df["age_days_true"] = df["age_weeks_true"].astype(float) * 7
        df["age_days_pred"] = df["age_weeks_pred"].astype(float) * 7
    elif {"age_days", "age_weeks_pred"}.issubset(df.columns):
        df["age_days_true"] = df["age_days"].astype(float)
        df["age_days_pred"] = df["age_weeks_pred"].astype(float) * 7
    elif {"age_days", "age_days_pred"}.issubset(df.columns):
        df = df.rename(columns={"age_days": "age_days_true"})
    else:
        raise ValueError(
            "Predictions must contain age_days_true/age_days_pred, "
            "age_weeks_true/age_weeks_pred, or legacy age_days/age_weeks_pred columns."
        )

    if "dataset_batch" not in df.columns:
        df["dataset_batch"] = "unknown"
    if "intervention" not in df.columns:
        df["intervention"] = "control"
    if "tissue" not in df.columns:
        df["tissue"] = "unknown"
    return df


def evaluate_predictions(pred_df: pd.DataFrame) -> dict:
    df = normalize_prediction_columns(pred_df)
    y_true = df["age_days_true"].values.astype(float)
    y_pred = df["age_days_pred"].values.astype(float)
    batches = df["dataset_batch"].fillna("unknown").astype(str).values
    interventions = df["intervention"].fillna("control").astype(str).values

    metrics = compute_regression_metrics(y_true, y_pred)
    metrics["n_samples"] = int(len(df))
    metrics["n_datasets"] = int(len(np.unique(batches)))

    acceleration = compute_age_acceleration(y_pred, y_true, batches)
    for intervention in ["CR", "rapamycin", "dwarfism", "castration", "diet_high_fat"]:
        if intervention in set(interventions):
            metrics.update(compute_intervention_metrics(acceleration, interventions, intervention))

    cross_maes = []
    for batch in np.unique(batches):
        mask = batches == batch
        if mask.sum() >= 3:
            cross_maes.append(float(np.mean(np.abs(y_true[mask] - y_pred[mask])) / 7))
    metrics["cross_dataset_mae"] = round(float(np.mean(cross_maes)), 3) if cross_maes else None
    return metrics


def check_targets(metrics: dict) -> dict:
    checks = {}
    for key, (op, threshold) in TARGETS.items():
        val = metrics.get(key)
        if val is None:
            checks[key] = "N/A"
        elif op == ">" and val > threshold:
            checks[key] = f"PASS {val} > {threshold}"
        elif op == "<" and val < threshold:
            checks[key] = f"PASS {val} < {threshold}"
        else:
            checks[key] = f"FAIL {val} target {op} {threshold}"
    return checks


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--predictions", help="predictions CSV with required columns")
    parser.add_argument("--result", help="existing result.json to validate")
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    if args.predictions:
        pred_df = pd.read_csv(args.predictions)
        metrics = evaluate_predictions(pred_df)
        checks = check_targets(metrics)
        payload = {**metrics, "threshold_checks": checks}

        print("\n=== Benchmark Metrics ===")
        for key, value in metrics.items():
            print(f"  {key}: {value}")
        print("\n=== Publication Threshold Checks ===")
        for key, value in checks.items():
            print(f"  {key}: {value}")

        out_path = Path(args.out) if args.out else Path(args.predictions).parent / "benchmark_result.json"
        out_path.write_text(json.dumps(payload, indent=2, default=as_jsonable), encoding="utf-8")
        print(f"Saved -> {out_path}")
        return

    if args.result:
        data = json.loads(Path(args.result).read_text(encoding="utf-8"))
        checks = check_targets(data)
        print("=== Threshold Checks ===")
        for key, value in checks.items():
            print(f"  {key}: {value}")
        return

    parser.print_help()


if __name__ == "__main__":
    main()

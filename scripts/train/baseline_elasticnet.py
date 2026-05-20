#!/usr/bin/env python3
"""
Phase 0 Baseline: ElasticNet clock on Thompson 2018 beta matrix.
Input:  results/phase0/beta_matrix_thompson.parquet
        metadata/unified_sample_metadata.csv
Output: results/phase0/result_elasticnet.json
        results/phase0/predictions_elasticnet.csv
"""
import json, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import ElasticNetCV
from sklearn.model_selection import GroupKFold, KFold
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline

ROOT = Path("/home/zdq-as/mouse_methyl_work")
PHASE0 = ROOT / "results" / "phase0"
OUT = PHASE0

def compute_metrics(y_true_weeks, y_pred_weeks):
    r, pval = stats.pearsonr(y_true_weeks, y_pred_weeks)
    mae = np.mean(np.abs(y_true_weeks - y_pred_weeks))
    medae = np.median(np.abs(y_true_weeks - y_pred_weeks))
    ss_res = np.sum((y_true_weeks - y_pred_weeks) ** 2)
    ss_tot = np.sum((y_true_weeks - y_true_weeks.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot
    return {"pearson_r": round(float(r), 4), "pearson_pval": float(pval),
            "mae_weeks": round(float(mae), 3), "medae_weeks": round(float(medae), 3),
            "r2": round(float(r2), 4)}

def main():
    print("=" * 60)
    print("Phase 0: ElasticNet Baseline Clock")
    print("=" * 60)

    beta_path = PHASE0 / "beta_matrix_thompson.parquet"
    meta_path = ROOT / "metadata" / "unified_sample_metadata.csv"

    if not beta_path.exists():
        print(f"[ERROR] Run 00_load_thompson_matrix.py first.")
        sys.exit(1)

    print("[Baseline] Loading beta matrix... (this takes a moment)")
    beta = pd.read_parquet(beta_path)   # CpGs × samples
    meta = pd.read_csv(meta_path)

    # Filter samples with known age
    meta = meta[meta["age_days"].notna()].copy()
    
    # Strip file extensions and extra text from beta matrix columns to match GSM IDs
    # e.g., GSM3394233_Mouse_Blood_SH009 -> GSM3394233
    beta.columns = [c.split("_")[0] for c in beta.columns]
    
    common = [s for s in meta["sample_id"] if s in beta.columns]
    meta = meta[meta["sample_id"].isin(common)].reset_index(drop=True)
    
    if len(common) == 0:
        print("[ERROR] No common samples between beta matrix and metadata.")
        sys.exit(1)
        
    X = beta[common].T.values.astype(np.float32)  # samples × CpGs
    y_days = meta["age_days"].values.astype(np.float32)
    y_log  = np.log1p(y_days)  # log(days+1) target
    groups = meta["dataset_batch"].values

    print(f"[Baseline] Samples: {X.shape[0]}, CpGs: {X.shape[1]}")
    print(f"[Baseline] Age range: {y_days.min():.0f}–{y_days.max():.0f} days")

    # Variance pre-filter: keep top 20k most variable CpGs
    var = np.nanvar(X, axis=0)
    top_idx = np.argsort(var)[-20_000:]
    X = X[:, top_idx]
    print(f"[Baseline] After variance filter: {X.shape[1]} CpGs")

    # CV — GroupKFold by dataset_batch if multiple, else standard KFold
    n_groups = len(np.unique(groups))
    if n_groups >= 3:
        cv = GroupKFold(n_splits=min(5, n_groups))
        cv_iter = list(cv.split(X, y_log, groups))
    else:
        cv = KFold(n_splits=5, shuffle=True, random_state=42)
        cv_iter = list(cv.split(X, y_log))

    # Pipeline: impute → scale → ElasticNet
    pipe = Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("scaler",  StandardScaler()),
        ("model",   ElasticNetCV(
            l1_ratio=[0.1, 0.5, 0.7, 0.9, 0.95, 0.99, 1.0],
            cv=3, max_iter=10000, n_jobs=-1, random_state=42
        ))
    ])

    # Collect OOF predictions
    y_pred_log = np.full_like(y_log, np.nan)
    for fold, (tr, te) in enumerate(cv_iter):
        pipe.fit(X[tr], y_log[tr])
        y_pred_log[te] = pipe.predict(X[te])
        print(f"  Fold {fold+1}: test n={len(te)}")

    # Convert back to weeks
    y_true_wk = y_days / 7
    y_pred_wk = np.expm1(y_pred_log) / 7

    metrics = compute_metrics(y_true_wk, y_pred_wk)
    print(f"\n[Baseline] OOF Metrics:")
    for k, v in metrics.items():
        print(f"  {k}: {v}")

    # Save predictions
    pred_df = meta[["sample_id", "age_days", "tissue", "dataset_batch", "intervention"]].copy()
    pred_df["age_weeks_true"] = y_true_wk
    pred_df["age_weeks_pred"] = y_pred_wk
    pred_df["residual_weeks"] = y_true_wk - y_pred_wk
    pred_df.to_csv(OUT / "predictions_elasticnet.csv", index=False)

    # Save result.json
    result = {
        "model_type": "elasticnet",
        "feature_type": "site_top20k",
        "n_features": int(X.shape[1]),
        **metrics,
        "n_samples": int(X.shape[0]),
        "dataset": "GSE120137_thompson_precomputed",
        "note": "Phase 0 prototype — no train/test split, OOF CV only"
    }
    (OUT / "result_elasticnet.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8")
    print(f"\n[Baseline] ✅ Saved → {OUT}")

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
CV-safe embedding-aware methylation clock training.

The benchmark path fits imputation, embedding, optional hybrid feature
selection, scaling, and model training inside each CV fold. The full-data
embedding artifacts written by this script are for interpretation only and are
not used for out-of-fold metrics.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.base import clone
from sklearn.decomposition import MiniBatchNMF, TruncatedSVD
from sklearn.impute import SimpleImputer
from sklearn.linear_model import ElasticNetCV, RidgeCV
from sklearn.model_selection import GroupKFold, KFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

try:
    import lightgbm as lgb
except ImportError:  # pragma: no cover - handled at runtime
    lgb = None

warnings.filterwarnings("ignore", message="An ill-conditioned matrix detected")
warnings.filterwarnings("ignore", message="Singular matrix in solving dual problem")
warnings.filterwarnings("ignore", message="X does not have valid feature names")
warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = Path("/home/zdq-as/mouse_methyl_work")
PHASE0 = ROOT / "results" / "phase0"


def compute_metrics(y_true_weeks: np.ndarray, y_pred_weeks: np.ndarray) -> dict:
    if np.std(y_pred_weeks) == 0 or np.std(y_true_weeks) == 0:
        r, pval = 0.0, 1.0
    else:
        r, pval = stats.pearsonr(y_true_weeks, y_pred_weeks)
        if np.isnan(r):
            r, pval = 0.0, 1.0
    mae = float(np.mean(np.abs(y_true_weeks - y_pred_weeks)))
    medae = float(np.median(np.abs(y_true_weeks - y_pred_weeks)))
    rmse = float(np.sqrt(np.mean((y_true_weeks - y_pred_weeks) ** 2)))
    ss_res = float(np.sum((y_true_weeks - y_pred_weeks) ** 2))
    ss_tot = float(np.sum((y_true_weeks - y_true_weeks.mean()) ** 2))
    r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0.0
    return {
        "pearson_r": round(float(r), 4),
        "pearson_pval": float(pval),
        "mae_weeks": round(mae, 3),
        "medae_weeks": round(medae, 3),
        "rmse_weeks": round(rmse, 3),
        "r2": round(float(r2), 4),
    }


def fast_pearson_correlations(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    X_mean = np.nanmean(X, axis=0)
    y_mean = float(np.mean(y))
    X_centered = X - X_mean
    y_centered = y - y_mean
    ss_X = np.nansum(X_centered**2, axis=0)
    ss_y = float(np.sum(y_centered**2))
    cov = np.nansum(X_centered * y_centered[:, np.newaxis], axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        corr = cov / np.sqrt(ss_X * ss_y)
    return np.nan_to_num(corr, nan=0.0, posinf=0.0, neginf=0.0)


def tissue_eta_squared(X: np.ndarray, labels: np.ndarray) -> np.ndarray:
    valid_labels = pd.Series(labels).fillna("unknown").astype(str).values
    unique = np.unique(valid_labels)
    overall = np.nanmean(X, axis=0)
    ss_between = np.zeros(X.shape[1], dtype=np.float64)
    for label in unique:
        mask = valid_labels == label
        if mask.sum() == 0:
            continue
        group_mean = np.nanmean(X[mask], axis=0)
        ss_between += mask.sum() * (group_mean - overall) ** 2
    ss_total = np.nansum((X - overall) ** 2, axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        eta = ss_between / ss_total
    return np.nan_to_num(eta, nan=0.0, posinf=0.0, neginf=0.0)


def select_features(X_train: np.ndarray, y_train: np.ndarray, n_features: int) -> np.ndarray:
    if n_features <= 0:
        return np.array([], dtype=int)
    n_features = min(int(n_features), X_train.shape[1])
    corrs = fast_pearson_correlations(X_train, y_train)
    return np.argsort(np.abs(corrs))[-n_features:]


def make_model(model_type: str, seed: int):
    if model_type == "elasticnet":
        return ElasticNetCV(
            l1_ratio=[0.1, 0.5, 0.7, 0.9, 0.95, 0.99],
            cv=3,
            n_jobs=-1,
            max_iter=10000,
            random_state=seed,
        )
    if model_type == "ridge":
        return RidgeCV(alphas=np.logspace(-2, 4, 25), cv=3)
    if model_type == "lgbm":
        if lgb is None:
            raise RuntimeError("lightgbm is not installed")
        return lgb.LGBMRegressor(
            n_estimators=200,
            learning_rate=0.03,
            num_leaves=31,
            subsample=0.9,
            colsample_bytree=0.9,
            random_state=seed,
            n_jobs=-1,
            verbose=-1,
        )
    raise ValueError(f"Unsupported model_type: {model_type}")


def make_embedder(embedder: str, n_components: int, seed: int):
    if embedder == "svd":
        return TruncatedSVD(n_components=n_components, random_state=seed)
    if embedder == "nmf":
        return MiniBatchNMF(
            n_components=n_components,
            init="nndsvda",
            random_state=seed,
            max_iter=200,
            batch_size=128,
        )
    raise ValueError(f"Unsupported embedder: {embedder}")


def build_cv(groups: np.ndarray, seed: int):
    unique_groups = np.unique(groups)
    if len(unique_groups) >= 2:
        cv = GroupKFold(n_splits=min(5, len(unique_groups)))
        return "GroupKFold_dataset_batch", list(cv.split(np.zeros(len(groups)), groups=groups))
    cv = KFold(n_splits=5, shuffle=True, random_state=seed)
    return "KFold_phase0_within_dataset", list(cv.split(np.zeros(len(groups))))


def load_data(matrix_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    meta_path = ROOT / "metadata" / "unified_sample_metadata.csv"
    if not matrix_path.exists():
        print(f"[ERROR] Matrix missing: {matrix_path}")
        sys.exit(1)
    if not meta_path.exists():
        print(f"[ERROR] Metadata missing: {meta_path}")
        sys.exit(1)
    matrix = pd.read_parquet(matrix_path)
    matrix.columns = [str(c).split("_")[0] for c in matrix.columns]
    meta = pd.read_csv(meta_path)
    meta = meta[meta["age_days"].notna()].copy()
    common = [sample for sample in meta["sample_id"] if sample in matrix.columns]
    meta = meta[meta["sample_id"].isin(common)].reset_index(drop=True)
    if not common:
        print("[ERROR] No common samples between matrix and metadata.")
        sys.exit(1)
    return matrix[common], meta


def model_feature_importance(fitted_pipe: Pipeline, n_embedding: int, n_hybrid: int) -> np.ndarray:
    model = fitted_pipe.named_steps["model"]
    total = n_embedding + n_hybrid
    if hasattr(model, "feature_importances_"):
        imp = np.asarray(model.feature_importances_, dtype=float)
    elif hasattr(model, "coef_"):
        imp = np.abs(np.ravel(model.coef_)).astype(float)
    else:
        imp = np.ones(total, dtype=float)
    if len(imp) != total:
        fixed = np.zeros(total, dtype=float)
        fixed[: min(total, len(imp))] = imp[: min(total, len(imp))]
        imp = fixed
    return imp


def write_embedding_artifacts(
    out_dir: Path,
    X_raw: np.ndarray,
    feature_ids: list[str],
    meta: pd.DataFrame,
    args: argparse.Namespace,
    embedding_dim_importance: np.ndarray,
    raw_feature_importance: np.ndarray,
    raw_feature_counts: np.ndarray,
) -> None:
    imputer = SimpleImputer(strategy=args.imputation)
    X_imp = imputer.fit_transform(X_raw).astype(np.float32)
    embedder = make_embedder(args.embedder, args.n_components, args.random_seed)
    sample_embeddings = embedder.fit_transform(X_imp).astype(np.float32)
    components = np.asarray(embedder.components_, dtype=np.float32).T

    emb_cols = [f"emb_{i:03d}" for i in range(args.n_components)]
    sample_df = pd.DataFrame(sample_embeddings, columns=emb_cols)
    sample_df.insert(0, "sample_id", meta["sample_id"].values)
    sample_df["dataset_batch"] = meta["dataset_batch"].values
    sample_df["tissue"] = meta["tissue"].values
    sample_df["age_days"] = meta["age_days"].values
    sample_df.to_parquet(out_dir / "sample_embeddings.parquet", index=False)

    feature_df = pd.DataFrame(components, columns=emb_cols)
    feature_df.insert(0, "region_id", feature_ids)
    top_idx = np.argmax(np.abs(components), axis=1)
    feature_df["top_latent_factor"] = [emb_cols[i] for i in top_idx]
    feature_df["top_latent_loading"] = components[np.arange(len(feature_ids)), top_idx]
    feature_df.to_parquet(out_dir / "feature_embeddings.parquet", index=False)

    age_corr = fast_pearson_correlations(X_raw, np.log1p(meta["age_days"].values.astype(float)))
    tissue_eta = tissue_eta_squared(X_imp, meta["tissue"].values)
    embedding_region_importance = np.abs(components) @ embedding_dim_importance
    denom = max(float(raw_feature_counts.max()), 1.0)
    stability = raw_feature_counts / denom

    interp = pd.DataFrame(
        {
            "region_id": feature_ids,
            "age_abs_corr": np.abs(age_corr),
            "age_corr": age_corr,
            "tissue_eta2": tissue_eta,
            "fold_selection_count": raw_feature_counts.astype(int),
            "fold_stability": stability,
            "raw_model_importance": raw_feature_importance,
            "embedding_model_importance": embedding_region_importance,
            "combined_importance": embedding_region_importance + raw_feature_importance,
            "top_latent_factor": feature_df["top_latent_factor"].values,
            "top_latent_loading": feature_df["top_latent_loading"].values,
        }
    )

    stats_path = PHASE0 / "region_stats_5kb.csv"
    if stats_path.exists():
        stats_df = pd.read_csv(stats_path)
        interp = interp.merge(stats_df, on="region_id", how="left")
    interp.sort_values("combined_importance", ascending=False).to_csv(
        out_dir / "feature_interpretation.csv", index=False
    )

    cluster_summary = (
        interp.groupby("top_latent_factor", dropna=False)
        .agg(
            n_regions=("region_id", "count"),
            mean_combined_importance=("combined_importance", "mean"),
            max_combined_importance=("combined_importance", "max"),
            mean_age_abs_corr=("age_abs_corr", "mean"),
            mean_tissue_eta2=("tissue_eta2", "mean"),
            mean_fold_stability=("fold_stability", "mean"),
        )
        .sort_values("max_combined_importance", ascending=False)
        .reset_index()
    )
    cluster_summary.to_csv(out_dir / "cluster_summary.csv", index=False)

    top_regions = interp.sort_values("combined_importance", ascending=False).head(20)
    report_lines = [
        "# Embedding Feature Report",
        "",
        "Full-data embeddings in this directory are interpretation artifacts only; benchmark metrics use fold-internal embeddings.",
        "",
        f"- Embedder: `{args.embedder}`",
        f"- Components: `{args.n_components}`",
        f"- Hybrid top regions: `{args.hybrid_top_regions}`",
        f"- Imputation: `{args.imputation}`",
        "",
        "## Top Regions",
        "",
        "| Rank | Region | Combined importance | Age abs corr | Tissue eta2 | Fold stability | Top factor |",
        "|---:|---|---:|---:|---:|---:|---|",
    ]
    for rank, row in enumerate(top_regions.itertuples(index=False), start=1):
        report_lines.append(
            f"| {rank} | {row.region_id} | {row.combined_importance:.6g} | "
            f"{row.age_abs_corr:.4f} | {row.tissue_eta2:.4f} | "
            f"{row.fold_stability:.3f} | {row.top_latent_factor} |"
        )
    report_lines.extend(
        [
            "",
            "## Top Latent Factors",
            "",
            "| Factor | Regions | Max importance | Mean age abs corr | Mean tissue eta2 |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for row in cluster_summary.head(10).itertuples(index=False):
        report_lines.append(
            f"| {row.top_latent_factor} | {row.n_regions} | {row.max_combined_importance:.6g} | "
            f"{row.mean_age_abs_corr:.4f} | {row.mean_tissue_eta2:.4f} |"
        )
    (out_dir / "embedding_feature_report.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix_path", default=str(PHASE0 / "region_matrix_5kb.parquet"))
    parser.add_argument("--embedder", default="svd", choices=["svd", "nmf"])
    parser.add_argument("--n_components", type=int, default=32)
    parser.add_argument("--model_type", default="ridge", choices=["ridge", "elasticnet", "lgbm"])
    parser.add_argument("--hybrid_top_regions", type=int, default=0)
    parser.add_argument("--imputation", default="median", choices=["median", "mean"])
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--random_seed", type=int, default=42)
    parser.add_argument("--randomize_labels", action="store_true")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    matrix, meta = load_data(Path(args.matrix_path))
    feature_ids = [str(idx) for idx in matrix.index]
    X_raw = matrix.T.values.astype(np.float32)
    y_days_true = meta["age_days"].values.astype(np.float32)
    y_log_fit = np.log1p(y_days_true).astype(np.float32)
    if args.randomize_labels:
        rng = np.random.default_rng(args.random_seed)
        y_log_fit = rng.permutation(y_log_fit)

    groups = meta["dataset_batch"].fillna("unknown").astype(str).values
    cv_strategy, cv_iter = build_cv(groups, args.random_seed)
    y_pred_log = np.full(len(meta), np.nan, dtype=np.float32)
    raw_feature_counts = np.zeros(X_raw.shape[1], dtype=np.int32)
    raw_feature_importance = np.zeros(X_raw.shape[1], dtype=np.float64)
    embedding_dim_importance = np.zeros(args.n_components, dtype=np.float64)

    base_pipe = Pipeline(
        [
            ("scaler", StandardScaler()),
            ("model", make_model(args.model_type, args.random_seed)),
        ]
    )

    for fold, (train_idx, test_idx) in enumerate(cv_iter, start=1):
        imputer = SimpleImputer(strategy=args.imputation)
        X_train_imp = imputer.fit_transform(X_raw[train_idx]).astype(np.float32)
        X_test_imp = imputer.transform(X_raw[test_idx]).astype(np.float32)

        embedder = make_embedder(args.embedder, args.n_components, args.random_seed + fold)
        X_train_emb = embedder.fit_transform(X_train_imp).astype(np.float32)
        X_test_emb = embedder.transform(X_test_imp).astype(np.float32)

        selected = select_features(X_raw[train_idx], y_log_fit[train_idx], args.hybrid_top_regions)
        raw_feature_counts[selected] += 1
        if len(selected) > 0:
            X_train_model = np.concatenate([X_train_emb, X_train_imp[:, selected]], axis=1)
            X_test_model = np.concatenate([X_test_emb, X_test_imp[:, selected]], axis=1)
        else:
            X_train_model = X_train_emb
            X_test_model = X_test_emb

        pipe = clone(base_pipe)
        pipe.fit(X_train_model, y_log_fit[train_idx])
        y_pred_log[test_idx] = pipe.predict(X_test_model)

        imp = model_feature_importance(pipe, args.n_components, len(selected))
        embedding_dim_importance += imp[: args.n_components]
        if len(selected) > 0:
            raw_feature_importance[selected] += imp[args.n_components :]
        print(
            f"  Fold {fold}: train n={len(train_idx)} test n={len(test_idx)} "
            f"embed={args.embedder}{args.n_components} hybrid_regions={len(selected)}"
        )

    y_pred_log = np.clip(y_pred_log, 0, np.log1p(max(float(y_days_true.max()) * 1.5, 1.0)))
    y_true_weeks = y_days_true / 7
    y_pred_days = np.expm1(y_pred_log)
    y_pred_weeks = y_pred_days / 7
    metrics = compute_metrics(y_true_weeks, y_pred_weeks)

    embedding_dim_importance /= max(len(cv_iter), 1)
    raw_feature_importance /= max(len(cv_iter), 1)

    write_embedding_artifacts(
        out_dir=out_dir,
        X_raw=X_raw,
        feature_ids=feature_ids,
        meta=meta,
        args=args,
        embedding_dim_importance=embedding_dim_importance,
        raw_feature_importance=raw_feature_importance,
        raw_feature_counts=raw_feature_counts,
    )

    result = {
        "model_type": args.model_type,
        "feature_type": "embedding_hybrid" if args.hybrid_top_regions > 0 else "embedding",
        "matrix_path": str(Path(args.matrix_path)),
        "embedder": args.embedder,
        "n_components": int(args.n_components),
        "hybrid_top_regions": int(args.hybrid_top_regions),
        "imputation": args.imputation,
        "cv_strategy": cv_strategy,
        "random_seed": int(args.random_seed),
        "randomize_labels": bool(args.randomize_labels),
        "n_samples": int(X_raw.shape[0]),
        "n_features_total": int(X_raw.shape[1]),
        "dataset_scope": "phase0_GSE120137_only_embedding",
        "leakage_controls": [
            "imputer_fit_train_fold_only",
            "embedder_fit_train_fold_only",
            "hybrid_feature_selection_fit_train_fold_only",
            "scaler_fit_train_fold_only",
            "model_fit_train_fold_only",
            "full_data_embedding_artifacts_for_interpretation_only",
        ],
        "exec_time_sec": round(time.time() - t0, 1),
        **metrics,
    }
    (out_dir / "result.json").write_text(json.dumps(result, indent=2), encoding="utf-8")

    pred_df = pd.DataFrame(
        {
            "sample_id": meta["sample_id"].values,
            "dataset_batch": meta.get("dataset_batch", pd.Series(["unknown"] * len(meta))).values,
            "tissue": meta.get("tissue", pd.Series(["unknown"] * len(meta))).values,
            "intervention": meta.get("intervention", pd.Series(["control"] * len(meta))).values,
            "age_days_true": y_days_true,
            "age_days_pred": y_pred_days,
            "age_weeks_true": y_true_weeks,
            "age_weeks_pred": y_pred_weeks,
            "residual_weeks": y_true_weeks - y_pred_weeks,
        }
    )
    pred_df.to_csv(out_dir / "predictions.csv", index=False)
    print(f"[{args.embedder}{args.n_components}+{args.model_type}] R:{metrics['pearson_r']} MAE:{metrics['mae_weeks']}wk R2:{metrics['r2']}")


if __name__ == "__main__":
    main()

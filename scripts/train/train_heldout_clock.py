#!/usr/bin/env python3
"""
Train on one methylation matrix and predict a held-out dataset matrix.

This is the v5 GSE120137 -> GSE80672 validation path. Feature selection,
imputation, scaling, and model fitting are all trained on the training dataset
only; the held-out dataset is transformed using those train-fitted objects.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.ensemble import RandomForestRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import ElasticNetCV, RidgeCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import QuantileTransformer, RobustScaler, StandardScaler

try:
    import lightgbm as lgb
except ImportError:  # pragma: no cover
    lgb = None

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


def select_features(X_train: np.ndarray, y_train: np.ndarray, n_features: int | str) -> np.ndarray:
    if str(n_features).lower() == "all":
        n_features = X_train.shape[1]
    n_features = min(int(n_features), X_train.shape[1])
    corrs = fast_pearson_correlations(X_train, y_train)
    return np.argsort(np.abs(corrs))[-n_features:]


def train_dataset_stability_mask(
    X_train: np.ndarray,
    groups_train: np.ndarray,
    min_group_presence: float | None,
    max_group_mean_shift: float | None,
) -> np.ndarray:
    """Return a train-only feature mask for dataset-level coverage/shift stability."""
    mask = np.ones(X_train.shape[1], dtype=bool)
    unique_groups = [group for group in pd.unique(groups_train) if pd.notna(group)]
    if len(unique_groups) < 2 and max_group_mean_shift is None:
        return mask

    group_means = []
    for group in unique_groups:
        group_X = X_train[groups_train == group]
        finite = np.isfinite(group_X)
        if min_group_presence is not None:
            mask &= finite.mean(axis=0) >= min_group_presence
        if max_group_mean_shift is not None:
            with np.errstate(invalid="ignore"):
                group_mean = np.nanmean(group_X, axis=0)
            group_means.append(group_mean)

    if max_group_mean_shift is not None and group_means:
        means = np.vstack(group_means)
        finite_means = np.isfinite(means).all(axis=0)
        with np.errstate(invalid="ignore"):
            shift = np.nanmax(means, axis=0) - np.nanmin(means, axis=0)
        mask &= finite_means & (shift <= max_group_mean_shift)
    return mask


def train_agebin_presence_mask(
    X_train: np.ndarray,
    age_weeks_train: np.ndarray,
    min_agebin_presence: float | None,
    min_agebin_samples: int,
) -> np.ndarray:
    """Return a train-only feature mask requiring coverage across populated age bins."""
    mask = np.ones(X_train.shape[1], dtype=bool)
    if min_agebin_presence is None:
        return mask
    bins = np.array([0, 4, 13, 26, 52, 104, np.inf], dtype=float)
    bin_ids = np.digitize(age_weeks_train, bins[1:-1], right=False)
    populated = [bin_id for bin_id in np.unique(bin_ids) if np.sum(bin_ids == bin_id) >= min_agebin_samples]
    if len(populated) < 2:
        return mask
    finite = np.isfinite(X_train)
    for bin_id in populated:
        mask &= finite[bin_ids == bin_id].mean(axis=0) >= min_agebin_presence
    return mask


def make_preprocess_steps(preprocess: str, imputation: str, seed: int, n_train: int) -> list[tuple]:
    steps: list[tuple] = [("imputer", SimpleImputer(strategy=imputation))]
    if preprocess == "standard":
        steps.append(("scaler", StandardScaler()))
    elif preprocess == "robust":
        steps.append(("scaler", RobustScaler()))
    elif preprocess == "quantile_uniform":
        steps.append(
            (
                "scaler",
                QuantileTransformer(
                    n_quantiles=max(10, min(1000, n_train)),
                    output_distribution="uniform",
                    random_state=seed,
                ),
            )
        )
    elif preprocess == "none":
        pass
    else:
        raise ValueError(f"Unsupported preprocess: {preprocess}")
    return steps


def build_target(age_days: np.ndarray, target_transform: str) -> np.ndarray:
    if target_transform == "log1p_days":
        return np.log1p(age_days).astype(np.float32)
    if target_transform == "linear_weeks":
        return (age_days / 7).astype(np.float32)
    raise ValueError(f"Unsupported target_transform: {target_transform}")


def inverse_target(pred: np.ndarray, target_transform: str, max_age_days: float) -> tuple[np.ndarray, np.ndarray]:
    if target_transform == "log1p_days":
        pred = np.clip(pred, 0, np.log1p(max(max_age_days * 1.5, 1.0)))
        pred_days = np.expm1(pred)
        return pred_days, pred_days / 7
    if target_transform == "linear_weeks":
        pred_weeks = np.clip(pred, 0, max(max_age_days / 7 * 1.5, 1.0))
        return pred_weeks * 7, pred_weeks
    raise ValueError(f"Unsupported target_transform: {target_transform}")


def make_model(
    model_type: str,
    seed: int,
    *,
    lgbm_n_estimators: int = 200,
    lgbm_learning_rate: float = 0.03,
    lgbm_num_leaves: int = 31,
    lgbm_subsample: float = 0.9,
    lgbm_colsample_bytree: float = 0.9,
    rf_n_estimators: int = 150,
    rf_max_depth: int | None = 10,
    rf_min_samples_leaf: int = 2,
    rf_max_features: float | str | None = 1.0,
):
    if model_type == "elasticnet":
        return ElasticNetCV(
            l1_ratio=[0.1, 0.5, 0.9],
            alphas=np.logspace(-2, 2, 12),
            cv=3,
            n_jobs=-1,
            max_iter=3000,
            tol=1e-3,
            random_state=seed,
        )
    if model_type == "ridge":
        return RidgeCV(alphas=np.logspace(-2, 4, 25), cv=3)
    if model_type == "lgbm":
        if lgb is None:
            raise RuntimeError("lightgbm is not installed")
        return lgb.LGBMRegressor(
            n_estimators=lgbm_n_estimators,
            learning_rate=lgbm_learning_rate,
            num_leaves=lgbm_num_leaves,
            subsample=lgbm_subsample,
            colsample_bytree=lgbm_colsample_bytree,
            random_state=seed,
            n_jobs=-1,
            verbose=-1,
        )
    if model_type == "rf":
        return RandomForestRegressor(
            n_estimators=rf_n_estimators,
            max_depth=rf_max_depth,
            min_samples_leaf=rf_min_samples_leaf,
            max_features=rf_max_features,
            random_state=seed,
            n_jobs=-1,
        )
    raise ValueError(f"Unsupported model_type: {model_type}")


def compute_metrics(y_true_weeks: np.ndarray, y_pred_weeks: np.ndarray) -> dict:
    if len(y_true_weeks) < 3 or np.std(y_true_weeks) == 0 or np.std(y_pred_weeks) == 0:
        r, pval = 0.0, 1.0
    else:
        r, pval = stats.pearsonr(y_true_weeks, y_pred_weeks)
    mae = float(np.mean(np.abs(y_true_weeks - y_pred_weeks)))
    medae = float(np.median(np.abs(y_true_weeks - y_pred_weeks)))
    rmse = float(np.sqrt(np.mean((y_true_weeks - y_pred_weeks) ** 2)))
    ss_res = float(np.sum((y_true_weeks - y_pred_weeks) ** 2))
    ss_tot = float(np.sum((y_true_weeks - y_true_weeks.mean()) ** 2))
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return {
        "pearson_r": round(float(r), 4),
        "pearson_pval": float(pval),
        "mae_weeks": round(mae, 3),
        "medae_weeks": round(medae, 3),
        "rmse_weeks": round(rmse, 3),
        "r2": round(float(r2), 4),
    }


def dataset_filter(meta: pd.DataFrame, selector: str) -> pd.Series:
    selector = selector.strip()
    if selector.startswith("all_except:"):
        excluded = {item.strip() for item in selector.split(":", 1)[1].split(",") if item.strip()}
        return ~meta["dataset_batch"].isin(excluded)
    if selector in {"all", "*"}:
        return pd.Series(True, index=meta.index)
    allowed = {item.strip() for item in selector.split(",") if item.strip()}
    return meta["dataset_batch"].isin(allowed)


def load_matrix_and_meta(matrix_path: Path, dataset_selector: str, metadata_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    if not matrix_path.exists():
        print(f"[ERROR] Matrix missing: {matrix_path}")
        sys.exit(1)
    matrix = pd.read_parquet(matrix_path)
    matrix.columns = [resolve_matrix_sample_id(col) for col in matrix.columns]
    if pd.Index(matrix.columns).duplicated().any():
        duplicates = sorted(pd.Index(matrix.columns)[pd.Index(matrix.columns).duplicated()].unique())
        raise SystemExit(f"Duplicate resolved matrix sample ids: {duplicates[:10]}")
    meta = pd.read_csv(metadata_path)
    meta = meta[dataset_filter(meta, dataset_selector) & meta["age_days"].notna()].copy()
    common = [sample for sample in meta["sample_id"] if sample in matrix.columns]
    meta = meta[meta["sample_id"].isin(common)].reset_index(drop=True)
    if not common:
        print(f"[ERROR] No common samples between {matrix_path} and metadata for {dataset_selector}")
        sys.exit(1)
    return matrix[common], meta


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--train_matrix", default=str(ROOT / "results" / "phase0" / "region_matrix_5kb.parquet"))
    parser.add_argument("--test_matrix", default=str(ROOT / "results" / "multidataset" / "GSE80672_region_matrix_5kb.parquet"))
    parser.add_argument("--train_dataset", default="GSE120137")
    parser.add_argument("--test_dataset", default="GSE80672")
    parser.add_argument("--train_datasets", default=None, help="Comma list, 'all', or 'all_except:GSE80672'.")
    parser.add_argument("--test_datasets", default=None, help="Comma list or single dataset. Defaults to --test_dataset.")
    parser.add_argument("--metadata_path", default=str(META_FILE))
    parser.add_argument("--model_type", default="lgbm", choices=["ridge", "elasticnet", "lgbm", "rf"])
    parser.add_argument("--n_feature_prefilter", default="500")
    parser.add_argument("--imputation", default="median", choices=["median", "mean"])
    parser.add_argument("--preprocess", default="standard", choices=["standard", "robust", "quantile_uniform", "none"])
    parser.add_argument("--min_train_feature_presence", type=float, default=0.0)
    parser.add_argument("--min_train_group_feature_presence", type=float, default=None)
    parser.add_argument("--max_train_dataset_mean_shift", type=float, default=None)
    parser.add_argument("--min_train_agebin_feature_presence", type=float, default=None)
    parser.add_argument("--min_train_agebin_samples", type=int, default=8)
    parser.add_argument("--target_transform", default="log1p_days", choices=["log1p_days", "linear_weeks"])
    parser.add_argument("--lgbm_n_estimators", type=int, default=200)
    parser.add_argument("--lgbm_learning_rate", type=float, default=0.03)
    parser.add_argument("--lgbm_num_leaves", type=int, default=31)
    parser.add_argument("--lgbm_subsample", type=float, default=0.9)
    parser.add_argument("--lgbm_colsample_bytree", type=float, default=0.9)
    parser.add_argument("--rf_n_estimators", type=int, default=150)
    parser.add_argument("--rf_max_depth", default="10", help="Integer depth or 'none'.")
    parser.add_argument("--rf_min_samples_leaf", type=int, default=2)
    parser.add_argument("--rf_max_features", default="1.0", help="Float, 'sqrt', 'log2', or 'none'.")
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--random_seed", type=int, default=42)
    parser.add_argument("--randomize_train_labels", action="store_true")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    started = time.time()

    train_selector = args.train_datasets or args.train_dataset
    test_selector = args.test_datasets or args.test_dataset
    metadata_path = Path(args.metadata_path)
    train_matrix, train_meta = load_matrix_and_meta(Path(args.train_matrix), train_selector, metadata_path)
    test_matrix, test_meta = load_matrix_and_meta(Path(args.test_matrix), test_selector, metadata_path)
    common_features = train_matrix.index.intersection(test_matrix.index)
    if len(common_features) == 0:
        print("[ERROR] No common features between train and held-out matrices.")
        sys.exit(1)
    train_matrix = train_matrix.loc[common_features]
    test_matrix = test_matrix.loc[common_features]

    X_train_raw = train_matrix.T.values.astype(np.float32)
    X_test_raw = test_matrix.T.values.astype(np.float32)
    y_train_days = train_meta["age_days"].values.astype(np.float32)
    y_test_days = test_meta["age_days"].values.astype(np.float32)
    y_train_fit = build_target(y_train_days, args.target_transform)
    y_train_days_for_filters = y_train_days.copy()
    if args.randomize_train_labels:
        rng = np.random.default_rng(args.random_seed)
        perm = rng.permutation(len(y_train_fit))
        y_train_fit = y_train_fit[perm]
        y_train_days_for_filters = y_train_days_for_filters[perm]

    rf_max_depth = None if str(args.rf_max_depth).lower() in {"none", "null"} else int(args.rf_max_depth)
    rf_max_features: float | str | None
    if str(args.rf_max_features).lower() in {"none", "null"}:
        rf_max_features = None
    elif str(args.rf_max_features).lower() in {"sqrt", "log2"}:
        rf_max_features = str(args.rf_max_features).lower()
    else:
        rf_max_features = float(args.rf_max_features)

    train_presence = np.isfinite(X_train_raw).mean(axis=0) >= args.min_train_feature_presence
    stability_count = None
    if args.min_train_group_feature_presence is not None or args.max_train_dataset_mean_shift is not None:
        stability_mask = train_dataset_stability_mask(
            X_train_raw,
            train_meta["dataset_batch"].fillna("unknown").astype(str).values,
            args.min_train_group_feature_presence,
            args.max_train_dataset_mean_shift,
        )
        train_presence &= stability_mask
        stability_count = int(stability_mask.sum())
    agebin_presence_count = None
    if args.min_train_agebin_feature_presence is not None:
        agebin_mask = train_agebin_presence_mask(
            X_train_raw,
            y_train_days_for_filters / 7,
            args.min_train_agebin_feature_presence,
            args.min_train_agebin_samples,
        )
        train_presence &= agebin_mask
        agebin_presence_count = int(agebin_mask.sum())
    presence_idx = np.flatnonzero(train_presence)
    if len(presence_idx) == 0:
        raise RuntimeError(f"No features pass min_train_feature_presence={args.min_train_feature_presence}")
    top_local_idx = select_features(X_train_raw[:, presence_idx], y_train_fit, args.n_feature_prefilter)
    top_idx = presence_idx[top_local_idx]
    pipe = Pipeline(
        [
            *make_preprocess_steps(args.preprocess, args.imputation, args.random_seed, X_train_raw.shape[0]),
            (
                "model",
                make_model(
                    args.model_type,
                    args.random_seed,
                    lgbm_n_estimators=args.lgbm_n_estimators,
                    lgbm_learning_rate=args.lgbm_learning_rate,
                    lgbm_num_leaves=args.lgbm_num_leaves,
                    lgbm_subsample=args.lgbm_subsample,
                    lgbm_colsample_bytree=args.lgbm_colsample_bytree,
                    rf_n_estimators=args.rf_n_estimators,
                    rf_max_depth=rf_max_depth,
                    rf_min_samples_leaf=args.rf_min_samples_leaf,
                    rf_max_features=rf_max_features,
                ),
            ),
        ]
    )
    pipe.fit(X_train_raw[:, top_idx], y_train_fit)
    y_pred_fit = pipe.predict(X_test_raw[:, top_idx])
    y_pred_days, y_pred_weeks = inverse_target(
        y_pred_fit,
        args.target_transform,
        max(float(y_train_days.max()), float(y_test_days.max())),
    )
    y_true_weeks = y_test_days / 7
    metrics = compute_metrics(y_true_weeks, y_pred_weeks)

    pred_df = pd.DataFrame(
        {
            "sample_id": test_meta["sample_id"].values,
            "dataset_batch": test_meta["dataset_batch"].values,
            "tissue": test_meta["tissue"].values,
            "intervention": test_meta["intervention"].values,
            "age_days_true": y_test_days,
            "age_days_pred": y_pred_days,
            "age_weeks_true": y_true_weeks,
            "age_weeks_pred": y_pred_weeks,
            "residual_weeks": y_true_weeks - y_pred_weeks,
        }
    )
    for optional_column in [
        "sex",
        "sample_id_source",
        "metadata_source",
        "condition_detail",
        "condition_family",
        "matrix_sample_id",
        "source_column",
    ]:
        if optional_column in test_meta.columns:
            pred_df[optional_column] = test_meta[optional_column].values
    pred_df.to_csv(out_dir / "predictions.csv", index=False)

    result = {
        "model_type": args.model_type,
        "feature_type": "heldout_region",
        "train_dataset": train_selector,
        "test_dataset": test_selector,
        "train_matrix": str(Path(args.train_matrix)),
        "test_matrix": str(Path(args.test_matrix)),
        "metadata_path": str(metadata_path),
        "sample_id_resolver": "gsm_or_gse60012_tile_short_id",
        "n_train_samples": int(X_train_raw.shape[0]),
        "n_test_samples": int(X_test_raw.shape[0]),
        "n_common_features": int(len(common_features)),
        "n_feature_prefilter": args.n_feature_prefilter,
        "n_features_passing_train_presence": int(len(presence_idx)),
        "n_features_selected": int(len(top_idx)),
        "imputation": args.imputation,
        "preprocess": args.preprocess,
        "min_train_feature_presence": float(args.min_train_feature_presence),
        "min_train_group_feature_presence": (
            None if args.min_train_group_feature_presence is None else float(args.min_train_group_feature_presence)
        ),
        "max_train_dataset_mean_shift": (
            None if args.max_train_dataset_mean_shift is None else float(args.max_train_dataset_mean_shift)
        ),
        "min_train_agebin_feature_presence": (
            None
            if args.min_train_agebin_feature_presence is None
            else float(args.min_train_agebin_feature_presence)
        ),
        "min_train_agebin_samples": int(args.min_train_agebin_samples),
        "n_features_stability": stability_count,
        "n_features_agebin_presence": agebin_presence_count,
        "target_transform": args.target_transform,
        "lgbm_n_estimators": int(args.lgbm_n_estimators),
        "lgbm_learning_rate": float(args.lgbm_learning_rate),
        "lgbm_num_leaves": int(args.lgbm_num_leaves),
        "lgbm_subsample": float(args.lgbm_subsample),
        "lgbm_colsample_bytree": float(args.lgbm_colsample_bytree),
        "rf_n_estimators": int(args.rf_n_estimators),
        "rf_max_depth": rf_max_depth,
        "rf_min_samples_leaf": int(args.rf_min_samples_leaf),
        "rf_max_features": rf_max_features,
        "random_seed": int(args.random_seed),
        "randomize_train_labels": bool(args.randomize_train_labels),
        "leakage_controls": [
            "feature_intersection_before_train_only_selection",
            "feature_selection_fit_train_dataset_only",
            "imputer_fit_train_dataset_only",
            "scaler_fit_train_dataset_only",
            "model_fit_train_dataset_only",
        ],
        "exec_time_sec": round(time.time() - started, 1),
        **metrics,
    }
    (out_dir / "result.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(
        f"[Heldout] {train_selector}->{test_selector} {args.model_type} "
        f"R:{metrics['pearson_r']} MAE:{metrics['mae_weeks']}wk R2:{metrics['r2']}"
    )


if __name__ == "__main__":
    main()

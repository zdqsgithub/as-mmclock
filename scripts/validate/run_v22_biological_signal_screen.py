#!/usr/bin/env python3
"""Run v22 leakage-controlled biological signal screens on mouse methylation matrices."""
from __future__ import annotations

import argparse
import json
import math
import re
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression, RidgeCV
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    f1_score,
    mean_absolute_error,
    mean_squared_error,
    r2_score,
    roc_auc_score,
)
from sklearn.model_selection import GroupKFold, KFold, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import LabelEncoder, RobustScaler, StandardScaler

try:
    import lightgbm as lgb
except ImportError:  # pragma: no cover - runtime capability
    lgb = None


LOCAL_ROOT = Path("/home/zdq-as/mouse_methyl_work")
REMOTE_ROOT = Path("/root/autodl-tmp/mouse_methyl_work")


def path_exists_safe(path: Path) -> bool:
    try:
        return path.exists()
    except PermissionError:
        return False


ROOT = REMOTE_ROOT if path_exists_safe(REMOTE_ROOT) else LOCAL_ROOT
DEFAULT_MATRIX = ROOT / "results" / "multidataset_v8_3_ablation" / "all6" / "all_rrbs_region_matrix_5kb.parquet"
DEFAULT_METADATA = ROOT / "metadata" / "model_sample_metadata_v8.csv"
DEFAULT_TARGET_REGISTRY = ROOT / "metadata" / "v22_biological_signal_targets.csv"
DEFAULT_OUT = ROOT / "results" / "v22_biological_signal_screen"
DOC_REPORT = ROOT / "doc" / "20_analysis" / "56_20260521_v22_biological_signal_screen_report.md"
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
    return {
        "chrom": chrom,
        "start": start,
        "end": end,
        "cluster_5mb": f"{chrom}:{cluster_start}-{cluster_start + 4_999_999}",
    }


def safe_float(value: object, default: float = 0.0) -> float:
    try:
        if pd.isna(value):
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def safe_int(value: object, default: int = 0) -> int:
    try:
        if pd.isna(value):
            return default
        return int(float(value))
    except (TypeError, ValueError):
        return default


def safe_corr(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    valid = np.isfinite(x) & np.isfinite(y)
    if valid.sum() < 5 or np.nanstd(x[valid]) == 0 or np.nanstd(y[valid]) == 0:
        return 0.0, 1.0
    r, p = stats.pearsonr(x[valid], y[valid])
    if not np.isfinite(r):
        return 0.0, 1.0
    return float(r), float(p)


def cramers_v(labels: pd.Series, groups: pd.Series) -> float:
    table = pd.crosstab(labels.fillna("NA").astype(str), groups.fillna("NA").astype(str))
    if table.shape[0] < 2 or table.shape[1] < 2 or table.to_numpy().sum() == 0:
        return 0.0
    chi2 = stats.chi2_contingency(table, correction=False)[0]
    n = table.to_numpy().sum()
    denom = n * (min(table.shape) - 1)
    return float(math.sqrt(chi2 / denom)) if denom > 0 else 0.0


def eta_squared_numeric_by_label(values: pd.Series, labels: pd.Series) -> float:
    df = pd.DataFrame({"value": pd.to_numeric(values, errors="coerce"), "label": labels}).dropna()
    if df["label"].nunique() < 2 or len(df) < 5:
        return 0.0
    grand = df["value"].mean()
    ss_between = 0.0
    ss_total = float(((df["value"] - grand) ** 2).sum())
    for _, group in df.groupby("label"):
        ss_between += float(len(group) * (group["value"].mean() - grand) ** 2)
    return float(ss_between / ss_total) if ss_total > 0 else 0.0


def fast_pearson_scores(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    y = y.astype(float)
    X_mean = np.nanmean(X, axis=0)
    y_mean = float(np.mean(y))
    X_centered = X - X_mean
    y_centered = y - y_mean
    ss_x = np.nansum(X_centered**2, axis=0)
    ss_y = float(np.sum(y_centered**2))
    cov = np.nansum(X_centered * y_centered[:, np.newaxis], axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        scores = cov / np.sqrt(ss_x * ss_y)
    return np.nan_to_num(scores, nan=0.0, posinf=0.0, neginf=0.0)


def fast_multiclass_scores(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    classes = np.unique(y)
    if len(classes) < 2:
        return np.zeros(X.shape[1], dtype=float)
    overall = np.nanmean(X, axis=0)
    ss_between = np.zeros(X.shape[1], dtype=float)
    ss_within = np.zeros(X.shape[1], dtype=float)
    n_total = np.zeros(X.shape[1], dtype=float)
    for cls in classes:
        Xc = X[y == cls]
        valid = np.isfinite(Xc)
        n_c = valid.sum(axis=0).astype(float)
        with np.errstate(invalid="ignore"):
            mean_c = np.nanmean(Xc, axis=0)
        diff = np.nan_to_num(mean_c - overall, nan=0.0)
        ss_between += n_c * diff**2
        ss_within += np.nansum((Xc - mean_c) ** 2, axis=0)
        n_total += n_c
    df_between = max(1, len(classes) - 1)
    df_within = np.maximum(1, n_total - len(classes))
    with np.errstate(divide="ignore", invalid="ignore"):
        f_scores = (ss_between / df_between) / (ss_within / df_within)
    return np.nan_to_num(f_scores, nan=0.0, posinf=0.0, neginf=0.0)


def parse_values(text: object) -> set[str]:
    if pd.isna(text) or str(text).strip() == "":
        return set()
    return {item.strip() for item in str(text).split(";") if item.strip()}


@dataclass
class TargetData:
    target_id: str
    target_kind: str
    meta: pd.DataFrame
    y: np.ndarray
    y_labels: np.ndarray
    class_names: list[str]
    status: str
    reason: str


def load_matrix_and_metadata(matrix_path: Path, metadata_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    if not matrix_path.exists():
        raise SystemExit(f"Matrix missing: {matrix_path}")
    if not metadata_path.exists():
        raise SystemExit(f"Metadata missing: {metadata_path}")
    matrix = pd.read_parquet(matrix_path)
    matrix.columns = [resolve_matrix_sample_id(col) for col in matrix.columns]
    if pd.Index(matrix.columns).duplicated().any():
        dupes = sorted(pd.Index(matrix.columns)[pd.Index(matrix.columns).duplicated()].unique())
        raise SystemExit(f"Duplicate resolved matrix sample IDs: {dupes[:10]}")
    meta = pd.read_csv(metadata_path)
    common = [sample for sample in meta["sample_id"].astype(str) if sample in matrix.columns]
    meta = meta[meta["sample_id"].astype(str).isin(common)].copy().reset_index(drop=True)
    matrix = matrix[common]
    return matrix, meta


def apply_dataset_scope(meta: pd.DataFrame, scope: str) -> pd.Series:
    scope = str(scope).strip()
    if scope in {"", "all"}:
        return pd.Series(True, index=meta.index)
    if scope == "all_non_unknown":
        return pd.Series(True, index=meta.index)
    allowed = {item.strip() for item in scope.split(";") if item.strip()}
    return meta["dataset_batch"].astype(str).isin(allowed)


def build_target(row: pd.Series, meta: pd.DataFrame, *, smoke: bool) -> TargetData:
    target_id = str(row["target_id"])
    kind = str(row["target_kind"])
    label_col = str(row["label_column"])
    min_total = safe_int(row.get("min_total"), 20)
    min_class = max(2, safe_int(row.get("min_class_count"), 5))
    min_groups = max(1, safe_int(row.get("min_groups"), 1))

    scoped = meta[apply_dataset_scope(meta, str(row.get("dataset_scope", "all")))].copy()
    if kind == "regression":
        scoped = scoped[pd.to_numeric(scoped[label_col], errors="coerce").notna()].copy()
        if label_col == "age_days":
            y = pd.to_numeric(scoped[label_col], errors="coerce").to_numpy(dtype=float) / 7.0
            y_labels = np.array(["age_weeks"] * len(scoped), dtype=object)
        else:
            y = pd.to_numeric(scoped[label_col], errors="coerce").to_numpy(dtype=float)
            y_labels = np.array([label_col] * len(scoped), dtype=object)
        status, reason = gate_regression(scoped, y, min_total, min_groups)
        return TargetData(target_id, kind, scoped, y, y_labels, [], status, reason)

    include_values = parse_values(row.get("include_values"))
    if include_values:
        scoped = scoped[scoped[label_col].astype(str).isin(include_values)].copy()
    else:
        scoped = scoped[scoped[label_col].notna()].copy()
        scoped = scoped[~scoped[label_col].astype(str).isin({"", "nan", "unknown", "NA"})].copy()

    if kind == "binary":
        positive = str(row["positive_label"])
        negative = str(row["negative_label"])
        scoped = scoped[scoped[label_col].astype(str).isin({positive, negative})].copy()
        y_labels = scoped[label_col].astype(str).to_numpy()
        y = (y_labels == positive).astype(int)
        class_names = [negative, positive]
    elif kind == "multiclass":
        counts = scoped[label_col].astype(str).value_counts()
        keep = set(counts[counts >= min_class].index)
        scoped = scoped[scoped[label_col].astype(str).isin(keep)].copy()
        y_labels = scoped[label_col].astype(str).to_numpy()
        encoder = LabelEncoder()
        y = encoder.fit_transform(y_labels)
        class_names = [str(item) for item in encoder.classes_]
    else:
        return TargetData(target_id, kind, scoped, np.array([]), np.array([]), [], "skipped", f"unsupported target_kind={kind}")

    if smoke and len(scoped) > 320:
        scoped, y, y_labels = stratified_downsample(scoped, y, y_labels, max_n=320)
    status, reason = gate_classification(scoped, y_labels, min_total, min_class, min_groups)
    return TargetData(target_id, kind, scoped, y, y_labels, class_names, status, reason)


def stratified_downsample(
    meta: pd.DataFrame, y: np.ndarray, y_labels: np.ndarray, *, max_n: int
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(22)
    pieces = []
    indices = np.arange(len(meta))
    per_label = max(10, max_n // max(1, len(np.unique(y_labels))))
    for label in np.unique(y_labels):
        idx = indices[y_labels == label]
        if len(idx) > per_label:
            idx = rng.choice(idx, size=per_label, replace=False)
        pieces.extend(idx.tolist())
    if len(pieces) > max_n:
        pieces = rng.choice(np.array(pieces), size=max_n, replace=False).tolist()
    pieces = sorted(pieces)
    return meta.iloc[pieces].reset_index(drop=True), y[pieces], y_labels[pieces]


def gate_regression(meta: pd.DataFrame, y: np.ndarray, min_total: int, min_groups: int) -> tuple[str, str]:
    if len(meta) < min_total:
        return "skipped", f"n_samples={len(meta)} < min_total={min_total}"
    if np.isfinite(y).sum() < min_total:
        return "skipped", "insufficient numeric target coverage"
    if pd.Series(meta["dataset_batch"]).nunique() < min_groups:
        return "passed_within_dataset", "single-dataset or low-group regression; use within-scope CV"
    return "passed", "target passed sample and group gates"


def gate_classification(
    meta: pd.DataFrame, y_labels: np.ndarray, min_total: int, min_class: int, min_groups: int
) -> tuple[str, str]:
    if len(meta) < min_total:
        return "skipped", f"n_samples={len(meta)} < min_total={min_total}"
    counts = pd.Series(y_labels).value_counts()
    if len(counts) < 2:
        return "skipped", "fewer than two classes after filtering"
    if int(counts.min()) < min_class:
        return "skipped", f"minimum class count={int(counts.min())} < min_class_count={min_class}"
    if pd.Series(meta["dataset_batch"]).nunique() < min_groups:
        return "passed_within_dataset", "single-dataset classification; use stratified within-scope CV"
    return "passed", "target passed sample and class gates"


def make_cv(target: TargetData, seed: int) -> tuple[str, list[tuple[np.ndarray, np.ndarray]]]:
    n = len(target.meta)
    groups = target.meta["dataset_batch"].fillna("unknown").astype(str).to_numpy()
    if target.target_kind == "regression":
        n_groups = len(np.unique(groups))
        if n_groups >= 2:
            cv = GroupKFold(n_splits=min(5, n_groups))
            return "GroupKFold_dataset_batch", list(cv.split(np.zeros(n), groups=groups))
        n_splits = min(5, max(2, n // 8))
        cv = KFold(n_splits=n_splits, shuffle=True, random_state=seed)
        return "KFold_within_scope", list(cv.split(np.zeros(n)))

    labels = np.asarray(target.y_labels)
    counts = pd.Series(labels).value_counts()
    min_count = int(counts.min()) if not counts.empty else 0
    n_groups = len(np.unique(groups))
    group_class = pd.crosstab(pd.Series(groups), pd.Series(labels))
    class_group_counts = (group_class > 0).sum(axis=0)
    if n_groups >= 2 and int(class_group_counts.min()) >= 2:
        n_splits = min(5, n_groups, int(class_group_counts.min()))
        cv = GroupKFold(n_splits=n_splits)
        return "GroupKFold_dataset_batch_class_covered", list(cv.split(np.zeros(n), labels, groups=groups))
    n_splits = min(5, min_count)
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    return "StratifiedKFold_within_scope", list(cv.split(np.zeros(n), labels))


def build_model(target_kind: str, model_name: str, seed: int, n_classes: int):
    if target_kind == "regression":
        if model_name == "ridge":
            return RidgeCV(alphas=np.logspace(-2, 4, 16))
        if model_name == "rf":
            return RandomForestRegressor(
                n_estimators=160,
                max_depth=12,
                min_samples_leaf=3,
                max_features="sqrt",
                n_jobs=-1,
                random_state=seed,
            )
        if model_name == "lgbm":
            if lgb is None:
                raise RuntimeError("lightgbm is not installed")
            return lgb.LGBMRegressor(
                n_estimators=180,
                learning_rate=0.03,
                num_leaves=31,
                subsample=0.9,
                colsample_bytree=0.9,
                random_state=seed,
                n_jobs=-1,
                verbose=-1,
            )
    else:
        if model_name == "logistic":
            return LogisticRegression(
                max_iter=1200,
                class_weight="balanced",
                solver="lbfgs",
            )
        if model_name == "rf":
            return RandomForestClassifier(
                n_estimators=180,
                max_depth=12,
                min_samples_leaf=3,
                max_features="sqrt",
                class_weight="balanced",
                n_jobs=-1,
                random_state=seed,
            )
        if model_name == "lgbm":
            if lgb is None:
                raise RuntimeError("lightgbm is not installed")
            return lgb.LGBMClassifier(
                n_estimators=180,
                learning_rate=0.03,
                num_leaves=31,
                subsample=0.9,
                colsample_bytree=0.9,
                class_weight="balanced",
                random_state=seed,
                n_jobs=-1,
                verbose=-1,
            )
    raise ValueError(f"Unsupported model for {target_kind}: {model_name}")


def make_pipeline(target_kind: str, model_name: str, preprocess: str, seed: int, n_classes: int) -> Pipeline:
    steps: list[tuple[str, Any]] = [("imputer", SimpleImputer(strategy="median"))]
    if preprocess == "standard":
        steps.append(("scaler", StandardScaler()))
    elif preprocess == "robust":
        steps.append(("scaler", RobustScaler()))
    elif preprocess == "none":
        pass
    else:
        raise ValueError(f"Unsupported preprocess: {preprocess}")
    steps.append(("model", build_model(target_kind, model_name, seed, n_classes)))
    return Pipeline(steps)


def select_features(
    X_train: np.ndarray,
    y_train: np.ndarray,
    *,
    target_kind: str,
    n_features: int,
    min_presence: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    presence = np.isfinite(X_train).mean(axis=0)
    eligible = presence >= min_presence
    eligible_idx = np.flatnonzero(eligible)
    if len(eligible_idx) == 0:
        raise RuntimeError("No features passed train-fold presence filter")
    X_eligible = X_train[:, eligible_idx]
    if target_kind in {"regression", "binary"}:
        scores = fast_pearson_scores(X_eligible, y_train.astype(float))
    else:
        scores = fast_multiclass_scores(X_eligible, y_train)
    n_keep = min(n_features, len(eligible_idx))
    local = np.argsort(np.abs(scores))[-n_keep:]
    selected = eligible_idx[local]
    selected_scores = scores[local]
    order = np.argsort(np.abs(selected_scores))[::-1]
    return selected[order], selected_scores[order], presence[selected[order]]


def compute_regression_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, Any]:
    r, p = safe_corr(y_true, y_pred)
    return {
        "pearson_r": round(r, 4),
        "pearson_pval": p,
        "mae_weeks": round(float(mean_absolute_error(y_true, y_pred)), 3),
        "rmse_weeks": round(float(math.sqrt(mean_squared_error(y_true, y_pred))), 3),
        "r2": round(float(r2_score(y_true, y_pred)), 4),
    }


def compute_classification_metrics(
    y_true: np.ndarray, y_pred: np.ndarray, y_score: np.ndarray | None, class_names: list[str]
) -> dict[str, Any]:
    metrics: dict[str, Any] = {
        "accuracy": round(float(accuracy_score(y_true, y_pred)), 4),
        "balanced_accuracy": round(float(balanced_accuracy_score(y_true, y_pred)), 4),
        "macro_f1": round(float(f1_score(y_true, y_pred, average="macro", zero_division=0)), 4),
    }
    if y_score is not None and len(class_names) == 2:
        positive_score = y_score[:, 1] if y_score.ndim == 2 else y_score
        try:
            metrics["roc_auc"] = round(float(roc_auc_score(y_true, positive_score)), 4)
        except ValueError:
            metrics["roc_auc"] = None
        try:
            metrics["average_precision"] = round(float(average_precision_score(y_true, positive_score)), 4)
        except ValueError:
            metrics["average_precision"] = None
    elif y_score is not None and len(class_names) > 2:
        try:
            metrics["roc_auc_ovr_macro"] = round(float(roc_auc_score(y_true, y_score, multi_class="ovr", average="macro")), 4)
        except ValueError:
            metrics["roc_auc_ovr_macro"] = None
    return metrics


def run_cv(
    *,
    matrix: pd.DataFrame,
    target: TargetData,
    model_name: str,
    preprocess: str,
    n_features: int,
    min_presence: float,
    seed: int,
    randomize_labels: bool,
) -> tuple[dict[str, Any], pd.DataFrame, pd.DataFrame]:
    sample_ids = target.meta["sample_id"].astype(str).tolist()
    X_raw = matrix[sample_ids].T.to_numpy(dtype=np.float32)
    y = np.asarray(target.y).copy()
    y_labels = np.asarray(target.y_labels).copy()
    if randomize_labels:
        rng = np.random.default_rng(seed)
        perm = rng.permutation(len(y))
        y = y[perm]
        y_labels = y_labels[perm]
    cv_strategy, cv_iter = make_cv(
        TargetData(
            target.target_id,
            target.target_kind,
            target.meta,
            y,
            y_labels,
            target.class_names,
            target.status,
            target.reason,
        ),
        seed,
    )
    if target.target_kind == "regression":
        y_pred = np.full(len(y), np.nan, dtype=float)
        y_score = None
    else:
        y_pred = np.full(len(y), -1, dtype=int)
        y_score = np.full((len(y), len(target.class_names)), np.nan, dtype=float)

    selected_rows: list[dict[str, Any]] = []
    for fold, (train_idx, test_idx) in enumerate(cv_iter, start=1):
        top_idx, top_scores, top_presence = select_features(
            X_raw[train_idx],
            y[train_idx],
            target_kind=target.target_kind,
            n_features=n_features,
            min_presence=min_presence,
        )
        pipe = make_pipeline(target.target_kind, model_name, preprocess, seed + fold, len(target.class_names))
        pipe.fit(X_raw[train_idx][:, top_idx], y[train_idx])
        fold_pred = pipe.predict(X_raw[test_idx][:, top_idx])
        y_pred[test_idx] = fold_pred
        if target.target_kind != "regression" and hasattr(pipe, "predict_proba"):
            prob = pipe.predict_proba(X_raw[test_idx][:, top_idx])
            y_score[test_idx, : prob.shape[1]] = prob
        for rank, (idx, score, presence) in enumerate(zip(top_idx, top_scores, top_presence, strict=True), start=1):
            selected_rows.append(
                {
                    "target_id": target.target_id,
                    "model_name": model_name,
                    "randomize_labels": randomize_labels,
                    "fold": fold,
                    "feature_rank": rank,
                    "feature_id": str(matrix.index[idx]),
                    "selector_score": float(score),
                    "train_presence": float(presence),
                }
            )
    if target.target_kind == "regression":
        metrics = compute_regression_metrics(np.asarray(target.y, dtype=float), y_pred)
        pred_df = target.meta[["sample_id", "dataset_batch", "tissue", "intervention", "age_weeks"]].copy()
        pred_df["target_id"] = target.target_id
        pred_df["model_name"] = model_name
        pred_df["randomize_labels"] = randomize_labels
        pred_df["y_true"] = np.asarray(target.y, dtype=float)
        pred_df["y_pred"] = y_pred
        pred_df["residual"] = pred_df["y_true"] - pred_df["y_pred"]
    else:
        score_arg = y_score if np.isfinite(y_score).all() else None
        metrics = compute_classification_metrics(np.asarray(target.y), y_pred, score_arg, target.class_names)
        pred_df = target.meta[["sample_id", "dataset_batch", "tissue", "intervention", "age_weeks"]].copy()
        pred_df["target_id"] = target.target_id
        pred_df["model_name"] = model_name
        pred_df["randomize_labels"] = randomize_labels
        pred_df["label_true"] = np.asarray(target.y_labels)
        pred_df["label_pred"] = [target.class_names[int(i)] if 0 <= int(i) < len(target.class_names) else "" for i in y_pred]
        if score_arg is not None:
            for idx, cls in enumerate(target.class_names):
                pred_df[f"prob_{cls}"] = score_arg[:, idx]

    metrics.update(
        {
            "target_id": target.target_id,
            "target_kind": target.target_kind,
            "model_name": model_name,
            "preprocess": preprocess,
            "n_feature_prefilter": int(n_features),
            "min_train_feature_presence": float(min_presence),
            "randomize_labels": bool(randomize_labels),
            "cv_strategy": cv_strategy,
            "n_samples": int(len(target.meta)),
            "n_classes": int(len(target.class_names)) if target.class_names else None,
            "classes": ";".join(target.class_names),
        }
    )
    return metrics, pred_df, pd.DataFrame(selected_rows)


def target_audit(target: TargetData) -> dict[str, Any]:
    meta = target.meta
    if target.target_kind == "regression":
        labels = pd.Series(pd.cut(target.y, bins=5, duplicates="drop"), index=meta.index)
    else:
        labels = pd.Series(target.y_labels, index=meta.index)
    row = {
        "target_id": target.target_id,
        "target_kind": target.target_kind,
        "status": target.status,
        "reason": target.reason,
        "n_samples": int(len(meta)),
        "n_datasets": int(meta["dataset_batch"].nunique()) if not meta.empty else 0,
        "n_tissues": int(meta["tissue"].nunique()) if "tissue" in meta.columns and not meta.empty else 0,
        "class_counts": pd.Series(target.y_labels).value_counts().to_dict() if target.target_kind != "regression" else {},
        "dataset_counts": meta["dataset_batch"].value_counts().to_dict() if not meta.empty else {},
        "tissue_counts": meta["tissue"].value_counts().to_dict() if "tissue" in meta.columns and not meta.empty else {},
        "label_dataset_cramers_v": cramers_v(labels, meta["dataset_batch"]) if not meta.empty else 0.0,
        "label_tissue_cramers_v": cramers_v(labels, meta["tissue"]) if "tissue" in meta.columns and not meta.empty else 0.0,
        "label_age_eta2": eta_squared_numeric_by_label(meta.get("age_weeks", pd.Series(dtype=float)), labels)
        if not meta.empty
        else 0.0,
    }
    risk = "low"
    if row["label_dataset_cramers_v"] >= 0.85 or row["label_tissue_cramers_v"] >= 0.85:
        risk = "high_confounding"
    elif row["label_dataset_cramers_v"] >= 0.5 or row["label_tissue_cramers_v"] >= 0.5 or row["label_age_eta2"] >= 0.5:
        risk = "moderate_confounding"
    row["confounding_risk"] = risk
    return row


def feature_audit(matrix: pd.DataFrame, meta: pd.DataFrame, selected: pd.DataFrame, top_n: int) -> pd.DataFrame:
    if selected.empty:
        return pd.DataFrame()
    selected = selected[selected["randomize_labels"].eq(False)].copy()
    if selected.empty:
        return pd.DataFrame()
    grouped = (
        selected.groupby(["target_id", "model_name", "feature_id"], dropna=False)
        .agg(
            selection_count=("feature_id", "count"),
            mean_abs_selector_score=("selector_score", lambda s: float(np.mean(np.abs(s)))),
            mean_selector_score=("selector_score", "mean"),
            mean_feature_rank=("feature_rank", "mean"),
            mean_train_presence=("train_presence", "mean"),
        )
        .reset_index()
        .sort_values(["target_id", "selection_count", "mean_abs_selector_score"], ascending=[True, False, False])
    )
    rows: list[dict[str, Any]] = []
    aligned = [sample for sample in meta["sample_id"].astype(str) if sample in matrix.columns]
    meta = meta[meta["sample_id"].astype(str).isin(aligned)].reset_index(drop=True)
    age = pd.to_numeric(meta["age_weeks"], errors="coerce").to_numpy(dtype=float)
    datasets = meta["dataset_batch"].fillna("unknown").astype(str).to_numpy()
    tissues = meta["tissue"].fillna("unknown").astype(str).to_numpy() if "tissue" in meta.columns else np.array(["unknown"] * len(meta))

    for _, row in grouped.groupby("target_id", group_keys=False).head(top_n).iterrows():
        fid = str(row["feature_id"])
        if fid not in matrix.index:
            continue
        x = matrix.loc[fid, aligned].astype(float).to_numpy()
        age_r, age_p = safe_corr(x, age)
        dataset_means = [float(np.nanmean(x[datasets == dataset])) for dataset in sorted(set(datasets)) if np.isfinite(x[datasets == dataset]).any()]
        tissue_means = [float(np.nanmean(x[tissues == tissue])) for tissue in sorted(set(tissues)) if np.isfinite(x[tissues == tissue]).any()]
        dataset_shift = float(np.nanmax(dataset_means) - np.nanmin(dataset_means)) if dataset_means else 0.0
        tissue_shift = float(np.nanmax(tissue_means) - np.nanmin(tissue_means)) if tissue_means else 0.0
        missing_fraction = float(np.mean(~np.isfinite(x)))
        label = "candidate_signal"
        if missing_fraction > 0.25:
            label = "coverage_driven"
        elif dataset_shift > 0.25:
            label = "dataset_confounded"
        elif tissue_shift > 0.25:
            label = "tissue_confounded"
        parsed = parse_region(fid)
        rows.append(
            {
                **parsed,
                **row.to_dict(),
                "age_pearson_r": round(age_r, 6),
                "age_pearson_p": age_p,
                "missing_fraction": round(missing_fraction, 6),
                "dataset_mean_shift": round(dataset_shift, 6),
                "tissue_mean_shift": round(tissue_shift, 6),
                "confounding_label": label,
                "interpretation_note": "target-associated methylation signal; mechanism not inferred",
            }
        )
    return pd.DataFrame(rows)


def autodl_candidates(registry: pd.DataFrame, gate: pd.DataFrame) -> pd.DataFrame:
    rows = []
    gate_by_id = gate.set_index("target_id").to_dict(orient="index") if not gate.empty else {}
    for _, target in registry.iterrows():
        target_id = str(target["target_id"])
        status = gate_by_id.get(target_id, {}).get("status", "not_evaluated")
        if str(status).startswith("passed"):
            rows.append(
                {
                    "target_id": target_id,
                    "target_kind": target["target_kind"],
                    "architecture": "tabular_res_mlp_or_autoencoder_probe",
                    "n_feature_prefilter": "500;1000;2000",
                    "preprocess": "standard;robust",
                    "local_training_authorized": False,
                    "autodl_only": True,
                    "priority": target.get("priority", ""),
                    "notes": "Use same held-out and random-label gates as local ML; CpGPT/MethFormer only as optional embedding experiment.",
                }
            )
    return pd.DataFrame(rows)


def write_report(
    out_dir: Path,
    *,
    registry: pd.DataFrame,
    gate: pd.DataFrame,
    benchmark: pd.DataFrame,
    random_benchmark: pd.DataFrame,
    audit: pd.DataFrame,
    feature_table: pd.DataFrame,
    summary: dict[str, Any],
) -> None:
    passed = gate[gate["status"].astype(str).str.startswith("passed")] if not gate.empty else pd.DataFrame()
    skipped = gate[~gate["status"].astype(str).str.startswith("passed")] if not gate.empty else pd.DataFrame()
    best_rows = []
    if not benchmark.empty:
        for target_id, part in benchmark.groupby("target_id"):
            metric = str(registry.set_index("target_id").loc[target_id, "primary_metric"])
            sort_col = metric if metric in part.columns else (
                "mae_weeks" if "mae_weeks" in part.columns else "balanced_accuracy"
            )
            ascending = sort_col in {"mae_weeks", "rmse_weeks"}
            best = part.sort_values(sort_col, ascending=ascending).head(1).iloc[0].to_dict()
            best_rows.append(best)
    best = pd.DataFrame(best_rows)

    lines = [
        "# v22 Biological Signal Screen Report",
        "",
        f"Date: {utc_now()}",
        "",
        "## Scope",
        "",
        "This report implements the v22 pivot from age-only clocks to biological-signal discovery using local mouse methylation resources. Age remains the anchor task; CR, sex, tissue, and condition labels are screened as auxiliary biological signals with explicit confounding audits.",
        "",
        "## Inputs",
        "",
        f"- matrix: `{summary['matrix']}`",
        f"- metadata: `{summary['metadata']}`",
        f"- target registry: `{summary['target_registry']}`",
        f"- aligned samples: `{summary['n_aligned_samples']}`",
        f"- regions: `{summary['n_regions']}`",
        "",
        "## Target Gates",
        "",
        f"- passed targets: `{len(passed)}`",
        f"- skipped targets: `{len(skipped)}`",
        "",
    ]
    if not passed.empty:
        lines.extend(["Passed targets:", ""])
        for _, row in passed.iterrows():
            lines.append(f"- `{row['target_id']}`: n={row['n_samples']}, status={row['status']}, reason={row['reason']}")
        lines.append("")
    if not skipped.empty:
        lines.extend(["Skipped targets:", ""])
        for _, row in skipped.iterrows():
            lines.append(f"- `{row['target_id']}`: {row['reason']}")
        lines.append("")
    lines.extend(["## Best Local ML Results", ""])
    if best.empty:
        lines.append("No local ML benchmark rows were produced.")
    else:
        for _, row in best.iterrows():
            bits = [f"`{row['target_id']}`", f"model={row['model_name']}", f"cv={row['cv_strategy']}"]
            for metric in ["mae_weeks", "pearson_r", "balanced_accuracy", "roc_auc", "macro_f1"]:
                if metric in row and pd.notna(row[metric]):
                    bits.append(f"{metric}={row[metric]}")
            lines.append("- " + ", ".join(bits))
    lines.extend(
        [
            "",
            "## Interpretation Guardrails",
            "",
            "Feature tables report target-associated methylation signals only. Regions marked dataset-, tissue-, or coverage-confounded must not be described as mechanisms without follow-up validation.",
            "",
            "## Outputs",
            "",
            f"- `{out_dir / 'target_registry.csv'}`",
            f"- `{out_dir / 'target_gate_table.csv'}`",
            f"- `{out_dir / 'benchmark_summary.csv'}`",
            f"- `{out_dir / 'random_label_sanity_summary.csv'}`",
            f"- `{out_dir / 'target_confounding_audit.csv'}`",
            f"- `{out_dir / 'selected_region_signal_audit.csv'}`",
            f"- `{out_dir / 'autodl_signal_candidate_manifest.csv'}`",
            f"- `{out_dir / 'v22_biological_signal_screen_summary.json'}`",
            "",
            "## External Anchors",
            "",
            "- GSE80672/Petkovich: https://pubmed.ncbi.nlm.nih.gov/28380383/",
            "- GSE121141/Meer whole-lifespan clock: https://elifesciences.org/articles/40675",
            "- Region-based RRBS clock: https://www.pure.ed.ac.uk/ws/portalfiles/portal/347454803/Aging_Cell_2023_Simpson_Region_based_epigenetic_clock_design_improves_RRBS_based_age_prediction.pdf",
            "- scEpiAge code: https://github.com/EpigenomeClock/scEpiAge",
        ]
    )
    report_text = "\n".join(lines) + "\n"
    (out_dir / "v22_biological_signal_screen_report.md").write_text(report_text, encoding="utf-8")
    if out_dir.resolve() == DEFAULT_OUT.resolve():
        DOC_REPORT.parent.mkdir(parents=True, exist_ok=True)
        DOC_REPORT.write_text(report_text, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--target-registry", type=Path, default=DEFAULT_TARGET_REGISTRY)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--targets", default="all", help="Comma-separated target IDs or 'all'.")
    parser.add_argument("--models", default="auto", help="Comma-separated models or 'auto'.")
    parser.add_argument("--n-feature-prefilter", type=int, default=500)
    parser.add_argument("--min-train-feature-presence", type=float, default=0.8)
    parser.add_argument("--preprocess", default="standard", choices=["standard", "robust", "none"])
    parser.add_argument("--random-seed", type=int, default=22)
    parser.add_argument("--include-random-labels", action="store_true")
    parser.add_argument("--smoke", action="store_true", help="Limit large classification targets for quick validation.")
    args = parser.parse_args()

    t0 = time.time()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "predictions").mkdir(exist_ok=True)
    matrix, meta = load_matrix_and_metadata(args.matrix, args.metadata)
    registry = pd.read_csv(args.target_registry)
    if args.targets != "all":
        wanted = {item.strip() for item in args.targets.split(",") if item.strip()}
        registry = registry[registry["target_id"].astype(str).isin(wanted)].copy()

    registry.to_csv(args.out_dir / "target_registry.csv", index=False)
    targets = [build_target(row, meta, smoke=args.smoke) for _, row in registry.iterrows()]
    gate_rows = []
    audit_rows = []
    benchmark_rows = []
    random_rows = []
    prediction_frames = []
    selected_frames = []

    for target in targets:
        audit_rows.append(target_audit(target))
        gate_rows.append(
            {
                "target_id": target.target_id,
                "target_kind": target.target_kind,
                "status": target.status,
                "reason": target.reason,
                "n_samples": int(len(target.meta)),
                "n_classes": int(len(target.class_names)) if target.class_names else None,
                "classes": ";".join(target.class_names),
                "n_datasets": int(target.meta["dataset_batch"].nunique()) if not target.meta.empty else 0,
                "n_tissues": int(target.meta["tissue"].nunique()) if "tissue" in target.meta.columns and not target.meta.empty else 0,
            }
        )
        if not target.status.startswith("passed"):
            continue
        model_names = choose_models(registry, target, args.models)
        for model_name in model_names:
            print(f"[v22] {target.target_id} {model_name} n={len(target.meta)}", flush=True)
            metrics, pred_df, selected = run_cv(
                matrix=matrix,
                target=target,
                model_name=model_name,
                preprocess=args.preprocess,
                n_features=args.n_feature_prefilter,
                min_presence=args.min_train_feature_presence,
                seed=args.random_seed,
                randomize_labels=False,
            )
            benchmark_rows.append(metrics)
            prediction_frames.append(pred_df)
            selected_frames.append(selected)
            pred_df.to_csv(args.out_dir / "predictions" / f"{target.target_id}_{model_name}_predictions.csv", index=False)
            if args.include_random_labels:
                rand_metrics, rand_pred, rand_selected = run_cv(
                    matrix=matrix,
                    target=target,
                    model_name=model_name,
                    preprocess=args.preprocess,
                    n_features=args.n_feature_prefilter,
                    min_presence=args.min_train_feature_presence,
                    seed=args.random_seed + 1009,
                    randomize_labels=True,
                )
                random_rows.append(rand_metrics)
                selected_frames.append(rand_selected)
                rand_pred.to_csv(
                    args.out_dir / "predictions" / f"{target.target_id}_{model_name}_random_label_predictions.csv",
                    index=False,
                )

    gate = pd.DataFrame(gate_rows)
    audit = pd.DataFrame(audit_rows)
    benchmark = pd.DataFrame(benchmark_rows)
    random_benchmark = pd.DataFrame(random_rows)
    selected = pd.concat(selected_frames, ignore_index=True) if selected_frames else pd.DataFrame()
    features = feature_audit(matrix, meta, selected, top_n=100)
    autodl = autodl_candidates(registry, gate)

    gate.to_csv(args.out_dir / "target_gate_table.csv", index=False)
    audit.to_csv(args.out_dir / "target_confounding_audit.csv", index=False)
    benchmark.to_csv(args.out_dir / "benchmark_summary.csv", index=False)
    random_benchmark.to_csv(args.out_dir / "random_label_sanity_summary.csv", index=False)
    selected.to_csv(args.out_dir / "fold_selected_features.csv", index=False)
    features.to_csv(args.out_dir / "selected_region_signal_audit.csv", index=False)
    autodl.to_csv(args.out_dir / "autodl_signal_candidate_manifest.csv", index=False)

    summary = {
        "timestamp": utc_now(),
        "matrix": str(args.matrix),
        "metadata": str(args.metadata),
        "target_registry": str(args.target_registry),
        "out_dir": str(args.out_dir),
        "n_aligned_samples": int(meta.shape[0]),
        "n_regions": int(matrix.shape[0]),
        "n_targets": int(len(targets)),
        "n_targets_passed": int(gate["status"].astype(str).str.startswith("passed").sum()) if not gate.empty else 0,
        "n_benchmark_rows": int(len(benchmark)),
        "n_random_label_rows": int(len(random_benchmark)),
        "n_selected_feature_rows": int(len(selected)),
        "n_interpreted_feature_rows": int(len(features)),
        "exec_time_sec": round(time.time() - t0, 1),
        "local_dl_started": False,
        "autodl_manifest_only": True,
        "interpretation_guardrail": "target-associated methylation signal; mechanism not inferred",
    }
    write_json(args.out_dir / "v22_biological_signal_screen_summary.json", summary)
    write_report(
        args.out_dir,
        registry=registry,
        gate=gate,
        benchmark=benchmark,
        random_benchmark=random_benchmark,
        audit=audit,
        feature_table=features,
        summary=summary,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


def choose_models(registry: pd.DataFrame, target: TargetData, models_arg: str) -> list[str]:
    if models_arg != "auto":
        return [item.strip() for item in models_arg.split(",") if item.strip()]
    row = registry[registry["target_id"].astype(str).eq(target.target_id)].iloc[0]
    recommended = str(row.get("recommended_local_model", "")).strip()
    if recommended:
        return [recommended]
    return ["ridge"] if target.target_kind == "regression" else ["logistic"]


if __name__ == "__main__":
    main()

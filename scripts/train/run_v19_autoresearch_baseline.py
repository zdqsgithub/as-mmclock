#!/usr/bin/env python3
"""v19 classical-ML autoresearch and final baseline fit.

The search is centered on the v18 winner and expands LGBM plus RandomForest
architecture parameters. Candidate scoring uses GroupKFold first, then LODO for
top candidates, with old-age stress penalties and random-label sanity.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import joblib
import numpy as np
import pandas as pd


ROOT = Path("/home/zdq-as/mouse_methyl_work")
DEFAULT_MATRIX = ROOT / "results" / "multidataset_v8_3_ablation" / "all6" / "all_rrbs_region_matrix_5kb.parquet"
DEFAULT_METADATA = ROOT / "metadata" / "model_sample_metadata_v8.csv"
DEFAULT_OUT = ROOT / "results" / "autoresearch_v19_ml_baseline"
REPORT = ROOT / "doc" / "20_analysis" / "50_20260520_v19_autoresearch_ml_baseline_report.md"
PYTHON = "/home/zdq-as/as-ds-ops/.venv/bin/python"
DATASETS = ["GSE120137", "GSE80672", "GSE93957", "GSE121141", "GSE60012", "GSE213628"]


BASE = {
    "model_type": "lgbm",
    "n_feature_prefilter": "1000",
    "imputation": "median",
    "preprocess": "quantile_uniform",
    "min_train_feature_presence": 0.8,
    "min_train_group_feature_presence": 0.8,
    "max_train_dataset_mean_shift": 0.10,
    "min_train_agebin_feature_presence": 0.8,
    "target_transform": "log1p_days",
    "lgbm_n_estimators": 200,
    "lgbm_learning_rate": 0.03,
    "lgbm_num_leaves": 31,
    "lgbm_subsample": 0.9,
    "lgbm_colsample_bytree": 0.9,
    "rf_n_estimators": 150,
    "rf_max_depth": "10",
    "rf_min_samples_leaf": 2,
    "rf_max_features": "1.0",
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str), encoding="utf-8")


def run(cmd: list[str], cwd: Path = ROOT) -> None:
    print(" ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=cwd, check=True)


def config_key(config: dict[str, Any]) -> tuple[tuple[str, Any], ...]:
    return tuple(sorted(config.items()))


def slug_config(config: dict[str, Any], idx: int) -> str:
    parts = [
        f"cfg_{idx:03d}",
        str(config["model_type"]),
        f"top{config['n_feature_prefilter']}",
        str(config["preprocess"]),
        str(config["imputation"]),
        f"s{config['max_train_dataset_mean_shift']}",
        f"p{config['min_train_feature_presence']}",
    ]
    if config["model_type"] == "lgbm":
        parts.extend(
            [
                f"ne{config['lgbm_n_estimators']}",
                f"lr{config['lgbm_learning_rate']}",
                f"leaves{config['lgbm_num_leaves']}",
            ]
        )
    if config["model_type"] == "rf":
        parts.extend(
            [
                f"trees{config['rf_n_estimators']}",
                f"depth{config['rf_max_depth']}",
                f"leaf{config['rf_min_samples_leaf']}",
                f"mf{config['rf_max_features']}",
            ]
        )
    return "_".join(str(part).replace(".", "p").replace(" ", "") for part in parts)


def candidate_configs() -> list[dict[str, Any]]:
    configs: list[dict[str, Any]] = []

    for top in ["500", "750", "1000", "1500", "2000"]:
        for shift in [0.05, 0.075, 0.10, 0.125, 0.15]:
            configs.append({**BASE, "n_feature_prefilter": top, "max_train_dataset_mean_shift": shift})

    for preprocess in ["quantile_uniform", "standard", "robust", "none"]:
        for top in ["500", "1000", "1500"]:
            configs.append({**BASE, "preprocess": preprocess, "n_feature_prefilter": top})

    for presence in [0.7, 0.8, 0.9]:
        for shift in [0.075, 0.10, 0.125]:
            configs.append(
                {
                    **BASE,
                    "min_train_feature_presence": presence,
                    "min_train_group_feature_presence": presence,
                    "min_train_agebin_feature_presence": presence,
                    "max_train_dataset_mean_shift": shift,
                }
            )

    for lgbm_params in [
        {"lgbm_n_estimators": 150, "lgbm_learning_rate": 0.04, "lgbm_num_leaves": 15},
        {"lgbm_n_estimators": 250, "lgbm_learning_rate": 0.025, "lgbm_num_leaves": 31},
        {"lgbm_n_estimators": 350, "lgbm_learning_rate": 0.02, "lgbm_num_leaves": 31},
        {"lgbm_n_estimators": 250, "lgbm_learning_rate": 0.02, "lgbm_num_leaves": 63},
        {"lgbm_n_estimators": 300, "lgbm_learning_rate": 0.015, "lgbm_num_leaves": 15},
        {"lgbm_n_estimators": 200, "lgbm_learning_rate": 0.03, "lgbm_num_leaves": 63},
    ]:
        for top in ["750", "1000", "1500"]:
            configs.append({**BASE, **lgbm_params, "n_feature_prefilter": top})

    rf_base = {**BASE, "model_type": "rf", "preprocess": "quantile_uniform"}
    for top in ["500", "1000", "1500", "2000"]:
        for rf_params in [
            {"rf_n_estimators": 250, "rf_max_depth": "8", "rf_min_samples_leaf": 2, "rf_max_features": "sqrt"},
            {"rf_n_estimators": 250, "rf_max_depth": "12", "rf_min_samples_leaf": 2, "rf_max_features": "sqrt"},
            {"rf_n_estimators": 300, "rf_max_depth": "10", "rf_min_samples_leaf": 4, "rf_max_features": "0.5"},
            {"rf_n_estimators": 350, "rf_max_depth": "14", "rf_min_samples_leaf": 2, "rf_max_features": "0.5"},
            {"rf_n_estimators": 300, "rf_max_depth": "none", "rf_min_samples_leaf": 5, "rf_max_features": "sqrt"},
            {"rf_n_estimators": 250, "rf_max_depth": "10", "rf_min_samples_leaf": 1, "rf_max_features": "1.0"},
        ]:
            configs.append({**rf_base, **rf_params, "n_feature_prefilter": top})

    for preprocess in ["standard", "robust"]:
        for rf_params in [
            {"rf_n_estimators": 300, "rf_max_depth": "10", "rf_min_samples_leaf": 2, "rf_max_features": "sqrt"},
            {"rf_n_estimators": 300, "rf_max_depth": "14", "rf_min_samples_leaf": 4, "rf_max_features": "0.5"},
        ]:
            configs.append({**rf_base, **rf_params, "preprocess": preprocess, "n_feature_prefilter": "1000"})

    deduped = []
    seen = set()
    for config in configs:
        key = config_key(config)
        if key not in seen:
            deduped.append(config)
            seen.add(key)
    return deduped


def model_param_args(config: dict[str, Any]) -> list[str]:
    return [
        "--lgbm_n_estimators",
        str(config["lgbm_n_estimators"]),
        "--lgbm_learning_rate",
        str(config["lgbm_learning_rate"]),
        "--lgbm_num_leaves",
        str(config["lgbm_num_leaves"]),
        "--lgbm_subsample",
        str(config["lgbm_subsample"]),
        "--lgbm_colsample_bytree",
        str(config["lgbm_colsample_bytree"]),
        "--rf_n_estimators",
        str(config["rf_n_estimators"]),
        "--rf_max_depth",
        str(config["rf_max_depth"]),
        "--rf_min_samples_leaf",
        str(config["rf_min_samples_leaf"]),
        "--rf_max_features",
        str(config["rf_max_features"]),
    ]


def train_args(config: dict[str, Any], matrix: Path, metadata: Path, out_dir: Path, randomize: bool = False) -> list[str]:
    cmd = [
        PYTHON,
        "scripts/train/train_clock.py",
        "--matrix_path",
        str(matrix),
        "--metadata_path",
        str(metadata),
        "--feature_type",
        "region",
        "--model_type",
        str(config["model_type"]),
        "--n_feature_prefilter",
        str(config["n_feature_prefilter"]),
        "--imputation",
        str(config["imputation"]),
        "--preprocess",
        str(config["preprocess"]),
        "--min_train_feature_presence",
        str(config["min_train_feature_presence"]),
        "--min_train_group_feature_presence",
        str(config["min_train_group_feature_presence"]),
        "--max_train_dataset_mean_shift",
        str(config["max_train_dataset_mean_shift"]),
        "--min_train_agebin_feature_presence",
        str(config["min_train_agebin_feature_presence"]),
        "--target_transform",
        str(config["target_transform"]),
        *model_param_args(config),
        "--output_dir",
        str(out_dir),
    ]
    if randomize:
        cmd.append("--randomize_labels")
    return cmd


def heldout_args(config: dict[str, Any], matrix: Path, metadata: Path, dataset: str, out_dir: Path) -> list[str]:
    return [
        PYTHON,
        "scripts/train/train_heldout_clock.py",
        "--train_matrix",
        str(matrix),
        "--test_matrix",
        str(matrix),
        "--train_datasets",
        f"all_except:{dataset}",
        "--test_datasets",
        dataset,
        "--metadata_path",
        str(metadata),
        "--model_type",
        str(config["model_type"]),
        "--n_feature_prefilter",
        str(config["n_feature_prefilter"]),
        "--imputation",
        str(config["imputation"]),
        "--preprocess",
        str(config["preprocess"]),
        "--min_train_feature_presence",
        str(config["min_train_feature_presence"]),
        "--min_train_group_feature_presence",
        str(config["min_train_group_feature_presence"]),
        "--max_train_dataset_mean_shift",
        str(config["max_train_dataset_mean_shift"]),
        "--min_train_agebin_feature_presence",
        str(config["min_train_agebin_feature_presence"]),
        "--target_transform",
        str(config["target_transform"]),
        *model_param_args(config),
        "--output_dir",
        str(out_dir),
    ]


def benchmark(predictions: Path, out_path: Path) -> dict[str, Any]:
    if not out_path.exists():
        run([PYTHON, "scripts/validate/benchmark_metrics.py", "--predictions", str(predictions), "--out", str(out_path)])
    return read_json(out_path)


def group_score(metrics: dict[str, Any]) -> float:
    r = float(metrics.get("pearson_r") or 0.0)
    r2 = max(0.0, float(metrics.get("r2") or 0.0))
    mae = float(metrics.get("mae_weeks") or 999.0)
    cr_auc = float(metrics.get("cr_detection_auc") or 0.5)
    mae_component = max(0.0, 1.0 - min(mae, 60.0) / 60.0)
    return round((0.38 * r) + (0.24 * r2) + (0.20 * mae_component) + (0.18 * cr_auc), 6)


def lodo_score(row: dict[str, Any]) -> float:
    mean_lodo_mae = float(row.get("lodo_mean_mae_weeks") or 999.0)
    worst_lodo_mae = float(row.get("lodo_worst_mae_weeks") or 999.0)
    mean_lodo_r = float(row.get("lodo_mean_pearson_r") or 0.0)
    old_mae = row.get("old_104w_weighted_mae_weeks")
    old_penalty = 0.0 if old_mae is None or pd.isna(old_mae) else min(float(old_mae), 100.0) / 100.0
    mean_mae_component = max(0.0, 1.0 - min(mean_lodo_mae, 60.0) / 60.0)
    worst_mae_component = max(0.0, 1.0 - min(worst_lodo_mae, 90.0) / 90.0)
    return round((0.34 * mean_lodo_r) + (0.34 * mean_mae_component) + (0.20 * worst_mae_component) - (0.12 * old_penalty), 6)


def summarize_lodo(lodo_rows: list[dict[str, Any]]) -> dict[str, Any]:
    df = pd.DataFrame(lodo_rows)
    old = df[df["old_104w_n"] > 0].copy()
    old_weighted_mae = None if old.empty else float(np.average(old["old_104w_mae_weeks"], weights=old["old_104w_n"]))
    cr_rows = df[df["cr_detection_auc"].notna()]
    return {
        "lodo_mean_mae_weeks": round(float(df["mae_weeks"].mean()), 3),
        "lodo_worst_mae_weeks": round(float(df["mae_weeks"].max()), 3),
        "lodo_mean_pearson_r": round(float(df["pearson_r"].mean()), 4),
        "lodo_min_pearson_r": round(float(df["pearson_r"].min()), 4),
        "lodo_cr_auc": None if cr_rows.empty else round(float(cr_rows["cr_detection_auc"].mean()), 4),
        "old_104w_weighted_mae_weeks": None if old_weighted_mae is None else round(old_weighted_mae, 3),
    }


def import_train_clock():
    path = ROOT / "scripts" / "train" / "train_clock.py"
    spec = importlib.util.spec_from_file_location("v19_train_clock_helpers", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def parse_rf_max_depth(value: Any) -> int | None:
    return None if str(value).lower() in {"none", "null"} else int(value)


def parse_rf_max_features(value: Any) -> float | str | None:
    text = str(value).lower()
    if text in {"none", "null"}:
        return None
    if text in {"sqrt", "log2"}:
        return text
    return float(value)


def fit_final_all_data_model(config: dict[str, Any], matrix: Path, metadata: Path, out_dir: Path) -> dict[str, Any]:
    helper = import_train_clock()
    beta, meta = helper.load_clock_data(matrix, metadata)
    X_raw = beta.T.values.astype(np.float32)
    y_days = meta["age_days"].values.astype(np.float32)
    y_fit = helper.build_target(y_days, config["target_transform"])
    groups = meta["dataset_batch"].fillna("unknown").astype(str).values

    presence = np.isfinite(X_raw).mean(axis=0) >= float(config["min_train_feature_presence"])
    stability = helper.train_dataset_stability_mask(
        X_raw,
        groups,
        float(config["min_train_group_feature_presence"]),
        float(config["max_train_dataset_mean_shift"]),
    )
    agebin = helper.train_agebin_presence_mask(
        X_raw,
        y_days / 7,
        float(config["min_train_agebin_feature_presence"]),
        8,
    )
    feature_mask = presence & stability & agebin
    presence_idx = np.flatnonzero(feature_mask)
    if len(presence_idx) == 0:
        raise RuntimeError("No features pass final all-data filters.")

    top_local_idx = helper.select_features(X_raw[:, presence_idx], y_fit, config["n_feature_prefilter"])
    top_idx = presence_idx[top_local_idx]
    pipe = helper.Pipeline(
        [
            *helper.make_preprocess_steps(config["preprocess"], config["imputation"], 42, X_raw.shape[0]),
            (
                "model",
                helper.make_model(
                    config["model_type"],
                    42,
                    lgbm_n_estimators=int(config["lgbm_n_estimators"]),
                    lgbm_learning_rate=float(config["lgbm_learning_rate"]),
                    lgbm_num_leaves=int(config["lgbm_num_leaves"]),
                    lgbm_subsample=float(config["lgbm_subsample"]),
                    lgbm_colsample_bytree=float(config["lgbm_colsample_bytree"]),
                    rf_n_estimators=int(config["rf_n_estimators"]),
                    rf_max_depth=parse_rf_max_depth(config["rf_max_depth"]),
                    rf_min_samples_leaf=int(config["rf_min_samples_leaf"]),
                    rf_max_features=parse_rf_max_features(config["rf_max_features"]),
                ),
            ),
        ]
    )
    pipe.fit(X_raw[:, top_idx], y_fit)
    pred_fit = pipe.predict(X_raw[:, top_idx])
    y_pred_days, y_pred_weeks = helper.inverse_target(pred_fit, config["target_transform"], float(y_days.max()))
    metrics = helper.compute_metrics(y_days / 7, y_pred_weeks)

    out_dir.mkdir(parents=True, exist_ok=True)
    joblib.dump(
        {
            "config": config,
            "matrix_path": str(matrix),
            "metadata_path": str(metadata),
            "selected_feature_indices": top_idx,
            "selected_feature_ids": beta.index[top_idx].astype(str).tolist(),
            "pipeline": pipe,
            "warning": "Final all-data fit is for baseline deployment/demo only; use LODO metrics for evidence.",
        },
        out_dir / "final_all_data_ml_baseline_model.joblib",
    )
    pd.DataFrame({"feature_id": beta.index[top_idx].astype(str), "feature_index": top_idx}).to_csv(
        out_dir / "selected_features.csv",
        index=False,
    )
    pd.DataFrame(
        {
            "sample_id": meta["sample_id"].values,
            "dataset_batch": meta.get("dataset_batch", pd.Series(["unknown"] * len(meta))).values,
            "tissue": meta.get("tissue", pd.Series(["unknown"] * len(meta))).values,
            "intervention": meta.get("intervention", pd.Series(["control"] * len(meta))).values,
            "age_days_true": y_days,
            "age_days_pred": y_pred_days,
            "age_weeks_true": y_days / 7,
            "age_weeks_pred": y_pred_weeks,
            "residual_weeks": y_days / 7 - y_pred_weeks,
        }
    ).to_csv(out_dir / "apparent_train_predictions.csv", index=False)
    payload = {
        "warning": "Apparent all-data fit metrics are optimistic and not held-out evidence.",
        "n_samples": int(X_raw.shape[0]),
        "n_features_total": int(X_raw.shape[1]),
        "n_features_after_presence": int(presence.sum()),
        "n_features_after_stability": int(stability.sum()),
        "n_features_after_agebin": int(agebin.sum()),
        "n_features_after_all_filters": int(feature_mask.sum()),
        "n_features_selected": int(len(top_idx)),
        **metrics,
    }
    write_json(out_dir / "apparent_train_metrics.json", payload)
    return payload


def md_table(df: pd.DataFrame) -> str:
    if df.empty:
        return "No rows."
    view = df.copy()
    for col in view.columns:
        view[col] = view[col].map(lambda value: "" if pd.isna(value) else str(value))
    return "\n".join(
        [
            "| " + " | ".join(view.columns) + " |",
            "| " + " | ".join(["---"] * len(view.columns)) + " |",
            *["| " + " | ".join(row) + " |" for row in view.to_numpy(dtype=str)],
        ]
    )


def write_report(
    out_root: Path,
    summary: pd.DataFrame,
    lodo: pd.DataFrame,
    random_metrics: dict[str, Any] | None,
    final_metrics: dict[str, Any],
) -> None:
    best = summary.iloc[0].to_dict()
    lines = [
        "# v19 Autoresearch ML Baseline Report",
        "",
        f"Date: {utc_now()}",
        "",
        "## Scope",
        "",
        "Classical ML autoresearch before deep learning. Search includes LGBM and",
        "RandomForest architecture parameters, leakage-safe GroupKFold, LODO checks",
        "for top candidates, random-label sanity, and final all-data baseline fit.",
        "",
        "## Selected Baseline",
        "",
        f"- config: `{best['config_id']}`",
        f"- model: `{best['model_type']}`",
        f"- top regions: `{best['n_feature_prefilter']}`",
        f"- preprocess/imputation: `{best['preprocess']}` / `{best['imputation']}`",
        f"- filters: presence `{best['min_train_feature_presence']}`, shift `{best['max_train_dataset_mean_shift']}`",
        f"- GroupKFold r/MAE/R2: `{best['pearson_r']}` / `{best['mae_weeks']}` / `{best['r2']}`",
        f"- LODO mean/worst MAE: `{best.get('lodo_mean_mae_weeks')}` / `{best.get('lodo_worst_mae_weeks')}`",
        f"- LODO CR AUC: `{best.get('lodo_cr_auc')}`",
        f"- old104 weighted MAE: `{best.get('old_104w_weighted_mae_weeks')}`",
        f"- final score: `{best.get('final_score')}`",
        "",
        "## Final All-Data Fit",
        "",
        f"- apparent r/MAE/R2: `{final_metrics.get('pearson_r')}` / `{final_metrics.get('mae_weeks')}` / `{final_metrics.get('r2')}`",
        f"- selected features: `{final_metrics.get('n_features_selected')}`",
        "- apparent all-data metrics are optimistic; LODO is the evidence baseline.",
        "",
        "## Random-Label Sanity",
        "",
        "`not_run`" if random_metrics is None else f"`{random_metrics}`",
        "",
        "## Top Ranked Configs",
        "",
        md_table(summary.head(15)),
        "",
        "## LODO Rows",
        "",
        md_table(lodo),
        "",
        "## Outputs",
        "",
        f"- `{out_root.relative_to(ROOT)}/v19_autoresearch_summary.csv`",
        f"- `{out_root.relative_to(ROOT)}/v19_lodo_summary.csv`",
        f"- `{out_root.relative_to(ROOT)}/final_all_data_model/final_all_data_ml_baseline_model.joblib`",
    ]
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    text = "\n".join(lines) + "\n"
    REPORT.write_text(text, encoding="utf-8")
    (out_root / "v19_autoresearch_ml_baseline_report.md").write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", default=str(DEFAULT_MATRIX))
    parser.add_argument("--metadata", default=str(DEFAULT_METADATA))
    parser.add_argument("--out-root", default=str(DEFAULT_OUT))
    parser.add_argument("--max-configs", type=int, default=0, help="0 means all generated candidates.")
    parser.add_argument("--top-lodo", type=int, default=8)
    parser.add_argument("--skip-random-label", action="store_true")
    parser.add_argument("--skip-final-fit", action="store_true")
    args = parser.parse_args()

    matrix = Path(args.matrix)
    metadata = Path(args.metadata)
    out_root = Path(args.out_root)
    out_root.mkdir(parents=True, exist_ok=True)
    started = time.time()

    configs = candidate_configs()
    if args.max_configs > 0:
        configs = configs[: args.max_configs]
    write_json(
        out_root / "v19_run_manifest.json",
        {
            "timestamp": utc_now(),
            "matrix": str(matrix),
            "metadata": str(metadata),
            "n_configs": len(configs),
            "top_lodo": args.top_lodo,
            "candidate_configs": configs,
        },
    )

    rows = []
    id_to_config = {}
    for idx, config in enumerate(configs, start=1):
        config_id = slug_config(config, idx)
        id_to_config[config_id] = config
        cfg_dir = out_root / "groupkfold" / config_id
        cfg_dir.mkdir(parents=True, exist_ok=True)
        write_json(cfg_dir / "config.json", config)
        if not (cfg_dir / "result.json").exists():
            run(train_args(config, matrix, metadata, cfg_dir))
        metrics = benchmark(cfg_dir / "predictions.csv", cfg_dir / "benchmark_result.json")
        rows.append(
            {
                "config_id": config_id,
                **config,
                **{
                    key: metrics.get(key)
                    for key in [
                        "n_samples",
                        "pearson_r",
                        "mae_weeks",
                        "medae_weeks",
                        "rmse_weeks",
                        "r2",
                        "cross_dataset_mae",
                        "cr_detection_auc",
                    ]
                },
                "group_score": group_score(metrics),
            }
        )

    summary = pd.DataFrame(rows).sort_values(["group_score", "pearson_r"], ascending=False).reset_index(drop=True)
    lodo_rows = []
    for config_id in summary.head(args.top_lodo)["config_id"].tolist():
        config = id_to_config[config_id]
        config_lodo_rows = []
        for dataset in DATASETS:
            heldout_dir = out_root / "lodo" / config_id / dataset
            if not (heldout_dir / "result.json").exists():
                run(heldout_args(config, matrix, metadata, dataset, heldout_dir))
            metrics = benchmark(heldout_dir / "predictions.csv", heldout_dir / "benchmark_result.json")
            pred = pd.read_csv(heldout_dir / "predictions.csv")
            old = pred[pred["age_weeks_true"] >= 104]
            row = {
                "config_id": config_id,
                "dataset": dataset,
                **{key: metrics.get(key) for key in ["n_samples", "pearson_r", "mae_weeks", "rmse_weeks", "r2", "cr_detection_auc"]},
                "old_104w_n": int(len(old)),
                "old_104w_mae_weeks": None if old.empty else round(float(np.mean(np.abs(old["age_weeks_true"] - old["age_weeks_pred"]))), 3),
            }
            lodo_rows.append(row)
            config_lodo_rows.append(row)
        lodo_summary = summarize_lodo(config_lodo_rows)
        mask = summary["config_id"].eq(config_id)
        for key, value in lodo_summary.items():
            summary.loc[mask, key] = value
        merged = summary.loc[mask].iloc[0].to_dict()
        summary.loc[mask, "lodo_score"] = lodo_score({**merged, **lodo_summary})

    summary["lodo_score"] = summary["lodo_score"].fillna(-1.0)
    summary["final_score"] = summary["group_score"] + summary["lodo_score"].clip(lower=0.0)
    summary = summary.sort_values(["final_score", "group_score"], ascending=False).reset_index(drop=True)

    random_metrics = None
    best_config_id = str(summary.iloc[0]["config_id"])
    best_config = id_to_config[best_config_id]
    if not args.skip_random_label:
        random_dir = out_root / "random_label_sanity" / best_config_id
        if not (random_dir / "result.json").exists():
            run(train_args(best_config, matrix, metadata, random_dir, randomize=True))
        random_metrics = benchmark(random_dir / "predictions.csv", random_dir / "benchmark_result.json")

    final_metrics: dict[str, Any] = {}
    if not args.skip_final_fit:
        final_metrics = fit_final_all_data_model(best_config, matrix, metadata, out_root / "final_all_data_model")

    lodo_df = pd.DataFrame(lodo_rows)
    summary.to_csv(out_root / "v19_autoresearch_summary.csv", index=False)
    lodo_df.to_csv(out_root / "v19_lodo_summary.csv", index=False)
    payload = {
        "timestamp": utc_now(),
        "exec_time_sec": round(time.time() - started, 1),
        "best_config": summary.iloc[0].to_dict(),
        "random_label_best": random_metrics,
        "final_all_data": final_metrics,
        "summary": summary.to_dict(orient="records"),
        "lodo": lodo_df.to_dict(orient="records"),
    }
    write_json(out_root / "v19_autoresearch_summary.json", payload)
    write_report(out_root, summary, lodo_df, random_metrics, final_metrics)
    print(json.dumps(payload["best_config"], indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Leakage-safe v18 baseline parameter search before deep learning.

This runner extends the v17 locked demo baseline with a small, explicit
classical-ML search. It uses the existing fold-internal training code and
benchmarks top GroupKFold candidates with leave-one-dataset-out validation.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


ROOT = Path("/home/zdq-as/mouse_methyl_work")
DEFAULT_MATRIX = ROOT / "results" / "multidataset_v8_3_ablation" / "all6" / "all_rrbs_region_matrix_5kb.parquet"
DEFAULT_METADATA = ROOT / "metadata" / "model_sample_metadata_v8.csv"
DEFAULT_OUT = ROOT / "results" / "baseline_param_search_v18"
REPORT = ROOT / "doc" / "20_analysis" / "49_20260520_v18_baseline_parameter_search_report.md"
PYTHON = "/home/zdq-as/as-ds-ops/.venv/bin/python"
DATASETS = ["GSE120137", "GSE80672", "GSE93957", "GSE121141", "GSE60012", "GSE213628"]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str), encoding="utf-8")


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def run(cmd: list[str], cwd: Path = ROOT) -> None:
    print(" ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=cwd, check=True)


def slug_config(config: dict[str, Any], idx: int) -> str:
    parts = [
        f"cfg_{idx:03d}",
        str(config["model_type"]),
        f"top{config['n_feature_prefilter']}",
        str(config["preprocess"]),
        str(config["imputation"]),
        f"p{config['min_train_feature_presence']}",
        f"gp{config['min_train_group_feature_presence']}",
        f"s{config['max_train_dataset_mean_shift']}",
        f"ab{config['min_train_agebin_feature_presence']}",
        str(config["target_transform"]),
    ]
    return "_".join(part.replace(".", "p") for part in parts)


def candidate_configs() -> list[dict[str, Any]]:
    locked = {
        "model_type": "lgbm",
        "n_feature_prefilter": "1000",
        "imputation": "median",
        "preprocess": "quantile_uniform",
        "min_train_feature_presence": 0.8,
        "min_train_group_feature_presence": 0.8,
        "max_train_dataset_mean_shift": 0.15,
        "min_train_agebin_feature_presence": 0.8,
        "target_transform": "log1p_days",
    }
    configs: list[dict[str, Any]] = [locked]

    for top in ["500", "2000", "5000"]:
        configs.append({**locked, "n_feature_prefilter": top})
    for preprocess in ["standard", "robust", "none"]:
        configs.append({**locked, "preprocess": preprocess})
    for shift in [0.10, 0.20, 0.30]:
        configs.append({**locked, "max_train_dataset_mean_shift": shift})
    for presence in [0.6, 0.7, 0.9]:
        configs.append(
            {
                **locked,
                "min_train_feature_presence": presence,
                "min_train_group_feature_presence": presence,
                "min_train_agebin_feature_presence": presence,
            }
        )
    for model in ["ridge", "elasticnet", "rf"]:
        configs.append({**locked, "model_type": model})
        configs.append({**locked, "model_type": model, "n_feature_prefilter": "5000"})
    for top in ["500", "1000", "2000", "5000", "10000"]:
        configs.append(
            {
                **locked,
                "model_type": "ridge",
                "n_feature_prefilter": top,
                "preprocess": "standard",
                "imputation": "mean",
                "target_transform": "linear_weeks",
            }
        )
    configs.append({**locked, "target_transform": "linear_weeks"})
    configs.append({**locked, "imputation": "mean"})

    deduped = []
    seen = set()
    for config in configs:
        key = tuple(sorted(config.items()))
        if key not in seen:
            deduped.append(config)
            seen.add(key)
    return deduped


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
    return round((0.40 * r) + (0.25 * r2) + (0.20 * mae_component) + (0.15 * cr_auc), 6)


def lodo_score(row: dict[str, Any]) -> float:
    mean_lodo_mae = float(row.get("lodo_mean_mae_weeks") or 999.0)
    worst_lodo_mae = float(row.get("lodo_worst_mae_weeks") or 999.0)
    mean_lodo_r = float(row.get("lodo_mean_pearson_r") or 0.0)
    old_mae = row.get("old_104w_weighted_mae_weeks")
    old_penalty = 0.0 if old_mae is None or pd.isna(old_mae) else min(float(old_mae), 100.0) / 100.0
    mean_mae_component = max(0.0, 1.0 - min(mean_lodo_mae, 60.0) / 60.0)
    worst_mae_component = max(0.0, 1.0 - min(worst_lodo_mae, 90.0) / 90.0)
    return round((0.35 * mean_lodo_r) + (0.35 * mean_mae_component) + (0.20 * worst_mae_component) - (0.10 * old_penalty), 6)


def summarize_lodo(lodo_rows: list[dict[str, Any]]) -> dict[str, Any]:
    df = pd.DataFrame(lodo_rows)
    old = df[df["old_104w_n"] > 0].copy()
    if old.empty:
        old_weighted_mae = None
    else:
        old_weighted_mae = float(np.average(old["old_104w_mae_weeks"], weights=old["old_104w_n"]))
    return {
        "lodo_mean_mae_weeks": round(float(df["mae_weeks"].mean()), 3),
        "lodo_worst_mae_weeks": round(float(df["mae_weeks"].max()), 3),
        "lodo_mean_pearson_r": round(float(df["pearson_r"].mean()), 4),
        "lodo_min_pearson_r": round(float(df["pearson_r"].min()), 4),
        "old_104w_weighted_mae_weeks": None if old_weighted_mae is None else round(old_weighted_mae, 3),
    }


def write_report(out_root: Path, summary: pd.DataFrame, lodo: pd.DataFrame, random_metrics: dict[str, Any] | None) -> None:
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

    best = summary.iloc[0].to_dict()
    lines = [
        "# v18 Baseline Parameter Search Report",
        "",
        f"Date: {utc_now()}",
        "",
        "## Scope",
        "",
        "Classical ML baseline search before deep learning. Training uses existing",
        "fold-internal feature selection, train-only imputation/scaling, dataset",
        "GroupKFold, and leave-one-dataset-out checks for the top candidates.",
        "",
        "## Best Candidate",
        "",
        f"- config: `{best['config_id']}`",
        f"- model: `{best['model_type']}`",
        f"- top regions: `{best['n_feature_prefilter']}`",
        f"- preprocess/imputation: `{best['preprocess']}` / `{best['imputation']}`",
        f"- filters: presence `{best['min_train_feature_presence']}`, group `{best['min_train_group_feature_presence']}`, shift `{best['max_train_dataset_mean_shift']}`, age-bin `{best['min_train_agebin_feature_presence']}`",
        f"- GroupKFold r/MAE/R2: `{best['pearson_r']}` / `{best['mae_weeks']}` / `{best['r2']}`",
        f"- LODO mean/worst MAE: `{best.get('lodo_mean_mae_weeks')}` / `{best.get('lodo_worst_mae_weeks')}`",
        f"- old 104w weighted MAE: `{best.get('old_104w_weighted_mae_weeks')}`",
        f"- final score: `{best.get('final_score')}`",
        "",
        "## Random-Label Sanity",
        "",
        "`not_run`" if random_metrics is None else f"`{random_metrics}`",
        "",
        "## Top Ranked Configs",
        "",
        md_table(summary.head(10)),
        "",
        "## LODO Rows",
        "",
        md_table(lodo),
        "",
        "## Outputs",
        "",
        f"- `{out_root.relative_to(ROOT)}/v18_baseline_param_search_summary.csv`",
        f"- `{out_root.relative_to(ROOT)}/v18_lodo_summary.csv`",
        f"- `{out_root.relative_to(ROOT)}/v18_baseline_param_search_summary.json`",
    ]
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    text = "\n".join(lines) + "\n"
    REPORT.write_text(text, encoding="utf-8")
    (out_root / "v18_baseline_parameter_search_report.md").write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", default=str(DEFAULT_MATRIX))
    parser.add_argument("--metadata", default=str(DEFAULT_METADATA))
    parser.add_argument("--out-root", default=str(DEFAULT_OUT))
    parser.add_argument("--max-configs", type=int, default=0, help="0 means all explicit configs.")
    parser.add_argument("--top-lodo", type=int, default=5)
    parser.add_argument("--skip-random-label", action="store_true")
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
        out_root / "v18_run_manifest.json",
        {
            "timestamp": utc_now(),
            "matrix": str(matrix),
            "metadata": str(metadata),
            "n_configs": len(configs),
            "top_lodo": args.top_lodo,
            "candidate_configs": configs,
            "guardrails": {
                "download_authorized": False,
                "bismark_authorized": False,
                "deep_learning_authorized": False,
                "feature_selection": "fold_internal_or_train_dataset_only",
            },
        },
    )

    rows = []
    for idx, config in enumerate(configs, start=1):
        config_id = slug_config(config, idx)
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
                **{key: metrics.get(key) for key in ["n_samples", "pearson_r", "mae_weeks", "medae_weeks", "rmse_weeks", "r2", "cross_dataset_mae", "cr_detection_auc"]},
                "group_score": group_score(metrics),
            }
        )

    summary = pd.DataFrame(rows).sort_values(["group_score", "pearson_r"], ascending=False).reset_index(drop=True)
    top_ids = summary.head(args.top_lodo)["config_id"].tolist()
    config_by_id = {slug_config(config, idx): config for idx, config in enumerate(configs, start=1)}
    lodo_rows = []
    lodo_summary_by_id = {}
    for config_id in top_ids:
        config = config_by_id[config_id]
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
        lodo_summary_by_id[config_id] = summarize_lodo(config_lodo_rows)

    for config_id, lodo_summary in lodo_summary_by_id.items():
        mask = summary["config_id"].eq(config_id)
        for key, value in lodo_summary.items():
            summary.loc[mask, key] = value
        base = summary.loc[mask].iloc[0].to_dict()
        summary.loc[mask, "lodo_score"] = lodo_score({**base, **lodo_summary})
    summary["lodo_score"] = summary["lodo_score"].fillna(-1.0)
    summary["final_score"] = summary["group_score"] + summary["lodo_score"].clip(lower=0.0)
    summary = summary.sort_values(["final_score", "group_score"], ascending=False).reset_index(drop=True)

    random_metrics = None
    if not args.skip_random_label and not summary.empty:
        best_id = str(summary.iloc[0]["config_id"])
        random_dir = out_root / "random_label_sanity" / best_id
        if not (random_dir / "result.json").exists():
            run(train_args(config_by_id[best_id], matrix, metadata, random_dir, randomize=True))
        random_metrics = benchmark(random_dir / "predictions.csv", random_dir / "benchmark_result.json")

    lodo_df = pd.DataFrame(lodo_rows)
    summary.to_csv(out_root / "v18_baseline_param_search_summary.csv", index=False)
    lodo_df.to_csv(out_root / "v18_lodo_summary.csv", index=False)
    payload = {
        "timestamp": utc_now(),
        "exec_time_sec": round(time.time() - started, 1),
        "best_config": summary.iloc[0].to_dict() if not summary.empty else None,
        "random_label_best": random_metrics,
        "summary": summary.to_dict(orient="records"),
        "lodo": lodo_df.to_dict(orient="records"),
    }
    write_json(out_root / "v18_baseline_param_search_summary.json", payload)
    write_report(out_root, summary, lodo_df, random_metrics)
    print(json.dumps(payload["best_config"], indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

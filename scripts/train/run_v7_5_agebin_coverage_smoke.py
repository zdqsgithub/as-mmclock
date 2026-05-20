#!/usr/bin/env python3
"""Run v7.5 train-fold age-bin coverage smoke configs."""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("/home/zdq-as/mouse_methyl_work")

CONFIGS = [
    {
        "name": "01_v74_best_no_agebin",
        "model": "lgbm",
        "preprocess": "robust",
        "presence": "0.95",
        "group_presence": "0.95",
        "max_shift": "0.15",
        "agebin_presence": None,
        "top_n": "1000",
    },
    {
        "name": "02_lgbm_robust_p095_top1000_agebin08",
        "model": "lgbm",
        "preprocess": "robust",
        "presence": "0.95",
        "group_presence": "0.95",
        "max_shift": "0.15",
        "agebin_presence": "0.8",
        "top_n": "1000",
    },
    {
        "name": "03_lgbm_robust_p095_top1000_agebin095",
        "model": "lgbm",
        "preprocess": "robust",
        "presence": "0.95",
        "group_presence": "0.95",
        "max_shift": "0.15",
        "agebin_presence": "0.95",
        "top_n": "1000",
    },
    {
        "name": "04_lgbm_robust_p095_top2000_agebin08",
        "model": "lgbm",
        "preprocess": "robust",
        "presence": "0.95",
        "group_presence": "0.95",
        "max_shift": "0.15",
        "agebin_presence": "0.8",
        "top_n": "2000",
    },
    {
        "name": "05_lgbm_quantile_p08_top1000_agebin08",
        "model": "lgbm",
        "preprocess": "quantile_uniform",
        "presence": "0.8",
        "group_presence": "0.8",
        "max_shift": "0.15",
        "agebin_presence": "0.8",
        "top_n": "1000",
    },
]


def run(cmd: list[str]) -> None:
    print(" ".join(cmd))
    subprocess.run(cmd, cwd=ROOT, check=True)


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def per_dataset_metrics(predictions: Path, out_path: Path) -> pd.DataFrame:
    pred = pd.read_csv(predictions)
    rows = []
    for dataset, sub in pred.groupby("dataset_batch", dropna=False):
        true = sub["age_weeks_true"].to_numpy(float)
        fitted = sub["age_weeks_pred"].to_numpy(float)
        if len(sub) >= 3 and np.std(true) > 0 and np.std(fitted) > 0:
            r, pval = stats.pearsonr(true, fitted)
        else:
            r, pval = 0.0, 1.0
        mae = float(np.mean(np.abs(true - fitted)))
        rmse = float(np.sqrt(np.mean((true - fitted) ** 2)))
        ss_res = float(np.sum((true - fitted) ** 2))
        ss_tot = float(np.sum((true - true.mean()) ** 2))
        rows.append(
            {
                "dataset_batch": dataset,
                "n_samples": int(len(sub)),
                "pearson_r": round(float(r), 4),
                "pearson_pval": float(pval),
                "mae_weeks": round(mae, 3),
                "rmse_weeks": round(rmse, 3),
                "r2": round(float(1 - ss_res / ss_tot), 4) if ss_tot > 0 else 0.0,
            }
        )
    metrics = pd.DataFrame(rows).sort_values("dataset_batch")
    metrics.to_csv(out_path, index=False)
    return metrics


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--matrix",
        default=str(ROOT / "results" / "multidataset" / "all_rrbs_region_matrix_5kb.parquet"),
    )
    parser.add_argument(
        "--metadata_path",
        default=str(ROOT / "metadata" / "model_sample_metadata_v7_1.csv"),
    )
    parser.add_argument("--out_root", default=str(ROOT / "results" / "benchmark_v7_5_gse121141_harmonization"))
    args = parser.parse_args()

    out_root = Path(args.out_root)
    out_root.mkdir(parents=True, exist_ok=True)
    rows = []
    for config in CONFIGS:
        out_dir = out_root / config["name"]
        out_dir.mkdir(parents=True, exist_ok=True)
        if (out_dir / "result.json").exists() and (out_dir / "benchmark_result.json").exists():
            print(f"[skip] {config['name']}: existing result.json and benchmark_result.json")
        else:
            cmd = [
                sys.executable,
                "scripts/train/train_clock.py",
                "--matrix_path",
                args.matrix,
                "--metadata_path",
                args.metadata_path,
                "--feature_type",
                "region",
                "--model_type",
                config["model"],
                "--n_feature_prefilter",
                config["top_n"],
                "--imputation",
                "median",
                "--preprocess",
                config["preprocess"],
                "--min_train_feature_presence",
                config["presence"],
                "--min_train_group_feature_presence",
                config["group_presence"],
                "--max_train_dataset_mean_shift",
                config["max_shift"],
                "--target_transform",
                "log1p_days",
                "--output_dir",
                str(out_dir),
            ]
            if config["agebin_presence"] is not None:
                cmd.extend(["--min_train_agebin_feature_presence", config["agebin_presence"]])
            run(cmd)
            run(
                [
                    sys.executable,
                    "scripts/validate/benchmark_metrics.py",
                    "--predictions",
                    str(out_dir / "predictions.csv"),
                    "--out",
                    str(out_dir / "benchmark_result.json"),
                ]
            )
        result = load_json(out_dir / "result.json")
        benchmark = load_json(out_dir / "benchmark_result.json")
        dataset_metrics = per_dataset_metrics(out_dir / "predictions.csv", out_dir / "dataset_metrics.csv")
        gse121141 = dataset_metrics[dataset_metrics["dataset_batch"].eq("GSE121141")]
        rows.append(
            {
                "experiment": config["name"],
                "model_type": config["model"],
                "preprocess": config["preprocess"],
                "min_train_feature_presence": float(config["presence"]),
                "min_train_group_feature_presence": float(config["group_presence"]),
                "max_train_dataset_mean_shift": float(config["max_shift"]),
                "min_train_agebin_feature_presence": (
                    None if config["agebin_presence"] is None else float(config["agebin_presence"])
                ),
                "n_feature_prefilter": config["top_n"],
                "pearson_r": benchmark.get("pearson_r"),
                "mae_weeks": benchmark.get("mae_weeks"),
                "r2": benchmark.get("r2"),
                "cross_dataset_mae": benchmark.get("cross_dataset_mae"),
                "cr_detection_auc": benchmark.get("cr_detection_auc"),
                "gse121141_mae_weeks": None if gse121141.empty else float(gse121141["mae_weeks"].iloc[0]),
                "gse121141_r": None if gse121141.empty else float(gse121141["pearson_r"].iloc[0]),
                "n_features_presence_mean": result.get("n_features_presence_mean"),
                "n_features_stability_mean": result.get("n_features_stability_mean"),
                "n_features_agebin_presence_mean": result.get("n_features_agebin_presence_mean"),
                "n_features_mean": result.get("n_features_mean"),
                "exec_time_sec": result.get("exec_time_sec"),
            }
        )
    summary = pd.DataFrame(rows).sort_values(["gse121141_mae_weeks", "mae_weeks"], ascending=[True, True])
    summary.to_csv(out_root / "v7_5_agebin_coverage_summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()


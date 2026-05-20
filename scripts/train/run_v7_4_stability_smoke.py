#!/usr/bin/env python3
"""Run v7.4 fold-internal coverage/stability-filter smoke configs."""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")

CONFIGS = [
    {
        "name": "01_lgbm_standard_p05_top500_no_stability",
        "model": "lgbm",
        "preprocess": "standard",
        "presence": "0.5",
        "top_n": "500",
        "target": "log1p_days",
        "group_presence": None,
        "max_shift": None,
    },
    {
        "name": "02_lgbm_standard_p05_top500_shift015",
        "model": "lgbm",
        "preprocess": "standard",
        "presence": "0.5",
        "top_n": "500",
        "target": "log1p_days",
        "group_presence": "0.8",
        "max_shift": "0.15",
    },
    {
        "name": "03_lgbm_standard_p05_top500_shift010",
        "model": "lgbm",
        "preprocess": "standard",
        "presence": "0.5",
        "top_n": "500",
        "target": "log1p_days",
        "group_presence": "0.8",
        "max_shift": "0.10",
    },
    {
        "name": "04_lgbm_quantile_p08_top1000_shift015",
        "model": "lgbm",
        "preprocess": "quantile_uniform",
        "presence": "0.8",
        "top_n": "1000",
        "target": "log1p_days",
        "group_presence": "0.8",
        "max_shift": "0.15",
    },
    {
        "name": "05_lgbm_robust_p095_top1000_shift015",
        "model": "lgbm",
        "preprocess": "robust",
        "presence": "0.95",
        "top_n": "1000",
        "target": "log1p_days",
        "group_presence": "0.95",
        "max_shift": "0.15",
    },
]


def run(cmd: list[str]) -> None:
    print(" ".join(cmd))
    subprocess.run(cmd, cwd=ROOT, check=True)


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


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
    parser.add_argument("--out_root", default=str(ROOT / "results" / "benchmark_v7_4_stability"))
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
                "--target_transform",
                config["target"],
                "--output_dir",
                str(out_dir),
            ]
            if config["group_presence"] is not None:
                cmd.extend(["--min_train_group_feature_presence", config["group_presence"]])
            if config["max_shift"] is not None:
                cmd.extend(["--max_train_dataset_mean_shift", config["max_shift"]])
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
        rows.append(
            {
                "experiment": config["name"],
                "model_type": config["model"],
                "preprocess": config["preprocess"],
                "min_train_feature_presence": float(config["presence"]),
                "min_train_group_feature_presence": (
                    None if config["group_presence"] is None else float(config["group_presence"])
                ),
                "max_train_dataset_mean_shift": None if config["max_shift"] is None else float(config["max_shift"]),
                "n_feature_prefilter": config["top_n"],
                "target_transform": config["target"],
                "pearson_r": benchmark.get("pearson_r"),
                "mae_weeks": benchmark.get("mae_weeks"),
                "r2": benchmark.get("r2"),
                "cross_dataset_mae": benchmark.get("cross_dataset_mae"),
                "cr_detection_auc": benchmark.get("cr_detection_auc"),
                "n_features_presence_mean": result.get("n_features_presence_mean"),
                "n_features_stability_mean": result.get("n_features_stability_mean"),
                "n_features_mean": result.get("n_features_mean"),
                "exec_time_sec": result.get("exec_time_sec"),
            }
        )
    summary = pd.DataFrame(rows).sort_values(["mae_weeks", "pearson_r"], ascending=[True, False])
    summary.to_csv(out_root / "v7_4_stability_summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Run the fixed v7.3 five-dataset outer-join benchmark smoke grid."""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")

CONFIGS = [
    ("01_lgbm_standard_p05_top500_log1p", "lgbm", "standard", "0.5", "500", "log1p_days"),
    ("02_lgbm_quantile_p08_top1000_log1p", "lgbm", "quantile_uniform", "0.8", "1000", "log1p_days"),
    ("03_lgbm_robust_p095_top1000_log1p", "lgbm", "robust", "0.95", "1000", "log1p_days"),
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
        default=str(ROOT / "results" / "multidataset_v7_3_outer" / "all_rrbs_region_matrix_5kb.parquet"),
    )
    parser.add_argument(
        "--metadata_path",
        default=str(ROOT / "metadata" / "model_sample_metadata_v7_1.csv"),
    )
    parser.add_argument("--out_root", default=str(ROOT / "results" / "benchmark_v7_3_outer_join"))
    args = parser.parse_args()

    out_root = Path(args.out_root)
    out_root.mkdir(parents=True, exist_ok=True)
    rows = []
    for name, model, preprocess, presence, top_n, target_transform in CONFIGS:
        out_dir = out_root / name
        out_dir.mkdir(parents=True, exist_ok=True)
        if (out_dir / "result.json").exists() and (out_dir / "benchmark_result.json").exists():
            print(f"[skip] {name}: existing result.json and benchmark_result.json")
        else:
            run(
                [
                    sys.executable,
                    "scripts/train/train_clock.py",
                    "--matrix_path",
                    args.matrix,
                    "--metadata_path",
                    args.metadata_path,
                    "--feature_type",
                    "region",
                    "--model_type",
                    model,
                    "--n_feature_prefilter",
                    top_n,
                    "--imputation",
                    "median",
                    "--preprocess",
                    preprocess,
                    "--min_train_feature_presence",
                    presence,
                    "--target_transform",
                    target_transform,
                    "--output_dir",
                    str(out_dir),
                ]
            )
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
                "experiment": name,
                "model_type": model,
                "preprocess": preprocess,
                "min_train_feature_presence": float(presence),
                "n_feature_prefilter": top_n,
                "target_transform": target_transform,
                "pearson_r": benchmark.get("pearson_r"),
                "mae_weeks": benchmark.get("mae_weeks"),
                "r2": benchmark.get("r2"),
                "cross_dataset_mae": benchmark.get("cross_dataset_mae"),
                "cr_detection_auc": benchmark.get("cr_detection_auc"),
                "n_features_presence_mean": result.get("n_features_presence_mean"),
                "n_features_mean": result.get("n_features_mean"),
                "exec_time_sec": result.get("exec_time_sec"),
            }
        )
    summary = pd.DataFrame(rows).sort_values(["mae_weeks", "pearson_r"], ascending=[True, False])
    summary.to_csv(out_root / "v7_3_outer_join_summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Run v7.2 GSE60012 auxiliary held-out validation."""
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
        "--train_matrix",
        default=str(ROOT / "results" / "multidataset" / "all_rrbs_region_matrix_5kb.parquet"),
    )
    parser.add_argument(
        "--test_matrix",
        default=str(ROOT / "results" / "multidataset" / "GSE60012_region_matrix_5kb.parquet"),
    )
    parser.add_argument(
        "--metadata_path",
        default=str(ROOT / "metadata" / "model_sample_metadata_v7_1.csv"),
    )
    parser.add_argument("--out_root", default=str(ROOT / "results" / "validation_v7_2_gse60012_aux_heldout"))
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
                    "scripts/train/train_heldout_clock.py",
                    "--train_matrix",
                    args.train_matrix,
                    "--test_matrix",
                    args.test_matrix,
                    "--metadata_path",
                    args.metadata_path,
                    "--train_datasets",
                    "all_except:GSE60012",
                    "--test_datasets",
                    "GSE60012",
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
                "n_train_samples": result.get("n_train_samples"),
                "n_test_samples": result.get("n_test_samples"),
                "n_common_features": result.get("n_common_features"),
                "n_features_passing_train_presence": result.get("n_features_passing_train_presence"),
                "exec_time_sec": result.get("exec_time_sec"),
            }
        )
    summary = pd.DataFrame(rows).sort_values(["mae_weeks", "pearson_r"], ascending=[True, False])
    summary.to_csv(out_root / "v7_2_gse60012_aux_summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()

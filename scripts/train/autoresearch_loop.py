#!/usr/bin/env python3
"""
Autoresearch loop for leakage-safe Phase 0 clock experiments.

The fixed search space contains 48 unique configurations. When asked for 100
attempts, this script runs all unique configurations once and records the
remaining attempt slots as skipped_exhausted instead of duplicating runs.
"""
from __future__ import annotations

import argparse
import itertools
import json
import random
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path("/home/zdq-as/mouse_methyl_work")
PYTHON = "/home/zdq-as/as-ds-ops/.venv/bin/python"
TRAIN_SCRIPT = ROOT / "scripts" / "train" / "train_clock.py"
BENCHMARK_SCRIPT = ROOT / "scripts" / "validate" / "benchmark_metrics.py"

MODEL_TYPES = ["ridge", "elasticnet", "lgbm", "rf"]
N_CPG_PREFILTER = [1000, 2000, 5000, 10000, 20000, 50000]
N_REGION_PREFILTER = ["500", "1000", "2000", "5000", "10000", "all"]
IMPUTATIONS = ["median", "mean"]


def all_configs(seed: int, feature_type: str) -> list[dict]:
    feature_counts = N_REGION_PREFILTER if feature_type == "region" else [str(n) for n in N_CPG_PREFILTER]
    configs = [
        {
            "model_type": model_type,
            "n_feature_prefilter": n_feature,
            "imputation": imputation,
            "feature_selector": "age_correlation",
            "feature_type": feature_type,
        }
        for model_type, n_feature, imputation in itertools.product(MODEL_TYPES, feature_counts, IMPUTATIONS)
    ]
    rng = random.Random(seed)
    rng.shuffle(configs)
    return configs


def composite_score(metrics: dict) -> float:
    pearson = float(metrics.get("pearson_r") or 0.0)
    r2 = max(0.0, float(metrics.get("r2") or 0.0))
    mae = max(float(metrics.get("mae_weeks") or 1e9), 1e-9)
    cr_auc = metrics.get("cr_detection_auc")
    rapa_auc = metrics.get("rapamycin_detection_auc")
    intervention_auc = cr_auc if cr_auc is not None else rapa_auc
    if intervention_auc is not None:
        return (0.45 * pearson) + (0.20 * r2) + (0.25 * float(intervention_auc)) + (0.10 / mae)
    return (0.60 * pearson) + (0.30 * r2) + (0.10 / mae)


def write_header(path: Path) -> None:
    path.write_text(
        "\t".join(
            [
                "exp_id",
                "status",
                "feature_type",
                "model_type",
                "n_feature_prefilter",
                "imputation",
                "pearson_r",
                "mae_weeks",
                "medae_weeks",
                "rmse_weeks",
                "r2",
                "cr_detection_auc",
                "rapamycin_detection_auc",
                "composite_score",
                "exec_time_sec",
                "error",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def append_row(path: Path, row: dict) -> None:
    columns = [
        "exp_id",
        "status",
        "feature_type",
        "model_type",
        "n_feature_prefilter",
        "imputation",
        "pearson_r",
        "mae_weeks",
        "medae_weeks",
        "rmse_weeks",
        "r2",
        "cr_detection_auc",
        "rapamycin_detection_auc",
        "composite_score",
        "exec_time_sec",
        "error",
    ]
    with path.open("a", encoding="utf-8") as handle:
        handle.write("\t".join(str(row.get(col, "")) for col in columns) + "\n")


def run_experiment(exp_id: str, config: dict, results_dir: Path, matrix_path: Path) -> tuple[dict, dict]:
    exp_dir = results_dir / exp_id
    exp_dir.mkdir(parents=True, exist_ok=True)
    train_cmd = [
        PYTHON,
        str(TRAIN_SCRIPT),
        "--model_type",
        config["model_type"],
        "--n_feature_prefilter",
        str(config["n_feature_prefilter"]),
        "--feature_type",
        config["feature_type"],
        "--matrix_path",
        str(matrix_path),
        "--imputation",
        config["imputation"],
        "--output_dir",
        str(exp_dir),
    ]
    subprocess.run(train_cmd, check=True)

    benchmark_cmd = [
        PYTHON,
        str(BENCHMARK_SCRIPT),
        "--predictions",
        str(exp_dir / "predictions.csv"),
        "--out",
        str(exp_dir / "benchmark_result.json"),
    ]
    subprocess.run(benchmark_cmd, check=True)

    result = json.loads((exp_dir / "result.json").read_text(encoding="utf-8"))
    benchmark = json.loads((exp_dir / "benchmark_result.json").read_text(encoding="utf-8"))
    return result, benchmark


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n_attempts", type=int, default=100)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--feature_type", default="site", choices=["site", "region"])
    parser.add_argument("--results_dir", default=None)
    parser.add_argument("--matrix_path", default=None)
    args = parser.parse_args()

    if args.feature_type == "region":
        default_results_dir = ROOT / "results" / "autoresearch_v3_region"
        default_matrix_path = ROOT / "results" / "phase0" / "region_matrix_5kb.parquet"
        loop_label = "v3 region"
    else:
        default_results_dir = ROOT / "results" / "autoresearch_v2"
        default_matrix_path = ROOT / "results" / "phase0" / "beta_matrix_thompson.parquet"
        loop_label = "v2 site"
    results_dir = Path(args.results_dir) if args.results_dir else default_results_dir
    matrix_path = Path(args.matrix_path) if args.matrix_path else default_matrix_path

    results_dir.mkdir(parents=True, exist_ok=True)
    results_log = results_dir / "autoresearch_summary.tsv"
    if not args.resume or not results_log.exists():
        write_header(results_log)

    if not TRAIN_SCRIPT.exists() or not BENCHMARK_SCRIPT.exists():
        print("[ERROR] Missing train or benchmark script.")
        sys.exit(1)

    configs = all_configs(args.seed, args.feature_type)
    n_to_run = min(args.n_attempts, len(configs))
    best_score = -float("inf")

    print("=" * 60)
    print(f"AS-DS-Ops: Epigenetic Clock Autoresearch Loop {loop_label}")
    print(f"Unique configs available: {len(configs)} | attempt slots requested: {args.n_attempts}")
    print(f"Feature matrix: {matrix_path}")
    print("=" * 60)

    for i, config in enumerate(configs[:n_to_run], start=1):
        exp_id = (
            f"exp_{i:03d}_{config['model_type']}_"
            f"{config['n_feature_prefilter']}_{config['imputation']}"
        )
        print(f"\n[{i}/{args.n_attempts}] Launching {exp_id}: {config}")
        started = time.time()
        try:
            result, benchmark = run_experiment(exp_id, config, results_dir, matrix_path)
            score = composite_score(benchmark)
            elapsed = round(time.time() - started, 1)
            row = {
                "exp_id": exp_id,
                "status": "completed",
                "feature_type": config["feature_type"],
                "model_type": config["model_type"],
                "n_feature_prefilter": config["n_feature_prefilter"],
                "imputation": config["imputation"],
                "pearson_r": benchmark.get("pearson_r"),
                "mae_weeks": benchmark.get("mae_weeks"),
                "medae_weeks": benchmark.get("medae_weeks"),
                "rmse_weeks": benchmark.get("rmse_weeks"),
                "r2": benchmark.get("r2"),
                "cr_detection_auc": benchmark.get("cr_detection_auc"),
                "rapamycin_detection_auc": benchmark.get("rapamycin_detection_auc"),
                "composite_score": round(score, 6),
                "exec_time_sec": result.get("exec_time_sec", elapsed),
                "error": "",
            }
            append_row(results_log, row)
            print(
                f"  [Result] R:{row['pearson_r']} MAE:{row['mae_weeks']} "
                f"R2:{row['r2']} Comp:{row['composite_score']}"
            )
            if score > best_score:
                best_score = score
                print(f"  New best score: {best_score:.6f}")
        except subprocess.CalledProcessError as exc:
            append_row(
                results_log,
                {
                    "exp_id": exp_id,
                    "status": "failed",
                    "feature_type": config["feature_type"],
                    "model_type": config["model_type"],
                    "n_feature_prefilter": config["n_feature_prefilter"],
                    "imputation": config["imputation"],
                    "error": str(exc)[:500],
                },
            )
            print(f"  [ERROR] Experiment failed: {exc}")

    for slot in range(n_to_run + 1, args.n_attempts + 1):
        append_row(
            results_log,
            {
                "exp_id": f"slot_{slot:03d}",
                "status": "skipped_exhausted",
                "feature_type": args.feature_type,
                "error": f"Fixed {loop_label} search space has only {len(configs)} unique configurations.",
            },
        )

    print("\n" + "=" * 60)
    print(f"Autoresearch {loop_label} complete. Results logged to {results_log}")


if __name__ == "__main__":
    main()

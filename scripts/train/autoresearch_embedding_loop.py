#!/usr/bin/env python3
"""
Autoresearch loop for embedding-aware Phase 0 clock experiments.
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
TRAIN_SCRIPT = ROOT / "scripts" / "train" / "train_embedding_clock.py"
BENCHMARK_SCRIPT = ROOT / "scripts" / "validate" / "benchmark_metrics.py"

EMBEDDERS = ["svd", "nmf"]
COMPONENTS = [16, 32, 64, 128]
MODEL_TYPES = ["ridge", "elasticnet", "lgbm"]
HYBRID_TOP_REGIONS = [0, 500, 1000]
IMPUTATIONS = ["median", "mean"]


def parse_csv_values(value: str | None, cast=str):
    if value is None:
        return None
    return [cast(item.strip()) for item in value.split(",") if item.strip()]


def all_configs(
    seed: int,
    embedders: list[str] | None = None,
    components: list[int] | None = None,
    model_types: list[str] | None = None,
    hybrid_top_regions: list[int] | None = None,
    imputations: list[str] | None = None,
) -> list[dict]:
    embedders = embedders or EMBEDDERS
    components = components or COMPONENTS
    model_types = model_types or MODEL_TYPES
    hybrid_top_regions = hybrid_top_regions or HYBRID_TOP_REGIONS
    imputations = imputations or IMPUTATIONS
    configs = [
        {
            "embedder": embedder,
            "n_components": n_components,
            "model_type": model_type,
            "hybrid_top_regions": hybrid_top_regions,
            "imputation": imputation,
        }
        for embedder, n_components, model_type, hybrid_top_regions, imputation in itertools.product(
            embedders, components, model_types, hybrid_top_regions, imputations
        )
    ]
    rng = random.Random(seed)
    rng.shuffle(configs)
    return configs


def composite_score(metrics: dict) -> float:
    pearson = float(metrics.get("pearson_r") or 0.0)
    r2 = max(0.0, float(metrics.get("r2") or 0.0))
    mae = max(float(metrics.get("mae_weeks") or 1e9), 1e-9)
    return (0.60 * pearson) + (0.30 * r2) + (0.10 / mae)


def write_header(path: Path) -> None:
    path.write_text(
        "\t".join(
            [
                "exp_id",
                "status",
                "embedder",
                "n_components",
                "model_type",
                "hybrid_top_regions",
                "imputation",
                "pearson_r",
                "mae_weeks",
                "medae_weeks",
                "rmse_weeks",
                "r2",
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
        "embedder",
        "n_components",
        "model_type",
        "hybrid_top_regions",
        "imputation",
        "pearson_r",
        "mae_weeks",
        "medae_weeks",
        "rmse_weeks",
        "r2",
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
        "--matrix_path",
        str(matrix_path),
        "--embedder",
        config["embedder"],
        "--n_components",
        str(config["n_components"]),
        "--model_type",
        config["model_type"],
        "--hybrid_top_regions",
        str(config["hybrid_top_regions"]),
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
    parser.add_argument("--matrix_path", default=str(ROOT / "results" / "phase0" / "region_matrix_5kb.parquet"))
    parser.add_argument("--results_dir", default=str(ROOT / "results" / "autoresearch_v4_embedding"))
    parser.add_argument("--embedders", default=None, help="Optional comma-separated subset, e.g. svd or svd,nmf.")
    parser.add_argument("--components", default=None, help="Optional comma-separated subset, e.g. 32,64.")
    parser.add_argument("--model_types", default=None, help="Optional comma-separated subset, e.g. ridge,lgbm.")
    parser.add_argument("--hybrid_top_regions", default=None, help="Optional comma-separated subset, e.g. 0,500.")
    parser.add_argument("--imputations", default=None, help="Optional comma-separated subset, e.g. median.")
    args = parser.parse_args()

    if not TRAIN_SCRIPT.exists() or not BENCHMARK_SCRIPT.exists():
        print("[ERROR] Missing train or benchmark script.")
        sys.exit(1)

    results_dir = Path(args.results_dir)
    matrix_path = Path(args.matrix_path)
    results_dir.mkdir(parents=True, exist_ok=True)
    results_log = results_dir / "autoresearch_summary.tsv"
    if not args.resume or not results_log.exists():
        write_header(results_log)

    configs = all_configs(
        args.seed,
        embedders=parse_csv_values(args.embedders, str),
        components=parse_csv_values(args.components, int),
        model_types=parse_csv_values(args.model_types, str),
        hybrid_top_regions=parse_csv_values(args.hybrid_top_regions, int),
        imputations=parse_csv_values(args.imputations, str),
    )
    n_to_run = min(args.n_attempts, len(configs))
    best_score = -float("inf")

    print("=" * 60)
    print("AS-DS-Ops: Epigenetic Clock Autoresearch Loop v4 embedding")
    print(f"Unique configs available: {len(configs)} | attempt slots requested: {args.n_attempts}")
    print(f"Feature matrix: {matrix_path}")
    print("=" * 60)

    for i, config in enumerate(configs[:n_to_run], start=1):
        hybrid = f"hybrid{config['hybrid_top_regions']}"
        exp_id = (
            f"exp_{i:03d}_{config['embedder']}{config['n_components']}_"
            f"{config['model_type']}_{hybrid}_{config['imputation']}"
        )
        print(f"\n[{i}/{args.n_attempts}] Launching {exp_id}: {config}")
        started = time.time()
        try:
            result, benchmark = run_experiment(exp_id, config, results_dir, matrix_path)
            score = composite_score(benchmark)
            row = {
                "exp_id": exp_id,
                "status": "completed",
                "embedder": config["embedder"],
                "n_components": config["n_components"],
                "model_type": config["model_type"],
                "hybrid_top_regions": config["hybrid_top_regions"],
                "imputation": config["imputation"],
                "pearson_r": benchmark.get("pearson_r"),
                "mae_weeks": benchmark.get("mae_weeks"),
                "medae_weeks": benchmark.get("medae_weeks"),
                "rmse_weeks": benchmark.get("rmse_weeks"),
                "r2": benchmark.get("r2"),
                "composite_score": round(score, 6),
                "exec_time_sec": result.get("exec_time_sec", round(time.time() - started, 1)),
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
                    "embedder": config["embedder"],
                    "n_components": config["n_components"],
                    "model_type": config["model_type"],
                    "hybrid_top_regions": config["hybrid_top_regions"],
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
                "error": f"Fixed v4 search space has only {len(configs)} unique configurations.",
            },
        )

    print("\n" + "=" * 60)
    print(f"Autoresearch v4 embedding complete. Results logged to {results_log}")


if __name__ == "__main__":
    main()

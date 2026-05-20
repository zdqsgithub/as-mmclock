#!/usr/bin/env python3
"""Run the v8 GSE213628 RALPH smoke benchmark.

The runner intentionally evaluates only the locked v7.5/v7.4 configurations.
It decides whether a later constrained autoresearch search is justified.
"""
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
        "name": "01_v75_lgbm_quantile_p08_top1000_agebin08",
        "model": "lgbm",
        "preprocess": "quantile_uniform",
        "presence": "0.8",
        "group_presence": "0.8",
        "max_shift": "0.15",
        "agebin_presence": "0.8",
        "top_n": "1000",
    },
    {
        "name": "02_v74_lgbm_robust_p095_top1000_shift015",
        "model": "lgbm",
        "preprocess": "robust",
        "presence": "0.95",
        "group_presence": "0.95",
        "max_shift": "0.15",
        "agebin_presence": None,
        "top_n": "1000",
    },
]


def run(cmd: list[str]) -> None:
    print(" ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=ROOT, check=True)


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def safe_r(true: np.ndarray, pred: np.ndarray) -> tuple[float, float]:
    if len(true) < 3 or np.std(true) == 0 or np.std(pred) == 0:
        return 0.0, 1.0
    r, pval = stats.pearsonr(true, pred)
    return float(r), float(pval)


def metric_row(dataset: str, sub: pd.DataFrame) -> dict:
    true = sub["age_weeks_true"].to_numpy(float)
    pred = sub["age_weeks_pred"].to_numpy(float)
    r, pval = safe_r(true, pred)
    ss_res = float(np.sum((true - pred) ** 2))
    ss_tot = float(np.sum((true - true.mean()) ** 2))
    return {
        "dataset_batch": dataset,
        "n_samples": int(len(sub)),
        "pearson_r": round(r, 4),
        "pearson_pval": pval,
        "mae_weeks": round(float(np.mean(np.abs(true - pred))), 3),
        "rmse_weeks": round(float(np.sqrt(np.mean((true - pred) ** 2))), 3),
        "r2": round(float(1 - ss_res / ss_tot), 4) if ss_tot > 0 else 0.0,
        "old_104w_n": int((sub["age_weeks_true"] >= 104).sum()),
        "old_104w_mae_weeks": (
            None
            if (sub["age_weeks_true"] >= 104).sum() == 0
            else round(
                float(
                    np.mean(
                        np.abs(
                            sub.loc[sub["age_weeks_true"] >= 104, "age_weeks_true"]
                            - sub.loc[sub["age_weeks_true"] >= 104, "age_weeks_pred"]
                        )
                    )
                ),
                3,
            )
        ),
    }


def per_dataset_metrics(predictions: Path, out_path: Path) -> pd.DataFrame:
    pred = pd.read_csv(predictions)
    rows = [metric_row(dataset, sub) for dataset, sub in pred.groupby("dataset_batch", dropna=False)]
    metrics = pd.DataFrame(rows).sort_values("dataset_batch")
    metrics.to_csv(out_path, index=False)
    return metrics


def train_args(config: dict, matrix: str, metadata_path: str, out_dir: Path, randomize: bool = False) -> list[str]:
    cmd = [
        sys.executable,
        "scripts/train/train_clock.py",
        "--matrix_path",
        matrix,
        "--metadata_path",
        metadata_path,
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
    if randomize:
        cmd.append("--randomize_labels")
    return cmd


def heldout_args(config: dict, matrix: str, metadata_path: str, test_dataset: str, out_dir: Path) -> list[str]:
    cmd = [
        sys.executable,
        "scripts/train/train_heldout_clock.py",
        "--train_matrix",
        matrix,
        "--test_matrix",
        matrix,
        "--train_datasets",
        f"all_except:{test_dataset}",
        "--test_datasets",
        test_dataset,
        "--metadata_path",
        metadata_path,
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
    return cmd


def run_benchmark(predictions: Path, out: Path) -> dict:
    run([sys.executable, "scripts/validate/benchmark_metrics.py", "--predictions", str(predictions), "--out", str(out)])
    return load_json(out)


def shuffle_intervention(predictions: Path, out_csv: Path, seed: int = 42) -> None:
    pred = pd.read_csv(predictions)
    rng = np.random.default_rng(seed)
    shuffled = pred["intervention"].to_numpy(str).copy()
    rng.shuffle(shuffled)
    pred["intervention"] = shuffled
    pred.to_csv(out_csv, index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--matrix",
        default=str(ROOT / "results" / "multidataset_v8_gse213628" / "all_rrbs_region_matrix_5kb.parquet"),
    )
    parser.add_argument("--metadata_path", default=str(ROOT / "metadata" / "model_sample_metadata_v8.csv"))
    parser.add_argument("--out_root", default=str(ROOT / "results" / "benchmark_v8_gse213628"))
    parser.add_argument("--validation_root", default=str(ROOT / "results"))
    args = parser.parse_args()

    out_root = Path(args.out_root)
    out_root.mkdir(parents=True, exist_ok=True)
    rows = []
    for config in CONFIGS:
        out_dir = out_root / config["name"]
        out_dir.mkdir(parents=True, exist_ok=True)
        if not (out_dir / "result.json").exists():
            run(train_args(config, args.matrix, args.metadata_path, out_dir))
        if not (out_dir / "benchmark_result.json").exists():
            run_benchmark(out_dir / "predictions.csv", out_dir / "benchmark_result.json")
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
                "gse121141_old_104w_mae_weeks": (
                    None if gse121141.empty else gse121141["old_104w_mae_weeks"].iloc[0]
                ),
                "n_features_mean": result.get("n_features_mean"),
                "exec_time_sec": result.get("exec_time_sec"),
            }
        )
    summary = pd.DataFrame(rows).sort_values(
        ["gse121141_old_104w_mae_weeks", "gse121141_mae_weeks", "mae_weeks"],
        ascending=[True, True, True],
    )
    summary_path = out_root / "v8_gse213628_ralph_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    best_name = str(summary.iloc[0]["experiment"])
    best_config = next(config for config in CONFIGS if config["name"] == best_name)
    print(summary.to_string(index=False), flush=True)
    print(f"[v8] best_config={best_name}", flush=True)

    random_dir = out_root / f"{best_name}_random_labels"
    if not (random_dir / "result.json").exists():
        run(train_args(best_config, args.matrix, args.metadata_path, random_dir, randomize=True))
    random_benchmark = run_benchmark(random_dir / "predictions.csv", random_dir / "benchmark_result.json")

    validation_root = Path(args.validation_root)
    heldout_results = {}
    for test_dataset, out_name in [
        ("GSE121141", "validation_v8_gse121141_all_except"),
        ("GSE80672", "validation_v8_gse80672_cr_all_except"),
    ]:
        heldout_dir = validation_root / out_name / best_name
        heldout_dir.mkdir(parents=True, exist_ok=True)
        if not (heldout_dir / "result.json").exists():
            run(heldout_args(best_config, args.matrix, args.metadata_path, test_dataset, heldout_dir))
        heldout_benchmark = run_benchmark(heldout_dir / "predictions.csv", heldout_dir / "benchmark_result.json")
        per_dataset_metrics(heldout_dir / "predictions.csv", heldout_dir / "dataset_metrics.csv")
        heldout_results[test_dataset] = {
            "dir": str(heldout_dir),
            "result": load_json(heldout_dir / "result.json"),
            "benchmark": heldout_benchmark,
        }
        if test_dataset == "GSE80672":
            shuffled_csv = heldout_dir / "predictions_shuffled_intervention.csv"
            shuffle_intervention(heldout_dir / "predictions.csv", shuffled_csv)
            heldout_results[test_dataset]["shuffled_benchmark"] = run_benchmark(
                shuffled_csv,
                heldout_dir / "benchmark_shuffled_intervention.json",
            )

    decision = {
        "best_config": best_name,
        "summary_path": str(summary_path),
        "random_label_benchmark": random_benchmark,
        "heldout_results": heldout_results,
        "v7_5_baselines": {
            "groupkfold_mae_weeks": 23.767,
            "gse121141_heldout_mae_weeks": 38.033,
            "gse121141_old_104w_mae_weeks": 75.386,
        },
    }
    (out_root / "v8_gse213628_ralph_decision.json").write_text(json.dumps(decision, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()

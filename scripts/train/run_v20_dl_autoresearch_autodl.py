#!/usr/bin/env python3
"""v20 AutoDL deep-learning autoresearch for the mouse methylation clock."""
from __future__ import annotations

import argparse
import itertools
import json
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


ROOT = Path("/home/zdq-as/mouse_methyl_work")
REMOTE_ROOT = Path("/root/mouse_methyl_work")
REMOTE_DATA_ROOT = Path("/root/autodl-tmp/mouse_methyl_work")
if REMOTE_ROOT.exists():
    ROOT = REMOTE_ROOT

DEFAULT_MATRIX = ROOT / "results" / "multidataset_v8_3_ablation" / "all6" / "all_rrbs_region_matrix_5kb.parquet"
if not DEFAULT_MATRIX.exists() and REMOTE_DATA_ROOT.exists():
    DEFAULT_MATRIX = (
        REMOTE_DATA_ROOT
        / "results"
        / "multidataset_v8_3_ablation"
        / "all6"
        / "all_rrbs_region_matrix_5kb.parquet"
    )
DEFAULT_METADATA = ROOT / "metadata" / "model_sample_metadata_v8.csv"
DEFAULT_OUT = (REMOTE_DATA_ROOT if REMOTE_DATA_ROOT.exists() else ROOT) / "results" / "autoresearch_v20_deep_learning"
REPORT = ROOT / "doc" / "20_analysis" / "51_20260520_v20_autodl_deep_learning_autoresearch_report.md"
TRAIN_SCRIPT = ROOT / "scripts" / "train" / "train_deep_clock.py"
BENCHMARK_SCRIPT = ROOT / "scripts" / "validate" / "benchmark_metrics.py"
DATASETS = ["GSE120137", "GSE80672", "GSE93957", "GSE121141", "GSE60012", "GSE213628"]
PYTHON = sys.executable

V19_BASELINE = {
    "config_id": "cfg_014_lgbm_top1000_quantile_uniform_median_s0p125_p0p8_ne200_lr0p03_leaves31",
    "group_pearson_r": 0.6315,
    "group_mae_weeks": 24.227,
    "group_r2": 0.2605,
    "lodo_mean_mae_weeks": 21.93,
    "lodo_worst_mae_weeks": 36.086,
    "lodo_cr_auc": 0.7793,
    "old_104w_weighted_mae_weeks": 43.762,
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str), encoding="utf-8")


def result_matches_config(result_path: Path, config: dict[str, Any], eval_mode: str) -> bool:
    if not result_path.exists():
        return False
    try:
        result = read_json(result_path)
    except Exception:
        return False
    training = result.get("training") or {}
    checks = {
        "architecture": config["architecture"],
        "n_feature_prefilter": str(config["n_feature_prefilter"]),
        "preprocess": config["preprocess"],
        "target_transform": config["target_transform"],
        "eval_mode": eval_mode,
    }
    for key, expected in checks.items():
        if str(result.get(key)) != str(expected):
            return False
    if int(training.get("epochs") or -1) != int(config["epochs"]):
        return False
    if int(training.get("patience") or -1) != int(config["patience"]):
        return False
    if float(training.get("lr") or -1.0) != float(config["lr"]):
        return False
    return True


def run(cmd: list[str], cwd: Path = ROOT) -> None:
    print(" ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=cwd, check=True)


def slug_value(value: Any) -> str:
    return str(value).replace(".", "p").replace(" ", "").replace("/", "-")


def config_id(config: dict[str, Any], idx: int) -> str:
    return "_".join(
        [
            f"cfg_{idx:03d}",
            config["architecture"],
            f"top{config['n_feature_prefilter']}",
            config["preprocess"],
            f"w{config['width']}",
            f"d{config['depth']}",
            f"do{slug_value(config['dropout'])}",
            f"lr{slug_value(config['lr'])}",
        ]
    )


def candidate_configs() -> list[dict[str, Any]]:
    base = {
        "imputation": "median",
        "min_train_feature_presence": 0.8,
        "min_train_group_feature_presence": 0.8,
        "max_train_dataset_mean_shift": 0.125,
        "min_train_agebin_feature_presence": 0.8,
        "target_transform": "log1p_days",
        "epochs": 300,
        "patience": 30,
        "batch_size": 64,
        "weight_decay": 1e-3,
        "huber_beta": 0.5,
        "amp": True,
    }
    configs: list[dict[str, Any]] = []
    for top, preprocess, width, depth, dropout, lr in itertools.product(
        ["500", "1000", "2000", "5000"],
        ["standard", "robust", "quantile_uniform"],
        [128, 256, 512],
        [2, 3, 4],
        [0.1, 0.25, 0.4],
        [1e-3, 3e-4],
    ):
        if len(configs) >= 45:
            break
        configs.append(
            {
                **base,
                "architecture": "res_mlp",
                "n_feature_prefilter": top,
                "preprocess": preprocess,
                "width": width,
                "depth": depth,
                "dropout": dropout,
                "lr": lr,
                "bottleneck": 64,
                "cnn_channels": 64,
                "cnn_kernel": 15,
                "token_dim": 16,
                "transformer_heads": 4,
                "transformer_layers": 2,
                "dae_latent": 64,
                "dae_noise": 0.1,
                "dae_pretrain_epochs": 0,
            }
        )
    for top in ["500", "1000", "2000"]:
        for width, bottleneck, dropout, lr in itertools.product([256, 512], [32, 64, 128], [0.1, 0.25], [1e-3, 3e-4]):
            configs.append(
                {
                    **base,
                    "architecture": "wide_deep",
                    "n_feature_prefilter": top,
                    "preprocess": "standard",
                    "width": width,
                    "depth": 3,
                    "dropout": dropout,
                    "lr": lr,
                    "bottleneck": bottleneck,
                    "cnn_channels": 64,
                    "cnn_kernel": 15,
                    "token_dim": 16,
                    "transformer_heads": 4,
                    "transformer_layers": 2,
                    "dae_latent": 64,
                    "dae_noise": 0.1,
                    "dae_pretrain_epochs": 0,
                }
            )
    for top, channels, kernel, dropout in itertools.product(["1000", "2000", "5000"], [32, 64, 128], [5, 15, 31], [0.1, 0.25]):
        configs.append(
            {
                **base,
                "architecture": "cnn1d",
                "n_feature_prefilter": top,
                "preprocess": "standard",
                "width": 256,
                "depth": 3,
                "dropout": dropout,
                "lr": 1e-3,
                "bottleneck": 64,
                "cnn_channels": channels,
                "cnn_kernel": kernel,
                "token_dim": 16,
                "transformer_heads": 4,
                "transformer_layers": 2,
                "dae_latent": 64,
                "dae_noise": 0.1,
                "dae_pretrain_epochs": 0,
            }
        )
    for top, token_dim, layers, dropout in itertools.product(["500", "1000"], [16, 32], [1, 2], [0.1, 0.25]):
        configs.append(
            {
                **base,
                "architecture": "tab_transformer",
                "n_feature_prefilter": top,
                "preprocess": "standard",
                "width": 256,
                "depth": 3,
                "dropout": dropout,
                "lr": 3e-4,
                "bottleneck": 64,
                "cnn_channels": 64,
                "cnn_kernel": 15,
                "token_dim": token_dim,
                "transformer_heads": 4,
                "transformer_layers": layers,
                "dae_latent": 64,
                "dae_noise": 0.1,
                "dae_pretrain_epochs": 0,
            }
        )
    for top, latent, noise, dropout in itertools.product(["1000", "2000"], [32, 64, 128], [0.05, 0.1, 0.2], [0.1, 0.25]):
        configs.append(
            {
                **base,
                "architecture": "dae_mlp",
                "n_feature_prefilter": top,
                "preprocess": "standard",
                "width": 256,
                "depth": 3,
                "dropout": dropout,
                "lr": 1e-3,
                "bottleneck": 64,
                "cnn_channels": 64,
                "cnn_kernel": 15,
                "token_dim": 16,
                "transformer_heads": 4,
                "transformer_layers": 2,
                "dae_latent": latent,
                "dae_noise": noise,
                "dae_pretrain_epochs": 20,
            }
        )
    for arch, top, width, depth, dropout, lr in itertools.product(
        ["res_mlp", "wide_deep"],
        ["500", "1000", "2000"],
        [256, 512],
        [3, 4],
        [0.1, 0.25],
        [1e-3, 3e-4],
    ):
        configs.append(
            {
                **base,
                "architecture": arch,
                "n_feature_prefilter": top,
                "preprocess": "robust",
                "width": width,
                "depth": depth,
                "dropout": dropout,
                "lr": lr,
                "bottleneck": 64,
                "cnn_channels": 64,
                "cnn_kernel": 15,
                "token_dim": 16,
                "transformer_heads": 4,
                "transformer_layers": 2,
                "dae_latent": 64,
                "dae_noise": 0.1,
                "dae_pretrain_epochs": 0,
            }
        )
    for top, channels, kernel in itertools.product(["1000", "2000"], [64, 128], [5, 15]):
        configs.append(
            {
                **base,
                "architecture": "cnn1d",
                "n_feature_prefilter": top,
                "preprocess": "robust",
                "width": 256,
                "depth": 3,
                "dropout": 0.1,
                "lr": 1e-3,
                "bottleneck": 64,
                "cnn_channels": channels,
                "cnn_kernel": kernel,
                "token_dim": 16,
                "transformer_heads": 4,
                "transformer_layers": 2,
                "dae_latent": 64,
                "dae_noise": 0.1,
                "dae_pretrain_epochs": 0,
            }
        )
    deduped = []
    seen = set()
    for cfg in configs:
        key = tuple(sorted(cfg.items()))
        if key not in seen:
            seen.add(key)
            deduped.append(cfg)
    return deduped


def train_cmd(
    config: dict[str, Any],
    matrix: Path,
    metadata: Path,
    out_dir: Path,
    eval_mode: str,
    heldout_dataset: str | None,
    seed: int,
    randomize: bool = False,
    epochs_override: int | None = None,
) -> list[str]:
    cmd = [
        PYTHON,
        str(TRAIN_SCRIPT),
        "--matrix_path",
        str(matrix),
        "--metadata_path",
        str(metadata),
        "--output_dir",
        str(out_dir),
        "--eval_mode",
        eval_mode,
        "--architecture",
        config["architecture"],
        "--n_feature_prefilter",
        str(config["n_feature_prefilter"]),
        "--imputation",
        config["imputation"],
        "--preprocess",
        config["preprocess"],
        "--min_train_feature_presence",
        str(config["min_train_feature_presence"]),
        "--min_train_group_feature_presence",
        str(config["min_train_group_feature_presence"]),
        "--max_train_dataset_mean_shift",
        str(config["max_train_dataset_mean_shift"]),
        "--min_train_agebin_feature_presence",
        str(config["min_train_agebin_feature_presence"]),
        "--target_transform",
        config["target_transform"],
        "--width",
        str(config["width"]),
        "--depth",
        str(config["depth"]),
        "--dropout",
        str(config["dropout"]),
        "--bottleneck",
        str(config["bottleneck"]),
        "--cnn_channels",
        str(config["cnn_channels"]),
        "--cnn_kernel",
        str(config["cnn_kernel"]),
        "--token_dim",
        str(config["token_dim"]),
        "--transformer_heads",
        str(config["transformer_heads"]),
        "--transformer_layers",
        str(config["transformer_layers"]),
        "--dae_latent",
        str(config["dae_latent"]),
        "--dae_noise",
        str(config["dae_noise"]),
        "--dae_pretrain_epochs",
        str(config["dae_pretrain_epochs"]),
        "--epochs",
        str(epochs_override or config["epochs"]),
        "--patience",
        str(config["patience"]),
        "--batch_size",
        str(config["batch_size"]),
        "--lr",
        str(config["lr"]),
        "--weight_decay",
        str(config["weight_decay"]),
        "--huber_beta",
        str(config["huber_beta"]),
        "--random_seed",
        str(seed),
    ]
    if config.get("amp"):
        cmd.append("--amp")
    if heldout_dataset:
        cmd.extend(["--heldout_dataset", heldout_dataset])
    if randomize:
        cmd.append("--randomize_labels")
    return cmd


def benchmark(predictions: Path, out_path: Path) -> dict[str, Any]:
    if not out_path.exists():
        run([PYTHON, str(BENCHMARK_SCRIPT), "--predictions", str(predictions), "--out", str(out_path)])
    return read_json(out_path)


def run_group_config(
    idx: int,
    cfg: dict[str, Any],
    matrix: Path,
    metadata: Path,
    out_root: Path,
) -> dict[str, Any]:
    cid = config_id(cfg, idx)
    cfg_dir = out_root / "groupkfold" / cid
    cfg_dir.mkdir(parents=True, exist_ok=True)
    write_json(cfg_dir / "config.json", cfg)
    if not result_matches_config(cfg_dir / "result.json", cfg, "groupkfold"):
        try:
            run(train_cmd(cfg, matrix, metadata, cfg_dir, "groupkfold", None, 42))
        except subprocess.CalledProcessError as exc:
            write_json(
                cfg_dir / "failed.json",
                {
                    "config_id": cid,
                    "status": "failed",
                    "returncode": exc.returncode,
                    "cmd": exc.cmd,
                    "config": cfg,
                    "timestamp": utc_now(),
                },
            )
            return {
                "config_id": cid,
                "status": "failed",
                **cfg,
                "error": f"returncode={exc.returncode}",
                "group_score": -1.0,
            }
    metrics = benchmark(cfg_dir / "predictions.csv", cfg_dir / "benchmark_result.json")
    result = read_json(cfg_dir / "result.json")
    return {
        "config_id": cid,
        "status": "completed",
        **cfg,
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
        "exec_time_sec": result.get("exec_time_sec"),
        "group_score": group_score(metrics),
    }


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
    cr_auc = float(row.get("lodo_cr_auc") or 0.5)
    old_mae = row.get("old_104w_weighted_mae_weeks")
    old_penalty = 0.0 if old_mae is None or pd.isna(old_mae) else min(float(old_mae), 100.0) / 100.0
    mean_mae_component = max(0.0, 1.0 - min(mean_lodo_mae, 60.0) / 60.0)
    worst_mae_component = max(0.0, 1.0 - min(worst_lodo_mae, 90.0) / 90.0)
    return round(
        (0.30 * mean_lodo_r)
        + (0.28 * mean_mae_component)
        + (0.18 * worst_mae_component)
        + (0.14 * cr_auc)
        - (0.10 * old_penalty),
        6,
    )


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
    final_metrics: dict[str, Any] | None,
) -> None:
    best = summary.iloc[0].to_dict() if not summary.empty else {}
    promoted = False
    if best:
        promoted = (
            float(best.get("lodo_mean_mae_weeks") or 999) < V19_BASELINE["lodo_mean_mae_weeks"]
            and float(best.get("lodo_worst_mae_weeks") or 999) <= V19_BASELINE["lodo_worst_mae_weeks"] + 3
        )
    lines = [
        "# v20 AutoDL Deep Learning Autoresearch Report",
        "",
        f"Date: {utc_now()}",
        "",
        "## Scope",
        "",
        "Controlled deep-learning POC on the current all6 5kb region matrix.",
        "Promotion requires LODO/old-age evidence against the v19 LGBM baseline.",
        "",
        "## Selected DL Candidate",
        "",
        "`not_available`" if not best else f"- config: `{best.get('config_id')}`",
        "" if not best else f"- architecture: `{best.get('architecture')}`",
        "" if not best else f"- GroupKFold r/MAE/R2: `{best.get('pearson_r')}` / `{best.get('mae_weeks')}` / `{best.get('r2')}`",
        "" if not best else f"- LODO mean/worst MAE: `{best.get('lodo_mean_mae_weeks')}` / `{best.get('lodo_worst_mae_weeks')}`",
        "" if not best else f"- LODO CR AUC: `{best.get('lodo_cr_auc')}`",
        "" if not best else f"- old104 weighted MAE: `{best.get('old_104w_weighted_mae_weeks')}`",
        "" if not best else f"- promoted over v19: `{promoted}`",
        "",
        "## v19 Baseline Comparator",
        "",
        f"`{V19_BASELINE}`",
        "",
        "## Random-Label Sanity",
        "",
        "`not_run`" if random_metrics is None else f"`{random_metrics}`",
        "",
        "## Final All-Data Fit",
        "",
        "`not_run`" if final_metrics is None else f"`{final_metrics}`",
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
        f"- `{out_root}/v20_dl_summary.csv`",
        f"- `{out_root}/v20_lodo_summary.csv`",
        f"- `{out_root}/final_all_data_model/final_deep_clock_model.pt`",
    ]
    text = "\n".join(line for line in lines if line is not None) + "\n"
    out_root.mkdir(parents=True, exist_ok=True)
    (out_root / "v20_autodl_deep_learning_autoresearch_report.md").write_text(text, encoding="utf-8")
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    try:
        REPORT.write_text(text, encoding="utf-8")
    except OSError:
        pass


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", default=str(DEFAULT_MATRIX))
    parser.add_argument("--metadata", default=str(DEFAULT_METADATA))
    parser.add_argument("--out-root", default=str(DEFAULT_OUT))
    parser.add_argument("--budget-hours", type=float, default=24.0)
    parser.add_argument("--max-configs", type=int, default=0)
    parser.add_argument("--top-lodo", type=int, default=12)
    parser.add_argument("--top-seed-repeats", type=int, default=3)
    parser.add_argument("--workers", type=int, default=1, help="Concurrent GroupKFold configs on one GPU.")
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--skip-random-label", action="store_true")
    parser.add_argument("--skip-final-fit", action="store_true")
    args = parser.parse_args()

    matrix = Path(args.matrix)
    metadata = Path(args.metadata)
    out_root = Path(args.out_root)
    out_root.mkdir(parents=True, exist_ok=True)
    started = time.time()
    deadline = started + args.budget_hours * 3600
    configs = candidate_configs()
    if args.smoke:
        configs = configs[:2]
        for cfg in configs:
            cfg["epochs"] = 2
            cfg["patience"] = 1
            cfg["dae_pretrain_epochs"] = min(int(cfg["dae_pretrain_epochs"]), 1)
    if args.max_configs > 0:
        configs = configs[: args.max_configs]
    write_json(
        out_root / "v20_run_manifest.json",
        {
            "timestamp": utc_now(),
            "python": PYTHON,
            "matrix": str(matrix),
            "metadata": str(metadata),
            "budget_hours": args.budget_hours,
            "n_configs": len(configs),
            "configs": configs,
        },
    )

    rows = []
    id_to_config = {}
    for idx, cfg in enumerate(configs, start=1):
        cid = config_id(cfg, idx)
        id_to_config[cid] = cfg

    group_jobs = list(enumerate(configs, start=1))
    max_workers = max(1, int(args.workers))
    if max_workers == 1:
        for idx, cfg in group_jobs:
            if time.time() > deadline and not args.smoke:
                print("[Budget] stopping GroupKFold search.", flush=True)
                break
            rows.append(run_group_config(idx, cfg, matrix, metadata, out_root))
            pd.DataFrame(rows).sort_values(["group_score", "pearson_r"], ascending=False).to_csv(
                out_root / "v20_dl_summary_partial.csv", index=False
            )
    else:
        with ThreadPoolExecutor(max_workers=max_workers) as pool:
            pending = {}
            cursor = 0
            while cursor < len(group_jobs) or pending:
                while cursor < len(group_jobs) and len(pending) < max_workers:
                    if time.time() > deadline and not args.smoke:
                        break
                    idx, cfg = group_jobs[cursor]
                    future = pool.submit(run_group_config, idx, cfg, matrix, metadata, out_root)
                    pending[future] = (idx, cfg)
                    cursor += 1
                if not pending:
                    break
                for future in as_completed(list(pending.keys()), timeout=None):
                    pending.pop(future)
                    rows.append(future.result())
                    pd.DataFrame(rows).sort_values(["group_score", "pearson_r"], ascending=False).to_csv(
                        out_root / "v20_dl_summary_partial.csv", index=False
                    )
                    break
                if time.time() > deadline and not args.smoke and cursor >= len(group_jobs):
                    break

    summary = pd.DataFrame(rows).sort_values(["group_score", "pearson_r"], ascending=False).reset_index(drop=True)
    lodo_rows = []
    for cid in summary.head(2 if args.smoke else args.top_lodo)["config_id"].tolist():
        if time.time() > deadline and not args.smoke:
            print("[Budget] stopping LODO search.", flush=True)
            break
        cfg = id_to_config[cid]
        config_lodo_rows = []
        for dataset in DATASETS:
            heldout_dir = out_root / "lodo" / cid / dataset
            if not result_matches_config(heldout_dir / "result.json", cfg, "heldout"):
                run(train_cmd(cfg, matrix, metadata, heldout_dir, "heldout", dataset, 42))
            metrics = benchmark(heldout_dir / "predictions.csv", heldout_dir / "benchmark_result.json")
            pred = pd.read_csv(heldout_dir / "predictions.csv")
            old = pred[pred["age_weeks_true"] >= 104]
            row = {
                "config_id": cid,
                "dataset": dataset,
                **{key: metrics.get(key) for key in ["n_samples", "pearson_r", "mae_weeks", "rmse_weeks", "r2", "cr_detection_auc"]},
                "old_104w_n": int(len(old)),
                "old_104w_mae_weeks": None
                if old.empty
                else round(float(np.mean(np.abs(old["age_weeks_true"] - old["age_weeks_pred"]))), 3),
            }
            lodo_rows.append(row)
            config_lodo_rows.append(row)
        lodo_summary = summarize_lodo(config_lodo_rows)
        mask = summary["config_id"].eq(cid)
        for key, value in lodo_summary.items():
            summary.loc[mask, key] = value
        merged = summary.loc[mask].iloc[0].to_dict()
        summary.loc[mask, "lodo_score"] = lodo_score({**merged, **lodo_summary})

    summary["lodo_score"] = summary["lodo_score"].fillna(-1.0)
    summary["final_score"] = summary["group_score"] + summary["lodo_score"].clip(lower=0.0)
    summary = summary.sort_values(["final_score", "group_score"], ascending=False).reset_index(drop=True)

    if not args.smoke and args.top_seed_repeats > 1 and not summary.empty:
        for cid in summary.head(3)["config_id"].tolist():
            cfg = id_to_config[cid]
            seed_scores = []
            for seed in [43, 44][: max(0, args.top_seed_repeats - 1)]:
                seed_dir = out_root / "seed_repeats" / cid / f"seed_{seed}"
                if not result_matches_config(seed_dir / "result.json", cfg, "groupkfold"):
                    run(train_cmd(cfg, matrix, metadata, seed_dir, "groupkfold", None, seed))
                metrics = benchmark(seed_dir / "predictions.csv", seed_dir / "benchmark_result.json")
                seed_scores.append(group_score(metrics))
            mask = summary["config_id"].eq(cid)
            if seed_scores:
                summary.loc[mask, "seed_repeat_group_score_mean"] = round(float(np.mean(seed_scores)), 6)

    random_metrics = None
    final_metrics = None
    if not summary.empty:
        best_id = str(summary.iloc[0]["config_id"])
        best_cfg = id_to_config[best_id]
        if not args.skip_random_label:
            random_dir = out_root / "random_label_sanity" / best_id
            if not result_matches_config(random_dir / "result.json", best_cfg, "groupkfold"):
                run(train_cmd(best_cfg, matrix, metadata, random_dir, "groupkfold", None, 99, randomize=True))
            random_metrics = benchmark(random_dir / "predictions.csv", random_dir / "benchmark_result.json")
        if not args.skip_final_fit:
            final_dir = out_root / "final_all_data_model"
            if not result_matches_config(final_dir / "result.json", best_cfg, "final"):
                run(train_cmd(best_cfg, matrix, metadata, final_dir, "final", None, 42))
            final_metrics = read_json(final_dir / "result.json")

    lodo_df = pd.DataFrame(lodo_rows)
    summary.to_csv(out_root / "v20_dl_summary.csv", index=False)
    summary.to_json(out_root / "v20_dl_summary.json", orient="records", indent=2)
    lodo_df.to_csv(out_root / "v20_lodo_summary.csv", index=False)
    write_report(out_root, summary, lodo_df, random_metrics, final_metrics)
    print(json.dumps(summary.iloc[0].to_dict() if not summary.empty else {}, indent=2, default=str), flush=True)
    print(f"[Done] elapsed_sec={round(time.time() - started, 1)} out={out_root}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

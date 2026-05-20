#!/usr/bin/env python3
"""Robust v23 cloud batch for DL age optimization plus biological-signal sidecars."""
from __future__ import annotations

import argparse
import itertools
import json
import os
import socket
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


LOCAL_ROOT = Path("/home/zdq-as/mouse_methyl_work")
REMOTE_ROOT = Path("/root/autodl-tmp/mouse_methyl_work")
ROOT = REMOTE_ROOT if REMOTE_ROOT.exists() else LOCAL_ROOT
DEFAULT_MATRIX = ROOT / "results" / "multidataset_v8_3_ablation" / "all6" / "all_rrbs_region_matrix_5kb.parquet"
DEFAULT_METADATA = ROOT / "metadata" / "model_sample_metadata_v8.csv"
DEFAULT_TARGET_REGISTRY = ROOT / "metadata" / "v22_biological_signal_targets.csv"
DEFAULT_OUT = ROOT / "results" / "v23_cloud_max_batch"
REPORT = ROOT / "doc" / "20_analysis" / "58_20260521_v23_cloud_max_batch_report.md"
TRAIN_SCRIPT = ROOT / "scripts" / "train" / "train_deep_clock.py"
BENCHMARK_SCRIPT = ROOT / "scripts" / "validate" / "benchmark_metrics.py"
BIO_SCREEN_SCRIPT = ROOT / "scripts" / "validate" / "run_v22_biological_signal_screen.py"
FEATURE_SCRIPT = ROOT / "scripts" / "validate" / "run_v21_feature_interpretation.py"
DATASETS = ["GSE120137", "GSE80672", "GSE93957", "GSE121141", "GSE60012", "GSE213628"]
PYTHON = sys.executable


V20_BEST = {
    "config_id": "cfg_241_res_mlp_top1000_robust_w256_d3_do0p1_lr0p0003",
    "group_mae_weeks": 20.318,
    "lodo_mean_mae_weeks": 18.109,
    "lodo_worst_mae_weeks": 30.915,
    "lodo_cr_auc": 0.8459,
    "old_104w_weighted_mae_weeks": 27.07,
}

V22_BEST = {
    "config_id": "v22_007_res_mlp_top1000_robust_w192_d3_do0p1_lr0p0003",
    "group_mae_weeks": 20.45,
    "lodo_mean_mae_weeks": 18.81,
    "lodo_worst_mae_weeks": 27.043,
    "lodo_cr_auc": 0.8349,
    "old_104w_weighted_mae_weeks": 25.673,
}


def utc_now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str), encoding="utf-8")
    tmp.replace(path)


def append_jsonl(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(payload, sort_keys=True, default=str) + "\n")


def slug(value: Any) -> str:
    return str(value).replace(".", "p").replace("/", "-").replace(" ", "")


def config_id(cfg: dict[str, Any], idx: int) -> str:
    arch_bits = [f"v23_{idx:03d}", cfg["architecture"], f"top{cfg['n_feature_prefilter']}", cfg["preprocess"]]
    if cfg["architecture"] == "cnn1d":
        arch_bits.extend([f"ch{cfg['cnn_channels']}", f"k{cfg['cnn_kernel']}"])
    elif cfg["architecture"] == "tab_transformer":
        arch_bits.extend([f"tok{cfg['token_dim']}", f"tl{cfg['transformer_layers']}"])
    elif cfg["architecture"] == "dae_mlp":
        arch_bits.extend([f"w{cfg['width']}", f"lat{cfg['dae_latent']}"])
    elif cfg["architecture"] == "wide_deep":
        arch_bits.extend([f"w{cfg['width']}", f"bn{cfg['bottleneck']}"])
    else:
        arch_bits.extend([f"w{cfg['width']}", f"d{cfg['depth']}"])
    arch_bits.extend([f"do{slug(cfg['dropout'])}", f"lr{slug(cfg['lr'])}"])
    return "_".join(arch_bits)


def base_config() -> dict[str, Any]:
    return {
        "imputation": "median",
        "min_train_feature_presence": 0.8,
        "min_train_group_feature_presence": 0.8,
        "max_train_dataset_mean_shift": 0.125,
        "min_train_agebin_feature_presence": 0.8,
        "target_transform": "log1p_days",
        "batch_size": 64,
        "weight_decay": 1e-3,
        "huber_beta": 0.5,
        "amp": True,
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


def candidate_configs() -> list[dict[str, Any]]:
    base = base_config()
    configs: list[dict[str, Any]] = []

    for top, width, depth in itertools.product(["1000", "1500", "2000", "5000"], [128, 192, 256], [2, 3]):
        configs.append(
            {
                **base,
                "architecture": "res_mlp",
                "n_feature_prefilter": top,
                "preprocess": "robust",
                "width": width,
                "depth": depth,
                "dropout": 0.1,
                "lr": 3e-4,
            }
        )
    for top, preprocess, width in itertools.product(["1000", "2000", "5000"], ["standard", "quantile_uniform"], [192, 256]):
        configs.append(
            {
                **base,
                "architecture": "res_mlp",
                "n_feature_prefilter": top,
                "preprocess": preprocess,
                "width": width,
                "depth": 3,
                "dropout": 0.1,
                "lr": 3e-4,
            }
        )
    for top, dropout, lr in itertools.product(["1000", "2000"], [0.05, 0.15, 0.25], [3e-4, 1e-4]):
        configs.append(
            {
                **base,
                "architecture": "res_mlp",
                "n_feature_prefilter": top,
                "preprocess": "robust",
                "width": 192,
                "depth": 3,
                "dropout": dropout,
                "lr": lr,
            }
        )
    for top, preprocess, width, bottleneck in itertools.product(["1000", "2000", "5000"], ["standard", "robust"], [256, 512], [64]):
        configs.append(
            {
                **base,
                "architecture": "wide_deep",
                "n_feature_prefilter": top,
                "preprocess": preprocess,
                "width": width,
                "depth": 3,
                "dropout": 0.1,
                "lr": 3e-4,
                "bottleneck": bottleneck,
            }
        )
    for top, channels, kernel in itertools.product(["1000", "2000"], [32, 64], [5, 15]):
        configs.append(
            {
                **base,
                "architecture": "cnn1d",
                "n_feature_prefilter": top,
                "preprocess": "standard",
                "width": 256,
                "depth": 3,
                "dropout": 0.1,
                "lr": 1e-3,
                "cnn_channels": channels,
                "cnn_kernel": kernel,
            }
        )
    for top, latent in itertools.product(["1000", "2000", "5000"], [32, 64]):
        configs.append(
            {
                **base,
                "architecture": "dae_mlp",
                "n_feature_prefilter": top,
                "preprocess": "robust",
                "width": 256,
                "depth": 3,
                "dropout": 0.1,
                "lr": 3e-4,
                "dae_latent": latent,
                "dae_noise": 0.1,
                "dae_pretrain_epochs": 8,
            }
        )
    for top, token_dim in itertools.product(["500", "1000"], [16, 32]):
        configs.append(
            {
                **base,
                "architecture": "tab_transformer",
                "n_feature_prefilter": top,
                "preprocess": "standard",
                "width": 256,
                "depth": 3,
                "dropout": 0.1,
                "lr": 3e-4,
                "token_dim": token_dim,
                "transformer_heads": 4,
                "transformer_layers": 1,
            }
        )
    return configs


def train_cmd(
    cfg: dict[str, Any],
    matrix: Path,
    metadata: Path,
    out_dir: Path,
    eval_mode: str,
    *,
    heldout_dataset: str | None = None,
    randomize: bool = False,
    seed: int = 42,
    epochs: int = 300,
    patience: int = 30,
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
        str(cfg["architecture"]),
        "--n_feature_prefilter",
        str(cfg["n_feature_prefilter"]),
        "--imputation",
        str(cfg["imputation"]),
        "--preprocess",
        str(cfg["preprocess"]),
        "--min_train_feature_presence",
        str(cfg["min_train_feature_presence"]),
        "--min_train_group_feature_presence",
        str(cfg["min_train_group_feature_presence"]),
        "--max_train_dataset_mean_shift",
        str(cfg["max_train_dataset_mean_shift"]),
        "--min_train_agebin_feature_presence",
        str(cfg["min_train_agebin_feature_presence"]),
        "--target_transform",
        str(cfg["target_transform"]),
        "--width",
        str(cfg["width"]),
        "--depth",
        str(cfg["depth"]),
        "--dropout",
        str(cfg["dropout"]),
        "--bottleneck",
        str(cfg["bottleneck"]),
        "--cnn_channels",
        str(cfg["cnn_channels"]),
        "--cnn_kernel",
        str(cfg["cnn_kernel"]),
        "--token_dim",
        str(cfg["token_dim"]),
        "--transformer_heads",
        str(cfg["transformer_heads"]),
        "--transformer_layers",
        str(cfg["transformer_layers"]),
        "--dae_latent",
        str(cfg["dae_latent"]),
        "--dae_noise",
        str(cfg["dae_noise"]),
        "--dae_pretrain_epochs",
        str(cfg["dae_pretrain_epochs"]),
        "--epochs",
        str(epochs),
        "--patience",
        str(patience),
        "--batch_size",
        str(cfg["batch_size"]),
        "--lr",
        str(cfg["lr"]),
        "--weight_decay",
        str(cfg["weight_decay"]),
        "--huber_beta",
        str(cfg["huber_beta"]),
        "--random_seed",
        str(seed),
        "--device",
        "cuda",
    ]
    if cfg.get("amp"):
        cmd.append("--amp")
    if heldout_dataset:
        cmd.extend(["--heldout_dataset", heldout_dataset])
    if randomize:
        cmd.append("--randomize_labels")
    return cmd


def result_matches(result_path: Path, cfg: dict[str, Any], eval_mode: str, epochs: int, *, randomize: bool = False) -> bool:
    if not result_path.exists():
        return False
    try:
        result = read_json(result_path)
    except Exception:
        return False
    training = result.get("training") or {}
    checks = {
        "architecture": cfg["architecture"],
        "n_feature_prefilter": str(cfg["n_feature_prefilter"]),
        "preprocess": cfg["preprocess"],
        "target_transform": cfg["target_transform"],
        "eval_mode": eval_mode,
        "randomize_labels": randomize,
    }
    for key, expected in checks.items():
        if str(result.get(key)) != str(expected):
            return False
    return int(training.get("epochs") or -1) == int(epochs)


def run_logged(cmd: list[str], log_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.setdefault("PYTHONUNBUFFERED", "1")
    env.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")
    with log_path.open("a", encoding="utf-8") as log:
        log.write(f"\n[{utc_now()}] $ {' '.join(cmd)}\n")
        log.flush()
        subprocess.run(cmd, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)


def benchmark(predictions: Path, out_path: Path, log_path: Path) -> dict[str, Any]:
    if not out_path.exists() or out_path.stat().st_size == 0:
        run_logged([PYTHON, str(BENCHMARK_SCRIPT), "--predictions", str(predictions), "--out", str(out_path)], log_path)
    return read_json(out_path)


def group_score(metrics: dict[str, Any]) -> float:
    r = float(metrics.get("pearson_r") or 0.0)
    r2 = max(0.0, float(metrics.get("r2") or 0.0))
    mae = float(metrics.get("mae_weeks") or 999.0)
    cr_auc = float(metrics.get("cr_detection_auc") or 0.5)
    mae_component = max(0.0, 1.0 - min(mae, 60.0) / 60.0)
    return round(0.36 * r + 0.24 * r2 + 0.22 * mae_component + 0.18 * cr_auc, 6)


def lodo_score(row: dict[str, Any]) -> float:
    mean_mae = float(row.get("lodo_mean_mae_weeks") or 999.0)
    worst_mae = float(row.get("lodo_worst_mae_weeks") or 999.0)
    mean_r = float(row.get("lodo_mean_pearson_r") or 0.0)
    cr_auc = float(row.get("lodo_cr_auc") or 0.5)
    old_mae = row.get("old_104w_weighted_mae_weeks")
    old_penalty = 0.0 if old_mae is None or pd.isna(old_mae) else min(float(old_mae), 100.0) / 100.0
    return round(
        0.30 * mean_r
        + 0.30 * max(0.0, 1.0 - min(mean_mae, 60.0) / 60.0)
        + 0.18 * max(0.0, 1.0 - min(worst_mae, 90.0) / 90.0)
        + 0.14 * cr_auc
        - 0.08 * old_penalty,
        6,
    )


def summarize_lodo(rows: list[dict[str, Any]]) -> dict[str, Any]:
    df = pd.DataFrame(rows)
    if df.empty:
        return {}
    old = df[df["old_104w_n"] > 0].copy()
    cr = df[df["cr_detection_auc"].notna()]
    return {
        "lodo_mean_mae_weeks": round(float(df["mae_weeks"].mean()), 3),
        "lodo_worst_mae_weeks": round(float(df["mae_weeks"].max()), 3),
        "lodo_mean_pearson_r": round(float(df["pearson_r"].mean()), 4),
        "lodo_min_pearson_r": round(float(df["pearson_r"].min()), 4),
        "lodo_cr_auc": None if cr.empty else round(float(cr["cr_detection_auc"].mean()), 4),
        "old_104w_weighted_mae_weeks": None
        if old.empty
        else round(float(np.average(old["old_104w_mae_weeks"], weights=old["old_104w_n"])), 3),
    }


def compact_result_fields(metrics: dict[str, Any]) -> dict[str, Any]:
    keys = [
        "n_samples",
        "pearson_r",
        "mae_weeks",
        "medae_weeks",
        "rmse_weeks",
        "r2",
        "cross_dataset_mae",
        "cr_detection_auc",
        "cr_cohens_d",
    ]
    return {key: metrics.get(key) for key in keys if key in metrics}


def run_train_task(task: dict[str, Any], matrix: Path, metadata: Path, epochs: int, patience: int, max_retries: int) -> dict[str, Any]:
    out_dir = Path(task["out_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    cfg = task["config"]
    write_json(out_dir / "config.json", {"task": task, "config": cfg})
    log_path = out_dir / "task.log"
    result_path = out_dir / "result.json"
    bench_path = out_dir / "benchmark_result.json"
    eval_mode = str(task["eval_mode"])
    randomize = bool(task.get("randomize", False))
    status_path = out_dir / "task_status.json"

    try:
        if result_matches(result_path, cfg, eval_mode, epochs, randomize=randomize) and bench_path.exists():
            metrics = read_json(bench_path)
            status = "skipped_existing"
        else:
            status = "completed"
            for attempt in range(1, max_retries + 2):
                write_json(
                    status_path,
                    {
                        "status": "running",
                        "attempt": attempt,
                        "timestamp": utc_now(),
                        "pid": os.getpid(),
                        "task_id": task["task_id"],
                    },
                )
                try:
                    run_logged(
                        train_cmd(
                            cfg,
                            matrix,
                            metadata,
                            out_dir,
                            eval_mode,
                            heldout_dataset=task.get("heldout_dataset"),
                            randomize=randomize,
                            seed=int(task.get("seed", 42)),
                            epochs=epochs,
                            patience=patience,
                        ),
                        log_path,
                    )
                    metrics = benchmark(out_dir / "predictions.csv", bench_path, log_path)
                    break
                except subprocess.CalledProcessError as exc:
                    write_json(
                        out_dir / "failed.json",
                        {
                            "timestamp": utc_now(),
                            "attempt": attempt,
                            "returncode": exc.returncode,
                            "task": task,
                            "cmd": exc.cmd,
                        },
                    )
                    if attempt > max_retries:
                        raise
                    time.sleep(20)
        result = read_json(result_path) if result_path.exists() else {}
        row = {
            "task_id": task["task_id"],
            "config_id": task["config_id"],
            "task_kind": task["task_kind"],
            "status": status,
            "eval_mode": eval_mode,
            "heldout_dataset": task.get("heldout_dataset"),
            "randomize": randomize,
            "out_dir": str(out_dir),
            **cfg,
            **compact_result_fields(metrics),
            "exec_time_sec": result.get("exec_time_sec"),
        }
        write_json(status_path, {"status": status, "timestamp": utc_now(), "row": row})
        return row
    except Exception as exc:
        row = {
            "task_id": task["task_id"],
            "config_id": task["config_id"],
            "task_kind": task["task_kind"],
            "status": "failed",
            "error": repr(exc),
            "eval_mode": eval_mode,
            "heldout_dataset": task.get("heldout_dataset"),
            "randomize": randomize,
            "out_dir": str(out_dir),
            **cfg,
        }
        write_json(status_path, {"status": "failed", "timestamp": utc_now(), "row": row})
        return row


def run_tasks_parallel(
    tasks: list[dict[str, Any]],
    *,
    matrix: Path,
    metadata: Path,
    epochs: int,
    patience: int,
    workers: int,
    max_retries: int,
    event_log: Path,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if not tasks:
        return rows
    with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
        future_to_task = {
            pool.submit(run_train_task, task, matrix, metadata, epochs, patience, max_retries): task for task in tasks
        }
        for future in as_completed(future_to_task):
            task = future_to_task[future]
            try:
                row = future.result()
            except Exception as exc:  # pragma: no cover - defensive guard
                row = {"task_id": task["task_id"], "config_id": task["config_id"], "status": "failed", "error": repr(exc)}
            rows.append(row)
            append_jsonl(event_log, {"timestamp": utc_now(), "event": "task_finished", "row": row})
    return rows


def make_group_tasks(configs: list[dict[str, Any]], out_root: Path) -> tuple[list[dict[str, Any]], dict[str, dict[str, Any]]]:
    tasks = []
    config_lookup: dict[str, dict[str, Any]] = {}
    for idx, cfg in enumerate(configs, start=1):
        cid = config_id(cfg, idx)
        config_lookup[cid] = cfg
        tasks.append(
            {
                "task_id": f"group::{cid}",
                "task_kind": "group",
                "config_id": cid,
                "config": cfg,
                "eval_mode": "groupkfold",
                "out_dir": str(out_root / "groupkfold" / cid),
                "seed": 42,
            }
        )
    return tasks, config_lookup


def make_lodo_tasks(summary: pd.DataFrame, config_lookup: dict[str, dict[str, Any]], out_root: Path, top_lodo: int) -> list[dict[str, Any]]:
    tasks = []
    for cid in summary.head(top_lodo)["config_id"].astype(str).tolist():
        cfg = config_lookup[cid]
        for dataset in DATASETS:
            tasks.append(
                {
                    "task_id": f"lodo::{cid}::{dataset}",
                    "task_kind": "lodo",
                    "config_id": cid,
                    "config": cfg,
                    "eval_mode": "heldout",
                    "heldout_dataset": dataset,
                    "out_dir": str(out_root / "lodo" / cid / dataset),
                    "seed": 42,
                }
            )
    return tasks


def make_random_tasks(summary: pd.DataFrame, config_lookup: dict[str, dict[str, Any]], out_root: Path, random_top: int) -> list[dict[str, Any]]:
    tasks = []
    for cid in summary.head(random_top)["config_id"].astype(str).tolist():
        cfg = config_lookup[cid]
        tasks.append(
            {
                "task_id": f"random::{cid}",
                "task_kind": "random_label",
                "config_id": cid,
                "config": cfg,
                "eval_mode": "groupkfold",
                "out_dir": str(out_root / "random_label_sanity" / cid),
                "seed": 99,
                "randomize": True,
            }
        )
    return tasks


def make_final_task(summary: pd.DataFrame, config_lookup: dict[str, dict[str, Any]], out_root: Path) -> list[dict[str, Any]]:
    if summary.empty:
        return []
    cid = str(summary.iloc[0]["config_id"])
    return [
        {
            "task_id": f"final::{cid}",
            "task_kind": "final",
            "config_id": cid,
            "config": config_lookup[cid],
            "eval_mode": "final",
            "out_dir": str(out_root / "final_all_data_model" / cid),
            "seed": 42,
        }
    ]


def rows_to_group_summary(rows: list[dict[str, Any]]) -> pd.DataFrame:
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df["group_score"] = df.apply(lambda row: group_score(row.to_dict()) if row.get("status") != "failed" else -1.0, axis=1)
    return df.sort_values(["group_score", "pearson_r"], ascending=False, na_position="last").reset_index(drop=True)


def attach_lodo_summary(summary: pd.DataFrame, lodo_rows: list[dict[str, Any]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    if summary.empty:
        return summary, pd.DataFrame(lodo_rows)
    lodo = pd.DataFrame(lodo_rows)
    if lodo.empty:
        summary["lodo_score"] = -1.0
        summary["final_score"] = summary["group_score"]
        return summary, lodo
    enriched = summary.copy()
    for cid, group in lodo[lodo["status"].ne("failed")].groupby("config_id"):
        records = group.to_dict(orient="records")
        lodo_stats = summarize_lodo(records)
        mask = enriched["config_id"].eq(cid)
        for key, value in lodo_stats.items():
            enriched.loc[mask, key] = value
        merged = enriched.loc[mask].iloc[0].to_dict()
        enriched.loc[mask, "lodo_score"] = lodo_score({**merged, **lodo_stats})
    enriched["lodo_score"] = enriched["lodo_score"].fillna(-1.0)
    enriched["final_score"] = enriched["group_score"] + enriched["lodo_score"].clip(lower=0.0)
    enriched = enriched.sort_values(["final_score", "group_score"], ascending=False, na_position="last").reset_index(drop=True)
    return enriched, lodo


def add_old104_metrics(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    enriched = []
    for row in rows:
        item = dict(row)
        if item.get("status") == "failed" or item.get("task_kind") != "lodo":
            enriched.append(item)
            continue
        pred_path = Path(item.get("out_dir", "")) / "predictions.csv" if item.get("out_dir") else Path()
        if not pred_path.exists():
            pred_path = Path(DEFAULT_OUT) / "lodo" / str(item["config_id"]) / str(item.get("heldout_dataset")) / "predictions.csv"
        if pred_path.exists():
            pred = pd.read_csv(pred_path)
            old = pred[pred["age_weeks_true"] >= 104]
            item["old_104w_n"] = int(len(old))
            item["old_104w_mae_weeks"] = None if old.empty else round(float(np.mean(np.abs(old["age_weeks_true"] - old["age_weeks_pred"]))), 3)
        else:
            item["old_104w_n"] = 0
            item["old_104w_mae_weeks"] = None
        enriched.append(item)
    return enriched


def run_bio_sidecars(args: argparse.Namespace, out_root: Path, event_log: Path) -> list[dict[str, Any]]:
    if args.skip_bio_sidecars:
        return []
    sidecar_root = out_root / "biology_sidecars"
    sidecar_root.mkdir(parents=True, exist_ok=True)
    jobs = [
        {
            "name": "auto_robust2000",
            "cmd": [
                PYTHON,
                str(BIO_SCREEN_SCRIPT),
                "--matrix",
                str(args.matrix),
                "--metadata",
                str(args.metadata),
                "--target-registry",
                str(args.target_registry),
                "--out-dir",
                str(sidecar_root / "auto_robust2000"),
                "--models",
                "auto",
                "--preprocess",
                "robust",
                "--n-feature-prefilter",
                "2000",
                "--include-random-labels",
            ],
        },
        {
            "name": "lgbm_age_2000",
            "cmd": [
                PYTHON,
                str(BIO_SCREEN_SCRIPT),
                "--matrix",
                str(args.matrix),
                "--metadata",
                str(args.metadata),
                "--target-registry",
                str(args.target_registry),
                "--out-dir",
                str(sidecar_root / "lgbm_age_2000"),
                "--targets",
                "age_weeks",
                "--models",
                "lgbm",
                "--preprocess",
                "none",
                "--n-feature-prefilter",
                "2000",
                "--include-random-labels",
            ],
        },
        {
            "name": "lgbm_classifiers_2000",
            "cmd": [
                PYTHON,
                str(BIO_SCREEN_SCRIPT),
                "--matrix",
                str(args.matrix),
                "--metadata",
                str(args.metadata),
                "--target-registry",
                str(args.target_registry),
                "--out-dir",
                str(sidecar_root / "lgbm_classifiers_2000"),
                "--targets",
                "sex_binary,tissue_multiclass,cr_vs_control,old_castration_vs_control,condition_family_multiclass",
                "--models",
                "lgbm",
                "--preprocess",
                "none",
                "--n-feature-prefilter",
                "2000",
                "--include-random-labels",
            ],
        },
        {
            "name": "rf_classifiers_1000",
            "cmd": [
                PYTHON,
                str(BIO_SCREEN_SCRIPT),
                "--matrix",
                str(args.matrix),
                "--metadata",
                str(args.metadata),
                "--target-registry",
                str(args.target_registry),
                "--out-dir",
                str(sidecar_root / "rf_classifiers_1000"),
                "--targets",
                "sex_binary,tissue_multiclass,cr_vs_control,old_castration_vs_control,condition_family_multiclass",
                "--models",
                "rf",
                "--preprocess",
                "none",
                "--n-feature-prefilter",
                "1000",
                "--include-random-labels",
            ],
        },
    ]
    if args.smoke:
        jobs = jobs[:1]
        jobs[0]["cmd"].append("--smoke")

    rows = []
    for job in jobs:
        log_path = sidecar_root / f"{job['name']}.log"
        summary_path = sidecar_root / job["name"] / "v22_biological_signal_screen_summary.json"
        if summary_path.exists():
            summary = read_json(summary_path)
            row = {"name": job["name"], "status": "skipped_existing", **summary}
            rows.append(row)
            append_jsonl(event_log, {"timestamp": utc_now(), "event": "bio_sidecar_finished", "row": row})
            continue
        try:
            run_logged(job["cmd"], log_path)
            summary = read_json(summary_path) if summary_path.exists() else {}
            row = {"name": job["name"], "status": "completed", **summary}
        except Exception as exc:
            row = {"name": job["name"], "status": "failed", "error": repr(exc)}
        rows.append(row)
        append_jsonl(event_log, {"timestamp": utc_now(), "event": "bio_sidecar_finished", "row": row})
    pd.DataFrame(rows).to_csv(out_root / "biology_sidecar_summary.csv", index=False)
    return rows


def run_feature_interpretation(args: argparse.Namespace, out_root: Path, final_rows: list[dict[str, Any]], event_log: Path) -> dict[str, Any] | None:
    if args.skip_feature_interpretation or not final_rows:
        return None
    best = final_rows[0]
    if best.get("status") == "failed":
        return None
    selected = Path(best["out_dir"]) / "selected_features.csv"
    if not selected.exists():
        return None
    ml_features = ROOT / "results" / "autoresearch_v19_ml_baseline" / "final_all_data_model" / "selected_features.csv"
    if not ml_features.exists():
        ml_features = ROOT / "results" / "v22_biological_signal_screen" / "fold_selected_features.csv"
    if not ml_features.exists():
        ml_features = selected
    out_dir = out_root / "feature_interpretation"
    cmd = [
        PYTHON,
        str(FEATURE_SCRIPT),
        "--matrix",
        str(args.matrix),
        "--metadata",
        str(args.metadata),
        "--ml-selected-features",
        str(ml_features),
        "--dl-selected-features",
        str(selected),
        "--out-dir",
        str(out_dir),
        "--top-n",
        str(args.feature_top_n),
    ]
    try:
        run_logged(cmd, out_dir / "feature_interpretation.log")
        summary = read_json(out_dir / "v21_feature_interpretation_summary.json")
        append_jsonl(event_log, {"timestamp": utc_now(), "event": "feature_interpretation_finished", "row": summary})
        return summary
    except Exception as exc:
        summary = {"status": "failed", "error": repr(exc)}
        append_jsonl(event_log, {"timestamp": utc_now(), "event": "feature_interpretation_failed", "row": summary})
        return summary


def md_table(df: pd.DataFrame, cols: list[str], n: int = 12) -> str:
    if df.empty:
        return "No rows."
    use_cols = [col for col in cols if col in df.columns]
    view = df[use_cols].head(n).copy()
    for col in view.columns:
        view[col] = view[col].map(lambda value: "" if pd.isna(value) else str(value))
    return "\n".join(
        [
            "| " + " | ".join(view.columns) + " |",
            "| " + " | ".join(["---"] * len(view.columns)) + " |",
            *["| " + " | ".join(row) + " |" for row in view.to_numpy(dtype=str)],
        ]
    )


def collect_bio_best(out_root: Path) -> pd.DataFrame:
    rows = []
    for path in sorted((out_root / "biology_sidecars").glob("*/benchmark_summary.csv")):
        try:
            df = pd.read_csv(path)
        except Exception:
            continue
        if df.empty:
            continue
        df.insert(0, "sidecar", path.parent.name)
        rows.append(df)
    if not rows:
        return pd.DataFrame()
    merged = pd.concat(rows, ignore_index=True)
    sort_cols = []
    if "roc_auc" in merged.columns:
        merged["_auc_sort"] = pd.to_numeric(merged["roc_auc"], errors="coerce").fillna(-1.0)
        sort_cols.append("_auc_sort")
    if "balanced_accuracy" in merged.columns:
        merged["_bal_sort"] = pd.to_numeric(merged["balanced_accuracy"], errors="coerce").fillna(-1.0)
        sort_cols.append("_bal_sort")
    if "mae_weeks" in merged.columns:
        merged["_mae_sort"] = pd.to_numeric(merged["mae_weeks"], errors="coerce").fillna(999.0)
    if sort_cols:
        merged = merged.sort_values(sort_cols + ["_mae_sort"], ascending=[False] * len(sort_cols) + [True])
    return merged.drop(columns=[col for col in ["_auc_sort", "_bal_sort", "_mae_sort"] if col in merged.columns])


def write_report(
    out_root: Path,
    group_summary: pd.DataFrame,
    final_summary: pd.DataFrame,
    lodo: pd.DataFrame,
    random_df: pd.DataFrame,
    final_rows: list[dict[str, Any]],
    bio_rows: list[dict[str, Any]],
    feature_summary: dict[str, Any] | None,
    args: argparse.Namespace,
) -> None:
    best = final_summary.iloc[0].to_dict() if not final_summary.empty else {}
    beats_v20 = False
    if best:
        beats_v20 = (
            float(best.get("lodo_mean_mae_weeks") or 999.0) <= V20_BEST["lodo_mean_mae_weeks"]
            and float(best.get("lodo_cr_auc") or 0.0) >= V20_BEST["lodo_cr_auc"]
        )
    bio_best = collect_bio_best(out_root)
    lines = [
        "# v23 Cloud Max Batch Report",
        "",
        f"Date: {utc_now()}",
        "",
        "## Scope",
        "",
        "Robust cloud batch using the approved RTX 5090 host. The run expands v22 age DL beyond residual MLP width sweeps and runs CPU sidecar biological-signal screens without touching raw FASTQ.",
        "",
        "## Execution",
        "",
        f"- host: `{socket.gethostname()}`",
        f"- python: `{PYTHON}`",
        f"- matrix: `{args.matrix}`",
        f"- output: `{out_root}`",
        f"- group configs requested: `{args.max_configs}`",
        f"- group/lodo GPU workers: `{args.group_workers}` / `{args.lodo_workers}`",
        "",
        "## Selected Candidate",
        "",
        "`not_available`" if not best else f"- config: `{best.get('config_id')}`",
        "" if not best else f"- architecture: `{best.get('architecture')}`",
        "" if not best else f"- GroupKFold r/MAE/R2: `{best.get('pearson_r')}` / `{best.get('mae_weeks')}` / `{best.get('r2')}`",
        "" if not best else f"- LODO mean/worst MAE: `{best.get('lodo_mean_mae_weeks')}` / `{best.get('lodo_worst_mae_weeks')}`",
        "" if not best else f"- LODO CR AUC: `{best.get('lodo_cr_auc')}`",
        "" if not best else f"- old104 weighted MAE: `{best.get('old_104w_weighted_mae_weeks')}`",
        "" if not best else f"- beats v20 on both LODO mean MAE and CR AUC: `{beats_v20}`",
        "",
        "## Comparators",
        "",
        f"- v20: `{V20_BEST}`",
        f"- v22: `{V22_BEST}`",
        "",
        "## Top DL Configs",
        "",
        md_table(
            final_summary if not final_summary.empty else group_summary,
            [
                "config_id",
                "status",
                "architecture",
                "n_feature_prefilter",
                "preprocess",
                "width",
                "depth",
                "dropout",
                "lr",
                "pearson_r",
                "mae_weeks",
                "r2",
                "cr_detection_auc",
                "lodo_mean_mae_weeks",
                "lodo_worst_mae_weeks",
                "lodo_cr_auc",
                "old_104w_weighted_mae_weeks",
                "final_score",
            ],
            n=16,
        ),
        "",
        "## Random-Label Sanity",
        "",
        md_table(random_df, ["config_id", "status", "pearson_r", "mae_weeks", "r2", "cr_detection_auc"], n=8),
        "",
        "## Biological Sidecars",
        "",
        f"- sidecar jobs: `{len(bio_rows)}`",
        md_table(
            bio_best,
            ["sidecar", "target_id", "model_name", "cv_strategy", "n_samples", "balanced_accuracy", "roc_auc", "mae_weeks", "pearson_r"],
            n=16,
        ),
        "",
        "## Feature Interpretation",
        "",
        "`not_run`" if feature_summary is None else f"`{feature_summary}`",
        "",
        "## Outputs",
        "",
        f"- `{out_root / 'v23_group_summary.csv'}`",
        f"- `{out_root / 'v23_lodo_summary.csv'}`",
        f"- `{out_root / 'v23_final_summary.csv'}`",
        f"- `{out_root / 'v23_random_label_summary.csv'}`",
        f"- `{out_root / 'biology_sidecars'}`",
        f"- `{out_root / 'feature_interpretation'}`",
    ]
    text = "\n".join(line for line in lines if line is not None) + "\n"
    (out_root / "v23_cloud_max_batch_report.md").write_text(text, encoding="utf-8")
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(text, encoding="utf-8")


def heartbeat_loop(out_root: Path, state: dict[str, Any], stop_event: threading.Event) -> None:
    while not stop_event.is_set():
        write_json(
            out_root / "heartbeat.json",
            {
                "timestamp": utc_now(),
                "pid": os.getpid(),
                "host": socket.gethostname(),
                **state,
            },
        )
        stop_event.wait(30)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--target-registry", type=Path, default=DEFAULT_TARGET_REGISTRY)
    parser.add_argument("--out-root", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--max-configs", type=int, default=72)
    parser.add_argument("--top-lodo", type=int, default=8)
    parser.add_argument("--random-top", type=int, default=2)
    parser.add_argument("--epochs", type=int, default=300)
    parser.add_argument("--patience", type=int, default=30)
    parser.add_argument("--group-workers", type=int, default=2)
    parser.add_argument("--lodo-workers", type=int, default=2)
    parser.add_argument("--max-retries", type=int, default=1)
    parser.add_argument("--feature-top-n", type=int, default=200)
    parser.add_argument("--skip-bio-sidecars", action="store_true")
    parser.add_argument("--skip-feature-interpretation", action="store_true")
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()

    args.out_root.mkdir(parents=True, exist_ok=True)
    event_log = args.out_root / "v23_event_log.jsonl"
    state = {"phase": "init", "started": utc_now()}
    stop_event = threading.Event()
    heartbeat = threading.Thread(target=heartbeat_loop, args=(args.out_root, state, stop_event), daemon=True)
    heartbeat.start()
    started = time.time()

    try:
        configs = candidate_configs()
        if args.smoke:
            configs = configs[:2]
            args.top_lodo = min(args.top_lodo, 1)
            args.random_top = min(args.random_top, 1)
            args.epochs = min(args.epochs, 5)
            args.patience = min(args.patience, 2)
        else:
            configs = configs[: args.max_configs]
        write_json(
            args.out_root / "v23_run_manifest.json",
            {
                "timestamp": utc_now(),
                "root": str(ROOT),
                "python": PYTHON,
                "matrix": str(args.matrix),
                "metadata": str(args.metadata),
                "target_registry": str(args.target_registry),
                "out_root": str(args.out_root),
                "configs": configs,
                "args": vars(args),
                "comparators": {"v20": V20_BEST, "v22": V22_BEST},
            },
        )

        state["phase"] = "biology_sidecars_and_groupkfold"
        bio_future = None
        bio_rows: list[dict[str, Any]] = []
        with ThreadPoolExecutor(max_workers=1) as side_pool:
            bio_future = side_pool.submit(run_bio_sidecars, args, args.out_root, event_log)
            group_tasks, config_lookup = make_group_tasks(configs, args.out_root)
            pd.DataFrame(group_tasks).to_json(args.out_root / "v23_group_task_manifest.json", orient="records", indent=2)
            group_rows = run_tasks_parallel(
                group_tasks,
                matrix=args.matrix,
                metadata=args.metadata,
                epochs=args.epochs,
                patience=args.patience,
                workers=args.group_workers,
                max_retries=args.max_retries,
                event_log=event_log,
            )
            group_summary = rows_to_group_summary(group_rows)
            group_summary.to_csv(args.out_root / "v23_group_summary.csv", index=False)
            group_summary.to_json(args.out_root / "v23_group_summary.json", orient="records", indent=2)
            bio_rows = bio_future.result()

        state["phase"] = "lodo"
        lodo_tasks = make_lodo_tasks(group_summary[group_summary["status"].ne("failed")], config_lookup, args.out_root, args.top_lodo)
        pd.DataFrame(lodo_tasks).to_json(args.out_root / "v23_lodo_task_manifest.json", orient="records", indent=2)
        lodo_rows_raw = run_tasks_parallel(
            lodo_tasks,
            matrix=args.matrix,
            metadata=args.metadata,
            epochs=args.epochs,
            patience=args.patience,
            workers=args.lodo_workers,
            max_retries=args.max_retries,
            event_log=event_log,
        )
        lodo_rows = add_old104_metrics(lodo_rows_raw)
        final_summary, lodo = attach_lodo_summary(group_summary, lodo_rows)
        lodo.to_csv(args.out_root / "v23_lodo_summary.csv", index=False)
        final_summary.to_csv(args.out_root / "v23_final_summary.csv", index=False)
        final_summary.to_json(args.out_root / "v23_final_summary.json", orient="records", indent=2)

        state["phase"] = "random_label_sanity"
        random_tasks = make_random_tasks(final_summary[final_summary["status"].ne("failed")], config_lookup, args.out_root, args.random_top)
        random_rows = run_tasks_parallel(
            random_tasks,
            matrix=args.matrix,
            metadata=args.metadata,
            epochs=args.epochs,
            patience=args.patience,
            workers=min(args.lodo_workers, max(1, args.random_top)),
            max_retries=args.max_retries,
            event_log=event_log,
        )
        random_df = pd.DataFrame(random_rows)
        random_df.to_csv(args.out_root / "v23_random_label_summary.csv", index=False)

        state["phase"] = "final_fit"
        final_tasks = make_final_task(final_summary[final_summary["status"].ne("failed")], config_lookup, args.out_root)
        final_rows = run_tasks_parallel(
            final_tasks,
            matrix=args.matrix,
            metadata=args.metadata,
            epochs=args.epochs,
            patience=args.patience,
            workers=1,
            max_retries=args.max_retries,
            event_log=event_log,
        )
        pd.DataFrame(final_rows).to_csv(args.out_root / "v23_final_fit_summary.csv", index=False)

        state["phase"] = "feature_interpretation"
        feature_summary = run_feature_interpretation(args, args.out_root, final_rows, event_log)

        state["phase"] = "report"
        write_report(args.out_root, group_summary, final_summary, lodo, random_df, final_rows, bio_rows, feature_summary, args)
        write_json(
            args.out_root / "v23_completion_state.json",
            {
                "status": "completed",
                "timestamp": utc_now(),
                "elapsed_sec": round(time.time() - started, 1),
                "best_config": None if final_summary.empty else final_summary.iloc[0].to_dict(),
                "n_group_rows": int(len(group_summary)),
                "n_lodo_rows": int(len(lodo)),
                "n_random_rows": int(len(random_df)),
                "n_bio_sidecars": int(len(bio_rows)),
                "feature_summary": feature_summary,
            },
        )
        state["phase"] = "completed"
        print(f"[Done] elapsed_sec={round(time.time() - started, 1)} out={args.out_root}")
        return 0
    finally:
        stop_event.set()
        heartbeat.join(timeout=2)
        write_json(args.out_root / "heartbeat.json", {"timestamp": utc_now(), "pid": os.getpid(), **state})


if __name__ == "__main__":
    raise SystemExit(main())

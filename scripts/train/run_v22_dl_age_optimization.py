#!/usr/bin/env python3
"""v22 narrow DL optimization for age-anchored biological methylation signals."""
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


LOCAL_ROOT = Path("/home/zdq-as/mouse_methyl_work")
REMOTE_ROOT = Path("/root/autodl-tmp/mouse_methyl_work")
ROOT = REMOTE_ROOT if REMOTE_ROOT.exists() else LOCAL_ROOT
DEFAULT_MATRIX = ROOT / "results" / "multidataset_v8_3_ablation" / "all6" / "all_rrbs_region_matrix_5kb.parquet"
DEFAULT_METADATA = ROOT / "metadata" / "model_sample_metadata_v8.csv"
DEFAULT_OUT = ROOT / "results" / "v22_dl_optimization"
REPORT = ROOT / "doc" / "20_analysis" / "57_20260521_v22_dl_optimization_report.md"
TRAIN_SCRIPT = ROOT / "scripts" / "train" / "train_deep_clock.py"
BENCHMARK_SCRIPT = ROOT / "scripts" / "validate" / "benchmark_metrics.py"
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


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str), encoding="utf-8")


def slug(value: Any) -> str:
    return str(value).replace(".", "p").replace("/", "-")


def config_id(cfg: dict[str, Any], idx: int) -> str:
    return "_".join(
        [
            f"v22_{idx:03d}",
            cfg["architecture"],
            f"top{cfg['n_feature_prefilter']}",
            cfg["preprocess"],
            f"w{cfg['width']}",
            f"d{cfg['depth']}",
            f"do{slug(cfg['dropout'])}",
            f"lr{slug(cfg['lr'])}",
        ]
    )


def candidate_configs() -> list[dict[str, Any]]:
    base = {
        "architecture": "res_mlp",
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
    rows: list[dict[str, Any]] = []
    for top in ["500", "1000", "1500", "2000"]:
        for width in [192, 256, 384]:
            for depth in [3, 4]:
                rows.append(
                    {
                        **base,
                        "n_feature_prefilter": top,
                        "preprocess": "robust",
                        "width": width,
                        "depth": depth,
                        "dropout": 0.1,
                        "lr": 3e-4,
                    }
                )
    for top in ["500", "1000", "2000"]:
        for preprocess in ["standard", "quantile_uniform"]:
            rows.append(
                {
                    **base,
                    "n_feature_prefilter": top,
                    "preprocess": preprocess,
                    "width": 256,
                    "depth": 3,
                    "dropout": 0.1,
                    "lr": 3e-4,
                }
            )
            rows.append(
                {
                    **base,
                    "n_feature_prefilter": top,
                    "preprocess": preprocess,
                    "width": 256,
                    "depth": 4,
                    "dropout": 0.1,
                    "lr": 1e-3,
                }
            )
    # A small regularization sweep around the prior v20 winner.
    for dropout in [0.05, 0.15, 0.25]:
        rows.append(
            {
                **base,
                "n_feature_prefilter": "1000",
                "preprocess": "robust",
                "width": 256,
                "depth": 3,
                "dropout": dropout,
                "lr": 3e-4,
            }
        )
    deduped = []
    seen = set()
    for row in rows:
        key = tuple(sorted(row.items()))
        if key not in seen:
            seen.add(key)
            deduped.append(row)
    return deduped


def result_matches(path: Path, cfg: dict[str, Any], eval_mode: str, epochs: int) -> bool:
    if not path.exists():
        return False
    try:
        result = read_json(path)
    except Exception:
        return False
    training = result.get("training") or {}
    return (
        str(result.get("architecture")) == str(cfg["architecture"])
        and str(result.get("n_feature_prefilter")) == str(cfg["n_feature_prefilter"])
        and str(result.get("preprocess")) == str(cfg["preprocess"])
        and str(result.get("eval_mode")) == eval_mode
        and int(training.get("epochs") or -1) == int(epochs)
        and float(training.get("lr") or -1.0) == float(cfg["lr"])
    )


def run(cmd: list[str], cwd: Path = ROOT) -> None:
    print(" ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=cwd, check=True)


def train_cmd(
    cfg: dict[str, Any],
    matrix: Path,
    metadata: Path,
    out_dir: Path,
    eval_mode: str,
    *,
    heldout_dataset: str | None = None,
    seed: int = 42,
    randomize: bool = False,
    epochs: int | None = None,
    patience: int | None = None,
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
        cfg["architecture"],
        "--n_feature_prefilter",
        str(cfg["n_feature_prefilter"]),
        "--imputation",
        cfg["imputation"],
        "--preprocess",
        cfg["preprocess"],
        "--min_train_feature_presence",
        str(cfg["min_train_feature_presence"]),
        "--min_train_group_feature_presence",
        str(cfg["min_train_group_feature_presence"]),
        "--max_train_dataset_mean_shift",
        str(cfg["max_train_dataset_mean_shift"]),
        "--min_train_agebin_feature_presence",
        str(cfg["min_train_agebin_feature_presence"]),
        "--target_transform",
        cfg["target_transform"],
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
        str(epochs or cfg["epochs"]),
        "--patience",
        str(patience or cfg["patience"]),
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
    ]
    if cfg.get("amp"):
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


def run_group(idx: int, cfg: dict[str, Any], matrix: Path, metadata: Path, out_root: Path, epochs: int, patience: int) -> dict[str, Any]:
    cid = config_id(cfg, idx)
    out = out_root / "groupkfold" / cid
    out.mkdir(parents=True, exist_ok=True)
    write_json(out / "config.json", cfg)
    if not result_matches(out / "result.json", cfg, "groupkfold", epochs):
        try:
            run(train_cmd(cfg, matrix, metadata, out, "groupkfold", seed=42, epochs=epochs, patience=patience))
        except subprocess.CalledProcessError as exc:
            write_json(out / "failed.json", {"cmd": exc.cmd, "returncode": exc.returncode, "config": cfg, "timestamp": utc_now()})
            return {"config_id": cid, "status": "failed", **cfg, "group_score": -1.0, "error": f"returncode={exc.returncode}"}
    metrics = benchmark(out / "predictions.csv", out / "benchmark_result.json")
    result = read_json(out / "result.json")
    return {
        "config_id": cid,
        "status": "completed",
        **cfg,
        **{k: metrics.get(k) for k in ["n_samples", "pearson_r", "mae_weeks", "medae_weeks", "rmse_weeks", "r2", "cross_dataset_mae", "cr_detection_auc"]},
        "exec_time_sec": result.get("exec_time_sec"),
        "group_score": group_score(metrics),
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


def compact_dict(payload: dict[str, Any] | None, keys: list[str]) -> dict[str, Any] | None:
    if payload is None:
        return None
    return {key: payload.get(key) for key in keys if key in payload}


def compact_final_metrics(payload: dict[str, Any] | None) -> dict[str, Any] | None:
    if payload is None:
        return None
    view = compact_dict(
        payload,
        [
            "architecture",
            "cuda_device",
            "torch_version",
            "n_samples",
            "n_features_total",
            "n_feature_prefilter",
            "pearson_r",
            "mae_weeks",
            "medae_weeks",
            "rmse_weeks",
            "r2",
            "exec_time_sec",
            "warning",
        ],
    )
    folds = payload.get("fold_rows") or []
    if folds and isinstance(folds[0], dict):
        view["n_features_selected"] = folds[0].get("n_features_selected")
        view["epochs_ran"] = folds[0].get("epochs_ran")
    return view


def write_report(out_root: Path, summary: pd.DataFrame, lodo: pd.DataFrame, random_metrics: dict[str, Any] | None, final_metrics: dict[str, Any] | None) -> None:
    best = summary.iloc[0].to_dict() if not summary.empty else {}
    beats_v20 = False
    if best:
        beats_v20 = (
            float(best.get("lodo_mean_mae_weeks") or 999.0) <= V20_BEST["lodo_mean_mae_weeks"]
            and float(best.get("lodo_cr_auc") or 0.0) >= V20_BEST["lodo_cr_auc"]
        )
    summary_cols = [
        "config_id",
        "status",
        "n_feature_prefilter",
        "preprocess",
        "width",
        "depth",
        "dropout",
        "pearson_r",
        "mae_weeks",
        "r2",
        "cr_detection_auc",
        "lodo_mean_mae_weeks",
        "lodo_worst_mae_weeks",
        "lodo_cr_auc",
        "old_104w_weighted_mae_weeks",
        "final_score",
    ]
    lodo_cols = [
        "config_id",
        "dataset",
        "n_samples",
        "pearson_r",
        "mae_weeks",
        "r2",
        "cr_detection_auc",
        "old_104w_n",
        "old_104w_mae_weeks",
    ]
    random_view = compact_dict(
        random_metrics,
        [
            "pearson_r",
            "mae_weeks",
            "medae_weeks",
            "rmse_weeks",
            "r2",
            "cr_detection_auc",
            "cr_cohens_d",
            "cross_dataset_mae",
            "n_samples",
        ],
    )
    final_view = compact_final_metrics(final_metrics)
    lines = [
        "# v22 DL Optimization Report",
        "",
        f"Date: {utc_now()}",
        "",
        "## Scope",
        "",
        "Narrow RTX 5090 DL optimization around the prior v20 res_mlp winner. The biological readout is age prediction plus CR detection from held-out age predictions; no CNN/Transformer expansion was used.",
        "",
        "## Selected Candidate",
        "",
        "`not_available`" if not best else f"- config: `{best.get('config_id')}`",
        "" if not best else f"- GroupKFold r/MAE/R2: `{best.get('pearson_r')}` / `{best.get('mae_weeks')}` / `{best.get('r2')}`",
        "" if not best else f"- LODO mean/worst MAE: `{best.get('lodo_mean_mae_weeks')}` / `{best.get('lodo_worst_mae_weeks')}`",
        "" if not best else f"- LODO CR AUC: `{best.get('lodo_cr_auc')}`",
        "" if not best else f"- old104 weighted MAE: `{best.get('old_104w_weighted_mae_weeks')}`",
        "" if not best else f"- beats v20 on both LODO mean MAE and CR AUC: `{beats_v20}`",
        "",
        "## v20 Comparator",
        "",
        f"`{V20_BEST}`",
        "",
        "## Random-Label Sanity",
        "",
        "`not_run`" if random_view is None else f"`{random_view}`",
        "",
        "## Final All-Data Fit",
        "",
        "`not_run`" if final_view is None else f"`{final_view}`",
        "",
        "## Top Configs",
        "",
        md_table(summary[[col for col in summary_cols if col in summary.columns]].head(12)),
        "",
        "## LODO Rows",
        "",
        md_table(lodo[[col for col in lodo_cols if col in lodo.columns]]),
        "",
        "## Outputs",
        "",
        f"- `{out_root / 'v22_dl_summary.csv'}`",
        f"- `{out_root / 'v22_lodo_summary.csv'}`",
        f"- `{out_root / 'random_label_sanity'}`",
        f"- `{out_root / 'final_all_data_model'}`",
    ]
    text = "\n".join(line for line in lines if line is not None) + "\n"
    (out_root / "v22_dl_optimization_report.md").write_text(text, encoding="utf-8")
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--out-root", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--max-configs", type=int, default=12)
    parser.add_argument("--top-lodo", type=int, default=4)
    parser.add_argument("--epochs", type=int, default=300)
    parser.add_argument("--patience", type=int, default=30)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--skip-random-label", action="store_true")
    parser.add_argument("--skip-final-fit", action="store_true")
    args = parser.parse_args()

    started = time.time()
    args.out_root.mkdir(parents=True, exist_ok=True)
    configs = candidate_configs()
    if args.smoke:
        configs = configs[:2]
        args.epochs = min(args.epochs, 5)
        args.patience = min(args.patience, 2)
        args.top_lodo = min(args.top_lodo, 1)
    elif args.max_configs > 0:
        configs = configs[: args.max_configs]

    write_json(
        args.out_root / "v22_dl_run_manifest.json",
        {
            "timestamp": utc_now(),
            "python": PYTHON,
            "matrix": str(args.matrix),
            "metadata": str(args.metadata),
            "out_root": str(args.out_root),
            "n_configs": len(configs),
            "epochs": args.epochs,
            "patience": args.patience,
            "configs": configs,
        },
    )

    rows: list[dict[str, Any]] = []
    id_to_config: dict[str, dict[str, Any]] = {}
    for idx, cfg in enumerate(configs, start=1):
        cid = config_id(cfg, idx)
        id_to_config[cid] = cfg
        rows.append(run_group(idx, cfg, args.matrix, args.metadata, args.out_root, args.epochs, args.patience))
        pd.DataFrame(rows).sort_values(["group_score", "pearson_r"], ascending=False).to_csv(
            args.out_root / "v22_dl_summary_partial.csv", index=False
        )

    summary = pd.DataFrame(rows).sort_values(["group_score", "pearson_r"], ascending=False).reset_index(drop=True)
    lodo_rows: list[dict[str, Any]] = []
    for cid in summary.head(args.top_lodo)["config_id"].tolist():
        cfg = id_to_config[cid]
        config_lodo = []
        for dataset in DATASETS:
            out = args.out_root / "lodo" / cid / dataset
            if not result_matches(out / "result.json", cfg, "heldout", args.epochs):
                run(train_cmd(cfg, args.matrix, args.metadata, out, "heldout", heldout_dataset=dataset, seed=42, epochs=args.epochs, patience=args.patience))
            metrics = benchmark(out / "predictions.csv", out / "benchmark_result.json")
            pred = pd.read_csv(out / "predictions.csv")
            old = pred[pred["age_weeks_true"] >= 104]
            row = {
                "config_id": cid,
                "dataset": dataset,
                **{k: metrics.get(k) for k in ["n_samples", "pearson_r", "mae_weeks", "rmse_weeks", "r2", "cr_detection_auc"]},
                "old_104w_n": int(len(old)),
                "old_104w_mae_weeks": None
                if old.empty
                else round(float(np.mean(np.abs(old["age_weeks_true"] - old["age_weeks_pred"]))), 3),
            }
            lodo_rows.append(row)
            config_lodo.append(row)
        lodo_summary = summarize_lodo(config_lodo)
        mask = summary["config_id"].eq(cid)
        for key, value in lodo_summary.items():
            summary.loc[mask, key] = value
        merged = summary.loc[mask].iloc[0].to_dict()
        summary.loc[mask, "lodo_score"] = lodo_score({**merged, **lodo_summary})

    summary["lodo_score"] = summary["lodo_score"].fillna(-1.0)
    summary["final_score"] = summary["group_score"] + summary["lodo_score"].clip(lower=0.0)
    summary = summary.sort_values(["final_score", "group_score"], ascending=False).reset_index(drop=True)

    random_metrics = None
    final_metrics = None
    if not summary.empty:
        best_id = str(summary.iloc[0]["config_id"])
        best_cfg = id_to_config[best_id]
        if not args.skip_random_label:
            out = args.out_root / "random_label_sanity" / best_id
            if not result_matches(out / "result.json", best_cfg, "groupkfold", args.epochs):
                run(train_cmd(best_cfg, args.matrix, args.metadata, out, "groupkfold", seed=99, randomize=True, epochs=args.epochs, patience=args.patience))
            random_metrics = benchmark(out / "predictions.csv", out / "benchmark_result.json")
        if not args.skip_final_fit:
            out = args.out_root / "final_all_data_model"
            if not result_matches(out / "result.json", best_cfg, "final", args.epochs):
                run(train_cmd(best_cfg, args.matrix, args.metadata, out, "final", seed=42, epochs=args.epochs, patience=args.patience))
            final_metrics = read_json(out / "result.json")

    lodo = pd.DataFrame(lodo_rows)
    summary.to_csv(args.out_root / "v22_dl_summary.csv", index=False)
    summary.to_json(args.out_root / "v22_dl_summary.json", orient="records", indent=2)
    lodo.to_csv(args.out_root / "v22_lodo_summary.csv", index=False)
    write_report(args.out_root, summary, lodo, random_metrics, final_metrics)
    print(json.dumps(summary.iloc[0].to_dict() if not summary.empty else {}, indent=2, sort_keys=True, default=str), flush=True)
    print(f"[Done] elapsed_sec={round(time.time() - started, 1)} out={args.out_root}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

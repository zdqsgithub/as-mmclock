#!/usr/bin/env python3
"""v17.1 demo-only narrow high-confidence ML optimization.

This is intentionally not broad autoresearch. It evaluates a small, pre-
registered set of configurations around the best v17/v7.5 LightGBM baseline.
"""
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


ROOT = Path("/home/zdq-as/mouse_methyl_work")
MATRIX = ROOT / "results" / "multidataset_v8_3_ablation" / "all6" / "all_rrbs_region_matrix_5kb.parquet"
METADATA = ROOT / "metadata" / "model_sample_metadata_v8.csv"
OUT_ROOT = ROOT / "results" / "demo_ml_narrow_optimization_v17_1"
REPORT = ROOT / "doc" / "20_analysis" / "49_20260520_demo_ml_narrow_optimization_report.md"
PROJECT_INDEX = ROOT / "doc" / "00_meta" / "02_20260519_project_status_index_v13.md"
DATASETS = ["GSE120137", "GSE80672", "GSE93957", "GSE121141", "GSE60012", "GSE213628"]

CONFIGS = [
    {
        "name": "c01_q_p08_top1000_agebin08",
        "model_type": "lgbm",
        "preprocess": "quantile_uniform",
        "presence": "0.8",
        "group_presence": "0.8",
        "max_shift": "0.15",
        "agebin_presence": "0.8",
        "top_n": "1000",
    },
    {
        "name": "c02_q_p08_top2000_agebin08",
        "model_type": "lgbm",
        "preprocess": "quantile_uniform",
        "presence": "0.8",
        "group_presence": "0.8",
        "max_shift": "0.15",
        "agebin_presence": "0.8",
        "top_n": "2000",
    },
    {
        "name": "c03_q_p08_top500_agebin08",
        "model_type": "lgbm",
        "preprocess": "quantile_uniform",
        "presence": "0.8",
        "group_presence": "0.8",
        "max_shift": "0.15",
        "agebin_presence": "0.8",
        "top_n": "500",
    },
    {
        "name": "c04_robust_p095_top1000_agebin08",
        "model_type": "lgbm",
        "preprocess": "robust",
        "presence": "0.95",
        "group_presence": "0.95",
        "max_shift": "0.15",
        "agebin_presence": "0.8",
        "top_n": "1000",
    },
    {
        "name": "c05_robust_p095_top2000_agebin08",
        "model_type": "lgbm",
        "preprocess": "robust",
        "presence": "0.95",
        "group_presence": "0.95",
        "max_shift": "0.15",
        "agebin_presence": "0.8",
        "top_n": "2000",
    },
    {
        "name": "c06_robust_p08_top1000_agebin08",
        "model_type": "lgbm",
        "preprocess": "robust",
        "presence": "0.8",
        "group_presence": "0.8",
        "max_shift": "0.15",
        "agebin_presence": "0.8",
        "top_n": "1000",
    },
]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def run(cmd: list[str], cwd: Path = ROOT) -> None:
    print(" ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=cwd, check=True)


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str), encoding="utf-8")


def train_cmd(config: dict[str, str], out_dir: Path, randomize: bool = False) -> list[str]:
    cmd = [
        sys.executable,
        "scripts/train/train_clock.py",
        "--matrix_path",
        str(MATRIX),
        "--metadata_path",
        str(METADATA),
        "--feature_type",
        "region",
        "--model_type",
        config["model_type"],
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
        "--min_train_agebin_feature_presence",
        config["agebin_presence"],
        "--target_transform",
        "log1p_days",
        "--output_dir",
        str(out_dir),
    ]
    if randomize:
        cmd.append("--randomize_labels")
    return cmd


def heldout_cmd(config: dict[str, str], dataset: str, out_dir: Path) -> list[str]:
    return [
        sys.executable,
        "scripts/train/train_heldout_clock.py",
        "--train_matrix",
        str(MATRIX),
        "--test_matrix",
        str(MATRIX),
        "--train_datasets",
        f"all_except:{dataset}",
        "--test_datasets",
        dataset,
        "--metadata_path",
        str(METADATA),
        "--model_type",
        config["model_type"],
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
        "--min_train_agebin_feature_presence",
        config["agebin_presence"],
        "--target_transform",
        "log1p_days",
        "--output_dir",
        str(out_dir),
    ]


def benchmark(predictions: Path, out_path: Path) -> dict[str, Any]:
    if not out_path.exists():
        run(
            [
                sys.executable,
                "scripts/validate/benchmark_metrics.py",
                "--predictions",
                str(predictions),
                "--out",
                str(out_path),
            ]
        )
    return read_json(out_path)


def metric_block(df: pd.DataFrame) -> dict[str, Any]:
    true = df["age_weeks_true"].to_numpy(float)
    pred = df["age_weeks_pred"].to_numpy(float)
    if len(df) < 3 or np.std(true) == 0 or np.std(pred) == 0:
        r = 0.0
    else:
        r = float(np.corrcoef(true, pred)[0, 1])
    ss_res = float(np.sum((true - pred) ** 2))
    ss_tot = float(np.sum((true - true.mean()) ** 2))
    return {
        "n_samples": int(len(df)),
        "pearson_r": round(r, 4),
        "mae_weeks": round(float(np.mean(np.abs(true - pred))), 3),
        "rmse_weeks": round(float(np.sqrt(np.mean((true - pred) ** 2))), 3),
        "r2": None if ss_tot <= 1e-8 else round(float(1 - ss_res / ss_tot), 4),
    }


def support_metrics(predictions: Path, out_csv: Path | None = None) -> dict[str, Any]:
    pred = pd.read_csv(predictions)
    meta = pd.read_csv(METADATA)
    rows = []
    for _, row in pred.iterrows():
        dataset = str(row["dataset_batch"])
        tissue = str(row.get("tissue", "unknown"))
        train = meta[(meta["dataset_batch"].astype(str) != dataset) & meta["age_weeks"].notna()].copy()
        same = train[train["tissue"].astype(str).eq(tissue)]
        same_n = int(len(same))
        max_age = float(same["age_weeks"].max()) if same_n else np.nan
        gap = float(row["age_weeks_true"]) - max_age if same_n else np.inf
        support_covered = same_n >= 10 and gap <= 8.0
        stress_group = ""
        if dataset == "GSE121141" and float(row["age_weeks_true"]) >= 104:
            stress_group = "gse121141_old104_target_tissue_stress"
        rows.append(
            {
                **row.to_dict(),
                "same_tissue_train_n": same_n,
                "same_tissue_train_max_age_weeks": None if not np.isfinite(max_age) else round(max_age, 3),
                "age_support_gap_weeks": None if not np.isfinite(gap) else round(gap, 3),
                "support_covered": bool(support_covered),
                "stress_test_group": stress_group,
            }
        )
    ann = pd.DataFrame(rows)
    if out_csv is not None:
        ann.to_csv(out_csv, index=False)
    support = ann[ann["support_covered"]]
    unsupported = ann[~ann["support_covered"]]
    old104 = ann[ann["stress_test_group"].eq("gse121141_old104_target_tissue_stress")]
    return {
        "support_covered": metric_block(support) if not support.empty else None,
        "unsupported": metric_block(unsupported) if not unsupported.empty else None,
        "gse121141_old104_stress": metric_block(old104) if not old104.empty else None,
    }


def random_pass(metrics: dict[str, Any], real_mae: float) -> bool:
    return abs(float(metrics.get("pearson_r", 999))) < 0.2 and float(metrics.get("mae_weeks", 0)) > real_mae


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


def write_report(summary: dict[str, Any], config_table: pd.DataFrame, lodo: pd.DataFrame) -> None:
    best = summary["selected_config"]
    lines = [
        "# v17.1 Demo ML Narrow Optimization Report",
        "",
        f"Date: {utc_now()}",
        "",
        "## Scope",
        "",
        "This is a narrow, high-confidence, demonstration-only ML optimization. It",
        "does not change the v13 Route C model scope and is not broad autoresearch.",
        "",
        "## Selected Configuration",
        "",
        f"- config: `{best['name']}`",
        f"- preprocess: `{best['preprocess']}`",
        f"- top regions: `{best['top_n']}`",
        f"- support-covered MAE: `{summary['selected_metrics']['support_covered_mae_weeks']}` weeks",
        f"- GroupKFold MAE: `{summary['selected_metrics']['groupkfold_mae_weeks']}` weeks",
        f"- random-label pass: `{summary['selected_metrics']['random_label_pass']}`",
        "",
        "## Config Comparison",
        "",
        md_table(
            config_table[
                [
                    "config",
                    "preprocess",
                    "top_n",
                    "groupkfold_mae_weeks",
                    "groupkfold_r",
                    "support_covered_mae_weeks",
                    "unsupported_mae_weeks",
                    "gse121141_old104_mae_weeks",
                    "cr_auc",
                ]
            ]
        ),
        "",
        "## Selected LODO",
        "",
        md_table(lodo),
        "",
        "## Decision",
        "",
        "The optimization produced a demo-selected ML baseline, but the strict",
        "held-out metrics remain far from publication-grade clock thresholds.",
        "The result is useful for demonstration and experiment planning, not for",
        "new full-lifespan old target-tissue claims.",
        "",
        "## Guardrails",
        "",
        "- no FASTQ download;",
        "- no Bismark;",
        "- no deep learning;",
        "- no human clock CpG mapping;",
        "- no dummy AUC;",
        "- no broad autoresearch beyond the six pre-registered configs.",
    ]
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    (OUT_ROOT / "v17_1_demo_ml_narrow_optimization_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def update_index() -> None:
    link = "- v17.1 demo ML narrow optimization report: `doc/20_analysis/49_20260520_demo_ml_narrow_optimization_report.md`"
    if not PROJECT_INDEX.exists():
        return
    text = PROJECT_INDEX.read_text(encoding="utf-8")
    if link in text:
        return
    marker = "- v17 demo ML best baseline report: `doc/20_analysis/48_20260520_demo_ml_best_baseline_report.md`"
    if marker in text:
        text = text.replace(marker, marker + "\n" + link)
    else:
        text += "\n" + link + "\n"
    PROJECT_INDEX.write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-root", default=str(OUT_ROOT))
    args = parser.parse_args()
    out_root = Path(args.out_root)
    out_root.mkdir(parents=True, exist_ok=True)
    started = time.time()

    manifest = {
        "timestamp": utc_now(),
        "demo_only": True,
        "narrow_high_confidence": True,
        "interpreter": sys.executable,
        "interpreter_note": "Compatible AS-DS-Ops venv used because the project .venv was re-created during session setup.",
        "matrix": str(MATRIX),
        "metadata": str(METADATA),
        "configs": CONFIGS,
        "download_authorized": False,
        "bismark_authorized": False,
        "autoresearch_authorized": False,
    }
    write_json(out_root / "run_manifest.json", manifest)

    rows: list[dict[str, Any]] = []
    for config in CONFIGS:
        config_dir = out_root / "configs" / config["name"]
        group_dir = config_dir / "groupkfold"
        if not (group_dir / "result.json").exists():
            run(train_cmd(config, group_dir))
        bm = benchmark(group_dir / "predictions.csv", group_dir / "benchmark_result.json")
        support = support_metrics(group_dir / "predictions.csv")
        rows.append(
            {
                "config": config["name"],
                "preprocess": config["preprocess"],
                "top_n": config["top_n"],
                "presence": config["presence"],
                "group_presence": config["group_presence"],
                "agebin_presence": config["agebin_presence"],
                "groupkfold_mae_weeks": bm.get("mae_weeks"),
                "groupkfold_r": bm.get("pearson_r"),
                "groupkfold_r2": bm.get("r2"),
                "support_covered_mae_weeks": support["support_covered"]["mae_weeks"],
                "unsupported_mae_weeks": support["unsupported"]["mae_weeks"],
                "gse121141_old104_mae_weeks": support["gse121141_old104_stress"]["mae_weeks"],
                "cr_auc": bm.get("cr_detection_auc"),
            }
        )

    config_table = pd.DataFrame(rows).sort_values(
        ["support_covered_mae_weeks", "groupkfold_mae_weeks", "gse121141_old104_mae_weeks"],
        ascending=[True, True, True],
    )
    config_table.to_csv(out_root / "config_comparison.csv", index=False)

    selected = None
    selected_random = None
    for _, row in config_table.head(3).iterrows():
        config = next(item for item in CONFIGS if item["name"] == row["config"])
        random_dir = out_root / "configs" / config["name"] / "random_label_sanity"
        if not (random_dir / "result.json").exists():
            run(train_cmd(config, random_dir, randomize=True))
        random_bm = benchmark(random_dir / "predictions.csv", random_dir / "benchmark_result.json")
        if random_pass(random_bm, float(row["groupkfold_mae_weeks"])):
            selected = config
            selected_random = random_bm
            break
    if selected is None:
        selected = next(item for item in CONFIGS if item["name"] == config_table.iloc[0]["config"])
        selected_random = read_json(out_root / "configs" / selected["name"] / "random_label_sanity" / "benchmark_result.json")

    selected_group_dir = out_root / "configs" / selected["name"] / "groupkfold"
    selected_support = support_metrics(selected_group_dir / "predictions.csv", out_root / "support_annotations_selected.csv")
    selected_bm = benchmark(selected_group_dir / "predictions.csv", selected_group_dir / "benchmark_result.json")

    lodo_rows = []
    for dataset in DATASETS:
        heldout_dir = out_root / "selected_lodo" / dataset
        if not (heldout_dir / "result.json").exists():
            run(heldout_cmd(selected, dataset, heldout_dir))
        bm = benchmark(heldout_dir / "predictions.csv", heldout_dir / "benchmark_result.json")
        pred = pd.read_csv(heldout_dir / "predictions.csv")
        old = pred[pred["age_weeks_true"] >= 104]
        lodo_rows.append(
            {
                "dataset": dataset,
                "n_samples": bm.get("n_samples"),
                "pearson_r": bm.get("pearson_r"),
                "mae_weeks": bm.get("mae_weeks"),
                "rmse_weeks": bm.get("rmse_weeks"),
                "r2": bm.get("r2"),
                "cr_detection_auc": bm.get("cr_detection_auc"),
                "old_104w_n": int(len(old)),
                "old_104w_mae_weeks": None if old.empty else round(float(np.mean(np.abs(old["age_weeks_true"] - old["age_weeks_pred"]))), 3),
            }
        )
    lodo = pd.DataFrame(lodo_rows)
    lodo.to_csv(out_root / "selected_lodo_summary.csv", index=False)

    selected_metrics = {
        "groupkfold_mae_weeks": selected_bm.get("mae_weeks"),
        "groupkfold_r": selected_bm.get("pearson_r"),
        "support_covered_mae_weeks": selected_support["support_covered"]["mae_weeks"],
        "unsupported_mae_weeks": selected_support["unsupported"]["mae_weeks"],
        "gse121141_old104_mae_weeks": selected_support["gse121141_old104_stress"]["mae_weeks"],
        "random_label_r": selected_random.get("pearson_r"),
        "random_label_mae_weeks": selected_random.get("mae_weeks"),
        "random_label_pass": random_pass(selected_random, float(selected_bm.get("mae_weeks"))),
    }
    summary = {
        "timestamp": utc_now(),
        "demo_only": True,
        "selected_config": selected,
        "selected_metrics": selected_metrics,
        "config_comparison": rows,
        "lodo": lodo_rows,
        "exec_time_sec": round(time.time() - started, 1),
        "guardrails": {
            "download_authorized": False,
            "bismark_authorized": False,
            "deep_learning_authorized": False,
            "headline_claim_changed": False,
        },
    }
    write_json(out_root / "v17_1_narrow_optimization_summary.json", summary)
    pd.DataFrame([{"selected_config": selected["name"], **selected_metrics}]).to_csv(out_root / "selected_summary.csv", index=False)
    write_report(summary, config_table, lodo)
    update_index()
    print(json.dumps(summary, indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

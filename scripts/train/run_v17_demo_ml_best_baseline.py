#!/usr/bin/env python3
"""Run the best locked ML baseline for demonstration-only full-data testing.

This runner does not perform autoresearch. It reuses the documented best ML
configuration from v7.5/v8.3 and writes outputs to a v17 demo-only directory.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import joblib
import numpy as np
import pandas as pd


ROOT = Path("/home/zdq-as/mouse_methyl_work")
OUT_ROOT = ROOT / "results" / "demo_ml_best_baseline_v17"
REPORT = ROOT / "doc" / "20_analysis" / "48_20260520_demo_ml_best_baseline_report.md"
PROJECT_INDEX = ROOT / "doc" / "00_meta" / "02_20260519_project_status_index_v13.md"

DEFAULT_MATRIX = ROOT / "results" / "multidataset_v8_3_ablation" / "all6" / "all_rrbs_region_matrix_5kb.parquet"
DEFAULT_METADATA = ROOT / "metadata" / "model_sample_metadata_v8.csv"
DATASETS = ["GSE120137", "GSE80672", "GSE93957", "GSE121141", "GSE60012", "GSE213628"]

CONFIG = {
    "name": "v75_lgbm_quantile_p08_top1000_agebin08",
    "model_type": "lgbm",
    "preprocess": "quantile_uniform",
    "imputation": "median",
    "min_train_feature_presence": "0.8",
    "min_train_group_feature_presence": "0.8",
    "max_train_dataset_mean_shift": "0.15",
    "min_train_agebin_feature_presence": "0.8",
    "n_feature_prefilter": "1000",
    "target_transform": "log1p_days",
}


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


def import_train_clock():
    path = ROOT / "scripts" / "train" / "train_clock.py"
    spec = importlib.util.spec_from_file_location("demo_train_clock_helpers", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


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


def train_clock_cmd(matrix: Path, metadata: Path, out_dir: Path, randomize: bool = False) -> list[str]:
    cmd = [
        sys.executable,
        "scripts/train/train_clock.py",
        "--matrix_path",
        str(matrix),
        "--metadata_path",
        str(metadata),
        "--feature_type",
        "region",
        "--model_type",
        CONFIG["model_type"],
        "--n_feature_prefilter",
        CONFIG["n_feature_prefilter"],
        "--imputation",
        CONFIG["imputation"],
        "--preprocess",
        CONFIG["preprocess"],
        "--min_train_feature_presence",
        CONFIG["min_train_feature_presence"],
        "--min_train_group_feature_presence",
        CONFIG["min_train_group_feature_presence"],
        "--max_train_dataset_mean_shift",
        CONFIG["max_train_dataset_mean_shift"],
        "--min_train_agebin_feature_presence",
        CONFIG["min_train_agebin_feature_presence"],
        "--target_transform",
        CONFIG["target_transform"],
        "--output_dir",
        str(out_dir),
    ]
    if randomize:
        cmd.append("--randomize_labels")
    return cmd


def heldout_cmd(matrix: Path, metadata: Path, test_dataset: str, out_dir: Path) -> list[str]:
    return [
        sys.executable,
        "scripts/train/train_heldout_clock.py",
        "--train_matrix",
        str(matrix),
        "--test_matrix",
        str(matrix),
        "--train_datasets",
        f"all_except:{test_dataset}",
        "--test_datasets",
        test_dataset,
        "--metadata_path",
        str(metadata),
        "--model_type",
        CONFIG["model_type"],
        "--n_feature_prefilter",
        CONFIG["n_feature_prefilter"],
        "--imputation",
        CONFIG["imputation"],
        "--preprocess",
        CONFIG["preprocess"],
        "--min_train_feature_presence",
        CONFIG["min_train_feature_presence"],
        "--min_train_group_feature_presence",
        CONFIG["min_train_group_feature_presence"],
        "--max_train_dataset_mean_shift",
        CONFIG["max_train_dataset_mean_shift"],
        "--min_train_agebin_feature_presence",
        CONFIG["min_train_agebin_feature_presence"],
        "--target_transform",
        CONFIG["target_transform"],
        "--output_dir",
        str(out_dir),
    ]


def metric_block(df: pd.DataFrame) -> dict[str, Any]:
    true = df["age_weeks_true"].to_numpy(float)
    pred = df["age_weeks_pred"].to_numpy(float)
    if len(df) < 3 or np.std(true) == 0 or np.std(pred) == 0:
        r = 0.0
    else:
        r = float(np.corrcoef(true, pred)[0, 1])
    ss_res = float(np.sum((true - pred) ** 2))
    ss_tot = float(np.sum((true - true.mean()) ** 2))
    r2 = None if ss_tot <= 1e-8 else round(float(1 - ss_res / ss_tot), 4)
    return {
        "n_samples": int(len(df)),
        "pearson_r": round(r, 4),
        "mae_weeks": round(float(np.mean(np.abs(true - pred))), 3),
        "medae_weeks": round(float(np.median(np.abs(true - pred))), 3),
        "rmse_weeks": round(float(np.sqrt(np.mean((true - pred) ** 2))), 3),
        "r2": r2,
    }


def annotate_support(predictions: Path, metadata: Path, out_csv: Path) -> dict[str, Any]:
    pred = pd.read_csv(predictions)
    meta = pd.read_csv(metadata)
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
                "support_reason": "same_tissue_n>=10_and_age_gap<=8w" if support_covered else "unsupported_same_tissue_or_age_range",
                "stress_test_group": stress_group,
            }
        )
    annotated = pd.DataFrame(rows)
    annotated.to_csv(out_csv, index=False)
    support = annotated[annotated["support_covered"]]
    unsupported = annotated[~annotated["support_covered"]]
    old104 = annotated[annotated["stress_test_group"].eq("gse121141_old104_target_tissue_stress")]
    return {
        "support_covered": metric_block(support) if not support.empty else None,
        "unsupported": metric_block(unsupported) if not unsupported.empty else None,
        "gse121141_old104_stress": metric_block(old104) if not old104.empty else None,
    }


def fit_final_all_data_model(matrix: Path, metadata: Path, out_dir: Path) -> dict[str, Any]:
    helper = import_train_clock()
    beta, meta = helper.load_clock_data(matrix, metadata)
    X_raw = beta.T.values.astype(np.float32)
    y_days = meta["age_days"].values.astype(np.float32)
    y_fit = helper.build_target(y_days, CONFIG["target_transform"])
    groups = meta["dataset_batch"].fillna("unknown").astype(str).values

    presence = np.isfinite(X_raw).mean(axis=0) >= float(CONFIG["min_train_feature_presence"])
    stability = helper.train_dataset_stability_mask(
        X_raw,
        groups,
        float(CONFIG["min_train_group_feature_presence"]),
        float(CONFIG["max_train_dataset_mean_shift"]),
    )
    agebin = helper.train_agebin_presence_mask(
        X_raw,
        y_days / 7,
        float(CONFIG["min_train_agebin_feature_presence"]),
        8,
    )
    feature_mask = presence & stability & agebin
    presence_idx = np.flatnonzero(feature_mask)
    if len(presence_idx) == 0:
        raise RuntimeError("No features pass final all-data filters.")

    top_local_idx = helper.select_features(X_raw[:, presence_idx], y_fit, CONFIG["n_feature_prefilter"])
    top_idx = presence_idx[top_local_idx]
    pipe = helper.Pipeline(
        [
            *helper.make_preprocess_steps(CONFIG["preprocess"], CONFIG["imputation"], 42, X_raw.shape[0]),
            ("model", helper.make_model(CONFIG["model_type"], 42)),
        ]
    )
    pipe.fit(X_raw[:, top_idx], y_fit)
    pred_fit = pipe.predict(X_raw[:, top_idx])
    y_pred_days, y_pred_weeks = helper.inverse_target(pred_fit, CONFIG["target_transform"], float(y_days.max()))
    metrics = helper.compute_metrics(y_days / 7, y_pred_weeks)

    out_dir.mkdir(parents=True, exist_ok=True)
    model_payload = {
        "demo_only": True,
        "config": CONFIG,
        "matrix_path": str(matrix),
        "metadata_path": str(metadata),
        "selected_feature_indices": top_idx,
        "selected_feature_ids": beta.index[top_idx].astype(str).tolist(),
        "pipeline": pipe,
    }
    joblib.dump(model_payload, out_dir / "final_all_data_lgbm_model.joblib")
    pd.DataFrame({"feature_id": beta.index[top_idx].astype(str), "feature_index": top_idx}).to_csv(
        out_dir / "selected_features.csv",
        index=False,
    )
    pred_df = pd.DataFrame(
        {
            "sample_id": meta["sample_id"].values,
            "dataset_batch": meta.get("dataset_batch", pd.Series(["unknown"] * len(meta))).values,
            "tissue": meta.get("tissue", pd.Series(["unknown"] * len(meta))).values,
            "intervention": meta.get("intervention", pd.Series(["control"] * len(meta))).values,
            "age_days_true": y_days,
            "age_days_pred": y_pred_days,
            "age_weeks_true": y_days / 7,
            "age_weeks_pred": y_pred_weeks,
            "residual_weeks": y_days / 7 - y_pred_weeks,
        }
    )
    pred_df.to_csv(out_dir / "apparent_train_predictions.csv", index=False)
    payload = {
        "demo_only": True,
        "warning": "Apparent all-data fit metrics are optimistic and not a held-out benchmark.",
        "n_samples": int(X_raw.shape[0]),
        "n_features_total": int(X_raw.shape[1]),
        "n_features_after_presence": int(presence.sum()),
        "n_features_after_stability": int(stability.sum()),
        "n_features_after_agebin": int(agebin.sum()),
        "n_features_after_all_filters": int(feature_mask.sum()),
        "n_features_selected": int(len(top_idx)),
        **metrics,
    }
    write_json(out_dir / "apparent_train_metrics.json", payload)
    return payload


def write_report(out_root: Path, summary: dict[str, Any], lodo: pd.DataFrame, support: dict[str, Any]) -> None:
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

    lines = [
        "# v17 Demo ML Best Baseline Report",
        "",
        f"Date: {utc_now()}",
        "",
        "## Scope",
        "",
        "This is a demonstration-only ML baseline. It does not change the accepted",
        "v13 Route C model scope and does not claim full-lifespan old target-tissue",
        "generalization.",
        "",
        "## Locked Configuration",
        "",
        "- model: `lgbm`",
        "- preprocessing: `quantile_uniform`",
        "- feature filters: presence `0.8`, group presence `0.8`, dataset mean shift `0.15`, age-bin presence `0.8`",
        "- feature selector: fold-internal age correlation",
        "- top regions: `1000`",
        "- target transform: `log1p_days`",
        "- matrix: `results/multidataset_v8_3_ablation/all6/all_rrbs_region_matrix_5kb.parquet`",
        "",
        "## Main Metrics",
        "",
        f"- GroupKFold MAE: `{summary['groupkfold'].get('mae_weeks')}` weeks",
        f"- GroupKFold r: `{summary['groupkfold'].get('pearson_r')}`",
        f"- GroupKFold R2: `{summary['groupkfold'].get('r2')}`",
        f"- Random-label sanity r: `{summary['random_label'].get('pearson_r')}`",
        f"- Random-label sanity MAE: `{summary['random_label'].get('mae_weeks')}` weeks",
        "",
        "## Support-Covered Split",
        "",
        f"- support-covered: `{support.get('support_covered')}`",
        f"- unsupported: `{support.get('unsupported')}`",
        f"- GSE121141 old104+ stress: `{support.get('gse121141_old104_stress')}`",
        "",
        "## Leave-One-Dataset-Out",
        "",
        md_table(lodo),
        "",
        "## Full-Data Demo Fit",
        "",
        "The final all-data model is saved only for demonstration. Its apparent",
        "training metrics are optimistic and must not be used as held-out evidence.",
        "",
        f"- apparent all-data MAE: `{summary['apparent_all_data'].get('mae_weeks')}` weeks",
        f"- apparent all-data r: `{summary['apparent_all_data'].get('pearson_r')}`",
        f"- selected features: `{summary['apparent_all_data'].get('n_features_selected')}`",
        "",
        "## Guardrails",
        "",
        "- demo-only training was explicitly requested by the user;",
        "- no autoresearch was run;",
        "- no FASTQ download or Bismark was run;",
        "- no human clock CpG mapping was used;",
        "- no dummy AUC was generated;",
        "- this result is for demonstration and experiment-planning only.",
        "",
        "## Outputs",
        "",
        f"- `{out_root.relative_to(ROOT)}/groupkfold/predictions.csv`",
        f"- `{out_root.relative_to(ROOT)}/support_annotations.csv`",
        f"- `{out_root.relative_to(ROOT)}/lodo_summary.csv`",
        f"- `{out_root.relative_to(ROOT)}/final_all_data_model/final_all_data_lgbm_model.joblib`",
    ]
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    (out_root / "demo_ml_best_baseline_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def update_index() -> None:
    if not PROJECT_INDEX.exists():
        return
    link = "- v17 demo ML best baseline report: `doc/20_analysis/48_20260520_demo_ml_best_baseline_report.md`"
    text = PROJECT_INDEX.read_text(encoding="utf-8")
    if link in text:
        return
    marker = "- current data exhaustion and Route A vs deep learning assessment: `doc/20_analysis/47_20260519_current_data_exhaustion_route_a_vs_deep_learning_assessment.md`"
    if marker in text:
        text = text.replace(marker, marker + "\n" + link)
    else:
        text += "\n" + link + "\n"
    PROJECT_INDEX.write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", default=str(DEFAULT_MATRIX))
    parser.add_argument("--metadata", default=str(DEFAULT_METADATA))
    parser.add_argument("--out-root", default=str(OUT_ROOT))
    args = parser.parse_args()

    matrix = Path(args.matrix)
    metadata = Path(args.metadata)
    out_root = Path(args.out_root)
    out_root.mkdir(parents=True, exist_ok=True)
    started = time.time()

    manifest = {
        "timestamp": utc_now(),
        "demo_only": True,
        "training_authorized_by_user_request": True,
        "download_authorized": False,
        "bismark_authorized": False,
        "autoresearch_authorized": False,
        "matrix": str(matrix),
        "metadata": str(metadata),
        "config": CONFIG,
    }
    write_json(out_root / "demo_run_manifest.json", manifest)

    group_dir = out_root / "groupkfold"
    if not (group_dir / "result.json").exists():
        run(train_clock_cmd(matrix, metadata, group_dir))
    group_metrics = benchmark(group_dir / "predictions.csv", group_dir / "benchmark_result.json")

    random_dir = out_root / "random_label_sanity"
    if not (random_dir / "result.json").exists():
        run(train_clock_cmd(matrix, metadata, random_dir, randomize=True))
    random_metrics = benchmark(random_dir / "predictions.csv", random_dir / "benchmark_result.json")

    support = annotate_support(group_dir / "predictions.csv", metadata, out_root / "support_annotations.csv")
    write_json(out_root / "support_metrics.json", support)

    lodo_rows = []
    for dataset in DATASETS:
        heldout_dir = out_root / "lodo" / dataset
        if not (heldout_dir / "result.json").exists():
            run(heldout_cmd(matrix, metadata, dataset, heldout_dir))
        bm = benchmark(heldout_dir / "predictions.csv", heldout_dir / "benchmark_result.json")
        row = {"dataset": dataset, **{key: bm.get(key) for key in ["n_samples", "pearson_r", "mae_weeks", "rmse_weeks", "r2", "cr_detection_auc"]}}
        pred = pd.read_csv(heldout_dir / "predictions.csv")
        old = pred[pred["age_weeks_true"] >= 104]
        row["old_104w_n"] = int(len(old))
        row["old_104w_mae_weeks"] = None if old.empty else round(float(np.mean(np.abs(old["age_weeks_true"] - old["age_weeks_pred"]))), 3)
        lodo_rows.append(row)
    lodo = pd.DataFrame(lodo_rows)
    lodo.to_csv(out_root / "lodo_summary.csv", index=False)

    apparent = fit_final_all_data_model(matrix, metadata, out_root / "final_all_data_model")
    summary = {
        "timestamp": utc_now(),
        "demo_only": True,
        "config": CONFIG,
        "groupkfold": group_metrics,
        "random_label": random_metrics,
        "support_metrics": support,
        "lodo": lodo_rows,
        "apparent_all_data": apparent,
        "exec_time_sec": round(time.time() - started, 1),
        "guardrails": {
            "autoresearch_authorized": False,
            "download_authorized": False,
            "bismark_authorized": False,
            "headline_claim_changed": False,
        },
    }
    write_json(out_root / "demo_ml_best_baseline_summary.json", summary)
    pd.DataFrame(
        [
            {
                "scope": "groupkfold_all6",
                **{key: group_metrics.get(key) for key in ["n_samples", "pearson_r", "mae_weeks", "rmse_weeks", "r2", "cross_dataset_mae", "cr_detection_auc"]},
            },
            {
                "scope": "random_label_all6",
                **{key: random_metrics.get(key) for key in ["n_samples", "pearson_r", "mae_weeks", "rmse_weeks", "r2", "cross_dataset_mae", "cr_detection_auc"]},
            },
            {
                "scope": "apparent_all_data_fit",
                **{key: apparent.get(key) for key in ["n_samples", "pearson_r", "mae_weeks", "rmse_weeks", "r2"]},
            },
        ]
    ).to_csv(out_root / "demo_ml_best_baseline_summary.csv", index=False)
    write_report(out_root, summary, lodo, support)
    update_index()
    print(json.dumps(summary, indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Run v8.6 leakage-safe calibration POC on existing v8.3 LODO predictions.

The calibration parameters are fitted only from non-target LODO prediction
residuals. GSE121141 labels are never used to fit calibration parameters.
This is a diagnostic POC and does not train a new clock model.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")
DEFAULT_V8_3_ROOT = ROOT / "results" / "benchmark_v8_3_ablation"
DEFAULT_OUT_DIR = ROOT / "results" / "v8_6_calibration_poc"
DEFAULT_REPORT = ROOT / "doc" / "20_analysis" / "21_20260518_v8_6_leakage_safe_calibration_poc_report.md"
TARGET_DATASET = "GSE121141"
OLD_THRESHOLD_WEEKS = 104.0
MIN_GROUP_N = 5

PRED_AGE_BINS = [0, 13, 26, 52, 78, 104, 130, np.inf]
PRED_AGE_LABELS = ["0-13w", "13-26w", "26-52w", "52-78w", "78-104w", "104-130w", "130w+"]


def md_table(df: pd.DataFrame, max_rows: int | None = None) -> str:
    if df.empty:
        return "No rows."
    view = df.copy()
    if max_rows is not None:
        view = view.head(max_rows)
    for col in view.columns:
        view[col] = view[col].map(lambda value: "" if pd.isna(value) else str(value))
    header = "| " + " | ".join(view.columns) + " |"
    sep = "| " + " | ".join(["---"] * len(view.columns)) + " |"
    rows = ["| " + " | ".join(row) + " |" for row in view.to_numpy(dtype=str)]
    return "\n".join([header, sep, *rows])


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def safe_read(path: Path, **kwargs) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, **kwargs)


def selected_variant_configs(v8_3_root: Path) -> list[tuple[str, str]]:
    summary = safe_read(v8_3_root / "variant_summary.tsv", sep="\t")
    if summary.empty:
        rows = []
        for path in sorted(v8_3_root.glob("*/*/lodo/*/predictions.csv")):
            rows.append((path.parents[2].name, path.parents[1].name))
        return sorted(set(rows))
    selected = summary[summary["selected_best_config"].fillna(False).astype(bool)]
    return list(selected[["variant", "config"]].itertuples(index=False, name=None))


def load_lodo_predictions(v8_3_root: Path, variant: str, config: str) -> pd.DataFrame:
    frames = []
    for path in sorted((v8_3_root / variant / config / "lodo").glob("*/predictions.csv")):
        dataset = path.parent.name
        frame = pd.read_csv(path)
        frame.insert(0, "test_dataset", dataset)
        frame.insert(0, "config", config)
        frame.insert(0, "variant", variant)
        frames.append(frame)
    if not frames:
        return pd.DataFrame()
    data = pd.concat(frames, ignore_index=True)
    data["tissue"] = data["tissue"].fillna("unknown").astype(str)
    data["residual_weeks"] = data["age_weeks_true"].astype(float) - data["age_weeks_pred"].astype(float)
    data["pred_age_bin"] = pd.cut(
        data["age_weeks_pred"].astype(float).clip(lower=0),
        bins=PRED_AGE_BINS,
        labels=PRED_AGE_LABELS,
        include_lowest=True,
    ).astype(str)
    return data


def metric_block(df: pd.DataFrame, pred_col: str) -> dict:
    if df.empty:
        return {"n": 0, "mae_weeks": np.nan, "median_ae_weeks": np.nan, "bias_weeks": np.nan, "rmse_weeks": np.nan}
    true = df["age_weeks_true"].astype(float).to_numpy()
    pred = df[pred_col].astype(float).to_numpy()
    residual = true - pred
    return {
        "n": int(len(df)),
        "mae_weeks": round(float(np.mean(np.abs(residual))), 3),
        "median_ae_weeks": round(float(np.median(np.abs(residual))), 3),
        "bias_weeks": round(float(np.mean(residual)), 3),
        "rmse_weeks": round(float(np.sqrt(np.mean(residual**2))), 3),
    }


def group_offsets(calibration: pd.DataFrame, keys: list[str], min_n: int) -> pd.DataFrame:
    if calibration.empty:
        return pd.DataFrame(columns=keys + ["n", "offset_weeks"])
    grouped = (
        calibration.groupby(keys, dropna=False)
        .agg(n=("sample_id", "count"), offset_weeks=("residual_weeks", "median"))
        .reset_index()
    )
    grouped = grouped[grouped["n"].ge(min_n)].copy()
    grouped["offset_weeks"] = grouped["offset_weeks"].astype(float)
    return grouped


def fit_calibrator(
    calibration: pd.DataFrame,
    method: str,
    target_tissues: set[str],
    scope: str,
    min_group_n: int = MIN_GROUP_N,
) -> dict:
    pool = calibration.copy()
    if scope == "shared_target_tissues_only":
        pool = pool[pool["tissue"].isin(target_tissues)].copy()
    if pool.empty:
        pool = calibration.copy()
        scope = "fallback_all_non_target"
    global_offset = float(pool["residual_weeks"].median()) if not pool.empty else 0.0
    payload = {
        "method": method,
        "scope": scope,
        "global_offset_weeks": global_offset,
        "n_calibration_samples": int(len(pool)),
        "source_datasets": ",".join(sorted(pool["test_dataset"].dropna().astype(str).unique())) if not pool.empty else "",
        "target_labels_used": False,
        "offset_tables": {},
    }
    if method in {"tissue_offset", "tissue_predbin_offset"}:
        payload["offset_tables"]["tissue"] = group_offsets(pool, ["tissue"], min_group_n)
    if method in {"predbin_offset", "tissue_predbin_offset"}:
        payload["offset_tables"]["pred_age_bin"] = group_offsets(pool, ["pred_age_bin"], min_group_n)
    if method == "tissue_predbin_offset":
        payload["offset_tables"]["tissue_pred_age_bin"] = group_offsets(pool, ["tissue", "pred_age_bin"], min_group_n)
    return payload


def lookup_offset(row: pd.Series, fit: dict) -> tuple[float, str]:
    method = fit["method"]
    global_offset = float(fit["global_offset_weeks"])
    tables = fit.get("offset_tables", {})
    if method == "global_offset":
        return global_offset, "global"
    if method == "tissue_offset":
        table = tables.get("tissue", pd.DataFrame())
        match = table[table["tissue"].astype(str).eq(str(row["tissue"]))] if not table.empty else pd.DataFrame()
        if not match.empty:
            return float(match.iloc[0]["offset_weeks"]), "tissue"
        return global_offset, "global_fallback"
    if method == "predbin_offset":
        table = tables.get("pred_age_bin", pd.DataFrame())
        match = table[table["pred_age_bin"].astype(str).eq(str(row["pred_age_bin"]))] if not table.empty else pd.DataFrame()
        if not match.empty:
            return float(match.iloc[0]["offset_weeks"]), "pred_age_bin"
        return global_offset, "global_fallback"
    if method == "tissue_predbin_offset":
        table = tables.get("tissue_pred_age_bin", pd.DataFrame())
        if not table.empty:
            match = table[
                table["tissue"].astype(str).eq(str(row["tissue"]))
                & table["pred_age_bin"].astype(str).eq(str(row["pred_age_bin"]))
            ]
            if not match.empty:
                return float(match.iloc[0]["offset_weeks"]), "tissue_pred_age_bin"
        tissue_table = tables.get("tissue", pd.DataFrame())
        tissue_match = tissue_table[tissue_table["tissue"].astype(str).eq(str(row["tissue"]))] if not tissue_table.empty else pd.DataFrame()
        if not tissue_match.empty:
            return float(tissue_match.iloc[0]["offset_weeks"]), "tissue_fallback"
        bin_table = tables.get("pred_age_bin", pd.DataFrame())
        bin_match = bin_table[bin_table["pred_age_bin"].astype(str).eq(str(row["pred_age_bin"]))] if not bin_table.empty else pd.DataFrame()
        if not bin_match.empty:
            return float(bin_match.iloc[0]["offset_weeks"]), "pred_age_bin_fallback"
        return global_offset, "global_fallback"
    raise ValueError(f"Unknown method: {method}")


def apply_calibrator(target: pd.DataFrame, fit: dict) -> pd.DataFrame:
    out = target.copy()
    offsets = []
    sources = []
    for _, row in out.iterrows():
        offset, source = lookup_offset(row, fit)
        offsets.append(offset)
        sources.append(source)
    out["calibration_method"] = fit["method"]
    out["calibration_scope"] = fit["scope"]
    out["calibration_offset_weeks"] = offsets
    out["calibration_offset_source"] = sources
    out["age_weeks_pred_calibrated"] = out["age_weeks_pred"].astype(float) + out["calibration_offset_weeks"].astype(float)
    out["residual_weeks_calibrated"] = out["age_weeks_true"].astype(float) - out["age_weeks_pred_calibrated"].astype(float)
    out["target_labels_used_for_calibration"] = False
    return out


def offset_rows(variant: str, config: str, fit: dict) -> list[dict]:
    rows = [
        {
            "variant": variant,
            "config": config,
            "method": fit["method"],
            "scope": fit["scope"],
            "level": "global",
            "key": "all",
            "n": fit["n_calibration_samples"],
            "offset_weeks": round(float(fit["global_offset_weeks"]), 3),
            "source_datasets": fit["source_datasets"],
            "target_labels_used": False,
        }
    ]
    for level, table in fit.get("offset_tables", {}).items():
        if table.empty:
            continue
        key_cols = [col for col in table.columns if col not in {"n", "offset_weeks"}]
        for _, row in table.iterrows():
            rows.append(
                {
                    "variant": variant,
                    "config": config,
                    "method": fit["method"],
                    "scope": fit["scope"],
                    "level": level,
                    "key": "|".join(str(row[col]) for col in key_cols),
                    "n": int(row["n"]),
                    "offset_weeks": round(float(row["offset_weeks"]), 3),
                    "source_datasets": fit["source_datasets"],
                    "target_labels_used": False,
                }
            )
    return rows


def summarize_application(pred: pd.DataFrame, variant: str, config: str, method: str, scope: str) -> list[dict]:
    rows = []
    for subset_name, subset in [
        ("all", pred),
        ("under104", pred[pred["age_weeks_true"].astype(float) < OLD_THRESHOLD_WEEKS]),
        ("old104plus", pred[pred["age_weeks_true"].astype(float) >= OLD_THRESHOLD_WEEKS]),
    ]:
        metrics = metric_block(subset, "age_weeks_pred_calibrated")
        raw = metric_block(subset, "age_weeks_pred")
        rows.append(
            {
                "variant": variant,
                "config": config,
                "method": method,
                "scope": scope,
                "test_dataset": TARGET_DATASET,
                "subset": subset_name,
                **metrics,
                "raw_mae_weeks": raw["mae_weeks"],
                "delta_mae_vs_raw": None if pd.isna(metrics["mae_weeks"]) else round(float(raw["mae_weeks"] - metrics["mae_weeks"]), 3),
                "target_labels_used_for_calibration": False,
            }
        )
    old = pred[pred["age_weeks_true"].astype(float) >= OLD_THRESHOLD_WEEKS]
    for tissue, tissue_df in old.groupby("tissue", dropna=False):
        metrics = metric_block(tissue_df, "age_weeks_pred_calibrated")
        raw = metric_block(tissue_df, "age_weeks_pred")
        rows.append(
            {
                "variant": variant,
                "config": config,
                "method": method,
                "scope": scope,
                "test_dataset": TARGET_DATASET,
                "subset": f"old104plus_tissue:{tissue}",
                **metrics,
                "raw_mae_weeks": raw["mae_weeks"],
                "delta_mae_vs_raw": None if pd.isna(metrics["mae_weeks"]) else round(float(raw["mae_weeks"] - metrics["mae_weeks"]), 3),
                "target_labels_used_for_calibration": False,
            }
        )
    return rows


def raw_summary(target: pd.DataFrame, variant: str, config: str) -> list[dict]:
    pred = target.copy()
    pred["age_weeks_pred_calibrated"] = pred["age_weeks_pred"].astype(float)
    rows = summarize_application(pred, variant, config, "raw", "none")
    for row in rows:
        row["delta_mae_vs_raw"] = 0.0
    return rows


def nontarget_validation(
    data: pd.DataFrame,
    target_tissues: set[str],
    variant: str,
    config: str,
    methods: list[str],
    scopes: list[str],
) -> tuple[list[dict], list[pd.DataFrame]]:
    rows = []
    frames = []
    eval_datasets = [value for value in sorted(data["test_dataset"].unique()) if value != TARGET_DATASET]
    for eval_dataset in eval_datasets:
        eval_df = data[data["test_dataset"].eq(eval_dataset)].copy()
        raw_metrics = metric_block(eval_df.assign(age_weeks_pred_calibrated=eval_df["age_weeks_pred"]), "age_weeks_pred_calibrated")
        for scope in scopes:
            for method in methods:
                calibration = data[
                    (~data["test_dataset"].isin([TARGET_DATASET, eval_dataset]))
                ].copy()
                fit = fit_calibrator(calibration, method, target_tissues, scope)
                applied = apply_calibrator(eval_df, fit)
                metrics = metric_block(applied, "age_weeks_pred_calibrated")
                rows.append(
                    {
                        "variant": variant,
                        "config": config,
                        "eval_dataset": eval_dataset,
                        "method": method,
                        "scope": fit["scope"],
                        **metrics,
                        "raw_mae_weeks": raw_metrics["mae_weeks"],
                        "delta_mae_vs_raw": round(float(raw_metrics["mae_weeks"] - metrics["mae_weeks"]), 3),
                        "calibration_source_datasets": fit["source_datasets"],
                        "target_dataset_excluded_from_fit": True,
                        "eval_dataset_excluded_from_fit": True,
                    }
                )
                frames.append(applied.assign(eval_dataset=eval_dataset))
    return rows, frames


def run_v86(v8_3_root: Path, out_dir: Path) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    methods = ["global_offset", "tissue_offset", "predbin_offset", "tissue_predbin_offset"]
    scopes = ["all_non_target", "shared_target_tissues_only"]
    summary_rows = []
    target_frames = []
    offset_summary_rows = []
    nontarget_rows = []
    nontarget_frames = []
    raw_lodo_frames = []

    for variant, config in selected_variant_configs(v8_3_root):
        data = load_lodo_predictions(v8_3_root, variant, config)
        if data.empty or TARGET_DATASET not in set(data["test_dataset"]):
            continue
        raw_lodo_frames.append(data)
        target = data[data["test_dataset"].eq(TARGET_DATASET)].copy()
        target_tissues = set(target["tissue"].dropna().astype(str).unique())
        calibration_base = data[~data["test_dataset"].eq(TARGET_DATASET)].copy()
        summary_rows.extend(raw_summary(target, variant, config))
        for scope in scopes:
            for method in methods:
                fit = fit_calibrator(calibration_base, method, target_tissues, scope)
                applied = apply_calibrator(target, fit)
                target_frames.append(applied)
                summary_rows.extend(summarize_application(applied, variant, config, method, fit["scope"]))
                offset_summary_rows.extend(offset_rows(variant, config, fit))
        val_rows, val_frames = nontarget_validation(data, target_tissues, variant, config, methods, scopes)
        nontarget_rows.extend(val_rows)
        nontarget_frames.extend(val_frames)

    summary = pd.DataFrame(summary_rows)
    target_predictions = pd.concat(target_frames, ignore_index=True) if target_frames else pd.DataFrame()
    offsets = pd.DataFrame(offset_summary_rows)
    nontarget = pd.DataFrame(nontarget_rows)
    nontarget_predictions = pd.concat(nontarget_frames, ignore_index=True) if nontarget_frames else pd.DataFrame()
    raw_lodo = pd.concat(raw_lodo_frames, ignore_index=True) if raw_lodo_frames else pd.DataFrame()

    summary.to_csv(out_dir / "gse121141_calibration_summary.csv", index=False)
    target_predictions.to_csv(out_dir / "gse121141_calibrated_predictions.csv", index=False)
    offsets.to_csv(out_dir / "calibration_offsets.csv", index=False)
    nontarget.to_csv(out_dir / "nontarget_lodo_calibration_validation.csv", index=False)
    nontarget_predictions.to_csv(out_dir / "nontarget_calibrated_predictions.csv", index=False)
    raw_lodo.to_csv(out_dir / "raw_lodo_predictions_used.csv", index=False)

    best = pd.DataFrame()
    if not summary.empty:
        best = (
            summary[summary["subset"].eq("old104plus")]
            .sort_values(["delta_mae_vs_raw", "mae_weeks"], ascending=[False, True])
            .groupby(["variant", "config"], as_index=False)
            .head(3)
        )
        best.to_csv(out_dir / "best_old104plus_calibration_by_variant.csv", index=False)

    payload = {
        "status": "complete",
        "target_dataset": TARGET_DATASET,
        "target_labels_used_for_calibration": False,
        "n_variant_configs": int(len(selected_variant_configs(v8_3_root))),
        "best_old104plus_delta_mae": None
        if best.empty
        else round(float(best["delta_mae_vs_raw"].max()), 3),
        "outputs": {
            "target_summary": str(out_dir / "gse121141_calibration_summary.csv"),
            "target_predictions": str(out_dir / "gse121141_calibrated_predictions.csv"),
            "offsets": str(out_dir / "calibration_offsets.csv"),
            "nontarget_validation": str(out_dir / "nontarget_lodo_calibration_validation.csv"),
            "best_old104plus": str(out_dir / "best_old104plus_calibration_by_variant.csv"),
        },
    }
    write_json(out_dir / "v8_6_summary.json", payload)
    return payload


def write_report(out_dir: Path, report_path: Path) -> None:
    summary = safe_read(out_dir / "gse121141_calibration_summary.csv")
    nontarget = safe_read(out_dir / "nontarget_lodo_calibration_validation.csv")
    best = safe_read(out_dir / "best_old104plus_calibration_by_variant.csv")

    old_view = (
        summary[summary["subset"].eq("old104plus")]
        .sort_values(["variant", "delta_mae_vs_raw"], ascending=[True, False])
        [
            [
                "variant",
                "config",
                "method",
                "scope",
                "n",
                "mae_weeks",
                "raw_mae_weeks",
                "delta_mae_vs_raw",
                "bias_weeks",
            ]
        ]
        if not summary.empty
        else pd.DataFrame()
    )
    all_view = (
        summary[summary["subset"].eq("all")]
        .sort_values(["variant", "delta_mae_vs_raw"], ascending=[True, False])
        [
            [
                "variant",
                "method",
                "scope",
                "n",
                "mae_weeks",
                "raw_mae_weeks",
                "delta_mae_vs_raw",
                "bias_weeks",
            ]
        ]
        if not summary.empty
        else pd.DataFrame()
    )
    nontarget_view = pd.DataFrame()
    if not nontarget.empty:
        nontarget_view = (
            nontarget.groupby(["variant", "method", "scope"], dropna=False)
            .agg(
                eval_datasets=("eval_dataset", "nunique"),
                mean_delta_mae_vs_raw=("delta_mae_vs_raw", "mean"),
                median_delta_mae_vs_raw=("delta_mae_vs_raw", "median"),
                mean_mae_weeks=("mae_weeks", "mean"),
                mean_raw_mae_weeks=("raw_mae_weeks", "mean"),
            )
            .reset_index()
            .sort_values(["variant", "mean_delta_mae_vs_raw"], ascending=[True, False])
        )
        for col in ["mean_delta_mae_vs_raw", "median_delta_mae_vs_raw", "mean_mae_weeks", "mean_raw_mae_weeks"]:
            nontarget_view[col] = nontarget_view[col].round(3)

    best_delta = None if best.empty else float(best["delta_mae_vs_raw"].max())
    nontarget_mean_best = None
    if not best.empty and not nontarget.empty:
        row = best.sort_values("delta_mae_vs_raw", ascending=False).iloc[0]
        matched = nontarget[
            nontarget["variant"].eq(row["variant"])
            & nontarget["method"].eq(row["method"])
            & nontarget["scope"].eq(row["scope"])
        ]
        if not matched.empty:
            nontarget_mean_best = float(matched["delta_mae_vs_raw"].mean())

    lines = [
        "# v8.6 Leakage-Safe Calibration POC Report",
        "",
        "Date: 2026-05-18",
        "",
        "## Summary",
        "",
        "v8.6 tested calibration methods on existing v8.3 LODO predictions only. It did not retrain models, download data, or start autoresearch.",
        "",
        "- Calibration fit data: non-GSE121141 LODO residuals only.",
        "- Target labels used for calibration: False.",
        "- Methods: global residual offset, tissue offset, predicted-age-bin offset, tissue+predicted-age-bin offset.",
        "- Scopes: all non-target tissues and shared GSE121141 target tissues only.",
        f"- Best GSE121141 old104+ delta MAE: {'' if best_delta is None else round(best_delta, 3)} weeks.",
        f"- Non-target validation mean delta for that best method: {'' if nontarget_mean_best is None else round(nontarget_mean_best, 3)} weeks.",
        "",
        "## GSE121141 Old104+ Results",
        "",
        md_table(old_view, max_rows=36),
        "",
        "## GSE121141 All-Age Results",
        "",
        md_table(all_view, max_rows=24),
        "",
        "## Non-Target LODO Calibration Validation",
        "",
        "For each non-target dataset, calibration was fitted from all other non-target datasets, excluding both GSE121141 and the evaluated dataset.",
        "",
        md_table(nontarget_view, max_rows=40),
        "",
        "## Decision",
        "",
    ]
    if best_delta is not None and best_delta >= 10 and nontarget_mean_best is not None and nontarget_mean_best >= 0:
        lines.append(
            "Calibration shows a potentially transferable signal. Next step should be a strict training-integrated calibration experiment, still without autoresearch."
        )
    else:
        lines.append(
            "Calibration does not yet provide a robust transferable fix. Do not promote it to headline benchmark or run autoresearch. Keep the conclusion that old same-tissue data support remains the main blocker."
        )
    lines.extend(
        [
            "",
            "## Outputs",
            "",
            "- `results/v8_6_calibration_poc/gse121141_calibration_summary.csv`",
            "- `results/v8_6_calibration_poc/gse121141_calibrated_predictions.csv`",
            "- `results/v8_6_calibration_poc/calibration_offsets.csv`",
            "- `results/v8_6_calibration_poc/nontarget_lodo_calibration_validation.csv`",
            "- `results/v8_6_calibration_poc/best_old104plus_calibration_by_variant.csv`",
        ]
    )
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--v8_3_root", default=str(DEFAULT_V8_3_ROOT))
    parser.add_argument("--out_dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument("--report", default=str(DEFAULT_REPORT))
    args = parser.parse_args()
    out_dir = Path(args.out_dir)
    payload = run_v86(Path(args.v8_3_root), out_dir)
    write_report(out_dir, Path(args.report))
    payload["outputs"]["report"] = str(Path(args.report))
    write_json(out_dir / "v8_6_summary.json", payload)
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Apply the v12 benchmark redefinition to existing prediction files.

This script only reads existing predictions and result metadata. It does not
train models, download data, rebuild matrices, or run autoresearch.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.metrics import f1_score, roc_auc_score

ROOT = Path("/home/zdq-as/mouse_methyl_work")
DEFAULT_ROOTS = [
    ROOT / "results" / "benchmark_v8_3_ablation",
    ROOT / "results" / "benchmark_v8_2_prefilter_liftover",
    ROOT / "results" / "benchmark_v7_5_gse121141_harmonization",
    ROOT / "results" / "validation_v7_5_gse121141_all_except",
    ROOT / "results" / "validation_v7_5_gse80672_cr_all_except",
]
OUT_DIR = ROOT / "results" / "benchmark_v12_1_redefined"
REPORT_PATH = ROOT / "doc" / "20_analysis" / "31_20260519_v12_1_redefined_benchmark_report.md"
V13_RFC_PATH = ROOT / "doc" / "20_analysis" / "32_20260519_v13_data_strategy_rfc.md"

TARGET_TISSUES = {"brain_cortex", "cortex", "heart", "lung"}
SUPPORT_MIN_SAME_TISSUE_N = 10
SUPPORT_MAX_AGE_GAP_WEEKS = 8.0


def clean(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    return str(value)


def as_jsonable(value: Any) -> Any:
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    return value


def write_json(path: Path, payload: dict | list) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, default=as_jsonable), encoding="utf-8")


def load_json(path: Path) -> dict:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def safe_path(path_text: str | None) -> Path | None:
    if not path_text:
        return None
    path = Path(path_text)
    if not path.is_absolute():
        path = ROOT / path
    return path


def dataset_filter(meta: pd.DataFrame, selector: str) -> pd.Series:
    selector = selector.strip()
    if selector.startswith("all_except:"):
        excluded = {item.strip() for item in selector.split(":", 1)[1].split(",") if item.strip()}
        return ~meta["dataset_batch"].astype(str).isin(excluded)
    if selector in {"all", "*"}:
        return pd.Series(True, index=meta.index)
    allowed = {item.strip() for item in selector.split(",") if item.strip()}
    return meta["dataset_batch"].astype(str).isin(allowed)


def normalize_tissue(value: Any) -> str:
    text = clean(value).strip().lower().replace(" ", "_").replace("-", "_")
    aliases = {
        "cortex": "brain_cortex",
        "brain_cortex": "brain_cortex",
        "heart": "heart",
        "lung": "lung",
        "liver": "liver",
        "blood": "blood",
        "skeletal_muscle": "skeletal_muscle",
        "muscle": "skeletal_muscle",
        "intestine": "intestine",
    }
    return aliases.get(text, text or "unknown")


def normalize_predictions(pred: pd.DataFrame) -> pd.DataFrame:
    df = pred.copy()
    if {"age_days_true", "age_days_pred"}.issubset(df.columns):
        pass
    elif {"age_weeks_true", "age_weeks_pred"}.issubset(df.columns):
        df["age_days_true"] = df["age_weeks_true"].astype(float) * 7
        df["age_days_pred"] = df["age_weeks_pred"].astype(float) * 7
    elif {"age_days", "age_weeks_pred"}.issubset(df.columns):
        df["age_days_true"] = df["age_days"].astype(float)
        df["age_days_pred"] = df["age_weeks_pred"].astype(float) * 7
    elif {"age_days", "age_days_pred"}.issubset(df.columns):
        df = df.rename(columns={"age_days": "age_days_true"})
    else:
        raise ValueError("Prediction file lacks standard age columns")
    if "age_weeks_true" not in df.columns:
        df["age_weeks_true"] = df["age_days_true"].astype(float) / 7
    if "age_weeks_pred" not in df.columns:
        df["age_weeks_pred"] = df["age_days_pred"].astype(float) / 7
    for col, default in [
        ("sample_id", ""),
        ("dataset_batch", "unknown"),
        ("tissue", "unknown"),
        ("intervention", "control"),
    ]:
        if col not in df.columns:
            df[col] = default
    df["tissue"] = df["tissue"].map(normalize_tissue)
    return df


def safe_pearson(y_true: np.ndarray, y_pred: np.ndarray) -> tuple[float | None, float | None]:
    if len(y_true) < 3 or np.isclose(np.std(y_true), 0.0) or np.isclose(np.std(y_pred), 0.0):
        return None, None
    r, pval = stats.pearsonr(y_true, y_pred)
    if not np.isfinite(r):
        return None, None
    return float(r), float(pval)


def regression_metrics(df: pd.DataFrame) -> dict:
    if df.empty:
        return {
            "n_samples": 0,
            "n_datasets": 0,
            "pearson_r": None,
            "pearson_pval": None,
            "mae_weeks": None,
            "medae_weeks": None,
            "rmse_weeks": None,
            "bias_weeks": None,
            "r2": None,
        }
    y_true = df["age_weeks_true"].astype(float).to_numpy()
    y_pred = df["age_weeks_pred"].astype(float).to_numpy()
    residual = y_true - y_pred
    r, pval = safe_pearson(y_true, y_pred)
    ss_res = float(np.sum((y_true - y_pred) ** 2))
    ss_tot = float(np.sum((y_true - y_true.mean()) ** 2))
    return {
        "n_samples": int(len(df)),
        "n_datasets": int(df["dataset_batch"].nunique()) if "dataset_batch" in df.columns else 0,
        "pearson_r": None if r is None else round(r, 4),
        "pearson_pval": pval,
        "mae_weeks": round(float(np.mean(np.abs(residual))), 3),
        "medae_weeks": round(float(np.median(np.abs(residual))), 3),
        "rmse_weeks": round(float(np.sqrt(np.mean(residual**2))), 3),
        "bias_weeks": round(float(np.mean(residual)), 3),
        "r2": round(1 - ss_res / ss_tot, 4) if ss_tot > 1e-8 else None,
    }


def compute_age_acceleration(df: pd.DataFrame) -> np.ndarray:
    y_true = df["age_days_true"].astype(float).to_numpy()
    y_pred = df["age_days_pred"].astype(float).to_numpy()
    batches = df["dataset_batch"].fillna("unknown").astype(str).to_numpy()
    acceleration = np.full(len(df), np.nan, dtype=float)
    for batch in np.unique(batches):
        mask = batches == batch
        if mask.sum() < 3 or np.std(y_true[mask]) == 0:
            acceleration[mask] = y_pred[mask] - y_true[mask]
            continue
        lr = LinearRegression()
        lr.fit(y_true[mask].reshape(-1, 1), y_pred[mask])
        acceleration[mask] = y_pred[mask] - lr.predict(y_true[mask].reshape(-1, 1))
    return acceleration


def cr_metrics(df: pd.DataFrame) -> dict:
    if df.empty or "intervention" not in df.columns:
        return {
            "cr_detection_auc": None,
            "cr_detection_f1": None,
            "cr_cohens_d": None,
            "cr_mannwhitney_p": None,
            "cr_status": "no_predictions",
        }
    work = df[df["intervention"].astype(str).isin(["CR", "control"])].copy()
    if work.empty or work["intervention"].nunique() < 2:
        return {
            "cr_detection_auc": None,
            "cr_detection_f1": None,
            "cr_cohens_d": None,
            "cr_mannwhitney_p": None,
            "cr_status": "missing_cr_or_control",
        }
    acceleration = compute_age_acceleration(work)
    labels = (work["intervention"].astype(str).to_numpy() == "CR").astype(int)
    valid = np.isfinite(acceleration)
    acceleration = acceleration[valid]
    labels = labels[valid]
    if len(np.unique(labels)) < 2:
        return {
            "cr_detection_auc": None,
            "cr_detection_f1": None,
            "cr_cohens_d": None,
            "cr_mannwhitney_p": None,
            "cr_status": "invalid_acceleration",
        }
    score = -acceleration
    auc = float(roc_auc_score(labels, score))
    pred_binary = (acceleration < 0).astype(int)
    f1 = float(f1_score(labels, pred_binary, zero_division=0))
    cr_accel = acceleration[labels == 1]
    control_accel = acceleration[labels == 0]
    pooled = np.sqrt((cr_accel.std(ddof=1) ** 2 + control_accel.std(ddof=1) ** 2) / 2)
    cohens_d = None if not np.isfinite(pooled) or pooled == 0 else float((control_accel.mean() - cr_accel.mean()) / pooled)
    try:
        _, mw_p = stats.mannwhitneyu(control_accel, cr_accel, alternative="greater")
        mw_p_value = float(mw_p)
    except Exception:
        mw_p_value = None
    return {
        "cr_detection_auc": round(auc, 4),
        "cr_detection_f1": round(f1, 4),
        "cr_cohens_d": round(cohens_d, 4) if cohens_d is not None else None,
        "cr_mannwhitney_p": round(mw_p_value, 6) if mw_p_value is not None else None,
        "cr_status": "computed",
    }


def infer_prediction_kind(path: Path, result: dict) -> str:
    parts = [part.lower() for part in path.parts]
    if "groupkfold" in parts or clean(result.get("cv_strategy")).startswith("GroupKFold"):
        return "groupkfold"
    if "lodo" in parts:
        return "lodo"
    if "train_dataset" in result:
        return "heldout"
    if any(part.startswith("heldout_") for part in parts):
        return "heldout"
    return "unknown"


def infer_train_selector(path: Path, result: dict, pred: pd.DataFrame) -> str | None:
    if "train_dataset" in result:
        return clean(result.get("train_dataset"))
    parts = list(path.parts)
    if "lodo" in parts:
        idx = parts.index("lodo")
        if idx + 1 < len(parts):
            return f"all_except:{parts[idx + 1]}"
    for part in parts:
        lower = part.lower()
        if lower.startswith("heldout_gse"):
            return f"all_except:{part.split('_', 1)[1].upper()}"
    if len(pred["dataset_batch"].dropna().unique()) == 1:
        return f"all_except:{pred['dataset_batch'].dropna().astype(str).iloc[0]}"
    return None


def load_metadata_for_result(result: dict) -> pd.DataFrame:
    metadata_path = safe_path(result.get("metadata_path"))
    if metadata_path is None or not metadata_path.exists():
        return pd.DataFrame()
    meta = pd.read_csv(metadata_path)
    if "age_days" in meta.columns and "age_weeks" not in meta.columns:
        meta["age_weeks"] = pd.to_numeric(meta["age_days"], errors="coerce") / 7
    meta = meta[meta.get("age_weeks", pd.Series(dtype=float)).notna()].copy()
    if "tissue" in meta.columns:
        meta["tissue"] = meta["tissue"].map(normalize_tissue)
    return meta


def train_pool_for_row(
    row: pd.Series,
    pred: pd.DataFrame,
    result: dict,
    kind: str,
    train_selector: str | None,
    metadata_cache: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    if kind == "groupkfold" or train_selector is None:
        return pred[pred["dataset_batch"].astype(str) != str(row["dataset_batch"])].copy()
    metadata_key = clean(result.get("metadata_path"))
    if metadata_key not in metadata_cache:
        metadata_cache[metadata_key] = load_metadata_for_result(result)
    meta = metadata_cache[metadata_key]
    if meta.empty or "dataset_batch" not in meta.columns:
        return pred[pred["dataset_batch"].astype(str) != str(row["dataset_batch"])].copy()
    return meta[dataset_filter(meta, train_selector)].copy()


def annotate_predictions(path: Path, result: dict) -> pd.DataFrame:
    pred = normalize_predictions(pd.read_csv(path))
    kind = infer_prediction_kind(path, result)
    train_selector = infer_train_selector(path, result, pred)
    metadata_cache: dict[str, pd.DataFrame] = {}
    rows = []
    for _, row in pred.iterrows():
        train_pool = train_pool_for_row(row, pred, result, kind, train_selector, metadata_cache)
        tissue = normalize_tissue(row["tissue"])
        same_tissue = train_pool[train_pool.get("tissue", pd.Series(dtype=str)).map(normalize_tissue).eq(tissue)].copy()
        same_n = int(len(same_tissue))
        if same_n:
            age_col = "age_weeks" if "age_weeks" in same_tissue.columns else "age_weeks_true"
            max_age = float(pd.to_numeric(same_tissue[age_col], errors="coerce").max())
        else:
            max_age = np.nan
        age_true = float(row["age_weeks_true"])
        gap = age_true - max_age if np.isfinite(max_age) else np.nan
        support_covered = bool(same_n >= SUPPORT_MIN_SAME_TISSUE_N and (np.isfinite(gap) and gap <= SUPPORT_MAX_AGE_GAP_WEEKS))
        if same_n < SUPPORT_MIN_SAME_TISSUE_N:
            support_reason = "same_tissue_train_n_lt_10"
        elif not np.isfinite(gap):
            support_reason = "same_tissue_train_age_missing"
        elif gap > SUPPORT_MAX_AGE_GAP_WEEKS:
            support_reason = "heldout_age_exceeds_train_same_tissue_max_plus_8w"
        else:
            support_reason = "support_covered"
        dataset = clean(row["dataset_batch"])
        is_old_target = dataset == "GSE121141" and age_true >= 104 and tissue in TARGET_TISSUES
        if is_old_target:
            stress_group = "gse121141_old104_target_tissue"
        elif not support_covered:
            stress_group = "unsupported_other"
        else:
            stress_group = "support_covered"
        rows.append(
            {
                **row.to_dict(),
                "source_prediction": str(path.relative_to(ROOT)),
                "source_dir": str(path.parent.relative_to(ROOT)),
                "prediction_kind": kind,
                "train_selector": train_selector or "groupkfold_inferred_non_same_dataset",
                "same_tissue_train_n": same_n,
                "same_tissue_train_max_age_weeks": None if not np.isfinite(max_age) else round(max_age, 3),
                "age_support_gap_weeks": None if not np.isfinite(gap) else round(gap, 3),
                "support_covered": support_covered,
                "support_reason": support_reason,
                "stress_test_group": stress_group,
            }
        )
    return pd.DataFrame(rows)


def find_prediction_files(roots: list[Path], include_random_labels: bool) -> list[Path]:
    paths: list[Path] = []
    for root in roots:
        if not root.exists():
            continue
        for path in sorted(root.rglob("predictions.csv")):
            lower = str(path).lower()
            if not include_random_labels and ("random_label" in lower or "random_labels" in lower):
                continue
            if "shuffled" in lower:
                continue
            paths.append(path)
    return sorted(dict.fromkeys(paths))


def build_summary(annotations: pd.DataFrame) -> pd.DataFrame:
    rows = []
    subsets = {
        "all_predictions": annotations,
        "support_covered_headline": annotations[annotations["support_covered"].astype(bool)],
        "unsupported_stress": annotations[~annotations["support_covered"].astype(bool)],
        "gse121141_old104_stress": annotations[annotations["stress_test_group"].eq("gse121141_old104_target_tissue")],
    }
    for subset_name, df in subsets.items():
        metric = regression_metrics(df)
        rows.append({"level": "aggregate", "source_prediction": "all", "subset": subset_name, **metric})
    for (source, subset_name), group in annotations.groupby(["source_prediction", "stress_test_group"], dropna=False):
        rows.append({"level": "source_stress_group", "source_prediction": source, "subset": subset_name, **regression_metrics(group)})
    for source, group in annotations.groupby("source_prediction", dropna=False):
        rows.append({"level": "source_support_covered", "source_prediction": source, "subset": "support_covered", **regression_metrics(group[group["support_covered"].astype(bool)])})
        rows.append({"level": "source_unsupported", "source_prediction": source, "subset": "unsupported", **regression_metrics(group[~group["support_covered"].astype(bool)])})
    return pd.DataFrame(rows)


def load_sibling_shuffled(path: Path) -> pd.DataFrame | None:
    shuffled = path.parent / "predictions_shuffled_intervention.csv"
    if not shuffled.exists():
        return None
    return normalize_predictions(pd.read_csv(shuffled))


def build_cr_metrics(annotations: pd.DataFrame) -> dict:
    gse = annotations[annotations["dataset_batch"].astype(str).eq("GSE80672")].copy()
    if gse.empty:
        return {"status": "no_gse80672_predictions", "sources": []}
    rows = []
    for source, group in gse.groupby("source_prediction", dropna=False):
        source_path = ROOT / source
        real_metrics = cr_metrics(group)
        shuffled_df = load_sibling_shuffled(source_path)
        if shuffled_df is None:
            shuffled_metrics = {"shuffled_sanity_status": "not_available_not_recomputed"}
        else:
            shuffled_metrics = {
                **cr_metrics(shuffled_df[shuffled_df["dataset_batch"].astype(str).eq("GSE80672")]),
                "shuffled_sanity_status": "computed_from_existing_file",
            }
        prefixed_shuffled = {
            (key if key == "shuffled_sanity_status" else f"shuffled_{key}"): value
            for key, value in shuffled_metrics.items()
        }
        rows.append(
            {
                "source_prediction": source,
                "n_samples": int(len(group)),
                "support_covered_n": int(group["support_covered"].astype(bool).sum()),
                **real_metrics,
                **prefixed_shuffled,
            }
        )
    best = sorted(
        [row for row in rows if row.get("cr_detection_auc") is not None],
        key=lambda row: row.get("cr_detection_auc", -1),
        reverse=True,
    )
    return {
        "status": "computed",
        "best_source": best[0]["source_prediction"] if best else None,
        "best_cr_detection_auc": best[0]["cr_detection_auc"] if best else None,
        "sources": rows,
    }


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


def write_report(report_path: Path, summary: pd.DataFrame, cr_payload: dict, n_files: int) -> None:
    headline = summary[(summary["level"].eq("aggregate")) & (summary["subset"].eq("support_covered_headline"))]
    unsupported = summary[(summary["level"].eq("aggregate")) & (summary["subset"].eq("unsupported_stress"))]
    old104 = summary[(summary["level"].eq("aggregate")) & (summary["subset"].eq("gse121141_old104_stress"))]
    headline_mae = headline["mae_weeks"].iloc[0] if not headline.empty else None
    unsupported_mae = unsupported["mae_weeks"].iloc[0] if not unsupported.empty else None
    accepted = (
        pd.notna(headline_mae)
        and pd.notna(unsupported_mae)
        and float(headline_mae) < float(unsupported_mae)
    )
    lines = [
        "# v12.1 Redefined Benchmark Report",
        "",
        "Date: 2026-05-19",
        "",
        "## Summary",
        "",
        "v12.1 applied the v12 benchmark redefinition to existing prediction files only. It did not retrain models, download data, rebuild matrices, run FASTQ/Bismark, or start autoresearch.",
        "",
        f"- Prediction files scanned: {n_files}",
        f"- Benchmark redefinition accepted: `{str(bool(accepted)).lower()}`",
        f"- Support-covered headline MAE: `{headline_mae}` weeks",
        f"- Unsupported stress-test MAE: `{unsupported_mae}` weeks",
        "",
        "## Aggregate Metrics",
        "",
        md_table(summary[summary["level"].eq("aggregate")]),
        "",
        "## GSE121141 old104+ Stress Metric",
        "",
        md_table(old104),
        "",
        "## GSE80672 CR Metrics",
        "",
        md_table(pd.DataFrame(cr_payload.get("sources", [])).sort_values("cr_detection_auc", ascending=False, na_position="last"), max_rows=20),
        "",
        "## Interpretation",
        "",
    ]
    if accepted:
        lines.append("Support-covered rows perform better than unsupported stress-test rows. This supports adopting the v12 benchmark redefinition.")
    else:
        lines.append("Support-covered rows do not clearly outperform unsupported stress-test rows. Do not start autoresearch; inspect metadata, matrices, or move to v13 data strategy.")
    lines.extend(
        [
            "",
            "## Outputs",
            "",
            "- `results/benchmark_v12_1_redefined/prediction_support_annotations.csv`",
            "- `results/benchmark_v12_1_redefined/support_covered_headline_metrics.json`",
            "- `results/benchmark_v12_1_redefined/unsupported_stress_test_metrics.json`",
            "- `results/benchmark_v12_1_redefined/gse121141_old104_stress_metrics.json`",
            "- `results/benchmark_v12_1_redefined/gse80672_cr_redefined_metrics.json`",
            "- `results/benchmark_v12_1_redefined/redefined_benchmark_summary.csv`",
        ]
    )
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_v13_rfc(path: Path) -> None:
    lines = [
        "# v13 Data Strategy RFC",
        "",
        "Date: 2026-05-19",
        "",
        "## Summary",
        "",
        "v13 starts only after v12 benchmark redefinition and v12.1 support-covered evaluation. It is a data strategy decision, not model tuning.",
        "",
        "## Route A: Generate or Collaborate",
        "",
        "- Acquire new old bulk brain_cortex/heart/lung RRBS/WGBS data with sample-specific exact age.",
        "- Minimum useful design: same tissue young/mid/old, with old samples at or above 104 weeks.",
        "- This is the preferred route if the project needs a true full-lifespan target-tissue clock.",
        "",
        "## Route B: Minimal FASTQ/Bismark Pilot",
        "",
        "- Only allowed for raw candidates with BioSample/RunInfo-confirmed sample-specific age, bulk target tissue, bisulfite assay, and practical run size.",
        "- Pilot size: 2-3 samples only.",
        "- Pilot goal: verify alignment, methylation extraction, 5kb common regions, and metadata traceability before any full ETL.",
        "",
        "## Route C: Accept Redefined Benchmark",
        "",
        "- Use support-covered chronological benchmark as the headline result.",
        "- Report GSE121141 old104+ as a stress-test/blocker metric.",
        "- Keep biological-age/CR claims limited to real held-out predictions and shuffled-intervention sanity checks.",
        "",
        "## Non-Negotiable Rules",
        "",
        "- No autoresearch until a new data gate passes or the redefined benchmark is formally accepted.",
        "- No dummy AUC.",
        "- No human clock CpG mapping.",
        "- No raw FASTQ download without explicit minimal ETL approval.",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out_dir", default=str(OUT_DIR))
    parser.add_argument("--report_path", default=str(REPORT_PATH))
    parser.add_argument("--v13_rfc_path", default=str(V13_RFC_PATH))
    parser.add_argument("--include_random_labels", action="store_true")
    parser.add_argument("--prediction_roots", nargs="*", default=[str(path) for path in DEFAULT_ROOTS])
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    roots = [safe_path(item) or Path(item) for item in args.prediction_roots]
    prediction_files = find_prediction_files(roots, args.include_random_labels)
    annotated_frames = []
    skipped = []
    for path in prediction_files:
        result = load_json(path.parent / "result.json")
        try:
            annotated_frames.append(annotate_predictions(path, result))
        except Exception as exc:
            skipped.append({"source_prediction": str(path.relative_to(ROOT)), "error": str(exc)[:500]})
    annotations = pd.concat(annotated_frames, ignore_index=True, sort=False) if annotated_frames else pd.DataFrame()
    annotations.to_csv(out_dir / "prediction_support_annotations.csv", index=False)
    pd.DataFrame(skipped).to_csv(out_dir / "skipped_predictions.csv", index=False)

    summary = build_summary(annotations) if not annotations.empty else pd.DataFrame()
    summary.to_csv(out_dir / "redefined_benchmark_summary.csv", index=False)

    support_metrics = regression_metrics(annotations[annotations["support_covered"].astype(bool)]) if not annotations.empty else regression_metrics(pd.DataFrame())
    unsupported_metrics = regression_metrics(annotations[~annotations["support_covered"].astype(bool)]) if not annotations.empty else regression_metrics(pd.DataFrame())
    old104_metrics = regression_metrics(annotations[annotations["stress_test_group"].eq("gse121141_old104_target_tissue")]) if not annotations.empty else regression_metrics(pd.DataFrame())
    cr_payload = build_cr_metrics(annotations) if not annotations.empty else {"status": "no_predictions", "sources": []}
    write_json(out_dir / "support_covered_headline_metrics.json", support_metrics)
    write_json(out_dir / "unsupported_stress_test_metrics.json", unsupported_metrics)
    write_json(out_dir / "gse121141_old104_stress_metrics.json", old104_metrics)
    write_json(out_dir / "gse80672_cr_redefined_metrics.json", cr_payload)
    write_report(Path(args.report_path), summary, cr_payload, len(prediction_files))
    write_v13_rfc(Path(args.v13_rfc_path))
    state = {
        "status": "completed",
        "prediction_files_scanned": len(prediction_files),
        "prediction_files_skipped": len(skipped),
        "n_annotated_rows": int(len(annotations)),
        "support_covered_headline_metrics": support_metrics,
        "unsupported_stress_test_metrics": unsupported_metrics,
        "gse121141_old104_stress_metrics": old104_metrics,
        "cr_status": cr_payload.get("status"),
        "outputs": {
            "out_dir": str(out_dir),
            "report_path": str(Path(args.report_path)),
            "v13_rfc_path": str(Path(args.v13_rfc_path)),
        },
    }
    write_json(out_dir / "v12_1_redefined_benchmark_state.json", state)
    print(json.dumps(state, indent=2, default=as_jsonable))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Run v8.4 old-tissue candidate and calibration diagnostics.

This script is intentionally diagnostic-only. It does not download large GEO
supplements, does not build matrices, and does not run autoresearch.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import re
import sys
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")
DEFAULT_OUT_DIR = ROOT / "results" / "v8_4_old_tissue_diagnostics"
DEFAULT_REPORT = ROOT / "doc" / "20_analysis" / "19_20260518_v8_4_old_tissue_candidate_calibration_report.md"

CORE_TARGET_TISSUES = {"brain_cortex", "heart", "lung"}
REFERENCE_TISSUES = {"liver"}
TARGET_OLD_WEEKS = 104.0
MONTH_TO_WEEK = 30.42 / 7.0
INTEGRATED_DATASETS = {"GSE120137", "GSE80672", "GSE93957", "GSE121141", "GSE60012", "GSE213628"}

# Official GEO hits from targeted v8.4 web screening. These are not accepted
# blindly; the script runs the same SOFT/filelist preflight used by v7.9.
WEB_SCREENED_ACCESSIONS = [
    "GSE232547",  # aged brain areas + lung, IMPLICON-seq; likely targeted assay.
    "GSE171236",  # aged hippocampal dentate gyrus RRBS/environmental enrichment.
    "GSE138368",  # related aged hippocampus RRBS series referenced by GSE171236.
    "GSE151541",  # aged lung Treg RRBS; cell-type sorted and <104w.
    "GSE225166",  # single-cell/low-coverage blood/liver age prediction.
]


def read_csv(path: Path, **kwargs) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, **kwargs)


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


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


def canonical_tissues(raw: str | float | None) -> set[str]:
    if raw is None or (isinstance(raw, float) and np.isnan(raw)):
        return set()
    text = str(raw).lower().replace("-", "_")
    labels: set[str] = set()
    if "brain_cortex" in text or "cortex" in text or "cortical" in text:
        labels.add("brain_cortex")
    if "hippocampus" in text or "dentate" in text or "cerebell" in text or "prefrontal" in text:
        labels.add("brain_other")
    if "brain" in text and "brain_cortex" not in labels and "cortex" not in text:
        labels.add("brain_other")
    if "heart" in text or "cardiac" in text:
        labels.add("heart")
    if "lung" in text:
        labels.add("lung")
    if "liver" in text or "hepatic" in text:
        labels.add("liver")
    if "blood" in text or "pbmc" in text:
        labels.add("blood")
    if "muscle" in text or "tibialis" in text or "quadriceps" in text:
        labels.add("skeletal_muscle")
    if "intestin" in text or "colon" in text or "organoid" in text:
        labels.add("intestine")
    return labels


def conservative_sample_age_weeks(row: pd.Series) -> float:
    try:
        age = float(row.get("age_weeks"))
        if np.isfinite(age):
            return age
    except (TypeError, ValueError):
        pass
    blob = f"{row.get('title', '')} {row.get('source_name', '')} {row.get('characteristics', '')}".lower()
    patterns = [
        (r"age\s*weeks?\s*:\s*(\d+(?:\.\d+)?)", 1.0),
        (r"(\d+(?:\.\d+)?)\s*[- ]?(?:week|weeks|wk|w)\s*[- ]?old", 1.0),
        (r"(\d+(?:\.\d+)?)\s*[- ]?(?:month|months|mo)\s*[- ]?old", MONTH_TO_WEEK),
        (r"age\s*:\s*(\d+(?:\.\d+)?)\s*[- ]?(?:month|months|mo)\b", MONTH_TO_WEEK),
    ]
    candidates: list[float] = []
    for pattern, multiplier in patterns:
        for match in re.finditer(pattern, blob):
            value = float(match.group(1)) * multiplier
            if 0 <= value <= 220:
                candidates.append(value)
    return float(max(candidates)) if candidates else np.nan


def tissue_counts_for_dataset(samples: pd.DataFrame, dataset: str) -> dict:
    subset = samples[samples["dataset"].astype(str).eq(dataset)].copy()
    if subset.empty:
        return {
            "target_old_sample_count": 0,
            "target_old_tissues": "",
            "reference_old_tissues": "",
            "target_old_counts": "",
            "max_age_by_target_tissue": "",
        }

    rows = []
    for _, row in subset.iterrows():
        age = conservative_sample_age_weeks(row)
        labels = canonical_tissues(row.get("tissue_guess"))
        if not labels:
            labels = canonical_tissues(f"{row.get('source_name', '')} {row.get('title', '')} {row.get('characteristics', '')}")
        for label in labels:
            rows.append({"dataset": dataset, "sample_id": row.get("sample_id"), "age_weeks": age, "tissue": label})
    parsed = pd.DataFrame(rows)
    if parsed.empty:
        return {
            "target_old_sample_count": 0,
            "target_old_tissues": "",
            "reference_old_tissues": "",
            "target_old_counts": "",
            "max_age_by_target_tissue": "",
        }
    old = parsed[parsed["age_weeks"].ge(TARGET_OLD_WEEKS)]
    target_old = old[old["tissue"].isin(CORE_TARGET_TISSUES)]
    ref_old = old[old["tissue"].isin(REFERENCE_TISSUES)]

    def count_string(frame: pd.DataFrame) -> str:
        if frame.empty:
            return ""
        counts = frame.groupby("tissue")["sample_id"].nunique().sort_values(ascending=False)
        return ";".join(f"{idx}:{int(value)}" for idx, value in counts.items())

    def max_age_string(frame: pd.DataFrame) -> str:
        if frame.empty:
            return ""
        ages = frame.groupby("tissue")["age_weeks"].max().sort_index()
        return ";".join(f"{idx}:{value:.1f}" for idx, value in ages.items())

    return {
        "target_old_sample_count": int(target_old["sample_id"].nunique()) if not target_old.empty else 0,
        "target_old_tissues": ";".join(sorted(set(target_old["tissue"]))),
        "reference_old_tissues": ";".join(sorted(set(ref_old["tissue"]))),
        "target_old_counts": count_string(target_old),
        "max_age_by_target_tissue": max_age_string(old[old["tissue"].isin(CORE_TARGET_TISSUES | REFERENCE_TISSUES)]),
    }


def import_discovery_module():
    path = ROOT / "scripts" / "etl" / "16_geo_old_age_candidate_discovery.py"
    spec = importlib.util.spec_from_file_location("geo_old_age_candidate_discovery", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def run_official_screen(out_dir: Path, skip_network: bool) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    inventory_path = out_dir / "official_web_screened_inventory.csv"
    samples_path = out_dir / "official_web_screened_samples.csv"
    supplements_path = out_dir / "official_web_screened_supplements.csv"
    if skip_network:
        return read_csv(inventory_path), read_csv(samples_path), read_csv(supplements_path)
    if inventory_path.exists() and samples_path.exists() and supplements_path.exists():
        return read_csv(inventory_path), read_csv(samples_path), read_csv(supplements_path)

    module = import_discovery_module()
    candidate_rows: list[dict] = []
    sample_rows: list[dict] = []
    supplement_rows: list[dict] = []
    for accession in WEB_SCREENED_ACCESSIONS:
        try:
            row, supplements, samples = module.build_candidate(
                accession,
                source_terms=["v8_4_targeted_official_web_screen"],
                eutils_summary=None,
                refresh_soft=False,
                do_head=True,
            )
        except Exception as exc:  # pragma: no cover - defensive preflight path.
            row = {
                "dataset": accession,
                "geo_url": module.geo_url(accession),
                "priority_tier": "ERROR",
                "recommended_action": "inspect_manually",
                "blockers": str(exc)[:500],
                "source_queries": "v8_4_targeted_official_web_screen",
            }
            supplements = []
            samples = []
        candidate_rows.append(row)
        supplement_rows.extend(supplements)
        sample_rows.extend(samples)

    inventory = pd.DataFrame(candidate_rows)
    samples = pd.DataFrame(sample_rows)
    supplements = pd.DataFrame(supplement_rows)
    out_dir.mkdir(parents=True, exist_ok=True)
    inventory.to_csv(inventory_path, index=False)
    samples.to_csv(samples_path, index=False)
    supplements.to_csv(supplements_path, index=False)
    return inventory, samples, supplements


def score_candidate(row: pd.Series) -> tuple[int, str, str, str]:
    blockers: list[str] = []
    action = "do_not_download"
    role = "not_gse121141_like"
    score = 0
    dataset = str(row.get("dataset", ""))

    if dataset in INTEGRATED_DATASETS:
        blockers.append("already_integrated")
        action = "already_in_matrix_or_current_target"

    schema = str(row.get("processed_schema_guess", "")).lower()
    if "cov" in schema or "bismark" in schema:
        score += 2
    else:
        blockers.append("processed_cov_not_detected")

    target_count = int(row.get("target_old_sample_count", 0) or 0)
    target_tissues = set(str(row.get("target_old_tissues", "")).split(";")) - {""}
    if target_count > 0:
        score += 4
    else:
        blockers.append("no_104w_plus_brain_cortex_heart_lung_samples")
    if len(target_tissues) >= 2:
        score += 4
    elif len(target_tissues) == 1:
        score += 2
        blockers.append("single_target_tissue_only")
    else:
        blockers.append("no_target_tissue_overlap")

    if (row.get("max_age_weeks") or 0) >= TARGET_OLD_WEEKS:
        score += 2
    else:
        blockers.append("max_exact_age_below_104w_or_missing")
    if int(row.get("exact_age_sample_count", 0) or 0) >= 10:
        score += 1
    else:
        blockers.append("limited_sample_specific_age_parse")
    if int(row.get("n_samples_soft", 0) or 0) >= 30:
        score += 1

    text = f"{row.get('title','')} {row.get('summary','')} {row.get('tissues','')}".lower()
    if any(token in text for token in ["impl", "imprinting", "single-cell", "single cell", "treg", "organoid", "hsc"]):
        score -= 3
        blockers.append("targeted_celltype_or_non_bulk_context")
    if "superseries" in text:
        score -= 3
        blockers.append("superseries_use_subseries")

    hard_context_blocker = (
        "targeted_celltype_or_non_bulk_context" in blockers
        or "superseries_use_subseries" in blockers
    )

    if dataset not in INTEGRATED_DATASETS and score >= 9 and len(target_tissues) >= 2 and not hard_context_blocker:
        action = "parser_smoke_then_download_decision"
        role = "p1_candidate_if_schema_smoke_passes"
    elif dataset not in INTEGRATED_DATASETS and target_count > 0:
        action = "adapter_audit_only"
        role = "p2_single_tissue_or_context_limited"
    elif str(row.get("reference_old_tissues", "")):
        action = "auxiliary_validation_only"
        role = "p3_reference_tissue_not_target_gap"

    return score, role, action, ";".join(dict.fromkeys(blockers))


def build_candidate_ranking(existing_inventory: pd.DataFrame, existing_samples: pd.DataFrame, extra_inventory: pd.DataFrame, extra_samples: pd.DataFrame) -> pd.DataFrame:
    inventory = pd.concat([existing_inventory, extra_inventory], ignore_index=True)
    samples = pd.concat([existing_samples, extra_samples], ignore_index=True)
    if inventory.empty:
        return pd.DataFrame()
    inventory = inventory.drop_duplicates(subset=["dataset"], keep="first").copy()
    rows = []
    for _, row in inventory.iterrows():
        dataset = str(row.get("dataset", ""))
        tissue_info = tissue_counts_for_dataset(samples, dataset)
        enriched = row.to_dict()
        enriched.update(tissue_info)
        score, role, action, blockers = score_candidate(pd.Series(enriched))
        enriched.update(
            {
                "v8_4_score": score,
                "v8_4_role": role,
                "v8_4_recommended_action": action,
                "v8_4_blockers": blockers,
                "already_integrated": dataset in INTEGRATED_DATASETS,
            }
        )
        rows.append(enriched)
    ranking = pd.DataFrame(rows)
    preferred_cols = [
        "dataset",
        "v8_4_score",
        "v8_4_role",
        "v8_4_recommended_action",
        "already_integrated",
        "target_old_sample_count",
        "target_old_tissues",
        "target_old_counts",
        "reference_old_tissues",
        "max_age_by_target_tissue",
        "max_age_weeks",
        "n_samples_soft",
        "tissue_count",
        "tissues",
        "processed_schema_guess",
        "preferred_processed_file",
        "preferred_size_bytes",
        "priority_tier",
        "v8_4_blockers",
        "geo_url",
        "title",
    ]
    cols = [col for col in preferred_cols if col in ranking.columns]
    ranking = ranking[cols + [col for col in ranking.columns if col not in cols]]
    return ranking.sort_values(["v8_4_score", "target_old_sample_count", "dataset"], ascending=[False, False, True])


def metric_block(df: pd.DataFrame, pred_col: str) -> dict:
    if df.empty:
        return {"n": 0, "mae_weeks": np.nan, "median_ae_weeks": np.nan, "bias_weeks": np.nan}
    err = df["age_weeks_true"].astype(float) - df[pred_col].astype(float)
    return {
        "n": int(len(df)),
        "mae_weeks": round(float(np.mean(np.abs(err))), 3),
        "median_ae_weeks": round(float(np.median(np.abs(err))), 3),
        "bias_weeks": round(float(np.mean(err)), 3),
    }


def prediction_paths(v8_3_root: Path) -> list[Path]:
    return sorted(v8_3_root.glob("*/*/heldout_gse121141/predictions.csv"))


def load_variant_summary(v8_3_root: Path) -> pd.DataFrame:
    path = v8_3_root / "variant_summary.tsv"
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def is_selected_config(variant_summary: pd.DataFrame, variant: str, config: str) -> bool:
    if variant_summary.empty:
        return False
    subset = variant_summary[
        variant_summary["variant"].astype(str).eq(variant)
        & variant_summary["config"].astype(str).eq(config)
    ]
    if subset.empty or "selected_best_config" not in subset:
        return False
    return bool(subset["selected_best_config"].fillna(False).iloc[0])


def calibrate_predictions(pred: pd.DataFrame) -> pd.DataFrame:
    out = pred.copy()
    out["raw_age_weeks_pred"] = out["age_weeks_pred"].astype(float)
    train = out[out["age_weeks_true"].astype(float) < TARGET_OLD_WEEKS].copy()
    global_offset = float(train["residual_weeks"].median()) if not train.empty else 0.0
    out["global_under104_intercept_offset"] = global_offset
    out["global_under104_intercept_pred"] = out["raw_age_weeks_pred"] + global_offset

    tissue_offsets = train.groupby("tissue")["residual_weeks"].median().to_dict() if not train.empty else {}
    out["tissue_under104_intercept_offset"] = out["tissue"].map(tissue_offsets).fillna(global_offset).astype(float)
    out["tissue_under104_intercept_pred"] = out["raw_age_weeks_pred"] + out["tissue_under104_intercept_offset"]
    out["diagnostic_only"] = True
    out["calibration_fit_note"] = "Offsets fitted on GSE121141 under-104w labels; diagnostic only, not a benchmark."
    return out


def summarize_calibration(pred: pd.DataFrame, variant: str, config: str, selected: bool) -> list[dict]:
    calibrated = calibrate_predictions(pred)
    modes = [
        ("raw", "raw_age_weeks_pred"),
        ("global_under104_intercept", "global_under104_intercept_pred"),
        ("tissue_under104_intercept", "tissue_under104_intercept_pred"),
    ]
    rows = []
    for mode, pred_col in modes:
        for subset_name, subset in [
            ("all", calibrated),
            ("under104", calibrated[calibrated["age_weeks_true"].astype(float) < TARGET_OLD_WEEKS]),
            ("old104plus", calibrated[calibrated["age_weeks_true"].astype(float) >= TARGET_OLD_WEEKS]),
        ]:
            metrics = metric_block(subset, pred_col)
            rows.append(
                {
                    "variant": variant,
                    "config": config,
                    "selected_best_config": selected,
                    "calibration_mode": mode,
                    "subset": subset_name,
                    **metrics,
                }
            )
        old = calibrated[calibrated["age_weeks_true"].astype(float) >= TARGET_OLD_WEEKS]
        for tissue, tissue_df in old.groupby("tissue", dropna=False):
            metrics = metric_block(tissue_df, pred_col)
            rows.append(
                {
                    "variant": variant,
                    "config": config,
                    "selected_best_config": selected,
                    "calibration_mode": mode,
                    "subset": f"old104plus_tissue:{tissue}",
                    **metrics,
                }
            )
    return rows


def build_calibration_outputs(v8_3_root: Path, out_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    variant_summary = load_variant_summary(v8_3_root)
    summary_rows: list[dict] = []
    calibrated_frames: list[pd.DataFrame] = []
    for path in prediction_paths(v8_3_root):
        config = path.parents[1].name
        variant = path.parents[2].name
        pred = pd.read_csv(path)
        selected = is_selected_config(variant_summary, variant, config)
        summary_rows.extend(summarize_calibration(pred, variant, config, selected))
        calibrated = calibrate_predictions(pred)
        calibrated.insert(0, "config", config)
        calibrated.insert(0, "variant", variant)
        calibrated.insert(2, "selected_best_config", selected)
        calibrated_frames.append(calibrated)
    summary = pd.DataFrame(summary_rows)
    calibrated_all = pd.concat(calibrated_frames, ignore_index=True) if calibrated_frames else pd.DataFrame()
    out_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(out_dir / "gse121141_tissue_calibration_summary.csv", index=False)
    calibrated_all.to_csv(out_dir / "gse121141_tissue_calibrated_predictions.csv", index=False)
    return summary, calibrated_all


def build_shared_tissue_support(v8_3_root: Path, out_dir: Path) -> pd.DataFrame:
    support = read_csv(v8_3_root / "gse121141_tissue_age_support.csv")
    if support.empty:
        return pd.DataFrame()
    support["range_status"] = np.where(support["outside_same_tissue_train_age_range"].astype(bool), "outside_same_tissue_age_range", "inside_same_tissue_age_range")
    support["age_subset"] = np.where(support["is_104w_plus"].astype(bool), "old104plus", "under104")
    rows = []
    for keys, group in support.groupby(["variant", "config", "age_subset", "range_status"], dropna=False):
        variant, config, age_subset, range_status = keys
        metrics = metric_block(group.rename(columns={"age_weeks_pred": "pred"}), "pred")
        rows.append(
            {
                "variant": variant,
                "config": config,
                "age_subset": age_subset,
                "range_status": range_status,
                **metrics,
                "median_train_same_tissue_age_max": round(float(group["train_same_tissue_age_max"].median()), 3)
                if group["train_same_tissue_age_max"].notna().any()
                else np.nan,
            }
        )
    summary = pd.DataFrame(rows)

    old_tissue = (
        support[support["is_104w_plus"].astype(bool)]
        .groupby(["variant", "config", "tissue"], dropna=False)
        .agg(
            n=("sample_id", "count"),
            mae_weeks=("abs_error_weeks", "mean"),
            median_abs_error_weeks=("abs_error_weeks", "median"),
            train_same_tissue_n=("train_same_tissue_n", "median"),
            train_same_tissue_age_max=("train_same_tissue_age_max", "max"),
            outside_range_rate=("outside_same_tissue_train_age_range", "mean"),
        )
        .reset_index()
    )
    summary.to_csv(out_dir / "gse121141_shared_tissue_support_summary.csv", index=False)
    old_tissue.to_csv(out_dir / "gse121141_old104_tissue_support_summary.csv", index=False)
    return summary


def write_report(
    path: Path,
    candidate_ranking: pd.DataFrame,
    calibration_summary: pd.DataFrame,
    support_summary: pd.DataFrame,
    official_screen_inventory: pd.DataFrame,
) -> None:
    p1 = candidate_ranking[
        (candidate_ranking["v8_4_recommended_action"] == "parser_smoke_then_download_decision")
        & (~candidate_ranking["already_integrated"].astype(bool))
    ]
    target_view_cols = [
        "dataset",
        "v8_4_score",
        "v8_4_role",
        "v8_4_recommended_action",
        "target_old_counts",
        "reference_old_tissues",
        "max_age_by_target_tissue",
        "processed_schema_guess",
        "v8_4_blockers",
        "title",
    ]
    target_view = candidate_ranking[[col for col in target_view_cols if col in candidate_ranking.columns]].head(14)
    web_view = candidate_ranking[candidate_ranking["dataset"].isin(WEB_SCREENED_ACCESSIONS)].copy()
    web_view = web_view[
        [
            col
            for col in [
                "dataset",
                "v8_4_role",
                "v8_4_recommended_action",
                "target_old_counts",
                "processed_schema_guess",
                "max_age_weeks",
                "v8_4_blockers",
                "geo_url",
            ]
            if col in web_view.columns
        ]
    ]

    selected_cal = calibration_summary[
        (calibration_summary["selected_best_config"].astype(bool))
        & (calibration_summary["subset"].isin(["all", "old104plus"]))
    ].copy()
    selected_cal = selected_cal[
        [
            "variant",
            "config",
            "calibration_mode",
            "subset",
            "n",
            "mae_weeks",
            "median_ae_weeks",
            "bias_weeks",
        ]
    ]

    old_support = support_summary[support_summary["age_subset"].eq("old104plus")].copy()
    old_support = old_support[
        [
            "variant",
            "config",
            "range_status",
            "n",
            "mae_weeks",
            "median_ae_weeks",
            "median_train_same_tissue_age_max",
        ]
    ]

    lines = [
        "# v8.4 Old Tissue Candidate and Calibration Diagnostics Report",
        "",
        "Date: 2026-05-18",
        "",
        "## Summary",
        "",
        "v8.4 did not run autoresearch, deep learning, FASTQ/Bismark, or large GEO downloads. It screened official GEO/SOFT/filelist metadata for GSE121141-like old-age tissues and ran diagnostic-only calibration on existing v8.3 GSE121141 held-out predictions.",
        "",
        f"- Non-integrated P1 candidates found: {len(p1)}",
        f"- Official web-screened accessions checked: {', '.join(WEB_SCREENED_ACCESSIONS)}",
        "- Core target tissues: brain_cortex, heart, lung; liver is reference-only.",
        "",
        "## Candidate Ranking",
        "",
        md_table(target_view),
        "",
        "## Official Web Screen Additions",
        "",
        md_table(web_view, max_rows=10),
        "",
        "## Shared-Tissue Age-Range Diagnostic",
        "",
        "For GSE121141 104w+ samples, this table separates samples that are inside vs outside the same-tissue age range available in training. If old samples are outside range, model tuning is unlikely to fix the extrapolation problem by itself.",
        "",
        md_table(old_support),
        "",
        "## Tissue-Specific Calibration Diagnostic",
        "",
        "Calibration modes using GSE121141 labels are diagnostic-only. They test whether errors look like a tissue offset, not whether a deployable cross-dataset model has improved.",
        "",
        md_table(selected_cal),
        "",
        "## Decision",
        "",
    ]
    if p1.empty:
        lines.extend(
            [
                "No new non-integrated P1 dataset was identified for old brain_cortex/heart/lung bulk RRBS with a clean processed COV schema. Do not start v8.4 downloads or autoresearch from this screen.",
                "",
                "The next useful branch is a narrower official GEO refresh for cortex/heart/lung old-age methylation datasets, plus parser smoke only for candidates that have sample-specific age, bulk tissue, and processed COV/bedGraph-like methylation files. Current P2/P3 hits should remain auxiliary or method-specific validation data.",
            ]
        )
    else:
        lines.extend(
            [
                "At least one non-integrated P1 candidate was found. The next step is parser smoke only, followed by matrix construction if sample IDs, coordinates, and beta ranges validate.",
            ]
        )
    lines.extend(
        [
            "",
            "## Files",
            "",
            "- `results/v8_4_old_tissue_diagnostics/old_tissue_candidate_ranking.csv`",
            "- `results/v8_4_old_tissue_diagnostics/official_web_screened_inventory.csv`",
            "- `results/v8_4_old_tissue_diagnostics/gse121141_shared_tissue_support_summary.csv`",
            "- `results/v8_4_old_tissue_diagnostics/gse121141_tissue_calibration_summary.csv`",
            "- `results/v8_4_old_tissue_diagnostics/gse121141_tissue_calibrated_predictions.csv`",
            "",
            "## Official Source URLs",
            "",
            "- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121141",
            "- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213628",
            "- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE232547",
            "- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171236",
            "- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151541",
            "- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225166",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out_dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument("--v8_3_root", default=str(ROOT / "results" / "benchmark_v8_3_ablation"))
    parser.add_argument("--candidate_inventory", default=str(ROOT / "metadata" / "geo_old_age_candidate_inventory.csv"))
    parser.add_argument("--candidate_samples", default=str(ROOT / "metadata" / "geo_old_age_candidate_samples.csv"))
    parser.add_argument("--skip_network", action="store_true", help="Use cached v8.4 official screen files only.")
    parser.add_argument("--report", default=str(DEFAULT_REPORT))
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    v8_3_root = Path(args.v8_3_root)

    existing_inventory = read_csv(Path(args.candidate_inventory))
    existing_samples = read_csv(Path(args.candidate_samples))
    official_inventory, official_samples, official_supplements = run_official_screen(out_dir, skip_network=args.skip_network)

    candidate_ranking = build_candidate_ranking(existing_inventory, existing_samples, official_inventory, official_samples)
    candidate_ranking.to_csv(out_dir / "old_tissue_candidate_ranking.csv", index=False)
    official_supplements.to_csv(out_dir / "official_web_screened_supplements.csv", index=False)

    calibration_summary, calibrated_predictions = build_calibration_outputs(v8_3_root, out_dir)
    support_summary = build_shared_tissue_support(v8_3_root, out_dir)

    payload = {
        "status": "complete",
        "non_integrated_p1_candidates": int(
            (
                (candidate_ranking.get("v8_4_recommended_action") == "parser_smoke_then_download_decision")
                & (~candidate_ranking.get("already_integrated", pd.Series(False, index=candidate_ranking.index)).astype(bool))
            ).sum()
        )
        if not candidate_ranking.empty
        else 0,
        "web_screened_accessions": WEB_SCREENED_ACCESSIONS,
        "outputs": {
            "candidate_ranking": str(out_dir / "old_tissue_candidate_ranking.csv"),
            "calibration_summary": str(out_dir / "gse121141_tissue_calibration_summary.csv"),
            "calibrated_predictions": str(out_dir / "gse121141_tissue_calibrated_predictions.csv"),
            "shared_tissue_support": str(out_dir / "gse121141_shared_tissue_support_summary.csv"),
            "report": str(Path(args.report)),
        },
    }
    write_json(out_dir / "v8_4_summary.json", payload)
    write_report(Path(args.report), candidate_ranking, calibration_summary, support_summary, official_inventory)
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

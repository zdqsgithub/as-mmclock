#!/usr/bin/env python3
"""v9 RALPH controller for mouse methylation clock data strategy.

This controller turns the v8 bottleneck into a repeatable gate process.  It is
safe by default: without --refresh_candidates or --run_smoke it only reads cached
preflight/smoke outputs and writes v9 decision ledgers.  Network metadata
refresh and small-file adapter smoke are explicit opt-ins.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import re
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")
OUT_DIR = ROOT / "results" / "ralph_v9_loop"
REPORT = OUT_DIR / "v9_success_or_failure_report.md"
REFERENCE_MATRIX = ROOT / "results" / "multidataset_v8_2_prefilter_liftover" / "all_rrbs_region_matrix_5kb.parquet"

V85_DIR = ROOT / "results" / "v8_5_targeted_geo_refresh"
V86_DIR = ROOT / "results" / "v8_6_calibration_poc"
V87_DIR = ROOT / "results" / "ralph_v8_loop"

TARGET_TISSUES = {"brain_cortex", "heart", "lung"}
OLD_THRESHOLD_WEEKS = 104.0
COMMON_REGION_GATE = 50_000
MAX_SMOKE_FILES_PER_DATASET = 3
MAX_SMOKE_FILE_BYTES = 60_000_000

V75_BASELINES = {
    "groupkfold_mae_weeks": 23.767,
    "gse121141_all_age_mae_weeks": 38.033,
    "gse121141_old104_mae_weeks": 75.386,
}

OFFICIAL_SOURCE_DOCS = {
    "geo_download": "https://www.ncbi.nlm.nih.gov/geo/info/download.html",
    "geo_programmatic_access": "https://www.ncbi.nlm.nih.gov/geo/info/geo_paccess.html",
    "geo_soft": "https://www.ncbi.nlm.nih.gov/geo/info/soft.html",
    "sra_download": "https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/",
    "ena_browser_api": "https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/browser-api.html",
}

V9_SEARCH_TERMS = [
    'Mus musculus[ORGN] AND "Methylation profiling by high throughput sequencing" AND (RRBS OR WGBS OR bisulfite) AND ("brain cortex" OR cortex OR heart OR lung) AND (aging OR aged OR old OR "24 month" OR "26 month" OR "30 month") AND gse[ETYP]',
    'Mus musculus[ORGN] AND "DNA methylation" AND ("brain cortex" OR cortex OR heart OR lung) AND ("24 month" OR "26 month" OR "30 month" OR aged OR old) AND gse[ETYP]',
    'Mus musculus[ORGN] AND "reduced representation bisulfite" AND (heart OR lung OR cortex) AND aging AND gse[ETYP]',
]

V9_SEED_ACCESSIONS = [
    "GSE225166",
    "GSE83947",
    "GSE134398",
    "GSE232547",
    "GSE171236",
    "GSE138368",
    "GSE151541",
    "GSE169234",
    "GSE103249",
    "GSE108762",
    "GSE281602",
    "GSE156557",
    "GSE215310",
    "GSE221124",
]

INTEGRATED_DATASETS = {"GSE120137", "GSE80672", "GSE93957", "GSE121141", "GSE60012", "GSE213628"}
HARD_BLOCKER_PATTERNS = {
    "single_cell_or_targeted_celltype": ["single-cell", "single cell", "targeted_celltype_or_non_bulk_context"],
    "low_coverage_or_tagged": ["ipcrtag", "itag", "low coverage", "low_coverage_or_tagged_context"],
    "untraceable_superseries": ["superseries_use_subseries", "superseries_or_subseries_required"],
    "organoid_or_in_vitro": ["organoid", "in vitro modeling"],
    "celltype_specific_hsc": [" hsc", "hematopoietic stem"],
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def import_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def discovery_module():
    return import_module(ROOT / "scripts" / "etl" / "16_geo_old_age_candidate_discovery.py", "geo_old_age_candidate_discovery_v9")


def v84_module():
    return import_module(ROOT / "scripts" / "validate" / "v8_4_old_tissue_candidate_and_calibration.py", "v8_4_helpers_v9")


def v8_loop_module():
    return import_module(ROOT / "scripts" / "validate" / "run_v8_ralph_loop.py", "run_v8_ralph_loop_helpers_v9")


def read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path) if path.exists() else pd.DataFrame()


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def append_jsonl(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(payload, sort_keys=True) + "\n")


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


def numeric(value: object) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return float("nan")
    if np.isfinite(parsed):
        return parsed
    return float("nan")


def combine_frames(paths: list[Path], source_label: str) -> pd.DataFrame:
    frames = []
    for path in paths:
        df = read_csv(path)
        if not df.empty:
            df = df.copy()
            df["v9_source_table"] = source_label + ":" + str(path.relative_to(ROOT))
            frames.append(df)
    return pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()


def load_cached_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    inventory = combine_frames(
        [
            ROOT / "metadata" / "geo_old_age_candidate_inventory.csv",
            V85_DIR / "targeted_search_inventory.csv",
            OUT_DIR / "v9_research_inventory.csv",
        ],
        "cached_inventory",
    )
    samples = combine_frames(
        [
            ROOT / "metadata" / "geo_old_age_candidate_samples.csv",
            V85_DIR / "targeted_search_samples.csv",
            OUT_DIR / "v9_research_samples.csv",
        ],
        "cached_samples",
    )
    supplements = combine_frames(
        [
            ROOT / "metadata" / "geo_old_age_candidate_supplements.csv",
            V85_DIR / "targeted_search_supplements.csv",
            OUT_DIR / "v9_research_supplements.csv",
        ],
        "cached_supplements",
    )
    if not inventory.empty and "dataset" in inventory:
        inventory = inventory.drop_duplicates("dataset", keep="last").reset_index(drop=True)
    if not samples.empty and {"dataset", "sample_id"}.issubset(samples.columns):
        samples = samples.drop_duplicates(["dataset", "sample_id"], keep="last").reset_index(drop=True)
    if not supplements.empty and {"dataset", "supplement_name"}.issubset(supplements.columns):
        supplements = supplements.drop_duplicates(["dataset", "supplement_name"], keep="last").reset_index(drop=True)
    return inventory, samples, supplements


def refresh_candidates(module, out_dir: Path, retmax: int, max_candidates: int) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    accessions = set(V9_SEED_ACCESSIONS)
    source_terms: dict[str, list[str]] = {accession: ["manual_seed_v9"] for accession in accessions}
    summaries: dict[str, dict] = {}
    found, found_terms, found_summaries = module.discover_accessions(V9_SEARCH_TERMS, retmax=retmax)
    accessions.update(found)
    summaries.update(found_summaries)
    for accession, terms in found_terms.items():
        source_terms.setdefault(accession, []).extend(terms)
    for accession in sorted(accessions):
        if accession not in summaries:
            try:
                summaries[accession] = module.eutils_summary_for_accession(accession)
            except Exception as exc:
                summaries[accession] = {"dataset": accession, "eutils_error": str(exc)[:300]}
    ranked = sorted(
        accessions,
        key=lambda value: (
            -module.rough_esummary_score(value, summaries.get(value, {}), source_terms.get(value, [])),
            int(value.removeprefix("GSE")) if value.removeprefix("GSE").isdigit() else 999_999_999,
        ),
    )
    selected = list(dict.fromkeys(ranked[:max_candidates] + V9_SEED_ACCESSIONS))
    candidate_rows: list[dict] = []
    sample_rows: list[dict] = []
    supplement_rows: list[dict] = []
    for accession in selected:
        try:
            row, supplements, samples = module.build_candidate(
                accession,
                source_terms=source_terms.get(accession, ["v9_search"]),
                eutils_summary=summaries.get(accession),
                refresh_soft=False,
                do_head=True,
            )
        except Exception as exc:
            row = {
                "dataset": accession,
                "geo_url": module.geo_url(accession),
                "priority_tier": "ERROR",
                "recommended_action": "inspect_manually",
                "blockers": str(exc)[:500],
                "source_queries": " || ".join(source_terms.get(accession, [])),
            }
            supplements = []
            samples = []
        row["v9_source_table"] = "v9_refresh"
        candidate_rows.append(row)
        for sample in samples:
            sample["v9_source_table"] = "v9_refresh"
        for supplement in supplements:
            supplement["v9_source_table"] = "v9_refresh"
        sample_rows.extend(samples)
        supplement_rows.extend(supplements)
        time.sleep(0.05)
    inventory = pd.DataFrame(candidate_rows)
    samples = pd.DataFrame(sample_rows)
    supplements = pd.DataFrame(supplement_rows)
    inventory.to_csv(out_dir / "v9_research_inventory.csv", index=False)
    samples.to_csv(out_dir / "v9_research_samples.csv", index=False)
    supplements.to_csv(out_dir / "v9_research_supplements.csv", index=False)
    return inventory, samples, supplements


def normalize_tissue(v84, row: pd.Series) -> str:
    labels = v84.canonical_tissues(row.get("tissue_guess"))
    if not labels:
        labels = v84.canonical_tissues(f"{row.get('source_name', '')} {row.get('title', '')} {row.get('characteristics', '')}")
    if "brain_cortex" in labels:
        return "brain_cortex"
    if "brain_other" in labels:
        return "brain_other"
    if labels:
        return sorted(labels)[0]
    raw = str(row.get("tissue_guess") or row.get("source_name") or "unknown")
    return raw.strip() or "unknown"


def sample_age_weeks(v84, row: pd.Series) -> float:
    age = numeric(row.get("age_weeks"))
    if np.isfinite(age):
        return age
    try:
        age = float(v84.conservative_sample_age_weeks(row))
    except Exception:
        return float("nan")
    return age if np.isfinite(age) else float("nan")


def file_size(row: pd.Series) -> int:
    value = row.get("filelist_size_bytes")
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def match_supplement_rows(samples: pd.DataFrame, supplements: pd.DataFrame, v84) -> pd.DataFrame:
    if samples.empty or supplements.empty:
        return pd.DataFrame()
    file_type = supplements.get("filelist_type", pd.Series(dtype=str)).astype(str).str.upper()
    archive_or_file = supplements.get("archive_or_file", pd.Series(dtype=str)).astype(str).str.lower()
    file_rows = supplements[archive_or_file.eq("file") & file_type.isin(["COV", "TXT", "BEDGRAPH"])].copy()
    rows = []
    for _, sample in samples.iterrows():
        sample_id = str(sample.get("sample_id", ""))
        if not sample_id:
            continue
        matches = file_rows[file_rows["supplement_name"].astype(str).str.contains(sample_id, regex=False)] if not file_rows.empty else pd.DataFrame()
        age = sample_age_weeks(v84, sample)
        tissue = normalize_tissue(v84, sample)
        for _, file_row in matches.iterrows():
            size_value = file_size(file_row)
            rows.append(
                {
                    "dataset": sample.get("dataset"),
                    "sample_id": sample_id,
                    "title": sample.get("title", ""),
                    "source_name": sample.get("source_name", ""),
                    "tissue": tissue,
                    "age_weeks": age,
                    "is_old104": bool(np.isfinite(age) and age >= OLD_THRESHOLD_WEEKS),
                    "is_target_tissue": tissue in TARGET_TISSUES,
                    "supplement_name": file_row.get("supplement_name", ""),
                    "filelist_type": str(file_row.get("filelist_type", "")).upper(),
                    "file_size_bytes": size_value,
                    "supplement_url": file_row.get("supplement_url", ""),
                    "under_size_cap": bool(size_value and size_value <= MAX_SMOKE_FILE_BYTES),
                }
            )
    return pd.DataFrame(rows)


def select_smoke_files(matched: pd.DataFrame) -> pd.DataFrame:
    if matched.empty:
        return matched
    data = matched[matched["under_size_cap"].astype(bool)].copy()
    if data.empty:
        return data
    data["priority"] = 0
    data.loc[data["is_target_tissue"].astype(bool), "priority"] += 10
    data.loc[data["is_old104"].astype(bool), "priority"] += 10
    data.loc[data["tissue"].eq("brain_cortex"), "priority"] += 5
    data.loc[data["filelist_type"].eq("COV"), "priority"] += 3
    data.loc[data["filelist_type"].eq("TXT"), "priority"] += 1
    data["size_rank"] = data["file_size_bytes"].astype(float)
    return data.sort_values(["priority", "size_rank"], ascending=[False, True]).head(MAX_SMOKE_FILES_PER_DATASET)


def context_blockers_for_candidate(dataset: str, inventory_row: pd.Series | None, dataset_supplements: pd.DataFrame | None = None) -> list[str]:
    text = dataset.lower()
    if inventory_row is not None:
        text += " " + " ".join(
            str(inventory_row.get(field, ""))
            for field in ["title", "summary", "overall_design", "blockers", "v8_4_blockers", "processed_schema_guess"]
        ).lower()
    if dataset_supplements is not None and not dataset_supplements.empty:
        names = dataset_supplements.get("supplement_name", pd.Series(dtype=str)).fillna("").astype(str).tolist()
        text += " " + " ".join(names).lower()
    blockers: list[str] = []
    for label, needles in HARD_BLOCKER_PATTERNS.items():
        if any(needle in text for needle in needles):
            blockers.append(label)
    return list(dict.fromkeys(blockers))


def sample_support(samples: pd.DataFrame, v84) -> pd.DataFrame:
    if samples.empty:
        return pd.DataFrame(columns=["sample_id", "tissue", "age_weeks", "is_old104", "is_target_tissue"])
    rows = []
    for _, sample in samples.iterrows():
        age = sample_age_weeks(v84, sample)
        tissue = normalize_tissue(v84, sample)
        rows.append(
            {
                "sample_id": sample.get("sample_id", ""),
                "tissue": tissue,
                "age_weeks": age,
                "is_old104": bool(np.isfinite(age) and age >= OLD_THRESHOLD_WEEKS),
                "is_target_tissue": tissue in TARGET_TISSUES,
            }
        )
    return pd.DataFrame(rows)


def existing_smoke_cache(dataset: str, selected: pd.DataFrame) -> pd.DataFrame:
    path = V87_DIR / f"{dataset}_adapter_smoke.csv"
    if selected.empty or not path.exists():
        return pd.DataFrame()
    smoke = read_csv(path)
    if smoke.empty:
        return smoke
    key_cols = ["sample_id", "supplement_name"]
    if not set(key_cols).issubset(smoke.columns) or not set(key_cols).issubset(selected.columns):
        return pd.DataFrame()
    keys = set(zip(selected["sample_id"].astype(str), selected["supplement_name"].astype(str)))
    smoke_keys = list(zip(smoke["sample_id"].astype(str), smoke["supplement_name"].astype(str)))
    cached = smoke[[item in keys for item in smoke_keys]].copy()
    if cached.empty:
        return cached
    cached["status"] = "existing_smoke_cache"
    return cached


def load_reference_regions(v8) -> set[str]:
    if REFERENCE_MATRIX.exists():
        return v8.load_reference_regions()
    return set()


def run_or_plan_smoke(
    dataset: str,
    selected: pd.DataFrame,
    out_dir: Path,
    v8,
    reference_regions: set[str],
    run_smoke: bool,
    network_log: Path,
) -> pd.DataFrame:
    smoke_dir = out_dir / "smoke_files" / dataset
    smoke_log = out_dir / "smoke_download_log.jsonl"
    rows = []
    if selected.empty:
        return pd.DataFrame()
    if not run_smoke:
        cached = existing_smoke_cache(dataset, selected)
        if not cached.empty:
            return cached
        planned = selected.copy()
        planned["status"] = "planned_dry_run_no_download"
        planned["schema_guess"] = "not_assessed"
        planned["rows_parseable_primary_autosomes"] = np.nan
        planned["rows_pass_coverage_ge5"] = np.nan
        planned["beta_range_valid"] = False
        planned["n_smoke_regions"] = np.nan
        planned["n_smoke_common_regions_with_reference"] = np.nan
        return planned
    for _, selected_row in selected.iterrows():
        filename = Path(str(selected_row["supplement_name"])).name
        url_candidates = v8.smoke_url_candidates(selected_row, filename)
        append_jsonl(
            network_log,
            {
                "timestamp": utc_now(),
                "dataset": dataset,
                "blocker_type": "adapter_smoke_url_resolution",
                "attempted_url": " | ".join(url for _, url in url_candidates),
                "source_doc": OFFICIAL_SOURCE_DOCS["geo_download"],
                "resolution": "try_sample_supplement_then_inventory_url",
                "applied_rule": "individual_sample_supplement_before_series_inventory",
            },
        )
        local_path, download_payload = v8.download_smoke_file(selected_row, smoke_dir, smoke_log)
        parse_payload = {}
        if local_path is not None:
            try:
                parse_payload = v8.parse_smoke_file(local_path, reference_regions)
            except Exception as exc:
                parse_payload = {"local_path": str(local_path), "parse_error": str(exc)[:500]}
        if download_payload.get("status") == "failed":
            append_jsonl(
                network_log,
                {
                    "timestamp": utc_now(),
                    "dataset": dataset,
                    "blocker_type": "download_failed",
                    "attempted_url": json.dumps(download_payload.get("url_errors", []))[:1000],
                    "source_doc": OFFICIAL_SOURCE_DOCS["geo_download"],
                    "resolution": "record_failure_and_keep_candidate_at_audit_gate",
                    "applied_rule": "no_tar_download_from_adapter_smoke_failure",
                },
            )
        rows.append({**selected_row.to_dict(), **download_payload, **parse_payload})
    return pd.DataFrame(rows)


def schema_smoke_pass(smoke: pd.DataFrame) -> bool:
    if smoke.empty:
        return False
    parseable = pd.to_numeric(smoke.get("rows_parseable_primary_autosomes", pd.Series(0, index=smoke.index)), errors="coerce").fillna(0).gt(0)
    coverage = pd.to_numeric(smoke.get("rows_pass_coverage_ge5", pd.Series(0, index=smoke.index)), errors="coerce").fillna(0).gt(0)
    beta_valid = smoke.get("beta_range_valid", pd.Series(False, index=smoke.index)).fillna(False).astype(bool)
    status_ok = smoke.get("status", pd.Series("", index=smoke.index)).astype(str).isin(["completed", "skipped", "existing_smoke_cache"])
    return bool((parseable & coverage & beta_valid & status_ok).any())


def smoke_common_regions(smoke: pd.DataFrame) -> int:
    if smoke.empty or "n_smoke_common_regions_with_reference" not in smoke:
        return 0
    values = pd.to_numeric(smoke["n_smoke_common_regions_with_reference"], errors="coerce").fillna(0)
    return int(values.max()) if not values.empty else 0


def assembly_compatible(inventory_row: pd.Series | None, smoke: pd.DataFrame) -> bool:
    text = ""
    if inventory_row is not None:
        text += " ".join(str(inventory_row.get(field, "")) for field in ["title", "summary", "overall_design", "processed_schema_guess"])
    if not smoke.empty:
        text += " " + " ".join(smoke.get("supplement_name", pd.Series(dtype=str)).astype(str).tolist())
    text = text.lower()
    return any(token in text for token in ["grcm38", "mm10", "mm9", "ncbi37"])


def classify_candidate(
    old_target_tissues: set[str],
    official_metadata_pass: bool,
    has_processed_file: bool,
    hard_blockers: list[str],
    schema_pass: bool,
    assembly_pass: bool,
) -> tuple[str, str, bool]:
    if hard_blockers:
        return "P4_blocked_auxiliary", ";".join(hard_blockers), False
    if official_metadata_pass and has_processed_file and ("brain_cortex" in old_target_tissues or len(old_target_tissues) >= 2):
        if schema_pass and assembly_pass:
            return "P1_headline_candidate", "matrix_gate_required_before_headline", True
        return "P1_headline_candidate", "adapter_smoke_or_assembly_gate_pending", False
    if official_metadata_pass and has_processed_file and len(old_target_tissues) == 1:
        return "P2_single_target_auxiliary", "single_target_tissue_only", False
    if official_metadata_pass and not has_processed_file:
        return "P3_minimal_fastq_pilot", "processed_matrix_absent_consider_2_3_sample_fastq_pilot", False
    return "P4_blocked_auxiliary", "no_sample_specific_old_target_tissue", False


def summarize_candidate(
    dataset: str,
    inventory: pd.DataFrame,
    samples: pd.DataFrame,
    supplements: pd.DataFrame,
    v84,
    v8,
    reference_regions: set[str],
    out_dir: Path,
    run_smoke: bool,
    network_log: Path,
) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    inv_row = None
    if not inventory.empty and "dataset" in inventory and dataset in set(inventory["dataset"].astype(str)):
        inv_row = inventory[inventory["dataset"].astype(str).eq(dataset)].iloc[0]
    ds_samples = samples[samples["dataset"].astype(str).eq(dataset)].copy() if not samples.empty else pd.DataFrame()
    ds_supplements = supplements[supplements["dataset"].astype(str).eq(dataset)].copy() if not supplements.empty else pd.DataFrame()
    support = sample_support(ds_samples, v84)
    old_target = support[support["is_old104"].astype(bool) & support["is_target_tissue"].astype(bool)] if not support.empty else pd.DataFrame()
    old_target_tissues = set(old_target["tissue"].astype(str)) if not old_target.empty else set()
    matched = match_supplement_rows(ds_samples, ds_supplements, v84)
    selected = select_smoke_files(matched)
    smoke = run_or_plan_smoke(dataset, selected, out_dir, v8, reference_regions, run_smoke, network_log)
    matched.to_csv(out_dir / f"{dataset}_matched_supplement_files.csv", index=False)
    smoke.to_csv(out_dir / f"{dataset}_adapter_smoke.csv", index=False)

    blockers = context_blockers_for_candidate(dataset, inv_row, ds_supplements)
    has_processed_file = bool(not matched.empty)
    metadata_pass = bool(not old_target.empty)
    schema_pass = schema_smoke_pass(smoke)
    assembly_pass = assembly_compatible(inv_row, smoke)
    candidate_tier, auxiliary_reason, p1_ready = classify_candidate(
        old_target_tissues,
        metadata_pass,
        has_processed_file,
        blockers,
        schema_pass,
        assembly_pass,
    )
    matrix_gate_pass = False
    headline_allowed = False
    if p1_ready:
        headline_allowed = False
        auxiliary_reason = "full_matrix_common_region_gate_not_assessed"
    next_step = {
        "P1_headline_candidate": "run_adapter_smoke_or_build_full_matrix_gate",
        "P2_single_target_auxiliary": "keep_for_single_tissue_diagnostics_only",
        "P3_minimal_fastq_pilot": "plan_2_3_sample_fastq_bismark_pilot",
        "P4_blocked_auxiliary": "do_not_promote_to_headline",
    }[candidate_tier]
    row = {
        "iteration": "v9.0",
        "dataset": dataset,
        "candidate_tier": candidate_tier,
        "hypothesis": "candidate may repair GSE121141 old same-tissue support",
        "action": "research_audit_smoke_gate_no_matrix_training",
        "official_metadata_pass": metadata_pass,
        "download_allowed": bool(not selected.empty),
        "smoke_download_pass": bool(not smoke.empty and smoke.get("status", pd.Series(dtype=str)).astype(str).isin(["completed", "skipped", "existing_smoke_cache"]).any()),
        "n_samples": int(len(ds_samples)),
        "n_old_target_samples": int(len(old_target)),
        "old_target_tissues": ";".join(sorted(old_target_tissues)),
        "n_matched_processed_files": int(len(matched)),
        "n_selected_smoke_files": int(len(selected)),
        "schema_smoke_pass": schema_pass,
        "assembly_compatible": assembly_pass,
        "smoke_common_regions_estimate": smoke_common_regions(smoke),
        "common_regions_pass": matrix_gate_pass,
        "common_regions_required": COMMON_REGION_GATE,
        "headline_candidate": bool(candidate_tier == "P1_headline_candidate"),
        "headline_allowed": headline_allowed,
        "minimal_fastq_pilot_candidate": bool(candidate_tier == "P3_minimal_fastq_pilot"),
        "auxiliary_only_reason": auxiliary_reason,
        "context_blockers": ";".join(blockers),
        "next_step": next_step,
    }
    return row, matched, smoke


def candidate_score(row: pd.Series) -> float:
    tier = str(row.get("priority_tier", ""))
    score = {
        "P1_download_next": 100.0,
        "P2_adapter_audit": 80.0,
        "P3_auxiliary_or_intervention": 45.0,
        "P4_low_priority": 5.0,
        "ERROR": -10.0,
    }.get(tier, 0.0)
    score += min(numeric(row.get("max_age_months")) if np.isfinite(numeric(row.get("max_age_months"))) else 0.0, 36.0)
    score += min(numeric(row.get("exact_age_sample_count")) if np.isfinite(numeric(row.get("exact_age_sample_count"))) else 0.0, 50.0) / 5.0
    tissues = str(row.get("tissues", "")).lower()
    for tissue in TARGET_TISSUES:
        if tissue.replace("_", " ") in tissues or tissue in tissues:
            score += 6.0
    schema = str(row.get("processed_schema_guess", "")).lower()
    if "cov" in schema:
        score += 8.0
    blockers = str(row.get("blockers", "")).lower()
    if "superseries" in blockers:
        score -= 20.0
    if "single_tissue_or_celltype" in blockers:
        score -= 6.0
    if "celltype" in blockers or "organoid" in blockers:
        score -= 20.0
    return score


def candidate_order(inventory: pd.DataFrame, samples: pd.DataFrame) -> list[str]:
    scored: list[tuple[float, str]] = []
    if not inventory.empty and "dataset" in inventory:
        for _, row in inventory.iterrows():
            dataset = str(row.get("dataset", ""))
            if dataset.startswith("GSE") and dataset not in INTEGRATED_DATASETS:
                scored.append((candidate_score(row), dataset))
    if not samples.empty and "dataset" in samples:
        for dataset in samples["dataset"].dropna().astype(str).tolist():
            if dataset.startswith("GSE") and dataset not in INTEGRATED_DATASETS:
                scored.append((0.0, dataset))
    for dataset in V9_SEED_ACCESSIONS:
        if dataset not in INTEGRATED_DATASETS:
            scored.append((1.0, dataset))
    ordered = []
    for _, dataset in sorted(scored, key=lambda item: (-item[0], item[1])):
        if dataset not in ordered:
            ordered.append(dataset)
    return ordered


def load_v86_summary() -> dict:
    path = V86_DIR / "v8_6_summary.json"
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}


def decide_state(gate: pd.DataFrame, refresh_round: int, dry_run: bool) -> dict:
    headline_ready = gate[gate["headline_candidate"].astype(bool) & gate["schema_smoke_pass"].astype(bool)] if not gate.empty else pd.DataFrame()
    headline_pending = gate[gate["headline_candidate"].astype(bool)] if not gate.empty else pd.DataFrame()
    pilot_pending = gate[gate["minimal_fastq_pilot_candidate"].astype(bool)] if not gate.empty else pd.DataFrame()
    only_auxiliary = bool(
        not gate.empty
        and headline_pending.empty
        and pilot_pending.empty
        and gate["candidate_tier"].astype(str).isin(["P2_single_target_auxiliary", "P4_blocked_auxiliary"]).all()
    )
    failure = bool(refresh_round >= 2 and only_auxiliary)
    if not headline_ready.empty:
        decision = "continue_v9_matrix_gate"
        next_action = "Build candidate matrix, require common regions >=50000, then run fixed held-out benchmarks."
        status = "continue"
    elif not headline_pending.empty or not pilot_pending.empty:
        decision = "continue_v9_audit_or_pilot"
        next_action = "Run adapter smoke for P1 candidates or design minimal FASTQ pilot for P3 candidates."
        status = "continue"
    elif failure:
        decision = "v9_failure_no_headline_data_after_two_refresh_rounds"
        next_action = "Stop v9 loop and move to broader processed methylome acquisition or minimal ETL planning."
        status = "failure"
    else:
        decision = "continue_v9_research"
        next_action = "Run a second official candidate refresh before declaring v9 failure."
        status = "continue"
    return {
        "timestamp": utc_now(),
        "loop_version": "v9.0",
        "status": status,
        "dry_run": dry_run,
        "refresh_round": refresh_round,
        "loop_decision": decision,
        "next_action": next_action,
        "success_criteria": {
            "headline_dataset_or_minimal_etl_dataset": "sample-specific age, old target tissue, parseable methylation, assembly/liftover, common 5kb regions >=50000",
            "gse121141_old104_mae_improvement_weeks": ">=10 vs 75.386w",
            "gse121141_all_age_mae_max": "<=40.033w",
            "groupkfold_mae_max": "<=25.767w",
            "random_label_sanity": "abs(r)<0.2 and MAE random-like",
            "shuffled_cr_sanity": "AUC chance-like and Cohen_d decreased",
        },
        "failure_criteria": {
            "two_refresh_rounds_without_p1_or_p3": bool(refresh_round >= 2 and headline_pending.empty and pilot_pending.empty),
            "only_auxiliary_contexts_available": only_auxiliary,
            "calibration_not_transferable_from_v8_6": True,
            "common_regions_lt_50000_blocks_training": True,
        },
        "metrics": {
            "n_candidates": int(len(gate)),
            "n_p1_headline_candidates": int((gate["candidate_tier"].astype(str).eq("P1_headline_candidate")).sum()) if not gate.empty else 0,
            "n_p2_auxiliary": int((gate["candidate_tier"].astype(str).eq("P2_single_target_auxiliary")).sum()) if not gate.empty else 0,
            "n_p3_fastq_pilot_candidates": int((gate["candidate_tier"].astype(str).eq("P3_minimal_fastq_pilot")).sum()) if not gate.empty else 0,
            "n_p4_blocked": int((gate["candidate_tier"].astype(str).eq("P4_blocked_auxiliary")).sum()) if not gate.empty else 0,
            "v75_gse121141_old104_mae_weeks": V75_BASELINES["gse121141_old104_mae_weeks"],
            "v86_best_old104plus_delta_mae_weeks": load_v86_summary().get("best_old104plus_delta_mae"),
        },
    }


def write_report(report_path: Path, gate: pd.DataFrame, decision: dict, run_smoke: bool, refresh_candidates_flag: bool) -> None:
    cols = [
        "dataset",
        "candidate_tier",
        "official_metadata_pass",
        "n_old_target_samples",
        "old_target_tissues",
        "n_matched_processed_files",
        "n_selected_smoke_files",
        "schema_smoke_pass",
        "assembly_compatible",
        "smoke_common_regions_estimate",
        "headline_candidate",
        "headline_allowed",
        "minimal_fastq_pilot_candidate",
        "auxiliary_only_reason",
        "next_step",
    ]
    view = gate[[col for col in cols if col in gate.columns]].copy() if not gate.empty else pd.DataFrame()
    lines = [
        "# v9 RALPH Loop Success or Failure Report",
        "",
        f"Date: {utc_now()}",
        "",
        "## Summary",
        "",
        "v9 converts the v8 old-age same-tissue bottleneck into a controlled data strategy loop. This run did not build matrices, train models, or run autoresearch.",
        "",
        f"- Loop decision: `{decision['loop_decision']}`",
        f"- Status: `{decision['status']}`",
        f"- Next action: {decision['next_action']}",
        f"- Candidate refresh used: `{refresh_candidates_flag}`",
        f"- Adapter smoke download used: `{run_smoke}`",
        f"- P1 headline candidates: {decision['metrics']['n_p1_headline_candidates']}",
        f"- P3 minimal FASTQ pilot candidates: {decision['metrics']['n_p3_fastq_pilot_candidates']}",
        "",
        "## Official Network/API Rules",
        "",
        "E-Utils/SOFT are metadata sources. Supplement files are resolved through GEO FTP paths. Individual sample supplements use `/geo/samples/GSM.../<GSM>/suppl/<filename>`; series archives use `/geo/series/GSE.../<GSE>/suppl/`.",
        "",
        md_table(pd.DataFrame([{"source": key, "url": url} for key, url in OFFICIAL_SOURCE_DOCS.items()])),
        "",
        "## Candidate Gate Table",
        "",
        md_table(view, max_rows=30),
        "",
        "## Success and Failure Gates",
        "",
        "- Success requires a headline or minimal-ETL dataset with common 5kb regions >=50000 plus fixed benchmark improvement on GSE121141 old104+ and sanity checks.",
        "- Failure requires two official refresh rounds without P1/P3 data, or matrix/benchmark evidence that old104+ remains around 70w+ despite valid headline data.",
        "- Constrained autoresearch is forbidden until a headline matrix improves old104+ MAE by at least 5w and sanity checks pass.",
        "",
        "## Outputs",
        "",
        "- `results/ralph_v9_loop/ralph_iteration_log.jsonl`",
        "- `results/ralph_v9_loop/ralph_decision_state.json`",
        "- `results/ralph_v9_loop/candidate_gate_table.csv`",
        "- `results/ralph_v9_loop/network_resolution_log.jsonl`",
        "- `results/ralph_v9_loop/candidate_smoke_manifest.csv`",
    ]
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out_dir", default=str(OUT_DIR))
    parser.add_argument("--report", default=str(REPORT))
    parser.add_argument("--refresh_candidates", action="store_true", help="Run official E-Utils/SOFT/FTP candidate refresh.")
    parser.add_argument("--run_smoke", action="store_true", help="Download up to 2-3 small per-sample files per candidate for adapter smoke.")
    parser.add_argument("--retmax", type=int, default=35)
    parser.add_argument("--max_candidates", type=int, default=28)
    parser.add_argument("--refresh_round", type=int, default=1)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    network_log = out_dir / "network_resolution_log.jsonl"
    append_jsonl(
        network_log,
        {
            "timestamp": utc_now(),
            "blocker_type": "official_database_rules_loaded",
            "attempted_url": "",
            "source_doc": " | ".join(OFFICIAL_SOURCE_DOCS.values()),
            "resolution": "metadata_api_and_ftp_roles_are_separated",
            "applied_rule": "eutils_for_metadata_geo_ftp_for_files_sra_ena_for_raw",
        },
    )

    discovery = discovery_module()
    v84 = v84_module()
    v8 = v8_loop_module()
    cached_inventory, cached_samples, cached_supplements = load_cached_inputs()
    frames_inventory = [cached_inventory]
    frames_samples = [cached_samples]
    frames_supplements = [cached_supplements]
    if args.refresh_candidates:
        refreshed = refresh_candidates(discovery, out_dir, args.retmax, args.max_candidates)
        frames_inventory.append(refreshed[0])
        frames_samples.append(refreshed[1])
        frames_supplements.append(refreshed[2])
    inventory = pd.concat([df for df in frames_inventory if not df.empty], ignore_index=True, sort=False) if any(not df.empty for df in frames_inventory) else pd.DataFrame()
    samples = pd.concat([df for df in frames_samples if not df.empty], ignore_index=True, sort=False) if any(not df.empty for df in frames_samples) else pd.DataFrame()
    supplements = pd.concat([df for df in frames_supplements if not df.empty], ignore_index=True, sort=False) if any(not df.empty for df in frames_supplements) else pd.DataFrame()
    if not inventory.empty and "dataset" in inventory:
        inventory = inventory.drop_duplicates("dataset", keep="last").reset_index(drop=True)
    if not samples.empty and {"dataset", "sample_id"}.issubset(samples.columns):
        samples = samples.drop_duplicates(["dataset", "sample_id"], keep="last").reset_index(drop=True)
    if not supplements.empty and {"dataset", "supplement_name"}.issubset(supplements.columns):
        supplements = supplements.drop_duplicates(["dataset", "supplement_name"], keep="last").reset_index(drop=True)

    reference_regions = load_reference_regions(v8)
    gate_rows = []
    smoke_manifest_rows = []
    iteration_log = out_dir / "ralph_iteration_log.jsonl"
    for dataset in candidate_order(inventory, samples)[: args.max_candidates]:
        row, matched, smoke = summarize_candidate(
            dataset,
            inventory,
            samples,
            supplements,
            v84,
            v8,
            reference_regions,
            out_dir,
            args.run_smoke,
            network_log,
        )
        gate_rows.append(row)
        for _, smoke_row in smoke.iterrows():
            smoke_manifest_rows.append(
                {
                    "dataset": dataset,
                    "sample_id": smoke_row.get("sample_id", ""),
                    "supplement_name": smoke_row.get("supplement_name", ""),
                    "status": smoke_row.get("status", ""),
                    "schema_guess": smoke_row.get("schema_guess", ""),
                    "rows_parseable_primary_autosomes": smoke_row.get("rows_parseable_primary_autosomes", ""),
                    "n_smoke_common_regions_with_reference": smoke_row.get("n_smoke_common_regions_with_reference", ""),
                }
            )
        append_jsonl(
            iteration_log,
            {
                "timestamp": utc_now(),
                "iteration": row["iteration"],
                "dataset": dataset,
                "candidate_tier": row["candidate_tier"],
                "hypothesis": row["hypothesis"],
                "action": row["action"],
                "metrics": {
                    "n_old_target_samples": row["n_old_target_samples"],
                    "n_matched_processed_files": row["n_matched_processed_files"],
                    "schema_smoke_pass": row["schema_smoke_pass"],
                    "smoke_common_regions_estimate": row["smoke_common_regions_estimate"],
                },
                "gate_decision": row["next_step"],
                "next_step": row["next_step"],
                "auxiliary_only_reason": row["auxiliary_only_reason"],
            },
        )

    gate = pd.DataFrame(gate_rows)
    smoke_manifest = pd.DataFrame(smoke_manifest_rows)
    gate.to_csv(out_dir / "candidate_gate_table.csv", index=False)
    smoke_manifest.to_csv(out_dir / "candidate_smoke_manifest.csv", index=False)
    decision = decide_state(gate, args.refresh_round, dry_run=not args.run_smoke)
    write_json(out_dir / "ralph_decision_state.json", decision)
    write_report(Path(args.report), gate, decision, args.run_smoke, args.refresh_candidates)
    print(json.dumps(decision, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

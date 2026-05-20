#!/usr/bin/env python3
"""v16 Route B candidate refresh and priority pilot controller.

The controller refreshes public mouse methylation candidates from official
GEO/E-Utils/SOFT/FTP plus SRA/ENA run metadata, builds a gate table and a pilot
priority queue, and records RALPH decisions. It does not train, run
autoresearch, download FASTQ, or run Bismark during refresh/gate modes.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")
OUT_DIR = ROOT / "results" / "ralph_v16_route_b_candidate_refresh"
DOC_REPORT = ROOT / "doc" / "20_analysis" / "46_20260519_v16_route_b_candidate_refresh_report.md"
PROJECT_INDEX = ROOT / "doc" / "00_meta" / "02_20260519_project_status_index_v13.md"

COMMON_REGION_GATE = 50_000
AGE_COVERAGE_GATE = 0.95
OLD_THRESHOLD_WEEKS = 104.0
TARGET_TISSUES = {"brain_cortex", "cortex", "heart", "lung"}

BASE_SEED_ACCESSIONS = [
    "GSE286302",
    "GSE281602",
    "GSE304754",
    "GSE292803",
    "GSE313770",
    "GSE225166",
    "GSE83947",
    "GSE134398",
    "GSE232547",
    "GSE290397",
    "GSE184267",
    "GSE325694",
    "GSE160074",
    "GSE312263",
    "GSE153331",
    "GSE232546",
    "GSE213628",
    "GSE80672",
]

SEARCH_TERMS = [
    'Mus musculus[ORGN] AND "Methylation profiling by high throughput sequencing" AND (RRBS OR WGBS OR bisulfite) AND ("brain cortex" OR cortex OR heart OR lung) AND (aging OR aged OR old OR "24 month" OR "26 month" OR "28 month") AND gse[ETYP]',
    'Mus musculus[ORGN] AND "DNA methylation" AND ("brain cortex" OR cortex OR heart OR lung) AND ("24 month" OR "26 month" OR "28 month" OR aged OR old) AND gse[ETYP]',
    'Mus musculus[ORGN] AND "reduced representation bisulfite" AND (heart OR lung OR cortex OR "brain cortex") AND (aging OR aged OR old) AND gse[ETYP]',
    'Mus musculus[ORGN] AND "whole genome bisulfite" AND (heart OR lung OR cortex OR "brain cortex") AND (aging OR aged OR old) AND gse[ETYP]',
    'Mus musculus[ORGN] AND WGBS AND (heart OR lung OR cortex OR "brain cortex") AND ("24 month" OR "26 month" OR "28 month") AND gse[ETYP]',
]

OFFICIAL_SOURCE_DOCS = {
    "geo_download": "https://www.ncbi.nlm.nih.gov/geo/info/download.html",
    "geo_programmatic_access": "https://www.ncbi.nlm.nih.gov/geo/info/geo_paccess.html",
    "geo_soft": "https://www.ncbi.nlm.nih.gov/geo/info/soft.html",
    "sra_download": "https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/",
    "ena_file_reports": "https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/file-reports.html",
}

DO_NOT_REPILOT_HEADLINE = {
    "GSE83947": "processed and raw Route B pilots already failed common-region gate",
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def clean(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    return str(value)


def numeric(value: Any, default: float = 0.0) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return default
    return parsed if math.isfinite(parsed) else default


def boolish(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    if isinstance(value, (int, float)):
        return bool(value)
    return str(value).strip().lower() in {"true", "1", "yes", "y", "pass", "passed"}


def import_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def v14_module():
    return import_module(ROOT / "scripts" / "validate" / "run_v14_public_data_rescue.py", "route_b_v14_helpers_for_v16")


def discovery_module():
    module = v14_module()
    return module.discovery_module()


def read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str), encoding="utf-8")


def append_jsonl(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(payload, sort_keys=True, default=str) + "\n")


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = []
        for row in rows:
            for key in row:
                if key not in fieldnames:
                    fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def md_table(rows: list[dict[str, Any]], columns: list[str], limit: int | None = None) -> str:
    if not rows:
        return "No rows."
    view = rows[:limit] if limit else rows
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for row in view:
        lines.append("| " + " | ".join(str(row.get(col, "")) for col in columns) + " |")
    return "\n".join(lines)


def collect_seed_accessions(extra_accessions: list[str]) -> list[str]:
    accessions = list(dict.fromkeys([*BASE_SEED_ACCESSIONS, *extra_accessions]))
    for path in [
        ROOT / "results" / "ralph_v14_public_data_rescue" / "candidate_gate_table.csv",
        ROOT / "results" / "ralph_v15_matrix_gate" / "candidate_gate_table.csv",
        ROOT / "results" / "ralph_v12_loop" / "candidate_gate_table.csv",
    ]:
        for row in read_csv(path):
            dataset = row.get("dataset") or row.get("accession") or ""
            if dataset.startswith("GSE") and dataset not in accessions:
                accessions.append(dataset)
    return accessions


def run_dry_refresh(out_dir: Path, accessions: list[str]) -> dict[str, Any]:
    refresh_rows = []
    for accession in accessions:
        refresh_rows.append(
            {
                "dataset": accession,
                "refresh_mode": "dry_run",
                "official_geo_url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}",
                "source_terms": "manual_seed_or_prior_gate",
                "download_authorized": False,
                "bismark_authorized": False,
                "training_authorized": False,
                "autoresearch_authorized": False,
            }
        )
    for idx, term in enumerate(SEARCH_TERMS, 1):
        refresh_rows.append(
            {
                "dataset": f"SEARCH_TERM_{idx}",
                "refresh_mode": "dry_run_query",
                "official_geo_url": "https://www.ncbi.nlm.nih.gov/gds",
                "source_terms": term,
                "download_authorized": False,
                "bismark_authorized": False,
                "training_authorized": False,
                "autoresearch_authorized": False,
            }
        )
    write_csv(out_dir / "candidate_refresh_table.csv", refresh_rows)
    state = {
        "timestamp": utc_now(),
        "loop_version": "v16_route_b_candidate_refresh",
        "mode": "refresh",
        "dry_run": True,
        "decision": "dry_run_refresh_plan_created_no_network_no_download",
        "n_seed_accessions": len(accessions),
        "n_search_terms": len(SEARCH_TERMS),
        "download_authorized": False,
        "bismark_authorized": False,
        "training_authorized": False,
        "autoresearch_authorized": False,
    }
    write_json(out_dir / "ralph_decision_state.json", state)
    return state


def local_v15_attempts() -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    for row in read_csv(ROOT / "results" / "ralph_v15_matrix_gate" / "matrix_gate_attempts.csv"):
        dataset = row.get("dataset", "")
        if dataset:
            out[dataset] = row
    return out


def candidate_failed_prior(dataset: str) -> tuple[bool, str]:
    if dataset in DO_NOT_REPILOT_HEADLINE:
        return True, DO_NOT_REPILOT_HEADLINE[dataset]
    v15 = local_v15_attempts().get(dataset)
    if not v15:
        return False, ""
    if v15.get("headline_matrix_gate_pass") == "False" and v15.get("gate_status") in {
        "matrix_gate_failed",
        "research_or_audit_only_not_matrix_ready",
        "auxiliary_matrix_gate_passed_not_headline",
    }:
        return True, v15.get("blocker_type", "")
    return False, ""


def classify_v16(candidate: dict[str, Any]) -> tuple[str, str, str, bool]:
    dataset = clean(candidate.get("dataset"))
    age_fraction = numeric(candidate.get("age_known_fraction"), 0.0)
    age_pass = age_fraction >= AGE_COVERAGE_GATE
    old_exact_n = int(numeric(candidate.get("old_exact_target_n"), 0))
    old_label_n = int(numeric(candidate.get("old_label_target_n"), 0))
    exact_target = old_exact_n > 0
    label_target = old_label_n > 0
    bulk_pass = not clean(candidate.get("context_blockers"))
    parseable = boolish(candidate.get("parseable_processed_schema"))
    processed_files = boolish(candidate.get("has_sample_processed_files"))
    common_regions = numeric(candidate.get("local_common_regions"), 0)
    common_pass = common_regions >= COMMON_REGION_GATE
    failed_prior, prior_reason = candidate_failed_prior(dataset)
    old_tissues = clean(candidate.get("old_exact_target_tissues")) or clean(candidate.get("old_label_target_tissues"))

    if failed_prior and dataset == "GSE83947":
        return "P2_auxiliary", "do_not_repilot_headline", f"prior_failed_candidate:{prior_reason}", False

    if age_pass and exact_target and bulk_pass and parseable and common_pass:
        return "P1_processed_headline", "ready_for_ralph_learn_pending_explicit_training_approval", "", True
    if age_pass and exact_target and bulk_pass and processed_files and not common_pass:
        return "P1_processed_headline", "processed_adapter_smoke_or_matrix_gate_required", "processed_files_present_common_regions_unknown_or_failed", False
    if age_pass and exact_target and bulk_pass and not processed_files:
        return "P3_raw_pilot", "raw_pilot_candidate_pending_priority_review", "no_processed_methylation_supplement_detected", False

    if exact_target or label_target or old_tissues:
        reasons = []
        if not age_pass:
            reasons.append("sample_specific_exact_age_coverage_lt_95pct")
        if not exact_target and label_target:
            reasons.append("old_age_group_label_not_exact_age")
        if not bulk_pass:
            reasons.append("non_bulk_or_cell_context")
        if processed_files and not parseable:
            reasons.append("processed_schema_unknown")
        if common_regions and not common_pass:
            reasons.append("common_region_gate_failed")
        if failed_prior:
            reasons.append(f"prior_failed_candidate:{prior_reason}")
        return "P2_auxiliary", "keep_for_diagnostic_or_adapter_evidence_not_headline", ";".join(reasons), False

    reasons = []
    if not exact_target:
        reasons.append("no_exact_old_target_tissue_samples")
    if not label_target:
        reasons.append("no_old_target_tissue_label")
    if not bulk_pass:
        reasons.append("blocked_context")
    return "P4_blocked", "do_not_promote_to_headline", ";".join(reasons), False


def priority_score(row: dict[str, Any]) -> tuple:
    tier = clean(row.get("tier"))
    tissue_text = clean(row.get("old_exact_target_tissues")) or clean(row.get("old_label_target_tissues"))
    if "brain_cortex" in tissue_text or "cortex" in tissue_text:
        tissue_rank = 0
    elif "heart" in tissue_text:
        tissue_rank = 1
    elif "lung" in tissue_text:
        tissue_rank = 2
    else:
        tissue_rank = 9
    tier_rank = {"P1_processed_headline": 0, "P3_raw_pilot": 1}.get(tier, 9)
    failed_prior = 1 if clean(row.get("prior_failure_reason")) else 0
    processed_rank = 0 if boolish(row.get("has_sample_processed_files")) else 1
    bytes_value = numeric(row.get("ena_total_fastq_bytes"), 10**18)
    if bytes_value <= 0:
        bytes_value = 10**18
    return (
        tier_rank,
        tissue_rank,
        failed_prior,
        processed_rank,
        -int(numeric(row.get("old_exact_target_n"), 0)),
        -numeric(row.get("age_known_fraction"), 0),
        bytes_value,
        clean(row.get("dataset")),
    )


def select_processed_pilot_rows(module: Any, dataset: str, samples: pd.DataFrame, supplements: pd.DataFrame, tier: str) -> list[dict[str, Any]]:
    rows = module.select_processed_pilot_rows(dataset, samples, supplements, "P3_raw_pilot" if tier == "P3_raw_pilot" else tier)
    for row in rows:
        row["download_authorized"] = False
        row["bismark_authorized"] = False
        row["training_authorized"] = False
        row["autoresearch_authorized"] = False
    return rows


def build_refresh_and_gate(
    out_dir: Path,
    accessions: list[str],
    retmax: int,
    max_candidates: int,
    no_eutils_search: bool,
    refresh_soft: bool,
    no_head: bool,
    skip_ena: bool,
    skip_sra: bool,
) -> dict[str, Any]:
    helper = v14_module()
    module = helper.discovery_module()
    selected, source_terms, summaries = helper.select_candidate_accessions(
        module,
        accessions=accessions,
        no_eutils_search=no_eutils_search,
        retmax=retmax,
        max_candidates=max_candidates,
    )
    # The search result limit should not drop pre-registered seeds, especially
    # prior failed candidates such as GSE83947 that must remain explicitly
    # blocked from new headline pilots.
    for accession in accessions:
        if accession not in selected:
            selected.append(accession)
            source_terms.setdefault(accession, ["manual_seed_or_prior_gate"])
    network_log = out_dir / "network_resolution_log.jsonl"
    if network_log.exists():
        network_log.unlink()

    refresh_rows: list[dict[str, Any]] = []
    gate_rows: list[dict[str, Any]] = []
    pilot_rows: list[dict[str, Any]] = []
    smoke_rows: list[dict[str, Any]] = []

    for accession in selected:
        try:
            row, supplements, samples = module.build_candidate(
                accession,
                source_terms=source_terms.get(accession, ["v16_search"]),
                eutils_summary=summaries.get(accession),
                refresh_soft=refresh_soft,
                do_head=not no_head,
            )
        except Exception as exc:
            row = {
                "dataset": accession,
                "geo_url": module.geo_url(accession),
                "title": "",
                "summary": "",
                "overall_design": "",
                "processed_schema_guess": "",
                "blockers": str(exc)[:500],
                "recommended_action": "inspect_manually",
            }
            supplements = []
            samples = []

        supplement_df = pd.DataFrame(supplements)
        sample_df = pd.DataFrame(samples)
        support = helper.summarize_sample_support(sample_df)
        blockers = helper.context_blockers(row, sample_df)
        local = helper.local_conversion_state(accession)
        bioprojects, bioproject_error = helper.extract_bioprojects_from_soft(module, accession)
        parseable = helper.is_parseable_processed(clean(row.get("processed_schema_guess")), supplement_df)
        processed_files = helper.has_sample_processed_files(supplement_df)
        failed_prior, prior_reason = candidate_failed_prior(accession)
        candidate = {
            **row,
            **support,
            **local,
            "dataset": accession,
            "bioprojects": ";".join(bioprojects),
            "bioproject_error": bioproject_error,
            "context_blockers": ";".join(blockers),
            "bulk_context_pass": not blockers,
            "parseable_processed_schema": parseable,
            "has_sample_processed_files": processed_files,
            "official_geo_url": module.geo_url(accession),
            "v16_source_terms": " || ".join(source_terms.get(accession, [])),
            "prior_failure_reason": prior_reason if failed_prior else "",
            "do_not_repilot_headline": accession in DO_NOT_REPILOT_HEADLINE,
            "download_authorized": False,
            "training_authorized": False,
            "bismark_authorized": False,
            "autoresearch_authorized": False,
        }
        tier, gate_status, blocker_type, headline_allowed = classify_v16(candidate)
        candidate.update(
            {
                "tier": tier,
                "gate_status": gate_status,
                "blocker_type": blocker_type or candidate.get("context_blockers") or candidate.get("local_reason") or candidate.get("blockers", ""),
                "headline_allowed": headline_allowed,
            }
        )

        if tier == "P4_blocked":
            run_rows = []
            network_summary = {
                "ena_n_runs": 0,
                "ena_total_fastq_bytes": 0,
                "sra_n_runs": 0,
                "sra_total_size_mb": 0.0,
            }
            append_jsonl(
                network_log,
                {
                    "timestamp": utc_now(),
                    "dataset": accession,
                    "blocker_type": "p4_blocked_skip_runinfo",
                    "attempted_url": "",
                    "source_doc": OFFICIAL_SOURCE_DOCS["geo_soft"],
                    "resolution": "soft_gate_sufficient_p4_never_pilot_skip_sra_ena_runinfo",
                    "applied_rule": "no_raw_run_metadata_for_p4_blocked_candidates",
                    "error": "",
                },
            )
        else:
            run_rows, network_summary = helper.build_network_and_pilot_rows(
                accession,
                bioprojects,
                tier,
                out_dir,
                network_log,
                query_sra=not skip_sra,
                query_ena=not skip_ena,
            )
        processed_rows = select_processed_pilot_rows(helper, accession, sample_df, supplement_df, tier)
        candidate.update(network_summary)
        refresh_rows.append(candidate)

        gate_row = {
            "dataset": accession,
            "tier": tier,
            "gate_status": gate_status,
            "blocker_type": candidate["blocker_type"],
            "recommended_action": gate_status,
            "headline_allowed": headline_allowed,
            "sample_specific_age_pass": numeric(candidate.get("age_known_fraction"), 0) >= AGE_COVERAGE_GATE,
            "target_tissue_pass": int(numeric(candidate.get("old_exact_target_n"), 0)) > 0,
            "bulk_context_pass": candidate["bulk_context_pass"],
            "schema_smoke_pass": bool(candidate.get("local_conversion_status")),
            "assembly_traceable": bool(candidate.get("local_common_regions")),
            "common_regions_pass": bool(candidate.get("local_common_region_gate")),
            "common_regions_estimate": candidate.get("local_common_regions", ""),
            "age_known_fraction": candidate.get("age_known_fraction", ""),
            "old_exact_target_n": candidate.get("old_exact_target_n", ""),
            "old_label_target_n": candidate.get("old_label_target_n", ""),
            "old_exact_target_tissues": candidate.get("old_exact_target_tissues", ""),
            "old_label_target_tissues": candidate.get("old_label_target_tissues", ""),
            "target_tissue_n": candidate.get("target_tissue_n", ""),
            "tissue_counts": candidate.get("tissue_counts", ""),
            "context_blockers": candidate.get("context_blockers", ""),
            "processed_schema_guess": candidate.get("processed_schema_guess", ""),
            "parseable_processed_schema": parseable,
            "has_sample_processed_files": processed_files,
            "preferred_url": candidate.get("preferred_url", ""),
            "bioprojects": candidate.get("bioprojects", ""),
            "ena_n_runs": candidate.get("ena_n_runs", 0),
            "ena_total_fastq_bytes": candidate.get("ena_total_fastq_bytes", 0),
            "sra_n_runs": candidate.get("sra_n_runs", 0),
            "local_conversion_status": candidate.get("local_conversion_status", ""),
            "local_conversion_state_path": candidate.get("local_conversion_state_path", ""),
            "local_reason": candidate.get("local_reason", ""),
            "prior_failure_reason": candidate.get("prior_failure_reason", ""),
            "download_authorized": False,
            "bismark_authorized": False,
            "training_authorized": False,
            "autoresearch_authorized": False,
        }
        gate_rows.append(gate_row)

        candidate_pilot_rows = processed_rows or run_rows
        for pilot in candidate_pilot_rows:
            pilot["tier"] = tier
            pilot["gate_status"] = gate_status
            pilot["candidate_recommended_action"] = gate_status
            pilot["download_authorized"] = False
            pilot["bismark_authorized"] = False
            pilot["training_authorized"] = False
            pilot["autoresearch_authorized"] = False
        pilot_rows.extend(candidate_pilot_rows)
        smoke_rows.append(
            {
                "dataset": accession,
                "local_conversion_status": candidate.get("local_conversion_status", ""),
                "local_conversion_state_path": candidate.get("local_conversion_state_path", ""),
                "local_common_regions": candidate.get("local_common_regions", ""),
                "local_common_region_gate": candidate.get("local_common_region_gate", False),
                "local_beta_min": candidate.get("local_beta_min", ""),
                "local_beta_max": candidate.get("local_beta_max", ""),
                "local_n_samples": candidate.get("local_n_samples", ""),
                "local_matrix_path": candidate.get("local_matrix_path", ""),
                "schema_smoke_pass": bool(candidate.get("local_conversion_status")),
            }
        )
        time.sleep(0.05)

    priority_rows = [
        row
        for row in gate_rows
        if row["tier"] in {"P1_processed_headline", "P3_raw_pilot"}
        and row["dataset"] not in DO_NOT_REPILOT_HEADLINE
        and row.get("gate_status") != "ready_for_ralph_learn_pending_explicit_training_approval"
    ]
    priority_rows.sort(key=priority_score)
    for idx, row in enumerate(priority_rows, 1):
        row["pilot_priority_rank"] = idx
        row["pilot_scope"] = "2_3_samples_only"
        row["pilot_authorization_state"] = "pending_explicit_authorize_pilot"

    write_csv(out_dir / "candidate_refresh_table.csv", refresh_rows)
    write_csv(out_dir / "candidate_gate_table.csv", gate_rows)
    write_csv(out_dir / "pilot_priority_queue.csv", priority_rows)
    write_csv(out_dir / "pilot_run_manifest.csv", pilot_rows)
    write_csv(out_dir / "candidate_smoke_manifest.csv", smoke_rows)

    # Reuse v15 matrix-gate attempts as historical strike evidence for the new refresh.
    v15_attempts = read_csv(ROOT / "results" / "ralph_v15_matrix_gate" / "matrix_gate_attempts.csv")
    v15_strikes = read_csv(ROOT / "results" / "ralph_v15_matrix_gate" / "three_strike_table.csv")
    write_csv(out_dir / "matrix_gate_attempts.csv", v15_attempts)
    write_csv(out_dir / "three_strike_table.csv", v15_strikes)

    n_headline_ready = sum(1 for row in gate_rows if row["headline_allowed"])
    n_priority = len(priority_rows)
    n_p1 = sum(1 for row in gate_rows if row["tier"] == "P1_processed_headline")
    n_p3 = sum(1 for row in gate_rows if row["tier"] == "P3_raw_pilot")
    n_p2 = sum(1 for row in gate_rows if row["tier"] == "P2_auxiliary")
    n_p4 = sum(1 for row in gate_rows if row["tier"] == "P4_blocked")

    if n_headline_ready:
        decision = "ready_for_ralph_learn_pending_explicit_training_approval"
        next_action = "prepare_guarded_fixed_benchmark_manifest_pending_explicit_training_approval"
    elif n_priority:
        decision = "pilot_candidates_prioritized_pending_explicit_single_candidate_authorization"
        next_action = "run_mode_pilot_for_top_ranked_candidate_only_after_explicit_authorize_pilot"
    else:
        decision = "no_new_route_b_pilot_candidate_keep_route_a_or_redefine_search"
        next_action = "return_to_route_a_or_research_query_refinement_do_not_train"

    state = {
        "timestamp": utc_now(),
        "loop_version": "v16_route_b_candidate_refresh",
        "mode": "gate",
        "decision": decision,
        "next_action": next_action,
        "target_success_state": "ready_for_ralph_learn_pending_explicit_training_approval",
        "metrics": {
            "n_candidates": len(gate_rows),
            "n_p1_processed_headline": n_p1,
            "n_p3_raw_pilot": n_p3,
            "n_p2_auxiliary": n_p2,
            "n_p4_blocked": n_p4,
            "n_priority_pilot_candidates": n_priority,
            "n_headline_ready": n_headline_ready,
            "common_region_gate": COMMON_REGION_GATE,
            "old_threshold_weeks": OLD_THRESHOLD_WEEKS,
        },
        "priority_candidates": [row["dataset"] for row in priority_rows[:10]],
        "download_authorized": False,
        "bismark_authorized": False,
        "training_authorized": False,
        "autoresearch_authorized": False,
        "outputs": {
            "candidate_refresh_table": str((out_dir / "candidate_refresh_table.csv").relative_to(ROOT)),
            "candidate_gate_table": str((out_dir / "candidate_gate_table.csv").relative_to(ROOT)),
            "pilot_priority_queue": str((out_dir / "pilot_priority_queue.csv").relative_to(ROOT)),
            "pilot_run_manifest": str((out_dir / "pilot_run_manifest.csv").relative_to(ROOT)),
            "matrix_gate_attempts": str((out_dir / "matrix_gate_attempts.csv").relative_to(ROOT)),
            "three_strike_table": str((out_dir / "three_strike_table.csv").relative_to(ROOT)),
            "network_resolution_log": str((out_dir / "network_resolution_log.jsonl").relative_to(ROOT)),
        },
        "official_source_docs": OFFICIAL_SOURCE_DOCS,
    }
    write_json(out_dir / "ralph_decision_state.json", state)
    append_jsonl(
        out_dir / "ralph_iteration_log.jsonl",
        {
            "timestamp": utc_now(),
            "phase": "RAP",
            "action": "route_b_official_candidate_refresh_and_gate_no_download_no_training",
            "decision": decision,
            "metrics": state["metrics"],
            "download_authorized": False,
            "bismark_authorized": False,
            "training_authorized": False,
            "autoresearch_authorized": False,
        },
    )
    write_report(out_dir, gate_rows, priority_rows, state)
    update_project_index()
    return state


def write_report(out_dir: Path, gate_rows: list[dict[str, Any]], priority_rows: list[dict[str, Any]], state: dict[str, Any]) -> None:
    tier_counts: dict[str, int] = {}
    for row in gate_rows:
        tier_counts[row["tier"]] = tier_counts.get(row["tier"], 0) + 1
    tier_rows = [{"tier": key, "count": value} for key, value in sorted(tier_counts.items())]
    gate_columns = [
        "dataset",
        "tier",
        "age_known_fraction",
        "old_exact_target_n",
        "old_exact_target_tissues",
        "context_blockers",
        "common_regions_estimate",
        "gate_status",
        "blocker_type",
    ]
    priority_columns = [
        "pilot_priority_rank",
        "dataset",
        "tier",
        "old_exact_target_tissues",
        "old_exact_target_n",
        "has_sample_processed_files",
        "ena_total_fastq_bytes",
        "gate_status",
        "blocker_type",
    ]
    text = f"""# v16 Route B Candidate Refresh Report

Date: {utc_now()}

## Summary

v16 refreshed Route B public candidates using official GEO/E-Utils/SOFT/FTP,
SRA RunInfo, and ENA file report metadata. This run did not download FASTQ,
run Bismark, train models, or start autoresearch.

Decision: `{state['decision']}`.

Next action: `{state['next_action']}`.

## Tier Counts

{md_table(tier_rows, ['tier', 'count'])}

## Candidate Gate Table

{md_table(gate_rows, gate_columns, limit=80)}

## Pilot Priority Queue

{md_table(priority_rows, priority_columns, limit=30)}

## Rules Applied

- `GSE83947` remains `do_not_repilot_headline` because processed and raw pilots
  already failed the common-region gate.
- P2 candidates are diagnostic only and cannot enter headline pilot.
- P4 candidates are blocked and cannot enter pilot.
- Any pilot remains limited to one candidate and 2-3 samples after explicit
  `--authorize-pilot`.

## Guardrails

- Download authorized: `false`
- Bismark authorized: `false`
- Training authorized: `false`
- Autoresearch authorized: `false`
- Human clock CpG mapping: forbidden
- Dummy AUC: forbidden

## Outputs

- `{state['outputs']['candidate_refresh_table']}`
- `{state['outputs']['candidate_gate_table']}`
- `{state['outputs']['pilot_priority_queue']}`
- `{state['outputs']['pilot_run_manifest']}`
- `{state['outputs']['network_resolution_log']}`

"""
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "v16_route_b_candidate_refresh_report.md").write_text(text, encoding="utf-8")
    DOC_REPORT.parent.mkdir(parents=True, exist_ok=True)
    DOC_REPORT.write_text(text, encoding="utf-8")


def update_project_index() -> None:
    link = "- v16 Route B candidate refresh report: `doc/20_analysis/46_20260519_v16_route_b_candidate_refresh_report.md`"
    if not PROJECT_INDEX.exists():
        return
    text = PROJECT_INDEX.read_text(encoding="utf-8")
    if link in text:
        return
    marker = "- v15 matrix-gate RALPH loop report: `doc/20_analysis/45_20260519_v15_matrix_gate_ralph_loop_report.md`"
    if marker in text:
        text = text.replace(marker, marker + "\n" + link)
    else:
        text += "\n" + link + "\n"
    PROJECT_INDEX.write_text(text, encoding="utf-8")


def run_pilot_mode(out_dir: Path, candidate: str, authorize_pilot: bool) -> dict[str, Any]:
    gate_rows = read_csv(out_dir / "candidate_gate_table.csv")
    priority_rows = read_csv(out_dir / "pilot_priority_queue.csv")
    gate = next((row for row in gate_rows if row.get("dataset") == candidate), None)
    priority = next((row for row in priority_rows if row.get("dataset") == candidate), None)
    if gate is None:
        decision = "candidate_missing_run_gate_first"
        reason = "candidate_not_found_in_gate_table"
    elif candidate in DO_NOT_REPILOT_HEADLINE:
        decision = "pilot_blocked_prior_failed_candidate"
        reason = DO_NOT_REPILOT_HEADLINE[candidate]
    elif priority is None:
        decision = "pilot_blocked_candidate_not_in_priority_queue"
        reason = gate.get("blocker_type", "")
    elif not authorize_pilot:
        decision = "pilot_manifest_ready_pending_explicit_authorize_pilot"
        reason = "authorization_flag_missing"
    else:
        decision = "pilot_authorized_but_execution_requires_candidate_specific_adapter"
        reason = "v16 controller does not batch-download; run candidate-specific processed adapter or Route B pilot script after reviewing manifest"

    pilot_manifest = [
        row for row in read_csv(out_dir / "pilot_run_manifest.csv") if row.get("dataset") == candidate
    ][:3]
    for row in pilot_manifest:
        row["selected_for_pilot"] = bool(priority) and candidate not in DO_NOT_REPILOT_HEADLINE
        row["authorize_pilot_flag"] = authorize_pilot
        row["download_authorized"] = False
        row["bismark_authorized"] = False
        row["training_authorized"] = False
        row["autoresearch_authorized"] = False
    write_csv(out_dir / f"{candidate}_pilot_selection_manifest.csv", pilot_manifest)
    state = {
        "timestamp": utc_now(),
        "loop_version": "v16_route_b_candidate_refresh",
        "mode": "pilot",
        "candidate": candidate,
        "decision": decision,
        "reason": reason,
        "n_manifest_rows": len(pilot_manifest),
        "download_authorized": False,
        "bismark_authorized": False,
        "training_authorized": False,
        "autoresearch_authorized": False,
        "pilot_selection_manifest": str((out_dir / f"{candidate}_pilot_selection_manifest.csv").relative_to(ROOT)),
    }
    write_json(out_dir / f"{candidate}_pilot_decision_state.json", state)
    append_jsonl(out_dir / "ralph_iteration_log.jsonl", state)
    return state


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=["refresh", "gate", "pilot"], default="gate")
    parser.add_argument("--out-dir", default=str(OUT_DIR))
    parser.add_argument("--accessions", default="")
    parser.add_argument("--retmax", type=int, default=40)
    parser.add_argument("--max-candidates", type=int, default=60)
    parser.add_argument("--no-eutils-search", action="store_true")
    parser.add_argument("--refresh-soft", action="store_true")
    parser.add_argument("--no-head", action="store_true")
    parser.add_argument("--skip-ena", action="store_true")
    parser.add_argument("--skip-sra", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--candidate", default="")
    parser.add_argument("--authorize-pilot", action="store_true")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    extras = [item.strip() for item in args.accessions.split(",") if item.strip()]
    accessions = collect_seed_accessions(extras)

    if args.mode == "refresh" and args.dry_run:
        state = run_dry_refresh(out_dir, accessions)
    elif args.mode in {"refresh", "gate"}:
        state = build_refresh_and_gate(
            out_dir=out_dir,
            accessions=accessions,
            retmax=args.retmax,
            max_candidates=args.max_candidates,
            no_eutils_search=args.no_eutils_search,
            refresh_soft=args.refresh_soft,
            no_head=args.no_head,
            skip_ena=args.skip_ena,
            skip_sra=args.skip_sra,
        )
        if args.mode == "refresh":
            state["mode"] = "refresh"
            write_json(out_dir / "ralph_decision_state.json", state)
    else:
        if not args.candidate:
            raise SystemExit("--candidate is required for --mode pilot")
        state = run_pilot_mode(out_dir, args.candidate, args.authorize_pilot)

    print(json.dumps(state, indent=2, sort_keys=True, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

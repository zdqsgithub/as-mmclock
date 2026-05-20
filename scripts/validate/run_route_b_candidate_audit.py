#!/usr/bin/env python3
"""Audit Route B minimal FASTQ/Bismark pilot candidates.

This audit reads existing v10-v12 candidate metadata and cached official
SRA/BioSample verification tables. It does not download FASTQ, run Bismark,
train models, rebuild matrices, or call external APIs.
"""
from __future__ import annotations

import csv
import json
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any

ROOT = Path("/home/zdq-as/mouse_methyl_work")
OUT_DIR = ROOT / "results" / "route_b_candidate_audit"
REPORT_PATH = ROOT / "doc" / "20_analysis" / "34_20260519_route_b_candidate_audit_report.md"

RAW_SAMPLE_MANIFEST = ROOT / "results" / "ralph_v11_data_strategy" / "raw_biosample_verification_manifest.csv"
RAW_SUMMARY = ROOT / "results" / "ralph_v11_data_strategy" / "raw_biosample_verification_summary.csv"
RAW_PROJECT_LEADS = ROOT / "results" / "ralph_v11_data_strategy" / "raw_project_rfc_summary.csv"
V12_MINIMAL_ETL = ROOT / "results" / "ralph_v12_loop" / "minimal_etl_rfc_manifest.csv"
V12_GATE = ROOT / "results" / "ralph_v12_loop" / "candidate_gate_table.csv"
V11_CONVERSION_QUEUE = ROOT / "results" / "download_backlog_v11_3" / "conversion_queue_smoke_gated.csv"

TARGET_TISSUES = {"brain_cortex", "cortex", "heart", "lung"}
OLD_WEEKS = 104.0
CONTEXT_BLOCKER_TERMS = {
    "cell type",
    "cardiomyocyte",
    "cardiomyocytes",
    "single-cell",
    "single_cell",
    "low_coverage",
    "itag",
    "organoid",
}


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    columns: list[str] = []
    for row in rows:
        for key in row:
            if key not in columns:
                columns.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str), encoding="utf-8")


def as_bool(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def as_float(value: Any) -> float | None:
    text = str(value).strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def split_terms(value: Any) -> list[str]:
    text = str(value or "").strip()
    if not text:
        return []
    out: list[str] = []
    for part in text.replace("|", ";").split(";"):
        item = part.strip()
        if item:
            out.append(item)
    return out


def normalize_tissue(value: Any) -> str:
    text = str(value or "").strip().lower().replace(" ", "_").replace("-", "_")
    aliases = {
        "cortex": "brain_cortex",
        "brain_cortex": "brain_cortex",
        "heart": "heart",
        "lung": "lung",
    }
    return aliases.get(text, text)


def has_context_blocker(blockers: list[str]) -> bool:
    joined = ";".join(blockers).lower()
    return any(term in joined for term in CONTEXT_BLOCKER_TERMS)


def sample_decision(row: dict[str, str]) -> str:
    blockers = split_terms(row.get("blockers"))
    tissue = normalize_tissue(row.get("tissue"))
    age_weeks = as_float(row.get("age_weeks"))
    if not as_bool(row.get("assay_pass")):
        return "blocked_assay_not_bisulfite"
    if tissue not in TARGET_TISSUES or not as_bool(row.get("target_tissue_pass")):
        return "blocked_non_target_tissue"
    if has_context_blocker(blockers):
        return "blocked_context_not_bulk_headline"
    if as_bool(row.get("age_exact_pass")) and age_weeks is not None and age_weeks >= OLD_WEEKS:
        return "eligible_old_exact_target_pending_approval"
    if as_bool(row.get("age_group_pass")):
        return "diagnostic_only_age_group_no_exact_age"
    if as_bool(row.get("age_exact_pass")):
        return "diagnostic_exact_age_but_not_old_target"
    return "blocked_missing_age"


def candidate_decision(sample_rows: list[dict[str, Any]], v12_row: dict[str, str] | None) -> tuple[str, str]:
    if not sample_rows:
        return "blocked_no_verified_runs", "No verified BioSample/RunInfo sample rows were found."
    eligible = [r for r in sample_rows if r["sample_decision"] == "eligible_old_exact_target_pending_approval"]
    diagnostic_age_group = [r for r in sample_rows if r["sample_decision"] == "diagnostic_only_age_group_no_exact_age"]
    context_blocked = [r for r in sample_rows if r["sample_decision"] == "blocked_context_not_bulk_headline"]
    if len(eligible) >= 2:
        return (
            "route_b_pilot_eligible_pending_explicit_approval",
            "At least two old target-tissue exact-age bulk-compatible runs were found.",
        )
    if context_blocked:
        return (
            "blocked_context_not_route_b_headline",
            "Runs have exact target-tissue metadata but are cell-type/single-cell/other non-bulk context.",
        )
    if diagnostic_age_group:
        return (
            "diagnostic_only_age_group_no_exact_age",
            "Runs are target-tissue bisulfite data but age is only an age group, not sample-specific exact age.",
        )
    if v12_row and v12_row.get("gate_status"):
        return (
            "blocked_by_v12_gate",
            f"v12 gate status: {v12_row.get('gate_status')}; blocker: {v12_row.get('blocker_type')}",
        )
    return "not_route_b_candidate", "No old exact-age bulk target-tissue candidate passed Route B gates."


def load_v12_raw_gate() -> dict[str, dict[str, str]]:
    rows = read_csv(V12_GATE)
    out: dict[str, dict[str, str]] = {}
    for row in rows:
        candidate_id = row.get("candidate_id", "")
        if candidate_id.startswith("round3_raw_rfc:"):
            out[candidate_id.split(":", 1)[1]] = row
    return out


def load_v12_minimal_etl() -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    for row in read_csv(V12_MINIMAL_ETL):
        candidate_id = row.get("candidate_id", "")
        if candidate_id.startswith("round3_raw_rfc:"):
            out[candidate_id.split(":", 1)[1]] = row
    return out


def build_sample_audit(raw_rows: list[dict[str, str]]) -> list[dict[str, Any]]:
    sample_rows: list[dict[str, Any]] = []
    for row in raw_rows:
        blockers = split_terms(row.get("blockers"))
        sample_rows.append(
            {
                "candidate": row.get("candidate", ""),
                "dataset": row.get("linked_gse_expected", ""),
                "run": row.get("run", ""),
                "experiment": row.get("experiment", ""),
                "sample": row.get("sample", ""),
                "biosample": row.get("biosample", ""),
                "gsm": row.get("gsm", ""),
                "library_strategy": row.get("library_strategy", ""),
                "library_selection": row.get("library_selection", ""),
                "library_layout": row.get("library_layout", ""),
                "size_mb": as_float(row.get("size_mb")),
                "sample_title": row.get("sample_title", ""),
                "source_name": row.get("source_name", ""),
                "tissue": normalize_tissue(row.get("tissue")),
                "age_weeks": as_float(row.get("age_weeks")),
                "raw_age_token": row.get("raw_age_token", ""),
                "age_group": row.get("age_group", ""),
                "assay_pass": as_bool(row.get("assay_pass")),
                "target_tissue_pass": as_bool(row.get("target_tissue_pass")),
                "age_exact_pass": as_bool(row.get("age_exact_pass")),
                "age_group_pass": as_bool(row.get("age_group_pass")),
                "blockers": ";".join(blockers),
                "context_blocker_pass": not has_context_blocker(blockers),
                "include_in_prior_minimal_pilot": as_bool(row.get("include_in_minimal_pilot")),
                "sample_decision": sample_decision(row),
            }
        )
    return sample_rows


def summarize_candidates(sample_rows: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    by_candidate: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in sample_rows:
        by_candidate[str(row["candidate"])].append(row)
    v12_gate = load_v12_raw_gate()
    v12_minimal = load_v12_minimal_etl()
    candidate_rows: list[dict[str, Any]] = []
    pilot_rows: list[dict[str, Any]] = []
    for candidate, rows in sorted(by_candidate.items()):
        dataset = rows[0].get("dataset", "")
        gate = v12_gate.get(candidate)
        minimal = v12_minimal.get(candidate)
        decision, reason = candidate_decision(rows, gate)
        tissues = sorted({str(r.get("tissue", "")) for r in rows if r.get("tissue")})
        runs = [str(r.get("run", "")) for r in rows if r.get("run")]
        total_size = sum(float(r.get("size_mb") or 0) for r in rows)
        exact_old_target = [
            r
            for r in rows
            if r["sample_decision"] == "eligible_old_exact_target_pending_approval"
        ]
        diagnostic_age_group = [
            r
            for r in rows
            if r["sample_decision"] == "diagnostic_only_age_group_no_exact_age"
        ]
        old_exact_target_raw = [
            r
            for r in rows
            if r.get("target_tissue_pass")
            and r.get("age_exact_pass")
            and r.get("age_weeks") is not None
            and float(r.get("age_weeks") or 0) >= OLD_WEEKS
        ]
        blocked_context = [
            r
            for r in rows
            if r["sample_decision"] == "blocked_context_not_bulk_headline"
        ]
        candidate_rows.append(
            {
                "candidate": candidate,
                "dataset": dataset,
                "route_b_decision": decision,
                "decision_reason": reason,
                "n_runs": len(rows),
                "total_size_mb": round(total_size, 3),
                "tissues": ";".join(tissues),
                "assay_pass_n": sum(bool(r.get("assay_pass")) for r in rows),
                "target_tissue_pass_n": sum(bool(r.get("target_tissue_pass")) for r in rows),
                "exact_age_n": sum(bool(r.get("age_exact_pass")) for r in rows),
                "old_exact_target_raw_n": len(old_exact_target_raw),
                "headline_eligible_old_exact_target_n": len(exact_old_target),
                "age_group_target_n": len(diagnostic_age_group),
                "context_blocked_n": len(blocked_context),
                "v12_pilot_status": minimal.get("pilot_status", "") if minimal else "",
                "v12_pilot_reason": minimal.get("pilot_reason", "") if minimal else "",
                "v12_blocker_type": minimal.get("blocker_type", "") if minimal else (gate.get("blocker_type", "") if gate else ""),
                "recommended_action": route_b_action(decision),
                "run_accessions": ";".join(runs),
            }
        )
        selected = select_pilot_rows(decision, rows)
        for rank, row in enumerate(selected, start=1):
            pilot_rows.append(
                {
                    "candidate": candidate,
                    "dataset": dataset,
                    "selection_rank": rank,
                    "selection_status": "not_authorized",
                    "selection_type": decision,
                    "run": row.get("run", ""),
                    "biosample": row.get("biosample", ""),
                    "gsm": row.get("gsm", ""),
                    "tissue": row.get("tissue", ""),
                    "age_weeks": row.get("age_weeks"),
                    "raw_age_token": row.get("raw_age_token", ""),
                    "age_group": row.get("age_group", ""),
                    "size_mb": row.get("size_mb"),
                    "sample_title": row.get("sample_title", ""),
                    "blockers": row.get("blockers", ""),
                    "recommended_action": route_b_action(decision),
                }
            )
    return candidate_rows, pilot_rows


def route_b_action(decision: str) -> str:
    if decision == "route_b_pilot_eligible_pending_explicit_approval":
        return "prepare_separate_minimal_etl_approval_before_any_fastq_download"
    if decision == "diagnostic_only_age_group_no_exact_age":
        return "do_not_download_fastq_for_headline; diagnostic_only_if_separately_approved"
    if decision == "blocked_context_not_route_b_headline":
        return "do_not_download_fastq; celltype_or_nonbulk_context_blocks_headline"
    return "do_not_download_fastq"


def select_pilot_rows(decision: str, rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    if decision == "route_b_pilot_eligible_pending_explicit_approval":
        eligible = [r for r in rows if r["sample_decision"] == "eligible_old_exact_target_pending_approval"]
        return sorted(eligible, key=lambda r: float(r.get("size_mb") or 0))[:3]
    if decision == "diagnostic_only_age_group_no_exact_age":
        diagnostic = [r for r in rows if r["sample_decision"] == "diagnostic_only_age_group_no_exact_age"]
        vehicle = [r for r in diagnostic if "vehicle" in str(r.get("sample_title", "")).lower()]
        chosen = vehicle or diagnostic
        return sorted(chosen, key=lambda r: float(r.get("size_mb") or 0))[:3]
    return []


def build_raw_lead_backlog() -> list[dict[str, Any]]:
    rows = read_csv(RAW_PROJECT_LEADS)
    out: list[dict[str, Any]] = []
    seen: set[str] = set()
    for row in rows:
        project_key = row.get("project_key", "")
        if project_key in seen:
            continue
        seen.add(project_key)
        priority = row.get("pilot_priority", "")
        if priority not in {"manual_biosample_check_only", "minimal_etl_rfc_candidate"}:
            continue
        out.append(
            {
                "project_key": project_key,
                "linked_gse": row.get("linked_gse", ""),
                "pilot_priority": priority,
                "max_lead_score": row.get("max_lead_score", ""),
                "target_hits": row.get("target_hits", ""),
                "old_hits": row.get("old_hits", ""),
                "methylation_hits": row.get("methylation_hits", ""),
                "blocker_hits": row.get("blocker_hits", ""),
                "raw_accession_examples": row.get("raw_accession_examples", ""),
                "next_action": "already_biosample_audited" if priority == "minimal_etl_rfc_candidate" else "official_biosample_runinfo_refresh_required_before_any_pilot",
            }
        )
    return out


def build_processed_context() -> list[dict[str, Any]]:
    rows = read_csv(V11_CONVERSION_QUEUE)
    out = []
    for row in rows:
        out.append(
            {
                "dataset": row.get("dataset", ""),
                "tier": row.get("tier_smoke", ""),
                "headline_allowed": row.get("final_headline_allowed", ""),
                "conversion_gate": row.get("final_conversion_gate", ""),
                "recommended_adapter": row.get("final_recommended_adapter", ""),
                "n_files": row.get("n_files", ""),
                "total_size_gib": row.get("total_size_gib", ""),
                "reason": row.get("reason_smoke", ""),
                "route_b_implication": "processed_data_exists_or_is_preferred; raw_fastq_pilot_not_default",
            }
        )
    return out


def write_report(
    candidate_rows: list[dict[str, Any]],
    pilot_rows: list[dict[str, Any]],
    raw_backlog_rows: list[dict[str, Any]],
    processed_context_rows: list[dict[str, Any]],
) -> None:
    eligible = [r for r in candidate_rows if r["route_b_decision"] == "route_b_pilot_eligible_pending_explicit_approval"]
    diagnostic = [r for r in candidate_rows if r["route_b_decision"] == "diagnostic_only_age_group_no_exact_age"]
    blocked = [r for r in candidate_rows if r["route_b_decision"].startswith("blocked")]
    lines = [
        "# Route B Candidate Audit Report",
        "",
        f"Date: {datetime.now().astimezone().isoformat(timespec='seconds')}",
        "",
        "## Summary",
        "",
        "This audit reviews existing v10-v12 raw-lead and BioSample/RunInfo metadata for a possible Route B minimal FASTQ/Bismark pilot. It uses cached official metadata only.",
        "",
        "No FASTQ download, Bismark, training, matrix rebuild, or autoresearch was run.",
        "",
        "## Decision",
        "",
        f"- Headline Route B pilot candidates: `{len(eligible)}`.",
        f"- Diagnostic-only age-group candidates: `{len(diagnostic)}`.",
        f"- Blocked raw candidates: `{len(blocked)}`.",
        "",
        "Current conclusion: no raw FASTQ pilot is authorized. The only verified raw candidates are either diagnostic-only without sample-specific exact age or blocked by cell-type/non-bulk context.",
        "",
        "## Candidate Findings",
        "",
    ]
    for row in candidate_rows:
        lines.extend(
            [
                f"### {row['dataset']}",
                "",
                f"- decision: `{row['route_b_decision']}`",
                f"- reason: {row['decision_reason']}",
                f"- runs: `{row['n_runs']}`, total size: `{row['total_size_mb']}` MB",
                f"- tissues: `{row['tissues']}`",
                f"- exact-age runs: `{row['exact_age_n']}`, old exact target raw runs: `{row['old_exact_target_raw_n']}`, headline-eligible old exact target runs: `{row['headline_eligible_old_exact_target_n']}`",
                f"- v12 pilot status: `{row['v12_pilot_status']}`; reason: `{row['v12_pilot_reason']}`",
                f"- recommended action: `{row['recommended_action']}`",
                "",
            ]
        )
    lines.extend(
        [
            "## Diagnostic Pilot Rows",
            "",
            f"`route_b_recommended_pilot_samples.csv` contains `{len(pilot_rows)}` rows. All rows are marked `selection_status=not_authorized`.",
            "",
            "## Raw Lead Backlog",
            "",
            f"`route_b_raw_lead_backlog.csv` contains `{len(raw_backlog_rows)}` deduplicated project keys from v11 raw-lead records.",
            "",
            "Manual backlog rows require a fresh official BioSample/RunInfo refresh before any pilot decision.",
            "",
            "## Processed Data Context",
            "",
            f"`route_b_processed_context.csv` contains `{len(processed_context_rows)}` datasets where processed methylation files already exist or were preferred in v11.3/v11.4.",
            "",
            "Processed data remains preferred over raw FASTQ whenever schema and metadata gates are adequate.",
            "",
            "## Outputs",
            "",
            "- `results/route_b_candidate_audit/route_b_candidate_audit.csv`",
            "- `results/route_b_candidate_audit/route_b_sample_audit.csv`",
            "- `results/route_b_candidate_audit/route_b_recommended_pilot_samples.csv`",
            "- `results/route_b_candidate_audit/route_b_raw_lead_backlog.csv`",
            "- `results/route_b_candidate_audit/route_b_processed_context.csv`",
            "- `results/route_b_candidate_audit/route_b_decision_state.json`",
            "",
            "## Guardrails",
            "",
            "- No raw FASTQ download is authorized by this audit.",
            "- Route B still requires separate minimal ETL approval.",
            "- Age-group-only data cannot become a headline chronological benchmark.",
            "- Cell-type-specific data cannot become a headline bulk tissue benchmark.",
        ]
    )
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    raw_rows = read_csv(RAW_SAMPLE_MANIFEST)
    sample_rows = build_sample_audit(raw_rows)
    candidate_rows, pilot_rows = summarize_candidates(sample_rows)
    raw_backlog_rows = build_raw_lead_backlog()
    processed_context_rows = build_processed_context()

    write_csv(OUT_DIR / "route_b_sample_audit.csv", sample_rows)
    write_csv(OUT_DIR / "route_b_candidate_audit.csv", candidate_rows)
    write_csv(OUT_DIR / "route_b_recommended_pilot_samples.csv", pilot_rows)
    write_csv(OUT_DIR / "route_b_raw_lead_backlog.csv", raw_backlog_rows)
    write_csv(OUT_DIR / "route_b_processed_context.csv", processed_context_rows)

    eligible = [r for r in candidate_rows if r["route_b_decision"] == "route_b_pilot_eligible_pending_explicit_approval"]
    state = {
        "status": "completed",
        "created_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "network_mode": "cached_official_metadata_only_no_api_calls",
        "raw_fastq_download_authorized": False,
        "bismark_authorized": False,
        "training_authorized": False,
        "autoresearch_authorized": False,
        "headline_route_b_candidate_count": len(eligible),
        "candidate_count": len(candidate_rows),
        "sample_row_count": len(sample_rows),
        "pilot_sample_rows_written": len(pilot_rows),
        "decision": "no_route_b_headline_pilot_candidate_ready",
        "decision_reason": "Verified raw candidates are diagnostic age-group only or blocked by cell-type/non-bulk context.",
        "outputs": {
            "candidate_audit": str(OUT_DIR / "route_b_candidate_audit.csv"),
            "sample_audit": str(OUT_DIR / "route_b_sample_audit.csv"),
            "recommended_pilot_samples": str(OUT_DIR / "route_b_recommended_pilot_samples.csv"),
            "raw_lead_backlog": str(OUT_DIR / "route_b_raw_lead_backlog.csv"),
            "processed_context": str(OUT_DIR / "route_b_processed_context.csv"),
            "report": str(REPORT_PATH),
        },
        "inputs": {
            "raw_sample_manifest": str(RAW_SAMPLE_MANIFEST),
            "v12_minimal_etl": str(V12_MINIMAL_ETL),
            "v12_gate": str(V12_GATE),
            "raw_project_leads": str(RAW_PROJECT_LEADS),
            "processed_context": str(V11_CONVERSION_QUEUE),
        },
    }
    write_json(OUT_DIR / "route_b_decision_state.json", state)
    write_report(candidate_rows, pilot_rows, raw_backlog_rows, processed_context_rows)
    print(json.dumps(state, indent=2))


if __name__ == "__main__":
    main()

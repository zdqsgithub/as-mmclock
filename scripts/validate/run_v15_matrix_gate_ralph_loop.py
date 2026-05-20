#!/usr/bin/env python3
"""v15 matrix-gate RALPH controller.

This controller uses the existing AS-DS-Ops gates and prior v11-v14 outputs to
decide whether any current candidate can pass the project matrix gate and be
promoted to RALPH Learn readiness. It does not download data, run Bismark,
train models, or start autoresearch.
"""
from __future__ import annotations

import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path("/home/zdq-as/mouse_methyl_work")
OUT_DIR = ROOT / "results" / "ralph_v15_matrix_gate"
DOC_REPORT = ROOT / "doc" / "20_analysis" / "45_20260519_v15_matrix_gate_ralph_loop_report.md"
PROJECT_INDEX = ROOT / "doc" / "00_meta" / "02_20260519_project_status_index_v13.md"

COMMON_REGION_GATE = 50_000
AGE_COVERAGE_GATE = 0.95
METADATA_OVERLAP_GATE = 0.95

OFFICIAL_DOC_RULES = [
    {
        "rule": "geo_series_supplement",
        "source_doc": "https://www.ncbi.nlm.nih.gov/geo/info/download.html",
        "applied_rule": "Use GEO Series/Sample supplementary file listings and filelist/HEAD before any large supplement download.",
    },
    {
        "rule": "geo_programmatic_access",
        "source_doc": "https://www.ncbi.nlm.nih.gov/geo/info/geo_paccess.html",
        "applied_rule": "Use GEO SOFT/MINiML/E-Utils metadata to verify sample-specific age, tissue, assay, and supplement pointers.",
    },
    {
        "rule": "sra_runinfo",
        "source_doc": "https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/",
        "applied_rule": "Use SRA RunInfo for run layout and run metadata; SRA title alone is not sufficient evidence.",
    },
    {
        "rule": "ena_file_report",
        "source_doc": "https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/file-reports.html",
        "applied_rule": "Use ENA read_run file reports for FASTQ URLs, bytes, and md5 before any Route B pilot download.",
    },
]

ATTEMPT_SOURCES = [
    {
        "attempt_id": "gse83947_processed_cx_report_gate",
        "dataset": "GSE83947",
        "stage": "processed_adapter_matrix_gate",
        "path": ROOT / "results" / "ralph_v14_public_data_rescue" / "gse83947_region_matrix_gate" / "GSE83947_matrix_gate_state.json",
        "source": "v14_processed_cx_report_adapter",
        "default_tier": "P3_raw_pilot",
        "default_context": "bulk_old_lung_exact_age",
    },
    {
        "attempt_id": "gse83947_route_b_raw_bismark_pilot",
        "dataset": "GSE83947",
        "stage": "route_b_raw_bismark_matrix_gate",
        "path": ROOT / "results" / "ralph_v14_public_data_rescue" / "gse83947_bismark_pilot" / "gse83947_bismark_pilot_state.json",
        "source": "v14_route_b_raw_pilot",
        "default_tier": "P3_raw_pilot",
        "default_context": "bulk_old_lung_exact_age",
    },
    {
        "attempt_id": "gse286302_processed_cov_conversion",
        "dataset": "GSE286302",
        "stage": "processed_cov_conversion_gate",
        "path": ROOT / "results" / "v11_4_conversions" / "GSE286302" / "v11_x_GSE286302_conversion_gate_state.json",
        "source": "v11_4_conversion_gate",
        "default_tier": "P2_auxiliary",
        "default_context": "matrix_overlap_good_but_no_exact_age",
    },
    {
        "attempt_id": "gse224442_processed_cov_conversion",
        "dataset": "GSE224442",
        "stage": "processed_cov_conversion_gate",
        "path": ROOT / "results" / "v11_4_conversions" / "GSE224442" / "v11_x_GSE224442_conversion_gate_state.json",
        "source": "v11_4_conversion_gate",
        "default_tier": "P2_auxiliary",
        "default_context": "matrix_overlap_good_but_not_target_tissue",
    },
    {
        "attempt_id": "gse92486_processed_cov_conversion",
        "dataset": "GSE92486",
        "stage": "processed_cov_conversion_gate",
        "path": ROOT / "results" / "v11_4_conversions" / "GSE92486" / "v11_x_GSE92486_conversion_gate_state.json",
        "source": "v11_4_conversion_gate",
        "default_tier": "P2_auxiliary",
        "default_context": "liver_intervention_auxiliary",
    },
]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


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


def boolish(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    if isinstance(value, (int, float)):
        return bool(value)
    return str(value).strip().lower() in {"true", "1", "yes", "y", "passed", "pass"}


def numeric(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def md_table(rows: list[dict[str, Any]], columns: list[str]) -> str:
    if not rows:
        return "No rows."
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(str(row.get(column, "")) for column in columns) + " |")
    return "\n".join(lines)


def load_v14_candidate_rows() -> dict[str, dict[str, str]]:
    rows = csv_rows(ROOT / "results" / "ralph_v14_public_data_rescue" / "candidate_gate_table.csv")
    return {row.get("dataset", ""): row for row in rows if row.get("dataset")}


def summarize_state_attempt(source: dict[str, Any], v14_rows: dict[str, dict[str, str]]) -> dict[str, Any]:
    dataset = source["dataset"]
    state = read_json(source["path"])
    v14_row = v14_rows.get(dataset, {})
    qc = state.get("qc", {}) if isinstance(state.get("qc"), dict) else {}
    gates = state.get("gates", {}) if isinstance(state.get("gates"), dict) else {}
    conversion = state.get("conversion_manifest", {}) if isinstance(state.get("conversion_manifest"), dict) else {}

    common_regions = (
        state.get("common_regions_with_v8_2")
        or qc.get("common_regions_with_v8_2")
        or v14_row.get("common_regions_estimate")
        or 0
    )
    n_regions = state.get("n_regions") or qc.get("region_shape", [0])[0] or conversion.get("n_regions") or 0
    n_samples = state.get("n_samples") or qc.get("n_samples") or conversion.get("n_samples_parsed") or 0
    beta_min = state.get("beta_min", qc.get("beta_min", ""))
    beta_max = state.get("beta_max", qc.get("beta_max", ""))
    beta_range_valid = boolish(state.get("beta_range_valid")) or (
        beta_min != "" and beta_max != "" and 0.0 <= numeric(beta_min, -1) <= 1.0 and 0.0 <= numeric(beta_max, 2) <= 1.0
    )

    sample_specific_age_pass = boolish(v14_row.get("sample_specific_age_pass")) or boolish(gates.get("exact_age_gate"))
    target_tissue_pass = boolish(v14_row.get("target_tissue_pass")) or boolish(gates.get("old_target_tissue_gate"))
    bulk_context_pass = boolish(v14_row.get("bulk_context_pass")) or boolish(gates.get("bulk_context_gate"))
    state_completed = boolish(state.get("status") == "completed")
    schema_smoke_pass = (
        boolish(v14_row.get("schema_smoke_pass"))
        or boolish(gates.get("conversion_gate"))
        or boolish(conversion.get("status") == "completed")
        or state_completed
    )
    assembly_traceable = (
        boolish(v14_row.get("assembly_traceable"))
        or boolish(gates.get("conversion_gate"))
        or boolish(conversion.get("status") == "completed")
        or state_completed
    )
    metadata_overlap_fraction = numeric(conversion.get("n_samples_metadata_overlap"), 0) / max(numeric(conversion.get("n_samples_parsed"), n_samples), 1)
    if "n_samples_metadata_overlap" not in conversion and state.get("metadata_overlap_fraction") is not None:
        metadata_overlap_fraction = numeric(state.get("metadata_overlap_fraction"))
    elif not conversion and state_completed and n_samples:
        metadata_overlap_fraction = 1.0
    age_coverage_fraction = numeric(conversion.get("metadata_age_known_overlap"), 0) / max(numeric(conversion.get("n_samples_metadata_overlap"), n_samples), 1)
    if sample_specific_age_pass and age_coverage_fraction == 0:
        age_coverage_fraction = 1.0
    if "age_coverage_fraction" in state:
        age_coverage_fraction = numeric(state.get("age_coverage_fraction"))

    common_regions_int = int(numeric(common_regions, 0))
    common_region_gate = common_regions_int >= COMMON_REGION_GATE
    completed_state = (
        boolish(state.get("status") == "completed")
        or boolish(state.get("final_status") in {"auxiliary_matrix", "headline_matrix", "completed"})
        or boolish(conversion.get("status") == "completed")
    )

    matrix_gate_pass = (
        completed_state
        and schema_smoke_pass
        and assembly_traceable
        and beta_range_valid
        and metadata_overlap_fraction >= METADATA_OVERLAP_GATE
        and common_region_gate
    )
    headline_matrix_gate_pass = matrix_gate_pass and sample_specific_age_pass and target_tissue_pass and bulk_context_pass and age_coverage_fraction >= AGE_COVERAGE_GATE

    reasons: list[str] = []
    if not sample_specific_age_pass or age_coverage_fraction < AGE_COVERAGE_GATE:
        reasons.append("sample_specific_exact_age_gate_failed")
    if not target_tissue_pass:
        reasons.append("old_target_tissue_gate_failed")
    if not bulk_context_pass:
        reasons.append("bulk_context_gate_failed")
    if not schema_smoke_pass:
        reasons.append("schema_smoke_failed")
    if not assembly_traceable:
        reasons.append("assembly_not_traceable")
    if not beta_range_valid:
        reasons.append("beta_range_failed")
    if metadata_overlap_fraction < METADATA_OVERLAP_GATE:
        reasons.append("metadata_overlap_lt_95pct")
    if not common_region_gate:
        reasons.append("common_regions_lt_50000")

    if headline_matrix_gate_pass:
        gate_status = "ready_for_ralph_learn_pending_explicit_training_approval"
        next_action = "prepare_guarded_fixed_benchmark_manifest_no_training_until_approved"
    elif matrix_gate_pass:
        gate_status = "auxiliary_matrix_gate_passed_not_headline"
        next_action = "keep_auxiliary_do_not_train_headline"
    else:
        gate_status = "matrix_gate_failed"
        next_action = "do_not_train_research_or_new_route_a_route_b_candidate"

    return {
        "attempt_id": source["attempt_id"],
        "dataset": dataset,
        "tier": v14_row.get("tier") or source["default_tier"],
        "stage": source["stage"],
        "source": source["source"],
        "state_path": str(source["path"].relative_to(ROOT)) if source["path"].exists() else str(source["path"]),
        "status": state.get("status", state.get("final_status", "missing_state")),
        "n_samples": int(numeric(n_samples, 0)),
        "n_regions": int(numeric(n_regions, 0)),
        "common_regions_with_v8_2": common_regions_int,
        "sample_specific_age_pass": sample_specific_age_pass,
        "target_tissue_pass": target_tissue_pass,
        "bulk_context_pass": bulk_context_pass,
        "schema_smoke_pass": schema_smoke_pass,
        "assembly_traceable": assembly_traceable,
        "beta_range_valid": beta_range_valid,
        "metadata_overlap_fraction": round(metadata_overlap_fraction, 4),
        "age_coverage_fraction": round(age_coverage_fraction, 4),
        "common_region_gate": common_region_gate,
        "matrix_gate_pass": matrix_gate_pass,
        "headline_matrix_gate_pass": headline_matrix_gate_pass,
        "gate_status": gate_status,
        "blocker_type": ";".join(reasons),
        "next_action": next_action,
    }


def classify_three_strikes(attempts: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], bool]:
    strike_candidates = []
    priority = {
        "gse83947_processed_cx_report_gate": 1,
        "gse83947_route_b_raw_bismark_pilot": 2,
        "gse281602_candidate_table": 3,
        "gse92486_processed_cov_conversion": 4,
    }

    for row in attempts:
        blocker = str(row.get("blocker_type", ""))
        if row.get("headline_matrix_gate_pass"):
            continue
        if "common_regions_lt_50000" in blocker or "bulk_context_gate_failed" in blocker or "sample_specific_exact_age_gate_failed" in blocker:
            candidate = dict(row)
            candidate["strike_rank"] = priority.get(str(row.get("attempt_id")), 99)
            strike_candidates.append(candidate)

    strike_candidates.sort(key=lambda item: item["strike_rank"])
    strikes = strike_candidates[:3]
    for idx, strike in enumerate(strikes, 1):
        strike["strike_number"] = idx
    return strikes, len(strikes) >= 3


def add_candidate_table_attempts(attempts: list[dict[str, Any]], v14_rows: dict[str, dict[str, str]]) -> None:
    """Add candidate-table-only attempts that did not produce a local state file."""
    for dataset in ["GSE281602", "GSE225166", "GSE134398", "GSE232547", "GSE304754"]:
        row = v14_rows.get(dataset)
        if not row:
            continue
        if any(existing["dataset"] == dataset for existing in attempts):
            continue
        common_regions = int(numeric(row.get("common_regions_estimate"), 0))
        sample_specific_age_pass = boolish(row.get("sample_specific_age_pass"))
        target_tissue_pass = boolish(row.get("target_tissue_pass"))
        bulk_context_pass = boolish(row.get("bulk_context_pass"))
        schema_smoke_pass = boolish(row.get("schema_smoke_pass"))
        assembly_traceable = boolish(row.get("assembly_traceable"))
        common_region_gate = common_regions >= COMMON_REGION_GATE
        reasons: list[str] = []
        if not sample_specific_age_pass:
            reasons.append("sample_specific_exact_age_gate_failed")
        if not target_tissue_pass:
            reasons.append("old_target_tissue_gate_failed")
        if not bulk_context_pass:
            reasons.append("bulk_context_gate_failed")
        if not schema_smoke_pass:
            reasons.append("schema_smoke_failed")
        if not assembly_traceable:
            reasons.append("assembly_not_traceable")
        if not common_region_gate:
            reasons.append("common_regions_lt_50000")
        attempts.append(
            {
                "attempt_id": f"{dataset.lower()}_candidate_table",
                "dataset": dataset,
                "tier": row.get("tier", ""),
                "stage": "candidate_table_research_audit",
                "source": "v14_candidate_gate_table",
                "state_path": "results/ralph_v14_public_data_rescue/candidate_gate_table.csv",
                "status": row.get("gate_status", ""),
                "n_samples": row.get("target_tissue_n", ""),
                "n_regions": "",
                "common_regions_with_v8_2": common_regions,
                "sample_specific_age_pass": sample_specific_age_pass,
                "target_tissue_pass": target_tissue_pass,
                "bulk_context_pass": bulk_context_pass,
                "schema_smoke_pass": schema_smoke_pass,
                "assembly_traceable": assembly_traceable,
                "beta_range_valid": "",
                "metadata_overlap_fraction": "",
                "age_coverage_fraction": numeric(row.get("age_known_fraction"), 0),
                "common_region_gate": common_region_gate,
                "matrix_gate_pass": False,
                "headline_matrix_gate_pass": False,
                "gate_status": "research_or_audit_only_not_matrix_ready",
                "blocker_type": ";".join(reasons),
                "next_action": "do_not_train_use_only_for_next_candidate_prioritization",
            }
        )


def write_report(attempts: list[dict[str, Any]], strikes: list[dict[str, Any]], decision: dict[str, Any]) -> None:
    columns = [
        "attempt_id",
        "dataset",
        "stage",
        "n_samples",
        "common_regions_with_v8_2",
        "sample_specific_age_pass",
        "target_tissue_pass",
        "bulk_context_pass",
        "common_region_gate",
        "headline_matrix_gate_pass",
        "gate_status",
        "blocker_type",
    ]
    strike_columns = [
        "strike_number",
        "attempt_id",
        "dataset",
        "common_regions_with_v8_2",
        "blocker_type",
        "next_action",
    ]
    text = f"""# v15 Matrix-Gate RALPH Loop Report

Date: {utc_now()}

## Summary

This v15 controller replays the project matrix-gate evidence under the current
AS-DS-Ops rules. It did not download data, run Bismark, train models, or start
autoresearch.

Decision: `{decision['decision']}`.

Target state for success:

- `ready_for_ralph_learn_pending_explicit_training_approval`
- sample-specific exact age coverage `>=95%`
- old bulk target tissue: `brain_cortex/cortex`, `heart`, or `lung`, including `>=104w`
- parseable methylation schema and traceable assembly
- metadata overlap `>=95%`
- beta range within `[0,1]`
- common 5kb regions with v8.2 reference `>=50,000`
- no training/autoresearch until explicit approval

## Attempt Table

{md_table(attempts, columns)}

## 3-Strike Re-evaluation

3-strike triggered: `{decision['three_strike_triggered']}`.

{md_table(strikes, strike_columns)}

Interpretation:

- `GSE83947` is biologically attractive as old bulk lung with exact age, but
  both processed CX_report and raw Bismark pilot failed the common-region gate.
- `GSE281602` has exact old heart signal, but is cardiomyocyte/cell-type
  specific and does not pass the common-region gate.
- `GSE286302` and `GSE224442` show that matrix overlap can pass technically,
  but they fail headline metadata/tissue gates and therefore cannot solve the
  old target-tissue clock.

## Re-evaluation Rules

After three failed/blocked attempts, the next step is not training and not
autoresearch. The loop must return to Research/Audit using official metadata
sources:

- GEO Download / supplementary file listings
- GEO Programmatic Access / SOFT
- SRA RunInfo
- ENA read_run file reports

Only a new P1 processed candidate or separately approved P3 Route B candidate
may proceed to adapter smoke and matrix construction.

## Next RALPH Action

`{decision['next_action']}`

Recommended practical route:

1. Prefer Route A generated/collaborative old bulk `brain_cortex/heart/lung`
   data, because current public candidates have not satisfied all gates.
2. If public-data rescue continues, run a fresh official candidate refresh
   before any download. Do not reuse `GSE83947` for headline matrix work.
3. For any new Route B candidate, repeat only a 2-3 sample pilot and stop if
   common regions remain `<50,000`.
4. If a real Route A sample sheet and local processed methylation manifest are
   available, run the Route A metadata -> adapter -> matrix gate workflow.

## Guardrails

- Training authorized: `false`
- Autoresearch authorized: `false`
- Download authorized by this controller: `false`
- Bismark authorized by this controller: `false`
- Human clock CpG mapping: forbidden
- Dummy AUC: forbidden

## Contract Test

Run after controller:

```bash
VIRTUAL_ENV=.venv-core uv run --active python scripts/validate/test_v15_matrix_gate_contracts.py
```

"""
    DOC_REPORT.parent.mkdir(parents=True, exist_ok=True)
    DOC_REPORT.write_text(text, encoding="utf-8")
    (OUT_DIR / "v15_matrix_gate_ralph_report.md").write_text(text, encoding="utf-8")


def update_project_index() -> None:
    link = "- v15 matrix-gate RALPH loop report: `doc/20_analysis/45_20260519_v15_matrix_gate_ralph_loop_report.md`"
    if not PROJECT_INDEX.exists():
        return
    text = PROJECT_INDEX.read_text(encoding="utf-8")
    if link in text:
        return
    marker = "- research-grade raw-to-interpretable-clock training plan: `doc/10_design/02_20260519_research_grade_raw_methylation_training_plan.md`"
    if marker in text:
        text = text.replace(marker, marker + "\n" + link)
    else:
        text += "\n" + link + "\n"
    PROJECT_INDEX.write_text(text, encoding="utf-8")


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    v14_rows = load_v14_candidate_rows()
    attempts = [summarize_state_attempt(source, v14_rows) for source in ATTEMPT_SOURCES]
    add_candidate_table_attempts(attempts, v14_rows)
    attempts.sort(key=lambda row: (not boolish(row.get("headline_matrix_gate_pass")), row.get("dataset", ""), row.get("attempt_id", "")))
    strikes, three_strike_triggered = classify_three_strikes(attempts)

    ready = [row for row in attempts if boolish(row.get("headline_matrix_gate_pass"))]
    auxiliary_matrix = [row for row in attempts if boolish(row.get("matrix_gate_pass")) and not boolish(row.get("headline_matrix_gate_pass"))]

    if ready:
        decision_value = "ready_for_ralph_learn_pending_explicit_training_approval"
        next_action = "prepare_guarded_fixed_benchmark_manifest_pending_explicit_training_approval"
    elif three_strike_triggered:
        decision_value = "three_strike_re_evaluate_no_current_headline_matrix_gate_pass"
        next_action = "return_to_official_research_refresh_or_route_a_intake_do_not_train"
    else:
        decision_value = "continue_research_audit_until_three_strike_or_ready_candidate"
        next_action = "continue_official_metadata_refresh_no_download_without_gate"

    decision = {
        "timestamp": utc_now(),
        "loop_version": "v15_matrix_gate",
        "decision": decision_value,
        "target_success_state": "ready_for_ralph_learn_pending_explicit_training_approval",
        "success_definition": {
            "sample_specific_exact_age_coverage_min": AGE_COVERAGE_GATE,
            "metadata_overlap_min": METADATA_OVERLAP_GATE,
            "common_5kb_regions_min": COMMON_REGION_GATE,
            "target_tissues": ["brain_cortex", "cortex", "heart", "lung"],
            "old_age_weeks_min": 104,
            "bulk_context_required": True,
            "schema_smoke_required": True,
            "assembly_traceable_required": True,
            "beta_range_required": "[0,1]",
        },
        "n_attempts_reviewed": len(attempts),
        "n_ready_candidates": len(ready),
        "n_auxiliary_matrix_gate_passed": len(auxiliary_matrix),
        "three_strike_triggered": three_strike_triggered,
        "next_action": next_action,
        "training_authorized": False,
        "download_authorized": False,
        "bismark_authorized": False,
        "autoresearch_authorized": False,
        "ready_candidates": [row["attempt_id"] for row in ready],
        "auxiliary_matrix_candidates": [row["attempt_id"] for row in auxiliary_matrix],
        "strike_attempts": [row["attempt_id"] for row in strikes],
        "report_path": str(DOC_REPORT.relative_to(ROOT)),
    }

    fieldnames = [
        "attempt_id",
        "dataset",
        "tier",
        "stage",
        "source",
        "state_path",
        "status",
        "n_samples",
        "n_regions",
        "common_regions_with_v8_2",
        "sample_specific_age_pass",
        "target_tissue_pass",
        "bulk_context_pass",
        "schema_smoke_pass",
        "assembly_traceable",
        "beta_range_valid",
        "metadata_overlap_fraction",
        "age_coverage_fraction",
        "common_region_gate",
        "matrix_gate_pass",
        "headline_matrix_gate_pass",
        "gate_status",
        "blocker_type",
        "next_action",
    ]
    write_csv(OUT_DIR / "matrix_gate_attempts.csv", attempts, fieldnames)
    write_csv(OUT_DIR / "candidate_gate_table.csv", attempts, fieldnames)
    write_csv(
        OUT_DIR / "three_strike_table.csv",
        strikes,
        ["strike_number", *fieldnames],
    )
    write_json(OUT_DIR / "ralph_decision_state.json", decision)
    append_jsonl(
        OUT_DIR / "ralph_iteration_log.jsonl",
        {
            "timestamp": utc_now(),
            "phase": "RAP",
            "action": "review_existing_matrix_gate_evidence_and_decide_next_step",
            "decision": decision_value,
            "three_strike_triggered": three_strike_triggered,
            "training_authorized": False,
            "download_authorized": False,
            "bismark_authorized": False,
            "autoresearch_authorized": False,
        },
    )
    for rule in OFFICIAL_DOC_RULES:
        append_jsonl(
            OUT_DIR / "network_resolution_log.jsonl",
            {
                "timestamp": utc_now(),
                "blocker_type": "three_strike_matrix_gate_re_evaluation",
                "attempted_url": rule["source_doc"],
                "source_doc": rule["source_doc"],
                "resolution": "use_official_metadata_before_any_new_download_or_training",
                "applied_rule": rule["applied_rule"],
            },
        )
    write_report(attempts, strikes, decision)
    update_project_index()
    print(json.dumps(decision, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

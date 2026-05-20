#!/usr/bin/env python3
"""Generate the concrete Route A data-generation RFC package.

This is a planning and governance script only. It does not train models,
download data, run Bismark, rebuild matrices, or authorize wet-lab work.
"""
from __future__ import annotations

import csv
import json
from datetime import datetime
from pathlib import Path
from typing import Any

ROOT = Path("/home/zdq-as/mouse_methyl_work")
OUT_DIR = ROOT / "results" / "route_a_data_generation_rfc"
REPORT_PATH = ROOT / "doc" / "20_analysis" / "36_20260519_route_a_old_tissue_data_generation_rfc.md"
CHECKLIST_PATH = ROOT / "doc" / "30_protocols" / "06_20260519_route_a_submission_checklist.md"
INDEX_PATH = ROOT / "doc" / "00_meta" / "02_20260519_project_status_index_v13.md"

DESIGN_CSV = OUT_DIR / "route_a_cohort_design.csv"
SCHEMA_CSV = OUT_DIR / "route_a_required_metadata_schema.csv"
GATES_JSON = OUT_DIR / "route_a_qc_benchmark_gates.json"
STATE_JSON = OUT_DIR / "route_a_decision_state.json"

TARGET_TISSUES = ["brain_cortex", "heart", "lung"]
AGE_STRATA = [
    {
        "age_stratum": "young",
        "target_age_weeks": 12,
        "acceptable_window_weeks": "8-16",
        "purpose": "anchor young adult methylation baseline",
    },
    {
        "age_stratum": "mid",
        "target_age_weeks": 52,
        "acceptable_window_weeks": "48-56",
        "purpose": "midlife interpolation support",
    },
    {
        "age_stratum": "late_mid",
        "target_age_weeks": 78,
        "acceptable_window_weeks": "72-84",
        "purpose": "pre-old transition support",
    },
    {
        "age_stratum": "old",
        "target_age_weeks": 112,
        "acceptable_window_weeks": "104-128",
        "purpose": "old target-tissue support for GSE121141 stress-test gap",
    },
]
SEXES = ["female", "male"]
PREFERRED_REPLICATES_PER_CELL = 4
MINIMUM_REPLICATES_PER_CELL = 3


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def cohort_design() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for tissue in TARGET_TISSUES:
        for age in AGE_STRATA:
            for sex in SEXES:
                rows.append(
                    {
                        "tissue": tissue,
                        "sex": sex,
                        "age_stratum": age["age_stratum"],
                        "target_age_weeks": age["target_age_weeks"],
                        "acceptable_window_weeks": age["acceptable_window_weeks"],
                        "minimum_replicates": MINIMUM_REPLICATES_PER_CELL,
                        "preferred_replicates": PREFERRED_REPLICATES_PER_CELL,
                        "minimum_samples": MINIMUM_REPLICATES_PER_CELL,
                        "preferred_samples": PREFERRED_REPLICATES_PER_CELL,
                        "purpose": age["purpose"],
                    }
                )
    return rows


def metadata_schema() -> list[dict[str, Any]]:
    required = [
        ("sample_id", "string", "required", "Stable unique sample ID used in files and sample sheet."),
        ("mouse_id", "string", "required", "Animal identifier; keep separate from sample_id if multiple tissues are collected."),
        ("age_days", "number", "required", "Exact age in days at collection."),
        ("age_weeks", "number", "required", "Exact age in weeks at collection; should equal age_days / 7."),
        ("raw_age_token", "string", "required", "Original age label from source records or sample sheet."),
        ("tissue", "enum", "required", "One of brain_cortex, heart, lung."),
        ("tissue_region_detail", "string", "required", "For cortex, specify dissection region and laterality when available."),
        ("sex", "enum", "required", "female or male; unknown is not acceptable for headline Route A."),
        ("strain", "string", "required", "Mouse strain, preferably a single controlled background such as C57BL/6J."),
        ("intervention", "string", "required", "control for the headline cohort; non-control samples must be separated."),
        ("diet", "string", "required", "Diet description or chow product if known."),
        ("housing_site", "string", "required", "Facility/site identifier for batch traceability."),
        ("collection_date", "date", "required", "Date of tissue collection."),
        ("extraction_batch", "string", "required", "DNA extraction batch."),
        ("library_batch", "string", "required", "Library preparation batch."),
        ("sequencing_batch", "string", "required", "Sequencing batch or lane group."),
        ("assay", "enum", "required", "bulk RRBS, WGBS, or another traceable bisulfite methylation assay."),
        ("genome_assembly", "enum", "required", "mm10/GRCm38 or mm39/GRCm39; liftover path must be documented."),
        ("fastq_r1", "path_or_uri", "expected", "Raw FASTQ R1 path or URI if raw data is delivered."),
        ("fastq_r2", "path_or_uri", "expected", "Raw FASTQ R2 path or URI for paired-end data."),
        ("fastq_sha256", "string", "expected", "Checksum manifest for raw files."),
        ("processed_coverage_path", "path_or_uri", "required", "Coverage-like CpG file or beta matrix path."),
        ("processed_schema", "string", "required", "Column schema and beta calculation definition."),
        ("processed_sha256", "string", "required", "Checksum for processed methylation file."),
        ("notes", "string", "optional", "Known deviations, pathology, low input, or QC comments."),
    ]
    return [
        {
            "field": field,
            "type": typ,
            "requirement": requirement,
            "description": description,
        }
        for field, typ, requirement, description in required
    ]


def gates() -> dict[str, Any]:
    return {
        "route": "A",
        "status": "rfc_only_not_authorized",
        "forbidden_by_this_rfc": {
            "wet_lab_authorized": False,
            "raw_fastq_download_authorized": False,
            "bismark_authorized": False,
            "training_authorized": False,
            "autoresearch_authorized": False,
        },
        "study_design_gate": {
            "species": "Mus musculus",
            "context": "bulk tissue only",
            "target_tissues": TARGET_TISSUES,
            "minimum_total_samples": 72,
            "preferred_total_samples": 96,
            "minimum_per_tissue": 24,
            "preferred_per_tissue": 32,
            "minimum_old_samples_per_tissue": 6,
            "preferred_old_samples_per_tissue": 8,
            "old_definition_weeks": ">=104",
            "sample_specific_exact_age_required": True,
        },
        "assay_gate": {
            "allowed": ["bulk RRBS", "bulk WGBS", "traceable bulk bisulfite methylation"],
            "blocked": ["single-cell", "cell-type targeted", "organoid", "in vitro", "low-coverage/iTAG"],
            "assembly_required": True,
            "processed_coverage_or_beta_required": True,
            "raw_file_checksums_expected": True,
        },
        "matrix_gate": {
            "sex_mt_exclusion_required": True,
            "five_kb_region_framework_required": True,
            "common_regions_with_v8_2_reference_min": 50000,
            "metadata_overlap_min": 0.95,
            "age_coverage_min": 0.95,
            "beta_range": [0.0, 1.0],
        },
        "benchmark_gate": {
            "fixed_initial_configs": [
                "lgbm quantile_uniform p=0.8 top1000 agebin=0.8",
                "lgbm robust p=0.95 top1000 shift=0.15",
            ],
            "must_run": [
                "GroupKFold by dataset_batch",
                "leave-one-dataset-out",
                "all_except:GSE121141 -> GSE121141",
                "all_except:GSE80672 -> GSE80672",
                "random-label sanity",
                "shuffled CR sanity when CR metrics are computed",
            ],
            "scientific_success_thresholds": {
                "GSE121141_old104_mae_improvement_weeks_vs_75_386": ">=10",
                "GSE121141_all_age_heldout_mae": "<=40.033",
                "GroupKFold_mae": "<=25.767",
                "random_label_abs_r": "<0.2",
                "shuffled_cr_auc": "0.4-0.6 when available",
            },
        },
    }


def write_report(design_rows: list[dict[str, Any]], gate_payload: dict[str, Any]) -> None:
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    preferred_total = sum(int(row["preferred_samples"]) for row in design_rows)
    minimum_total = sum(int(row["minimum_samples"]) for row in design_rows)
    lines = [
        "# Route A Old Target-Tissue Data Generation RFC",
        "",
        f"Date: {datetime.now().isoformat(timespec='seconds')}",
        "",
        "## Summary",
        "",
        "Route A is the preferred path after v13 Route C acceptance and the Route B backlog refresh.",
        "The current public/raw backlog did not produce a headline-ready minimal FASTQ pilot candidate.",
        "This RFC defines the minimum data-generation or collaboration package needed to test true old target-tissue generalization.",
        "",
        "This document does not authorize wet-lab work, data download, Bismark, model training, or autoresearch.",
        "",
        "## Scientific Target",
        "",
        "- Species: `Mus musculus`.",
        "- Context: bulk tissue only.",
        "- Target tissues: `brain_cortex`, `heart`, `lung`.",
        "- Target old-age support: sample-specific exact age with old samples at `>=104w`.",
        "- Preferred assay: bulk RRBS or WGBS with explicit genome assembly and processed methylation output.",
        "",
        "## Cohort Design",
        "",
        f"- Minimum viable design: `{minimum_total}` samples.",
        f"- Preferred design: `{preferred_total}` samples.",
        "- Preferred cell structure: 3 tissues x 4 age strata x 2 sexes x 4 replicates.",
        "- Minimum cell structure: 3 tissues x 4 age strata x 2 sexes x 3 replicates.",
        "- Old stratum target: `112w`, acceptable window `104-128w`.",
        "",
        "The cohort design table is written to `results/route_a_data_generation_rfc/route_a_cohort_design.csv`.",
        "",
        "## Required Metadata",
        "",
        "The required schema is written to `results/route_a_data_generation_rfc/route_a_required_metadata_schema.csv`.",
        "The key non-negotiable fields are `sample_id`, `mouse_id`, exact `age_days/age_weeks`, `tissue`, `sex`, `strain`, `intervention`, `assay`, `genome_assembly`, processed methylation path, and checksums.",
        "",
        "## Gate Criteria",
        "",
        "- Metadata overlap must be `>=95%`.",
        "- Exact age coverage must be `>=95%`.",
        "- At least `50,000` common 5kb regions with the v8.2/v13 reference must be available.",
        "- Sex chromosomes and MT must be removable.",
        "- Beta values must stay in `[0,1]`.",
        "- Random-label sanity must pass before any model claim.",
        "- CR or other biological-age metrics remain research-level and require real held-out predictions.",
        "",
        "## Promotion Rule",
        "",
        "Route A data may enter headline benchmarking only after metadata, schema, matrix, and sanity gates pass.",
        "Constrained autoresearch remains forbidden unless the new data improves old104+ support by at least `5w` and sanity checks pass.",
        "",
        "## Outputs",
        "",
        "- `results/route_a_data_generation_rfc/route_a_cohort_design.csv`",
        "- `results/route_a_data_generation_rfc/route_a_required_metadata_schema.csv`",
        "- `results/route_a_data_generation_rfc/route_a_qc_benchmark_gates.json`",
        "- `results/route_a_data_generation_rfc/route_a_decision_state.json`",
        "- `doc/30_protocols/06_20260519_route_a_submission_checklist.md`",
        "",
        "## Current Decision",
        "",
        f"`{gate_payload['status']}`. Separate approval is required before any data generation, raw download, Bismark, training, or autoresearch.",
    ]
    REPORT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_checklist() -> None:
    CHECKLIST_PATH.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Route A Submission Checklist",
        "",
        "Date: 2026-05-19",
        "",
        "Use this checklist before accepting generated or collaborative old target-tissue methylation data.",
        "",
        "## Pre-Acceptance",
        "",
        "- Confirm mouse bulk tissue context.",
        "- Confirm exact sample-specific age for every sample.",
        "- Confirm target tissue labels are one of `brain_cortex`, `heart`, `lung`.",
        "- Confirm old samples include `>=104w` for each target tissue intended for headline claims.",
        "- Confirm sex, strain, diet, intervention, batch, and collection metadata are present.",
        "- Confirm assay is bulk RRBS/WGBS or traceable bisulfite methylation.",
        "- Confirm genome assembly and coordinate system are explicit.",
        "- Confirm processed methylation schema is documented.",
        "- Confirm file checksums are supplied.",
        "",
        "## Adapter Smoke",
        "",
        "- Parse the first 1,000-10,000 rows or 2-3 sample files.",
        "- Verify coordinates, beta range, sample IDs, and metadata join.",
        "- Verify sex/MT exclusion can be applied.",
        "- Estimate common 5kb region count against v8.2/v13 reference.",
        "- Stop before training if common regions are `<50,000`.",
        "",
        "## Benchmark Promotion",
        "",
        "- Run fixed v7/v8 configs before any search.",
        "- Run GroupKFold, LODO, GSE121141 held-out, GSE80672 CR held-out, random-label sanity, and shuffled CR sanity when applicable.",
        "- Keep GSE121141 old104+ stress-test reporting even if Route A improves support.",
        "- Do not run autoresearch unless new-data gates and sanity checks pass.",
        "",
        "## Not Authorized By This Checklist",
        "",
        "- Wet-lab execution.",
        "- Raw FASTQ download.",
        "- Bismark/FASTQ ETL.",
        "- Model training.",
        "- Autoresearch.",
    ]
    CHECKLIST_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def update_index() -> None:
    if not INDEX_PATH.exists():
        return
    text = INDEX_PATH.read_text(encoding="utf-8")
    entries = [
        "- Route A concrete data-generation RFC: `doc/20_analysis/36_20260519_route_a_old_tissue_data_generation_rfc.md`",
        "- Route A submission checklist: `doc/30_protocols/06_20260519_route_a_submission_checklist.md`",
    ]
    marker = "- Route B manual backlog refresh report: `doc/20_analysis/35_20260519_route_b_manual_backlog_runinfo_refresh_report.md`"
    for entry in entries:
        if entry in text:
            continue
        if marker in text:
            text = text.replace(marker, marker + "\n" + entry)
        else:
            text += "\n" + entry + "\n"
    INDEX_PATH.write_text(text, encoding="utf-8")


def main() -> None:
    design_rows = cohort_design()
    schema_rows = metadata_schema()
    gate_payload = gates()
    state = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "route": "A",
        "status": "rfc_only_not_authorized",
        "minimum_total_samples": sum(int(row["minimum_samples"]) for row in design_rows),
        "preferred_total_samples": sum(int(row["preferred_samples"]) for row in design_rows),
        "target_tissues": TARGET_TISSUES,
        "old_definition_weeks": ">=104",
        "wet_lab_authorized": False,
        "download_authorized": False,
        "fastq_download_authorized": False,
        "bismark_authorized": False,
        "training_authorized": False,
        "autoresearch_authorized": False,
        "next_action": "seek_explicit_route_a_approval_or_collaboration_scope_before_data_generation",
    }
    write_csv(DESIGN_CSV, design_rows)
    write_csv(SCHEMA_CSV, schema_rows)
    write_json(GATES_JSON, gate_payload)
    write_json(STATE_JSON, state)
    write_report(design_rows, gate_payload)
    write_checklist()
    update_index()
    print(json.dumps(state, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

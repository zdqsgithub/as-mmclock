#!/usr/bin/env python3
"""Build the Route A partner handoff package.

The package contains sample-sheet templates, file-manifest templates, and a
protocol document for external collaborators. It does not authorize data
generation, downloads, Bismark, training, or autoresearch.
"""
from __future__ import annotations

import csv
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

ROOT = Path("/home/zdq-as/mouse_methyl_work")
TEMPLATE_DIR = ROOT / "metadata" / "templates"
OUT_DIR = ROOT / "results" / "route_a_handoff_pack"
REPORT_PATH = ROOT / "doc" / "30_protocols" / "07_20260519_route_a_partner_handoff_package.md"
INDEX_PATH = ROOT / "doc" / "00_meta" / "02_20260519_project_status_index_v13.md"

DESIGN = ROOT / "results" / "route_a_data_generation_rfc" / "route_a_cohort_design.csv"
SCHEMA = ROOT / "results" / "route_a_data_generation_rfc" / "route_a_required_metadata_schema.csv"

MIN_SAMPLE_SHEET = TEMPLATE_DIR / "route_a_sample_sheet_minimum_72_template.csv"
PREF_SAMPLE_SHEET = TEMPLATE_DIR / "route_a_sample_sheet_preferred_96_template.csv"
MIN_FILE_MANIFEST = TEMPLATE_DIR / "route_a_file_manifest_minimum_72_template.csv"
PREF_FILE_MANIFEST = TEMPLATE_DIR / "route_a_file_manifest_preferred_96_template.csv"
HANDOFF_MANIFEST = OUT_DIR / "route_a_handoff_manifest.json"
VALIDATION_REPORT = OUT_DIR / "route_a_template_validation_report.json"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows and fieldnames is None:
        path.write_text("", encoding="utf-8")
        return
    if fieldnames is None:
        fieldnames = []
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


def metadata_fields() -> list[str]:
    rows = read_csv(SCHEMA)
    fields = [row["field"] for row in rows]
    for extra in ["age_stratum", "replicate_index", "planned_sample_role"]:
        if extra not in fields:
            fields.append(extra)
    return fields


def build_sample_sheet(replicates_per_cell: int, label: str) -> list[dict[str, Any]]:
    design_rows = read_csv(DESIGN)
    rows: list[dict[str, Any]] = []
    for cell in design_rows:
        tissue = cell["tissue"]
        sex = cell["sex"]
        age_stratum = cell["age_stratum"]
        age_weeks = float(cell["target_age_weeks"])
        age_days = int(round(age_weeks * 7))
        for rep in range(1, replicates_per_cell + 1):
            sample_id = f"ROUTEA_{label}_{tissue.upper()}_{sex[0].upper()}_{int(age_weeks):03d}W_R{rep:02d}"
            mouse_id = f"MOUSE_{label}_{tissue.upper()}_{sex[0].upper()}_{int(age_weeks):03d}W_R{rep:02d}"
            rows.append(
                {
                    "sample_id": sample_id,
                    "mouse_id": mouse_id,
                    "age_days": age_days,
                    "age_weeks": age_weeks,
                    "raw_age_token": f"{age_weeks:g} weeks",
                    "tissue": tissue,
                    "tissue_region_detail": "TO_BE_FILLED",
                    "sex": sex,
                    "strain": "TO_BE_FILLED",
                    "intervention": "control",
                    "diet": "TO_BE_FILLED",
                    "housing_site": "TO_BE_FILLED",
                    "collection_date": "TO_BE_FILLED",
                    "extraction_batch": "TO_BE_FILLED",
                    "library_batch": "TO_BE_FILLED",
                    "sequencing_batch": "TO_BE_FILLED",
                    "assay": "bulk RRBS or WGBS",
                    "genome_assembly": "TO_BE_FILLED",
                    "fastq_r1": "TO_BE_FILLED",
                    "fastq_r2": "TO_BE_FILLED",
                    "fastq_sha256": "TO_BE_FILLED",
                    "processed_coverage_path": "TO_BE_FILLED",
                    "processed_schema": "TO_BE_FILLED",
                    "processed_sha256": "TO_BE_FILLED",
                    "notes": "",
                    "age_stratum": age_stratum,
                    "replicate_index": rep,
                    "planned_sample_role": f"{label}_route_a_template",
                }
            )
    return rows


def build_file_manifest(sample_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for sample in sample_rows:
        rows.append(
            {
                "sample_id": sample["sample_id"],
                "file_role": "processed_coverage",
                "path_or_uri": "TO_BE_FILLED",
                "sha256": "TO_BE_FILLED",
                "file_size_bytes": "TO_BE_FILLED",
                "genome_assembly": "TO_BE_FILLED",
                "processed_schema": "TO_BE_FILLED",
                "notes": "",
            }
        )
    return rows


def run_template_validation() -> dict[str, Any]:
    cmd = [
        sys.executable,
        str(ROOT / "scripts" / "validate" / "validate_route_a_submission.py"),
        "--sample-sheet",
        str(MIN_SAMPLE_SHEET),
        "--file-manifest",
        str(MIN_FILE_MANIFEST),
        "--allow-placeholders",
        "--output",
        str(VALIDATION_REPORT),
    ]
    completed = subprocess.run(cmd, cwd=ROOT, check=False, text=True, capture_output=True)
    return {
        "command": " ".join(cmd),
        "returncode": completed.returncode,
        "stdout_tail": completed.stdout[-2000:],
        "stderr_tail": completed.stderr[-2000:],
    }


def write_report(manifest: dict[str, Any]) -> None:
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Route A Partner Handoff Package",
        "",
        f"Date: {datetime.now().isoformat(timespec='seconds')}",
        "",
        "## Purpose",
        "",
        "This package converts the Route A RFC into files that can be sent to a collaborator or experimental team before any data generation or transfer.",
        "It is a metadata and planning package only.",
        "",
        "## Included Files",
        "",
        "- `metadata/templates/route_a_sample_sheet_minimum_72_template.csv`",
        "- `metadata/templates/route_a_sample_sheet_preferred_96_template.csv`",
        "- `metadata/templates/route_a_file_manifest_minimum_72_template.csv`",
        "- `metadata/templates/route_a_file_manifest_preferred_96_template.csv`",
        "- `scripts/validate/validate_route_a_submission.py`",
        "- `results/route_a_handoff_pack/route_a_handoff_manifest.json`",
        "- `results/route_a_handoff_pack/route_a_template_validation_report.json`",
        "",
        "## How To Use",
        "",
        "1. Fill the sample sheet with exact age, tissue, sex, strain, intervention, batch, assay, genome assembly, and processed methylation metadata.",
        "2. Fill the file manifest with processed coverage/beta file paths and checksums.",
        "3. Run the validator before accepting the submission:",
        "",
        "```bash",
        "VIRTUAL_ENV=.venv-core uv run --active python scripts/validate/validate_route_a_submission.py \\",
        "  --sample-sheet metadata/templates/route_a_sample_sheet_minimum_72_template.csv \\",
        "  --file-manifest metadata/templates/route_a_file_manifest_minimum_72_template.csv",
        "```",
        "",
        "For a template dry check only, add `--allow-placeholders`.",
        "",
        "## Current Template Scope",
        "",
        f"- Minimum samples: `{manifest['minimum_template_samples']}`.",
        f"- Preferred samples: `{manifest['preferred_template_samples']}`.",
        "- Target tissues: `brain_cortex`, `heart`, `lung`.",
        "- Age strata: `12w`, `52w`, `78w`, `112w`.",
        "- Old-age support threshold: `>=104w`.",
        "",
        "## Guardrails",
        "",
        "- This package does not authorize wet-lab execution.",
        "- This package does not authorize raw FASTQ download.",
        "- This package does not authorize Bismark/FASTQ ETL.",
        "- This package does not authorize model training or autoresearch.",
        "- Any real submission must pass the validator and then adapter/matrix gates before benchmarking.",
    ]
    REPORT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def update_index() -> None:
    if not INDEX_PATH.exists():
        return
    text = INDEX_PATH.read_text(encoding="utf-8")
    entry = "- Route A partner handoff package: `doc/30_protocols/07_20260519_route_a_partner_handoff_package.md`"
    if entry in text:
        return
    marker = "- Route A submission checklist: `doc/30_protocols/06_20260519_route_a_submission_checklist.md`"
    if marker in text:
        text = text.replace(marker, marker + "\n" + entry)
    else:
        text += "\n" + entry + "\n"
    INDEX_PATH.write_text(text, encoding="utf-8")


def main() -> None:
    minimum_rows = build_sample_sheet(3, "MIN72")
    preferred_rows = build_sample_sheet(4, "PREF96")
    fields = metadata_fields()
    write_csv(MIN_SAMPLE_SHEET, minimum_rows, fields)
    write_csv(PREF_SAMPLE_SHEET, preferred_rows, fields)
    write_csv(MIN_FILE_MANIFEST, build_file_manifest(minimum_rows))
    write_csv(PREF_FILE_MANIFEST, build_file_manifest(preferred_rows))
    validation = run_template_validation()
    manifest = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "status": "handoff_package_created",
        "minimum_template_samples": len(minimum_rows),
        "preferred_template_samples": len(preferred_rows),
        "sample_sheet_minimum": str(MIN_SAMPLE_SHEET),
        "sample_sheet_preferred": str(PREF_SAMPLE_SHEET),
        "file_manifest_minimum": str(MIN_FILE_MANIFEST),
        "file_manifest_preferred": str(PREF_FILE_MANIFEST),
        "validator": "scripts/validate/validate_route_a_submission.py",
        "template_validation": validation,
        "wet_lab_authorized": False,
        "download_authorized": False,
        "bismark_authorized": False,
        "training_authorized": False,
        "autoresearch_authorized": False,
    }
    write_json(HANDOFF_MANIFEST, manifest)
    write_report(manifest)
    update_index()
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

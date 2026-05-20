#!/usr/bin/env python3
"""Validate a Route A incoming sample sheet and optional file manifest.

This validator checks metadata and manifest gates only. It does not read data
files, download data, run Bismark, build matrices, train models, or authorize a
submission.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any

ROOT = Path("/home/zdq-as/mouse_methyl_work")
DEFAULT_SCHEMA = ROOT / "results" / "route_a_data_generation_rfc" / "route_a_required_metadata_schema.csv"
DEFAULT_OUT = ROOT / "results" / "route_a_handoff_pack" / "route_a_submission_validation_report.json"

TARGET_TISSUES = {"brain_cortex", "heart", "lung"}
ALLOWED_SEX = {"female", "male"}
ALLOWED_ASSAY_TERMS = ("rrbs", "wgbs", "bisulfite")
PLACEHOLDER_TOKENS = {"", "to_be_filled", "tbd", "na", "n/a", "placeholder"}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def required_fields(schema_path: Path) -> list[str]:
    rows = read_csv(schema_path)
    return [row["field"] for row in rows if row.get("requirement") == "required"]


def is_placeholder(value: Any) -> bool:
    return str(value or "").strip().lower() in PLACEHOLDER_TOKENS


def safe_float(value: Any) -> float | None:
    try:
        text = str(value or "").strip()
        if not text:
            return None
        return float(text)
    except ValueError:
        return None


def add(errors: list[dict[str, Any]], severity: str, code: str, message: str, row: int | None = None) -> None:
    errors.append({"severity": severity, "code": code, "message": message, "row": row})


def validate_sample_sheet(rows: list[dict[str, str]], required: list[str], allow_placeholders: bool) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    issues: list[dict[str, Any]] = []
    if not rows:
        add(issues, "error", "empty_sample_sheet", "Sample sheet has no rows.")
        return issues, {}

    columns = set(rows[0].keys())
    for field in required:
        if field not in columns:
            add(issues, "error", "missing_required_column", f"Missing required column: {field}")

    sample_ids: list[str] = []
    by_tissue: Counter[str] = Counter()
    old_by_tissue: Counter[str] = Counter()
    by_cell: Counter[tuple[str, str, str]] = Counter()
    age_coverage_n = 0

    for idx, row in enumerate(rows, start=2):
        sample_id = str(row.get("sample_id", "")).strip()
        sample_ids.append(sample_id)
        for field in required:
            value = row.get(field, "")
            if is_placeholder(value) and not allow_placeholders:
                add(issues, "error", "missing_required_value", f"Missing required value for {field}.", idx)

        tissue = str(row.get("tissue", "")).strip()
        sex = str(row.get("sex", "")).strip().lower()
        assay = str(row.get("assay", "")).strip().lower()
        intervention = str(row.get("intervention", "")).strip().lower()
        age_days = safe_float(row.get("age_days"))
        age_weeks = safe_float(row.get("age_weeks"))

        if sample_id and sample_id.lower() not in PLACEHOLDER_TOKENS:
            pass
        elif not allow_placeholders:
            add(issues, "error", "invalid_sample_id", "sample_id is empty or placeholder.", idx)

        if tissue not in TARGET_TISSUES and not allow_placeholders:
            add(issues, "error", "invalid_tissue", f"Invalid tissue `{tissue}`; expected one of {sorted(TARGET_TISSUES)}.", idx)
        if sex not in ALLOWED_SEX and not allow_placeholders:
            add(issues, "error", "invalid_sex", f"Invalid sex `{sex}`; expected female or male.", idx)
        if assay and not any(term in assay for term in ALLOWED_ASSAY_TERMS) and not allow_placeholders:
            add(issues, "error", "invalid_assay", "Assay must be bulk RRBS, WGBS, or traceable bisulfite methylation.", idx)
        if intervention and intervention != "control" and not allow_placeholders:
            add(issues, "warning", "non_control_intervention", "Non-control samples cannot enter headline Route A without separate stratification.", idx)

        if age_days is None or age_weeks is None:
            if not allow_placeholders:
                add(issues, "error", "missing_exact_age", "age_days and age_weeks must both be numeric.", idx)
        else:
            age_coverage_n += 1
            if not math.isclose(age_days / 7.0, age_weeks, abs_tol=0.25):
                add(issues, "error", "age_unit_mismatch", "age_days / 7 must match age_weeks within 0.25 weeks.", idx)

        if tissue in TARGET_TISSUES:
            by_tissue[tissue] += 1
            if age_weeks is not None and age_weeks >= 104:
                old_by_tissue[tissue] += 1
            age_stratum = str(row.get("age_stratum", infer_age_stratum(age_weeks))).strip() or infer_age_stratum(age_weeks)
            by_cell[(tissue, sex, age_stratum)] += 1

    duplicated = sorted(sample for sample, count in Counter(sample_ids).items() if sample and count > 1)
    if duplicated:
        add(issues, "error", "duplicate_sample_id", f"Duplicate sample_id values: {duplicated[:20]}")

    age_coverage = age_coverage_n / len(rows)
    if age_coverage < 0.95 and not allow_placeholders:
        add(issues, "error", "age_coverage_below_95_percent", f"Exact age coverage is {age_coverage:.3f}; expected >=0.95.")

    if len(rows) < 72 and not allow_placeholders:
        add(issues, "error", "sample_count_below_minimum", f"Sample count is {len(rows)}; Route A minimum is 72.")

    for tissue in sorted(TARGET_TISSUES):
        if by_tissue[tissue] < 24 and not allow_placeholders:
            add(issues, "error", "tissue_count_below_minimum", f"{tissue} has {by_tissue[tissue]} samples; expected >=24.")
        if old_by_tissue[tissue] < 6 and not allow_placeholders:
            add(issues, "error", "old_tissue_count_below_minimum", f"{tissue} has {old_by_tissue[tissue]} old samples; expected >=6.")

    metrics = {
        "n_samples": len(rows),
        "age_coverage": round(age_coverage, 4),
        "tissue_counts": dict(sorted(by_tissue.items())),
        "old_tissue_counts": dict(sorted(old_by_tissue.items())),
        "cell_counts_min": min(by_cell.values()) if by_cell else 0,
        "cell_counts_max": max(by_cell.values()) if by_cell else 0,
    }
    return issues, metrics


def infer_age_stratum(age_weeks: float | None) -> str:
    if age_weeks is None:
        return ""
    if age_weeks < 24:
        return "young"
    if age_weeks < 64:
        return "mid"
    if age_weeks < 104:
        return "late_mid"
    return "old"


def validate_file_manifest(rows: list[dict[str, str]], sample_ids: set[str], allow_placeholders: bool) -> list[dict[str, Any]]:
    issues: list[dict[str, Any]] = []
    required_columns = {"sample_id", "file_role", "path_or_uri", "sha256", "file_size_bytes", "genome_assembly", "processed_schema"}
    if not rows:
        add(issues, "warning", "empty_file_manifest", "File manifest has no rows.")
        return issues
    columns = set(rows[0].keys())
    for field in sorted(required_columns - columns):
        add(issues, "error", "missing_manifest_column", f"Missing file manifest column: {field}")
    for idx, row in enumerate(rows, start=2):
        sample_id = str(row.get("sample_id", "")).strip()
        if sample_id not in sample_ids and not allow_placeholders:
            add(issues, "error", "manifest_sample_not_in_sheet", f"Manifest sample_id `{sample_id}` is not in sample sheet.", idx)
        for field in ["file_role", "path_or_uri", "sha256", "genome_assembly", "processed_schema"]:
            if is_placeholder(row.get(field)) and not allow_placeholders:
                add(issues, "error", "manifest_missing_required_value", f"Missing manifest value for {field}.", idx)
        size = safe_float(row.get("file_size_bytes"))
        if size is None and not allow_placeholders:
            add(issues, "error", "manifest_invalid_file_size", "file_size_bytes must be numeric.", idx)
    processed_sample_ids = {
        str(row.get("sample_id", "")).strip()
        for row in rows
        if "processed" in str(row.get("file_role", "")).lower()
        or "coverage" in str(row.get("file_role", "")).lower()
        or "beta" in str(row.get("file_role", "")).lower()
    }
    missing_processed = sorted(sample_ids - processed_sample_ids)
    if missing_processed and not allow_placeholders:
        add(
            issues,
            "error",
            "sample_without_processed_manifest_row",
            f"{len(missing_processed)} sample_id values have no processed methylation manifest row; examples: {missing_processed[:10]}",
        )
    return issues


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-sheet", required=True, type=Path)
    parser.add_argument("--file-manifest", type=Path)
    parser.add_argument("--schema", type=Path, default=DEFAULT_SCHEMA)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--allow-placeholders", action="store_true", help="Allow template placeholders for dry template checks.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    sample_rows = read_csv(args.sample_sheet)
    required = required_fields(args.schema)
    issues, metrics = validate_sample_sheet(sample_rows, required, args.allow_placeholders)
    if args.file_manifest:
        manifest_rows = read_csv(args.file_manifest)
        issues.extend(validate_file_manifest(manifest_rows, {row.get("sample_id", "") for row in sample_rows}, args.allow_placeholders))
    errors = [issue for issue in issues if issue["severity"] == "error"]
    warnings = [issue for issue in issues if issue["severity"] == "warning"]
    payload = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "sample_sheet": str(args.sample_sheet),
        "file_manifest": str(args.file_manifest) if args.file_manifest else "",
        "allow_placeholders": args.allow_placeholders,
        "status": "passed" if not errors else "failed",
        "error_count": len(errors),
        "warning_count": len(warnings),
        "metrics": metrics,
        "issues": issues,
        "download_authorized": False,
        "bismark_authorized": False,
        "training_authorized": False,
        "autoresearch_authorized": False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))
    raise SystemExit(0 if not errors else 2)


if __name__ == "__main__":
    main()

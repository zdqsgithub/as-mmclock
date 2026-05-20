#!/usr/bin/env python3
"""Build a Route A file manifest from local staged processed files.

This helper reads a filled Route A sample sheet, computes local file size and
sha256 for each processed methylation file, and writes a file manifest accepted
by the Route A gates. It never downloads data, runs Bismark, trains models, or
starts autoresearch.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime
from pathlib import Path
from typing import Any

ROOT = Path("/home/zdq-as/mouse_methyl_work")
DEFAULT_OUT_BASE = ROOT / "results" / "route_a_local_file_manifest"
REMOTE_PREFIXES = ("http://", "https://", "ftp://", "s3://", "gs://")
PLACEHOLDER_TOKENS = {"", "to_be_filled", "tbd", "na", "n/a", "placeholder", "not_provided_processed_only"}
SUPPORTED_SCHEMAS = {"bismark_cov_6col", "cpg_beta_table", "region_beta_matrix"}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["sample_id", "file_role", "path_or_uri", "sha256", "file_size_bytes", "genome_assembly", "processed_schema", "notes"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str), encoding="utf-8")


def is_placeholder(value: Any) -> bool:
    return str(value or "").strip().lower() in PLACEHOLDER_TOKENS


def is_remote(value: Any) -> bool:
    return str(value or "").strip().lower().startswith(REMOTE_PREFIXES)


def normalize_schema(value: Any) -> str:
    text = str(value or "").strip().lower().replace("-", "_").replace(" ", "_")
    aliases = {
        "bismark": "bismark_cov_6col",
        "bismark_cov": "bismark_cov_6col",
        "bismark_cov_6col": "bismark_cov_6col",
        "cpg_beta": "cpg_beta_table",
        "cpg_beta_table": "cpg_beta_table",
        "beta_table": "cpg_beta_table",
        "region_beta": "region_beta_matrix",
        "region_beta_matrix": "region_beta_matrix",
        "5kb_region_beta_matrix": "region_beta_matrix",
    }
    return aliases.get(text, text)


def resolve_local_path(value: str) -> Path:
    if is_placeholder(value):
        raise ValueError("processed_path_is_placeholder")
    if is_remote(value):
        raise ValueError("remote_uri_not_supported")
    path = Path(str(value))
    return path if path.is_absolute() else ROOT / path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build(args: argparse.Namespace) -> dict[str, Any]:
    rows = read_csv(args.sample_sheet)
    errors: list[dict[str, Any]] = []
    warnings: list[dict[str, Any]] = []
    manifest_rows: list[dict[str, Any]] = []
    file_cache: dict[str, dict[str, Any]] = {}
    seen_sample_ids: set[str] = set()

    for idx, row in enumerate(rows, start=2):
        sample_id = str(row.get(args.sample_id_column, "")).strip()
        if not sample_id or sample_id in seen_sample_ids:
            errors.append({"row": idx, "code": "missing_or_duplicate_sample_id", "sample_id": sample_id})
            continue
        seen_sample_ids.add(sample_id)

        schema = normalize_schema(row.get(args.schema_column, ""))
        if schema not in SUPPORTED_SCHEMAS:
            errors.append({"row": idx, "sample_id": sample_id, "code": "unsupported_or_missing_processed_schema", "processed_schema": row.get(args.schema_column, "")})
            continue

        assembly = str(row.get(args.assembly_column, "")).strip()
        if is_placeholder(assembly):
            errors.append({"row": idx, "sample_id": sample_id, "code": "missing_genome_assembly"})
            continue

        path_text = str(row.get(args.path_column, "")).strip()
        try:
            path = resolve_local_path(path_text)
        except ValueError as exc:
            errors.append({"row": idx, "sample_id": sample_id, "code": str(exc), "path_or_uri": path_text})
            continue
        if not path.exists():
            errors.append({"row": idx, "sample_id": sample_id, "code": "local_processed_file_not_found", "path_or_uri": str(path)})
            continue
        if not path.is_file():
            errors.append({"row": idx, "sample_id": sample_id, "code": "local_processed_path_not_file", "path_or_uri": str(path)})
            continue

        key = str(path)
        if key not in file_cache:
            file_cache[key] = {"file_size_bytes": path.stat().st_size, "sha256": sha256_file(path)}
        file_info = file_cache[key]
        manifest_rows.append(
            {
                "sample_id": sample_id,
                "file_role": args.file_role,
                "path_or_uri": str(path),
                "sha256": file_info["sha256"],
                "file_size_bytes": int(file_info["file_size_bytes"]),
                "genome_assembly": assembly,
                "processed_schema": schema,
                "notes": f"generated_from_sample_sheet:{args.path_column}",
            }
        )

        existing_sha = str(row.get("processed_sha256", "")).strip().lower()
        if existing_sha and not is_placeholder(existing_sha) and existing_sha != str(file_info["sha256"]).lower():
            warnings.append({"row": idx, "sample_id": sample_id, "code": "sample_sheet_processed_sha256_differs_from_computed_manifest_sha256"})

    write_csv(args.output, manifest_rows)
    status = "passed" if not errors else "failed"
    report = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "status": status,
        "sample_sheet": str(args.sample_sheet),
        "output_manifest": str(args.output),
        "path_column": args.path_column,
        "schema_column": args.schema_column,
        "assembly_column": args.assembly_column,
        "file_role": args.file_role,
        "n_sample_sheet_rows": len(rows),
        "n_manifest_rows": len(manifest_rows),
        "n_unique_files": len(file_cache),
        "error_count": len(errors),
        "warning_count": len(warnings),
        "errors": errors[:100],
        "warnings": warnings[:100],
        "download_authorized": False,
        "bismark_authorized": False,
        "training_authorized": False,
        "autoresearch_authorized": False,
    }
    write_json(args.report, report)
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-sheet", required=True, type=Path)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--report", type=Path, default=None)
    parser.add_argument("--submission-id", default="route_a_submission")
    parser.add_argument("--sample-id-column", default="sample_id")
    parser.add_argument("--path-column", default="processed_coverage_path")
    parser.add_argument("--schema-column", default="processed_schema")
    parser.add_argument("--assembly-column", default="genome_assembly")
    parser.add_argument("--file-role", default="processed_methylation")
    args = parser.parse_args()
    out_dir = DEFAULT_OUT_BASE / args.submission_id
    if args.output is None:
        args.output = out_dir / "route_a_file_manifest.csv"
    if args.report is None:
        args.report = out_dir / "route_a_file_manifest_build_report.json"
    return args


def main() -> None:
    args = parse_args()
    report = build(args)
    print(json.dumps(report, indent=2, sort_keys=True, default=str))
    raise SystemExit(0 if report.get("status") == "passed" else 2)


if __name__ == "__main__":
    main()

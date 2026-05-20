#!/usr/bin/env python3
"""Run the Route A gated intake workflow.

This controller enforces the Route A sequence:

1. metadata gate with validate_route_a_submission.py
2. adapter smoke only after metadata passes
3. matrix gate only after adapter smoke passes
4. RALPH Learn/Benchmark readiness only after matrix gate passes

It never downloads data, runs Bismark, trains models, or starts autoresearch.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

ROOT = Path("/home/zdq-as/mouse_methyl_work")
DEFAULT_OUT = ROOT / "results" / "route_a_gate_workflow"
VALIDATOR = ROOT / "scripts" / "validate" / "validate_route_a_submission.py"
SOP_PATH = ROOT / "doc" / "30_protocols" / "08_20260519_route_a_gate_workflow_sop.md"
DRY_REPORT_PATH = ROOT / "doc" / "20_analysis" / "37_20260519_route_a_gate_workflow_dry_run_report.md"
INTAKE_REPORT_PATH = ROOT / "doc" / "20_analysis" / "38_20260519_route_a_intake_gate_readiness_report.md"
INDEX_PATH = ROOT / "doc" / "00_meta" / "02_20260519_project_status_index_v13.md"

PLACEHOLDER_TOKENS = {"", "to_be_filled", "tbd", "na", "n/a", "placeholder"}
LOCAL_URI_PREFIXES = ("http://", "https://", "ftp://", "s3://", "gs://")
COORD_CHROM_COLUMNS = {"chr", "chrom", "chromosome", "seqnames"}
COORD_POS_COLUMNS = {"pos", "position", "start", "cpg_pos", "coordinate"}
BETA_COLUMNS = {"beta", "methylation", "methylation_level", "percent_methylation", "meth_percent"}
SUPPORTED_SCHEMAS = {"bismark_cov_6col", "cpg_beta_table", "region_beta_matrix"}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def is_placeholder(value: Any) -> bool:
    return str(value or "").strip().lower() in PLACEHOLDER_TOKENS


def run_metadata_gate(args: argparse.Namespace, out_dir: Path) -> dict[str, Any]:
    output = out_dir / "metadata_gate_report.json"
    cmd = [
        sys.executable,
        str(VALIDATOR),
        "--sample-sheet",
        str(args.sample_sheet),
        "--output",
        str(output),
    ]
    if args.file_manifest:
        cmd.extend(["--file-manifest", str(args.file_manifest)])
    if args.allow_placeholders:
        cmd.append("--allow-placeholders")
    completed = subprocess.run(cmd, cwd=ROOT, check=False, text=True, capture_output=True)
    if output.exists():
        payload = json.loads(output.read_text(encoding="utf-8"))
    else:
        payload = {
            "status": "failed",
            "error_count": 1,
            "issues": [
                {
                    "severity": "error",
                    "code": "metadata_validator_did_not_write_output",
                    "message": completed.stderr[-1000:] or completed.stdout[-1000:],
                }
            ],
        }
    payload["command"] = " ".join(cmd)
    payload["returncode"] = completed.returncode
    payload["stdout_tail"] = completed.stdout[-2000:]
    payload["stderr_tail"] = completed.stderr[-2000:]
    write_json(output, payload)
    return payload


def adapter_smoke(args: argparse.Namespace, metadata_gate: dict[str, Any], out_dir: Path) -> dict[str, Any]:
    if metadata_gate.get("status") != "passed":
        payload = {
            "status": "blocked",
            "reason": "metadata_gate_failed",
            "adapter_smoke_pass": False,
            "download_authorized": False,
            "bismark_authorized": False,
        }
        write_json(out_dir / "adapter_smoke_report.json", payload)
        return payload

    if args.allow_placeholders:
        payload = {
            "status": "pending_real_submission",
            "reason": "template_or_placeholder_mode",
            "adapter_smoke_pass": False,
            "next_action": "replace placeholders with real processed methylation paths and rerun without --allow-placeholders",
            "download_authorized": False,
            "bismark_authorized": False,
        }
        write_json(out_dir / "adapter_smoke_report.json", payload)
        return payload

    if not args.file_manifest:
        payload = {
            "status": "pending_file_manifest",
            "reason": "no_file_manifest_provided",
            "adapter_smoke_pass": False,
            "download_authorized": False,
            "bismark_authorized": False,
        }
        write_json(out_dir / "adapter_smoke_report.json", payload)
        return payload

    manifest = read_csv(args.file_manifest)
    processed = [
        row
        for row in manifest
        if "processed" in str(row.get("file_role", "")).lower()
        or "coverage" in str(row.get("file_role", "")).lower()
        or "beta" in str(row.get("file_role", "")).lower()
    ]
    file_reports: list[dict[str, Any]] = []
    local_checked = 0
    passed_files = 0
    for row in processed[: args.max_smoke_files]:
        report = smoke_file(row, args.max_smoke_rows)
        file_reports.append(report)
        if report["status"] != "remote_or_missing_path":
            local_checked += 1
        if report.get("file_smoke_pass"):
            passed_files += 1

    if not processed:
        status = "pending_processed_methylation_manifest_rows"
        reason = "file_manifest_has_no_processed_coverage_or_beta_rows"
    elif local_checked == 0:
        status = "pending_local_or_staged_processed_files"
        reason = "processed paths are remote, missing, or placeholders; stage files before adapter smoke"
    elif passed_files == 0:
        status = "failed"
        reason = "no_processed_file_passed_schema_smoke"
    else:
        status = "passed"
        reason = "at_least_one_processed_file_passed_schema_smoke"

    payload = {
        "status": status,
        "reason": reason,
        "adapter_smoke_pass": status == "passed",
        "processed_manifest_rows": len(processed),
        "local_checked_files": local_checked,
        "passed_files": passed_files,
        "file_reports": file_reports,
        "download_authorized": False,
        "bismark_authorized": False,
    }
    write_json(out_dir / "adapter_smoke_report.json", payload)
    return payload


def smoke_file(row: dict[str, str], max_rows: int) -> dict[str, Any]:
    path_text = str(row.get("path_or_uri", "")).strip()
    if is_placeholder(path_text) or path_text.startswith(LOCAL_URI_PREFIXES):
        return {
            "sample_id": row.get("sample_id", ""),
            "path_or_uri": path_text,
            "status": "remote_or_missing_path",
            "file_smoke_pass": False,
            "reason": "path_is_placeholder_or_remote_uri",
        }
    path = Path(path_text)
    if not path.is_absolute():
        path = ROOT / path
    if not path.exists():
        return {
            "sample_id": row.get("sample_id", ""),
            "path_or_uri": path_text,
            "status": "remote_or_missing_path",
            "file_smoke_pass": False,
            "reason": "local_path_not_found",
        }
    schema = normalize_schema(row.get("processed_schema"))
    if schema == "bismark_cov_6col":
        return smoke_bismark_cov(row, path, max_rows)
    if schema and schema not in SUPPORTED_SCHEMAS:
        return {
            "sample_id": row.get("sample_id", ""),
            "path_or_uri": path_text,
            "status": "failed",
            "file_smoke_pass": False,
            "reason": f"unsupported_processed_schema:{schema}",
        }
    try:
        header, records = read_head_records(path, max_rows)
    except Exception as exc:  # noqa: BLE001 - smoke must report blockers and continue
        return {
            "sample_id": row.get("sample_id", ""),
            "path_or_uri": path_text,
            "status": "failed",
            "file_smoke_pass": False,
            "reason": f"read_failed:{type(exc).__name__}:{str(exc)[:200]}",
        }

    lower_header = [col.strip().lower() for col in header]
    has_region = "region_id" in lower_header
    has_coordinate = bool(COORD_CHROM_COLUMNS & set(lower_header)) and bool(COORD_POS_COLUMNS & set(lower_header))
    beta_values = extract_beta_values(header, records)
    beta_range_pass = bool(beta_values) and min(beta_values) >= 0.0 and max(beta_values) <= 1.0
    schema_pass = (has_region or has_coordinate) and beta_range_pass
    return {
        "sample_id": row.get("sample_id", ""),
        "path_or_uri": path_text,
        "status": "passed" if schema_pass else "failed",
        "file_smoke_pass": schema_pass,
        "n_rows_checked": len(records),
        "n_columns": len(header),
        "has_region_id": has_region,
        "has_coordinate_columns": has_coordinate,
        "beta_values_checked": len(beta_values),
        "beta_min": min(beta_values) if beta_values else "",
        "beta_max": max(beta_values) if beta_values else "",
        "reason": "schema_and_beta_range_passed" if schema_pass else "missing_coordinate_or_beta_range_failed",
    }


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


def smoke_bismark_cov(row: dict[str, str], path: Path, max_rows: int) -> dict[str, Any]:
    beta_values = []
    rows_checked = 0
    parseable = 0
    opener = gzip.open if str(path).endswith(".gz") else open
    try:
        with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
            for _, line in zip(range(max_rows), handle):
                rows_checked += 1
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 6:
                    parts = line.rstrip("\n").split(",")
                if len(parts) < 6:
                    continue
                chrom = str(parts[0]).strip()
                if chrom in {"chrX", "chrY", "chrM", "X", "Y", "M", "MT"}:
                    continue
                try:
                    meth = float(parts[4])
                    unmeth = float(parts[5])
                except ValueError:
                    continue
                total = meth + unmeth
                if total <= 0:
                    continue
                parseable += 1
                beta_values.append(meth / total)
    except Exception as exc:  # noqa: BLE001
        return {
            "sample_id": row.get("sample_id", ""),
            "path_or_uri": str(path),
            "status": "failed",
            "file_smoke_pass": False,
            "reason": f"read_failed:{type(exc).__name__}:{str(exc)[:200]}",
        }
    beta_range_pass = bool(beta_values) and min(beta_values) >= 0.0 and max(beta_values) <= 1.0
    return {
        "sample_id": row.get("sample_id", ""),
        "path_or_uri": str(path),
        "status": "passed" if beta_range_pass else "failed",
        "file_smoke_pass": beta_range_pass,
        "n_rows_checked": rows_checked,
        "n_columns": 6,
        "has_region_id": False,
        "has_coordinate_columns": True,
        "beta_values_checked": len(beta_values),
        "parseable_rows": parseable,
        "beta_min": min(beta_values) if beta_values else "",
        "beta_max": max(beta_values) if beta_values else "",
        "reason": "bismark_cov_schema_and_beta_range_passed" if beta_range_pass else "bismark_cov_beta_range_failed",
    }


def read_head_records(path: Path, max_rows: int) -> tuple[list[str], list[dict[str, str]]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        delimiter = "\t" if sample.count("\t") >= sample.count(",") else ","
        reader = csv.DictReader(handle, delimiter=delimiter)
        records = []
        for _, row in zip(range(max_rows), reader):
            records.append(row)
        return list(reader.fieldnames or []), records


def extract_beta_values(header: list[str], records: list[dict[str, str]]) -> list[float]:
    lower_to_original = {col.strip().lower(): col for col in header}
    beta_cols = [lower_to_original[col] for col in BETA_COLUMNS if col in lower_to_original]
    if not beta_cols:
        numeric_cols = []
        blocked = COORD_CHROM_COLUMNS | COORD_POS_COLUMNS | {"region_id", "chrom", "chr"}
        for col in header:
            if col.strip().lower() not in blocked:
                numeric_cols.append(col)
        beta_cols = numeric_cols[: min(20, len(numeric_cols))]
    values: list[float] = []
    for row in records:
        for col in beta_cols:
            text = str(row.get(col, "")).strip()
            if text in {"", "NA", "NaN", "nan"}:
                continue
            try:
                value = float(text)
            except ValueError:
                continue
            values.append(value)
            if len(values) >= 1000:
                return values
    return values


def matrix_gate(args: argparse.Namespace, adapter: dict[str, Any], out_dir: Path) -> dict[str, Any]:
    if adapter.get("status") != "passed":
        payload = {
            "status": "blocked",
            "reason": "adapter_smoke_not_passed",
            "matrix_gate_pass": False,
        }
        write_json(out_dir / "matrix_gate_report.json", payload)
        return payload
    if not args.matrix_manifest:
        payload = {
            "status": "pending_matrix_manifest",
            "reason": "adapter_passed_but_no_matrix_manifest_provided",
            "matrix_gate_pass": False,
        }
        write_json(out_dir / "matrix_gate_report.json", payload)
        return payload
    manifest = json.loads(args.matrix_manifest.read_text(encoding="utf-8"))
    common_regions = number_from_manifest(manifest, ["common_regions", "common_region_count", "n_common_regions"])
    metadata_overlap = number_from_manifest(manifest, ["metadata_overlap", "metadata_overlap_rate"])
    age_coverage = number_from_manifest(manifest, ["age_coverage", "age_coverage_rate"])
    beta_min = number_from_manifest(manifest, ["beta_min", "min_beta"])
    beta_max = number_from_manifest(manifest, ["beta_max", "max_beta"])
    issues = []
    if manifest.get("status") != "completed":
        issues.append(f"manifest_status_{manifest.get('status', 'missing')}")
    if common_regions is None or common_regions < 50000:
        issues.append("common_regions_below_50000")
    if metadata_overlap is None or metadata_overlap < 0.95:
        issues.append("metadata_overlap_below_95_percent")
    if age_coverage is None or age_coverage < 0.95:
        issues.append("age_coverage_below_95_percent")
    if beta_min is not None and beta_min < 0:
        issues.append("beta_min_below_0")
    if beta_max is not None and beta_max > 1:
        issues.append("beta_max_above_1")
    payload = {
        "status": "passed" if not issues else "failed",
        "matrix_gate_pass": not issues,
        "issues": issues,
        "metrics": {
            "manifest_status": manifest.get("status"),
            "common_regions": common_regions,
            "metadata_overlap": metadata_overlap,
            "age_coverage": age_coverage,
            "beta_min": beta_min,
            "beta_max": beta_max,
        },
    }
    write_json(out_dir / "matrix_gate_report.json", payload)
    return payload


def number_from_manifest(payload: dict[str, Any], keys: list[str]) -> float | None:
    stack = [payload]
    while stack:
        item = stack.pop()
        if isinstance(item, dict):
            for key in keys:
                if key in item:
                    try:
                        return float(item[key])
                    except (TypeError, ValueError):
                        return None
            stack.extend(value for value in item.values() if isinstance(value, dict))
    return None


def learn_readiness(matrix: dict[str, Any], out_dir: Path) -> dict[str, Any]:
    if matrix.get("status") != "passed":
        payload = {
            "status": "blocked",
            "reason": "matrix_gate_not_passed",
            "ralph_learn_ready": False,
            "training_authorized": False,
            "autoresearch_authorized": False,
        }
    else:
        payload = {
            "status": "ready_pending_explicit_training_approval",
            "reason": "metadata_adapter_matrix_gates_passed",
            "ralph_learn_ready": True,
            "training_authorized": False,
            "autoresearch_authorized": False,
            "next_fixed_benchmarks": [
                "GroupKFold by dataset_batch",
                "leave-one-dataset-out",
                "all_except:GSE121141 -> GSE121141",
                "all_except:GSE80672 -> GSE80672",
                "random-label sanity",
                "shuffled CR sanity when CR metrics are computed",
            ],
        }
    write_json(out_dir / "ralph_learn_readiness.json", payload)
    return payload


def write_sop() -> None:
    SOP_PATH.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Route A Gate Workflow SOP",
        "",
        "Date: 2026-05-19",
        "",
        "This SOP controls how generated or collaborative Route A data enters the project.",
        "",
        "## State Machine",
        "",
        "1. Metadata gate: run `validate_route_a_submission.py` on the sample sheet and file manifest.",
        "2. Adapter smoke: only after metadata passes, inspect 1-3 local processed methylation files or the first 1,000-10,000 rows.",
        "3. Matrix gate: only after adapter smoke passes, check local path existence, manifest `file_size_bytes`, manifest `sha256`, common 5kb regions, metadata overlap, age coverage, sex/MT exclusion, and beta range.",
        "4. RALPH Learn readiness: only after matrix gate passes, prepare fixed benchmark commands. Training still requires explicit approval.",
        "",
        "## Default Command",
        "",
        "Build a local file manifest from a filled sample sheet after files are staged:",
        "",
        "```bash",
        "VIRTUAL_ENV=.venv-core uv run --active python scripts/validate/build_route_a_local_file_manifest.py \\",
        "  --sample-sheet <filled_sample_sheet.csv> \\",
        "  --submission-id <submission_id> \\",
        "  --output results/route_a_local_file_manifest/<submission_id>/route_a_file_manifest.csv",
        "```",
        "",
        "Preferred one-command runner for a real filled submission. If `--file-manifest` is omitted, the runner builds a local manifest from the sample sheet first:",
        "",
        "```bash",
        "VIRTUAL_ENV=.venv-core uv run --active python scripts/validate/run_route_a_intake_to_learn_ready.py \\",
        "  --sample-sheet <filled_sample_sheet.csv> \\",
        "  --submission-id <submission_id>",
        "```",
        "",
        "Pass `--file-manifest <filled_file_manifest.csv>` only when using a pre-built manifest. This runner executes the same gate order below and stops before training. It only prepares a guarded RALPH Learn command package when all gates pass.",
        "",
        "Metadata and adapter gate:",
        "",
        "```bash",
        "VIRTUAL_ENV=.venv-core uv run --active python scripts/validate/run_route_a_gate_workflow.py \\",
        "  --sample-sheet <filled_sample_sheet.csv> \\",
        "  --file-manifest <filled_file_manifest.csv>",
        "```",
        "",
        "Build matrix manifest after adapter smoke passes:",
        "",
        "```bash",
        "VIRTUAL_ENV=.venv-core uv run --active python scripts/etl/22_build_route_a_matrix.py \\",
        "  --sample-sheet <filled_sample_sheet.csv> \\",
        "  --file-manifest <filled_file_manifest.csv> \\",
        "  --submission-id <submission_id> \\",
        "  --reference-matrix results/multidataset_v8_2_prefilter_liftover/all_rrbs_region_matrix_5kb.parquet \\",
        "  --out-dir results/route_a_intake/<submission_id>",
        "```",
        "",
        "Matrix gate and RALPH Learn readiness:",
        "",
        "```bash",
        "VIRTUAL_ENV=.venv-core uv run --active python scripts/validate/run_route_a_gate_workflow.py \\",
        "  --sample-sheet <filled_sample_sheet.csv> \\",
        "  --file-manifest <filled_file_manifest.csv> \\",
        "  --matrix-manifest results/route_a_intake/<submission_id>/route_a_matrix_manifest.json \\",
        "  --output-dir results/route_a_intake/<submission_id>",
        "```",
        "",
        "Prepare a guarded RALPH Learn command package after readiness:",
        "",
        "```bash",
        "VIRTUAL_ENV=.venv-core uv run --active python scripts/validate/run_route_a_ralph_learn_package.py \\",
        "  --intake-dir results/route_a_intake/<submission_id> \\",
        "  --base-matrix results/multidataset_v8_2_prefilter_liftover/all_rrbs_region_matrix_5kb.parquet \\",
        "  --base-manifest results/multidataset_v8_2_prefilter_liftover/all_rrbs_matrix_manifest.json \\",
        "  --out-dir results/route_a_ralph_learn_package/<submission_id>",
        "```",
        "",
        "The generated shell script must remain guarded and must not run training until explicit approval is given.",
        "",
        "Use `--allow-placeholders` only for template dry-runs.",
        "",
        "Contract regression test:",
        "",
        "```bash",
        "VIRTUAL_ENV=.venv-core uv run --active python scripts/validate/test_route_a_gate_contracts.py",
        "```",
        "",
        "This test reruns a synthetic gate smoke, verifies the local manifest generator and runner auto-manifest mode, verifies the matrix thresholds, checks that remote paths, bad beta values, and checksum mismatches are blocked, and confirms the generated command script exits before any training command.",
        "",
        "## Guardrails",
        "",
        "- No download is authorized by this workflow.",
        "- Remote URIs are recorded as pending; files must be staged separately before adapter smoke.",
        "- No Bismark/FASTQ ETL is authorized.",
        "- No training or autoresearch is authorized.",
        "- If any gate fails, stop and fix data/metadata before moving forward.",
    ]
    SOP_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_dry_report(state: dict[str, Any]) -> None:
    DRY_REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Route A Gate Workflow Dry-Run Report",
        "",
        f"Date: {datetime.now().isoformat(timespec='seconds')}",
        "",
        "## Summary",
        "",
        "The Route A gate workflow was executed on the 72-sample template in placeholder mode.",
        "This validates the state machine and metadata validator plumbing only.",
        "",
        "## Decision",
        "",
        f"- workflow decision: `{state['decision']}`",
        f"- metadata status: `{state['metadata_gate_status']}`",
        f"- adapter status: `{state['adapter_smoke_status']}`",
        f"- matrix status: `{state['matrix_gate_status']}`",
        f"- RALPH Learn status: `{state['ralph_learn_status']}`",
        "",
        "No download, Bismark, training, or autoresearch was run.",
        "",
        "## Outputs",
        "",
        "- `results/route_a_gate_workflow/route_a_gate_state.json`",
        "- `results/route_a_gate_workflow/metadata_gate_report.json`",
        "- `results/route_a_gate_workflow/adapter_smoke_report.json`",
        "- `results/route_a_gate_workflow/matrix_gate_report.json`",
        "- `results/route_a_gate_workflow/ralph_learn_readiness.json`",
    ]
    DRY_REPORT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_intake_report(state: dict[str, Any], metadata: dict[str, Any], adapter: dict[str, Any], matrix: dict[str, Any], learn: dict[str, Any], out_dir: Path) -> None:
    INTAKE_REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Route A Intake Gate Readiness Report",
        "",
        f"Date: {datetime.now().isoformat(timespec='seconds')}",
        "",
        "## Summary",
        "",
        "This report records the latest Route A intake gate workflow state.",
        "It does not authorize download, Bismark, model training, or autoresearch.",
        "",
        "## Decision",
        "",
        f"- workflow decision: `{state['decision']}`",
        f"- metadata gate: `{state['metadata_gate_status']}`",
        f"- adapter smoke: `{state['adapter_smoke_status']}`",
        f"- matrix gate: `{state['matrix_gate_status']}`",
        f"- RALPH Learn readiness: `{state['ralph_learn_status']}`",
        "",
        "## Gate Details",
        "",
        f"- metadata errors: `{metadata.get('error_count', 'N/A')}`; warnings: `{metadata.get('warning_count', 'N/A')}`",
        f"- adapter reason: `{adapter.get('reason', '')}`",
        f"- adapter files passed: `{adapter.get('passed_files', 0)}` / `{adapter.get('local_checked_files', 0)}`",
        f"- matrix reason/issues: `{matrix.get('reason', '') or ';'.join(matrix.get('issues', []))}`",
        f"- matrix metrics: `{matrix.get('metrics', {})}`",
        f"- learn reason: `{learn.get('reason', '')}`",
        "",
        "## Outputs",
        "",
        f"- gate output directory: `{out_dir}`",
        "- `metadata_gate_report.json`",
        "- `adapter_smoke_report.json`",
        "- `matrix_gate_report.json`",
        "- `ralph_learn_readiness.json`",
        "- `route_a_gate_state.json`",
        "",
        "## Guardrails",
        "",
        f"- download_authorized: `{state.get('download_authorized')}`",
        f"- bismark_authorized: `{state.get('bismark_authorized')}`",
        f"- training_authorized: `{state.get('training_authorized')}`",
        f"- autoresearch_authorized: `{state.get('autoresearch_authorized')}`",
    ]
    INTAKE_REPORT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def update_index() -> None:
    if not INDEX_PATH.exists():
        return
    text = INDEX_PATH.read_text(encoding="utf-8")
    entries = [
        "- Route A gate workflow SOP: `doc/30_protocols/08_20260519_route_a_gate_workflow_sop.md`",
        "- Route A gate workflow dry-run: `doc/20_analysis/37_20260519_route_a_gate_workflow_dry_run_report.md`",
        "- Route A intake gate readiness report: `doc/20_analysis/38_20260519_route_a_intake_gate_readiness_report.md`",
    ]
    marker = "- Route A partner handoff package: `doc/30_protocols/07_20260519_route_a_partner_handoff_package.md`"
    for entry in entries:
        if entry in text:
            continue
        if marker in text:
            text = text.replace(marker, marker + "\n" + entry)
        else:
            text += "\n" + entry + "\n"
    INDEX_PATH.write_text(text, encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-sheet", required=True, type=Path)
    parser.add_argument("--file-manifest", type=Path)
    parser.add_argument("--matrix-manifest", type=Path)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--allow-placeholders", action="store_true")
    parser.add_argument("--max-smoke-files", type=int, default=3)
    parser.add_argument("--max-smoke-rows", type=int, default=10000)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    out_dir = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    metadata = run_metadata_gate(args, out_dir)
    adapter = adapter_smoke(args, metadata, out_dir)
    matrix = matrix_gate(args, adapter, out_dir)
    learn = learn_readiness(matrix, out_dir)
    if metadata.get("status") != "passed":
        decision = "blocked_metadata_gate_failed"
    elif adapter.get("status") != "passed":
        decision = "metadata_passed_waiting_for_adapter_smoke"
    elif matrix.get("status") != "passed":
        decision = "adapter_passed_waiting_for_matrix_gate"
    elif learn.get("ralph_learn_ready"):
        decision = "ready_for_ralph_learn_pending_explicit_training_approval"
    else:
        decision = "blocked_unknown_state"
    state = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "decision": decision,
        "metadata_gate_status": metadata.get("status"),
        "adapter_smoke_status": adapter.get("status"),
        "matrix_gate_status": matrix.get("status"),
        "ralph_learn_status": learn.get("status"),
        "sample_sheet": str(args.sample_sheet),
        "file_manifest": str(args.file_manifest) if args.file_manifest else "",
        "matrix_manifest": str(args.matrix_manifest) if args.matrix_manifest else "",
        "download_authorized": False,
        "bismark_authorized": False,
        "training_authorized": False,
        "autoresearch_authorized": False,
    }
    write_json(out_dir / "route_a_gate_state.json", state)
    write_sop()
    write_dry_report(state)
    write_intake_report(state, metadata, adapter, matrix, learn, out_dir)
    update_index()
    print(json.dumps(state, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

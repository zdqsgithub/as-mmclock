#!/usr/bin/env python3
"""Run Route A intake through RALPH Learn readiness without training.

This orchestrator is for real generated/collaborative Route A submissions. It
wraps the existing gated scripts in the only allowed order:

1. metadata gate + adapter smoke
2. matrix build
3. matrix gate + readiness decision
4. guarded RALPH Learn command package

It never downloads data, runs Bismark, trains models, or starts autoresearch.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Any

ROOT = Path("/home/zdq-as/mouse_methyl_work")
DEFAULT_REFERENCE = ROOT / "results" / "multidataset_v8_2_prefilter_liftover" / "all_rrbs_region_matrix_5kb.parquet"
DEFAULT_BASE_MATRIX = DEFAULT_REFERENCE
DEFAULT_BASE_MANIFEST = ROOT / "results" / "multidataset_v8_2_prefilter_liftover" / "all_rrbs_matrix_manifest.json"
DEFAULT_INTAKE_BASE = ROOT / "results" / "route_a_intake"
DEFAULT_PACKAGE_BASE = ROOT / "results" / "route_a_ralph_learn_package"
REPORT_PATH = ROOT / "doc" / "20_analysis" / "40_20260519_route_a_intake_to_learn_ready_runner_report.md"
INDEX_PATH = ROOT / "doc" / "00_meta" / "02_20260519_project_status_index_v13.md"

GATE_SCRIPT = ROOT / "scripts" / "validate" / "run_route_a_gate_workflow.py"
MANIFEST_BUILDER_SCRIPT = ROOT / "scripts" / "validate" / "build_route_a_local_file_manifest.py"
MATRIX_SCRIPT = ROOT / "scripts" / "etl" / "22_build_route_a_matrix.py"
PACKAGE_SCRIPT = ROOT / "scripts" / "validate" / "run_route_a_ralph_learn_package.py"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str), encoding="utf-8")


def run_stage(name: str, cmd: list[str], out_dir: Path, stages: list[dict[str, Any]]) -> subprocess.CompletedProcess[str]:
    started = time.time()
    result = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True, check=False)
    stage = {
        "stage": name,
        "returncode": result.returncode,
        "command": " ".join(cmd),
        "exec_time_sec": round(time.time() - started, 3),
        "stdout_log": str(out_dir / f"{name}.stdout.log"),
        "stderr_log": str(out_dir / f"{name}.stderr.log"),
    }
    (out_dir / f"{name}.stdout.log").write_text(result.stdout, encoding="utf-8")
    (out_dir / f"{name}.stderr.log").write_text(result.stderr, encoding="utf-8")
    stages.append(stage)
    return result


def summarize_gate(intake_dir: Path) -> dict[str, Any]:
    state_path = intake_dir / "route_a_gate_state.json"
    if not state_path.exists():
        return {"decision": "missing_route_a_gate_state"}
    state = read_json(state_path)
    return {
        "decision": state.get("decision"),
        "metadata_gate_status": state.get("metadata_gate_status"),
        "adapter_smoke_status": state.get("adapter_smoke_status"),
        "matrix_gate_status": state.get("matrix_gate_status"),
        "ralph_learn_status": state.get("ralph_learn_status"),
    }


def summarize_matrix(intake_dir: Path) -> dict[str, Any]:
    manifest_path = intake_dir / "route_a_matrix_manifest.json"
    blocker_path = intake_dir / "route_a_matrix_blocker.json"
    if manifest_path.exists():
        manifest = read_json(manifest_path)
        return {
            "status": manifest.get("status"),
            "common_regions": manifest.get("common_regions"),
            "metadata_overlap": manifest.get("metadata_overlap"),
            "age_coverage": manifest.get("age_coverage"),
            "beta_min": manifest.get("beta_min"),
            "beta_max": manifest.get("beta_max"),
            "blockers": manifest.get("blockers", []),
        }
    if blocker_path.exists():
        return read_json(blocker_path)
    return {"status": "missing_matrix_manifest"}


def summarize_package(package_dir: Path) -> dict[str, Any]:
    state_path = package_dir / "route_a_ralph_learn_package_state.json"
    if not state_path.exists():
        return {"status": "not_prepared"}
    state = read_json(state_path)
    return {
        "status": state.get("status"),
        "combined_matrix_path": state.get("combined_matrix_path"),
        "model_metadata_path": state.get("model_metadata_path"),
        "command_manifest_csv": state.get("command_manifest_csv"),
        "training_authorized": state.get("training_authorized"),
        "download_authorized": state.get("download_authorized"),
        "bismark_authorized": state.get("bismark_authorized"),
        "autoresearch_authorized": state.get("autoresearch_authorized"),
    }


def write_report(state: dict[str, Any]) -> None:
    lines = [
        "# Route A Intake To Learn Ready Runner Report",
        "",
        f"Date: {state['timestamp']}",
        "",
        "## Summary",
        "",
        "This runner executes the Route A gate sequence up to RALPH Learn readiness.",
        "It does not train, download, run Bismark, or start autoresearch.",
        "",
        "## State",
        "",
        f"- submission_id: `{state['submission_id']}`",
        f"- status: `{state['status']}`",
        f"- intake_dir: `{state['intake_dir']}`",
        f"- package_dir: `{state.get('package_dir', '')}`",
        f"- gate decision: `{state.get('gate_summary', {}).get('decision')}`",
        f"- matrix status: `{state.get('matrix_summary', {}).get('status')}`",
        f"- package status: `{state.get('package_summary', {}).get('status')}`",
        "",
        "## Guardrails",
        "",
        f"- training_authorized: `{state['training_authorized']}`",
        f"- download_authorized: `{state['download_authorized']}`",
        f"- bismark_authorized: `{state['bismark_authorized']}`",
        f"- autoresearch_authorized: `{state['autoresearch_authorized']}`",
        "",
        "## Stage Return Codes",
        "",
    ]
    for stage in state.get("stages", []):
        lines.append(f"- `{stage['stage']}`: returncode `{stage['returncode']}`")
    lines.extend(
        [
            "",
            "## Next Step",
            "",
            "If status is `ready_for_ralph_learn_pending_explicit_training_approval`,",
            "review the guarded command manifest and request explicit training approval",
            "before running any benchmark command.",
        ]
    )
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def update_index() -> None:
    entry = "- Route A intake-to-learn-ready runner report: `doc/20_analysis/40_20260519_route_a_intake_to_learn_ready_runner_report.md`"
    if not INDEX_PATH.exists():
        return
    text = INDEX_PATH.read_text(encoding="utf-8")
    if entry in text:
        return
    marker = "- Route A RALPH Learn package report: `doc/20_analysis/39_20260519_route_a_ralph_learn_package_report.md`"
    if marker in text:
        text = text.replace(marker, marker + "\n" + entry)
    else:
        text += "\n" + entry + "\n"
    INDEX_PATH.write_text(text, encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-sheet", required=True, type=Path)
    parser.add_argument("--file-manifest", type=Path)
    parser.add_argument("--submission-id", required=True)
    parser.add_argument("--reference-matrix", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--base-matrix", type=Path, default=DEFAULT_BASE_MATRIX)
    parser.add_argument("--base-manifest", type=Path, default=DEFAULT_BASE_MANIFEST)
    parser.add_argument("--intake-dir", type=Path, default=None)
    parser.add_argument("--package-dir", type=Path, default=None)
    parser.add_argument("--max-smoke-files", type=int, default=3)
    parser.add_argument("--max-smoke-rows", type=int, default=10000)
    parser.add_argument("--build-file-manifest", action="store_true", help="Build a local file manifest from the sample sheet before gating.")
    parser.add_argument("--file-manifest-output", type=Path, default=None, help="Output path when building a local file manifest.")
    parser.add_argument("--skip-package", action="store_true", help="Stop after readiness and do not prepare the guarded command package.")
    parser.add_argument("--skip-report-update", action="store_true", help="Do not update the shared report or project index; useful for contract tests.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    started = time.time()
    intake_dir = args.intake_dir or DEFAULT_INTAKE_BASE / args.submission_id
    package_dir = args.package_dir or DEFAULT_PACKAGE_BASE / args.submission_id
    intake_dir.mkdir(parents=True, exist_ok=True)
    stages: list[dict[str, Any]] = []

    if args.build_file_manifest or args.file_manifest is None:
        manifest_out = args.file_manifest_output or ROOT / "results" / "route_a_local_file_manifest" / args.submission_id / "route_a_file_manifest.csv"
        manifest_report = manifest_out.parent / "route_a_file_manifest_build_report.json"
        manifest_cmd = [
            sys.executable,
            str(MANIFEST_BUILDER_SCRIPT),
            "--sample-sheet",
            str(args.sample_sheet),
            "--submission-id",
            str(args.submission_id),
            "--output",
            str(manifest_out),
            "--report",
            str(manifest_report),
        ]
        manifest_stage = run_stage("00_build_local_file_manifest", manifest_cmd, intake_dir, stages)
        if manifest_stage.returncode != 0:
            status = "blocked_file_manifest_build_failed"
            args.file_manifest = manifest_out
            state = build_state(args, intake_dir, package_dir, stages, status, started, summarize_gate(intake_dir), summarize_matrix(intake_dir), summarize_package(package_dir))
            finish(state, exit_code=2)
        args.file_manifest = manifest_out

    pre_gate_cmd = [
        sys.executable,
        str(GATE_SCRIPT),
        "--sample-sheet",
        str(args.sample_sheet),
        "--file-manifest",
        str(args.file_manifest),
        "--output-dir",
        str(intake_dir),
        "--max-smoke-files",
        str(args.max_smoke_files),
        "--max-smoke-rows",
        str(args.max_smoke_rows),
    ]
    pre_gate = run_stage("01_metadata_adapter_gate", pre_gate_cmd, intake_dir, stages)
    gate_summary = summarize_gate(intake_dir)
    if pre_gate.returncode != 0 or gate_summary.get("decision") not in {
        "adapter_passed_waiting_for_matrix_gate",
        "ready_for_ralph_learn_pending_explicit_training_approval",
    }:
        status = gate_summary.get("decision", "blocked_pre_matrix_gate")
        state = build_state(args, intake_dir, package_dir, stages, status, started, gate_summary, summarize_matrix(intake_dir), summarize_package(package_dir))
        finish(state, exit_code=2)

    matrix_cmd = [
        sys.executable,
        str(MATRIX_SCRIPT),
        "--sample-sheet",
        str(args.sample_sheet),
        "--file-manifest",
        str(args.file_manifest),
        "--submission-id",
        str(args.submission_id),
        "--reference-matrix",
        str(args.reference_matrix),
        "--out-dir",
        str(intake_dir),
    ]
    matrix = run_stage("02_build_matrix", matrix_cmd, intake_dir, stages)
    matrix_summary = summarize_matrix(intake_dir)
    if matrix.returncode != 0 or matrix_summary.get("status") != "completed":
        status = "blocked_matrix_gate"
        state = build_state(args, intake_dir, package_dir, stages, status, started, gate_summary, matrix_summary, summarize_package(package_dir))
        finish(state, exit_code=2)

    final_gate_cmd = pre_gate_cmd + ["--matrix-manifest", str(intake_dir / "route_a_matrix_manifest.json")]
    final_gate = run_stage("03_matrix_gate_readiness", final_gate_cmd, intake_dir, stages)
    gate_summary = summarize_gate(intake_dir)
    if final_gate.returncode != 0 or gate_summary.get("decision") != "ready_for_ralph_learn_pending_explicit_training_approval":
        status = gate_summary.get("decision", "blocked_readiness_gate")
        state = build_state(args, intake_dir, package_dir, stages, status, started, gate_summary, matrix_summary, summarize_package(package_dir))
        finish(state, exit_code=2)

    package_summary = summarize_package(package_dir)
    if not args.skip_package:
        package_cmd = [
            sys.executable,
            str(PACKAGE_SCRIPT),
            "--intake-dir",
            str(intake_dir),
            "--submission-id",
            str(args.submission_id),
            "--base-matrix",
            str(args.base_matrix),
            "--base-manifest",
            str(args.base_manifest),
            "--out-dir",
            str(package_dir),
        ]
        package = run_stage("04_prepare_guarded_learn_package", package_cmd, intake_dir, stages)
        package_summary = summarize_package(package_dir)
        if package.returncode != 0 or package_summary.get("status") != "prepared_pending_explicit_training_approval":
            status = "blocked_learn_package_preparation"
            state = build_state(args, intake_dir, package_dir, stages, status, started, gate_summary, matrix_summary, package_summary)
            finish(state, exit_code=2)

    state = build_state(
        args,
        intake_dir,
        package_dir,
        stages,
        "ready_for_ralph_learn_pending_explicit_training_approval",
        started,
        gate_summary,
        matrix_summary,
        package_summary,
    )
    finish(state, exit_code=0)


def build_state(
    args: argparse.Namespace,
    intake_dir: Path,
    package_dir: Path,
    stages: list[dict[str, Any]],
    status: str,
    started: float,
    gate_summary: dict[str, Any],
    matrix_summary: dict[str, Any],
    package_summary: dict[str, Any],
) -> dict[str, Any]:
    return {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "submission_id": args.submission_id,
        "status": status,
        "sample_sheet": str(args.sample_sheet),
        "file_manifest": str(args.file_manifest),
        "intake_dir": str(intake_dir),
        "package_dir": str(package_dir),
        "reference_matrix": str(args.reference_matrix),
        "base_matrix": str(args.base_matrix),
        "base_manifest": str(args.base_manifest),
        "gate_summary": gate_summary,
        "matrix_summary": matrix_summary,
        "package_summary": package_summary,
        "stages": stages,
        "training_authorized": False,
        "download_authorized": False,
        "bismark_authorized": False,
        "autoresearch_authorized": False,
        "report_update_enabled": not args.skip_report_update,
        "file_manifest_built_by_runner": any(stage["stage"] == "00_build_local_file_manifest" for stage in stages),
        "exec_time_sec": round(time.time() - started, 3),
    }


def finish(state: dict[str, Any], exit_code: int) -> None:
    intake_dir = Path(state["intake_dir"])
    write_json(intake_dir / "route_a_intake_to_learn_ready_state.json", state)
    if state.get("report_update_enabled", True):
        write_report(state)
        update_index()
    print(json.dumps(state, indent=2, sort_keys=True, default=str))
    raise SystemExit(exit_code)


if __name__ == "__main__":
    main()

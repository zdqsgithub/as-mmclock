#!/usr/bin/env python3
"""Route A gate/readiness regression contracts.

These tests exercise only Route A intake gates and command-package guardrails.
They do not train models, download data, run FASTQ/Bismark ETL, or start
autoresearch.
"""
from __future__ import annotations

import io
import csv
import json
import subprocess
import sys
import unittest
from datetime import datetime
from pathlib import Path
from typing import Any

ROOT = Path("/home/zdq-as/mouse_methyl_work")
OUT_DIR = ROOT / "results" / "route_a_contract_tests"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


class RouteAGateContracts(unittest.TestCase):
    def test_intake_runner_reaches_ready_without_package_or_training(self) -> None:
        intake_dir = OUT_DIR / "synthetic_contract_intake"
        package_dir = OUT_DIR / "synthetic_contract_package"
        completed = subprocess.run(
            [
                sys.executable,
                "scripts/validate/run_route_a_intake_to_learn_ready.py",
                "--sample-sheet",
                "results/route_a_intake/synthetic_ready/sample_sheet.csv",
                "--file-manifest",
                "results/route_a_intake/synthetic_ready/file_manifest.csv",
                "--submission-id",
                "synthetic_contract",
                "--intake-dir",
                str(intake_dir),
                "--package-dir",
                str(package_dir),
                "--skip-package",
                "--skip-report-update",
            ],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        state = load_json(intake_dir / "route_a_intake_to_learn_ready_state.json")
        self.assertEqual(state.get("status"), "ready_for_ralph_learn_pending_explicit_training_approval")
        self.assertEqual(state.get("gate_summary", {}).get("metadata_gate_status"), "passed")
        self.assertEqual(state.get("gate_summary", {}).get("adapter_smoke_status"), "passed")
        self.assertEqual(state.get("gate_summary", {}).get("matrix_gate_status"), "passed")
        self.assertFalse(state.get("training_authorized"))
        self.assertFalse(state.get("download_authorized"))
        self.assertFalse(state.get("bismark_authorized"))
        self.assertFalse(state.get("autoresearch_authorized"))

    def test_intake_runner_can_auto_build_local_file_manifest(self) -> None:
        intake_dir = OUT_DIR / "auto_manifest_intake"
        manifest_path = OUT_DIR / "auto_manifest" / "route_a_file_manifest.csv"
        completed = subprocess.run(
            [
                sys.executable,
                "scripts/validate/run_route_a_intake_to_learn_ready.py",
                "--sample-sheet",
                "results/route_a_intake/synthetic_ready/sample_sheet.csv",
                "--submission-id",
                "synthetic_auto_manifest",
                "--intake-dir",
                str(intake_dir),
                "--package-dir",
                str(OUT_DIR / "auto_manifest_package"),
                "--file-manifest-output",
                str(manifest_path),
                "--skip-package",
                "--skip-report-update",
            ],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertTrue(manifest_path.exists())
        state = load_json(intake_dir / "route_a_intake_to_learn_ready_state.json")
        self.assertEqual(state.get("status"), "ready_for_ralph_learn_pending_explicit_training_approval")
        self.assertTrue(state.get("file_manifest_built_by_runner"))
        self.assertFalse(state.get("training_authorized"))
        self.assertFalse(state.get("download_authorized"))

    def test_synthetic_ready_matrix_gate_thresholds(self) -> None:
        manifest = load_json(ROOT / "results" / "route_a_intake" / "synthetic_ready" / "route_a_matrix_manifest.json")
        self.assertEqual(manifest.get("status"), "completed")
        self.assertGreaterEqual(int(manifest.get("common_regions", 0)), 50000)
        self.assertGreaterEqual(float(manifest.get("metadata_overlap", 0)), 0.95)
        self.assertGreaterEqual(float(manifest.get("age_coverage", 0)), 0.95)
        self.assertGreaterEqual(float(manifest.get("beta_min", -1)), 0.0)
        self.assertLessEqual(float(manifest.get("beta_max", 2)), 1.0)
        self.assertTrue(manifest.get("file_integrity_checked"))
        self.assertIn("sha256", str(manifest.get("file_integrity_checks", "")))
        self.assertEqual(manifest.get("blockers"), [])
        self.assertFalse(manifest.get("training_authorized"))
        self.assertFalse(manifest.get("download_authorized"))
        self.assertFalse(manifest.get("bismark_authorized"))
        self.assertFalse(manifest.get("autoresearch_authorized"))

    def test_guarded_command_manifest_cannot_execute_accidentally(self) -> None:
        script_path = ROOT / "results" / "route_a_ralph_learn_package" / "synthetic_ready" / "ralph_learn_commands.sh"
        text = script_path.read_text(encoding="utf-8")
        exit_pos = text.find("exit 2")
        train_pos = text.find("train_clock.py")
        self.assertNotEqual(exit_pos, -1)
        self.assertNotEqual(train_pos, -1)
        self.assertLess(exit_pos, train_pos)
        state = load_json(ROOT / "results" / "route_a_ralph_learn_package" / "synthetic_ready" / "route_a_ralph_learn_package_state.json")
        self.assertEqual(state.get("status"), "prepared_pending_explicit_training_approval")
        self.assertFalse(state.get("training_authorized"))
        self.assertFalse(state.get("download_authorized"))
        self.assertFalse(state.get("bismark_authorized"))
        self.assertFalse(state.get("autoresearch_authorized"))

    def test_bad_beta_fixture_is_blocked_at_adapter_smoke(self) -> None:
        out_dir = OUT_DIR / "bad_beta_gate"
        completed = subprocess.run(
            [
                sys.executable,
                "scripts/validate/run_route_a_gate_workflow.py",
                "--sample-sheet",
                "results/route_a_intake/synthetic_schema_smoke/sample_sheet.csv",
                "--file-manifest",
                "results/route_a_intake/synthetic_schema_smoke/file_manifest_bad_beta.csv",
                "--output-dir",
                str(out_dir),
            ],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        state = load_json(out_dir / "route_a_gate_state.json")
        self.assertNotEqual(state.get("decision"), "ready_for_ralph_learn_pending_explicit_training_approval")
        self.assertIn(state.get("adapter_smoke_status"), {"failed", "blocked"})
        self.assertFalse(state.get("training_authorized"))
        self.assertFalse(state.get("download_authorized"))

    def test_file_checksum_mismatch_blocks_matrix_builder(self) -> None:
        out_dir = OUT_DIR / "bad_checksum_matrix"
        manifest_path = OUT_DIR / "file_manifest_bad_checksum.csv"
        source_manifest = ROOT / "results" / "route_a_intake" / "synthetic_ready" / "file_manifest.csv"
        with source_manifest.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
            fieldnames = list(rows[0].keys())
        rows[0]["sha256"] = "0" * 64
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        with manifest_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

        completed = subprocess.run(
            [
                sys.executable,
                "scripts/etl/22_build_route_a_matrix.py",
                "--sample-sheet",
                "results/route_a_intake/synthetic_ready/sample_sheet.csv",
                "--file-manifest",
                str(manifest_path),
                "--submission-id",
                "synthetic_bad_checksum",
                "--out-dir",
                str(out_dir),
            ],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 2, completed.stdout)
        blocker = load_json(out_dir / "route_a_matrix_blocker.json")
        self.assertEqual(blocker.get("reason"), "local_processed_path_gate_failed")
        issues = blocker.get("path_issues", [])
        self.assertTrue(any(issue.get("issue") == "file_sha256_mismatch" for issue in issues))

    def test_local_file_manifest_builder_outputs_usable_manifest(self) -> None:
        manifest_path = OUT_DIR / "generated_from_sample_sheet_file_manifest.csv"
        report_path = OUT_DIR / "generated_from_sample_sheet_file_manifest_report.json"
        intake_dir = OUT_DIR / "generated_manifest_intake"
        package_dir = OUT_DIR / "generated_manifest_package"
        completed = subprocess.run(
            [
                sys.executable,
                "scripts/validate/build_route_a_local_file_manifest.py",
                "--sample-sheet",
                "results/route_a_intake/synthetic_ready/sample_sheet.csv",
                "--submission-id",
                "synthetic_generated_manifest",
                "--output",
                str(manifest_path),
                "--report",
                str(report_path),
            ],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        report = load_json(report_path)
        self.assertEqual(report.get("status"), "passed")
        self.assertEqual(int(report.get("n_manifest_rows", 0)), 72)
        self.assertEqual(int(report.get("n_unique_files", 0)), 1)

        ready = subprocess.run(
            [
                sys.executable,
                "scripts/validate/run_route_a_intake_to_learn_ready.py",
                "--sample-sheet",
                "results/route_a_intake/synthetic_ready/sample_sheet.csv",
                "--file-manifest",
                str(manifest_path),
                "--submission-id",
                "synthetic_generated_manifest",
                "--intake-dir",
                str(intake_dir),
                "--package-dir",
                str(package_dir),
                "--skip-package",
                "--skip-report-update",
            ],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(ready.returncode, 0, ready.stderr)
        state = load_json(intake_dir / "route_a_intake_to_learn_ready_state.json")
        self.assertEqual(state.get("status"), "ready_for_ralph_learn_pending_explicit_training_approval")
        self.assertFalse(state.get("training_authorized"))

    def test_local_file_manifest_builder_blocks_remote_uri(self) -> None:
        source = ROOT / "results" / "route_a_intake" / "synthetic_ready" / "sample_sheet.csv"
        bad_sheet = OUT_DIR / "sample_sheet_remote_path.csv"
        with source.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
            fieldnames = list(rows[0].keys())
        rows[0]["processed_coverage_path"] = "https://example.org/not_local.cov.gz"
        bad_sheet.parent.mkdir(parents=True, exist_ok=True)
        with bad_sheet.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        completed = subprocess.run(
            [
                sys.executable,
                "scripts/validate/build_route_a_local_file_manifest.py",
                "--sample-sheet",
                str(bad_sheet),
                "--submission-id",
                "synthetic_remote_blocked",
                "--output",
                str(OUT_DIR / "remote_blocked_manifest.csv"),
                "--report",
                str(OUT_DIR / "remote_blocked_manifest_report.json"),
            ],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 2, completed.stdout)
        report = load_json(OUT_DIR / "remote_blocked_manifest_report.json")
        self.assertEqual(report.get("status"), "failed")
        self.assertTrue(any(error.get("code") == "remote_uri_not_supported" for error in report.get("errors", [])))
        self.assertFalse(report.get("training_authorized"))


def run_suite() -> tuple[unittest.result.TestResult, str]:
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(RouteAGateContracts)
    stream = io.StringIO()
    runner = unittest.TextTestRunner(stream=stream, verbosity=2)
    result = runner.run(suite)
    return result, stream.getvalue()


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    result, text = run_suite()
    payload = {
        "created_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "status": "passed" if result.wasSuccessful() else "failed",
        "tests_run": result.testsRun,
        "failures": len(result.failures),
        "errors": len(result.errors),
        "failure_details": [detail for _, detail in result.failures],
        "error_details": [detail for _, detail in result.errors],
        "training_authorized": False,
        "download_authorized": False,
        "bismark_authorized": False,
        "autoresearch_authorized": False,
    }
    (OUT_DIR / "test_report.txt").write_text(text, encoding="utf-8")
    (OUT_DIR / "test_report.json").write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    print(text)
    print(json.dumps(payload, indent=2, sort_keys=True))
    sys.exit(0 if result.wasSuccessful() else 1)


if __name__ == "__main__":
    main()

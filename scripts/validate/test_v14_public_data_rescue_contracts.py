#!/usr/bin/env python3
"""Regression contracts for v14 public-data rescue outputs.

The tests only read v14 outputs. They do not query networks, download files,
run Bismark, train models, or start autoresearch.
"""
from __future__ import annotations

import csv
import io
import json
import sys
import unittest
from datetime import datetime
from pathlib import Path
from typing import Any

ROOT = Path("/home/zdq-as/mouse_methyl_work")
OUT_DIR = ROOT / "results" / "ralph_v14_public_data_rescue"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


class V14PublicDataRescueContracts(unittest.TestCase):
    def test_guardrails_are_disabled(self) -> None:
        state = read_json(OUT_DIR / "ralph_decision_state.json")
        self.assertFalse(state.get("training_authorized"))
        self.assertFalse(state.get("download_authorized"))
        self.assertFalse(state.get("bismark_authorized"))
        self.assertFalse(state.get("autoresearch_authorized"))
        self.assertFalse(state.get("human_clock_mapping_authorized"))
        self.assertFalse(state.get("dummy_auc_authorized"))

    def test_required_outputs_exist(self) -> None:
        for name in [
            "candidate_refresh_table.csv",
            "candidate_gate_table.csv",
            "pilot_run_manifest.csv",
            "candidate_smoke_manifest.csv",
            "network_resolution_log.jsonl",
            "v14_public_data_rescue_report.md",
        ]:
            self.assertTrue((OUT_DIR / name).exists(), name)

    def test_no_headline_candidate_is_promoted(self) -> None:
        rows = read_rows(OUT_DIR / "candidate_gate_table.csv")
        self.assertGreater(len(rows), 0)
        promoted = [row for row in rows if row.get("headline_allowed") == "True"]
        self.assertEqual(promoted, [])

    def test_gse83947_is_processed_adapter_pilot_only(self) -> None:
        rows = read_rows(OUT_DIR / "candidate_gate_table.csv")
        gse83947 = [row for row in rows if row.get("dataset") == "GSE83947"]
        self.assertEqual(len(gse83947), 1)
        self.assertEqual(gse83947[0].get("tier"), "P3_raw_pilot")
        self.assertIn("processed_adapter_smoke", gse83947[0].get("gate_status", ""))
        pilots = [row for row in read_rows(OUT_DIR / "pilot_run_manifest.csv") if row.get("dataset") == "GSE83947"]
        self.assertGreaterEqual(len(pilots), 1)
        self.assertTrue(any(row.get("processed_supplement_url") for row in pilots))
        self.assertTrue(all(row.get("download_authorized") == "False" for row in pilots))

    def test_existing_auxiliary_evidence_is_not_promoted(self) -> None:
        rows = {row["dataset"]: row for row in read_rows(OUT_DIR / "candidate_gate_table.csv")}
        self.assertEqual(rows["GSE286302"]["tier"], "P2_auxiliary")
        self.assertEqual(rows["GSE281602"]["tier"], "P2_auxiliary")
        self.assertEqual(rows["GSE213628"]["tier"], "P2_auxiliary")

    def test_gse83947_matrix_gate_does_not_authorize_learn(self) -> None:
        state_path = OUT_DIR / "gse83947_region_matrix_gate" / "GSE83947_matrix_gate_state.json"
        self.assertTrue(state_path.exists())
        state = read_json(state_path)
        self.assertEqual(state.get("status"), "completed")
        self.assertFalse(state.get("common_region_gate_50000"))
        self.assertLess(int(state.get("common_regions_with_v8_2", 0)), 50000)
        self.assertFalse(state.get("training_authorized"))
        self.assertFalse(state.get("bismark_authorized"))
        self.assertFalse(state.get("autoresearch_authorized"))
        self.assertFalse(state.get("headline_allowed"))

    def test_gse83947_route_b_preflight_is_ready_after_environment_gate(self) -> None:
        state_path = OUT_DIR / "gse83947_raw_pilot_preflight" / "route_b_gse83947_raw_pilot_preflight_state.json"
        self.assertTrue(state_path.exists())
        state = read_json(state_path)
        self.assertTrue(state.get("metadata_gate"))
        self.assertTrue(state.get("tool_gate"))
        self.assertEqual(state.get("decision"), "pending_explicit_download_and_compute_approval")
        self.assertFalse(state.get("raw_fastq_download_authorized"))
        self.assertFalse(state.get("bismark_authorized"))
        self.assertFalse(state.get("training_authorized"))
        self.assertGreaterEqual(int(state.get("total_fastq_bytes", 0)), 1_000_000_000)

    def test_route_b_environment_and_fastq_download_are_traceable(self) -> None:
        env_state = read_json(OUT_DIR / "route_b_environment" / "route_b_bismark_environment_manifest.json")
        self.assertEqual(env_state.get("status"), "ready")
        self.assertTrue(env_state.get("reference", {}).get("bismark_index_ready"))
        self.assertIn("GRCm38", env_state.get("reference", {}).get("assembly", ""))
        self.assertFalse(env_state.get("authorizations", {}).get("training_authorized"))
        self.assertFalse(env_state.get("authorizations", {}).get("autoresearch_authorized"))

        download_state = read_json(OUT_DIR / "gse83947_raw_fastq_pilot" / "raw_fastq_download_state.json")
        self.assertEqual(download_state.get("status"), "completed")
        self.assertTrue(download_state.get("raw_fastq_download_authorized"))
        self.assertTrue(download_state.get("all_size_match"))
        self.assertTrue(download_state.get("all_md5_match"))
        self.assertEqual(int(download_state.get("n_samples", 0)), 3)
        self.assertFalse(download_state.get("bismark_authorized"))
        self.assertFalse(download_state.get("training_authorized"))
        self.assertFalse(download_state.get("autoresearch_authorized"))

    def test_route_b_bismark_pilot_fails_matrix_gate_without_training(self) -> None:
        state_path = OUT_DIR / "gse83947_bismark_pilot" / "gse83947_bismark_pilot_state.json"
        self.assertTrue(state_path.exists())
        state = read_json(state_path)
        self.assertEqual(state.get("status"), "completed")
        self.assertEqual(int(state.get("n_samples", 0)), 3)
        self.assertTrue(state.get("bismark_authorized"))
        self.assertTrue(state.get("beta_range_valid"))
        self.assertFalse(state.get("common_region_gate_50000"))
        self.assertLess(int(state.get("common_regions_with_v8_2", 0)), 50000)
        self.assertFalse(state.get("headline_allowed"))
        self.assertEqual(state.get("decision"), "matrix_gate_failed_auxiliary_only")
        self.assertFalse(state.get("training_authorized"))
        self.assertFalse(state.get("autoresearch_authorized"))


def run_suite() -> tuple[unittest.result.TestResult, str]:
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(V14PublicDataRescueContracts)
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream, verbosity=2).run(suite)
    return result, stream.getvalue()


def main() -> None:
    result, text = run_suite()
    payload = {
        "created_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "status": "passed" if result.wasSuccessful() else "failed",
        "tests_run": result.testsRun,
        "failures": len(result.failures),
        "errors": len(result.errors),
        "failure_details": [detail for _, detail in result.failures],
        "error_details": [detail for _, detail in result.errors],
    }
    (OUT_DIR / "test_report.txt").write_text(text, encoding="utf-8")
    (OUT_DIR / "test_report.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(text)
    print(json.dumps(payload, indent=2))
    sys.exit(0 if result.wasSuccessful() else 1)


if __name__ == "__main__":
    main()

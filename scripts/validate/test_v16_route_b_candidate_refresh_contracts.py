#!/usr/bin/env python3
"""Contract tests for the v16 Route B candidate refresh controller."""
from __future__ import annotations

import csv
import json
import unittest
from pathlib import Path


ROOT = Path("/home/zdq-as/mouse_methyl_work")
OUT_DIR = ROOT / "results" / "ralph_v16_route_b_candidate_refresh"


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def as_bool(value: str) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


class V16RouteBCandidateRefreshContracts(unittest.TestCase):
    def test_required_outputs_exist(self) -> None:
        for name in [
            "candidate_refresh_table.csv",
            "candidate_gate_table.csv",
            "pilot_priority_queue.csv",
            "pilot_run_manifest.csv",
            "matrix_gate_attempts.csv",
            "three_strike_table.csv",
            "network_resolution_log.jsonl",
            "ralph_decision_state.json",
        ]:
            self.assertTrue((OUT_DIR / name).exists(), name)

    def test_decision_state_guardrails(self) -> None:
        state = read_json(OUT_DIR / "ralph_decision_state.json")
        self.assertEqual(state["loop_version"], "v16_route_b_candidate_refresh")
        self.assertIn(
            state["decision"],
            {
                "ready_for_ralph_learn_pending_explicit_training_approval",
                "pilot_candidates_prioritized_pending_explicit_single_candidate_authorization",
                "no_new_route_b_pilot_candidate_keep_route_a_or_redefine_search",
            },
        )
        self.assertFalse(state["download_authorized"])
        self.assertFalse(state["bismark_authorized"])
        self.assertFalse(state["training_authorized"])
        self.assertFalse(state["autoresearch_authorized"])
        self.assertEqual(state["metrics"]["common_region_gate"], 50000)

    def test_gse83947_is_explicitly_not_repiloted(self) -> None:
        gate = read_csv(OUT_DIR / "candidate_gate_table.csv")
        rows = [row for row in gate if row["dataset"] == "GSE83947"]
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["tier"], "P2_auxiliary")
        self.assertEqual(rows[0]["gate_status"], "do_not_repilot_headline")
        self.assertIn("prior_failed_candidate", rows[0]["blocker_type"])
        queue = read_csv(OUT_DIR / "pilot_priority_queue.csv")
        self.assertNotIn("GSE83947", {row["dataset"] for row in queue})

    def test_priority_queue_only_contains_p1_or_p3_and_no_authorization(self) -> None:
        queue = read_csv(OUT_DIR / "pilot_priority_queue.csv")
        for row in queue:
            self.assertIn(row["tier"], {"P1_processed_headline", "P3_raw_pilot"})
            self.assertEqual(row["pilot_authorization_state"], "pending_explicit_authorize_pilot")
            self.assertFalse(as_bool(row["download_authorized"]))
            self.assertFalse(as_bool(row["bismark_authorized"]))
            self.assertFalse(as_bool(row["training_authorized"]))
            self.assertFalse(as_bool(row["autoresearch_authorized"]))

    def test_blocked_candidates_are_not_headline_allowed(self) -> None:
        gate = read_csv(OUT_DIR / "candidate_gate_table.csv")
        for row in gate:
            if row["tier"] in {"P2_auxiliary", "P4_blocked"}:
                self.assertFalse(as_bool(row["headline_allowed"]))
            if row["tier"] == "P4_blocked":
                self.assertNotIn(row["dataset"], {item["dataset"] for item in read_csv(OUT_DIR / "pilot_priority_queue.csv")})

    def test_network_log_records_official_sources_and_no_fastq_download(self) -> None:
        lines = [
            json.loads(line)
            for line in (OUT_DIR / "network_resolution_log.jsonl").read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
        self.assertGreater(len(lines), 0)
        docs = {line.get("source_doc", "") for line in lines}
        self.assertTrue(any("ncbi.nlm.nih.gov/geo/info/soft" in doc for doc in docs))
        self.assertTrue(
            any("sradownload" in doc for doc in docs)
            or any("ena-docs.readthedocs.io" in doc for doc in docs)
        )
        for line in lines:
            self.assertNotEqual(line.get("resolution"), "fastq_downloaded")


if __name__ == "__main__":
    unittest.main(verbosity=2)

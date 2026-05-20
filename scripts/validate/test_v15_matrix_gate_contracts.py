#!/usr/bin/env python3
"""Contract tests for the v15 matrix-gate RALPH controller."""
from __future__ import annotations

import csv
import json
import unittest
from pathlib import Path


ROOT = Path("/home/zdq-as/mouse_methyl_work")
OUT_DIR = ROOT / "results" / "ralph_v15_matrix_gate"


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


class V15MatrixGateContracts(unittest.TestCase):
    def test_decision_state_guardrails(self) -> None:
        state = read_json(OUT_DIR / "ralph_decision_state.json")
        self.assertEqual(state["decision"], "three_strike_re_evaluate_no_current_headline_matrix_gate_pass")
        self.assertTrue(state["three_strike_triggered"])
        self.assertEqual(state["ready_candidates"], [])
        self.assertFalse(state["training_authorized"])
        self.assertFalse(state["download_authorized"])
        self.assertFalse(state["bismark_authorized"])
        self.assertFalse(state["autoresearch_authorized"])

    def test_success_definition_is_matrix_gate_specific(self) -> None:
        state = read_json(OUT_DIR / "ralph_decision_state.json")
        success = state["success_definition"]
        self.assertEqual(success["common_5kb_regions_min"], 50000)
        self.assertEqual(success["sample_specific_exact_age_coverage_min"], 0.95)
        self.assertIn("brain_cortex", success["target_tissues"])
        self.assertIn("heart", success["target_tissues"])
        self.assertIn("lung", success["target_tissues"])

    def test_auxiliary_matrix_is_not_headline(self) -> None:
        state = read_json(OUT_DIR / "ralph_decision_state.json")
        self.assertIn("gse286302_processed_cov_conversion", state["auxiliary_matrix_candidates"])
        self.assertIn("gse224442_processed_cov_conversion", state["auxiliary_matrix_candidates"])
        rows = read_csv(OUT_DIR / "matrix_gate_attempts.csv")
        for row in rows:
            if row["attempt_id"] in state["auxiliary_matrix_candidates"]:
                self.assertEqual(row["matrix_gate_pass"], "True")
                self.assertEqual(row["headline_matrix_gate_pass"], "False")

    def test_three_strike_table_blocks_training(self) -> None:
        rows = read_csv(OUT_DIR / "three_strike_table.csv")
        self.assertEqual(len(rows), 3)
        datasets = {row["dataset"] for row in rows}
        self.assertIn("GSE83947", datasets)
        for row in rows:
            self.assertNotEqual(row["headline_matrix_gate_pass"], "True")

    def test_network_resolution_log_uses_official_docs(self) -> None:
        log_path = OUT_DIR / "network_resolution_log.jsonl"
        lines = [json.loads(line) for line in log_path.read_text(encoding="utf-8").splitlines() if line.strip()]
        docs = {line["source_doc"] for line in lines}
        self.assertIn("https://www.ncbi.nlm.nih.gov/geo/info/download.html", docs)
        self.assertIn("https://www.ncbi.nlm.nih.gov/geo/info/geo_paccess.html", docs)
        self.assertIn("https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/", docs)
        self.assertIn("https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/file-reports.html", docs)


if __name__ == "__main__":
    unittest.main(verbosity=2)


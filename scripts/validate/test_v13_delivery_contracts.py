#!/usr/bin/env python3
"""Smoke/regression contracts for the v13 delivery freeze.

The tests only read v12.1/v13 outputs, rerun the v13 decision script, and write
test reports under results/v13_delivery_freeze. They do not train models,
download data, run FASTQ/Bismark, rebuild matrices, or start autoresearch.
"""
from __future__ import annotations

import csv
import importlib.util
import io
import json
import subprocess
import sys
import unittest
from datetime import datetime
from pathlib import Path
from typing import Any

ROOT = Path("/home/zdq-as/mouse_methyl_work")
OUT_DIR = ROOT / "results" / "v13_delivery_freeze"


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


class V13DeliveryContracts(unittest.TestCase):
    def test_core_environment_does_not_include_training_stack(self) -> None:
        self.assertIsNone(importlib.util.find_spec("torch"))
        self.assertIsNone(importlib.util.find_spec("lightning"))

    def test_v13_strategy_decision_reproduces_with_uv(self) -> None:
        completed = subprocess.run(
            ["uv", "run", "--active", "python", "scripts/validate/run_v13_strategy_decision.py"],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        payload = json.loads(completed.stdout)
        self.assertEqual(payload.get("selected_route"), "route_c_accept_redefined_benchmark")
        self.assertFalse(payload.get("raw_fastq_download_authorized"))
        self.assertFalse(payload.get("training_authorized"))
        self.assertFalse(payload.get("autoresearch_authorized"))

    def test_decision_state_guardrails(self) -> None:
        decision = load_json(ROOT / "results" / "ralph_v13_strategy" / "v13_route_decision_state.json")
        self.assertEqual(decision.get("selected_route"), "route_c_accept_redefined_benchmark")
        self.assertFalse(decision.get("raw_fastq_download_authorized"))
        self.assertFalse(decision.get("training_authorized"))
        self.assertFalse(decision.get("autoresearch_authorized"))

    def test_support_covered_rows_obey_definition(self) -> None:
        rows = read_csv_rows(ROOT / "results" / "benchmark_v12_1_redefined" / "prediction_support_annotations.csv")
        covered_count = 0
        for row in rows:
            if row.get("support_covered") != "True":
                continue
            covered_count += 1
            same_tissue_n = int(float(row["same_tissue_train_n"]))
            gap = float(row["age_support_gap_weeks"])
            self.assertGreaterEqual(same_tissue_n, 10)
            self.assertLessEqual(gap, 8.0)
        self.assertGreater(covered_count, 0)

    def test_route_table_selects_only_route_c(self) -> None:
        rows = read_csv_rows(ROOT / "results" / "ralph_v13_strategy" / "v13_route_table.csv")
        selected = [row for row in rows if row.get("selected_now") == "True"]
        self.assertEqual([row["route"] for row in selected], ["C_accept_redefined_benchmark"])

    def test_cr_metrics_are_real_not_dummy(self) -> None:
        cr_metrics = load_json(ROOT / "results" / "benchmark_v12_1_redefined" / "gse80672_cr_redefined_metrics.json")
        self.assertEqual(cr_metrics.get("status"), "computed")
        self.assertNotEqual(cr_metrics.get("status"), "dummy")
        self.assertNotEqual(cr_metrics.get("best_source"), "dummy")
        self.assertIsNotNone(cr_metrics.get("best_source"))
        self.assertIsNotNone(cr_metrics.get("best_cr_detection_auc"))
        self.assertNotEqual(float(cr_metrics["best_cr_detection_auc"]), 0.5)


def run_suite() -> tuple[unittest.result.TestResult, str]:
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(V13DeliveryContracts)
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
    }
    (OUT_DIR / "test_report.txt").write_text(text, encoding="utf-8")
    (OUT_DIR / "test_report.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(text)
    print(json.dumps(payload, indent=2))
    sys.exit(0 if result.wasSuccessful() else 1)


if __name__ == "__main__":
    main()

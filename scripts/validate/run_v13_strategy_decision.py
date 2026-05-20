#!/usr/bin/env python3
"""Create the v13 data-strategy decision package.

v13 is a strategy decision after v12/v12.1, not an experiment runner. This
script reads the v12.1 redefined benchmark state and writes a route decision,
model scope statement, and next-action table. It does not download data, run
FASTQ/Bismark, train models, or start autoresearch.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

ROOT = Path("/home/zdq-as/mouse_methyl_work")
DEFAULT_STATE = ROOT / "results" / "benchmark_v12_1_redefined" / "v12_1_redefined_benchmark_state.json"
DEFAULT_CR = ROOT / "results" / "benchmark_v12_1_redefined" / "gse80672_cr_redefined_metrics.json"
OUT_DIR = ROOT / "results" / "ralph_v13_strategy"
REPORT = OUT_DIR / "v13_data_strategy_decision_report.md"
SCOPE = OUT_DIR / "model_scope_statement.md"


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str), encoding="utf-8")


def metric(payload: dict[str, Any], section: str, key: str) -> Any:
    value = payload.get(section, {})
    return value.get(key) if isinstance(value, dict) else None


def decide_route(state: dict[str, Any]) -> tuple[str, str]:
    headline_mae = metric(state, "support_covered_headline_metrics", "mae_weeks")
    unsupported_mae = metric(state, "unsupported_stress_test_metrics", "mae_weeks")
    old104_mae = metric(state, "gse121141_old104_stress_metrics", "mae_weeks")
    if headline_mae is not None and unsupported_mae is not None and float(headline_mae) < float(unsupported_mae):
        return (
            "route_c_accept_redefined_benchmark",
            "support_covered_metrics_outperform_unsupported_stress_metrics",
        )
    if old104_mae is not None and float(old104_mae) > 65.386:
        return (
            "route_a_or_b_required_before_model_work",
            "support_covered_benchmark_not_sufficient_or_old104_stress_remains_high",
        )
    return ("manual_review_required", "metrics_do_not_match_v13_default_decision_rule")


def build_route_table(state: dict[str, Any]) -> list[dict[str, Any]]:
    selected_route, reason = decide_route(state)
    return [
        {
            "route": "A_generate_or_collaborate",
            "status": "open_not_started",
            "selected_now": selected_route == "route_a_or_b_required_before_model_work",
            "entry_gate": "Need true full-lifespan old bulk brain_cortex/heart/lung clock.",
            "next_action": "Design or source old target-tissue RRBS/WGBS cohort.",
        },
        {
            "route": "B_minimal_fastq_bismark_pilot",
            "status": "blocked_until_explicit_approval",
            "selected_now": False,
            "entry_gate": "BioSample/RunInfo proves sample-specific exact age, bulk target tissue, assay, and practical run size.",
            "next_action": "Open separate minimal ETL approval before any FASTQ download.",
        },
        {
            "route": "C_accept_redefined_benchmark",
            "status": "accepted" if selected_route == "route_c_accept_redefined_benchmark" else "manual_review",
            "selected_now": selected_route == "route_c_accept_redefined_benchmark",
            "entry_gate": reason,
            "next_action": "Use support-covered chronological benchmark as current model headline.",
        },
    ]


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def md_table(rows: list[dict[str, Any]]) -> str:
    if not rows:
        return "No rows."
    columns = list(rows[0].keys())
    header = "| " + " | ".join(columns) + " |"
    sep = "| " + " | ".join(["---"] * len(columns)) + " |"
    body = [
        "| " + " | ".join("" if row.get(col) is None else str(row.get(col)) for col in columns) + " |"
        for row in rows
    ]
    return "\n".join([header, sep, *body])


def write_scope_statement(path: Path, state: dict[str, Any], cr_payload: dict[str, Any]) -> None:
    headline = state.get("support_covered_headline_metrics", {})
    unsupported = state.get("unsupported_stress_test_metrics", {})
    old104 = state.get("gse121141_old104_stress_metrics", {})
    best_cr_auc = cr_payload.get("best_cr_detection_auc")
    best_source = cr_payload.get("best_source")
    lines = [
        "# Mouse Methylation Clock Model Scope Statement",
        "",
        "Date: 2026-05-19",
        "",
        "## Accepted Scope",
        "",
        "The current model benchmark is valid only for support-covered chronological-age prediction. A prediction is support-covered when the training data contains at least 10 same-tissue samples and the held-out age is no more than 8 weeks above the training same-tissue maximum age.",
        "",
        f"- Support-covered headline MAE: {headline.get('mae_weeks')} weeks",
        f"- Support-covered samples evaluated: {headline.get('n_samples')}",
        f"- Unsupported stress-test MAE: {unsupported.get('mae_weeks')} weeks",
        "",
        "## Required Stress-Test Reporting",
        "",
        "`GSE121141 old104+ brain_cortex/heart/lung` must continue to be reported as a stress-test/blocker metric, not as the headline pass/fail benchmark.",
        "",
        f"- GSE121141 old104+ stress MAE: {old104.get('mae_weeks')} weeks",
        f"- GSE121141 old104+ stress samples: {old104.get('n_samples')}",
        "",
        "## Biological-Age / CR Scope",
        "",
        "CR metrics are research-level validation only and must use real held-out predictions plus shuffled-intervention sanity when available.",
        "",
        f"- Best observed CR AUC in v12.1: {best_cr_auc}",
        f"- Best CR source: {best_source}",
        "",
        "## Out of Scope",
        "",
        "- Claims of full-lifespan old target-tissue generalization.",
        "- Claims that old104+ GSE121141 failure is solved.",
        "- Autoresearch on unsupported old104+ metrics.",
        "- Raw FASTQ/Bismark execution without a separate minimal ETL approval.",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_report(path: Path, state: dict[str, Any], cr_payload: dict[str, Any], route_table: list[dict[str, Any]]) -> None:
    selected = [row for row in route_table if bool(row.get("selected_now"))]
    selected_route = selected[0]["route"] if selected else "manual_review_required"
    headline = state.get("support_covered_headline_metrics", {})
    unsupported = state.get("unsupported_stress_test_metrics", {})
    old104 = state.get("gse121141_old104_stress_metrics", {})
    lines = [
        "# v13 Data Strategy Decision Report",
        "",
        "Date: 2026-05-19",
        "",
        "## Summary",
        "",
        "v13 formalizes the post-v12 decision. It does not download data, run FASTQ/Bismark, train models, or start autoresearch.",
        "",
        f"- Selected route: `{selected_route}`",
        f"- Support-covered headline MAE: `{headline.get('mae_weeks')}` weeks",
        f"- Unsupported stress-test MAE: `{unsupported.get('mae_weeks')}` weeks",
        f"- GSE121141 old104+ stress MAE: `{old104.get('mae_weeks')}` weeks",
        f"- Best CR AUC: `{cr_payload.get('best_cr_detection_auc')}`",
        "",
        "## Route Decision Table",
        "",
        md_table(route_table),
        "",
        "## Decision",
        "",
        "Accept Route C for the current project state: support-covered chronological benchmark is the model headline, and old target-tissue full-lifespan validation moves to a future data project.",
        "",
        "Route A remains the preferred path if the project requires a true full-lifespan target-tissue mouse clock. Route B remains blocked until a separate minimal ETL approval identifies a valid raw candidate.",
        "",
        "## Outputs",
        "",
        "- `results/ralph_v13_strategy/v13_route_decision_state.json`",
        "- `results/ralph_v13_strategy/v13_route_table.csv`",
        "- `results/ralph_v13_strategy/model_scope_statement.md`",
        "- `results/ralph_v13_strategy/v13_data_strategy_decision_report.md`",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--state", default=str(DEFAULT_STATE))
    parser.add_argument("--cr_metrics", default=str(DEFAULT_CR))
    parser.add_argument("--out_dir", default=str(OUT_DIR))
    args = parser.parse_args()

    state = load_json(Path(args.state))
    cr_payload = load_json(Path(args.cr_metrics))
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    route, reason = decide_route(state)
    route_table = build_route_table(state)
    write_csv(out_dir / "v13_route_table.csv", route_table)
    write_scope_statement(out_dir / "model_scope_statement.md", state, cr_payload)
    decision_state = {
        "status": "completed",
        "selected_route": route,
        "decision_reason": reason,
        "raw_fastq_download_authorized": False,
        "training_authorized": False,
        "autoresearch_authorized": False,
        "support_covered_headline_metrics": state.get("support_covered_headline_metrics", {}),
        "unsupported_stress_test_metrics": state.get("unsupported_stress_test_metrics", {}),
        "gse121141_old104_stress_metrics": state.get("gse121141_old104_stress_metrics", {}),
        "cr_metrics_summary": {
            "status": cr_payload.get("status"),
            "best_source": cr_payload.get("best_source"),
            "best_cr_detection_auc": cr_payload.get("best_cr_detection_auc"),
        },
        "outputs": {
            "route_table": str(out_dir / "v13_route_table.csv"),
            "model_scope_statement": str(out_dir / "model_scope_statement.md"),
            "decision_report": str(out_dir / "v13_data_strategy_decision_report.md"),
        },
    }
    write_json(out_dir / "v13_route_decision_state.json", decision_state)
    write_report(out_dir / "v13_data_strategy_decision_report.md", state, cr_payload, route_table)
    print(json.dumps(decision_state, indent=2, default=str))


if __name__ == "__main__":
    main()

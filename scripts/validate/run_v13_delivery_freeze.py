#!/usr/bin/env python3
"""Freeze the v13 delivery package for reproducibility.

This script only reads existing v12.1/v13 outputs and writes v13 delivery
manifests. It does not train models, download data, run FASTQ/Bismark, rebuild
matrices, or start autoresearch.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

ROOT = Path("/home/zdq-as/mouse_methyl_work")
AS_DS_OPS_PROGRAM = Path("/home/zdq-as/as-ds-ops/program.md")
OUT_DIR = ROOT / "results" / "v13_delivery_freeze"

ARTIFACTS: list[tuple[Path, str]] = [
    (ROOT / "requirements.txt", "environment_input"),
    (ROOT / "requirements-core.txt", "environment_input_core"),
    (ROOT / "requirements-viz.txt", "environment_input_viz"),
    (ROOT / "requirements-full.txt", "environment_input_full"),
    (ROOT / "scripts" / "validate" / "run_v12_1_redefined_benchmark.py", "reproduction_script"),
    (ROOT / "scripts" / "validate" / "run_v13_strategy_decision.py", "reproduction_script"),
    (ROOT / "scripts" / "validate" / "run_v13_delivery_freeze.py", "delivery_script"),
    (ROOT / "scripts" / "validate" / "test_v13_delivery_contracts.py", "delivery_test"),
    (ROOT / "results" / "benchmark_v12_1_redefined" / "v12_1_redefined_benchmark_state.json", "v12_1_state"),
    (ROOT / "results" / "benchmark_v12_1_redefined" / "support_covered_headline_metrics.json", "v12_1_metric"),
    (ROOT / "results" / "benchmark_v12_1_redefined" / "unsupported_stress_test_metrics.json", "v12_1_metric"),
    (ROOT / "results" / "benchmark_v12_1_redefined" / "gse121141_old104_stress_metrics.json", "v12_1_metric"),
    (ROOT / "results" / "benchmark_v12_1_redefined" / "gse80672_cr_redefined_metrics.json", "v12_1_metric"),
    (ROOT / "results" / "benchmark_v12_1_redefined" / "prediction_support_annotations.csv", "v12_1_annotation"),
    (ROOT / "results" / "benchmark_v12_1_redefined" / "redefined_benchmark_summary.csv", "v12_1_summary"),
    (ROOT / "results" / "ralph_v13_strategy" / "v13_route_decision_state.json", "v13_decision"),
    (ROOT / "results" / "ralph_v13_strategy" / "v13_route_table.csv", "v13_decision"),
    (ROOT / "results" / "ralph_v13_strategy" / "model_scope_statement.md", "v13_scope"),
    (ROOT / "results" / "ralph_v13_strategy" / "v13_data_strategy_decision_report.md", "v13_report"),
    (ROOT / "doc" / "20_analysis" / "31_20260519_v12_1_redefined_benchmark_report.md", "analysis_doc"),
    (ROOT / "doc" / "20_analysis" / "32_20260519_v13_data_strategy_rfc.md", "analysis_doc"),
    (ROOT / "doc" / "20_analysis" / "33_20260519_v13_delivery_freeze_report.md", "analysis_doc"),
    (ROOT / "doc" / "00_meta" / "01_20260519_v13_model_card.md", "model_card"),
    (ROOT / "doc" / "00_meta" / "02_20260519_project_status_index_v13.md", "project_index"),
    (
        ROOT / "doc" / "30_protocols" / "04_20260519_route_a_old_tissue_data_generation_rfc_template.md",
        "route_a_template",
    ),
    (
        ROOT / "doc" / "30_protocols" / "05_20260519_route_b_minimal_fastq_bismark_pilot_rfc_template.md",
        "route_b_template",
    ),
    (ROOT / "results" / "v13_delivery_freeze" / "test_report.json", "delivery_test_output"),
    (ROOT / "results" / "v13_delivery_freeze" / "test_report.txt", "delivery_test_output"),
    (AS_DS_OPS_PROGRAM, "project_norm"),
]


def run_command(command: list[str], cwd: Path = ROOT) -> dict[str, Any]:
    completed = subprocess.run(command, cwd=cwd, text=True, capture_output=True, check=False)
    return {
        "command": command,
        "returncode": completed.returncode,
        "stdout": completed.stdout.strip(),
        "stderr": completed.stderr.strip(),
    }


def sha256_file(path: Path) -> str | None:
    if not path.exists() or not path.is_file():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str), encoding="utf-8")


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def artifact_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path, role in ARTIFACTS:
        rows.append(
            {
                "path": rel(path),
                "role": role,
                "exists": path.exists(),
                "size_bytes": path.stat().st_size if path.exists() and path.is_file() else None,
                "sha256": sha256_file(path),
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def summarize_dependencies(freeze_output: str) -> dict[str, Any]:
    packages = [line for line in freeze_output.splitlines() if line.strip()]
    wanted = {}
    for line in packages:
        name = line.split("==", 1)[0].lower()
        if name in {"pandas", "numpy", "scipy", "scikit-learn", "pyarrow", "fastparquet", "torch", "lightning"}:
            wanted[name] = line
    training_stack = any(
        line.lower().startswith(("torch==", "lightning==", "pytorch-lightning==", "nvidia-", "triton=="))
        for line in packages
    )
    return {
        "package_count": len(packages),
        "key_packages": wanted,
        "training_stack_present": training_stack,
        "freeze": packages,
    }


def write_summary(
    path: Path,
    state: dict[str, Any],
    decision: dict[str, Any],
    cr_payload: dict[str, Any],
    test_report: dict[str, Any],
) -> None:
    headline = state.get("support_covered_headline_metrics", {})
    unsupported = state.get("unsupported_stress_test_metrics", {})
    old104 = state.get("gse121141_old104_stress_metrics", {})
    lines = [
        "# v13 Delivery Freeze Summary",
        "",
        f"Date: {datetime.now().astimezone().isoformat(timespec='seconds')}",
        "",
        "## Decision",
        "",
        "- Selected route: `route_c_accept_redefined_benchmark`.",
        "- Current model headline is support-covered chronological-age prediction.",
        "- `GSE121141 old104+ brain_cortex/heart/lung` remains a stress-test/blocker metric.",
        "- This freeze does not authorize training, downloads, FASTQ/Bismark, or autoresearch.",
        "",
        "## Metrics",
        "",
        f"- Support-covered headline MAE: `{headline.get('mae_weeks')}` weeks over `{headline.get('n_samples')}` rows.",
        f"- Unsupported stress-test MAE: `{unsupported.get('mae_weeks')}` weeks over `{unsupported.get('n_samples')}` rows.",
        f"- GSE121141 old104+ stress MAE: `{old104.get('mae_weeks')}` weeks over `{old104.get('n_samples')}` rows.",
        f"- Best observed CR AUC: `{cr_payload.get('best_cr_detection_auc')}`.",
        "",
        "## Reproduction Outputs",
        "",
        "- `results/v13_delivery_freeze/environment_manifest.json`",
        "- `results/v13_delivery_freeze/artifact_inventory.csv`",
        "- `results/v13_delivery_freeze/reproducibility_manifest.json`",
        "- `results/v13_delivery_freeze/v13_delivery_summary.md`",
        "- `results/v13_delivery_freeze/test_report.json`",
        "- `results/v13_delivery_freeze/test_report.txt`",
        "- `doc/00_meta/01_20260519_v13_model_card.md`",
        "- `doc/00_meta/02_20260519_project_status_index_v13.md`",
        "- `doc/30_protocols/04_20260519_route_a_old_tissue_data_generation_rfc_template.md`",
        "- `doc/30_protocols/05_20260519_route_b_minimal_fastq_bismark_pilot_rfc_template.md`",
        "",
        "## Dependency Mode",
        "",
        "- Default delivery reproduction uses `requirements-core.txt`.",
        "- `requirements-full.txt` is reserved for training/deep-learning/full historical environments.",
        "",
        "## Guardrails",
        "",
        f"- raw_fastq_download_authorized: `{decision.get('raw_fastq_download_authorized')}`",
        f"- training_authorized: `{decision.get('training_authorized')}`",
        f"- autoresearch_authorized: `{decision.get('autoresearch_authorized')}`",
        "",
        "## Delivery Tests",
        "",
        f"- status: `{test_report.get('status', 'not_run')}`",
        f"- tests_run: `{test_report.get('tests_run')}`",
        f"- failures: `{test_report.get('failures')}`",
        f"- errors: `{test_report.get('errors')}`",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out_dir", default=str(OUT_DIR))
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    uv_version = run_command(["uv", "--version"])
    freeze = run_command(["uv", "pip", "freeze"])
    rows = artifact_rows()
    missing = [row["path"] for row in rows if not row["exists"]]
    write_csv(out_dir / "artifact_inventory.csv", rows)

    state_path = ROOT / "results" / "benchmark_v12_1_redefined" / "v12_1_redefined_benchmark_state.json"
    decision_path = ROOT / "results" / "ralph_v13_strategy" / "v13_route_decision_state.json"
    cr_path = ROOT / "results" / "benchmark_v12_1_redefined" / "gse80672_cr_redefined_metrics.json"
    test_report_path = out_dir / "test_report.json"
    state = load_json(state_path)
    decision = load_json(decision_path)
    cr_payload = load_json(cr_path)
    test_report = load_json(test_report_path)

    env_manifest = {
        "created_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "cwd": str(ROOT),
        "uv_version": uv_version,
        "python_executable": sys.executable,
        "python_version": sys.version,
        "python_platform": platform.platform(),
        "dependency_snapshot": summarize_dependencies(freeze["stdout"]),
        "dependency_snapshot_returncode": freeze["returncode"],
        "guardrails": {
            "training_authorized": False,
            "download_authorized": False,
            "raw_fastq_download_authorized": False,
            "bismark_authorized": False,
            "autoresearch_authorized": False,
        },
    }
    write_json(out_dir / "environment_manifest.json", env_manifest)

    repro_manifest = {
        "status": "completed" if not missing else "completed_with_missing_artifacts",
        "missing_artifacts": missing,
        "commands": [
            "uv venv .venv-core",
            "VIRTUAL_ENV=.venv-core uv pip install -r requirements-core.txt",
            "VIRTUAL_ENV=.venv-core uv run --active python scripts/validate/run_v12_1_redefined_benchmark.py",
            "VIRTUAL_ENV=.venv-core uv run --active python scripts/validate/run_v13_strategy_decision.py",
            "VIRTUAL_ENV=.venv-core uv run --active python scripts/validate/run_v13_delivery_freeze.py",
            "VIRTUAL_ENV=.venv-core uv run --active python scripts/validate/test_v13_delivery_contracts.py",
        ],
        "dependency_files": {
            "core_delivery_default": "requirements-core.txt",
            "visualization_optional": "requirements-viz.txt",
            "full_training_or_historical": "requirements-full.txt",
            "compatibility_entrypoint": "requirements.txt",
        },
        "input_artifacts": rows,
        "decision": {
            "selected_route": decision.get("selected_route"),
            "decision_reason": decision.get("decision_reason"),
            "raw_fastq_download_authorized": decision.get("raw_fastq_download_authorized"),
            "training_authorized": decision.get("training_authorized"),
            "autoresearch_authorized": decision.get("autoresearch_authorized"),
        },
        "metrics": {
            "support_covered_headline_metrics": state.get("support_covered_headline_metrics", {}),
            "unsupported_stress_test_metrics": state.get("unsupported_stress_test_metrics", {}),
            "gse121141_old104_stress_metrics": state.get("gse121141_old104_stress_metrics", {}),
            "cr_metrics_summary": {
                "status": cr_payload.get("status"),
                "best_source": cr_payload.get("best_source"),
                "best_cr_detection_auc": cr_payload.get("best_cr_detection_auc"),
            },
        },
        "delivery_tests": test_report or {"status": "not_run"},
    }
    write_json(out_dir / "reproducibility_manifest.json", repro_manifest)
    write_summary(out_dir / "v13_delivery_summary.md", state, decision, cr_payload, test_report)
    print(json.dumps({"status": repro_manifest["status"], "out_dir": str(out_dir), "missing_artifacts": missing}, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Prepare Route A RALPH Learn fixed-benchmark package without training.

This script consumes a Route A intake directory that has reached
``ready_for_ralph_learn_pending_explicit_training_approval`` and writes the
combined matrix, model metadata, and command manifest needed for a later
explicitly-authorized benchmark run.
"""
from __future__ import annotations

import argparse
import csv
import json
import shlex
import time
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")
DEFAULT_BASE_MATRIX = ROOT / "results" / "multidataset_v8_2_prefilter_liftover" / "all_rrbs_region_matrix_5kb.parquet"
DEFAULT_BASE_MANIFEST = ROOT / "results" / "multidataset_v8_2_prefilter_liftover" / "all_rrbs_matrix_manifest.json"
DEFAULT_OUT_BASE = ROOT / "results" / "route_a_ralph_learn_package"
REPORT_PATH = ROOT / "doc" / "20_analysis" / "39_20260519_route_a_ralph_learn_package_report.md"
INDEX_PATH = ROOT / "doc" / "00_meta" / "02_20260519_project_status_index_v13.md"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str), encoding="utf-8")


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def require_ready(intake_dir: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    gate = read_json(intake_dir / "route_a_gate_state.json")
    manifest = read_json(intake_dir / "route_a_matrix_manifest.json")
    if gate.get("decision") != "ready_for_ralph_learn_pending_explicit_training_approval":
        raise SystemExit(f"Route A intake is not ready: {gate.get('decision')}")
    if manifest.get("status") != "completed":
        raise SystemExit(f"Route A matrix manifest is not completed: {manifest.get('status')}")
    return gate, manifest


def route_a_metadata(sample_sheet: Path, dataset_batch: str) -> pd.DataFrame:
    df = pd.read_csv(sample_sheet)
    df["dataset_batch"] = dataset_batch
    df["sample_id_source"] = "route_a_sample_sheet"
    df["metadata_source"] = "route_a_generated_or_collaborative_submission"
    df["assay"] = df.get("assay", "route_a_processed_methylation")
    for column in ["intervention", "tissue", "sex", "strain"]:
        if column not in df.columns:
            df[column] = ""
    return df


def build_model_metadata(base_metadata_path: Path, route_a_sample_sheet: Path, dataset_batch: str, out_path: Path) -> pd.DataFrame:
    base = pd.read_csv(base_metadata_path)
    route = route_a_metadata(route_a_sample_sheet, dataset_batch)
    all_columns = list(dict.fromkeys(list(base.columns) + list(route.columns)))
    base = base.reindex(columns=all_columns)
    route = route.reindex(columns=all_columns)
    combined = pd.concat([base, route], ignore_index=True)
    combined.to_csv(out_path, index=False)
    return combined


def build_combined_matrix(base_matrix_path: Path, route_a_matrix_path: Path, out_path: Path) -> dict[str, Any]:
    base = pd.read_parquet(base_matrix_path)
    route = pd.read_parquet(route_a_matrix_path)
    common = sorted(set(base.index.astype(str)) & set(route.index.astype(str)))
    if len(common) < 50_000:
        raise SystemExit(f"Combined Route A benchmark blocked: common regions {len(common)} < 50000")
    combined = pd.concat([base.loc[common], route.loc[common]], axis=1).astype("float32")
    combined.to_parquet(out_path, compression="zstd")
    return {
        "n_base_regions": int(base.shape[0]),
        "n_route_a_regions": int(route.shape[0]),
        "n_common_regions": int(len(common)),
        "n_base_samples": int(base.shape[1]),
        "n_route_a_samples": int(route.shape[1]),
        "n_combined_samples": int(combined.shape[1]),
        "matrix_path": str(out_path),
    }


def shell_join(parts: list[str]) -> str:
    return " ".join(shlex.quote(str(part)) for part in parts)


def train_command(
    *,
    matrix_path: Path,
    metadata_path: Path,
    out_dir: Path,
    model_type: str,
    preprocess: str,
    presence: float,
    top: str,
    target_transform: str,
    imputation: str = "median",
    agebin: float | None = None,
    shift: float | None = None,
    randomize: bool = False,
) -> str:
    cmd = [
        "VIRTUAL_ENV=.venv",
        "uv",
        "run",
        "--active",
        "python",
        "scripts/train/train_clock.py",
        "--matrix_path",
        str(matrix_path),
        "--feature_type",
        "region",
        "--metadata_path",
        str(metadata_path),
        "--model_type",
        model_type,
        "--preprocess",
        preprocess,
        "--imputation",
        imputation,
        "--min_train_feature_presence",
        str(presence),
        "--n_feature_prefilter",
        str(top),
        "--target_transform",
        target_transform,
        "--output_dir",
        str(out_dir),
    ]
    if agebin is not None:
        cmd.extend(["--min_train_agebin_feature_presence", str(agebin)])
    if shift is not None:
        cmd.extend(["--max_train_dataset_mean_shift", str(shift)])
    if randomize:
        cmd.append("--randomize_labels")
    return shell_join(cmd)


def benchmark_command(run_dir: Path) -> str:
    return shell_join(
        [
            "VIRTUAL_ENV=.venv-core",
            "uv",
            "run",
            "--active",
            "python",
            "scripts/validate/benchmark_metrics.py",
            "--predictions",
            str(run_dir / "predictions.csv"),
            "--out",
            str(run_dir / "benchmark_result.json"),
        ]
    )


def command_manifest(matrix_path: Path, metadata_path: Path, out_dir: Path) -> list[dict[str, Any]]:
    configs = [
        {
            "config_id": "lgbm_quantile_p08_top1000_agebin08",
            "model_type": "lgbm",
            "preprocess": "quantile_uniform",
            "presence": 0.8,
            "top": "1000",
            "target_transform": "log1p_days",
            "agebin": 0.8,
            "shift": None,
        },
        {
            "config_id": "lgbm_robust_p095_top1000_shift015",
            "model_type": "lgbm",
            "preprocess": "robust",
            "presence": 0.95,
            "top": "1000",
            "target_transform": "log1p_days",
            "agebin": None,
            "shift": 0.15,
        },
    ]
    rows: list[dict[str, Any]] = []
    for config in configs:
        run_dir = out_dir / config["config_id"]
        random_dir = out_dir / f"{config['config_id']}_random_label"
        train = train_command(
            matrix_path=matrix_path,
            metadata_path=metadata_path,
            out_dir=run_dir,
            model_type=config["model_type"],
            preprocess=config["preprocess"],
            presence=config["presence"],
            top=config["top"],
            target_transform=config["target_transform"],
            agebin=config["agebin"],
            shift=config["shift"],
        )
        random_train = train_command(
            matrix_path=matrix_path,
            metadata_path=metadata_path,
            out_dir=random_dir,
            model_type=config["model_type"],
            preprocess=config["preprocess"],
            presence=config["presence"],
            top=config["top"],
            target_transform=config["target_transform"],
            agebin=config["agebin"],
            shift=config["shift"],
            randomize=True,
        )
        rows.append(
            {
                **config,
                "run_type": "groupkfold_fixed_benchmark",
                "training_authorized": False,
                "command": train,
                "benchmark_command": benchmark_command(run_dir),
                "expected_predictions": str(run_dir / "predictions.csv"),
                "expected_result": str(run_dir / "result.json"),
                "expected_benchmark_result": str(run_dir / "benchmark_result.json"),
            }
        )
        rows.append(
            {
                **config,
                "run_type": "random_label_sanity",
                "training_authorized": False,
                "command": random_train,
                "benchmark_command": benchmark_command(random_dir),
                "expected_predictions": str(random_dir / "predictions.csv"),
                "expected_result": str(random_dir / "result.json"),
                "expected_benchmark_result": str(random_dir / "benchmark_result.json"),
            }
        )
    return rows


def write_shell_script(path: Path, rows: list[dict[str, Any]]) -> None:
    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "echo 'This script is a command manifest only. Review and explicitly authorize training before running.'",
        "exit 2",
        "",
    ]
    for row in rows:
        lines.extend(
            [
                f"# {row['run_type']} {row['config_id']}",
                row["command"],
                row["benchmark_command"],
                "",
            ]
        )
    path.write_text("\n".join(lines), encoding="utf-8")


def write_report(state: dict[str, Any], rows: list[dict[str, Any]]) -> None:
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Route A RALPH Learn Package Report",
        "",
        f"Date: {datetime.now().isoformat(timespec='seconds')}",
        "",
        "## Summary",
        "",
        "This package prepares fixed Route A benchmark inputs and commands after the intake gate reaches readiness.",
        "No training, download, Bismark, or autoresearch was run.",
        "",
        "## Package State",
        "",
        f"- status: `{state['status']}`",
        f"- submission_id: `{state['submission_id']}`",
        f"- combined matrix: `{state['combined_matrix_path']}`",
        f"- model metadata: `{state['model_metadata_path']}`",
        f"- command count: `{len(rows)}`",
        "",
        "## Commands",
        "",
        "Commands are written to `ralph_learn_command_manifest.csv/json` and `ralph_learn_commands.sh`.",
        "The shell script exits before executing commands until training is explicitly authorized.",
        "",
        "## Guardrails",
        "",
        f"- training_authorized: `{state['training_authorized']}`",
        f"- autoresearch_authorized: `{state['autoresearch_authorized']}`",
        f"- download_authorized: `{state['download_authorized']}`",
        f"- bismark_authorized: `{state['bismark_authorized']}`",
    ]
    REPORT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def update_index() -> None:
    if not INDEX_PATH.exists():
        return
    text = INDEX_PATH.read_text(encoding="utf-8")
    entry = "- Route A RALPH Learn package report: `doc/20_analysis/39_20260519_route_a_ralph_learn_package_report.md`"
    if entry in text:
        return
    marker = "- Route A intake gate readiness report: `doc/20_analysis/38_20260519_route_a_intake_gate_readiness_report.md`"
    if marker in text:
        text = text.replace(marker, marker + "\n" + entry)
    else:
        text += "\n" + entry + "\n"
    INDEX_PATH.write_text(text, encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--intake-dir", required=True, type=Path)
    parser.add_argument("--submission-id", default=None)
    parser.add_argument("--base-matrix", type=Path, default=DEFAULT_BASE_MATRIX)
    parser.add_argument("--base-manifest", type=Path, default=DEFAULT_BASE_MANIFEST)
    parser.add_argument("--out-dir", type=Path, default=None)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    started = time.time()
    intake_dir = args.intake_dir
    gate, route_manifest = require_ready(intake_dir)
    submission_id = args.submission_id or route_manifest.get("submission_id") or intake_dir.name
    out_dir = args.out_dir or DEFAULT_OUT_BASE / str(submission_id)
    out_dir.mkdir(parents=True, exist_ok=True)

    base_manifest = read_json(args.base_manifest)
    base_metadata_path = ROOT / base_manifest["metadata_path"]
    model_metadata_path = out_dir / "route_a_model_metadata.csv"
    model_meta = build_model_metadata(
        base_metadata_path=base_metadata_path,
        route_a_sample_sheet=ROOT / route_manifest["sample_sheet"],
        dataset_batch=route_manifest["dataset_batch"],
        out_path=model_metadata_path,
    )
    combined_matrix_path = out_dir / "route_a_plus_v8_2_region_matrix_5kb.parquet"
    matrix_metrics = build_combined_matrix(
        base_matrix_path=args.base_matrix,
        route_a_matrix_path=ROOT / route_manifest["region_matrix_path"],
        out_path=combined_matrix_path,
    )
    rows = command_manifest(combined_matrix_path, model_metadata_path, out_dir)
    write_csv(out_dir / "ralph_learn_command_manifest.csv", rows)
    write_json(out_dir / "ralph_learn_command_manifest.json", {"commands": rows})
    write_shell_script(out_dir / "ralph_learn_commands.sh", rows)
    state = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "status": "prepared_pending_explicit_training_approval",
        "submission_id": submission_id,
        "intake_dir": str(intake_dir),
        "gate_decision": gate.get("decision"),
        "combined_matrix_path": str(combined_matrix_path),
        "model_metadata_path": str(model_metadata_path),
        "base_matrix": str(args.base_matrix),
        "base_metadata_path": str(base_metadata_path),
        "route_a_matrix_manifest": str(intake_dir / "route_a_matrix_manifest.json"),
        "matrix_metrics": matrix_metrics,
        "metadata_rows": int(len(model_meta)),
        "route_a_dataset_batch": route_manifest["dataset_batch"],
        "command_manifest_csv": str(out_dir / "ralph_learn_command_manifest.csv"),
        "command_manifest_json": str(out_dir / "ralph_learn_command_manifest.json"),
        "command_shell_script": str(out_dir / "ralph_learn_commands.sh"),
        "training_authorized": False,
        "autoresearch_authorized": False,
        "download_authorized": False,
        "bismark_authorized": False,
        "exec_time_sec": round(time.time() - started, 3),
    }
    write_json(out_dir / "route_a_ralph_learn_package_state.json", state)
    write_report(state, rows)
    update_index()
    print(json.dumps(state, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

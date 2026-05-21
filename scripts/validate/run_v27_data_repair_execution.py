#!/usr/bin/env python3
"""Run the v27 lightweight data repair control checks.

This controller intentionally avoids starting heavy downloads or Bismark jobs.
It records whether the current v25 ETL is stable enough to allow a later raw
repair phase, and it aggregates outputs from the v27 lightweight repair audits.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")
DEFAULT_OUT = ROOT / "results" / "v27_data_repair_execution"
WATCHDOG_STATE = Path("/data/mouse_methyl/processed_v21_raw_etl/v21_raw_etl_watchdog_state.json")
SAMPLE_STATUS = Path("/data/mouse_methyl/processed_v21_raw_etl/v21_raw_etl_sample_status.csv")


def utc_now_dt() -> datetime:
    return datetime.now(timezone.utc)


def utc_now() -> str:
    return utc_now_dt().isoformat()


def read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str), encoding="utf-8")


def run_step(name: str, cmd: list[str], out_dir: Path) -> dict[str, Any]:
    started = utc_now()
    proc = subprocess.run(cmd, cwd=str(ROOT), text=True, capture_output=True, check=False)
    log_prefix = out_dir / "logs" / name
    log_prefix.parent.mkdir(parents=True, exist_ok=True)
    log_prefix.with_suffix(".stdout.log").write_text(proc.stdout, encoding="utf-8")
    log_prefix.with_suffix(".stderr.log").write_text(proc.stderr, encoding="utf-8")
    return {
        "name": name,
        "cmd": " ".join(cmd),
        "started_at": started,
        "finished_at": utc_now(),
        "returncode": proc.returncode,
        "stdout_log": str(log_prefix.with_suffix(".stdout.log")),
        "stderr_log": str(log_prefix.with_suffix(".stderr.log")),
    }


def process_alive(pid: Any) -> bool:
    try:
        pid_int = int(pid)
    except (TypeError, ValueError):
        return False
    return Path(f"/proc/{pid_int}").exists()


def resource_snapshot() -> dict[str, Any]:
    data_usage = shutil.disk_usage("/data")
    root_usage = shutil.disk_usage("/")
    load1, load5, load15 = os.getloadavg()
    return {
        "nproc": os.cpu_count(),
        "load1": round(load1, 3),
        "load5": round(load5, 3),
        "load15": round(load15, 3),
        "data_free_gib": round(data_usage.free / 1024**3, 3),
        "root_free_gib": round(root_usage.free / 1024**3, 3),
    }


def v25_stability() -> dict[str, Any]:
    state = read_json(WATCHDOG_STATE)
    snapshot = resource_snapshot()
    timestamp = state.get("timestamp")
    heartbeat_age_sec = None
    if timestamp:
        try:
            heartbeat_age_sec = (utc_now_dt() - datetime.fromisoformat(str(timestamp))).total_seconds()
        except ValueError:
            heartbeat_age_sec = None
    samples_failed = int(state.get("samples_failed") or 0)
    status_rows = int(state.get("status_rows") or 0)
    status_failed_rows = 0
    if SAMPLE_STATUS.exists():
        try:
            status_df = pd.read_csv(SAMPLE_STATUS)
            status_failed_rows = int(status_df["status"].astype(str).eq("failed").sum()) if "status" in status_df else 0
        except Exception:
            status_failed_rows = 0
    alive = process_alive(state.get("pid"))
    stable = (
        state.get("status") == "running"
        and alive
        and heartbeat_age_sec is not None
        and heartbeat_age_sec < 300
        and samples_failed == 0
        and status_failed_rows == 0
        and snapshot["data_free_gib"] > 10_000
        and snapshot["load15"] < (snapshot["nproc"] or 1) * 1.2
    )
    return {
        "watchdog_state_path": str(WATCHDOG_STATE),
        "state": state,
        "resources": snapshot,
        "heartbeat_age_sec": None if heartbeat_age_sec is None else round(float(heartbeat_age_sec), 3),
        "runner_alive": alive,
        "sample_status_rows": status_rows,
        "sample_status_failed_rows": status_failed_rows,
        "stable_for_heavy_download_or_raw_repair": stable,
        "stability_policy": "heavy raw repair waits for heartbeat<300s, runner alive, no failed rows, /data>10TiB free, load15<1.2*nproc",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--skip-steps", action="store_true", help="Only write the controller state.")
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    stability = v25_stability()
    write_json(args.out_dir / "v27_v25_stability_gate.json", stability)
    step_results: list[dict[str, Any]] = []
    if not args.skip_steps:
        py = str(ROOT / ".venv" / "bin" / "python")
        steps = [
            ("gse60012_mapping", [py, "scripts/etl/26_repair_gse60012_gsm_tile_mapping.py"]),
            ("gse93957_root_cause", [py, "scripts/validate/run_v27_gse93957_root_cause_audit.py"]),
            ("processed_candidate_manifest", [py, "scripts/etl/27_build_v27_processed_candidate_backlog.py"]),
            ("raw_missing_repair_queue", [py, "scripts/etl/28_build_v27_raw_missing_repair_queue.py"]),
            ("processed_pilot_manifest", [py, "scripts/etl/30_build_v27_processed_pilot_manifest.py", "--include-secondary-new"]),
            ("v11_3_hard_verify", [py, "scripts/etl/21_build_v11_3_conversion_queue.py"]),
            ("v11_3_schema_smoke", [py, "scripts/validate/v11_3_processed_schema_smoke.py"]),
        ]
        for name, cmd in steps:
            step_results.append(run_step(name, cmd, args.out_dir))
            pd.DataFrame(step_results).to_csv(args.out_dir / "v27_step_status.csv", index=False)

    gse60012 = read_json(args.out_dir / "gse60012_mapping" / "gse60012_mapping_summary.json")
    gse93957 = read_json(args.out_dir / "gse93957_root_cause" / "gse93957_root_cause_summary.json")
    processed = read_json(args.out_dir / "processed_candidates" / "v27_processed_candidate_manifest_summary.json")
    raw_repair = read_json(args.out_dir / "raw_repair" / "v27_raw_missing_repair_queue_summary.json")
    processed_pilot = read_json(args.out_dir / "processed_pilot" / "v27_processed_pilot_manifest_summary.json")
    v11_schema = read_json(ROOT / "results" / "download_backlog_v11_3" / "schema_smoke_summary.json")

    decision = {
        "timestamp": utc_now(),
        "status": "completed_with_failed_steps" if any(step.get("returncode") for step in step_results) else "completed",
        "v25_stable_for_heavy_raw_repair": stability["stable_for_heavy_download_or_raw_repair"],
        "raw_repair_action": "eligible_to_build_repair_queue" if stability["stable_for_heavy_download_or_raw_repair"] else "defer_until_v25_stable",
        "gse60012_mapping_status": gse60012.get("status", "not_run"),
        "gse93957_status": gse93957.get("status", "not_run"),
        "processed_candidate_manifest_status": processed.get("status", "not_run"),
        "raw_missing_repair_queue_status": raw_repair.get("status", "not_run"),
        "raw_missing_repair_download_allowed_runs": raw_repair.get("n_download_allowed", 0),
        "processed_pilot_manifest_status": processed_pilot.get("status", "not_run"),
        "processed_pilot_planned_download_rows": processed_pilot.get("n_planned_download_rows", 0),
        "v11_3_schema_smoke_status": v11_schema.get("status", "not_run"),
        "step_status_path": str(args.out_dir / "v27_step_status.csv"),
        "stability_gate_path": str(args.out_dir / "v27_v25_stability_gate.json"),
        "policy": "no_new_raw_download_or_bismark_started_by_this_controller",
    }
    write_json(args.out_dir / "v27_decision_state.json", decision)
    print(json.dumps(decision, indent=2, sort_keys=True, default=str))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Watchdog for the v21 raw Bismark ETL queue.

Default mode does not start the runner. Use both `--start-runner` and
`--authorize-bismark` after explicit pilot approval.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path("/home/zdq-as/mouse_methyl_work")
DEFAULT_QUEUE = ROOT / "results" / "ralph_v21_raid_raw_clock" / "raw_etl_queue.csv"
DEFAULT_OUT_ROOT = Path("/data/mouse_methyl/processed_v21_raw_etl")
RUNNER = ROOT / "scripts" / "etl" / "23_run_raw_bismark_queue.py"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str), encoding="utf-8")
    tmp.replace(path)


def append_log(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(payload, sort_keys=True, default=str) + "\n")


def recommend_parallel(requested: int) -> int:
    cpu_count = os.cpu_count() or 1
    try:
        load1, _, _ = os.getloadavg()
    except OSError:
        load1 = 0.0
    if load1 > cpu_count * 0.85 and requested > 2:
        return 2
    return max(1, requested)


def read_status(status_csv: Path) -> dict[str, Any]:
    if not status_csv.exists():
        return {"completed": 0, "failed": 0, "current_sample": None, "rows": 0}
    try:
        df = pd.read_csv(status_csv)
    except Exception:
        return {"completed": 0, "failed": 0, "current_sample": None, "rows": 0}
    completed = int(df.get("status", pd.Series(dtype=str)).eq("completed").sum())
    failed = int(df.get("status", pd.Series(dtype=str)).eq("failed").sum())
    current_sample = None
    if not df.empty and {"dataset", "sample_id"}.issubset(df.columns):
        last = df.iloc[-1]
        current_sample = f"{last.get('dataset')}/{last.get('sample_id')}"
    return {"completed": completed, "failed": failed, "current_sample": current_sample, "rows": int(len(df))}


def state_payload(
    *,
    status: str,
    pid: int | None,
    queue: Path,
    out_root: Path,
    status_csv: Path,
    requested_parallel: int,
    recommended_parallel: int,
    exit_code: int | None,
    started_at: str | None,
    message: str,
) -> dict[str, Any]:
    progress = read_status(status_csv)
    try:
        queue_rows = len(pd.read_csv(queue)) if queue.exists() else 0
    except Exception:
        queue_rows = 0
    return {
        "timestamp": utc_now(),
        "status": status,
        "pid": pid,
        "queue": str(queue),
        "out_root": str(out_root),
        "requested_parallel_samples": requested_parallel,
        "recommended_parallel_samples": recommended_parallel,
        "queue_rows": queue_rows,
        "samples_completed": progress["completed"],
        "samples_failed": progress["failed"],
        "status_rows": progress["rows"],
        "current_sample": progress["current_sample"],
        "retry_count": 0,
        "exit_code": exit_code,
        "started_at": started_at,
        "message": message,
        "training_authorized": False,
        "autoresearch_authorized": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--queue", type=Path, default=DEFAULT_QUEUE)
    parser.add_argument("--out-root", type=Path, default=DEFAULT_OUT_ROOT)
    parser.add_argument("--poll-seconds", type=int, default=60)
    parser.add_argument("--parallel-samples", type=int, default=3)
    parser.add_argument("--bismark-threads", type=int, default=8)
    parser.add_argument("--extract-threads", type=int, default=4)
    parser.add_argument("--max-samples", type=int, default=0)
    parser.add_argument("--start-runner", action="store_true")
    parser.add_argument("--authorize-bismark", action="store_true")
    args = parser.parse_args()

    args.out_root.mkdir(parents=True, exist_ok=True)
    state_path = args.out_root / "v21_raw_etl_watchdog_state.json"
    log_path = args.out_root / "v21_raw_etl_watchdog.log"
    status_csv = args.out_root / "v21_raw_etl_sample_status.csv"
    recommended = recommend_parallel(args.parallel_samples)

    if not args.start_runner or not args.authorize_bismark:
        status = "blocked"
        reason = "start_runner_and_authorize_bismark_required"
        payload = state_payload(
            status=status,
            pid=None,
            queue=args.queue,
            out_root=args.out_root,
            status_csv=status_csv,
            requested_parallel=args.parallel_samples,
            recommended_parallel=recommended,
            exit_code=2,
            started_at=None,
            message=reason,
        )
        write_json(state_path, payload)
        append_log(log_path, payload)
        print(json.dumps(payload, indent=2, sort_keys=True))
        raise SystemExit(2)

    cmd = [
        sys.executable,
        str(RUNNER),
        "--queue",
        str(args.queue),
        "--out-root",
        str(args.out_root),
        "--parallel-samples",
        str(recommended),
        "--bismark-threads",
        str(args.bismark_threads),
        "--extract-threads",
        str(args.extract_threads),
        "--authorize-bismark",
    ]
    if args.max_samples:
        cmd.extend(["--max-samples", str(args.max_samples)])
    started_at = utc_now()
    stdout_log = args.out_root / "v21_raw_etl_runner.stdout.log"
    stderr_log = args.out_root / "v21_raw_etl_runner.stderr.log"
    with stdout_log.open("a", encoding="utf-8") as out_handle, stderr_log.open("a", encoding="utf-8") as err_handle:
        proc = subprocess.Popen(cmd, cwd=str(ROOT), stdout=out_handle, stderr=err_handle, text=True)
        while proc.poll() is None:
            recommended_now = recommend_parallel(args.parallel_samples)
            payload = state_payload(
                status="running",
                pid=proc.pid,
                queue=args.queue,
                out_root=args.out_root,
                status_csv=status_csv,
                requested_parallel=args.parallel_samples,
                recommended_parallel=recommended_now,
                exit_code=None,
                started_at=started_at,
                message="runner_active",
            )
            write_json(state_path, payload)
            append_log(log_path, payload)
            time.sleep(max(5, args.poll_seconds))
        exit_code = proc.returncode

    final_status = "completed" if exit_code == 0 else "failed"
    payload = state_payload(
        status=final_status,
        pid=proc.pid,
        queue=args.queue,
        out_root=args.out_root,
        status_csv=status_csv,
        requested_parallel=args.parallel_samples,
        recommended_parallel=recommend_parallel(args.parallel_samples),
        exit_code=exit_code,
        started_at=started_at,
        message=f"runner_exit_{exit_code}",
    )
    write_json(state_path, payload)
    append_log(log_path, payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    raise SystemExit(exit_code)


if __name__ == "__main__":
    main()

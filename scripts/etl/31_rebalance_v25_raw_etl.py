#!/usr/bin/env python3
"""Gracefully rebalance a running v25 raw ETL watchdog.

The script waits until the old watchdog process group has no active Bismark,
Bowtie2, Samtools, or methylation-extractor subprocesses, then terminates the
old group and starts a fresh watchdog with the requested resource profile.
"""
from __future__ import annotations

import argparse
import json
import os
import signal
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path("/home/zdq-as/mouse_methyl_work")
WATCHDOG = ROOT / "scripts" / "etl" / "24_raw_etl_watchdog.py"
DEFAULT_QUEUE = ROOT / "results" / "ralph_v21_raid_raw_clock" / "raw_etl_queue_v25_authorized_full_passed_gate.csv"
DEFAULT_OUT_ROOT = Path("/data/mouse_methyl/processed_v21_raw_etl")
ACTIVE_PATTERNS = (
    "bismark --genome",
    "bowtie2 ",
    "bowtie2-align",
    "samtools ",
    "bismark_methylation_extractor",
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str), encoding="utf-8")
    tmp.replace(path)


def append_jsonl(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(payload, sort_keys=True, default=str) + "\n")


def process_group(pid: int) -> int | None:
    try:
        return os.getpgid(pid)
    except ProcessLookupError:
        return None


def ps_rows() -> list[tuple[int, int, str, str]]:
    completed = subprocess.run(
        ["ps", "-eo", "pid=,pgid=,stat=,cmd="],
        text=True,
        capture_output=True,
        check=False,
    )
    rows: list[tuple[int, int, str, str]] = []
    for line in completed.stdout.splitlines():
        parts = line.strip().split(maxsplit=3)
        if len(parts) < 4:
            continue
        try:
            rows.append((int(parts[0]), int(parts[1]), parts[2], parts[3]))
        except ValueError:
            continue
    return rows


def active_children(pgid: int) -> list[dict[str, Any]]:
    active = []
    for pid, row_pgid, stat, cmd in ps_rows():
        if row_pgid != pgid:
            continue
        if any(pattern in cmd for pattern in ACTIVE_PATTERNS):
            active.append({"pid": pid, "stat": stat, "cmd": cmd[:300]})
    return active


def terminate_group(pgid: int, grace_seconds: int) -> None:
    try:
        os.killpg(pgid, signal.SIGTERM)
    except ProcessLookupError:
        return
    deadline = time.time() + max(1, grace_seconds)
    while time.time() < deadline:
        if not any(row_pgid == pgid for _, row_pgid, _, _ in ps_rows()):
            return
        time.sleep(1)
    try:
        os.killpg(pgid, signal.SIGKILL)
    except ProcessLookupError:
        return


def start_watchdog(args: argparse.Namespace, log_path: Path) -> int:
    cmd = [
        sys.executable,
        "-u",
        str(WATCHDOG),
        "--queue",
        str(args.queue),
        "--out-root",
        str(args.out_root),
        "--poll-seconds",
        str(args.watchdog_poll_seconds),
        "--parallel-samples",
        str(args.parallel_samples),
        "--bismark-threads",
        str(args.bismark_threads),
        "--extract-threads",
        str(args.extract_threads),
        "--start-runner",
        "--authorize-bismark",
    ]
    stdout_log = args.out_root / "v25_rebalanced_watchdog.stdout.log"
    stderr_log = args.out_root / "v25_rebalanced_watchdog.stderr.log"
    with stdout_log.open("a", encoding="utf-8") as out_handle, stderr_log.open("a", encoding="utf-8") as err_handle:
        proc = subprocess.Popen(
            cmd,
            cwd=str(ROOT),
            stdout=out_handle,
            stderr=err_handle,
            text=True,
            start_new_session=True,
        )
    append_jsonl(log_path, {"timestamp": utc_now(), "status": "new_watchdog_started", "pid": proc.pid, "cmd": " ".join(cmd)})
    return int(proc.pid)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--old-watchdog-pid", type=int, required=True)
    parser.add_argument("--queue", type=Path, default=DEFAULT_QUEUE)
    parser.add_argument("--out-root", type=Path, default=DEFAULT_OUT_ROOT)
    parser.add_argument("--parallel-samples", type=int, default=3)
    parser.add_argument("--bismark-threads", type=int, default=6)
    parser.add_argument("--extract-threads", type=int, default=4)
    parser.add_argument("--poll-seconds", type=int, default=5)
    parser.add_argument("--watchdog-poll-seconds", type=int, default=60)
    parser.add_argument("--terminate-grace-seconds", type=int, default=30)
    args = parser.parse_args()

    state_path = args.out_root / "v25_rebalance_state.json"
    log_path = args.out_root / "v25_rebalance_events.jsonl"
    started_at = utc_now()
    append_jsonl(
        log_path,
        {
            "timestamp": started_at,
            "status": "rebalance_wait_started",
            "old_watchdog_pid": args.old_watchdog_pid,
            "target_parallel_samples": args.parallel_samples,
            "target_bismark_threads": args.bismark_threads,
        },
    )

    while True:
        pgid = process_group(args.old_watchdog_pid)
        if pgid is None:
            new_pid = start_watchdog(args, log_path)
            payload = {
                "timestamp": utc_now(),
                "status": "started_new_watchdog_after_old_exit",
                "old_watchdog_pid": args.old_watchdog_pid,
                "new_watchdog_pid": new_pid,
                "started_at": started_at,
            }
            write_json(state_path, payload)
            print(json.dumps(payload, indent=2, sort_keys=True))
            return
        active = active_children(pgid)
        payload = {
            "timestamp": utc_now(),
            "status": "waiting_for_idle" if active else "old_group_idle",
            "old_watchdog_pid": args.old_watchdog_pid,
            "old_pgid": pgid,
            "active_process_count": len(active),
            "active_processes": active[:12],
            "started_at": started_at,
        }
        write_json(state_path, payload)
        append_jsonl(log_path, payload)
        if not active:
            terminate_group(pgid, args.terminate_grace_seconds)
            new_pid = start_watchdog(args, log_path)
            final = {
                "timestamp": utc_now(),
                "status": "rebalanced",
                "old_watchdog_pid": args.old_watchdog_pid,
                "old_pgid": pgid,
                "new_watchdog_pid": new_pid,
                "target_parallel_samples": args.parallel_samples,
                "target_bismark_threads": args.bismark_threads,
                "target_extract_threads": args.extract_threads,
                "started_at": started_at,
            }
            write_json(state_path, final)
            print(json.dumps(final, indent=2, sort_keys=True))
            return
        time.sleep(max(1, args.poll_seconds))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Background watchdog downloader for v11.3 GEO processed supplements.

The downloader is manifest-driven and safe to restart. It downloads files with
HTTP Range chunks, resumes incomplete chunks, writes structured logs, and emits
a watchdog state file with current speed and ETA.
"""
from __future__ import annotations

import argparse
import json
import math
import shutil
import signal
import subprocess
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def format_seconds(seconds: float | int | None) -> str:
    if seconds is None or seconds < 0 or math.isinf(seconds):
        return "unknown"
    seconds = int(seconds)
    hours, rem = divmod(seconds, 3600)
    minutes, sec = divmod(rem, 60)
    if hours:
        return f"{hours}h{minutes:02d}m{sec:02d}s"
    if minutes:
        return f"{minutes}m{sec:02d}s"
    return f"{sec}s"


def append_jsonl(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(payload, sort_keys=True) + "\n")


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    tmp.replace(path)


def file_progress(path: Path, expected: int) -> int:
    if expected <= 0:
        return 0
    if path.exists() and path.stat().st_size == expected:
        return expected
    parts_dir = path.with_name(path.name + ".parts")
    if parts_dir.exists():
        total = 0
        for part in parts_dir.iterdir():
            if part.is_file() and (part.name.startswith("part_") or part.name.endswith(".resume")):
                total += part.stat().st_size
        return min(total, expected)
    if path.exists():
        return min(path.stat().st_size, expected)
    return 0


class DownloadState:
    def __init__(self, rows: list[dict[str, Any]], out_dir: Path):
        self.rows = rows
        self.out_dir = out_dir
        self.lock = threading.Lock()
        self.status: dict[int, dict[str, Any]] = {}
        for row in rows:
            row_id = int(row["row_id"])
            self.status[row_id] = {
                "row_id": row_id,
                "dataset": row["dataset"],
                "supplement_name": row["supplement_name"],
                "status": row.get("status", "planned"),
                "retry_count": int(row.get("retry_count") or 0),
                "last_error": str(row.get("last_error") or ""),
                "started_at": "",
                "finished_at": "",
            }

    def update(self, row_id: int, **kwargs: Any) -> None:
        with self.lock:
            self.status.setdefault(row_id, {"row_id": row_id}).update(kwargs)

    def snapshot_status(self) -> list[dict[str, Any]]:
        with self.lock:
            return [dict(v) for _, v in sorted(self.status.items())]

    def row_status(self, row_id: int) -> str:
        with self.lock:
            return str(self.status.get(row_id, {}).get("status", "unknown"))


def download_chunk(
    url: str,
    parts_dir: Path,
    idx: int,
    start: int,
    end: int,
    *,
    chunk_attempts: int,
    speed_limit: int,
    speed_time: int,
) -> tuple[int, bool, str]:
    chunk_path = parts_dir / f"part_{idx:03d}"
    expected_len = end - start + 1
    if chunk_path.exists() and chunk_path.stat().st_size == expected_len:
        return idx, True, "already_verified"

    details: list[str] = []
    no_progress = 0
    for _ in range(chunk_attempts):
        current = chunk_path.stat().st_size if chunk_path.exists() else 0
        if current == expected_len:
            return idx, True, "verified"
        if current > expected_len:
            chunk_path.unlink()
            current = 0
        resume_start = start + current
        tmp_path = parts_dir / f"part_{idx:03d}.resume"
        tmp_path.unlink(missing_ok=True)
        cmd = [
            "curl",
            "-L",
            "--fail",
            "--retry",
            "4",
            "--retry-delay",
            "5",
            "--connect-timeout",
            "30",
            "--speed-limit",
            str(speed_limit),
            "--speed-time",
            str(speed_time),
            "--range",
            f"{resume_start}-{end}",
            "-o",
            str(tmp_path),
            url,
        ]
        proc = subprocess.run(cmd, text=True, capture_output=True)
        added = tmp_path.stat().st_size if tmp_path.exists() else 0
        if added:
            with chunk_path.open("ab") as out_handle, tmp_path.open("rb") as in_handle:
                shutil.copyfileobj(in_handle, out_handle)
            no_progress = 0
        else:
            no_progress += 1
        tmp_path.unlink(missing_ok=True)
        details.append(proc.stderr[-500:])
        if chunk_path.exists() and chunk_path.stat().st_size == expected_len:
            return idx, True, "verified"
        if no_progress >= 3:
            break
    actual = chunk_path.stat().st_size if chunk_path.exists() else 0
    return idx, False, f"chunk_incomplete:{actual}:{expected_len}:{''.join(details[-2:])}"


def download_file_with_ranges(
    row: dict[str, Any],
    *,
    max_chunks: int,
    chunk_attempts: int,
    speed_limit: int,
    speed_time: int,
) -> tuple[bool, str]:
    url = str(row["official_url"])
    path = Path(str(row["local_path"]))
    expected = int(row["expected_size_bytes"])
    if not url or expected <= 0:
        return False, "missing_url_or_expected_size"
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.stat().st_size == expected:
        return True, "already_verified"

    parts_dir = path.with_name(path.name + ".parts")
    parts_dir.mkdir(parents=True, exist_ok=True)
    chunks = max(1, min(max_chunks, expected))
    chunk_size = math.ceil(expected / chunks)
    ranges: list[tuple[int, int, int]] = []
    for idx in range(chunks):
        start = idx * chunk_size
        end = min(expected - 1, (idx + 1) * chunk_size - 1)
        if start <= end:
            ranges.append((idx, start, end))

    failures: list[dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=min(chunks, len(ranges))) as executor:
        futures = {
            executor.submit(
                download_chunk,
                url,
                parts_dir,
                idx,
                start,
                end,
                chunk_attempts=chunk_attempts,
                speed_limit=speed_limit,
                speed_time=speed_time,
            ): idx
            for idx, start, end in ranges
        }
        for future in as_completed(futures):
            idx, ok, detail = future.result()
            if not ok:
                failures.append({"chunk": idx, "detail": detail})
    if failures:
        return False, json.dumps(failures[:3], sort_keys=True)

    tmp_path = path.with_name(path.name + ".tmp")
    with tmp_path.open("wb") as out_handle:
        for idx, _, _ in ranges:
            chunk_path = parts_dir / f"part_{idx:03d}"
            with chunk_path.open("rb") as in_handle:
                shutil.copyfileobj(in_handle, out_handle)
    actual = tmp_path.stat().st_size
    if actual != expected:
        tmp_path.unlink(missing_ok=True)
        return False, f"assembled_size_mismatch:{actual}:{expected}"
    tmp_path.replace(path)
    shutil.rmtree(parts_dir, ignore_errors=True)
    return True, "downloaded_range_verified"


def monitor_loop(
    state: DownloadState,
    stop_event: threading.Event,
    *,
    state_path: Path,
    status_csv_path: Path,
    poll_seconds: int,
    started_at: float,
) -> None:
    last_t = time.time()
    last_bytes = 0
    while not stop_event.is_set():
        now = time.time()
        rows = state.rows
        total_expected = sum(int(r["expected_size_bytes"]) for r in rows)
        downloaded = sum(file_progress(Path(str(r["local_path"])), int(r["expected_size_bytes"])) for r in rows)
        verified = sum(
            int(r["expected_size_bytes"])
            for r in rows
            if Path(str(r["local_path"])).exists()
            and Path(str(r["local_path"])).stat().st_size == int(r["expected_size_bytes"])
        )
        dt = max(0.001, now - last_t)
        speed = max(0.0, (downloaded - last_bytes) / dt)
        remaining = max(0, total_expected - downloaded)
        eta = remaining / speed if speed > 0 else None
        statuses = state.snapshot_status()
        counts = pd.Series([s.get("status", "unknown") for s in statuses]).value_counts().to_dict() if statuses else {}
        active = [s for s in statuses if s.get("status") in {"downloading", "retrying"}]
        payload = {
            "timestamp": utc_now(),
            "started_at_epoch": started_at,
            "runtime_seconds": round(now - started_at, 1),
            "total_expected_bytes": total_expected,
            "downloaded_bytes_est": downloaded,
            "verified_bytes": verified,
            "remaining_bytes_est": remaining,
            "speed_bps_window": speed,
            "speed_mib_s_window": speed / 1024 / 1024,
            "eta_seconds": eta,
            "eta_human": format_seconds(eta) if eta is not None else "unknown",
            "status_counts": counts,
            "active": active[:10],
        }
        atomic_write_text(state_path, json.dumps(payload, indent=2, sort_keys=True))
        status_rows = []
        status_by_id = {int(s["row_id"]): s for s in statuses}
        for row in rows:
            row_id = int(row["row_id"])
            merged = dict(row)
            merged.update(status_by_id.get(row_id, {}))
            merged["bytes_downloaded_est"] = file_progress(Path(str(row["local_path"])), int(row["expected_size_bytes"]))
            status_rows.append(merged)
        pd.DataFrame(status_rows).to_csv(status_csv_path, index=False)
        last_t = now
        last_bytes = downloaded
        stop_event.wait(poll_seconds)


def process_dataset(
    dataset: str,
    rows: list[dict[str, Any]],
    state: DownloadState,
    *,
    log_path: Path,
    max_retries: int,
    max_file_chunks: int,
    chunk_attempts: int,
    speed_limit: int,
    speed_time: int,
    stop_event: threading.Event,
) -> None:
    for row in rows:
        if stop_event.is_set():
            return
        row_id = int(row["row_id"])
        path = Path(str(row["local_path"]))
        expected = int(row["expected_size_bytes"])
        if path.exists() and path.stat().st_size == expected:
            state.update(row_id, status="verified", finished_at=utc_now(), last_error="")
            append_jsonl(log_path, {"timestamp": utc_now(), "row_id": row_id, "dataset": dataset, "status": "verified_existing", "path": str(path)})
            continue
        state.update(row_id, status="downloading", started_at=utc_now(), last_error="")
        append_jsonl(
            log_path,
            {
                "timestamp": utc_now(),
                "row_id": row_id,
                "dataset": dataset,
                "sample_id": row.get("sample_id"),
                "supplement_name": row.get("supplement_name"),
                "status": "download_started",
                "url": row.get("official_url"),
                "expected_size_bytes": expected,
            },
        )
        last_error = ""
        for retry in range(max_retries + 1):
            if stop_event.is_set():
                return
            if retry:
                state.update(row_id, status="retrying", retry_count=retry, last_error=last_error)
                time.sleep(min(60, 5 * retry))
            ok, detail = download_file_with_ranges(
                row,
                max_chunks=max_file_chunks,
                chunk_attempts=chunk_attempts,
                speed_limit=speed_limit,
                speed_time=speed_time,
            )
            if ok:
                state.update(row_id, status="verified", finished_at=utc_now(), retry_count=retry, last_error="")
                append_jsonl(
                    log_path,
                    {
                        "timestamp": utc_now(),
                        "row_id": row_id,
                        "dataset": dataset,
                        "sample_id": row.get("sample_id"),
                        "supplement_name": row.get("supplement_name"),
                        "status": "verified",
                        "retry_count": retry,
                        "path": str(path),
                        "actual_size": path.stat().st_size if path.exists() else 0,
                    },
                )
                break
            last_error = detail
            append_jsonl(
                log_path,
                {
                    "timestamp": utc_now(),
                    "row_id": row_id,
                    "dataset": dataset,
                    "sample_id": row.get("sample_id"),
                    "supplement_name": row.get("supplement_name"),
                    "status": "download_attempt_failed",
                    "retry_count": retry,
                    "error": detail[-2000:],
                },
            )
        else:
            state.update(row_id, status="failed", finished_at=utc_now(), last_error=last_error)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--out_dir", type=Path, default=ROOT / "results" / "download_backlog_v11_3")
    parser.add_argument("--max_dataset_workers", type=int, default=3)
    parser.add_argument("--max_file_chunks", type=int, default=8)
    parser.add_argument("--max_retries", type=int, default=8)
    parser.add_argument("--chunk_attempts", type=int, default=20)
    parser.add_argument("--poll_seconds", type=int, default=30)
    parser.add_argument("--speed_limit", type=int, default=1024)
    parser.add_argument("--speed_time", type=int, default=180)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    state_path = args.out_dir / "watchdog_state.json"
    status_csv_path = args.out_dir / "download_status.csv"
    log_path = args.out_dir / "download_log.jsonl"

    df = pd.read_csv(args.manifest)
    rows = [
        row
        for row in df.to_dict(orient="records")
        if bool(row.get("download_allowed", True))
    ]
    state = DownloadState(rows, args.out_dir)
    stop_event = threading.Event()

    def handle_signal(signum: int, _frame: Any) -> None:
        append_jsonl(log_path, {"timestamp": utc_now(), "status": "signal_received", "signum": signum})
        stop_event.set()

    signal.signal(signal.SIGTERM, handle_signal)
    signal.signal(signal.SIGINT, handle_signal)

    started_at = time.time()
    append_jsonl(
        log_path,
        {
            "timestamp": utc_now(),
            "status": "watchdog_started",
            "manifest": str(args.manifest),
            "n_rows": len(rows),
            "max_dataset_workers": args.max_dataset_workers,
            "max_file_chunks": args.max_file_chunks,
        },
    )
    monitor = threading.Thread(
        target=monitor_loop,
        args=(state, stop_event),
        kwargs={
            "state_path": state_path,
            "status_csv_path": status_csv_path,
            "poll_seconds": args.poll_seconds,
            "started_at": started_at,
        },
        daemon=True,
    )
    monitor.start()

    by_dataset: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        by_dataset.setdefault(str(row["dataset"]), []).append(row)

    with ThreadPoolExecutor(max_workers=args.max_dataset_workers) as executor:
        futures = [
            executor.submit(
                process_dataset,
                dataset,
                dataset_rows,
                state,
                log_path=log_path,
                max_retries=args.max_retries,
                max_file_chunks=args.max_file_chunks,
                chunk_attempts=args.chunk_attempts,
                speed_limit=args.speed_limit,
                speed_time=args.speed_time,
                stop_event=stop_event,
            )
            for dataset, dataset_rows in sorted(by_dataset.items())
        ]
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as exc:  # noqa: BLE001 - log and continue other datasets
                append_jsonl(log_path, {"timestamp": utc_now(), "status": "dataset_worker_exception", "error": repr(exc)})

    stop_event.set()
    monitor.join(timeout=5)
    append_jsonl(log_path, {"timestamp": utc_now(), "status": "watchdog_finished"})


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Resume-safe ENA-first downloader for v27 raw missing-run repair.

The runner writes status/log/state files and only downloads when
`--authorize-download` is supplied. It does not rename or overwrite existing
FASTQ files. If a target file already exists with an unexpected size or MD5, the
run is blocked for manual review.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")
DEFAULT_QUEUE = ROOT / "results" / "v27_data_repair_execution" / "raw_repair" / "v27_raw_missing_repair_queue.csv"
DEFAULT_OUT = ROOT / "results" / "v27_data_repair_execution" / "raw_repair"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def append_jsonl(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(payload, sort_keys=True, default=str) + "\n")


def atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str), encoding="utf-8")
    tmp.replace(path)


def split_semicolon(value: Any) -> list[str]:
    text = str(value or "").strip()
    if not text or text.lower() == "nan":
        return []
    return [item.strip() for item in text.split(";") if item.strip()]


def expected_sizes(value: Any) -> list[int]:
    sizes = []
    for item in split_semicolon(value):
        try:
            sizes.append(int(float(item)))
        except ValueError:
            sizes.append(0)
    return sizes


def md5sum(path: Path) -> str:
    digest = hashlib.md5()  # noqa: S324 - checksum verification, not security
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_file(path: Path, expected_size: int, expected_md5: str) -> tuple[bool, str]:
    if not path.exists():
        return False, "missing"
    actual = path.stat().st_size
    if expected_size > 0 and actual != expected_size:
        return False, f"size_mismatch:{actual}:{expected_size}"
    if expected_md5 and md5sum(path).lower() != expected_md5.lower():
        return False, "md5_mismatch"
    return True, "verified"


def file_name_from_url(url: str) -> str:
    name = Path(urlparse(url).path).name
    if not name:
        raise ValueError(f"cannot infer filename from URL: {url}")
    return name


def interleave_by_dataset(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Keep per-dataset priority order while avoiding one-source bursts."""
    groups: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        groups.setdefault(str(row.get("dataset", "")), []).append(row)
    ordered_keys = list(groups)
    interleaved = []
    while any(groups.values()):
        for key in ordered_keys:
            if groups[key]:
                interleaved.append(groups[key].pop(0))
    return interleaved


def observed_target_bytes(rows: list[dict[str, Any]]) -> tuple[int, dict[str, int]]:
    total = 0
    by_dataset: dict[str, int] = {}
    for row in rows:
        dataset = str(row.get("dataset", ""))
        target_dir = Path(str(row.get("target_fastq_dir", "")))
        for url in split_semicolon(row.get("fastq_urls", "")):
            try:
                name = file_name_from_url(url)
            except ValueError:
                continue
            final_path = target_dir / name
            part_path = target_dir / f"{name}.v27part"
            path = final_path if final_path.exists() else part_path
            if not path.exists():
                continue
            try:
                size = path.stat().st_size
            except FileNotFoundError:
                continue
            total += size
            by_dataset[dataset] = by_dataset.get(dataset, 0) + size
    return total, by_dataset


def run_curl(url: str, part_path: Path, args: argparse.Namespace) -> tuple[bool, str]:
    part_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "curl",
        "-L",
        "--fail",
        "--retry",
        str(args.curl_retries),
        "--retry-connrefused",
        "--retry-all-errors",
        "--retry-max-time",
        str(args.retry_max_time),
        "--retry-delay",
        "5",
        "--connect-timeout",
        "30",
        "--speed-limit",
        str(args.speed_limit),
        "--speed-time",
        str(args.speed_time),
        "--continue-at",
        "-",
        "-o",
        str(part_path),
        url,
    ]
    proc = subprocess.run(cmd, text=True, capture_output=True)
    if proc.returncode != 0:
        return False, proc.stderr[-2000:] or f"curl_exit_{proc.returncode}"
    return True, "curl_completed"


def download_one_file(
    *,
    url: str,
    target_dir: Path,
    expected_size: int,
    expected_md5: str,
    args: argparse.Namespace,
) -> tuple[bool, str, str, int]:
    name = file_name_from_url(url)
    final_path = target_dir / name
    ok, detail = verify_file(final_path, expected_size, expected_md5)
    if ok:
        return True, str(final_path), "already_verified", final_path.stat().st_size
    if final_path.exists() and detail != "missing":
        return False, str(final_path), f"blocked_existing_{detail}", final_path.stat().st_size

    part_path = target_dir / f"{name}.v27part"
    if expected_size > 0 and part_path.exists() and part_path.stat().st_size > expected_size:
        part_path.unlink()
    ok, detail = run_curl(url, part_path, args)
    if not ok:
        size = part_path.stat().st_size if part_path.exists() else 0
        return False, str(final_path), detail, size
    ok, detail = verify_file(part_path, expected_size, expected_md5)
    if not ok:
        size = part_path.stat().st_size if part_path.exists() else 0
        return False, str(final_path), f"post_download_{detail}", size
    if final_path.exists():
        return False, str(final_path), "blocked_target_appeared_during_download", final_path.stat().st_size
    part_path.replace(final_path)
    return True, str(final_path), "downloaded_verified", final_path.stat().st_size


class Status:
    def __init__(self, rows: list[dict[str, Any]], out_dir: Path):
        self.rows = rows
        self.out_dir = out_dir
        self.lock = threading.Lock()
        self.status: dict[str, dict[str, Any]] = {
            str(row["run_accession"]): {
                "dataset": row["dataset"],
                "run_accession": row["run_accession"],
                "gsm": row.get("gsm", ""),
                "status": "planned",
                "started_at": "",
                "finished_at": "",
                "downloaded_bytes": 0,
                "last_error": "",
            }
            for row in rows
        }

    def update(self, run: str, **kwargs: Any) -> None:
        with self.lock:
            self.status.setdefault(run, {"run_accession": run}).update(kwargs)

    def snapshot(self) -> list[dict[str, Any]]:
        with self.lock:
            return [dict(v) for _, v in sorted(self.status.items())]


def monitor_loop(status: Status, stop: threading.Event, out_dir: Path, poll_seconds: int, started: float) -> None:
    state_path = out_dir / "download_watchdog_state.json"
    status_path = out_dir / "download_status.csv"
    while not stop.is_set():
        rows = status.snapshot()
        counts = pd.Series([r.get("status", "unknown") for r in rows]).value_counts().to_dict() if rows else {}
        active = [r for r in rows if r.get("status") == "downloading"]
        total_bytes = sum(int(r.get("downloaded_bytes", 0) or 0) for r in rows)
        observed_file_bytes, observed_by_dataset = observed_target_bytes(status.rows)
        atomic_write_json(
            state_path,
            {
                "timestamp": utc_now(),
                "runtime_seconds": round(time.time() - started, 1),
                "status_counts": counts,
                "active": active[:12],
                "downloaded_bytes_observed": total_bytes,
                "downloaded_gib_observed": round(total_bytes / 1024**3, 6),
                "observed_file_bytes": observed_file_bytes,
                "observed_file_gib": round(observed_file_bytes / 1024**3, 6),
                "observed_file_gib_by_dataset": {
                    key: round(value / 1024**3, 6) for key, value in sorted(observed_by_dataset.items())
                },
            },
        )
        pd.DataFrame(rows).to_csv(status_path, index=False)
        stop.wait(max(5, poll_seconds))


def process_run(row: dict[str, Any], status: Status, args: argparse.Namespace, log_path: Path) -> None:
    run = str(row["run_accession"])
    if not bool(row.get("download_allowed", True)):
        reason = str(row.get("blocked_reason") or "download_not_allowed")
        status.update(run, status="blocked", finished_at=utc_now(), last_error=reason)
        append_jsonl(log_path, {"timestamp": utc_now(), "run_accession": run, "status": "blocked", "reason": reason})
        return
    urls = split_semicolon(row.get("fastq_urls", ""))
    sizes = expected_sizes(row.get("fastq_bytes", ""))
    md5s = split_semicolon(row.get("fastq_md5s", ""))
    if not urls:
        status.update(run, status="blocked", finished_at=utc_now(), last_error="missing_fastq_urls")
        return
    target_dir = Path(str(row["target_fastq_dir"]))
    status.update(run, status="downloading", started_at=utc_now(), n_files=len(urls))
    append_jsonl(log_path, {"timestamp": utc_now(), "run_accession": run, "dataset": row.get("dataset"), "status": "run_started", "n_files": len(urls)})
    downloaded_bytes = 0
    outputs = []
    for idx, url in enumerate(urls):
        expected_size = sizes[idx] if idx < len(sizes) else 0
        expected_md5 = md5s[idx] if idx < len(md5s) else ""
        ok, path, detail, observed = download_one_file(
            url=url,
            target_dir=target_dir,
            expected_size=expected_size,
            expected_md5=expected_md5,
            args=args,
        )
        downloaded_bytes += observed
        outputs.append(path)
        append_jsonl(
            log_path,
            {
                "timestamp": utc_now(),
                "run_accession": run,
                "url": url,
                "path": path,
                "status": "file_verified" if ok else "file_failed",
                "detail": detail,
                "observed_bytes": observed,
                "expected_bytes": expected_size,
            },
        )
        status.update(run, downloaded_bytes=downloaded_bytes, output_paths=";".join(outputs), last_error="" if ok else detail)
        if not ok:
            status.update(run, status="failed", finished_at=utc_now())
            return
    status.update(run, status="verified", finished_at=utc_now(), downloaded_bytes=downloaded_bytes, output_paths=";".join(outputs), last_error="")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--queue", type=Path, default=DEFAULT_QUEUE)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--max-workers", type=int, default=4)
    parser.add_argument("--max-runs", type=int, default=0)
    parser.add_argument("--poll-seconds", type=int, default=30)
    parser.add_argument("--curl-retries", type=int, default=8)
    parser.add_argument("--retry-max-time", type=int, default=1800)
    parser.add_argument("--speed-limit", type=int, default=1024)
    parser.add_argument("--speed-time", type=int, default=300)
    parser.add_argument("--interleave-datasets", action="store_true")
    parser.add_argument("--authorize-download", action="store_true")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    log_path = args.out_dir / "download_log.jsonl"
    queue = pd.read_csv(args.queue)
    rows = queue.sort_values("repair_priority").to_dict(orient="records")
    if args.interleave_datasets:
        rows = interleave_by_dataset(rows)
    if args.max_runs:
        rows = rows[: args.max_runs]
    if not args.authorize_download:
        payload = {
            "timestamp": utc_now(),
            "status": "blocked_authorization_required",
            "queue": str(args.queue),
            "n_rows_selected": len(rows),
            "message": "rerun with --authorize-download to start ENA FASTQ repair",
        }
        atomic_write_json(args.out_dir / "download_watchdog_state.json", payload)
        print(json.dumps(payload, indent=2, sort_keys=True))
        raise SystemExit(2)

    status = Status(rows, args.out_dir)
    stop = threading.Event()
    started = time.time()
    monitor = threading.Thread(target=monitor_loop, args=(status, stop, args.out_dir, args.poll_seconds, started), daemon=True)
    monitor.start()
    append_jsonl(log_path, {"timestamp": utc_now(), "status": "watchdog_started", "queue": str(args.queue), "n_rows": len(rows), "max_workers": args.max_workers})
    with ThreadPoolExecutor(max_workers=max(1, args.max_workers)) as executor:
        futures = [executor.submit(process_run, row, status, args, log_path) for row in rows]
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as exc:  # noqa: BLE001 - keep other runs moving
                append_jsonl(log_path, {"timestamp": utc_now(), "status": "run_worker_exception", "error": repr(exc)})
    stop.set()
    monitor.join(timeout=5)
    rows_out = status.snapshot()
    pd.DataFrame(rows_out).to_csv(args.out_dir / "download_status.csv", index=False)
    counts = pd.Series([r.get("status", "unknown") for r in rows_out]).value_counts().to_dict() if rows_out else {}
    final = {
        "timestamp": utc_now(),
        "status": "completed",
        "queue": str(args.queue),
        "n_rows_selected": len(rows),
        "status_counts": counts,
        "runtime_seconds": round(time.time() - started, 1),
    }
    atomic_write_json(args.out_dir / "download_watchdog_state.json", final)
    append_jsonl(log_path, {"timestamp": utc_now(), "status": "watchdog_finished", "status_counts": counts})
    print(json.dumps(final, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

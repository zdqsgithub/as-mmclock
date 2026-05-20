#!/usr/bin/env python3
"""
Download GEO supplementary files with resume and structured JSONL logging.

Files are downloaded to a .part path first, then atomically renamed after the
downloaded size matches the inventory remote size. Existing verified files are
not overwritten.
"""
from __future__ import annotations

import argparse
import csv
import json
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path("/home/zdq-as/mouse_methyl_work")
META_DIR = ROOT / "metadata"
DEFAULT_INVENTORY = META_DIR / "geo_supplement_inventory.csv"
DEFAULT_LOG = META_DIR / "download_log.jsonl"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def append_log(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps({"timestamp": utc_now(), **payload}) + "\n")


def read_inventory(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def download_with_resume(row: dict, log_path: Path, timeout: int = 300, retries: int = 5) -> dict:
    dataset = row["dataset"]
    url = row["supplement_url"]
    target = Path(row["local_path"])
    remote_size = int(row["remote_size_bytes"] or 0)
    target.parent.mkdir(parents=True, exist_ok=True)
    part = target.with_suffix(target.suffix + ".part")

    if target.exists() and remote_size and target.stat().st_size == remote_size:
        payload = {
            "dataset": dataset,
            "file_accession": target.name,
            "status": "skipped",
            "reason": "already_verified",
            "bytes_downloaded": target.stat().st_size,
            "verified": True,
            "failure_count": 0,
            "error": "",
        }
        append_log(log_path, payload)
        return payload

    for attempt in range(1, retries + 1):
        existing = part.stat().st_size if part.exists() else 0
        headers = {}
        if existing > 0:
            headers["Range"] = f"bytes={existing}-"
        request = urllib.request.Request(url, headers=headers)
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response, part.open("ab") as handle:
                while True:
                    chunk = response.read(1024 * 1024)
                    if not chunk:
                        break
                    handle.write(chunk)

            downloaded = part.stat().st_size
            if remote_size and downloaded != remote_size:
                raise RuntimeError(f"downloaded size {downloaded} != remote size {remote_size}")
            part.replace(target)
            payload = {
                "dataset": dataset,
                "file_accession": target.name,
                "status": "completed",
                "bytes_downloaded": target.stat().st_size,
                "verified": bool(remote_size and target.stat().st_size == remote_size),
                "failure_count": attempt - 1,
                "error": "",
            }
            append_log(log_path, payload)
            return payload
        except urllib.error.HTTPError as exc:
            # Some GEO endpoints ignore Range and return 200. Restart once.
            if exc.code == 416 and remote_size and existing == remote_size:
                part.replace(target)
                payload = {
                    "dataset": dataset,
                    "file_accession": target.name,
                    "status": "completed",
                    "bytes_downloaded": target.stat().st_size,
                    "verified": True,
                    "failure_count": attempt - 1,
                    "error": "",
                }
                append_log(log_path, payload)
                return payload
            error = str(exc)[:500]
        except Exception as exc:
            error = str(exc)[:500]

        payload = {
            "dataset": dataset,
            "file_accession": target.name,
            "status": "failed",
            "bytes_downloaded": part.stat().st_size if part.exists() else 0,
            "verified": False,
            "failure_count": attempt,
            "error": error,
        }
        append_log(log_path, payload)
        if attempt < retries:
            time.sleep(min(2**attempt, 300))

    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inventory", default=str(DEFAULT_INVENTORY))
    parser.add_argument("--dataset", default="GSE80672")
    parser.add_argument("--log", default=str(DEFAULT_LOG))
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument("--retries", type=int, default=5)
    args = parser.parse_args()

    rows = [row for row in read_inventory(Path(args.inventory)) if row["dataset"] == args.dataset]
    if not rows:
        raise SystemExit(f"No inventory rows for {args.dataset}. Run 08_geo_supplement_inventory.py first.")

    for row in rows:
        result = download_with_resume(row, Path(args.log), timeout=args.timeout, retries=args.retries)
        print(
            f"[Download] {row['dataset']} {row['preferred_processed_file']} "
            f"status={result['status']} verified={result['verified']} bytes={result['bytes_downloaded']}"
        )


if __name__ == "__main__":
    main()

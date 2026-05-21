#!/usr/bin/env python3
"""Build the v27 raw missing-run repair queue for priority datasets.

This script is read-only with respect to FASTQ files. It narrows the v26 local
missing-run audit to the v27 priority datasets, checks whether runs are still
missing locally, fetches ENA direct FASTQ metadata when possible, and writes a
resume-safe queue for the downloader. It does not download.
"""
from __future__ import annotations

import argparse
import csv
import json
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")
DEFAULT_IN = ROOT / "results" / "v26_data_repair_research" / "local_raw_missing_fastq_repair_queue.csv"
DEFAULT_OUT = ROOT / "results" / "v27_data_repair_execution" / "raw_repair"
RAW_ROOT = Path("/data/mouse_methyl/raw")
PRIORITY_DATASETS = ["GSE80672", "GSE121141"]

ERROR_RANK = {
    "orphan_lock_retry": 0,
    "fasterq_timeout_retry_or_ena": 0,
    "not_attempted_or_log_missing": 1,
}
DATASET_RANK = {"GSE80672": 0, "GSE121141": 1}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def fastq_paths_for_run(dataset: str, run: str) -> list[Path]:
    fastq_dir = RAW_ROOT / dataset / "fastq"
    if not fastq_dir.exists():
        return []
    return sorted(
        path
        for pattern in (f"{run}*.fastq", f"{run}*.fastq.gz", f"{run}*.fq", f"{run}*.fq.gz")
        for path in fastq_dir.glob(pattern)
        if path.is_file()
    )


def local_complete(dataset: str, run: str, layout: str) -> tuple[bool, int, str]:
    paths = fastq_paths_for_run(dataset, run)
    if str(layout).upper() == "PAIRED":
        names = [p.name for p in paths]
        has_r1 = any(token in name for name in names for token in ("_1.fastq", "_R1", "_1.fq"))
        has_r2 = any(token in name for name in names for token in ("_2.fastq", "_R2", "_2.fq"))
        complete = len(paths) >= 2 and has_r1 and has_r2
    else:
        complete = len(paths) >= 1
    return complete, len(paths), ";".join(str(p) for p in paths)


def split_ena_values(value: str) -> list[str]:
    value = str(value or "").strip()
    if not value or value.lower() == "nan":
        return []
    return [item.strip() for item in value.split(";") if item.strip()]


def fetch_ena_run(run: str, timeout: int = 45) -> dict[str, Any]:
    fields = "run_accession,fastq_ftp,fastq_md5,fastq_bytes,library_layout"
    query = urllib.parse.urlencode(
        {
            "accession": run,
            "result": "read_run",
            "fields": fields,
            "format": "tsv",
        }
    )
    url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?{query}"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "as-mmclock-v27/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as response:
            text = response.read().decode("utf-8", errors="replace")
    except Exception as exc:  # noqa: BLE001 - reported in manifest
        return {"run_accession": run, "ena_error": str(exc)[:500]}
    lines = [line for line in text.splitlines() if line.strip()]
    if len(lines) < 2:
        return {"run_accession": run, "ena_error": "no_ena_filereport_row"}
    row = next(csv.DictReader(lines, delimiter="\t"), None)
    if not row:
        return {"run_accession": run, "ena_error": "ena_parse_empty"}
    ftp_values = split_ena_values(row.get("fastq_ftp", ""))
    urls = [
        value if value.startswith(("http://", "https://")) else f"https://{value}"
        for value in ftp_values
    ]
    md5s = split_ena_values(row.get("fastq_md5", ""))
    sizes = split_ena_values(row.get("fastq_bytes", ""))
    return {
        "run_accession": run,
        "ena_error": "",
        "ena_library_layout": row.get("library_layout", ""),
        "fastq_urls": ";".join(urls),
        "fastq_md5s": ";".join(md5s),
        "fastq_bytes": ";".join(sizes),
        "ena_n_fastq": len(urls),
    }


def build_queue(args: argparse.Namespace) -> tuple[pd.DataFrame, dict[str, Any]]:
    missing = pd.read_csv(args.input)
    datasets = set(args.datasets)
    priority = missing[missing["dataset"].astype(str).isin(datasets)].copy()
    if priority.empty:
        queue = pd.DataFrame()
    else:
        with ThreadPoolExecutor(max_workers=args.ena_workers) as executor:
            futures = {
                executor.submit(fetch_ena_run, str(row.run_accession)): str(row.run_accession)
                for row in priority.itertuples(index=False)
            }
            ena_rows = []
            for future in as_completed(futures):
                ena_rows.append(future.result())
                time.sleep(args.ena_spacing_seconds)
        ena = pd.DataFrame(ena_rows)
        queue = priority.merge(ena, on="run_accession", how="left")
        local_checks = [
            local_complete(str(row.dataset), str(row.run_accession), str(row.library_layout))
            for row in queue.itertuples(index=False)
        ]
        queue["current_local_fastq_complete"] = [item[0] for item in local_checks]
        queue["current_local_fastq_count"] = [item[1] for item in local_checks]
        queue["current_local_fastq_paths"] = [item[2] for item in local_checks]
        queue["target_fastq_dir"] = [str(RAW_ROOT / str(ds) / "fastq") for ds in queue["dataset"]]
        queue["dataset_rank"] = queue["dataset"].astype(str).map(DATASET_RANK).fillna(99).astype(int)
        queue["repair_error_rank"] = (
            queue["latest_log_error_class"].fillna("not_attempted_or_log_missing").astype(str).map(ERROR_RANK).fillna(2).astype(int)
        )
        queue["download_strategy"] = "ena_direct_first_then_sra_if_needed"
        queue["download_allowed"] = ~queue["current_local_fastq_complete"].astype(bool)
        queue["blocked_reason"] = ""
        queue.loc[queue["current_local_fastq_complete"], "blocked_reason"] = "already_complete_local_fastq"
        queue.loc[queue["fastq_urls"].fillna("").eq(""), "download_allowed"] = False
        queue.loc[queue["fastq_urls"].fillna("").eq(""), "blocked_reason"] = "missing_ena_fastq_urls"
        queue["repair_priority"] = range(1, len(queue) + 1)
        queue = queue.sort_values(
            ["repair_error_rank", "dataset_rank", "run_total_bases", "run_accession"],
            ascending=[True, True, False, True],
        ).reset_index(drop=True)
        queue["repair_priority"] = queue.index + 1

    summary = {
        "timestamp": utc_now(),
        "status": "completed",
        "source_path": str(args.input),
        "datasets_requested": list(args.datasets),
        "n_rows": int(len(queue)),
        "n_download_allowed": int(queue["download_allowed"].sum()) if not queue.empty else 0,
        "n_already_complete": int(queue["current_local_fastq_complete"].sum()) if not queue.empty else 0,
        "n_missing_ena": int(queue["fastq_urls"].fillna("").eq("").sum()) if not queue.empty else 0,
        "download_started": False,
    }
    return queue, summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_IN)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--datasets", nargs="+", default=PRIORITY_DATASETS)
    parser.add_argument("--ena-workers", type=int, default=8)
    parser.add_argument("--ena-spacing-seconds", type=float, default=0.0)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    queue, summary = build_queue(args)
    queue_path = args.out_dir / "v27_raw_missing_repair_queue.csv"
    summary_path = args.out_dir / "v27_raw_missing_repair_queue_summary.json"
    queue.to_csv(queue_path, index=False)
    summary["queue_path"] = str(queue_path)
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

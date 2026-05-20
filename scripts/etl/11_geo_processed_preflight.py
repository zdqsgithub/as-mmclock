#!/usr/bin/env python3
"""Preflight GEO processed methylation supplements before large downloads.

The preflight records official GEO FTP facts only: supplement directory listing,
``filelist.txt`` rows when present, HEAD metadata, schema guesses, and whether a
file should enter the v6 benchmark. It intentionally does not download large
supplement files.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path("/home/zdq-as/mouse_methyl_work")
META_DIR = ROOT / "metadata"

OFFICIAL_DOCS = [
    "https://www.ncbi.nlm.nih.gov/geo/info/download.html",
    "https://www.ncbi.nlm.nih.gov/geo/info/geo_paccess.html",
    "https://www.ncbi.nlm.nih.gov/geo/info/soft.html",
]

RULES = {
    "GSE213628": [
        ("GSE213628_RAW.tar", "bismark_cov_per_sample_tar", 1, False),
    ],
    "GSE80672": [
        ("GSE80672_RAW.tar", "gse80672_overlap_percentage_coverage", 1, True),
    ],
    "GSE93957": [
        ("GSE93957_RAW.tar", "bismark_cov_per_sample_tar", 1, True),
        ("GSE93957_exchange.tar.gz", "bismark_cov_per_sample_tar_alt_archive", 2, False),
    ],
    "GSE121141": [
        ("GSE121141_RAW.tar", "bismark_cov_per_sample_tar", 1, True),
    ],
    "GSE60012": [
        ("GSE60012_100bpTiles_RRBS_Mouse.txt.gz", "gse60012_100bp_tile_matrix", 1, True),
    ],
    "GSE52266": [
        ("GSE52266_SerreRRBSPercentsNov2013.txt.gz", "gse52266_percent_matrix_low_priority", 5, False),
        ("GSE52266_SerreRRBSCountsNov2013.txt.gz", "gse52266_count_matrix_low_priority", 6, False),
    ],
}


def series_bucket(gse: str) -> str:
    match = re.fullmatch(r"GSE(\d+)", gse)
    if not match:
        raise ValueError(f"Invalid GSE accession: {gse}")
    digits = match.group(1)
    return f"GSE{digits[:-3]}nnn"


def suppl_dir_url(gse: str) -> str:
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{series_bucket(gse)}/{gse}/suppl/"


def fetch_text(url: str, timeout: int = 30) -> str:
    try:
        with urllib.request.urlopen(url, timeout=timeout) as response:
            return response.read().decode("utf-8", errors="replace")
    except Exception:
        return ""


def head_url(url: str, timeout: int = 30) -> dict:
    request = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            headers = response.headers
            return {
                "http_status": int(response.status),
                "remote_size_bytes": int(headers.get("Content-Length") or 0) or None,
                "remote_last_modified": headers.get("Last-Modified"),
                "content_type": headers.get("Content-Type"),
                "available": True,
                "head_error": "",
            }
    except urllib.error.HTTPError as exc:
        return {
            "http_status": int(exc.code),
            "remote_size_bytes": None,
            "remote_last_modified": None,
            "content_type": None,
            "available": False,
            "head_error": str(exc)[:300],
        }
    except Exception as exc:
        return {
            "http_status": None,
            "remote_size_bytes": None,
            "remote_last_modified": None,
            "content_type": None,
            "available": False,
            "head_error": str(exc)[:300],
        }


def parse_filelist(text: str) -> dict[str, dict]:
    rows: dict[str, dict] = {}
    reader = csv.DictReader(text.splitlines(), delimiter="\t")
    for row in reader:
        name = row.get("Name")
        if name:
            rows[name] = row
    return rows


def directory_names(html: str) -> set[str]:
    names = set(re.findall(r'href="([^"]+)"', html))
    names.update(re.findall(r"\s([A-Za-z0-9_.+-]+\.(?:tar|gz|txt|csv|tsv|bed)(?:\.gz)?)\s", html))
    return {name for name in names if not name.startswith("?")}


def build_rows(gses: list[str]) -> list[dict]:
    rows = []
    timestamp = datetime.now(timezone.utc).isoformat()
    for gse in gses:
        base = suppl_dir_url(gse)
        listing_text = fetch_text(base)
        filelist_text = fetch_text(base + "filelist.txt")
        filelist = parse_filelist(filelist_text)
        listing_names = directory_names(listing_text)
        for filename, schema, priority, include in RULES.get(gse, [(f"{gse}_RAW.tar", "unknown_raw_tar", 9, False)]):
            url = base + filename
            head = head_url(url)
            filelist_row = filelist.get(filename, {})
            listed = filename in listing_names or filename in filelist
            blocker = ""
            if not listed:
                blocker = "not_in_ftp_listing"
            if not head["available"]:
                blocker = "head_not_available"
            if "low_priority" in schema:
                blocker = "low_priority_not_in_v6_main_benchmark"
            rows.append(
                {
                    "dataset": gse,
                    "supplement_url": url,
                    "supplement_dir_url": base,
                    "preferred_processed_file": filename,
                    "schema_guess": schema,
                    "download_priority": priority,
                    "include_in_v6_benchmark": include,
                    "listed_in_ftp": listed,
                    "filelist_type": filelist_row.get("Type", ""),
                    "filelist_size": filelist_row.get("Size", ""),
                    "filelist_time": filelist_row.get("Time", ""),
                    "preflight_blocker": blocker,
                    "preflight_timestamp_utc": timestamp,
                    **head,
                }
            )
    return rows


def write_outputs(rows: list[dict], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "dataset",
        "supplement_url",
        "supplement_dir_url",
        "preferred_processed_file",
        "schema_guess",
        "download_priority",
        "include_in_v6_benchmark",
        "listed_in_ftp",
        "available",
        "http_status",
        "remote_size_bytes",
        "remote_last_modified",
        "content_type",
        "filelist_type",
        "filelist_size",
        "filelist_time",
        "preflight_blocker",
        "preflight_timestamp_utc",
        "head_error",
    ]
    with out_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    payload = {"official_docs": OFFICIAL_DOCS, "rows": rows}
    out_path.with_suffix(".json").write_text(json.dumps(payload, indent=2), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--datasets", default="GSE80672,GSE93957,GSE121141,GSE60012,GSE52266,GSE213628")
    parser.add_argument("--out", default=str(META_DIR / "geo_processed_preflight.csv"))
    args = parser.parse_args()
    gses = [item.strip() for item in args.datasets.split(",") if item.strip()]
    rows = build_rows(gses)
    write_outputs(rows, Path(args.out))
    print(f"[Preflight] Wrote {args.out}")
    for row in rows:
        print(
            f"  {row['dataset']} {row['preferred_processed_file']} "
            f"schema={row['schema_guess']} available={row['available']} "
            f"size={row['remote_size_bytes']} blocker={row['preflight_blocker'] or 'none'}"
        )


if __name__ == "__main__":
    main()

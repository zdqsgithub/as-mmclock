#!/usr/bin/env python3
"""Build an inventory of GEO supplementary methylation files.

The v6 inventory is dataset-aware. Some GEO series expose a ``GSE*_RAW.tar`` of
per-sample coverage files, while older intervention datasets expose aggregated
matrix files. This script records the preferred processed supplement, schema
guess, and local verification status without downloading data.
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
from typing import Iterable

ROOT = Path("/home/zdq-as/mouse_methyl_work")
META_DIR = ROOT / "metadata"
RAW_SUPP_DIR = ROOT / "raw_downloads" / "geo_supplements"


def series_bucket(gse: str) -> str:
    match = re.fullmatch(r"GSE(\d+)", gse)
    if not match:
        raise ValueError(f"Invalid GSE accession: {gse}")
    digits = match.group(1)
    prefix = digits[:-3] + "nnn"
    return f"GSE{prefix}"


SUPPLEMENT_RULES = {
    "GSE213628": [
        {
            "filename": "GSE213628_RAW.tar",
            "file_type": "GEO_RAW_TAR_OF_COV",
            "schema_guess": "bismark_cov_per_sample_tar",
            "download_priority": 1,
            "include_in_v6_benchmark": False,
        }
    ],
    "GSE80672": [
        {
            "filename": "GSE80672_RAW.tar",
            "file_type": "GEO_RAW_TAR_OF_TXT",
            "schema_guess": "gse80672_overlap_percentage_coverage",
            "download_priority": 1,
            "include_in_v6_benchmark": True,
        }
    ],
    "GSE93957": [
        {
            "filename": "GSE93957_RAW.tar",
            "file_type": "GEO_RAW_TAR_OF_COV_TXT",
            "schema_guess": "bismark_cov_per_sample_tar",
            "download_priority": 1,
            "include_in_v6_benchmark": True,
        }
    ],
    "GSE121141": [
        {
            "filename": "GSE121141_RAW.tar",
            "file_type": "GEO_RAW_TAR_OF_BISMARK_COV",
            "schema_guess": "bismark_cov_per_sample_tar",
            "download_priority": 1,
            "include_in_v6_benchmark": True,
        }
    ],
    "GSE60012": [
        {
            "filename": "GSE60012_100bpTiles_RRBS_Mouse.txt.gz",
            "file_type": "GEO_AGGREGATED_TILE_MATRIX",
            "schema_guess": "gse60012_100bp_tile_matrix",
            "download_priority": 1,
            "include_in_v6_benchmark": True,
        }
    ],
    "GSE52266": [
        {
            "filename": "GSE52266_SerreRRBSPercentsNov2013.txt.gz",
            "file_type": "GEO_AGGREGATED_PERCENT_MATRIX",
            "schema_guess": "gse52266_percent_matrix_low_priority",
            "download_priority": 5,
            "include_in_v6_benchmark": False,
        },
        {
            "filename": "GSE52266_SerreRRBSCountsNov2013.txt.gz",
            "file_type": "GEO_AGGREGATED_COUNT_MATRIX",
            "schema_guess": "gse52266_count_matrix_low_priority",
            "download_priority": 6,
            "include_in_v6_benchmark": False,
        },
    ],
}


def supplement_url(gse: str, filename: str) -> str:
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{series_bucket(gse)}/{gse}/suppl/{filename}"


def supplement_rules(gse: str) -> list[dict]:
    if gse in SUPPLEMENT_RULES:
        return SUPPLEMENT_RULES[gse]
    return [
        {
            "filename": f"{gse}_RAW.tar",
            "file_type": "GEO_RAW_TAR_UNKNOWN_SCHEMA",
            "schema_guess": "unknown_raw_tar",
            "download_priority": 9,
            "include_in_v6_benchmark": False,
        }
    ]


def head_url(url: str, timeout: int = 30) -> dict:
    request = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            headers = response.headers
            return {
                "http_status": int(response.status),
                "remote_size_bytes": int(headers.get("Content-Length") or 0) or None,
                "remote_last_modified": headers.get("Last-Modified"),
                "available": True,
                "error": "",
            }
    except urllib.error.HTTPError as exc:
        return {
            "http_status": int(exc.code),
            "remote_size_bytes": None,
            "remote_last_modified": None,
            "available": False,
            "error": str(exc)[:500],
        }
    except Exception as exc:
        return {
            "http_status": None,
            "remote_size_bytes": None,
            "remote_last_modified": None,
            "available": False,
            "error": str(exc)[:500],
        }


def local_status(gse: str, url: str, remote_size: int | None) -> dict:
    filename = Path(urllib.request.urlparse(url).path).name
    local_path = RAW_SUPP_DIR / gse / filename
    if not local_path.exists():
        return {
            "local_path": str(local_path),
            "local_size_bytes": 0,
            "local_status": "missing",
            "verified": False,
        }
    local_size = local_path.stat().st_size
    verified = bool(remote_size and local_size == remote_size)
    if verified:
        status = "verified"
    elif local_size > 0:
        status = "partial_or_size_mismatch"
    else:
        status = "empty"
    return {
        "local_path": str(local_path),
        "local_size_bytes": local_size,
        "local_status": status,
        "verified": verified,
    }


def build_inventory(gses: Iterable[str]) -> list[dict]:
    rows = []
    for gse in gses:
        for rule in supplement_rules(gse):
            url = supplement_url(gse, rule["filename"])
            remote = head_url(url)
            local = local_status(gse, url, remote.get("remote_size_bytes"))
            rows.append(
                {
                    "dataset": gse,
                    "supplement_url": url,
                    "preferred_processed_file": Path(urllib.request.urlparse(url).path).name,
                    "file_type": rule["file_type"],
                    "schema_guess": rule["schema_guess"],
                    "download_priority": rule["download_priority"],
                    "include_in_v6_benchmark": rule["include_in_v6_benchmark"],
                    "inventory_timestamp_utc": datetime.now(timezone.utc).isoformat(),
                    **remote,
                    **local,
                }
            )
    return rows


def write_inventory(rows: list[dict], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "dataset",
        "supplement_url",
        "preferred_processed_file",
        "file_type",
        "schema_guess",
        "download_priority",
        "include_in_v6_benchmark",
        "available",
        "http_status",
        "remote_size_bytes",
        "remote_last_modified",
        "local_path",
        "local_size_bytes",
        "local_status",
        "verified",
        "inventory_timestamp_utc",
        "error",
    ]
    with out_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    (out_path.with_suffix(".json")).write_text(json.dumps(rows, indent=2), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--datasets",
        default="GSE80672,GSE93957,GSE121141,GSE60012,GSE52266,GSE213628",
        help="Comma-separated GEO series accessions.",
    )
    parser.add_argument("--out", default=str(META_DIR / "geo_supplement_inventory.csv"))
    args = parser.parse_args()

    gses = [item.strip() for item in args.datasets.split(",") if item.strip()]
    rows = build_inventory(gses)
    out_path = Path(args.out)
    write_inventory(rows, out_path)
    print(f"[Inventory] Wrote {out_path}")
    for row in rows:
        print(
            f"  {row['dataset']}: available={row['available']} "
            f"file={row['preferred_processed_file']} schema={row['schema_guess']} "
            f"remote={row['remote_size_bytes']} local_status={row['local_status']}"
        )


if __name__ == "__main__":
    main()

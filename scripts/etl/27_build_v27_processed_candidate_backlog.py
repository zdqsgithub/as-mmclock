#!/usr/bin/env python3
"""Build the v27 processed-supplement candidate manifest.

This is a planning/preflight manifest builder. It discovers official GEO
processed supplement URLs for new P0/P1 candidates, records expected sizes when
available, and marks files already present locally. It does not download.
"""
from __future__ import annotations

import argparse
import html.parser
import json
import re
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")
DEFAULT_OUT = ROOT / "results" / "v27_data_repair_execution" / "processed_candidates"
DEFAULT_DOWNLOAD_ROOT = ROOT / "raw_downloads" / "geo_supplements_v27"
EXISTING_ROOTS = [ROOT / "raw_downloads" / "geo_supplements_v11_3", ROOT / "raw_downloads"]

DATASETS = {
    "GSE225166": {
        "priority": "P0_processed_parser_pilot",
        "allowed_patterns": [r"GSE225166_RAW\.tar$", r"\.cov\.gz$"],
        "reason": "mouse sparse/single-cell low-coverage age methylation; processed COV tar first",
    },
    "GSE225173": {
        "priority": "P0_processed_parser_pilot",
        "allowed_patterns": [r"GSE225173_RAW\.tar$", r"\.cov\.gz$"],
        "reason": "sibling GEO record for GSE225166 methylation supplements; track for duplicate-aware parser pilot",
    },
    "GSE233734": {
        "priority": "P0_processed_parser_pilot",
        "allowed_patterns": [r"GSE233734_RAW\.tar$", r"\.bedGraph(?:\.gz)?$", r"\.bedgraph(?:\.gz)?$"],
        "reason": "102 colon RRBS aging/stage trajectory; processed BEDGRAPH first",
    },
    "GSE304754": {
        "priority": "P1_processed_preflight_then_raw_pilot_if_needed",
        "allowed_patterns": [r"GSE304754_RAW\.tar$", r"\.cov\.gz$", r"\.bedGraph(?:\.gz)?$", r"\.bedgraph(?:\.gz)?$"],
        "reason": "hippocampus aging/cognitive rank candidate; use processed if present before raw ENA pilot",
    },
}


class LinkParser(html.parser.HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.hrefs: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag.lower() != "a":
            return
        for key, value in attrs:
            if key.lower() == "href" and value:
                self.hrefs.append(value)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def geo_series_prefix(dataset: str) -> str:
    digits = dataset.replace("GSE", "")
    return f"GSE{digits[:-3]}nnn"


def geo_sample_supplement_url(filename: str) -> str:
    sample_id = filename.split("_", 1)[0].split(".", 1)[0]
    digits = sample_id.replace("GSM", "")
    if not sample_id.startswith("GSM") or not digits.isdigit() or len(digits) <= 3:
        return ""
    group = f"GSM{digits[:-3]}nnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{group}/{sample_id}/suppl/{filename}"


def fetch_text(url: str, timeout: int = 45) -> tuple[str, str]:
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "as-mmclock-v27/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as response:
            return response.read().decode("utf-8", errors="replace"), ""
    except Exception as exc:  # noqa: BLE001 - report all preflight failures
        return "", str(exc)[:500]


def parse_filelist(dataset: str, text: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for line in text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        # GEO filelist usually has: archive/file, name, timestamp, size, type.
        name = ""
        size = 0
        file_type = ""
        for part in parts:
            if re.search(r"\.(tar|gz|txt|bedGraph|bedgraph|cov)(?:\.gz)?$", part):
                name = part.strip()
            if re.fullmatch(r"\d+(?:\.\d+)?(?:e[+-]?\d+)?", part.strip(), flags=re.IGNORECASE):
                try:
                    size = int(float(part))
                except ValueError:
                    pass
        if len(parts) >= 5:
            file_type = parts[-1].strip()
        if name:
            rows.append({"dataset": dataset, "supplement_name": name, "expected_size_bytes": size, "file_type": file_type, "source": "filelist"})
    return rows


def parse_html_listing(dataset: str, text: str) -> list[dict[str, Any]]:
    parser = LinkParser()
    parser.feed(text)
    rows = []
    for href in parser.hrefs:
        name = href.split("/")[-1]
        if not name or name in {"../", "."}:
            continue
        if re.search(r"\.(tar|gz|txt|bedGraph|bedgraph|cov)(?:\.gz)?$", name):
            rows.append({"dataset": dataset, "supplement_name": name, "expected_size_bytes": 0, "file_type": "", "source": "html_listing"})
    return rows


def discover_dataset(dataset: str) -> tuple[list[dict[str, Any]], str]:
    prefix = geo_series_prefix(dataset)
    base = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{dataset}/suppl"
    filelist_url = f"{base}/filelist.txt"
    text, error = fetch_text(filelist_url)
    if text:
        rows = parse_filelist(dataset, text)
        if rows:
            for row in rows:
                row["inventory_url"] = filelist_url
            return rows, ""
    html_text, html_error = fetch_text(base + "/")
    if html_text:
        rows = parse_html_listing(dataset, html_text)
        for row in rows:
            row["inventory_url"] = base + "/"
        return rows, ""
    return [], error or html_error or "no_supplement_listing"


def allowed(dataset: str, filename: str) -> bool:
    patterns = DATASETS[dataset]["allowed_patterns"]
    return any(re.search(pattern, filename, flags=re.IGNORECASE) for pattern in patterns)


def find_existing(filename: str, expected_size: int) -> tuple[str, str]:
    for root in EXISTING_ROOTS:
        if not root.exists():
            continue
        for path in root.rglob(filename):
            if path.is_file():
                if expected_size <= 0 or path.stat().st_size == expected_size:
                    return str(path), "already_verified_existing_path" if expected_size > 0 else "already_present_size_unknown"
    return "", "planned"


def build_manifest(download_root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, Any]] = []
    skipped: list[dict[str, Any]] = []
    row_id = 0
    for dataset, spec in DATASETS.items():
        discovered, error = discover_dataset(dataset)
        if error:
            skipped.append({"dataset": dataset, "reason": "discover_failed", "error": error})
        for item in discovered:
            filename = str(item["supplement_name"])
            if not allowed(dataset, filename):
                skipped.append({"dataset": dataset, "supplement_name": filename, "reason": "not_allowed_for_v27"})
                continue
            expected = int(item.get("expected_size_bytes") or 0)
            prefix = geo_series_prefix(dataset)
            official_url = geo_sample_supplement_url(filename) or f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{dataset}/suppl/{filename}"
            existing_path, status = find_existing(filename, expected)
            local_path = existing_path or str(download_root / dataset / filename)
            rows.append(
                {
                    "row_id": row_id,
                    "dataset": dataset,
                    "priority": spec["priority"],
                    "reason": spec["reason"],
                    "supplement_name": filename,
                    "file_type": item.get("file_type", ""),
                    "expected_size_bytes": expected,
                    "official_url": official_url,
                    "inventory_url": item.get("inventory_url", ""),
                    "local_path": local_path,
                    "status": status,
                    "download_allowed": True,
                }
            )
            row_id += 1
    return pd.DataFrame(rows), pd.DataFrame(skipped)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--download-root", type=Path, default=DEFAULT_DOWNLOAD_ROOT)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.download_root.mkdir(parents=True, exist_ok=True)

    manifest, skipped = build_manifest(args.download_root)
    manifest_path = args.out_dir / "v27_processed_candidate_download_manifest.csv"
    skipped_path = args.out_dir / "v27_processed_candidate_skipped.csv"
    manifest.to_csv(manifest_path, index=False)
    skipped.to_csv(skipped_path, index=False)
    planned = manifest[manifest["status"].eq("planned")] if not manifest.empty else pd.DataFrame()
    summary = {
        "timestamp": utc_now(),
        "status": "completed",
        "manifest_path": str(manifest_path),
        "skipped_path": str(skipped_path),
        "n_manifest_rows": int(len(manifest)),
        "n_planned_rows": int(len(planned)),
        "planned_known_bytes": int(planned["expected_size_bytes"].sum()) if not planned.empty else 0,
        "datasets": sorted(manifest["dataset"].unique().tolist()) if not manifest.empty else [],
        "download_started": False,
    }
    (args.out_dir / "v27_processed_candidate_manifest_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

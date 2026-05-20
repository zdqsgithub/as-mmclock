#!/usr/bin/env python3
"""Build the v11.3 opportunistic processed supplement download backlog.

This script does not download data. It turns the existing official GEO
inventory into a manifest for the watchdog downloader. The queue is limited to
processed methylation supplements that are useful as auxiliary/diagnostic data;
it intentionally excludes known P4, superseries duplicates, already integrated
headline inputs, and raw-only FASTQ/SRA work.
"""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")
DEFAULT_OUT = ROOT / "results" / "download_backlog_v11_3"
SUPPLEMENTS = ROOT / "metadata" / "geo_old_age_candidate_supplements.csv"
INVENTORY = ROOT / "metadata" / "geo_old_age_candidate_inventory.csv"
DEFAULT_DOWNLOAD_ROOT = ROOT / "raw_downloads" / "geo_supplements_v11_3"
EXISTING_DOWNLOAD_ROOT = ROOT / "raw_downloads"

DOWNLOAD_DATASETS = {
    "GSE175410": {
        "tier": "P2_diagnostic",
        "reason": "24mo skeletal_muscle processed COV; useful auxiliary muscle aging parser/overlap data",
        "allowed_types": {"COV"},
        "headline_allowed": False,
    },
    "GSE281602": {
        "tier": "P3_celltype_auxiliary",
        "reason": "exact 4mo/28mo heart labels but cardiomyocyte/cell-type specific; auxiliary only",
        "allowed_types": {"COV"},
        "headline_allowed": False,
    },
    "GSE129712": {
        "tier": "P2_diagnostic",
        "reason": "26mo intestine TXT processed methylation; schema smoke/auxiliary only",
        "allowed_types": {"TXT"},
        "headline_allowed": False,
    },
    "GSE286302": {
        "tier": "P3_diagnostic",
        "reason": "old lung/liver/muscle processed COV; no exact age_weeks, diagnostic only",
        "allowed_types": {"COV"},
        "headline_allowed": False,
    },
    "GSE92486": {
        "tier": "P2_intervention_auxiliary",
        "reason": "liver DR/AL old-age processed TXT/COV-like files; intervention auxiliary",
        "allowed_types": {"TXT"},
        "headline_allowed": False,
    },
    "GSE224442": {
        "tier": "P2_intervention_auxiliary",
        "reason": "blood/liver parabiosis processed COV; intervention auxiliary",
        "allowed_types": {"COV"},
        "headline_allowed": False,
    },
}

SKIP_DATASETS = {
    "GSE213628": "already integrated in v8; do not redownload",
    "GSE213723": "superseries duplicate; use GSE213628 subseries",
    "GSE233879": "HSC/cell-type context and not target old tissue",
    "GSE276335": "HSC/cell-type context and low priority",
    "GSE47815": "HSC/cell-type context and unknown raw schema",
    "GSE231658": "skin/BW schema; does not address target brain_cortex/heart/lung gap",
    "GSE266961": "large mixed-species/low-age non-target dataset",
    "GSE295059": "organoid/in-vitro context",
    "GSE304754": "no parseable processed methylation supplement detected in current inventory",
    "GSE312263": "young hippocampus only and no old target age support",
    "GSE313655": "no parseable processed methylation supplement detected in current inventory",
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def geo_sample_supplement_url(sample_id: str, filename: str) -> str:
    digits = sample_id.replace("GSM", "")
    if not sample_id.startswith("GSM") or not digits.isdigit() or len(digits) <= 3:
        return ""
    group = f"GSM{digits[:-3]}nnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{group}/{sample_id}/suppl/{filename}"


def sample_id_from_filename(filename: str) -> str:
    token = filename.split("_", 1)[0].split(".", 1)[0]
    return token if token.startswith("GSM") and token[3:].isdigit() else ""


def normalized_size(row: dict) -> int:
    for key in ("remote_size_bytes", "filelist_size_bytes", "size"):
        value = row.get(key)
        if pd.notna(value):
            try:
                return int(float(value))
            except (TypeError, ValueError):
                continue
    return 0


def find_existing_file(filename: str, expected_size: int) -> tuple[str, str]:
    if not filename or expected_size <= 0:
        return "", ""
    for path in EXISTING_DOWNLOAD_ROOT.rglob(filename):
        if path.is_file() and path.stat().st_size == expected_size:
            return str(path), "already_verified_existing_path"
    return "", ""


def load_inventory_notes() -> dict[str, dict]:
    inventory = pd.read_csv(INVENTORY)
    return {
        str(row["dataset"]): row
        for row in inventory.to_dict(orient="records")
        if pd.notna(row.get("dataset"))
    }


def build_manifest(download_root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    supp = pd.read_csv(SUPPLEMENTS)
    notes = load_inventory_notes()
    rows: list[dict] = []
    skipped: list[dict] = []
    row_id = 0

    for dataset, spec in DOWNLOAD_DATASETS.items():
        ds = supp[
            (supp["dataset"].astype(str) == dataset)
            & (supp["archive_or_file"].astype(str).str.lower() == "file")
        ].copy()
        if ds.empty:
            skipped.append({"dataset": dataset, "reason": "no_file_level_supplements_in_inventory"})
            continue
        for item in ds.to_dict(orient="records"):
            file_type = str(item.get("filelist_type", "")).upper()
            filename = str(item.get("supplement_name", ""))
            if file_type not in spec["allowed_types"]:
                continue
            if not filename.endswith((".gz", ".txt.gz", ".cov.gz", ".bedGraph.gz", ".bedgraph.gz")):
                continue
            expected = normalized_size(item)
            if expected <= 0:
                skipped.append({"dataset": dataset, "supplement_name": filename, "reason": "missing_size"})
                continue
            sample_id = sample_id_from_filename(filename)
            inventory_url = str(item.get("supplement_url") or "")
            sample_url = geo_sample_supplement_url(sample_id, filename) if sample_id else ""
            url = sample_url or inventory_url
            existing_path, existing_status = find_existing_file(filename, expected)
            if existing_path:
                local_path = existing_path
                initial_status = existing_status
            else:
                local_path = str(download_root / dataset / filename)
                initial_status = "planned"
            inv = notes.get(dataset, {})
            rows.append(
                {
                    "row_id": row_id,
                    "dataset": dataset,
                    "sample_id": sample_id,
                    "supplement_name": filename,
                    "file_type": file_type,
                    "schema_guess": inv.get("processed_schema_guess", ""),
                    "tier": spec["tier"],
                    "reason": spec["reason"],
                    "headline_allowed": spec["headline_allowed"],
                    "download_allowed": True,
                    "expected_size_bytes": expected,
                    "official_url": url,
                    "sample_level_url": sample_url,
                    "inventory_url": inventory_url,
                    "local_path": local_path,
                    "status": initial_status,
                    "retry_count": 0,
                    "last_error": "",
                    "priority_tier_source": inv.get("priority_tier", ""),
                    "recommended_action_source": inv.get("recommended_action", ""),
                    "tissues_source": inv.get("tissues", ""),
                    "age_values_weeks_source": inv.get("age_values_weeks", ""),
                }
            )
            row_id += 1

    for dataset, reason in SKIP_DATASETS.items():
        inv = notes.get(dataset, {})
        skipped.append(
            {
                "dataset": dataset,
                "reason": reason,
                "priority_tier_source": inv.get("priority_tier", ""),
                "recommended_action_source": inv.get("recommended_action", ""),
                "preferred_size_bytes": inv.get("preferred_size_bytes", ""),
                "supplement_types": inv.get("supplement_types", ""),
            }
        )

    manifest = pd.DataFrame(rows)
    skipped_df = pd.DataFrame(skipped)
    return manifest, skipped_df


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out_dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--download_root", type=Path, default=DEFAULT_DOWNLOAD_ROOT)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.download_root.mkdir(parents=True, exist_ok=True)

    manifest, skipped = build_manifest(args.download_root)
    manifest_path = args.out_dir / "download_manifest.csv"
    skipped_path = args.out_dir / "download_skipped_manifest.csv"
    manifest.to_csv(manifest_path, index=False)
    skipped.to_csv(skipped_path, index=False)

    planned = manifest[manifest["status"] == "planned"]
    verified = manifest[manifest["status"].astype(str).str.startswith("already_verified")]
    summary = {
        "timestamp": utc_now(),
        "status": "completed",
        "manifest_path": str(manifest_path),
        "skipped_manifest_path": str(skipped_path),
        "download_root": str(args.download_root),
        "n_rows_total": int(len(manifest)),
        "n_rows_planned": int(len(planned)),
        "n_rows_already_verified": int(len(verified)),
        "planned_bytes": int(planned["expected_size_bytes"].sum()) if not planned.empty else 0,
        "already_verified_bytes": int(verified["expected_size_bytes"].sum()) if not verified.empty else 0,
        "datasets_planned": sorted(planned["dataset"].unique().tolist()) if not planned.empty else [],
        "datasets_all": sorted(manifest["dataset"].unique().tolist()) if not manifest.empty else [],
    }
    (args.out_dir / "download_manifest_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

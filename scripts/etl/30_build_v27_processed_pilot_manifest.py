#!/usr/bin/env python3
"""Build a duplicate-aware v27 processed parser pilot manifest.

The manifest is compatible with scripts/etl/20_download_backlog_watchdog.py.
It selects a bounded P0 pilot instead of downloading every newly discovered
supplement immediately, and it reuses already verified v11.3 processed files
for secondary candidates.
"""
from __future__ import annotations

import argparse
import json
import re
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")
DEFAULT_CANDIDATES = ROOT / "results" / "v27_data_repair_execution" / "processed_candidates" / "v27_processed_candidate_download_manifest.csv"
DEFAULT_V11 = ROOT / "results" / "download_backlog_v11_3" / "download_manifest.csv"
DEFAULT_OUT = ROOT / "results" / "v27_data_repair_execution" / "processed_pilot"

SECONDARY_DATASETS = {"GSE175410", "GSE224442", "GSE92486"}
TISSUE_RE = re.compile(r"_(Liver|Heart|Cortex|Lung)_", re.IGNORECASE)
AGE_RE = re.compile(r"_(\d+)_M_")


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def head_size(url: str, timeout: int = 30) -> int:
    if not url:
        return 0
    try:
        req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "as-mmclock-v27/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as response:
            value = response.headers.get("Content-Length", "")
        return int(value) if value.isdigit() else 0
    except Exception:
        return 0


def infer_tissue(filename: str) -> str:
    match = TISSUE_RE.search(filename)
    return match.group(1).lower() if match else ""


def infer_age_weeks(filename: str) -> float | None:
    match = AGE_RE.search(filename)
    if not match:
        return None
    return float(match.group(1))


def hard_status(row: pd.Series) -> tuple[str, int]:
    path = Path(str(row.get("local_path", "")))
    expected = int(row.get("expected_size_bytes") or 0)
    actual = path.stat().st_size if path.exists() else 0
    if path.exists() and (expected <= 0 or actual == expected):
        return "already_verified_existing_path", actual
    return str(row.get("status", "planned") or "planned"), actual


def select_gse225166(df: pd.DataFrame, max_files: int) -> pd.DataFrame:
    cov = df[(df["dataset"].eq("GSE225166")) & (df["supplement_name"].str.endswith(".cov.gz"))].copy()
    if cov.empty:
        return cov
    cov["pilot_tissue"] = cov["supplement_name"].map(infer_tissue)
    cov["pilot_age_weeks"] = cov["supplement_name"].map(infer_age_weeks)
    selected = []
    # First cover each tissue/age pair, then fill by file size/name.
    for _, row in cov.sort_values(["pilot_age_weeks", "pilot_tissue", "supplement_name"]).groupby(
        ["pilot_age_weeks", "pilot_tissue"],
        dropna=False,
        sort=True,
    ).head(1).iterrows():
        selected.append(row)
        if len(selected) >= max_files:
            break
    if len(selected) < max_files:
        used = {row["supplement_name"] for row in selected}
        for _, row in cov.sort_values(["expected_size_bytes", "supplement_name"], ascending=[False, True]).iterrows():
            if row["supplement_name"] in used:
                continue
            selected.append(row)
            if len(selected) >= max_files:
                break
    return pd.DataFrame(selected)


def normalize_candidate_rows(df: pd.DataFrame, out_dir: Path) -> pd.DataFrame:
    rows = df.copy()
    if "sample_id" not in rows.columns:
        rows["sample_id"] = rows["supplement_name"].astype(str).str.extract(r"(GSM\d+)", expand=False).fillna("")
    rows["schema_guess"] = rows["supplement_name"].astype(str).map(
        lambda name: "bismark_cov_per_sample_file" if name.endswith(".cov.gz") else "bedgraph_pair_or_combined"
    )
    rows["tier"] = rows.get("priority", "v27_processed_pilot")
    rows["headline_allowed"] = False
    rows["retry_count"] = 0
    rows["last_error"] = ""
    rows["pilot_source"] = "v27_new_processed_candidates"
    return rows


def build_manifest(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    candidates = pd.read_csv(args.candidates)
    selected_parts = []
    skipped = []

    gse225166 = select_gse225166(candidates, args.gse225166_files)
    selected_parts.append(gse225166)
    duplicate_rows = candidates[candidates["dataset"].eq("GSE225173")].copy()
    if not duplicate_rows.empty:
        skipped.extend(
            {
                "dataset": row.dataset,
                "supplement_name": row.supplement_name,
                "reason": "duplicate_filename_with_gse225166_pilot_or_sibling_record",
            }
            for row in duplicate_rows.itertuples(index=False)
            if str(row.supplement_name) in set(candidates[candidates["dataset"].eq("GSE225166")]["supplement_name"].astype(str))
        )

    p0_bed = candidates[candidates["dataset"].isin(["GSE233734"])].copy()
    selected_parts.append(p0_bed)
    if args.include_secondary_new:
        selected_parts.append(candidates[candidates["dataset"].isin(["GSE304754"])].copy())

    selected = pd.concat([part for part in selected_parts if not part.empty], ignore_index=True) if selected_parts else pd.DataFrame()
    selected = normalize_candidate_rows(selected, args.out_dir) if not selected.empty else selected

    for idx, row in selected.iterrows():
        expected = int(row.get("expected_size_bytes") or 0)
        if expected <= 0:
            selected.at[idx, "expected_size_bytes"] = head_size(str(row.get("official_url", "")))
    if not selected.empty:
        selected["download_allowed"] = selected["expected_size_bytes"].astype(int).gt(0)
        selected.loc[~selected["download_allowed"], "status"] = "blocked_missing_content_length"

    if args.include_v11_secondary and args.v11_manifest.exists():
        v11 = pd.read_csv(args.v11_manifest)
        v11 = v11[v11["dataset"].isin(SECONDARY_DATASETS)].copy()
        for idx, row in v11.iterrows():
            status, actual = hard_status(row)
            v11.at[idx, "status"] = status
            v11.at[idx, "actual_size_bytes"] = actual
        v11["pilot_source"] = "v11_3_verified_secondary_candidate"
        selected = pd.concat([selected, v11], ignore_index=True, sort=False)

    if selected.empty:
        manifest = pd.DataFrame()
    else:
        selected = selected.drop_duplicates(["dataset", "supplement_name"], keep="first").copy()
        selected["row_id"] = range(len(selected))
        required = [
            "row_id",
            "dataset",
            "sample_id",
            "supplement_name",
            "file_type",
            "schema_guess",
            "tier",
            "reason",
            "headline_allowed",
            "download_allowed",
            "expected_size_bytes",
            "official_url",
            "local_path",
            "status",
            "retry_count",
            "last_error",
            "pilot_source",
        ]
        for col in required:
            if col not in selected.columns:
                selected[col] = ""
        manifest = selected[required].copy()
        manifest["download_allowed"] = manifest["download_allowed"].fillna(True).astype(bool)
        manifest.loc[manifest["expected_size_bytes"].fillna(0).astype(int).le(0), "download_allowed"] = False

    skipped_df = pd.DataFrame(skipped)
    planned = manifest[manifest["download_allowed"] & manifest["status"].eq("planned")] if not manifest.empty else pd.DataFrame()
    verified = manifest[manifest["status"].astype(str).str.startswith("already_verified")] if not manifest.empty else pd.DataFrame()
    summary = {
        "timestamp": utc_now(),
        "status": "completed",
        "n_manifest_rows": int(len(manifest)),
        "n_planned_download_rows": int(len(planned)),
        "planned_download_bytes": int(planned["expected_size_bytes"].sum()) if not planned.empty else 0,
        "n_already_verified_rows": int(len(verified)),
        "already_verified_bytes": int(verified["expected_size_bytes"].sum()) if not verified.empty else 0,
        "datasets": sorted(manifest["dataset"].unique().tolist()) if not manifest.empty else [],
        "download_started": False,
    }
    return manifest, skipped_df, summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, default=DEFAULT_CANDIDATES)
    parser.add_argument("--v11-manifest", type=Path, default=DEFAULT_V11)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--gse225166-files", type=int, default=32)
    parser.add_argument("--include-secondary-new", action="store_true")
    parser.add_argument("--include-v11-secondary", action="store_true", default=True)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    manifest, skipped, summary = build_manifest(args)
    manifest_path = args.out_dir / "v27_processed_pilot_download_manifest.csv"
    skipped_path = args.out_dir / "v27_processed_pilot_skipped.csv"
    summary_path = args.out_dir / "v27_processed_pilot_manifest_summary.json"
    manifest.to_csv(manifest_path, index=False)
    skipped.to_csv(skipped_path, index=False)
    summary["manifest_path"] = str(manifest_path)
    summary["skipped_path"] = str(skipped_path)
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

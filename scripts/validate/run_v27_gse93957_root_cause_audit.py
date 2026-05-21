#!/usr/bin/env python3
"""Audit why the current GSE93957 raw ETL branch fails usable-region gates."""
from __future__ import annotations

import argparse
import gzip
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")
DEFAULT_ROOT = Path("/data/mouse_methyl/processed_v21_raw_etl/GSE93957")
DEFAULT_OUT = ROOT / "results" / "v27_data_repair_execution" / "gse93957_root_cause"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def extract_number(pattern: str, text: str, *, integer: bool = False) -> float | int | None:
    match = re.search(pattern, text, flags=re.IGNORECASE)
    if not match:
        return None
    value = match.group(1).replace(",", "")
    try:
        return int(value) if integer else float(value)
    except ValueError:
        return None


def read_text(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8", errors="replace")
    except Exception:
        return ""


def count_cov_rows(path: Path, max_rows: int = 1_000_000) -> dict[str, Any]:
    rows = 0
    coverage_ge_1 = 0
    coverage_ge_5 = 0
    coverage_values: list[float] = []
    try:
        with gzip.open(path, "rt", errors="replace") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#") or line.startswith("track"):
                    continue
                rows += 1
                if rows > max_rows:
                    break
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 6:
                    try:
                        total = float(parts[4]) + float(parts[5])
                    except ValueError:
                        continue
                    coverage_values.append(total)
                    if total >= 1:
                        coverage_ge_1 += 1
                    if total >= 5:
                        coverage_ge_5 += 1
    except Exception as exc:
        return {"cov_read_error": str(exc)[:500]}
    return {
        "cov_rows_checked": int(rows),
        "cov_rows_coverage_ge_1": int(coverage_ge_1),
        "cov_rows_coverage_ge_5": int(coverage_ge_5),
        "cov_max_coverage_checked": max(coverage_values) if coverage_values else None,
        "cov_median_coverage_checked": float(pd.Series(coverage_values).median()) if coverage_values else None,
    }


def audit_sample(sample_dir: Path) -> dict[str, Any]:
    sample_id = sample_dir.name
    stats_path = sample_dir / f"{sample_id}_raw_parse_stats.json"
    stats = json.loads(stats_path.read_text(encoding="utf-8")) if stats_path.exists() else {}
    run_dirs = sorted([path for path in (sample_dir / "runs").glob("*") if path.is_dir()]) if (sample_dir / "runs").exists() else []
    report_text = "\n".join(read_text(path) for run in run_dirs for path in run.glob("*_PE_report.txt"))
    split_text = "\n".join(read_text(path) for run in run_dirs for path in run.glob("*splitting_report.txt"))
    cov_paths = [path for run in run_dirs for path in run.glob("*.bismark.cov.gz")]
    cov_stats = count_cov_rows(cov_paths[0]) if cov_paths else {}

    analysed = extract_number(r"Sequence pairs analysed in total:\s*([\d,]+)", report_text, integer=True)
    unique = extract_number(r"unique best hit:\s*([\d,]+)", report_text, integer=True)
    mapping_eff = extract_number(r"Mapping efficiency:\s*([\d.]+)%", report_text)
    processed_lines = extract_number(r"Processed\s+([\d,]+)\s+lines in total", split_text, integer=True)
    return {
        "sample_id": sample_id,
        "run_accessions": ";".join(run.name for run in run_dirs),
        "bismark_sequence_pairs": analysed,
        "bismark_unique_alignments": unique,
        "bismark_mapping_efficiency_pct": mapping_eff,
        "extractor_processed_lines": processed_lines,
        "parse_rows_total": stats.get("rows_total"),
        "parse_rows_primary_autosomes": stats.get("rows_parseable_primary_autosomes"),
        "parse_rows_pass_coverage": stats.get("rows_pass_coverage"),
        "parse_n_regions": stats.get("n_regions"),
        "parse_beta_min": stats.get("beta_min"),
        "parse_beta_max": stats.get("beta_max"),
        "n_cov_files": len(cov_paths),
        "first_cov_path": str(cov_paths[0]) if cov_paths else "",
        **cov_stats,
    }


def classify(rows: pd.DataFrame) -> tuple[str, str]:
    if rows.empty:
        return "blocked_no_completed_samples", "No completed GSE93957 ETL sample directories with parse stats were found."
    median_mapping = rows["bismark_mapping_efficiency_pct"].dropna().median()
    median_cov_ge5 = rows["cov_rows_coverage_ge_5"].dropna().median() if "cov_rows_coverage_ge_5" in rows else 0
    median_regions = rows["parse_n_regions"].dropna().median()
    if pd.notna(median_mapping) and median_mapping < 1.0:
        return (
            "blocked_alignment_or_pretrim_issue",
            "Bismark mapping efficiency is near zero, so more download is not the next fix. Recheck trimming/adapters/library settings and assembly before more ETL.",
        )
    if pd.notna(median_cov_ge5) and median_cov_ge5 < 50_000:
        return (
            "blocked_low_coverage_after_alignment",
            "Alignment exists but COV coverage above 5 is too sparse for 5kb region matrix gates.",
        )
    if pd.notna(median_regions) and median_regions < 50_000:
        return (
            "blocked_parser_or_region_gate_issue",
            "Coverage is not obviously absent, but parsed region counts fail the sample region gate.",
        )
    return "candidate_for_parser_rebuild", "Existing outputs may be sufficient for a parser-only matrix rebuild."


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--etl-root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    sample_dirs = sorted([path for path in args.etl_root.glob("GSM*") if (path / f"{path.name}_raw_parse_stats.json").exists()])
    rows = pd.DataFrame([audit_sample(path) for path in sample_dirs])
    status, recommendation = classify(rows)
    audit_path = args.out_dir / "gse93957_completed_sample_audit.csv"
    rows.to_csv(audit_path, index=False)
    summary = {
        "timestamp": utc_now(),
        "status": status,
        "recommendation": recommendation,
        "n_samples_audited": int(len(rows)),
        "median_mapping_efficiency_pct": None if rows.empty else rows["bismark_mapping_efficiency_pct"].dropna().median(),
        "median_parse_regions": None if rows.empty else rows["parse_n_regions"].dropna().median(),
        "median_cov_rows_coverage_ge_5": None if rows.empty or "cov_rows_coverage_ge_5" not in rows else rows["cov_rows_coverage_ge_5"].dropna().median(),
        "audit_path": str(audit_path),
        "policy": "no_more_download_or_etl_until_root_cause_fixed",
    }
    summary_path = args.out_dir / "gse93957_root_cause_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True, default=str), encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True, default=str))


if __name__ == "__main__":
    main()

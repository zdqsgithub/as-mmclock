#!/usr/bin/env python3
"""Verify v11 raw pilot candidates with SRA RunInfo, GEO GSM, and BioSample metadata.

This script does not download FASTQ files. It builds a sample/run verification
manifest that decides whether a later 2-3 run FASTQ/Bismark pilot is justified.
"""
from __future__ import annotations

import csv
import gzip
import io
import json
import re
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")
OUT_DIR = ROOT / "results" / "ralph_v11_data_strategy"
OUT = OUT_DIR / "raw_biosample_verification_manifest.csv"
SUMMARY_OUT = OUT_DIR / "raw_biosample_verification_summary.csv"
STATE_OUT = OUT_DIR / "raw_biosample_decision_state.json"
RAW_PROJECT_SUMMARY = OUT_DIR / "raw_project_rfc_summary.csv"

CANDIDATES = [
    {"candidate": "GSE286302_PRJNA1208582_SRP556367", "study": "SRP556367", "linked_gse": "GSE286302"},
    {"candidate": "GSE281602_PRJNA1184779_SRP544506", "study": "SRP544506", "linked_gse": "GSE281602"},
]
MAX_RUNS_PER_CANDIDATE = 12

TARGET_PATTERNS = {
    "brain_cortex": ["brain cortex", "cortex"],
    "heart": ["heart", "cardiac"],
    "lung": ["lung"],
}
BLOCKER_TERMS = [
    "single-cell",
    "single cell",
    "organoid",
    "cell line",
    "in vitro",
    "ipcrtag",
    "itag",
    "rna-seq",
    "cell type",
    "cardiomyocyte",
    "cardiomyocytes",
]
AGE_RE = re.compile(r"(?P<value>\d+(?:\.\d+)?)\s*[- ]?(?P<unit>day|days|week|weeks|wk|wks|mo|mos|month|months|m|year|years|yr|yrs)", re.I)


def fetch_text(url: str, timeout: int = 45) -> str:
    request = urllib.request.Request(url, headers={"User-Agent": "mouse-methyl-v11-biosample-verify/1.0"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        data = response.read()
        if url.endswith(".gz"):
            data = gzip.decompress(data)
        return data.decode("utf-8", errors="replace")


def runinfo(study: str) -> pd.DataFrame:
    url = f"https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc={study}"
    text = fetch_text(url)
    return pd.read_csv(io.StringIO(text))


def priority_runs(candidate: dict) -> set[str]:
    if not RAW_PROJECT_SUMMARY.exists():
        return set()
    summary = pd.read_csv(RAW_PROJECT_SUMMARY)
    keys = {candidate["study"], candidate["linked_gse"]}
    rows = summary[summary["project_key"].astype(str).isin(keys)]
    runs = set()
    for value in rows.get("raw_accession_examples", pd.Series(dtype=str)).fillna(""):
        runs.update(re.findall(r"SRR\d+", str(value)))
    return runs


def geo_sample_text(gsm: str) -> str:
    if not str(gsm).startswith("GSM"):
        return ""
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}&targ=self&view=full&form=text"
    try:
        return fetch_text(url)
    except Exception:
        return ""


def biosample_text(samn: str) -> str:
    if not str(samn).startswith("SAMN"):
        return ""
    url = f"https://www.ncbi.nlm.nih.gov/biosample/{samn}?report=full&format=text"
    try:
        return fetch_text(url)
    except Exception:
        return ""


def parse_geo_field(text: str, field: str) -> str:
    values = re.findall(rf"^!Sample_{re.escape(field)}\s*=\s*(.*)$", text, flags=re.M)
    return " | ".join(value.strip() for value in values if value.strip())


def parse_age_weeks(text: str) -> tuple[float | None, str]:
    match = AGE_RE.search(text or "")
    if not match:
        return None, ""
    value = float(match.group("value"))
    unit = match.group("unit").lower()
    if unit.startswith("day"):
        weeks = value / 7.0
    elif unit.startswith(("week", "wk")):
        weeks = value
    elif unit in {"m", "mo", "mos", "month", "months"}:
        weeks = value * 30.42 / 7.0
    else:
        weeks = value * 365.0 / 7.0
    return weeks, match.group(0)


def age_group(text: str) -> str:
    lower = (text or "").lower()
    if re.search(r"\b(aged|old)\b", lower):
        return "old"
    if re.search(r"\b(young|juvenile)\b", lower):
        return "young"
    return ""


def detect_tissues(text: str) -> list[str]:
    lower = (text or "").lower()
    hits = []
    for label, needles in TARGET_PATTERNS.items():
        if any(needle in lower for needle in needles):
            hits.append(label)
    return sorted(set(hits))


def has_blocker(text: str) -> str:
    lower = (text or "").lower()
    return ";".join(term for term in BLOCKER_TERMS if term in lower)


def main() -> None:
    rows = []
    bio_cache: dict[str, str] = {}
    geo_cache: dict[str, str] = {}
    for candidate in CANDIDATES:
        frame = runinfo(candidate["study"])
        selected_runs = priority_runs(candidate)
        if selected_runs:
            frame = frame[frame["Run"].astype(str).isin(selected_runs)].copy()
        if len(frame) > MAX_RUNS_PER_CANDIDATE:
            frame = frame.head(MAX_RUNS_PER_CANDIDATE).copy()
        for _, run in frame.iterrows():
            gsm = str(run.get("LibraryName", "") or run.get("SampleName", ""))
            samn = str(run.get("BioSample", ""))
            if gsm not in geo_cache:
                geo_cache[gsm] = geo_sample_text(gsm)
                time.sleep(0.1)
            if samn not in bio_cache:
                bio_cache[samn] = biosample_text(samn)
                time.sleep(0.1)
            geo_text = geo_cache.get(gsm, "")
            bio_text = bio_cache.get(samn, "")
            sample_title = parse_geo_field(geo_text, "title")
            source_name = parse_geo_field(geo_text, "source_name_ch1")
            characteristics = parse_geo_field(geo_text, "characteristics_ch1")
            series = parse_geo_field(geo_text, "series_id")
            combined = " ".join(
                str(value)
                for value in [
                    run.get("Run", ""),
                    run.get("Experiment", ""),
                    run.get("LibraryName", ""),
                    run.get("LibraryStrategy", ""),
                    run.get("LibrarySelection", ""),
                    run.get("SampleName", ""),
                    sample_title,
                    source_name,
                    characteristics,
                    bio_text,
                ]
            )
            age_weeks, raw_age_token = parse_age_weeks(combined)
            tissues = detect_tissues(combined)
            blockers = has_blocker(combined)
            assay_pass = str(run.get("LibraryStrategy", "")).lower() == "bisulfite-seq"
            target_tissue_pass = bool(tissues)
            age_exact_pass = age_weeks is not None
            group = age_group(combined)
            age_group_pass = group in {"old", "young"}
            size_mb = float(run.get("size_MB", 0) or 0)
            pilot_role = ""
            if assay_pass and target_tissue_pass and not blockers and age_exact_pass:
                pilot_role = "chronological_pilot_candidate"
            elif assay_pass and target_tissue_pass and not blockers and age_group_pass:
                pilot_role = "diagnostic_age_group_pilot_candidate"
            rows.append(
                {
                    "candidate": candidate["candidate"],
                    "linked_gse_expected": candidate["linked_gse"],
                    "run": run.get("Run", ""),
                    "experiment": run.get("Experiment", ""),
                    "sample": run.get("Sample", ""),
                    "biosample": samn,
                    "gsm": gsm,
                    "geo_series_id": series,
                    "library_strategy": run.get("LibraryStrategy", ""),
                    "library_selection": run.get("LibrarySelection", ""),
                    "library_layout": run.get("LibraryLayout", ""),
                    "size_mb": size_mb,
                    "sample_title": sample_title,
                    "source_name": source_name,
                    "characteristics": characteristics,
                    "biosample_text_compact": " ".join(bio_text.split())[:500],
                    "age_weeks": round(age_weeks, 3) if age_weeks is not None else "",
                    "raw_age_token": raw_age_token,
                    "age_group": group,
                    "tissue": ";".join(tissues),
                    "assay_pass": bool(assay_pass),
                    "target_tissue_pass": bool(target_tissue_pass),
                    "age_exact_pass": bool(age_exact_pass),
                    "age_group_pass": bool(age_group_pass),
                    "blockers": blockers,
                    "pilot_role": pilot_role,
                    "include_in_minimal_pilot": bool(pilot_role),
                }
            )
    manifest = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(OUT, index=False, quoting=csv.QUOTE_MINIMAL)
    summary = (
        manifest.groupby(["candidate", "pilot_role"], dropna=False)
        .agg(
            n_runs=("run", "nunique"),
            n_samples=("gsm", "nunique"),
            total_size_mb=("size_mb", "sum"),
            tissues=("tissue", lambda values: ";".join(sorted(set(v for value in values for v in str(value).split(";") if v)))),
            age_groups=("age_group", lambda values: ";".join(sorted(set(v for v in values if str(v))))),
            exact_age_n=("age_exact_pass", "sum"),
        )
        .reset_index()
    )
    summary.to_csv(SUMMARY_OUT, index=False)
    chronological_n = int(manifest["pilot_role"].eq("chronological_pilot_candidate").sum())
    diagnostic_n = int(manifest["pilot_role"].eq("diagnostic_age_group_pilot_candidate").sum())
    if chronological_n >= 2:
        status = "pilot_ready_chronological"
        decision = "v11_1_chronological_minimal_etl_pilot_allowed_after_resource_check"
        next_action = "Prepare a 2-3 run FASTQ/Bismark pilot manifest; do not start full ETL."
    elif diagnostic_n >= 2:
        status = "diagnostic_only"
        decision = "v11_1_only_diagnostic_old_lung_pilot_candidate"
        next_action = "Optional 2-3 run diagnostic pilot for GSE286302 old lung only; not a headline chronological benchmark."
    else:
        status = "failure"
        decision = "v11_1_no_verified_minimal_etl_candidate"
        next_action = "Stop raw-pilot path and move to external data generation/acquisition or benchmark redefinition."
    state = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "loop_version": "v11.1",
        "status": status,
        "loop_decision": decision,
        "metrics": {
            "n_manifest_rows": int(len(manifest)),
            "n_chronological_pilot_rows": chronological_n,
            "n_diagnostic_pilot_rows": diagnostic_n,
            "n_candidates": int(manifest["candidate"].nunique()),
        },
        "next_action": next_action,
    }
    STATE_OUT.write_text(json.dumps(state, indent=2, sort_keys=True), encoding="utf-8")
    print(summary.to_string(index=False))
    print(json.dumps(state, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

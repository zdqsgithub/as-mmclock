#!/usr/bin/env python3
"""v14 public-data rescue RALPH controller.

This script implements the safe public-data branch after v13 Route C. It uses
official GEO/E-Utils/SOFT/FTP metadata plus optional SRA/ENA run-file metadata
to classify old mouse methylome candidates. It does not download FASTQ files,
run Bismark, train models, or start autoresearch.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import importlib.util
import json
import math
import re
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")
OUT_DIR = ROOT / "results" / "ralph_v14_public_data_rescue"
DOC_REPORT = ROOT / "doc" / "20_analysis" / "43_20260519_v14_public_data_rescue_report.md"

COMMON_REGION_GATE = 50_000
OLD_THRESHOLD_WEEKS = 104.0
TARGET_TISSUES = {"brain_cortex", "cortex", "heart", "lung"}
PREVIOUSLY_INTEGRATED_OR_EVALUATED = {
    "GSE120137",
    "GSE80672",
    "GSE93957",
    "GSE121141",
    "GSE60012",
    "GSE213628",
}

SEED_ACCESSIONS = [
    "GSE286302",
    "GSE281602",
    "GSE304754",
    "GSE292803",
    "GSE313770",
    "GSE225166",
    "GSE83947",
    "GSE134398",
]

SEARCH_TERMS = [
    'Mus musculus[ORGN] AND "Methylation profiling by high throughput sequencing" AND (RRBS OR WGBS OR bisulfite) AND ("brain cortex" OR cortex OR heart OR lung OR hippocampus) AND (aging OR aged OR old OR "24 month" OR "28 month") AND gse[ETYP]',
    'Mus musculus[ORGN] AND "DNA methylation" AND ("brain cortex" OR cortex OR heart OR lung) AND ("24 month" OR "28 month" OR aged OR old) AND gse[ETYP]',
    'Mus musculus[ORGN] AND "reduced representation bisulfite" AND (heart OR lung OR cortex OR hippocampus) AND aging AND gse[ETYP]',
]

OFFICIAL_SOURCE_DOCS = {
    "geo_download": "https://www.ncbi.nlm.nih.gov/geo/info/download.html",
    "geo_programmatic_access": "https://www.ncbi.nlm.nih.gov/geo/info/geo_paccess.html",
    "geo_soft": "https://www.ncbi.nlm.nih.gov/geo/info/soft.html",
    "sra_download": "https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/",
    "ena_file_reports": "https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/file-reports.html",
}

MANUAL_NOTES = {
    "GSE286302": "old lung RRBS/COV diagnostic; sample titles use young/aged labels but exact age_weeks is not sample-specific",
    "GSE281602": "28mo heart RRBS but FACS cardiomyocyte nuclei; cell-type context blocks bulk heart headline",
    "GSE304754": "28mo hippocampus/liver/colon/fecal RRBS; brain-adjacent diagnostic, not cortex headline",
    "GSE292803": "single-cell multi-omic aging brain atlas; not bulk and max age appears below old104 threshold",
    "GSE313770": "male single-cell multi-omic aging brain atlas; not bulk and max age appears below old104 threshold",
    "GSE225166": "single-cell/low-coverage/iTAG age methylation; adapter audit/reference only",
}

HARD_BLOCKER_PATTERNS = {
    "single_cell_or_single_nucleus": [
        "single-cell",
        "single cell",
        "single-nucleus",
        "single nucleus",
        "snmc",
        "snm3c",
        "single cell multi-omic",
    ],
    "targeted_cell_type_or_sorted": [
        "cardiomyocyte",
        "macrophage",
        "microglia",
        "hsc",
        "hematopoietic stem",
        "facs",
        "sorted",
        "cd45",
        "siglecf",
        "nuclei",
    ],
    "low_coverage_or_itag": ["low-coverage", "low coverage", "itag", "ipcrtag"],
    "organoid_or_in_vitro": ["organoid", "in vitro"],
    "superseries_or_untraceable": ["superseries"],
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def clean(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    return str(value)


def import_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def discovery_module():
    return import_module(ROOT / "scripts" / "etl" / "16_geo_old_age_candidate_discovery.py", "geo_old_age_candidate_discovery_v14")


def read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str), encoding="utf-8")


def append_jsonl(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(payload, sort_keys=True, default=str) + "\n")


def md_table(df: pd.DataFrame, max_rows: int | None = None) -> str:
    if df.empty:
        return "No rows."
    view = df.copy()
    if max_rows is not None:
        view = view.head(max_rows)
    for col in view.columns:
        view[col] = view[col].map(lambda value: "" if pd.isna(value) else str(value))
    header = "| " + " | ".join(view.columns) + " |"
    sep = "| " + " | ".join(["---"] * len(view.columns)) + " |"
    rows = ["| " + " | ".join(row) + " |" for row in view.to_numpy(dtype=str)]
    return "\n".join([header, sep, *rows])


def fetch_text(url: str, timeout: int = 45) -> tuple[str, str]:
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "mouse-methyl-v14-public-data-rescue/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as response:
            data = response.read()
        if url.endswith(".gz"):
            data = gzip.decompress(data)
        return data.decode("utf-8", errors="replace"), ""
    except Exception as exc:
        return "", str(exc)[:500]


def extract_bioprojects_from_soft(module, accession: str) -> tuple[list[str], str]:
    soft = module.cached_soft_text(accession, refresh=False)
    if not soft.text:
        return [], soft.error or "soft_unavailable"
    projects = sorted(set(re.findall(r"PRJ[DEN]A?\d+", soft.text)))
    if not projects:
        projects = sorted(set(re.findall(r"PRJ[DEN]\d+", soft.text)))
    return projects, ""


def ena_file_report(accession: str) -> tuple[list[dict[str, str]], str]:
    fields = "run_accession,fastq_ftp,fastq_md5,fastq_bytes,library_strategy,library_layout,scientific_name,sample_accession"
    params = {
        "accession": accession,
        "result": "read_run",
        "fields": fields,
        "format": "tsv",
    }
    url = "https://www.ebi.ac.uk/ena/portal/api/filereport?" + urllib.parse.urlencode(params)
    text, error = fetch_text(url, timeout=60)
    if error:
        return [], error
    reader = csv.DictReader(text.splitlines(), delimiter="\t")
    return list(reader), ""


def sra_runinfo(accession: str) -> tuple[list[dict[str, str]], str]:
    url = "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?" + urllib.parse.urlencode({"acc": accession})
    text, error = fetch_text(url, timeout=60)
    if error:
        return [], error
    reader = csv.DictReader(text.splitlines())
    return list(reader), ""


def numeric(value: Any) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return parsed if math.isfinite(parsed) else float("nan")


def parse_age_weeks_v14(text: str) -> float:
    patterns = [
        re.compile(
            r"age\s*[:=]\s*(?P<value>\d+(?:\.\d+)?)\s*(?P<unit>d|day|days|w|wk|wks|week|weeks|m/o|mo|mos|month|months|m|y|yr|yrs|year|years)",
            re.IGNORECASE,
        ),
        re.compile(
            r"(?<![A-Za-z0-9])(?P<value>\d+(?:\.\d+)?)[\s-]*(?P<unit>d|day|days|w|wk|wks|week|weeks|m/o|mo|mos|month|months|m|y|yr|yrs|year|years)(?![A-Za-z0-9])",
            re.IGNORECASE,
        ),
    ]
    candidates: list[float] = []
    for pattern in patterns:
        for match in pattern.finditer(text):
            value = float(match.group("value"))
            unit = match.group("unit").lower()
            if unit in {"d", "day", "days"}:
                weeks = value / 7.0
            elif unit in {"w", "wk", "wks", "week", "weeks"}:
                weeks = value
            elif unit in {"m/o", "mo", "mos", "month", "months", "m"}:
                weeks = value * 30.42 / 7.0
            elif unit in {"y", "yr", "yrs", "year", "years"}:
                weeks = value * 365.0 / 7.0
            else:
                continue
            if 0 <= weeks <= 220:
                candidates.append(weeks)
    return max(candidates) if candidates else float("nan")


def canonical_tissue(text: str) -> str:
    lowered = f" {text.lower().replace('-', ' ').replace('_', ' ')} "
    if " brain cortex " in lowered or " cortex " in lowered:
        return "brain_cortex"
    if " hippocampus " in lowered or " hippocampal " in lowered:
        return "brain_adjacent_hippocampus"
    if " pvn " in lowered or " paraventricular " in lowered or " hypothalamus " in lowered:
        return "brain_adjacent_hypothalamus"
    if " brain " in lowered:
        return "brain_other"
    if " heart " in lowered or " cardiac " in lowered or " cardiomyocyte " in lowered:
        return "heart"
    if " lung " in lowered:
        return "lung"
    if " liver " in lowered:
        return "liver"
    if " muscle " in lowered or " skeletal muscle " in lowered:
        return "skeletal_muscle"
    if " colon " in lowered or " intestine " in lowered:
        return "intestine"
    return "unknown"


def old_label(text: str) -> bool:
    lowered = f" {text.lower()} "
    return any(token in lowered for token in [" aged", " old", "28-mo", "28 mo", "24-mo", "24 mo", "26-mo", "26 mo", "30-mo", "30 mo"])


def context_blockers(row: dict[str, Any], samples: pd.DataFrame) -> list[str]:
    text = " ".join(clean(row.get(field)) for field in ["dataset", "title", "summary", "overall_design", "blockers"])
    if not samples.empty:
        text += " " + " ".join(
            " ".join(clean(sample.get(field)) for field in ["title", "source_name", "characteristics"])
            for _, sample in samples.head(80).iterrows()
        )
    lowered = text.lower()
    blockers: list[str] = []
    for label, needles in HARD_BLOCKER_PATTERNS.items():
        if any(needle in lowered for needle in needles):
            blockers.append(label)
    return blockers


def summarize_sample_support(samples: pd.DataFrame) -> dict[str, Any]:
    if samples.empty:
        return {
            "n_samples": 0,
            "age_known_n": 0,
            "age_known_fraction": 0.0,
            "old_exact_target_n": 0,
            "old_label_target_n": 0,
            "target_tissue_n": 0,
            "old_exact_target_tissues": "",
            "old_label_target_tissues": "",
            "max_age_weeks": "",
            "tissue_counts": "",
        }
    rows = []
    for _, sample in samples.iterrows():
        blob = " ".join(clean(sample.get(field)) for field in ["title", "source_name", "characteristics"])
        tissue = canonical_tissue(blob)
        age = numeric(sample.get("age_weeks"))
        if not math.isfinite(age):
            age = parse_age_weeks_v14(blob)
        rows.append(
            {
                "tissue": tissue,
                "age_weeks": age,
                "age_known": math.isfinite(age),
                "old_exact": math.isfinite(age) and age >= OLD_THRESHOLD_WEEKS,
                "old_label": old_label(blob),
                "target": tissue in TARGET_TISSUES,
            }
        )
    data = pd.DataFrame(rows)
    n = len(data)
    age_known_n = int(data["age_known"].sum())
    old_exact_target = data[data["old_exact"] & data["target"]]
    old_label_target = data[data["old_label"] & data["target"]]
    target = data[data["target"]]
    ages = pd.to_numeric(data["age_weeks"], errors="coerce").dropna()
    tissue_counts = Counter(data["tissue"].astype(str))
    return {
        "n_samples": int(n),
        "age_known_n": age_known_n,
        "age_known_fraction": round(age_known_n / n, 4) if n else 0.0,
        "old_exact_target_n": int(len(old_exact_target)),
        "old_label_target_n": int(len(old_label_target)),
        "target_tissue_n": int(len(target)),
        "old_exact_target_tissues": ";".join(sorted(set(old_exact_target["tissue"].astype(str)))),
        "old_label_target_tissues": ";".join(sorted(set(old_label_target["tissue"].astype(str)))),
        "max_age_weeks": round(float(ages.max()), 3) if not ages.empty else "",
        "tissue_counts": ";".join(f"{key}:{value}" for key, value in tissue_counts.most_common()),
    }


def is_parseable_processed(schema: str, supplements: pd.DataFrame) -> bool:
    text = f"{schema} " + " ".join(supplements.get("supplement_name", pd.Series(dtype=str)).astype(str).tolist())
    lowered = text.lower()
    return any(token in lowered for token in ["cov", "bismark.cov", ".cov", "bedgraph"])


def has_sample_processed_files(supplements: pd.DataFrame) -> bool:
    if supplements.empty:
        return False
    archive = supplements.get("archive_or_file", pd.Series(dtype=str)).astype(str).str.lower()
    file_type = supplements.get("filelist_type", pd.Series(dtype=str)).astype(str).str.upper()
    names = supplements.get("supplement_name", pd.Series(dtype=str)).astype(str)
    return bool((archive.eq("file") & (file_type.isin(["COV", "TXT", "BEDGRAPH"]) | names.str.contains("CX_report|bismark|cov|bedgraph", case=False, regex=True))).any())


def select_processed_pilot_rows(dataset: str, samples: pd.DataFrame, supplements: pd.DataFrame, tier: str) -> list[dict[str, Any]]:
    if samples.empty or supplements.empty:
        return []
    file_rows = supplements[
        supplements.get("archive_or_file", pd.Series(dtype=str)).astype(str).str.lower().eq("file")
    ].copy()
    if file_rows.empty:
        return []
    selected: list[dict[str, Any]] = []
    for _, sample in samples.iterrows():
        sample_id = clean(sample.get("sample_id"))
        if not sample_id:
            continue
        blob = " ".join(clean(sample.get(field)) for field in ["title", "source_name", "characteristics"])
        tissue = canonical_tissue(blob)
        age = numeric(sample.get("age_weeks"))
        if not math.isfinite(age):
            age = parse_age_weeks_v14(blob)
        is_target = tissue in TARGET_TISSUES
        is_old = math.isfinite(age) and age >= OLD_THRESHOLD_WEEKS
        if not (is_target and (is_old or old_label(blob))):
            continue
        matches = file_rows[file_rows.get("supplement_name", pd.Series(dtype=str)).astype(str).str.contains(sample_id, regex=False)]
        for _, file_row in matches.iterrows():
            selected.append(
                {
                    "dataset": dataset,
                    "bioproject": "",
                    "run_accession": "",
                    "sample_accession": sample_id,
                    "library_strategy": "processed_methylation_supplement",
                    "library_layout": "",
                    "scientific_name": "Mus musculus",
                    "fastq_bytes": "",
                    "fastq_ftp": "",
                    "fastq_md5": "",
                    "processed_supplement_name": file_row.get("supplement_name", ""),
                    "processed_supplement_url": file_row.get("supplement_url", ""),
                    "processed_supplement_bytes": file_row.get("filelist_size_bytes", ""),
                    "processed_supplement_type": file_row.get("filelist_type", ""),
                    "sample_tissue": tissue,
                    "sample_age_weeks": round(float(age), 3) if math.isfinite(age) else "",
                    "selected_for_download": False,
                    "download_authorized": False,
                    "pilot_authorization_state": "metadata_only_no_download",
                    "pilot_role": "candidate_processed_adapter_pilot" if tier == "P3_raw_pilot" else "auxiliary_or_blocked_reference",
                }
            )
    selected = sorted(
        selected,
        key=lambda row: (
            0 if clean(row.get("sample_tissue")) == "brain_cortex" else 1 if clean(row.get("sample_tissue")) == "heart" else 2,
            numeric(row.get("processed_supplement_bytes")) if math.isfinite(numeric(row.get("processed_supplement_bytes"))) else 10**18,
            clean(row.get("sample_accession")),
        ),
    )
    return selected[:3]


def local_conversion_state(dataset: str) -> dict[str, Any]:
    paths = [
        ROOT / "results" / "v11_4_conversions" / dataset / f"v11_x_{dataset}_conversion_gate_state.json",
        ROOT / "results" / "v11_4_conversions" / dataset / f"v11_4_{dataset.lower()}_conversion_smoke_state.json",
        ROOT / "results" / "v11_2_gse286302_cov_pilot" / "v11_2_decision_state.json" if dataset == "GSE286302" else Path("__missing__"),
    ]
    for path in paths:
        payload = read_json(path)
        if not payload:
            continue
        qc = payload.get("qc") or payload.get("overlap") or {}
        manifest = payload.get("conversion_manifest") or {}
        return {
            "local_conversion_state_path": str(path.relative_to(ROOT)),
            "local_conversion_status": payload.get("final_status") or payload.get("status") or manifest.get("status") or "",
            "local_headline_allowed": bool(payload.get("headline_allowed", False)),
            "local_reason": payload.get("reason") or payload.get("reason_not_headline") or payload.get("decision") or "",
            "local_common_regions": qc.get("common_regions_with_v8_2") or qc.get("common_regions") or "",
            "local_common_region_gate": bool(qc.get("common_region_gate_50000", False) or (numeric(qc.get("common_regions_with_v8_2") or qc.get("common_regions")) >= COMMON_REGION_GATE)),
            "local_beta_min": qc.get("beta_min", ""),
            "local_beta_max": qc.get("beta_max", ""),
            "local_n_samples": qc.get("n_samples") or manifest.get("n_samples_parsed") or "",
            "local_matrix_path": manifest.get("region_matrix_path") or qc.get("pilot_region_matrix") or "",
        }
    return {
        "local_conversion_state_path": "",
        "local_conversion_status": "",
        "local_headline_allowed": False,
        "local_reason": "",
        "local_common_regions": "",
        "local_common_region_gate": False,
        "local_beta_min": "",
        "local_beta_max": "",
        "local_n_samples": "",
        "local_matrix_path": "",
    }


def classify_candidate(row: dict[str, Any]) -> tuple[str, str, str, bool]:
    age_pass = float(row.get("age_known_fraction") or 0.0) >= 0.95
    exact_target = int(row.get("old_exact_target_n") or 0) > 0
    label_target = int(row.get("old_label_target_n") or 0) > 0
    bulk_pass = not clean(row.get("context_blockers"))
    parseable = bool(row.get("parseable_processed_schema"))
    processed_files = bool(row.get("has_sample_processed_files"))
    common_regions = numeric(row.get("local_common_regions"))
    common_pass = math.isfinite(common_regions) and common_regions >= COMMON_REGION_GATE
    dataset = clean(row.get("dataset"))

    if dataset in PREVIOUSLY_INTEGRATED_OR_EVALUATED:
        reason = "previously_integrated_or_evaluated_dataset"
        if clean(row.get("local_reason")):
            reason += ";" + clean(row.get("local_reason"))
        return "P2_auxiliary", "keep_as_existing_evidence_not_new_v14_headline", reason, False

    if age_pass and exact_target and bulk_pass and parseable and common_pass:
        return "P1_headline", "ready_for_fixed_benchmark_pending_explicit_training_approval", "", True
    if age_pass and exact_target and bulk_pass and processed_files and not parseable:
        return "P3_raw_pilot", "prepare_2_3_sample_processed_adapter_smoke_no_download_in_this_run", "processed_supplement_schema_unknown_requires_adapter_smoke", False
    if age_pass and exact_target and bulk_pass and not parseable:
        return "P3_raw_pilot", "prepare_2_3_sample_raw_pilot_manifest_no_download_in_this_run", "no_parseable_processed_matrix", False
    if age_pass and exact_target and bulk_pass and parseable and not common_pass:
        return "P3_raw_pilot", "processed_or_raw_pilot_matrix_gate_required", "common_regions_not_passed_or_unknown", False
    if exact_target or label_target or row.get("dataset") in {"GSE286302", "GSE281602", "GSE304754"}:
        reasons = []
        if not age_pass:
            reasons.append("sample_specific_exact_age_coverage_lt_95pct")
        if not bulk_pass:
            reasons.append("non_bulk_or_cell_context")
        if not exact_target and label_target:
            reasons.append("old_age_group_label_not_exact_age")
        if parseable and not common_pass and clean(row.get("local_common_regions")):
            reasons.append("common_region_gate_failed")
        return "P2_auxiliary", "keep_for_diagnostic_or_adapter_evidence_not_headline", ";".join(reasons), False
    reasons = []
    if not exact_target:
        reasons.append("no_exact_old_target_tissue_samples")
    if not label_target:
        reasons.append("no_old_target_tissue_label")
    if not bulk_pass:
        reasons.append("blocked_context")
    return "P4_blocked", "do_not_promote_to_headline", ";".join(reasons), False


def select_candidate_accessions(module, accessions: list[str], no_eutils_search: bool, retmax: int, max_candidates: int) -> tuple[list[str], dict[str, list[str]], dict[str, dict[str, Any]]]:
    source_terms: dict[str, list[str]] = {acc: ["manual_seed_v14"] for acc in accessions}
    summaries: dict[str, dict[str, Any]] = {}
    selected = set(accessions)
    if not no_eutils_search:
        found, found_terms, found_summaries = module.discover_accessions(SEARCH_TERMS, retmax=retmax)
        selected.update(found)
        summaries.update(found_summaries)
        for accession, terms in found_terms.items():
            source_terms.setdefault(accession, []).extend(terms)
    for accession in list(selected):
        if accession not in summaries:
            try:
                summaries[accession] = module.eutils_summary_for_accession(accession)
            except Exception as exc:
                summaries[accession] = {"dataset": accession, "eutils_error": str(exc)[:300]}
    ranked = sorted(
        selected,
        key=lambda value: (
            0 if value in accessions else 1,
            -module.rough_esummary_score(value, summaries.get(value, {}), source_terms.get(value, [])),
            int(value.removeprefix("GSE")) if value.removeprefix("GSE").isdigit() else 999_999_999,
        ),
    )
    return ranked[:max_candidates], source_terms, summaries


def build_network_and_pilot_rows(
    dataset: str,
    bioprojects: list[str],
    tier: str,
    out_dir: Path,
    network_log: Path,
    query_sra: bool,
    query_ena: bool,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    run_rows: list[dict[str, Any]] = []
    summary: dict[str, Any] = {
        "bioprojects": ";".join(bioprojects),
        "ena_n_runs": 0,
        "ena_total_fastq_bytes": 0,
        "sra_n_runs": 0,
        "sra_total_size_mb": 0.0,
        "network_errors": "",
    }
    errors: list[str] = []
    for project in bioprojects[:3]:
        if query_ena:
            rows, error = ena_file_report(project)
            append_jsonl(
                network_log,
                {
                    "timestamp": utc_now(),
                    "dataset": dataset,
                    "blocker_type": "ena_file_report",
                    "attempted_url": f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={project}&result=read_run",
                    "source_doc": OFFICIAL_SOURCE_DOCS["ena_file_reports"],
                    "resolution": "metadata_only_run_file_report",
                    "applied_rule": "no_fastq_download",
                    "error": error,
                },
            )
            if error:
                errors.append(f"ENA:{project}:{error}")
            else:
                summary["ena_n_runs"] += len(rows)
                for run in rows[:3]:
                    bytes_text = clean(run.get("fastq_bytes"))
                    total = 0
                    for item in bytes_text.split(";"):
                        try:
                            total += int(item)
                        except ValueError:
                            pass
                    summary["ena_total_fastq_bytes"] += total
                    run_rows.append(
                        {
                            "dataset": dataset,
                            "bioproject": project,
                            "run_accession": run.get("run_accession", ""),
                            "sample_accession": run.get("sample_accession", ""),
                            "library_strategy": run.get("library_strategy", ""),
                            "library_layout": run.get("library_layout", ""),
                            "scientific_name": run.get("scientific_name", ""),
                            "fastq_bytes": total,
                            "fastq_ftp": run.get("fastq_ftp", ""),
                            "fastq_md5": run.get("fastq_md5", ""),
                            "selected_for_download": False,
                            "download_authorized": False,
                            "pilot_authorization_state": "metadata_only_no_download",
                            "pilot_role": "candidate_raw_pilot" if tier == "P3_raw_pilot" else "auxiliary_or_blocked_reference",
                        }
                    )
        if query_sra:
            rows, error = sra_runinfo(project)
            append_jsonl(
                network_log,
                {
                    "timestamp": utc_now(),
                    "dataset": dataset,
                    "blocker_type": "sra_runinfo",
                    "attempted_url": f"https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc={project}",
                    "source_doc": OFFICIAL_SOURCE_DOCS["sra_download"],
                    "resolution": "metadata_only_runinfo",
                    "applied_rule": "no_fastq_download",
                    "error": error,
                },
            )
            if error:
                errors.append(f"SRA:{project}:{error}")
            else:
                summary["sra_n_runs"] += len(rows)
                for run in rows:
                    summary["sra_total_size_mb"] += numeric(run.get("size_MB")) if math.isfinite(numeric(run.get("size_MB"))) else 0.0
    summary["network_errors"] = " | ".join(errors)
    if not run_rows:
        run_rows.append(
            {
                "dataset": dataset,
                "bioproject": ";".join(bioprojects),
                "run_accession": "",
                "sample_accession": "",
                "library_strategy": "",
                "library_layout": "",
                "scientific_name": "",
                "fastq_bytes": "",
                "fastq_ftp": "",
                "fastq_md5": "",
                "selected_for_download": False,
                "download_authorized": False,
                "pilot_authorization_state": "no_run_metadata_or_not_queried",
                "pilot_role": "candidate_raw_pilot" if tier == "P3_raw_pilot" else "auxiliary_or_blocked_reference",
            }
        )
    return run_rows, summary


def write_report(out_dir: Path, doc_report: Path, refresh: pd.DataFrame, gates: pd.DataFrame, state: dict[str, Any]) -> None:
    gate_cols = [
        "dataset",
        "tier",
        "age_known_fraction",
        "old_exact_target_n",
        "old_label_target_n",
        "old_exact_target_tissues",
        "context_blockers",
        "local_common_regions",
        "headline_allowed",
        "gate_status",
        "blocker_type",
    ]
    gate_view = gates[[col for col in gate_cols if col in gates.columns]].copy()
    tier_counts = gates["tier"].value_counts().rename_axis("tier").reset_index(name="count") if not gates.empty else pd.DataFrame()
    lines = [
        "# v14 Public-Data Rescue Report",
        "",
        f"Date: {utc_now()}",
        "",
        "## Summary",
        "",
        "v14 ran a public-data RALPH Research/Audit pass against official GEO/E-Utils/SOFT/FTP metadata and SRA/ENA run-file metadata. This run did not download FASTQ, run Bismark, train models, or start autoresearch.",
        "",
        f"- Decision: `{state['loop_decision']}`",
        f"- Next action: {state['next_action']}",
        f"- Headline candidates allowed now: {state['metrics']['n_headline_allowed']}",
        f"- P3 raw-pilot candidates: {state['metrics']['n_p3_raw_pilot']}",
        "",
        "## Tier Counts",
        "",
        md_table(tier_counts),
        "",
        "## Candidate Gate Table",
        "",
        md_table(gate_view, max_rows=80),
        "",
        "## Interpretation",
        "",
        "- `GSE286302` remains useful as old-lung processed-COV compatibility evidence, but it lacks sample-specific exact age and is lung-only.",
        "- `GSE281602` has exact old heart samples, but the context is cardiomyocyte/cell-type specific and its common-region overlap failed the project gate.",
        "- `GSE304754` is brain-adjacent hippocampus rather than brain_cortex/cortex headline support.",
        "- Single-cell/low-coverage/iTAG brain or methylation-age datasets are reference or adapter-audit material only.",
        "",
        "## Guardrails",
        "",
        "- Training authorized: `false`",
        "- Download authorized: `false` in this run; only run metadata was queried.",
        "- Bismark authorized: `false`",
        "- Autoresearch authorized: `false`",
        "- Human clock CpG mapping: forbidden",
        "- Dummy AUC: forbidden",
        "",
        "## Outputs",
        "",
        f"- `{(out_dir / 'candidate_refresh_table.csv').relative_to(ROOT)}`",
        f"- `{(out_dir / 'candidate_gate_table.csv').relative_to(ROOT)}`",
        f"- `{(out_dir / 'pilot_run_manifest.csv').relative_to(ROOT)}`",
        f"- `{(out_dir / 'network_resolution_log.jsonl').relative_to(ROOT)}`",
        f"- `{(out_dir / 'ralph_decision_state.json').relative_to(ROOT)}`",
        "",
        "## Sources",
        "",
        "- GEO Download / Programmatic Access / SOFT official documentation.",
        "- NCBI SRA download documentation and ENA file report API documentation.",
        "- GEO accession pages and cached SOFT metadata for candidate accessions.",
    ]
    out_report = out_dir / "v14_public_data_rescue_report.md"
    out_report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    doc_report.parent.mkdir(parents=True, exist_ok=True)
    doc_report.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", default=str(OUT_DIR))
    parser.add_argument("--doc-report", default=str(DOC_REPORT))
    parser.add_argument("--accessions", default=",".join(SEED_ACCESSIONS))
    parser.add_argument("--retmax", type=int, default=25)
    parser.add_argument("--max-candidates", type=int, default=24)
    parser.add_argument("--no-eutils-search", action="store_true")
    parser.add_argument("--refresh-soft", action="store_true")
    parser.add_argument("--no-head", action="store_true")
    parser.add_argument("--skip-ena", action="store_true")
    parser.add_argument("--skip-sra", action="store_true")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    network_log = out_dir / "network_resolution_log.jsonl"
    if network_log.exists():
        network_log.unlink()

    module = discovery_module()
    accessions = [item.strip() for item in args.accessions.split(",") if item.strip()]
    selected, source_terms, summaries = select_candidate_accessions(
        module,
        accessions=accessions,
        no_eutils_search=args.no_eutils_search,
        retmax=args.retmax,
        max_candidates=args.max_candidates,
    )

    refresh_rows: list[dict[str, Any]] = []
    gate_rows: list[dict[str, Any]] = []
    pilot_rows: list[dict[str, Any]] = []
    smoke_rows: list[dict[str, Any]] = []

    for accession in selected:
        try:
            row, supplements, samples = module.build_candidate(
                accession,
                source_terms=source_terms.get(accession, ["v14_search"]),
                eutils_summary=summaries.get(accession),
                refresh_soft=args.refresh_soft,
                do_head=not args.no_head,
            )
        except Exception as exc:
            row = {
                "dataset": accession,
                "geo_url": module.geo_url(accession),
                "title": "",
                "summary": "",
                "overall_design": "",
                "processed_schema_guess": "",
                "blockers": str(exc)[:500],
                "recommended_action": "inspect_manually",
            }
            supplements = []
            samples = []
        supplement_df = pd.DataFrame(supplements)
        sample_df = pd.DataFrame(samples)
        support = summarize_sample_support(sample_df)
        blockers = context_blockers(row, sample_df)
        local = local_conversion_state(accession)
        bioprojects, bioproject_error = extract_bioprojects_from_soft(module, accession)
        parseable = is_parseable_processed(clean(row.get("processed_schema_guess")), supplement_df)
        processed_files = has_sample_processed_files(supplement_df)
        candidate = {
            **row,
            **support,
            **local,
            "dataset": accession,
            "manual_note": MANUAL_NOTES.get(accession, ""),
            "bioprojects": ";".join(bioprojects),
            "bioproject_error": bioproject_error,
            "context_blockers": ";".join(blockers),
            "bulk_context_pass": not blockers,
            "parseable_processed_schema": parseable,
            "has_sample_processed_files": processed_files,
            "official_geo_url": module.geo_url(accession),
            "v14_source_terms": " || ".join(source_terms.get(accession, [])),
        }
        tier, gate_status, blocker_type, headline_allowed = classify_candidate(candidate)
        candidate.update(
            {
                "tier": tier,
                "gate_status": gate_status,
                "blocker_type": blocker_type or candidate.get("context_blockers") or candidate.get("local_reason") or candidate.get("blockers", ""),
                "headline_allowed": headline_allowed,
                "download_authorized": False,
                "training_authorized": False,
                "bismark_authorized": False,
                "autoresearch_authorized": False,
            }
        )
        run_rows, network_summary = build_network_and_pilot_rows(
            accession,
            bioprojects,
            tier,
            out_dir,
            network_log,
            query_sra=not args.skip_sra,
            query_ena=not args.skip_ena,
        )
        processed_rows = select_processed_pilot_rows(accession, sample_df, supplement_df, tier)
        candidate.update(network_summary)
        refresh_rows.append(candidate)
        gate_rows.append(
            {
                "dataset": accession,
                "tier": tier,
                "gate_status": gate_status,
                "blocker_type": candidate["blocker_type"],
                "recommended_action": gate_status,
                "headline_allowed": headline_allowed,
                "sample_specific_age_pass": float(candidate.get("age_known_fraction") or 0.0) >= 0.95,
                "target_tissue_pass": int(candidate.get("old_exact_target_n") or 0) > 0,
                "bulk_context_pass": candidate["bulk_context_pass"],
                "schema_smoke_pass": bool(candidate.get("local_conversion_status")),
                "assembly_traceable": bool(candidate.get("local_common_regions")),
                "common_regions_pass": bool(candidate.get("local_common_region_gate")),
                "common_regions_estimate": candidate.get("local_common_regions", ""),
                "age_known_fraction": candidate.get("age_known_fraction", ""),
                "old_exact_target_n": candidate.get("old_exact_target_n", ""),
                "old_label_target_n": candidate.get("old_label_target_n", ""),
                "old_exact_target_tissues": candidate.get("old_exact_target_tissues", ""),
                "old_label_target_tissues": candidate.get("old_label_target_tissues", ""),
                "target_tissue_n": candidate.get("target_tissue_n", ""),
                "tissue_counts": candidate.get("tissue_counts", ""),
                "context_blockers": candidate.get("context_blockers", ""),
                "processed_schema_guess": candidate.get("processed_schema_guess", ""),
                "has_sample_processed_files": candidate.get("has_sample_processed_files", False),
                "preferred_url": candidate.get("preferred_url", ""),
                "bioprojects": candidate.get("bioprojects", ""),
                "ena_n_runs": candidate.get("ena_n_runs", 0),
                "sra_n_runs": candidate.get("sra_n_runs", 0),
                "local_conversion_status": candidate.get("local_conversion_status", ""),
                "local_conversion_state_path": candidate.get("local_conversion_state_path", ""),
                "local_reason": candidate.get("local_reason", ""),
            }
        )
        pilot_rows.extend(processed_rows or run_rows)
        smoke_rows.append(
            {
                "dataset": accession,
                "local_conversion_status": candidate.get("local_conversion_status", ""),
                "local_conversion_state_path": candidate.get("local_conversion_state_path", ""),
                "local_common_regions": candidate.get("local_common_regions", ""),
                "local_common_region_gate": candidate.get("local_common_region_gate", False),
                "local_beta_min": candidate.get("local_beta_min", ""),
                "local_beta_max": candidate.get("local_beta_max", ""),
                "local_n_samples": candidate.get("local_n_samples", ""),
                "local_matrix_path": candidate.get("local_matrix_path", ""),
                "schema_smoke_pass": bool(candidate.get("local_conversion_status")),
            }
        )
        time.sleep(0.05)

    refresh_df = pd.DataFrame(refresh_rows)
    gate_df = pd.DataFrame(gate_rows)
    pilot_df = pd.DataFrame(pilot_rows)
    smoke_df = pd.DataFrame(smoke_rows)

    refresh_df.to_csv(out_dir / "candidate_refresh_table.csv", index=False)
    gate_df.to_csv(out_dir / "candidate_gate_table.csv", index=False)
    pilot_df.to_csv(out_dir / "pilot_run_manifest.csv", index=False)
    smoke_df.to_csv(out_dir / "candidate_smoke_manifest.csv", index=False)

    n_headline = int(gate_df["headline_allowed"].astype(bool).sum()) if not gate_df.empty else 0
    n_p3 = int(gate_df["tier"].astype(str).eq("P3_raw_pilot").sum()) if not gate_df.empty else 0
    n_p2 = int(gate_df["tier"].astype(str).eq("P2_auxiliary").sum()) if not gate_df.empty else 0
    if n_headline:
        decision = "v14_public_data_headline_candidate_ready_pending_explicit_training_approval"
        next_action = "Run fixed RALPH Learn benchmark only after explicit training approval."
    elif n_p3:
        decision = "v14_public_data_raw_pilot_candidate_found_no_download_started"
        next_action = "Review pilot_run_manifest and authorize only 2-3 sample raw/processed pilot if the candidate is scientifically acceptable."
    else:
        decision = "v14_no_public_headline_or_raw_pilot_candidate_keep_v13_route_c"
        next_action = "Keep v13 support-covered benchmark boundary; public candidates remain auxiliary/blocked."

    state = {
        "timestamp": utc_now(),
        "loop_version": "v14.0",
        "loop_decision": decision,
        "next_action": next_action,
        "training_authorized": False,
        "download_authorized": False,
        "bismark_authorized": False,
        "autoresearch_authorized": False,
        "human_clock_mapping_authorized": False,
        "dummy_auc_authorized": False,
        "metrics": {
            "n_candidates": int(len(gate_df)),
            "n_headline_allowed": n_headline,
            "n_p1_headline": int(gate_df["tier"].astype(str).eq("P1_headline").sum()) if not gate_df.empty else 0,
            "n_p2_auxiliary": n_p2,
            "n_p3_raw_pilot": n_p3,
            "n_p4_blocked": int(gate_df["tier"].astype(str).eq("P4_blocked").sum()) if not gate_df.empty else 0,
            "common_region_gate": COMMON_REGION_GATE,
            "old_threshold_weeks": OLD_THRESHOLD_WEEKS,
        },
        "official_source_docs": OFFICIAL_SOURCE_DOCS,
        "search_terms": SEARCH_TERMS,
        "seed_accessions": accessions,
        "outputs": {
            "candidate_refresh_table": str((out_dir / "candidate_refresh_table.csv").relative_to(ROOT)),
            "candidate_gate_table": str((out_dir / "candidate_gate_table.csv").relative_to(ROOT)),
            "pilot_run_manifest": str((out_dir / "pilot_run_manifest.csv").relative_to(ROOT)),
            "candidate_smoke_manifest": str((out_dir / "candidate_smoke_manifest.csv").relative_to(ROOT)),
            "network_resolution_log": str((out_dir / "network_resolution_log.jsonl").relative_to(ROOT)),
        },
    }
    write_json(out_dir / "ralph_decision_state.json", state)
    append_jsonl(
        out_dir / "ralph_iteration_log.jsonl",
        {
            "timestamp": utc_now(),
            "phase": "RAP",
            "hypothesis": "public data can repair or falsify old target tissue support gap",
            "action": "official_metadata_refresh_candidate_gate_no_download_no_training",
            "decision": decision,
            "next_action": next_action,
            "metrics": state["metrics"],
        },
    )
    write_report(out_dir, Path(args.doc_report), refresh_df, gate_df, state)
    print(json.dumps(state, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

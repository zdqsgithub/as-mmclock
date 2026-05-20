#!/usr/bin/env python3
"""Refresh Route B manual raw-lead backlog with official run metadata only.

This script queries NCBI RunInfo/SRA XML and ENA read_run metadata for backlog
items that require manual BioSample/RunInfo refresh. It does not download FASTQ,
run Bismark, build matrices, train models, or authorize a pilot.
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

ROOT = Path("/home/zdq-as/mouse_methyl_work")
OUT_DIR = ROOT / "results" / "route_b_candidate_audit"
BACKLOG = OUT_DIR / "route_b_raw_lead_backlog.csv"
REPORT_PATH = ROOT / "doc" / "20_analysis" / "35_20260519_route_b_manual_backlog_runinfo_refresh_report.md"
INDEX_PATH = ROOT / "doc" / "00_meta" / "02_20260519_project_status_index_v13.md"

ROW_OUT = OUT_DIR / "route_b_manual_backlog_runinfo_refresh.csv"
SUMMARY_OUT = OUT_DIR / "route_b_manual_backlog_candidate_summary.csv"
NETWORK_LOG = OUT_DIR / "route_b_manual_backlog_network_log.jsonl"
STATE_OUT = OUT_DIR / "route_b_manual_backlog_decision_state.json"

OLD_WEEKS = 104.0
TARGET_TISSUES = {"brain_cortex", "heart", "lung"}
MANUAL_ACTION = "official_biosample_runinfo_refresh_required_before_any_pilot"

AGE_RE = re.compile(
    r"(?P<value>\d+(?:\.\d+)?)\s*[- ]?"
    r"(?P<unit>day|days|d|week|weeks|wk|wks|w|month|months|mo|mos|m|year|years|yr|yrs|y)\b",
    re.I,
)

ACCESSION_RE = {
    "prjna": re.compile(r"PRJNA\d+", re.I),
    "srp": re.compile(r"SRP\d+", re.I),
    "srr": re.compile(r"SRR\d+", re.I),
    "prjeb": re.compile(r"PRJEB\d+", re.I),
    "erp": re.compile(r"ERP\d+", re.I),
    "err": re.compile(r"ERR\d+", re.I),
}

ASSAY_TERMS = {
    "bisulfite",
    "bisulfite-seq",
    "bs-seq",
    "rrbs",
    "wgbs",
    "methyl-seq",
    "methylseq",
    "em-seq",
    "enzymatic methyl",
}
ASSAY_BLOCKERS = {"rna-seq", "chip-seq", "atac-seq", "rna seq", "rnaseq"}
CONTEXT_BLOCKERS = {
    "single-cell",
    "single cell",
    "single_cell",
    "cell type",
    "cell-type",
    "cardiomyocyte",
    "cardiomyocytes",
    "organoid",
    "in vitro",
    "cell line",
    "ipcrtag",
    "itag",
    "low-coverage",
}


@dataclass
class Group:
    group_id: str
    backend: str
    query_accession: str
    row_count: int
    source_project_keys: list[str]
    raw_accessions: list[str]
    evidence_snippets: list[str]


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str), encoding="utf-8")


def append_network_log(entry: dict[str, Any]) -> None:
    NETWORK_LOG.parent.mkdir(parents=True, exist_ok=True)
    with NETWORK_LOG.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(entry, sort_keys=True, default=str) + "\n")


def fetch_text(url: str, *, endpoint: str, timeout: int = 60) -> str:
    started = datetime.now().isoformat(timespec="seconds")
    request = urllib.request.Request(url, headers={"User-Agent": "mouse-methyl-route-b-refresh/1.0"})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            data = response.read()
            text = data.decode("utf-8", errors="replace")
            append_network_log(
                {
                    "timestamp": started,
                    "endpoint": endpoint,
                    "url": url,
                    "status": "ok",
                    "bytes": len(data),
                    "source_doc": source_doc(endpoint),
                    "resolution": "metadata_only_request_succeeded",
                    "applied_rule": "no_fastq_download",
                }
            )
            return text
    except Exception as exc:  # noqa: BLE001 - network failure must be logged and non-fatal
        status = "error"
        if isinstance(exc, urllib.error.HTTPError):
            status = f"http_{exc.code}"
        append_network_log(
            {
                "timestamp": started,
                "endpoint": endpoint,
                "url": url,
                "status": status,
                "blocker_type": type(exc).__name__,
                "error": str(exc)[:500],
                "source_doc": source_doc(endpoint),
                "resolution": "record_blocker_and_continue",
                "applied_rule": "metadata_refresh_is_non_fatal",
            }
        )
        return ""


def source_doc(endpoint: str) -> str:
    if endpoint.startswith("ncbi"):
        return "NCBI SRA RunInfo/E-Utilities official metadata endpoint"
    if endpoint.startswith("ena"):
        return "ENA Portal API official read_run filereport endpoint"
    return "official_metadata_endpoint"


def extract_accessions(text: str, kind: str) -> list[str]:
    return sorted({item.upper() for item in ACCESSION_RE[kind].findall(text or "")})


def first_accession(text: str, kinds: list[str]) -> str:
    for kind in kinds:
        values = extract_accessions(text, kind)
        if values:
            return values[0]
    return ""


def make_groups(rows: list[dict[str, str]]) -> list[Group]:
    manual_rows = [row for row in rows if row.get("next_action") == MANUAL_ACTION]
    buckets: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in manual_rows:
        text = row.get("raw_accession_examples", "")
        err = first_accession(text, ["err"])
        if err:
            key = ("ena_run", err)
        else:
            key = ("ncbi_group", first_accession(text, ["prjna", "srp", "srr"]) or row.get("project_key", ""))
        buckets[key].append(row)

    groups: list[Group] = []
    for (backend, query_accession), bucket in sorted(buckets.items(), key=lambda item: item[0][1]):
        raw_values: list[str] = []
        snippets: list[str] = []
        keys: list[str] = []
        for row in bucket:
            keys.append(row.get("project_key", ""))
            raw_values.extend([item for item in row.get("raw_accession_examples", "").split(";") if item])
            snippet = row.get("evidence_snippet", "")
            if snippet:
                snippets.append(snippet[:350])
        groups.append(
            Group(
                group_id=query_accession,
                backend=backend,
                query_accession=query_accession,
                row_count=len(bucket),
                source_project_keys=sorted(set(keys)),
                raw_accessions=sorted(set(raw_values)),
                evidence_snippets=snippets[:4],
            )
        )
    return groups


def ncbi_runinfo(accession: str) -> list[dict[str, str]]:
    url = "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=" + urllib.parse.quote(accession)
    text = fetch_text(url, endpoint="ncbi_runinfo")
    if not text.strip() or text.startswith("Error"):
        return []
    try:
        return list(csv.DictReader(io.StringIO(text)))
    except csv.Error:
        append_network_log(
            {
                "timestamp": datetime.now().isoformat(timespec="seconds"),
                "endpoint": "ncbi_runinfo",
                "url": url,
                "status": "schema_mismatch",
                "blocker_type": "csv_parse_failed",
                "source_doc": source_doc("ncbi_runinfo"),
                "resolution": "record_blocker_and_continue",
                "applied_rule": "metadata_refresh_is_non_fatal",
            }
        )
        return []


def ncbi_sra_xml(run: str) -> dict[str, str]:
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=" + urllib.parse.quote(run) + "&retmode=xml"
    text = fetch_text(url, endpoint="ncbi_sra_efetch")
    if not text.strip():
        return {}
    try:
        root = ET.fromstring(text)
    except ET.ParseError:
        return {}
    package = root.find(".//EXPERIMENT_PACKAGE")
    if package is None:
        return {}
    experiment = package.find(".//EXPERIMENT")
    library = package.find(".//LIBRARY_DESCRIPTOR")
    sample_descriptor = package.find(".//SAMPLE_DESCRIPTOR")
    title = text_of(experiment.find("TITLE") if experiment is not None else None)
    study_ref = package.find(".//STUDY_REF")
    run_node = package.find(".//RUN")
    return {
        "experiment": attr(experiment, "accession"),
        "run": attr(run_node, "accession") or run,
        "study": attr(study_ref, "accession"),
        "sample": attr(sample_descriptor, "accession"),
        "biosample": external_id(sample_descriptor, "BioSample"),
        "gsm": external_id(sample_descriptor, "GEO") or text_of(library.find("LIBRARY_NAME") if library is not None else None),
        "title": title,
        "library_strategy": text_of(library.find("LIBRARY_STRATEGY") if library is not None else None),
        "library_source": text_of(library.find("LIBRARY_SOURCE") if library is not None else None),
        "library_selection": text_of(library.find("LIBRARY_SELECTION") if library is not None else None),
        "library_layout": library_layout(library),
        "scientific_name": text_of(package.find(".//SCIENTIFIC_NAME")),
    }


def text_of(node: ET.Element | None) -> str:
    if node is None or node.text is None:
        return ""
    return " ".join(node.text.split())


def attr(node: ET.Element | None, name: str) -> str:
    if node is None:
        return ""
    return str(node.attrib.get(name, ""))


def external_id(node: ET.Element | None, namespace: str) -> str:
    if node is None:
        return ""
    for child in node.findall(".//EXTERNAL_ID"):
        if child.attrib.get("namespace") == namespace and child.text:
            return child.text.strip()
    return ""


def library_layout(library: ET.Element | None) -> str:
    if library is None:
        return ""
    layout = library.find("LIBRARY_LAYOUT")
    if layout is None:
        return ""
    if list(layout):
        return list(layout)[0].tag
    return ""


def ena_runinfo(run: str) -> list[dict[str, str]]:
    fields = [
        "run_accession",
        "experiment_accession",
        "study_accession",
        "sample_accession",
        "secondary_sample_accession",
        "scientific_name",
        "instrument_platform",
        "library_strategy",
        "library_source",
        "library_selection",
        "library_layout",
        "read_count",
        "base_count",
        "fastq_bytes",
        "submitted_ftp",
        "experiment_title",
        "sample_title",
    ]
    params = urllib.parse.urlencode(
        {
            "accession": run,
            "result": "read_run",
            "fields": ",".join(fields),
            "format": "tsv",
        }
    )
    url = "https://www.ebi.ac.uk/ena/portal/api/filereport?" + params
    text = fetch_text(url, endpoint="ena_filereport_read_run")
    if not text.strip():
        return []
    return list(csv.DictReader(io.StringIO(text), delimiter="\t"))


def choose_ncbi_runs(group: Group, rows: list[dict[str, str]], max_runs: int) -> list[dict[str, str]]:
    priority = extract_accessions(";".join(group.raw_accessions), "srr")
    if priority and rows:
        by_run = {str(row.get("Run", "")).upper(): row for row in rows}
        selected = [by_run[run] for run in priority if run in by_run]
        if selected:
            return selected[:max_runs]
    if priority:
        return [{"Run": run} for run in priority[:max_runs]]
    return rows[:max_runs]


def refresh_group(group: Group, max_runs_per_group: int, request_delay: float) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if group.backend == "ena_run":
        for item in ena_runinfo(group.query_accession):
            rows.append(classify_ena_row(group, item))
        time.sleep(request_delay)
        return rows

    runinfo_rows = ncbi_runinfo(group.query_accession)
    time.sleep(request_delay)
    selected_rows = choose_ncbi_runs(group, runinfo_rows, max_runs_per_group)
    runinfo_by_run = {str(row.get("Run", "")).upper(): row for row in runinfo_rows}
    for runinfo_row in selected_rows:
        run = str(runinfo_row.get("Run", "")).upper()
        xml_meta = ncbi_sra_xml(run) if run else {}
        merged = {**runinfo_by_run.get(run, {}), **runinfo_row, **xml_meta}
        rows.append(classify_ncbi_row(group, merged))
        time.sleep(request_delay)
    return rows


def classify_ncbi_row(group: Group, row: dict[str, str]) -> dict[str, Any]:
    title = row.get("title", "") or row.get("LibraryName", "") or row.get("SampleName", "")
    metadata_text = " ".join(
        str(value)
        for value in [
            title,
            row.get("LibraryStrategy", ""),
            row.get("library_strategy", ""),
            row.get("LibrarySelection", ""),
            row.get("library_selection", ""),
            row.get("LibrarySource", ""),
            row.get("library_source", ""),
            row.get("ScientificName", ""),
            row.get("scientific_name", ""),
        ]
        if value
    )
    age_weeks, raw_age_token = parse_age_weeks(title)
    tissue = detect_tissue(title)
    assay_pass, assay_reason = detect_assay(metadata_text)
    species = row.get("ScientificName", "") or row.get("scientific_name", "")
    return build_classified_row(
        group=group,
        metadata_source="ncbi_runinfo_sra_xml",
        run=row.get("Run", "") or row.get("run", ""),
        experiment=row.get("Experiment", "") or row.get("experiment", ""),
        study=row.get("SRAStudy", "") or row.get("study", ""),
        bioproject=row.get("BioProject", ""),
        sample=row.get("Sample", "") or row.get("sample", ""),
        biosample=row.get("BioSample", "") or row.get("biosample", ""),
        gsm=row.get("LibraryName", "") or row.get("gsm", ""),
        library_strategy=row.get("LibraryStrategy", "") or row.get("library_strategy", ""),
        library_source=row.get("LibrarySource", "") or row.get("library_source", ""),
        library_selection=row.get("LibrarySelection", "") or row.get("library_selection", ""),
        library_layout=row.get("LibraryLayout", "") or row.get("library_layout", ""),
        size_mb=safe_float(row.get("size_MB", "")),
        scientific_name=species,
        title_text=title,
        age_weeks=age_weeks,
        raw_age_token=raw_age_token,
        tissue=tissue,
        assay_pass=assay_pass,
        assay_reason=assay_reason,
    )


def classify_ena_row(group: Group, row: dict[str, str]) -> dict[str, Any]:
    title = " ".join(value for value in [row.get("experiment_title", ""), row.get("sample_title", "")] if value)
    metadata_text = " ".join(
        str(value)
        for value in [
            title,
            row.get("library_strategy", ""),
            row.get("library_selection", ""),
            row.get("library_source", ""),
            row.get("scientific_name", ""),
        ]
        if value
    )
    age_weeks, raw_age_token = parse_age_weeks(title)
    tissue = detect_tissue(title)
    assay_pass, assay_reason = detect_assay(metadata_text)
    size_mb = size_mb_from_ena(row.get("fastq_bytes", ""))
    return build_classified_row(
        group=group,
        metadata_source="ena_read_run_filereport",
        run=row.get("run_accession", ""),
        experiment=row.get("experiment_accession", ""),
        study=row.get("study_accession", ""),
        bioproject="",
        sample=row.get("sample_accession", ""),
        biosample=row.get("secondary_sample_accession", ""),
        gsm="",
        library_strategy=row.get("library_strategy", ""),
        library_source=row.get("library_source", ""),
        library_selection=row.get("library_selection", ""),
        library_layout=row.get("library_layout", ""),
        size_mb=size_mb,
        scientific_name=row.get("scientific_name", ""),
        title_text=title,
        age_weeks=age_weeks,
        raw_age_token=raw_age_token,
        tissue=tissue,
        assay_pass=assay_pass,
        assay_reason=assay_reason,
    )


def build_classified_row(
    *,
    group: Group,
    metadata_source: str,
    run: str,
    experiment: str,
    study: str,
    bioproject: str,
    sample: str,
    biosample: str,
    gsm: str,
    library_strategy: str,
    library_source: str,
    library_selection: str,
    library_layout: str,
    size_mb: float | None,
    scientific_name: str,
    title_text: str,
    age_weeks: float | None,
    raw_age_token: str,
    tissue: str,
    assay_pass: bool,
    assay_reason: str,
) -> dict[str, Any]:
    species_pass = "mus musculus" in scientific_name.lower()
    target_tissue_pass = tissue in TARGET_TISSUES
    exact_age_pass = age_weeks is not None
    age_group = detect_age_group(title_text)
    old_age_pass = bool(age_weeks is not None and age_weeks >= OLD_WEEKS) or age_group == "old"
    blockers = detect_context_blockers(title_text + " " + library_strategy + " " + library_source)
    bulk_context_pass = not blockers
    sample_decision = decide_sample(
        species_pass=species_pass,
        assay_pass=assay_pass,
        target_tissue_pass=target_tissue_pass,
        exact_age_pass=exact_age_pass,
        old_age_pass=old_age_pass,
        bulk_context_pass=bulk_context_pass,
        age_group=age_group,
    )
    return {
        "group_id": group.group_id,
        "backend": group.backend,
        "query_accession": group.query_accession,
        "source_project_keys": ";".join(group.source_project_keys),
        "metadata_source": metadata_source,
        "run": run,
        "experiment": experiment,
        "study": study,
        "bioproject": bioproject,
        "sample": sample,
        "biosample": biosample,
        "gsm": gsm,
        "scientific_name": scientific_name,
        "library_strategy": library_strategy,
        "library_source": library_source,
        "library_selection": library_selection,
        "library_layout": library_layout,
        "size_mb": round(size_mb, 3) if size_mb is not None else "",
        "title_text": title_text,
        "tissue": tissue,
        "age_weeks": round(age_weeks, 3) if age_weeks is not None else "",
        "raw_age_token": raw_age_token,
        "age_group": age_group,
        "species_pass": species_pass,
        "assay_pass": assay_pass,
        "assay_reason": assay_reason,
        "target_tissue_pass": target_tissue_pass,
        "exact_age_pass": exact_age_pass,
        "old_age_pass": old_age_pass,
        "bulk_context_pass": bulk_context_pass,
        "context_blockers": ";".join(blockers),
        "sample_decision": sample_decision,
        "fastq_download_authorized": False,
        "bismark_authorized": False,
    }


def parse_age_weeks(text: str) -> tuple[float | None, str]:
    match = AGE_RE.search(text or "")
    if not match:
        return None, ""
    value = float(match.group("value"))
    unit = match.group("unit").lower()
    if unit in {"d", "day", "days"}:
        weeks = value / 7.0
    elif unit in {"w", "wk", "wks", "week", "weeks"}:
        weeks = value
    elif unit in {"m", "mo", "mos", "month", "months"}:
        weeks = value * 30.42 / 7.0
    else:
        weeks = value * 365.0 / 7.0
    return weeks, match.group(0)


def detect_age_group(text: str) -> str:
    lower = (text or "").lower()
    if re.search(r"\b(aged|old)\b", lower):
        return "old"
    if re.search(r"\b(young|juvenile)\b", lower):
        return "young"
    return ""


def detect_tissue(text: str) -> str:
    lower = (text or "").lower().replace("_", " ")
    if "brain cortex" in lower or re.search(r"\bcortex\b", lower):
        return "brain_cortex"
    if re.search(r"\bheart\b|\bcardiac\b", lower):
        return "heart"
    if re.search(r"\blung\b", lower):
        return "lung"
    if "small intestinal" in lower or "intestine" in lower or "intestinal" in lower or "crypt" in lower:
        return "intestine"
    if "liver" in lower:
        return "liver"
    if "kidney" in lower:
        return "kidney"
    if "spleen" in lower:
        return "spleen"
    return ""


def detect_assay(text: str) -> tuple[bool, str]:
    lower = (text or "").lower()
    for blocker in ASSAY_BLOCKERS:
        if blocker in lower:
            return False, f"blocked_{blocker.replace(' ', '_')}"
    if any(term in lower for term in ASSAY_TERMS):
        return True, "bisulfite_or_methylation_term_found"
    return False, "no_bisulfite_methylation_term"


def detect_context_blockers(text: str) -> list[str]:
    lower = (text or "").lower()
    return sorted(term for term in CONTEXT_BLOCKERS if term in lower)


def decide_sample(
    *,
    species_pass: bool,
    assay_pass: bool,
    target_tissue_pass: bool,
    exact_age_pass: bool,
    old_age_pass: bool,
    bulk_context_pass: bool,
    age_group: str,
) -> str:
    if not species_pass:
        return "blocked_non_mouse"
    if not assay_pass:
        return "blocked_non_bisulfite_methylation_assay"
    if not target_tissue_pass:
        return "blocked_non_target_tissue"
    if not bulk_context_pass:
        return "blocked_context_not_bulk_headline"
    if exact_age_pass and old_age_pass:
        return "eligible_old_exact_target_pending_approval"
    if age_group == "old":
        return "diagnostic_only_age_group_no_exact_age"
    if exact_age_pass:
        return "diagnostic_exact_age_but_not_old_target"
    return "blocked_missing_sample_specific_age"


def safe_float(value: Any) -> float | None:
    try:
        text = str(value).strip()
        if not text:
            return None
        return float(text)
    except (TypeError, ValueError):
        return None


def size_mb_from_ena(value: str) -> float | None:
    parts = [safe_float(part) for part in str(value or "").split(";") if str(part).strip()]
    numeric = [part for part in parts if part is not None]
    if not numeric:
        return None
    return sum(numeric) / 1024.0 / 1024.0


def summarize(rows: list[dict[str, Any]], groups: list[Group]) -> list[dict[str, Any]]:
    by_group: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        by_group[str(row["group_id"])].append(row)
    group_by_id = {group.group_id: group for group in groups}
    summary_rows: list[dict[str, Any]] = []
    for group_id in sorted(group_by_id):
        group = group_by_id[group_id]
        group_rows = by_group.get(group_id, [])
        decisions = sorted({str(row.get("sample_decision", "")) for row in group_rows if row.get("sample_decision")})
        eligible = [row for row in group_rows if row.get("sample_decision") == "eligible_old_exact_target_pending_approval"]
        diagnostic_age = [row for row in group_rows if row.get("sample_decision") == "diagnostic_only_age_group_no_exact_age"]
        decision, action, reason = decide_group(group_rows, eligible, diagnostic_age)
        tissues = sorted({str(row.get("tissue", "")) for row in group_rows if row.get("tissue")})
        age_values = [safe_float(row.get("age_weeks")) for row in group_rows]
        ages = [age for age in age_values if age is not None]
        summary_rows.append(
            {
                "group_id": group_id,
                "backend": group.backend,
                "query_accession": group.query_accession,
                "source_project_keys": ";".join(group.source_project_keys),
                "backlog_row_count": group.row_count,
                "metadata_rows": len(group_rows),
                "route_b_decision": decision,
                "decision_reason": reason,
                "recommended_action": action,
                "n_species_pass": count_true(group_rows, "species_pass"),
                "n_assay_pass": count_true(group_rows, "assay_pass"),
                "n_target_tissue": count_true(group_rows, "target_tissue_pass"),
                "n_exact_age": count_true(group_rows, "exact_age_pass"),
                "n_old_age": count_true(group_rows, "old_age_pass"),
                "n_bulk_context_pass": count_true(group_rows, "bulk_context_pass"),
                "n_headline_eligible": len(eligible),
                "n_diagnostic_age_group": len(diagnostic_age),
                "tissues": ";".join(tissues),
                "min_age_weeks": round(min(ages), 3) if ages else "",
                "max_age_weeks": round(max(ages), 3) if ages else "",
                "sample_decisions": ";".join(decisions),
                "raw_accession_examples": ";".join(group.raw_accessions[:12]),
            }
        )
    return summary_rows


def count_true(rows: list[dict[str, Any]], field: str) -> int:
    return sum(bool(row.get(field)) for row in rows)


def decide_group(
    rows: list[dict[str, Any]],
    eligible: list[dict[str, Any]],
    diagnostic_age: list[dict[str, Any]],
) -> tuple[str, str, str]:
    if len(eligible) >= 2:
        return (
            "route_b_pilot_eligible_pending_explicit_approval",
            "prepare_separate_minimal_etl_approval_before_any_fastq_download",
            "At least two old exact-age target-tissue bulk-compatible methylation runs were found.",
        )
    if eligible:
        return (
            "insufficient_replicates_for_route_b_pilot",
            "do_not_download_fastq; seek_more_replicates_or_route_a",
            "Only one old exact-age target-tissue bulk-compatible methylation run was found.",
        )
    if diagnostic_age:
        return (
            "diagnostic_only_age_group_no_exact_age",
            "do_not_download_fastq_for_headline; diagnostic_only_if_separately_approved",
            "Target-tissue methylation data have old/aged labels but no sample-specific exact age.",
        )
    if not rows:
        return (
            "blocked_no_official_metadata_rows",
            "do_not_download_fastq",
            "Official metadata refresh returned no parseable run rows.",
        )
    if any(row.get("sample_decision") == "blocked_context_not_bulk_headline" for row in rows):
        return (
            "blocked_context_not_bulk_headline",
            "do_not_download_fastq",
            "Candidate is cell-type/single-cell/other non-bulk context.",
        )
    if not any(row.get("assay_pass") for row in rows):
        return (
            "blocked_non_bisulfite_methylation_assay",
            "do_not_download_fastq",
            "Official metadata does not indicate bisulfite/methylation sequencing.",
        )
    if not any(row.get("target_tissue_pass") for row in rows):
        return (
            "blocked_non_target_tissue",
            "do_not_download_fastq",
            "Official metadata does not include brain_cortex/heart/lung target tissue.",
        )
    return (
        "not_route_b_headline_candidate",
        "do_not_download_fastq",
        "Candidate failed at least one Route B headline gate.",
    )


def write_report(summary_rows: list[dict[str, Any]], row_count: int, state: dict[str, Any]) -> None:
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    headline = [row for row in summary_rows if row.get("route_b_decision") == "route_b_pilot_eligible_pending_explicit_approval"]
    diagnostic = [
        row
        for row in summary_rows
        if row.get("route_b_decision") in {"diagnostic_only_age_group_no_exact_age", "insufficient_replicates_for_route_b_pilot"}
    ]
    blocked = [row for row in summary_rows if row not in headline and row not in diagnostic]
    def report_rank(row: dict[str, Any]) -> tuple[int, str]:
        decision = str(row.get("route_b_decision", ""))
        if decision == "route_b_pilot_eligible_pending_explicit_approval":
            return (0, str(row.get("group_id", "")))
        if str(row.get("backend")) == "ncbi_group":
            return (1, str(row.get("group_id", "")))
        if decision != "blocked_non_bisulfite_methylation_assay":
            return (2, str(row.get("group_id", "")))
        return (3, str(row.get("group_id", "")))

    top_rows = sorted(summary_rows, key=report_rank)
    lines = [
        "# Route B Manual Backlog Official Metadata Refresh Report",
        "",
        f"Date: {datetime.now().isoformat(timespec='seconds')}",
        "",
        "## Summary",
        "",
        "This refresh checks the manual Route B raw-lead backlog using official NCBI RunInfo/SRA XML and ENA read_run metadata only.",
        "",
        "No FASTQ download, Bismark, training, matrix rebuild, or autoresearch was run.",
        "",
        "## Decision",
        "",
        f"- Metadata groups refreshed: `{len(summary_rows)}`.",
        f"- Run rows audited: `{row_count}`.",
        f"- Headline Route B pilot-ready groups: `{len(headline)}`.",
        f"- Diagnostic-only or insufficient-replicate groups: `{len(diagnostic)}`.",
        f"- Blocked groups: `{len(blocked)}`.",
        f"- Overall decision: `{state['decision']}`.",
        "",
        "No raw FASTQ pilot is authorized by this report.",
        "",
        "## Key Findings",
        "",
    ]
    for row in top_rows[:12]:
        lines.extend(
            [
                f"### {row.get('group_id')}",
                "",
                f"- decision: `{row.get('route_b_decision')}`",
                f"- reason: {row.get('decision_reason')}",
                f"- runs audited: `{row.get('metadata_rows')}`, tissues: `{row.get('tissues') or 'N/A'}`",
                f"- age range weeks: `{row.get('min_age_weeks') or 'N/A'}` to `{row.get('max_age_weeks') or 'N/A'}`",
                f"- assay-pass / target / exact-age / old-age: `{row.get('n_assay_pass')}` / `{row.get('n_target_tissue')}` / `{row.get('n_exact_age')}` / `{row.get('n_old_age')}`",
                f"- recommended action: `{row.get('recommended_action')}`",
                "",
            ]
        )
    lines.extend(
        [
            "## Outputs",
            "",
            "- `results/route_b_candidate_audit/route_b_manual_backlog_runinfo_refresh.csv`",
            "- `results/route_b_candidate_audit/route_b_manual_backlog_candidate_summary.csv`",
            "- `results/route_b_candidate_audit/route_b_manual_backlog_network_log.jsonl`",
            "- `results/route_b_candidate_audit/route_b_manual_backlog_decision_state.json`",
            "",
            "## Guardrails",
            "",
            "- This is metadata refresh only.",
            "- Route B pilot still requires a separate explicit minimal ETL approval.",
            "- Age-group-only candidates cannot become headline chronological benchmarks.",
            "- Cell-type-specific, RNA-Seq, or non-target tissue candidates remain blocked.",
        ]
    )
    REPORT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def update_index() -> None:
    if not INDEX_PATH.exists():
        return
    text = INDEX_PATH.read_text(encoding="utf-8")
    bullet = "- Route B manual backlog refresh report: `doc/20_analysis/35_20260519_route_b_manual_backlog_runinfo_refresh_report.md`"
    if bullet in text:
        return
    marker = "- Route B candidate audit report: `doc/20_analysis/34_20260519_route_b_candidate_audit_report.md`"
    if marker in text:
        text = text.replace(marker, marker + "\n" + bullet)
    else:
        text += "\n" + bullet + "\n"
    INDEX_PATH.write_text(text, encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-runs-per-group", type=int, default=40)
    parser.add_argument("--request-delay", type=float, default=0.34)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if NETWORK_LOG.exists():
        NETWORK_LOG.unlink()
    backlog_rows = read_csv(BACKLOG)
    groups = make_groups(backlog_rows)
    all_rows: list[dict[str, Any]] = []
    for group in groups:
        all_rows.extend(refresh_group(group, args.max_runs_per_group, args.request_delay))
    summary_rows = summarize(all_rows, groups)
    write_csv(ROW_OUT, all_rows)
    write_csv(SUMMARY_OUT, summary_rows)
    headline_count = sum(
        row.get("route_b_decision") == "route_b_pilot_eligible_pending_explicit_approval" for row in summary_rows
    )
    diagnostic_count = sum(
        row.get("route_b_decision") in {"diagnostic_only_age_group_no_exact_age", "insufficient_replicates_for_route_b_pilot"}
        for row in summary_rows
    )
    state = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "loop_version": "route_b_manual_backlog_refresh_v1",
        "input_backlog": str(BACKLOG),
        "groups_refreshed": len(groups),
        "metadata_rows": len(all_rows),
        "headline_route_b_candidate_count": headline_count,
        "diagnostic_or_insufficient_candidate_count": diagnostic_count,
        "decision": "route_b_headline_pilot_ready_pending_explicit_approval" if headline_count else "no_route_b_headline_pilot_candidate_ready",
        "download_authorized": False,
        "fastq_download_authorized": False,
        "bismark_authorized": False,
        "training_authorized": False,
        "autoresearch_authorized": False,
        "next_action": (
            "prepare_separate_minimal_etl_approval_for_selected_headline_candidate"
            if headline_count
            else "do_not_download_fastq; keep_route_a_or_benchmark_route_c_as_current_baseline"
        ),
    }
    write_json(STATE_OUT, state)
    write_report(summary_rows, len(all_rows), state)
    update_index()
    print(json.dumps(state, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

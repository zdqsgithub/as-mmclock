#!/usr/bin/env python3
"""Discover old-age mouse methylation GEO candidates without large downloads.

This v7.9 preflight uses official NCBI/GEO endpoints only:

- E-Utils for searchable GEO metadata.
- GEO family SOFT files for sample/series metadata.
- GEO FTP ``suppl/filelist.txt`` plus HEAD requests for supplement facts.

It intentionally does not download large ``*_RAW.tar`` archives. The goal is to
rank candidate datasets that can fix the v7.8 age-range blocker.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

ROOT = Path("/home/zdq-as/mouse_methyl_work")
META_DIR = ROOT / "metadata"
SOFT_CACHE_DIR = META_DIR / "geo_candidate_soft"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

OFFICIAL_SOURCE_URLS = [
    "https://www.ncbi.nlm.nih.gov/geo/info/download.html",
    "https://www.ncbi.nlm.nih.gov/geo/info/geo_paccess.html",
    "https://www.ncbi.nlm.nih.gov/geo/info/soft.html",
]

DEFAULT_SEARCH_TERMS = [
    'Mus musculus[ORGN] AND (RRBS OR "reduced representation bisulfite") AND (aging OR age OR old OR month) AND gse[ETYP]',
    'Mus musculus[ORGN] AND "Methylation profiling by high throughput sequencing" AND aging AND gse[ETYP]',
    'Mus musculus[ORGN] AND "DNA methylation" AND "24 month" AND gse[ETYP]',
    'Mus musculus[ORGN] AND "DNA methylation" AND "26 month" AND gse[ETYP]',
]

SEED_ACCESSIONS = [
    "GSE213628",  # RRBS, multi-tissue, 4/12/18/24 month and old heart samples.
    "GSE224442",  # RRBS, parabiosis/recovery, liver and blood, young/old.
    "GSE92486",  # WGBS, liver, AL/DR, 5/26 month.
    "GSE233879",  # RRBS subset, old HSC/blood aging/rejuvenation.
    "GSE295059",  # RRBS, 4/24 month intestinal organoids.
    "GSE286302",  # RRBS, multi-tissue circadian intervention in aging mice.
    "GSE175410",  # RRBS, late-life exercise skeletal muscle.
    "GSE276335",  # WGBS/RRBS-like COV, HSC aging perturbation.
]

MONTH_TO_WEEK = 30.42 / 7.0
EXACT_AGE_RE = re.compile(
    r"(?<![A-Za-z0-9])(?P<value>\d+(?:\.\d+)?)\s*"
    r"(?P<unit>d|day|days|w|wk|wks|week|weeks|m/o|mo|mos|month|months|m|y|yr|yrs|year|years)"
    r"(?![A-Za-z0-9])",
    re.IGNORECASE,
)
AGE_LABEL_UNIT_RE = re.compile(
    r"age\s*\((?P<unit>days?|weeks?|months?|years?|d|w|m|y|yrs?|mos?)\)\s*:\s*(?P<value>\d+(?:\.\d+)?)",
    re.IGNORECASE,
)
AGE_LABEL_VALUE_RE = re.compile(
    r"age\s*:\s*(?P<value>\d+(?:\.\d+)?)\s*(?P<unit>d|day|days|w|wk|wks|week|weeks|m/o|mo|mos|month|months|m|y|yr|yrs|year|years)",
    re.IGNORECASE,
)

TISSUE_PATTERNS = {
    "blood": ["whole blood", "blood", "pbmc"],
    "liver": ["liver", "hepatic"],
    "heart": ["heart", "cardiac"],
    "lung": ["lung"],
    "kidney": ["kidney"],
    "spleen": ["spleen"],
    "skeletal_muscle": ["skeletal muscle", "quadriceps", "muscle", "myonuclei"],
    "brain": ["brain", "cortex", "hippocampus", "pvn", "hypothalamus"],
    "intestine": ["intestine", "colon", "crypt", "organoid"],
    "hsc_bone_marrow": ["hsc", "hematopoietic", "bone marrow", "cd150"],
    "adipose": ["adipose"],
}

INTERVENTION_PATTERNS = {
    "dietary_restriction": ["dietary restriction", "calorie restriction", "caloric restriction", "diet: dr", "group: dr"],
    "parabiosis": ["parabiosis", "heterochronic", "isochronic"],
    "recovery": ["recovery", "rec"],
    "exercise": ["exercise", "running", "wheel", "treadmill"],
    "circadian": ["circadian", "3da", "chemogenetic", "pvn"],
    "reprogramming": ["reprogramming", "oskm", "yamanaka"],
    "drug": ["decitabine", "elamipretide", "glp-1", "rapamycin"],
    "genetic": ["knock-out", "knockout", "ko", "gadd45", "sirt6"],
}


@dataclass(frozen=True)
class UrlFetch:
    text: str
    error: str = ""


def series_bucket(gse: str) -> str:
    match = re.fullmatch(r"GSE(\d+)", gse)
    if not match:
        raise ValueError(f"Invalid GSE accession: {gse}")
    digits = match.group(1)
    return f"GSE{digits[:-3]}nnn"


def geo_url(gse: str) -> str:
    return f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}"


def soft_url(gse: str) -> str:
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{series_bucket(gse)}/{gse}/soft/{gse}_family.soft.gz"


def suppl_dir_url(gse: str) -> str:
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{series_bucket(gse)}/{gse}/suppl/"


def fetch_text(url: str, timeout: int = 45) -> UrlFetch:
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "mouse-methyl-v7.9-preflight/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as response:
            data = response.read()
            if url.endswith(".gz"):
                data = gzip.decompress(data)
            return UrlFetch(data.decode("utf-8", errors="replace"))
    except Exception as exc:
        return UrlFetch("", str(exc)[:500])


def head_url(url: str, timeout: int = 30) -> dict:
    request = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "mouse-methyl-v7.9-preflight/1.0"})
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


def eutils_get(endpoint: str, params: dict, timeout: int = 45) -> str:
    url = EUTILS + endpoint + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"User-Agent": "mouse-methyl-v7.9-preflight/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as response:
        return response.read().decode("utf-8", errors="replace")


def eutils_search(term: str, retmax: int) -> list[str]:
    xml = eutils_get("esearch.fcgi", {"db": "gds", "term": term, "retmax": retmax})
    root = ET.fromstring(xml)
    return [node.text for node in root.findall(".//Id") if node.text]


def parse_esummary(xml: str) -> list[dict]:
    root = ET.fromstring(xml)
    rows = []
    for doc in root.findall(".//DocumentSummary"):
        samples = doc.findall("Samples/Sample")
        rows.append(
            {
                "uid": doc.attrib.get("uid", ""),
                "dataset": doc.findtext("Accession") or "",
                "title": doc.findtext("title") or "",
                "summary": doc.findtext("summary") or "",
                "taxon": doc.findtext("taxon") or "",
                "entry_type": doc.findtext("entryType") or "",
                "experiment_type": doc.findtext("gdsType") or "",
                "public_date": doc.findtext("PDAT") or "",
                "supp_file": doc.findtext("suppFile") or "",
                "n_samples_esummary": len(samples),
                "sample_titles_esummary": " | ".join((sample.findtext("Title") or "") for sample in samples[:12]),
            }
        )
    return rows


def eutils_summaries_for_uids(uids: list[str]) -> list[dict]:
    rows = []
    for start in range(0, len(uids), 50):
        chunk = uids[start : start + 50]
        if not chunk:
            continue
        xml = eutils_get("esummary.fcgi", {"db": "gds", "version": "2.0", "id": ",".join(chunk)})
        rows.extend(parse_esummary(xml))
        time.sleep(0.34)
    return rows


def eutils_summary_for_accession(accession: str) -> dict:
    ids = eutils_search(f"{accession}[ACCN]", retmax=5)
    if not ids:
        return {"dataset": accession, "eutils_error": "no_uid_for_accession"}
    rows = eutils_summaries_for_uids(ids)
    for row in rows:
        if row.get("dataset") == accession:
            return row
    return rows[0] if rows else {"dataset": accession, "eutils_error": "summary_not_found"}


def discover_accessions(search_terms: list[str], retmax: int) -> tuple[set[str], dict[str, list[str]], dict[str, dict]]:
    accessions: set[str] = set()
    source_terms: dict[str, list[str]] = defaultdict(list)
    summaries: dict[str, dict] = {}
    for term in search_terms:
        try:
            uids = eutils_search(term, retmax=retmax)
            rows = eutils_summaries_for_uids(uids)
        except Exception as exc:
            print(f"[WARN] E-Utils search failed for {term!r}: {exc}")
            continue
        for row in rows:
            acc = row.get("dataset")
            if not acc or not acc.startswith("GSE"):
                continue
            accessions.add(acc)
            source_terms[acc].append(term)
            summaries[acc] = row
    return accessions, source_terms, summaries


def cached_soft_text(gse: str, refresh: bool = False) -> UrlFetch:
    SOFT_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache = SOFT_CACHE_DIR / f"{gse}_family.soft.gz"
    if cache.exists() and not refresh:
        try:
            return UrlFetch(gzip.decompress(cache.read_bytes()).decode("utf-8", errors="replace"))
        except Exception as exc:
            return UrlFetch("", f"cache_read_failed: {exc}")
    url = soft_url(gse)
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "mouse-methyl-v7.9-preflight/1.0"})
        with urllib.request.urlopen(req, timeout=60) as response:
            data = response.read()
        cache.write_bytes(data)
        return UrlFetch(gzip.decompress(data).decode("utf-8", errors="replace"))
    except Exception as exc:
        return UrlFetch("", str(exc)[:500])


def parse_soft(text: str) -> tuple[dict[str, str], list[dict[str, str]]]:
    series: dict[str, list[str]] = defaultdict(list)
    samples: list[dict[str, list[str]]] = []
    current_type = ""
    current_sample: dict[str, list[str]] | None = None
    for raw_line in text.splitlines():
        if raw_line.startswith("^SERIES"):
            current_type = "series"
            current_sample = None
            if "=" in raw_line:
                series["series_accession"].append(raw_line.split("=", 1)[1].strip())
            continue
        if raw_line.startswith("^SAMPLE"):
            current_type = "sample"
            current_sample = defaultdict(list)
            samples.append(current_sample)
            if "=" in raw_line:
                current_sample["sample_id"].append(raw_line.split("=", 1)[1].strip())
            continue
        if not raw_line.startswith("!") or "=" not in raw_line:
            continue
        label, value = raw_line.split("=", 1)
        label = label.lstrip("!").strip()
        value = value.strip()
        if current_type == "series":
            series[label].append(value)
        elif current_type == "sample" and current_sample is not None:
            current_sample[label].append(value)

    flat_series = {key: " | ".join(values) for key, values in series.items()}
    flat_samples: list[dict[str, str]] = []
    for sample in samples:
        flat_samples.append({key: " | ".join(values) for key, values in sample.items()})
    return flat_series, flat_samples


SAMPLE_SPECIFIC_FIELDS = [
    "Sample_title",
    "Sample_source_name_ch1",
    "Sample_characteristics_ch1",
    "Sample_description",
]


def text_blob(sample: dict[str, str], fields: Iterable[str] | None = None) -> str:
    if fields is None:
        fields = SAMPLE_SPECIFIC_FIELDS
    return " ".join(str(sample.get(field, "")) for field in fields)


def age_to_weeks(value: float, unit: str) -> float | None:
    unit = unit.lower()
    if unit in {"d", "day", "days"}:
        return value / 7.0
    if unit in {"w", "wk", "wks", "week", "weeks"}:
        return value
    if unit in {"m/o", "mo", "mos", "month", "months", "m"}:
        return value * MONTH_TO_WEEK
    if unit in {"y", "yr", "yrs", "year", "years"}:
        return value * 365.0 / 7.0
    return None


def parse_age_weeks(text: str) -> tuple[float | None, str, str]:
    candidates = []
    for pattern in (AGE_LABEL_UNIT_RE, AGE_LABEL_VALUE_RE):
        for match in pattern.finditer(text):
            weeks = age_to_weeks(float(match.group("value")), match.group("unit"))
            if weeks is not None and 0 <= weeks <= 220:
                candidates.append((weeks, match.group(0), match.group("unit").lower(), 0))
    for match in EXACT_AGE_RE.finditer(text):
        value = float(match.group("value"))
        unit = match.group("unit").lower()
        weeks = age_to_weeks(value, unit)
        if weeks is None:
            continue
        if 0 <= weeks <= 220:
            candidates.append((weeks, match.group(0), unit, 1))
    if not candidates:
        return None, "", ""
    candidates.sort(key=lambda item: (item[3], -item[0]))
    weeks, raw, unit, _ = candidates[0]
    return round(float(weeks), 3), raw, unit


def detect_terms(text: str, patterns: dict[str, list[str]]) -> set[str]:
    lowered = f" {text.lower().replace('-', ' ')} "
    found = set()
    for label, needles in patterns.items():
        for needle in needles:
            normalized = f" {needle.lower().replace('-', ' ')} "
            if normalized.strip() in lowered:
                found.add(label)
                break
    return found


def parse_filelist(text: str) -> list[dict]:
    if not text.strip():
        return []
    reader = csv.DictReader(text.splitlines(), delimiter="\t")
    rows = []
    for row in reader:
        name = row.get("Name")
        if name:
            rows.append(row)
    return rows


def fetch_supplements(gse: str, do_head: bool = True) -> list[dict]:
    base = suppl_dir_url(gse)
    filelist = fetch_text(base + "filelist.txt")
    rows = parse_filelist(filelist.text)
    supplements = []
    for row in rows:
        name = row.get("Name", "")
        if not name:
            continue
        url = base + urllib.parse.quote(name)
        head = head_url(url) if do_head and row.get("#Archive/File") == "Archive" else {}
        supplements.append(
            {
                "dataset": gse,
                "supplement_name": name,
                "archive_or_file": row.get("#Archive/File", ""),
                "filelist_time": row.get("Time", ""),
                "filelist_size_bytes": int(row.get("Size") or 0) or None,
                "filelist_type": row.get("Type", ""),
                "supplement_url": url,
                **head,
            }
        )
    if not supplements and filelist.error:
        supplements.append(
            {
                "dataset": gse,
                "supplement_name": "",
                "archive_or_file": "",
                "filelist_time": "",
                "filelist_size_bytes": None,
                "filelist_type": "",
                "supplement_url": base + "filelist.txt",
                "filelist_error": filelist.error,
            }
        )
    return supplements


def infer_schema(summary: dict, supplements: list[dict]) -> tuple[str, str]:
    names = " ".join(row.get("supplement_name", "") for row in supplements).lower()
    types = {str(row.get("filelist_type", "")).upper() for row in supplements}
    title_summary = f"{summary.get('title', '')} {summary.get('summary', '')}".lower()
    if "cov" in types or ".cov" in names or "bismark.cov" in names:
        if "wgbs" in title_summary or "whole genome" in title_summary:
            return "wgbs_cov_per_sample_tar", ""
        return "bismark_cov_per_sample_tar", ""
    if "bedgraph" in types or "bedgraph" in names:
        return "bedgraph_low_priority_adapter", "bedgraph_schema_not_currently_supported"
    if "raw.tar" in names and "methylation" in title_summary:
        return "unknown_methylation_raw_tar", "requires_parser_smoke_before_download"
    return "unknown_or_non_methylation", "no_parseable_processed_methylation_schema_detected"


def summarize_samples(gse: str, samples: list[dict]) -> tuple[dict, list[dict]]:
    sample_rows = []
    ages = []
    tissue_counter: Counter[str] = Counter()
    intervention_counter: Counter[str] = Counter()
    exact_age_count = 0
    old_label_count = 0
    for sample in samples:
        sample_blob = text_blob(sample, SAMPLE_SPECIFIC_FIELDS)
        protocol_blob = text_blob(
            sample,
            [
                "Sample_title",
                "Sample_source_name_ch1",
                "Sample_characteristics_ch1",
                "Sample_treatment_protocol_ch1",
            ],
        )
        age_weeks, raw_age, age_unit = parse_age_weeks(sample_blob)
        tissues = sorted(detect_terms(sample_blob, TISSUE_PATTERNS))
        interventions = sorted(detect_terms(protocol_blob, INTERVENTION_PATTERNS))
        if age_weeks is not None:
            ages.append(age_weeks)
            exact_age_count += 1
        if "old" in sample_blob.lower() or "aged" in sample_blob.lower():
            old_label_count += 1
        for tissue in tissues or ["unknown"]:
            tissue_counter[tissue] += 1
        for intervention in interventions or ["none_detected"]:
            intervention_counter[intervention] += 1
        sample_rows.append(
            {
                "dataset": gse,
                "sample_id": sample.get("sample_id", ""),
                "title": sample.get("Sample_title", ""),
                "age_weeks": age_weeks,
                "raw_age_token": raw_age,
                "raw_age_unit": age_unit,
                "tissue_guess": ";".join(tissues) if tissues else "",
                "intervention_guess": ";".join(interventions) if interventions else "",
                "source_name": sample.get("Sample_source_name_ch1", ""),
                "characteristics": sample.get("Sample_characteristics_ch1", "")[:500],
            }
        )
    summary = {
        "n_samples_soft": len(samples),
        "exact_age_sample_count": exact_age_count,
        "old_label_sample_count": old_label_count,
        "min_age_weeks": round(float(min(ages)), 3) if ages else None,
        "max_age_weeks": round(float(max(ages)), 3) if ages else None,
        "max_age_months": round(float(max(ages)) / MONTH_TO_WEEK, 3) if ages else None,
        "age_values_weeks": ";".join(str(round(v, 3)) for v in sorted(set(ages))) if ages else "",
        "tissue_count": len([key for key in tissue_counter if key != "unknown"]),
        "tissues": ";".join(f"{key}:{value}" for key, value in tissue_counter.most_common()),
        "intervention_terms": ";".join(f"{key}:{value}" for key, value in intervention_counter.most_common()),
    }
    return summary, sample_rows


def summarize_supplements(supplements: list[dict]) -> dict:
    archive_rows = [row for row in supplements if row.get("archive_or_file") == "Archive"]
    raw_archive = next((row for row in archive_rows if "RAW.tar" in row.get("supplement_name", "")), None)
    type_counter = Counter(str(row.get("filelist_type", "")) for row in supplements if row.get("filelist_type"))
    return {
        "supplement_count": len(supplements),
        "supplement_types": ";".join(f"{key}:{value}" for key, value in type_counter.most_common()),
        "preferred_processed_file": (raw_archive or archive_rows[0]).get("supplement_name", "") if archive_rows else "",
        "preferred_size_bytes": (raw_archive or archive_rows[0]).get("filelist_size_bytes") if archive_rows else None,
        "preferred_url": (raw_archive or archive_rows[0]).get("supplement_url", "") if archive_rows else "",
    }


def rank_candidate(row: dict) -> tuple[str, str, str]:
    blockers = []
    score = 0
    if row.get("taxon") == "Mus musculus":
        score += 2
    if "Methylation profiling" in str(row.get("experiment_type", "")):
        score += 3
    if row.get("exact_age_sample_count", 0) >= 10:
        score += 3
    else:
        blockers.append("limited_exact_age_parse")
    if (row.get("max_age_months") or 0) >= 24:
        score += 4
    else:
        blockers.append("no_24_month_plus_exact_age")
    if row.get("tissue_count", 0) >= 4:
        score += 4
    elif row.get("tissue_count", 0) >= 2:
        score += 2
    else:
        blockers.append("single_tissue_or_celltype")
    if "COV" in str(row.get("supplement_types", "")) or "cov" in str(row.get("processed_schema_guess", "")):
        score += 3
    else:
        blockers.append("processed_cov_not_detected")
    if row.get("n_samples_soft", 0) >= 50:
        score += 2
    if row.get("intervention_terms") and "none_detected" not in str(row.get("intervention_terms")):
        score += 1
    text = f"{row.get('title','')} {row.get('summary','')}".lower()
    title_text = str(row.get("title", "")).lower()
    if "organoid" in title_text or "in vitro modeling" in title_text:
        score -= 3
        blockers.append("in_vitro_or_organoid_context")
    if "hsc" in text or "hematopoietic stem" in text:
        score -= 2
        blockers.append("celltype_specific_hsc_context")
    if "superseries" in text:
        score -= 3
        blockers.append("superseries_use_subseries_when_available")

    has_cov = "cov" in str(row.get("processed_schema_guess", "")).lower()
    has_old_age = (row.get("max_age_months") or 0) >= 24
    has_intervention = row.get("intervention_terms") and "none_detected" not in str(row.get("intervention_terms"))
    is_multitissue = row.get("tissue_count", 0) >= 4
    if "superseries_use_subseries_when_available" in blockers:
        tier = "P3_auxiliary_or_intervention"
        action = "use_subseries_not_superseries_for_download"
    elif score >= 16 and is_multitissue and has_old_age and has_cov:
        tier = "P1_download_next"
        action = "add_to_v8_download_and_conversion_plan"
    elif score >= 12 or (has_old_age and has_cov and has_intervention):
        tier = "P2_adapter_audit"
        action = "parser_smoke_then_decide_download"
    elif score >= 8:
        tier = "P3_auxiliary_or_intervention"
        action = "keep_for_targeted_validation_not_primary_clock_training"
    else:
        tier = "P4_low_priority"
        action = "do_not_download_until_main_gap_is_fixed"
    return tier, action, ";".join(dict.fromkeys(blockers))


def rough_esummary_score(accession: str, row: dict, source_terms_for_accession: list[str]) -> int:
    if accession in SEED_ACCESSIONS:
        return 10_000
    text = f"{row.get('title','')} {row.get('summary','')} {row.get('experiment_type','')} {row.get('supp_file','')}".lower()
    score = 0
    if row.get("taxon") == "Mus musculus":
        score += 3
    if "methylation profiling" in text:
        score += 4
    if "rrbs" in text or "reduced representation" in text:
        score += 4
    if "bisulfite" in text or "bs-seq" in text or "wgbs" in text:
        score += 2
    if "cov" in text:
        score += 3
    for needle in ["aging", "aged", " old ", "24 month", "26 month", "30 month", "rejuvenation", "epigenetic clock"]:
        if needle in f" {text} ":
            score += 2
    if "month" in text:
        score += 1
    if source_terms_for_accession:
        score += min(3, len(source_terms_for_accession))
    if "mus musculus" not in str(row.get("taxon", "")).lower():
        score -= 5
    if "methylation" not in text:
        score -= 8
    return score


def build_candidate(
    gse: str,
    source_terms: list[str],
    eutils_summary: dict | None,
    refresh_soft: bool,
    do_head: bool,
) -> tuple[dict, list[dict], list[dict]]:
    eutils_summary = eutils_summary or eutils_summary_for_accession(gse)
    soft = cached_soft_text(gse, refresh=refresh_soft)
    soft_series: dict[str, str] = {}
    soft_samples: list[dict[str, str]] = []
    if soft.text:
        soft_series, soft_samples = parse_soft(soft.text)
    sample_summary, sample_rows = summarize_samples(gse, soft_samples)
    supplements = fetch_supplements(gse, do_head=do_head)
    supp_summary = summarize_supplements(supplements)
    schema_guess, schema_blocker = infer_schema(eutils_summary, supplements)

    summary_text = soft_series.get("Series_summary") or eutils_summary.get("summary", "")
    title = soft_series.get("Series_title") or eutils_summary.get("title", "")
    row = {
        "dataset": gse,
        "geo_url": geo_url(gse),
        "soft_url": soft_url(gse),
        "supplement_dir_url": suppl_dir_url(gse),
        "title": title,
        "summary": summary_text[:1000],
        "overall_design": soft_series.get("Series_overall_design", "")[:1000],
        "taxon": eutils_summary.get("taxon", ""),
        "entry_type": eutils_summary.get("entry_type", ""),
        "experiment_type": eutils_summary.get("experiment_type", soft_series.get("Series_type", "")),
        "public_date": eutils_summary.get("public_date", ""),
        "supp_file_eutils": eutils_summary.get("supp_file", ""),
        "source_queries": " || ".join(source_terms),
        "soft_fetch_error": soft.error,
        "preflight_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        **sample_summary,
        **supp_summary,
        "processed_schema_guess": schema_guess,
        "schema_blocker": schema_blocker,
    }
    tier, action, blockers = rank_candidate(row)
    row["priority_tier"] = tier
    row["recommended_action"] = action
    row["blockers"] = ";".join(item for item in [blockers, schema_blocker] if item)
    return row, supplements, sample_rows


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = sorted({key for row in rows for key in row.keys()})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--retmax", type=int, default=35)
    parser.add_argument("--out_prefix", default=str(META_DIR / "geo_old_age_candidate"))
    parser.add_argument("--accessions", default=",".join(SEED_ACCESSIONS))
    parser.add_argument("--max_candidates", type=int, default=28)
    parser.add_argument("--no_eutils_search", action="store_true")
    parser.add_argument("--refresh_soft", action="store_true")
    parser.add_argument("--no_head", action="store_true")
    args = parser.parse_args()

    source_terms: dict[str, list[str]] = defaultdict(list)
    summaries: dict[str, dict] = {}
    accessions = {item.strip() for item in args.accessions.split(",") if item.strip()}
    for accession in accessions:
        source_terms[accession].append("manual_seed_v7_9")

    if not args.no_eutils_search:
        found, found_terms, found_summaries = discover_accessions(DEFAULT_SEARCH_TERMS, args.retmax)
        accessions.update(found)
        summaries.update(found_summaries)
        for accession, terms in found_terms.items():
            source_terms[accession].extend(terms)

    for accession in list(accessions):
        if accession not in summaries:
            try:
                summaries[accession] = eutils_summary_for_accession(accession)
            except Exception:
                summaries[accession] = {"dataset": accession}

    ranked_accessions = sorted(
        accessions,
        key=lambda value: (
            -rough_esummary_score(value, summaries.get(value, {}), source_terms.get(value, [])),
            int(value.removeprefix("GSE")),
        ),
    )
    if args.max_candidates and len(ranked_accessions) > args.max_candidates:
        kept = set(ranked_accessions[: args.max_candidates])
        kept.update(SEED_ACCESSIONS)
        ranked_accessions = [accession for accession in ranked_accessions if accession in kept]
    print(f"[v7.9] Detailed SOFT/filelist preflight for {len(ranked_accessions)} candidates", flush=True)

    candidate_rows = []
    supplement_rows = []
    sample_rows = []
    for gse in ranked_accessions:
        try:
            row, supplements, samples = build_candidate(
                gse,
                source_terms=source_terms.get(gse, []),
                eutils_summary=summaries.get(gse),
                refresh_soft=args.refresh_soft,
                do_head=not args.no_head,
            )
        except Exception as exc:
            row = {
                "dataset": gse,
                "geo_url": geo_url(gse),
                "preflight_timestamp_utc": datetime.now(timezone.utc).isoformat(),
                "priority_tier": "ERROR",
                "recommended_action": "inspect_manually",
                "blockers": str(exc)[:500],
                "source_queries": " || ".join(source_terms.get(gse, [])),
            }
            supplements = []
            samples = []
        candidate_rows.append(row)
        supplement_rows.extend(supplements)
        sample_rows.extend(samples)
        print(
            f"[v7.9] {gse}: tier={row.get('priority_tier')} "
            f"max_age_months={row.get('max_age_months')} tissues={row.get('tissue_count')} "
            f"schema={row.get('processed_schema_guess')} action={row.get('recommended_action')}"
            ,
            flush=True,
        )
        time.sleep(0.1)

    out_prefix = Path(args.out_prefix)
    write_csv(out_prefix.with_name(out_prefix.name + "_inventory.csv"), candidate_rows)
    write_csv(out_prefix.with_name(out_prefix.name + "_supplements.csv"), supplement_rows)
    write_csv(out_prefix.with_name(out_prefix.name + "_samples.csv"), sample_rows)
    payload = {
        "official_source_urls": OFFICIAL_SOURCE_URLS,
        "search_terms": DEFAULT_SEARCH_TERMS,
        "seed_accessions": SEED_ACCESSIONS,
        "candidate_count": len(candidate_rows),
        "candidates": candidate_rows,
    }
    out_prefix.with_name(out_prefix.name + "_inventory.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[v7.9] Wrote {out_prefix.with_name(out_prefix.name + '_inventory.csv')}")


if __name__ == "__main__":
    main()

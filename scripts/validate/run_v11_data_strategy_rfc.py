#!/usr/bin/env python3
"""v11 data-strategy RFC controller for mouse methylation clocks.

This controller starts after v10 fails to find P1/P3 candidates in the current
public GEO processed/seed space. It broadens discovery through official PubMed,
GEO, and SRA E-Utils metadata, extracts accession leads, verifies new GSE leads
through GEO/SOFT/FTP preflight, and writes an RFC decision. It does not download
large files, build matrices, run FASTQ/Bismark, train models, or start
autoresearch.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")
OUT_DIR = ROOT / "results" / "ralph_v11_data_strategy"
REPORT = OUT_DIR / "v11_data_strategy_rfc_report.md"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

TARGET_TISSUES = {"brain_cortex", "heart", "lung"}
OLD_THRESHOLD_WEEKS = 104.0

OFFICIAL_SOURCE_DOCS = {
    "geo_download": "https://www.ncbi.nlm.nih.gov/geo/info/download.html",
    "geo_programmatic_access": "https://www.ncbi.nlm.nih.gov/geo/info/geo_paccess.html",
    "geo_soft": "https://www.ncbi.nlm.nih.gov/geo/info/soft.html",
    "sra_download": "https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/",
    "ena_browser_api": "https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/browser-api.html",
    "ncbi_eutils": "https://www.ncbi.nlm.nih.gov/books/NBK25501/",
}

PUBMED_QUERIES = [
    '("Mus musculus" OR mouse OR mice) AND (aging OR aged OR old) AND ("DNA methylation" OR methylome) AND (RRBS OR WGBS OR bisulfite) AND (cortex OR "brain cortex" OR heart OR lung)',
    '("Mus musculus" OR mouse OR mice) AND ("epigenetic clock" OR "methylation clock") AND (RRBS OR WGBS OR bisulfite OR methylome)',
    '("Mus musculus" OR mouse OR mice) AND ("24 month" OR "26 month" OR "28 month" OR "30 month") AND ("DNA methylation" OR methylome)',
    '("Mus musculus" OR mouse OR mice) AND (rejuvenation OR parabiosis OR senescence) AND ("DNA methylation" OR methylome) AND (cortex OR heart OR lung OR brain)',
]

SRA_QUERIES = [
    'Mus musculus[Organism] AND ("bisulfite" OR "RRBS" OR "WGBS" OR "methylation") AND ("aging" OR "aged" OR "old") AND ("brain cortex" OR cortex OR heart OR lung)',
    'Mus musculus[Organism] AND ("DNA methylation" OR methylome) AND ("24 month" OR "26 month" OR "28 month" OR "30 month")',
    'Mus musculus[Organism] AND ("reduced representation bisulfite" OR "whole genome bisulfite") AND (aging OR aged OR old)',
]

ACCESSION_RE = re.compile(
    r"\b("
    r"GSE\d+|GSM\d+|GPL\d+|"
    r"SRP\d+|SRX\d+|SRS\d+|SRR\d+|SRA\d+|"
    r"ERP\d+|ERX\d+|ERS\d+|ERR\d+|"
    r"DRP\d+|DRX\d+|DRS\d+|DRR\d+|"
    r"PRJNA\d+|PRJEB\d+|PRJDB\d+"
    r")\b",
    re.IGNORECASE,
)

METHYLATION_TERMS = [
    "rrbs",
    "wgbs",
    "bisulfite",
    "methylome",
    "dna methylation",
    "methylation profiling",
]
OLD_TERMS = ["aging", "aged", "old", "24 month", "26 month", "28 month", "30 month", "lifespan"]
TARGET_TERMS = ["brain cortex", "cortex", "heart", "lung"]
BLOCKER_TERMS = ["single-cell", "single cell", "organoid", "cell line", "in vitro", "ipcrtag", "itag", "hsc"]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def import_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def discovery_module():
    return import_module(ROOT / "scripts" / "etl" / "16_geo_old_age_candidate_discovery.py", "geo_old_age_candidate_discovery_v11")


def v9_module():
    return import_module(ROOT / "scripts" / "validate" / "run_v9_ralph_loop.py", "run_v9_ralph_loop_helpers_v11")


def v10_module():
    return import_module(ROOT / "scripts" / "validate" / "run_v10_ralph_loop.py", "run_v10_ralph_loop_helpers_v11")


def eutils_get(endpoint: str, params: dict, timeout: int = 45) -> str:
    query = urllib.parse.urlencode(params)
    request = urllib.request.Request(
        EUTILS + endpoint + "?" + query,
        headers={"User-Agent": "mouse-methyl-v11-data-strategy/1.0"},
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read().decode("utf-8", errors="replace")


def fetch_url_text(url: str, timeout: int = 45) -> str:
    request = urllib.request.Request(url, headers={"User-Agent": "mouse-methyl-v11-data-strategy/1.0"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read().decode("utf-8", errors="replace")


def eutils_search(db: str, term: str, retmax: int) -> list[str]:
    xml = eutils_get("esearch.fcgi", {"db": db, "term": term, "retmax": retmax, "usehistory": "n"})
    root = ET.fromstring(xml)
    return [node.text for node in root.findall(".//Id") if node.text]


def eutils_fetch_pubmed(ids: list[str]) -> list[dict]:
    rows = []
    for start in range(0, len(ids), 80):
        chunk = ids[start : start + 80]
        if not chunk:
            continue
        xml = eutils_get("efetch.fcgi", {"db": "pubmed", "id": ",".join(chunk), "retmode": "xml"})
        root = ET.fromstring(xml)
        for article in root.findall(".//PubmedArticle"):
            pmid = article.findtext(".//PMID") or ""
            title_node = article.find(".//ArticleTitle")
            title = " ".join("".join(title_node.itertext()).split()) if title_node is not None else ""
            abstract_parts = []
            for abstract_node in article.findall(".//AbstractText"):
                abstract_text = " ".join("".join(abstract_node.itertext()).split())
                if abstract_text:
                    abstract_parts.append(abstract_text)
            abstract = " ".join(abstract_parts)
            journal = article.findtext(".//Journal/Title") or ""
            year = (
                article.findtext(".//JournalIssue/PubDate/Year")
                or article.findtext(".//DateCompleted/Year")
                or article.findtext(".//PubDate/Year")
                or ""
            )
            rows.append(
                {
                    "source_db": "pubmed",
                    "source_id": pmid,
                    "title": title,
                    "abstract_or_summary": abstract,
                    "journal": journal,
                    "year": year,
                    "source_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else "",
                }
            )
        time.sleep(0.34)
    return rows


def eutils_summary_sra(ids: list[str]) -> list[dict]:
    rows = []
    for start in range(0, len(ids), 80):
        chunk = ids[start : start + 80]
        if not chunk:
            continue
        xml = eutils_get("esummary.fcgi", {"db": "sra", "id": ",".join(chunk), "retmode": "xml"})
        root = ET.fromstring(xml)
        for docsum in root.findall(".//DocSum"):
            row = defaultdict(str)
            row["source_db"] = "sra"
            row["source_id"] = docsum.findtext("Id") or ""
            parts = []
            for item in docsum.findall("Item"):
                name = item.attrib.get("Name", "")
                text = "".join(item.itertext())
                if name:
                    row[name] = text
                parts.append(text)
            title = row.get("Title") or row.get("Study", "")
            rows.append(
                {
                    "source_db": "sra",
                    "source_id": row["source_id"],
                    "title": " ".join(str(title).split()),
                    "abstract_or_summary": " ".join(" ".join(parts).split())[:5000],
                    "journal": "",
                    "year": "",
                    "source_url": f"https://www.ncbi.nlm.nih.gov/sra/?term={row['source_id']}" if row["source_id"] else "",
                }
            )
        time.sleep(0.34)
    return rows


def extract_accessions(text: str) -> list[str]:
    return sorted({match.group(1).upper() for match in ACCESSION_RE.finditer(text or "")})


def keyword_hits(text: str, terms: list[str]) -> list[str]:
    lower = (text or "").lower()
    return [term for term in terms if term in lower]


def score_lead(text: str, accessions: list[str]) -> tuple[int, list[str]]:
    lower = (text or "").lower()
    reasons = []
    score = 0
    methyl_hits = keyword_hits(lower, METHYLATION_TERMS)
    old_hits = keyword_hits(lower, OLD_TERMS)
    target_hits = keyword_hits(lower, TARGET_TERMS)
    blocker_hits = keyword_hits(lower, BLOCKER_TERMS)
    if methyl_hits:
        score += 3
        reasons.append("methylation_context")
    if old_hits:
        score += 3
        reasons.append("old_age_context")
    if target_hits:
        score += 3
        reasons.append("target_tissue_context:" + ",".join(target_hits[:4]))
    if any(acc.startswith("GSE") for acc in accessions):
        score += 2
        reasons.append("geo_accession")
    if any(acc.startswith(("SRP", "SRX", "SRR", "SRA", "PRJ")) for acc in accessions):
        score += 1
        reasons.append("raw_accession")
    if blocker_hits:
        score -= 3
        reasons.append("blocker_terms:" + ",".join(blocker_hits[:4]))
    return score, reasons


def literature_search(retmax_pubmed: int, retmax_sra: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    source_rows = []
    network_rows = []
    seen_pubmed_ids = set()
    for query in PUBMED_QUERIES:
        try:
            ids = [uid for uid in eutils_search("pubmed", query, retmax_pubmed) if uid not in seen_pubmed_ids]
            seen_pubmed_ids.update(ids)
            source_rows.extend(eutils_fetch_pubmed(ids))
            network_rows.append({"timestamp": utc_now(), "source": "pubmed", "query": query, "n_ids": len(ids), "error": ""})
        except Exception as exc:
            network_rows.append({"timestamp": utc_now(), "source": "pubmed", "query": query, "n_ids": 0, "error": str(exc)[:500]})
    seen_sra_ids = set()
    for query in SRA_QUERIES:
        try:
            ids = [uid for uid in eutils_search("sra", query, retmax_sra) if uid not in seen_sra_ids]
            seen_sra_ids.update(ids)
            source_rows.extend(eutils_summary_sra(ids))
            network_rows.append({"timestamp": utc_now(), "source": "sra", "query": query, "n_ids": len(ids), "error": ""})
        except Exception as exc:
            network_rows.append({"timestamp": utc_now(), "source": "sra", "query": query, "n_ids": 0, "error": str(exc)[:500]})
    return pd.DataFrame(source_rows), pd.DataFrame(network_rows)


def build_lead_table(source_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in source_df.iterrows():
        text = " ".join(str(row.get(field, "")) for field in ["title", "abstract_or_summary"])
        accessions = extract_accessions(text)
        score, reasons = score_lead(text, accessions)
        rows.append(
            {
                "source_db": row.get("source_db", ""),
                "source_id": row.get("source_id", ""),
                "title": row.get("title", ""),
                "year": row.get("year", ""),
                "source_url": row.get("source_url", ""),
                "accessions": ";".join(accessions),
                "gse_accessions": ";".join(acc for acc in accessions if acc.startswith("GSE")),
                "raw_accessions": ";".join(
                    acc
                    for acc in accessions
                    if acc.startswith(("SRP", "SRX", "SRR", "SRA", "ERP", "ERX", "ERR", "DRP", "DRX", "DRR", "PRJ"))
                ),
                "lead_score": score,
                "lead_reasons": ";".join(reasons),
                "target_hits": ";".join(keyword_hits(text, TARGET_TERMS)),
                "old_hits": ";".join(keyword_hits(text, OLD_TERMS)),
                "methylation_hits": ";".join(keyword_hits(text, METHYLATION_TERMS)),
                "blocker_hits": ";".join(keyword_hits(text, BLOCKER_TERMS)),
                "snippet": text[:800],
            }
        )
    leads = pd.DataFrame(rows)
    if leads.empty:
        return leads
    leads = leads.sort_values(["lead_score", "source_db", "source_id"], ascending=[False, True, True]).reset_index(drop=True)
    return leads


def known_gse_accessions() -> set[str]:
    paths = [
        ROOT / "results" / "ralph_v10_loop" / "candidate_gate_table.csv",
        ROOT / "results" / "ralph_v9_loop" / "candidate_gate_table.csv",
        ROOT / "metadata" / "geo_old_age_candidate_inventory.csv",
    ]
    known = set()
    for path in paths:
        if not path.exists():
            continue
        df = pd.read_csv(path)
        if "dataset" in df:
            known.update(str(value).upper() for value in df["dataset"].dropna())
    return known


def resolve_gsm_to_gse(gsm_accessions: list[str], max_resolve: int) -> pd.DataFrame:
    rows = []
    for gsm in gsm_accessions[:max_resolve]:
        url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}&targ=self&view=full&form=text"
        try:
            text = fetch_url_text(url)
            gse_values = sorted(set(re.findall(r"Sample_series_id\s*=\s*(GSE\d+)", text)))
            error = ""
        except Exception as exc:
            gse_values = []
            error = str(exc)[:500]
        rows.append(
            {
                "gsm": gsm,
                "resolved_gse": ";".join(gse_values),
                "n_resolved_gse": len(gse_values),
                "source_url": url,
                "error": error,
            }
        )
        time.sleep(0.34)
    return pd.DataFrame(rows)


def verify_new_gse_leads(gse_leads: list[str], max_verify: int, run_smoke: bool) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if not gse_leads:
        empty = pd.DataFrame()
        return empty, empty, empty, empty
    discovery = discovery_module()
    v9 = v9_module()
    v10 = v10_module()
    v84 = v9.v84_module()
    v8 = v9.v8_loop_module()

    inventory_rows = []
    sample_rows = []
    supplement_rows = []
    for gse in gse_leads[:max_verify]:
        try:
            row, supplements, samples = discovery.build_candidate(
                gse,
                source_terms=["v11_pubmed_or_sra_lead"],
                eutils_summary=None,
                refresh_soft=False,
                do_head=True,
            )
        except Exception as exc:
            row = {
                "dataset": gse,
                "geo_url": discovery.geo_url(gse),
                "priority_tier": "ERROR",
                "recommended_action": "inspect_manually",
                "blockers": str(exc)[:500],
            }
            supplements = []
            samples = []
        inventory_rows.append(row)
        sample_rows.extend(samples)
        supplement_rows.extend(supplements)
        time.sleep(0.1)

    inventory = pd.DataFrame(inventory_rows)
    samples = pd.DataFrame(sample_rows)
    supplements = pd.DataFrame(supplement_rows)
    gate_rows = []
    smoke_rows = []
    network_log = OUT_DIR / "network_resolution_log.jsonl"
    for dataset in inventory.get("dataset", pd.Series(dtype=str)).dropna().astype(str).unique():
        row, matched, smoke = v10.summarize_candidate(
            dataset,
            inventory,
            samples,
            supplements,
            v9,
            v84,
            v8,
            OUT_DIR,
            run_smoke,
            network_log,
        )
        row["v11_verification_source"] = "pubmed_or_sra_accession_lead"
        gate_rows.append(row)
        if not smoke.empty:
            smoke_rows.extend(smoke.to_dict(orient="records"))
        matched.to_csv(OUT_DIR / f"{dataset}_v11_matched_supplement_files.csv", index=False)
    return inventory, samples, supplements, pd.DataFrame(gate_rows), pd.DataFrame(smoke_rows)


def decide_state(leads: pd.DataFrame, gate: pd.DataFrame, raw_leads: pd.DataFrame) -> dict:
    n_p1 = int(gate.get("candidate_tier", pd.Series(dtype=str)).astype(str).str.startswith("P1").sum()) if not gate.empty else 0
    n_p3 = int(gate.get("candidate_tier", pd.Series(dtype=str)).astype(str).str.startswith("P3").sum()) if not gate.empty else 0
    n_raw_high = int((raw_leads.get("lead_score", pd.Series(dtype=float)) >= 7).sum()) if not raw_leads.empty else 0
    if n_p1 or n_p3:
        status = "continue"
        decision = "v11_candidate_found_promote_to_smoke_or_minimal_etl_plan"
        next_action = "Run explicit adapter smoke for P1/P3 leads, then matrix gate if smoke passes."
    elif n_raw_high:
        status = "rfc_required"
        decision = "v11_raw_literature_leads_need_manual_biosample_verification"
        next_action = "Open a minimal ETL RFC for high-scoring raw accessions; verify BioSample age/tissue before FASTQ pilot."
    else:
        status = "failure"
        decision = "v11_no_actionable_public_or_raw_target_tissue_leads"
        next_action = "Stop model work; pursue external data acquisition, new experiment, or benchmark redefinition."
    return {
        "timestamp": utc_now(),
        "loop_version": "v11.0",
        "status": status,
        "loop_decision": decision,
        "metrics": {
            "n_literature_or_sra_records": int(len(leads)),
            "n_high_scoring_raw_leads": n_raw_high,
            "n_verified_new_gse": int(len(gate)),
            "n_p1": n_p1,
            "n_p3": n_p3,
        },
        "next_action": next_action,
        "success_gate": {
            "candidate_gate": "P1 or P3 with target old tissue and traceable methylation data",
            "matrix_gate": "common 5kb regions >=50000",
            "benchmark_gate": "old104+ MAE improves >=10w and sanity checks pass",
        },
        "failure_gate": {
            "no_actionable_leads": "No P1/P3 and no high-scoring raw leads after PubMed/SRA official refresh",
            "raw_lead_requires_rfc": "Raw-only leads require separate BioSample verification and pilot plan",
        },
    }


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


def write_report(state: dict, leads: pd.DataFrame, gate: pd.DataFrame, raw_leads: pd.DataFrame) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    tier_counts = {}
    if not gate.empty and "candidate_tier" in gate:
        tier_counts = Counter(gate["candidate_tier"].fillna("")).most_common()
    gate_view_cols = [
        "dataset",
        "candidate_tier",
        "n_old_target_samples",
        "old_target_tissues",
        "n_old_adjacent_samples",
        "old_adjacent_tissues",
        "n_matched_processed_files",
        "schema_smoke_pass",
        "assembly_compatible",
        "smoke_common_regions_estimate",
        "minimal_fastq_pilot_candidate",
        "auxiliary_or_blocker_reason",
        "next_step",
    ]
    lead_view_cols = [
        "source_db",
        "source_id",
        "lead_score",
        "gse_accessions",
        "raw_accessions",
        "target_hits",
        "old_hits",
        "blocker_hits",
        "title",
    ]
    lines = [
        "# v11 Data Strategy RFC Report",
        "",
        f"Date: {utc_now()}",
        "",
        "## Summary",
        "",
        "v11 broadened discovery beyond the v10 GEO processed/seed space using official PubMed and SRA E-Utils metadata. It did not download large files, build matrices, run FASTQ/Bismark, train models, or start autoresearch.",
        "",
        f"- Status: `{state['status']}`",
        f"- Decision: `{state['loop_decision']}`",
        f"- Next action: {state['next_action']}",
        "",
        "## Metrics",
        "",
        md_table(pd.DataFrame([state["metrics"]])),
        "",
        "## Verified New GSE Gate Counts",
        "",
        md_table(pd.DataFrame(tier_counts, columns=["candidate_tier", "count"])),
        "",
        "## Verified New GSE Candidates",
        "",
        md_table(gate[[col for col in gate_view_cols if col in gate.columns]], max_rows=40),
        "",
        "## Top Literature/SRA Leads",
        "",
        md_table(leads[[col for col in lead_view_cols if col in leads.columns]], max_rows=30),
        "",
        "## Raw Leads Requiring RFC",
        "",
        md_table(raw_leads[[col for col in lead_view_cols if col in raw_leads.columns]], max_rows=30),
        "",
        "## Official Source Rules",
        "",
        md_table(pd.DataFrame([{"source": key, "url": value} for key, value in OFFICIAL_SOURCE_DOCS.items()])),
        "",
        "## Interpretation",
        "",
    ]
    if state["status"] == "continue":
        lines.extend(
            [
                "At least one P1/P3 candidate was found. The next step is explicit adapter smoke only; full downloads and training remain blocked until the matrix gate passes.",
            ]
        )
    elif state["status"] == "rfc_required":
        lines.extend(
            [
                "No verified new GSE P1/P3 candidate was found, but high-scoring raw accession leads exist. These are not enough to start FASTQ/Bismark. A minimal ETL RFC must first verify BioSample-level age, tissue, strain, assay, and run count.",
            ]
        )
    else:
        lines.extend(
            [
                "No actionable public or raw target-tissue leads were found in this official metadata refresh. Continuing model tuning on the same matrix is not justified; the project needs external data acquisition, new experiments, or benchmark redefinition.",
            ]
        )
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def append_jsonl(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(payload, sort_keys=True) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--retmax_pubmed", type=int, default=80)
    parser.add_argument("--retmax_sra", type=int, default=80)
    parser.add_argument("--max_verify_gse", type=int, default=25)
    parser.add_argument("--max_resolve_gsm", type=int, default=80)
    parser.add_argument("--run_smoke", action="store_true", help="Allow small adapter smoke downloads for verified new GSE leads.")
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    append_jsonl(
        OUT_DIR / "ralph_iteration_log.jsonl",
        {
            "timestamp": utc_now(),
            "loop": "v11_data_strategy",
            "action": "start",
            "run_smoke": bool(args.run_smoke),
            "retmax_pubmed": args.retmax_pubmed,
            "retmax_sra": args.retmax_sra,
            "max_resolve_gsm": args.max_resolve_gsm,
        },
    )
    sources, network = literature_search(args.retmax_pubmed, args.retmax_sra)
    sources.to_csv(OUT_DIR / "literature_sra_source_records.csv", index=False)
    network.to_csv(OUT_DIR / "network_query_log.csv", index=False)
    for record in network.to_dict(orient="records"):
        append_jsonl(OUT_DIR / "network_resolution_log.jsonl", record)

    leads = build_lead_table(sources)
    leads.to_csv(OUT_DIR / "literature_sra_accession_leads.csv", index=False)

    known = known_gse_accessions()
    extracted_gse = []
    extracted_gsm = []
    for values in leads.get("gse_accessions", pd.Series(dtype=str)).fillna(""):
        extracted_gse.extend([value for value in str(values).split(";") if value])
    for values in leads.get("accessions", pd.Series(dtype=str)).fillna(""):
        extracted_gsm.extend([value for value in str(values).split(";") if value.startswith("GSM")])

    gsm_resolution = resolve_gsm_to_gse(sorted(set(extracted_gsm)), args.max_resolve_gsm)
    gsm_resolution.to_csv(OUT_DIR / "gsm_to_gse_resolution.csv", index=False)
    for values in gsm_resolution.get("resolved_gse", pd.Series(dtype=str)).fillna(""):
        extracted_gse.extend([value for value in str(values).split(";") if value])

    new_gse = sorted({value for value in extracted_gse if value not in known})

    inventory, samples, supplements, gate, smoke = verify_new_gse_leads(new_gse, args.max_verify_gse, args.run_smoke)
    inventory.to_csv(OUT_DIR / "verified_new_gse_inventory.csv", index=False)
    samples.to_csv(OUT_DIR / "verified_new_gse_samples.csv", index=False)
    supplements.to_csv(OUT_DIR / "verified_new_gse_supplements.csv", index=False)
    gate.to_csv(OUT_DIR / "verified_new_gse_gate_table.csv", index=False)
    pd.DataFrame(smoke).to_csv(OUT_DIR / "candidate_smoke_manifest.csv", index=False)

    raw_leads = leads[
        (leads.get("raw_accessions", pd.Series(dtype=str)).fillna("").astype(str) != "")
        & (leads.get("lead_score", pd.Series(dtype=float)) >= 7)
    ].copy()
    raw_leads.to_csv(OUT_DIR / "raw_accession_rfc_leads.csv", index=False)

    state = decide_state(leads, gate, raw_leads)
    write_json(OUT_DIR / "ralph_decision_state.json", state)
    write_report(state, leads, gate, raw_leads)
    append_jsonl(
        OUT_DIR / "ralph_iteration_log.jsonl",
        {
            "timestamp": utc_now(),
            "loop": "v11_data_strategy",
            "action": "decision",
            "status": state["status"],
            "decision": state["loop_decision"],
            "metrics": state["metrics"],
        },
    )
    print(json.dumps(state, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Summarize v11 raw accession leads into minimal-ETL RFC project candidates."""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")
V11_DIR = ROOT / "results" / "ralph_v11_data_strategy"
OUT = V11_DIR / "raw_project_rfc_summary.csv"

INTEGRATED_GSE = {"GSE120137", "GSE80672", "GSE93957", "GSE121141", "GSE60012", "GSE213628"}
TARGET_TISSUES = {"brain cortex", "cortex", "heart", "lung"}
BLOCKER_TERMS = {"single-cell", "single cell", "organoid", "cell line", "in vitro", "ipcrtag", "itag", "hsc"}


def split_values(value: object) -> list[str]:
    if pd.isna(value):
        return []
    return [item for item in str(value).split(";") if item]


def extract(pattern: str, text: str) -> list[str]:
    return sorted(set(re.findall(pattern, text or "", flags=re.IGNORECASE)))


def load_gsm_resolution() -> dict[str, list[str]]:
    path = V11_DIR / "gsm_to_gse_resolution.csv"
    if not path.exists():
        return {}
    df = pd.read_csv(path)
    mapping: dict[str, list[str]] = {}
    for _, row in df.iterrows():
        mapping[str(row.get("gsm", "")).upper()] = split_values(row.get("resolved_gse", ""))
    return mapping


def resolved_gse_for_lead(row: pd.Series, gsm_map: dict[str, list[str]]) -> list[str]:
    text = " ".join(str(row.get(field, "")) for field in ["accessions", "raw_accessions", "gse_accessions", "snippet"])
    gsm_values = extract(r"GSM\d+", text)
    resolved = []
    for gsm in gsm_values:
        resolved.extend(gsm_map.get(gsm.upper(), []))
    return sorted(set(resolved))


def project_keys(row: pd.Series, gsm_map: dict[str, list[str]]) -> list[str]:
    text = " ".join(str(row.get(field, "")) for field in ["accessions", "raw_accessions", "gse_accessions", "snippet"])
    gse = sorted(set(split_values(row.get("gse_accessions")) + resolved_gse_for_lead(row, gsm_map)))
    prjna = extract(r"PRJNA\d+", text)
    srp = extract(r"SRP\d+", text)
    keys = gse + prjna + srp
    return keys or [str(row.get("source_id", ""))]


def decide_priority(row: dict) -> str:
    linked_gse = set(split_values(row.get("linked_gse", "")))
    blockers = set(split_values(row.get("blocker_hits", "")))
    tissues = set(split_values(row.get("target_hits", "")))
    old_hits = split_values(row.get("old_hits", ""))
    if linked_gse & INTEGRATED_GSE:
        return "already_integrated_or_previously_tested"
    if blockers:
        return "blocked_context_manual_only"
    if tissues & TARGET_TISSUES and old_hits and int(row.get("n_sra_records", 0)) >= 2:
        return "minimal_etl_rfc_candidate"
    return "manual_biosample_check_only"


def main() -> None:
    leads_path = V11_DIR / "raw_accession_rfc_leads.csv"
    if not leads_path.exists():
        raise SystemExit(f"Missing {leads_path}")
    leads = pd.read_csv(leads_path)
    gsm_map = load_gsm_resolution()
    rows = []
    for _, lead in leads.iterrows():
        lead_gse = sorted(set(split_values(lead.get("gse_accessions")) + resolved_gse_for_lead(lead, gsm_map)))
        for key in project_keys(lead, gsm_map):
            rows.append(
                {
                    "project_key": key.upper(),
                    "source_id": lead.get("source_id", ""),
                    "linked_gse": ";".join(lead_gse),
                    "raw_accessions": lead.get("raw_accessions", ""),
                    "target_hits": lead.get("target_hits", ""),
                    "old_hits": lead.get("old_hits", ""),
                    "methylation_hits": lead.get("methylation_hits", ""),
                    "blocker_hits": lead.get("blocker_hits", ""),
                    "lead_score": lead.get("lead_score", 0),
                    "snippet": lead.get("snippet", ""),
                }
            )
    expanded = pd.DataFrame(rows)
    summary_rows = []
    for key, group in expanded.groupby("project_key", dropna=False):
        linked_gse = sorted(set(item for values in group["linked_gse"].fillna("") for item in split_values(values)))
        raw_accessions = sorted(set(item for values in group["raw_accessions"].fillna("") for item in split_values(values)))
        target_hits = sorted(set(item for values in group["target_hits"].fillna("") for item in split_values(values)))
        old_hits = sorted(set(item for values in group["old_hits"].fillna("") for item in split_values(values)))
        methylation_hits = sorted(set(item for values in group["methylation_hits"].fillna("") for item in split_values(values)))
        blocker_hits = sorted(set(item for values in group["blocker_hits"].fillna("") for item in split_values(values)))
        snippets = " || ".join(str(value)[:250] for value in group["snippet"].dropna().head(4))
        row = {
            "project_key": key,
            "n_sra_records": int(group["source_id"].nunique()),
            "max_lead_score": float(pd.to_numeric(group["lead_score"], errors="coerce").max()),
            "linked_gse": ";".join(linked_gse),
            "n_raw_accessions": len(raw_accessions),
            "raw_accession_examples": ";".join(raw_accessions[:10]),
            "target_hits": ";".join(target_hits),
            "old_hits": ";".join(old_hits),
            "methylation_hits": ";".join(methylation_hits),
            "blocker_hits": ";".join(blocker_hits),
            "evidence_snippet": snippets,
        }
        row["pilot_priority"] = decide_priority(row)
        summary_rows.append(row)
    summary = pd.DataFrame(summary_rows)
    if not summary.empty:
        priority_order = {
            "minimal_etl_rfc_candidate": 0,
            "manual_biosample_check_only": 1,
            "blocked_context_manual_only": 2,
            "already_integrated_or_previously_tested": 3,
        }
        summary["_rank"] = summary["pilot_priority"].map(priority_order).fillna(99)
        summary = summary.sort_values(["_rank", "max_lead_score", "n_sra_records"], ascending=[True, False, False]).drop(columns=["_rank"])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(OUT, index=False)
    print(summary.head(30).to_string(index=False))


if __name__ == "__main__":
    main()

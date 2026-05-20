#!/usr/bin/env python3
"""Run v8.5 targeted GEO refresh and GSE225166 adapter audit.

This is a preflight/audit step only. It uses official GEO/E-Utils/SOFT/filelist
metadata and does not download large supplement archives or run model training.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")
OUT_DIR = ROOT / "results" / "v8_5_targeted_geo_refresh"
REPORT = ROOT / "doc" / "20_analysis" / "20_20260518_v8_5_targeted_geo_refresh_gse225166_audit_report.md"

TARGETED_SEARCH_TERMS = [
    'Mus musculus[ORGN] AND "Methylation profiling by high throughput sequencing" AND (RRBS OR "reduced representation" OR bisulfite) AND (aged OR aging OR old OR "24 month" OR "26 month" OR "30 month") AND cortex AND gse[ETYP]',
    'Mus musculus[ORGN] AND "Methylation profiling by high throughput sequencing" AND (RRBS OR "reduced representation" OR bisulfite) AND (aged OR aging OR old OR "24 month" OR "26 month" OR "30 month") AND heart AND gse[ETYP]',
    'Mus musculus[ORGN] AND "Methylation profiling by high throughput sequencing" AND (RRBS OR "reduced representation" OR bisulfite) AND (aged OR aging OR old OR "24 month" OR "26 month" OR "30 month") AND lung AND gse[ETYP]',
    'Mus musculus[ORGN] AND "DNA methylation" AND ("brain cortex" OR "frontal cortex" OR "prefrontal cortex" OR hippocampus) AND (aged OR aging OR "24 month" OR "26 month" OR old) AND gse[ETYP]',
    'Mus musculus[ORGN] AND "DNA methylation" AND ("heart" OR "cardiac" OR "lung") AND (aged OR aging OR "24 month" OR "26 month" OR old) AND gse[ETYP]',
]

SEED_ACCESSIONS = [
    "GSE225166",
    "GSE232547",
    "GSE171236",
    "GSE138368",
    "GSE151541",
    "GSE169234",
    "GSE103249",
    "GSE108762",
    "GSE281602",
]

INTEGRATED_DATASETS = {"GSE120137", "GSE80672", "GSE93957", "GSE121141", "GSE60012", "GSE213628"}


def import_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def discovery_module():
    return import_module(ROOT / "scripts" / "etl" / "16_geo_old_age_candidate_discovery.py", "geo_old_age_candidate_discovery_v85")


def v84_module():
    return import_module(ROOT / "scripts" / "validate" / "v8_4_old_tissue_candidate_and_calibration.py", "v8_4_helpers_v85")


def read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path) if path.exists() else pd.DataFrame()


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


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


def discover_targeted_accessions(module, retmax: int) -> tuple[list[str], dict[str, list[str]], dict[str, dict]]:
    accessions: set[str] = set(SEED_ACCESSIONS)
    source_terms: dict[str, list[str]] = {accession: ["manual_seed_v8_5"] for accession in SEED_ACCESSIONS}
    summaries: dict[str, dict] = {}
    for term in TARGETED_SEARCH_TERMS:
        try:
            uids = module.eutils_search(term, retmax=retmax)
            rows = module.eutils_summaries_for_uids(uids)
        except Exception as exc:
            print(f"[v8.5][WARN] E-Utils search failed: {term!r}: {exc}", flush=True)
            continue
        for row in rows:
            accession = row.get("dataset", "")
            if not accession.startswith("GSE"):
                continue
            accessions.add(accession)
            source_terms.setdefault(accession, []).append(term)
            summaries[accession] = row
        time.sleep(0.34)
    for accession in sorted(accessions):
        if accession not in summaries:
            try:
                summaries[accession] = module.eutils_summary_for_accession(accession)
            except Exception as exc:
                summaries[accession] = {"dataset": accession, "eutils_error": str(exc)[:300]}
    ranked = sorted(
        accessions,
        key=lambda value: (
            -module.rough_esummary_score(value, summaries.get(value, {}), source_terms.get(value, [])),
            int(value.removeprefix("GSE")) if value.removeprefix("GSE").isdigit() else 999999999,
        ),
    )
    return ranked, source_terms, summaries


def preflight_candidates(module, accessions: list[str], source_terms: dict[str, list[str]], summaries: dict[str, dict], max_candidates: int) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    selected = list(dict.fromkeys(accessions[:max_candidates] + SEED_ACCESSIONS))
    candidate_rows: list[dict] = []
    sample_rows: list[dict] = []
    supplement_rows: list[dict] = []
    for accession in selected:
        try:
            row, supplements, samples = module.build_candidate(
                accession,
                source_terms=source_terms.get(accession, ["v8_5_targeted_search"]),
                eutils_summary=summaries.get(accession),
                refresh_soft=False,
                do_head=True,
            )
        except Exception as exc:
            row = {
                "dataset": accession,
                "geo_url": module.geo_url(accession),
                "priority_tier": "ERROR",
                "recommended_action": "inspect_manually",
                "blockers": str(exc)[:500],
                "source_queries": " || ".join(source_terms.get(accession, [])),
            }
            supplements = []
            samples = []
        candidate_rows.append(row)
        sample_rows.extend(samples)
        supplement_rows.extend(supplements)
        print(
            f"[v8.5] {accession}: tier={row.get('priority_tier')} "
            f"schema={row.get('processed_schema_guess')} max_age={row.get('max_age_weeks')}",
            flush=True,
        )
        time.sleep(0.05)
    return pd.DataFrame(candidate_rows), pd.DataFrame(sample_rows), pd.DataFrame(supplement_rows)


def enrich_candidate_ranking(v84, inventory: pd.DataFrame, samples: pd.DataFrame) -> pd.DataFrame:
    existing_inventory = read_csv(ROOT / "metadata" / "geo_old_age_candidate_inventory.csv")
    existing_samples = read_csv(ROOT / "metadata" / "geo_old_age_candidate_samples.csv")
    ranking = v84.build_candidate_ranking(existing_inventory, existing_samples, inventory, samples)
    if ranking.empty:
        return ranking
    ranking["v8_5_direct_mainline_ready"] = (
        ranking["v8_4_recommended_action"].eq("parser_smoke_then_download_decision")
        & (~ranking["already_integrated"].astype(bool))
        & (~ranking["v8_4_blockers"].astype(str).str.contains("targeted_celltype_or_non_bulk_context|superseries_use_subseries", na=False))
    )
    return ranking


def first_nonempty(*values: str) -> str:
    for value in values:
        if isinstance(value, str) and value.strip():
            return value.strip()
    return ""


def audit_gse225166(v84, samples: pd.DataFrame, supplements: pd.DataFrame, inventory: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    sample = samples[samples["dataset"].astype(str).eq("GSE225166")].copy()
    supp = supplements[supplements["dataset"].astype(str).eq("GSE225166")].copy()
    if sample.empty:
        cached = read_csv(OUT_DIR / "targeted_search_samples.csv")
        sample = cached[cached["dataset"].astype(str).eq("GSE225166")].copy()
    rows = []
    file_rows = supp[supp.get("filelist_type", pd.Series(dtype=str)).astype(str).eq("COV")].copy() if not supp.empty else pd.DataFrame()
    for _, row in sample.iterrows():
        sample_id = str(row.get("sample_id", ""))
        age = v84.conservative_sample_age_weeks(row)
        tissue_labels = v84.canonical_tissues(row.get("tissue_guess"))
        if not tissue_labels:
            tissue_labels = v84.canonical_tissues(
                f"{row.get('source_name', '')} {row.get('title', '')} {row.get('characteristics', '')}"
            )
        tissue = first_nonempty(";".join(sorted(tissue_labels)), str(row.get("tissue_guess", "")), str(row.get("source_name", "")))
        matched = file_rows[file_rows["supplement_name"].astype(str).str.contains(sample_id, regex=False)] if not file_rows.empty else pd.DataFrame()
        file_name = "" if matched.empty else str(matched.iloc[0].get("supplement_name", ""))
        file_size = np.nan if matched.empty else matched.iloc[0].get("filelist_size_bytes")
        text = f"{row.get('title', '')} {row.get('source_name', '')} {row.get('characteristics', '')} {file_name}".lower()
        rows.append(
            {
                "sample_id": sample_id,
                "title": row.get("title", ""),
                "source_name": row.get("source_name", ""),
                "tissue": tissue,
                "age_weeks": round(float(age), 3) if np.isfinite(age) else np.nan,
                "age_months": round(float(age) * 7.0 / 30.42, 3) if np.isfinite(age) else np.nan,
                "is_old_104w_plus": bool(np.isfinite(age) and age >= 104.0),
                "is_core_target_tissue": bool(set(str(tissue).split(";")) & {"brain_cortex", "heart", "lung"}),
                "has_cov_file": bool(not matched.empty),
                "cov_file_size_bytes": file_size,
                "supplement_name": file_name,
                "low_coverage_or_tag_context": bool("ipcrtag" in text or "itag" in text or "low coverage" in text),
                "single_cell_context": bool("single-cell" in text or "single cell" in text),
            }
        )
    audit = pd.DataFrame(rows)
    if audit.empty:
        summary = pd.DataFrame()
    else:
        summary = (
            audit.groupby(["tissue", "age_weeks"], dropna=False)
            .agg(
                n_samples=("sample_id", "count"),
                n_cov_files=("has_cov_file", "sum"),
                low_coverage_or_tag_rate=("low_coverage_or_tag_context", "mean"),
                single_cell_context_rate=("single_cell_context", "mean"),
            )
            .reset_index()
            .sort_values(["tissue", "age_weeks"])
        )
    inv = inventory[inventory["dataset"].astype(str).eq("GSE225166")]
    summary_text = " ".join(inv.get("summary", pd.Series(dtype=str)).dropna().astype(str).tolist()).lower()
    title_text = " ".join(inv.get("title", pd.Series(dtype=str)).dropna().astype(str).tolist()).lower()
    n_samples = int(len(audit))
    n_cov = int(audit["has_cov_file"].sum()) if not audit.empty else 0
    n_old_target = int((audit["is_old_104w_plus"] & audit["is_core_target_tissue"]).sum()) if not audit.empty else 0
    n_old_cortex = int((audit["is_old_104w_plus"] & audit["tissue"].astype(str).str.contains("brain_cortex")).sum()) if not audit.empty else 0
    n_old_heart = int((audit["is_old_104w_plus"] & audit["tissue"].astype(str).str.contains("heart")).sum()) if not audit.empty else 0
    n_old_lung = int((audit["is_old_104w_plus"] & audit["tissue"].astype(str).str.contains("lung")).sum()) if not audit.empty else 0
    context_blockers = []
    if "single-cell" in summary_text or "single cell" in summary_text or "single-cell" in title_text or "single cell" in title_text:
        context_blockers.append("series_single_cell_or_low_coverage_age_clock_context")
    if audit["low_coverage_or_tag_context"].mean() > 0.5 if not audit.empty else False:
        context_blockers.append("iPCRtag_or_iTAG_files_indicate_low_coverage_tagged_context")
    if n_old_cortex == 0:
        context_blockers.append("no_104w_plus_cortex_samples")
    if n_cov != n_samples:
        context_blockers.append("not_all_samples_have_cov_file")
    recommendation = "auxiliary_adapter_smoke_only"
    if not context_blockers and n_old_target >= 20 and n_old_cortex > 0:
        recommendation = "parser_smoke_candidate"
    payload = {
        "dataset": "GSE225166",
        "n_samples": n_samples,
        "n_cov_files_matched": n_cov,
        "n_old_core_target_samples": n_old_target,
        "n_old_cortex": n_old_cortex,
        "n_old_heart": n_old_heart,
        "n_old_lung": n_old_lung,
        "max_age_weeks": None if audit.empty or audit["age_weeks"].isna().all() else round(float(audit["age_weeks"].max()), 3),
        "raw_archive_size_bytes": None
        if supp.empty or supp[ "archive_or_file"].astype(str).ne("Archive").all()
        else int(float(supp[supp["archive_or_file"].astype(str).eq("Archive")].iloc[0].get("filelist_size_bytes") or 0)),
        "context_blockers": context_blockers,
        "recommendation": recommendation,
        "headline_benchmark_allowed": False,
    }
    return audit, summary, payload


def write_report(path: Path, ranking: pd.DataFrame, gse225166_summary: pd.DataFrame, audit_payload: dict, search_meta: dict) -> None:
    p1 = ranking[ranking.get("v8_5_direct_mainline_ready", pd.Series(False, index=ranking.index)).astype(bool)] if not ranking.empty else pd.DataFrame()
    candidate_cols = [
        "dataset",
        "v8_5_direct_mainline_ready",
        "v8_4_role",
        "v8_4_recommended_action",
        "target_old_counts",
        "max_age_by_target_tissue",
        "processed_schema_guess",
        "v8_4_blockers",
        "title",
        "geo_url",
    ]
    candidate_view = ranking[[col for col in candidate_cols if col in ranking.columns]].head(18) if not ranking.empty else pd.DataFrame()
    audit_view = gse225166_summary[
        (gse225166_summary["age_weeks"].ge(95, fill_value=False))
        | (gse225166_summary["tissue"].astype(str).str.contains("brain_cortex|heart|lung", na=False))
    ].copy() if not gse225166_summary.empty else pd.DataFrame()

    lines = [
        "# v8.5 Targeted GEO Refresh and GSE225166 Adapter Audit Report",
        "",
        "Date: 2026-05-18",
        "",
        "## Summary",
        "",
        "v8.5 refreshed official GEO candidates using targeted E-Utils queries and audited GSE225166 without downloading large supplements or running training.",
        "",
        f"- Targeted search candidates preflighted: {search_meta.get('candidate_count')}",
        f"- Direct mainline-ready new candidates: {len(p1)}",
        f"- GSE225166 recommendation: {audit_payload.get('recommendation')}",
        f"- GSE225166 old core-target samples: {audit_payload.get('n_old_core_target_samples')} "
        f"(cortex={audit_payload.get('n_old_cortex')}, heart={audit_payload.get('n_old_heart')}, lung={audit_payload.get('n_old_lung')})",
        "",
        "## Candidate Ranking",
        "",
        md_table(candidate_view),
        "",
        "## GSE225166 Adapter Audit",
        "",
        f"- Matched COV files: {audit_payload.get('n_cov_files_matched')}/{audit_payload.get('n_samples')}",
        f"- RAW archive size: {audit_payload.get('raw_archive_size_bytes')} bytes",
        f"- Context blockers: {', '.join(audit_payload.get('context_blockers', [])) or 'none'}",
        "- Headline benchmark allowed: False",
        "",
        md_table(audit_view, max_rows=30),
        "",
        "## Decision",
        "",
    ]
    if p1.empty:
        lines.extend(
            [
                "No new direct mainline-ready old brain_cortex/heart/lung bulk RRBS dataset was found. Do not download a new dataset for v8.5 and do not start autoresearch.",
                "",
                "GSE225166 has useful processed COV files and old heart/lung support, but its series context is single-cell/low-coverage age prediction and it lacks 104w+ cortex samples. It should remain an auxiliary adapter-smoke candidate, not a chronological headline benchmark dataset.",
            ]
        )
    else:
        lines.append("At least one direct candidate passed the metadata gate. Only those candidates should proceed to parser smoke; no training until parser/matrix manifests pass.")
    lines.extend(
        [
            "",
            "## Outputs",
            "",
            "- `results/v8_5_targeted_geo_refresh/targeted_search_inventory.csv`",
            "- `results/v8_5_targeted_geo_refresh/targeted_search_samples.csv`",
            "- `results/v8_5_targeted_geo_refresh/targeted_search_supplements.csv`",
            "- `results/v8_5_targeted_geo_refresh/v8_5_candidate_ranking.csv`",
            "- `results/v8_5_targeted_geo_refresh/gse225166_adapter_audit.csv`",
            "- `results/v8_5_targeted_geo_refresh/gse225166_tissue_age_summary.csv`",
            "- `results/v8_5_targeted_geo_refresh/gse225166_adapter_audit.json`",
            "",
            "## Official Sources",
            "",
            "- https://www.ncbi.nlm.nih.gov/geo/info/geo_paccess.html",
            "- https://www.ncbi.nlm.nih.gov/geo/info/download.html",
            "- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225166",
            "- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE232547",
            "- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171236",
            "- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138368",
            "- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151541",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out_dir", default=str(OUT_DIR))
    parser.add_argument("--retmax", type=int, default=25)
    parser.add_argument("--max_candidates", type=int, default=35)
    parser.add_argument("--skip_network", action="store_true")
    parser.add_argument("--report", default=str(REPORT))
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    discovery = discovery_module()
    v84 = v84_module()

    inventory_path = out_dir / "targeted_search_inventory.csv"
    samples_path = out_dir / "targeted_search_samples.csv"
    supplements_path = out_dir / "targeted_search_supplements.csv"

    if args.skip_network and inventory_path.exists() and samples_path.exists() and supplements_path.exists():
        inventory = pd.read_csv(inventory_path)
        samples = pd.read_csv(samples_path)
        supplements = pd.read_csv(supplements_path)
        accessions = sorted(inventory["dataset"].dropna().astype(str).unique())
    else:
        accessions, source_terms, summaries = discover_targeted_accessions(discovery, args.retmax)
        inventory, samples, supplements = preflight_candidates(
            discovery,
            accessions,
            source_terms,
            summaries,
            max_candidates=args.max_candidates,
        )
        inventory.to_csv(inventory_path, index=False)
        samples.to_csv(samples_path, index=False)
        supplements.to_csv(supplements_path, index=False)

    ranking = enrich_candidate_ranking(v84, inventory, samples)
    ranking.to_csv(out_dir / "v8_5_candidate_ranking.csv", index=False)

    audit, tissue_summary, audit_payload = audit_gse225166(v84, samples, supplements, inventory)
    audit.to_csv(out_dir / "gse225166_adapter_audit.csv", index=False)
    tissue_summary.to_csv(out_dir / "gse225166_tissue_age_summary.csv", index=False)
    write_json(out_dir / "gse225166_adapter_audit.json", audit_payload)

    search_meta = {
        "targeted_search_terms": TARGETED_SEARCH_TERMS,
        "seed_accessions": SEED_ACCESSIONS,
        "candidate_count": int(inventory["dataset"].nunique()) if not inventory.empty else 0,
        "direct_mainline_ready_count": int(ranking.get("v8_5_direct_mainline_ready", pd.Series(dtype=bool)).sum()) if not ranking.empty else 0,
    }
    write_json(out_dir / "v8_5_search_summary.json", search_meta)
    write_report(Path(args.report), ranking, tissue_summary, audit_payload, search_meta)
    print(
        json.dumps(
            {
                "status": "complete",
                "direct_mainline_ready_count": search_meta["direct_mainline_ready_count"],
                "gse225166_recommendation": audit_payload["recommendation"],
                "outputs": {
                    "ranking": str(out_dir / "v8_5_candidate_ranking.csv"),
                    "gse225166_audit": str(out_dir / "gse225166_adapter_audit.csv"),
                    "report": str(Path(args.report)),
                },
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()

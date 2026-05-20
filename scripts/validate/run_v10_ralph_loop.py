#!/usr/bin/env python3
"""v10 RALPH loop for broader data acquisition and minimal ETL decisions.

The controller is safe by default. It can refresh official GEO metadata and
write candidate gates, but it does not download large files, build matrices, or
train models. Adapter smoke downloads require --run_smoke. FASTQ/Bismark ETL is
not executed here; P3 candidates are written to a pilot manifest for a separate
minimal ETL plan.
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
OUT_DIR = ROOT / "results" / "ralph_v10_loop"
REPORT = OUT_DIR / "v10_success_or_failure_report.md"
REFERENCE_MATRIX = ROOT / "results" / "multidataset_v8_2_prefilter_liftover" / "all_rrbs_region_matrix_5kb.parquet"

TARGET_TISSUES = {"brain_cortex", "heart", "lung"}
ADJACENT_TISSUES = {
    "brain_other",
    "intestine",
    "colon",
    "liver",
    "spleen",
    "kidney",
    "blood",
    "skeletal_muscle",
    "adipose",
}
OLD_THRESHOLD_WEEKS = 104.0
COMMON_REGION_GATE = 50_000

V75_BASELINES = {
    "groupkfold_mae_weeks": 23.767,
    "gse121141_all_age_mae_weeks": 38.033,
    "gse121141_old104_mae_weeks": 75.386,
}

OFFICIAL_SOURCE_DOCS = {
    "geo_download": "https://www.ncbi.nlm.nih.gov/geo/info/download.html",
    "geo_programmatic_access": "https://www.ncbi.nlm.nih.gov/geo/info/geo_paccess.html",
    "geo_soft": "https://www.ncbi.nlm.nih.gov/geo/info/soft.html",
    "sra_download": "https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/",
    "ena_browser_api": "https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/browser-api.html",
}

V10_SEARCH_TERMS = [
    'Mus musculus[ORGN] AND "Methylation profiling by high throughput sequencing" AND (RRBS OR WGBS OR bisulfite OR "reduced representation") AND ("brain cortex" OR cortex OR heart OR lung) AND (aging OR aged OR old OR "24 month" OR "26 month" OR "28 month" OR "30 month") AND gse[ETYP]',
    'Mus musculus[ORGN] AND "DNA methylation" AND (RRBS OR WGBS OR bisulfite) AND ("brain" OR cortex OR heart OR lung OR colon OR liver) AND ("24 month" OR "26 month" OR "28 month" OR "30 month" OR aged OR old) AND gse[ETYP]',
    'Mus musculus[ORGN] AND "reduced representation bisulfite" AND (aging OR aged OR old OR lifespan OR rejuvenation) AND gse[ETYP]',
    'Mus musculus[ORGN] AND "whole genome bisulfite" AND (aging OR aged OR old) AND gse[ETYP]',
]

V10_SEED_ACCESSIONS = [
    "GSE233734",  # old colon RRBS, adjacent old-age support
    "GSE171604",  # old satellite cells WGBS, likely auxiliary
    "GSE171236",  # aged hippocampus/enrichment, adjacent brain-old
    "GSE156557",  # old reprogramming tissues, adjacent/intervention
    "GSE215310",
    "GSE221124",
    "GSE225166",
    "GSE83947",
    "GSE134398",
    "GSE134238",
    "GSE134397",
    "GSE232547",
    "GSE224442",
    "GSE92486",
    "GSE175410",
    "GSE295059",
    "GSE286302",
    "GSE281602",
]

INTEGRATED_DATASETS = {"GSE120137", "GSE80672", "GSE93957", "GSE121141", "GSE60012", "GSE213628"}
HARD_BLOCKERS = {
    "single_cell_or_targeted_celltype": ["single-cell", "single cell", "targeted_celltype_or_non_bulk_context"],
    "low_coverage_or_tagged": ["ipcrtag", "itag", "low coverage", "low_coverage_or_tagged_context", "padlock"],
    "untraceable_superseries": ["superseries_use_subseries", "superseries_or_subseries_required"],
    "organoid_or_in_vitro": ["organoid", "in vitro modeling"],
    "celltype_specific_hsc": [" hsc", "hematopoietic stem"],
}


def import_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def discovery_module():
    return import_module(ROOT / "scripts" / "etl" / "16_geo_old_age_candidate_discovery.py", "geo_old_age_candidate_discovery_v10")


def v9_module():
    return import_module(ROOT / "scripts" / "validate" / "run_v9_ralph_loop.py", "run_v9_ralph_loop_helpers_v10")


def read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path) if path.exists() else pd.DataFrame()


def combine_frames(paths: list[Path]) -> pd.DataFrame:
    frames = []
    for path in paths:
        df = read_csv(path)
        if not df.empty:
            df = df.copy()
            df["v10_source_table"] = str(path.relative_to(ROOT))
            frames.append(df)
    return pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()


def load_cached_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    inventory = combine_frames(
        [
            ROOT / "metadata" / "geo_old_age_candidate_inventory.csv",
            ROOT / "results" / "v8_5_targeted_geo_refresh" / "targeted_search_inventory.csv",
            ROOT / "results" / "ralph_v9_loop" / "v9_research_inventory.csv",
            OUT_DIR / "v10_research_inventory.csv",
        ]
    )
    samples = combine_frames(
        [
            ROOT / "metadata" / "geo_old_age_candidate_samples.csv",
            ROOT / "results" / "v8_5_targeted_geo_refresh" / "targeted_search_samples.csv",
            ROOT / "results" / "ralph_v9_loop" / "v9_research_samples.csv",
            OUT_DIR / "v10_research_samples.csv",
        ]
    )
    supplements = combine_frames(
        [
            ROOT / "metadata" / "geo_old_age_candidate_supplements.csv",
            ROOT / "results" / "v8_5_targeted_geo_refresh" / "targeted_search_supplements.csv",
            ROOT / "results" / "ralph_v9_loop" / "v9_research_supplements.csv",
            OUT_DIR / "v10_research_supplements.csv",
        ]
    )
    if not inventory.empty and "dataset" in inventory:
        inventory = inventory.drop_duplicates("dataset", keep="last").reset_index(drop=True)
    if not samples.empty and {"dataset", "sample_id"}.issubset(samples.columns):
        samples = samples.drop_duplicates(["dataset", "sample_id"], keep="last").reset_index(drop=True)
    if not supplements.empty and {"dataset", "supplement_name"}.issubset(supplements.columns):
        supplements = supplements.drop_duplicates(["dataset", "supplement_name"], keep="last").reset_index(drop=True)
    return inventory, samples, supplements


def refresh_candidates(discovery, out_dir: Path, retmax: int, max_candidates: int) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    accessions = set(V10_SEED_ACCESSIONS)
    source_terms = {accession: ["manual_seed_v10"] for accession in accessions}
    summaries: dict[str, dict] = {}
    found, found_terms, found_summaries = discovery.discover_accessions(V10_SEARCH_TERMS, retmax=retmax)
    accessions.update(found)
    summaries.update(found_summaries)
    for accession, terms in found_terms.items():
        source_terms.setdefault(accession, []).extend(terms)
    for accession in sorted(accessions):
        if accession not in summaries:
            try:
                summaries[accession] = discovery.eutils_summary_for_accession(accession)
            except Exception as exc:
                summaries[accession] = {"dataset": accession, "eutils_error": str(exc)[:300]}
    ranked = sorted(
        accessions,
        key=lambda value: (
            -discovery.rough_esummary_score(value, summaries.get(value, {}), source_terms.get(value, [])),
            int(value.removeprefix("GSE")) if value.removeprefix("GSE").isdigit() else 999_999_999,
        ),
    )
    selected = list(dict.fromkeys(ranked[:max_candidates] + V10_SEED_ACCESSIONS))
    candidate_rows: list[dict] = []
    sample_rows: list[dict] = []
    supplement_rows: list[dict] = []
    for accession in selected:
        try:
            row, supplements, samples = discovery.build_candidate(
                accession,
                source_terms=source_terms.get(accession, ["v10_search"]),
                eutils_summary=summaries.get(accession),
                refresh_soft=False,
                do_head=True,
            )
        except Exception as exc:
            row = {
                "dataset": accession,
                "geo_url": discovery.geo_url(accession),
                "priority_tier": "ERROR",
                "recommended_action": "inspect_manually",
                "blockers": str(exc)[:500],
                "source_queries": " || ".join(source_terms.get(accession, [])),
            }
            supplements = []
            samples = []
        row["v10_source_table"] = "v10_refresh"
        candidate_rows.append(row)
        for sample in samples:
            sample["v10_source_table"] = "v10_refresh"
        for supplement in supplements:
            supplement["v10_source_table"] = "v10_refresh"
        sample_rows.extend(samples)
        supplement_rows.extend(supplements)
        time.sleep(0.05)
    inventory = pd.DataFrame(candidate_rows)
    samples = pd.DataFrame(sample_rows)
    supplements = pd.DataFrame(supplement_rows)
    inventory.to_csv(out_dir / "v10_research_inventory.csv", index=False)
    samples.to_csv(out_dir / "v10_research_samples.csv", index=False)
    supplements.to_csv(out_dir / "v10_research_supplements.csv", index=False)
    return inventory, samples, supplements


def candidate_text(inventory_row: pd.Series | None, supplement_rows: pd.DataFrame) -> str:
    text = ""
    if inventory_row is not None:
        text += " ".join(
            str(inventory_row.get(field, ""))
            for field in ["title", "summary", "overall_design", "blockers", "processed_schema_guess", "tissues"]
        )
    if not supplement_rows.empty:
        text += " " + " ".join(supplement_rows.get("supplement_name", pd.Series(dtype=str)).fillna("").astype(str).tolist())
    return text.lower()


def hard_blockers(text: str) -> list[str]:
    blockers = []
    for label, needles in HARD_BLOCKERS.items():
        if any(needle in text for needle in needles):
            blockers.append(label)
    return list(dict.fromkeys(blockers))


def has_processed_methylation(matches: pd.DataFrame, inventory_row: pd.Series | None) -> bool:
    if not matches.empty:
        return True
    schema = "" if inventory_row is None else str(inventory_row.get("processed_schema_guess", "")).lower()
    return any(token in schema for token in ["cov", "bedgraph", "methylation_raw_tar"])


def raw_pilot_possible(inventory_row: pd.Series | None, supplement_rows: pd.DataFrame) -> bool:
    text = candidate_text(inventory_row, supplement_rows)
    if "sra" in text or "fastq" in text or "raw.tar" in text:
        return True
    preferred = "" if inventory_row is None else str(inventory_row.get("supp_file_eutils", "")).lower()
    return "sra" in preferred or "fastq" in preferred


def classify(
    old_target_tissues: set[str],
    old_adjacent_tissues: set[str],
    old_target_samples: int,
    old_adjacent_samples: int,
    processed: bool,
    raw_possible: bool,
    blockers: list[str],
    schema_pass: bool,
    assembly_pass: bool,
) -> tuple[str, str, bool]:
    if old_target_samples > 0 and processed and not blockers:
        if "brain_cortex" in old_target_tissues or len(old_target_tissues) >= 2:
            return "P1_headline_exact_target", "matrix_gate_required" if schema_pass and assembly_pass else "adapter_or_assembly_gate_pending", True
        return "P2_adjacent_old_age_auxiliary", "single_target_tissue_only", False
    if old_adjacent_samples > 0 and processed and not blockers:
        return "P2_adjacent_old_age_auxiliary", "old_non_target_tissue_only", False
    if old_target_samples > 0 and raw_possible and not blockers:
        return "P3_raw_exact_target_pilot", "plan_2_3_sample_fastq_bismark_pilot", False
    if blockers:
        return "P4_low_coverage_or_targeted_bridge", ";".join(blockers), False
    return "P5_blocked", "no_sample_specific_old_target_or_adjacent_tissue", False


def score_inventory(v9, row: pd.Series) -> float:
    score = v9.candidate_score(row)
    text = " ".join(str(row.get(field, "")) for field in ["title", "summary", "tissues", "blockers"]).lower()
    if "colon" in text or "intestine" in text:
        score += 8
    if "hippocampus" in text or "brain" in text:
        score += 8
    if "24" in text or "28" in text or "30" in text:
        score += 4
    return score


def candidate_order(v9, inventory: pd.DataFrame, samples: pd.DataFrame) -> list[str]:
    scored: list[tuple[float, str]] = []
    if not inventory.empty and "dataset" in inventory:
        for _, row in inventory.iterrows():
            dataset = str(row.get("dataset", ""))
            if dataset.startswith("GSE") and dataset not in INTEGRATED_DATASETS:
                scored.append((score_inventory(v9, row), dataset))
    if not samples.empty and "dataset" in samples:
        for dataset in samples["dataset"].dropna().astype(str).tolist():
            if dataset.startswith("GSE") and dataset not in INTEGRATED_DATASETS:
                scored.append((0.0, dataset))
    for dataset in V10_SEED_ACCESSIONS:
        if dataset not in INTEGRATED_DATASETS:
            scored.append((1.0, dataset))
    ordered = []
    for _, dataset in sorted(scored, key=lambda item: (-item[0], item[1])):
        if dataset not in ordered:
            ordered.append(dataset)
    return ordered


def summarize_candidate(dataset: str, inventory: pd.DataFrame, samples: pd.DataFrame, supplements: pd.DataFrame, v9, v84, v8, out_dir: Path, run_smoke: bool, network_log: Path) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    inv_row = None
    if not inventory.empty and "dataset" in inventory and dataset in set(inventory["dataset"].astype(str)):
        inv_row = inventory[inventory["dataset"].astype(str).eq(dataset)].iloc[0]
    ds_samples = samples[samples["dataset"].astype(str).eq(dataset)].copy() if not samples.empty else pd.DataFrame()
    ds_supplements = supplements[supplements["dataset"].astype(str).eq(dataset)].copy() if not supplements.empty else pd.DataFrame()
    support = v9.sample_support(ds_samples, v84)
    old = support[support["is_old104"].astype(bool)] if not support.empty else pd.DataFrame()
    old_target = old[old["tissue"].astype(str).isin(TARGET_TISSUES)] if not old.empty else pd.DataFrame()
    old_adjacent = old[old["tissue"].astype(str).isin(ADJACENT_TISSUES)] if not old.empty else pd.DataFrame()
    matched = v9.match_supplement_rows(ds_samples, ds_supplements, v84)
    selected = v9.select_smoke_files(matched)
    reference_regions = v9.load_reference_regions(v8)
    smoke = v9.run_or_plan_smoke(dataset, selected, out_dir, v8, reference_regions, run_smoke, network_log)
    matched.to_csv(out_dir / f"{dataset}_matched_supplement_files.csv", index=False)
    smoke.to_csv(out_dir / f"{dataset}_adapter_smoke.csv", index=False)
    text = candidate_text(inv_row, ds_supplements)
    blockers = hard_blockers(text)
    schema_pass = v9.schema_smoke_pass(smoke)
    assembly_pass = v9.assembly_compatible(inv_row, smoke)
    processed = has_processed_methylation(matched, inv_row)
    raw_possible = raw_pilot_possible(inv_row, ds_supplements)
    tier, reason, headline_candidate = classify(
        set(old_target["tissue"].astype(str)) if not old_target.empty else set(),
        set(old_adjacent["tissue"].astype(str)) if not old_adjacent.empty else set(),
        int(len(old_target)),
        int(len(old_adjacent)),
        processed,
        raw_possible,
        blockers,
        schema_pass,
        assembly_pass,
    )
    next_step = {
        "P1_headline_exact_target": "run_adapter_smoke_or_build_full_matrix_gate",
        "P2_adjacent_old_age_auxiliary": "keep_for_adjacent_or_single_tissue_diagnostics",
        "P3_raw_exact_target_pilot": "write_minimal_fastq_bismark_pilot_plan",
        "P4_low_coverage_or_targeted_bridge": "compatibility_audit_only",
        "P5_blocked": "do_not_promote",
    }[tier]
    row = {
        "iteration": "v10.0",
        "dataset": dataset,
        "candidate_tier": tier,
        "hypothesis": "broader acquisition or minimal ETL may repair old target tissue support",
        "official_metadata_pass": bool(len(old_target) > 0 or len(old_adjacent) > 0),
        "n_samples": int(len(ds_samples)),
        "n_old_target_samples": int(len(old_target)),
        "old_target_tissues": ";".join(sorted(set(old_target["tissue"].astype(str)))) if not old_target.empty else "",
        "n_old_adjacent_samples": int(len(old_adjacent)),
        "old_adjacent_tissues": ";".join(sorted(set(old_adjacent["tissue"].astype(str)))) if not old_adjacent.empty else "",
        "n_matched_processed_files": int(len(matched)),
        "n_selected_smoke_files": int(len(selected)),
        "schema_smoke_pass": schema_pass,
        "assembly_compatible": assembly_pass,
        "smoke_common_regions_estimate": v9.smoke_common_regions(smoke),
        "common_regions_required": COMMON_REGION_GATE,
        "headline_candidate": bool(headline_candidate),
        "minimal_fastq_pilot_candidate": tier == "P3_raw_exact_target_pilot",
        "raw_pilot_possible": raw_possible,
        "auxiliary_or_blocker_reason": reason,
        "context_blockers": ";".join(blockers),
        "next_step": next_step,
    }
    return row, matched, smoke


def decide_state(gate: pd.DataFrame, refresh_round: int) -> dict:
    p1 = gate[gate["candidate_tier"].astype(str).eq("P1_headline_exact_target")] if not gate.empty else pd.DataFrame()
    p3 = gate[gate["candidate_tier"].astype(str).eq("P3_raw_exact_target_pilot")] if not gate.empty else pd.DataFrame()
    if not p1.empty:
        decision = "continue_v10_p1_matrix_or_smoke_gate"
        status = "continue"
        next_action = "Run adapter smoke for P1 if needed, then build matrix only if common regions >=50000."
    elif not p3.empty:
        decision = "continue_v10_minimal_etl_pilot"
        status = "continue"
        next_action = "Prepare a 2-3 sample FASTQ/Bismark pilot; do not run full ETL."
    elif refresh_round >= 2:
        decision = "v10_failure_no_p1_or_p3_after_broadened_search"
        status = "failure"
        next_action = "Stop data-candidate loop; plan broader acquisition outside current public GEO processed/seed space."
    else:
        decision = "continue_v10_research"
        status = "continue"
        next_action = "Run a second broadened official refresh before declaring v10 failure."
    return {
        "timestamp": v9_module().utc_now(),
        "loop_version": "v10.0",
        "status": status,
        "refresh_round": refresh_round,
        "loop_decision": decision,
        "next_action": next_action,
        "success_criteria": {
            "dataset_gate": "P1 or P3 full dataset with common 5kb regions >=50000",
            "gse121141_old104_mae_improvement_weeks": ">=10 vs 75.386w",
            "gse121141_all_age_mae_max": "<=40.033w",
            "groupkfold_mae_max": "<=25.767w",
            "random_label_sanity": "abs(r)<0.2 and MAE random-like",
            "shuffled_cr_sanity": "AUC chance-like and Cohen_d decreased",
        },
        "failure_criteria": {
            "two_broadened_refresh_rounds_without_p1_or_p3": bool(refresh_round >= 2 and p1.empty and p3.empty),
            "only_auxiliary_or_blocked_data": bool(not gate.empty and p1.empty and p3.empty),
            "minimal_etl_pilot_failure_blocks_full_etl": True,
        },
        "metrics": {
            "n_candidates": int(len(gate)),
            "n_p1": int((gate["candidate_tier"].astype(str).eq("P1_headline_exact_target")).sum()) if not gate.empty else 0,
            "n_p2": int((gate["candidate_tier"].astype(str).eq("P2_adjacent_old_age_auxiliary")).sum()) if not gate.empty else 0,
            "n_p3": int((gate["candidate_tier"].astype(str).eq("P3_raw_exact_target_pilot")).sum()) if not gate.empty else 0,
            "n_p4": int((gate["candidate_tier"].astype(str).eq("P4_low_coverage_or_targeted_bridge")).sum()) if not gate.empty else 0,
            "n_p5": int((gate["candidate_tier"].astype(str).eq("P5_blocked")).sum()) if not gate.empty else 0,
            "v75_gse121141_old104_mae_weeks": V75_BASELINES["gse121141_old104_mae_weeks"],
        },
    }


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def write_report(path: Path, gate: pd.DataFrame, decision: dict, refreshed: bool, run_smoke: bool) -> None:
    v9 = v9_module()
    cols = [
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
    view = gate[[col for col in cols if col in gate.columns]].copy() if not gate.empty else pd.DataFrame()
    lines = [
        "# v10 RALPH Loop Success or Failure Report",
        "",
        f"Date: {v9.utc_now()}",
        "",
        "## Summary",
        "",
        "v10 broadened the data strategy beyond current v9 candidates. It did not download large files, build matrices, run FASTQ/Bismark, train models, or start autoresearch.",
        "",
        f"- Loop decision: `{decision['loop_decision']}`",
        f"- Status: `{decision['status']}`",
        f"- Next action: {decision['next_action']}",
        f"- Candidate refresh used: `{refreshed}`",
        f"- Adapter smoke download used: `{run_smoke}`",
        f"- P1 exact target candidates: {decision['metrics']['n_p1']}",
        f"- P3 raw exact target pilot candidates: {decision['metrics']['n_p3']}",
        "",
        "## Candidate Gate Table",
        "",
        v9.md_table(view, max_rows=40),
        "",
        "## Official Network/API Rules",
        "",
        "Use E-Utils/SOFT for metadata, GEO FTP for supplements, and SRA/ENA only for raw pilot decisions. Non-official web search can add seeds, but every seed must be verified through official records.",
        "",
        v9.md_table(pd.DataFrame([{"source": key, "url": value} for key, value in OFFICIAL_SOURCE_DOCS.items()])),
        "",
        "## Outputs",
        "",
        "- `results/ralph_v10_loop/ralph_iteration_log.jsonl`",
        "- `results/ralph_v10_loop/ralph_decision_state.json`",
        "- `results/ralph_v10_loop/candidate_gate_table.csv`",
        "- `results/ralph_v10_loop/network_resolution_log.jsonl`",
        "- `results/ralph_v10_loop/candidate_smoke_manifest.csv`",
        "- `results/ralph_v10_loop/minimal_etl_pilot_manifest.csv`",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out_dir", default=str(OUT_DIR))
    parser.add_argument("--report", default=str(REPORT))
    parser.add_argument("--refresh_candidates", action="store_true")
    parser.add_argument("--run_smoke", action="store_true")
    parser.add_argument("--retmax", type=int, default=50)
    parser.add_argument("--max_candidates", type=int, default=45)
    parser.add_argument("--refresh_round", type=int, default=1)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    v9 = v9_module()
    discovery = discovery_module()
    v84 = v9.v84_module()
    v8 = v9.v8_loop_module()
    network_log = out_dir / "network_resolution_log.jsonl"
    v9.append_jsonl(
        network_log,
        {
            "timestamp": v9.utc_now(),
            "blocker_type": "official_database_rules_loaded",
            "attempted_url": "",
            "source_doc": " | ".join(OFFICIAL_SOURCE_DOCS.values()),
            "resolution": "v10_broader_search_uses_official_metadata_before_download",
            "applied_rule": "no_large_download_no_training_until_p1_or_p3_gate",
        },
    )

    cached_inventory, cached_samples, cached_supplements = load_cached_inputs()
    frames_inventory = [cached_inventory]
    frames_samples = [cached_samples]
    frames_supplements = [cached_supplements]
    if args.refresh_candidates:
        refreshed = refresh_candidates(discovery, out_dir, args.retmax, args.max_candidates)
        frames_inventory.append(refreshed[0])
        frames_samples.append(refreshed[1])
        frames_supplements.append(refreshed[2])

    inventory = pd.concat([df for df in frames_inventory if not df.empty], ignore_index=True, sort=False) if any(not df.empty for df in frames_inventory) else pd.DataFrame()
    samples = pd.concat([df for df in frames_samples if not df.empty], ignore_index=True, sort=False) if any(not df.empty for df in frames_samples) else pd.DataFrame()
    supplements = pd.concat([df for df in frames_supplements if not df.empty], ignore_index=True, sort=False) if any(not df.empty for df in frames_supplements) else pd.DataFrame()
    if not inventory.empty and "dataset" in inventory:
        inventory = inventory.drop_duplicates("dataset", keep="last").reset_index(drop=True)
    if not samples.empty and {"dataset", "sample_id"}.issubset(samples.columns):
        samples = samples.drop_duplicates(["dataset", "sample_id"], keep="last").reset_index(drop=True)
    if not supplements.empty and {"dataset", "supplement_name"}.issubset(supplements.columns):
        supplements = supplements.drop_duplicates(["dataset", "supplement_name"], keep="last").reset_index(drop=True)

    gate_rows = []
    smoke_rows = []
    pilot_rows = []
    iteration_log = out_dir / "ralph_iteration_log.jsonl"
    for dataset in candidate_order(v9, inventory, samples)[: args.max_candidates]:
        row, matched, smoke = summarize_candidate(dataset, inventory, samples, supplements, v9, v84, v8, out_dir, args.run_smoke, network_log)
        gate_rows.append(row)
        for _, smoke_row in smoke.iterrows():
            smoke_rows.append(
                {
                    "dataset": dataset,
                    "sample_id": smoke_row.get("sample_id", ""),
                    "supplement_name": smoke_row.get("supplement_name", ""),
                    "status": smoke_row.get("status", ""),
                    "schema_guess": smoke_row.get("schema_guess", ""),
                    "n_smoke_common_regions_with_reference": smoke_row.get("n_smoke_common_regions_with_reference", ""),
                }
            )
        if row["minimal_fastq_pilot_candidate"]:
            pilot_rows.append(
                {
                    "dataset": dataset,
                    "n_old_target_samples": row["n_old_target_samples"],
                    "old_target_tissues": row["old_target_tissues"],
                    "pilot_sample_count": min(3, row["n_old_target_samples"]),
                    "pilot_rule": "select oldest target-tissue samples with sample-specific age and raw availability",
                    "status": "planned_not_executed",
                }
            )
        v9.append_jsonl(
            iteration_log,
            {
                "timestamp": v9.utc_now(),
                "iteration": row["iteration"],
                "dataset": dataset,
                "candidate_tier": row["candidate_tier"],
                "metrics": {
                    "n_old_target_samples": row["n_old_target_samples"],
                    "n_old_adjacent_samples": row["n_old_adjacent_samples"],
                    "n_matched_processed_files": row["n_matched_processed_files"],
                    "schema_smoke_pass": row["schema_smoke_pass"],
                },
                "gate_decision": row["next_step"],
            },
        )

    gate = pd.DataFrame(gate_rows)
    smoke_manifest = pd.DataFrame(smoke_rows)
    pilot_manifest = pd.DataFrame(pilot_rows)
    gate.to_csv(out_dir / "candidate_gate_table.csv", index=False)
    smoke_manifest.to_csv(out_dir / "candidate_smoke_manifest.csv", index=False)
    pilot_manifest.to_csv(out_dir / "minimal_etl_pilot_manifest.csv", index=False)
    decision = decide_state(gate, args.refresh_round)
    write_json(out_dir / "ralph_decision_state.json", decision)
    write_report(Path(args.report), gate, decision, args.refresh_candidates, args.run_smoke)
    print(json.dumps(decision, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

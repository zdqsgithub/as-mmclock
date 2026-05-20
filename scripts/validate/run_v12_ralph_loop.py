#!/usr/bin/env python3
"""v12 RALPH loop for data strategy and benchmark redefinition.

This controller is intentionally safe by default. It reuses v9-v11 official
refreshes, raw RFC evidence, and v11.x conversion gates to close the three-round
v12 data-strategy loop. It does not download FASTQ, build new matrices, run
benchmarks, or start autoresearch.
"""
from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")
OUT_DIR = ROOT / "results" / "ralph_v12_loop"
REPORT = OUT_DIR / "v12_success_or_failure_report.md"
REDEFINITION_RFC = OUT_DIR / "benchmark_redefinition_rfc.md"

TARGET_TISSUES = {"brain_cortex", "cortex", "heart", "lung"}
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
    "ncbi_eutils": "https://www.ncbi.nlm.nih.gov/books/NBK25501/",
}

STANDARD_COLUMNS = [
    "candidate_id",
    "refresh_round",
    "dataset",
    "accession_type",
    "source_query",
    "official_url",
    "title",
    "sample_id",
    "biosample",
    "run_accession",
    "age_weeks",
    "raw_age_token",
    "tissue",
    "tissue_source",
    "sex",
    "strain",
    "assay",
    "processed_schema_guess",
    "supplement_url",
    "run_size_bytes",
    "tier",
    "gate_status",
    "blocker_type",
    "recommended_action",
    "sample_specific_age_pass",
    "target_tissue_pass",
    "bulk_context_pass",
    "schema_smoke_pass",
    "assembly_traceable",
    "common_regions_estimate",
    "headline_allowed",
    "auxiliary_only_reason",
]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def clean(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    return str(value)


def boolish(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    text = clean(value).strip().lower()
    return text in {"1", "true", "yes", "y", "pass", "passed"}


def number_or_blank(value: Any) -> Any:
    text = clean(value)
    if not text:
        return ""
    try:
        value_f = float(text)
    except ValueError:
        return text
    if math.isnan(value_f):
        return ""
    if value_f.is_integer():
        return int(value_f)
    return value_f


def positive_number(value: Any) -> bool:
    parsed = number_or_blank(value)
    if parsed == "":
        return False
    try:
        return float(parsed) > 0
    except (TypeError, ValueError):
        return False


def read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path) if path.exists() else pd.DataFrame()


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


def combine_inventory() -> pd.DataFrame:
    paths = [
        ROOT / "metadata" / "geo_old_age_candidate_inventory.csv",
        ROOT / "results" / "ralph_v9_loop" / "v9_research_inventory.csv",
        ROOT / "results" / "ralph_v10_loop" / "v10_research_inventory.csv",
        ROOT / "results" / "ralph_v11_data_strategy" / "verified_new_gse_inventory.csv",
    ]
    frames = []
    for path in paths:
        df = read_csv(path)
        if not df.empty:
            df = df.copy()
            df["inventory_source"] = str(path.relative_to(ROOT))
            frames.append(df)
    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, ignore_index=True, sort=False)
    if "dataset" in out.columns:
        out = out.drop_duplicates("dataset", keep="last")
    return out.reset_index(drop=True)


def inventory_lookup(inventory: pd.DataFrame, dataset: str, field: str) -> str:
    if inventory.empty or "dataset" not in inventory.columns or field not in inventory.columns:
        return ""
    hit = inventory[inventory["dataset"].astype(str).eq(dataset)]
    if hit.empty:
        return ""
    return clean(hit.iloc[0].get(field, ""))


def map_v10_tier(value: Any) -> str:
    text = clean(value)
    if text.startswith("P1"):
        return "P1_processed_headline"
    if text.startswith("P2"):
        return "P2_auxiliary"
    if text.startswith("P3"):
        return "P3_minimal_etl_rfc"
    if text.startswith("P4"):
        return "P4_blocked_context"
    return "P5_redefinition_evidence"


def normalize_v10_like(path: Path, refresh_round: str, inventory: pd.DataFrame) -> list[dict[str, Any]]:
    df = read_csv(path)
    rows: list[dict[str, Any]] = []
    if df.empty:
        return rows
    for _, row in df.iterrows():
        dataset = clean(row.get("dataset"))
        if not dataset:
            continue
        old_target = clean(row.get("old_target_tissues"))
        old_adjacent = clean(row.get("old_adjacent_tissues"))
        tier = map_v10_tier(row.get("candidate_tier"))
        common = number_or_blank(row.get("smoke_common_regions_estimate"))
        context_blockers = clean(row.get("context_blockers"))
        reason = clean(row.get("auxiliary_or_blocker_reason"))
        rows.append(
            {
                "candidate_id": f"{refresh_round}:{dataset}",
                "refresh_round": refresh_round,
                "dataset": dataset,
                "accession_type": "GEO_series",
                "source_query": inventory_lookup(inventory, dataset, "source_queries") or clean(row.get("hypothesis")),
                "official_url": inventory_lookup(inventory, dataset, "geo_url") or f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={dataset}",
                "title": inventory_lookup(inventory, dataset, "title"),
                "sample_id": "",
                "biosample": "",
                "run_accession": "",
                "age_weeks": inventory_lookup(inventory, dataset, "age_values_weeks"),
                "raw_age_token": "",
                "tissue": old_target or old_adjacent or inventory_lookup(inventory, dataset, "tissues"),
                "tissue_source": "GEO_SOFT_candidate_samples",
                "sex": "",
                "strain": "",
                "assay": inventory_lookup(inventory, dataset, "experiment_type"),
                "processed_schema_guess": inventory_lookup(inventory, dataset, "processed_schema_guess"),
                "supplement_url": inventory_lookup(inventory, dataset, "preferred_url"),
                "run_size_bytes": inventory_lookup(inventory, dataset, "preferred_size_bytes"),
                "tier": tier,
                "gate_status": clean(row.get("next_step")),
                "blocker_type": context_blockers or reason,
                "recommended_action": clean(row.get("next_step")),
                "sample_specific_age_pass": boolish(row.get("official_metadata_pass")),
                "target_tissue_pass": positive_number(row.get("n_old_target_samples")),
                "bulk_context_pass": not bool(context_blockers),
                "schema_smoke_pass": boolish(row.get("schema_smoke_pass")),
                "assembly_traceable": boolish(row.get("assembly_compatible")),
                "common_regions_estimate": common,
                "headline_allowed": boolish(row.get("headline_candidate")),
                "auxiliary_only_reason": reason or context_blockers,
            }
        )
    return rows


def normalize_v11x(path: Path, inventory: pd.DataFrame) -> list[dict[str, Any]]:
    df = read_csv(path)
    rows: list[dict[str, Any]] = []
    if df.empty:
        return rows
    for _, row in df.iterrows():
        dataset = clean(row.get("dataset"))
        if not dataset:
            continue
        reason = clean(row.get("reason"))
        exact_age = boolish(row.get("exact_age_gate"))
        target = boolish(row.get("old_target_tissue_gate"))
        bulk = boolish(row.get("bulk_context_gate"))
        common = number_or_blank(row.get("common_regions_with_v8_2"))
        if clean(row.get("final_status")) == "auxiliary_matrix":
            tier = "P2_auxiliary_matrix"
        elif "adapter_required" in clean(row.get("decision")):
            tier = "P2_adapter_required_auxiliary"
        else:
            tier = "P5_redefinition_evidence"
        rows.append(
            {
                "candidate_id": f"round3_downloaded_backlog:{dataset}",
                "refresh_round": "round3_downloaded_backlog",
                "dataset": dataset,
                "accession_type": "GEO_processed_backlog",
                "source_query": "v11.x downloaded processed-data backlog",
                "official_url": inventory_lookup(inventory, dataset, "geo_url") or f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={dataset}",
                "title": inventory_lookup(inventory, dataset, "title"),
                "sample_id": "",
                "biosample": "",
                "run_accession": "",
                "age_weeks": "",
                "raw_age_token": "",
                "tissue": inventory_lookup(inventory, dataset, "tissues"),
                "tissue_source": "v11_x_conversion_gate",
                "sex": "",
                "strain": "",
                "assay": inventory_lookup(inventory, dataset, "experiment_type"),
                "processed_schema_guess": inventory_lookup(inventory, dataset, "processed_schema_guess"),
                "supplement_url": inventory_lookup(inventory, dataset, "preferred_url"),
                "run_size_bytes": "",
                "tier": tier,
                "gate_status": clean(row.get("decision")),
                "blocker_type": reason,
                "recommended_action": "do_not_train_keep_auxiliary",
                "sample_specific_age_pass": exact_age,
                "target_tissue_pass": target,
                "bulk_context_pass": bulk,
                "schema_smoke_pass": clean(row.get("final_status")) == "auxiliary_matrix",
                "assembly_traceable": boolish(row.get("common_region_gate")),
                "common_regions_estimate": common,
                "headline_allowed": boolish(row.get("headline_allowed")),
                "auxiliary_only_reason": reason,
            }
        )
    return rows


def normalize_raw_biosample(summary_path: Path, manifest_path: Path) -> list[dict[str, Any]]:
    summary = read_csv(summary_path)
    manifest = read_csv(manifest_path)
    rows: list[dict[str, Any]] = []
    if summary.empty:
        return rows
    for _, row in summary.iterrows():
        candidate = clean(row.get("candidate"))
        if not candidate:
            continue
        hits = manifest[manifest["candidate"].astype(str).eq(candidate)] if not manifest.empty and "candidate" in manifest.columns else pd.DataFrame()
        include = hits[hits.get("include_in_minimal_pilot", pd.Series(dtype=bool)).map(boolish)] if not hits.empty and "include_in_minimal_pilot" in hits.columns else pd.DataFrame()
        blockers = sorted(set(clean(value) for value in hits.get("blockers", pd.Series(dtype=str)).dropna() if clean(value))) if not hits.empty else []
        exact_n = int(number_or_blank(row.get("exact_age_n")) or 0)
        n_samples = int(number_or_blank(row.get("n_samples")) or 0)
        target_pass = bool(clean(row.get("tissues")) in TARGET_TISSUES or any(t in clean(row.get("tissues")) for t in TARGET_TISSUES))
        exact_pass = exact_n >= max(1, int(0.95 * n_samples)) if n_samples else False
        if not include.empty and target_pass and exact_pass and not blockers:
            tier = "P3_minimal_etl_rfc"
            recommended = "prepare_2_3_sample_fastq_bismark_pilot_rfc"
        elif not include.empty:
            tier = "P3_diagnostic_raw_rfc"
            recommended = "diagnostic_only_do_not_promote_without_exact_age_bulk_context"
        else:
            tier = "P5_redefinition_evidence"
            recommended = "do_not_download_fastq"
        run_accessions = ";".join(sorted(set(clean(v) for v in hits.get("run", pd.Series(dtype=str)).dropna() if clean(v)))) if not hits.empty else ""
        biosamples = ";".join(sorted(set(clean(v) for v in hits.get("biosample", pd.Series(dtype=str)).dropna() if clean(v)))) if not hits.empty else ""
        rows.append(
            {
                "candidate_id": f"round3_raw_rfc:{candidate}",
                "refresh_round": "round3_raw_rfc",
                "dataset": candidate.split("_")[0],
                "accession_type": "SRA_BioSample_raw_rfc",
                "source_query": "v11 raw BioSample/RunInfo verification",
                "official_url": "",
                "title": candidate,
                "sample_id": "",
                "biosample": biosamples,
                "run_accession": run_accessions,
                "age_weeks": "",
                "raw_age_token": clean(row.get("age_groups")),
                "tissue": clean(row.get("tissues")),
                "tissue_source": "BioSample_RunInfo",
                "sex": "",
                "strain": "",
                "assay": "raw_fastq_bisulfite_candidate",
                "processed_schema_guess": "raw_fastq_requires_bismark",
                "supplement_url": "",
                "run_size_bytes": int(float(row.get("total_size_mb", 0)) * 1024 * 1024) if clean(row.get("total_size_mb")) else "",
                "tier": tier,
                "gate_status": recommended,
                "blocker_type": ";".join(blockers),
                "recommended_action": recommended,
                "sample_specific_age_pass": exact_pass,
                "target_tissue_pass": target_pass,
                "bulk_context_pass": not blockers,
                "schema_smoke_pass": False,
                "assembly_traceable": False,
                "common_regions_estimate": "",
                "headline_allowed": False,
                "auxiliary_only_reason": ";".join(blockers) or clean(row.get("pilot_role")),
            }
        )
    return rows


def aggregate_smoke_manifests() -> pd.DataFrame:
    frames = []
    for path in [
        ROOT / "results" / "ralph_v9_loop" / "candidate_smoke_manifest.csv",
        ROOT / "results" / "ralph_v10_loop" / "candidate_smoke_manifest.csv",
        ROOT / "results" / "ralph_v11_data_strategy" / "candidate_smoke_manifest.csv",
    ]:
        df = read_csv(path)
        if not df.empty:
            df = df.copy()
            df["source_manifest"] = str(path.relative_to(ROOT))
            frames.append(df)
    return pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()


def build_minimal_etl_manifest(candidate_table: pd.DataFrame) -> pd.DataFrame:
    raw = candidate_table[candidate_table["accession_type"].astype(str).eq("SRA_BioSample_raw_rfc")].copy()
    if raw.empty:
        return pd.DataFrame(
            columns=[
                "candidate_id",
                "dataset",
                "pilot_status",
                "pilot_reason",
                "pilot_sample_count",
                "run_accession",
                "blocker_type",
            ]
        )
    raw["pilot_status"] = raw["tier"].map(
        lambda value: "authorized_rfc_only_not_downloaded" if str(value) == "P3_minimal_etl_rfc" else "not_authorized"
    )
    raw["pilot_reason"] = raw["recommended_action"]
    raw["pilot_sample_count"] = raw["run_accession"].map(lambda value: min(3, len([item for item in str(value).split(";") if item])))
    return raw[["candidate_id", "dataset", "pilot_status", "pilot_reason", "pilot_sample_count", "run_accession", "blocker_type"]]


def decide_state(candidate_table: pd.DataFrame, refresh_rounds: int) -> dict[str, Any]:
    p1 = candidate_table[candidate_table["tier"].astype(str).eq("P1_processed_headline")]
    p3 = candidate_table[candidate_table["tier"].astype(str).eq("P3_minimal_etl_rfc")]
    headline = candidate_table[candidate_table["headline_allowed"].map(boolish)]
    n_scientific = int(len(p1) + len(p3) + len(headline))
    if n_scientific:
        status = "continue"
        decision = "v12_candidate_found_requires_matrix_or_minimal_etl_gate"
        next_action = "Run adapter smoke or minimal ETL RFC approval; no training until matrix gate passes."
    elif refresh_rounds >= 3:
        status = "benchmark_redefinition_success"
        decision = "v12_benchmark_redefinition_success_no_headline_data_after_three_rounds"
        next_action = "Adopt benchmark redefinition RFC; do not run autoresearch on unsupported old104+ headline."
    else:
        status = "continue"
        decision = "v12_continue_candidate_refresh"
        next_action = "Complete three refresh rounds before declaring benchmark redefinition."
    tier_counts = Counter(candidate_table["tier"].fillna("").astype(str))
    return {
        "timestamp": utc_now(),
        "loop_version": "v12.0",
        "status": status,
        "loop_decision": decision,
        "refresh_rounds_evaluated": int(refresh_rounds),
        "next_action": next_action,
        "success_gates": {
            "scientific_data_success": {
                "dataset_gate": "P1 processed headline or P3 minimal-ETL-authorized old bulk target tissue",
                "gse121141_old104_mae_weeks_max": 65.386,
                "gse121141_all_age_mae_weeks_max": 40.033,
                "groupkfold_mae_weeks_max": 25.767,
                "random_label_sanity": "abs(r)<0.2 and random-like MAE",
                "shuffled_cr_sanity": "AUC 0.4-0.6 and Cohen_d decreases",
            },
            "benchmark_redefinition_success": "Three refresh rounds produce no P1/P3 headline data and an auditable RFC is generated.",
        },
        "failure_gates": {
            "no_redefinition_rfc": "No P1/P3 and no auditable benchmark redefinition",
            "three_headline_matrices_still_fail": "Three valid headline matrices leave old104+ MAE near 70w+",
            "sanity_failure": "random-label or shuffled-CR sanity fails",
            "unresolved_network_or_schema_blocker": "Cannot close official metadata/API provenance",
        },
        "metrics": {
            "n_candidate_rows": int(len(candidate_table)),
            "n_p1_processed_headline": int(tier_counts.get("P1_processed_headline", 0)),
            "n_p3_minimal_etl_rfc": int(tier_counts.get("P3_minimal_etl_rfc", 0)),
            "n_p2_auxiliary": int(sum(count for tier, count in tier_counts.items() if tier.startswith("P2"))),
            "n_p4_blocked": int(sum(count for tier, count in tier_counts.items() if tier.startswith("P4"))),
            "n_p5_redefinition_evidence": int(tier_counts.get("P5_redefinition_evidence", 0)),
            "n_headline_allowed": int(len(headline)),
            "v75_groupkfold_mae_weeks": V75_BASELINES["groupkfold_mae_weeks"],
            "v75_gse121141_all_age_mae_weeks": V75_BASELINES["gse121141_all_age_mae_weeks"],
            "v75_gse121141_old104_mae_weeks": V75_BASELINES["gse121141_old104_mae_weeks"],
        },
    }


def write_benchmark_redefinition_rfc(path: Path, candidate_table: pd.DataFrame, state: dict[str, Any]) -> None:
    blocker_cols = ["tier", "blocker_type"]
    blockers = (
        candidate_table[blocker_cols]
        .fillna("")
        .groupby(blocker_cols, dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values("count", ascending=False)
        if not candidate_table.empty
        else pd.DataFrame()
    )
    lines = [
        "# v12 Benchmark Redefinition RFC",
        "",
        f"Date: {utc_now()}",
        "",
        "## Decision",
        "",
        "Redefine the headline benchmark if no P1/P3 old bulk brain_cortex/heart/lung dataset is found after three v12 refresh rounds.",
        "",
        "## Evidence",
        "",
        f"- Refresh rounds evaluated: {state['refresh_rounds_evaluated']}",
        f"- Candidate rows evaluated: {state['metrics']['n_candidate_rows']}",
        f"- P1 processed headline candidates: {state['metrics']['n_p1_processed_headline']}",
        f"- P3 minimal ETL RFC candidates: {state['metrics']['n_p3_minimal_etl_rfc']}",
        f"- Headline-allowed candidates: {state['metrics']['n_headline_allowed']}",
        "",
        "## Blocker Summary",
        "",
        md_table(blockers, max_rows=40),
        "",
        "## New Headline Benchmark",
        "",
        "- Primary chronological benchmark is limited to samples whose held-out fold has same-tissue and age-range support in training data.",
        "- A held-out sample is support-covered when train data contains at least 10 same-tissue samples and train same-tissue max age is within 8 weeks of the held-out age.",
        "- GSE121141 old104+ brain_cortex/heart/lung remains reported, but is downgraded to a stress-test/blocker metric rather than pass/fail headline.",
        "- GSE80672 CR/rapamycin biological-age metrics remain allowed only from real held-out predictions and shuffled-intervention sanity checks.",
        "- Old target-tissue data acquisition moves to v13 as new data generation, collaboration, or explicitly approved minimal ETL.",
        "",
        "## Non-Negotiable Constraints",
        "",
        "- No autoresearch on unsupported old104+ headline data.",
        "- No human clock CpG mapping.",
        "- No dummy AUC.",
        "- No full-data feature selection or preprocessing outside train fold/train dataset.",
        "- No raw FASTQ download without a separate minimal ETL approval.",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_report(path: Path, candidate_table: pd.DataFrame, state: dict[str, Any], rfc_path: Path) -> None:
    view_cols = [
        "refresh_round",
        "dataset",
        "tier",
        "tissue",
        "processed_schema_guess",
        "common_regions_estimate",
        "gate_status",
        "blocker_type",
        "recommended_action",
    ]
    view = candidate_table[[col for col in view_cols if col in candidate_table.columns]].copy()
    lines = [
        "# v12 RALPH Loop Success or Failure Report",
        "",
        f"Date: {utc_now()}",
        "",
        "## Summary",
        "",
        "v12 closed the three-round data strategy review using cached official v9-v11 refreshes, raw RFC verification, and v11.x conversion gates. It did not download FASTQ, build new matrices, run benchmarks, or start autoresearch.",
        "",
        f"- Status: `{state['status']}`",
        f"- Decision: `{state['loop_decision']}`",
        f"- Next action: {state['next_action']}",
        f"- Benchmark redefinition RFC: `{rfc_path.relative_to(ROOT)}`",
        "",
        "## Metrics",
        "",
        md_table(pd.DataFrame([state["metrics"]])),
        "",
        "## Candidate Gate Table",
        "",
        md_table(view, max_rows=80),
        "",
        "## Success / Failure Gate",
        "",
        "- Scientific data success was not reached because no P1 processed headline or P3 authorized minimal-ETL dataset passed the data gate.",
        "- Benchmark redefinition success is reached when three refresh rounds are complete and the RFC is auditable.",
        "- If benchmark redefinition is rejected, the next project phase is v13 data generation/collaboration/minimal ETL, not autoresearch.",
        "",
        "## Outputs",
        "",
        "- `results/ralph_v12_loop/ralph_iteration_log.jsonl`",
        "- `results/ralph_v12_loop/ralph_decision_state.json`",
        "- `results/ralph_v12_loop/candidate_gate_table.csv`",
        "- `results/ralph_v12_loop/network_resolution_log.jsonl`",
        "- `results/ralph_v12_loop/candidate_smoke_manifest.csv`",
        "- `results/ralph_v12_loop/minimal_etl_rfc_manifest.csv`",
        "- `results/ralph_v12_loop/benchmark_redefinition_rfc.md`",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out_dir", default=str(OUT_DIR))
    parser.add_argument("--refresh_rounds", type=int, default=3)
    parser.add_argument(
        "--cached_only",
        action="store_true",
        default=True,
        help="Use existing official v9-v11 artifacts only. This is the default safe behavior.",
    )
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    inventory = combine_inventory()
    rows: list[dict[str, Any]] = []
    rows.extend(normalize_v10_like(ROOT / "results" / "ralph_v9_loop" / "candidate_gate_table.csv", "round1_cached_v9_official_geo", inventory))
    rows.extend(normalize_v10_like(ROOT / "results" / "ralph_v10_loop" / "candidate_gate_table.csv", "round1_cached_v10_broadened_geo", inventory))
    rows.extend(normalize_v10_like(ROOT / "results" / "ralph_v11_data_strategy" / "verified_new_gse_gate_table.csv", "round2_cached_pubmed_sra_verified_gse", inventory))
    rows.extend(normalize_v11x(ROOT / "results" / "ralph_v11_x_loop" / "candidate_gate_table.csv", inventory))
    rows.extend(
        normalize_raw_biosample(
            ROOT / "results" / "ralph_v11_data_strategy" / "raw_biosample_verification_summary.csv",
            ROOT / "results" / "ralph_v11_data_strategy" / "raw_biosample_verification_manifest.csv",
        )
    )

    candidate_table = pd.DataFrame(rows)
    if candidate_table.empty:
        candidate_table = pd.DataFrame(columns=STANDARD_COLUMNS)
    for col in STANDARD_COLUMNS:
        if col not in candidate_table.columns:
            candidate_table[col] = ""
    candidate_table = candidate_table[STANDARD_COLUMNS].drop_duplicates("candidate_id", keep="last").reset_index(drop=True)
    candidate_table.to_csv(out_dir / "candidate_gate_table.csv", index=False)

    smoke = aggregate_smoke_manifests()
    smoke.to_csv(out_dir / "candidate_smoke_manifest.csv", index=False)
    minimal = build_minimal_etl_manifest(candidate_table)
    minimal.to_csv(out_dir / "minimal_etl_rfc_manifest.csv", index=False)

    network_log = out_dir / "network_resolution_log.jsonl"
    append_jsonl(
        network_log,
        {
            "timestamp": utc_now(),
            "blocker_type": "official_database_rules_loaded",
            "attempted_url": "",
            "source_doc": " | ".join(OFFICIAL_SOURCE_DOCS.values()),
            "resolution": "v12_cached_three_round_review_uses_prior_official_outputs",
            "applied_rule": "no_large_download_no_training_until_p1_or_p3_gate",
        },
    )
    append_jsonl(
        out_dir / "ralph_iteration_log.jsonl",
        {
            "timestamp": utc_now(),
            "loop": "v12_data_strategy_benchmark_redefinition",
            "action": "aggregate_cached_three_rounds",
            "cached_only": bool(args.cached_only),
            "n_candidate_rows": int(len(candidate_table)),
        },
    )

    state = decide_state(candidate_table, args.refresh_rounds)
    write_json(out_dir / "ralph_decision_state.json", state)
    write_benchmark_redefinition_rfc(out_dir / "benchmark_redefinition_rfc.md", candidate_table, state)
    write_report(out_dir / "v12_success_or_failure_report.md", candidate_table, state, out_dir / "benchmark_redefinition_rfc.md")
    append_jsonl(
        out_dir / "ralph_iteration_log.jsonl",
        {
            "timestamp": utc_now(),
            "loop": "v12_data_strategy_benchmark_redefinition",
            "action": "decision",
            "status": state["status"],
            "decision": state["loop_decision"],
            "metrics": state["metrics"],
        },
    )
    print(json.dumps(state, indent=2, sort_keys=True, default=str))


if __name__ == "__main__":
    main()

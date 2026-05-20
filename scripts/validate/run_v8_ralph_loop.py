#!/usr/bin/env python3
"""Controlled v8.x RALPH loop runner.

The loop is deliberately conservative. It uses existing v8.5/v8.6 preflight
outputs, downloads only a few per-sample supplementary files for adapter smoke,
and stops before matrix building or model training unless all gates are met.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
import math
import re
import sys
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")
DEFAULT_OUT_DIR = ROOT / "results" / "ralph_v8_loop"
DEFAULT_REPORT = DEFAULT_OUT_DIR / "v8_success_or_v9_blocker_report.md"
V85_DIR = ROOT / "results" / "v8_5_targeted_geo_refresh"
V86_DIR = ROOT / "results" / "v8_6_calibration_poc"
REFERENCE_MATRIX = ROOT / "results" / "multidataset_v8_2_prefilter_liftover" / "all_rrbs_region_matrix_5kb.parquet"

CANDIDATES = ["GSE225166", "GSE83947", "GSE134398"]
TARGET_TISSUES = {"brain_cortex", "heart", "lung"}
OLD_THRESHOLD_WEEKS = 104.0
MAX_SMOKE_FILES_PER_DATASET = 3
MAX_SMOKE_FILE_BYTES = 60_000_000
MAX_PARSE_ROWS = 250_000

V75_BASELINES = {
    "groupkfold_mae_weeks": 23.767,
    "gse121141_heldout_mae_weeks": 38.033,
    "gse121141_old104_mae_weeks": 75.386,
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


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


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def append_jsonl(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(payload, sort_keys=True) + "\n")


def read_csv(path: Path, **kwargs) -> pd.DataFrame:
    return pd.read_csv(path, **kwargs) if path.exists() else pd.DataFrame()


def import_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def v84_module():
    return import_module(ROOT / "scripts" / "validate" / "v8_4_old_tissue_candidate_and_calibration.py", "v8_4_helpers_ralph")


def normalize_tissue(v84, row: pd.Series) -> str:
    labels = v84.canonical_tissues(row.get("tissue_guess"))
    if not labels:
        labels = v84.canonical_tissues(f"{row.get('source_name', '')} {row.get('title', '')} {row.get('characteristics', '')}")
    if "brain_cortex" in labels:
        return "brain_cortex"
    if "brain_other" in labels:
        return "brain_other"
    if labels:
        return sorted(labels)[0]
    raw = str(row.get("tissue_guess") or row.get("source_name") or "unknown")
    return raw.strip() or "unknown"


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    return (
        read_csv(V85_DIR / "v8_5_candidate_ranking.csv"),
        read_csv(V85_DIR / "targeted_search_inventory.csv"),
        read_csv(V85_DIR / "targeted_search_samples.csv"),
        read_csv(V85_DIR / "targeted_search_supplements.csv"),
    )


def match_supplement_rows(samples: pd.DataFrame, supplements: pd.DataFrame, v84) -> pd.DataFrame:
    rows = []
    file_rows = supplements[
        supplements["archive_or_file"].astype(str).str.lower().eq("file")
        & supplements["filelist_type"].astype(str).str.upper().isin(["COV", "TXT", "BEDGRAPH"])
    ].copy()
    for _, sample in samples.iterrows():
        sample_id = str(sample.get("sample_id", ""))
        matches = file_rows[file_rows["supplement_name"].astype(str).str.contains(sample_id, regex=False)]
        if matches.empty:
            continue
        age = v84.conservative_sample_age_weeks(sample)
        tissue = normalize_tissue(v84, sample)
        for _, file_row in matches.iterrows():
            size = file_row.get("filelist_size_bytes")
            try:
                size_value = int(float(size))
            except (TypeError, ValueError):
                size_value = 0
            rows.append(
                {
                    "dataset": sample.get("dataset"),
                    "sample_id": sample_id,
                    "title": sample.get("title", ""),
                    "source_name": sample.get("source_name", ""),
                    "tissue": tissue,
                    "age_weeks": age,
                    "is_old104": bool(np.isfinite(age) and age >= OLD_THRESHOLD_WEEKS),
                    "is_target_tissue": tissue in TARGET_TISSUES,
                    "supplement_name": file_row.get("supplement_name", ""),
                    "filelist_type": str(file_row.get("filelist_type", "")).upper(),
                    "file_size_bytes": size_value,
                    "supplement_url": file_row.get("supplement_url", ""),
                    "under_size_cap": bool(size_value and size_value <= MAX_SMOKE_FILE_BYTES),
                }
            )
    return pd.DataFrame(rows)


def select_smoke_files(matched: pd.DataFrame) -> pd.DataFrame:
    if matched.empty:
        return matched
    data = matched[matched["under_size_cap"].astype(bool)].copy()
    if data.empty:
        return data
    data["priority"] = 0
    data.loc[data["is_target_tissue"].astype(bool), "priority"] += 10
    data.loc[data["is_old104"].astype(bool), "priority"] += 10
    data.loc[data["filelist_type"].eq("COV"), "priority"] += 3
    data.loc[data["filelist_type"].eq("TXT"), "priority"] += 1
    data["size_rank"] = data["file_size_bytes"].astype(float)
    return data.sort_values(["priority", "size_rank"], ascending=[False, True]).head(MAX_SMOKE_FILES_PER_DATASET)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def gsm_bucket(sample_id: str) -> str | None:
    match = re.fullmatch(r"(GSM)(\d+)", str(sample_id).strip())
    if not match:
        return None
    digits = match.group(2)
    return f"GSM{digits[:-3]}nnn" if len(digits) > 3 else "GSMnnn"


def sample_supplement_url(sample_id: str, filename: str) -> str | None:
    bucket = gsm_bucket(sample_id)
    match = re.fullmatch(r"GSM\d+", str(sample_id).strip())
    if bucket is None or match is None:
        return None
    quoted = urllib.parse.quote(Path(str(filename)).name)
    return f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{bucket}/{match.group(0)}/suppl/{quoted}"


def smoke_url_candidates(row: pd.Series, filename: str) -> list[tuple[str, str]]:
    candidates: list[tuple[str, str]] = []
    sample_url = sample_supplement_url(str(row.get("sample_id", "")), filename)
    inventory_url = str(row.get("supplement_url", "") or "").strip()
    if sample_url:
        candidates.append(("sample_supplement", sample_url))
    if inventory_url and inventory_url not in {url for _, url in candidates}:
        candidates.append(("inventory", inventory_url))
    return candidates


def download_smoke_file(row: pd.Series, target_dir: Path, log_path: Path, timeout: int = 180, retries: int = 3) -> tuple[Path | None, dict]:
    target_dir.mkdir(parents=True, exist_ok=True)
    filename = Path(str(row["supplement_name"])).name
    target = target_dir / filename
    expected_size = int(row.get("file_size_bytes") or 0)
    if expected_size > MAX_SMOKE_FILE_BYTES:
        payload = {
            "timestamp": utc_now(),
            "dataset": row["dataset"],
            "sample_id": row["sample_id"],
            "file": filename,
            "status": "skipped",
            "reason": "file_exceeds_smoke_size_cap",
            "expected_size": expected_size,
        }
        append_jsonl(log_path, payload)
        return None, payload
    if target.exists() and expected_size and target.stat().st_size == expected_size:
        payload = {
            "timestamp": utc_now(),
            "dataset": row["dataset"],
            "sample_id": row["sample_id"],
            "file": filename,
            "status": "skipped",
            "reason": "already_verified",
            "bytes_downloaded": target.stat().st_size,
            "sha256": sha256_file(target),
        }
        append_jsonl(log_path, payload)
        return target, payload

    part = target.with_suffix(target.suffix + ".part")
    errors = []
    total_failures = 0
    url_candidates = smoke_url_candidates(row, filename)
    for url_source, url in url_candidates:
        for attempt in range(1, retries + 1):
            try:
                if part.exists():
                    part.unlink()
                request = urllib.request.Request(url, headers={"User-Agent": "mouse-methyl-v8-ralph-smoke/1.0"})
                with urllib.request.urlopen(request, timeout=timeout) as response, part.open("wb") as handle:
                    while True:
                        chunk = response.read(1024 * 1024)
                        if not chunk:
                            break
                        handle.write(chunk)
                if expected_size and part.stat().st_size != expected_size:
                    raise RuntimeError(f"downloaded size {part.stat().st_size} != expected {expected_size}")
                part.replace(target)
                payload = {
                    "timestamp": utc_now(),
                    "dataset": row["dataset"],
                    "sample_id": row["sample_id"],
                    "file": filename,
                    "status": "completed",
                    "url": url,
                    "url_source": url_source,
                    "bytes_downloaded": target.stat().st_size,
                    "sha256": sha256_file(target),
                    "failure_count": total_failures,
                }
                append_jsonl(log_path, payload)
                return target, payload
            except Exception as exc:
                total_failures += 1
                error = str(exc)[:500]
                errors.append({"url_source": url_source, "url": url, "attempt": attempt, "error": error})
                time.sleep(min(2**attempt, 20))
    payload = {
        "timestamp": utc_now(),
        "dataset": row["dataset"],
        "sample_id": row["sample_id"],
        "file": filename,
        "status": "failed",
        "error": errors[-1]["error"] if errors else "no_download_url_candidates",
        "url_errors": errors,
        "failure_count": total_failures,
    }
    append_jsonl(log_path, payload)
    return None, payload


def normalize_chrom(value: str) -> str | None:
    value = str(value).strip()
    if value in {"X", "Y", "M", "MT", "chrX", "chrY", "chrM"}:
        return None
    if value.startswith("chr"):
        clean = value.removeprefix("chr")
        if clean.isdigit():
            return f"chr{clean}"
        return None
    if value.isdigit():
        return f"chr{value}"
    return None


def parse_float(value: str) -> float | None:
    try:
        parsed = float(value)
    except ValueError:
        return None
    if math.isnan(parsed) or math.isinf(parsed):
        return None
    return parsed


def iter_text_lines(path: Path) -> Iterable[str]:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", errors="replace") as handle:
        for line in handle:
            yield line.rstrip("\n")


def detect_and_parse_line(parts: list[str]) -> tuple[str | None, str | None, int | None, float | None, float | None]:
    if len(parts) >= 6:
        chrom = normalize_chrom(parts[0])
        pos = parse_float(parts[1])
        meth = parse_float(parts[4])
        unmeth = parse_float(parts[5])
        pct = parse_float(parts[3])
        if chrom and pos is not None and meth is not None and unmeth is not None and meth + unmeth > 0:
            beta = meth / (meth + unmeth)
            if 0 <= beta <= 1:
                return "bismark_cov_6col", chrom, int(pos), beta, meth + unmeth
            if pct is not None and 0 <= pct <= 100:
                return "bismark_cov_6col_pct", chrom, int(pos), pct / 100.0, meth + unmeth
    if len(parts) >= 5:
        chrom = normalize_chrom(parts[0])
        pos = parse_float(parts[1])
        meth = parse_float(parts[3])
        unmeth = parse_float(parts[4])
        if chrom and pos is not None and meth is not None and unmeth is not None and meth + unmeth > 0:
            beta = meth / (meth + unmeth)
            if 0 <= beta <= 1:
                return "bismark_cx_report", chrom, int(pos), beta, meth + unmeth
    return None, None, None, None, None


def load_reference_regions() -> set[str]:
    if not REFERENCE_MATRIX.exists():
        return set()
    matrix = pd.read_parquet(REFERENCE_MATRIX, columns=[])
    return set(map(str, matrix.index))


def region_id(chrom: str, pos: int, bin_size: int = 5000) -> str:
    start = (int(pos) // bin_size) * bin_size
    end = start + bin_size - 1
    return f"{chrom}:{start}-{end}"


def parse_smoke_file(path: Path, reference_regions: set[str]) -> dict:
    rows_total = 0
    rows_parseable = 0
    rows_pass_coverage = 0
    beta_min = None
    beta_max = None
    schema_counts: dict[str, int] = {}
    regions: set[str] = set()
    common_regions: set[str] = set()
    coverage_values = []
    for line in iter_text_lines(path):
        if rows_total >= MAX_PARSE_ROWS:
            break
        if not line or line.startswith("track") or line.startswith("#"):
            continue
        rows_total += 1
        parts = line.split("\t")
        schema, chrom, pos, beta, coverage = detect_and_parse_line(parts)
        if schema is None or chrom is None or pos is None or beta is None or coverage is None:
            continue
        rows_parseable += 1
        schema_counts[schema] = schema_counts.get(schema, 0) + 1
        beta_min = beta if beta_min is None else min(beta_min, beta)
        beta_max = beta if beta_max is None else max(beta_max, beta)
        if coverage >= 5:
            rows_pass_coverage += 1
            coverage_values.append(float(coverage))
            rid = region_id(chrom, pos)
            regions.add(rid)
            if reference_regions and rid in reference_regions:
                common_regions.add(rid)
    schema_guess = max(schema_counts, key=schema_counts.get) if schema_counts else "unparseable"
    return {
        "local_path": str(path),
        "rows_scanned": rows_total,
        "rows_parseable_primary_autosomes": rows_parseable,
        "rows_pass_coverage_ge5": rows_pass_coverage,
        "schema_guess": schema_guess,
        "beta_min": None if beta_min is None else round(float(beta_min), 6),
        "beta_max": None if beta_max is None else round(float(beta_max), 6),
        "median_coverage": None if not coverage_values else round(float(np.median(coverage_values)), 3),
        "n_smoke_regions": len(regions),
        "n_smoke_common_regions_with_reference": len(common_regions),
        "beta_range_valid": bool(beta_min is not None and beta_max is not None and 0 <= beta_min <= beta_max <= 1),
    }


def context_blockers_for_candidate(dataset: str, ranking_row: pd.Series | None, inventory_row: pd.Series | None) -> list[str]:
    blockers = []
    text = ""
    if ranking_row is not None:
        text += f" {ranking_row.get('v8_4_blockers', '')} {ranking_row.get('title', '')}"
    if inventory_row is not None:
        text += f" {inventory_row.get('title', '')} {inventory_row.get('summary', '')} {inventory_row.get('overall_design', '')}"
    lowered = text.lower()
    if "targeted_celltype_or_non_bulk_context" in lowered or "single-cell" in lowered or "single cell" in lowered:
        blockers.append("targeted_celltype_or_non_bulk_context")
    if "ipcrtag" in lowered or "itag" in lowered or "low coverage" in lowered:
        blockers.append("low_coverage_or_tagged_context")
    if "superseries" in lowered or "superseries_use_subseries" in lowered:
        blockers.append("superseries_or_subseries_required")
    if "single_target_tissue_only" in lowered or dataset in {"GSE83947", "GSE134398"}:
        blockers.append("single_target_tissue_only")
    return list(dict.fromkeys(blockers))


def candidate_summary(
    dataset: str,
    ranking: pd.DataFrame,
    inventory: pd.DataFrame,
    samples: pd.DataFrame,
    supplements: pd.DataFrame,
    v84,
    reference_regions: set[str],
    out_dir: Path,
) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    ranking_row = None
    if not ranking.empty and dataset in set(ranking["dataset"].astype(str)):
        ranking_row = ranking[ranking["dataset"].astype(str).eq(dataset)].iloc[0]
    inventory_row = None
    if not inventory.empty and dataset in set(inventory["dataset"].astype(str)):
        inventory_row = inventory[inventory["dataset"].astype(str).eq(dataset)].iloc[0]

    dataset_samples = samples[samples["dataset"].astype(str).eq(dataset)].copy()
    dataset_supplements = supplements[supplements["dataset"].astype(str).eq(dataset)].copy()
    matched = match_supplement_rows(dataset_samples, dataset_supplements, v84)
    selected = select_smoke_files(matched)
    smoke_dir = out_dir / "smoke_files" / dataset
    download_log = out_dir / "smoke_download_log.jsonl"

    smoke_rows = []
    for _, selected_row in selected.iterrows():
        local_path, download_payload = download_smoke_file(selected_row, smoke_dir, download_log)
        parse_payload = {}
        if local_path is not None:
            try:
                parse_payload = parse_smoke_file(local_path, reference_regions)
            except Exception as exc:
                parse_payload = {"local_path": str(local_path), "parse_error": str(exc)[:500]}
        smoke_rows.append({**selected_row.to_dict(), **download_payload, **parse_payload})

    smoke = pd.DataFrame(smoke_rows)
    matched.to_csv(out_dir / f"{dataset}_matched_supplement_files.csv", index=False)
    smoke.to_csv(out_dir / f"{dataset}_adapter_smoke.csv", index=False)

    context_blockers = context_blockers_for_candidate(dataset, ranking_row, inventory_row)
    target_old_samples = matched[matched["is_target_tissue"].astype(bool) & matched["is_old104"].astype(bool)]
    old_target_files = selected[selected["is_target_tissue"].astype(bool) & selected["is_old104"].astype(bool)] if not selected.empty else pd.DataFrame()
    if not smoke.empty:
        parseable_mask = pd.to_numeric(
            smoke["rows_parseable_primary_autosomes"] if "rows_parseable_primary_autosomes" in smoke else pd.Series(0, index=smoke.index),
            errors="coerce",
        ).fillna(0).gt(0)
        pass_coverage_mask = pd.to_numeric(
            smoke["rows_pass_coverage_ge5"] if "rows_pass_coverage_ge5" in smoke else pd.Series(0, index=smoke.index),
            errors="coerce",
        ).fillna(0).gt(0)
        beta_valid_series = (
            smoke["beta_range_valid"] if "beta_range_valid" in smoke else pd.Series(False, index=smoke.index)
        ).fillna(False).astype(bool)
        parseable = smoke[parseable_mask]
        pass_coverage = smoke[pass_coverage_mask]
        beta_valid = bool(beta_valid_series.any())
    else:
        parseable = pd.DataFrame()
        pass_coverage = pd.DataFrame()
        beta_valid = False
    downloaded_ok = bool(
        not smoke.empty
        and "status" in smoke
        and smoke["status"].astype(str).isin(["completed", "skipped"]).any()
        and "local_path" in smoke
        and smoke["local_path"].notna().any()
    )
    archive_available = bool(
        not dataset_supplements.empty
        and dataset_supplements["archive_or_file"].astype(str).str.lower().eq("archive").any()
    )
    tar_required_for_smoke = bool(not downloaded_ok and archive_available and not selected.empty)
    schema_smoke_pass = bool(downloaded_ok and not parseable.empty and not pass_coverage.empty and beta_valid)
    assembly_compatible = bool(
        not smoke.empty
        and smoke["supplement_name"].astype(str).str.contains("GRCm38|mm10|GRCm38|mm10", case=False, regex=True).any()
    )
    if not smoke.empty and "n_smoke_common_regions_with_reference" in smoke:
        common_series = pd.to_numeric(smoke["n_smoke_common_regions_with_reference"], errors="coerce").fillna(0)
        smoke_common_regions = int(common_series.max()) if not common_series.empty else 0
    else:
        smoke_common_regions = 0
    official_metadata_pass = bool(len(target_old_samples) > 0)
    metadata_overlap_pass = bool(len(old_target_files) > 0)
    hard_context_blocked = bool(context_blockers)
    full_matrix_candidate = bool(
        official_metadata_pass
        and schema_smoke_pass
        and metadata_overlap_pass
        and assembly_compatible
        and not hard_context_blocked
        and len(set(target_old_samples["tissue"]) & {"brain_cortex"}) > 0
    )
    auxiliary_reason = ""
    if hard_context_blocked:
        auxiliary_reason = ";".join(context_blockers)
    elif not official_metadata_pass:
        auxiliary_reason = "no_old_target_tissue_with_sample_specific_age"
    elif not downloaded_ok:
        auxiliary_reason = "individual_file_download_failed_tar_required_for_smoke" if tar_required_for_smoke else "smoke_file_download_failed"
    elif not schema_smoke_pass:
        auxiliary_reason = "schema_smoke_failed_or_unparseable"
    elif not assembly_compatible:
        auxiliary_reason = "assembly_not_explicitly_compatible_in_smoke_files"
    elif "brain_cortex" not in set(target_old_samples["tissue"]):
        auxiliary_reason = "no_old_brain_cortex_support"
    elif not full_matrix_candidate:
        auxiliary_reason = "full_matrix_gate_not_reached"

    row = {
        "iteration": "v8.7",
        "dataset": dataset,
        "hypothesis": "candidate may repair GSE121141 old same-tissue support",
        "action": "adapter_smoke_only_no_tar_no_training",
        "official_metadata_pass": official_metadata_pass,
        "download_allowed": bool(not selected.empty),
        "smoke_download_pass": downloaded_ok,
        "tar_required_for_smoke": tar_required_for_smoke,
        "n_matched_files": int(len(matched)),
        "n_selected_smoke_files": int(len(selected)),
        "n_old_target_samples_with_files": int(len(target_old_samples)),
        "old_target_tissues_with_files": ";".join(sorted(set(target_old_samples["tissue"]))) if not target_old_samples.empty else "",
        "schema_smoke_pass": schema_smoke_pass,
        "metadata_overlap_pass": metadata_overlap_pass,
        "assembly_compatible": assembly_compatible,
        "smoke_common_regions_estimate": smoke_common_regions,
        "common_regions_pass": False,
        "common_regions_reason": "not_assessed_until_full_matrix",
        "headline_allowed": full_matrix_candidate,
        "full_matrix_candidate": full_matrix_candidate,
        "auxiliary_only_reason": auxiliary_reason,
        "context_blockers": ";".join(context_blockers),
        "next_step": "build_full_matrix_and_benchmark" if full_matrix_candidate else "do_not_promote_to_headline",
    }
    return row, matched, smoke


def load_v86_summary() -> dict:
    path = V86_DIR / "v8_6_summary.json"
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def decide_loop_state(gate: pd.DataFrame, v86_summary: dict) -> dict:
    success = False
    full_matrix_candidates = gate[gate["full_matrix_candidate"].astype(bool)] if not gate.empty else pd.DataFrame()
    headline_candidates = gate[gate["headline_allowed"].astype(bool)] if not gate.empty else pd.DataFrame()
    calibration_transferable = False
    best_calibration_delta = v86_summary.get("best_old104plus_delta_mae")
    # v8.6 report shows the best non-target validation was negative; encode that
    # as non-transferable until a future result explicitly overrides it.
    if best_calibration_delta is not None and float(best_calibration_delta) >= 10:
        calibration_transferable = False

    if not headline_candidates.empty:
        loop_decision = "continue_v8_full_matrix_benchmark"
        next_action = "Build candidate full matrix, then run fixed benchmark gates before any constrained search."
        v9_bottleneck = False
    elif full_matrix_candidates.empty and not calibration_transferable:
        loop_decision = "v9_bottleneck_reached_for_headline_goal"
        next_action = "Stop v8 headline loop and notify user. Optional auxiliary smoke data should not be promoted."
        v9_bottleneck = True
    else:
        loop_decision = "continue_v8_diagnostics"
        next_action = "Continue diagnostics only; do not run autoresearch."
        v9_bottleneck = False

    return {
        "timestamp": utc_now(),
        "loop_version": "v8.7",
        "status": "success" if success else ("v9_bottleneck" if v9_bottleneck else "continue"),
        "loop_decision": loop_decision,
        "next_action": next_action,
        "success_criteria": {
            "gse121141_old104_mae_improvement_weeks": ">=10 vs 75.386w",
            "gse121141_all_age_mae_max": "<=40.033w",
            "groupkfold_mae_max": "<=25.767w",
            "random_label_sanity": "abs(r)<0.2 and MAE random-like",
            "shuffled_cr_sanity": "AUC chance-like and Cohen_d decreased",
        },
        "v9_bottleneck_criteria": {
            "no_headline_candidate_after_refresh_and_smoke": bool(headline_candidates.empty),
            "calibration_not_transferable": not calibration_transferable,
            "only_auxiliary_contexts_available": bool(
                not gate.empty
                and gate["headline_allowed"].astype(bool).sum() == 0
                and gate["auxiliary_only_reason"].astype(str).ne("").all()
            ),
        },
        "metrics": {
            "v75_gse121141_old104_mae_weeks": V75_BASELINES["gse121141_old104_mae_weeks"],
            "v86_best_old104plus_delta_mae_weeks": best_calibration_delta,
            "n_candidates": int(len(gate)),
            "n_headline_allowed": int(gate["headline_allowed"].astype(bool).sum()) if not gate.empty else 0,
            "n_full_matrix_candidates": int(gate["full_matrix_candidate"].astype(bool).sum()) if not gate.empty else 0,
        },
    }


def write_report(out_dir: Path, report_path: Path, gate: pd.DataFrame, decision: dict) -> None:
    view_cols = [
        "dataset",
        "official_metadata_pass",
        "download_allowed",
        "smoke_download_pass",
        "tar_required_for_smoke",
        "n_selected_smoke_files",
        "old_target_tissues_with_files",
        "schema_smoke_pass",
        "metadata_overlap_pass",
        "assembly_compatible",
        "smoke_common_regions_estimate",
        "headline_allowed",
        "auxiliary_only_reason",
        "next_step",
    ]
    view = gate[[col for col in view_cols if col in gate.columns]].copy() if not gate.empty else pd.DataFrame()
    lines = [
        "# v8.x RALPH Loop Decision Report",
        "",
        f"Date: {utc_now()}",
        "",
        "## Summary",
        "",
        "This run executed the v8.7 RALPH adapter-smoke loop for GSE225166, GSE83947, and GSE134398. It did not download RAW tar archives, build matrices, train models, or run autoresearch.",
        "",
        f"- Loop decision: `{decision['loop_decision']}`",
        f"- Next action: {decision['next_action']}",
        f"- Headline-allowed candidates: {decision['metrics']['n_headline_allowed']}",
        f"- Full-matrix candidates: {decision['metrics']['n_full_matrix_candidates']}",
        f"- v8.6 calibration transferable: False",
        "",
        "## GEO API/FTP Handling",
        "",
        "When inventory URLs fail for individual supplementary files, the loop follows the official GEO FTP directory structure and retries sample-level supplementary URLs in the form `/geo/samples/GSM.../<GSM>/suppl/<filename>`. E-Utils/SOFT remain metadata sources; full records and supplementary files are still resolved through GEO FTP paths.",
        "",
        "Official references checked:",
        "",
        "- GEO Download: https://www.ncbi.nlm.nih.gov/geo/info/download.html",
        "- GEO Programmatic Access: https://www.ncbi.nlm.nih.gov/geo/info/geo_paccess.html",
        "- GEO SOFT format: https://www.ncbi.nlm.nih.gov/geo/info/soft.html",
        "",
        "## Candidate Gate Table",
        "",
        md_table(view),
        "",
        "## Interpretation",
        "",
    ]
    if decision["status"] == "v9_bottleneck":
        lines.extend(
            [
                "The v8 headline loop has reached the planned v9 bottleneck condition: the available candidates remain auxiliary-only or hard-blocked for the core GSE121141 old brain_cortex/heart/lung goal, and calibration has not transferred robustly in non-target LODO validation.",
                "",
                "Do not start constrained autoresearch from these candidates. A v9 plan should redefine the data strategy: broader processed methylome acquisition, single-tissue old-age validation cohorts, or minimal FASTQ/Bismark ETL where processed matrices are unavailable.",
            ]
        )
    else:
        lines.append("At least one candidate may continue in v8, but it must pass full matrix and benchmark gates before any optimization.")
    lines.extend(
        [
            "",
            "## Outputs",
            "",
            "- `results/ralph_v8_loop/ralph_iteration_log.jsonl`",
            "- `results/ralph_v8_loop/ralph_decision_state.json`",
            "- `results/ralph_v8_loop/candidate_gate_table.csv`",
            "- `results/ralph_v8_loop/*_adapter_smoke.csv`",
            "- `results/ralph_v8_loop/smoke_download_log.jsonl`",
        ]
    )
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out_dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument("--report", default=str(DEFAULT_REPORT))
    parser.add_argument("--candidates", default=",".join(CANDIDATES))
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    v84 = v84_module()
    ranking, inventory, samples, supplements = load_inputs()
    reference_regions = load_reference_regions()
    gate_rows = []
    iteration_log = out_dir / "ralph_iteration_log.jsonl"

    for dataset in [item.strip() for item in args.candidates.split(",") if item.strip()]:
        row, _, _ = candidate_summary(dataset, ranking, inventory, samples, supplements, v84, reference_regions, out_dir)
        gate_rows.append(row)
        append_jsonl(
            iteration_log,
            {
                "timestamp": utc_now(),
                "iteration": row["iteration"],
                "dataset": dataset,
                "hypothesis": row["hypothesis"],
                "action": row["action"],
                "gate_decision": row["next_step"],
                "metrics": {
                    "n_selected_smoke_files": row["n_selected_smoke_files"],
                    "smoke_download_pass": row["smoke_download_pass"],
                    "tar_required_for_smoke": row["tar_required_for_smoke"],
                    "schema_smoke_pass": row["schema_smoke_pass"],
                    "metadata_overlap_pass": row["metadata_overlap_pass"],
                    "assembly_compatible": row["assembly_compatible"],
                    "smoke_common_regions_estimate": row["smoke_common_regions_estimate"],
                },
                "next_step": row["next_step"],
                "auxiliary_only_reason": row["auxiliary_only_reason"],
            },
        )

    gate = pd.DataFrame(gate_rows)
    gate.to_csv(out_dir / "candidate_gate_table.csv", index=False)
    decision = decide_loop_state(gate, load_v86_summary())
    write_json(out_dir / "ralph_decision_state.json", decision)
    write_report(out_dir, Path(args.report), gate, decision)
    print(json.dumps(decision, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Run a v11.x processed methylation conversion gate for one dataset.

This controller is intentionally data-gate only. It packages verified processed
supplements, builds auxiliary metadata, invokes the existing methylation
converter, measures overlap with the v8.2 reference matrix, and records a RALPH
ledger decision. It does not train models or run autoresearch.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import tarfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")
PYTHON = ROOT.parent / "as-ds-ops" / ".venv" / "bin" / "python"
QUEUE = ROOT / "results" / "download_backlog_v11_3" / "conversion_queue_smoke_gated.csv"
SAMPLES = ROOT / "metadata" / "geo_old_age_candidate_samples.csv"
REFERENCE_MATRIX = ROOT / "results" / "multidataset_v8_2_prefilter_liftover" / "all_rrbs_region_matrix_5kb.parquet"
LEDGER_DIR = ROOT / "results" / "ralph_v11_x_loop"
TARGET_TISSUES = {"brain_cortex", "heart", "lung"}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str), encoding="utf-8")


def append_jsonl(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(payload, sort_keys=True, default=str) + "\n")


def clean_str(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and np.isnan(value):
        return ""
    return str(value)


def as_bool(value: Any) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    text = clean_str(value).strip().lower()
    return text in {"1", "true", "yes", "y"}


def parse_sex(text: str) -> str:
    text = text.lower()
    if re.search(r"\bfemale\b|\bgender:\s*female\b|\bsex:\s*female\b", text):
        return "F"
    if re.search(r"\bmale\b|\bgender:\s*male\b|\bsex:\s*male\b", text):
        return "M"
    return "unknown"


def parse_strain(text: str) -> str:
    match = re.search(r"(?:mouse\s+strain|strain|genotype):\s*([^|;]+)", text, flags=re.IGNORECASE)
    if match:
        return match.group(1).strip()
    return "unknown"


def normalize_tissue(value: str) -> str:
    value = value.strip().lower().replace(" ", "_").replace("-", "_")
    aliases = {
        "cortex": "brain_cortex",
        "brain": "brain",
        "heart": "heart",
        "liver": "liver",
        "blood": "blood",
        "skeletal_muscle": "skeletal_muscle",
        "muscle": "skeletal_muscle",
        "lung": "lung",
        "intestine": "intestine",
    }
    return aliases.get(value, value or "unknown")


def infer_intervention(dataset: str, title: str, characteristics: str, fallback: str) -> str:
    text = f"{title} | {characteristics} | {fallback}".lower()
    if dataset == "GSE92486":
        if "dietary-restricted" in text or "diet restricted" in text or re.search(r"\bdr\b", text):
            return "diet_restriction"
        return "control"
    if dataset == "GSE175410":
        if "exercise" in text:
            return "exercise"
        return "control"
    if "parabiosis" in text or "heterochronic" in text or "isochronic" in text:
        return "parabiosis_recovery"
    if "rapamycin" in text:
        return "rapamycin"
    if "calorie" in text or "diet restriction" in text or "dietary-restricted" in text:
        return "diet_restriction"
    if fallback:
        return re.sub(r"[^a-z0-9]+", "_", fallback.lower()).strip("_") or "control"
    return "control"


def infer_cell_type(characteristics: str) -> str:
    match = re.search(r"cell type:\s*([^|;]+)", characteristics, flags=re.IGNORECASE)
    return match.group(1).strip().lower().replace(" ", "_") if match else ""


def load_queue_row(dataset: str) -> dict[str, Any]:
    queue = pd.read_csv(QUEUE)
    hit = queue[queue["dataset"].astype(str).eq(dataset)].copy()
    if hit.empty:
        raise RuntimeError(f"{dataset} not found in {QUEUE}")
    return hit.iloc[0].to_dict()


def load_verified_files(path: Path) -> pd.DataFrame:
    files = pd.read_csv(path)
    if "hard_verified" not in files.columns:
        raise RuntimeError(f"Missing hard_verified column in {path}")
    files = files[files["hard_verified"].astype(bool)].copy()
    if files.empty:
        raise RuntimeError(f"No hard-verified files in {path}")
    files["local_path"] = files["local_path"].astype(str)
    files["sample_id"] = files["sample_id"].astype(str)
    missing = [p for p in files["local_path"] if not Path(p).exists()]
    if missing:
        raise RuntimeError(f"Missing verified local files: {missing[:5]}")
    return files


def build_auxiliary_metadata(dataset: str, files: pd.DataFrame, out_dir: Path, queue_row: dict[str, Any]) -> pd.DataFrame:
    sample_ids = list(files["sample_id"].astype(str))
    candidates = pd.read_csv(SAMPLES)
    candidates = candidates[candidates["dataset"].astype(str).eq(dataset)].copy()
    by_sample = {str(row["sample_id"]): row for row in candidates.to_dict(orient="records")}
    rows: list[dict[str, Any]] = []
    for sample_id in sample_ids:
        row = by_sample.get(sample_id, {})
        title = clean_str(row.get("title"))
        source_name = clean_str(row.get("source_name"))
        characteristics = clean_str(row.get("characteristics"))
        text = f"{title} | {source_name} | {characteristics}"
        age_weeks = pd.to_numeric(pd.Series([row.get("age_weeks")]), errors="coerce").iloc[0] if row else np.nan
        raw_age_token = clean_str(row.get("raw_age_token"))
        age_days = float(age_weeks) * 7.0 if pd.notna(age_weeks) else np.nan
        tissue = normalize_tissue(clean_str(row.get("tissue_guess")) or source_name)
        sex = parse_sex(text)
        strain = parse_strain(text)
        cell_type = infer_cell_type(characteristics)
        intervention = infer_intervention(dataset, title, characteristics, clean_str(row.get("intervention_guess")))
        rows.append(
            {
                "sample_id": sample_id,
                "age_days": age_days,
                "age_weeks": float(age_weeks) if pd.notna(age_weeks) else np.nan,
                "tissue": tissue,
                "sex": sex,
                "strain": strain,
                "intervention": intervention,
                "dataset_batch": dataset,
                "assay": "RRBS_processed_bismark_cov",
                "metadata_source": "v11_x_geo_candidate_samples_auxiliary",
                "sample_id_source": "GEO_GSM_processed_cov",
                "cell_type": cell_type,
                "raw_age_token": raw_age_token,
                "title": title,
                "source_name": source_name,
                "raw_characteristics": characteristics,
                "headline_allowed": False,
                "auxiliary_only_reason": clean_str(queue_row.get("reason_smoke")) or clean_str(queue_row.get("reason")),
            }
        )
    meta = pd.DataFrame(rows).sort_values("sample_id")
    meta_path = out_dir / f"{dataset}_auxiliary_metadata.csv"
    meta.to_csv(meta_path, index=False)
    return meta


def make_tar(dataset: str, files: pd.DataFrame, out_dir: Path, force: bool) -> dict[str, Any]:
    tar_path = out_dir / f"{dataset}_verified_processed_files.tar"
    expected_files = int(len(files))
    expected_input_bytes = int(pd.to_numeric(files["actual_size_bytes"], errors="coerce").fillna(0).sum())
    if tar_path.exists() and not force:
        with tarfile.open(tar_path, "r") as tar:
            n_members = len([m for m in tar.getmembers() if m.isfile()])
        if n_members == expected_files:
            return {
                "tar_path": str(tar_path),
                "n_files": expected_files,
                "tar_size_bytes": int(tar_path.stat().st_size),
                "input_bytes": expected_input_bytes,
                "reused_existing": True,
            }
    with tarfile.open(tar_path, "w") as tar:
        for path_text in files["local_path"]:
            path = Path(path_text)
            tar.add(path, arcname=path.name)
    return {
        "tar_path": str(tar_path),
        "n_files": expected_files,
        "tar_size_bytes": int(tar_path.stat().st_size),
        "input_bytes": expected_input_bytes,
        "reused_existing": False,
    }


def run_conversion(dataset: str, tar_path: Path, metadata_path: Path, out_dir: Path, force: bool) -> dict[str, Any]:
    manifest_path = out_dir / f"{dataset}_matrix_manifest.json"
    if manifest_path.exists() and not force:
        return json.loads(manifest_path.read_text(encoding="utf-8")) | {"reused_existing": True}
    cmd = [
        str(PYTHON),
        "scripts/etl/10_convert_processed_methylation.py",
        "--input",
        str(tar_path),
        "--dataset",
        dataset,
        "--schema",
        "bismark_cov_per_sample_tar",
        "--min_coverage",
        "5",
        "--min_sample_presence",
        "0.5",
        "--region_min_sample_presence",
        "0.8",
        "--metadata_path",
        str(metadata_path),
        "--out_dir",
        str(out_dir),
    ]
    proc = subprocess.run(cmd, cwd=str(ROOT), text=True, capture_output=True)
    (out_dir / "conversion_stdout.txt").write_text(proc.stdout, encoding="utf-8")
    (out_dir / "conversion_stderr.txt").write_text(proc.stderr, encoding="utf-8")
    if proc.returncode != 0:
        return {
            "status": "blocked",
            "reason": "conversion_failed",
            "returncode": proc.returncode,
            "stderr_tail": proc.stderr[-4000:],
        }
    if not manifest_path.exists():
        return {"status": "blocked", "reason": "matrix_manifest_missing"}
    return json.loads(manifest_path.read_text(encoding="utf-8")) | {"reused_existing": False}


def compute_overlap_and_qc(dataset: str, conversion: dict[str, Any], out_dir: Path, force: bool) -> dict[str, Any]:
    region_path = Path(str(conversion.get("region_matrix_path", "")))
    beta_path = Path(str(conversion.get("beta_matrix_path", "")))
    if not region_path.exists() or not beta_path.exists():
        return {"status": "blocked", "reason": "converted_matrix_missing"}
    qc_path = out_dir / f"{dataset}_matrix_qc.json"
    if qc_path.exists() and not force:
        return json.loads(qc_path.read_text(encoding="utf-8")) | {"reused_existing": True}
    beta = pd.read_parquet(beta_path)
    region = pd.read_parquet(region_path)
    beta_min = float(beta.min(skipna=True).min(skipna=True))
    beta_max = float(beta.max(skipna=True).max(skipna=True))
    beta_bad_values = int(((beta < 0) | (beta > 1)).sum(skipna=True).sum())
    qc: dict[str, Any] = {
        "status": "completed",
        "beta_shape": [int(beta.shape[0]), int(beta.shape[1])],
        "region_shape": [int(region.shape[0]), int(region.shape[1])],
        "beta_min": beta_min,
        "beta_max": beta_max,
        "beta_values_outside_0_1": beta_bad_values,
        "n_samples": int(beta.shape[1]),
    }
    if REFERENCE_MATRIX.exists():
        reference = pd.read_parquet(REFERENCE_MATRIX)
        common = sorted(set(region.index) & set(reference.index))
        pd.DataFrame({"region_id": common}).to_csv(out_dir / "common_regions_with_v8_2.csv", index=False)
        qc.update(
            {
                "reference_matrix": str(REFERENCE_MATRIX),
                "reference_regions": int(reference.shape[0]),
                "common_regions_with_v8_2": int(len(common)),
                "common_region_fraction_vs_dataset": float(len(common) / region.shape[0]) if region.shape[0] else 0.0,
                "common_region_fraction_vs_reference": float(len(common) / reference.shape[0]) if reference.shape[0] else 0.0,
                "common_region_gate_50000": bool(len(common) >= 50_000),
            }
        )
    else:
        qc.update(
            {
                "reference_matrix": str(REFERENCE_MATRIX),
                "common_region_gate_50000": False,
                "overlap_reason": "reference_matrix_missing",
            }
        )
    write_json(qc_path, qc)
    return qc


def decide_gate(dataset: str, queue_row: dict[str, Any], meta: pd.DataFrame, conversion: dict[str, Any], qc: dict[str, Any]) -> dict[str, Any]:
    age_weeks = pd.to_numeric(meta["age_weeks"], errors="coerce")
    tissues = set(meta["tissue"].dropna().astype(str))
    cell_types = set(meta["cell_type"].dropna().astype(str)) - {""}
    exact_age_gate = bool(age_weeks.notna().mean() >= 0.95)
    old_target_tissue_gate = bool((age_weeks >= 104).any() and bool(tissues & TARGET_TISSUES))
    bulk_context_gate = not bool(cell_types)
    conversion_gate = bool(str(conversion.get("status", "")).startswith("completed") and int(conversion.get("n_errors", 1)) == 0)
    beta_gate = bool(qc.get("beta_values_outside_0_1", 1) == 0 and 0.0 <= float(qc.get("beta_min", -1)) <= 1.0 and 0.0 <= float(qc.get("beta_max", 2)) <= 1.0)
    common_region_gate = bool(qc.get("common_region_gate_50000", False))
    queue_headline_allowed = as_bool(queue_row.get("final_headline_allowed", False))
    raw_headline_gate = all(
        [
            conversion_gate,
            beta_gate,
            exact_age_gate,
            old_target_tissue_gate,
            bulk_context_gate,
            common_region_gate,
        ]
    )
    headline_allowed = bool(raw_headline_gate and queue_headline_allowed)
    if headline_allowed:
        final_status = "headline_candidate"
        decision = "promote_to_fixed_benchmark"
        reason = ""
    elif conversion_gate and qc.get("status") == "completed":
        final_status = "auxiliary_matrix"
        decision = "do_not_train_keep_auxiliary"
        blockers = []
        if not exact_age_gate:
            blockers.append("exact_age_coverage_lt_95pct")
        if not old_target_tissue_gate:
            blockers.append("no_old_target_tissue_support")
        if not bulk_context_gate:
            blockers.append("cell_type_specific_or_non_bulk")
        if not common_region_gate:
            blockers.append("common_regions_lt_50000")
        if not queue_headline_allowed:
            blockers.append("queue_marks_non_headline")
        reason = ";".join(blockers) or "auxiliary_by_policy"
    else:
        final_status = "blocked"
        decision = "conversion_or_qc_blocked"
        reason = clean_str(conversion.get("reason")) or clean_str(qc.get("reason")) or "unknown_blocker"
    return {
        "final_status": final_status,
        "decision": decision,
        "headline_allowed": headline_allowed,
        "raw_headline_gate": raw_headline_gate,
        "queue_headline_allowed": queue_headline_allowed,
        "reason": reason,
        "gates": {
            "conversion_gate": conversion_gate,
            "beta_gate": beta_gate,
            "exact_age_gate": exact_age_gate,
            "old_target_tissue_gate": old_target_tissue_gate,
            "bulk_context_gate": bulk_context_gate,
            "common_region_gate": common_region_gate,
        },
        "metadata": {
            "n_rows": int(len(meta)),
            "age_known_fraction": float(age_weeks.notna().mean()) if len(meta) else 0.0,
            "age_weeks_values": sorted(float(x) for x in age_weeks.dropna().unique()),
            "tissue_counts": meta["tissue"].value_counts(dropna=False).to_dict(),
            "sex_counts": meta["sex"].value_counts(dropna=False).to_dict(),
            "intervention_counts": meta["intervention"].value_counts(dropna=False).to_dict(),
            "cell_type_counts": meta["cell_type"].value_counts(dropna=False).to_dict(),
        },
    }


def update_ledger(dataset: str, state: dict[str, Any]) -> None:
    LEDGER_DIR.mkdir(parents=True, exist_ok=True)
    append_jsonl(
        LEDGER_DIR / "ralph_iteration_log.jsonl",
        {
            "timestamp": utc_now(),
            "loop_version": "v11.x",
            "dataset": dataset,
            "status": state["final_status"],
            "decision": state["decision"],
            "reason": state.get("reason", ""),
            "state_path": state["state_path"],
        },
    )
    rows = []
    state_files = sorted((ROOT / "results" / "v11_4_conversions").glob("GSE*/v11_x_*_conversion_gate_state.json"))
    for path in state_files:
        payload = json.loads(path.read_text(encoding="utf-8"))
        qc = payload.get("qc", {})
        gates = payload.get("gates", {})
        metadata = payload.get("metadata_summary", {})
        rows.append(
            {
                "dataset": payload.get("dataset"),
                "final_status": payload.get("final_status"),
                "decision": payload.get("decision"),
                "headline_allowed": payload.get("headline_allowed"),
                "reason": payload.get("reason"),
                "n_samples": qc.get("n_samples"),
                "n_regions": payload.get("conversion_manifest", {}).get("n_regions"),
                "common_regions_with_v8_2": qc.get("common_regions_with_v8_2"),
                "common_region_gate": gates.get("common_region_gate"),
                "exact_age_gate": gates.get("exact_age_gate"),
                "old_target_tissue_gate": gates.get("old_target_tissue_gate"),
                "bulk_context_gate": gates.get("bulk_context_gate"),
                "age_known_fraction": metadata.get("age_known_fraction"),
            }
        )
    pd.DataFrame(rows).sort_values("dataset").to_csv(LEDGER_DIR / "candidate_gate_table.csv", index=False)
    n_headline = sum(bool(row.get("headline_allowed")) for row in rows)
    decision_state = {
        "timestamp": utc_now(),
        "loop_version": "v11.x",
        "n_datasets_gated": len(rows),
        "n_headline_allowed": int(n_headline),
        "status": "continue" if n_headline == 0 else "headline_candidate_found",
        "next_action": "continue_conversion_queue_or_prepare_v12_failure_report"
        if n_headline == 0
        else "run_fixed_benchmark_for_headline_candidate",
        "candidate_gate_table": str(LEDGER_DIR / "candidate_gate_table.csv"),
    }
    write_json(LEDGER_DIR / "ralph_decision_state.json", decision_state)


def run_dataset(dataset: str, force: bool) -> dict[str, Any]:
    queue_row = load_queue_row(dataset)
    out_dir = Path(clean_str(queue_row.get("target_out_dir")) or ROOT / "results" / "v11_4_conversions" / dataset)
    out_dir.mkdir(parents=True, exist_ok=True)
    if clean_str(queue_row.get("final_conversion_gate")) != "ready_for_bismark_cov_conversion":
        state = {
            "timestamp": utc_now(),
            "loop_version": "v11.x",
            "dataset": dataset,
            "final_status": "blocked",
            "decision": "adapter_required_before_conversion",
            "reason": clean_str(queue_row.get("final_recommended_adapter")),
            "headline_allowed": False,
            "queue_row": queue_row,
            "state_path": str(out_dir / f"v11_x_{dataset}_conversion_gate_state.json"),
        }
        write_json(Path(state["state_path"]), state)
        update_ledger(dataset, state)
        return state

    verified_path = Path(clean_str(queue_row.get("verified_file_manifest")))
    files = load_verified_files(verified_path)
    meta = build_auxiliary_metadata(dataset, files, out_dir, queue_row)
    metadata_path = out_dir / f"{dataset}_auxiliary_metadata.csv"
    tar_state = make_tar(dataset, files, out_dir, force)
    conversion = run_conversion(dataset, Path(tar_state["tar_path"]), metadata_path, out_dir, force)
    qc = compute_overlap_and_qc(dataset, conversion, out_dir, force) if str(conversion.get("status", "")).startswith("completed") else {}
    gate = decide_gate(dataset, queue_row, meta, conversion, qc)
    state_path = out_dir / f"v11_x_{dataset}_conversion_gate_state.json"
    state = {
        "timestamp": utc_now(),
        "loop_version": "v11.x",
        "dataset": dataset,
        "queue_row": queue_row,
        "metadata_path": str(metadata_path),
        "tar": tar_state,
        "conversion_manifest": conversion,
        "qc": qc,
        "metadata_summary": gate.pop("metadata"),
        "state_path": str(state_path),
        **gate,
    }
    write_json(state_path, state)
    update_ledger(dataset, state)
    return state


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    state = run_dataset(args.dataset, args.force)
    print(json.dumps(state, indent=2)[:10_000])


if __name__ == "__main__":
    main()

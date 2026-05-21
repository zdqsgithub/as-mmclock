#!/usr/bin/env python3
"""Prepare gate-passed v24 raw-derived matrices for local biological-signal ML.

This script does not touch FASTQ/BAM/COV files. It combines only datasets that
passed the raw matrix gate and writes a complete-region matrix plus an
exploratory target registry for small-n local ML screens.
"""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path("/home/zdq-as/mouse_methyl_work")
DEFAULT_GSE121141_MATRIX = Path(
    "/data/mouse_methyl/processed_v21_raw_etl/GSE121141/GSE121141_raw_region_matrix_5kb.parquet"
)
DEFAULT_GSE80672_MATRIX = Path(
    "/data/mouse_methyl/processed_v21_raw_etl/GSE80672/GSE80672_raw_region_matrix_5kb.parquet"
)
DEFAULT_METADATA = ROOT / "metadata" / "model_sample_metadata_v8.csv"
DEFAULT_GATE = ROOT / "results" / "ralph_v21_raid_raw_clock" / "matrix_gate_table.csv"
DEFAULT_OUT_DIR = ROOT / "results" / "v24_raw_signal_screen"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str), encoding="utf-8")


def read_gate(path: Path) -> set[str]:
    if not path.exists():
        raise SystemExit(f"Matrix gate table missing: {path}")
    gate = pd.read_csv(path)
    required = {"dataset", "matrix_gate_passed"}
    missing = required - set(gate.columns)
    if missing:
        raise SystemExit(f"Matrix gate table missing columns: {sorted(missing)}")
    passed = gate[gate["matrix_gate_passed"].astype(str).str.lower().isin({"true", "1"})]
    return set(passed["dataset"].astype(str))


def target_registry() -> pd.DataFrame:
    rows = [
        {
            "target_id": "age_weeks_raw_v24",
            "target_kind": "regression",
            "label_column": "age_days",
            "positive_label": "",
            "negative_label": "",
            "dataset_scope": "all",
            "include_values": "",
            "min_total": 24,
            "min_class_count": 0,
            "min_groups": 2,
            "recommended_local_model": "ridge",
            "primary_metric": "mae_weeks",
            "priority": "primary_exploratory",
            "biological_question": "Chronological age signal in the gate-passed raw-derived methylation matrix.",
        },
        {
            "target_id": "cr_vs_control_raw_v24",
            "target_kind": "binary",
            "label_column": "intervention",
            "positive_label": "CR",
            "negative_label": "control",
            "dataset_scope": "GSE80672",
            "include_values": "CR;control",
            "min_total": 20,
            "min_class_count": 8,
            "min_groups": 1,
            "recommended_local_model": "logistic",
            "primary_metric": "roc_auc",
            "priority": "primary_aux_exploratory",
            "biological_question": "Calorie-restriction-associated methylation signal in raw-derived GSE80672 blood samples.",
        },
        {
            "target_id": "dataset_batch_raw_v24",
            "target_kind": "binary",
            "label_column": "dataset_batch",
            "positive_label": "GSE80672",
            "negative_label": "GSE121141",
            "dataset_scope": "all",
            "include_values": "GSE121141;GSE80672",
            "min_total": 20,
            "min_class_count": 5,
            "min_groups": 1,
            "recommended_local_model": "logistic",
            "primary_metric": "balanced_accuracy",
            "priority": "confounding_probe",
            "biological_question": "Dataset/batch separability probe used to flag confounding, not a biological endpoint.",
        },
        {
            "target_id": "tissue_blood_vs_brain_raw_v24",
            "target_kind": "binary",
            "label_column": "tissue",
            "positive_label": "blood",
            "negative_label": "brain_cortex",
            "dataset_scope": "all",
            "include_values": "blood;brain_cortex",
            "min_total": 20,
            "min_class_count": 5,
            "min_groups": 1,
            "recommended_local_model": "logistic",
            "primary_metric": "balanced_accuracy",
            "priority": "confounding_probe",
            "biological_question": "Tissue separability probe used to flag age/CR model confounding.",
        },
    ]
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gse121141-matrix", type=Path, default=DEFAULT_GSE121141_MATRIX)
    parser.add_argument("--gse80672-matrix", type=Path, default=DEFAULT_GSE80672_MATRIX)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--matrix-gate-table", type=Path, default=DEFAULT_GATE)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()

    passed = read_gate(args.matrix_gate_table)
    needed = {"GSE121141", "GSE80672"}
    if not needed.issubset(passed):
        raise SystemExit(f"Required datasets have not passed matrix gate: {sorted(needed - passed)}")
    if not args.gse121141_matrix.exists():
        raise SystemExit(f"Missing matrix: {args.gse121141_matrix}")
    if not args.gse80672_matrix.exists():
        raise SystemExit(f"Missing matrix: {args.gse80672_matrix}")
    if not args.metadata.exists():
        raise SystemExit(f"Missing metadata: {args.metadata}")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    gse121141 = pd.read_parquet(args.gse121141_matrix)
    gse80672 = pd.read_parquet(args.gse80672_matrix)
    common_regions = gse121141.index.intersection(gse80672.index)
    combined = pd.concat([gse121141.loc[common_regions], gse80672.loc[common_regions]], axis=1)
    complete = combined[combined.notna().all(axis=1)].copy()
    if complete.shape[0] < 50_000:
        raise SystemExit(f"Complete common regions below gate: {complete.shape[0]} < 50000")

    metadata = pd.read_csv(args.metadata)
    samples = [sample for sample in complete.columns.astype(str) if sample in set(metadata["sample_id"].astype(str))]
    complete = complete[samples]
    metadata = metadata[metadata["sample_id"].astype(str).isin(samples)].copy()
    metadata["raw_v24_signal_screen"] = True
    metadata = metadata.sort_values(["dataset_batch", "age_weeks", "sample_id"]).reset_index(drop=True)
    complete = complete[metadata["sample_id"].astype(str).tolist()]

    matrix_path = args.out_dir / "raw_v24_gse121141_gse80672_common_complete_region_matrix_5kb.parquet"
    metadata_path = args.out_dir / "raw_v24_signal_metadata.csv"
    registry_path = args.out_dir / "raw_v24_biological_signal_targets.csv"
    complete.to_parquet(matrix_path)
    metadata.to_csv(metadata_path, index=False)
    target_registry().to_csv(registry_path, index=False)

    summary = {
        "timestamp": utc_now(),
        "status": "prepared",
        "matrix": str(matrix_path),
        "metadata": str(metadata_path),
        "target_registry": str(registry_path),
        "datasets": metadata["dataset_batch"].value_counts().to_dict(),
        "interventions": metadata["intervention"].value_counts(dropna=False).to_dict(),
        "tissues": metadata["tissue"].value_counts(dropna=False).to_dict(),
        "sexes": metadata["sex"].value_counts(dropna=False).to_dict(),
        "n_samples": int(metadata.shape[0]),
        "n_regions_intersection": int(len(common_regions)),
        "n_regions_complete_all_samples": int(complete.shape[0]),
        "matrix_gate_table": str(args.matrix_gate_table),
        "interpretation_guardrail": "small-n exploratory raw-derived screen; target-associated signal only",
    }
    write_json(args.out_dir / "raw_v24_signal_input_summary.json", summary)
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Plan a balanced raw ETL expansion queue after v21 pilot gates pass.

The queue remains read-only with respect to FASTQ inputs. It is intended for
`scripts/etl/24_raw_etl_watchdog.py` and includes already-completed pilot
samples so dataset-level matrices can be rebuilt with old plus new samples.
"""
from __future__ import annotations

import argparse
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path("/home/zdq-as/mouse_methyl_work")
DEFAULT_PAIRING = ROOT / "results" / "ralph_v21_raid_raw_clock" / "pairing_qc_manifest.csv"
DEFAULT_METADATA = ROOT / "metadata" / "model_sample_metadata_v8.csv"
DEFAULT_ETL_ROOT = Path("/data/mouse_methyl/processed_v21_raw_etl")
DEFAULT_OUT = ROOT / "results" / "ralph_v21_raid_raw_clock" / "raw_etl_queue_v24_balanced_expansion.csv"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def infer_mate(path: str) -> int | None:
    name = Path(path).name
    for pattern in (r"(?:^|[_\-.])R?([12])(?:[_\-.]|$)", r"_R([12])_", r"_([12])\.f(?:ast)?q(?:\.gz)?$"):
        match = re.search(pattern, name, flags=re.IGNORECASE)
        if match:
            return int(match.group(1))
    return None


def load_metadata(path: Path) -> pd.DataFrame:
    meta = pd.read_csv(path)
    if "dataset" not in meta.columns:
        meta["dataset"] = meta["dataset_batch"].astype(str).str.extract(r"^(GSE\d+)")[0].fillna(meta["dataset_batch"].astype(str))
    return meta


def sample_level_pairing(pairing: pd.DataFrame) -> pd.DataFrame:
    complete = pairing[pairing["pairing_status"].eq("complete_pair")].copy()
    rows: list[dict[str, Any]] = []
    for (dataset, sample_id), group in complete.groupby(["dataset", "sample_id"], dropna=False):
        r1_paths: list[str] = []
        r2_paths: list[str] = []
        run_accessions: list[str] = []
        for _, row in group.sort_values("run_accession").iterrows():
            paths = [item for item in str(row.get("local_paths") or "").split(";") if item]
            r1 = sorted(path for path in paths if infer_mate(path) == 1)
            r2 = sorted(path for path in paths if infer_mate(path) == 2)
            if not r1 or not r2:
                continue
            r1_paths.extend(r1)
            r2_paths.extend(r2)
            run_accessions.append(str(row["run_accession"]))
        if not r1_paths or len(r1_paths) != len(r2_paths):
            continue
        rows.append(
            {
                "dataset": dataset,
                "sample_id": sample_id,
                "run_accessions": ";".join(run_accessions),
                "r1_paths": ";".join(r1_paths),
                "r2_paths": ";".join(r2_paths),
                "n_runs": len(run_accessions),
                "n_fastq_files": len(r1_paths) + len(r2_paths),
                "fastq_gib": float(group["total_size_gib"].sum()),
            }
        )
    return pd.DataFrame(rows)


def existing_completed_samples(etl_root: Path) -> set[tuple[str, str]]:
    completed: set[tuple[str, str]] = set()
    if not etl_root.exists():
        return completed
    for path in etl_root.glob("GSE*/GSM*/*_raw_region_beta.parquet"):
        completed.add((path.parents[1].name, path.parent.name))
    return completed


def round_robin_by_group(df: pd.DataFrame, group_cols: list[str], limit: int, required: set[str]) -> list[str]:
    if df.empty or limit <= 0:
        return []
    selected: list[str] = []
    seen: set[str] = set()

    required_rows = df[df["sample_id"].isin(required)].sort_values(["sample_id"])
    for sample_id in required_rows["sample_id"].tolist():
        if sample_id not in seen:
            selected.append(sample_id)
            seen.add(sample_id)
    groups = []
    for key, group in df.sort_values(group_cols + ["fastq_gib", "sample_id"]).groupby(group_cols, dropna=False):
        groups.append((key, group[~group["sample_id"].isin(seen)]["sample_id"].tolist()))
    remaining_slots = max(0, limit - len(selected))
    if remaining_slots and len(groups) > remaining_slots:
        if remaining_slots == 1:
            indexes = [len(groups) // 2]
        else:
            indexes = sorted({round(i * (len(groups) - 1) / (remaining_slots - 1)) for i in range(remaining_slots)})
        groups = [groups[idx] for idx in indexes]
    while len(selected) < limit:
        progressed = False
        for _, sample_ids in groups:
            while sample_ids and sample_ids[0] in seen:
                sample_ids.pop(0)
            if sample_ids and len(selected) < limit:
                sample_id = sample_ids.pop(0)
                selected.append(sample_id)
                seen.add(sample_id)
                progressed = True
        if not progressed:
            break
    return selected[:limit]


def build_queue(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame]:
    pairing = pd.read_csv(args.pairing)
    meta = load_metadata(args.metadata)
    samples = sample_level_pairing(pairing)
    merged = samples.merge(
        meta[
            [
                "sample_id",
                "dataset",
                "age_weeks",
                "tissue",
                "sex",
                "intervention",
                "title",
                "condition_family",
            ]
        ],
        on=["dataset", "sample_id"],
        how="left",
    )
    merged = merged[merged["age_weeks"].notna()].copy()
    existing = existing_completed_samples(args.etl_root)
    existing_by_dataset: dict[str, set[str]] = {}
    for dataset, sample_id in existing:
        existing_by_dataset.setdefault(dataset, set()).add(sample_id)

    selected: list[str] = []
    selection_rows: list[dict[str, Any]] = []

    def add(dataset: str, sample_ids: list[str], reason: str) -> None:
        for rank, sample_id in enumerate(sample_ids, start=1):
            key = f"{dataset}/{sample_id}"
            if key in selected:
                continue
            selected.append(key)
            row = merged[(merged["dataset"].eq(dataset)) & (merged["sample_id"].eq(sample_id))].iloc[0].to_dict()
            row["selection_reason"] = reason
            row["selection_rank_within_reason"] = rank
            row["already_completed_raw_region_beta"] = (dataset, sample_id) in existing
            selection_rows.append(row)

    gse121141 = merged[merged["dataset"].eq("GSE121141")].sort_values(["age_weeks", "sample_id"])
    add("GSE121141", gse121141["sample_id"].tolist(), "all_complete_paired_metadata_mapped")

    gse80672 = merged[merged["dataset"].eq("GSE80672")].copy()
    cr = gse80672[gse80672["intervention"].eq("CR")]
    control = gse80672[gse80672["intervention"].eq("control")]
    add(
        "GSE80672",
        round_robin_by_group(cr, ["age_weeks"], args.gse80672_cr_samples, existing_by_dataset.get("GSE80672", set())),
        "cr_age_balanced",
    )
    add(
        "GSE80672",
        round_robin_by_group(control, ["age_weeks"], args.gse80672_control_samples, existing_by_dataset.get("GSE80672", set())),
        "control_age_balanced",
    )

    gse93957 = merged[merged["dataset"].eq("GSE93957")].copy()
    add(
        "GSE93957",
        round_robin_by_group(gse93957, ["age_weeks", "tissue"], args.gse93957_samples, existing_by_dataset.get("GSE93957", set())),
        "age_tissue_balanced",
    )

    selected_df = pd.DataFrame(selection_rows)
    if selected_df.empty:
        return pd.DataFrame(), selected_df
    selected_df["pilot_priority"] = selected_df["dataset"].map({"GSE121141": 0, "GSE80672": 1, "GSE93957": 2}).fillna(99).astype(int)
    selected_df = selected_df.sort_values(
        ["pilot_priority", "already_completed_raw_region_beta", "selection_reason", "age_weeks", "sample_id"],
        ascending=[True, False, True, True, True],
    ).reset_index(drop=True)
    queue_cols = [
        "dataset",
        "sample_id",
        "run_accessions",
        "r1_paths",
        "r2_paths",
        "n_runs",
        "n_fastq_files",
        "pilot_priority",
    ]
    queue = selected_df[queue_cols].copy()
    queue["queue_status"] = "v24_balanced_expansion"
    queue["bismark_authorized"] = False
    queue["training_authorized"] = False
    return queue, selected_df


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairing", type=Path, default=DEFAULT_PAIRING)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--etl-root", type=Path, default=DEFAULT_ETL_ROOT)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--gse80672-cr-samples", type=int, default=12)
    parser.add_argument("--gse80672-control-samples", type=int, default=12)
    parser.add_argument("--gse93957-samples", type=int, default=12)
    args = parser.parse_args()

    queue, selected = build_queue(args)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    queue.to_csv(args.out, index=False)
    selected_path = args.out.with_name(args.out.stem + "_selection_audit.csv")
    selected.to_csv(selected_path, index=False)
    summary = {
        "timestamp": utc_now(),
        "queue_path": str(args.out),
        "selection_audit_path": str(selected_path),
        "n_queue_rows": int(len(queue)),
        "datasets": queue["dataset"].value_counts().to_dict() if not queue.empty else {},
        "already_completed_rows": int(selected.get("already_completed_raw_region_beta", pd.Series(dtype=bool)).fillna(False).sum()) if not selected.empty else 0,
        "new_rows": int((~selected.get("already_completed_raw_region_beta", pd.Series(dtype=bool)).fillna(False)).sum()) if not selected.empty else 0,
        "policy": "read_only_fastq_balanced_expansion_after_v21_gate",
    }
    summary_path = args.out.with_name(args.out.stem + "_summary.json")
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

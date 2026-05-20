#!/usr/bin/env python3
"""Build the v11.3 post-download conversion queue.

The watchdog status file can be stale if the final monitor poll happens before
the last worker finishes. This script therefore performs hard size verification
against the manifest and only queues fully verified datasets for later schema
smoke/conversion.
"""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")
DEFAULT_IN = ROOT / "results" / "download_backlog_v11_3" / "download_manifest.csv"
DEFAULT_OUT = ROOT / "results" / "download_backlog_v11_3"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def hard_verify(row: dict) -> dict:
    path = Path(str(row["local_path"]))
    expected = int(row["expected_size_bytes"])
    exists = path.exists()
    actual = path.stat().st_size if exists else 0
    ok = exists and actual == expected
    status = "verified" if ok else "missing"
    if exists and actual != expected:
        status = "size_mismatch"
    return {
        **row,
        "hard_verify_status": status,
        "hard_verified": ok,
        "actual_size_bytes": actual,
        "size_delta_bytes": actual - expected,
    }


def conversion_schema(dataset: str, file_types: set[str], schema_guess: str) -> tuple[str, str]:
    if file_types == {"COV"} and schema_guess == "bismark_cov_per_sample_tar":
        return "bismark_cov_per_sample_tar", "ready_for_bismark_cov_tar"
    if dataset in {"GSE92486", "GSE129712"} and "TXT" in file_types:
        return "txt_cov_like_pending_smoke", "requires_txt_adapter_smoke"
    if "BEDGRAPH" in file_types:
        return "bedgraph_pending_smoke", "requires_bedgraph_adapter_smoke"
    return schema_guess or "unknown", "requires_schema_smoke"


def build_queue(verified: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for dataset, ds in verified.groupby("dataset", sort=True):
        file_types = set(ds["file_type"].astype(str))
        schema_guess = str(ds["schema_guess"].dropna().iloc[0]) if ds["schema_guess"].notna().any() else ""
        target_schema, gate_status = conversion_schema(dataset, file_types, schema_guess)
        target_dir = ROOT / "results" / "v11_4_conversions" / dataset
        target_tar = target_dir / f"{dataset}_verified_processed_files.tar"
        rows.append(
            {
                "dataset": dataset,
                "n_files": int(len(ds)),
                "total_size_bytes": int(ds["expected_size_bytes"].sum()),
                "total_size_gib": round(float(ds["expected_size_bytes"].sum()) / 1024**3, 6),
                "file_types": ";".join(sorted(file_types)),
                "source_schema_guess": schema_guess,
                "conversion_schema": target_schema,
                "gate_status": gate_status,
                "headline_allowed": bool(ds["headline_allowed"].fillna(False).astype(bool).any()),
                "tier": ";".join(sorted(set(ds["tier"].astype(str)))),
                "reason": " | ".join(sorted(set(ds["reason"].astype(str))))[:1000],
                "verified_file_manifest": str(target_dir / f"{dataset}_verified_files.csv"),
                "target_tar_path": str(target_tar),
                "target_out_dir": str(target_dir),
                "conversion_allowed_after_smoke": gate_status == "ready_for_bismark_cov_tar",
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=DEFAULT_IN)
    parser.add_argument("--out_dir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(args.manifest)
    verified_rows = [hard_verify(row) for row in manifest.to_dict(orient="records")]
    completion = pd.DataFrame(verified_rows)
    completion_path = args.out_dir / "download_completion_summary.csv"
    completion.to_csv(completion_path, index=False)

    verified = completion[completion["hard_verified"]].copy()
    failed = completion[~completion["hard_verified"]].copy()
    dataset_counts = (
        completion.groupby(["dataset", "hard_verify_status"])
        .size()
        .unstack(fill_value=0)
        .reset_index()
    )
    dataset_counts.to_csv(args.out_dir / "download_completion_by_dataset.csv", index=False)

    queue = build_queue(verified)
    queue_path = args.out_dir / "conversion_queue.csv"
    queue.to_csv(queue_path, index=False)

    for dataset, ds in verified.groupby("dataset", sort=True):
        target_dir = ROOT / "results" / "v11_4_conversions" / dataset
        target_dir.mkdir(parents=True, exist_ok=True)
        ds.sort_values("supplement_name").to_csv(target_dir / f"{dataset}_verified_files.csv", index=False)

    summary = {
        "timestamp": utc_now(),
        "status": "completed" if failed.empty else "blocked_incomplete_downloads",
        "manifest_path": str(args.manifest),
        "completion_summary_path": str(completion_path),
        "completion_by_dataset_path": str(args.out_dir / "download_completion_by_dataset.csv"),
        "conversion_queue_path": str(queue_path),
        "n_manifest_rows": int(len(manifest)),
        "n_hard_verified": int(len(verified)),
        "n_failed_or_missing": int(len(failed)),
        "verified_bytes": int(verified["expected_size_bytes"].sum()),
        "failed_or_missing_rows": failed[["dataset", "supplement_name", "hard_verify_status", "actual_size_bytes", "expected_size_bytes"]].to_dict(orient="records")[:20],
        "datasets_queued": sorted(queue["dataset"].tolist()) if not queue.empty else [],
    }
    (args.out_dir / "download_completion_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Run the v14 GSE83947 processed-supplement adapter smoke.

This is a narrow 2-3 file processed TXT/CX_report smoke. It downloads only the
sample-level GEO processed supplement files listed in the v14 pilot manifest,
records sha256/size, parses a bounded number of rows, and estimates overlap
with the v8.2 5kb reference regions. It does not download FASTQ, run Bismark,
train models, or start autoresearch.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")
V14_DIR = ROOT / "results" / "ralph_v14_public_data_rescue"
OUT_DIR = V14_DIR / "gse83947_processed_adapter_smoke"


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


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", default="GSE83947")
    parser.add_argument("--pilot-manifest", default=str(V14_DIR / "pilot_run_manifest.csv"))
    parser.add_argument("--out-dir", default=str(OUT_DIR))
    parser.add_argument("--max-files", type=int, default=3)
    parser.add_argument("--authorize-download", action="store_true", help="Required to download the 2-3 small processed supplement files.")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    state_path = out_dir / "gse83947_processed_adapter_smoke_state.json"

    if not args.authorize_download:
        state = {
            "timestamp": utc_now(),
            "status": "blocked",
            "reason": "download_not_authorized",
            "training_authorized": False,
            "bismark_authorized": False,
            "autoresearch_authorized": False,
        }
        write_json(state_path, state)
        print(json.dumps(state, indent=2))
        sys.exit(2)

    manifest = pd.read_csv(args.pilot_manifest)
    data = manifest[
        manifest["dataset"].astype(str).eq(args.dataset)
        & manifest.get("processed_supplement_url", pd.Series("", index=manifest.index)).fillna("").astype(str).ne("")
    ].copy()
    data = data.head(args.max_files)
    if data.empty:
        state = {
            "timestamp": utc_now(),
            "status": "blocked",
            "reason": "no_processed_supplement_rows_in_pilot_manifest",
            "training_authorized": False,
            "bismark_authorized": False,
            "autoresearch_authorized": False,
        }
        write_json(state_path, state)
        print(json.dumps(state, indent=2))
        sys.exit(1)

    v8 = import_module(ROOT / "scripts" / "validate" / "run_v8_ralph_loop.py", "run_v8_helpers_for_v14_gse83947")
    reference_regions = v8.load_reference_regions()
    rows = []
    download_dir = out_dir / "downloaded_processed_supplements"
    download_log = out_dir / "download_log.jsonl"

    for _, row in data.iterrows():
        download_row = pd.Series(
            {
                "dataset": row["dataset"],
                "sample_id": row["sample_accession"],
                "supplement_name": row["processed_supplement_name"],
                "supplement_url": row["processed_supplement_url"],
                "file_size_bytes": int(float(row["processed_supplement_bytes"])),
            }
        )
        local_path, payload = v8.download_smoke_file(download_row, download_dir, download_log, timeout=180, retries=3)
        parse_payload = {}
        if local_path is not None:
            parse_payload = v8.parse_smoke_file(local_path, reference_regions)
        rows.append(
            {
                **row.to_dict(),
                **{f"download_{key}": value for key, value in payload.items()},
                **parse_payload,
            }
        )

    smoke = pd.DataFrame(rows)
    smoke_path = out_dir / "gse83947_processed_adapter_smoke_manifest.csv"
    smoke.to_csv(smoke_path, index=False)
    rows_parseable = pd.to_numeric(smoke.get("rows_parseable_primary_autosomes", pd.Series(dtype=float)), errors="coerce").fillna(0)
    rows_coverage = pd.to_numeric(smoke.get("rows_pass_coverage_ge5", pd.Series(dtype=float)), errors="coerce").fillna(0)
    common = pd.to_numeric(smoke.get("n_smoke_common_regions_with_reference", pd.Series(dtype=float)), errors="coerce").fillna(0)
    beta_valid = smoke.get("beta_range_valid", pd.Series(False, index=smoke.index)).fillna(False).astype(bool)
    completed = smoke.get("download_status", pd.Series("", index=smoke.index)).astype(str).isin(["completed", "skipped"])
    status = "passed" if bool((completed & rows_parseable.gt(0) & rows_coverage.gt(0) & beta_valid).all()) else "failed"
    state = {
        "timestamp": utc_now(),
        "status": status,
        "dataset": args.dataset,
        "n_files_requested": int(len(data)),
        "n_files_completed_or_cached": int(completed.sum()),
        "n_files_parseable": int(rows_parseable.gt(0).sum()),
        "n_files_beta_range_valid": int(beta_valid.sum()),
        "max_smoke_common_regions_with_reference": int(common.max()) if not common.empty else 0,
        "smoke_manifest": str(smoke_path.relative_to(ROOT)),
        "download_log": str(download_log.relative_to(ROOT)),
        "training_authorized": False,
        "bismark_authorized": False,
        "autoresearch_authorized": False,
        "headline_allowed": False,
        "next_action": "If this passes, implement a GSE83947 CX_report-to-region matrix adapter and run matrix gate before any Learn/Benchmark approval.",
    }
    write_json(state_path, state)
    print(json.dumps(state, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

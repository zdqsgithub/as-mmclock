#!/usr/bin/env python3
"""Build a pilot 5kb region matrix from GSE83947 CX_report smoke files.

This is a local processed-supplement matrix gate. It consumes only the three
downloaded CX_report files from the v14 adapter smoke, builds a small 5kb region
matrix, and checks overlap with the v8.2 reference matrix. It does not download
data, run Bismark, train models, or start autoresearch.
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")
V14_DIR = ROOT / "results" / "ralph_v14_public_data_rescue"
SMOKE_DIR = V14_DIR / "gse83947_processed_adapter_smoke"
OUT_DIR = V14_DIR / "gse83947_region_matrix_gate"
REFERENCE_MATRIX = ROOT / "results" / "multidataset_v8_2_prefilter_liftover" / "all_rrbs_region_matrix_5kb.parquet"
COMMON_REGION_GATE = 50_000
MIN_COVERAGE = 5
REGION_BIN_SIZE = 5_000


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str), encoding="utf-8")


def normalize_chrom(value: str) -> str | None:
    value = str(value).strip()
    if value.startswith("chr"):
        clean = value.removeprefix("chr")
    else:
        clean = value
    if clean in {"X", "Y", "M", "MT"}:
        return None
    if clean.isdigit():
        return f"chr{clean}"
    return None


def region_id(chrom: str, pos: int, bin_size: int = REGION_BIN_SIZE) -> str:
    start = (int(pos) // bin_size) * bin_size
    return f"{chrom}:{start}-{start + bin_size - 1}"


def iter_lines(path: Path) -> Iterable[str]:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", errors="replace") as handle:
        for line in handle:
            yield line.rstrip("\n")


def parse_cx_report_to_regions(path: Path, sample_id: str, min_coverage: int) -> tuple[pd.Series, dict]:
    """Parse Bismark CX_report: chrom, pos, strand, methylated, unmethylated, context, trinuc."""
    beta_sum: dict[str, float] = defaultdict(float)
    beta_n: dict[str, int] = defaultdict(int)
    coverage_values: list[float] = []
    stats = {
        "sample_id": sample_id,
        "local_path": str(path),
        "rows_total": 0,
        "rows_parseable_primary_autosomes": 0,
        "rows_pass_coverage_ge5": 0,
        "beta_min": None,
        "beta_max": None,
        "schema": "bismark_cx_report",
    }
    for line in iter_lines(path):
        if not line or line.startswith("#") or line.startswith("track"):
            continue
        stats["rows_total"] += 1
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        chrom = normalize_chrom(parts[0])
        if chrom is None:
            continue
        try:
            pos = int(float(parts[1]))
            methylated = float(parts[3])
            unmethylated = float(parts[4])
        except ValueError:
            continue
        total = methylated + unmethylated
        if total <= 0:
            continue
        beta = methylated / total
        if not (0.0 <= beta <= 1.0):
            continue
        stats["rows_parseable_primary_autosomes"] += 1
        stats["beta_min"] = beta if stats["beta_min"] is None else min(float(stats["beta_min"]), beta)
        stats["beta_max"] = beta if stats["beta_max"] is None else max(float(stats["beta_max"]), beta)
        if total < min_coverage:
            continue
        rid = region_id(chrom, pos)
        beta_sum[rid] += beta
        beta_n[rid] += 1
        coverage_values.append(total)
        stats["rows_pass_coverage_ge5"] += 1
    values = {rid: beta_sum[rid] / beta_n[rid] for rid in beta_sum}
    series = pd.Series(values, dtype=np.float32, name=sample_id)
    stats["n_regions"] = int(len(series))
    stats["median_coverage"] = None if not coverage_values else round(float(np.median(coverage_values)), 3)
    stats["beta_min"] = None if stats["beta_min"] is None else round(float(stats["beta_min"]), 6)
    stats["beta_max"] = None if stats["beta_max"] is None else round(float(stats["beta_max"]), 6)
    return series, stats


def load_reference_regions(path: Path) -> set[str]:
    if not path.exists():
        return set()
    ref = pd.read_parquet(path, columns=[])
    return set(map(str, ref.index))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--smoke-manifest", default=str(SMOKE_DIR / "gse83947_processed_adapter_smoke_manifest.csv"))
    parser.add_argument("--reference-matrix", default=str(REFERENCE_MATRIX))
    parser.add_argument("--out-dir", default=str(OUT_DIR))
    parser.add_argument("--min-coverage", type=int, default=MIN_COVERAGE)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    smoke = pd.read_csv(args.smoke_manifest)
    required = smoke[smoke.get("local_path", pd.Series("", index=smoke.index)).fillna("").astype(str).ne("")].copy()
    series_list: list[pd.Series] = []
    parse_stats: list[dict] = []
    for _, row in required.iterrows():
        sample_id = str(row["sample_accession"])
        path = Path(str(row["local_path"]))
        series, stats = parse_cx_report_to_regions(path, sample_id, min_coverage=args.min_coverage)
        stats.update(
            {
                "sample_age_weeks": row.get("sample_age_weeks", ""),
                "sample_tissue": row.get("sample_tissue", ""),
                "processed_supplement_name": row.get("processed_supplement_name", ""),
                "processed_supplement_bytes": row.get("processed_supplement_bytes", ""),
            }
        )
        series_list.append(series)
        parse_stats.append(stats)

    matrix = pd.concat(series_list, axis=1).sort_index() if series_list else pd.DataFrame()
    reference_regions = load_reference_regions(Path(args.reference_matrix))
    common_regions = sorted(set(map(str, matrix.index)).intersection(reference_regions))
    parse_stats_df = pd.DataFrame(parse_stats)
    parse_stats_path = out_dir / "GSE83947_sample_parse_stats.csv"
    region_matrix_path = out_dir / "GSE83947_region_matrix_5kb.parquet"
    region_stats_path = out_dir / "GSE83947_region_stats_5kb.csv"
    common_path = out_dir / "common_regions_with_v8_2.csv"
    manifest_path = out_dir / "GSE83947_matrix_gate_state.json"

    if not matrix.empty:
        matrix.to_parquet(region_matrix_path)
        present = matrix.notna().sum(axis=1)
        stats = pd.DataFrame(
            {
                "region_id": matrix.index,
                "n_samples_present": present.to_numpy(),
                "mean_beta": matrix.mean(axis=1, skipna=True).to_numpy(),
                "std_beta": matrix.std(axis=1, skipna=True).fillna(0.0).to_numpy(),
            }
        )
        chrom_start_end = stats["region_id"].str.extract(r"^(chr[^:]+):(\d+)-(\d+)$")
        stats["chrom"] = chrom_start_end[0]
        stats["start"] = pd.to_numeric(chrom_start_end[1], errors="coerce")
        stats["end"] = pd.to_numeric(chrom_start_end[2], errors="coerce")
        stats[["region_id", "chrom", "start", "end", "n_samples_present", "mean_beta", "std_beta"]].to_csv(region_stats_path, index=False)
    else:
        pd.DataFrame().to_parquet(region_matrix_path)
        pd.DataFrame(columns=["region_id", "chrom", "start", "end", "n_samples_present", "mean_beta", "std_beta"]).to_csv(region_stats_path, index=False)
    parse_stats_df.to_csv(parse_stats_path, index=False)
    pd.DataFrame({"region_id": common_regions}).to_csv(common_path, index=False)

    beta_min = float(matrix.min(skipna=True).min()) if not matrix.empty else math.nan
    beta_max = float(matrix.max(skipna=True).max()) if not matrix.empty else math.nan
    common_region_gate = len(common_regions) >= COMMON_REGION_GATE
    state = {
        "timestamp": utc_now(),
        "status": "completed",
        "dataset": "GSE83947",
        "schema": "bismark_cx_report",
        "min_coverage": args.min_coverage,
        "region_bin_size": REGION_BIN_SIZE,
        "n_samples": int(matrix.shape[1]) if not matrix.empty else 0,
        "n_regions": int(matrix.shape[0]) if not matrix.empty else 0,
        "reference_regions": int(len(reference_regions)),
        "common_regions_with_v8_2": int(len(common_regions)),
        "common_region_gate_50000": bool(common_region_gate),
        "beta_min": None if math.isnan(beta_min) else beta_min,
        "beta_max": None if math.isnan(beta_max) else beta_max,
        "beta_range_valid": bool(not matrix.empty and 0.0 <= beta_min <= beta_max <= 1.0),
        "region_matrix_path": str(region_matrix_path.relative_to(ROOT)),
        "region_stats_path": str(region_stats_path.relative_to(ROOT)),
        "sample_parse_stats_path": str(parse_stats_path.relative_to(ROOT)),
        "common_regions_path": str(common_path.relative_to(ROOT)),
        "training_authorized": False,
        "bismark_authorized": False,
        "autoresearch_authorized": False,
        "headline_allowed": False,
        "decision": "matrix_gate_failed_auxiliary_only" if not common_region_gate else "matrix_gate_passed_pending_metadata_and_training_approval",
        "next_action": "Do not train unless common-region gate and metadata/headline gates are approved explicitly.",
    }
    write_json(manifest_path, state)
    print(json.dumps(state, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

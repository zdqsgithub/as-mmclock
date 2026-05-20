#!/usr/bin/env python3
"""
Build a region-averaged methylation matrix from a CpG beta matrix.

The default is a Simpson/Meer-style 5 kb non-overlapping window summary for the
current Phase 0 GSE120137 matrix. Rows are genomic regions, columns are samples,
and values are mean beta across CpGs in the region, ignoring missing CpGs.
"""
from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")
PHASE0 = ROOT / "results" / "phase0"
EXCLUDE_CHROMS = {"chrX", "chrY", "chrM", "X", "Y", "MT", "M"}


def parse_cpg_id(cpg_id: object, bin_size: int) -> tuple[str, int, int, str] | None:
    text = str(cpg_id)
    if "_" not in text:
        return None
    chrom, pos_text = text.rsplit("_", 1)
    if chrom in EXCLUDE_CHROMS:
        return None
    try:
        pos = int(pos_text)
    except ValueError:
        return None
    start = (pos // bin_size) * bin_size
    end = start + bin_size - 1
    region_id = f"{chrom}:{start}-{end}"
    return chrom, start, end, region_id


def chrom_sort_value(chrom: str) -> tuple[int, str]:
    clean = chrom.removeprefix("chr")
    if clean.isdigit():
        return int(clean), ""
    return 10_000, clean


def build_region_matrix(
    input_matrix: Path,
    output_matrix: Path,
    stats_path: Path,
    bin_size: int,
    min_cpgs_per_region: int,
    min_sample_presence: float,
) -> dict:
    started = time.time()
    beta = pd.read_parquet(input_matrix)
    beta = beta.astype(np.float32)

    parsed_rows = []
    keep_positions = []
    for row_number, cpg_id in enumerate(beta.index):
        parsed = parse_cpg_id(cpg_id, bin_size)
        if parsed is None:
            continue
        chrom, start, end, region_id = parsed
        keep_positions.append(row_number)
        parsed_rows.append(
            {
                "cpg_id": str(cpg_id),
                "chrom": chrom,
                "start": start,
                "end": end,
                "region_id": region_id,
            }
        )

    if not parsed_rows:
        raise RuntimeError(f"No parseable autosomal CpG IDs found in {input_matrix}")

    cpg_meta = pd.DataFrame(parsed_rows)
    beta = beta.iloc[keep_positions]
    region_labels = pd.Series(cpg_meta["region_id"].values, index=beta.index)

    region_matrix = beta.groupby(region_labels, sort=False).mean().astype(np.float32)
    region_counts = cpg_meta.groupby("region_id", sort=False).size().rename("n_cpgs")
    region_coords = (
        cpg_meta[["region_id", "chrom", "start", "end"]]
        .drop_duplicates("region_id")
        .set_index("region_id")
    )
    stats = region_coords.join(region_counts)
    stats["n_samples_present"] = region_matrix.notna().sum(axis=1).astype(int)
    stats["mean_beta"] = region_matrix.mean(axis=1, skipna=True).astype(float)
    stats["std_beta"] = region_matrix.std(axis=1, skipna=True).astype(float)

    min_samples = int(math.ceil(region_matrix.shape[1] * min_sample_presence))
    keep = (stats["n_cpgs"] >= min_cpgs_per_region) & (stats["n_samples_present"] >= min_samples)
    stats = stats.loc[keep].copy()
    stats["_chrom_sort"] = [chrom_sort_value(chrom) for chrom in stats["chrom"]]
    stats = stats.sort_values(["_chrom_sort", "start", "end"]).drop(columns=["_chrom_sort"])
    region_matrix = region_matrix.loc[stats.index]

    output_matrix.parent.mkdir(parents=True, exist_ok=True)
    stats_path.parent.mkdir(parents=True, exist_ok=True)
    region_matrix.to_parquet(output_matrix, compression="zstd")
    stats.reset_index().to_csv(stats_path, index=False)

    manifest = {
        "source_matrix": str(input_matrix),
        "output_matrix": str(output_matrix),
        "stats_path": str(stats_path),
        "bin_size": int(bin_size),
        "min_cpgs_per_region": int(min_cpgs_per_region),
        "min_sample_presence": float(min_sample_presence),
        "min_samples_present": int(min_samples),
        "n_cpg_input": int(len(beta)),
        "n_samples": int(region_matrix.shape[1]),
        "n_regions": int(region_matrix.shape[0]),
        "exec_time_sec": round(time.time() - started, 1),
    }
    manifest_path = output_matrix.with_suffix(".manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input_matrix",
        default=str(PHASE0 / "beta_matrix_thompson.parquet"),
        help="CpG beta matrix parquet, rows=CpGs and columns=samples.",
    )
    parser.add_argument("--bin_size", type=int, default=5000)
    parser.add_argument("--min_cpgs_per_region", type=int, default=3)
    parser.add_argument("--min_sample_presence", type=float, default=0.8)
    parser.add_argument("--output_matrix", default=None)
    parser.add_argument("--stats_out", default=None)
    args = parser.parse_args()

    input_matrix = Path(args.input_matrix)
    suffix = f"{args.bin_size // 1000}kb" if args.bin_size % 1000 == 0 else str(args.bin_size)
    output_matrix = Path(args.output_matrix) if args.output_matrix else PHASE0 / f"region_matrix_{suffix}.parquet"
    stats_path = Path(args.stats_out) if args.stats_out else PHASE0 / f"region_stats_{suffix}.csv"

    manifest = build_region_matrix(
        input_matrix=input_matrix,
        output_matrix=output_matrix,
        stats_path=stats_path,
        bin_size=args.bin_size,
        min_cpgs_per_region=args.min_cpgs_per_region,
        min_sample_presence=args.min_sample_presence,
    )
    print("[Region] Complete")
    for key, value in manifest.items():
        print(f"  {key}: {value}")


if __name__ == "__main__":
    main()

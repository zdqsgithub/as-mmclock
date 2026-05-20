#!/usr/bin/env python3
"""Build coverage-weighted 5kb RRBS region matrices from processed supplements.

This v7.6 ETL bypasses the site-level beta matrix and aggregates methylated and
total coverage counts directly into 5kb windows per sample. It is intended for
datasets where processed supplements contain usable coverage/count columns.
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
import re
import tarfile
import time
from collections import defaultdict
from pathlib import Path
from typing import BinaryIO, Iterable

import numpy as np
import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")
OUT_DIR = ROOT / "results" / "multidataset_v7_6_weighted"

NC_TO_CHROM = {
    "NC_000067.6": "chr1",
    "NC_000068.7": "chr2",
    "NC_000069.6": "chr3",
    "NC_000070.6": "chr4",
    "NC_000071.6": "chr5",
    "NC_000072.6": "chr6",
    "NC_000073.6": "chr7",
    "NC_000074.6": "chr8",
    "NC_000075.6": "chr9",
    "NC_000076.6": "chr10",
    "NC_000077.6": "chr11",
    "NC_000078.6": "chr12",
    "NC_000079.6": "chr13",
    "NC_000080.6": "chr14",
    "NC_000081.6": "chr15",
    "NC_000082.6": "chr16",
    "NC_000083.6": "chr17",
    "NC_000084.6": "chr18",
    "NC_000085.6": "chr19",
    "NC_000086.7": "chrX",
    "NC_000087.7": "chrY",
    "NC_005089.1": "chrM",
}
EXCLUDE_CHROMS = {"chrX", "chrY", "chrM", "X", "Y", "M", "MT"}
COORD_RE = re.compile(r"\|ref\|([^|]+)\|:(\d+)")
GSM_RE = re.compile(r"(GSM\d+)")
GSE_RE = re.compile(r"(GSE\d+)")


def chrom_sort_value(chrom: str) -> tuple[int, str]:
    clean = chrom.removeprefix("chr")
    if clean.isdigit():
        return int(clean), ""
    return 10_000, clean


def sample_id_from_name(name: str) -> str:
    match = GSM_RE.search(Path(name).name)
    if not match:
        raise ValueError(f"Cannot infer GSM sample id from {name}")
    return match.group(1)


def infer_dataset(input_path: Path, explicit: str | None) -> str:
    if explicit:
        return explicit
    match = GSE_RE.search(str(input_path))
    if match:
        return match.group(1)
    raise ValueError(f"Cannot infer dataset from {input_path}")


def infer_schema(dataset: str, explicit: str) -> str:
    if explicit != "auto":
        return explicit
    if dataset == "GSE80672":
        return "gse80672_overlap_percentage_coverage"
    if dataset in {"GSE93957", "GSE121141"}:
        return "bismark_cov_per_sample_tar"
    raise ValueError(f"No automatic coverage-weighted schema for {dataset}")


def normalize_plain_chrom(chrom: str) -> str | None:
    chrom = str(chrom).strip()
    if chrom in EXCLUDE_CHROMS:
        return None
    if chrom.startswith("chr"):
        return chrom
    if chrom.isdigit():
        return f"chr{chrom}"
    return None


def normalize_refseq_coord(value: str) -> tuple[str, int] | None:
    match = COORD_RE.search(value)
    if not match:
        return None
    accession, pos_text = match.groups()
    chrom = NC_TO_CHROM.get(accession)
    if chrom is None or chrom in EXCLUDE_CHROMS:
        return None
    return chrom, int(pos_text)


def region_id_for(chrom: str, pos: int, bin_size: int) -> tuple[str, int, int]:
    start = (pos // bin_size) * bin_size
    end = start + bin_size - 1
    return f"{chrom}:{start}-{end}", start, end


def iter_tar_members(input_path: Path, sample_limit: int | None = None) -> Iterable[tuple[str, BinaryIO]]:
    with tarfile.open(input_path, "r") as tar:
        members = [
            member
            for member in tar.getmembers()
            if member.isfile() and member.name.endswith((".txt.gz", ".cov.gz", ".bismark.cov.gz"))
        ]
        for idx, member in enumerate(sorted(members, key=lambda item: item.name)):
            if sample_limit is not None and idx >= sample_limit:
                break
            extracted = tar.extractfile(member)
            if extracted is not None:
                yield member.name, extracted


def parse_bismark_sample(
    handle: BinaryIO,
    *,
    sample_id: str,
    min_coverage: int,
    bin_size: int,
    max_rows: int | None,
) -> tuple[pd.Series, pd.Series, pd.Series, dict]:
    meth_sum: dict[str, float] = defaultdict(float)
    total_sum: dict[str, float] = defaultdict(float)
    cpg_count: dict[str, int] = defaultdict(int)
    coord_map: dict[str, tuple[str, int, int]] = {}
    stats = {
        "sample_id": sample_id,
        "rows_total": 0,
        "rows_parseable_primary_autosomes": 0,
        "rows_pass_coverage": 0,
        "total_coverage_pass": 0.0,
        "schema": "bismark_cov_per_sample_tar",
    }
    with gzip.open(handle, "rt") as text:
        for row_idx, line in enumerate(text, start=1):
            if max_rows is not None and row_idx > max_rows:
                break
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            stats["rows_total"] += 1
            chrom = normalize_plain_chrom(parts[0])
            if chrom is None:
                continue
            try:
                pos = int(float(parts[1]))
                methylated = float(parts[4])
                unmethylated = float(parts[5])
            except ValueError:
                continue
            stats["rows_parseable_primary_autosomes"] += 1
            total = methylated + unmethylated
            if total < min_coverage or total <= 0:
                continue
            region_id, start, end = region_id_for(chrom, pos, bin_size)
            coord_map[region_id] = (chrom, start, end)
            meth_sum[region_id] += methylated
            total_sum[region_id] += total
            cpg_count[region_id] += 1
            stats["rows_pass_coverage"] += 1
            stats["total_coverage_pass"] += total
    beta = {region: meth_sum[region] / total_sum[region] for region in meth_sum if total_sum[region] > 0}
    return (
        pd.Series(beta, dtype=np.float32, name=sample_id),
        pd.Series(total_sum, dtype=np.float32, name=sample_id),
        pd.Series(cpg_count, dtype=np.float32, name=sample_id),
        stats | {"n_regions_pass": len(beta), "coord_map": coord_map},
    )


def parse_overlap_sample(
    handle: BinaryIO,
    *,
    sample_id: str,
    min_coverage: int,
    bin_size: int,
    max_rows: int | None,
) -> tuple[pd.Series, pd.Series, pd.Series, dict]:
    meth_sum: dict[str, float] = defaultdict(float)
    total_sum: dict[str, float] = defaultdict(float)
    cpg_count: dict[str, int] = defaultdict(int)
    coord_map: dict[str, tuple[str, int, int]] = {}
    stats = {
        "sample_id": sample_id,
        "rows_total": 0,
        "rows_parseable_primary_autosomes": 0,
        "rows_pass_coverage": 0,
        "total_coverage_pass": 0.0,
        "schema": "gse80672_overlap_percentage_coverage",
    }
    with gzip.open(handle, "rt") as text:
        header = text.readline().rstrip("\n").split("\t")
        stats["header"] = header
        if len(header) < 3 or header[0] != "index" or "Percentage" not in header[1] or "Coverage" not in header[2]:
            raise ValueError(f"Unsupported overlap schema for {sample_id}: {header}")
        for row_idx, line in enumerate(text, start=1):
            if max_rows is not None and row_idx > max_rows:
                break
            stats["rows_total"] += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            coord = normalize_refseq_coord(parts[0])
            if coord is None:
                continue
            try:
                pct = float(parts[1])
                total = float(parts[2])
            except ValueError:
                continue
            stats["rows_parseable_primary_autosomes"] += 1
            if total < min_coverage or total <= 0:
                continue
            chrom, pos = coord
            region_id, start, end = region_id_for(chrom, pos, bin_size)
            coord_map[region_id] = (chrom, start, end)
            meth_sum[region_id] += (pct / 100.0) * total
            total_sum[region_id] += total
            cpg_count[region_id] += 1
            stats["rows_pass_coverage"] += 1
            stats["total_coverage_pass"] += total
    beta = {region: meth_sum[region] / total_sum[region] for region in meth_sum if total_sum[region] > 0}
    return (
        pd.Series(beta, dtype=np.float32, name=sample_id),
        pd.Series(total_sum, dtype=np.float32, name=sample_id),
        pd.Series(cpg_count, dtype=np.float32, name=sample_id),
        stats | {"n_regions_pass": len(beta), "coord_map": coord_map},
    )


def build_stats(
    beta_matrix: pd.DataFrame,
    coverage_matrix: pd.DataFrame,
    cpg_count_matrix: pd.DataFrame,
    coord_map: dict[str, tuple[str, int, int]],
    min_sample_presence: float,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    min_samples = int(math.ceil(beta_matrix.shape[1] * min_sample_presence))
    n_present = beta_matrix.notna().sum(axis=1).astype(int)
    keep = n_present >= min_samples
    beta_matrix = beta_matrix.loc[keep].copy()
    coverage_matrix = coverage_matrix.reindex(beta_matrix.index)
    cpg_count_matrix = cpg_count_matrix.reindex(beta_matrix.index)
    stats = pd.DataFrame(index=beta_matrix.index)
    stats["chrom"] = [coord_map[region][0] for region in beta_matrix.index]
    stats["start"] = [coord_map[region][1] for region in beta_matrix.index]
    stats["end"] = [coord_map[region][2] for region in beta_matrix.index]
    stats["n_samples_present"] = beta_matrix.notna().sum(axis=1).astype(int)
    stats["mean_beta"] = beta_matrix.mean(axis=1, skipna=True).astype(float)
    stats["std_beta"] = beta_matrix.std(axis=1, skipna=True).astype(float)
    stats["mean_total_coverage"] = coverage_matrix.mean(axis=1, skipna=True).astype(float)
    stats["median_total_coverage"] = coverage_matrix.median(axis=1, skipna=True).astype(float)
    stats["mean_cpgs_per_sample"] = cpg_count_matrix.mean(axis=1, skipna=True).astype(float)
    stats["_chrom_sort"] = [chrom_sort_value(chrom) for chrom in stats["chrom"]]
    stats = stats.sort_values(["_chrom_sort", "start", "end"]).drop(columns=["_chrom_sort"])
    beta_matrix = beta_matrix.loc[stats.index]
    coverage_matrix = coverage_matrix.loc[stats.index]
    cpg_count_matrix = cpg_count_matrix.loc[stats.index]
    return beta_matrix, coverage_matrix, cpg_count_matrix, stats.reset_index(names="region_id")


def convert(input_path: Path, args: argparse.Namespace) -> dict:
    started = time.time()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    dataset = infer_dataset(input_path, args.dataset)
    schema = infer_schema(dataset, args.schema)
    parser = parse_overlap_sample if schema == "gse80672_overlap_percentage_coverage" else parse_bismark_sample

    beta_series: dict[str, pd.Series] = {}
    coverage_series: dict[str, pd.Series] = {}
    cpg_count_series: dict[str, pd.Series] = {}
    sample_stats = []
    coord_map: dict[str, tuple[str, int, int]] = {}
    errors = []
    for name, handle in iter_tar_members(input_path, args.sample_limit):
        try:
            sample_id = sample_id_from_name(name)
            beta, coverage, cpg_count, stats = parser(
                handle,
                sample_id=sample_id,
                min_coverage=args.min_coverage,
                bin_size=args.bin_size,
                max_rows=args.max_rows,
            )
            if not beta.empty:
                beta_series[sample_id] = beta
                coverage_series[sample_id] = coverage
                cpg_count_series[sample_id] = cpg_count
                coord_map.update(stats.pop("coord_map"))
            sample_stats.append(stats)
            if len(sample_stats) % 25 == 0:
                print(f"  Parsed {len(sample_stats)} samples for {dataset}")
        except Exception as exc:
            errors.append({"member": name, "error": str(exc)[:500]})

    if not beta_series:
        raise SystemExit(f"No parseable coverage-weighted samples for {dataset}; first errors={errors[:3]}")
    beta_matrix = pd.DataFrame(beta_series).astype(np.float32)
    coverage_matrix = pd.DataFrame(coverage_series).astype(np.float32)
    cpg_count_matrix = pd.DataFrame(cpg_count_series).astype(np.float32)
    n_regions_before_presence = int(beta_matrix.shape[0])
    beta_matrix, coverage_matrix, cpg_count_matrix, region_stats = build_stats(
        beta_matrix,
        coverage_matrix,
        cpg_count_matrix,
        coord_map,
        args.region_min_sample_presence,
    )

    matrix_path = out_dir / f"{dataset}_coverage_weighted_region_matrix_5kb.parquet"
    stats_path = out_dir / f"{dataset}_coverage_weighted_region_stats_5kb.csv"
    coverage_path = out_dir / f"{dataset}_coverage_weighted_region_total_coverage_5kb.parquet"
    cpg_count_path = out_dir / f"{dataset}_coverage_weighted_region_cpg_counts_5kb.parquet"
    beta_matrix.to_parquet(matrix_path, compression="zstd")
    region_stats.to_csv(stats_path, index=False)
    if args.write_coverage_matrix:
        coverage_matrix.to_parquet(coverage_path, compression="zstd")
        cpg_count_matrix.to_parquet(cpg_count_path, compression="zstd")

    manifest = {
        "status": "completed",
        "dataset": dataset,
        "schema": schema,
        "input_path": str(input_path),
        "aggregation": "coverage_weighted_5kb_region",
        "min_coverage": int(args.min_coverage),
        "region_bin_size": int(args.bin_size),
        "region_min_sample_presence": float(args.region_min_sample_presence),
        "n_samples_parsed": int(len(beta_series)),
        "n_regions_before_presence_filter": n_regions_before_presence,
        "n_regions": int(beta_matrix.shape[0]),
        "matrix_path": str(matrix_path),
        "stats_path": str(stats_path),
        "coverage_matrix_path": str(coverage_path) if args.write_coverage_matrix else None,
        "cpg_count_matrix_path": str(cpg_count_path) if args.write_coverage_matrix else None,
        "n_errors": int(len(errors)),
        "errors": errors[:20],
        "exec_time_sec": round(time.time() - started, 1),
    }
    pd.DataFrame(sample_stats).to_csv(out_dir / f"{dataset}_coverage_weighted_sample_parse_stats.csv", index=False)
    (out_dir / f"{dataset}_coverage_weighted_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(f"[CoverageWeighted] {dataset} -> {matrix_path} shape={beta_matrix.shape}")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--dataset", default=None)
    parser.add_argument(
        "--schema",
        default="auto",
        choices=["auto", "gse80672_overlap_percentage_coverage", "bismark_cov_per_sample_tar"],
    )
    parser.add_argument("--min_coverage", type=int, default=5)
    parser.add_argument("--region_min_sample_presence", type=float, default=0.8)
    parser.add_argument("--bin_size", type=int, default=5000)
    parser.add_argument("--out_dir", default=str(OUT_DIR))
    parser.add_argument("--sample_limit", type=int, default=None)
    parser.add_argument("--max_rows", type=int, default=None)
    parser.add_argument("--write_coverage_matrix", action="store_true")
    args = parser.parse_args()

    manifest = convert(Path(args.input), args)
    print(json.dumps(manifest, indent=2)[:4000])


if __name__ == "__main__":
    main()


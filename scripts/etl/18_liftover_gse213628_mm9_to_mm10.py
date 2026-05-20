#!/usr/bin/env python3
"""Lift mm9-style mouse methylation matrices to mm10 and rebuild 5kb regions.

The script uses a UCSC chain file and performs point liftover on row positions.
For CpG matrices this maps CpG coordinates. For tile matrices this maps the
tile start coordinate, which is acceptable for overlap diagnostics but should
remain lower priority than base-resolution CpG inputs. UCSC ``mm9ToMm10``
stores mm9 in the chain target fields and mm10 in the query fields, so this
script maps source target fields to destination query fields.
"""
from __future__ import annotations

import argparse
import bisect
import gzip
import json
import math
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")
DEFAULT_BETA = ROOT / "results" / "multidataset_v8_gse213628" / "GSE213628_beta_matrix.parquet"
DEFAULT_OUT = ROOT / "results" / "multidataset_v8_1_gse213628_lifted"
DEFAULT_CHAIN = ROOT / "resources" / "liftover" / "mm9ToMm10.over.chain.gz"
EXCLUDE_CHROMS = {"chrX", "chrY", "chrM", "X", "Y", "M", "MT"}


@dataclass(frozen=True)
class ChainInterval:
    source_start: int
    source_end: int
    target_chrom: str
    target_block_start: int
    target_size: int
    target_strand: str
    score: int

    def map_pos0(self, source_pos0: int) -> int:
        offset = source_pos0 - self.source_start
        target_strand_coord = self.target_block_start + offset
        if self.target_strand == "+":
            return target_strand_coord
        return self.target_size - target_strand_coord - 1


def chrom_sort_value(chrom: str) -> tuple[int, str]:
    clean = chrom.removeprefix("chr")
    if clean.isdigit():
        return int(clean), ""
    return 10_000, clean


def parse_chain(chain_path: Path) -> dict[str, tuple[list[int], list[ChainInterval]]]:
    chains: dict[str, list[ChainInterval]] = {}
    opener = gzip.open if str(chain_path).endswith(".gz") else open
    with opener(chain_path, "rt") as handle:
        header = None
        source_pos = target_pos = 0
        for line in handle:
            line = line.strip()
            if not line:
                header = None
                continue
            if line.startswith("chain "):
                parts = line.split()
                header = {
                    "score": int(parts[1]),
                    "source_chrom": parts[2],
                    "source_strand": parts[4],
                    "source_start": int(parts[5]),
                    "target_chrom": parts[7],
                    "target_size": int(parts[8]),
                    "target_strand": parts[9],
                    "target_start": int(parts[10]),
                }
                if header["source_strand"] != "+":
                    header = None
                    continue
                source_pos = int(header["source_start"])
                target_pos = int(header["target_start"])
                continue
            if header is None:
                continue
            parts = [int(item) for item in line.split()]
            size = parts[0]
            interval = ChainInterval(
                source_start=source_pos,
                source_end=source_pos + size,
                target_chrom=header["target_chrom"],
                target_block_start=target_pos,
                target_size=header["target_size"],
                target_strand=header["target_strand"],
                score=header["score"],
            )
            chains.setdefault(header["source_chrom"], []).append(interval)
            source_pos += size
            target_pos += size
            if len(parts) == 3:
                source_pos += parts[1]
                target_pos += parts[2]
    indexed = {}
    for chrom, intervals in chains.items():
        intervals = sorted(intervals, key=lambda item: (item.source_start, -item.score))
        indexed[chrom] = ([item.source_start for item in intervals], intervals)
    return indexed


def parse_cpg_id(cpg_id: object) -> tuple[str, int] | None:
    try:
        chrom, pos_text = str(cpg_id).rsplit("_", 1)
        pos = int(pos_text)
    except Exception:
        return None
    if chrom in EXCLUDE_CHROMS or pos <= 0:
        return None
    return chrom, pos


def liftover_point(
    chain_index: dict[str, tuple[list[int], list[ChainInterval]]],
    chrom: str,
    pos_1based: int,
) -> tuple[str | None, int | None, str]:
    if chrom not in chain_index:
        return None, None, "missing_chrom_chain"
    source_pos0 = pos_1based - 1
    starts, intervals = chain_index[chrom]
    idx = bisect.bisect_right(starts, source_pos0)
    best = None
    scan = idx - 1
    while scan >= 0 and intervals[scan].source_start <= source_pos0:
        interval = intervals[scan]
        if interval.source_start <= source_pos0 < interval.source_end:
            if best is None or interval.score > best.score:
                best = interval
        if scan > 0 and intervals[scan - 1].source_end <= source_pos0:
            break
        scan -= 1
    if best is None:
        return None, None, "unmapped"
    target_pos0 = best.map_pos0(source_pos0)
    return best.target_chrom, target_pos0 + 1, best.query_strand


def build_region_matrix(beta_df: pd.DataFrame, bin_size: int, min_sample_presence: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    parsed_rows = []
    keep_positions = []
    for row_number, cpg_id in enumerate(beta_df.index):
        parsed = parse_cpg_id(cpg_id)
        if parsed is None:
            continue
        chrom, pos = parsed
        start = (pos // bin_size) * bin_size
        end = start + bin_size - 1
        parsed_rows.append(
            {
                "cpg_id": str(cpg_id),
                "chrom": chrom,
                "start": start,
                "end": end,
                "region_id": f"{chrom}:{start}-{end}",
            }
        )
        keep_positions.append(row_number)
    cpg_meta = pd.DataFrame(parsed_rows)
    region_labels = pd.Series(cpg_meta["region_id"].values, index=beta_df.index[keep_positions])
    region_matrix = beta_df.iloc[keep_positions].groupby(region_labels, sort=False).mean().astype(np.float32)
    region_counts = cpg_meta.groupby("region_id", sort=False).size().rename("n_cpgs")
    region_coords = cpg_meta[["region_id", "chrom", "start", "end"]].drop_duplicates("region_id").set_index("region_id")
    stats = region_coords.join(region_counts)
    stats["n_samples_present"] = region_matrix.notna().sum(axis=1).astype(int)
    stats["mean_beta"] = region_matrix.mean(axis=1, skipna=True).astype(float)
    stats["std_beta"] = region_matrix.std(axis=1, skipna=True).astype(float)
    min_samples = int(math.ceil(region_matrix.shape[1] * min_sample_presence))
    keep = stats["n_samples_present"] >= min_samples
    stats = stats.loc[keep].copy()
    stats["_chrom_sort"] = [chrom_sort_value(chrom) for chrom in stats["chrom"]]
    stats = stats.sort_values(["_chrom_sort", "start", "end"]).drop(columns=["_chrom_sort"])
    return region_matrix.loc[stats.index], stats.reset_index()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--beta_matrix", default=str(DEFAULT_BETA))
    parser.add_argument("--chain", default=str(DEFAULT_CHAIN))
    parser.add_argument("--out_dir", default=str(DEFAULT_OUT))
    parser.add_argument("--dataset", default="GSE213628")
    parser.add_argument("--source_schema", default="bismark_cov_per_sample_tar")
    parser.add_argument("--source_assembly", default="mm9")
    parser.add_argument("--target_assembly", default="mm10")
    parser.add_argument("--bin_size", type=int, default=5000)
    parser.add_argument("--region_min_sample_presence", type=float, default=0.8)
    parser.add_argument("--max_cpgs", type=int, default=None, help="Debug limit; omit for full matrix.")
    args = parser.parse_args()

    started = time.time()
    beta_path = Path(args.beta_matrix)
    chain_path = Path(args.chain)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if not beta_path.exists():
        raise SystemExit(f"Missing beta matrix: {beta_path}")
    if not chain_path.exists():
        raise SystemExit(f"Missing chain file: {chain_path}")

    print(f"[liftover] loading chain {chain_path}", flush=True)
    chain_index = parse_chain(chain_path)
    print(f"[liftover] chain chromosomes={len(chain_index)}", flush=True)
    beta = pd.read_parquet(beta_path)
    if args.max_cpgs is not None:
        beta = beta.iloc[: args.max_cpgs].copy()
    print(f"[liftover] beta shape={beta.shape}", flush=True)

    mapping_rows = []
    target_ids = []
    keep_positions = []
    failure_counts: dict[str, int] = {}
    for idx, cpg_id in enumerate(beta.index):
        parsed = parse_cpg_id(cpg_id)
        if parsed is None:
            reason = "bad_cpg_id_or_excluded_chrom"
            failure_counts[reason] = failure_counts.get(reason, 0) + 1
            continue
        chrom, pos = parsed
        target_chrom, target_pos, chain_strand = liftover_point(chain_index, chrom, pos)
        if target_chrom is None or target_pos is None or target_chrom in EXCLUDE_CHROMS:
            reason = chain_strand if target_chrom is None else "target_excluded_chrom"
            failure_counts[reason] = failure_counts.get(reason, 0) + 1
            continue
        target_id = f"{target_chrom}_{target_pos}"
        target_ids.append(target_id)
        keep_positions.append(idx)
        mapping_rows.append(
            {
                "source_cpg_id": str(cpg_id),
                "source_chrom": chrom,
                "source_pos": pos,
                "target_cpg_id": target_id,
                "target_chrom": target_chrom,
                "target_pos": int(target_pos),
                "chain_target_strand": chain_strand,
            }
        )
        if len(mapping_rows) % 250_000 == 0:
            print(f"[liftover] mapped {len(mapping_rows)} CpGs", flush=True)

    if not keep_positions:
        raise SystemExit("No CpGs lifted over.")
    lifted = beta.iloc[keep_positions].copy()
    lifted.index = target_ids
    before_duplicate_cpgs = int(lifted.index.duplicated().sum())
    if before_duplicate_cpgs:
        lifted = lifted.groupby(lifted.index, sort=False).mean().astype(np.float32)
    else:
        lifted = lifted.astype(np.float32)

    lifted_beta_path = out_dir / f"{args.dataset}_beta_matrix.parquet"
    lifted.to_parquet(lifted_beta_path, compression="zstd")
    mapping_df = pd.DataFrame(mapping_rows)
    mapping_path = out_dir / f"{args.dataset}_liftover_cpg_map.parquet"
    mapping_df.to_parquet(mapping_path, compression="zstd")

    region_matrix, region_stats = build_region_matrix(
        lifted,
        bin_size=args.bin_size,
        min_sample_presence=args.region_min_sample_presence,
    )
    region_path = out_dir / f"{args.dataset}_region_matrix_5kb.parquet"
    stats_path = out_dir / f"{args.dataset}_region_stats_5kb.csv"
    region_matrix.to_parquet(region_path, compression="zstd")
    region_stats.to_csv(stats_path, index=False)
    manifest = {
        "status": "completed",
        "dataset": args.dataset,
        "source_assembly": args.source_assembly,
        "target_assembly": args.target_assembly,
        "chain_file": str(chain_path),
        "input_beta_matrix": str(beta_path),
        "n_cpg_input": int(beta.shape[0]),
        "n_cpg_lifted_before_duplicate_collapse": int(len(keep_positions)),
        "n_cpg_duplicate_targets": before_duplicate_cpgs,
        "n_cpg_lifted_after_duplicate_collapse": int(lifted.shape[0]),
        "liftover_failure_counts": failure_counts,
        "lifted_beta_matrix_path": str(lifted_beta_path),
        "liftover_map_path": str(mapping_path),
        "region_matrix_path": str(region_path),
        "region_stats_path": str(stats_path),
        "region_bin_size": int(args.bin_size),
        "region_min_sample_presence": float(args.region_min_sample_presence),
        "n_regions": int(region_matrix.shape[0]),
        "sex_mt_region_rows": int(region_stats["chrom"].isin(["chrX", "chrY", "chrM", "X", "Y", "M", "MT"]).sum()),
        "exec_time_sec": round(time.time() - started, 1),
    }
    (out_dir / f"{args.dataset}_liftover_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    (out_dir / f"{args.dataset}_matrix_manifest.json").write_text(
        json.dumps(
            {
                "status": "completed",
                "dataset": args.dataset,
                "schema": f"{args.source_schema}_lifted_{args.source_assembly}_to_{args.target_assembly}",
                "source_manifest": str(out_dir / f"{args.dataset}_liftover_manifest.json"),
                "n_samples_parsed": int(lifted.shape[1]),
                "n_samples_metadata_overlap": int(lifted.shape[1]),
                "metadata_age_known_overlap": int(lifted.shape[1]),
                "n_cpg_after_presence_filter": int(lifted.shape[0]),
                "n_regions": int(region_matrix.shape[0]),
                "sex_mt_region_rows": int(manifest["sex_mt_region_rows"]),
                "beta_matrix_path": str(lifted_beta_path),
                "region_matrix_path": str(region_path),
                "region_stats_path": str(stats_path),
                "assembly_harmonization": "mm9_to_mm10_ucsc_liftover",
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(json.dumps(manifest, indent=2)[:4000])


if __name__ == "__main__":
    main()

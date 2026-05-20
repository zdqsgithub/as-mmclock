#!/usr/bin/env python3
"""Convert processed GEO RRBS methylation supplements into model matrices.

Supported schemas:
- ``gse80672_overlap_percentage_coverage``: per-sample overlap tables with
  RefSeq coordinates, methylation percentage, and coverage.
- ``bismark_cov_per_sample_tar``: GEO tar archives of per-sample Bismark
  coverage files with ``chrom start end pct methylated unmethylated`` columns.
- ``gse60012_100bp_tile_matrix``: aggregated 100 bp tile methylation matrix.

Outputs are dataset-aware: ``{dataset}_beta_matrix.parquet``,
``{dataset}_region_matrix_5kb.parquet``, and ``{dataset}_matrix_manifest.json``.
"""
from __future__ import annotations

import argparse
import bisect
import gzip
import json
import math
import re
import tarfile
import tempfile
import time
from pathlib import Path
from dataclasses import dataclass
from typing import BinaryIO, Iterable

import numpy as np
import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")
MULTI_DIR = ROOT / "results" / "multidataset"
META_FILE = ROOT / "metadata" / "unified_sample_metadata.csv"

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


@dataclass(frozen=True)
class LiftOverInterval:
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


class LiftOverMapper:
    """Point liftover using UCSC chain files.

    UCSC ``mm9ToMm10`` chain files store the source assembly in the target
    fields and the destination assembly in the query fields. Query coordinates
    on ``-`` chains are reverse-strand coordinates, so they are converted back
    to forward coordinates during point mapping.
    """

    def __init__(self, chain_path: Path):
        self.chain_path = chain_path
        self.index = self._parse_chain(chain_path)
        self.cache: dict[tuple[str, int], tuple[str | None, int | None, str]] = {}

    @staticmethod
    def _parse_chain(chain_path: Path) -> dict[str, tuple[list[int], list[LiftOverInterval]]]:
        chains: dict[str, list[LiftOverInterval]] = {}
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
                    if len(parts) < 13:
                        header = None
                        continue
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
                interval = LiftOverInterval(
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

    def map_point(self, chrom: str, pos_1based: int) -> tuple[str | None, int | None, str]:
        key = (chrom, int(pos_1based))
        cached = self.cache.get(key)
        if cached is not None:
            return cached
        if chrom not in self.index:
            result = (None, None, "missing_chrom_chain")
            self.cache[key] = result
            return result
        source_pos0 = int(pos_1based) - 1
        starts, intervals = self.index[chrom]
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
            result = (None, None, "unmapped")
        else:
            result = (best.target_chrom, best.map_pos0(source_pos0) + 1, best.target_strand)
        self.cache[key] = result
        return result

    def map_interval(self, chrom: str, start_1based: int, end_1based: int) -> tuple[str | None, int | None, str]:
        points = [int(start_1based), int((int(start_1based) + int(end_1based)) // 2), int(end_1based)]
        mapped = [self.map_point(chrom, pos) for pos in points]
        ok = [(m_chrom, m_pos) for m_chrom, m_pos, _ in mapped if m_chrom is not None and m_pos is not None]
        if not ok:
            return None, None, "interval_unmapped"
        chrom_counts: dict[str, int] = {}
        for m_chrom, _ in ok:
            chrom_counts[m_chrom] = chrom_counts.get(m_chrom, 0) + 1
        best_chrom = max(chrom_counts, key=chrom_counts.get)
        positions = sorted(int(m_pos) for m_chrom, m_pos in ok if m_chrom == best_chrom)
        if len(positions) < 2 and len(ok) > 1:
            return None, None, "interval_split_mapping"
        return best_chrom, positions[len(positions) // 2], "interval_mapped"


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
    return "GSE80672"


def infer_schema(dataset: str, input_path: Path, explicit: str) -> str:
    if explicit != "auto":
        return explicit
    if dataset == "GSE80672":
        return "gse80672_overlap_percentage_coverage"
    if dataset in {"GSE93957", "GSE121141", "GSE213628"}:
        return "bismark_cov_per_sample_tar"
    if dataset == "GSE60012":
        return "gse60012_100bp_tile_matrix"
    name = input_path.name
    if name.endswith(".tar"):
        return "bismark_cov_per_sample_tar"
    return "unknown"


def normalize_coord(value: str) -> str | None:
    match = COORD_RE.search(value)
    if not match:
        return None
    accession, pos_text = match.groups()
    chrom = NC_TO_CHROM.get(accession)
    if chrom is None or chrom in EXCLUDE_CHROMS:
        return None
    return f"{chrom}_{int(pos_text)}"


def normalize_plain_chrom(chrom: str) -> str | None:
    chrom = str(chrom).strip()
    if chrom in EXCLUDE_CHROMS:
        return None
    if chrom.startswith("chr"):
        return chrom
    if chrom.isdigit():
        return f"chr{chrom}"
    return None


def open_text_member(path: Path | None = None, fileobj: BinaryIO | None = None):
    if path is not None:
        if str(path).endswith(".gz"):
            return gzip.open(path, "rt")
        return path.open("rt", encoding="utf-8")
    if fileobj is not None:
        return gzip.open(fileobj, "rt")
    raise ValueError("Either path or fileobj is required")


def parse_overlap_table(
    *,
    sample_id: str,
    path: Path | None = None,
    fileobj: BinaryIO | None = None,
    min_coverage: int,
    max_rows: int | None = None,
    coordinate_mapper: LiftOverMapper | None = None,
) -> tuple[pd.Series, dict]:
    records: dict[str, float] = {}
    stats = {
        "sample_id": sample_id,
        "rows_total": 0,
        "rows_parseable_primary_autosomes": 0,
        "rows_pass_coverage": 0,
        "schema": "gse80672_overlap_percentage_coverage",
    }
    with open_text_member(path=path, fileobj=fileobj) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        stats["header"] = header
        if len(header) < 3 or header[0] != "index" or "Percentage" not in header[1] or "Coverage" not in header[2]:
            raise ValueError(f"Unsupported processed methylation schema for {sample_id}: {header}")
        for row_idx, line in enumerate(handle, start=1):
            if max_rows is not None and row_idx > max_rows:
                break
            stats["rows_total"] += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            coord = normalize_coord(parts[0])
            if coord is None:
                continue
            stats["rows_parseable_primary_autosomes"] += 1
            try:
                pct = float(parts[1])
                coverage = float(parts[2])
            except ValueError:
                continue
            if coverage < min_coverage:
                continue
            stats["rows_pass_coverage"] += 1
            records[coord] = pct / 100.0
    return pd.Series(records, dtype=np.float32, name=sample_id), stats


def parse_bismark_cov_table(
    *,
    sample_id: str,
    path: Path | None = None,
    fileobj: BinaryIO | None = None,
    min_coverage: int,
    max_rows: int | None = None,
    coordinate_mapper: LiftOverMapper | None = None,
) -> tuple[pd.Series, dict]:
    records: dict[str, float] = {}
    stats = {
        "sample_id": sample_id,
        "rows_total": 0,
        "rows_parseable_primary_autosomes": 0,
        "rows_pass_coverage": 0,
        "schema": "bismark_cov_per_sample_tar",
    }
    if coordinate_mapper is not None:
        stats["rows_liftover_mapped"] = 0
        stats["liftover_failure_counts"] = {}
    with open_text_member(path=path, fileobj=fileobj) as handle:
        for row_idx, line in enumerate(handle, start=1):
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
            stats["rows_pass_coverage"] += 1
            if coordinate_mapper is not None:
                mapped_chrom, mapped_pos, reason = coordinate_mapper.map_point(chrom, pos)
                if mapped_chrom is None or mapped_pos is None or mapped_chrom in EXCLUDE_CHROMS:
                    failures = stats["liftover_failure_counts"]
                    failures[reason] = failures.get(reason, 0) + 1
                    continue
                chrom = mapped_chrom
                pos = mapped_pos
                stats["rows_liftover_mapped"] += 1
            records[f"{chrom}_{pos}"] = methylated / total
    return pd.Series(records, dtype=np.float32, name=sample_id), stats


def iter_input_members(input_path: Path, sample_limit: int | None = None) -> Iterable[tuple[str, Path | None, bytes | None]]:
    if input_path.suffix == ".tar":
        with tarfile.open(input_path, "r") as tar:
            members = [
                m
                for m in tar.getmembers()
                if m.isfile() and m.name.endswith((".txt.gz", ".cov.gz", ".bismark.cov.gz"))
            ]
            members = sorted(members, key=lambda m: m.name)
            for idx, member in enumerate(members):
                if sample_limit is not None and idx >= sample_limit:
                    break
                extracted = tar.extractfile(member)
                if extracted is None:
                    continue
                yield member.name, None, extracted.read()
    else:
        yield input_path.name, input_path, None


def sanitize_sample_label(value: str) -> str:
    value = re.sub(r"\s+", "_", str(value).strip())
    value = re.sub(r"[^A-Za-z0-9_.+-]+", "_", value)
    return value.strip("_") or "sample"


def parse_gse60012_tile_matrix(
    input_path: Path,
    max_rows: int | None = None,
    coordinate_mapper: LiftOverMapper | None = None,
) -> tuple[pd.DataFrame, list[dict], dict]:
    opener = gzip.open if str(input_path).endswith(".gz") else open
    with opener(input_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
    if len(header) < 4 or header[0] not in {"#Chr", "Chr", "chr"}:
        raise ValueError(f"Unsupported GSE60012 tile matrix header: {header[:8]}")

    raw_sample_labels = header[3:]
    sample_ids = [
        f"GSE60012_tile_{idx:03d}_{sanitize_sample_label(label)}"
        for idx, label in enumerate(raw_sample_labels, start=1)
    ]
    df = pd.read_csv(
        input_path,
        sep="\t",
        na_values=["NA", "NaN", ""],
        nrows=max_rows,
        compression="infer",
        low_memory=False,
    )
    chrom_col = df.columns[0]
    start_col = df.columns[1]
    end_col = df.columns[2]
    chroms = df[chrom_col].map(normalize_plain_chrom)
    keep = chroms.notna()
    df = df.loc[keep].copy()
    chroms = chroms.loc[keep]
    liftover_failure_counts: dict[str, int] = {}
    if coordinate_mapper is None:
        cpg_ids = chroms.astype(str) + "_" + df[start_col].astype(int).astype(str)
        values = df.iloc[:, 3:].apply(pd.to_numeric, errors="coerce").astype(np.float32)
    else:
        start_values = pd.to_numeric(df[start_col], errors="coerce")
        end_values = pd.to_numeric(df[end_col], errors="coerce")
        target_ids = []
        keep_positions = []
        for row_number, (chrom, start_value, end_value) in enumerate(
            zip(chroms.to_numpy(), start_values.to_numpy(), end_values.to_numpy(), strict=True)
        ):
            if pd.isna(start_value) or pd.isna(end_value):
                liftover_failure_counts["bad_interval"] = liftover_failure_counts.get("bad_interval", 0) + 1
                continue
            mapped_chrom, mapped_pos, reason = coordinate_mapper.map_interval(
                str(chrom),
                int(start_value),
                int(end_value),
            )
            if mapped_chrom is None or mapped_pos is None or mapped_chrom in EXCLUDE_CHROMS:
                liftover_failure_counts[reason] = liftover_failure_counts.get(reason, 0) + 1
                continue
            target_ids.append(f"{mapped_chrom}_{int(mapped_pos)}")
            keep_positions.append(row_number)
        values = df.iloc[keep_positions, 3:].apply(pd.to_numeric, errors="coerce").astype(np.float32)
        cpg_ids = pd.Series(target_ids, index=values.index)
    values.columns = sample_ids
    values.index = cpg_ids.values
    stats = [
        {
            "sample_id": sample_id,
            "source_column": source_column,
            "rows_total": int(len(df)),
            "rows_parseable_primary_autosomes": int(len(df)),
            "rows_pass_coverage": None,
            "schema": "gse60012_100bp_tile_matrix",
        }
        for sample_id, source_column in zip(sample_ids, raw_sample_labels, strict=True)
    ]
    meta = {
        "raw_sample_columns": len(raw_sample_labels),
        "metadata_mapping_status": "blocked_duplicate_non_gsm_tile_columns",
    }
    if coordinate_mapper is not None:
        meta.update(
            {
                "assembly_harmonization": "prefilter_liftover",
                "liftover_rows_input": int(len(df)),
                "liftover_rows_mapped": int(len(values)),
                "liftover_failure_counts": liftover_failure_counts,
            }
        )
    return values, stats, meta


def build_region_matrix(beta_df: pd.DataFrame, bin_size: int, min_sample_presence: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    parsed_rows = []
    keep_positions = []
    for row_number, cpg_id in enumerate(beta_df.index):
        chrom, pos_text = str(cpg_id).rsplit("_", 1)
        pos = int(pos_text)
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
    region_labels = pd.Series(cpg_meta["region_id"].values, index=beta_df.index)
    region_matrix = beta_df.iloc[keep_positions].groupby(region_labels, sort=False).mean().astype(np.float32)
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
    keep = stats["n_samples_present"] >= min_samples
    stats = stats.loc[keep].copy()
    stats["_chrom_sort"] = [chrom_sort_value(chrom) for chrom in stats["chrom"]]
    stats = stats.sort_values(["_chrom_sort", "start", "end"]).drop(columns=["_chrom_sort"])
    return region_matrix.loc[stats.index], stats.reset_index()


def write_blocker(out_dir: Path, dataset: str, payload: dict) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / f"{dataset}_conversion_blocker.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    lines = [
        f"# {dataset} Processed Methylation Conversion Blocker",
        "",
        f"- Reason: {payload.get('reason')}",
        f"- Input: `{payload.get('input_path')}`",
        "",
        "Training and CR validation were not run because a usable matrix was not created.",
    ]
    (out_dir / f"{dataset}_conversion_blocker.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def convert(input_path: Path, args: argparse.Namespace) -> dict:
    started = time.time()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    dataset = infer_dataset(input_path, args.dataset)
    schema = infer_schema(dataset, input_path, args.schema)
    if not input_path.exists():
        payload = {
            "status": "blocked",
            "reason": "input_supplement_missing",
            "input_path": str(input_path),
            "dataset": dataset,
            "schema": schema,
        }
        write_blocker(out_dir, dataset, payload)
        return payload

    coordinate_mapper = None
    if getattr(args, "liftover_chain", None):
        chain_path = Path(args.liftover_chain)
        if not chain_path.exists():
            payload = {
                "status": "blocked",
                "reason": "liftover_chain_missing",
                "input_path": str(input_path),
                "dataset": dataset,
                "schema": schema,
                "liftover_chain": str(chain_path),
            }
            write_blocker(out_dir, dataset, payload)
            return payload
        print(f"  Loading liftover chain {chain_path}")
        coordinate_mapper = LiftOverMapper(chain_path)

    sample_series = {}
    sample_stats = []
    schema_extra: dict = {}
    errors = []
    if schema == "gse60012_100bp_tile_matrix":
        try:
            beta_df, sample_stats, schema_extra = parse_gse60012_tile_matrix(
                input_path,
                max_rows=args.max_rows,
                coordinate_mapper=coordinate_mapper,
            )
            sample_series = {column: beta_df[column] for column in beta_df.columns}
        except Exception as exc:
            errors.append({"member": str(input_path), "error": str(exc)[:500]})
    else:
        for name, path, blob in iter_input_members(input_path, args.sample_limit):
            try:
                sample_id = sample_id_from_name(name)
                parser = parse_overlap_table if schema == "gse80672_overlap_percentage_coverage" else parse_bismark_cov_table
                if blob is not None:
                    with tempfile.TemporaryFile() as handle:
                        handle.write(blob)
                        handle.seek(0)
                        series, stats = parser(
                            sample_id=sample_id,
                            fileobj=handle,
                            min_coverage=args.min_coverage,
                            max_rows=args.max_rows,
                            coordinate_mapper=coordinate_mapper if parser is parse_bismark_cov_table else None,
                        )
                else:
                    series, stats = parser(
                        sample_id=sample_id,
                        path=path,
                        min_coverage=args.min_coverage,
                        max_rows=args.max_rows,
                        coordinate_mapper=coordinate_mapper if parser is parse_bismark_cov_table else None,
                    )
                if not series.empty:
                    sample_series[sample_id] = series
                sample_stats.append(stats)
                if len(sample_stats) % 25 == 0:
                    print(f"  Parsed {len(sample_stats)} samples")
            except Exception as exc:
                errors.append({"member": name, "error": str(exc)[:500]})

    if not sample_series:
        payload = {
            "status": "blocked",
            "reason": "no_parseable_sample_series",
            "input_path": str(input_path),
            "dataset": dataset,
            "schema": schema,
            "errors": errors[:20],
            "n_errors": len(errors),
        }
        write_blocker(out_dir, dataset, payload)
        return payload

    if schema != "gse60012_100bp_tile_matrix":
        beta_df = pd.DataFrame(sample_series).astype(np.float32)
    n_feature_rows_before_duplicate_collapse = int(beta_df.shape[0])
    n_duplicate_feature_rows_collapsed = int(beta_df.index.duplicated().sum())
    if n_duplicate_feature_rows_collapsed:
        beta_df = beta_df.groupby(beta_df.index, sort=False).mean().astype(np.float32)
    min_samples = int(math.ceil(beta_df.shape[1] * args.min_sample_presence))
    n_cpg_before = int(beta_df.shape[0])
    beta_df = beta_df.dropna(thresh=min_samples)
    beta_path = out_dir / f"{dataset}_beta_matrix.parquet"
    beta_df.to_parquet(beta_path, compression="zstd")

    cpg_stats = pd.DataFrame(
        {
            "cpg_id": beta_df.index,
            "n_samples_present": beta_df.notna().sum(axis=1).values,
            "mean_beta": beta_df.mean(axis=1).values,
            "std_beta": beta_df.std(axis=1).values,
        }
    )
    cpg_stats.to_csv(out_dir / f"{dataset}_cpg_coverage_stats.csv", index=False)

    region_matrix, region_stats = build_region_matrix(beta_df, args.bin_size, args.region_min_sample_presence)
    region_path = out_dir / f"{dataset}_region_matrix_5kb.parquet"
    region_stats_path = out_dir / f"{dataset}_region_stats_5kb.csv"
    region_matrix.to_parquet(region_path, compression="zstd")
    region_stats.to_csv(region_stats_path, index=False)

    metadata_path = Path(getattr(args, "metadata_path", META_FILE))
    meta = pd.read_csv(metadata_path)
    gse_meta = meta[meta["dataset_batch"] == dataset].copy()
    common = sorted(set(gse_meta["sample_id"]) & set(beta_df.columns))
    label_counts = gse_meta[gse_meta["sample_id"].isin(common)]["intervention"].value_counts().to_dict()
    status = "completed" if common else "completed_unmapped"
    liftover_failure_counts: dict[str, int] = {}
    liftover_rows_mapped = 0
    if coordinate_mapper is not None:
        for stats in sample_stats:
            failures = stats.get("liftover_failure_counts")
            if isinstance(failures, dict):
                for reason, count in failures.items():
                    liftover_failure_counts[reason] = liftover_failure_counts.get(reason, 0) + int(count)
            if "rows_liftover_mapped" in stats:
                liftover_rows_mapped += int(stats["rows_liftover_mapped"])

    manifest = {
        "status": status,
        "dataset": dataset,
        "input_path": str(input_path),
        "schema": schema,
        "min_coverage": int(args.min_coverage),
        "min_sample_presence": float(args.min_sample_presence),
        "region_bin_size": int(args.bin_size),
        "region_min_sample_presence": float(args.region_min_sample_presence),
        "n_samples_parsed": int(len(sample_series)),
        "n_samples_metadata_overlap": int(len(common)),
        "metadata_intervention_counts": label_counts,
        "metadata_age_known_overlap": int(gse_meta[gse_meta["sample_id"].isin(common)]["age_days"].notna().sum()),
        "n_feature_rows_before_duplicate_collapse": n_feature_rows_before_duplicate_collapse,
        "n_duplicate_feature_rows_collapsed": n_duplicate_feature_rows_collapsed,
        "n_cpg_before_presence_filter": n_cpg_before,
        "n_cpg_after_presence_filter": int(beta_df.shape[0]),
        "n_regions": int(region_matrix.shape[0]),
        "sex_mt_region_rows": int(region_stats["chrom"].isin(["chrX", "chrY", "chrM", "X", "Y", "M", "MT"]).sum()),
        "beta_matrix_path": str(beta_path),
        "region_matrix_path": str(region_path),
        "region_stats_path": str(region_stats_path),
        "metadata_path": str(metadata_path),
        "n_errors": len(errors),
        "errors": errors[:20],
        "exec_time_sec": round(time.time() - started, 1),
        **schema_extra,
    }
    if coordinate_mapper is not None:
        manifest.update(
            {
                "assembly_harmonization": "prefilter_liftover",
                "source_assembly": args.source_assembly,
                "target_assembly": args.target_assembly,
                "liftover_chain": str(Path(args.liftover_chain)),
                "liftover_chain_direction": "source_target_fields_to_query_fields",
                "liftover_rows_mapped": int(liftover_rows_mapped or schema_extra.get("liftover_rows_mapped", 0)),
                "liftover_failure_counts_total": liftover_failure_counts
                or schema_extra.get("liftover_failure_counts", {}),
                "liftover_unique_coordinate_cache_size": int(len(coordinate_mapper.cache)),
            }
        )
    (out_dir / f"{dataset}_matrix_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    if dataset == "GSE80672":
        (out_dir / "matrix_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    pd.DataFrame(sample_stats).to_csv(out_dir / f"{dataset}_sample_parse_stats.csv", index=False)
    if not common:
        write_blocker(
            out_dir,
            dataset,
            {
                "status": "blocked",
                "reason": "metadata_sample_mapping_unresolved",
                "input_path": str(input_path),
                "dataset": dataset,
                "schema": schema,
                "n_samples_parsed": int(len(sample_series)),
                "n_samples_metadata_overlap": 0,
                **schema_extra,
            },
        )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="GEO processed supplement file or tar archive.")
    parser.add_argument("--dataset", default=None, help="GSE accession. Inferred from input path when omitted.")
    parser.add_argument(
        "--schema",
        default="auto",
        choices=[
            "auto",
            "gse80672_overlap_percentage_coverage",
            "bismark_cov_per_sample_tar",
            "gse60012_100bp_tile_matrix",
        ],
    )
    parser.add_argument("--min_coverage", type=int, default=5)
    parser.add_argument("--min_sample_presence", type=float, default=0.5)
    parser.add_argument("--region_min_sample_presence", type=float, default=0.8)
    parser.add_argument("--bin_size", type=int, default=5000)
    parser.add_argument("--out_dir", default=str(MULTI_DIR))
    parser.add_argument("--metadata_path", default=str(META_FILE))
    parser.add_argument("--liftover_chain", default=None, help="Optional UCSC chain file for pre-filter coordinate liftover.")
    parser.add_argument("--source_assembly", default=None)
    parser.add_argument("--target_assembly", default=None)
    parser.add_argument("--sample_limit", type=int, default=None)
    parser.add_argument("--max_rows", type=int, default=None, help="Parser smoke limit per sample.")
    args = parser.parse_args()

    manifest = convert(Path(args.input), args)
    print(f"[Convert] status={manifest['status']}")
    print(json.dumps(manifest, indent=2)[:4000])


if __name__ == "__main__":
    main()

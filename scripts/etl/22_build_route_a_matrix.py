#!/usr/bin/env python3
"""Build Route A beta/region matrices from staged processed methylation files.

Supported local processed schemas:
- bismark_cov_6col: per-sample Bismark coverage files
- cpg_beta_table: CpG beta table, per-sample long/wide or per-sample file
- region_beta_matrix: 5kb region beta matrix with region_id plus sample columns

This script reads local processed methylation files only. It does not download
data, run Bismark, train models, or authorize benchmarking.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import re
import time
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")
DEFAULT_REFERENCE = ROOT / "results" / "multidataset_v8_2_prefilter_liftover" / "all_rrbs_region_matrix_5kb.parquet"
DEFAULT_OUT_BASE = ROOT / "results" / "route_a_intake"

PLACEHOLDER_TOKENS = {"", "to_be_filled", "tbd", "na", "n/a", "placeholder"}
REMOTE_PREFIXES = ("http://", "https://", "ftp://", "s3://", "gs://")
EXCLUDE_CHROMS = {"chrx", "chry", "chrm", "x", "y", "m", "mt"}
CHROM_COLUMNS = ["chr", "chrom", "chromosome", "seqnames"]
POS_COLUMNS = ["pos", "position", "start", "cpg_pos", "coordinate"]
BETA_COLUMNS = ["beta", "methylation", "methylation_level", "meth_percent", "percent_methylation"]
REGION_RE = re.compile(r"^(chr)?(?P<chrom>[0-9]+|X|Y|M|MT):(?P<start>\d+)-(?P<end>\d+)$", re.I)
SHA256_RE = re.compile(r"^[0-9a-f]{64}$", re.I)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str), encoding="utf-8")


def is_placeholder(value: Any) -> bool:
    return str(value or "").strip().lower() in PLACEHOLDER_TOKENS


def is_remote(value: Any) -> bool:
    return str(value or "").strip().lower().startswith(REMOTE_PREFIXES)


def normalize_schema(value: Any) -> str:
    text = str(value or "").strip().lower().replace("-", "_").replace(" ", "_")
    aliases = {
        "bismark": "bismark_cov_6col",
        "bismark_cov": "bismark_cov_6col",
        "bismark_cov_6col": "bismark_cov_6col",
        "cpg_beta": "cpg_beta_table",
        "cpg_beta_table": "cpg_beta_table",
        "beta_table": "cpg_beta_table",
        "region_beta": "region_beta_matrix",
        "region_beta_matrix": "region_beta_matrix",
        "5kb_region_beta_matrix": "region_beta_matrix",
    }
    return aliases.get(text, text)


def normalize_chrom(value: Any) -> str | None:
    text = str(value or "").strip()
    if not text:
        return None
    lower = text.lower()
    if lower in EXCLUDE_CHROMS:
        return None
    if lower.startswith("chr"):
        chrom = "chr" + text[3:]
    elif text.isdigit():
        chrom = "chr" + text
    else:
        return None
    if chrom.lower() in EXCLUDE_CHROMS:
        return None
    return chrom


def cpg_id(chrom: Any, pos: Any) -> str | None:
    norm = normalize_chrom(chrom)
    if norm is None:
        return None
    try:
        pos_int = int(float(str(pos)))
    except ValueError:
        return None
    return f"{norm}_{pos_int}"


def region_from_cpg_id(value: str, bin_size: int) -> tuple[str, str, int, int] | None:
    try:
        chrom, pos_text = str(value).rsplit("_", 1)
        chrom = normalize_chrom(chrom) or ""
        pos = int(float(pos_text))
    except ValueError:
        return None
    if not chrom:
        return None
    start = (pos // bin_size) * bin_size
    end = start + bin_size - 1
    return f"{chrom}:{start}-{end}", chrom, start, end


def parse_region_id(value: Any) -> tuple[str, str, int, int] | None:
    text = str(value or "").strip()
    match = REGION_RE.match(text)
    if not match:
        return None
    chrom = normalize_chrom(match.group("chrom"))
    if chrom is None:
        return None
    start = int(match.group("start"))
    end = int(match.group("end"))
    return f"{chrom}:{start}-{end}", chrom, start, end


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if str(path).endswith(".gz") else path.open("rt", encoding="utf-8", errors="replace")


def delimiter_for(path: Path) -> str:
    with open_text(path) as handle:
        sample = handle.read(4096)
    return "\t" if sample.count("\t") >= sample.count(",") else ","


def resolve_local_path(value: str) -> Path:
    if is_placeholder(value):
        raise ValueError("path_or_uri_is_placeholder")
    if is_remote(value):
        raise ValueError("remote_uri_not_supported")
    path = Path(str(value))
    return path if path.is_absolute() else ROOT / path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def metadata_metrics(sample_rows: list[dict[str, str]]) -> dict[str, Any]:
    n = len(sample_rows)
    age_known = 0
    target_tissue_rows = 0
    noncontrol = 0
    tissues: dict[str, int] = defaultdict(int)
    for row in sample_rows:
        if str(row.get("age_days", "")).strip() and str(row.get("age_weeks", "")).strip():
            age_known += 1
        tissue = str(row.get("tissue", "")).strip()
        tissues[tissue] += 1
        if tissue in {"brain_cortex", "heart", "lung"}:
            target_tissue_rows += 1
        if str(row.get("intervention", "control")).strip().lower() != "control":
            noncontrol += 1
    return {
        "n_sample_sheet_rows": n,
        "age_coverage": age_known / n if n else 0.0,
        "target_tissue_fraction": target_tissue_rows / n if n else 0.0,
        "tissue_counts": dict(sorted(tissues.items())),
        "noncontrol_samples": noncontrol,
        "non_headline_stratification_risk": bool(noncontrol),
    }


def file_manifest_rows(file_manifest: Path) -> list[dict[str, str]]:
    rows = read_csv(file_manifest)
    processed = [
        row
        for row in rows
        if "processed" in str(row.get("file_role", "")).lower()
        or "coverage" in str(row.get("file_role", "")).lower()
        or "beta" in str(row.get("file_role", "")).lower()
    ]
    return processed


def validate_paths(rows: list[dict[str, str]]) -> list[dict[str, Any]]:
    issues = []
    actual_cache: dict[str, dict[str, Any]] = {}
    for idx, row in enumerate(rows, start=2):
        path_text = row.get("path_or_uri", "")
        try:
            path = resolve_local_path(path_text)
        except ValueError as exc:
            issues.append({"row": idx, "sample_id": row.get("sample_id", ""), "issue": str(exc), "path_or_uri": path_text})
            continue
        if not path.exists():
            issues.append({"row": idx, "sample_id": row.get("sample_id", ""), "issue": "local_path_not_found", "path_or_uri": str(path)})
            continue

        path_key = str(path)
        if path_key not in actual_cache:
            actual_cache[path_key] = {
                "file_size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        actual = actual_cache[path_key]

        expected_size_text = str(row.get("file_size_bytes", "")).strip()
        try:
            expected_size = int(float(expected_size_text))
        except ValueError:
            issues.append({"row": idx, "sample_id": row.get("sample_id", ""), "issue": "file_size_bytes_invalid", "path_or_uri": str(path), "expected": expected_size_text})
            continue
        if expected_size != int(actual["file_size_bytes"]):
            issues.append(
                {
                    "row": idx,
                    "sample_id": row.get("sample_id", ""),
                    "issue": "file_size_bytes_mismatch",
                    "path_or_uri": str(path),
                    "expected": expected_size,
                    "actual": int(actual["file_size_bytes"]),
                }
            )

        expected_sha = str(row.get("sha256", "")).strip().lower()
        if not SHA256_RE.match(expected_sha):
            issues.append({"row": idx, "sample_id": row.get("sample_id", ""), "issue": "file_sha256_invalid_or_placeholder", "path_or_uri": str(path), "expected": expected_sha})
            continue
        if expected_sha != str(actual["sha256"]).lower():
            issues.append(
                {
                    "row": idx,
                    "sample_id": row.get("sample_id", ""),
                    "issue": "file_sha256_mismatch",
                    "path_or_uri": str(path),
                    "expected": expected_sha,
                    "actual": actual["sha256"],
                }
            )
    return issues


def parse_bismark_cov(path: Path, sample_id: str, min_coverage: int) -> tuple[pd.Series, dict[str, Any]]:
    records: dict[str, float] = {}
    stats = {
        "sample_id": sample_id,
        "path": str(path),
        "schema": "bismark_cov_6col",
        "rows_total": 0,
        "rows_parseable_primary_autosomes": 0,
        "rows_pass_coverage": 0,
    }
    with open_text(path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                parts = line.rstrip("\n").split(",")
            if len(parts) < 6:
                continue
            stats["rows_total"] += 1
            cid = cpg_id(parts[0], parts[1])
            if cid is None:
                continue
            try:
                meth = float(parts[4])
                unmeth = float(parts[5])
            except ValueError:
                continue
            stats["rows_parseable_primary_autosomes"] += 1
            total = meth + unmeth
            if total < min_coverage or total <= 0:
                continue
            stats["rows_pass_coverage"] += 1
            records[cid] = meth / total
    return pd.Series(records, name=sample_id, dtype=np.float32), stats


def load_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep=delimiter_for(path), compression="infer", low_memory=False)


def lower_column_map(columns: Iterable[Any]) -> dict[str, str]:
    return {str(col).strip().lower(): str(col) for col in columns}


def first_present(columns: dict[str, str], candidates: list[str]) -> str | None:
    for candidate in candidates:
        if candidate in columns:
            return columns[candidate]
    return None


def beta_from_value(value: Any) -> float | None:
    try:
        beta = float(value)
    except (TypeError, ValueError):
        return None
    if math.isnan(beta):
        return None
    return beta


def parse_cpg_beta_table(path: Path, manifest_sample_id: str | None, sample_ids: set[str]) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    df = load_table(path)
    colmap = lower_column_map(df.columns)
    sample_col = colmap.get("sample_id")
    chrom_col = first_present(colmap, CHROM_COLUMNS)
    pos_col = first_present(colmap, POS_COLUMNS)
    beta_col = first_present(colmap, BETA_COLUMNS)
    stats: list[dict[str, Any]] = []
    if chrom_col is None or pos_col is None:
        raise ValueError(f"cpg_beta_table_missing_coordinate_columns:{path}")
    ids = [cpg_id(chrom, pos) for chrom, pos in zip(df[chrom_col], df[pos_col], strict=False)]
    keep = [idx for idx, value in enumerate(ids) if value is not None]
    if not keep:
        raise ValueError(f"cpg_beta_table_no_autosomal_coordinates:{path}")
    cpg_ids = [str(ids[idx]) for idx in keep]
    if sample_col and beta_col:
        long = df.iloc[keep][[sample_col, beta_col]].copy()
        long["_cpg_id"] = cpg_ids
        long[beta_col] = long[beta_col].map(beta_from_value)
        long = long.dropna(subset=[beta_col])
        long = long[long[sample_col].astype(str).isin(sample_ids)]
        matrix = long.pivot_table(index="_cpg_id", columns=sample_col, values=beta_col, aggfunc="mean").astype(np.float32)
        for sample_id in matrix.columns:
            stats.append({"sample_id": sample_id, "path": str(path), "schema": "cpg_beta_table_long", "rows_pass_coverage": int(matrix[sample_id].notna().sum())})
        return matrix, stats
    if beta_col:
        sample_id = manifest_sample_id
        if not sample_id:
            raise ValueError(f"cpg_beta_table_per_sample_requires_manifest_sample_id:{path}")
        values = df.iloc[keep][beta_col].map(beta_from_value)
        series = pd.Series(values.to_numpy(), index=cpg_ids, name=sample_id, dtype=np.float32).dropna()
        stats.append({"sample_id": sample_id, "path": str(path), "schema": "cpg_beta_table_per_sample", "rows_pass_coverage": int(series.notna().sum())})
        return series.to_frame(), stats

    blocked = {chrom_col, pos_col}
    candidate_cols = [col for col in df.columns if col not in blocked and str(col) in sample_ids]
    if not candidate_cols:
        raise ValueError(f"cpg_beta_table_no_beta_columns_or_matching_sample_columns:{path}")
    values = df.iloc[keep][candidate_cols].apply(lambda col: col.map(beta_from_value)).astype(np.float32)
    values.index = cpg_ids
    for sample_id in values.columns:
        stats.append({"sample_id": sample_id, "path": str(path), "schema": "cpg_beta_table_wide", "rows_pass_coverage": int(values[sample_id].notna().sum())})
    return values, stats


def parse_region_beta_matrix(path: Path, sample_ids: set[str], manifest_sample_id: str | None = None) -> tuple[pd.DataFrame, pd.DataFrame, list[dict[str, Any]]]:
    df = load_table(path)
    colmap = lower_column_map(df.columns)
    region_col = colmap.get("region_id")
    beta_col = first_present(colmap, BETA_COLUMNS)
    stats: list[dict[str, Any]] = []
    if region_col is None:
        chrom_col = first_present(colmap, CHROM_COLUMNS)
        start_col = colmap.get("start")
        end_col = colmap.get("end")
        if chrom_col is None or start_col is None or end_col is None:
            raise ValueError(f"region_beta_matrix_missing_region_or_coordinate_columns:{path}")
        regions = []
        for chrom, start, end in zip(df[chrom_col], df[start_col], df[end_col], strict=False):
            norm = normalize_chrom(chrom)
            if norm is None:
                regions.append(None)
                continue
            try:
                regions.append(f"{norm}:{int(float(start))}-{int(float(end))}")
            except ValueError:
                regions.append(None)
        region_ids = regions
    else:
        region_ids = []
        for value in df[region_col]:
            parsed = parse_region_id(value)
            region_ids.append(parsed[0] if parsed else None)

    keep = [idx for idx, value in enumerate(region_ids) if value is not None]
    if not keep:
        raise ValueError(f"region_beta_matrix_no_autosomal_regions:{path}")
    kept_region_ids = [str(region_ids[idx]) for idx in keep]
    if beta_col:
        if not manifest_sample_id:
            raise ValueError(f"region_beta_matrix_per_sample_requires_manifest_sample_id:{path}")
        values = df.iloc[keep][beta_col].map(beta_from_value)
        matrix = pd.Series(values.to_numpy(), index=kept_region_ids, name=manifest_sample_id, dtype=np.float32).dropna().to_frame()
    else:
        blocked_cols = {region_col, "chr", "chrom", "chromosome", "seqnames", "start", "end"}
        candidate_cols = [col for col in df.columns if str(col) in sample_ids and col not in blocked_cols]
        if not candidate_cols:
            raise ValueError(f"region_beta_matrix_no_matching_sample_columns:{path}")
        matrix = df.iloc[keep][candidate_cols].apply(lambda col: col.map(beta_from_value)).astype(np.float32)
        matrix.index = kept_region_ids
    matrix = matrix.groupby(matrix.index, sort=False).mean().astype(np.float32)
    region_stats = region_stats_from_region_matrix(matrix)
    for sample_id in matrix.columns:
        stats.append({"sample_id": sample_id, "path": str(path), "schema": "region_beta_matrix", "rows_pass_coverage": int(matrix[sample_id].notna().sum())})
    return matrix, region_stats, stats


def build_region_matrix(beta_df: pd.DataFrame, bin_size: int, min_sample_presence: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    parsed = []
    keep_index = []
    for cpg in beta_df.index:
        region = region_from_cpg_id(str(cpg), bin_size)
        if region is None:
            continue
        region_id, chrom, start, end = region
        parsed.append({"cpg_id": str(cpg), "region_id": region_id, "chrom": chrom, "start": start, "end": end})
        keep_index.append(cpg)
    if not parsed:
        raise ValueError("no_cpg_rows_can_be_mapped_to_regions")
    meta = pd.DataFrame(parsed)
    labels = pd.Series(meta["region_id"].values, index=keep_index)
    region_matrix = beta_df.loc[keep_index].groupby(labels, sort=False).mean().astype(np.float32)
    stats = region_stats_from_region_matrix(region_matrix, cpg_meta=meta)
    min_samples = int(math.ceil(region_matrix.shape[1] * min_sample_presence))
    keep = stats["n_samples_present"] >= min_samples
    stats = stats.loc[keep].copy()
    return region_matrix.loc[stats.index], stats.reset_index()


def region_stats_from_region_matrix(region_matrix: pd.DataFrame, cpg_meta: pd.DataFrame | None = None) -> pd.DataFrame:
    rows = []
    for region_id in region_matrix.index:
        parsed = parse_region_id(region_id)
        if parsed is None:
            continue
        _, chrom, start, end = parsed
        rows.append({"region_id": region_id, "chrom": chrom, "start": start, "end": end})
    coords = pd.DataFrame(rows).drop_duplicates("region_id").set_index("region_id")
    if cpg_meta is not None:
        counts = cpg_meta.groupby("region_id", sort=False).size().rename("n_cpgs")
    else:
        counts = pd.Series(1, index=region_matrix.index, name="n_cpgs")
    stats = coords.join(counts)
    stats["n_samples_present"] = region_matrix.notna().sum(axis=1).astype(int)
    stats["mean_beta"] = region_matrix.mean(axis=1, skipna=True).astype(float)
    stats["std_beta"] = region_matrix.std(axis=1, skipna=True).astype(float)
    stats["_chrom_sort"] = [chrom_sort_value(chrom) for chrom in stats["chrom"]]
    stats = stats.sort_values(["_chrom_sort", "start", "end"]).drop(columns=["_chrom_sort"])
    return stats


def chrom_sort_value(chrom: str) -> tuple[int, str]:
    clean = str(chrom).removeprefix("chr")
    if clean.isdigit():
        return int(clean), ""
    return 10_000, clean


def beta_bounds(df: pd.DataFrame) -> tuple[float | None, float | None, int]:
    if df.empty:
        return None, None, 0
    numeric = df.apply(pd.to_numeric, errors="coerce")
    min_value = float(numeric.min(skipna=True).min(skipna=True))
    max_value = float(numeric.max(skipna=True).max(skipna=True))
    bad = int(((numeric < 0) | (numeric > 1)).sum(skipna=True).sum())
    return min_value, max_value, bad


def merge_matrices(matrices: list[pd.DataFrame]) -> pd.DataFrame:
    if not matrices:
        raise ValueError("no_matrices_to_merge")
    merged = pd.concat(matrices, axis=1)
    merged = merged.loc[:, ~merged.columns.duplicated()]
    return merged.groupby(merged.index, sort=False).mean().astype(np.float32)


def load_reference_regions(path: Path) -> set[str]:
    if not path.exists():
        return set()
    ref = pd.read_parquet(path, columns=[])
    return set(ref.index.astype(str))


def write_blocker(out_dir: Path, payload: dict[str, Any]) -> None:
    write_json(out_dir / "route_a_matrix_blocker.json", payload)
    lines = [
        "# Route A Matrix Build Blocker",
        "",
        f"- Reason: {payload.get('reason')}",
        f"- Status: `{payload.get('status')}`",
        "",
        "No training, Bismark, download, or autoresearch was run.",
    ]
    (out_dir / "route_a_matrix_blocker.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def build(args: argparse.Namespace) -> dict[str, Any]:
    started = time.time()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sample_rows = read_csv(args.sample_sheet)
    sample_ids = {row["sample_id"] for row in sample_rows}
    processed_rows = file_manifest_rows(args.file_manifest)
    path_issues = validate_paths(processed_rows)
    if path_issues:
        payload = {
            "status": "blocked",
            "reason": "local_processed_path_gate_failed",
            "path_issues": path_issues[:50],
            "n_path_issues": len(path_issues),
            "file_integrity_checked": True,
            "file_integrity_checks": "local_path_exists; file_size_bytes; sha256",
        }
        write_blocker(out_dir, payload)
        return payload

    by_schema: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in processed_rows:
        schema = normalize_schema(row.get("processed_schema"))
        by_schema[schema].append(row)
    supported = {"bismark_cov_6col", "cpg_beta_table", "region_beta_matrix"}
    unsupported = sorted(set(by_schema) - supported)
    if unsupported:
        payload = {"status": "blocked", "reason": "unsupported_processed_schema", "unsupported_schemas": unsupported}
        write_blocker(out_dir, payload)
        return payload

    cpg_matrices: list[pd.DataFrame] = []
    region_matrices: list[pd.DataFrame] = []
    sample_stats: list[dict[str, Any]] = []
    errors: list[dict[str, Any]] = []
    seen_path_schema: set[tuple[str, str]] = set()
    for schema, rows in by_schema.items():
        for row in rows:
            sample_id = row.get("sample_id", "")
            path = resolve_local_path(row.get("path_or_uri", ""))
            key = (str(path), schema)
            if schema in {"cpg_beta_table", "region_beta_matrix"} and key in seen_path_schema:
                continue
            seen_path_schema.add(key)
            try:
                if schema == "bismark_cov_6col":
                    if sample_id not in sample_ids:
                        raise ValueError(f"manifest_sample_id_not_in_sample_sheet:{sample_id}")
                    series, stats = parse_bismark_cov(path, sample_id, args.min_coverage)
                    cpg_matrices.append(series.to_frame())
                    sample_stats.append(stats)
                elif schema == "cpg_beta_table":
                    matrix, stats = parse_cpg_beta_table(path, sample_id if sample_id in sample_ids else None, sample_ids)
                    cpg_matrices.append(matrix)
                    sample_stats.extend(stats)
                elif schema == "region_beta_matrix":
                    matrix, _, stats = parse_region_beta_matrix(path, sample_ids, sample_id if sample_id in sample_ids else None)
                    region_matrices.append(matrix)
                    sample_stats.extend(stats)
            except Exception as exc:  # noqa: BLE001 - all schema blockers are reported in manifest
                errors.append({"sample_id": sample_id, "path": str(path), "schema": schema, "error": str(exc)[:500]})

    if errors:
        payload = {"status": "blocked", "reason": "processed_schema_parse_errors", "errors": errors[:50], "n_errors": len(errors)}
        write_blocker(out_dir, payload)
        return payload
    if not cpg_matrices and not region_matrices:
        payload = {"status": "blocked", "reason": "no_parseable_processed_methylation"}
        write_blocker(out_dir, payload)
        return payload

    beta_df = merge_matrices(cpg_matrices) if cpg_matrices else pd.DataFrame()
    if not beta_df.empty:
        min_samples = int(math.ceil(beta_df.shape[1] * args.min_sample_presence))
        beta_df = beta_df.dropna(thresh=min_samples)
        region_from_cpg, region_stats = build_region_matrix(beta_df, args.bin_size, args.region_min_sample_presence)
        region_matrices.append(region_from_cpg)
    region_matrix = merge_matrices(region_matrices)
    region_stats = region_stats_from_region_matrix(region_matrix).reset_index()
    min_region_samples = int(math.ceil(region_matrix.shape[1] * args.region_min_sample_presence))
    keep_regions = region_stats["n_samples_present"] >= min_region_samples
    region_stats = region_stats.loc[keep_regions].copy()
    region_matrix = region_matrix.loc[region_stats["region_id"].astype(str)]

    beta_path = out_dir / "route_a_beta_matrix.parquet"
    region_path = out_dir / "route_a_region_matrix_5kb.parquet"
    region_stats_path = out_dir / "route_a_region_stats_5kb.csv"
    if beta_df.empty:
        region_matrix.to_parquet(beta_path, compression="zstd")
    else:
        beta_df.to_parquet(beta_path, compression="zstd")
    region_matrix.to_parquet(region_path, compression="zstd")
    region_stats.to_csv(region_stats_path, index=False)
    pd.DataFrame(sample_stats).to_csv(out_dir / "route_a_sample_parse_stats.csv", index=False)

    metadata = metadata_metrics(sample_rows)
    matrix_sample_ids = set(region_matrix.columns.astype(str))
    metadata_overlap_count = len(sample_ids & matrix_sample_ids)
    metadata_overlap = metadata_overlap_count / len(sample_ids) if sample_ids else 0.0
    beta_min, beta_max, beta_bad = beta_bounds(region_matrix)
    reference_regions = load_reference_regions(args.reference_matrix)
    common_regions = sorted(set(region_matrix.index.astype(str)) & reference_regions) if reference_regions else []
    pd.DataFrame({"region_id": common_regions}).to_csv(out_dir / "route_a_common_regions_with_reference.csv", index=False)

    status = "completed"
    blockers = []
    if metadata_overlap < 0.95:
        blockers.append("metadata_overlap_below_95_percent")
    if metadata["age_coverage"] < 0.95:
        blockers.append("age_coverage_below_95_percent")
    if len(common_regions) < args.min_common_regions:
        blockers.append("common_regions_below_threshold")
    if beta_bad or beta_min is None or beta_max is None or beta_min < 0 or beta_max > 1:
        blockers.append("beta_range_failed")
    if blockers:
        status = "blocked"

    manifest = {
        "status": status,
        "submission_id": args.submission_id,
        "dataset_batch": f"ROUTEA_{args.submission_id}",
        "sample_sheet": str(args.sample_sheet),
        "file_manifest": str(args.file_manifest),
        "reference_matrix": str(args.reference_matrix),
        "supported_schemas": sorted(supported),
        "schemas_seen": sorted(by_schema),
        "file_integrity_checked": True,
        "file_integrity_checks": "local_path_exists; file_size_bytes; sha256",
        "min_coverage": int(args.min_coverage),
        "min_sample_presence": float(args.min_sample_presence),
        "region_min_sample_presence": float(args.region_min_sample_presence),
        "region_bin_size": int(args.bin_size),
        "n_sample_sheet_rows": len(sample_rows),
        "n_matrix_samples": int(region_matrix.shape[1]),
        "metadata_overlap_count": int(metadata_overlap_count),
        "metadata_overlap": float(metadata_overlap),
        "age_coverage": float(metadata["age_coverage"]),
        "target_tissue_fraction": float(metadata["target_tissue_fraction"]),
        "noncontrol_samples": int(metadata["noncontrol_samples"]),
        "non_headline_stratification_risk": bool(metadata["non_headline_stratification_risk"]),
        "tissue_counts": metadata["tissue_counts"],
        "n_cpg_rows": int(beta_df.shape[0]) if not beta_df.empty else 0,
        "n_regions": int(region_matrix.shape[0]),
        "common_regions": int(len(common_regions)),
        "common_regions_with_reference": int(len(common_regions)),
        "min_common_regions_required": int(args.min_common_regions),
        "beta_min": beta_min,
        "beta_max": beta_max,
        "beta_values_outside_0_1": int(beta_bad),
        "blockers": blockers,
        "beta_matrix_path": str(beta_path),
        "region_matrix_path": str(region_path),
        "region_stats_path": str(region_stats_path),
        "common_regions_path": str(out_dir / "route_a_common_regions_with_reference.csv"),
        "sample_parse_stats_path": str(out_dir / "route_a_sample_parse_stats.csv"),
        "exec_time_sec": round(time.time() - started, 3),
        "download_authorized": False,
        "bismark_authorized": False,
        "training_authorized": False,
        "autoresearch_authorized": False,
        "created_at": datetime.now().isoformat(timespec="seconds"),
    }
    write_json(out_dir / "route_a_matrix_manifest.json", manifest)
    blocker_json = out_dir / "route_a_matrix_blocker.json"
    blocker_md = out_dir / "route_a_matrix_blocker.md"
    if status == "completed":
        blocker_json.unlink(missing_ok=True)
        blocker_md.unlink(missing_ok=True)
    else:
        write_blocker(out_dir, manifest | {"reason": ";".join(blockers) or "unknown_matrix_gate_blocker"})
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-sheet", required=True, type=Path)
    parser.add_argument("--file-manifest", required=True, type=Path)
    parser.add_argument("--submission-id", required=True)
    parser.add_argument("--reference-matrix", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--min-coverage", type=int, default=5)
    parser.add_argument("--min-sample-presence", type=float, default=0.5)
    parser.add_argument("--region-min-sample-presence", type=float, default=0.8)
    parser.add_argument("--bin-size", type=int, default=5000)
    parser.add_argument("--min-common-regions", type=int, default=50000)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.out_dir is None:
        args.out_dir = DEFAULT_OUT_BASE / args.submission_id
    manifest = build(args)
    print(json.dumps(manifest, indent=2, sort_keys=True, default=str)[:6000])
    raise SystemExit(0 if manifest.get("status") == "completed" else 2)


if __name__ == "__main__":
    main()

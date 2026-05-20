#!/usr/bin/env python3
"""Run the approved GSE83947 Route B minimal Bismark pilot and matrix gate."""
from __future__ import annotations

import argparse
import gzip
import json
import math
import os
import subprocess
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")
V14_DIR = ROOT / "results" / "ralph_v14_public_data_rescue"
RAW_DIR = V14_DIR / "gse83947_raw_fastq_pilot"
OUT_DIR = V14_DIR / "gse83947_bismark_pilot"
ENV_DIR = ROOT / "tools" / "route_b_bismark_env"
REF_DIR = ROOT / "references" / "GRCm38_ensembl102"
REFERENCE_MATRIX = ROOT / "results" / "multidataset_v8_2_prefilter_liftover" / "all_rrbs_region_matrix_5kb.parquet"
COMMON_REGION_GATE = 50_000


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str), encoding="utf-8")


def env() -> dict[str, str]:
    out = os.environ.copy()
    out["PATH"] = f"{ENV_DIR / 'bin'}:{out.get('PATH', '')}"
    return out


def run(cmd: list[str], cwd: Path, log_prefix: Path) -> subprocess.CompletedProcess[str]:
    log_prefix.parent.mkdir(parents=True, exist_ok=True)
    completed = subprocess.run(cmd, cwd=str(cwd), env=env(), text=True, capture_output=True, check=False)
    log_prefix.with_suffix(".stdout.log").write_text(completed.stdout, encoding="utf-8")
    log_prefix.with_suffix(".stderr.log").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(f"Command failed ({completed.returncode}): {' '.join(cmd)}\n{completed.stderr[-2000:]}")
    return completed


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


def region_id(chrom: str, pos: int, bin_size: int = 5000) -> str:
    start = (int(pos) // bin_size) * bin_size
    return f"{chrom}:{start}-{start + bin_size - 1}"


def iter_lines(path: Path) -> Iterable[str]:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", errors="replace") as handle:
        for line in handle:
            yield line.rstrip("\n")


def parse_cov_to_regions(path: Path, sample_id: str, min_coverage: int) -> tuple[pd.Series, dict]:
    beta_sum: dict[str, float] = defaultdict(float)
    beta_n: dict[str, int] = defaultdict(int)
    coverage_values: list[float] = []
    stats = {
        "sample_id": sample_id,
        "coverage_path": str(path),
        "rows_total": 0,
        "rows_parseable_primary_autosomes": 0,
        "rows_pass_coverage_ge5": 0,
        "beta_min": None,
        "beta_max": None,
        "schema": "bismark_cov_6col",
    }
    for line in iter_lines(path):
        if not line or line.startswith("track") or line.startswith("#"):
            continue
        stats["rows_total"] += 1
        parts = line.split("\t")
        if len(parts) < 6:
            continue
        chrom = normalize_chrom(parts[0])
        if chrom is None:
            continue
        try:
            pos = int(float(parts[1]))
            pct = float(parts[3])
            methylated = float(parts[4])
            unmethylated = float(parts[5])
        except ValueError:
            continue
        total = methylated + unmethylated
        if total <= 0:
            continue
        beta = methylated / total
        if not (0.0 <= beta <= 1.0):
            if 0.0 <= pct <= 100.0:
                beta = pct / 100.0
            else:
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
    series = pd.Series({rid: beta_sum[rid] / beta_n[rid] for rid in beta_sum}, dtype=np.float32, name=sample_id)
    stats["n_regions"] = int(len(series))
    stats["median_coverage"] = None if not coverage_values else round(float(np.median(coverage_values)), 3)
    stats["beta_min"] = None if stats["beta_min"] is None else round(float(stats["beta_min"]), 6)
    stats["beta_max"] = None if stats["beta_max"] is None else round(float(stats["beta_max"]), 6)
    return series, stats


def reference_regions() -> set[str]:
    ref = pd.read_parquet(REFERENCE_MATRIX, columns=[])
    return set(map(str, ref.index))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download-manifest", default=str(RAW_DIR / "raw_fastq_download_manifest.csv"))
    parser.add_argument("--sample-manifest", default=str(V14_DIR / "gse83947_raw_pilot_preflight" / "route_b_gse83947_raw_pilot_manifest.csv"))
    parser.add_argument("--out-dir", default=str(OUT_DIR))
    parser.add_argument("--parallel", type=int, default=4)
    parser.add_argument("--extract-parallel", type=int, default=4)
    parser.add_argument("--min-coverage", type=int, default=5)
    parser.add_argument("--authorize-bismark", action="store_true")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if not args.authorize_bismark:
        state = {
            "timestamp": utc_now(),
            "status": "blocked",
            "reason": "bismark_not_authorized",
            "training_authorized": False,
            "autoresearch_authorized": False,
        }
        write_json(out_dir / "gse83947_bismark_pilot_state.json", state)
        print(json.dumps(state, indent=2))
        raise SystemExit(2)

    download = pd.read_csv(args.download_manifest)
    samples = pd.read_csv(args.sample_manifest)
    command_rows = []
    cov_paths: dict[str, Path] = {}
    alignments: list[dict] = []
    for sample_id, sample_rows in download.groupby("sample_id"):
        mates = sample_rows.sort_values("mate")
        r1 = Path(str(mates[mates["mate"].eq(1)].iloc[0]["local_path"]))
        r2 = Path(str(mates[mates["mate"].eq(2)].iloc[0]["local_path"]))
        sample_out = out_dir / "bismark" / sample_id
        sample_out.mkdir(parents=True, exist_ok=True)
        bam_candidates = list(sample_out.glob("*.bam"))
        if not bam_candidates:
            cmd = [
                str(ENV_DIR / "bin" / "bismark"),
                "--genome_folder",
                str(REF_DIR),
                "--bowtie2",
                "--path_to_bowtie2",
                str(ENV_DIR / "bin"),
                "--parallel",
                str(args.parallel),
                "-o",
                str(sample_out),
                "-1",
                str(r1),
                "-2",
                str(r2),
            ]
            run(cmd, cwd=ROOT, log_prefix=sample_out / "01_bismark_alignment")
            command_rows.append({"sample_id": sample_id, "step": "bismark_alignment", "command": " ".join(cmd)})
            bam_candidates = list(sample_out.glob("*.bam"))
        else:
            command_rows.append(
                {
                    "sample_id": sample_id,
                    "step": "bismark_alignment_reused",
                    "command": f"existing_bam={bam_candidates[0]}",
                }
            )
        if not bam_candidates:
            raise RuntimeError(f"No BAM output found for {sample_id}")
        bam = bam_candidates[0]
        cov_candidates = [p for p in sample_out.glob("*.bismark.cov.gz") if p.stat().st_size > 0]
        if not cov_candidates:
            cmd = [
                str(ENV_DIR / "bin" / "bismark_methylation_extractor"),
                "--paired-end",
                "--gzip",
                "--bedGraph",
                "--parallel",
                str(args.extract_parallel),
                "--buffer_size",
                "4G",
                "--output",
                str(sample_out),
                str(bam),
            ]
            run(cmd, cwd=ROOT, log_prefix=sample_out / "02_bismark_methylation_extractor")
            command_rows.append({"sample_id": sample_id, "step": "methylation_extractor", "command": " ".join(cmd)})
            cov_candidates = [p for p in sample_out.glob("*.bismark.cov.gz") if p.stat().st_size > 0]
        else:
            command_rows.append(
                {
                    "sample_id": sample_id,
                    "step": "methylation_extractor_reused",
                    "command": f"existing_coverage={cov_candidates[0]}",
                }
            )
        if not cov_candidates:
            raise RuntimeError(f"No bismark.cov.gz output found for {sample_id}")
        cov_paths[sample_id] = cov_candidates[0]
        alignments.append({"sample_id": sample_id, "bam": str(bam), "coverage": str(cov_candidates[0])})

    series_list: list[pd.Series] = []
    parse_stats: list[dict] = []
    for sample_id, cov in cov_paths.items():
        series, stats = parse_cov_to_regions(cov, sample_id, args.min_coverage)
        sample_meta = samples[samples["sample_id"].astype(str).eq(sample_id)]
        if not sample_meta.empty:
            stats["age_weeks"] = sample_meta.iloc[0].get("age_weeks", "")
            stats["tissue"] = sample_meta.iloc[0].get("tissue", "")
            stats["run_accession"] = sample_meta.iloc[0].get("run_accession", "")
        series_list.append(series)
        parse_stats.append(stats)
    matrix = pd.concat(series_list, axis=1).sort_index() if series_list else pd.DataFrame()
    ref_regions = reference_regions()
    common = sorted(set(map(str, matrix.index)).intersection(ref_regions))
    region_matrix_path = out_dir / "GSE83947_bismark_pilot_region_matrix_5kb.parquet"
    region_stats_path = out_dir / "GSE83947_bismark_pilot_region_stats_5kb.csv"
    parse_stats_path = out_dir / "GSE83947_bismark_pilot_sample_parse_stats.csv"
    common_path = out_dir / "common_regions_with_v8_2.csv"
    command_manifest_path = out_dir / "bismark_command_manifest.csv"
    matrix.to_parquet(region_matrix_path)
    pd.DataFrame(parse_stats).to_csv(parse_stats_path, index=False)
    pd.DataFrame({"region_id": common}).to_csv(common_path, index=False)
    if not matrix.empty:
        stats = pd.DataFrame(
            {
                "region_id": matrix.index,
                "n_samples_present": matrix.notna().sum(axis=1).to_numpy(),
                "mean_beta": matrix.mean(axis=1, skipna=True).to_numpy(),
                "std_beta": matrix.std(axis=1, skipna=True).fillna(0.0).to_numpy(),
            }
        )
        chrom_start_end = stats["region_id"].str.extract(r"^(chr[^:]+):(\d+)-(\d+)$")
        stats["chrom"] = chrom_start_end[0]
        stats["start"] = pd.to_numeric(chrom_start_end[1], errors="coerce")
        stats["end"] = pd.to_numeric(chrom_start_end[2], errors="coerce")
        stats[["region_id", "chrom", "start", "end", "n_samples_present", "mean_beta", "std_beta"]].to_csv(region_stats_path, index=False)
    pd.DataFrame(command_rows).to_csv(command_manifest_path, index=False)
    beta_min = float(matrix.min(skipna=True).min()) if not matrix.empty else math.nan
    beta_max = float(matrix.max(skipna=True).max()) if not matrix.empty else math.nan
    state = {
        "timestamp": utc_now(),
        "status": "completed",
        "dataset": "GSE83947",
        "n_samples": int(matrix.shape[1]),
        "n_regions": int(matrix.shape[0]),
        "reference_regions": int(len(ref_regions)),
        "common_regions_with_v8_2": int(len(common)),
        "common_region_gate_50000": bool(len(common) >= COMMON_REGION_GATE),
        "beta_min": None if math.isnan(beta_min) else beta_min,
        "beta_max": None if math.isnan(beta_max) else beta_max,
        "beta_range_valid": bool(not matrix.empty and 0.0 <= beta_min <= beta_max <= 1.0),
        "alignment_outputs": alignments,
        "region_matrix_path": str(region_matrix_path.relative_to(ROOT)),
        "region_stats_path": str(region_stats_path.relative_to(ROOT)),
        "sample_parse_stats_path": str(parse_stats_path.relative_to(ROOT)),
        "common_regions_path": str(common_path.relative_to(ROOT)),
        "command_manifest_path": str(command_manifest_path.relative_to(ROOT)),
        "bismark_authorized": True,
        "training_authorized": False,
        "autoresearch_authorized": False,
        "headline_allowed": bool(len(common) >= COMMON_REGION_GATE),
        "decision": "matrix_gate_passed_pending_explicit_training_approval" if len(common) >= COMMON_REGION_GATE else "matrix_gate_failed_auxiliary_only",
    }
    write_json(out_dir / "gse83947_bismark_pilot_state.json", state)
    print(json.dumps(state, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

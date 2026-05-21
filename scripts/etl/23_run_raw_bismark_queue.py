#!/usr/bin/env python3
"""Run a guarded v21 raw FASTQ to Bismark COV to 5kb matrix queue.

This script is intentionally blocked by default. Pass `--authorize-bismark`
only after the v21 RALPH decision state approves a small raw pilot.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import subprocess
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd


ROOT = Path("/home/zdq-as/mouse_methyl_work")
DEFAULT_QUEUE = ROOT / "results" / "ralph_v21_raid_raw_clock" / "raw_etl_queue.csv"
DEFAULT_OUT_ROOT = Path("/data/mouse_methyl/processed_v21_raw_etl")
ENV_DIR = ROOT / "tools" / "route_b_bismark_env"
REF_DIR = ROOT / "references" / "GRCm38_ensembl102"
REGION_BIN_SIZE = 5_000


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str), encoding="utf-8")


def append_jsonl(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(payload, sort_keys=True, default=str) + "\n")


def env(env_dir: Path) -> dict[str, str]:
    out = os.environ.copy()
    out["PATH"] = f"{env_dir / 'bin'}:{out.get('PATH', '')}"
    return out


def run_command(cmd: list[str], cwd: Path, log_prefix: Path, env_dir: Path) -> None:
    log_prefix.parent.mkdir(parents=True, exist_ok=True)
    completed = subprocess.run(cmd, cwd=str(cwd), env=env(env_dir), text=True, capture_output=True, check=False)
    log_prefix.with_suffix(".stdout.log").write_text(completed.stdout, encoding="utf-8")
    log_prefix.with_suffix(".stderr.log").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(f"Command failed ({completed.returncode}): {' '.join(cmd)}\n{completed.stderr[-2000:]}")


def normalize_chrom(value: str) -> str | None:
    value = str(value).strip()
    clean = value.removeprefix("chr") if value.startswith("chr") else value
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


def parse_cov_paths_to_regions(paths: list[Path], sample_id: str, min_coverage: int) -> tuple[pd.Series, dict[str, Any]]:
    methylated_by_region: dict[str, float] = defaultdict(float)
    total_by_region: dict[str, float] = defaultdict(float)
    coverage_values: list[float] = []
    stats: dict[str, Any] = {
        "sample_id": sample_id,
        "coverage_paths": ";".join(str(path) for path in paths),
        "rows_total": 0,
        "rows_parseable_primary_autosomes": 0,
        "rows_pass_coverage": 0,
        "beta_min": None,
        "beta_max": None,
        "schema": "bismark_cov_6col",
    }
    for path in paths:
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
                    methylated = beta * total
                else:
                    continue
            stats["rows_parseable_primary_autosomes"] += 1
            stats["beta_min"] = beta if stats["beta_min"] is None else min(float(stats["beta_min"]), beta)
            stats["beta_max"] = beta if stats["beta_max"] is None else max(float(stats["beta_max"]), beta)
            if total < min_coverage:
                continue
            rid = region_id(chrom, pos)
            methylated_by_region[rid] += methylated
            total_by_region[rid] += total
            coverage_values.append(total)
            stats["rows_pass_coverage"] += 1
    values = {
        rid: methylated_by_region[rid] / total_by_region[rid]
        for rid in methylated_by_region
        if total_by_region[rid] > 0
    }
    series = pd.Series(values, dtype=np.float32, name=sample_id)
    stats["n_regions"] = int(len(series))
    stats["median_coverage"] = None if not coverage_values else round(float(np.median(coverage_values)), 3)
    stats["beta_min"] = None if stats["beta_min"] is None else round(float(stats["beta_min"]), 8)
    stats["beta_max"] = None if stats["beta_max"] is None else round(float(stats["beta_max"]), 8)
    return series, stats


def sample_region_gate_error(series: pd.Series, stats: dict[str, Any], min_sample_regions: int) -> str | None:
    n_regions = int(len(series))
    stats["sample_region_gate_min_regions"] = int(min_sample_regions)
    stats["sample_region_gate_status"] = "pass" if n_regions >= min_sample_regions else "fail"
    if n_regions < min_sample_regions:
        return f"Sample region QC failed: n_regions={n_regions} < min_sample_regions={min_sample_regions}"
    return None


def split_paths(value: Any) -> list[Path]:
    return [Path(item) for item in str(value or "").split(";") if item]


def existing_covs(run_dir: Path) -> list[Path]:
    return sorted([path for path in run_dir.glob("*.bismark.cov.gz") if path.stat().st_size > 0])


def bam_passes_quickcheck(path: Path, env_dir: Path) -> bool:
    samtools = env_dir / "bin" / "samtools"
    if not samtools.exists():
        return True
    completed = subprocess.run(
        [str(samtools), "quickcheck", "-q", str(path)],
        text=True,
        capture_output=True,
        check=False,
    )
    return completed.returncode == 0


def existing_bams(run_dir: Path, env_dir: Path) -> list[Path]:
    candidates: list[Path] = []
    reports = list(run_dir.glob("*_PE_report.txt")) + list(run_dir.glob("*_SE_report.txt"))
    has_report = any(path.stat().st_size > 0 for path in reports)
    if not has_report:
        return candidates
    for path in sorted(run_dir.glob("*.bam")):
        if path.stat().st_size <= 0 or ".temp." in path.name:
            continue
        if bam_passes_quickcheck(path, env_dir):
            candidates.append(path)
    return candidates


def process_sample(row: dict[str, Any], args: argparse.Namespace) -> tuple[str, str, pd.Series | None, dict[str, Any], list[dict[str, Any]]]:
    dataset = str(row["dataset"])
    sample_id = str(row["sample_id"])
    sample_dir = Path(args.out_root) / dataset / sample_id
    sample_dir.mkdir(parents=True, exist_ok=True)
    sample_matrix_path = sample_dir / f"{sample_id}_raw_region_beta.parquet"
    sample_stats_path = sample_dir / f"{sample_id}_raw_parse_stats.json"
    commands: list[dict[str, Any]] = []

    if sample_matrix_path.exists() and sample_stats_path.exists():
        frame = pd.read_parquet(sample_matrix_path)
        series = frame.iloc[:, 0].rename(sample_id)
        stats = json.loads(sample_stats_path.read_text(encoding="utf-8"))
        stats["resume_status"] = "reused_sample_region_beta"
        if error := sample_region_gate_error(series, stats, args.min_sample_regions):
            raise RuntimeError(error)
        return dataset, sample_id, series, stats, commands

    r1_paths = split_paths(row["r1_paths"])
    r2_paths = split_paths(row["r2_paths"])
    run_accessions = [item for item in str(row.get("run_accessions") or "").split(";") if item]
    if len(r1_paths) != len(r2_paths):
        raise RuntimeError(f"{dataset}/{sample_id} has mismatched R1/R2 counts: {len(r1_paths)} != {len(r2_paths)}")
    if not run_accessions or len(run_accessions) != len(r1_paths):
        run_accessions = [f"run_{idx + 1:03d}" for idx in range(len(r1_paths))]

    cov_paths: list[Path] = []
    for run_accession, r1, r2 in zip(run_accessions, r1_paths, r2_paths, strict=True):
        run_dir = sample_dir / "runs" / run_accession
        run_dir.mkdir(parents=True, exist_ok=True)
        bam_candidates = existing_bams(run_dir, Path(args.env_dir))
        if not bam_candidates:
            cmd = [
                str(Path(args.env_dir) / "bin" / "bismark"),
                "--genome_folder",
                str(args.reference_dir),
                "--bowtie2",
                "--path_to_bowtie2",
                str(Path(args.env_dir) / "bin"),
                "--parallel",
                str(args.bismark_threads),
                "-o",
                str(run_dir),
                "-1",
                str(r1),
                "-2",
                str(r2),
            ]
            run_command(cmd, run_dir, run_dir / "01_bismark_alignment", Path(args.env_dir))
            commands.append({"dataset": dataset, "sample_id": sample_id, "run_accession": run_accession, "step": "bismark_alignment", "command": " ".join(cmd)})
            bam_candidates = existing_bams(run_dir, Path(args.env_dir))
        else:
            commands.append({"dataset": dataset, "sample_id": sample_id, "run_accession": run_accession, "step": "bismark_alignment_reused", "command": f"existing_bam={bam_candidates[0]}"})
        if not bam_candidates:
            raise RuntimeError(f"No BAM output found for {dataset}/{sample_id}/{run_accession}")

        cov_candidates = existing_covs(run_dir)
        if not cov_candidates:
            cmd = [
                str(Path(args.env_dir) / "bin" / "bismark_methylation_extractor"),
                "--paired-end",
                "--gzip",
                "--bedGraph",
                "--parallel",
                str(args.extract_threads),
                "--buffer_size",
                "4G",
                "--output",
                str(run_dir),
                str(bam_candidates[0]),
            ]
            run_command(cmd, run_dir, run_dir / "02_bismark_methylation_extractor", Path(args.env_dir))
            commands.append({"dataset": dataset, "sample_id": sample_id, "run_accession": run_accession, "step": "methylation_extractor", "command": " ".join(cmd)})
            cov_candidates = existing_covs(run_dir)
        else:
            commands.append({"dataset": dataset, "sample_id": sample_id, "run_accession": run_accession, "step": "methylation_extractor_reused", "command": f"existing_cov={cov_candidates[0]}"})
        if not cov_candidates:
            raise RuntimeError(f"No Bismark COV output found for {dataset}/{sample_id}/{run_accession}")
        cov_paths.extend(cov_candidates)

    series, stats = parse_cov_paths_to_regions(cov_paths, sample_id, args.min_coverage)
    pd.DataFrame({sample_id: series}).to_parquet(sample_matrix_path)
    gate_error = sample_region_gate_error(series, stats, args.min_sample_regions)
    write_json(sample_stats_path, stats)
    if gate_error:
        raise RuntimeError(gate_error)
    return dataset, sample_id, series, stats, commands


def finalize_dataset_matrices(out_root: Path, sample_results: list[tuple[str, str, pd.Series | None, dict[str, Any], list[dict[str, Any]]]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    by_dataset: dict[str, list[pd.Series]] = defaultdict(list)
    for dataset, _, series, _, _ in sample_results:
        if series is not None and not series.empty:
            by_dataset[dataset].append(series)
    for dataset, series_list in by_dataset.items():
        matrix = pd.concat(series_list, axis=1).sort_index()
        dataset_dir = out_root / dataset
        dataset_dir.mkdir(parents=True, exist_ok=True)
        matrix_path = dataset_dir / f"{dataset}_raw_region_matrix_5kb.parquet"
        stats_path = dataset_dir / f"{dataset}_raw_region_stats_5kb.csv"
        matrix.to_parquet(matrix_path)
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
        stats[["region_id", "chrom", "start", "end", "n_samples_present", "mean_beta", "std_beta"]].to_csv(stats_path, index=False)
        rows.append(
            {
                "dataset": dataset,
                "matrix_path": str(matrix_path),
                "region_stats_path": str(stats_path),
                "n_samples": int(matrix.shape[1]),
                "n_regions": int(matrix.shape[0]),
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--queue", type=Path, default=DEFAULT_QUEUE)
    parser.add_argument("--out-root", type=Path, default=DEFAULT_OUT_ROOT)
    parser.add_argument("--env-dir", type=Path, default=ENV_DIR)
    parser.add_argument("--reference-dir", type=Path, default=REF_DIR)
    parser.add_argument("--parallel-samples", type=int, default=3)
    parser.add_argument("--bismark-threads", type=int, default=8)
    parser.add_argument("--extract-threads", type=int, default=4)
    parser.add_argument("--min-coverage", type=int, default=5)
    parser.add_argument("--min-sample-regions", type=int, default=50_000)
    parser.add_argument("--max-samples", type=int, default=0, help="0 means all queued rows.")
    parser.add_argument("--state-tag", default="", help="Optional suffix for state/status/log files when running a supplemental queue.")
    parser.add_argument("--skip-dataset-finalize", action="store_true", help="Only write per-sample outputs; do not rewrite dataset-level matrices.")
    parser.add_argument("--authorize-bismark", action="store_true")
    args = parser.parse_args()

    args.out_root.mkdir(parents=True, exist_ok=True)
    tag = f"_{args.state_tag.strip().replace('/', '_')}" if args.state_tag else ""
    state_path = args.out_root / f"v21_raw_etl_state{tag}.json"
    status_csv = args.out_root / f"v21_raw_etl_sample_status{tag}.csv"
    command_manifest_path = args.out_root / f"v21_raw_etl_command_manifest{tag}.csv"
    log_path = args.out_root / f"v21_raw_etl_log{tag}.jsonl"

    if not args.authorize_bismark:
        state = {
            "timestamp": utc_now(),
            "status": "blocked",
            "reason": "bismark_not_authorized",
            "queue": str(args.queue),
            "out_root": str(args.out_root),
            "training_authorized": False,
            "autoresearch_authorized": False,
        }
        write_json(state_path, state)
        print(json.dumps(state, indent=2, sort_keys=True))
        raise SystemExit(2)

    queue = pd.read_csv(args.queue)
    if args.max_samples and args.max_samples > 0:
        queue = queue.head(args.max_samples).copy()
    queue_rows = queue.to_dict(orient="records")
    sample_status: list[dict[str, Any]] = []
    sample_results: list[tuple[str, str, pd.Series | None, dict[str, Any], list[dict[str, Any]]]] = []
    all_commands: list[dict[str, Any]] = []
    started = time.time()

    state = {
        "timestamp": utc_now(),
        "status": "running",
        "pid": os.getpid(),
        "queue": str(args.queue),
        "out_root": str(args.out_root),
        "parallel_samples": int(args.parallel_samples),
        "bismark_threads": int(args.bismark_threads),
        "extract_threads": int(args.extract_threads),
        "n_queue_rows": int(len(queue_rows)),
        "training_authorized": False,
        "autoresearch_authorized": False,
    }
    write_json(state_path, state)

    with ThreadPoolExecutor(max_workers=max(1, args.parallel_samples)) as executor:
        futures = {executor.submit(process_sample, row, args): row for row in queue_rows}
        for future in as_completed(futures):
            row = futures[future]
            dataset = str(row["dataset"])
            sample_id = str(row["sample_id"])
            try:
                result = future.result()
                sample_results.append(result)
                _, _, series, stats, commands = result
                all_commands.extend(commands)
                sample_status.append(
                    {
                        "timestamp": utc_now(),
                        "dataset": dataset,
                        "sample_id": sample_id,
                        "status": "completed",
                        "n_regions": int(0 if series is None else len(series)),
                        "rows_total": stats.get("rows_total"),
                        "rows_pass_coverage": stats.get("rows_pass_coverage"),
                        "beta_min": stats.get("beta_min"),
                        "beta_max": stats.get("beta_max"),
                        "error": "",
                    }
                )
                append_jsonl(log_path, {"timestamp": utc_now(), "event": "sample_completed", "dataset": dataset, "sample_id": sample_id})
            except Exception as exc:
                sample_status.append(
                    {
                        "timestamp": utc_now(),
                        "dataset": dataset,
                        "sample_id": sample_id,
                        "status": "failed",
                        "n_regions": 0,
                        "error": str(exc)[-2000:],
                    }
                )
                append_jsonl(log_path, {"timestamp": utc_now(), "event": "sample_failed", "dataset": dataset, "sample_id": sample_id, "error": str(exc)[-2000:]})
            pd.DataFrame(sample_status).to_csv(status_csv, index=False)

    dataset_rows = [] if args.skip_dataset_finalize else finalize_dataset_matrices(args.out_root, sample_results)
    pd.DataFrame(all_commands).to_csv(command_manifest_path, index=False)
    status_df = pd.DataFrame(sample_status)
    completed = int(status_df["status"].eq("completed").sum()) if not status_df.empty else 0
    failed = int(status_df["status"].eq("failed").sum()) if not status_df.empty else 0
    final_state = {
        "timestamp": utc_now(),
        "status": "completed_with_failures" if failed else "completed",
        "pid": os.getpid(),
        "queue": str(args.queue),
        "out_root": str(args.out_root),
        "n_queue_rows": int(len(queue_rows)),
        "n_samples_completed": completed,
        "n_samples_failed": failed,
        "dataset_matrices": dataset_rows,
        "sample_status_csv": str(status_csv),
        "command_manifest": str(command_manifest_path),
        "exec_time_sec": round(time.time() - started, 1),
        "training_authorized": False,
        "autoresearch_authorized": False,
    }
    write_json(state_path, final_state)
    print(json.dumps(final_state, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

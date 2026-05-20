#!/usr/bin/env python3
"""
Phase 5: Beta Matrix Builder
==============================
Loads all .bismark.cov.gz files, applies coverage filter (>=5x),
builds sample x CpG beta matrix, exports to Parquet.

Usage:
    python scripts/etl/05_build_beta_matrix.py [--coverage 5] [--min-samples 0.5]
"""
import argparse, gzip, sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import numpy as np
import pandas as pd

ROOT     = Path("/home/zdq-as/mouse_methyl_work")
COV_DIR  = ROOT / "cov_files"
OUT_DIR  = ROOT / "results" / "beta_matrices"
META_FILE = ROOT / "metadata" / "unified_sample_metadata.csv"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# CpGs on sex chromosomes and MT — excluded from training
EXCLUDE_CHROMS = {"chrX", "chrY", "chrM", "X", "Y", "MT", "M"}


def load_cov_gz(path: Path, min_cov: int) -> pd.Series:
    """
    Load Bismark coverage file.
    Format: chr  start  end  methylation%  count_M  count_U
    Returns pd.Series: index=chr_pos string, value=beta [0,1] or NaN
    """
    records = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue
            chrom, start, _, pct, m, u = parts[:6]
            if chrom in EXCLUDE_CHROMS:
                continue
            total = int(m) + int(u)
            if total < min_cov:
                continue
            key = f"{chrom}_{start}"
            records[key] = float(m) / total
    return pd.Series(records, dtype=np.float32)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--coverage", type=int, default=5,
                    help="Minimum CpG coverage (default 5)")
    ap.add_argument("--min-samples", type=float, default=0.5,
                    help="Min fraction of samples a CpG must have coverage in (default 0.5)")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--gse", default=None, help="Process only this GEO dataset")
    args = ap.parse_args()

    if not META_FILE.exists():
        print(f"[ERROR] Run 01_build_metadata.py first.")
        sys.exit(1)

    meta = pd.read_csv(META_FILE)
    meta = meta[meta["fastq_on_disk"] == True]
    if args.gse:
        meta = meta[meta["dataset_batch"] == args.gse]

    # Find all .cov.gz files
    cov_files = {}
    for _, row in meta.drop_duplicates("sample_id").iterrows():
        gse, srr = row.dataset_batch, row.sample_id
        for pat in [f"{srr}*.bismark.cov.gz", f"{srr}*.cov.gz",
                    f"{srr}*_CpG_report.txt.gz"]:
            matches = list((COV_DIR / gse / srr).glob(pat)) if (COV_DIR / gse / srr).exists() else []
            if matches:
                cov_files[srr] = matches[0]
                break

    if not cov_files:
        print(f"[ERROR] No .cov.gz files found in {COV_DIR}")
        print("        Run phases 3-4 (Bismark + methylation extractor) first.")
        sys.exit(1)

    print(f"[Beta] Found {len(cov_files)} coverage files | min_cov={args.coverage}x")

    # Load in parallel
    sample_series = {}
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(load_cov_gz, path, args.coverage): srr
                   for srr, path in cov_files.items()}
        for i, fut in enumerate(as_completed(futures)):
            srr = futures[fut]
            try:
                s = fut.result()
                sample_series[srr] = s
                if (i + 1) % 50 == 0:
                    print(f"  Loaded {i+1}/{len(cov_files)}: {srr} ({len(s)} CpGs)")
            except Exception as e:
                print(f"  [WARN] Failed {srr}: {e}")

    print(f"[Beta] Loaded {len(sample_series)} samples")

    # Build DataFrame: CpGs × samples
    beta_df = pd.DataFrame(sample_series).astype(np.float32)  # rows=CpGs

    # Filter: CpGs present in >= min_samples fraction of samples
    min_n = int(len(sample_series) * args.min_samples)
    n_before = len(beta_df)
    beta_df = beta_df.dropna(thresh=min_n)
    print(f"[Beta] CpG filter ({args.min_samples*100:.0f}% samples): "
          f"{n_before} → {len(beta_df)} CpGs")

    # Export full matrix
    full_path = OUT_DIR / "beta_matrix_full.parquet"
    beta_df.to_parquet(full_path, compression="zstd")
    print(f"[Beta] Full matrix → {full_path} ({full_path.stat().st_size/1e6:.0f} MB)")

    # Union-filtered (80% samples — model-ready)
    min_n_strict = int(len(sample_series) * 0.8)
    beta_filtered = beta_df.dropna(thresh=min_n_strict)
    filtered_path = OUT_DIR / "beta_matrix_union_filtered.parquet"
    beta_filtered.to_parquet(filtered_path, compression="zstd")
    print(f"[Beta] Union-filtered matrix → {filtered_path} "
          f"({len(beta_filtered)} CpGs × {beta_filtered.shape[1]} samples)")

    # CpG coverage stats
    stats = pd.DataFrame({
        "cpg_id": beta_df.index,
        "n_samples_present": beta_df.notna().sum(axis=1).values,
        "mean_beta": beta_df.mean(axis=1).values,
        "std_beta": beta_df.std(axis=1).values,
    })
    stats.to_csv(OUT_DIR / "cpg_coverage_stats.csv", index=False)
    print(f"[Beta] CpG stats → {OUT_DIR / 'cpg_coverage_stats.csv'}")
    print(f"\n[Beta] ✅ Complete")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Phase 0: Thompson 2018 Beta Matrix Builder
==========================================
Parses the extracted .CGfinal.txt.gz files from GSE120137_RAW.tar.
Format of CGfinal: chr context pos count_M count_total beta ...
Example: chr1 cG 3014928 2 4 0.500 0.754

Input:  /home/zdq-as/mouse_methyl_work/cov_files/GSE120137/
Output: /home/zdq-as/mouse_methyl_work/results/phase0/beta_matrix_thompson.parquet
"""
import gzip, sys, multiprocessing
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd
import numpy as np

ROOT = Path("/home/zdq-as/mouse_methyl_work")
IN_DIR = ROOT / "cov_files" / "GSE120137"
OUT_DIR = ROOT / "results" / "phase0"
OUT_DIR.mkdir(parents=True, exist_ok=True)

EXCLUDE_CHROMS = {"chrX", "chrY", "chrM", "X", "Y", "M"}

def process_cgfinal(path: Path) -> pd.Series:
    records = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 6: continue
            chrom, _, pos, _, count_total, beta = parts[:6]
            if chrom in EXCLUDE_CHROMS: continue
            if int(count_total) < 5: continue  # coverage filter
            key = f"{chrom}_{pos}"
            records[key] = float(beta)
    return pd.Series(records, dtype=np.float32)

def main():
    print("=" * 60)
    print("Phase 0: Thompson Beta Matrix Builder from CGfinal")
    print("=" * 60)
    
    files = list(IN_DIR.glob("*.CGfinal.txt.gz"))
    if not files:
        print(f"[ERROR] No files found in {IN_DIR}")
        sys.exit(1)
        
    print(f"Found {len(files)} coverage files. Processing...")
    
    sample_series = {}
    with ProcessPoolExecutor(max_workers=min(16, multiprocessing.cpu_count())) as pool:
        futures = {pool.submit(process_cgfinal, p): p.name.split('.')[0] for p in files}
        for i, fut in enumerate(as_completed(futures)):
            srr = futures[fut]
            try:
                s = fut.result()
                sample_series[srr] = s
                if (i+1) % 50 == 0:
                    print(f"  Processed {i+1}/{len(files)}: {srr} ({len(s)} CpGs >=5x)")
            except Exception as e:
                print(f"  [ERROR] {srr}: {e}")
                
    print(f"Building matrix...")
    df = pd.DataFrame(sample_series).astype(np.float32)
    print(f"Raw shape: {df.shape}")
    
    # Filter missing
    min_samples = int(df.shape[1] * 0.8)
    df_filt = df.dropna(thresh=min_samples)
    print(f"Union-filtered shape (80% samples): {df_filt.shape}")
    
    out_path = OUT_DIR / "beta_matrix_thompson.parquet"
    df_filt.to_parquet(out_path, compression="zstd")
    print(f"Saved → {out_path} ({out_path.stat().st_size/1e6:.1f} MB)")
    print("✅ Complete.")

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Phase 2: Trim Galore Batch Processor (RRBS mode)
=================================================
Reads metadata/unified_sample_metadata.csv, finds FASTQs on disk,
and runs Trim Galore with --rrbs flag in parallel.

Usage:
    python scripts/etl/02_trim_galore_batch.py [--workers N] [--gse GSE120137]
    --workers: parallel jobs (default: 4, adjust to CPU count)
    --gse: process only this dataset (default: all)
    --dry-run: print commands without executing

Output: /home/zdq-as/mouse_methyl_work/trimmed/{GSE}/{SRR}/
"""
import argparse, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import pandas as pd

ROOT       = Path("/home/zdq-as/mouse_methyl_work")
RAW_DATA   = Path("/home/zdq-as/as-ds-ops/data/mouse_methyl/raw")
TRIM_OUT   = ROOT / "trimmed"
META_FILE  = ROOT / "metadata" / "unified_sample_metadata.csv"
LOG_DIR    = ROOT / "results" / "qc"
LOG_DIR.mkdir(parents=True, exist_ok=True)
TRIM_OUT.mkdir(parents=True, exist_ok=True)

TRIM_GALORE = "trim_galore"  # must be in PATH


def find_fastqs(gse: str, srr: str) -> tuple[Path | None, Path | None]:
    """Return (R1, R2) FASTQ paths or (single, None) for SE."""
    fastq_dir = RAW_DATA / gse / "fastq"
    r1 = fastq_dir / f"{srr}_1.fastq"
    r2 = fastq_dir / f"{srr}_2.fastq"
    se = fastq_dir / f"{srr}.fastq"
    # gzipped variants
    r1gz = fastq_dir / f"{srr}_1.fastq.gz"
    r2gz = fastq_dir / f"{srr}_2.fastq.gz"
    segz = fastq_dir / f"{srr}.fastq.gz"
    if r1.exists() and r2.exists(): return r1, r2
    if r1gz.exists() and r2gz.exists(): return r1gz, r2gz
    if se.exists(): return se, None
    if segz.exists(): return segz, None
    return None, None


def is_trimmed(gse: str, srr: str, paired: bool) -> bool:
    out_dir = TRIM_OUT / gse / srr
    if not out_dir.exists(): return False
    if paired:
        return (out_dir / f"{srr}_1_val_1.fq").exists() or \
               (out_dir / f"{srr}_1_val_1.fq.gz").exists()
    return (out_dir / f"{srr}_trimmed.fq").exists() or \
           (out_dir / f"{srr}_trimmed.fq.gz").exists()


def trim_sample(gse: str, srr: str, dry_run: bool) -> dict:
    r1, r2 = find_fastqs(gse, srr)
    if r1 is None:
        return {"srr": srr, "gse": gse, "status": "skip_no_fastq"}

    paired = r2 is not None
    if is_trimmed(gse, srr, paired):
        return {"srr": srr, "gse": gse, "status": "skip_done"}

    out_dir = TRIM_OUT / gse / srr
    out_dir.mkdir(parents=True, exist_ok=True)

    cmd = [TRIM_GALORE, "--rrbs", "--quality", "20", "--length", "20",
           "--cores", "2", "--fastqc",
           "--output_dir", str(out_dir)]
    if paired:
        cmd += ["--paired", str(r1), str(r2)]
    else:
        cmd += [str(r1)]

    if dry_run:
        print("  [DRY-RUN]", " ".join(cmd))
        return {"srr": srr, "gse": gse, "status": "dry_run"}

    t0 = time.time()
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        elapsed = time.time() - t0
        if result.returncode != 0:
            return {"srr": srr, "gse": gse, "status": "failed",
                    "error": result.stderr[-500:], "seconds": elapsed}
        return {"srr": srr, "gse": gse, "status": "done",
                "paired": paired, "seconds": round(elapsed, 1)}
    except subprocess.TimeoutExpired:
        return {"srr": srr, "gse": gse, "status": "timeout", "seconds": 3600}
    except Exception as e:
        return {"srr": srr, "gse": gse, "status": "error", "error": str(e)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--gse", default=None)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    if not META_FILE.exists():
        print(f"[ERROR] Run 01_build_metadata.py first: {META_FILE}")
        sys.exit(1)

    meta = pd.read_csv(META_FILE)
    meta = meta[meta["fastq_on_disk"] == True].copy()
    if args.gse:
        meta = meta[meta["dataset_batch"] == args.gse]

    tasks = [(row.dataset_batch, row.sample_id)
             for _, row in meta.drop_duplicates("sample_id").iterrows()]
    print(f"[Trim] {len(tasks)} samples to process | {args.workers} workers")

    results = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(trim_sample, gse, srr, args.dry_run): (gse, srr)
                   for gse, srr in tasks}
        for fut in as_completed(futures):
            r = fut.result()
            status = r.get("status", "?")
            print(f"  {r['gse']}/{r['srr']}: {status}")
            results.append(r)

    import json
    log_path = LOG_DIR / "trim_galore_log.jsonl"
    with open(log_path, "w") as f:
        for r in results:
            f.write(json.dumps(r) + "\n")

    done    = sum(1 for r in results if r["status"] == "done")
    skipped = sum(1 for r in results if r["status"].startswith("skip"))
    failed  = sum(1 for r in results if r["status"] in ("failed", "timeout", "error"))
    print(f"\n[Trim] Done: {done} | Skipped: {skipped} | Failed: {failed}")
    print(f"[Trim] Log → {log_path}")


if __name__ == "__main__":
    main()

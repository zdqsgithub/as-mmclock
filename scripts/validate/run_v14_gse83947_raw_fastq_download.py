#!/usr/bin/env python3
"""Download and verify the approved GSE83947 Route B raw FASTQ pilot files."""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")
V14_DIR = ROOT / "results" / "ralph_v14_public_data_rescue"
PREFLIGHT_DIR = V14_DIR / "gse83947_raw_pilot_preflight"
OUT_DIR = V14_DIR / "gse83947_raw_fastq_pilot"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def md5_file(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str), encoding="utf-8")


def https_url(raw: str) -> str:
    raw = raw.strip()
    if raw.startswith("http://") or raw.startswith("https://") or raw.startswith("ftp://"):
        return raw
    return "https://" + raw


def download(url: str, target: Path, log_dir: Path) -> subprocess.CompletedProcess[str]:
    target.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["wget", "-c", "--tries=5", "--timeout=60", "-O", str(target), url]
    completed = subprocess.run(cmd, text=True, capture_output=True, check=False)
    safe_name = target.name.replace("/", "_")
    (log_dir / f"{safe_name}.wget.stdout.log").write_text(completed.stdout, encoding="utf-8")
    (log_dir / f"{safe_name}.wget.stderr.log").write_text(completed.stderr, encoding="utf-8")
    return completed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", default=str(PREFLIGHT_DIR / "route_b_gse83947_raw_pilot_manifest.csv"))
    parser.add_argument("--out-dir", default=str(OUT_DIR))
    parser.add_argument("--authorize-download", action="store_true")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    fastq_dir = out_dir / "fastq"
    log_dir = out_dir / "download_logs"
    out_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    if not args.authorize_download:
        state = {
            "timestamp": utc_now(),
            "status": "blocked",
            "reason": "download_not_authorized",
            "raw_fastq_download_authorized": False,
            "bismark_authorized": False,
            "training_authorized": False,
            "autoresearch_authorized": False,
        }
        write_json(out_dir / "raw_fastq_download_state.json", state)
        print(json.dumps(state, indent=2))
        raise SystemExit(2)

    manifest = pd.read_csv(args.manifest)
    rows = []
    for _, row in manifest.iterrows():
        urls = str(row["fastq_ftp"]).split(";")
        md5s = str(row["fastq_md5"]).split(";")
        sizes = [int(x) for x in str(row["fastq_bytes"]).split(";") if str(x).isdigit()]
        for idx, url in enumerate(urls, start=1):
            url = https_url(url)
            target = fastq_dir / Path(url).name
            expected_md5 = md5s[idx - 1] if idx - 1 < len(md5s) else ""
            expected_size = sizes[idx - 1] if idx - 1 < len(sizes) else None
            status = "existing"
            returncode = 0
            if not target.exists() or (expected_size and target.stat().st_size != expected_size):
                completed = download(url, target, log_dir)
                returncode = completed.returncode
                status = "downloaded" if returncode == 0 else "failed"
            observed_size = target.stat().st_size if target.exists() else 0
            observed_md5 = md5_file(target) if target.exists() else ""
            rows.append(
                {
                    "dataset": row["dataset"],
                    "sample_id": row["sample_id"],
                    "run_accession": row["run_accession"],
                    "mate": idx,
                    "url": url,
                    "local_path": str(target),
                    "expected_size": expected_size,
                    "observed_size": observed_size,
                    "size_match": bool(expected_size is not None and observed_size == expected_size),
                    "expected_md5": expected_md5,
                    "observed_md5": observed_md5,
                    "md5_match": bool(expected_md5 and observed_md5 == expected_md5),
                    "sha256": sha256_file(target) if target.exists() else "",
                    "status": status,
                    "returncode": returncode,
                }
            )
    df = pd.DataFrame(rows)
    download_manifest = out_dir / "raw_fastq_download_manifest.csv"
    df.to_csv(download_manifest, index=False)
    all_pass = bool(not df.empty and df["size_match"].all() and df["md5_match"].all() and df["returncode"].eq(0).all())
    state = {
        "timestamp": utc_now(),
        "status": "completed" if all_pass else "failed",
        "n_files": int(len(df)),
        "n_samples": int(df["sample_id"].nunique()) if not df.empty else 0,
        "total_bytes": int(df["observed_size"].sum()) if not df.empty else 0,
        "all_size_match": bool(df["size_match"].all()) if not df.empty else False,
        "all_md5_match": bool(df["md5_match"].all()) if not df.empty else False,
        "download_manifest": str(download_manifest.relative_to(ROOT)),
        "raw_fastq_download_authorized": True,
        "bismark_authorized": False,
        "training_authorized": False,
        "autoresearch_authorized": False,
    }
    write_json(out_dir / "raw_fastq_download_state.json", state)
    print(json.dumps(state, indent=2, sort_keys=True))
    raise SystemExit(0 if all_pass else 1)


if __name__ == "__main__":
    main()

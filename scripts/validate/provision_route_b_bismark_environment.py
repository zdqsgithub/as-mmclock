#!/usr/bin/env python3
"""Provision the Route B local Bismark environment and GRCm38 index.

This script is intentionally project-local. It installs micromamba under
``tools/micromamba`` if needed, creates a conda-style environment under
``tools/route_b_bismark_env`` with Bismark/Bowtie2/Samtools, downloads the
Ensembl GRCm38 primary assembly FASTA, and builds a Bismark Bowtie2 bisulfite
index.

It does not download FASTQ, train models, or run autoresearch.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import tarfile
import tempfile
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path("/home/zdq-as/mouse_methyl_work")
OUT_DIR = ROOT / "results" / "ralph_v14_public_data_rescue" / "route_b_environment"
TOOLS_DIR = ROOT / "tools"
MAMBA_DIR = TOOLS_DIR / "micromamba"
MAMBA_BIN = MAMBA_DIR / "bin" / "micromamba"
MAMBA_ROOT = TOOLS_DIR / "micromamba-root"
ENV_DIR = TOOLS_DIR / "route_b_bismark_env"
REF_DIR = ROOT / "references" / "GRCm38_ensembl102"
REF_URL = "https://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz"
REF_GZ = REF_DIR / "Mus_musculus.GRCm38.dna.primary_assembly.fa.gz"
REF_FA = REF_DIR / "Mus_musculus.GRCm38.dna.primary_assembly.fa"
MANIFEST = OUT_DIR / "route_b_bismark_environment_manifest.json"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def run(cmd: list[str], *, env: dict[str, str] | None = None, cwd: Path | None = None, log_name: str | None = None) -> subprocess.CompletedProcess[str]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    merged_env = os.environ.copy()
    if env:
        merged_env.update(env)
    completed = subprocess.run(cmd, cwd=str(cwd or ROOT), env=merged_env, text=True, capture_output=True, check=False)
    if log_name:
        (OUT_DIR / f"{log_name}.stdout.log").write_text(completed.stdout, encoding="utf-8")
        (OUT_DIR / f"{log_name}.stderr.log").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(f"Command failed ({completed.returncode}): {' '.join(cmd)}\n{completed.stderr[-2000:]}")
    return completed


def sha256_file(path: Path) -> str:
    import hashlib

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_file(url: str, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    part = target.with_suffix(target.suffix + ".part")
    if target.exists() and target.stat().st_size > 0:
        return
    req = urllib.request.Request(url, headers={"User-Agent": "mouse-methyl-routeb-provision/1.0"})
    with urllib.request.urlopen(req, timeout=180) as response, part.open("wb") as handle:
        shutil.copyfileobj(response, handle)
    part.replace(target)


def ensure_micromamba() -> None:
    if MAMBA_BIN.exists():
        return
    MAMBA_DIR.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        tar_path = tmp_path / "micromamba.tar.bz2"
        download_file("https://micro.mamba.pm/api/micromamba/linux-64/latest", tar_path)
        with tarfile.open(tar_path, "r:bz2") as tar:
            member = next(m for m in tar.getmembers() if m.name.endswith("bin/micromamba"))
            member.name = "micromamba"
            tar.extract(member, path=tmp_path)
        (MAMBA_DIR / "bin").mkdir(parents=True, exist_ok=True)
        shutil.move(str(tmp_path / "micromamba"), MAMBA_BIN)
        MAMBA_BIN.chmod(0o755)


def mamba_env() -> dict[str, str]:
    return {
        "MAMBA_ROOT_PREFIX": str(MAMBA_ROOT),
        "PATH": f"{ENV_DIR / 'bin'}:{MAMBA_DIR / 'bin'}:{os.environ.get('PATH', '')}",
    }


def ensure_tool_env() -> None:
    bismark = ENV_DIR / "bin" / "bismark"
    bowtie2 = ENV_DIR / "bin" / "bowtie2"
    samtools = ENV_DIR / "bin" / "samtools"
    if bismark.exists() and bowtie2.exists() and samtools.exists():
        return
    ensure_micromamba()
    cmd = [
        str(MAMBA_BIN),
        "create",
        "-y",
        "-p",
        str(ENV_DIR),
        "-c",
        "conda-forge",
        "-c",
        "bioconda",
        "bismark=0.24.2",
        "bowtie2",
        "samtools",
    ]
    run(cmd, env=mamba_env(), log_name="01_micromamba_create_route_b_env")


def ensure_reference() -> None:
    REF_DIR.mkdir(parents=True, exist_ok=True)
    download_file(REF_URL, REF_GZ)
    if not REF_FA.exists():
        run(["gzip", "-dk", str(REF_GZ)], log_name="02_decompress_grcm38_fasta")


def bismark_index_ready() -> bool:
    bisulfite = REF_DIR / "Bisulfite_Genome"
    return bisulfite.exists() and bool(list(bisulfite.rglob("*.bt2")) or list(bisulfite.rglob("*.bt2l")))


def ensure_bismark_index(threads: int) -> None:
    if bismark_index_ready():
        return
    env = mamba_env()
    cmd = [
        str(ENV_DIR / "bin" / "bismark_genome_preparation"),
        "--bowtie2",
        "--parallel",
        str(threads),
        str(REF_DIR),
    ]
    run(cmd, env=env, cwd=ROOT, log_name="03_bismark_genome_preparation_grcm38")


def tool_version(cmd: list[str], env: dict[str, str]) -> str:
    completed = subprocess.run(cmd, env={**os.environ, **env}, text=True, capture_output=True, check=False)
    return (completed.stdout or completed.stderr).strip().splitlines()[0] if (completed.stdout or completed.stderr).strip() else ""


def build_manifest() -> dict[str, Any]:
    env = mamba_env()
    return {
        "timestamp": utc_now(),
        "status": "ready" if bismark_index_ready() else "blocked",
        "micromamba": str(MAMBA_BIN),
        "micromamba_version": tool_version([str(MAMBA_BIN), "--version"], env) if MAMBA_BIN.exists() else "",
        "environment_dir": str(ENV_DIR),
        "tools": {
            "bismark": str(ENV_DIR / "bin" / "bismark"),
            "bismark_version": tool_version([str(ENV_DIR / "bin" / "bismark"), "--version"], env) if (ENV_DIR / "bin" / "bismark").exists() else "",
            "bowtie2": str(ENV_DIR / "bin" / "bowtie2"),
            "bowtie2_version": tool_version([str(ENV_DIR / "bin" / "bowtie2"), "--version"], env) if (ENV_DIR / "bin" / "bowtie2").exists() else "",
            "samtools": str(ENV_DIR / "bin" / "samtools"),
            "samtools_version": tool_version([str(ENV_DIR / "bin" / "samtools"), "--version"], env) if (ENV_DIR / "bin" / "samtools").exists() else "",
        },
        "reference": {
            "assembly": "GRCm38/mm10",
            "source": "Ensembl release 102 primary assembly",
            "url": REF_URL,
            "fasta_gz": str(REF_GZ),
            "fasta_gz_sha256": sha256_file(REF_GZ) if REF_GZ.exists() else "",
            "fasta": str(REF_FA),
            "fasta_sha256": sha256_file(REF_FA) if REF_FA.exists() else "",
            "bismark_index_dir": str(REF_DIR / "Bisulfite_Genome"),
            "bismark_index_ready": bismark_index_ready(),
        },
        "authorizations": {
            "raw_fastq_download_authorized": False,
            "bismark_pilot_authorized": False,
            "training_authorized": False,
            "autoresearch_authorized": False,
        },
    }


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threads", type=int, default=12)
    parser.add_argument("--skip-index-build", action="store_true")
    args = parser.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ensure_micromamba()
    ensure_tool_env()
    ensure_reference()
    if not args.skip_index_build:
        ensure_bismark_index(args.threads)
    manifest = build_manifest()
    write_json(MANIFEST, manifest)
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

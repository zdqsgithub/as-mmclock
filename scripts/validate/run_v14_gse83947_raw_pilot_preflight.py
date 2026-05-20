#!/usr/bin/env python3
"""Preflight a Route B minimal raw FASTQ/Bismark pilot for GSE83947.

This script maps selected GSE83947 old-lung GSM samples to official
BioSample/SRA/ENA run metadata and checks whether local Bismark tooling and
reference indexes are available. It does not download FASTQ, run Bismark, train
models, or start autoresearch.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import shutil
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")
V14_DIR = ROOT / "results" / "ralph_v14_public_data_rescue"
OUT_DIR = V14_DIR / "gse83947_raw_pilot_preflight"
DOC_PATH = ROOT / "doc" / "20_analysis" / "44_20260519_v14_gse83947_route_b_raw_pilot_preflight.md"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str), encoding="utf-8")


def fetch_text(url: str, timeout: int = 60) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": "mouse-methyl-v14-routeb-preflight/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as response:
        data = response.read()
    if url.endswith(".gz"):
        data = gzip.decompress(data)
    return data.decode("utf-8", errors="replace")


def series_bucket(gse: str) -> str:
    digits = gse.removeprefix("GSE")
    return f"GSE{digits[:-3]}nnn"


def soft_url(gse: str) -> str:
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{series_bucket(gse)}/{gse}/soft/{gse}_family.soft.gz"


def parse_soft_blocks(gse: str) -> dict[str, dict[str, Any]]:
    cache_dir = ROOT / "metadata" / "geo_candidate_soft"
    cache = cache_dir / f"{gse}_family.soft.gz"
    if cache.exists():
        text = gzip.decompress(cache.read_bytes()).decode("utf-8", errors="replace")
    else:
        text = fetch_text(soft_url(gse))
        cache_dir.mkdir(parents=True, exist_ok=True)
        cache.write_bytes(gzip.compress(text.encode("utf-8")))
    blocks: dict[str, dict[str, Any]] = {}
    for raw_block in re.split(r"\n(?=\^SAMPLE = )", text):
        if not raw_block.startswith("^SAMPLE = "):
            continue
        sample_id = raw_block.splitlines()[0].split("=", 1)[1].strip()
        payload: dict[str, Any] = {
            "sample_id": sample_id,
            "title": "",
            "source_name": "",
            "characteristics": [],
            "biosample": "",
            "srx": "",
        }
        for line in raw_block.splitlines():
            if line.startswith("!Sample_title = "):
                payload["title"] = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_source_name_ch1 = "):
                payload["source_name"] = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_characteristics_ch1 = "):
                payload["characteristics"].append(line.split("=", 1)[1].strip())
            elif line.startswith("!Sample_relation = BioSample:"):
                match = re.search(r"(SAMN\d+)", line)
                if match:
                    payload["biosample"] = match.group(1)
            elif line.startswith("!Sample_relation = SRA:"):
                match = re.search(r"(SRX\d+)", line)
                if match:
                    payload["srx"] = match.group(1)
        blocks[sample_id] = payload
    return blocks


def sample_age_weeks(meta: dict[str, Any]) -> float | None:
    text = " | ".join(meta.get("characteristics", []))
    match = re.search(r"age\s*\(in months\)\s*:\s*(\d+(?:\.\d+)?)M", text, re.IGNORECASE)
    if not match:
        return None
    return round(float(match.group(1)) * 30.42 / 7.0, 3)


def sample_tissue(meta: dict[str, Any]) -> str:
    text = f"{meta.get('title', '')} {meta.get('source_name', '')} {' | '.join(meta.get('characteristics', []))}".lower()
    if "lung" in text:
        return "lung"
    if "heart" in text:
        return "heart"
    if "cortex" in text:
        return "brain_cortex"
    return "unknown"


def sample_treatment(meta: dict[str, Any]) -> str:
    text = " | ".join(meta.get("characteristics", []))
    match = re.search(r"treatment \(bisulfite/oxybisulfite\):\s*([^|]+)", text, re.IGNORECASE)
    return match.group(1).strip().lower() if match else ""


def ena_file_report(accession: str) -> list[dict[str, str]]:
    fields = "run_accession,sample_accession,secondary_sample_accession,experiment_accession,fastq_ftp,fastq_md5,fastq_bytes,library_strategy,library_layout,read_count,base_count"
    url = "https://www.ebi.ac.uk/ena/portal/api/filereport?" + urllib.parse.urlencode(
        {"accession": accession, "result": "read_run", "fields": fields, "format": "tsv", "limit": "0"}
    )
    text = fetch_text(url)
    return list(csv.DictReader(text.splitlines(), delimiter="\t"))


def command_path(name: str) -> str:
    project_tool = ROOT / "tools" / "route_b_bismark_env" / "bin" / name
    if project_tool.exists():
        return str(project_tool)
    return shutil.which(name) or ""


def find_bismark_indexes(search_roots: list[Path]) -> list[str]:
    hits: list[str] = []
    for root in search_roots:
        if not root.exists():
            continue
        for path in root.rglob("Bisulfite_Genome"):
            if path.is_dir():
                bt2 = list(path.rglob("*.bt2")) + list(path.rglob("*.bt2l"))
                if bt2:
                    hits.append(str(path))
        for path in root.rglob("*.bt2"):
            hits.append(str(path.parent))
            if len(hits) > 20:
                break
    return sorted(set(hits))[:20]


def write_report(path: Path, state: dict[str, Any], manifest: pd.DataFrame) -> None:
    view_cols = [
        "sample_id",
        "biosample",
        "srx",
        "run_accession",
        "age_weeks",
        "tissue",
        "fastq_total_bytes",
        "library_strategy",
        "library_layout",
        "download_authorized",
    ]
    view = manifest[[col for col in view_cols if col in manifest.columns]].copy()
    if view.empty:
        table = "No rows."
    else:
        for col in view.columns:
            view[col] = view[col].map(lambda value: "" if pd.isna(value) else str(value))
        header = "| " + " | ".join(view.columns) + " |"
        sep = "| " + " | ".join(["---"] * len(view.columns)) + " |"
        rows = ["| " + " | ".join(row) + " |" for row in view.to_numpy(dtype=str)]
        table = "\n".join([header, sep, *rows])
    lines = [
        "# v14 GSE83947 Route B Raw Pilot Preflight",
        "",
        f"Date: {utc_now()}",
        "",
        "## Summary",
        "",
        "This preflight mapped selected old-lung `GSE83947` GSM samples to official GEO SOFT, BioSample, SRA experiment, ENA run and FASTQ metadata. It did not download FASTQ, run Bismark, train models, or start autoresearch.",
        "",
        f"- Status: `{state['status']}`",
        f"- Decision: `{state['decision']}`",
        f"- Reason: {state['reason']}",
        "",
        "## Selected Runs",
        "",
        table,
        "",
        "## Environment Gate",
        "",
        f"- bismark: `{state['environment']['bismark'] or 'missing'}`",
        f"- bowtie2: `{state['environment']['bowtie2'] or 'missing'}`",
        f"- samtools: `{state['environment']['samtools'] or 'missing'}`",
        f"- Bismark index candidates: `{'; '.join(state['environment']['bismark_index_candidates']) or 'missing'}`",
        f"- free disk GB: `{state['environment']['free_disk_gb']}`",
        "",
        "## Guardrails",
        "",
        "- raw FASTQ download authorized: `false`",
        "- Bismark authorized: `false`",
        "- training authorized: `false`",
        "- autoresearch authorized: `false`",
        "",
        "## Next Action",
        "",
        "Provision a project-local Bismark/Bowtie2/Samtools environment and a traceable mm10/GRCm38 Bismark index, then rerun this preflight. Do not download FASTQ before the environment gate passes.",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", default="GSE83947")
    parser.add_argument("--pilot-manifest", default=str(V14_DIR / "pilot_run_manifest.csv"))
    parser.add_argument("--out-dir", default=str(OUT_DIR))
    parser.add_argument("--doc-path", default=str(DOC_PATH))
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    soft = parse_soft_blocks(args.dataset)
    ena_rows = ena_file_report("PRJNA327577")
    ena_by_biosample = {row.get("sample_accession", ""): row for row in ena_rows}
    ena_by_srx = {row.get("experiment_accession", ""): row for row in ena_rows}
    rows = []
    candidates = []
    for sample_id, meta in soft.items():
        ena = ena_by_biosample.get(meta.get("biosample", "")) or ena_by_srx.get(meta.get("srx", "")) or {}
        fastq_bytes = [int(x) for x in str(ena.get("fastq_bytes", "")).split(";") if x.isdigit()]
        age = sample_age_weeks(meta)
        tissue = sample_tissue(meta)
        treatment = sample_treatment(meta)
        if tissue != "lung" or age is None or age < 104.0 or treatment != "bisulfite":
            continue
        candidates.append((sum(fastq_bytes), int(ena.get("read_count") or 0), sample_id, meta, ena, age, tissue, treatment))
    candidates = sorted(candidates, reverse=True)[:3]
    for _, _, sample_id, meta, ena, age, tissue, treatment in candidates:
        fastq_bytes = [int(x) for x in str(ena.get("fastq_bytes", "")).split(";") if x.isdigit()]
        rows.append(
            {
                "dataset": args.dataset,
                "sample_id": sample_id,
                "title": meta.get("title", ""),
                "biosample": meta.get("biosample", ""),
                "srx": meta.get("srx", ""),
                "run_accession": ena.get("run_accession", ""),
                "library_strategy": ena.get("library_strategy", ""),
                "library_layout": ena.get("library_layout", ""),
                "read_count": ena.get("read_count", ""),
                "base_count": ena.get("base_count", ""),
                "fastq_ftp": ena.get("fastq_ftp", ""),
                "fastq_md5": ena.get("fastq_md5", ""),
                "fastq_bytes": ena.get("fastq_bytes", ""),
                "fastq_total_bytes": sum(fastq_bytes),
                "tissue": tissue,
                "age_weeks": age,
                "treatment": treatment,
                "selection_rule": "top_3_old_lung_bisulfite_by_fastq_bytes",
                "processed_supplement_name": "",
                "download_authorized": False,
                "bismark_authorized": False,
                "training_authorized": False,
                "autoresearch_authorized": False,
            }
        )
    manifest = pd.DataFrame(rows)
    manifest_path = out_dir / "route_b_gse83947_raw_pilot_manifest.csv"
    manifest.to_csv(manifest_path, index=False)

    usage = shutil.disk_usage(ROOT)
    env = {
        "bismark": command_path("bismark"),
        "bowtie2": command_path("bowtie2"),
        "samtools": command_path("samtools"),
        "bismark_index_candidates": find_bismark_indexes([ROOT / "references", ROOT / "genomes", ROOT / "data", ROOT]),
        "free_disk_gb": round(usage.free / (1024**3), 2),
    }
    metadata_gate = bool(
        len(manifest) == 3
        and manifest["run_accession"].fillna("").astype(str).ne("").all()
        and manifest["fastq_ftp"].fillna("").astype(str).ne("").all()
        and manifest["age_weeks"].fillna("").astype(str).ne("").all()
        and manifest["tissue"].astype(str).str.lower().eq("lung").all()
    )
    tool_gate = bool(env["bismark"] and env["bowtie2"] and env["samtools"] and env["bismark_index_candidates"])
    status = "ready_for_explicit_download_and_bismark_approval" if metadata_gate and tool_gate else "blocked"
    reason = []
    if not metadata_gate:
        reason.append("metadata_or_run_mapping_gate_failed")
    if not env["bismark"]:
        reason.append("missing_bismark")
    if not env["bowtie2"]:
        reason.append("missing_bowtie2")
    if not env["samtools"]:
        reason.append("missing_samtools")
    if not env["bismark_index_candidates"]:
        reason.append("missing_bismark_reference_index")
    state = {
        "timestamp": utc_now(),
        "dataset": args.dataset,
        "status": status,
        "decision": "do_not_download_fastq_until_environment_gate_passes" if status == "blocked" else "pending_explicit_download_and_compute_approval",
        "reason": ";".join(reason) or "metadata_and_environment_gates_passed",
        "n_selected_samples": int(len(manifest)),
        "total_fastq_bytes": int(manifest["fastq_total_bytes"].sum()) if not manifest.empty else 0,
        "metadata_gate": metadata_gate,
        "tool_gate": tool_gate,
        "environment": env,
        "manifest_path": str(manifest_path.relative_to(ROOT)),
        "raw_fastq_download_authorized": False,
        "bismark_authorized": False,
        "training_authorized": False,
        "autoresearch_authorized": False,
    }
    state_path = out_dir / "route_b_gse83947_raw_pilot_preflight_state.json"
    write_json(state_path, state)
    write_report(Path(args.doc_path), state, manifest)
    write_report(out_dir / "route_b_gse83947_raw_pilot_preflight_report.md", state, manifest)
    print(json.dumps(state, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

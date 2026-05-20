#!/usr/bin/env python3
"""v21 RAID raw FASTQ RALPH inventory, environment, and matrix gate.

Default mode is non-destructive: it inventories local RAID FASTQ files, builds a
small pilot ETL queue, validates local tooling, checks for any existing v21 raw
matrices, and writes a decision state. It does not run Bismark, download data,
train models, or start autoresearch.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


ROOT = Path("/home/zdq-as/mouse_methyl_work")
RAW_ROOT = Path("/data/mouse_methyl/raw")
DEFAULT_OUT = ROOT / "results" / "ralph_v21_raid_raw_clock"
DEFAULT_ETL_ROOT = Path("/data/mouse_methyl/processed_v21_raw_etl")
DEFAULT_METADATA = ROOT / "metadata" / "model_sample_metadata_v8.csv"
DEFAULT_REFERENCE = ROOT / "results" / "multidataset_v8_3_ablation" / "all6" / "all_rrbs_region_matrix_5kb.parquet"
DEFAULT_PROCESSED_ROOTS = [
    ROOT / "results" / "multidataset_v8_3_ablation" / "all6",
    ROOT / "results" / "multidataset_v8_2_prefilter_liftover",
    ROOT / "results" / "multidataset",
]
ENV_DIR = ROOT / "tools" / "route_b_bismark_env"
REF_DIR = ROOT / "references" / "GRCm38_ensembl102"
PRIORITY_DATASETS = ["GSE121141", "GSE80672", "GSE93957", "GSE60012", "GSE120137"]
FASTQ_SUFFIXES = (".fastq", ".fastq.gz", ".fq", ".fq.gz")
COMMON_REGION_GATE = 50_000
METADATA_OVERLAP_GATE = 0.95
AGE_COVERAGE_GATE = 0.95
RUN_RE = re.compile(r"((?:SRR|ERR|DRR)\d+)", re.IGNORECASE)
GSM_RE = re.compile(r"(GSM\d+)", re.IGNORECASE)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str), encoding="utf-8")


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def parse_priority(value: str | None) -> list[str]:
    if not value:
        return PRIORITY_DATASETS
    return [item.strip() for item in value.split(",") if item.strip()]


def infer_dataset(path: Path, raw_root: Path) -> str:
    try:
        return path.relative_to(raw_root).parts[0]
    except Exception:
        return "unknown"


def infer_run_accession(path: Path) -> str:
    match = RUN_RE.search(path.name)
    return match.group(1).upper() if match else path.stem.split("_")[0]


def infer_mate(path: Path) -> int | None:
    name = path.name
    patterns = [
        r"(?:^|[_\-.])R?([12])(?:[_\-.]|$)",
        r"_R([12])_",
        r"_([12])\.f(?:ast)?q(?:\.gz)?$",
    ]
    for pattern in patterns:
        match = re.search(pattern, name, flags=re.IGNORECASE)
        if match:
            return int(match.group(1))
    return None


def sample_id_from_text(*values: object) -> str:
    for value in values:
        match = GSM_RE.search(str(value or ""))
        if match:
            return match.group(1).upper()
    return ""


def inventory_fastq(raw_root: Path, datasets: list[str]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    dataset_dirs = [raw_root / ds for ds in datasets if (raw_root / ds).exists()]
    if not dataset_dirs:
        dataset_dirs = [path for path in raw_root.iterdir() if path.is_dir()] if raw_root.exists() else []
    for dataset_dir in dataset_dirs:
        dataset = dataset_dir.name
        for path in sorted(dataset_dir.rglob("*")):
            if not path.is_file():
                continue
            lower = path.name.lower()
            if not lower.endswith(FASTQ_SUFFIXES):
                continue
            stat = path.stat()
            rows.append(
                {
                    "dataset": dataset,
                    "run_accession": infer_run_accession(path),
                    "mate": infer_mate(path),
                    "file_name": path.name,
                    "local_path": str(path),
                    "size_bytes": int(stat.st_size),
                    "size_gib": round(float(stat.st_size) / 1024**3, 6),
                    "mtime_utc": datetime.fromtimestamp(stat.st_mtime, tz=timezone.utc).isoformat(),
                    "is_gzip": lower.endswith(".gz"),
                    "raw_root": str(raw_root),
                }
            )
    return pd.DataFrame(rows)


def load_sra_metadata(raw_root: Path, datasets: list[str]) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for dataset in datasets:
        path = raw_root / dataset / "metadata" / "sra_metadata.tsv"
        if not path.exists():
            continue
        try:
            df = pd.read_csv(path, sep=None, engine="python")
        except Exception:
            df = pd.read_csv(path, sep="\t")
        df["dataset"] = dataset
        df["sra_metadata_path"] = str(path)
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    meta = pd.concat(frames, ignore_index=True)
    if "run_accession" in meta.columns:
        meta["run_accession"] = meta["run_accession"].astype(str).str.upper()
    for col in ["experiment_title", "experiment_desc", "sample_title"]:
        if col not in meta.columns:
            meta[col] = ""
    meta["sample_id"] = [
        sample_id_from_text(row.get("experiment_title"), row.get("experiment_desc"), row.get("sample_title"))
        for _, row in meta.iterrows()
    ]
    return meta


def build_sample_run_manifest(inventory: pd.DataFrame, sra_meta: pd.DataFrame) -> pd.DataFrame:
    if inventory.empty:
        return pd.DataFrame()
    grouped = (
        inventory.groupby(["dataset", "run_accession"], dropna=False)
        .agg(
            local_fastq_files=("local_path", lambda s: ";".join(sorted(map(str, s)))),
            local_fastq_count=("local_path", "count"),
            local_size_bytes=("size_bytes", "sum"),
            mates_present=("mate", lambda s: ";".join(str(int(value)) for value in sorted(s.dropna().unique()))),
        )
        .reset_index()
    )
    if sra_meta.empty or "run_accession" not in sra_meta.columns:
        grouped["sample_id"] = ""
        grouped["library_layout"] = ""
        grouped["library_strategy"] = ""
        grouped["sample_accession"] = ""
        grouped["experiment_accession"] = ""
        grouped["biosample"] = ""
        grouped["metadata_status"] = "missing_sra_metadata"
        return grouped
    keep_cols = [
        col
        for col in [
            "dataset",
            "run_accession",
            "sample_id",
            "library_layout",
            "library_strategy",
            "library_selection",
            "sample_accession",
            "experiment_accession",
            "experiment_title",
            "biosample",
            "bioproject",
            "run_total_bases",
        ]
        if col in sra_meta.columns
    ]
    merged = grouped.merge(sra_meta[keep_cols].drop_duplicates(), on=["dataset", "run_accession"], how="left")
    merged["metadata_status"] = np.where(merged.get("sample_id", "").fillna("").astype(str).ne(""), "mapped", "missing_sample_id")
    return merged


def build_pairing_qc(inventory: pd.DataFrame, sample_run: pd.DataFrame) -> pd.DataFrame:
    if inventory.empty:
        return pd.DataFrame()
    rows: list[dict[str, Any]] = []
    sample_lookup = sample_run.set_index(["dataset", "run_accession"]).to_dict(orient="index") if not sample_run.empty else {}
    for (dataset, run), group in inventory.groupby(["dataset", "run_accession"], dropna=False):
        mates = set(int(v) for v in group["mate"].dropna().astype(int).tolist())
        missing = []
        if mates and 1 not in mates:
            missing.append("1")
        if mates and 2 not in mates:
            missing.append("2")
        layout = "paired" if mates == {1, 2} else ("single_or_unknown" if not mates else "incomplete_paired")
        meta = sample_lookup.get((dataset, run), {})
        rows.append(
            {
                "dataset": dataset,
                "run_accession": run,
                "sample_id": meta.get("sample_id", ""),
                "library_layout": meta.get("library_layout", ""),
                "fastq_count": int(len(group)),
                "mates_present": ";".join(str(m) for m in sorted(mates)),
                "missing_mate": ";".join(missing),
                "pairing_status": "complete_pair" if layout == "paired" else layout,
                "total_size_bytes": int(group["size_bytes"].sum()),
                "total_size_gib": round(float(group["size_bytes"].sum()) / 1024**3, 6),
                "local_paths": ";".join(sorted(group["local_path"].astype(str).tolist())),
            }
        )
    return pd.DataFrame(rows)


def build_raw_etl_queue(
    pairing_qc: pd.DataFrame,
    priority: list[str],
    pilot_samples: int,
    *,
    queue_all_priority_pilots: bool,
) -> pd.DataFrame:
    if pairing_qc.empty:
        return pd.DataFrame()
    complete = pairing_qc[pairing_qc["pairing_status"].eq("complete_pair")].copy()
    if complete.empty:
        return pd.DataFrame()
    rows: list[dict[str, Any]] = []
    priority_rank = {dataset: idx for idx, dataset in enumerate(priority)}
    complete["priority_rank"] = complete["dataset"].map(priority_rank).fillna(len(priority)).astype(int)
    for dataset in priority:
        ds = complete[complete["dataset"].eq(dataset)].copy()
        if ds.empty:
            continue
        ds["sample_group"] = ds["sample_id"].fillna("").astype(str)
        ds.loc[ds["sample_group"].eq(""), "sample_group"] = ds.loc[ds["sample_group"].eq(""), "run_accession"]
        for sample_group, group in ds.sort_values(["sample_group", "run_accession"]).groupby("sample_group"):
            r1_paths = []
            r2_paths = []
            run_accessions = []
            for _, row in group.iterrows():
                parts = [Path(p) for p in str(row["local_paths"]).split(";") if p]
                r1 = [str(p) for p in parts if infer_mate(p) == 1]
                r2 = [str(p) for p in parts if infer_mate(p) == 2]
                if not r1 or not r2:
                    continue
                r1_paths.extend(sorted(r1))
                r2_paths.extend(sorted(r2))
                run_accessions.append(str(row["run_accession"]))
            if not r1_paths or not r2_paths:
                continue
            rows.append(
                {
                    "dataset": dataset,
                    "sample_id": sample_group,
                    "run_accessions": ";".join(run_accessions),
                    "r1_paths": ";".join(r1_paths),
                    "r2_paths": ";".join(r2_paths),
                    "n_runs": len(run_accessions),
                    "n_fastq_files": len(r1_paths) + len(r2_paths),
                    "pilot_priority": priority_rank[dataset],
                    "queue_status": "planned_pilot",
                    "bismark_authorized": False,
                    "training_authorized": False,
                }
            )
            if len([r for r in rows if r["dataset"] == dataset]) >= pilot_samples:
                break
        if rows and not queue_all_priority_pilots:
            break
    return pd.DataFrame(rows)


def command_version(cmd: list[str]) -> dict[str, Any]:
    try:
        completed = subprocess.run(cmd, text=True, capture_output=True, timeout=30, check=False)
        text = (completed.stdout or completed.stderr or "").strip().splitlines()
        return {
            "command": " ".join(cmd),
            "exists": True,
            "returncode": completed.returncode,
            "version_head": text[:3],
        }
    except FileNotFoundError:
        return {"command": " ".join(cmd), "exists": False, "returncode": None, "version_head": []}
    except Exception as exc:
        return {"command": " ".join(cmd), "exists": True, "returncode": None, "error": str(exc), "version_head": []}


def validate_environment(reference_dir: Path, env_dir: Path, raw_root: Path, etl_root: Path) -> dict[str, Any]:
    uv = command_version(["uv", "--version"])
    core_python = ROOT / ".venv-core" / "bin" / "python"
    import_cmd = [str(core_python), "-c", "import pandas, numpy, sklearn, pyarrow; print('imports_ok')"]
    if not core_python.exists():
        import_cmd = ["uv", "run", "python", "-c", "import pandas, numpy, sklearn, pyarrow; print('imports_ok')"]
    imports = command_version(import_cmd)
    bismark_bin = env_dir / "bin" / "bismark"
    bowtie_bin = env_dir / "bin" / "bowtie2"
    samtools_bin = env_dir / "bin" / "samtools"
    ct = reference_dir / "Bisulfite_Genome" / "CT_conversion"
    ga = reference_dir / "Bisulfite_Genome" / "GA_conversion"
    ct_indexes = list(ct.glob("*.bt2*")) if ct.exists() else []
    ga_indexes = list(ga.glob("*.bt2*")) if ga.exists() else []
    disks = {}
    for name, path in {"raw_root": raw_root, "etl_root": etl_root, "project_root": ROOT}.items():
        target = path if path.exists() else path.parent
        usage = shutil.disk_usage(target)
        disks[name] = {
            "path": str(path),
            "total_gib": round(usage.total / 1024**3, 3),
            "used_gib": round(usage.used / 1024**3, 3),
            "free_gib": round(usage.free / 1024**3, 3),
        }
    checks = {
        "timestamp": utc_now(),
        "uv": uv,
        "project_imports": imports,
        "bismark": command_version([str(bismark_bin), "--version"]),
        "bowtie2": command_version([str(bowtie_bin), "--version"]),
        "samtools": command_version([str(samtools_bin), "--version"]),
        "reference_dir": str(reference_dir),
        "reference_fasta_exists": (reference_dir / "Mus_musculus.GRCm38.dna.primary_assembly.fa").exists(),
        "ct_index_count": len(ct_indexes),
        "ga_index_count": len(ga_indexes),
        "disk": disks,
    }
    checks["environment_gate_passed"] = bool(
        uv.get("exists")
        and imports.get("returncode") == 0
        and bismark_bin.exists()
        and bowtie_bin.exists()
        and samtools_bin.exists()
        and checks["reference_fasta_exists"]
        and len(ct_indexes) > 0
        and len(ga_indexes) > 0
        and disks["etl_root"]["free_gib"] >= 100
    )
    return checks


def read_reference_regions(path: Path) -> set[str]:
    if not path.exists():
        return set()
    try:
        ref = pd.read_parquet(path, columns=[])
    except Exception:
        ref = pd.read_parquet(path)
    return set(map(str, ref.index))


def find_processed_matrix(dataset: str, roots: list[Path]) -> Path | None:
    candidates = []
    for root in roots:
        if root.exists():
            candidates.extend(root.glob(f"{dataset}*_region_matrix_5kb.parquet"))
    return sorted(candidates)[0] if candidates else None


def matrix_dataset_from_path(path: Path) -> str:
    for part in path.parts:
        if part.startswith("GSE"):
            return part
    match = re.search(r"(GSE\d+)", path.name)
    return match.group(1) if match else "unknown"


def consistency_audit(raw_matrix: pd.DataFrame, processed_path: Path | None) -> dict[str, Any]:
    if processed_path is None or not processed_path.exists() or raw_matrix.empty:
        return {"processed_matrix_path": "", "status": "not_available"}
    try:
        processed = pd.read_parquet(processed_path)
    except Exception as exc:
        return {"processed_matrix_path": str(processed_path), "status": "read_failed", "error": str(exc)}
    common_regions = raw_matrix.index.intersection(processed.index)
    common_samples = raw_matrix.columns.intersection(processed.columns)
    if len(common_regions) == 0 or len(common_samples) == 0:
        return {
            "processed_matrix_path": str(processed_path),
            "status": "no_common_region_or_sample",
            "common_regions": int(len(common_regions)),
            "common_samples": int(len(common_samples)),
        }
    corrs = []
    for sample in common_samples[: min(10, len(common_samples))]:
        a = raw_matrix.loc[common_regions, sample].astype(float)
        b = processed.loc[common_regions, sample].astype(float)
        valid = a.notna() & b.notna()
        if valid.sum() >= 100:
            corrs.append(float(a[valid].corr(b[valid])))
    return {
        "processed_matrix_path": str(processed_path),
        "status": "completed" if corrs else "insufficient_overlap_for_correlation",
        "common_regions": int(len(common_regions)),
        "common_samples": int(len(common_samples)),
        "median_sample_corr": None if not corrs else round(float(np.nanmedian(corrs)), 6),
    }


def build_matrix_gate_table(
    etl_root: Path,
    reference_matrix: Path,
    metadata_path: Path,
    priority: list[str],
    processed_roots: list[Path],
) -> pd.DataFrame:
    reference_regions = read_reference_regions(reference_matrix)
    meta = pd.read_csv(metadata_path) if metadata_path.exists() else pd.DataFrame(columns=["sample_id", "age_days"])
    meta_ids = set(meta.get("sample_id", pd.Series(dtype=str)).astype(str))
    age_ids = set(meta[meta.get("age_days", pd.Series(dtype=float)).notna()]["sample_id"].astype(str)) if not meta.empty else set()
    paths = sorted(etl_root.rglob("*_raw_region_matrix_5kb.parquet")) if etl_root.exists() else []
    if not paths:
        return pd.DataFrame(
            [
                {
                    "dataset": dataset,
                    "matrix_path": "",
                    "status": "pending_raw_pilot",
                    "n_samples": 0,
                    "n_regions": 0,
                    "common_regions": 0,
                    "metadata_overlap_fraction": 0.0,
                    "age_coverage_fraction": 0.0,
                    "beta_range_valid": False,
                    "sex_mt_excluded": False,
                    "matrix_gate_passed": False,
                    "decision": "raw_inventory_ready_pilot_pending",
                }
                for dataset in priority
            ]
        )
    rows: list[dict[str, Any]] = []
    for path in paths:
        dataset = matrix_dataset_from_path(path)
        try:
            matrix = pd.read_parquet(path)
            matrix.index = matrix.index.astype(str)
            matrix.columns = matrix.columns.astype(str)
            n_samples = int(matrix.shape[1])
            n_regions = int(matrix.shape[0])
            common_regions = len(set(matrix.index).intersection(reference_regions)) if reference_regions else 0
            beta_min = float(matrix.min(skipna=True).min()) if not matrix.empty else np.nan
            beta_max = float(matrix.max(skipna=True).max()) if not matrix.empty else np.nan
            metadata_overlap = sum(col in meta_ids for col in matrix.columns) / n_samples if n_samples else 0.0
            age_coverage = sum(col in age_ids for col in matrix.columns) / n_samples if n_samples else 0.0
            sex_mt_excluded = not any(
                str(idx).startswith(("chrX:", "chrY:", "chrM:", "chrMT:", "MT:", "X:", "Y:")) for idx in matrix.index
            )
            beta_valid = bool(n_samples and n_regions and 0.0 <= beta_min <= beta_max <= 1.0)
            processed_path = find_processed_matrix(dataset, processed_roots)
            audit = consistency_audit(matrix, processed_path)
            gate = bool(
                metadata_overlap >= METADATA_OVERLAP_GATE
                and age_coverage >= AGE_COVERAGE_GATE
                and beta_valid
                and sex_mt_excluded
                and common_regions >= COMMON_REGION_GATE
            )
            decision = (
                "raw_matrix_ready_for_learn_pending_training_approval"
                if gate
                else "raw_matrix_gate_failed_or_pending_repair"
            )
            rows.append(
                {
                    "dataset": dataset,
                    "matrix_path": str(path),
                    "status": "completed",
                    "n_samples": n_samples,
                    "n_regions": n_regions,
                    "common_regions": int(common_regions),
                    "metadata_overlap_fraction": round(float(metadata_overlap), 6),
                    "age_coverage_fraction": round(float(age_coverage), 6),
                    "beta_min": None if np.isnan(beta_min) else round(beta_min, 8),
                    "beta_max": None if np.isnan(beta_max) else round(beta_max, 8),
                    "beta_range_valid": beta_valid,
                    "sex_mt_excluded": bool(sex_mt_excluded),
                    "processed_consistency_status": audit.get("status"),
                    "processed_consistency_median_corr": audit.get("median_sample_corr"),
                    "processed_common_regions": audit.get("common_regions"),
                    "processed_common_samples": audit.get("common_samples"),
                    "matrix_gate_passed": gate,
                    "decision": decision,
                }
            )
        except Exception as exc:
            rows.append(
                {
                    "dataset": dataset,
                    "matrix_path": str(path),
                    "status": "read_failed",
                    "error": str(exc),
                    "matrix_gate_passed": False,
                    "decision": "raw_matrix_gate_failed_or_pending_repair",
                }
            )
    return pd.DataFrame(rows)


def write_ml_command_manifest(out_dir: Path, matrix_gate: pd.DataFrame, metadata_path: Path) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    passed = matrix_gate[matrix_gate.get("matrix_gate_passed", False).astype(bool)] if not matrix_gate.empty else pd.DataFrame()
    for _, row in passed.iterrows():
        matrix = row["matrix_path"]
        for model in ["ridge", "elasticnet", "lgbm"]:
            for top in ["500", "1000", "2000"]:
                cmd = [
                    "uv",
                    "run",
                    "python",
                    "scripts/train/train_clock.py",
                    "--matrix_path",
                    str(matrix),
                    "--metadata_path",
                    str(metadata_path),
                    "--feature_type",
                    "region",
                    "--model_type",
                    model,
                    "--n_feature_prefilter",
                    top,
                    "--preprocess",
                    "standard",
                    "--output_dir",
                    f"results/ralph_v21_raid_raw_clock/ml/{row['dataset']}/{model}_top{top}",
                ]
                rows.append(
                    {
                        "dataset": row["dataset"],
                        "model_type": model,
                        "top_regions": top,
                        "command": " ".join(cmd),
                        "training_authorized": False,
                        "autoresearch_authorized": False,
                    }
                )
    manifest = pd.DataFrame(rows)
    manifest.to_csv(out_dir / "local_ml_command_manifest.csv", index=False)
    return manifest


def write_report(
    out_dir: Path,
    inventory: pd.DataFrame,
    sample_run: pd.DataFrame,
    pairing_qc: pd.DataFrame,
    queue: pd.DataFrame,
    env_state: dict[str, Any],
    matrix_gate: pd.DataFrame,
    decision: dict[str, Any],
) -> None:
    def markdown_table(df: pd.DataFrame) -> str:
        if df.empty:
            return ""
        headers = list(df.columns)
        lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
        for _, row in df.iterrows():
            values = [str(row.get(col, "")) for col in headers]
            lines.append("| " + " | ".join(value.replace("|", "\\|") for value in values) + " |")
        return "\n".join(lines)

    lines = [
        "# v21 RAID Raw RALPH Report",
        "",
        f"Date: {utc_now()}",
        "",
        "## Scope",
        "",
        "Read-only inventory and gate preflight for local RAID FASTQ. This run did not start Bismark, downloads, training, or autoresearch.",
        "",
        "## Raw Inventory",
        "",
        f"- FASTQ files: `{len(inventory)}`",
        f"- FASTQ size GiB: `{round(float(inventory['size_bytes'].sum()) / 1024**3, 3) if not inventory.empty else 0}`",
        f"- datasets with FASTQ: `{', '.join(sorted(inventory['dataset'].unique())) if not inventory.empty else 'none'}`",
        f"- sample/run rows: `{len(sample_run)}`",
        f"- complete paired runs: `{int(pairing_qc['pairing_status'].eq('complete_pair').sum()) if not pairing_qc.empty else 0}`",
        f"- pilot queue rows: `{len(queue)}`",
        "",
        "## Environment Gate",
        "",
        f"- passed: `{env_state.get('environment_gate_passed')}`",
        f"- CT index files: `{env_state.get('ct_index_count')}`",
        f"- GA index files: `{env_state.get('ga_index_count')}`",
        "",
        "## Matrix Gate",
        "",
    ]
    if matrix_gate.empty:
        lines.append("No matrix gate rows were produced.")
    else:
        cols = [
            "dataset",
            "status",
            "n_samples",
            "n_regions",
            "common_regions",
            "metadata_overlap_fraction",
            "age_coverage_fraction",
            "beta_range_valid",
            "sex_mt_excluded",
            "matrix_gate_passed",
            "decision",
        ]
        lines.append(markdown_table(matrix_gate[[c for c in cols if c in matrix_gate.columns]]))
    lines.extend(
        [
            "",
            "## Decision",
            "",
            f"- state: `{decision['state']}`",
            f"- training authorized: `{decision['training_authorized']}`",
            f"- autoresearch authorized: `{decision['autoresearch_authorized']}`",
            f"- next action: {decision['next_action']}",
            "",
            "## Outputs",
            "",
            "- `raw_inventory.csv`",
            "- `sample_run_manifest.csv`",
            "- `pairing_qc_manifest.csv`",
            "- `raw_etl_queue.csv`",
            "- `environment_gate.json`",
            "- `matrix_gate_table.csv`",
            "- `watchdog_state.json`",
            "- `ralph_decision_state.json`",
        ]
    )
    write_text(out_dir / "v21_raid_raw_clock_report.md", "\n".join(lines) + "\n")


def decide(env_state: dict[str, Any], matrix_gate: pd.DataFrame, queue: pd.DataFrame) -> dict[str, Any]:
    gate_passed = bool((not matrix_gate.empty) and matrix_gate.get("matrix_gate_passed", pd.Series(dtype=bool)).fillna(False).astype(bool).any())
    failed_critical = 0
    if not matrix_gate.empty:
        pending = matrix_gate.get("status", pd.Series(dtype=str)).astype(str).eq("pending_raw_pilot")
        failed = ~matrix_gate.get("matrix_gate_passed", pd.Series(False, index=matrix_gate.index)).fillna(False).astype(bool)
        failed_critical = int((failed & ~pending).sum())
    if failed_critical >= 3:
        state = "raw_re_evaluation_required_no_training"
        next_action = "Stop raw ETL expansion and re-evaluate Bismark/GEO/SRA/ENA/metadata before any training."
    elif gate_passed:
        state = "raw_matrix_ready_for_learn_pending_training_approval"
        next_action = "Review matrix gate and record explicit training approval before local ML."
    elif not env_state.get("environment_gate_passed"):
        state = "blocked_environment_provisioning_required"
        next_action = "Provision uv/Bismark/Bowtie2/Samtools/reference index before raw pilot."
    elif queue.empty:
        state = "blocked_no_complete_paired_pilot_fastq"
        next_action = "Resolve raw FASTQ pairing or SRA metadata before Bismark."
    else:
        state = "raw_inventory_ready_pilot_pending_bismark_approval"
        next_action = "Run only the guarded 2-3 sample raw pilot after explicit Bismark approval."
    return {
        "timestamp": utc_now(),
        "state": state,
        "training_authorized": False,
        "autoresearch_authorized": False,
        "download_authorized": False,
        "bismark_authorized": False,
        "matrix_gate_passed": gate_passed,
        "critical_failed_matrix_attempts": failed_critical,
        "next_action": next_action,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-root", type=Path, default=RAW_ROOT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--etl-root", type=Path, default=DEFAULT_ETL_ROOT)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--reference-matrix", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--reference-dir", type=Path, default=REF_DIR)
    parser.add_argument("--env-dir", type=Path, default=ENV_DIR)
    parser.add_argument("--priority-datasets", default=",".join(PRIORITY_DATASETS))
    parser.add_argument("--pilot-samples", type=int, default=3)
    parser.add_argument(
        "--queue-all-priority-pilots",
        action="store_true",
        help="Queue pilot rows for every priority dataset. Default queues only the first available priority dataset.",
    )
    args = parser.parse_args()

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    priority = parse_priority(args.priority_datasets)

    inventory = inventory_fastq(args.raw_root, priority)
    sra_meta = load_sra_metadata(args.raw_root, priority)
    sample_run = build_sample_run_manifest(inventory, sra_meta)
    pairing_qc = build_pairing_qc(inventory, sample_run)
    queue = build_raw_etl_queue(
        pairing_qc,
        priority,
        args.pilot_samples,
        queue_all_priority_pilots=args.queue_all_priority_pilots,
    )
    env_state = validate_environment(args.reference_dir, args.env_dir, args.raw_root, args.etl_root)
    matrix_gate = build_matrix_gate_table(
        args.etl_root,
        args.reference_matrix,
        args.metadata,
        priority,
        DEFAULT_PROCESSED_ROOTS,
    )
    ml_manifest = write_ml_command_manifest(out_dir, matrix_gate, args.metadata)
    decision = decide(env_state, matrix_gate, queue)

    inventory.to_csv(out_dir / "raw_inventory.csv", index=False)
    sample_run.to_csv(out_dir / "sample_run_manifest.csv", index=False)
    pairing_qc.to_csv(out_dir / "pairing_qc_manifest.csv", index=False)
    queue.to_csv(out_dir / "raw_etl_queue.csv", index=False)
    matrix_gate.to_csv(out_dir / "matrix_gate_table.csv", index=False)
    write_json(out_dir / "environment_gate.json", env_state)
    watchdog_state = {
        "timestamp": utc_now(),
        "status": "not_started",
        "pid": None,
        "current_sample": None,
        "retry_count": 0,
        "exit_code": None,
        "recommended_parallel_samples": 3,
        "queue_path": rel(out_dir / "raw_etl_queue.csv"),
    }
    write_json(out_dir / "watchdog_state.json", watchdog_state)
    write_json(out_dir / "ralph_decision_state.json", decision)
    write_report(out_dir, inventory, sample_run, pairing_qc, queue, env_state, matrix_gate, decision)

    summary = {
        **decision,
        "n_fastq_files": int(len(inventory)),
        "n_queue_rows": int(len(queue)),
        "environment_gate_passed": bool(env_state.get("environment_gate_passed")),
        "n_ml_commands_prepared": int(len(ml_manifest)),
        "output_dir": str(out_dir),
    }
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Run the v11.4 GSE281602 processed-COV conversion smoke.

This is a data-gate step, not a benchmark. GSE281602 is cardiomyocyte /
cell-type specific, so the resulting matrix remains auxiliary even if parsing
and region overlap pass.
"""
from __future__ import annotations

import json
import re
import subprocess
import tarfile
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")
PYTHON = ROOT.parent / "as-ds-ops" / ".venv" / "bin" / "python"
OUT_DIR = ROOT / "results" / "v11_4_conversions" / "GSE281602"
VERIFIED_FILES = OUT_DIR / "GSE281602_verified_files.csv"
TAR_PATH = OUT_DIR / "GSE281602_verified_processed_files.tar"
METADATA_PATH = OUT_DIR / "GSE281602_auxiliary_metadata.csv"
REFERENCE_MATRIX = ROOT / "results" / "multidataset_v8_2_prefilter_liftover" / "all_rrbs_region_matrix_5kb.parquet"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def parse_age_sex_from_title(title: str) -> tuple[float | None, float | None, str, str]:
    title = str(title)
    age_match = re.search(r"(\d+(?:\.\d+)?)\s*-\s*mo", title, flags=re.IGNORECASE)
    sex = "unknown"
    if re.search(r"\bfemale\b", title, flags=re.IGNORECASE):
        sex = "F"
    elif re.search(r"\bmale\b", title, flags=re.IGNORECASE):
        sex = "M"
    if not age_match:
        return None, None, sex, ""
    months = float(age_match.group(1))
    age_days = months * 30.42
    return age_days, age_days / 7.0, sex, f"{age_match.group(1)}mo"


def build_auxiliary_metadata() -> pd.DataFrame:
    samples = pd.read_csv(ROOT / "metadata" / "geo_old_age_candidate_samples.csv")
    rows = samples[samples["dataset"].astype(str).eq("GSE281602")].copy()
    if rows.empty:
        raise RuntimeError("No GSE281602 rows found in geo_old_age_candidate_samples.csv")
    parsed = []
    for row in rows.to_dict(orient="records"):
        age_days, age_weeks, sex, raw_age_token = parse_age_sex_from_title(str(row.get("title", "")))
        parsed.append(
            {
                "sample_id": row["sample_id"],
                "age_days": age_days,
                "age_weeks": age_weeks,
                "tissue": "heart",
                "sex": sex,
                "strain": "wildtype",
                "intervention": "control",
                "dataset_batch": "GSE281602",
                "assay": "RRBS_processed_bismark_cov",
                "metadata_source": "v11_4_geo_soft_header_auxiliary",
                "sample_id_source": "GEO_GSM_processed_cov",
                "cell_type": "cardiomyocytes",
                "raw_age_token": raw_age_token,
                "title": row.get("title", ""),
                "source_name": row.get("source_name", ""),
                "raw_characteristics": row.get("characteristics", ""),
                "headline_allowed": False,
                "auxiliary_only_reason": "cardiomyocyte_cell_type_specific_not_bulk_heart",
            }
        )
    meta = pd.DataFrame(parsed).sort_values("sample_id")
    meta.to_csv(METADATA_PATH, index=False)
    return meta


def make_tar() -> dict:
    files = pd.read_csv(VERIFIED_FILES)
    files = files[files["hard_verified"].astype(bool)].copy()
    if files.empty:
        raise RuntimeError(f"No hard-verified files in {VERIFIED_FILES}")
    missing = [path for path in files["local_path"].astype(str) if not Path(path).exists()]
    if missing:
        raise RuntimeError(f"Missing verified local files: {missing[:5]}")
    with tarfile.open(TAR_PATH, "w") as tar:
        for path_text in files["local_path"].astype(str):
            path = Path(path_text)
            tar.add(path, arcname=path.name)
    return {
        "tar_path": str(TAR_PATH),
        "n_files": int(len(files)),
        "tar_size_bytes": int(TAR_PATH.stat().st_size),
        "input_bytes": int(files["actual_size_bytes"].sum()),
    }


def run_conversion() -> dict:
    cmd = [
        str(PYTHON),
        "scripts/etl/10_convert_processed_methylation.py",
        "--input",
        str(TAR_PATH),
        "--dataset",
        "GSE281602",
        "--schema",
        "bismark_cov_per_sample_tar",
        "--min_coverage",
        "5",
        "--min_sample_presence",
        "0.5",
        "--region_min_sample_presence",
        "0.8",
        "--metadata_path",
        str(METADATA_PATH),
        "--out_dir",
        str(OUT_DIR),
    ]
    proc = subprocess.run(cmd, cwd=str(ROOT), text=True, capture_output=True)
    (OUT_DIR / "conversion_stdout.txt").write_text(proc.stdout, encoding="utf-8")
    (OUT_DIR / "conversion_stderr.txt").write_text(proc.stderr, encoding="utf-8")
    if proc.returncode != 0:
        return {
            "status": "blocked",
            "reason": "conversion_failed",
            "returncode": proc.returncode,
            "stderr_tail": proc.stderr[-4000:],
        }
    manifest_path = OUT_DIR / "GSE281602_matrix_manifest.json"
    if not manifest_path.exists():
        return {"status": "blocked", "reason": "matrix_manifest_missing"}
    return json.loads(manifest_path.read_text(encoding="utf-8"))


def compute_overlap_and_qc(conversion: dict) -> dict:
    region_path = Path(str(conversion.get("region_matrix_path", "")))
    beta_path = Path(str(conversion.get("beta_matrix_path", "")))
    if not region_path.exists() or not beta_path.exists():
        return {"status": "blocked", "reason": "converted_matrix_missing"}
    beta = pd.read_parquet(beta_path)
    region = pd.read_parquet(region_path)
    beta_min = float(beta.min(skipna=True).min(skipna=True))
    beta_max = float(beta.max(skipna=True).max(skipna=True))
    beta_bad_values = int(((beta < 0) | (beta > 1)).sum(skipna=True).sum())
    qc = {
        "status": "completed",
        "beta_shape": [int(beta.shape[0]), int(beta.shape[1])],
        "region_shape": [int(region.shape[0]), int(region.shape[1])],
        "beta_min": beta_min,
        "beta_max": beta_max,
        "beta_values_outside_0_1": beta_bad_values,
        "n_samples": int(beta.shape[1]),
    }
    if REFERENCE_MATRIX.exists():
        reference = pd.read_parquet(REFERENCE_MATRIX)
        common = sorted(set(region.index) & set(reference.index))
        pd.DataFrame({"region_id": common}).to_csv(OUT_DIR / "common_regions_with_v8_2.csv", index=False)
        qc.update(
            {
                "reference_matrix": str(REFERENCE_MATRIX),
                "reference_regions": int(reference.shape[0]),
                "common_regions_with_v8_2": int(len(common)),
                "common_region_fraction_vs_gse281602": float(len(common) / region.shape[0]) if region.shape[0] else 0.0,
                "common_region_fraction_vs_reference": float(len(common) / reference.shape[0]) if reference.shape[0] else 0.0,
                "common_region_gate_50000": bool(len(common) >= 50_000),
            }
        )
    else:
        qc.update(
            {
                "reference_matrix": str(REFERENCE_MATRIX),
                "common_region_gate_50000": False,
                "overlap_reason": "reference_matrix_missing",
            }
        )
    return qc


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    meta = build_auxiliary_metadata()
    tar_state = make_tar()
    conversion = run_conversion()
    qc = compute_overlap_and_qc(conversion) if conversion.get("status", "").startswith("completed") else {}
    state = {
        "timestamp": utc_now(),
        "loop_version": "v11.4",
        "dataset": "GSE281602",
        "status": "conversion_smoke_completed" if qc.get("status") == "completed" else "blocked",
        "decision": "auxiliary_matrix_only_not_headline",
        "metadata_path": str(METADATA_PATH),
        "n_metadata_rows": int(len(meta)),
        "age_weeks_values": sorted(float(x) for x in meta["age_weeks"].dropna().unique()),
        "tissue_counts": meta["tissue"].value_counts().to_dict(),
        "sex_counts": meta["sex"].value_counts().to_dict(),
        "tar": tar_state,
        "conversion_manifest": conversion,
        "qc": qc,
        "headline_allowed": False,
        "reason_not_headline": "GSE281602 is cardiomyocyte/cell-type specific, not bulk heart; use as auxiliary conversion/overlap gate only.",
    }
    write_json(OUT_DIR / "v11_4_gse281602_conversion_smoke_state.json", state)
    print(json.dumps(state, indent=2)[:8000])


if __name__ == "__main__":
    main()

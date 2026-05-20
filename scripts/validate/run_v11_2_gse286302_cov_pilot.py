#!/usr/bin/env python3
"""Run the v11.2 GSE286302 processed-COV diagnostic pilot.

This is not a headline benchmark. It downloads a tiny processed-COV subset,
builds a 5kb region matrix, and measures overlap with the current v8.2
multidataset region matrix. It avoids FASTQ/Bismark because the local alignment
toolchain is not currently available and processed COV supplements exist.
"""
from __future__ import annotations

import json
import math
import shutil
import subprocess
import tarfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")
OUT_DIR = ROOT / "results" / "v11_2_gse286302_cov_pilot"
DOWNLOAD_DIR = ROOT / "raw_downloads" / "geo_supplements" / "GSE286302_v11_2"
SUPPLEMENTS = ROOT / "metadata" / "geo_old_age_candidate_supplements.csv"
REFERENCE_MATRIX = ROOT / "results" / "multidataset_v8_2_prefilter_liftover" / "all_rrbs_region_matrix_5kb.parquet"
DOWNLOAD_CHUNKS = 8

SELECTED_SAMPLES = [
    {"sample_id": "GSM8723148", "age_group": "young", "condition_detail": "lung, young, rep1"},
    {"sample_id": "GSM8723152", "age_group": "old", "condition_detail": "lung, aged, vehicle, rep1"},
]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def geo_sample_supplement_url(sample_id: str, filename: str) -> str:
    """Build the official GEO sample-level supplement URL.

    GEO stores sample supplements under a grouped GSM directory where the last
    three digits are replaced by nnn, e.g. GSM8723148 -> GSM8723nnn.
    """
    digits = sample_id.replace("GSM", "")
    if not digits.isdigit() or len(digits) <= 3:
        return ""
    group = f"GSM{digits[:-3]}nnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{group}/{sample_id}/suppl/{filename}"


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def append_jsonl(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(payload, sort_keys=True) + "\n")


def build_manifest() -> pd.DataFrame:
    supp = pd.read_csv(SUPPLEMENTS)
    rows = []
    for selected in SELECTED_SAMPLES:
        sample_id = selected["sample_id"]
        hit = supp[
            (supp["dataset"].astype(str) == "GSE286302")
            & (supp["supplement_name"].astype(str).str.startswith(sample_id))
            & (supp["supplement_name"].astype(str).str.endswith(".bismark.cov.gz"))
        ].copy()
        if hit.empty:
            rows.append({**selected, "status": "blocked_missing_supplement"})
            continue
        item = hit.iloc[0].to_dict()
        remote_size = item.get("remote_size_bytes")
        if pd.isna(remote_size):
            remote_size = item.get("filelist_size_bytes")
        if pd.isna(remote_size):
            remote_size = item.get("size")
        if pd.isna(remote_size):
            remote_size = 0
        supplement_name = str(item.get("supplement_name", ""))
        inventory_url = item.get("supplement_url") or item.get("url") or ""
        url = geo_sample_supplement_url(sample_id, supplement_name) or inventory_url
        rows.append(
            {
                **selected,
                "dataset_batch": "GSE286302",
                "tissue": "lung",
                "intervention": "control",
                "assay": "RRBS_processed_bismark_cov",
                "metadata_source": "v11_2_processed_cov_diagnostic_pilot",
                "supplement_name": supplement_name,
                "inventory_url": inventory_url,
                "url": url,
                "remote_size_bytes": int(float(remote_size or 0)),
                "download_path": str(DOWNLOAD_DIR / str(item.get("supplement_name", ""))),
                "status": "planned",
            }
        )
    manifest = pd.DataFrame(rows)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(OUT_DIR / "gse286302_cov_pilot_manifest.csv", index=False)
    metadata = manifest[manifest["status"] == "planned"].copy()
    metadata["age_days"] = pd.NA
    metadata["age_weeks"] = pd.NA
    metadata["sex"] = "unknown"
    metadata["strain"] = "unknown"
    metadata["sample_id_source"] = "GEO_GSM_processed_cov"
    metadata[
        [
            "sample_id",
            "age_days",
            "age_weeks",
            "tissue",
            "sex",
            "strain",
            "intervention",
            "dataset_batch",
            "assay",
            "metadata_source",
            "age_group",
            "condition_detail",
            "sample_id_source",
        ]
    ].to_csv(OUT_DIR / "gse286302_cov_pilot_metadata.csv", index=False)
    return manifest


def download_with_curl_ranges(url: str, path: Path, expected: int, chunks: int = DOWNLOAD_CHUNKS) -> tuple[bool, str]:
    if expected <= 0 or not url:
        return False, "missing_expected_size_or_url"
    parts_dir = path.with_name(path.name + ".parts")
    parts_dir.mkdir(parents=True, exist_ok=True)
    chunk_ranges = []
    chunk_size = math.ceil(expected / chunks)
    for idx in range(chunks):
        start = idx * chunk_size
        end = min(expected - 1, (idx + 1) * chunk_size - 1)
        if start > end:
            continue
        chunk_ranges.append((idx, start, end))

    def fetch_chunk(item: tuple[int, int, int]) -> tuple[int, bool, str]:
        idx, start, end = item
        chunk_path = parts_dir / f"part_{idx:03d}"
        expected_len = end - start + 1
        if chunk_path.exists() and chunk_path.stat().st_size == expected_len:
            return idx, True, "already_verified"
        details = []
        stalled_attempts = 0
        for _ in range(12):
            current = chunk_path.stat().st_size if chunk_path.exists() else 0
            if current == expected_len:
                return idx, True, "resumed_verified"
            if current > expected_len:
                chunk_path.unlink()
                current = 0
            resume_start = start + current
            tmp_path = parts_dir / f"part_{idx:03d}.resume"
            tmp_path.unlink(missing_ok=True)
            cmd = [
                "curl",
                "-L",
                "--fail",
                "--retry",
                "3",
                "--retry-delay",
                "5",
                "--range",
                f"{resume_start}-{end}",
                "-o",
                str(tmp_path),
                url,
            ]
            proc = subprocess.run(cmd, text=True, capture_output=True)
            added = tmp_path.stat().st_size if tmp_path.exists() else 0
            if added:
                with chunk_path.open("ab") as out_handle, tmp_path.open("rb") as in_handle:
                    shutil.copyfileobj(in_handle, out_handle)
                stalled_attempts = 0
            else:
                stalled_attempts += 1
            tmp_path.unlink(missing_ok=True)
            details.append(proc.stderr[-400:])
            if chunk_path.exists() and chunk_path.stat().st_size == expected_len:
                return idx, True, "resumed_verified"
            if proc.returncode == 0 and not added:
                break
            if stalled_attempts >= 3:
                break
        actual = chunk_path.stat().st_size if chunk_path.exists() else 0
        return idx, False, f"chunk_incomplete:{actual}:{expected_len}:{''.join(details[-2:])}"

    failures = []
    with ThreadPoolExecutor(max_workers=min(chunks, len(chunk_ranges))) as executor:
        future_to_chunk = {executor.submit(fetch_chunk, item): item for item in chunk_ranges}
        for future in as_completed(future_to_chunk):
            idx, ok, detail = future.result()
            if not ok:
                failures.append({"chunk": idx, "detail": detail})
    if failures:
        return False, json.dumps(failures[:3], sort_keys=True)

    tmp_path = path.with_name(path.name + ".tmp")
    with tmp_path.open("wb") as out_handle:
        for idx, _, _ in chunk_ranges:
            chunk_path = parts_dir / f"part_{idx:03d}"
            with chunk_path.open("rb") as in_handle:
                shutil.copyfileobj(in_handle, out_handle)
    assembled_size = tmp_path.stat().st_size
    if assembled_size != expected:
        tmp_path.unlink(missing_ok=True)
        return False, f"assembled_size_mismatch:{assembled_size}:{expected}"
    tmp_path.replace(path)
    shutil.rmtree(parts_dir, ignore_errors=True)
    return True, "range_download_completed"


def download_files(manifest: pd.DataFrame) -> list[Path]:
    DOWNLOAD_DIR.mkdir(parents=True, exist_ok=True)
    paths = []
    log_path = OUT_DIR / "download_log.jsonl"
    for row in manifest.to_dict(orient="records"):
        if row.get("status") != "planned":
            append_jsonl(log_path, {"timestamp": utc_now(), "sample_id": row.get("sample_id"), "status": row.get("status")})
            continue
        path = Path(str(row["download_path"]))
        expected = int(row.get("remote_size_bytes") or 0)
        if path.exists() and expected and path.stat().st_size == expected:
            status = "already_verified"
        else:
            ok, detail = download_with_curl_ranges(str(row["url"]), path, expected)
            status = "downloaded_range" if ok else "download_failed"
            if not ok:
                append_jsonl(
                    log_path,
                    {
                        "timestamp": utc_now(),
                        "sample_id": row.get("sample_id"),
                        "status": status,
                        "url": row.get("url"),
                        "stderr": detail,
                    },
                )
                continue
        size = path.stat().st_size if path.exists() else 0
        verified = bool(expected and size == expected)
        append_jsonl(
            log_path,
            {
                "timestamp": utc_now(),
                "sample_id": row.get("sample_id"),
                "status": status,
                "path": str(path),
                "expected_size": expected,
                "actual_size": size,
                "verified_size": verified,
            },
        )
        if verified or path.exists():
            paths.append(path)
    return paths


def make_tar(paths: list[Path]) -> Path:
    tar_path = OUT_DIR / "GSE286302_v11_2_cov_pilot.tar"
    with tarfile.open(tar_path, "w") as tar:
        for path in paths:
            tar.add(path, arcname=path.name)
    return tar_path


def run_conversion(tar_path: Path) -> dict:
    cmd = [
        str(ROOT.parent / "as-ds-ops" / ".venv" / "bin" / "python"),
        "scripts/etl/10_convert_processed_methylation.py",
        "--input",
        str(tar_path),
        "--dataset",
        "GSE286302",
        "--schema",
        "bismark_cov_per_sample_tar",
        "--min_coverage",
        "5",
        "--min_sample_presence",
        "0.5",
        "--region_min_sample_presence",
        "0.67",
        "--metadata_path",
        str(OUT_DIR / "gse286302_cov_pilot_metadata.csv"),
        "--out_dir",
        str(OUT_DIR),
    ]
    proc = subprocess.run(cmd, cwd=str(ROOT), text=True, capture_output=True)
    (OUT_DIR / "conversion_stdout.txt").write_text(proc.stdout, encoding="utf-8")
    (OUT_DIR / "conversion_stderr.txt").write_text(proc.stderr, encoding="utf-8")
    if proc.returncode != 0:
        return {"status": "blocked", "reason": "conversion_failed", "stderr": proc.stderr[-2000:]}
    manifest_path = OUT_DIR / "GSE286302_matrix_manifest.json"
    return json.loads(manifest_path.read_text()) if manifest_path.exists() else {"status": "blocked", "reason": "manifest_missing"}


def compute_overlap(conversion_manifest: dict) -> dict:
    region_path = Path(conversion_manifest.get("region_matrix_path", ""))
    if not region_path.exists() or not REFERENCE_MATRIX.exists():
        return {"status": "blocked", "reason": "region_or_reference_matrix_missing"}
    pilot = pd.read_parquet(region_path)
    reference = pd.read_parquet(REFERENCE_MATRIX)
    common = sorted(set(pilot.index) & set(reference.index))
    overlap = {
        "status": "completed",
        "pilot_regions": int(pilot.shape[0]),
        "reference_regions": int(reference.shape[0]),
        "common_regions": int(len(common)),
        "common_region_fraction_vs_pilot": float(len(common) / pilot.shape[0]) if pilot.shape[0] else 0.0,
        "common_region_fraction_vs_reference": float(len(common) / reference.shape[0]) if reference.shape[0] else 0.0,
        "pilot_samples": list(pilot.columns),
        "reference_matrix": str(REFERENCE_MATRIX),
        "pilot_region_matrix": str(region_path),
    }
    pd.DataFrame({"region_id": common}).to_csv(OUT_DIR / "common_regions_with_v8_2.csv", index=False)
    return overlap


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    tool_status = {
        "wget": shutil.which("wget"),
        "bismark": shutil.which("bismark"),
        "bowtie2": shutil.which("bowtie2"),
        "samtools": shutil.which("samtools"),
        "mode": "processed_cov_no_alignment",
    }
    manifest = build_manifest()
    paths = download_files(manifest)
    if len(paths) < 2:
        state = {"status": "blocked", "reason": "fewer_than_two_cov_files_downloaded", "tool_status": tool_status}
        write_json(OUT_DIR / "v11_2_decision_state.json", state)
        print(json.dumps(state, indent=2))
        return
    tar_path = make_tar(paths)
    conversion = run_conversion(tar_path)
    overlap = compute_overlap(conversion) if conversion.get("status", "").startswith("completed") else {}
    state = {
        "timestamp": utc_now(),
        "loop_version": "v11.2",
        "status": "diagnostic_completed" if overlap.get("status") == "completed" else "blocked",
        "decision": "processed_cov_diagnostic_overlap_only_not_headline",
        "tool_status": tool_status,
        "n_files_downloaded": len(paths),
        "tar_path": str(tar_path),
        "conversion_manifest": conversion,
        "overlap": overlap,
        "headline_allowed": False,
        "reason_not_headline": "GSE286302 pilot has age groups only and no exact age_weeks; it is old-lung diagnostic only.",
    }
    write_json(OUT_DIR / "v11_2_decision_state.json", state)
    print(json.dumps(state, indent=2)[:6000])


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Lightweight schema smoke test for v11.3 verified processed supplements."""
from __future__ import annotations

import argparse
import gzip
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")
DEFAULT_COMPLETION = ROOT / "results" / "download_backlog_v11_3" / "download_completion_summary.csv"
DEFAULT_OUT = ROOT / "results" / "download_backlog_v11_3"
EXCLUDE_CHROMS = {"chrX", "chrY", "chrM", "X", "Y", "M", "MT"}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def read_head_lines(path: Path, max_lines: int) -> tuple[list[str], str]:
    opener = gzip.open if str(path).endswith(".gz") else open
    lines: list[str] = []
    try:
        with opener(path, "rt", errors="replace") as handle:
            for _ in range(max_lines):
                line = handle.readline()
                if not line:
                    break
                lines.append(line.rstrip("\n"))
        return lines, ""
    except Exception as exc:  # noqa: BLE001 - smoke must report any parser failure
        return lines, str(exc)[:500]


def normalize_chrom(value: str) -> str | None:
    value = str(value).strip()
    if value in EXCLUDE_CHROMS:
        return None
    if value.startswith("chr"):
        return value
    if value.isdigit():
        return f"chr{value}"
    return None


def is_accession_like_chrom(value: str) -> bool:
    value = str(value).strip()
    prefixes = ("CM", "GL", "JH", "KB", "NC_")
    return value.startswith(prefixes) and any(char.isdigit() for char in value)


def is_float(value: str) -> bool:
    try:
        float(value)
        return True
    except ValueError:
        return False


def detect_schema(lines: list[str]) -> dict:
    data_lines = [line for line in lines if line and not line.startswith("#")]
    if not data_lines:
        return {"schema_detected": "empty_or_comment_only", "schema_pass": False}
    first = data_lines[0]
    sep = "\t" if "\t" in first else None
    if sep is None:
        return {"schema_detected": "unknown_no_tabs", "schema_pass": False, "first_data_line": first[:200]}
    header = first.split("\t")
    if header[:12] == [
        "chr",
        "pos",
        "strand",
        "context",
        "ratio",
        "eff_CT_count",
        "C_count",
        "CT_count",
        "rev_G_count",
        "rev_GA_count",
        "CI_lower",
        "CI_upper",
    ]:
        body = [line.split("\t") for line in data_lines[1:201]]
        chrom_ok = 0
        numeric_rows = 0
        ratio_values = []
        coverage_values = []
        for row in body:
            if len(row) < 8:
                continue
            if normalize_chrom(row[0]) is not None:
                chrom_ok += 1
            if all(is_float(row[idx]) for idx in [1, 4, 6, 7]):
                numeric_rows += 1
                ratio_values.append(float(row[4]))
                coverage_values.append(float(row[7]))
        schema_pass = bool(
            body
            and numeric_rows
            and ratio_values
            and min(ratio_values) >= 0
            and max(ratio_values) <= 1
        )
        return {
            "schema_detected": "methylratio_cg_12col",
            "schema_pass": schema_pass,
            "n_data_lines_checked": len(body),
            "min_cols": min((len(row) for row in body), default=0),
            "max_cols": max((len(row) for row in body), default=0),
            "chrom_ok_rows": chrom_ok,
            "numeric_position_rows": numeric_rows,
            "ratio_min": min(ratio_values) if ratio_values else None,
            "ratio_max": max(ratio_values) if ratio_values else None,
            "coverage_min": min(coverage_values) if coverage_values else None,
            "coverage_max": max(coverage_values) if coverage_values else None,
            "first_data_line": first[:200],
        }
    rows = [line.split("\t") for line in data_lines[:200]]
    ncols = [len(row) for row in rows]
    chrom_ok = 0
    accession_chrom_rows = 0
    numeric_positions = 0
    beta_like = 0
    cov_like = 0
    pct_values = []
    beta_values = []
    coverage_values = []
    for row in rows:
        if len(row) >= 1 and normalize_chrom(row[0]) is not None:
            chrom_ok += 1
        if len(row) >= 1 and is_accession_like_chrom(row[0]):
            accession_chrom_rows += 1
        if len(row) >= 3 and is_float(row[1]) and is_float(row[2]):
            numeric_positions += 1
        if len(row) >= 6 and all(is_float(row[idx]) for idx in [1, 2, 3, 4, 5]):
            pct = float(row[3])
            methylated = float(row[4])
            unmethylated = float(row[5])
            pct_values.append(pct)
            coverage_values.append(methylated + unmethylated)
            cov_like += 1
        if len(row) >= 4 and is_float(row[3]):
            beta = float(row[3])
            if 0 <= beta <= 1:
                beta_like += 1
                beta_values.append(beta)
    min_cols = min(ncols)
    max_cols = max(ncols)
    if cov_like and (chrom_ok or accession_chrom_rows) and numeric_positions:
        pct_min = min(pct_values) if pct_values else None
        pct_max = max(pct_values) if pct_values else None
        coverage_min = min(coverage_values) if coverage_values else None
        coverage_max = max(coverage_values) if coverage_values else None
        schema_pass = pct_min is not None and 0 <= pct_min <= 100 and pct_max is not None and 0 <= pct_max <= 100 and coverage_max is not None and coverage_max > 0
        schema_detected = "bismark_cov_6col" if chrom_ok else "bismark_cov_6col_accession_chrom"
        return {
            "schema_detected": schema_detected,
            "schema_pass": bool(schema_pass),
            "n_data_lines_checked": len(rows),
            "min_cols": min_cols,
            "max_cols": max_cols,
            "chrom_ok_rows": chrom_ok,
            "accession_chrom_rows": accession_chrom_rows,
            "numeric_position_rows": numeric_positions,
            "cov_like_rows": cov_like,
            "pct_min": pct_min,
            "pct_max": pct_max,
            "coverage_min": coverage_min,
            "coverage_max": coverage_max,
            "first_data_line": first[:200],
        }
    if beta_like and chrom_ok:
        return {
            "schema_detected": "bed_or_matrix_beta_like",
            "schema_pass": True,
            "n_data_lines_checked": len(rows),
            "min_cols": min_cols,
            "max_cols": max_cols,
            "chrom_ok_rows": chrom_ok,
            "beta_like_rows": beta_like,
            "beta_min": min(beta_values) if beta_values else None,
            "beta_max": max(beta_values) if beta_values else None,
            "first_data_line": first[:200],
        }
    return {
        "schema_detected": "unknown_tabular",
        "schema_pass": False,
        "n_data_lines_checked": len(rows),
        "min_cols": min_cols,
        "max_cols": max_cols,
        "chrom_ok_rows": chrom_ok,
        "accession_chrom_rows": accession_chrom_rows,
        "numeric_position_rows": numeric_positions,
        "first_data_line": first[:200],
    }


def smoke_file(row: dict, max_lines: int) -> dict:
    path = Path(str(row["local_path"]))
    lines, error = read_head_lines(path, max_lines=max_lines)
    result = {
        "dataset": row["dataset"],
        "sample_id": row.get("sample_id", ""),
        "supplement_name": row["supplement_name"],
        "file_type": row["file_type"],
        "path": str(path),
        "read_error": error,
        "n_lines_read": len(lines),
        "gzip_readable": not bool(error) and len(lines) > 0,
    }
    if error:
        result.update({"schema_detected": "read_error", "schema_pass": False})
        return result
    result.update(detect_schema(lines))
    return result


def dataset_gate(smoke: pd.DataFrame, completion: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for dataset, ds in completion.groupby("dataset", sort=True):
        smoke_ds = smoke[smoke["dataset"] == dataset]
        schema_counts = smoke_ds["schema_detected"].value_counts().to_dict() if not smoke_ds.empty else {}
        schema_pass = bool(smoke_ds["schema_pass"].fillna(False).all()) if not smoke_ds.empty else False
        schemas = set(smoke_ds["schema_detected"].astype(str))
        if schemas == {"bismark_cov_6col"} and schema_pass:
            gate = "ready_for_bismark_cov_conversion"
            adapter = "bismark_cov_per_sample_tar"
        elif schemas == {"bismark_cov_6col_accession_chrom"} and schema_pass:
            gate = "adapter_required_before_conversion"
            adapter = "bismark_cov_with_accession_chrom_map"
        elif schemas == {"methylratio_cg_12col"} and schema_pass:
            gate = "adapter_required_before_conversion"
            adapter = "methylratio_cg_12col_to_cov"
        elif schema_pass:
            gate = "adapter_required_before_conversion"
            adapter = ";".join(sorted(schemas))
        else:
            gate = "blocked_schema_smoke_failed"
            adapter = ";".join(sorted(schemas)) if schemas else "none"
        rows.append(
            {
                "dataset": dataset,
                "n_verified_files": int(len(ds)),
                "file_types": ";".join(sorted(set(ds["file_type"].astype(str)))),
                "smoke_files_checked": int(len(smoke_ds)),
                "schema_counts": json.dumps(schema_counts, sort_keys=True),
                "all_smoke_pass": schema_pass,
                "conversion_gate": gate,
                "recommended_adapter": adapter,
                "headline_allowed_now": bool(ds["headline_allowed"].fillna(False).astype(bool).any()),
                "tier": ";".join(sorted(set(ds["tier"].astype(str)))),
                "reason": " | ".join(sorted(set(ds["reason"].astype(str))))[:1000],
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--completion", type=Path, default=DEFAULT_COMPLETION)
    parser.add_argument("--out_dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--samples_per_dataset", type=int, default=3)
    parser.add_argument("--max_lines", type=int, default=1000)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    completion = pd.read_csv(args.completion)
    verified = completion[completion["hard_verified"]].copy()
    sampled_rows = []
    for dataset, ds in verified.groupby("dataset", sort=True):
        sample_n = min(args.samples_per_dataset, len(ds))
        sampled_rows.extend(ds.sort_values("supplement_name").head(sample_n).to_dict(orient="records"))
    smoke_rows = [smoke_file(row, args.max_lines) for row in sampled_rows]
    smoke = pd.DataFrame(smoke_rows)
    smoke_path = args.out_dir / "schema_smoke_summary.csv"
    smoke.to_csv(smoke_path, index=False)
    gate = dataset_gate(smoke, verified)
    gate_path = args.out_dir / "dataset_download_gate_summary.csv"
    gate.to_csv(gate_path, index=False)
    summary = {
        "timestamp": utc_now(),
        "status": "completed",
        "completion_path": str(args.completion),
        "schema_smoke_summary_path": str(smoke_path),
        "dataset_gate_summary_path": str(gate_path),
        "datasets_ready_for_bismark_cov_conversion": gate.loc[
            gate["conversion_gate"] == "ready_for_bismark_cov_conversion",
            "dataset",
        ].tolist(),
        "datasets_requiring_adapter": gate.loc[
            gate["conversion_gate"] == "adapter_required_before_conversion",
            "dataset",
        ].tolist(),
        "datasets_blocked": gate.loc[
            gate["conversion_gate"] == "blocked_schema_smoke_failed",
            "dataset",
        ].tolist(),
    }
    (args.out_dir / "schema_smoke_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

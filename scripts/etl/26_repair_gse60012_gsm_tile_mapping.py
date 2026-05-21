#!/usr/bin/env python3
"""Repair the GSE60012 GSM-to-tile metadata crosswalk.

The existing GSE60012 training matrix uses synthetic tile IDs parsed from the
processed matrix header, while the local raw holdings are keyed by GEO GSM/SRR.
This script builds a deterministic crosswalk by matching tissue/sex/age and
condition tokens, then ordering duplicated samples within each matched group.
"""
from __future__ import annotations

import argparse
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")
DEFAULT_HEADER_METADATA = ROOT / "metadata" / "gse60012_header_sample_metadata.csv"
DEFAULT_SRA_METADATA = Path("/data/mouse_methyl/raw/GSE60012/metadata/sra_metadata.tsv")
DEFAULT_FASTQ_DIR = Path("/data/mouse_methyl/raw/GSE60012/fastq")
DEFAULT_OUT_DIR = ROOT / "results" / "v27_data_repair_execution" / "gse60012_mapping"

RRBS_MOUSE_RE = re.compile(
    r"(?P<gsm>GSM\d+):\s*RRBS_(?P<rep>[^_]+)_Mouse_(?P<tissue>[^_]+)_(?P<sex>[FM])_(?P<age>\d+(?:\.\d+)?w)_(?P<condition>.+?);",
    flags=re.IGNORECASE,
)
AGE_RE = re.compile(r"^(?P<value>\d+(?:\.\d+)?)w$", flags=re.IGNORECASE)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def clean_condition(value: object) -> str:
    text = "" if pd.isna(value) else str(value).strip()
    text = text.replace("vehacle", "vehicle")
    text = text.replace("castaration", "castration")
    text = re.sub(r"\s*\+\s*", "_+_", text)
    text = re.sub(r"\s+", "_", text)
    text = re.sub(r"_+", "_", text)
    return text.strip("_")


def age_text_to_weeks(value: str) -> float | None:
    match = AGE_RE.match(str(value).strip())
    if not match:
        return None
    return float(match.group("value"))


def key_for(tissue: object, sex: object, age_weeks: object, condition: object) -> str:
    try:
        age = float(age_weeks)
    except (TypeError, ValueError):
        age_token = ""
    else:
        age_token = f"{age:g}w"
    return "|".join(
        [
            str(tissue or "").strip().lower(),
            str(sex or "").strip().upper(),
            age_token,
            clean_condition(condition).lower(),
        ]
    )


def run_fastqs_complete(run_accession: str, fastq_dir: Path, layout: str) -> bool:
    layout = str(layout or "").upper()
    if layout == "PAIRED":
        return bool(list(fastq_dir.glob(f"{run_accession}_1.fastq*")) and list(fastq_dir.glob(f"{run_accession}_2.fastq*")))
    return bool(list(fastq_dir.glob(f"{run_accession}*.fastq*")))


def load_raw_metadata(path: Path, fastq_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    sra = pd.read_csv(path)
    parsed: list[dict[str, Any]] = []
    excluded: list[dict[str, Any]] = []
    for row in sra.to_dict(orient="records"):
        title = str(row.get("experiment_title") or row.get("experiment_desc") or "")
        match = RRBS_MOUSE_RE.search(title)
        if not match:
            excluded.append(
                {
                    "run_accession": row.get("run_accession", ""),
                    "experiment_title": title,
                    "library_strategy": row.get("library_strategy", ""),
                    "organism_name": row.get("organism_name", ""),
                    "exclude_reason": "not_mouse_rrbs_title_schema",
                }
            )
            continue
        groups = match.groupdict()
        age_weeks = age_text_to_weeks(groups["age"])
        if age_weeks is None:
            excluded.append(
                {
                    "run_accession": row.get("run_accession", ""),
                    "experiment_title": title,
                    "exclude_reason": "unparseable_age",
                }
            )
            continue
        parsed.append(
            {
                "gsm_sample_id": groups["gsm"],
                "run_accession": row.get("run_accession", ""),
                "experiment_accession": row.get("experiment_accession", ""),
                "sample_accession": row.get("sample_accession", ""),
                "biosample": row.get("biosample", ""),
                "replicate_token": groups["rep"],
                "replicate_order": int(groups["rep"]) if str(groups["rep"]).isdigit() else 10_000,
                "raw_tissue": groups["tissue"],
                "sex": groups["sex"].upper(),
                "age_weeks": age_weeks,
                "condition_detail": clean_condition(groups["condition"]),
                "experiment_title": title,
                "library_layout": row.get("library_layout", ""),
                "local_fastq_complete": run_fastqs_complete(str(row.get("run_accession", "")), fastq_dir, str(row.get("library_layout", ""))),
            }
        )
    raw = pd.DataFrame(parsed)
    if raw.empty:
        return raw, pd.DataFrame(excluded)
    raw["match_key"] = [
        key_for(row["raw_tissue"], row["sex"], row["age_weeks"], row["condition_detail"])
        for row in raw.to_dict(orient="records")
    ]
    return raw, pd.DataFrame(excluded)


def load_header_metadata(path: Path) -> pd.DataFrame:
    meta = pd.read_csv(path)
    rows: list[dict[str, Any]] = []
    for row in meta.to_dict(orient="records"):
        rows.append(
            {
                **row,
                "header_tissue_raw": row.get("raw_source_name", ""),
                "header_condition_clean": clean_condition(row.get("condition_detail", "")),
                "match_key": key_for(row.get("raw_source_name", ""), row.get("sex", ""), row.get("age_weeks", ""), row.get("condition_detail", "")),
            }
        )
    return pd.DataFrame(rows)


def build_crosswalk(header: pd.DataFrame, raw: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    mapped: list[dict[str, Any]] = []
    unresolved: list[dict[str, Any]] = []
    raw_by_key = {
        key: group.sort_values(["replicate_order", "gsm_sample_id", "run_accession"]).reset_index(drop=True)
        for key, group in raw.groupby("match_key", dropna=False)
    }
    for key, header_group in header.sort_values(["sample_id", "matrix_sample_id"]).groupby("match_key", dropna=False):
        raw_group = raw_by_key.get(key, pd.DataFrame())
        header_group = header_group.reset_index(drop=True)
        if raw_group.empty:
            for row in header_group.to_dict(orient="records"):
                unresolved.append({**row, "mapping_status": "unresolved_no_raw_group", "raw_group_count": 0})
            continue
        raw_count = len(raw_group)
        header_count = len(header_group)
        if raw_count < header_count:
            for row in header_group.to_dict(orient="records"):
                unresolved.append({**row, "mapping_status": "unresolved_raw_count_lt_header_count", "raw_group_count": raw_count})
            continue
        confidence = "exact_group_count_ordered" if raw_count == header_count else "raw_group_has_extra_rows_ordered"
        for idx, row in header_group.iterrows():
            raw_row = raw_group.iloc[idx].to_dict()
            mapped.append(
                {
                    "tile_sample_id": row["sample_id"],
                    "matrix_sample_id": row.get("matrix_sample_id", ""),
                    "source_column": row.get("source_column", ""),
                    "match_key": key,
                    "gsm_sample_id": raw_row["gsm_sample_id"],
                    "run_accession": raw_row["run_accession"],
                    "experiment_accession": raw_row["experiment_accession"],
                    "sample_accession": raw_row["sample_accession"],
                    "biosample": raw_row["biosample"],
                    "replicate_token": raw_row["replicate_token"],
                    "replicate_order": raw_row["replicate_order"],
                    "local_fastq_complete": raw_row["local_fastq_complete"],
                    "mapping_status": "mapped",
                    "mapping_confidence": confidence,
                    "raw_group_count": raw_count,
                    "header_group_count": header_count,
                    "age_weeks": row.get("age_weeks", ""),
                    "tissue": row.get("tissue", ""),
                    "sex": row.get("sex", ""),
                    "condition_family": row.get("condition_family", ""),
                    "condition_detail": row.get("condition_detail", ""),
                    "raw_experiment_title": raw_row.get("experiment_title", ""),
                }
            )
    return pd.DataFrame(mapped), pd.DataFrame(unresolved)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--header-metadata", type=Path, default=DEFAULT_HEADER_METADATA)
    parser.add_argument("--sra-metadata", type=Path, default=DEFAULT_SRA_METADATA)
    parser.add_argument("--fastq-dir", type=Path, default=DEFAULT_FASTQ_DIR)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--min-map-fraction", type=float, default=0.95)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    header = load_header_metadata(args.header_metadata)
    raw, raw_excluded = load_raw_metadata(args.sra_metadata, args.fastq_dir)
    crosswalk, unresolved = build_crosswalk(header, raw)

    crosswalk_path = args.out_dir / "gse60012_gsm_tile_crosswalk.csv"
    unresolved_path = args.out_dir / "gse60012_unresolved_tile_rows.csv"
    raw_path = args.out_dir / "gse60012_parsed_raw_sra_rows.csv"
    raw_excluded_path = args.out_dir / "gse60012_raw_sra_exclusions.csv"
    repaired_metadata_path = args.out_dir / "gse60012_header_metadata_with_gsm_crosswalk.csv"
    summary_path = args.out_dir / "gse60012_mapping_summary.json"

    crosswalk.to_csv(crosswalk_path, index=False)
    unresolved.to_csv(unresolved_path, index=False)
    raw.to_csv(raw_path, index=False)
    raw_excluded.to_csv(raw_excluded_path, index=False)

    repaired = header.merge(crosswalk, left_on="sample_id", right_on="tile_sample_id", how="left", suffixes=("", "_crosswalk"))
    repaired["gsm_crosswalk_status"] = repaired["mapping_status"].fillna("unmapped")
    repaired.to_csv(repaired_metadata_path, index=False)

    map_fraction = 0.0 if header.empty else len(crosswalk) / len(header)
    summary = {
        "timestamp": utc_now(),
        "status": "passed" if map_fraction >= args.min_map_fraction else "blocked",
        "header_rows": int(len(header)),
        "mapped_rows": int(len(crosswalk)),
        "unresolved_rows": int(len(unresolved)),
        "raw_parsed_rows": int(len(raw)),
        "raw_excluded_rows": int(len(raw_excluded)),
        "map_fraction": round(float(map_fraction), 6),
        "min_map_fraction": float(args.min_map_fraction),
        "mapped_local_fastq_complete_rows": int(crosswalk["local_fastq_complete"].fillna(False).sum()) if not crosswalk.empty else 0,
        "crosswalk_path": str(crosswalk_path),
        "repaired_metadata_path": str(repaired_metadata_path),
        "unresolved_path": str(unresolved_path),
        "raw_excluded_path": str(raw_excluded_path),
        "policy": "mapping_only_no_training_authorization",
    }
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

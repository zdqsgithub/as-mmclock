#!/usr/bin/env python3
"""Build v7.1 model metadata with GSE60012 synthetic header-derived samples."""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")
BASE_METADATA = ROOT / "metadata" / "unified_sample_metadata.csv"
GSE60012_STATS = ROOT / "results" / "multidataset" / "GSE60012_sample_parse_stats.csv"
OUT_METADATA = ROOT / "metadata" / "gse60012_header_sample_metadata.csv"
OUT_MODEL_METADATA = ROOT / "metadata" / "model_sample_metadata_v7_1.csv"
OUT_EXCLUSIONS = ROOT / "results" / "multidataset_v7_1" / "GSE60012_header_mapping_exclusions.csv"

SHORT_TILE_RE = re.compile(r"^(GSE60012_tile_\d{3})")
AGE_RE = re.compile(r"^(\d+(?:\.\d+)?)w$", re.IGNORECASE)

TISSUE_MAP = {
    "cerebellum": "brain_cerebellum",
    "liver": "liver",
    "muscle": "muscle",
    "spleen": "spleen",
}


def clean_source_label(value: object) -> str:
    text = "" if pd.isna(value) else str(value).strip()
    text = text.replace("vehacle", "vehicle")
    text = text.replace("castaration", "castration")
    text = re.sub(r"\s*\+\s*", "_+_", text)
    text = re.sub(r"\s+", "_", text)
    text = re.sub(r"_+", "_", text)
    return text.strip("_")


def short_sample_id(value: str) -> str | None:
    match = SHORT_TILE_RE.match(str(value))
    return match.group(1) if match else None


def parse_source_column(row: pd.Series) -> tuple[dict | None, dict | None]:
    source_column = clean_source_label(row.get("source_column"))
    short_id = short_sample_id(str(row.get("sample_id")))
    base = {
        "sample_id": short_id or row.get("sample_id"),
        "matrix_sample_id": row.get("sample_id"),
        "source_column": source_column,
        "dataset_batch": "GSE60012",
    }
    if not short_id:
        return None, {**base, "exclude_reason": "missing_tile_sample_id"}
    if not source_column:
        return None, {**base, "exclude_reason": "empty_source_column"}
    parts = source_column.split("_")
    if len(parts) < 4:
        return None, {**base, "exclude_reason": "not_enough_header_tokens"}

    tissue_raw, sex_raw, age_raw = parts[0], parts[1], parts[2]
    tissue = TISSUE_MAP.get(tissue_raw.lower())
    if tissue is None:
        return None, {**base, "exclude_reason": f"unsupported_tissue:{tissue_raw}"}
    if sex_raw not in {"F", "M"}:
        return None, {**base, "exclude_reason": f"unsupported_sex:{sex_raw}"}
    age_match = AGE_RE.match(age_raw)
    if not age_match:
        return None, {**base, "exclude_reason": f"unparseable_age:{age_raw}"}

    age_weeks = float(age_match.group(1))
    condition_detail = "_".join(parts[3:]).strip("_")
    if not condition_detail:
        return None, {**base, "exclude_reason": "missing_condition_detail"}
    condition_lower = condition_detail.lower()
    intervention = "control" if condition_lower == "normal" else "other_intervention"
    condition_family = "normal" if condition_lower == "normal" else condition_lower

    record = {
        "sample_id": short_id,
        "title": f"GSE60012 synthetic {source_column}",
        "age_days": age_weeks * 7,
        "age_weeks": age_weeks,
        "raw_age_value": age_weeks,
        "raw_age_unit": "weeks",
        "normalized_age_unit": "weeks",
        "age_unit_source": "matrix_header_synthetic",
        "age_parse_warning": "",
        "tissue": tissue,
        "strain": "C57/bl",
        "sex": sex_raw,
        "intervention": intervention,
        "dataset_batch": "GSE60012",
        "assay": "RRBS_tile_matrix",
        "fastq_on_disk": False,
        "parse_status": "ok_header_synthetic_non_gsm",
        "raw_source_name": tissue_raw,
        "raw_characteristics": (
            f"matrix_sample_id: {row.get('sample_id')} | source_column: {source_column} | "
            f"condition_detail: {condition_detail}"
        ),
        "sample_id_source": "matrix_header_synthetic",
        "metadata_source": "matrix_header_synthetic_non_gsm",
        "matrix_sample_id": row.get("sample_id"),
        "source_column": source_column,
        "condition_detail": condition_detail,
        "condition_family": condition_family,
    }
    return record, None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base_metadata", default=str(BASE_METADATA))
    parser.add_argument("--parse_stats", default=str(GSE60012_STATS))
    parser.add_argument("--out_metadata", default=str(OUT_METADATA))
    parser.add_argument("--out_model_metadata", default=str(OUT_MODEL_METADATA))
    parser.add_argument("--out_exclusions", default=str(OUT_EXCLUSIONS))
    args = parser.parse_args()

    stats = pd.read_csv(args.parse_stats)
    records: list[dict] = []
    exclusions: list[dict] = []
    for _, row in stats.iterrows():
        record, exclusion = parse_source_column(row)
        if record is not None:
            records.append(record)
        if exclusion is not None:
            exclusions.append(exclusion)

    gse_meta = pd.DataFrame(records)
    if gse_meta.empty:
        raise SystemExit("No GSE60012 header-derived samples were parseable.")
    if gse_meta["sample_id"].duplicated().any():
        dupes = sorted(gse_meta.loc[gse_meta["sample_id"].duplicated(), "sample_id"].unique())
        raise SystemExit(f"Duplicate synthetic sample_id values: {dupes[:10]}")

    out_metadata = Path(args.out_metadata)
    out_model_metadata = Path(args.out_model_metadata)
    out_exclusions = Path(args.out_exclusions)
    out_metadata.parent.mkdir(parents=True, exist_ok=True)
    out_model_metadata.parent.mkdir(parents=True, exist_ok=True)
    out_exclusions.parent.mkdir(parents=True, exist_ok=True)
    gse_meta.to_csv(out_metadata, index=False)
    pd.DataFrame(exclusions).to_csv(out_exclusions, index=False)

    base = pd.read_csv(args.base_metadata)
    base = base[base["dataset_batch"] != "GSE60012"].copy()
    for column in gse_meta.columns:
        if column not in base.columns:
            base[column] = ""
    for column in base.columns:
        if column not in gse_meta.columns:
            gse_meta[column] = ""
    model_meta = pd.concat([base, gse_meta[base.columns]], ignore_index=True)
    if model_meta["sample_id"].duplicated().any():
        dupes = sorted(model_meta.loc[model_meta["sample_id"].duplicated(), "sample_id"].unique())
        raise SystemExit(f"Duplicate model sample_id values: {dupes[:10]}")
    model_meta.to_csv(out_model_metadata, index=False)

    print(
        "[GSE60012 metadata] "
        f"parsed={len(gse_meta)} excluded={len(exclusions)} "
        f"model_metadata={len(model_meta)}"
    )
    print(f"  wrote {out_metadata}")
    print(f"  wrote {out_model_metadata}")
    print(f"  wrote {out_exclusions}")


if __name__ == "__main__":
    main()

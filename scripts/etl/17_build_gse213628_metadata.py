#!/usr/bin/env python3
"""Build v8 model metadata rows for GSE213628.

GSE213628 is introduced as an old-age multi-tissue chronological-age support
dataset. The builder uses sample-specific fields only and intentionally does
not infer sample ages from protocol-level age lists.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")
DEFAULT_BASE = ROOT / "metadata" / "model_sample_metadata_v7_1.csv"
DEFAULT_CANDIDATES = ROOT / "metadata" / "geo_old_age_candidate_samples.csv"
DEFAULT_GSE_OUT = ROOT / "metadata" / "gse213628_sample_metadata.csv"
DEFAULT_MODEL_OUT = ROOT / "metadata" / "model_sample_metadata_v8.csv"

MONTH_TO_WEEK = 30.42 / 7.0
AGE_TITLE_RE = re.compile(r"(?P<value>\d+(?:\.\d+)?)\s*month", re.IGNORECASE)


def parse_age_from_title(title: str) -> tuple[float | None, float | None, str, str]:
    match = AGE_TITLE_RE.search(str(title))
    if not match:
        return None, None, "", "missing_sample_title_age"
    months = float(match.group("value"))
    weeks = round(months * MONTH_TO_WEEK, 3)
    return months, weeks, match.group(0), ""


def normalize_tissue(source_name: str, tissue_guess: str) -> str:
    text = f"{source_name} {tissue_guess}".lower()
    if "blood" in text:
        return "blood"
    if "heart" in text:
        return "heart"
    if "kidney" in text:
        return "kidney"
    if "liver" in text:
        return "liver"
    if "lung" in text:
        return "lung"
    if "spleen" in text:
        return "spleen"
    if "muscle" in text:
        return "skeletal_muscle"
    if "intestine" in text or "colon" in text:
        return "intestine"
    return str(tissue_guess or "unknown").split(";")[0] or "unknown"


def parse_characteristic(characteristics: str, key: str, default: str = "unknown") -> str:
    pattern = re.compile(rf"(?:^|\|)\s*{re.escape(key)}\s*:\s*([^|]+)", re.IGNORECASE)
    match = pattern.search(str(characteristics))
    if not match:
        return default
    return match.group(1).strip()


def build_gse_rows(candidate_path: Path) -> pd.DataFrame:
    candidates = pd.read_csv(candidate_path)
    gse = candidates[candidates["dataset"].eq("GSE213628")].copy()
    if gse.empty:
        raise SystemExit(f"No GSE213628 rows found in {candidate_path}")

    rows = []
    for _, row in gse.sort_values("sample_id").iterrows():
        raw_age_value, age_weeks, raw_age_token, warning = parse_age_from_title(row.get("title", ""))
        age_days = None if age_weeks is None else round(age_weeks * 7, 3)
        characteristics = str(row.get("characteristics", ""))
        source_name = str(row.get("source_name", ""))
        tissue = normalize_tissue(source_name, str(row.get("tissue_guess", "")))
        rows.append(
            {
                "sample_id": row["sample_id"],
                "title": row.get("title", ""),
                "age_days": age_days,
                "age_weeks": age_weeks,
                "raw_age_value": raw_age_value,
                "raw_age_unit": "months" if raw_age_value is not None else "",
                "normalized_age_unit": "months" if raw_age_value is not None else "",
                "age_unit_source": "sample_title",
                "age_parse_warning": warning,
                "tissue": tissue,
                "strain": parse_characteristic(characteristics, "strain"),
                "sex": "unknown",
                "intervention": "control",
                "dataset_batch": "GSE213628",
                "assay": "RRBS",
                "fastq_on_disk": False,
                "parse_status": "ok" if age_weeks is not None and tissue != "unknown" else "needs_review",
                "raw_source_name": source_name,
                "raw_characteristics": characteristics,
                "sample_id_source": "GEO_SOFT_GSM",
                "metadata_source": "official_geo_soft_v8",
                "matrix_sample_id": row["sample_id"],
                "source_column": np.nan,
                "condition_detail": "none",
                "condition_family": "control",
            }
        )
    return pd.DataFrame(rows)


def validate(gse_meta: pd.DataFrame) -> dict:
    duplicate_count = int(gse_meta["sample_id"].duplicated().sum())
    age_known = int(gse_meta["age_weeks"].notna().sum())
    age_coverage = age_known / max(len(gse_meta), 1)
    tissue_count = int(gse_meta["tissue"].nunique())
    max_age = float(gse_meta["age_weeks"].max())
    status = "passed"
    failures = []
    if len(gse_meta) < 110:
        failures.append("sample_count_lt_110")
    if age_coverage < 0.95:
        failures.append("age_coverage_lt_95pct")
    if not (120 <= max_age <= 123):
        failures.append("max_age_not_around_121_7w")
    if tissue_count < 6:
        failures.append("tissue_count_lt_6")
    if duplicate_count:
        failures.append("duplicate_sample_id")
    if failures:
        status = "failed"
    return {
        "status": status,
        "failures": failures,
        "n_samples": int(len(gse_meta)),
        "age_known": age_known,
        "age_coverage": round(float(age_coverage), 4),
        "max_age_weeks": round(max_age, 3),
        "tissue_count": tissue_count,
        "tissue_counts": gse_meta["tissue"].value_counts().to_dict(),
        "duplicate_sample_ids": duplicate_count,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base_metadata", default=str(DEFAULT_BASE))
    parser.add_argument("--candidate_samples", default=str(DEFAULT_CANDIDATES))
    parser.add_argument("--gse_out", default=str(DEFAULT_GSE_OUT))
    parser.add_argument("--model_out", default=str(DEFAULT_MODEL_OUT))
    args = parser.parse_args()

    base = pd.read_csv(args.base_metadata)
    gse_meta = build_gse_rows(Path(args.candidate_samples))
    validation = validate(gse_meta)
    if validation["status"] != "passed":
        raise SystemExit(f"GSE213628 metadata validation failed: {validation}")

    all_columns = list(dict.fromkeys([*base.columns, *gse_meta.columns]))
    model_meta = pd.concat([base.reindex(columns=all_columns), gse_meta.reindex(columns=all_columns)], ignore_index=True)
    if model_meta["sample_id"].duplicated().any():
        dupes = model_meta.loc[model_meta["sample_id"].duplicated(), "sample_id"].tolist()
        raise SystemExit(f"Duplicate sample_id after v8 metadata merge: {dupes[:10]}")

    Path(args.gse_out).parent.mkdir(parents=True, exist_ok=True)
    gse_meta.to_csv(args.gse_out, index=False)
    model_meta.to_csv(args.model_out, index=False)
    print(f"[GSE213628 metadata] wrote {args.gse_out} rows={len(gse_meta)}")
    print(f"[v8 metadata] wrote {args.model_out} rows={len(model_meta)}")
    print(validation)


if __name__ == "__main__":
    main()

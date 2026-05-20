#!/usr/bin/env python3
"""Audit GSE121141 metadata/source/batch coupling against prediction errors."""
from __future__ import annotations

import argparse
import gzip
import json
import re
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("/home/zdq-as/mouse_methyl_work")
GSM_RE = re.compile(r"(GSM\d+)")
TITLE_RE = re.compile(r"^([A-Za-z]+)(\d{2})(\d+)$")
FLOWCELL_RE = re.compile(r"^H[A-Z0-9]{5,}(?:XX|XY)$")
LANE_RE = re.compile(r"^L\d{1,2}$")
BATCH_RE = re.compile(r"^(DO|KD)\d+")


def parse_characteristics(values: list[str]) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for value in values:
        if ":" not in value:
            continue
        key, raw = value.split(":", 1)
        parsed[key.strip().lower().replace(" ", "_")] = raw.strip()
    return parsed


def parse_soft(path: Path, dataset: str) -> pd.DataFrame:
    rows = []
    current: dict[str, object] | None = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current is not None:
                    rows.append(current)
                current = {
                    "sample_id": line.split("=", 1)[1].strip(),
                    "characteristics": [],
                    "relations": [],
                }
                continue
            if current is None:
                continue
            if line.startswith("^"):
                rows.append(current)
                current = None
                continue
            if " = " not in line:
                continue
            key, value = line.split(" = ", 1)
            key = key.lstrip("!")
            if key == "Sample_series_id" and value != dataset:
                current["not_target_dataset"] = True
            elif key == "Sample_characteristics_ch1":
                current.setdefault("characteristics", []).append(value)
            elif key == "Sample_relation":
                current.setdefault("relations", []).append(value)
            elif key in {
                "Sample_title",
                "Sample_source_name_ch1",
                "Sample_supplementary_file_1",
                "Sample_instrument_model",
                "Sample_library_selection",
                "Sample_library_source",
                "Sample_library_strategy",
                "Sample_extract_protocol_ch1",
                "Sample_data_processing",
            }:
                current[key.removeprefix("Sample_").lower()] = value
        if current is not None:
            rows.append(current)

    normalized = []
    for row in rows:
        if row.get("not_target_dataset"):
            continue
        fields = parse_characteristics(row.get("characteristics", []))
        relations = row.get("relations", [])
        relation_text = " | ".join(relations)
        sra = None
        biosample = None
        for item in relations:
            if "sra?term=" in item:
                sra = item.rsplit("=", 1)[-1]
            if "biosample/" in item:
                biosample = item.rsplit("/", 1)[-1]
        normalized.append(
            {
                "sample_id": row.get("sample_id"),
                "soft_title": row.get("title"),
                "soft_source_name": row.get("source_name_ch1"),
                "soft_tissue": fields.get("tissue"),
                "soft_gender": fields.get("gender"),
                "soft_age_months": fields.get("age,months"),
                "soft_strain": fields.get("strain"),
                "soft_genotype": fields.get("genotype"),
                "soft_sra": sra,
                "soft_biosample": biosample,
                "soft_relations": relation_text,
                "supplementary_file": row.get("supplementary_file_1"),
                "instrument_model": row.get("instrument_model"),
                "library_selection": row.get("library_selection"),
                "library_source": row.get("library_source"),
                "library_strategy": row.get("library_strategy"),
            }
        )
    return pd.DataFrame(normalized)


def parse_title(title: object) -> dict[str, object]:
    text = str(title or "")
    match = TITLE_RE.match(text)
    if not match:
        return {"title_prefix": None, "title_age_code": None, "title_animal_code": None}
    prefix, age_code, animal_code = match.groups()
    return {
        "title_prefix": prefix,
        "title_age_code": int(age_code),
        "title_animal_code": animal_code,
    }


def parse_file_tokens(value: object) -> dict[str, object]:
    text = str(value or "")
    basename = Path(text).name
    match = GSM_RE.search(basename)
    sample_id = match.group(1) if match else None
    core = basename
    for suffix in [".bismark.cov.gz", ".cov.txt.gz", ".cov.gz", ".txt.gz"]:
        if core.endswith(suffix):
            core = core[: -len(suffix)]
            break
    tokens = core.split("_")
    title_token = tokens[1] if len(tokens) > 1 else None
    batch_tokens = [tok for tok in tokens if BATCH_RE.match(tok)]
    flowcells = [tok for tok in tokens if FLOWCELL_RE.match(tok)]
    lanes = [tok for tok in tokens if LANE_RE.match(tok)]
    return {
        "file_basename": basename,
        "file_sample_id": sample_id,
        "file_title_token": title_token,
        "file_batch_token": batch_tokens[0] if batch_tokens else None,
        "file_batch_family": batch_tokens[0][:2] if batch_tokens else None,
        "file_flowcell": flowcells[0] if flowcells else None,
        "file_lane": lanes[0] if lanes else None,
        "file_combined": "combined" in basename.lower(),
        "file_val_trimmed": "_val_" in basename,
        "file_token_count": len(tokens),
    }


def tar_member_table(path: Path) -> pd.DataFrame:
    rows = []
    with tarfile.open(path, "r") as tar:
        for member in tar.getmembers():
            if not member.isfile():
                continue
            parsed = parse_file_tokens(member.name)
            parsed["tar_member"] = member.name
            parsed["tar_size_bytes"] = int(member.size)
            rows.append(parsed)
    return pd.DataFrame(rows)


def age_bins(age_weeks: pd.Series) -> pd.Series:
    return pd.cut(
        age_weeks,
        bins=[0, 4, 13, 26, 52, 104, np.inf],
        labels=["0-4w", "4-13w", "13-26w", "26-52w", "52-104w", "104w+"],
        include_lowest=True,
    ).astype(str)


def add_prediction_columns(base: pd.DataFrame, predictions: Path, label: str) -> pd.DataFrame:
    if not predictions.exists():
        return base
    pred = pd.read_csv(predictions)
    pred = pred[pred["dataset_batch"].eq("GSE121141")].copy()
    pred[f"{label}_age_pred_weeks"] = pred["age_weeks_pred"].astype(float)
    pred[f"{label}_residual_weeks"] = pred["age_weeks_true"].astype(float) - pred["age_weeks_pred"].astype(float)
    pred[f"{label}_abs_error_weeks"] = pred[f"{label}_residual_weeks"].abs()
    keep = [
        "sample_id",
        f"{label}_age_pred_weeks",
        f"{label}_residual_weeks",
        f"{label}_abs_error_weeks",
    ]
    return base.merge(pred[keep], on="sample_id", how="left")


def group_error_summary(df: pd.DataFrame, variables: list[str], labels: list[str]) -> pd.DataFrame:
    rows = []
    for variable in variables:
        if variable not in df.columns:
            continue
        for label in labels:
            error_col = f"{label}_abs_error_weeks"
            residual_col = f"{label}_residual_weeks"
            if error_col not in df.columns:
                continue
            grouped = df.groupby(variable, dropna=False)
            for value, sub in grouped:
                valid = sub[sub[error_col].notna()]
                if valid.empty:
                    continue
                rows.append(
                    {
                        "grouping_variable": variable,
                        "grouping_value": value,
                        "prediction_run": label,
                        "n_samples": int(len(valid)),
                        "age_min": float(valid["age_weeks"].min()),
                        "age_median": float(valid["age_weeks"].median()),
                        "age_max": float(valid["age_weeks"].max()),
                        "tissues": ";".join(sorted(valid["tissue"].dropna().astype(str).unique())),
                        "mean_abs_error_weeks": round(float(valid[error_col].mean()), 3),
                        "median_abs_error_weeks": round(float(valid[error_col].median()), 3),
                        "mean_residual_weeks": round(float(valid[residual_col].mean()), 3),
                    }
                )
    return pd.DataFrame(rows)


def contingency_tables(df: pd.DataFrame, out_dir: Path) -> pd.DataFrame:
    pairs = [
        ("age_months", "tissue"),
        ("age_months", "file_batch_family"),
        ("age_months", "file_batch_token"),
        ("age_months", "file_flowcell"),
        ("age_months", "file_lane"),
        ("tissue", "file_flowcell"),
        ("tissue", "file_batch_family"),
        ("age_bin", "tissue"),
    ]
    summary_rows = []
    for left, right in pairs:
        if left not in df.columns or right not in df.columns:
            continue
        table = pd.crosstab(df[left].fillna("NA"), df[right].fillna("NA"))
        table.to_csv(out_dir / f"contingency_{left}_by_{right}.csv")
        total = int(table.to_numpy().sum())
        row_purity = float((table.max(axis=1).sum() / total)) if total else np.nan
        col_purity = float((table.max(axis=0).sum() / total)) if total else np.nan
        summary_rows.append(
            {
                "left": left,
                "right": right,
                "n": total,
                "row_weighted_purity": round(row_purity, 4) if np.isfinite(row_purity) else None,
                "col_weighted_purity": round(col_purity, 4) if np.isfinite(col_purity) else None,
                "n_left_levels": int(table.shape[0]),
                "n_right_levels": int(table.shape[1]),
            }
        )
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(out_dir / "contingency_coupling_summary.csv", index=False)
    return summary


def correlation_summary(df: pd.DataFrame, labels: list[str]) -> pd.DataFrame:
    numeric_cols = [
        "rows_total",
        "rows_parseable_primary_autosomes",
        "rows_pass_coverage",
        "coverage_pass_fraction",
        "tar_size_bytes",
        "mean_region_total_coverage",
        "median_region_total_coverage",
        "mean_region_cpg_count",
        "region_presence_fraction",
        "age_weeks",
    ]
    rows = []
    for label in labels:
        error_col = f"{label}_abs_error_weeks"
        if error_col not in df.columns:
            continue
        for col in numeric_cols:
            if col not in df.columns:
                continue
            sub = df[[col, error_col]].dropna()
            if len(sub) < 3 or sub[col].std() == 0 or sub[error_col].std() == 0:
                continue
            r, pval = stats.pearsonr(sub[col].astype(float), sub[error_col].astype(float))
            rows.append(
                {
                    "prediction_run": label,
                    "variable": col,
                    "pearson_r_abs_error": round(float(r), 4),
                    "pearson_p": float(pval),
                    "n": int(len(sub)),
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata_path", default=str(ROOT / "metadata" / "model_sample_metadata_v7_1.csv"))
    parser.add_argument("--soft", default=str(ROOT / "metadata" / "geo_downloads" / "GSE121141_family.soft.gz"))
    parser.add_argument("--tar", default=str(ROOT / "raw_downloads" / "geo_supplements" / "GSE121141" / "GSE121141_RAW.tar"))
    parser.add_argument("--parse_stats", default=str(ROOT / "results" / "multidataset" / "GSE121141_sample_parse_stats.csv"))
    parser.add_argument(
        "--coverage_qc",
        default=str(
            ROOT
            / "results"
            / "benchmark_v7_6_coverage_weighted"
            / "coverage_qc"
            / "GSE121141_sample_coverage_count_qc.csv"
        ),
    )
    parser.add_argument(
        "--v74_predictions",
        default=str(
            ROOT
            / "results"
            / "benchmark_v7_4_stability"
            / "05_lgbm_robust_p095_top1000_shift015"
            / "predictions.csv"
        ),
    )
    parser.add_argument(
        "--v75_predictions",
        default=str(
            ROOT
            / "results"
            / "validation_v7_5_gse121141_all_except"
            / "05_lgbm_quantile_p08_top1000_agebin08"
            / "predictions.csv"
        ),
    )
    parser.add_argument(
        "--v76_predictions",
        default=str(
            ROOT
            / "results"
            / "validation_v7_6_gse121141_all_except"
            / "02_weighted_v75_quantile_p08_top1000_agebin08"
            / "predictions.csv"
        ),
    )
    parser.add_argument("--out_dir", default=str(ROOT / "results" / "benchmark_v7_7_gse121141_metadata_batch_audit"))
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(args.metadata_path)
    meta = meta[meta["dataset_batch"].eq("GSE121141") & meta["age_days"].notna()].copy()
    meta["age_months"] = (meta["age_weeks"].astype(float) / 4.34524).round(1)
    meta["age_bin"] = age_bins(meta["age_weeks"])

    soft = parse_soft(Path(args.soft), "GSE121141")
    tar_members = tar_member_table(Path(args.tar))
    parse_stats = pd.read_csv(args.parse_stats)
    coverage_qc = pd.read_csv(args.coverage_qc) if Path(args.coverage_qc).exists() else pd.DataFrame()

    df = meta.merge(soft, on="sample_id", how="left")
    file_tokens = df["supplementary_file"].map(parse_file_tokens).apply(pd.Series)
    df = pd.concat([df, file_tokens], axis=1)
    df = df.merge(
        tar_members,
        left_on="sample_id",
        right_on="file_sample_id",
        how="left",
        suffixes=("", "_tar"),
    )
    # If SOFT supplementary parsing already provided file tokens, prefer those and use tar fields for size/member only.
    for col in ["tar_member", "tar_size_bytes"]:
        if f"{col}_tar" in df.columns and col not in df.columns:
            df[col] = df[f"{col}_tar"]
    if "tar_size_bytes_tar" in df.columns:
        df["tar_size_bytes"] = df.get("tar_size_bytes", df["tar_size_bytes_tar"]).fillna(df["tar_size_bytes_tar"])
    if "tar_member_tar" in df.columns:
        df["tar_member"] = df.get("tar_member", df["tar_member_tar"]).fillna(df["tar_member_tar"])

    title_parts = df["title"].map(parse_title).apply(pd.Series)
    df = pd.concat([df, title_parts], axis=1)
    df = df.merge(parse_stats, on="sample_id", how="left")
    if "rows_pass_coverage" in df.columns and "rows_parseable_primary_autosomes" in df.columns:
        df["coverage_pass_fraction"] = df["rows_pass_coverage"] / df["rows_parseable_primary_autosomes"]
    if not coverage_qc.empty:
        keep_cols = [
            "sample_id",
            "mean_region_total_coverage",
            "median_region_total_coverage",
            "mean_region_cpg_count",
            "median_region_cpg_count",
            "region_presence_fraction",
        ]
        df = df.merge(coverage_qc[[col for col in keep_cols if col in coverage_qc.columns]], on="sample_id", how="left")

    labels = ["v74_groupkfold", "v75_heldout", "v76_weighted_heldout"]
    df = add_prediction_columns(df, Path(args.v74_predictions), labels[0])
    df = add_prediction_columns(df, Path(args.v75_predictions), labels[1])
    df = add_prediction_columns(df, Path(args.v76_predictions), labels[2])

    priority_cols = [
        "sample_id",
        "title",
        "age_weeks",
        "age_months",
        "age_bin",
        "tissue",
        "sex",
        "soft_sra",
        "soft_biosample",
        "file_batch_family",
        "file_batch_token",
        "file_flowcell",
        "file_lane",
        "file_combined",
        "tar_size_bytes",
        "rows_pass_coverage",
        "coverage_pass_fraction",
        "mean_region_total_coverage",
        "mean_region_cpg_count",
        "v74_groupkfold_abs_error_weeks",
        "v75_heldout_abs_error_weeks",
        "v76_weighted_heldout_abs_error_weeks",
        "v75_heldout_residual_weeks",
        "v76_weighted_heldout_residual_weeks",
    ]
    remaining = [col for col in df.columns if col not in priority_cols]
    df[[col for col in priority_cols if col in df.columns] + remaining].to_csv(
        out_dir / "gse121141_sample_source_error_audit.csv", index=False
    )

    old = df[df["age_weeks"] >= 104].copy()
    old.to_csv(out_dir / "gse121141_old_age_source_error_table.csv", index=False)

    group_vars = [
        "age_months",
        "age_bin",
        "tissue",
        "title_prefix",
        "file_batch_family",
        "file_batch_token",
        "file_flowcell",
        "file_lane",
        "file_combined",
    ]
    group_error_summary(df, group_vars, labels).to_csv(out_dir / "error_by_metadata_batch_group.csv", index=False)
    contingency = contingency_tables(df, out_dir)
    corr = correlation_summary(df, labels)
    corr.to_csv(out_dir / "numeric_qc_error_correlations.csv", index=False)

    summary = {
        "dataset": "GSE121141",
        "n_samples": int(len(df)),
        "n_old_age_ge_104w": int((df["age_weeks"] >= 104).sum()),
        "age_levels": sorted(df["age_months"].dropna().unique().tolist()),
        "tissues": sorted(df["tissue"].dropna().unique().tolist()),
        "batch_families": sorted(df["file_batch_family"].dropna().unique().tolist()),
        "flowcells": sorted(df["file_flowcell"].dropna().unique().tolist()),
        "lanes": sorted(df["file_lane"].dropna().unique().tolist()),
        "contingency_summary_rows": int(len(contingency)),
        "numeric_correlation_rows": int(len(corr)),
        "v75_old_age_mean_abs_error": (
            None
            if "v75_heldout_abs_error_weeks" not in old.columns
            else round(float(old["v75_heldout_abs_error_weeks"].mean()), 3)
        ),
        "v75_non_old_mean_abs_error": (
            None
            if "v75_heldout_abs_error_weeks" not in df.columns
            else round(float(df[df["age_weeks"] < 104]["v75_heldout_abs_error_weeks"].mean()), 3)
        ),
        "diagnostic_only": True,
    }
    (out_dir / "gse121141_metadata_batch_audit_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(f"[v7.7 GSE121141 audit] Wrote outputs -> {out_dir}")


if __name__ == "__main__":
    main()

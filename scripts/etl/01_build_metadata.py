#!/usr/bin/env python3
"""
Phase 1: Unified Sample Metadata Builder

Builds a normalized sample metadata table from local GEO SOFT files via GEOparse.
The parser is intentionally conservative and keeps raw fields for auditability.
"""
from __future__ import annotations

import re
import warnings
from pathlib import Path

import GEOparse
import pandas as pd

try:
    import yaml
except ImportError:  # pragma: no cover - depends on local environment
    yaml = None

warnings.filterwarnings("ignore", category=UserWarning)

ROOT = Path("/home/zdq-as/mouse_methyl_work")
CONFIG_DIR = ROOT / "configs"
META_OUT = ROOT / "metadata"
GEO_DIR = META_OUT / "geo_downloads"
META_OUT.mkdir(exist_ok=True)
GEO_DIR.mkdir(exist_ok=True)

DATASET_PROFILES = {
    "GSE120137": {"assay": "RRBS"},
    "GSE80672": {"assay": "RRBS"},
    "GSE93957": {"assay": "RRBS"},
    "GSE121141": {"assay": "RRBS"},
    "GSE60012": {"assay": "RRBS"},
    "GSE52266": {"assay": "RRBS"},
    "GSE80761": {"assay": "RRBS"},
    "GSE45361": {"assay": "RRBS"},
}

WORD_NUMBERS = {
    "zero": 0,
    "one": 1,
    "two": 2,
    "three": 3,
    "four": 4,
    "five": 5,
    "six": 6,
    "seven": 7,
    "eight": 8,
    "nine": 9,
    "ten": 10,
    "eleven": 11,
    "twelve": 12,
    "thirteen": 13,
    "fourteen": 14,
    "fifteen": 15,
    "sixteen": 16,
    "seventeen": 17,
    "eighteen": 18,
    "nineteen": 19,
    "twenty": 20,
    "thirty": 30,
}

NUMBER_PAT = (
    r"\d+(?:\.\d+)?|zero|one|two|three|four|five|six|seven|eight|nine|ten|"
    r"eleven|twelve|thirteen|fourteen|fifteen|sixteen|seventeen|eighteen|"
    r"nineteen|twenty|thirty"
)
UNIT_PAT = r"d|day|days|w|wk|wks|week|weeks|m|mo|mos|month|months|y|yr|yrs|year|years"


def load_dataset_overrides() -> dict:
    path = CONFIG_DIR / "metadata_overrides.yaml"
    if not path.exists():
        return {}
    if yaml is None:
        datasets: dict[str, dict[str, str]] = {}
        current_dataset = None
        for line in path.read_text(encoding="utf-8").splitlines():
            stripped = line.strip()
            if not stripped or stripped.startswith("#") or stripped == "datasets:":
                continue
            if line.startswith("  ") and not line.startswith("    ") and stripped.endswith(":"):
                current_dataset = stripped[:-1]
                datasets[current_dataset] = {}
                continue
            if current_dataset and line.startswith("    ") and ":" in stripped and not stripped.startswith("-"):
                key, value = stripped.split(":", 1)
                datasets[current_dataset][key] = value.strip().strip("\"'")
        return datasets
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    return data.get("datasets", {})


def normalize_text(value: object) -> str:
    return re.sub(r"\s+", " ", str(value).replace("\xa0", " ")).strip()


def parse_characteristics(chars: list[str]) -> tuple[dict[str, str], str]:
    fields: dict[str, str] = {}
    raw_parts: list[str] = []
    for char in chars:
        clean = normalize_text(char)
        raw_parts.append(clean)
        if ":" not in clean:
            continue
        key, value = clean.split(":", 1)
        key_norm = re.sub(r"[^a-z0-9]+", "_", key.lower()).strip("_")
        fields[key_norm] = normalize_text(value)
    return fields, " | ".join(raw_parts)


def parse_number(token: str) -> float | None:
    token = token.lower().strip()
    if re.fullmatch(r"\d+(?:\.\d+)?", token):
        return float(token)
    if token in WORD_NUMBERS:
        return float(WORD_NUMBERS[token])
    if "-" in token:
        total = 0.0
        for part in token.split("-"):
            if part not in WORD_NUMBERS:
                return None
            total += WORD_NUMBERS[part]
        return total
    return None


def unit_multiplier(unit: str | None) -> float | None:
    if not unit:
        return None
    unit = unit.lower()
    if unit in {"d", "day", "days"}:
        return 1.0
    if unit in {"w", "wk", "wks", "week", "weeks"}:
        return 7.0
    if unit in {"m", "mo", "mos", "month", "months"}:
        return 30.42
    if unit in {"y", "yr", "yrs", "year", "years"}:
        return 365.0
    return None


def canonical_unit(unit: str | None) -> str | None:
    if unit is None:
        return None
    unit = unit.lower()
    if unit in {"d", "day", "days"}:
        return "days"
    if unit in {"w", "wk", "wks", "week", "weeks"}:
        return "weeks"
    if unit in {"m", "mo", "mos", "month", "months"}:
        return "months"
    if unit in {"y", "yr", "yrs", "year", "years"}:
        return "years"
    return unit


def unit_from_key(key: str) -> str | None:
    key = key.lower()
    if any(x in key for x in ["year", "yr"]):
        return "years"
    if any(x in key for x in ["month", "mo"]):
        return "months"
    if any(x in key for x in ["week", "wk"]):
        return "weeks"
    if "day" in key:
        return "days"
    return None


def age_parse_result(
    value: float,
    raw_unit: str | None,
    default_unit: str | None,
    override_unit: str | None,
    override_source: str | None,
) -> dict | None:
    raw_unit = canonical_unit(raw_unit or default_unit)
    normalized_unit = canonical_unit(override_unit or raw_unit)
    mult = unit_multiplier(normalized_unit)
    if mult is None:
        return None
    warning = None
    age_unit_source = "parsed"
    if override_unit:
        age_unit_source = override_source or "dataset_override"
        if raw_unit and raw_unit != normalized_unit:
            warning = f"age_unit_override:{raw_unit}->{normalized_unit}"
    return {
        "age_days": value * mult,
        "raw_age_value": value,
        "raw_age_unit": raw_unit,
        "normalized_age_unit": normalized_unit,
        "age_unit_source": age_unit_source,
        "age_parse_warning": warning,
    }


def parse_age_expression_detail(
    text: str,
    default_unit: str | None = None,
    override_unit: str | None = None,
    override_source: str | None = None,
) -> dict | None:
    text = normalize_text(text).lower()
    explicit = re.search(rf"\b({NUMBER_PAT})\s*({UNIT_PAT})\b", text)
    if explicit:
        val = parse_number(explicit.group(1))
        if val is not None:
            return age_parse_result(val, explicit.group(2), default_unit, override_unit, override_source)

    compact = re.search(rf"\b({NUMBER_PAT})(d|w|m|y)\b", text)
    if compact:
        val = parse_number(compact.group(1))
        if val is not None:
            return age_parse_result(val, compact.group(2), default_unit, override_unit, override_source)

    labelled = re.search(rf"\bage[^:;=,]*[:=,]\s*({NUMBER_PAT})\s*({UNIT_PAT})?\b", text)
    if labelled:
        val = parse_number(labelled.group(1))
        if val is not None:
            return age_parse_result(val, labelled.group(2), default_unit, override_unit, override_source)

    val = parse_number(text)
    if val is not None:
        return age_parse_result(val, None, default_unit, override_unit, override_source)
    return None


def parse_age_expression(text: str, default_unit: str | None = None) -> float | None:
    result = parse_age_expression_detail(text, default_unit)
    return result["age_days"] if result else None


def parse_age_detail(
    fields: dict[str, str],
    fallback_text: str,
    dataset_override: dict | None = None,
) -> dict:
    dataset_override = dataset_override or {}
    override_unit = dataset_override.get("age_unit_override")
    override_source = dataset_override.get("age_unit_source")

    for key, value in fields.items():
        if "age" not in key:
            continue
        result = parse_age_expression_detail(value, unit_from_key(key), override_unit, override_source)
        if result is not None:
            return result
        result = parse_age_expression_detail(
            f"{key}: {value}", unit_from_key(key), override_unit, override_source
        )
        if result is not None:
            return result

    # Fallback is deliberately limited to strings that mention age to avoid
    # extracting unrelated ages from growth protocols.
    for chunk in fallback_text.split("|"):
        if "age" not in chunk.lower():
            continue
        result = parse_age_expression_detail(chunk, None, override_unit, override_source)
        if result is not None:
            return result
    return {
        "age_days": None,
        "raw_age_value": None,
        "raw_age_unit": None,
        "normalized_age_unit": None,
        "age_unit_source": None,
        "age_parse_warning": None,
    }


def parse_age_days(fields: dict[str, str], fallback_text: str) -> float | None:
    result = parse_age_detail(fields, fallback_text)
    return result["age_days"]


def normalize_tissue(fields: dict[str, str], title: str, source: str) -> str:
    text = " ".join(
        [
            fields.get("tissue", ""),
            fields.get("organ", ""),
            source,
            title,
        ]
    ).lower()
    tissue_map = [
        ("whole blood", "blood"),
        ("blood", "blood"),
        ("cerebellum", "brain_cerebellum"),
        ("hippocampus", "brain_hippo"),
        ("cortex", "brain_cortex"),
        ("brain", "brain"),
        ("fibroblast", "fibroblast"),
        ("adipose", "adipose"),
        ("liver", "liver"),
        ("kidney", "kidney"),
        ("lung", "lung"),
        ("heart", "heart"),
        ("muscle", "muscle"),
        ("spleen", "spleen"),
        ("pancreas", "pancreas"),
        ("adrenal", "adrenal"),
        ("endometrium", "endometrium"),
    ]
    for needle, tissue in tissue_map:
        if needle in text:
            return tissue
    return "unknown"


def normalize_sex(fields: dict[str, str]) -> str:
    value = normalize_text(fields.get("sex", fields.get("gender", ""))).lower()
    if value in {"m", "male"}:
        return "M"
    if value in {"f", "female"}:
        return "F"
    return "unknown"


def normalize_strain(fields: dict[str, str]) -> str:
    value = normalize_text(fields.get("strain", ""))
    return value if value else "unknown"


def normalize_intervention(fields: dict[str, str], title: str, source: str, treatment: str) -> str:
    text = " | ".join(
        [
            title,
            source,
            treatment,
            fields.get("diet", ""),
            fields.get("adult_diet", ""),
            fields.get("maternal_diet", ""),
            fields.get("condition", ""),
            fields.get("genetic_condition", ""),
            fields.get("genotype", ""),
        ]
    ).lower()

    if "rapamycin" in text:
        return "rapamycin"
    if (
        "calorie restricted" in text
        or "caloric restricted" in text
        or "diet restriction" in text
        or "dietary restriction" in text
        or re.search(r"\bcr\b", text)
    ):
        return "CR"
    if "castrat" in text:
        return "castration"
    if "high fat" in text or "high-fat" in text or "obesogenic" in text:
        return "diet_high_fat"
    if any(x in text for x in ["standard", "low fat", "low-fat", "wild type", "normal", "control"]):
        return "control"
    return "control"


def build_dataset(gse: str, dataset_overrides: dict) -> list[dict]:
    print(f"\n[Meta] Parsing {gse} via GEOparse...")
    try:
        gse_data = GEOparse.get_GEO(geo=gse, destdir=str(GEO_DIR), silent=True)
    except Exception as exc:
        print(f"  [Error] Failed to get {gse}: {exc}")
        return []

    records: list[dict] = []
    for gsm_name, gsm in sorted(gse_data.gsms.items()):
        title = normalize_text(gsm.metadata.get("title", [""])[0])
        source = normalize_text(gsm.metadata.get("source_name_ch1", [""])[0])
        treatment = normalize_text(gsm.metadata.get("treatment_protocol_ch1", [""])[0])
        chars = gsm.metadata.get("characteristics_ch1", [])
        fields, raw_characteristics = parse_characteristics(chars)

        fallback_text = " | ".join([title, source, treatment, raw_characteristics])
        age_detail = parse_age_detail(fields, fallback_text, dataset_overrides.get(gse, {}))
        age_days = age_detail["age_days"]
        intervention = normalize_intervention(fields, title, source, treatment)
        parse_status = []
        if age_days is None:
            parse_status.append("missing_age")
        if normalize_tissue(fields, title, source) == "unknown":
            parse_status.append("missing_tissue")

        records.append(
            {
                "sample_id": gsm_name,
                "title": title,
                "age_days": round(float(age_days), 3) if age_days is not None else None,
                "age_weeks": round(float(age_days) / 7, 3) if age_days is not None else None,
                "raw_age_value": age_detail["raw_age_value"],
                "raw_age_unit": age_detail["raw_age_unit"],
                "normalized_age_unit": age_detail["normalized_age_unit"],
                "age_unit_source": age_detail["age_unit_source"],
                "age_parse_warning": age_detail["age_parse_warning"],
                "tissue": normalize_tissue(fields, title, source),
                "strain": normalize_strain(fields),
                "sex": normalize_sex(fields),
                "intervention": intervention,
                "dataset_batch": gse,
                "assay": DATASET_PROFILES[gse]["assay"],
                "fastq_on_disk": True,
                "parse_status": "ok" if not parse_status else ";".join(parse_status),
                "raw_source_name": source,
                "raw_characteristics": raw_characteristics,
            }
        )

    age_known = sum(1 for record in records if record["age_days"] is not None)
    interventions = sorted({record["intervention"] for record in records})
    print(f"  Found {len(records)} samples. {age_known} have age info. Interventions: {interventions}")
    return records


def write_completeness(df: pd.DataFrame) -> None:
    rows = []
    for dataset, group in df.groupby("dataset_batch"):
        rows.append(
            {
                "dataset": dataset,
                "n_samples": len(group),
                "n_age_known": int(group["age_days"].notna().sum()),
                "age_known_pct": round(float(group["age_days"].notna().mean() * 100), 1),
                "n_fastq_on_disk": int(group["fastq_on_disk"].sum()),
                "tissues": ",".join(sorted(group["tissue"].dropna().unique())),
                "intervention": ",".join(sorted(group["intervention"].dropna().unique())),
            }
        )
    out = META_OUT / "dataset_completeness.csv"
    pd.DataFrame(rows).sort_values("dataset").to_csv(out, index=False)
    print(f"[Meta] Completeness summary saved -> {out}")


def main() -> None:
    print("=" * 60)
    print("Phase 1: Unified GEO Metadata Builder")
    print("=" * 60)

    dataset_overrides = load_dataset_overrides()
    all_records: list[dict] = []
    for gse in DATASET_PROFILES:
        all_records.extend(build_dataset(gse, dataset_overrides))

    df = pd.DataFrame(all_records)
    out = META_OUT / "unified_sample_metadata.csv"
    df.to_csv(out, index=False)
    write_completeness(df)

    print(f"\n[Meta] Total records: {len(df)}")
    print(f"[Meta] Age known: {df['age_days'].notna().sum()}")
    print(f"[Meta] Non-control interventions: {(df['intervention'] != 'control').sum()}")
    print(f"[Meta] Saved -> {out}")


if __name__ == "__main__":
    main()

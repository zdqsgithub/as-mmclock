#!/usr/bin/env python3
"""Build a unified RRBS 5kb region matrix across converted datasets."""
from __future__ import annotations

import argparse
import json
import re
import time
from pathlib import Path

import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")
MULTI_DIR = ROOT / "results" / "multidataset"
META_FILE = ROOT / "metadata" / "unified_sample_metadata.csv"
GSM_RE = re.compile(r"(GSM\d+)")
GSE60012_TILE_RE = re.compile(r"^(GSE60012_tile_\d{3})")


def resolve_matrix_sample_id(value: object) -> str:
    text = str(value)
    tile = GSE60012_TILE_RE.match(text)
    if tile:
        return tile.group(1)
    gsm = GSM_RE.search(text)
    if gsm:
        return gsm.group(1)
    return text


def existing_path(candidates: list[Path]) -> Path:
    for path in candidates:
        if path.exists():
            return path
    return candidates[0]


def dataset_matrix_path(dataset: str, input_dirs: list[Path]) -> Path:
    if dataset == "GSE120137":
        return ROOT / "results" / "phase0" / "region_matrix_5kb.parquet"
    return existing_path([input_dir / f"{dataset}_region_matrix_5kb.parquet" for input_dir in input_dirs])


def dataset_manifest_path(dataset: str, input_dirs: list[Path]) -> Path:
    if dataset == "GSE120137":
        return ROOT / "results" / "phase0" / "region_stats_5kb.csv"
    path = existing_path([input_dir / f"{dataset}_matrix_manifest.json" for input_dir in input_dirs])
    if dataset == "GSE80672" and not path.exists():
        return existing_path([input_dir / "matrix_manifest.json" for input_dir in input_dirs])
    return path


def load_manifest(dataset: str, matrix: pd.DataFrame, meta: pd.DataFrame, input_dirs: list[Path]) -> dict:
    if dataset == "GSE120137":
        columns = [resolve_matrix_sample_id(col) for col in matrix.columns]
        overlap = sorted(set(columns) & set(meta.loc[meta["dataset_batch"] == dataset, "sample_id"]))
        return {
            "dataset": dataset,
            "status": "completed",
            "source": "phase0_region_matrix",
            "n_samples_parsed": int(matrix.shape[1]),
            "n_samples_metadata_overlap": int(len(overlap)),
            "metadata_age_known_overlap": int(
                meta[(meta["dataset_batch"] == dataset) & (meta["sample_id"].isin(overlap))]["age_days"].notna().sum()
            ),
            "n_regions": int(matrix.shape[0]),
        }
    path = dataset_manifest_path(dataset, input_dirs)
    if not path.exists():
        return {"dataset": dataset, "status": "missing_manifest", "n_samples_metadata_overlap": 0}
    manifest = json.loads(path.read_text(encoding="utf-8"))
    manifest.setdefault("dataset", dataset)
    return manifest


def load_dataset_matrix(dataset: str, meta: pd.DataFrame, input_dirs: list[Path]) -> tuple[pd.DataFrame | None, dict]:
    path = dataset_matrix_path(dataset, input_dirs)
    if not path.exists():
        return None, {"dataset": dataset, "status": "missing_matrix", "matrix_path": str(path)}
    matrix = pd.read_parquet(path)
    matrix.columns = [resolve_matrix_sample_id(col) for col in matrix.columns]
    manifest = load_manifest(dataset, matrix, meta, input_dirs)
    allowed_statuses = {"completed"}
    if dataset == "GSE60012" and manifest.get("status") == "completed_unmapped":
        allowed_statuses.add("completed_unmapped")
        manifest["metadata_source"] = "matrix_header_synthetic_non_gsm"
    if manifest.get("status") not in allowed_statuses:
        manifest["skip_reason"] = f"manifest_status_{manifest.get('status')}"
        return None, manifest
    gse_meta = meta[(meta["dataset_batch"] == dataset) & meta["age_days"].notna()].copy()
    common = [sample for sample in gse_meta["sample_id"] if sample in matrix.columns]
    if not common:
        manifest["skip_reason"] = "no_metadata_age_overlap"
        return None, manifest
    matrix = matrix[common]
    manifest.update(
        {
            "matrix_path": str(path),
            "n_samples_included": int(matrix.shape[1]),
            "n_regions_loaded": int(matrix.shape[0]),
        }
    )
    return matrix, manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--datasets", default="GSE120137,GSE80672,GSE93957,GSE121141,GSE60012")
    parser.add_argument("--join", default="inner", choices=["inner", "outer"])
    parser.add_argument("--out_dir", default=str(MULTI_DIR))
    parser.add_argument(
        "--input_dirs",
        default=str(MULTI_DIR),
        help="Comma-separated directories containing per-dataset region matrices; first existing path wins.",
    )
    parser.add_argument("--metadata_path", default=str(META_FILE))
    parser.add_argument("--min_regions", type=int, default=0)
    args = parser.parse_args()

    started = time.time()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(args.metadata_path)
    datasets = [item.strip() for item in args.datasets.split(",") if item.strip()]
    input_dirs = [Path(item.strip()) for item in args.input_dirs.split(",") if item.strip()]

    matrices: list[pd.DataFrame] = []
    manifests: list[dict] = []
    for dataset in datasets:
        matrix, manifest = load_dataset_matrix(dataset, meta, input_dirs)
        manifests.append(manifest)
        if matrix is not None:
            matrices.append(matrix)
            print(f"  Include {dataset}: {matrix.shape[0]} regions x {matrix.shape[1]} samples")
        else:
            print(f"  Skip {dataset}: {manifest.get('skip_reason') or manifest.get('status')}")

    if len(matrices) < 2:
        raise SystemExit("Need at least two metadata-mapped matrices to build a multidataset benchmark matrix.")

    region_sets = [set(matrix.index) for matrix in matrices]
    if args.join == "inner":
        regions = sorted(set.intersection(*region_sets))
    else:
        regions = sorted(set.union(*region_sets))
    if len(regions) < args.min_regions:
        blocker = {
            "status": "blocked",
            "reason": "common_region_count_below_threshold",
            "join": args.join,
            "min_regions": int(args.min_regions),
            "n_regions": int(len(regions)),
            "datasets_requested": datasets,
            "dataset_manifests": manifests,
            "metadata_path": str(Path(args.metadata_path)),
        }
        (out_dir / "all_rrbs_matrix_blocker.json").write_text(json.dumps(blocker, indent=2), encoding="utf-8")
        raise SystemExit(
            f"Common region count {len(regions)} is below --min_regions={args.min_regions}; wrote blocker."
        )
    aligned = [matrix.reindex(regions) for matrix in matrices]
    combined = pd.concat(aligned, axis=1).astype("float32")
    out_path = out_dir / "all_rrbs_region_matrix_5kb.parquet"
    combined.to_parquet(out_path, compression="zstd")

    included_datasets = [m.get("dataset", "unknown") for m in manifests if m.get("n_samples_included")]
    common_meta = meta[meta["sample_id"].isin([resolve_matrix_sample_id(col) for col in combined.columns])]
    manifest = {
        "status": "completed",
        "join": args.join,
        "metadata_path": str(Path(args.metadata_path)),
        "input_dirs": [str(path) for path in input_dirs],
        "sample_id_resolver": "gsm_or_gse60012_tile_short_id",
        "datasets_requested": datasets,
        "datasets_included": included_datasets,
        "dataset_manifests": manifests,
        "n_regions": int(combined.shape[0]),
        "n_samples": int(combined.shape[1]),
        "dataset_sample_counts": common_meta["dataset_batch"].value_counts().to_dict(),
        "intervention_counts": common_meta["intervention"].value_counts().to_dict(),
        "age_known_samples": int(common_meta["age_days"].notna().sum()),
        "matrix_path": str(out_path),
        "exec_time_sec": round(time.time() - started, 1),
    }
    (out_dir / "all_rrbs_matrix_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[Multidataset] Wrote {out_path} shape={combined.shape}")
    print(f"[Multidataset] Included: {', '.join(included_datasets)}")


if __name__ == "__main__":
    main()

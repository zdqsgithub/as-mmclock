#!/usr/bin/env python3
"""Build v7.6 mixed coverage-weighted multidataset region matrix."""
from __future__ import annotations

import argparse
import json
import re
import time
from pathlib import Path

import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")
WEIGHTED_DIR = ROOT / "results" / "multidataset_v7_6_weighted"
MULTI_DIR = ROOT / "results" / "multidataset"
META_FILE = ROOT / "metadata" / "model_sample_metadata_v7_1.csv"
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


def weighted_matrix_path(dataset: str, weighted_dir: Path) -> Path:
    return weighted_dir / f"{dataset}_coverage_weighted_region_matrix_5kb.parquet"


def fallback_matrix_path(dataset: str) -> Path:
    if dataset == "GSE120137":
        return ROOT / "results" / "phase0" / "region_matrix_5kb.parquet"
    return MULTI_DIR / f"{dataset}_region_matrix_5kb.parquet"


def manifest_path(dataset: str, weighted_dir: Path) -> Path:
    return weighted_dir / f"{dataset}_coverage_weighted_manifest.json"


def load_dataset(dataset: str, meta: pd.DataFrame, weighted_dir: Path, allow_fallback: bool) -> tuple[pd.DataFrame | None, dict]:
    weighted_path = weighted_matrix_path(dataset, weighted_dir)
    source = "coverage_weighted"
    path = weighted_path
    if not path.exists():
        if not allow_fallback:
            return None, {"dataset": dataset, "status": "missing_weighted_matrix", "path": str(path)}
        source = "unweighted_fallback"
        path = fallback_matrix_path(dataset)
    if not path.exists():
        return None, {"dataset": dataset, "status": "missing_matrix", "path": str(path), "source": source}
    matrix = pd.read_parquet(path)
    matrix.columns = [resolve_matrix_sample_id(col) for col in matrix.columns]
    gse_meta = meta[(meta["dataset_batch"] == dataset) & meta["age_days"].notna()].copy()
    common = [sample for sample in gse_meta["sample_id"] if sample in matrix.columns]
    if not common:
        return None, {"dataset": dataset, "status": "no_metadata_overlap", "path": str(path), "source": source}
    manifest = {
        "dataset": dataset,
        "status": "completed",
        "source": source,
        "path": str(path),
        "n_regions_loaded": int(matrix.shape[0]),
        "n_samples_loaded": int(matrix.shape[1]),
        "n_samples_included": int(len(common)),
    }
    if source == "coverage_weighted" and manifest_path(dataset, weighted_dir).exists():
        manifest["weighted_manifest"] = json.loads(manifest_path(dataset, weighted_dir).read_text(encoding="utf-8"))
    return matrix[common], manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--datasets", default="GSE120137,GSE80672,GSE93957,GSE121141")
    parser.add_argument("--weighted_dir", default=str(WEIGHTED_DIR))
    parser.add_argument("--metadata_path", default=str(META_FILE))
    parser.add_argument("--out_dir", default=str(WEIGHTED_DIR))
    parser.add_argument("--join", default="inner", choices=["inner", "outer"])
    parser.add_argument("--min_regions", type=int, default=50000)
    parser.add_argument("--allow_fallback", action="store_true")
    args = parser.parse_args()

    started = time.time()
    weighted_dir = Path(args.weighted_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(args.metadata_path)
    datasets = [item.strip() for item in args.datasets.split(",") if item.strip()]

    matrices = []
    manifests = []
    for dataset in datasets:
        matrix, manifest = load_dataset(dataset, meta, weighted_dir, args.allow_fallback)
        manifests.append(manifest)
        if matrix is None:
            print(f"  Skip {dataset}: {manifest.get('status')}")
            continue
        matrices.append(matrix)
        print(f"  Include {dataset}: {matrix.shape[0]} regions x {matrix.shape[1]} samples ({manifest['source']})")
    if len(matrices) < 2:
        raise SystemExit("Need at least two datasets to build v7.6 matrix.")

    region_sets = [set(matrix.index) for matrix in matrices]
    regions = sorted(set.intersection(*region_sets) if args.join == "inner" else set.union(*region_sets))
    if len(regions) < args.min_regions:
        blocker = {
            "status": "blocked",
            "reason": "region_count_below_threshold",
            "join": args.join,
            "n_regions": int(len(regions)),
            "min_regions": int(args.min_regions),
            "datasets": datasets,
            "dataset_manifests": manifests,
            "metadata_path": str(Path(args.metadata_path)),
        }
        (out_dir / "all_rrbs_weighted_matrix_blocker.json").write_text(
            json.dumps(blocker, indent=2), encoding="utf-8"
        )
        raise SystemExit(f"Region count {len(regions)} below min_regions={args.min_regions}; wrote blocker.")

    combined = pd.concat([matrix.reindex(regions) for matrix in matrices], axis=1).astype("float32")
    out_path = out_dir / "all_rrbs_region_matrix_5kb.parquet"
    combined.to_parquet(out_path, compression="zstd")
    common_meta = meta[meta["sample_id"].isin([resolve_matrix_sample_id(col) for col in combined.columns])]
    manifest = {
        "status": "completed",
        "matrix_family": "v7_6_coverage_weighted_where_available",
        "join": args.join,
        "metadata_path": str(Path(args.metadata_path)),
        "sample_id_resolver": "gsm_or_gse60012_tile_short_id",
        "allow_fallback": bool(args.allow_fallback),
        "datasets_requested": datasets,
        "dataset_manifests": manifests,
        "n_regions": int(combined.shape[0]),
        "n_samples": int(combined.shape[1]),
        "dataset_sample_counts": common_meta["dataset_batch"].value_counts().to_dict(),
        "intervention_counts": common_meta["intervention"].value_counts().to_dict(),
        "matrix_path": str(out_path),
        "exec_time_sec": round(time.time() - started, 1),
    }
    (out_dir / "all_rrbs_matrix_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[v7.6 matrix] Wrote {out_path} shape={combined.shape}")


if __name__ == "__main__":
    main()


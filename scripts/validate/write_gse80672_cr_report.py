#!/usr/bin/env python3
"""Write the v5 GSE80672 held-out CR validation report."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")
DEFAULT_MATRIX_DIR = ROOT / "results" / "multidataset"
DEFAULT_VALIDATION_DIR = ROOT / "results" / "validation_gse80672_cr"
DEFAULT_DOC = ROOT / "doc" / "20_analysis" / "03_20260518_multidataset_gse80672_cr_validation_report.md"


def load_json(path: Path) -> dict | None:
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def fmt(value, digits: int = 3) -> str:
    if value is None:
        return "N/A"
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def build_blocker_report(blocker: dict, inventory: pd.DataFrame | None) -> str:
    inventory_lines = []
    if inventory is not None and not inventory.empty:
        row = inventory.iloc[0].to_dict()
        inventory_lines = [
            "| Dataset | Available | Remote size | Local status | Local path |",
            "|---|---:|---:|---|---|",
            f"| {row.get('dataset')} | {row.get('available')} | {row.get('remote_size_bytes')} | {row.get('local_status')} | `{row.get('local_path')}` |",
            "",
        ]
    return "\n".join(
        [
            "# Multidataset GSE80672 CR Validation v5 Report",
            "",
            "> Date: 2026-05-18  ",
            "> Scope: GSE80672 processed supplement inventory, conversion, held-out CR validation",
            "",
            "## Status",
            "",
            "Blocked before training. A GSE80672 matrix was not created, so held-out age prediction and CR statistics were not computed.",
            "",
            "## Inventory",
            "",
            *inventory_lines,
            "## Blocker",
            "",
            f"- Reason: `{blocker.get('reason')}`",
            f"- Input: `{blocker.get('input_path')}`",
            f"- Errors recorded: {blocker.get('n_errors', 0)}",
            "",
            "## Decision",
            "",
            "Do not run biological-age validation without a real GSE80672 matrix. The next step is to inspect the processed supplement schema and, if it cannot be converted, move to a minimal FASTQ-to-Bismark ETL plan for GSE80672 only.",
            "",
        ]
    )


def build_success_report(
    manifest: dict,
    heldout: dict,
    benchmark: dict,
    shuffled: dict | None,
    predictions_path: Path,
) -> str:
    cr_counts = manifest.get("metadata_intervention_counts", {})
    shuffled_auc = None if shuffled is None else shuffled.get("cr_detection_auc")
    shuffled_f1 = None if shuffled is None else shuffled.get("cr_detection_f1")
    cr_auc = benchmark.get("cr_detection_auc")
    cr_f1 = benchmark.get("cr_detection_f1")
    cr_d = benchmark.get("cr_cohens_d")
    cr_p = benchmark.get("cr_mannwhitney_p")

    return "\n".join(
        [
            "# Multidataset GSE80672 CR Validation v5 Report",
            "",
            "> Date: 2026-05-18  ",
            "> Scope: GSE80672 processed supplement conversion, GSE120137 -> GSE80672 held-out age prediction, real CR residual validation",
            "",
            "## Decision",
            "",
            "GSE80672 processed methylation data were converted into beta and 5kb region matrices and aligned to project metadata. Held-out validation now uses real GSE80672 predictions rather than metadata-only labels or dummy CR metrics.",
            "",
            "## Matrix Build",
            "",
            "| Artifact | Value |",
            "|---|---:|",
            f"| Parsed samples | {manifest.get('n_samples_parsed')} |",
            f"| Metadata-overlap samples | {manifest.get('n_samples_metadata_overlap')} |",
            f"| Intervention counts | {cr_counts} |",
            f"| CpGs after presence filter | {manifest.get('n_cpg_after_presence_filter')} |",
            f"| 5kb regions | {manifest.get('n_regions')} |",
            f"| Conversion seconds | {manifest.get('exec_time_sec')} |",
            "",
            f"- Beta matrix: `{manifest.get('beta_matrix_path')}`",
            f"- Region matrix: `{manifest.get('region_matrix_path')}`",
            f"- Region stats: `{manifest.get('region_stats_path')}`",
            "",
            "## Held-Out Age Benchmark",
            "",
            "| Train -> Test | Model | Common regions | Selected regions | Pearson r | MAE weeks | R2 |",
            "|---|---|---:|---:|---:|---:|---:|",
            f"| {heldout.get('train_dataset')} -> {heldout.get('test_dataset')} | {heldout.get('model_type')} | {heldout.get('n_common_features')} | {heldout.get('n_features_selected')} | {fmt(benchmark.get('pearson_r'), 4)} | {fmt(benchmark.get('mae_weeks'))} | {fmt(benchmark.get('r2'), 4)} |",
            "",
            f"- Predictions: `{predictions_path}`",
            "- Leakage controls: feature selection, imputation, scaling, and model fit use the training dataset only.",
            "",
            "Phase 0 reference points:",
            "",
            "| Benchmark | Pearson r | MAE weeks | R2 |",
            "|---|---:|---:|---:|",
            "| v3 region best, GSE120137 CV | 0.8530 | 11.550 | 0.7051 |",
            "| v4 embedding POC best, GSE120137 CV | 0.8594 | 11.106 | 0.7221 |",
            f"| v5 held-out GSE80672 | {fmt(benchmark.get('pearson_r'), 4)} | {fmt(benchmark.get('mae_weeks'))} | {fmt(benchmark.get('r2'), 4)} |",
            "",
            "## CR Biological-Age Validation",
            "",
            "CR statistics are computed only from real held-out predictions. Positive Cohen's d means controls have higher age acceleration than CR after chronological-age adjustment.",
            "",
            "| Metric | Value |",
            "|---|---:|",
            f"| CR AUC | {fmt(cr_auc, 4)} |",
            f"| CR F1 | {fmt(cr_f1, 4)} |",
            f"| CR Cohen's d | {fmt(cr_d, 4)} |",
            f"| CR Mann-Whitney p | {fmt(cr_p, 6)} |",
            "",
            "## Sanity Check",
            "",
            "| Check | CR AUC | CR F1 | Interpretation |",
            "|---|---:|---:|---|",
            f"| Shuffled intervention labels | {fmt(shuffled_auc, 4)} | {fmt(shuffled_f1, 4)} | should collapse toward chance |",
            "",
            "## Next Step",
            "",
            "If held-out MAE is high or CR statistics are unstable, prioritize adding GSE93957/GSE121141/GSE60012 processed matrices and switch benchmark selection to true `GroupKFold(dataset_batch)`. Do not expand embedding or deep learning until the multidataset baseline stabilizes.",
            "",
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix_dir", default=str(DEFAULT_MATRIX_DIR))
    parser.add_argument("--validation_dir", default=str(DEFAULT_VALIDATION_DIR))
    parser.add_argument("--inventory", default=str(ROOT / "metadata" / "geo_supplement_inventory.csv"))
    parser.add_argument("--doc_out", default=str(DEFAULT_DOC))
    args = parser.parse_args()

    matrix_dir = Path(args.matrix_dir)
    validation_dir = Path(args.validation_dir)
    inventory = pd.read_csv(args.inventory) if Path(args.inventory).exists() else None

    blocker = load_json(matrix_dir / "GSE80672_conversion_blocker.json")
    manifest = load_json(matrix_dir / "matrix_manifest.json")
    heldout = load_json(validation_dir / "result.json")
    benchmark = load_json(validation_dir / "benchmark_result.json")
    shuffled = load_json(validation_dir / "benchmark_shuffled_intervention.json")

    if blocker is not None and (manifest is None or manifest.get("status") != "completed"):
        text = build_blocker_report(blocker, inventory)
    elif manifest is None or heldout is None or benchmark is None:
        text = build_blocker_report(
            {
                "reason": "missing_required_validation_artifacts",
                "input_path": str(matrix_dir),
                "n_errors": 0,
            },
            inventory,
        )
    else:
        text = build_success_report(manifest, heldout, benchmark, shuffled, validation_dir / "predictions.csv")

    validation_dir.mkdir(parents=True, exist_ok=True)
    (validation_dir / "cr_validation_report.md").write_text(text, encoding="utf-8")
    doc_out = Path(args.doc_out)
    doc_out.parent.mkdir(parents=True, exist_ok=True)
    doc_out.write_text(text, encoding="utf-8")
    print(f"[Report] Wrote {validation_dir / 'cr_validation_report.md'}")
    print(f"[Report] Wrote {doc_out}")


if __name__ == "__main__":
    main()

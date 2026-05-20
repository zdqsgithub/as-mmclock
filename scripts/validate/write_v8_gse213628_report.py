#!/usr/bin/env python3
"""Write the v8 GSE213628 RALPH report."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}


def fmt(value, digits: int = 3) -> str:
    if value is None or value != value:
        return "N/A"
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def metric_from_dataset_csv(path: Path, dataset: str, key: str):
    if not path.exists():
        return None
    df = pd.read_csv(path)
    sub = df[df["dataset_batch"].eq(dataset)]
    if sub.empty or key not in sub.columns:
        return None
    return sub[key].iloc[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark_root", default=str(ROOT / "results" / "benchmark_v8_gse213628"))
    parser.add_argument("--matrix_dir", default=str(ROOT / "results" / "multidataset_v8_gse213628"))
    parser.add_argument("--metadata", default=str(ROOT / "metadata" / "gse213628_sample_metadata.csv"))
    parser.add_argument(
        "--out",
        default=str(ROOT / "doc" / "20_analysis" / "15_20260518_v8_gse213628_ralph_report.md"),
    )
    args = parser.parse_args()

    benchmark_root = Path(args.benchmark_root)
    matrix_dir = Path(args.matrix_dir)
    metadata = pd.read_csv(args.metadata) if Path(args.metadata).exists() else pd.DataFrame()
    summary_path = benchmark_root / "v8_gse213628_ralph_summary.tsv"
    summary = pd.read_csv(summary_path, sep="\t") if summary_path.exists() else pd.DataFrame()
    decision = load_json(benchmark_root / "v8_gse213628_ralph_decision.json")
    all_manifest = load_json(matrix_dir / "all_rrbs_matrix_manifest.json")
    gse_manifest = load_json(matrix_dir / "GSE213628_matrix_manifest.json")

    best_name = decision.get("best_config")
    best = summary[summary["experiment"].eq(best_name)].iloc[0].to_dict() if best_name and not summary.empty else {}
    random_benchmark = decision.get("random_label_benchmark", {})
    heldouts = decision.get("heldout_results", {})
    gse121141 = heldouts.get("GSE121141", {})
    gse80672 = heldouts.get("GSE80672", {})
    gse121141_dir = Path(gse121141.get("dir", ""))
    gse80672_dir = Path(gse80672.get("dir", ""))
    gse121141_old_mae = metric_from_dataset_csv(gse121141_dir / "dataset_metrics.csv", "GSE121141", "old_104w_mae_weeks")

    v75_group_mae = 23.767
    v75_old_mae = 75.386
    group_improvement = None if not best else v75_group_mae - float(best.get("mae_weeks"))
    old_improvement = None if gse121141_old_mae is None else v75_old_mae - float(gse121141_old_mae)
    autoresearch_ready = bool(
        group_improvement is not None
        and old_improvement is not None
        and group_improvement >= 2
        and old_improvement >= 10
        and abs(float(random_benchmark.get("pearson_r", 999))) < 0.2
    )

    lines = [
        "# v8 GSE213628 老龄多组织数据闭环报告",
        "",
        "Date: 2026-05-18",
        "",
        "## Summary",
        "",
        (
            "v8 将 `GSE213628` 加入多数据集 5kb region benchmark，用来检验 v7.8 发现的 "
            "old-age same-tissue support blocker 是否被缓解。"
        ),
        "",
        f"- Best config: `{best_name or 'N/A'}`",
        f"- GroupKFold MAE: {fmt(best.get('mae_weeks'))}w; v7.5 baseline 23.767w; improvement {fmt(group_improvement)}w",
        f"- GSE121141 held-out 104w+ MAE: {fmt(gse121141_old_mae)}w; v7.5 diagnostic baseline 75.386w; improvement {fmt(old_improvement)}w",
        f"- Random-label Pearson r: {fmt(random_benchmark.get('pearson_r'), 4)}",
        f"- Constrained autoresearch ready: `{autoresearch_ready}`",
        "",
        "## Data Artifacts",
        "",
        f"- GSE213628 metadata rows: {len(metadata)}",
        f"- GSE213628 parsed samples: {fmt(gse_manifest.get('n_samples_parsed'), 0)}",
        f"- GSE213628 metadata overlap: {fmt(gse_manifest.get('n_samples_metadata_overlap'), 0)}",
        f"- GSE213628 5kb regions: {fmt(gse_manifest.get('n_regions'), 0)}",
        f"- v8 matrix shape: {fmt(all_manifest.get('n_regions'), 0)} regions x {fmt(all_manifest.get('n_samples'), 0)} samples",
        f"- Included datasets: {', '.join(all_manifest.get('datasets_included', []))}",
        "",
        "## Benchmark Summary",
        "",
    ]
    if not summary.empty:
        display_cols = [
            "experiment",
            "pearson_r",
            "mae_weeks",
            "r2",
            "cross_dataset_mae",
            "cr_detection_auc",
            "gse121141_mae_weeks",
            "gse121141_old_104w_mae_weeks",
        ]
        lines.append(summary[[c for c in display_cols if c in summary.columns]].to_markdown(index=False))
    else:
        lines.append("No benchmark summary found.")
    lines.extend(
        [
            "",
            "## Held-Out Validation",
            "",
            "| Test dataset | Pearson r | MAE weeks | R2 | CR AUC | Shuffled CR AUC |",
            "|---|---:|---:|---:|---:|---:|",
            (
                f"| GSE121141 | {fmt(gse121141.get('benchmark', {}).get('pearson_r'), 4)} | "
                f"{fmt(gse121141.get('benchmark', {}).get('mae_weeks'))} | "
                f"{fmt(gse121141.get('benchmark', {}).get('r2'), 4)} | N/A | N/A |"
            ),
            (
                f"| GSE80672 | {fmt(gse80672.get('benchmark', {}).get('pearson_r'), 4)} | "
                f"{fmt(gse80672.get('benchmark', {}).get('mae_weeks'))} | "
                f"{fmt(gse80672.get('benchmark', {}).get('r2'), 4)} | "
                f"{fmt(gse80672.get('benchmark', {}).get('cr_detection_auc'), 4)} | "
                f"{fmt(gse80672.get('shuffled_benchmark', {}).get('cr_detection_auc'), 4)} |"
            ),
            "",
            "## Decision",
            "",
        ]
    )
    if autoresearch_ready:
        lines.append(
            "v8 passed the preconditions for constrained autoresearch. Run 30-50 configs within the locked "
            "`lgbm/ridge/elasticnet` + region-filter search space."
        )
    else:
        lines.append(
            "v8 did not meet all thresholds for constrained autoresearch. Keep the result as a RALPH-loop data "
            "diagnostic and prioritize schema/coverage harmonization or targeted P2 validation before larger search."
        )
    lines.extend(
        [
            "",
            "## Files",
            "",
            f"- Matrix manifest: `{matrix_dir / 'all_rrbs_matrix_manifest.json'}`",
            f"- Benchmark summary: `{summary_path}`",
            f"- Decision JSON: `{benchmark_root / 'v8_gse213628_ralph_decision.json'}`",
        ]
    )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[Report] wrote {out}")


if __name__ == "__main__":
    main()

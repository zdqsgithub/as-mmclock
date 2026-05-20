#!/usr/bin/env python3
"""Run v8.3 matrix ablation and tissue-support diagnostics.

This runner is intentionally narrow: it compares fixed dataset combinations and
the two locked v8 RALPH configs. It does not perform autoresearch.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("/home/zdq-as/mouse_methyl_work")

VARIANTS = {
    "core4": ["GSE120137", "GSE80672", "GSE93957", "GSE121141"],
    "core4_plus_gse213628": ["GSE120137", "GSE80672", "GSE93957", "GSE121141", "GSE213628"],
    "core4_plus_gse60012": ["GSE120137", "GSE80672", "GSE93957", "GSE121141", "GSE60012"],
    "all6": ["GSE120137", "GSE80672", "GSE93957", "GSE121141", "GSE60012", "GSE213628"],
}

CONFIGS = [
    {
        "name": "01_v75_lgbm_quantile_p08_top1000_agebin08",
        "model": "lgbm",
        "preprocess": "quantile_uniform",
        "presence": "0.8",
        "group_presence": "0.8",
        "max_shift": "0.15",
        "agebin_presence": "0.8",
        "top_n": "1000",
    },
    {
        "name": "02_v74_lgbm_robust_p095_top1000_shift015",
        "model": "lgbm",
        "preprocess": "robust",
        "presence": "0.95",
        "group_presence": "0.95",
        "max_shift": "0.15",
        "agebin_presence": None,
        "top_n": "1000",
    },
]

V75_BASELINES = {
    "groupkfold_mae_weeks": 23.767,
    "gse121141_old_104w_mae_weeks": 75.386,
}


def run(cmd: list[str], cwd: Path = ROOT) -> None:
    print(" ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=cwd, check=True)


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}


def safe_r(true: np.ndarray, pred: np.ndarray) -> tuple[float, float]:
    if len(true) < 3 or np.std(true) == 0 or np.std(pred) == 0:
        return 0.0, 1.0
    r, pval = stats.pearsonr(true, pred)
    if np.isnan(r):
        return 0.0, 1.0
    return float(r), float(pval)


def metric_row(dataset: str, pred: pd.DataFrame) -> dict:
    true = pred["age_weeks_true"].to_numpy(float)
    yhat = pred["age_weeks_pred"].to_numpy(float)
    r, pval = safe_r(true, yhat)
    ss_res = float(np.sum((true - yhat) ** 2))
    ss_tot = float(np.sum((true - true.mean()) ** 2))
    old = pred[pred["age_weeks_true"] >= 104]
    return {
        "test_dataset": dataset,
        "n_samples": int(len(pred)),
        "pearson_r": round(r, 4),
        "pearson_pval": pval,
        "mae_weeks": round(float(np.mean(np.abs(true - yhat))), 3),
        "rmse_weeks": round(float(np.sqrt(np.mean((true - yhat) ** 2))), 3),
        "r2": round(float(1 - ss_res / ss_tot), 4) if ss_tot > 0 else 0.0,
        "old_104w_n": int(len(old)),
        "old_104w_mae_weeks": (
            None
            if old.empty
            else round(float(np.mean(np.abs(old["age_weeks_true"] - old["age_weeks_pred"]))), 3)
        ),
    }


def run_benchmark(predictions: Path, out_path: Path) -> dict:
    if not out_path.exists():
        run([sys.executable, "scripts/validate/benchmark_metrics.py", "--predictions", str(predictions), "--out", str(out_path)])
    return load_json(out_path)


def train_args(config: dict, matrix: Path, metadata_path: Path, out_dir: Path, randomize: bool = False) -> list[str]:
    cmd = [
        sys.executable,
        "scripts/train/train_clock.py",
        "--matrix_path",
        str(matrix),
        "--metadata_path",
        str(metadata_path),
        "--feature_type",
        "region",
        "--model_type",
        config["model"],
        "--n_feature_prefilter",
        config["top_n"],
        "--imputation",
        "median",
        "--preprocess",
        config["preprocess"],
        "--min_train_feature_presence",
        config["presence"],
        "--min_train_group_feature_presence",
        config["group_presence"],
        "--max_train_dataset_mean_shift",
        config["max_shift"],
        "--target_transform",
        "log1p_days",
        "--output_dir",
        str(out_dir),
    ]
    if config["agebin_presence"] is not None:
        cmd.extend(["--min_train_agebin_feature_presence", config["agebin_presence"]])
    if randomize:
        cmd.append("--randomize_labels")
    return cmd


def heldout_args(config: dict, matrix: Path, metadata_path: Path, test_dataset: str, out_dir: Path) -> list[str]:
    cmd = [
        sys.executable,
        "scripts/train/train_heldout_clock.py",
        "--train_matrix",
        str(matrix),
        "--test_matrix",
        str(matrix),
        "--train_datasets",
        f"all_except:{test_dataset}",
        "--test_datasets",
        test_dataset,
        "--metadata_path",
        str(metadata_path),
        "--model_type",
        config["model"],
        "--n_feature_prefilter",
        config["top_n"],
        "--imputation",
        "median",
        "--preprocess",
        config["preprocess"],
        "--min_train_feature_presence",
        config["presence"],
        "--min_train_group_feature_presence",
        config["group_presence"],
        "--max_train_dataset_mean_shift",
        config["max_shift"],
        "--target_transform",
        "log1p_days",
        "--output_dir",
        str(out_dir),
    ]
    if config["agebin_presence"] is not None:
        cmd.extend(["--min_train_agebin_feature_presence", config["agebin_presence"]])
    return cmd


def build_variant_matrix(
    variant: str,
    datasets: list[str],
    matrix_root: Path,
    metadata_path: Path,
    input_dirs: str,
    min_regions: int,
) -> tuple[Path | None, dict]:
    out_dir = matrix_root / variant
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = out_dir / "all_rrbs_matrix_manifest.json"
    blocker_path = out_dir / "all_rrbs_matrix_blocker.json"
    if not manifest_path.exists() and not blocker_path.exists():
        cmd = [
            sys.executable,
            "scripts/etl/12_build_multidataset_region_matrix.py",
            "--datasets",
            ",".join(datasets),
            "--metadata_path",
            str(metadata_path),
            "--input_dirs",
            input_dirs,
            "--out_dir",
            str(out_dir),
            "--min_regions",
            str(min_regions),
        ]
        print(" ".join(cmd), flush=True)
        completed = subprocess.run(cmd, cwd=ROOT, check=False)
        if completed.returncode != 0 and not blocker_path.exists():
            completed.check_returncode()
    if blocker_path.exists() and not manifest_path.exists():
        return None, load_json(blocker_path)
    manifest = load_json(manifest_path)
    return out_dir / "all_rrbs_region_matrix_5kb.parquet", manifest


def shuffle_intervention(predictions: Path, out_csv: Path, seed: int = 42) -> None:
    pred = pd.read_csv(predictions)
    rng = np.random.default_rng(seed)
    shuffled = pred["intervention"].to_numpy(str).copy()
    rng.shuffle(shuffled)
    pred["intervention"] = shuffled
    pred.to_csv(out_csv, index=False)


def add_age_bin(df: pd.DataFrame, col: str) -> pd.Series:
    return pd.cut(
        df[col],
        bins=[0, 4, 13, 26, 52, 104, np.inf],
        labels=["0-4w", "4-13w", "13-26w", "26-52w", "52-104w", "104w+"],
        include_lowest=True,
    ).astype(str)


def write_residual_summary(predictions: Path, variant: str, config_name: str) -> pd.DataFrame:
    pred = pd.read_csv(predictions)
    pred["variant"] = variant
    pred["config"] = config_name
    pred["age_bin"] = add_age_bin(pred, "age_weeks_true")
    pred["abs_error_weeks"] = (pred["age_weeks_true"] - pred["age_weeks_pred"]).abs()
    return (
        pred.groupby(["variant", "config", "dataset_batch", "tissue", "age_bin"], dropna=False)
        .agg(
            n_samples=("sample_id", "count"),
            mae_weeks=("abs_error_weeks", "mean"),
            median_abs_error_weeks=("abs_error_weeks", "median"),
            residual_weeks_mean=("residual_weeks", "mean"),
            residual_weeks_median=("residual_weeks", "median"),
        )
        .reset_index()
    )


def tissue_support_audit(
    predictions: Path,
    metadata_path: Path,
    train_datasets: list[str],
    variant: str,
    config_name: str,
) -> pd.DataFrame:
    pred = pd.read_csv(predictions)
    meta = pd.read_csv(metadata_path)
    train = meta[meta["dataset_batch"].isin(train_datasets) & meta["age_weeks"].notna()].copy()
    rows = []
    for _, row in pred.iterrows():
        tissue = row.get("tissue")
        same = train[train["tissue"].astype(str).eq(str(tissue))]
        rows.append(
            {
                "variant": variant,
                "config": config_name,
                "sample_id": row["sample_id"],
                "dataset_batch": row["dataset_batch"],
                "tissue": tissue,
                "age_weeks_true": float(row["age_weeks_true"]),
                "age_weeks_pred": float(row["age_weeks_pred"]),
                "residual_weeks": float(row["residual_weeks"]),
                "abs_error_weeks": abs(float(row["age_weeks_true"]) - float(row["age_weeks_pred"])),
                "is_104w_plus": bool(float(row["age_weeks_true"]) >= 104),
                "train_same_tissue_n": int(len(same)),
                "train_same_tissue_age_min": None if same.empty else float(same["age_weeks"].min()),
                "train_same_tissue_age_max": None if same.empty else float(same["age_weeks"].max()),
                "outside_same_tissue_train_age_range": bool(
                    same.empty
                    or float(row["age_weeks_true"]) < float(same["age_weeks"].min())
                    or float(row["age_weeks_true"]) > float(same["age_weeks"].max())
                ),
                "train_datasets": ",".join(train_datasets),
            }
        )
    return pd.DataFrame(rows)


def run_heldout(
    config: dict,
    matrix: Path,
    metadata_path: Path,
    test_dataset: str,
    out_dir: Path,
    shuffle_cr: bool = False,
) -> tuple[dict, pd.DataFrame]:
    out_dir.mkdir(parents=True, exist_ok=True)
    if not (out_dir / "result.json").exists():
        run(heldout_args(config, matrix, metadata_path, test_dataset, out_dir))
    benchmark = run_benchmark(out_dir / "predictions.csv", out_dir / "benchmark_result.json")
    if shuffle_cr:
        shuffled_csv = out_dir / "predictions_shuffled_intervention.csv"
        shuffled_json = out_dir / "benchmark_shuffled_intervention.json"
        if not shuffled_csv.exists():
            shuffle_intervention(out_dir / "predictions.csv", shuffled_csv)
        benchmark["shuffled_intervention"] = run_benchmark(shuffled_csv, shuffled_json)
    return benchmark, pd.read_csv(out_dir / "predictions.csv")


def config_by_name(name: str) -> dict:
    for config in CONFIGS:
        if config["name"] == name:
            return config
    raise KeyError(name)


def write_report(
    out_path: Path,
    variant_summary: pd.DataFrame,
    lodo_summary: pd.DataFrame,
    matrix_rows: list[dict],
    random_rows: list[dict],
    gse121141_support: pd.DataFrame,
) -> None:
    def md_table(df: pd.DataFrame) -> str:
        if df.empty:
            return "No rows."
        view = df.copy()
        for col in view.columns:
            view[col] = view[col].map(lambda value: "" if pd.isna(value) else str(value))
        header = "| " + " | ".join(view.columns) + " |"
        sep = "| " + " | ".join(["---"] * len(view.columns)) + " |"
        rows = ["| " + " | ".join(row) + " |" for row in view.to_numpy(dtype=str)]
        return "\n".join([header, sep, *rows])

    ready = variant_summary[
        (variant_summary["groupkfold_mae_weeks"] <= V75_BASELINES["groupkfold_mae_weeks"])
        & (variant_summary["gse121141_old_104w_improvement_vs_v8_2"] >= 10)
        & (variant_summary["random_label_pass"] == True)
        & (variant_summary["gse80672_shuffled_cr_pass"] == True)
    ]
    lines = [
        "# v8.3 Matrix Ablation 与 Tissue-Shared Diagnostics 报告",
        "",
        "Date: 2026-05-18",
        "",
        "## Summary",
        "",
        "v8.3 比较了 core4、加入 GSE213628、加入 GSE60012、all6 四种 strict matrix variant，并只运行两个锁定配置。没有执行 autoresearch、deep learning 或新下载。",
        "",
        f"- Autoresearch-ready variants: {len(ready)}",
        f"- Best GroupKFold MAE: {variant_summary['groupkfold_mae_weeks'].min():.3f}w",
        f"- Best GSE121141 104w+ MAE: {variant_summary['gse121141_old_104w_mae_weeks'].min():.3f}w",
        "",
        "## Matrix Variants",
        "",
        md_table(pd.DataFrame(matrix_rows)),
        "",
        "## Variant Summary",
        "",
        md_table(
            variant_summary[
                [
                    "variant",
                    "config",
                    "groupkfold_mae_weeks",
                    "groupkfold_r",
                    "gse121141_mae_weeks",
                    "gse121141_old_104w_mae_weeks",
                    "gse80672_mae_weeks",
                    "gse80672_cr_auc",
                    "gse80672_shuffled_cr_auc",
                    "random_label_r",
                    "random_label_pass",
                    "gse80672_shuffled_cr_pass",
                    "autoresearch_ready",
                ]
            ]
        ),
        "",
        "## LODO Summary",
        "",
        md_table(lodo_summary),
        "",
        "## GSE121141 104w+ Tissue Support",
        "",
    ]
    old = gse121141_support[gse121141_support["is_104w_plus"]].copy()
    if old.empty:
        lines.append("No GSE121141 104w+ rows found.")
    else:
        tissue_summary = (
            old.groupby(["variant", "config", "tissue"], dropna=False)
            .agg(
                n_samples=("sample_id", "count"),
                mae_weeks=("abs_error_weeks", "mean"),
                train_same_tissue_n=("train_same_tissue_n", "median"),
                train_same_tissue_age_max=("train_same_tissue_age_max", "max"),
                outside_range_rate=("outside_same_tissue_train_age_range", "mean"),
            )
            .reset_index()
        )
        lines.append(md_table(tissue_summary))
    lines.extend(
        [
            "",
            "## Decision",
            "",
        ]
    )
    if ready.empty:
        lines.append(
            "No variant met all v8.3 acceptance criteria. Do not start v8.3.1 autoresearch. "
            "Use the ablation outputs to decide whether GSE60012/GSE213628 should remain auxiliary and prioritize "
            "GSE121141-like same-tissue old-age data or finer tissue/schema harmonization."
        )
    else:
        lines.append(
            "At least one variant met all criteria. Use only those variants for v8.3.1 constrained autoresearch."
        )
    lines.extend(
        [
            "",
            "## Files",
            "",
            "- `results/benchmark_v8_3_ablation/variant_summary.tsv`",
            "- `results/benchmark_v8_3_ablation/lodo_summary_by_variant.csv`",
            "- `results/benchmark_v8_3_ablation/gse121141_tissue_age_support.csv`",
            "- `results/benchmark_v8_3_ablation/residual_by_dataset_tissue_agebin.csv`",
        ]
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata_path", default=str(ROOT / "metadata" / "model_sample_metadata_v8.csv"))
    parser.add_argument("--input_dirs", default="results/multidataset_v8_2_prefilter_liftover,results/multidataset")
    parser.add_argument("--matrix_root", default=str(ROOT / "results" / "multidataset_v8_3_ablation"))
    parser.add_argument("--out_root", default=str(ROOT / "results" / "benchmark_v8_3_ablation"))
    parser.add_argument("--min_regions", type=int, default=50_000)
    parser.add_argument(
        "--report",
        default=str(ROOT / "doc" / "20_analysis" / "18_20260518_v8_3_matrix_ablation_tissue_diagnostics_report.md"),
    )
    args = parser.parse_args()

    metadata_path = Path(args.metadata_path)
    matrix_root = Path(args.matrix_root)
    out_root = Path(args.out_root)
    out_root.mkdir(parents=True, exist_ok=True)

    variant_rows = []
    lodo_rows = []
    matrix_rows = []
    random_rows = []
    support_frames = []
    residual_frames = []

    for variant, datasets in VARIANTS.items():
        print(f"[v8.3] variant={variant}", flush=True)
        matrix, manifest = build_variant_matrix(
            variant,
            datasets,
            matrix_root,
            metadata_path,
            args.input_dirs,
            args.min_regions,
        )
        matrix_rows.append(
            {
                "variant": variant,
                "datasets": ",".join(datasets),
                "status": manifest.get("status"),
                "n_regions": manifest.get("n_regions"),
                "n_samples": manifest.get("n_samples"),
                "blocker_reason": manifest.get("reason"),
            }
        )
        if matrix is None:
            continue

        config_rows = []
        for config in CONFIGS:
            config_dir = out_root / variant / config["name"]
            group_dir = config_dir / "groupkfold"
            group_dir.mkdir(parents=True, exist_ok=True)
            if not (group_dir / "result.json").exists():
                run(train_args(config, matrix, metadata_path, group_dir))
            group_benchmark = run_benchmark(group_dir / "predictions.csv", group_dir / "benchmark_result.json")
            residual_frames.append(write_residual_summary(group_dir / "predictions.csv", variant, config["name"]))

            gse121_dir = config_dir / "heldout_gse121141"
            gse121_benchmark, _ = run_heldout(config, matrix, metadata_path, "GSE121141", gse121_dir)
            support_frames.append(
                tissue_support_audit(
                    gse121_dir / "predictions.csv",
                    metadata_path,
                    [dataset for dataset in datasets if dataset != "GSE121141"],
                    variant,
                    config["name"],
                )
            )
            gse121_pred = pd.read_csv(gse121_dir / "predictions.csv")
            old = gse121_pred[gse121_pred["age_weeks_true"] >= 104]
            gse121_old_mae = None if old.empty else round(float(np.mean(np.abs(old["age_weeks_true"] - old["age_weeks_pred"]))), 3)

            gse80672_dir = config_dir / "heldout_gse80672"
            gse80672_benchmark, _ = run_heldout(
                config,
                matrix,
                metadata_path,
                "GSE80672",
                gse80672_dir,
                shuffle_cr=True,
            )
            shuffled = gse80672_benchmark.get("shuffled_intervention", {})

            row = {
                "variant": variant,
                "config": config["name"],
                "datasets": ",".join(datasets),
                "matrix_regions": manifest.get("n_regions"),
                "matrix_samples": manifest.get("n_samples"),
                "groupkfold_r": group_benchmark.get("pearson_r"),
                "groupkfold_mae_weeks": group_benchmark.get("mae_weeks"),
                "groupkfold_r2": group_benchmark.get("r2"),
                "cross_dataset_mae": group_benchmark.get("cross_dataset_mae"),
                "gse121141_r": gse121_benchmark.get("pearson_r"),
                "gse121141_mae_weeks": gse121_benchmark.get("mae_weeks"),
                "gse121141_r2": gse121_benchmark.get("r2"),
                "gse121141_old_104w_mae_weeks": gse121_old_mae,
                "gse121141_old_104w_improvement_vs_v8_2": None if gse121_old_mae is None else round(84.425 - gse121_old_mae, 3),
                "gse80672_r": gse80672_benchmark.get("pearson_r"),
                "gse80672_mae_weeks": gse80672_benchmark.get("mae_weeks"),
                "gse80672_cr_auc": gse80672_benchmark.get("cr_detection_auc"),
                "gse80672_cr_cohens_d": gse80672_benchmark.get("cr_cohens_d"),
                "gse80672_shuffled_cr_auc": shuffled.get("cr_detection_auc"),
                "gse80672_shuffled_cr_cohens_d": shuffled.get("cr_cohens_d"),
                "gse80672_shuffled_cr_pass": (
                    shuffled.get("cr_detection_auc") is not None
                    and abs(float(shuffled.get("cr_detection_auc")) - 0.5) <= 0.1
                    and abs(float(shuffled.get("cr_cohens_d", 999))) < 0.3
                ),
            }
            config_rows.append(row)
            variant_rows.append(row)

        config_df = pd.DataFrame(config_rows).sort_values(
            ["gse121141_old_104w_mae_weeks", "gse121141_mae_weeks", "groupkfold_mae_weeks"],
            ascending=[True, True, True],
        )
        best_name = str(config_df.iloc[0]["config"])
        best_config = config_by_name(best_name)

        random_dir = out_root / variant / best_name / "random_labels"
        random_dir.mkdir(parents=True, exist_ok=True)
        if not (random_dir / "result.json").exists():
            run(train_args(best_config, matrix, metadata_path, random_dir, randomize=True))
        random_benchmark = run_benchmark(random_dir / "predictions.csv", random_dir / "benchmark_result.json")
        random_pass = abs(float(random_benchmark.get("pearson_r", 999))) < 0.2 and float(random_benchmark.get("mae_weeks", 0)) > float(config_df.iloc[0]["groupkfold_mae_weeks"])
        random_rows.append(
            {
                "variant": variant,
                "config": best_name,
                "random_label_r": random_benchmark.get("pearson_r"),
                "random_label_mae_weeks": random_benchmark.get("mae_weeks"),
                "random_label_pass": random_pass,
            }
        )
        variant_rows[-len(config_rows) + int(config_df.index[0] - min(config_df.index))]["selected_best_config"] = True

        for test_dataset in datasets:
            lodo_dir = out_root / variant / best_name / "lodo" / test_dataset
            lodo_benchmark, lodo_pred = run_heldout(best_config, matrix, metadata_path, test_dataset, lodo_dir)
            row = metric_row(test_dataset, lodo_pred)
            row.update(
                {
                    "variant": variant,
                    "config": best_name,
                    "cr_detection_auc": lodo_benchmark.get("cr_detection_auc"),
                }
            )
            lodo_rows.append(row)

    variant_summary = pd.DataFrame(variant_rows)
    random_summary = pd.DataFrame(random_rows)
    if not random_summary.empty:
        variant_summary = variant_summary.merge(random_summary, on=["variant", "config"], how="left")
    variant_summary["selected_best_config"] = variant_summary.get("selected_best_config", False).fillna(False)
    variant_summary["random_label_pass"] = variant_summary["random_label_pass"].fillna(False)
    variant_summary["autoresearch_ready"] = (
        (variant_summary["groupkfold_mae_weeks"] <= V75_BASELINES["groupkfold_mae_weeks"])
        & (variant_summary["gse121141_old_104w_improvement_vs_v8_2"] >= 10)
        & (variant_summary["random_label_pass"])
        & (variant_summary["gse80672_shuffled_cr_pass"])
    )
    variant_summary.to_csv(out_root / "variant_summary.tsv", sep="\t", index=False)

    lodo_summary = pd.DataFrame(lodo_rows)
    lodo_summary.to_csv(out_root / "lodo_summary_by_variant.csv", index=False)

    support = pd.concat(support_frames, ignore_index=True) if support_frames else pd.DataFrame()
    support.to_csv(out_root / "gse121141_tissue_age_support.csv", index=False)

    residual = pd.concat(residual_frames, ignore_index=True) if residual_frames else pd.DataFrame()
    residual.to_csv(out_root / "residual_by_dataset_tissue_agebin.csv", index=False)

    write_report(
        Path(args.report),
        variant_summary,
        lodo_summary,
        matrix_rows,
        random_rows,
        support,
    )
    print(f"[v8.3] wrote {out_root / 'variant_summary.tsv'}", flush=True)
    print(f"[v8.3] wrote {args.report}", flush=True)


if __name__ == "__main__":
    main()

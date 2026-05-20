# v7 Harmonization QC 与分层 Benchmark 报告

Date: 2026-05-18

## Summary

本轮基于 v6 已构建的 4 数据集统一 5kb region matrix，没有继续下载新数据、没有启动 FASTQ/Bismark、没有扩大 embedding/deep learning。实现重点是：

- 在 `train_clock.py` / `train_heldout_clock.py` 中加入 fold-internal `preprocess`、`min_train_feature_presence`、`target_transform`。
- 新增 multidataset QC audit，诊断 dataset/tissue/age-bin 的 missingness、mean beta shift、PCA batch separation 和残差分层。
- 跑完固定 8 个 v7 harmonization smoke 配置，并对最佳配置做 random-label sanity、GSE80672 held-out CR validation、shuffled intervention sanity。

结论：v7 最佳 GroupKFold MAE 为 24.392 周，略差于 v6 基线 24.351 周，未达到改善 >=2 周的进入 v7 autoresearch 条件。当前误差主要来自 dataset/tissue/age-range shift 和 schema/coverage harmonization 不充分，而不是简单 scaler 或 presence 阈值可以解决的问题。

## Inputs And Outputs

Input matrix:

- `results/multidataset/all_rrbs_region_matrix_5kb.parquet`
- 947 samples x 72,079 common 5kb autosomal regions
- Datasets: GSE120137=549, GSE80672=255, GSE121141=81, GSE93957=62

New/updated scripts:

- `scripts/validate/multidataset_qc_audit.py`
- `scripts/train/run_v7_harmonization_smoke.py`
- `scripts/train/train_clock.py`
- `scripts/train/train_heldout_clock.py`

Main outputs:

- `results/benchmark_v7_harmonization/qc/`
- `results/benchmark_v7_harmonization/v7_smoke_summary.tsv`
- `results/benchmark_v7_harmonization/01_lgbm_standard_p05_top500_log1p/lodo_dataset_metrics.csv`
- `results/benchmark_v7_harmonization/random_label_01_lgbm_standard_p05_top500_log1p/`
- `results/validation_v7_gse80672_cr_all_except/01_lgbm_standard_p05_top500_log1p/`

## QC Audit

Dataset feature shift:

| dataset | n | mean missing | mean abs region shift | median abs shift | top shift abs beta |
|---|---:|---:|---:|---:|---:|
| GSE120137 | 549 | 0.0123 | 0.0603 | 0.0567 | 0.3739 |
| GSE121141 | 81 | 0.0070 | 0.0804 | 0.0713 | 0.6274 |
| GSE80672 | 255 | 0.0039 | 0.0876 | 0.0740 | 0.6490 |
| GSE93957 | 62 | 0.0014 | 0.1202 | 0.1180 | 0.7740 |

PCA separation:

| component | variance | dataset eta2 | tissue eta2 | age-bin eta2 |
|---|---:|---:|---:|---:|
| PC1 | 0.4417 | 0.9448 | 0.4921 | 0.2780 |
| PC2 | 0.1151 | 0.3572 | 0.8001 | 0.0761 |
| PC4 | 0.0491 | 0.1997 | 0.8728 | 0.0747 |
| PC5 | 0.0349 | 0.0400 | 0.8214 | 0.0112 |
| PC7 | 0.0140 | 0.0121 | 0.8810 | 0.0023 |
| PC8 | 0.0110 | 0.7081 | 0.1268 | 0.1002 |

Interpretation:

- PC1 is dominated by dataset identity, not age.
- Several later PCs are dominated by tissue identity.
- This is consistent with the poor GroupKFold generalization: the common region matrix still carries strong batch/tissue structure.

Largest residual strata:

| stratum | n | MAE weeks | median abs error |
|---|---:|---:|---:|
| GSE121141 heart 104w+ | 5 | 100.1 | 97.8 |
| GSE121141 brain_cortex 104w+ | 5 | 88.6 | 88.1 |
| GSE121141 lung 104w+ | 5 | 76.6 | 73.8 |
| GSE80672 blood 104w+ | 74 | 48.1 | 49.2 |
| GSE120137 liver 52-104w | 20 | 43.2 | 45.6 |
| GSE120137 lung 52-104w | 20 | 38.8 | 44.3 |

GSE60012 mapping audit:

- Status: blocked.
- Reason: official tile matrix columns are not unique GSM accessions.
- Tile columns: 154.
- Duplicate source columns: 119.
- Metadata overlap: 0.
- Decision: do not force GSE60012 into v7 training; keep it for v7.1 mapping/schema work.

## v7 Smoke Grid

| rank | config | r | MAE w | R2 | cross-dataset MAE | CR AUC |
|---:|---|---:|---:|---:|---:|---:|
| 1 | lgbm standard p=0.5 top500 log1p | 0.5798 | 24.392 | 0.3079 | 26.115 | 0.6547 |
| 2 | lgbm quantile_uniform p=0.8 top1000 log1p | 0.5570 | 24.626 | 0.2730 | 26.576 | 0.7147 |
| 3 | lgbm robust p=0.95 top1000 log1p | 0.5645 | 24.674 | 0.2652 | 26.600 | 0.7064 |
| 4 | lgbm robust p=0.8 top1000 log1p | 0.5553 | 24.697 | 0.2548 | 26.685 | 0.6959 |
| 5 | lgbm standard p=0.8 top1000 log1p | 0.5491 | 24.939 | 0.2328 | 26.828 | 0.6624 |
| 6 | lgbm robust p=0.8 top1000 linear_weeks | 0.4532 | 30.896 | 0.0153 | 31.519 | 0.4920 |
| 7 | ridge robust p=0.8 top500 log1p | 0.3657 | 36.057 | -0.9202 | 31.888 | 0.6468 |
| 8 | elasticnet robust p=0.8 top1000 log1p | 0.2915 | 37.897 | -1.1081 | 34.748 | 0.4557 |

ElasticNetCV was unstable under the initial broad grid and was narrowed for smoke feasibility. Even after narrowing, it emitted convergence warnings and performed poorly; it should not be promoted without a more careful regularization study.

## v6 vs v7 Best

| run | config | r | MAE w | R2 | cross-dataset MAE | CR AUC |
|---|---|---:|---:|---:|---:|---:|
| v6 baseline | lgbm top500 median | 0.5812 | 24.351 | 0.3096 | 26.012 | 0.6561 |
| v7 best | lgbm standard p=0.5 top500 log1p | 0.5798 | 24.392 | 0.3079 | 26.115 | 0.6547 |

Delta v7-v6:

- MAE: +0.041 weeks, worse.
- R2: -0.0017, worse.
- cross-dataset MAE: +0.103 weeks, worse.
- CR AUC: -0.0014, worse.

Therefore v7 smoke does not meet the >=2 week improvement threshold for 30-50 config v7 autoresearch.

## Leave-One-Dataset-Out Metrics

The GroupKFold folds are dataset-held-out folds.

| held-out dataset | n | r | MAE w | RMSE w | R2 |
|---|---:|---:|---:|---:|---:|
| GSE120137 | 549 | 0.4987 | 21.390 | 25.680 | 0.2096 |
| GSE121141 | 81 | 0.3564 | 37.746 | 48.512 | -0.5772 |
| GSE80672 | 255 | 0.7322 | 28.441 | 34.851 | 0.4331 |
| GSE93957 | 62 | 0.8064 | 16.883 | 18.995 | -0.6928 |

GSE121141 is the clearest age-range/tissue mismatch failure. GSE80672 remains useful for CR validation but held-out chronological error is still too high for strong biological-age claims.

## Sanity Checks

Random-label GroupKFold using v7 best configuration:

| check | r | MAE w | R2 | cross-dataset MAE |
|---|---:|---:|---:|---:|
| real labels | 0.5798 | 24.392 | 0.3079 | 26.115 |
| randomized labels | -0.0411 | 35.197 | -0.4661 | 36.009 |

Result: passed. Pearson r collapses and MAE moves toward random-label behavior.

GSE80672 held-out CR validation:

| run | r | MAE w | R2 | CR AUC | CR F1 | Cohen d | Mann-Whitney p |
|---|---:|---:|---:|---:|---:|---:|---:|
| real intervention | 0.7322 | 28.441 | 0.4331 | 0.6616 | 0.3019 | 0.4813 | 0.001569 |
| shuffled intervention | 0.7322 | 28.441 | 0.4331 | 0.5245 | 0.2138 | 0.0492 | 0.327352 |

Result: shuffled intervention removes the CR signal, so the weak CR signal is not caused by the benchmark code alone. However, CR AUC and Cohen d remain below acceptance thresholds.

## Decision

Do not run a larger v7 autoresearch yet. Conservative sklearn preprocessing did not improve cross-dataset age prediction and did not improve GSE80672 CR validation. The strongest diagnosis is:

1. Dataset shift remains severe: PC1 dataset eta2 is 0.9448.
2. Tissue composition remains severe: multiple PCs have tissue eta2 >0.80.
3. Age-range mismatch is visible in residuals, especially GSE121141 old tissues and GSE80672 old blood.
4. Feature presence thresholds and quantile/robust scaling do not solve this mismatch.

Recommended next step:

- Prioritize v7.1 schema/coverage harmonization and sample mapping, especially GSE60012 mapping repair.
- Revisit region definitions and coverage-aware aggregation before expanding model search.
- Keep embedding as interpretation layer only after a stable multi-dataset benchmark improves.
- Do not move to MLP/CNN/Transformer until held-out dataset MAE and CR validation improve under simpler models.

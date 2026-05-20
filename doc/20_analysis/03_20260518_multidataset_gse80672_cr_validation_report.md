# Multidataset GSE80672 CR Validation v5 Report

> Date: 2026-05-18  
> Scope: GSE80672 processed supplement conversion, GSE120137 -> GSE80672 held-out age prediction, real CR residual validation

## Decision

GSE80672 processed methylation data were converted into beta and 5kb region matrices and aligned to project metadata. Held-out validation now uses real GSE80672 predictions rather than metadata-only labels or dummy CR metrics.

## Matrix Build

| Artifact | Value |
|---|---:|
| Parsed samples | 255 |
| Metadata-overlap samples | 255 |
| Intervention counts | {'control': 223, 'CR': 32} |
| CpGs after presence filter | 1940728 |
| 5kb regions | 128739 |
| Conversion seconds | 958.2 |

- Beta matrix: `results/multidataset/GSE80672_beta_matrix.parquet`
- Region matrix: `results/multidataset/GSE80672_region_matrix_5kb.parquet`
- Region stats: `results/multidataset/GSE80672_region_stats_5kb.csv`

## Held-Out Age Benchmark

| Train -> Test | Model | Common regions | Selected regions | Pearson r | MAE weeks | R2 |
|---|---|---:|---:|---:|---:|---:|
| GSE120137 -> GSE80672 | lgbm | 73074 | 500 | 0.7295 | 33.408 | 0.3438 |

- Predictions: `/home/zdq-as/mouse_methyl_work/results/validation_gse80672_cr/predictions.csv`
- Leakage controls: feature selection, imputation, scaling, and model fit use the training dataset only.

Phase 0 reference points:

| Benchmark | Pearson r | MAE weeks | R2 |
|---|---:|---:|---:|
| v3 region best, GSE120137 CV | 0.8530 | 11.550 | 0.7051 |
| v4 embedding POC best, GSE120137 CV | 0.8594 | 11.106 | 0.7221 |
| v5 held-out GSE80672 | 0.7295 | 33.408 | 0.3438 |

## CR Biological-Age Validation

CR statistics are computed only from real held-out predictions. Positive Cohen's d means controls have higher age acceleration than CR after chronological-age adjustment.

| Metric | Value |
|---|---:|
| CR AUC | 0.6717 |
| CR F1 | 0.3014 |
| CR Cohen's d | 0.4177 |
| CR Mann-Whitney p | 0.000849 |

## Sanity Check

| Check | CR AUC | CR F1 | Interpretation |
|---|---:|---:|---|
| Shuffled intervention labels | 0.5074 | 0.1781 | should collapse toward chance |

## Next Step

If held-out MAE is high or CR statistics are unstable, prioritize adding GSE93957/GSE121141/GSE60012 processed matrices and switch benchmark selection to true `GroupKFold(dataset_batch)`. Do not expand embedding or deep learning until the multidataset baseline stabilizes.

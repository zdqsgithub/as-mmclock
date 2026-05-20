# v7.2 GSE60012 Auxiliary Held-Out 报告

Date: 2026-05-18

## Summary

v7.2 没有继续强行构建 5-dataset strict inner-join matrix。主 benchmark 仍保留 v7 的 4-dataset matrix，GSE60012 只作为 auxiliary held-out dataset：训练用 GSE120137/GSE80672/GSE93957/GSE121141，测试用 GSE60012 header-derived synthetic samples。

结果：GSE60012 held-out 最佳 MAE 为 12.216 周，优于 train-median naive baseline 27.043 周；但 random-label 训练 MAE 反而更低，为 7.971 周，且真实模型 R2 为 -2.4807。因此 GSE60012 不能作为模型泛化改善证据，只能保留为 schema/coverage stress test 和辅助数据集。

## Implemented Changes

- Added `scripts/train/run_v7_2_gse60012_aux_heldout.py`.
- Extended `train_heldout_clock.py`:
  - supports `--randomize_train_labels`;
  - preserves optional synthetic metadata columns in predictions, including `condition_family`, `condition_detail`, `metadata_source`, and `matrix_sample_id`.
- Outputs:
  - `results/validation_v7_2_gse60012_aux_heldout/`
  - `results/validation_v7_2_gse80672_cr_all_except/03_lgbm_robust_p095_top1000_log1p/`

## GSE60012 Held-Out Results

Training matrix:

- `results/multidataset/all_rrbs_region_matrix_5kb.parquet`
- Training datasets: all except GSE60012
- Training samples: 947

Test matrix:

- `results/multidataset/GSE60012_region_matrix_5kb.parquet`
- Test samples: 152
- Common train/test regions: 25,374

Fixed auxiliary configs:

| config | r | MAE w | R2 | train presence features |
|---|---:|---:|---:|---:|
| lgbm robust p=0.95 top1000 log1p | 0.2563 | 12.216 | -2.4807 | 24,364 |
| lgbm standard p=0.5 top500 log1p | 0.3055 | 14.520 | -3.8669 | 25,374 |
| lgbm quantile_uniform p=0.8 top1000 log1p | 0.2500 | 15.884 | -4.4542 | 25,374 |

Naive baselines:

| baseline | pred weeks | MAE w | R2 |
|---|---:|---:|---:|
| train median age | 43.457 | 27.043 | -9.6870 |
| test median age oracle | 20.000 | 5.743 | -0.1703 |

Interpretation:

- Best model improves 54.8% vs train-median baseline.
- But it is worse than the test-median oracle and still has negative R2.
- This means GSE60012 age distribution/calibration dominates MAE; MAE alone is not adequate evidence of transfer.

## Sanity Check

Randomized train labels using the best GSE60012 config:

| run | r | MAE w | R2 |
|---|---:|---:|---:|
| real labels | 0.2563 | 12.216 | -2.4807 |
| randomized train labels | 0.0007 | 7.971 | -0.4028 |

Result: sanity check fails if MAE is required to worsen under random labels. Pearson r collapses, but random-label MAE is lower than real-label MAE. This is consistent with a calibration/age-distribution artifact rather than a robust age clock signal.

## Condition Residual Audit

Largest GSE60012 residual strata under best model:

| stratum | n | MAE w | mean residual w |
|---|---:|---:|---:|
| normal muscle F 3w | 5 | 33.122 | -33.122 |
| normal spleen M 1w | 1 | 31.971 | -31.971 |
| normal muscle M 3w | 5 | 29.483 | -29.483 |
| normal spleen M 20w | 4 | 23.534 | -23.534 |
| normal muscle M 20w | 9 | 23.165 | -23.165 |
| young_castrated muscle M 20w | 5 | 20.062 | -20.062 |

Full audit:

- `results/validation_v7_2_gse60012_aux_heldout/gse60012_best_condition_residual_audit.csv`

The model generally over-predicts young GSE60012 tissue ages. GSE60012 condition labels remain exploratory and are not promoted to castration biological-age validation.

## GSE80672 CR Check With GSE60012-Best Config

The best GSE60012 auxiliary config was also tested on GSE80672:

| run | r | MAE w | R2 | CR AUC | CR F1 | Cohen d |
|---|---:|---:|---:|---:|---:|---:|
| v7 best CR baseline | 0.7322 | 28.441 | 0.4331 | 0.6616 | 0.3019 | 0.4813 |
| v7.2 GSE60012-best config | 0.6904 | 30.110 | 0.3670 | 0.7064 | 0.2949 | 0.7512 |
| shuffled intervention | 0.6904 | 30.110 | 0.3670 | 0.5174 | 0.2051 | 0.0686 |

The CR signal disappears after label shuffling, but the chronological MAE is worse than v7 best. This config should not replace the v7 CR baseline.

## Decision

Do not promote GSE60012 to the main benchmark yet.

Keep:

- GSE60012 header-derived metadata as a valid auxiliary asset.
- GSE60012 held-out scripts/results as a stress test for schema and coverage transfer.

Do not claim:

- 5-dataset GroupKFold improvement.
- GSE60012 biological-age validation.
- Castration/hormone intervention effect.

Recommended next step:

1. Preserve v7 4-dataset matrix as the main benchmark.
2. If continuing with GSE60012, test an explicit outer-join benchmark family with train-fold feature presence filtering and report it separately from strict-common-region GroupKFold.
3. Prioritize coverage-aware region harmonization before model search: GSE60012 has too little strict overlap with the current 4-dataset common matrix.

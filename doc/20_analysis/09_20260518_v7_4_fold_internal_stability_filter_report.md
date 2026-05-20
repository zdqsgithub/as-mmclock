# v7.4 Fold-Internal Stability Filter 报告

Date: 2026-05-18

## Summary

v7.4 按 v7.3 的结论回到 4-dataset strict-common 主矩阵，不继续 outer-join，不扩大 deep learning/embedding 搜索。核心改动是在每个 CV train fold 内增加 dataset-level coverage/mean-shift stability filter：只允许在训练 fold 内各 dataset 覆盖稳定、beta 均值漂移较小的 region 进入后续 age-correlation feature selection。

结果：v7.4 best 为 `lgbm robust p=0.95 top1000 log1p + group_presence=0.95 + max_shift=0.15`。它相对 v7 best 的 GroupKFold MAE 改善 0.599 周，Pearson r 从 0.5798 提升到 0.6240，GSE80672 CR held-out MAE 从 28.441 周改善到 26.543 周，CR AUC 从 0.6616 提升到 0.7649。random-label 和 shuffled-intervention sanity 均通过。

结论：stability filter 值得保留为下一轮候选基线，但不满足 MAE 改善 >=2 周的阈值，因此不启动 30-50 config autoresearch。下一步应做更细的 coverage-aware/region definition harmonization，而不是扩大模型复杂度。

## Implemented Changes

- Updated `scripts/train/train_clock.py`
  - added `--min_train_group_feature_presence`
  - added `--max_train_dataset_mean_shift`
  - stability filtering is fit only on each train fold
  - feature selection still runs after the train-only stability/presence mask
- Updated `scripts/train/train_heldout_clock.py`
  - added the same train-only stability filter for held-out validation
  - output records `n_features_stability`
- Added v7.4 runner:
  - `scripts/train/run_v7_4_stability_smoke.py`
- Outputs:
  - `results/benchmark_v7_4_stability/`
  - `results/validation_v7_4_gse80672_cr_all_except/`

## Smoke Grid

Matrix: `results/multidataset/all_rrbs_region_matrix_5kb.parquet`

Metadata: `metadata/model_sample_metadata_v7_1.csv`

CV: `GroupKFold(dataset_batch)`

| config | r | MAE w | R2 | cross-dataset MAE | CR AUC | train-fold stable features mean |
|---|---:|---:|---:|---:|---:|---:|
| lgbm robust p=0.95 top1000 shift=0.15 | 0.6240 | 23.793 | 0.3170 | 24.707 | 0.7885 | 32,659 |
| lgbm quantile p=0.8 top1000 shift=0.15 | 0.6144 | 23.939 | 0.3027 | 24.721 | 0.8428 | 34,770 |
| lgbm standard p=0.5 top500 no stability | 0.5798 | 24.392 | 0.3079 | 26.115 | 0.6547 | N/A |
| lgbm standard p=0.5 top500 shift=0.15 | 0.6026 | 24.561 | 0.3009 | 24.977 | 0.6572 | 34,770 |
| lgbm standard p=0.5 top500 shift=0.10 | 0.5958 | 24.991 | 0.2523 | 24.622 | 0.7494 | 24,087 |

Comparison to v7:

| benchmark | matrix family | r | MAE w | R2 | cross-dataset MAE | CR AUC |
|---|---|---:|---:|---:|---:|---:|
| v7 best | 4-dataset strict-common | 0.5798 | 24.392 | 0.3079 | 26.115 | 0.6547 |
| v7.4 best | 4-dataset strict-common + train-fold stability | 0.6240 | 23.793 | 0.3170 | 24.707 | 0.7885 |

Delta: MAE improves by 0.599 weeks, r improves by 0.0442, cross-dataset MAE improves by 1.408 weeks, and CR AUC improves by 0.1338.

## Per-Dataset GroupKFold Metrics

Best v7.4 config: `lgbm robust p=0.95 top1000 shift=0.15`.

| held-out dataset | n | r | MAE w | RMSE w | R2 |
|---|---:|---:|---:|---:|---:|
| GSE120137 | 549 | 0.5157 | 21.636 | 26.459 | 0.1609 |
| GSE121141 | 81 | 0.3162 | 39.239 | 50.355 | -0.6993 |
| GSE80672 | 255 | 0.7635 | 26.543 | 32.704 | 0.5008 |
| GSE93957 | 62 | 0.8648 | 11.410 | 13.183 | 0.1846 |

GSE121141 remains the largest failure mode. The stability filter helps overall signal and GSE80672, but it does not solve dataset-specific age/schema/tissue mismatch.

## Sanity Checks

Random-label GroupKFold on v7.4 best config:

| run | r | MAE w | R2 | cross-dataset MAE |
|---|---:|---:|---:|---:|
| real labels | 0.6240 | 23.793 | 0.3170 | 24.707 |
| randomized labels | 0.0021 | 36.026 | -0.4952 | 36.540 |

Pearson r collapses and MAE rises sharply, so this passes the leakage sanity check.

GSE80672 held-out CR validation:

| run | r | MAE w | R2 | CR AUC | CR F1 | Cohen d | Mann-Whitney p |
|---|---:|---:|---:|---:|---:|---:|---:|
| v7 CR baseline | 0.7322 | 28.441 | 0.4331 | 0.6616 | 0.3019 | 0.4813 | N/A |
| v7.4 best | 0.7635 | 26.543 | 0.5008 | 0.7649 | 0.3373 | 1.0137 | 0.000001 |
| v7.4 shuffled intervention | 0.7635 | 26.543 | 0.5008 | 0.4802 | 0.1807 | -0.0554 | 0.641571 |

The CR signal improves and collapses after label shuffling. This supports retaining v7.4 as the current CR validation candidate, but it is still research-grade only and not a biological-age claim.

## Decision

Promote v7.4 stability filtering as a candidate preprocessing option, but do not start broad v7.4 autoresearch.

Reasons:

- GroupKFold MAE improves only 0.599 weeks, below the pre-set 2-week threshold.
- GSE80672 held-out MAE and CR AUC improve, and sanity checks pass.
- GSE121141 remains poor, indicating unresolved dataset/schema/tissue mismatch rather than simple model underfitting.
- The best improvements come from train-fold feature stability, not from adding model complexity.

Recommended next step:

1. Build dataset-pair feature overlap and region stability diagnostics for GSE121141 vs the other datasets.
2. Test stricter region definitions or coverage-aware region aggregation, still under train-fold-only evaluation.
3. Keep v7.4 best as the current CR held-out candidate while avoiding biological-age claims until more held-out intervention datasets are available.


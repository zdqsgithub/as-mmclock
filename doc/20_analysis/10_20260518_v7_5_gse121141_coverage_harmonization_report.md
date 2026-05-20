# v7.5 GSE121141 Coverage-Aware Harmonization 报告

Date: 2026-05-18

## Summary

v7.5 针对 v7.4 的主要失败点 GSE121141 做了两件事：第一，生成 GSE121141-focused coverage/shift/age-range 诊断；第二，在严格无泄漏的 GroupKFold 框架内测试 train-fold-only age-bin coverage filter。

结果：GSE121141 的问题不是 strict-common 矩阵中的简单缺失率过高。它在 72,079 个 strict-common regions 上覆盖率很高，`GSE121141_presence_median=1.0`，但预测误差和年龄高度相关，尤其 52 周以上和 104 周以上样本被系统性预测过年轻。v7.5 best 只把 GSE121141 MAE 从 39.239 周降到 38.033 周，整体 MAE 从 23.793 周降到 23.767 周，改善很小。

结论：age-bin coverage filter 通过 sanity，但不能解决 GSE121141 泛化失败。下一步不应扩大 LGBM/embedding/deep learning 搜索，而应转向更基础的数据层 harmonization：region aggregation 方式、coverage/count 层面的跨 schema 标准化、以及 GSE121141 老年段的组织/批次混杂诊断。

## Implemented Changes

- Updated `scripts/train/train_clock.py`
  - added `--min_train_agebin_feature_presence`
  - added `--min_train_agebin_samples`
  - age-bin coverage mask is fit only on train fold
  - random-label sanity also permutes the age labels used by target-dependent filters
- Updated `scripts/train/train_heldout_clock.py`
  - added the same train-only age-bin coverage filter
- Added diagnostics:
  - `scripts/validate/v7_5_gse121141_harmonization_audit.py`
- Added smoke runner:
  - `scripts/train/run_v7_5_agebin_coverage_smoke.py`
- Outputs:
  - `results/benchmark_v7_5_gse121141_harmonization/`
  - `results/benchmark_v7_5_gse121141_harmonization/qc/`
  - `results/validation_v7_5_gse121141_all_except/`
  - `results/validation_v7_5_gse80672_cr_all_except/`

## GSE121141 Diagnostics

Strict-common matrix:

| metric | value |
|---|---:|
| samples | 947 |
| regions | 72,079 |
| GSE121141 samples | 81 |
| all-data stable regions, p>=0.95 shift<=0.15 | 21,561 |
| all-data stable regions, p>=0.80 shift<=0.15 | 23,092 |
| GSE121141 mean region presence | 0.9930 |
| GSE121141 median region presence | 1.0000 |
| GSE121141 regions with abs mean shift >0.15 vs others | 2,695 |
| GSE121141 regions with abs mean shift >0.10 vs others | 9,532 |

Pairwise raw region overlap:

| pair | intersection regions | jaccard |
|---|---:|---:|
| GSE120137 x GSE121141 | 77,227 | 0.5233 |
| GSE121141 x GSE80672 | 118,304 | 0.7515 |
| GSE121141 x GSE93957 | 137,872 | 0.7616 |

Interpretation: GSE121141 has enough region overlap and high strict-common coverage. The failure is more consistent with age/tissue/batch/schema mismatch than with simple missingness.

Sample-level error correlations under v7.4 best:

| dataset | missing vs abs error r | age vs abs error r | MAE w |
|---|---:|---:|---:|
| GSE120137 | -0.0345 | 0.5558 | 21.636 |
| GSE121141 | 0.0223 | 0.8174 | 39.239 |
| GSE80672 | 0.0256 | 0.5329 | 26.543 |
| GSE93957 | -0.1044 | -0.2456 | 11.410 |

For GSE121141, sample missingness is not associated with error, but age is strongly associated with error.

## Smoke Grid

Matrix: `results/multidataset/all_rrbs_region_matrix_5kb.parquet`

Metadata: `metadata/model_sample_metadata_v7_1.csv`

CV: `GroupKFold(dataset_batch)`

| config | r | MAE w | R2 | cross-dataset MAE | CR AUC | GSE121141 MAE w |
|---|---:|---:|---:|---:|---:|---:|
| lgbm quantile p=0.8 top1000 agebin=0.8 | 0.6170 | 23.767 | 0.3135 | 24.328 | 0.8180 | 38.033 |
| lgbm robust p=0.95 top2000 agebin=0.8 | 0.6198 | 24.631 | 0.2682 | 24.855 | 0.7485 | 38.963 |
| v7.4 best, no agebin | 0.6240 | 23.793 | 0.3170 | 24.707 | 0.7885 | 39.239 |
| lgbm robust p=0.95 top1000 agebin=0.8 | 0.6240 | 23.793 | 0.3170 | 24.707 | 0.7885 | 39.239 |
| lgbm robust p=0.95 top1000 agebin=0.95 | 0.6122 | 24.069 | 0.3062 | 24.935 | 0.8062 | 39.300 |

The best v7.5 config improves GSE121141 MAE by 1.206 weeks and total MAE by 0.026 weeks versus v7.4 best. This is below the threshold for a larger search.

## GSE121141 Held-Out

Best v7.5 config: `lgbm quantile_uniform p=0.8 top1000 agebin=0.8`.

`all_except:GSE121141 -> GSE121141`:

| r | MAE w | RMSE w | R2 |
|---:|---:|---:|---:|
| 0.3580 | 38.033 | 48.778 | -0.5945 |

GSE121141 remains poor. Age-stratified errors show systematic underprediction for old samples:

| tissue | age bin | n | MAE w | residual mean |
|---|---|---:|---:|---:|
| brain_cortex | 104w+ | 5 | 86.108 | 86.108 |
| heart | 104w+ | 5 | 96.501 | 96.501 |
| liver | 104w+ | 4 | 35.647 | 35.005 |
| lung | 104w+ | 5 | 75.342 | 75.342 |
| brain_cortex | 52-104w | 10 | 39.935 | 39.935 |
| heart | 52-104w | 10 | 34.126 | 34.126 |
| liver | 52-104w | 10 | 29.256 | 21.805 |
| lung | 52-104w | 10 | 38.655 | 30.808 |

Positive residual means true age is older than predicted age. The model is compressing older GSE121141 samples toward younger predictions.

## Sanity And CR Validation

Random-label sanity on v7.5 best:

| run | r | MAE w | R2 | cross-dataset MAE |
|---|---:|---:|---:|---:|
| real labels | 0.6170 | 23.767 | 0.3135 | 24.328 |
| randomized labels | -0.0215 | 36.562 | -0.5346 | 37.339 |

This passes leakage sanity.

GSE80672 held-out CR validation:

| run | r | MAE w | R2 | CR AUC | CR F1 | Cohen d |
|---|---:|---:|---:|---:|---:|---:|
| v7.4 best | 0.7635 | 26.543 | 0.5008 | 0.7649 | 0.3373 | 1.0137 |
| v7.5 best | 0.7303 | 26.706 | 0.4856 | 0.8072 | 0.3373 | 1.2558 |
| v7.5 shuffled intervention | 0.7303 | 26.706 | 0.4856 | 0.5284 | 0.2289 | 0.1179 |

v7.5 improves CR AUC and Cohen's d, and shuffled labels collapse toward chance, but chronological held-out MAE is slightly worse than v7.4.

## Decision

Do not start v7.5 autoresearch.

Keep:

- v7.4 best as the main chronological held-out baseline.
- v7.5 best as an alternate CR-oriented configuration, because CR AUC improves and shuffled-intervention sanity passes.
- v7.5 diagnostics as evidence that GSE121141 failure is not primarily simple coverage missingness.

Do not claim:

- age-bin coverage filtering solves cross-dataset generalization;
- GSE121141 generalization is acceptable;
- CR AUC alone is enough to promote v7.5 over v7.4.

Recommended next step:

1. Inspect GSE121141 old-age samples and source metadata for tissue/age/batch coupling.
2. Build coverage/count-level harmonization diagnostics from raw per-sample `.cov.gz`, not only beta regions.
3. Test alternative region aggregation such as coverage-weighted 5kb means or age/tissue-balanced region stability, still fit only inside train folds.
4. Keep model class fixed; the current failure is a data/schema harmonization problem, not evidence for deep learning.


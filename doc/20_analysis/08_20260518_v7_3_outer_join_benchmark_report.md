# v7.3 Outer-Join Benchmark Family 报告

Date: 2026-05-18

## Summary

v7.3 按 v7.2 建议建立了独立的 outer-join benchmark family。它不替代 v7 strict-common 4-dataset 主 benchmark，也不用于证明 GSE60012 已成为主线数据集。目标是测试：在保留所有 dataset union regions，并依赖 train-fold feature presence filtering 的情况下，是否能缓解 GSE60012 strict-overlap 太低的问题。

结果：outer-join 可以运行 5-dataset GroupKFold，但性能比 v7 主线更差。v7.3 best MAE 为 27.006 周，低于 v7 best 的 24.392 周；GSE80672 CR held-out 也明显退化，CR AUC 接近 chance。因此 outer-join 不应进入 autoresearch，也不应替代 v7 主 benchmark。

## Implemented Changes

- Added v7.3 runner:
  - `scripts/train/run_v7_3_outer_join_smoke.py`
- Built outer-join matrix:
  - `results/multidataset_v7_3_outer/all_rrbs_region_matrix_5kb.parquet`
  - shape: 243,026 regions x 1,099 samples
  - datasets: GSE120137, GSE80672, GSE93957, GSE121141, GSE60012
- Outputs:
  - `results/benchmark_v7_3_outer_join/`
  - `results/benchmark_v7_3_outer_join/qc/`
  - `results/validation_v7_3_outer_gse80672_cr_all_except/`

## Outer Matrix QC

The outer matrix keeps many dataset-specific regions, but this creates high missingness:

| dataset | n | sample missing fraction | mean abs region shift |
|---|---:|---:|---:|
| GSE120137 | 549 | 0.6846 | 0.0664 |
| GSE121141 | 81 | 0.4105 | 0.0771 |
| GSE60012 | 152 | 0.5682 | 0.1011 |
| GSE80672 | 255 | 0.4754 | 0.0705 |
| GSE93957 | 62 | 0.3022 | 0.0972 |

PCA remains dataset-dominated:

| component | variance | dataset eta2 | tissue eta2 | age-bin eta2 |
|---|---:|---:|---:|---:|
| PC1 | 0.1536 | 0.9158 | 0.4738 | 0.1085 |
| PC2 | 0.0801 | 0.9379 | 0.4367 | 0.5737 |
| PC3 | 0.0628 | 0.9064 | 0.5982 | 0.1748 |
| PC4 | 0.0472 | 0.6819 | 0.3687 | 0.0274 |

This confirms outer-join mainly converts low overlap into high missingness rather than removing batch structure.

## Smoke Grid

| config | r | MAE w | R2 | cross-dataset MAE | CR AUC | train-fold presence mean |
|---|---:|---:|---:|---:|---:|---:|
| lgbm quantile p=0.8 top1000 log1p | 0.4990 | 27.006 | 0.0712 | 26.147 | 0.5020 | 64,338 |
| lgbm robust p=0.95 top1000 log1p | 0.4304 | 27.928 | 0.0160 | 27.080 | 0.5215 | 35,429 |
| lgbm standard p=0.5 top500 log1p | 0.3278 | 29.550 | -0.0807 | 28.766 | 0.5148 | 87,852 |

Comparison:

| benchmark | datasets | matrix family | r | MAE w | R2 | CR AUC |
|---|---:|---|---:|---:|---:|---:|
| v7 best | 4 | strict-common | 0.5798 | 24.392 | 0.3079 | 0.6547 |
| v7.3 best | 5 | outer-join | 0.4990 | 27.006 | 0.0712 | 0.5020 |

v7.3 is worse by 2.614 weeks MAE and loses the weak CR signal seen in v7.

## Leave-One-Dataset-Out

Best v7.3 model: `lgbm quantile_uniform p=0.8 top1000 log1p`.

| held-out dataset | n | r | MAE w | RMSE w | R2 |
|---|---:|---:|---:|---:|---:|
| GSE120137 | 549 | 0.3769 | 25.632 | 31.387 | -0.1808 |
| GSE121141 | 81 | 0.3667 | 44.359 | 55.633 | -1.0742 |
| GSE60012 | 152 | 0.1919 | 15.375 | 18.711 | -3.6376 |
| GSE80672 | 255 | 0.6477 | 35.570 | 45.078 | 0.0516 |
| GSE93957 | 62 | 0.7620 | 9.798 | 11.863 | 0.3398 |

GSE60012 improves relative to v7.2 train/test strict intersection MAE for some configs, but R2 remains strongly negative. GSE80672 and GSE121141 degrade substantially.

## Sanity And CR Validation

Random-label sanity on v7.3 best:

| run | r | MAE w | R2 |
|---|---:|---:|---:|
| real labels | 0.4990 | 27.006 | 0.0712 |
| randomized labels | 0.0966 | 29.779 | -0.0701 |

Pearson r collapses under random labels and MAE worsens modestly. This passes a minimal leakage check, but the real-label signal remains weak.

GSE80672 CR held-out with v7.3 best:

| run | r | MAE w | R2 | CR AUC | CR F1 | Cohen d |
|---|---:|---:|---:|---:|---:|---:|
| v7 best CR baseline | 0.7322 | 28.441 | 0.4331 | 0.6616 | 0.3019 | 0.4813 |
| v7.3 outer best | 0.6477 | 35.570 | 0.0516 | 0.5077 | 0.2209 | 0.0495 |
| v7.3 shuffled intervention | 0.6477 | 35.570 | 0.0516 | 0.4418 | 0.1963 | -0.1747 |

CR validation is worse than v7 and should not be used as a biological-age baseline.

## Decision

Do not promote v7.3 outer-join to the main benchmark and do not start outer-join autoresearch.

Keep:

- v7 strict-common 4-dataset matrix as the current main benchmark.
- GSE60012 header-derived metadata and auxiliary held-out results as stress-test assets.

Do not claim:

- outer-join improves cross-dataset generalization;
- GSE60012 improves CR validation;
- 5-dataset benchmark supersedes v7.

Recommended next step:

1. Move to coverage-aware region harmonization instead of broader model search.
2. Build dataset-pair/leave-one-dataset feature overlap diagnostics to identify stable region subsets before training.
3. Consider region definitions beyond fixed 5kb bins, but keep evaluation under v7 strict-common and held-out datasets.

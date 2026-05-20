# v7.6 Coverage-Weighted Region Harmonization 报告

Date: 2026-05-18

## Summary

v7.6 按 v7.5 的建议进入 raw coverage/count 层，不扩大模型、不做 embedding/deep learning。目标是测试：把 processed supplements 中的 methylated/total coverage 直接聚合成 coverage-weighted 5kb region beta，是否能缓解 GSE121141 老年样本系统性预测过年轻的问题。

结果是否定的。coverage-weighted 聚合成功构建并通过 random-label sanity，但没有改善 GSE121141，也明显退化 GSE80672 CR held-out。v7.6 best overall 为 `lgbm quantile_uniform p=0.8 top1000 agebin=0.8`，GroupKFold MAE 为 24.985 周，差于 v7.4 的 23.793 周和 v7.5 的 23.767 周；GSE121141 held-out MAE 为 39.038 周，差于 v7.5 的 38.033 周；GSE80672 CR held-out MAE 退化到 30.203 周，CR AUC 降到 0.7354。

结论：不要推广 coverage-weighted matrix，也不要基于 v7.6 开始 autoresearch。GSE121141 的主要问题仍不是 coverage depth 或 CpG count，而是年龄/组织/批次/schema 的系统性偏移，尤其老年样本被压缩到年轻预测。

## Implemented Changes

- Added coverage-weighted ETL:
  - `scripts/etl/14_build_coverage_weighted_region_matrix.py`
  - reads local GEO supplement tar files
  - aggregates methylated and total coverage directly to 5kb windows
  - outputs coverage-weighted region beta, stats, and optional coverage/count matrices
- Added v7.6 multidataset builder:
  - `scripts/etl/15_build_v7_6_weighted_multidataset_matrix.py`
  - uses weighted matrices where available
  - marks GSE120137 as `unweighted_fallback` because raw counts are not available in the current project state
- Added coverage/count audit:
  - `scripts/validate/v7_6_coverage_count_audit.py`
- Added v7.6 benchmark runner:
  - `scripts/train/run_v7_6_weighted_smoke.py`

## Matrix Outputs

Coverage-weighted datasets:

| dataset | schema | regions | samples | source |
|---|---|---:|---:|---|
| GSE80672 | percentage + coverage | 128,739 | 255 | coverage-weighted |
| GSE93957 | Bismark cov | 172,468 | 62 | coverage-weighted |
| GSE121141 | Bismark cov | 148,312 | 81 | coverage-weighted |
| GSE120137 | phase0 region matrix | 77,826 | 549 | unweighted fallback |

Combined v7.6 matrix:

| metric | value |
|---|---:|
| matrix | `results/multidataset_v7_6_weighted/all_rrbs_region_matrix_5kb.parquet` |
| samples | 947 |
| common regions | 72,122 |
| datasets | GSE120137, GSE80672, GSE93957, GSE121141 |

The common region count is essentially unchanged from v7/v7.5, so the benchmark is directly comparable.

## Coverage/Count Audit

GSE121141 coverage/count audit:

| metric | value |
|---|---:|
| samples | 81 |
| weighted regions | 148,312 |
| old-age samples >=104w | 19 |
| mean region total coverage | 1128.29 |
| mean region CpG count | 23.71 |
| regions with abs log2 old/non-old coverage ratio >1 | 35 |

Coverage/error correlations:

| variable | abs-error r | p |
|---|---:|---:|
| mean region total coverage | -0.1557 | 0.1651 |
| median region total coverage | -0.1255 | 0.2643 |
| mean region CpG count | -0.0301 | 0.7898 |
| region presence fraction | -0.0085 | 0.9400 |
| age weeks | 0.8164 | 1.59e-20 |

This confirms the v7.5 diagnosis: GSE121141 error is strongly age-associated and not explained by sample coverage depth, CpG count, or region presence.

## Smoke Grid

Matrix: `results/multidataset_v7_6_weighted/all_rrbs_region_matrix_5kb.parquet`

CV: `GroupKFold(dataset_batch)`

| config | r | MAE w | R2 | cross-dataset MAE | CR AUC | GSE121141 MAE w |
|---|---:|---:|---:|---:|---:|---:|
| weighted standard p=0.5 top500 | 0.5386 | 25.412 | 0.1876 | 26.637 | 0.7614 | 38.855 |
| weighted v7.5 quantile p=0.8 top1000 agebin=0.8 | 0.6004 | 24.985 | 0.2639 | 25.489 | 0.7231 | 39.038 |
| weighted robust p=0.95 top1000 agebin=0.95 | 0.5910 | 25.038 | 0.2459 | 25.917 | 0.7091 | 39.434 |
| weighted v7.4 robust p=0.95 top1000 | 0.5886 | 25.363 | 0.2273 | 26.240 | 0.6662 | 39.965 |

Comparison to previous baselines:

| benchmark | matrix family | r | MAE w | GSE121141 MAE w | CR AUC |
|---|---|---:|---:|---:|---:|
| v7.4 best | unweighted strict-common | 0.6240 | 23.793 | 39.239 | 0.7885 |
| v7.5 best | unweighted strict-common + agebin filter | 0.6170 | 23.767 | 38.033 | 0.8180 |
| v7.6 best overall | coverage-weighted where available | 0.6004 | 24.985 | 39.038 | 0.7231 |

v7.6 is worse than both v7.4 and v7.5 on the headline chronological benchmark.

## Held-Out Sanity

Random-label sanity on v7.6 best overall:

| run | r | MAE w | R2 | cross-dataset MAE |
|---|---:|---:|---:|---:|
| real labels | 0.6004 | 24.985 | 0.2639 | 25.489 |
| randomized labels | 0.0574 | 38.523 | -0.7624 | 35.948 |

Leakage sanity passes.

GSE121141 held-out:

| run | r | MAE w | R2 |
|---|---:|---:|---:|
| v7.5 best | 0.3580 | 38.033 | -0.5945 |
| v7.6 best overall | 0.2759 | 39.038 | -0.6839 |

GSE121141 remains poor and gets worse under coverage-weighted aggregation.

GSE80672 CR held-out:

| run | r | MAE w | R2 | CR AUC | CR F1 | Cohen d |
|---|---:|---:|---:|---:|---:|---:|
| v7.4 best | 0.7635 | 26.543 | 0.5008 | 0.7649 | 0.3373 | 1.0137 |
| v7.5 best | 0.7303 | 26.706 | 0.4856 | 0.8072 | 0.3373 | 1.2558 |
| v7.6 best overall | 0.7282 | 30.203 | 0.4062 | 0.7354 | 0.3484 | 0.7984 |
| v7.6 shuffled intervention | 0.7282 | 30.203 | 0.4062 | 0.5349 | 0.2452 | 0.0393 |

The shuffled-intervention sanity behaves correctly, but the real v7.6 CR signal is weaker than v7.5 and chronological error is much worse.

## Decision

Do not promote v7.6 coverage-weighted matrix to the main benchmark.

Keep:

- v7.4 best as the chronological held-out baseline.
- v7.5 best as the CR-oriented alternate.
- v7.6 coverage/count audit outputs as evidence that GSE121141 failure is not primarily coverage depth or region presence.

Do not do:

- v7.6 autoresearch;
- broader LGBM/RF/ElasticNet tuning on coverage-weighted matrix;
- deep learning or embedding expansion based on this result;
- biological-age claims from v7.6 CR metrics.

Recommended next step:

1. Audit GSE121141 source metadata and sample structure for age/tissue/batch coupling, especially old-age samples.
2. Compare GSE121141 old-age beta distributions against same-tissue younger samples and against other datasets at stable regions.
3. Test a sample-level calibration or dataset-level residual correction only as a diagnostic, not as a benchmark claim, because unseen-dataset calibration semantics are risky.
4. If intervention validation remains the priority, prioritize adding another real intervention dataset over further GSE121141 tuning.


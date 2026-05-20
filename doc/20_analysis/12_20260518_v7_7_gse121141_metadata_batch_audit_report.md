# v7.7 GSE121141 Metadata/Batch Coupling Audit 报告

Date: 2026-05-18

## Summary

v7.7 不训练新模型，专门审计 GSE121141 的 metadata、SOFT source、supplement filename token、coverage/count QC 与预测误差。目标是回答 v7.5/v7.6 留下的问题：GSE121141 老年样本为什么被系统性预测过年轻。

结果：GSE121141 失败主要是年龄外推/年龄分布问题，而不是简单 coverage、CpG count、tar 文件大小、flowcell、lane 或单一 tissue 问题。20 月和 30 月样本在所有组织中都被预测过年轻；30 月样本 v7.5 held-out 平均绝对误差为 75.386 周，非 30 月样本为 26.586 周。数值 QC 与误差相关性很弱，只有 age_weeks 与误差强相关。

结论：不要继续围绕 GSE121141 做 coverage/matrix 聚合微调，也不要上深度学习。下一步应回到科学建模层面：要么显式限制当前 benchmark 的年龄适用范围，要么增加更多老年、多组织、跨数据集训练样本；如果继续研究 GSE121141，则应做诊断性的 age-range calibration，而不能把它作为无校准跨数据集泛化成功证据。

## Implemented Changes

- Added audit script:
  - `scripts/validate/v7_7_gse121141_metadata_batch_audit.py`
- Parsed:
  - `metadata/geo_downloads/GSE121141_family.soft.gz`
  - `raw_downloads/geo_supplements/GSE121141/GSE121141_RAW.tar`
  - `results/multidataset/GSE121141_sample_parse_stats.csv`
  - `results/benchmark_v7_6_coverage_weighted/coverage_qc/GSE121141_sample_coverage_count_qc.csv`
  - v7.4/v7.5/v7.6 predictions
- Outputs:
  - `results/benchmark_v7_7_gse121141_metadata_batch_audit/gse121141_sample_source_error_audit.csv`
  - `results/benchmark_v7_7_gse121141_metadata_batch_audit/gse121141_old_age_source_error_table.csv`
  - `results/benchmark_v7_7_gse121141_metadata_batch_audit/error_by_metadata_batch_group.csv`
  - `results/benchmark_v7_7_gse121141_metadata_batch_audit/contingency_*`
  - `results/benchmark_v7_7_gse121141_metadata_batch_audit/numeric_qc_error_correlations.csv`

## Dataset Structure

GSE121141 has 81 model-aligned samples:

| variable | levels |
|---|---|
| ages | 6, 10, 12, 20, 30 months |
| tissues | brain_cortex, heart, liver, lung |
| batch families | DO, KD |
| flowcells | H3V3MCCXY, H53LFCCXY, HFLMFALXX |
| lanes | L1, L2, L3, L4, L5, L6, L8 |
| old-age samples >=104w | 19 |

Age x tissue is balanced except the 10-month group is only 2 samples:

| age months | brain_cortex | heart | liver | lung |
|---:|---:|---:|---:|---:|
| 6 | 5 | 5 | 5 | 5 |
| 10 | 0 | 1 | 1 | 0 |
| 12 | 5 | 5 | 5 | 5 |
| 20 | 5 | 5 | 5 | 5 |
| 30 | 5 | 5 | 4 | 5 |

## Coupling Diagnostics

Contingency purity:

| left | right | row purity | column purity | comment |
|---|---|---:|---:|---|
| age_months | tissue | 0.2593 | 0.2469 | age is not confounded with tissue |
| age_months | file_batch_family | 0.8025 | 0.3210 | DO/KD is age-skewed but not uniquely old-age |
| age_months | file_flowcell | 0.5432 | 0.3086 | weak/moderate coupling |
| age_months | file_lane | 0.3210 | 0.3333 | weak coupling |
| tissue | file_flowcell | 0.5432 | 0.3457 | weak/moderate coupling |
| tissue | file_batch_family | 0.7778 | 0.3457 | tissue distribution differs by DO/KD |

Age x batch family:

| age months | DO | KD |
|---:|---:|---:|
| 6 | 5 | 15 |
| 10 | 2 | 0 |
| 12 | 0 | 20 |
| 20 | 5 | 15 |
| 30 | 6 | 13 |

Batch family is not sufficient to explain failure: 30-month samples occur in both DO and KD, and 20-month samples also occur in both.

## Error Pattern

v7.5 held-out error by age:

| age months | n | MAE w | residual mean w |
|---:|---:|---:|---:|
| 6 | 20 | 10.658 | -6.339 |
| 10 | 2 | 7.723 | 7.723 |
| 12 | 20 | 21.887 | 14.238 |
| 20 | 20 | 49.099 | 49.099 |
| 30 | 19 | 75.386 | 75.251 |

Positive residual means true age is older than predicted age. The model systematically underpredicts older GSE121141 samples.

v7.5 held-out error by tissue:

| tissue | n | MAE w | residual mean w |
|---|---:|---:|---:|
| brain_cortex | 20 | 42.611 | 41.599 |
| heart | 21 | 42.779 | 36.679 |
| liver | 20 | 24.288 | 16.319 |
| lung | 20 | 42.217 | 32.827 |

Liver is less severe, but every tissue shows strong old-age underprediction.

Old-age v7.5 held-out error by tissue:

| age months | tissue | n | MAE w | residual mean w |
|---:|---|---:|---:|---:|
| 20 | brain_cortex | 5 | 56.264 | 56.264 |
| 20 | heart | 5 | 49.833 | 49.833 |
| 20 | liver | 5 | 40.852 | 40.852 |
| 20 | lung | 5 | 49.448 | 49.448 |
| 30 | brain_cortex | 5 | 86.108 | 86.108 |
| 30 | heart | 5 | 96.501 | 96.501 |
| 30 | liver | 4 | 35.647 | 35.005 |
| 30 | lung | 5 | 75.342 | 75.342 |

This is an age-range compression problem across tissues, not a single tissue artifact.

## QC/Error Correlations

Correlations with absolute error:

| run | variable | r | interpretation |
|---|---|---:|---|
| v74_groupkfold | age_weeks | 0.8174 | dominant association |
| v75_heldout | age_weeks | 0.8164 | dominant association |
| v76_weighted_heldout | age_weeks | 0.7842 | dominant association |
| v75_heldout | rows_pass_coverage | -0.0241 | no meaningful association |
| v75_heldout | coverage_pass_fraction | -0.0369 | no meaningful association |
| v75_heldout | tar_size_bytes | -0.0514 | no meaningful association |
| v75_heldout | mean_region_total_coverage | -0.1557 | weak, not significant |
| v75_heldout | mean_region_cpg_count | -0.0301 | no meaningful association |
| v75_heldout | region_presence_fraction | -0.0085 | no meaningful association |

The coverage/count evidence agrees with v7.6: coverage depth is not the primary driver.

## Decision

Do not run another matrix/feature/filter search targeting GSE121141 right now.

Keep:

- v7.4 best as the chronological baseline.
- v7.5 best as the CR-oriented alternate.
- v7.7 audit tables as the source-of-truth diagnosis for GSE121141 failure.

Do not claim:

- GSE121141 has acceptable held-out age generalization.
- coverage/count harmonization fixed old-age prediction.
- batch family, flowcell, lane, or tissue alone explains the failure.

Recommended next step:

1. Add an explicit age-range diagnostic to benchmark reports, including MAE by age bin and a warning when held-out age exceeds the reliable range learned from other datasets.
2. Run a diagnostic-only age calibration experiment for GSE121141 residuals, clearly labeled as post-hoc and not a valid unseen-dataset benchmark.
3. Prioritize adding more old-age multi-tissue processed datasets or intervention datasets over more model tuning.
4. If model development continues before new data, test conservative target handling such as age-bin stratified loss/sample weighting inside train folds, with random-label sanity.


# v7.8 Age-Range Diagnostics and Calibration Probe 报告

Date: 2026-05-18

## Summary

v7.8 把年龄适用范围诊断正式加入 benchmark 层，不训练新模型，也不扩大搜索。它读取 v7.4/v7.5/v7.6 的现有 predictions，生成 age-range flags、同组织训练年龄支持、年龄段误差、预测年龄压缩斜率，并对 GSE121141 做明确标注为 diagnostic-only 的 post-hoc calibration probe。

结果：GSE121141 的失败来自同组织老年样本外推和预测年龄压缩。全局训练集中有老年血液样本，但 GSE121141 held-out 的 brain/heart 同组织训练最大年龄只有 41 周，liver/lung 最大 86.9 周；GSE121141 130.4 周样本全部超过同组织训练年龄支持。v7.5 held-out 的 GSE121141 预测斜率只有 0.1997，104w+ MAE 为 75.386 周。

结论：当前模型不能被声明为覆盖完整 mouse lifespan 的跨数据集 clock。v7.4/v7.5 可保留为研究级 baseline，但报告必须带 age-range warning。下一步如果继续建模，应优先补充 old-age multi-tissue 训练数据，或者把 benchmark 明确限制在有同组织训练支持的年龄范围内。

## Implemented Changes

- Added diagnostic script:
  - `scripts/validate/v7_8_age_range_diagnostics.py`
- Inputs:
  - v7.4 GroupKFold best predictions
  - v7.5 GSE121141 held-out predictions
  - v7.5 GSE80672 CR held-out predictions
  - v7.6 GSE121141 weighted held-out predictions
  - `metadata/model_sample_metadata_v7_1.csv`
- Outputs:
  - `results/benchmark_v7_8_age_range_diagnostics/predictions_with_age_range_flags.csv`
  - `age_range_support_by_dataset_tissue_agebin.csv`
  - `prediction_age_compression_summary.csv`
  - `age_range_warning_flags.csv`
  - `posthoc_oracle_linear_calibration_metrics.csv`
  - `posthoc_leave_age_level_out_calibration_metrics.csv`

## Warning Summary

v7.8 generated 21 warning rows:

| warning | count |
|---|---:|
| test_age_exceeds_same_tissue_train_support | 12 |
| prediction_age_range_compression | 6 |
| high_old_age_error | 3 |

These warnings are diagnostic flags; they do not alter model metrics.

## Prediction Compression

| run | dataset | slope pred~true | true range w | pred range w | MAE w |
|---|---|---:|---:|---:|---:|
| v7.4 GroupKFold | GSE121141 | 0.1777 | 104.297 | 108.305 | 39.239 |
| v7.5 held-out | GSE121141 | 0.1997 | 104.297 | 115.847 | 38.033 |
| v7.6 weighted held-out | GSE121141 | 0.1642 | 104.297 | 109.965 | 39.038 |
| v7.5 CR held-out | GSE80672 | 0.4047 | 149.188 | 111.246 | 26.706 |

The GSE121141 slope is far below 1, meaning old samples are compressed toward younger predictions.

## Same-Tissue Age Support

For v7.5 GSE121141 held-out:

| tissue | test age bin | test max w | same-tissue train max w | pct above same-tissue support | MAE w |
|---|---|---:|---:|---:|---:|
| brain_cortex | 104w+ | 130.371 | 41.000 | 1.000 | 86.108 |
| heart | 104w+ | 130.371 | 41.000 | 1.000 | 96.501 |
| liver | 104w+ | 130.371 | 86.914 | 1.000 | 35.647 |
| lung | 104w+ | 130.371 | 86.914 | 1.000 | 75.342 |
| brain_cortex | 52-104w | 86.914 | 41.000 | 1.000 | 39.935 |
| heart | 52-104w | 86.914 | 41.000 | 1.000 | 34.126 |
| liver | 52-104w | 86.914 | 86.914 | 0.500 | 29.256 |
| lung | 52-104w | 86.914 | 86.914 | 0.500 | 38.655 |

This is the key finding: global training age range is misleading because old support is mostly blood, not matched multi-tissue samples.

GSE80672 has a related warning:

| run | dataset | age bin | same-tissue train max w | warning |
|---|---|---|---:|---|
| v7.5 CR held-out | GSE80672 blood | 104w+ | 92.462 | all 104w+ blood samples exceed same-tissue train support |

This means CR biological-age residuals should remain research-grade, with age-range caveats.

## Post-Hoc Calibration Probe

These calibration outputs are invalid as benchmark metrics because they use GSE121141 labels after prediction. They are only used to test whether the failure is a simple linear scale issue.

Oracle in-sample linear calibration:

| run | MAE before w | calibrated MAE w | calibration slope |
|---|---:|---:|---:|
| v7.4 GroupKFold GSE121141 | 39.239 | 31.569 | 0.5626 |
| v7.5 held-out GSE121141 | 38.033 | 30.808 | 0.6415 |
| v7.6 weighted GSE121141 | 39.038 | 32.329 | 0.4636 |

Even oracle in-sample calibration leaves MAE around 31 weeks, so the problem is not only a simple intercept/slope correction.

Leave-one-age-level-out calibration for v7.5 GSE121141:

| held-out age months | n | calibrated MAE w |
|---:|---:|---:|
| 6 | 20 | 56.444 |
| 10 | 2 | 26.595 |
| 12 | 20 | 24.367 |
| 20 | 20 | 21.759 |
| 30 | 19 | 71.696 |

Leaving out 30-month samples still fails badly. Post-hoc age calibration does not solve old-age extrapolation.

## Decision

Promote v7.8 diagnostics into the standard reporting layer for future runs.

Keep:

- v7.4 best as chronological baseline, with age-range warning.
- v7.5 best as CR-oriented alternate, with age-range warning.
- v7.8 age-range flags as mandatory context for GSE121141/GSE80672 claims.

Do not claim:

- full-lifespan cross-dataset mouse age prediction;
- successful GSE121141 old-age generalization;
- validated biological-age effects at old ages without age-range caveats.

Recommended next step:

1. Add old-age multi-tissue processed datasets before further model tuning.
2. If no new data is added, restrict headline benchmark to age/tissue ranges with train support and report old-age extrapolation separately.
3. For CR validation, keep GSE80672 residual analysis but explicitly flag 104w+ blood samples as beyond same-tissue training support.
4. Do not use post-hoc calibration as a benchmark; use it only to demonstrate the failure mode.


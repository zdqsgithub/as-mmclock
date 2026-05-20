# v8.2 Pre-Filter Liftover Harmonization + RALPH Report

Date: 2026-05-18

## Summary

v8.2 修正了 v8.1 的核心问题：UCSC `mm9ToMm10` chain 中 source assembly 位于 chain target fields，destination assembly 位于 query fields；v8.1 的后过滤 point-liftover 不能作为正式判断。v8.2 改为在 raw COV / 100bp tile 解析阶段先做 mm9->mm10 liftover，再做 sample presence filtering 和 5kb region aggregation。

结果分两层：

- 数据层成功：6-dataset strict common regions 从 v8.1 的 13,357 恢复到 65,870，超过 50,000 门禁。
- 模型层未成功：锁定 RALPH smoke 没有改善 v7.5 baseline；GSE121141 old-age held-out 仍很差，random-label sanity 也处于边界失败。

因此 v8.2 不进入 30-50 config constrained autoresearch。

UCSC chain 方向依据：UCSC chain format 定义 `tName/tStart/tEnd` 为 target/reference，`qName/qStart/qEnd` 为 query；当 query strand 为 `-` 时需要用 chromosome size 转回 forward coordinates。参考官方文档：https://genome.ucsc.edu/goldenpath/help/chain.html

## Code Changes

- Updated `scripts/etl/10_convert_processed_methylation.py`
  - Added `--liftover_chain`, `--source_assembly`, `--target_assembly`.
  - Added pre-filter liftover for Bismark COV rows.
  - Added interval-aware liftover for GSE60012 100bp tile rows.
  - Added duplicate coordinate collapse before sample presence filtering.
  - Added assembly/liftover fields to matrix manifests.
- Updated `scripts/etl/18_liftover_gse213628_mm9_to_mm10.py`
  - Corrected UCSC chain direction for future diagnostic use.
  - Kept this as post-filter diagnostic; official v8.2 matrix uses `10_convert_processed_methylation.py`.

## Data Outputs

- `results/multidataset_v8_2_prefilter_liftover/GSE213628_beta_matrix.parquet`
- `results/multidataset_v8_2_prefilter_liftover/GSE213628_region_matrix_5kb.parquet`
- `results/multidataset_v8_2_prefilter_liftover/GSE213628_matrix_manifest.json`
- `results/multidataset_v8_2_prefilter_liftover/GSE60012_beta_matrix.parquet`
- `results/multidataset_v8_2_prefilter_liftover/GSE60012_region_matrix_5kb.parquet`
- `results/multidataset_v8_2_prefilter_liftover/GSE60012_matrix_manifest.json`
- `results/multidataset_v8_2_prefilter_liftover/all_rrbs_region_matrix_5kb.parquet`
- `results/multidataset_v8_2_prefilter_liftover/all_rrbs_matrix_manifest.json`
- `results/multidataset_v8_2_prefilter_liftover/region_overlap_diagnostics.csv`

## Liftover Validation

| Dataset | Input rows | Rows after duplicate collapse | Rows after presence filter | 5kb regions | Notes |
|---|---:|---:|---:|---:|---|
| GSE213628 | 8,820,394 | 8,820,394 | 1,753,535 | 124,695 | raw COV CpG pre-filter liftover |
| GSE60012 | 817,799 | 817,791 | 462,470 | 109,805 | 100bp tile interval pre-filter liftover |

GSE213628 liftover failure rate was very low after fixing chain direction. GSE60012 remains lower confidence because it is tile-level, synthetic-header metadata rather than GSM-level CpG calls.

## Matrix Gate

| Matrix set | Common 5kb regions | Threshold | Status |
|---|---:|---:|---|
| v8 unlifted all | 22,822 | 50,000 | blocked |
| v8.1 post-filter lifted all | 13,357 | 50,000 | blocked/superseded |
| v8.2 pre-filter lifted core without GSE60012 | 70,302 | 50,000 | pass |
| v8.2 pre-filter lifted all six datasets | 65,870 | 50,000 | pass |

Final v8.2 matrix:

- Shape: 65,870 regions x 1,219 samples
- Datasets: GSE120137, GSE80672, GSE93957, GSE121141, GSE60012, GSE213628
- Age-known samples: 1,219
- Intervention labels: control 1,126; CR 32; other_intervention 61

## RALPH Smoke Results

Locked configs only:

| Config | r | MAE weeks | R2 | cross-dataset MAE | CR AUC | GSE121141 MAE | GSE121141 104w+ MAE |
|---|---:|---:|---:|---:|---:|---:|---:|
| lgbm quantile p0.8 top1000 agebin0.8 | 0.6253 | 24.843 | 0.2091 | 23.359 | 0.7567 | 41.616 | 84.611 |
| lgbm robust p0.95 top1000 shift0.15 | 0.6137 | 25.358 | 0.1819 | 24.017 | 0.7151 | 43.195 | 84.425 |

v7.5 baselines:

- GroupKFold MAE: 23.767w
- GSE121141 held-out MAE: 38.033w
- GSE121141 104w+ MAE: 75.386w

v8.2 did not improve the chronological benchmark or the GSE121141 old-age target.

## Sanity Checks

Random-label sanity for selected robust config:

- Pearson r: 0.2005
- MAE: 31.582w
- R2: -0.2728

This is worse than real labels, but `abs(r)<0.2` was the pre-registered threshold and the observed value is 0.2005. Treat this as a boundary failure rather than a clean pass.

GSE80672 shuffled intervention sanity:

- True held-out CR AUC: 0.7424
- Shuffled CR AUC: 0.5083
- True Cohen's d: 0.8891
- Shuffled Cohen's d: -0.0865

The CR signal collapses under shuffled labels, but true held-out CR AUC is below the 0.8 target and below the v7.5 CR-oriented result.

## LODO Metrics

Selected robust config, leave-one-dataset-out:

| Held-out dataset | n | r | MAE weeks | R2 | 104w+ n | 104w+ MAE | CR AUC |
|---|---:|---:|---:|---:|---:|---:|---:|
| GSE120137 | 549 | 0.5710 | 27.315 | -0.4135 | 0 | N/A | N/A |
| GSE80672 | 255 | 0.7605 | 30.092 | 0.3753 | 74 | 53.927 | 0.7424 |
| GSE93957 | 62 | 0.8568 | 7.505 | 0.6074 | 0 | N/A | N/A |
| GSE121141 | 81 | 0.4248 | 37.485 | -0.6373 | 19 | 79.921 | N/A |
| GSE60012 | 152 | 0.7864 | 7.434 | -0.2530 | 0 | N/A | N/A |
| GSE213628 | 120 | 0.6062 | 23.747 | 0.3670 | 42 | 30.532 | N/A |

Interpretation:

- GSE213628 itself is predictable in old-age held-out mode, but adding it does not transfer well to GSE121141 old-age samples.
- GSE121141 remains the main extrapolation blocker.
- GSE120137 and GSE80672 also show substantial held-out error, so the issue is broader than one old-age dataset.

## Decision

Do not run constrained autoresearch.

Reasons:

- GroupKFold MAE worsened versus v7.5.
- GSE121141 104w+ MAE worsened versus v7.5.
- GSE80672 held-out CR AUC is below target and below the v7.5 CR-oriented baseline.
- Random-label sanity is borderline at r=0.2005.

v8.2 solved the assembly-aware matrix blocker, but did not solve the modeling/generalization problem.

## Recommended Next Step

Run v8.3 as a narrow ablation/diagnostic, not autoresearch:

1. Compare v8.2 matrices with and without GSE60012 and with/without GSE213628.
2. Test whether GSE213628 helps only when GSE121141 is held out, or whether it is diluted by tissue/schema differences.
3. Add tissue-aware LODO summaries for shared tissues only.
4. Revisit feature harmonization only after identifying which dataset/tissue pair causes the remaining residual shift.

Deep learning, embedding search, and broad autoresearch should remain paused until this diagnostic explains the residual mismatch.

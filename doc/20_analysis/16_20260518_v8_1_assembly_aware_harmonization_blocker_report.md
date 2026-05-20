# v8.1 Assembly-Aware Harmonization Blocker Report

Date: 2026-05-18

## Summary

v8.1 按 AS-DS-Ops 的无泄漏、可追溯、先过数据门禁再训练原则，完成了 assembly-aware harmonization 预检、mm9->mm10 liftover、5kb region matrix 重建和 strict inner-join matrix 尝试。

结论：v8.1 不能进入 RALPH smoke、GroupKFold benchmark 或 constrained autoresearch。原因是 liftover 后 common 5kb regions 没有恢复，反而低于 v8 unlifted 结果。v8.1-core strict common regions 只有 16,678；包含 lifted GSE60012 后只有 13,357，均低于 50,000 准入阈值。

这说明当前 blocker 不是单纯的“把 GSE213628 从 mm9 lift 到 mm10”就能解决。继续训练会把 coordinate/schema mismatch 当作模型目标优化，不符合项目规范。

## Assembly Evidence

本轮只使用本地 cached GEO SOFT 和官方 UCSC chain，不下载新大数据集。

| Dataset | SOFT/processing assembly evidence | v8.1 action |
|---|---|---|
| GSE120137 | `Genome_build: mm10` | keep native |
| GSE80672 | `Genome_build: GRCm38.p2` | keep native |
| GSE93957 | mapped to mouse `GRCm38`; `Genome_build: GRCm38` | keep native |
| GSE121141 | `GRCm38/mm10`; `Genome_build: mm10` | keep native |
| GSE60012 | mouse samples `mm9`; human samples `hg19` | lifted mm9->mm10 for diagnostic |
| GSE213628 | `Assembly: mm9`; SOFT also contains contradictory `hg19` text | lifted mm9->mm10 using `Assembly: mm9` as source |

The chain file used was:

- `resources/liftover/mm9ToMm10.over.chain.gz`
- Source: UCSC `mm9ToMm10.over.chain.gz`

## Implemented Changes

New/updated code:

- `scripts/etl/18_liftover_gse213628_mm9_to_mm10.py`
  - Parses UCSC chain directly; no new dependency.
  - Performs point liftover on CpG/tile row coordinates.
  - Rebuilds 5kb region matrix after liftover.
  - Supports dataset/source schema arguments, so it can handle both GSE213628 CpG rows and GSE60012 tile-start rows.

Generated v8.1 artifacts:

- `results/multidataset_v8_1_gse213628_lifted/GSE213628_beta_matrix.parquet`
- `results/multidataset_v8_1_gse213628_lifted/GSE213628_liftover_cpg_map.parquet`
- `results/multidataset_v8_1_gse213628_lifted/GSE213628_region_matrix_5kb.parquet`
- `results/multidataset_v8_1_gse213628_lifted/GSE213628_liftover_manifest.json`
- `results/multidataset_v8_1_gse213628_lifted/GSE60012_beta_matrix.parquet`
- `results/multidataset_v8_1_gse213628_lifted/GSE60012_liftover_cpg_map.parquet`
- `results/multidataset_v8_1_gse213628_lifted/GSE60012_region_matrix_5kb.parquet`
- `results/multidataset_v8_1_gse213628_lifted/GSE60012_liftover_manifest.json`
- `results/multidataset_v8_1_gse213628_lifted/region_overlap_diagnostics.csv`
- `results/multidataset_v8_1_gse213628_lifted/all_rrbs_matrix_blocker.json`
- `results/multidataset_v8_1_all_mm10_lifted/all_rrbs_matrix_blocker.json`

## Liftover Validation

| Dataset | Input rows | Lifted rows | Lifted rate | 5kb regions | sex/MT region rows |
|---|---:|---:|---:|---:|---:|
| GSE213628 | 1,753,608 | 995,024 | 56.7% | 70,499 | 0 |
| GSE60012 | 462,500 | 261,177 | 56.5% | 61,598 | 0 |

Important note: GSE60012 is a 100bp tile matrix, so v8.1 maps tile-start coordinates. This is valid as a diagnostic but remains lower confidence than base-resolution CpG liftover.

## Matrix Gate Results

Strict inner-join attempts:

| Matrix set | Common 5kb regions | Threshold | Status |
|---|---:|---:|---|
| v8 unlifted core: GSE120137+GSE80672+GSE93957+GSE121141+GSE213628 | 28,545 | 50,000 | blocked in v8 |
| v8.1 lifted core: GSE120137+GSE80672+GSE93957+GSE121141+GSE213628 | 16,678 | 50,000 | blocked |
| v8 unlifted all: core+GSE60012 | 22,822 | 50,000 | blocked in v8 |
| v8.1 all lifted-mm10: core+GSE60012+GSE213628 | 13,357 | 50,000 | blocked |

Selected overlap diagnostics after liftover:

| Dataset set | Common 5kb regions |
|---|---:|
| GSE80672+GSE93957+GSE121141 | 112,958 |
| GSE120137+GSE80672+GSE93957+GSE121141 | 72,079 |
| GSE60012_lifted_mm10+GSE213628_lifted_mm10 | 51,676 |
| GSE93957+GSE213628_lifted_mm10 | 34,313 |
| GSE93957+GSE121141+GSE213628_lifted_mm10 | 28,814 |
| GSE80672+GSE93957+GSE121141+GSE213628_lifted_mm10 | 24,297 |
| GSE120137+GSE80672+GSE93957+GSE121141+GSE213628_lifted_mm10 | 16,678 |

The main benchmark datasets already have strong mutual overlap without GSE213628. The large drop appears specifically when the lifted old-age GSE213628 matrix is included. This points to remaining schema/library/coordinate incompatibility rather than a model-capacity issue.

## Decision

Stop v8.1 before training.

Skipped by design:

- v8 RALPH smoke configs
- GroupKFold benchmark
- leave-one-dataset-out benchmark
- GSE80672 CR held-out validation
- random-label sanity
- constrained autoresearch

Reason:

- The pre-registered matrix gate was `common_regions >= 50,000`.
- Both v8.1 strict matrices failed the gate.
- Liftover decreased overlap versus the v8 native-coordinate attempt.
- Running autoresearch here would optimize a coordinate harmonization blocker, not a real biological/chronological age objective.

## Interpretation

The likely causes are:

1. GSE213628 SOFT metadata is internally inconsistent: it reports `Assembly: mm9`, but also contains text saying reads were aligned to `hg19`. Treating it as mm9 and lifting to mm10 did not improve common-region overlap.
2. RRBS library design, coverage filtering, and processed-file schemas differ enough that strict 5kb inner join remains too small even after coordinate conversion.
3. Point liftover after matrix-level filtering may lose regions that would be recoverable if liftover were applied earlier, before CpG/coverage filtering and region aggregation.
4. GSE60012 is tile-level mm9 data, not base-resolution CpG data; it should remain auxiliary unless interval-aware tile liftover is implemented.

## Recommended v8.2

Do not do autoresearch yet.

Recommended next step:

1. Add explicit `assembly` and `coordinate_source` fields to every per-dataset matrix manifest.
2. Rebuild mm9 datasets with interval-aware harmonization before presence filtering:
   - GSE213628: liftover raw COV CpG coordinates before CpG presence filtering and 5kb aggregation.
   - GSE60012: liftover 100bp tile intervals, not only tile-start points.
3. Compare three controlled matrix variants:
   - native-coordinate diagnostic,
   - point-lifted diagnostic,
   - pre-filter interval/CpG-lifted rebuilt matrix.
4. Only if the strict common-region count recovers to `>=50,000`, rerun the two v8 RALPH smoke configs.
5. If strict common regions still fail, move to a separate diagnostic track for outer-join/union matrix with train-fold presence filters, clearly marked non-headline until leakage and missingness semantics are validated.

No biological-age or CR conclusion is made from v8.1.

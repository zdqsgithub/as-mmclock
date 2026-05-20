# v8 GSE213628 老龄多组织数据闭环 Blocker 报告

Date: 2026-05-18

## Summary

v8 完成了 `GSE213628` 的官方 GEO 下载、metadata 构建、parser smoke 和全量 5kb region matrix 转换，但在多数据集 strict inner-join matrix 阶段按准入规则停止。

核心结果：`GSE213628` 本身质量可用，metadata 也可追溯；但它与现有 benchmark 矩阵的共同 5kb regions 过少。包含 `GSE60012` 时只有 22,822 个 common regions；排除已知会压低交集的 `GSE60012` 后，v8-core 也只有 28,545 个 common regions，低于计划规定的 50,000 阈值。因此没有继续运行 GroupKFold、held-out、random-label 或 CR validation。

这不是模型搜索问题，不能通过 autoresearch 解决。下一步应先做 assembly/coordinate harmonization，尤其是确认 `GSE213628` 的 SOFT 中 `Assembly: mm9` 与当前 mm10/GRCm38 风格 region 坐标之间的兼容性。

## Implemented Artifacts

- Added `GSE213628` support to GEO inventory/preflight:
  - `scripts/etl/08_geo_supplement_inventory.py`
  - `scripts/etl/11_geo_processed_preflight.py`
- Added v8 metadata builder:
  - `scripts/etl/17_build_gse213628_metadata.py`
- Extended conversion and multidataset builder compatibility:
  - `scripts/etl/10_convert_processed_methylation.py`
  - `scripts/etl/12_build_multidataset_region_matrix.py`
- Added v8 benchmark/report runners for later use after matrix blocker is fixed:
  - `scripts/train/run_v8_gse213628_ralph_smoke.py`
  - `scripts/validate/write_v8_gse213628_report.py`

Generated data artifacts:

- `metadata/gse213628_sample_metadata.csv`
- `metadata/model_sample_metadata_v8.csv`
- `metadata/geo_supplement_inventory_v8_gse213628.csv`
- `metadata/geo_processed_preflight_v8_gse213628.csv`
- `raw_downloads/geo_supplements/GSE213628/GSE213628_RAW.tar`
- `results/parser_smoke_gse213628/`
- `results/multidataset_v8_gse213628/GSE213628_region_matrix_5kb.parquet`
- `results/multidataset_v8_gse213628/all_rrbs_matrix_blocker.json`
- `results/multidataset_v8_gse213628/region_overlap_diagnostics.csv`

## Validation Results

GSE213628 metadata validation passed:

| Check | Result |
|---|---:|
| metadata rows | 120 |
| age-known samples | 120 |
| age coverage | 100% |
| max age | 121.68 weeks |
| tissue count | 8 |
| duplicate sample_id | 0 |

GEO download validation passed:

| File | Status | Verified bytes |
|---|---|---:|
| `GSE213628_RAW.tar` | completed | 3,260,723,200 |

Parser smoke passed:

| Check | Result |
|---|---:|
| smoke samples | 3 |
| metadata overlap | 3 |
| beta min/max | 0.0 / 1.0 |
| parser errors | 0 |

Full conversion passed:

| Check | Result |
|---|---:|
| parsed samples | 120 |
| metadata overlap | 120 |
| CpGs before presence filter | 8,820,735 |
| CpGs after presence filter | 1,753,608 |
| 5kb regions | 124,683 |
| sex/MT region rows | 0 |
| parser errors | 0 |

## Matrix Blocker

Strict inner-join matrix attempts:

| Dataset set | Common 5kb regions | Threshold | Status |
|---|---:|---:|---|
| GSE120137 + GSE80672 + GSE93957 + GSE121141 + GSE60012 + GSE213628 | 22,822 | 50,000 | blocked |
| GSE120137 + GSE80672 + GSE93957 + GSE121141 + GSE213628 | 28,545 | 50,000 | blocked |

Important overlap diagnostics:

| Dataset set | Common 5kb regions |
|---|---:|
| GSE93957 + GSE213628 | 59,717 |
| GSE121141 + GSE213628 | 51,802 |
| GSE93957 + GSE121141 + GSE213628 | 49,878 |
| GSE80672 + GSE213628 | 45,541 |
| GSE120137 + GSE213628 | 30,573 |
| GSE120137 + GSE80672 + GSE93957 + GSE121141 + GSE213628 | 28,545 |

The main loss is not caused only by `GSE60012`. `GSE213628` has limited intersection with `GSE120137`, and the SOFT metadata reports `Assembly: mm9`. This makes a direct strict-coordinate merge with existing matrices scientifically risky.

## Decision

Stop v8 before training. Do not run constrained autoresearch yet.

Reason:

- The v8 acceptance rule required common 5kb regions `>=50,000`.
- The best strict inner-join v8-core matrix has only 28,545 regions.
- Running GroupKFold on this matrix would answer a different question and would not be comparable to v7.4/v7.5.
- Autoresearch would optimize a coordinate/schema blocker rather than a real model objective.

Recommended next step:

1. Implement v8.1 assembly-aware harmonization:
   - confirm assembly for each matrix source;
   - liftover or rebuild `GSE213628` regions into the same mm10/GRCm38 coordinate frame as the main benchmark;
   - rebuild 5kb regions after coordinate harmonization.
2. Re-run strict inner-join validation before any training.
3. Only if common regions recover to `>=50,000`, rerun the v8 RALPH smoke configs.
4. Keep outer-join or pairwise-intersection runs diagnostic-only until coordinate harmonization is resolved.

No biological-age claim or CR validation was made in v8.

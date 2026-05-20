# v7.1 GSE60012 Header Mapping 与 Harmonization 报告

Date: 2026-05-18

## Summary

v7.1 实现了 GSE60012 official tile matrix header-derived metadata，不再要求 tile columns 唯一映射到 GSM。GSE60012 样本被明确标记为 `matrix_header_synthetic_non_gsm`，只能用于 age-prediction benchmark 设计，不能做 GSM 级 GEO 样本声明。

Mapping 阶段通过，但 5-dataset inner-join matrix 阶段按准入规则停止：加入 GSE60012 后共同 5kb regions 只有 25,374，低于计划规定的 50,000 阈值。因此没有继续运行 v7.1 smoke benchmark、random-label、GSE80672 CR held-out 或 GSE60012 held-out age validation。

## Implemented Changes

- Added GSE60012 header-derived metadata builder:
  - `scripts/etl/13_build_gse60012_header_metadata.py`
  - Outputs:
    - `metadata/gse60012_header_sample_metadata.csv`
    - `metadata/model_sample_metadata_v7_1.csv`
    - `results/multidataset_v7_1/GSE60012_header_mapping_exclusions.csv`
- Updated matrix/training/QC entry points with metadata path and sample-id resolver:
  - GSM columns resolve to `GSM...`
  - GSE60012 tile columns resolve to short IDs such as `GSE60012_tile_001`
  - Full tile labels remain traceable in `matrix_sample_id/source_column`
- Added v7.1 fixed smoke runner:
  - `scripts/train/run_v7_1_harmonization_smoke.py`
  - Not executed because matrix validation stopped first.

## Mapping Validation

GSE60012 header mapping result:

| metric | value |
|---|---:|
| Raw tile columns | 154 |
| Parsed synthetic samples | 152 |
| Excluded columns | 2 |
| Duplicate model sample IDs | 0 |
| Parsed control samples | 91 |
| Parsed other_intervention samples | 61 |

Excluded columns:

| sample | source_column | reason |
|---|---|---|
| GSE60012_tile_080 | `Liver_M` | not enough header tokens |
| GSE60012_tile_154 | empty | empty source column |

Required GSE60012 metadata fields are complete for all 152 parsed samples:

- `age_days`
- `age_weeks`
- `tissue`
- `sex`
- `dataset_batch`
- `sample_id_source`

## Matrix Validation Blocker

Attempted v7.1 matrix build:

```bash
python scripts/etl/12_build_multidataset_region_matrix.py \
  --metadata_path metadata/model_sample_metadata_v7_1.csv \
  --out_dir results/multidataset_v7_1 \
  --min_regions 50000
```

Included dataset matrix sizes before intersection:

| dataset | regions | samples |
|---|---:|---:|
| GSE120137 | 77,826 | 549 |
| GSE80672 | 128,739 | 255 |
| GSE93957 | 171,918 | 62 |
| GSE121141 | 146,979 | 81 |
| GSE60012 | 109,777 | 152 |

Blocker:

| check | threshold | observed | status |
|---|---:|---:|---|
| 5-dataset common 5kb regions | >=50,000 | 25,374 | blocked |

The blocker was written to:

- `results/multidataset_v7_1/all_rrbs_matrix_blocker.json`

## Region Overlap Diagnostics

The previous v7 four-dataset intersection remains healthy:

| datasets | common regions |
|---|---:|
| GSE120137 + GSE80672 + GSE93957 + GSE121141 | 72,079 |

Any 4-dataset intersection containing GSE60012 falls sharply:

| datasets | common regions |
|---|---:|
| GSE120137 + GSE80672 + GSE121141 + GSE60012 | 25,560 |
| GSE120137 + GSE80672 + GSE93957 + GSE60012 | 25,506 |
| GSE120137 + GSE93957 + GSE121141 + GSE60012 | 26,728 |
| GSE80672 + GSE93957 + GSE121141 + GSE60012 | 36,746 |
| all 5 datasets | 25,374 |

Full diagnostics:

- `results/multidataset_v7_1/region_overlap_diagnostics.csv`

## GSE60012 Condition Audit

GSE60012 non-NORMAL conditions were retained only as exploratory labels. They were not promoted to castration biological-age validation.

Condition summary:

| condition_family | n |
|---|---:|
| normal | 91 |
| old_castration | 10 |
| old_castration_control | 10 |
| young_castrated | 10 |
| young_control_castrated | 8 |
| testosterone_+_castration | 6 |
| vehicle_+_castration | 6 |
| testosterone | 4 |
| testosterone_control | 4 |
| castration | 3 |

Full audit:

- `results/multidataset_v7_1/GSE60012_condition_audit.csv`

## Decision

Do not run v7.1 smoke benchmark yet. The v7.1 acceptance criteria explicitly says to stop if the 5-dataset common region count is below 50,000. Running GroupKFold on 25,374 strict-common regions would answer a different question and would make v7 vs v7.1 comparison hard to interpret.

Recommended next step:

1. Keep GSE60012 header-derived metadata as valid auxiliary metadata.
2. Do not merge GSE60012 by strict inner join for the main benchmark.
3. Design v7.2 around coverage-aware matrix strategies:
   - 4-dataset v7 baseline matrix remains the main benchmark.
   - GSE60012 can be tested as a held-out auxiliary dataset using feature intersection with train data, but not as a strict 5-dataset inner-join matrix.
   - Alternatively evaluate an outer-join matrix with train-fold presence filtering, explicitly labeled as a different benchmark family.

No biological-age claim was made from GSE60012 in v7.1.

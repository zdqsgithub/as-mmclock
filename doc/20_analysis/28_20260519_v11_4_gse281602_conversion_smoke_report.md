# v11.4 GSE281602 Conversion Smoke Report

Date: 2026-05-19 00:05 CST

Project standard: AS-DS-Ops v3.3. This is a data-gate conversion smoke, not a
training or biological-age benchmark.

## Summary

`GSE281602` was used as the first v11.4 direct-conversion smoke because it is
the smallest verified Bismark-COV dataset in the v11.3 backlog.

Result:

- conversion succeeded;
- all 8 samples were parsed;
- all 8 samples mapped to auxiliary metadata;
- beta values were valid in `[0, 1]`;
- 5kb region matrix was created;
- common-region gate failed against the v8.2 reference matrix.

Decision:

`GSE281602` remains auxiliary only. It validates the packaging and Bismark-COV
conversion path, but it must not enter headline chronological benchmark or
autoresearch because it is cardiomyocyte/cell-type specific and common 5kb
region overlap is below the project gate.

## Outputs

Script:

- `scripts/validate/run_v11_4_gse281602_conversion_smoke.py`

Output directory:

- `results/v11_4_conversions/GSE281602/`

Key files:

- `GSE281602_auxiliary_metadata.csv`
- `GSE281602_verified_processed_files.tar`
- `GSE281602_beta_matrix.parquet`
- `GSE281602_region_matrix_5kb.parquet`
- `GSE281602_region_stats_5kb.csv`
- `GSE281602_matrix_manifest.json`
- `GSE281602_sample_parse_stats.csv`
- `common_regions_with_v8_2.csv`
- `v11_4_gse281602_conversion_smoke_state.json`

## Metadata Gate

| Field | Value |
|---|---:|
| Metadata rows | 8 |
| Metadata overlap | 8 |
| Age-known overlap | 8 |
| Tissue | heart |
| Cell type | cardiomyocytes |
| Sex counts | F=4, M=4 |
| Age weeks | 17.383, 121.680 |
| Intervention | control=8 |
| Headline allowed | false |

Reason not headline:

`GSE281602` is cardiomyocyte/cell-type specific, not bulk heart. It can be used
as auxiliary matrix/parser evidence, but not as a GSE121141-like headline
chronological validation dataset.

## Conversion Gate

| Metric | Value |
|---|---:|
| Input files | 8 |
| Input bytes | 233,585,452 |
| Tar bytes | 233,605,120 |
| Schema | `bismark_cov_per_sample_tar` |
| Min coverage | 5 |
| Sample presence filter | 0.5 |
| Region bin size | 5,000 |
| Region presence filter | 0.8 |
| Samples parsed | 8 |
| Parser errors | 0 |
| CpG rows before presence filter | 3,324,410 |
| CpG rows after presence filter | 2,779,291 |
| Regions | 121,333 |
| sex/MT regions | 0 |
| Conversion runtime | 57.4 sec |

Sample parse stats were consistent across all eight samples: each sample had
approximately 3.8M-4.4M input COV rows, 3.7M-4.2M autosomal parseable rows, and
2.5M-2.9M rows passing coverage.

## Matrix QC

| Metric | Value |
|---|---:|
| Beta matrix shape | 2,779,291 x 8 |
| Region matrix shape | 121,333 x 8 |
| beta min | 0.0 |
| beta max | 1.0 |
| beta values outside `[0,1]` | 0 |

## Region Overlap Gate

Reference matrix:

- `results/multidataset_v8_2_prefilter_liftover/all_rrbs_region_matrix_5kb.parquet`

| Metric | Value |
|---|---:|
| GSE281602 regions | 121,333 |
| Reference regions | 65,870 |
| Common regions | 25,321 |
| Fraction vs GSE281602 | 0.209 |
| Fraction vs reference | 0.384 |
| Common-region gate `>=50,000` | false |

This fails the matrix gate for headline training. The failure is consistent with
`GSE281602` being a cell-type/cardiomyocyte dataset rather than a bulk tissue
dataset designed to align with the current multi-dataset RRBS benchmark.

## Next Step

Use this result to validate the direct COV conversion path, then process the next
ready direct-COV dataset only as a data gate:

1. `GSE224442`: intervention auxiliary, COV schema ready.
2. `GSE92486`: TXT filename but content is COV-like; COV schema ready.
3. `GSE286302`: larger diagnostic COV set; no exact age weeks, diagnostic only.

No training, GroupKFold, autoresearch, or biological-age statistics should be
run from `GSE281602`.


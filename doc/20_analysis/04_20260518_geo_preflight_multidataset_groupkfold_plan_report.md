# GEO Preflight and Multidataset GroupKFold v6 Report

> Date: 2026-05-18  
> Scope: official GEO preflight, processed supplement conversion, unified RRBS 5kb region matrix, dataset-heldout benchmark

## Decision

v6 moved the project from metadata-only supplement assumptions to official GEO preflight plus dataset-aware conversion. The main benchmark now uses four metadata-mapped RRBS datasets: GSE120137, GSE80672, GSE93957, and GSE121141. GSE60012 was converted but excluded from the benchmark because the official tile matrix columns are not unique GSM sample accessions.

Official GEO references used:

- GEO Download: https://www.ncbi.nlm.nih.gov/geo/info/download.html
- GEO Programmatic Access: https://www.ncbi.nlm.nih.gov/geo/info/geo_paccess.html
- GEO SOFT format: https://www.ncbi.nlm.nih.gov/geo/info/soft.html

## Implemented Artifacts

- Added `scripts/etl/11_geo_processed_preflight.py` for official FTP/filelist/HEAD preflight without large downloads.
- Extended `scripts/etl/08_geo_supplement_inventory.py` for dataset-aware supplement priority and schema guesses.
- Refactored `scripts/etl/10_convert_processed_methylation.py` to support GSE80672 overlap tables, Bismark coverage tar archives, and GSE60012 100bp tile matrices.
- Added `scripts/etl/12_build_multidataset_region_matrix.py` for unified inner-join 5kb region matrices.
- Extended `scripts/train/train_heldout_clock.py` with `--train_datasets all_except:<GSE>` and multi-dataset selectors.
- Updated `scripts/train/train_clock.py` dataset scope reporting for true multidataset GroupKFold.

## Data and Schema Results

| Dataset | Supplement | Schema | Status | Metadata overlap | Regions |
|---|---|---|---|---:|---:|
| GSE120137 | existing phase0 matrix | 5kb region | included | 549 | 77,826 |
| GSE80672 | `GSE80672_RAW.tar` | overlap percentage/coverage | included | 255 | 128,739 |
| GSE93957 | `GSE93957_RAW.tar` | Bismark coverage tar | included | 62 | 171,918 |
| GSE121141 | `GSE121141_RAW.tar` | Bismark coverage tar | included | 81 | 146,979 |
| GSE60012 | `GSE60012_100bpTiles_RRBS_Mouse.txt.gz` | 100bp tile matrix | excluded | 0 | 109,777 |
| GSE52266 | percent/count matrices | low-priority adapter | inventory only | N/A | N/A |

Unified matrix:

| Artifact | Value |
|---|---:|
| Matrix | `results/multidataset/all_rrbs_region_matrix_5kb.parquet` |
| Join mode | inner region intersection |
| Included datasets | GSE120137, GSE80672, GSE93957, GSE121141 |
| Samples | 947 |
| Common 5kb regions | 72,079 |
| Intervention labels | 915 control, 32 CR |

GSE60012 blocker:

- Matrix columns are phenotype labels rather than unique GSM IDs.
- Many labels are duplicated, so mapping by order would be unsafe.
- The data were converted and documented but excluded from training.

## Benchmark Results

GroupKFold by `dataset_batch`, model `lgbm`, 500 train-fold-selected regions, median imputation:

| Benchmark | Pearson r | MAE weeks | R2 | Cross-dataset MAE |
|---|---:|---:|---:|---:|
| v3 region best, GSE120137 CV | 0.8530 | 11.550 | 0.7051 | N/A |
| v4 embedding POC best, GSE120137 CV | 0.8594 | 11.106 | 0.7221 | N/A |
| v5 GSE120137 -> GSE80672 | 0.7295 | 33.408 | 0.3438 | 33.408 |
| v6 GroupKFold, 4 datasets | 0.5812 | 24.351 | 0.3096 | 26.012 |
| v6 all-except-GSE80672 -> GSE80672 | 0.7322 | 28.441 | 0.4331 | 28.441 |

LODO fold metrics from the GroupKFold predictions:

| Held-out dataset | n | Pearson r | MAE weeks | R2 |
|---|---:|---:|---:|---:|
| GSE120137 | 549 | 0.4987 | 21.377 | 0.2115 |
| GSE80672 | 255 | 0.7322 | 28.441 | 0.4331 |
| GSE93957 | 62 | 0.8064 | 16.883 | -0.6928 |
| GSE121141 | 81 | 0.3647 | 37.347 | -0.5673 |

Random-label sanity check:

| Check | Pearson r | MAE weeks | R2 |
|---|---:|---:|---:|
| v6 randomized labels | -0.0610 | 35.344 | -0.4799 |

The random-label collapse supports that the v6 GroupKFold path is not obviously leaking labels.

## Intervention Validation

CR validation on GSE80672 uses real held-out predictions only.

| Run | CR AUC | CR F1 | Cohen's d | Mann-Whitney p |
|---|---:|---:|---:|---:|
| v6 all-except-GSE80672 -> GSE80672 | 0.6616 | 0.3019 | 0.4813 | 0.001569 |
| shuffled intervention labels | 0.4865 | 0.2013 | -0.0155 | 0.597672 |

The shuffled-label collapse is good, but the true CR effect remains below the preset threshold. This is a weak research signal, not a biological-age claim.

## Scientific Judgment

v6 improves the data architecture substantially, but it also shows that the current 5kb region clock is still not robust enough across processed RRBS datasets. The next step should be schema/data harmonization, not deep learning:

1. Audit age/tissue/strain composition and methylation value distributions per dataset to find the largest batch shifts.
2. Try conservative normalization variants fitted inside train folds only, such as dataset-wise coverage/presence filters and tissue-stratified reporting.
3. Resolve GSE60012 sample mapping only if an official mapping from tile columns to GSM accessions can be established; otherwise keep it excluded.
4. Keep embedding as interpretation layer after a stable multidataset baseline exists.

# v21 Raw Pilot Execution Report

Date: 2026-05-21

## Scope

This run executed the approved v21 raw pilot path on local RAID FASTQ. It did not run downloads, local deep learning, AutoDL jobs, or broad autoresearch.

## Environment

- `uv` and `.venv-core` Python imports were validated.
- Bismark/Bowtie2/Samtools and the GRCm38/mm10 Bismark index were available.
- Raw FASTQ inventory found 1,457 FASTQ files across `GSE120137`, `GSE121141`, `GSE60012`, `GSE80672`, and `GSE93957`.
- The ETL runner used local Bismark with guarded watchdog execution.

## Raw ETL Completed

### GSE121141 pilot

Queue: `results/ralph_v21_raid_raw_clock/raw_etl_queue.csv`

| sample_id | run_accessions | status | n_regions | rows_total | rows_pass_coverage | beta_min | beta_max |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| GSM3426754 | SRR8032927;SRR8032928 | completed | 87046 | 3344280 | 1760944 | 0.0 | 1.0 |
| GSM3426755 | SRR8032929;SRR8032930 | completed | 109993 | 3800509 | 2602574 | 0.0 | 1.0 |
| GSM3426757 | SRR8032932;SRR8032933 | completed | 131821 | 6499145 | 4400399 | 0.0 | 1.0 |

Matrix: `/data/mouse_methyl/processed_v21_raw_etl/GSE121141/GSE121141_raw_region_matrix_5kb.parquet`

### GSE80672 pilot

Queue: `results/ralph_v21_raid_raw_clock/raw_etl_queue_gse80672_pilot.csv`

| sample_id | title | intervention | age_weeks | run_accessions | status | n_regions | rows_total | rows_pass_coverage | beta_min | beta_max |
| --- | --- | --- | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: |
| GSM2132852 | M0806 | control | 34.766 | SRR3440750 | completed | 131645 | 3342011 | 1630770 | 0.0 | 1.0 |
| GSM2132905 | M2404 | control | 104.297 | SRR3440835 | completed | 144913 | 3329037 | 2217133 | 0.0 | 1.0 |
| GSM2132736 | B6_CR2303 | CR | 99.951 | SRR3440596 | completed | 126023 | 3004232 | 1574651 | 0.0 | 1.0 |

Matrix: `/data/mouse_methyl/processed_v21_raw_etl/GSE80672/GSE80672_raw_region_matrix_5kb.parquet`

## Matrix Gate

`results/ralph_v21_raid_raw_clock/matrix_gate_table.csv`

| dataset | n_samples | n_regions | common_regions | metadata_overlap_fraction | age_coverage_fraction | processed_consistency_median_corr | matrix_gate_passed |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| GSE121141 | 3 | 135778 | 58280 | 1.0 | 1.0 | 0.923432 | True |
| GSE80672 | 3 | 150576 | 65793 | 1.0 | 1.0 | 0.914328 | True |

Decision state: `raw_matrix_ready_for_learn_pending_training_approval`

No critical raw pilot strikes were recorded.

## ML Smoke Check

A 6-sample combined common-region pilot matrix was created at:

`results/ralph_v21_raid_raw_clock/ml/raw_pilot_combined/GSE121141_GSE80672_raw_pilot_common_region_matrix_5kb.parquet`

Combined matrix:

- samples: 6
- datasets: `GSE121141`, `GSE80672`
- common 5kb regions between both raw pilots: 122213

Ridge top500 GroupKFold smoke:

| run | n_samples | cv_strategy | pearson_r | mae_weeks | r2 | status |
| --- | ---: | --- | ---: | ---: | ---: | --- |
| real labels | 6 | GroupKFold_dataset_batch | -0.6870 | 59.313 | -2.5246 | path executed, not a valid success metric |
| random labels | 6 | GroupKFold_dataset_batch | 0.8213 | 59.024 | -2.2934 | sanity not interpretable at n=6 |

This ML smoke confirms that the raw-derived matrix can enter the leakage-controlled training path. It does not meet the v21 ML success gate and should not be used for model selection or biological interpretation.

## Current Decision

Raw ETL and matrix gates passed for the first two priority pilots. Formal ML/autoresearch remains blocked by sample count and support coverage, not by ETL.

Next required action is to expand raw ETL with age/condition-balanced samples before any narrow autoresearch or AutoDL packaging. Priority expansion should add enough samples per dataset and age stratum to make GroupKFold, CR held-out checks, and random-label sanity meaningful.

## Key Outputs

- `results/ralph_v21_raid_raw_clock/raw_inventory.csv`
- `results/ralph_v21_raid_raw_clock/sample_run_manifest.csv`
- `results/ralph_v21_raid_raw_clock/pairing_qc_manifest.csv`
- `results/ralph_v21_raid_raw_clock/matrix_gate_table.csv`
- `results/ralph_v21_raid_raw_clock/ralph_decision_state.json`
- `results/ralph_v21_raid_raw_clock/etl_archives/v21_raw_etl_sample_status_gse121141_pilot.csv`
- `results/ralph_v21_raid_raw_clock/etl_archives/v21_raw_etl_sample_status_gse80672_pilot.csv`
- `results/ralph_v21_raid_raw_clock/ml/raw_pilot_combined/combined_matrix_manifest.json`

# v24 Raw ETL Balanced Expansion and Signal Screen Report

Date: 2026-05-21

## Scope

This report records the first post-pilot raw ETL expansion after the v21 RAID
raw gates. The work stayed on the local workstation for FASTQ, Bismark, COV,
and matrix construction. GPU/cloud DL was not started for this raw subset.

## Raw ETL Execution

- raw inventory: `1457` FASTQ files, `6198.314` GiB, read-only under `/data/mouse_methyl/raw`
- full v24 queue: `43` rows (`GSE121141=7`, `GSE80672=24`, `GSE93957=12`)
- effective queue after the GSE93957 stop: `31` rows (`GSE121141=7`, `GSE80672=24`)
- local ETL output root: `/data/mouse_methyl/processed_v21_raw_etl`
- Bismark runner patch: existing BAM reuse now requires a Bismark report plus `samtools quickcheck`
- sample gate patch: samples with fewer than `50000` 5kb regions fail sample QC instead of entering dataset matrices

## Matrix Gates

| dataset | valid samples | regions | common regions | metadata overlap | age coverage | processed median corr | gate |
| --- | --- | --- | --- | --- | --- | --- | --- |
| GSE121141 | 7 | 167244 | 60278 | 1.0 | 1.0 | 0.939409 | pass |
| GSE80672 | 23 | 184647 | 65862 | 1.0 | 1.0 | 0.920575 | pass |

Both promoted matrices have beta values in `[0, 1]`, primary autosome-only
5kb regions, traceable metadata overlap, age coverage, and raw-vs-processed
consistency audits.

`GSE80672/GSM2132900` was excluded by the new sample gate:
`n_regions=680 < min_sample_regions=50000`.

## GSE93957 Stop

GSE93957 was stopped before expansion after four consecutive samples produced
zero usable primary-autosome 5kb regions under the current mm10/Bismark COV
parser:

| sample | rows total | rows pass coverage | usable regions |
| --- | --- | --- | --- |
| GSM2465647 | 969 | 0 | 0 |
| GSM2465650 | 373 | 0 | 0 |
| GSM2465662 | 789 | 0 | 0 |
| GSM2465655 | 495 | 0 | 0 |

This is a re-evaluation condition, not a training input. Next work for GSE93957
should check library type, assembly/schema compatibility, Bismark parameters,
and metadata/run mapping before any more ETL.

## Raw v24 Signal Screen

The two gate-passed raw-derived matrices were combined on complete common 5kb
regions:

- samples: `30` (`GSE121141=7`, `GSE80672=23`)
- complete common regions: `58215`
- interventions: `CR=12`, `control=18`
- tissues: `blood=22`, `brain_cortex=7`, `fibroblast=1`

Leakage-controlled local ML was run with fold-internal feature selection,
`500` regions, `0.9` train-fold feature presence, robust scaling, and
random-label sanity checks.

| target | CV | model | main result | random-label check | interpretation |
| --- | --- | --- | --- | --- | --- |
| age_weeks_raw_v24 | GroupKFold by dataset | ridge | MAE `43.603w`, r `0.2921` | r `-0.1953`, MAE `42.583w` | not a useful raw age clock yet |
| cr_vs_control_raw_v24 | StratifiedKFold within GSE80672 | logistic | AUC `0.8636`, balanced accuracy `0.822` | AUC `0.4848`, balanced accuracy `0.4394` | promising CR-associated methylation signal |
| dataset_batch_raw_v24 | StratifiedKFold | logistic | AUC `1.0` | AUC `0.0` | strong batch/dataset separability |
| tissue_blood_vs_brain_raw_v24 | StratifiedKFold | logistic | AUC `1.0` | AUC `0.0` | strong tissue/dataset confounding probe |

The CR screen is the biologically meaningful result from this raw subset:
held-out predictions beat random labels and the target confounding audit is
low. It should still be reported as a CR-associated methylation signal, not a
clinical or mechanistic claim. The age result is not strong enough to justify
DL expansion on this raw subset.

## DL Decision

DL optimization has already been run on the larger processed matrix track in
v23. For the raw-derived v24 subset, DL should wait until raw ETL expands to a
substantially larger sample count. The current raw matrix is valuable for ETL
validation and CR signal discovery, but `n=30` is not enough for robust neural
network optimization.

## Decision State

- raw matrix state: `raw_matrix_ready_for_learn_pending_training_approval`
- local ML on gate-passed raw matrices: completed for an exploratory v24 signal screen
- broad autoresearch: not started for raw v24
- raw DL: not started; wait for more raw samples
- GSE93957: stop and re-evaluate

## Key Outputs

- `results/ralph_v21_raid_raw_clock/raw_etl_queue_v24_balanced_expansion.csv`
- `results/ralph_v21_raid_raw_clock/raw_etl_queue_v24_effective_gse121141_gse80672.csv`
- `results/ralph_v21_raid_raw_clock/matrix_gate_table.csv`
- `/data/mouse_methyl/processed_v21_raw_etl/v24_full_queue_sample_status_before_gse93957_stop.csv`
- `/data/mouse_methyl/processed_v21_raw_etl/v21_raw_etl_sample_status.csv`
- `results/v24_raw_signal_screen/raw_v24_signal_input_summary.json`
- `results/v24_raw_signal_screen/benchmark_summary.csv`
- `results/v24_raw_signal_screen/random_label_sanity_summary.csv`
- `results/v24_raw_signal_screen/target_confounding_audit.csv`
- `results/v24_raw_signal_screen/selected_region_signal_audit.csv`
- `results/v24_raw_signal_screen/v22_biological_signal_screen_report.md`

## Next Work

1. Expand raw ETL within the already-passed families, prioritizing more
   GSE80672 CR/control and more GSE121141 ages/tissues.
2. Re-evaluate GSE93957 before touching more of its FASTQ files.
3. Rerun raw signal ML after the next expansion and require CR random-label
   sanity to remain near random.
4. Only package raw-derived DL for AutoDL after raw sample count is large enough
   for validation beyond small-n memorization.

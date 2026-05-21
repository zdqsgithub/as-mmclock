# v21 RAID Raw RALPH Report

Date: 2026-05-21T04:35:41.547092+00:00

## Scope

Read-only inventory and gate preflight for local RAID FASTQ. This run did not start Bismark, downloads, training, or autoresearch.

## Raw Inventory

- FASTQ files: `1457`
- FASTQ size GiB: `6198.314`
- datasets with FASTQ: `GSE120137, GSE121141, GSE60012, GSE80672, GSE93957`
- sample/run rows: `1019`
- complete paired runs: `438`
- pilot queue rows: `12`

## Environment Gate

- passed: `True`
- CT index files: `6`
- GA index files: `6`

## Matrix Gate

| dataset | status | n_samples | n_regions | common_regions | metadata_overlap_fraction | age_coverage_fraction | beta_range_valid | sex_mt_excluded | matrix_gate_passed | decision |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE121141 | completed | 7 | 167244 | 60278 | 1.0 | 1.0 | True | True | True | raw_matrix_ready_for_learn_pending_training_approval |
| GSE80672 | completed | 23 | 184647 | 65862 | 1.0 | 1.0 | True | True | True | raw_matrix_ready_for_learn_pending_training_approval |

## Decision

- state: `raw_matrix_ready_for_learn_pending_training_approval`
- training authorized: `False`
- autoresearch authorized: `False`
- next action: Review matrix gate and record explicit training approval before local ML.

## Outputs

- `raw_inventory.csv`
- `sample_run_manifest.csv`
- `pairing_qc_manifest.csv`
- `raw_etl_queue.csv`
- `environment_gate.json`
- `matrix_gate_table.csv`
- `watchdog_state.json`
- `ralph_decision_state.json`

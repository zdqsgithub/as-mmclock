# v13 Delivery Freeze Summary

Date: 2026-05-19T09:13:47+08:00

## Decision

- Selected route: `route_c_accept_redefined_benchmark`.
- Current model headline is support-covered chronological-age prediction.
- `GSE121141 old104+ brain_cortex/heart/lung` remains a stress-test/blocker metric.
- This freeze does not authorize training, downloads, FASTQ/Bismark, or autoresearch.

## Metrics

- Support-covered headline MAE: `20.286` weeks over `15573` rows.
- Unsupported stress-test MAE: `36.659` weeks over `7620` rows.
- GSE121141 old104+ stress MAE: `88.565` weeks over `420` rows.
- Best observed CR AUC: `0.8744`.

## Reproduction Outputs

- `results/v13_delivery_freeze/environment_manifest.json`
- `results/v13_delivery_freeze/artifact_inventory.csv`
- `results/v13_delivery_freeze/reproducibility_manifest.json`
- `results/v13_delivery_freeze/v13_delivery_summary.md`
- `results/v13_delivery_freeze/test_report.json`
- `results/v13_delivery_freeze/test_report.txt`
- `doc/00_meta/01_20260519_v13_model_card.md`
- `doc/00_meta/02_20260519_project_status_index_v13.md`
- `doc/30_protocols/04_20260519_route_a_old_tissue_data_generation_rfc_template.md`
- `doc/30_protocols/05_20260519_route_b_minimal_fastq_bismark_pilot_rfc_template.md`

## Dependency Mode

- Default delivery reproduction uses `requirements-core.txt`.
- `requirements-full.txt` is reserved for training/deep-learning/full historical environments.

## Guardrails

- raw_fastq_download_authorized: `False`
- training_authorized: `False`
- autoresearch_authorized: `False`

## Delivery Tests

- status: `passed`
- tests_run: `6`
- failures: `0`
- errors: `0`

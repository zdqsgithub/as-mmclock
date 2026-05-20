# Route A Intake Gate Readiness Report

Date: 2026-05-19T09:58:47

## Summary

This report records the latest Route A intake gate workflow state.
It does not authorize download, Bismark, model training, or autoresearch.

## Decision

- workflow decision: `ready_for_ralph_learn_pending_explicit_training_approval`
- metadata gate: `passed`
- adapter smoke: `passed`
- matrix gate: `passed`
- RALPH Learn readiness: `ready_pending_explicit_training_approval`

## Gate Details

- metadata errors: `0`; warnings: `0`
- adapter reason: `at_least_one_processed_file_passed_schema_smoke`
- adapter files passed: `3` / `3`
- matrix reason/issues: ``
- matrix metrics: `{'manifest_status': 'completed', 'common_regions': 50010.0, 'metadata_overlap': 1.0, 'age_coverage': 1.0, 'beta_min': 0.050000112503767014, 'beta_max': 0.9499999284744263}`
- learn reason: `metadata_adapter_matrix_gates_passed`

## Outputs

- gate output directory: `/home/zdq-as/mouse_methyl_work/results/route_a_contract_tests/generated_manifest_intake`
- `metadata_gate_report.json`
- `adapter_smoke_report.json`
- `matrix_gate_report.json`
- `ralph_learn_readiness.json`
- `route_a_gate_state.json`

## Guardrails

- download_authorized: `False`
- bismark_authorized: `False`
- training_authorized: `False`
- autoresearch_authorized: `False`

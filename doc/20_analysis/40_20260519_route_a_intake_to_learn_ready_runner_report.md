# Route A Intake To Learn Ready Runner Report

Date: 2026-05-19T09:40:53

## Summary

This runner executes the Route A gate sequence up to RALPH Learn readiness.
It does not train, download, run Bismark, or start autoresearch.

## State

- submission_id: `synthetic_ready`
- status: `ready_for_ralph_learn_pending_explicit_training_approval`
- intake_dir: `/home/zdq-as/mouse_methyl_work/results/route_a_intake/synthetic_ready`
- package_dir: `/home/zdq-as/mouse_methyl_work/results/route_a_ralph_learn_package/synthetic_ready`
- gate decision: `ready_for_ralph_learn_pending_explicit_training_approval`
- matrix status: `completed`
- package status: `prepared_pending_explicit_training_approval`

## Guardrails

- training_authorized: `False`
- download_authorized: `False`
- bismark_authorized: `False`
- autoresearch_authorized: `False`

## Stage Return Codes

- `01_metadata_adapter_gate`: returncode `0`
- `02_build_matrix`: returncode `0`
- `03_matrix_gate_readiness`: returncode `0`
- `04_prepare_guarded_learn_package`: returncode `0`

## Next Step

If status is `ready_for_ralph_learn_pending_explicit_training_approval`,
review the guarded command manifest and request explicit training approval
before running any benchmark command.

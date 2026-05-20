# Route A RALPH Learn Package Report

Date: 2026-05-19T09:40:53

## Summary

This package prepares fixed Route A benchmark inputs and commands after the intake gate reaches readiness.
No training, download, Bismark, or autoresearch was run.

## Package State

- status: `prepared_pending_explicit_training_approval`
- submission_id: `synthetic_ready`
- combined matrix: `/home/zdq-as/mouse_methyl_work/results/route_a_ralph_learn_package/synthetic_ready/route_a_plus_v8_2_region_matrix_5kb.parquet`
- model metadata: `/home/zdq-as/mouse_methyl_work/results/route_a_ralph_learn_package/synthetic_ready/route_a_model_metadata.csv`
- command count: `4`

## Commands

Commands are written to `ralph_learn_command_manifest.csv/json` and `ralph_learn_commands.sh`.
The shell script exits before executing commands until training is explicitly authorized.

## Guardrails

- training_authorized: `False`
- autoresearch_authorized: `False`
- download_authorized: `False`
- bismark_authorized: `False`

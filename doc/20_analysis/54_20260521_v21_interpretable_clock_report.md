# v21 Interpretable Clock Report

Date: 2026-05-20T16:24:31.815089+00:00

## Scope

Region-level interpretation of selected v19 ML and v20 DL clock features. The labels describe age-associated methylation signals and confounding risk; they do not establish aging mechanisms.

## Summary

- top regions interpreted: `100`
- cross-tissue stable: `86`
- tissue-specific: `1`
- dataset-confounded: `0`
- coverage-driven: `0`

## Important Caveat

Permutation, DL gradient, and DL occlusion columns are deterministic interpretation proxies unless a downstream run supplies model-native importance tensors. They are adequate for triage and reporting guardrails, not mechanistic claims.

## Outputs

- `/home/zdq-as/mouse_methyl_work/results/ralph_v21_raid_raw_clock/feature_interpretation/top_region_interpretation.csv`
- `/home/zdq-as/mouse_methyl_work/results/ralph_v21_raid_raw_clock/feature_interpretation/region_cluster_summary.csv`
- `/home/zdq-as/mouse_methyl_work/results/ralph_v21_raid_raw_clock/feature_interpretation/feature_confounding_audit.csv`
- `/home/zdq-as/mouse_methyl_work/results/ralph_v21_raid_raw_clock/feature_interpretation/v21_feature_interpretation_summary.json`

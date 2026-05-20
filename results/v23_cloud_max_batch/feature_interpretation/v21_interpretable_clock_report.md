# v21 Interpretable Clock Report

Date: 2026-05-20T22:27:43.025933+00:00

## Scope

Region-level interpretation of selected v19 ML and v20 DL clock features. The labels describe age-associated methylation signals and confounding risk; they do not establish aging mechanisms.

## Summary

- top regions interpreted: `250`
- cross-tissue stable: `180`
- tissue-specific: `2`
- dataset-confounded: `0`
- coverage-driven: `0`

## Important Caveat

Permutation, DL gradient, and DL occlusion columns are deterministic interpretation proxies unless a downstream run supplies model-native importance tensors. They are adequate for triage and reporting guardrails, not mechanistic claims.

## Outputs

- `results/v23_cloud_max_batch/feature_interpretation/top_region_interpretation.csv`
- `results/v23_cloud_max_batch/feature_interpretation/region_cluster_summary.csv`
- `results/v23_cloud_max_batch/feature_interpretation/feature_confounding_audit.csv`
- `results/v23_cloud_max_batch/feature_interpretation/v21_feature_interpretation_summary.json`

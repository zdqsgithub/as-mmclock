# v22 Biological Signal Screen Report

Date: 2026-05-21T04:42:34.742308+00:00

## Scope

This report implements the v22 pivot from age-only clocks to biological-signal discovery using local mouse methylation resources. Age remains the anchor task; CR, sex, tissue, and condition labels are screened as auxiliary biological signals with explicit confounding audits.

## Inputs

- matrix: `results/v24_raw_signal_screen/raw_v24_gse121141_gse80672_common_complete_region_matrix_5kb.parquet`
- metadata: `results/v24_raw_signal_screen/raw_v24_signal_metadata.csv`
- target registry: `results/v24_raw_signal_screen/raw_v24_biological_signal_targets.csv`
- aligned samples: `30`
- regions: `58215`

## Target Gates

- passed targets: `4`
- skipped targets: `0`

Passed targets:

- `age_weeks_raw_v24`: n=30, status=passed, reason=target passed sample and group gates
- `cr_vs_control_raw_v24`: n=23, status=passed, reason=target passed sample and class gates
- `dataset_batch_raw_v24`: n=30, status=passed, reason=target passed sample and class gates
- `tissue_blood_vs_brain_raw_v24`: n=29, status=passed, reason=target passed sample and class gates

## Best Local ML Results

- `age_weeks_raw_v24`, model=ridge, cv=GroupKFold_dataset_batch, mae_weeks=43.603, pearson_r=0.2921
- `cr_vs_control_raw_v24`, model=logistic, cv=StratifiedKFold_within_scope, balanced_accuracy=0.822, roc_auc=0.8636, macro_f1=0.8231
- `dataset_batch_raw_v24`, model=logistic, cv=StratifiedKFold_within_scope, balanced_accuracy=1.0, roc_auc=1.0, macro_f1=1.0
- `tissue_blood_vs_brain_raw_v24`, model=logistic, cv=StratifiedKFold_within_scope, balanced_accuracy=1.0, roc_auc=1.0, macro_f1=1.0

## Interpretation Guardrails

Feature tables report target-associated methylation signals only. Regions marked dataset-, tissue-, or coverage-confounded must not be described as mechanisms without follow-up validation.

## Outputs

- `results/v24_raw_signal_screen/target_registry.csv`
- `results/v24_raw_signal_screen/target_gate_table.csv`
- `results/v24_raw_signal_screen/benchmark_summary.csv`
- `results/v24_raw_signal_screen/random_label_sanity_summary.csv`
- `results/v24_raw_signal_screen/target_confounding_audit.csv`
- `results/v24_raw_signal_screen/selected_region_signal_audit.csv`
- `results/v24_raw_signal_screen/autodl_signal_candidate_manifest.csv`
- `results/v24_raw_signal_screen/v22_biological_signal_screen_summary.json`

## External Anchors

- GSE80672/Petkovich: https://pubmed.ncbi.nlm.nih.gov/28380383/
- GSE121141/Meer whole-lifespan clock: https://elifesciences.org/articles/40675
- Region-based RRBS clock: https://www.pure.ed.ac.uk/ws/portalfiles/portal/347454803/Aging_Cell_2023_Simpson_Region_based_epigenetic_clock_design_improves_RRBS_based_age_prediction.pdf
- scEpiAge code: https://github.com/EpigenomeClock/scEpiAge

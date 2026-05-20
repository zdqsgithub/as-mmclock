# v22 Biological Signal Screen Report

Date: 2026-05-20T19:49:50.360089+00:00

## Scope

This report implements the v22 pivot from age-only clocks to biological-signal discovery using local mouse methylation resources. Age remains the anchor task; CR, sex, tissue, and condition labels are screened as auxiliary biological signals with explicit confounding audits.

## Inputs

- matrix: `/home/zdq-as/mouse_methyl_work/results/multidataset_v8_3_ablation/all6/all_rrbs_region_matrix_5kb.parquet`
- metadata: `/home/zdq-as/mouse_methyl_work/metadata/model_sample_metadata_v8.csv`
- target registry: `/home/zdq-as/mouse_methyl_work/metadata/v22_biological_signal_targets.csv`
- aligned samples: `1219`
- regions: `65870`

## Target Gates

- passed targets: `6`
- skipped targets: `1`

Passed targets:

- `age_weeks`: n=1219, status=passed, reason=target passed sample and group gates
- `sex_binary`: n=1037, status=passed, reason=target passed sample and class gates
- `tissue_multiclass`: n=1193, status=passed, reason=target passed sample and class gates
- `cr_vs_control`: n=255, status=passed, reason=target passed sample and class gates
- `old_castration_vs_control`: n=20, status=passed, reason=target passed sample and class gates
- `condition_family_multiclass`: n=141, status=passed, reason=target passed sample and class gates

Skipped targets:

- `gse52266_hfd_vs_control`: n_samples=0 < min_total=30

## Best Local ML Results

- `age_weeks`, model=ridge, cv=GroupKFold_dataset_batch, mae_weeks=33.214, pearson_r=0.4217
- `condition_family_multiclass`, model=logistic, cv=StratifiedKFold_within_scope, balanced_accuracy=0.5794, macro_f1=0.5547
- `cr_vs_control`, model=logistic, cv=StratifiedKFold_within_scope, balanced_accuracy=0.9308, roc_auc=0.9937, macro_f1=0.9366
- `old_castration_vs_control`, model=logistic, cv=StratifiedKFold_within_scope, balanced_accuracy=0.65, roc_auc=0.72, macro_f1=0.6419
- `sex_binary`, model=logistic, cv=GroupKFold_dataset_batch_class_covered, balanced_accuracy=0.7181, roc_auc=0.7667, macro_f1=0.6899
- `tissue_multiclass`, model=logistic, cv=StratifiedKFold_within_scope, balanced_accuracy=0.9927, macro_f1=0.9815

## Interpretation Guardrails

Feature tables report target-associated methylation signals only. Regions marked dataset-, tissue-, or coverage-confounded must not be described as mechanisms without follow-up validation.

## Outputs

- `results/v22_biological_signal_screen/target_registry.csv`
- `results/v22_biological_signal_screen/target_gate_table.csv`
- `results/v22_biological_signal_screen/benchmark_summary.csv`
- `results/v22_biological_signal_screen/random_label_sanity_summary.csv`
- `results/v22_biological_signal_screen/target_confounding_audit.csv`
- `results/v22_biological_signal_screen/selected_region_signal_audit.csv`
- `results/v22_biological_signal_screen/autodl_signal_candidate_manifest.csv`
- `results/v22_biological_signal_screen/v22_biological_signal_screen_summary.json`

## External Anchors

- GSE80672/Petkovich: https://pubmed.ncbi.nlm.nih.gov/28380383/
- GSE121141/Meer whole-lifespan clock: https://elifesciences.org/articles/40675
- Region-based RRBS clock: https://www.pure.ed.ac.uk/ws/portalfiles/portal/347454803/Aging_Cell_2023_Simpson_Region_based_epigenetic_clock_design_improves_RRBS_based_age_prediction.pdf
- scEpiAge code: https://github.com/EpigenomeClock/scEpiAge

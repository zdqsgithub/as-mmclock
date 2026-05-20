# v22 Biological Signal Screen Report

Date: 2026-05-20T21:34:40.440597+00:00

## Scope

This report implements the v22 pivot from age-only clocks to biological-signal discovery using local mouse methylation resources. Age remains the anchor task; CR, sex, tissue, and condition labels are screened as auxiliary biological signals with explicit confounding audits.

## Inputs

- matrix: `/root/autodl-tmp/mouse_methyl_work/results/multidataset_v8_3_ablation/all6/all_rrbs_region_matrix_5kb.parquet`
- metadata: `/root/autodl-tmp/mouse_methyl_work/metadata/model_sample_metadata_v8.csv`
- target registry: `/root/autodl-tmp/mouse_methyl_work/metadata/v22_biological_signal_targets.csv`
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

- `age_weeks`, model=ridge, cv=GroupKFold_dataset_batch, mae_weeks=33.286, pearson_r=0.4162
- `condition_family_multiclass`, model=logistic, cv=StratifiedKFold_within_scope, balanced_accuracy=0.3913, macro_f1=0.3922
- `cr_vs_control`, model=logistic, cv=StratifiedKFold_within_scope, balanced_accuracy=0.9464, roc_auc=0.9851, macro_f1=0.9464
- `old_castration_vs_control`, model=logistic, cv=StratifiedKFold_within_scope, balanced_accuracy=0.7, roc_auc=0.68, macro_f1=0.697
- `sex_binary`, model=logistic, cv=GroupKFold_dataset_batch_class_covered, balanced_accuracy=0.6963, roc_auc=0.7591, macro_f1=0.6627
- `tissue_multiclass`, model=logistic, cv=StratifiedKFold_within_scope, balanced_accuracy=0.9948, macro_f1=0.9885

## Interpretation Guardrails

Feature tables report target-associated methylation signals only. Regions marked dataset-, tissue-, or coverage-confounded must not be described as mechanisms without follow-up validation.

## Outputs

- `results/v23_cloud_max_batch/biology_sidecars/auto_robust2000/target_registry.csv`
- `results/v23_cloud_max_batch/biology_sidecars/auto_robust2000/target_gate_table.csv`
- `results/v23_cloud_max_batch/biology_sidecars/auto_robust2000/benchmark_summary.csv`
- `results/v23_cloud_max_batch/biology_sidecars/auto_robust2000/random_label_sanity_summary.csv`
- `results/v23_cloud_max_batch/biology_sidecars/auto_robust2000/target_confounding_audit.csv`
- `results/v23_cloud_max_batch/biology_sidecars/auto_robust2000/selected_region_signal_audit.csv`
- `results/v23_cloud_max_batch/biology_sidecars/auto_robust2000/autodl_signal_candidate_manifest.csv`
- `results/v23_cloud_max_batch/biology_sidecars/auto_robust2000/v22_biological_signal_screen_summary.json`

## External Anchors

- GSE80672/Petkovich: https://pubmed.ncbi.nlm.nih.gov/28380383/
- GSE121141/Meer whole-lifespan clock: https://elifesciences.org/articles/40675
- Region-based RRBS clock: https://www.pure.ed.ac.uk/ws/portalfiles/portal/347454803/Aging_Cell_2023_Simpson_Region_based_epigenetic_clock_design_improves_RRBS_based_age_prediction.pdf
- scEpiAge code: https://github.com/EpigenomeClock/scEpiAge

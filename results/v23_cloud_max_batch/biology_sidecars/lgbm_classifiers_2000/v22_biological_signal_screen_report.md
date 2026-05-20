# v22 Biological Signal Screen Report

Date: 2026-05-20T21:39:17.175650+00:00

## Scope

This report implements the v22 pivot from age-only clocks to biological-signal discovery using local mouse methylation resources. Age remains the anchor task; CR, sex, tissue, and condition labels are screened as auxiliary biological signals with explicit confounding audits.

## Inputs

- matrix: `/root/autodl-tmp/mouse_methyl_work/results/multidataset_v8_3_ablation/all6/all_rrbs_region_matrix_5kb.parquet`
- metadata: `/root/autodl-tmp/mouse_methyl_work/metadata/model_sample_metadata_v8.csv`
- target registry: `/root/autodl-tmp/mouse_methyl_work/metadata/v22_biological_signal_targets.csv`
- aligned samples: `1219`
- regions: `65870`

## Target Gates

- passed targets: `5`
- skipped targets: `0`

Passed targets:

- `sex_binary`: n=1037, status=passed, reason=target passed sample and class gates
- `tissue_multiclass`: n=1193, status=passed, reason=target passed sample and class gates
- `cr_vs_control`: n=255, status=passed, reason=target passed sample and class gates
- `old_castration_vs_control`: n=20, status=passed, reason=target passed sample and class gates
- `condition_family_multiclass`: n=141, status=passed, reason=target passed sample and class gates

## Best Local ML Results

- `condition_family_multiclass`, model=lgbm, cv=StratifiedKFold_within_scope, balanced_accuracy=0.5179, macro_f1=0.5461
- `cr_vs_control`, model=lgbm, cv=StratifiedKFold_within_scope, balanced_accuracy=0.7901, roc_auc=0.9742, macro_f1=0.8343
- `old_castration_vs_control`, model=lgbm, cv=StratifiedKFold_within_scope, balanced_accuracy=0.5, roc_auc=0.5, macro_f1=0.3333
- `sex_binary`, model=lgbm, cv=GroupKFold_dataset_batch_class_covered, balanced_accuracy=0.7204, roc_auc=0.7214, macro_f1=0.6578
- `tissue_multiclass`, model=lgbm, cv=StratifiedKFold_within_scope, balanced_accuracy=0.9817, macro_f1=0.9824

## Interpretation Guardrails

Feature tables report target-associated methylation signals only. Regions marked dataset-, tissue-, or coverage-confounded must not be described as mechanisms without follow-up validation.

## Outputs

- `results/v23_cloud_max_batch/biology_sidecars/lgbm_classifiers_2000/target_registry.csv`
- `results/v23_cloud_max_batch/biology_sidecars/lgbm_classifiers_2000/target_gate_table.csv`
- `results/v23_cloud_max_batch/biology_sidecars/lgbm_classifiers_2000/benchmark_summary.csv`
- `results/v23_cloud_max_batch/biology_sidecars/lgbm_classifiers_2000/random_label_sanity_summary.csv`
- `results/v23_cloud_max_batch/biology_sidecars/lgbm_classifiers_2000/target_confounding_audit.csv`
- `results/v23_cloud_max_batch/biology_sidecars/lgbm_classifiers_2000/selected_region_signal_audit.csv`
- `results/v23_cloud_max_batch/biology_sidecars/lgbm_classifiers_2000/autodl_signal_candidate_manifest.csv`
- `results/v23_cloud_max_batch/biology_sidecars/lgbm_classifiers_2000/v22_biological_signal_screen_summary.json`

## External Anchors

- GSE80672/Petkovich: https://pubmed.ncbi.nlm.nih.gov/28380383/
- GSE121141/Meer whole-lifespan clock: https://elifesciences.org/articles/40675
- Region-based RRBS clock: https://www.pure.ed.ac.uk/ws/portalfiles/portal/347454803/Aging_Cell_2023_Simpson_Region_based_epigenetic_clock_design_improves_RRBS_based_age_prediction.pdf
- scEpiAge code: https://github.com/EpigenomeClock/scEpiAge

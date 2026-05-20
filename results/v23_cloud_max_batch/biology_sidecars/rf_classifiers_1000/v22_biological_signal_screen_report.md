# v22 Biological Signal Screen Report

Date: 2026-05-20T21:40:41.710452+00:00

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

- `condition_family_multiclass`, model=rf, cv=StratifiedKFold_within_scope, balanced_accuracy=0.2143, macro_f1=0.2275
- `cr_vs_control`, model=rf, cv=StratifiedKFold_within_scope, balanced_accuracy=0.7165, roc_auc=0.9762, macro_f1=0.7774
- `old_castration_vs_control`, model=rf, cv=StratifiedKFold_within_scope, balanced_accuracy=0.55, roc_auc=0.55, macro_f1=0.5396
- `sex_binary`, model=rf, cv=GroupKFold_dataset_batch_class_covered, balanced_accuracy=0.7021, roc_auc=0.7563, macro_f1=0.667
- `tissue_multiclass`, model=rf, cv=StratifiedKFold_within_scope, balanced_accuracy=0.9867, macro_f1=0.9856

## Interpretation Guardrails

Feature tables report target-associated methylation signals only. Regions marked dataset-, tissue-, or coverage-confounded must not be described as mechanisms without follow-up validation.

## Outputs

- `results/v23_cloud_max_batch/biology_sidecars/rf_classifiers_1000/target_registry.csv`
- `results/v23_cloud_max_batch/biology_sidecars/rf_classifiers_1000/target_gate_table.csv`
- `results/v23_cloud_max_batch/biology_sidecars/rf_classifiers_1000/benchmark_summary.csv`
- `results/v23_cloud_max_batch/biology_sidecars/rf_classifiers_1000/random_label_sanity_summary.csv`
- `results/v23_cloud_max_batch/biology_sidecars/rf_classifiers_1000/target_confounding_audit.csv`
- `results/v23_cloud_max_batch/biology_sidecars/rf_classifiers_1000/selected_region_signal_audit.csv`
- `results/v23_cloud_max_batch/biology_sidecars/rf_classifiers_1000/autodl_signal_candidate_manifest.csv`
- `results/v23_cloud_max_batch/biology_sidecars/rf_classifiers_1000/v22_biological_signal_screen_summary.json`

## External Anchors

- GSE80672/Petkovich: https://pubmed.ncbi.nlm.nih.gov/28380383/
- GSE121141/Meer whole-lifespan clock: https://elifesciences.org/articles/40675
- Region-based RRBS clock: https://www.pure.ed.ac.uk/ws/portalfiles/portal/347454803/Aging_Cell_2023_Simpson_Region_based_epigenetic_clock_design_improves_RRBS_based_age_prediction.pdf
- scEpiAge code: https://github.com/EpigenomeClock/scEpiAge

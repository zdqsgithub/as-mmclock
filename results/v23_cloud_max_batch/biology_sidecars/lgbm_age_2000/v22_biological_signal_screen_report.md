# v22 Biological Signal Screen Report

Date: 2026-05-20T21:35:30.724781+00:00

## Scope

This report implements the v22 pivot from age-only clocks to biological-signal discovery using local mouse methylation resources. Age remains the anchor task; CR, sex, tissue, and condition labels are screened as auxiliary biological signals with explicit confounding audits.

## Inputs

- matrix: `/root/autodl-tmp/mouse_methyl_work/results/multidataset_v8_3_ablation/all6/all_rrbs_region_matrix_5kb.parquet`
- metadata: `/root/autodl-tmp/mouse_methyl_work/metadata/model_sample_metadata_v8.csv`
- target registry: `/root/autodl-tmp/mouse_methyl_work/metadata/v22_biological_signal_targets.csv`
- aligned samples: `1219`
- regions: `65870`

## Target Gates

- passed targets: `1`
- skipped targets: `0`

Passed targets:

- `age_weeks`: n=1219, status=passed, reason=target passed sample and group gates

## Best Local ML Results

- `age_weeks`, model=lgbm, cv=GroupKFold_dataset_batch, mae_weeks=25.115, pearson_r=0.5704

## Interpretation Guardrails

Feature tables report target-associated methylation signals only. Regions marked dataset-, tissue-, or coverage-confounded must not be described as mechanisms without follow-up validation.

## Outputs

- `results/v23_cloud_max_batch/biology_sidecars/lgbm_age_2000/target_registry.csv`
- `results/v23_cloud_max_batch/biology_sidecars/lgbm_age_2000/target_gate_table.csv`
- `results/v23_cloud_max_batch/biology_sidecars/lgbm_age_2000/benchmark_summary.csv`
- `results/v23_cloud_max_batch/biology_sidecars/lgbm_age_2000/random_label_sanity_summary.csv`
- `results/v23_cloud_max_batch/biology_sidecars/lgbm_age_2000/target_confounding_audit.csv`
- `results/v23_cloud_max_batch/biology_sidecars/lgbm_age_2000/selected_region_signal_audit.csv`
- `results/v23_cloud_max_batch/biology_sidecars/lgbm_age_2000/autodl_signal_candidate_manifest.csv`
- `results/v23_cloud_max_batch/biology_sidecars/lgbm_age_2000/v22_biological_signal_screen_summary.json`

## External Anchors

- GSE80672/Petkovich: https://pubmed.ncbi.nlm.nih.gov/28380383/
- GSE121141/Meer whole-lifespan clock: https://elifesciences.org/articles/40675
- Region-based RRBS clock: https://www.pure.ed.ac.uk/ws/portalfiles/portal/347454803/Aging_Cell_2023_Simpson_Region_based_epigenetic_clock_design_improves_RRBS_based_age_prediction.pdf
- scEpiAge code: https://github.com/EpigenomeClock/scEpiAge

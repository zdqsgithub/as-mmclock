# v12.1 Redefined Benchmark Report

Date: 2026-05-19

## Summary

v12.1 applied the v12 benchmark redefinition to existing prediction files only. It did not retrain models, download data, rebuild matrices, run FASTQ/Bismark, or start autoresearch.

- Prediction files scanned: 53
- Benchmark redefinition accepted: `true`
- Support-covered headline MAE: `20.286` weeks
- Unsupported stress-test MAE: `36.659` weeks

## Aggregate Metrics

| level | source_prediction | subset | n_samples | n_datasets | pearson_r | pearson_pval | mae_weeks | medae_weeks | rmse_weeks | bias_weeks | r2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| aggregate | all | all_predictions | 23193 | 6 | 0.6161 | 0.0 | 25.666 | 19.992 | 33.411 | 14.344 | 0.239 |
| aggregate | all | support_covered_headline | 15573 | 6 | 0.5441 | 0.0 | 20.286 | 15.389 | 26.356 | 7.287 | 0.2259 |
| aggregate | all | unsupported_stress | 7620 | 5 | 0.6435 | 0.0 | 36.659 | 33.596 | 44.475 | 28.766 | -0.0084 |
| aggregate | all | gse121141_old104_stress | 420 | 1 |  |  | 88.565 | 94.452 | 90.776 | 88.565 |  |

## GSE121141 old104+ Stress Metric

| level | source_prediction | subset | n_samples | n_datasets | pearson_r | pearson_pval | mae_weeks | medae_weeks | rmse_weeks | bias_weeks | r2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| aggregate | all | gse121141_old104_stress | 420 | 1 |  |  | 88.565 | 94.452 | 90.776 | 88.565 |  |

## GSE80672 CR Metrics

| source_prediction | n_samples | support_covered_n | cr_detection_auc | cr_detection_f1 | cr_cohens_d | cr_mannwhitney_p | cr_status | shuffled_sanity_status | shuffled_cr_detection_auc | shuffled_cr_detection_f1 | shuffled_cr_cohens_d | shuffled_cr_mannwhitney_p | shuffled_cr_status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| results/benchmark_v8_3_ablation/core4_plus_gse213628/02_v74_lgbm_robust_p095_top1000_shift015/groupkfold/predictions.csv | 255 | 178 | 0.8744 | 0.3529 | 1.5393 | 0.0 | computed | not_available_not_recomputed |  |  |  |  |  |
| results/benchmark_v8_3_ablation/core4_plus_gse213628/02_v74_lgbm_robust_p095_top1000_shift015/heldout_gse80672/predictions.csv | 255 | 178 | 0.8744 | 0.3529 | 1.5393 | 0.0 | computed | computed_from_existing_file | 0.527 | 0.2118 | 0.0783 | 0.310875 | computed |
| results/benchmark_v7_5_gse121141_harmonization/05_lgbm_quantile_p08_top1000_agebin08/predictions.csv | 255 | 169 | 0.8072 | 0.3373 | 1.2558 | 0.0 | computed | not_available_not_recomputed |  |  |  |  |  |
| results/benchmark_v8_3_ablation/core4/01_v75_lgbm_quantile_p08_top1000_agebin08/groupkfold/predictions.csv | 255 | 169 | 0.8072 | 0.3373 | 1.2558 | 0.0 | computed | not_available_not_recomputed |  |  |  |  |  |
| results/benchmark_v8_3_ablation/core4/01_v75_lgbm_quantile_p08_top1000_agebin08/heldout_gse80672/predictions.csv | 255 | 178 | 0.8072 | 0.3373 | 1.2558 | 0.0 | computed | computed_from_existing_file | 0.5284 | 0.2289 | 0.1179 | 0.30188 | computed |
| results/benchmark_v8_3_ablation/core4/01_v75_lgbm_quantile_p08_top1000_agebin08/lodo/GSE80672/predictions.csv | 255 | 178 | 0.8072 | 0.3373 | 1.2558 | 0.0 | computed | not_available_not_recomputed |  |  |  |  |  |
| results/validation_v7_5_gse80672_cr_all_except/05_lgbm_quantile_p08_top1000_agebin08/predictions.csv | 255 | 169 | 0.8072 | 0.3373 | 1.2558 | 0.0 | computed | computed_from_existing_file | 0.5284 | 0.2289 | 0.1179 | 0.30188 | computed |
| results/benchmark_v8_3_ablation/core4_plus_gse213628/01_v75_lgbm_quantile_p08_top1000_agebin08/lodo/GSE80672/predictions.csv | 255 | 178 | 0.7997 | 0.325 | 1.1727 | 0.0 | computed | not_available_not_recomputed |  |  |  |  |  |
| results/benchmark_v8_3_ablation/core4_plus_gse213628/01_v75_lgbm_quantile_p08_top1000_agebin08/groupkfold/predictions.csv | 255 | 178 | 0.7997 | 0.325 | 1.1727 | 0.0 | computed | not_available_not_recomputed |  |  |  |  |  |
| results/benchmark_v8_3_ablation/core4_plus_gse213628/01_v75_lgbm_quantile_p08_top1000_agebin08/heldout_gse80672/predictions.csv | 255 | 178 | 0.7997 | 0.325 | 1.1727 | 0.0 | computed | computed_from_existing_file | 0.5392 | 0.225 | 0.1372 | 0.236888 | computed |
| results/benchmark_v7_5_gse121141_harmonization/03_lgbm_robust_p095_top1000_agebin095/predictions.csv | 255 | 169 | 0.7904 | 0.3275 | 1.132 | 0.0 | computed | not_available_not_recomputed |  |  |  |  |  |
| results/benchmark_v8_2_prefilter_liftover/01_v75_lgbm_quantile_p08_top1000_agebin08/predictions.csv | 255 | 178 | 0.7797 | 0.3355 | 1.0603 | 0.0 | computed | not_available_not_recomputed |  |  |  |  |  |
| results/benchmark_v8_3_ablation/all6/01_v75_lgbm_quantile_p08_top1000_agebin08/groupkfold/predictions.csv | 255 | 178 | 0.7797 | 0.3355 | 1.0603 | 0.0 | computed | not_available_not_recomputed |  |  |  |  |  |
| results/benchmark_v8_3_ablation/all6/01_v75_lgbm_quantile_p08_top1000_agebin08/lodo/GSE80672/predictions.csv | 255 | 178 | 0.7797 | 0.3355 | 1.0603 | 0.0 | computed | not_available_not_recomputed |  |  |  |  |  |
| results/benchmark_v8_3_ablation/all6/01_v75_lgbm_quantile_p08_top1000_agebin08/heldout_gse80672/predictions.csv | 255 | 178 | 0.7797 | 0.3355 | 1.0603 | 0.0 | computed | computed_from_existing_file | 0.5137 | 0.2452 | -0.0401 | 0.401336 | computed |
| results/benchmark_v8_3_ablation/core4/02_v74_lgbm_robust_p095_top1000_shift015/groupkfold/predictions.csv | 255 | 169 | 0.7649 | 0.3373 | 1.0137 | 1e-06 | computed | not_available_not_recomputed |  |  |  |  |  |
| results/benchmark_v7_5_gse121141_harmonization/02_lgbm_robust_p095_top1000_agebin08/predictions.csv | 255 | 169 | 0.7649 | 0.3373 | 1.0137 | 1e-06 | computed | not_available_not_recomputed |  |  |  |  |  |
| results/benchmark_v7_5_gse121141_harmonization/01_v74_best_no_agebin/predictions.csv | 255 | 169 | 0.7649 | 0.3373 | 1.0137 | 1e-06 | computed | not_available_not_recomputed |  |  |  |  |  |
| results/benchmark_v8_3_ablation/core4/02_v74_lgbm_robust_p095_top1000_shift015/heldout_gse80672/predictions.csv | 255 | 178 | 0.7649 | 0.3373 | 1.0137 | 1e-06 | computed | computed_from_existing_file | 0.4802 | 0.1807 | -0.0554 | 0.641571 | computed |
| results/benchmark_v8_3_ablation/all6/02_v74_lgbm_robust_p095_top1000_shift015/groupkfold/predictions.csv | 255 | 178 | 0.7424 | 0.3038 | 0.8891 | 5e-06 | computed | not_available_not_recomputed |  |  |  |  |  |

## Interpretation

Support-covered rows perform better than unsupported stress-test rows. This supports adopting the v12 benchmark redefinition.

## Outputs

- `results/benchmark_v12_1_redefined/prediction_support_annotations.csv`
- `results/benchmark_v12_1_redefined/support_covered_headline_metrics.json`
- `results/benchmark_v12_1_redefined/unsupported_stress_test_metrics.json`
- `results/benchmark_v12_1_redefined/gse121141_old104_stress_metrics.json`
- `results/benchmark_v12_1_redefined/gse80672_cr_redefined_metrics.json`
- `results/benchmark_v12_1_redefined/redefined_benchmark_summary.csv`

# v8.6 Leakage-Safe Calibration POC Report

Date: 2026-05-18

## Summary

v8.6 tested calibration methods on existing v8.3 LODO predictions only. It did not retrain models, download data, or start autoresearch.

- Calibration fit data: non-GSE121141 LODO residuals only.
- Target labels used for calibration: False.
- Methods: global residual offset, tissue offset, predicted-age-bin offset, tissue+predicted-age-bin offset.
- Scopes: all non-target tissues and shared GSE121141 target tissues only.
- Best GSE121141 old104+ delta MAE: 13.48 weeks.
- Non-target validation mean delta for that best method: -5.733 weeks.

## GSE121141 Old104+ Results

| variant | config | method | scope | n | mae_weeks | raw_mae_weeks | delta_mae_vs_raw | bias_weeks |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| all6 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | predbin_offset | all_non_target | 19 | 68.736 | 78.09 | 9.354 | 68.736 |
| all6 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | tissue_predbin_offset | all_non_target | 19 | 68.754 | 78.09 | 9.336 | 68.754 |
| all6 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | tissue_predbin_offset | shared_target_tissues_only | 19 | 68.754 | 78.09 | 9.336 | 68.754 |
| all6 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | predbin_offset | shared_target_tissues_only | 19 | 69.225 | 78.09 | 8.865 | 69.225 |
| all6 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | global_offset | all_non_target | 19 | 70.788 | 78.09 | 7.302 | 70.788 |
| all6 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | tissue_offset | all_non_target | 19 | 72.52 | 78.09 | 5.57 | 72.52 |
| all6 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | tissue_offset | shared_target_tissues_only | 19 | 72.52 | 78.09 | 5.57 | 72.52 |
| all6 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | global_offset | shared_target_tissues_only | 19 | 77.58 | 78.09 | 0.51 | 77.58 |
| all6 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | raw | none | 19 | 78.09 | 78.09 | 0.0 | 78.09 |
| core4 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | predbin_offset | all_non_target | 19 | 70.455 | 75.386 | 4.931 | 68.659 |
| core4 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | global_offset | all_non_target | 19 | 72.687 | 75.386 | 2.699 | 72.234 |
| core4 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | raw | none | 19 | 75.386 | 75.386 | 0.0 | 75.251 |
| core4 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | predbin_offset | shared_target_tissues_only | 19 | 75.512 | 75.386 | -0.126 | 75.512 |
| core4 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | global_offset | shared_target_tissues_only | 19 | 76.914 | 75.386 | -1.528 | 76.914 |
| core4 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | tissue_predbin_offset | all_non_target | 19 | 78.947 | 75.386 | -3.561 | 78.69 |
| core4 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | tissue_predbin_offset | shared_target_tissues_only | 19 | 78.947 | 75.386 | -3.561 | 78.69 |
| core4 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | tissue_offset | all_non_target | 19 | 80.753 | 75.386 | -5.367 | 80.496 |
| core4 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | tissue_offset | shared_target_tissues_only | 19 | 80.753 | 75.386 | -5.367 | 80.496 |
| core4_plus_gse213628 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | global_offset | all_non_target | 19 | 63.49 | 75.648 | 12.158 | 63.387 |
| core4_plus_gse213628 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | predbin_offset | shared_target_tissues_only | 19 | 64.394 | 75.648 | 11.254 | 64.394 |
| core4_plus_gse213628 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | tissue_predbin_offset | all_non_target | 19 | 67.953 | 75.648 | 7.695 | 67.953 |
| core4_plus_gse213628 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | tissue_predbin_offset | shared_target_tissues_only | 19 | 67.953 | 75.648 | 7.695 | 67.953 |
| core4_plus_gse213628 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | predbin_offset | all_non_target | 19 | 68.285 | 75.648 | 7.363 | 68.285 |
| core4_plus_gse213628 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | global_offset | shared_target_tissues_only | 19 | 72.082 | 75.648 | 3.566 | 72.082 |
| core4_plus_gse213628 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | tissue_offset | all_non_target | 19 | 72.663 | 75.648 | 2.985 | 72.663 |
| core4_plus_gse213628 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | tissue_offset | shared_target_tissues_only | 19 | 72.663 | 75.648 | 2.985 | 72.663 |
| core4_plus_gse213628 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | raw | none | 19 | 75.648 | 75.648 | 0.0 | 75.648 |
| core4_plus_gse60012 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | predbin_offset | all_non_target | 19 | 68.73 | 82.21 | 13.48 | 68.543 |
| core4_plus_gse60012 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | global_offset | all_non_target | 19 | 72.867 | 82.21 | 9.343 | 72.867 |
| core4_plus_gse60012 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | tissue_predbin_offset | all_non_target | 19 | 77.024 | 82.21 | 5.186 | 77.024 |
| core4_plus_gse60012 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | tissue_predbin_offset | shared_target_tissues_only | 19 | 77.024 | 82.21 | 5.186 | 77.024 |
| core4_plus_gse60012 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | predbin_offset | shared_target_tissues_only | 19 | 78.37 | 82.21 | 3.84 | 78.37 |
| core4_plus_gse60012 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | tissue_offset | all_non_target | 19 | 78.568 | 82.21 | 3.642 | 78.568 |
| core4_plus_gse60012 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | tissue_offset | shared_target_tissues_only | 19 | 78.568 | 82.21 | 3.642 | 78.568 |
| core4_plus_gse60012 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | global_offset | shared_target_tissues_only | 19 | 80.006 | 82.21 | 2.204 | 80.006 |
| core4_plus_gse60012 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | raw | none | 19 | 82.21 | 82.21 | 0.0 | 82.21 |

## GSE121141 All-Age Results

| variant | method | scope | n | mae_weeks | raw_mae_weeks | delta_mae_vs_raw | bias_weeks |
| --- | --- | --- | --- | --- | --- | --- | --- |
| all6 | tissue_predbin_offset | all_non_target | 81 | 32.896 | 37.499 | 4.603 | 23.881 |
| all6 | tissue_predbin_offset | shared_target_tissues_only | 81 | 32.896 | 37.499 | 4.603 | 23.881 |
| all6 | global_offset | all_non_target | 81 | 33.602 | 37.499 | 3.897 | 26.319 |
| all6 | predbin_offset | shared_target_tissues_only | 81 | 34.401 | 37.499 | 3.098 | 25.796 |
| all6 | predbin_offset | all_non_target | 81 | 34.57 | 37.499 | 2.929 | 24.992 |
| all6 | tissue_offset | all_non_target | 81 | 34.852 | 37.499 | 2.647 | 28.405 |
| all6 | tissue_offset | shared_target_tissues_only | 81 | 34.852 | 37.499 | 2.647 | 28.405 |
| all6 | global_offset | shared_target_tissues_only | 81 | 37.179 | 37.499 | 0.32 | 33.11 |
| all6 | raw | none | 81 | 37.499 | 37.499 | 0.0 | 33.62 |
| core4 | global_offset | all_non_target | 81 | 36.191 | 38.033 | 1.842 | 28.899 |
| core4 | predbin_offset | all_non_target | 81 | 37.402 | 38.033 | 0.631 | 29.245 |
| core4 | raw | none | 81 | 38.033 | 38.033 | 0.0 | 31.916 |
| core4 | global_offset | shared_target_tissues_only | 81 | 39.131 | 38.033 | -1.098 | 33.579 |
| core4 | predbin_offset | shared_target_tissues_only | 81 | 39.841 | 38.033 | -1.808 | 33.399 |
| core4 | tissue_predbin_offset | all_non_target | 81 | 42.152 | 38.033 | -4.119 | 36.22 |
| core4 | tissue_predbin_offset | shared_target_tissues_only | 81 | 42.152 | 38.033 | -4.119 | 36.22 |
| core4 | tissue_offset | all_non_target | 81 | 42.413 | 38.033 | -4.38 | 36.953 |
| core4 | tissue_offset | shared_target_tissues_only | 81 | 42.413 | 38.033 | -4.38 | 36.953 |
| core4_plus_gse213628 | global_offset | all_non_target | 81 | 30.791 | 35.433 | 4.642 | 16.884 |
| core4_plus_gse213628 | predbin_offset | shared_target_tissues_only | 81 | 31.602 | 35.433 | 3.831 | 19.059 |
| core4_plus_gse213628 | predbin_offset | all_non_target | 81 | 33.275 | 35.433 | 2.158 | 21.106 |
| core4_plus_gse213628 | global_offset | shared_target_tissues_only | 81 | 33.685 | 35.433 | 1.748 | 25.579 |
| core4_plus_gse213628 | tissue_predbin_offset | all_non_target | 81 | 34.808 | 35.433 | 0.625 | 19.641 |
| core4_plus_gse213628 | tissue_predbin_offset | shared_target_tissues_only | 81 | 34.808 | 35.433 | 0.625 | 19.641 |

## Non-Target LODO Calibration Validation

For each non-target dataset, calibration was fitted from all other non-target datasets, excluding both GSE121141 and the evaluated dataset.

| variant | method | scope | eval_datasets | mean_delta_mae_vs_raw | median_delta_mae_vs_raw | mean_mae_weeks | mean_raw_mae_weeks |
| --- | --- | --- | --- | --- | --- | --- | --- |
| all6 | global_offset | shared_target_tissues_only | 5 | -2.398 | -0.05 | 21.354 | 18.956 |
| all6 | predbin_offset | shared_target_tissues_only | 5 | -2.592 | -2.078 | 21.549 | 18.956 |
| all6 | global_offset | all_non_target | 5 | -3.444 | -1.216 | 22.401 | 18.956 |
| all6 | tissue_offset | shared_target_tissues_only | 5 | -3.465 | -1.291 | 22.421 | 18.956 |
| all6 | tissue_offset | all_non_target | 5 | -3.552 | -2.212 | 22.508 | 18.956 |
| all6 | tissue_predbin_offset | shared_target_tissues_only | 5 | -3.862 | -3.649 | 22.818 | 18.956 |
| all6 | predbin_offset | all_non_target | 5 | -5.441 | -2.058 | 24.397 | 18.956 |
| all6 | tissue_predbin_offset | all_non_target | 5 | -6.684 | -4.336 | 25.64 | 18.956 |
| core4 | global_offset | all_non_target | 3 | -1.564 | -1.011 | 21.323 | 19.759 |
| core4 | predbin_offset | all_non_target | 3 | -1.757 | -1.773 | 21.516 | 19.759 |
| core4 | global_offset | shared_target_tissues_only | 3 | -2.002 | -2.87 | 21.761 | 19.759 |
| core4 | predbin_offset | shared_target_tissues_only | 3 | -2.003 | -1.942 | 21.762 | 19.759 |
| core4 | tissue_predbin_offset | all_non_target | 3 | -2.032 | -2.421 | 21.791 | 19.759 |
| core4 | tissue_predbin_offset | shared_target_tissues_only | 3 | -2.27 | -2.501 | 22.029 | 19.759 |
| core4 | tissue_offset | shared_target_tissues_only | 3 | -2.271 | -3.206 | 22.03 | 19.759 |
| core4 | tissue_offset | all_non_target | 3 | -2.434 | -1.901 | 22.193 | 19.759 |
| core4_plus_gse213628 | tissue_offset | all_non_target | 4 | -3.766 | -3.064 | 24.708 | 20.941 |
| core4_plus_gse213628 | global_offset | shared_target_tissues_only | 4 | -4.272 | -2.118 | 25.212 | 20.941 |
| core4_plus_gse213628 | global_offset | all_non_target | 4 | -4.297 | -2.283 | 25.238 | 20.941 |
| core4_plus_gse213628 | tissue_offset | shared_target_tissues_only | 4 | -4.44 | -2.87 | 25.382 | 20.941 |
| core4_plus_gse213628 | tissue_predbin_offset | all_non_target | 4 | -4.895 | -4.512 | 25.836 | 20.941 |
| core4_plus_gse213628 | predbin_offset | shared_target_tissues_only | 4 | -5.201 | -4.825 | 26.142 | 20.941 |
| core4_plus_gse213628 | predbin_offset | all_non_target | 4 | -5.907 | -5.802 | 26.848 | 20.941 |
| core4_plus_gse213628 | tissue_predbin_offset | shared_target_tissues_only | 4 | -6.009 | -6.176 | 26.95 | 20.941 |
| core4_plus_gse60012 | global_offset | shared_target_tissues_only | 4 | -1.064 | -0.852 | 20.258 | 19.194 |
| core4_plus_gse60012 | tissue_offset | shared_target_tissues_only | 4 | -3.732 | -3.968 | 22.926 | 19.194 |
| core4_plus_gse60012 | global_offset | all_non_target | 4 | -5.688 | -4.677 | 24.882 | 19.194 |
| core4_plus_gse60012 | predbin_offset | all_non_target | 4 | -5.733 | -5.26 | 24.927 | 19.194 |
| core4_plus_gse60012 | tissue_offset | all_non_target | 4 | -5.889 | -5.736 | 25.083 | 19.194 |
| core4_plus_gse60012 | tissue_predbin_offset | all_non_target | 4 | -6.881 | -6.699 | 26.075 | 19.194 |
| core4_plus_gse60012 | predbin_offset | shared_target_tissues_only | 4 | -7.334 | -4.333 | 26.528 | 19.194 |
| core4_plus_gse60012 | tissue_predbin_offset | shared_target_tissues_only | 4 | -7.42 | -3.912 | 26.614 | 19.194 |

## Decision

Calibration does not yet provide a robust transferable fix. Do not promote it to headline benchmark or run autoresearch. Keep the conclusion that old same-tissue data support remains the main blocker.

## Outputs

- `results/v8_6_calibration_poc/gse121141_calibration_summary.csv`
- `results/v8_6_calibration_poc/gse121141_calibrated_predictions.csv`
- `results/v8_6_calibration_poc/calibration_offsets.csv`
- `results/v8_6_calibration_poc/nontarget_lodo_calibration_validation.csv`
- `results/v8_6_calibration_poc/best_old104plus_calibration_by_variant.csv`

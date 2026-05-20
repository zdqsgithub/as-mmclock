# v19 Autoresearch ML Baseline Report

Date: 2026-05-20T04:39:07.559961+00:00

## Scope

Classical ML autoresearch before deep learning. Search includes LGBM and
RandomForest architecture parameters, leakage-safe GroupKFold, LODO checks
for top candidates, random-label sanity, and final all-data baseline fit.

## Selected Baseline

- config: `cfg_014_lgbm_top1000_quantile_uniform_median_s0p125_p0p8_ne200_lr0p03_leaves31`
- model: `lgbm`
- top regions: `1000`
- preprocess/imputation: `quantile_uniform` / `median`
- filters: presence `0.8`, shift `0.125`
- GroupKFold r/MAE/R2: `0.6315` / `24.227` / `0.2605`
- LODO mean/worst MAE: `21.93` / `36.086`
- LODO CR AUC: `0.7793`
- old104 weighted MAE: `43.762`
- final score: `1.0667309999999999`

## Final All-Data Fit

- apparent r/MAE/R2: `0.9944` / `2.319` / `0.9868`
- selected features: `1000`
- apparent all-data metrics are optimistic; LODO is the evidence baseline.

## Random-Label Sanity

`{'pearson_r': 0.1199, 'pearson_pval': 2.7093383976549727e-05, 'mae_weeks': 33.026, 'medae_weeks': 22.705, 'rmse_weeks': 43.212, 'r2': -0.3763, 'mae_pct_lifespan': 0.3303, 'n_samples': 1219, 'n_datasets': 6, 'cr_detection_auc': 0.76, 'cr_detection_f1': 0.0828, 'cr_cohens_d': 0.9898, 'cr_mannwhitney_p': 0.0, 'cross_dataset_mae': 31.835, 'threshold_checks': {'pearson_r': 'FAIL 0.1199 target > 0.9', 'mae_weeks': 'FAIL 33.026 target < 3.5', 'medae_weeks': 'FAIL 22.705 target < 3.0', 'r2': 'FAIL -0.3763 target > 0.8', 'cr_detection_auc': 'FAIL 0.76 target > 0.8', 'cr_detection_f1': 'FAIL 0.0828 target > 0.75', 'cr_cohens_d': 'PASS 0.9898 > 0.8', 'cross_dataset_mae': 'FAIL 31.835 target < 5.0'}}`

## Top Ranked Configs

| config_id | model_type | n_feature_prefilter | imputation | preprocess | min_train_feature_presence | min_train_group_feature_presence | max_train_dataset_mean_shift | min_train_agebin_feature_presence | target_transform | lgbm_n_estimators | lgbm_learning_rate | lgbm_num_leaves | lgbm_subsample | lgbm_colsample_bytree | rf_n_estimators | rf_max_depth | rf_min_samples_leaf | rf_max_features | n_samples | pearson_r | mae_weeks | medae_weeks | rmse_weeks | r2 | cross_dataset_mae | cr_detection_auc | group_score | lodo_mean_mae_weeks | lodo_worst_mae_weeks | lodo_mean_pearson_r | lodo_min_pearson_r | lodo_cr_auc | old_104w_weighted_mae_weeks | lodo_score | final_score |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| cfg_014_lgbm_top1000_quantile_uniform_median_s0p125_p0p8_ne200_lr0p03_leaves31 | lgbm | 1000 | median | quantile_uniform | 0.8 | 0.8 | 0.125 | 0.8 | log1p_days | 200 | 0.03 | 31 | 0.9 | 0.9 | 150 | 10 | 2 | 1.0 | 1219 | 0.6315 | 24.227 | 20.134 | 31.675 | 0.2605 | 23.185 | 0.7524 | 0.557165 | 21.93 | 36.086 | 0.6663 | 0.4475 | 0.7793 | 43.762 | 0.509566 | 1.0667309999999999 |
| cfg_034_lgbm_top1500_none_median_s0p1_p0p8_ne200_lr0p03_leaves31 | lgbm | 1500 | median | none | 0.8 | 0.8 | 0.1 | 0.8 | log1p_days | 200 | 0.03 | 31 | 0.9 | 0.9 | 150 | 10 | 2 | 1.0 | 1219 | 0.649 | 23.994 | 18.778 | 31.887 | 0.2506 | 23.591 | 0.7248 | 0.557248 | 21.519 | 36.188 | 0.6696 | 0.4421 | 0.7592 | 49.806 | 0.505538 | 1.062786 |
| cfg_036_lgbm_top1000_quantile_uniform_median_s0p1_p0p7_ne200_lr0p03_leaves31 | lgbm | 1000 | median | quantile_uniform | 0.7 | 0.7 | 0.1 | 0.7 | log1p_days | 200 | 0.03 | 31 | 0.9 | 0.9 | 150 | 10 | 2 | 1.0 | 1219 | 0.63 | 23.943 | 17.998 | 31.426 | 0.2721 | 23.924 | 0.7155 | 0.553684 | 22.08 | 34.506 | 0.6583 | 0.4881 | 0.7472 | 49.541 | 0.502573 | 1.056257 |
| cfg_049_lgbm_top1500_quantile_uniform_median_s0p1_p0p8_ne350_lr0p02_leaves31 | lgbm | 1500 | median | quantile_uniform | 0.8 | 0.8 | 0.1 | 0.8 | log1p_days | 350 | 0.02 | 31 | 0.9 | 0.9 | 150 | 10 | 2 | 1.0 | 1219 | 0.6558 | 23.868 | 17.939 | 31.851 | 0.2523 | 23.56 | 0.6789 | 0.552398 | 21.568 | 36.53 | 0.6718 | 0.4657 | 0.7021 | 51.635 | 0.503054 | 1.0554519999999998 |
| cfg_022_lgbm_top2000_quantile_uniform_median_s0p075_p0p8_ne200_lr0p03_leaves31 | lgbm | 2000 | median | quantile_uniform | 0.8 | 0.8 | 0.075 | 0.8 | log1p_days | 200 | 0.03 | 31 | 0.9 | 0.9 | 150 | 10 | 2 | 1.0 | 1219 | 0.633 | 24.182 | 18.642 | 31.813 | 0.2541 | 23.394 | 0.7697 | 0.559463 | 22.038 | 35.247 | 0.6314 | 0.4698 | 0.808 | 48.124 | 0.493719 | 1.053182 |
| cfg_033_lgbm_top1000_none_median_s0p1_p0p8_ne200_lr0p03_leaves31 | lgbm | 1000 | median | none | 0.8 | 0.8 | 0.1 | 0.8 | log1p_days | 200 | 0.03 | 31 | 0.9 | 0.9 | 150 | 10 | 2 | 1.0 | 1219 | 0.6255 | 24.024 | 18.33 | 31.552 | 0.2663 | 23.99 | 0.7282 | 0.552598 | 22.216 | 35.968 | 0.6634 | 0.4895 | 0.7609 | 49.816 | 0.499957 | 1.052555 |
| cfg_030_lgbm_top1000_robust_median_s0p1_p0p8_ne200_lr0p03_leaves31 | lgbm | 1000 | median | robust | 0.8 | 0.8 | 0.1 | 0.8 | log1p_days | 200 | 0.03 | 31 | 0.9 | 0.9 | 150 | 10 | 2 | 1.0 | 1219 | 0.6238 | 24.006 | 18.376 | 31.55 | 0.2663 | 24.1 | 0.7306 | 0.552444 | 22.262 | 36.092 | 0.6602 | 0.4729 | 0.7777 | 49.546 | 0.498657 | 1.051101 |
| cfg_038_lgbm_top1000_quantile_uniform_median_s0p075_p0p9_ne200_lr0p03_leaves31 | lgbm | 1000 | median | quantile_uniform | 0.9 | 0.9 | 0.075 | 0.9 | log1p_days | 200 | 0.03 | 31 | 0.9 | 0.9 | 150 | 10 | 2 | 1.0 | 1219 | 0.6248 | 24.578 | 18.53 | 32.413 | 0.2257 | 24.226 | 0.7914 | 0.552117 | 22.26 | 34.64 | 0.6425 | 0.5069 | 0.8104 | 49.804 | 0.495567 | 1.0476839999999998 |
| cfg_013_lgbm_top1000_quantile_uniform_median_s0p1_p0p8_ne200_lr0p03_leaves31 | lgbm | 1000 | median | quantile_uniform | 0.8 | 0.8 | 0.1 | 0.8 | log1p_days | 200 | 0.03 | 31 | 0.9 | 0.9 | 150 | 10 | 2 | 1.0 | 1219 | 0.627 | 23.956 | 17.952 | 31.505 | 0.2684 | 23.993 | 0.7162 | 0.551739 |  |  |  |  |  |  | -1.0 | 0.551739 |
| cfg_045_lgbm_top1000_quantile_uniform_median_s0p1_p0p8_ne250_lr0p025_leaves31 | lgbm | 1000 | median | quantile_uniform | 0.8 | 0.8 | 0.1 | 0.8 | log1p_days | 250 | 0.025 | 31 | 0.9 | 0.9 | 150 | 10 | 2 | 1.0 | 1219 | 0.6266 | 23.98 | 18.361 | 31.471 | 0.27 | 24.144 | 0.7119 | 0.551117 |  |  |  |  |  |  | -1.0 | 0.551117 |
| cfg_005_lgbm_top500_quantile_uniform_median_s0p15_p0p8_ne200_lr0p03_leaves31 | lgbm | 500 | median | quantile_uniform | 0.8 | 0.8 | 0.15 | 0.8 | log1p_days | 200 | 0.03 | 31 | 0.9 | 0.9 | 150 | 10 | 2 | 1.0 | 1219 | 0.6166 | 24.484 | 18.544 | 32.191 | 0.2362 | 23.967 | 0.787 | 0.551043 |  |  |  |  |  |  | -1.0 | 0.551043 |
| cfg_019_lgbm_top1500_quantile_uniform_median_s0p125_p0p8_ne200_lr0p03_leaves31 | lgbm | 1500 | median | quantile_uniform | 0.8 | 0.8 | 0.125 | 0.8 | log1p_days | 200 | 0.03 | 31 | 0.9 | 0.9 | 150 | 10 | 2 | 1.0 | 1219 | 0.6276 | 25.249 | 20.178 | 33.109 | 0.192 | 23.865 | 0.8343 | 0.550579 |  |  |  |  |  |  | -1.0 | 0.550579 |
| cfg_021_lgbm_top2000_quantile_uniform_median_s0p05_p0p8_ne200_lr0p03_leaves31 | lgbm | 2000 | median | quantile_uniform | 0.8 | 0.8 | 0.05 | 0.8 | log1p_days | 200 | 0.03 | 31 | 0.9 | 0.9 | 150 | 10 | 2 | 1.0 | 1219 | 0.614 | 24.957 | 19.961 | 32.421 | 0.2253 | 22.617 | 0.8132 | 0.550578 |  |  |  |  |  |  | -1.0 | 0.550578 |
| cfg_017_lgbm_top1500_quantile_uniform_median_s0p075_p0p8_ne200_lr0p03_leaves31 | lgbm | 1500 | median | quantile_uniform | 0.8 | 0.8 | 0.075 | 0.8 | log1p_days | 200 | 0.03 | 31 | 0.9 | 0.9 | 150 | 10 | 2 | 1.0 | 1219 | 0.6302 | 24.244 | 18.305 | 32.001 | 0.2452 | 23.791 | 0.7265 | 0.548281 |  |  |  |  |  |  | -1.0 | 0.548281 |
| cfg_023_lgbm_top2000_quantile_uniform_median_s0p1_p0p8_ne200_lr0p03_leaves31 | lgbm | 2000 | median | quantile_uniform | 0.8 | 0.8 | 0.1 | 0.8 | log1p_days | 200 | 0.03 | 31 | 0.9 | 0.9 | 150 | 10 | 2 | 1.0 | 1219 | 0.6451 | 24.055 | 18.127 | 31.549 | 0.2664 | 23.707 | 0.6631 | 0.548249 |  |  |  |  |  |  | -1.0 | 0.548249 |

## LODO Rows

| config_id | dataset | n_samples | pearson_r | mae_weeks | rmse_weeks | r2 | cr_detection_auc | old_104w_n | old_104w_mae_weeks |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| cfg_022_lgbm_top2000_quantile_uniform_median_s0p075_p0p8_ne200_lr0p03_leaves31 | GSE120137 | 549 | 0.5497 | 25.036 | 31.958 | -0.2241 |  | 0 |  |
| cfg_022_lgbm_top2000_quantile_uniform_median_s0p075_p0p8_ne200_lr0p03_leaves31 | GSE80672 | 255 | 0.7848 | 28.141 | 33.953 | 0.462 | 0.808 | 74 | 47.725 |
| cfg_022_lgbm_top2000_quantile_uniform_median_s0p075_p0p8_ne200_lr0p03_leaves31 | GSE93957 | 62 | 0.7995 | 7.701 | 9.234 | 0.6 |  | 0 |  |
| cfg_022_lgbm_top2000_quantile_uniform_median_s0p075_p0p8_ne200_lr0p03_leaves31 | GSE121141 | 81 | 0.5069 | 35.247 | 47.361 | -0.5032 |  | 19 | 77.386 |
| cfg_022_lgbm_top2000_quantile_uniform_median_s0p075_p0p8_ne200_lr0p03_leaves31 | GSE60012 | 152 | 0.678 | 9.274 | 12.589 | -1.0993 |  | 0 |  |
| cfg_022_lgbm_top2000_quantile_uniform_median_s0p075_p0p8_ne200_lr0p03_leaves31 | GSE213628 | 120 | 0.4698 | 26.828 | 33.17 | 0.1484 |  | 42 | 35.588 |
| cfg_034_lgbm_top1500_none_median_s0p1_p0p8_ne200_lr0p03_leaves31 | GSE120137 | 549 | 0.6189 | 24.202 | 31.09 | -0.1585 |  | 0 |  |
| cfg_034_lgbm_top1500_none_median_s0p1_p0p8_ne200_lr0p03_leaves31 | GSE80672 | 255 | 0.7867 | 29.54 | 36.81 | 0.3676 | 0.7592 | 74 | 55.037 |
| cfg_034_lgbm_top1500_none_median_s0p1_p0p8_ne200_lr0p03_leaves31 | GSE93957 | 62 | 0.8418 | 6.766 | 8.182 | 0.6859 |  | 0 |  |
| cfg_034_lgbm_top1500_none_median_s0p1_p0p8_ne200_lr0p03_leaves31 | GSE121141 | 81 | 0.4421 | 36.188 | 46.942 | -0.4768 |  | 19 | 75.692 |
| cfg_034_lgbm_top1500_none_median_s0p1_p0p8_ne200_lr0p03_leaves31 | GSE60012 | 152 | 0.7177 | 9.084 | 11.501 | -0.752 |  | 0 |  |
| cfg_034_lgbm_top1500_none_median_s0p1_p0p8_ne200_lr0p03_leaves31 | GSE213628 | 120 | 0.6105 | 23.335 | 28.604 | 0.3667 |  | 42 | 28.879 |
| cfg_014_lgbm_top1000_quantile_uniform_median_s0p125_p0p8_ne200_lr0p03_leaves31 | GSE120137 | 549 | 0.5758 | 25.998 | 32.927 | -0.2995 |  | 0 |  |
| cfg_014_lgbm_top1000_quantile_uniform_median_s0p125_p0p8_ne200_lr0p03_leaves31 | GSE80672 | 255 | 0.8032 | 26.12 | 32.384 | 0.5105 | 0.7793 | 74 | 44.778 |
| cfg_014_lgbm_top1000_quantile_uniform_median_s0p125_p0p8_ne200_lr0p03_leaves31 | GSE93957 | 62 | 0.8074 | 7.797 | 9.124 | 0.6094 |  | 0 |  |
| cfg_014_lgbm_top1000_quantile_uniform_median_s0p125_p0p8_ne200_lr0p03_leaves31 | GSE121141 | 81 | 0.4475 | 36.086 | 47.125 | -0.4883 |  | 19 | 75.545 |
| cfg_014_lgbm_top1000_quantile_uniform_median_s0p125_p0p8_ne200_lr0p03_leaves31 | GSE60012 | 152 | 0.7593 | 11.639 | 14.645 | -1.8412 |  | 0 |  |
| cfg_014_lgbm_top1000_quantile_uniform_median_s0p125_p0p8_ne200_lr0p03_leaves31 | GSE213628 | 120 | 0.6048 | 23.94 | 28.965 | 0.3507 |  | 42 | 27.595 |
| cfg_036_lgbm_top1000_quantile_uniform_median_s0p1_p0p7_ne200_lr0p03_leaves31 | GSE120137 | 549 | 0.6011 | 23.001 | 29.295 | -0.0286 |  | 0 |  |
| cfg_036_lgbm_top1000_quantile_uniform_median_s0p1_p0p7_ne200_lr0p03_leaves31 | GSE80672 | 255 | 0.7712 | 29.858 | 37.359 | 0.3486 | 0.7472 | 74 | 55.63 |
| cfg_036_lgbm_top1000_quantile_uniform_median_s0p1_p0p7_ne200_lr0p03_leaves31 | GSE93957 | 62 | 0.8254 | 7.455 | 8.975 | 0.6221 |  | 0 |  |
| cfg_036_lgbm_top1000_quantile_uniform_median_s0p1_p0p7_ne200_lr0p03_leaves31 | GSE121141 | 81 | 0.4881 | 34.506 | 45.047 | -0.3599 |  | 19 | 72.238 |
| cfg_036_lgbm_top1000_quantile_uniform_median_s0p1_p0p7_ne200_lr0p03_leaves31 | GSE60012 | 152 | 0.7077 | 11.312 | 14.474 | -1.7751 |  | 0 |  |
| cfg_036_lgbm_top1000_quantile_uniform_median_s0p1_p0p7_ne200_lr0p03_leaves31 | GSE213628 | 120 | 0.556 | 26.347 | 31.724 | 0.221 |  | 42 | 28.546 |
| cfg_033_lgbm_top1000_none_median_s0p1_p0p8_ne200_lr0p03_leaves31 | GSE120137 | 549 | 0.5959 | 23.144 | 29.478 | -0.0415 |  | 0 |  |
| cfg_033_lgbm_top1000_none_median_s0p1_p0p8_ne200_lr0p03_leaves31 | GSE80672 | 255 | 0.7647 | 30.124 | 37.399 | 0.3472 | 0.7609 | 74 | 55.752 |
| cfg_033_lgbm_top1000_none_median_s0p1_p0p8_ne200_lr0p03_leaves31 | GSE93957 | 62 | 0.8479 | 7.168 | 8.494 | 0.6615 |  | 0 |  |
| cfg_033_lgbm_top1000_none_median_s0p1_p0p8_ne200_lr0p03_leaves31 | GSE121141 | 81 | 0.4895 | 35.968 | 46.415 | -0.4438 |  | 19 | 74.874 |
| cfg_033_lgbm_top1000_none_median_s0p1_p0p8_ne200_lr0p03_leaves31 | GSE60012 | 152 | 0.7129 | 10.781 | 14.053 | -1.616 |  | 0 |  |
| cfg_033_lgbm_top1000_none_median_s0p1_p0p8_ne200_lr0p03_leaves31 | GSE213628 | 120 | 0.5694 | 26.11 | 31.395 | 0.2371 |  | 42 | 28.02 |
| cfg_030_lgbm_top1000_robust_median_s0p1_p0p8_ne200_lr0p03_leaves31 | GSE120137 | 549 | 0.5875 | 23.18 | 29.459 | -0.0401 |  | 0 |  |
| cfg_030_lgbm_top1000_robust_median_s0p1_p0p8_ne200_lr0p03_leaves31 | GSE80672 | 255 | 0.776 | 29.496 | 36.962 | 0.3624 | 0.7777 | 74 | 55.204 |
| cfg_030_lgbm_top1000_robust_median_s0p1_p0p8_ne200_lr0p03_leaves31 | GSE93957 | 62 | 0.8495 | 7.253 | 8.589 | 0.6539 |  | 0 |  |
| cfg_030_lgbm_top1000_robust_median_s0p1_p0p8_ne200_lr0p03_leaves31 | GSE121141 | 81 | 0.4729 | 36.092 | 46.455 | -0.4462 |  | 19 | 74.346 |
| cfg_030_lgbm_top1000_robust_median_s0p1_p0p8_ne200_lr0p03_leaves31 | GSE60012 | 152 | 0.7292 | 10.54 | 13.532 | -1.4257 |  | 0 |  |
| cfg_030_lgbm_top1000_robust_median_s0p1_p0p8_ne200_lr0p03_leaves31 | GSE213628 | 120 | 0.5458 | 27.011 | 32.029 | 0.206 |  | 42 | 28.359 |
| cfg_049_lgbm_top1500_quantile_uniform_median_s0p1_p0p8_ne350_lr0p02_leaves31 | GSE120137 | 549 | 0.6224 | 23.799 | 30.688 | -0.1288 |  | 0 |  |
| cfg_049_lgbm_top1500_quantile_uniform_median_s0p1_p0p8_ne350_lr0p02_leaves31 | GSE80672 | 255 | 0.7818 | 29.958 | 37.64 | 0.3388 | 0.7021 | 74 | 56.84 |
| cfg_049_lgbm_top1500_quantile_uniform_median_s0p1_p0p8_ne350_lr0p02_leaves31 | GSE93957 | 62 | 0.822 | 6.965 | 8.423 | 0.6672 |  | 0 |  |
| cfg_049_lgbm_top1500_quantile_uniform_median_s0p1_p0p8_ne350_lr0p02_leaves31 | GSE121141 | 81 | 0.4657 | 36.53 | 47.434 | -0.5079 |  | 19 | 77.103 |
| cfg_049_lgbm_top1500_quantile_uniform_median_s0p1_p0p8_ne350_lr0p02_leaves31 | GSE60012 | 152 | 0.7294 | 8.786 | 11.023 | -0.6094 |  | 0 |  |
| cfg_049_lgbm_top1500_quantile_uniform_median_s0p1_p0p8_ne350_lr0p02_leaves31 | GSE213628 | 120 | 0.6093 | 23.368 | 28.539 | 0.3696 |  | 42 | 30.944 |
| cfg_038_lgbm_top1000_quantile_uniform_median_s0p075_p0p9_ne200_lr0p03_leaves31 | GSE120137 | 549 | 0.5983 | 24.579 | 31.607 | -0.1973 |  | 0 |  |
| cfg_038_lgbm_top1000_quantile_uniform_median_s0p075_p0p9_ne200_lr0p03_leaves31 | GSE80672 | 255 | 0.7479 | 29.912 | 36.326 | 0.3841 | 0.8104 | 74 | 51.531 |
| cfg_038_lgbm_top1000_quantile_uniform_median_s0p075_p0p9_ne200_lr0p03_leaves31 | GSE93957 | 62 | 0.7928 | 8.601 | 10.402 | 0.4924 |  | 0 |  |
| cfg_038_lgbm_top1000_quantile_uniform_median_s0p075_p0p9_ne200_lr0p03_leaves31 | GSE121141 | 81 | 0.5254 | 34.64 | 46.453 | -0.4461 |  | 19 | 76.698 |
| cfg_038_lgbm_top1000_quantile_uniform_median_s0p075_p0p9_ne200_lr0p03_leaves31 | GSE60012 | 152 | 0.6839 | 9.527 | 12.908 | -1.2071 |  | 0 |  |
| cfg_038_lgbm_top1000_quantile_uniform_median_s0p075_p0p9_ne200_lr0p03_leaves31 | GSE213628 | 120 | 0.5069 | 26.301 | 32.349 | 0.1901 |  | 42 | 34.595 |

## Outputs

- `results/autoresearch_v19_ml_baseline/v19_autoresearch_summary.csv`
- `results/autoresearch_v19_ml_baseline/v19_lodo_summary.csv`
- `results/autoresearch_v19_ml_baseline/final_all_data_model/final_all_data_ml_baseline_model.joblib`

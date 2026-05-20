# v18 Baseline Parameter Search Report

Date: 2026-05-20T03:17:54.712962+00:00

## Scope

Classical ML baseline search before deep learning. Training uses existing
fold-internal feature selection, train-only imputation/scaling, dataset
GroupKFold, and leave-one-dataset-out checks for the top candidates.

## Best Candidate

- config: `cfg_008_lgbm_top1000_quantile_uniform_median_p0p8_gp0p8_s0p1_ab0p8_log1p_days`
- model: `lgbm`
- top regions: `1000`
- preprocess/imputation: `quantile_uniform` / `median`
- filters: presence `0.8`, group `0.8`, shift `0.1`, age-bin `0.8`
- GroupKFold r/MAE/R2: `0.627` / `23.956` / `0.2684`
- LODO mean/worst MAE: `22.344` / `36.099`
- old 104w weighted MAE: `49.738`
- final score: `1.065829`

## Random-Label Sanity

`{'pearson_r': 0.1279, 'pearson_pval': 7.482512739013961e-06, 'mae_weeks': 32.832, 'medae_weeks': 22.757, 'rmse_weeks': 43.637, 'r2': -0.4035, 'mae_pct_lifespan': 0.3283, 'n_samples': 1219, 'n_datasets': 6, 'cr_detection_auc': 0.6191, 'cr_detection_f1': 0.0651, 'cr_cohens_d': 0.5048, 'cr_mannwhitney_p': 0.010726, 'cross_dataset_mae': 30.898, 'threshold_checks': {'pearson_r': 'FAIL 0.1279 target > 0.9', 'mae_weeks': 'FAIL 32.832 target < 3.5', 'medae_weeks': 'FAIL 22.757 target < 3.0', 'r2': 'FAIL -0.4035 target > 0.8', 'cr_detection_auc': 'FAIL 0.6191 target > 0.8', 'cr_detection_f1': 'FAIL 0.0651 target > 0.75', 'cr_cohens_d': 'FAIL 0.5048 target > 0.8', 'cross_dataset_mae': 'FAIL 30.898 target < 5.0'}}`

## Top Ranked Configs

| config_id | model_type | n_feature_prefilter | imputation | preprocess | min_train_feature_presence | min_train_group_feature_presence | max_train_dataset_mean_shift | min_train_agebin_feature_presence | target_transform | n_samples | pearson_r | mae_weeks | medae_weeks | rmse_weeks | r2 | cross_dataset_mae | cr_detection_auc | group_score | lodo_mean_mae_weeks | lodo_worst_mae_weeks | lodo_mean_pearson_r | lodo_min_pearson_r | old_104w_weighted_mae_weeks | lodo_score | final_score |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| cfg_008_lgbm_top1000_quantile_uniform_median_p0p8_gp0p8_s0p1_ab0p8_log1p_days | lgbm | 1000 | median | quantile_uniform | 0.8 | 0.8 | 0.1 | 0.8 | log1p_days | 1219 | 0.627 | 23.956 | 17.952 | 31.505 | 0.2684 | 23.993 | 0.7162 | 0.545477 | 22.344 | 36.099 | 0.659 | 0.4742 | 49.738 | 0.520352 | 1.065829 |
| cfg_005_lgbm_top1000_standard_median_p0p8_gp0p8_s0p15_ab0p8_log1p_days | lgbm | 1000 | median | standard | 0.8 | 0.8 | 0.15 | 0.8 | log1p_days | 1219 | 0.6252 | 24.85 | 20.299 | 32.483 | 0.2223 | 23.488 | 0.8004 | 0.542882 | 22.409 | 38.579 | 0.658 | 0.3708 | 47.228 | 0.516622 | 1.059504 |
| cfg_006_lgbm_top1000_robust_median_p0p8_gp0p8_s0p15_ab0p8_log1p_days | lgbm | 1000 | median | robust | 0.8 | 0.8 | 0.15 | 0.8 | log1p_days | 1219 | 0.6256 | 24.758 | 19.728 | 32.604 | 0.2165 | 23.644 | 0.7644 | 0.536498 | 22.132 | 38.01 | 0.6625 | 0.3869 | 49.341 | 0.518964 | 1.055462 |
| cfg_026_lgbm_top1000_quantile_uniform_mean_p0p8_gp0p8_s0p15_ab0p8_log1p_days | lgbm | 1000 | mean | quantile_uniform | 0.8 | 0.8 | 0.15 | 0.8 | log1p_days | 1219 | 0.6209 | 24.875 | 20.353 | 32.642 | 0.2147 | 23.44 | 0.7831 | 0.536583 | 22.258 | 38.386 | 0.6473 | 0.3743 | 46.264 | 0.51515 | 1.051733 |
| cfg_002_lgbm_top500_quantile_uniform_median_p0p8_gp0p8_s0p15_ab0p8_log1p_days | lgbm | 500 | median | quantile_uniform | 0.8 | 0.8 | 0.15 | 0.8 | log1p_days | 1219 | 0.6166 | 24.484 | 18.544 | 32.191 | 0.2362 | 23.967 | 0.787 | 0.542127 | 22.57 | 40.111 | 0.6389 | 0.3283 | 47.799 | 0.505022 | 1.047149 |
| cfg_001_lgbm_top1000_quantile_uniform_median_p0p8_gp0p8_s0p15_ab0p8_log1p_days | lgbm | 1000 | median | quantile_uniform | 0.8 | 0.8 | 0.15 | 0.8 | log1p_days | 1219 | 0.6253 | 24.843 | 19.912 | 32.757 | 0.2091 | 23.359 | 0.7567 | 0.53309 |  |  |  |  |  | -1.0 | 0.53309 |
| cfg_013_lgbm_top1000_quantile_uniform_median_p0p9_gp0p9_s0p15_ab0p9_log1p_days | lgbm | 1000 | median | quantile_uniform | 0.9 | 0.9 | 0.15 | 0.9 | log1p_days | 1219 | 0.6214 | 24.831 | 19.531 | 32.748 | 0.2095 | 23.587 | 0.7634 | 0.532675 |  |  |  |  |  | -1.0 | 0.532675 |
| cfg_007_lgbm_top1000_none_median_p0p8_gp0p8_s0p15_ab0p8_log1p_days | lgbm | 1000 | median | none | 0.8 | 0.8 | 0.15 | 0.8 | log1p_days | 1219 | 0.6207 | 25.009 | 19.873 | 32.738 | 0.21 | 23.496 | 0.7658 | 0.532287 |  |  |  |  |  | -1.0 | 0.532287 |
| cfg_011_lgbm_top1000_quantile_uniform_median_p0p6_gp0p6_s0p15_ab0p6_log1p_days | lgbm | 1000 | median | quantile_uniform | 0.6 | 0.6 | 0.15 | 0.6 | log1p_days | 1219 | 0.6238 | 24.861 | 19.828 | 32.686 | 0.2126 | 23.372 | 0.7378 | 0.53047 |  |  |  |  |  | -1.0 | 0.53047 |
| cfg_012_lgbm_top1000_quantile_uniform_median_p0p7_gp0p7_s0p15_ab0p7_log1p_days | lgbm | 1000 | median | quantile_uniform | 0.7 | 0.7 | 0.15 | 0.7 | log1p_days | 1219 | 0.6238 | 24.861 | 19.828 | 32.686 | 0.2126 | 23.372 | 0.7378 | 0.53047 |  |  |  |  |  | -1.0 | 0.53047 |

## LODO Rows

| config_id | dataset | n_samples | pearson_r | mae_weeks | rmse_weeks | r2 | cr_detection_auc | old_104w_n | old_104w_mae_weeks |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| cfg_008_lgbm_top1000_quantile_uniform_median_p0p8_gp0p8_s0p1_ab0p8_log1p_days | GSE120137 | 549 | 0.6011 | 23.001 | 29.295 | -0.0286 |  | 0 |  |
| cfg_008_lgbm_top1000_quantile_uniform_median_p0p8_gp0p8_s0p1_ab0p8_log1p_days | GSE80672 | 255 | 0.7712 | 29.858 | 37.359 | 0.3486 | 0.7472 | 74 | 55.63 |
| cfg_008_lgbm_top1000_quantile_uniform_median_p0p8_gp0p8_s0p1_ab0p8_log1p_days | GSE93957 | 62 | 0.8315 | 7.697 | 9.216 | 0.6015 |  | 0 |  |
| cfg_008_lgbm_top1000_quantile_uniform_median_p0p8_gp0p8_s0p1_ab0p8_log1p_days | GSE121141 | 81 | 0.4742 | 36.099 | 46.925 | -0.4756 |  | 19 | 75.87 |
| cfg_008_lgbm_top1000_quantile_uniform_median_p0p8_gp0p8_s0p1_ab0p8_log1p_days | GSE60012 | 152 | 0.709 | 11.103 | 14.439 | -1.7616 |  | 0 |  |
| cfg_008_lgbm_top1000_quantile_uniform_median_p0p8_gp0p8_s0p1_ab0p8_log1p_days | GSE213628 | 120 | 0.5668 | 26.304 | 31.589 | 0.2277 |  | 42 | 27.534 |
| cfg_005_lgbm_top1000_standard_median_p0p8_gp0p8_s0p15_ab0p8_log1p_days | GSE120137 | 549 | 0.5904 | 26.766 | 33.771 | -0.367 |  | 0 |  |
| cfg_005_lgbm_top1000_standard_median_p0p8_gp0p8_s0p15_ab0p8_log1p_days | GSE80672 | 255 | 0.7565 | 28.975 | 34.947 | 0.43 | 0.8103 | 74 | 48.977 |
| cfg_005_lgbm_top1000_standard_median_p0p8_gp0p8_s0p15_ab0p8_log1p_days | GSE93957 | 62 | 0.8594 | 7.795 | 9.43 | 0.5828 |  | 0 |  |
| cfg_005_lgbm_top1000_standard_median_p0p8_gp0p8_s0p15_ab0p8_log1p_days | GSE121141 | 81 | 0.3708 | 38.579 | 50.199 | -0.6888 |  | 19 | 79.596 |
| cfg_005_lgbm_top1000_standard_median_p0p8_gp0p8_s0p15_ab0p8_log1p_days | GSE60012 | 152 | 0.7696 | 8.433 | 11.11 | -0.6352 |  | 0 |  |
| cfg_005_lgbm_top1000_standard_median_p0p8_gp0p8_s0p15_ab0p8_log1p_days | GSE213628 | 120 | 0.6011 | 23.903 | 28.77 | 0.3593 |  | 42 | 29.503 |
| cfg_002_lgbm_top500_quantile_uniform_median_p0p8_gp0p8_s0p15_ab0p8_log1p_days | GSE120137 | 549 | 0.5992 | 25.205 | 31.925 | -0.2216 |  | 0 |  |
| cfg_002_lgbm_top500_quantile_uniform_median_p0p8_gp0p8_s0p15_ab0p8_log1p_days | GSE80672 | 255 | 0.7454 | 28.774 | 35.367 | 0.4162 | 0.7944 | 74 | 49.858 |
| cfg_002_lgbm_top500_quantile_uniform_median_p0p8_gp0p8_s0p15_ab0p8_log1p_days | GSE93957 | 62 | 0.7787 | 7.654 | 9.236 | 0.5998 |  | 0 |  |
| cfg_002_lgbm_top500_quantile_uniform_median_p0p8_gp0p8_s0p15_ab0p8_log1p_days | GSE121141 | 81 | 0.3283 | 40.111 | 51.87 | -0.803 |  | 19 | 82.702 |
| cfg_002_lgbm_top500_quantile_uniform_median_p0p8_gp0p8_s0p15_ab0p8_log1p_days | GSE60012 | 152 | 0.7904 | 9.926 | 12.506 | -1.0719 |  | 0 |  |
| cfg_002_lgbm_top500_quantile_uniform_median_p0p8_gp0p8_s0p15_ab0p8_log1p_days | GSE213628 | 120 | 0.5917 | 23.753 | 29.372 | 0.3323 |  | 42 | 28.381 |
| cfg_026_lgbm_top1000_quantile_uniform_mean_p0p8_gp0p8_s0p15_ab0p8_log1p_days | GSE120137 | 549 | 0.5813 | 27.112 | 34.2 | -0.4019 |  | 0 |  |
| cfg_026_lgbm_top1000_quantile_uniform_mean_p0p8_gp0p8_s0p15_ab0p8_log1p_days | GSE80672 | 255 | 0.7563 | 28.687 | 34.429 | 0.4468 | 0.7913 | 74 | 47.259 |
| cfg_026_lgbm_top1000_quantile_uniform_mean_p0p8_gp0p8_s0p15_ab0p8_log1p_days | GSE93957 | 62 | 0.8105 | 7.735 | 9.187 | 0.604 |  | 0 |  |
| cfg_026_lgbm_top1000_quantile_uniform_mean_p0p8_gp0p8_s0p15_ab0p8_log1p_days | GSE121141 | 81 | 0.3743 | 38.386 | 50.101 | -0.6822 |  | 19 | 79.761 |
| cfg_026_lgbm_top1000_quantile_uniform_mean_p0p8_gp0p8_s0p15_ab0p8_log1p_days | GSE60012 | 152 | 0.7717 | 7.566 | 9.908 | -0.3004 |  | 0 |  |
| cfg_026_lgbm_top1000_quantile_uniform_mean_p0p8_gp0p8_s0p15_ab0p8_log1p_days | GSE213628 | 120 | 0.5895 | 24.06 | 29.2 | 0.3401 |  | 42 | 29.359 |
| cfg_006_lgbm_top1000_robust_median_p0p8_gp0p8_s0p15_ab0p8_log1p_days | GSE120137 | 549 | 0.6009 | 26.238 | 33.197 | -0.3209 |  | 0 |  |
| cfg_006_lgbm_top1000_robust_median_p0p8_gp0p8_s0p15_ab0p8_log1p_days | GSE80672 | 255 | 0.7541 | 29.601 | 35.83 | 0.4008 | 0.7843 | 74 | 51.606 |
| cfg_006_lgbm_top1000_robust_median_p0p8_gp0p8_s0p15_ab0p8_log1p_days | GSE93957 | 62 | 0.8601 | 7.362 | 9.127 | 0.6092 |  | 0 |  |
| cfg_006_lgbm_top1000_robust_median_p0p8_gp0p8_s0p15_ab0p8_log1p_days | GSE121141 | 81 | 0.3869 | 38.01 | 49.919 | -0.67 |  | 19 | 79.484 |
| cfg_006_lgbm_top1000_robust_median_p0p8_gp0p8_s0p15_ab0p8_log1p_days | GSE60012 | 152 | 0.7763 | 7.914 | 10.561 | -0.4773 |  | 0 |  |
| cfg_006_lgbm_top1000_robust_median_p0p8_gp0p8_s0p15_ab0p8_log1p_days | GSE213628 | 120 | 0.5966 | 23.669 | 28.854 | 0.3556 |  | 42 | 31.714 |

## Outputs

- `results/baseline_param_search_v18/v18_baseline_param_search_summary.csv`
- `results/baseline_param_search_v18/v18_lodo_summary.csv`
- `results/baseline_param_search_v18/v18_baseline_param_search_summary.json`

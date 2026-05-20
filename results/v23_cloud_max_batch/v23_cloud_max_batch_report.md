# v23 Cloud Max Batch Report

Date: 2026-05-20T22:27:43Z

## Scope

Robust cloud batch using the approved RTX 5090 host. The run expands v22 age DL beyond residual MLP width sweeps and runs CPU sidecar biological-signal screens without touching raw FASTQ.

## Execution

- host: `autodl-container-31eb41bf76-0d5d2e15`
- python: `/root/autodl-tmp/venvs/mmclock/bin/python`
- matrix: `/root/autodl-tmp/mouse_methyl_work/results/multidataset_v8_3_ablation/all6/all_rrbs_region_matrix_5kb.parquet`
- output: `results/v23_cloud_max_batch`
- group configs requested: `72`
- group/lodo GPU workers: `4` / `4`

## Selected Candidate

- config: `v23_007_res_mlp_top1500_robust_w128_d2_do0p1_lr0p0003`
- architecture: `res_mlp`
- GroupKFold r/MAE/R2: `0.7088` / `20.361` / `0.3851`
- LODO mean/worst MAE: `17.464` / `29.15`
- LODO CR AUC: `0.8159`
- old104 weighted MAE: `24.523`
- beats v20 on both LODO mean MAE and CR AUC: `False`

## Comparators

- v20: `{'config_id': 'cfg_241_res_mlp_top1000_robust_w256_d3_do0p1_lr0p0003', 'group_mae_weeks': 20.318, 'lodo_mean_mae_weeks': 18.109, 'lodo_worst_mae_weeks': 30.915, 'lodo_cr_auc': 0.8459, 'old_104w_weighted_mae_weeks': 27.07}`
- v22: `{'config_id': 'v22_007_res_mlp_top1000_robust_w192_d3_do0p1_lr0p0003', 'group_mae_weeks': 20.45, 'lodo_mean_mae_weeks': 18.81, 'lodo_worst_mae_weeks': 27.043, 'lodo_cr_auc': 0.8349, 'old_104w_weighted_mae_weeks': 25.673}`

## Top DL Configs

| config_id | status | architecture | n_feature_prefilter | preprocess | width | depth | dropout | lr | pearson_r | mae_weeks | r2 | cr_detection_auc | lodo_mean_mae_weeks | lodo_worst_mae_weeks | lodo_cr_auc | old_104w_weighted_mae_weeks | final_score |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| v23_007_res_mlp_top1500_robust_w128_d2_do0p1_lr0p0003 | skipped_existing | res_mlp | 1500 | robust | 128 | 2 | 0.1 | 0.0003 | 0.7088 | 20.361 | 0.3851 | 0.8071 | 17.464 | 29.15 | 0.8159 | 24.523 | 1.280351 |
| v23_004_res_mlp_top1000_robust_w192_d3_do0p1_lr0p0003 | skipped_existing | res_mlp | 1000 | robust | 192 | 3 | 0.1 | 0.0003 | 0.7099 | 20.45 | 0.4004 | 0.7832 | 18.81 | 27.043 | 0.8349 | 25.673 | 1.2780550000000002 |
| v23_046_res_mlp_top2000_robust_w192_d3_do0p15_lr0p0001 | skipped_existing | res_mlp | 2000 | robust | 192 | 3 | 0.15 | 0.0001 | 0.7071 | 20.32 | 0.3972 | 0.7218 | 19.449 | 29.898 | 0.8024 | 33.181 | 1.246631 |
| v23_012_res_mlp_top1500_robust_w256_d3_do0p1_lr0p0003 | skipped_existing | res_mlp | 1500 | robust | 256 | 3 | 0.1 | 0.0003 | 0.6869 | 20.942 | 0.3432 | 0.7502 | 18.448 | 28.854 | 0.8114 | 26.274 | 1.2444600000000001 |
| v23_013_res_mlp_top2000_robust_w128_d2_do0p1_lr0p0003 | skipped_existing | res_mlp | 2000 | robust | 128 | 2 | 0.1 | 0.0003 | 0.7025 | 20.602 | 0.3883 | 0.7123 | 18.603 | 28.305 | 0.7934 | 28.451 | 1.240645 |
| v23_037_res_mlp_top1000_robust_w192_d3_do0p05_lr0p0003 | skipped_existing | res_mlp | 1000 | robust | 192 | 3 | 0.05 | 0.0003 | 0.6969 | 21.588 | 0.347 | 0.7731 | 19.222 | 31.086 | 0.8761 | 33.942 | 1.240364 |
| v23_014_res_mlp_top2000_robust_w128_d3_do0p1_lr0p0003 | skipped_existing | res_mlp | 2000 | robust | 128 | 3 | 0.1 | 0.0003 | 0.7066 | 20.49 | 0.4161 | 0.7355 | 19.179 | 30.515 | 0.7399 | 34.865 | 1.238289 |
| v23_040_res_mlp_top1000_robust_w192_d3_do0p15_lr0p0001 | skipped_existing | res_mlp | 1000 | robust | 192 | 3 | 0.15 | 0.0001 | 0.6842 | 21.313 | 0.3665 | 0.7375 | 20.268 | 30.796 | 0.8494 | 27.855 | 1.224414 |
| v23_016_res_mlp_top2000_robust_w192_d3_do0p1_lr0p0003 | skipped_existing | res_mlp | 2000 | robust | 192 | 3 | 0.1 | 0.0003 | 0.6916 | 21.332 | 0.3616 | 0.7074 |  |  |  |  | 0.604875 |
| v23_015_res_mlp_top2000_robust_w192_d2_do0p1_lr0p0003 | skipped_existing | res_mlp | 2000 | robust | 192 | 2 | 0.1 | 0.0003 | 0.6855 | 21.869 | 0.3254 | 0.7748 |  |  |  |  | 0.604154 |
| v23_009_res_mlp_top1500_robust_w192_d2_do0p1_lr0p0003 | skipped_existing | res_mlp | 1500 | robust | 192 | 2 | 0.1 | 0.0003 | 0.6956 | 21.585 | 0.3438 | 0.7222 |  |  |  |  | 0.603779 |
| v23_006_res_mlp_top1000_robust_w256_d3_do0p1_lr0p0003 | skipped_existing | res_mlp | 1000 | robust | 256 | 3 | 0.1 | 0.0003 | 0.6819 | 21.495 | 0.3424 | 0.7488 |  |  |  |  | 0.603629 |
| v23_044_res_mlp_top2000_robust_w192_d3_do0p05_lr0p0001 | skipped_existing | res_mlp | 2000 | robust | 192 | 3 | 0.05 | 0.0001 | 0.7027 | 21.74 | 0.3248 | 0.7273 |  |  |  |  | 0.602125 |
| v23_038_res_mlp_top1000_robust_w192_d3_do0p05_lr0p0001 | skipped_existing | res_mlp | 1000 | robust | 192 | 3 | 0.05 | 0.0001 | 0.6905 | 22.022 | 0.3008 | 0.7811 |  |  |  |  | 0.600623 |
| v23_045_res_mlp_top2000_robust_w192_d3_do0p15_lr0p0003 | skipped_existing | res_mlp | 2000 | robust | 192 | 3 | 0.15 | 0.0003 | 0.686 | 21.302 | 0.3471 | 0.7109 |  |  |  |  | 0.600119 |
| v23_043_res_mlp_top2000_robust_w192_d3_do0p05_lr0p0003 | skipped_existing | res_mlp | 2000 | robust | 192 | 3 | 0.05 | 0.0003 | 0.6888 | 21.688 | 0.3541 | 0.6771 |  |  |  |  | 0.595307 |

## Random-Label Sanity

| config_id | status | pearson_r | mae_weeks | r2 | cr_detection_auc |
| --- | --- | --- | --- | --- | --- |
| v23_004_res_mlp_top1000_robust_w192_d3_do0p1_lr0p0003 | completed | 0.1672 | 32.546 | -0.2617 | 0.6319 |
| v23_007_res_mlp_top1500_robust_w128_d2_do0p1_lr0p0003 | completed | 0.2451 | 34.234 | -0.4088 | 0.3739 |

## Biological Sidecars

- sidecar jobs: `4`
| sidecar | target_id | model_name | cv_strategy | n_samples | balanced_accuracy | roc_auc | mae_weeks | pearson_r |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| auto_robust2000 | cr_vs_control | logistic | StratifiedKFold_within_scope | 255 | 0.9464 | 0.9851 |  |  |
| rf_classifiers_1000 | cr_vs_control | rf | StratifiedKFold_within_scope | 255 | 0.7165 | 0.9762 |  |  |
| lgbm_classifiers_2000 | cr_vs_control | lgbm | StratifiedKFold_within_scope | 255 | 0.7901 | 0.9742 |  |  |
| auto_robust2000 | sex_binary | logistic | GroupKFold_dataset_batch_class_covered | 1037 | 0.6963 | 0.7591 |  |  |
| rf_classifiers_1000 | sex_binary | rf | GroupKFold_dataset_batch_class_covered | 1037 | 0.7021 | 0.7563 |  |  |
| lgbm_classifiers_2000 | sex_binary | lgbm | GroupKFold_dataset_batch_class_covered | 1037 | 0.7204 | 0.7214 |  |  |
| auto_robust2000 | old_castration_vs_control | logistic | StratifiedKFold_within_scope | 20 | 0.7 | 0.68 |  |  |
| rf_classifiers_1000 | old_castration_vs_control | rf | StratifiedKFold_within_scope | 20 | 0.55 | 0.55 |  |  |
| lgbm_classifiers_2000 | old_castration_vs_control | lgbm | StratifiedKFold_within_scope | 20 | 0.5 | 0.5 |  |  |
| auto_robust2000 | tissue_multiclass | logistic | StratifiedKFold_within_scope | 1193 | 0.9948 |  |  |  |
| rf_classifiers_1000 | tissue_multiclass | rf | StratifiedKFold_within_scope | 1193 | 0.9867 |  |  |  |
| lgbm_classifiers_2000 | tissue_multiclass | lgbm | StratifiedKFold_within_scope | 1193 | 0.9817 |  |  |  |
| lgbm_classifiers_2000 | condition_family_multiclass | lgbm | StratifiedKFold_within_scope | 141 | 0.5179 |  |  |  |
| auto_robust2000 | condition_family_multiclass | logistic | StratifiedKFold_within_scope | 141 | 0.3913 |  |  |  |
| rf_classifiers_1000 | condition_family_multiclass | rf | StratifiedKFold_within_scope | 141 | 0.2143 |  |  |  |
| lgbm_age_2000 | age_weeks | lgbm | GroupKFold_dataset_batch | 1219 |  |  | 25.115 | 0.5704 |

## Feature Interpretation

`{'labels': {'cross_tissue_stable': 180, 'mixed_or_low_support': 68, 'tissue_specific': 2}, 'matrix': '/root/autodl-tmp/mouse_methyl_work/results/multidataset_v8_3_ablation/all6/all_rrbs_region_matrix_5kb.parquet', 'metadata': '/root/autodl-tmp/mouse_methyl_work/metadata/model_sample_metadata_v8.csv', 'n_clusters': 208, 'n_confounding_audit_rows': 2, 'n_regions_interpreted': 250, 'n_samples': 1219, 'timestamp': '2026-05-20T22:27:43.025445+00:00'}`

## Outputs

- `results/v23_cloud_max_batch/v23_group_summary.csv`
- `results/v23_cloud_max_batch/v23_lodo_summary.csv`
- `results/v23_cloud_max_batch/v23_final_summary.csv`
- `results/v23_cloud_max_batch/v23_random_label_summary.csv`
- `results/v23_cloud_max_batch/biology_sidecars`
- `results/v23_cloud_max_batch/feature_interpretation`

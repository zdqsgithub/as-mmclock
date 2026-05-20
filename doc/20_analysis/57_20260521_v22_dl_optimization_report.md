# v22 DL Optimization Report

Date: 2026-05-20T21:00:43.464462+00:00

## Scope

Narrow RTX 5090 DL optimization around the prior v20 res_mlp winner. The biological readout is age prediction plus CR detection from held-out age predictions; no CNN/Transformer expansion was used.

## Selected Candidate

- config: `v22_007_res_mlp_top1000_robust_w192_d3_do0p1_lr0p0003`
- GroupKFold r/MAE/R2: `0.7099` / `20.45` / `0.4004`
- LODO mean/worst MAE: `18.81` / `27.043`
- LODO CR AUC: `0.8349`
- old104 weighted MAE: `25.673`
- beats v20 on both LODO mean MAE and CR AUC: `False`

## v20 Comparator

`{'config_id': 'cfg_241_res_mlp_top1000_robust_w256_d3_do0p1_lr0p0003', 'group_mae_weeks': 20.318, 'lodo_mean_mae_weeks': 18.109, 'lodo_worst_mae_weeks': 30.915, 'lodo_cr_auc': 0.8459, 'old_104w_weighted_mae_weeks': 27.07}`

## Random-Label Sanity

`{'pearson_r': 0.1672, 'mae_weeks': 32.546, 'medae_weeks': 26.363, 'rmse_weeks': 41.373, 'r2': -0.2617, 'cr_detection_auc': 0.6319, 'cr_cohens_d': 0.4112, 'cross_dataset_mae': 32.793, 'n_samples': 1219}`

## Final All-Data Fit

`{'architecture': 'res_mlp', 'cuda_device': 'NVIDIA GeForce RTX 5090', 'torch_version': '2.11.0+cu128', 'n_samples': 1219, 'n_features_total': 65870, 'n_feature_prefilter': '1000', 'pearson_r': 0.9748, 'mae_weeks': 10.536, 'medae_weeks': 7.101, 'rmse_weeks': 14.985, 'r2': 0.8345, 'exec_time_sec': 16.1, 'warning': 'Apparent all-data fit metrics are optimistic and not held-out evidence.', 'n_features_selected': 1000, 'epochs_ran': 82}`

## Top Configs

| config_id | status | n_feature_prefilter | preprocess | width | depth | dropout | pearson_r | mae_weeks | r2 | cr_detection_auc | lodo_mean_mae_weeks | lodo_worst_mae_weeks | lodo_cr_auc | old_104w_weighted_mae_weeks | final_score |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| v22_007_res_mlp_top1000_robust_w192_d3_do0p1_lr0p0003 | completed | 1000 | robust | 192 | 3 | 0.1 | 0.7099 | 20.45 | 0.4004 | 0.7832 | 18.81 | 27.043 | 0.8349 | 25.673 | 1.2780550000000002 |
| v22_009_res_mlp_top1000_robust_w256_d3_do0p1_lr0p0003 | completed | 1000 | robust | 256 | 3 | 0.1 | 0.6819 | 21.495 | 0.3424 | 0.7488 | 18.276 | 28.696 | 0.8156 | 28.116 | 1.2351079999999999 |
| v22_011_res_mlp_top1000_robust_w384_d3_do0p1_lr0p0003 | completed | 1000 | robust | 384 | 3 | 0.1 | 0.6635 | 22.078 | 0.3207 | 0.7325 | 19.442 | 29.186 | 0.8222 | 29.166 | 1.210638 |
| v22_008_res_mlp_top1000_robust_w192_d4_do0p1_lr0p0003 | completed | 1000 | robust | 192 | 4 | 0.1 | 0.6719 | 21.619 | 0.3306 | 0.7569 | 19.654 | 31.892 | 0.8132 | 31.616 | 1.209181 |
| v22_010_res_mlp_top1000_robust_w256_d4_do0p1_lr0p0003 | completed | 1000 | robust | 256 | 4 | 0.1 | 0.6655 | 21.986 | 0.3196 | 0.6892 |  |  |  |  | 0.579725 |
| v22_005_res_mlp_top500_robust_w384_d3_do0p1_lr0p0003 | completed | 500 | robust | 384 | 3 | 0.1 | 0.6562 | 22.716 | 0.2751 | 0.7806 |  |  |  |  | 0.579472 |
| v22_006_res_mlp_top500_robust_w384_d4_do0p1_lr0p0003 | completed | 500 | robust | 384 | 4 | 0.1 | 0.6651 | 22.724 | 0.2849 | 0.7059 |  |  |  |  | 0.571553 |
| v22_003_res_mlp_top500_robust_w256_d3_do0p1_lr0p0003 | completed | 500 | robust | 256 | 3 | 0.1 | 0.6456 | 23.106 | 0.2742 | 0.7646 |  |  |  |  | 0.57113 |
| v22_012_res_mlp_top1000_robust_w384_d4_do0p1_lr0p0003 | completed | 1000 | robust | 384 | 4 | 0.1 | 0.6447 | 22.794 | 0.2892 | 0.6672 |  |  |  |  | 0.558018 |
| v22_004_res_mlp_top500_robust_w256_d4_do0p1_lr0p0003 | completed | 500 | robust | 256 | 4 | 0.1 | 0.6147 | 24.159 | 0.2142 | 0.7757 |  |  |  |  | 0.543743 |
| v22_001_res_mlp_top500_robust_w192_d3_do0p1_lr0p0003 | completed | 500 | robust | 192 | 3 | 0.1 | 0.6298 | 23.438 | 0.2438 | 0.6902 |  |  |  |  | 0.543537 |
| v22_002_res_mlp_top500_robust_w192_d4_do0p1_lr0p0003 | completed | 500 | robust | 192 | 4 | 0.1 | 0.6046 | 24.647 | 0.2066 | 0.7374 |  |  |  |  | 0.5296 |

## LODO Rows

| config_id | dataset | n_samples | pearson_r | mae_weeks | r2 | cr_detection_auc | old_104w_n | old_104w_mae_weeks |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| v22_007_res_mlp_top1000_robust_w192_d3_do0p1_lr0p0003 | GSE120137 | 549 | 0.6321 | 21.093 | 0.0677 |  | 0 |  |
| v22_007_res_mlp_top1000_robust_w192_d3_do0p1_lr0p0003 | GSE80672 | 255 | 0.9165 | 15.426 | 0.8125 | 0.8349 | 74 | 14.15 |
| v22_007_res_mlp_top1000_robust_w192_d3_do0p1_lr0p0003 | GSE93957 | 62 | 0.7832 | 8.02 | 0.4829 |  | 0 |  |
| v22_007_res_mlp_top1000_robust_w192_d3_do0p1_lr0p0003 | GSE121141 | 81 | 0.6068 | 27.043 | 0.0794 |  | 19 | 58.685 |
| v22_007_res_mlp_top1000_robust_w192_d3_do0p1_lr0p0003 | GSE60012 | 152 | 0.7301 | 15.072 | -4.5777 |  | 0 |  |
| v22_007_res_mlp_top1000_robust_w192_d3_do0p1_lr0p0003 | GSE213628 | 120 | 0.575 | 26.205 | 0.2183 |  | 42 | 31.043 |
| v22_009_res_mlp_top1000_robust_w256_d3_do0p1_lr0p0003 | GSE120137 | 549 | 0.6361 | 21.529 | 0.02 |  | 0 |  |
| v22_009_res_mlp_top1000_robust_w256_d3_do0p1_lr0p0003 | GSE80672 | 255 | 0.8921 | 17.447 | 0.7605 | 0.8156 | 74 | 20.005 |
| v22_009_res_mlp_top1000_robust_w256_d3_do0p1_lr0p0003 | GSE93957 | 62 | 0.741 | 7.855 | 0.5345 |  | 0 |  |
| v22_009_res_mlp_top1000_robust_w256_d3_do0p1_lr0p0003 | GSE121141 | 81 | 0.5741 | 28.696 | -0.0533 |  | 19 | 63.871 |
| v22_009_res_mlp_top1000_robust_w256_d3_do0p1_lr0p0003 | GSE60012 | 152 | 0.7003 | 10.931 | -2.5654 |  | 0 |  |
| v22_009_res_mlp_top1000_robust_w256_d3_do0p1_lr0p0003 | GSE213628 | 120 | 0.6274 | 23.201 | 0.3569 |  | 42 | 26.231 |
| v22_008_res_mlp_top1000_robust_w192_d4_do0p1_lr0p0003 | GSE120137 | 549 | 0.615 | 23.342 | -0.1336 |  | 0 |  |
| v22_008_res_mlp_top1000_robust_w192_d4_do0p1_lr0p0003 | GSE80672 | 255 | 0.8939 | 16.32 | 0.7741 | 0.8132 | 74 | 18.407 |
| v22_008_res_mlp_top1000_robust_w192_d4_do0p1_lr0p0003 | GSE93957 | 62 | 0.7347 | 8.408 | 0.46 |  | 0 |  |
| v22_008_res_mlp_top1000_robust_w192_d4_do0p1_lr0p0003 | GSE121141 | 81 | 0.5409 | 31.892 | -0.2651 |  | 19 | 70.054 |
| v22_008_res_mlp_top1000_robust_w192_d4_do0p1_lr0p0003 | GSE60012 | 152 | 0.7312 | 12.835 | -3.0452 |  | 0 |  |
| v22_008_res_mlp_top1000_robust_w192_d4_do0p1_lr0p0003 | GSE213628 | 120 | 0.5739 | 25.127 | 0.2659 |  | 42 | 37.501 |
| v22_011_res_mlp_top1000_robust_w384_d3_do0p1_lr0p0003 | GSE120137 | 549 | 0.6323 | 22.904 | -0.114 |  | 0 |  |
| v22_011_res_mlp_top1000_robust_w384_d3_do0p1_lr0p0003 | GSE80672 | 255 | 0.8502 | 16.548 | 0.7171 | 0.8222 | 74 | 20.178 |
| v22_011_res_mlp_top1000_robust_w384_d3_do0p1_lr0p0003 | GSE93957 | 62 | 0.7817 | 7.249 | 0.6048 |  | 0 |  |
| v22_011_res_mlp_top1000_robust_w384_d3_do0p1_lr0p0003 | GSE121141 | 81 | 0.5626 | 29.186 | -0.066 |  | 19 | 62.909 |
| v22_011_res_mlp_top1000_robust_w384_d3_do0p1_lr0p0003 | GSE60012 | 152 | 0.6919 | 17.559 | -6.2813 |  | 0 |  |
| v22_011_res_mlp_top1000_robust_w384_d3_do0p1_lr0p0003 | GSE213628 | 120 | 0.6359 | 23.207 | 0.3308 |  | 42 | 29.737 |

## Outputs

- `results/v22_dl_optimization/v22_dl_summary.csv`
- `results/v22_dl_optimization/v22_lodo_summary.csv`
- `results/v22_dl_optimization/random_label_sanity`
- `results/v22_dl_optimization/final_all_data_model`

# v8.3 Matrix Ablation 与 Tissue-Shared Diagnostics 报告

Date: 2026-05-18

## Summary

v8.3 比较了 core4、加入 GSE213628、加入 GSE60012、all6 四种 strict matrix variant，并只运行两个锁定配置。没有执行 autoresearch、deep learning 或新下载。

- Autoresearch-ready variants: 0
- Best GroupKFold MAE: 23.767w
- Best GSE121141 104w+ MAE: 75.386w

## Matrix Variants

| variant | datasets | status | n_regions | n_samples | blocker_reason |
| --- | --- | --- | --- | --- | --- |
| core4 | GSE120137,GSE80672,GSE93957,GSE121141 | completed | 72079 | 947 |  |
| core4_plus_gse213628 | GSE120137,GSE80672,GSE93957,GSE121141,GSE213628 | completed | 70302 | 1067 |  |
| core4_plus_gse60012 | GSE120137,GSE80672,GSE93957,GSE121141,GSE60012 | completed | 67435 | 1099 |  |
| all6 | GSE120137,GSE80672,GSE93957,GSE121141,GSE60012,GSE213628 | completed | 65870 | 1219 |  |

## Variant Summary

| variant | config | groupkfold_mae_weeks | groupkfold_r | gse121141_mae_weeks | gse121141_old_104w_mae_weeks | gse80672_mae_weeks | gse80672_cr_auc | gse80672_shuffled_cr_auc | random_label_r | random_label_pass | gse80672_shuffled_cr_pass | autoresearch_ready |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| core4 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | 23.767 | 0.617 | 38.033 | 75.386 | 26.706 | 0.8072 | 0.5284 | -0.0215 | True | True | False |
| core4 | 02_v74_lgbm_robust_p095_top1000_shift015 | 23.793 | 0.624 | 39.239 | 77.948 | 26.543 | 0.7649 | 0.4802 |  | False | True | False |
| core4_plus_gse213628 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | 25.585 | 0.6299 | 35.433 | 75.648 | 25.257 | 0.7997 | 0.5392 | -0.1139 | True | True | False |
| core4_plus_gse213628 | 02_v74_lgbm_robust_p095_top1000_shift015 | 25.015 | 0.6414 | 35.716 | 76.934 | 27.274 | 0.8744 | 0.527 |  | False | True | False |
| core4_plus_gse60012 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | 26.652 | 0.5618 | 40.919 | 82.21 | 34.081 | 0.7216 | 0.5048 | 0.135 | True | True | False |
| core4_plus_gse60012 | 02_v74_lgbm_robust_p095_top1000_shift015 | 26.025 | 0.582 | 40.766 | 82.888 | 31.939 | 0.6962 | 0.5233 |  | False | True | False |
| all6 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | 24.843 | 0.6253 | 37.499 | 78.09 | 28.86 | 0.7797 | 0.5137 | 0.0401 | True | True | False |
| all6 | 02_v74_lgbm_robust_p095_top1000_shift015 | 25.358 | 0.6137 | 37.485 | 79.921 | 30.092 | 0.7424 | 0.5083 |  | False | True | False |

## LODO Summary

| test_dataset | n_samples | pearson_r | pearson_pval | mae_weeks | rmse_weeks | r2 | old_104w_n | old_104w_mae_weeks | variant | config | cr_detection_auc |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE120137 | 549 | 0.5038 | 1.1392583470073409e-36 | 21.761 | 26.797 | 0.1393 | 0 |  | core4 | 01_v75_lgbm_quantile_p08_top1000_agebin08 |  |
| GSE80672 | 255 | 0.7303 | 9.295257883002063e-44 | 26.706 | 33.2 | 0.4856 | 74 | 41.463 | core4 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | 0.8072 |
| GSE93957 | 62 | 0.8124 | 1.1045819962953267e-15 | 10.81 | 12.457 | 0.2719 | 0 |  | core4 | 01_v75_lgbm_quantile_p08_top1000_agebin08 |  |
| GSE121141 | 81 | 0.358 | 0.0010342320551891748 | 38.033 | 48.778 | -0.5945 | 19 | 75.386 | core4 | 01_v75_lgbm_quantile_p08_top1000_agebin08 |  |
| GSE120137 | 549 | 0.537 | 2.471756905010821e-42 | 26.711 | 33.909 | -0.3781 | 0 |  | core4_plus_gse213628 | 01_v75_lgbm_quantile_p08_top1000_agebin08 |  |
| GSE80672 | 255 | 0.8341 | 2.662607777850446e-67 | 25.257 | 32.17 | 0.517 | 74 | 47.037 | core4_plus_gse213628 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | 0.7997 |
| GSE93957 | 62 | 0.8289 | 8.955328335891764e-17 | 8.47 | 9.981 | 0.5326 | 0 |  | core4_plus_gse213628 | 01_v75_lgbm_quantile_p08_top1000_agebin08 |  |
| GSE121141 | 81 | 0.3558 | 0.0011138908253111124 | 35.433 | 46.692 | -0.4611 | 19 | 75.648 | core4_plus_gse213628 | 01_v75_lgbm_quantile_p08_top1000_agebin08 |  |
| GSE213628 | 120 | 0.6367 | 5.4362663805845895e-15 | 23.326 | 27.849 | 0.3997 | 42 | 30.135 | core4_plus_gse213628 | 01_v75_lgbm_quantile_p08_top1000_agebin08 |  |
| GSE120137 | 549 | 0.4993 | 6.009057511920623e-36 | 28.914 | 36.325 | -0.5815 | 0 |  | core4_plus_gse60012 | 01_v75_lgbm_quantile_p08_top1000_agebin08 |  |
| GSE80672 | 255 | 0.5871 | 5.066072141156998e-25 | 34.081 | 40.834 | 0.2218 | 74 | 58.68 | core4_plus_gse60012 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | 0.7216 |
| GSE93957 | 62 | 0.8233 | 2.1622563560555291e-16 | 7.596 | 9.665 | 0.5618 | 0 |  | core4_plus_gse60012 | 01_v75_lgbm_quantile_p08_top1000_agebin08 |  |
| GSE121141 | 81 | 0.3718 | 0.0006313939459531661 | 40.919 | 52.005 | -0.8125 | 19 | 82.21 | core4_plus_gse60012 | 01_v75_lgbm_quantile_p08_top1000_agebin08 |  |
| GSE60012 | 152 | 0.7034 | 5.286090261962145e-24 | 6.186 | 8.419 | 0.0612 | 0 |  | core4_plus_gse60012 | 01_v75_lgbm_quantile_p08_top1000_agebin08 |  |
| GSE120137 | 549 | 0.5862 | 5.5606184545284015e-52 | 27.102 | 34.113 | -0.3948 | 0 |  | all6 | 01_v75_lgbm_quantile_p08_top1000_agebin08 |  |
| GSE80672 | 255 | 0.7652 | 2.700169647032373e-50 | 28.86 | 35.209 | 0.4214 | 74 | 50.319 | all6 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | 0.7797 |
| GSE93957 | 62 | 0.8418 | 1.0370166689990876e-17 | 7.902 | 9.303 | 0.5939 | 0 |  | all6 | 01_v75_lgbm_quantile_p08_top1000_agebin08 |  |
| GSE121141 | 81 | 0.3922 | 0.00029357157712094084 | 37.499 | 49.07 | -0.6137 | 19 | 78.09 | all6 | 01_v75_lgbm_quantile_p08_top1000_agebin08 |  |
| GSE60012 | 152 | 0.7851 | 5.062692398278156e-33 | 7.558 | 9.769 | -0.264 | 0 |  | all6 | 01_v75_lgbm_quantile_p08_top1000_agebin08 |  |
| GSE213628 | 120 | 0.6067 | 2.069860258385904e-13 | 23.359 | 28.684 | 0.3632 | 42 | 28.375 | all6 | 01_v75_lgbm_quantile_p08_top1000_agebin08 |  |

## GSE121141 104w+ Tissue Support

| variant | config | tissue | n_samples | mae_weeks | train_same_tissue_n | train_same_tissue_age_max | outside_range_rate |
| --- | --- | --- | --- | --- | --- | --- | --- |
| all6 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | brain_cortex | 5 | 92.44202306208452 | 16.0 | 41.0 | 1.0 |
| all6 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | heart | 5 | 100.4449531148728 | 37.0 | 121.68 | 1.0 |
| all6 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | liver | 4 | 40.98668340113731 | 181.0 | 104.297 | 1.0 |
| all6 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | lung | 5 | 71.0643709797526 | 79.0 | 104.297 | 1.0 |
| all6 | 02_v74_lgbm_robust_p095_top1000_shift015 | brain_cortex | 5 | 91.4891212645123 | 16.0 | 41.0 | 1.0 |
| all6 | 02_v74_lgbm_robust_p095_top1000_shift015 | heart | 5 | 98.85422091173 | 37.0 | 121.68 | 1.0 |
| all6 | 02_v74_lgbm_robust_p095_top1000_shift015 | liver | 4 | 43.81549429736222 | 181.0 | 104.297 | 1.0 |
| all6 | 02_v74_lgbm_robust_p095_top1000_shift015 | lung | 5 | 78.30371253881901 | 79.0 | 104.297 | 1.0 |
| core4 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | brain_cortex | 5 | 86.10802925284693 | 16.0 | 41.0 | 1.0 |
| core4 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | heart | 5 | 96.50110477566723 | 15.0 | 41.0 | 1.0 |
| core4 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | liver | 4 | 35.647327538540196 | 75.0 | 86.914 | 1.0 |
| core4 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | lung | 5 | 75.34159976127074 | 73.0 | 86.914 | 1.0 |
| core4 | 02_v74_lgbm_robust_p095_top1000_shift015 | brain_cortex | 5 | 87.77775372699669 | 16.0 | 41.0 | 1.0 |
| core4 | 02_v74_lgbm_robust_p095_top1000_shift015 | heart | 5 | 102.04474671586323 | 15.0 | 41.0 | 1.0 |
| core4 | 02_v74_lgbm_robust_p095_top1000_shift015 | liver | 4 | 39.28111225875671 | 75.0 | 86.914 | 1.0 |
| core4 | 02_v74_lgbm_robust_p095_top1000_shift015 | lung | 5 | 74.95604866261661 | 73.0 | 86.914 | 1.0 |
| core4_plus_gse213628 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | brain_cortex | 5 | 87.58990424647058 | 16.0 | 41.0 | 1.0 |
| core4_plus_gse213628 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | heart | 5 | 97.69546478076958 | 37.0 | 121.68 | 1.0 |
| core4_plus_gse213628 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | liver | 4 | 38.59922338943153 | 81.0 | 104.297 | 1.0 |
| core4_plus_gse213628 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | lung | 5 | 71.29952151580252 | 79.0 | 104.297 | 1.0 |
| core4_plus_gse213628 | 02_v74_lgbm_robust_p095_top1000_shift015 | brain_cortex | 5 | 89.66803931713608 | 16.0 | 41.0 | 1.0 |
| core4_plus_gse213628 | 02_v74_lgbm_robust_p095_top1000_shift015 | heart | 5 | 95.90005296163227 | 37.0 | 121.68 | 1.0 |
| core4_plus_gse213628 | 02_v74_lgbm_robust_p095_top1000_shift015 | liver | 4 | 42.63744698206679 | 81.0 | 104.297 | 1.0 |
| core4_plus_gse213628 | 02_v74_lgbm_robust_p095_top1000_shift015 | lung | 5 | 72.67219688176002 | 79.0 | 104.297 | 1.0 |
| core4_plus_gse60012 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | brain_cortex | 5 | 89.90775356942667 | 16.0 | 41.0 | 1.0 |
| core4_plus_gse60012 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | heart | 5 | 100.5593817244841 | 15.0 | 41.0 | 1.0 |
| core4_plus_gse60012 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | liver | 4 | 55.99245764597579 | 175.0 | 86.914 | 1.0 |
| core4_plus_gse60012 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | lung | 5 | 77.13533045942336 | 73.0 | 86.914 | 1.0 |
| core4_plus_gse60012 | 02_v74_lgbm_robust_p095_top1000_shift015 | brain_cortex | 5 | 88.85903148793827 | 16.0 | 41.0 | 1.0 |
| core4_plus_gse60012 | 02_v74_lgbm_robust_p095_top1000_shift015 | heart | 5 | 103.17416120060265 | 15.0 | 41.0 | 1.0 |
| core4_plus_gse60012 | 02_v74_lgbm_robust_p095_top1000_shift015 | liver | 4 | 56.107176567047865 | 175.0 | 86.914 | 1.0 |
| core4_plus_gse60012 | 02_v74_lgbm_robust_p095_top1000_shift015 | lung | 5 | 78.05624068206468 | 73.0 | 86.914 | 1.0 |

## Decision

No variant met all v8.3 acceptance criteria. Do not start v8.3.1 autoresearch. Use the ablation outputs to decide whether GSE60012/GSE213628 should remain auxiliary and prioritize GSE121141-like same-tissue old-age data or finer tissue/schema harmonization.

## Files

- `results/benchmark_v8_3_ablation/variant_summary.tsv`
- `results/benchmark_v8_3_ablation/lodo_summary_by_variant.csv`
- `results/benchmark_v8_3_ablation/gse121141_tissue_age_support.csv`
- `results/benchmark_v8_3_ablation/residual_by_dataset_tissue_agebin.csv`

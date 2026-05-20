# v17.1 Demo ML Narrow Optimization Report

Date: 2026-05-20T15:51:42.569867+00:00

## Scope

This is a narrow, high-confidence, demonstration-only ML optimization. It
does not change the v13 Route C model scope and is not broad autoresearch.

## Selected Configuration

- config: `c03_q_p08_top500_agebin08`
- preprocess: `quantile_uniform`
- top regions: `500`
- support-covered MAE: `19.888` weeks
- GroupKFold MAE: `24.484` weeks
- random-label pass: `True`

## Config Comparison

| config | preprocess | top_n | groupkfold_mae_weeks | groupkfold_r | support_covered_mae_weeks | unsupported_mae_weeks | gse121141_old104_mae_weeks | cr_auc |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| c03_q_p08_top500_agebin08 | quantile_uniform | 500 | 24.484 | 0.6166 | 19.888 | 40.259 | 84.681 | 0.787 |
| c06_robust_p08_top1000_agebin08 | robust | 1000 | 24.758 | 0.6256 | 20.215 | 40.354 | 83.933 | 0.7644 |
| c01_q_p08_top1000_agebin08 | quantile_uniform | 1000 | 24.843 | 0.6253 | 20.394 | 40.118 | 84.611 | 0.7567 |
| c04_robust_p095_top1000_agebin08 | robust | 1000 | 25.358 | 0.6137 | 20.739 | 41.214 | 84.425 | 0.7151 |
| c05_robust_p095_top2000_agebin08 | robust | 2000 | 25.507 | 0.6334 | 20.914 | 41.272 | 84.026 | 0.7789 |
| c02_q_p08_top2000_agebin08 | quantile_uniform | 2000 | 25.631 | 0.634 | 21.22 | 40.775 | 85.399 | 0.7597 |

## Selected LODO

| dataset | n_samples | pearson_r | mae_weeks | rmse_weeks | r2 | cr_detection_auc | old_104w_n | old_104w_mae_weeks |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE120137 | 549 | 0.5992 | 25.205 | 31.925 | -0.2216 |  | 0 |  |
| GSE80672 | 255 | 0.7454 | 28.774 | 35.367 | 0.4162 | 0.7944 | 74 | 49.858 |
| GSE93957 | 62 | 0.7787 | 7.654 | 9.236 | 0.5998 |  | 0 |  |
| GSE121141 | 81 | 0.3283 | 40.111 | 51.87 | -0.803 |  | 19 | 82.702 |
| GSE60012 | 152 | 0.7904 | 9.926 | 12.506 | -1.0719 |  | 0 |  |
| GSE213628 | 120 | 0.5917 | 23.753 | 29.372 | 0.3323 |  | 42 | 28.381 |

## Decision

The optimization produced a demo-selected ML baseline, but the strict
held-out metrics remain far from publication-grade clock thresholds.
The result is useful for demonstration and experiment planning, not for
new full-lifespan old target-tissue claims.

## Guardrails

- no FASTQ download;
- no Bismark;
- no deep learning;
- no human clock CpG mapping;
- no dummy AUC;
- no broad autoresearch beyond the six pre-registered configs.

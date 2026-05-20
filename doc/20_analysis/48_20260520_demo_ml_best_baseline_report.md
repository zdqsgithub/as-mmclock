# v17 Demo ML Best Baseline Report

Date: 2026-05-20T01:35:03.343230+00:00

## Scope

This is a demonstration-only ML baseline. It does not change the accepted
v13 Route C model scope and does not claim full-lifespan old target-tissue
generalization.

## Locked Configuration

- model: `lgbm`
- preprocessing: `quantile_uniform`
- feature filters: presence `0.8`, group presence `0.8`, dataset mean shift `0.15`, age-bin presence `0.8`
- feature selector: fold-internal age correlation
- top regions: `1000`
- target transform: `log1p_days`
- matrix: `results/multidataset_v8_3_ablation/all6/all_rrbs_region_matrix_5kb.parquet`

## Main Metrics

- GroupKFold MAE: `24.843` weeks
- GroupKFold r: `0.6253`
- GroupKFold R2: `0.2091`
- Random-label sanity r: `0.0401`
- Random-label sanity MAE: `32.841` weeks

## Support-Covered Split

- support-covered: `{'n_samples': 944, 'pearson_r': 0.5873, 'mae_weeks': 20.394, 'medae_weeks': 14.736, 'rmse_weeks': 26.706, 'r2': 0.2009}`
- unsupported: `{'n_samples': 275, 'pearson_r': 0.5681, 'mae_weeks': 40.118, 'medae_weeks': 35.097, 'rmse_weeks': 48.043, 'r2': -0.2025}`
- GSE121141 old104+ stress: `{'n_samples': 19, 'pearson_r': 0.0, 'mae_weeks': 84.611, 'medae_weeks': 95.892, 'rmse_weeks': 89.634, 'r2': None}`

## Leave-One-Dataset-Out

| dataset | n_samples | pearson_r | mae_weeks | rmse_weeks | r2 | cr_detection_auc | old_104w_n | old_104w_mae_weeks |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE120137 | 549 | 0.5862 | 27.102 | 34.113 | -0.3948 |  | 0 |  |
| GSE80672 | 255 | 0.7652 | 28.86 | 35.209 | 0.4214 | 0.7797 | 74 | 50.319 |
| GSE93957 | 62 | 0.8418 | 7.902 | 9.303 | 0.5939 |  | 0 |  |
| GSE121141 | 81 | 0.3922 | 37.499 | 49.07 | -0.6137 |  | 19 | 78.09 |
| GSE60012 | 152 | 0.7851 | 7.558 | 9.769 | -0.264 |  | 0 |  |
| GSE213628 | 120 | 0.6067 | 23.359 | 28.684 | 0.3632 |  | 42 | 28.375 |

## Full-Data Demo Fit

The final all-data model is saved only for demonstration. Its apparent
training metrics are optimistic and must not be used as held-out evidence.

- apparent all-data MAE: `2.264` weeks
- apparent all-data r: `0.9946`
- selected features: `1000`

## Guardrails

- demo-only training was explicitly requested by the user;
- no autoresearch was run;
- no FASTQ download or Bismark was run;
- no human clock CpG mapping was used;
- no dummy AUC was generated;
- this result is for demonstration and experiment-planning only.

## Outputs

- `results/demo_ml_best_baseline_v17/groupkfold/predictions.csv`
- `results/demo_ml_best_baseline_v17/support_annotations.csv`
- `results/demo_ml_best_baseline_v17/lodo_summary.csv`
- `results/demo_ml_best_baseline_v17/final_all_data_model/final_all_data_lgbm_model.joblib`

# Autoresearch v2 Architecture Fix Report

> Date: 2026-05-17  
> Scope: Phase 0 GSE120137-only beta matrix, leakage-safe CV, standardized benchmark

## What Changed

- Metadata parsing was rebuilt from local GEO SOFT fields. Age coverage improved from 779/1187 to 1155/1187 samples.
- GSE80672 now has 255/255 age labels and CR/control intervention labels (32 CR, 223 control).
- Fold-internal age-correlation CpG selection replaced the previous whole-dataset prefilter, removing the known leakage path.
- Prediction files now use the benchmark schema: `age_days_true`, `age_days_pred`, `dataset_batch`, `tissue`, `intervention`.
- Benchmark metrics now report CR/rapamycin metrics only when real intervention labels exist in the evaluated prediction set.
- Autoresearch v2 writes to `results/autoresearch_v2/` and leaves v1 untouched.

## Validation Results

Metadata validation:

| Dataset | Samples | Age known | Non-control labels |
|---|---:|---:|---:|
| GSE120137 | 549 | 549 | 0 |
| GSE80672 | 255 | 255 | 32 |
| GSE121141 | 81 | 81 | 0 |
| GSE52266 | 40 | 40 | 30 |
| GSE60012 | 173 | 164 | 47 |
| GSE93957 | 62 | 62 | 0 |
| GSE80761 | 4 | 4 | 0 |
| GSE45361 | 23 | 0 | 0 |

Smoke tests:

| Config | Pearson r | MAE weeks | R2 |
|---|---:|---:|---:|
| ridge 10k mean | 0.8217 | 12.697 | 0.6467 |
| ridge 10k median | 0.8214 | 12.704 | 0.6464 |
| lgbm 5k median | 0.8332 | 13.327 | 0.6340 |
| elasticnet 50k mean | 0.8173 | 12.665 | 0.6390 |
| rf 5k mean | 0.7930 | 16.227 | 0.4956 |

Random-label sanity check:

| Config | Pearson r | MAE weeks | R2 |
|---|---:|---:|---:|
| randomized ridge 10k mean | 0.0160 | 26.305 | -0.2029 |

The random-label result collapses as expected, supporting that the repaired fold-internal feature selection removed the obvious leakage path.

## Autoresearch v2 Results

The fixed v2 search space contains 48 unique configurations (`4 model types x 6 feature counts x 2 imputations`). The run requested 100 attempt slots; all 48 unique configurations completed and the remaining 52 slots were recorded as `skipped_exhausted` rather than duplicated.

Best by composite score:

| Rank | Experiment | Model | CpGs | Imputation | Pearson r | MAE weeks | R2 |
|---:|---|---|---:|---|---:|---:|---:|
| 1 | exp_024_lgbm_1000_mean | lgbm | 1000 | mean | 0.8404 | 12.750 | 0.6516 |
| 2 | exp_003_lgbm_1000_median | lgbm | 1000 | median | 0.8361 | 12.870 | 0.6460 |
| 3 | exp_014_lgbm_2000_median | lgbm | 2000 | median | 0.8365 | 13.055 | 0.6385 |
| 4 | exp_037_lgbm_2000_mean | lgbm | 2000 | mean | 0.8343 | 12.994 | 0.6394 |
| 5 | exp_017_lgbm_5000_median | lgbm | 5000 | median | 0.8332 | 13.327 | 0.6340 |

Best by MAE:

| Rank | Experiment | Model | CpGs | Imputation | Pearson r | MAE weeks | R2 |
|---:|---|---|---:|---|---:|---:|---:|
| 1 | exp_043_elasticnet_2000_median | elasticnet | 2000 | median | 0.8129 | 12.471 | 0.6326 |
| 2 | exp_044_elasticnet_2000_mean | elasticnet | 2000 | mean | 0.8127 | 12.494 | 0.6324 |
| 3 | exp_025_elasticnet_20000_median | elasticnet | 20000 | median | 0.8173 | 12.550 | 0.6446 |
| 4 | exp_012_ridge_5000_median | ridge | 5000 | median | 0.8163 | 12.554 | 0.6396 |
| 5 | exp_005_elasticnet_20000_mean | elasticnet | 20000 | mean | 0.8178 | 12.555 | 0.6444 |

Best by model family:

| Model | Runs | Best r | Best MAE weeks | Best composite |
|---|---:|---:|---:|---:|
| lgbm | 12 | 0.8404 | 12.750 | 0.707563 |
| ridge | 12 | 0.8217 | 12.554 | 0.694906 |
| elasticnet | 12 | 0.8178 | 12.471 | 0.691965 |
| rf | 12 | 0.8148 | 15.148 | 0.661082 |

## Interpretation

The v1 best run reported `ridge 10k mean, r=0.8364, MAE=10.107w`, but that run used whole-dataset CpG prefiltering and a dummy `cr_auc=0.5`. After the v2 leakage fix, the best Pearson result is slightly higher (`lgbm 1000 mean, r=0.8404`) but MAE is still much worse than the publication target (`12.750w` vs target `<3.5w`). The lowest MAE is `12.471w`.

This supports the planned scientific decision: do not invest in MLP/CNN/Transformer yet. The next meaningful step is to build region-based features and complete multi-dataset beta matrices so the model can reduce tissue/batch sensitivity and enable real CR/rapamycin biological validation.

## Notes

- Current v2 training is still Phase 0 GSE120137-only. GSE80672 intervention labels are now available in metadata, but no GSE80672 beta matrix is available in the current model input.
- GSE80672 SOFT uses the source label `age (years)` with values that look biologically unusual for mice. The parser follows the specified unit rule for now; before using GSE80672 as a final biological validation set, the age unit should be cross-checked against the publication/sample table.

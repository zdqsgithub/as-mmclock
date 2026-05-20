# Embedding-Aware Clock v4 POC Report

> Date: 2026-05-18  
> Scope: Phase 0 GSE120137 5kb region matrix, CV-safe endogenous embeddings, feature interpretation, external foundation-model audit

## Decision

Embedding was implemented as both a training input and an interpretation layer. The benchmark path uses only project-native mouse embeddings fitted inside each CV fold. External MethylGPT/CpGPT integration remains a gated feasibility POC because the current project input is mouse RRBS 5kb region bins, not the human CpG/probe vocabularies used by those models.

References:
- MethylGPT: https://github.com/albert-ying/MethylGPT
- CpGPT: https://spacefrontiers.org/r/10.1101/2024.10.24.619766

## Implemented Artifacts

- `scripts/train/train_embedding_clock.py`: CV-safe SVD/NMF embedding clock trainer.
- `scripts/train/autoresearch_embedding_loop.py`: v4 embedding autoresearch loop with optional search filters.
- `scripts/audit/external_embedding_audit.py`: MethylGPT/CpGPT install/runtime/compatibility audit.
- Per-run outputs: `sample_embeddings.parquet`, `feature_embeddings.parquet`, `feature_interpretation.csv`, `cluster_summary.csv`, and `embedding_feature_report.md`.

Leakage controls in the trainer:

- Imputation is fitted on train folds only.
- SVD/NMF embedder is fitted on train folds only for metrics.
- Hybrid top-region selection is fitted on train folds only.
- Scaling and model fitting are train-fold only.
- Full-data embedding artifacts are labeled as interpretation-only and are not used for out-of-fold predictions.

## Smoke Results

| Config | Pearson r | MAE weeks | R2 | Interpretation |
|---|---:|---:|---:|---|
| svd32 + ridge | 0.6478 | 19.254 | 0.3132 | embedding-only linear is weak |
| svd64 + lgbm | 0.7752 | 14.173 | 0.5821 | embedding-only LGBM improves but underperforms v3 |
| nmf32 + elasticnet | 0.2048 | 25.583 | -0.1045 | NMF is slow and weak in this setup |
| svd64 + top500 regions + lgbm | 0.8594 | 11.106 | 0.7221 | best v4 POC result |
| randomized svd64 + top500 + lgbm | -0.0089 | 27.044 | -0.2696 | leakage sanity check passed |

The best smoke configuration improves over v3 best MAE `11.550w` by `0.444w`, but this is below the pre-set `>=1w` threshold for a larger v4 search.

## Autoresearch v4 POC

NMF was gated after smoke because it took several minutes for one run and produced poor accuracy. A focused SVD-only autoresearch POC ran 12 configurations from the v4 search interface.

Best SVD-only autoresearch runs:

| Rank | Experiment | Pearson r | MAE weeks | R2 |
|---:|---|---:|---:|---:|
| 1 | exp_008_svd32_lgbm_hybrid500_mean | 0.8596 | 11.288 | 0.7172 |
| 2 | exp_007_svd16_lgbm_hybrid500_mean | 0.8553 | 11.391 | 0.7097 |
| 3 | exp_009_svd64_lgbm_hybrid1000_median | 0.8551 | 11.447 | 0.7097 |

The POC confirms the pattern: embedding-only loses information, while embedding plus fold-selected raw regions is slightly better than v3 but not enough to justify large embedding search before adding datasets.

## Feature Interpretation

The best interpretation artifact is in `results/smoke_v4_embedding/svd64_top500_lgbm/`.

Top regions show stable fold selection but many have substantial tissue eta-squared values, meaning some apparent age signal is still mixed with tissue structure in the GSE120137-only matrix. Examples:

| Region | Combined importance | Age abs corr | Tissue eta2 | Fold stability | Top factor |
|---|---:|---:|---:|---:|---|
| chr2:18685000-18689999 | 47.8071 | 0.5464 | 0.1451 | 1.000 | emb_013 |
| chr10:45485000-45489999 | 44.9051 | 0.3942 | 0.1233 | 1.000 | emb_009 |
| chr11:21990000-21994999 | 43.6159 | 0.3880 | 0.4381 | 1.000 | emb_009 |
| chr2:147080000-147084999 | 42.6259 | 0.3427 | 0.6188 | 1.000 | emb_010 |

Interpretation conclusion: embedding is useful for grouping and ranking region clusters, but the current Phase 0 feature signals are not cleanly age-only. Multi-dataset validation is required before biological interpretation.

## External Model Audit

Audit output: `results/external_embedding_audit/`.

| Model | Installed | Compatibility | Decision |
|---|---:|---|---|
| MethylGPT | false | blocked_for_benchmark | not used in v4 benchmark |
| CpGPT | false | blocked_for_benchmark | not used in v4 benchmark |

Reason: current input is mouse 5kb RRBS regions; public MethylGPT/CpGPT paths require compatible methylation CpG/probe vocabularies and validated feature ordering. AS-DS-Ops forbids unvalidated human-to-mouse CpG transfer for model features.

## Next Step

Do not expand v4 embedding search yet. The best validated gain is under 1 week, and NMF is not promising. The next scientifically useful step is to build additional dataset matrices, starting with GSE80672 processed methylation data for true CR validation, then GSE93957/GSE121141/GSE60012 for GroupKFold cross-dataset evaluation.

Embedding should be retained as an interpretation layer and revisited after multi-dataset matrices exist.

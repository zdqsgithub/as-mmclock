# v8.4 Old Tissue Candidate and Calibration Diagnostics Report

Date: 2026-05-18

## Summary

v8.4 did not run autoresearch, deep learning, FASTQ/Bismark, or large GEO downloads. It screened official GEO/SOFT/filelist metadata for GSE121141-like old-age tissues and ran diagnostic-only calibration on existing v8.3 GSE121141 held-out predictions.

- Non-integrated P1 candidates found: 0
- Official web-screened accessions checked: GSE232547, GSE171236, GSE138368, GSE151541, GSE225166
- Core target tissues: brain_cortex, heart, lung; liver is reference-only.

## Candidate Ranking

| dataset | v8_4_score | v8_4_role | v8_4_recommended_action | target_old_counts | reference_old_tissues | max_age_by_target_tissue | processed_schema_guess | v8_4_blockers | title |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE213628 | 14 | p3_reference_tissue_not_target_gap | auxiliary_validation_only | heart:10;lung:2 | liver | heart:121.7;liver:104.3;lung:104.3 | bismark_cov_per_sample_tar | already_integrated | DNA methylation entropy as a measure of stem cell replication and aging [BiSulfite-seq] |
| GSE225166 | 11 | p2_single_tissue_or_context_limited | adapter_audit_only | lung:13;heart:11 | liver | heart:104.3;liver:104.3;lung:104.3 | bismark_cov_per_sample_tar | targeted_celltype_or_non_bulk_context | Predicting age in single cells and low coverage DNA methylation data [RRBS] |
| GSE213723 | 11 | p2_single_tissue_or_context_limited | adapter_audit_only | heart:10;lung:2 | liver | heart:121.7;liver:104.3;lung:104.3 | bismark_cov_per_sample_tar | superseries_use_subseries | DNA methylation entropy as a measure of stem cell replication and aging |
| GSE224442 | 6 | p3_reference_tissue_not_target_gap | auxiliary_validation_only |  | liver | liver:113.0 | bismark_cov_per_sample_tar | no_104w_plus_brain_cortex_heart_lung_samples;no_target_tissue_overlap | DNA methylation changes in mice induced by parabiosis and recovery [RRBS] |
| GSE92486 | 6 | p3_reference_tissue_not_target_gap | auxiliary_validation_only |  | liver | liver:113.0 | bismark_cov_per_sample_tar | no_104w_plus_brain_cortex_heart_lung_samples;no_target_tissue_overlap | Dietary restriction protects from age-associated DNA methylation and induces epigenetic reprogramming of lipid metabolism |
| GSE175410 | 5 | not_gse121141_like | do_not_download |  |  |  | bismark_cov_per_sample_tar | no_104w_plus_brain_cortex_heart_lung_samples;no_target_tissue_overlap | Exercise Mitigates Skeletal Muscle Epigenetic Aging |
| GSE266961 | 4 | not_gse121141_like | do_not_download |  |  |  | bismark_cov_per_sample_tar | no_104w_plus_brain_cortex_heart_lung_samples;no_target_tissue_overlap;max_exact_age_below_104w_or_missing | Novel enzyme-based reduced representation method for DNA methylation profiling with low inputs |
| GSE304754 | 4 | not_gse121141_like | do_not_download |  |  |  | unknown_or_non_methylation | processed_cov_not_detected;no_104w_plus_brain_cortex_heart_lung_samples;no_target_tissue_overlap | Multi-omics analysis highlights the link of aging-related cognitive decline with systemic inflammation and alterations of tissue-maintenance [RRBS] |
| GSE319289 | 4 | not_gse121141_like | do_not_download |  |  |  | bismark_cov_per_sample_tar | no_104w_plus_brain_cortex_heart_lung_samples;no_target_tissue_overlap;max_exact_age_below_104w_or_missing | The Age-Dependent Resident Myonuclear Multi-Omic Response to an Acute Skeletal Muscle Hypertrophic Stimulus in Mice |
| GSE129712 | 3 | not_gse121141_like | do_not_download |  |  |  | unknown_or_non_methylation | processed_cov_not_detected;no_104w_plus_brain_cortex_heart_lung_samples;no_target_tissue_overlap | Transcriptional and epigenetic landscape of the mouse intestine during aging identifies key molecular drivers of aging -associated dysfunctions and diseases [RRBS] |
| GSE286302 | 3 | not_gse121141_like | do_not_download |  |  |  | bismark_cov_per_sample_tar | no_104w_plus_brain_cortex_heart_lung_samples;no_target_tissue_overlap;max_exact_age_below_104w_or_missing;limited_sample_specific_age_parse | Restoring Circadian Rhythms in the Paraventricular Nucleus Reverses the Epigenetic Clock of Multi Tissues in Aging Mice. [RRBS] |
| GSE295059 | 3 | not_gse121141_like | do_not_download |  |  |  | bismark_cov_per_sample_tar | no_104w_plus_brain_cortex_heart_lung_samples;no_target_tissue_overlap;targeted_celltype_or_non_bulk_context | In Vitro Modeling of Age-Related Methylation Changes with Intestinal Organoids [RRBS] |
| GSE138368 | 2 | not_gse121141_like | do_not_download |  |  |  | unknown_or_non_methylation | processed_cov_not_detected;no_104w_plus_brain_cortex_heart_lung_samples;no_target_tissue_overlap;max_exact_age_below_104w_or_missing | Environmental enrichment restores a young DNA methylation landscape in the aged mouse hippocampus |
| GSE281602 | 2 | not_gse121141_like | do_not_download |  |  |  | bismark_cov_per_sample_tar | no_104w_plus_brain_cortex_heart_lung_samples;no_target_tissue_overlap;max_exact_age_below_104w_or_missing;limited_sample_specific_age_parse | Transcriptional profiling of aged and young mouse hearts [RRBS] |

## Official Web Screen Additions

| dataset | v8_4_role | v8_4_recommended_action | target_old_counts | processed_schema_guess | max_age_weeks | v8_4_blockers | geo_url |
| --- | --- | --- | --- | --- | --- | --- | --- |
| GSE225166 | p2_single_tissue_or_context_limited | adapter_audit_only | lung:13;heart:11 | bismark_cov_per_sample_tar | 104.297 | targeted_celltype_or_non_bulk_context | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225166 |
| GSE138368 | not_gse121141_like | do_not_download |  | unknown_or_non_methylation | 73.877 | processed_cov_not_detected;no_104w_plus_brain_cortex_heart_lung_samples;no_target_tissue_overlap;max_exact_age_below_104w_or_missing | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138368 |
| GSE171236 | not_gse121141_like | do_not_download |  | unknown_or_non_methylation | 73.877 | processed_cov_not_detected;no_104w_plus_brain_cortex_heart_lung_samples;no_target_tissue_overlap;max_exact_age_below_104w_or_missing | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171236 |
| GSE151541 | not_gse121141_like | do_not_download |  | bismark_cov_per_sample_tar |  | no_104w_plus_brain_cortex_heart_lung_samples;no_target_tissue_overlap;max_exact_age_below_104w_or_missing;limited_sample_specific_age_parse;targeted_celltype_or_non_bulk_context | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151541 |
| GSE232547 | not_gse121141_like | do_not_download |  | unknown_or_non_methylation |  | processed_cov_not_detected;no_104w_plus_brain_cortex_heart_lung_samples;no_target_tissue_overlap;max_exact_age_below_104w_or_missing;limited_sample_specific_age_parse;targeted_celltype_or_non_bulk_context | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE232547 |

## Shared-Tissue Age-Range Diagnostic

For GSE121141 104w+ samples, this table separates samples that are inside vs outside the same-tissue age range available in training. If old samples are outside range, model tuning is unlikely to fix the extrapolation problem by itself.

| variant | config | range_status | n | mae_weeks | median_ae_weeks | median_train_same_tissue_age_max |
| --- | --- | --- | --- | --- | --- | --- |
| all6 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | outside_same_tissue_age_range | 19 | 78.09 | 84.852 | 104.297 |
| all6 | 02_v74_lgbm_robust_p095_top1000_shift015 | outside_same_tissue_age_range | 19 | 79.921 | 89.175 | 104.297 |
| core4 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | outside_same_tissue_age_range | 19 | 75.386 | 90.41 | 41.0 |
| core4 | 02_v74_lgbm_robust_p095_top1000_shift015 | outside_same_tissue_age_range | 19 | 77.948 | 93.522 | 41.0 |
| core4_plus_gse213628 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | outside_same_tissue_age_range | 19 | 75.648 | 86.53 | 104.297 |
| core4_plus_gse213628 | 02_v74_lgbm_robust_p095_top1000_shift015 | outside_same_tissue_age_range | 19 | 76.934 | 87.973 | 104.297 |
| core4_plus_gse60012 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | outside_same_tissue_age_range | 19 | 82.21 | 96.038 | 41.0 |
| core4_plus_gse60012 | 02_v74_lgbm_robust_p095_top1000_shift015 | outside_same_tissue_age_range | 19 | 82.888 | 92.98 | 41.0 |

## Tissue-Specific Calibration Diagnostic

Calibration modes using GSE121141 labels are diagnostic-only. They test whether errors look like a tissue offset, not whether a deployable cross-dataset model has improved.

| variant | config | calibration_mode | subset | n | mae_weeks | median_ae_weeks | bias_weeks |
| --- | --- | --- | --- | --- | --- | --- | --- |
| all6 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | raw | all | 81 | 37.499 | 26.325 | 33.62 |
| all6 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | raw | old104plus | 19 | 78.09 | 84.852 | 78.09 |
| all6 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | global_under104_intercept | all | 81 | 29.638 | 24.136 | 14.124 |
| all6 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | global_under104_intercept | old104plus | 19 | 60.181 | 65.356 | 58.594 |
| all6 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | tissue_under104_intercept | all | 81 | 29.157 | 24.413 | 16.057 |
| all6 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | tissue_under104_intercept | old104plus | 19 | 59.966 | 70.488 | 59.966 |
| core4 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | raw | all | 81 | 38.033 | 32.241 | 31.916 |
| core4 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | raw | old104plus | 19 | 75.386 | 90.41 | 75.251 |
| core4 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | global_under104_intercept | all | 81 | 30.259 | 23.165 | 13.059 |
| core4 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | global_under104_intercept | old104plus | 19 | 58.895 | 71.554 | 56.394 |
| core4 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | tissue_under104_intercept | all | 81 | 29.493 | 23.758 | 14.012 |
| core4 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | tissue_under104_intercept | old104plus | 19 | 58.105 | 68.923 | 56.92 |
| core4_plus_gse213628 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | raw | all | 81 | 35.433 | 27.118 | 29.145 |
| core4_plus_gse213628 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | raw | old104plus | 19 | 75.648 | 86.53 | 75.648 |
| core4_plus_gse213628 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | global_under104_intercept | all | 81 | 30.966 | 22.433 | 17.794 |
| core4_plus_gse213628 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | global_under104_intercept | old104plus | 19 | 64.304 | 75.179 | 64.298 |
| core4_plus_gse213628 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | tissue_under104_intercept | all | 81 | 29.307 | 24.512 | 15.265 |
| core4_plus_gse213628 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | tissue_under104_intercept | old104plus | 19 | 61.138 | 65.872 | 61.138 |
| core4_plus_gse60012 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | raw | all | 81 | 40.919 | 34.633 | 37.519 |
| core4_plus_gse60012 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | raw | old104plus | 19 | 82.21 | 96.038 | 82.21 |
| core4_plus_gse60012 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | global_under104_intercept | all | 81 | 30.865 | 26.314 | 16.185 |
| core4_plus_gse60012 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | global_under104_intercept | old104plus | 19 | 60.876 | 74.705 | 60.876 |
| core4_plus_gse60012 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | tissue_under104_intercept | all | 81 | 30.513 | 25.716 | 15.392 |
| core4_plus_gse60012 | 01_v75_lgbm_quantile_p08_top1000_agebin08 | tissue_under104_intercept | old104plus | 19 | 59.859 | 72.3 | 59.859 |

## Decision

No new non-integrated P1 dataset was identified for old brain_cortex/heart/lung bulk RRBS with a clean processed COV schema. Do not start v8.4 downloads or autoresearch from this screen.

The next useful branch is a narrower official GEO refresh for cortex/heart/lung old-age methylation datasets, plus parser smoke only for candidates that have sample-specific age, bulk tissue, and processed COV/bedGraph-like methylation files. Current P2/P3 hits should remain auxiliary or method-specific validation data.

## Files

- `results/v8_4_old_tissue_diagnostics/old_tissue_candidate_ranking.csv`
- `results/v8_4_old_tissue_diagnostics/official_web_screened_inventory.csv`
- `results/v8_4_old_tissue_diagnostics/gse121141_shared_tissue_support_summary.csv`
- `results/v8_4_old_tissue_diagnostics/gse121141_tissue_calibration_summary.csv`
- `results/v8_4_old_tissue_diagnostics/gse121141_tissue_calibrated_predictions.csv`

## Official Source URLs

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121141
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213628
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE232547
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171236
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151541
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225166

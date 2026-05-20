# v8.5 Targeted GEO Refresh and GSE225166 Adapter Audit Report

Date: 2026-05-18

## Summary

v8.5 refreshed official GEO candidates using targeted E-Utils queries and audited GSE225166 without downloading large supplements or running training.

- Targeted search candidates preflighted: 38
- Direct mainline-ready new candidates: 0
- GSE225166 recommendation: auxiliary_adapter_smoke_only
- GSE225166 old core-target samples: 24 (cortex=0, heart=11, lung=13)

## Candidate Ranking

| dataset | v8_5_direct_mainline_ready | v8_4_role | v8_4_recommended_action | target_old_counts | max_age_by_target_tissue | processed_schema_guess | v8_4_blockers | title | geo_url |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE213628 | False | p3_reference_tissue_not_target_gap | auxiliary_validation_only | heart:10;lung:2 | heart:121.7;liver:104.3;lung:104.3 | bismark_cov_per_sample_tar | already_integrated | DNA methylation entropy as a measure of stem cell replication and aging [BiSulfite-seq] | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213628 |
| GSE225166 | False | p2_single_tissue_or_context_limited | adapter_audit_only | lung:13;heart:11 | heart:104.3;liver:104.3;lung:104.3 | bismark_cov_per_sample_tar | targeted_celltype_or_non_bulk_context | Predicting age in single cells and low coverage DNA methylation data [RRBS] | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225166 |
| GSE213723 | False | p2_single_tissue_or_context_limited | adapter_audit_only | heart:10;lung:2 | heart:121.7;liver:104.3;lung:104.3 | bismark_cov_per_sample_tar | superseries_use_subseries | DNA methylation entropy as a measure of stem cell replication and aging | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213723 |
| GSE83947 | False | p2_single_tissue_or_context_limited | adapter_audit_only | lung:180 | liver:108.6;lung:108.6 | unknown_or_non_methylation | processed_cov_not_detected;single_target_tissue_only | Cytosine modifications exhibit circadian oscillations that are involved in epigenetic diversity and aging [Bisulfite-Seq] | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE83947 |
| GSE80672 | False | not_gse121141_like | already_in_matrix_or_current_target | lung:6 | lung:130.4 | unknown_methylation_raw_tar | already_integrated;processed_cov_not_detected;single_target_tissue_only | Using DNA methylation profiling to evaluate biological age and longevity interventions. | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80672 |
| GSE134398 | False | p2_single_tissue_or_context_limited | adapter_audit_only | lung:8 | lung:104.3 | bismark_cov_per_sample_tar | single_target_tissue_only;superseries_use_subseries | The Aging Microenvironment Drives Alveolar Macrophage Changes in Aging | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134398 |
| GSE224442 | False | p3_reference_tissue_not_target_gap | auxiliary_validation_only |  | liver:113.0 | bismark_cov_per_sample_tar | no_104w_plus_brain_cortex_heart_lung_samples;no_target_tissue_overlap | DNA methylation changes in mice induced by parabiosis and recovery [RRBS] | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224442 |
| GSE92486 | False | p3_reference_tissue_not_target_gap | auxiliary_validation_only |  | liver:113.0 | bismark_cov_per_sample_tar | no_104w_plus_brain_cortex_heart_lung_samples;no_target_tissue_overlap | Dietary restriction protects from age-associated DNA methylation and induces epigenetic reprogramming of lipid metabolism | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92486 |
| GSE175410 | False | not_gse121141_like | do_not_download |  |  | bismark_cov_per_sample_tar | no_104w_plus_brain_cortex_heart_lung_samples;no_target_tissue_overlap | Exercise Mitigates Skeletal Muscle Epigenetic Aging | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE175410 |
| GSE153331 | False | not_gse121141_like | do_not_download |  |  | bismark_cov_per_sample_tar | no_104w_plus_brain_cortex_heart_lung_samples;no_target_tissue_overlap;max_exact_age_below_104w_or_missing | Glial Reactivity and Cognitive Decline Follow Chronic Heterochromatin Loss in Neurons | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153331 |
| GSE266961 | False | not_gse121141_like | do_not_download |  |  | bismark_cov_per_sample_tar | no_104w_plus_brain_cortex_heart_lung_samples;no_target_tissue_overlap;max_exact_age_below_104w_or_missing | Novel enzyme-based reduced representation method for DNA methylation profiling with low inputs | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266961 |
| GSE281062 | False | not_gse121141_like | do_not_download |  |  | unknown_methylation_raw_tar | processed_cov_not_detected;no_104w_plus_brain_cortex_heart_lung_samples;no_target_tissue_overlap | Methylation Array Signals are Predictive of Chronological Age Without Bisulfite Conversion | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE281062 |
| GSE292344 | False | not_gse121141_like | do_not_download |  |  | bismark_cov_per_sample_tar | no_104w_plus_brain_cortex_heart_lung_samples;no_target_tissue_overlap;max_exact_age_below_104w_or_missing | Hypercapnia induces long-term DNA methylation changes in skeletal muscle | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE292344 |
| GSE304754 | False | not_gse121141_like | do_not_download |  |  | unknown_or_non_methylation | processed_cov_not_detected;no_104w_plus_brain_cortex_heart_lung_samples;no_target_tissue_overlap | Multi-omics analysis highlights the link of aging-related cognitive decline with systemic inflammation and alterations of tissue-maintenance [RRBS] | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE304754 |
| GSE319289 | False | not_gse121141_like | do_not_download |  |  | bismark_cov_per_sample_tar | no_104w_plus_brain_cortex_heart_lung_samples;no_target_tissue_overlap;max_exact_age_below_104w_or_missing | The Age-Dependent Resident Myonuclear Multi-Omic Response to an Acute Skeletal Muscle Hypertrophic Stimulus in Mice | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE319289 |
| GSE93957 | False | not_gse121141_like | already_in_matrix_or_current_target |  |  | wgbs_cov_per_sample_tar | already_integrated;no_104w_plus_brain_cortex_heart_lung_samples;no_target_tissue_overlap;max_exact_age_below_104w_or_missing | Multi-tissue DNA methylation age predictor in mouse | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE93957 |
| GSE121141 | False | not_gse121141_like | already_in_matrix_or_current_target |  |  | bismark_cov_per_sample_tar | already_integrated;no_104w_plus_brain_cortex_heart_lung_samples;no_target_tissue_overlap;max_exact_age_below_104w_or_missing;limited_sample_specific_age_parse | A whole lifespan mouse multi-tissue DNA methylation clock | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121141 |
| GSE129712 | False | not_gse121141_like | do_not_download |  |  | unknown_or_non_methylation | processed_cov_not_detected;no_104w_plus_brain_cortex_heart_lung_samples;no_target_tissue_overlap | Transcriptional and epigenetic landscape of the mouse intestine during aging identifies key molecular drivers of aging -associated dysfunctions and diseases [RRBS] | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129712 |

## GSE225166 Adapter Audit

- Matched COV files: 232/232
- RAW archive size: 5042442240 bytes
- Context blockers: series_single_cell_or_low_coverage_age_clock_context, iPCRtag_or_iTAG_files_indicate_low_coverage_tagged_context, no_104w_plus_cortex_samples
- Headline benchmark allowed: False

| tissue | age_weeks | n_samples | n_cov_files | low_coverage_or_tag_rate | single_cell_context_rate |
| --- | --- | --- | --- | --- | --- |
| B_cells | 104.297 | 2 | 2 | 1.0 | 0.0 |
| CD4_T_Cells | 104.297 | 2 | 2 | 1.0 | 0.0 |
| CD8_T_Cells | 104.297 | 2 | 2 | 1.0 | 0.0 |
| Skin | 104.297 | 2 | 2 | 1.0 | 0.0 |
| Testes | 104.297 | 1 | 1 | 1.0 | 0.0 |
| brain_other | 104.297 | 13 | 13 | 1.0 | 0.0 |
| heart | 13.037 | 10 | 10 | 1.0 | 0.0 |
| heart | 65.186 | 12 | 12 | 1.0 | 0.0 |
| heart | 104.297 | 11 | 11 | 1.0 | 0.0 |
| heart |  | 1 | 1 | 1.0 | 0.0 |
| intestine | 104.297 | 2 | 2 | 1.0 | 0.0 |
| kidney | 104.297 | 2 | 2 | 1.0 | 0.0 |
| liver | 104.297 | 11 | 11 | 1.0 | 0.0 |
| lung | 13.037 | 12 | 12 | 1.0 | 0.0 |
| lung | 26.074 | 2 | 2 | 1.0 | 0.0 |
| lung | 39.111 | 2 | 2 | 1.0 | 0.0 |
| lung | 65.186 | 14 | 14 | 1.0 | 0.0 |
| lung | 104.297 | 13 | 13 | 1.0 | 0.0 |
| skeletal_muscle | 104.297 | 2 | 2 | 1.0 | 0.0 |

## Decision

No new direct mainline-ready old brain_cortex/heart/lung bulk RRBS dataset was found. Do not download a new dataset for v8.5 and do not start autoresearch.

GSE225166 has useful processed COV files and old heart/lung support, but its series context is single-cell/low-coverage age prediction and it lacks 104w+ cortex samples. It should remain an auxiliary adapter-smoke candidate, not a chronological headline benchmark dataset.

## Outputs

- `results/v8_5_targeted_geo_refresh/targeted_search_inventory.csv`
- `results/v8_5_targeted_geo_refresh/targeted_search_samples.csv`
- `results/v8_5_targeted_geo_refresh/targeted_search_supplements.csv`
- `results/v8_5_targeted_geo_refresh/v8_5_candidate_ranking.csv`
- `results/v8_5_targeted_geo_refresh/gse225166_adapter_audit.csv`
- `results/v8_5_targeted_geo_refresh/gse225166_tissue_age_summary.csv`
- `results/v8_5_targeted_geo_refresh/gse225166_adapter_audit.json`

## Official Sources

- https://www.ncbi.nlm.nih.gov/geo/info/geo_paccess.html
- https://www.ncbi.nlm.nih.gov/geo/info/download.html
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225166
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE232547
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171236
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138368
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151541

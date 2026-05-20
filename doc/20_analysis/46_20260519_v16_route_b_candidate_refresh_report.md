# v16 Route B Candidate Refresh Report

Date: 2026-05-19T06:16:39.980568+00:00

## Summary

v16 refreshed Route B public candidates using official GEO/E-Utils/SOFT/FTP,
SRA RunInfo, and ENA file report metadata. This run did not download FASTQ,
run Bismark, train models, or start autoresearch.

Decision: `pilot_candidates_prioritized_pending_explicit_single_candidate_authorization`.

Next action: `run_mode_pilot_for_top_ranked_candidate_only_after_explicit_authorize_pilot`.

## Tier Counts

| tier | count |
| --- | --- |
| P1_processed_headline | 1 |
| P2_auxiliary | 13 |
| P4_blocked | 56 |

## Candidate Gate Table

| dataset | tier | age_known_fraction | old_exact_target_n | old_exact_target_tissues | context_blockers | common_regions_estimate | gate_status | blocker_type |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE92486 | P4_blocked | 1.0 | 0 |  |  | 6500 | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label |
| GSE175410 | P4_blocked | 1.0 | 0 |  |  |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label |
| GSE213628 | P1_processed_headline | 0.975 | 12 | heart;lung |  |  | processed_adapter_smoke_or_matrix_gate_required | processed_files_present_common_regions_unknown_or_failed |
| GSE224442 | P4_blocked | 1.0 | 0 |  | organoid_or_in_vitro | 65827 | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label;blocked_context |
| GSE233879 | P4_blocked | 0.35 | 0 |  | single_cell_or_single_nucleus;targeted_cell_type_or_sorted |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label;blocked_context |
| GSE286302 | P2_auxiliary | 0.0 | 0 |  |  | 63524 | keep_for_diagnostic_or_adapter_evidence_not_headline | sample_specific_exact_age_coverage_lt_95pct;old_age_group_label_not_exact_age;prior_failed_candidate:sample_specific_exact_age_gate_failed;old_target_tissue_gate_failed |
| GSE295059 | P4_blocked | 0.9714 | 0 |  | organoid_or_in_vitro |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label;blocked_context |
| GSE156557 | P4_blocked | 1.0 | 0 |  |  |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label |
| GSE150670 | P4_blocked | 0.4702 | 0 |  |  |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label |
| GSE281602 | P2_auxiliary | 1.0 | 4 | heart | targeted_cell_type_or_sorted | 25321 | keep_for_diagnostic_or_adapter_evidence_not_headline | non_bulk_or_cell_context;common_region_gate_failed;prior_failed_candidate:bulk_context_gate_failed;common_regions_lt_50000 |
| GSE80672 | P2_auxiliary | 0.1647 | 6 | lung |  |  | keep_for_diagnostic_or_adapter_evidence_not_headline | sample_specific_exact_age_coverage_lt_95pct;processed_schema_unknown |
| GSE179880 | P4_blocked | 0.0 | 0 |  |  |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label |
| GSE304754 | P4_blocked | 1.0 | 0 |  |  |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label |
| GSE184267 | P4_blocked | 1.0 | 0 |  |  |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label |
| GSE232547 | P2_auxiliary | 0.0 | 0 |  |  |  | keep_for_diagnostic_or_adapter_evidence_not_headline | sample_specific_exact_age_coverage_lt_95pct;old_age_group_label_not_exact_age;prior_failed_candidate:sample_specific_exact_age_gate_failed;old_target_tissue_gate_failed;schema_smoke_failed;assembly_not_traceable;common_regions_lt_50000 |
| GSE290397 | P4_blocked | 0.0 | 0 |  |  |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label |
| GSE44117 | P4_blocked | 0.48 | 0 |  | targeted_cell_type_or_sorted |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label;blocked_context |
| GSE47815 | P4_blocked | 1.0 | 0 |  | targeted_cell_type_or_sorted |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label;blocked_context |
| GSE134267 | P4_blocked | 1.0 | 0 |  |  |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label |
| GSE134636 | P4_blocked | 1.0 | 0 |  |  |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label |
| GSE151541 | P2_auxiliary | 1.0 | 0 |  | targeted_cell_type_or_sorted |  | keep_for_diagnostic_or_adapter_evidence_not_headline | old_age_group_label_not_exact_age;non_bulk_or_cell_context |
| GSE157553 | P4_blocked | 0.0 | 0 |  |  |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label |
| GSE180433 | P4_blocked | 1.0 | 0 |  | targeted_cell_type_or_sorted |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label;blocked_context |
| GSE197045 | P4_blocked | 1.0 | 0 |  | targeted_cell_type_or_sorted |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label;blocked_context |
| GSE225305 | P4_blocked | 1.0 | 0 |  |  |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label |
| GSE281062 | P4_blocked | 1.0 | 0 |  |  |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label |
| GSE129712 | P4_blocked | 1.0 | 0 |  |  |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label |
| GSE135208 | P4_blocked | 1.0 | 0 |  |  |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label |
| GSE143890 | P2_auxiliary | 0.3478 | 0 |  |  |  | keep_for_diagnostic_or_adapter_evidence_not_headline | sample_specific_exact_age_coverage_lt_95pct;old_age_group_label_not_exact_age |
| GSE225166 | P2_auxiliary | 0.9914 | 35 | brain_cortex;heart;lung | single_cell_or_single_nucleus;low_coverage_or_itag |  | keep_for_diagnostic_or_adapter_evidence_not_headline | non_bulk_or_cell_context;prior_failed_candidate:bulk_context_gate_failed;schema_smoke_failed;assembly_not_traceable;common_regions_lt_50000 |
| GSE231658 | P4_blocked | 1.0 | 0 |  |  |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label |
| GSE266961 | P4_blocked | 0.8108 | 0 |  | targeted_cell_type_or_sorted |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label;blocked_context |
| GSE319289 | P4_blocked | 1.0 | 0 |  | targeted_cell_type_or_sorted |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label;blocked_context |
| GSE325694 | P4_blocked | 0.0 | 0 |  |  |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label |
| GSE134238 | P4_blocked | 0.0 | 0 |  | targeted_cell_type_or_sorted |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label;blocked_context |
| GSE134872 | P4_blocked | 1.0 | 0 |  |  |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label |
| GSE140125 | P4_blocked | 1.0 | 0 |  |  |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label |
| GSE160074 | P2_auxiliary | 1.0 | 0 |  |  |  | keep_for_diagnostic_or_adapter_evidence_not_headline | old_age_group_label_not_exact_age |
| GSE215310 | P4_blocked | 1.0 | 0 |  |  |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label |
| GSE218648 | P4_blocked | 1.0 | 0 |  | targeted_cell_type_or_sorted |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label;blocked_context |
| GSE233734 | P4_blocked | 1.0 | 0 |  |  |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label |
| GSE312263 | P4_blocked | 0.0 | 0 |  |  |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label |
| GSE77079 | P4_blocked | 1.0 | 0 |  |  |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label |
| GSE85772 | P4_blocked | 0.9474 | 0 |  |  |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label |
| GSE103886 | P4_blocked | 1.0 | 0 |  |  |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label |
| GSE153331 | P4_blocked | 0.4507 | 0 |  | targeted_cell_type_or_sorted |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label;blocked_context |
| GSE154949 | P4_blocked | 1.0 | 0 |  | organoid_or_in_vitro |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label;blocked_context |
| GSE158455 | P4_blocked | 1.0 | 0 |  |  |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label |
| GSE163037 | P4_blocked | 0.0 | 0 |  |  |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label |
| GSE213723 | P2_auxiliary | 0.9792 | 12 | heart;lung | superseries_or_untraceable |  | keep_for_diagnostic_or_adapter_evidence_not_headline | non_bulk_or_cell_context |
| GSE232546 | P4_blocked | 0.0 | 0 |  |  |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label |
| GSE236776 | P4_blocked | 1.0 | 0 |  | single_cell_or_single_nucleus;targeted_cell_type_or_sorted |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label;blocked_context |
| GSE236784 | P4_blocked | 1.0 | 0 |  | single_cell_or_single_nucleus;targeted_cell_type_or_sorted |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label;blocked_context |
| GSE236789 | P4_blocked | 1.0 | 0 |  | single_cell_or_single_nucleus;targeted_cell_type_or_sorted |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label;blocked_context |
| GSE292344 | P4_blocked | 1.0 | 0 |  |  |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label |
| GSE120132 | P4_blocked | 1.0 | 0 |  |  |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label |
| GSE134398 | P2_auxiliary | 0.6173 | 8 | lung | targeted_cell_type_or_sorted;superseries_or_untraceable |  | keep_for_diagnostic_or_adapter_evidence_not_headline | sample_specific_exact_age_coverage_lt_95pct;non_bulk_or_cell_context;prior_failed_candidate:sample_specific_exact_age_gate_failed;bulk_context_gate_failed;schema_smoke_failed;assembly_not_traceable;common_regions_lt_50000 |
| GSE138368 | P4_blocked | 1.0 | 0 |  |  |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label |
| GSE171236 | P4_blocked | 1.0 | 0 |  |  |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label |
| GSE249018 | P4_blocked | 1.0 | 0 |  | targeted_cell_type_or_sorted |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label;blocked_context |
| GSE292803 | P4_blocked | 1.0 | 0 |  | single_cell_or_single_nucleus;targeted_cell_type_or_sorted |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label;blocked_context |
| GSE313770 | P4_blocked | 1.0 | 0 |  | single_cell_or_single_nucleus;targeted_cell_type_or_sorted |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label;blocked_context |
| GSE83947 | P2_auxiliary | 1.0 | 180 | lung |  |  | do_not_repilot_headline | prior_failed_candidate:processed and raw Route B pilots already failed common-region gate |
| GSE221124 | P4_blocked | 1.0 | 0 |  |  |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label |
| GSE134397 | P2_auxiliary | 0.6264 | 8 | lung | targeted_cell_type_or_sorted |  | keep_for_diagnostic_or_adapter_evidence_not_headline | sample_specific_exact_age_coverage_lt_95pct;non_bulk_or_cell_context |
| GSE52266 | P4_blocked | 1.0 | 0 |  |  |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label |
| GSE232548 | P2_auxiliary | 0.0 | 0 |  | superseries_or_untraceable |  | keep_for_diagnostic_or_adapter_evidence_not_headline | sample_specific_exact_age_coverage_lt_95pct;old_age_group_label_not_exact_age;non_bulk_or_cell_context |
| GSE290999 | P4_blocked | 0.6385 | 0 |  | single_cell_or_single_nucleus |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label;blocked_context |
| GSE49191 | P4_blocked | 0.4167 | 0 |  | targeted_cell_type_or_sorted;superseries_or_untraceable |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label;blocked_context |
| GSE53742 | P4_blocked | 0.5 | 0 |  | targeted_cell_type_or_sorted |  | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label;blocked_context |

## Pilot Priority Queue

| pilot_priority_rank | dataset | tier | old_exact_target_tissues | old_exact_target_n | has_sample_processed_files | ena_total_fastq_bytes | gate_status | blocker_type |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | GSE213628 | P1_processed_headline | heart;lung | 12 | True | 5175141456 | processed_adapter_smoke_or_matrix_gate_required | processed_files_present_common_regions_unknown_or_failed |

## Rules Applied

- `GSE83947` remains `do_not_repilot_headline` because processed and raw pilots
  already failed the common-region gate.
- P2 candidates are diagnostic only and cannot enter headline pilot.
- P4 candidates are blocked and cannot enter pilot.
- Any pilot remains limited to one candidate and 2-3 samples after explicit
  `--authorize-pilot`.

## Guardrails

- Download authorized: `false`
- Bismark authorized: `false`
- Training authorized: `false`
- Autoresearch authorized: `false`
- Human clock CpG mapping: forbidden
- Dummy AUC: forbidden

## Outputs

- `results/ralph_v16_route_b_candidate_refresh/candidate_refresh_table.csv`
- `results/ralph_v16_route_b_candidate_refresh/candidate_gate_table.csv`
- `results/ralph_v16_route_b_candidate_refresh/pilot_priority_queue.csv`
- `results/ralph_v16_route_b_candidate_refresh/pilot_run_manifest.csv`
- `results/ralph_v16_route_b_candidate_refresh/network_resolution_log.jsonl`


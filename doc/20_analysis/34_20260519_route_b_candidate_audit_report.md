# Route B Candidate Audit Report

Date: 2026-05-19T08:15:53+08:00

## Summary

This audit reviews existing v10-v12 raw-lead and BioSample/RunInfo metadata for a possible Route B minimal FASTQ/Bismark pilot. It uses cached official metadata only.

No FASTQ download, Bismark, training, matrix rebuild, or autoresearch was run.

## Decision

- Headline Route B pilot candidates: `0`.
- Diagnostic-only age-group candidates: `1`.
- Blocked raw candidates: `1`.

Current conclusion: no raw FASTQ pilot is authorized. The only verified raw candidates are either diagnostic-only without sample-specific exact age or blocked by cell-type/non-bulk context.

## Candidate Findings

### GSE281602

- decision: `blocked_context_not_route_b_headline`
- reason: Runs have exact target-tissue metadata but are cell-type/single-cell/other non-bulk context.
- runs: `7`, total size: `34479.0` MB
- tissues: `heart`
- exact-age runs: `7`, old exact target raw runs: `4`, headline-eligible old exact target runs: `0`
- v12 pilot status: `not_authorized`; reason: `do_not_download_fastq`
- recommended action: `do_not_download_fastq; celltype_or_nonbulk_context_blocks_headline`

### GSE286302

- decision: `diagnostic_only_age_group_no_exact_age`
- reason: Runs are target-tissue bisulfite data but age is only an age group, not sample-specific exact age.
- runs: `7`, total size: `19707.0` MB
- tissues: `lung`
- exact-age runs: `0`, old exact target raw runs: `0`, headline-eligible old exact target runs: `0`
- v12 pilot status: `not_authorized`; reason: `diagnostic_only_do_not_promote_without_exact_age_bulk_context`
- recommended action: `do_not_download_fastq_for_headline; diagnostic_only_if_separately_approved`

## Diagnostic Pilot Rows

`route_b_recommended_pilot_samples.csv` contains `3` rows. All rows are marked `selection_status=not_authorized`.

## Raw Lead Backlog

`route_b_raw_lead_backlog.csv` contains `42` deduplicated project keys from v11 raw-lead records.

Manual backlog rows require a fresh official BioSample/RunInfo refresh before any pilot decision.

## Processed Data Context

`route_b_processed_context.csv` contains `6` datasets where processed methylation files already exist or were preferred in v11.3/v11.4.

Processed data remains preferred over raw FASTQ whenever schema and metadata gates are adequate.

## Outputs

- `results/route_b_candidate_audit/route_b_candidate_audit.csv`
- `results/route_b_candidate_audit/route_b_sample_audit.csv`
- `results/route_b_candidate_audit/route_b_recommended_pilot_samples.csv`
- `results/route_b_candidate_audit/route_b_raw_lead_backlog.csv`
- `results/route_b_candidate_audit/route_b_processed_context.csv`
- `results/route_b_candidate_audit/route_b_decision_state.json`

## Guardrails

- No raw FASTQ download is authorized by this audit.
- Route B still requires separate minimal ETL approval.
- Age-group-only data cannot become a headline chronological benchmark.
- Cell-type-specific data cannot become a headline bulk tissue benchmark.

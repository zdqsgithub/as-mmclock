# Route B Manual Backlog Official Metadata Refresh Report

Date: 2026-05-19T08:25:56

## Summary

This refresh checks the manual Route B raw-lead backlog using official NCBI RunInfo/SRA XML and ENA read_run metadata only.

No FASTQ download, Bismark, training, matrix rebuild, or autoresearch was run.

## Decision

- Metadata groups refreshed: `32`.
- Run rows audited: `46`.
- Headline Route B pilot-ready groups: `0`.
- Diagnostic-only or insufficient-replicate groups: `0`.
- Blocked groups: `32`.
- Overall decision: `no_route_b_headline_pilot_candidate_ready`.

No raw FASTQ pilot is authorized by this report.

## Key Findings

### PRJNA1231360

- decision: `blocked_non_target_tissue`
- reason: Official metadata does not include brain_cortex/heart/lung target tissue.
- runs audited: `7`, tissues: `intestine`
- age range weeks: `28.0` to `28.0`
- assay-pass / target / exact-age / old-age: `7` / `0` / `7` / `0`
- recommended action: `do_not_download_fastq`

### PRJNA1254502

- decision: `blocked_non_bisulfite_methylation_assay`
- reason: Official metadata does not indicate bisulfite/methylation sequencing.
- runs audited: `3`, tissues: `N/A`
- age range weeks: `N/A` to `N/A`
- assay-pass / target / exact-age / old-age: `0` / `0` / `0` / `0`
- recommended action: `do_not_download_fastq`

### PRJNA492485

- decision: `blocked_non_target_tissue`
- reason: Official metadata does not include brain_cortex/heart/lung target tissue.
- runs audited: `7`, tissues: `N/A`
- age range weeks: `N/A` to `N/A`
- assay-pass / target / exact-age / old-age: `7` / `0` / `0` / `7`
- recommended action: `do_not_download_fastq`

### ERR11479077

- decision: `blocked_non_bisulfite_methylation_assay`
- reason: Official metadata does not indicate bisulfite/methylation sequencing.
- runs audited: `1`, tissues: `N/A`
- age range weeks: `N/A` to `N/A`
- assay-pass / target / exact-age / old-age: `0` / `0` / `0` / `0`
- recommended action: `do_not_download_fastq`

### ERR11479078

- decision: `blocked_non_bisulfite_methylation_assay`
- reason: Official metadata does not indicate bisulfite/methylation sequencing.
- runs audited: `1`, tissues: `N/A`
- age range weeks: `N/A` to `N/A`
- assay-pass / target / exact-age / old-age: `0` / `0` / `0` / `0`
- recommended action: `do_not_download_fastq`

### ERR11479079

- decision: `blocked_non_bisulfite_methylation_assay`
- reason: Official metadata does not indicate bisulfite/methylation sequencing.
- runs audited: `1`, tissues: `N/A`
- age range weeks: `N/A` to `N/A`
- assay-pass / target / exact-age / old-age: `0` / `0` / `0` / `0`
- recommended action: `do_not_download_fastq`

### ERR11479080

- decision: `blocked_non_bisulfite_methylation_assay`
- reason: Official metadata does not indicate bisulfite/methylation sequencing.
- runs audited: `1`, tissues: `N/A`
- age range weeks: `N/A` to `N/A`
- assay-pass / target / exact-age / old-age: `0` / `0` / `0` / `0`
- recommended action: `do_not_download_fastq`

### ERR11479081

- decision: `blocked_non_bisulfite_methylation_assay`
- reason: Official metadata does not indicate bisulfite/methylation sequencing.
- runs audited: `1`, tissues: `N/A`
- age range weeks: `N/A` to `N/A`
- assay-pass / target / exact-age / old-age: `0` / `0` / `0` / `0`
- recommended action: `do_not_download_fastq`

### ERR11479082

- decision: `blocked_non_bisulfite_methylation_assay`
- reason: Official metadata does not indicate bisulfite/methylation sequencing.
- runs audited: `1`, tissues: `N/A`
- age range weeks: `N/A` to `N/A`
- assay-pass / target / exact-age / old-age: `0` / `0` / `0` / `0`
- recommended action: `do_not_download_fastq`

### ERR11479083

- decision: `blocked_non_bisulfite_methylation_assay`
- reason: Official metadata does not indicate bisulfite/methylation sequencing.
- runs audited: `1`, tissues: `N/A`
- age range weeks: `N/A` to `N/A`
- assay-pass / target / exact-age / old-age: `0` / `0` / `0` / `0`
- recommended action: `do_not_download_fastq`

### ERR11479084

- decision: `blocked_non_bisulfite_methylation_assay`
- reason: Official metadata does not indicate bisulfite/methylation sequencing.
- runs audited: `1`, tissues: `N/A`
- age range weeks: `N/A` to `N/A`
- assay-pass / target / exact-age / old-age: `0` / `0` / `0` / `0`
- recommended action: `do_not_download_fastq`

### ERR11479085

- decision: `blocked_non_bisulfite_methylation_assay`
- reason: Official metadata does not indicate bisulfite/methylation sequencing.
- runs audited: `1`, tissues: `N/A`
- age range weeks: `N/A` to `N/A`
- assay-pass / target / exact-age / old-age: `0` / `0` / `0` / `0`
- recommended action: `do_not_download_fastq`

## Outputs

- `results/route_b_candidate_audit/route_b_manual_backlog_runinfo_refresh.csv`
- `results/route_b_candidate_audit/route_b_manual_backlog_candidate_summary.csv`
- `results/route_b_candidate_audit/route_b_manual_backlog_network_log.jsonl`
- `results/route_b_candidate_audit/route_b_manual_backlog_decision_state.json`

## Guardrails

- This is metadata refresh only.
- Route B pilot still requires a separate explicit minimal ETL approval.
- Age-group-only candidates cannot become headline chronological benchmarks.
- Cell-type-specific, RNA-Seq, or non-target tissue candidates remain blocked.

# v14 Public-Data Rescue Report

Date: 2026-05-19T02:37:59.873609+00:00

## Summary

v14 ran a public-data RALPH Research/Audit pass against official GEO/E-Utils/SOFT/FTP metadata and SRA/ENA run-file metadata. This run did not download FASTQ, run Bismark, train models, or start autoresearch.

- Decision: `v14_public_data_raw_pilot_candidate_found_no_download_started`
- Next action: Review pilot_run_manifest and authorize only 2-3 sample raw/processed pilot if the candidate is scientifically acceptable.
- Headline candidates allowed now: 0
- P3 raw-pilot candidates: 1

## Tier Counts

| tier | count |
| --- | --- |
| P2_auxiliary | 9 |
| P4_blocked | 8 |
| P3_raw_pilot | 1 |

## Candidate Gate Table

| dataset | tier | age_known_fraction | old_exact_target_n | old_label_target_n | old_exact_target_tissues | context_blockers | headline_allowed | gate_status | blocker_type |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE286302 | P2_auxiliary | 0.0 | 0 | 8 |  |  | False | keep_for_diagnostic_or_adapter_evidence_not_headline | sample_specific_exact_age_coverage_lt_95pct;old_age_group_label_not_exact_age |
| GSE281602 | P2_auxiliary | 1.0 | 4 | 4 | heart | targeted_cell_type_or_sorted | False | keep_for_diagnostic_or_adapter_evidence_not_headline | non_bulk_or_cell_context;common_region_gate_failed |
| GSE304754 | P2_auxiliary | 1.0 | 0 | 0 |  |  | False | keep_for_diagnostic_or_adapter_evidence_not_headline | single_tissue_or_celltype;processed_cov_not_detected;no_parseable_processed_methylation_schema_detected |
| GSE225166 | P2_auxiliary | 0.9914 | 35 | 0 | brain_cortex;heart;lung | single_cell_or_single_nucleus;low_coverage_or_itag | False | keep_for_diagnostic_or_adapter_evidence_not_headline | non_bulk_or_cell_context |
| GSE134398 | P2_auxiliary | 0.6173 | 8 | 52 | lung | targeted_cell_type_or_sorted;superseries_or_untraceable | False | keep_for_diagnostic_or_adapter_evidence_not_headline | sample_specific_exact_age_coverage_lt_95pct;non_bulk_or_cell_context |
| GSE292803 | P4_blocked | 1.0 | 0 | 0 |  | single_cell_or_single_nucleus;targeted_cell_type_or_sorted | False | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label;blocked_context |
| GSE313770 | P4_blocked | 1.0 | 0 | 0 |  | single_cell_or_single_nucleus;targeted_cell_type_or_sorted | False | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label;blocked_context |
| GSE83947 | P3_raw_pilot | 1.0 | 180 | 0 | lung |  | False | prepare_2_3_sample_processed_adapter_smoke_no_download_in_this_run | processed_supplement_schema_unknown_requires_adapter_smoke |
| GSE213628 | P2_auxiliary | 0.975 | 12 | 10 | heart;lung |  | False | keep_as_existing_evidence_not_new_v14_headline | previously_integrated_or_evaluated_dataset |
| GSE80672 | P2_auxiliary | 0.1647 | 6 | 0 | lung |  | False | keep_as_existing_evidence_not_new_v14_headline | previously_integrated_or_evaluated_dataset |
| GSE232547 | P2_auxiliary | 0.0 | 0 | 12 |  |  | False | keep_for_diagnostic_or_adapter_evidence_not_headline | sample_specific_exact_age_coverage_lt_95pct;old_age_group_label_not_exact_age |
| GSE290397 | P4_blocked | 0.0 | 0 | 0 |  |  | False | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label |
| GSE184267 | P4_blocked | 1.0 | 0 | 0 |  |  | False | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label |
| GSE325694 | P4_blocked | 0.0 | 0 | 0 |  |  | False | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label |
| GSE160074 | P2_auxiliary | 1.0 | 0 | 12 |  |  | False | keep_for_diagnostic_or_adapter_evidence_not_headline | old_age_group_label_not_exact_age |
| GSE312263 | P4_blocked | 0.0 | 0 | 0 |  |  | False | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label |
| GSE153331 | P4_blocked | 0.4507 | 0 | 0 |  | targeted_cell_type_or_sorted | False | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label;blocked_context |
| GSE232546 | P4_blocked | 0.0 | 0 | 0 |  |  | False | do_not_promote_to_headline | no_exact_old_target_tissue_samples;no_old_target_tissue_label |

## Interpretation

- `GSE286302` remains useful as old-lung processed-COV compatibility evidence, but it lacks sample-specific exact age and is lung-only.
- `GSE281602` has exact old heart samples, but the context is cardiomyocyte/cell-type specific and its common-region overlap failed the project gate.
- `GSE304754` is brain-adjacent hippocampus rather than brain_cortex/cortex headline support.
- Single-cell/low-coverage/iTAG brain or methylation-age datasets are reference or adapter-audit material only.

## Guardrails

- Training authorized: `false`
- Download authorized: `false` in this run; only run metadata was queried.
- Bismark authorized: `false`
- Autoresearch authorized: `false`
- Human clock CpG mapping: forbidden
- Dummy AUC: forbidden

## Pilot Candidate

The only v14 P3 candidate is `GSE83947`.

- Gate: `prepare_2_3_sample_processed_adapter_smoke_no_download_in_this_run`
- Rationale: exact-age old lung samples exist, but the public supplement is
  `CX_report.txt.gz` / TXT-style processed methylation, so it needs a dedicated
  adapter smoke before any matrix or benchmark.
- Pilot manifest now selects sample-level processed supplement URLs, not raw
  FASTQ, because GEO exposes individual small TXT files in addition to the large
  series TAR.
- Download remains unauthorized in this run.

The narrow `GSE83947` processed TXT/CX_report adapter smoke was then executed
for three old lung files. This was a small GEO sample-level processed
supplement smoke, not FASTQ, not Bismark, and not training.

Smoke result:

| Metric | Value |
| --- | ---: |
| Files downloaded | 3 |
| Total bytes | 17,644,422 |
| Schema detected | `bismark_cx_report` |
| Beta range valid | 3/3 |
| Parseable files | 3/3 |
| Max smoke common 5kb regions with v8.2 reference | 121 |

Key outputs:

- `results/ralph_v14_public_data_rescue/gse83947_processed_adapter_smoke/gse83947_processed_adapter_smoke_state.json`
- `results/ralph_v14_public_data_rescue/gse83947_processed_adapter_smoke/gse83947_processed_adapter_smoke_manifest.csv`
- `results/ralph_v14_public_data_rescue/gse83947_processed_adapter_smoke/download_log.jsonl`

Interpretation: the sample-level TXT files are parseable and beta values are
valid, but the bounded smoke only scanned 250k rows per file and does not yet
prove the full 5kb common-region matrix gate. The next step is a dedicated
`GSE83947` CX_report-to-region matrix adapter/gate before any Learn/Benchmark
approval.

## GSE83947 Matrix Gate

The dedicated CX_report-to-5kb adapter was implemented and run on the same three
processed old lung files. It parses CX_report columns as:
`chrom, pos, strand, methylated, unmethylated, context, trinucleotide`, applies
`coverage >= 5`, excludes sex/MT chromosomes, aggregates mean beta per 5kb
region, and checks overlap with the v8.2 reference matrix.

| Metric | Value |
| --- | ---: |
| Samples | 3 |
| 5kb regions | 1,500 |
| v8.2 reference regions | 65,870 |
| Common regions with v8.2 | 539 |
| Common-region gate `>=50,000` | false |
| Beta range | 0.0-0.252459 |

Decision: `matrix_gate_failed_auxiliary_only`.

This closes the v14 public-data rescue loop for `GSE83947`: it is technically
parseable, but far too sparse relative to the project 5kb reference matrix to
enter RALPH Learn/Benchmark. It remains auxiliary adapter evidence only.

Additional outputs:

- `scripts/validate/run_v14_gse83947_region_matrix_gate.py`
- `results/ralph_v14_public_data_rescue/gse83947_region_matrix_gate/GSE83947_matrix_gate_state.json`
- `results/ralph_v14_public_data_rescue/gse83947_region_matrix_gate/GSE83947_region_matrix_5kb.parquet`
- `results/ralph_v14_public_data_rescue/gse83947_region_matrix_gate/GSE83947_sample_parse_stats.csv`

## Route B Raw Pilot Preflight

Because the processed CX_report matrix gate failed, a Route B raw pilot preflight
was prepared for the highest-depth old lung bisulfite runs in `GSE83947`. The
preflight uses official GEO SOFT to map GSM -> BioSample/SRX and ENA read_run
metadata to map SRX/BioSample -> SRR/FASTQ. It did not download FASTQ.

Selected raw pilot candidates:

| sample_id | run | age_weeks | tissue | FASTQ bytes |
| --- | --- | ---: | --- | ---: |
| GSM2223569 | SRR3738487 | 108.643 | lung | 408,118,908 |
| GSM2223613 | SRR3738531 | 108.643 | lung | 366,111,345 |
| GSM2223583 | SRR3738501 | 108.643 | lung | 293,772,100 |

Metadata gate: passed. Environment gate: failed.

Blockers:

- `bismark` missing
- `bowtie2` missing
- `samtools` missing
- mm10/GRCm38 Bismark bisulfite genome index missing

Decision: `do_not_download_fastq_until_environment_gate_passes`.

Additional outputs:

- `scripts/validate/run_v14_gse83947_raw_pilot_preflight.py`
- `results/ralph_v14_public_data_rescue/gse83947_raw_pilot_preflight/route_b_gse83947_raw_pilot_manifest.csv`
- `results/ralph_v14_public_data_rescue/gse83947_raw_pilot_preflight/route_b_gse83947_raw_pilot_preflight_state.json`
- `doc/20_analysis/44_20260519_v14_gse83947_route_b_raw_pilot_preflight.md`

## Route B Environment Provisioning

After explicit approval, the Route B environment gate was provisioned locally
inside the project. This step installed/registered project-local tools and built
a traceable GRCm38/mm10 Bismark index; it did not train models or start
autoresearch.

| Item | Value |
| --- | --- |
| Environment status | `ready` |
| Bismark | `tools/route_b_bismark_env/bin/bismark` |
| Bowtie2 | `tools/route_b_bismark_env/bin/bowtie2` |
| Samtools | `tools/route_b_bismark_env/bin/samtools` |
| Reference | Ensembl release 102 `Mus_musculus.GRCm38.dna.primary_assembly.fa.gz` |
| Assembly | `GRCm38/mm10` |
| Bismark index | `references/GRCm38_ensembl102/Bisulfite_Genome` |
| FASTA SHA256 | `a3628d245170e72f2f7c40354c4c496052f46808264d73c94720bfd546322184` |
| FASTA.GZ SHA256 | `285bc481d583ab65b13d91853bf743acf950710afb3302264a4b4f116b6049c1` |

Output:

- `scripts/validate/provision_route_b_bismark_environment.py`
- `results/ralph_v14_public_data_rescue/route_b_environment/route_b_bismark_environment_manifest.json`

## Route B FASTQ Download

After the environment gate passed and explicit approval was given, the selected
three `GSE83947` old-lung runs were downloaded as the minimal raw pilot.

| Metric | Value |
| --- | ---: |
| Samples | 3 |
| FASTQ files | 6 |
| Total bytes | 1,068,002,353 |
| Size checks | passed |
| MD5 checks | passed |
| Training authorized | false |
| Autoresearch authorized | false |

Output:

- `scripts/validate/run_v14_gse83947_raw_fastq_download.py`
- `results/ralph_v14_public_data_rescue/gse83947_raw_fastq_pilot/raw_fastq_download_manifest.csv`
- `results/ralph_v14_public_data_rescue/gse83947_raw_fastq_pilot/raw_fastq_download_state.json`

## Route B Minimal Bismark Pilot

The three downloaded FASTQ pairs were aligned with local Bismark/Bowtie2 against
the project GRCm38/mm10 Bismark index. Methylation extraction was run as a
CpG-oriented coverage path, and 5kb region aggregation used the same
sex/MT-exclusion and beta-range checks as the existing matrix gates.

| Metric | Value |
| --- | ---: |
| Samples aligned | 3 |
| Pilot 5kb regions | 5,796 |
| v8.2 reference regions | 65,870 |
| Common regions with v8.2 | 1,814 |
| Common-region gate `>=50,000` | false |
| Beta range | 0.0-1.0 |
| Decision | `matrix_gate_failed_auxiliary_only` |

Interpretation: the raw FASTQ/Bismark pilot is technically executable and
traceable, but `GSE83947` remains far too sparse relative to the project
multi-dataset 5kb region reference. It cannot enter RALPH Learn/Benchmark and
does not change the v13 Route C boundary. It remains auxiliary old-lung
diagnostic evidence only.

Outputs:

- `scripts/validate/run_v14_gse83947_bismark_pilot.py`
- `results/ralph_v14_public_data_rescue/gse83947_bismark_pilot/gse83947_bismark_pilot_state.json`
- `results/ralph_v14_public_data_rescue/gse83947_bismark_pilot/GSE83947_bismark_pilot_region_matrix_5kb.parquet`
- `results/ralph_v14_public_data_rescue/gse83947_bismark_pilot/GSE83947_bismark_pilot_sample_parse_stats.csv`
- `results/ralph_v14_public_data_rescue/gse83947_bismark_pilot/common_regions_with_v8_2.csv`

## Contract Tests

- `scripts/validate/test_v14_public_data_rescue_contracts.py`: 9/9 passed.
- `scripts/validate/test_v13_delivery_contracts.py`: 6/6 passed.

## Outputs

- `results/ralph_v14_public_data_rescue/candidate_refresh_table.csv`
- `results/ralph_v14_public_data_rescue/candidate_gate_table.csv`
- `results/ralph_v14_public_data_rescue/pilot_run_manifest.csv`
- `results/ralph_v14_public_data_rescue/network_resolution_log.jsonl`
- `results/ralph_v14_public_data_rescue/ralph_decision_state.json`

## Sources

- GEO Download / Programmatic Access / SOFT official documentation.
- NCBI SRA download documentation and ENA file report API documentation.
- GEO accession pages and cached SOFT metadata for candidate accessions.

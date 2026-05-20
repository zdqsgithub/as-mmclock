# v11 Data Strategy and Raw Lead RFC Report

Date: 2026-05-18

## Summary

v11 broadened discovery beyond the v10 GEO processed/seed space using official
PubMed, SRA E-Utils, GEO GSM-to-GSE resolution, and GEO/SOFT/FTP verification.
It did not download large files, build matrices, run FASTQ/Bismark, train
models, or start autoresearch.

Final state: `rfc_required`

Decision: `v11_raw_literature_leads_need_manual_biosample_verification`

## Outputs

- `scripts/validate/run_v11_data_strategy_rfc.py`
- `scripts/validate/summarize_v11_raw_pilot_leads.py`
- `results/ralph_v11_data_strategy/literature_sra_source_records.csv`
- `results/ralph_v11_data_strategy/literature_sra_accession_leads.csv`
- `results/ralph_v11_data_strategy/gsm_to_gse_resolution.csv`
- `results/ralph_v11_data_strategy/verified_new_gse_gate_table.csv`
- `results/ralph_v11_data_strategy/raw_accession_rfc_leads.csv`
- `results/ralph_v11_data_strategy/raw_project_rfc_summary.csv`
- `results/ralph_v11_data_strategy/raw_biosample_verification_manifest.csv`
- `results/ralph_v11_data_strategy/raw_biosample_verification_summary.csv`
- `results/ralph_v11_data_strategy/raw_biosample_decision_state.json`
- `results/ralph_v11_data_strategy/v11_data_strategy_rfc_report.md`

## Gate Results

| Metric | Value |
| --- | ---: |
| PubMed/SRA metadata records | 363 |
| High-scoring raw leads | 198 |
| New GSE verified through GEO/SOFT/FTP | 4 |
| New P1 candidates | 0 |
| New P3 candidates | 0 |

The newly verified GSE accessions were all P4:

| Dataset | Decision |
| --- | --- |
| GSE232548 | superseries traceability risk; not headline |
| GSE290999 | single-cell/targeted context; not headline |
| GSE49191 | superseries/HSC context; not headline |
| GSE53742 | HSC/cell-type context; not headline |

## Raw Project RFC Candidates

The raw-project summary identified two first-pass candidates that justify
BioSample/RunInfo verification, but not FASTQ download yet:

| Candidate | Evidence | Required next check |
| --- | --- | --- |
| GSE286302 / PRJNA1208582 / SRP556367 | SRA titles include old lung bisulfite runs with aged vehicle and aged 3dA replicates | verify BioSample age/tissue, control status, release/accessibility, run layout, size |
| GSE281602 / PRJNA1184779 / SRP544506 | GSM resolves to GSE281602; SRA titles include heart bisulfite/RRBS-like records with age labels | verify whether aged heart samples are sample-specific and not only summary-text artifacts |

PRJNA882035/SRP398068 resolves to GSE213628/GSE213723 and is lower priority
because GSE213628 has already been integrated/diagnosed in v8.

## Interpretation

v11 did not find a public processed dataset that can be promoted directly into
the matrix. It did find raw-level leads that are plausible enough for a minimal
ETL RFC. This changes the next step from "search more GEO processed files" to
"verify raw sample metadata at BioSample/RunInfo level."

The project still must not run autoresearch, embedding search, deep learning, or
full FASTQ ETL. The only justified next action is a v11.1 BioSample/RunInfo
verification manifest for the top raw candidates.

## Recommended v11.1

Build `results/ralph_v11_data_strategy/raw_biosample_verification_manifest.csv`
for:

1. `GSE286302 / PRJNA1208582 / SRP556367`
2. `GSE281602 / PRJNA1184779 / SRP544506`

Only if that manifest passes should the project request a 2-3 run minimal
FASTQ/Bismark pilot.

## v11.1 Verification Result

The v11.1 RunInfo/GEO/BioSample verification manifest was executed. It found no
headline chronological pilot candidate.

| Candidate | Verified result | Decision |
| --- | --- | --- |
| GSE286302 / PRJNA1208582 / SRP556367 | 7 old lung RRBS/bisulfite runs; BioSample/GEO confirm lung and age group `aged`, but no exact age_weeks | optional diagnostic old-lung pilot only |
| GSE281602 / PRJNA1184779 / SRP544506 | 7 heart/cardiomyocyte RRBS/bisulfite runs with exact 4mo/28mo labels | blocked from headline because samples are cardiomyocyte/cell-type specific, not bulk heart |

Final v11.1 state:

- Status: `diagnostic_only`
- Decision: `v11_1_only_diagnostic_old_lung_pilot_candidate`
- Chronological pilot rows: 0
- Diagnostic pilot rows: 7

Therefore the raw path still does not provide a headline chronological-age fix.
The only technically allowed pilot is an optional 2-3 run GSE286302 old-lung
diagnostic pilot. It cannot be used as the headline GSE121141 old104+ solution
because it lacks exact sample age and only covers old lung.

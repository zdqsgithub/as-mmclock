# v12 Data Strategy and Benchmark Redefinition SOP

Date: 2026-05-19

## Purpose

v12 starts after v11.x closes the downloaded processed-data backlog without a
headline old target-tissue candidate. v12 is a RALPH data-strategy loop, not a
model-optimization loop. Its purpose is to decide whether public traceable data
can still repair the GSE121141 old104+ brain_cortex/heart/lung support gap, or
whether the benchmark must be formally redefined.

## State Machine

1. Research
   - Round 1: reuse official GEO/E-Utils/SOFT/FTP/SRA/ENA refresh evidence from v9-v10.
   - Round 2: reuse PubMed/SRA accession evidence and verified GSE leads from v11.
   - Round 3: reuse raw BioSample/RunInfo verification and v11.x downloaded backlog conversion gates.
   - Non-official sources may only add clues; every candidate must return to official GEO/SRA/ENA metadata before promotion.

2. Audit
   - Allowed: HEAD/filelist/SOFT/RunInfo/BioSample checks and 2-3 small schema-smoke files.
   - Forbidden by default: large tar downloads, full FASTQ downloads, full Bismark ETL, model training, and autoresearch.
   - Every candidate must be classified as P1, P2, P3, P4, or P5.

3. Learn
   - Only P1 processed headline data or approved P3 minimal-ETL data may build a matrix.
   - Matrix gate requires metadata overlap, beta range `[0,1]`, sex/MT exclusion, assembly traceability, and common 5kb regions `>=50000`.
   - Fixed benchmark only: GroupKFold, LODO, GSE121141 held-out, GSE80672 held-out, random-label sanity, shuffled-CR sanity.

4. Promote or Stop
   - Promote to constrained autoresearch only after old104+ MAE improves at least 5 weeks and sanity checks pass.
   - Scientific success requires old104+ MAE improvement at least 10 weeks.
   - If three refresh rounds produce no P1/P3 headline data, generate benchmark redefinition RFC.

## Success Gates

Scientific data success requires all of:

- P1 processed headline or P3 minimal-ETL-authorized dataset.
- Mouse bulk RRBS/WGBS/bisulfite methylation.
- Sample-specific age coverage `>=95%`.
- Old target tissue: brain_cortex/cortex, heart, or lung at `>=104w`.
- No single-cell, cell-type targeted, organoid/in-vitro, low-coverage/iTAG, or untraceable superseries context.
- Parseable methylation schema and traceable assembly/liftover.
- Common 5kb regions `>=50000`.
- GSE121141 old104+ MAE `<=65.386w`.
- GSE121141 all-age held-out MAE `<=40.033w`.
- GroupKFold MAE `<=25.767w`.
- Random-label and shuffled-CR sanity checks pass.

Benchmark redefinition success requires:

- three refresh rounds are complete;
- all candidates have gate decisions and blockers;
- no P1/P3 headline data remains;
- a benchmark redefinition RFC is generated and auditable.

## Failure Gates

v12 fails if:

- no P1/P3 data exists and no auditable benchmark redefinition can be formed;
- three valid headline matrices still leave GSE121141 old104+ MAE near `70w+`;
- all candidates are auxiliary or blocked by context;
- common regions are `<50000` or sanity checks fail;
- calibration only improves the target set and does not transfer in non-target LODO;
- official network/API/schema provenance cannot be closed.

## Outputs

- `results/ralph_v12_loop/ralph_iteration_log.jsonl`
- `results/ralph_v12_loop/ralph_decision_state.json`
- `results/ralph_v12_loop/candidate_gate_table.csv`
- `results/ralph_v12_loop/network_resolution_log.jsonl`
- `results/ralph_v12_loop/candidate_smoke_manifest.csv`
- `results/ralph_v12_loop/minimal_etl_rfc_manifest.csv`
- `results/ralph_v12_loop/benchmark_redefinition_rfc.md`
- `results/ralph_v12_loop/v12_success_or_failure_report.md`

## Non-Negotiable Rules

- RALPH before autoresearch.
- No dummy CR-AUC.
- No human clock CpG mapping.
- No full-data feature selection or preprocessing.
- No biological-age claim without real held-out predictions.
- No raw FASTQ download without a separate minimal ETL approval.


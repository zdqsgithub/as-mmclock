# v12 RALPH Success and Failure Plan

Date: 2026-05-19

## Summary

v12 follows the v11.x conversion gate conclusion: the downloaded backlog did not
produce a headline old bulk brain_cortex/heart/lung mouse methylome dataset. The
next step is not training. v12 is a three-round data-strategy and benchmark
redefinition loop.

The loop has two possible successful exits:

- scientific data success, if a real P1/P3 headline dataset fixes the old
  GSE121141 target-tissue support gap;
- benchmark redefinition success, if three documented refresh rounds show that
  public traceable data is insufficient and the project defines a new auditable
  headline benchmark.

## Inputs

v12 reuses existing official evidence:

- v9/v10 GEO/E-Utils/SOFT/FTP/SRA/ENA refresh outputs;
- v11 PubMed/SRA accession and BioSample/RunInfo verification outputs;
- v11.x downloaded processed-data conversion gate outputs.

Default execution does not download FASTQ, build matrices, train models, or run
autoresearch.

## Success Definition

Scientific data success requires a P1 processed headline or P3 minimal-ETL
authorized dataset that passes sample-specific age, old target tissue, bulk
context, schema, assembly/liftover, and common-region gates. It must then pass
fixed benchmark gates:

- GSE121141 old104+ MAE `<=65.386w`;
- GSE121141 all-age held-out MAE `<=40.033w`;
- GroupKFold MAE `<=25.767w`;
- random-label sanity passes;
- GSE80672 shuffled-CR sanity passes.

Benchmark redefinition success is allowed when three refresh rounds produce no
P1/P3 headline data and the RFC explains why GSE121141 old104+
brain_cortex/heart/lung should become a stress-test rather than the headline
pass/fail benchmark.

## Failure Definition

v12 fails if:

- no headline data exists and no auditable benchmark redefinition can be formed;
- three valid headline matrices still leave old104+ MAE near `70w+`;
- all candidates are auxiliary/blocked and the project refuses benchmark
  redefinition;
- matrix/sanity gates fail;
- provenance cannot be closed through official GEO/SRA/ENA metadata.

## Implementation

Controller:

- `scripts/validate/run_v12_ralph_loop.py`

Primary outputs:

- `results/ralph_v12_loop/candidate_gate_table.csv`
- `results/ralph_v12_loop/ralph_decision_state.json`
- `results/ralph_v12_loop/benchmark_redefinition_rfc.md`
- `results/ralph_v12_loop/v12_success_or_failure_report.md`

## Expected Decision

Based on v9-v11 evidence, the expected first v12 run is benchmark redefinition
success rather than scientific data success. That is a valid endpoint: it stops
unsupported model tuning and moves the project to v13 data generation,
collaboration, or explicitly approved minimal ETL.

## Execution Result

The v12 cached-only controller was run on 2026-05-19 and produced the expected
benchmark redefinition endpoint.

Outputs:

- `results/ralph_v12_loop/candidate_gate_table.csv`
- `results/ralph_v12_loop/ralph_decision_state.json`
- `results/ralph_v12_loop/benchmark_redefinition_rfc.md`
- `results/ralph_v12_loop/v12_success_or_failure_report.md`

Decision:

`v12_benchmark_redefinition_success_no_headline_data_after_three_rounds`

Key metrics:

| Metric | Value |
|---|---:|
| Candidate rows | 107 |
| P1 processed headline candidates | 0 |
| P3 minimal ETL RFC candidates | 0 |
| Headline-allowed candidates | 0 |
| P2 auxiliary candidates | 14 |
| P4 blocked candidates | 49 |
| P5 redefinition evidence rows | 43 |

Conclusion:

Do not run autoresearch or training against the unsupported old104+
GSE121141 headline. Adopt the benchmark redefinition RFC or move to v13 data
generation/collaboration/minimal ETL.

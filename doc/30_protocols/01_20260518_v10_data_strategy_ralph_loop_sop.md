# v10 Data Strategy RALPH Loop SOP

Date: 2026-05-18

## Purpose

v10 is a data-strategy loop for mouse methylation clocks. It starts only after
v8/v9 show that model tuning, leakage-safe calibration, and the current
processed GEO candidate space do not repair GSE121141 old-age
brain_cortex/heart/lung generalization.

The loop is designed to answer one question:

Can a traceable public or minimal-ETL mouse methylome dataset provide old target
tissue support that passes matrix gates and improves GSE121141 104w+ held-out
MAE?

## Default Rule

Do not run autoresearch, deep learning, embedding search, or large downloads
until the data gate is passed. v10 is allowed to refresh official metadata,
build candidate gate tables, and plan adapter smoke tests. It is not allowed to
claim biological age performance from metadata, proxy labels, or dummy metrics.

## State Machine

1. Research
   - Query official GEO/E-Utils/SOFT and GEO FTP metadata.
   - Use non-official sources only as accession seeds.
   - Do not download large archives in this step.

2. Audit
   - Classify candidates into P1-P5 tiers.
   - Download at most 2-3 small per-sample COV/bedGraph/TXT files only when
     `--run_smoke` is explicitly used.
   - Verify sample-specific age, target tissue, assembly, coordinate schema,
     beta range, and estimated common 5kb region overlap.

3. Learn
   - Run only if a P1 or P3 candidate passes the matrix gate.
   - Build dataset matrix and multidataset matrix with traceable assembly or
     liftover provenance.
   - Run fixed benchmark suite: GroupKFold, LODO, GSE121141 held-out,
     GSE80672 CR held-out, random-label sanity, and shuffled CR sanity.

4. Promote or Stop
   - Promote only if old104+ MAE improves and all sanity gates pass.
   - Stop with failure if two broadened refresh rounds find no P1/P3 candidate,
     or if headline matrices repeatedly fail to reduce old104+ MAE.

## Candidate Tiers

| Tier | Meaning | Headline allowed |
| --- | --- | --- |
| P1 | Bulk RRBS/WGBS processed methylation, old brain_cortex/heart/lung, sample-specific age | Yes |
| P2 | Old adjacent tissue processed methylation, e.g. liver, blood, intestine, muscle | No; diagnostic only |
| P3 | Old target tissue with raw FASTQ/SRA only, suitable for minimal 2-3 sample ETL pilot | Pilot only |
| P4 | Low-coverage, tagged, single-cell, targeted cell type, organoid, or superseries context | No |
| P5 | Missing sample-specific age, missing target/adjacent old tissue, no parseable methylation evidence | No |

## Success Gate

v10 success requires all conditions below:

- At least one P1 or P3 dataset passes sample-specific age, target tissue,
  parseable methylation schema, assembly/liftover, and common 5kb regions
  `>=50000`.
- GSE121141 104w+ MAE improves by `>=10w` versus v7.5 core4 baseline
  `75.386w`.
- GSE121141 all-age held-out MAE is `<=40.033w`.
- GroupKFold MAE is `<=25.767w`.
- Random-label sanity passes: `abs(r)<0.2` and MAE returns to random-like
  level.
- GSE80672 shuffled CR sanity passes: AUC near chance and Cohen's d decreases.

## Failure Gate

v10 failure is declared when either condition holds:

- Two broadened official refresh rounds find no P1 or P3 candidate.
- Only P2/P4/P5 datasets are available, so the current public GEO
  processed/seed space cannot supply headline old target tissue support.

If failure is declared, do not continue tuning the same matrix. The next phase
must be broader acquisition outside the current seed space, a deliberate
minimal ETL design, or a redefined benchmark.

## Required Outputs

- `results/ralph_v10_loop/ralph_iteration_log.jsonl`
- `results/ralph_v10_loop/ralph_decision_state.json`
- `results/ralph_v10_loop/candidate_gate_table.csv`
- `results/ralph_v10_loop/network_resolution_log.jsonl`
- `results/ralph_v10_loop/candidate_smoke_manifest.csv`
- `results/ralph_v10_loop/minimal_etl_pilot_manifest.csv`
- `results/ralph_v10_loop/v10_success_or_failure_report.md`

## Official Source Rules

- E-Utils and SOFT are metadata/search sources.
- Supplement files must be resolved through official GEO FTP paths.
- Sample-level supplements use `/geo/samples/GSM.../<GSM>/suppl/`.
- Series-level archives use `/geo/series/GSE.../<GSE>/suppl/`.
- Raw pilot decisions may use SRA/ENA metadata, but FASTQ/Bismark ETL is a
  separate plan and must not be mixed into the metadata-only loop.


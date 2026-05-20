# v11.x RALPH Conversion Gate Report

Date: 2026-05-19 01:16 CST

Project standard: AS-DS-Ops v3.3. This report records a data-gate RALPH loop,
not a training benchmark, biological-age validation, or autoresearch run.

## Summary

The v11.x execution branch closed the downloaded processed-data conversion gate
for the current backlog.

Result:

- all direct Bismark-COV candidates that were allowed after smoke were converted;
- no candidate passed the headline data gate;
- no GroupKFold, held-out benchmark, CR validation, or autoresearch was run;
- two remaining candidates require dataset-specific adapters, but both are
  diagnostic/auxiliary by tissue context and are not expected to resolve the
  GSE121141 old brain_cortex/heart/lung blocker.

Decision:

`v11.x` did not find a new headline dataset from the downloaded backlog. The
current branch should stop before training. The next scientific step is a v12
data-strategy decision, unless we intentionally implement the two low-priority
adapters for auxiliary matrix coverage.

## RALPH Ledger Outputs

Controller:

- `scripts/validate/run_v11_x_conversion_gate.py`

Ledger directory:

- `results/ralph_v11_x_loop/`

Key files:

- `ralph_iteration_log.jsonl`
- `ralph_decision_state.json`
- `candidate_gate_table.csv`

Per-dataset state files:

- `results/v11_4_conversions/GSE224442/v11_x_GSE224442_conversion_gate_state.json`
- `results/v11_4_conversions/GSE286302/v11_x_GSE286302_conversion_gate_state.json`
- `results/v11_4_conversions/GSE92486/v11_x_GSE92486_conversion_gate_state.json`
- `results/v11_4_conversions/GSE129712/v11_x_GSE129712_conversion_gate_state.json`
- `results/v11_4_conversions/GSE175410/v11_x_GSE175410_conversion_gate_state.json`
- `results/v11_4_conversions/GSE281602/v11_4_gse281602_conversion_smoke_state.json`

## Candidate Gate Table

| Dataset | Status | Samples | Regions | Common with v8.2 | Exact age | Old target tissue | Decision |
|---|---:|---:|---:|---:|---|---|---|
| GSE224442 | auxiliary matrix | 78 | 146,915 | 65,827 | pass | fail | keep auxiliary |
| GSE281602 | auxiliary matrix | 8 | 121,333 | 25,321 | pass | fail | keep auxiliary |
| GSE286302 | auxiliary matrix | 49 | 111,808 | 63,524 | fail | fail | keep auxiliary |
| GSE92486 | auxiliary matrix | 18 | 123,457 | 6,500 | pass | fail | keep auxiliary |
| GSE129712 | blocked | - | - | - | not evaluated | fail by context | adapter required |
| GSE175410 | blocked | - | - | - | not evaluated | fail by context | adapter required |

## Direct Conversion Results

### GSE224442

Conversion succeeded.

- beta matrix: 3,486,147 CpGs x 78 samples
- 5kb region matrix: 146,915 regions x 78 samples
- common 5kb regions with v8.2 reference: 65,827
- beta range: `[0, 1]`
- metadata age coverage: 100%
- tissues: liver and blood
- intervention: parabiosis recovery

Gate decision:

`GSE224442` is technically useful because common-region overlap is high, but it
does not provide old brain_cortex/heart/lung support. It is auxiliary only and
must not enter the v11 headline chronological benchmark.

### GSE281602

Conversion succeeded in the first v11.4 smoke.

- beta matrix: 2,779,291 CpGs x 8 samples
- 5kb region matrix: 121,333 regions x 8 samples
- common 5kb regions with v8.2 reference: 25,321
- beta range: `[0, 1]`
- metadata age coverage: 100%
- tissue: heart
- context: cardiomyocyte/cell-type specific

Gate decision:

`GSE281602` is not bulk heart and fails the common-region threshold. It remains
auxiliary parser evidence only.

### GSE286302

Conversion succeeded after the large matrix merge.

- beta matrix: 15,820,781 CpGs x 49 samples
- 5kb region matrix: 111,808 regions x 49 samples
- common 5kb regions with v8.2 reference: 63,524
- beta range: `[0, 1]`
- metadata age coverage: 0%
- tissues: liver, skeletal_muscle, lung
- intervention/context: circadian recovery

Gate decision:

`GSE286302` has adequate technical overlap, but lacks sample-specific exact age
and does not provide old target-tissue support for the GSE121141 blocker. It is
auxiliary only.

### GSE92486

Conversion succeeded.

- beta matrix: 1,089,436 CpGs x 18 samples
- 5kb region matrix: 123,457 regions x 18 samples
- common 5kb regions with v8.2 reference: 6,500
- beta range: `[0, 1]`
- metadata age coverage: 100%
- tissue: liver
- intervention: control and diet restriction

Gate decision:

`GSE92486` is useful as a small liver diet-restriction auxiliary dataset, but it
fails the common-region threshold and does not target brain_cortex/heart/lung.

## Adapter-Required Candidates

### GSE129712

Smoke schema:

`methylratio_cg_12col`

Required adapter:

`methylratio_cg_12col_to_cov`

Gate decision:

This is a 26-month intestine diagnostic dataset. It is not a headline
brain_cortex/heart/lung candidate. Adapter implementation can be deferred unless
we need intestine auxiliary coverage.

### GSE175410

Smoke schema:

`bismark_cov_6col_accession_chrom`

Required adapter:

`bismark_cov_with_accession_chrom_map`

Gate decision:

This is a 24-month skeletal_muscle diagnostic dataset. It can be useful for
muscle aging coverage, but it is not a headline candidate for the GSE121141
brain_cortex/heart/lung failure.

## Success / Failure Gate Assessment

v11.x headline success requires a downloaded backlog dataset to satisfy:

- bulk old target tissue: brain_cortex, heart, or lung;
- sample-specific exact age;
- parseable methylation matrix;
- assembly-compatible 5kb region matrix;
- common regions `>= 50,000`;
- no leakage-risk training shortcut.

No candidate satisfied all gates.

Controlled v11.x failure condition is met for the downloaded backlog:

- direct-conversion candidates either lack old target tissue support, fail common
  regions, or lack exact age;
- adapter-required candidates are auxiliary by tissue context;
- no candidate justifies GroupKFold, GSE121141 held-out, CR validation, or
  autoresearch.

## Next Step

Do not run training on these v11.x matrices as headline data.

Recommended next action:

1. Start v12 as a data acquisition strategy loop focused only on official,
   sample-specific, old bulk brain_cortex/heart/lung mouse methylome candidates.
2. Keep `GSE224442` and `GSE286302` as auxiliary matrix resources because their
   common-region overlap is technically strong.
3. Defer `GSE129712` and `GSE175410` adapters unless auxiliary intestine/muscle
   coverage becomes a concrete v12 need.
4. If no official processed old target-tissue dataset is found in v12, define a
   minimal FASTQ/Bismark pilot for a small P3 candidate instead of tuning models
   on the current matrix.


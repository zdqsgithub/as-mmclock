# Current Data Exhaustion and Route A vs Deep Learning Assessment

Date: 2026-05-19

## Executive Judgment

Within the project-registered public data search space, the existing data have
been deeply explored enough to make a strategic decision:

- More training, autoresearch, or deep learning on the current matrices is not
  expected to solve the target scientific problem.
- The highest-probability path to the original goal is Route A: obtain or
  generate traceable old bulk `brain_cortex/heart/lung` mouse RRBS/WGBS data
  with exact sample-specific age and processed methylation files.
- Route B remains useful as a low-cost public-data rescue path, but v16 found
  only one priority candidate (`GSE213628`) and it still requires a processed
  adapter/matrix gate. It is not currently Learn-ready.

This is not a claim that every possible public dataset worldwide has been
exhausted. It is a claim that the current official GEO/E-Utils/SOFT/FTP,
SRA/ENA, v9-v16 RALPH search space no longer supports blind model optimization.

## Evidence From Project Gates

### v13 Benchmark Redefinition

The accepted v13 state is Route C: model claims are limited to support-covered
chronological-age prediction.

Key metrics:

- support-covered headline MAE: `20.286w`
- unsupported stress-test MAE: `36.659w`
- `GSE121141 old104+ brain_cortex/heart/lung` stress MAE: `88.565w`
- best observed CR AUC: `0.8744`, research-level held-out validation only

The old104+ target-tissue failure remains a stress-test/blocker, not a headline
pass/fail benchmark for the current model.

### v15 Matrix-Gate Rescue

v15 reviewed existing processed/raw rescue attempts and triggered the 3-strike
rule:

- no ready candidates;
- two auxiliary matrices (`GSE224442`, `GSE286302`) could not become headline;
- `GSE83947` processed and raw Route B pilot failed matrix gate;
- next action was to return to official refresh or Route A intake, not train.

The key blocker was not model architecture; it was lack of a candidate that
simultaneously passed exact age, old target tissue, bulk context, assembly,
beta range, and common 5kb region gates.

### v16 Route B Candidate Refresh

v16 refreshed `70` public candidates through official metadata/API routes:

- `P1_processed_headline`: `1`
- `P2_auxiliary`: `13`
- `P3_raw_pilot`: `0`
- `P4_blocked`: `56`
- Learn-ready candidates: `0`
- priority pilot candidates: `1`

The only priority candidate is:

| candidate | tier | status | old target support | schema | note |
| --- | --- | --- | --- | --- | --- |
| `GSE213628` | `P1_processed_headline` | `processed_adapter_smoke_or_matrix_gate_required` | `heart;lung`, `12` old exact samples | `bismark_cov_per_sample_tar` | needs processed adapter/matrix gate; not training-ready |

`GSE83947` is explicitly locked as `P2_auxiliary / do_not_repilot_headline`
because both processed and raw pilots already failed the common-region gate.

## Route A vs Deep Learning

### Route A: Highest Probability Path

Route A directly addresses the dominant blocker: absent old bulk
same-target-tissue support.

Required incoming data:

- mouse bulk RRBS/WGBS/bisulfite methylation;
- exact sample-specific age, preferably young/mid/old with old `>=104w`;
- `brain_cortex/cortex`, `heart`, and/or `lung`;
- sex, strain, intervention, assay, assembly, and sample provenance;
- processed methylation files or coverage/beta tables with checksums;
- enough overlap to pass `>=50,000` common 5kb regions.

Estimated success probability relative to project goal: high, if the data are
designed against the gate rather than collected opportunistically.

### Route B: Still Worth a Small Gate, Not a Main Strategy

Route B should continue only as one-candidate-at-a-time gated pilot work.
Current next candidate is `GSE213628` processed COV smoke/matrix gate.

Route B is lower probability because public datasets tend to fail one of:

- no exact sample-specific age;
- wrong tissue or adjacent tissue only;
- sorted/cell-type/single-cell context;
- assembly or schema mismatch;
- common 5kb region overlap below gate;
- age/tissue support not aligned to `GSE121141 old104+`.

### Deep Learning: Not the Right Next Main Move

Deep learning can be useful later for representation learning or feature
interpretation, but under the current data state it is not expected to fix the
core target:

- It cannot synthesize missing old bulk `brain_cortex/heart/lung` support.
- It is more likely to learn dataset, tissue, coverage, or assembly artifacts
  when the held-out stress set is outside support.
- It would violate the current AS-DS-Ops rule unless a new RFC explicitly
  authorizes it after data gates pass.
- Existing v4 embedding work already showed that representation methods should
  be gated and interpreted, not treated as a shortcut around missing data.

Current recommended role for deep learning:

- not headline training;
- possible future POC after Route A or a successful Route B matrix gate;
- use only with fixed held-out benchmarks, random-label sanity, and feature
  interpretation requirements;
- keep region-level interpretable features as the default scientific baseline.

## Decision

The project should not spend more cycles on broad autoresearch or deep learning
against the current matrix set.

Priority order:

1. Run the `GSE213628` processed COV adapter/matrix gate only if explicitly
   authorized as the next Route B pilot.
2. In parallel or next, pursue Route A data generation/collaboration using the
   existing Route A submission checklist and gate workflow.
3. Treat deep learning as a later RFC, only after a new headline matrix passes
   metadata, schema, assembly, common-region, and sanity gates.

## Current Authorized State

- download authorized: `false`
- Bismark authorized: `false`
- training authorized: `false`
- autoresearch authorized: `false`
- current model claim: support-covered chronological-age prediction only
- `GSE121141 old104+ brain_cortex/heart/lung`: stress-test/blocker metric

## Next Work Package

The most useful next engineering work is to prepare the `GSE213628` processed
adapter/matrix pilot as a guarded Route B action:

1. select 2-3 processed `.bismark.cov.gz` files from the v16 manifest;
2. download only those small processed files after explicit approval;
3. validate schema, beta range, coordinates, assembly, and sex/MT exclusion;
4. build a tiny region overlap estimate against the v8.2 reference;
5. stop immediately if common 5kb regions are `<50,000`;
6. if it passes, promote to full processed matrix gate, still without training
   until explicit Learn/Benchmark approval.


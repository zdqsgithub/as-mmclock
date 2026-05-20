# Route A Old Target-Tissue Data Generation RFC

Date: 2026-05-19T08:37:32

## Summary

Route A is the preferred path after v13 Route C acceptance and the Route B backlog refresh.
The current public/raw backlog did not produce a headline-ready minimal FASTQ pilot candidate.
This RFC defines the minimum data-generation or collaboration package needed to test true old target-tissue generalization.

This document does not authorize wet-lab work, data download, Bismark, model training, or autoresearch.

## Scientific Target

- Species: `Mus musculus`.
- Context: bulk tissue only.
- Target tissues: `brain_cortex`, `heart`, `lung`.
- Target old-age support: sample-specific exact age with old samples at `>=104w`.
- Preferred assay: bulk RRBS or WGBS with explicit genome assembly and processed methylation output.

## Cohort Design

- Minimum viable design: `72` samples.
- Preferred design: `96` samples.
- Preferred cell structure: 3 tissues x 4 age strata x 2 sexes x 4 replicates.
- Minimum cell structure: 3 tissues x 4 age strata x 2 sexes x 3 replicates.
- Old stratum target: `112w`, acceptable window `104-128w`.

The cohort design table is written to `results/route_a_data_generation_rfc/route_a_cohort_design.csv`.

## Required Metadata

The required schema is written to `results/route_a_data_generation_rfc/route_a_required_metadata_schema.csv`.
The key non-negotiable fields are `sample_id`, `mouse_id`, exact `age_days/age_weeks`, `tissue`, `sex`, `strain`, `intervention`, `assay`, `genome_assembly`, processed methylation path, and checksums.

## Gate Criteria

- Metadata overlap must be `>=95%`.
- Exact age coverage must be `>=95%`.
- At least `50,000` common 5kb regions with the v8.2/v13 reference must be available.
- Sex chromosomes and MT must be removable.
- Beta values must stay in `[0,1]`.
- Random-label sanity must pass before any model claim.
- CR or other biological-age metrics remain research-level and require real held-out predictions.

## Promotion Rule

Route A data may enter headline benchmarking only after metadata, schema, matrix, and sanity gates pass.
Constrained autoresearch remains forbidden unless the new data improves old104+ support by at least `5w` and sanity checks pass.

## Outputs

- `results/route_a_data_generation_rfc/route_a_cohort_design.csv`
- `results/route_a_data_generation_rfc/route_a_required_metadata_schema.csv`
- `results/route_a_data_generation_rfc/route_a_qc_benchmark_gates.json`
- `results/route_a_data_generation_rfc/route_a_decision_state.json`
- `doc/30_protocols/06_20260519_route_a_submission_checklist.md`

## Current Decision

`rfc_only_not_authorized`. Separate approval is required before any data generation, raw download, Bismark, training, or autoresearch.

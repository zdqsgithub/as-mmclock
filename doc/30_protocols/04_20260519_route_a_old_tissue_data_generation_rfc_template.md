# Route A RFC Template: Old Target-Tissue Data Generation / Collaboration

Date: 2026-05-19

## Summary

This RFC is for generating or collaborating on new old bulk mouse methylation
data to support full-lifespan brain_cortex/heart/lung chronological-age
benchmarking.

Completing this RFC does not authorize downloads, training, or ETL execution.

## Required Study Design

- Species: mouse.
- Assay: bulk RRBS, WGBS, or traceable bisulfite methylation.
- Tissues: at least one of `brain_cortex/cortex`, `heart`, or `lung`.
- Ages: young, mid, and old samples; old samples must include `>=104w`.
- Age metadata: sample-specific exact age, not age group only.
- Metadata: sample_id, age_days, age_weeks, tissue, sex, strain, intervention,
  dataset_batch, assay, source, and processing batch.
- Context: bulk tissue; not single-cell, organoid, cell-type targeted, or
  low-coverage/iTAG unless separately justified as auxiliary only.

## Matrix Requirements

- Assembly must be explicit and traceable.
- CpG coordinates must be convertible to the project's 5kb region framework.
- sex chromosomes and MT must be removable.
- Common 5kb regions with the v8.2/v13 reference must be expected to be
  `>=50000`.
- Coverage and missingness thresholds must be documented before training.

## Benchmark Gate

Promote to headline benchmark only if:

- metadata overlap is `>=95%`;
- age coverage is `>=95%`;
- target-tissue old samples are present;
- schema conversion passes adapter smoke;
- random-label sanity passes;
- no dummy CR/biological-age metric is used.

## Proposed Cohort

Fill in before approval:

- candidate owner or collaborator:
- tissue set:
- age design:
- sample count per tissue/age/sex:
- assay:
- expected assembly:
- expected data delivery format:
- ethical/data-use constraints:
- expected timeline:

## Approval Decision

- approved:
- reviewer:
- date:
- allowed actions:
- explicitly forbidden actions:

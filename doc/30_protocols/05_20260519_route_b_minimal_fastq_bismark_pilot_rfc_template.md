# Route B RFC Template: Minimal FASTQ/Bismark Pilot

Date: 2026-05-19

## Summary

This RFC is for a strictly limited 2-3 sample FASTQ/Bismark pilot when a raw
candidate may supply old bulk brain_cortex/heart/lung support but lacks usable
processed methylation files.

Completing this template does not authorize raw FASTQ download or Bismark
execution. Approval must be explicit.

## Pre-Approval Gate

All must be true before any pilot:

- BioSample/RunInfo confirms sample-specific exact age.
- Tissue is bulk `brain_cortex/cortex`, `heart`, or `lung`.
- Assay is bisulfite methylation compatible with the project.
- Run layout, read length, and total size are practical for a 2-3 sample pilot.
- Candidate is not single-cell, organoid, cell-type targeted, low-coverage/iTAG,
  or an untraceable superseries headline.
- SRA title alone is not used as proof of age or tissue.

## Pilot Scope

- Maximum samples: `3`.
- Goal: verify metadata, alignment feasibility, methylation extraction, assembly,
  coordinate compatibility, 5kb region overlap, and matrix join feasibility.
- Output is diagnostic only and cannot be used as headline benchmark evidence.

## Required Pilot Outputs

- download manifest with URLs, run accessions, sizes, and checksums;
- BioSample/RunInfo metadata snapshot;
- Bismark command manifest;
- alignment/extraction QC;
- CpG/region overlap estimate;
- blocker or promotion decision.

## Stop Conditions

Stop immediately if:

- age is age-group only;
- tissue is not target tissue;
- assay or layout is incompatible;
- common 5kb region estimate is below `50000`;
- pilot resource use exceeds the approved scope.

## Candidate Details

Fill in before approval:

- accession:
- BioProject:
- BioSample IDs:
- run accessions:
- tissue:
- age token and exact age:
- assay:
- layout/read length:
- expected size:
- selected 2-3 pilot samples:

## Approval Decision

- approved:
- reviewer:
- date:
- allowed download scope:
- allowed compute scope:
- explicitly forbidden actions:

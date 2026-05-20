# v13 Data Strategy RFC

Date: 2026-05-19

## Summary

v13 starts only after v12 benchmark redefinition and v12.1 support-covered evaluation. It is a data strategy decision, not model tuning.

## Route A: Generate or Collaborate

- Acquire new old bulk brain_cortex/heart/lung RRBS/WGBS data with sample-specific exact age.
- Minimum useful design: same tissue young/mid/old, with old samples at or above 104 weeks.
- This is the preferred route if the project needs a true full-lifespan target-tissue clock.

## Route B: Minimal FASTQ/Bismark Pilot

- Only allowed for raw candidates with BioSample/RunInfo-confirmed sample-specific age, bulk target tissue, bisulfite assay, and practical run size.
- Pilot size: 2-3 samples only.
- Pilot goal: verify alignment, methylation extraction, 5kb common regions, and metadata traceability before any full ETL.

## Route C: Accept Redefined Benchmark

- Use support-covered chronological benchmark as the headline result.
- Report GSE121141 old104+ as a stress-test/blocker metric.
- Keep biological-age/CR claims limited to real held-out predictions and shuffled-intervention sanity checks.

## Non-Negotiable Rules

- No autoresearch until a new data gate passes or the redefined benchmark is formally accepted.
- No dummy AUC.
- No human clock CpG mapping.
- No raw FASTQ download without explicit minimal ETL approval.

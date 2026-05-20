# Route A Submission Checklist

Date: 2026-05-19

Use this checklist before accepting generated or collaborative old target-tissue methylation data.

## Pre-Acceptance

- Confirm mouse bulk tissue context.
- Confirm exact sample-specific age for every sample.
- Confirm target tissue labels are one of `brain_cortex`, `heart`, `lung`.
- Confirm old samples include `>=104w` for each target tissue intended for headline claims.
- Confirm sex, strain, diet, intervention, batch, and collection metadata are present.
- Confirm assay is bulk RRBS/WGBS or traceable bisulfite methylation.
- Confirm genome assembly and coordinate system are explicit.
- Confirm processed methylation schema is documented.
- Confirm file checksums are supplied.

## Adapter Smoke

- Parse the first 1,000-10,000 rows or 2-3 sample files.
- Verify coordinates, beta range, sample IDs, and metadata join.
- Verify sex/MT exclusion can be applied.
- Estimate common 5kb region count against v8.2/v13 reference.
- Stop before training if common regions are `<50,000`.

## Benchmark Promotion

- Run fixed v7/v8 configs before any search.
- Run GroupKFold, LODO, GSE121141 held-out, GSE80672 CR held-out, random-label sanity, and shuffled CR sanity when applicable.
- Keep GSE121141 old104+ stress-test reporting even if Route A improves support.
- Do not run autoresearch unless new-data gates and sanity checks pass.

## Not Authorized By This Checklist

- Wet-lab execution.
- Raw FASTQ download.
- Bismark/FASTQ ETL.
- Model training.
- Autoresearch.

# Route A Partner Handoff Package

Date: 2026-05-19T08:41:58

## Purpose

This package converts the Route A RFC into files that can be sent to a collaborator or experimental team before any data generation or transfer.
It is a metadata and planning package only.

## Included Files

- `metadata/templates/route_a_sample_sheet_minimum_72_template.csv`
- `metadata/templates/route_a_sample_sheet_preferred_96_template.csv`
- `metadata/templates/route_a_file_manifest_minimum_72_template.csv`
- `metadata/templates/route_a_file_manifest_preferred_96_template.csv`
- `scripts/validate/validate_route_a_submission.py`
- `results/route_a_handoff_pack/route_a_handoff_manifest.json`
- `results/route_a_handoff_pack/route_a_template_validation_report.json`

## How To Use

1. Fill the sample sheet with exact age, tissue, sex, strain, intervention, batch, assay, genome assembly, and processed methylation metadata.
2. Fill the file manifest with processed coverage/beta file paths and checksums.
3. Run the validator before accepting the submission:

```bash
VIRTUAL_ENV=.venv-core uv run --active python scripts/validate/validate_route_a_submission.py \
  --sample-sheet metadata/templates/route_a_sample_sheet_minimum_72_template.csv \
  --file-manifest metadata/templates/route_a_file_manifest_minimum_72_template.csv
```

For a template dry check only, add `--allow-placeholders`.

## Current Template Scope

- Minimum samples: `72`.
- Preferred samples: `96`.
- Target tissues: `brain_cortex`, `heart`, `lung`.
- Age strata: `12w`, `52w`, `78w`, `112w`.
- Old-age support threshold: `>=104w`.

## Guardrails

- This package does not authorize wet-lab execution.
- This package does not authorize raw FASTQ download.
- This package does not authorize Bismark/FASTQ ETL.
- This package does not authorize model training or autoresearch.
- Any real submission must pass the validator and then adapter/matrix gates before benchmarking.

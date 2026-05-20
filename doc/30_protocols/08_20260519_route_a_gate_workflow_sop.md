# Route A Gate Workflow SOP

Date: 2026-05-19

This SOP controls how generated or collaborative Route A data enters the project.

## State Machine

1. Metadata gate: run `validate_route_a_submission.py` on the sample sheet and file manifest.
2. Adapter smoke: only after metadata passes, inspect 1-3 local processed methylation files or the first 1,000-10,000 rows.
3. Matrix gate: only after adapter smoke passes, check local path existence, manifest `file_size_bytes`, manifest `sha256`, common 5kb regions, metadata overlap, age coverage, sex/MT exclusion, and beta range.
4. RALPH Learn readiness: only after matrix gate passes, prepare fixed benchmark commands. Training still requires explicit approval.

## Default Command

Build a local file manifest from a filled sample sheet after files are staged:

```bash
VIRTUAL_ENV=.venv-core uv run --active python scripts/validate/build_route_a_local_file_manifest.py \
  --sample-sheet <filled_sample_sheet.csv> \
  --submission-id <submission_id> \
  --output results/route_a_local_file_manifest/<submission_id>/route_a_file_manifest.csv
```

Preferred one-command runner for a real filled submission. If `--file-manifest` is omitted, the runner builds a local manifest from the sample sheet first:

```bash
VIRTUAL_ENV=.venv-core uv run --active python scripts/validate/run_route_a_intake_to_learn_ready.py \
  --sample-sheet <filled_sample_sheet.csv> \
  --submission-id <submission_id>
```

Pass `--file-manifest <filled_file_manifest.csv>` only when using a pre-built manifest. This runner executes the same gate order below and stops before training. It only prepares a guarded RALPH Learn command package when all gates pass.

Metadata and adapter gate:

```bash
VIRTUAL_ENV=.venv-core uv run --active python scripts/validate/run_route_a_gate_workflow.py \
  --sample-sheet <filled_sample_sheet.csv> \
  --file-manifest <filled_file_manifest.csv>
```

Build matrix manifest after adapter smoke passes:

```bash
VIRTUAL_ENV=.venv-core uv run --active python scripts/etl/22_build_route_a_matrix.py \
  --sample-sheet <filled_sample_sheet.csv> \
  --file-manifest <filled_file_manifest.csv> \
  --submission-id <submission_id> \
  --reference-matrix results/multidataset_v8_2_prefilter_liftover/all_rrbs_region_matrix_5kb.parquet \
  --out-dir results/route_a_intake/<submission_id>
```

Matrix gate and RALPH Learn readiness:

```bash
VIRTUAL_ENV=.venv-core uv run --active python scripts/validate/run_route_a_gate_workflow.py \
  --sample-sheet <filled_sample_sheet.csv> \
  --file-manifest <filled_file_manifest.csv> \
  --matrix-manifest results/route_a_intake/<submission_id>/route_a_matrix_manifest.json \
  --output-dir results/route_a_intake/<submission_id>
```

Prepare a guarded RALPH Learn command package after readiness:

```bash
VIRTUAL_ENV=.venv-core uv run --active python scripts/validate/run_route_a_ralph_learn_package.py \
  --intake-dir results/route_a_intake/<submission_id> \
  --base-matrix results/multidataset_v8_2_prefilter_liftover/all_rrbs_region_matrix_5kb.parquet \
  --base-manifest results/multidataset_v8_2_prefilter_liftover/all_rrbs_matrix_manifest.json \
  --out-dir results/route_a_ralph_learn_package/<submission_id>
```

The generated shell script must remain guarded and must not run training until explicit approval is given.

Use `--allow-placeholders` only for template dry-runs.

Contract regression test:

```bash
VIRTUAL_ENV=.venv-core uv run --active python scripts/validate/test_route_a_gate_contracts.py
```

This test reruns a synthetic gate smoke, verifies the local manifest generator and runner auto-manifest mode, verifies the matrix thresholds, checks that remote paths, bad beta values, and checksum mismatches are blocked, and confirms the generated command script exits before any training command.

## Guardrails

- No download is authorized by this workflow.
- Remote URIs are recorded as pending; files must be staged separately before adapter smoke.
- No Bismark/FASTQ ETL is authorized.
- No training or autoresearch is authorized.
- If any gate fails, stop and fix data/metadata before moving forward.

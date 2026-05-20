# v13 Delivery Freeze Report

Date: 2026-05-19

## Summary

This report freezes the current v13 project state as a reproducible delivery
baseline. The freeze uses `uv` for Python environment management and does not
run model training, downloads, FASTQ/Bismark ETL, matrix rebuilds, or
autoresearch.

The accepted project route is **Route C: Accept Redefined Benchmark**.

## Environment

- Environment manager: `uv 0.11.7`
- Project interpreter: `.venv/bin/python3`
- Python version: `3.13.13`
- Dependency snapshot: `results/v13_delivery_freeze/environment_manifest.json`
- Full artifact inventory: `results/v13_delivery_freeze/artifact_inventory.csv`

The default v13 delivery reproduction now uses `requirements-core.txt`, not the
full training stack. This avoids installing torch/lightning/CUDA packages for
benchmark redefinition and decision-package reproduction.

Dependency split:

- `requirements-core.txt`: v12.1/v13 delivery reproduction.
- `requirements-viz.txt`: optional plotting/reporting dependencies.
- `requirements-full.txt`: training, deep-learning, and full historical
  environment.
- `requirements.txt`: compatibility entrypoint that references
  `requirements-full.txt`.

## Reproduction Commands

```bash
uv venv .venv
uv venv .venv-core
VIRTUAL_ENV=.venv-core uv pip install -r requirements-core.txt
VIRTUAL_ENV=.venv-core uv run --active python scripts/validate/run_v12_1_redefined_benchmark.py
VIRTUAL_ENV=.venv-core uv run --active python scripts/validate/run_v13_strategy_decision.py
VIRTUAL_ENV=.venv-core uv run --active python scripts/validate/run_v13_delivery_freeze.py
VIRTUAL_ENV=.venv-core uv run --active python scripts/validate/test_v13_delivery_contracts.py
```

## Frozen Outputs

- `results/benchmark_v12_1_redefined/`
- `results/ralph_v13_strategy/`
- `results/v13_delivery_freeze/environment_manifest.json`
- `results/v13_delivery_freeze/artifact_inventory.csv`
- `results/v13_delivery_freeze/reproducibility_manifest.json`
- `results/v13_delivery_freeze/v13_delivery_summary.md`
- `results/v13_delivery_freeze/test_report.json`
- `results/v13_delivery_freeze/test_report.txt`

Core environment validation:

- `.venv-core` size: approximately `245M`.
- Full `.venv` size from the previous full install: approximately `4.9G`.
- Core environment contains pandas/numpy/scipy/scikit-learn and does not contain
  torch/lightning.

## Model Scope

Current headline claims are limited to support-covered chronological-age
prediction. A prediction is support-covered when train data contains at least
10 same-tissue samples and held-out age is no more than 8 weeks above the train
same-tissue maximum age.

Current headline metrics:

- Support-covered headline MAE: `20.286` weeks over `15573` rows.
- Unsupported stress-test MAE: `36.659` weeks over `7620` rows.
- GSE121141 old104+ stress MAE: `88.565` weeks over `420` rows.
- Best observed CR AUC: `0.8744`; research-level held-out validation only.

## Stress-Test Boundary

`GSE121141 old104+ brain_cortex/heart/lung` remains unresolved. It must be
reported as a stress-test/blocker metric and cannot be used as the current
headline pass/fail benchmark.

Unsupported rows remain visible in reports, but they must not drive
autoresearch or model claims.

## Delivery Tests

The v13 delivery contract tests passed:

- `run_v13_strategy_decision.py` reproduces through `uv run python`.
- v13 decision state selects only Route C.
- raw FASTQ download, training, and autoresearch are all unauthorized.
- support-covered rows obey same-tissue sample count and age-gap rules.
- CR metrics come from real held-out prediction files and are not dummy AUC.

Test outputs:

- `results/v13_delivery_freeze/test_report.json`
- `results/v13_delivery_freeze/test_report.txt`

## Next Routes

Route A remains the preferred path for a true full-lifespan target-tissue mouse
clock: generate or collaborate on old bulk brain_cortex/heart/lung RRBS/WGBS
data with sample-specific age.

Route B remains blocked until a separate minimal FASTQ/Bismark pilot approval
identifies a valid raw candidate with BioSample/RunInfo-confirmed sample age,
target tissue, assay, layout, and practical run size.

Route C is accepted for the current project state: future model use must stay
inside the support-covered benchmark boundary unless Route A or Route B supplies
new valid data.

## Guardrails

- No autoresearch from unsupported old104+ metrics.
- No raw FASTQ download without separate minimal ETL approval.
- No Bismark/FASTQ ETL in this delivery freeze.
- No dummy AUC.
- No human clock CpG mapping.
- No claims of full-lifespan old target-tissue generalization.

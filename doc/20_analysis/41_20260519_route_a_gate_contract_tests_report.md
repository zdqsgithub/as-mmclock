# Route A Gate Contract Tests Report

Date: 2026-05-19

## Summary

Route A gate/readiness behavior is now covered by a dedicated unittest smoke
suite. The suite does not train, download, run Bismark, or start autoresearch.

## Test Entry

```bash
VIRTUAL_ENV=.venv-core uv run --active python scripts/validate/test_route_a_gate_contracts.py
```

## Result

- status: `passed`
- tests_run: `8`
- failures: `0`
- errors: `0`

Output files:

- `results/route_a_contract_tests/test_report.json`
- `results/route_a_contract_tests/test_report.txt`

## Covered Contracts

- The intake-to-learn-ready runner reaches
  `ready_for_ralph_learn_pending_explicit_training_approval` on the synthetic
  fixture without authorizing training.
- The Route A matrix manifest requires `common_regions >= 50000`, metadata
  overlap `>=0.95`, age coverage `>=0.95`, and beta range `[0,1]`.
- The local file manifest generator computes `file_size_bytes` and `sha256`
  from staged files and can feed the intake runner.
- The intake runner can auto-build the local file manifest when `--file-manifest`
  is omitted.
- Remote URIs are blocked by the local file manifest generator.
- Local processed methylation files must match manifest `file_size_bytes` and
  `sha256`; checksum mismatch blocks matrix construction.
- The guarded RALPH Learn shell script exits before the first `train_clock.py`
  command.
- A bad beta fixture is blocked at adapter smoke and cannot reach readiness.

## Guardrails

- training_authorized: `False`
- download_authorized: `False`
- bismark_authorized: `False`
- autoresearch_authorized: `False`

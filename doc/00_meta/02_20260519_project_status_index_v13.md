# Mouse Methylation Clock Project Status Index v13

Date: 2026-05-19

## Current Status

The current accepted project state is v13 Route C: accept the redefined
support-covered chronological benchmark.

Current headline benchmark:

- support-covered chronological-age prediction only.

Required stress-test:

- `GSE121141 old104+ brain_cortex/heart/lung`, reported as stress-test/blocker,
  not as headline pass/fail.

## Phase Summary

- v1-v4: prototype clock, architecture repair, region-based features, and
  embedding-aware proof of concept.
- v5-v8: multidataset matrix work, GSE80672 CR validation, GSE60012 synthetic
  metadata, GSE213628, assembly-aware harmonization, and ablation diagnostics.
- v9-v12: RALPH data strategy loops. These confirmed that available public
  processed data did not supply traceable old bulk brain_cortex/heart/lung
  headline support.
- v12.1-v13: benchmark redefinition and Route C acceptance. The model scope is
  now bounded to support-covered chronological-age prediction.
- v14: public-data rescue and Route B minimal raw pilot. `GSE83947` old-lung
  FASTQ/Bismark was executable, but common 5kb overlap with v8.2 was only
  `1,814`, so it failed the `>=50,000` matrix gate and remains auxiliary only.

## Current Metrics

- Support-covered headline MAE: `20.286` weeks.
- Unsupported stress-test MAE: `36.659` weeks.
- GSE121141 old104+ stress MAE: `88.565` weeks.
- Best observed CR AUC: `0.8744`, research-level only.

## Entry Points

- AS-DS-Ops program: `/home/zdq-as/as-ds-ops/program.md`
- v13 model card: `doc/00_meta/01_20260519_v13_model_card.md`
- v12.1 benchmark report: `doc/20_analysis/31_20260519_v12_1_redefined_benchmark_report.md`
- v13 data strategy RFC: `doc/20_analysis/32_20260519_v13_data_strategy_rfc.md`
- v13 delivery freeze report: `doc/20_analysis/33_20260519_v13_delivery_freeze_report.md`
- Route B candidate audit report: `doc/20_analysis/34_20260519_route_b_candidate_audit_report.md`
- Route B manual backlog refresh report: `doc/20_analysis/35_20260519_route_b_manual_backlog_runinfo_refresh_report.md`
- Route A submission checklist: `doc/30_protocols/06_20260519_route_a_submission_checklist.md`
- Route A partner handoff package: `doc/30_protocols/07_20260519_route_a_partner_handoff_package.md`
- Route A intake gate readiness report: `doc/20_analysis/38_20260519_route_a_intake_gate_readiness_report.md`
- Route A RALPH Learn package report: `doc/20_analysis/39_20260519_route_a_ralph_learn_package_report.md`
- Route A intake-to-learn-ready runner report: `doc/20_analysis/40_20260519_route_a_intake_to_learn_ready_runner_report.md`
- Route A gate contract tests report: `doc/20_analysis/41_20260519_route_a_gate_contract_tests_report.md`
- Route A P0-P2 execution status: `doc/20_analysis/42_20260519_route_a_p0_p2_execution_status.md`
- Route A local file manifest builder: `scripts/validate/build_route_a_local_file_manifest.py`
- Route A gate workflow dry-run: `doc/20_analysis/37_20260519_route_a_gate_workflow_dry_run_report.md`
- Route A gate workflow SOP: `doc/30_protocols/08_20260519_route_a_gate_workflow_sop.md`
- Route A concrete data-generation RFC: `doc/20_analysis/36_20260519_route_a_old_tissue_data_generation_rfc.md`
- v14 public-data rescue report: `doc/20_analysis/43_20260519_v14_public_data_rescue_report.md`
- v14 GSE83947 Route B raw pilot preflight/report: `doc/20_analysis/44_20260519_v14_gse83947_route_b_raw_pilot_preflight.md`
- research-grade raw-to-interpretable-clock training plan: `doc/10_design/02_20260519_research_grade_raw_methylation_training_plan.md`
- v15 matrix-gate RALPH loop report: `doc/20_analysis/45_20260519_v15_matrix_gate_ralph_loop_report.md`
- v16 Route B candidate refresh report: `doc/20_analysis/46_20260519_v16_route_b_candidate_refresh_report.md`
- current data exhaustion and Route A vs deep learning assessment: `doc/20_analysis/47_20260519_current_data_exhaustion_route_a_vs_deep_learning_assessment.md`
- v17 demo ML best baseline report: `doc/20_analysis/48_20260520_demo_ml_best_baseline_report.md`
- v17.1 demo ML narrow optimization report: `doc/20_analysis/49_20260520_demo_ml_narrow_optimization_report.md`
- v20 AutoDL DL autoresearch report: `doc/20_analysis/51_20260520_v20_autodl_deep_learning_autoresearch_report.md`
- v21 raw-to-interpretable clock plan: `doc/20_analysis/52_20260521_v21_raw_to_interpretable_clock_plan.md`
- v21 raw-to-interpretable clock SOP: `doc/30_protocols/09_20260521_raw_to_interpretable_clock_sop.md`
- v21 raw pilot execution report: `doc/20_analysis/55_20260521_v21_raw_pilot_execution_report.md`
- v21 RAID raw RALPH controller: `scripts/validate/run_v21_raid_raw_ralph_loop.py`
- v21 raw ETL queue/watchdog: `scripts/etl/23_run_raw_bismark_queue.py`, `scripts/etl/24_raw_etl_watchdog.py`
- v21 AutoDL export/import: `scripts/train/export_v21_autodl_package.py`, `scripts/train/import_v21_autodl_results.py`
- v21 feature interpretation: `scripts/validate/run_v21_feature_interpretation.py`
- v13 delivery outputs: `results/v13_delivery_freeze/`
- v13 route decision: `results/ralph_v13_strategy/v13_route_decision_state.json`

## Environment

Use `uv` for Python environment management.

Default delivery reproduction:

```bash
uv venv .venv-core
VIRTUAL_ENV=.venv-core uv pip install -r requirements-core.txt
VIRTUAL_ENV=.venv-core uv run --active python scripts/validate/run_v12_1_redefined_benchmark.py
VIRTUAL_ENV=.venv-core uv run --active python scripts/validate/run_v13_strategy_decision.py
VIRTUAL_ENV=.venv-core uv run --active python scripts/validate/run_v13_delivery_freeze.py
VIRTUAL_ENV=.venv-core uv run --active python scripts/validate/test_v13_delivery_contracts.py
```

Use `requirements-full.txt` only for training/deep-learning/full historical
environments.

v21 raw and interpretation tooling uses the project `uv` environment. Local ML
is CPU-only (`ridge`, `elasticnet`, `lgbm`) after matrix gate and explicit
training approval. DL runs are packaged for AutoDL only.

## Next Strategy Routes

Route A remains the preferred path for a true full-lifespan target-tissue mouse
clock: generate or collaborate on old bulk brain_cortex/heart/lung RRBS/WGBS
data with sample-specific age.

Route B was attempted for `GSE83947` after explicit approval. The environment,
FASTQ download, checksum, Bismark alignment, and coverage extraction passed, but
the matrix gate failed (`1,814` common 5kb regions with v8.2, below `50,000`).
Future Route B work needs a different candidate or a revised data-generation RFC;
`GSE83947` does not authorize training or autoresearch.

v21 raw pilots have now passed matrix gates for `GSE121141` and `GSE80672`.
This validates the local FASTQ->Bismark COV->5kb matrix ETL path for two
priority datasets, but the current raw pilot matrix is still too small for
formal ML, CR, random-label, AutoDL, or biological-interpretation success
claims. The next v21 step is age/condition-balanced raw expansion before
autoresearch.

## Guardrails

- No autoresearch from unsupported old104+ metrics.
- No raw FASTQ download without separate minimal ETL approval.
- No Bismark/FASTQ ETL without approval.
- No full-scale v21 raw expansion before a 2-3 sample pilot passes metadata,
  age, beta, sex/MT, and common-region gates.
- No local GPU/DL execution; v21 DL artifacts must come from an AutoDL package
  or another explicitly approved GPU host.
- No biological mechanism claim from a top region unless confounding and
  coverage audits support that interpretation.
- No dummy AUC.
- No human clock CpG mapping.
- No claims of full-lifespan old target-tissue generalization.

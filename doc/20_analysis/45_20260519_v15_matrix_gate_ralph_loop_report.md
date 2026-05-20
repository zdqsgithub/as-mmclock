# v15 Matrix-Gate RALPH Loop Report

Date: 2026-05-19T05:01:57.330013+00:00

## Summary

This v15 controller replays the project matrix-gate evidence under the current
AS-DS-Ops rules. It did not download data, run Bismark, train models, or start
autoresearch.

Decision: `three_strike_re_evaluate_no_current_headline_matrix_gate_pass`.

Target state for success:

- `ready_for_ralph_learn_pending_explicit_training_approval`
- sample-specific exact age coverage `>=95%`
- old bulk target tissue: `brain_cortex/cortex`, `heart`, or `lung`, including `>=104w`
- parseable methylation schema and traceable assembly
- metadata overlap `>=95%`
- beta range within `[0,1]`
- common 5kb regions with v8.2 reference `>=50,000`
- no training/autoresearch until explicit approval

## Attempt Table

| attempt_id | dataset | stage | n_samples | common_regions_with_v8_2 | sample_specific_age_pass | target_tissue_pass | bulk_context_pass | common_region_gate | headline_matrix_gate_pass | gate_status | blocker_type |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gse134398_candidate_table | GSE134398 | candidate_table_research_audit | 158 | 0 | False | True | False | False | False | research_or_audit_only_not_matrix_ready | sample_specific_exact_age_gate_failed;bulk_context_gate_failed;schema_smoke_failed;assembly_not_traceable;common_regions_lt_50000 |
| gse224442_processed_cov_conversion | GSE224442 | processed_cov_conversion_gate | 78 | 65827 | True | False | True | True | False | auxiliary_matrix_gate_passed_not_headline | old_target_tissue_gate_failed |
| gse225166_candidate_table | GSE225166 | candidate_table_research_audit | 111 | 0 | True | True | False | False | False | research_or_audit_only_not_matrix_ready | bulk_context_gate_failed;schema_smoke_failed;assembly_not_traceable;common_regions_lt_50000 |
| gse232547_candidate_table | GSE232547 | candidate_table_research_audit | 24 | 0 | False | False | True | False | False | research_or_audit_only_not_matrix_ready | sample_specific_exact_age_gate_failed;old_target_tissue_gate_failed;schema_smoke_failed;assembly_not_traceable;common_regions_lt_50000 |
| gse281602_candidate_table | GSE281602 | candidate_table_research_audit | 8 | 25321 | True | True | False | False | False | research_or_audit_only_not_matrix_ready | bulk_context_gate_failed;common_regions_lt_50000 |
| gse286302_processed_cov_conversion | GSE286302 | processed_cov_conversion_gate | 49 | 63524 | False | False | True | True | False | auxiliary_matrix_gate_passed_not_headline | sample_specific_exact_age_gate_failed;old_target_tissue_gate_failed |
| gse304754_candidate_table | GSE304754 | candidate_table_research_audit | 0 | 0 | True | False | True | False | False | research_or_audit_only_not_matrix_ready | old_target_tissue_gate_failed;schema_smoke_failed;assembly_not_traceable;common_regions_lt_50000 |
| gse83947_processed_cx_report_gate | GSE83947 | processed_adapter_matrix_gate | 3 | 539 | True | True | True | False | False | matrix_gate_failed | common_regions_lt_50000 |
| gse83947_route_b_raw_bismark_pilot | GSE83947 | route_b_raw_bismark_matrix_gate | 3 | 1814 | True | True | True | False | False | matrix_gate_failed | common_regions_lt_50000 |
| gse92486_processed_cov_conversion | GSE92486 | processed_cov_conversion_gate | 18 | 6500 | True | False | True | False | False | matrix_gate_failed | old_target_tissue_gate_failed;common_regions_lt_50000 |

## 3-Strike Re-evaluation

3-strike triggered: `True`.

| strike_number | attempt_id | dataset | common_regions_with_v8_2 | blocker_type | next_action |
| --- | --- | --- | --- | --- | --- |
| 1 | gse83947_processed_cx_report_gate | GSE83947 | 539 | common_regions_lt_50000 | do_not_train_research_or_new_route_a_route_b_candidate |
| 2 | gse83947_route_b_raw_bismark_pilot | GSE83947 | 1814 | common_regions_lt_50000 | do_not_train_research_or_new_route_a_route_b_candidate |
| 3 | gse281602_candidate_table | GSE281602 | 25321 | bulk_context_gate_failed;common_regions_lt_50000 | do_not_train_use_only_for_next_candidate_prioritization |

Interpretation:

- `GSE83947` is biologically attractive as old bulk lung with exact age, but
  both processed CX_report and raw Bismark pilot failed the common-region gate.
- `GSE281602` has exact old heart signal, but is cardiomyocyte/cell-type
  specific and does not pass the common-region gate.
- `GSE286302` and `GSE224442` show that matrix overlap can pass technically,
  but they fail headline metadata/tissue gates and therefore cannot solve the
  old target-tissue clock.

## Re-evaluation Rules

After three failed/blocked attempts, the next step is not training and not
autoresearch. The loop must return to Research/Audit using official metadata
sources:

- GEO Download / supplementary file listings
- GEO Programmatic Access / SOFT
- SRA RunInfo
- ENA read_run file reports

Only a new P1 processed candidate or separately approved P3 Route B candidate
may proceed to adapter smoke and matrix construction.

## Next RALPH Action

`return_to_official_research_refresh_or_route_a_intake_do_not_train`

Recommended practical route:

1. Prefer Route A generated/collaborative old bulk `brain_cortex/heart/lung`
   data, because current public candidates have not satisfied all gates.
2. If public-data rescue continues, run a fresh official candidate refresh
   before any download. Do not reuse `GSE83947` for headline matrix work.
3. For any new Route B candidate, repeat only a 2-3 sample pilot and stop if
   common regions remain `<50,000`.
4. If a real Route A sample sheet and local processed methylation manifest are
   available, run the Route A metadata -> adapter -> matrix gate workflow.

## Guardrails

- Training authorized: `false`
- Autoresearch authorized: `false`
- Download authorized by this controller: `false`
- Bismark authorized by this controller: `false`
- Human clock CpG mapping: forbidden
- Dummy AUC: forbidden

## Contract Test

Run after controller:

```bash
VIRTUAL_ENV=.venv-core uv run --active python scripts/validate/test_v15_matrix_gate_contracts.py
```


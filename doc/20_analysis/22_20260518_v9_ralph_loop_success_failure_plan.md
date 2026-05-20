# v9 RALPH Loop Success and Failure Plan

Date: 2026-05-18

## Summary

v1-v8 established that the current bottleneck is not model complexity. The
remaining failure mode is old-age same-tissue support for GSE121141
brain_cortex/heart/lung. v9 therefore changes the project from optimization to a
data-gated RALPH loop.

The first v9 dry run produced `results/ralph_v9_loop/` and confirmed the loop can
read existing v8.4-v8.7 outputs, generate a gate table, and decide the next step
without downloading large files or training models. A first official refresh was
then run with E-Utils/SOFT/FTP filelist/HEAD only; no supplement archives, matrices,
or training runs were started.

## Lessons from v1-v8

- v1-v2: metric integrity matters more than search size; fold-internal feature selection and true units are mandatory.
- v3: region features are the correct RRBS baseline because site overlap is fragile across studies.
- v4: embedding is useful for interpretation, but it cannot substitute for valid held-out data.
- v5-v6: real CR validation requires actual held-out predictions, not metadata or dummy AUC.
- v7: preprocessing/harmonization helps locally but does not fix old-age tissue mismatch.
- v8: adding available old-age datasets and calibration did not solve GSE121141 old104+ generalization.
- v8.7: GSE225166, GSE83947, and GSE134398 are downloadable and parseable in smoke tests, but remain auxiliary or blocked for headline use.

## v9 Decision Logic

v9 should continue only while it can find either:

- P1 headline data: bulk RRBS/WGBS processed methylation with sample-specific old target tissue.
- P3 minimal ETL data: metadata passes, but processed methylation is absent and a 2-3 sample FASTQ/Bismark pilot is justified.

v9 should stop when two official refresh rounds find no P1/P3 candidate, or when
valid headline matrices still leave GSE121141 old104+ MAE near `70w+`.

## Final v9.0 Result

- Output directory: `results/ralph_v9_loop/`
- Decision: `v9_failure_no_headline_data_after_two_refresh_rounds`
- Reason: two official refresh passes plus web-scout seed expansion produced no P1 headline candidate and no P3 minimal FASTQ pilot candidate.
- P1 headline candidates: 0
- P3 minimal FASTQ pilot candidates: 0
- P2 auxiliary candidates: 1 (`GSE83947`, lung-only)
- Important blocked candidates: `GSE225166` remains P4 because smoke/parser pass is outweighed by targeted cell-type/low-coverage tagged context; `GSE134398` remains P4 because of superseries/subseries traceability.
- Additional official GEO web scouting found no new P1 candidate. `GSE171236` is already in the seed list but is hippocampus/brain-other and lacks 24-month-plus target support in the current parser. `GSE215310`, `GSE156557`, and `GSE221124` were added as low-priority seeds for future official refresh, but their visible GEO summaries are liver/pancreas/spleen/cerebellum oriented rather than GSE121141-like brain_cortex/heart/lung headline data.
- Next action: stop v9 and move to v10/v9-follow-on data strategy: broader processed methylome acquisition or a minimal ETL plan for target-tissue old-age cohorts.

## Why Not Autoresearch Yet

Autoresearch is only meaningful after a data candidate passes matrix and sanity
gates. Running 30-50 configs on the current data would optimize against a known
support gap and risks selecting models that exploit batch or tissue artifacts.

## Required Outputs for Future Iterations

- `ralph_iteration_log.jsonl`
- `ralph_decision_state.json`
- `candidate_gate_table.csv`
- `network_resolution_log.jsonl`
- `candidate_smoke_manifest.csv`
- `v9_success_or_failure_report.md`

## Recommended Next Work

Do not run constrained autoresearch from the current v9 candidate set. The next
plan should define a broader acquisition strategy or minimal FASTQ/Bismark ETL
pilot for old target tissues, with the same no-leakage benchmark gates preserved.

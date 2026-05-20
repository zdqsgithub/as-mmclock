# v9 RALPH Loop SOP

Date: 2026-05-18

## Purpose

v9 is a data-strategy loop for the mouse methylation clock project. It exists to
resolve the v8 bottleneck: poor old-age GSE121141 brain_cortex/heart/lung
generalization. It is not an autoresearch or deep-learning loop.

## State Machine

1. R / Research
   - Use official E-Utils, GEO SOFT, GEO FTP filelist, HEAD requests, SRA and ENA metadata.
   - Do not download large archives.
   - Record network/API issues in `network_resolution_log.jsonl`.

2. A / Audit
   - Download only 2-3 small per-sample COV/bedGraph/TXT files when `--run_smoke` is explicitly enabled.
   - Verify sample-specific age, tissue, assembly, beta range, coverage, and smoke common-region overlap.
   - Stop at audit if age or tissue is inferred only from global protocol text.

3. L / Learn
   - Build a candidate matrix only after P1/P3 gates pass.
   - Require common 5kb regions `>=50000` before training.
   - Run fixed benchmarks only: GroupKFold, LODO, GSE121141 held-out, GSE80672 held-out, random-label sanity, shuffled CR sanity.

4. P / Promote or Stop
   - Promote to constrained autoresearch only after the old104+ MAE improves by at least 5 weeks and sanity checks pass.
   - Stop v9 if two official refresh rounds produce no P1/P3 candidate, or valid headline matrices still leave old104+ MAE near 70 weeks.

## Candidate Tiers

- P1 headline candidate: bulk RRBS/WGBS or equivalent processed methylation, sample-specific age, old target tissue, parseable schema, assembly/liftover path.
- P2 single-target auxiliary: one old target tissue only, useful for diagnostics but not headline proof.
- P3 minimal FASTQ pilot: metadata and tissue are valid but processed methylation is absent; only a 2-3 sample pilot is allowed before any larger ETL.
- P4 blocked auxiliary: single-cell, targeted cell type, low-coverage/iTAG, organoid/in vitro, HSC-only, or untraceable superseries contexts.

## Network and API Rules

- E-Utils/SOFT are metadata sources, not file download authorities.
- GEO sample-level supplements use `/geo/samples/GSM.../<GSM>/suppl/<filename>`.
- GEO series archives use `/geo/series/GSE.../<GSE>/suppl/`.
- SRA/ENA raw data is considered only after processed supplements fail and the P3 pilot gate is approved.
- Non-official web search may suggest candidates, but every candidate must be verified against official GEO/SRA/ENA metadata before entering the ledger.

## Success Definition

v9 succeeds only when all conditions are met:

- A headline or minimal-ETL dataset passes sample-specific age, old target tissue, parseable methylation, assembly/liftover, and common region gates.
- GSE121141 old104+ MAE improves by `>=10w` versus the v7.5 `75.386w` baseline.
- GSE121141 all-age held-out MAE is `<=40.033w`.
- GroupKFold MAE is `<=25.767w`.
- Random-label sanity has `abs(r)<0.2` and random-like MAE.
- GSE80672 shuffled CR sanity is chance-like and Cohen's d decreases.

## Failure Definition

v9 fails when any condition is met:

- Two official refresh rounds produce no P1/P3 candidate.
- Three valid headline candidates still leave GSE121141 old104+ MAE near `70w+`.
- Only auxiliary single-cell, targeted, low-coverage/iTAG, organoid, HSC, or untraceable superseries data is available.
- Calibration improves only the target dataset and fails non-target LODO transfer.
- Matrix common regions are `<50000` or sanity checks fail.

## Command Templates

Dry-run, no network downloads:

```bash
/home/zdq-as/as-ds-ops/.venv/bin/python scripts/validate/run_v9_ralph_loop.py --max_candidates 28
```

Official refresh, no supplement download:

```bash
/home/zdq-as/as-ds-ops/.venv/bin/python scripts/validate/run_v9_ralph_loop.py --refresh_candidates --max_candidates 28
```

Adapter smoke for candidates selected by the controller:

```bash
/home/zdq-as/as-ds-ops/.venv/bin/python scripts/validate/run_v9_ralph_loop.py --refresh_candidates --run_smoke --max_candidates 28
```

## Non-Negotiable Rules

- No RRBS deduplication.
- No human clock CpG mapping.
- No dummy CR-AUC.
- No full-data feature selection.
- No biological-age claim without real held-out predictions.
- No autoresearch until the data gate and sanity checks pass.

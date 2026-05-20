# v10 Broader Data Acquisition and Minimal ETL RALPH Report

Date: 2026-05-18

## Summary

v10 executed a broader RALPH data-strategy loop after v9 failed to find a
headline public processed dataset for old GSE121141-like target tissues. The
loop expanded official GEO/E-Utils/SOFT/FTP discovery and reclassified
candidates into P1-P5 tiers, but it did not download large archives, build new
matrices, run FASTQ/Bismark, train models, or start autoresearch.

Final decision: `v10_failure_no_p1_or_p3_after_broadened_search`.

## Execution

Controller:

- `scripts/validate/run_v10_ralph_loop.py`

Main command pattern:

```bash
/home/zdq-as/as-ds-ops/.venv/bin/python scripts/validate/run_v10_ralph_loop.py \
  --refresh_candidates \
  --refresh_round 2 \
  --max_candidates 60 \
  --retmax 90
```

Outputs:

- `results/ralph_v10_loop/ralph_iteration_log.jsonl`
- `results/ralph_v10_loop/ralph_decision_state.json`
- `results/ralph_v10_loop/candidate_gate_table.csv`
- `results/ralph_v10_loop/network_resolution_log.jsonl`
- `results/ralph_v10_loop/candidate_smoke_manifest.csv`
- `results/ralph_v10_loop/minimal_etl_pilot_manifest.csv`
- `results/ralph_v10_loop/v10_success_or_failure_report.md`

## Gate Counts

From `ralph_decision_state.json`:

| Metric | Value |
| --- | ---: |
| Total candidates evaluated | 60 |
| P1 headline exact target | 0 |
| P2 adjacent old-age auxiliary | 7 |
| P3 raw exact target pilot | 0 |
| P4 low-coverage/targeted bridge | 11 |
| P5 blocked | 42 |

Because `P1=0` and `P3=0` after two broadened official refresh rounds, v10
meets the predefined failure gate.

## Notable Candidates

| Dataset | Tier | Evidence | Decision |
| --- | --- | --- | --- |
| GSE225166 | P4 | old heart/lung present, processed files present, but single-cell/targeted and low-coverage context; common-region smoke estimate 5088 | compatibility audit only |
| GSE83947 | P4 | old lung support, but low-coverage/tagged context and poor assembly/common-region signal | compatibility audit only |
| GSE134398 | P4/P5 context | old lung signal but superseries traceability risk | do not promote |
| GSE224442 | P2 | old blood/liver processed support | adjacent diagnostic only |
| GSE92486 | P2 | old liver processed support | adjacent diagnostic only |
| GSE129712 | P2 | old intestine support | adjacent diagnostic only |
| GSE281062 | P2 | old brain_other support | adjacent diagnostic only |
| GSE197045 | P2 | old skeletal_muscle support | adjacent diagnostic only |
| GSE175410 | P2 | old skeletal_muscle support | adjacent diagnostic only |
| GSE85772 | P2 | old liver support | adjacent diagnostic only |

None of these candidates satisfies the v10 headline requirement: bulk RRBS/WGBS
or minimal ETL-ready old brain_cortex/heart/lung with sample-specific age and
matrix path likely to yield common 5kb regions `>=50000`.

## Interpretation

The current bottleneck is no longer a training-architecture question. v1-v8
already fixed leakage, metadata, region aggregation, assembly/liftover, matrix
ablation, and calibration diagnostics. v9/v10 then tested whether public
candidate discovery can provide old target-tissue support. It cannot, within the
official GEO processed/seed space exercised here.

The correct conclusion is not to run more autoresearch on the same matrices.
The next phase must change data strategy:

- broader acquisition outside the current public GEO processed candidate set;
- a deliberate minimal FASTQ/Bismark ETL pilot for manually justified old target
  tissue datasets;
- collaboration or generation of new old brain_cortex/heart/lung mouse RRBS/WGBS
  data;
- or a benchmark redefinition that stops treating GSE121141 old104+
  brain_cortex/heart/lung as the headline bottleneck.

## Next Recommendation

Start a v10.1 or v11 data-acquisition RFC instead of model tuning. The RFC
should decide whether to:

1. perform manual non-GEO and literature-backed candidate acquisition with
   official accession verification;
2. run a minimal FASTQ/Bismark pilot on a small number of old target-tissue raw
   samples if a justified accession is found;
3. purchase or generate old target-tissue methylome data;
4. redefine the headline benchmark to a tissue/data regime that has sufficient
   public support.

Until one of these is approved, constrained autoresearch, embedding search, and
deep learning are not scientifically meaningful for the current target.


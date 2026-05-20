# v21 Raw-to-Interpretable Mouse Methylation Clock Plan

Date: 2026-05-21

## Scope

v21 moves the project from a processed-matrix demo baseline toward a traceable
raw FASTQ to Bismark COV to 5kb region matrix workflow, then guarded local ML,
AutoDL-only DL, and feature/region interpretation.

This plan does not authorize full 6.5T FASTQ processing. The first execution
unit is a small raw pilot and matrix gate.

## Implementation

New files:

- `scripts/validate/run_v21_raid_raw_ralph_loop.py`
- `scripts/etl/23_run_raw_bismark_queue.py`
- `scripts/etl/24_raw_etl_watchdog.py`
- `scripts/train/export_v21_autodl_package.py`
- `scripts/train/import_v21_autodl_results.py`
- `scripts/validate/run_v21_feature_interpretation.py`
- `doc/30_protocols/09_20260521_raw_to_interpretable_clock_sop.md`

Primary output directory:

- `results/ralph_v21_raid_raw_clock/`

Raw ETL working directory:

- `/data/mouse_methyl/processed_v21_raw_etl/`

## Dataset Priority

1. `GSE121141`: core stress-test, first 2-3 sample raw pilot.
2. `GSE80672`: CR validation pilot after first gate passes.
3. `GSE93957` / `GSE60012`: cross-dataset support.
4. `GSE120137`: raw-vs-processed consistency audit.

## Success Criteria

Raw ETL gate success:

- at least one priority dataset produces Bismark COV and a 5kb region matrix;
- metadata overlap and age coverage are both `>=95%`;
- beta range is valid;
- sex/MT regions are excluded;
- common 5kb regions are `>=50000`;
- raw-derived and processed matrices are audited when a processed baseline
  exists for that dataset.

Local ML demo success:

- GroupKFold MAE `<=24.5w`;
- support-covered MAE `<=20.5w`;
- random-label sanity has `abs(r)<0.2` and MAE at random level;
- GSE80672 CR uses real held-out predictions only.

AutoDL DL demo success:

- `res_mlp` centered search reproduces or exceeds v20;
- GroupKFold MAE `<=20.5w`;
- LODO mean MAE `<=18.5w`;
- no local GPU/DL run is used.

Interpretation success:

- at least 100 top regions have stable importance, age direction,
  confounding labels, and traceable region IDs;
- the report separates cross-tissue stable, tissue-specific,
  dataset-confounded, and coverage-driven signals.

## Current Guarded Commands

```bash
uv run python scripts/validate/run_v21_raid_raw_ralph_loop.py
uv run python scripts/etl/23_run_raw_bismark_queue.py --queue results/ralph_v21_raid_raw_clock/raw_etl_queue.csv
uv run python scripts/etl/24_raw_etl_watchdog.py --queue results/ralph_v21_raid_raw_clock/raw_etl_queue.csv
uv run python scripts/train/export_v21_autodl_package.py --matrix results/multidataset_v8_3_ablation/all6/all_rrbs_region_matrix_5kb.parquet --metadata metadata/model_sample_metadata_v8.csv
uv run python scripts/validate/run_v21_feature_interpretation.py
```

The ETL and watchdog commands are blocked by default unless
`--authorize-bismark` is provided.

## Decision Boundary

If three pilots fail by common-region overlap, metadata, assembly, schema, beta
range, or age coverage, v21 raw ETL expansion stops. The next action is
re-evaluation of Bismark configuration, GEO/SRA/ENA documentation, assay type,
assembly, and local metadata rather than training or wider model search.

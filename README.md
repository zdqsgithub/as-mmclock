# as-mmclock

Mouse methylation clock and biological-signal discovery workflows.

This repository tracks the reusable code, protocols, metadata tables, and lightweight reports for the local mouse methylation project. Large local assets are intentionally excluded from git, including FASTQ files, Bismark outputs, reference indexes, parquet matrices, virtual environments, and downloaded GEO/SRA archives.

## Current Focus

- Build traceable FASTQ-to-region-matrix ETL for local RAID mouse methylation resources.
- Train leakage-controlled ML baselines for mouse methylation age and auxiliary biological signals.
- Keep deep learning runs packaged for AutoDL rather than local GPU execution.
- Interpret selected 5kb methylation regions with explicit dataset, tissue, and coverage confounding labels.

## Important Local Paths

- Raw FASTQ lives outside this repo under `/data/mouse_methyl/raw`.
- v21 raw ETL pilot outputs live under `/data/mouse_methyl/processed_v21_raw_etl`.
- Model-ready processed matrices are generated locally under `results/` and are ignored by git.
- Lightweight v22 screening outputs are tracked under `results/v22_biological_signal_screen/`.

## Reproduce v22 Screen

```bash
uv run python scripts/validate/run_v22_biological_signal_screen.py \
  --include-random-labels \
  --n-feature-prefilter 500 \
  --out-dir results/v22_biological_signal_screen
```

The command expects the local processed 5kb matrix at:

```text
results/multidataset_v8_3_ablation/all6/all_rrbs_region_matrix_5kb.parquet
```

That matrix is intentionally not committed because it is a generated data artifact.

## Repository Boundary

GitHub is used for source control and lightweight scientific reporting. Full raw data, derived matrices, Bismark indexes, BAM/COV files, AutoDL packages, and broad benchmark outputs remain local or should be stored in an artifact/data system instead of git.

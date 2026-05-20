# v21 Raw-to-Interpretable Clock SOP

Date: 2026-05-21

## Purpose

This SOP controls the v21 path from local RAID FASTQ to a traceable 5kb
methylation region matrix, guarded ML benchmarks, AutoDL-only DL runs, and
feature/region interpretation.

## State Machine

1. Raw inventory: scan `/data/mouse_methyl/raw` read-only and write FASTQ,
   metadata, run, and pairing manifests.
2. Environment gate: verify `uv`, project imports, Bismark, Bowtie2, Samtools,
   GRCm38/mm10 Bismark index files, and disk space.
3. Raw pilot: run only a 2-3 sample pilot for `GSE121141` unless a decision
   state explicitly selects another priority dataset.
4. Matrix gate: build Bismark COV and 5kb region matrix, then check metadata
   overlap, exact age coverage, beta range, sex/MT exclusion, and common 5kb
   regions.
5. Learn readiness: if the gate passes, emit fixed local ML commands with
   `training_authorized=false` until explicit approval is recorded.
6. AutoDL package: export matrix, metadata, scripts, and a narrow `res_mlp`
   command manifest for GPU execution outside the local workstation.
7. Interpretation: produce top-region, cluster, and confounding audit tables
   before any v21 delivery report.

## Default Commands

Inventory and gate preflight. This command does not run Bismark, download data,
train, or start autoresearch:

```bash
uv run python scripts/validate/run_v21_raid_raw_ralph_loop.py
```

Build the guarded raw ETL queue from the latest inventory:

```bash
uv run python scripts/validate/run_v21_raid_raw_ralph_loop.py \
  --pilot-samples 3 \
  --priority-datasets GSE121141,GSE80672,GSE93957,GSE60012,GSE120137
```

Run a raw Bismark queue only after explicit approval. Keep outputs under
`/data/mouse_methyl/processed_v21_raw_etl/`:

```bash
uv run python scripts/etl/23_run_raw_bismark_queue.py \
  --queue results/ralph_v21_raid_raw_clock/raw_etl_queue.csv \
  --out-root /data/mouse_methyl/processed_v21_raw_etl \
  --parallel-samples 3 \
  --bismark-threads 8 \
  --authorize-bismark
```

Use the watchdog for background/resume execution after approval:

```bash
uv run python scripts/etl/24_raw_etl_watchdog.py \
  --queue results/ralph_v21_raid_raw_clock/raw_etl_queue.csv \
  --out-root /data/mouse_methyl/processed_v21_raw_etl \
  --start-runner \
  --authorize-bismark
```

After a raw matrix gate passes and training approval is recorded, run local ML
only with `uv run` and CPU-friendly models. The v21 controller writes command
manifests; do not hand-edit them into a broader search space.

Package DL work for AutoDL only:

```bash
uv run python scripts/train/export_v21_autodl_package.py \
  --matrix <raw_or_approved_combined_matrix.parquet> \
  --metadata metadata/model_sample_metadata_v8.csv
```

Import AutoDL results after the remote run:

```bash
uv run python scripts/train/import_v21_autodl_results.py \
  --source <autodl_results_dir_or_tgz>
```

Run feature interpretation:

```bash
uv run python scripts/validate/run_v21_feature_interpretation.py
```

## Gates

Raw matrix promotion requires:

- metadata overlap `>=95%`;
- age coverage `>=95%`;
- beta values within `[0,1]`;
- sex chromosome and mitochondrial region exclusion;
- common autosomal 5kb regions `>=50000`;
- raw-vs-processed consistency audit when the same dataset has a processed
  baseline.

The only successful raw state is
`raw_matrix_ready_for_learn_pending_training_approval`. Training and
autoresearch remain unauthorized until a separate approval is recorded.

## Failure Handling

Apply a 3-strike rule. If three raw pilots fail or are blocked by common-region
overlap, metadata, assembly, schema, beta range, or age coverage, stop raw ETL
expansion and write a re-evaluation state. Do not continue training,
autoresearch, or DL packaging from the blocked candidate set.

## Guardrails

- RAID FASTQ is read-only input.
- No downloads start from this workflow by default.
- No local GPU/DL run is authorized.
- No dummy CR or intervention AUC is allowed.
- All preprocessing, imputation, feature selection, and model fitting must be
  fit inside the train fold or train dataset only.
- Feature interpretation must distinguish cross-tissue stable,
  tissue-specific, dataset-confounded, and coverage-driven signals.

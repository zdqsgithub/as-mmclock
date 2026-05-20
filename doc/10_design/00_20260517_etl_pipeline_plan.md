# Mouse Methylation ETL Pipeline — Local Machine Processing Plan

> **Date:** 2026-05-17  
> **Scope:** FASTQ → Beta Matrix (all processing on local machine, no GPU required)  
> **Input:** `/home/zdq-as/as-ds-ops/data/mouse_methyl/raw/` (6.5 TB, 8 GEO datasets)  
> **Output:** `/home/zdq-as/mouse_methyl_work/` (beta matrices, metadata, QC reports)  
> **Conforms to:** AS-DS-Ops SKILL.md v2.5.0 methylation extension rules

---

## Phase 0: Rapid Prototyping with Thompson Pre-Computed Matrix

**Timeline:** Day 1–2 (immediate)  
**Goal:** Get a working ElasticNet baseline running in < 2 days without any alignment

### Available Data
- `GSE120137_MM10_EpigeneticAgeClock_RidgeRegression.csv.gz` (2 MB, already on disk)
  - Contains: pre-processed beta values for ~193k CpGs × ~548 samples
  - Includes: Ridge regression coefficients from Thompson 2018

### Steps

```python
# scripts/etl/00_load_thompson_matrix.py
# 1. Decompress and load the pre-built beta matrix
# 2. Parse sample metadata from GEO (age_days, tissue, strain, sex, intervention)
# 3. Apply coverage filter: keep CpGs with < 20% missing values
# 4. Build metadata CSV: sample_id, age_days, tissue, strain, sex, dataset_batch
# 5. Export: processed_beta_matrix.parquet + sample_metadata.csv
```

**Output:** `results/phase0/` — baseline model results

---

## Phase 1: Metadata Integration & Standardization

**Timeline:** Day 2–5  
**Goal:** Create a unified, publication-ready sample metadata table across all 8 datasets

### Required Fields (per SKILL.md rules)

| Field | Standardization Rule | Source |
|-------|---------------------|--------|
| `sample_id` | SRR accession | SRA metadata |
| `age_days` | weeks×7, months×30.42 | GEO sample attributes |
| `tissue` | Ontology: blood/liver/lung/heart/brain_cortex/brain_hippo/brain_cerebellum/muscle/spleen/fibroblast | GEO |
| `strain` | C57BL/6 / BALB/c / DBA/2 / CAST/EiJ / hybrid | GEO |
| `sex` | M / F / unknown | GEO |
| `intervention` | control / CR / rapamycin / dwarfism / castration / diet | GEO |
| `dataset_batch` | GEO series ID (e.g., GSE80672) | derived |
| `assay` | RRBS (all; WGBS excluded from training) | GEO |
| `paired_end` | True/False | SRA metadata |

### Scripts

```bash
# scripts/etl/01_build_metadata.py
# Downloads GEO Series Matrix files, parses sample attributes
# Cross-references with SRA Run Selector JSON from existing metadata dir
# Outputs: metadata/unified_sample_metadata.csv
```

### QC Checks
- [ ] No sample with missing `age_days` enters training
- [ ] No sample with `assay == WGBS` enters training set (external validation only)
- [ ] Age distribution plot by dataset (flag imbalance > 18 months)

---

## Phase 2: FASTQ Quality Control & Trimming

**Timeline:** Day 5–14 (parallel with metadata, batch processing)  
**Tool:** Trim Galore ≥ 0.6.0  
**Constraint:** RRBS-specific trimming (`--rrbs` flag); DO NOT use generic adapter trimming

### Trim Galore Command (RRBS PE)
```bash
trim_galore \
  --rrbs \
  --paired \
  --quality 20 \
  --length 20 \
  --cores 8 \
  --fastqc \
  --output_dir /home/zdq-as/mouse_methyl_work/trimmed/{GSE_ID}/{SRR}/ \
  {SRR}_1.fastq {SRR}_2.fastq
```

### Script
```bash
# scripts/etl/02_trim_galore_batch.py
# Reads unified_sample_metadata.csv
# Generates SLURM-style job array (or local parallel via subprocess)
# 8 workers × local machine
# Logs: results/qc/trimming_summary.tsv
```

### Hard QC Filters (Red Lines)
- FastQC GC content deviation > ±5% → flag for manual review
- Adapter content > 30% after trimming → discard (re-check if single-end was run as PE)
- Very short reads (< 20 bp after trimming) > 20% → discard

---

## Phase 3: Bismark Alignment

**Timeline:** Day 10–30 (resource-intensive; GPU NOT needed; multi-core CPU)  
**Tool:** Bismark ≥ 0.24.0  
**Reference:** GRCm38 / mm10 (must be indexed once locally)

### One-Time: Build Bismark Index
```bash
# scripts/etl/03a_build_bismark_index.sh
# Reference: GRCm38.primary_assembly.fa (~2.7 GB download from Ensembl)
# Index location: /home/zdq-as/mouse_methyl_work/reference/mm10_bismark_index/
bismark-genome-preparation \
  --parallel 8 \
  /home/zdq-as/mouse_methyl_work/reference/mm10/
```

### Per-Sample Alignment
```bash
bismark \
  --genome /home/zdq-as/mouse_methyl_work/reference/mm10/ \
  --parallel 4 \
  --non_directional \          # for RRBS: MspI cuts are non-directional
  -1 {SRR}_1_val_1.fq \
  -2 {SRR}_2_val_2.fq \
  --output_dir /home/zdq-as/mouse_methyl_work/bismark_bam/{GSE}/{SRR}/
```

### ⚠️ CRITICAL: NO DEDUPLICATION FOR RRBS
```bash
# DO NOT RUN: bismark_deduplicate (destroys signal for RRBS)
# MspI cut sites produce identical fragment starts by design
```

### QC Hard Filters (SKILL.md Red Lines)
- Bisulfite conversion rate < 95% → discard sample (flag in metadata)
- Mapping rate < 50% → investigate; may indicate wrong reference or data corruption
- Mean CpG coverage < 5× → discard

---

## Phase 4: Methylation Extraction

**Timeline:** Day 25–40  
**Tool:** Bismark Methylation Extractor

```bash
bismark_methylation_extractor \
  --paired-end \
  --comprehensive \
  --CpG_context \
  --cytosine_report \
  --genome_folder /home/zdq-as/mouse_methyl_work/reference/mm10/ \
  --output /home/zdq-as/mouse_methyl_work/cov_files/{GSE}/{SRR}/ \
  {SRR}_bismark_bt2_pe.bam
```

**Output per sample:** `{SRR}.bismark.cov.gz`  
Format: `chr  start  end  methylation%  count_M  count_U`

---

## Phase 5: Beta Matrix Construction

**Timeline:** Day 40–50  
**Script:** `scripts/etl/05_build_beta_matrix.py`

### Algorithm

```python
# 1. Load all .cov.gz files
# 2. For each CpG: compute beta = count_M / (count_M + count_U)
# 3. Filter: CpGs with coverage < 5x → set to NaN
# 4. Filter: CpGs present in < 50% of samples within a dataset → drop
# 5. Union CpGs across datasets: only keep sites present in >= 80% of ALL samples
# 6. Exclude: sex chromosomes (chrX, chrY), mitochondrial (chrM)
# 7. Output:
#    - beta_matrix_full.parquet (all CpGs)
#    - beta_matrix_union_filtered.parquet (model-ready)
#    - cpg_coverage_stats.csv (per-CpG coverage summary)
```

### Storage Format
- **Parquet** (not CSV) for efficient memory-mapped access during training
- Rows: CpG sites (`chr_pos` string key e.g. `chr1_3000827`)
- Columns: SRR sample IDs
- Values: float32 beta [0.0, 1.0], NaN for low-coverage

---

## Phase 6: Feature Engineering Options

### Option A: Site-Based (ElasticNet baseline)
- Direct beta values per CpG
- ~18k–50k sites after union filtering
- Input to ElasticNet/Ridge/MLP

### Option B: Region-Based (Meer 2023 style)
```python
# Tile genome into 5 kb non-overlapping windows
# Compute mean beta per window per sample
# ~50k regions (vs. ~millions of raw CpGs)
# More robust to coverage variation
```

### Option C: Epiallele / Heterogeneity Features (novel DL input)
```python
# From BAM files: compute per-read methylation strings
# WSH (within-sample heterogeneity) metrics:
#   - PDR (proportion of discordant reads)
#   - Epipolymorphism
#   - Shannon entropy of methylation patterns
# These are age-correlated beyond mean beta (unexplored for clocks)
```

---

## Phase 7: QC Report

**Script:** `scripts/etl/07_qc_report.py`

**Outputs:**
- `results/qc/sample_qc_summary.csv` — per-sample pass/fail
- `results/qc/conversion_rate_distribution.png`
- `results/qc/coverage_distribution.png`
- `results/qc/age_distribution_by_dataset.png`
- `results/qc/pca_batch_check.png` — detect batch effects
- `results/qc/tissue_heatmap.png` — sample distribution

---

## Resource Requirements (Local Machine)

| Step | CPU | Memory | Disk I/O | Estimated Time |
|------|-----|--------|----------|----------------|
| Phase 0 (prototype) | 4 cores | 16 GB | Low | 2–4 hours |
| Phase 1 (metadata) | 2 cores | 4 GB | Low | 1–2 days |
| Phase 2 (trim) | 32 cores | 16 GB | High | 3–7 days |
| Phase 3 (Bismark) | 32 cores | 64 GB | Very High | 14–21 days |
| Phase 4 (extract) | 16 cores | 32 GB | High | 5–10 days |
| Phase 5 (beta matrix) | 8 cores | 64 GB | Medium | 2–3 days |

**Critical:** Bismark alignment requires ~10 GB RAM per 4-thread job. Run 8 parallel jobs = 80 GB RAM. Check available RAM before starting.

---

## Data Flow Diagram

```
/as-ds-ops/data/mouse_methyl/raw/{GSE}/fastq/*.fastq
          ↓ Phase 2: Trim Galore
/mouse_methyl_work/trimmed/{GSE}/{SRR}/*_val_*.fq
          ↓ Phase 3: Bismark Align
/mouse_methyl_work/bismark_bam/{GSE}/{SRR}/*.bam
          ↓ Phase 4: Methylation Extractor
/mouse_methyl_work/cov_files/{GSE}/{SRR}/*.bismark.cov.gz
          ↓ Phase 5: Beta Matrix Builder
/mouse_methyl_work/results/beta_matrices/
  ├── beta_matrix_full.parquet
  ├── beta_matrix_union_filtered.parquet
  └── sample_metadata_unified.csv
          ↓ Upload to Cloud GPU
          ↓ Phase: Model Training
```

---

## Deliverables for Cloud Training Hand-off

The following files must be ready before cloud GPU training begins:

- [ ] `metadata/unified_sample_metadata.csv` — all samples with age, tissue, strain, sex, intervention, batch
- [ ] `results/beta_matrices/beta_matrix_union_filtered.parquet` — coverage-filtered beta values
- [ ] `results/qc/sample_qc_summary.csv` — QC pass/fail per sample
- [ ] `configs/train_test_split.yaml` — predefined train/val/test assignments
- [ ] `results/qc/pca_batch_check.png` — batch effect assessment

---

## References

1. Krueger & Andrews (2011). Bismark. Bioinformatics.
2. Martin (2011). Trim Galore / Cutadapt. EMBnet Journal.
3. Meer et al. (2023). Region-based epigenetic clock design. Aging Cell.
4. Thompson et al. (2018). GSE120137 processed matrix.

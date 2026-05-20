# Mouse Methylation Biological Age Prediction — Project Overview

> **Project Root:** `/home/zdq-as/mouse_methyl_work`  
> **Spec Root:** `/home/zdq-as/as-ds-ops`  
> **Created:** 2026-05-17  
> **Updated:** 2026-05-17 (Phase 0 Autoresearch Execution)
> **Status:** 🟢 Phase 0 Executing (Autoresearch 100+ Models)  
> **Conforms to:** AS-DS-Ops SKILL.md v3.0.0 + Methyl-Clock Extension v1.0

---

## Mission

Build a **production-grade mouse epigenetic clock** using RRBS bisulfite sequencing data that:

1. **Predicts chronological age** from raw methylation data with MAE < 3.5 weeks and Pearson r > 0.90
2. **Measures biological age acceleration** validated against calorie restriction (CR) and rapamycin interventions
3. **Benchmarks multiple architectures** from classical linear models to deep learning (CNN/Transformer)
4. **Uses autoresearch loops** (karpathy-inspired) for autonomous hyperparameter optimization
5. Outputs **all standard epigenetic clock benchmark metrics** required for publication

---

## Data Summary (Current Status: 2026-05-17)

| Dataset | FASTQs on disk | Size | Status |
|---------|----------------|------|--------|
| GSE120137 (Thompson 2018) | 573 files | 2.6 TB | ✅ Mostly complete |
| GSE80672 (Petkovich 2017) | 532 files | 1.8 TB | ✅ Mostly complete |
| GSE60012 (Reizel 2015) | 210 files | 736 GB | ✅ Mostly complete |
| GSE93957 (Stubbs 2017) | 118 files | 843 GB | ✅ Mostly complete |
| GSE52266 (Cannon 2013) | 128 files | 182 GB | ✅ Complete |
| GSE121141 (Meer 2018) | 24 files | 176 GB | ⚠️ Partial |
| GSE45361 (Schillebeeckx) | 22 files | 148 GB | ✅ Complete |
| GSE80761 (Zhang) | 4 files | 47 GB | ✅ Complete |

**Total raw data on disk:** ~6.5 TB (in `/home/zdq-as/as-ds-ops/data/mouse_methyl/raw/`)  

**Processing Milestones:**
- **Phase 1 Metadata:** ✅ All 8 datasets metadata successfully extracted via `GEOparse` with clinical age (days/weeks) and interventions dynamically parsed. (`metadata/unified_sample_metadata.csv`)
- **Phase 0 Matrix:** ✅ Prototype Bismark coverage files from GSE120137 parsed into high-performance `beta_matrix_thompson.parquet` (921k CpGs × 549 samples).
- **Phase 0 Autoresearch:** 🔄 Currently orchestrating 100+ deep learning experiments on the prototype matrix.

---

## Architecture Philosophy

```
Raw FASTQ → ETL Pipeline → Beta Matrix → Model Training (Cloud GPU) → Benchmark
    ↑                                           ↑
Local machine                           Cloud (GPU rental)
```

### Local Machine (ETL)
- Trim Galore → Bismark alignment → Methylation extraction → Coverage filtering
- Beta matrix construction, QC, metadata integration
- Feature engineering (site-based, region-based, epiallele-level)

### Cloud GPU (Training)
- ElasticNet / Ridge (sklearn, quick baseline)
- CNN-based models (PyTorch)
- Transformer / MethylBERT-style (PyTorch)
- Autoresearch optimization loops

---

## Primary Benchmark Targets

| Metric | Category | Target | Reference |
|--------|----------|--------|-----------|
| Pearson r | Regression | > 0.90 | Stubbs 2017, Thompson 2018 |
| MAE (weeks) | Regression | < 3.5 weeks | Meer 2018 |
| MedAE (weeks) | Regression | < 3.0 weeks | Thompson 2018 |
| R² | Regression | > 0.80 | standard |
| AUC-ROC (CR vs. control) | Classification | > 0.80 | Petkovich 2017 |
| F1 (CR detection) | Classification | > 0.75 | Petkovich 2017 |
| Age acceleration (CR) | Biological validation | Significant negative | Meer 2023 |
| Age acceleration (rapamycin) | Biological validation | Significant negative | GSE67507 |
| Cross-dataset MAE | Generalization | < 5.0 weeks | Meer 2023 |

---

## Directory Structure

```
mouse_methyl_work/
├── doc/                        # All project documents (AS-DS-Ops convention)
│   ├── 00_meta/                # Project specs, overviews
│   ├── 10_design/              # Architecture plans, RFCs
│   ├── 20_analysis/            # Analysis reports, data surveys
│   ├── 30_protocols/           # ETL SOPs, training protocols
│   ├── 40_meeting/             # Decision logs
│   └── 99_archive/             # Superseded docs
├── scripts/                    # All executable code
│   ├── etl/                    # Data processing (local machine)
│   ├── train/                  # Model training (cloud GPU)
│   ├── validate/               # Evaluation and benchmarking
│   └── utils/                  # Shared utilities
├── configs/                    # YAML/JSON configuration files
├── notebooks/                  # Exploratory analysis Jupyter notebooks
├── metadata/                   # Sample metadata CSVs
├── models/                     # Saved model checkpoints (small files only)
├── results/                    # Benchmark results and autoresearch logs
│   ├── benchmark/
│   └── autoresearch_log/
└── raw_downloads/              # (currently empty; data in as-ds-ops/data/)
```

---

## Key References

1. Stubbs et al. (2017). *Multi-tissue DNA methylation age predictor in mouse.* Genome Biol. GEO: GSE93957.
2. Petkovich et al. (2017). *Using DNA methylation profiling to evaluate biological age.* Cell Metab. GEO: GSE80672.
3. Thompson et al. (2018). *A multi-tissue full developmental time course epigenetic mouse clock.* Aging. GEO: GSE120137.
4. Meer et al. (2018). *A whole lifespan mouse multi-tissue DNA methylation clock.* eLife. GEO: GSE121141.
5. Meer et al. (2023). *Region-based epigenetic clock design improves RRBS-based age prediction.* Aging Cell.
6. karpathy/autoresearch — Autonomous research loop framework.

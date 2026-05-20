# Mouse Epigenetic Clock — Model Architecture & Optimization Plan

> **Date:** 2026-05-17  
> **Scope:** Model training strategies, cloud GPU architecture, autoresearch loop  
> **Target:** Best-in-class mouse age prediction with biological validation  
> **Conforms to:** AS-DS-Ops SKILL.md v2.5.0 + karpathy/autoresearch pattern

---

## Overview

This document describes the **model training and optimization strategy** for the mouse epigenetic clock project. All training runs on **cloud GPU** (no local GPU required). The autoresearch loop automatically explores the hypothesis space and logs results.

---

## Tier 0: Baseline Models (Required Before DL)

> All baselines MUST be implemented and benchmarked before any deep learning

### B1: Thompson Ridge Coefficients (Published Baseline)
- Use the pre-computed `GSE120137_MM10_EpigeneticAgeClock_RidgeRegression.csv.gz`
- Apply published Ridge coefficients directly to held-out test set
- **Purpose:** Sanity check — our pipeline must match or beat this
- **Expected:** Pearson r ≈ 0.95, MAE ≈ 3.3 weeks (on training data)

### B2: ElasticNet (L1+L2)
```python
from sklearn.linear_model import ElasticNetCV
model = ElasticNetCV(
    l1_ratio=[0.1, 0.5, 0.7, 0.9, 0.95, 0.99, 1.0],
    cv=GroupKFold(n_splits=5),
    max_iter=10000,
    n_jobs=-1
)
```
- **Features:** Union-filtered beta matrix (~18k–50k CpGs)
- **Target:** log(age_days + 1)
- **CV:** GroupKFold by dataset_batch (prevents dataset leakage)
- **Expected N_sites selected:** 300–500 CpGs

### B3: Ridge (All-CpG)
```python
from sklearn.linear_model import RidgeCV
model = RidgeCV(alphas=np.logspace(-3, 3, 100), cv=5)
```
- **Advantage:** More robust for mixed-strain cohorts
- **Advantage:** Better detects subtle intervention effects (CR, rapamycin)

### B4: Random Forest (Ensemble Baseline)
```python
from sklearn.ensemble import RandomForestRegressor
model = RandomForestRegressor(n_estimators=500, max_depth=6, n_jobs=-1)
```
- Use on **region-averaged beta** (not full CpG matrix — too wide for RF)
- **Purpose:** Capture non-linear interactions; validate if DL is needed

---

## Tier 1: Regional Clock (Meer 2023 Style)

> Key insight: Regional averaging (5 kb windows) improves cross-dataset transfer

```python
# Region definition: non-overlapping 5 kb bins across mm10
# ~500k bins total; after filtering: ~50k informative regions

# For each sample, for each region:
# region_beta = mean(site_betas within region, ignore NaN)

# Then apply ElasticNet/Ridge on region_beta matrix
# Expected: fewer sites needed (200–300 regions), better transferability
```

**Benchmark target:**
- Cross-dataset MAE: < 5 weeks (vs. site-based MAE which can be 8+ weeks)
- CR detection AUC: > 0.80

---

## Tier 2: Deep Learning — 1D-CNN

> Minimal DL; validates whether nonlinearity helps beyond linear models

### Architecture
```python
class MethylCNN(nn.Module):
    def __init__(self, n_cpgs):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(n_cpgs, 1024), nn.ReLU(), nn.Dropout(0.3),
            nn.Linear(1024, 256), nn.ReLU(), nn.Dropout(0.2),
            nn.Linear(256, 64), nn.ReLU(),
            nn.Linear(64, 1)  # log(age_days + 1)
        )
    
    def forward(self, x):
        return self.net(x)
```

**Training config:**
```yaml
optimizer: AdamW
lr: 1e-3
weight_decay: 1e-4
batch_size: 32
max_epochs: 200
early_stopping: patience=20, monitor=val_mae
scheduler: CosineAnnealingLR
loss: HuberLoss(delta=1.0)  # robust to age outliers
```

---

## Tier 3: Deep Learning — Sequence-Aware CNN (Novel Innovation)

> Combines CpG beta values with genomic sequence context

### Input Representation
```
For each CpG site i:
  - beta_i: scalar [0.0, 1.0]
  - sequence_context: 101 bp window centered on CpG → one-hot encoded → (101, 5)
    [A, T, C, G, CpG_flag]

Combined per-CpG input: (101, 5) + scalar beta
```

### Architecture (Sequence-Aware Multi-Scale CNN)
```python
class SequenceAwareMethylClock(nn.Module):
    """
    Innovation: treats each CpG as a sequence context + methylation pair.
    Instead of a flat vector, learns which sequence motifs are age-informative.
    """
    def __init__(self, n_cpgs=5000, seq_len=101):
        super().__init__()
        # Per-CpG sequence encoder (shared weights across CpGs)
        self.seq_encoder = nn.Sequential(
            nn.Conv1d(5, 32, kernel_size=9, padding=4), nn.ReLU(),
            nn.Conv1d(32, 64, kernel_size=5, padding=2), nn.ReLU(),
            nn.AdaptiveAvgPool1d(8),
            nn.Flatten()  # → 512-dim per CpG
        )
        # Methylation context fusion
        self.fuse = nn.Linear(512 + 1, 64)  # 512 seq + 1 beta
        # Across-CpG aggregation (attention)
        self.attn = nn.MultiheadAttention(64, num_heads=4, batch_first=True)
        # Regression head
        self.head = nn.Sequential(
            nn.Linear(64, 32), nn.ReLU(),
            nn.Linear(32, 1)
        )
```

**This is the main innovation** — no published mouse clock has combined sequence context with methylation signal end-to-end.

---

## Tier 4: Transformer Foundation Model Fine-tuning

> Fine-tune CpGPT or MethylGPT on our dataset (if pre-trained weights available)

- **CpGPT** (2024): Pre-trained on 100M+ human CpG sites; requires species transfer evaluation
- **MethylGPT** (2024): Broader epigenomic tasks; check availability of mouse pre-trained weights
- **Strategy:** Fine-tune last 2 layers only on our 1000+ mouse samples

---

## Cross-Validation Protocol (NON-NEGOTIABLE)

```python
# GroupKFold stratified by BOTH dataset_batch AND tissue
# Never let samples from the same dataset appear in both train and test
from sklearn.model_selection import GroupKFold

# Primary CV: 5-fold with dataset-level groups
cv = GroupKFold(n_splits=5)
groups = sample_metadata['dataset_batch']  # GSE80672, GSE120137, etc.

# For final evaluation only:
# Train: GSE120137 + GSE121141 + GSE60012 (N ≈ 700–800)
# Val:   GSE93957 (N = 62, held out during hyperparam tuning)
# Test:  GSE80672 (N ≈ 255, includes CR/rapamycin — biological validation)
```

---

## Benchmark Metrics Implementation

```python
# scripts/validate/benchmark_metrics.py

import numpy as np
from scipy import stats
from sklearn.metrics import roc_auc_score, f1_score

def compute_regression_metrics(y_true, y_pred, log_scale=True):
    """Standard epigenetic clock regression metrics."""
    if log_scale:
        # Convert from log(days+1) back to days
        y_true = np.expm1(y_true) / 7  # → weeks
        y_pred = np.expm1(y_pred) / 7
    
    r, pval = stats.pearsonr(y_true, y_pred)
    mae = np.mean(np.abs(y_true - y_pred))
    medae = np.median(np.abs(y_true - y_pred))
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - y_true.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot
    
    return {'pearson_r': r, 'pearson_pval': pval, 
            'mae_weeks': mae, 'medae_weeks': medae, 'r2': r2}

def compute_biological_age_acceleration(y_pred, y_true, group_labels):
    """Age acceleration per dataset to avoid batch confounding."""
    from sklearn.linear_model import LinearRegression
    accelerations = {}
    for group in np.unique(group_labels):
        mask = group_labels == group
        lr = LinearRegression()
        lr.fit(y_pred[mask].reshape(-1, 1), y_true[mask])
        residuals = y_true[mask] - lr.predict(y_pred[mask].reshape(-1, 1))
        accelerations[group] = residuals
    return accelerations

def compute_intervention_classification(accelerations, intervention_labels):
    """
    AUC and F1 for binary CR detection.
    'CR' samples should have significantly negative age acceleration.
    """
    # acceleration is predictor; intervention (CR vs control) is label
    y_score = -accelerations  # negative acceleration → higher CR score
    y_true_binary = (intervention_labels == 'CR').astype(int)
    
    auc = roc_auc_score(y_true_binary, y_score)
    # Threshold at 0 (negative acceleration = CR-like)
    y_pred_binary = (accelerations < 0).astype(int)
    f1 = f1_score(y_true_binary, y_pred_binary)
    
    return {'cr_detection_auc': auc, 'cr_detection_f1': f1}
```

---

## Autoresearch Loop (Implemented & Executing)

### Hypothesis Space (Updated 2026-05-17)

| Parameter | Range | Default | Notes |
|-----------|-------|---------|-------|
| `model_type` | elasticnet / ridge / mlp / lgbm / rf | elasticnet | Active in Phase 0 |
| `feature_type` | age_correlation | age_correlation | Replaced 'variance' due to tissue bias |
| `n_cpg_prefilter` | 5k / 10k / 20k / 50k / 100k | 20k | Pearson correlation with age |
| `imputation` | median / mean | median | Handled inside CV pipeline |
| `batch_size` (NN) | 16 / 32 / 64 | 32 | |

> **Critical Discovery (Tissue Variance Trap):** 
> Initial baselines using maximum variance (`np.nanvar`) for feature selection failed catastrophically (R² < 0) because the most variable CpGs in a multi-tissue dataset are driven by cell-type differences, not chronological age. We successfully pivoted to **Age-Correlation Filtering** (computing Pearson correlation with age across all samples), which immediately yielded $R > 0.78$ even with simple linear models.

### Loop Configuration
```yaml
# configs/autoresearch_config.yaml
metric_name: pearson_r
maximize: true
budget_seconds: 1800  # 30 min per experiment on cloud GPU
n_iterations: 30
results_tsv: results/autoresearch_log/results.tsv
git_dir: /home/zdq-as/mouse_methyl_work
```

### Loop Protocol
```
LOOP for n_iterations:
  1. git commit current state (save checkpoint)
  2. Sample next hypothesis from param_space (hill-climb from best + random mutation)
  3. Run training script: python scripts/train/train_clock.py --config {hypothesis_config}
  4. Read result from results/last_run.json (pearson_r, mae_weeks, cr_auc)
  5. If improved → KEEP commit, update best
  6. If not improved → git reset --hard (discard changes)
  7. Log to results/autoresearch_log/results.tsv
```

---

## Cloud GPU Training Setup

### Recommended Cloud Instance
- **Provider:** Lambda Labs / Vast.ai / RunPod
- **Instance:** A100 40GB (preferred) or RTX 4090 24GB
- **Cost estimate:** ~$1–2/hr; budget 20 hrs = $20–40 per architecture tier

### Data Transfer Plan
1. Compress beta matrix: `parquet → zstd` (~50–200 MB)
2. rsync to cloud instance scratch
3. Training run → download only `results.json` + model checkpoint
4. Local machine stores all source data permanently

### Environment (cloud)
```bash
# configs/environment_cloud.yml
name: mouse-methyl-clock
channels:
  - conda-forge
  - pytorch
dependencies:
  - python=3.11
  - pytorch>=2.1
  - torchvision
  - scikit-learn>=1.4
  - pandas>=2.0
  - numpy>=1.26
  - scipy>=1.12
  - pyarrow>=14.0
  - optuna>=3.0   # hyperparameter search
  - pip:
    - lightning>=2.2
    - wandb        # optional: experiment tracking
```

---

## Model Registry & Reproducibility

Every training run MUST produce:
```
results/runs/{timestamp}_{model_type}/
├── result.json          # all metrics
├── config.yaml          # exact hyperparameters used
├── model_checkpoint.pt  # or model_coefficients.csv for linear
├── predictions.csv      # y_true, y_pred, sample_id, dataset
└── reproducibility/
    ├── commands.sh
    ├── environment.yml
    └── checksums.sha256
```

---

## Publication-Ready Benchmark Table

Final benchmark report must include ALL of the following:

| Clock | Architecture | N_features | Train MAE (wks) | Val MAE (wks) | Test MAE (wks) | Pearson r (test) | CR AUC | CR F1 |
|-------|-------------|------------|-----------------|---------------|----------------|------------------|--------|-------|
| Thompson 2018 (published) | Ridge | 193k | — | — | ~3.3 | ~0.95 | — | — |
| Ours: ElasticNet | Linear | ~400 | | | | | | |
| Ours: Regional-Ridge | Linear | ~200 regions | | | | | | |
| Ours: MLP | DL | ~20k | | | | | | |
| Ours: SeqCNN | DL+Sequence | ~5k | | | | | | |
| Ours: Best (autoresearch) | Auto | Auto | | | | | | |

---

## Red Flags — Stop Conditions

- Any model reports test Pearson r > 0.99 → **data leakage — audit CV split**
- Age acceleration not significant for CR mice → **clock is not biologically valid**
- Training MAE << Test MAE (gap > 2×) → **overfitting — add regularization**
- NaN loss at epoch 1 → **check age transform and NaN handling in beta matrix**

---

## References

1. Thompson et al. (2018). Aging. GEO: GSE120137. (Ridge baseline)
2. Petkovich et al. (2017). Cell Metab. (CR validation, GSE80672)
3. Meer et al. (2023). Aging Cell. (Regional clock design)
4. karpathy/autoresearch — Autonomous research loop
5. CpGPT (2024). Transformer foundation model for CpG methylation.
6. DeepAge (2024). TCN architecture for biological age. bioRxiv.

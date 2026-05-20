# Literature-Constrained Region Clock v3 Report

> Date: 2026-05-17  
> Scope: GSE80672 age-unit provenance, Phase 0 GSE120137 5 kb region matrix, leakage-safe v3 autoresearch

## Literature Decision

The next step was region-based modeling rather than deep learning or full FASTQ ETL. This follows the project AS-DS-Ops constraints and the RRBS clock literature:

- Simpson/Meer 2023 reports that region-based RRBS clock design improves robustness over single-CpG models in mouse methylation clocks.
- Thompson 2018/GSE120137 remains a Phase 0 within-dataset training source, so its KFold metrics are not a cross-dataset benchmark.
- GSE80672 is the key Petkovich intervention dataset for CR/biological-age validation, but local inventory currently contains only metadata and the GEO SOFT file, not a usable beta/coverage matrix.
- Stubbs 2017 and Bell 2019 support separating chronological age accuracy from biological age/intervention claims.

References:
- Simpson/Meer 2023: https://mayoclinic.elsevierpure.com/en/publications/region-based-epigenetic-clock-design-improves-rrbs-based-age-pred/
- Thompson 2018: https://www.aging-us.com/article/101590/text
- GSE80672 GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80672
- Stubbs 2017: https://pmc.ncbi.nlm.nih.gov/articles/PMC5389178/
- Bell 2019: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1824-y

## Implemented Changes

- Added `configs/metadata_overrides.yaml` and applied a GSE80672 `age_unit_override: months` in `scripts/etl/01_build_metadata.py`.
- Metadata now keeps `raw_age_value`, `raw_age_unit`, `normalized_age_unit`, `age_unit_source`, and `age_parse_warning`.
- Added `scripts/etl/06_build_region_matrix.py` for fixed 5 kb non-overlapping region aggregation from the existing GSE120137 beta matrix.
- Extended `scripts/train/train_clock.py` with `--matrix_path`, `--feature_type site|region`, and `--n_feature_prefilter`, while keeping `--n_cpg_prefilter` compatible.
- Extended `scripts/train/autoresearch_loop.py` to run region v3 into `results/autoresearch_v3_region/` without overwriting v1/v2.

## Validation Results

Metadata:

| Dataset | Check | Result |
|---|---|---:|
| GSE80672 | Samples with age | 255/255 |
| GSE80672 | Age unit source | 255 `dataset_override_literature_geo` |
| GSE80672 | Age weeks range | 2.912-152.100 |
| GSE80672 | Intervention labels | 223 control, 32 CR |

Region matrix:

| Artifact | Result |
|---|---:|
| Matrix | `results/phase0/region_matrix_5kb.parquet` |
| Stats | `results/phase0/region_stats_5kb.csv` |
| Shape | 77,826 regions x 549 samples |
| Filters | min 3 CpGs/region, present in >=440 samples |
| Sex/MT chromosomes | excluded |

Smoke tests:

| Config | Pearson r | MAE weeks | R2 |
|---|---:|---:|---:|
| ridge region 1000 mean | 0.6606 | 17.897 | 0.3657 |
| elasticnet region 2000 median | 0.7769 | 14.050 | 0.5711 |
| lgbm region 1000 mean | 0.8507 | 11.832 | 0.6976 |
| rf region 1000 mean | 0.8166 | 14.253 | 0.6030 |

Random-label sanity check:

| Config | Pearson r | MAE weeks | R2 |
|---|---:|---:|---:|
| randomized lgbm region 1000 mean | -0.0097 | 27.059 | -0.2740 |

The random-label collapse supports that the region training path does not introduce obvious label leakage.

## Autoresearch v3 Results

The v3 search space contained 48 unique configurations. All 48 completed; the remaining 52 requested attempt slots were recorded as `skipped_exhausted`.

Best v3 by composite score and MAE:

| Rank | Experiment | Model | Regions | Imputation | Pearson r | MAE weeks | R2 |
|---:|---|---|---:|---|---:|---:|---:|
| 1 | exp_003_lgbm_500_median | lgbm | 500 | median | 0.8530 | 11.550 | 0.7051 |
| 2 | exp_024_lgbm_500_mean | lgbm | 500 | mean | 0.8524 | 11.616 | 0.7041 |
| 3 | exp_037_lgbm_1000_mean | lgbm | 1000 | mean | 0.8507 | 11.832 | 0.6976 |
| 4 | exp_014_lgbm_1000_median | lgbm | 1000 | median | 0.8512 | 11.908 | 0.6967 |
| 5 | exp_017_lgbm_2000_median | lgbm | 2000 | median | 0.8447 | 12.019 | 0.6825 |

v2 vs v3:

| Benchmark | Best experiment | Pearson r | MAE weeks | R2 |
|---|---|---:|---:|---:|
| v2 site best composite | exp_024_lgbm_1000_mean | 0.8404 | 12.750 | 0.6516 |
| v2 site best MAE | exp_043_elasticnet_2000_median | 0.8129 | 12.471 | 0.6326 |
| v3 region best | exp_003_lgbm_500_median | 0.8530 | 11.550 | 0.7051 |

Region v3 improves the repaired Phase 0 benchmark modestly, but the best MAE remains far above the publication target of <3.5 weeks.

## Scientific Judgment

Region aggregation helps, but only modestly on the current GSE120137-only Phase 0 matrix. This is not evidence to move to MLP/CNN/Transformer. The limiting factor is still biological/data coverage: no true multi-dataset beta matrix and no usable local GSE80672 matrix for CR validation.

Recommended next step:

1. Download or convert the GSE80672 GEO processed TXT supplement before full FASTQ ETL.
2. Build a GSE80672-compatible beta/region matrix and run held-out CR validation with real age acceleration metrics.
3. Add GSE93957/GSE121141/GSE60012 matrices to enable GroupKFold by `dataset_batch`.
4. Consider DBSCAN/conserved CpG cluster regions after the fixed 5 kb baseline, not before.

Deep learning remains out of scope until the cross-dataset matrix exists and classical site/region baselines stabilize.

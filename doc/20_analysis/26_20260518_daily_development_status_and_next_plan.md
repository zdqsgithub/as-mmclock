# 2026-05-18 Daily Development Record, Current Status, and Next Plan

Date: 2026-05-18

Project: mouse methylation chronological-age clock from raw/processed RRBS/WGBS
methylation data.

Method standard: AS-DS-Ops v3.3, RALPH-first. Data gates must pass before
autoresearch. All benchmark metrics remain research-grade only.

## Executive Summary

Today the project moved from model-search thinking into a data-gated acquisition
and diagnostic loop. The main scientific blocker is unchanged: GSE121141 old
brain_cortex/heart/lung extrapolation is not fixed by current model tuning,
region features, embedding POC, calibration, or already integrated public
datasets.

The practical progress today was:

- public processed/matrix candidates were audited through official GEO/E-Utils,
  SOFT, filelist, sample supplement paths, and SRA/BioSample checks;
- v11.2 proved that GSE286302 processed COV files are technically compatible
  with the current 5kb region representation;
- v11.3 started a robust background download backlog for auxiliary/diagnostic
  processed methylation datasets;
- the project now has a reproducible path to build auxiliary matrices after the
  background downloads finish.

Current decision:

- Do not run autoresearch yet.
- Do not start FASTQ/Bismark ETL yet.
- Do not promote GSE286302/GSE281602/GSE224442/GSE92486/GSE129712/GSE175410 to
  headline chronological benchmark without gates.
- Continue the v11.3 processed-download backlog and prepare schema smoke /
  conversion queues.

## AS-DS-Ops Rules Applied Today

The following project rules were enforced:

- RALPH before autoresearch.
- Official metadata first: GEO/E-Utils/SOFT/filelist/FTP/SRA/BioSample.
- E-Utils is metadata only; actual GEO supplements must use official series or
  sample supplement URLs.
- Sample-level supplement fallback is required when series-level individual
  file URLs return `404`.
- RRBS is not deduplicated.
- No human clock CpG mapping.
- No dummy CR/rapamycin AUC.
- No full-data feature selection.
- No biological-age claim without real held-out predictions.
- Diagnostic-only data must stay out of headline chronological benchmarks.

## Development Record

### v2-v4: Architecture and Feature Foundations

Earlier today the project had already stabilized the core no-leakage training
architecture:

- metadata parsing was repaired and made traceable;
- feature selection was moved inside CV folds;
- prediction schema and benchmark units were standardized;
- region-based 5kb matrix construction was added;
- embedding-aware POC was defined as two tracks:
  sample-level prediction and region/feature-level interpretation.

Outcome:

- region features became the main feature representation;
- embedding was retained as an interpretation layer, not as the main way to fix
  cross-dataset generalization.

### v5-v7.1: Multi-Dataset Matrix and Harmonization

The project then moved to processed GEO supplements and multi-dataset matrices:

- GSE80672, GSE93957, GSE121141, GSE60012, and later GSE213628 were connected or
  audited through processed data paths where possible;
- GroupKFold and held-out validation replaced Phase 0-style KFold as the main
  cross-dataset benchmark;
- v7/v7.1 tested fold-internal preprocessing, feature presence, target
  transforms, GSE60012 synthetic header metadata, and GSE80672 CR validation.

Important baseline metrics retained in project reports:

- v6 GroupKFold: `r=0.5812`, `MAE=24.351w`, `R2=0.3096`.
- v6 GSE80672 CR held-out: `MAE=28.441w`, `CR AUC=0.6616`.
- v7.5 core4 baseline: `GroupKFold MAE=23.767w`.
- GSE121141 old104+ remained the main failure mode, with old104+ MAE near the
  `70w+` range in multiple diagnostics.

Outcome:

- harmonization helped some metrics but did not solve GSE121141 old-age
  same-tissue extrapolation.

### v8.x: Same-Tissue Old-Age Diagnostics

v8.x focused on the hypothesis that GSE121141 old brain_cortex/heart/lung
failure was caused by insufficient same-tissue old-age support.

Actions:

- GSE213628 was integrated and then assembly-aware/liftover issues were
  diagnosed and repaired;
- matrix ablations tested core4, core4+GSE213628, core4+GSE60012, and all6;
- tissue/age-support diagnostics showed that adding available auxiliary data
  did not reliably repair GSE121141 old104+ extrapolation;
- calibration was tested only as a diagnostic layer and was not promoted.

Outcome:

- v8.x did not justify deep learning or broader autoresearch;
- the bottleneck remained data support, not local model tuning.

### v9-v10: Public Candidate Search and Failure Gates

v9/v10 broadened the official public-data search.

Actions:

- candidate discovery used official GEO/E-Utils/SOFT/FTP/filelist rules;
- candidates were tiered as P1/P2/P3/P4;
- no broad downloads were allowed before data gates;
- success/failure definitions were formalized:
  GSE121141 old104+ improvement, all-age held-out MAE, GroupKFold MAE, random
  label sanity, and shuffled CR sanity.

Outcome:

- no new P1/P3 public processed headline dataset was found;
- v10 concluded that another blind public processed search was unlikely to fix
  the main blocker.

### v11.0-v11.1: Raw-Lead RFC and BioSample Verification

v11 moved from processed-only discovery to raw-lead RFC without immediately
starting raw FASTQ/Bismark.

Outputs:

- `results/ralph_v11_data_strategy/literature_sra_source_records.csv`
- `results/ralph_v11_data_strategy/raw_accession_rfc_leads.csv`
- `results/ralph_v11_data_strategy/raw_project_rfc_summary.csv`
- `results/ralph_v11_data_strategy/raw_biosample_verification_manifest.csv`
- `results/ralph_v11_data_strategy/raw_biosample_decision_state.json`

Result:

- `GSE286302 / PRJNA1208582 / SRP556367`: valid old-lung diagnostic lead, but no
  exact age_weeks; not headline.
- `GSE281602 / PRJNA1184779 / SRP544506`: exact 4mo/28mo heart/cardiomyocyte
  labels, but cell-type specific; not bulk heart headline.

Decision:

- v11.1 state: `diagnostic_only`.
- no chronological pilot rows passed.
- no FASTQ/Bismark pilot was authorized.

### v11.2: GSE286302 Processed-COV Diagnostic Pilot

v11.2 used processed COV files instead of FASTQ/Bismark because official GEO
sample supplements were available and local Bismark/Bowtie2/Samtools were not
installed.

Key operational finding:

- the series supplement directory exposed only `GSE286302_RAW.tar` and
  `filelist.txt`;
- direct series-level individual COV URLs returned `404`;
- the correct official path was sample-level:
  `/geo/samples/GSM8723nnn/<GSM>/suppl/<filename>`.

Outputs:

- `scripts/validate/run_v11_2_gse286302_cov_pilot.py`
- `results/v11_2_gse286302_cov_pilot/GSE286302_beta_matrix.parquet`
- `results/v11_2_gse286302_cov_pilot/GSE286302_region_matrix_5kb.parquet`
- `results/v11_2_gse286302_cov_pilot/common_regions_with_v8_2.csv`
- `results/v11_2_gse286302_cov_pilot/v11_2_decision_state.json`
- `doc/20_analysis/25_20260518_v11_2_gse286302_processed_cov_diagnostic_pilot_report.md`

Result:

| Metric | Value |
| --- | ---: |
| Downloaded COV files | 2 |
| Parsed samples | 2 |
| Known age overlap | 0 |
| CpG rows after presence filter | 17,206,857 |
| 5kb regions | 113,275 |
| Common regions with v8.2 | 63,591 |

Decision:

- technical matrix/overlap gate passed;
- scientific chronological benchmark gate failed because exact age_weeks is
  unavailable;
- GSE286302 remains diagnostic-only.

### v11.3: Background Processed Download Backlog

v11.3 was started to opportunistically download auxiliary/diagnostic processed
methylation files while preserving the RALPH gate discipline.

Scripts added:

- `scripts/etl/19_build_download_backlog_v11_3.py`
- `scripts/etl/20_download_backlog_watchdog.py`

Outputs:

- `results/download_backlog_v11_3/download_manifest.csv`
- `results/download_backlog_v11_3/download_skipped_manifest.csv`
- `results/download_backlog_v11_3/download_manifest_summary.json`
- `results/download_backlog_v11_3/download_log.jsonl`
- `results/download_backlog_v11_3/download_status.csv`
- `results/download_backlog_v11_3/watchdog_state.json`
- `results/download_backlog_v11_3/watchdog.pid`

Datasets in download backlog:

| Dataset | Role | Headline allowed now |
| --- | --- | --- |
| GSE129712 | intestine TXT adapter/smoke auxiliary | no |
| GSE175410 | skeletal_muscle old-age auxiliary | no |
| GSE224442 | parabiosis blood/liver intervention auxiliary | no |
| GSE281602 | cardiomyocyte/heart cell-type auxiliary | no |
| GSE286302 | old lung/liver/muscle diagnostic-only | no |
| GSE92486 | liver dietary restriction auxiliary | no |

Datasets intentionally skipped:

- GSE213628: already integrated/diagnosed.
- GSE213723: superseries duplicate; use subseries.
- GSE233879/GSE276335/GSE47815: HSC/cell-type context.
- GSE266961: large mixed-species/non-target low-age dataset.
- GSE295059: organoid/in-vitro context.
- GSE231658: skin/BW schema, not target gap.
- GSE312263/GSE313655: insufficient old target/schema support.

Background process:

- PID: recorded in `results/download_backlog_v11_3/watchdog.pid`.
- Started with `setsid -f` so the process survives the interactive shell.
- Uses 3 dataset workers and 8 HTTP Range chunks per active file.
- Uses structured logs, retry, resume, and size verification.

Snapshot at approximately 2026-05-18 15:09 Asia/Shanghai:

| Metric | Value |
| --- | ---: |
| Total expected bytes | 12,780,766,370 |
| Downloaded estimate | 1,862,388,197 |
| Verified bytes | 1,736,559,077 |
| Already verified existing files | 2 |
| Newly verified files | 64 |
| Active downloads | 3 |
| Remaining planned files | 123 |
| Window speed | 0.75 MiB/s |
| ETA | 3h52m |

The downloader is working, but NCBI occasionally returns `503`. The watchdog
records these attempts and resumes chunks rather than restarting entire files.

## Current Project Status

### What Is Working

- Fold-internal no-leakage training architecture exists.
- Standard prediction/benchmark schema exists.
- Region-based 5kb features are supported.
- Multi-dataset matrices and assembly-aware/liftover workflows exist.
- GSE80672 CR validation exists as held-out research-grade validation.
- Background GEO processed download is robust and resumable.
- Sample-level supplement fallback is now encoded in project SOP.

### What Remains Blocked

- No current dataset fixes GSE121141 old brain_cortex/heart/lung old104+
  extrapolation.
- GSE286302 lacks exact age_weeks.
- GSE281602 is cardiomyocyte/cell-type specific, not bulk heart.
- GSE224442/GSE92486/GSE129712/GSE175410 are auxiliary or single-tissue and do
  not directly solve the headline target.
- FASTQ/Bismark is blocked because no candidate currently passes the
  chronological P3 gate and local Bismark/Bowtie2/Samtools are unavailable.
- Autoresearch remains blocked because data gates are unresolved.

### Current Scientific Interpretation

The best explanation remains:

1. Model architecture was not the primary bottleneck after the v2-v4 fixes.
2. Region-based features are appropriate and technically stable.
3. Cross-dataset old-age failure is dominated by data/tissue/age-range support,
   schema/coverage shift, or study-specific biology.
4. More tuning on the current matrix is unlikely to produce a defensible
   headline solution.
5. The next useful work is controlled data acquisition, schema smoke, matrix
   gate checks, and only then fixed benchmarks.

## Immediate Next Plan

### Step 1: Let v11.3 Downloads Finish

Do not interrupt the watchdog unless disk or process health fails.

Monitor:

```bash
cat results/download_backlog_v11_3/watchdog_state.json
tail -f results/download_backlog_v11_3/download_log.jsonl
ps -p $(cat results/download_backlog_v11_3/watchdog.pid)
```

Stop only if:

- disk free space becomes unsafe;
- watchdog exits unexpectedly;
- a dataset repeatedly fails all retries;
- NCBI blocks too many connections and progress stalls for hours.

### Step 2: Build Conversion Queue

After downloads finish or a dataset completes:

- create `scripts/etl/21_build_v11_3_conversion_queue.py`;
- read `results/download_backlog_v11_3/download_status.csv`;
- group only fully verified datasets;
- output `results/download_backlog_v11_3/conversion_queue.csv`;
- do not convert partially downloaded datasets.

### Step 3: Run Schema Smoke

Create or run a light schema smoke script:

- input: verified files per dataset;
- check gzip readability;
- inspect first rows;
- infer COV/TXT/BEDGRAPH-like schema;
- verify coordinate columns and methylation value range;
- report whether current converter can handle it.

Expected output:

- `results/download_backlog_v11_3/schema_smoke_summary.csv`
- `results/download_backlog_v11_3/dataset_download_gate_summary.csv`

### Step 4: Convert Only Gated Datasets

Allowed first conversion order after full verification:

1. GSE175410: small COV, skeletal muscle auxiliary.
2. GSE281602: small COV, cardiomyocyte auxiliary.
3. GSE286302: COV diagnostic-only.
4. GSE224442: COV intervention auxiliary.
5. GSE92486: TXT adapter after smoke.
6. GSE129712: TXT adapter after smoke.

Each conversion must produce:

- `{dataset}_beta_matrix.parquet`
- `{dataset}_region_matrix_5kb.parquet`
- `{dataset}_region_stats_5kb.csv`
- `{dataset}_matrix_manifest.json`
- region overlap summary against v8.2/reference matrix.

### Step 5: Gate Before Benchmark

Only if a dataset passes these checks should it be considered for any fixed
benchmark:

- sample metadata can be joined without synthetic overclaim;
- exact age_weeks exists for chronological training/held-out use;
- tissue context is compatible with the benchmark question;
- common 5kb regions are `>=50,000`;
- schema and assembly are traceable;
- no intervention metric is computed without real held-out predictions.

### Step 6: Decide v12 Direction

After v11.3 conversion gates:

- If no exact-age same-tissue old brain_cortex/heart/lung dataset is found,
  declare public processed auxiliary expansion insufficient.
- Then choose one v12 path:
  1. targeted manual acquisition/contact/data-generation RFC;
  2. minimal FASTQ/Bismark pilot only if a P3 candidate finally passes;
  3. benchmark redefinition to separate chronological clock from intervention
     diagnostic auxiliary analyses.

## Do Not Do Next

- Do not start broad autoresearch.
- Do not run MLP/CNN/Transformer.
- Do not claim biological age effects from metadata labels.
- Do not merge auxiliary/cell-type/organoid data into headline benchmark without
  gate approval.
- Do not launch full FASTQ/Bismark ETL while current candidate gate is blocked.

## Final State for Today

The project is in an active v11.3 acquisition state:

- background downloads are running;
- current work should be light, planning/schema-oriented, and IO-aware;
- the next decisive event is not a model result but a dataset gate result after
  verified downloads and schema smoke.

## 23:52 CST Addendum: v11.4 Download Gate Completed

The v11.3 background backlog has now completed and v11.4 schema smoke has been
run.

Download hard verification:

- `192/192` manifest rows verified.
- `0` failed or missing rows.
- total verified size: `12,780,766,370` bytes (`11.903 GiB`).

Authoritative v11.4 outputs:

- `results/download_backlog_v11_3/download_completion_summary.json`
- `results/download_backlog_v11_3/dataset_download_gate_summary.csv`
- `results/download_backlog_v11_3/conversion_queue_smoke_gated.csv`
- `doc/20_analysis/27_20260518_v11_4_download_gate_schema_smoke_report.md`

Current conversion gate:

- Ready for current Bismark-COV conversion:
  `GSE224442`, `GSE281602`, `GSE286302`, `GSE92486`.
- Adapter required before conversion:
  `GSE129712` requires `methylratio_cg_12col_to_cov`;
  `GSE175410` requires accession chromosome normalization.

Updated next step:

1. Run the smallest direct-conversion smoke on `GSE281602`.
2. If valid, process `GSE224442`, `GSE92486`, and `GSE286302`.
3. Implement adapters for `GSE129712` and `GSE175410`.
4. Only after matrix, metadata, and common-region gates pass, decide whether any
   dataset can enter a fixed benchmark.

# v27 Data Repair Execution Report

Date: 2026-05-21

## Scope

This v27 pass implements the approved data-repair plan:

- repair GSE80672/GSE121141 missing raw runs after v25 ETL stability;
- pause GSE93957 as an ETL/schema/assembly/parser problem rather than a download problem;
- repair GSE60012 GSM-to-tile metadata mapping;
- start processed parser pilots for GSE225166/GSE225173, GSE233734, GSE304754, GSE175410, GSE224442, and GSE92486;
- keep modeling/DL downstream of matrix gates and schema smoke checks.

## Implemented

New scripts:

- `scripts/etl/26_repair_gse60012_gsm_tile_mapping.py`
- `scripts/validate/run_v27_gse93957_root_cause_audit.py`
- `scripts/etl/27_build_v27_processed_candidate_backlog.py`
- `scripts/etl/28_build_v27_raw_missing_repair_queue.py`
- `scripts/etl/29_run_v27_raw_missing_repair_queue.py`
- `scripts/etl/30_build_v27_processed_pilot_manifest.py`
- `scripts/validate/run_v27_data_repair_execution.py`

Key outputs:

- `results/v27_data_repair_execution/v27_decision_state.json`
- `results/v27_data_repair_execution/gse60012_mapping/gse60012_mapping_summary.json`
- `results/v27_data_repair_execution/gse93957_root_cause/gse93957_root_cause_summary.json`
- `results/v27_data_repair_execution/raw_repair/v27_raw_missing_repair_queue.csv`
- `results/v27_data_repair_execution/processed_pilot/v27_processed_pilot_download_manifest.csv`
- `results/v27_data_repair_execution/v27_background_launch_state.json`

## Results

GSE60012 mapping repair passed:

- header metadata rows: 152
- mapped rows: 152
- mapping fraction: 1.0
- locally complete FASTQ-mapped rows: 100

GSE93957 is paused:

- audited completed samples: 4
- median Bismark mapping efficiency: 0.0%
- median parsed regions: 0
- conclusion: do not download more first; recheck trimming/adapters/library/assembly/parser.

Processed candidate discovery:

- discovered datasets: GSE225166, GSE225173, GSE233734, GSE304754
- full candidate manifest rows: 471
- known planned bytes: 23.32 GB
- GSE225166/GSE225173 are sibling records with duplicated individual COV names; pilot uses sample-level GEO URLs for individual COV files.

Processed pilot:

- manifest rows: 148
- already verified rows from v11.3: 111
- planned new download rows: 37
- planned new bytes: 2.32 GB
- datasets covered: GSE175410, GSE224442, GSE225166, GSE233734, GSE304754, GSE92486

Raw missing repair:

- queued priority missing runs: 218
- GSE80672: 114 runs, about 207 GiB compressed FASTQ
- GSE121141: 104 runs, about 424 GiB compressed FASTQ
- all queued runs have ENA FASTQ URLs, MD5, and byte sizes.

## Background Jobs

v25 raw ETL remained stable when launch decisions were made:

- status: running
- completed samples: 6
- failed samples: 0
- current queue: `raw_etl_queue_v25_authorized_full_passed_gate.csv`

Started low-priority background jobs:

- processed pilot download v3: `results/v27_data_repair_execution/processed_pilot_download_v3`
- raw missing repair: `results/v27_data_repair_execution/raw_repair_download`

The first processed-pilot launch exposed a real GEO URL issue: GSE225166 filelist entries 404 at series-level individual-file URLs. The manifest builder was fixed to use sample-level GEO supplement URLs, generated partial chunks were cleaned, and the processed-pilot downloader was restarted as v3.

## Gates

Training and DL optimization are still gated. They should start only after:

- processed pilot downloads hard-verify and pass schema smoke;
- raw-derived matrices keep beta in [0, 1], primary autosome traceability, and common 5kb region gates;
- GSE80672 CR and GSE121141 stress-test support are computed from real held-out predictions.

Current biological interpretation status: not yet upgraded. Existing v22/v23 models show useful age-associated methylation signal, but v27 is still in data expansion/repair. Do not claim new biological mechanism until v27 matrices and feature-confounding audits pass.

# v23 Data Resource Maximization SOP

Date: 2026-05-21

## Purpose

Use all available mouse methylation resources to discover biologically meaningful methylation signals, with biological-age prediction as a priority when the data support it. This SOP clarifies that raw FASTQ gates authorize staged expansion after success; they are not a permanent restriction on raw-data use.

## Resource Split

- **Local workstation / RAID**: FASTQ inventory, checksum/size audit, sample-run mapping, trimming, Bismark/Bowtie2/Samtools, COV generation, 5kb matrix construction, raw-vs-processed consistency checks, metadata repair, and CPU ML.
- **Cloud GPU host**: DL training, DL architecture comparisons, final all-data DL fits, random-label DL sanity, and large neural-search batches.
- **Notify user before GPU change**: only if a non-DL raw-data step truly requires GPU acceleration. Standard FASTQ, FASTA indexing, Bismark, Bowtie2, and Samtools are CPU/IO workflows.

## Authorization Ladder

1. Read-only raw inventory is authorized.
2. 2-3 sample raw pilot is authorized when sample/run metadata and assay/reference context are mapped.
3. Dataset-level raw ETL expansion is authorized after the pilot matrix gate passes.
4. Local ML and cloud DL on raw-derived matrices are authorized after matrix gate success and sanity-check planning.
5. If three pilots fail on critical gates, stop expansion and enter re-evaluation.

## Required Raw Matrix Gates

- metadata overlap `>=95%`;
- exact age coverage `>=95%`;
- beta values in `[0,1]`;
- assembly and region schema traceable;
- sex chromosomes and mitochondrial regions excluded or traceably flagged;
- common autosomal 5kb regions `>=50000`;
- raw-vs-processed consistency audit when a same-dataset processed baseline exists.

## Robust Long-Run Rules

Every long raw ETL or cloud DL run must write:

- run manifest;
- task manifest;
- heartbeat;
- event log;
- per-task logs;
- status JSON with PID, retry count, current task, exit code, and timestamp;
- summary CSV/JSON and markdown report.

Use background execution (`nohup`, watchdog, or scheduler). Interactive SSH is for smoke tests only.

## Reporting Guardrails

Reports may claim age-associated, CR-associated, tissue-associated, or condition-associated methylation signals when supported by held-out or leakage-controlled predictions. They must not claim aging mechanisms from predictive regions alone, especially when regions are labeled tissue-, dataset-, or coverage-confounded.

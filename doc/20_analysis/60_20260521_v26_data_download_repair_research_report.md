# v26 Data Download Provenance, Repair, and Expansion Research Report

Date: 2026-05-21

## Scope

This report audits how the current local mouse methylation datasets were
downloaded, identifies repairable gaps in the existing raw FASTQ holdings, and
prioritizes additional public datasets for biological-signal ML/DL work. No
new FASTQ, COV, BEDGRAPH, or GEO supplement downloads were started in this
work.

The active local raw ETL watchdog was left running. At the time of this audit,
the v25 queue was processing the already gate-passed families (`GSE80672` and
`GSE121141`) under `/data/mouse_methyl/processed_v21_raw_etl`.

## Local Download Provenance

The existing raw FASTQ datasets under `/data/mouse_methyl/raw/GSE*/fastq` were
downloaded by the older AS-DS-Ops SRA downloader in
`/home/zdq-as/as-project1/pipeline/geo_sra_downloader.py`, with orchestration
settings in `/home/zdq-as/as-project1/pipeline/config.json`.

The observed workflow was:

1. Query GEO/SRA run metadata with `pysradb`.
2. Persist `metadata/sra_metadata.json` and `metadata/sra_metadata.tsv`.
3. Run SRA Toolkit `prefetch --max-size u -O <ssd_cache> <SRR>`.
4. Run `fasterq-dump --split-files -e 12 -O <run_cache> -t <run_cache> <SRR>`.
5. Keep structured attempt state in `metadata/download_log.jsonl`.
6. Use SSD cache `/mouse_methyl_work/raw_downloads/<GSE>/<SRR>/` and final raw
   storage `/data/mouse_methyl/raw/<GSE>/fastq/`.

This matches NCBI SRA Toolkit guidance that `prefetch` + `fasterq-dump` is the
standard fast path and that failed `prefetch` runs can resume rather than start
from scratch.

Processed GEO supplements were handled separately by the repo v11.3 downloader:

- manifest builder: `scripts/etl/19_build_download_backlog_v11_3.py`
- watchdog downloader: `scripts/etl/20_download_backlog_watchdog.py`
- method: official GEO HTTP URLs, HTTP Range chunking, resume, size checks, and
  structured watchdog/status files

## Current Local Raw Holdings

| dataset | SRA runs | complete local runs | missing / incomplete | FASTQ GiB | logged completed | logged failed |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| GSE120137 | 739 | 573 | 166 | 2611.194 | 572 | 4 |
| GSE121141 | 116 | 12 | 104 | 175.848 | 8 | 2 |
| GSE45361 | 23 | 22 | 1 | 147.987 | 22 | 0 |
| GSE52266 | 128 | 128 | 0 | 181.425 | 128 | 0 |
| GSE60012 | 173 | 109 | 64 | 735.219 | 109 | 0 |
| GSE80672 | 380 | 266 | 114 | 1833.056 | 265 | 2 |
| GSE80761 | 4 | 4 | 0 | 46.848 | 4 | 0 |
| GSE93957 | 62 | 59 | 3 | 842.997 | 58 | 1 |

Detailed outputs:

- `results/v26_data_repair_research/local_raw_download_provenance_summary.csv`
- `results/v26_data_repair_research/local_raw_download_log_flat.csv`
- `results/v26_data_repair_research/local_raw_run_fastq_completeness.csv`
- `results/v26_data_repair_research/local_raw_missing_fastq_repair_queue.csv`
- `results/v26_data_repair_research/local_raw_fastq_completeness_summary.csv`

## Repair Findings

The existing gaps are mostly repairable, but should not all be repaired at once.
The active bottleneck is not download tooling alone; it is whether each dataset
can pass metadata, assembly, beta-range, and common-region gates after ETL.

| dataset | priority | finding | action |
| --- | --- | --- | --- |
| GSE80672 | P0 | 23 raw-derived samples already passed matrix gate; v25 full local queue is running; 114 SRA runs are still missing/incomplete. | Let v25 ETL finish first. Then retry metadata-relevant missing CR/control runs; use one SRA retry, then ENA direct FASTQ + md5 fallback. |
| GSE121141 | P0 | 7 raw-derived brain cortex samples passed gate; 104 SRA runs missing/incomplete; more age/tissue coverage is needed before raw DL. | Fill older ages/tissues after current ETL load drops. Retry `SRR8032931` and `SRR8032940` with lower `fasterq-dump` threads or ENA fallback. |
| GSE93957 | P1 stop/re-evaluate | Four ETL attempts produced fewer than 1000 COV rows and zero usable 5kb autosomal regions. | Do not download more yet. Inspect original processing notes, Bismark reports, adapters, library type, assembly, and COV parser assumptions. |
| GSE60012 | P1 metadata repair | 109 raw GSMs exist locally, but current training metadata uses tile-style IDs, so exact overlap is zero. | Repair GSM-to-tile/sample mapping before raw expansion or training. |
| GSE120137 | P2 | 573 single-end runs exist and 166 are missing; processed baseline exists; current raw runner is paired-first. | Add/test single-end Bismark queue support on 2-3 samples, then raw-vs-processed consistency audit before filling missing runs. |
| GSE52266 | P2 | No missing FASTQ. Single-end/MBD-Seq metadata requires careful filtering. | Use later as a single-end ETL validation dataset. |
| GSE45361 | P3 | One missing run, but age labels are absent in current metadata. | Treat as non-age auxiliary unless phenotype labels are curated. |
| GSE80761 | P3 | Complete but only four runs. | Use only as small schema/sanity support. |

Repair action table:

- `results/v26_data_repair_research/data_repair_action_plan.csv`

## New Data Expansion Candidates

The best immediate expansion is processed COV/BEDGRAPH data first, not raw
FASTQ. Processed downloads are small enough to parse quickly and decide whether
raw ETL is worth the CPU/IO cost.

| priority | dataset | why useful | recommended mode |
| --- | --- | --- | --- |
| P0 | `GSE225166/GSE225173` | 232 methylation samples in the subseries, multi-tissue mouse ages 10-101 weeks, designed around sparse/single-cell age prediction (`scEpiAge`). | Download processed `GSE225166_RAW.tar` COV first; build sparse/low-coverage parser. |
| P0 | `GSE233734` | 102 colon RRBS samples, 3-28 month lifespan trajectory, nonlinear stage-of-aging target, processed methylation/coverage BEDGRAPH available. | Download processed BEDGRAPH first; build parser and stage target. |
| P1 | `GSE304754/PRJEB73981` | Brain/hippocampus/cognitive aging methylation, linked to multi-omics and cognitive performance; ENA lists FASTQ files. | Search for processed matrix first; otherwise ENA 2-3 sample raw pilot. |
| P1 | `GSE175410` | Late-life exercise skeletal muscle epigenetic aging; processed Bismark COV exists; high biological intervention value. | Use existing v11.3 processed COV backlog. |
| P1 | `GSE224442` | Blood/liver parabiosis and recovery, 78 RRBS samples, processed COV tar. | Use existing v11.3 processed COV backlog. |
| P1 | `GSE92486` | Liver dietary restriction/ad lib, young/old BS-seq methylation plus RNA-seq; strong CR-like validation target. | Download processed TXT/COV-like tar; filter BS-seq samples. |
| P2 | `GSE129712` | Intestinal aging RRBS, 24 samples, 2/18/26 months, processed TXT tar. | Processed parser pilot. |
| P2 | `GSE286302` | Very recent circadian/PVN/3dA healthspan intervention with raw and processed RRBS deposited. | Processed COV first if accessible; intervention validation, not primary age training. |
| P3 | `GSE295059` | Intestinal organoid aging/passaging/decitabine, 70 RRBS samples, processed COV tar. | Auxiliary only; mm9/liftover and organoid confounding required. |

Candidate table:

- `results/v26_data_repair_research/additional_data_candidate_plan.csv`

## Open Source / Model Scan

There are relevant open projects, but most consume methylation matrices rather
than raw FASTQ. The project should keep its local FASTQ-to-region ETL, then
adapt selected model ideas.

| project | relevance |
| --- | --- |
| MouseEpigeneticClock | Direct mouse age-clock comparator; 329 CpG sites; original paper processed FASTQ with Trim Galore/Bismark. |
| scEpiAge | Mouse sparse/single-cell DNA methylation age model for liver/lung/blood; useful after `GSE225166/GSE225173` parser. |
| MethylNet | PyTorch VAE + classifier/regressor + SHAP for methylation matrices; useful DL/interpretable architecture reference. |
| MethFormer | HuggingFace transformer for binned methylation masked regression; closest public pattern for self-supervised methylation representation once raw matrices scale. |
| XAI-AGE | Pathway-aware explainable age DL; human-focused but useful as an interpretation design pattern. |
| AltumAge / methylclock / dnaMethyAge | Useful human-array comparators; not immediately usable for mouse RRBS raw matrices. |
| YAME | Efficient sequence-level methylome storage/summarization; possible future performance layer, not a clock model. |

Project scan table:

- `results/v26_data_repair_research/open_source_methylation_signal_project_scan.csv`

## Recommended Execution Order

1. Let the current v25 local ETL queue finish or reach a stable checkpoint.
2. In parallel, implement light parser pilots for processed `GSE225166` and
   `GSE233734`; they are high-value and do not require raw FASTQ ETL first.
3. Run existing v11.3 processed-supplement backlog for `GSE175410`,
   `GSE224442`, `GSE92486`, and `GSE129712` once local IO is not saturated.
4. Repair `GSE60012` metadata mapping before any GSE60012 training claim.
5. Add single-end ETL support, then pilot `GSE120137` and `GSE52266` raw
   samples.
6. Re-evaluate `GSE93957` before any more download/ETL.
7. Only after sample count increases and sanity checks pass, package raw-derived
   matrices for cloud DL.

## Decision

Yes, more useful public methylation data can be downloaded, and yes, existing
data can be repaired. The safest high-yield path is:

- processed COV/BEDGRAPH expansion first (`GSE225166`, `GSE233734`, existing
  v11.3 backlog);
- targeted raw FASTQ repair second (`GSE80672`, `GSE121141`);
- single-end ETL support third (`GSE120137`, `GSE52266`);
- no raw DL expansion until sample count and gate-passed matrix coverage are
  larger.

The current biologically meaningful signal remains the gate-passed GSE80672
CR-associated methylation signal. The age-clock track needs more raw-derived
samples or the new processed age-rich datasets before DL optimization is worth
another cloud run.

## Sources

- NCBI SRA Toolkit `prefetch` + `fasterq-dump` workflow:
  https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump
- ENA direct FASTQ and md5 download structure:
  https://ena-docs.readthedocs.io/en/latest/retrieval/file-download.html
- ENA SRA FTP structure:
  https://ena-docs.readthedocs.io/en/latest/retrieval/file-download/sra-ftp-structure.html
- `GSE225166`: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225166
- `GSE233734`: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233734
- `GSE304754/PRJEB73981`: https://www.omicsdi.org/dataset/project/PRJEB73981
- `GSE175410` sample evidence: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5333080
- `GSE224442`: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224442
- `GSE92486`: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92486
- `GSE129712`: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129712
- `GSE286302` Cell article: https://www.sciencedirect.com/science/article/pii/S0092867426001030
- `GSE295059`: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE295059
- MouseEpigeneticClock: https://github.com/EpigenomeClock/MouseEpigeneticClock
- scEpiAge: https://github.com/EpigenomeClock/scEpiAge
- MethylNet: https://github.com/Christensen-Lab-Dartmouth/MethylNet
- MethFormer: https://huggingface.co/CChahrour/Methformer
- XAI-AGE: https://github.com/Paureel/XAI-AGE
- ComputAge Bench: https://huggingface.co/datasets/computage/computage_bench
- YAME: https://github.com/zhou-lab/YAME

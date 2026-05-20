# v11 Minimal ETL RFC SOP

Date: 2026-05-18

## Purpose

v11 starts after v10 fails to find a directly usable P1/P3 public GEO
processed dataset. v11 may consider raw FASTQ/Bismark only as a controlled
minimal pilot, and only after official metadata confirms sample-specific age,
target tissue, assay, and traceable run accessions.

This SOP does not approve full ETL. It defines the gate for deciding whether a
2-3 sample pilot is scientifically justified.

## Allowed Inputs

- `results/ralph_v11_data_strategy/raw_accession_rfc_leads.csv`
- `results/ralph_v11_data_strategy/raw_project_rfc_summary.csv`
- Official NCBI SRA RunInfo or E-Utils metadata
- Official GEO GSM/GSE text records
- Official BioSample metadata

## Candidate Priority

Only candidates with all of the following may enter a pilot RFC:

- mouse bulk RRBS/WGBS/bisulfite sequencing;
- old target tissue: brain_cortex/cortex, heart, or lung;
- sample-specific age or age group that can be mapped to weeks;
- sample-specific tissue;
- run accessions with public FASTQ availability;
- no single-cell, targeted cell type, organoid, cell-line, iTAG, or RNA-seq
  context;
- enough same-tissue young/old or old-control samples to make the pilot
  interpretable.

## Current v11 Pilot RFC Candidates

| Candidate | Evidence | Status |
| --- | --- | --- |
| GSE286302 / PRJNA1208582 / SRP556367 | SRA titles include old lung bisulfite runs, e.g. aged vehicle and aged 3dA replicates | BioSample/RunInfo verification required |
| GSE281602 / PRJNA1184779 / SRP544506 | SRA/GSM records resolve to GSE281602 and include heart bisulfite/RRBS-like records with age labels | BioSample/RunInfo verification required |

Previously tested or lower-priority raw leads, such as PRJNA882035 linked to
GSE213628/GSE213723, are not first-choice pilots because they overlap with
datasets already integrated or diagnosed in v8.

## Pilot Gate

Before any FASTQ download, create a verification manifest with:

- project accession, study accession, experiment accession, run accession;
- sample accession and BioSample accession;
- sample title/source;
- age token and normalized age_weeks;
- tissue;
- sex and strain if available;
- assay/library strategy;
- layout and expected FASTQ count;
- total size per run;
- reason for inclusion or exclusion.

Pilot may proceed only if at least 2-3 samples pass verification and the design
contains interpretable age/tissue support. Preferred design is same-tissue
young/old. Old-only same-tissue support can be used only as a diagnostic, not as
a headline chronological benchmark.

## Minimal Pilot Scope

If the verification manifest passes:

1. Download only 2-3 runs with `prefetch`/ENA fallback and resume logging.
2. Convert to FASTQ with `fasterq-dump`.
3. Run Bismark using the existing project reference/assembly policy.
4. Build coverage and 5kb region matrix for pilot samples only.
5. Estimate overlap with v8.2/v10 matrices before any full download.

## Stop Conditions

Stop before FASTQ download if:

- BioSample lacks sample-specific age or tissue;
- the assay is RNA-seq, single-cell, tagged/low-coverage, cell-type-specific,
  organoid, or in vitro;
- no same-tissue old target support remains after verification;
- run sizes make a 2-3 sample pilot impractical without a separate resource
  approval;
- common 5kb region estimate is unlikely to reach `>=50000`.

## v11.1 Outcome

BioSample/RunInfo verification demoted the two current candidates:

- `GSE286302 / PRJNA1208582 / SRP556367`: valid old lung RRBS/bisulfite
  diagnostic-only candidate, but no exact age_weeks and not a chronological
  headline dataset.
- `GSE281602 / PRJNA1184779 / SRP544506`: exact 4mo/28mo heart labels, but
  BioSample/GEO identifies cardiomyocyte/cell-type samples, not bulk heart.

Thus no v11.1 candidate currently authorizes a headline chronological pilot.
Only an optional GSE286302 diagnostic old-lung pilot remains scientifically
defensible.

## v11.2 Processed-COV Pilot Outcome

The optional GSE286302 pilot was run as a processed-COV-first diagnostic, not
as FASTQ/Bismark ETL. This was the correct route because official processed
`.bismark.cov.gz` sample supplements exist and the local Bismark/Bowtie2/Samtools
toolchain is not installed.

Operational rule added:

- If GEO `filelist.txt` lists sample files but direct series supplement URLs
  return `404`, resolve individual processed files through official sample-level
  supplement paths:
  `/geo/samples/GSMxxxxnnn/<GSM>/suppl/<filename>`.
- Keep the failed series-level URL attempts in the download log for traceability.
- Prefer resumable range downloads for large processed files when single
  connections are unstable or slow.
- Do not download the full series `RAW.tar` when the needed sample-level
  supplements are available.

v11.2 result:

- 2 GSE286302 lung COV files downloaded and size-verified.
- Existing converter produced `113,275` 5kb regions with `63,591` common
  regions against the v8.2 reference matrix.
- The dataset remains diagnostic-only because sample-level exact `age_weeks`
  is unavailable.

# v11.4 Download Gate and Schema Smoke Report

Date: 2026-05-18 23:52 CST

Project standard: AS-DS-Ops v3.3. RALPH-first; no autoresearch while data gates
remain unresolved.

## Summary

The v11.3 processed supplement download backlog has completed and passed hard
verification. All `192/192` manifest rows are present and size-verified, with
`12,780,766,370` verified bytes (`11.903 GiB`).

v11.4 then ran schema smoke over the verified files. No dataset is blocked at
the download/schema-smoke layer, but only four datasets can enter the existing
Bismark-COV conversion path immediately. Two datasets require explicit adapters
before conversion.

This is still a data-gate stage. None of these datasets is promoted to headline
chronological benchmark or biological-age validation yet.

## Outputs

Scripts:

- `scripts/etl/21_build_v11_3_conversion_queue.py`
- `scripts/validate/v11_3_processed_schema_smoke.py`

Primary outputs:

- `results/download_backlog_v11_3/download_completion_summary.json`
- `results/download_backlog_v11_3/download_completion_by_dataset.csv`
- `results/download_backlog_v11_3/conversion_queue.csv`
- `results/download_backlog_v11_3/schema_smoke_summary.csv`
- `results/download_backlog_v11_3/dataset_download_gate_summary.csv`
- `results/download_backlog_v11_3/conversion_queue_smoke_gated.csv`

Per-dataset verified file manifests:

- `results/v11_4_conversions/GSE129712/GSE129712_verified_files.csv`
- `results/v11_4_conversions/GSE175410/GSE175410_verified_files.csv`
- `results/v11_4_conversions/GSE224442/GSE224442_verified_files.csv`
- `results/v11_4_conversions/GSE281602/GSE281602_verified_files.csv`
- `results/v11_4_conversions/GSE286302/GSE286302_verified_files.csv`
- `results/v11_4_conversions/GSE92486/GSE92486_verified_files.csv`

## Download Verification

| Dataset | Verified files | Size GiB | File type | Current role |
|---|---:|---:|---|---|
| GSE129712 | 24 | 0.325 | TXT | P2 diagnostic |
| GSE175410 | 15 | 0.336 | COV | P2 diagnostic |
| GSE224442 | 78 | 2.260 | COV | P2 intervention auxiliary |
| GSE281602 | 8 | 0.218 | COV | P3 cell-type auxiliary |
| GSE286302 | 49 | 5.911 | COV | P3 diagnostic |
| GSE92486 | 18 | 2.853 | TXT | P2 intervention auxiliary |

Hard verification result:

- Manifest rows: `192`
- Verified rows: `192`
- Failed or missing rows: `0`
- Status: `completed`

## Schema Smoke Gates

| Dataset | Smoke schema | Gate | Adapter |
|---|---|---|---|
| GSE224442 | `bismark_cov_6col` | ready for current COV conversion | `bismark_cov_per_sample_tar` |
| GSE281602 | `bismark_cov_6col` | ready for current COV conversion | `bismark_cov_per_sample_tar` |
| GSE286302 | `bismark_cov_6col` | ready for current COV conversion | `bismark_cov_per_sample_tar` |
| GSE92486 | `bismark_cov_6col` | ready for current COV conversion | `bismark_cov_per_sample_tar` |
| GSE129712 | `methylratio_cg_12col` | adapter required before conversion | `methylratio_cg_12col_to_cov` |
| GSE175410 | `bismark_cov_6col_accession_chrom` | adapter required before conversion | `bismark_cov_with_accession_chrom_map` |

Interpretation:

- `GSE92486` uses `.txt` filenames, but smoke shows the contents are Bismark-like
  six-column coverage records. It can use the same COV converter after packaging.
- `GSE129712` has methylratio-style twelve-column records and needs a dedicated
  parser that derives beta from `ratio` or `C_count / CT_count`, while preserving
  coverage fields.
- `GSE175410` has Bismark-like six-column records, but chromosome labels are
  assembly accession IDs such as `CM...`. It must not be treated as standard
  `chrN` coordinates until a traceable accession-to-chromosome map or liftover
  rule is added.

## Conversion Queue

Use `conversion_queue_smoke_gated.csv` as the authoritative v11.4 queue. The
older `conversion_queue.csv` was generated before content-level schema smoke and
is now superseded for conversion decisions.

Immediate conversion candidates:

1. `GSE281602`: smallest COV set; useful to test the packaging/conversion path.
2. `GSE224442`: COV intervention auxiliary.
3. `GSE92486`: TXT filename but COV-like content; liver DR/AL auxiliary.
4. `GSE286302`: largest ready COV set; old lung/liver/muscle diagnostic only.

Adapter queue:

1. `GSE129712`: implement `methylratio_cg_12col_to_cov`.
2. `GSE175410`: implement accession chromosome normalization before conversion.

## Gates Before Training

After conversion, each dataset still must pass:

- metadata join with traceable sample IDs;
- exact age coverage if used for chronological benchmark;
- tissue compatibility with the benchmark question;
- sex/MT exclusion;
- common 5kb regions `>=50,000`;
- region overlap manifest against the v8.2/v11 reference matrix;
- no intervention metrics unless real held-out predictions exist.

If a dataset is cell-type, intervention-specific, synthetic, or lacks exact
age, it remains auxiliary/diagnostic even if conversion succeeds.

## Next Step

Run a minimal direct-conversion smoke on `GSE281602` first, because it is the
smallest ready set. If packaging and conversion produce valid beta, region, and
manifest outputs, process `GSE224442`, `GSE92486`, and `GSE286302` in that order.

Do not run autoresearch or GroupKFold benchmark until converted matrices pass
the metadata and common-region gates.

## 2026-05-19 Addendum: GSE281602 Smoke Result

`GSE281602` conversion smoke completed successfully as a parser and matrix
construction test.

Key result:

- samples parsed: `8`
- metadata overlap: `8`
- beta matrix: `2,779,291 x 8`
- region matrix: `121,333 x 8`
- beta range: `[0.0, 1.0]`
- parser errors: `0`
- common regions with v8.2 reference: `25,321`
- common-region gate `>=50,000`: `false`

Decision:

`GSE281602` remains auxiliary only. It validates the direct Bismark-COV
conversion path, but it is not eligible for headline training because it is
cardiomyocyte/cell-type specific and fails the common-region gate.

Detailed report:

- `doc/20_analysis/28_20260519_v11_4_gse281602_conversion_smoke_report.md`

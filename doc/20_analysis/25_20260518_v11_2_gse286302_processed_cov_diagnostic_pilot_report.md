# v11.2 GSE286302 Processed-COV Diagnostic Pilot Report

Date: 2026-05-18

## Summary

v11.2 executed the optional GSE286302 old-lung diagnostic pilot without
starting FASTQ/Bismark. The project used official GEO sample-level processed
supplements, built a 2-sample Bismark COV pilot, converted it into 5kb region
features, and measured overlap with the current v8.2 lifted multidataset
region matrix.

Final state: `diagnostic_completed`

Decision: `processed_cov_diagnostic_overlap_only_not_headline`

## Inputs

- Young lung: `GSM8723148_GY1L.bismark.cov.gz`
- Old vehicle lung: `GSM8723152_GV3L.bismark.cov.gz`
- Official sample supplement pattern:
  `/geo/samples/GSM8723nnn/<GSM>/suppl/<filename>`
- Reference overlap matrix:
  `results/multidataset_v8_2_prefilter_liftover/all_rrbs_region_matrix_5kb.parquet`

The series supplement directory only exposes `GSE286302_RAW.tar` and
`filelist.txt`; direct per-file series URLs returned `404`. The pilot therefore
used the official sample-level GEO supplement path and recorded the failed
series-path attempts in `download_log.jsonl`.

## Outputs

- `scripts/validate/run_v11_2_gse286302_cov_pilot.py`
- `results/v11_2_gse286302_cov_pilot/gse286302_cov_pilot_manifest.csv`
- `results/v11_2_gse286302_cov_pilot/gse286302_cov_pilot_metadata.csv`
- `results/v11_2_gse286302_cov_pilot/download_log.jsonl`
- `results/v11_2_gse286302_cov_pilot/GSE286302_v11_2_cov_pilot.tar`
- `results/v11_2_gse286302_cov_pilot/GSE286302_beta_matrix.parquet`
- `results/v11_2_gse286302_cov_pilot/GSE286302_region_matrix_5kb.parquet`
- `results/v11_2_gse286302_cov_pilot/GSE286302_region_stats_5kb.csv`
- `results/v11_2_gse286302_cov_pilot/common_regions_with_v8_2.csv`
- `results/v11_2_gse286302_cov_pilot/v11_2_decision_state.json`

## Conversion Result

| Metric | Value |
| --- | ---: |
| Downloaded COV files | 2 |
| Parsed samples | 2 |
| Metadata overlap | 2 |
| Known age overlap | 0 |
| CpG rows after presence filter | 17,206,857 |
| 5kb regions | 113,275 |
| sex/MT region rows | 0 |
| Converter errors | 0 |
| Conversion runtime | 126.7 sec |

## Region Overlap

| Metric | Value |
| --- | ---: |
| Pilot regions | 113,275 |
| v8.2 reference regions | 65,870 |
| Common regions | 63,591 |
| Fraction of pilot regions shared | 0.561 |
| Fraction of v8.2 reference regions shared | 0.965 |

The diagnostic pilot passes the technical overlap gate (`common_regions >=
50,000`). This means GSE286302 processed COV files are compatible with the
project's current 5kb region representation and assembly/liftover policy at a
matrix level.

## Gate Decision

GSE286302 is not promoted to headline chronological benchmark:

- It provides age group labels but no exact sample-level `age_weeks`.
- The pilot covers lung only.
- It can help diagnose old-lung feature compatibility, but it cannot answer
  the GSE121141 old104+ chronological-age target across brain_cortex/heart/lung.

Allowed next use:

- auxiliary/diagnostic old-lung matrix compatibility evidence;
- optional old-lung condition audit if more GSE286302 processed COV samples are
  downloaded later.

Disallowed next use:

- headline GroupKFold/LODO chronological age benchmark;
- biological-age or intervention claim;
- autoresearch trigger.

## Next Step

v11.2 removes the technical processed-COV blocker for GSE286302 but does not
solve the data gate. The next scientifically useful step remains a v12/v11.3
data acquisition decision: find or produce exact-age, bulk, same-tissue old
mouse methylome data for brain_cortex/heart/lung, or formally declare the
public-data route failed and move to a minimal ETL/design-RFC path.

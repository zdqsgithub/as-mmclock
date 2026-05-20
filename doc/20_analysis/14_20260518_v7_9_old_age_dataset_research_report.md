# v7.9 Old-Age Dataset Research and GEO Preflight Report

Date: 2026-05-18

## Summary

v7.9 executed the research/preflight step recommended by v7.8. The goal was not to tune another model, but to find official GEO processed methylation datasets that can repair the current age-range blocker: GSE121141 and GSE80672 old-age samples exceed same-tissue training support, so the current v7.4/v7.5 clocks are not full-lifespan cross-dataset clocks.

The main result is clear: `GSE213628` is the best next dataset to add. It is mouse RRBS/Bisulfite-seq, has 120 samples, multiple tissues, exact ages through 24 months plus heart samples up to 28 months, and an official processed `GSE213628_RAW.tar` of COV files. It directly targets the old-age multi-tissue gap found in v7.8.

Do not download a broad set yet. The next execution step should be a narrow v8: download and convert `GSE213628` first, then rerun GroupKFold/LODO and GSE121141 old-age held-out diagnostics. P2 datasets should wait for parser smoke tests and should not become the primary benchmark until GSE213628 proves useful.

## Official Method

Official sources used:

- GEO download documentation: https://www.ncbi.nlm.nih.gov/geo/info/download.html
- GEO programmatic access/E-Utils documentation: https://www.ncbi.nlm.nih.gov/geo/info/geo_paccess.html
- GEO SOFT format documentation: https://www.ncbi.nlm.nih.gov/geo/info/soft.html
- GSE213628 official GEO page: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213628
- GSE224442 official GEO page: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224442
- GSE92486 official GEO page: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92486
- GSE233879 official GEO page: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233879
- GSE295059 official GEO page: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE295059

Implemented script:

- `scripts/etl/16_geo_old_age_candidate_discovery.py`

Generated artifacts:

- `metadata/geo_old_age_candidate_inventory.csv`
- `metadata/geo_old_age_candidate_inventory.json`
- `metadata/geo_old_age_candidate_supplements.csv`
- `metadata/geo_old_age_candidate_samples.csv`
- `metadata/geo_candidate_soft/*.soft.gz`

The script uses only official E-Utils summaries, GEO family SOFT files, and GEO FTP `suppl/filelist.txt` plus `HEAD`. It does not download large `*_RAW.tar` archives.

## Candidate Ranking

| Dataset | Tier | Role | Samples | Max age months | Tissues | Processed schema | Size |
|---|---|---|---:|---:|---:|---|---:|
| GSE213628 | P1 | next primary old-age multi-tissue training data | 120 | 28.0 | 8 | COV tar | 3.26 GB |
| GSE224442 | P2 | parabiosis/recovery validation, liver+blood | 78 | 26.0 | 2 | COV tar | 2.43 GB |
| GSE92486 | P2 | DR/AL old liver validation, WGBS-like cov txt | 30 | 26.0 | 1 | COV/TXT tar | 3.06 GB |
| GSE175410 | P2 | exercise old skeletal muscle auxiliary validation | 15 | 24.0 | 1 | COV tar | 0.36 GB |
| GSE295059 | P2 | organoid aging, not primary in vivo benchmark | 70 | 24.0 | 1 | COV tar | 5.13 GB |
| GSE233879 | P2 | HSC/blood aging/rejuvenation, cell-type limited | 80 | 10.0 parsed | 2 | COV tar | 2.12 GB |
| GSE286302 | P3 | circadian intervention, metadata needs age resolver | 49 | not sample-specific | 3 | COV tar | 6.35 GB |
| GSE213723 | P3 | SuperSeries; use GSE213628 subseries instead | 144 | 28.0 | 8 | COV tar | 3.33 GB |

Notes:

- `GSE213628` should be treated as the next primary dataset, not `GSE213723`, because `GSE213723` is a SuperSeries and duplicates subseries context.
- `GSE224442` has valid old/young parabiosis/recovery structure, but it is an intervention dataset with only liver and blood; it is better as validation than as the first fix for multi-tissue old-age support.
- `GSE92486` is scientifically important for dietary restriction, but it is liver-only and WGBS/BS-seq rather than the current RRBS mainline. It needs adapter smoke testing before benchmark use.
- `GSE286302` is interesting but the SOFT records describe 6- and 15-month mice at the protocol level, not a clean sample-specific age field. It should not enter the benchmark until a conservative age resolver is implemented.

## Decision

Proceed with a narrow v8:

1. Extend supplement inventory/download for `GSE213628` only.
2. Download `GSE213628_RAW.tar` with resume and structured `metadata/download_log.jsonl`.
3. Run a parser smoke on 2-3 COV members before full conversion.
4. Convert to 5kb region matrix with the existing Bismark COV converter path.
5. Build a v8 multidataset matrix including GSE120137, GSE80672, GSE93957, GSE121141, GSE60012, and GSE213628.
6. Rerun the fixed v7.5/v7.8 diagnostic configuration:
   - GroupKFold by `dataset_batch`
   - LODO for every dataset
   - `all_except:GSE121141 -> GSE121141`
   - `all_except:GSE80672 -> GSE80672` CR residual validation
   - random-label and shuffled-intervention sanity checks

Do not expand deep learning or embedding search until GSE213628 shows whether the old-age same-tissue support blocker is reduced.

## Acceptance Criteria for v8

- `GSE213628` matrix builds with metadata-mapped samples >= 110.
- Common 5kb regions with current v7.1/v7.5 matrix are >= 50,000.
- GSE121141 held-out 104w+ MAE improves by at least 10 weeks or the report shows that the remaining failure is not age-range support.
- GroupKFold MAE improves by at least 2 weeks over v7.5 best (`23.767w`) without random-label leakage.
- GSE80672 CR held-out metrics are recalculated only from real predictions; no dummy AUC.

If GSE213628 does not improve old-age support, the next priority is not deep learning. The next priority is schema harmonization and careful inclusion of P2 intervention/single-tissue datasets as validation-only cohorts.

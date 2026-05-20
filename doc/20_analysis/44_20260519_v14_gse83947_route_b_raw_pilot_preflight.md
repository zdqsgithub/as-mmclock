# v14 GSE83947 Route B Raw Pilot Preflight

Date: 2026-05-19T03:44:25.773624+00:00

## Summary

This preflight mapped selected old-lung `GSE83947` GSM samples to official GEO SOFT, BioSample, SRA experiment, ENA run and FASTQ metadata. It did not download FASTQ, run Bismark, train models, or start autoresearch.

- Status: `ready_for_explicit_download_and_bismark_approval`
- Decision: `pending_explicit_download_and_compute_approval`
- Reason: metadata_and_environment_gates_passed

## Selected Runs

| sample_id | biosample | srx | run_accession | age_weeks | tissue | fastq_total_bytes | library_strategy | library_layout | download_authorized |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSM2223569 | SAMN05341816 | SRX1893806 | SRR3738487 | 108.643 | lung | 408118908 | Bisulfite-Seq | PAIRED | False |
| GSM2223613 | SAMN05341784 | SRX1893850 | SRR3738531 | 108.643 | lung | 366111345 | Bisulfite-Seq | PAIRED | False |
| GSM2223583 | SAMN05341848 | SRX1893820 | SRR3738501 | 108.643 | lung | 293772100 | Bisulfite-Seq | PAIRED | False |

## Environment Gate

- bismark: `/home/zdq-as/mouse_methyl_work/tools/route_b_bismark_env/bin/bismark`
- bowtie2: `/home/zdq-as/mouse_methyl_work/tools/route_b_bismark_env/bin/bowtie2`
- samtools: `/home/zdq-as/mouse_methyl_work/tools/route_b_bismark_env/bin/samtools`
- Bismark index candidates: `/home/zdq-as/mouse_methyl_work/references/GRCm38_ensembl102/Bisulfite_Genome; /home/zdq-as/mouse_methyl_work/references/GRCm38_ensembl102/Bisulfite_Genome/CT_conversion; /home/zdq-as/mouse_methyl_work/references/GRCm38_ensembl102/Bisulfite_Genome/GA_conversion`
- free disk GB: `4404.61`

## Guardrails

- raw FASTQ download authorized: `false`
- Bismark authorized: `false`
- training authorized: `false`
- autoresearch authorized: `false`

## Next Action

Environment provisioning was completed after this preflight, followed by an
explicitly approved minimal FASTQ/Bismark pilot for the three listed runs. The
pilot aligned all three samples and produced valid methylation coverage, but the
5kb region matrix had only `1,814` common regions with the v8.2 reference, below
the `>=50,000` gate.

Final Route B decision for this candidate:
`matrix_gate_failed_auxiliary_only`.

This does not authorize model training, RALPH Learn/Benchmark, or autoresearch.
`GSE83947` remains old-lung auxiliary evidence only.

Follow-up outputs:

- `results/ralph_v14_public_data_rescue/route_b_environment/route_b_bismark_environment_manifest.json`
- `results/ralph_v14_public_data_rescue/gse83947_raw_fastq_pilot/raw_fastq_download_state.json`
- `results/ralph_v14_public_data_rescue/gse83947_bismark_pilot/gse83947_bismark_pilot_state.json`
- `doc/20_analysis/43_20260519_v14_public_data_rescue_report.md`

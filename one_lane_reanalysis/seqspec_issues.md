# Seqspec Issues Summary

Generated from the one-lane reanalysis outputs in `one_lane_reanalysis/`.

## Inputs used

- `charles_gemx_v3`: `charles_gemx_v3/sample_metadata_gcp_2026_02_15.csv`, group `2_1`
- `charles_htv2`: `charles_htv2/sample_metadata_gcp_2026_02_15.csv`, group `2_1`
- `gary`: `gary_sample_metadata_gcp_2026_02_26.csv`, group `1_sample`
- `huangfu`: `/Users/lf588/Downloads/crispr_validator/huangfu_sample_metadata_gcp_2026_02_15.csv`, group `1_5`

## Observed seqspec-side problems

### charles_gemx_v3

- `gRNA guide`: seqspec says `R2 reverse 63-82` while the predictor found `R2 reverse 63-81`.
- The guide metadata file is all 19 bp spacers, so the remaining 1 bp mismatch is consistent with a 20 bp guide interval in seqspec versus 19 bp guide references in metadata.
- `scRNA`: `seqspec index` only exposed the RNA interval in this run, so seqspec-side `barcode` and `umi` were missing from the comparison.

### charles_htv2

- `gRNA guide`: seqspec says `R2 reverse 63-82` while the predictor found `R2 reverse 63-81`.
- `gRNA umi`: seqspec says `R1 forward 16-27` while the predictor found `R1 forward 16-25`.
- The guide metadata file is all 19 bp spacers, so the guide mismatch is again consistent with a 20 bp guide interval in seqspec versus 19 bp guide references in metadata.
- `scRNA`: `seqspec index` only exposed the RNA interval in this run, so seqspec-side `barcode` and `umi` were missing from the comparison.

### gary

- `seqspec index` returned no indexed regions for the `gRNA`, `scRNA`, or `hash` seqspec files in this run.
- Because of that, every region in the comparison is flagged as `missing_seqspec`.
- This means the current Gary seqspec files are present and readable, but they are not yielding tabular intervals through the `seqspec index` call used by this pipeline.

### huangfu

- `seqspec index` returned no indexed regions for both the `gRNA` and `scRNA` seqspec files in this run.
- Because of that, every region in the comparison is flagged as `missing_seqspec`.
- This is separate from the parser bug that was fixed earlier for Huangfu file ids; after that fix, the run completed successfully, but the seqspec files still did not emit indexed regions in the comparison step.

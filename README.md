



# Seqspec Validator CLI

## Install

Install directly from GitHub:

```bash
pip install "git+https://github.com/LucasSilvaFerreira/crispr_validator.git"
```

After install, the console scripts are:

- `seqspec-parser`
- `seqspec-check`
- `downloading-from-samplesheet`

Example:

```bash
seqspec-parser igvf \
  --accession IGVFDS9445RJOU \
  --analysis-root igvf_analysis \
  --igvf-keypair igvf_key.json \
  --one-lane
```

Python dependencies are installed by `pip`, but a few runtime tools are still external:

- `seqspec` must be on `PATH` for comparison runs
- `gsutil` is still required when samplesheets reference `gs://` assets

`seqspec_parser.py` compares `seqspec index` output against the interval predictor from `seqspec_check.py` and writes:

- `seqspec_report.html`
- `prediction_report.html`
- `comparison.html`
- `analysis_summary.json`

It supports two execution modes.

# Quick DACC Guide to Check Already Deposited Samples

Use this workflow to check samples that have already been deposited in DACC. (for not samplesheet check the mode3: single-lane)



## Download an IGVF portal analysis set and run validation

The command below downloads an IGVF portal analysis set, including the sample sheet and related files, and then runs validation after the download completes.

```bash
python3 seqspec_parser.py igvf \
  --accession IGVFDS9445RJOU \
  --analysis-root igvf_analysis \
  --igvf-keypair igvf_key.json \
  --one-lane
```

This script downloads:

- Seqspec files
- FASTQ sample information
- Metadata, including hash and guide information

## Credentials format

Provide your credentials in a JSON file like this:

```json
{
  "Key": "xxxxxx",
  "secret": "xxxxxxxxx"
}
```

## Validation output

The validation step returns a JSON file.

In that output, you can find the predicted flag for each modality, similar to what is reported when using SeqSpec indexing.

Example fields include:

```json
"technology": "0,0,16:0,16,26:1,20,41"

or 

"flag": "perfect_match"
```

Where:

- `technology` describes the predicted region structure in the cb/umi/feature format (kb)
- `flag` indicates whether the seqspec and the predicted regions match

## HTML reports

The workflow also generates HTML reports in the same output directory.

These reports can be helpful for spotting differences between the seqspec definition and the predicted regions.





## Mode 1: Samplesheet

Use this when you have a metadata table with `R1_path`, `R2_path`, `file_modality`, `seqspec`, `barcode_onlist`, `guide_design`, and optional `barcode_hashtag_map`.

Example:

```bash
python3 seqspec_parser.py samplesheet \
  --samplesheet charles_htv2/sample_metadata_gcp_2026_02_15.csv \
  --analysis-root on_lane_more_reads/charles_htv2 \
  --filter sequencing_run=2 \
  --filter lane=1
```

This groups rows by `sequencing_run,lane` by default and runs one comparison per complete modality set.

The samplesheet can reference:

- local files
- `gs://` objects
- authenticated IGVF portal download URLs (`https://api.data.igvf.org/...`)

If the sheet contains IGVF portal URLs, pass credentials with `--igvf-keypair` or `--igvf-api-key` plus `--igvf-api-secret` (or set `IGVF_API_KEY` and `IGVF_SECRET_KEY`).

## Mode 2: IGVF Analysis Set

Use this when you want the validator to generate the samplesheet directly from an IGVF analysis-set accession and then optionally download the needed reads/reference files through the IGVF API.

If `vendor/CRISPR_Pipeline/download_development/generate_per_sample.py` is present, `igvf` mode uses that upstream script first and falls back to the local generator when the upstream script is absent or fails.

Write only the samplesheet:

```bash
python3 seqspec_parser.py igvf \
  --accession IGVFDS9445RJOU \
  --analysis-root igvf_analysis \
  --samplesheet-only \
  --igvf-keypair igvf_key.json
```

Generate the samplesheet and continue into validation:

```bash
python3 seqspec_parser.py igvf \
  --accession IGVFDS9445RJOU \
  --analysis-root igvf_analysis \
  --igvf-keypair igvf_key.json
```

Generate the upstream samplesheet and stage only minimal read chunks plus reference files:

```bash
python3 seqspec_parser.py igvf \
  --accession IGVFDS9445RJOU \
  --analysis-root igvf_analysis \
  --igvf-keypair igvf_key.json \
  --download-only
```

Generate the upstream samplesheet but process only one representative lane:

```bash
python3 seqspec_parser.py igvf \
  --accession IGVFDS9445RJOU \
  --analysis-root igvf_analysis \
  --igvf-keypair igvf_key.json \
  --one-lane
```

Generate the upstream samplesheet, stage one representative lane, and stop before validation:

```bash
python3 seqspec_parser.py igvf \
  --accession IGVFDS9445RJOU \
  --analysis-root igvf_analysis \
  --igvf-keypair igvf_key.json \
  --one-lane \
  --download-only
```

The generated samplesheet is validator-compatible and includes:

- `R1_path` / `R2_path` as portal download URLs
- `R1_md5sum` / `R2_md5sum`
- portal-backed `seqspec`, barcode whitelist, guide metadata, and optional hash metadata URLs

`--chunk-bytes` still applies in IGVF mode, so FASTQ downloads can be staged as byte-range chunks for quick validation runs, or downloaded in full with `--chunk-bytes 0`.
`--one-lane` is optional and only affects the processing step after the samplesheet is generated; by default, `igvf` mode still processes every lane.

## Mode 3: Single Lane

Use this when you are outside the samplesheet context and want to provide one lane directly from explicit files.

Required inputs:

- scRNA `R1`, `R2`, and seqspec
- guide `R1`, `R2`, and seqspec
- barcode whitelist
- guide metadata

Optional hash inputs:

- hash `R1`, `R2`, seqspec, and hash metadata

Example without hash:

```bash
python3 seqspec_parser.py single-lane \
  --analysis-root single_lane_demo \
  --group-label charles_htv2_lane1 \
  --scrna-r1 on_lane_more_reads/charles_htv2/2_1/downloads/scRNA_R1_10mb.fastq.gz \
  --scrna-r2 on_lane_more_reads/charles_htv2/2_1/downloads/scRNA_R2_10mb.fastq.gz \
  --scrna-seqspec on_lane_more_reads/charles_htv2/2_1/seqspec/scRNA_seqspec.yaml \
  --guide-r1 on_lane_more_reads/charles_htv2/2_1/downloads/gRNA_R1_10mb.fastq.gz \
  --guide-r2 on_lane_more_reads/charles_htv2/2_1/downloads/gRNA_R2_10mb.fastq.gz \
  --guide-seqspec on_lane_more_reads/charles_htv2/2_1/seqspec/gRNA_seqspec.yaml \
  --barcode-whitelist on_lane_more_reads/charles_htv2/2_1/downloads/scRNA_barcodeWhitelist.tsv \
  --guide-metadata on_lane_more_reads/charles_htv2/2_1/downloads/gRNA_guideMetadata.tsv
```

Example with hash:

```bash
python3 seqspec_parser.py single-lane \
  --analysis-root single_lane_demo_gary \
  --group-label gary_lane \
  --scrna-r1 on_lane_more_reads/gary/1_sample/downloads/scRNA_R1_10mb.fastq \
  --scrna-r2 on_lane_more_reads/gary/1_sample/downloads/scRNA_R2_10mb.fastq \
  --scrna-seqspec on_lane_more_reads/gary/1_sample/seqspec/scRNA_seqspec.yaml \
  --guide-r1 on_lane_more_reads/gary/1_sample/downloads/gRNA_R1_10mb.fastq \
  --guide-r2 on_lane_more_reads/gary/1_sample/downloads/gRNA_R2_10mb.fastq \
  --guide-seqspec on_lane_more_reads/gary/1_sample/seqspec/gRNA_seqspec.yaml \
  --hash-r1 on_lane_more_reads/gary/1_sample/downloads/hash_R1_10mb.fastq \
  --hash-r2 on_lane_more_reads/gary/1_sample/downloads/hash_R2_10mb.fastq \
  --hash-seqspec on_lane_more_reads/gary/1_sample/seqspec/hash_seqspec.yaml \
  --hash-metadata on_lane_more_reads/gary/1_sample/downloads/hash_hashMetadata.tsv \
  --barcode-whitelist on_lane_more_reads/gary/1_sample/downloads/scRNA_barcodeWhitelist.tsv \
  --guide-metadata on_lane_more_reads/gary/1_sample/downloads/gRNA_guideMetadata.tsv
```

This command has already been run in this repository. The resulting example reports are in:

- `single_lane_demo_gary/index.html`
- `single_lane_demo_gary/gary_lane/comparison.html`
- `single_lane_demo_gary/gary_lane/seqspec_report.html`
- `single_lane_demo_gary/gary_lane/prediction_report.html`

## How The Predictor Works

The predictor is deliberately simple and seqspec-oriented. It tries to recover one shared barcode and one shared UMI for the lane, then adds modality-specific RNA, guide, and optional hash intervals.

### Barcode

- A whitelist is loaded from `barcode_whitelist`, and the dominant barcode length in that file becomes the barcode length used for the scan.
- Only a barcode subsample is used for placement. By default that is the first `10,000` reads from the preferred source modality.
- The source modality is chosen in priority order `rna -> guide -> hash`, unless `--barcode-source` forces one.
- For each candidate configuration, the code scans `R1` and `R2`, on both `forward` and `reverse`, and counts exact whitelist hits at every possible start position.
- The winning barcode configuration is the region with the highest score:

```text
barcode_score = 4 * dominant_read_fraction + 2 * position_purity + hit_ratio
```

- `hit_ratio` is the fraction of sampled reads that contain at least one whitelist barcode hit.
- `position_purity` is the fraction of all whitelist hits that land at the dominant start position.
- `dominant_read_fraction` is the fraction of all sampled reads supporting the dominant start position.
- The barcode is considered well supported only when hits cluster strongly at one interval. The predictor does not assume `R1` a priori.

### UMI

- The UMI is not searched globally. It is searched only after the barcode interval is fixed.
- Reads are reoriented into the winning barcode strand, and only barcode-supporting reads are kept for the UMI step.
- Barcode placement still uses the first `10,000` reads by default, but UMI inference uses the full sampled source pool already loaded for that modality. With current defaults that means up to `100,000` reads.
- The code only tests the two windows immediately after the barcode:
  - `barcode|UMI(10)`
  - `barcode|UMI(12)`
- Each candidate window is scored from the extracted strings:

```text
umi_score =
  1.5 * unique_fraction
  + 1.1 * (mean_entropy / 2)
  - 0.8 * top_sequence_fraction
  - 0.5 * n_fraction
```

- `unique_fraction` is the fraction of unique UMI strings among the barcode-supporting reads.
- `mean_entropy` is the average per-base Shannon entropy across the window.
- `top_sequence_fraction` penalizes windows dominated by one repeated sequence.
- `n_fraction` penalizes windows containing ambiguous bases.
- After scoring, the code explicitly compares the `10 bp` and `12 bp` windows anchored immediately after the barcode.
- A `12 bp` UMI is accepted only when the extra `2 bp` are genuinely informative:
  - they add enough entropy in the 2 bp tail
  - they improve uniqueness or collision reduction enough relative to the 10 bp anchor
- If the extra 2 bases add little evidence, the shorter `10 bp` interval is kept. This is why increasing the number of reads can change a borderline `10 vs 12` decision.

### RNA

- RNA is not sequence-matched in the current predictor.
- Once the shared barcode orientation is known, the RNA region is placed on the full opposite read and opposite strand.
- Example: if the barcode is `R1 forward`, RNA is emitted as the full `R2 reverse`.

### Guide And Hash

- Guide and hash inference use the same exact-match strategy.
- Reference sequences come from:
  - `guide_metadata` for guide
  - `hash_metadata` for hash
- The code loads the feature sequences, deduplicates them, and scans `R1` and `R2` on both strands.
- Matching is exact sequence search against the metadata strings and their reverse complements.
- For each read/strand configuration, the winning region is the dominant exact-match interval scored as:

```text
feature_score = 3 * hit_ratio + 2 * position_purity + flank_purity + (1 - gini)
```

- `hit_ratio` is the fraction of sampled reads with any exact feature hit.
- `position_purity` is the concentration of hits at the dominant start position.
- `flank_purity` is the consistency of the 12 bp immediately upstream of the dominant interval. This stabilizes calls like guides that tend to occur in one structural context.
- `gini` measures how uneven the hits are across reference sequences. Lower gini means the matched reference set is less dominated by one sequence and is therefore more compatible with a real library rather than one artifact sequence.
- For guides, the report can also emit the `12 bp` immediately before the winning guide interval as `guide_prefix`.
- Hash uses the hash metadata table from the samplesheet or `--hash-metadata` input. It does not reuse the guide metadata.

### Shared Design Assumption

- The barcode and UMI are inferred once per lane from the chosen source modality and then reused for RNA, guide, and optional hash.
- This follows the assumption that the modalities in the lane share the same bead or cell barcode design.

### Comparison Against Seqspec

- `seqspec index` is run separately for each modality seqspec, and the extracted seqspec intervals are compared with the predictor intervals.
- Interval comparison uses read identity and coordinate distance.
- If seqspec and prediction are on different reads, the comparison is still flagged as wrong.
- If they are on the same read but opposite strands, the interval flag ignores the strand mismatch and a separate strand-status column records that difference.
- Flags are:
  - `perfect_match`: same read and exact interval
  - `close_enough`: same read and at most `2 bp` away
  - `very_distant`: larger coordinate disagreement or different reads
  - `missing_seqspec`: seqspec did not yield a matching region for that modality
  - `missing_prediction`: predictor did not emit that region

## Practical Reading Of The Reports

- `seqspec_report.html` shows what `seqspec index` extracted from the seqspec YAML and which indexed region was used for comparison.
- `prediction_report.html` shows the inferred intervals, barcode position histogram, UMI candidates, and guide/hash exact-match diagnostics.
- `comparison.html` is the direct seqspec-vs-prediction table with interval and strand flags.

## Notes

- `samplesheet` mode still supports `--downloads-dir` when filters reduce the run to one group.
- `igvf` mode writes a local TSV first, then reuses the same grouped validator flow as `samplesheet` mode.
- `--download-only` stages chunked FASTQs and the full seqspec / whitelist / metadata files without running seqspec or predictor reports.
- `single-lane` mode treats the inputs as one isolated lane and does not try to infer multi-lane grouping.
- Local paths, `gs://` paths, and authenticated IGVF portal URLs are supported.
- If hash is used in `single-lane` mode, all four hash flags must be provided together.

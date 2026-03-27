# seqspec index Commands And Errors

Generated from the one-lane reanalysis outputs.

## charles_gemx_v3 / 2_1

### gRNA

#### tab

Command:
```bash
seqspec index -m crispr -s file -i dummy_guide_R1.fq,dummy_guide_R2.fq one_lane_reanalysis/charles_gemx_v3/2_1/seqspec/gRNA_seqspec.yaml
```

Exit code: `0`

stdout:
```text
dummy_guide_R1.fq	barcode	barcode	0	16
dummy_guide_R1.fq	umi	umi	16	26
dummy_guide_R2.fq	feature_primer	custom_primer	0	26
dummy_guide_R2.fq	sgrna_scaffold	linker	26	63
dummy_guide_R2.fq	guide	sgrna_target	63	83
dummy_guide_R2.fq	TSO	linker	83	90

```

stderr:
```text
2026-03-19 17:58:00.467 system_profiler[44140:41626615] XType: Using static font registry.
Matplotlib is building the font cache; this may take a moment.

```

#### kb

Command:
```bash
seqspec index -m crispr -s file -i dummy_guide_R1.fq,dummy_guide_R2.fq -t kb one_lane_reanalysis/charles_gemx_v3/2_1/seqspec/gRNA_seqspec.yaml
```

Exit code: `0`

stdout:
```text
0,0,16:0,16,26:1,63,83

```

stderr:
```text
<empty>
```

### scRNA

#### tab

Command:
```bash
seqspec index -m rna -s file -i dummy_rna_R1,dummy_rna_R2.fq one_lane_reanalysis/charles_gemx_v3/2_1/seqspec/scRNA_seqspec.yaml
```

Exit code: `0`

stdout:
```text
dummy_rna_R2.fq	rna	cdna	0	90

```

stderr:
```text
<empty>
```

#### kb

Command:
```bash
seqspec index -m rna -s file -i dummy_rna_R1,dummy_rna_R2.fq -t kb one_lane_reanalysis/charles_gemx_v3/2_1/seqspec/scRNA_seqspec.yaml
```

Exit code: `0`

stdout:
```text
-1,-1,-1:-1,-1,-1:0,0,90

```

stderr:
```text
<empty>
```


## charles_htv2 / 2_1

### gRNA

#### tab

Command:
```bash
seqspec index -m crispr -s file -i dummy_guide_R1.fq,dummy_guide_R2.fq one_lane_reanalysis/charles_htv2/2_1/seqspec/gRNA_seqspec.yaml
```

Exit code: `0`

stdout:
```text
dummy_guide_R1.fq	barcode	barcode	0	16
dummy_guide_R1.fq	umi	umi	16	28
dummy_guide_R2.fq	feature_primer	custom_primer	0	26
dummy_guide_R2.fq	sgrna_scaffold	linker	26	63
dummy_guide_R2.fq	guide	sgrna_target	63	83
dummy_guide_R2.fq	TSO	linker	83	90

```

stderr:
```text
<empty>
```

#### kb

Command:
```bash
seqspec index -m crispr -s file -i dummy_guide_R1.fq,dummy_guide_R2.fq -t kb one_lane_reanalysis/charles_htv2/2_1/seqspec/gRNA_seqspec.yaml
```

Exit code: `0`

stdout:
```text
0,0,16:0,16,28:1,63,83

```

stderr:
```text
<empty>
```

### scRNA

#### tab

Command:
```bash
seqspec index -m rna -s file -i dummy_rna_R1,dummy_rna_R2.fq one_lane_reanalysis/charles_htv2/2_1/seqspec/scRNA_seqspec.yaml
```

Exit code: `0`

stdout:
```text
dummy_rna_R2.fq	rna	cdna	0	90

```

stderr:
```text
<empty>
```

#### kb

Command:
```bash
seqspec index -m rna -s file -i dummy_rna_R1,dummy_rna_R2.fq -t kb one_lane_reanalysis/charles_htv2/2_1/seqspec/scRNA_seqspec.yaml
```

Exit code: `0`

stdout:
```text
-1,-1,-1:-1,-1,-1:0,0,90

```

stderr:
```text
<empty>
```


## gary / 1_sample

### gRNA

#### tab

Command:
```bash
seqspec index -m crispr -s file -i dummy_guide_R1,dummy_guide_R2 one_lane_reanalysis/gary/1_sample/seqspec/gRNA_seqspec.yaml
```

Exit code: `0`

stdout:
```text


```

stderr:
```text
<empty>
```

#### kb

Command:
```bash
seqspec index -m crispr -s file -i dummy_guide_R1,dummy_guide_R2 -t kb one_lane_reanalysis/gary/1_sample/seqspec/gRNA_seqspec.yaml
```

Exit code: `0`

stdout:
```text
-1,-1,-1:-1,-1,-1:

```

stderr:
```text
<empty>
```

### hash

#### tab

Command:
```bash
seqspec index -m tag -s file -i dummy_hash_R1,dummy_hash_R2 one_lane_reanalysis/gary/1_sample/seqspec/hash_seqspec.yaml
```

Exit code: `0`

stdout:
```text


```

stderr:
```text
<empty>
```

#### kb

Command:
```bash
seqspec index -m tag -s file -i dummy_hash_R1,dummy_hash_R2 -t kb one_lane_reanalysis/gary/1_sample/seqspec/hash_seqspec.yaml
```

Exit code: `0`

stdout:
```text
-1,-1,-1:-1,-1,-1:

```

stderr:
```text
<empty>
```

### scRNA

#### tab

Command:
```bash
seqspec index -m rna -s file -i dummy_rna_R1,dummy_rna_R2 one_lane_reanalysis/gary/1_sample/seqspec/scRNA_seqspec.yaml
```

Exit code: `0`

stdout:
```text


```

stderr:
```text
<empty>
```

#### kb

Command:
```bash
seqspec index -m rna -s file -i dummy_rna_R1,dummy_rna_R2 -t kb one_lane_reanalysis/gary/1_sample/seqspec/scRNA_seqspec.yaml
```

Exit code: `0`

stdout:
```text
-1,-1,-1:-1,-1,-1:

```

stderr:
```text
<empty>
```


## huangfu / 1_5

### gRNA

#### tab

Command:
```bash
seqspec index -m crispr -s file -i IGVFFI7570VNWD,IGVFFI4496IORJ,IGVFFI9535CMTO,IGVFFI3874BDJV one_lane_reanalysis/huangfu/1_5/seqspec/gRNA_seqspec.yaml
```

Exit code: `0`

stdout:
```text


```

stderr:
```text
<empty>
```

#### kb

Command:
```bash
seqspec index -m crispr -s file -i IGVFFI7570VNWD,IGVFFI4496IORJ,IGVFFI9535CMTO,IGVFFI3874BDJV -t kb one_lane_reanalysis/huangfu/1_5/seqspec/gRNA_seqspec.yaml
```

Exit code: `0`

stdout:
```text
-1,-1,-1:-1,-1,-1:

```

stderr:
```text
<empty>
```

### scRNA

#### tab

Command:
```bash
seqspec index -m rna -s file -i IGVFFI6630MDRJ,IGVFFI8015BUKO,IGVFFI3359QDOY,IGVFFI0934OGKX one_lane_reanalysis/huangfu/1_5/seqspec/scRNA_seqspec.yaml
```

Exit code: `0`

stdout:
```text


```

stderr:
```text
<empty>
```

#### kb

Command:
```bash
seqspec index -m rna -s file -i IGVFFI6630MDRJ,IGVFFI8015BUKO,IGVFFI3359QDOY,IGVFFI0934OGKX -t kb one_lane_reanalysis/huangfu/1_5/seqspec/scRNA_seqspec.yaml
```

Exit code: `0`

stdout:
```text
-1,-1,-1:-1,-1,-1:

```

stderr:
```text
<empty>
```


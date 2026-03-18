# fgumi TAPs/Illumina 5-Base Pipeline Guide

This guide describes how to process TAPs (TET-assisted pyridine borane sequencing) and Illumina 5-base chemistry data through fgumi's consensus pipeline. It covers both simplex and duplex consensus calling workflows.

## Background

TAPs converts **methylated** cytosines (5mC and 5hmC) to thymine through TET oxidation followed by pyridine borane reduction. This is the opposite of EM-Seq:

- **Methylated C** (5mC/5hmC): converted to T by TET + pyridine borane
- **Unmethylated C**: remains as C in sequencing reads

| | EM-Seq | TAPs |
|---|---|---|
| **What gets converted** | Unmethylated C → T | Methylated C → T |
| **C in read at ref-C** | Methylated (protected) | Unmethylated (not a target) |
| **T in read at ref-C** | Unmethylated (converted) | Methylated (converted) |

### Impact on UMI Processing

TAPs does **not** require special UMI handling — UMI sequences contain unmethylated cytosines (synthetic oligonucleotides), and since TAPs only converts *methylated* C, UMIs are unaffected.

At the consensus calling stage, reads showing T at a reference C position represent *methylated* cytosines (opposite of EM-Seq). The `--methylation-mode taps` flag tells fgumi to interpret the MM/ML methylation probability accordingly.

---

## Pipeline Overview

The TAPs pipeline follows the same structure as the [standard consensus pipeline](best-practice-consensus-pipeline.md), with `--methylation-mode taps` at the consensus and filter steps:

```
Phase 1: FASTQ → Grouped BAM
  extract → [correct] → fastq | bwa mem | zipper → sort → group

Phase 2: Grouped BAM → Filtered Consensus
  simplex/duplex → fastq | bwa mem | zipper → filter → sort
```

Steps with TAPs-specific flags:

| Step | Flag | Purpose |
|------|------|---------|
| `simplex` | `--methylation-mode taps --ref` | Methylation-aware consensus (TAPs interpretation) |
| `duplex` | `--methylation-mode taps --ref` | Methylation-aware consensus (TAPs interpretation) |
| `filter` | `--methylation-mode taps`, `--min-conversion-fraction` | TAPs-specific conversion fraction check |

---

## Workflow A: Random UMIs (No Fixed UMI Set)

### Step 1: UMI Extraction

Extract UMIs from FASTQ. No TAPs-specific flags needed.

**Simplex (single UMI per read pair):**

```bash
fgumi extract \
  --inputs r1.fq.gz r2.fq.gz \
  --read-structures 8M+T +T \
  --sample "sample_name" \
  --library "library_name" \
  --output unmapped.bam \
  --threads 4
```

**Duplex (UMI from both ends):**

```bash
fgumi extract \
  --inputs r1.fq.gz r2.fq.gz \
  --read-structures 8M+T 8M+T \
  --sample "sample_name" \
  --library "library_name" \
  --output unmapped.bam \
  --threads 4
```

### Step 2: Alignment

```bash
fgumi fastq --input unmapped.bam \
  | bwa mem -t 16 -p -K 150000000 -Y ref.fa - \
  | fgumi zipper --unmapped unmapped.bam --reference ref.fa --output aligned.bam
```

### Step 3: Sort

```bash
fgumi sort \
  --input aligned.bam \
  --output sorted.bam \
  --order template-coordinate \
  --threads 8 \
  --max-memory 4G
```

### Step 4: UMI Grouping

Standard grouping — no TAPs-specific flags needed.

**Simplex:**

```bash
fgumi group \
  --input sorted.bam \
  --output grouped.bam \
  --strategy adjacency \
  --edits 1 \
  --family-size-histogram fam_sizes.txt \
  --threads 8
```

**Duplex:**

```bash
fgumi group \
  --input sorted.bam \
  --output grouped.bam \
  --strategy paired \
  --edits 1 \
  --family-size-histogram fam_sizes.txt \
  --threads 8
```

### Step 5: Consensus Calling

Use `--methylation-mode taps` and `--ref` to enable TAPs methylation-aware consensus.

**Simplex:**

```bash
fgumi simplex \
  --input grouped.bam \
  --output consensus.bam \
  --min-reads 1 \
  --min-input-base-quality 20 \
  --output-per-base-tags \
  --methylation-mode taps \
  --ref ref.fa \
  --threads 8
```

**Duplex:**

```bash
fgumi duplex \
  --input grouped.bam \
  --output consensus.bam \
  --min-reads 1 \
  --min-input-base-quality 20 \
  --output-per-base-tags \
  --methylation-mode taps \
  --ref ref.fa \
  --threads 8
```

### Step 6: Re-alignment

Consensus reads are unmapped and must be re-aligned:

```bash
fgumi fastq --input consensus.bam \
  | bwa mem -t 16 -p -K 150000000 -Y ref.fa - \
  | fgumi zipper --unmapped consensus.bam --reference ref.fa --output consensus.mapped.bam
```

### Step 7: Filtering

Use `--methylation-mode taps` with `--min-conversion-fraction` to check TAPs specificity. In TAPs, non-CpG cytosines should remain unconverted (they are unmethylated). A high unconverted fraction at non-CpG sites indicates good TAPs specificity.

**Simplex filtering:**

```bash
fgumi filter \
  --input consensus.mapped.bam \
  --output filtered.bam \
  --ref ref.fa \
  --min-reads 3 \
  --max-base-error-rate 0.1 \
  --max-no-call-fraction 0.2 \
  --min-methylation-depth 3 \
  --methylation-mode taps \
  --min-conversion-fraction 0.9 \
  --reverse-per-base-tags \
  --threads 8
```

**Duplex filtering:**

```bash
fgumi filter \
  --input consensus.mapped.bam \
  --output filtered.bam \
  --ref ref.fa \
  --min-reads 10,5,3 \
  --max-base-error-rate 0.1 \
  --max-no-call-fraction 0.2 \
  --min-methylation-depth 10,5,3 \
  --require-single-strand-agreement \
  --require-strand-methylation-agreement \
  --methylation-mode taps \
  --min-conversion-fraction 0.9 \
  --reverse-per-base-tags \
  --threads 8
```

### Step 8: Final Sort

```bash
fgumi sort \
  --input filtered.bam \
  --output final.bam \
  --order coordinate \
  --threads 8
```

---

## Workflow B: Fixed UMIs (Known UMI Set)

When UMIs come from a fixed set, add a correction step before alignment.

### Step 1: UMI Extraction

Same as Workflow A.

### Step 2: UMI Correction

Correct UMIs against the known whitelist:

```bash
fgumi correct \
  --input unmapped.bam \
  --output corrected.bam \
  --umi-files known_umis.txt \
  --max-mismatches 1 \
  --min-distance 1 \
  --metrics correction_metrics.txt \
  --threads 8
```

### Steps 3-8: Alignment through Final Sort

Same as Workflow A (steps 2-8).

---

## TAPs Output Tags

When `--methylation-mode taps` is enabled, consensus reads carry the same BAM tags as EM-Seq mode. The tag names and formats are identical — only the interpretation of MM/ML probabilities differs.

### Simplex Output Tags

| Tag | Type | Description |
|-----|------|-------------|
| `MM` | `Z` | SAM-spec methylation modification calls (sparse format) |
| `ML` | `B:C` | Methylation modification probabilities (companion to MM) |
| `cu` | `B:s` | Per-base unconverted count (reads showing C at ref-C) |
| `ct` | `B:s` | Per-base converted count (reads showing T at ref-C) |

### Duplex Output Tags

All simplex tags above (combined from both strands), plus per-strand tags:

| Tag | Type | Description |
|-----|------|-------------|
| `am` | `Z` | AB strand methylation calls (MM format, no ML companion) |
| `bm` | `Z` | BA strand methylation calls (MM format, no ML companion) |
| `au` | `B:s` | AB strand unconverted count |
| `at` | `B:s` | AB strand converted count |
| `bu` | `B:s` | BA strand unconverted count |
| `bt` | `B:s` | BA strand converted count |

### MM/ML Interpretation Differences

The `cu` and `ct` count tags have the same meaning in both chemistries:
- `cu`: reads showing C (unconverted) at a reference C position
- `ct`: reads showing T (converted) at a reference C position

The **MM/ML probability** differs:
- **EM-Seq**: `prob = cu / (cu + ct)` — higher probability = more methylated (C stayed as C because it was methylated)
- **TAPs**: `prob = ct / (cu + ct)` — higher probability = more methylated (C was converted to T because it was methylated)

The `MM`/`ML` tags follow the [SAM-spec methylation format](https://samtools.github.io/hts-specs/SAMtags.pdf) and are compatible with downstream methylation analysis tools.

---

## TAPs Filter Options

The `filter` command provides methylation-specific options that work with both EM-Seq and TAPs data. These operate on the `cu`/`ct`/`au`/`at`/`bu`/`bt` count tags.

### `--min-methylation-depth`

Identical to EM-Seq. Per-base masking based on methylation evidence depth. Bases where `cu[i] + ct[i]` is below the threshold are masked to N.

Accepts 1-3 comma-delimited values for duplex reads:

| Values | Meaning |
|--------|---------|
| `5` | 5 for all levels |
| `10,5` | 10 for duplex combined, 5 for each strand |
| `10,5,3` | 10 for duplex combined, 5 for AB strand, 3 for BA strand |

### `--require-strand-methylation-agreement`

Identical to EM-Seq. Duplex-only, per-base masking. Requires `--ref`.

### `--min-conversion-fraction` with `--methylation-mode`

Read-level filter. Requires `--ref`. Accepts a value between 0.0 and 1.0.

**With `--methylation-mode em-seq`:** Checks `ct / (cu + ct) >= threshold` at non-CpG ref-C positions. High conversion = good enzymatic conversion efficiency.

**With `--methylation-mode taps`:** Checks `cu / (cu + ct) >= threshold` at non-CpG ref-C positions. High non-conversion at non-CpG = good TAPs specificity (non-CpG Cs are unmethylated and should not be converted).

CpG positions are excluded from both calculations because they may have variable methylation status.

---

## Recommended Parameters

### TAPs Simplex (Moderate Stringency)

```bash
fgumi simplex --min-reads 1 --min-input-base-quality 20 --output-per-base-tags --methylation-mode taps --ref ref.fa
fgumi filter --ref ref.fa --min-reads 3 --max-base-error-rate 0.1 --min-methylation-depth 3 \
  --methylation-mode taps --min-conversion-fraction 0.9
```

### TAPs Duplex (High Specificity)

```bash
fgumi duplex --min-reads 1 --min-input-base-quality 20 --output-per-base-tags --methylation-mode taps --ref ref.fa
fgumi filter --ref ref.fa --min-reads 10,5,3 --max-base-error-rate 0.1 --min-methylation-depth 10,5,3 \
  --require-single-strand-agreement --require-strand-methylation-agreement \
  --methylation-mode taps --min-conversion-fraction 0.9
```

### TAPs Deduplication (No Consensus)

For workflows that mark duplicates without consensus calling:

```bash
fgumi dedup \
  --input sorted.bam \
  --output deduped.bam \
  --metrics metrics.txt
```

---

## Key Differences from EM-Seq Pipeline

| Aspect | EM-Seq | TAPs |
|--------|--------|------|
| Consensus calling flag | `--methylation-mode em-seq --ref` | `--methylation-mode taps --ref` |
| Filter flag | `--methylation-mode em-seq` | `--methylation-mode taps` |
| Conversion fraction meaning | ct/(cu+ct) >= threshold | cu/(cu+ct) >= threshold |

---

## Troubleshooting

### Missing MM/ML Tags on Output

Ensure both `--methylation-mode taps` and `--ref` are provided to the consensus caller. The reference FASTA must have an accompanying `.dict` file (generate with `samtools dict` if missing).

### Unexpected Masking from Strand Methylation Agreement

`--require-strand-methylation-agreement` only applies to duplex reads at CpG sites. If you see excessive masking:

1. Check that your library has adequate duplex coverage at CpG sites
2. Consider whether strand-specific methylation differences are biologically expected (e.g., imprinted regions)
3. This filter requires both strands to have evidence — positions with zero evidence on either strand are not masked

### Reads Filtered by Conversion Fraction

If many reads fail `--min-conversion-fraction` with `--methylation-mode taps`:

1. A high failure rate indicates that non-CpG cytosines are being converted, suggesting the TAPs chemistry is not sufficiently specific to methylated C
2. Try lowering the threshold (e.g., 0.8 instead of 0.9)
3. Reads with no non-CpG cytosine positions automatically pass this filter

### Accidentally Using EM-Seq Flags

If you use `--methylation-mode em-seq` instead of `--methylation-mode taps`, the methylation probabilities will be inverted — methylated positions will show low probability and vice versa. If downstream analysis shows unexpected methylation patterns, verify you used the correct mode for your chemistry.

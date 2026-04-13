# Methylation Pipeline Guide

This guide describes how to process methylation sequencing data through fgumi's consensus pipeline. It covers EM-Seq and TAPs/Illumina 5-base chemistries, for both simplex and duplex consensus calling workflows.

## Background

Both EM-Seq and TAPs detect cytosine methylation by converting one class of cytosines to thymine, but they target opposite classes:

| | EM-Seq | TAPs |
|---|---|---|
| **Chemistry** | TET2 + APOBEC | TET oxidation + pyridine borane |
| **What gets converted** | Unmethylated C → T | Methylated C → T |
| **C in read at ref-C** | Methylated (protected) | Unmethylated (not a target) |
| **T in read at ref-C** | Unmethylated (converted) | Methylated (converted) |

### Impact on UMI Processing

C→T conversion affects consensus calling: at a reference C position, reads showing T are not errors — they represent conversion events. Standard consensus calling would treat C/T disagreements as sequencing errors and penalize quality. Methylation mode recognizes these as conversion events and tracks per-base evidence through consensus calling.

**UMI sequences:**
- **EM-Seq**: UMIs should be synthesized with methylated cytosines (5mC) to protect them from enzymatic conversion. Unmethylated C in UMIs is a library prep issue.
- **TAPs**: UMIs are unaffected — synthetic oligonucleotides contain unmethylated cytosines, which TAPs does not convert.

---

## Pipeline Overview

The methylation pipeline follows the same structure as the [standard consensus pipeline](best-practices.md), with additional flags at the consensus, re-alignment, and filter steps. Methylation mode is supported by `simplex` and `duplex` consensus callers. The `codec` caller does not support methylation mode.

```text
Phase 1: FASTQ → Grouped BAM
  extract → [correct] → fastq | aligner | zipper → sort → group

Phase 2: Grouped BAM → Filtered Consensus
  simplex/duplex → fastq | aligner | zipper → filter → sort
```

### Chemistry-Specific Steps

| Step | EM-Seq | TAPs |
|------|--------|------|
| Alignment | `bwameth` (bisulfite-aware) | `bwa mem` (standard) |
| Consensus | `--methylation-mode em-seq --ref` | `--methylation-mode taps --ref` |
| Re-alignment zipper | `--restore-unconverted-bases` | (no additional flags) |
| Filter | `--methylation-mode em-seq` | `--methylation-mode taps` |

---

## Workflow A: Random UMIs (No Fixed UMI Set)

This is the simpler case. Random UMIs (e.g., random 8-mers ligated during library prep) do not need correction against a whitelist.

### Step 1: UMI Extraction

Extract UMIs from FASTQ. No methylation-specific flags needed here.

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

**EM-Seq** — use a bisulfite-aware aligner (bwameth) because unmethylated C→T conversion looks like bisulfite conversion:

```bash
fgumi fastq --input unmapped.bam --no-read-suffix \
  | bwameth.py --reference ref.fa --threads 16 --interleaved /dev/stdin \
  | samtools view -b \
  | fgumi zipper --unmapped unmapped.bam --reference ref.fa --output aligned.bam
```

**TAPs** — use a standard aligner (bwa mem) because only methylated Cs are converted, leaving most Cs intact:

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

Use `--methylation-mode` and `--ref` to enable methylation-aware consensus.

**Simplex:**

```bash
fgumi simplex \
  --input grouped.bam \
  --output consensus.bam \
  --min-reads 1 \
  --min-input-base-quality 20 \
  --output-per-base-tags \
  --methylation-mode <em-seq|taps> \
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
  --methylation-mode <em-seq|taps> \
  --ref ref.fa \
  --threads 8
```

### Step 6: Re-alignment

Consensus reads are unmapped and must be re-aligned.

**EM-Seq** — use `--restore-unconverted-bases` so that bases normalized during consensus (T→C at ref-C positions) are restored before bisulfite-aware re-alignment:

```bash
fgumi fastq --input consensus.bam --no-read-suffix \
  | bwameth.py --reference ref.fa --threads 16 --interleaved /dev/stdin \
  | samtools view -b \
  | fgumi zipper --unmapped consensus.bam --reference ref.fa --restore-unconverted-bases --output consensus.mapped.bam
```

**TAPs:**

```bash
fgumi fastq --input consensus.bam \
  | bwa mem -t 16 -p -K 150000000 -Y ref.fa - \
  | fgumi zipper --unmapped consensus.bam --reference ref.fa --output consensus.mapped.bam
```

### Step 7: Filtering

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
  --methylation-mode <em-seq|taps> \
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
  --methylation-mode <em-seq|taps> \
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

When UMIs come from a fixed set (e.g., a synthesized pool of known sequences), add a correction step before alignment. This maps observed UMIs back to the correct whitelist entry.

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

If your UMI design includes unmethylated cytosines, add `--allow-c-to-t`. This flag applies uniformly across all UMI segments regardless of read-pair index, since both R1 and R2 UMI segments are in forward orientation. Only C-to-T tolerance is needed; G-to-A tolerance is not required.

### Steps 3-8: Alignment through Final Sort

After correction, the remaining steps are the same as Workflow A (steps 2-8).

---

## Output Tags

When methylation mode is enabled, consensus reads carry additional BAM tags for methylation evidence.

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

### MM/ML Probability Interpretation

The `cu` and `ct` count tags have the same meaning in both chemistries:
- `cu`: reads showing C (unconverted) at a reference C position
- `ct`: reads showing T (converted) at a reference C position

The **MM/ML probability** differs:
- **EM-Seq**: `prob = cu / (cu + ct)` — higher probability = more methylated (C stayed as C because it was protected)
- **TAPs**: `prob = ct / (cu + ct)` — higher probability = more methylated (C was converted to T because it was methylated)

The `MM`/`ML` tags follow the [SAM-spec methylation format](https://samtools.github.io/hts-specs/SAMtags.pdf) and are compatible with downstream methylation analysis tools.

---

## Filter Options

The `filter` command provides methylation-specific options. These operate on the `cu`/`ct`/`au`/`at`/`bu`/`bt` count tags emitted by methylation-aware consensus calling.

### `--min-methylation-depth`

Per-base masking based on methylation evidence depth. Bases where `cu[i] + ct[i]` is below the threshold are masked to N.

Accepts 1-3 comma-delimited values for duplex reads, following the same convention as `--min-reads`:

| Values | Meaning |
|--------|---------|
| `5` | 5 for all levels |
| `10,5` | 10 for duplex combined, 5 for each strand |
| `10,5,3` | 10 for duplex combined, 5 for AB strand, 3 for BA strand |

For simplex reads, only the first value is used.

### `--require-strand-methylation-agreement`

Duplex-only, per-base masking. Requires `--ref`.

At each CpG dinucleotide in the reference, compares the methylation call from the top strand (AB: `au`/`at` at the C position) with the call from the bottom strand (BA: `bu`/`bt` at the G position). If one strand calls methylated and the other calls unmethylated, both positions of the CpG are masked to N.

This is analogous to `--require-single-strand-agreement` but specific to methylation status at CpG sites rather than raw base identity.

### `--min-conversion-fraction`

Read-level filter. Requires `--ref` and `--methylation-mode`. Accepts a value between 0.0 and 1.0.

Computes the conversion fraction at non-CpG reference cytosine positions across the read:

- **EM-Seq** (`--methylation-mode em-seq`): checks `ct / (cu + ct) >= threshold`. Non-CpG cytosines are expected to be unmethylated and therefore converted. High conversion = good enzymatic conversion efficiency.
- **TAPs** (`--methylation-mode taps`): checks `cu / (cu + ct) >= threshold`. Non-CpG cytosines are expected to be unmethylated and therefore *not* converted. High non-conversion at non-CpG = good TAPs specificity.

CpG positions are excluded from both calculations because they may have variable methylation status.

---

## Recommended Parameters

### Simplex (Moderate Stringency)

```bash
fgumi simplex --min-reads 1 --min-input-base-quality 20 --output-per-base-tags \
  --methylation-mode <em-seq|taps> --ref ref.fa
fgumi filter --ref ref.fa --min-reads 3 --max-base-error-rate 0.1 --min-methylation-depth 3 \
  --methylation-mode <em-seq|taps> --min-conversion-fraction 0.9
```

### Duplex (High Specificity)

```bash
fgumi duplex --min-reads 1 --min-input-base-quality 20 --output-per-base-tags \
  --methylation-mode <em-seq|taps> --ref ref.fa
fgumi filter --ref ref.fa --min-reads 10,5,3 --max-base-error-rate 0.1 --min-methylation-depth 10,5,3 \
  --require-single-strand-agreement --require-strand-methylation-agreement \
  --methylation-mode <em-seq|taps> --min-conversion-fraction 0.9
```

### Deduplication (No Consensus)

For workflows that mark duplicates without consensus calling:

```bash
fgumi dedup \
  --input sorted.bam \
  --output deduped.bam \
  --metrics metrics.txt
```

---

## Troubleshooting

### Low Family Sizes / Too Many UMI Groups

If family size histograms show many singletons:

1. Check that `--edits` is appropriate for your UMI length
2. For fixed UMIs, review correction metrics to see how many UMIs are being corrected vs rejected
3. **EM-Seq only:** verify that UMI sequences are synthesized with methylated cytosines to protect them from enzymatic conversion

### Missing MM/ML Tags on Output

Ensure both `--methylation-mode` and `--ref` are provided to the consensus caller. The reference FASTA must have an accompanying `.dict` file (generate with `samtools dict` if missing).

### Unexpected Masking from Strand Methylation Agreement

`--require-strand-methylation-agreement` only applies to duplex reads at CpG sites. If you see excessive masking:

1. Check that your library has adequate duplex coverage at CpG sites
2. Consider whether strand-specific methylation differences are biologically expected (e.g., imprinted regions)
3. This filter requires both strands to have evidence — positions with zero evidence on either strand are not masked

### Reads Filtered by Conversion Fraction

If many reads fail `--min-conversion-fraction`:

1. **EM-Seq:** this indicates potential issues with enzymatic conversion efficiency
2. **TAPs:** this indicates non-CpG cytosines are being converted, suggesting insufficient TAPs specificity
3. Try lowering the threshold (e.g., 0.8 instead of 0.9)
4. Check the overall conversion rate in your library QC metrics
5. Reads with no non-CpG cytosine positions (e.g., very short reads aligned to AT-rich regions) automatically pass this filter

### Using the Wrong Methylation Mode

If you use `--methylation-mode em-seq` for TAPs data (or vice versa), the methylation probabilities will be inverted — methylated positions will show low probability and vice versa. If downstream analysis shows unexpected methylation patterns, verify you used the correct mode for your chemistry.

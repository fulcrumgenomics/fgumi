# fgumi EM-Seq Pipeline Guide

This guide describes how to process enzymatic methyl-seq (EM-Seq) data through fgumi's consensus pipeline. It covers both workflows with and without fixed (known) UMI sets, for simplex and duplex consensus calling.

## Background

EM-Seq uses enzymatic conversion (TET2 + APOBEC) to distinguish methylated from unmethylated cytosines:

- **Methylated C** (5mC): protected by TET2, remains as C in sequencing reads
- **Unmethylated C**: deaminated by APOBEC to uracil, read as T after PCR

Conversion happens **before PCR**, so all reads in a UMI family agree on conversion status at each position. fgumi's EM-Seq mode tracks this per-base evidence through consensus calling and provides methylation-aware filtering.

### Impact on UMI Processing

C→T conversion affects consensus calling: at a reference C position, reads showing T are not errors — they represent unmethylated cytosines. Standard consensus calling would treat C/T disagreements as sequencing errors and penalize quality. EM-Seq mode recognizes these as conversion events.

**Note:** UMI sequences should be synthesized with methylated cytosines (5mC) to protect them from enzymatic conversion. If UMIs contain unmethylated cytosines, that is a library preparation issue.

---

## Pipeline Overview

The EM-Seq pipeline follows the same structure as the [standard consensus pipeline](best-practices.md), with additional flags at key steps:

```text
Phase 1: FASTQ → Grouped BAM
  extract → [correct] → fastq | bwameth | zipper → sort → group

Phase 2: Grouped BAM → Filtered Consensus
  simplex/duplex → fastq | bwameth | zipper → filter → sort
```

Steps marked with EM-Seq-specific flags:

| Step | Flag | Purpose |
|------|------|---------|
| `simplex` | `--em-seq --ref` | Methylation-aware consensus |
| `duplex` | `--em-seq --ref` | Methylation-aware consensus |
| `filter` | `--min-methylation-depth`, `--require-strand-methylation-agreement`, `--min-conversion-fraction` | Methylation-specific filtering |

---

## Workflow A: Random UMIs (No Fixed UMI Set)

This is the simpler case. Random UMIs (e.g., random 8-mers ligated during library prep) do not need correction against a whitelist.

### Step 1: UMI Extraction

Extract UMIs from FASTQ. No EM-Seq-specific flags needed here.

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
fgumi fastq --input unmapped.bam --no-read-suffix \
  | bwameth.py --reference ref.fa --threads 16 --interleaved /dev/stdin \
  | samtools view -b \
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

Use `--em-seq` and `--ref` to enable methylation-aware consensus.

**Simplex:**

```bash
fgumi simplex \
  --input grouped.bam \
  --output consensus.bam \
  --min-reads 1 \
  --min-input-base-quality 20 \
  --output-per-base-tags \
  --em-seq \
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
  --em-seq \
  --ref ref.fa \
  --threads 8
```

### Step 6: Re-alignment

Consensus reads are unmapped and must be re-aligned:

```bash
fgumi fastq --input consensus.bam --no-read-suffix \
  | bwameth.py --reference ref.fa --threads 16 --interleaved /dev/stdin \
  | samtools view -b \
  | fgumi zipper --unmapped consensus.bam --reference ref.fa --restore-unconverted-bases --output consensus.mapped.bam
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

Same as Workflow A:

```bash
fgumi extract \
  --inputs r1.fq.gz r2.fq.gz \
  --read-structures 8M+T +T \
  --sample "sample_name" \
  --library "library_name" \
  --output unmapped.bam \
  --threads 4
```

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

### Steps 3–8: Alignment through Final Sort

After correction, the remaining steps are the same as Workflow A (steps 2–8).

---

## EM-Seq Output Tags

When `--em-seq` is enabled, consensus reads carry additional BAM tags for methylation evidence.

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

The `MM`/`ML` tags follow the [SAM-spec methylation format](https://samtools.github.io/hts-specs/SAMtags.pdf) and are compatible with downstream methylation analysis tools.

---

## EM-Seq Filter Options

The `filter` command provides three methylation-specific options. These operate on the `cu`/`ct`/`au`/`at`/`bu`/`bt` count tags emitted by `--em-seq` consensus calling.

### `--min-methylation-depth`

Per-base masking based on methylation evidence depth. Bases where `cu[i] + ct[i]` is below the threshold are masked to N.

Accepts 1–3 comma-delimited values for duplex reads, following the same convention as `--min-reads`:

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

Read-level filter. Requires `--ref`. Accepts a value between 0.0 and 1.0.

Computes the C→T conversion fraction at non-CpG reference cytosine positions across the read. Non-CpG cytosines are expected to be unmethylated and therefore converted. A low conversion rate suggests incomplete enzymatic conversion, and the entire read is discarded.

CpG positions are excluded from this calculation because they may be methylated (and therefore correctly unconverted).

---

## Recommended Parameters

### EM-Seq Simplex (Moderate Stringency)

```bash
fgumi simplex --min-reads 1 --min-input-base-quality 20 --output-per-base-tags --em-seq --ref ref.fa
fgumi filter --ref ref.fa --min-reads 3 --max-base-error-rate 0.1 --min-methylation-depth 3 --min-conversion-fraction 0.9
```

### EM-Seq Duplex (High Specificity)

```bash
fgumi duplex --min-reads 1 --min-input-base-quality 20 --output-per-base-tags --em-seq --ref ref.fa
fgumi filter --ref ref.fa --min-reads 10,5,3 --max-base-error-rate 0.1 --min-methylation-depth 10,5,3 \
  --require-single-strand-agreement --require-strand-methylation-agreement --min-conversion-fraction 0.9
```

### EM-Seq Deduplication (No Consensus)

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
3. Verify that UMI sequences are synthesized with methylated cytosines to protect them from enzymatic conversion

### Missing MM/ML Tags on Output

Ensure both `--em-seq` and `--ref` are provided to the consensus caller. The reference FASTA must have an accompanying `.dict` file (generate with `samtools dict` if missing).

### Unexpected Masking from Strand Methylation Agreement

`--require-strand-methylation-agreement` only applies to duplex reads at CpG sites. If you see excessive masking:

1. Check that your library has adequate duplex coverage at CpG sites
2. Consider whether strand-specific methylation differences are biologically expected (e.g., imprinted regions)
3. This filter requires both strands to have evidence — positions with zero evidence on either strand are not masked

### Reads Filtered by Conversion Fraction

If many reads fail `--min-conversion-fraction`:

1. This indicates potential issues with enzymatic conversion efficiency
2. Try lowering the threshold (e.g., 0.8 instead of 0.9)
3. Check the overall conversion rate in your library QC metrics
4. Reads with no non-CpG cytosine positions (e.g., very short reads aligned to AT-rich regions) automatically pass this filter

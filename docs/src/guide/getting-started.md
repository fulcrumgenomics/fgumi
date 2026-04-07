# Getting Started

This guide walks through a basic fgumi workflow from FASTQ files to filtered consensus reads.

## Prerequisites

- fgumi installed (see [Installation](../index.md#installation))
- A reference genome FASTA (with BWA index)
- Paired-end FASTQ files with UMI sequences

## Basic Workflow

### 1. Extract UMIs from FASTQ

Extract UMIs from FASTQ reads and create an unmapped BAM. The `--read-structures` argument tells fgumi where UMI bases are located in each read. See [Read Structures](read-structures.md) for details.

```bash
fgumi extract \
  --inputs R1.fastq.gz R2.fastq.gz \
  --read-structures +T +M \
  --output unaligned.bam \
  --sample MySample \
  --library MyLibrary
```

### 2. (Optional) Correct UMIs

If using a fixed set of known UMIs, correct sequencing errors:

```bash
fgumi correct \
  --input unaligned.bam \
  --output corrected.bam \
  --umi-files umis.txt \
  --min-distance 1
```

### 3. Align and Sort

Use fgumi's streaming pipeline to align with BWA and sort into template-coordinate order in a single pass:

```bash
fgumi fastq --input unaligned.bam \
  | bwa mem -p ref.fa - \
  | fgumi zipper --unmapped unaligned.bam \
  | fgumi sort --output sorted.bam --order template-coordinate
```

This pipes reads through:
1. `fastq` — converts unmapped BAM to interleaved FASTQ
2. `bwa mem` — aligns reads to the reference
3. `zipper` — merges aligned reads with original unmapped BAM to restore UMI tags
4. `sort` — sorts into template-coordinate order for grouping

> **Note:** `fgumi zipper` accepts SAM input (piped from the aligner) or a BAM file via `--input`.
> Piping SAM directly from the aligner is preferred for best performance; BAM input is
> functional but involves an extra decode step.

For single-cell data, pass `--cell-tag CB` to `sort` to include the cell barcode in the
template-coordinate sort key, keeping templates from different cells at the same locus separate:

```bash
fgumi fastq --input unaligned.bam \
  | bwa mem -p ref.fa - \
  | fgumi zipper --unmapped unaligned.bam \
  | fgumi sort --output sorted.bam --order template-coordinate --cell-tag CB
```

### 3b. (Optional) Merge Multiple BAMs

If processing multiple lanes or flowcells separately, merge the sorted BAMs before grouping:

```bash
fgumi merge \
  --order template-coordinate \
  --output merged.bam \
  lane1_sorted.bam lane2_sorted.bam lane3_sorted.bam
```

All inputs must be sorted in the same order. For large numbers of files, use `--input-list`:

```bash
fgumi merge \
  --order template-coordinate \
  --input-list bam_paths.txt \
  --output merged.bam
```

For single-cell data, pass `--cell-tag CB` to include the cell barcode in the merge key.

### 4. Group Reads by UMI

Group reads from the same original molecule together.

For duplex workflows, use `paired` strategy:

```bash
fgumi group \
  --input sorted.bam \
  --output grouped.bam \
  --strategy paired
```

For simplex/codec workflows, use `adjacency` strategy:

```bash
fgumi group \
  --input sorted.bam \
  --output grouped.bam \
  --strategy adjacency
```

To collect all grouping QC metrics under a single prefix:

```bash
fgumi group \
  --input sorted.bam \
  --output grouped.bam \
  --strategy adjacency \
  --metrics group_metrics
```

This writes `group_metrics.family_sizes.txt`, `group_metrics.grouping_metrics.txt`, and
`group_metrics.position_group_sizes.txt` in one step.

See [UMI Grouping](umi-grouping.md) for details on grouping strategies.

### 5. Call Consensus Reads

Choose the consensus calling method based on your library preparation:

**Simplex consensus** (single-strand):
```bash
fgumi simplex \
  --input grouped.bam \
  --output consensus.bam
```

**Duplex consensus** (double-strand):
```bash
fgumi duplex \
  --input grouped.bam \
  --output duplex.bam
```

**CODEC consensus**:
```bash
fgumi codec \
  --input grouped.bam \
  --output codec_consensus.bam
```

See [Consensus Calling](consensus-calling.md) and [Duplex Consensus Calling](duplex-consensus-calling.md) for details.

### 6. (Optional) Collect QC Metrics

Collect QC metrics before filtering to understand your library.

**For simplex libraries**, use `simplex-metrics` on the grouped BAM:

```bash
fgumi simplex-metrics \
  --input grouped.bam \
  --output simplex_metrics
```

**For duplex libraries**, use `duplex-metrics` on the grouped BAM:

```bash
fgumi duplex-metrics \
  --input grouped.bam \
  --output duplex_metrics
```

Both commands write a set of metrics files under the given output prefix. See
[Working with Metrics](working-with-metrics.md) for details on interpreting the output.

### 7. Filter Consensus Reads

Filter consensus reads based on quality metrics. The `--min-reads` format depends on the
consensus type:

**For simplex consensus** (single integer):
```bash
fgumi filter \
  --input consensus.bam \
  --output filtered.bam \
  --ref ref.fa \
  --min-reads 1
```

**For duplex consensus** (three comma-separated values: duplex,AB,BA):
```bash
fgumi filter \
  --input duplex.bam \
  --output filtered.bam \
  --ref ref.fa \
  --min-reads 1,1,1
```

### 8. (Optional) Clip Overlapping Reads

Clip overlapping bases in read pairs to avoid double-counting evidence:

```bash
fgumi clip \
  --input filtered.bam \
  --output clipped.bam \
  --ref ref.fa
```

## What's Next

- [Best Practices](best-practices.md) — recommended parameter settings and pipeline configuration
- [Performance Tuning](performance-tuning.md) — threading, memory, and compression optimization
- [Working with Metrics](working-with-metrics.md) — understanding fgumi's output metrics

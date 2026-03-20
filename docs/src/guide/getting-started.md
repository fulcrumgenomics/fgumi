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

### 4. Group Reads by UMI

Group reads from the same original molecule together:

```bash
fgumi group \
  --input sorted.bam \
  --output grouped.bam \
  --strategy paired     # for duplex workflows
  # or --strategy adjacency for simplex/codec workflows
```

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

### 6. (Optional) Collect Duplex Metrics

For duplex libraries, collect QC metrics before filtering:

```bash
fgumi duplex-metrics \
  --input grouped.bam \
  --output metrics
```

### 7. Filter Consensus Reads

Filter consensus reads based on quality metrics:

```bash
fgumi filter \
  --input consensus.bam \
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

# Getting Started

This guide walks through a basic fgumi workflow from FASTQ files to filtered consensus reads.

## Prerequisites

- fgumi installed (see [Installation](../index.md#installation))
- A reference genome FASTA (with BWA index)
- Paired-end FASTQ files with UMI sequences

## Basic Workflow

This guide runs each stage as its own command so you can see the whole pipeline. Once you are comfortable with the stages, you can collapse supported runs of them into a single [`fgumi runall`](running-pipelines.md) command — same output, no intermediate files. We point out where along the way.

### 1. Extract UMIs from FASTQ

Extract UMIs from the FASTQ reads and create an unmapped BAM. The `--read-structures` argument tells fgumi where the UMI bases sit in each read; here `8M+T` means the first 8 bases of R1 are the UMI followed by template, and `+T` means R2 is all template. Adjust these to match your library — see [Read Structures](read-structures.md).

```bash
fgumi extract \
  --inputs R1.fastq.gz R2.fastq.gz \
  --read-structures 8M+T +T \
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

Aligners read and write plain reads, so they drop the UMI tags carried on the unmapped BAM. The workflow therefore aligns the reads, then uses `fgumi zipper` to copy the UMI (and other) tags from the original unmapped BAM back onto the aligned reads — you cannot pipe the aligner straight into `sort` without losing the UMIs. fgumi runs all of this as one streaming pass:

```bash
fgumi fastq --input unaligned.bam \
  | bwa mem -p ref.fa - \
  | fgumi zipper --unmapped unaligned.bam --reference ref.fa \
  | fgumi sort --output sorted.bam --order template-coordinate
```

This pipes reads through four stages:

1. `fastq` — converts the unmapped BAM to interleaved FASTQ
2. `bwa mem` — aligns reads to the reference
3. `zipper` — merges aligned reads with the original unmapped BAM to restore UMI tags
4. `sort` — sorts into template-coordinate order for grouping

> **Note:** `fgumi zipper` accepts SAM or BAM input, on stdin or via `--input`. For best
> performance, pipe uncompressed BAM from the aligner (e.g. `bwa-mem3 mem --bam=0`) — this
> skips both the SAM text formatting on the aligner side and the SAM parsing on the zipper
> side. SAM is fine for aligners that can't emit BAM; compressed BAM on a pipe is not
> recommended (wasted CPU on both ends).

For single-cell data, `fgumi sort` automatically adds the `CB` cell barcode tag to the
template-coordinate sort key, so reads from different cells at the same genomic position stay
separate. No extra flags are needed.

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

For single-cell data, the `CB` cell barcode tag is automatically included in the merge key.

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

Choose the consensus method based on your library preparation: `simplex` for single-strand UMI
libraries, `duplex` for double-stranded duplex libraries, and `codec` for CODEC. If you are
unsure, see [Choosing a Consensus Method](consensus-calling.md#choosing-a-consensus-method).

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

> **Do it in one command.** Steps 3–5 — align, sort, group, and call consensus — can run as a
> single streaming `fgumi runall` invocation that starts from the unmapped BAM produced in step 1,
> with no intermediate BAMs and byte-for-byte identical output:
>
> ```bash
> fgumi runall \
>   --start-from align --stop-after consensus --consensus simplex \
>   --input unaligned.bam --ref ref.fa --output consensus.bam \
>   --aligner::preset bwa \
>   --group::strategy adjacency --simplex::min-reads 1 \
>   --threads 8
> ```
>
> (`runall` runs the aligner for you via `--aligner::preset` and inserts the sort step
> automatically.) See
> [Running Pipelines](running-pipelines.md) for the full stage model and more examples. The
> step-by-step commands remain useful for inspecting intermediates or collecting per-stage metrics.

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

- [Running Pipelines](running-pipelines.md) — fuse these steps with `fgumi runall`, and when to prefer it over individual commands
- [Best Practices](best-practices.md) — recommended parameter settings and pipeline configuration
- [Performance Tuning](performance-tuning.md) — threading, memory, and compression optimization
- [Working with Metrics](working-with-metrics.md) — understanding fgumi's output metrics
- [Glossary](glossary.md) and [Troubleshooting](troubleshooting.md) — terminology reference and fixes for common errors

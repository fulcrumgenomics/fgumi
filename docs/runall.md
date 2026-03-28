# fgumi runall

`fgumi runall` runs any or all of the consensus-calling workflow in one command.
It streams data between stages in memory, eliminating intermediate BAM files and
significantly reducing wall-clock time compared to running each step separately.
The output is equivalent to running the individual commands in sequence.

## When to use runall

Use `runall` when you want to go from raw data to consensus reads without managing
intermediate files. The `--start-from` option lets you enter the workflow at any
stage, depending on where your data is:

| `--start-from` | Input | Stages run |
|---|---|---|
| `extract` | FASTQ files | extract, align, zipper, sort, group, consensus, filter |
| `correct` | Unmapped BAM with UMIs | correct, fastq, align, zipper, sort, group, consensus, filter |
| `fastq` | Unmapped BAM | fastq, align, zipper, sort, group, consensus, filter |
| `align` | FASTQ file | align, zipper, sort, group, consensus, filter |
| `zipper` | Mapped + unmapped BAM | zipper, sort, group, consensus, filter |
| `sort` | Mapped BAM | sort, group, consensus, filter |
| `group` | Grouped BAM (MI tags) | consensus, filter |

## Quick start

```bash
# Full pipeline from FASTQ to filtered consensus reads
fgumi runall -i R1.fq.gz R2.fq.gz -o consensus.bam -r ref.fa \
  --start-from extract \
  --extract::sample SampleA --extract::library LibA \
  --extract::read-structures 3M2S146T 3M2S146T \
  --aligner::preset bwa-mem2 -t 8

# From a mapped BAM (e.g. output of fgumi zipper)
fgumi runall -i mapped.bam -o consensus.bam -r ref.fa --start-from sort -t 8

# From unmapped + mapped BAMs
fgumi runall -i mapped.bam --unmapped unmapped.bam -o consensus.bam -r ref.fa \
  --start-from zipper -t 8

# From a grouped BAM (MI tags already assigned)
fgumi runall -i grouped.bam -o consensus.bam -r ref.fa --start-from group
```

### Early stopping

Use `--stop-after` to exit the workflow early. For example, to sort and group
without calling consensus:

```bash
fgumi runall -i mapped.bam -o grouped.bam -r ref.fa \
  --start-from sort --stop-after group
```

### Consensus modes

```bash
# Simplex (default)
fgumi runall -i grouped.bam -o consensus.bam -r ref.fa --consensus simplex

# Duplex
fgumi runall -i grouped.bam -o consensus.bam -r ref.fa --consensus duplex

# CODEC
fgumi runall -i grouped.bam -o consensus.bam -r ref.fa --consensus codec
```

## Performance

With 4 or more threads (`-t 4`), `runall` automatically uses a parallel
work-stealing scheduler for the group, consensus, filter, compress, and write
stages. With fewer threads it runs a simpler serial path.

The main performance advantage over running individual commands is eliminating
intermediate BAM I/O. A typical sequential workflow writes and re-reads ~3
intermediate BAMs; `runall` streams everything in memory.

## Scoped options

Stage-specific options use `::` scoping to avoid name collisions. For example,
`--sort::memory-limit` controls the sort stage while `--filter::min-reads`
controls filtering. Run `fgumi runall --help` for the full list.

Common options:

| Option | Default | Description |
|---|---|---|
| `--sort::memory-limit` | 805306368 | Memory limit in bytes for external sort |
| `--group::strategy` | `adjacency` | UMI assignment strategy |
| `--group::max-edits` | 1 | Maximum UMI mismatches for grouping |
| `--filter::min-reads` | 1 | Minimum supporting reads per consensus |
| `--consensus::error-rate-pre-umi` | 45 | Pre-UMI error rate (Phred) |
| `--consensus::error-rate-post-umi` | 40 | Post-UMI error rate (Phred) |
| `--extract::read-structures` | -- | Read structures, one per input FASTQ |
| `--aligner::preset` | -- | Aligner preset: `bwa-mem`, `bwa-mem2` |
| `--correct::max-mismatches` | 2 | Maximum mismatches for UMI correction |

## Equivalence guarantee

`fgumi runall` produces equivalent output to running the individual commands in
sequence. Integration tests verify this property across all `--start-from` modes
and consensus algorithms.

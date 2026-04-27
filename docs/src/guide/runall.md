# The `runall` Command

`fgumi runall` runs the full UMI-processing pipeline end-to-end in a single
process, streaming data between stages in memory. It replaces the need to
invoke each step (`extract`, `fastq`, aligner, `zipper`, `sort`, `group`,
consensus, `filter`) as separate subprocesses with intermediate BAM files.

## When to use `runall`

Use `runall` when:
- You're running the standard UMI pipeline end-to-end.
- You want to skip writing intermediate BAMs to disk.
- You want a single progress display and single metrics file.

Use the step-by-step chain when:
- You're debugging a specific stage in isolation.
- You need to customize a step with flags that `runall` doesn't expose.
- You're integrating into a workflow manager (Nextflow/Snakemake) where each
  step is a separate task.

## Entry points: `--start-from`

`runall` can start from several points depending on what you already have:

- `extract` (default) — start from raw FASTQ files; runs the full pipeline.
- `correct` — start from an unmapped BAM with UMIs; runs
  correct -> fastq -> aligner -> zipper -> sort -> group -> consensus -> filter.
- `fastq` — start from an unmapped BAM (UMIs already extracted); runs
  fastq -> aligner -> zipper -> sort -> ...
- `sort` — start from a mapped BAM; runs sort -> group -> consensus -> filter.
- `group` — start from a template-coordinate sorted BAM with MI tags; runs
  consensus -> filter.

## Stopping early: `--stop-after`

Stop after a specific stage to produce intermediate output. Valid values:
`extract`, `correct`, `fastq`, `align`, `zipper`, `sort`, `group`,
`consensus`, `filter` (default). The value must be at or after the
`--start-from` stage.

For example, `--start-from extract --stop-after sort` produces a
template-coordinate sorted BAM with UMIs extracted and reads aligned — useful
as input to other tools.

## Consensus modes: `--consensus`

- `simplex` (default) — single-strand consensus.
- `duplex` — duplex consensus (requires paired UMIs; see
  [UMI Grouping](umi-grouping.md#paired)).
- `codec` — CODEC consensus (R1 and R2 from opposite strands of the same
  molecule).

## Aligner configuration

`--aligner::preset` selects a built-in aligner preset:

- `bwa-mem` — invokes `bwa mem -p -K <chunk-size> -t <threads> <ref> /dev/stdin`.
- `bwa-mem2` — invokes `bwa-mem2 mem -p -K <chunk-size> -t <threads> <ref> /dev/stdin`.

`--aligner::command` lets you provide a custom aligner command template:

```bash
--aligner::command "minimap2 -ax sr {ref} -t {threads} - -"
```

Use `{ref}` for the reference path and `{threads}` for the thread count.
`--aligner::preset` and `--aligner::command` are mutually exclusive.

Additional aligner knobs:

- `--aligner::threads` — threads for the aligner subprocess (default: 4).
- `--aligner::chunk-size` — processing chunk size in bases (BWA `-K` flag).

## Examples

### Full pipeline from FASTQ to filtered consensus BAM

```bash
fgumi runall \
    --start-from extract \
    --input r1.fastq.gz r2.fastq.gz \
    --extract::sample SAMPLE --extract::library LIB \
    --extract::read-structures 8M+T 8M+T \
    --reference hs38DH.fa \
    --output consensus.bam
```

### Consensus only (input is already grouped)

```bash
fgumi runall \
    --start-from group \
    --input grouped.bam \
    --output consensus.bam \
    --consensus simplex
```

### Duplex consensus end-to-end

```bash
fgumi runall \
    --start-from extract \
    --input r1.fq.gz r2.fq.gz \
    --extract::sample S --extract::library L \
    --extract::read-structures 8M+T 8M+T \
    --reference hs38DH.fa \
    --consensus duplex \
    --group::strategy paired \
    --output duplex.bam
```

### Stop after sort to produce an aligned, sorted BAM

```bash
fgumi runall \
    --start-from extract --stop-after sort \
    --input r1.fq.gz r2.fq.gz \
    --extract::sample S --extract::library L \
    --extract::read-structures 8M+T 8M+T \
    --reference hs38DH.fa \
    --output aligned.sorted.bam
```

### Inspect the resolved plan without running

```bash
fgumi runall --start-from extract --explain \
    --input r1.fq.gz r2.fq.gz \
    --extract::sample S --extract::library L \
    --extract::read-structures 8M+T 8M+T \
    --reference hs38DH.fa \
    --output out.bam
```

Prints the full pipeline plan (source, stages with parameters, sink, runtime
limits) and exits without executing.

### Capture per-stage metrics

```bash
fgumi runall ... --metrics metrics.tsv
```

Writes a tab-separated file:

```
stage	wall_time_secs	records_in	records_out
GroupAssign	0.850	12	12
Consensus	0.003	12	3
Filter	0.001	3	3
```

## Performance

See the [Performance Tuning](performance-tuning.md) chapter for guidance on
thread counts, compression level, and sort memory.

## Signal handling

`Ctrl-C` during `runall` triggers clean shutdown: all workers observe the
cancel, stages stop, and the partial `.tmp` output file is discarded. No
corrupt BAM is left behind.

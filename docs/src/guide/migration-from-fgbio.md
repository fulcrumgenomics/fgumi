# Migration from fgbio

fgumi is the Rust successor to [fgbio](https://github.com/fulcrumgenomics/fgbio) for UMI-based tools. This guide maps fgbio tools to their fgumi equivalents and highlights key differences.

## Command Mapping

| fgbio Tool | fgumi Command | Notes |
|-----------|---------------|-------|
| `ExtractUmisFromBam` | `extract` | Extracts directly from FASTQ (not BAM) |
| `CorrectUmis` | `correct` | |
| `ZipperBams` | `zipper` | Also replaces `picard MergeBamAlignment`; accepts SAM or BAM input |
| `SortBam` | `sort` | Adds template-coordinate sort order with optional cell barcode key |
| `GroupReadsByUmi` | `group` | Same strategies: identity, edit, adjacency, paired |
| `CallMolecularConsensusReads` | `simplex` | |
| `CallDuplexConsensusReads` | `duplex` | |
| `CallCodecConsensusReads` | `codec` | |
| `FilterConsensusReads` | `filter` | |
| `ClipBam` | `clip` | |
| `CollectDuplexSeqMetrics` | `duplex-metrics` | |
| *(no equivalent)* | `simplex-metrics` | New: simplex QC metrics (yield, family sizes, UMI counts) |
| *(samtools merge)* | `merge` | k-way merge of pre-sorted BAMs; supports all sort orders |
| `ReviewConsensusVariants` | `review` | |

## Key Differences

### Input Format

fgbio's `ExtractUmisFromBam` takes an unmapped BAM as input. fgumi's `extract` takes FASTQ files directly, which is more common in practice and avoids an unnecessary BAM conversion step.

### Streaming Pipeline

fgumi supports Unix pipe-based streaming for the alignment workflow:

```bash
fgumi fastq --input unaligned.bam \
  | bwa mem -p -K 150000000 -Y ref.fa - \
  | fgumi zipper --unmapped unaligned.bam \
  | fgumi sort --output sorted.bam --order template-coordinate
```

This replaces multiple separate fgbio/picard steps (SortBam, ZipperBams/MergeBamAlignment) with a single streaming pass. `fgumi zipper` accepts SAM piped from the aligner (preferred) or a BAM file via `--input`.

### Merging Multiple BAMs

fgbio users who relied on `samtools merge` to combine per-lane BAMs before grouping should use
`fgumi merge` instead. It performs an equivalent k-way merge and correctly handles
template-coordinate order with cell barcodes:

```bash
# fgbio/samtools workflow
samtools merge -n merged.bam lane1.bam lane2.bam lane3.bam

# fgumi equivalent (also supports template-coordinate and queryname sort orders)
fgumi merge --order template-coordinate --output merged.bam \
  lane1.bam lane2.bam lane3.bam
```

### Simplex QC Metrics

fgbio has no equivalent to `fgumi simplex-metrics`. This command provides yield curves,
family size distributions, and UMI frequency statistics specifically for simplex sequencing
experiments, analogous to what `duplex-metrics` provides for duplex experiments.

### Threading Model

fgumi uses a multi-threaded pipeline architecture where reading, processing, and writing happen
concurrently. Most commands accept `--threads` to control parallelism. See
[Performance Tuning](performance-tuning.md) for details.

### Grouping Strategies

fgumi supports the same four UMI assignment strategies as fgbio:
- `identity` — exact UMI matching only
- `edit` — edit-distance clustering
- `adjacency` — directional adjacency (recommended for most use cases)
- `paired` — paired adjacency for duplex workflows

The algorithms are equivalent but fgumi's implementations are optimized for throughput.

### Group Metrics

fgumi's `group` command now produces a third metrics file beyond family sizes and grouping
metrics: `position_group_sizes.txt`, a histogram of how many UMI families appear at each
genomic position. This has no fgbio equivalent but is useful for detecting UMI exhaustion or
abnormal duplication patterns.

Use the `--metrics PREFIX` flag to write all three files in one step.

### Metrics Compatibility

fgumi's `simplex` and `duplex` stats output uses the same three-column key-value format as
fgbio's `CallMolecularConsensusReads`, allowing direct comparison with `fgumi compare metrics`.

### Sort Orders

fgumi's `sort` command supports the same sort orders as fgbio:
- `coordinate` — standard genomic coordinate sort
- `queryname` — sort by read name
- `template-coordinate` — sort by template 5' positions (required input for `group`)

For single-cell data, `fgumi sort --order template-coordinate --cell-tag CB` includes the cell
barcode in the sort key so that templates from different cells at the same locus are not
interleaved. fgbio's template-coordinate sort does not support this.

### Boolean Flag Values

fgumi boolean flags (e.g. `--output-per-base-tags`, `--trim`, `--require-single-strand-agreement`)
accept the following values: `true`/`false`, `yes`/`no`, `y`/`n`, `t`/`f` (case-insensitive).
fgbio uses standard true/false only.

### Removed Options

The `--sort-order` flag has been removed from `simplex` and `codec`. Output sort order for
consensus reads is determined by the downstream pipeline step (`zipper` + `sort`), not by the
consensus caller itself.

## What fgumi Does Not Replace

fgumi focuses on UMI-based tools. The following fgbio tools do **not** have fgumi equivalents:

- Non-UMI tools (e.g., `TrimFastq`, `ErrorRateByReadPosition`, `EstimatePoolingFractions`)
- VCF tools (e.g., `FilterSomaticVcf`, `HapTyper`)
- FASTQ/FASTA utilities (e.g., `FastqToBam`, `HardMaskFasta`)

Continue using fgbio for these tools.

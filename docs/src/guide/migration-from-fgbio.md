# Migration from fgbio

fgumi is the Rust successor to [fgbio](https://github.com/fulcrumgenomics/fgbio) for UMI-based tools. This guide maps fgbio tools to their fgumi equivalents and highlights key differences.

## Command Mapping

| fgbio Tool | fgumi Command | Notes |
|-----------|---------------|-------|
| `ExtractUmisFromBam` | `extract` | Extracts directly from FASTQ (not BAM) |
| `CorrectUmis` | `correct` | |
| `ZipperBams` | `zipper` | Also replaces `picard MergeBamAlignment` |
| `SortBam` | `sort` | Adds template-coordinate sort order |
| `GroupReadsByUmi` | `group` | Same strategies: identity, edit, adjacency, paired |
| `CallMolecularConsensusReads` | `simplex` | |
| `CallDuplexConsensusReads` | `duplex` | |
| `CallCodecConsensusReads` | `codec` | |
| `FilterConsensusReads` | `filter` | |
| `ClipBam` | `clip` | |
| `CollectDuplexSeqMetrics` | `duplex-metrics` | |
| `ReviewConsensusVariants` | `review` | |

## Key Differences

### Input Format

fgbio's `ExtractUmisFromBam` takes an unmapped BAM as input. fgumi's `extract` takes FASTQ files directly, which is more common in practice and avoids an unnecessary BAM conversion step.

### Streaming Pipeline

fgumi supports Unix pipe-based streaming for the alignment workflow:

```bash
fgumi fastq --input unaligned.bam \
  | bwa mem -p ref.fa - \
  | fgumi zipper --unmapped unaligned.bam \
  | fgumi sort --output sorted.bam --order template-coordinate
```

This replaces multiple separate fgbio/picard steps (SortBam, ZipperBams/MergeBamAlignment) with a single streaming pass.

### Threading Model

fgumi uses a multi-threaded pipeline architecture where reading, processing, and writing happen concurrently. Most commands accept `--threads` to control parallelism. See [Performance Tuning](performance-tuning.md) for details.

### Grouping Strategies

fgumi supports the same four UMI assignment strategies as fgbio:
- `identity` — exact UMI matching only
- `edit` — edit-distance clustering
- `adjacency` — directional adjacency (recommended for most use cases)
- `paired` — paired adjacency for duplex workflows

The algorithms are equivalent but fgumi's implementations are optimized for throughput.

### Metrics Compatibility

fgumi's `simplex` and `duplex` stats output uses the same three-column key-value format as fgbio's `CallMolecularConsensusReads`, allowing direct comparison with `fgumi compare metrics`.

### Sort Orders

fgumi's `sort` command supports the same sort orders as fgbio:
- `coordinate` — standard genomic coordinate sort
- `queryname` — sort by read name
- `template-coordinate` — sort by template 5' positions (required input for `group`)

## What fgumi Does Not Replace

fgumi focuses on UMI-based tools. The following fgbio tools do **not** have fgumi equivalents:

- Non-UMI tools (e.g., `TrimFastq`, `ErrorRateByReadPosition`, `EstimatePoolingFractions`)
- VCF tools (e.g., `FilterSomaticVcf`, `HapTyper`)
- FASTQ/FASTA utilities (e.g., `FastqToBam`, `HardMaskFasta`)

Continue using fgbio for these tools.

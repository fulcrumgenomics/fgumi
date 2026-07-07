<h1><span style="color:#26a8e0">fg</span><span style="color:#38b44a">umi</span></h1>

High-performance tools for UMI-tagged sequencing data: UMI extraction and correction, alignment support, grouping, deduplication, simplex/duplex/CODEC consensus calling, filtering, and QC metrics.

![fgumi Pipeline](images/fgumi_subway.png)

The diagram shows the workflow from FASTQ files to filtered consensus reads:
- **Red**: Simplex (single-strand) consensus
- **Blue**: Duplex (double-strand) consensus
- **Green**: CODEC consensus
- **Orange**: Optional UMI correction for fixed UMI sets

## Where to Use fgumi

### Command Line

Install and run fgumi directly on your data. See the [Getting Started](guide/getting-started.md) guide.

### Nextflow Pipeline

Use [fastquorum](https://github.com/fulcrumgenomics/fastquorum) for an end-to-end Nextflow workflow from FASTQ to consensus reads using fgumi.

### Latch.bio

Run fgumi in the cloud with a point-and-click interface via [Latch.bio](https://latch.bio) — no installation required.

## Installation

### Pre-built Binaries

Pre-built binaries for common operating systems and architectures are attached to each [release](https://github.com/fulcrumgenomics/fgumi/releases/latest).

### Cargo

```bash
cargo install fgumi
```

### Bioconda

```bash
conda install -c bioconda fgumi
```

### From Source

```bash
git clone https://github.com/fulcrumgenomics/fgumi
cd fgumi
cargo build --release
```

## Available Commands

> **Running more than one stage?** Use [`fgumi runall`](guide/running-pipelines.md) to fuse the whole pipeline into a single streaming command — no intermediate BAMs, identical output. Use the individual commands below for a single stage, to inspect an intermediate, or to collect per-stage metrics.

The commands below are grouped by pipeline stage.

| Stage | Command | Description |
|-------|---------|-------------|
| UMI extraction | `extract` | Extract UMIs from FASTQ and create an unmapped BAM |
| UMI extraction | `correct` | Correct UMIs in a BAM to a fixed set of UMIs |
| Alignment support | `fastq` | Convert BAM to FASTQ (interleaved or paired), optionally embedding the UMI in the read name |
| Alignment support | `zipper` | Merge an aligned BAM with its unmapped BAM, restoring tags |
| Alignment support | `sort` | Sort by coordinate, queryname, or template-coordinate |
| Alignment support | `merge` | Merge pre-sorted BAM files into one sorted BAM |
| Grouping & dedup | `group` | Group reads by UMI to identify reads from the same molecule |
| Grouping & dedup | `dedup` | Mark or remove PCR duplicates using UMI information |
| Consensus | `simplex` | Call single-strand consensus reads |
| Consensus | `duplex` | Call duplex (double-strand) consensus reads |
| Consensus | `codec` | Call CODEC duplex consensus reads |
| Post-consensus | `filter` | Filter and mask consensus reads by quality |
| Post-consensus | `clip` | Clip overlapping bases in read pairs |
| QC metrics | `simplex-metrics` | Collect QC metrics for simplex data |
| QC metrics | `duplex-metrics` | Collect QC metrics for duplex data |
| QC metrics | `review` | Extract data to review consensus variant calls |
| Utilities | `downsample` | Downsample a BAM by UMI family |
| Utilities | `runall` | Fused multi-stage pipeline (extract → … → filter, no intermediate BAMs) |

See the [Tool Reference](tools/README.md) for detailed documentation of each command, and [Running Pipelines](guide/running-pipelines.md) for `runall`.

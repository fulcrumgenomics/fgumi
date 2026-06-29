![Build](https://github.com/fulcrumgenomics/fgumi/actions/workflows/check.yml/badge.svg)
[![Latest release](https://img.shields.io/github/downloads-pre/fulcrumgenomics/fgumi/latest/total)](https://github.com/fulcrumgenomics/fgumi/releases)
[![Version at crates.io](https://img.shields.io/crates/v/fgumi)](https://crates.io/crates/fgumi)
[![Documentation at docs.rs](https://img.shields.io/docsrs/fgumi)](https://docs.rs/fgumi)
[![codecov](https://codecov.io/gh/fulcrumgenomics/fgumi/graph/badge.svg)](https://codecov.io/gh/fulcrumgenomics/fgumi)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/fgumi.svg?label=bioconda)](https://bioconda.github.io/recipes/fgumi/README.html)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18702470.svg)](https://doi.org/10.5281/zenodo.18702470)

# fgumi

**⚠️ RESEARCH PREVIEW - USE AT YOUR OWN RISK**

This software is currently a **research preview**. While we have extensively
tested these tools across a wide variety of vendor-provided data, we make **no
guarantees** regarding correctness or stability.

High-performance tools for UMI-tagged sequencing data, from Fulcrum Genomics.

<p>
<a href="https://fulcrumgenomics.com">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/fulcrumgenomics/fgumi/main/.github/logos/fulcrumgenomics-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="https://raw.githubusercontent.com/fulcrumgenomics/fgumi/main/.github/logos/fulcrumgenomics-light.svg">
  <img alt="Fulcrum Genomics" src="https://raw.githubusercontent.com/fulcrumgenomics/fgumi/main/.github/logos/fulcrumgenomics-light.svg" height="100">
</picture>
</a>
</p>

[Visit us at Fulcrum Genomics](https://www.fulcrumgenomics.com) to learn more about how we can power your Bioinformatics with fgumi and beyond.

<a href="mailto:contact@fulcrumgenomics.com?subject=[GitHub inquiry]"><img src="https://img.shields.io/badge/Email_us-%2338b44a.svg?&style=for-the-badge&logo=gmail&logoColor=white"/></a>
<a href="https://www.fulcrumgenomics.com"><img src="https://img.shields.io/badge/Visit_Us-%2326a8e0.svg?&style=for-the-badge&logo=wordpress&logoColor=white"/></a>

## Overview

fgumi provides comprehensive functionality for:
- **UMI extraction** from FASTQ files
- **Read grouping** by UMI with multiple assignment strategies
- **UMI-aware deduplication** for marking/removing PCR duplicates
- **Consensus calling** (simplex, duplex, and CODEC)
- **Quality filtering** of consensus reads
- **Read clipping** for overlapping pairs
- **Metrics collection** for QC and analysis

## Pipeline Overview

<p align="center">
  <img src="https://raw.githubusercontent.com/fulcrumgenomics/fgumi/main/docs/images/fgumi_subway.png" alt="fgumi Pipeline" width="800"/>
</p>

The diagram shows the workflow from FASTQ files to filtered consensus reads:
- **Red**: Simplex (single-strand) consensus
- **Blue**: Duplex (double-strand) consensus
- **Green**: CODEC consensus
- **Orange**: Optional UMI correction for fixed UMI sets

## Where to Use fgumi

- **Command line** — install and run fgumi directly on your data; see the [Getting Started guide](docs/src/guide/getting-started.md).
- **Nextflow** — [fastquorum](https://github.com/fulcrumgenomics/fastquorum) provides an end-to-end FASTQ-to-consensus workflow built on fgumi.
- **Latch.bio** — run fgumi in the cloud through a point-and-click interface via [Latch.bio](https://latch.bio), no installation required.

## Resources

* [Documentation](https://docs.rs/fgumi)
* [Getting Started](docs/src/guide/getting-started.md): Tutorial workflow from FASTQ to consensus
* [Running Pipelines](docs/src/guide/running-pipelines.md): `fgumi runall` and when to use it
* [Best Practice Pipeline](docs/src/guide/best-practices.md): Recommended workflow from FASTQ to consensus
* [Performance Tuning Guide](docs/src/guide/performance-tuning.md): Threading, memory, and compression optimization
* [Snakemake Pipeline](docs/FastqToConsensus-RnD.smk): Reference implementation
* [Metrics](docs/src/guide/working-with-metrics.md): Output metrics documentation
* [Developing](docs/DEVELOPING.md): Developer guide
* [Compare CLI](docs/compare-cli.md): Compare command documentation (feature-gated)
* [Simulate CLI](docs/simulate-cli.md): Simulate command documentation (feature-gated)
* [Releases](https://github.com/fulcrumgenomics/fgumi/releases)
* [Issues](https://github.com/fulcrumgenomics/fgumi/issues): Report a bug or request a feature
* [Pull requests](https://github.com/fulcrumgenomics/fgumi/pulls): Submit a patch or new feature
* [Discussions](https://github.com/fulcrumgenomics/fgumi/discussions): Ask a question
* [Contributors guide](CONTRIBUTING.md)
* [License](LICENSE): Released under the MIT license

## Installation

### Downloading a pre-built binary

Pre-built binaries for the most common operating systems and CPU architectures are attached to each [release](https://github.com/fulcrumgenomics/fgumi/releases/latest) for this project.

### Installing with Cargo

```
cargo install fgumi
```

### Building from source

Clone the repository:

```
git clone https://github.com/fulcrumgenomics/fgumi
```

Build the `release` version:

```
cd fgumi
cargo build --release
```

### Optional Features

| Feature | Description |
|---------|-------------|
| `compare` | Developer tools for comparing BAMs and metrics |
| `simulate` | Commands for generating synthetic test data |
| `profile-adjacency` | Enable profiling output for adjacency UMI assigner |
Enable with: `cargo build --release --features <feature>`

## Available Tools

| Command | Description | Equivalent Tool(s) |
|---------|-------------|------------------|
| `extract` | Extract UMIs from FASTQ and create an unmapped BAM | `fgbio FastqToBam` + `fgbio ExtractUmisFromBam` |
| `correct` | Correct UMIs in a BAM to a fixed set of UMIs | `fgbio CorrectUmis` |
| `fastq` | Convert BAM to interleaved FASTQ for an aligner | `samtools fastq` |
| `zipper` | Merge an aligned BAM with its unmapped BAM, restoring tags | `fgbio ZipperBams`, `picard MergeBamAlignment` |
| `sort` | Sort by coordinate, queryname, or template-coordinate | `samtools sort` |
| `merge` | Merge pre-sorted BAM files into one sorted BAM | `samtools merge` |
| `group` | Group reads by UMI | `fgbio GroupReadsByUmi` |
| `dedup` | Mark/remove UMI-aware duplicates | `gatk UmiAwareMarkDuplicatesWithMateCigar`, `umi-tools dedup` |
| `simplex` | Call single-strand consensus reads | `fgbio CallMolecularConsensusReads` |
| `duplex` | Call duplex consensus reads | `fgbio CallDuplexConsensusReads` |
| `codec` | Call CODEC consensus | `fgbio CallCodecConsensusReads` |
| `filter` | Filter consensus reads | `fgbio FilterConsensusReads` |
| `clip` | Clip overlapping read pairs | `fgbio ClipBam` |
| `simplex-metrics` | Collect simplex QC metrics | — |
| `duplex-metrics` | Collect duplex QC metrics | `fgbio CollectDuplexSeqMetrics` |
| `review` | Review consensus variants | `fgbio ReviewConsensusVariants` |
| `downsample` | Downsample BAM by UMI family | N/A |
| `runall` | Fused multi-stage pipeline (extract → … → filter, no intermediate BAMs) | N/A |
| `compare <cmd>` | Compare files (feature-gated) | N/A |
| `simulate <cmd>` | Generate test data (feature-gated) | N/A |

## Usage

For detailed usage of each command, run:
```bash
fgumi <command> --help
```

### Basic Workflow

> **Tip:** the multi-stage steps below can run as a single `fgumi runall` command — same output, no intermediate BAMs. See [Running Pipelines](docs/src/guide/running-pipelines.md). The individual commands are shown here so each stage is visible and tunable.

1. **Extract UMIs** from FASTQ:
```bash
fgumi extract \
  --inputs R1.fastq.gz R2.fastq.gz \
  --read-structures +T +M \
  --output unaligned.bam \
  --sample MySample \
  --library MyLibrary
```

2. **(Optional) Correct UMIs** for fixed UMI sets:
```bash
fgumi correct \
  --input unaligned.bam \
  --output corrected.bam \
  --umi-files umis.txt \
  --min-distance 1
```

3. **Align and sort** reads using fgumi fastq + zipper + sort pipeline:
```bash
fgumi fastq --input unaligned.bam \
  | bwa mem -p ref.fa - \
  | fgumi zipper --unmapped unaligned.bam --reference ref.fa \
  | fgumi sort --output sorted.bam --order template-coordinate
```

4. **Group** reads by UMI:
```bash
fgumi group \
  --input sorted.bam \
  --output grouped.bam \
  --strategy paired   # for duplex workflows
  # or --strategy adjacency for simplex/codec workflows
```

5. **Call consensus** reads:
```bash
# Simplex consensus
fgumi simplex \
  --input grouped.bam \
  --output consensus.bam

# Or duplex consensus
fgumi duplex \
  --input grouped.bam \
  --output duplex.bam

# Or codec consensus
fgumi codec \
  --input grouped.bam \
  --output codec_consensus.bam
```

6. **(Optional) Collect duplex metrics**:
```bash
fgumi duplex-metrics \
  --input grouped.bam \
  --output metrics
```

7. **Filter** consensus reads:
```bash
fgumi filter \
  --input consensus.bam \
  --output filtered.bam \
  --ref ref.fa \
  --min-reads 1,1,1
```

## Performance Options

fgumi supports multi-threading and memory management for optimal performance:

- **Threading**: `--threads N` for parallel processing
- **Memory**: `--max-memory 768` (plain numbers are MiB; supports human-readable formats like `2GB`; pass `auto` to size to the host)
- **Compression**: `--compression-level 1-12` for speed vs size trade-offs

> **Memory Model — Important for fgbio users:**
> Unlike fgbio's JVM `-Xmx` which sets a hard ceiling on total process memory,
> fgumi's `--max-memory` controls pipeline queue backpressure only. Two things
> to be aware of:
>
> 1. **Per-thread scaling (default):** `--max-memory 768 --threads 8` allocates
>    768 MiB **per thread** = ~6 GB total queue memory. Use `--memory-per-thread false`
>    for a fixed total budget, or `--max-memory auto` to size the budget to the host
>    (cgroup-aware, minus `--memory-reserve`) so the command self-throttles instead of OOM-ing.
> 2. **Queue memory < total process memory:** Actual RSS will be higher due to
>    UMI data structures, decompressors, thread stacks, and working buffers.
>
> See the [Performance Tuning Guide](docs/src/guide/performance-tuning.md)
> for detailed guidance, including scenario-based configurations and troubleshooting.

## Performance

fgumi is written in Rust and tuned for high-throughput BAM processing.

### Command-Level Optimizations

| Command | Key Optimizations |
|---------|-------------------|
| `extract` | Work-stealing thread pool, streaming I/O |
| `correct` | N-gram indexing with pigeonhole principle, BK-tree for k>1 |
| `group` | 2-bit UMI encoding, N-gram/BK-tree indexing, directed adjacency graph |
| `simplex` | Fast-path for unanimous consensus, parallel processing |
| `duplex` | Parallel duplex calling, efficient strand matching |
| `codec` | Parallel CODEC consensus |
| `filter` | Streaming filter with parallel processing |
| `clip` | Parallel overlap detection and clipping |
| `sort` | External merge sort, configurable memory limit |

### General Optimizations

- **2-bit DNA encoding**: 4 bases in 1 byte, 32 bases in u64
- **CPU intrinsics**: XOR + popcount for Hamming distance
- **Work-stealing scheduler**: Unified pipeline with dynamic load balancing
- **libdeflate**: Fast BGZF compression

## Acknowledgements

fgumi's UMI grouping algorithms are inspired by:

- [UMI-tools](https://github.com/CGATOxford/UMI-tools) (Smith et al. 2017) - The directed adjacency
  method for UMI deduplication with count gradient constraints.
- [UMICollapse](https://github.com/Daniel-Liu-c0deb0t/UMICollapse) (Liu 2019) - N-gram and BK-tree
  indexing strategies for efficient similarity search in UMI deduplication.

## Authors

- [Nils Homer](https://github.com/nh13)
- [Tim Fennell](https://github.com/tfenne)

## Sponsors

Development of fgumi is supported by [Fulcrum Genomics](https://www.fulcrumgenomics.com).

[Become a sponsor](https://github.com/sponsors/fulcrumgenomics)

## Disclaimer

This software is under active development.
While we make a best effort to test this software and to fix issues as they are reported, this software is provided as-is without any warranty (see the [license](LICENSE) for details).
Please submit an [issue](https://github.com/fulcrumgenomics/fgumi/issues), and better yet a [pull request](https://github.com/fulcrumgenomics/fgumi/pulls) as well, if you discover a bug or identify a missing feature.
Please contact [Fulcrum Genomics](https://www.fulcrumgenomics.com) if you are considering using this software or are interested in sponsoring its development.

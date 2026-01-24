# fgumi

> [!CAUTION]
> **ALPHA SOFTWARE - USE AT YOUR OWN RISK**
>
> This software is currently in **ALPHA**. While we have extensively tested these
> tools across a wide variety of vendor-provided data, **no guarantees are made**
> regarding correctness or stability.
>
> We are targeting **June 1, 2026** to recommend fgumi over
> [fgbio](https://github.com/fulcrumgenomics/fgbio) for production use.

<!-- Add your logo here -->

[![GitHub workflow build status](https://img.shields.io/github/actions/workflow/status/fulcrumgenomics/fgumi/check.yml)](https://github.com/fulcrumgenomics/fgumi/actions?query=workflow%3A%22Check+and+Test%22)
[![Latest release](https://img.shields.io/github/downloads-pre/fulcrumgenomics/fgumi/latest/total)](https://github.com/fulcrumgenomics/fgumi/releases)
[![Version at crates.io](https://img.shields.io/crates/v/fgumi)](https://crates.io/crates/fgumi)
[![Documentation at docs.rs](https://img.shields.io/docsrs/fgumi)](https://docs.rs/fgumi)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

Fulcrum Genomics Unique Molecular Indexing (UMI) Tools - a suite of high-performance tools for working with UMI-tagged sequencing data.

fgumi provides comprehensive functionality for:
- **UMI extraction** from FASTQ files
- **Read grouping** by UMI with multiple assignment strategies
- **Consensus calling** (simplex, duplex, and vanilla)
- **Quality filtering** of consensus reads
- **Read clipping** for overlapping pairs
- **Metrics collection** for QC and analysis

## Resources

* [Documentation](https://docs.rs/fgumi)
* [Releases](https://github.com/fulcrumgenomics/fgumi/releases)
* [Issues](https://github.com/fulcrumgenomics/fgumi/issues): report a bug or request a feature
* [Pull requests](https://github.com/fulcrumgenomics/fgumi/pulls): submit a patch or new feature
* [Discussions](https://github.com/fulcrumgenomics/fgumi/discussions): as a question
* [Contributors guide](docs/CONTRIBUTING.md)
* [License](LICENSE): Released under the MIT license

## Installation

### Downloading a pre-built binary

Pre-built binaries for the most common operating systems and CPU architectures are attached to each [release](https://github.com/jdidion/pest-test/releases/latest) for this project.

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
| `isal` | Use ISA-L/igzip for BGZF compression (benchmarking) |

Enable with: `cargo build --release --features <feature>`

**Note:** The `isal` feature requires the vendored submodule. Clone with submodules:

```
git clone --recursive https://github.com/fulcrumgenomics/fgumi
```

Or initialize after cloning: `git submodule update --init`

## Available Tools

| Command | Description |
|---------|-------------|
| `extract` | Extract UMIs from FASTQ files and write to unaligned BAM |
| `group` | Group reads by UMI using various assignment strategies |
| `simplex` | Call single-strand consensus reads from UMI groups |
| `duplex` | Call duplex consensus from paired single-strand consensuses |
| `codec` | Call CODEC consensus (duplex from single read pairs) |
| `filter` | Filter and mask consensus reads based on quality metrics |
| `clip` | Clip overlapping read pairs to avoid double-counting |
| `correct` | Correct UMIs based on sequence similarity |
| `review` | Review and analyze consensus variants |
| `duplex-metrics` | Collect detailed metrics on duplex consensus |
| `zipper` | Restore original FASTQ from unaligned BAM |

## Usage

For detailed usage of each command, run:
```bash
fgumi <command> --help
```

### Basic Workflow

1. **Extract UMIs** from FASTQ:
```bash
fgumi extract \
  --inputs R1.fastq.gz R2.fastq.gz \
  --read-structures +T +M \
  --output unaligned.bam
```

2. **Align** reads (using your favorite aligner):
```bash
bwa mem ref.fa unaligned.bam | samtools sort -o aligned.bam
```

3. **Group** reads by UMI:
```bash
fgumi group \
  --input aligned.bam \
  --output grouped.bam \
  --strategy adjacency
```

4. **Call consensus** reads:
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

5. **Filter** consensus reads:
```bash
fgumi filter \
  --input consensus.bam \
  --output filtered.bam \
  --reference ref.fa
```


## Threading Options

fgumi provides two mutually exclusive threading modes to accommodate different computing environments:

### `--threads N` (Strict Mode)

**Use this on HPC clusters or shared systems** where you need predictable resource usage.

```bash
fgumi group --input in.bam --output out.bam --strategy adjacency --threads 8
```

- **Hard cap**: Guarantees the tool spawns at most N total threads
- **Profile-aware allocation**: Thread distribution is optimized for each command's workload
- **Example**: `--threads 8` spawns exactly 8 threads, distributed based on command profile

This mode is ideal when:
- Running on SLURM/PBS with allocated CPU limits
- Sharing a machine with other users
- You need reproducible resource usage for benchmarking

### Thread Allocation Summary

#### Profile-Aware Allocation (`--threads` mode)

Each command has a workload profile that determines how threads are distributed. This is based on profiling data showing where each command spends its CPU time:

| Profile | Commands | Distribution | Why |
|---------|----------|--------------|-----|
| **Work-Stealing** | `extract` | Dynamic | Threads balance read/process/write work automatically |
| **Writer-Heavy** | `clip` | ~60% writer, ~25% worker, ~15% reader | BGZF compression dominates |
| **Worker-Heavy** | `simplex`, `duplex` | ~60% worker, ~25% writer, ~15% reader | Consensus calling dominates |
| **Reader-Heavy** | `correct` | ~60% reader, ~25% worker, ~15% writer | BGZF decompression dominates |
| **Balanced** | `group`, `codec`, `filter` | ~35% each | Evenly distributed workload |

**Writer-Heavy** (`clip`):

| --threads | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 |
|-----------|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| Reader    | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 2 | 2  | 2  | 2  | 3  | 3  | 3  | 3  |
| Workers   | 1 | 1 | 1 | 1 | 1 | 1 | 2 | 2 | 2 | 2  | 2  | 3  | 3  | 3  | 3  | 4  |
| Writer    | 1 | 1 | 1 | 2 | 3 | 4 | 4 | 5 | 5 | 6  | 6  | 7  | 7  | 8  | 9  | 9  |

*Note: The `extract` command uses a work-stealing thread pool instead of fixed allocation, where all N threads dynamically balance between reading, processing, and writing work. This provides ~15-20% better performance than fixed allocation.*

**Worker-Heavy** (`simplex`, `duplex`):

| --threads | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 |
|-----------|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| Reader    | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 2 | 2  | 2  | 2  | 3  | 3  | 3  | 3  |
| Workers   | 1 | 1 | 1 | 2 | 3 | 3 | 4 | 5 | 5 | 6  | 6  | 7  | 7  | 8  | 9  | 9  |
| Writer    | 1 | 1 | 1 | 1 | 1 | 2 | 2 | 2 | 2 | 2  | 2  | 3  | 3  | 3  | 3  | 4  |

**Reader-Heavy** (`correct`):

| --threads | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 |
|-----------|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| Reader    | 1 | 1 | 1 | 2 | 3 | 3 | 4 | 5 | 5 | 6  | 6  | 7  | 7  | 8  | 9  | 9  |
| Workers   | 1 | 1 | 1 | 1 | 1 | 2 | 2 | 2 | 2 | 2  | 2  | 3  | 3  | 3  | 3  | 4  |
| Writer    | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 2 | 2  | 2  | 2  | 3  | 3  | 3  | 3  |

**Balanced** (`group`, `codec`, `filter`):

| --threads | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 |
|-----------|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| Reader    | 1 | 1 | 1 | 1 | 1 | 2 | 2 | 2 | 3 | 3  | 3  | 4  | 4  | 4  | 5  | 5  |
| Workers   | 1 | 1 | 1 | 2 | 2 | 2 | 3 | 3 | 3 | 3  | 3  | 4  | 4  | 4  | 5  | 5  |
| Writer    | 1 | 1 | 1 | 1 | 2 | 2 | 2 | 3 | 3 | 4  | 4  | 4  | 5  | 5  | 5  | 6  |

*Note: For N ≤ 3, each pool requires at least 1 thread, so the minimum is 3 threads total.*

### Default Behavior

If no threading option is specified, fgumi runs in **single-threaded mode**. For parallel processing, specify `--threads N`.

### Recommendations

| Scenario | Recommendation |
|----------|----------------|
| SLURM job with 8 CPUs allocated | `--threads 8` |
| Dedicated 16-core workstation | `--threads 16` |
| Shared login node | `--threads 2` (or avoid heavy jobs) |
| Snakemake/Nextflow with resource limits | `--threads {threads}` |

## Performance

fgumi is written in Rust for maximum performance. Key optimizations include:

- **N-gram indexing**: For the adjacency and paired grouping strategies, fgumi uses an N-gram index
  (with pigeonhole principle filtering) to achieve sub-linear candidate search when there are many
  unique UMIs per genomic position. This optimization is inspired by
  [UMICollapse](https://github.com/Daniel-Liu-c0deb0t/UMICollapse) (Liu 2019).
- **Multi-threading**: The `group`, `simplex`, `duplex`, `codec`, `filter`, and `extract` commands support parallel
  processing via `--threads N`. See [Threading Options](#threading-options) above.
- **Efficient bit encoding**: UMI sequences are encoded as 2-bit per base for fast Hamming distance
  computation using CPU intrinsics.

## Acknowledgements

fgumi's UMI grouping algorithms are inspired by:

- [UMI-tools](https://github.com/CGATOxford/UMI-tools) (Smith et al. 2017) - The directed adjacency
  method for UMI deduplication with count gradient constraints.
- [UMICollapse](https://github.com/Daniel-Liu-c0deb0t/UMICollapse) (Liu 2019) - N-gram and BK-tree
  indexing strategies for efficient similarity search in UMI deduplication.

## Disclaimer

This software is under active development.
While we make a best effort to test this software and to fix issues as they are reported, this software is provided as-is without any warranty (see the [license](LICENSE) for details).
Please submit an [issue](https://github.com/fulcrumgenomics/fgumi/issues), and better yet a [pull request](https://github.com/fulcrumgenomics/fgumi/pulls) as well, if you discover a bug or identify a missing feature.
Please contact [Fulcrum Genomics](https://www.fulcrumgenomics.com) if you are considering using this software or are interested in sponsoring its development.

<h1><span style="color:#26a8e0">fg</span><span style="color:#38b44a">umi</span></h1>

High-performance tools for UMI-tagged sequencing data: extraction, grouping, and consensus calling.

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

### Cargo

```bash
cargo install fgumi
```

This builds `fgumi` without the in-process aligner, so it needs no C++ toolchain and no AVX2. To include the aligner, enable the `aligner-bwa-mem3` feature. It needs a C++17 compiler and the zlib headers, and on x86_64 the resulting binary needs a CPU with AVX2:

```bash
cargo install fgumi --features aligner-bwa-mem3
```

Building the reference index still needs the `bwa-mem3` command-line tool (`bwa-mem3 index ref.fa`).

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

| Command | Description |
|---------|-------------|
| `extract` | Extract UMIs from FASTQ files |
| `correct` | Correct UMIs based on sequence similarity |
| `fastq` | Convert BAM to FASTQ (interleaved or paired split), optionally embedding the UMI in the read name |
| `zipper` | Merge alignments back onto the unmapped BAM, restoring its tags |
| `sort` | Sort BAM by coordinate/queryname/template |
| `group` | Group reads by UMI |
| `dedup` | Mark/remove UMI-aware duplicates |
| `simplex` | Call single-strand consensus reads |
| `duplex` | Call duplex consensus reads |
| `codec` | Call CODEC consensus |
| `filter` | Filter consensus reads |
| `runall` | Fuse two or more pipeline stages into a single streaming pass (see [Running Pipelines](guide/running-pipelines.md)) |
| `clip` | Clip overlapping read pairs |
| `duplex-metrics` | Collect duplex metrics |
| `review` | Review consensus variants |
| `downsample` | Downsample BAM by UMI family |
| `retag` | Rewrite SAM tags (copy/move/delete/pair) |
| `simplex-metrics` | Collect simplex metrics |
| `merge` | Merge sorted BAM files |

See the [Tool Reference](tools/index.md) for detailed documentation of each command.

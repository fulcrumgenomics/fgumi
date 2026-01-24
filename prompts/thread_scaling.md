# Thread Scaling Issues - Debugging Prompt

## Overview

Benchmark analysis of fgumi commit `b698874654e15f589776f9ea379fa6160f2d112d` revealed several operations with poor thread scaling or performance regressions at higher thread counts (t4/t8). These need investigation and optimization.

## Environment

- **Benchmarks repo**: `/Users/nhomer/work/fg-nils-00005/work/git/fgumi-benchmarks`
- **fgumi repo**: `/Users/nhomer/work/fg-nils-00005/work/git/fgumi`
- **fgumi binary**: `target/release/fgumi`
- **fgumi version**: `0.1.0-b698874654e15f589776f9ea379fa6160f2d112d-dirty`

```bash
# Set up environment variables for commands below
export BENCHMARKS=/Users/nhomer/work/fg-nils-00005/work/git/fgumi-benchmarks
cd /Users/nhomer/work/fg-nils-00005/work/git/fgumi
cargo build --release
```

## Related Issues

See also `prompts/memory_investigation.md` for memory tracking discrepancies that may contribute to thread scaling issues.

---

## Issue 1: duplex-metrics - Regression at t8

### Symptoms
| Sample | fgbio(s) | t1x | t4x | t8x | t8/t1 |
|--------|----------|-----|-----|-----|-------|
| SRR6109255 | 22.4 | 1.5 | 1.5 | **0.6** | 0.40 |
| SRR6109273 | 14.8 | 1.8 | 1.6 | 1.7 | 0.94 |
| synthetic-pipeline-xlarge | 378.5 | 2.0 | 1.8 | 1.9 | 0.95 |

**Critical**: On SRR6109255, t8 is **slower than fgbio** (0.6x) while t1 is faster (1.5x). This is a clear regression.

### Commands to reproduce

```bash
# Single-threaded (baseline)
/usr/bin/time -l target/release/fgumi duplex-metrics \
  --input $BENCHMARKS/data/processed/SRR6109255.group.paired.fgbio.bam \
  --output /tmp/test-duplex-metrics-t1

# 8 threads (regression)
/usr/bin/time -l target/release/fgumi duplex-metrics \
  --input $BENCHMARKS/data/processed/SRR6109255.group.paired.fgbio.bam \
  --output /tmp/test-duplex-metrics-t8
```

### Investigation areas
1. Check if `duplex-metrics` even uses threading (it may be single-threaded by design)
2. If threaded, look for lock contention in hash maps or counters
3. Profile with `samply` or `instruments` to find hotspots
4. Check for false sharing in data structures

### Relevant source files
- `src/commands/duplex_metrics.rs`

---

## Issue 2: duplex - Regression from t4 to t8

### Symptoms
| Sample | fgbio(s) | t1x | t4x | t8x | t8/t4 |
|--------|----------|-----|-----|-----|-------|
| synthetic-pipeline-xlarge | 207.7 | 2.3 | 5.4 | **2.2** | 0.41 |

The speedup **drops from 5.4x at t4 to 2.2x at t8** - adding more threads makes it slower.

### Commands to reproduce

```bash
# 4 threads (good)
/usr/bin/time -l target/release/fgumi duplex \
  --input $BENCHMARKS/data/processed/synthetic-pipeline-xlarge.group.paired.fgbio.bam \
  --output /tmp/test-duplex-t4.bam \
  --threads 4 \
  --min-reads 1 \
  --min-input-base-quality 10

# 8 threads (regression)
/usr/bin/time -l target/release/fgumi duplex \
  --input $BENCHMARKS/data/processed/synthetic-pipeline-xlarge.group.paired.fgbio.bam \
  --output /tmp/test-duplex-t8.bam \
  --threads 8 \
  --min-reads 1 \
  --min-input-base-quality 10
```

### Investigation areas
1. Lock contention in the unified pipeline
2. False sharing between threads (check struct alignment/padding)
3. Memory bandwidth saturation
4. Check queue sizes and backpressure in the pipeline

### Relevant source files
- `src/commands/duplex.rs`
- `src/lib/unified_pipeline/`

---

## Issue 3: zipper - Regression from t4 to t8

### Symptoms
| Sample | fgbio(s) | t1x | t4x | t8x | t8/t4 |
|--------|----------|-----|-----|-----|-------|
| SRR6109255 | 32.0 | 1.7 | 5.1 | **1.7** | 0.33 |
| synthetic-pipeline-xlarge | 882.6 | 1.8 | 3.0 | 3.8 | 1.27 |

On SRR6109255, performance **drops from 5.1x at t4 to 1.7x at t8** (same as t1).

### Commands to reproduce

```bash
# 4 threads (good)
samtools view -h $BENCHMARKS/data/processed/SRR6109255.aligned.bam | \
/usr/bin/time -l target/release/fgumi zipper \
  --unmapped $BENCHMARKS/data/processed/SRR6109255.extract.fgbio.bam \
  --reference $BENCHMARKS/data/resources/hs38DH.fa \
  --output /tmp/test-zipper-t4.bam \
  --threads 4

# 8 threads (regression)
samtools view -h $BENCHMARKS/data/processed/SRR6109255.aligned.bam | \
/usr/bin/time -l target/release/fgumi zipper \
  --unmapped $BENCHMARKS/data/processed/SRR6109255.extract.fgbio.bam \
  --reference $BENCHMARKS/data/resources/hs38DH.fa \
  --output /tmp/test-zipper-t8.bam \
  --threads 8
```

### Investigation areas
1. The zipper reads from stdin (piped from samtools) - check if input parsing is the bottleneck
2. Reference FASTA lookup may have contention
3. Check thread pool behavior with stdin input

### Relevant source files
- `src/commands/zipper.rs`

---

## Issue 4: extract - Plateaus at t4

### Symptoms
| Sample | fgbio(s) | t1x | t4x | t8x | t8/t4 |
|--------|----------|-----|-----|-----|-------|
| synthetic-pipeline-xlarge | 675.0 | 3.0 | 8.7 | **8.8** | 1.01 |

No improvement from t4 to t8 (8.7x → 8.8x).

### Commands to reproduce

```bash
# 4 threads
/usr/bin/time -l target/release/fgumi extract \
  --inputs $BENCHMARKS/data/raw/synthetic-pipeline-xlarge_1.fastq.gz \
           $BENCHMARKS/data/raw/synthetic-pipeline-xlarge_2.fastq.gz \
  --output /tmp/test-extract-t4.bam \
  --read-structures 8M+T 8M+T \
  --sample test --library lib \
  --threads 4

# 8 threads
/usr/bin/time -l target/release/fgumi extract \
  --inputs $BENCHMARKS/data/raw/synthetic-pipeline-xlarge_1.fastq.gz \
           $BENCHMARKS/data/raw/synthetic-pipeline-xlarge_2.fastq.gz \
  --output /tmp/test-extract-t8.bam \
  --read-structures 8M+T 8M+T \
  --sample test --library lib \
  --threads 8
```

### Investigation areas
1. FASTQ decompression may be I/O bound or single-threaded
2. Check if gzip decompression is parallelized (consider libdeflate or pigz integration)
3. BAM output compression may be the bottleneck
4. Profile to identify whether reading, processing, or writing is limiting

### Relevant source files
- `src/commands/extract.rs`

---

## General Debugging Tools

### Profiling with samply (recommended for macOS)
```bash
cargo install samply
samply record -- /path/to/fgumi <command> <args>
```

### Check pipeline stats
Add `--pipeline-stats` flag to fgumi commands to see internal timing breakdown.

### Memory profiling
```bash
# Use /usr/bin/time -l for peak memory
/usr/bin/time -l /path/to/fgumi <command>
# Look for "maximum resident set size" and "peak memory footprint"
```

---

## Summary Table

| Priority | Command | Issue | Impact |
|----------|---------|-------|--------|
| **Critical** | duplex-metrics | t8 slower than t1 | 0.4x regression |
| **High** | duplex | t8 slower than t4 | 0.41x t4→t8 |
| **High** | zipper | t8 slower than t4 | 0.33x t4→t8 |
| **Medium** | extract | No t4→t8 scaling | 1.01x plateau |

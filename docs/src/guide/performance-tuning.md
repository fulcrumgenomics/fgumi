# Performance Tuning

This chapter covers runtime characteristics of the fgumi pipeline with real
measured numbers, followed by guidance on threading, memory, and compression.

## Overview of runtime costs

When running the full pipeline, alignment (`bwa mem`) dominates. The following
numbers are from a 29 M paired-read Agilent HS2 run (`agilent-hs2`, ~16 GB of
gzipped FASTQ input) aligned to `hs38DH.fa` with `-t 8` throughout. See the
[Reproducing these numbers](#reproducing-these-numbers) section below for an
end-to-end script.

| stage                     | wall-clock (-t 8) | fraction |
|---------------------------|------------------:|---------:|
| bwa mem                   |             1152s |      81% |
| extract                   |              109s |       8% |
| zipper                    |               78s |       5% |
| sort                      |               27s |       2% |
| group                     |               14s |       1% |
| simplex                   |               31s |       2% |
| filter                    |                6s |      <1% |
| **total (step-by-step)**  |     **1416s**     |          |
| **total (runall)**        |     **1036s**     |          |
| **speedup (runall)**      |     **~27% (6:20 saved)** |  |

Output BAM sizes matched closely: 916 MB (step-by-step) vs. 908 MB (`runall`),
a ~1% difference attributable to non-determinism in `bwa mem -t 8`.

Key observation: `bwa` is the Amdahl's-law bottleneck for extract-based runs.
`runall` saves time by overlapping extract, zipper, sort, group, consensus,
and filter with `bwa`, not by making `bwa` faster.

### IDT cfDNA (smaller dataset)

On a smaller run (IDT cfDNA, ~1 GB FASTQ input), the same pattern held but with
smaller absolute savings:

- step-by-step at `-t 8`: ~66s
- `runall` at `-t 8`: ~48s
- speedup: ~12-15%

Relative savings scale with how much non-alignment work there is to overlap.

## When `runall` wins biggest

`runall` wins most when:
- Post-alignment stages (sort, group, consensus, filter) are non-trivial.
- You don't already have intermediate BAMs on disk — most of the savings come
  from eliminating intermediate BAM writes and reads.
- You have CPU headroom for the non-`bwa` stages to run concurrently with `bwa`.

`runall` wins less when:
- Alignment is already >90% of wall-clock — speedup is bounded by
  `1 / (bwa fraction)` (Amdahl).
- You already have intermediate BAMs and just need consensus + filter. In that
  case use `--start-from group` to skip the alignment portion entirely.

## Memory

Peak queue-held memory on the agilent-hs2 run was ~65 MiB. **This figure is
memory inside inter-stage queues only.** It does NOT include:
- the `bwa` subprocess (hundreds of MiB on a human reference),
- per-stage scratch buffers (BGZF decompression arenas, SAM parser state,
  ring buffers),
- items currently being processed inside a stage's `process()` call,
- payloads held in sink/coalesce reorder buffers after being popped,
- the allocator's overhead and fragmentation, thread stacks, I/O buffers.

See `src/lib/pipeline/memory.rs` for the exact scope of the tracker. Treat
queue-held memory as a **lower bound** on process RSS, not an estimate of it.

For a process-level memory figure, measure with `/usr/bin/time -v` (or
`/usr/bin/time -l` on macOS):

```bash
/usr/bin/time -v target/release/fgumi runall ...
```

and read `Maximum resident set size` from the output.

## Threads

`--threads` controls the fgumi worker pool used for BAM I/O and pipeline
stages. For `bwa`-dominated runs, match or slightly exceed your physical core
count — past that, contention with `bwa` hurts more than parallelism helps.

`--aligner::threads` controls `bwa`'s own thread pool (passed through as
`bwa mem -t N`). Usually set this to the machine's core count (independent of
the fgumi worker count; `bwa` is a subprocess and has its own scheduler).

## Coming from fgbio?

If you're used to fgbio's JVM-based memory model (`java -Xmx4g`), there are
important differences in how fgumi manages memory:

| | fgbio (JVM) | fgumi |
|---|---|---|
| **Memory control** | `-Xmx` sets a hard ceiling on the entire process | `--queue-memory` controls pipeline queue backpressure |
| **Enforcement** | Hard limit — JVM throws `OutOfMemoryError` at the ceiling | Soft limit — triggers backpressure to slow producers |
| **Scope** | Total process memory (heap + off-heap) | Queue memory only; does not cover UMI data structures, decompressors, thread stacks, or working buffers |
| **Scaling** | Fixed regardless of threads | Per-thread by default (`--queue-memory 768 --threads 8` = ~6 GB) |
| **Recommendation** | Set once and forget | Monitor RSS and adjust; use `--queue-memory-per-thread false` for a fixed total budget |

**Key takeaway:** fgumi's actual process memory (RSS) will be higher than the
`--queue-memory` value. When estimating memory needs, account for queue memory,
UMI grouping data structures, per-thread decompressor/compressor instances,
thread stacks, and I/O buffers.

For memory-constrained environments, start with `--queue-memory-per-thread false`
and a conservative total budget, then increase if throughput is too low.

## Threading options (detailed)

### No-flag Fast Path (default)
- **Usage**: Omit `--threads` entirely
- **Behavior**: Uses optimized single-threaded fast path with minimal overhead
- **Best for**: Small files, memory-constrained systems, debugging

### Explicit Single-threaded Mode
- **Usage**: `--threads 1`
- **Behavior**: Uses the unified pipeline with a single worker thread — same pipeline as `--threads N` but with N=1; does **not** use the no-flag fast path
- **Best for**: Isolating pipeline behavior in a single-threaded context

### Multi-threaded Mode
- **Usage**: `--threads N` where N > 1
- **Behavior**: Uses unified 7-step pipeline with work-stealing scheduler
- **Best for**: Large files, high-performance systems, production workloads

## Memory management

fgumi's unified memory management controls pipeline queue memory to prevent
out-of-memory conditions while maintaining throughput.

### Queue Memory Options

```bash
# Basic usage (768MB per thread - default)
fgumi filter --queue-memory 768 --queue-memory-per-thread true

# Human-readable formats
fgumi filter --queue-memory 2GB
fgumi filter --queue-memory 1024MiB

# Fixed total memory (no per-thread scaling)
fgumi filter --queue-memory 4096 --queue-memory-per-thread false
```

### Memory Scaling Behavior

| Threads | Per-thread Mode | Fixed Mode |
|---------|----------------|------------|
| 1       | 768MB          | 768MB      |
| 4       | 3GB            | 768MB      |
| 8       | 6GB            | 768MB      |
| 16      | 12GB           | 768MB      |

### Memory Validation

- **System check**: Warns if requesting >90% of available system memory
- **Overflow protection**: Prevents integer overflow with checked arithmetic
- **Decimal support**: Accepts formats like `1.5GB` in addition to integers

## Compression options

### Compression Level
- **Range**: 1 (fastest) to 12 (best compression)
- **Default**: 1 (fastest) for most commands; `fgumi merge` defaults to 6
- **Usage**: `--compression-level N`

### Compression Threading
- **Default**: Matches `--threads` setting
- **Override**: `--compression-threads N`
- **Best practice**: Usually leave at default

## Scenario-based configurations

### High-Throughput Server
**Goal**: Maximum processing speed for large datasets

```bash
fgumi filter \
  --threads 16 \
  --queue-memory 1GB \
  --compression-level 3 \
  --input large_dataset.bam \
  --output filtered.bam
```

### Memory-Constrained Node
**Goal**: Minimize memory usage while maintaining reasonable performance

```bash
fgumi filter \
  --threads 8 \
  --queue-memory 512 \
  --queue-memory-per-thread false \
  --compression-level 6 \
  --input dataset.bam \
  --output filtered.bam
```

### Fast Local SSD
**Goal**: Optimize for fast I/O with minimal compression overhead

```bash
fgumi filter \
  --threads 8 \
  --queue-memory 2GB \
  --compression-level 1 \
  --input dataset.bam \
  --output filtered.bam
```

### Network Storage
**Goal**: Minimize network I/O with maximum compression

```bash
fgumi filter \
  --threads 4 \
  --queue-memory 512 \
  --compression-level 9 \
  --input dataset.bam \
  --output filtered.bam
```

### Development/Testing
**Goal**: Fast iteration with minimal resource usage

```bash
fgumi filter \
  --queue-memory 256 \
  --compression-level 1 \
  --input small_test.bam \
  --output test_output.bam
```

## Verbose logging

Use `--verbose` (or `-v`) to enable debug-level logging for any command:

```bash
fgumi group --verbose --input reads.bam --output grouped.bam
```

This is equivalent to setting `RUST_LOG=debug`. If `RUST_LOG` is explicitly set,
it takes precedence over `--verbose`.

## Advanced pipeline options

The following options are available on all multi-threaded pipeline commands.
They are hidden from the default help text but can be useful for debugging and
performance analysis.

### Pipeline Statistics

```bash
fgumi group --pipeline-stats --input reads.bam --output grouped.bam
```

Prints detailed per-step timing, throughput, contention metrics, and per-thread
work distribution at completion.

### Scheduler Strategy

```bash
fgumi group --scheduler balanced-chase-drain --input reads.bam --output grouped.bam
```

Controls which scheduling strategy threads use for work assignment. The default
(`balanced-chase-drain`) is recommended for most workloads. Available
strategies:

| Strategy | Description |
|----------|-------------|
| `balanced-chase-drain` | Default. Balanced work distribution with output drain mode. |
| `fixed-priority` | Static thread roles (reader, writer, workers). Simple baseline. |
| `chase-bottleneck` | Threads dynamically follow work through the pipeline. |

Other experimental strategies are available (`thompson-sampling`, `ucb`,
`epsilon-greedy`, etc.) but are not recommended for production use.

### Deadlock Detection

```bash
# Adjust timeout (default: 10 seconds, 0 to disable)
fgumi group --deadlock-timeout 30 --input reads.bam --output grouped.bam

# Enable automatic recovery (default: detection only)
fgumi group --deadlock-recover --input reads.bam --output grouped.bam
```

The pipeline monitors for progress stalls. When no queue operations succeed for
the timeout duration, diagnostic information is logged (queue depths, memory
usage, per-queue timestamps).

With `--deadlock-recover`, the pipeline progressively doubles queue memory
limits (2x, 4x, up to 8x) to resolve backpressure deadlocks, then restores
original limits after 30 seconds of sustained progress.

## Command-specific considerations

### Extract
- Benefits from high memory (large FASTQ processing)
- Compression level affects output size significantly

### Zipper
- For best throughput, pipe SAM directly from the aligner rather than writing
  and re-reading a BAM
- BAM file input is supported via `--input` but adds an extra decode step
- The zipper pipeline uses raw-byte merging internally: aligned records are not
  fully decoded and re-encoded unless the record actually needs modification,
  which eliminates a significant CPU bottleneck on high-throughput runs

### Sort
- Uses an internal LoserTree (tournament tree) for k-way merging, which
  performs significantly better than a simple heap merge when the number of
  sorted runs is large
- `--max-memory` controls how much RAM is used for sort buffers; increase for
  large files to reduce the number of intermediate merge passes
- For template-coordinate sort with single-cell data, the `CB` tag is included
  automatically

### Merge
- `fgumi merge` performs a k-way merge using a LoserTree for efficient
  multi-file merging
- Thread count (`--threads`) controls compression parallelism, not merge
  concurrency
- For template-coordinate merges with single-cell data, the `CB` tag is
  included automatically

### Group/Dedup
- Memory usage scales with UMI diversity and the number of reads at any given
  position
- Higher thread counts improve UMI processing
- The `--metrics PREFIX` flag writes all grouping metrics in one step with
  minimal overhead

### Simplex/Duplex Metrics
- Both `simplex-metrics` and `duplex-metrics` are single-threaded; they do not
  benefit from `--threads`
- Memory usage is proportional to the number of unique genomic positions in
  the input

### Consensus (Simplex/Duplex/CODEC)
- Memory proportional to family sizes
- Benefits from balanced threading and memory

### Filter
- Streaming operation benefits from pipeline memory
- Compression affects final output size

## Reproducing these numbers

See `scripts/bench-repro.sh` to reproduce the agilent-hs2 measurements.
Requires:
- `bwa` and `samtools` on `PATH`
- Agilent HS2 FASTQ files (from vendor data or fgumi-benchmarks)
- `hs38DH.fa` reference with accompanying `bwa` index and samtools `.dict`

Usage:

```bash
./scripts/bench-repro.sh /tmp/fgumi-bench ./target/release/fgumi
```

The script runs the full step-by-step chain and the `runall` equivalent with
`/usr/bin/time` around each invocation, so you can compare per-stage and
end-to-end wall clock.

## Advanced: profile yourself

To get per-function flamegraphs of a real run:

```bash
cargo build --release
sudo perf record -F 99 -g -- target/release/fgumi runall ...
sudo perf report
```

Or via `flamegraph-rs`:

```bash
cargo install flamegraph
cargo flamegraph --bin fgumi -- runall --start-from extract ...
```

## Performance monitoring

### Memory Usage
- Monitor system memory usage during execution
- Watch for "exceeds available memory" warnings
- Adjust `--queue-memory` if seeing swap activity

### Thread Utilization
- Use `htop` or similar to monitor CPU usage
- All threads should show activity during processing
- Consider reducing threads if not fully utilized

### I/O Patterns
- Monitor disk I/O with `iotop`
- Network storage may benefit from lower thread counts
- SSD storage can handle higher thread counts

## Troubleshooting

### Out of Memory Errors
1. Reduce `--queue-memory`
2. Set `--queue-memory-per-thread false` for fixed limits
3. Reduce `--threads`

### Poor Performance
1. Increase `--threads` if CPU usage is low
2. Increase `--queue-memory` if I/O bound
3. Reduce `--compression-level` if CPU bound

### Pipeline Appears Stuck
If a command hangs without producing output:
1. Check if a deadlock warning appears in the log (default timeout: 10 seconds)
2. Run with `--verbose` to see detailed pipeline activity
3. Run with `--pipeline-stats` to see per-step metrics at completion
4. Try `--deadlock-recover` to allow automatic recovery from backpressure deadlocks
5. Reduce `--threads` — fewer threads means simpler scheduling and less contention

### System Memory Warnings
```text
Requested memory 16GB exceeds 90% of system memory (14.4GB)
```
- Reduce memory allocation or add more RAM
- Consider using `--queue-memory-per-thread false`

## Migration from legacy parameters

If using deprecated `--queue-memory-limit-mb`:

```bash
# Old (deprecated)
fgumi group --queue-memory-limit-mb 4096

# New (recommended)
fgumi group --queue-memory 4096 --queue-memory-per-thread false
```

The new parameters provide better control and human-readable formats while
maintaining backward compatibility.

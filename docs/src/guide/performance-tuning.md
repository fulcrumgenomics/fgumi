# Performance Tuning Guide

fgumi provides three key options to optimize performance for your system: threading, memory management, and compression. This guide explains how to configure these options for different scenarios.

## Coming from fgbio?

If you're used to fgbio's JVM-based memory model (`java -Xmx4g`), there are important differences in how fgumi manages memory:

| | fgbio (JVM) | fgumi |
|---|---|---|
| **Memory control** | `-Xmx` sets a hard ceiling on the entire process | `--queue-memory` controls pipeline queue backpressure |
| **Enforcement** | Hard limit — JVM throws `OutOfMemoryError` at the ceiling | Soft limit — triggers backpressure to slow producers |
| **Scope** | Total process memory (heap + off-heap) | Queue memory only; does not cover UMI data structures, decompressors, thread stacks, or working buffers |
| **Scaling** | Fixed regardless of threads | Per-thread by default (`--queue-memory 768 --threads 8` = ~6 GB) |
| **Recommendation** | Set once and forget | Monitor RSS and adjust; use `--queue-memory-per-thread false` for a fixed total budget |

**Key takeaway:** fgumi's actual process memory (RSS) will be higher than the `--queue-memory` value. When estimating memory needs, account for:
- Queue memory (controlled by `--queue-memory`)
- UMI grouping data structures (scales with UMI diversity and position depth)
- Per-thread decompressor and compressor instances
- Thread stacks and I/O buffers

For memory-constrained environments, start with `--queue-memory-per-thread false` and a conservative total budget, then increase if throughput is too low.

## Threading Options

### Single-threaded Mode
- **Usage**: `--threads 1` or omit the parameter
- **Behavior**: Uses optimized fast path with minimal overhead
- **Best for**: Small files, memory-constrained systems, debugging

### Multi-threaded Mode
- **Usage**: `--threads N` where N > 1
- **Behavior**: Uses unified 7-step pipeline with work-stealing scheduler
- **Best for**: Large files, high-performance systems, production workloads

## Memory Management

fgumi's unified memory management controls pipeline queue memory to prevent out-of-memory conditions while maintaining throughput.

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

## Compression Options

### Compression Level
- **Range**: 1 (fastest) to 12 (best compression)
- **Default**: 1 (fastest)
- **Usage**: `--compression-level N`

### Compression Threading
- **Default**: Matches `--threads` setting
- **Override**: `--compression-threads N`
- **Best practice**: Usually leave at default

## Scenario-Based Configurations

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

**Rationale**:
- High thread count for parallel processing
- Generous memory for pipeline buffers
- Lower compression for speed

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

**Rationale**:
- Moderate thread count
- Fixed memory limit (512MB total)
- Default compression for balance

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

**Rationale**:
- High memory for large pipeline buffers
- Minimal compression (I/O not bottleneck)

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

**Rationale**:
- Moderate threading to avoid overwhelming network
- Conservative memory usage
- Maximum compression to reduce network transfer

### Development/Testing
**Goal**: Fast iteration with minimal resource usage

```bash
fgumi filter \
  --queue-memory 256 \
  --compression-level 1 \
  --input small_test.bam \
  --output test_output.bam
```

**Rationale**:
- Single-threaded for simplicity
- Minimal memory footprint
- Fast compression for quick turnaround

## Verbose Logging

Use `--verbose` (or `-v`) to enable debug-level logging for any command:

```bash
fgumi group --verbose --input reads.bam --output grouped.bam
```

This is equivalent to setting `RUST_LOG=debug`. If `RUST_LOG` is explicitly set, it takes precedence over `--verbose`.

## Advanced Pipeline Options

The following options are available on all multi-threaded pipeline commands. They are hidden from the default help text but can be useful for debugging and performance analysis.

### Pipeline Statistics

```bash
fgumi group --pipeline-stats --input reads.bam --output grouped.bam
```

Prints detailed per-step timing, throughput, contention metrics, and per-thread work distribution at completion.

### Scheduler Strategy

```bash
fgumi group --scheduler balanced-chase-drain --input reads.bam --output grouped.bam
```

Controls which scheduling strategy threads use for work assignment. The default (`balanced-chase-drain`) is recommended for most workloads. Available strategies:

| Strategy | Description |
|----------|-------------|
| `balanced-chase-drain` | Default. Balanced work distribution with output drain mode. |
| `fixed-priority` | Static thread roles (reader, writer, workers). Simple baseline. |
| `chase-bottleneck` | Threads dynamically follow work through the pipeline. |

Other experimental strategies are available (`thompson-sampling`, `ucb`, `epsilon-greedy`, etc.) but are not recommended for production use.

### Deadlock Detection

```bash
# Adjust timeout (default: 10 seconds, 0 to disable)
fgumi group --deadlock-timeout 30 --input reads.bam --output grouped.bam

# Enable automatic recovery (default: detection only)
fgumi group --deadlock-recover --input reads.bam --output grouped.bam
```

The pipeline monitors for progress stalls. When no queue operations succeed for the timeout duration, diagnostic information is logged (queue depths, memory usage, per-queue timestamps).

With `--deadlock-recover`, the pipeline progressively doubles queue memory limits (2x, 4x, up to 8x) to resolve backpressure deadlocks, then restores original limits after 30 seconds of sustained progress.

## Performance Monitoring

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

## Command-Specific Considerations

### Extract
- Benefits from high memory (large FASTQ processing)
- Compression level affects output size significantly

### Group/Dedup
- Memory usage scales with UMI diversity
- Higher thread counts improve UMI processing

### Consensus (Simplex/Duplex/CODEC)
- Memory proportional to family sizes
- Benefits from balanced threading and memory

### Filter
- Streaming operation benefits from pipeline memory
- Compression affects final output size

## Migration from Legacy Parameters

If using deprecated `--queue-memory-limit-mb`:

```bash
# Old (deprecated)
fgumi group --queue-memory-limit-mb 4096

# New (recommended)
fgumi group --queue-memory 4096 --queue-memory-per-thread false
```

The new parameters provide better control and human-readable formats while maintaining backward compatibility.

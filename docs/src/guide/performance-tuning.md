# Performance Tuning Guide

fgumi provides three key options to optimize performance for your system: threading, memory management, and compression. This guide explains how to configure these options for different scenarios.

## Coming from fgbio?

If you're used to fgbio's JVM-based memory model (`java -Xmx4g`), there are important differences in how fgumi manages memory:

| | fgbio (JVM) | fgumi |
|---|---|---|
| **Memory control** | `-Xmx` sets a hard ceiling on the entire process | `--max-memory` controls pipeline queue backpressure |
| **Enforcement** | Hard limit — JVM throws `OutOfMemoryError` at the ceiling | Soft limit — triggers backpressure to slow producers |
| **Scope** | Total process memory (heap + off-heap) | Queue memory only; does not cover UMI data structures, decompressors, thread stacks, or working buffers |
| **Scaling** | Fixed regardless of threads | Per-thread by default (`--max-memory 768 --threads 8` = ~6 GB) |
| **Recommendation** | Set once and forget | Monitor RSS and adjust; use `--max-memory auto` (or `--memory-per-thread false`) for a host-aware/fixed total budget |

**Key takeaway:** fgumi's actual process memory (RSS) will be higher than the `--max-memory` value. When estimating memory needs, account for:
- Queue memory (controlled by `--max-memory`)
- UMI grouping data structures (scales with UMI diversity and position depth)
- Per-thread decompressor and compressor instances
- Thread stacks and I/O buffers

For memory-constrained environments, pass `--max-memory auto` to detect (cgroup-aware) host memory and subtract `--memory-reserve`, so the budget shrinks to fit the host. Alternatively, start with `--memory-per-thread false` and a conservative total budget, then increase if throughput is too low.

## Threading Options

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

## Memory Management

fgumi's unified memory management controls pipeline queue memory to prevent out-of-memory conditions while maintaining throughput.

### Memory Options

```bash
# Basic usage (768 MiB per thread - default; plain numbers are MiB)
fgumi filter --max-memory 768 --memory-per-thread true

# Human-readable formats
fgumi filter --max-memory 2GB
fgumi filter --max-memory 1024MiB

# Fixed total memory (no per-thread scaling)
fgumi filter --max-memory 4096 --memory-per-thread false

# Host-aware: detect (cgroup-aware) RAM and subtract --memory-reserve
fgumi filter --max-memory auto
fgumi filter --max-memory auto --memory-reserve 12GiB
```

This is the same `--max-memory` surface as `fgumi sort`. The default is opt-in:
the budget stays 768 MiB/thread unless you pass `auto`, so on a fixed-RAM host
(e.g. a 30 GiB container at `--threads 16`) prefer `--max-memory auto` or a fixed
total budget to avoid OOM.

### Memory Scaling Behavior

| Threads | Per-thread Mode | Fixed Mode |
|---------|----------------|------------|
| 1       | 768 MiB        | 768 MiB    |
| 4       | 3 GiB          | 768 MiB    |
| 8       | 6 GiB          | 768 MiB    |
| 16      | 12 GiB         | 768 MiB    |

With `--max-memory auto`, the total budget is instead `host RAM − reserve`
(divided across threads when per-thread), so it never scales past the host.

### Deprecated flags

`--queue-memory`, `--queue-memory-per-thread`, and `--queue-memory-limit-mb`
remain accepted as hidden aliases of `--max-memory` / `--memory-per-thread` for
backward compatibility, but new scripts should use the `--max-memory` flags.

## Compression Options

### Compression Level
- **Range**: 1 (fastest) to 12 (best compression)
- **Default**: 1 (fastest) for most commands; `fgumi merge` defaults to 6
- **Usage**: `--compression-level N`

### Compression Threading
- **Default**: Matches `--threads` setting
- **Override**: `--compression-threads N`
- **Best practice**: Usually leave at default

## I/O and Storage Tuning

For sequential workloads like BAM and FASTQ processing, I/O throughput is often the
bottleneck — not CPU. Two areas to check: OS readahead and volume throughput.

### OS Readahead

The Linux kernel prefetches file data into the page cache ahead of the application.
The default readahead window is typically 128 KB, which fgumi's decompression threads
can easily outpace. When that happens the processing thread stalls waiting on disk.

Check the current readahead (in 512-byte sectors):

```bash
blockdev --getra /dev/nvme1n1    # e.g. 256 = 128 KB
```

For sequential BAM/FASTQ workloads, increasing to 4 MB eliminates most I/O stalls:

```bash
# 4 MB = 8192 sectors (requires root)
sudo blockdev --setra 8192 /dev/nvme1n1
```

This setting does not persist across reboots. Add it to a startup script or udev rule
if needed.

### `--async-reader` (Experimental)

When you cannot tune OS readahead — containers, managed cloud instances, network
mounts — `--async-reader` provides a similar benefit from userspace. It spawns a
dedicated I/O thread that reads raw bytes into a bounded queue ahead of the
decompression step, so processing threads do not block on disk.

```bash
fgumi group \
  --async-reader \
  --threads 8 \
  --input reads.bam \
  --output grouped.bam
```

`--async-reader` works with all input types: BAM files, BGZF/gzip/plain FASTQs,
and piped stdin. It is supported by all commands that read BAM/FASTQ input,
including `sort`. It is most effective when I/O latency is high (network storage,
cold page cache, small OS readahead). On systems where you can already set 4 MB+
readahead, the additional benefit is modest.

### AWS EBS Volume Throughput

On AWS, `gp3` volumes default to 125 MB/s throughput regardless of size. For BAM
processing this is often the binding constraint. Increasing to 300-500 MB/s is
inexpensive and has a large impact:

```bash
# Increase throughput on an existing volume (takes effect within minutes)
aws ec2 modify-volume \
  --volume-id vol-0123456789abcdef0 \
  --throughput 500
```

For sustained sequential I/O, also consider increasing IOPS (default 3000) if your
reads are small. Monitor with `iostat -x 1` to confirm the volume is the bottleneck
before spending on higher provisioned throughput.

## Scenario-Based Configurations

### High-Throughput Server
**Goal**: Maximum processing speed for large datasets

```bash
fgumi filter \
  --threads 16 \
  --max-memory 1GB \
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
  --max-memory 512 \
  --memory-per-thread false \
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
  --max-memory 2GB \
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
  --async-reader \
  --threads 4 \
  --max-memory 512 \
  --compression-level 9 \
  --input dataset.bam \
  --output filtered.bam
```

**Rationale**:
- `--async-reader` hides network I/O latency (see [I/O and Storage Tuning](#io-and-storage-tuning))
- Moderate threading to avoid overwhelming network
- Conservative memory usage
- Maximum compression to reduce network transfer

### Development/Testing
**Goal**: Fast iteration with minimal resource usage

```bash
fgumi filter \
  --max-memory 256 \
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
- Adjust `--max-memory` if seeing swap activity (or pass `--max-memory auto`)

### Thread Utilization
- Use `htop` or similar to monitor CPU usage
- All threads should show activity during processing
- Consider reducing threads if not fully utilized

### I/O Patterns
- Monitor disk I/O with `iotop` or `iostat -x 1`
- If threads are idle waiting on I/O, increase OS readahead or try `--async-reader` (see [I/O and Storage Tuning](#io-and-storage-tuning))
- Network storage may benefit from lower thread counts
- SSD storage can handle higher thread counts

## Troubleshooting

### Out of Memory Errors
1. Pass `--max-memory auto` to size the budget to the host, or reduce `--max-memory`
2. Set `--memory-per-thread false` for a fixed total budget
3. Reduce `--threads`

### Poor Performance
1. Increase `--threads` if CPU usage is low
2. Increase `--max-memory` if I/O bound
3. Reduce `--compression-level` if CPU bound
4. Check OS readahead and EBS throughput if disk I/O is the bottleneck (see [I/O and Storage Tuning](#io-and-storage-tuning))

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
- Consider using `--memory-per-thread false` (or `--max-memory auto`)

## Command-Specific Considerations

### Extract
- Benefits from high memory (large FASTQ processing)
- Compression level affects output size significantly

### Zipper
- For best throughput, pipe uncompressed BAM from the aligner (e.g. `bwa-mem3 mem --bam=0`).
  Uncompressed BAM skips SAM text formatting on the aligner side and SAM parsing on the zipper
  side, and adds only ~26 bytes of BGZF framing per ~64 KiB block
- SAM input is fine for aligners that can't emit BAM; compressed BAM on a pipe wastes CPU on
  both ends for data the sort step will re-compress anyway
- The zipper pipeline uses raw-byte merging internally: aligned records are not fully decoded and
  re-encoded unless the record actually needs modification, which eliminates a significant CPU
  bottleneck on high-throughput runs

### Sort
- Uses an internal LoserTree (tournament tree) for k-way merging, which performs significantly
  better than a simple heap merge when the number of sorted runs is large
- `--max-memory` controls how much RAM is used for sort buffers; increase for large files to
  reduce the number of intermediate merge passes
- For template-coordinate sort with single-cell data, the `CB` tag is included automatically
- `--async-reader` is supported and can improve Phase 1 (input reading) throughput when disk
  latency is high or the OS page cache readahead is small

### Merge
- `fgumi merge` performs a k-way merge using a LoserTree for efficient multi-file merging
- Thread count (`--threads`) controls compression parallelism, not merge concurrency
- For template-coordinate merges with single-cell data, the `CB` tag is included automatically

### Group/Dedup
- Memory usage scales with UMI diversity and the number of reads at any given position
- Higher thread counts improve UMI processing
- The `--metrics PREFIX` flag writes all grouping metrics in one step with minimal overhead

### Simplex/Duplex Metrics
- Both `simplex-metrics` and `duplex-metrics` are single-threaded; they do not benefit from `--threads`
- Memory usage is proportional to the number of unique genomic positions in the input

### Consensus (Simplex/Duplex/CODEC)
- Memory proportional to family sizes
- Benefits from balanced threading and memory

### Filter
- Streaming operation benefits from pipeline memory
- Compression affects final output size

## Migration from Legacy Parameters

The `--queue-memory`, `--queue-memory-per-thread`, and `--queue-memory-limit-mb`
flags are deprecated (but still accepted as hidden aliases). Migrate to the
`--max-memory` family, which matches `fgumi sort`:

```bash
# Old (deprecated)
fgumi group --queue-memory-limit-mb 4096
fgumi group --queue-memory 4096 --queue-memory-per-thread false

# New (recommended)
fgumi group --max-memory 4096 --memory-per-thread false

# Or let fgumi size the budget to the host:
fgumi group --max-memory auto
```

The new parameters provide host-aware (`auto`) sizing and human-readable formats while maintaining backward compatibility.

# Performance Tuning Guide

fgumi provides three key options to optimize performance for your system: threading, memory management, and compression. This guide explains how to configure these options for different scenarios.

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
- **Range**: 1 (fastest) to 9 (best compression)
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

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
- **Behavior**: Command-dependent. Commands that special-case `ThreadingMode::SingleThreaded` run an optimized fast path with minimal pipeline overhead; others (e.g. `correct`) always run the declarative chain and simply execute it with a single worker. Either way, omitting `--threads` never runs more than one worker.
- **Best for**: Small files, memory-constrained systems, debugging

### Explicit Single-threaded Mode
- **Usage**: `--threads 1`
- **Behavior**: Uses the chain pipeline with a single worker thread — same pipeline as `--threads N` but with N=1; does **not** use the no-flag fast path
- **Best for**: Isolating pipeline behavior in a single-threaded context

### Multi-threaded Mode
- **Usage**: `--threads N` where N > 1
- **Behavior**: Uses the declarative chain pipeline with N worker threads (the same engine as `--threads 1`, scaled out); the legacy `--scheduler` flag is inert on the chain — it selects no dispatch policy
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

### What the Budget Bounds

`--max-memory` is a **total**: the BAM pipeline stops admitting input once the
bytes queued between its stages reach it, so a slow or contended output device
backs pressure up to the reader instead of filling every queue.

It is not a per-stage allowance. The reorder-buffered stages back off at
512 MiB and the processed queue at 256 MiB, and those marks do not rise with
the budget. A budget below one of them pulls it down; a budget above it leaves
it where it is. Raising `--max-memory` therefore raises how much the pipeline
may hold in total, not how much any one stage may hold, and the startup line
reports both numbers:

```
Queue memory budget: 6.0 GiB total (768.0 MiB/thread × 8 threads); per-stage high-water marks 512.0 MiB (reorder-buffered stages), 256.0 MiB (processed queue)
```

Note also that `fgumi extract`'s FASTQ pipeline shares this option and logs the
same line, but does not yet enforce the total — there only the per-stage marks
bind. That gap is tracked in
[fgumi#766](https://github.com/fulcrumgenomics/fgumi/issues/766).

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

When any step woke another thread or held an item, a second `wake` table follows the step table, one row per step:

- `notify` / `nowait` / `suppressed`: pool event-count notifies issued on the step's `Progress` (or, in a chain with a dedicated-thread step, on a dispatch that flushed a held item); how many of them found no parked worker; and how many found one but left it parked because the pool's concurrency ceiling was reached (0 unless a ceiling is set). `nowait` and `suppressed` are subsets of `notify`: `notify − nowait − suppressed` notifies woke a worker.
- `unpark`: direct wakes of a dedicated driver thread or a pinned pool worker that consumes the step's output, or of the step's own thread, for a same-thread consumer of a flushed held item.
- `reverse`: threads woken because this step's pops made room for an item they were holding.
- `fallback`: pool wakes that found no parked event-count waiter and woke a worker idling on its own timer instead: one that is pinned or holding an item, else one parked only because a phase cap refused it (when the woken step's own cap, if any, has a free permit).
- `gated_off`: output branches of a dedicated-thread step that pushed nothing on a `Progress` dispatch (or, in a chain with a dedicated-thread step, on a dispatch that flushed a held item), so no wake was sent for them; a dispatch counts once per such branch.
- `req_woken` / `req_awake` / `req_pending`: outcomes of a step's explicit request for a pool worker: another worker was woken; no other worker was parked where a request reaches it; or the request was dropped because an earlier one is still unanswered (at most one such event-count request is in flight per run).
- `holds` / `retries`: items whose push into one of the step's memory-bounded outputs was refused back to the step and later went through, and refused retries in between; `held_p50_ms` / `held_p99_ms` / `held_max_ms` are the wait from the first refusal to the push (p50/p99 are histogram upper bounds). On an ordered output a refusal is the queue's or the reorder buffer's cap; an item the reorder buffer accepted while the queue was full is not a hold, and a held item it accepted ends its hold. A parallel step's holds are timed per worker; any other step's per step, so a held item another worker pushes ends the hold.

Each pool worker row also reports `waits notified% timed_out%`: how many times it waited on the shared event-count, the share of those waits a peer ended while the worker was blocked (`notified%`), and the share its own timer ended (`timed_out%`). The rest of its waits never blocked: work was published between the worker arming its wait and blocking, so it re-polled at once. A worker that idles on its own timer instead (a pinned worker, the only worker, one holding an item, or one a phase cap refused, which that cap's release wakes) adds `timer_parks=N (unparked U, timed_out T)`. Each dedicated-thread (`detached`) line reports its parks as `N parks (unparked U, timed out T, avg …µs)`; a high timed-out share means the thread is being fed by its timer rather than by wakes.

### Sort Statistics

```bash
fgumi sort --sort-stats --input reads.bam --output sorted.bam
```

`fgumi sort` runs through the same shared `ChainBuilder` pipeline as every other command, but it carries a dedicated `--sort-stats` flag instead of `--pipeline-stats` for its own merge-loop diagnostic. Whenever the k-way merge runs -- any sort that spills, or a no-spill sort that still holds more than one in-memory chunk -- it prints a `Sort merge diag: stalls=... contention=... output_full=... progress_dispatches=...` line reporting merge-loop stalls (waiting on decompress), contention (dispatches that produced nothing), and output backpressure. `stalls=` counts only stalls the merge could not resolve at once: when the merge finds the block it is waiting on has landed by the time it registers to be woken for it, it keeps merging and the stall is not counted, so compare `stalls=` only between runs of the same release.

When the merge reads spill files, the diag line is followed by the merge-demand lines, in this order:

- `Merge demand:` stall episodes (one per wait for a block, however many times the merge re-checks within it) and their exact total time; then `wakes delivered W of P parking registrations (R registrations)`: of the `R` times the merge registered to be woken for the spill file it waits on, `P` were followed by a park, and `W` parks were ended by the delivery of that file's block (the rest ended on the merge's idle timer).
- `Awaited slot at stall:` the share of stall episodes whose awaited spill file had no block being decompressed (`starved`) or had blocks being decompressed (`decompressing`). Under `--file-granularity` decompression runs inline and is never tracked per block, so every stall reports `starved` and the line says the stalls are not classified.
- `Pool at stall:` the merge asks the pool for one sleeping worker once per stall episode, before it parks; the outcomes of those requests (`woken`, `all-awake`: no worker was asleep where a request reaches it, `pending`: an earlier request was still unanswered, `unavailable`: no pool to ask), the share of requests that found a sleeping worker (parked on the pool's shared wait, or idling on its own timer), and the stall time of the episodes whose request did. The same outcomes, per step, are the `req_*` columns of the `--pipeline-stats` wake table.
- `Merge output:` how many output batches the merge flushed early, short of their target size, because it stalled with records already merged.

Only when the sort spills nothing *and* fits in a single in-memory chunk does no k-way merge run; there it instead prints one `Sort fast-path diag: ...` line noting the single-chunk in-memory fast path was taken. Off by default; it is instrumentation for performance work, read from a log with a `grep`.

### Scheduler Strategy (legacy, inert)

The `--scheduler` flag selects a legacy scheduler strategy that the typed-step chain engine does **not** consume, so it has no effect on any chain-backed command. It is retained for backward compatibility only; setting it to a non-default value logs a warning that the requested strategy is ignored. The strategy names still parse, and the common ones once meant:

| Strategy | Former behavior |
|----------|-----------------|
| `balanced-chase-drain` | Default. Balanced work distribution with output drain mode. |
| `fixed-priority` | Static thread roles (reader, writer, workers). Simple baseline. |
| `chase-bottleneck` | Threads dynamically follow work through the pipeline. |

Additional legacy strategy names (`thompson-sampling`, `ucb`, `epsilon-greedy`, and others) also parse but are equally inert.

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
- `--max-temp-files` sets how many spilled runs may be live at once; when the limit is reached,
  adjacent runs are merged — the smallest first — until the count is back under it. The final k-way merge opens every remaining run at once, so
  this limit is what bounds the sort's open file descriptors — and consolidation rewrites data
  that is already sorted, making it pure overhead whenever the descriptor budget could have
  carried the runs. The default, `auto`, sizes the limit to the process's soft open-file limit
  (`ulimit -n`), less a reserve for the input, output and index handles, and capped at a tested
  maximum. On a host with a low `ulimit -n`, raising it lets the sort avoid consolidation
  entirely. Every spilling sort reports `Merge sources:` in its summary (equal to `Spill runs:`
  unless it consolidated); a sort that consolidated also reports `Consolidations:`, and its
  phase-timing roll-up has a consolidation bucket. Pass
  an explicit value (`--max-temp-files 64`) to pin it instead — if that value exceeds the
  open-file budget, the sort says so at startup rather than failing partway through with "Too
  many open files". Must be at least 2; the output is unchanged either way
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

# fgumi Memory Investigation Prompt

## Overview

During benchmarking of fgumi commit `b698874654e15f589776f9ea379fa6160f2d112d`, we discovered a significant discrepancy between the internal memory tracking and actual memory usage reported by the OS. This prompt documents the issue and provides guidance for investigation.

## The Problem

### Memory Tracking Discrepancy

For the `group` command on `synthetic-pipeline-xlarge` (90M records, 16.8M unique molecules):

| Metric | Reported Value |
|--------|----------------|
| Pipeline internal tracking | 2,265 MB |
| Peak Q3 Reorder Buffer | 598 MB |
| Queue memory limit | 4,096 MB |
| **OS maximum resident set size** | **20,600 MB** |
| OS peak memory footprint | 11,400 MB |

**The internal tracking shows ~2.3 GB but the OS reports 20.6 GB - an 8-9x discrepancy.**

## Root Cause Analysis (COMPLETED)

### Investigation Summary

The 18+ GB of untracked memory comes primarily from **parallel `process_fn` executions** in the Process step. With 8 threads, each thread independently creates:

1. **UMI Assigner data structures** (`AdjacencyUmiAssigner`):
   - `AHashMap<Umi, usize>` for counting UMIs
   - `Vec<(Umi, usize)>` for sorted UMI list
   - Graph node structures for adjacency
   - Children vectors for graph traversal

2. **MI grouping HashMap** (`by_mi: AHashMap<String, Vec<Template>>`):
   - Groups templates by assigned molecule ID
   - 16.8M unique MIs × ~80 bytes overhead

3. **String allocations**:
   - MI tag strings created with `format!`
   - MI string pool (`mi_pool: AHashMap<&str, Arc<BString>>`)

### Memory Math

For the xlarge dataset (16.8M unique molecules, ~18.7M templates after filtering):

| Component | Per-Thread | 8 Threads |
|-----------|------------|-----------|
| UMI count HashMap | ~1.3 GB | ~10 GB |
| Sorted UMI vector | ~750 MB | ~6 GB |
| Graph nodes | ~670 MB | ~5.4 GB |
| by_mi HashMap | ~1.3 GB | ~10 GB |
| **Total (overlapping)** | **~2-3 GB** | **~16-24 GB** |

Note: Not all threads hit peak memory simultaneously, but with large position groups,
multiple threads can be processing multi-million template groups concurrently.

### Why Pipeline Tracking Misses This

The pipeline tracks memory at queue boundaries:
- `q3_reorder_heap_bytes`: Before Group step
- `q4_groups_heap_bytes`: After Group, before Process
- `q5_processed_heap_bytes`: After Process, before Serialize

But allocations **inside** `process_fn` are never tracked:
```rust
// In group.rs process_fn - ALL OF THIS IS UNTRACKED:
let assigner = strategy.new_assigner_full(effective_edits, 1, index_threshold);
assign_umi_groups_impl(&mut templates, assigner.as_ref(), ...);  // Creates HashMaps
let mut by_mi: AHashMap<String, Vec<Template>> = AHashMap::with_capacity(estimated_groups);
// ... grouping and sorting ...
```

## Solution Requirements

**Hard constraints:**
1. **Memory target: ≤ 8 GB peak** (currently 20.6 GB)
2. **NO runtime regression** (currently 62 seconds for xlarge)

## Implementation Plan

### Phase 1: Process Step Memory Tracking (Optional)

Track allocations during process_fn and apply backpressure:

```rust
// Add to process_fn signature:
fn process_fn(group: G, memory_tracker: &AtomicU64) -> P {
    // Track UMI assigner memory
    let assigner_mem = estimate_assigner_memory(group.templates.len());
    memory_tracker.fetch_add(assigner_mem, Ordering::Relaxed);

    // ... do work ...

    // Track by_mi HashMap memory
    let by_mi_mem = estimate_hashmap_memory(by_mi.len());
    memory_tracker.fetch_add(by_mi_mem, Ordering::Relaxed);

    // ... finish and release ...
    memory_tracker.fetch_sub(total_mem, Ordering::Relaxed);
}

// In try_step_process - check before starting:
if process_memory_tracker.load() > PROCESS_MEMORY_LIMIT {
    return StepResult::OutputFull;  // Apply backpressure
}
```

### Phase 2: Efficient UMI Data Structures

**Goal:** Reduce per-thread memory for UMI assignment.

#### 2a: Bit-Encoded UMIs

Replace `String` UMIs with compact `BitEnc` representation:

```rust
// Current: Umi = String (24 bytes + heap allocation per UMI)
// Proposed: Umi = BitEnc (8 bytes inline for UMIs up to 32bp)

// In assigner.rs:
pub struct CompactUmi(u64);  // 2 bits per base, up to 32bp

impl CompactUmi {
    fn from_bytes(bytes: &[u8]) -> Option<Self> { ... }
    fn hamming_distance(&self, other: &Self) -> u32 { ... }
}
```

**Memory savings:** ~60-70% reduction in UMI storage

#### 2b: Reusable Assigner Pools

Pool UMI assigner data structures instead of reallocating:

```rust
// Thread-local assigner with reusable buffers:
thread_local! {
    static ASSIGNER_BUFFERS: RefCell<AssignerBuffers> = RefCell::new(AssignerBuffers::new());
}

struct AssignerBuffers {
    counts: AHashMap<CompactUmi, usize>,
    sorted_umis: Vec<(CompactUmi, usize)>,
    nodes: Vec<Node>,
    // Clear and reuse between position groups
}
```

## Implementation Order

| Priority | Phase | Memory Impact | Runtime Impact | Effort |
|----------|-------|---------------|----------------|--------|
| 1 | **2a: BitEnc UMIs** | ~60-70% UMI reduction | **5-15% faster** (better cache, less alloc) | Medium |
| 2 | **2b: Assigner pools** | ~20-30% further reduction | **5-10% faster** (reduced allocations) | Medium |
| 3 | **1b: Process tracking** | Better visibility | Neutral (<1% overhead) | Low |

**Key insight:** Phases 2a and 2b actually *improve* runtime by reducing allocation overhead and improving cache locality. The net effect of implementing 2a + 2b should be **memory reduced to ~6-8 GB with ~5-10% faster runtime**.

## Validation Plan

### Memory Testing

```bash
# Run with memory tracking
/usr/bin/time -l target/release/fgumi group \
    --input $BENCHMARKS/data/processed/synthetic-pipeline-xlarge.sorted.bam \
    --output /tmp/test-group.bam \
    --strategy Identity \
    --edits 1 \
    --threads 8 \
    --queue-memory-limit 4096 \
    --pipeline-stats \
    2>&1 | tee /tmp/group-memory-test.log

# Verify max RSS ≤ 8 GB (8589934592 bytes)
grep "maximum resident set size" /tmp/group-memory-test.log
```

### Performance Testing

```bash
# Must NOT regress from baseline 62 seconds
hyperfine --warmup 1 --runs 3 \
    'target/release/fgumi group --input $BENCHMARKS/data/processed/synthetic-pipeline-xlarge.sorted.bam --output /dev/null --strategy Identity --edits 1 --threads 8'
```

### Correctness Testing

```bash
# Output must be identical to baseline
diff <(samtools view baseline.bam | sort) <(samtools view new.bam | sort)
```

## Key Files to Modify

| File | Changes |
|------|---------|
| `src/lib/unified_pipeline/bam.rs` | Add process memory semaphore/tracking |
| `src/commands/group.rs` | Integrate memory tracking in process_fn |
| `src/lib/umi/assigner.rs` | Add CompactUmi, pooled buffers |
| `src/lib/unified_pipeline/base.rs` | Add memory tracking to ProcessPipelineState |

## Success Criteria

- [ ] Peak RSS ≤ 8 GB on xlarge dataset
- [ ] Runtime ≤ 62 seconds (no regression)
- [ ] All existing tests pass
- [ ] Output byte-for-byte identical to baseline

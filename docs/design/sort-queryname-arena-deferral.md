# Queryname sort: arena unification deferred

## Status: DEFERRED (not a feat-runall blocker)

The in-memory Phase-1 buffer for the **queryname** sort orders stores
`Vec<(K, RawRecord)>` — one owned record-byte allocation *and* one owned key
allocation per record. The **coordinate** and **template-coordinate** orders do
not: the inline-buffer work (#432) copies each record's bytes once into a shared
segmented byte arena (`RecordBuffer` / `SegmentedBuf` in `inline.rs`) and stores
only `RecordRef { sort_key, offset, len }` entries, which are radix/comparison
-sorted and then written by slicing the arena. The owned `RecordSource::Iterator`
`Direct` arm (`read_ahead.rs`) — and its per-record `RawRecord` allocation — is
used **only** by the queryname path.

This is the substance of finding **S3-015**. The reviewer's suggested fix
(hold a scratch `RawRecord` in the `Direct` variant and `mem::take` it out per
`next()`) is a **no-op**: `mem::take` hands the caller the scratch's buffer and
leaves an empty (capacity-0) `Vec` behind, so the next `read_raw_record`
`resize` reallocates exactly as `RawRecord::default()` does today. No reuse is
possible on owned iteration because `sort_queryname_keyed` *retains* every
yielded record in its `entries` buffer until spill — each live record genuinely
needs its own buffer.

## Why the queryname path is the holdout

Coordinate/template keys are fixed-width packed integers (`u64`/`u96`/…), so the
`RecordRef`-style entry carries the key inline and the arena holds bytes only.
Queryname keys are **variable-length byte strings** (the read name), so the
naive entry can't be a fixed-width field — which is why this path fell back to
owning both the record and the key.

## The principled solution

Extend the inline arena to the queryname path so all three orders share one
zero-per-record-allocation representation:

1. **Record bytes** → the shared `SegmentedBuf` arena (offset/len ref), exactly
   like coordinate/template. Records are copied in once via
   `RecordSource::next_record_borrowed` (the borrow-in-place path), so the owned
   `Iterator::next` `Direct` arm becomes dead and can be removed.
2. **Keys** → one of:
   - a parallel **key-arena**: name bytes pooled with offset/len refs (no
     per-key `Vec`), entries become `(KeyRef, RecordRef)`; or
   - a **small inline key**: a name hash + tiebreak stored inline, with the full
     name re-derived from the record slice (via the arena) at comparison time
     for the rare hash-collision path.

The sort then sorts `(key, ref)` entries; spill and final write slice the arena.
This collapses coordinate, template-coordinate, and queryname onto one
`KeyedRecordBuffer<K>`-style abstraction.

## Prerequisites / why it is deferred

- A `KeyedRecordBuffer<K>` (or generalization of `RecordBuffer`) that holds a
  variable-length key alongside each `RecordRef`, plus the key-arena or
  inline-key comparator.
- The spill writer (`PooledChunkWriter`) and merge packing (`MemorySources`)
  must consume `(key, range)` slices instead of owned `RawRecord`s on the
  queryname path (the coordinate/template paths already do the arena form).
- A **sort bench matrix** (queryname-lex and queryname-natural, in-memory and
  multi-spill, single- and multi-threaded) to confirm the arena form does not
  regress comparison-sort throughput — the keyed comparison is hotter than the
  packed-integer radix, so the key representation must be benched, not assumed.

Until those land, the per-record allocation on the queryname path is an
accepted, documented tradeoff (see the `sort_queryname_keyed` comment noting the
owned-`RawRecord` accumulation, and the `RecordSource::Iterator` `Direct` arm
PERF NOTE). It is correctness-neutral; only allocation churn is at stake, which
`force_mi_collect()` after each spill already mitigates.

# Sort Phase-2 merge: two drivers, unification deferred

## Status: DEFERRED (not a feat-runall blocker)

fgumi's external sort has **two** Phase-2 (merge) drivers, and they exist by
design:

- **Pool / standalone** (`crates/fgumi-sort/src/worker_pool.rs`,
  `external.rs::merge_chunks_generic`) — drives standalone file-to-file
  `fgumi sort`. Work-stealing pool reads + decompresses spill blocks in
  parallel; the merge consumer is a **parked OS thread**
  (`advance_to_next_block`). Keeps `--write-index`, `--verify`, zstd spill.
- **Streaming / runall** (`crates/fgumi-sort/src/merge_slots.rs`,
  `external.rs::MergeDriver`) — drives the fused in-pipeline `runall` sort. A
  cooperative `try_run` state machine that yields (`Stalled`/`WouldBlock`)
  instead of blocking, so it composes with the typed-step pipeline framework.

Unifying these onto a single driver is a deliberate **future** refactor with its
own design + bench gate + deadlock re-analysis. It is intentionally NOT done in
the feat-runall line.

## What is already shared (so the duplication is smaller than it looks)

- `LoserTree<K>` (`loser_tree.rs`) — the actual k-way merge / tie-break
  primitive. Both drivers use the one implementation.
- `RawSortKey` comparators and the `keys.rs` key types.
- `ReorderBuffer` (`fgumi-bam-io`) — also used by Phase-1 read-ahead, so it is a
  shared utility, not pool-specific duplication.

What genuinely differs is the **driver** (the outer drive loop and the
read/decompress topology), not the merge core.

## Why the streaming path dropped the pool's machinery

The pre-v4 streaming slot also carried a `raw_blocks` FIFO, a `decomp_in_flight`
counter, a reorder buffer, and a cap+gap-filler admission rule. That design
**deadlocked at production scale**: the framework's drain protocol could Skip
workers during a transient "all slots at cap simultaneously" window. v4 collapsed
read-and-decompress into one inline op per worker per slot, which made a plain
bounded FIFO correct (blocks decompress strictly in read order) and removed the
deadlock surface.

The pool path keeps all of that machinery and is **immune** to the same
deadlock, because (a) its consumer is a parked OS thread the engine cannot Skip,
and (b) it is single-*reader* but **multi-*decompressor*** — workers decompress
in parallel, so completion order ≠ pop order, and the reorder buffer +
in-flight counter are a *correctness* dependency (they re-order out-of-order
completions and close the `is_drained` truncation race), not redundancy.

## Prerequisites before unification is even attemptable

1. The streaming path gains a true file-to-file drive (BAM source +
   bgzf-decompress + parse → SortAndSpill → SortSpillDecompress → SortMerge →
   serialize → WriteBgzfFile).
2. `--write-index` (`.bai`) wired onto/alongside that chain (the pool path has
   it via `IndexBamFinalizeHook`, PR #393; the streaming step does not).
3. `--verify` ported (today a separate read-and-check path, designed to stay
   legacy).
4. A deadlock-safe consumer for the file-to-file case — the streaming consumer
   is the cooperative `try_run` model that hit the v4 production deadlock.

## Explicit risk

Routing standalone sort through the streaming path **re-opens the v4
cooperative-consumer deadlock** that standalone sort is currently immune to
(parked OS thread). Until all four prerequisites land, it is a net capability
LOSS + added risk and buys no shipped feature.

## Safe-to-pursue-independently sub-step (does NOT require the above)

Factor the shared **per-winner merge-emit core**: both drivers already share
`LoserTree`, so extract a `MergeSource<K>` trait (`try_next_record`) + a single
`merge_emit(tree, sources, emit_fn)` core that the pool's blocking loop and the
streaming state machine both call one-winner-at-a-time. Keep the
blocking-vs-cooperative **outer** drive separate — that is the deadlock-relevant
axis and must NOT be merged. MED risk (hot-loop perf; needs a sort bench
matrix), but carries no deadlock surface.

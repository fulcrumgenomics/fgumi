//! `ReadBlocks` serial arena-admit step, `InflateToArena` parallel step, and
//! `FindBoundariesAndSort` serial step.
//!
//! `ReadBlocks`: serial step that consumes raw `BgzfBlock`s, uses an [`ArenaPool`]
//! (capacity 1) to acquire/reuse arenas, reserves a fixed front region
//! `[0, FRONT_REGION)` per run for the straddler carry (Task 2), reserves an
//! uninit slot per block via `grow_uninit` at offsets `>= FRONT_REGION`, and
//! seals a run when its cumulative uncompressed bytes reach `run_cap = memory_limit`.
//! Sealed mid-stream runs emit [`ArenaBlock`]s with `seals_to_spill = true`;
//! the final residual run emits with `seals_to_spill = false`.
//! This is the second `unsafe` site in `fgumi-pipeline-io`; see CLAUDE.md
//! §"Approved hot-path unsafe (parallel-inflate sort ingest, fgumi-pipeline-io)"
//! for the full justification.
//!
//! `InflateToArena`: each worker decompresses one BGZF block directly into its
//! disjoint arena slot via [`fgumi_bgzf::decompress_into_slice`].  This is the first
//! `unsafe` site in `fgumi-pipeline-io`; see CLAUDE.md §"Approved hot-path
//! unsafe (parallel-inflate sort ingest, fgumi-pipeline-io)" for the
//! full justification and SAFETY invariant.
//!
//! `FindBoundariesAndSort`: serial step that scans a run's contiguous arena span
//! directly (no copy), skips the BAM header via [`crate::boundaries::bam_header_len`],
//! builds `(body_offset, block_size)` refs, sorts via
//! [`fgumi_sort::coordinate_chunk_from_arena_refs`], and emits
//! `SortChunkEvent::Residual` + `AllAnnounced`.  Single-run / no-spill scope;
//! multi-run seal and straddler handling are 3b.4.
//!
//! Together these three steps form the BAM sort ingest front for the block-input
//! path.  The chain builder's `add_sort` (`src/lib/pipeline/chains/builder.rs`)
//! wires them as `ReadBlocks → InflateToArena → FindBoundariesAndSort` (feeding
//! the Phase-1 spill/merge tail) when the sort source is raw BGZF blocks.

use std::collections::VecDeque;
use std::io;
use std::sync::Arc;

use fgumi_bgzf::{Decompressor, decompress_into_slice};
use fgumi_sort::{
    ArenaPool, InMemoryChunk, PooledSegmentedBuf, RawSortKey, RecordRef,
    coordinate_chunk_from_refs, extract_coordinate_key_inline, queryname_chunk_from_arena_refs,
};

use crate::types::BgzfBlock;

use fgumi_pipeline_core::held::HeldSlot;
use fgumi_pipeline_core::item::{HeapSize, Ordered};
use fgumi_pipeline_core::outputs::{OrderedBytesSingle, Single};
use fgumi_pipeline_core::queues::QueueSpec;
use fgumi_pipeline_core::reorder::BranchOrdering;
use fgumi_pipeline_core::step::{DetachedGroup, Step, StepCtx, StepKind, StepOutcome, StepProfile};
use fgumi_pipeline_core::{HeldRetry, Unpushed};

use crate::boundaries::bam_header_len;
use crate::sort::protocol::{MemoryChunkErased, SortChunkEvent};

// ============================================================================
// Constants
// ============================================================================

/// Reserved front region at the start of every run's arena, in bytes.
///
/// Block 0 always lands at offset `FRONT_REGION`.  The region `[0, FRONT_REGION)`
/// is left uninitialized by `ReadBlocks` and filled by `FindBoundariesAndSort`
/// when it writes a carried straddler from the previous run (Task 2, 3b.4).
/// Must be larger than any single BAM record; 8 MiB covers all realistic records.
pub const FRONT_REGION: usize = 8 * 1024 * 1024;

/// Maximum uncompressed size of a single BGZF block: the BGZF spec bounds an
/// uncompressed block at 2^16 = 65536 bytes.  Used to cap a run's seal budget
/// (`run_cap`) so the single-segment-per-run invariant
/// `FRONT_REGION + Σ(block ISIZE) ≤ segment_size` holds with a true upper bound
/// of margin (see `ReadBlocks::new`).  This must be `>=` any real block's ISIZE;
/// 65536 is the exact spec maximum, so the headroom is a genuine bound rather
/// than relying on the `< run_cap` cumulative check being off by one.
const MAX_BGZF_BLOCK: usize = 1 << 16;

/// Bytes ahead of the current scan cursor to software-prefetch in the
/// `FindBoundariesAndSort` boundary+key scan.  Chosen from a microbench over a
/// cold ~2.6 GiB arena (matching the ~220 B/record production density): 2 KiB
/// gave the best speedup (~15%); ≤1 KiB was negligible (too little lead time),
/// 4 KiB matched 2 KiB.  At ~220 B/record this is ~9 records of lead.
const SCAN_PREFETCH_DISTANCE: usize = 2048;

/// Max blocks `ReadBlocks` admits per `try_run` dispatch.
///
/// `ReadBlocks` runs on the coordination driver (`StepKind::Detached`), whose
/// drain-first `round_robin_dispatch` restarts the whole downstream walk on every
/// `Progress` (see `runtime::driver`): admitting one block per dispatch made this
/// serial admit the input-side throughput bottleneck — the raw-block reader
/// upstream backed up while the parallel `InflateToArena` pool workers starved
/// (empty pops), leaving the pool short of full CPU occupancy. Admitting a bounded
/// batch per dispatch amortises the per-dispatch queue/reorder overhead so the
/// admit keeps the inflaters saturated, while the cap bounds how long the driver
/// dwells on admit (fairness vs the group's other steps) and how far it runs
/// ahead of the byte/count-bounded output queue (which backpressures early anyway).
const ADMIT_BATCH: usize = 64;

/// Software-prefetch (read, into L1, temporal) the cache line containing `byte`.
/// The `FindBoundariesAndSort` scan walks the run's arena cold (it was written by
/// the parallel `InflateToArena` workers long before this serial scan runs), so it
/// is latency-bound on cache misses; prefetching [`SCAN_PREFETCH_DISTANCE`] ahead
/// hides them.  `cfg`-gated to the supported architectures; a no-op elsewhere.
///
/// SAFETY note: this is the only place `fgumi-pipeline-io` uses an architecture
/// intrinsic.  Both `prfm` (`aarch64`) and `_mm_prefetch` (`x86_64`) are
/// *non-faulting hints* — they never read or write observable memory and never
/// trap, even on an unmapped address.  `byte` is a live `&u8` (the caller
/// bounds-checks the index), so the pointer is valid to name.  See CLAUDE.md
/// §"Approved hot-path unsafe (parallel-inflate sort ingest, fgumi-pipeline-io)".
#[inline]
fn prefetch_read_l1(byte: &u8) {
    let ptr: *const u8 = byte;
    #[cfg(target_arch = "aarch64")]
    #[allow(unsafe_code)]
    // SAFETY: `prfm pldl1keep` is a non-faulting prefetch hint over a valid pointer.
    unsafe {
        core::arch::asm!(
            "prfm pldl1keep, [{p}]",
            p = in(reg) ptr,
            options(nostack, readonly, preserves_flags),
        );
    }
    #[cfg(target_arch = "x86_64")]
    #[allow(unsafe_code)]
    // SAFETY: `_mm_prefetch` is a non-faulting prefetch hint over a valid pointer.
    unsafe {
        core::arch::x86_64::_mm_prefetch::<{ core::arch::x86_64::_MM_HINT_T0 }>(ptr.cast());
    }
    #[cfg(not(any(target_arch = "aarch64", target_arch = "x86_64")))]
    {
        let _ = ptr; // no portable stable prefetch; hint is a no-op on other arches
    }
}

// ============================================================================
// ReadBlocks: serial arena-admit step
// ============================================================================

/// `Serial` step that consumes raw [`BgzfBlock`]s, acquires arenas from an
/// [`ArenaPool`] (capacity 1), grows the whole arena segment to full capacity
/// once at acquire (so block slots can be sliced by inflate workers without
/// further `&mut` growth), reserves a front region `[0, FRONT_REGION)` per run,
/// assigns each block a slot at an arithmetic offset (`FRONT_REGION` + prefix-sum
/// of ISIZE), seals a run when cumulative uncompressed bytes reach `memory_limit`,
/// and emits one [`ArenaBlock`] per block.
///
/// **Eager (streaming) emission.** Blocks are emitted *as they are read*, not
/// buffered until the run seals: the `Arc<PooledSegmentedBuf>` is created at
/// acquire and the slot offset is known immediately (from the block's ISIZE,
/// without inflating), so an inflate worker can start on block 0 while
/// `ReadBlocks` is still reading block 5000 — the Read‖Inflate overlap.  Exactly
/// ONE block is withheld at a time (`last_block`): `is_last_of_run` and
/// `seals_to_spill` are run-level facts known only when the run seals, so the
/// most-recently-admitted block is held until either a successor in the same run
/// arrives (confirming it is *not* last → emit it tagged `false`) or the run
/// seals (it *is* last → stamp `is_last_of_run`/`seals_to_spill` and emit).  This
/// single held block coincides with the deferred-seal block, so the two are one
/// mechanism.
///
/// Mid-stream seals set `seals_to_spill = true`; the final (residual) run sets
/// `seals_to_spill = false`.  Only the LAST block of a run carries a meaningful
/// `is_last_of_run`/`seals_to_spill`; downstream (`FindBoundariesAndSort`) reads
/// those fields only on the `is_last_of_run` block, so non-last blocks emit
/// `is_last_of_run = false, seals_to_spill = false`.
///
/// With pool capacity 1, the next run's arena cannot be acquired until the
/// prior run's [`Arc<PooledSegmentedBuf>`] drops (i.e. after `CompressSpill`
/// consumes and releases it).  If acquisition fails, `try_run` returns
/// `NoProgress` (backpressure) and retries on the next call.
pub struct ReadBlocks {
    /// Bounded pool (capacity 1) that owns arena storage and reuses it across runs.
    pool: Arc<ArenaPool>,
    /// Size of each arena segment (`FRONT_REGION + run_cap + MAX_BGZF_BLOCK`); the
    /// whole segment is grown live once at acquire.
    segment_size: usize,
    /// `Arc`-shared arena for the current run, grown to full capacity at acquire so
    /// blocks (and their inflate slots) can be handed out eagerly without further
    /// `&mut` growth; `None` until the first block of the current run is admitted,
    /// and again after the run seals (dropping this step's handle so the arena can
    /// return to the pool once all emitted clones release it).
    arena: Option<Arc<PooledSegmentedBuf>>,
    /// Next slot offset within the current run's arena (set to `FRONT_REGION` on a
    /// fresh run, then advanced by each block's ISIZE).
    run_offset: u64,
    /// The most-recently-admitted block of the current run, withheld from emission
    /// until we learn whether another block joins this run, so `is_last_of_run` /
    /// `seals_to_spill` (run-level facts known only at seal) can be stamped on it.
    last_block: Option<ArenaBlock>,
    /// A block popped from the input but not yet admitted because a deferred
    /// seal fired first and the pool is momentarily exhausted (the just-sealed
    /// run's arena has not returned yet).  Re-admitted on a later `try_run` once
    /// the pool yields an arena.  Holding it here prevents losing the block.
    deferred_block: Option<BgzfBlock>,
    /// Emitted `ArenaBlock`s ready to push to the output queue.
    emit: VecDeque<ArenaBlock>,
    /// Held output slot for the backpressure retry path.
    held: HeldSlot<Unpushed<ArenaBlock>>,
    /// Monotonically-increasing ordinal assigned to each admitted block (across all runs).
    next_ordinal: u64,
    /// Current run sequence number.
    run_seq: u32,
    /// Cumulative uncompressed bytes in the current run (resets on each seal).
    run_cumulative: usize,
    /// Budget per run: seal when `run_cumulative >= run_cap`.
    run_cap: usize,
    /// Armed when the current run reached `run_cap`, but the seal is DEFERRED
    /// until the next block arrives.  Deferring by one block guarantees that a
    /// spill seal is only ever fired when at least one more block follows (so
    /// the straddler carry from a spilled run always has a subsequent run to
    /// complete it).  If the input drains while a seal is armed, the run is
    /// sealed as the RESIDUAL instead (the final run's last record always
    /// completes at EOF, so the carry is empty).  Without this, the budget seal
    /// could fire on the very last input block, emitting a trailing *spill*
    /// whose carry has no following run — surfacing as a spurious "truncated
    /// BAM" error.
    seal_pending: bool,
    /// Set to `true` after the final run has been emitted.
    finished: bool,
    /// Byte limit for the output queue.
    output_byte_limit: u64,
}

impl ReadBlocks {
    /// Create a new `ReadBlocks` step.
    ///
    /// - `memory_limit`: bytes of record data a run holds before sealing — the
    ///   FULL in-memory budget (`--max-memory × threads`, via
    ///   `resolve_memory_budget`).  This is the SAME spill trigger legacy uses:
    ///   the front spills to disk ONLY when the data exceeds this budget.  For
    ///   data that fits the budget the result is exactly ONE in-memory run, ZERO
    ///   disk spills, and one radix sort over the full budget — legacy's
    ///   in-memory algorithm and footprint, but with the no-copy arena ingest.
    /// - `output_byte_limit`: byte-bounds the output queue.
    ///
    /// The arena segment is sized so one full run always fits within ONE segment
    /// (`FRONT_REGION` + the run budget + one deferred-seal tipping block), which
    /// keeps the `FindBoundariesAndSort` contiguous-slice scan valid WITHOUT
    /// capping the run below the budget.  (The earlier design fixed the segment at
    /// 256 MiB and capped `run_cap` to fit it, which forced the front to spill
    /// even when the data fit `--max-memory` — a wall-clock regression vs legacy.
    /// Sizing the segment to the budget instead removes that forced spill.)
    #[must_use]
    pub fn new(memory_limit: usize, output_byte_limit: u64) -> Self {
        // Seal a run only when its record data reaches the full budget (legacy's
        // trigger); never cap below it.
        let run_cap = memory_limit;
        // One run = front region + up to `run_cap` of data + the single block that
        // tips `run_cumulative` over `run_cap` (the deferred seal holds the NEXT
        // block for the following run, so the current run overshoots by at most one
        // block).  Sizing the segment to that upper bound guarantees a run is
        // gap-free within one segment, so the scan's single `arena.slice(..)` over
        // `[scan_start, run_end)` never spans a segment boundary.
        let segment_size = FRONT_REGION + run_cap + MAX_BGZF_BLOCK;
        let pool = ArenaPool::new(1, segment_size);
        Self {
            pool,
            segment_size,
            arena: None,
            run_offset: 0,
            last_block: None,
            deferred_block: None,
            emit: VecDeque::new(),
            held: HeldSlot::new(),
            next_ordinal: 0,
            run_seq: 0,
            run_cumulative: 0,
            run_cap,
            seal_pending: false,
            finished: false,
            output_byte_limit,
        }
    }

    /// Ensure the current run's arena is acquired, grown to full capacity, and
    /// `Arc`-shared, with `run_offset` reset to `FRONT_REGION`.
    ///
    /// Returns `true` if the arena is ready, `false` if the pool is exhausted
    /// (backpressure: the previous run's chunk has not been consumed yet).
    ///
    /// # Side effects
    ///
    /// Acquires an arena from the pool, calls `reserve_full_capacity`, then grows
    /// the ENTIRE segment to `segment_size` live bytes via a single `unsafe`
    /// `grow_uninit` call (see the `// SAFETY:` comment inside), and wraps it in an
    /// `Arc<PooledSegmentedBuf>`.  This `grow_uninit` is the second `unsafe` site in
    /// `fgumi-pipeline-io`; see CLAUDE.md §"Approved hot-path unsafe (parallel-inflate
    /// sort ingest, fgumi-pipeline-io)" — the `ReadBlocks::ensure_arena` grow-once
    /// bullet — for the full justification.  Growing once, before sharing, is the
    /// soundness keystone: every block's slot offset is then computed arithmetically
    /// (`FRONT_REGION` + prefix-sum of ISIZE) so no further `&mut` growth is needed
    /// while inflate workers hold disjoint `slice_mut` views.
    fn ensure_arena(&mut self) -> bool {
        if self.arena.is_some() {
            return true;
        }
        let Some(mut arena) = self.pool.try_acquire() else {
            return false;
        };
        arena.reserve_full_capacity();
        // Grow the whole segment to `segment_size` live bytes in ONE call, here on
        // the serial admit path BEFORE the arena is wrapped in an `Arc` or any slice
        // is handed out — so no concurrent borrow can exist during the grow.
        // SAFETY: sole writer, no live borrow (the `Arc` is created only after this
        // returns).  `segment_size` is exactly the pool's segment size and the
        // capacity `reserve_full_capacity` just guaranteed, so `grow_uninit` performs
        // no realloc and the slot stays in one segment.  `u8` has no validity
        // invariant; every live byte is written exactly once — block slots by their
        // inflate worker, the `[0, FRONT_REGION)` carry region by `FindBoundariesAndSort`
        // — before any read.  Bytes never read (the unused tail past a run's end, and
        // the front region when there is no straddler) are never observed: `FBS` only
        // ever slices `[scan_start, run_end)`.
        #[allow(unsafe_code)]
        let _all = unsafe { arena.grow_uninit(self.segment_size) };
        self.arena = Some(Arc::new(PooledSegmentedBuf::pooled(arena, Arc::clone(&self.pool))));
        self.run_offset = FRONT_REGION as u64;
        true
    }

    /// Admit one [`BgzfBlock`]: ensure the arena exists, assign the block an
    /// arithmetic slot offset, and emit it eagerly (the PRIOR withheld block, now
    /// confirmed to have a same-run successor, is pushed to the output queue; THIS
    /// block becomes the new withheld `last_block`).
    ///
    /// Returns `Ok(())` on success, or `Err(b)` handing the block back if the pool
    /// is exhausted (the arena could not be acquired) — so the caller never loses
    /// the block.
    ///
    /// No `unsafe` here: the arena was grown to full capacity once by `ensure_arena`,
    /// so the slot `(offset, len)` is already a live region of the shared arena;
    /// the inflate worker writes it via `slice_mut`.
    pub(crate) fn admit_block(&mut self, b: BgzfBlock) -> Result<(), BgzfBlock> {
        if !self.ensure_arena() {
            return Err(b);
        }
        let arena = self.arena.as_ref().expect("arena present after ensure_arena");

        let len = b.uncompressed_size;
        let offset = self.run_offset;
        self.run_offset += u64::from(len);

        let ordinal = self.next_ordinal;
        self.next_ordinal += 1;
        self.run_cumulative += len as usize;

        // `is_last_of_run` / `seals_to_spill` are run-level facts known only at seal,
        // so the just-admitted block is withheld as the tentative last block.  The
        // PRIOR withheld block now has a same-run successor → it is NOT last → emit it.
        let block = ArenaBlock {
            arena: Arc::clone(arena),
            ordinal,
            offset,
            len,
            block: b.bytes,
            is_last_of_run: false,
            run_seq: self.run_seq,
            seals_to_spill: false,
        };
        if let Some(prev) = self.last_block.replace(block) {
            self.emit.push_back(prev);
        }
        Ok(())
    }

    /// Seal the current run: stamp the withheld `last_block` as the run's final
    /// block (`is_last_of_run = true` plus the run's `seals_to_spill`), emit it,
    /// drop this step's arena handle (so the arena can return to the pool once all
    /// emitted clones release it), and advance `run_seq`.
    ///
    /// All non-last blocks of the run were already emitted eagerly by `admit_block`;
    /// only the single withheld block remains.
    ///
    /// `seals_to_spill`: `true` for mid-stream seals (the run becomes a disk spill),
    /// `false` for the final residual run.
    ///
    /// A `seal_run` call always follows at least one admit for the run, so
    /// `last_block` is `Some`; the `try_run` EOF guard returns before sealing an
    /// arena-less, blockless state.
    fn seal_run(&mut self, seals_to_spill: bool) {
        if let Some(mut last) = self.last_block.take() {
            last.is_last_of_run = true;
            last.seals_to_spill = seals_to_spill;
            self.emit.push_back(last);
        }
        // Drop our handle to the sealed run's arena; the emitted blocks (and their
        // downstream consumers) keep it alive until they release their `Arc` clones,
        // at which point `PooledSegmentedBuf::drop` returns it to the capacity-1 pool.
        self.arena = None;
        self.run_seq += 1;
        self.run_cumulative = 0;
        self.run_offset = 0;
    }

    /// Test seam: admit blocks, then freeze the arena and drain all emitted
    /// [`ArenaBlock`]s into a `Vec` for direct inspection, bypassing the pipeline
    /// framework.  Uses `seals_to_spill = false` (residual) for the single flush.
    #[cfg(test)]
    pub(crate) fn seal_and_drain_for_test(&mut self) -> Vec<ArenaBlock> {
        self.seal_run(false);
        self.emit.drain(..).collect()
    }
}

impl Step for ReadBlocks {
    type Input = BgzfBlock;
    type Outputs = OrderedBytesSingle<ArenaBlock>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "ReadBlocks",
            // Off-pool on the coordination driver (N+2): the serial arena-admit
            // runs on a dedicated thread instead of stealing a pool worker slot
            // from the parallel inflaters. Detached collapses the `ByItemOrdinal`
            // output to `None` exactly as `Serial` did (transport-identical).
            kind: StepKind::Detached,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn detached_group(&self) -> DetachedGroup {
        DetachedGroup::Shared(crate::sort::SORT_COORD_GROUP)
    }

    // A single cohesive deferred-seal state machine: held-output retry, seal
    // arming/execution, batched admit, and staged-emit drain are tightly coupled
    // by `seal_pending` / `deferred_block` / `held` and are clearer read top to
    // bottom than split across helpers that would each need the same state.
    #[allow(clippy::too_many_lines)]
    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        // 1. Retry any held output first (backpressure path).
        if let Some(unpushed) = self.held.take() {
            match ctx.outputs.retry(unpushed) {
                Ok(()) => {}
                Err(again) => {
                    self.held.put(again);
                    return Ok(StepOutcome::Contention);
                }
            }
        }

        // 2. Drain any staged emitted ArenaBlocks before accepting new input.
        if let Some(block) = self.emit.pop_front() {
            match ctx.outputs.push(block) {
                Ok(()) => return Ok(StepOutcome::Progress),
                Err(unpushed) => {
                    self.held.put(unpushed);
                    return Ok(StepOutcome::Progress);
                }
            }
        }

        // 3. If finished and emit is empty, we are done.
        if self.finished {
            return Ok(StepOutcome::Finished);
        }

        // 4. Re-admit a block held over from a prior deferred seal, if any.  The
        //    seal already fired; we just need an arena to land this block in.
        if self.deferred_block.is_some() {
            // Drain any still-staged emit blocks first so the pool can recycle.
            if let Some(out) = self.emit.pop_front() {
                match ctx.outputs.push(out) {
                    Ok(()) => return Ok(StepOutcome::Progress),
                    Err(unpushed) => {
                        self.held.put(unpushed);
                        return Ok(StepOutcome::Progress);
                    }
                }
            }
            // Acquire the arena BEFORE consuming the held block, so a pool
            // exhaustion does not drop it.
            if !self.ensure_arena() {
                // Pool still exhausted (the just-sealed run's arena has not
                // returned yet): keep holding the block and backpressure.
                return Ok(StepOutcome::NoProgress);
            }
            let block = self.deferred_block.take().expect("deferred_block present");
            match self.admit_block(block) {
                Ok(()) => {}
                Err(b) => {
                    // ensure_arena just returned true, so admit cannot fail; restore
                    // the block defensively rather than drop it.
                    self.deferred_block = Some(b);
                    debug_assert!(false, "admit_block must succeed after ensure_arena");
                    return Ok(StepOutcome::NoProgress);
                }
            }
            if self.run_cumulative >= self.run_cap {
                self.seal_pending = true;
            }
            return Ok(StepOutcome::Progress);
        }

        // 5. Pop and admit BgzfBlocks — BATCHED (up to `ADMIT_BATCH`) to keep the
        //    parallel inflaters fed.  See `ADMIT_BATCH`: one-per-dispatch made this
        //    admit the input-side bottleneck under the coordination driver's
        //    drain-first restart.  Only the clean admit path is batched; a seal
        //    boundary and
        //    output backpressure both break the batch and return, preserving the
        //    exact single-dispatch semantics of the delicate deferred-seal state
        //    machine (a seal is rare — ~one per spill).
        let mut admitted_any = false;
        for _ in 0..ADMIT_BATCH {
            let Some(block) = ctx.input.pop() else { break };
            // A seal is armed from a previous block hitting the budget.  The
            // arrival of THIS block proves a following run exists, so it is safe
            // to seal the previous run as a disk spill (its straddler carry will
            // be completed by the run this block opens).  End the batch here: the
            // just-sealed run's arena must be consumed before the next admit.
            if self.seal_pending {
                self.seal_run(true); // deferred mid-stream seal → spill
                self.seal_pending = false;
                // The just-sealed run's arena is held (as an Arc) in `emit` until
                // it is pushed downstream and consumed, so the capacity-1 pool is
                // momentarily empty.  Hold this block and re-admit it once the
                // arena returns (step 4 on a later call).  Drain `emit` now.
                self.deferred_block = Some(block);
                if let Some(out) = self.emit.pop_front() {
                    match ctx.outputs.push(out) {
                        Ok(()) => return Ok(StepOutcome::Progress),
                        Err(unpushed) => {
                            self.held.put(unpushed);
                            return Ok(StepOutcome::Progress);
                        }
                    }
                }
                return Ok(StepOutcome::Progress);
            }
            match self.admit_block(block) {
                Ok(()) => {}
                Err(b) => {
                    // Pool exhausted on a fresh run's first block: hold the block and
                    // backpressure until the prior run's chunk is consumed and its
                    // arena returns.  (Mid-run admits never fail — the arena is
                    // already acquired — so this is the cross-run boundary case.)
                    // If we already admitted this batch, that IS progress.
                    self.deferred_block = Some(b);
                    return Ok(if admitted_any {
                        StepOutcome::Progress
                    } else {
                        StepOutcome::NoProgress
                    });
                }
            }
            admitted_any = true;
            // Check if we hit the run budget; arm the deferred seal.  The next
            // loop iteration observes `seal_pending` and runs the seal path above.
            if self.run_cumulative >= self.run_cap {
                self.seal_pending = true;
            }
            // Drain staged emit eagerly within the batch so the inflaters see
            // blocks promptly and `emit` stays bounded.  On output backpressure,
            // hold and return — the batch ends; the rest is picked up next
            // dispatch (step 1 retries the held item first).
            while let Some(out) = self.emit.pop_front() {
                match ctx.outputs.push(out) {
                    Ok(()) => {}
                    Err(unpushed) => {
                        self.held.put(unpushed);
                        return Ok(StepOutcome::Progress);
                    }
                }
            }
        }
        if admitted_any {
            return Ok(StepOutcome::Progress);
        }

        // 6. No input available.
        if !ctx.input.is_drained() {
            return Ok(StepOutcome::NoProgress);
        }

        // 7. Input fully drained: seal the final (residual) run.  A pending seal
        //    is DOWNGRADED to a residual here — there is no following run to
        //    complete a straddler carry, and the final run's last record always
        //    completes at EOF, so the carry is empty.  `seal_run(false)` emits
        //    the residual.
        self.seal_pending = false;
        // Only seal if a block is withheld; if no block is held and the arena is
        // None, there was no data at all after the previous seal — i.e. truly
        // empty input (no blocks ever admitted).  NOTE: the deferred-seal design
        // guarantees the *final* run is always a residual (step 5 only seals a
        // spill once a following block proves a next run exists; otherwise this
        // step downgrades the pending seal to a residual), so this branch is NOT
        // reached after any data — it is the empty-input case only.
        if self.arena.is_none() && self.last_block.is_none() {
            self.finished = true;
            return Ok(StepOutcome::Finished);
        }
        self.seal_run(false); // final run → residual
        self.finished = true;

        if let Some(block) = self.emit.pop_front() {
            match ctx.outputs.push(block) {
                Ok(()) => return Ok(StepOutcome::Progress),
                Err(unpushed) => {
                    self.held.put(unpushed);
                    return Ok(StepOutcome::Progress);
                }
            }
        }

        Ok(StepOutcome::Finished)
    }

    fn new_worker_copy(&self) -> Self {
        // Serial steps are never cloned by the framework; this is unreachable.
        panic!("ReadBlocks is Serial — new_worker_copy should never be called")
    }
}

// ============================================================================
// Item types
// ============================================================================

/// Input to `InflateToArena`.
///
/// Carries a grown-but-uninit slot `(offset, len)` in `arena` plus the full
/// raw BGZF `block` bytes to inflate into it.  The slot was reserved by
/// `grow_uninit` on the serial admit path; it must be fully written by the
/// inflate worker before any read.
pub struct ArenaBlock {
    /// The shared arena into which this block's bytes will be inflated.
    pub arena: Arc<PooledSegmentedBuf>,
    /// Global ordinal, used by `ByItemOrdinal` reordering so downstream steps
    /// receive blocks in the original file order.
    pub ordinal: u64,
    /// Byte offset of this block's slot within the arena (returned by
    /// `grow_uninit`).
    pub offset: u64,
    /// Uncompressed size (ISIZE from the BGZF footer == slot length).
    pub len: u32,
    /// Complete raw BGZF block bytes (header + deflate payload + footer).
    pub block: Vec<u8>,
    /// `true` if this is the last block of the current run (e.g. a BAM file
    /// segment); used by downstream steps to detect run boundaries.
    pub is_last_of_run: bool,
    /// Run sequence number; increments with each new run.
    pub run_seq: u32,
    /// `true` if this run seals to a disk spill (mid-stream seal); `false` if
    /// this is the final residual run (in-memory, not spilled).
    pub seals_to_spill: bool,
}

impl HeapSize for ArenaBlock {
    fn heap_size(&self) -> usize {
        self.block.len()
    }
}

impl Ordered for ArenaBlock {
    fn ordinal(&self) -> u64 {
        self.ordinal
    }
}

/// Completion token emitted by `InflateToArena` after successfully inflating
/// one block.  The decompressed bytes now live in `arena` at the byte range
/// `offset..offset + len`.  The token carries no heap data of its own
/// (`heap_size == 0`).
pub struct InflatedBlock {
    /// The arena holding the decompressed bytes.
    pub arena: Arc<PooledSegmentedBuf>,
    /// Global ordinal (same value as the originating `ArenaBlock`).
    pub ordinal: u64,
    /// Byte offset of the decompressed data within the arena.
    pub offset: u64,
    /// Byte length of the decompressed data.
    pub len: u32,
    /// Forwarded from `ArenaBlock`.
    pub is_last_of_run: bool,
    /// Forwarded from `ArenaBlock`.
    pub run_seq: u32,
    /// Forwarded from `ArenaBlock`: `true` if this run seals to a disk spill.
    pub seals_to_spill: bool,
}

impl HeapSize for InflatedBlock {
    fn heap_size(&self) -> usize {
        0
    }
}

impl Ordered for InflatedBlock {
    fn ordinal(&self) -> u64 {
        self.ordinal
    }
}

// ============================================================================
// Step
// ============================================================================

/// `Parallel + ByItemOrdinal` step that decompresses each `ArenaBlock`'s BGZF
/// bytes directly into its pre-reserved arena slot.
///
/// Each parallel worker clone holds its own [`Decompressor`] (allocated via
/// `new_worker_copy`) so there is no contention on the decompression state.
pub struct InflateToArena {
    decompressor: Decompressor,
    held: HeldSlot<Unpushed<InflatedBlock>>,
    output_byte_limit: u64,
}

impl InflateToArena {
    /// Create a new `InflateToArena` step.
    ///
    /// `output_byte_limit` bounds the byte-counted output queue (tokens carry
    /// `heap_size == 0`, so this mainly controls queue depth via the framework's
    /// backpressure mechanism).
    #[must_use]
    pub fn new(output_byte_limit: u64) -> Self {
        Self { decompressor: Decompressor::new(), held: HeldSlot::new(), output_byte_limit }
    }

    /// Decompress `item.block` into the arena slot `(item.offset, item.len)`,
    /// returning an [`InflatedBlock`] completion token on success.
    ///
    /// # Errors
    ///
    /// Returns an `io::Error` if BGZF decompression fails or if the
    /// decompressed length does not match `item.len`.
    ///
    /// # Safety (caller contract)
    ///
    /// See the `#[allow(unsafe_code)]` block inside — the caller must have
    /// reserved the slot via `grow_uninit` before constructing the `ArenaBlock`.
    // `pub(crate)` so the unit test below can call it directly without wiring up
    // the full pipeline framework.
    pub(crate) fn inflate_one(&mut self, item: ArenaBlock) -> io::Result<InflatedBlock> {
        let ArenaBlock {
            arena,
            ordinal,
            offset,
            len,
            block,
            is_last_of_run,
            run_seq,
            seals_to_spill,
        } = item;

        // SAFETY: `(offset, len)` was reserved by `grow_uninit` on the serial
        // ReadBlocks admit path before this `ArenaBlock` was enqueued, so the
        // slot is live within the arena's allocated storage.  The ISIZE
        // prefix-sum partitions the arena into non-overlapping slots, so this
        // `&mut [u8]` aliases no other concurrent inflate worker's slice.
        // `u8` has no validity invariant, so writing before reading is the only
        // required contract — `decompress_into_slice` below fills every byte.
        // `offset` is a u64 byte offset into the arena; on a 64-bit platform
        // this always fits in usize — the arena itself cannot exceed
        // `isize::MAX` bytes, which is the Rust allocation bound.
        let offset_usize =
            usize::try_from(offset).expect("arena offset must fit in usize on this platform");
        let len_usize = len as usize; // u32 always fits in usize

        #[allow(unsafe_code)]
        let slot = unsafe { arena.slice_mut(offset_usize, len_usize) };

        let n = decompress_into_slice(&block, &mut self.decompressor, slot)?;
        debug_assert_eq!(n, len_usize, "decompressed length must match ISIZE");

        Ok(InflatedBlock { arena, ordinal, offset, len, is_last_of_run, run_seq, seals_to_spill })
    }
}

impl Step for InflateToArena {
    type Input = ArenaBlock;
    type Outputs = OrderedBytesSingle<InflatedBlock>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "InflateToArena",
            kind: StepKind::Parallel,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        // Retry any held output first (backpressure path — mirror BgzfDecompress).
        if let Some(unpushed) = self.held.take() {
            match ctx.outputs.retry(unpushed) {
                Ok(()) => {}
                Err(again) => {
                    self.held.put(again);
                    // `Contention` keeps the worker alive for retry — `NoProgress`
                    // would let the framework silently drop the held item if input
                    // is also drained.
                    return Ok(StepOutcome::Contention);
                }
            }
        }

        let Some(block) = ctx.input.pop() else {
            // No input this call.  If upstream is fully drained the held slot was
            // already flushed by the Contention preamble above, so every item has
            // been processed.  For a Parallel step only the last clone to finish
            // closes the shared output (gated by StepDrainCounter in the driver).
            if ctx.input.is_drained() {
                return Ok(StepOutcome::Finished);
            }
            return Ok(StepOutcome::NoProgress);
        };

        let inflated = self.inflate_one(block)?;

        match ctx.outputs.push(inflated) {
            Ok(()) => Ok(StepOutcome::Progress),
            Err(unpushed) => {
                self.held.put(unpushed);
                Ok(StepOutcome::Progress)
            }
        }
    }

    fn new_worker_copy(&self) -> Self {
        Self::new(self.output_byte_limit)
    }
}

// ============================================================================
// Arena sort strategy (per-order key extraction + seal)
// ============================================================================

/// Per-sort-order policy for the arena-front [`FindBoundariesAndSort`] scan.
///
/// The strategy owns the accumulated per-record refs and knows how to (a) extract
/// a record's sort key from its arena-resident body and (b) at run seal, sort
/// those refs and wrap the shared arena into the correctly-typed
/// [`MemoryChunkErased`] — with zero record-body copies. `FindBoundariesAndSort`
/// is generic over this trait and monomorphised per order, so the per-record
/// [`push_record`](Self::push_record) call inlines into the boundary scan hot
/// loop with no dynamic dispatch.
pub trait ArenaSortStrategy: Send + 'static {
    /// Reserve capacity for approximately `est_records` refs at the start of a
    /// run, so the incremental per-record pushes do not reallocate a multi-GB ref
    /// buffer mid-run.
    fn reserve_for_run(&mut self, est_records: usize);

    /// Extract the sort key from `body` (the record's BAM body, `block_size`
    /// prefix excluded) and accumulate a ref pointing at `(body_off, len)` in the
    /// shared inflate arena. Called once per record on the boundary scan hot path.
    ///
    /// # Errors
    ///
    /// Returns an error if the record is invalid for this order (e.g. a
    /// template-coordinate dropped-lane violation). Coordinate never errors.
    fn push_record(&mut self, body: &[u8], body_off: u64, len: u32) -> io::Result<()>;

    /// Sort the refs accumulated for this run and wrap `arena` into the erased
    /// chunk (zero record copies), resetting the accumulator for the next run.
    fn seal(&mut self, arena: Arc<PooledSegmentedBuf>, sort_threads: usize) -> MemoryChunkErased;

    /// A fresh, empty strategy carrying the same configuration — used to build a
    /// worker copy of the step. `FindBoundariesAndSort` is `Serial`, so this only
    /// ever constructs its single working instance.
    #[must_use]
    fn fresh(&self) -> Self
    where
        Self: Sized;
}

/// Coordinate-order strategy: extracts the fixed `u64` coordinate key inline and
/// accumulates plain [`RecordRef`]s; seals to [`MemoryChunkErased::Coordinate`]
/// via [`coordinate_chunk_from_refs`].
pub struct CoordinateStrategy {
    /// BAM header reference-sequence count, used by [`extract_coordinate_key_inline`].
    n_ref: u32,
    /// Coordinate-key refs accumulated across the current run's blocks. Filled
    /// incrementally by the scan and consumed (via `mem::take`) at seal, leaving
    /// an empty `Vec` for the next run.
    refs: Vec<RecordRef>,
}

impl CoordinateStrategy {
    /// Create a coordinate strategy for a header with `n_ref` reference sequences.
    #[must_use]
    pub fn new(n_ref: u32) -> Self {
        Self { n_ref, refs: Vec::new() }
    }
}

impl ArenaSortStrategy for CoordinateStrategy {
    #[inline]
    fn reserve_for_run(&mut self, est_records: usize) {
        self.refs.reserve(est_records);
    }

    #[inline]
    fn push_record(&mut self, body: &[u8], body_off: u64, len: u32) -> io::Result<()> {
        let sort_key = extract_coordinate_key_inline(body, self.n_ref);
        self.refs.push(RecordRef::new(sort_key, body_off, len));
        Ok(())
    }

    fn seal(&mut self, arena: Arc<PooledSegmentedBuf>, sort_threads: usize) -> MemoryChunkErased {
        let refs = std::mem::take(&mut self.refs);
        MemoryChunkErased::Coordinate(coordinate_chunk_from_refs(arena, refs, sort_threads))
    }

    fn fresh(&self) -> Self {
        Self::new(self.n_ref)
    }
}

/// Template-coordinate strategy: wraps a [`fgumi_sort::TemplateArenaAccumulator`],
/// which owns the library / cell-barcode / MI and `--key-types` narrowed-lane
/// machinery the template key needs; accumulates arena-pointing refs and seals to
/// [`MemoryChunkErased::TemplateCoordinate`] (an arena-backed
/// `InMemoryChunk<TemplateKey>`), byte-identical to the owned `TemplateChunkSorter`.
pub struct TemplateStrategy {
    acc: fgumi_sort::TemplateArenaAccumulator,
}

impl TemplateStrategy {
    /// Wrap a template accumulator (built from the header via
    /// [`fgumi_sort::TemplateArenaAccumulator::from_header`]).
    #[must_use]
    pub fn new(acc: fgumi_sort::TemplateArenaAccumulator) -> Self {
        Self { acc }
    }
}

impl ArenaSortStrategy for TemplateStrategy {
    #[inline]
    fn reserve_for_run(&mut self, est_records: usize) {
        self.acc.reserve(est_records);
    }

    #[inline]
    fn push_record(&mut self, body: &[u8], body_off: u64, len: u32) -> io::Result<()> {
        self.acc
            .push(body, body_off, len)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("{e:#}")))
    }

    fn seal(&mut self, arena: Arc<PooledSegmentedBuf>, sort_threads: usize) -> MemoryChunkErased {
        MemoryChunkErased::TemplateCoordinate(self.acc.seal(arena, sort_threads))
    }

    fn fresh(&self) -> Self {
        Self { acc: self.acc.fresh() }
    }
}

/// Queryname-order strategy (lexicographic or natural, selected by the key type
/// `K` + `wrap`): extracts the embedded read-name key from each record and
/// accumulates `(key, offset, len)` refs pointing into the shared arena; at seal
/// it comparator-sorts the refs (variable-length names are not radix-able) into
/// an arena-backed [`InMemoryChunk<K>`], wrapped into the matching erased variant.
///
/// The record bodies stay in the arena (zero-copy); only the small name bytes are
/// owned by each key, exactly as the legacy queryname sort. Queryname's
/// tie order is unspecified (the integration parity gate is name-order, not
/// byte-identity), so one globally-sorted chunk per run is correct.
pub struct QuerynameStrategy<K: RawSortKey> {
    /// `(key, body_offset, len)` refs accumulated across the current run's blocks.
    refs: Vec<(K, u64, u32)>,
    /// Bounded rayon pool (sized to `sort_threads`) the per-run comparator sort
    /// installs into, so the parallel sort does not oversubscribe the pipeline's
    /// worker pool on a spill. Built once on the first seal, reused; `None` on a
    /// fresh copy. Mirrors [`TemplateArenaAccumulator`](fgumi_sort::TemplateArenaAccumulator)'s pool.
    sort_pool: Option<rayon::ThreadPool>,
    /// Erases the sorted `InMemoryChunk<K>` into the correct `MemoryChunkErased`
    /// arm (`QuerynameLex` for the lex key, `QuerynameNatural` for the natural
    /// key). A `fn` pointer so [`fresh`](ArenaSortStrategy::fresh) can copy it.
    wrap: fn(InMemoryChunk<K>) -> MemoryChunkErased,
}

impl<K: RawSortKey> QuerynameStrategy<K> {
    /// Build a queryname strategy that erases its sealed chunk via `wrap`.
    #[must_use]
    pub fn new(wrap: fn(InMemoryChunk<K>) -> MemoryChunkErased) -> Self {
        Self { refs: Vec::new(), sort_pool: None, wrap }
    }

    /// Test-only observation of the bounded sort pool's thread count, so a test
    /// can assert the Phase-1 `sort_threads` value actually SIZED the worker pool
    /// (the runtime effect), not merely that it was plumbed. `None` until the
    /// first [`seal`](ArenaSortStrategy::seal) builds the pool.
    #[cfg(test)]
    pub(crate) fn sort_pool_threads(&self) -> Option<usize> {
        self.sort_pool.as_ref().map(rayon::ThreadPool::current_num_threads)
    }
}

impl<K: RawSortKey + 'static> ArenaSortStrategy for QuerynameStrategy<K> {
    #[inline]
    fn reserve_for_run(&mut self, est_records: usize) {
        self.refs.reserve(est_records);
    }

    #[inline]
    fn push_record(&mut self, body: &[u8], body_off: u64, len: u32) -> io::Result<()> {
        // Queryname keys are EMBEDDED_IN_RECORD: the name lives in `body`, so the
        // key is reconstructed straight from it (no SortContext needed).
        self.refs.push((K::extract_from_record(body), body_off, len));
        Ok(())
    }

    fn seal(&mut self, arena: Arc<PooledSegmentedBuf>, sort_threads: usize) -> MemoryChunkErased {
        let refs = std::mem::take(&mut self.refs);
        let pool = self.sort_pool.get_or_insert_with(|| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(sort_threads.max(1))
                .thread_name(|i| format!("qname-sort-{i}"))
                .build()
                .expect("build bounded queryname-sort rayon pool")
        });
        let wrap = self.wrap;
        pool.install(move || wrap(queryname_chunk_from_arena_refs(arena, refs)))
    }

    fn fresh(&self) -> Self {
        Self::new(self.wrap)
    }
}

// ============================================================================
// FindBoundariesAndSort step
// ============================================================================

/// `Serial` step that scans each run's contiguous arena span directly (no copy),
/// skips the BAM header on run 0 via [`bam_header_len`], builds `(body_offset,
/// block_size)` refs, sorts via [`coordinate_chunk_from_arena_refs`], and emits
/// `SortChunkEvent::Spill` (for sealed mid-stream runs) or `SortChunkEvent::Residual`
/// (for the final run), followed by `SortChunkEvent::AllAnnounced`.
///
/// **Straddler handling (3b.4 front-region carry):** A BAM record that spans the
/// boundary between run k and run k+1 is handled as follows:
/// - While scanning run k's span, the trailing partial record (a 4-byte
///   `block_size` prefix whose body extends past the run's end) is detected when
///   `cur + 4 <= span.len()` but `cur + 4 + bs > span.len()`.  Those partial bytes
///   (`span[cur..]`) are copied into `self.carry`.
/// - When the first block of run k+1 arrives (arena already frozen), the carry is
///   written right-aligned into that arena's front region at
///   `[FRONT_REGION - carry.len(), FRONT_REGION)` via `unsafe { arena.slice_mut(...) }`.
///   The scan then starts at `FRONT_REGION - carry.len()` so that `carry ++ block-0-head`
///   is physically contiguous, forming one complete record ref that is naturally
///   record 0 after sorting.
///
/// **Parity invariant:** sealed mid-stream runs emit `Spill{seq: 0, 1, ...}`;
/// exactly one final run emits `Residual`.  `AllAnnounced` carries
/// `slot_count = spilled_run_count, memory_chunk_count = 1`.
pub struct FindBoundariesAndSort<S: ArenaSortStrategy> {
    /// Per-order key extraction + seal policy; owns the accumulated sort refs.
    strategy: S,
    /// Threads handed to the per-chunk coordinate sort (`coordinate_chunk_from_refs`);
    /// `>1` enables the parallel radix on large chunks.  Matches the pipeline's
    /// configured thread count.
    sort_threads: usize,
    output_byte_limit: u64,
    /// Held slot for backpressure retry on emitting events.
    held: HeldSlot<Unpushed<SortChunkEvent>>,
    /// The shared arena that holds the current run's decompressed bytes (set on
    /// the first block of each run, cleared after the run is scanned).
    arena: Option<Arc<PooledSegmentedBuf>>,
    /// Arena offset of the scan start for the current run.
    ///
    /// - Run 0 with empty carry: `FRONT_REGION + bam_header_len(...)`.
    /// - Run k>0 with non-empty carry: `FRONT_REGION - carry.len()` (after
    ///   writing carry into the front region).
    /// - Run k>0 with empty carry (no straddler): `FRONT_REGION`.
    scan_start: u64,
    /// Arena offset one past the last byte of the current run span (updated with
    /// each block).
    run_end: u64,
    /// Arena offset of the next un-parsed record's first byte within the current
    /// run.  The incremental scan advances this as in-order blocks extend `run_end`,
    /// parking it at the start of the first record that overruns the bytes inflated
    /// so far; it resumes from here when the next block arrives.  At seal,
    /// `[scan_cursor, run_end)` is exactly the trailing partial record (the straddler
    /// carry for a spill, or empty for a clean residual).
    scan_cursor: u64,
    /// `true` once run 0's BAM header has been skipped.  Set `true` immediately for
    /// runs `> 0` (no header).  Stays `false` on run 0 until enough blocks have been
    /// inflated for [`bam_header_len`] to parse the full header — a many-reference
    /// header can exceed a single 64 KiB block, so the skip may take several blocks.
    header_skipped: bool,
    /// Pending events staged after each run scan; drained in subsequent `try_run`
    /// calls.  Holds `Spill`/`Residual` then `AllAnnounced` (only on the final run).
    pending: VecDeque<SortChunkEvent>,
    /// `true` once the final run has been scanned and `AllAnnounced` staged.
    finalized: bool,
    /// Trailing partial record bytes from the previous run (the carry buffer).
    /// Non-empty when run k ended mid-record; written into run k+1's front region.
    carry: Vec<u8>,
    /// Sequence number for the next run's `Spill` event (monotonically increasing).
    next_seq: u32,
    /// Cumulative record count across all runs.
    total_records: u64,
    /// Number of runs that sealed to a disk spill (all runs except the final one).
    spilled_run_count: u32,
}

impl<S: ArenaSortStrategy> FindBoundariesAndSort<S> {
    /// Test-only borrow of the per-order strategy, so a test can observe
    /// strategy-owned state (e.g. the bounded sort pool built at seal) after
    /// driving the step — confirming the `sort_threads` handed to [`new`](Self::new)
    /// is forwarded into `strategy.seal`.
    #[cfg(test)]
    pub(crate) fn strategy(&self) -> &S {
        &self.strategy
    }

    /// Create a new `FindBoundariesAndSort` step over the given per-order
    /// `strategy` (e.g. [`CoordinateStrategy`], which carries the header's
    /// reference-sequence count for key extraction). `sort_threads` is the thread
    /// count handed to the per-chunk sort (the pipeline's configured threads).
    /// `output_byte_limit` byte-bounds the output queue.
    #[must_use]
    pub fn new(strategy: S, sort_threads: usize, output_byte_limit: u64) -> Self {
        Self {
            strategy,
            sort_threads,
            output_byte_limit,
            held: HeldSlot::new(),
            arena: None,
            scan_start: 0,
            run_end: 0,
            scan_cursor: 0,
            header_skipped: false,
            pending: VecDeque::new(),
            finalized: false,
            carry: Vec::new(),
            next_seq: 0,
            total_records: 0,
            spilled_run_count: 0,
        }
    }

    /// Ingest one `InflatedBlock`, extending the current run's contiguous arena span
    /// and then parsing every record made complete by it via
    /// [`scan_available`](Self::scan_available) — so the scan overlaps the run's
    /// still-inflating tail (Inflate‖Scan).  When the block is the last of a run
    /// (`is_last_of_run`), the run is sealed: the trailing partial record (if any)
    /// becomes the carry, the accumulated refs are sorted, and a `Spill` or
    /// `Residual` event is staged in `self.pending`.  After the final run's
    /// `Residual`, an `AllAnnounced` is also staged.
    ///
    /// On the first block of each run, the straddler carry from the previous run (if
    /// any) is written right-aligned into the new arena's front region via an
    /// `unsafe` `slice_mut` call (3rd `fgumi-pipeline-io` unsafe site — see
    /// CLAUDE.md §"Approved hot-path unsafe (parallel-inflate sort ingest,
    /// fgumi-pipeline-io)" for the full justification).
    ///
    /// # Errors
    ///
    /// Returns an `io::Error` if:
    /// - The carry to write into the front region is longer than `FRONT_REGION`.
    /// - A `block_size` value in the record stream overflows `u32`.
    /// - At seal: run 0's header never fully arrived, the carry exceeds `FRONT_REGION`,
    ///   or the final residual run ends mid-record (truncated BAM).
    pub(crate) fn ingest_block(&mut self, block: &InflatedBlock) -> io::Result<()> {
        let block_end = block.offset + u64::from(block.len);

        if self.arena.is_none() {
            // ----------------------------------------------------------------
            // First block of this run.
            // ----------------------------------------------------------------

            // Straddler carry: write carry right-aligned into the front region.
            let scan_start = if !self.carry.is_empty() {
                let l = self.carry.len();
                if l > FRONT_REGION {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "FindBoundariesAndSort: carry ({l} bytes) exceeds FRONT_REGION \
                             ({FRONT_REGION}); record exceeds max straddler size — \
                             raise --max-memory"
                        ),
                    ));
                }
                let write_offset = FRONT_REGION - l;
                // SAFETY: `[FRONT_REGION - l, FRONT_REGION)` is the front region reserved
                // by `ReadBlocks::ensure_arena` via `grow_uninit(FRONT_REGION)` before any
                // block slot was allocated.  Soundness rests on DISJOINTNESS, not on
                // ordering: this range lies strictly below `FRONT_REGION`, while every
                // inflate block slot lies at or above `FRONT_REGION` (block 0 starts there
                // by construction in ReadBlocks), so the `&mut [u8]` synthesized here
                // aliases no slot an `InflateToArena` worker is concurrently writing — and
                // run k+1's later blocks MAY still be inflating when this carry write runs
                // (FBS is serial but does not wait for the whole run to inflate).  The
                // arena's backing segment was frozen with `reserve_full_capacity` before
                // being shared, so no `grow_uninit`/realloc can move these bytes underneath
                // a live worker slice.  This write happens exactly once per straddler per
                // run; `u8` has no validity invariant.
                #[allow(unsafe_code)]
                let dst = unsafe { block.arena.slice_mut(write_offset, l) };
                dst.copy_from_slice(&self.carry);
                self.carry.clear();
                #[allow(clippy::cast_possible_truncation)]
                let start = write_offset as u64;
                start
            } else if self.next_seq == 0 {
                // Run 0: BAM header starts at block.offset (the first byte of the
                // run's data).  In the real pipeline block.offset == FRONT_REGION
                // (since ReadBlocks reserves the front region first), but in unit
                // tests the arena may be laid out without the front region.
                block.offset
            } else {
                // Non-first run, no carry: records start at FRONT_REGION.
                #[allow(clippy::cast_possible_truncation)]
                let start = FRONT_REGION as u64;
                start
            };

            self.arena = Some(Arc::clone(&block.arena));
            self.scan_start = scan_start;
            self.scan_cursor = scan_start;
            // Runs > 0 have no BAM header (it lives only at the very start of input);
            // their first record is the straddler/`FRONT_REGION` data, so nothing to
            // skip.  Run 0 must skip the header, which `scan_available` does once
            // enough blocks have arrived.
            self.header_skipped = self.next_seq != 0;
            // Pre-size the ref buffer to the run's upper bound (the arena holds at
            // most `segment_size` bytes of record data) so the incremental pushes do
            // not trigger doubling reallocations of a multi-GB `Vec` mid-run.
            let est = (block.arena.len() / 192).max(1024);
            self.strategy.reserve_for_run(est);
        } else {
            debug_assert_eq!(
                block.offset, self.run_end,
                "FindBoundariesAndSort: non-contiguous arena block (expected offset {}, got {})",
                self.run_end, block.offset
            );
        }
        self.run_end = block_end;

        // Parse every record made complete by the bytes inflated so far, overlapping
        // the scan with the still-inflating tail of the run (Inflate‖Scan).
        self.scan_available(&block.arena)?;

        if block.is_last_of_run {
            self.seal_run(block.seals_to_spill, block.run_seq)?;
        }

        Ok(())
    }

    /// Parse every record that is fully present in `[scan_cursor, run_end)`, pushing
    /// a [`RecordRef`] (with its coordinate key extracted inline) for each, and park
    /// `scan_cursor` at the first record that overruns the bytes inflated so far.
    ///
    /// Records are gap-free across block slots (the arena is one contiguous segment
    /// with prefix-summed offsets), so the cursor walks straight through block
    /// boundaries — a record that started in block k and continues into block k+1 is
    /// simply parsed once k+1 has extended `run_end`.  A trailing partial record at
    /// the *run* boundary is not handled here; it is left in `[scan_cursor, run_end)`
    /// for `seal_run` to carry (spill) or reject (truncated residual).
    ///
    /// For run 0, the BAM header is skipped first; if the header is not yet fully
    /// inflated (`bam_header_len` returns `None`), the scan returns and retries on the
    /// next block.
    ///
    /// # Errors
    ///
    /// Returns an `io::Error` on a `block_size` value that overflows `u32`.
    fn scan_available(&mut self, arena: &PooledSegmentedBuf) -> io::Result<()> {
        let scan_start = self.scan_start;
        let scan_start_usize = usize::try_from(scan_start).expect("scan_start must fit in usize");
        let run_end_usize = usize::try_from(self.run_end).expect("run_end must fit in usize");
        let avail_len = run_end_usize.checked_sub(scan_start_usize).expect("run_end >= scan_start");

        // Borrow the bytes inflated so far for this run directly from the arena — NO
        // copy.  The whole run lives in one segment, so this single slice spans
        // `[scan_start, run_end)` and is re-taken (cheaply) as `run_end` grows.
        let span = arena.slice(scan_start_usize, avail_len);

        // Run 0: skip the BAM header once it is fully present.  A many-reference
        // header can exceed one BGZF block, so `None` means "wait for more blocks",
        // NOT an error (the error is raised at seal if the header never completes).
        if !self.header_skipped {
            // `?` surfaces a wrong-magic stream as an error here; `Ok(None)` still
            // means "header not fully inflated yet — wait for more blocks".
            match bam_header_len(span)? {
                Some(h) => {
                    self.scan_cursor = scan_start + h as u64;
                    self.header_skipped = true;
                }
                None => return Ok(()),
            }
        }

        // Parse complete `[block_size(4)][body]` frames from the cursor, extracting
        // the coordinate key inline (the body is already resident in `span`).
        let mut cur = usize::try_from(self.scan_cursor - scan_start)
            .expect("cursor offset within run fits in usize");
        let mut run_records: u64 = 0;
        loop {
            if cur + 4 > span.len() {
                // Not even a full 4-byte length prefix present yet — wait.
                break;
            }
            let bs = u32::from_le_bytes([span[cur], span[cur + 1], span[cur + 2], span[cur + 3]])
                as usize;
            if cur + 4 + bs > span.len() {
                // Record body not fully inflated yet — wait for the next block.
                break;
            }
            // Software-prefetch a few records ahead to hide the cold-arena miss
            // latency of this forward scan.
            let pf = cur + SCAN_PREFETCH_DISTANCE;
            if pf < span.len() {
                prefetch_read_l1(&span[pf]);
            }
            #[allow(clippy::cast_possible_truncation)]
            let body_arena_off = scan_start + cur as u64 + 4;
            let bs_u32 = u32::try_from(bs).map_err(|_| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("FindBoundariesAndSort: block_size {bs} overflows u32"),
                )
            })?;
            self.strategy.push_record(&span[cur + 4..cur + 4 + bs], body_arena_off, bs_u32)?;
            cur += 4 + bs;
            run_records += 1;
        }
        self.scan_cursor = scan_start + cur as u64;
        self.total_records += run_records;
        Ok(())
    }

    /// Seal the current run: finalize the trailing carry, sort the refs accumulated
    /// incrementally by [`scan_available`](Self::scan_available), and stage a `Spill`
    /// or `Residual` event.
    ///
    /// By the time this is called, `scan_available` has already run for the run's last
    /// block (in `ingest_block`), so every complete record is in `self.refs` and
    /// `scan_cursor` is parked at the start of the trailing partial record (if any).
    /// `[scan_cursor, run_end)` is therefore exactly that partial record:
    /// - for a mid-stream (spill) seal it becomes the straddler carry for the next run;
    /// - for the final (residual) seal it MUST be empty — a non-empty tail means the
    ///   BAM ended mid-record (truncated/malformed), surfaced as a hard error.
    ///
    /// # Errors
    ///
    /// Returns an `io::Error` if run 0's header never fully arrived, on a truncated
    /// final record, or if the carry exceeds `FRONT_REGION`.
    fn seal_run(&mut self, seals_to_spill: bool, run_seq: u32) -> io::Result<()> {
        let arena = self.arena.take().expect("seal_run called with no arena");

        // Run 0's header must have been skipped by now; if `scan_available` never
        // managed it, the input is too short to contain a valid BAM header.
        if !self.header_skipped {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "FindBoundariesAndSort: input too short to contain a valid BAM header",
            ));
        }

        // The trailing partial record, if any, is `[scan_cursor, run_end)`.
        let scan_cursor_usize =
            usize::try_from(self.scan_cursor).expect("scan_cursor must fit in usize");
        let run_end_usize = usize::try_from(self.run_end).expect("run_end must fit in usize");
        let tail_len =
            run_end_usize.checked_sub(scan_cursor_usize).expect("run_end >= scan_cursor");

        if seals_to_spill {
            // Mid-stream seal → the tail is the straddler carry for the next run.
            if tail_len > 0 {
                let tail = arena.slice(scan_cursor_usize, tail_len);
                if tail_len > FRONT_REGION {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "FindBoundariesAndSort: carry ({tail_len} bytes) exceeds FRONT_REGION \
                             ({FRONT_REGION}); record exceeds max straddler size — \
                             raise --max-memory"
                        ),
                    ));
                }
                self.carry.extend_from_slice(tail);
            }
        } else if tail_len > 0 {
            // Final (residual) run with a leftover tail → truncated/malformed BAM.
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "truncated BAM: final record incomplete at end of input",
            ));
        }

        // Sort the refs accumulated across this run's blocks and wrap the arena into
        // the erased chunk — no record copies (keys were extracted during the scan).
        let chunk = self.strategy.seal(Arc::clone(&arena), self.sort_threads);

        // Verify ReadBlocks and FindBoundariesAndSort are in lockstep on run
        // sequencing.  A desync here indicates a bug in the pipeline wiring (e.g.
        // ReadBlocks emitting the wrong seq on an is_last_of_run block).
        debug_assert_eq!(
            run_seq, self.next_seq,
            "ReadBlocks run_seq desynced from FindBoundariesAndSort seq counter"
        );

        let seq = self.next_seq;
        self.next_seq += 1;

        if seals_to_spill {
            // Mid-stream sealed run → disk spill.
            debug_assert!(
                !self.finalized,
                "FindBoundariesAndSort: Spill emitted after finalization"
            );
            self.pending.push_back(SortChunkEvent::Spill {
                seq,
                chunk,
                records_ingested_so_far: self.total_records,
            });
            self.spilled_run_count += 1;
        } else {
            // Final (residual) run.  A non-empty carry here means the BAM ends
            // mid-record — the file is truncated or malformed.  We must surface
            // this as a hard error in all build profiles: a debug_assert! would
            // silently drop the partial record in release builds.
            if !self.carry.is_empty() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "truncated BAM: final record incomplete at end of input",
                ));
            }
            self.pending.push_back(SortChunkEvent::Residual {
                chunk,
                records_ingested_so_far: self.total_records,
            });
            self.pending.push_back(SortChunkEvent::AllAnnounced {
                slot_count: self.spilled_run_count,
                memory_chunk_count: 1,
                total_records: self.total_records,
            });
            self.finalized = true;
        }

        // Reset per-run span tracking.  `self.refs` was emptied by `mem::take`; the
        // next run's first block re-initializes `scan_cursor`/`header_skipped`.
        self.scan_start = 0;
        self.run_end = 0;
        self.scan_cursor = 0;
        self.header_skipped = false;

        Ok(())
    }

    /// Called when input is fully drained and no `is_last_of_run` block was received
    /// (empty pipeline — no arena was ingested at all).  Stages a well-formed
    /// sentinel so downstream steps can complete.
    ///
    /// The "all runs sealed to spills, no residual" branch below is DEFENSIVE: the
    /// `ReadBlocks` deferred-seal design guarantees the final run is always a
    /// residual (it seals a spill only once a following block proves a next run
    /// exists, and downgrades a pending seal to a residual at EOF), so for valid
    /// input that branch is unreachable.  It is retained as a hard error on a
    /// non-empty carry so a genuinely truncated BAM (or a future regression in the
    /// seal logic) surfaces loudly instead of silently dropping a record.
    ///
    /// Returns `Some(first_event)` popped from `self.pending`, or `None` if no
    /// arena was ingested.
    pub(crate) fn finalize(&mut self) -> io::Result<Option<SortChunkEvent>> {
        if self.finalized {
            // Already handled by the last is_last_of_run block — nothing to do.
            return Ok(self.pending.pop_front());
        }
        if self.arena.is_none() && self.next_seq == 0 {
            // No arena was ever ingested (empty input).
            return Ok(None);
        }
        if self.arena.is_none() && self.next_seq > 0 {
            // DEFENSIVE / unreachable for valid input: the deferred-seal design
            // always makes the final run a residual, so we should never finalize
            // with spills-only-and-no-residual. If we somehow do, a non-empty carry
            // means the BAM ends mid-record (truncated input or a seal-logic
            // regression) — surface it as a hard error rather than drop a record.
            if !self.carry.is_empty() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "truncated BAM: final record incomplete at end of input (no-residual path)",
                ));
            }
            // Emit AllAnnounced with memory_chunk_count=0 so SortMerge completes.
            self.pending.push_back(SortChunkEvent::AllAnnounced {
                slot_count: self.spilled_run_count,
                memory_chunk_count: 0,
                total_records: self.total_records,
            });
            self.finalized = true;
            return Ok(self.pending.pop_front());
        }
        // Edge case: arena present but no is_last_of_run seen — treat as residual.
        if self.arena.is_some() {
            self.seal_run(false, self.next_seq)?;
        }
        Ok(self.pending.pop_front())
    }

    fn flush_held(&mut self, ctx: &mut StepCtx<'_, Self>) -> bool {
        // `true` once the slot is clear (was empty, or the held event flushed);
        // `false` while it's still held under backpressure. Uses the canonical
        // re-hold helper so the put-back-on-reject invariant lives in one place.
        !matches!(ctx.outputs.retry_held(&mut self.held), HeldRetry::StillHeld)
    }

    fn emit_pending(&mut self, ctx: &mut StepCtx<'_, Self>) -> StepOutcome {
        let Some(event) = self.pending.pop_front() else {
            return StepOutcome::NoProgress;
        };
        match ctx.outputs.push(event) {
            Ok(()) => StepOutcome::Progress,
            Err(unpushed) => {
                self.held.put(unpushed);
                StepOutcome::Progress
            }
        }
    }
}

impl<S: ArenaSortStrategy> Step for FindBoundariesAndSort<S> {
    type Input = InflatedBlock;
    type Outputs = Single<SortChunkEvent>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "FindBoundariesAndSort",
            // Off-pool on the coordination driver (N+2): the serial boundary scan
            // + radix sort + seal/spill framing runs on the dedicated coordination
            // thread, keeping the pool on pure parallel inflate/compress.
            kind: StepKind::Detached,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::None],
        }
    }

    fn detached_group(&self) -> DetachedGroup {
        DetachedGroup::Shared(crate::sort::SORT_COORD_GROUP)
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        // Retry any held output first (backpressure path).
        if !self.flush_held(ctx) {
            return Ok(StepOutcome::Contention);
        }

        // Drain any staged pending events before processing new input.
        if !self.pending.is_empty() {
            return Ok(self.emit_pending(ctx));
        }

        // If already finalized, drain pending or finish.
        if self.finalized {
            return Ok(StepOutcome::Finished);
        }

        // Try to pop and ingest one InflatedBlock.
        if let Some(block) = ctx.input.pop() {
            self.ingest_block(&block)?;
            // If ingest_block sealed a run it may have staged events; drain them.
            if !self.pending.is_empty() {
                return Ok(self.emit_pending(ctx));
            }
            return Ok(StepOutcome::Progress);
        }

        // No input available.
        if !ctx.input.is_drained() {
            return Ok(StepOutcome::NoProgress);
        }

        // Input fully drained: finalize (handles empty-input edge case).
        if let Some(first_event) = self.finalize()? {
            match ctx.outputs.push(first_event) {
                Ok(()) => {
                    // The push landed, so `held` is empty: it is safe to drain the
                    // next staged event this call (e.g. `AllAnnounced` after
                    // `Residual`).  `emit_pending` will itself hold on a full queue.
                    if !self.pending.is_empty() {
                        return Ok(self.emit_pending(ctx));
                    }
                }
                Err(unpushed) => {
                    // The output queue is full: hold `first_event` and retry it on
                    // the next `try_run` via `flush_held`.  Do NOT drain `pending`
                    // here — pushing a later event now would place it ahead of the
                    // still-held one, and a subsequent hold would clobber the held
                    // event (silent loss).  Keep the one-event-per-call invariant
                    // exactly like the normal path above.
                    self.held.put(unpushed);
                }
            }
            return Ok(StepOutcome::Progress);
        }

        // No arena was ingested (empty input).
        Ok(StepOutcome::Finished)
    }

    fn new_worker_copy(&self) -> Self {
        Self::new(self.strategy.fresh(), self.sort_threads, self.output_byte_limit)
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use fgumi_sort::{CoordinateChunkSorter, PooledSegmentedBuf, SegmentedBuf};
    use std::sync::Arc;

    // -----------------------------------------------------------------------
    // Test helpers shared across ReadBlocks and InflateToArena tests.
    // -----------------------------------------------------------------------

    /// Build a single BGZF block containing `payload` using `InlineBgzfCompressor`
    /// at compression level 0 (stored / uncompressed).
    /// Returns the raw block bytes (header + deflate payload + footer).
    fn make_test_bgzf_block(payload: &[u8]) -> Vec<u8> {
        let mut compressor = fgumi_bgzf::writer::InlineBgzfCompressor::new(0);
        compressor.write_all(payload).expect("compress payload");
        compressor.flush().expect("flush compressor");
        let mut blocks = compressor.take_blocks();
        assert_eq!(blocks.len(), 1, "payload must fit in one BGZF block");
        blocks.remove(0).data
    }

    /// Read the BGZF footer ISIZE field (last 4 bytes of the block), which is
    /// the uncompressed size mod 2^32.  Returns the value as `usize`.
    fn uncompressed_size_of(block: &[u8]) -> usize {
        assert!(block.len() >= 4, "block too short to contain BGZF footer");
        let n = block.len();
        u32::from_le_bytes([block[n - 4], block[n - 3], block[n - 2], block[n - 1]]) as usize
    }

    /// Build a single BGZF block containing `payload` using `InlineBgzfCompressor`
    /// at level 6, returning raw block bytes.  Used by `InflateToArena` tests.
    fn compress_one_block(payload: &[u8]) -> Vec<u8> {
        let mut compressor = fgumi_bgzf::writer::InlineBgzfCompressor::new(6);
        compressor.write_all(payload).expect("compress payload");
        compressor.flush().expect("flush compressor");
        let mut blocks = compressor.take_blocks();
        assert_eq!(blocks.len(), 1, "payload must fit in one BGZF block");
        blocks.remove(0).data
    }

    // -----------------------------------------------------------------------
    // ReadBlocks unit tests
    // -----------------------------------------------------------------------

    /// Two synthetic BGZF blocks are admitted; after freeze the emitted
    /// `ArenaBlock`s must carry contiguous offsets starting at `FRONT_REGION`,
    /// correct `len`/`ordinal` values, and share the same arena `Arc`.
    #[allow(unsafe_code)]
    #[test]
    fn read_blocks_admits_grows_and_emits_arena_blocks() {
        let blk0 = make_test_bgzf_block(b"first-block-decompressed-payload");
        let blk1 = make_test_bgzf_block(b"second-block-payload-2");
        let isz0 = uncompressed_size_of(&blk0);
        let isz1 = uncompressed_size_of(&blk1);

        // memory_limit larger than both blocks combined → single run
        let mut step = ReadBlocks::new(64 * 1024 * 1024, 64 * 1024 * 1024);
        assert!(
            step.admit_block(BgzfBlock {
                batch_serial: 0,
                bytes: blk0.clone(),
                uncompressed_size: u32::try_from(isz0).unwrap(),
            })
            .is_ok()
        );
        assert!(
            step.admit_block(BgzfBlock {
                batch_serial: 1,
                bytes: blk1.clone(),
                uncompressed_size: u32::try_from(isz1).unwrap(),
            })
            .is_ok()
        );
        let emitted: Vec<ArenaBlock> = step.seal_and_drain_for_test();

        assert_eq!(emitted.len(), 2);
        // Block 0 lands at FRONT_REGION (after the reserved front region).
        assert_eq!(emitted[0].offset, FRONT_REGION as u64, "block 0 must start at FRONT_REGION");
        // Block 1 is contiguous after block 0.
        assert_eq!(
            usize::try_from(emitted[1].offset).unwrap(),
            FRONT_REGION + isz0,
            "offsets must be contiguous after FRONT_REGION"
        );
        assert!(emitted[1].is_last_of_run, "last block must carry is_last_of_run=true");
        assert!(!emitted[0].is_last_of_run, "first block must not carry is_last_of_run");
        assert_eq!(emitted[0].len as usize, isz0);
        assert_eq!(emitted[1].len as usize, isz1);
        assert_eq!(emitted[0].ordinal, 0);
        assert_eq!(emitted[1].ordinal, 1);
        assert_eq!(emitted[0].run_seq, 0);
        assert_eq!(emitted[1].run_seq, 0);
        assert!(!emitted[0].seals_to_spill, "residual run must not seal to spill");
        assert!(!emitted[1].seals_to_spill, "residual run must not seal to spill");
        // All blocks share the same arena Arc.
        assert!(Arc::ptr_eq(&emitted[0].arena, &emitted[1].arena));
    }

    /// Feed enough blocks to force TWO runs.
    ///
    /// Asserts:
    /// - Run 0 blocks: `run_seq == 0`, offsets start at `FRONT_REGION`,
    ///   last block has `is_last_of_run`, all have `seals_to_spill = true`.
    /// - Run 1 blocks: `run_seq == 1`, offsets restart at `FRONT_REGION`
    ///   (fresh/reused arena from the pool), `seals_to_spill = false`.
    /// - Ordinals are monotonically increasing across both runs.
    /// - Run 0 and run 1 arena `Arc`s differ (or are the same reused physical
    ///   buffer returned via the pool after run 0's Arc drops).
    #[allow(unsafe_code)]
    #[allow(clippy::too_many_lines)] // exhaustive two-run state-machine assertions
    #[test]
    fn read_blocks_two_run_seal() {
        // Make two blocks; each has isz bytes of uncompressed data.
        let payload0 = b"run0-block0-payload".as_ref();
        let payload1 = b"run0-block1-payload".as_ref();
        let payload2 = b"run1-block0-payload".as_ref();
        let payload3 = b"run1-block1-payload".as_ref();

        let blk0 = make_test_bgzf_block(payload0);
        let blk1 = make_test_bgzf_block(payload1);
        let blk2 = make_test_bgzf_block(payload2);
        let blk3 = make_test_bgzf_block(payload3);

        let isz0 = uncompressed_size_of(&blk0);
        let isz1 = uncompressed_size_of(&blk1);
        let isz2 = uncompressed_size_of(&blk2);
        let isz3 = uncompressed_size_of(&blk3);

        // memory_limit set to force a seal after run 0 (blk0+blk1):
        // use isz0+isz1 as the cap so that after admitting blk1, run_cumulative
        // >= run_cap and a seal fires.  The arena segment is derived from this
        // budget (FRONT_REGION + run_cap + one block), so all tiny test blocks
        // still fit in one segment per run.
        let run_cap = isz0 + isz1;
        let mut step = ReadBlocks::new(run_cap, 64 * 1024 * 1024);

        // Admit run 0 blocks.
        assert!(
            step.admit_block(BgzfBlock {
                batch_serial: 0,
                bytes: blk0.clone(),
                uncompressed_size: u32::try_from(isz0).unwrap(),
            })
            .is_ok()
        );
        assert!(
            step.admit_block(BgzfBlock {
                batch_serial: 1,
                bytes: blk1.clone(),
                uncompressed_size: u32::try_from(isz1).unwrap(),
            })
            .is_ok()
        );

        // After run 0's blocks fill run_cap, seal run 0 explicitly and collect its
        // blocks (simulating the seal that happens when run_cumulative >= run_cap).
        // In the pipeline try_run, the seal fires after admit; here we do it manually
        // via a dedicated seal so we can inspect both runs independently.
        step.seal_run(true); // mid-stream → seals_to_spill = true

        // Collect run 0 blocks.
        let run0_blocks: Vec<ArenaBlock> = step.emit.drain(..).collect();
        assert_eq!(run0_blocks.len(), 2, "run 0 must emit 2 blocks");

        // Run 0 invariants.  Eager emission stamps `is_last_of_run` /
        // `seals_to_spill` only on the run's LAST block (the withheld one); non-last
        // blocks carry `false` (downstream reads these fields only on the last block).
        for blk in &run0_blocks {
            assert_eq!(blk.run_seq, 0, "run 0 blocks must have run_seq=0");
        }
        assert!(!run0_blocks[0].is_last_of_run, "run 0 block 0 must not be last_of_run");
        assert!(run0_blocks[1].is_last_of_run, "run 0 block 1 must be last_of_run");
        assert!(!run0_blocks[0].seals_to_spill, "non-last block carries seals_to_spill=false");
        assert!(run0_blocks[1].seals_to_spill, "run 0 last block must seal to spill");
        assert_eq!(
            run0_blocks[0].offset, FRONT_REGION as u64,
            "run 0 block 0 must start at FRONT_REGION"
        );
        assert_eq!(
            usize::try_from(run0_blocks[1].offset).unwrap(),
            FRONT_REGION + isz0,
            "run 0 block 1 must be contiguous after block 0"
        );
        assert_eq!(run0_blocks[0].ordinal, 0, "global ordinal monotonic");
        assert_eq!(run0_blocks[1].ordinal, 1, "global ordinal monotonic");

        // Verify pool exhaustion: while run 0's Arc is still alive, admitting a new
        // block must fail (pool capacity 1 → all arenas in flight).
        {
            let probe_blk = make_test_bgzf_block(b"pool-exhaustion-probe");
            let probe_isz = uncompressed_size_of(&probe_blk);
            let result = step.admit_block(BgzfBlock {
                batch_serial: 99,
                bytes: probe_blk,
                uncompressed_size: u32::try_from(probe_isz).unwrap(),
            });
            assert!(result.is_err(), "pool must be exhausted while run 0 Arc is still in flight");
            // A failed admit (ensure_arena returned false) must not leak state: the
            // block is handed back in the Err, no block is withheld, arena stays None.
            assert!(step.last_block.is_none(), "failed admit must not withhold a block");
            assert!(step.arena.is_none(), "failed admit must leave arena as None");
        }

        // Drop run 0 blocks → their Arc<PooledSegmentedBuf> refcount falls to 0 →
        // PooledSegmentedBuf::drop releases the arena back to the pool.
        drop(run0_blocks);

        // Now admit run 1 blocks — the pool has one free (reused) arena.
        assert!(
            step.admit_block(BgzfBlock {
                batch_serial: 2,
                bytes: blk2.clone(),
                uncompressed_size: u32::try_from(isz2).unwrap(),
            })
            .is_ok(),
            "run 1 admit must succeed after pool release"
        );
        assert!(
            step.admit_block(BgzfBlock {
                batch_serial: 3,
                bytes: blk3.clone(),
                uncompressed_size: u32::try_from(isz3).unwrap(),
            })
            .is_ok(),
            "run 1 second admit must succeed"
        );

        // Seal run 1 as the residual.
        step.seal_run(false);
        let run1_blocks: Vec<ArenaBlock> = step.emit.drain(..).collect();
        assert_eq!(run1_blocks.len(), 2, "run 1 must emit 2 blocks");

        // Run 1 invariants.
        for blk in &run1_blocks {
            assert_eq!(blk.run_seq, 1, "run 1 blocks must have run_seq=1");
            assert!(!blk.seals_to_spill, "run 1 (residual) blocks must have seals_to_spill=false");
        }
        assert!(!run1_blocks[0].is_last_of_run, "run 1 block 0 must not be last_of_run");
        assert!(run1_blocks[1].is_last_of_run, "run 1 block 1 must be last_of_run");
        // Run 1 offsets restart at FRONT_REGION (fresh/reused arena).
        assert_eq!(
            run1_blocks[0].offset, FRONT_REGION as u64,
            "run 1 block 0 must restart at FRONT_REGION"
        );
        assert_eq!(
            usize::try_from(run1_blocks[1].offset).unwrap(),
            FRONT_REGION + isz2,
            "run 1 block 1 must be contiguous after block 0"
        );
        // Ordinals are globally monotonic across runs.
        assert_eq!(run1_blocks[0].ordinal, 2, "global ordinal monotonic across runs");
        assert_eq!(run1_blocks[1].ordinal, 3, "global ordinal monotonic across runs");

        // All run 1 blocks share the same arena Arc.
        assert!(Arc::ptr_eq(&run1_blocks[0].arena, &run1_blocks[1].arena));

        // Pool behavioral check: the run 1 arena was acquired after the pool round-
        // tripped run 0's arena back.  Verify the underlying storage was reused
        // (not a fresh allocation) by comparing the allocated_capacity of the run 1
        // arena: the reused arena retains its segment allocations (reset_for_reuse
        // keeps the Vec capacity), so its allocated_capacity is >= the derived
        // segment size (FRONT_REGION + run_cap + one block, pre-reserved by
        // reserve_full_capacity in ensure_arena).
        let expected_segment_size = FRONT_REGION + run_cap + MAX_BGZF_BLOCK;
        let run1_alloc_cap = run1_blocks[0].arena.allocated_capacity();
        assert!(
            run1_alloc_cap >= expected_segment_size,
            "reused arena must retain its segment capacity (got {run1_alloc_cap}, expected >= {expected_segment_size})"
        );
    }

    /// Build a minimal valid BAM record body (no `block_size` prefix) with the given
    /// ref\_id/pos and a one-character read name `name`.  Returns the body bytes.
    ///
    /// Copied verbatim from `crates/fgumi-sort/src/ref_sort.rs`'s test module
    /// (the `(i32, i32, u8)` version) so tests in this crate can construct
    /// identical records without a cross-crate dependency on a `#[cfg(test)]`
    /// helper.
    fn coord_body(ref_id: i32, pos: i32, name: u8) -> Vec<u8> {
        let mut b = Vec::new();
        b.extend_from_slice(&ref_id.to_le_bytes());
        b.extend_from_slice(&pos.to_le_bytes());
        b.push(2); // l_read_name (incl NUL): "<name>\0"
        b.push(0); // mapq
        b.extend_from_slice(&0u16.to_le_bytes()); // bin
        b.extend_from_slice(&0u16.to_le_bytes()); // n_cigar_op
        b.extend_from_slice(&0u16.to_le_bytes()); // flag (forward)
        b.extend_from_slice(&0u32.to_le_bytes()); // l_seq
        b.extend_from_slice(&(-1i32).to_le_bytes()); // next_refID
        b.extend_from_slice(&(-1i32).to_le_bytes()); // next_pos
        b.extend_from_slice(&0i32.to_le_bytes()); // tlen
        b.push(name); // read_name char
        b.push(0); // read_name NUL terminator
        b
    }

    /// Build a minimal BAM binary header parseable by `bam_header_len`.
    ///
    /// Layout: `magic(4) + l_text=0(4) + n_ref(4)` followed by `n_ref` entries of
    /// `l_name=2(4) + "r\0"(2) + l_ref=100(4)`.
    fn minimal_bam_header(n_ref: u32) -> Vec<u8> {
        let mut h = Vec::new();
        h.extend_from_slice(b"BAM\x01"); // magic
        h.extend_from_slice(&0u32.to_le_bytes()); // l_text = 0
        h.extend_from_slice(&n_ref.to_le_bytes()); // n_ref
        for _ in 0..n_ref {
            h.extend_from_slice(&2u32.to_le_bytes()); // l_name = 2
            h.push(b'r');
            h.push(0); // name "r\0"
            h.extend_from_slice(&100u32.to_le_bytes()); // l_ref = 100
        }
        h
    }

    /// Oracle test: `FindBoundariesAndSort` emitted chunk must be byte-identical to
    /// the copy-based `CoordinateChunkSorter` over the same record bodies.
    #[allow(unsafe_code)]
    #[test]
    fn find_boundaries_and_sort_single_run_matches_oracle() {
        let n_ref = 4u32;
        // Records deliberately out of coordinate order; distinct names witness stability.
        let recs = vec![
            coord_body(2, 100, b'a'),
            coord_body(0, 50, b'b'),
            coord_body(2, 10, b'c'),
            coord_body(0, 50, b'd'), // coordinate tie with b'b' → stable sort keeps b before d
            coord_body(1, 999, b'e'),
        ];

        // Build an arena: [BAM header][[block_size][body]...] as one contiguous run.
        let header = minimal_bam_header(n_ref);
        let mut arena = SegmentedBuf::with_capacity(0, 1 << 20);
        arena.reserve_full_capacity();

        // SAFETY: every slot fully written before any read.
        let h_off = unsafe { arena.grow_uninit(header.len()) };
        unsafe { arena.slice_mut(h_off, header.len()) }.copy_from_slice(&header);
        // h_off is usize; on 64-bit targets usize fits in u64.
        #[allow(clippy::cast_possible_truncation)]
        let run_start = h_off as u64;

        for r in &recs {
            let bs = u32::try_from(r.len()).unwrap();
            let po = unsafe { arena.grow_uninit(4) };
            unsafe { arena.slice_mut(po, 4) }.copy_from_slice(&bs.to_le_bytes());
            let bo = unsafe { arena.grow_uninit(r.len()) };
            unsafe { arena.slice_mut(bo, r.len()) }.copy_from_slice(r);
        }
        // arena.len() is usize; usize fits in u64 on 64-bit targets.
        #[allow(clippy::cast_possible_truncation)]
        let run_len = arena.len() as u64 - run_start;
        let arena = Arc::new(PooledSegmentedBuf::unpooled(arena));

        // Drive ingest_block + finalize directly (the test seam).
        let mut step =
            FindBoundariesAndSort::new(CoordinateStrategy::new(n_ref), 1, 64 * 1024 * 1024);
        step.ingest_block(&InflatedBlock {
            arena: Arc::clone(&arena),
            ordinal: 0,
            offset: run_start,
            len: u32::try_from(run_len).unwrap(),
            is_last_of_run: true,
            run_seq: 0,
            seals_to_spill: false,
        })
        .expect("ingest_block must succeed");
        let ev = step.finalize().expect("finalize must succeed").expect("residual event");

        let chunk_bytes: Vec<Vec<u8>> = match ev {
            SortChunkEvent::Residual { chunk: MemoryChunkErased::Coordinate(c), .. } => {
                (0..c.len()).map(|i| c.record_bytes(i).to_vec()).collect()
            }
            _ => panic!("expected Residual coordinate chunk, got something else"),
        };

        // Oracle: copy-based in-memory coordinate sort of the same record bodies.
        let mut oracle = CoordinateChunkSorter::for_test(usize::MAX, n_ref);
        for r in &recs {
            let _ = oracle.push(r).unwrap();
        }
        let oc = oracle.take_sorted_chunk();
        let oracle_bytes: Vec<Vec<u8>> =
            (0..oc.len()).map(|i| oc.record_bytes(i).to_vec()).collect();

        assert_eq!(
            chunk_bytes, oracle_bytes,
            "arena-scan single-run sort must be byte-identical to the copy-based sorter"
        );
    }

    /// Two-run straddler test: a record that spans the boundary between run 0 and
    /// run 1 must appear as record 0 of run 1's chunk with bytes equal to
    /// `carry ++ head` (the full original record), and the merged union of both
    /// runs (stable, lower-run first) must equal the oracle (`CoordinateChunkSorter`
    /// over the same records in input order).
    ///
    /// Test layout:
    /// - Run 0 arena (no `FRONT_REGION` prefix — built manually like the oracle test):
    ///   `[header][rec0_prefix+body][rec1_prefix+body][straddler_prefix+partial_body]`
    ///   The straddler record's `block_size` prefix PLUS the first few body bytes land
    ///   in run 0; the rest of the body starts at `FRONT_REGION` in run 1's arena.
    /// - Run 1 arena: `FRONT_REGION` uninit bytes reserved, then `straddler_tail +
    ///   rec2_prefix + body`.
    ///   `FindBoundariesAndSort` writes the carry into `[FRONT_REGION - L, FRONT_REGION)`,
    ///   so `carry ++ straddler_tail` is contiguous and forms the complete record.
    #[allow(unsafe_code)]
    #[allow(clippy::too_many_lines)]
    #[test]
    fn find_boundaries_and_sort_two_run_straddler() {
        let n_ref = 2u32;
        let seg_size = 256 * 1024 * 1024usize;

        // Records: rec0 + rec1 go into run 0 (complete); straddler goes across
        // the boundary; rec2 goes into run 1 (complete).
        let rec0 = coord_body(0, 10, b'a');
        let rec1 = coord_body(1, 20, b'b');
        let straddler = coord_body(0, 5, b'c'); // will be straddled across runs
        let rec2 = coord_body(1, 5, b'd');

        let all_recs = vec![&rec0, &rec1, &straddler, &rec2];

        // -----------------------------------------------------------------------
        // Build run 0 arena: [header][rec0][rec1][straddler_prefix+partial_body]
        // We use the same manual layout as the oracle test (no FRONT_REGION prefix).
        // -----------------------------------------------------------------------
        let header = minimal_bam_header(n_ref);

        let mut arena0 = SegmentedBuf::with_capacity(0, seg_size);
        arena0.reserve_full_capacity();

        // Write header.
        let h_off = unsafe { arena0.grow_uninit(header.len()) };
        unsafe { arena0.slice_mut(h_off, header.len()) }.copy_from_slice(&header);

        // Write rec0 (prefix + body) — complete.
        let bs0 = u32::try_from(rec0.len()).unwrap();
        let po = unsafe { arena0.grow_uninit(4) };
        unsafe { arena0.slice_mut(po, 4) }.copy_from_slice(&bs0.to_le_bytes());
        let bo = unsafe { arena0.grow_uninit(rec0.len()) };
        unsafe { arena0.slice_mut(bo, rec0.len()) }.copy_from_slice(&rec0);

        // Write rec1 (prefix + body) — complete.
        let bs1 = u32::try_from(rec1.len()).unwrap();
        let po = unsafe { arena0.grow_uninit(4) };
        unsafe { arena0.slice_mut(po, 4) }.copy_from_slice(&bs1.to_le_bytes());
        let bo = unsafe { arena0.grow_uninit(rec1.len()) };
        unsafe { arena0.slice_mut(bo, rec1.len()) }.copy_from_slice(&rec1);

        // Write straddler prefix (4-byte block_size) + partial body.
        // Split: put the 4-byte prefix + first 3 bytes of body in run 0.
        let bs_str = u32::try_from(straddler.len()).unwrap();
        let partial_len = 3usize; // bytes of straddler body in run 0
        assert!(partial_len < straddler.len(), "partial_len must be < straddler body size");
        let po = unsafe { arena0.grow_uninit(4) };
        unsafe { arena0.slice_mut(po, 4) }.copy_from_slice(&bs_str.to_le_bytes());
        let bp = unsafe { arena0.grow_uninit(partial_len) };
        unsafe { arena0.slice_mut(bp, partial_len) }.copy_from_slice(&straddler[..partial_len]);

        // run 0 span: from header offset (0) to end of arena.
        let run0_start = h_off as u64;
        let run0_end = arena0.len() as u64;
        let run0_len = usize::try_from(run0_end - run0_start).unwrap();

        let arena0 = Arc::new(PooledSegmentedBuf::unpooled(arena0));

        // -----------------------------------------------------------------------
        // Build run 1 arena: [FRONT_REGION uninit][straddler_tail][rec2 prefix+body]
        // -----------------------------------------------------------------------
        let mut arena1 = SegmentedBuf::with_capacity(0, seg_size);
        arena1.reserve_full_capacity();

        // Reserve the FRONT_REGION prefix (uninit — the carry will be written here
        // by FindBoundariesAndSort).
        let _front = unsafe { arena1.grow_uninit(FRONT_REGION) };
        assert_eq!(arena1.len(), FRONT_REGION);

        // Write straddler tail (remaining body bytes) at FRONT_REGION.
        let tail_len = straddler.len() - partial_len;
        let st_off = unsafe { arena1.grow_uninit(tail_len) };
        assert_eq!(st_off, FRONT_REGION, "straddler tail must start at FRONT_REGION");
        unsafe { arena1.slice_mut(st_off, tail_len) }.copy_from_slice(&straddler[partial_len..]);

        // Write rec2 (prefix + body) — complete.
        let bs2 = u32::try_from(rec2.len()).unwrap();
        let po = unsafe { arena1.grow_uninit(4) };
        unsafe { arena1.slice_mut(po, 4) }.copy_from_slice(&bs2.to_le_bytes());
        let bo = unsafe { arena1.grow_uninit(rec2.len()) };
        unsafe { arena1.slice_mut(bo, rec2.len()) }.copy_from_slice(&rec2);

        // run 1 spans from FRONT_REGION (the inflate data, not the front region).
        let run1_data_start = FRONT_REGION as u64;
        let run1_data_end = arena1.len() as u64;
        let run1_data_len = usize::try_from(run1_data_end - run1_data_start).unwrap();

        let arena1 = Arc::new(PooledSegmentedBuf::unpooled(arena1));

        // -----------------------------------------------------------------------
        // Drive FindBoundariesAndSort across both runs.
        // -----------------------------------------------------------------------
        let mut fbs =
            FindBoundariesAndSort::new(CoordinateStrategy::new(n_ref), 1, 64 * 1024 * 1024);

        // Run 0: one block, is_last_of_run = true, seals_to_spill = true.
        fbs.ingest_block(&InflatedBlock {
            arena: Arc::clone(&arena0),
            ordinal: 0,
            offset: run0_start,
            len: u32::try_from(run0_len).unwrap(),
            is_last_of_run: true,
            run_seq: 0,
            seals_to_spill: true,
        })
        .expect("run 0 ingest must succeed");

        // After run 0's last block, a Spill event must be staged.
        assert_eq!(fbs.pending.len(), 1, "run 0 must stage exactly one Spill event");
        let ev0 = fbs.pending.pop_front().unwrap();
        let run0_chunk = match ev0 {
            SortChunkEvent::Spill { seq, chunk: MemoryChunkErased::Coordinate(c), .. } => {
                assert_eq!(seq, 0, "first Spill must have seq=0");
                c
            }
            _ => panic!("expected Spill(Coordinate) for run 0, got something else"),
        };

        // Assert run 0's carry: 4 bytes (prefix) + partial_len body bytes.
        let expected_carry_len = 4 + partial_len;
        assert_eq!(
            fbs.carry.len(),
            expected_carry_len,
            "carry after run 0 must be {expected_carry_len} bytes (prefix + partial body)"
        );

        // Run 0 chunk must NOT contain the straddler — only rec0 and rec1.
        assert_eq!(
            run0_chunk.len(),
            2,
            "run 0 chunk must have exactly 2 complete records (rec0, rec1)"
        );

        // Run 1: one block, is_last_of_run = true, seals_to_spill = false.
        fbs.ingest_block(&InflatedBlock {
            arena: Arc::clone(&arena1),
            ordinal: 1,
            offset: run1_data_start,
            len: u32::try_from(run1_data_len).unwrap(),
            is_last_of_run: true,
            run_seq: 1,
            seals_to_spill: false,
        })
        .expect("run 1 ingest must succeed");

        // After run 1's last block, pending has Residual + AllAnnounced.
        assert_eq!(fbs.pending.len(), 2, "run 1 must stage Residual + AllAnnounced");
        assert!(fbs.carry.is_empty(), "carry must be empty after the final run");

        let ev_residual = fbs.pending.pop_front().unwrap();
        let SortChunkEvent::Residual { chunk: MemoryChunkErased::Coordinate(run1_chunk), .. } =
            ev_residual
        else {
            panic!("expected Residual(Coordinate) for run 1, got something else")
        };

        let ev_announced = fbs.pending.pop_front().unwrap();
        match ev_announced {
            SortChunkEvent::AllAnnounced { slot_count, memory_chunk_count, total_records } => {
                assert_eq!(slot_count, 1, "AllAnnounced: slot_count must be 1 (one spill)");
                assert_eq!(memory_chunk_count, 1, "AllAnnounced: memory_chunk_count must be 1");
                assert_eq!(total_records, 4, "AllAnnounced: total_records must be 4");
            }
            _ => panic!("expected AllAnnounced, got something else"),
        }

        // Run 1 chunk: straddler (record 0) + rec2 (record 1) — 2 records.
        assert_eq!(run1_chunk.len(), 2, "run 1 chunk must have 2 records (straddler + rec2)");

        // Assert straddler is record 0 of run 1's sorted chunk.
        // The straddler has key (ref_id=0, pos=5) and rec2 has (ref_id=1, pos=5);
        // coordinate sort orders by ref_id first, so straddler (ref_id=0) comes before
        // rec2 (ref_id=1).
        let straddler_bytes = run1_chunk.record_bytes(0).to_vec();
        assert_eq!(
            straddler_bytes, straddler,
            "run 1 record 0 must be the straddler (carry ++ head = full original record)"
        );
        let rec2_bytes = run1_chunk.record_bytes(1).to_vec();
        assert_eq!(rec2_bytes, rec2, "run 1 record 1 must be rec2");

        // -----------------------------------------------------------------------
        // Oracle check: stable merge of run0 + run1 by coordinate key, lower
        // source index first (run0 < run1) must equal the oracle over all 4 records.
        //
        // Oracle sort order: (ref_id=0,pos=5)=straddler, (ref_id=0,pos=10)=rec0,
        //   (ref_id=1,pos=5)=rec2, (ref_id=1,pos=20)=rec1.
        //
        // Merge gives: run0 records sorted = [rec0(0,10), rec1(1,20)];
        //              run1 records sorted = [straddler(0,5), rec2(1,5)].
        // Interleaved lower-run-first stable merge:
        //   Compare run0[0]=(0,10) vs run1[0]=(0,5): run1 wins → straddler
        //   Compare run0[0]=(0,10) vs run1[1]=(1,5): run0 wins → rec0
        //   Compare run0[1]=(1,20) vs run1[1]=(1,5): run1 wins → rec2
        //   run1 exhausted → run0[1] = rec1
        // Merge result: [straddler, rec0, rec2, rec1]
        // -----------------------------------------------------------------------
        let mut oracle = CoordinateChunkSorter::for_test(usize::MAX, n_ref);
        for r in &all_recs {
            let _ = oracle.push(r).unwrap();
        }
        let oracle_chunk = oracle.take_sorted_chunk();
        let oracle_bytes: Vec<Vec<u8>> =
            (0..oracle_chunk.len()).map(|i| oracle_chunk.record_bytes(i).to_vec()).collect();
        // Expected oracle order: [(0,5)=straddler, (0,10)=rec0, (1,5)=rec2, (1,20)=rec1]
        assert_eq!(oracle_bytes[0], straddler, "oracle[0] must be straddler");
        assert_eq!(oracle_bytes[1], rec0, "oracle[1] must be rec0");
        assert_eq!(oracle_bytes[2], rec2, "oracle[2] must be rec2");
        assert_eq!(oracle_bytes[3], rec1, "oracle[3] must be rec1");

        // Collect sorted bytes from both run chunks (run0 first, then run1) via a
        // manual stable merge that mirrors the merge engine's lower-source-index-first
        // tie-break.
        //
        // run0 sorted: [rec0(0,10), rec1(1,20)]
        // run1 sorted: [straddler(0,5), rec2(1,5)]
        let run0_sorted: Vec<Vec<u8>> =
            (0..run0_chunk.len()).map(|i| run0_chunk.record_bytes(i).to_vec()).collect();
        let run1_sorted: Vec<Vec<u8>> =
            (0..run1_chunk.len()).map(|i| run1_chunk.record_bytes(i).to_vec()).collect();

        // Build a merged sequence via coordinate key comparison.
        // We use the oracle order as the expected merged order (they must match).
        let mut merged: Vec<Vec<u8>> = Vec::with_capacity(4);
        let mut i0 = 0usize;
        let mut i1 = 0usize;
        while i0 < run0_sorted.len() || i1 < run1_sorted.len() {
            if i0 >= run0_sorted.len() {
                merged.push(run1_sorted[i1].clone());
                i1 += 1;
            } else if i1 >= run1_sorted.len() {
                merged.push(run0_sorted[i0].clone());
                i0 += 1;
            } else {
                // Extract coordinate key: ref_id (first i32) and pos (second i32).
                let key0 = {
                    let b = &run0_sorted[i0];
                    let r = i32::from_le_bytes(b[0..4].try_into().unwrap());
                    let p = i32::from_le_bytes(b[4..8].try_into().unwrap());
                    (r, p)
                };
                let key1 = {
                    let b = &run1_sorted[i1];
                    let r = i32::from_le_bytes(b[0..4].try_into().unwrap());
                    let p = i32::from_le_bytes(b[4..8].try_into().unwrap());
                    (r, p)
                };
                // Stable: on tie, run0 (lower source index) wins.
                if key1 < key0 {
                    merged.push(run1_sorted[i1].clone());
                    i1 += 1;
                } else {
                    merged.push(run0_sorted[i0].clone());
                    i0 += 1;
                }
            }
        }

        assert_eq!(
            merged, oracle_bytes,
            "stable merge of run0 + run1 chunks must equal the oracle over all 4 records"
        );
    }

    #[allow(unsafe_code)]
    #[test]
    fn inflate_writes_decompressed_bytes_into_the_slot() {
        // Build a payload that is comfortably under one BGZF block (< 64 KiB).
        let payload = b"PARALLEL-INFLATE-ARENA-TEST".repeat(50);
        let block = compress_one_block(&payload);
        let isize = u32::try_from(payload.len()).unwrap();

        // Construct an arena with enough capacity for the payload, reserve it,
        // then carve out an uninit slot for the inflate worker to fill.
        let mut arena = SegmentedBuf::with_capacity(0, 1 << 20);
        arena.reserve_full_capacity();
        // SAFETY: slot is fully written by `inflate_one` below before any read.
        let offset = unsafe { arena.grow_uninit(payload.len()) } as u64;
        let arena = Arc::new(PooledSegmentedBuf::unpooled(arena));

        let item = ArenaBlock {
            arena: Arc::clone(&arena),
            ordinal: 0,
            offset,
            len: isize,
            block,
            is_last_of_run: true,
            run_seq: 0,
            seals_to_spill: false,
        };

        let mut step = InflateToArena::new(64 * 1024 * 1024);
        let inflated = step.inflate_one(item).expect("inflate must succeed");

        assert_eq!(inflated.offset, offset, "offset must be forwarded unchanged");
        assert_eq!(inflated.len, isize, "len must be forwarded unchanged");
        assert!(!inflated.seals_to_spill, "seals_to_spill must be forwarded");
        // Confirm the decompressed bytes are in the arena at the reserved slot.
        assert_eq!(
            arena.slice(
                usize::try_from(offset).unwrap(),
                isize as usize, // u32 always fits in usize
            ),
            &payload[..],
            "arena slot must contain the original payload after inflate"
        );
    }

    /// Regression coverage for the `FindBoundariesAndSort` `try_run`
    /// finalize/drained-branch held-overwrite bug.
    ///
    /// The bug: on the drained path `try_run` did `push(first_event)` and then,
    /// in the SAME call, unconditionally drained the rest of `pending` via
    /// `emit_pending` WITHOUT re-checking `held`.  If the first push is rejected
    /// (full output queue → `first_event` goes into `held`), the follow-on
    /// `emit_pending` would push a LATER event (`AllAnnounced`) ahead of the
    /// still-held `Residual` and then `held.put(...)` it — clobbering the held
    /// `Residual` (`HeldSlot::put` asserts on a double-put, so it would panic /
    /// lose the record).  The fix drains `pending` only when the first push
    /// SUCCEEDED, preserving the one-event-per-call discipline of the normal path.
    ///
    /// This test pins the ordering contract `finalize()` hands `try_run`, which is
    /// exactly what the one-event-per-call fix must emit in order: `Residual`
    /// FIRST, then `AllAnnounced`, each staged exactly once.  If `finalize()`
    /// staged them in the wrong order (or duplicated one), the held/retry path in
    /// `try_run` could not emit a correct stream no matter how careful it is.
    ///
    /// Coverage limit: the held-overwrite lived inside `try_run`, which needs a
    /// backpressuring `StepCtx` (an `OutputHandles<Single<_>>` whose queue rejects
    /// the first push).  `OutputHandles::new` is `pub(crate)` to
    /// `fgumi-pipeline-core`, so a rejecting output context cannot be built from
    /// this crate without editing another crate; the full `try_run` held/retry path
    /// is exercised end-to-end by the framework-driven tests in `sort/tests.rs`.
    /// Here we assert the finalize-branch invariant as directly as the in-crate
    /// seams (`finalize`, `pending`) allow.
    #[allow(unsafe_code)]
    #[test]
    fn finalize_stages_residual_before_all_announced_each_once() {
        let n_ref = 2u32;
        // A complete run (header + two whole records, no trailing partial), so
        // sealing it as a residual leaves an empty carry.
        let recs = [coord_body(0, 10, b'a'), coord_body(1, 20, b'b')];
        let header = minimal_bam_header(n_ref);
        let mut arena = SegmentedBuf::with_capacity(0, 1 << 20);
        arena.reserve_full_capacity();
        // SAFETY: every slot fully written before any read.
        let h_off = unsafe { arena.grow_uninit(header.len()) };
        unsafe { arena.slice_mut(h_off, header.len()) }.copy_from_slice(&header);
        #[allow(clippy::cast_possible_truncation)]
        let run_start = h_off as u64;
        for r in &recs {
            let bs = u32::try_from(r.len()).unwrap();
            let po = unsafe { arena.grow_uninit(4) };
            unsafe { arena.slice_mut(po, 4) }.copy_from_slice(&bs.to_le_bytes());
            let bo = unsafe { arena.grow_uninit(r.len()) };
            unsafe { arena.slice_mut(bo, r.len()) }.copy_from_slice(r);
        }
        #[allow(clippy::cast_possible_truncation)]
        let run_len = arena.len() as u64 - run_start;
        let arena = Arc::new(PooledSegmentedBuf::unpooled(arena));

        let mut step =
            FindBoundariesAndSort::new(CoordinateStrategy::new(n_ref), 1, 64 * 1024 * 1024);
        // Ingest the run WITHOUT `is_last_of_run`, so `ingest_block` does NOT seal:
        // the seal (and the `Residual` + `AllAnnounced` staging) then happens in
        // `finalize()` — the path `try_run`'s drained branch drives.
        step.ingest_block(&InflatedBlock {
            arena: Arc::clone(&arena),
            ordinal: 0,
            offset: run_start,
            len: u32::try_from(run_len).unwrap(),
            is_last_of_run: false,
            run_seq: 0,
            seals_to_spill: false,
        })
        .expect("ingest_block must succeed");
        assert!(step.pending.is_empty(), "no events staged before finalize");

        // finalize() returns the FIRST event (Residual) and leaves the remainder
        // staged in `pending`.  This is the one-event-at-a-time hand-off that
        // `try_run`'s drained branch must respect: emit `first_event`, and only
        // drain `pending` if that push landed.
        let first = step.finalize().expect("finalize must succeed").expect("first event present");
        assert!(
            matches!(first, SortChunkEvent::Residual { .. }),
            "finalize must return Residual as the first event"
        );
        // Exactly one event remains staged, and it is AllAnnounced (never emitted
        // ahead of the Residual).
        assert_eq!(step.pending.len(), 1, "exactly one event remains staged after the Residual");
        match step.pending.front().expect("AllAnnounced staged") {
            SortChunkEvent::AllAnnounced { slot_count, memory_chunk_count, .. } => {
                assert_eq!(*slot_count, 0, "no spills → slot_count 0");
                assert_eq!(*memory_chunk_count, 1, "one in-memory residual chunk");
            }
            _ => panic!("second staged event must be AllAnnounced"),
        }

        // A subsequent finalize() (mirroring a later drained `try_run` after the
        // Residual flushed) hands back AllAnnounced, then nothing — proving each
        // event is produced exactly once and in order.
        let second = step.finalize().expect("finalize must succeed").expect("AllAnnounced present");
        assert!(
            matches!(second, SortChunkEvent::AllAnnounced { .. }),
            "second finalize must return AllAnnounced"
        );
        assert!(step.pending.is_empty(), "no further events staged");
        assert!(
            step.finalize().expect("finalize must succeed").is_none(),
            "no more events after Residual + AllAnnounced"
        );
    }

    /// Runtime proof that `--sort-threads` (Phase 1) controls the actual sort
    /// worker count, not just that it parses.
    ///
    /// The chain builder resolves `--sort-threads` (falling back to `--threads`)
    /// into `num_phase1_threads` and hands it to `FindBoundariesAndSort::new`,
    /// which forwards it to `strategy.seal(arena, self.sort_threads)`. The
    /// queryname strategy builds a bounded rayon pool sized to that value and runs
    /// the per-run comparator sort inside it — so the pool's thread count IS the
    /// effective Phase-1 concurrency. (The D1.1 regression handed the strategy the
    /// raw global `--threads` instead, silently ignoring `--sort-threads`.)
    ///
    /// `ThreadPool::current_num_threads()` is the fixed pool size, so this is
    /// deterministic — unlike counting how many workers a given input happens to
    /// keep busy. Sort output is byte-identical across thread counts, so this pool
    /// observation is the only way to assert the flag's runtime effect.
    #[allow(unsafe_code)]
    #[test]
    fn phase1_sort_threads_sizes_the_queryname_worker_pool() {
        use fgumi_sort::RawQuerynameLexKey;

        for sort_threads in [1usize, 3] {
            let n_ref = 2u32;
            let recs = [coord_body(0, 10, b'a'), coord_body(1, 20, b'b')];
            let header = minimal_bam_header(n_ref);
            let mut arena = SegmentedBuf::with_capacity(0, 1 << 20);
            arena.reserve_full_capacity();
            // SAFETY: every slot fully written before any read.
            let h_off = unsafe { arena.grow_uninit(header.len()) };
            unsafe { arena.slice_mut(h_off, header.len()) }.copy_from_slice(&header);
            #[allow(clippy::cast_possible_truncation)]
            let run_start = h_off as u64;
            for r in &recs {
                let bs = u32::try_from(r.len()).unwrap();
                let po = unsafe { arena.grow_uninit(4) };
                unsafe { arena.slice_mut(po, 4) }.copy_from_slice(&bs.to_le_bytes());
                let bo = unsafe { arena.grow_uninit(r.len()) };
                unsafe { arena.slice_mut(bo, r.len()) }.copy_from_slice(r);
            }
            #[allow(clippy::cast_possible_truncation)]
            let run_len = arena.len() as u64 - run_start;
            let arena = Arc::new(PooledSegmentedBuf::unpooled(arena));

            let mut step = FindBoundariesAndSort::new(
                QuerynameStrategy::<RawQuerynameLexKey>::new(MemoryChunkErased::QuerynameLex),
                sort_threads,
                64 * 1024 * 1024,
            );
            assert_eq!(
                step.strategy().sort_pool_threads(),
                None,
                "sort pool is not built until the first seal"
            );
            step.ingest_block(&InflatedBlock {
                arena: Arc::clone(&arena),
                ordinal: 0,
                offset: run_start,
                len: u32::try_from(run_len).unwrap(),
                is_last_of_run: false,
                run_seq: 0,
                seals_to_spill: false,
            })
            .expect("ingest_block must succeed");
            // finalize() seals the residual run, which builds + installs the
            // bounded sort pool sized to `sort_threads`.
            let _ = step.finalize().expect("finalize must succeed");
            assert_eq!(
                step.strategy().sort_pool_threads(),
                Some(sort_threads),
                "Phase-1 queryname sort pool must be sized to sort_threads={sort_threads} \
                 (the value FindBoundariesAndSort was constructed with)"
            );
        }
    }
}

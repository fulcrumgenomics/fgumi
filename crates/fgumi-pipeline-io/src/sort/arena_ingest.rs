//! `ReadBlocks` serial arena-admit step, `InflateToArena` parallel step, and
//! the per-order [`ArenaSortStrategy`]s the record-input `SortBuffer` uses.
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
//! disjoint arena slot via [`fgumi_bgzf::decompress_into_slice_with_crc`] under
//! the step's CRC policy.  This is the first
//! `unsafe` site in `fgumi-pipeline-io`; see CLAUDE.md §"Approved hot-path
//! unsafe (parallel-inflate sort ingest, fgumi-pipeline-io)" for the
//! full justification and SAFETY invariant.
//!
//! Together `ReadBlocks` and `InflateToArena` form the head of the BAM sort's
//! block-input front. The chain builder's `add_sort`
//! (`src/lib/pipeline/chains/builder.rs`) follows them with the order's split
//! front in [`arena_split`](super::arena_split) (boundaries → key extraction →
//! run sort), feeding the Phase-1 spill/merge tail, when the sort source is raw
//! BGZF blocks.

use std::collections::VecDeque;
use std::io;
use std::sync::Arc;
use std::time::Instant;

use fgumi_bam_io::pread::RawFrame;
use fgumi_bgzf::{Decompressor, decompress_into_slice_with_crc};
use fgumi_sort::{
    ArenaPool, BoundedSortPool, PooledSegmentedBuf, QN_CLASS_BATCH_RECORDS, QnBatchLcp,
    QnComparator, QnRef, QnRefBuilder, QnSortState, RecordRef, coordinate_chunk_from_refs,
    extract_coordinate_key_inline, queryname_chunk_from_qn_refs,
};

use crate::types::BgzfBlock;

use fgumi_pipeline_core::Unpushed;
use fgumi_pipeline_core::held::HeldSlot;
use fgumi_pipeline_core::item::{HeapSize, Ordered};
use fgumi_pipeline_core::outputs::OrderedBytesSingle;
use fgumi_pipeline_core::queues::QueueSpec;
use fgumi_pipeline_core::reorder::BranchOrdering;
use fgumi_pipeline_core::step::{
    CounterSpec, DetachedGroup, Step, StepCtx, StepKind, StepOutcome, StepProfile,
};
use fgumi_pipeline_core::{PhaseCap, WholeCap, admit_input};

use crate::sort::arena_split::QnChunkOrder;
use crate::sort::protocol::MemoryChunkErased;

// ============================================================================
// Constants
// ============================================================================

/// Reserved front region at the start of every run's arena, in bytes.
///
/// Block 0 always lands at offset `FRONT_REGION`.  The region `[0, FRONT_REGION)`
/// is left uninitialized by `ReadBlocks` and filled by the split front's boundary
/// step (`FindRecordBoundaries`) when it writes a carried straddler from the previous run (Task 2, 3b.4).
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

/// `ReadBlocks` counter slot index: BGZF blocks consumed from `ctx.input` this call.
const READ_BLOCKS_BLOCKS: usize = 0;
/// `ReadBlocks` counter slot index: compressed bytes consumed from `ctx.input` this call.
const READ_BLOCKS_BYTES_READ: usize = 1;

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
/// `is_last_of_run`/`seals_to_spill`; downstream (`FindRecordBoundaries`) reads
/// those fields only on the `is_last_of_run` block, so non-last blocks emit
/// `is_last_of_run = false, seals_to_spill = false`.
///
/// With pool capacity 1, the next run's arena cannot be acquired until the
/// prior run's [`Arc<PooledSegmentedBuf>`] drops (i.e. after the spill steps
/// and the merge release its chunk).  If acquisition fails, `try_run` returns
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
    /// When the current wait for a free arena began (its first refusal);
    /// `None` while not waiting.
    wait_started: Option<Instant>,
    /// Waits for a free arena that ended with an acquire.
    arena_waits: u64,
    /// Total time from each wait's first refusal to its acquire.
    arena_wait_ns: u64,
    /// Total time from the release that ended each wait to its acquire — the
    /// wake latency, counted only for a release after the wait began.
    arena_wake_ns: u64,
    /// The level the arena-wait line is still to be logged at — `info` under
    /// `--sort-stats`, `debug` otherwise; `None` once logged.
    arena_line_level: Option<log::Level>,
    /// Test hook run between the refused acquire and the waiter registration.
    #[cfg(test)]
    before_register: Option<Box<dyn FnOnce() + Send>>,
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
    /// keeps the boundary step's contiguous-slice scan valid WITHOUT
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
            wait_started: None,
            arena_waits: 0,
            arena_wait_ns: 0,
            arena_wake_ns: 0,
            arena_line_level: Some(log::Level::Debug),
            #[cfg(test)]
            before_register: None,
        }
    }

    /// Log the arena-wait line at `info` (`--sort-stats`) instead of `debug`.
    #[must_use]
    pub fn with_sort_stats(mut self, on: bool) -> Self {
        self.arena_line_level = Some(if on { log::Level::Info } else { log::Level::Debug });
        self
    }

    /// Acquire the next run's arena, registering as the pool's waiter when it
    /// is exhausted so the release — on whichever thread drops the previous
    /// run's last reference — unparks this step's driver. Accounts each wait
    /// and its wake latency.
    fn try_acquire_arena(&mut self) -> Option<PooledSegmentedBuf> {
        let acquired = self.pool.try_acquire().or_else(|| {
            #[cfg(test)]
            if let Some(hook) = self.before_register.take() {
                hook();
            }
            // Register, then retry: a release after the registration unparks
            // us; one between the refused acquire and the registration is
            // caught by the retry.
            self.pool.register_waiter();
            self.pool.try_acquire()
        });
        match acquired {
            None => {
                self.wait_started.get_or_insert_with(Instant::now);
                None
            }
            Some(arena) => {
                if let Some(started) = self.wait_started.take() {
                    let now = Instant::now();
                    self.arena_waits += 1;
                    self.arena_wait_ns += nanos(now - started);
                    if let Some(released) = self.pool.released_at().filter(|&r| r >= started) {
                        self.arena_wake_ns += nanos(now.saturating_duration_since(released));
                    }
                }
                Some(arena)
            }
        }
    }

    /// The arena-wait summary line.
    fn arena_wait_line(&self) -> String {
        #[allow(clippy::cast_precision_loss)] // display seconds
        let wait_s = self.arena_wait_ns as f64 / 1e9;
        format!(
            "Phase 1: arena waits {} ({wait_s:.1} s waiting for a free arena; {} ms from return \
             to reuse)",
            self.arena_waits,
            self.arena_wake_ns / 1_000_000
        )
    }

    /// Log [`Self::arena_wait_line`] once.
    fn log_arena_waits(&mut self) {
        if let Some(level) = self.arena_line_level.take() {
            log::log!(level, "{}", self.arena_wait_line());
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
        let Some(mut arena) = self.try_acquire_arena() else {
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
        // inflate worker, the `[0, FRONT_REGION)` carry region by the boundary step
        // (`FindRecordBoundaries`) — before any read.  Bytes never read (the unused
        // tail past a run's end, and the front region when there is no straddler)
        // are never observed: the boundary step only ever slices
        // `[scan_start, run_end)`.
        #[allow(unsafe_code)]
        let _all = unsafe { arena.grow_uninit(self.segment_size) };
        // `ArenaPool::try_acquire` already returns the `PooledSegmentedBuf` wrapper
        // (it wraps internally so a caller cannot orphan a pooled arena — the hang
        // that motivated making `PooledSegmentedBuf::pooled` private), so the arena
        // is wrapped by construction here and only needs sharing.
        self.arena = Some(Arc::new(arena));
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

    /// Test seam: acquire (or wait for) the run's arena.
    #[cfg(test)]
    pub(crate) fn ensure_arena_for_test(&mut self) -> bool {
        self.ensure_arena()
    }

    /// Test seam: hand the current run's arena to the caller (as if the run
    /// had sealed and its chunk were still downstream).
    #[cfg(test)]
    pub(crate) fn take_arena_for_test(&mut self) -> Option<Arc<PooledSegmentedBuf>> {
        self.arena.take()
    }

    /// Test seam: run `hook` between the next refused acquire and the waiter
    /// registration.
    #[cfg(test)]
    pub(crate) fn before_register_for_test(&mut self, hook: Box<dyn FnOnce() + Send>) {
        self.before_register = Some(hook);
    }

    /// Test seam: `(arena_waits, arena_wait_ns, arena_wake_ns)`.
    #[cfg(test)]
    pub(crate) fn arena_wait_stats(&self) -> (u64, u64, u64) {
        (self.arena_waits, self.arena_wait_ns, self.arena_wake_ns)
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

/// A duration in whole nanoseconds, saturating.
fn nanos(d: std::time::Duration) -> u64 {
    u64::try_from(d.as_nanos()).unwrap_or(u64::MAX)
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

    fn counters(&self) -> &'static [CounterSpec] {
        const SPECS: &[CounterSpec] =
            &[CounterSpec::new("blocks", "blocks"), CounterSpec::new("bytes_read", "bytes")];
        SPECS
    }

    // `try_run` bumps the `blocks`/`bytes_read` counters exactly ONCE per call
    // with the batch total (never per-record), after delegating the actual
    // (long, multi-early-return) state machine to `try_run_inner` below.
    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        let mut admitted_blocks: u64 = 0;
        let mut admitted_bytes: u64 = 0;
        let outcome = self.try_run_inner(ctx, &mut admitted_blocks, &mut admitted_bytes);
        ctx.counters.add(READ_BLOCKS_BLOCKS, admitted_blocks);
        ctx.counters.add(READ_BLOCKS_BYTES_READ, admitted_bytes);
        if matches!(outcome, Ok(StepOutcome::Finished)) {
            self.log_arena_waits();
        }
        outcome
    }

    fn new_worker_copy(&self) -> Self {
        // Serial steps are never cloned by the framework; this is unreachable.
        panic!("ReadBlocks is Serial — new_worker_copy should never be called")
    }
}

impl ReadBlocks {
    /// Body of [`Step::try_run`], factored out so the counter bump in the trait
    /// method stays a single pair of calls regardless of which early-return path
    /// below fires. `admitted_blocks`/`admitted_bytes` accumulate the batch total
    /// for this call as blocks are popped from `ctx.input` (see step 5); the
    /// caller bumps `ctx.counters` with the final tally after this returns.
    ///
    /// A single cohesive deferred-seal state machine: held-output retry, seal
    /// arming/execution, batched admit, and staged-emit drain are tightly coupled
    /// by `seal_pending` / `deferred_block` / `held` and are clearer read top to
    /// bottom than split across helpers that would each need the same state —
    /// hence the line count.
    ///
    /// Returns `io::Result` (rather than a bare `StepOutcome`) to stay
    /// signature-symmetric with the `Step::try_run` it was extracted from, even
    /// though no path here currently returns `Err` — a future error path (e.g. a
    /// new `?`-fallible admit step) should be free to add one without having to
    /// first change this function's return type.
    #[allow(clippy::too_many_lines, clippy::unnecessary_wraps)]
    fn try_run_inner(
        &mut self,
        ctx: &mut StepCtx<'_, Self>,
        admitted_blocks: &mut u64,
        admitted_bytes: &mut u64,
    ) -> io::Result<StepOutcome> {
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
            // Count the block as consumed the moment it is popped — regardless
            // of whether it is admitted immediately, deferred to a seal, or
            // handed back on pool exhaustion below — since input consumption,
            // not admission outcome, is what the counter tracks.
            *admitted_blocks += 1;
            *admitted_bytes += block.bytes.len() as u64;
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
                    // The pop took an input item, so this is `Progress` (the
                    // `StepOutcome::Progress` contract): it freed a slot on the
                    // input edge, and only a `Progress` delivers the reverse wake
                    // to an upstream producer holding an item for that slot. The
                    // next call finds the deferred block and reports `NoProgress`
                    // until the arena returns, without popping.
                    self.deferred_block = Some(b);
                    return Ok(StepOutcome::Progress);
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
    /// Complete raw BGZF block bytes (header + deflate payload + footer),
    /// borrowed from its read slice on the native input path.
    pub block: RawFrame,
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
    /// The frame's charge ([`crate::types::raw_frame_heap_size`]): an owned
    /// (carried) frame's allocation, a borrowed frame's range.
    fn heap_size(&self) -> usize {
        crate::types::raw_frame_heap_size(&self.block)
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
    /// Whether each block's BGZF footer CRC32 is compared after inflate. The
    /// ISIZE / exact-fill checks always run; see
    /// [`fgumi_bgzf::decompress_into_slice_with_crc`].
    verify_crc: bool,
    /// Admission cap shared with the other phase-1 pool steps (set by the
    /// chain builder via `with_phase_cap`; `None` = uncapped).
    cap: Option<Arc<PhaseCap>>,
}

impl InflateToArena {
    /// Create a new `InflateToArena` step that verifies every block's CRC32.
    ///
    /// Equivalent to [`Self::new_with_crc`]`(output_byte_limit, true)`.
    #[must_use]
    pub fn new(output_byte_limit: u64) -> Self {
        Self::new_with_crc(output_byte_limit, true)
    }

    /// Create a new `InflateToArena` step with an explicit CRC32 policy.
    ///
    /// `output_byte_limit` bounds the byte-counted output queue (tokens carry
    /// `heap_size == 0`, so this mainly controls queue depth via the framework's
    /// backpressure mechanism). `verify_crc = false` skips only the footer
    /// CRC32 compare, for trusted input; the ISIZE / exact-fill checks always
    /// run.
    #[must_use]
    pub fn new_with_crc(output_byte_limit: u64, verify_crc: bool) -> Self {
        Self {
            decompressor: Decompressor::new(),
            held: HeldSlot::new(),
            output_byte_limit,
            verify_crc,
            cap: None,
        }
    }

    /// Share a phase admission cap with the other pool steps of this phase. At
    /// most `cap.max()` workers are inside any of those steps at once. `None`
    /// leaves the step uncapped.
    #[must_use]
    pub fn with_phase_cap(mut self, cap: Option<Arc<PhaseCap>>) -> Self {
        self.cap = cap;
        self
    }

    /// One `try_run` with the phase cap moved out of `self` (see `try_run`).
    fn run_once(
        &mut self,
        ctx: &mut StepCtx<'_, Self>,
        cap: Option<&PhaseCap>,
    ) -> io::Result<StepOutcome> {
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

        // Admission (see the pipeline-core `admission` module): a capped step
        // takes its phase permit BEFORE popping (an empty poll takes none and
        // is recorded as an empty pop), so a refused clone leaves the item for
        // a clone that may run and reports `Capped`. Uncapped, this is a no-op
        // and the pop below runs exactly as without a cap. The permit is held
        // to the end of this call.
        let _permit = match admit_input(ctx.input, cap) {
            Ok(permit) => permit,
            Err(outcome) => return Ok(outcome),
        };

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
        // required contract — `decompress_into_slice_with_crc` below fills every
        // byte (the exact-fill guarantee is unconditional; `verify_crc` gates
        // only the CRC32 compare).
        // `offset` is a u64 byte offset into the arena; on a 64-bit platform
        // this always fits in usize — the arena itself cannot exceed
        // `isize::MAX` bytes, which is the Rust allocation bound.
        let offset_usize =
            usize::try_from(offset).expect("arena offset must fit in usize on this platform");
        let len_usize = len as usize; // u32 always fits in usize

        #[allow(unsafe_code)]
        let slot = unsafe { arena.slice_mut(offset_usize, len_usize) };

        let n =
            decompress_into_slice_with_crc(&block, &mut self.decompressor, slot, self.verify_crc)?;
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
        // Moved out for the call: the permit borrows the cap while `inflate_one`
        // needs `&mut self`. A move, so no reference-count traffic.
        let cap = self.cap.take();
        let outcome = self.run_once(ctx, cap.as_deref());
        self.cap = cap;
        outcome
    }

    fn new_worker_copy(&self) -> Self {
        Self {
            cap: self.cap.clone(),
            ..Self::new_with_crc(self.output_byte_limit, self.verify_crc)
        }
    }

    fn phase_cap(&self) -> Option<&PhaseCap> {
        self.cap.as_deref()
    }
}

// ============================================================================
// Arena sort strategy (per-order key extraction + seal)
// ============================================================================

/// Per-sort-order policy for the record-input front (`SortBuffer`), which
/// pushes records one at a time.
///
/// The strategy owns the accumulated per-record refs and knows how to (a) extract
/// a record's sort key from its arena-resident body and (b) at run seal, sort
/// those refs and wrap the shared arena into the correctly-typed
/// [`MemoryChunkErased`] — with zero record-body copies. `SortBuffer` is generic
/// over this trait and monomorphised per order, so the per-record
/// [`push_record`](Self::push_record) call inlines with no dynamic dispatch. The
/// block-input front extracts in batches instead
/// ([`ArenaKeyExtractor`](super::arena_split::ArenaKeyExtractor)).
pub trait ArenaSortStrategy: Send + 'static {
    /// Extract the sort key from `body` (the record's BAM body, `block_size`
    /// prefix excluded) and accumulate a ref pointing at `(body_off, len)` in the
    /// shared inflate arena. Called once per record on the boundary scan hot path.
    ///
    /// # Errors
    ///
    /// Returns an error if the record is invalid for this order (e.g. a
    /// template-coordinate dropped-lane violation). Coordinate never errors.
    fn push_record(&mut self, body: &[u8], body_off: u64, len: u32) -> io::Result<()>;

    /// Sort the refs accumulated for this run on `sort_threads` threads and wrap
    /// `arena` into the erased chunk (zero record copies), resetting the
    /// accumulator for the next run.
    fn seal(&mut self, arena: Arc<PooledSegmentedBuf>, sort_threads: usize) -> MemoryChunkErased;

    /// Whether sealing the refs accumulated so far would sort on more than one
    /// thread at `sort_threads` — what decides whether the seal must take the
    /// whole phase cap (see `seal_under_cap`).
    fn sorts_in_parallel(&self, sort_threads: usize) -> bool;
}

/// Reserve the whole phase-1 cap for an off-pool run sort: the one reservation
/// path of every run seal ([`seal_under_cap`] and the split front's `SortRun`).
///
/// Returns `Ok(None)` when there is no cap (an unset `--sort-threads`) or the
/// sort runs on one thread — a serial sort is the coordination thread's own
/// work, which the cap does not count. Otherwise registers a pending
/// reservation (new pool admissions are refused while it waits), parks until
/// every permit is held, and returns the guard; the sort runs while it is held.
///
/// # Errors
///
/// Returns [`io::ErrorKind::Interrupted`] when the run is cancelled (or failed)
/// while waiting; the caller must not sort.
pub(crate) fn reserve_for_run_sort(
    cap: Option<&PhaseCap>,
    sorts_in_parallel: bool,
) -> io::Result<Option<WholeCap<'_>>> {
    match cap {
        Some(cap) if sorts_in_parallel => cap.acquire_whole().map(Some).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::Interrupted,
                "sort cancelled while waiting for its phase cap",
            )
        }),
        _ => Ok(None),
    }
}

/// Seal one run under the phase-1 cap.
///
/// Uncapped: the sort runs at its full `sort_threads` width. Capped (an
/// explicit `--sort-threads`): when the run sorts in parallel, the per-run
/// sort is phase-1 work on its own rayon threads beside the pool steps that
/// share the cap, so it takes the whole cap ([`reserve_for_run_sort`]) and
/// sorts on exactly the cap's width — pool admission pauses until it is done,
/// and sort threads + inflate + spill-compression workers never exceed the cap
/// together. A run that sorts on one thread takes no permit.
///
/// # Errors
///
/// Returns [`io::ErrorKind::Interrupted`] without sorting when the run is
/// cancelled (or failed) while the seal waits for the cap.
pub(crate) fn seal_under_cap<S: ArenaSortStrategy>(
    strategy: &mut S,
    cap: Option<&PhaseCap>,
    arena: Arc<PooledSegmentedBuf>,
    sort_threads: usize,
) -> io::Result<MemoryChunkErased> {
    let sort_threads = sort_threads.max(1);
    let whole = reserve_for_run_sort(cap, strategy.sorts_in_parallel(sort_threads))?;
    let width = whole.as_ref().map_or(sort_threads, WholeCap::width);
    Ok(strategy.seal(arena, width))
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
    fn push_record(&mut self, body: &[u8], body_off: u64, len: u32) -> io::Result<()> {
        let sort_key = extract_coordinate_key_inline(body, self.n_ref);
        self.refs.push(RecordRef::new(sort_key, body_off, len));
        Ok(())
    }

    fn seal(&mut self, arena: Arc<PooledSegmentedBuf>, sort_threads: usize) -> MemoryChunkErased {
        let refs = std::mem::take(&mut self.refs);
        MemoryChunkErased::Coordinate(coordinate_chunk_from_refs(arena, refs, sort_threads))
    }

    fn sorts_in_parallel(&self, sort_threads: usize) -> bool {
        fgumi_sort::radix_sorts_in_parallel(self.refs.len(), sort_threads)
    }
}

/// Template-coordinate strategy: wraps a [`fgumi_sort::TemplateArenaAccumulator`],
/// which owns the library / cell-barcode / MI and `--key-types` narrowed-lane
/// machinery the template key needs; accumulates arena-pointing refs and seals to
/// [`MemoryChunkErased::TemplateCoordinate`] (an arena-backed
/// `InMemoryChunk<TemplateKey>`), byte-identical to the owned `TemplateChunkSorter`.
///
/// Its production user is `SortBuffer`, the record-input front; the block-input
/// arena front sorts template order through [`super::arena_split`].
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
    fn push_record(&mut self, body: &[u8], body_off: u64, len: u32) -> io::Result<()> {
        self.acc
            .push(body, body_off, len)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("{e:#}")))
    }

    fn seal(&mut self, arena: Arc<PooledSegmentedBuf>, sort_threads: usize) -> MemoryChunkErased {
        MemoryChunkErased::TemplateCoordinate(self.acc.seal(arena, sort_threads))
    }

    fn sorts_in_parallel(&self, sort_threads: usize) -> bool {
        fgumi_sort::radix_sorts_in_parallel(self.acc.pending_len(), sort_threads)
    }
}

/// Queryname-order strategy: builds one 40-byte [`QnRef`] per record — a
/// window of the read name behind a learned common prefix, the flag class and
/// the record's arena position — and at seal sorts the refs with
/// [`queryname_chunk_from_qn_refs`] into an arena-backed `InMemoryChunk<QnKey>`,
/// wrapped into the matching erased variant. Nothing is copied or allocated
/// per record.
///
/// The arena front's queryname orders run on the split front's
/// [`QuerynameKeys`](super::arena_split::QuerynameKeys), which extracts through
/// the same [`QnRefBuilder`]; this strategy is the record-input path.
///
/// The prefix is learned per batch of [`QN_CLASS_BATCH_RECORDS`] records. This
/// front pushes records one at a time with no arena to read back until the
/// seal, so a batch whose guessed prefix failed is recorded at the batch's
/// close and re-windowed at seal, before the sort.
///
/// Equal keys keep input order (ties break by arena offset), so one sorted
/// chunk per run is the stable sort of the run.
pub struct QuerynameStrategy<C: QnComparator> {
    /// Refs accumulated across the current run.
    refs: Vec<QnRef>,
    /// The sort's shared class table and counters.
    keys: QnSortState,
    /// Carries the class guess across batches.
    builder: QnRefBuilder,
    /// The common prefix of the open batch's names.
    lcp: QnBatchLcp,
    /// Index in `refs` where the open batch starts.
    batch_start: usize,
    /// Closed batches to re-window at seal: `(refs range, class)`.
    pending: Vec<(std::ops::Range<usize>, u8)>,
    /// Bounded rayon pool (sized to `sort_threads` on the first seal) the run's
    /// comparator sort runs on, so it does not oversubscribe the pipeline's
    /// workers.
    sort_pool: BoundedSortPool,
    /// The comparator, which also picks the erased chunk variant
    /// ([`QnChunkOrder::wrap`]).
    _order: std::marker::PhantomData<C>,
}

impl<C: QnChunkOrder> QuerynameStrategy<C> {
    /// A queryname strategy over the sort's shared key state `keys`; `C` picks
    /// the sealed chunk's erased variant.
    #[must_use]
    pub fn new(keys: QnSortState) -> Self {
        Self {
            refs: Vec::new(),
            keys,
            builder: QnRefBuilder::default(),
            lcp: QnBatchLcp::default(),
            batch_start: 0,
            pending: Vec::new(),
            sort_pool: BoundedSortPool::new("qname-sort"),
            _order: std::marker::PhantomData,
        }
    }

    /// Replace the key state (before any record is pushed): the chain builder
    /// shares one state across a sort's front and merge.
    pub fn set_keys(&mut self, keys: QnSortState) {
        debug_assert!(self.refs.is_empty(), "the key state is set before ingest");
        self.keys = keys;
    }

    /// Test-only: the run's refs as the seal would sort them — the open batch
    /// closed and every adopted batch re-windowed.
    #[cfg(test)]
    pub(crate) fn refs_for_test(&mut self, arena: &PooledSegmentedBuf) -> Vec<QnRef> {
        self.close_batch();
        for (range, class) in self.pending.drain(..) {
            QnRefBuilder::rewindow(&self.keys.table, arena, &mut self.refs[range], class);
        }
        self.refs.clone()
    }

    /// Close the open batch: adopt a learned class if its guess failed (recorded
    /// for re-windowing at seal) and count it.
    fn close_batch(&mut self) {
        let batch = &self.refs[self.batch_start..];
        if batch.is_empty() {
            return;
        }
        let adopted = self.builder.adopt_batch_class(&self.keys.table, batch, &self.lcp);
        let class0 =
            if adopted.is_some() { 0 } else { batch.iter().filter(|r| r.0.class == 0).count() };
        self.keys.stats.note_batch(batch.len(), class0, adopted.is_some());
        if let Some(class) = adopted {
            self.pending.push((self.batch_start..self.refs.len(), class));
        }
        self.batch_start = self.refs.len();
        self.lcp = QnBatchLcp::default();
    }
}

impl<C: QnChunkOrder> ArenaSortStrategy for QuerynameStrategy<C> {
    #[inline]
    fn push_record(&mut self, body: &[u8], body_off: u64, len: u32) -> io::Result<()> {
        let r = self.builder.build(&self.keys.table, &mut self.lcp, body, body_off, len);
        self.refs.push(r);
        if self.refs.len() - self.batch_start == QN_CLASS_BATCH_RECORDS {
            self.close_batch();
        }
        Ok(())
    }

    fn seal(&mut self, arena: Arc<PooledSegmentedBuf>, sort_threads: usize) -> MemoryChunkErased {
        self.close_batch();
        for (range, class) in self.pending.drain(..) {
            QnRefBuilder::rewindow(&self.keys.table, &arena, &mut self.refs[range], class);
        }
        let refs = std::mem::take(&mut self.refs);
        self.batch_start = 0;
        let keys = &self.keys;
        self.sort_pool.bounded(sort_threads).install(move || {
            C::wrap(queryname_chunk_from_qn_refs::<C>(
                &keys.table,
                arena,
                refs,
                keys.compare_stats(),
            ))
        })
    }

    /// The comparator sort always runs inside the per-sort rayon pool, so it is
    /// parallel whenever that pool has more than one thread.
    fn sorts_in_parallel(&self, sort_threads: usize) -> bool {
        sort_threads > 1
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use fgumi_sort::{PooledSegmentedBuf, SegmentedBuf};
    use rstest::rstest;
    use std::sync::Arc;
    use std::sync::atomic::{AtomicBool, Ordering};

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

    /// Acquire a pooled arena for constructing block tokens in tests.
    fn test_arena() -> Arc<PooledSegmentedBuf> {
        let pool = ArenaPool::new(1, 1024);
        Arc::new(pool.try_acquire().expect("a fresh pool always has one arena"))
    }

    // -----------------------------------------------------------------------
    // Block token accounting (HeapSize / Ordered)
    // -----------------------------------------------------------------------

    #[test]
    fn arena_block_charges_only_its_compressed_bytes_to_the_heap() {
        // The arena is shared via `Arc`, so it must NOT be counted per block —
        // byte-bounded queues would otherwise over-charge by the arena size and
        // stall the pipeline.
        let block = ArenaBlock {
            arena: test_arena(),
            ordinal: 7,
            offset: 64,
            len: 128,
            block: vec![0xAB; 40].into(),
            is_last_of_run: true,
            run_seq: 2,
            seals_to_spill: false,
        };
        assert_eq!(block.heap_size(), 40, "only the owned BGZF bytes are charged");
        assert_eq!(block.ordinal(), 7, "ordinal drives ByItemOrdinal reordering");
    }

    /// An owned (carried) frame is charged its allocation, like `BgzfBlock`'s.
    #[test]
    fn arena_block_charges_an_owned_frame_its_capacity() {
        let mut bytes = Vec::with_capacity(64);
        bytes.extend_from_slice(&[0xAB; 40]);
        let block = ArenaBlock {
            arena: test_arena(),
            ordinal: 0,
            offset: 0,
            len: 128,
            block: bytes.into(),
            is_last_of_run: false,
            run_seq: 0,
            seals_to_spill: false,
        };
        assert_eq!(block.heap_size(), 64);
    }

    #[test]
    fn inflated_block_is_heap_free_because_its_bytes_live_in_the_arena() {
        let token = InflatedBlock {
            arena: test_arena(),
            ordinal: 11,
            offset: 0,
            len: 256,
            is_last_of_run: false,
            run_seq: 1,
            seals_to_spill: true,
        };
        assert_eq!(token.heap_size(), 0, "the token owns no bytes of its own");
        assert_eq!(token.ordinal(), 11);
    }

    // -----------------------------------------------------------------------
    // Step wiring
    // -----------------------------------------------------------------------

    #[test]
    fn read_blocks_runs_off_pool_on_the_coordination_driver() {
        let step = ReadBlocks::new(1 << 20, 4096);
        let profile = step.profile();
        assert_eq!(profile.name, "ReadBlocks");
        // Detached keeps the serial arena-admit off a pool worker so it cannot
        // steal a slot from the parallel inflaters.
        assert_eq!(profile.kind, StepKind::Detached);
        assert!(!profile.sticky);
        assert_eq!(profile.branch_ordering, vec![BranchOrdering::ByItemOrdinal]);
        match profile.output_queues.as_slice() {
            [QueueSpec::ByteBounded { limit_bytes }] => assert_eq!(*limit_bytes, 4096),
            other => panic!("expected one byte-bounded queue, got {other:?}"),
        }
        assert_eq!(step.detached_group(), DetachedGroup::Shared(crate::sort::SORT_COORD_GROUP));
    }

    #[test]
    fn inflate_to_arena_is_parallel_and_worker_copies_share_only_config() {
        let step = InflateToArena::new(8192);
        let profile = step.profile();
        assert_eq!(profile.name, "InflateToArena");
        assert_eq!(profile.kind, StepKind::Parallel, "inflation fans out across workers");
        assert!(!profile.sticky);
        assert_eq!(profile.branch_ordering, vec![BranchOrdering::ByItemOrdinal]);

        // Each worker gets an independent copy carrying the same queue bound.
        let worker = step.new_worker_copy();
        match worker.profile().output_queues.as_slice() {
            [QueueSpec::ByteBounded { limit_bytes }] => assert_eq!(*limit_bytes, 8192),
            other => panic!("expected one byte-bounded queue, got {other:?}"),
        }
    }

    // -----------------------------------------------------------------------
    // Sort strategies
    // -----------------------------------------------------------------------

    #[test]
    fn coordinate_strategy_accumulates_refs_and_seal_drains_them() {
        let mut strategy = CoordinateStrategy::new(4);

        // Two synthetic record bodies; only the leading tid/pos fields are read
        // by the inline key extractor, so a minimal fixed-size body suffices.
        let body = vec![0u8; 36];
        strategy.push_record(&body, 0, 36).unwrap();
        strategy.push_record(&body, 36, 36).unwrap();
        assert_eq!(strategy.refs.len(), 2, "each pushed record yields one arena ref");

        let chunk = strategy.seal(test_arena(), 1);
        assert!(matches!(chunk, MemoryChunkErased::Coordinate(_)), "coordinate seal");
        assert!(strategy.refs.is_empty(), "seal must drain the accumulator for the next run");
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

    /// A dispatch that pops a block reports `Progress` even when the arena pool
    /// is exhausted and the block is deferred: the pop freed a slot on the
    /// input edge, and the plan's reverse wake to an upstream producer holding
    /// an item for it runs only on `Progress`. The next dispatch pops nothing
    /// and is idle until the arena returns, then admits the deferred block.
    #[test]
    fn a_popped_block_deferred_on_an_exhausted_pool_is_progress() {
        use fgumi_pipeline_core::testing::StepProbe;
        let payload = b"deferred-block-payload";
        let bytes = make_test_bgzf_block(payload);
        let isz = u32::try_from(uncompressed_size_of(&bytes)).unwrap();
        let mut step = ReadBlocks::new(64 * 1024 * 1024, 64 * 1024 * 1024);
        let probe = StepProbe::new(&step);
        let stolen = step.pool.try_acquire().expect("the pool's one arena");
        probe.push_input(BgzfBlock {
            batch_serial: 0,
            bytes: bytes.into(),
            uncompressed_size: isz,
            index: None,
        });
        assert_eq!(probe.try_run(&mut step).expect("try_run"), StepOutcome::Progress);
        assert!(probe.input_is_empty(), "the block was popped");
        assert!(step.deferred_block.is_some(), "...and deferred, not dropped");
        assert_eq!(
            probe.try_run(&mut step).expect("try_run"),
            StepOutcome::NoProgress,
            "the arena is still out: idle, and nothing popped"
        );
        drop(stolen);
        assert_eq!(probe.try_run(&mut step).expect("try_run"), StepOutcome::Progress);
        assert!(step.deferred_block.is_none(), "the arena returned: the block is admitted");
    }

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
                bytes: blk0.clone().into(),
                uncompressed_size: u32::try_from(isz0).unwrap(),
                index: None,
            })
            .is_ok()
        );
        assert!(
            step.admit_block(BgzfBlock {
                batch_serial: 1,
                bytes: blk1.clone().into(),
                uncompressed_size: u32::try_from(isz1).unwrap(),
                index: None,
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
                bytes: blk0.clone().into(),
                uncompressed_size: u32::try_from(isz0).unwrap(),
                index: None,
            })
            .is_ok()
        );
        assert!(
            step.admit_block(BgzfBlock {
                batch_serial: 1,
                bytes: blk1.clone().into(),
                uncompressed_size: u32::try_from(isz1).unwrap(),
                index: None,
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
                bytes: probe_blk.into(),
                uncompressed_size: u32::try_from(probe_isz).unwrap(),
                index: None,
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
                bytes: blk2.clone().into(),
                uncompressed_size: u32::try_from(isz2).unwrap(),
                index: None,
            })
            .is_ok(),
            "run 1 admit must succeed after pool release"
        );
        assert!(
            step.admit_block(BgzfBlock {
                batch_serial: 3,
                bytes: blk3.clone().into(),
                uncompressed_size: u32::try_from(isz3).unwrap(),
                index: None,
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
            block: block.into(),
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

    /// Lost-wake guard: a release between the refused acquire and the waiter
    /// registration is caught by the retry, deterministically — the step gets
    /// its arena without ever waiting.
    #[test]
    fn register_then_retry_catches_a_release_in_between() {
        let mut step = ReadBlocks::new(1 << 16, 1 << 20);
        assert!(step.ensure_arena_for_test(), "the pool's one arena");
        let held = step.take_arena_for_test().expect("arena");
        step.before_register_for_test(Box::new(move || drop(held)));
        assert!(
            step.ensure_arena_for_test(),
            "the retry after registering finds the released arena"
        );
        assert_eq!(step.arena_wait_stats().0, 0, "no wait was recorded");
    }

    /// Wait vs wake latency: the wait spans first refusal → acquire; the wake
    /// latency spans release → acquire. Refuse, wait ≥ 100 ms, release, wait
    /// ≥ 30 ms, then acquire.
    #[test]
    fn read_blocks_separates_arena_wait_from_wake_latency() {
        use std::time::Duration;
        let mut step = ReadBlocks::new(1 << 16, 1 << 20);
        assert!(step.ensure_arena_for_test());
        let held = step.take_arena_for_test().expect("arena");
        assert!(!step.ensure_arena_for_test(), "refused while the arena is held");
        std::thread::sleep(Duration::from_millis(100));
        drop(held);
        std::thread::sleep(Duration::from_millis(30));
        assert!(step.ensure_arena_for_test(), "acquired after the release");
        let (waits, wait_ns, wake_ns) = step.arena_wait_stats();
        assert_eq!(waits, 1);
        assert!(wait_ns >= 130_000_000, "wait {wait_ns} ns");
        assert!(wake_ns >= 30_000_000, "wake {wake_ns} ns");
        assert!(wake_ns <= wait_ns - 100_000_000, "wake {wake_ns} ns within wait {wait_ns} ns");
    }

    /// The record-input front cannot re-window a batch when it closes (it has
    /// no arena yet): a failed guess is recorded and applied at seal. The first
    /// batch is built in class 0, learns `A00132:53:HFHJ:1:`, and comes out of
    /// the seal in that class; the next batch is built in it directly.
    #[test]
    fn queryname_strategy_rewindows_pending_batches_at_seal() {
        use fgumi_raw_bam::testutil::make_bam_bytes;
        use fgumi_sort::{QN_CLASS_BATCH_RECORDS, QnLex, RawQuerynameLexKey, RawSortKey};

        let n = QN_CLASS_BATCH_RECORDS + 100;
        let bodies: Vec<Vec<u8>> = (0..n)
            .map(|i| {
                let name = format!("A00132:53:HFHJ:1:{}:{}", 1101 + (i * 7919) % 1000, i);
                make_bam_bytes(0, 0, 0x41, name.as_bytes(), &[], 4, -1, -1, &[])
            })
            .collect();
        let mut arena = SegmentedBuf::new();
        let keys = QnSortState::new(false);
        let mut strategy = QuerynameStrategy::<QnLex>::new(keys.clone());
        for b in &bodies {
            let off = arena.extend_from_slice(b) as u64;
            strategy.push_record(b, off, u32::try_from(b.len()).unwrap()).unwrap();
        }
        let MemoryChunkErased::QuerynameLex(chunk) =
            strategy.seal(Arc::new(PooledSegmentedBuf::unpooled(arena)), 2)
        else {
            panic!("a lex strategy seals a QuerynameLex chunk");
        };
        assert_eq!(keys.table.len(), 2, "the first batch learned one class");
        assert_eq!(keys.table.ecp(1), b"A00132:53:HFHJ:1:");
        assert!(
            (0..chunk.len()).all(|i| chunk.key_at(i).class == 1),
            "every ref, the first batch's included, is in the learned class"
        );
        let s = keys.stats.snapshot();
        assert_eq!((s.batches, s.reclassified, s.class0_refs), (2, 1, 0));
        let mut want = bodies.clone();
        want.sort_by(|a, b| {
            RawQuerynameLexKey::extract_from_record(a)
                .cmp(&RawQuerynameLexKey::extract_from_record(b))
        });
        let got: Vec<Vec<u8>> = (0..chunk.len()).map(|i| chunk.record_bytes(i).to_vec()).collect();
        assert!(got == want, "the sealed chunk is in owned-key order");
    }

    /// Build a one-block `ArenaBlock` over a fresh unpooled arena with its slot
    /// reserved, for `inflate_one` tests. Returns the token and the arena so the
    /// caller can read the slot back.
    fn arena_block_for(
        block: Vec<u8>,
        payload_len: usize,
    ) -> (ArenaBlock, Arc<PooledSegmentedBuf>) {
        let mut arena = SegmentedBuf::with_capacity(0, 1 << 20);
        arena.reserve_full_capacity();
        // Reserve the slot with the safe API (zero placeholder); `inflate_one`
        // fully overwrites it before any read, so the placeholder is never
        // observed. `extend_from_slice` returns the write offset to use.
        let offset = arena.extend_from_slice(&vec![0u8; payload_len]) as u64;
        let arena = Arc::new(PooledSegmentedBuf::unpooled(arena));
        let item = ArenaBlock {
            arena: Arc::clone(&arena),
            ordinal: 0,
            offset,
            len: u32::try_from(payload_len).unwrap(),
            block: block.into(),
            is_last_of_run: true,
            run_seq: 0,
            seals_to_spill: false,
        };
        (item, arena)
    }

    /// Flip one footer-CRC32 bit of a raw BGZF block (the payload is untouched).
    fn with_flipped_crc(mut block: Vec<u8>) -> Vec<u8> {
        let crc_off = block.len() - fgumi_bgzf::BGZF_FOOTER_SIZE;
        block[crc_off] ^= 0x01;
        block
    }

    /// `InflateToArena::new_with_crc` decides whether a CRC-only-corrupted block
    /// is rejected or inflated into the slot. Covered on both decode branches
    /// (level 6 → libdeflater; level 0 → stored fast path) because the CRC is
    /// checked at two separate call sites in fgumi-bgzf. The plain `new`
    /// constructor must behave as `verify_crc = true`.
    #[rstest]
    #[case::deflate_verify_rejects(6, Some(true), false)]
    #[case::deflate_skip_accepts(6, Some(false), true)]
    #[case::stored_verify_rejects(0, Some(true), false)]
    #[case::stored_skip_accepts(0, Some(false), true)]
    #[case::plain_new_verifies(6, None, false)]
    fn inflate_honors_verify_crc(
        #[case] level: u32,
        #[case] verify_crc: Option<bool>,
        #[case] expect_ok: bool,
    ) {
        let payload = b"CRC-POLICY-ARENA-TEST".repeat(40);
        let block = match level {
            0 => make_test_bgzf_block(&payload),
            _ => compress_one_block(&payload),
        };
        let (item, arena) = arena_block_for(with_flipped_crc(block), payload.len());
        let mut step = match verify_crc {
            Some(v) => InflateToArena::new_with_crc(64 * 1024 * 1024, v),
            None => InflateToArena::new(64 * 1024 * 1024),
        };
        let result = step.inflate_one(item);
        if expect_ok {
            let inflated = result.expect("verify_crc=false must inflate past a CRC-only fault");
            assert_eq!(usize::try_from(inflated.len).unwrap(), payload.len());
            assert_eq!(arena.slice(0, payload.len()), &payload[..], "slot holds the payload");
        } else {
            // `InflatedBlock` is not `Debug`, so match rather than `expect_err`.
            let Err(err) = result else {
                panic!("verify_crc=true must reject a CRC32 mismatch");
            };
            assert!(err.to_string().contains("CRC32"), "got: {err}");
        }
    }

    /// Parallel workers are cloned via `new_worker_copy`; the policy must ride
    /// along or only the seed worker would honor `--no-check-crc`.
    #[rstest]
    #[case::skip(false)]
    #[case::verify(true)]
    fn new_worker_copy_propagates_verify_crc(#[case] verify_crc: bool) {
        let seed = InflateToArena::new_with_crc(1024, verify_crc);
        let mut copy = seed.new_worker_copy();
        let payload = b"WORKER-COPY".repeat(20);
        let (item, _arena) =
            arena_block_for(with_flipped_crc(compress_one_block(&payload)), payload.len());
        assert_eq!(copy.inflate_one(item).is_ok(), !verify_crc);
    }

    /// What a probe strategy saw while it sorted.
    #[derive(Debug, PartialEq, Eq)]
    struct SealSeen {
        /// Threads the seal was asked to sort on.
        threads: usize,
        /// `cap.active()` during the sort.
        active: usize,
        /// Whether a pool step's `try_acquire` was refused during the sort.
        pool_refused: bool,
        /// Whether the pool holders had already dropped their permits.
        pool_released: bool,
    }

    /// A strategy that records, at each seal, what [`SealSeen`] reports.
    /// `parallel` is what it answers to `sorts_in_parallel`.
    struct ActiveProbe {
        cap: Arc<PhaseCap>,
        parallel: bool,
        released: Arc<AtomicBool>,
        seen: Vec<SealSeen>,
    }

    impl ActiveProbe {
        fn new(cap: &Arc<PhaseCap>, parallel: bool) -> Self {
            Self {
                cap: Arc::clone(cap),
                parallel,
                released: Arc::new(AtomicBool::new(true)),
                seen: Vec::new(),
            }
        }
    }

    impl ArenaSortStrategy for ActiveProbe {
        fn push_record(&mut self, _body: &[u8], _body_off: u64, _len: u32) -> io::Result<()> {
            Ok(())
        }
        fn seal(
            &mut self,
            _arena: Arc<PooledSegmentedBuf>,
            sort_threads: usize,
        ) -> MemoryChunkErased {
            let pool_released = self.released.load(Ordering::Acquire);
            let active = self.cap.active();
            let pool_refused = self.cap.try_acquire().is_none();
            self.seen.push(SealSeen { threads: sort_threads, active, pool_refused, pool_released });
            MemoryChunkErased::Coordinate(fgumi_sort::InMemoryChunk::default())
        }
        fn sorts_in_parallel(&self, _sort_threads: usize) -> bool {
            self.parallel
        }
    }

    /// A capped parallel per-run sort holds the WHOLE phase-1 cap while it
    /// sorts at the cap's width: it starts only after the pool holders dropped
    /// their permits (the dropper sets `released` BEFORE releasing, so the
    /// sort can only see `true` if it waited), sees every permit held, and a
    /// pool step's `try_acquire` is refused while it sorts. Everything is
    /// released afterwards.
    #[test]
    fn capped_parallel_seal_holds_the_whole_cap_while_it_sorts() {
        let cap = PhaseCap::new("sort-phase1", 3);
        // `Pipeline::run` binds every reported cap; this drives the seal directly.
        let signal = fgumi_pipeline_core::signal::PipelineSignal::new();
        cap.bind_signal_for_test(&signal);
        let mut probe = ActiveProbe::new(&cap, true);
        let whole =
            |pool_released| SealSeen { threads: 3, active: 3, pool_refused: true, pool_released };

        // An idle cap: the whole cap, at once.
        seal_under_cap(&mut probe, Some(&cap), test_arena(), 3).unwrap();
        assert_eq!(probe.seen.pop(), Some(whole(true)));
        assert_eq!(cap.active(), 0, "released after the sort");

        // Pool steps hold permits: the sort waits for both to be released,
        // then holds all three.
        probe.released.store(false, Ordering::Release);
        let released = Arc::clone(&probe.released);
        std::thread::scope(|scope| {
            let held: Vec<_> = (0..2).map(|_| cap.try_acquire().expect("pool")).collect();
            scope.spawn(move || {
                std::thread::sleep(std::time::Duration::from_millis(20));
                released.store(true, Ordering::Release);
                drop(held);
            });
            seal_under_cap(&mut probe, Some(&cap), test_arena(), 3).unwrap();
        });
        assert_eq!(probe.seen.pop(), Some(whole(true)), "waited for the pool, then held all");
        assert_eq!(cap.active(), 0);
    }

    /// A run that sorts on one thread does not take the whole cap — it would
    /// stall every pool step for serial work — and neither does an uncapped
    /// seal. Both sort at the configured width with no permit taken.
    #[rstest]
    #[case::capped_serial(true, false)]
    #[case::uncapped_parallel(false, true)]
    fn seal_takes_no_permit_when_serial_or_uncapped(#[case] capped: bool, #[case] parallel: bool) {
        let cap = PhaseCap::new("sort-phase1", 3);
        let mut probe = ActiveProbe::new(&cap, parallel);
        let _pool = cap.try_acquire().expect("a pool step holds a permit");
        seal_under_cap(&mut probe, capped.then_some(&*cap), test_arena(), 3).unwrap();
        assert_eq!(
            probe.seen.pop(),
            Some(SealSeen { threads: 3, active: 1, pool_refused: false, pool_released: true }),
            "did not wait for the pool holder, held nothing, did not block admission"
        );
    }

    /// Cancel while a capped seal waits for the whole cap: the seal gives up
    /// with `Interrupted` and the strategy never sorts.
    #[test]
    fn cancel_while_a_seal_waits_skips_the_sort() {
        let cap = PhaseCap::new("sort-phase1", 2);
        let signal = fgumi_pipeline_core::signal::PipelineSignal::new();
        cap.bind_signal_for_test(&signal);
        let held = cap.try_acquire().expect("a pool step holds a permit");
        let (tx, rx) = std::sync::mpsc::channel();
        let seal_cap = Arc::clone(&cap);
        let waiter = std::thread::spawn(move || {
            let mut probe = ActiveProbe::new(&seal_cap, true);
            let result = seal_under_cap(&mut probe, Some(&seal_cap), test_arena(), 2);
            let _ = tx.send(result.map(|_| ()).map_err(|e| e.kind()));
            probe.seen.len()
        });
        while cap.pending_reservations() == 0 {
            std::thread::yield_now();
        }
        signal.cancel();
        let result = rx.recv_timeout(std::time::Duration::from_secs(10)).expect("woken by cancel");
        assert_eq!(result, Err(io::ErrorKind::Interrupted));
        assert_eq!(waiter.join().unwrap(), 0, "no sort ran");
        assert_eq!(cap.pending_reservations(), 0);
        drop(held);
    }

    /// Inflate under `--sort-threads`: the phase-1 admission contract. The
    /// arena backing the block stays alive for the whole check.
    #[test]
    fn capped_step_follows_the_admission_contract() {
        let payload = b"CAPPED-INFLATE".repeat(20);
        let (item, _arena) = arena_block_for(compress_one_block(&payload), payload.len());
        fgumi_pipeline_core::testing::assert_admission_contract(
            "sort-phase1",
            |cap| InflateToArena::new(1 << 20).with_phase_cap(cap),
            item,
        );
    }
}

//! `SortSpillDecompress` — Parallel typed step that reads spill chunk
//! files, decompresses their blocks, and pushes the decompressed bytes
//! into per-slot bounded queues on `SortMergeSlot`.
//!
//! # Two decompression granularities
//!
//! The step supports two strategies, selected by [`SortDecompressTuning`]:
//!
//! - **file-granularity (`file_granularity == true`, the fallback):** one worker
//!   owns a file's decompression at a time. Under the per-slot reader lock it
//!   reads AND decompresses up to `block_batch` blocks inline, in read order,
//!   and pushes them to the slot's FIFO. No reorder buffer is needed — a plain
//!   FIFO suffices because read-and-decompress is a single inline operation. This
//!   is the proven path (see the `merge_slots` module header, "What used to live
//!   here, and why it's gone (v4 vs v3.1)").
//!
//! - **block-parallel (`file_granularity == false`):** multiple workers
//!   decompress different blocks of the SAME file concurrently. Each `try_run`
//!   holds the reader lock only for the READ (sequence-tagging each raw block via
//!   `SortMergeReader::next_seq`), releases it, then decompresses its own batch
//!   OUTSIDE the lock and reassembles via the slot's `ReorderBuffer`. The read
//!   and decompression of a given block still happen within a single `try_run`
//!   of a single worker — the lock is merely released between them. Parallelism
//!   comes from multiple workers each grabbing the lock briefly, reading their
//!   own batch, and decompressing concurrently — NOT from splitting read and
//!   decompress across dispatches (which would re-open the v3 Skip-wedge
//!   deadlock).
//!
//! # HARD INVARIANT
//!
//! A spill block must be read AND decompressed within a single `try_run` by a
//! single worker. Both paths uphold this.
//!
//! # Memory note
//!
//! `--max-memory` does NOT bound Phase-2 decompressed memory today: the FIFO is
//! count-bounded (`PHASE2_DECOMP_CAP`). The block-parallel path's reorder window
//! is the additional decompressed-memory surface a slow straggler could grow, so
//! it is explicitly bounded per-slot by `window_budget` (derived from the step's
//! `output_byte_limit`) via [`SortMergeSlot::bp_reorder_admits`].

use std::io;
use std::sync::Arc;

use fgumi_sort::{SortMergeSlot, SpillBlockDecompressor};
use parking_lot::Mutex;

use crate::sort::protocol::{SortPhase1Event, SortPhase2Event};
use fgumi_pipeline_core::{
    PhaseCap, Unpushed,
    held::HeldSlot,
    outputs::Single,
    queues::QueueSpec,
    reorder::BranchOrdering,
    step::{Step, StepCtx, StepKind, StepOutcome, StepProfile},
};

/// Default reorder-window byte budget substituted when the caller's
/// `output_byte_limit` is `0`. A `0` memory limit is normalized to a fixed cap
/// rather than treated as "unlimited": `SortMergeSlot::bp_reorder_admits`
/// (via `ReorderBuffer::would_accept`) reads `window_budget == 0` as *no bound*,
/// so passing a resolved-to-zero budget straight through would remove the only
/// byte cap on decompressed stragglers. Matches
/// `fgumi_pipeline_core::reorder::DEFAULT_REORDER_OVERFLOW_BYTES` (256 MiB).
const DEFAULT_REORDER_WINDOW_BYTES: u64 = 256 * 1024 * 1024;

/// Tuning for the Phase-2 spill decompression granularity.
///
/// Threaded from `SortOptions` (`--sort::file-granularity` /
/// `--sort::block-batch`) through `ChainBuilder::add_sort` into
/// [`SortSpillDecompress::new`].
#[derive(Debug, Clone, Copy)]
pub struct SortDecompressTuning {
    /// `false` (default) — block-level parallel: hold the reader lock only for
    /// the READ (sequence-tagged), release, decompress OUTSIDE the lock; multiple
    /// workers decompress one file's blocks concurrently, reassembled by a
    /// `ReorderBuffer`. The hardened production default (loom + soak matrix).
    ///
    /// `true` — one worker owns a file's decompression at a time (inline under
    /// the reader lock, in-order, plain FIFO): the older single-worker-per-file
    /// fallback.
    pub file_granularity: bool,
    /// Number of raw blocks claimed per reader-lock acquisition (replaces the
    /// formerly-hardcoded batch size). Default `4` (restores the original
    /// `MAX_BATCH_PER_CALL`; a fleet decompress-throughput bench will pick the
    /// final value). Must be `>= 1`; [`SortSpillDecompress::new`] clamps lower
    /// values (the single construction chokepoint, so every entry point —
    /// standalone `--block-batch`, `runall`, direct construction — is defended).
    pub block_batch: usize,
}

impl Default for SortDecompressTuning {
    fn default() -> Self {
        // Block-parallel is the production default: it cleared the hardening gate
        // (loom over the real `SortMergeSlot` + the external-watchdog soak matrix
        // + the reorder-window byte cap). `file_granularity = true` is the
        // single-worker-per-file fallback. `block_batch` default is 4 (the
        // original `MAX_BATCH_PER_CALL`); a fleet bench will tune it. Matches
        // `SortOptions::default`.
        Self { file_granularity: false, block_batch: 4 }
    }
}

/// What [`SortSpillDecompress::scan_slots`] found.
#[derive(Clone, Copy)]
struct SlotScan {
    /// Some registered slot has not reached its queue EOF.
    alive: bool,
    /// Some live slot has FIFO room, so a fill could do work.
    fillable: bool,
}

/// What one admission-controlled fill attempt did.
enum Fill {
    /// A slot was filled.
    Filled,
    /// The phase cap refused this worker.
    Refused,
    /// Nothing to fill right now (or no slot is live).
    Nothing,
}

struct RegisteredSpill {
    slot: Arc<SortMergeSlot>,
}

/// Parallel step that reads + decompresses spill chunk files and
/// pushes results into per-slot queues. Forwards `SortPhase1Event`s
/// verbatim to `SortMerge`.
pub struct SortSpillDecompress {
    registry: Arc<Mutex<Vec<RegisteredSpill>>>,
    block_dec: SpillBlockDecompressor,
    held: HeldSlot<Unpushed<SortPhase2Event>>,
    output_byte_limit: u64,
    /// Phase-2 admission cap shared with the other merge-phase pool steps
    /// (`None` = uncapped until the chain builder sets it).
    cap: Option<Arc<PhaseCap>>,
    tuning: SortDecompressTuning,
    /// Per-slot reorder-window byte budget for the block-parallel path. Derived
    /// from `output_byte_limit` (one per-step byte budget per slot). Bounds the
    /// reorder buffer so a slow straggler can't balloon decompressed memory.
    window_budget: u64,
    /// The sort's merge-wide demand. After every delivery, EOF finalize or
    /// failure on a slot this step calls `notify_delivered(slot.file_id)`, which
    /// unparks the merge iff it is waiting on that file. `None` when the step
    /// runs without a merge to wake (unit tests that drive the step alone).
    demand: Option<Arc<fgumi_sort::MergeDemand>>,
}

impl SortSpillDecompress {
    /// Construct a fresh step with an empty registry.
    ///
    /// `output_byte_limit` byte-bounds the forwarded-event output queue.
    /// The forwarded `SortPhase2Event::MemoryChunk` variant retains sorted
    /// record chunks, so this queue must budget on bytes (`HeapSize`), not
    /// event count, to keep retained memory a function of configuration.
    ///
    /// `tuning` selects the decompression granularity (see
    /// [`SortDecompressTuning`]). The block-parallel path derives its per-slot
    /// reorder-window budget from `output_byte_limit`.
    #[must_use]
    pub fn new(output_byte_limit: u64, tuning: SortDecompressTuning) -> Self {
        // Clamp `block_batch` to >= 1. A value of 0 reads zero blocks per
        // acquisition, which on the inline path declares a phantom EOF after
        // reading nothing (silent record loss) and on the block-parallel path
        // never sets `reader_eof`/`queue_eof` (the merge livelocks). 0 is
        // nonsensical for a "blocks per read" knob, so we normalize rather than
        // propagate it. This is the single construction chokepoint for the step,
        // so it defends every entry point (CLI, runall, direct construction).
        let tuning = SortDecompressTuning { block_batch: tuning.block_batch.max(1), ..tuning };
        // Normalize a zero reorder-window budget to a sane default rather than
        // propagating it: `bp_reorder_admits` treats `window_budget == 0` as
        // "unlimited", so a resolved-to-zero budget (e.g. `--max-memory 0`) would
        // silently remove the byte cap on decompressed stragglers. This mirrors
        // the legacy `effective_limit` 0-normalization and is applied at the same
        // construction chokepoint as the `block_batch` clamp, defending every
        // entry point (CLI, runall, direct construction).
        let window_budget =
            if output_byte_limit == 0 { DEFAULT_REORDER_WINDOW_BYTES } else { output_byte_limit };
        Self {
            registry: Arc::new(Mutex::new(Vec::new())),
            block_dec: SpillBlockDecompressor::new(),
            held: HeldSlot::new(),
            output_byte_limit,
            cap: None,
            tuning,
            window_budget,
            demand: None,
        }
    }

    /// Share the sort's [`fgumi_sort::MergeDemand`] (one per sort, built by the
    /// chain builder's `add_sort`), so deliveries wake the merge awaiting them.
    /// On the inline (file-granularity) path it also marks the demand's stalls
    /// as unclassified (`MergeDemandStats::mark_inline_decompress`).
    #[must_use]
    pub fn with_merge_demand(mut self, demand: Arc<fgumi_sort::MergeDemand>) -> Self {
        if self.tuning.file_granularity {
            demand.stats().mark_inline_decompress();
        }
        self.demand = Some(demand);
        self
    }

    /// The merge demand this step notifies (tests: a clone must carry it).
    #[cfg(test)]
    #[must_use]
    pub(crate) fn merge_demand_for_test(&self) -> Option<&Arc<fgumi_sort::MergeDemand>> {
        self.demand.as_ref()
    }

    /// One pass over the registry, under its lock and without cloning it:
    /// whether any registered slot is still live (has not reached its queue
    /// EOF), and whether any live slot has FIFO room — the only state in which
    /// a fill (either path) can do work.
    fn scan_slots(&self) -> SlotScan {
        use std::sync::atomic::Ordering;
        let mut scan = SlotScan { alive: false, fillable: false };
        for entry in self.registry.lock().iter() {
            let slot = &entry.slot;
            if slot.queue_eof.load(Ordering::Acquire) {
                continue;
            }
            scan.alive = true;
            if slot.bp_fifo_room() > 0 {
                scan.fillable = true;
                break;
            }
        }
        scan
    }

    /// Whether any registered slot has not reached its queue EOF — the
    /// uncapped path's liveness check after a failed fill: one atomic load per
    /// slot under the registry lock, no FIFO inspection.
    fn any_slot_alive(&self) -> bool {
        self.registry
            .lock()
            .iter()
            .any(|entry| !entry.slot.queue_eof.load(std::sync::atomic::Ordering::Acquire))
    }

    /// One fill attempt under `cap` (moved out of `self` by the caller). A
    /// poll with no fillable slot — none live, or every live slot's FIFO full —
    /// takes no permit and is not a refusal, so idle clones never crowd the
    /// output compressor sharing the phase-2 cap; otherwise the permit is held
    /// across the fill. The uncapped path does not come here: it fills
    /// directly and pays no scan.
    fn fill_under_cap(&mut self, cap: &PhaseCap, scan: SlotScan) -> io::Result<Fill> {
        if !scan.fillable {
            return Ok(Fill::Nothing);
        }
        let Some(_permit) = cap.try_acquire() else {
            return Ok(Fill::Refused);
        };
        Ok(if self.try_fill_some_slot()? { Fill::Filled } else { Fill::Nothing })
    }

    /// Share the phase-2 admission cap with the other merge-phase pool steps
    /// (`--merge-threads`): at most `cap.max()` workers decompress spill
    /// blocks / compress output at once.
    #[must_use]
    pub fn with_phase_cap(mut self, cap: Option<Arc<PhaseCap>>) -> Self {
        self.cap = cap;
        self
    }

    fn flush_held(&mut self, ctx: &mut StepCtx<'_, Self>) -> bool {
        let Some(unpushed) = self.held.take() else {
            return true;
        };
        match ctx.outputs.retry(unpushed) {
            Ok(()) => true,
            Err(again) => {
                self.held.put(again);
                false
            }
        }
    }

    fn push_or_hold(&mut self, ctx: &mut StepCtx<'_, Self>, event: SortPhase2Event) -> bool {
        match ctx.outputs.push(event) {
            Ok(()) => true,
            Err(unpushed) => {
                self.held.put(unpushed);
                false
            }
        }
    }

    fn snapshot_registry(&self) -> Vec<Arc<SortMergeSlot>> {
        let registry = self.registry.lock();
        registry.iter().map(|e| Arc::clone(&e.slot)).collect()
    }

    /// Slot indices ordered by ascending FIFO block count (most-starved first) — the
    /// emptiest-first refill forecaster (see the budget-refill design spec §4.1).
    /// Snapshots each slot's [`SortMergeSlot::fifo_len`] once (O(N) brief locks), then
    /// sorts the indices, so the per-slot lock is taken exactly once per dispatch — not
    /// inside the sort comparator. For the typical spill count (tens) this is negligible
    /// next to the decompression work; if N grows into the hundreds, profile (a relaxed
    /// cached length on the slot is the fallback) before keeping it.
    ///
    /// Slots that have already signalled `queue_eof` are dropped before the FIFO
    /// lengths are read. The registry is append-only, so a drained slot would
    /// otherwise stay in the scan for the rest of the run — taking its FIFO lock
    /// on every dispatch only to be rejected immediately by
    /// `try_fill_inline_slot` / `try_fill_block_parallel_slot`. With many spill
    /// files that drained tail dominates the scan. Filtering is scheduling-only:
    /// an EOF slot can never make progress, so skipping it changes no output.
    #[must_use]
    pub(crate) fn emptiest_first_order(slots: &[Arc<SortMergeSlot>]) -> Vec<usize> {
        use std::sync::atomic::Ordering;

        let mut order: Vec<usize> =
            (0..slots.len()).filter(|&i| !slots[i].queue_eof.load(Ordering::Acquire)).collect();
        // `sort_by_cached_key`, not `sort_by_key`: the key takes the slot's FIFO
        // lock, and `sort_by_key` would re-take it O(n log n) times. This way each
        // surviving slot is locked exactly once, and EOF slots not at all.
        order.sort_by_cached_key(|&i| slots[i].fifo_len());
        order
    }

    fn try_fill_some_slot(&mut self) -> io::Result<bool> {
        let slots = self.snapshot_registry();
        // Refill the most-starved slot first so a free worker tops up the slot the merge
        // will exhaust soonest, rather than the first in registry order. Scheduling-only:
        // admission and the read-and-decompress-in-one-`try_run` invariant are unchanged,
        // so this cannot affect output or wedge progress (a non-progressing slot returns
        // `false` fast and the loop falls through to the next).
        for i in Self::emptiest_first_order(&slots) {
            let slot = &slots[i];
            let progressed = if self.tuning.file_granularity {
                self.try_fill_inline_slot(slot)?
            } else {
                self.try_fill_block_parallel_slot(slot)?
            };
            if progressed {
                return Ok(true);
            }
        }
        Ok(false)
    }

    /// Inline (file-granularity) fill: read AND decompress up to `block_batch`
    /// blocks under the reader lock, push them to the FIFO in read order. One
    /// worker owns a slot at a time; no reorder buffer needed.
    fn try_fill_inline_slot(&mut self, slot: &Arc<SortMergeSlot>) -> io::Result<bool> {
        use std::sync::atomic::Ordering;

        if slot.queue_eof.load(Ordering::Acquire) {
            return Ok(false);
        }

        let mut reader_guard = match slot.reader.try_lock() {
            Ok(guard) => guard,
            // Contended (another worker owns the slot): a normal skip.
            Err(std::sync::TryLockError::WouldBlock) => return Ok(false),
            // Poisoned: a fill worker panicked mid-read. Fail the slot CLOSED so
            // `SortMerge` surfaces the failure; swallowing it as a skip would
            // leave `queue_eof` unset and spin `Contention` forever (deadlock).
            Err(std::sync::TryLockError::Poisoned(_)) => {
                self.mark_slot_failed(slot);
                return Err(io::Error::other(
                    "spill reader mutex poisoned: a decompress fill worker panicked",
                ));
            }
        };

        let room = {
            let dec = slot.decompressed.lock().expect("decompressed mutex poisoned");
            fgumi_sort::PHASE2_DECOMP_CAP.saturating_sub(dec.len())
        };
        if room == 0 {
            return Ok(false);
        }
        let want = room.min(self.tuning.block_batch);

        let decompressed_batch =
            match self.block_dec.read_blocks(&mut reader_guard.inner, slot.codec, want) {
                Ok(b) => b,
                Err(e) => {
                    // Centralized in `mark_slot_failed` so failure semantics stay
                    // in one place (see the block-parallel path's use of it).
                    self.mark_slot_failed(slot);
                    drop(reader_guard);
                    return Err(e);
                }
            };
        let got = decompressed_batch.len();
        let hit_eof = got < want;

        if got == 0 {
            {
                let _g = slot.decompressed.lock().expect("decompressed mutex poisoned");
                slot.queue_eof.store(true, Ordering::Release);
            }
            drop(reader_guard);
            self.notify(slot);
            return Ok(true);
        }

        {
            let mut dec = slot.decompressed.lock().expect("decompressed mutex poisoned");
            for b in decompressed_batch {
                dec.push_back(b);
            }
            if hit_eof {
                slot.queue_eof.store(true, Ordering::Release);
            }
        }
        drop(reader_guard);
        self.notify(slot);
        Ok(true)
    }

    /// Block-parallel fill: under the reader lock read (only) up to `block_batch`
    /// raw blocks, sequence-tag them, release the lock, decompress OUTSIDE the
    /// lock, then reassemble via the slot's reorder buffer and drain in-order
    /// blocks into the FIFO. Multiple workers run this concurrently on the same
    /// slot.
    fn try_fill_block_parallel_slot(&mut self, slot: &Arc<SortMergeSlot>) -> io::Result<bool> {
        use std::sync::atomic::Ordering;

        if slot.queue_eof.load(Ordering::Acquire) {
            return Ok(false);
        }

        // Phase A: read a fresh batch if the reader is still live and the
        // FIFO / reorder window admit more.
        if !slot.reader_eof.load(Ordering::Acquire) {
            let acquired = match slot.reader.try_lock() {
                Ok(guard) => Some(guard),
                // Contended: fall through to the drain-only phase below.
                Err(std::sync::TryLockError::WouldBlock) => None,
                // Poisoned: a fill worker panicked mid-read. Fail closed so
                // `SortMerge` surfaces it instead of spinning forever.
                Err(std::sync::TryLockError::Poisoned(_)) => {
                    self.mark_slot_failed(slot);
                    return Err(io::Error::other(
                        "spill reader mutex poisoned: a decompress fill worker panicked",
                    ));
                }
            };
            if let Some(mut reader_guard) = acquired {
                // Re-check under the lock: another worker may have hit EOF.
                if !slot.reader_eof.load(Ordering::Acquire) {
                    let next_seq = reader_guard.next_seq;
                    let fifo_room = slot.bp_fifo_room();
                    let admit =
                        fifo_room > 0 && slot.bp_reorder_admits(next_seq, self.window_budget);
                    // NB: the reorder-window budget is checked once here (for
                    // `next_seq`), then up to `want` (≤ `block_batch`) blocks are
                    // inserted below without a per-block re-check. So the reorder
                    // window can transiently exceed `window_budget` by up to
                    // `block_batch - 1` blocks. This overshoot is bounded and by
                    // design: `block_batch` is small (default 4) and configurable,
                    // so worst-case resident bytes stay `O(window_budget +
                    // block_batch × block_size)` — not the unbounded growth the
                    // window guards against. Per-block admission is intentionally
                    // avoided to keep the reader-lock hold short (read the whole
                    // batch, release, decompress outside the lock).
                    if admit {
                        // Bound the read by FIFO room (as the inline path does):
                        // reading `block_batch` when only `fifo_room < block_batch`
                        // slots can drain would over-admit the surplus into the
                        // reorder window. `want >= 1` since `fifo_room > 0`.
                        let want = self.tuning.block_batch.min(fifo_room);
                        let start_seq = reader_guard.next_seq;
                        let raw = match self.block_dec.read_raw(
                            &mut reader_guard.inner,
                            slot.codec,
                            want,
                        ) {
                            Ok(r) => r,
                            Err(e) => {
                                self.mark_slot_failed(slot);
                                drop(reader_guard);
                                return Err(e);
                            }
                        };
                        let got = raw.len();
                        // EOF only when the reader returned fewer than we asked
                        // for (`want`); a FIFO-limited short read is not EOF.
                        let hit_eof = got < want;
                        // Stamp the read range and account for it BEFORE releasing
                        // the lock, so a concurrent worker observing EOF cannot
                        // race ahead of this batch's in-flight accounting. The
                        // publish order (reserve `in_flight` before setting
                        // `reader_eof`) is the loom-verified protocol; it lives in
                        // `SortMergeSlot::bp_commit_read` as the single source of
                        // truth, so this call site and the loom model share it (see
                        // that method's doc and fgumi-sort tests/loom_merge_slots.rs).
                        reader_guard.next_seq += got as u64;
                        slot.bp_commit_read(got, hit_eof);
                        drop(reader_guard);

                        // Decompress OUTSIDE the reader lock (still this try_run).
                        let mut blocks = Vec::with_capacity(got);
                        for raw_block in &raw {
                            match self.block_dec.decompress_one(slot.codec, raw_block) {
                                Ok(d) => blocks.push(d),
                                Err(e) => {
                                    self.mark_slot_failed(slot);
                                    return Err(e);
                                }
                            }
                        }
                        // A batch that only filled the reorder window (its front
                        // block is still another worker's) delivers nothing and
                        // so wakes nobody; the worker that closes the gap
                        // notifies.
                        if slot.bp_insert_drain_finalize(start_seq, blocks, got) {
                            self.notify(slot);
                        }
                        return Ok(true);
                    }
                }
            }
        }

        // Phase B: drain-only. Flush any now-in-order blocks the FIFO can accept
        // (it may have freed up, or another worker delivered a straggler) and
        // finalize EOF if fully delivered.
        let progressed = slot.bp_drain_and_finalize();
        if progressed {
            self.notify(slot);
        }
        Ok(progressed)
    }

    /// Mark a slot as failed (decompression / read error): set `decomp_error`
    /// and `queue_eof` under the `decompressed` mutex so the consumer surfaces
    /// the error in preference to a clean EOF, then wake the merge if it awaits
    /// this slot.
    fn mark_slot_failed(&self, slot: &Arc<SortMergeSlot>) {
        use std::sync::atomic::Ordering;
        {
            let _g = slot.decompressed.lock().expect("decompressed mutex poisoned");
            slot.decomp_error.store(true, Ordering::Release);
            slot.queue_eof.store(true, Ordering::Release);
        }
        self.notify(slot);
    }

    /// Tell the merge that `slot` received a block, EOF or failure. Called
    /// after the slot's `decompressed` mutex is released (the ordering the
    /// wake protocol in `fgumi_sort::MergeDemand` relies on); unparks the merge
    /// only when it awaits this slot.
    #[inline]
    fn notify(&self, slot: &SortMergeSlot) {
        if let Some(d) = &self.demand {
            d.notify_delivered(slot.file_id);
        }
    }
}

impl Clone for SortSpillDecompress {
    fn clone(&self) -> Self {
        Self {
            registry: Arc::clone(&self.registry),
            block_dec: SpillBlockDecompressor::new(),
            held: HeldSlot::new(),
            output_byte_limit: self.output_byte_limit,
            // Share the SAME cap so it is global across all worker clones.
            cap: self.cap.clone(),
            tuning: self.tuning,
            window_budget: self.window_budget,
            // Every clone notifies the one demand the merge awaits on.
            demand: self.demand.clone(),
        }
    }
}

impl Step for SortSpillDecompress {
    type Input = SortPhase1Event;
    type Outputs = Single<SortPhase2Event>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "SortSpillDecompress",
            kind: StepKind::Parallel,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::None],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        // 1. Drain held output first.
        if !self.flush_held(ctx) {
            return Ok(StepOutcome::Contention);
        }

        // 2. Pop one input event, register if SpillReady, forward all.
        if let Some(event) = ctx.input.pop() {
            let forwarded = match event {
                SortPhase1Event::SpillReady { slot, path, records_ingested_so_far } => {
                    self.registry.lock().push(RegisteredSpill { slot: Arc::clone(&slot) });
                    SortPhase2Event::SpillReady { slot, path, records_ingested_so_far }
                }
                SortPhase1Event::MemoryChunk { chunk, records_ingested_so_far } => {
                    SortPhase2Event::MemoryChunk { chunk, records_ingested_so_far }
                }
                SortPhase1Event::AllAnnounced { slot_count, memory_chunk_count, total_records } => {
                    SortPhase2Event::AllAnnounced { slot_count, memory_chunk_count, total_records }
                }
            };
            let _ = self.push_or_hold(ctx, forwarded);
            return Ok(StepOutcome::Progress);
        }

        // 3. Greedy slot-fill, admission-controlled by --merge-threads. The
        // permit covers only the fill (spill read + decompress); the event
        // forwarding above is bookkeeping and stays uncapped. Only a capped
        // poll scans the registry (once per call, to decide whether to take a
        // permit); the uncapped path fills directly, as before the cap existed.
        //
        // The cap is moved out for the fill: the permit borrows it while the
        // fill needs `&mut self` (the decompressor), which a field borrow
        // cannot split. A move, so no reference-count traffic.
        let cap = self.cap.take();
        let (fill, scan) = match cap.as_deref() {
            None => {
                let filled = self.try_fill_some_slot();
                (filled.map(|f| if f { Fill::Filled } else { Fill::Nothing }), None)
            }
            Some(cap) => {
                let scan = self.scan_slots();
                (self.fill_under_cap(cap, scan), Some(scan))
            }
        };
        self.cap = cap;
        match fill? {
            Fill::Filled => return Ok(StepOutcome::Progress),
            Fill::Refused => return Ok(StepOutcome::Capped),
            Fill::Nothing => {}
        }

        // 4. No fill work. A capped poll that found every live slot full is
        // idle, not contended; uncapped, a live slot is contention.
        let outcome_if_alive = match scan {
            Some(scan) if scan.alive => {
                Some(if scan.fillable { StepOutcome::Contention } else { StepOutcome::NoProgress })
            }
            Some(_) => None,
            None => self.any_slot_alive().then_some(StepOutcome::Contention),
        };
        if let Some(outcome) = outcome_if_alive {
            return Ok(outcome);
        }

        if ctx.input.is_drained() {
            return Ok(StepOutcome::Finished);
        }
        Ok(StepOutcome::NoProgress)
    }

    fn new_worker_copy(&self) -> Self {
        self.clone()
    }

    fn phase_cap(&self) -> Option<&PhaseCap> {
        self.cap.as_deref()
    }
}

#[cfg(test)]
mod tests;

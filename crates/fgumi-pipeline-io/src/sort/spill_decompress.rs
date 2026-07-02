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
use std::sync::atomic::{AtomicUsize, Ordering};

use fgumi_sort::{SortMergeSlot, SpillBlockDecompressor};
use parking_lot::Mutex;

use crate::sort::protocol::{SortPhase1Event, SortPhase2Event};
use fgumi_pipeline_core::{
    Unpushed,
    held::HeldSlot,
    outputs::Single,
    queues::QueueSpec,
    reorder::BranchOrdering,
    step::{Step, StepCtx, StepKind, StepOutcome, StepProfile},
};

/// Default reorder-window byte budget substituted when the caller's
/// `output_byte_limit` is `0`. Mirrors the legacy pipeline's `effective_limit`
/// (`unified_pipeline`), which normalizes a `0` memory limit to a fixed cap
/// rather than treating it as "unlimited": `SortMergeSlot::bp_reorder_admits`
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
    /// values and the CLI rejects them.
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

/// CAS-acquire one decompress permit. `max == None` ⇒ unbounded (always succeeds).
///
/// Returns the [`DecompressPermit`] on success so ownership of the slot is
/// encoded in the type: the count is released exactly once, when the returned
/// permit drops, and acquisition cannot be separated from release. `None` ⇒ the
/// cap is full and no slot was taken.
#[must_use]
fn try_acquire(active: &Arc<AtomicUsize>, max: Option<usize>) -> Option<DecompressPermit> {
    let Some(max) = max else {
        active.fetch_add(1, Ordering::AcqRel);
        return Some(DecompressPermit { active: Arc::clone(active) });
    };
    let max = max.max(1);
    let mut cur = active.load(Ordering::Acquire);
    loop {
        if cur >= max {
            return None;
        }
        match active.compare_exchange_weak(cur, cur + 1, Ordering::AcqRel, Ordering::Acquire) {
            Ok(_) => return Some(DecompressPermit { active: Arc::clone(active) }),
            Err(observed) => cur = observed,
        }
    }
}

/// RAII release of one decompress permit. Holds an OWNED `Arc` clone so it does
/// not borrow `self` while `try_fill_some_slot(&mut self)` runs.
struct DecompressPermit {
    active: Arc<AtomicUsize>,
}
impl Drop for DecompressPermit {
    fn drop(&mut self) {
        self.active.fetch_sub(1, Ordering::AcqRel);
    }
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
    /// Shared admission counter: how many worker clones are currently inside the
    /// decompress branch. Shared across all clones so the cap is global.
    active: Arc<AtomicUsize>,
    /// Maximum concurrent decompress workers. `None` ⇒ unbounded. Set from
    /// `--merge-threads` (Phase-2 CPU control).
    max_decompress: Option<usize>,
    tuning: SortDecompressTuning,
    /// Per-slot reorder-window byte budget for the block-parallel path. Derived
    /// from `output_byte_limit` (one per-step byte budget per slot). Bounds the
    /// reorder buffer so a slow straggler can't balloon decompressed memory.
    window_budget: u64,
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
            active: Arc::new(AtomicUsize::new(0)),
            max_decompress: None,
            tuning,
            window_budget,
        }
    }

    /// Cap the number of concurrent spill-decompression workers (Phase-2 CPU
    /// control via `--merge-threads`). `None` (default) leaves it unbounded.
    #[must_use]
    pub fn with_max_concurrency(mut self, max: Option<usize>) -> Self {
        self.max_decompress = max;
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
    #[must_use]
    pub(crate) fn emptiest_first_order(slots: &[Arc<SortMergeSlot>]) -> Vec<usize> {
        let lens: Vec<usize> = slots.iter().map(|s| s.fifo_len()).collect();
        let mut order: Vec<usize> = (0..slots.len()).collect();
        order.sort_by_key(|&i| lens[i]);
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

        let Ok(mut reader_guard) = slot.reader.try_lock() else {
            return Ok(false);
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
                    Self::mark_slot_failed(slot);
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
            if let Ok(mut reader_guard) = slot.reader.try_lock() {
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
                                Self::mark_slot_failed(slot);
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
                                    Self::mark_slot_failed(slot);
                                    return Err(e);
                                }
                            }
                        }
                        slot.bp_insert_drain_finalize(start_seq, blocks, got);
                        return Ok(true);
                    }
                }
            }
        }

        // Phase B: drain-only. Flush any now-in-order blocks the FIFO can accept
        // (it may have freed up, or another worker delivered a straggler) and
        // finalize EOF if fully delivered.
        Ok(slot.bp_drain_and_finalize())
    }

    /// Mark a slot as failed (decompression / read error): set `decomp_error`
    /// and `queue_eof` under the `decompressed` mutex so the consumer surfaces
    /// the error in preference to a clean EOF.
    fn mark_slot_failed(slot: &Arc<SortMergeSlot>) {
        use std::sync::atomic::Ordering;
        let _g = slot.decompressed.lock().expect("decompressed mutex poisoned");
        slot.decomp_error.store(true, Ordering::Release);
        slot.queue_eof.store(true, Ordering::Release);
    }
}

impl Clone for SortSpillDecompress {
    fn clone(&self) -> Self {
        Self {
            registry: Arc::clone(&self.registry),
            block_dec: SpillBlockDecompressor::new(),
            held: HeldSlot::new(),
            output_byte_limit: self.output_byte_limit,
            // Share the SAME counter so the cap is global across all worker clones.
            active: Arc::clone(&self.active),
            max_decompress: self.max_decompress,
            tuning: self.tuning,
            window_budget: self.window_budget,
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

        // 3. Greedy slot-fill, admission-controlled by --merge-threads.
        if let Some(_permit) = try_acquire(&self.active, self.max_decompress) {
            if self.try_fill_some_slot()? {
                return Ok(StepOutcome::Progress);
            }
            // permit drops here (also on the `?` early return), releasing the count
        }

        // 4. No fill work.
        let any_alive = self
            .snapshot_registry()
            .iter()
            .any(|slot| !slot.queue_eof.load(std::sync::atomic::Ordering::Acquire));
        if any_alive {
            return Ok(StepOutcome::Contention);
        }

        if ctx.input.is_drained() {
            return Ok(StepOutcome::Finished);
        }
        Ok(StepOutcome::NoProgress)
    }

    fn new_worker_copy(&self) -> Self {
        self.clone()
    }
}

#[cfg(test)]
mod tests;

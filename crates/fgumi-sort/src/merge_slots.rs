//! Per-spill-file shared state for the unified-pipeline merge phase.
//!
//! `SortMergeSlot` is the per-spill-file shared state shuttled between
//! `SortSpillDecompress` (producer side: reads raw BGZF blocks from
//! disk, decompresses inline, pushes to the slot's bounded queue) and
//! `SortMerge` (consumer side: pops decompressed blocks from the
//! queue, parses records, drives the k-way merge) via
//! `Arc<SortMergeSlot>` clones.
//!
//! # Per-slot bounded queue design (v4 — see commit `9c39dea` / PR #389 and
//! `docs/design/sort-phase2-unification-deferral.md`)
//!
//! Each slot carries a bounded queue of decompressed BGZF blocks
//! (`PHASE2_DECOMP_CAP` entries). Backpressure lives here — the
//! producer is **non-blocking**: pushes only when the queue has
//! space; otherwise skips this slot and tries the next. The consumer
//! is also **non-blocking**: when the queue is empty but the slot is
//! not yet `queue_eof` it reports `WouldBlock` (see
//! `external.rs::slot_try_load_block`) and the cooperative `SortMerge`
//! step yields, so the framework re-dispatches it after the producer
//! has refilled the slot. No condvar, no parked thread.
//!
//! ## Why "queue has space OR consumer has a block OR slot EOF" is
//! the full state space
//!
//! For any slot at any wall-clock time, one of the following is true:
//!
//! 1. `decompressed.len() < PHASE2_DECOMP_CAP` — producer can push.
//! 2. `decompressed.len() > 0` — consumer can pop.
//! 3. `queue_eof == true` — slot is done; consumer returns EOF.
//!
//! The "consumer would block" path is reachable only when
//! `decompressed.len() == 0 && !queue_eof`, in which case the
//! producer will eventually flip the state to either (1) (push more
//! blocks) or (3) (set `queue_eof` on the read returning fewer
//! bytes than asked), and the next consumer dispatch observes it. No
//! "transient cap with all-workers-Skip" window.
//!
//! ## Atomic ordering (`decomp_error` / `queue_eof`)
//!
//! Producers' BOTH success and error paths MUST hold the
//! `decompressed` mutex while storing `queue_eof` (and the error
//! path additionally stores `decomp_error`). The consumer always
//! acquires the same mutex at the top of its poll loop. The
//! mutex's release-acquire chain establishes happens-before for
//! BOTH atomics simultaneously — regardless of which the consumer
//! loads first. Without this discipline, a stale `decomp_error`
//! load can race a fresh `queue_eof` load and produce silent
//! truncation.
//!
//! ## What used to live here, and why it's gone (v4 vs v3.1)
//!
//! Pre-v4 the slot also carried a `raw_blocks: Mutex<VecDeque>`
//! queue and a `decomp_in_flight: AtomicUsize` counter. Pre-v4
//! workers split work into a separate "read raw" step and a "claim
//! and decompress" step, with cap+gap-filler admission. That
//! design deadlocked at production scale because the framework's
//! drain protocol could Skip workers during a transient "all slots
//! at cap simultaneously" window. v4 collapses read-and-decompress
//! into one inline operation per worker per slot per `try_run`,
//! eliminating the `raw_blocks` queue, the in-flight counter, the
//! reorder buffer (FIFO suffices because one worker reads per slot
//! at a time via the reader lock), and the gap-filler escape.
//!
//! ## The OTHER Phase-2 implementation (`worker_pool.rs`) still has all of
//! this — by design.
//!
//! Since the P6/P7 unification (#395) THIS module is the production Phase-2 for
//! both standalone `fgumi sort` and the fused `runall` sort. `worker_pool::
//! Phase2FileState` (plus the gap-filler this module dropped) no longer drives
//! a production sort; it is retained as the `RawExternalSorter::sort` library
//! path and the `#[cfg(test)]` parity oracle. It keeps the reorder buffer +
//! in-flight counter because its single-reader/**multi-decompressor** topology
//! genuinely needs them (a plain FIFO suffices HERE only because read-and-
//! decompress is one inline op, so blocks decompress strictly in read order).
//! (History: commit `9d6d7e9` / PR #395 and
//! `docs/design/sort-phase2-unification-deferral.md`.)

use std::collections::VecDeque;
use std::fs::File;
use std::io::BufReader;

// Concurrency primitives are sourced from `loom` under `--cfg loom` so the
// model-checking test (`tests/loom_merge_slots.rs`) exercises the REAL
// `SortMergeSlot` atomics/mutexes — every interleaving and memory reordering of
// the block-parallel EOF/in-flight/finalize protocol — instead of a hand-copied
// re-implementation. Under a normal build these are the `std` types verbatim.
#[cfg(loom)]
use loom::sync::Mutex;
#[cfg(loom)]
use loom::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
#[cfg(not(loom))]
use std::sync::Mutex;
#[cfg(not(loom))]
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

use fgumi_bam_io::reorder::ReorderBuffer;

use crate::codec::SpillCodec;

/// Per-slot decompressed-block queue cap. Bounds in-flight
/// decompressed memory: `num_slots × PHASE2_DECOMP_CAP ×` per-entry size,
/// where each entry is a decompressed BGZF block or zstd frame ranging from
/// ~64 KB (BGZF) up to 256 KB (zstd worst case).
///
/// Raised from 8 to 32 (increment 1a) to give the work-stealing decompressor
/// more runway ahead of the `Detached` merge, which was input-starved
/// (`MergeDiag stalls`) at the spill-heavy operating point. The `worker_pool.rs`
/// copy (the `cfg(test)`/library oracle path) is intentionally left at 8 — these
/// two are no longer equal.
///
/// Hard cap — there is no admission escape. Producers skip a slot
/// when its queue is at this cap, returning to it on a later
/// `try_run` after the consumer has drained. It stays a per-slot, independent
/// bound (deadlock-safety is unchanged — see the module header).
pub const PHASE2_DECOMP_CAP: usize = 32;

/// Disk reader state for a single spill file. Mutex'd separately
/// from `decompressed` so the producer can hold the reader during
/// disk I/O without blocking the consumer's pop.
pub struct SortMergeReader {
    /// Buffered file handle.
    pub inner: BufReader<File>,
    /// Next per-slot sequence number to assign to a raw block read from this
    /// file. Used ONLY by the block-parallel `SortSpillDecompress` path: each
    /// raw block read under the reader lock is stamped with a monotonically
    /// increasing sequence number so the out-of-order decompression results can
    /// be reassembled in read order by the slot's [`SortMergeSlot::reorder`]
    /// buffer. Because every read is serialized under the reader lock, the
    /// sequence numbers are dense and assigned in strict read order. Unused by
    /// the inline (file-granularity) path, which never reorders.
    pub next_seq: u64,
}

/// Per-spill-file shared state.
///
/// Constructed by `SortAndSpill` (one per closed spill file via
/// `slots_for_chunk_files`), embedded in
/// `SortPhase1Event::SpillReady` as `Arc<SortMergeSlot>`, forwarded
/// verbatim by `SortSpillDecompress`, and finally installed in
/// `SortMerge`'s slot table. Drops when the last `Arc` is released
/// — typically after `SortMerge`'s merge driver is exhausted.
pub struct SortMergeSlot {
    /// Stable identifier — the index this slot occupies in the
    /// merge driver's source list. `SortMerge` orders sources by
    /// `file_id` so the `LoserTree` tie-break for equal sort keys
    /// is deterministic and matches the legacy chunk-files order.
    pub file_id: u32,
    /// Spill codec of this chunk's file, detected from the file magic when the
    /// slot is opened (`slots_for_chunk_files`). The `SortSpillDecompress` step
    /// reads it to decompress BGZF blocks or zstd frames; `reader` is already
    /// positioned past any codec file-magic.
    pub codec: SpillCodec,
    /// Disk reader. Held only while reading raw bytes from disk.
    /// One worker at a time per slot via `try_lock`.
    pub reader: Mutex<SortMergeReader>,
    /// Bounded queue of decompressed blocks (each a BGZF block or a zstd
    /// frame, per this slot's `codec`), FIFO. Pushed by the producer
    /// (`SortSpillDecompress`) after read + inline decompress; popped by the
    /// consumer (`SortMerge` via `slot_try_load_block`). Bounded at
    /// `PHASE2_DECOMP_CAP`; producer skips when at cap.
    pub decompressed: Mutex<VecDeque<Vec<u8>>>,
    /// Set true once the producer detects EOF on the disk reader
    /// AND has pushed any final batch of decompressed blocks to
    /// `decompressed`. After this transition, the slot will never
    /// receive another push. Consumer surfaces this as a clean EOF
    /// when its poll finds `decompressed.is_empty() && queue_eof`.
    ///
    /// **Atomic ordering:** producer must hold the `decompressed`
    /// mutex while storing this. Consumer reads it under the same
    /// mutex. Mutex release-acquire creates happens-before.
    pub queue_eof: AtomicBool,
    /// Set true if BGZF decompression of a raw block fails. Consumer
    /// surfaces this as `Err` rather than the silent `Ok(false)` that
    /// an empty queue + `queue_eof` would look like.
    ///
    /// **Atomic ordering:** identical to `queue_eof`. Producer
    /// stores while holding `decompressed`; consumer reads under
    /// the same lock.
    pub decomp_error: AtomicBool,

    // ── Block-parallel decompression state (file_granularity == false) ───────
    //
    // These three fields are inert in the inline (file-granularity) path. In
    // the block-parallel path multiple workers decompress one file's blocks
    // concurrently: the READ is serialized under `reader` (sequence-tagged via
    // `SortMergeReader::next_seq`), but the decompression happens outside the
    // lock, so results complete out of order and are reassembled here.
    /// Number of raw blocks that have been READ (under the reader lock) but not
    /// yet inserted into `reorder`. Incremented under the reader lock when a
    /// batch is read; decremented after the worker inserts that batch into
    /// `reorder`. The slot may declare EOF only once this reaches zero — a
    /// worker that observes reader-EOF must not truncate the merge while another
    /// worker still holds an in-flight (read-but-undelivered) block.
    pub in_flight: AtomicUsize,
    /// Set true once a read returns fewer raw blocks than requested, i.e. the
    /// disk reader reached a clean EOF. Distinct from `queue_eof`: `reader_eof`
    /// means "no more blocks will be read", whereas `queue_eof` means "every
    /// block has been read, decompressed, reassembled, and delivered to the
    /// FIFO". `queue_eof` is set only when `reader_eof && in_flight == 0 &&
    /// reorder.is_empty()`. Stored under the reader lock (Release), read in the
    /// finalize path (Acquire); being an atomic it does not participate in lock
    /// ordering.
    pub reader_eof: AtomicBool,
    /// Per-slot reorder buffer that reassembles out-of-order decompression
    /// results back into read (sequence) order before they are drained into the
    /// FIFO. Lock order: acquire `reorder` BEFORE `decompressed` (never the
    /// reverse); the reader lock, when held, is outermost. The consumer never
    /// touches this — it only pops the in-order FIFO.
    pub reorder: Mutex<ReorderBuffer<Vec<u8>>>,
}

impl SortMergeSlot {
    /// Construct an empty slot for `file_id` backed by `reader` (positioned
    /// past any codec file-magic) with the detected `codec`.
    #[must_use]
    pub fn new(file_id: u32, reader: BufReader<File>, codec: SpillCodec) -> Self {
        Self {
            file_id,
            codec,
            reader: Mutex::new(SortMergeReader { inner: reader, next_seq: 0 }),
            decompressed: Mutex::new(VecDeque::with_capacity(PHASE2_DECOMP_CAP)),
            queue_eof: AtomicBool::new(false),
            decomp_error: AtomicBool::new(false),
            in_flight: AtomicUsize::new(0),
            reader_eof: AtomicBool::new(false),
            reorder: Mutex::new(ReorderBuffer::new()),
        }
    }

    /// Returns `true` when this slot has cleanly delivered all of its
    /// output: `queue_eof` is set, the decompressed queue is empty, and
    /// no decompression error was recorded.
    ///
    /// A slot whose `decomp_error` flag is set is **never** reported as
    /// drained — an errored slot must surface as an error to its
    /// consumer, not be mistaken for clean EOF. Callers using
    /// `is_drained()` as a completion predicate must check
    /// [`Self::has_error`] (or the source-level error flags) to
    /// distinguish "still producing" from "failed".
    ///
    /// Both `queue_eof` and `decomp_error` are read while holding the
    /// `decompressed` mutex, honoring the file's read-under-lock
    /// discipline so the release-acquire chain establishes
    /// happens-before for both flags simultaneously.
    ///
    /// # Panics
    ///
    /// Panics if `decompressed` mutex is poisoned.
    #[must_use]
    pub fn is_drained(&self) -> bool {
        // Fast path: not yet at EOF — no need to take the lock.
        if !self.queue_eof.load(Ordering::Acquire) {
            return false;
        }
        let guard = self.decompressed.lock().expect("SortMergeSlot decompressed mutex poisoned");
        // Read `decomp_error` under the same lock the producer held when
        // storing it. An errored slot is not a clean drain.
        if self.decomp_error.load(Ordering::Acquire) {
            return false;
        }
        guard.is_empty()
    }

    /// Returns `true` if this slot recorded a decompression error.
    ///
    /// Read under the `decompressed` mutex to honor the file's
    /// read-under-lock discipline. Use alongside [`Self::is_drained`] to
    /// distinguish a clean EOF (`is_drained() == true`) from a failed
    /// slot (`has_error() == true`), since `is_drained()` returns
    /// `false` in both the "still producing" and "errored" cases.
    ///
    /// # Panics
    ///
    /// Panics if `decompressed` mutex is poisoned.
    #[must_use]
    pub fn has_error(&self) -> bool {
        let _guard = self.decompressed.lock().expect("SortMergeSlot decompressed mutex poisoned");
        self.decomp_error.load(Ordering::Acquire)
    }

    /// Gather probe statistics for this slot: `(pending_blocks,
    /// pending_bytes, active)`.
    ///
    /// `pending_blocks` is the count of decompressed blocks waiting
    /// for the consumer. `pending_bytes` is the sum of their byte
    /// lengths. `active` is `!queue_eof` (the slot is still being
    /// fed by the producer).
    ///
    /// # Panics
    ///
    /// Panics if `decompressed` mutex is poisoned.
    #[must_use]
    pub fn probe_stats(&self) -> (u64, u64, bool) {
        let dec = self.decompressed.lock().expect("SortMergeSlot decompressed mutex poisoned");
        #[allow(clippy::cast_possible_truncation)]
        let pending_blocks = dec.len() as u64;
        let pending_bytes: u64 = dec.iter().map(|buf| buf.len() as u64).sum();
        drop(dec);
        // `Acquire` for uniformity with every other `queue_eof` reader, even
        // though this is a best-effort diagnostics probe.
        let active = !self.queue_eof.load(Ordering::Acquire);
        (pending_blocks, pending_bytes, active)
    }

    /// Current number of decompressed blocks resident in the FIFO. Locks
    /// `decompressed` for an O(1) `len()` read — used by the decompressor's
    /// emptiest-first refill order (most-starved slot first). Cheaper than
    /// [`Self::probe_stats`], which also sums per-block byte lengths.
    ///
    /// # Panics
    ///
    /// Panics if the `decompressed` mutex is poisoned.
    #[must_use]
    pub fn fifo_len(&self) -> usize {
        self.decompressed.lock().expect("SortMergeSlot decompressed mutex poisoned").len()
    }

    // ── Block-parallel decompression helpers (file_granularity == false) ─────

    /// Block-parallel admission control: may a worker read another batch of raw
    /// blocks (whose first block would be tagged `next_seq`) into this slot's
    /// reorder window?
    ///
    /// Combines [`ReorderBuffer::would_accept`] (the deadlock-free predicate the
    /// pipeline uses elsewhere) with a hard `heap_bytes < window_budget`
    /// backstop. The backstop is what actually bounds memory: `would_accept`
    /// alone returns `true` (accept-all) while the buffer is *stuck* (the front
    /// sequence not yet decompressed), which would let a slow straggler balloon
    /// the window. The backstop is deadlock-safe **in this topology** because
    /// reads are serialized and densely sequenced, so the front gap is ALWAYS an
    /// already-read in-flight block (some worker is decompressing it) — never an
    /// unread block that a new read would be required to fetch. Refusing new
    /// reads therefore cannot wedge progress.
    ///
    /// `window_budget == 0` means unlimited.
    ///
    /// # Panics
    ///
    /// Panics if the `reorder` mutex is poisoned.
    #[must_use]
    pub fn bp_reorder_admits(&self, next_seq: u64, window_budget: u64) -> bool {
        let rb = self.reorder.lock().expect("reorder mutex poisoned");
        if !rb.would_accept(next_seq, window_budget) {
            return false;
        }
        window_budget == 0 || rb.heap_bytes() < window_budget
    }

    /// Reserve `count` in-flight blocks just read under the reader lock.
    ///
    /// **Ordering requirement:** the caller MUST call this BEFORE
    /// [`Self::bp_set_reader_eof`] for the EOF-carrying batch. Both stores are
    /// `Release`; the lock-free finalizer in [`Self::drain_locked_and_finalize`]
    /// reads `reader_eof` (Acquire) then `in_flight` (Acquire) WITHOUT the
    /// reader lock, so only this publish order guarantees that observing
    /// `reader_eof == true` also makes this batch's `in_flight` increment
    /// visible — otherwise the finalizer can declare a clean EOF that truncates
    /// the EOF read's own still-in-flight block (loom-verified; see
    /// `tests/loom_merge_slots.rs`).
    pub fn bp_add_in_flight(&self, count: usize) {
        self.in_flight.fetch_add(count, Ordering::Release);
    }

    /// Mark the disk reader as having reached a clean EOF (called under the
    /// reader lock when a read returns fewer blocks than requested).
    ///
    /// **Ordering requirement:** call this AFTER [`Self::bp_add_in_flight`] has
    /// reserved the current batch — see that method's note. Publishing
    /// `reader_eof` before the in-flight reservation reopens a truncation race.
    pub fn bp_set_reader_eof(&self) {
        self.reader_eof.store(true, Ordering::Release);
    }

    /// Publish the accounting for a batch of `count` raw blocks just read under
    /// the reader lock, in the one order that is correct: reserve the in-flight
    /// blocks FIRST, then (on a short read) set `reader_eof`.
    ///
    /// This is the single source of truth for the publish order — both the
    /// production worker (`SortSpillDecompress::try_fill_block_parallel_slot`)
    /// and the loom model (`tests/loom_merge_slots.rs`) call it, so the ordering
    /// is model-checked against the real code and the two cannot drift.
    ///
    /// **Why this order (loom-verified).** The lock-free finalizer in
    /// [`Self::drain_locked_and_finalize`] reads `reader_eof` (Acquire) then
    /// `in_flight` (Acquire) WITHOUT the reader lock, so the reader lock does not
    /// order this batch's accounting against it. Both stores are `Release`; a
    /// finalizer that observes `reader_eof == true` therefore also observes this
    /// (possibly EOF-carrying) batch's `in_flight` increment, and so cannot
    /// finalize `queue_eof` while the batch is still in flight. Setting
    /// `reader_eof` first would let a finalizer see `reader_eof == true` with
    /// `in_flight == 0` after earlier blocks drained, finalizing a clean EOF that
    /// silently truncates this batch's own block. Swapping the two lines makes
    /// `tests/loom_merge_slots.rs` fail (as it did before the fix in `ac0a2ad`).
    pub fn bp_commit_read(&self, count: usize, hit_eof: bool) {
        self.bp_add_in_flight(count);
        if hit_eof {
            self.bp_set_reader_eof();
        }
    }

    /// Number of additional decompressed blocks the FIFO can accept before it
    /// hits [`PHASE2_DECOMP_CAP`].
    ///
    /// # Panics
    ///
    /// Panics if the `decompressed` mutex is poisoned.
    #[must_use]
    pub fn bp_fifo_room(&self) -> usize {
        let dec = self.decompressed.lock().expect("decompressed mutex poisoned");
        PHASE2_DECOMP_CAP.saturating_sub(dec.len())
    }

    /// Tracked reorder-window heap bytes (for tests / diagnostics).
    ///
    /// # Panics
    ///
    /// Panics if the `reorder` mutex is poisoned.
    #[must_use]
    pub fn bp_reorder_heap_bytes(&self) -> u64 {
        self.reorder.lock().expect("reorder mutex poisoned").heap_bytes()
    }

    /// Insert a freshly-decompressed batch `[start_seq, start_seq + count)` into
    /// the reorder buffer, release the in-flight reservation, drain any now-ready
    /// (in-order) blocks into the FIFO (bounded by [`PHASE2_DECOMP_CAP`]), and
    /// finalize `queue_eof` if the slot is fully delivered.
    ///
    /// `blocks.len()` need not equal `count`: an empty final read (`count == 0`)
    /// still calls through here so the EOF can be finalized once the last
    /// in-flight block drains. Returns `true` if it made progress (drained at
    /// least one block or finalized EOF).
    ///
    /// # Panics
    ///
    /// Panics if the `reorder` or `decompressed` mutex is poisoned.
    pub fn bp_insert_drain_finalize(
        &self,
        start_seq: u64,
        blocks: Vec<Vec<u8>>,
        count: usize,
    ) -> bool {
        let mut rb = self.reorder.lock().expect("reorder mutex poisoned");
        for (i, block) in blocks.into_iter().enumerate() {
            let size = block.len();
            // `start_seq + i` cannot exceed the number of blocks read from one
            // spill file, which is far below u64::MAX.
            rb.insert_with_size(start_seq + i as u64, block, size);
        }
        // The reservation now lives in `reorder`, so release it. AcqRel so the
        // finalize load below sees a consistent count.
        self.in_flight.fetch_sub(count, Ordering::AcqRel);
        self.drain_locked_and_finalize(&mut rb)
    }

    /// Drain-only block-parallel pass: move any in-order ready blocks from the
    /// reorder buffer into the FIFO (used when the FIFO had no room earlier, or
    /// post-reader-EOF to flush stragglers other workers inserted) and finalize
    /// `queue_eof`. Returns `true` if it made progress.
    ///
    /// # Panics
    ///
    /// Panics if the `reorder` or `decompressed` mutex is poisoned.
    pub fn bp_drain_and_finalize(&self) -> bool {
        let mut rb = self.reorder.lock().expect("reorder mutex poisoned");
        self.drain_locked_and_finalize(&mut rb)
    }

    /// Shared drain + finalize, called with the `reorder` lock held. Acquires
    /// `decompressed` (lock order: reorder → decompressed) so `queue_eof` is
    /// stored under the same mutex the consumer reads it under.
    fn drain_locked_and_finalize(&self, rb: &mut ReorderBuffer<Vec<u8>>) -> bool {
        let mut dec = self.decompressed.lock().expect("decompressed mutex poisoned");
        let mut room = PHASE2_DECOMP_CAP.saturating_sub(dec.len());
        let mut drained = 0usize;
        while room > 0 {
            let Some(block) = rb.try_pop_next() else { break };
            dec.push_back(block);
            room -= 1;
            drained += 1;
        }
        // Finalize EOF: reader is exhausted, every read block has been inserted
        // (in_flight == 0), and the reorder buffer is fully drained. Stored
        // under `decompressed` so the consumer's release-acquire chain sees it.
        let mut finalized = false;
        if !self.queue_eof.load(Ordering::Acquire)
            && self.reader_eof.load(Ordering::Acquire)
            && self.in_flight.load(Ordering::Acquire) == 0
            && rb.is_empty()
        {
            self.queue_eof.store(true, Ordering::Release);
            finalized = true;
        }
        drained > 0 || finalized
    }
}

// These exercise `SortMergeSlot` outside `loom::model`, which is illegal once
// the primitives are loom's, so they compile only in a normal (non-loom) build.
// The `--cfg loom` invocation runs `tests/loom_merge_slots.rs` exclusively.
#[cfg(all(test, not(loom)))]
mod tests {
    use super::*;

    fn empty_reader() -> BufReader<File> {
        BufReader::new(tempfile::tempfile().expect("create tempfile"))
    }

    #[test]
    fn new_slot_starts_empty_and_not_drained() {
        let slot = SortMergeSlot::new(0, empty_reader(), SpillCodec::Bgzf);
        assert_eq!(slot.file_id, 0);
        assert!(!slot.is_drained(), "fresh slot not drained until queue_eof");
        let (pending_blocks, pending_bytes, active) = slot.probe_stats();
        assert_eq!(pending_blocks, 0);
        assert_eq!(pending_bytes, 0);
        assert!(active, "active until queue_eof flips");
    }

    #[test]
    fn drained_when_queue_eof_and_empty() {
        let slot = SortMergeSlot::new(0, empty_reader(), SpillCodec::Bgzf);

        // Push a decompressed block; not drained even if queue_eof.
        slot.decompressed.lock().unwrap().push_back(vec![0xAB, 0xCD]);
        slot.queue_eof.store(true, Ordering::Release);
        assert!(!slot.is_drained(), "not drained while decompressed non-empty");

        // Consume the block.
        let popped = slot.decompressed.lock().unwrap().pop_front();
        assert_eq!(popped, Some(vec![0xAB, 0xCD]));
        assert!(slot.is_drained(), "drained once queue empty AND queue_eof");
    }

    #[test]
    fn fifo_len_reports_block_count() {
        let slot = SortMergeSlot::new(0, empty_reader(), SpillCodec::Bgzf);
        assert_eq!(slot.fifo_len(), 0);
        slot.decompressed.lock().unwrap().push_back(vec![0u8; 10]);
        slot.decompressed.lock().unwrap().push_back(vec![0u8; 20]);
        assert_eq!(slot.fifo_len(), 2, "counts blocks, not bytes");
    }

    #[test]
    fn probe_stats_counts_decompressed_only() {
        let slot = SortMergeSlot::new(0, empty_reader(), SpillCodec::Bgzf);
        slot.decompressed.lock().unwrap().push_back(vec![0u8; 1024]);
        slot.decompressed.lock().unwrap().push_back(vec![0u8; 512]);

        let (pending_blocks, pending_bytes, active) = slot.probe_stats();
        assert_eq!(pending_blocks, 2);
        assert_eq!(pending_bytes, 1024 + 512);
        assert!(active, "not yet queue_eof");
    }

    #[test]
    fn fifo_order() {
        let slot = SortMergeSlot::new(7, empty_reader(), SpillCodec::Bgzf);
        // Push in order; pop must yield in same order.
        slot.decompressed.lock().unwrap().push_back(vec![0]);
        slot.decompressed.lock().unwrap().push_back(vec![1]);
        slot.decompressed.lock().unwrap().push_back(vec![2]);
        slot.decompressed.lock().unwrap().push_back(vec![3]);

        let mut popped = Vec::new();
        while let Some(b) = slot.decompressed.lock().unwrap().pop_front() {
            popped.push(b[0]);
        }
        assert_eq!(popped, vec![0, 1, 2, 3]);
    }

    #[test]
    fn decomp_error_flag_set_and_read() {
        let slot = SortMergeSlot::new(0, empty_reader(), SpillCodec::Bgzf);
        assert!(!slot.decomp_error.load(Ordering::Acquire));
        slot.decomp_error.store(true, Ordering::Release);
        assert!(slot.decomp_error.load(Ordering::Acquire));
    }

    #[test]
    fn errored_slot_is_not_drained() {
        let slot = SortMergeSlot::new(0, empty_reader(), SpillCodec::Bgzf);

        // Mark EOF with an empty queue but a decompression error set: this
        // must NOT be reported as a clean drain, otherwise a caller using
        // is_drained() as a completion predicate would treat the failure as
        // successful EOF and silently truncate the merge.
        slot.queue_eof.store(true, Ordering::Release);
        slot.decomp_error.store(true, Ordering::Release);
        assert!(!slot.is_drained(), "errored slot must not be reported as drained");
        assert!(slot.has_error(), "has_error must surface the recorded decomp error");
    }

    #[test]
    fn clean_eof_reports_no_error() {
        let slot = SortMergeSlot::new(0, empty_reader(), SpillCodec::Bgzf);
        slot.queue_eof.store(true, Ordering::Release);
        assert!(slot.is_drained(), "clean empty + queue_eof is drained");
        assert!(!slot.has_error(), "clean drain has no error");
    }

    // ── Block-parallel helper tests ─────────────────────────────────────────

    /// A single in-order batch drains straight into the FIFO and (because the
    /// reader is at EOF with no other in-flight blocks) finalizes `queue_eof`.
    #[test]
    fn bp_single_batch_drains_and_finalizes() {
        let slot = SortMergeSlot::new(0, empty_reader(), SpillCodec::Bgzf);
        slot.bp_add_in_flight(2);
        slot.bp_set_reader_eof();
        let progressed = slot.bp_insert_drain_finalize(0, vec![vec![1u8], vec![2u8]], 2);
        assert!(progressed);
        assert!(slot.queue_eof.load(Ordering::Acquire), "reader_eof + drained ⇒ queue_eof");
        let mut dec = slot.decompressed.lock().unwrap();
        assert_eq!(dec.pop_front(), Some(vec![1u8]));
        assert_eq!(dec.pop_front(), Some(vec![2u8]));
    }

    /// EOF-with-stragglers: a worker reads the FINAL (short) batch and hits
    /// reader-EOF while another worker still holds an earlier in-flight batch.
    /// The EOF-observing worker must NOT finalize `queue_eof` (which would
    /// truncate the merge); only after the straggler is delivered, in order,
    /// does the slot finalize. No block is dropped or reordered.
    #[test]
    fn bp_eof_with_straggler_does_not_truncate() {
        let slot = SortMergeSlot::new(0, empty_reader(), SpillCodec::Bgzf);

        // Worker A reserved seqs 0,1; worker B reserved seqs 2,3 on the final
        // (short) read and set reader_eof.
        slot.bp_add_in_flight(2); // A: seq 0,1
        slot.bp_add_in_flight(2); // B: seq 2,3 (final batch)
        slot.bp_set_reader_eof();

        // Worker B finishes decompressing FIRST and delivers seqs 2,3. Gap at
        // 0,1 ⇒ nothing drains, and in_flight (A's 2) is non-zero ⇒ no EOF.
        let b_progress = slot.bp_insert_drain_finalize(2, vec![vec![2u8], vec![3u8]], 2);
        assert!(!b_progress, "straggler ahead-of-gap insert drains nothing and can't finalize");
        assert!(!slot.queue_eof.load(Ordering::Acquire), "must NOT declare EOF with A in flight");
        assert!(slot.decompressed.lock().unwrap().is_empty(), "nothing delivered yet");

        // Worker A finishes and delivers seqs 0,1 ⇒ all four drain in order and
        // EOF finalizes.
        let a_progress = slot.bp_insert_drain_finalize(0, vec![vec![0u8], vec![1u8]], 2);
        assert!(a_progress);
        assert!(slot.queue_eof.load(Ordering::Acquire), "EOF finalizes once straggler delivered");

        let mut dec = slot.decompressed.lock().unwrap();
        let order: Vec<u8> = std::iter::from_fn(|| dec.pop_front()).map(|b| b[0]).collect();
        assert_eq!(order, vec![0, 1, 2, 3], "blocks delivered in read order, none truncated");
    }

    /// The reorder window stays bounded when the front sequence straggles: a
    /// worker keeps decompressing ahead-of-gap blocks, but `bp_reorder_admits`
    /// refuses new reads once `heap_bytes` reaches the window budget, so the
    /// buffer never balloons past `budget + one batch`.
    #[test]
    fn bp_reorder_window_is_bounded_under_straggler() {
        const BLOCK: usize = 1024;
        const BUDGET: u64 = 4 * BLOCK as u64; // 4 blocks
        let slot = SortMergeSlot::new(0, empty_reader(), SpillCodec::Bgzf);

        // Seq 0 is the straggler — it is reserved but never decompressed/inserted,
        // so the buffer can never pop and keeps a permanent front gap.
        slot.bp_add_in_flight(1); // seq 0 in flight forever (the straggler)

        // Workers race ahead delivering seqs 1,2,3,… as long as admission allows.
        let mut next = 1u64;
        let mut admitted = 0;
        while slot.bp_reorder_admits(next, BUDGET) {
            slot.bp_add_in_flight(1);
            slot.bp_insert_drain_finalize(next, vec![vec![0u8; BLOCK]], 1);
            // Front gap at seq 0 ⇒ nothing drains.
            assert!(slot.decompressed.lock().unwrap().is_empty());
            assert!(!slot.queue_eof.load(Ordering::Acquire), "never EOF while seq 0 in flight");
            assert!(
                slot.bp_reorder_heap_bytes() <= BUDGET,
                "reorder window must stay within budget, got {} > {BUDGET}",
                slot.bp_reorder_heap_bytes(),
            );
            next += 1;
            admitted += 1;
            assert!(admitted < 1000, "admission must eventually backpressure, not loop forever");
        }
        assert!(admitted > 0, "should admit at least some ahead-of-gap blocks");
        assert!(
            slot.bp_reorder_heap_bytes() >= BUDGET.saturating_sub(BLOCK as u64),
            "should have filled the window before backpressuring",
        );
    }
}

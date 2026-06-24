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
//! Standalone file-to-file `fgumi sort` uses `worker_pool::Phase2FileState`
//! plus the gap-filler this module dropped. That path is NOT broken: its merge
//! consumer is a parked OS thread (no framework Skip), so the v4 deadlock
//! cannot occur there, and its single-reader/**multi-decompressor** topology
//! genuinely needs the reorder buffer + in-flight counter (a plain FIFO
//! suffices HERE only because read-and-decompress is one inline op, so blocks
//! decompress strictly in read order). Unifying the two drivers is a deferred
//! goal — see commit `9d6d7e9` / PR #395 and
//! `docs/design/sort-phase2-unification-deferral.md` — gated on this streaming
//! path gaining file-to-file drive + `--write-index` + `--verify` + a
//! deadlock-safe consumer.
// TODO(#395): unify the two Phase-2 sort drivers — see docs/design/sort-phase2-unification-deferral.md

use std::collections::VecDeque;
use std::fs::File;
use std::io::BufReader;
use std::sync::Mutex;
use std::sync::atomic::{AtomicBool, Ordering};

use crate::codec::SpillCodec;

/// Per-slot decompressed-block queue cap. Bounds in-flight
/// decompressed memory: `num_slots × PHASE2_DECOMP_CAP ×` per-entry size,
/// where each entry is a decompressed BGZF block or zstd frame ranging from
/// ~64 KB (BGZF) up to 256 KB (zstd worst case). Matches
/// `worker_pool::PHASE2_DECOMP_CAP` numerically.
///
/// Hard cap — there is no admission escape. Producers skip a slot
/// when its queue is at this cap, returning to it on a later
/// `try_run` after the consumer has drained.
pub const PHASE2_DECOMP_CAP: usize = 8;

/// Disk reader state for a single spill file. Mutex'd separately
/// from `decompressed` so the producer can hold the reader during
/// disk I/O without blocking the consumer's pop.
pub struct SortMergeReader {
    /// Buffered file handle.
    pub inner: BufReader<File>,
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
}

impl SortMergeSlot {
    /// Construct an empty slot for `file_id` backed by `reader` (positioned
    /// past any codec file-magic) with the detected `codec`.
    #[must_use]
    pub fn new(file_id: u32, reader: BufReader<File>, codec: SpillCodec) -> Self {
        Self {
            file_id,
            codec,
            reader: Mutex::new(SortMergeReader { inner: reader }),
            decompressed: Mutex::new(VecDeque::with_capacity(PHASE2_DECOMP_CAP)),
            queue_eof: AtomicBool::new(false),
            decomp_error: AtomicBool::new(false),
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
}

#[cfg(test)]
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
}

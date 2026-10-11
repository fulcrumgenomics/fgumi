//! Per-spill-file shared state for the chain merge phase.
//!
//! `SortMergeSlot` is the per-spill-file shared state shuttled between
//! `SortSpillDecompress` (producer side: reads the spill file
//! positionally, parses and decompresses its frames, pushes to the slot's
//! bounded queue) and
//! `SortMerge` (consumer side: pops decompressed blocks from the
//! queue, parses records, drives the k-way merge) via
//! `Arc<SortMergeSlot>` clones.
//!
//! # Per-slot bounded queue design (v4 — see commit `9c39dea` / PR #389)
//!
//! Each slot carries a bounded queue of decompressed BGZF blocks
//! (`PHASE2_DECOMP_CAP` entries). Backpressure lives here — the
//! producer is **non-blocking**: pushes only when the queue has
//! space; otherwise skips this slot and tries the next. The consumer
//! is also **non-blocking**: when the queue is empty but the slot is
//! not yet `queue_eof` it reports `WouldBlock` (see
//! `external.rs::slot_try_load_block`) and the cooperative `SortMerge`
//! step yields: it declares the slot awaited
//! ([`crate::MergeDemand::await_slot`]) and its driver thread parks until
//! the producer's next delivery, EOF or failure on that slot unparks it
//! ([`crate::MergeDemand::notify_delivered`]). No condvar: the producer
//! never blocks, and only a delivery to the awaited slot wakes the merge.
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
//! ## Reads, the raw stash, and who owns a raw block
//!
//! The slot carries its spill file as a positional handle (`source`,
//! `body_start`, `len`) and no read buffer. Bytes arrive as leased read slices;
//! [`SortMergeSlot::bp_ingest_slice`] reorders them by their per-slot slice
//! sequence, parses every frame of each in-order slice under the `reader` lock
//! (a [`SpillFrameParser`] carries a frame that straddles two slices) and
//! appends the frames to `raw_stash`, sequence-stamped. A raw block is owned by
//! **the slot** from its parse ([`SortMergeSlot::bp_commit_read`] reserves it
//! as in flight) until its decompressed bytes are inserted
//! ([`SortMergeSlot::bp_insert_drain_finalize`] releases the reservation): any
//! worker — or the merge itself — may claim it in between. `queue_eof`
//! finalizes only when the reader is at EOF, nothing is in flight and the
//! reorder buffer is drained, so a stashed block can never be truncated away.
//!
//! The reorder buffer ([`SortMergeSlot::reorder`]) reassembles decompression
//! results that complete out of order, and [`SortMergeSlot::in_flight`] counts
//! parsed-but-undelivered blocks. Lock order, outermost first: `reader` →
//! `raw_stash` → `reorder` → `decompressed`; the consumer takes only
//! `raw_stash` (a claim) and `decompressed` (a pop), never `reader`. The
//! lock-free mirrors (`stash_len`, `stash_bytes`, `fifo_len_mirror`,
//! `reorder_len_mirror`, `issued_bytes`) are stored under the lock they mirror
//! (or by the single planner) so schedulers can read them without locking.
//!
//! ## Who drives it, and the OTHER Phase-2 implementation
//!
//! `SortMergeSlot` is the Phase-2 of the chain sort — standalone `fgumi sort`
//! and the fused `runall` sort. The spill writer opens each spill file as a
//! slot (`external.rs::open_spill_slot`), `fgumi sort`'s `SortSpillDecompress`
//! fills it (the `bp_*` methods), and `SortMerge` merges the set through
//! `MergeDriver::from_slots`.
//!
//! `worker_pool::Phase2FileState` is the other Phase-2: `RawExternalSorter`'s
//! `sort_records` drives it (`fgumi simulate`, and the `#[cfg(test)]` parity
//! oracle). `fgumi merge` uses neither: `RawExternalSorter::merge_bams` merges
//! its sorted inputs with `run_merge_loop`, which has no spill files and no
//! Phase-2 slots. `Phase2FileState` keeps
//! its own reorder buffer and in-flight counter because its
//! single-reader/**multi-decompressor** topology needs them, and it retains
//! the gap-filler this module dropped. Which implementation each command runs
//! is stated here, on [`PHASE2_DECOMP_CAP`], and on
//! `worker_pool::Phase2FileState`; keep the three in sync. (History: commit
//! `9d6d7e9` / PR #395.)

use std::collections::VecDeque;
use std::fs::File;
use std::io;
use std::sync::Arc;

// Concurrency primitives are sourced from `loom` under `--cfg loom` so the
// model-checking test (`tests/loom_merge_slots.rs`) exercises the REAL
// `SortMergeSlot` atomics/mutexes — every interleaving and memory reordering of
// the block-parallel EOF/in-flight/finalize protocol — instead of a hand-copied
// re-implementation. Under a normal build these are the `std` types verbatim.
#[cfg(loom)]
use loom::sync::Mutex;
#[cfg(loom)]
use loom::sync::atomic::{AtomicBool, AtomicU8, AtomicU32, AtomicU64, AtomicUsize, Ordering};
#[cfg(not(loom))]
use std::sync::Mutex;
#[cfg(not(loom))]
use std::sync::atomic::{AtomicBool, AtomicU8, AtomicU32, AtomicU64, AtomicUsize, Ordering};

use fgumi_bam_io::pread::{RawFrame, SliceLease};
use fgumi_bam_io::reorder::ReorderBuffer;

use crate::codec::SpillCodec;
use crate::spill_block_reader::{RawBlock, SpillFrameParser};

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
///
/// # Not the only bound of that name
///
/// `worker_pool.rs` declares a *different* constant of the same name
/// (`pub(crate) const PHASE2_DECOMP_CAP: usize = 8`) for its own Phase-2, and
/// that is the one governing `RawExternalSorter::sort_records`' merge
/// (`fgumi simulate`, and the test parity oracle). This constant bounds the
/// chain sort's slots (`fgumi sort`, the `runall` sort). So
/// `fgumi_sort::PHASE2_DECOMP_CAP` resolves to 32 while `fgumi simulate`'s
/// sort runs bounded at 8 (`fgumi merge`'s `merge_bams` has no Phase-2 slots,
/// so neither bounds it). The two are deliberately unequal (see this constant's history
/// above); the confusable part is only which command each one bounds.
pub const PHASE2_DECOMP_CAP: usize = 32;

/// Parse state for a single spill file: the slices that arrived out of order,
/// the frame parser, and the raw-block sequence. Mutex'd separately from
/// `decompressed` so parsing never blocks the consumer's pop. The slot holds no
/// read buffer of its own: the bytes arrive as leased read slices.
pub struct SortMergeReader {
    /// Slices that landed ahead of their predecessors, keyed by their per-slot
    /// slice sequence (the read edge is unordered), with each slice's `last`
    /// flag.
    pub slices: ReorderBuffer<(SliceLease, bool)>,
    /// Frame parser over the in-order slices (carries a frame that straddles
    /// two slices).
    pub parser: SpillFrameParser,
    /// Next per-slot sequence number to stamp on a parsed raw block: dense and
    /// in file order, because parsing is serialized under this lock. The
    /// slot's [`SortMergeSlot::reorder`] buffer reassembles the out-of-order
    /// decompression results back into this order.
    pub next_seq: u64,
    /// Heap bytes of the parser's carried partial frame, as last charged to
    /// [`SortMergeSlot::stash_bytes`].
    carry_charge: u64,
    /// File cursor of the positional window reader (`SortSpillDecompress`'s
    /// read-on-demand path).
    pub next_offset: u64,
    /// Frames the window reader parsed past what its caller asked for, handed
    /// out first on its next read.
    pub pending: VecDeque<RawFrame>,
}

impl SortMergeReader {
    fn new(codec: SpillCodec, body_start: u64) -> Self {
        Self {
            slices: ReorderBuffer::new(),
            parser: SpillFrameParser::new(codec),
            next_seq: 0,
            carry_charge: 0,
            next_offset: body_start,
            pending: VecDeque::new(),
        }
    }
}

/// What one [`SortMergeSlot::bp_ingest_slice`] did.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IngestOutcome {
    /// Frames parsed into the raw stash (from every slice that became in order).
    pub frames: usize,
    /// Net change in the slot's charged stash bytes ([`SortMergeSlot::stash_bytes`]):
    /// the resident slices the new frames pin, their owned allocations, and
    /// the change in the parser's carried partial frame.
    pub stashed_bytes: i64,
    /// The file's last slice was parsed with no partial frame left.
    pub hit_eof: bool,
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
    /// slot is opened (`open_spill_slot`). The `SortSpillDecompress` step
    /// reads it to decompress BGZF blocks or zstd frames.
    pub codec: SpillCodec,
    /// Positional handle on the spill file (no buffer: reads are positional
    /// slices into pooled buffers).
    source: Arc<File>,
    /// Offset of the first frame (past any codec file magic).
    body_start: u64,
    /// File length in bytes.
    len: u64,
    /// Parse state. Lock order (outermost first): `reader` → `raw_stash` →
    /// `reorder` → `decompressed`; the consumer takes only `raw_stash` and
    /// `decompressed`.
    pub reader: Mutex<SortMergeReader>,
    /// Parsed, not yet claimed raw frames, in sequence order. A raw block is
    /// owned by the slot from its parse (`bp_commit_read`) until its
    /// decompressed bytes are inserted (`bp_insert_drain_finalize`).
    pub raw_stash: Mutex<VecDeque<RawBlock>>,
    /// Lock-free mirror: resident bytes the stash holds — every read slice with
    /// a frame still in `raw_stash` (a borrowed frame pins its whole slice),
    /// every owned frame's allocation, and the parser's carried partial frame.
    /// Charged as frames are stashed and released as they are claimed (the
    /// slice with its last stashed frame; see [`RawBlock::charge`]), so it is
    /// what the slot holds, not the frames' lengths.
    pub stash_bytes: AtomicU64,
    /// Lock-free mirror: frames in `raw_stash`.
    pub stash_len: AtomicU32,
    /// Bytes requested from the file but not yet parsed (the read planner adds
    /// on issue, [`Self::bp_ingest_slice`] subtracts on parse).
    pub issued_bytes: AtomicU64,
    /// Lock-free mirror of `decompressed.len()`, stored under that lock.
    pub fifo_len_mirror: AtomicU32,
    /// Lock-free mirror of the blocks waiting in `reorder`, stored under that
    /// lock.
    pub reorder_len_mirror: AtomicU32,
    /// Soft cap on the decompressed FIFO, at most [`PHASE2_DECOMP_CAP`].
    pub fifo_cap: AtomicU32,
    /// The slot's read-ahead class (diagnostics).
    pub class: AtomicU8,
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
    /// Construct an empty slot for `file_id` over `source`, whose frames span
    /// `[body_start, len)`, with the detected `codec`. Allocates no read
    /// buffer.
    #[must_use]
    pub fn new(
        file_id: u32,
        source: Arc<File>,
        body_start: u64,
        len: u64,
        codec: SpillCodec,
    ) -> Self {
        Self {
            file_id,
            codec,
            source,
            body_start,
            len,
            reader: Mutex::new(SortMergeReader::new(codec, body_start)),
            raw_stash: Mutex::new(VecDeque::new()),
            stash_bytes: AtomicU64::new(0),
            stash_len: AtomicU32::new(0),
            issued_bytes: AtomicU64::new(0),
            fifo_len_mirror: AtomicU32::new(0),
            reorder_len_mirror: AtomicU32::new(0),
            fifo_cap: AtomicU32::new(mirror(PHASE2_DECOMP_CAP)),
            class: AtomicU8::new(0),
            decompressed: Mutex::new(VecDeque::with_capacity(PHASE2_DECOMP_CAP)),
            queue_eof: AtomicBool::new(false),
            decomp_error: AtomicBool::new(false),
            in_flight: AtomicUsize::new(0),
            reader_eof: AtomicBool::new(false),
            reorder: Mutex::new(ReorderBuffer::new()),
        }
    }

    /// An empty slot over an empty anonymous file (test support).
    ///
    /// # Panics
    /// If the anonymous file cannot be created.
    #[cfg(any(test, loom, feature = "test-utils"))]
    #[must_use]
    pub fn for_test(file_id: u32, codec: SpillCodec) -> Self {
        Self::new(file_id, Arc::new(tempfile::tempfile().expect("tempfile")), 0, 0, codec)
    }

    /// The positional handle on the spill file.
    #[must_use]
    pub fn source(&self) -> &Arc<File> {
        &self.source
    }

    /// Offset of the first frame.
    #[must_use]
    pub fn body_start(&self) -> u64 {
        self.body_start
    }

    /// File length in bytes.
    #[must_use]
    pub fn len(&self) -> u64 {
        self.len
    }

    /// Whether the file holds no frame bytes.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len <= self.body_start
    }

    /// Fail the slot: set `decomp_error` and `queue_eof` under the
    /// `decompressed` mutex (the consumer reads both under it), so the merge
    /// surfaces the error in preference to a clean EOF.
    ///
    /// # Panics
    /// If the `decompressed` mutex is poisoned.
    pub fn mark_failed(&self) {
        let _g = lock_ranked(&self.decompressed, LockRank::Decompressed);
        self.decomp_error.store(true, Ordering::Release);
        self.queue_eof.store(true, Ordering::Release);
    }

    /// Pop the FIFO head from `dec` (this slot's locked `decompressed` queue)
    /// and store the lock-free length mirror under that lock — the one consumer
    /// pop, so the mirror the claim admission reads is the consumer's.
    pub fn pop_locked(&self, dec: &mut VecDeque<Vec<u8>>) -> Option<Vec<u8>> {
        let block = dec.pop_front()?;
        self.fifo_len_mirror.store(mirror(dec.len()), Ordering::Relaxed);
        Some(block)
    }

    /// Lock the FIFO and [`Self::pop_locked`] its head.
    ///
    /// # Panics
    /// If the `decompressed` mutex is poisoned.
    pub fn pop_decompressed(&self) -> Option<Vec<u8>> {
        let mut dec = lock_ranked(&self.decompressed, LockRank::Decompressed);
        self.pop_locked(&mut dec)
    }

    /// Decompressed blocks in the FIFO, read lock-free (exact as of the last
    /// push or pop).
    #[must_use]
    pub fn fifo_len_relaxed(&self) -> usize {
        self.fifo_len_mirror.load(Ordering::Relaxed) as usize
    }

    /// Decompressed blocks waiting in `reorder`, read lock-free.
    #[must_use]
    pub fn reorder_len_relaxed(&self) -> usize {
        self.reorder_len_mirror.load(Ordering::Relaxed) as usize
    }

    /// Record `bytes` requested from the file (the read planner's side of
    /// `issued_bytes`).
    pub fn bp_note_issued(&self, bytes: u64) {
        self.issued_bytes.fetch_add(bytes, Ordering::Relaxed);
    }

    /// Ingest one read slice of this file (slice sequence `seq`; `last` when it
    /// covers the file's final byte). Under `reader`: reorder the slice by
    /// `seq`, then for every slice now in order parse ALL its frames, publish
    /// them as in flight (`bp_commit_read`) and append them to `raw_stash`
    /// stamped with dense sequence numbers. When the last slice is parsed with
    /// no partial frame left the read side is at EOF, and the slot is drained
    /// and finalized here, so a last slice that carries no frame (only the BGZF
    /// EOF marker) still sets `queue_eof`.
    ///
    /// This crate has no handle on the merge's wake, so the **caller must wake
    /// a merge awaiting this slot** (`MergeDemand::notify_delivered`) when the
    /// outcome reports `hit_eof` (the slot may now be `queue_eof`) and when this
    /// returns `Err` (the slot is failed); otherwise a parked merge sleeps until
    /// its timer.
    ///
    /// # Errors
    /// A malformed frame, or a partial frame at the end of the last slice
    /// (`UnexpectedEof`); the slot is marked failed first.
    ///
    /// # Panics
    /// If a slot mutex is poisoned.
    pub fn bp_ingest_slice(
        &self,
        seq: u32,
        lease: SliceLease,
        last: bool,
    ) -> io::Result<IngestOutcome> {
        let mut out = IngestOutcome { frames: 0, stashed_bytes: 0, hit_eof: false };
        let mut r = lock_ranked(&self.reader, LockRank::Reader);
        let size = lease.len();
        r.slices.insert_with_size(u64::from(seq), (lease, last), size);
        let mut frames = Vec::new();
        while let Some((lease, last)) = r.slices.try_pop_next() {
            let parsed = r
                .parser
                .push(&lease, &mut frames)
                .and_then(|_| if last { r.parser.finish() } else { Ok(()) });
            if let Err(e) = parsed {
                drop(r);
                self.mark_failed();
                return Err(e);
            }
            let prev = self.issued_bytes.fetch_sub(lease.len() as u64, Ordering::Relaxed);
            debug_assert!(prev >= lease.len() as u64, "a slice landed that was never issued");
            // The slice stays resident while any frame borrowing it is stashed;
            // its last borrowed frame carries the slice's charge.
            let slice_charge = lease.capacity() as u64;
            let last_borrowed = frames
                .iter()
                .rposition(|f| matches!(f, fgumi_bam_io::pread::RawFrame::Borrowed { .. }));
            drop(lease);
            if last {
                out.hit_eof = true;
            }
            let carry = r.parser.carry_capacity() as u64;
            let carry_before = std::mem::replace(&mut r.carry_charge, carry);
            let mut stash = lock_ranked(&self.raw_stash, LockRank::RawStash);
            let n = frames.len();
            let mut charged = 0u64;
            for (i, frame) in frames.drain(..).enumerate() {
                let charge =
                    match &frame {
                        fgumi_bam_io::pread::RawFrame::Owned(v) => v.capacity() as u64,
                        fgumi_bam_io::pread::RawFrame::Borrowed { .. } => {
                            if Some(i) == last_borrowed { slice_charge } else { 0 }
                        }
                    };
                out.frames += 1;
                charged += charge;
                stash.push_back(RawBlock { seq: r.next_seq, frame, charge });
                r.next_seq += 1;
            }
            // Charge before publishing (a claim releases only what was charged)
            // and the carry's change with it: the carried frame moved into an
            // owned frame charged above, or a new tail was carried.
            let delta = i64::try_from(charged + carry).unwrap_or(i64::MAX)
                - i64::try_from(carry_before).unwrap_or(i64::MAX);
            if delta >= 0 {
                self.stash_bytes.fetch_add(delta.unsigned_abs(), Ordering::Relaxed);
            } else {
                self.stash_bytes.fetch_sub(delta.unsigned_abs(), Ordering::Relaxed);
            }
            out.stashed_bytes += delta;
            self.publish_stashed(stash.len(), n, last);
        }
        if out.hit_eof {
            self.bp_drain_and_finalize();
        }
        Ok(out)
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
    /// `decomp_error` is read while holding the `decompressed` mutex, honoring
    /// the file's read-under-lock discipline. `queue_eof` is read *before* the
    /// lock, on the fast path, and that is safe for the reverse reason: a
    /// producer only stores `queue_eof` while holding this mutex, so observing
    /// it `true` means the storing critical section has already released, and
    /// the lock taken immediately below then establishes happens-before for
    /// `decomp_error`. Observing it `false` unlocked can only be stale in the
    /// direction that returns `false` early, which is the answer a
    /// not-yet-finalized slot warrants anyway.
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
        let guard = lock_ranked(&self.decompressed, LockRank::Decompressed);
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
        let _guard = lock_ranked(&self.decompressed, LockRank::Decompressed);
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
        let dec = lock_ranked(&self.decompressed, LockRank::Decompressed);
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
        lock_ranked(&self.decompressed, LockRank::Decompressed).len()
    }

    /// The merge can make progress on this slot: a decompressed block is queued,
    /// or the slot reached EOF (clean or failed). Read under `decompressed`, so
    /// it orders against a producer's push (see `merge_demand`'s module doc).
    ///
    /// # Panics
    ///
    /// Panics if the `decompressed` mutex is poisoned.
    #[must_use]
    pub fn has_block_or_eof(&self) -> bool {
        let d = lock_ranked(&self.decompressed, LockRank::Decompressed);
        !d.is_empty() || self.queue_eof.load(Ordering::Acquire)
    }

    /// What this slot is doing while the merge stalls on it, from lock-free
    /// state only: the merge must not take a slot lock, which would make a
    /// decompress worker skip the very slot the merge is waiting on.
    /// `Decompressing` when some claimed block is in flight, else `Stashed`
    /// when frames wait in the stash, else `Issued` when reads are
    /// outstanding, else `Starved`. Only the block-parallel path raises these
    /// mirrors, so under
    /// `--sort::file-granularity` (inline decompress) every stall reports
    /// [`AwaitedSlotState::Starved`]; the `--sort-stats` awaited-slot line says
    /// the buckets are not classified there
    /// ([`crate::MergeDemandStats::mark_inline_decompress`]).
    ///
    /// [`AwaitedSlotState::Starved`]: crate::AwaitedSlotState::Starved
    #[must_use]
    pub fn awaited_state(&self) -> crate::AwaitedSlotState {
        use crate::AwaitedSlotState as S;
        // `in_flight` first: an ingest stores `stash_len` before raising
        // `in_flight` (`publish_stashed`), so this pair never reads a stashed
        // batch as decompressing.
        let in_flight = self.in_flight.load(Ordering::Acquire);
        let stashed = self.stash_len.load(Ordering::Acquire) as usize;
        if in_flight > stashed {
            S::Decompressing
        } else if stashed > 0 {
            S::Stashed
        } else if self.issued_bytes.load(Ordering::Relaxed) > 0 {
            S::Issued
        } else {
            S::Starved
        }
    }

    // ── Block-parallel decompression helpers (file_granularity == false) ─────

    /// Block-parallel admission control for the fused read path: may a worker
    /// read another batch of raw blocks (whose first block would be tagged
    /// `next_seq`) into this slot's reorder window? The same rule
    /// [`Self::bp_claim_raw`] applies to a claim (`window_admits`).
    ///
    /// # Panics
    ///
    /// Panics if the `reorder` mutex is poisoned.
    #[must_use]
    pub fn bp_reorder_admits(&self, next_seq: u64, window_budget: u64) -> bool {
        let rb = lock_ranked(&self.reorder, LockRank::Reorder);
        window_admits(&rb, next_seq, window_budget)
    }

    /// Reserve `count` in-flight blocks just read under the reader lock.
    ///
    /// **Ordering requirement:** the caller MUST call this BEFORE
    /// [`Self::bp_set_reader_eof`] for the EOF-carrying batch. Both stores are
    /// `Release`; the lock-free finalizer in `drain_locked_and_finalize`
    /// reads `reader_eof` (Acquire) then `in_flight` (Acquire) WITHOUT the
    /// reader lock, so only this publish order guarantees that observing
    /// `reader_eof == true` also makes this batch's `in_flight` increment
    /// visible — otherwise the finalizer can declare a clean EOF that truncates
    /// the EOF read's own still-in-flight block (loom-verified; see
    /// `tests/loom_merge_slots.rs`).
    pub(crate) fn bp_add_in_flight(&self, count: usize) {
        self.in_flight.fetch_add(count, Ordering::Release);
    }

    /// Mark the disk reader as having reached a clean EOF (called under the
    /// reader lock when a read returns fewer blocks than requested).
    ///
    /// **Ordering requirement:** call this AFTER [`Self::bp_add_in_flight`] has
    /// reserved the current batch — see that method's note. Publishing
    /// `reader_eof` before the in-flight reservation reopens a truncation race.
    pub(crate) fn bp_set_reader_eof(&self) {
        self.reader_eof.store(true, Ordering::Release);
    }

    /// Publish `count` frames just appended to the stash (now `stash_len`
    /// long), called with the `raw_stash` lock held. The stash mirror is stored
    /// before the in-flight reservation is raised (Release), so a reader that
    /// loads `in_flight` (Acquire) before `stash_len` — [`Self::awaited_state`] —
    /// never sees the batch in flight but not stashed. The reservation still
    /// exists before any claim can take one of the frames: claims need the
    /// `raw_stash` lock this caller holds.
    fn publish_stashed(&self, stash_len: usize, count: usize, hit_eof: bool) {
        self.stash_len.store(mirror(stash_len), Ordering::Relaxed);
        self.bp_commit_read(count, hit_eof);
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
    /// `drain_locked_and_finalize` reads `reader_eof` (Acquire) then
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
        let dec = lock_ranked(&self.decompressed, LockRank::Decompressed);
        PHASE2_DECOMP_CAP.saturating_sub(dec.len())
    }

    /// Tracked reorder-window heap bytes (for tests / diagnostics).
    ///
    /// # Panics
    ///
    /// Panics if the `reorder` mutex is poisoned.
    #[must_use]
    pub fn bp_reorder_heap_bytes(&self) -> u64 {
        lock_ranked(&self.reorder, LockRank::Reorder).heap_bytes()
    }

    /// Insert a freshly-decompressed batch `[start_seq, start_seq + count)` into
    /// the reorder buffer, release the in-flight reservation, drain any now-ready
    /// (in-order) blocks into the FIFO (bounded by [`PHASE2_DECOMP_CAP`]), and
    /// finalize `queue_eof` if the slot is fully delivered.
    ///
    /// `count` is the reservation being released and **must** equal
    /// `blocks.len()`. The one exception is the empty final read
    /// (`count == 0`, `blocks` empty), which still calls through here so the
    /// EOF can be finalized once the last in-flight block drains. Returns
    /// `true` if it made progress (drained at least one block or finalized
    /// EOF).
    ///
    /// Delivering fewer blocks than reserved is **silent truncation**, not a
    /// tolerated shape: the short batch releases the reservation for sequence
    /// numbers that were never inserted, so `in_flight` reaches zero and
    /// `reorder` reports empty with those sequences simply absent. The slot
    /// then finalizes a *clean* `queue_eof` — `is_drained()` true,
    /// `has_error()` false — and the consumer treats a truncated spill file as
    /// a complete one. Releasing *more* than was reserved is worse still: the
    /// `fetch_sub` below wraps `in_flight` to `usize::MAX`, after which the
    /// finalize predicate's `in_flight == 0` can never hold and the slot wedges
    /// — the merge consumer polls `WouldBlock` forever rather than failing.
    ///
    /// Both are caller contract violations with no in-band signal, so they are
    /// asserted here. A caller that decompresses fewer blocks than it reserved
    /// must store `true` into [`SortMergeSlot::decomp_error`] while holding the
    /// `decompressed` mutex — see the module header's ordering rules — rather
    /// than short-batching this call.
    ///
    /// # Panics
    ///
    /// Panics if `blocks.len() != count`, or if `count` exceeds the outstanding
    /// in-flight reservation. Both checks are **always on**, release included:
    /// the failure they catch is silent record loss or a wedged merge, neither
    /// of which has an in-band signal, and release is the configuration that
    /// ships. The cost is two compares per *batch* — not per record — on a path
    /// that decompresses up to [`PHASE2_DECOMP_CAP`] blocks per call. Both run
    /// before any lock is taken, so a violation cannot poison the `reorder` or
    /// `decompressed` mutex.
    ///
    /// Also panics if the `reorder` or `decompressed` mutex is already poisoned.
    pub fn bp_insert_drain_finalize(
        &self,
        start_seq: u64,
        blocks: Vec<Vec<u8>>,
        count: usize,
    ) -> bool {
        assert_eq!(
            blocks.len(),
            count,
            "bp_insert_drain_finalize would release a reservation of {count} while delivering \
             {} blocks; a short batch finalizes a clean EOF over the missing sequences",
            blocks.len()
        );
        // Load once: the message must report the value that actually failed the
        // check, not a re-read that another worker may have moved in between.
        let in_flight = self.in_flight.load(Ordering::Acquire);
        assert!(
            in_flight >= count,
            "bp_insert_drain_finalize would release {count} in-flight blocks with only \
             {in_flight} reserved; the fetch_sub below would wrap and wedge the slot",
        );
        let mut rb = lock_ranked(&self.reorder, LockRank::Reorder);
        for (i, block) in blocks.into_iter().enumerate() {
            let size = block.len();
            // `start_seq + i` cannot exceed the number of blocks read from one
            // spill file, which is far below u64::MAX.
            rb.insert_with_size(start_seq + i as u64, block, size);
        }
        self.reorder_len_mirror.store(mirror(rb.len()), Ordering::Relaxed);
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
        let mut rb = lock_ranked(&self.reorder, LockRank::Reorder);
        self.drain_locked_and_finalize(&mut rb)
    }

    /// Claim the raw stash head for decompression, or `None` when it may not
    /// be admitted yet. The head is admitted when it is the reorder buffer's
    /// next-expected sequence — always — or when the FIFO has room
    /// (`fifo_len < fifo_cap`) and the reorder window admits it
    /// (`window_admits`). The escape exists for the FIFO
    /// cap: when the head is the front, every lower sequence has been drained
    /// and the reorder buffer is empty, so the window always admits it; but the
    /// FIFO may be at its cap, and refusing the front there would leave a
    /// consumer that is about to pop its last queued block with nothing
    /// claimable behind it. The stash is in sequence order and
    /// claims pop its head, so the head is the lowest unclaimed sequence:
    /// either it is the front, or the front is already in flight. Lock order
    /// `raw_stash` → `reorder`. The claimed block stays in flight until its
    /// decompressed bytes are inserted (`bp_insert_drain_finalize`).
    ///
    /// # Panics
    /// If a slot mutex is poisoned.
    pub fn bp_claim_raw(&self, window_budget: u64) -> Option<RawBlock> {
        let mut stash = lock_ranked(&self.raw_stash, LockRank::RawStash);
        let head = stash.front()?.seq;
        let admitted = {
            let rb = lock_ranked(&self.reorder, LockRank::Reorder);
            head == rb.next_seq()
                || (self.fifo_len_relaxed() < self.fifo_cap.load(Ordering::Relaxed) as usize
                    && window_admits(&rb, head, window_budget))
        };
        if !admitted {
            return None;
        }
        let block = stash.pop_front()?;
        self.stash_len.store(mirror(stash.len()), Ordering::Relaxed);
        self.stash_bytes.fetch_sub(block.charge, Ordering::Relaxed);
        Some(block)
    }

    /// Stash `frames` (owned) exactly as [`Self::bp_ingest_slice`] does after
    /// parsing a slice — in flight first, then appended with dense sequence
    /// numbers; `last` marks the read side at EOF and drains/finalizes. Test and
    /// loom seam: loom cannot model the lease's `Arc`.
    ///
    /// # Panics
    /// If a slot mutex is poisoned.
    #[cfg(any(test, loom, feature = "test-utils"))]
    pub fn bp_stash_frames_for_test(&self, frames: Vec<Vec<u8>>, last: bool) {
        let mut r = lock_ranked(&self.reader, LockRank::Reader);
        let mut stash = lock_ranked(&self.raw_stash, LockRank::RawStash);
        let n = frames.len();
        for frame in frames {
            let charge = frame.capacity() as u64;
            self.stash_bytes.fetch_add(charge, Ordering::Relaxed);
            stash.push_back(RawBlock { seq: r.next_seq, frame: RawFrame::Owned(frame), charge });
            r.next_seq += 1;
        }
        self.publish_stashed(stash.len(), n, last);
        drop(stash);
        if last {
            self.bp_drain_and_finalize();
        }
    }

    /// Shared drain + finalize, called with the `reorder` lock held. Acquires
    /// `decompressed` (lock order: reorder → decompressed) so `queue_eof` is
    /// stored under the same mutex the consumer reads it under.
    fn drain_locked_and_finalize(&self, rb: &mut ReorderBuffer<Vec<u8>>) -> bool {
        let mut dec = lock_ranked(&self.decompressed, LockRank::Decompressed);
        let mut room = PHASE2_DECOMP_CAP.saturating_sub(dec.len());
        let mut drained = 0usize;
        while room > 0 {
            let Some(block) = rb.try_pop_next() else { break };
            dec.push_back(block);
            room -= 1;
            drained += 1;
        }
        self.fifo_len_mirror.store(mirror(dec.len()), Ordering::Relaxed);
        self.reorder_len_mirror.store(mirror(rb.len()), Ordering::Relaxed);
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

/// The slot's locks, in the one order they may be taken (outermost first).
#[derive(Clone, Copy, Debug)]
enum LockRank {
    Reader = 0,
    RawStash = 1,
    Reorder = 2,
    Decompressed = 3,
}

/// Debug-build lock-order checker: a per-thread bitmask of the slot lock ranks
/// held; taking a lock asserts that no lock of the same or a later rank is
/// held. Compiled only under `debug_assertions` (and not under loom, whose
/// model threads share one OS thread's locals).
#[cfg(all(debug_assertions, not(loom)))]
mod lock_order {
    use std::cell::Cell;

    thread_local! {
        static HELD: Cell<u8> = const { Cell::new(0) };
    }

    /// One held rank; clears it on drop.
    pub(super) struct Token(u8);

    pub(super) fn enter(rank: super::LockRank) -> Token {
        let bit = rank as u8;
        HELD.with(|held| {
            let mask = held.get();
            assert!(
                mask >> bit == 0,
                "slot lock order violated: taking {rank:?} while holding ranks {mask:#06b} \
                 (order: reader -> raw_stash -> reorder -> decompressed)"
            );
            held.set(mask | (1 << bit));
        });
        Token(bit)
    }

    impl Drop for Token {
        fn drop(&mut self) {
            HELD.with(|held| held.set(held.get() & !(1 << self.0)));
        }
    }
}

/// Release and loom builds: no checker.
#[cfg(not(all(debug_assertions, not(loom))))]
mod lock_order {
    pub(super) struct Token;

    pub(super) fn enter(_rank: super::LockRank) -> Token {
        Token
    }
}

#[cfg(loom)]
type Guard<'a, T> = loom::sync::MutexGuard<'a, T>;
#[cfg(not(loom))]
type Guard<'a, T> = std::sync::MutexGuard<'a, T>;

/// A slot lock's guard paired with its lock-order token. The guard is
/// declared first, so it is released before the token clears the rank.
struct Ranked<'a, T> {
    guard: Guard<'a, T>,
    _token: lock_order::Token,
}

impl<T> std::ops::Deref for Ranked<'_, T> {
    type Target = T;
    fn deref(&self) -> &T {
        &self.guard
    }
}

impl<T> std::ops::DerefMut for Ranked<'_, T> {
    fn deref_mut(&mut self) -> &mut T {
        &mut self.guard
    }
}

/// Take one of a slot's locks at `rank` (checked in debug builds).
///
/// # Panics
/// If the mutex is poisoned, or (debug builds) the lock order is violated.
fn lock_ranked<T>(m: &Mutex<T>, rank: LockRank) -> Ranked<'_, T> {
    let token = lock_order::enter(rank);
    let guard = m.lock().unwrap_or_else(|_| panic!("SortMergeSlot {rank:?} mutex poisoned"));
    Ranked { guard, _token: token }
}

/// The reorder window's admission rule for a claim of sequence `seq`:
/// [`ReorderBuffer::would_accept`] (the deadlock-free predicate the pipeline
/// uses elsewhere) plus a hard `heap_bytes < window_budget` backstop
/// (`window_budget == 0` = unlimited). The backstop is what bounds memory:
/// `would_accept` alone accepts everything while the buffer is *stuck* (its
/// front sequence not yet decompressed), which would let a straggler balloon
/// the window. It is deadlock-safe here because parsing is serialized and
/// densely sequenced, so the front gap is always a claimed, in-flight block —
/// never an unread one a refused claim would have had to fetch — and the
/// front itself bypasses this rule ([`SortMergeSlot::bp_claim_raw`]).
fn window_admits<T>(rb: &ReorderBuffer<T>, seq: u64, window_budget: u64) -> bool {
    rb.would_accept(seq, window_budget) && (window_budget == 0 || rb.heap_bytes() < window_budget)
}

/// A queue length as a `u32` mirror value (a slot's queues hold far fewer
/// than `u32::MAX` blocks).
fn mirror(len: usize) -> u32 {
    u32::try_from(len).unwrap_or(u32::MAX)
}

// These exercise `SortMergeSlot` outside `loom::model`, which is illegal once
// the primitives are loom's, so they compile only in a normal (non-loom) build.
// The `--cfg loom` invocation runs `tests/loom_merge_slots.rs` exclusively.
#[cfg(all(test, not(loom)))]
mod tests {
    use rstest::rstest;

    use super::*;

    fn frames(n: u8) -> Vec<Vec<u8>> {
        (0..n).map(|i| vec![i; 100]).collect()
    }

    /// A head that is not the reorder front is refused while
    /// the window holds bytes at its budget, and admitted once the front lands
    /// and it becomes the front. (The front itself never meets a full window:
    /// when the head is the front the reorder buffer is empty. The escape's
    /// case is the FIFO cap — `claim_admits_the_front_over_a_full_fifo`.)
    #[test]
    fn claim_refuses_a_later_head_over_the_window_budget() {
        let slot = SortMergeSlot::for_test(0, SpillCodec::Bgzf);
        slot.bp_stash_frames_for_test(frames(3), true);
        let b0 = slot.bp_claim_raw(1).expect("the front is always admitted");
        assert_eq!(b0.seq, 0);
        // Seq 1 is not the front (0 is in flight): the 1-byte window refuses
        // it only once it holds bytes, so fill the window with a later seq.
        let b1 = slot.bp_claim_raw(1).expect("an empty window admits");
        assert!(!slot.bp_insert_drain_finalize(1, vec![vec![1; 64]], 1), "seq 1 waits for 0");
        assert_eq!(slot.bp_reorder_heap_bytes(), 64, "seq 1 waits behind the front");
        assert!(slot.bp_claim_raw(1).is_none(), "seq 2 is not the front and the window is full");
        assert!(slot.bp_insert_drain_finalize(0, vec![vec![0; 64]], 1));
        let b2 = slot.bp_claim_raw(1).expect("seq 2 is now the front");
        assert_eq!((b1.seq, b2.seq), (1, 2));
    }

    /// Claim admission reads the FIFO mirror the consumer's pop stores: with
    /// the FIFO at its cap and the front in flight, a later head is refused
    /// until the consumer pops, then admitted.
    #[test]
    fn claim_admission_follows_the_consumer_pop() {
        let slot = SortMergeSlot::for_test(0, SpillCodec::Bgzf);
        slot.fifo_cap.store(1, Ordering::Relaxed);
        slot.bp_stash_frames_for_test(frames(3), true);
        let b0 = slot.bp_claim_raw(u64::MAX).unwrap();
        slot.bp_insert_drain_finalize(b0.seq, vec![vec![0]], 1);
        let b1 = slot.bp_claim_raw(u64::MAX).expect("the front");
        assert!(slot.bp_claim_raw(u64::MAX).is_none(), "seq 2: FIFO at its cap");
        assert_eq!(slot.pop_decompressed(), Some(vec![0]));
        let b2 = slot.bp_claim_raw(u64::MAX).expect("the pop made room");
        assert_eq!((b1.seq, b2.seq), (1, 2));
    }

    /// The front escape: with the FIFO at its cap the head that is the reorder
    /// front is still admitted; a later head is refused.
    #[test]
    fn claim_admits_the_front_over_a_full_fifo() {
        let slot = SortMergeSlot::for_test(0, SpillCodec::Bgzf);
        slot.fifo_cap.store(1, Ordering::Relaxed);
        slot.bp_stash_frames_for_test(frames(3), true);
        let b0 = slot.bp_claim_raw(u64::MAX).unwrap();
        slot.bp_insert_drain_finalize(b0.seq, vec![vec![0]], 1);
        assert_eq!(slot.fifo_len_relaxed(), 1, "the FIFO is at its cap");
        let b1 = slot.bp_claim_raw(u64::MAX).expect("seq 1 is the front: admitted over the cap");
        assert_eq!(b1.seq, 1);
        assert!(slot.bp_claim_raw(u64::MAX).is_none(), "seq 2: not the front, FIFO at cap");
    }

    #[test]
    fn claim_pops_in_seq_order() {
        let slot = SortMergeSlot::for_test(0, SpillCodec::Bgzf);
        slot.bp_stash_frames_for_test(frames(5), false);
        let seqs: Vec<u64> =
            std::iter::from_fn(|| slot.bp_claim_raw(u64::MAX)).map(|b| b.seq).collect();
        assert_eq!(seqs, vec![0, 1, 2, 3, 4]);
        assert_eq!(slot.stash_len.load(Ordering::Relaxed), 0);
        assert_eq!(slot.stash_bytes.load(Ordering::Relaxed), 0);
    }

    #[rstest]
    #[case::starved(0, 0, 0, crate::AwaitedSlotState::Starved)]
    #[case::issued(4096, 0, 0, crate::AwaitedSlotState::Issued)]
    #[case::stashed(0, 2, 0, crate::AwaitedSlotState::Stashed)]
    #[case::decompressing(0, 2, 1, crate::AwaitedSlotState::Decompressing)]
    fn awaited_state_reads_the_mirrors(
        #[case] issued: u64,
        #[case] stashed: u8,
        #[case] claimed: usize,
        #[case] want: crate::AwaitedSlotState,
    ) {
        let slot = SortMergeSlot::for_test(0, SpillCodec::Bgzf);
        slot.bp_note_issued(issued);
        slot.bp_stash_frames_for_test(frames(stashed), false);
        let _claimed: Vec<_> = (0..claimed).filter_map(|_| slot.bp_claim_raw(u64::MAX)).collect();
        assert_eq!(slot.awaited_state(), want);
    }

    /// Debug builds check the lock order: taking an outer lock while holding
    /// an inner one panics.
    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "slot lock order violated")]
    fn lock_order_violation_panics_in_debug() {
        let slot = SortMergeSlot::for_test(0, SpillCodec::Bgzf);
        let _dec = lock_ranked(&slot.decompressed, LockRank::Decompressed);
        let _stash = lock_ranked(&slot.raw_stash, LockRank::RawStash);
    }

    /// `mark_failed` sets both flags, read under the lock the consumer uses:
    /// an error, never a clean drain.
    #[test]
    fn mark_failed_sets_both_flags_under_the_lock() {
        let slot = SortMergeSlot::for_test(0, SpillCodec::Bgzf);
        slot.mark_failed();
        assert!(slot.has_error());
        assert!(slot.queue_eof.load(Ordering::Acquire));
        assert!(!slot.is_drained(), "a failed slot is not a clean drain");
    }

    #[test]
    fn new_slot_starts_empty_and_not_drained() {
        let slot = SortMergeSlot::for_test(0, SpillCodec::Bgzf);
        assert_eq!(slot.file_id, 0);
        assert!(!slot.is_drained(), "fresh slot not drained until queue_eof");
        let (pending_blocks, pending_bytes, active) = slot.probe_stats();
        assert_eq!(pending_blocks, 0);
        assert_eq!(pending_bytes, 0);
        assert!(active, "active until queue_eof flips");
    }

    #[test]
    fn drained_when_queue_eof_and_empty() {
        let slot = SortMergeSlot::for_test(0, SpillCodec::Bgzf);

        // Push a decompressed block; not drained even if queue_eof.
        slot.decompressed.lock().unwrap().push_back(vec![0xAB, 0xCD]);
        slot.queue_eof.store(true, Ordering::Release);
        assert!(!slot.is_drained(), "not drained while decompressed non-empty");

        // Consume the block.
        let popped = slot.pop_decompressed();
        assert_eq!(popped, Some(vec![0xAB, 0xCD]));
        assert!(slot.is_drained(), "drained once queue empty AND queue_eof");
    }

    #[test]
    fn fifo_len_reports_block_count() {
        let slot = SortMergeSlot::for_test(0, SpillCodec::Bgzf);
        assert_eq!(slot.fifo_len(), 0);
        slot.decompressed.lock().unwrap().push_back(vec![0u8; 10]);
        slot.decompressed.lock().unwrap().push_back(vec![0u8; 20]);
        assert_eq!(slot.fifo_len(), 2, "counts blocks, not bytes");
    }

    /// The merge's re-check must treat EOF as "can progress" — otherwise the
    /// merge parks at a run end where no further delivery will wake it.
    #[test]
    fn has_block_or_eof_is_true_for_a_block_or_eof_only() {
        let slot = SortMergeSlot::for_test(0, SpillCodec::Bgzf);
        assert!(!slot.has_block_or_eof(), "fresh slot: no block, no EOF");
        slot.decompressed.lock().unwrap().push_back(vec![1]);
        assert!(slot.has_block_or_eof(), "a queued block");
        let eof = SortMergeSlot::for_test(1, SpillCodec::Bgzf);
        eof.queue_eof.store(true, Ordering::Release);
        assert!(eof.has_block_or_eof(), "EOF with an empty FIFO");
    }

    #[test]
    fn probe_stats_counts_decompressed_only() {
        let slot = SortMergeSlot::for_test(0, SpillCodec::Bgzf);
        slot.decompressed.lock().unwrap().push_back(vec![0u8; 1024]);
        slot.decompressed.lock().unwrap().push_back(vec![0u8; 512]);

        let (pending_blocks, pending_bytes, active) = slot.probe_stats();
        assert_eq!(pending_blocks, 2);
        assert_eq!(pending_bytes, 1024 + 512);
        assert!(active, "not yet queue_eof");
    }

    /// Blocks reach the FIFO in READ order even when decompression completes
    /// out of order — the reorder buffer's entire purpose. Hand-pushing into
    /// `decompressed` and popping it back would assert that `VecDeque` is a
    /// queue, which is true of `VecDeque` and says nothing about the slot.
    #[test]
    fn fifo_order() {
        let slot = SortMergeSlot::for_test(7, SpillCodec::Bgzf);
        slot.bp_commit_read(4, true);

        // The later half of the file finishes decompressing first. Nothing may
        // drain past the gap at seq 0, however ready seqs 2 and 3 are.
        slot.bp_insert_drain_finalize(2, vec![vec![2], vec![3]], 2);
        assert_eq!(slot.fifo_len(), 0, "nothing drains while seq 0 is missing");

        // The front of the file lands; now the whole run drains, in read order.
        slot.bp_insert_drain_finalize(0, vec![vec![0], vec![1]], 2);
        let popped: Vec<u8> = slot.decompressed.lock().unwrap().drain(..).map(|b| b[0]).collect();
        assert_eq!(popped, vec![0, 1, 2, 3], "read order, not completion order");
        assert!(slot.queue_eof.load(Ordering::Acquire), "all 4 delivered ⇒ EOF finalizes");
    }

    /// `has_error` must report the flag through its own lock-and-load, not
    /// merely echo the atomic. Asserting on `decomp_error` directly would test
    /// `AtomicBool` rather than anything belonging to `SortMergeSlot`.
    #[test]
    fn has_error_reports_the_decomp_error_flag() {
        let slot = SortMergeSlot::for_test(0, SpillCodec::Bgzf);
        assert!(!slot.has_error(), "a fresh slot has no error");
        slot.decomp_error.store(true, Ordering::Release);
        assert!(slot.has_error(), "has_error must observe a stored decomp_error");
    }

    #[test]
    fn errored_slot_is_not_drained() {
        let slot = SortMergeSlot::for_test(0, SpillCodec::Bgzf);

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
        let slot = SortMergeSlot::for_test(0, SpillCodec::Bgzf);
        slot.queue_eof.store(true, Ordering::Release);
        assert!(slot.is_drained(), "clean empty + queue_eof is drained");
        assert!(!slot.has_error(), "clean drain has no error");
    }

    // ── Block-parallel helper tests ─────────────────────────────────────────

    /// A single in-order batch drains straight into the FIFO and (because the
    /// reader is at EOF with no other in-flight blocks) finalizes `queue_eof`.
    #[test]
    fn bp_single_batch_drains_and_finalizes() {
        let slot = SortMergeSlot::for_test(0, SpillCodec::Bgzf);
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
        let slot = SortMergeSlot::for_test(0, SpillCodec::Bgzf);

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
        let slot = SortMergeSlot::for_test(0, SpillCodec::Bgzf);

        // Stamp each block with its sequence number. Asserting only a count at
        // the end would pass for a drain that dropped one block and duplicated
        // another, or that delivered the window out of order.
        let block_for = |seq: u64| {
            let mut block = vec![0u8; BLOCK];
            block[0] = u8::try_from(seq).expect("the window budget keeps seqs far below 256");
            block
        };

        // Seq 0 is the straggler — it is reserved but never decompressed/inserted,
        // so the buffer can never pop and keeps a permanent front gap.
        slot.bp_add_in_flight(1); // seq 0 in flight forever (the straggler)

        // Workers race ahead delivering seqs 1,2,3,… as long as admission allows.
        let mut next = 1u64;
        let mut admitted = 0;
        while slot.bp_reorder_admits(next, BUDGET) {
            slot.bp_add_in_flight(1);
            slot.bp_insert_drain_finalize(next, vec![block_for(next)], 1);
            // Front gap at seq 0 ⇒ nothing drains.
            assert!(
                slot.decompressed.lock().unwrap().is_empty(),
                "no block may reach the FIFO while the seq-0 gap is open",
            );
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

        // Now the part that makes the EOF check mean something. Until here
        // `reader_eof` was never set, so asserting `!queue_eof` inside the loop
        // was vacuous — the finalize predicate could not fire for want of
        // `reader_eof`, whatever the straggler did. Set it, and the predicate
        // becomes falsifiable: the only thing still holding EOF back is seq 0,
        // which is in flight and absent from `reorder`.
        slot.bp_set_reader_eof();
        assert!(!slot.bp_drain_and_finalize(), "nothing can drain past the seq-0 gap");
        assert!(
            !slot.queue_eof.load(Ordering::Acquire),
            "must not finalize EOF while the straggler at seq 0 is still in flight",
        );

        // Deliver the straggler; now the whole window drains in order and EOF
        // finalizes, confirming the blocks held behind the gap were retained
        // rather than dropped.
        slot.bp_insert_drain_finalize(0, vec![block_for(0)], 1);
        assert!(
            slot.queue_eof.load(Ordering::Acquire),
            "EOF must finalize once the straggler lands and the window drains",
        );
        let delivered: Vec<u8> =
            slot.decompressed.lock().unwrap().drain(..).map(|block| block[0]).collect();
        let expected: Vec<u8> = (0..=u8::try_from(admitted)
            .expect("admitted is bounded by the window budget"))
            .collect();
        assert_eq!(
            delivered, expected,
            "the straggler and every admitted block must reach the FIFO exactly once, in read \
             order — not merely in the right quantity",
        );
    }

    /// The drain stops at [`PHASE2_DECOMP_CAP`], so a slot can sit at
    /// reader-EOF with `in_flight == 0` and *still* owe blocks that did not fit
    /// in the FIFO. Finalizing `queue_eof` there would strand them: the
    /// consumer stops at `is_drained()`, and the blocks left in `reorder` are
    /// never popped. Completion instead depends on the driver calling
    /// [`SortMergeSlot::bp_drain_and_finalize`] again once the consumer has made
    /// room — the one drain path whose correctness lives outside this module.
    ///
    /// Two boundaries: a FIFO exactly at cap (`room == 0`, the loop never runs)
    /// and one block short of it (`room` hits zero *inside* the loop, after a
    /// partial delivery). Both must defer the finalize.
    #[rstest]
    #[case::fifo_exactly_at_cap(PHASE2_DECOMP_CAP, 0, false)]
    #[case::fifo_one_short_of_cap(PHASE2_DECOMP_CAP - 1, 1, true)]
    fn bp_full_fifo_defers_finalize_until_the_consumer_makes_room(
        #[case] prefill: usize,
        #[case] drained_now: usize,
        #[case] progressed_now: bool,
    ) {
        const DEFERRED: [u8; 2] = [0xAA, 0xBB];
        let slot = SortMergeSlot::for_test(0, SpillCodec::Bgzf);

        // Occupy the FIFO so the incoming batch cannot fully land. The prefill
        // bytes stay below 0xAA, so they never alias the delivered blocks.
        // A constant byte, distinct from DEFERRED: the prefill is only there to
        // occupy room, and keying it to `i` would make the test panic the day
        // PHASE2_DECOMP_CAP is tuned above 256.
        {
            let mut dec = slot.decompressed.lock().unwrap();
            for _ in 0..prefill {
                dec.push_back(vec![0u8]);
            }
        }

        // The final read: two blocks, reader at EOF, nothing else in flight.
        // Every finalize precondition except `reorder.is_empty()` now holds.
        slot.bp_commit_read(2, true);
        let progressed =
            slot.bp_insert_drain_finalize(0, DEFERRED.iter().map(|&b| vec![b]).collect(), 2);

        assert_eq!(progressed, progressed_now, "progress is reported iff a block actually drained");
        assert_eq!(slot.fifo_len(), prefill + drained_now, "the FIFO fills only to the cap");
        assert!(
            !slot.queue_eof.load(Ordering::Acquire),
            "{} block(s) are still owed from `reorder`; finalizing EOF here strands them",
            DEFERRED.len() - drained_now,
        );
        assert!(!slot.is_drained(), "a slot that still owes blocks is not drained");
        assert_eq!(slot.bp_fifo_room(), 0, "the drain ran until the FIFO had no room left");

        // The consumer pops; the driver's follow-up drain must deliver the
        // remainder in read order and only then finalize.
        {
            let mut dec = slot.decompressed.lock().unwrap();
            for _ in 0..DEFERRED.len() {
                dec.pop_front();
            }
        }
        assert!(slot.bp_drain_and_finalize(), "room freed ⇒ the deferred blocks drain");
        assert!(
            slot.queue_eof.load(Ordering::Acquire),
            "EOF finalizes once `reorder` empties and nothing is in flight",
        );

        let delivered: Vec<Vec<u8>> = slot.decompressed.lock().unwrap().drain(..).collect();
        assert_eq!(
            delivered.len(),
            prefill,
            "prefill + {} delivered, less the {} the consumer popped: a dropped or duplicated \
             block moves this count",
            DEFERRED.len(),
            DEFERRED.len(),
        );
        assert_eq!(
            delivered[delivered.len() - DEFERRED.len()..],
            DEFERRED.map(|b| vec![b]),
            "the deferred blocks arrive behind the prefill, in read order",
        );
    }

    /// The two caller-contract violations the assertions exist to catch: a
    /// short batch (silent truncation) and an over-release (wraps `in_flight`
    /// and wedges the slot). The assertions are always on, so these tests run
    /// in release too — which is the configuration whose behaviour they pin.
    #[test]
    #[should_panic(expected = "while delivering")]
    fn bp_short_batch_is_rejected() {
        let slot = SortMergeSlot::for_test(0, SpillCodec::Bgzf);
        slot.bp_commit_read(2, true);
        // Deliver 1 of the 2 reserved blocks: without the assertion this
        // finalizes a clean queue_eof over the missing seq 1.
        slot.bp_insert_drain_finalize(0, vec![vec![0xAAu8; 8]], 2);
    }

    #[test]
    #[should_panic(expected = "would wrap")]
    fn bp_over_release_is_rejected() {
        let slot = SortMergeSlot::for_test(0, SpillCodec::Bgzf);
        slot.bp_commit_read(1, true);
        // Release 2 with only 1 reserved: the fetch_sub would wrap in_flight
        // to usize::MAX and the slot could never finalize.
        slot.bp_insert_drain_finalize(0, vec![vec![1u8; 8], vec![2u8; 8]], 2);
    }
}

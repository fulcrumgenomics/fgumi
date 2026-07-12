//! `PairRawFastq` step: a cheap two-input chunk-level pairing step that joins
//! the `FastqRawChunk` streams emitted by two concurrent per-stream
//! [`ReadFastqInputs`] readers (R1 + R2) into a single
//! [`PairedRawFastqBatch`] stream for the downstream `Parallel`
//! [`ParseAndZipFastq`] step.
//!
//! `Serial` + `Affinity::None`. This step mirrors the legacy pipeline's
//! `FindBoundaries` stage: it does the cheap chunk-level pairing serially
//! (matching the R1 chunk and R2 chunk at the same `chunk_serial`) and assigns
//! a fresh, globally-unique, monotonically-increasing `ordinal` to each emitted
//! pair. The expensive per-record FASTQ parse and template build then run in
//! parallel in [`ParseAndZipFastq`].
//!
//! ## Why serial ordinal assignment is load-bearing
//!
//! [`ParseAndZipFastq`] is `Parallel` + `ByItemOrdinal`, so the framework
//! reorders its output by each item's `ordinal()`. For that reorder to be
//! correct the ordinals MUST be globally unique and dense. Assigning them here,
//! serially, from a single `next_ordinal` counter on this step guarantees both
//! properties. (The per-`FastqRawChunk` ordinal from the readers is NOT reused
//! — two readers each carry their own ordinal sequence; we mint a fresh one for
//! the paired batch.)
//!
//! ## Completion
//!
//! `try_run` opportunistically buffers whichever branch has a chunk ready and
//! emits a pair as soon as both streams' chunks for the lowest `chunk_serial`
//! are present. Once **both** input branches are drained, `finalize_pairs`
//! flushes the remaining complete pairs in `chunk_serial` order (one per call)
//! and returns `StepOutcome::Finished`, bailing with an error if either stream
//! ran short.

use std::collections::BTreeMap;
use std::io;

use super::read_fastq::FastqRawChunk;
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::item::{HeapSize, Ordered};
use crate::pipeline::core::step::Affinity;
use crate::pipeline::core::{
    BranchOrdering, OrderedBytesSingle, QueueSpec, Step2, StepCtx2, StepKind, StepOutcome,
    StepProfile, Unpushed,
};

/// Backpressure cap on the total bytes buffered across the two partial-pair
/// slots. Mirrors `ZipFastqRecords::DEFAULT_PENDING_BACKPRESSURE_BYTES`.
const DEFAULT_PENDING_BACKPRESSURE_BYTES: usize = 256 * 1024 * 1024;

// ─────────────────────────────────────────────────────────────────────────────
// PairedRawFastqBatch — the raw bytes of one R1 chunk + one R2 chunk that share
// a chunk_serial, plus a freshly-minted globally-unique ordinal.
// ─────────────────────────────────────────────────────────────────────────────

/// A matched pair of raw (decompressed, unparsed) FASTQ byte chunks — one from
/// stream A (R1) and one from stream B (R2) — that share a `chunk_serial`.
///
/// `ordinal` is minted serially by [`PairRawFastq`] and is the key the
/// framework uses to reorder the `Parallel` [`ParseAndZipFastq`] output.
/// `chunk_serial` is carried through for diagnostics / desync messages only.
///
/// [`ParseAndZipFastq`]: super::parse_zip_fastq::ParseAndZipFastq
pub struct PairedRawFastqBatch {
    /// Globally-unique monotonic ordinal, minted by [`PairRawFastq`].
    pub ordinal: u64,
    /// The shared per-stream round-robin cycle serial (for desync diagnostics).
    pub chunk_serial: u64,
    /// Raw bytes of the stream-A (R1) chunk, whole-record-aligned.
    pub data_a: Vec<u8>,
    /// Raw bytes of the stream-B (R2) chunk, whole-record-aligned.
    pub data_b: Vec<u8>,
}

impl HeapSize for PairedRawFastqBatch {
    fn heap_size(&self) -> usize {
        self.data_a.len() + self.data_b.len()
    }
}

impl Ordered for PairedRawFastqBatch {
    fn ordinal(&self) -> u64 {
        self.ordinal
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PairRawFastq — Serial two-input pairing step
// ─────────────────────────────────────────────────────────────────────────────

/// `Serial` two-input chunk-level pairing step. Buffers the most recent
/// unmatched chunk from each stream keyed by `chunk_serial`; once both streams
/// have produced their chunk for the lowest pending `chunk_serial`, emits a
/// [`PairedRawFastqBatch`] with a fresh `ordinal`.
pub struct PairRawFastq {
    /// Partial pairs keyed by `chunk_serial`. Each value is `(Option<a>,
    /// Option<b>)`; a complete pair is emitted (and removed) as soon as both
    /// slots are filled and it is the lowest pending serial.
    pending: BTreeMap<u64, (Option<FastqRawChunk>, Option<FastqRawChunk>)>,
    pending_total_bytes: usize,
    pending_backpressure_bytes: usize,
    /// Serial, monotonic ordinal assignment point for the downstream parallel
    /// reorder. Incremented once per emitted pair.
    next_ordinal: u64,
    held: HeldSlot<Unpushed<PairedRawFastqBatch>>,
    output_byte_limit: u64,
}

impl PairRawFastq {
    #[must_use]
    pub fn new(output_byte_limit: u64) -> Self {
        Self {
            pending: BTreeMap::new(),
            pending_total_bytes: 0,
            pending_backpressure_bytes: DEFAULT_PENDING_BACKPRESSURE_BYTES,
            next_ordinal: 0,
            held: HeldSlot::new(),
            output_byte_limit,
        }
    }

    /// Override the soft pending-backpressure limit (test-only).
    ///
    /// Production wires the default 256 MiB limit; tests shrink it to a few
    /// bytes so a skewed-completion scenario can cross it with a handful of
    /// small chunks rather than hundreds of megabytes of data.
    #[cfg(test)]
    fn with_backpressure_bytes(mut self, bytes: usize) -> Self {
        self.pending_backpressure_bytes = bytes;
        self
    }

    /// Buffer a chunk from stream A (slot 0) or stream B (slot 1).
    ///
    /// Each reader mints exactly one chunk per `(chunk_serial, stream)`, so a
    /// target slot is always `None` here. The byte accounting (X5-004) depends
    /// on that invariant: it adds `chunk.data.len()` to `pending_total_bytes`,
    /// which feeds the catastrophic-desync guard and the backpressure decision.
    /// A *hard* `assert!` enforces single-chunk-per-serial in release builds
    /// too: overwriting an already-filled slot would silently change emitted
    /// read identity (and drift the byte total), so a retry/replay/re-chunk
    /// regression must fail fast rather than corrupt output.
    fn buffer_chunk(&mut self, chunk: FastqRawChunk, is_a: bool) {
        self.pending_total_bytes += chunk.data.len();
        let entry = self.pending.entry(chunk.chunk_serial).or_insert((None, None));
        let slot = if is_a { &mut entry.0 } else { &mut entry.1 };
        assert!(
            slot.is_none(),
            "duplicate chunk for serial {} stream {}: single-chunk-per-serial invariant violated",
            chunk.chunk_serial,
            if is_a { "A" } else { "B" },
        );
        *slot = Some(chunk);
    }

    /// If the lowest pending `chunk_serial` has both stream slots filled, remove
    /// it and produce a `PairedRawFastqBatch` with a fresh ordinal.
    fn try_emit_complete(&mut self) -> Option<PairedRawFastqBatch> {
        let lowest_serial = *self.pending.keys().next()?;
        let entry = self.pending.get(&lowest_serial)?;
        if entry.0.is_none() || entry.1.is_none() {
            return None;
        }

        let (a, b) = self.pending.remove(&lowest_serial).unwrap();
        let a = a.unwrap();
        let b = b.unwrap();
        self.pending_total_bytes =
            self.pending_total_bytes.saturating_sub(a.data.len() + b.data.len());

        let ordinal = self.next_ordinal;
        self.next_ordinal += 1;
        Some(PairedRawFastqBatch {
            ordinal,
            chunk_serial: lowest_serial,
            data_a: a.data,
            data_b: b.data,
        })
    }
}

/// Decide which input branches to pull this dispatch, given backpressure state.
///
/// Returns `(pull_a, pull_b)`. The pending-pair buffer must be byte-bounded, but
/// bounding it by globally refusing to pull deadlocks: when one FASTQ stream
/// races ahead, the buffer fills with its unmatched chunks, and the chunk needed
/// to complete the lowest pending pair sits unpulled in the *behind* stream's
/// queue. So under backpressure we still pull whichever stream is missing from
/// the front (lowest pending) pair — the pull that lets the front pair emit and
/// the buffer drain — and throttle only the stream that is already ahead.
///
/// - Not over the limit → pull both (normal streaming).
/// - Over the limit, front pair missing a slot → pull only that (behind) stream.
/// - Over the limit, front pair already complete (or no pending) → pull neither;
///   the caller drains via `try_emit_complete` instead of accumulating more.
fn pull_decision(
    over_backpressure: bool,
    front_missing_a: bool,
    front_missing_b: bool,
) -> (bool, bool) {
    if !over_backpressure {
        return (true, true);
    }
    (front_missing_a, front_missing_b)
}

impl Step2 for PairRawFastq {
    type InputA = FastqRawChunk;
    type InputB = FastqRawChunk;
    type Outputs = OrderedBytesSingle<PairedRawFastqBatch>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "PairRawFastq",
            kind: StepKind::Serial,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            // FIFO (no reorder stage) is BOTH correct and required for
            // throughput here. Each emitted `PairedRawFastqBatch` already
            // carries a dense, globally-unique `ordinal` minted serially by
            // this step, and the downstream `ParseAndZipFastq` is Parallel +
            // `ByItemOrdinal`, so it re-establishes input order by that
            // ordinal regardless of arrival order. Imposing `ByItemOrdinal`
            // on THIS edge would insert a `ReorderStage` that releases items
            // strictly one-at-a-time in ordinal order, serializing the
            // Parallel consumer (its workers would contend on a single
            // `next_serial` slot) and defeating the whole point of moving the
            // parse + template build into a Parallel step. This mirrors the
            // legacy `MergeFastqRawChunks` → `ParseFastqChunks` edge, which
            // was FIFO for the same reason.
            branch_ordering: vec![BranchOrdering::None],
        }
    }

    fn affinity(&self) -> Affinity {
        // Let any free worker drive the pairing; pinning it to a reader's
        // worker would serialize read + pair on one thread.
        Affinity::None
    }

    fn try_run(&mut self, ctx: &mut StepCtx2<'_, Self>) -> io::Result<StepOutcome> {
        // 1. Drain held output slot first.
        if let Some(unpushed) = self.held.take() {
            match ctx.outputs.retry(unpushed) {
                Ok(()) => {}
                Err(again) => {
                    self.held.put(again);
                    return Ok(StepOutcome::Contention);
                }
            }
        }

        // 2. Catastrophic-desync guard: one stream producing far more data
        //    than its counterpart balloons the partial-pair buffer. With the
        //    per-stream backpressure below this should never trip (the ahead
        //    stream is throttled at 1x), but keep it as a safety net.
        if self.pending_total_bytes > 2 * self.pending_backpressure_bytes {
            return Err(io::Error::other(format!(
                "PairRawFastq: catastrophic stream desync — pending buffer ({} bytes) \
                 exceeds 2x backpressure limit ({} bytes). One FASTQ stream is producing \
                 far more data than its counterpart.",
                self.pending_total_bytes, self.pending_backpressure_bytes
            )));
        }

        // 3. Drain available chunks and emit every pair we can in this dispatch.
        //    Doing this in a loop (rather than one chunk per call) amortizes the
        //    Serial mutex / dispatch overhead and keeps the downstream Parallel
        //    `ParseAndZipFastq` workers fed.
        //
        //    Backpressure is PER-STREAM, not global: once the pending buffer is
        //    over the limit we keep pulling whichever stream is missing from the
        //    lowest pending pair (the pull that completes it and drains the
        //    buffer) and throttle only the stream that is ahead. Globally
        //    refusing to pull here deadlocks — the completing chunk sits unpulled
        //    in the behind stream's full queue while the buffer never drains
        //    (a fast aligner makes the front-end run fast enough to desync the
        //    R1/R2 readers past the limit and trip this).
        let mut pulled = false;
        let mut emitted = false;
        loop {
            let over_backpressure = self.pending_total_bytes > self.pending_backpressure_bytes;
            let (front_missing_a, front_missing_b) = match self.pending.values().next() {
                Some((a, b)) => (a.is_none(), b.is_none()),
                None => (false, false),
            };
            let (pull_a, pull_b) =
                pull_decision(over_backpressure, front_missing_a, front_missing_b);

            let mut pulled_this_iter = false;
            if pull_a {
                if let Some(chunk) = ctx.a.pop() {
                    self.buffer_chunk(chunk, /* is_a */ true);
                    pulled_this_iter = true;
                }
            }
            if pull_b {
                if let Some(chunk) = ctx.b.pop() {
                    self.buffer_chunk(chunk, /* is_a */ false);
                    pulled_this_iter = true;
                }
            }
            pulled |= pulled_this_iter;

            // Emit every complete pair currently at the front of the buffer.
            while let Some(pair) = self.try_emit_complete() {
                emitted = true;
                if let Err(unpushed) = ctx.outputs.push(pair) {
                    // Output queue is full: hold the rejected pair and yield so
                    // the downstream consumer can drain it.
                    self.held.put(unpushed);
                    return Ok(StepOutcome::Progress);
                }
            }

            // Stop when neither branch yielded a new chunk this iteration.
            if !pulled_this_iter {
                break;
            }
        }

        if pulled || emitted {
            return Ok(StepOutcome::Progress);
        }

        // Nothing pulled or emitted this call. When over backpressure, pulling
        // of the surplus (ahead) stream has stopped, so if the front pending pair
        // is missing a record from a stream that is already drained, that stream
        // can never complete the pair AND the ahead stream can never drain —
        // `finalize_pairs` (both-drained) is unreachable, and `pending_total_bytes`
        // plateaus below the 2x catastrophic-desync abort so that guard never
        // fires either. This is the unequal-record-count case (e.g. a truncated
        // R1): fail fast instead of spinning `NoProgress` forever (a silent hang).
        // Under normal backpressure both streams still drain and the richer
        // `finalize_pairs` desync report handles the mismatch.
        if self.pending_total_bytes > self.pending_backpressure_bytes {
            if let Some((&serial, (a_slot, b_slot))) = self.pending.iter().next() {
                let a_short = a_slot.is_none() && ctx.a.is_drained();
                let b_short = b_slot.is_none() && ctx.b.is_drained();
                if a_short || b_short {
                    // The stream whose front slot is missing is the one that ended
                    // early. Report through the shared helper so this mid-stream
                    // fail-fast names the short stream and orphaned serial exactly
                    // like the both-drained `finalize_pairs` path.
                    let short_stream = usize::from(!a_short);
                    return Err(out_of_sync_error(short_stream, serial));
                }
            }
        }

        // If both input branches are drained, the pairing is complete — flush
        // remaining pairs and finish.
        if ctx.a.is_drained() && ctx.b.is_drained() {
            return self.finalize_pairs(ctx);
        }
        Ok(StepOutcome::NoProgress)
    }
}

/// Build the canonical FASTQ out-of-sync error identifying the stream that ended
/// early. `short_stream` is the 0-based index of the stream that ran out of
/// records first (0 = R1/A, 1 = R2/B) and `chunk_serial` is the orphaned serial
/// still holding only the other stream's chunk.
///
/// Shared by both fail paths — the both-drained [`PairRawFastq::finalize_pairs`]
/// desync report and the over-backpressure mid-stream fail-fast — so the two
/// surface byte-for-byte identical diagnostics regardless of which path detects
/// the mismatch.
fn out_of_sync_error(short_stream: usize, chunk_serial: u64) -> io::Error {
    io::Error::other(format!(
        "FASTQ sources out of sync: stream {short_stream} ended before chunk_serial \
         {chunk_serial} while the other stream had more records"
    ))
}

impl PairRawFastq {
    /// Both input branches are drained: emit the next remaining complete pair
    /// (in `chunk_serial` order), or report `Finished` once `pending` is empty.
    /// Any serial holding only one stream's chunk means the two FASTQ streams
    /// desynchronized (unequal record counts) — a fatal error. One pair per
    /// call (re-dispatched); a bounced push is parked in `held` and retried by
    /// `try_run`'s held-drain preamble next pass.
    fn finalize_pairs(&mut self, ctx: &mut StepCtx2<'_, Self>) -> io::Result<StepOutcome> {
        let Some(&lowest_serial) = self.pending.keys().next() else {
            return Ok(StepOutcome::Finished);
        };
        let entry = self.pending.get(&lowest_serial).unwrap();
        if entry.0.is_none() || entry.1.is_none() {
            let short = usize::from(entry.0.is_some());
            return Err(out_of_sync_error(short, lowest_serial));
        }
        let pair = self.try_emit_complete().expect("entry is complete");
        if let Err(unpushed) = ctx.outputs.push(pair) {
            self.held.put(unpushed);
        }
        Ok(StepOutcome::Progress)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn chunk(stream_idx: usize, chunk_serial: u64, data: &[u8]) -> FastqRawChunk {
        // The pairing reads only `chunk_serial` and `data`; ordinal/stream_idx
        // are carried by the readers but irrelevant to the pair join. Derive a
        // deterministic ordinal for clarity.
        let ordinal = chunk_serial * 2 + stream_idx as u64;
        FastqRawChunk { ordinal, stream_idx, chunk_serial, data: data.to_vec() }
    }

    #[test]
    fn pull_decision_pulls_both_when_not_backpressured() {
        assert_eq!(pull_decision(false, true, true), (true, true));
        assert_eq!(pull_decision(false, false, false), (true, true));
    }

    #[test]
    fn pull_decision_under_backpressure_pulls_only_the_behind_stream() {
        // The deadlock fix: when over the limit and the front pair is missing
        // B, we MUST still pull B (the behind stream) so the pair can complete
        // and the buffer drain — never refuse the completing pull.
        assert_eq!(pull_decision(true, false, true), (false, true), "missing B → pull B only");
        assert_eq!(pull_decision(true, true, false), (true, false), "missing A → pull A only");
    }

    #[test]
    fn pull_decision_under_backpressure_with_complete_front_pulls_neither() {
        // Front pair is complete; draining (emit) is the way forward, not pulling
        // more — so throttle both streams.
        assert_eq!(pull_decision(true, false, false), (false, false));
    }

    #[test]
    fn paired_batch_heap_size_and_ordinal() {
        let batch = PairedRawFastqBatch {
            ordinal: 9,
            chunk_serial: 3,
            data_a: vec![0u8; 10],
            data_b: vec![0u8; 7],
        };
        assert_eq!(batch.heap_size(), 17);
        assert_eq!(batch.ordinal(), 9);
    }

    #[test]
    fn profile_advertises_serial_fifo() {
        let step = PairRawFastq::new(1024);
        let p = step.profile();
        assert_eq!(p.name, "PairRawFastq");
        assert_eq!(p.kind, StepKind::Serial);
        assert!(!p.sticky);
        assert_eq!(step.affinity(), Affinity::None);
        // FIFO on this edge: the items carry a dense ordinal and the Parallel
        // downstream reorders, so a reorder stage here would needlessly
        // serialize the consumer.
        assert_eq!(p.branch_ordering, vec![BranchOrdering::None]);
        assert!(matches!(p.output_queues[0], QueueSpec::ByteBounded { .. }));
    }

    /// X5-004: a duplicate chunk for an already-filled `(chunk_serial, stream)`
    /// slot violates the single-chunk-per-serial invariant the byte accounting
    /// depends on. The hard `assert!` in `buffer_chunk` must catch it in all
    /// builds rather than letting `pending_total_bytes` silently drift.
    #[test]
    #[should_panic(expected = "single-chunk-per-serial invariant violated")]
    fn buffer_chunk_rejects_duplicate_slot() {
        let mut step = PairRawFastq::new(64 * 1024);
        step.buffer_chunk(chunk(0, 0, b"raw_a0"), true);
        // Second chunk for the same (serial 0, stream A) slot — must trip the guard.
        step.buffer_chunk(chunk(0, 0, b"raw_a0_again"), true);
    }

    /// Buffering a-then-b for a serial completes the pair; out-of-order serial
    /// arrival still emits in `chunk_serial` order with dense ordinals.
    #[test]
    fn pairs_in_chunk_serial_order_with_fresh_dense_ordinals() {
        let mut step = PairRawFastq::new(64 * 1024 * 1024);

        // Buffer serial 0 (A then B) and serial 1 (B then A) interleaved.
        step.buffer_chunk(chunk(0, 0, b"raw_a0"), true);
        assert!(step.try_emit_complete().is_none(), "serial 0 missing B");
        step.buffer_chunk(chunk(1, 1, b"raw_b1"), false);
        assert!(step.try_emit_complete().is_none(), "lowest serial 0 still missing B");
        step.buffer_chunk(chunk(1, 0, b"raw_b0"), false);

        // Now serial 0 is complete and is the lowest → emits first, ordinal 0.
        let p0 = step.try_emit_complete().unwrap();
        assert_eq!(p0.ordinal, 0);
        assert_eq!(p0.chunk_serial, 0);
        assert_eq!(p0.data_a, b"raw_a0");
        assert_eq!(p0.data_b, b"raw_b0");

        // serial 1 still missing A.
        assert!(step.try_emit_complete().is_none());
        step.buffer_chunk(chunk(0, 1, b"raw_a1"), true);
        let p1 = step.try_emit_complete().unwrap();
        assert_eq!(p1.ordinal, 1, "ordinals are dense and monotonic");
        assert_eq!(p1.chunk_serial, 1);
        assert_eq!(p1.data_a, b"raw_a1");
        assert_eq!(p1.data_b, b"raw_b1");

        assert_eq!(step.pending_total_bytes, 0);
    }

    // ─────────────────────────────────────────────────────────────────────────
    // X3-003: drive the paired-FASTQ source through a REAL fused chain under
    // back-pressure with skewed stream completion, on a watchdog thread.
    //
    // The in-module `pull_decision` tests above only check the decision
    // predicate in isolation. This test exercises the actual `PairRawFastq`
    // step inside a `Pipeline` with a tiny byte limit and one stream (R1)
    // racing far ahead of the other (R2): R1 emits ALL its chunks before R2
    // emits any. The partial-pair buffer crosses the limit while the front
    // pair is missing R2; the deadlock fix must keep pulling R2 (the behind
    // stream) so the front pair completes and the buffer drains. A regression
    // that globally refuses to pull under back-pressure wedges the step, and
    // the watchdog fires.
    // ─────────────────────────────────────────────────────────────────────────

    use std::sync::Arc;
    use std::sync::Mutex;
    use std::sync::atomic::{AtomicUsize, Ordering as AtomicOrd};

    use crate::pipeline::core::{PipelineBuilder, PipelineConfig, Step, StepCtx};

    /// Bytes per emitted raw chunk. Sized so that buffering a SINGLE unmatched
    /// leader chunk already crosses the pairing step's soft byte limit — that
    /// way the partial-pair buffer is over-limit with the front pair missing
    /// its R2 slot, which is exactly the state where a "refuse all pulls under
    /// back-pressure" regression deadlocks (it will not pull the available R2
    /// chunk that would complete the front pair and drain the buffer).
    const X3_CHUNK_BYTES: usize = 1024;
    /// Number of chunk serials each stream emits.
    const X3_N_SERIALS: u64 = 16;
    /// Soft backpressure limit for the pairing step under test. Sized so a
    /// single buffered unmatched chunk crosses the SOFT limit (1x) but stays
    /// under the HARD desync limit (2x): `768 < 1024 (one chunk) < 1536`. This
    /// puts the step in the over-soft-limit / front-pair-incomplete regime that
    /// the deadlock fix targets, WITHOUT tripping the catastrophic-desync abort.
    const X3_SOFT_LIMIT_BYTES: usize = 768;
    /// Per-source output-queue byte limit. Small (a few chunks) so a
    /// back-pressured pairing step leaves a source's output queue full and that
    /// source cannot drain — the precondition for the deadlock to manifest (if
    /// the streams always drained, `finalize_pairs` would flush everything even
    /// under a refuse-to-pull regression, masking the bug).
    const X3_SOURCE_QUEUE_BYTES: u64 = 4 * 1024;
    /// R2 holds back until R1 has emitted this many chunks, so the front pair
    /// (serial 0) is reliably buffered R1-first and over the limit with its R2
    /// slot still missing. Kept minimal (just enough to order the front pair) so
    /// it does not also stall the FIXED path, which over the limit pulls ONLY
    /// the behind stream.
    const X3_LAG_UNTIL: usize = 1;

    /// Serial source emitting `FastqRawChunk`s for ONE stream, serials `0..N`,
    /// at the given `stream_idx`. Mirrors the per-stream `ReadFastqInputs`
    /// reader: one chunk per dispatch, each carrying a globally-unique ordinal.
    ///
    /// To force the skewed-completion regime the deadlock fix targets, the
    /// stream is optionally gated on a shared "leader progress" counter: the
    /// leader (R1) increments it on every emit, and the laggard (R2) refuses to
    /// emit (`NoProgress`) until the leader has raced `lag_until` chunks ahead.
    /// This piles the leader's unmatched chunks in `PairRawFastq`'s partial-pair
    /// buffer past the byte limit while the front pair is still missing R2 —
    /// exactly the state where a "refuse all pulls under back-pressure"
    /// regression deadlocks (the completing R2 chunk is available but never
    /// pulled).
    struct OneStreamSource {
        stream_idx: usize,
        next_serial: u64,
        n_serials: u64,
        held: HeldSlot<Unpushed<FastqRawChunk>>,
        output_byte_limit: u64,
        /// Shared leader-progress counter (chunks the leader has emitted).
        leader_progress: Arc<AtomicUsize>,
        /// If `Some(n)`, this is the laggard: do not emit until the leader has
        /// emitted at least `n` chunks. If `None`, this is the leader.
        lag_until: Option<usize>,
    }

    impl OneStreamSource {
        fn leader(
            stream_idx: usize,
            n_serials: u64,
            output_byte_limit: u64,
            leader_progress: Arc<AtomicUsize>,
        ) -> Self {
            Self {
                stream_idx,
                next_serial: 0,
                n_serials,
                held: HeldSlot::new(),
                output_byte_limit,
                leader_progress,
                lag_until: None,
            }
        }

        fn laggard(
            stream_idx: usize,
            n_serials: u64,
            output_byte_limit: u64,
            leader_progress: Arc<AtomicUsize>,
            lag_until: usize,
        ) -> Self {
            Self {
                stream_idx,
                next_serial: 0,
                n_serials,
                held: HeldSlot::new(),
                output_byte_limit,
                leader_progress,
                lag_until: Some(lag_until),
            }
        }
    }

    impl Step for OneStreamSource {
        type Input = ();
        type Outputs = OrderedBytesSingle<FastqRawChunk>;

        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "OneStreamSource",
                kind: StepKind::Serial,
                sticky: false,
                output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
                branch_ordering: vec![BranchOrdering::ByItemOrdinal],
            }
        }

        fn affinity(&self) -> Affinity {
            // Distinct workers so R1 and R2 can race independently (mirrors the
            // production per-stream reader affinities).
            Affinity::Worker(self.stream_idx)
        }

        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            if let Some(unpushed) = self.held.take() {
                match ctx.outputs.retry(unpushed) {
                    Ok(()) => {}
                    Err(again) => {
                        self.held.put(again);
                        return Ok(StepOutcome::Contention);
                    }
                }
            }

            if self.next_serial >= self.n_serials {
                return Ok(StepOutcome::Finished);
            }

            // Laggard: stall until the leader has raced far enough ahead to
            // pile its chunks past the byte limit with the front pair missing.
            if let Some(lag_until) = self.lag_until {
                if self.leader_progress.load(AtomicOrd::Acquire) < lag_until {
                    return Ok(StepOutcome::NoProgress);
                }
            }

            let chunk_serial = self.next_serial;
            self.next_serial += 1;

            // Globally-unique ordinal across both streams (N == 2):
            // serial*2 + stream_idx.
            let ordinal = chunk_serial * 2 + self.stream_idx as u64;
            let chunk = FastqRawChunk {
                ordinal,
                stream_idx: self.stream_idx,
                chunk_serial,
                data: vec![b'A'; X3_CHUNK_BYTES],
            };
            match ctx.outputs.push(chunk) {
                Ok(()) => {
                    if self.lag_until.is_none() {
                        self.leader_progress.fetch_add(1, AtomicOrd::Release);
                    }
                    Ok(StepOutcome::Progress)
                }
                Err(unpushed) => {
                    self.held.put(unpushed);
                    if self.lag_until.is_none() {
                        self.leader_progress.fetch_add(1, AtomicOrd::Release);
                    }
                    Ok(StepOutcome::Progress)
                }
            }
        }
    }

    /// Serial sink that records each paired batch's `chunk_serial` in the order
    /// it is popped. Declared `Serial` (not `Parallel`) so the recorded order is
    /// exactly the emission order of the upstream `Serial` `PairRawFastq` step —
    /// a `Parallel` sink would let multiple workers race on `pop()` and reorder
    /// the observations, masking any ordering regression in the pairing step.
    #[derive(Clone)]
    struct PairSink {
        seen_serials: Arc<Mutex<Vec<u64>>>,
        count: Arc<AtomicUsize>,
    }

    impl Step for PairSink {
        type Input = PairedRawFastqBatch;
        type Outputs = ();

        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "PairSink",
                kind: StepKind::Serial,
                sticky: false,
                output_queues: vec![],
                branch_ordering: vec![],
            }
        }

        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            match ctx.input.pop() {
                Some(batch) => {
                    // Each emitted pair must carry equal-length R1/R2 data.
                    assert_eq!(batch.data_a.len(), X3_CHUNK_BYTES);
                    assert_eq!(batch.data_b.len(), X3_CHUNK_BYTES);
                    self.seen_serials.lock().unwrap().push(batch.chunk_serial);
                    self.count.fetch_add(1, AtomicOrd::Relaxed);
                    Ok(StepOutcome::Progress)
                }
                None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
                None => Ok(StepOutcome::NoProgress),
            }
        }

        fn new_worker_copy(&self) -> Self {
            self.clone()
        }
    }

    /// End-to-end: two per-stream sources (R1 races ahead, R2 lags) feed
    /// `PairRawFastq` (tiny byte limit) which feeds a collecting sink. Driven on
    /// a watchdog thread so the deadlock regression surfaces as a timeout. The
    /// run must complete and emit every pair in `chunk_serial` order.
    #[test]
    fn pair_fastq_does_not_deadlock_under_backpressure_with_skewed_streams() {
        use std::time::{Duration, Instant};

        let seen = Arc::new(Mutex::new(Vec::<u64>::new()));
        let count = Arc::new(AtomicUsize::new(0));

        // R1 (leader) races ahead; R2 (laggard) stalls until R1 has emitted
        // `X3_LAG_UNTIL` chunks, forcing the partial-pair buffer over the limit
        // with the front pair missing R2.
        let leader_progress = Arc::new(AtomicUsize::new(0));

        let seen_for_thread = Arc::clone(&seen);
        let count_for_thread = Arc::clone(&count);
        let leader_for_thread = Arc::clone(&leader_progress);
        let worker = std::thread::Builder::new()
            .name("pair-fastq-x3-003".into())
            .spawn(move || -> io::Result<()> {
                let r1 = OneStreamSource::leader(
                    0,
                    X3_N_SERIALS,
                    X3_SOURCE_QUEUE_BYTES,
                    Arc::clone(&leader_for_thread),
                );
                let r2 = OneStreamSource::laggard(
                    1,
                    X3_N_SERIALS,
                    X3_SOURCE_QUEUE_BYTES,
                    leader_for_thread,
                    X3_LAG_UNTIL,
                );
                let pair =
                    PairRawFastq::new(64 * 1024).with_backpressure_bytes(X3_SOFT_LIMIT_BYTES);
                let sink = PairSink { seen_serials: seen_for_thread, count: count_for_thread };

                let builder = PipelineBuilder::new();
                let r1_tail = builder.append_source(r1);
                let r2_tail = builder.append_source(r2);
                let pair_tail = builder.append_step2(pair, r1_tail, r2_tail);
                builder.append_step(sink, pair_tail);

                let pipeline = builder.build().map_err(io::Error::other)?;
                // 3 workers: R1, R2, and a free worker for pairing/sink.
                pipeline
                    .run(PipelineConfig { threads: 3, ..Default::default() })
                    .map_err(io::Error::other)
            })
            .expect("spawn pipeline worker");

        let deadline = Instant::now() + Duration::from_secs(30);
        while !worker.is_finished() {
            assert!(
                Instant::now() < deadline,
                "PairRawFastq pipeline did not finish within 30s — likely the X3-003 \
                 refuse-to-pull-under-backpressure deadlock regression"
            );
            std::thread::sleep(Duration::from_millis(25));
        }
        worker.join().expect("pipeline thread panicked").expect("pipeline run returned Err");

        // Full output: one pair per serial, every serial present AND observed in
        // ascending `chunk_serial` order. The Serial `PairSink` records in pop
        // order, so we assert the observed vector directly (no sorting) — an
        // out-of-order emission from `PairRawFastq` would now fail the test.
        let serials = seen.lock().unwrap().clone();
        assert_eq!(
            serials.len(),
            usize::try_from(X3_N_SERIALS).unwrap(),
            "expected one emitted pair per serial, got {}",
            serials.len()
        );
        let expected: Vec<u64> = (0..X3_N_SERIALS).collect();
        assert_eq!(
            serials, expected,
            "every chunk_serial must be paired exactly once and emitted in order"
        );
    }

    // ─────────────────────────────────────────────────────────────────────────
    // S5a2-011: when the two FASTQ streams have UNEQUAL record counts, one
    // stream's reader drains while the other still has a chunk for the lowest
    // pending serial. `finalize_pairs` must abort with an "out of sync" error
    // that names the stream index that ENDED EARLY (the missing one) — and it
    // must do so symmetrically, whichever of the two streams is short. The
    // `run_desync_harness` helper below drives both directions; the two tests
    // assert the stream-1-short and stream-0-short cases respectively.
    //
    // Wiring mirrors the X3-003 harness above: two `OneStreamSource` leaders
    // (no lag — we want clean end-of-stream desync at drain, not mid-stream
    // backpressure), one feeding 3 serials and one feeding 2. Serial 2 buffers
    // only the long stream's chunk; both readers then drain, `finalize_pairs`
    // finds serial 2 holding one side present and the other missing, computes
    // `short` from which side is absent, and reports "stream <short> ended
    // before chunk_serial 2". The backpressure limit is kept generous so the
    // catastrophic-desync hard-abort (2x) never fires first — the one-serial
    // skew keeps the partial-pair buffer tiny.
    // ─────────────────────────────────────────────────────────────────────────

    /// The long stream emits this many serials.
    const DESYNC_LONG_SERIALS: u64 = 3;
    /// The short stream emits this many serials (it ends one serial early).
    const DESYNC_SHORT_SERIALS: u64 = 2;

    /// Run the two-leader desync harness with `serials0` chunks on stream 0 and
    /// `serials1` on stream 1, returning the pipeline run result. Whichever
    /// stream is short drains first, leaving the lowest pending serial holding
    /// only the long stream's chunk, so `finalize_pairs` must abort with an
    /// "out of sync" error naming the SHORT (missing) stream index. Shared by
    /// both desync tests so the stream-0-short and stream-1-short cases exercise
    /// byte-for-byte the same wiring (only the lengths swap).
    fn run_desync_harness(serials0: u64, serials1: u64) -> Result<(), String> {
        run_desync_harness_with_backpressure(serials0, serials1, None)
    }

    fn run_desync_harness_with_backpressure(
        serials0: u64,
        serials1: u64,
        backpressure_bytes: Option<usize>,
    ) -> Result<(), String> {
        use std::time::{Duration, Instant};

        // Both sources are leaders with no lag gate; the shared counter is
        // unused for ordering here but required by the constructor signature.
        let unused_progress = Arc::new(AtomicUsize::new(0));

        let progress_for_thread = Arc::clone(&unused_progress);
        let worker = std::thread::Builder::new()
            .name("pair-fastq-s5a2-011".into())
            .spawn(move || -> Result<(), String> {
                let stream0 = OneStreamSource::leader(
                    0,
                    serials0,
                    X3_SOURCE_QUEUE_BYTES,
                    Arc::clone(&progress_for_thread),
                );
                let stream1 = OneStreamSource::leader(
                    1,
                    serials1,
                    X3_SOURCE_QUEUE_BYTES,
                    progress_for_thread,
                );
                // Default (generous) backpressure lets the streams reach the
                // end-of-stream `finalize_pairs` desync path; a small cap instead
                // exercises the mid-stream fail-fast (surplus over 1x stops
                // pulling before both branches drain).
                let pair = match backpressure_bytes {
                    Some(bytes) => PairRawFastq::new(64 * 1024).with_backpressure_bytes(bytes),
                    None => PairRawFastq::new(64 * 1024),
                };
                let sink = PairSink {
                    seen_serials: Arc::new(Mutex::new(Vec::new())),
                    count: Arc::new(AtomicUsize::new(0)),
                };

                let builder = PipelineBuilder::new();
                let tail0 = builder.append_source(stream0);
                let tail1 = builder.append_source(stream1);
                let pair_tail = builder.append_step2(pair, tail0, tail1);
                builder.append_step(sink, pair_tail);

                let pipeline = builder.build().map_err(|e| e.to_string())?;
                pipeline
                    .run(PipelineConfig { threads: 3, ..Default::default() })
                    .map_err(|e| e.to_string())
            })
            .expect("spawn pipeline worker");

        let deadline = Instant::now() + Duration::from_secs(30);
        while !worker.is_finished() {
            assert!(
                Instant::now() < deadline,
                "PairRawFastq desync pipeline did not finish within 30s — likely a deadlock"
            );
            std::thread::sleep(Duration::from_millis(25));
        }
        worker.join().expect("pipeline thread panicked")
    }

    #[test]
    fn finalize_pairs_reports_desync_with_short_stream_index() {
        // Stream 0 long (3 serials), stream 1 short (2): serial 2 buffers only
        // stream 0's chunk, so `short = usize::from(entry.0.is_some()) = 1` and
        // the error must name the MISSING (short) stream, index 1.
        let result = run_desync_harness(DESYNC_LONG_SERIALS, DESYNC_SHORT_SERIALS);
        let err = result.expect_err("unequal stream lengths must abort the pipeline run");
        assert!(
            err.contains("out of sync"),
            "desync abort must report an out-of-sync error, got: {err}"
        );
        assert!(
            err.contains("stream 1 ended before chunk_serial 2"),
            "desync error must name the short stream (index 1) and the orphaned \
             serial (2), got: {err}"
        );
    }

    #[test]
    fn finalize_pairs_reports_desync_when_stream_zero_is_short() {
        // Mirror of the above with the lengths swapped: stream 0 short (2
        // serials), stream 1 long (3). Serial 2 now buffers only stream 1's
        // chunk, so the missing-stream index is 0 — guards the symmetric branch
        // that a stream-1-only test would leave unexercised.
        let result = run_desync_harness(DESYNC_SHORT_SERIALS, DESYNC_LONG_SERIALS);
        let err = result.expect_err("unequal stream lengths must abort the pipeline run");
        assert!(
            err.contains("out of sync"),
            "desync abort must report an out-of-sync error, got: {err}"
        );
        assert!(
            err.contains("stream 0 ended before chunk_serial 2"),
            "desync error must name the short stream (index 0) and the orphaned \
             serial (2), got: {err}"
        );
    }

    #[test]
    fn truncated_stream_over_backpressure_fails_fast_instead_of_hanging() {
        // Stream 0 (a) has 3 chunks (X3_CHUNK_BYTES each), stream 1 (b) has 1.
        // With a small backpressure cap, after the first pair the surplus from
        // stream 0 crosses 1x the cap, so per-stream backpressure STOPS pulling
        // stream 0 — leaving its remaining chunks unpulled (so it never drains),
        // and pending plateaus BELOW the 2x catastrophic-desync abort, which
        // therefore never fires. The front pair is missing stream 1, which is
        // drained. Pre-fix this spun `NoProgress` forever (the harness's 30s
        // deadline would trip); the fix detects the unpairable record and fails
        // fast. This also proves the run terminates (no hang).
        let result = run_desync_harness_with_backpressure(3, 1, Some(X3_SOFT_LIMIT_BYTES));
        let err = result.expect_err("truncated stream over backpressure must abort, not hang");
        assert!(
            err.contains("out of sync") && err.contains("stream 1 ended before chunk_serial 1"),
            "backpressure fail-fast must reuse the shared out-of-sync error naming the \
             short stream (index 1) and the orphaned serial (1), got: {err}"
        );
    }

    #[test]
    fn truncated_stream_zero_over_backpressure_fails_fast_instead_of_hanging() {
        // Mirror of the above with the lengths swapped: stream 0 (a) has 1 chunk,
        // stream 1 (b) has 3. Now the surplus is stream 1, and the front pair is
        // missing stream 0, which is drained — so this drives the symmetric
        // `a_slot.is_none() && ctx.a.is_drained()` branch that the stream-1-short
        // test leaves unexercised. The abort must name the short stream (index 0).
        let result = run_desync_harness_with_backpressure(1, 3, Some(X3_SOFT_LIMIT_BYTES));
        let err = result.expect_err("truncated stream over backpressure must abort, not hang");
        assert!(
            err.contains("out of sync") && err.contains("stream 0 ended before chunk_serial 1"),
            "backpressure fail-fast must reuse the shared out-of-sync error naming the \
             short stream (index 0) and the orphaned serial (1), got: {err}"
        );
    }
}

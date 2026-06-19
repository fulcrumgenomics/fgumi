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

    /// Buffer a chunk from stream A (slot 0) or stream B (slot 1).
    fn buffer_chunk(&mut self, chunk: FastqRawChunk, is_a: bool) {
        self.pending_total_bytes += chunk.data.len();
        let entry = self.pending.entry(chunk.chunk_serial).or_insert((None, None));
        if is_a {
            entry.0 = Some(chunk);
        } else {
            entry.1 = Some(chunk);
        }
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

        // Nothing pulled or emitted this call. If both input branches are
        // drained, the pairing is complete — flush remaining pairs and finish.
        if ctx.a.is_drained() && ctx.b.is_drained() {
            return self.finalize_pairs(ctx);
        }
        Ok(StepOutcome::NoProgress)
    }
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
            return Err(io::Error::other(format!(
                "FASTQ sources out of sync: stream {short} ended before chunk_serial \
                 {lowest_serial} while the other stream had more records"
            )));
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
}

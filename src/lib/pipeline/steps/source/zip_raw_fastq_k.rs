//! `ZipRawFastqK` — the `Serial`, homogeneous K-input lockstep aligner for the
//! parallel FASTQ decode front.
//!
//! Sits at the convergence of K independent per-stream decode sub-chains (each
//! `ReadFastqBlocks(i) → FastqDecompress → FindFastqBoundaries(i)` for bgzip
//! input, or a single-stream `ReadFastqInputs` for gzip), each emitting
//! [`FastqRawChunk`](super::read_fastq::FastqRawChunk) tagged with its
//! `chunk_serial`. `ZipRawFastqK` pulls one chunk from every stream for the
//! current row (by `chunk_serial`), and once **all K** are present emits an
//! [`NRawFastqBatch`](super::fastq_zip::NRawFastqBatch) carrying the K aligned
//! raw byte-chunks plus a freshly-minted dense `ordinal`. The expensive parse
//! and template build then run in parallel in
//! [`ParseAndZipFastqN`](super::parse_zip_fastq::ParseAndZipFastqN).
//!
//! ## Lockstep, not merge
//!
//! This is a positional zip: row `N` = chunk `N` from every stream, joined by
//! `chunk_serial`. It is NOT a key-ordered merge, so there is no per-stream
//! ahead/behind pull arithmetic (contrast the former 2-input `PairRawFastq`'s
//! `pull_decision`). The backpressure rule is a one-liner: **pull only the
//! streams whose slot for the lowest pending row is still empty.** A stream that
//! already contributed its chunk to the front row is not pulled again until the
//! row completes and drains, so a fast stream cannot run ahead and balloon the
//! buffer past one row's worth of chunks.
//!
//! ## Ordinal
//!
//! `ZipRawFastqK` mints its own dense, monotonic `ordinal` (one per emitted
//! `NRawFastqBatch`), which the `Parallel` `ParseAndZipFastqN` reorders by. Its
//! own output edge is therefore FIFO (`BranchOrdering::None`): a reorder stage
//! here would serialize that parallel consumer (same rationale documented on
//! the former `PairRawFastq`).

use std::collections::BTreeMap;
use std::io;

use super::fastq_zip::NRawFastqBatch;
use super::read_fastq::FastqRawChunk;
use crate::pipeline::core::Unpushed;
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::outputs::OrderedBytesSingle;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{
    Affinity, Step, StepCtx, StepCtxK, StepK, StepKind, StepOutcome, StepProfile,
};

/// Backpressure cap on total bytes buffered across partial rows. Mirrors the
/// former `PairRawFastq`/`ZipFastqRecords` limit. With the lockstep pull rule a
/// single row of unmatched chunks is the steady-state buffer, so this only trips
/// on a genuine unbounded desync (unequal stream lengths).
const DEFAULT_PENDING_BACKPRESSURE_BYTES: usize = 256 * 1024 * 1024;

/// `Serial` homogeneous K-input lockstep aligner. Buffers, per `chunk_serial`,
/// the K per-stream chunks; emits an `NRawFastqBatch` once the lowest pending
/// serial has all K.
pub struct ZipRawFastqK {
    k: usize,
    /// Partial rows keyed by `chunk_serial`. Each value is a `Vec` of K slots
    /// (`slots[stream_idx]`); a row is emitted (and removed) once every slot is
    /// filled and it is the lowest pending serial.
    pending: BTreeMap<u64, Vec<Option<Vec<u8>>>>,
    pending_total_bytes: usize,
    pending_backpressure_bytes: usize,
    /// Self-minted dense ordinal for the downstream parallel reorder.
    next_ordinal: u64,
    held: HeldSlot<Unpushed<NRawFastqBatch>>,
    output_byte_limit: u64,
}

impl ZipRawFastqK {
    /// Construct a K-stream lockstep aligner. `k` must be >= 2 (a single stream
    /// needs no cross-stream zip; use the single-input parse path).
    ///
    /// # Panics
    ///
    /// Panics if `k < 2`; a single stream must use the single-input parse path
    /// (`WrapRawFastq1`) instead.
    #[must_use]
    pub fn new(k: usize, output_byte_limit: u64) -> Self {
        assert!(k >= 2, "ZipRawFastqK requires k >= 2, got {k}");
        Self {
            k,
            pending: BTreeMap::new(),
            pending_total_bytes: 0,
            pending_backpressure_bytes: DEFAULT_PENDING_BACKPRESSURE_BYTES,
            next_ordinal: 0,
            held: HeldSlot::new(),
            output_byte_limit,
        }
    }

    /// Override the soft backpressure limit (test-only).
    #[cfg(test)]
    fn with_backpressure_bytes(mut self, bytes: usize) -> Self {
        self.pending_backpressure_bytes = bytes;
        self
    }

    /// Buffer a chunk from `stream_idx` into its slot at `chunk_serial`. The
    /// slot must be empty (one chunk per `(chunk_serial, stream)` — the lockstep
    /// pull rule guarantees it); a duplicate is a wiring/replay bug that would
    /// silently change emitted identity, so it fails hard in all builds.
    fn buffer(&mut self, chunk: FastqRawChunk) {
        let k = self.k;
        self.pending_total_bytes += chunk.data.capacity();
        let row = self
            .pending
            .entry(chunk.chunk_serial)
            .or_insert_with(|| (0..k).map(|_| None).collect());
        let slot = &mut row[chunk.stream_idx];
        assert!(
            slot.is_none(),
            "duplicate chunk for serial {} stream {}: one-chunk-per-(serial,stream) invariant \
             violated",
            chunk.chunk_serial,
            chunk.stream_idx,
        );
        *slot = Some(chunk.data);
    }

    /// If the lowest pending serial has all K slots filled, remove it and build
    /// an `NRawFastqBatch` with a fresh ordinal.
    fn try_emit(&mut self) -> Option<NRawFastqBatch> {
        let lowest = *self.pending.keys().next()?;
        if !self.pending.get(&lowest)?.iter().all(Option::is_some) {
            return None;
        }
        let slots = self.pending.remove(&lowest).unwrap();
        let streams: Vec<Vec<u8>> = slots.into_iter().map(Option::unwrap).collect();
        self.pending_total_bytes =
            self.pending_total_bytes.saturating_sub(streams.iter().map(Vec::capacity).sum());
        let ordinal = self.next_ordinal;
        self.next_ordinal += 1;
        Some(NRawFastqBatch { ordinal, chunk_serial: lowest, streams })
    }

    /// Push a batch, holding it on rejection.
    fn emit(&mut self, batch: NRawFastqBatch, ctx: &mut StepCtxK<'_, Self>) {
        if let Err(unpushed) = ctx.outputs.push(batch) {
            self.held.put(unpushed);
        }
    }
}

impl StepK for ZipRawFastqK {
    type Input = FastqRawChunk;
    type Outputs = OrderedBytesSingle<NRawFastqBatch>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "ZipRawFastqK",
            kind: StepKind::Serial,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            // FIFO: we mint a dense ordinal that the Parallel ParseAndZipFastqN
            // reorders by; a reorder stage here would serialize that consumer.
            branch_ordering: vec![BranchOrdering::None],
        }
    }

    fn input_count(&self) -> usize {
        self.k
    }

    fn affinity(&self) -> Affinity {
        // Any free worker may drive the (cheap) alignment; pinning it to a
        // reader's worker would serialize read + zip on one thread.
        Affinity::None
    }

    fn try_run(&mut self, ctx: &mut StepCtxK<'_, Self>) -> io::Result<StepOutcome> {
        // 1. Drain a held emit first.
        if let Some(unpushed) = self.held.take() {
            match ctx.outputs.retry(unpushed) {
                Ok(()) => {}
                Err(again) => {
                    self.held.put(again);
                    return Ok(StepOutcome::Contention);
                }
            }
        }

        // 2. Catastrophic-desync guard: one stream producing far more than its
        // peers balloons the buffer. With the lockstep pull rule (step 3) this
        // should never trip on well-formed equal-length input; it is the safety
        // net for unbounded divergence.
        if self.pending_total_bytes > 2 * self.pending_backpressure_bytes {
            return Err(io::Error::other(format!(
                "ZipRawFastqK: catastrophic stream desync — pending buffer ({} bytes) exceeds 2x \
                 backpressure limit ({} bytes). One FASTQ stream is producing far more data than \
                 its peers.",
                self.pending_total_bytes, self.pending_backpressure_bytes
            )));
        }

        // 3. Lockstep pull + emit. Pull ONLY the streams missing from the lowest
        // pending row (the ones needed to complete it); a stream already present
        // at the front is not re-pulled, so it cannot run ahead. Loop so a
        // dispatch amortizes the Serial mutex and keeps the Parallel consumer fed.
        let mut did_work = false;
        loop {
            // Determine which streams are missing from the lowest pending row.
            // `None` pending row (empty buffer) means every stream is missing.
            let front_missing: Vec<bool> = match self.pending.values().next() {
                Some(row) => row.iter().map(Option::is_none).collect(),
                None => vec![true; self.k],
            };

            let mut pulled_this_iter = false;
            for (i, &missing) in front_missing.iter().enumerate() {
                if missing && let Some(chunk) = ctx.inputs[i].pop() {
                    self.buffer(chunk);
                    pulled_this_iter = true;
                    did_work = true;
                }
            }

            // Emit every complete row now at the front.
            while let Some(batch) = self.try_emit() {
                did_work = true;
                self.emit(batch, ctx);
                if self.held.is_held() {
                    // Output full — yield so the consumer drains; the held-slot
                    // preamble retries next dispatch.
                    return Ok(StepOutcome::Progress);
                }
            }

            if !pulled_this_iter {
                break;
            }
        }

        if did_work {
            return Ok(StepOutcome::Progress);
        }

        // 4. No progress this call. If every input is drained, the pairing is
        // complete: flush any remaining complete row, else (a partial row with a
        // drained stream) report the unequal-length desync; once empty, Finished.
        if ctx.inputs.iter().all(|h| h.is_drained()) {
            if let Some(batch) = self.try_emit() {
                self.emit(batch, ctx);
                return Ok(StepOutcome::Progress);
            }
            if let Some((&serial, row)) = self.pending.iter().next() {
                // Some stream ended early: name a missing stream and a present
                // one, matching the fgbio-consistent "out of sync" wording.
                let short = row.iter().position(Option::is_none).unwrap_or(0);
                let present = row.iter().position(Option::is_some).map_or(usize::MAX, |p| p);
                return Err(io::Error::other(format!(
                    "FASTQ sources out of sync: R{} ended before R{} at chunk_serial {serial} \
                     while other streams had more records",
                    short + 1,
                    present.wrapping_add(1),
                )));
            }
            return Ok(StepOutcome::Finished);
        }

        // 5. Mid-stream fail-fast: not every input is drained, but a stream that
        // the front (lowest pending) row is still WAITING on has already drained.
        // That row can never complete (the drained stream will never produce its
        // chunk), and the streams that ARE still live are held at the front —
        // ZipRawFastqK does not re-pull a stream already present in the front
        // row — so their queues never drain and the all-drained arm above is
        // unreachable. This is the unequal-length desync (a short stream ended
        // early). Without this check the step spins `NoProgress` forever — a
        // silent hang. Mirrors `finalize`'s wording so both fail paths report
        // the same "out of sync" diagnostic. (A live stream missing from the
        // front row is a transient — its chunk is simply still upstream — so
        // only a *drained* missing stream trips this.)
        if let Some((&serial, row)) = self.pending.iter().next()
            && let Some(short) = row
                .iter()
                .enumerate()
                .find_map(|(i, slot)| (slot.is_none() && ctx.inputs[i].is_drained()).then_some(i))
        {
            let present = row.iter().position(Option::is_some).map_or(usize::MAX, |p| p);
            return Err(io::Error::other(format!(
                "FASTQ sources out of sync: R{} ended before R{} at chunk_serial {serial} \
                 while other streams had more records",
                short + 1,
                present.wrapping_add(1),
            )));
        }
        Ok(StepOutcome::NoProgress)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// WrapRawFastq1 — single-input adapter into the unified K-stream decode front.
// ─────────────────────────────────────────────────────────────────────────────

/// `Serial` single-input adapter that wraps each [`FastqRawChunk`] from one
/// stream into a one-element [`NRawFastqBatch`], minting its own dense ordinal.
///
/// [`ZipRawFastqK`] requires `k >= 2` (a lone stream needs no cross-stream
/// alignment), so the K == 1 case cannot use it. This step is the degenerate
/// `k == 1` shim: it puts the single stream's bytes into `streams[0]` and mints
/// a dense, monotonic `ordinal` so the downstream `Parallel`
/// [`ParseAndZipFastqN`](super::parse_zip_fastq::ParseAndZipFastqN) reorders by
/// it — exactly the ordinal contract `ZipRawFastqK` provides for `k >= 2`.
///
/// FIFO output (`BranchOrdering::None`) for the same reason as `ZipRawFastqK`:
/// the emitted batch already carries a dense ordinal that the Parallel consumer
/// reorders by, so a reorder stage here would needlessly serialize it.
pub struct WrapRawFastq1 {
    next_ordinal: u64,
    held: HeldSlot<Unpushed<NRawFastqBatch>>,
    output_byte_limit: u64,
}

impl WrapRawFastq1 {
    #[must_use]
    pub fn new(output_byte_limit: u64) -> Self {
        Self { next_ordinal: 0, held: HeldSlot::new(), output_byte_limit }
    }
}

impl Step for WrapRawFastq1 {
    type Input = FastqRawChunk;
    type Outputs = OrderedBytesSingle<NRawFastqBatch>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "WrapRawFastq1",
            kind: StepKind::Serial,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            // FIFO: we mint a dense ordinal that the Parallel ParseAndZipFastqN
            // reorders by; a reorder stage here would serialize that consumer.
            branch_ordering: vec![BranchOrdering::None],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        // Drain a held emit first (single-occupancy held slot).
        if let Some(unpushed) = self.held.take() {
            match ctx.outputs.retry(unpushed) {
                Ok(()) => {}
                Err(again) => {
                    self.held.put(again);
                    return Ok(StepOutcome::Contention);
                }
            }
        }

        let Some(chunk) = ctx.input.pop() else {
            if ctx.input.is_drained() {
                return Ok(StepOutcome::Finished);
            }
            return Ok(StepOutcome::NoProgress);
        };

        let ordinal = self.next_ordinal;
        self.next_ordinal += 1;
        let batch =
            NRawFastqBatch { ordinal, chunk_serial: chunk.chunk_serial, streams: vec![chunk.data] };

        match ctx.outputs.push(batch) {
            Ok(()) => Ok(StepOutcome::Progress),
            Err(unpushed) => {
                self.held.put(unpushed);
                Ok(StepOutcome::Progress)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pipeline::core::item::Ordered;

    fn chunk(stream_idx: usize, chunk_serial: u64, data: &[u8]) -> FastqRawChunk {
        FastqRawChunk { ordinal: 0, stream_idx, chunk_serial, data: data.to_vec() }
    }

    #[test]
    fn profile_is_serial_fifo() {
        let s = ZipRawFastqK::new(3, 1 << 20);
        let p = s.profile();
        assert_eq!(p.name, "ZipRawFastqK");
        assert_eq!(p.kind, StepKind::Serial);
        assert_eq!(p.branch_ordering, vec![BranchOrdering::None]);
        assert_eq!(s.input_count(), 3);
    }

    /// A row completes only when all K slots are filled, and emits in
    /// `chunk_serial` order with dense ordinals.
    #[test]
    fn emits_when_row_complete_with_dense_ordinals() {
        let mut s = ZipRawFastqK::new(3, 1 << 20);
        // Serial 0: fill streams 0 and 2, still missing 1.
        s.buffer(chunk(0, 0, b"a0"));
        s.buffer(chunk(2, 0, b"c0"));
        assert!(s.try_emit().is_none(), "serial 0 missing stream 1");
        // Serial 1 arrives before serial 0 completes — must not jump the queue.
        s.buffer(chunk(0, 1, b"a1"));
        assert!(s.try_emit().is_none(), "lowest serial 0 still incomplete");
        // Complete serial 0.
        s.buffer(chunk(1, 0, b"b0"));
        let r0 = s.try_emit().expect("serial 0 complete");
        assert_eq!(r0.ordinal, 0);
        assert_eq!(r0.chunk_serial, 0);
        assert_eq!(r0.streams, vec![b"a0".to_vec(), b"b0".to_vec(), b"c0".to_vec()]);
        // Serial 1 still missing 1 and 2.
        assert!(s.try_emit().is_none());
        s.buffer(chunk(1, 1, b"b1"));
        s.buffer(chunk(2, 1, b"c1"));
        let r1 = s.try_emit().expect("serial 1 complete");
        assert_eq!(r1.ordinal, 1, "dense monotonic ordinal");
        assert_eq!(r1.chunk_serial, 1);
        assert_eq!(s.pending_total_bytes, 0, "buffer drains to zero");
    }

    #[test]
    #[should_panic(expected = "one-chunk-per-(serial,stream) invariant violated")]
    fn duplicate_slot_panics() {
        let mut s = ZipRawFastqK::new(2, 1 << 20);
        s.buffer(chunk(0, 0, b"a"));
        s.buffer(chunk(0, 0, b"a-again"));
    }

    /// A finalize with a partial front row (one stream ran short) reports the
    /// fgbio-consistent "out of sync" desync rather than emitting or hanging.
    /// (Drives `try_emit`/`pending` directly; the full drained-input path is
    /// covered by the chain-level `mismatched_stream_lengths_*` test.)
    #[test]
    fn partial_front_row_is_out_of_sync() {
        let mut s = ZipRawFastqK::new(2, 1 << 20);
        // Stream 0 produced serial 0; stream 1 never did.
        s.buffer(chunk(0, 0, b"a0"));
        assert!(s.try_emit().is_none(), "incomplete row must not emit");
        // The row is genuinely stuck: stream 0 present, stream 1 missing.
        let (_, row) = s.pending.iter().next().unwrap();
        assert!(row[0].is_some() && row[1].is_none(), "front row waits on stream 1");
    }

    /// The test-only backpressure override lowers the soft limit so the
    /// pending-byte accounting can be exercised with small chunks: buffering a
    /// single stream's chunk grows `pending_total_bytes`, and it drains back to
    /// zero once the row completes and emits.
    #[test]
    fn with_backpressure_bytes_tracks_pending_total() {
        let mut s = ZipRawFastqK::new(2, 1 << 20).with_backpressure_bytes(8);
        s.buffer(chunk(0, 0, b"aaaa"));
        assert!(s.pending_total_bytes >= 4, "buffered bytes are accounted");
        assert!(s.try_emit().is_none(), "row incomplete until both streams present");
        s.buffer(chunk(1, 0, b"bbbb"));
        let row = s.try_emit().expect("row complete");
        assert_eq!(row.chunk_serial, 0);
        assert_eq!(s.pending_total_bytes, 0, "buffer drains to zero after emit");
    }

    #[test]
    fn wrap_raw_fastq1_profile_is_serial_fifo() {
        let s = WrapRawFastq1::new(1 << 20);
        let p = s.profile();
        assert_eq!(p.name, "WrapRawFastq1");
        assert_eq!(p.kind, StepKind::Serial);
        assert_eq!(p.branch_ordering, vec![BranchOrdering::None]);
    }

    /// Wrapping preserves the chunk's `chunk_serial` and bytes into a
    /// single-stream `NRawFastqBatch`, minting a fresh dense ordinal.
    #[test]
    fn wrap_raw_fastq1_wraps_with_dense_ordinal() {
        let mut s = WrapRawFastq1::new(1 << 20);
        // Simulate two successive wraps by driving `next_ordinal` the way
        // `try_run` does (the StepCtx plumbing is exercised by integration tests).
        let c0 = chunk(0, 3, b"row0");
        let ord0 = s.next_ordinal;
        s.next_ordinal += 1;
        let b0 =
            NRawFastqBatch { ordinal: ord0, chunk_serial: c0.chunk_serial, streams: vec![c0.data] };
        assert_eq!(b0.ordinal(), 0);
        assert_eq!(b0.chunk_serial, 3);
        assert_eq!(b0.streams, vec![b"row0".to_vec()]);

        let c1 = chunk(0, 4, b"row1");
        let ord1 = s.next_ordinal;
        s.next_ordinal += 1;
        let b1 =
            NRawFastqBatch { ordinal: ord1, chunk_serial: c1.chunk_serial, streams: vec![c1.data] };
        assert_eq!(b1.ordinal(), 1, "dense monotonic ordinal");
        assert_eq!(b1.streams.len(), 1);
    }
}

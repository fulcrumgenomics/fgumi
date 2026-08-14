//! `ZipFastqRecords` step: joins per-stream `FastqChunkBatch` items into
//! `FastqTemplateBatch` items by matching chunk serials across streams.
//!
//! `Serial` + `Affinity::None`. Receives round-robin chunks from
//! `ReadFastqInputs` and buffers them in a `BTreeMap` keyed by
//! `chunk_serial`. When all `n_streams` slots for the lowest serial are
//! filled, the records are zipped into `FastqTemplate`s and emitted.

use std::collections::BTreeMap;
use std::io;

use crate::fastq_parse::{FastqRecord, strip_read_suffix};
use crate::grouper::FastqTemplate;
use crate::pipeline::core::Unpushed;
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::item::{HeapSize, Ordered};
use crate::pipeline::core::outputs::OrderedBytesSingle;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};

use super::read_fastq::FastqChunkBatch;

// ─────────────────────────────────────────────────────────────────────────────
// FastqTemplateBatch
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug)]
pub struct FastqTemplateBatch {
    pub batch_serial: u64,
    pub templates: Vec<FastqTemplate>,
    total_bytes: usize,
}

impl FastqTemplateBatch {
    #[must_use]
    pub fn new(batch_serial: u64, templates: Vec<FastqTemplate>) -> Self {
        let total_bytes: usize = templates.iter().map(HeapSize::heap_size).sum();
        Self { batch_serial, templates, total_bytes }
    }
}

impl HeapSize for FastqTemplateBatch {
    fn heap_size(&self) -> usize {
        self.total_bytes
    }
}

impl Ordered for FastqTemplateBatch {
    fn ordinal(&self) -> u64 {
        self.batch_serial
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ZipFastqRecords — Serial step
// ─────────────────────────────────────────────────────────────────────────────

const DEFAULT_PENDING_BACKPRESSURE_BYTES: usize = 256 * 1024 * 1024;

/// Outcome of a full egress drain attempt ([`ZipFastqRecords::drain_complete`]).
enum DrainStop {
    /// The lowest serial is incomplete (or `pending` is empty): nothing more to
    /// emit right now. `emitted` records whether at least one batch left.
    Exhausted { emitted: bool },
    /// A downstream push could not complete; the unpushed batch is now in
    /// `self.held`. The caller must stop draining immediately (no further serial
    /// may be emitted while the single held slot is occupied) and return
    /// `Progress` — parking a complete serial in `held` is real forward progress,
    /// and the held-slot retry preamble on the next dispatch flushes it before
    /// any more work runs (see the `try_run` call site).
    Held,
}

pub struct ZipFastqRecords {
    n_streams: usize,
    pending: BTreeMap<u64, Vec<Option<FastqChunkBatch>>>,
    pending_total_bytes: usize,
    pending_backpressure_bytes: usize,
    held: HeldSlot<Unpushed<FastqTemplateBatch>>,
    output_byte_limit: u64,
    next_batch_serial: u64,
    /// Whether the soft-limit warning has already fired for the current
    /// over-limit episode. This step is `Serial` and the driver re-dispatches it
    /// in a tight loop, so an unguarded warning emits one identical line per
    /// dispatch for the whole time a legitimately-ahead stream stays over the
    /// limit — the exact case the branch exists to admit. Cleared when the
    /// buffer drops back under the limit, so each episode warns once.
    over_soft_limit_warned: bool,
}

impl ZipFastqRecords {
    /// Construct a zipper for `n_streams` FASTQ streams.
    ///
    /// # Preconditions
    ///
    /// `n_streams >= 1`. The count is wired from the FASTQ-input count, which
    /// is always at least one, and `try_emit_complete` indexes
    /// `stream_records[0]`, so a zero-stream configuration would panic on that
    /// index. A `debug_assert!` makes an accidental zero-stream wiring fail
    /// loudly in debug/test builds.
    #[must_use]
    pub fn new(n_streams: usize, output_byte_limit: u64) -> Self {
        debug_assert!(n_streams >= 1, "ZipFastqRecords requires n_streams >= 1, got {n_streams}");
        Self {
            n_streams,
            pending: BTreeMap::new(),
            pending_total_bytes: 0,
            pending_backpressure_bytes: DEFAULT_PENDING_BACKPRESSURE_BYTES,
            held: HeldSlot::new(),
            output_byte_limit,
            next_batch_serial: 0,
            over_soft_limit_warned: false,
        }
    }

    /// Latch the over-soft-limit episode, returning `true` only on the dispatch
    /// that *enters* it.
    ///
    /// The caller warns on `true`. Because this step is `Serial` and the driver
    /// re-dispatches it in a tight loop, warning unconditionally would emit one
    /// identical line per dispatch for as long as a legitimately-ahead stream
    /// stays over the limit.
    fn enter_over_soft_limit(&mut self) -> bool {
        let first = !self.over_soft_limit_warned;
        self.over_soft_limit_warned = true;
        first
    }

    /// Clear the episode latch, so the next crossing warns again.
    fn leave_over_soft_limit(&mut self) {
        self.over_soft_limit_warned = false;
    }

    /// Override the soft pending-backpressure limit (test-only).
    ///
    /// Production wires the default 256 MiB limit; tests shrink it to a few KiB
    /// so a lag-then-burst scenario can cross the soft/hard thresholds with a
    /// handful of small records rather than hundreds of megabytes of data.
    #[cfg(test)]
    fn with_backpressure_bytes(mut self, bytes: usize) -> Self {
        self.pending_backpressure_bytes = bytes;
        self
    }

    /// Egress: emit every currently-complete consecutive lowest serial.
    ///
    /// Loops [`Self::try_emit_complete`] until it yields `None` (the lowest
    /// serial is incomplete or `pending` is empty) or a downstream push fails.
    /// On a failed push the unpushed batch is stashed in `self.held` and the
    /// loop stops — the single held slot can hold only one item, so no further
    /// serial may be emitted until it flushes on a later dispatch.
    ///
    /// Backpressure never gates this: draining is the only action that relieves
    /// pending-byte pressure, so it must always run. Bounded: each successful
    /// emit removes one key from `pending` (in `try_emit_complete`), so the loop
    /// runs at most `pending.len()` times and strictly shrinks `pending`.
    ///
    /// Precondition: `self.held` is empty (the held-slot preamble in `try_run`
    /// guarantees this by returning `Contention` on a failed retry).
    fn drain_complete(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<DrainStop> {
        debug_assert!(
            !self.held.is_held(),
            "drain_complete requires an empty held slot — the try_run preamble must clear it"
        );
        let mut emitted = false;
        loop {
            match self.try_emit_complete()? {
                None => return Ok(DrainStop::Exhausted { emitted }),
                Some(template_batch) => match ctx.outputs.push(template_batch) {
                    Ok(()) => emitted = true,
                    Err(unpushed) => {
                        self.held.put(unpushed);
                        return Ok(DrainStop::Held);
                    }
                },
            }
        }
    }

    fn try_emit_complete(&mut self) -> io::Result<Option<FastqTemplateBatch>> {
        let Some(lowest_serial) = self.pending.keys().next().copied() else {
            return Ok(None);
        };

        let slots = self.pending.get(&lowest_serial).unwrap();
        let all_filled = slots.iter().all(Option::is_some);
        if !all_filled {
            return Ok(None);
        }

        let chunks: Vec<FastqChunkBatch> = self
            .pending
            .remove(&lowest_serial)
            .unwrap()
            .into_iter()
            .map(|opt| opt.expect("every slot is filled — checked by `all_filled` above"))
            .collect();

        // Subtract using each batch's own `HeapSize`, which is the measure the
        // byte-bounded queues budget and the same one the insert side adds. A
        // hand-rolled sum over `records` would be a *second* measure of the same
        // chunk: it happens to agree today because the records are moved through
        // unchanged, but `FastqChunkBatch::new` takes `total_bytes` as a
        // parameter rather than deriving it, so the two can diverge — and then
        // `pending_total_bytes` would drift against the queue budget and the
        // soft/hard thresholds would fire on a different measure than the queues
        // enforce.
        let chunk_bytes: usize = chunks.iter().map(HeapSize::heap_size).sum();
        self.pending_total_bytes = self.pending_total_bytes.saturating_sub(chunk_bytes);

        let mut stream_records: Vec<Vec<FastqRecord>> =
            chunks.into_iter().map(|chunk| chunk.records).collect();

        let n_records = stream_records[0].len();

        // Paired/N-stream FASTQ is positionally parallel: every stream must
        // contribute the same number of records for this chunk. Reject a
        // mismatch explicitly (matches the N==2 ParseAndZipFastq path and
        // fgbio's FastqSource.zipped "out of sync" check) rather than letting
        // it surface as an incidental name mismatch from the back-pop below.
        for (stream_idx, recs) in stream_records.iter().enumerate().skip(1) {
            if recs.len() != n_records {
                return Err(io::Error::other(format!(
                    "FASTQ sources out of sync at chunk_serial {lowest_serial}: \
                     stream 0 has {n_records} records, stream {stream_idx} has {}",
                    recs.len(),
                )));
            }
        }

        // Zip across streams: record i from each stream forms one template.
        // Drain from the back (via pop) to avoid O(n^2) from front removal;
        // reverse in place afterwards.
        let mut templates = Vec::with_capacity(n_records);
        for _ in 0..n_records {
            let mut records = Vec::with_capacity(self.n_streams);
            let mut base_name: Option<Vec<u8>> = None;

            for (stream_idx, stream) in stream_records.iter_mut().enumerate() {
                // Unreachable after the equal-count check above; kept defensive.
                let Some(record) = stream.pop() else {
                    return Err(io::Error::other(format!(
                        "FASTQ sources out of sync: stream {stream_idx} ran short at chunk_serial {lowest_serial}",
                    )));
                };

                let stripped = strip_read_suffix(record.name());
                match &base_name {
                    None => {
                        base_name = Some(stripped.to_vec());
                    }
                    Some(expected) => {
                        if stripped != expected.as_slice() {
                            return Err(io::Error::new(
                                io::ErrorKind::InvalidData,
                                format!(
                                    "FASTQ read name mismatch at chunk_serial {}: stream 0 has '{}', stream {} has '{}'",
                                    lowest_serial,
                                    String::from_utf8_lossy(expected),
                                    stream_idx,
                                    String::from_utf8_lossy(stripped),
                                ),
                            ));
                        }
                    }
                }
                records.push(record);
            }

            templates.push(FastqTemplate { name: base_name.unwrap(), records });
        }

        // We popped from the back, so reverse to restore original order.
        templates.reverse();

        let serial = self.next_batch_serial;
        self.next_batch_serial += 1;
        Ok(Some(FastqTemplateBatch::new(serial, templates)))
    }
}

impl Step for ZipFastqRecords {
    type Input = FastqChunkBatch;
    type Outputs = OrderedBytesSingle<FastqTemplateBatch>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "ZipFastqRecords",
            kind: StepKind::Serial,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        // ── (a) Held-slot retry preamble ────────────────────────────────────
        // Single-occupancy: if the retry still can't push, re-hold and return
        // Contention so control never reaches egress/ingress with a full slot.
        if let Some(unpushed) = self.held.take() {
            match ctx.outputs.retry(unpushed) {
                Ok(()) => {}
                Err(again) => {
                    self.held.put(again);
                    return Ok(StepOutcome::Contention);
                }
            }
        }

        // ── (b) Full egress drain — ALWAYS runs, before any backpressure gate.
        // Draining the complete lowest serials is the ONLY action that relieves
        // pending-byte pressure, so it must never be throttled (the previous
        // code returned NoProgress/Err on backpressure *before* emitting, which
        // wedged the step whenever the lowest serial was already complete but
        // over the soft limit — see issue notes for S5a1-001).
        match self.drain_complete(ctx)? {
            // `Held` means a complete serial was removed from `pending` (bytes
            // subtracted) and parked in the held slot for retry — that is real
            // forward progress ("held one for later"), not mutex contention.
            // Misreporting it as `Contention` on this liveness-critical
            // backpressure path would understate progress to the driver.
            DrainStop::Held => return Ok(StepOutcome::Progress),
            DrainStop::Exhausted { emitted } => {
                // If anything was emitted, report Progress now: pressure dropped
                // and the next dispatch re-evaluates ingress against the new
                // (lower) byte total. We deliberately do NOT also pop input this
                // dispatch — that keeps the held-slot single-occupancy reasoning
                // trivial and gives the driver a clean re-dispatch point.
                if emitted {
                    return Ok(StepOutcome::Progress);
                }
                // Nothing emitted ⇒ the lowest serial is incomplete (or pending
                // is empty). Fall through to the desync check / ingress gate.
            }
        }

        // Reaching here: `pending` is empty OR its lowest serial is incomplete,
        // and the held slot is empty.

        // ── (c) Hard desync abort — predicated on a genuinely-stuck serial ───
        // Abort ONLY when the lowest serial is still incomplete after a full
        // drain AND the buffer is over the hard limit. A complete lowest serial
        // is drained in (b) before this check, so a non-empty `pending` here
        // means the lowest serial is incomplete; if bytes are also over 2x, one
        // stream is producing unboundedly more than its counterpart — a real
        // desync. (The old predicate fired on bytes alone, false-aborting an
        // in-sync but ahead stream whose lowest serial was emittable.)
        let lowest_incomplete = !self.pending.is_empty();
        if lowest_incomplete && self.pending_total_bytes > 2 * self.pending_backpressure_bytes {
            let (&serial, _) = self.pending.iter().next().expect("pending non-empty");
            return Err(io::Error::other(format!(
                "ZipFastqRecords: catastrophic stream desync — lowest chunk_serial {serial} \
                 is still incomplete while pending buffer ({} bytes) exceeds 2x backpressure \
                 limit ({} bytes). One FASTQ stream is producing far more data than its \
                 counterpart.",
                self.pending_total_bytes, self.pending_backpressure_bytes
            )));
        }

        // ── (c') Advisory soft limit — warn, never refuse ────────────────────
        // Being over the soft limit does *not* throttle ingress, and cannot: the
        // chunk that completes the lowest serial (and so lets egress drain) is
        // itself still queued in `ctx.input`, so refusing input here would stall
        // the very drain that relieves the pressure. That is the S5a1-001
        // deadlock. The hard 2x abort above is therefore the sole bound on
        // `pending`, and this arm only reports.
        //
        // An earlier revision also refused ingress when `pending` was empty, as a
        // precautionary throttle. That branch could never run: `pending_total_bytes`
        // moves only with `pending` — bytes are added on insert and the same sum
        // is subtracted when the entry is removed — so an empty `pending` always
        // means zero bytes, which is never over a positive limit. The assert below
        // pins that invariant, so a future accounting change that breaks it is
        // loud in debug builds rather than silently resurrecting dead code.
        debug_assert!(
            !self.pending.is_empty() || self.pending_total_bytes == 0,
            "pending accounting desync: empty pending must imply zero bytes, got {}",
            self.pending_total_bytes
        );
        if self.pending_total_bytes > self.pending_backpressure_bytes {
            // Once per over-limit episode, not once per dispatch.
            if self.enter_over_soft_limit() {
                log::warn!(
                    "ZipFastqRecords: pending buffer exceeds backpressure limit ({} > {}); lowest \
                     serial incomplete, continuing ingress to complete it",
                    self.pending_total_bytes,
                    self.pending_backpressure_bytes
                );
            }
        } else {
            self.leave_over_soft_limit();
        }

        // ── (d) Ingress: pop one input chunk + insert, then drain again ──────
        let Some(batch) = ctx.input.pop() else {
            // No input this call. Egress already ran in (b) and emitted nothing,
            // so the lowest serial is incomplete. If upstream is drained, this
            // is the residual-completion / EOF path: every remaining entry must
            // be complete (all streams ended cleanly) — a `None` slot means a
            // stream ended early (out of sync).
            if ctx.input.is_drained() {
                if let Some((&serial, slots)) = self.pending.iter().next() {
                    for (idx, slot) in slots.iter().enumerate() {
                        if slot.is_none() {
                            return Err(io::Error::other(format!(
                                "FASTQ sources out of sync: stream {idx} ended before \
                                 chunk_serial {serial} while other streams had more records",
                            )));
                        }
                    }
                    // All slots present but the full drain in (b) emitted
                    // nothing: unreachable (a complete lowest entry always
                    // emits). Fail loud — returning NoProgress here would hang
                    // (input is drained, so the step is never re-fed).
                    return Err(io::Error::other(format!(
                        "ZipFastqRecords: chunk_serial {serial} has all streams present but \
                         did not emit — internal invariant violated",
                    )));
                }
                return Ok(StepOutcome::Finished);
            }
            return Ok(StepOutcome::NoProgress);
        };

        let chunk_serial = batch.chunk_serial;
        let stream_idx = batch.stream_idx;
        // The batch's own `HeapSize`, matching the removal side and the measure
        // the byte-bounded queues budget — see the note in `try_emit_complete`.
        let batch_bytes: usize = batch.heap_size();

        let n = self.n_streams;
        let slots =
            self.pending.entry(chunk_serial).or_insert_with(|| (0..n).map(|_| None).collect());
        // `stream_idx` is set by the upstream splitter and is an invariant, not
        // input — but indexing on a violated invariant would panic inside a
        // worker. Fail closed with the step's own error path instead, keeping a
        // `debug_assert!` so a violation is loud in development.
        debug_assert!(stream_idx < n, "stream_idx {stream_idx} out of range for {n} streams");
        let Some(slot) = slots.get_mut(stream_idx) else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("ZipFastqRecords: stream_idx {stream_idx} out of range for {n} streams"),
            ));
        };
        *slot = Some(batch);
        self.pending_total_bytes += batch_bytes;

        // The new chunk may have completed the lowest serial. Run the full drain
        // again so a completion is emitted immediately (one unified egress
        // path). Held single-occupancy holds: the preamble guaranteed the slot
        // empty and this drain only `put`s after a `take`-confirmed empty slot.
        // Both drain outcomes are forward progress here — this dispatch consumed
        // an input chunk, and `Held` additionally parked a freshly dequeued
        // serial for retry — so report `Progress` either way (never `Contention`).
        self.drain_complete(ctx)?;
        Ok(StepOutcome::Progress)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fastq_parse::FastqRecord;
    use crate::pipeline::core::item::Ordered;

    fn make_record(name: &str, seq: &str) -> FastqRecord {
        let qual: String = std::iter::repeat_n('I', seq.len()).collect();
        let data = format!("@{name}\n{seq}\n+\n{qual}\n");
        FastqRecord::from_slice(data.as_bytes()).unwrap()
    }

    fn make_chunk(
        stream_idx: usize,
        chunk_serial: u64,
        records: Vec<FastqRecord>,
    ) -> FastqChunkBatch {
        let total_bytes: usize = records.iter().map(HeapSize::heap_size).sum();
        // ZipFastqRecords joins by (stream_idx, chunk_serial) and does not read
        // `ordinal`, so the per-chunk ordinal is irrelevant to the join. Derive
        // a deterministic unique value for clarity.
        let ordinal = chunk_serial * 8 + stream_idx as u64;
        FastqChunkBatch::new(ordinal, stream_idx, chunk_serial, records, total_bytes)
    }

    #[test]
    fn fastq_template_batch_heap_size_and_ordinal() {
        let batch = FastqTemplateBatch::new(42, vec![]);
        assert_eq!(batch.heap_size(), 0);
        assert_eq!(batch.ordinal(), 42);
    }

    #[test]
    fn test_zip_1_stream() {
        let mut zip = ZipFastqRecords::new(1, 64 * 1024 * 1024);

        let r1 = make_record("read1", "ACGT");
        let r2 = make_record("read2", "TGCA");
        let chunk = make_chunk(0, 0, vec![r1, r2]);

        // Insert the chunk.
        zip.pending.insert(0, vec![Some(chunk)]);

        let batch = zip.try_emit_complete().unwrap().unwrap();
        assert_eq!(batch.templates.len(), 2);
        assert_eq!(batch.templates[0].name, b"read1");
        assert_eq!(batch.templates[0].records.len(), 1);
        assert_eq!(batch.templates[1].name, b"read2");
        assert_eq!(batch.templates[1].records.len(), 1);
        assert_eq!(batch.ordinal(), 0);
    }

    #[test]
    fn test_zip_2_streams() {
        let mut zip = ZipFastqRecords::new(2, 64 * 1024 * 1024);

        let r1_s0 = make_record("read1/1", "ACGT");
        let r2_s0 = make_record("read2/1", "TGCA");
        let chunk_s0 = make_chunk(0, 0, vec![r1_s0, r2_s0]);

        let r1_s1 = make_record("read1/2", "GGGG");
        let r2_s1 = make_record("read2/2", "CCCC");
        let chunk_s1 = make_chunk(1, 0, vec![r1_s1, r2_s1]);

        zip.pending.insert(0, vec![Some(chunk_s0), Some(chunk_s1)]);

        let batch = zip.try_emit_complete().unwrap().unwrap();
        assert_eq!(batch.templates.len(), 2);
        assert_eq!(batch.templates[0].name, b"read1");
        assert_eq!(batch.templates[0].records.len(), 2);
        assert_eq!(batch.templates[0].records[0].sequence(), b"ACGT");
        assert_eq!(batch.templates[0].records[1].sequence(), b"GGGG");
        assert_eq!(batch.templates[1].name, b"read2");
        assert_eq!(batch.templates[1].records.len(), 2);
    }

    #[test]
    fn test_zip_3_streams() {
        let mut zip = ZipFastqRecords::new(3, 64 * 1024 * 1024);

        let mk = |suffix: &str, seq: &str| make_record(&format!("read1{suffix}"), seq);
        let chunk_s0 = make_chunk(0, 0, vec![mk("/1", "AAAA")]);
        let chunk_s1 = make_chunk(1, 0, vec![mk("/2", "CCCC")]);
        let chunk_s2 = make_chunk(2, 0, vec![mk("", "GGGG")]);

        zip.pending.insert(0, vec![Some(chunk_s0), Some(chunk_s1), Some(chunk_s2)]);

        let batch = zip.try_emit_complete().unwrap().unwrap();
        assert_eq!(batch.templates.len(), 1);
        assert_eq!(batch.templates[0].name, b"read1");
        assert_eq!(batch.templates[0].records.len(), 3);
    }

    #[test]
    fn test_zip_4_streams() {
        let mut zip = ZipFastqRecords::new(4, 64 * 1024 * 1024);

        let mk = |suffix: &str, seq: &str| make_record(&format!("readA{suffix}"), seq);
        let chunk_s0 = make_chunk(0, 0, vec![mk("/1", "AAAA")]);
        let chunk_s1 = make_chunk(1, 0, vec![mk("/2", "CCCC")]);
        let chunk_s2 = make_chunk(2, 0, vec![mk("", "GGGG")]);
        let chunk_s3 = make_chunk(3, 0, vec![mk("", "TTTT")]);

        zip.pending.insert(0, vec![Some(chunk_s0), Some(chunk_s1), Some(chunk_s2), Some(chunk_s3)]);

        let batch = zip.try_emit_complete().unwrap().unwrap();
        assert_eq!(batch.templates.len(), 1);
        assert_eq!(batch.templates[0].name, b"readA");
        assert_eq!(batch.templates[0].records.len(), 4);
    }

    #[test]
    fn test_zip_desync_bails() {
        let mut zip = ZipFastqRecords::new(2, 64 * 1024 * 1024);

        // Stream 0 has chunk_serial 0, stream 1 is missing.
        let r1 = make_record("read1/1", "ACGT");
        let chunk = make_chunk(0, 0, vec![r1]);
        zip.pending.insert(0, vec![Some(chunk), None]);

        // The drained-completion path in `try_run` walks `pending` and errors
        // on a half-filled serial. Since we can't build a `StepCtx` here, test
        // the pending state detection directly.
        let slots = zip.pending.get(&0).unwrap();
        let all_filled = slots.iter().all(Option::is_some);
        assert!(!all_filled, "incomplete slot should be detected");

        // Verify that try_emit_complete returns None for partial entries.
        let result = zip.try_emit_complete().unwrap();
        assert!(result.is_none());
    }

    #[test]
    fn test_zip_name_mismatch_errors() {
        let mut zip = ZipFastqRecords::new(2, 64 * 1024 * 1024);

        let r1_s0 = make_record("readA/1", "ACGT");
        let r1_s1 = make_record("readB/2", "GGGG");
        let chunk_s0 = make_chunk(0, 0, vec![r1_s0]);
        let chunk_s1 = make_chunk(1, 0, vec![r1_s1]);

        zip.pending.insert(0, vec![Some(chunk_s0), Some(chunk_s1)]);

        let err = zip.try_emit_complete().unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("mismatch"), "expected name mismatch error, got: {msg}");
    }

    /// Within a complete chunk, unequal per-stream record counts are rejected
    /// by the explicit strict count check (fgbio-consistent), not surfaced as a
    /// name mismatch. Matches the `N==2` `ParseAndZipFastq` path.
    #[test]
    fn test_zip_unequal_counts_out_of_sync() {
        let mut zip = ZipFastqRecords::new(2, 64 * 1024 * 1024);

        // Same chunk_serial, both streams present, but stream 0 has 2 records
        // and stream 1 has 1.
        let s0 =
            make_chunk(0, 0, vec![make_record("read1/1", "ACGT"), make_record("read2/1", "TT")]);
        let s1 = make_chunk(1, 0, vec![make_record("read1/2", "GGGG")]);
        zip.pending.insert(0, vec![Some(s0), Some(s1)]);

        let err = zip.try_emit_complete().unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("out of sync"), "expected count out-of-sync error, got: {msg}");
    }

    #[test]
    fn test_zip_multiple_serials_in_order() {
        let mut zip = ZipFastqRecords::new(1, 64 * 1024 * 1024);

        let chunk0 = make_chunk(0, 0, vec![make_record("read1", "AAAA")]);
        let chunk1 = make_chunk(0, 1, vec![make_record("read2", "CCCC")]);

        zip.pending.insert(0, vec![Some(chunk0)]);
        zip.pending.insert(1, vec![Some(chunk1)]);

        let batch0 = zip.try_emit_complete().unwrap().unwrap();
        assert_eq!(batch0.ordinal(), 0);
        assert_eq!(batch0.templates[0].name, b"read1");

        let batch1 = zip.try_emit_complete().unwrap().unwrap();
        assert_eq!(batch1.ordinal(), 1);
        assert_eq!(batch1.templates[0].name, b"read2");
    }

    #[test]
    fn profile_advertises_serial_byordinal() {
        let step = ZipFastqRecords::new(2, 1024);
        let p = step.profile();
        assert_eq!(p.name, "ZipFastqRecords");
        assert_eq!(p.kind, StepKind::Serial);
        assert!(!p.sticky);
        assert_eq!(p.branch_ordering, vec![BranchOrdering::ByItemOrdinal]);
    }

    // ─────────────────────────────────────────────────────────────────────────
    // S5a1-001 regression: lag-then-burst over the soft backpressure limit.
    //
    // The buggy code returned NoProgress/Err on backpressure *before* draining
    // the (already-complete) lowest serial, so once `pending` crossed the soft
    // limit with the lowest serial incomplete and input still queued, the step
    // wedged forever (NoProgress on every re-dispatch, never popping the chunk
    // that would complete the serial). The test below drives the REAL `try_run`
    // through a `Pipeline` with a tiny soft limit and a lag-then-burst stream
    // pattern; under the bug it deadlocks and the watchdog fires, under the fix
    // it completes and emits every zipped template in order.
    // ─────────────────────────────────────────────────────────────────────────

    use std::sync::Arc;
    use std::sync::Mutex;
    use std::sync::atomic::{AtomicUsize, Ordering as AtomicOrd};

    use crate::pipeline::core::step::Affinity;
    use crate::pipeline::core::{PipelineBuilder, PipelineConfig};

    /// `(batch_serial, template names)` recorded by the collecting sink.
    type SeenBatches = Arc<Mutex<Vec<(u64, Vec<Vec<u8>>)>>>;

    /// Bytes per FASTQ record this source emits. Sized so a small run of
    /// records crosses the few-KiB soft limit the test sets on the zipper.
    const TEST_RECORD_SEQ_LEN: usize = 200;
    /// Number of chunk serials emitted per stream.
    const TEST_N_SERIALS: u64 = 12;
    /// Soft backpressure limit for the zipper under test. Sized so the lagging
    /// streams (streams 0 and 1, `TEST_N_SERIALS` records each at
    /// `TEST_RECORD_SEQ_LEN` bytes) cross the soft limit before the stream-2
    /// burst arrives, yet stay UNDER the 2x hard-desync limit — the lag is a
    /// legitimately-ahead in-sync stream, not a genuine desync, so the hard arm
    /// must NOT fire. (~12 records/stream * 2 streams * ~280 B ≈ 6.7 KiB, which
    /// is > soft (5 KiB) and < hard (10 KiB).)
    const TEST_SOFT_LIMIT_BYTES: usize = 5 * 1024;

    fn big_record(name: &str) -> FastqRecord {
        let seq: String = std::iter::repeat_n('A', TEST_RECORD_SEQ_LEN).collect();
        make_record(name, &seq)
    }

    /// Serial source that emits `FastqChunkBatch` items in a deliberately
    /// skewed order: all of streams 0 and 1 (serials `0..N`) FIRST, then all of
    /// stream 2 in a burst. The leading streams pile up incomplete lowest
    /// serials in the zipper (stream 2's slot is missing) and push it over the
    /// soft limit before the completing burst arrives.
    struct LagThenBurstSource {
        /// Emission plan: each entry is `(stream_idx, chunk_serial)`.
        plan: Vec<(usize, u64)>,
        next: usize,
        ordinals: crate::pipeline::steps::source::read_fastq::FastqOrdinalSequence,
        held: HeldSlot<Unpushed<FastqChunkBatch>>,
        output_byte_limit: u64,
    }

    impl LagThenBurstSource {
        fn new(n_serials: u64, n_streams_total: usize, output_byte_limit: u64) -> Self {
            let mut plan = Vec::new();
            // Streams 0..n-1 first (the "lag"): every serial of every leading
            // stream before any of the final stream.
            for stream_idx in 0..n_streams_total - 1 {
                for serial in 0..n_serials {
                    plan.push((stream_idx, serial));
                }
            }
            // Final stream in a burst: a run of consecutive lowest serials
            // becomes complete at once.
            for serial in 0..n_serials {
                plan.push((n_streams_total - 1, serial));
            }
            Self {
                plan,
                next: 0,
                ordinals: crate::pipeline::steps::source::read_fastq::FastqOrdinalSequence::new(),
                held: HeldSlot::new(),
                output_byte_limit,
            }
        }

        /// A plan in which only stream 0 ever emits, so the lowest serial can
        /// never complete and `pending` grows without bound — a genuine desync,
        /// as opposed to the in-sync-but-ahead lag the burst plan models.
        fn with_only_the_first_stream(n_serials: u64, output_byte_limit: u64) -> Self {
            let plan = (0..n_serials).map(|serial| (0usize, serial)).collect();
            Self {
                plan,
                next: 0,
                ordinals: crate::pipeline::steps::source::read_fastq::FastqOrdinalSequence::new(),
                held: HeldSlot::new(),
                output_byte_limit,
            }
        }

        /// The same plan with the final stream stopped one serial short, so a
        /// `pending` entry still has an empty slot when input drains — the
        /// truncated-input shape the EOF-desync guard exists to diagnose.
        fn with_final_stream_truncated(
            n_serials: u64,
            n_streams_total: usize,
            output_byte_limit: u64,
        ) -> Self {
            let mut source = Self::new(n_serials, n_streams_total, output_byte_limit);
            source.plan.pop().expect("plan is non-empty");
            source
        }
    }

    impl Step for LagThenBurstSource {
        type Input = ();
        type Outputs = OrderedBytesSingle<FastqChunkBatch>;

        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "LagThenBurstSource",
                kind: StepKind::Serial,
                sticky: false,
                output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
                branch_ordering: vec![BranchOrdering::ByItemOrdinal],
            }
        }

        fn affinity(&self) -> Affinity {
            Affinity::Reader
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

            let Some(&(stream_idx, chunk_serial)) = self.plan.get(self.next) else {
                return Ok(StepOutcome::Finished);
            };
            self.next += 1;

            // Trailing non-digit so `strip_read_suffix` does not eat the serial
            // (it strips a `_`/`/`/`.`/`:` + `1`/`2` pair-suffix).
            let name = format!("read{chunk_serial}x");
            let record = big_record(&name);
            let total_bytes = record.heap_size();
            // Draw from the same shared sequence production uses, so the
            // harness cannot drift from what the reader emits. Ordinals follow
            // emission order; the zipper keys on `chunk_serial`, so the skewed
            // plan this source exists to produce is unaffected.
            let ordinal = self.ordinals.next();
            let chunk =
                FastqChunkBatch::new(ordinal, stream_idx, chunk_serial, vec![record], total_bytes);

            match ctx.outputs.push(chunk) {
                Ok(()) => Ok(StepOutcome::Progress),
                Err(unpushed) => {
                    self.held.put(unpushed);
                    Ok(StepOutcome::Progress)
                }
            }
        }
    }

    /// Serial (single-consumer) sink that records, per emitted batch, its
    /// `batch_serial` and the names of the templates it carries, so the test can
    /// assert full drain and correct zipping/order. A `Serial` sink drains the
    /// queue in FIFO order on one worker, so the recorded sequence is the true
    /// emission order — letting the test assert ordering directly instead of
    /// sorting (which would only prove set equality).
    #[derive(Clone)]
    struct CollectingSink {
        seen: SeenBatches,
        count: Arc<AtomicUsize>,
    }

    impl Step for CollectingSink {
        type Input = FastqTemplateBatch;
        type Outputs = ();

        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "CollectingSink",
                kind: StepKind::Serial,
                sticky: false,
                output_queues: vec![],
                branch_ordering: vec![],
            }
        }

        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            match ctx.input.pop() {
                Some(batch) => {
                    let names: Vec<Vec<u8>> =
                        batch.templates.iter().map(|t| t.name.clone()).collect();
                    self.seen.lock().unwrap().push((batch.batch_serial, names));
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

    /// The soft-limit warning must fire once per over-limit *episode*, not once
    /// per dispatch. The step is `Serial` and the driver re-dispatches it in a
    /// tight loop, so an unlatched warning floods the log for the whole time an
    /// ahead-but-in-sync stream sits over the limit — burying everything else.
    #[test]
    fn soft_limit_warning_latches_once_per_episode() {
        let mut zip = ZipFastqRecords::new(1, 64 * 1024);

        assert!(zip.enter_over_soft_limit(), "the dispatch that enters the episode must warn");
        assert!(!zip.enter_over_soft_limit(), "a later dispatch in the same episode must not");
        assert!(!zip.enter_over_soft_limit(), "nor any after that");

        zip.leave_over_soft_limit();

        assert!(zip.enter_over_soft_limit(), "a fresh episode must warn again");
    }

    /// Sized so a stream-0-only plan clears the 2x hard-desync ceiling with room
    /// to spare: 12 serials * ~280 B/record is ~3.4 KB, well past 2 * 512 B.
    /// Deliberately not `TEST_SOFT_LIMIT_BYTES` — that one is tuned so the
    /// lag-then-burst case stays *under* the hard ceiling.
    const TEST_HARD_ABORT_LIMIT_BYTES: usize = 512;

    /// A genuinely desynced stream must abort, and this is the only bound on
    /// `pending`.
    ///
    /// The soft limit deliberately does not refuse ingress — doing so deadlocks,
    /// since the chunk that completes the lowest serial is itself still queued
    /// upstream — so this 2x abort is what stops an unbounded buffer. Its
    /// predicate is a conjunction (lowest serial still incomplete after a full
    /// drain AND bytes past 2x), and a regression weakening either half, or
    /// reverting to a bytes-only predicate, would leave the suite green while
    /// removing the ceiling entirely.
    ///
    /// This pins the direction `zip_does_not_deadlock_on_lag_then_burst_over_soft_limit`
    /// cannot: that test proves the abort does *not* fire for an ahead-but-in-sync
    /// stream. Together they fix both sides of the predicate.
    #[test]
    fn zip_aborts_on_a_genuine_stream_desync() {
        let n_streams = 2usize;
        let seen = Arc::new(Mutex::new(Vec::<(u64, Vec<Vec<u8>>)>::new()));
        let count = Arc::new(AtomicUsize::new(0));

        // Only stream 0 emits, so serial 0 never completes and every chunk stays
        // resident in `pending`.
        let source = LagThenBurstSource::with_only_the_first_stream(TEST_N_SERIALS, 64 * 1024);
        let zip = ZipFastqRecords::new(n_streams, 64 * 1024)
            .with_backpressure_bytes(TEST_HARD_ABORT_LIMIT_BYTES);
        let sink = CollectingSink { seen: Arc::clone(&seen), count: Arc::clone(&count) };

        let builder = PipelineBuilder::new();
        builder.chain(source).chain(zip).chain(sink).into_sink_marker();
        let pipeline = builder.build().expect("chain builds");

        let err = pipeline
            .run(PipelineConfig { threads: 2, ..Default::default() })
            .expect_err("a one-sided stream must abort rather than buffer without bound");

        let msg = err.to_string();
        assert!(
            msg.contains("catastrophic stream desync"),
            "expected the hard 2x desync abort, got: {msg}"
        );
    }

    /// A stream that ends before its counterparts must fail the run with a
    /// diagnosable message rather than silently dropping the unmatched records.
    ///
    /// This is the guard that turns a truncated FASTQ input — one file cut short,
    /// a partial download — into an error naming the stream and serial. It sits
    /// on the input-drained path inside `try_run`, so it needs a real `Pipeline`
    /// to reach: the `try_emit_complete` unit tests above cannot drive it.
    #[test]
    fn zip_reports_a_stream_that_ends_early() {
        let n_streams = 2usize;
        let seen = Arc::new(Mutex::new(Vec::<(u64, Vec<Vec<u8>>)>::new()));
        let count = Arc::new(AtomicUsize::new(0));

        let source =
            LagThenBurstSource::with_final_stream_truncated(TEST_N_SERIALS, n_streams, 64 * 1024);
        let zip = ZipFastqRecords::new(n_streams, 64 * 1024);
        let sink = CollectingSink { seen: Arc::clone(&seen), count: Arc::clone(&count) };

        let builder = PipelineBuilder::new();
        builder.chain(source).chain(zip).chain(sink).into_sink_marker();
        let pipeline = builder.build().expect("chain builds");

        let err = pipeline
            .run(PipelineConfig { threads: 2, ..Default::default() })
            .expect_err("a stream that ends early must fail the run");

        let msg = err.to_string();
        assert!(
            msg.contains("ended before"),
            "expected the EOF-desync diagnosis naming the stream and serial, got: {msg}"
        );
    }

    /// Drive `ZipFastqRecords` through a real `Pipeline` under a lag-then-burst
    /// stream pattern that crosses its (shrunk) soft backpressure limit, on a
    /// watchdog thread so a deadlock regression surfaces as a timeout rather
    /// than hanging the test runner.
    ///
    /// Pre-fix this deadlocks: the zipper crosses the soft limit while the
    /// lowest serial is incomplete, returns `NoProgress` before draining/popping,
    /// and never admits the completing burst → the watchdog fires. Post-fix the
    /// advisory soft limit admits the completing chunks, egress fully drains,
    /// and the run completes with every template zipped in order.
    #[test]
    fn zip_does_not_deadlock_on_lag_then_burst_over_soft_limit() {
        use std::time::{Duration, Instant};

        let n_streams = 3usize;
        let seen = Arc::new(Mutex::new(Vec::<(u64, Vec<Vec<u8>>)>::new()));
        let count = Arc::new(AtomicUsize::new(0));

        let seen_for_thread = Arc::clone(&seen);
        let count_for_thread = Arc::clone(&count);
        let worker = std::thread::Builder::new()
            .name("zip-lag-burst".into())
            .spawn(move || -> io::Result<()> {
                let source = LagThenBurstSource::new(TEST_N_SERIALS, n_streams, 64 * 1024);
                let zip = ZipFastqRecords::new(n_streams, 64 * 1024)
                    .with_backpressure_bytes(TEST_SOFT_LIMIT_BYTES);
                let sink = CollectingSink { seen: seen_for_thread, count: count_for_thread };

                let builder = PipelineBuilder::new();
                builder.chain(source).chain(zip).chain(sink).into_sink_marker();
                let pipeline = builder.build().map_err(io::Error::other)?;
                // A modest thread count keeps the Serial zipper single-driven
                // while letting the source/sink run on other workers.
                pipeline
                    .run(PipelineConfig { threads: 2, ..Default::default() })
                    .map_err(io::Error::other)
            })
            .expect("spawn pipeline worker");

        // Generous timeout: the run completes in well under a second on a
        // healthy system; 30 s leaves headroom for slow CI while still failing
        // fast on the deadlock regression.
        let deadline = Instant::now() + Duration::from_secs(30);
        while !worker.is_finished() {
            assert!(
                Instant::now() < deadline,
                "ZipFastqRecords pipeline did not finish within 30s — likely the S5a1-001 \
                 backpressure-starves-egress deadlock regression"
            );
            std::thread::sleep(Duration::from_millis(25));
        }
        worker.join().expect("pipeline thread panicked").expect("pipeline run returned Err");

        // Full drain: one batch per serial, exactly TEST_N_SERIALS batches.
        let seen = seen.lock().unwrap().clone();
        assert_eq!(
            seen.len(),
            usize::try_from(TEST_N_SERIALS).unwrap(),
            "expected one emitted batch per serial (full drain), got {}",
            seen.len()
        );

        // Correctness + order: the single-consumer sink records the true emission
        // order, so assert `batch_serial`s are contiguous 0..N directly (no sort)
        // — proving `ZipFastqRecords` emits in order, not merely as a set. Each
        // batch carries the single expected zipped template `read_<serial>`.
        for (i, (serial, names)) in seen.iter().enumerate() {
            assert_eq!(*serial, i as u64, "batch serials must be contiguous and in order");
            assert_eq!(names.len(), 1, "each chunk carried exactly one record");
            assert_eq!(
                names[0],
                format!("read{i}x").into_bytes(),
                "template name must match the source record for serial {i}"
            );
        }
    }
}

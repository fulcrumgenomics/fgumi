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

pub struct ZipFastqRecords {
    n_streams: usize,
    pending: BTreeMap<u64, Vec<Option<FastqChunkBatch>>>,
    pending_total_bytes: usize,
    pending_backpressure_bytes: usize,
    held: HeldSlot<Unpushed<FastqTemplateBatch>>,
    output_byte_limit: u64,
    next_batch_serial: u64,
}

impl ZipFastqRecords {
    #[must_use]
    pub fn new(n_streams: usize, output_byte_limit: u64) -> Self {
        Self {
            n_streams,
            pending: BTreeMap::new(),
            pending_total_bytes: 0,
            pending_backpressure_bytes: DEFAULT_PENDING_BACKPRESSURE_BYTES,
            held: HeldSlot::new(),
            output_byte_limit,
            next_batch_serial: 0,
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

        let chunks = self.pending.remove(&lowest_serial).unwrap();
        let mut stream_records: Vec<Vec<FastqRecord>> =
            chunks.into_iter().map(|opt| opt.unwrap().records).collect();

        // Subtract the bytes from pending accounting.
        let chunk_bytes: usize =
            stream_records.iter().flat_map(|recs| recs.iter()).map(HeapSize::heap_size).sum();
        self.pending_total_bytes = self.pending_total_bytes.saturating_sub(chunk_bytes);

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

        let mut templates = Vec::with_capacity(n_records);

        // Zip across streams: record i from each stream forms one template.
        // Drain from the back (via pop) to avoid O(n^2) from front removal;
        // reverse after.
        let mut reversed_templates = Vec::with_capacity(n_records);
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

            reversed_templates.push(FastqTemplate { name: base_name.unwrap(), records });
        }

        // We popped from the back, so reverse to restore original order.
        reversed_templates.reverse();
        templates.extend(reversed_templates);

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
        // 1. Drain held slot.
        if let Some(unpushed) = self.held.take() {
            match ctx.outputs.retry(unpushed) {
                Ok(()) => {}
                Err(again) => {
                    self.held.put(again);
                    return Ok(StepOutcome::Contention);
                }
            }
        }

        // 2. Check pending backpressure (D3).
        if self.pending_total_bytes > 2 * self.pending_backpressure_bytes {
            return Err(io::Error::other(format!(
                "ZipFastqRecords: catastrophic stream desync — pending buffer ({} bytes) \
                 exceeds 2x backpressure limit ({} bytes). One FASTQ stream is producing \
                 far more data than its counterpart.",
                self.pending_total_bytes, self.pending_backpressure_bytes
            )));
        }
        if self.pending_total_bytes > self.pending_backpressure_bytes {
            log::warn!(
                "ZipFastqRecords: pending buffer exceeds backpressure limit ({} > {})",
                self.pending_total_bytes,
                self.pending_backpressure_bytes
            );
            return Ok(StepOutcome::NoProgress);
        }

        // 3. Pop a chunk from input.
        let Some(batch) = ctx.input.pop() else {
            // No input this call. Emit any pending complete entry first.
            if let Some(template_batch) = self.try_emit_complete()? {
                if let Err(unpushed) = ctx.outputs.push(template_batch) {
                    self.held.put(unpushed);
                }
                return Ok(StepOutcome::Progress);
            }
            // Nothing emittable. If upstream is drained, finalize. `held` is
            // empty here (step 1 returned `Contention` otherwise), so the
            // remaining work is the residual `pending` walk: every remaining
            // entry must be complete (all streams ended cleanly) — a `None`
            // slot means a stream ended early (out of sync). `Finished` once
            // `pending` is fully drained.
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
                    // All slots present but `try_emit_complete` did not emit:
                    // unreachable (a complete lowest entry always emits). Fail
                    // loud — returning `NoProgress` here would hang (input is
                    // drained, so the step would never be re-fed and never
                    // reach `Finished`).
                    return Err(io::Error::other(format!(
                        "ZipFastqRecords: chunk_serial {serial} has all streams present but \
                         did not emit — internal invariant violated",
                    )));
                }
                return Ok(StepOutcome::Finished);
            }
            return Ok(StepOutcome::NoProgress);
        };

        // 4. Insert into pending.
        let chunk_serial = batch.chunk_serial;
        let stream_idx = batch.stream_idx;
        let batch_bytes: usize = batch.records.iter().map(HeapSize::heap_size).sum();

        let n = self.n_streams;
        let slots =
            self.pending.entry(chunk_serial).or_insert_with(|| (0..n).map(|_| None).collect());
        slots[stream_idx] = Some(batch);
        self.pending_total_bytes += batch_bytes;

        // 5. Try to emit if the lowest serial is complete.
        if let Some(template_batch) = self.try_emit_complete()? {
            if let Err(unpushed) = ctx.outputs.push(template_batch) {
                self.held.put(unpushed);
            }
        }

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
}

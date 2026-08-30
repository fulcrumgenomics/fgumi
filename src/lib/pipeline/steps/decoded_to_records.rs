//! `DecodedRecordBatchToRecordBatch` adapter step.
//!
//! Bridges the decoded-record view (`DecodedRecordBatch`, emitted by the SAM
//! parse preamble's `ParseSamChunk`) to the flat-record view (`RecordBatch`,
//! consumed by the sort ingest `SortBuffer`). `SortBuffer` operates on raw
//! record bytes and never reads a `DecodedRecord`'s pre-computed `GroupKey`
//! (sort computes its own sort keys straight off the raw bytes), so this step
//! discards the key and repacks each record's raw bytes — mirroring
//! `TemplatesToRecordBatch`'s flatten of `BamTemplateBatch`, the other
//! non-arena source shape `add_sort` accepts.
//!
//! Exists specifically so a standalone (or sort-first-intermediate) chain over
//! a SAM source can reach the sort stage: SAM has no BGZF framing for the
//! parallel-inflate arena front (`ReadBlocks → InflateToArena →
//! FindBoundariesAndSort`), so it always decodes through `ParseSamChunk`
//! first, same as every other stage's SAM ingest.
//!
//! `Parallel + ByItemOrdinal` — the conversion is per-batch stateless, so the
//! framework can clone this step across workers.

use std::io;

use crate::pipeline::core::Unpushed;
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::outputs::OrderedBytesSingle;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};
use crate::pipeline::steps::types::{DecodedRecordBatch, RecordBatch, RecordBatchBuilder};

/// `Parallel + ByItemOrdinal` adapter step.
///
/// Each `try_run` pops one [`DecodedRecordBatch`], appends every record's raw
/// bytes into a [`RecordBatchBuilder`], and emits the resulting [`RecordBatch`]
/// downstream. Backpressure is handled via a [`HeldSlot`] preserving the input
/// batch's serial ordinal.
pub struct DecodedRecordBatchToRecordBatch {
    held: HeldSlot<Unpushed<RecordBatch>>,
    output_byte_limit: u64,
}

impl DecodedRecordBatchToRecordBatch {
    #[must_use]
    pub fn new(output_byte_limit: u64) -> Self {
        Self { held: HeldSlot::new(), output_byte_limit }
    }
}

impl Clone for DecodedRecordBatchToRecordBatch {
    fn clone(&self) -> Self {
        Self { held: HeldSlot::new(), output_byte_limit: self.output_byte_limit }
    }
}

impl Step for DecodedRecordBatchToRecordBatch {
    type Input = DecodedRecordBatch;
    type Outputs = OrderedBytesSingle<RecordBatch>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "DecodedRecordBatchToRecordBatch",
            kind: StepKind::Parallel,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if let Some(unpushed) = self.held.take() {
            match ctx.outputs.retry(unpushed) {
                Ok(()) => {}
                Err(again) => {
                    self.held.put(again);
                    // `Contention` (not `NoProgress`) — see `TemplatesToRecordBatch`
                    // (and `BgzfDecompress`) for the rationale: `NoProgress` from a
                    // Parallel step on a held slot risks the framework marking this
                    // worker `Skip` if the upstream input is also drained, silently
                    // dropping the held item.
                    return Ok(StepOutcome::Contention);
                }
            }
        }

        let Some(batch) = ctx.input.pop() else {
            if ctx.input.is_drained() {
                return Ok(StepOutcome::Finished);
            }
            return Ok(StepOutcome::NoProgress);
        };

        let serial = batch.batch_serial();
        // Tight capacity from the exact sum of record body bytes, matching
        // `TemplatesToRecordBatch`'s approach (rather than `total_bytes()`,
        // which over-estimates via each record's allocated capacity).
        let bytes_cap: usize = batch.records().iter().map(|r| r.raw_bytes().len()).sum();
        let records_cap = batch.records().len();
        let mut builder = RecordBatchBuilder::with_capacity(serial, bytes_cap, records_cap);
        for record in batch.into_records() {
            builder.push_record_bytes(record.raw_bytes());
        }
        let rb = builder.build();

        match ctx.outputs.push(rb) {
            Ok(()) => Ok(StepOutcome::Progress),
            Err(unpushed) => {
                self.held.put(unpushed);
                Ok(StepOutcome::Progress)
            }
        }
    }

    fn new_worker_copy(&self) -> Self {
        self.clone()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fgumi_bam_io::{DecodedRecord, GroupKey};
    use fgumi_raw_bam::SamBuilder;

    fn make_decoded_record(qname: &[u8]) -> DecodedRecord {
        let mut b = SamBuilder::new();
        b.read_name(qname).flags(0).sequence(b"ACGT").qualities(b"IIII");
        DecodedRecord::from_raw_bytes(b.build(), GroupKey::default())
    }

    #[test]
    fn profile_advertises_parallel_byordinal() {
        let s = DecodedRecordBatchToRecordBatch::new(1024);
        let p = s.profile();
        assert_eq!(p.name, "DecodedRecordBatchToRecordBatch");
        assert_eq!(p.kind, StepKind::Parallel);
        assert!(!p.sticky);
        assert_eq!(p.branch_ordering, vec![BranchOrdering::ByItemOrdinal]);
    }

    #[test]
    fn clone_resets_held_slot() {
        let mut s = DecodedRecordBatchToRecordBatch::new(1024);
        let _ = s.held.take(); // empty
        let cloned = s.clone();
        assert!(!cloned.held.is_held());
    }

    /// Flattening preserves record count, serial, and each record's raw bytes
    /// verbatim (sort needs byte-identical records, just reordered).
    #[test]
    fn flatten_preserves_record_count_serial_and_bytes() {
        let batch = DecodedRecordBatch::new(
            7,
            vec![make_decoded_record(b"q1"), make_decoded_record(b"q2")],
        );
        let expected_bytes: Vec<Vec<u8>> =
            batch.records().iter().map(|r| r.raw_bytes().to_vec()).collect();

        let bytes_cap: usize = batch.records().iter().map(|r| r.raw_bytes().len()).sum();
        let records_cap = batch.records().len();
        let serial = batch.batch_serial();
        let mut builder = RecordBatchBuilder::with_capacity(serial, bytes_cap, records_cap);
        for record in batch.into_records() {
            builder.push_record_bytes(record.raw_bytes());
        }
        let rb = builder.build();

        assert_eq!(rb.batch_serial(), 7);
        assert_eq!(rb.len(), 2, "2 records flattened from 2 decoded records");
        let actual_bytes: Vec<Vec<u8>> = rb.iter_record_bytes().map(<[u8]>::to_vec).collect();
        assert_eq!(actual_bytes, expected_bytes, "record bytes must survive verbatim");
    }
}

//! `TemplatesToRecordBatch` adapter step.
//!
//! Bridges the queryname-template view (`BamTemplateBatch`,
//! emitted by `AlignAndMergeStep` and `ZipperMergeStep`) to the
//! flat-record view (`RecordBatch`, consumed by the sort ingest
//! (`SortBuffer`), `DecodeRecords`, etc.). Flattens each `Template`'s records into
//! a single `RecordBatch` carrying the input batch's serial
//! ordinal, so `BranchOrdering::ByItemOrdinal` consumers downstream
//! see a monotonic sequence.
//!
//! `Parallel + ByItemOrdinal` — the conversion is per-batch
//! stateless, so the framework can clone this step across workers.
//! Used to fuse AAM-fed runall pipelines (eliminating the
//! tempfile bridge between AAM and the downstream sort).

use std::io;

use crate::pipeline::core::Unpushed;
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::outputs::OrderedBytesSingle;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};
use crate::pipeline::steps::types::{BamTemplateBatch, RecordBatch, RecordBatchBuilder};

/// `Parallel + ByItemOrdinal` adapter step.
///
/// Each `try_run` pops one [`BamTemplateBatch`], appends every
/// record from every contained template into a
/// [`RecordBatchBuilder`], and emits the resulting [`RecordBatch`]
/// downstream. Backpressure is handled via a [`HeldSlot`]
/// preserving the input batch's serial ordinal.
pub struct TemplatesToRecordBatch {
    held: HeldSlot<Unpushed<RecordBatch>>,
    output_byte_limit: u64,
}

impl TemplatesToRecordBatch {
    #[must_use]
    pub fn new(output_byte_limit: u64) -> Self {
        Self { held: HeldSlot::new(), output_byte_limit }
    }
}

impl Clone for TemplatesToRecordBatch {
    fn clone(&self) -> Self {
        Self { held: HeldSlot::new(), output_byte_limit: self.output_byte_limit }
    }
}

impl Step for TemplatesToRecordBatch {
    type Input = BamTemplateBatch;
    type Outputs = OrderedBytesSingle<RecordBatch>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "TemplatesToRecordBatch",
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
                    // `Contention` (not `NoProgress`) — see
                    // `BgzfDecompress` for the rationale. Returning
                    // `NoProgress` from a Parallel step on a held
                    // slot risks the framework marking this worker
                    // `Skip` if the upstream input is also drained,
                    // silently dropping the held item.
                    return Ok(StepOutcome::Contention);
                }
            }
        }

        let Some(batch) = ctx.input.pop() else {
            // No input this call. If upstream is drained, every item has been
            // processed (held output was flushed by the Contention preamble
            // above) and this step will never push again — report Finished.
            // For a Parallel step only the last clone to finish closes the
            // shared output (gated by the StepDrainCounter in the driver).
            if ctx.input.is_drained() {
                return Ok(StepOutcome::Finished);
            }
            return Ok(StepOutcome::NoProgress);
        };

        let serial = batch.batch_serial;
        // Pre-size the backing buffer from the exact sum of record
        // body bytes. `batch.total_bytes` over-estimates because
        // `Template::heap_size` includes the queryname bytes, but
        // those don't appear in the `RecordBatch.backing` (only the
        // body bytes do). Computing the tight number once is cheap
        // and keeps the builder's allocation in a single mimalloc
        // size class.
        let bytes_cap: usize = batch
            .templates
            .iter()
            .flat_map(|t| t.records.iter().map(fgumi_raw_bam::RawRecord::len))
            .sum();
        let records_cap: usize = batch.templates.iter().map(|t| t.records.len()).sum();
        let mut builder = RecordBatchBuilder::with_capacity(serial, bytes_cap, records_cap);
        for template in batch.templates {
            for record in template.records {
                builder.push_record_bytes(record.as_ref());
            }
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
    use crate::template::Template;
    use fgumi_raw_bam::SamBuilder;

    fn make_record(qname: &[u8], flags: u16) -> fgumi_raw_bam::RawRecord {
        let mut b = SamBuilder::new();
        b.read_name(qname).flags(flags).sequence(b"ACGT").qualities(b"IIII");
        b.build()
    }

    /// Build a paired-end template: one R1 primary + one R2 primary
    /// (the minimum that satisfies `Template::from_records`).
    fn make_paired_template(qname: &[u8]) -> Template {
        use fgumi_raw_bam::flags::{FIRST_SEGMENT, LAST_SEGMENT, PAIRED};
        let r1 = make_record(qname, PAIRED | FIRST_SEGMENT);
        let r2 = make_record(qname, PAIRED | LAST_SEGMENT);
        Template::from_records(vec![r1, r2]).expect("paired template")
    }

    #[test]
    fn profile_advertises_parallel_byordinal() {
        let s = TemplatesToRecordBatch::new(1024);
        let p = s.profile();
        assert_eq!(p.name, "TemplatesToRecordBatch");
        assert_eq!(p.kind, StepKind::Parallel);
        assert!(!p.sticky);
        assert_eq!(p.branch_ordering, vec![BranchOrdering::ByItemOrdinal]);
    }

    #[test]
    fn clone_resets_held_slot() {
        let mut s = TemplatesToRecordBatch::new(1024);
        // We can't easily put a real Unpushed into held without a
        // live framework, so just assert the clone has a fresh slot
        // (no panic in HeldSlot::new(), which would fire if we
        // accidentally cloned a held item).
        let _ = s.held.take(); // empty
        let cloned = s.clone();
        assert!(!cloned.held.is_held());
    }

    /// Verify the flattening logic: 2 paired templates → 4 records,
    /// preserving the input batch's serial ordinal.
    #[test]
    fn flatten_preserves_record_count_and_serial() {
        let batch = BamTemplateBatch::new(
            42,
            vec![make_paired_template(b"q1"), make_paired_template(b"q2")],
        );
        let total_bytes = batch.total_bytes;
        let serial = batch.batch_serial;
        let records_cap: usize = batch.templates.iter().map(|t| t.records.len()).sum();
        let mut builder = RecordBatchBuilder::with_capacity(serial, total_bytes, records_cap);
        for template in batch.templates {
            for record in template.records {
                builder.push_record_bytes(record.as_ref());
            }
        }
        let rb = builder.build();
        assert_eq!(rb.batch_serial(), 42);
        assert_eq!(rb.len(), 4, "4 records flattened from 2 paired templates");
    }
}

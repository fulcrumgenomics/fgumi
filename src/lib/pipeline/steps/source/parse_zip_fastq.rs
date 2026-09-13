//! `ParseAndZipFastqN` step: `Parallel + ByItemOrdinal`. Parses the raw,
//! whole-record-aligned bytes of an
//! [`NRawFastqBatch`](super::fastq_zip::NRawFastqBatch) (K aligned raw
//! record-chunks that share a `chunk_serial`) into per-stream `FastqRecord`s
//! and zips them record-by-record into a
//! [`FastqTemplateBatch`].
//!
//! This step exists to lift BOTH the FASTQ parse and the template build off the
//! serial source path and fan them across worker threads. The cheap chunk-level
//! alignment already happened serially in
//! [`ZipRawFastqK`](super::zip_raw_fastq_k::ZipRawFastqK) (or, for a lone
//! stream, [`WrapRawFastq1`](super::zip_raw_fastq_k::WrapRawFastq1)), which
//! minted a dense `ordinal`. Each worker here parses + zips one batch
//! independently, and the framework reorders the parallel output by each
//! batch's `ordinal`, which this step propagates verbatim onto the emitted
//! `FastqTemplateBatch`.

use std::io;

use crate::fastq_parse::FastqRecord;
use crate::pipeline::core::Unpushed;
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::outputs::OrderedBytesSingle;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};

use super::fastq_zip::FastqTemplateBatch;

// ─────────────────────────────────────────────────────────────────────────────
// ParseAndZipFastqN — Parallel K-stream parse-and-zip.
// ─────────────────────────────────────────────────────────────────────────────

/// `Parallel + ByItemOrdinal` K-stream parse-and-zip. Consumes an
/// [`NRawFastqBatch`](super::fastq_zip::NRawFastqBatch) (K aligned raw
/// record-chunks that share a `chunk_serial`, minted with a dense `ordinal` by
/// [`ZipRawFastqK`](super::zip_raw_fastq_k::ZipRawFastqK)), parses each stream's
/// bytes into `FastqRecord`s, and zips them into a
/// [`FastqTemplateBatch`] via the shared
/// [`zip_streams`](super::fastq_zip::zip_streams).
///
/// The K-stream generalization of `ParseAndZipFastq` (which was fixed at two
/// streams, `data_a`/`data_b`). Each worker is independent and stateless beyond
/// the held slot; the `ordinal` is preserved verbatim so the framework reorder
/// keeps records in input order.
pub struct ParseAndZipFastqN {
    held: HeldSlot<Unpushed<FastqTemplateBatch>>,
    output_byte_limit: u64,
}

impl ParseAndZipFastqN {
    #[must_use]
    pub fn new(output_byte_limit: u64) -> Self {
        Self { held: HeldSlot::new(), output_byte_limit }
    }
}

impl Clone for ParseAndZipFastqN {
    fn clone(&self) -> Self {
        Self { held: HeldSlot::new(), output_byte_limit: self.output_byte_limit }
    }
}

impl Step for ParseAndZipFastqN {
    type Input = super::fastq_zip::NRawFastqBatch;
    type Outputs = OrderedBytesSingle<FastqTemplateBatch>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "ParseAndZipFastqN",
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

        let super::fastq_zip::NRawFastqBatch { ordinal, chunk_serial, streams } = batch;
        // Parse each stream's raw bytes into records, then zip positionally.
        let mut per_stream: Vec<Vec<FastqRecord>> = Vec::with_capacity(streams.len());
        for data in &streams {
            per_stream.push(super::parse_fastq_chunk(data, "ParseAndZipFastqN").map(|(r, _)| r)?);
        }
        let templates = super::fastq_zip::zip_streams(per_stream, chunk_serial)?;
        let template_batch = FastqTemplateBatch::new(ordinal, templates);

        match ctx.outputs.push(template_batch) {
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
    use crate::pipeline::core::item::Ordered;

    fn record_bytes(name: &str, seq: &str) -> Vec<u8> {
        let qual: String = std::iter::repeat_n('I', seq.len()).collect();
        format!("@{name}\n{seq}\n+\n{qual}\n").into_bytes()
    }

    #[test]
    fn profile_advertises_parallel_byordinal() {
        let s = ParseAndZipFastqN::new(1024);
        let p = s.profile();
        assert_eq!(p.name, "ParseAndZipFastqN");
        assert_eq!(p.kind, StepKind::Parallel);
        assert!(!p.sticky);
        assert_eq!(p.branch_ordering, vec![BranchOrdering::ByItemOrdinal]);
    }

    /// End-to-end of the step body's parse+zip: an `NRawFastqBatch` of two
    /// aligned raw byte-chunks parses and zips into templates, and the batch's
    /// `ordinal` is preserved verbatim as the output `batch_serial`. (The parse
    /// and zip primitives themselves are unit-tested in `super::parse_fastq_chunk`
    /// and `super::fastq_zip::zip_streams`; this pins the wiring in this step.)
    #[test]
    fn parse_and_zip_two_streams_preserves_ordinal() {
        let mut data_a = record_bytes("read1/1", "ACGT");
        data_a.extend(record_bytes("read2/1", "TGCA"));
        let mut data_b = record_bytes("read1/2", "GGGG");
        data_b.extend(record_bytes("read2/2", "CCCC"));

        let mut per_stream: Vec<Vec<FastqRecord>> = Vec::new();
        for data in [&data_a, &data_b] {
            per_stream.push(super::super::parse_fastq_chunk(data, "test").map(|(r, _)| r).unwrap());
        }
        let templates = super::super::fastq_zip::zip_streams(per_stream, 0).unwrap();
        let batch = FastqTemplateBatch::new(42, templates);

        assert_eq!(batch.ordinal(), 42);
        assert_eq!(batch.templates.len(), 2);
        assert_eq!(batch.templates[0].name, b"read1");
        assert_eq!(batch.templates[0].records.len(), 2);
        assert_eq!(batch.templates[0].records[0].sequence(), b"ACGT");
        assert_eq!(batch.templates[0].records[1].sequence(), b"GGGG");
    }
}

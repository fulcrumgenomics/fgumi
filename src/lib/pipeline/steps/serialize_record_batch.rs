//! `SerializeRecordBatch`: `Parallel + ByItemOrdinal` step that frames
//! a [`RecordBatch`] (sorted record bodies from `SortMerge`) as
//! on-disk BAM record framing (4-byte little-endian `block_size`
//! prefix + record body) inside a single [`DecompressedBlock`]
//! ready for [`BgzfCompress`].
//!
//! Sibling of [`crate::pipeline::steps::serialize::SerializeBamRecords`]
//! (which takes [`BamTemplateBatch`] as input — useful when the
//! upstream is template-grouped). Both steps share the same
//! per-record BAM framing logic via [`frame_record_into`]; the
//! difference is purely the input shape. This variant exists so the
//! runall fused `AAM → … → SortMerge → … → Write` chain can pipe
//! sorted records straight to BGZF compression without going
//! through a `RecordBatch → BamTemplateBatch` re-grouping shim.

use std::io;

use crate::pipeline::core::Unpushed;
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::outputs::OrderedBytesSingle;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};
use crate::pipeline::steps::types::{DecompressedBlock, RecordBatch};

/// Per-worker serialization scratch capacity. Matches
/// [`crate::pipeline::steps::serialize::SerializeBamRecords`]'s
/// 1 MiB so the two paths share a single mimalloc size class.
const SERIALIZE_SCRATCH_CAPACITY: usize = 1024 * 1024;

/// Append one BAM-framed record (4-byte little-endian
/// `block_size` + body) into `scratch`. Shared by both
/// [`SerializeRecordBatch`] (this module) and
/// [`crate::pipeline::steps::serialize::SerializeBamRecords`]
/// (sibling module) so the framing logic doesn't drift between
/// the two emitters.
///
/// # Errors
///
/// Returns `io::ErrorKind::InvalidData` if `body.len()` exceeds
/// `u32::MAX` — unreachable in practice (BAM body size is u32-
/// bounded and `tuning.per_step_byte_limit` is 4 MiB by default).
pub(crate) fn frame_record_into(scratch: &mut Vec<u8>, body: &[u8]) -> io::Result<()> {
    let block_size = u32::try_from(body.len()).map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("record exceeds u32 BAM block_size: {}", body.len()),
        )
    })?;
    scratch.extend_from_slice(&block_size.to_le_bytes());
    scratch.extend_from_slice(body);
    Ok(())
}

pub struct SerializeRecordBatch {
    held: HeldSlot<Unpushed<DecompressedBlock>>,
    output_byte_limit: u64,
    /// Per-worker output scratch — see `BgzfDecompress::output_scratch`
    /// for the `mem::replace` + fixed-size-class rationale.
    output_scratch: Vec<u8>,
}

impl SerializeRecordBatch {
    #[must_use]
    pub fn new(output_byte_limit: u64) -> Self {
        Self {
            held: HeldSlot::new(),
            output_byte_limit,
            output_scratch: Vec::with_capacity(SERIALIZE_SCRATCH_CAPACITY),
        }
    }
}

impl Clone for SerializeRecordBatch {
    fn clone(&self) -> Self {
        Self {
            held: HeldSlot::new(),
            output_byte_limit: self.output_byte_limit,
            output_scratch: Vec::with_capacity(SERIALIZE_SCRATCH_CAPACITY),
        }
    }
}

impl Step for SerializeRecordBatch {
    type Input = RecordBatch;
    type Outputs = OrderedBytesSingle<DecompressedBlock>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "SerializeRecordBatch",
            kind: StepKind::Parallel,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn new_worker_copy(&self) -> Self {
        self.clone()
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

        for body in batch.iter_record_bytes() {
            frame_record_into(&mut self.output_scratch, body)?;
        }
        let bytes = std::mem::replace(
            &mut self.output_scratch,
            Vec::with_capacity(SERIALIZE_SCRATCH_CAPACITY),
        );

        // `batch_serial` here is whatever `SortMerge` assigned (a
        // fresh sequence starting at 0 per merge driver), NOT
        // AAM's serial. That's the intended ordering for downstream
        // `BranchOrdering::ByItemOrdinal` consumers — BGZF blocks
        // need to be written in sort-output order, not
        // input-batch order.
        let out = DecompressedBlock { batch_serial: batch.batch_serial(), bytes };
        match ctx.outputs.push(out) {
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
    use crate::pipeline::steps::parse::bam::parse_records;
    use crate::pipeline::steps::types::RecordBatchBuilder;

    #[test]
    fn profile_advertises_parallel_byordinal() {
        let s = SerializeRecordBatch::new(1024);
        let p = s.profile();
        assert_eq!(p.name, "SerializeRecordBatch");
        assert_eq!(p.kind, StepKind::Parallel);
        assert_eq!(p.branch_ordering, vec![BranchOrdering::ByItemOrdinal]);
    }

    /// Round-trip: build a `RecordBatch` with two records, manually
    /// run the serializer's framing logic, and verify `parse_records`
    /// recovers the original record bodies in order.
    #[test]
    fn serialize_then_parse_roundtrip() {
        let r1: Vec<u8> = vec![0xAA; 8];
        let r2: Vec<u8> = vec![0xBB; 12];
        let mut builder = RecordBatchBuilder::with_capacity(7, r1.len() + r2.len(), 2);
        builder.push_record_bytes(&r1);
        builder.push_record_bytes(&r2);
        let rb = builder.build();

        // Drive the serializer's framing logic manually (mirrors the
        // step body; the step-trait integration is exercised by the
        // higher-level pipeline tests).
        let mut bytes = Vec::new();
        for body in rb.iter_record_bytes() {
            let block_size = u32::try_from(body.len()).unwrap();
            bytes.extend_from_slice(&block_size.to_le_bytes());
            bytes.extend_from_slice(body);
        }
        let parsed = parse_records(&bytes).expect("parse");
        assert_eq!(parsed.len(), 2);
        assert_eq!(parsed[0].as_ref(), &r1[..]);
        assert_eq!(parsed[1].as_ref(), &r2[..]);
    }
}

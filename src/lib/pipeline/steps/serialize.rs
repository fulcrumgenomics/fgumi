//! `SerializeBamRecords`: `Parallel + ByItemOrdinal` step that frames
//! `RawRecord`s as on-disk BAM record framing (4-byte little-endian
//! `block_size` prefix + record body) into a single
//! `DecompressedBlock`, ready for `BgzfCompress`. Input is
//! `BamTemplateBatch` (templates of records). Inverse of
//! `ParseBamRecords + GroupBam`. Used by the group/consensus chains.
//!
//! Output `batch_serial` mirrors the input's.

use std::io;

use crate::pipeline::core::Unpushed;
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::outputs::OrderedBytesSingle;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};
use crate::pipeline::steps::types::{BamTemplateBatch, DecompressedBlock};

/// Per-worker serialization scratch capacity. Sized to the typical
/// 500-template batch (~ 256 B/record × ~2 records/template + framing).
/// Mirrors legacy `bam.rs:1561` `SERIALIZATION_BUFFER_CAPACITY` (64 KiB)
/// scaled up to match the new framework's larger batch sizes.
const SERIALIZE_SCRATCH_CAPACITY: usize = 1024 * 1024;

/// Append one BAM-framed record (4-byte little-endian `block_size` + body)
/// into `scratch`. This is the canonical on-disk BAM record framing used by
/// the output path; the terminal standalone-sort path inlines the identical
/// layout in `SortMerge`'s `BlockBuilder` (`fgumi-pipeline-io`), which must
/// stay byte-for-byte in sync with this helper.
///
/// # Errors
///
/// Returns `io::ErrorKind::InvalidData` if `body.len()` exceeds `u32::MAX` —
/// unreachable in practice (BAM body size is u32-bounded and
/// `tuning.per_step_byte_limit` is 4 MiB by default).
fn frame_record_into(scratch: &mut Vec<u8>, body: &[u8]) -> io::Result<()> {
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

/// `Parallel + ByItemOrdinal` serializer. Templates → record-aligned bytes.
pub struct SerializeBamRecords {
    held: HeldSlot<Unpushed<DecompressedBlock>>,
    output_byte_limit: u64,
    /// Per-worker output scratch — see `BgzfDecompress::output_scratch`
    /// for the `mem::replace` + fixed-size-class rationale.
    output_scratch: Vec<u8>,
}

impl SerializeBamRecords {
    #[must_use]
    pub fn new(output_byte_limit: u64) -> Self {
        Self {
            held: HeldSlot::new(),
            output_byte_limit,
            output_scratch: Vec::with_capacity(SERIALIZE_SCRATCH_CAPACITY),
        }
    }
}

impl Clone for SerializeBamRecords {
    fn clone(&self) -> Self {
        Self {
            held: HeldSlot::new(),
            output_byte_limit: self.output_byte_limit,
            output_scratch: Vec::with_capacity(SERIALIZE_SCRATCH_CAPACITY),
        }
    }
}

impl Step for SerializeBamRecords {
    type Input = BamTemplateBatch;
    type Outputs = OrderedBytesSingle<DecompressedBlock>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "SerializeBamRecords",
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

        // Serialize into the per-worker scratch buffer; `mem::replace` it
        // out at the end. Same-size-class recycling via mimalloc's
        // thread-local cache is the goal — see `BgzfDecompress` for the
        // rationale. Per-record framing shared with `SerializeRecordBatch`
        // via [`frame_record_into`].
        for template in &batch.templates {
            for record in &template.records {
                frame_record_into(&mut self.output_scratch, record.as_ref())?;
            }
        }
        let bytes = std::mem::replace(
            &mut self.output_scratch,
            Vec::with_capacity(SERIALIZE_SCRATCH_CAPACITY),
        );

        let out = DecompressedBlock { batch_serial: batch.batch_serial, bytes };
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
    use crate::template::Template;
    use fgumi_raw_bam::RawRecord;

    #[test]
    fn profile_advertises_parallel_byordinal() {
        let s = SerializeBamRecords::new(1024);
        let p = s.profile();
        assert_eq!(p.name, "SerializeBamRecords");
        assert_eq!(p.kind, StepKind::Parallel);
        assert_eq!(p.branch_ordering, vec![BranchOrdering::ByItemOrdinal]);
    }

    /// Build a minimal `RawRecord` whose bytes match what `parse_records`
    /// would emit (a body of `len` arbitrary bytes; serialization
    /// prepends the 4-byte `block_size`).
    fn raw_of_len(len: usize, fill: u8) -> RawRecord {
        vec![fill; len].into()
    }

    #[test]
    fn serialize_then_parse_roundtrip() {
        // Two templates, one record each. After serializing we expect
        // [block_size=8, 8 bytes][block_size=12, 12 bytes].
        let r1 = raw_of_len(8, 0xAA);
        let r2 = raw_of_len(12, 0xBB);
        let mut t1 = Template::new(b"q1".to_vec());
        t1.records.push(r1.clone());
        let mut t2 = Template::new(b"q2".to_vec());
        t2.records.push(r2.clone());
        let input = BamTemplateBatch::new(7, vec![t1, t2]);

        // Drive the serializer manually (mirrors the closure path; the
        // step trait integration is exercised by `tests.rs` via `run`).
        let record_count: usize = input.templates.iter().map(|t| t.records.len()).sum();
        let mut bytes = Vec::with_capacity(input.total_bytes + 4 * record_count);
        for template in &input.templates {
            for record in &template.records {
                let block_size = u32::try_from(record.len()).unwrap();
                bytes.extend_from_slice(&block_size.to_le_bytes());
                bytes.extend_from_slice(record.as_ref());
            }
        }
        let parsed = parse_records(&bytes).expect("parse");
        assert_eq!(parsed.len(), 2);
        assert_eq!(parsed[0].as_ref(), r1.as_ref());
        assert_eq!(parsed[1].as_ref(), r2.as_ref());
    }

    /// `frame_record_into` writes `[u32 LE block_size][body]` and is the
    /// canonical layout the terminal-sort `BlockBuilder` must reproduce
    /// byte-for-byte; verify the exact prefix + body layout here.
    #[test]
    fn frame_record_into_writes_le_prefix_then_body() {
        let mut scratch = Vec::new();
        frame_record_into(&mut scratch, &[0xAA; 8]).unwrap();
        frame_record_into(&mut scratch, &[0xBB; 12]).unwrap();
        let mut expected = Vec::new();
        expected.extend_from_slice(&8u32.to_le_bytes());
        expected.extend_from_slice(&[0xAA; 8]);
        expected.extend_from_slice(&12u32.to_le_bytes());
        expected.extend_from_slice(&[0xBB; 12]);
        assert_eq!(scratch, expected);
        // Round-trips through the BAM record parser.
        let parsed = parse_records(&scratch).expect("parse");
        assert_eq!(parsed.len(), 2);
        assert_eq!(parsed[0].as_ref(), &[0xAA; 8]);
        assert_eq!(parsed[1].as_ref(), &[0xBB; 12]);
    }
}

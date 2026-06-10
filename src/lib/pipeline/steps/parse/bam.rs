//! `ParseBamRecords` mid-step. `Parallel + ByItemOrdinal`. Walks the
//! record-aligned bytes from `FindBamBoundaries` and produces a
//! `RecordBatch` that shares the input block's bytes (zero per-record
//! allocation — the block's `Vec<u8>` is moved into the
//! `RecordBatch`'s backing buffer, plus one `Vec<(u32, u32)>` of
//! `(start, end)` ranges).

use std::io;

use fgumi_raw_bam::RawRecord;

use crate::pipeline::core::Unpushed;
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::outputs::OrderedBytesSingle;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};
use crate::pipeline::steps::types::{DecompressedBlock, RecordBatch};

/// `Parallel + ByItemOrdinal` parser. Walks record-aligned bytes (one
/// record = `block_size: u32 LE` + body) and emits a `RecordBatch`.
///
/// The bytes coming in MUST be record-aligned (post-`FindBamBoundaries`).
/// If they aren't, parsing returns an error.
pub struct ParseBamRecords {
    held: HeldSlot<Unpushed<RecordBatch>>,
    output_byte_limit: u64,
}

impl ParseBamRecords {
    #[must_use]
    pub fn new(output_byte_limit: u64) -> Self {
        Self { held: HeldSlot::new(), output_byte_limit }
    }
}

impl Clone for ParseBamRecords {
    fn clone(&self) -> Self {
        Self { held: HeldSlot::new(), output_byte_limit: self.output_byte_limit }
    }
}

impl Step for ParseBamRecords {
    type Input = DecompressedBlock;
    type Outputs = OrderedBytesSingle<RecordBatch>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "ParseBamRecords",
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
                    // Return `Contention` (not `NoProgress`) so the framework
                    // doesn't observe drain on this worker while we still
                    // have a held item to push — observing drain here
                    // would mark the worker `Skip` and silently drop the
                    // held item.
                    return Ok(StepOutcome::Contention);
                }
            }
        }

        let Some(block) = ctx.input.pop() else {
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
        let DecompressedBlock { batch_serial, bytes } = block;

        // Zero per-record allocation: parse only the body byte ranges,
        // then hand the block's `Vec<u8>` to the `RecordBatch` as the
        // shared backing buffer. The ranges Vec is the only new
        // allocation per batch (plus the move of `bytes` into the
        // batch).
        let ranges = parse_record_ranges(&bytes)?;
        let batch = RecordBatch::from_parsed(batch_serial, bytes, ranges);

        match ctx.outputs.push(batch) {
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

/// Walk record-aligned bytes and produce `(body_start, body_end)` ranges
/// per record. Each record is preceded by a 4-byte little-endian
/// `block_size` prefix; the body is `block_size` bytes that follow.
/// Range values are byte offsets into the caller's buffer.
///
/// Used by [`ParseBamRecords`]'s zero-alloc path: the input block's
/// `Vec<u8>` is moved into the output `RecordBatch` as a shared backing
/// buffer, with these ranges identifying each record's body.
///
/// # Errors
///
/// Returns `Err` if a record runs past the buffer end or if there are
/// trailing partial-record bytes.
pub(crate) fn parse_record_ranges(bytes: &[u8]) -> io::Result<Vec<(u32, u32)>> {
    let mut ranges = Vec::new();
    let mut cursor = 0usize;
    while cursor + 4 <= bytes.len() {
        let block_size = u32::from_le_bytes([
            bytes[cursor],
            bytes[cursor + 1],
            bytes[cursor + 2],
            bytes[cursor + 3],
        ]) as usize;
        let body_start = cursor + 4;
        let record_end = body_start + block_size;
        if record_end > bytes.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "parse_record_ranges: record extends past buffer end \
                     (cursor={cursor}, block_size={block_size}, buffer_len={})",
                    bytes.len()
                ),
            ));
        }
        let start = u32::try_from(body_start).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("parse_record_ranges: body_start {body_start} exceeds u32::MAX"),
            )
        })?;
        let end = u32::try_from(record_end).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("parse_record_ranges: body_end {record_end} exceeds u32::MAX"),
            )
        })?;
        ranges.push((start, end));
        cursor = record_end;
    }
    if cursor != bytes.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "parse_record_ranges: trailing partial record \
                 (cursor={cursor}, buffer_len={})",
                bytes.len()
            ),
        ));
    }
    Ok(ranges)
}

/// Walk record-aligned bytes and produce `Vec<RawRecord>` (one heap
/// allocation per record). Convenience wrapper for test paths and the
/// `SerializeBamRecords` round-trip; the production
/// `ParseBamRecords` step uses [`parse_record_ranges`] for the
/// zero-alloc path.
///
/// # Errors
///
/// Returns `Err` if a record runs past the buffer end or
/// `block_size` is too small for a valid BAM record.
pub(crate) fn parse_records(bytes: &[u8]) -> io::Result<Vec<RawRecord>> {
    let mut records = Vec::new();
    let mut cursor = 0;
    while cursor + 4 <= bytes.len() {
        let block_size = u32::from_le_bytes([
            bytes[cursor],
            bytes[cursor + 1],
            bytes[cursor + 2],
            bytes[cursor + 3],
        ]) as usize;
        let record_end = cursor + 4 + block_size;
        if record_end > bytes.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "ParseBamRecords: record extends past buffer end \
                     (cursor={cursor}, block_size={block_size}, buffer_len={})",
                    bytes.len()
                ),
            ));
        }
        let record: RawRecord = bytes[cursor + 4..record_end].to_vec().into();
        records.push(record);
        cursor = record_end;
    }
    if cursor != bytes.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "ParseBamRecords: trailing partial record \
                 (cursor={cursor}, buffer_len={})",
                bytes.len()
            ),
        ));
    }
    Ok(records)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn profile_advertises_parallel_byordinal() {
        let s = ParseBamRecords::new(1024);
        let p = s.profile();
        assert_eq!(p.name, "ParseBamRecords");
        assert_eq!(p.kind, StepKind::Parallel);
        assert_eq!(p.branch_ordering, vec![BranchOrdering::ByItemOrdinal]);
    }

    #[test]
    fn parse_records_decodes_block_size_prefix() {
        // Build two synthetic "records": [block_size=8, 8 bytes of payload], twice.
        let mut bytes = Vec::new();
        for &p in &[1u8, 2u8] {
            let payload = vec![p; 8];
            bytes.extend_from_slice(&u32::try_from(payload.len()).unwrap().to_le_bytes());
            bytes.extend_from_slice(&payload);
        }
        let records = parse_records(&bytes).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].len(), 8);
        assert_eq!(records[1].len(), 8);
        assert_eq!(records[0].as_ref(), &[1u8; 8]);
        assert_eq!(records[1].as_ref(), &[2u8; 8]);
    }
}

//! `ParseBamRecords` mid-step. `Parallel + ByItemOrdinal`. Walks the
//! record-aligned bytes from `FindBamBoundaries` and produces a
//! `RecordBatch` that shares the input block's bytes (zero per-record
//! allocation — the block's `Vec<u8>` is moved into the
//! `RecordBatch`'s backing buffer, plus one `Vec<(u32, u32)>` of
//! `(start, end)` ranges).

use std::io;
use std::sync::Arc;

use fgumi_bam_io::ProgressTracker;
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
    /// Records parsed across all Parallel workers (shared `Arc`), logged every
    /// 1M under `RUST_LOG=info`. Compared against `SortBuffer`'s ingest rate,
    /// this splits the read-supply (decode chain) cost from `SortBuffer.push`.
    parse_progress: Arc<ProgressTracker>,
}

impl ParseBamRecords {
    #[must_use]
    pub fn new(output_byte_limit: u64) -> Self {
        Self {
            held: HeldSlot::new(),
            output_byte_limit,
            parse_progress: Arc::new(
                ProgressTracker::new("Parse records").with_interval(1_000_000),
            ),
        }
    }
}

impl Clone for ParseBamRecords {
    fn clone(&self) -> Self {
        // Share the parse counter across workers so the count is global.
        Self {
            held: HeldSlot::new(),
            output_byte_limit: self.output_byte_limit,
            parse_progress: Arc::clone(&self.parse_progress),
        }
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
        self.parse_progress.log_if_needed(ranges.len() as u64);
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

/// Smallest number of bytes a single framed BAM record can occupy: the 4-byte
/// little-endian `block_size` prefix plus the 32-byte fixed-field header
/// ([`fgumi_raw_bam::MIN_BAM_RECORD_LEN`]). Dividing a buffer length by this
/// yields an upper bound on the records it can contain, which is what the two
/// parsers below use to pre-size their output.
const MIN_FRAMED_RECORD_LEN: usize = 4 + fgumi_raw_bam::MIN_BAM_RECORD_LEN;

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
    // Pre-size from the smallest record a well-formed stream can hold (a 4-byte
    // `block_size` prefix plus the 32-byte fixed header), so
    // `len / MIN_FRAMED_RECORD_LEN` is an upper bound on the record count and
    // the per-record pushes below never regrow.
    //
    // That bound over-reserves several times over for real records (60-200
    // bytes), and this Vec is RETAINED: `RecordBatch::from_parsed` stores it
    // verbatim and `RecordBatch::heap_size` bills `ranges.capacity()` against
    // the pipeline's byte budget. Leaving the slack in place would inflate every
    // batch's claim (~16% for a 64 KiB block of 150-byte records) and cut how
    // many batches fit in flight. So the slack is released once at the end —
    // see the `shrink_to_fit` below.
    let mut ranges = Vec::with_capacity(bytes.len() / MIN_FRAMED_RECORD_LEN + 1);
    let mut cursor = 0usize;
    while cursor + 4 <= bytes.len() {
        let block_size = u32::from_le_bytes([
            bytes[cursor],
            bytes[cursor + 1],
            bytes[cursor + 2],
            bytes[cursor + 3],
        ]) as usize;
        let body_start = cursor + 4;
        // `block_size` is read straight out of the stream, so fold it in with
        // checked arithmetic before it is compared against the buffer or used to
        // slice. Unreachable on 64-bit (a `u32` widened to `usize` cannot
        // overflow this sum), but the framing rule is that every length-derived
        // range is checked before use, not that it happens to fit today.
        let Some(record_end) = body_start.checked_add(block_size) else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "parse_record_ranges: record end overflows \
                     (cursor={cursor}, block_size={block_size})"
                ),
            ));
        };
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
    // Return the table at exactly its populated size. This costs one copy of
    // `len * 8` bytes, against the ~9 regrow-and-copy rounds the pre-size above
    // avoided — and it keeps the batch's `heap_size` honest, so byte-bounded
    // admission bills only memory that is actually needed.
    ranges.shrink_to_fit();
    Ok(ranges)
}

/// Walk record-aligned bytes and produce `Vec<RawRecord>` (one heap
/// allocation per record, via `to_vec()`).
///
/// Used by:
/// - the `SerializeBamRecords` round-trip and test paths;
/// - the production [`DecodeRecords`](super::decode::DecodeRecords) step, whose
///   downstream `DecodedRecord::from_raw_bytes` takes an **owned** `RawRecord`
///   (it stores the bytes, unlike `ParseBamRecords` which only emits ranges into
///   the shared block) — so the production decode path does allocate one `Vec`
///   per record here. This is inherent to the owned-per-record `DecodedRecord`
///   representation; sharing the block buffer across a batch's records (e.g.
///   `Arc<[u8]>` + per-record range) is a `DecodedRecord` representation change
///   deferred behind a benchmark.
///
/// The other production parser, the `ParseBamRecords` step, instead uses
/// [`parse_record_ranges`] for its zero-alloc path.
///
/// # Errors
///
/// Returns `Err` if a record runs past the buffer end, or if the buffer ends
/// with a trailing partial record (the final `block_size` does not consume the
/// remaining bytes exactly). It does NOT enforce a per-record minimum body
/// size, so a `block_size` below the BAM fixed-field minimum is accepted here
/// and rejected (if at all) by later decode of the record body.
pub(crate) fn parse_records(bytes: &[u8]) -> io::Result<Vec<RawRecord>> {
    // Same upper-bound pre-size as `parse_record_ranges` — this sibling walks
    // the identical framing and pushes once per record, so it has the identical
    // regrowth cost.
    //
    // Deliberately NOT shrunk to fit, unlike that sibling. This `Vec` is
    // consumed immediately by the caller (`DecodeRecords` maps it into
    // `Vec<DecodedRecord>` and drops it), so its capacity never reaches a queued
    // item and is never billed against the pipeline's byte budget. Paying a copy
    // to shrink a Vec that is about to be dropped would be pure cost.
    let mut records = Vec::with_capacity(bytes.len() / MIN_FRAMED_RECORD_LEN + 1);
    let mut cursor = 0;
    while cursor + 4 <= bytes.len() {
        let block_size = u32::from_le_bytes([
            bytes[cursor],
            bytes[cursor + 1],
            bytes[cursor + 2],
            bytes[cursor + 3],
        ]) as usize;
        // Checked for the same reason as `parse_record_ranges` above.
        let Some(record_end) = cursor.checked_add(4).and_then(|s| s.checked_add(block_size)) else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "parse_records: record end overflows \
                     (cursor={cursor}, block_size={block_size})"
                ),
            ));
        };
        if record_end > bytes.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "parse_records: record extends past buffer end \
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
                "parse_records: trailing partial record \
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
    use rstest::rstest;

    #[test]
    fn profile_advertises_parallel_byordinal() {
        let s = ParseBamRecords::new(1024);
        let p = s.profile();
        assert_eq!(p.name, "ParseBamRecords");
        assert_eq!(p.kind, StepKind::Parallel);
        assert_eq!(p.branch_ordering, vec![BranchOrdering::ByItemOrdinal]);
    }

    /// Frame `payload` as one record: `block_size: u32 LE` + body.
    fn framed(payload: &[u8]) -> Vec<u8> {
        let mut out = u32::try_from(payload.len()).expect("fits u32").to_le_bytes().to_vec();
        out.extend_from_slice(payload);
        out
    }

    /// `parse_record_ranges` is the zero-alloc sibling of `parse_records`: it
    /// must identify the same record bodies, expressed as offsets into the
    /// caller's buffer rather than as copies.
    #[test]
    fn parse_record_ranges_brackets_each_record_body() {
        let payloads: [&[u8]; 3] = [&[1u8; 8], &[2u8; 0], &[3u8; 20]];
        let bytes: Vec<u8> = payloads.iter().flat_map(|p| framed(p)).collect();

        let ranges = parse_record_ranges(&bytes).expect("record-aligned bytes parse");

        assert_eq!(ranges, vec![(4, 12), (16, 16), (20, 40)]);
        for (i, (start, end)) in ranges.iter().enumerate() {
            assert_eq!(
                &bytes[*start as usize..*end as usize],
                payloads[i],
                "range {i} must bracket exactly that record's body",
            );
        }
    }

    /// The returned table is retained inside `RecordBatch`, whose `heap_size`
    /// bills `ranges.capacity()` against the pipeline's byte budget. So the
    /// upper-bound pre-size must not survive into the returned Vec: excess
    /// reserved slots would inflate every batch's claim and cut how many fit in
    /// flight, for memory no record needs.
    #[test]
    fn parse_record_ranges_returns_a_table_with_no_excess_capacity() {
        // Records well above the 36-byte minimum, so the pre-size over-reserves
        // by several times and the slack is real rather than incidental.
        let payloads: Vec<Vec<u8>> = (0..40).map(|_| vec![0u8; 160]).collect();
        let refs: Vec<&[u8]> = payloads.iter().map(Vec::as_slice).collect();
        let bytes: Vec<u8> = refs.iter().flat_map(|p| framed(p)).collect();

        // Sanity: the pre-size really would over-reserve here.
        let upper_bound = bytes.len() / MIN_FRAMED_RECORD_LEN + 1;
        assert!(
            upper_bound > payloads.len(),
            "fixture must over-reserve for this test to mean anything \
             ({upper_bound} <= {})",
            payloads.len(),
        );

        let ranges = parse_record_ranges(&bytes).expect("parses");
        assert_eq!(ranges.len(), payloads.len());
        assert_eq!(
            ranges.capacity(),
            ranges.len(),
            "the retained ranges table must carry no reserved slack — it is \
             billed by RecordBatch::heap_size against the byte budget",
        );
    }

    /// Empty input holds no records — a legal batch, not an error.
    #[test]
    fn parse_record_ranges_accepts_an_empty_buffer() {
        assert!(parse_record_ranges(&[]).expect("empty is legal").is_empty());
        assert!(parse_records(&[]).expect("empty is legal").is_empty());
    }

    /// Both parsers must reject malformed input rather than silently truncating
    /// it — a record whose declared `block_size` runs past the buffer, and a
    /// trailing fragment too short to frame. Silently dropping either would turn
    /// stream corruption into missing records.
    #[rstest]
    #[case::block_size_overruns_the_buffer(&[0x40, 0, 0, 0, 1, 2, 3])]
    #[case::trailing_bytes_cannot_frame_a_record(&[1, 2, 3])]
    #[case::single_trailing_byte(&[9])]
    fn both_parsers_reject_a_malformed_tail(#[case] tail: &[u8]) {
        let mut bytes = framed(&[1u8; 8]);
        bytes.extend_from_slice(tail);

        let ranges_err = parse_record_ranges(&bytes).expect_err("must reject a malformed tail");
        assert_eq!(ranges_err.kind(), io::ErrorKind::InvalidData);

        let records_err = parse_records(&bytes).expect_err("must reject a malformed tail");
        assert_eq!(records_err.kind(), io::ErrorKind::InvalidData);
    }

    /// The two parsers must agree on where records are: same count, and the
    /// bytes `parse_records` copies must equal the slices `parse_record_ranges`
    /// points at. They feed different steps over the same stream, so a
    /// disagreement would mean the two decode paths see different records.
    #[test]
    fn the_two_parsers_agree_on_record_boundaries() {
        let payloads: Vec<Vec<u8>> = (0..16u8).map(|i| vec![i; usize::from(i) + 1]).collect();
        let bytes: Vec<u8> = payloads.iter().flat_map(|p| framed(p)).collect();

        let ranges = parse_record_ranges(&bytes).expect("parses");
        let records = parse_records(&bytes).expect("parses");

        assert_eq!(ranges.len(), records.len());
        for (range, record) in ranges.iter().zip(&records) {
            assert_eq!(&bytes[range.0 as usize..range.1 as usize], record.as_ref());
        }
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

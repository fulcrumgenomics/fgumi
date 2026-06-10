//! `ParseFastqChunks` mid-step. `Parallel + ByItemOrdinal`. Parses the raw,
//! whole-record-aligned FASTQ bytes of a `FastqRawChunk` into a
//! `FastqChunkBatch` of `FastqRecord`s.
//!
//! This step exists to lift the FASTQ-parse cost off the single
//! `ReadFastqInputs` reader thread: `ReadFastqInputs` now only reads and
//! gzip-decompresses bytes (serial, I/O-bound), and one `ParseFastqChunks`
//! worker per pipeline thread parses the chunks concurrently. The framework
//! reorders the parallel output by each chunk's globally-unique `ordinal`,
//! which `ParseFastqChunks` faithfully propagates from input to output.

use std::io;

use crate::fastq_parse::FastqRecord;
use crate::pipeline::core::Unpushed;
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::item::HeapSize;
use crate::pipeline::core::outputs::OrderedBytesSingle;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};

use super::read_fastq::{FastqChunkBatch, FastqRawChunk};

/// Parse the whole-record-aligned raw FASTQ bytes of one chunk into a
/// `Vec<FastqRecord>`, returning the records and their total heap size.
///
/// `data` is expected to contain a whole number of 4-line FASTQ records, each
/// terminated by `\n` (as produced by `read_fastq_raw_bytes_from_bufread`).
/// Record boundaries are located with the SIMD lexer (`find_record_offsets`),
/// matching the legacy pipeline — a single vectorized pass over the chunk
/// instead of a per-byte scalar scan.
fn parse_fastq_records(data: &[u8]) -> io::Result<(Vec<FastqRecord>, usize)> {
    // `find_record_offsets` returns `[0, end_1, .., end_n]` where `end_n` is one
    // byte past the final *complete* record; trailing incomplete bytes are
    // excluded. The reader guarantees whole-record-aligned chunks, so the last
    // offset must equal `data.len()` — anything else is a truncated record.
    let offsets = fgumi_simd_fastq::find_record_offsets(data);
    let last = offsets.last().copied().unwrap_or(0);
    if last != data.len() {
        return Err(io::Error::other(format!(
            "ParseFastqChunks: truncated FASTQ record ({} trailing bytes after the last complete record)",
            data.len() - last
        )));
    }

    let num_records = offsets.len().saturating_sub(1);
    let mut records = Vec::with_capacity(num_records);
    let mut total_bytes = 0usize;
    for window in offsets.windows(2) {
        let record = FastqRecord::from_slice(&data[window[0]..window[1]])?;
        total_bytes += record.heap_size();
        records.push(record);
    }
    Ok((records, total_bytes))
}

/// `Parallel + ByItemOrdinal` FASTQ parser. Each worker is independent and
/// stateless beyond the held slot — a `FastqRawChunk` in, a `FastqChunkBatch`
/// out, with the globally-unique `ordinal` (and the `stream_idx` /
/// `chunk_serial` join keys) preserved verbatim.
pub struct ParseFastqChunks {
    held: HeldSlot<Unpushed<FastqChunkBatch>>,
    output_byte_limit: u64,
}

impl ParseFastqChunks {
    #[must_use]
    pub fn new(output_byte_limit: u64) -> Self {
        Self { held: HeldSlot::new(), output_byte_limit }
    }
}

impl Clone for ParseFastqChunks {
    fn clone(&self) -> Self {
        Self { held: HeldSlot::new(), output_byte_limit: self.output_byte_limit }
    }
}

impl Step for ParseFastqChunks {
    type Input = FastqRawChunk;
    type Outputs = OrderedBytesSingle<FastqChunkBatch>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "ParseFastqChunks",
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
                    // `Contention` (not `NoProgress`) keeps the worker alive
                    // for retry — `NoProgress` would let the framework mark
                    // this worker `Skip` if input is also drained, silently
                    // dropping the held item.
                    return Ok(StepOutcome::Contention);
                }
            }
        }

        let Some(chunk) = ctx.input.pop() else {
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

        let (records, total_bytes) = parse_fastq_records(&chunk.data)?;
        let batch = FastqChunkBatch::new(
            chunk.ordinal,
            chunk.stream_idx,
            chunk.chunk_serial,
            records,
            total_bytes,
        );

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pipeline::core::item::Ordered;

    #[test]
    fn profile_advertises_parallel_byordinal() {
        let s = ParseFastqChunks::new(1024);
        let p = s.profile();
        assert_eq!(p.name, "ParseFastqChunks");
        assert_eq!(p.kind, StepKind::Parallel);
        assert!(!p.sticky);
        assert_eq!(p.branch_ordering, vec![BranchOrdering::ByItemOrdinal]);
    }

    #[test]
    fn parse_fastq_records_parses_all_records() {
        let data = b"@read1\nACGT\n+\nIIII\n@read2\nTGCA\n+\nJJJJ\n";
        let (records, total_bytes) = parse_fastq_records(data).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name(), b"read1");
        assert_eq!(records[0].sequence(), b"ACGT");
        assert_eq!(records[1].name(), b"read2");
        assert_eq!(records[1].sequence(), b"TGCA");
        assert!(total_bytes > 0);
    }

    #[test]
    fn parse_fastq_records_empty_input() {
        let (records, total_bytes) = parse_fastq_records(b"").unwrap();
        assert!(records.is_empty());
        assert_eq!(total_bytes, 0);
    }

    #[test]
    fn parse_fastq_records_truncated_errors() {
        // Only 2 of 4 lines present.
        let data = b"@read1\nACGT\n";
        let err = parse_fastq_records(data).unwrap_err();
        assert!(err.to_string().contains("truncated"), "got: {err}");
    }

    /// The core correctness gate: parsing a `FastqRawChunk` must preserve the
    /// globally-unique `ordinal` and the `(stream_idx, chunk_serial)` join
    /// keys verbatim on the produced `FastqChunkBatch`.
    #[test]
    fn parse_preserves_ordinal_stream_and_serial() {
        let data = b"@read1\nACGT\n+\nIIII\n".to_vec();
        let chunk = FastqRawChunk { ordinal: 17, stream_idx: 1, chunk_serial: 4, data };

        let (records, total_bytes) = parse_fastq_records(&chunk.data).unwrap();
        let batch = FastqChunkBatch::new(
            chunk.ordinal,
            chunk.stream_idx,
            chunk.chunk_serial,
            records,
            total_bytes,
        );

        assert_eq!(batch.ordinal(), 17, "framework reorder key must be the unique ordinal");
        assert_eq!(batch.stream_idx, 1);
        assert_eq!(batch.chunk_serial, 4, "zip join key must be preserved");
        assert_eq!(batch.records.len(), 1);
        assert_eq!(batch.records[0].name(), b"read1");
    }
}

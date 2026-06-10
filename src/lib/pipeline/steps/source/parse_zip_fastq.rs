//! `ParseAndZipFastq` step: `Parallel + ByItemOrdinal`. Parses the raw,
//! whole-record-aligned bytes of a [`PairedRawFastqBatch`] (one R1 chunk + one
//! R2 chunk) into per-stream `FastqRecord`s and zips them record-by-record into
//! a [`FastqTemplateBatch`].
//!
//! This step exists to lift BOTH the FASTQ parse and the template build off the
//! serial source path and fan them across worker threads. It mirrors the legacy
//! pipeline's `Decode` stage: the cheap chunk-level pairing already happened
//! serially in [`PairRawFastq`] (which minted a globally-unique `ordinal`), and
//! each worker here parses + zips one paired batch independently. The framework
//! reorders the parallel output by each batch's `ordinal`, which this step
//! propagates verbatim onto the emitted `FastqTemplateBatch`.
//!
//! Used only for the common paired-end (N == 2) case. The N == 1 and N >= 3
//! paths keep the separate `ParseFastqChunks` → `ZipFastqRecords` structure.

use std::io;

use crate::fastq_parse::{FastqRecord, strip_read_suffix};
use crate::grouper::FastqTemplate;
use crate::pipeline::core::Unpushed;
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::outputs::OrderedBytesSingle;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};

use super::pair_fastq::PairedRawFastqBatch;
use super::zip_fastq::FastqTemplateBatch;

/// Parse the whole-record-aligned raw FASTQ bytes of one chunk into a
/// `Vec<FastqRecord>`.
///
/// `data` is expected to contain a whole number of 4-line FASTQ records, each
/// terminated by `\n` (as produced by `read_fastq_raw_bytes_from_bufread`).
/// Record boundaries are located with the SIMD lexer (`find_record_offsets`),
/// matching the legacy pipeline — a single vectorized pass over the chunk
/// instead of a per-byte scalar scan.
fn parse_fastq_records(data: &[u8]) -> io::Result<Vec<FastqRecord>> {
    // `find_record_offsets` returns `[0, end_1, .., end_n]` where `end_n` is one
    // byte past the final *complete* record; trailing incomplete bytes are
    // excluded. The reader guarantees whole-record-aligned chunks, so the last
    // offset must equal `data.len()` — anything else is a truncated record.
    let offsets = fgumi_simd_fastq::find_record_offsets(data);
    let last = offsets.last().copied().unwrap_or(0);
    if last != data.len() {
        return Err(io::Error::other(format!(
            "ParseAndZipFastq: truncated FASTQ record ({} trailing bytes after the last complete record)",
            data.len() - last
        )));
    }

    let num_records = offsets.len().saturating_sub(1);
    let mut records = Vec::with_capacity(num_records);
    for window in offsets.windows(2) {
        records.push(FastqRecord::from_slice(&data[window[0]..window[1]])?);
    }
    Ok(records)
}

/// Zip the per-stream records of one paired batch into `FastqTemplate`s.
///
/// Record `i` from stream A and record `i` from stream B form one template; the
/// names must match after [`strip_read_suffix`]. This mirrors the legacy
/// `ZipFastqRecords::try_emit_complete` logic verbatim — including its error
/// semantics — so output (and failure behavior) is byte-identical to the
/// `ParseFastqChunks` → `ZipFastqRecords` path the N != 2 chains still use.
///
/// Strict, fgbio-consistent mismatch handling (see `FastqSource.zipped`): a
/// record-count difference between the streams is an explicit "out of sync"
/// error, and a read name that disagrees after suffix-stripping is a "read name
/// mismatch" error. The count check is the same one `ZipFastqRecords` applies
/// on the N>=3 path, so both paths report identically.
fn zip_records(
    records_a: Vec<FastqRecord>,
    records_b: Vec<FastqRecord>,
    chunk_serial: u64,
) -> io::Result<Vec<FastqTemplate>> {
    // Paired FASTQ is positionally parallel: every stream must contribute the
    // same number of records for this chunk. Reject a mismatch explicitly
    // rather than letting it surface as an incidental name mismatch from the
    // back-pop below.
    if records_a.len() != records_b.len() {
        return Err(io::Error::other(format!(
            "FASTQ sources out of sync at chunk_serial {chunk_serial}: \
             stream 0 has {} records, stream 1 has {}",
            records_a.len(),
            records_b.len(),
        )));
    }

    let n_records = records_a.len();
    let mut streams: [Vec<FastqRecord>; 2] = [records_a, records_b];

    // Pop from the back of each stream (avoids O(n^2) front removal) and
    // reverse at the end to restore the original record order.
    let mut reversed_templates = Vec::with_capacity(n_records);
    for _ in 0..n_records {
        let mut records = Vec::with_capacity(2);
        let mut base_name: Option<Vec<u8>> = None;

        for (stream_idx, stream) in streams.iter_mut().enumerate() {
            // Unreachable after the equal-count check above; kept defensive.
            let Some(record) = stream.pop() else {
                return Err(io::Error::other(format!(
                    "FASTQ sources out of sync: stream {stream_idx} ran short at \
                     chunk_serial {chunk_serial}",
                )));
            };

            let stripped = strip_read_suffix(record.name());
            match &base_name {
                None => base_name = Some(stripped.to_vec()),
                Some(expected) => {
                    if stripped != expected.as_slice() {
                        return Err(io::Error::other(format!(
                            "FASTQ read name mismatch at chunk_serial {}: stream 0 has '{}', \
                             stream {} has '{}'",
                            chunk_serial,
                            String::from_utf8_lossy(expected),
                            stream_idx,
                            String::from_utf8_lossy(stripped),
                        )));
                    }
                }
            }
            records.push(record);
        }

        reversed_templates.push(FastqTemplate { name: base_name.unwrap(), records });
    }

    reversed_templates.reverse();
    Ok(reversed_templates)
}

/// `Parallel + ByItemOrdinal` FASTQ parse-and-zip. Each worker is independent
/// and stateless beyond the held slot — a `PairedRawFastqBatch` in, a
/// `FastqTemplateBatch` out, with the globally-unique `ordinal` preserved
/// verbatim as the output `batch_serial`.
pub struct ParseAndZipFastq {
    held: HeldSlot<Unpushed<FastqTemplateBatch>>,
    output_byte_limit: u64,
}

impl ParseAndZipFastq {
    #[must_use]
    pub fn new(output_byte_limit: u64) -> Self {
        Self { held: HeldSlot::new(), output_byte_limit }
    }
}

impl Clone for ParseAndZipFastq {
    fn clone(&self) -> Self {
        Self { held: HeldSlot::new(), output_byte_limit: self.output_byte_limit }
    }
}

impl Step for ParseAndZipFastq {
    type Input = PairedRawFastqBatch;
    type Outputs = OrderedBytesSingle<FastqTemplateBatch>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "ParseAndZipFastq",
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

        let records_a = parse_fastq_records(&batch.data_a)?;
        let records_b = parse_fastq_records(&batch.data_b)?;
        let templates = zip_records(records_a, records_b, batch.chunk_serial)?;
        // Preserve the serially-minted ordinal so the framework reorder keeps
        // the FASTQ records in input order.
        let template_batch = FastqTemplateBatch::new(batch.ordinal, templates);

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
        let s = ParseAndZipFastq::new(1024);
        let p = s.profile();
        assert_eq!(p.name, "ParseAndZipFastq");
        assert_eq!(p.kind, StepKind::Parallel);
        assert!(!p.sticky);
        assert_eq!(p.branch_ordering, vec![BranchOrdering::ByItemOrdinal]);
    }

    #[test]
    fn parse_fastq_records_parses_all_records() {
        let mut data = record_bytes("read1", "ACGT");
        data.extend(record_bytes("read2", "TGCA"));
        let records = parse_fastq_records(&data).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name(), b"read1");
        assert_eq!(records[1].name(), b"read2");
    }

    #[test]
    fn parse_fastq_records_truncated_errors() {
        let data = b"@read1\nACGT\n";
        let err = parse_fastq_records(data).unwrap_err();
        assert!(err.to_string().contains("truncated"), "got: {err}");
    }

    /// The core zip correctness: record i from A and record i from B form one
    /// template, the base read name (suffix-stripped) is preserved, and R1/R2
    /// land in records[0]/records[1] in input order.
    #[test]
    fn zip_pairs_records_and_strips_suffixes() {
        let mut data_a = record_bytes("read1/1", "ACGT");
        data_a.extend(record_bytes("read2/1", "TGCA"));
        let mut data_b = record_bytes("read1/2", "GGGG");
        data_b.extend(record_bytes("read2/2", "CCCC"));

        let recs_a = parse_fastq_records(&data_a).unwrap();
        let recs_b = parse_fastq_records(&data_b).unwrap();
        let templates = zip_records(recs_a, recs_b, 0).unwrap();

        assert_eq!(templates.len(), 2);
        assert_eq!(templates[0].name, b"read1");
        assert_eq!(templates[0].records.len(), 2);
        assert_eq!(templates[0].records[0].sequence(), b"ACGT");
        assert_eq!(templates[0].records[1].sequence(), b"GGGG");
        assert_eq!(templates[1].name, b"read2");
        assert_eq!(templates[1].records[0].sequence(), b"TGCA");
        assert_eq!(templates[1].records[1].sequence(), b"CCCC");
    }

    #[test]
    fn zip_name_mismatch_errors() {
        let recs_a = parse_fastq_records(&record_bytes("readA/1", "ACGT")).unwrap();
        let recs_b = parse_fastq_records(&record_bytes("readB/2", "GGGG")).unwrap();
        let err = zip_records(recs_a, recs_b, 0).unwrap_err();
        assert!(err.to_string().contains("mismatch"), "expected mismatch error, got: {err}");
    }

    /// Unequal record counts are rejected by an explicit, strict count check
    /// (fgbio-consistent) — distinct from a name mismatch. A has 2 records, B
    /// has 1, so the streams are out of sync regardless of read names.
    #[test]
    fn zip_unequal_counts_errors() {
        let mut data_a = record_bytes("read1/1", "ACGT");
        data_a.extend(record_bytes("read2/1", "TGCA"));
        let recs_a = parse_fastq_records(&data_a).unwrap();
        let recs_b = parse_fastq_records(&record_bytes("read1/2", "GGGG")).unwrap();
        let err = zip_records(recs_a, recs_b, 7).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("out of sync"), "expected count out-of-sync error, got: {msg}");
        assert!(msg.contains("chunk_serial 7"), "should name the chunk_serial, got: {msg}");
    }

    /// A shorter than B is symmetric: the explicit count check rejects it the
    /// same way, naming the differing counts rather than relying on back-pop
    /// name misalignment.
    #[test]
    fn zip_a_shorter_than_b_errors() {
        let recs_a = parse_fastq_records(&record_bytes("read1/1", "ACGT")).unwrap();
        let mut data_b = record_bytes("read1/2", "GGGG");
        data_b.extend(record_bytes("read2/2", "CCCC"));
        let recs_b = parse_fastq_records(&data_b).unwrap();
        let err = zip_records(recs_a, recs_b, 0).unwrap_err();
        assert!(err.to_string().contains("out of sync"), "got: {err}");
    }

    /// The output `FastqTemplateBatch` must carry the input batch's ordinal as
    /// its `batch_serial` so the framework reorder preserves input order.
    #[test]
    fn output_ordinal_is_input_ordinal() {
        let data_a = record_bytes("r/1", "ACGT");
        let data_b = record_bytes("r/2", "GGGG");
        let recs_a = parse_fastq_records(&data_a).unwrap();
        let recs_b = parse_fastq_records(&data_b).unwrap();
        let templates = zip_records(recs_a, recs_b, 0).unwrap();
        let batch = FastqTemplateBatch::new(42, templates);
        assert_eq!(batch.ordinal(), 42);
    }
}

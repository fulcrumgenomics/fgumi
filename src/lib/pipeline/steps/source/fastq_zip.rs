//! Shared K-way FASTQ zip primitives: the [`NRawFastqBatch`] item that carries
//! K aligned raw record-chunks, and the [`zip_streams`] function that turns K
//! per-stream parsed record lists into [`FastqTemplate`]s.
//!
//! These are consumed by the parallel K-way decode front:
//! [`ZipRawFastqK`](super::zip_raw_fastq_k::ZipRawFastqK) (Serial) aligns the K
//! streams' raw chunks by `chunk_serial` and emits an [`NRawFastqBatch`];
//! [`ParseAndZipFastqN`](super::parse_zip_fastq::ParseAndZipFastqN) (Parallel)
//! parses each stream's bytes and calls [`zip_streams`] to build the templates.
//!
//! The zip logic is the K-general form lifted verbatim (behavior-preserving)
//! from the former `ZipFastqRecords::try_emit_complete` — same positional
//! join, same strict record-count / read-name concordance checks, same
//! fgbio-consistent "out of sync" diagnostics — so output is byte-identical to
//! the path it replaces.

use std::io;

use crate::fastq_parse::{FastqRecord, strip_read_suffix};
use crate::grouper::FastqTemplate;
use crate::pipeline::core::item::{HeapSize, Ordered};

/// A batch of zipped [`FastqTemplate`]s carrying its ordering serial, the final
/// product of the FASTQ decode front (consumed by `add_extract`'s `ExtractStep`).
///
/// `batch_serial` is the ordinal the framework reorders by; `total_bytes` is
/// cached at construction for O(1) `HeapSize`.
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

/// K aligned raw (decompressed, unparsed) FASTQ byte chunks — one per input
/// stream, all sharing a `chunk_serial` — plus a freshly-minted globally-unique
/// `ordinal`. The K-stream generalization of the former `PairedRawFastqBatch`
/// (which carried exactly `data_a` + `data_b`).
///
/// `ordinal` is minted serially by
/// [`ZipRawFastqK`](super::zip_raw_fastq_k::ZipRawFastqK) and is the key the
/// framework uses to reorder the `Parallel`
/// [`ParseAndZipFastqN`](super::parse_zip_fastq::ParseAndZipFastqN) output.
/// `chunk_serial` is carried for desync diagnostics only. `streams[i]` is the
/// whole-record-aligned bytes of stream `i` for this row.
pub struct NRawFastqBatch {
    /// Globally-unique monotonic ordinal, minted by `ZipRawFastqK`.
    pub ordinal: u64,
    /// The shared per-stream cycle serial (for desync diagnostics).
    pub chunk_serial: u64,
    /// Raw bytes per stream, index = `stream_idx`, whole-record-aligned.
    pub streams: Vec<Vec<u8>>,
}

impl HeapSize for NRawFastqBatch {
    fn heap_size(&self) -> usize {
        // Allocated capacity, not logical length — the chunk buffers arrive from
        // the readers built with `Vec::with_capacity`, so capacity exceeds length
        // on any short (final/partial) chunk. Matches `FastqRawChunk::heap_size`
        // and the byte-bounded queue budget; measuring `len()` would under-count
        // the memory the output queue actually holds.
        self.streams.iter().map(Vec::capacity).sum::<usize>()
            + self.streams.capacity() * std::mem::size_of::<Vec<u8>>()
    }
}

impl Ordered for NRawFastqBatch {
    fn ordinal(&self) -> u64 {
        self.ordinal
    }
}

/// Zip K per-stream parsed record lists into templates: record `i` from every
/// stream forms one [`FastqTemplate`], with the base read name (suffix-stripped)
/// taken from stream 0 and required to match across all streams.
///
/// K-general lift of the former `ZipFastqRecords::try_emit_complete` /
/// `ParseAndZipFastq::zip_records` logic — identical positional join and
/// identical strict, fgbio-consistent checks:
///
/// - **Unequal per-stream record counts** are an explicit "out of sync" error
///   (naming the stream that ran short), not a silent truncation.
/// - **A read name that disagrees** across streams after
///   [`strip_read_suffix`] is a "read name mismatch" `InvalidData` error.
///
/// `chunk_serial` is only used to label the diagnostics.
///
/// # Errors
///
/// Returns `Other` ("out of sync") on a per-stream record-count mismatch, or
/// `InvalidData` ("read name mismatch") when a template's segments disagree on
/// the base name.
///
/// # Panics
///
/// Panics only on an internal invariant violation: a template row is assembled
/// with no segments (so no base name is established). This cannot happen for
/// `n_streams >= 1`, which is asserted on entry.
pub fn zip_streams(
    mut stream_records: Vec<Vec<FastqRecord>>,
    chunk_serial: u64,
) -> io::Result<Vec<FastqTemplate>> {
    let n_streams = stream_records.len();
    debug_assert!(n_streams >= 1, "zip_streams needs at least one stream");
    let n_records = stream_records[0].len();

    // Positional-parallel: every stream must contribute the same number of
    // records for this chunk. Reject a mismatch explicitly (fgbio parity)
    // rather than letting it surface as an incidental name mismatch below.
    for (stream_idx, recs) in stream_records.iter().enumerate().skip(1) {
        if recs.len() != n_records {
            // Name streams R1/R2/… (1-based) and the shorter as "ended before"
            // the other, matching the serial oracle's wording (issue #773).
            let (ended, before) =
                if recs.len() < n_records { (stream_idx, 0) } else { (0, stream_idx) };
            return Err(io::Error::other(format!(
                "FASTQ sources out of sync at chunk_serial {chunk_serial}: R{} ended before \
                 R{} (R1 has {n_records} record(s), R{} has {} in this chunk)",
                ended + 1,
                before + 1,
                stream_idx + 1,
                recs.len(),
            )));
        }
    }

    // Zip across streams: record i from each stream forms one template. Drain
    // from the back (pop) to avoid O(n^2) front removal; reverse at the end.
    let mut templates = Vec::with_capacity(n_records);
    for _ in 0..n_records {
        let mut records = Vec::with_capacity(n_streams);
        let mut base_name: Option<Vec<u8>> = None;

        for (stream_idx, stream) in stream_records.iter_mut().enumerate() {
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
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!(
                                "FASTQ read name mismatch at chunk_serial {}: stream 0 has \
                                 '{}', stream {} has '{}'",
                                chunk_serial,
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

        templates.push(FastqTemplate { name: base_name.unwrap(), records });
    }

    // Popped from the back, so reverse to restore original record order.
    templates.reverse();
    Ok(templates)
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    fn rec(name: &str, seq: &str) -> FastqRecord {
        let qual: String = std::iter::repeat_n('I', seq.len()).collect();
        FastqRecord::from_slice(format!("@{name}\n{seq}\n+\n{qual}\n").as_bytes()).unwrap()
    }

    #[test]
    fn nraw_batch_heap_size_counts_capacity_and_ordinal() {
        let batch = NRawFastqBatch {
            ordinal: 7,
            chunk_serial: 3,
            streams: vec![vec![0u8; 10], vec![0u8; 7]],
        };
        assert_eq!(batch.ordinal(), 7);
        assert!(batch.heap_size() >= 17);
    }

    /// K=3 positional zip: record i from each of the three streams forms one
    /// template, base name from stream 0, in input order.
    #[test]
    fn zip_streams_joins_three_streams_in_order() {
        let s0 = vec![rec("read0/1", "AAAA"), rec("read1/1", "CCCC")];
        let s1 = vec![rec("read0/2", "GGGG"), rec("read1/2", "TTTT")];
        let s2 = vec![rec("read0", "ACAC"), rec("read1", "GTGT")];
        let templates = zip_streams(vec![s0, s1, s2], 0).unwrap();
        assert_eq!(templates.len(), 2);
        assert_eq!(templates[0].name, b"read0");
        assert_eq!(templates[0].records.len(), 3);
        assert_eq!(templates[0].records[0].sequence(), b"AAAA");
        assert_eq!(templates[0].records[1].sequence(), b"GGGG");
        assert_eq!(templates[0].records[2].sequence(), b"ACAC");
        assert_eq!(templates[1].name, b"read1");
        assert_eq!(templates[1].records.len(), 3);
    }

    /// The out-of-sync diagnostic must name WHICH stream ran short (issue #773),
    /// not merely that a mismatch occurred: a reversed `(ended, before)`
    /// direction is a silent regression that a "contains out of sync" check
    /// would miss. Both branches of the direction conditional are covered.
    #[rstest]
    #[case::second_stream_short(2, 1, "R2 ended before R1")]
    #[case::first_stream_short(1, 2, "R1 ended before R2")]
    fn zip_streams_unequal_counts_names_the_short_stream(
        #[case] n0: usize,
        #[case] n1: usize,
        #[case] expected: &str,
    ) {
        let s0: Vec<_> = (0..n0).map(|i| rec(&format!("read{i}/1"), "AAAA")).collect();
        let s1: Vec<_> = (0..n1).map(|i| rec(&format!("read{i}/2"), "GGGG")).collect();
        let err = zip_streams(vec![s0, s1], 5).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("out of sync"), "got: {msg}");
        assert!(msg.contains("chunk_serial 5"), "got: {msg}");
        assert!(msg.contains(expected), "expected {expected:?} in: {msg}");
    }

    #[test]
    fn zip_streams_name_mismatch_errors() {
        let s0 = vec![rec("readA/1", "AAAA")];
        let s1 = vec![rec("readB/2", "GGGG")];
        let err = zip_streams(vec![s0, s1], 0).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert!(err.to_string().contains("mismatch"), "got: {err}");
    }

    /// K=1 is a valid degenerate case: one record per template, no cross-stream
    /// check to fail.
    #[test]
    fn zip_streams_single_stream() {
        let s0 = vec![rec("read0", "AAAA"), rec("read1", "CCCC")];
        let templates = zip_streams(vec![s0], 0).unwrap();
        assert_eq!(templates.len(), 2);
        assert_eq!(templates[0].records.len(), 1);
        assert_eq!(templates[0].name, b"read0");
    }
}

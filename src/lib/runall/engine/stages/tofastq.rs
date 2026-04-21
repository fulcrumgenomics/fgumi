//! `ToFastq` stage: converts serialized BAM records to interleaved FASTQ
//! bytes.
//!
//! Parallel pool stage. Byte-for-byte compatible with `fgumi fastq`
//! standalone command output: same record filtering (excludes
//! `SECONDARY` + `SUPPLEMENTARY` by default), same `/1` `/2` suffix
//! rule (by `FIRST_SEGMENT` / `LAST_SEGMENT` flag bits, not by record
//! position), same reverse-complement for `REVERSE`-flagged reads,
//! same LF-only line endings, trailing newline on the last record.
//!
//! ## Input
//!
//! [`SerializedBatch`] — length-prefixed BAM records (as emitted by
//! `ExtractStage` / `CorrectStage` / `BamBatchedSource`). The
//! `secondary` stream is intentionally ignored (`ToFastq` only emits
//! the primary stream to FASTQ).
//!
//! ## Output
//!
//! [`SerializedFastqBatch`] — plain (uncompressed) interleaved FASTQ
//! bytes with the record count and the input `ordinal` preserved
//! verbatim.
//!
//! ## Ordering guarantees
//!
//! Record order within the batch is preserved. Cross-batch ordering
//! is the upstream `ordinal` sequence, preserved verbatim.
//!
//! ## Memory model
//!
//! 1:1 with input; output buffer preallocated to
//! `primary.data.len() * 2` (FASTQ is text with name+seq+qual
//! redundancy — ~2.5× BAM bytes in practice for 150 bp reads).
//!
//! ## Determinism
//!
//! Byte-identical per batch: conversion is a pure function of the
//! input record bytes and the three filter flags.

use anyhow::{Context, Result};
use fgumi_raw_bam as bam_fields;
use fgumi_raw_bam::RawRecordView;

use crate::runall::engine::grouping_types::iter_length_prefixed;
use crate::runall::engine::output_types::{SerializedBatch, SerializedFastqBatch};
use crate::runall::engine::stage::{Parallelism, Stage};

/// Default flags to exclude: `SECONDARY` (0x100) + `SUPPLEMENTARY` (0x800).
/// Matches `fgumi fastq`'s default.
pub const DEFAULT_EXCLUDE_FLAGS: u16 =
    bam_fields::flags::SECONDARY | bam_fields::flags::SUPPLEMENTARY;

// Phred→ASCII and base-complement tables live on `fgumi fastq` and are
// reused verbatim — BAM sequences only contain uppercase ACGTN, so the two
// tables produce identical output for every record this stage ever sees.
use crate::commands::fastq::{COMPLEMENT, QUAL_TO_ASCII};

/// Parallel stage that converts BAM records to FASTQ bytes.
#[derive(Debug, Clone, Copy)]
pub struct ToFastq {
    /// When `true`, do not append `/1` or `/2` to read names.
    no_read_suffix: bool,
    /// Records whose flags intersect this mask are skipped. Default filters
    /// out SECONDARY and SUPPLEMENTARY alignments (matches `fgumi fastq`).
    exclude_flags: u16,
    /// Records must have ALL bits in this mask set to be included. Default
    /// `0` means no additional requirement.
    require_flags: u16,
}

impl ToFastq {
    /// Construct a `ToFastq` stage with the `fgumi fastq` default filtering
    /// (exclude SECONDARY + SUPPLEMENTARY, no required flags).
    #[must_use]
    pub fn with_defaults(no_read_suffix: bool) -> Self {
        Self { no_read_suffix, exclude_flags: DEFAULT_EXCLUDE_FLAGS, require_flags: 0 }
    }

    /// Construct with explicit flag filtering.
    #[must_use]
    pub fn new(no_read_suffix: bool, exclude_flags: u16, require_flags: u16) -> Self {
        Self { no_read_suffix, exclude_flags, require_flags }
    }
}

/// Convert length-prefixed BAM record bytes into interleaved FASTQ bytes.
///
/// Iterates over `bam_data` (a concatenation of `u32 LE block_size` + payload
/// records, as emitted by `ExtractStage` / `CorrectStage` / `BamBatchedSource`)
/// and appends FASTQ-formatted text to `out`. Returns the number of records
/// written (after filtering).
///
/// # Arguments
///
/// * `bam_data` - Concatenated length-prefixed BAM record bytes.
/// * `out` - Output buffer; FASTQ bytes are appended here.
/// * `no_read_suffix` - When `true`, suppress `/1` / `/2` suffixes on read names.
/// * `exclude_flags` - Skip records whose flags intersect this mask.
/// * `require_flags` - Skip records that do not have ALL of these flag bits set
///   (pass `0` to disable).
///
/// # Panics
///
/// Does not panic in practice: the `try_into().unwrap()` for the 4-byte block-size
/// prefix is preceded by a length `ensure!` that returns an error first.
///
/// # Errors
///
/// Returns an error if `bam_data` is truncated or any record has an empty name.
pub fn bam_records_to_fastq(
    bam_data: &[u8],
    out: &mut Vec<u8>,
    no_read_suffix: bool,
    exclude_flags: u16,
    require_flags: u16,
) -> Result<u64> {
    let mut record_count: u64 = 0;

    let mut seq_buf: Vec<u8> = Vec::with_capacity(512);
    let mut qual_buf: Vec<u8> = Vec::with_capacity(512);

    for record in iter_length_prefixed(bam_data) {
        let record = record.context("ToFastq: framed-BAM parse error")?;

        let flags = RawRecordView::new(record).flags();

        // Filter: exclude any record hitting exclude_flags, or missing
        // required_flags.
        if flags & exclude_flags != 0 {
            continue;
        }
        if require_flags != 0 && (flags & require_flags) != require_flags {
            continue;
        }

        // Name (null terminator already stripped by read_name).
        let name = bam_fields::read_name(record);
        anyhow::ensure!(!name.is_empty(), "ToFastq: record with empty name");

        // Suffix: by FIRST_SEGMENT / LAST_SEGMENT flag bits, not position.
        let is_first = flags & bam_fields::flags::FIRST_SEGMENT != 0;
        let is_last = flags & bam_fields::flags::LAST_SEGMENT != 0;
        let suffix: &[u8] = if no_read_suffix {
            b""
        } else if is_first && !is_last {
            b"/1"
        } else if is_last && !is_first {
            b"/2"
        } else {
            b""
        };

        // Sequence (decoded to ASCII) + quality (Phred → ASCII+33).
        seq_buf.clear();
        seq_buf.extend(bam_fields::extract_sequence(record));

        let raw_qual = bam_fields::quality_scores_slice(record);
        qual_buf.clear();
        qual_buf.reserve(raw_qual.len());
        for &q in raw_qual {
            qual_buf.push(QUAL_TO_ASCII[q as usize]);
        }

        // Reverse-complement for REVERSE-flagged reads (matches fgumi fastq).
        if flags & bam_fields::flags::REVERSE != 0 {
            seq_buf.reverse();
            for base in &mut seq_buf {
                *base = COMPLEMENT[*base as usize];
            }
            qual_buf.reverse();
        }

        // @name{/1|/2|}\n seq\n +\n qual\n
        out.push(b'@');
        out.extend_from_slice(name);
        out.extend_from_slice(suffix);
        out.push(b'\n');
        out.extend_from_slice(&seq_buf);
        out.extend_from_slice(b"\n+\n");
        out.extend_from_slice(&qual_buf);
        out.push(b'\n');
        record_count += 1;
    }

    Ok(record_count)
}

impl Stage for ToFastq {
    type Input = SerializedBatch;
    type Output = SerializedFastqBatch;

    #[tracing::instrument(name = "tofastq", skip_all)]
    fn process(&mut self, input: Self::Input, output: &mut dyn FnMut(Self::Output)) -> Result<()> {
        let SerializedBatch { primary, secondary: _, ordinal } = input;

        // Heuristic initial capacity: ~2.5x BAM bytes, since FASTQ is text
        // with name+seq+qual redundancy. This over-allocates slightly for
        // short reads but avoids reallocation for typical 150bp records.
        let mut out: Vec<u8> = Vec::with_capacity(primary.data.len() * 2);
        let record_count = bam_records_to_fastq(
            &primary.data,
            &mut out,
            self.no_read_suffix,
            self.exclude_flags,
            self.require_flags,
        )?;

        output(SerializedFastqBatch { data: out, record_count, ordinal });
        Ok(())
    }

    fn parallelism(&self) -> Parallelism {
        Parallelism::Parallel
    }

    fn output_memory_estimate(&self, out: &Self::Output) -> usize {
        out.data.capacity()
    }

    fn name(&self) -> &'static str {
        "ToFastq"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::runall::engine::output_types::{RawBytes, SerializedBatch};

    /// Build a minimal BAM record with the given name/seq/quality/flags and
    /// no aux tags. Uses fgumi-raw-bam's test builder.
    fn make_record(name: &[u8], flag: u16, seq: &[u8], qual_phred: &[u8]) -> Vec<u8> {
        // Build the fixed-part record. We'll manually inject seq/qual since
        // testutil::make_bam_bytes zero-fills them.
        let rec =
            bam_fields::testutil::make_bam_bytes(-1, -1, flag, name, &[], seq.len(), -1, -1, &[]);
        // Pack sequence into rec at seq_offset, write quality at qual_offset.
        let mut rec = rec;
        let seq_off = bam_fields::seq_offset(&rec);
        let mut packed: Vec<u8> = Vec::new();
        bam_fields::pack_sequence_into(&mut packed, seq);
        rec[seq_off..seq_off + packed.len()].copy_from_slice(&packed);
        let qual_off = bam_fields::qual_offset(&rec);
        rec[qual_off..qual_off + qual_phred.len()].copy_from_slice(qual_phred);
        rec
    }

    fn batch(records: &[Vec<u8>], ordinal: u64) -> SerializedBatch {
        let mut data = Vec::new();
        for r in records {
            let size = u32::try_from(r.len()).unwrap();
            data.extend_from_slice(&size.to_le_bytes());
            data.extend_from_slice(r);
        }
        SerializedBatch {
            primary: RawBytes { data, record_count: records.len() as u64 },
            secondary: None,
            ordinal,
        }
    }

    #[test]
    fn test_tofastq_single_record_basic_format() {
        let mut stage = ToFastq::with_defaults(false);
        let r = make_record(b"read1", 0, b"ACGT", &[30, 31, 32, 33]);
        let input = batch(&[r], 7);

        let mut captured = None;
        stage.process(input, &mut |v| captured = Some(v)).unwrap();
        let out = captured.unwrap();

        assert_eq!(out.ordinal, 7);
        assert_eq!(out.record_count, 1);
        // @read1\nACGT\n+\n[?@AB]\n   (30+33=63=?, 31+33=64=@, 32+33=65=A, 33+33=66=B)
        assert_eq!(out.data, b"@read1\nACGT\n+\n?@AB\n");
    }

    #[test]
    fn test_tofastq_paired_first_segment_suffix() {
        let mut stage = ToFastq::with_defaults(false);
        // Flags: PAIRED (0x1) + FIRST_SEGMENT (0x40) = 0x41
        let r1 = make_record(b"readA", 0x41, b"AA", &[30, 30]);
        // Flags: PAIRED (0x1) + LAST_SEGMENT (0x80) = 0x81
        let r2 = make_record(b"readA", 0x81, b"CC", &[31, 31]);
        let input = batch(&[r1, r2], 0);

        let mut captured = None;
        stage.process(input, &mut |v| captured = Some(v)).unwrap();
        let out = captured.unwrap();

        assert_eq!(out.record_count, 2);
        let expected = b"@readA/1\nAA\n+\n??\n@readA/2\nCC\n+\n@@\n";
        assert_eq!(out.data, expected);
    }

    #[test]
    fn test_tofastq_no_suffix_disables_slash_1_slash_2() {
        let mut stage = ToFastq::with_defaults(true);
        let r1 = make_record(b"readB", 0x41, b"GG", &[20, 20]);
        let input = batch(&[r1], 0);

        let mut captured = None;
        stage.process(input, &mut |v| captured = Some(v)).unwrap();
        let out = captured.unwrap();

        assert_eq!(out.data, b"@readB\nGG\n+\n55\n");
    }

    #[test]
    fn test_tofastq_excludes_secondary_and_supplementary() {
        let mut stage = ToFastq::with_defaults(false);
        // SECONDARY = 0x100, SUPPLEMENTARY = 0x800
        let kept = make_record(b"keep", 0, b"A", &[25]);
        let secondary = make_record(b"sec", 0x100, b"A", &[25]);
        let supplementary = make_record(b"supp", 0x800, b"A", &[25]);
        let input = batch(&[kept, secondary, supplementary], 0);

        let mut captured = None;
        stage.process(input, &mut |v| captured = Some(v)).unwrap();
        let out = captured.unwrap();

        assert_eq!(out.record_count, 1);
        assert_eq!(out.data, b"@keep\nA\n+\n:\n");
    }

    #[test]
    fn test_tofastq_reverse_flag_reverse_complements() {
        let mut stage = ToFastq::with_defaults(false);
        // REVERSE = 0x10
        let r = make_record(b"rev", 0x10, b"AACG", &[30, 31, 32, 33]);
        let input = batch(&[r], 0);

        let mut captured = None;
        stage.process(input, &mut |v| captured = Some(v)).unwrap();
        let out = captured.unwrap();

        // revcomp(AACG) = CGTT; reverse(quals 30,31,32,33) = 33,32,31,30
        // ASCII: 33+33=66=B, 32+33=65=A, 31+33=64=@, 30+33=63=?
        assert_eq!(out.data, b"@rev\nCGTT\n+\nBA@?\n");
    }

    #[test]
    fn test_tofastq_empty_input_emits_empty() {
        let mut stage = ToFastq::with_defaults(false);
        let input = SerializedBatch {
            primary: RawBytes { data: Vec::new(), record_count: 0 },
            secondary: None,
            ordinal: 42,
        };

        let mut captured = None;
        stage.process(input, &mut |v| captured = Some(v)).unwrap();
        let out = captured.unwrap();

        assert_eq!(out.ordinal, 42);
        assert_eq!(out.record_count, 0);
        assert!(out.data.is_empty());
    }

    #[test]
    fn test_tofastq_preserves_ordinal() {
        let mut stage = ToFastq::with_defaults(false);
        for ord in [0u64, 1, 17, u64::MAX - 1] {
            let r = make_record(b"r", 0, b"A", &[20]);
            let input = batch(&[r], ord);
            let mut captured = None;
            stage.process(input, &mut |v| captured = Some(v)).unwrap();
            assert_eq!(captured.unwrap().ordinal, ord);
        }
    }

    #[test]
    fn test_tofastq_parallelism_and_name() {
        let stage = ToFastq::with_defaults(false);
        assert_eq!(stage.parallelism(), Parallelism::Parallel);
        assert_eq!(stage.name(), "ToFastq");
    }
}

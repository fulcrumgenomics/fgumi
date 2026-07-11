//! FASTQ record parsing utilities.
//!
//! Provides a compact `FastqRecord` type backed by a single heap allocation
//! and helpers for parsing one or more records out of a byte buffer. Used by
//! the `extract` command's BGZF/gzip FASTQ ingestion path.

use std::io;

use crate::unified_pipeline::MemoryEstimate;

/// A parsed FASTQ record stored as a single heap allocation.
///
/// The raw bytes are stored as `@name\nseq\n+\nqual\n` in `data`.
/// Boundary indices allow zero-copy slicing of name, sequence, and quality.
///
/// # Memory layout
///
/// ```text
/// data: [@][name...]['\n'][seq...]['\n'][+...]['\n'][qual...]['\n']
///         1..name_end      name_end+1..seq_end             qual_start..data.len()-1
/// ```
#[derive(Debug, Clone)]
pub struct FastqRecord {
    /// Raw bytes: `@name\nseq\n+\nqual\n`
    data: Vec<u8>,
    /// End of name slice: `data[1..name_end]` (skip leading `@`)
    name_end: u32,
    /// End of sequence slice: `data[name_end+1..seq_end]`
    seq_end: u32,
    /// Start of quality slice: `data[qual_start..data.len()-1]`
    qual_start: u32,
}

impl FastqRecord {
    /// Construct a `FastqRecord` from a raw FASTQ record slice.
    ///
    /// The slice must contain exactly one complete FASTQ record in the form
    /// `@name\nseq\n+\nqual`, where the terminating newline after the quality line
    /// is optional (a final record at EOF may omit it). CRLF (`\r\n`) line endings
    /// are accepted: the `\r` before each `\n` is a line terminator and is stripped
    /// from the name, sequence, and quality.
    ///
    /// # Errors
    ///
    /// Returns an error if the slice is empty, does not start with `@`, has a
    /// missing `+` separator line, or has mismatched sequence and quality lengths.
    pub fn from_slice(data: &[u8]) -> io::Result<Self> {
        if data.is_empty() || data[0] != b'@' {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "FASTQ record must start with @",
            ));
        }

        // Find the first three newlines to locate name, seq, and plus lines.
        // The quality ends at the fourth newline (last byte of a complete record).
        let mut newline_positions = [0usize; 3];
        let mut count = 0;
        for (i, &byte) in data.iter().enumerate() {
            if byte == b'\n' {
                if count < 3 {
                    newline_positions[count] = i;
                    count += 1;
                } else {
                    break;
                }
            }
        }

        if count < 3 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "FASTQ record must have at least 3 internal newlines",
            ));
        }

        let name_end = newline_positions[0]; // data[1..name_end] = name
        let seq_end = newline_positions[1]; // data[name_end+1..seq_end] = seq
        let plus_end = newline_positions[2]; // data[seq_end+1..plus_end] = "+"
        let qual_start = plus_end + 1;

        // Validate separator line starts with '+'
        if data[seq_end + 1] != b'+' {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "FASTQ separator line must start with +",
            ));
        }

        // Validate sequence and quality lengths match. Compare the CRLF-stripped
        // fields, not the raw spans: on a CRLF file every interior line keeps a
        // trailing '\r', but an unterminated final quality line does not, so the
        // raw spans would differ by one and spuriously reject a valid record.
        let seq_len = strip_trailing_cr(&data[name_end + 1..seq_end]).len();
        // Quality ends at the trailing '\n' (excluded) or at data.len() if no trailing newline.
        let qual_end = if data.last() == Some(&b'\n') { data.len() - 1 } else { data.len() };
        // A quality-less record (`@name\nseq\n+\n`) leaves qual_start past qual_end;
        // the slice below would panic. Reject it with a clean Err, matching the
        // sibling parser (fgumi-simd-fastq try_parse_single_record). Reachable via
        // stitch_cross_block_record / the middle-block path, which feed raw record
        // windows straight to from_slice.
        if qual_start > qual_end {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "FASTQ record is missing its quality line",
            ));
        }
        let qual_len = strip_trailing_cr(&data[qual_start..qual_end]).len();

        if seq_len != qual_len {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Sequence length ({seq_len}) != quality length ({qual_len})"),
            ));
        }

        let to_u32 = |v: usize, field: &str| {
            u32::try_from(v).map_err(|_| {
                io::Error::new(io::ErrorKind::InvalidData, format!("{field} overflows u32"))
            })
        };
        Ok(Self {
            data: data.to_vec(),
            name_end: to_u32(name_end, "name_end")?,
            seq_end: to_u32(seq_end, "seq_end")?,
            qual_start: to_u32(qual_start, "qual_start")?,
        })
    }

    /// Returns the read name bytes (excludes the leading `@`).
    #[inline]
    #[must_use]
    pub fn name(&self) -> &[u8] {
        strip_trailing_cr(&self.data[1..self.name_end as usize])
    }

    /// Returns the sequence bytes.
    #[inline]
    #[must_use]
    pub fn sequence(&self) -> &[u8] {
        strip_trailing_cr(&self.data[self.name_end as usize + 1..self.seq_end as usize])
    }

    /// Returns the quality bytes (Phred+33 encoded ASCII).
    #[inline]
    #[must_use]
    pub fn quality(&self) -> &[u8] {
        let qual_end =
            if self.data.last() == Some(&b'\n') { self.data.len() - 1 } else { self.data.len() };
        strip_trailing_cr(&self.data[self.qual_start as usize..qual_end])
    }
}

/// Strip a single trailing carriage return so CRLF (`\r\n`)-terminated FASTQ
/// files parse identically to LF-terminated ones. Record boundaries here are
/// found by scanning for `\n` only, so on a CRLF file the `\r` preceding each
/// `\n` is left as the last byte of the name, sequence, and quality fields; it
/// is a line terminator, not data, and must be dropped (matching fgbio's
/// `Source.getLines`). A stray `\r` in the quality string is byte 13, below the
/// printable-ASCII floor, so leaving it in would derail quality-encoding
/// detection. Only one trailing `\r` is removed.
#[inline]
#[must_use]
fn strip_trailing_cr(field: &[u8]) -> &[u8] {
    if let [rest @ .., b'\r'] = field { rest } else { field }
}

impl MemoryEstimate for FastqRecord {
    fn estimate_heap_size(&self) -> usize {
        self.data.capacity()
    }
}

/// Result of parsing a FASTQ record.
#[derive(Debug)]
enum FastqParseResult {
    /// Record incomplete - need more data.
    Incomplete,
    /// Parse error.
    Error(io::Error),
}

/// Parse FASTQ records from a byte buffer.
///
/// FASTQ format:
/// ```text
/// @read_name
/// SEQUENCE
/// +
/// QUALITY
/// ```
///
/// Returns parsed records and leftover bytes for incomplete records.
///
/// `at_eof` marks `data` as the complete, final input rather than one chunk of a
/// stream. It is passed through to `parse_single_fastq_record`: only when
/// `true` is a final record without a terminating newline accepted; when `false`
/// such a tail is returned as leftover so an incremental caller (e.g.
/// `FastqGrouper::add_bytes_for_stream`) can complete it with the next chunk.
///
/// # Errors
///
/// Returns an error if a record has mismatched sequence and quality lengths.
pub fn parse_fastq_records(data: &[u8], at_eof: bool) -> io::Result<(Vec<FastqRecord>, Vec<u8>)> {
    let mut records = Vec::new();
    let mut pos = 0;

    while pos < data.len() {
        // Find the start of a record (@ character at start of line)
        if data[pos] != b'@' {
            // Skip until we find @ at start of line or run out of data
            while pos < data.len() && data[pos] != b'@' {
                pos += 1;
            }
            if pos >= data.len() {
                return Ok((records, Vec::new()));
            }
        }

        // Try to parse a complete record
        match parse_single_fastq_record(&data[pos..], at_eof) {
            Ok((record, consumed)) => {
                records.push(record);
                pos += consumed;
            }
            Err(FastqParseResult::Incomplete) => {
                // Not enough data for a complete record
                return Ok((records, data[pos..].to_vec()));
            }
            Err(FastqParseResult::Error(e)) => {
                return Err(e);
            }
        }
    }

    Ok((records, Vec::new()))
}

/// Parse a single FASTQ record from the beginning of a buffer.
///
/// `at_eof` marks `data` as the complete, final input. When `true`, a record
/// whose quality line lacks a terminating newline is accepted (a final record at
/// EOF may omit it). When `false`, `data` is an incremental chunk that more bytes
/// may follow, so an unterminated final quality line is reported as
/// [`FastqParseResult::Incomplete`] and held for the next chunk rather than being
/// mistaken for a complete (shorter) record.
///
/// Returns (record, `bytes_consumed`) or an error.
fn parse_single_fastq_record(
    data: &[u8],
    at_eof: bool,
) -> Result<(FastqRecord, usize), FastqParseResult> {
    let mut pos = 0;

    // Line 1: @name
    if data.is_empty() || data[0] != b'@' {
        return Err(FastqParseResult::Error(io::Error::new(
            io::ErrorKind::InvalidData,
            "FASTQ record must start with @",
        )));
    }

    let name_end_rel = find_newline(&data[pos..]).ok_or(FastqParseResult::Incomplete)?;
    pos += name_end_rel + 1; // advance past name line + newline

    // Line 2: sequence
    if pos >= data.len() {
        return Err(FastqParseResult::Incomplete);
    }
    let seq_end_rel = find_newline(&data[pos..]).ok_or(FastqParseResult::Incomplete)?;
    // Compare CRLF-stripped lengths below: on a CRLF file every interior line keeps
    // a trailing '\r' that an unterminated final quality line would not, so raw
    // spans would differ by one and spuriously reject a valid record.
    let seq_len = strip_trailing_cr(&data[pos..pos + seq_end_rel]).len();
    pos += seq_end_rel + 1;

    // Line 3: + (separator)
    if pos >= data.len() {
        return Err(FastqParseResult::Incomplete);
    }
    if data[pos] != b'+' {
        return Err(FastqParseResult::Error(io::Error::new(
            io::ErrorKind::InvalidData,
            "FASTQ separator line must start with +",
        )));
    }
    let plus_end_rel = find_newline(&data[pos..]).ok_or(FastqParseResult::Incomplete)?;
    pos += plus_end_rel + 1;

    // Line 4: quality (may lack trailing newline at EOF)
    if pos >= data.len() {
        return Err(FastqParseResult::Incomplete);
    }
    let (qual_raw_len, advance) = if let Some(rel) = find_newline(&data[pos..]) {
        (rel, rel + 1)
    } else if at_eof {
        // Final input with no trailing newline — treat the remaining bytes as the
        // quality line. Only valid when `data` is known to be the complete input:
        // on an incremental chunk the tail may be a partial quality line that
        // continues in the next chunk, so it must be held as leftover instead.
        let remaining = data.len() - pos;
        (remaining, remaining)
    } else {
        return Err(FastqParseResult::Incomplete);
    };
    let qual_len = strip_trailing_cr(&data[pos..pos + qual_raw_len]).len();
    pos += advance;

    // Validate lengths match
    if seq_len != qual_len {
        return Err(FastqParseResult::Error(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Sequence length ({seq_len}) != quality length ({qual_len})"),
        )));
    }

    let record = FastqRecord::from_slice(&data[..pos]).map_err(FastqParseResult::Error)?;
    Ok((record, pos))
}

/// Find the position of the next newline character.
fn find_newline(data: &[u8]) -> Option<usize> {
    data.iter().position(|&b| b == b'\n')
}

/// Strip a mate-indicator suffix and CASAVA comment from a read name for comparison.
///
/// First truncates at the first ASCII space (removing CASAVA-style comments like
/// `read1 1:N:0:ATCACG`), then strips a trailing `/` followed by a single ASCII
/// digit (`0`-`9`), e.g. `read/1` -> `read`.
///
/// This mirrors fgbio's `FastqSource` read-name canonicalization (see
/// `com/fulcrumgenomics/fastq/FastqSource.scala`), which strips only a trailing
/// `/` + single digit. It deliberately does **not** strip `.`/`_`/`:` separators
/// or multi-digit runs (`read/12` is left intact): those separators are not a
/// reliable mate indicator, and treating them as one makes genuinely mismatched
/// pairs such as `read_1` (R1) vs `read_2` (R2) collapse to the same base name
/// and be silently accepted as "in sync". Used to validate that reads at the same
/// position across paired/interleaved FASTQ streams have matching names.
#[must_use]
pub fn strip_read_suffix(name: &[u8]) -> &[u8] {
    // Truncate at the first space (CASAVA comment separator).
    let name = match name.iter().position(|&b| b == b' ') {
        Some(space_pos) => &name[..space_pos],
        None => name,
    };

    // Strip a trailing '/' + single ASCII digit (fgbio FastqSource parity).
    if name.len() >= 2 {
        let last = name[name.len() - 1];
        let sep = name[name.len() - 2];
        if sep == b'/' && last.is_ascii_digit() {
            return &name[..name.len() - 2];
        }
    }
    name
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[test]
    fn test_parse_single_fastq_record() {
        let data = b"@read1\nACGT\n+\nIIII\n";
        let (record, consumed) =
            parse_single_fastq_record(data, true).expect("parse single FASTQ record");
        assert_eq!(record.name(), b"read1");
        assert_eq!(record.sequence(), b"ACGT");
        assert_eq!(record.quality(), b"IIII");
        assert_eq!(consumed, data.len());
    }

    #[test]
    fn test_parse_fastq_records_multiple() {
        let data = b"@read1\nACGT\n+\nIIII\n@read2\nTGCA\n+\nJJJJ\n";
        let (records, leftover) =
            parse_fastq_records(data, true).expect("failed to parse FASTQ records");
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name(), b"read1");
        assert_eq!(records[1].name(), b"read2");
        assert!(leftover.is_empty());
    }

    /// EXT3-06: CRLF (`\r\n`)-terminated records must expose the same
    /// name/sequence/quality as LF records — the `\r` is a line terminator, not
    /// data, and must not leak into any field (a `\r` in the quality would be
    /// byte 13, below the printable-ASCII floor, and derail encoding detection).
    #[test]
    fn from_slice_strips_crlf_terminators() {
        let rec = FastqRecord::from_slice(b"@r1\r\nACGT\r\n+\r\nIIII\r\n")
            .expect("CRLF record should parse");
        assert_eq!(rec.name(), b"r1");
        assert_eq!(rec.sequence(), b"ACGT");
        assert_eq!(rec.quality(), b"IIII");

        // Also through the streaming entry point, including a CRLF final record
        // with no trailing newline.
        let (records, leftover) =
            parse_fastq_records(b"@r1\r\nACGT\r\n+\r\nIIII\r\n@r2\r\nTGCA\r\n+\r\nJJJJ", true)
                .expect("CRLF stream should parse");
        assert_eq!(records.len(), 2);
        assert_eq!(records[1].name(), b"r2");
        assert_eq!(records[1].sequence(), b"TGCA");
        assert_eq!(records[1].quality(), b"JJJJ");
        assert!(leftover.is_empty(), "the unterminated CRLF final record must be consumed");
    }

    /// A quality-less record (`@name\nseq\n+\n`, no quality line) yields
    /// `qual_start > qual_end`. `from_slice` must reject it with a clean `Err`
    /// rather than panicking on the `data[qual_start..qual_end]` slice. This is
    /// reachable in production: `stitch_cross_block_record` and the middle-block
    /// path hand raw record windows straight to `from_slice`, so a quality-less
    /// record straddling a BGZF block boundary reconstructs to exactly this
    /// shape. The sibling parser (`fgumi-simd-fastq::try_parse_single_record`)
    /// already guards the identical case; both parsers must agree (reject, not
    /// panic).
    #[rstest]
    #[case::lf(b"@r\nA\n+\n")]
    #[case::crlf(b"@r\r\nA\r\n+\r\n")]
    #[case::longer_lf(b"@read1\nACGT\n+\n")]
    #[case::plus_only_no_trailing_newline(b"@r\nA\n+")]
    fn from_slice_rejects_quality_less_record_without_panicking(#[case] data: &[u8]) {
        let result = FastqRecord::from_slice(data);
        assert!(
            result.is_err(),
            "quality-less record {data:?} must return Err, not panic or succeed"
        );
    }

    #[test]
    fn test_parse_fastq_incomplete_record() {
        let data = b"@read1\nACGT\n+\n";
        let (records, leftover) =
            parse_fastq_records(data, true).expect("failed to parse FASTQ records");
        assert!(records.is_empty());
        assert_eq!(leftover, data);
    }

    #[test]
    fn test_parse_fastq_eof_without_trailing_newline() {
        // Quality line has no trailing newline — should still parse at EOF.
        let data = b"@read1\nACGT\n+\nIIII";
        let (records, leftover) =
            parse_fastq_records(data, true).expect("failed to parse FASTQ records");
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].name(), b"read1");
        assert_eq!(records[0].sequence(), b"ACGT");
        assert_eq!(records[0].quality(), b"IIII");
        assert!(leftover.is_empty());
    }

    /// Mirror of `test_parse_fastq_eof_without_trailing_newline` for the
    /// incremental path: on a non-final chunk (`at_eof = false`) an unterminated
    /// final record must NOT be consumed — its bytes are returned as leftover so
    /// the caller can complete it with the next chunk. Otherwise a chunk boundary
    /// landing mid-quality would be mistaken for a complete (shorter) record.
    #[test]
    fn test_parse_fastq_not_eof_holds_unterminated_record_as_leftover() {
        let data = b"@read1\nACGT\n+\nIIII";
        let (records, leftover) =
            parse_fastq_records(data, false).expect("non-final chunk must not error");
        assert!(records.is_empty(), "an unterminated final record must not be yielded mid-stream");
        assert_eq!(leftover, data, "the whole record must be held as leftover for the next chunk");
    }

    #[test]
    fn test_parse_fastq_eof_no_newline_seq_qual_mismatch() {
        // Quality shorter than sequence at EOF — should error.
        let data = b"@read1\nACGT\n+\nIII";
        let result = parse_fastq_records(data, true);
        assert!(
            result.is_err() || {
                let (recs, leftover) = result.unwrap();
                recs.is_empty() && !leftover.is_empty()
            }
        );
    }

    #[test]
    fn test_parse_fastq_with_leftover() {
        let data = b"@read1\nACGT\n+\nIIII\n@read2\nTG";
        let (records, leftover) =
            parse_fastq_records(data, true).expect("failed to parse FASTQ records");
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].name(), b"read1");
        assert_eq!(leftover, b"@read2\nTG");
    }

    #[rstest]
    // A trailing '/' + single digit is stripped (any digit, fgbio parity).
    #[case::slash_one(b"read/1", b"read")]
    #[case::slash_two(b"read/2", b"read")]
    #[case::slash_three(b"read/3", b"read")]
    #[case::slash_zero(b"read/0", b"read")]
    // Non-'/' separators are preserved: they are not a reliable mate indicator,
    // and stripping them collapses genuinely different reads (see the collision
    // test below). This narrows the previous over-lenient behavior to match
    // fgbio's FastqSource.
    #[case::dot_one(b"read.1", b"read.1")]
    #[case::dot_two(b"read.2", b"read.2")]
    #[case::underscore_one(b"read_1", b"read_1")]
    #[case::underscore_two(b"read_2", b"read_2")]
    #[case::colon_one(b"read:1", b"read:1")]
    #[case::colon_two(b"read:2", b"read:2")]
    // Multi-digit runs are preserved (only a single trailing digit is stripped).
    #[case::slash_multidigit(b"read/12", b"read/12")]
    // Non-digit after '/' is preserved.
    #[case::slash_nondigit(b"read/x", b"read/x")]
    // A trailing space comment is stripped first, then the '/'+digit suffix.
    #[case::comment_and_slash(b"read/1 1:N:0:ATCACG", b"read")]
    // A comment is stripped even when there is no mate suffix.
    #[case::comment_only(b"read 1:N:0:ATCACG", b"read")]
    // No suffix at all.
    #[case::no_suffix(b"read", b"read")]
    #[case::single_char(b"a", b"a")]
    #[case::empty(b"", b"")]
    fn test_strip_read_suffix(#[case] input: &[u8], #[case] expected: &[u8]) {
        assert_eq!(strip_read_suffix(input), expected);
    }

    #[test]
    fn test_strip_read_suffix_rejects_mismatched_underscore_pair() {
        // Behavior tradeoff / regression guard: `read_1` (R1) and `read_2` (R2)
        // are genuinely different reads. The previous over-lenient rule stripped
        // the `_1`/`_2` suffix, collapsing both to `read` so a mismatched pair of
        // streams was silently accepted as "in sync". The fgbio-parity rule now
        // preserves them, so the two no longer collide and sync-validation
        // correctly rejects the pair. This also means non-standard paired FASTQs
        // that legitimately use `_1`/`_2` (or `.1`/`.2`/`:1`/`:2`) suffixes are
        // now rejected — matching fgbio, which rejects them too.
        assert_eq!(strip_read_suffix(b"read_1"), b"read_1");
        assert_eq!(strip_read_suffix(b"read_2"), b"read_2");
        assert_ne!(strip_read_suffix(b"read_1"), strip_read_suffix(b"read_2"));
    }
}

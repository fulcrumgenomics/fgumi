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
    /// The slice must contain exactly one complete FASTQ record in the form:
    /// `@name\nseq\n+\nqual\n`
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

        // Validate sequence and quality lengths match.
        let seq_len = seq_end - (name_end + 1);
        // Quality ends at the trailing '\n' (excluded) or at data.len() if no trailing newline.
        let qual_end = if data.last() == Some(&b'\n') { data.len() - 1 } else { data.len() };
        let qual_len = qual_end - qual_start;

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
        &self.data[1..self.name_end as usize]
    }

    /// Returns the sequence bytes.
    #[inline]
    #[must_use]
    pub fn sequence(&self) -> &[u8] {
        &self.data[self.name_end as usize + 1..self.seq_end as usize]
    }

    /// Returns the quality bytes (Phred+33 encoded ASCII).
    #[inline]
    #[must_use]
    pub fn quality(&self) -> &[u8] {
        let qual_end =
            if self.data.last() == Some(&b'\n') { self.data.len() - 1 } else { self.data.len() };
        &self.data[self.qual_start as usize..qual_end]
    }
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
/// # Errors
///
/// Returns an error if a record has mismatched sequence and quality lengths.
pub fn parse_fastq_records(data: &[u8]) -> io::Result<(Vec<FastqRecord>, Vec<u8>)> {
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
        match parse_single_fastq_record(&data[pos..]) {
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
/// Returns (record, `bytes_consumed`) or an error.
fn parse_single_fastq_record(data: &[u8]) -> Result<(FastqRecord, usize), FastqParseResult> {
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
    let seq_len = seq_end_rel;
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
    let (qual_len, advance) = if let Some(rel) = find_newline(&data[pos..]) {
        (rel, rel + 1)
    } else {
        // EOF without trailing newline — treat remaining bytes as the quality line.
        let remaining = data.len() - pos;
        (remaining, remaining)
    };
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

/// Strip mate-indicator suffix and CASAVA comment from a read name for comparison.
///
/// First truncates at the first ASCII space (removing CASAVA-style comments like
/// `read1 1:N:0:ATCACG`), then strips common pair suffixes: `/1`, `/2`, `.1`,
/// `.2`, `_1`, `_2`, `:1`, `:2`.
#[must_use]
pub fn strip_read_suffix(name: &[u8]) -> &[u8] {
    // Truncate at the first space (CASAVA comment separator).
    let name = match name.iter().position(|&b| b == b' ') {
        Some(space_pos) => &name[..space_pos],
        None => name,
    };

    // Strip common pair suffixes: separator + digit where separator is
    // one of '/', '.', '_', ':' and digit is '1' or '2'.
    if name.len() >= 2 {
        let last = name[name.len() - 1];
        let sep = name[name.len() - 2];
        if (last == b'1' || last == b'2')
            && (sep == b'/' || sep == b'.' || sep == b'_' || sep == b':')
        {
            return &name[..name.len() - 2];
        }
    }
    name
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_single_fastq_record() {
        let data = b"@read1\nACGT\n+\nIIII\n";
        let (record, consumed) =
            parse_single_fastq_record(data).expect("parse single FASTQ record");
        assert_eq!(record.name(), b"read1");
        assert_eq!(record.sequence(), b"ACGT");
        assert_eq!(record.quality(), b"IIII");
        assert_eq!(consumed, data.len());
    }

    #[test]
    fn test_parse_fastq_records_multiple() {
        let data = b"@read1\nACGT\n+\nIIII\n@read2\nTGCA\n+\nJJJJ\n";
        let (records, leftover) = parse_fastq_records(data).expect("failed to parse FASTQ records");
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name(), b"read1");
        assert_eq!(records[1].name(), b"read2");
        assert!(leftover.is_empty());
    }

    #[test]
    fn test_parse_fastq_incomplete_record() {
        let data = b"@read1\nACGT\n+\n";
        let (records, leftover) = parse_fastq_records(data).expect("failed to parse FASTQ records");
        assert!(records.is_empty());
        assert_eq!(leftover, data);
    }

    #[test]
    fn test_parse_fastq_eof_without_trailing_newline() {
        // Quality line has no trailing newline — should still parse.
        let data = b"@read1\nACGT\n+\nIIII";
        let (records, leftover) = parse_fastq_records(data).expect("failed to parse FASTQ records");
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].name(), b"read1");
        assert_eq!(records[0].sequence(), b"ACGT");
        assert_eq!(records[0].quality(), b"IIII");
        assert!(leftover.is_empty());
    }

    #[test]
    fn test_parse_fastq_eof_no_newline_seq_qual_mismatch() {
        // Quality shorter than sequence at EOF — should error.
        let data = b"@read1\nACGT\n+\nIII";
        let result = parse_fastq_records(data);
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
        let (records, leftover) = parse_fastq_records(data).expect("failed to parse FASTQ records");
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].name(), b"read1");
        assert_eq!(leftover, b"@read2\nTG");
    }

    #[test]
    fn test_strip_read_suffix() {
        // Slash separators (original behavior).
        assert_eq!(strip_read_suffix(b"read1/1"), b"read1");
        assert_eq!(strip_read_suffix(b"read1/2"), b"read1");
        // Dot, underscore, colon separators.
        assert_eq!(strip_read_suffix(b"read1.1"), b"read1");
        assert_eq!(strip_read_suffix(b"read1.2"), b"read1");
        assert_eq!(strip_read_suffix(b"read1_1"), b"read1");
        assert_eq!(strip_read_suffix(b"read1_2"), b"read1");
        assert_eq!(strip_read_suffix(b"read1:1"), b"read1");
        assert_eq!(strip_read_suffix(b"read1:2"), b"read1");
        // CASAVA-style headers with space-separated comments.
        assert_eq!(strip_read_suffix(b"read1/1 1:N:0:ATCACG"), b"read1");
        // No pair suffix, but CASAVA comment is still stripped.
        assert_eq!(strip_read_suffix(b"read1 1:N:0:ATCACG"), b"read1");
        // No suffix.
        assert_eq!(strip_read_suffix(b"read1"), b"read1");
        assert_eq!(strip_read_suffix(b"a"), b"a");
        assert_eq!(strip_read_suffix(b""), b"" as &[u8]);
    }
}

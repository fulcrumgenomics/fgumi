//! FASTQ parser: uses SIMD-produced newline bitmasks to find record boundaries.

use crate::lexer;

/// A borrowed FASTQ record with zero-copy slices into the input buffer.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastqRecord<'a> {
    /// Read name (without leading `@` or trailing newline).
    pub name: &'a [u8],
    /// Sequence bases.
    pub sequence: &'a [u8],
    /// Quality scores (Phred-encoded ASCII).
    pub quality: &'a [u8],
}

/// Find FASTQ record boundary offsets in a byte buffer.
///
/// Returns a `Vec` of byte offsets marking the start of each complete record.
/// The first offset is always 0 (if there are any records). The last offset
/// points one byte past the final complete record.
///
/// Incomplete trailing records (fewer than 4 newlines) are excluded; the caller
/// should treat `data[last_offset..]` as leftover bytes.
///
/// # Algorithm
///
/// Processes 64-byte blocks through the SIMD lexer to produce newline bitmasks.
/// Walks each bitmask with `trailing_zeros()` to locate newlines. Every 4th
/// newline marks a record boundary.
///
/// # Example
///
/// ```
/// use fgumi_simd_fastq::find_record_offsets;
///
/// let data = b"@r1\nACGT\n+\nIIII\n@r2\nTT\n+\nJJ\n";
/// let offsets = find_record_offsets(data);
/// assert_eq!(offsets, vec![0, 16, 28]);
/// ```
#[must_use]
#[allow(clippy::missing_panics_doc)]
pub fn find_record_offsets(data: &[u8]) -> Vec<usize> {
    if data.is_empty() {
        return vec![0];
    }

    let mut offsets = Vec::with_capacity(data.len() / 300 + 2);
    offsets.push(0);

    let mut newline_count: u32 = 0;
    let num_full_blocks = data.len() / 64;
    let remainder = data.len() % 64;

    // Process full 64-byte blocks
    for block_idx in 0..num_full_blocks {
        let block_start = block_idx * 64;
        // Indexing is safe: block_start + 64 <= num_full_blocks * 64 <= data.len()
        let block: &[u8; 64] =
            data[block_start..block_start + 64].try_into().expect("slice is exactly 64 bytes");
        let newlines = lexer::lex_block(block);
        process_bitmask(newlines, block_start, &mut newline_count, &mut offsets);
    }

    // Process the remaining bytes (< 64) with zero-padding
    if remainder > 0 {
        let block_start = num_full_blocks * 64;
        let mut padded = [0u8; 64];
        padded[..remainder].copy_from_slice(&data[block_start..]);
        let newlines = lexer::lex_block(&padded);
        // Mask out padding positions — only consider bits 0..remainder
        let valid_mask = if remainder < 64 { (1u64 << remainder) - 1 } else { u64::MAX };
        let masked_newlines = newlines & valid_mask;
        process_bitmask(masked_newlines, block_start, &mut newline_count, &mut offsets);
    }

    offsets
}

/// Process a newline bitmask from one 64-byte block, updating the record offsets.
///
/// Every 4th newline marks the end of a FASTQ record (the byte after the newline
/// is the start of the next record).
#[inline]
fn process_bitmask(
    mut newlines: u64,
    block_start: usize,
    newline_count: &mut u32,
    offsets: &mut Vec<usize>,
) {
    while newlines != 0 {
        let bit_pos = newlines.trailing_zeros() as usize;
        *newline_count += 1;

        if (*newline_count).is_multiple_of(4) {
            // This newline ends a FASTQ record. The next record starts at bit_pos + 1.
            offsets.push(block_start + bit_pos + 1);
        }

        // Clear the lowest set bit
        newlines &= newlines - 1;
    }
}

/// Iterator over zero-copy FASTQ records parsed from a byte buffer.
///
/// Uses [`find_record_offsets`] to locate record boundaries, then slices
/// each record into name, sequence, and quality fields.
///
/// # Example
///
/// ```
/// use fgumi_simd_fastq::parse_records;
///
/// let data = b"@read1\nACGT\n+\nIIII\n";
/// let records: Vec<_> = parse_records(data).collect();
/// assert_eq!(records[0].name, b"read1");
/// assert_eq!(records[0].sequence, b"ACGT");
/// assert_eq!(records[0].quality, b"IIII");
/// ```
pub fn parse_records(data: &[u8]) -> impl Iterator<Item = FastqRecord<'_>> {
    let offsets = find_record_offsets(data);
    RecordIter { data, offsets, idx: 0 }
}

/// Iterator over FASTQ records in a byte buffer.
struct RecordIter<'a> {
    data: &'a [u8],
    offsets: Vec<usize>,
    idx: usize,
}

impl<'a> Iterator for RecordIter<'a> {
    type Item = FastqRecord<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.idx + 1 >= self.offsets.len() {
            return None;
        }
        let start = self.offsets[self.idx];
        let end = self.offsets[self.idx + 1];
        self.idx += 1;
        Some(parse_single_record(&self.data[start..end]))
    }
}

/// Parse a single FASTQ record from a byte slice.
///
/// Expects format: `@name\nsequence\n+\nquality\n`
///
/// # Panics
///
/// Panics if the record does not contain exactly 4 newline-delimited lines
/// or does not start with `@`.
pub(crate) fn parse_single_record(record: &[u8]) -> FastqRecord<'_> {
    assert!(!record.is_empty() && record[0] == b'@', "FASTQ record must start with @");

    // Find the 3 internal newlines (the 4th newline is the last byte)
    let mut newline_positions = [0usize; 3];
    let mut count = 0;

    // We only need to find the first 3 newlines; the 4th is at end-1
    for (i, &byte) in record.iter().enumerate() {
        if byte == b'\n' {
            if count < 3 {
                newline_positions[count] = i;
                count += 1;
            } else {
                break;
            }
        }
    }

    assert_eq!(count, 3, "FASTQ record must have at least 3 internal newlines");

    let name_end = newline_positions[0];
    let seq_end = newline_positions[1];
    let plus_end = newline_positions[2];

    // name: skip '@', up to first newline
    let name = &record[1..name_end];
    // sequence: after first newline, up to second newline
    let sequence = &record[name_end + 1..seq_end];
    // quality: after third newline, up to end (excluding trailing newline if present)
    let qual_start = plus_end + 1;
    let qual_end = if record.last() == Some(&b'\n') { record.len() - 1 } else { record.len() };
    let quality = &record[qual_start..qual_end];

    FastqRecord { name, sequence, quality }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_input() {
        let offsets = find_record_offsets(b"");
        assert_eq!(offsets, vec![0]);
    }

    #[test]
    fn test_single_record() {
        let data = b"@r1\nACGT\n+\nIIII\n";
        let offsets = find_record_offsets(data);
        assert_eq!(offsets, vec![0, 16]);
    }

    #[test]
    fn test_two_records() {
        let data = b"@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nJJJJ\n";
        let offsets = find_record_offsets(data);
        assert_eq!(offsets, vec![0, 16, 32]);
    }

    #[test]
    fn test_incomplete_trailing_record() {
        let data = b"@r1\nACGT\n+\nIIII\n@r2\nTT";
        let offsets = find_record_offsets(data);
        assert_eq!(offsets, vec![0, 16]);
    }

    #[test]
    fn test_single_base_reads() {
        // @r\nA\n+\nI\n = 9 bytes: @(0) r(1) \n(2) A(3) \n(4) +(5) \n(6) I(7) \n(8)
        let data = b"@r\nA\n+\nI\n";
        assert_eq!(data.len(), 9);
        let offsets = find_record_offsets(data);
        assert_eq!(offsets, vec![0, 9]);
    }

    #[test]
    fn test_parse_single_record() {
        let data = b"@read1\nACGT\n+\nIIII\n";
        let records: Vec<_> = parse_records(data).collect();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].name, b"read1");
        assert_eq!(records[0].sequence, b"ACGT");
        assert_eq!(records[0].quality, b"IIII");
    }

    #[test]
    fn test_parse_multiple_records() {
        let data = b"@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nJJJJ\n";
        let records: Vec<_> = parse_records(data).collect();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name, b"r1");
        assert_eq!(records[0].sequence, b"ACGT");
        assert_eq!(records[0].quality, b"IIII");
        assert_eq!(records[1].name, b"r2");
        assert_eq!(records[1].sequence, b"TTTT");
        assert_eq!(records[1].quality, b"JJJJ");
    }

    #[test]
    fn test_record_spanning_block_boundary() {
        // Create a record that starts in one 64-byte block and ends in the next
        // "@" + 60 chars name + "\n" = 62 bytes for header line
        // "ACGT\n+\nIIII\n" = 13 bytes for remaining lines
        // Total: 75 bytes — spans two 64-byte blocks
        let name = "X".repeat(60);
        let data = format!("@{name}\nACGT\n+\nIIII\n");
        let offsets = find_record_offsets(data.as_bytes());
        assert_eq!(offsets, vec![0, data.len()]);

        let records: Vec<_> = parse_records(data.as_bytes()).collect();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].name, name.as_bytes());
        assert_eq!(records[0].sequence, b"ACGT");
        assert_eq!(records[0].quality, b"IIII");
    }

    #[test]
    fn test_long_sequence_spanning_multiple_blocks() {
        // 200-base sequence spans 4 blocks
        let seq = "A".repeat(200);
        let qual = "I".repeat(200);
        let data = format!("@r1\n{seq}\n+\n{qual}\n");
        let offsets = find_record_offsets(data.as_bytes());
        assert_eq!(offsets, vec![0, data.len()]);

        let records: Vec<_> = parse_records(data.as_bytes()).collect();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].sequence.len(), 200);
    }

    #[test]
    fn test_n_bases_and_mixed_case() {
        let data = b"@r1\nAcGtNn\n+\nIIIIII\n";
        let records: Vec<_> = parse_records(data).collect();
        assert_eq!(records[0].sequence, b"AcGtNn");
    }

    #[test]
    fn test_many_small_records() {
        // 10 minimal records
        let mut data = Vec::new();
        for i in 0..10 {
            data.extend_from_slice(format!("@r{i}\nA\n+\nI\n").as_bytes());
        }
        let offsets = find_record_offsets(&data);
        assert_eq!(offsets.len(), 11); // 10 records + initial 0
        let records: Vec<_> = parse_records(&data).collect();
        assert_eq!(records.len(), 10);
        for (i, rec) in records.iter().enumerate() {
            assert_eq!(rec.name, format!("r{i}").as_bytes());
        }
    }

    #[test]
    fn test_exactly_64_bytes() {
        // Craft a record that is exactly 64 bytes
        // @name(55 chars)\nA\n+\nI\n = 55 + 1 + 1 + 1 + 1 + 1 + 1 + 1 + 1 + 1 = hmm
        // @(1) + name(53) + \n(1) + A(1) + \n(1) + +(1) + \n(1) + I(1) + \n(1) = 61
        // Need 64: @(1) + name(56) + \n(1) + A(1) + \n(1) + +(1) + \n(1) + I(1) + \n(1) = 64
        let name = "X".repeat(56);
        let data = format!("@{name}\nA\n+\nI\n");
        assert_eq!(data.len(), 64);
        let offsets = find_record_offsets(data.as_bytes());
        assert_eq!(offsets, vec![0, 64]);
    }

    #[test]
    fn test_no_complete_records() {
        let data = b"@r1\nACGT\n+\n";
        let offsets = find_record_offsets(data);
        assert_eq!(offsets, vec![0]);
    }
}

//! Buffered FASTQ reader using SIMD-accelerated parsing.
//!
//! [`SimdFastqReader`] wraps a `BufRead` and yields owned FASTQ records,
//! serving as a drop-in replacement for `seq_io::fastq::Reader`.

use std::io::{self, BufRead};

use crate::parser::{self, parse_single_record};

/// Default internal buffer size (1 MiB), matching `seq_io`'s default.
const DEFAULT_BUFFER_SIZE: usize = 1 << 20;

/// An owned FASTQ record with heap-allocated name, sequence, and quality.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OwnedFastqRecord {
    /// Read name (without leading `@`).
    pub name: Vec<u8>,
    /// Sequence bases.
    pub sequence: Vec<u8>,
    /// Quality scores (Phred-encoded ASCII).
    pub quality: Vec<u8>,
}

/// Buffered FASTQ reader that uses SIMD-accelerated record boundary detection.
///
/// Reads chunks from the underlying `BufRead`, finds record boundaries with
/// [`find_record_offsets`](crate::find_record_offsets), and yields owned records.
///
/// # Example
///
/// ```
/// use fgumi_simd_fastq::SimdFastqReader;
/// use std::io::Cursor;
///
/// let data = b"@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nJJJJ\n";
/// let mut reader = SimdFastqReader::new(Cursor::new(&data[..]));
///
/// let rec = reader.next().unwrap().unwrap();
/// assert_eq!(rec.name, b"r1");
/// assert_eq!(rec.sequence, b"ACGT");
/// ```
pub struct SimdFastqReader<R: BufRead> {
    inner: R,
    /// Internal buffer holding data read from the source.
    buffer: Vec<u8>,
    /// Pre-computed record boundary offsets within `buffer[..valid]`.
    offsets: Vec<usize>,
    /// Index into `offsets` for the next record to yield.
    next_record_idx: usize,
    /// Number of valid bytes in `buffer`.
    valid: usize,
    /// True when the underlying reader has returned 0 bytes.
    at_eof: bool,
}

impl<R: BufRead> SimdFastqReader<R> {
    /// Create a new reader with the default buffer size (1 MiB).
    pub fn new(inner: R) -> Self {
        Self::with_capacity(inner, DEFAULT_BUFFER_SIZE)
    }

    /// Create a new reader with a custom buffer capacity.
    pub fn with_capacity(inner: R, capacity: usize) -> Self {
        Self {
            inner,
            buffer: Vec::with_capacity(capacity),
            offsets: Vec::new(),
            next_record_idx: 0,
            valid: 0,
            at_eof: false,
        }
    }

    /// Fill the internal buffer, preserving any leftover bytes from incomplete records.
    ///
    /// Returns `true` if new data was read (or there are still records to yield).
    fn fill_buffer(&mut self) -> io::Result<bool> {
        // Determine leftover: bytes from the last complete record offset to end of valid data.
        let leftover_start = if self.offsets.is_empty() {
            0
        } else {
            // The last offset in `offsets` is the start of the first incomplete record
            // (or the end of the last complete record, which is the same thing).
            self.offsets.last().copied().unwrap_or(0)
        };

        // Move leftover bytes to the front of the buffer
        if leftover_start > 0 && leftover_start < self.valid {
            self.buffer.copy_within(leftover_start..self.valid, 0);
            self.valid -= leftover_start;
        } else if leftover_start >= self.valid {
            self.valid = 0;
        }

        // If the buffer is full of leftover (no complete records found), grow it
        // so we can read more data and find the end of the current record.
        if self.valid >= self.buffer.capacity() {
            self.buffer.reserve(self.buffer.capacity().max(4096));
        }
        self.buffer.resize(self.buffer.capacity(), 0);

        // Read new data into the buffer after the leftover
        let mut total_read = 0;
        while self.valid + total_read < self.buffer.len() {
            let buf = &mut self.buffer[self.valid + total_read..];
            if buf.is_empty() {
                break;
            }
            match self.inner.read(buf) {
                Ok(0) => {
                    self.at_eof = true;
                    break;
                }
                Ok(n) => total_read += n,
                Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
                Err(e) => return Err(e),
            }
        }

        self.valid += total_read;
        self.buffer.truncate(self.valid);

        // Find record boundaries in the buffer
        self.offsets = parser::find_record_offsets(&self.buffer[..self.valid]);
        self.next_record_idx = 0;

        // We have data if there are any complete records
        Ok(self.offsets.len() > 1 || (!self.at_eof && self.valid > 0))
    }
}

impl<R: BufRead> Iterator for SimdFastqReader<R> {
    type Item = io::Result<OwnedFastqRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if self.next_record_idx + 1 < self.offsets.len() {
                let start = self.offsets[self.next_record_idx];
                let end = self.offsets[self.next_record_idx + 1];
                self.next_record_idx += 1;

                let borrowed = parse_single_record(&self.buffer[start..end]);
                return Some(Ok(OwnedFastqRecord {
                    name: borrowed.name.to_vec(),
                    sequence: borrowed.sequence.to_vec(),
                    quality: borrowed.quality.to_vec(),
                }));
            }

            if self.at_eof {
                // Check for leftover bytes that form an incomplete record
                let leftover_start = self.offsets.last().copied().unwrap_or(0);
                if leftover_start < self.valid {
                    return Some(Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "Truncated FASTQ record at EOF ({} leftover bytes)",
                            self.valid - leftover_start
                        ),
                    )));
                }
                return None;
            }

            match self.fill_buffer() {
                Ok(true) => {}
                Ok(false) => return None,
                Err(e) => return Some(Err(e)),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_reader_single_record() {
        let data = b"@r1\nACGT\n+\nIIII\n";
        let mut reader = SimdFastqReader::new(Cursor::new(&data[..]));

        let rec = reader
            .next()
            .expect("reader should yield a record")
            .expect("record should parse successfully");
        assert_eq!(rec.name, b"r1");
        assert_eq!(rec.sequence, b"ACGT");
        assert_eq!(rec.quality, b"IIII");

        assert!(reader.next().is_none());
    }

    #[test]
    fn test_reader_multiple_records() {
        let data = b"@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nJJJJ\n";
        let mut reader = SimdFastqReader::new(Cursor::new(&data[..]));

        let rec1 = reader
            .next()
            .expect("reader should yield a record")
            .expect("record should parse successfully");
        assert_eq!(rec1.name, b"r1");

        let rec2 = reader
            .next()
            .expect("reader should yield a record")
            .expect("record should parse successfully");
        assert_eq!(rec2.name, b"r2");

        assert!(reader.next().is_none());
    }

    #[test]
    fn test_reader_tiny_buffer() {
        // Use a very small buffer to force multiple refills
        let data = b"@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nJJJJ\n";
        let mut reader = SimdFastqReader::with_capacity(Cursor::new(&data[..]), 20);

        let rec1 = reader
            .next()
            .expect("reader should yield a record")
            .expect("record should parse successfully");
        assert_eq!(rec1.name, b"r1");
        assert_eq!(rec1.sequence, b"ACGT");

        let rec2 = reader
            .next()
            .expect("reader should yield a record")
            .expect("record should parse successfully");
        assert_eq!(rec2.name, b"r2");
        assert_eq!(rec2.sequence, b"TTTT");

        assert!(reader.next().is_none());
    }

    #[test]
    fn test_reader_empty_input() {
        let data = b"";
        let mut reader = SimdFastqReader::new(Cursor::new(&data[..]));
        assert!(reader.next().is_none());
    }

    #[test]
    fn test_reader_long_records() {
        let seq = "A".repeat(500);
        let qual = "I".repeat(500);
        let data = format!("@longread\n{seq}\n+\n{qual}\n");
        let mut reader = SimdFastqReader::with_capacity(Cursor::new(data.as_bytes()), 256);

        let rec = reader
            .next()
            .expect("reader should yield a record")
            .expect("record should parse successfully");
        assert_eq!(rec.name, b"longread");
        assert_eq!(rec.sequence.len(), 500);
        assert_eq!(rec.quality.len(), 500);

        assert!(reader.next().is_none());
    }
}

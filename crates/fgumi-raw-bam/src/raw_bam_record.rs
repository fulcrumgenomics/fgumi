//! Raw BAM record reader that bypasses noodles' Record type.
//!
//! This module provides [`RawRecord`], a zero-overhead wrapper around raw BAM bytes
//! that enables high-performance sorting without full record decoding.
//!
//! # Why Not Use `noodles::bam::Record`?
//!
//! The standard `noodles::bam::Record` type wraps raw bytes but doesn't expose them
//! via `AsRef<[u8]>`. This module provides direct access to raw bytes for:
//! - Zero-copy field extraction using fixed byte offsets
//! - Direct writes to output without re-encoding
//! - 3-4x lower memory usage than decoded `RecordBuf`
//!
//! # BAM Record Binary Layout
//!
//! ```text
//! Offset  Size  Field
//! ------  ----  -----
//! 0-3     4     refID (i32) - reference sequence ID
//! 4-7     4     pos (i32) - 0-based leftmost position
//! 8       1     l_read_name (u8) - length including null
//! 9       1     mapq (u8) - mapping quality
//! 10-11   2     bin (u16) - BAM bin
//! 12-13   2     n_cigar_op (u16) - CIGAR operation count
//! 14-15   2     flag (u16) - bitwise flags
//! 16-19   4     l_seq (u32) - sequence length
//! 20-23   4     next_refID (i32) - mate reference ID
//! 24-27   4     next_pos (i32) - mate position
//! 28-31   4     tlen (i32) - template length
//! 32+     var   read_name, CIGAR, sequence, quality, aux data
//! ```
//!
//! This module will be removed once noodles exposes raw bytes from Record.

use std::io::{self, Read};

/// A raw BAM record stored as bytes.
///
/// This is a zero-overhead wrapper that provides `AsRef<[u8]>` access to the
/// underlying BAM record bytes, enabling high-performance field extraction
/// and direct output writes.
#[derive(Clone, Default, Eq, PartialEq)]
pub struct RawRecord(Vec<u8>);

impl RawRecord {
    /// Creates a new empty raw record.
    #[inline]
    #[must_use]
    pub fn new() -> Self {
        Self(Vec::new())
    }

    /// Creates a raw record with the given capacity.
    #[inline]
    #[must_use]
    pub fn with_capacity(capacity: usize) -> Self {
        Self(Vec::with_capacity(capacity))
    }

    /// Returns the length of the record in bytes.
    #[inline]
    #[must_use]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Returns true if the record is empty.
    #[inline]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Clears the record, removing all bytes.
    #[inline]
    pub fn clear(&mut self) {
        self.0.clear();
    }

    /// Returns the inner bytes, consuming the record.
    #[inline]
    #[must_use]
    pub fn into_inner(self) -> Vec<u8> {
        self.0
    }
}

impl AsRef<[u8]> for RawRecord {
    #[inline]
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

impl std::ops::Deref for RawRecord {
    type Target = [u8];

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl From<Vec<u8>> for RawRecord {
    #[inline]
    fn from(buf: Vec<u8>) -> Self {
        Self(buf)
    }
}

/// Reads a single raw BAM record from the given reader.
///
/// This function reads:
/// 1. The 4-byte `block_size` prefix (little-endian u32)
/// 2. The record data (`block_size` bytes)
///
/// Returns the number of bytes read (excluding the 4-byte prefix), or 0 at EOF.
///
/// # Errors
///
/// Returns an error if:
/// - The reader returns an error
/// - EOF is reached in the middle of a record
pub fn read_raw_record<R>(reader: &mut R, record: &mut RawRecord) -> io::Result<usize>
where
    R: Read,
{
    let block_size = match read_block_size(reader)? {
        0 => return Ok(0),
        n => n,
    };

    record.0.resize(block_size, 0);
    reader.read_exact(&mut record.0)?;

    Ok(block_size)
}

/// Reads the 4-byte block size prefix.
///
/// Returns 0 at EOF (no bytes available).
fn read_block_size<R>(reader: &mut R) -> io::Result<usize>
where
    R: Read,
{
    let mut buf = [0u8; 4];

    // Try to read the first byte to detect EOF
    match reader.read(&mut buf[..1]) {
        Ok(0) => return Ok(0), // EOF
        Ok(1) => {}
        Ok(_) => unreachable!(),
        Err(ref e) if e.kind() == io::ErrorKind::Interrupted => {
            return read_block_size(reader);
        }
        Err(e) => return Err(e),
    }

    // Read the remaining 3 bytes
    reader.read_exact(&mut buf[1..])?;

    let n = u32::from_le_bytes(buf);
    usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

/// A reader for raw BAM records.
///
/// This wraps any `Read` type (typically a BGZF reader) and provides
/// methods for reading raw BAM records as bytes.
pub struct RawBamReader<R> {
    inner: R,
}

impl<R: Read> RawBamReader<R> {
    /// Creates a new raw BAM reader wrapping the given reader.
    ///
    /// Note: The caller is responsible for reading/skipping the BAM header
    /// before reading records. Use `skip_header` or `read_header_with_noodles`
    /// to handle the header.
    #[inline]
    pub fn new(inner: R) -> Self {
        Self { inner }
    }

    /// Returns a reference to the inner reader.
    #[inline]
    pub fn get_ref(&self) -> &R {
        &self.inner
    }

    /// Returns a mutable reference to the inner reader.
    #[inline]
    pub fn get_mut(&mut self) -> &mut R {
        &mut self.inner
    }

    /// Consumes the reader and returns the inner reader.
    #[inline]
    pub fn into_inner(self) -> R {
        self.inner
    }

    /// Reads a single raw BAM record.
    ///
    /// Returns the number of bytes read, or 0 at EOF.
    ///
    /// # Errors
    ///
    /// Returns an error if reading fails or if EOF is encountered mid-record.
    #[inline]
    pub fn read_record(&mut self, record: &mut RawRecord) -> io::Result<usize> {
        read_raw_record(&mut self.inner, record)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_raw_record_new() {
        let record = RawRecord::new();
        assert!(record.is_empty());
        assert_eq!(record.len(), 0);
    }

    #[test]
    fn test_raw_record_as_ref() {
        let record = RawRecord::from(vec![1, 2, 3, 4]);
        assert_eq!(record.as_ref(), &[1, 2, 3, 4]);
    }

    #[test]
    fn test_raw_record_deref() {
        let record = RawRecord::from(vec![1, 2, 3, 4]);
        assert_eq!(&*record, &[1, 2, 3, 4]);
    }

    #[test]
    fn test_read_raw_record_success() {
        // block_size = 8, followed by 8 bytes of data
        let data = [
            0x08, 0x00, 0x00, 0x00, // block_size = 8
            0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, // record data
        ];
        let mut reader = &data[..];
        let mut record = RawRecord::new();

        let n = read_raw_record(&mut reader, &mut record).unwrap();
        assert_eq!(n, 8);
        assert_eq!(record.as_ref(), &[1, 2, 3, 4, 5, 6, 7, 8]);
    }

    #[test]
    fn test_read_raw_record_eof() {
        let data: &[u8] = &[];
        let mut reader = data;
        let mut record = RawRecord::new();

        let n = read_raw_record(&mut reader, &mut record).unwrap();
        assert_eq!(n, 0);
    }

    #[test]
    fn test_read_raw_record_truncated_size() {
        let data = [0x08, 0x00]; // Incomplete block_size
        let mut reader = &data[..];
        let mut record = RawRecord::new();

        let result = read_raw_record(&mut reader, &mut record);
        assert!(result.is_err());
    }

    #[test]
    fn test_read_raw_record_truncated_data() {
        let data = [
            0x08, 0x00, 0x00, 0x00, // block_size = 8
            0x01, 0x02, 0x03, // Only 3 bytes of data
        ];
        let mut reader = &data[..];
        let mut record = RawRecord::new();

        let result = read_raw_record(&mut reader, &mut record);
        assert!(result.is_err());
    }

    #[test]
    fn test_read_multiple_records() {
        let data = [
            0x04, 0x00, 0x00, 0x00, // block_size = 4
            0x01, 0x02, 0x03, 0x04, // record 1
            0x02, 0x00, 0x00, 0x00, // block_size = 2
            0x05, 0x06, // record 2
        ];
        let mut reader = &data[..];
        let mut record = RawRecord::new();

        let n = read_raw_record(&mut reader, &mut record).unwrap();
        assert_eq!(n, 4);
        assert_eq!(record.as_ref(), &[1, 2, 3, 4]);

        let n = read_raw_record(&mut reader, &mut record).unwrap();
        assert_eq!(n, 2);
        assert_eq!(record.as_ref(), &[5, 6]);

        let n = read_raw_record(&mut reader, &mut record).unwrap();
        assert_eq!(n, 0); // EOF
    }
}

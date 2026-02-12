//! Direct raw BAM record reading for high-performance sorting.
//!
//! This module provides a reader that bypasses noodles Record objects,
//! reading raw BAM bytes directly from BGZF blocks. This eliminates
//! per-record allocation and deserialization overhead.
//!
//! # Architecture
//!
//! ```text
//! ┌─────────────────┐    ┌─────────────────┐    ┌─────────────────┐
//! │  BGZF Blocks    │───>│   Decompress    │───>│  Parse Records  │
//! │   (batched)     │    │  (libdeflater)  │    │   (zero-copy)   │
//! └─────────────────┘    └─────────────────┘    └─────────────────┘
//! ```
//!
//! # Performance
//!
//! Compared to noodles-based reading:
//! - Eliminates Record object allocation per record
//! - Batched I/O reduces syscall overhead
//! - Direct buffer management avoids extra copies
//!
//! # Note
//!
//! This module is currently not integrated into the sort pipeline.
//! It provides infrastructure for future optimization.

use crate::bgzf_reader::{decompress_block_into, read_raw_blocks};
use libdeflater::Decompressor;
use std::io::{self, BufReader, Read};

/// Number of BGZF blocks to read per batch.
const BLOCKS_PER_BATCH: usize = 64;

/// BAM magic bytes (appears in decompressed data, not raw file).
const BAM_MAGIC: &[u8; 4] = b"BAM\x01";

/// A raw BAM record reader that reads directly from BGZF blocks.
///
/// This reader bypasses noodles and reads raw BAM record bytes directly,
/// eliminating per-record allocation overhead.
pub struct RawBamRecordReader<R: Read> {
    /// Buffered reader for the input.
    reader: BufReader<R>,
    /// Decompressor for BGZF blocks.
    decompressor: Decompressor,
    /// Buffer for decompressed data.
    decompressed: Vec<u8>,
    /// Current position in decompressed buffer.
    position: usize,
    /// Whether we've reached EOF.
    eof: bool,
    /// Whether the BAM header has been skipped.
    header_skipped: bool,
}

impl<R: Read> RawBamRecordReader<R> {
    /// Create a new raw BAM record reader.
    ///
    /// Decompresses initial BGZF blocks and verifies BAM magic.
    /// Call `skip_header()` before reading records.
    ///
    /// # Errors
    ///
    /// Returns an error if the input is not a valid BAM file.
    pub fn new(reader: R) -> io::Result<Self> {
        let reader = BufReader::with_capacity(256 * 1024, reader);

        let mut this = Self {
            reader,
            decompressor: Decompressor::new(),
            decompressed: Vec::with_capacity(64 * 65536), // ~4MB for 64 blocks
            position: 0,
            eof: false,
            header_skipped: false,
        };

        // Decompress initial blocks to access BAM header
        this.refill_buffer()?;

        // Verify BAM magic (in decompressed data)
        if this.decompressed.len() < 4 {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "File too small to contain BAM magic",
            ));
        }
        if &this.decompressed[0..4] != BAM_MAGIC {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "Not a BAM file: expected magic {:?}, got {:?}",
                    BAM_MAGIC,
                    &this.decompressed[0..4]
                ),
            ));
        }
        this.position = 4; // Skip past magic

        Ok(this)
    }

    /// Create a new raw BAM record reader from a `BufReader`.
    ///
    /// The reader should be positioned at the start of a BGZF-compressed BAM file.
    ///
    /// # Errors
    ///
    /// Returns an error if the input is not a valid BAM file.
    pub fn from_buf_reader(reader: BufReader<R>) -> io::Result<Self> {
        let mut this = Self {
            reader,
            decompressor: Decompressor::new(),
            decompressed: Vec::with_capacity(64 * 65536),
            position: 0,
            eof: false,
            header_skipped: false,
        };

        // Decompress initial blocks
        this.refill_buffer()?;

        // Verify BAM magic
        if this.decompressed.len() < 4 || &this.decompressed[0..4] != BAM_MAGIC {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "Not a BAM file"));
        }
        this.position = 4;

        Ok(this)
    }

    /// Skip the BAM header and reference sequences.
    ///
    /// Must be called before reading records. Returns the raw header bytes
    /// so they can be parsed if needed.
    ///
    /// # Errors
    ///
    /// Returns an error if the header is already skipped or the header is truncated.
    #[allow(clippy::cast_possible_truncation)]
    pub fn skip_header(&mut self) -> io::Result<Vec<u8>> {
        if self.header_skipped {
            return Err(io::Error::other("Header already skipped"));
        }

        // Ensure we have data
        self.ensure_data()?;

        let mut header_bytes = Vec::new();

        // Read header text length (4 bytes)
        let l_text = self.read_u32()? as usize;
        header_bytes.extend_from_slice(&(l_text as u32).to_le_bytes());

        // Skip header text
        if !self.ensure_bytes(l_text)? {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "Truncated BAM header (text)",
            ));
        }
        header_bytes.extend_from_slice(&self.decompressed[self.position..self.position + l_text]);
        self.position += l_text;

        // Read number of references (4 bytes)
        let n_ref = self.read_u32()? as usize;
        header_bytes.extend_from_slice(&(n_ref as u32).to_le_bytes());

        // Skip reference sequences
        for _ in 0..n_ref {
            // Reference name length
            let l_name = self.read_u32()? as usize;
            header_bytes.extend_from_slice(&(l_name as u32).to_le_bytes());

            // Reference name
            if !self.ensure_bytes(l_name)? {
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "Truncated BAM header (ref name)",
                ));
            }
            header_bytes
                .extend_from_slice(&self.decompressed[self.position..self.position + l_name]);
            self.position += l_name;

            // Reference length
            let l_ref = self.read_u32()?;
            header_bytes.extend_from_slice(&l_ref.to_le_bytes());
        }

        self.header_skipped = true;
        Ok(header_bytes)
    }

    /// Read the next raw BAM record.
    ///
    /// Returns `Some(record_bytes)` if a record is available, `None` at EOF.
    /// The returned bytes exclude the 4-byte `block_size` prefix (matching noodles format).
    ///
    /// # Errors
    ///
    /// Returns an error if the header has not been skipped or the record is truncated.
    pub fn next_record(&mut self) -> io::Result<Option<Vec<u8>>> {
        if !self.header_skipped {
            return Err(io::Error::other("Must call skip_header() first"));
        }

        // Ensure we have at least 4 bytes for block_size
        if !self.ensure_bytes(4)? {
            return Ok(None); // EOF
        }

        // Read block_size (4 bytes, little-endian)
        let block_size = u32::from_le_bytes([
            self.decompressed[self.position],
            self.decompressed[self.position + 1],
            self.decompressed[self.position + 2],
            self.decompressed[self.position + 3],
        ]) as usize;

        // Validate block_size
        if block_size < 32 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Invalid BAM block_size: {block_size}"),
            ));
        }

        // Ensure we have the full record
        let total_size = 4 + block_size;
        if !self.ensure_bytes(total_size)? {
            return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "Truncated BAM record"));
        }

        // Copy record bytes (without the 4-byte block_size prefix, matching noodles format)
        let record = self.decompressed[self.position + 4..self.position + total_size].to_vec();
        self.position += total_size;

        Ok(Some(record))
    }

    /// Read a u32 from the decompressed buffer.
    fn read_u32(&mut self) -> io::Result<u32> {
        if !self.ensure_bytes(4)? {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "Unexpected EOF while reading u32",
            ));
        }
        let val = u32::from_le_bytes([
            self.decompressed[self.position],
            self.decompressed[self.position + 1],
            self.decompressed[self.position + 2],
            self.decompressed[self.position + 3],
        ]);
        self.position += 4;
        Ok(val)
    }

    /// Ensure at least `n` bytes are available in the decompressed buffer.
    ///
    /// Returns `false` if EOF is reached before `n` bytes are available.
    fn ensure_bytes(&mut self, n: usize) -> io::Result<bool> {
        while self.position + n > self.decompressed.len() {
            if self.eof {
                return Ok(false);
            }
            self.refill_buffer()?;
        }
        Ok(true)
    }

    /// Ensure we have some data in the buffer.
    fn ensure_data(&mut self) -> io::Result<()> {
        if self.position >= self.decompressed.len() && !self.eof {
            self.refill_buffer()?;
        }
        Ok(())
    }

    /// Refill the decompressed buffer by reading more BGZF blocks.
    fn refill_buffer(&mut self) -> io::Result<()> {
        // Compact buffer: move remaining data to front
        if self.position > 0 {
            let remaining = self.decompressed.len() - self.position;
            if remaining > 0 {
                self.decompressed.copy_within(self.position.., 0);
            }
            self.decompressed.truncate(remaining);
            self.position = 0;
        }

        // Read more BGZF blocks
        let blocks = read_raw_blocks(&mut self.reader, BLOCKS_PER_BATCH)?;
        if blocks.is_empty() {
            self.eof = true;
            return Ok(());
        }

        // Decompress blocks
        for block in &blocks {
            decompress_block_into(block, &mut self.decompressor, &mut self.decompressed)?;
        }

        Ok(())
    }
}

/// Iterator adapter for `RawBamRecordReader`.
impl<R: Read> Iterator for RawBamRecordReader<R> {
    type Item = io::Result<Vec<u8>>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.next_record() {
            Ok(Some(record)) => Some(Ok(record)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

/// A batched raw BAM record reader for reduced channel overhead.
///
/// Similar to `RawBamRecordReader` but yields batches of records.
pub struct BatchedRawBamReader<R: Read> {
    inner: RawBamRecordReader<R>,
    batch_size: usize,
}

impl<R: Read> BatchedRawBamReader<R> {
    /// Create a new batched reader.
    ///
    /// # Errors
    ///
    /// Returns an error if the input is not a valid BAM file.
    pub fn new(reader: R, batch_size: usize) -> io::Result<Self> {
        Ok(Self { inner: RawBamRecordReader::new(reader)?, batch_size })
    }

    /// Skip the BAM header.
    ///
    /// # Errors
    ///
    /// Returns an error if the header is already skipped or the header is truncated.
    pub fn skip_header(&mut self) -> io::Result<Vec<u8>> {
        self.inner.skip_header()
    }

    /// Read the next batch of records.
    ///
    /// Returns `None` at EOF, otherwise returns a vector of records.
    /// The vector may be smaller than `batch_size` at the end of the file.
    ///
    /// # Errors
    ///
    /// Returns an error if a record is truncated or invalid.
    pub fn next_batch(&mut self) -> io::Result<Option<Vec<Vec<u8>>>> {
        let mut batch = Vec::with_capacity(self.batch_size);

        for _ in 0..self.batch_size {
            match self.inner.next_record()? {
                Some(record) => batch.push(record),
                None => break,
            }
        }

        if batch.is_empty() { Ok(None) } else { Ok(Some(batch)) }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    /// Build a minimal BGZF-compressed BAM with the given records.
    /// Each record is just the raw bytes (without the 4-byte `block_size` prefix).
    #[allow(clippy::cast_possible_truncation)]
    fn build_test_bam(header_text: &str, refs: &[(&str, u32)], records: &[Vec<u8>]) -> Vec<u8> {
        let mut raw_bam = Vec::new();

        // Magic
        raw_bam.extend_from_slice(b"BAM\x01");

        // Header text
        let text_bytes = header_text.as_bytes();
        raw_bam.extend_from_slice(&(text_bytes.len() as u32).to_le_bytes());
        raw_bam.extend_from_slice(text_bytes);

        // References
        raw_bam.extend_from_slice(&(refs.len() as u32).to_le_bytes());
        for (name, length) in refs {
            let name_with_null = format!("{name}\0");
            raw_bam.extend_from_slice(&(name_with_null.len() as u32).to_le_bytes());
            raw_bam.extend_from_slice(name_with_null.as_bytes());
            raw_bam.extend_from_slice(&length.to_le_bytes());
        }

        // Records
        for record in records {
            raw_bam.extend_from_slice(&(record.len() as u32).to_le_bytes());
            raw_bam.extend_from_slice(record);
        }

        // Compress with BGZF
        let mut compressed = Vec::new();
        {
            let mut writer = noodles_bgzf::io::Writer::new(&mut compressed);
            writer.write_all(&raw_bam).unwrap();
            writer.finish().unwrap();
        }
        compressed
    }

    /// Build a minimal BAM record with the given read name.
    #[allow(clippy::cast_possible_truncation)]
    fn make_minimal_record(name: &[u8]) -> Vec<u8> {
        let l_read_name = (name.len() + 1) as u8; // +1 for null terminator
        let total = 32 + l_read_name as usize;
        let mut rec = vec![0u8; total];
        // tid at offset 0-3: -1 (unmapped)
        rec[0..4].copy_from_slice(&(-1i32).to_le_bytes());
        // pos at offset 4-7: -1
        rec[4..8].copy_from_slice(&(-1i32).to_le_bytes());
        // l_read_name at offset 8
        rec[8] = l_read_name;
        // Copy name + null terminator (null is already 0 from vec init)
        rec[32..32 + name.len()].copy_from_slice(name);
        rec
    }

    #[test]
    fn test_blocks_per_batch() {
        assert_eq!(BLOCKS_PER_BATCH, 64);
    }

    #[test]
    fn test_new_valid_bam() {
        let data = build_test_bam("@HD\tVN:1.6\n", &[], &[]);
        let reader = RawBamRecordReader::new(io::Cursor::new(data));
        assert!(reader.is_ok(), "Expected valid BAM to succeed: {:?}", reader.err());
    }

    #[test]
    fn test_new_invalid_magic() {
        // Write non-BAM data into a BGZF stream
        let mut compressed = Vec::new();
        {
            let mut writer = noodles_bgzf::io::Writer::new(&mut compressed);
            writer.write_all(b"NOT_BAM!").unwrap();
            writer.finish().unwrap();
        }
        let result = RawBamRecordReader::new(io::Cursor::new(compressed));
        let Err(err) = result else { panic!("Expected error for invalid magic") };
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert!(err.to_string().contains("Not a BAM file"), "Unexpected error message: {err}");
    }

    #[test]
    fn test_new_empty_input() {
        // An empty BGZF stream: just the EOF block
        let mut compressed = Vec::new();
        {
            let writer = noodles_bgzf::io::Writer::new(&mut compressed);
            // Write nothing, just finish to produce the EOF block
            writer.finish().unwrap();
        }
        let result = RawBamRecordReader::new(io::Cursor::new(compressed));
        let Err(err) = result else { panic!("Expected error for empty input") };
        assert!(
            err.to_string().contains("File too small")
                || err.to_string().contains("Not a BAM file"),
            "Unexpected error message: {err}"
        );
    }

    #[test]
    fn test_skip_header_returns_bytes() {
        let header_text = "@HD\tVN:1.6\n";
        let refs = vec![("chr1", 1000u32), ("chr2", 2000u32)];
        let data = build_test_bam(header_text, &refs, &[]);
        let mut reader = RawBamRecordReader::new(io::Cursor::new(data)).unwrap();
        let header_bytes = reader.skip_header().unwrap();

        // The header_bytes should contain l_text + header_text + n_ref + ref data
        // Verify we can parse it back
        assert!(!header_bytes.is_empty());

        // First 4 bytes: l_text
        let l_text = u32::from_le_bytes(header_bytes[0..4].try_into().unwrap()) as usize;
        assert_eq!(l_text, header_text.len());

        // Then header text
        let parsed_text = &header_bytes[4..4 + l_text];
        assert_eq!(parsed_text, header_text.as_bytes());

        // Then n_ref
        let offset = 4 + l_text;
        let n_ref = u32::from_le_bytes(header_bytes[offset..offset + 4].try_into().unwrap());
        assert_eq!(n_ref, 2);
    }

    #[test]
    fn test_skip_header_twice_errors() {
        let data = build_test_bam("@HD\tVN:1.6\n", &[], &[]);
        let mut reader = RawBamRecordReader::new(io::Cursor::new(data)).unwrap();
        reader.skip_header().unwrap();
        let result = reader.skip_header();
        assert!(result.is_err());
        assert!(
            result.unwrap_err().to_string().contains("already skipped"),
            "Expected 'already skipped' error"
        );
    }

    #[test]
    fn test_next_record_without_skip_header_errors() {
        let data = build_test_bam("@HD\tVN:1.6\n", &[], &[]);
        let mut reader = RawBamRecordReader::new(io::Cursor::new(data)).unwrap();
        let result = reader.next_record();
        assert!(result.is_err());
        assert!(
            result.unwrap_err().to_string().contains("skip_header"),
            "Expected error about skip_header"
        );
    }

    #[test]
    fn test_read_single_record() {
        let rec = make_minimal_record(b"R");
        let data = build_test_bam("@HD\tVN:1.6\n", &[], std::slice::from_ref(&rec));
        let mut reader = RawBamRecordReader::new(io::Cursor::new(data)).unwrap();
        reader.skip_header().unwrap();

        let record = reader.next_record().unwrap();
        assert!(record.is_some(), "Expected one record");
        let record = record.unwrap();
        assert_eq!(record, rec, "Record bytes should match");

        // No more records
        let eof = reader.next_record().unwrap();
        assert!(eof.is_none(), "Expected EOF after single record");
    }

    #[test]
    fn test_read_multiple_records() {
        let rec_a = make_minimal_record(b"A");
        let rec_b = make_minimal_record(b"B");
        let rec_c = make_minimal_record(b"C");
        let data =
            build_test_bam("@HD\tVN:1.6\n", &[], &[rec_a.clone(), rec_b.clone(), rec_c.clone()]);
        let mut reader = RawBamRecordReader::new(io::Cursor::new(data)).unwrap();
        reader.skip_header().unwrap();

        let r1 = reader.next_record().unwrap().expect("record 1");
        let r2 = reader.next_record().unwrap().expect("record 2");
        let r3 = reader.next_record().unwrap().expect("record 3");

        assert_eq!(r1, rec_a);
        assert_eq!(r2, rec_b);
        assert_eq!(r3, rec_c);

        assert!(reader.next_record().unwrap().is_none(), "Expected EOF");
    }

    #[test]
    fn test_iterator_adapter() {
        let rec_a = make_minimal_record(b"X");
        let rec_b = make_minimal_record(b"Y");
        let data = build_test_bam("@HD\tVN:1.6\n", &[], &[rec_a.clone(), rec_b.clone()]);
        let mut reader = RawBamRecordReader::new(io::Cursor::new(data)).unwrap();
        reader.skip_header().unwrap();

        let records: Vec<Vec<u8>> = reader.map(|r| r.unwrap()).collect();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0], rec_a);
        assert_eq!(records[1], rec_b);
    }

    #[test]
    fn test_batched_reader() {
        let rec_a = make_minimal_record(b"A");
        let rec_b = make_minimal_record(b"B");
        let rec_c = make_minimal_record(b"C");
        let data =
            build_test_bam("@HD\tVN:1.6\n", &[], &[rec_a.clone(), rec_b.clone(), rec_c.clone()]);
        let mut reader = BatchedRawBamReader::new(io::Cursor::new(data), 2).unwrap();
        reader.skip_header().unwrap();

        // First batch: 2 records
        let batch1 = reader.next_batch().unwrap().expect("batch 1");
        assert_eq!(batch1.len(), 2);
        assert_eq!(batch1[0], rec_a);
        assert_eq!(batch1[1], rec_b);

        // Second batch: 1 record (remainder)
        let batch2 = reader.next_batch().unwrap().expect("batch 2");
        assert_eq!(batch2.len(), 1);
        assert_eq!(batch2[0], rec_c);

        // No more batches
        assert!(reader.next_batch().unwrap().is_none(), "Expected no more batches");
    }

    #[test]
    fn test_from_buf_reader() {
        let data = build_test_bam("@HD\tVN:1.6\n", &[("chr1", 500)], &[]);
        let buf_reader = BufReader::new(io::Cursor::new(data));
        let mut reader = RawBamRecordReader::from_buf_reader(buf_reader).unwrap();
        let header_bytes = reader.skip_header().unwrap();
        assert!(!header_bytes.is_empty());

        // Verify no records
        assert!(reader.next_record().unwrap().is_none());
    }

    #[test]
    fn test_batched_reader_eof() {
        let data = build_test_bam("@HD\tVN:1.6\n", &[], &[]);
        let mut reader = BatchedRawBamReader::new(io::Cursor::new(data), 10).unwrap();
        reader.skip_header().unwrap();

        // No records, should return None immediately
        assert!(reader.next_batch().unwrap().is_none(), "Expected None for empty BAM");

        // Calling again should still return None
        assert!(reader.next_batch().unwrap().is_none(), "Expected None on repeated call");
    }
}

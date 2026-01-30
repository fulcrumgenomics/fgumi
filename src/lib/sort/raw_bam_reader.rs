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
        self.ensure_bytes(l_text)?;
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
            self.ensure_bytes(l_name)?;
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
    /// The returned bytes include the full BAM record (`block_size` + data).
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
        self.ensure_bytes(4)?;
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
    pub fn new(reader: R, batch_size: usize) -> io::Result<Self> {
        Ok(Self { inner: RawBamRecordReader::new(reader)?, batch_size })
    }

    /// Skip the BAM header.
    pub fn skip_header(&mut self) -> io::Result<Vec<u8>> {
        self.inner.skip_header()
    }

    /// Read the next batch of records.
    ///
    /// Returns `None` at EOF, otherwise returns a vector of records.
    /// The vector may be smaller than `batch_size` at the end of the file.
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

    #[test]
    fn test_blocks_per_batch() {
        assert_eq!(BLOCKS_PER_BATCH, 64);
    }
}

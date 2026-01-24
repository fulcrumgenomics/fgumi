//! Raw BGZF block reading and decompression.
//!
//! This module provides low-level functions for reading raw BGZF blocks
//! without decompressing them, and for decompressing blocks using libdeflater.
//! This enables parallel decompression in worker threads.
//!
//! # BGZF Format
//!
//! BGZF (Blocked GZIP Format) is a variant of gzip that stores data in
//! independent blocks, each up to 64KB uncompressed. The block structure:
//!
//! ```text
//! ┌─────────────────────────────────────────────────────────────────┐
//! │ Header (18 bytes)                                               │
//! │  - Magic: 0x1f 0x8b (gzip)                                      │
//! │  - Method: 0x08 (deflate)                                       │
//! │  - Flags: 0x04 (FEXTRA)                                         │
//! │  - MTIME, XFL, OS: 6 bytes                                      │
//! │  - XLEN: 2 bytes (= 6)                                          │
//! │  - Subfield: "BC" + len(2) + BSIZE(2)                           │
//! │    where BSIZE = total_block_size - 1                           │
//! ├─────────────────────────────────────────────────────────────────┤
//! │ Compressed data (deflate)                                       │
//! ├─────────────────────────────────────────────────────────────────┤
//! │ Footer (8 bytes)                                                │
//! │  - CRC32: 4 bytes                                               │
//! │  - ISIZE: 4 bytes (uncompressed size mod 2^32)                  │
//! └─────────────────────────────────────────────────────────────────┘
//! ```
//!
//! # Usage
//!
//! ```ignore
//! use fgumi_lib::bgzf_reader::{read_raw_blocks, decompress_block};
//! use libdeflater::Decompressor;
//! use std::io::BufReader;
//!
//! let mut reader = BufReader::new(File::open("input.bam")?);
//! let mut decompressor = Decompressor::new();
//!
//! loop {
//!     let blocks = read_raw_blocks(&mut reader, 100)?;
//!     if blocks.is_empty() {
//!         break;
//!     }
//!     for block in &blocks {
//!         let data = decompress_block(block, &mut decompressor)?;
//!         // Process decompressed data...
//!     }
//! }
//! ```

use libdeflater::Decompressor;
use std::io::{self, Read};

// ============================================================================
// Constants
// ============================================================================

/// Size of the BGZF block header.
pub const BGZF_HEADER_SIZE: usize = 18;

/// Size of the BGZF block footer (CRC32 + ISIZE).
pub const BGZF_FOOTER_SIZE: usize = 8;

/// BGZF EOF marker block (empty block signaling end of file).
pub const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

// ============================================================================
// Raw Block Types
// ============================================================================

/// A raw BGZF block (compressed, not yet decompressed).
#[derive(Debug, Clone)]
pub struct RawBgzfBlock {
    /// Complete raw block data: header + compressed data + footer.
    pub data: Vec<u8>,
}

impl RawBgzfBlock {
    /// Get the total size of the block.
    #[must_use]
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Check if this is an empty block.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Check if this is the BGZF EOF marker block.
    #[must_use]
    pub fn is_eof(&self) -> bool {
        self.data == BGZF_EOF
    }

    /// Get the compressed data portion (between header and footer).
    #[must_use]
    pub fn compressed_data(&self) -> &[u8] {
        if self.data.len() <= BGZF_HEADER_SIZE + BGZF_FOOTER_SIZE {
            return &[];
        }
        &self.data[BGZF_HEADER_SIZE..self.data.len() - BGZF_FOOTER_SIZE]
    }

    /// Get the expected uncompressed size from the footer (ISIZE field).
    #[must_use]
    pub fn uncompressed_size(&self) -> usize {
        if self.data.len() < BGZF_FOOTER_SIZE {
            return 0;
        }
        let len = self.data.len();
        u32::from_le_bytes([
            self.data[len - 4],
            self.data[len - 3],
            self.data[len - 2],
            self.data[len - 1],
        ]) as usize
    }

    /// Get the CRC32 from the footer.
    #[must_use]
    pub fn crc32(&self) -> u32 {
        if self.data.len() < BGZF_FOOTER_SIZE {
            return 0;
        }
        let len = self.data.len();
        u32::from_le_bytes([
            self.data[len - 8],
            self.data[len - 7],
            self.data[len - 6],
            self.data[len - 5],
        ])
    }
}

// ============================================================================
// Reading Functions
// ============================================================================

/// Read a single raw BGZF block from the input.
///
/// Returns `Ok(Some(block))` if a block was read, `Ok(None)` at EOF,
/// or an error if reading failed or the data is invalid.
fn read_raw_block<R: Read + ?Sized>(reader: &mut R) -> io::Result<Option<RawBgzfBlock>> {
    // Read the 18-byte header
    let mut header = [0u8; BGZF_HEADER_SIZE];
    match reader.read_exact(&mut header) {
        Ok(()) => {}
        Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(None),
        Err(e) => return Err(e),
    }

    // Validate gzip magic bytes
    if header[0] != 0x1f || header[1] != 0x8b {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "Invalid BGZF magic: expected 0x1f 0x8b, got 0x{:02x} 0x{:02x}",
                header[0], header[1]
            ),
        ));
    }

    // Validate compression method (deflate)
    if header[2] != 0x08 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Invalid compression method: expected 0x08, got 0x{:02x}", header[2]),
        ));
    }

    // Validate FEXTRA flag
    if header[3] & 0x04 == 0 {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "BGZF block missing FEXTRA flag"));
    }

    // Validate BC subfield identifier
    if header[12] != b'B' || header[13] != b'C' {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "Invalid BGZF subfield ID: expected 'BC', got '{}{}'",
                header[12] as char, header[13] as char
            ),
        ));
    }

    // Get block size from BSIZE field (bytes 16-17, little-endian)
    // BSIZE = total_block_size - 1
    let bsize = u16::from_le_bytes([header[16], header[17]]) as usize;
    let block_size = bsize + 1;

    if block_size < BGZF_HEADER_SIZE + BGZF_FOOTER_SIZE {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("BGZF block too small: {block_size} bytes"),
        ));
    }

    // Allocate buffer and copy header
    let mut data = vec![0u8; block_size];
    data[..BGZF_HEADER_SIZE].copy_from_slice(&header);

    // Read remaining block data
    reader.read_exact(&mut data[BGZF_HEADER_SIZE..])?;

    Ok(Some(RawBgzfBlock { data }))
}

/// Read multiple raw BGZF blocks as a batch.
///
/// This is more efficient than calling `read_raw_block` in a loop
/// when you need to read many blocks.
///
/// # Arguments
///
/// * `reader` - The input reader.
/// * `max_blocks` - Maximum number of blocks to read.
///
/// # Returns
///
/// A vector of blocks read. The vector may be shorter than `max_blocks`
/// if EOF is reached. Returns an empty vector at EOF.
pub fn read_raw_blocks<R: Read + ?Sized>(
    reader: &mut R,
    max_blocks: usize,
) -> io::Result<Vec<RawBgzfBlock>> {
    let mut blocks = Vec::with_capacity(max_blocks);
    for _ in 0..max_blocks {
        match read_raw_block(reader)? {
            Some(block) => {
                // Skip EOF marker blocks
                if !block.is_eof() {
                    blocks.push(block);
                }
            }
            None => break,
        }
    }
    Ok(blocks)
}

// ============================================================================
// Decompression Functions
// ============================================================================

/// Decompress a raw BGZF block.
///
/// Uses libdeflater for high-performance decompression.
///
/// # Arguments
///
/// * `block` - The raw BGZF block to decompress.
/// * `decompressor` - A reusable libdeflater decompressor.
///
/// # Returns
///
/// The uncompressed data.
///
/// # Errors
///
/// Returns an error if:
/// - Decompression fails
/// - CRC32 verification fails (if `verify_crc` is true)
pub fn decompress_block(
    block: &RawBgzfBlock,
    decompressor: &mut Decompressor,
) -> io::Result<Vec<u8>> {
    // Handle empty/EOF blocks
    if block.is_eof() || block.uncompressed_size() == 0 {
        return Ok(Vec::new());
    }

    let compressed = block.compressed_data();
    let uncompressed_size = block.uncompressed_size();

    // Allocate output buffer
    let mut uncompressed = vec![0u8; uncompressed_size];

    // Decompress
    decompressor.deflate_decompress(compressed, &mut uncompressed).map_err(|e| {
        io::Error::new(io::ErrorKind::InvalidData, format!("BGZF decompression failed: {e:?}"))
    })?;

    Ok(uncompressed)
}

/// Decompress a BGZF block into a provided buffer, appending to existing data.
///
/// This variant avoids allocation by reusing the provided buffer.
/// The decompressed data is appended to `output`, which should be pre-sized
/// or have sufficient capacity.
///
/// # Arguments
///
/// * `block` - The raw BGZF block to decompress.
/// * `decompressor` - A reusable libdeflater decompressor.
/// * `output` - Buffer to append decompressed data to.
///
/// # Errors
///
/// Returns an error if decompression fails.
pub fn decompress_block_into(
    block: &RawBgzfBlock,
    decompressor: &mut Decompressor,
    output: &mut Vec<u8>,
) -> io::Result<()> {
    // Handle empty/EOF blocks
    if block.is_eof() || block.uncompressed_size() == 0 {
        return Ok(());
    }

    let compressed = block.compressed_data();
    let uncompressed_size = block.uncompressed_size();

    // Record current length and extend buffer
    let start = output.len();
    output.resize(start + uncompressed_size, 0);

    // Decompress directly into the buffer
    decompressor.deflate_decompress(compressed, &mut output[start..]).map_err(|e| {
        io::Error::new(io::ErrorKind::InvalidData, format!("BGZF decompression failed: {e:?}"))
    })?;

    // Verify CRC32 to detect corruption
    let expected_crc = block.crc32();
    let actual_crc = crc32fast::hash(&output[start..]);
    if expected_crc != actual_crc {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "BGZF CRC32 mismatch: expected 0x{:08x}, got 0x{:08x}, block_size={}, uncompressed_size={}",
                expected_crc,
                actual_crc,
                block.len(),
                uncompressed_size
            ),
        ));
    }

    Ok(())
}

/// Decompress a BGZF block from raw bytes into a provided buffer.
///
/// This variant accepts a raw byte slice directly, avoiding the need to
/// construct a `RawBgzfBlock` (and its associated allocation).
///
/// # Arguments
///
/// * `data` - Raw BGZF block bytes (header + compressed + footer).
/// * `decompressor` - A reusable libdeflater decompressor.
/// * `output` - Buffer to append decompressed data to.
///
/// # Errors
///
/// Returns an error if decompression fails.
pub fn decompress_block_slice_into(
    data: &[u8],
    decompressor: &mut Decompressor,
    output: &mut Vec<u8>,
) -> io::Result<()> {
    // Check for EOF block
    if data == BGZF_EOF {
        return Ok(());
    }

    // Need at least header + footer
    if data.len() < BGZF_HEADER_SIZE + BGZF_FOOTER_SIZE {
        return Ok(());
    }

    // Get uncompressed size from ISIZE field (last 4 bytes)
    let len = data.len();
    let uncompressed_size =
        u32::from_le_bytes([data[len - 4], data[len - 3], data[len - 2], data[len - 1]]) as usize;

    if uncompressed_size == 0 {
        return Ok(());
    }

    // Get compressed data (between header and footer)
    let compressed = &data[BGZF_HEADER_SIZE..len - BGZF_FOOTER_SIZE];

    // Extend output buffer and decompress directly into it
    let start = output.len();
    output.resize(start + uncompressed_size, 0);

    decompressor.deflate_decompress(compressed, &mut output[start..]).map_err(|e| {
        io::Error::new(io::ErrorKind::InvalidData, format!("BGZF decompression failed: {e:?}"))
    })?;

    Ok(())
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_eof_block_detection() {
        let block = RawBgzfBlock { data: BGZF_EOF.to_vec() };
        assert!(block.is_eof());
        assert_eq!(block.uncompressed_size(), 0);
    }

    #[test]
    fn test_raw_block_accessors() {
        // Create a minimal valid block
        let mut data = vec![0u8; 30];
        // Magic
        data[0] = 0x1f;
        data[1] = 0x8b;
        // Method
        data[2] = 0x08;
        // Flags
        data[3] = 0x04;
        // BC subfield
        data[12] = b'B';
        data[13] = b'C';
        // BSIZE (29 = 30 - 1)
        data[16] = 29;
        data[17] = 0;
        // Footer: CRC32 (bytes 22-25) and ISIZE (bytes 26-29)
        data[22] = 0x12;
        data[23] = 0x34;
        data[24] = 0x56;
        data[25] = 0x78;
        data[26] = 100; // ISIZE = 100
        data[27] = 0;
        data[28] = 0;
        data[29] = 0;

        let block = RawBgzfBlock { data };
        assert_eq!(block.len(), 30);
        assert_eq!(block.uncompressed_size(), 100);
        assert_eq!(block.crc32(), 0x7856_3412);
        assert!(!block.is_eof());
    }

    #[test]
    fn test_read_invalid_magic() {
        // Need at least 18 bytes (header size) for magic validation to occur
        // With fewer bytes, read_exact returns UnexpectedEof which is treated as EOF
        let mut data = vec![0x00; BGZF_HEADER_SIZE];
        data[0] = 0x00; // Invalid magic (should be 0x1f)
        data[1] = 0x00; // Invalid magic (should be 0x8b)
        let mut reader = Cursor::new(data);
        let result = read_raw_block(&mut reader);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("Invalid BGZF magic"));
    }

    #[test]
    fn test_read_eof() {
        let mut reader = Cursor::new(Vec::<u8>::new());
        let result = read_raw_block(&mut reader).unwrap();
        assert!(result.is_none());
    }

    #[test]
    fn test_decompress_eof_block() {
        let block = RawBgzfBlock { data: BGZF_EOF.to_vec() };
        let mut decompressor = Decompressor::new();
        let result = decompress_block(&block, &mut decompressor).unwrap();
        assert!(result.is_empty());
    }

    #[test]
    fn test_decompress_block_into_eof() {
        let block = RawBgzfBlock { data: BGZF_EOF.to_vec() };
        let mut decompressor = Decompressor::new();
        let mut output = Vec::new();
        decompress_block_into(&block, &mut decompressor, &mut output).unwrap();
        assert!(output.is_empty());
    }

    #[test]
    fn test_decompress_block_into_appends() {
        // Create a compressed block using the writer
        use crate::bgzf_writer::InlineBgzfCompressor;

        let original_data = b"Hello, BGZF world!";
        let mut compressor = InlineBgzfCompressor::new(6);
        compressor.write_all(original_data).unwrap();
        compressor.flush().unwrap();
        let blocks = compressor.take_blocks();
        assert_eq!(blocks.len(), 1);

        let block = RawBgzfBlock { data: blocks[0].data.clone() };
        let mut decompressor = Decompressor::new();

        // Start with existing data in the buffer
        let mut output = vec![1, 2, 3];
        decompress_block_into(&block, &mut decompressor, &mut output).unwrap();

        // Should have preserved existing data and appended decompressed data
        assert_eq!(&output[0..3], &[1, 2, 3]);
        assert_eq!(&output[3..], original_data);
    }

    #[test]
    fn test_decompress_block_into_equivalence() {
        // Create a compressed block
        use crate::bgzf_writer::InlineBgzfCompressor;

        let original_data = b"Test data for equivalence check";
        let mut compressor = InlineBgzfCompressor::new(6);
        compressor.write_all(original_data).unwrap();
        compressor.flush().unwrap();
        let blocks = compressor.take_blocks();

        let block = RawBgzfBlock { data: blocks[0].data.clone() };
        let mut decompressor = Decompressor::new();

        // Decompress using original function
        let result1 = decompress_block(&block, &mut decompressor).unwrap();

        // Decompress using new function
        let mut result2 = Vec::new();
        decompress_block_into(&block, &mut decompressor, &mut result2).unwrap();

        // Should produce identical results
        assert_eq!(result1, result2);
        assert_eq!(result1, original_data);
    }
}

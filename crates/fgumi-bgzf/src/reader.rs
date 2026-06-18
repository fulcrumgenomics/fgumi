//! Raw BGZF block reading and decompression.
//!
//! This module provides low-level functions for reading raw BGZF blocks
//! without decompressing them, and for decompressing blocks using libdeflater.
//! For deflate **stored** blocks (`BTYPE = 00`), the payload is copied
//! directly into the output buffer without invoking libdeflater — see
//! [`decompress_block_into`]. This separation also enables parallel
//! decompression in worker threads.
//!
//! # BGZF Format
//!
//! BGZF (Blocked GZIP Format) is a variant of gzip that stores data in
//! independent blocks, each up to 64KB uncompressed. The block structure:
//!
//! ```text
//! +-------------------------------------------------------------------+
//! | Header (18 bytes)                                                 |
//! |  - Magic: 0x1f 0x8b (gzip)                                       |
//! |  - Method: 0x08 (deflate)                                         |
//! |  - Flags: 0x04 (FEXTRA)                                          |
//! |  - MTIME, XFL, OS: 6 bytes                                       |
//! |  - XLEN: 2 bytes (= 6)                                           |
//! |  - Subfield: "BC" + len(2) + BSIZE(2)                            |
//! |    where BSIZE = total_block_size - 1                             |
//! +-------------------------------------------------------------------+
//! | Compressed data (deflate)                                         |
//! +-------------------------------------------------------------------+
//! | Footer (8 bytes)                                                  |
//! |  - CRC32: 4 bytes                                                 |
//! |  - ISIZE: 4 bytes (uncompressed size mod 2^32)                    |
//! +-------------------------------------------------------------------+
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
// Footer-parsing helpers (operate on raw byte slices)
// ============================================================================

/// Get the compressed data portion from a raw BGZF block slice (between header and footer).
#[must_use]
fn compressed_data_from_slice(data: &[u8]) -> &[u8] {
    if data.len() <= BGZF_HEADER_SIZE + BGZF_FOOTER_SIZE {
        return &[];
    }
    &data[BGZF_HEADER_SIZE..data.len() - BGZF_FOOTER_SIZE]
}

/// Get the expected uncompressed size (ISIZE field) from the footer of a raw BGZF block slice.
///
/// Returns 0 if the slice is too short to contain a footer.
#[must_use]
fn uncompressed_size_from_slice(data: &[u8]) -> usize {
    if data.len() < BGZF_FOOTER_SIZE {
        return 0;
    }
    let len = data.len();
    // ISIZE is always < 64KB for BGZF blocks, fits in usize on all platforms
    u32::from_le_bytes([data[len - 4], data[len - 3], data[len - 2], data[len - 1]]) as usize
}

/// Get the CRC32 from the footer of a raw BGZF block slice.
///
/// Returns 0 if the slice is too short to contain a footer.
#[must_use]
fn crc32_from_slice(data: &[u8]) -> u32 {
    if data.len() < BGZF_FOOTER_SIZE {
        return 0;
    }
    let len = data.len();
    u32::from_le_bytes([data[len - 8], data[len - 7], data[len - 6], data[len - 5]])
}

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
        compressed_data_from_slice(&self.data)
    }

    /// Get the expected uncompressed size from the footer (ISIZE field).
    #[must_use]
    pub fn uncompressed_size(&self) -> usize {
        uncompressed_size_from_slice(&self.data)
    }

    /// Get the CRC32 from the footer.
    #[must_use]
    pub fn crc32(&self) -> u32 {
        crc32_from_slice(&self.data)
    }
}

// ============================================================================
// Reading Functions
// ============================================================================

/// Read a single raw BGZF block from the input.
///
/// Returns `Ok(Some(block))` if a block was read, `Ok(None)` at EOF,
/// or an error if reading failed or the data is invalid.
///
/// # Errors
///
/// Returns an error if the block header is invalid or reading fails.
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
    // BSIZE is u16, so block_size fits in usize on all platforms
    let bsize = usize::from(u16::from_le_bytes([header[16], header[17]]));
    let block_size = bsize + 1;

    if block_size < BGZF_HEADER_SIZE + BGZF_FOOTER_SIZE {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("BGZF block too small: {block_size} bytes"),
        ));
    }

    // Reserve the exact block size and append the header, then read the
    // remaining bytes via `read_to_end` into the uninitialized spare capacity.
    // This avoids the `vec![0u8; block_size]` zero-fill that `read_exact` would
    // immediately overwrite (std fills the spare capacity without a memset). The
    // length check below restores `read_exact`'s "error on a short block"
    // semantics, since `read_to_end` itself returns `Ok` on early EOF.
    let mut data = Vec::with_capacity(block_size);
    data.extend_from_slice(&header);
    let remaining = block_size - BGZF_HEADER_SIZE;
    let read = reader.take(remaining as u64).read_to_end(&mut data)?;
    if read != remaining {
        return Err(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            format!("truncated BGZF block: expected {remaining} more bytes, got {read}"),
        ));
    }

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
///
/// # Errors
///
/// Returns an error if reading or validation fails.
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
// Decompression Helpers
// ============================================================================

/// Verify that decompressed data matches the expected size and CRC32 checksum.
///
/// # Arguments
///
/// * `decompressed` - The decompressed data to verify.
/// * `expected_size` - The expected uncompressed size from the BGZF footer.
/// * `expected_crc` - The expected CRC32 checksum from the BGZF footer.
/// * `block_len` - The total block length (for error messages).
///
/// # Errors
///
/// Returns an error if the size or CRC32 does not match.
fn verify_decompression(
    decompressed: &[u8],
    expected_size: usize,
    expected_crc: u32,
    block_len: usize,
) -> io::Result<()> {
    let actual_size = decompressed.len();
    if actual_size != expected_size {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "BGZF decompression size mismatch: expected {expected_size}, got {actual_size}"
            ),
        ));
    }

    let actual_crc = crc32fast::hash(decompressed);
    if expected_crc != actual_crc {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "BGZF CRC32 mismatch: expected 0x{expected_crc:08x}, got 0x{actual_crc:08x}, \
                 block_size={block_len}, uncompressed_size={expected_size}",
            ),
        ));
    }

    Ok(())
}

// ============================================================================
// Decompression Functions
// ============================================================================

/// Decompress a raw BGZF block.
///
/// Uses libdeflater for high-performance decompression, with a fast path that
/// skips libdeflater for deflate stored blocks (`BTYPE = 00`) — see
/// [`decompress_block_into`] for details.
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
/// Returns an error if decompression fails or CRC32 verification fails.
pub fn decompress_block(
    block: &RawBgzfBlock,
    decompressor: &mut Decompressor,
) -> io::Result<Vec<u8>> {
    // Pre-size to one BGZF block. Avoids the multi-step `RawVec::finish_grow`
    // we see in profiles (~1.9% of decompress CPU on the merge hot path):
    // starting from capacity 0, the resize inside `decompress_and_verify` has
    // to grow through several capacity classes before reaching the final
    // ~64 KiB. A single `with_capacity` lands at the right size on the first
    // allocation.
    let mut output = Vec::with_capacity(crate::BGZF_MAX_BLOCK_SIZE);
    decompress_block_into(block, decompressor, &mut output)?;
    Ok(output)
}

/// Decompress a BGZF block into a provided buffer, appending to existing data.
///
/// This variant avoids allocation by reusing the provided buffer.
/// The decompressed data is appended to `output`, which should be pre-sized
/// or have sufficient capacity.
///
/// Blocks whose deflate payload is a single **stored** sub-block
/// (`BTYPE = 00`) take a fast path: the LEN bytes of payload are copied
/// straight into `output`, skipping the libdeflater call. Real producers
/// that emit stored blocks include `samtools view -u`, htsjdk's level-0
/// writer, and our own [`InlineBgzfCompressor::new(0)`](crate::writer::InlineBgzfCompressor::new).
/// CRC32 verification against the BGZF footer still runs on the copied bytes;
/// a CRC or size mismatch leaves `output` rolled back to its original length.
///
/// EOF marker blocks and blocks with `ISIZE == 0` are no-ops.
///
/// # Arguments
///
/// * `block` - The raw BGZF block to decompress.
/// * `decompressor` - A reusable libdeflater decompressor.
/// * `output` - Buffer to append decompressed data to.
///
/// # Errors
///
/// Returns an error if decompression fails or CRC32 verification fails.
pub fn decompress_block_into(
    block: &RawBgzfBlock,
    decompressor: &mut Decompressor,
    output: &mut Vec<u8>,
) -> io::Result<()> {
    if block.is_eof() || block.uncompressed_size() == 0 {
        return Ok(());
    }

    decompress_and_verify(
        block.compressed_data(),
        block.uncompressed_size(),
        block.crc32(),
        block.len(),
        decompressor,
        output,
    )
}

/// Decompress a BGZF block from raw bytes into a provided buffer.
///
/// This variant accepts a raw byte slice directly, avoiding the need to
/// construct a `RawBgzfBlock` (and its associated allocation). Like
/// [`decompress_block_into`], it transparently picks the stored-block fast
/// path for deflate `BTYPE = 00` blocks.
///
/// Slices that are too short to contain a BGZF header + footer, or whose
/// footer reports `ISIZE == 0`, are no-ops.
///
/// # Arguments
///
/// * `data` - Raw BGZF block bytes (header + compressed + footer).
/// * `decompressor` - A reusable libdeflater decompressor.
/// * `output` - Buffer to append decompressed data to.
///
/// # Errors
///
/// Returns an error if decompression fails or CRC32 verification fails.
pub fn decompress_block_slice_into(
    data: &[u8],
    decompressor: &mut Decompressor,
    output: &mut Vec<u8>,
) -> io::Result<()> {
    if data.len() < BGZF_HEADER_SIZE + BGZF_FOOTER_SIZE {
        return Ok(());
    }

    let uncompressed_size = uncompressed_size_from_slice(data);
    if uncompressed_size == 0 {
        return Ok(());
    }

    decompress_and_verify(
        compressed_data_from_slice(data),
        uncompressed_size,
        crc32_from_slice(data),
        data.len(),
        decompressor,
        output,
    )
}

/// Decompress (or copy, for deflate stored blocks) BGZF block data into the
/// output buffer and verify the result.
///
/// This is the shared implementation for both [`decompress_block_into`] and
/// [`decompress_block_slice_into`], consolidating the decompress + resize +
/// verify logic.
///
/// Includes a fast path for deflate **stored** blocks (`BTYPE = 00`): the
/// uncompressed payload is copied straight into `output` (delegated to
/// [`copy_stored_and_verify`]), bypassing libdeflater. Producers like
/// `samtools view -u`, htsjdk's level-0 writer, and our own
/// [`InlineBgzfCompressor::new(0)`](crate::writer::InlineBgzfCompressor::new)
/// emit these. CRC32 verification still runs on the copied bytes.
///
/// On any error (decompress failure, size mismatch, CRC mismatch, malformed
/// stored framing), `output` is rolled back to its original length so partial
/// data never leaks to callers.
///
/// # Arguments
///
/// * `compressed` - The deflate-framed payload between the BGZF header and
///   footer (i.e. what [`RawBgzfBlock::compressed_data`] returns).
/// * `uncompressed_size` - The expected decompressed length, from the BGZF
///   footer's ISIZE field.
/// * `expected_crc` - The expected CRC32, from the BGZF footer.
/// * `block_len` - Total BGZF block length, used only for error messages.
/// * `decompressor` - Reusable libdeflater decompressor (unused on the
///   stored-block path).
/// * `output` - Buffer that the decompressed bytes are appended to.
fn decompress_and_verify(
    compressed: &[u8],
    uncompressed_size: usize,
    expected_crc: u32,
    block_len: usize,
    decompressor: &mut Decompressor,
    output: &mut Vec<u8>,
) -> io::Result<()> {
    // Stored-block fast path. The deflate frame for a stored block is:
    //   byte 0   : BFINAL | BTYPE | (5 padding bits to next byte boundary)
    //   bytes 1-2: LEN  (little-endian, u16)
    //   bytes 3-4: NLEN (one's complement of LEN; not checked here — the
    //              BGZF footer's CRC32/ISIZE are authoritative)
    //   bytes 5..: LEN bytes of uncompressed payload
    // BTYPE lives in bits 1-2 of the first byte, so `b & 0b110 == 0` means
    // "stored". `payload_len == LEN + 5` is the structural guarantee that
    // there's exactly one stored sub-block spanning the BGZF payload — the
    // form every real level-0 producer emits.
    if !compressed.is_empty() && compressed[0] & 0b110 == 0 {
        return copy_stored_and_verify(
            compressed,
            uncompressed_size,
            expected_crc,
            block_len,
            output,
        );
    }

    let start = output.len();
    output.resize(start + uncompressed_size, 0);

    let result = (|| {
        let bytes_written =
            decompressor.deflate_decompress(compressed, &mut output[start..]).map_err(|e| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("BGZF decompression failed: {e:?}"),
                )
            })?;

        verify_decompression(
            &output[start..start + bytes_written],
            uncompressed_size,
            expected_crc,
            block_len,
        )
    })();

    if result.is_err() {
        output.truncate(start);
    }
    result
}

/// Decompress a deflate **stored** sub-block by copying its payload directly
/// into `output`, skipping libdeflater.
///
/// The caller is responsible for confirming that `compressed[0]` has
/// `BTYPE = 00` before calling — this function then validates the rest of
/// the stored framing:
///
/// * `compressed.len() >= 5` (room for the BFINAL/BTYPE byte + LEN + NLEN).
/// * `LEN + 5 == compressed.len()` — exactly one stored sub-block fills the
///   BGZF payload. Every real level-0 producer (`samtools view -u`, htsjdk,
///   [`InlineBgzfCompressor`](crate::writer::InlineBgzfCompressor)) emits
///   this shape. We intentionally do **not** fall back to libdeflater if
///   this check fails: the input is either malformed or uses a multi-sub-
///   block stored stream we have no real-world reason to support, and a
///   loud error beats silently papering over corruption.
/// * `LEN == uncompressed_size` — the deflate frame's LEN agrees with the
///   BGZF footer's ISIZE field.
///
/// NLEN (one's complement of LEN) is not checked. NLEN doesn't cover the
/// payload bytes, so a corrupt NLEN with an intact payload would pass the
/// framing check anyway; the BGZF footer's CRC32 is the authoritative
/// integrity check on the data itself.
///
/// CRC32 verification against the BGZF footer is still performed on the
/// copied bytes (the BGZF spec mandates it, and the bypass is intended to
/// skip the libdeflater call, not the integrity check).
///
/// On any framing or verification failure, `output` is rolled back to its
/// original length.
fn copy_stored_and_verify(
    compressed: &[u8],
    uncompressed_size: usize,
    expected_crc: u32,
    block_len: usize,
    output: &mut Vec<u8>,
) -> io::Result<()> {
    // Deflate framing is 5 bytes (BFINAL/BTYPE byte + LEN + NLEN).
    if compressed.len() < 5 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "BGZF stored block too small: {} bytes (need >= 5 for deflate framing)",
                compressed.len()
            ),
        ));
    }
    let len = usize::from(u16::from_le_bytes([compressed[1], compressed[2]]));
    // For a well-formed level-0 BGZF block, the entire payload is one stored
    // sub-block: 5 framing bytes + LEN payload bytes. If this doesn't hold,
    // the file is malformed (or contains a multi-sub-block stream, which no
    // real producer emits at level 0).
    if len + 5 != compressed.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("BGZF stored block size mismatch: LEN={len}, payload={}", compressed.len()),
        ));
    }
    if len != uncompressed_size {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("BGZF stored block ISIZE mismatch: footer={uncompressed_size}, LEN={len}"),
        ));
    }

    let start = output.len();
    // `len + 5 == compressed.len()` is checked above, so `&compressed[5..]`
    // is exactly the LEN payload bytes.
    output.extend_from_slice(&compressed[5..]);
    let result = verify_decompression(
        &output[start..start + len],
        uncompressed_size,
        expected_crc,
        block_len,
    );
    if result.is_err() {
        output.truncate(start);
    }
    result
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
        let mut data = vec![0x00; BGZF_HEADER_SIZE];
        data[0] = 0x00; // Invalid magic (should be 0x1f)
        data[1] = 0x00; // Invalid magic (should be 0x8b)
        let mut reader = Cursor::new(data);
        let result = read_raw_block(&mut reader);
        assert!(result.is_err());
        assert!(
            result
                .expect_err("should fail with invalid BGZF magic")
                .to_string()
                .contains("Invalid BGZF magic")
        );
    }

    #[test]
    fn test_read_eof() {
        let mut reader = Cursor::new(Vec::<u8>::new());
        let result = read_raw_block(&mut reader).expect("reading raw BGZF block should succeed");
        assert!(result.is_none());
    }

    #[test]
    fn test_decompress_eof_block() {
        let block = RawBgzfBlock { data: BGZF_EOF.to_vec() };
        let mut decompressor = Decompressor::new();
        let result = decompress_block(&block, &mut decompressor)
            .expect("decompressing block should succeed");
        assert!(result.is_empty());
    }

    #[test]
    fn test_decompress_block_into_eof() {
        let block = RawBgzfBlock { data: BGZF_EOF.to_vec() };
        let mut decompressor = Decompressor::new();
        let mut output = Vec::new();
        decompress_block_into(&block, &mut decompressor, &mut output)
            .expect("decompressing block into buffer should succeed");
        assert!(output.is_empty());
    }

    #[test]
    fn test_decompress_block_into_appends() {
        // Create a compressed block using the writer
        use crate::writer::InlineBgzfCompressor;

        let original_data = b"Hello, BGZF world!";
        let mut compressor = InlineBgzfCompressor::new(6);
        compressor.write_all(original_data).expect("writing original data should succeed");
        compressor.flush().expect("flushing compressor should succeed");
        let blocks = compressor.take_blocks();
        assert_eq!(blocks.len(), 1);

        let block = RawBgzfBlock { data: blocks[0].data.clone() };
        let mut decompressor = Decompressor::new();

        // Start with existing data in the buffer
        let mut output = vec![1, 2, 3];
        decompress_block_into(&block, &mut decompressor, &mut output)
            .expect("decompressing block into buffer should succeed");

        // Should have preserved existing data and appended decompressed data
        assert_eq!(&output[0..3], &[1, 2, 3]);
        assert_eq!(&output[3..], original_data);
    }

    #[test]
    fn test_decompress_block_into_equivalence() {
        // Create a compressed block
        use crate::writer::InlineBgzfCompressor;

        let original_data = b"Test data for equivalence check";
        let mut compressor = InlineBgzfCompressor::new(6);
        compressor.write_all(original_data).expect("writing original data should succeed");
        compressor.flush().expect("flushing compressor should succeed");
        let blocks = compressor.take_blocks();

        let block = RawBgzfBlock { data: blocks[0].data.clone() };
        let mut decompressor = Decompressor::new();

        // Decompress using original function
        let result1 = decompress_block(&block, &mut decompressor)
            .expect("decompressing block should succeed");

        // Decompress using new function
        let mut result2 = Vec::new();
        decompress_block_into(&block, &mut decompressor, &mut result2)
            .expect("decompressing block into result2 should succeed");

        // Should produce identical results
        assert_eq!(result1, result2);
        assert_eq!(result1, original_data);
    }

    #[test]
    fn test_decompress_and_verify_truncates_output_on_error() {
        // Create a valid compressed block
        use crate::writer::InlineBgzfCompressor;

        let original_data = b"Test data for truncation check";
        let mut compressor = InlineBgzfCompressor::new(6);
        compressor.write_all(original_data).expect("writing original data should succeed");
        compressor.flush().expect("flushing compressor should succeed");
        let blocks = compressor.take_blocks();
        let block = RawBgzfBlock { data: blocks[0].data.clone() };

        let mut decompressor = Decompressor::new();
        let mut output = vec![1, 2, 3];

        // Call decompress_and_verify with a wrong CRC to trigger failure
        let result = decompress_and_verify(
            block.compressed_data(),
            block.uncompressed_size(),
            block.crc32().wrapping_add(1), // wrong CRC
            block.len(),
            &mut decompressor,
            &mut output,
        );

        assert!(result.is_err());
        // The output buffer should be rolled back to its original length
        assert_eq!(output, vec![1, 2, 3]);
    }

    /// Level-0 BGZF blocks contain a deflate **stored** sub-block — the reader
    /// must round-trip these without invoking libdeflater.
    #[test]
    fn test_decompress_stored_block_roundtrip() {
        use crate::writer::InlineBgzfCompressor;

        // Picked to span more than one BGZF block so we exercise multi-block decode.
        let original: Vec<u8> = (0u32..(100 * 1024)).map(|i| (i % 251) as u8).collect();
        let mut compressor = InlineBgzfCompressor::new(0);
        compressor.write_all(&original).expect("write");
        compressor.flush().expect("flush");
        let blocks = compressor.take_blocks();
        assert!(blocks.len() >= 2, "expected at least two blocks for >64 KiB at level 0");

        // Sanity: each block's payload should look like a stored deflate frame —
        // first byte's BTYPE bits (1-2) are 00.
        for cb in &blocks {
            let block = RawBgzfBlock { data: cb.data.clone() };
            let payload = block.compressed_data();
            assert!(!payload.is_empty(), "stored block payload must be non-empty");
            assert_eq!(payload[0] & 0b110, 0, "stored block should have BTYPE=00");
        }

        let mut decompressor = Decompressor::new();
        let mut out = Vec::with_capacity(original.len());
        for cb in &blocks {
            let block = RawBgzfBlock { data: cb.data.clone() };
            decompress_block_into(&block, &mut decompressor, &mut out).expect("decompress");
        }
        assert_eq!(out, original);
    }

    /// Same round-trip via [`decompress_block_slice_into`], which is on the
    /// hot path for fastq pipelines (avoids materialising a `RawBgzfBlock`).
    #[test]
    fn test_decompress_stored_block_slice_into_roundtrip() {
        use crate::writer::InlineBgzfCompressor;

        let original = b"stored slice payload bytes here";
        let mut compressor = InlineBgzfCompressor::new(0);
        compressor.write_all(original).expect("write");
        compressor.flush().expect("flush");
        let blocks = compressor.take_blocks();
        assert_eq!(blocks.len(), 1);

        let mut decompressor = Decompressor::new();
        let mut out = Vec::new();
        decompress_block_slice_into(&blocks[0].data, &mut decompressor, &mut out)
            .expect("decompress");
        assert_eq!(out.as_slice(), original.as_slice());
    }

    /// A stored block whose `LEN` field disagrees with the BGZF payload size
    /// is malformed; the fast path should reject it (rather than silently
    /// truncating or falling through to libdeflater). This exercises the
    /// `len + 5 != compressed.len()` branch in `copy_stored_and_verify`.
    #[test]
    fn test_decompress_stored_block_rejects_bad_len() {
        use crate::writer::InlineBgzfCompressor;

        let original = b"stored len mismatch case";
        let mut compressor = InlineBgzfCompressor::new(0);
        compressor.write_all(original).expect("write");
        compressor.flush().expect("flush");
        let blocks = compressor.take_blocks();
        let mut data = blocks[0].data.clone();

        // Corrupt LEN inside the deflate frame (compressed payload starts at
        // BGZF_HEADER_SIZE; the LEN field is bytes 1-2 of the deflate frame).
        let len_off = BGZF_HEADER_SIZE + 1;
        data[len_off] = data[len_off].wrapping_add(1);

        let mut decompressor = Decompressor::new();
        let mut out = vec![0xAA, 0xBB];
        let err = decompress_block_slice_into(&data, &mut decompressor, &mut out)
            .expect_err("expected malformed stored block to be rejected");
        assert!(
            err.to_string().contains("stored block"),
            "error should mention stored block: {err}"
        );
        // Output buffer must be rolled back on failure.
        assert_eq!(out, vec![0xAA, 0xBB]);
    }

    /// A stored block whose payload bytes round-trip cleanly through the
    /// framing checks but whose BGZF footer CRC32 is wrong must still be
    /// rejected — the bypass skips libdeflater, not the integrity check.
    /// Exercises `verify_decompression` failure on the stored-block path,
    /// including the output-buffer rollback in `copy_stored_and_verify`.
    #[test]
    fn test_decompress_stored_block_rejects_bad_crc() {
        use crate::writer::InlineBgzfCompressor;

        let original = b"stored crc check payload bytes";
        let mut compressor = InlineBgzfCompressor::new(0);
        compressor.write_all(original).expect("write");
        compressor.flush().expect("flush");
        let blocks = compressor.take_blocks();
        let mut data = blocks[0].data.clone();

        // Corrupt the CRC32 field in the footer (last 8 bytes; CRC32 is the
        // first 4 of those, ISIZE is the last 4). Flipping a single bit is
        // enough to make verification fail.
        let crc_off = data.len() - BGZF_FOOTER_SIZE;
        data[crc_off] ^= 0x01;

        let mut decompressor = Decompressor::new();
        let mut out = vec![0xDE, 0xAD];
        let err = decompress_block_slice_into(&data, &mut decompressor, &mut out)
            .expect_err("expected CRC mismatch to be rejected");
        assert!(err.to_string().contains("CRC32"), "error should mention CRC32: {err}");
        // Output buffer must be rolled back on failure.
        assert_eq!(out, vec![0xDE, 0xAD]);
    }

    /// A stored block whose deflate LEN agrees with the wire payload size
    /// but disagrees with the BGZF footer's ISIZE must be rejected.
    /// Exercises the `len != uncompressed_size` branch in
    /// `copy_stored_and_verify`.
    #[test]
    fn test_decompress_stored_block_rejects_isize_mismatch() {
        use crate::writer::InlineBgzfCompressor;

        let original = b"stored isize mismatch case";
        let mut compressor = InlineBgzfCompressor::new(0);
        compressor.write_all(original).expect("write");
        compressor.flush().expect("flush");
        let blocks = compressor.take_blocks();
        let mut data = blocks[0].data.clone();

        // Corrupt the low byte of ISIZE in the footer (last 4 bytes of the
        // 8-byte footer). LEN inside the deflate frame and the wire payload
        // size remain consistent, so we land in the `len != uncompressed_size`
        // check rather than the framing-size check.
        let isize_off = data.len() - 4;
        data[isize_off] = data[isize_off].wrapping_add(1);

        let mut decompressor = Decompressor::new();
        let mut out = Vec::new();
        let err = decompress_block_slice_into(&data, &mut decompressor, &mut out)
            .expect_err("expected ISIZE mismatch to be rejected");
        assert!(
            err.to_string().contains("ISIZE mismatch"),
            "error should mention ISIZE mismatch: {err}"
        );
        assert!(out.is_empty(), "output should be rolled back on failure");
    }

    /// A BGZF block whose compressed payload is fewer than 5 bytes can't
    /// contain a valid deflate stored frame (BFINAL/BTYPE byte + LEN + NLEN
    /// = 5 bytes minimum). The bypass guards against this. We can't
    /// construct such an input via `InlineBgzfCompressor` (which always
    /// emits valid frames), so we synthesise one byte by byte.
    #[test]
    fn test_decompress_stored_block_rejects_truncated_framing() {
        // Build a BGZF block with a 4-byte compressed payload (zeros — first
        // byte has BTYPE=00 so we enter the stored fast path), an ISIZE of 1
        // in the footer so the wrapper doesn't short-circuit on
        // `uncompressed_size == 0`, and a CRC of 0 (we never reach the CRC
        // check). Layout: 18-byte header + 4-byte payload + 8-byte footer
        // = 30 bytes total.
        const BLOCK_SIZE: usize = BGZF_HEADER_SIZE + 4 + BGZF_FOOTER_SIZE;
        let mut data = vec![0u8; BLOCK_SIZE];
        // gzip magic + deflate method + FEXTRA flag.
        data[0] = 0x1f;
        data[1] = 0x8b;
        data[2] = 0x08;
        data[3] = 0x04;
        // BC subfield ID.
        data[12] = b'B';
        data[13] = b'C';
        // BSIZE = total_block_size - 1, little-endian u16 at offset 16.
        let bsize_bytes = u16::try_from(BLOCK_SIZE - 1).expect("block fits in u16").to_le_bytes();
        data[16] = bsize_bytes[0];
        data[17] = bsize_bytes[1];
        // Payload bytes 18..22 stay zero from the vec init — first byte's
        // BTYPE bits (1-2) are 00 so we enter the stored-block fast path.
        // Footer: CRC32 = 0 (unused — we error before checking it). ISIZE = 1
        // so `decompress_block_into` doesn't skip on `isize == 0`.
        data[BLOCK_SIZE - 4] = 1;

        let mut decompressor = Decompressor::new();
        let mut out = Vec::new();
        let err = decompress_block_slice_into(&data, &mut decompressor, &mut out)
            .expect_err("expected truncated stored framing to be rejected");
        assert!(
            err.to_string().contains("stored block too small"),
            "error should mention truncated stored framing: {err}"
        );
        assert!(out.is_empty(), "output should be rolled back on failure");
    }
}

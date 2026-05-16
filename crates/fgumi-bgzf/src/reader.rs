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

    /// Take ownership of the block's raw byte buffer, leaving the
    /// block with an empty `Vec<u8>`.
    ///
    /// Use this to return the underlying allocation to a recycling
    /// buffer pool once the block has been decompressed.
    #[must_use]
    pub fn take_data(&mut self) -> Vec<u8> {
        std::mem::take(&mut self.data)
    }
}

// ============================================================================
// Reading Functions
// ============================================================================

/// Read a single raw BGZF block from the input, sourcing its byte
/// storage from `get_buf` (only invoked once a block header is
/// successfully read).
///
/// This is the buffer-recycling variant of [`read_raw_block`]: callers
/// can pass a closure that pulls from a `Vec<u8>` recycling pool to
/// avoid the per-block `mi_malloc` that the convenience entry point
/// pays. On clean EOF (`Ok(None)`) the closure is *not* called, so the
/// pool is never charged for a wasted checkout.
///
/// # Errors
///
/// Returns an error if the block header is invalid or reading fails.
fn read_raw_block_into<R, F>(reader: &mut R, get_buf: F) -> io::Result<Option<RawBgzfBlock>>
where
    R: Read + ?Sized,
    F: FnOnce() -> Vec<u8>,
{
    // Read the 18-byte header (or detect clean EOF).
    let mut header = [0u8; BGZF_HEADER_SIZE];
    match reader.read_exact(&mut header) {
        Ok(()) => {}
        Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(None),
        Err(e) => return Err(e),
    }

    // Validate header (same checks as `read_raw_block`).
    if header[0] != 0x1f || header[1] != 0x8b {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "Invalid BGZF magic: expected 0x1f 0x8b, got 0x{:02x} 0x{:02x}",
                header[0], header[1]
            ),
        ));
    }
    if header[2] != 0x08 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Invalid compression method: expected 0x08, got 0x{:02x}", header[2]),
        ));
    }
    if header[3] & 0x04 == 0 {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "BGZF block missing FEXTRA flag"));
    }
    if header[12] != b'B' || header[13] != b'C' {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "Invalid BGZF subfield ID: expected 'BC', got '{}{}'",
                header[12] as char, header[13] as char
            ),
        ));
    }

    let bsize = usize::from(u16::from_le_bytes([header[16], header[17]]));
    let block_size = bsize + 1;
    if block_size < BGZF_HEADER_SIZE + BGZF_FOOTER_SIZE {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("BGZF block too small: {block_size} bytes"),
        ));
    }

    // Acquire a buffer from the caller-supplied source only now —
    // after a valid header has been read, so a checkout never charges
    // for a wasted-on-EOF cycle.
    let mut data = get_buf();
    data.clear();
    data.resize(block_size, 0);
    data[..BGZF_HEADER_SIZE].copy_from_slice(&header);
    reader.read_exact(&mut data[BGZF_HEADER_SIZE..])?;

    Ok(Some(RawBgzfBlock { data }))
}

/// Read a single raw BGZF block from the input, allocating a fresh
/// `Vec<u8>` per block. Convenience wrapper around
/// [`read_raw_block_into`] for callers that do not want to manage a
/// recycling buffer pool.
///
/// Returns `Ok(Some(block))` if a block was read, `Ok(None)` at EOF,
/// or an error if reading failed or the data is invalid.
///
/// # Errors
///
/// Returns an error if the block header is invalid or reading fails.
#[cfg(test)]
fn read_raw_block<R: Read + ?Sized>(reader: &mut R) -> io::Result<Option<RawBgzfBlock>> {
    read_raw_block_into(reader, Vec::new)
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
    // No recycling pool to feed; discard the EOF marker's buffer.
    read_raw_blocks_into(reader, max_blocks, Vec::new, |_| {})
}

/// Read multiple raw BGZF blocks as a batch, sourcing each block's
/// byte storage from `get_buf` and routing any discarded blocks'
/// buffers through `recycle`.
///
/// The buffer-recycling variant of [`read_raw_blocks`]: pass a
/// closure that pulls `Vec<u8>` from a recycling pool to avoid the
/// per-block `mi_malloc` that the convenience entry point pays.
///
/// `get_buf` is invoked only when a block header is successfully
/// read, so clean EOF detection does not waste a checkout.
/// `recycle` is invoked for the BGZF EOF marker block that closes
/// every stream — its buffer is consumed via `get_buf` to read the
/// trailer bytes but is not returned to the caller in the result
/// `Vec`. Without the `recycle` callback the EOF marker would
/// silently lose one pooled buffer per file (a 1-buffer-per-stream
/// pool drain across Phase 2's K spill files).
///
/// # Errors
///
/// Returns an error if reading or validation fails. On a mid-block
/// I/O error the freshly-checked-out buffer for the failing block
/// is dropped (not routed through `recycle`); this is acceptable
/// because the pool's `try_send`-based `checkin` already drops
/// excess buffers, and a mid-block error aborts the whole sort.
pub fn read_raw_blocks_into<R, F, G>(
    reader: &mut R,
    max_blocks: usize,
    mut get_buf: F,
    mut recycle: G,
) -> io::Result<Vec<RawBgzfBlock>>
where
    R: Read + ?Sized,
    F: FnMut() -> Vec<u8>,
    G: FnMut(Vec<u8>),
{
    let mut blocks = Vec::with_capacity(max_blocks);
    for _ in 0..max_blocks {
        match read_raw_block_into(reader, &mut get_buf)? {
            Some(mut block) => {
                if block.is_eof() {
                    // Skip the EOF marker block; return its buffer
                    // through the recycle callback so the pool
                    // does not lose one slot per stream.
                    recycle(block.take_data());
                } else {
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
/// Returns an error if decompression fails.
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

/// Decompress raw deflate data into the output buffer and verify the result.
///
/// This is the shared implementation for both [`decompress_block_into`] and
/// [`decompress_block_slice_into`], consolidating the decompress + resize + verify logic.
fn decompress_and_verify(
    compressed: &[u8],
    uncompressed_size: usize,
    expected_crc: u32,
    block_len: usize,
    decompressor: &mut Decompressor,
    output: &mut Vec<u8>,
) -> io::Result<()> {
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
    fn test_read_raw_block_into_reuses_buffer() {
        // Build a single-block BGZF stream by compressing a payload
        // and concatenating its raw bytes.
        use crate::writer::InlineBgzfCompressor;
        let payload = b"hello, recycled buffer";
        let mut compressor = InlineBgzfCompressor::new(6);
        compressor.write_all(payload).expect("write");
        compressor.flush().expect("flush");
        let blocks = compressor.take_blocks();
        assert_eq!(blocks.len(), 1, "small payload fits in one BGZF block");
        let stream = blocks.into_iter().next().unwrap().data;

        // Pre-warm a buffer with capacity well above the block size
        // so we can prove the buffer is reused, not freshly allocated.
        let recycled = Vec::with_capacity(64 * 1024);
        let recycled_ptr = recycled.as_ptr();
        let recycled_cap = recycled.capacity();

        let mut reader = std::io::Cursor::new(stream);
        let mut buf_holder = Some(recycled);
        let block = read_raw_block_into(&mut reader, || buf_holder.take().expect("buf consumed"))
            .expect("read should succeed")
            .expect("one block before EOF");

        // The returned block's data Vec is exactly the one we supplied:
        // same allocation (ptr equal at offset 0) and capacity preserved.
        assert!(block.data.capacity() >= recycled_cap);
        assert!(std::ptr::eq(block.data.as_ptr(), recycled_ptr));

        // Decompresses back to the original payload.
        let mut decompressor = Decompressor::new();
        let out = decompress_block(&block, &mut decompressor).expect("decompress");
        assert_eq!(out.as_slice(), payload);
    }

    #[test]
    fn test_read_raw_block_into_skips_callback_on_eof() {
        // get_buf must not be invoked when the reader is already at
        // clean EOF (no header bytes available).
        let mut reader = std::io::Cursor::new(Vec::<u8>::new());
        let called = std::cell::Cell::new(false);
        let result = read_raw_block_into(&mut reader, || {
            called.set(true);
            Vec::new()
        })
        .expect("EOF should not be an error");
        assert!(result.is_none(), "clean EOF returns None");
        assert!(!called.get(), "get_buf should not be invoked on EOF");
    }

    #[test]
    fn test_take_data_extracts_owned_buffer() {
        let mut block = RawBgzfBlock { data: vec![1, 2, 3, 4, 5] };
        let data = block.take_data();
        assert_eq!(data, vec![1, 2, 3, 4, 5]);
        assert!(block.data.is_empty(), "block's data should be drained");
    }

    #[test]
    fn test_read_raw_blocks_into_recycles_eof_marker_and_iterates_get_buf() {
        // Build a 3-block stream + EOF marker by compressing three
        // payloads separately. Each `flush` emits the buffered bytes
        // as its own BGZF block, so three flushes produce three blocks
        // even though the payloads are small.
        use crate::writer::InlineBgzfCompressor;
        let payloads: [&[u8]; 3] = [b"alpha", b"beta!", b"gamma!!"];
        let mut stream: Vec<u8> = Vec::new();
        for p in payloads {
            let mut compressor = InlineBgzfCompressor::new(6);
            compressor.write_all(p).expect("write");
            compressor.flush().expect("flush");
            let blocks = compressor.take_blocks();
            assert_eq!(blocks.len(), 1);
            stream.extend_from_slice(&blocks[0].data);
        }
        // Append the BGZF EOF marker so the reader produces an
        // EOF-marker block that `read_raw_blocks_into` must recycle.
        stream.extend_from_slice(&BGZF_EOF);

        let mut reader = std::io::Cursor::new(stream);
        let get_calls = std::cell::Cell::new(0usize);
        let recycled = std::cell::RefCell::new(Vec::<Vec<u8>>::new());

        let blocks = read_raw_blocks_into(
            &mut reader,
            8, // max_blocks
            || {
                get_calls.set(get_calls.get() + 1);
                Vec::with_capacity(64 * 1024)
            },
            |buf| recycled.borrow_mut().push(buf),
        )
        .expect("read");

        // 3 payload blocks returned (EOF marker skipped).
        assert_eq!(blocks.len(), 3, "3 payload blocks returned");
        // `get_buf` was called once per block including the EOF marker = 4 times.
        assert_eq!(get_calls.get(), 4, "get_buf called once per real block");
        // The EOF marker's buffer was routed through `recycle`.
        assert_eq!(recycled.borrow().len(), 1, "EOF marker buffer recycled once");
        assert!(
            recycled.borrow()[0].capacity() >= 64 * 1024,
            "recycled buffer preserves original capacity"
        );

        // Payload bytes round-trip via decompress_block.
        let mut decompressor = Decompressor::new();
        for (block, expected) in blocks.iter().zip(payloads.iter()) {
            let out = decompress_block(block, &mut decompressor).expect("decompress");
            assert_eq!(out.as_slice(), *expected);
        }
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
}

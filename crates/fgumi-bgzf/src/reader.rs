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
///
/// Re-exported from [`crate::header`], which owns the layout this reader frames
/// blocks with, so the constant and the field offsets cannot drift apart.
pub use crate::header::BGZF_HEADER_SIZE;

/// Size of the BGZF block footer (CRC32 + ISIZE).
pub const BGZF_FOOTER_SIZE: usize = 8;

/// Largest uncompressed payload a single BGZF block can hold.
///
/// A block's `BSIZE` is a `u16` holding the total size minus one, which caps a
/// block at 64 KiB; the SAM/BGZF spec and htslib's `BGZF_MAX_BLOCK_SIZE` both
/// use this bound, and the `bgzf` crate refuses to *write* a larger one
/// (`BlockSizeExceeded`). The footer's ISIZE is a `u32`, so a corrupt or
/// hostile block can claim up to 4 GiB; anything a caller sizes from that value
/// has to be bounded first.
///
/// Note this is deliberately 64 KiB rather than [`crate::writer::BGZF_MAX_BLOCK_SIZE`]
/// (65280), which is the size *we* fill blocks to. Other writers legitimately
/// emit up to the spec limit, so validating against our own chunk size would
/// reject valid input.
pub const MAX_UNCOMPRESSED_BLOCK_SIZE: usize = 64 * 1024;

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
/// Returns `Ok(Some(block))` if a block was read, or `Ok(None)` at a clean EOF
/// -- meaning the stream ended *before* the first header byte. Input that ends
/// part-way through a block, header included, is truncated rather than finished
/// and reports `UnexpectedEof`.
///
/// # Errors
///
/// Returns an error if the block header is invalid, the block is truncated, or
/// reading fails.
fn read_raw_block<R: Read + ?Sized>(reader: &mut R) -> io::Result<Option<RawBgzfBlock>> {
    // Read the 18-byte header. Only a stream that ends *before* the first header
    // byte is a clean EOF; a header that starts and then runs out is truncated
    // input and must surface as `UnexpectedEof`, the same contract the block-body
    // length check below enforces. Reading the header with a single `read_exact`
    // would collapse both cases into `Ok(None)` and silently drop a partial block.
    // Splitting the header read costs nothing on the paths that matter: every
    // production caller wraps its input in a `BufReader`, so the one-byte probe is
    // a buffer memcpy rather than a second syscall.
    let mut header = [0u8; BGZF_HEADER_SIZE];
    loop {
        match reader.read(&mut header[..1]) {
            Ok(0) => return Ok(None),
            Ok(_) => break,
            // `read` (unlike `read_exact`) surfaces `Interrupted` to the caller.
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    reader.read_exact(&mut header[1..])?;

    // One predicate for every reader, so a header this accepts is one the input
    // classifier and `noodles_bgzf` accept too. This is stricter than the checks
    // that used to live here, which took `FLG & FEXTRA` and never looked at XLEN
    // or SLEN at all. A `BC`-first header with a trailing subfield (XLEN 10) was
    // framed happily — `BSIZE` is still at 16-17 in that layout — and the four
    // extra bytes then rode into the payload, which every consumer takes as
    // `[18 .. len - 8]`. The failure surfaced downstream as `BGZF stored block
    // size mismatch: LEN=89, payload=16`, a message about stored-block internals
    // that says nothing about the header being the problem.
    // The rejection is carried, not stringified: `HeaderRejection` is an
    // `std::error::Error`, so passing it whole keeps it reachable through
    // `source()` / `downcast_ref` for a caller that wants the failing field
    // programmatically. `io::Error`'s `Display` delegates to it, so the rendered
    // message is unchanged.
    crate::header::validate(&header)
        .map_err(|rejection| io::Error::new(io::ErrorKind::InvalidData, rejection))?;

    // BSIZE is the total block size minus one, and must be big enough to describe
    // a header plus a footer — see `block_size_checked`, which is the one
    // definition of that floor.
    let Some(block_size) = crate::header::block_size_checked(&header) else {
        let stored = crate::header::block_size(&header)
            .expect("a validated header is long enough to carry BSIZE");
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("BGZF block too small: {stored} bytes"),
        ));
    };

    // Reserve exactly what the block needs and fill it, rather than
    // `vec![0u8; block_size]`, which memsets up to 64 KiB that is immediately
    // overwritten.
    let mut data = Vec::with_capacity(block_size);
    data.extend_from_slice(&header);

    // Read remaining block data. `read_to_end` on a bounded `take` avoids
    // pre-zeroing, but unlike `read_exact` it returns `Ok` on early EOF, so the
    // length must be checked explicitly -- without this a truncated BAM would
    // silently short-read instead of erroring.
    let remaining = block_size - BGZF_HEADER_SIZE;
    reader.take(remaining as u64).read_to_end(&mut data)?;
    if data.len() != block_size {
        return Err(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            format!("truncated BGZF block: expected {block_size} bytes, got {}", data.len()),
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

/// The uncompressed size a BGZF block's footer claims, validated against the
/// per-block maximum.
///
/// This is the size a caller must give [`decompress_into_slice`]'s `out`, and
/// it is the number a caller allocates from, so it is checked rather than
/// returned raw: ISIZE is a `u32` sitting in the file, and a corrupt or hostile
/// footer claiming 4 GiB would otherwise become a 4 GiB allocation in the
/// caller before this crate ever saw the block.
///
/// Takes the block as a `&[u8]` so reading the footer costs nothing —
/// [`RawBgzfBlock::uncompressed_size`] answers the same question but needs an
/// owned `Vec`, which would mean copying the block to size a slot for it.
///
/// # Errors
///
/// Returns `io::ErrorKind::InvalidData` if `block` is too short to hold a BGZF
/// header + footer, or if the claimed size exceeds
/// [`MAX_UNCOMPRESSED_BLOCK_SIZE`].
pub fn uncompressed_size(block: &[u8]) -> io::Result<usize> {
    if block.len() < BGZF_HEADER_SIZE + BGZF_FOOTER_SIZE {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "BGZF block too short to contain header + footer",
        ));
    }
    let claimed = uncompressed_size_from_slice(block);
    if claimed > MAX_UNCOMPRESSED_BLOCK_SIZE {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "BGZF block claims an uncompressed size of {claimed} bytes, above the \
                 {MAX_UNCOMPRESSED_BLOCK_SIZE}-byte maximum for a single block"
            ),
        ));
    }
    Ok(claimed)
}

/// Whether a BGZF block's deflate payload is a **stored** (uncompressed) frame,
/// i.e. `BTYPE = 00`.
///
/// `BTYPE` lives in bits 1-2 of the deflate stream's first byte, so masking with
/// `0b110` isolates it. An empty payload is not a stored frame — there is no
/// first byte to read, and the framing parser needs five.
///
/// Shared by every decompress entry point so the dispatch predicate cannot
/// drift between them; [`parse_stored_frame`] validates the rest of the framing
/// once this returns `true`.
#[must_use]
fn is_stored_block(compressed: &[u8]) -> bool {
    !compressed.is_empty() && compressed[0] & 0b110 == 0
}

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

/// Decompress a full BGZF block's DEFLATE payload into a caller-provided,
/// pre-sized slice. `out.len()` must equal the block's ISIZE (uncompressed
/// size), which callers should take from [`uncompressed_size`] — it reads the
/// footer straight off the same `&[u8]` and bounds the claim, so the slot a
/// caller allocates can never come from an unvalidated `u32`.
///
/// This is the fixed-slice analogue of [`decompress_block_slice_into`] (which
/// appends to a `Vec`); it lets a caller decompress straight into an arena
/// slot rather than into a buffer it then has to copy out of. The two names are
/// close and their `slice` refers to opposite operands: in
/// [`decompress_block_slice_into`] it is the *input* block, given as a slice
/// rather than a [`RawBgzfBlock`]; here it is the *output*.
///
/// Like the sibling decompressors ([`decompress_block_into`] /
/// [`decompress_block_slice_into`]), the decompressed payload is verified
/// against the BGZF footer's ISIZE (it must exactly fill `out`) and CRC32, so
/// a short fill or a silently-corrupt-but-decodable block is caught here
/// rather than fed into the arena.
///
/// # Returns
///
/// The number of bytes written, which on success is always `out.len()`: the
/// decompressing paths are held to it by the exact-fill check, and a block with
/// an ISIZE of zero (the BGZF EOF marker) returns `0` against the zero-length
/// slot its ISIZE requires. It is returned so the call reads like the
/// `Read`-style APIs it sits beside, not because a short result is reachable.
///
/// # Errors
///
/// Returns an `io::Error` if `block` is shorter than a BGZF header + footer, if
/// the footer's ISIZE exceeds [`MAX_UNCOMPRESSED_BLOCK_SIZE`], if `out` is not
/// sized to that ISIZE, if the DEFLATE stream is invalid, if it does not
/// exactly fill `out`, or if the CRC32 does not match the footer.
///
/// Note this validates the block's *framing*, not its full header: a payload
/// whose header [`crate::header::validate`] would reject can still reach the
/// decompressor and fail there. Callers reading from a stream get the header
/// check from [`read_raw_blocks`]; this entry point is for callers that already
/// hold a framed block.
///
/// On error `out` is left clobbered — unlike the `Vec` siblings, which roll
/// their output back. There is nothing to roll back to here: the buffer belongs
/// to the caller, who must treat its contents as undefined unless this returns
/// `Ok`.
pub fn decompress_into_slice(
    block: &[u8],
    decompressor: &mut Decompressor,
    out: &mut [u8],
) -> io::Result<usize> {
    // Same accessor the caller sizes `out` with, so the two cannot disagree
    // about either the value or the bound. It carries the length and ISIZE
    // checks, which is why neither is repeated here.
    let uncompressed_size = uncompressed_size(block)?;
    // Check the caller's slot against the footer up front rather than letting a
    // mis-sized `out` surface downstream. Both paths below compare against
    // `out.len()` on the assumption it *is* the ISIZE, so without this a wrongly
    // sized slot is reported as a fault in the block: the stored path would say
    // `ISIZE mismatch: footer=<out.len()>`, naming a value the footer does not
    // contain, and send a reader looking for file corruption instead of at the
    // arena that sized the slot.
    if out.len() != uncompressed_size {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "BGZF output slice is {} bytes, but the block's ISIZE is {uncompressed_size}",
                out.len()
            ),
        ));
    }
    // A block that carries nothing — the BGZF EOF marker is the one every
    // reader meets — is done, and both `Vec` siblings short-circuit here too.
    // libdeflater does return `Ok(0)` for the EOF marker's payload against a
    // zero-length slice, so this is not repairing a failure; it is making the
    // answer ours rather than resting on an undocumented edge of the C library.
    //
    // Placed *after* the size check on purpose: ahead of it, a caller passing a
    // wrongly-sized slot for a zero-ISIZE block would get a silent `Ok(0)`
    // instead of being told the slot is wrong.
    if uncompressed_size == 0 {
        return Ok(0);
    }
    let compressed = compressed_data_from_slice(block);
    // Stored-block fast path — level-0 blocks (`samtools view -u`, htsjdk's
    // level-0 writer, [`InlineBgzfCompressor::new(0)`]) skip the libdeflater
    // round-trip and get the stored-framing-specific LEN/ISIZE diagnostics.
    if is_stored_block(compressed) {
        return copy_stored_and_verify_slice(compressed, out, crc32_from_slice(block), block.len());
    }
    deflate_into_slice_and_verify(
        compressed,
        crc32_from_slice(block),
        block.len(),
        decompressor,
        out,
    )
}

/// Inflate `compressed` into the whole of `out` and verify the result against
/// the BGZF footer, returning the number of bytes written.
///
/// `out.len()` is taken as the expected uncompressed size, so
/// [`verify_decompression`] checks both the exact fill
/// (`bytes_written == out.len()`) and the CRC32. Shared by
/// [`decompress_into_slice`] and [`decompress_and_verify`]'s non-stored branch
/// so the inflate-then-verify invariant lives in one place — the same reason
/// [`parse_stored_frame`] exists for the stored branch.
fn deflate_into_slice_and_verify(
    compressed: &[u8],
    expected_crc: u32,
    block_len: usize,
    decompressor: &mut Decompressor,
    out: &mut [u8],
) -> io::Result<usize> {
    let bytes_written = decompressor.deflate_decompress(compressed, out).map_err(|e| {
        io::Error::new(io::ErrorKind::InvalidData, format!("BGZF decompression failed: {e:?}"))
    })?;
    verify_decompression(&out[..bytes_written], out.len(), expected_crc, block_len)?;
    Ok(bytes_written)
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
    // `payload_len == LEN + 5` is the structural guarantee that there's exactly
    // one stored sub-block spanning the BGZF payload — the form every real
    // level-0 producer emits.
    if is_stored_block(compressed) {
        return copy_stored_and_verify(
            compressed,
            uncompressed_size,
            expected_crc,
            block_len,
            output,
        );
    }

    let start = output.len();
    // Sizing the tail to `uncompressed_size` is what lets the shared helper take
    // `out.len()` as the expected size, so the exact-fill check it performs is
    // the same one this function used to do inline.
    output.resize(start + uncompressed_size, 0);

    let result = deflate_into_slice_and_verify(
        compressed,
        expected_crc,
        block_len,
        decompressor,
        &mut output[start..],
    );

    if result.is_err() {
        output.truncate(start);
    }
    result.map(|_| ())
}

/// Parse and validate a deflate **stored** sub-block frame, returning the
/// LEN-byte payload slice (`&compressed[5..]`). Shared by
/// [`copy_stored_and_verify`] and [`copy_stored_and_verify_slice`] so the
/// stored-framing checks cannot drift between the `Vec` and fixed-slice
/// decompress entry points.
///
/// The caller is responsible for confirming that `compressed[0]` has
/// `BTYPE = 00` before calling. This validates the rest of the framing:
///
/// * `compressed.len() >= 5` (room for the BFINAL/BTYPE byte + LEN + NLEN).
/// * `LEN + 5 == compressed.len()` — exactly one stored sub-block fills the
///   BGZF payload. Every real level-0 producer (`samtools view -u`, htsjdk,
///   [`InlineBgzfCompressor`](crate::writer::InlineBgzfCompressor)) emits
///   this shape. We intentionally do **not** fall back to libdeflater if
///   this check fails: the input is either malformed or uses a multi-sub-
///   block stored stream we have no real-world reason to support, and a
///   loud error beats silently papering over corruption.
/// * `LEN == expected_len` — the deflate frame's LEN agrees with the BGZF
///   footer's ISIZE (the caller passes the footer ISIZE, or the caller-sized
///   output length that equals it).
///
/// NLEN (one's complement of LEN) is not checked. NLEN doesn't cover the
/// payload bytes, so a corrupt NLEN with an intact payload would pass the
/// framing check anyway; the BGZF footer's CRC32 is the authoritative
/// integrity check on the data itself.
///
/// # Errors
///
/// Returns `io::ErrorKind::InvalidData` if the stored framing is truncated,
/// spans more than one sub-block, or its LEN disagrees with `expected_len`.
fn parse_stored_frame(compressed: &[u8], expected_len: usize) -> io::Result<&[u8]> {
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
    if len != expected_len {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("BGZF stored block ISIZE mismatch: footer={expected_len}, LEN={len}"),
        ));
    }
    // `len + 5 == compressed.len()`, so this is exactly the LEN payload bytes.
    Ok(&compressed[5..])
}

/// Decompress a deflate **stored** sub-block by copying its payload directly
/// into `output`, skipping libdeflater. The stored framing is validated by
/// [`parse_stored_frame`]; CRC32 verification against the BGZF footer is still
/// performed on the copied bytes (the BGZF spec mandates it, and the bypass is
/// intended to skip the libdeflater call, not the integrity check).
///
/// The caller must have confirmed `compressed[0]` has `BTYPE = 00` first.
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
    let payload = parse_stored_frame(compressed, uncompressed_size)?;
    let start = output.len();
    output.extend_from_slice(payload);
    let result = verify_decompression(&output[start..], uncompressed_size, expected_crc, block_len);
    if result.is_err() {
        output.truncate(start);
    }
    result
}

/// Slice-writing sibling of [`copy_stored_and_verify`]: copy a deflate
/// **stored** sub-block's payload straight into a caller-sized `out` slice
/// (whose length is the block's ISIZE), skipping libdeflater, then verify the
/// BGZF footer CRC32 over the copied bytes. Returns the number of bytes written
/// (`== out.len()` on success). Used by [`decompress_into_slice`]'s fast path.
///
/// The caller must have confirmed `compressed[0]` has `BTYPE = 00` first.
/// Unlike the `Vec` sibling there is nothing to roll back: `out` is the
/// caller's buffer, and — matching [`decompress_into_slice`]'s non-stored path
/// — a verification failure leaves `out` holding the (CRC-rejected) bytes,
/// which the caller discards along with the returned error.
fn copy_stored_and_verify_slice(
    compressed: &[u8],
    out: &mut [u8],
    expected_crc: u32,
    block_len: usize,
) -> io::Result<usize> {
    // `out` is caller-sized to the footer's ISIZE, so passing `out.len()` also
    // enforces LEN == ISIZE; the returned payload is then exactly `out.len()`.
    let payload = parse_stored_frame(compressed, out.len())?;
    out.copy_from_slice(payload);
    verify_decompression(out, out.len(), expected_crc, block_len)?;
    Ok(out.len())
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
    use std::io::Cursor;

    /// Compress `payload` into a single BGZF block and return its raw bytes.
    fn single_block_bytes(payload: &[u8]) -> Vec<u8> {
        let mut compressor = crate::writer::InlineBgzfCompressor::new(1);
        compressor.write_all(payload).expect("buffer payload");
        compressor.flush().expect("flush to a block");
        let mut bytes = Vec::new();
        compressor.write_blocks_to(&mut bytes).expect("emit block bytes");
        bytes
    }

    /// A `BC`-first header with a trailing extra subfield must be rejected at the
    /// header, not carried into the payload.
    ///
    /// RFC 1952 permits the second subfield, and `BSIZE` is still at bytes 16-17,
    /// so the block used to *frame* correctly — the four extra bytes then rode
    /// into the deflate payload, which this reader takes as `[18 .. len - 8]`,
    /// and the run died downstream on `BGZF stored block size mismatch`, a
    /// message about stored-block internals that named nothing about the header.
    #[test]
    fn test_read_raw_block_rejects_a_trailing_extra_subfield() {
        let full = single_block_bytes(b"the quick brown fox");
        let mut block = full[..10].to_vec(); // magic, CM, FLG, MTIME, XFL, OS
        block.extend_from_slice(&10u16.to_le_bytes()); // XLEN = 10, was 6
        block.extend_from_slice(&full[12..18]); // BC, SLEN, BSIZE — still first
        block.extend_from_slice(b"XY"); // a foreign subfield behind `BC`
        block.extend_from_slice(&0u16.to_le_bytes()); // ...of zero length
        block.extend_from_slice(&full[18..]); // payload and footer

        let err = read_raw_block(&mut Cursor::new(block))
            .expect_err("a header whose extra field is not exactly `BC` must be rejected");
        assert_eq!(err.kind(), io::ErrorKind::InvalidData, "unexpected kind: {err}");
        assert!(
            err.to_string().contains("extra-field length"),
            "the error must name the offending field, got: {err}"
        );
    }

    /// `read_raw_block` reserves capacity and fills with `read_to_end` rather
    /// than `vec![0u8; block_size]`, avoiding a memset of up to 64 KiB that is
    /// immediately overwritten. Unlike `read_exact`, `read_to_end` returns `Ok`
    /// on early EOF, so the explicit length check is what keeps a truncated BAM
    /// from silently short-reading. Without it this test would see `Ok(Some(..))`
    /// carrying a short block.
    #[test]
    fn test_read_raw_block_rejects_a_truncated_block() {
        // Build a well-formed block, then cut bytes off the end.
        let full = single_block_bytes(b"the quick brown fox jumps over the lazy dog");
        assert!(full.len() > BGZF_HEADER_SIZE + BGZF_FOOTER_SIZE, "need a real block");

        // Full block round-trips.
        let mut cursor = Cursor::new(full.clone());
        let block = read_raw_block(&mut cursor).expect("full block should read").expect("a block");
        let full_len = block.data.len();

        // Dropping even one byte must be an error, not a short block.
        let truncated = full[..full_len - 1].to_vec();
        let mut cursor = Cursor::new(truncated);
        let err = read_raw_block(&mut cursor)
            .expect_err("a truncated block must error rather than short-read");
        assert_eq!(
            err.kind(),
            io::ErrorKind::UnexpectedEof,
            "expected UnexpectedEof for a truncated block, got: {err}",
        );
    }

    /// A stream that ends before the first header byte is a clean EOF -- the
    /// only case that may report `Ok(None)`.
    #[test]
    fn test_read_raw_block_reports_eof_on_an_empty_stream() {
        let mut cursor = Cursor::new(Vec::new());
        assert!(
            read_raw_block(&mut cursor).expect("empty stream is a clean EOF").is_none(),
            "an empty stream must read as EOF, not a block",
        );
    }

    /// A header that starts and then runs out is truncated input, not EOF. Every
    /// non-empty prefix shorter than a full header must surface `UnexpectedEof`,
    /// the same contract the block-body length check enforces. The value list is
    /// exhaustive over `1..BGZF_HEADER_SIZE`; the assertion below fails loudly if
    /// the constant ever changes out from under it.
    #[rstest]
    fn test_read_raw_block_rejects_a_partial_header(
        #[values(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17)] prefix_len: usize,
    ) {
        assert_eq!(BGZF_HEADER_SIZE, 18, "prefix_len values must cover 1..BGZF_HEADER_SIZE");

        let full = single_block_bytes(b"payload behind a header that gets cut short");
        let mut cursor = Cursor::new(full[..prefix_len].to_vec());
        let err = read_raw_block(&mut cursor)
            .expect_err("a partial header must error rather than read as EOF");
        assert_eq!(
            err.kind(),
            io::ErrorKind::UnexpectedEof,
            "expected UnexpectedEof for a {prefix_len}-byte header prefix, got: {err}",
        );
    }

    /// The block bytes must be byte-identical to what the pre-optimization
    /// `vec![0u8; n]` + `read_exact` path produced, header included.
    #[test]
    fn test_read_raw_block_bytes_match_a_prezeroed_read() {
        let full = single_block_bytes(b"payload bytes for comparison");

        let mut cursor = Cursor::new(full.clone());
        let block = read_raw_block(&mut cursor).expect("read").expect("a block");

        // Independent oracle: the block is a prefix of the stream of exactly
        // its own length, so compare against the raw bytes directly.
        assert_eq!(
            block.data.as_slice(),
            &full[..block.data.len()],
            "block bytes must match the stream verbatim",
        );
    }

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
    /// `len + 5 != compressed.len()` branch in `parse_stored_frame`.
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
    /// Exercises the `len != expected_len` branch in `parse_stored_frame`.
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

    // ── `uncompressed_size` ─────────────────────────────────────────────────

    /// The public accessor callers size their slot with must agree with the
    /// footer for a real block, and must reject the two inputs that would
    /// otherwise become a bad allocation: a block too short to have a footer,
    /// and a footer claiming more than a block can hold.
    #[rstest]
    #[case::level_0_stored(0)]
    #[case::level_6_deflate(6)]
    fn uncompressed_size_reports_the_footer(#[case] level: u32) {
        let payload = one_block_payload();
        let block = first_block_at_level(&payload, level);
        assert_eq!(uncompressed_size(&block).expect("valid block"), payload.len());
    }

    /// A claim above the per-block maximum is refused rather than handed back
    /// for a caller to allocate from. ISIZE is a `u32`, so the untrusted range
    /// runs to 4 GiB.
    #[test]
    fn uncompressed_size_rejects_a_claim_above_the_block_maximum() {
        let mut block = first_block_at_level(&one_block_payload(), 6);
        let len = block.len();
        let claimed = u32::try_from(MAX_UNCOMPRESSED_BLOCK_SIZE + 1).expect("fits");
        block[len - 4..].copy_from_slice(&claimed.to_le_bytes());
        let err = uncompressed_size(&block).expect_err("an oversized claim must be rejected");
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert!(err.to_string().contains("above the"), "got: {err}");
    }

    /// Exactly the per-block maximum is legal — other writers may fill a block
    /// to the spec limit even though this crate stops at `BGZF_MAX_BLOCK_SIZE`,
    /// so the bound must not be off by one against them.
    #[test]
    fn uncompressed_size_accepts_the_block_maximum() {
        let mut block = first_block_at_level(&one_block_payload(), 6);
        let len = block.len();
        let claimed = u32::try_from(MAX_UNCOMPRESSED_BLOCK_SIZE).expect("fits");
        block[len - 4..].copy_from_slice(&claimed.to_le_bytes());
        assert_eq!(
            uncompressed_size(&block).expect("the maximum itself is valid"),
            MAX_UNCOMPRESSED_BLOCK_SIZE
        );
    }

    /// Too short to hold a footer at all — must error rather than read the
    /// ISIZE out of whatever bytes happen to be there.
    #[test]
    fn uncompressed_size_rejects_a_block_without_a_footer() {
        let err = uncompressed_size(&[0u8; 10]).expect_err("a short block must be rejected");
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert!(err.to_string().contains("too short to contain"), "got: {err}");
    }

    // ── `decompress_into_slice` ─────────────────────────────────────────────

    /// Compress `payload` into BGZF blocks at `level` and return the first
    /// block's raw bytes. Level 0 emits deflate **stored** blocks (the fast
    /// path); any other level emits a real deflate stream.
    fn first_block_at_level(payload: &[u8], level: u32) -> Vec<u8> {
        let mut compressor = crate::writer::InlineBgzfCompressor::new(level);
        compressor.write_all(payload).expect("write payload");
        compressor.flush().expect("flush");
        compressor.take_blocks().into_iter().next().expect("payload fits in one block").data
    }

    /// A payload that fits comfortably in one BGZF block (< 64 KiB).
    fn one_block_payload() -> Vec<u8> {
        b"the quick brown fox jumps over the lazy dog".repeat(100)
    }

    /// Both decompress paths must fill a pre-sized `&mut [u8]` slot with
    /// exactly the original bytes: level 0 through the stored-block copy, any
    /// other level through libdeflater. `expect_stored` pins which path the
    /// case actually takes, so a producer change that stops emitting stored
    /// blocks fails here rather than silently reducing coverage to one path.
    #[rstest]
    #[case::level_0_stored(0, true)]
    #[case::level_6_deflate(6, false)]
    fn decompress_into_slice_fills_presized_slot(#[case] level: u32, #[case] expect_stored: bool) {
        let payload = one_block_payload();
        let block = first_block_at_level(&payload, level);

        let compressed = compressed_data_from_slice(&block);
        assert_eq!(
            is_stored_block(compressed),
            expect_stored,
            "level {level} took the unexpected decompress path"
        );

        let mut out = vec![0u8; uncompressed_size_from_slice(&block)];
        let n = decompress_into_slice(&block, &mut Decompressor::new(), &mut out)
            .expect("decompress into slot");
        assert_eq!(n, payload.len());
        assert_eq!(&out[..n], payload.as_slice());
    }

    /// The BGZF EOF marker carries no payload, and every stream ends with one,
    /// so a caller decompressing each block of a stream in turn hits it. It has
    /// to succeed against a zero-length slot, as the two `Vec` siblings do.
    ///
    /// Its payload is `03 00` — a fixed-Huffman empty frame — so
    /// `is_stored_block` is false and it would otherwise reach libdeflater with
    /// a zero-length output slice. That returns `Ok(0)` today; this pins the
    /// result to the crate rather than to that behaviour.
    #[test]
    fn decompress_into_slice_accepts_the_eof_marker() {
        let compressed = compressed_data_from_slice(&BGZF_EOF);
        assert!(
            !is_stored_block(compressed),
            "the EOF payload is a fixed-Huffman frame, so the stored path must not claim it"
        );
        assert_eq!(uncompressed_size(&BGZF_EOF).expect("the EOF marker is a valid block"), 0);

        let n = decompress_into_slice(&BGZF_EOF, &mut Decompressor::new(), &mut [])
            .expect("the EOF marker must decompress to nothing");
        assert_eq!(n, 0);
    }

    /// A zero-ISIZE block must still reject a wrongly-sized slot. The obvious
    /// place to short-circuit on `uncompressed_size == 0` is *before* the
    /// slot-size check, which would turn this into a silent `Ok(0)` and hide the
    /// caller's sizing bug — the same misreport the up-front check exists to
    /// prevent.
    #[test]
    fn decompress_into_slice_rejects_a_sized_slot_for_an_empty_block() {
        let mut out = [0u8; 16];
        let err = decompress_into_slice(&BGZF_EOF, &mut Decompressor::new(), &mut out)
            .expect_err("a 16-byte slot for a 0-byte block must be rejected");
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert!(err.to_string().contains("output slice is"), "got: {err}");
    }

    /// A block at `level` with a flipped footer CRC32 bit: it still produces
    /// the right number of bytes, but the checksum no longer matches. Run at
    /// both levels because the CRC is verified at two separate call sites —
    /// `copy_stored_and_verify_slice` for level 0 and
    /// `deflate_into_slice_and_verify` for the rest — and deleting either one
    /// used to leave the whole suite green.
    fn block_bad_crc(level: u32) -> (Vec<u8>, Vec<u8>) {
        let mut block = first_block_at_level(&one_block_payload(), level);
        let len = block.len();
        block[len - 8] ^= 0x01; // footer CRC32 is bytes [len - 8 .. len - 4]
        let out = vec![0u8; uncompressed_size_from_slice(&block)];
        (block, out)
    }

    /// A level-6 block whose footer ISIZE is 16 too high, with `out` sized to
    /// that inflated ISIZE. libdeflater reports the true (shorter) length via
    /// `Ok(n)` rather than erroring, so this is what reaches the exact-fill
    /// check — the guarantee the docs single out, and the one branch that
    /// stops the arena being handed a slot with a stale tail.
    fn block_short_fill() -> (Vec<u8>, Vec<u8>) {
        let mut block = first_block_at_level(&one_block_payload(), 6);
        let len = block.len();
        let inflated = uncompressed_size_from_slice(&block) + 16;
        block[len - 4..].copy_from_slice(&u32::try_from(inflated).expect("fits").to_le_bytes());
        (block, vec![0u8; inflated])
    }

    /// A footer claiming more than a single BGZF block can hold. ISIZE is a
    /// `u32`, so an unbounded contract would have a caller allocate up to 4 GiB
    /// from a corrupt footer before anything validated it.
    fn block_isize_above_max() -> (Vec<u8>, Vec<u8>) {
        let mut block = first_block_at_level(&one_block_payload(), 6);
        let len = block.len();
        let claimed = u32::try_from(MAX_UNCOMPRESSED_BLOCK_SIZE + 1).expect("fits");
        block[len - 4..].copy_from_slice(&claimed.to_le_bytes());
        // Deliberately a *correctly* sized slot for the claim: the bound must
        // fire on the claim itself, not merely because `out` disagrees with it.
        (block, vec![0u8; MAX_UNCOMPRESSED_BLOCK_SIZE + 1])
    }

    /// A caller-oversized `out` (ISIZE + 16). The error must name the slice as
    /// the thing that is wrong: before the up-front size check this surfaced as
    /// a complaint about the block's own ISIZE, sending a reader after file
    /// corruption that isn't there. The level is irrelevant — this check fires
    /// before the stored/deflate dispatch — so one case covers it.
    fn block_oversized_out() -> (Vec<u8>, Vec<u8>) {
        let block = first_block_at_level(&one_block_payload(), 6);
        let out = vec![0u8; uncompressed_size_from_slice(&block) + 16];
        (block, out)
    }

    /// Synthesised stored block whose 4-byte payload can't hold the 5-byte
    /// deflate stored frame (BFINAL/BTYPE + LEN + NLEN), so the stored fast
    /// path errors before copying. `InlineBgzfCompressor` only emits valid
    /// frames, so this is built byte by byte: 18-byte header + 4-byte payload
    /// + 8-byte footer.
    ///
    /// The header is a **well-formed** one — `XLEN = 6` and the `BC` subfield's
    /// `SLEN = 2` are set, not left zero — so the fixture isolates the stored
    /// framing defect it is named for. An invalid header here would be rejected
    /// for a different reason by anything that validates one, making the case
    /// pass for the wrong cause.
    fn block_truncated_stored_framing() -> (Vec<u8>, Vec<u8>) {
        const BLOCK_SIZE: usize = BGZF_HEADER_SIZE + 4 + BGZF_FOOTER_SIZE;
        let mut data = vec![0u8; BLOCK_SIZE];
        data[0] = 0x1f; // gzip magic + deflate method + FEXTRA flag
        data[1] = 0x8b;
        data[2] = 0x08;
        data[3] = 0x04;
        data[10] = 0x06; // XLEN = 6: the extra field is exactly the BC subfield
        data[11] = 0x00;
        data[12] = b'B'; // BC subfield ID
        data[13] = b'C';
        data[14] = 0x02; // SLEN = 2: BC holds a two-byte BSIZE
        data[15] = 0x00;
        let bsize_bytes = u16::try_from(BLOCK_SIZE - 1).expect("block fits in u16").to_le_bytes();
        data[16] = bsize_bytes[0];
        data[17] = bsize_bytes[1];
        // Payload bytes 18..22 stay zero → BTYPE bits (1-2) are 00, the stored
        // fast path. Footer ISIZE = 1 matches the caller-sized `out`; the CRC
        // is never reached.
        data[BLOCK_SIZE - 4] = 1;
        debug_assert!(crate::header::is_bgzf_header(&data), "fixture header must be well-formed");
        (data, vec![0u8; 1])
    }

    /// A block shorter than the 26-byte header + footer minimum — must error
    /// rather than panic on an out-of-bounds slice or a subtract overflow.
    fn block_too_short() -> (Vec<u8>, Vec<u8>) {
        (vec![0u8; 10], vec![0u8; 16])
    }

    /// Every malformed input is rejected as `InvalidData` with a message that
    /// names what went wrong, so a corrupt block is diagnosable rather than a
    /// bare "decompression failed".
    #[rstest]
    #[case::bad_crc_deflate(block_bad_crc(6), "CRC32")]
    #[case::bad_crc_stored(block_bad_crc(0), "CRC32")]
    #[case::short_fill(block_short_fill(), "size mismatch")]
    #[case::isize_above_max(block_isize_above_max(), "above the")]
    #[case::oversized_out(block_oversized_out(), "output slice is")]
    #[case::truncated_stored_framing(block_truncated_stored_framing(), "stored block too small")]
    #[case::too_short_block(block_too_short(), "too short to contain")]
    fn decompress_into_slice_rejects_invalid(
        #[case] block_and_out: (Vec<u8>, Vec<u8>),
        #[case] expect_substr: &str,
    ) {
        let (block, mut out) = block_and_out;
        let err = decompress_into_slice(&block, &mut Decompressor::new(), &mut out)
            .expect_err("malformed block must be rejected");
        assert_eq!(err.kind(), io::ErrorKind::InvalidData, "expected InvalidData, got {err:?}");
        assert!(
            err.to_string().contains(expect_substr),
            "error should contain {expect_substr:?}, got: {err}"
        );
    }
}

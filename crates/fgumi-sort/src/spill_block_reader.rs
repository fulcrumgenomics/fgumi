//! Codec-aware streaming block decompressor for the typed-step
//! `SortSpillDecompress` pipeline step.
//!
//! A spill chunk is either a BGZF block-stream (self-framed `1f 8b` blocks) or
//! a zstd "ZSP1" stream (`[u32 LE frame-len][zstd frame]` records after the
//! 4-byte file magic). Both standalone `fgumi sort` (via the worker pool) and
//! the fused streaming sort write spill chunks with whichever
//! [`SpillCodec`] the sorter is configured for; this
//! decompressor lets the streaming `SortSpillDecompress` step read either.
//!
//! Each decompressed BGZF block (or zstd frame) is returned as one `Vec<u8>`.
//! The block/frame boundaries are immaterial: the downstream `MergeDriver`
//! parses `[len][record]` records out of a slot's decompressed-block queue and
//! reassembles records that span block boundaries, so any chunking of the
//! decompressed byte stream is correct.

use std::io::{self, Read};

use fgumi_bam_io::pread::{RawFrame, SliceLease};
use fgumi_bgzf::reader::{BgzfSliceFramer, SliceFrame};
use libdeflater::Decompressor as BgzfDecompressor;
use zstd::bulk::Decompressor as ZstdDecompressor;

use crate::codec::SpillCodec;
use crate::worker_pool::read_length_prefix;

/// Hard cap on the `u32 LE` length prefix of any zstd spill frame, shared by
/// every spill reader. Frames are produced one per ~64 KiB of input; even
/// pathological expansion can't reach this. Beyond it, we treat the value as
/// corruption rather than allocate gigabytes.
pub(crate) const MAX_ZSTD_FRAME_BYTES: usize = 2 * 1024 * 1024;

/// Size of a zstd spill frame's `u32 LE` length prefix.
const ZSTD_LEN_PREFIX: usize = 4;

/// One parsed, not yet decompressed spill frame and its per-slot sequence
/// number (dense, in file order).
#[derive(Debug)]
pub struct RawBlock {
    /// Per-slot sequence number of the frame.
    pub seq: u64,
    /// The frame's bytes: a BGZF block, or a zstd frame body without its
    /// length prefix.
    pub frame: RawFrame,
    /// The resident bytes the slot's stash stops holding when this frame is
    /// claimed: an owned frame's allocation; for a borrowed frame, its whole
    /// read slice if it is the slice's last frame in the stash, else nothing
    /// (a slice is resident while any of its frames is stashed).
    pub charge: u64,
}

/// Incremental frame parser over a spill file's read slices.
///
/// Cuts every complete frame out of `carry ++ slice`: a frame wholly inside
/// the slice is [`RawFrame::Borrowed`] from its lease (zero copy), the frame
/// that completes the carry is [`RawFrame::Owned`], and the trailing partial
/// frame is copied into the carry (at most one frame), so no slice is pinned
/// by the next slice's frames. BGZF uses [`BgzfSliceFramer`] (EOF-marker
/// blocks are skipped, as `read_raw` skips them); zstd parses
/// `[u32 LE len][frame]` records, rejecting a length over
/// `MAX_ZSTD_FRAME_BYTES` with the same error the sequential reader raises.
#[derive(Debug)]
pub struct SpillFrameParser {
    codec: SpillCodec,
    /// zstd: the length prefix of an incomplete `[len][frame]` record (the
    /// first `prefix_len` bytes).
    prefix: [u8; ZSTD_LEN_PREFIX],
    /// zstd: bytes of `prefix` held.
    prefix_len: usize,
    /// zstd: the incomplete frame body, allocated at its full length once the
    /// prefix is complete (it moves out as the emitted frame).
    body: Vec<u8>,
    /// zstd: the incomplete frame's body length (valid once `prefix_len` is 4).
    body_len: usize,
    /// BGZF: the framer (which keeps its own carry).
    bgzf: BgzfSliceFramer,
    /// BGZF: reused framer output.
    frames: Vec<SliceFrame>,
}

impl SpillFrameParser {
    /// A parser at the start of a spill body (after any file magic).
    #[must_use]
    pub fn new(codec: SpillCodec) -> Self {
        Self {
            codec,
            prefix: [0; ZSTD_LEN_PREFIX],
            prefix_len: 0,
            body: Vec::new(),
            body_len: 0,
            bgzf: BgzfSliceFramer::new(),
            frames: Vec::new(),
        }
    }

    /// Bytes of an incomplete frame held from earlier slices.
    #[must_use]
    pub fn carry_len(&self) -> usize {
        match self.codec {
            SpillCodec::Bgzf => self.bgzf.carry_len(),
            SpillCodec::Zstd => self.prefix_len + self.body.len(),
        }
    }

    /// Heap bytes the carried partial frame holds (its allocation, which is
    /// the frame's full size once that is known).
    #[must_use]
    pub fn carry_capacity(&self) -> usize {
        match self.codec {
            SpillCodec::Bgzf => self.bgzf.carry_capacity(),
            SpillCodec::Zstd => self.body.capacity(),
        }
    }

    /// Push every complete frame of `carry ++ slice` onto `out`; returns the
    /// number pushed.
    ///
    /// # Errors
    /// `InvalidData` for a malformed BGZF header or an oversized zstd length.
    pub fn push(&mut self, slice: &SliceLease, out: &mut Vec<RawFrame>) -> io::Result<usize> {
        match self.codec {
            SpillCodec::Bgzf => self.push_bgzf(slice, out),
            SpillCodec::Zstd => self.push_zstd(slice, out),
        }
    }

    /// End of the spill body.
    ///
    /// # Errors
    /// `UnexpectedEof` when a partial frame is still carried.
    pub fn finish(&mut self) -> io::Result<()> {
        match self.codec {
            SpillCodec::Bgzf => self.bgzf.finish(),
            SpillCodec::Zstd if self.carry_len() == 0 => Ok(()),
            SpillCodec::Zstd => Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                format!("truncated zstd spill frame: {} trailing bytes", self.carry_len()),
            )),
        }
    }

    fn push_bgzf(&mut self, slice: &SliceLease, out: &mut Vec<RawFrame>) -> io::Result<usize> {
        self.frames.clear();
        let n = self.bgzf.push(slice, &mut self.frames)?;
        out.extend(self.frames.drain(..).map(|f| match f {
            SliceFrame::Within(r) => RawFrame::borrowed(slice, r),
            SliceFrame::Carried(v) => RawFrame::Owned(v),
        }));
        Ok(n)
    }

    fn push_zstd(&mut self, slice: &SliceLease, out: &mut Vec<RawFrame>) -> io::Result<usize> {
        let before = out.len();
        let mut pos = 0usize;
        if self.prefix_len > 0 {
            if self.prefix_len < ZSTD_LEN_PREFIX {
                let take = (ZSTD_LEN_PREFIX - self.prefix_len).min(slice.len());
                self.prefix[self.prefix_len..self.prefix_len + take]
                    .copy_from_slice(&slice[..take]);
                self.prefix_len += take;
                pos = take;
                if self.prefix_len < ZSTD_LEN_PREFIX {
                    return Ok(0);
                }
                self.start_body(checked_zstd_len(&self.prefix)?);
            }
            let need = self.body_len - self.body.len();
            if slice.len() - pos < need {
                self.body.extend_from_slice(&slice[pos..]);
                return Ok(0);
            }
            self.body.extend_from_slice(&slice[pos..pos + need]);
            pos += need;
            self.prefix_len = 0;
            out.push(RawFrame::Owned(std::mem::take(&mut self.body)));
        }
        while slice.len() - pos >= ZSTD_LEN_PREFIX {
            let len = checked_zstd_len(&slice[pos..pos + ZSTD_LEN_PREFIX])?;
            let body = pos + ZSTD_LEN_PREFIX;
            if slice.len() - body < len {
                break;
            }
            out.push(RawFrame::borrowed(slice, body..body + len));
            pos = body + len;
        }
        let tail = &slice[pos..];
        if tail.len() < ZSTD_LEN_PREFIX {
            self.prefix[..tail.len()].copy_from_slice(tail);
            self.prefix_len = tail.len();
        } else {
            self.prefix.copy_from_slice(&tail[..ZSTD_LEN_PREFIX]);
            self.prefix_len = ZSTD_LEN_PREFIX;
            self.start_body(checked_zstd_len(&self.prefix)?);
            self.body.extend_from_slice(&tail[ZSTD_LEN_PREFIX..]);
        }
        Ok(out.len() - before)
    }

    /// Begin carrying a frame body of `len` bytes, allocated once at that
    /// size (it moves out as the emitted frame, so it cannot be reused).
    fn start_body(&mut self, len: usize) {
        self.body_len = len;
        self.body = Vec::with_capacity(len);
    }
}

/// A zstd frame's length from its `u32 LE` prefix, rejected over
/// [`MAX_ZSTD_FRAME_BYTES`] (the one copy of that check and its message,
/// shared with the sequential reader's `read_length_prefix`).
pub(crate) fn checked_zstd_len(prefix: &[u8]) -> io::Result<usize> {
    let len = u32::from_le_bytes(prefix.try_into().expect("a 4-byte prefix")) as usize;
    if len > MAX_ZSTD_FRAME_BYTES {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "zstd spill frame length {len} exceeds MAX_ZSTD_FRAME_BYTES ({MAX_ZSTD_FRAME_BYTES}): file likely corrupted",
            ),
        ));
    }
    Ok(len)
}

/// Output staging-buffer capacity for a single decompressed block/frame. Mirrors
/// the worker pool's `ZSTD_FRAME_DECOMP_CAP` / `BgzfDecompress` scratch sizing so
/// the mimalloc size-class reuse pattern matches.
const SCRATCH_CAP: usize = 256 * 1024;

/// Ceiling for sizing the scratch buffer up from [`SCRATCH_CAP`].
///
/// A single BAM record can legitimately exceed `SCRATCH_CAP` (long reads), so
/// the reader sizes up to fit rather than refusing — but not without limit: the
/// size comes from the frame header, which a corrupt frame can claim to be
/// anything, so an unbounded allocation on garbage would take the process out.
///
/// This is deliberately the *same* constant the writer enforces, imported rather
/// than redeclared: if the reader capped lower than the writer accepted, a block
/// in between would write successfully and be unreadable.
use crate::spill_block::MAX_SPILL_BLOCK_LEN as MAX_SCRATCH_CAP;

/// Per-worker codec-aware decompressor. Holds a libdeflate decompressor (BGZF)
/// and a zstd decompressor plus reusable scratch buffers; one instance per
/// `SortSpillDecompress` worker copy.
pub struct SpillBlockDecompressor {
    bgzf: BgzfDecompressor,
    zstd: ZstdDecompressor<'static>,
    /// Reused output buffer for the BGZF path (`mem::replace`d into the result).
    bgzf_scratch: Vec<u8>,
    /// Reused output buffer for the zstd path (decompressed-into, then copied).
    zstd_buf: Vec<u8>,
    /// Reused input buffer for the zstd path (one compressed frame at a time).
    zstd_frame: Vec<u8>,
}

impl SpillBlockDecompressor {
    /// Construct a fresh decompressor.
    ///
    /// # Panics
    ///
    /// Panics if the zstd decompressor cannot be created (allocation failure).
    #[must_use]
    pub fn new() -> Self {
        Self {
            bgzf: BgzfDecompressor::new(),
            zstd: ZstdDecompressor::new().expect("zstd decompressor init"),
            bgzf_scratch: Vec::with_capacity(SCRATCH_CAP),
            zstd_buf: Vec::new(),
            zstd_frame: Vec::new(),
        }
    }

    /// Read and decompress up to `max` blocks from `reader` using `codec`,
    /// returning the decompressed block bytes. A result with fewer than `max`
    /// entries (including an empty `Vec`) signals that the reader reached a
    /// clean EOF.
    ///
    /// The reader must be positioned at a block boundary — for zstd, at the
    /// start of a `[len][frame]` record (i.e. the `ZSP1` file magic already
    /// consumed); for BGZF, at the start of a block. `slots_for_chunk_files`
    /// detects the codec from the file magic and positions the reader
    /// accordingly when opening each slot.
    ///
    /// # Errors
    ///
    /// Propagates I/O errors, BGZF/zstd decompression failures, and truncation
    /// (a partial length prefix or frame body at EOF).
    pub fn read_blocks<R: Read + ?Sized>(
        &mut self,
        reader: &mut R,
        codec: SpillCodec,
        max: usize,
    ) -> io::Result<Vec<Vec<u8>>> {
        match codec {
            SpillCodec::Bgzf => self.read_bgzf_blocks(reader, max),
            SpillCodec::Zstd => self.read_zstd_frames(reader, max),
        }
    }

    fn read_bgzf_blocks<R: Read + ?Sized>(
        &mut self,
        reader: &mut R,
        max: usize,
    ) -> io::Result<Vec<Vec<u8>>> {
        let raw_blocks = fgumi_bgzf::reader::read_raw_blocks(reader, max)?;
        let mut out = Vec::with_capacity(raw_blocks.len());
        for raw in raw_blocks {
            fgumi_bgzf::reader::decompress_block_slice_into(
                &raw.data,
                &mut self.bgzf,
                &mut self.bgzf_scratch,
            )?;
            // mem::replace the filled scratch out and re-allocate a fresh one,
            // so the consumer owns the bytes and the next decompress reuses a
            // same-size-class allocation (matches `BgzfDecompress`).
            out.push(std::mem::replace(&mut self.bgzf_scratch, Vec::with_capacity(SCRATCH_CAP)));
        }
        Ok(out)
    }

    /// Decompress one zstd frame into the reusable scratch buffer, sizing it up
    /// first if the output would not fit, and return the decompressed length.
    ///
    /// Both read paths funnel through here so the sizing policy cannot drift
    /// between them — they previously each held their own copy of the
    /// resize-then-decompress sequence, and only one of the two would have been
    /// fixed by a change to either.
    ///
    /// The size comes from the frame header rather than from retrying on a
    /// too-small-destination error. `decompress_to_buffer` surfaces zstd's
    /// `get_error_name` text through `io::Error`, and that wording is not an
    /// API contract — matching on it means a zstd bump could silently turn
    /// "size up and succeed" into a hard failure on a spill file this very
    /// process wrote. Reading the header also allocates exactly once instead of
    /// doubling into place.
    fn decompress_zstd_frame(&mut self, frame: &[u8]) -> io::Result<usize> {
        // Both spill writers compress each block in a single
        // `Compressor::compress` call, which records the content size in the
        // frame header, so a missing size means the frame did not come from us.
        let content_size = zstd::zstd_safe::get_frame_content_size(frame)
            .map_err(|e| io::Error::other(format!("zstd spill frame header: {e}")))?
            .ok_or_else(|| {
                io::Error::other("zstd spill frame carries no content size (not written by fgumi)")
            })?;
        let needed = usize::try_from(content_size).map_err(|_| {
            io::Error::other(format!("zstd spill frame content size {content_size} exceeds usize"))
        })?;
        if needed > MAX_SCRATCH_CAP {
            return Err(io::Error::other(format!(
                "zstd spill frame decompresses to {needed} bytes, over the \
                 {MAX_SCRATCH_CAP}-byte scratch ceiling (or the frame is corrupt)"
            )));
        }
        // Reuse by CAPACITY, not by length. `decompress_to_buffer` sets the
        // vector's length to the decompressed byte count, so after any frame
        // `len` is the size of that frame's output, not the buffer size — a
        // `resize(want, 0)` guarded on `len` would therefore zero-fill the
        // difference on essentially every frame, immediately before zstd
        // overwrites it. `clear` + `reserve` keeps the allocation and zeroes
        // nothing; `reserve` is a no-op once capacity is high enough, so the
        // buffer still never shrinks.
        let want = needed.max(SCRATCH_CAP);
        self.zstd_buf.clear();
        self.zstd_buf.reserve(want);
        self.zstd
            .decompress_to_buffer(frame, &mut self.zstd_buf)
            .map_err(|e| io::Error::other(format!("zstd spill frame decompress: {e}")))
    }

    fn read_zstd_frames<R: Read + ?Sized>(
        &mut self,
        reader: &mut R,
        max: usize,
    ) -> io::Result<Vec<Vec<u8>>> {
        let mut out = Vec::with_capacity(max);
        for _ in 0..max {
            // `read_length_prefix` returns `Ok(None)` only at a clean frame
            // boundary EOF; a 1–3 byte partial prefix surfaces as an error.
            let Some(frame_len) = read_length_prefix(reader)? else {
                break;
            };
            // Read into spare capacity rather than `resize(frame_len, 0)` +
            // `read_exact`: the resize zero-fills the whole frame on every
            // iteration only to overwrite it immediately. `read_to_end` appends,
            // so the reuse is by capacity and nothing is zeroed. It also stops
            // short instead of erroring, so the short read is checked here to
            // keep `read_exact`'s `UnexpectedEof` semantics.
            self.zstd_frame.clear();
            self.zstd_frame.reserve(frame_len);
            Read::take(&mut *reader, frame_len as u64).read_to_end(&mut self.zstd_frame)?;
            if self.zstd_frame.len() != frame_len {
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    format!(
                        "truncated zstd spill frame: expected {frame_len} bytes, got {}",
                        self.zstd_frame.len(),
                    ),
                ));
            }
            // Move the frame out so `decompress_zstd_frame` can take `&mut self`,
            // then put it back to keep its allocation across frames. An error
            // here abandons the whole read, so not restoring it on that path
            // costs nothing.
            let frame = std::mem::take(&mut self.zstd_frame);
            let n = self.decompress_zstd_frame(&frame)?;
            self.zstd_frame = frame;
            out.push(self.zstd_buf[..n].to_vec());
        }
        Ok(out)
    }

    /// Read up to `max` *raw* (still-compressed) blocks from `reader` using
    /// `codec`, returning the compressed payloads **without** decompressing them
    /// (test oracle: the sequential reader the slice parser,
    /// [`SpillFrameParser`], is checked against).
    ///
    /// For BGZF each returned `Vec<u8>` is one complete raw block (header +
    /// compressed data + footer), exactly what [`Self::decompress_one`] expects;
    /// EOF-marker blocks are skipped. For zstd each is one raw frame body (the
    /// `[u32 LE len]` prefix is consumed here). A result shorter than `max`
    /// (including empty) signals a clean EOF, matching [`Self::read_blocks`].
    ///
    /// # Errors
    ///
    /// Propagates I/O errors and truncation (a partial BGZF block, or a partial
    /// zstd length prefix / frame body at EOF).
    #[cfg(any(test, feature = "test-utils"))]
    pub fn read_raw<R: Read + ?Sized>(
        &mut self,
        reader: &mut R,
        codec: SpillCodec,
        max: usize,
    ) -> io::Result<Vec<Vec<u8>>> {
        match codec {
            SpillCodec::Bgzf => {
                let raw_blocks = fgumi_bgzf::reader::read_raw_blocks(reader, max)?;
                Ok(raw_blocks.into_iter().map(|b| b.data).collect())
            }
            SpillCodec::Zstd => {
                let mut out = Vec::with_capacity(max);
                for _ in 0..max {
                    let Some(frame_len) = read_length_prefix(reader)? else {
                        break;
                    };
                    // Same zero-fill removal as `read_zstd_frames`: `vec![0u8; n]`
                    // zeroes the frame immediately before the read overwrites it.
                    // The *reuse* half of that fix does not transfer — `read_raw`
                    // hands each frame out as an owned `Vec` — but the zero-fill
                    // half does. `read_to_end` stops short instead of erroring, so
                    // the short read is checked to keep `read_exact`'s
                    // `UnexpectedEof` semantics.
                    let mut frame = Vec::with_capacity(frame_len);
                    Read::take(&mut *reader, frame_len as u64).read_to_end(&mut frame)?;
                    if frame.len() != frame_len {
                        return Err(io::Error::new(
                            io::ErrorKind::UnexpectedEof,
                            format!(
                                "truncated zstd spill frame: expected {frame_len} bytes, got {}",
                                frame.len(),
                            ),
                        ));
                    }
                    out.push(frame);
                }
                Ok(out)
            }
        }
    }

    /// Decompress a single raw block/frame (one parsed by [`SpillFrameParser`]).
    ///
    /// `raw` is one BGZF raw block (header + compressed data + footer) or one
    /// zstd frame body, per `codec`. Returns the decompressed bytes. Uses the
    /// worker's reusable scratch buffers, so this is cheap to call in a loop over
    /// a freshly-read batch.
    ///
    /// # Errors
    ///
    /// Propagates BGZF/zstd decompression failures (bad CRC, size mismatch, or a
    /// malformed frame).
    pub fn decompress_one(&mut self, codec: SpillCodec, raw: &[u8]) -> io::Result<Vec<u8>> {
        match codec {
            SpillCodec::Bgzf => {
                fgumi_bgzf::reader::decompress_block_slice_into(
                    raw,
                    &mut self.bgzf,
                    &mut self.bgzf_scratch,
                )?;
                Ok(std::mem::replace(&mut self.bgzf_scratch, Vec::with_capacity(SCRATCH_CAP)))
            }
            SpillCodec::Zstd => {
                let n = self.decompress_zstd_frame(raw)?;
                Ok(self.zstd_buf[..n].to_vec())
            }
        }
    }
}

impl Default for SpillBlockDecompressor {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pooled_chunk_writer::PooledChunkWriter;
    use crate::worker_pool::SortWorkerPool;
    use std::io::Cursor;
    use std::sync::Arc;

    /// A length prefix that promises more bytes than follow must surface
    /// `UnexpectedEof` from BOTH zstd read paths.
    ///
    /// `read_zstd_frames` and `read_raw` each read the frame body via
    /// `Read::take(..).read_to_end(..)`, which — unlike the `read_exact` it
    /// replaced — *stops short* instead of erroring. The explicit length check is
    /// what restores the loud failure, and it exists in two places, so this is a
    /// guard-set-parity case: drop the check in either path and a truncated spill
    /// gets silently decompressed as a short frame with every test still green.
    #[rstest::rstest]
    #[case::read_zstd_frames(false)]
    #[case::read_raw(true)]
    fn a_truncated_zstd_frame_is_unexpected_eof(#[case] via_read_raw: bool) {
        // `[u32 LE 64]` promising 64 body bytes, followed by only 8.
        let mut bytes = 64u32.to_le_bytes().to_vec();
        bytes.extend_from_slice(&[0xEE; 8]);

        let mut d = SpillBlockDecompressor::new();
        let mut cursor = Cursor::new(bytes);
        let err = if via_read_raw {
            d.read_raw(&mut cursor, SpillCodec::Zstd, 4).expect_err("truncated frame must error")
        } else {
            d.read_zstd_frames(&mut cursor, 4).expect_err("truncated frame must error")
        };
        assert_eq!(
            err.kind(),
            io::ErrorKind::UnexpectedEof,
            "a short frame must be UnexpectedEof, not a silent short read; got: {err}",
        );
        assert!(
            err.to_string().contains("truncated zstd spill frame"),
            "the error must name the truncation; got: {err}",
        );
    }

    /// A zstd frame with no recorded content size must be rejected, not
    /// decompressed into an unsized buffer.
    ///
    /// The reader sizes its scratch from the frame header, so a frame that omits
    /// the content size has nothing to size from. Both spill writers compress
    /// each block in one `bulk::compress` call, which always records it — so a
    /// size-less frame did not come from fgumi, and guessing at its size is the
    /// wrong response. `stream::encode_all` produces exactly that shape (verified
    /// against zstd 0.13: bulk records the size, the streaming APIs do not), so
    /// this is a real frame rather than a hand-forged header.
    #[test]
    fn a_zstd_frame_without_a_content_size_is_rejected() {
        let payload = vec![9u8; 5000];
        let sizeless = zstd::stream::encode_all(&payload[..], 1).expect("stream-compress");
        assert!(
            zstd::zstd_safe::get_frame_content_size(&sizeless).expect("readable header").is_none(),
            "precondition: the streaming API must omit the content size",
        );

        let mut d = SpillBlockDecompressor::new();
        let err = d
            .decompress_one(SpillCodec::Zstd, &sizeless)
            .expect_err("a size-less frame must be rejected");
        assert!(
            err.to_string().contains("carries no content size"),
            "the error must say the size is missing; got: {err}",
        );
    }

    /// A single zstd block larger than `SCRATCH_CAP` must survive the round
    /// trip. Before the reader sized its scratch to the frame, such a block
    /// wrote successfully and then failed on read — a spill file that could be
    /// produced but not consumed. Refusing the write was never the answer: a
    /// long-read BAM record alone can exceed 256 KiB.
    ///
    /// `compress_block` does now bound what it accepts, but at
    /// `MAX_SPILL_BLOCK_LEN` (64 MiB), which is the reader's ceiling too — the
    /// two agree by construction. This block is 768 KiB, above `SCRATCH_CAP`
    /// and far below that bound, so it exercises scratch growth specifically
    /// and not the writer's limit (which
    /// `compress_block_accepts_exactly_what_the_reader_can_decompress` covers).
    #[test]
    fn oversized_zstd_block_round_trips_through_the_reader() {
        let raw: Vec<u8> = (0..(SCRATCH_CAP * 3)).map(|i| u8::try_from(i % 251).unwrap()).collect();
        assert!(raw.len() > SCRATCH_CAP, "precondition: larger than the fixed scratch");

        let mut compressor =
            crate::spill_block::SpillBlockCompressor::new(SpillCodec::Zstd, 3).expect("compressor");
        let framed = compressor.compress_block(&raw).expect("compress accepts an oversized block");
        // `compress_block` prefixes the zstd frame with a u32 length; the
        // block-level entry point takes the bare frame.
        let frame = &framed[4..];

        let mut decompressor = SpillBlockDecompressor::new();
        let out = decompressor
            .decompress_one(SpillCodec::Zstd, frame)
            .expect("reader must grow its scratch rather than fail");
        assert_eq!(out, raw, "oversized block round-trips byte-for-byte");

        // The grown buffer is reused, so a normal-sized block still works after.
        let small = vec![7u8; 1024];
        let small_framed = compressor.compress_block(&small).expect("compress");
        let small_frame = &small_framed[4..];
        let small_out =
            decompressor.decompress_one(SpillCodec::Zstd, small_frame).expect("decompress");
        assert_eq!(small_out, small, "a later small block is unaffected by the growth");
    }

    /// A spill chunk written with codec X, read back through
    /// `SpillBlockDecompressor`, must yield the exact `[len][record]` byte
    /// stream that was written — for BOTH codecs. This is the round-trip the
    /// streaming merge depends on.
    fn roundtrip(codec: SpillCodec) {
        use crate::keys::{RawCoordinateKey, RawSortKey};
        // Build framed records whose total size (~440 KB) far exceeds the
        // writer's ~64 KB block cap, so the chunk is written as MANY blocks /
        // frames. That is what actually exercises the per-block decompress and
        // the downstream record reassembly across block boundaries — the whole
        // reason `read_blocks` returns per-block `Vec<u8>`s. Each record is
        // stamped with its index at both ends so a misaligned reassembly is
        // caught, not just a stream of zeros.
        let records: Vec<Vec<u8>> = (0usize..80)
            .map(|i| {
                let size = 3000 + (i % 13) * 500; // 3000..9000 bytes
                let mut r = vec![0u8; size];
                let stamp = u8::try_from(i % 251).expect("i % 251 fits u8");
                r[0] = stamp;
                let last = r.len() - 1;
                r[last] = stamp;
                r
            })
            .collect();

        let dir = tempfile::tempdir().expect("tempdir");
        let path = dir.path().join("chunk.spill");
        let pool = Arc::new(SortWorkerPool::new(1, 1, 6, codec));
        {
            let mut w = PooledChunkWriter::<RawCoordinateKey>::new(Arc::clone(&pool), &path, codec)
                .expect("writer");
            for rec in &records {
                let key = RawCoordinateKey::extract_from_record(rec);
                w.write_record(&key, rec).expect("write");
            }
            w.start_finish().expect("start_finish").wait().expect("finish");
        }

        // Re-open + position past the codec magic exactly like
        // slots_for_chunk_files does.
        let mut file = std::fs::File::open(&path).expect("open");
        let mut magic = [0u8; 4];
        let filled = crate::external::read_exact_or_eof(&mut file, &mut magic).expect("magic");
        let detected = if filled {
            SpillCodec::from_magic(&magic).unwrap_or(SpillCodec::Bgzf)
        } else {
            SpillCodec::Bgzf
        };
        assert_eq!(detected, codec, "codec must be detected from the file magic");
        if matches!(detected, SpillCodec::Bgzf) {
            use std::io::Seek;
            file.seek(std::io::SeekFrom::Start(0)).expect("seek");
        }

        let mut reader = std::io::BufReader::new(file);
        let mut dec = SpillBlockDecompressor::new();
        let mut all = Vec::new();
        loop {
            let blocks = dec.read_blocks(&mut reader, detected, 4).expect("read_blocks");
            let got = blocks.len();
            for b in blocks {
                all.extend_from_slice(&b);
            }
            if got < 4 {
                break;
            }
        }

        // The decompressed stream is `[u32 LE len][record]` per record, in order.
        let mut cursor = Cursor::new(&all);
        let mut read_back: Vec<Vec<u8>> = Vec::new();
        let mut len_buf = [0u8; 4];
        while Read::read(&mut cursor, &mut len_buf).map(|n| n == 4).unwrap_or(false) {
            let len = u32::from_le_bytes(len_buf) as usize;
            let mut rec = vec![0u8; len];
            cursor.read_exact(&mut rec).expect("record body");
            read_back.push(rec);
        }
        assert_eq!(read_back, records, "round-trip mismatch for {codec:?}");
    }

    #[test]
    fn roundtrip_bgzf() {
        roundtrip(SpillCodec::Bgzf);
    }

    #[test]
    fn roundtrip_zstd() {
        roundtrip(SpillCodec::Zstd);
    }

    /// The split `read_raw` + `decompress_one` path (used by the block-parallel
    /// `SortSpillDecompress`) must yield the exact same per-block bytes as the
    /// inline `read_blocks` path. We read the same chunk twice and compare.
    fn read_raw_matches_read_blocks(codec: SpillCodec) {
        use crate::keys::{RawCoordinateKey, RawSortKey};

        let records: Vec<Vec<u8>> = (0usize..80)
            .map(|i| {
                let size = 3000 + (i % 13) * 500;
                let mut r = vec![0u8; size];
                let stamp = u8::try_from(i % 251).expect("i % 251 fits u8");
                r[0] = stamp;
                let last = r.len() - 1;
                r[last] = stamp;
                r
            })
            .collect();

        let dir = tempfile::tempdir().expect("tempdir");
        let path = dir.path().join("chunk.spill");
        let pool = Arc::new(SortWorkerPool::new(1, 1, 6, codec));
        {
            let mut w = PooledChunkWriter::<RawCoordinateKey>::new(Arc::clone(&pool), &path, codec)
                .expect("writer");
            for rec in &records {
                let key = RawCoordinateKey::extract_from_record(rec);
                w.write_record(&key, rec).expect("write");
            }
            w.start_finish().expect("start_finish").wait().expect("finish");
        }

        // Helper: open the chunk positioned past any codec magic.
        let open = || {
            let mut file = std::fs::File::open(&path).expect("open");
            let mut magic = [0u8; 4];
            let filled = crate::external::read_exact_or_eof(&mut file, &mut magic).expect("magic");
            let detected = if filled {
                SpillCodec::from_magic(&magic).unwrap_or(SpillCodec::Bgzf)
            } else {
                SpillCodec::Bgzf
            };
            if matches!(detected, SpillCodec::Bgzf) {
                use std::io::Seek;
                file.seek(std::io::SeekFrom::Start(0)).expect("seek");
            }
            (std::io::BufReader::new(file), detected)
        };

        // Path A: inline read_blocks.
        let (mut reader_a, detected) = open();
        let mut dec_a = SpillBlockDecompressor::new();
        let mut blocks_a: Vec<Vec<u8>> = Vec::new();
        loop {
            let blocks = dec_a.read_blocks(&mut reader_a, detected, 4).expect("read_blocks");
            let got = blocks.len();
            blocks_a.extend(blocks);
            if got < 4 {
                break;
            }
        }

        // Path B: read_raw + decompress_one.
        let (mut reader_b, _) = open();
        let mut dec_b = SpillBlockDecompressor::new();
        let mut blocks_b: Vec<Vec<u8>> = Vec::new();
        loop {
            let raws = dec_b.read_raw(&mut reader_b, detected, 4).expect("read_raw");
            let got = raws.len();
            for raw in &raws {
                blocks_b.push(dec_b.decompress_one(detected, raw).expect("decompress_one"));
            }
            if got < 4 {
                break;
            }
        }

        assert_eq!(
            blocks_a, blocks_b,
            "split read_raw path must match inline read_blocks ({codec:?})"
        );
    }

    #[test]
    fn read_raw_matches_read_blocks_bgzf() {
        read_raw_matches_read_blocks(SpillCodec::Bgzf);
    }

    #[test]
    fn read_raw_matches_read_blocks_zstd() {
        read_raw_matches_read_blocks(SpillCodec::Zstd);
    }

    // ---- SpillFrameParser ----

    /// A spill file of `n` blocks of pseudo-random (poorly compressible)
    /// records-sized payloads in `codec`: magic, blocks, trailer.
    fn spill_stream(codec: SpillCodec, n: usize) -> Vec<u8> {
        let mut c = crate::spill_block::SpillBlockCompressor::new(codec, 1).unwrap();
        let mut out = crate::spill_block::spill_magic(codec).to_vec();
        let mut state = 0x2545_f491_4f6c_dd1du64 ^ n as u64;
        for i in 0..n {
            let len = 1_000 + (i * 7_919) % 90_000;
            let raw: Vec<u8> = (0..len)
                .map(|_| {
                    state ^= state << 13;
                    state ^= state >> 7;
                    state ^= state << 17;
                    state.to_le_bytes()[2] & 0x3f
                })
                .collect();
            out.extend_from_slice(&c.compress_block(&raw).unwrap());
        }
        out.extend_from_slice(crate::spill_block::spill_trailer(codec));
        out
    }

    fn body_start(codec: SpillCodec) -> usize {
        crate::spill_block::spill_magic(codec).len()
    }

    proptest::proptest! {
        /// Both codecs: frames from random slice cuts (including
        /// tiny ones that end inside a length prefix or a header) equal
        /// `read_raw` over the whole stream; whole-in-slice frames are borrowed
        /// within bounds; the carry never exceeds one frame; every slice returns
        /// to its pool.
        #[test]
        fn parser_matches_read_raw_over_random_cuts(
            zstd in proptest::bool::ANY,
            n in 1usize..40,
            cuts in proptest::collection::vec(1usize..300_000, 1..50),
            tiny in proptest::bool::ANY,
        ) {
            use fgumi_bam_io::pread::SliceBufferPool;
            let codec = if zstd { SpillCodec::Zstd } else { SpillCodec::Bgzf };
            let stream = spill_stream(codec, n);
            // `read_raw` pre-allocates `max`, so bound it by the stream length.
            let oracle = SpillBlockDecompressor::new()
                .read_raw(&mut &stream[body_start(codec)..], codec, stream.len())
                .unwrap();
            let pool = SliceBufferPool::new(64);
            let mut p = SpillFrameParser::new(codec);
            let mut got = Vec::new();
            let mut pos = body_start(codec);
            let mut i = 0;
            while pos < stream.len() {
                let cut = if tiny { 1 + cuts[i % cuts.len()] % 23 } else { cuts[i % cuts.len()] };
                let len = cut.min(stream.len() - pos);
                i += 1;
                let lease = pool.lease(stream[pos..pos + len].to_vec());
                let mut out = Vec::new();
                p.push(&lease, &mut out).unwrap();
                // Only the frame completing the carry may be owned: every frame
                // wholly inside the slice is borrowed.
                let owned = out.iter().filter(|f| matches!(f, RawFrame::Owned(_))).count();
                proptest::prop_assert!(owned <= 1, "{owned} owned frames from one slice");
                for f in out {
                    if let RawFrame::Borrowed { range, .. } = &f {
                        proptest::prop_assert!((range.end as usize) <= lease.len());
                    }
                    got.push(f.bytes().to_vec());
                }
                proptest::prop_assert!(p.carry_len() <= MAX_ZSTD_FRAME_BYTES + 18);
                pos += len;
            }
            p.finish().unwrap();
            proptest::prop_assert_eq!(got, oracle);
            drop(p);
            proptest::prop_assert_eq!(pool.resident_bytes(), 0);
        }
    }

    /// A `[u32 len]` over the cap is `InvalidData` with the existing message.
    #[test]
    fn zstd_frame_length_over_cap_is_invalid_data() {
        let pool = fgumi_bam_io::pread::SliceBufferPool::new(1);
        let mut bytes = u32::try_from(MAX_ZSTD_FRAME_BYTES + 1).unwrap().to_le_bytes().to_vec();
        bytes.extend_from_slice(&[0u8; 16]);
        let mut p = SpillFrameParser::new(SpillCodec::Zstd);
        let err = p.push(&pool.lease(bytes), &mut Vec::new()).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert!(err.to_string().contains("exceeds MAX_ZSTD_FRAME_BYTES"), "{err}");
    }

    /// The same rejection when the length prefix is split across slices.
    #[test]
    fn zstd_frame_length_over_cap_split_across_slices_is_invalid_data() {
        let pool = fgumi_bam_io::pread::SliceBufferPool::new(2);
        let len = u32::try_from(MAX_ZSTD_FRAME_BYTES + 1).unwrap().to_le_bytes();
        let mut p = SpillFrameParser::new(SpillCodec::Zstd);
        assert_eq!(p.push(&pool.lease(len[..2].to_vec()), &mut Vec::new()).unwrap(), 0);
        let err = p.push(&pool.lease(len[2..].to_vec()), &mut Vec::new()).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
    }

    /// A slice holding only the BGZF EOF marker yields no frame and leaves no
    /// carry.
    #[test]
    fn eof_marker_only_slice_yields_no_frame() {
        let pool = fgumi_bam_io::pread::SliceBufferPool::new(1);
        let mut p = SpillFrameParser::new(SpillCodec::Bgzf);
        let mut out = Vec::new();
        assert_eq!(p.push(&pool.lease(fgumi_bgzf::BGZF_EOF.to_vec()), &mut out).unwrap(), 0);
        assert_eq!(p.carry_len(), 0);
        p.finish().unwrap();
    }

    /// A frame carried across slices is allocated once at its full size —
    /// whether the cut falls inside the zstd length prefix / BGZF header or
    /// after it — and the carry reports that allocation while partial, so a
    /// stash charging the carry charges exactly what the emitted frame holds.
    #[rstest::rstest]
    #[case::bgzf_after_header(SpillCodec::Bgzf, &[0.6])]
    #[case::bgzf_inside_header(SpillCodec::Bgzf, &[0.001, 0.6])]
    #[case::zstd_after_prefix(SpillCodec::Zstd, &[0.6])]
    #[case::zstd_inside_prefix(SpillCodec::Zstd, &[0.0, 0.6])]
    fn a_carried_frame_is_allocated_once_at_its_size(
        #[case] codec: SpillCodec,
        #[case] cuts: &[f64],
    ) {
        let stream = spill_stream(codec, 1);
        let body = &stream[body_start(codec)..];
        let frame = &SpillBlockDecompressor::new().read_raw(&mut &body[..], codec, 4).unwrap()[0];
        let record = if codec == SpillCodec::Zstd { frame.len() + 4 } else { frame.len() };
        let pool = fgumi_bam_io::pread::SliceBufferPool::new(4);
        let mut p = SpillFrameParser::new(codec);
        let mut out = Vec::new();
        let mut from = 0;
        #[allow(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            clippy::cast_precision_loss,
            reason = "test cut points"
        )]
        for to in cuts.iter().map(|c| ((record as f64 * c) as usize).max(2)) {
            p.push(&pool.lease(body[from..to].to_vec()), &mut out).unwrap();
            assert!(out.is_empty());
            if to >= 18 {
                assert_eq!(p.carry_capacity(), frame.len(), "the carry holds the frame's size");
            }
            from = to;
        }
        p.push(&pool.lease(body[from..record].to_vec()), &mut out).unwrap();
        let [RawFrame::Owned(v)] = &out[..] else { panic!("one carried frame: {out:?}") };
        assert_eq!((v.len(), v.capacity()), (frame.len(), frame.len()));
        assert_eq!(p.carry_capacity(), 0);
    }

    /// A stream that ends mid-frame fails `finish` with `UnexpectedEof`.
    #[rstest::rstest]
    #[case::bgzf(SpillCodec::Bgzf)]
    #[case::zstd(SpillCodec::Zstd)]
    fn truncated_stream_fails_finish(#[case] codec: SpillCodec) {
        let stream = spill_stream(codec, 3);
        let pool = fgumi_bam_io::pread::SliceBufferPool::new(1);
        let mut p = SpillFrameParser::new(codec);
        p.push(&pool.lease(stream[body_start(codec)..stream.len() - 40].to_vec()), &mut Vec::new())
            .unwrap();
        assert_eq!(p.finish().unwrap_err().kind(), io::ErrorKind::UnexpectedEof);
    }
}

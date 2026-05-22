//! Streaming reader for "ZSP1"-format spill files.
//!
//! On-disk format:
//! ```text
//! [4-byte ASCII "ZSP1"]
//! [u32 LE compressed-len][zstd frame] x N
//! ```
//!
//! This module presents a [`Read`]-compatible view over the decompressed
//! content of such a file, so the legacy
//! [`GenericKeyedChunkReader`](crate::external::GenericKeyedChunkReader) record
//! parser can consume zstd-compressed spills without further changes. The
//! worker-pool reader takes the faster, per-frame parallel-decompress path;
//! this stream reader is the fallback used by the consolidation merge and the
//! indexing path.

use crate::codec::ZSPILL_MAGIC;
use std::io::{self, Read};

/// Cap on the uncompressed size of a single zstd spill frame. Production
/// frames are bounded by the producer's staging buffer
/// (`BGZF_MAX_BLOCK_SIZE` + padding ~= 68 KB), so 256 KiB leaves comfortable
/// slack. If a frame ever decompresses to more bytes than this,
/// `zstd::bulk::Decompressor::decompress_to_buffer` surfaces a clear error
/// rather than silently truncating. Kept in sync with `ZSTD_FRAME_DECOMP_CAP`
/// in `worker_pool.rs`.
const FRAME_DECOMP_CAP: usize = 256 * 1024;

/// Streaming decompressor for "ZSP1" spill files.
pub struct ZspillStreamReader<R: Read> {
    inner: R,
    decompressor: zstd::bulk::Decompressor<'static>,
    buf: Vec<u8>,
    pos: usize,
}

impl<R: Read> ZspillStreamReader<R> {
    /// Open a streaming zstd reader, consuming the four-byte file magic.
    pub fn new(mut inner: R) -> io::Result<Self> {
        let mut magic = [0u8; 4];
        inner.read_exact(&mut magic)?;
        if magic != ZSPILL_MAGIC {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("missing ZSPILL_MAGIC: got {magic:?}, expected {ZSPILL_MAGIC:?}"),
            ));
        }
        let decompressor = zstd::bulk::Decompressor::new()
            .map_err(|e| io::Error::other(format!("zstd decompressor init: {e}")))?;
        Ok(Self { inner, decompressor, buf: Vec::with_capacity(FRAME_DECOMP_CAP), pos: 0 })
    }

    /// Fill the internal buffer with the next decompressed frame. Returns
    /// `false` on a clean EOF at a frame boundary; propagates an error if the
    /// length prefix is truncated, exceeds `MAX_ZSTD_FRAME_BYTES`, or the
    /// frame body is short.
    fn fill_next_frame(&mut self) -> io::Result<bool> {
        let Some(frame_len) = crate::worker_pool::read_length_prefix(&mut self.inner)? else {
            return Ok(false);
        };
        let mut frame = vec![0u8; frame_len];
        self.inner.read_exact(&mut frame)?;

        // Resize to the cap so we have writable room; then truncate to the
        // actual decompressed length. Capacity stays at FRAME_DECOMP_CAP for
        // the next call.
        self.buf.resize(FRAME_DECOMP_CAP, 0);
        let n = self
            .decompressor
            .decompress_to_buffer(&frame, &mut self.buf[..])
            .map_err(|e| io::Error::other(format!("zstd decompress: {e}")))?;
        self.buf.truncate(n);
        self.pos = 0;
        Ok(true)
    }
}

impl<R: Read> Read for ZspillStreamReader<R> {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        if out.is_empty() {
            return Ok(0);
        }
        if self.pos >= self.buf.len() && !self.fill_next_frame()? {
            return Ok(0);
        }
        let avail = self.buf.len() - self.pos;
        let n = avail.min(out.len());
        out[..n].copy_from_slice(&self.buf[self.pos..self.pos + n]);
        self.pos += n;
        Ok(n)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
    use std::io::Cursor;

    /// Round-trip a sequence of independent frames through the stream reader.
    #[test]
    fn roundtrip_three_frames() {
        let frames: Vec<Vec<u8>> = vec![
            b"hello, world. this is frame zero.".to_vec(),
            b"another payload for frame one, slightly longer than the first.".to_vec(),
            b"and a final frame.".to_vec(),
        ];

        let mut compressor = zstd::bulk::Compressor::new(1).expect("compressor");
        let mut buf: Vec<u8> = ZSPILL_MAGIC.to_vec();
        for f in &frames {
            let frame = compressor.compress(f).expect("compress");
            let frame_len = u32::try_from(frame.len()).expect("test frame fits in u32");
            buf.extend_from_slice(&frame_len.to_le_bytes());
            buf.extend_from_slice(&frame);
        }

        let mut reader = ZspillStreamReader::new(Cursor::new(buf)).expect("open");
        let mut out = Vec::new();
        reader.read_to_end(&mut out).expect("read_to_end");
        let expected: Vec<u8> = frames.into_iter().flatten().collect();
        assert_eq!(out, expected);
    }

    /// Empty body (magic only) reads as EOF.
    #[test]
    fn roundtrip_empty_body() {
        let buf: Vec<u8> = ZSPILL_MAGIC.to_vec();
        let mut reader = ZspillStreamReader::new(Cursor::new(buf)).expect("open");
        let mut out = Vec::new();
        reader.read_to_end(&mut out).expect("read_to_end");
        assert!(out.is_empty());
    }

    /// Missing magic surfaces as `InvalidData`.
    #[test]
    fn missing_magic_errors() {
        let Err(err) = ZspillStreamReader::new(Cursor::new(b"NOTZ".to_vec())) else {
            panic!("expected error");
        };
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
    }

    /// A truncated length prefix (1–3 bytes after magic) must surface as
    /// `UnexpectedEof`, not be silently treated as a clean stream end.
    #[rstest]
    #[case(1)]
    #[case(2)]
    #[case(3)]
    fn truncated_length_prefix_errors(#[case] partial: usize) {
        let mut buf: Vec<u8> = ZSPILL_MAGIC.to_vec();
        buf.extend_from_slice(&[0xAB; 4][..partial]);
        let mut reader = ZspillStreamReader::new(Cursor::new(buf)).expect("open");
        let mut out = Vec::new();
        let err = reader.read_to_end(&mut out).expect_err("expected truncation error");
        assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof, "partial={partial} got {err:?}",);
    }

    /// A truncated frame body (length prefix says N, body is < N) must surface
    /// as `UnexpectedEof`.
    #[test]
    fn truncated_frame_body_errors() {
        let mut compressor = zstd::bulk::Compressor::new(1).expect("compressor");
        let frame = compressor.compress(b"hello").expect("compress");
        let frame_len = u32::try_from(frame.len()).expect("fits");

        let mut buf: Vec<u8> = ZSPILL_MAGIC.to_vec();
        buf.extend_from_slice(&frame_len.to_le_bytes());
        // Push only the first half of the frame.
        buf.extend_from_slice(&frame[..frame.len() / 2]);

        let mut reader = ZspillStreamReader::new(Cursor::new(buf)).expect("open");
        let mut out = Vec::new();
        let err = reader.read_to_end(&mut out).expect_err("expected truncation error");
        assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
    }

    /// A length prefix larger than `MAX_ZSTD_FRAME_BYTES` must be rejected as
    /// corruption rather than triggering an unbounded allocation.
    #[test]
    fn oversized_length_prefix_errors() {
        let bogus_len: u32 =
            u32::try_from(crate::worker_pool::MAX_ZSTD_FRAME_BYTES).expect("fits") + 1;
        let mut buf: Vec<u8> = ZSPILL_MAGIC.to_vec();
        buf.extend_from_slice(&bogus_len.to_le_bytes());

        let mut reader = ZspillStreamReader::new(Cursor::new(buf)).expect("open");
        let mut out = Vec::new();
        let err = reader.read_to_end(&mut out).expect_err("expected length-cap error");
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
    }
}

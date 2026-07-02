//! Codec-aware streaming block decompressor for the typed-step
//! `SortSpillDecompress` pipeline step.
//!
//! A spill chunk is either a BGZF block-stream (self-framed `1f 8b` blocks) or
//! a zstd "ZSP1" stream (`[u32 LE frame-len][zstd frame]` records after the
//! 4-byte file magic). Both standalone `fgumi sort` (via the worker pool) and
//! the fused streaming sort write spill chunks with whichever
//! [`SpillCodec`](crate::codec::SpillCodec) the sorter is configured for; this
//! decompressor lets the streaming `SortSpillDecompress` step read either.
//!
//! Each decompressed BGZF block (or zstd frame) is returned as one `Vec<u8>`.
//! The block/frame boundaries are immaterial: the downstream `MergeDriver`
//! parses `[len][record]` records out of a slot's decompressed-block queue and
//! reassembles records that span block boundaries, so any chunking of the
//! decompressed byte stream is correct.

use std::io::{self, Read};

use libdeflater::Decompressor as BgzfDecompressor;
use zstd::bulk::Decompressor as ZstdDecompressor;

use crate::codec::SpillCodec;
use crate::worker_pool::read_length_prefix;

/// Output staging-buffer capacity for a single decompressed block/frame. Mirrors
/// the worker pool's `ZSTD_FRAME_DECOMP_CAP` / `BgzfDecompress` scratch sizing so
/// the mimalloc size-class reuse pattern matches.
const SCRATCH_CAP: usize = 256 * 1024;

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

    fn read_zstd_frames<R: Read + ?Sized>(
        &mut self,
        reader: &mut R,
        max: usize,
    ) -> io::Result<Vec<Vec<u8>>> {
        if self.zstd_buf.len() < SCRATCH_CAP {
            self.zstd_buf.resize(SCRATCH_CAP, 0);
        }
        let mut out = Vec::with_capacity(max);
        for _ in 0..max {
            // `read_length_prefix` returns `Ok(None)` only at a clean frame
            // boundary EOF; a 1–3 byte partial prefix surfaces as an error.
            let Some(frame_len) = read_length_prefix(reader)? else {
                break;
            };
            self.zstd_frame.clear();
            self.zstd_frame.resize(frame_len, 0);
            reader.read_exact(&mut self.zstd_frame)?;
            let n = self
                .zstd
                .decompress_to_buffer(&self.zstd_frame, &mut self.zstd_buf)
                .map_err(|e| io::Error::other(format!("zstd spill frame decompress: {e}")))?;
            out.push(self.zstd_buf[..n].to_vec());
        }
        Ok(out)
    }

    /// Read up to `max` *raw* (still-compressed) blocks from `reader` using
    /// `codec`, returning the compressed payloads **without** decompressing them.
    ///
    /// This is the read half of the block-parallel `SortSpillDecompress` path:
    /// the caller acquires the per-slot reader lock, calls `read_raw` to pull a
    /// batch of compressed blocks (which it sequence-tags), releases the lock,
    /// and then decompresses each block via [`Self::decompress_one`] *outside*
    /// the lock so multiple workers decompress one file's blocks concurrently.
    /// The read and the decompression of a given block still happen within a
    /// single `try_run` of a single worker — only the lock is released between
    /// them — which preserves the read-and-decompress-together invariant that
    /// the FIFO inline path also upholds.
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
                    let mut frame = vec![0u8; frame_len];
                    reader.read_exact(&mut frame)?;
                    out.push(frame);
                }
                Ok(out)
            }
        }
    }

    /// Decompress a single raw block/frame previously read by [`Self::read_raw`].
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
                if self.zstd_buf.len() < SCRATCH_CAP {
                    self.zstd_buf.resize(SCRATCH_CAP, 0);
                }
                let n = self
                    .zstd
                    .decompress_to_buffer(raw, &mut self.zstd_buf)
                    .map_err(|e| io::Error::other(format!("zstd spill frame decompress: {e}")))?;
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
}

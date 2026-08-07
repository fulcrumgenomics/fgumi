//! Block-granular spill compression kernel shared by the Phase-1 spill steps
//! (`SpillGather` → `SpillCompress` → `SpillWrite`).
//!
//! The legacy [`SyncSpillWriter`](crate::sync_spill_writer) streams a whole
//! sorted chunk through one stateful compressor on a single worker. To compress
//! a spill across the framework pool instead, the chunk is cut into raw blocks by
//! the (Serial) serialize step, each block is compressed **independently** by a
//! (Parallel) compress worker, and the writer concatenates the results per file.
//! This module holds the three pieces that keep that byte-compatible with the
//! existing spill readers:
//!
//! - [`frame_keyed_record_into`] — append one `[key?][u32 LE len][record]` record
//!   to a raw block buffer (the identical on-disk record layout
//!   `SyncSpillWriter::write_record` produces).
//! - [`SpillBlockCompressor`] — a per-worker, codec-fixed compressor whose
//!   [`compress_block`](SpillBlockCompressor::compress_block) turns one raw block
//!   into a self-contained, independently-decodable unit: framed BGZF block(s)
//!   for bgzf, or a `[u32 LE len][zstd frame]` for zstd.
//! - [`spill_magic`] / [`spill_trailer`] — the per-codec file prologue/epilogue
//!   the writer brackets the compressed blocks with.
//!
//! Because a BGZF/zstd spill stream is just independent blocks concatenated and
//! the readers stream across boundaries, **any** blocking of the same raw
//! `[key?][len][record]…` byte stream reads back to the identical records (see
//! the round-trip + re-blocking-independence tests below and the format notes in
//! [`sync_spill_writer`](crate::sync_spill_writer)). Block boundaries are
//! therefore the serialize step's free choice, and compression is a pure function
//! of each block — so the result is byte-for-byte independent of how many workers
//! compress concurrently. The one bound on that freedom is
//! [`MAX_SPILL_BLOCK_LEN`].

use std::io;

use fgumi_bgzf::BGZF_EOF;
use fgumi_bgzf::writer::InlineBgzfCompressor;
use zstd::bulk::Compressor as ZstdCompressor;

/// Largest raw block [`SpillBlockCompressor::compress_block`] will accept, and
/// the ceiling the reader sizes its decompression scratch to.
///
/// This is one contract with two enforcement points, which is why it lives here
/// (the writer side) and is imported by `spill_block_reader`, rather than being
/// declared next to the buffer it bounds. The reader has to cap *something*: the
/// decompressed size it reads out of a zstd frame header is attacker- (or
/// corruption-) controlled, so an uncapped allocation would take the process
/// out. Once the reader caps, the writer must reject the same size, or a block
/// above the cap writes successfully and is unreadable — a spill file this
/// process produced and cannot consume.
///
/// 64 MiB is far above any single BAM record in practice (long reads are ~1 MiB
/// at the extreme), so it bounds corruption without bounding real data.
pub(crate) const MAX_SPILL_BLOCK_LEN: usize = 64 * 1024 * 1024;

use crate::codec::{SpillCodec, ZSPILL_MAGIC};
use crate::keys::RawSortKey;

/// Append one keyed record to `block` in the spill record framing:
/// `[key bytes if !K::EMBEDDED_IN_RECORD][u32 LE record-len][record body]`.
///
/// This is the exact layout `SyncSpillWriter::write_record` writes; the
/// serialize step accumulates records into a raw block with this helper, and the
/// merge / Phase-2 readers parse it back unchanged.
///
/// # Errors
///
/// Returns an error if the record is longer than `u32::MAX` (cannot fit the
/// length prefix) or key serialization fails.
pub fn frame_keyed_record_into<K: RawSortKey>(
    block: &mut Vec<u8>,
    key: &K,
    record: &[u8],
) -> io::Result<()> {
    // Validate the length BEFORE writing the key, so an oversized record fails
    // loud without leaving partial (key) bytes in the block.
    let record_len = u32::try_from(record.len()).map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("BAM record too large ({} bytes) for a u32 length prefix", record.len()),
        )
    })?;
    if !K::EMBEDDED_IN_RECORD {
        // `Vec<u8>: Write`, so the key serializes straight into the block buffer.
        key.write_to(block)?;
    }
    block.extend_from_slice(&record_len.to_le_bytes());
    block.extend_from_slice(record);
    Ok(())
}

/// Per-worker, codec-fixed block compressor. Constructed once per compress
/// worker (like `BgzfCompress`'s `InlineBgzfCompressor`) and reused across the
/// blocks that worker handles.
pub enum SpillBlockCompressor {
    /// bgzf: each block becomes one or more framed BGZF blocks (header + deflate
    /// + footer). No EOF marker — that is the file trailer (see [`spill_trailer`]).
    Bgzf(InlineBgzfCompressor),
    /// zstd: each block becomes one `[u32 LE frame-len][zstd frame]` unit.
    Zstd(ZstdCompressor<'static>),
}

impl SpillBlockCompressor {
    /// Build a compressor for `codec` at `compression` (bgzf level — `0` writes
    /// framed *stored* blocks; or zstd level, which must be ≥ 1).
    ///
    /// # Errors
    ///
    /// Returns an error if the zstd compressor cannot be initialized.
    pub fn new(codec: SpillCodec, compression: u32) -> io::Result<Self> {
        match codec {
            SpillCodec::Bgzf => Ok(Self::Bgzf(InlineBgzfCompressor::new(compression))),
            SpillCodec::Zstd => {
                #[allow(clippy::cast_possible_wrap)]
                let compressor = ZstdCompressor::new(compression as i32).map_err(|e| {
                    io::Error::other(format!("zstd compressor init (level {compression}): {e}"))
                })?;
                Ok(Self::Zstd(compressor))
            }
        }
    }

    /// Compress one raw block into a self-contained, independently-decodable
    /// unit. Concatenating these units (bracketed by [`spill_magic`] /
    /// [`spill_trailer`]) reproduces the same decompressed stream as the
    /// streaming `SyncSpillWriter`, regardless of where the block boundaries fall.
    ///
    /// An empty block yields an empty `Vec` (no empty BGZF block / zstd frame is
    /// emitted), matching the streaming writer's "no empty frames" invariant.
    ///
    /// # Errors
    ///
    /// Returns an error if compression fails, if a zstd frame exceeds `u32::MAX`,
    /// or if `raw` exceeds `MAX_SPILL_BLOCK_LEN` — see that constant for why
    /// the write side is where that has to fail.
    pub fn compress_block(&mut self, raw: &[u8]) -> io::Result<Vec<u8>> {
        if raw.is_empty() {
            return Ok(Vec::new());
        }
        // Fail here rather than let the reader fail later. The block boundary is
        // the serialize step's free choice (see the module doc), so nothing in
        // the type system stops it picking one the reader cannot decompress —
        // and a spill file that writes successfully and cannot be read back is
        // the worst shape this failure can take. Same reason the length check in
        // `frame_keyed_record_into` runs before any bytes are written.
        if raw.len() > MAX_SPILL_BLOCK_LEN {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "spill block of {} bytes exceeds the {MAX_SPILL_BLOCK_LEN}-byte limit the \
                     reader can decompress",
                    raw.len(),
                ),
            ));
        }
        match self {
            Self::Bgzf(compressor) => {
                compressor.write_all(raw)?;
                // Flush forces a block boundary at the end of this raw block, so
                // the unit is independently decodable.
                compressor.flush()?;
                let mut out = Vec::new();
                for block in compressor.take_blocks() {
                    out.extend_from_slice(&block.data);
                }
                Ok(out)
            }
            Self::Zstd(compressor) => {
                let frame = compressor
                    .compress(raw)
                    .map_err(|e| io::Error::other(format!("zstd compress: {e}")))?;
                let frame_len = u32::try_from(frame.len()).map_err(|_| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        "zstd frame larger than 4 GiB cannot fit a u32 length prefix",
                    )
                })?;
                let mut out = Vec::with_capacity(4 + frame.len());
                out.extend_from_slice(&frame_len.to_le_bytes());
                out.extend_from_slice(&frame);
                Ok(out)
            }
        }
    }
}

/// File prologue for `codec`: `ZSPILL_MAGIC` for zstd, empty for bgzf (a BGZF
/// stream needs no prologue — each block self-identifies via the gzip magic).
#[must_use]
pub fn spill_magic(codec: SpillCodec) -> &'static [u8] {
    match codec {
        SpillCodec::Bgzf => &[],
        SpillCodec::Zstd => &ZSPILL_MAGIC,
    }
}

/// File epilogue for `codec`: the `BGZF_EOF` empty-block terminator for bgzf,
/// empty for zstd (the zstd spill format has no trailing marker).
#[must_use]
pub fn spill_trailer(codec: SpillCodec) -> &'static [u8] {
    match codec {
        SpillCodec::Bgzf => &BGZF_EOF,
        SpillCodec::Zstd => &[],
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::external::GenericKeyedChunkReader;
    use crate::inline::TemplateKey;
    use std::fs::File;
    use std::io::{BufWriter, Write as _};
    use std::path::Path;
    use tempfile::TempDir;

    #[allow(clippy::cast_possible_truncation)]
    fn make_key(i: u64) -> TemplateKey {
        TemplateKey::new(
            i as i32,
            i as i32,
            false,
            i32::MAX,
            i32::MAX,
            false,
            0,
            0,
            (0, false),
            i,
            false,
        )
    }

    #[allow(clippy::cast_possible_truncation)]
    fn sample_records(n: u64) -> Vec<(TemplateKey, Vec<u8>)> {
        (0..n).map(|i| (make_key(i), vec![(i % 256) as u8; 200 + (i as usize % 50)])).collect()
    }

    /// The writer must reject exactly what the reader cannot decompress.
    ///
    /// These two limits are one contract with two enforcement points, and the
    /// failure mode when they disagree is asymmetric: a writer that accepts more
    /// than the reader accepts produces a spill file this process wrote and
    /// cannot read back. Rejecting at write time turns that into a loud error
    /// on the block that caused it.
    ///
    /// Both cases allocate a buffer of the size under test, so the at-limit case
    /// really does build and compress 64 MiB. That is deliberate — it is the
    /// only way to show the limit is inclusive rather than off by one — and it
    /// is affordable because the buffer is zeroed, which zstd collapses almost
    /// instantly. The over-limit case is the cheap one: `compress_block` checks
    /// the length before touching the compressor, so it short-circuits.
    #[rstest::rstest]
    #[case::at_the_limit(MAX_SPILL_BLOCK_LEN, true)]
    #[case::one_over_the_limit(MAX_SPILL_BLOCK_LEN + 1, false)]
    fn compress_block_accepts_exactly_what_the_reader_can_decompress(
        #[case] len: usize,
        #[case] accepted: bool,
    ) {
        // `compress_block` checks the length before touching the compressor, so
        // a zeroed buffer of the boundary size exercises the branch without
        // needing meaningful content.
        let raw = vec![0u8; len];
        let mut compressor = SpillBlockCompressor::new(SpillCodec::Zstd, 1).expect("compressor");
        let result = compressor.compress_block(&raw);
        assert_eq!(
            result.is_ok(),
            accepted,
            "a {len}-byte block should {} at the {MAX_SPILL_BLOCK_LEN}-byte limit",
            if accepted { "be accepted" } else { "be rejected" },
        );
        if !accepted {
            let err = result.unwrap_err();
            assert_eq!(err.kind(), io::ErrorKind::InvalidInput);
            assert!(
                err.to_string().contains("exceeds the"),
                "the error must name the limit; got: {err}",
            );
        }
    }

    /// A zero-length body must frame as `[key][0u32]` and read back as an empty
    /// record, on BOTH key branches.
    ///
    /// The failure this guards is not an error but a silent desync: the length
    /// prefix is what the reader advances by, so a framer or reader that
    /// mishandles a zero length leaves the stream misaligned and every
    /// *subsequent* record decodes as garbage. The neighbouring record proves
    /// alignment survived — asserting only that the empty record round-trips
    /// would miss exactly that.
    #[rstest::rstest]
    #[case::non_embedded_key(false)]
    #[case::embedded_key(true)]
    fn zero_length_record_body_frames_and_keeps_the_stream_aligned(#[case] embedded: bool) {
        let mut block = Vec::new();
        if embedded {
            let key = crate::keys::RawCoordinateKey { sort_key: 7 };
            frame_keyed_record_into(&mut block, &key, &[]).unwrap();
            // `[u32 LE 0]` and nothing else — no key prefix, no body.
            assert_eq!(block, 0u32.to_le_bytes().to_vec(), "embedded zero-length framing");
            frame_keyed_record_into(&mut block, &key, &[0xAB, 0xCD]).unwrap();
            let mut expected = 0u32.to_le_bytes().to_vec();
            expected.extend_from_slice(&2u32.to_le_bytes());
            expected.extend_from_slice(&[0xAB, 0xCD]);
            assert_eq!(block, expected, "the next record must start right after the empty one");
        } else {
            let key = make_key(3);
            frame_keyed_record_into(&mut block, &key, &[]).unwrap();
            let mut key_bytes = Vec::new();
            key.write_to(&mut key_bytes).unwrap();
            let mut expected = key_bytes.clone();
            expected.extend_from_slice(&0u32.to_le_bytes());
            assert_eq!(block, expected, "non-embedded zero-length framing");
            frame_keyed_record_into(&mut block, &key, &[0xAB, 0xCD]).unwrap();
            expected.extend_from_slice(&key_bytes);
            expected.extend_from_slice(&2u32.to_le_bytes());
            expected.extend_from_slice(&[0xAB, 0xCD]);
            assert_eq!(block, expected, "the next record must start right after the empty one");
        }
    }

    /// An embedded key (`RawCoordinateKey`, `EMBEDDED_IN_RECORD == true`) writes
    /// no key prefix — the framed block is exactly `[u32 LE len][record]`. Locks
    /// the embedded layout at the kernel boundary (the streaming-writer tests use
    /// the non-embedded `TemplateKey`, which exercises the key-prefix branch).
    #[test]
    fn frame_omits_key_prefix_for_embedded_keys() {
        let mut block = Vec::new();
        let key = crate::keys::RawCoordinateKey { sort_key: 0x0102_0304_0506_0708 };
        frame_keyed_record_into(&mut block, &key, &[0xAB, 0xCD]).unwrap();
        // Just the u32 LE length (2) + the 2 record bytes — no key prefix.
        assert_eq!(block, [2, 0, 0, 0, 0xAB, 0xCD]);
    }

    /// Write a spill file via the block kernel: frame records into raw blocks
    /// (cutting a new block every `records_per_block` records, to force many
    /// independent compressed units), compress each block, and bracket with the
    /// codec magic/trailer.
    fn write_via_kernel(
        path: &Path,
        codec: SpillCodec,
        compression: u32,
        records: &[(TemplateKey, Vec<u8>)],
        records_per_block: usize,
    ) {
        let mut comp = SpillBlockCompressor::new(codec, compression).unwrap();
        let mut out = BufWriter::with_capacity(256 * 1024, File::create(path).unwrap());
        out.write_all(spill_magic(codec)).unwrap();
        let mut block = Vec::new();
        let mut in_block = 0usize;
        let flush =
            |block: &mut Vec<u8>, comp: &mut SpillBlockCompressor, out: &mut BufWriter<File>| {
                out.write_all(&comp.compress_block(block).unwrap()).unwrap();
                block.clear();
            };
        for (k, r) in records {
            frame_keyed_record_into(&mut block, k, r).unwrap();
            in_block += 1;
            if in_block >= records_per_block {
                flush(&mut block, &mut comp, &mut out);
                in_block = 0;
            }
        }
        flush(&mut block, &mut comp, &mut out);
        out.write_all(spill_trailer(codec)).unwrap();
        out.flush().unwrap();
    }

    fn read_back(path: &Path) -> Vec<(TemplateKey, Vec<u8>)> {
        let mut reader =
            GenericKeyedChunkReader::<TemplateKey>::open(path, None).expect("open reader");
        let mut buf = Vec::new();
        let mut out = Vec::new();
        while let Some(key) = reader.next_record(&mut buf).expect("read record") {
            out.push((key, buf.clone()));
        }
        out
    }

    /// Kernel-built files read back to the original records for every codec/level
    /// — including a small block size that forces many independent compressed
    /// units, proving block-independent compression is reader-compatible.
    #[test]
    fn kernel_blocks_read_back_records() {
        let dir = TempDir::new().unwrap();
        let records = sample_records(400);
        for (codec, level, name) in
            [(SpillCodec::Zstd, 3, "z3"), (SpillCodec::Bgzf, 1, "b1"), (SpillCodec::Bgzf, 0, "b0")]
        {
            let path = dir.path().join(format!("kernel_{name}.keyed"));
            write_via_kernel(&path, codec, level, &records, 7);
            assert_eq!(read_back(&path), records, "kernel round-trip mismatch for {name}");
        }
    }

    /// Kernel output matches the streaming `write_sorted_chunk` oracle, read back
    /// through the production reader — locks block-granular ≡ streaming.
    #[test]
    fn kernel_matches_streaming_oracle_read_back() {
        let dir = TempDir::new().unwrap();
        let records = sample_records(500);
        let keyed: Vec<(TemplateKey, fgumi_raw_bam::RawRecord)> =
            records.iter().map(|(k, r)| (*k, fgumi_raw_bam::RawRecord::from(r.clone()))).collect();
        for (codec, level) in [(SpillCodec::Zstd, 3), (SpillCodec::Bgzf, 1)] {
            let kernel_path = dir.path().join("kernel.keyed");
            let oracle_path = dir.path().join("oracle.keyed");
            write_via_kernel(&kernel_path, codec, level, &records, 11);
            crate::write_sorted_chunk(&oracle_path, codec, level, &keyed).unwrap();
            assert_eq!(
                read_back(&kernel_path),
                read_back(&oracle_path),
                "kernel vs streaming-oracle read-back differ ({codec:?})"
            );
        }
    }

    /// Different block boundaries produce files that read back identically — the
    /// invariant that lets the (Serial) serialize step choose boundaries freely.
    #[test]
    fn reblocking_is_boundary_independent() {
        let dir = TempDir::new().unwrap();
        let records = sample_records(300);
        for codec in [SpillCodec::Zstd, SpillCodec::Bgzf] {
            let a = dir.path().join("a.keyed");
            let b = dir.path().join("b.keyed");
            write_via_kernel(&a, codec, 1, &records, 3);
            write_via_kernel(&b, codec, 1, &records, 64);
            assert_eq!(read_back(&a), read_back(&b), "re-blocking changed readback ({codec:?})");
        }
    }
}

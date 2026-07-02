//! Synchronous keyed-chunk spill writer (inline compression, no worker pool).
//!
//! The P6 Phase-1 decomposition runs the spill **compress** as a `Parallel`
//! `CompressSpill` step on the framework work-stealing pool, replacing the async
//! [`PooledChunkWriter`] → `SortWorkerPool` path so the sort no longer needs a
//! second, private thread pool. Each `CompressSpill` `try_run` already executes
//! on a framework worker, so the per-chunk write compresses **inline** on that
//! thread.
//!
//! The on-disk format is byte-compatible with [`PooledChunkWriter`], so
//! `SortSpillDecompress` (the Phase-2 streaming reader) and
//! `GenericKeyedChunkReader` (the merge reader) read these files unchanged:
//!
//! - **bgzf** (any `compression`, including `0`): framed BGZF blocks (header +
//!   deflate + footer) via [`InlineBgzfCompressor`] then a trailing `BGZF_EOF`,
//!   exactly like `PooledChunkWriter`'s bgzf path. Level 0 still produces *framed*
//!   (stored) BGZF blocks — NOT raw bytes — because the bgzf reader requires the
//!   `0x1f 0x8b` block magic; an unframed raw spill would fail the reader's
//!   magic check.
//! - **zstd** (the production default): [`ZSPILL_MAGIC`] then
//!   `[u32 LE frame-len][zstd frame]` per ≤`BGZF_MAX_BLOCK_SIZE` raw block, with
//!   no trailing marker — mirroring `PooledChunkWriter`'s zstd path and the
//!   `ZspillStreamReader` format.
//!
//! In both cases block/frame *boundaries* may differ from the pooled writer, but
//! the decompressed byte stream (and therefore every record read back) is
//! identical, because the reader streams across boundaries.
//!
//! [`PooledChunkWriter`]: crate::pooled_chunk_writer::PooledChunkWriter

use std::fs::File;
use std::io::{BufWriter, Write};
use std::marker::PhantomData;
use std::path::Path;

use anyhow::Result;
use fgumi_bgzf::writer::InlineBgzfCompressor;
use fgumi_bgzf::{BGZF_EOF, BGZF_MAX_BLOCK_SIZE};
use fgumi_raw_bam::RawRecord;
use zstd::bulk::Compressor as ZstdCompressor;

use crate::codec::{SpillCodec, ZSPILL_MAGIC};
use crate::keys::RawSortKey;

/// Compress and write a fully-sorted in-memory chunk to `path` as a keyed spill
/// file, inline on the calling thread (no `SortWorkerPool`).
///
/// This is the single-chunk entry point the P6 `CompressSpill` step calls per
/// sorted chunk: it creates a [`SyncSpillWriter`] for `codec`/`compression`,
/// writes every `(key, record)` pair in `records` order, and closes the file.
/// The on-disk format is byte-compatible with `PooledChunkWriter`, so
/// [`open_spill_slot`](crate::open_spill_slot) → `SortSpillDecompress` / the
/// merge reader consume it unchanged.
///
/// # Errors
///
/// Returns an error if the file cannot be created, the zstd compressor cannot be
/// initialized, or any record write/flush fails.
pub fn write_sorted_chunk<K: RawSortKey>(
    path: &Path,
    codec: SpillCodec,
    compression: u32,
    records: &[(K, RawRecord)],
) -> Result<()> {
    let mut writer = SyncSpillWriter::<K>::create(path, codec, compression)?;
    for (key, record) in records {
        writer.write_record(key, record.as_ref())?;
    }
    writer.finish()
}

/// Write an already-sorted [`InMemoryChunk`](crate::InMemoryChunk) to a spill
/// file — the zero-copy analogue of [`write_sorted_chunk`]. The chunk's records
/// share an `Arc<SegmentedBuf>` backing store, so this iterates by index and
/// writes each record's bytes directly (no owned `RawRecord`s), avoiding the
/// per-record copy the buffer chain previously paid to materialise
/// `Vec<(K, RawRecord)>`.
///
/// # Errors
///
/// Propagates I/O / compression errors from the underlying writer.
pub fn write_sorted_chunk_inmem<K: RawSortKey>(
    path: &Path,
    codec: SpillCodec,
    compression: u32,
    chunk: &crate::InMemoryChunk<K>,
) -> Result<()> {
    let mut writer = SyncSpillWriter::<K>::create(path, codec, compression)?;
    for i in 0..chunk.len() {
        writer.write_record(chunk.key_at(i), chunk.record_bytes(i))?;
    }
    writer.finish()
}

/// Synchronous keyed-chunk spill writer that compresses inline on the calling
/// (framework-worker) thread — no `SortWorkerPool`.
pub(crate) enum SyncSpillWriter<K: RawSortKey> {
    /// bgzf: framed BGZF blocks (any level, including a stored level-0 block) +
    /// trailing `BGZF_EOF`.
    Bgzf(BgzfSpillWriter<K>),
    /// zstd: inline-compressed `[u32 len][frame]` blocks after `ZSPILL_MAGIC`.
    Zstd(ZstdSpillWriter<K>),
}

impl<K: RawSortKey> SyncSpillWriter<K> {
    /// Create a writer for `path` using `codec` at `compression` level.
    ///
    /// For zstd, `compression` is the zstd level (must be ≥ 1 — level 0 is
    /// rejected up front by `SortOptions::validate`, since zstd has no
    /// uncompressed mode). For bgzf, `compression == 0` writes *framed* stored
    /// (uncompressed) BGZF blocks and `> 0` writes deflate-compressed BGZF blocks
    /// at that level — in both cases valid, reader-consumable BGZF.
    ///
    /// # Errors
    ///
    /// Returns an error if the output file cannot be created or the zstd
    /// compressor cannot be initialized.
    pub(crate) fn create(path: &Path, codec: SpillCodec, compression: u32) -> Result<Self> {
        match codec {
            SpillCodec::Bgzf => Ok(Self::Bgzf(BgzfSpillWriter::create(path, compression)?)),
            SpillCodec::Zstd => Ok(Self::Zstd(ZstdSpillWriter::create(path, compression)?)),
        }
    }

    /// Write one keyed record in the spill frame format.
    ///
    /// # Errors
    ///
    /// Returns an error if key serialization or the underlying write fails.
    pub(crate) fn write_record(&mut self, key: &K, record: &[u8]) -> Result<()> {
        match self {
            Self::Bgzf(w) => w.write_record(key, record),
            Self::Zstd(w) => w.write_record(key, record),
        }
    }

    /// Flush and close the chunk file.
    ///
    /// # Errors
    ///
    /// Returns an error if a final flush/compress fails.
    pub(crate) fn finish(self) -> Result<()> {
        match self {
            Self::Bgzf(w) => w.finish(),
            Self::Zstd(w) => w.finish(),
        }
    }
}

/// Inline BGZF spill writer: framed BGZF blocks (via [`InlineBgzfCompressor`])
/// followed by a trailing `BGZF_EOF`, matching `PooledChunkWriter`'s bgzf output.
/// Records are framed `[key.write_to() if !EMBEDDED][u32 LE record-len][record]`
/// into the compressor's byte stream, identical to the zstd arm.
pub(crate) struct BgzfSpillWriter<K: RawSortKey> {
    writer: BufWriter<File>,
    compressor: InlineBgzfCompressor,
    _marker: PhantomData<K>,
}

impl<K: RawSortKey> BgzfSpillWriter<K> {
    fn create(path: &Path, level: u32) -> Result<Self> {
        let file = File::create(path)?;
        let writer = BufWriter::with_capacity(256 * 1024, file);
        Ok(Self { writer, compressor: InlineBgzfCompressor::new(level), _marker: PhantomData })
    }

    /// Drain any completed (full-block) compressed output to disk, bounding the
    /// compressor's retained-block memory during a large chunk write.
    fn drain_blocks(&mut self) -> Result<()> {
        for block in self.compressor.take_blocks() {
            self.writer.write_all(&block.data)?;
        }
        Ok(())
    }

    fn write_record(&mut self, key: &K, record: &[u8]) -> Result<()> {
        // Validate the record length BEFORE any compressor write, so an oversized
        // record fails loud without leaving partial (key) bytes in the stream.
        let record_len = u32::try_from(record.len())
            .map_err(|_| anyhow::anyhow!("BAM record too large ({} bytes)", record.len()))?;
        if !K::EMBEDDED_IN_RECORD {
            let mut key_bytes = Vec::new();
            key.write_to(&mut key_bytes)?;
            self.compressor.write_all(&key_bytes)?;
        }
        self.compressor.write_all(&record_len.to_le_bytes())?;
        // Drain between block-sized chunks so retained compressed-block memory
        // stays bounded even for a very large (long-read) record, rather than
        // scaling with the record size.
        for chunk in record.chunks(BGZF_MAX_BLOCK_SIZE) {
            self.compressor.write_all(chunk)?;
            self.drain_blocks()?;
        }
        // Catch the key/length bytes when `record` is empty (the loop is a no-op).
        self.drain_blocks()
    }

    fn finish(mut self) -> Result<()> {
        self.compressor.flush()?;
        self.drain_blocks()?;
        // BGZF stream terminator (the empty-block EOF marker), as the pooled
        // writer and noodles bgzf writer both emit.
        self.writer.write_all(&BGZF_EOF)?;
        self.writer.flush()?;
        Ok(())
    }
}

/// Inline zstd spill writer: `ZSPILL_MAGIC` then length-prefixed zstd frames,
/// one frame per ≤`BGZF_MAX_BLOCK_SIZE` raw block.
pub(crate) struct ZstdSpillWriter<K: RawSortKey> {
    writer: BufWriter<File>,
    compressor: ZstdCompressor<'static>,
    /// Raw (uncompressed) staging block; flushed as a zstd frame at
    /// `BGZF_MAX_BLOCK_SIZE`.
    block: Vec<u8>,
    _marker: PhantomData<K>,
}

impl<K: RawSortKey> ZstdSpillWriter<K> {
    fn create(path: &Path, level: u32) -> Result<Self> {
        let file = File::create(path)?;
        let mut writer = BufWriter::with_capacity(256 * 1024, file);
        writer.write_all(&ZSPILL_MAGIC)?;
        #[allow(clippy::cast_possible_wrap)]
        let compressor = ZstdCompressor::new(level as i32)
            .map_err(|e| anyhow::anyhow!("zstd compressor init (level {level}): {e}"))?;
        Ok(Self {
            writer,
            compressor,
            block: Vec::with_capacity(BGZF_MAX_BLOCK_SIZE + 1024),
            _marker: PhantomData,
        })
    }

    /// Compress the staged block as one zstd frame and write `[u32 len][frame]`.
    /// No-op when the block is empty (no empty frames in the stream).
    fn flush_block(&mut self) -> Result<()> {
        if self.block.is_empty() {
            return Ok(());
        }
        let frame = self
            .compressor
            .compress(&self.block)
            .map_err(|e| anyhow::anyhow!("zstd compress: {e}"))?;
        let frame_len = u32::try_from(frame.len())
            .map_err(|_| anyhow::anyhow!("zstd frame larger than 4 GiB cannot fit a u32 prefix"))?;
        self.writer.write_all(&frame_len.to_le_bytes())?;
        self.writer.write_all(&frame)?;
        self.block.clear();
        Ok(())
    }

    /// Append `data` to the staging block, flushing a frame whenever the block
    /// fills. A record (or key) larger than one block therefore spans frames —
    /// safe because the reader streams across frame boundaries.
    fn append(&mut self, mut data: &[u8]) -> Result<()> {
        while !data.is_empty() {
            let space = BGZF_MAX_BLOCK_SIZE.saturating_sub(self.block.len());
            let n = data.len().min(space);
            self.block.extend_from_slice(&data[..n]);
            data = &data[n..];
            if self.block.len() >= BGZF_MAX_BLOCK_SIZE {
                self.flush_block()?;
            }
        }
        Ok(())
    }

    fn write_record(&mut self, key: &K, record: &[u8]) -> Result<()> {
        // Validate the record length BEFORE any append, so an oversized record
        // fails loud without leaving partial (key) bytes in the staging block.
        let record_len = u32::try_from(record.len())
            .map_err(|_| anyhow::anyhow!("BAM record too large ({} bytes)", record.len()))?;
        if !K::EMBEDDED_IN_RECORD {
            // Keys are small and fixed-size; a transient buffer keeps the
            // borrow simple. The record path below dominates the spill cost.
            let mut key_bytes = Vec::new();
            key.write_to(&mut key_bytes)?;
            self.append(&key_bytes)?;
        }
        self.append(&record_len.to_le_bytes())?;
        self.append(record)?;
        Ok(())
    }

    fn finish(mut self) -> Result<()> {
        self.flush_block()?;
        self.writer.flush()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::external::GenericKeyedChunkReader;
    use crate::inline::TemplateKey;
    use crate::pooled_chunk_writer::PooledChunkWriter;
    use crate::worker_pool::SortWorkerPool;
    use std::sync::Arc;
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

    /// The sync zstd writer round-trips every record back through the same reader
    /// the production merge/Phase-2 path uses, and the file carries `ZSPILL_MAGIC`.
    #[test]
    fn sync_zstd_round_trips() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("sync_zstd.keyed");
        let records = sample_records(300);

        let mut w = SyncSpillWriter::<TemplateKey>::create(&path, SpillCodec::Zstd, 3).unwrap();
        for (k, r) in &records {
            w.write_record(k, r).unwrap();
        }
        w.finish().unwrap();

        let bytes = std::fs::read(&path).unwrap();
        assert_eq!(&bytes[..ZSPILL_MAGIC.len()], &ZSPILL_MAGIC[..], "missing ZSPILL_MAGIC");
        assert_eq!(read_back(&path), records, "sync zstd round-trip mismatch");
    }

    /// Cross-check: the sync zstd writer and the pooled (async) zstd writer
    /// produce files that read back to the *same* records — proving the formats
    /// are interchangeable for `SortSpillDecompress` / the merge reader.
    #[test]
    fn sync_zstd_matches_pooled_zstd_read_back() {
        let dir = TempDir::new().unwrap();
        let records = sample_records(500);

        let sync_path = dir.path().join("sync.keyed");
        let mut w =
            SyncSpillWriter::<TemplateKey>::create(&sync_path, SpillCodec::Zstd, 3).unwrap();
        for (k, r) in &records {
            w.write_record(k, r).unwrap();
        }
        w.finish().unwrap();

        let pooled_path = dir.path().join("pooled.keyed");
        let pool = Arc::new(SortWorkerPool::new(2, 3, 6, SpillCodec::Zstd));
        {
            let mut pw = PooledChunkWriter::<TemplateKey>::new(
                Arc::clone(&pool),
                &pooled_path,
                SpillCodec::Zstd,
            )
            .unwrap();
            for (k, r) in &records {
                pw.write_record(k, r).unwrap();
            }
            pw.finish().unwrap();
        }
        if let Ok(p) = Arc::try_unwrap(pool) {
            p.shutdown();
        }

        assert_eq!(
            read_back(&sync_path),
            read_back(&pooled_path),
            "sync vs pooled read-back differ"
        );
    }

    /// Records larger than one block must span frames and still round-trip.
    #[test]
    fn sync_zstd_spans_frames_for_large_records() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("large.keyed");
        // Two records each ~3x the block size — forces multi-frame spanning.
        let big = vec![0xABu8; BGZF_MAX_BLOCK_SIZE * 3 + 17];
        let records = vec![(make_key(1), big.clone()), (make_key(2), big)];

        let mut w = SyncSpillWriter::<TemplateKey>::create(&path, SpillCodec::Zstd, 1).unwrap();
        for (k, r) in &records {
            w.write_record(k, r).unwrap();
        }
        w.finish().unwrap();

        assert_eq!(read_back(&path), records, "large-record frame spanning round-trip mismatch");
    }

    /// Uncompressed bgzf (level 0) must write *framed* (stored) BGZF blocks —
    /// starting with the `0x1f 0x8b` gzip magic — NOT raw bytes, or the bgzf
    /// reader's magic check rejects the spill. Regression for the
    /// `--temp-compression 0 --temp-codec bgzf` path. Also cross-checks that the
    /// records read back match the pooled writer at level 0.
    #[test]
    fn sync_bgzf_level0_is_framed_and_matches_pooled() {
        let dir = TempDir::new().unwrap();
        let records = sample_records(300);

        let sync_path = dir.path().join("sync_l0.keyed");
        let mut w =
            SyncSpillWriter::<TemplateKey>::create(&sync_path, SpillCodec::Bgzf, 0).unwrap();
        for (k, r) in &records {
            w.write_record(k, r).unwrap();
        }
        w.finish().unwrap();

        let bytes = std::fs::read(&sync_path).unwrap();
        assert_eq!(
            &bytes[..2],
            &[0x1f, 0x8b],
            "level-0 bgzf spill must be framed BGZF (gzip magic), not raw"
        );
        assert_eq!(read_back(&sync_path), records, "level-0 bgzf round-trip mismatch");

        // Cross-check vs the pooled writer at level 0 (its pool is built with
        // temp_compression = 0).
        let pooled_path = dir.path().join("pooled_l0.keyed");
        let pool = Arc::new(SortWorkerPool::new(2, 0, 6, SpillCodec::Bgzf));
        {
            let mut pw = PooledChunkWriter::<TemplateKey>::new(
                Arc::clone(&pool),
                &pooled_path,
                SpillCodec::Bgzf,
            )
            .unwrap();
            for (k, r) in &records {
                pw.write_record(k, r).unwrap();
            }
            pw.finish().unwrap();
        }
        if let Ok(p) = Arc::try_unwrap(pool) {
            p.shutdown();
        }
        assert_eq!(
            read_back(&sync_path),
            read_back(&pooled_path),
            "level-0 bgzf: sync vs pooled read-back differ"
        );
    }

    /// `write_sorted_chunk` (the per-chunk entry the `CompressSpill` step calls)
    /// round-trips every record through the merge reader, and `open_spill_slot`
    /// opens the result with the requested `file_id` and the detected codec.
    #[test]
    fn write_sorted_chunk_round_trips_and_open_spill_slot_sets_file_id() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("chunk_0007.keyed");
        let records = sample_records(400);
        // `write_sorted_chunk` takes `(K, RawRecord)`; `read_back` yields
        // `(K, Vec<u8>)`. Build the keyed form once and compare read-back vs the
        // original `Vec<u8>` payloads.
        let keyed: Vec<(TemplateKey, RawRecord)> =
            records.iter().map(|(k, r)| (*k, RawRecord::from(r.clone()))).collect();

        crate::write_sorted_chunk(&path, SpillCodec::Zstd, 3, &keyed).unwrap();
        assert_eq!(read_back(&path), records, "write_sorted_chunk round-trip mismatch");

        let slot = crate::open_spill_slot(&path, 7).expect("open spill slot");
        assert_eq!(slot.file_id, 7, "open_spill_slot must honor the requested file_id");
        assert_eq!(slot.codec, SpillCodec::Zstd, "codec must be detected from the file magic");

        // The bgzf arm: codec detection must fall back to Bgzf (no ZSPILL magic).
        let bgzf_path = dir.path().join("chunk_0008.keyed");
        crate::write_sorted_chunk(&bgzf_path, SpillCodec::Bgzf, 1, &keyed).unwrap();
        let bgzf_slot = crate::open_spill_slot(&bgzf_path, 8).expect("open bgzf spill slot");
        assert_eq!(bgzf_slot.file_id, 8);
        assert_eq!(bgzf_slot.codec, SpillCodec::Bgzf);
    }

    /// The bgzf arm delegates to the proven `GenericKeyedChunkWriter`; confirm
    /// the unified type round-trips there too.
    #[test]
    fn sync_bgzf_round_trips() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("sync_bgzf.keyed");
        let records = sample_records(200);

        let mut w = SyncSpillWriter::<TemplateKey>::create(&path, SpillCodec::Bgzf, 1).unwrap();
        for (k, r) in &records {
            w.write_record(k, r).unwrap();
        }
        w.finish().unwrap();

        assert_eq!(read_back(&path), records, "sync bgzf round-trip mismatch");
    }
}

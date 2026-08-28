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

use std::io::{BufWriter, Write};
use std::marker::PhantomData;
use std::path::{Path, PathBuf};

use anyhow::Result;
use fgumi_bgzf::writer::InlineBgzfCompressor;
use fgumi_bgzf::{BGZF_EOF, BGZF_MAX_BLOCK_SIZE};
use fgumi_raw_bam::RawRecord;
use tempfile::NamedTempFile;
use zstd::bulk::Compressor as ZstdCompressor;

use crate::codec::{SpillCodec, ZSPILL_MAGIC};
use crate::keys::RawSortKey;

/// Open a staging temp file in the same directory as `path`.
///
/// The spill is written to this temp file and only renamed into place by
/// [`persist_spill`] once `finish()` succeeds, so a crashed or errored write
/// never leaves a partial file at `path` (a `ZspillStreamReader` would otherwise
/// accept a partial file that happens to end on a frame boundary as a valid,
/// shorter spill and silently drop the trailing records). The temp file shares
/// `path`'s directory so the rename is same-filesystem and therefore atomic, and
/// `NamedTempFile` removes it on drop if `finish()` is never reached.
fn stage_spill(path: &Path) -> Result<(BufWriter<NamedTempFile>, PathBuf)> {
    let dir = match path.parent() {
        Some(parent) if !parent.as_os_str().is_empty() => parent,
        _ => Path::new("."),
    };
    let temp = NamedTempFile::new_in(dir)?;
    Ok((BufWriter::with_capacity(256 * 1024, temp), path.to_path_buf()))
}

/// Atomically publish a finished spill temp file at `final_path`.
///
/// Uses `persist_noclobber`, so it fails closed if a file already exists at
/// `final_path` rather than truncating a prior spill — the same "never reuse an
/// existing path" guarantee `SpillWrite::open_file`'s `create_new` provides,
/// preserved now that the destination is opened by rename rather than directly.
fn persist_spill(writer: BufWriter<NamedTempFile>, final_path: &Path) -> Result<()> {
    let temp = writer.into_inner().map_err(|e| anyhow::anyhow!("flushing spill temp file: {e}"))?;
    temp.persist_noclobber(final_path).map_err(|e| {
        anyhow::anyhow!("publishing spill file {}: {}", final_path.display(), e.error)
    })?;
    Ok(())
}

/// Compress and write a fully-sorted in-memory chunk to `path` as a keyed spill
/// file, inline on the calling thread (no `SortWorkerPool`).
///
/// This is the single-chunk entry point the P6 `CompressSpill` step calls per
/// sorted chunk: it creates a `SyncSpillWriter` for `codec`/`compression`,
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
    writer: BufWriter<NamedTempFile>,
    /// Final destination the temp file is renamed to by `finish()`.
    final_path: PathBuf,
    compressor: InlineBgzfCompressor,
    /// Reused staging buffer for serialized keys (non-embedded keys only).
    ///
    /// Kept on the writer rather than allocated per record: this is the
    /// queryname and template-coordinate spill path, so a per-record `Vec`
    /// costs one malloc/free for every record in the sort.
    key_scratch: Vec<u8>,
    _marker: PhantomData<K>,
}

impl<K: RawSortKey> BgzfSpillWriter<K> {
    fn create(path: &Path, level: u32) -> Result<Self> {
        // Stage into a temp file and rename into place only after `finish()`
        // succeeds (see `stage_spill` / `persist_spill`): a partial write never
        // reaches `path`, and the final rename still fails closed on a stale
        // path, matching `SpillWrite::open_file`.
        let (writer, final_path) = stage_spill(path)?;
        Ok(Self {
            writer,
            final_path,
            compressor: InlineBgzfCompressor::new(level),
            key_scratch: Vec::new(),
            _marker: PhantomData,
        })
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
            // `compressor` and `key_scratch` are disjoint fields, so the scratch
            // can be filled and written without moving it out of `self`.
            self.key_scratch.clear();
            key.write_to(&mut self.key_scratch)?;
            self.compressor.write_all(&self.key_scratch)?;
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
        persist_spill(self.writer, &self.final_path)
    }
}

/// Inline zstd spill writer: `ZSPILL_MAGIC` then length-prefixed zstd frames,
/// one frame per ≤`BGZF_MAX_BLOCK_SIZE` raw block.
pub(crate) struct ZstdSpillWriter<K: RawSortKey> {
    writer: BufWriter<NamedTempFile>,
    /// Final destination the temp file is renamed to by `finish()`.
    final_path: PathBuf,
    compressor: ZstdCompressor<'static>,
    /// Raw (uncompressed) staging block; flushed as a zstd frame at
    /// `BGZF_MAX_BLOCK_SIZE`.
    block: Vec<u8>,
    /// Reused destination for each compressed zstd frame, so `flush_block` does
    /// not allocate a fresh `Vec` per frame (`compress` returns an owned one).
    frame_buf: Vec<u8>,
    /// Reused staging buffer for serialized keys — see
    /// [`BgzfSpillWriter::key_scratch`].
    key_scratch: Vec<u8>,
    _marker: PhantomData<K>,
}

impl<K: RawSortKey> ZstdSpillWriter<K> {
    fn create(path: &Path, level: u32) -> Result<Self> {
        // See `BgzfSpillWriter::create`: staged temp file, atomic rename on
        // finish, fail-closed (no clobber) on the destination.
        let (mut writer, final_path) = stage_spill(path)?;
        writer.write_all(&ZSPILL_MAGIC)?;
        #[allow(clippy::cast_possible_wrap)]
        let compressor = ZstdCompressor::new(level as i32)
            .map_err(|e| anyhow::anyhow!("zstd compressor init (level {level}): {e}"))?;
        Ok(Self {
            writer,
            final_path,
            compressor,
            block: Vec::with_capacity(BGZF_MAX_BLOCK_SIZE + 1024),
            frame_buf: Vec::new(),
            key_scratch: Vec::new(),
            _marker: PhantomData,
        })
    }

    /// Compress the staged block as one zstd frame and write `[u32 len][frame]`.
    /// No-op when the block is empty (no empty frames in the stream).
    fn flush_block(&mut self) -> Result<()> {
        if self.block.is_empty() {
            return Ok(());
        }
        // Reuse `frame_buf` across frames rather than letting `compress` return a
        // freshly-allocated `Vec` for every non-empty spill frame.
        self.frame_buf.clear();
        self.frame_buf.reserve(zstd::zstd_safe::compress_bound(self.block.len()));
        self.compressor
            .compress_to_buffer(&self.block, &mut self.frame_buf)
            .map_err(|e| anyhow::anyhow!("zstd compress: {e}"))?;
        let frame_len = u32::try_from(self.frame_buf.len())
            .map_err(|_| anyhow::anyhow!("zstd frame larger than 4 GiB cannot fit a u32 prefix"))?;
        self.writer.write_all(&frame_len.to_le_bytes())?;
        self.writer.write_all(&self.frame_buf)?;
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
            // `append` takes `&mut self`, so unlike the bgzf arm the scratch
            // cannot stay borrowed from `self` across the call — move it out and
            // put it back, which reuses the allocation across records all the
            // same. An error here abandons the whole spill, so not restoring the
            // buffer on that path costs nothing.
            let mut key_bytes = std::mem::take(&mut self.key_scratch);
            key_bytes.clear();
            key.write_to(&mut key_bytes)?;
            self.append(&key_bytes)?;
            self.key_scratch = key_bytes;
        }
        self.append(&record_len.to_le_bytes())?;
        self.append(record)?;
        Ok(())
    }

    fn finish(mut self) -> Result<()> {
        self.flush_block()?;
        self.writer.flush()?;
        persist_spill(self.writer, &self.final_path)
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

    /// A minimal BAM fixed-block body with the coordinate fields set: `ref_id` at
    /// 0..4, `pos` at 4..8, `flag` at 14..16 — the three
    /// `RawCoordinateKey::extract_from_record` reads. 32 bytes is the fixed block,
    /// plus a tail byte to make each record distinguishable.
    fn coordinate_body(tid: i32, pos: i32, tail: u8) -> Vec<u8> {
        let mut b = vec![0u8; 32];
        b[0..4].copy_from_slice(&tid.to_le_bytes());
        b[4..8].copy_from_slice(&pos.to_le_bytes());
        b[14..16].copy_from_slice(&0u16.to_le_bytes()); // flags: mapped, forward
        b.push(tail);
        b
    }

    fn read_back_coordinate(path: &Path) -> Vec<Vec<u8>> {
        let mut reader = GenericKeyedChunkReader::<crate::keys::RawCoordinateKey>::open(path, None)
            .expect("open reader");
        let mut buf = Vec::new();
        let mut out = Vec::new();
        while reader.next_record(&mut buf).expect("read record").is_some() {
            out.push(buf.clone());
        }
        out
    }

    /// The embedded-key writer path, on both codecs.
    ///
    /// Every other writer test uses `TemplateKey`, whose `EMBEDDED_IN_RECORD` is
    /// `false` — so they all take the `if !K::EMBEDDED_IN_RECORD` branch and write
    /// a key prefix. The *skip* side of that branch, which the coordinate and
    /// queryname sorts use, was never exercised at this layer: those records carry
    /// their key inside the body and the spill frame is `[u32 len][record]` with no
    /// prefix at all. A writer that emitted a prefix anyway would desync the
    /// reader, and no existing test would notice.
    #[rstest::rstest]
    #[case::zstd(SpillCodec::Zstd)]
    #[case::bgzf(SpillCodec::Bgzf)]
    fn sync_writers_round_trip_embedded_keys(#[case] codec: SpillCodec) {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("embedded.keyed");
        // Ascending positions so the spill is a sorted run, as production writes.
        let bodies: Vec<Vec<u8>> =
            (0..16i32).map(|i| coordinate_body(0, i * 10, u8::try_from(i).unwrap())).collect();

        let mut w =
            SyncSpillWriter::<crate::keys::RawCoordinateKey>::create(&path, codec, 1).unwrap();
        for b in &bodies {
            let key = crate::keys::RawCoordinateKey::extract_from_record(b);
            w.write_record(&key, b).unwrap();
        }
        w.finish().unwrap();

        assert_eq!(
            read_back_coordinate(&path),
            bodies,
            "embedded-key records must round-trip with no key prefix ({codec:?})",
        );
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

        // Anchor one side to the original `records` first: a comparison of only
        // sync-vs-pooled would still pass if BOTH writers dropped every record
        // (or the reader returned nothing for both).
        let sync_back = read_back(&sync_path);
        assert_eq!(sync_back, records, "sync zstd writer must round-trip the original records");
        assert_eq!(sync_back, read_back(&pooled_path), "sync vs pooled read-back differ");
    }

    /// A zero-length record body must round-trip through the writer, on BOTH
    /// codecs, without desynchronising the records around it.
    ///
    /// This is the sibling of the kernel-level
    /// `zero_length_record_body_frames_and_keeps_the_stream_aligned`: that one
    /// pins `frame_keyed_record_into`, this one pins the streaming writers.
    ///
    /// The empty body is the degenerate shape of both arms' inner loops —
    /// `record.chunks(BGZF_MAX_BLOCK_SIZE)` and the zstd `append` both iterate
    /// zero times — so the record reaches the stream as `[key][0u32]` and
    /// nothing else, entirely via the surrounding flush paths. The empty record
    /// sits between two non-empty ones so a misframed zero length surfaces as
    /// its neighbours decoding wrong, rather than as one silently-absent record.
    ///
    /// Note on scope: the trailing `drain_blocks()` in the bgzf arm is NOT what
    /// makes this pass. Removing it keeps every test green, because `finish()`
    /// flushes the compressor before draining — that call bounds retained
    /// block memory between records, which is what its own comment claims, and
    /// is not load-bearing for correctness.
    #[rstest::rstest]
    #[case::zstd(SpillCodec::Zstd)]
    #[case::bgzf(SpillCodec::Bgzf)]
    fn sync_writers_round_trip_a_zero_length_record_body(#[case] codec: SpillCodec) {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("empty-body.keyed");
        let records: Vec<(TemplateKey, Vec<u8>)> = vec![
            (make_key(1), vec![0xAA; 64]),
            (make_key(2), Vec::new()),
            (make_key(3), vec![0xBB; 96]),
        ];

        let mut w = SyncSpillWriter::<TemplateKey>::create(&path, codec, 1).unwrap();
        for (k, r) in &records {
            w.write_record(k, r).unwrap();
        }
        w.finish().unwrap();

        assert_eq!(
            read_back(&path),
            records,
            "a zero-length body must round-trip and keep its neighbours aligned ({codec:?})",
        );
    }

    /// Records larger than one block must span blocks and still round-trip —
    /// on BOTH codecs.
    ///
    /// The two arms interleave differently and neither can vouch for the other:
    /// zstd fills a staging block and flushes a frame when it is full, while
    /// bgzf pushes `record.chunks(BGZF_MAX_BLOCK_SIZE)` into the compressor and
    /// drains completed blocks between chunks. That chunk/drain interleaving is
    /// exactly where a boundary defect would live, and it was untested — the
    /// case ran zstd only.
    #[rstest::rstest]
    #[case::zstd(SpillCodec::Zstd)]
    #[case::bgzf(SpillCodec::Bgzf)]
    fn sync_writers_span_blocks_for_large_records(#[case] codec: SpillCodec) {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("large.keyed");
        // Two records each ~3x the block size — forces multi-block spanning.
        // The `+ 17` keeps the tail from landing on a block boundary.
        let big = vec![0xABu8; BGZF_MAX_BLOCK_SIZE * 3 + 17];
        let records = vec![(make_key(1), big.clone()), (make_key(2), big)];

        let mut w = SyncSpillWriter::<TemplateKey>::create(&path, codec, 1).unwrap();
        for (k, r) in &records {
            w.write_record(k, r).unwrap();
        }
        w.finish().unwrap();

        assert_eq!(
            read_back(&path),
            records,
            "large-record block spanning round-trip mismatch for {codec:?}",
        );
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

    /// `write_sorted_chunk_inmem` must produce a byte-identical spill file to
    /// `write_sorted_chunk` for the same records.
    ///
    /// It is a separate public entry point that walks the chunk by index
    /// (`key_at(i)` / `record_bytes(i)`) instead of iterating owned pairs, and
    /// nothing covered it. An off-by-one or a key/record mismatch in that
    /// indexing still yields a spill file that parses cleanly, so the merge
    /// would emit wrong records rather than fail — the parity assertion is what
    /// makes that visible.
    #[test]
    fn write_sorted_chunk_inmem_matches_write_sorted_chunk() {
        let dir = TempDir::new().unwrap();
        let records = sample_records(400);

        let keyed: Vec<(TemplateKey, RawRecord)> =
            records.iter().map(|(k, r)| (*k, RawRecord::from(r.clone()))).collect();
        let owned: Vec<(TemplateKey, Vec<u8>)> =
            records.iter().map(|(k, r)| (*k, r.clone())).collect();
        let chunk = crate::InMemoryChunk::from_owned_records(owned);

        let from_pairs = dir.path().join("pairs.keyed");
        let from_chunk = dir.path().join("chunk.keyed");
        crate::write_sorted_chunk(&from_pairs, SpillCodec::Zstd, 3, &keyed).unwrap();
        crate::write_sorted_chunk_inmem(&from_chunk, SpillCodec::Zstd, 3, &chunk).unwrap();

        assert_eq!(read_back(&from_chunk), records, "in-mem writer round-trip mismatch");
        assert_eq!(
            std::fs::read(&from_chunk).unwrap(),
            std::fs::read(&from_pairs).unwrap(),
            "the two writers must produce byte-identical spill files",
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

    /// The bgzf arm frames records itself and drives `InlineBgzfCompressor`
    /// directly; confirm the unified type round-trips there too.
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

    /// A spill file is published atomically: it appears at its final path only
    /// after `finish()`, a writer dropped without `finish()` leaves nothing
    /// there, and `finish()` refuses to clobber an existing spill (preserving
    /// the fail-closed guarantee `create_new` gave), on BOTH codecs.
    ///
    /// The partial-file hazard this guards is real: `ZspillStreamReader` treats
    /// EOF after a complete frame as a clean end, so a partial file ending on a
    /// frame boundary would otherwise read back as a valid, shorter spill.
    #[rstest::rstest]
    #[case::zstd(SpillCodec::Zstd)]
    #[case::bgzf(SpillCodec::Bgzf)]
    fn spill_is_published_atomically_and_never_clobbers(#[case] codec: SpillCodec) {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("chunk_atomic.keyed");

        // A writer dropped before `finish()` must leave no file at the final
        // path (the temp file is removed on drop).
        {
            let mut w = SyncSpillWriter::<TemplateKey>::create(&path, codec, 1).unwrap();
            w.write_record(&make_key(1), &[0xAB; 64]).unwrap();
            // dropped here without `finish()`
        }
        assert!(
            !path.exists(),
            "an unfinished spill must not appear at its final path ({codec:?})"
        );

        // `finish()` publishes the file with exactly the written records.
        let records = sample_records(24);
        let mut w = SyncSpillWriter::<TemplateKey>::create(&path, codec, 1).unwrap();
        for (k, r) in &records {
            w.write_record(k, r).unwrap();
        }
        w.finish().unwrap();
        assert!(path.exists(), "a finished spill must appear at its final path ({codec:?})");
        assert_eq!(
            read_back(&path),
            records,
            "published spill must hold the written records ({codec:?})"
        );

        // A second writer targeting the same path fails closed at `finish()`
        // rather than truncating the existing spill.
        let mut w2 = SyncSpillWriter::<TemplateKey>::create(&path, codec, 1).unwrap();
        w2.write_record(&make_key(99), &[0xCD; 32]).unwrap();
        assert!(
            w2.finish().is_err(),
            "finishing onto an existing spill must fail closed, not clobber ({codec:?})",
        );
        assert_eq!(
            read_back(&path),
            records,
            "the pre-existing spill must be left intact ({codec:?})"
        );
    }
}

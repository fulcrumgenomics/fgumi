//! BAM writer factories and related types.
//!
//! This module provides writer factories for BAM output, including standard BAM writers,
//! raw-bytes BAM writers, an incremental-indexing BAM writer, and BAI sidecar
//! helpers (`bai_sidecar_path` for the naming convention, `write_bai_index` for
//! a prebuilt index, and `write_bai_sidecar` to index a finished BAM and write
//! its sidecar in one call).

use anyhow::{Context, Result};
use bgzf::CompressionLevel;
use noodles::bam::bai;
use noodles::bgzf::VirtualPosition;
use noodles::core::Position;
use noodles::sam::Header;
use noodles_bgzf::io::Writer as BgzfWriter;
use noodles_csi::binning_index::Indexer;
use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;
use noodles_csi::binning_index::index::reference_sequence::index::LinearIndex;
use std::collections::{HashMap, VecDeque};
use std::fs::File;
use std::hash::{BuildHasherDefault, Hasher};
use std::io::{self, Write};
use std::num::NonZero;
use std::path::{Path, PathBuf};

use crate::paths::is_stdout_path;
use crate::vendored::{BlockInfoRx, MultithreadedWriter, MultithreadedWriterBuilder};

/// Maximum uncompressed BGZF block size.
///
/// [`compute_end_vpos`] serves two writers, so this must match the block size
/// both of them flush at:
///
/// - The sort engine's pooled writer flushes at `fgumi_bgzf::BGZF_MAX_BLOCK_SIZE`,
///   which is itself this same `bgzf::BGZF_BLOCK_SIZE` — deriving the constant
///   here (rather than hardcoding 65280) means that path's virtual-offset math
///   and block layout cannot silently diverge.
/// - [`IndexingBamWriter`] flushes at the vendored `MultithreadedWriter`'s own
///   `BGZF_BLOCK_SIZE` (`crate::vendored::bgzf_multithreaded`), a deliberately
///   independent copy that keeps the vendored module self-contained. It is equal
///   to this value today; keep the two in step, since nothing enforces it.
const MAX_BLOCK_SIZE: usize = bgzf::BGZF_BLOCK_SIZE;

/// Fast, non-cryptographic hasher for the dense sequential `u64` block numbers
/// used as `block_positions` keys.
///
/// One multiplicative (Fibonacci) step spreads sequential keys across the full
/// 64-bit space — enough entropy for both hashbrown's bucket index and its
/// control byte — at the cost of a single multiply and no DoS-resistant mixing,
/// which this internal, non-adversarial map does not need. A standard hasher
/// (`ahash`, `FxHash`) would also serve, but each costs several
/// multiply-equivalents per `u64`; this map is hot in the merge, so the single
/// multiply is deliberate.
#[derive(Default)]
struct BlockHasher(u64);

impl Hasher for BlockHasher {
    #[inline]
    fn write_u64(&mut self, n: u64) {
        // 2^64 / golden ratio; odd, so the multiply is a bijection over u64.
        self.0 = (self.0 ^ n).wrapping_mul(0x9E37_79B9_7F4A_7C15);
    }

    #[inline]
    fn finish(&self) -> u64 {
        self.0
    }

    fn write(&mut self, bytes: &[u8]) {
        // `u64` keys hash via `write_u64`; this path is unused but required by
        // the trait. Fold bytes in (FNV-style) so any stray key still hashes.
        for &b in bytes {
            self.0 = (self.0 ^ u64::from(b)).wrapping_mul(0x0100_0000_01B3);
        }
    }
}

/// `block_number -> compressed start offset`, hashed with [`BlockHasher`].
type BlockPositions = HashMap<u64, u64, BuildHasherDefault<BlockHasher>>;

/// Resolve the end virtual position for a BAM record that begins in
/// `start_block` at offset `start_block_start` in the compressed stream.
///
/// `end_offset` is the uncompressed offset (relative to `start_block`) one byte
/// past the record's last byte. If the record fits in `start_block`, returns
/// `(start_block_start, end_offset)`. Otherwise walks forward through the
/// `block_positions` map `MAX_BLOCK_SIZE` bytes at a time until the remaining
/// overflow fits in a single block. Returns `None` if any required subsequent
/// block hasn't been completed yet (caller should defer the entry).
#[allow(clippy::cast_possible_truncation)]
fn compute_end_vpos(
    block_positions: &BlockPositions,
    start_block: u64,
    start_block_start: u64,
    end_offset: usize,
) -> Option<VirtualPosition> {
    if end_offset <= MAX_BLOCK_SIZE {
        return Some(
            VirtualPosition::try_from((start_block_start, end_offset as u16))
                .unwrap_or(VirtualPosition::MIN),
        );
    }

    let mut overflow = end_offset - MAX_BLOCK_SIZE;
    let mut block = start_block + 1;
    while overflow > MAX_BLOCK_SIZE {
        overflow -= MAX_BLOCK_SIZE;
        block += 1;
    }

    let &block_start = block_positions.get(&block)?;
    Some(VirtualPosition::try_from((block_start, overflow as u16)).unwrap_or(VirtualPosition::MIN))
}

/// Opaque wrapper over either a single- or multi-threaded BGZF writer.
///
/// Uses `Box<dyn Write + Send>` to support both file and stdout output. The
/// internal representation is intentionally private so that the choice of
/// underlying multi-threaded BGZF writer is not part of the public API.
pub struct BgzfWriterEnum {
    inner: BgzfWriterImpl,
}

enum BgzfWriterImpl {
    SingleThreaded(BgzfWriter<Box<dyn Write + Send>>),
    MultiThreaded(MultithreadedWriter<Box<dyn Write + Send>>),
}

impl BgzfWriterEnum {
    pub(crate) fn single_threaded(writer: BgzfWriter<Box<dyn Write + Send>>) -> Self {
        Self { inner: BgzfWriterImpl::SingleThreaded(writer) }
    }

    pub(crate) fn multi_threaded(writer: MultithreadedWriter<Box<dyn Write + Send>>) -> Self {
        Self { inner: BgzfWriterImpl::MultiThreaded(writer) }
    }
}

impl Write for BgzfWriterEnum {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match &mut self.inner {
            BgzfWriterImpl::SingleThreaded(w) => w.write(buf),
            BgzfWriterImpl::MultiThreaded(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match &mut self.inner {
            BgzfWriterImpl::SingleThreaded(w) => w.flush(),
            BgzfWriterImpl::MultiThreaded(w) => w.flush(),
        }
    }
}

impl BgzfWriterEnum {
    /// Finish writing and close the writer properly.
    /// This is especially important for multi-threaded writers to ensure
    /// all data is flushed and the EOF marker is written.
    ///
    /// # Errors
    /// Returns an error if flushing or finalizing the writer fails.
    pub fn finish(self) -> io::Result<()> {
        match self.inner {
            BgzfWriterImpl::SingleThreaded(mut w) => {
                w.flush()?;
                // Single-threaded writer writes EOF on drop
                Ok(())
            }
            BgzfWriterImpl::MultiThreaded(mut w) => {
                w.finish().map_err(|e| io::Error::other(e.to_string()))?;
                Ok(())
            }
        }
    }
}

/// Type alias for a BAM writer that supports both single and multi-threaded BGZF.
pub type BamWriter = noodles::bam::io::Writer<BgzfWriterEnum>;

/// Create a [`BgzfWriterEnum`] from a writer, selecting single- or multi-threaded based on
/// `threads`.
#[allow(clippy::cast_possible_truncation)]
pub(crate) fn make_bgzf_writer(
    output: Box<dyn Write + Send>,
    threads: usize,
    compression_level: u32,
) -> BgzfWriterEnum {
    if threads > 1 {
        let worker_count = NonZero::new(threads).expect("threads > 1 checked above");
        let mut builder = MultithreadedWriterBuilder::default().set_worker_count(worker_count);

        if let Ok(cl) = CompressionLevel::new(compression_level as u8) {
            builder = builder.set_compression_level(cl);
        }

        BgzfWriterEnum::multi_threaded(builder.build_from_writer(output))
    } else {
        let level = noodles_bgzf::io::writer::CompressionLevel::new(compression_level as u8)
            .unwrap_or_default();
        let writer = noodles_bgzf::io::writer::Builder::default()
            .set_compression_level(level)
            .build_from_writer(output);
        BgzfWriterEnum::single_threaded(writer)
    }
}

/// Write a BAM header (magic, SAM text, reference sequences) to any writer.
///
/// Shared implementation used by both [`RawBamWriter`] and [`IndexingBamWriter`].
///
/// Note: we use a manual implementation rather than delegating to
/// `noodles::bam::io::Writer::write_header` because the noodles encoder
/// produces subtly different output that causes read-back failures with
/// some writer backends (e.g. `MultithreadedWriter`).
///
/// # Errors
///
/// Returns an [`io::Error`] if writing to `writer` fails or if noodles
/// cannot serialize the SAM header text.
#[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
pub fn write_bam_header(writer: &mut impl Write, header: &Header) -> io::Result<()> {
    // BAM magic
    writer.write_all(fgumi_raw_bam::BAM_MAGIC)?;

    // Header text (SAM header serialized using noodles)
    let mut sam_writer = noodles::sam::io::Writer::new(Vec::new());
    sam_writer.write_header(header)?;
    let header_bytes = sam_writer.into_inner();
    let l_text = header_bytes.len() as i32;
    writer.write_all(&l_text.to_le_bytes())?;
    writer.write_all(&header_bytes)?;

    // Reference sequences
    let n_ref = header.reference_sequences().len() as i32;
    writer.write_all(&n_ref.to_le_bytes())?;

    for (name, map) in header.reference_sequences() {
        // l_name: length of name + null terminator
        let l_name = (name.len() + 1) as u32;
        writer.write_all(&l_name.to_le_bytes())?;
        writer.write_all(name)?;
        writer.write_all(&[0u8])?; // null terminator

        // l_ref: reference length
        let l_ref = map.length().get() as i32;
        writer.write_all(&l_ref.to_le_bytes())?;
    }

    Ok(())
}

/// Raw BAM writer for writing raw record bytes directly.
///
/// This bypasses the noodles Record encoding path for maximum performance
/// when you already have raw BAM bytes. Writes records as:
/// - 4-byte `block_size` (little-endian)
/// - raw BAM record bytes
pub struct RawBamWriter {
    inner: BgzfWriterEnum,
}

impl RawBamWriter {
    /// Create a new raw BAM writer from a BGZF writer.
    #[must_use]
    pub fn new(inner: BgzfWriterEnum) -> Self {
        Self { inner }
    }

    /// Write the BAM header.
    ///
    /// # Errors
    /// Returns an error if writing to the underlying writer fails.
    pub fn write_header(&mut self, header: &Header) -> io::Result<()> {
        write_bam_header(&mut self.inner, header)
    }

    /// Write a raw BAM record.
    ///
    /// The bytes should be the raw BAM record data (without the 4-byte `block_size` prefix).
    ///
    /// # Errors
    /// Returns an error if writing to the underlying writer fails.
    #[inline]
    #[allow(clippy::cast_possible_truncation)]
    pub fn write_raw_record(&mut self, record_bytes: &[u8]) -> io::Result<()> {
        use std::io::Write;

        // Write block_size (4 bytes, little-endian)
        let block_size = record_bytes.len() as u32;
        self.inner.write_all(&block_size.to_le_bytes())?;

        // Write raw record bytes
        self.inner.write_all(record_bytes)
    }

    /// Write pre-formatted BAM record bytes directly to the BGZF stream.
    ///
    /// Unlike `write_raw_record`, this does NOT add a length prefix — the data
    /// must already contain properly formatted BAM records (4-byte LE size + record bytes).
    /// Used by the pipeline for bulk secondary output.
    ///
    /// # Errors
    /// Returns an error if writing to the underlying writer fails.
    pub fn write_raw_bytes(&mut self, data: &[u8]) -> io::Result<()> {
        use std::io::Write;
        self.inner.write_all(data)
    }

    /// Finish writing and close the writer.
    ///
    /// # Errors
    /// Returns an error if finalizing the writer fails.
    pub fn finish(self) -> io::Result<()> {
        self.inner.finish()
    }
}

/// Create a raw BAM writer for writing raw record bytes directly.
///
/// This is more efficient than `create_bam_writer` when you already have raw BAM bytes
/// because it bypasses the noodles Record encoding path.
///
/// # Errors
/// Returns an error if the file cannot be created or the header cannot be written.
///
/// # Panics
/// Panics if `threads > 1` but `NonZero::new` fails (should not happen).
pub fn create_raw_bam_writer<P: AsRef<Path>>(
    path: P,
    header: &Header,
    threads: usize,
    compression_level: u32,
) -> Result<RawBamWriter> {
    let path_ref = path.as_ref();
    let output = open_output_writer(path_ref)?;

    let bgzf_writer = make_bgzf_writer(output, threads, compression_level);

    let mut writer = RawBamWriter::new(bgzf_writer);
    writer
        .write_header(header)
        .with_context(|| format!("Failed to write header to: {}", path_ref.display()))?;
    Ok(writer)
}

/// Cached index entry awaiting block position resolution.
///
/// When using multi-threaded compression, we don't know the exact compressed
/// position until the block is actually written. This struct caches the
/// information needed to build the index entry once the position is known.
#[derive(Debug)]
struct CachedIndexEntry {
    /// Block number this record starts in.
    block_number: u64,
    /// Uncompressed offset within the block.
    offset_in_block: usize,
    /// Length of record (including 4-byte size prefix).
    record_len: usize,
    /// Alignment context: (reference ID, start, end, mapped flag).
    alignment_context: Option<(usize, Position, Position, bool)>,
}

/// Extract alignment context from raw BAM bytes.
///
/// Returns `Some((ref_id, start, end, is_mapped))` for mapped reads, or `None`
/// for unmapped reads (the indexer still counts them via a `None` context).
#[inline]
#[allow(clippy::cast_sign_loss, clippy::cast_possible_wrap)]
pub(crate) fn extract_alignment_context(bam: &[u8]) -> Option<(usize, Position, Position, bool)> {
    let v = fgumi_raw_bam::RawRecordView::new(bam);
    let tid = v.ref_id();
    let pos = v.pos();
    let flags = v.flags();

    let is_unmapped = (flags & fgumi_raw_bam::flags::UNMAPPED) != 0;

    // Unmapped reads: return None (indexer handles them specially).
    if tid < 0 || is_unmapped {
        return None;
    }

    // Calculate alignment end from CIGAR using byte-safe access
    // (doesn't require 4-byte alignment like get_cigar_ops).
    let ref_len = fgumi_raw_bam::reference_length_from_raw_bam(bam);

    // BAM positions are 0-based, Position is 1-based.
    let start = Position::try_from((pos + 1) as usize).ok()?;
    let end = Position::try_from((pos + ref_len) as usize).ok()?;

    Some((tid as usize, start, end, true))
}

/// Incrementally builds a BAI index from per-record positions plus per-block
/// compressed offsets, independent of the concrete BGZF writer.
///
/// A writer feeds each written record's `(block_number, offset_in_block,
/// record_len)` via [`BaiBuilder::record`], and each finalized block's
/// `(block_number, compressed_start)` via [`BaiBuilder::note_block`].
/// [`BaiBuilder::resolve`] converts every record whose block position is now
/// known into an indexed chunk (records spanning multiple blocks are deferred
/// until their terminal block is known). [`BaiBuilder::build`] produces the
/// final [`bai::Index`] once all records are resolved.
///
/// This is the shared core behind both [`IndexingBamWriter`] (over the vendored
/// multi-threaded writer) and the sort engine's pool-integrated indexing writer,
/// so the virtual-offset accounting lives in exactly one place.
pub struct BaiBuilder {
    indexer: Indexer<LinearIndex>,
    /// Records awaiting block-position resolution, in write order (which is
    /// non-decreasing block order) — resolved strictly from the front.
    entry_cache: VecDeque<CachedIndexEntry>,
    /// `block_number` -> compressed start offset in the file.
    block_positions: BlockPositions,
    /// Watermark: every block number below this has been pruned from
    /// `block_positions` once its records resolved. Keeps the map bounded by the
    /// in-flight block window instead of growing one entry per output block for
    /// the whole (WGS-scale) file.
    pruned_below: u64,
}

impl BaiBuilder {
    /// Create an empty builder.
    #[must_use]
    pub fn new() -> Self {
        Self {
            indexer: Indexer::default(),
            entry_cache: VecDeque::new(),
            block_positions: BlockPositions::default(),
            pruned_below: 0,
        }
    }

    /// Record a written BAM record's position for later index resolution.
    ///
    /// `block_number`/`offset_in_block` locate the record's first byte in the
    /// uncompressed stream; `record_len` includes the 4-byte length prefix.
    /// The alignment context is extracted from `record_bytes` eagerly.
    pub fn record(
        &mut self,
        block_number: u64,
        offset_in_block: usize,
        record_len: usize,
        record_bytes: &[u8],
    ) {
        let alignment_context = extract_alignment_context(record_bytes);
        self.entry_cache.push_back(CachedIndexEntry {
            block_number,
            offset_in_block,
            record_len,
            alignment_context,
        });
    }

    /// Note the compressed start offset of a finalized block.
    pub fn note_block(&mut self, block_number: u64, compressed_start: u64) {
        self.block_positions.insert(block_number, compressed_start);
    }

    /// Resolve pending records whose block positions are now known.
    ///
    /// Records are cached in write order (non-decreasing block number) and
    /// blocks are noted in increasing order, so resolution proceeds strictly
    /// from the front — which is also the order the indexer must receive records
    /// in. It stops at the first record whose start (or, for a record spanning
    /// multiple blocks, terminal) block has not yet been noted, leaving it and
    /// everything after it cached. This is `O(records resolved)` per call, not
    /// `O(cache length)`, so per-record resolution stays amortized `O(1)`.
    ///
    /// # Errors
    /// Returns an error if the indexer rejects a record.
    #[allow(clippy::cast_possible_truncation)]
    pub fn resolve(&mut self) -> io::Result<()> {
        while let Some(entry) = self.entry_cache.front() {
            // Copy out the fields so the front() borrow ends before pop_front().
            let block_number = entry.block_number;
            let offset_in_block = entry.offset_in_block;
            let record_len = entry.record_len;
            let alignment_context = entry.alignment_context;

            let Some(&block_start) = self.block_positions.get(&block_number) else {
                break;
            };

            // End position may fall in this block or a later one (long reads /
            // large aux payloads can span more than one BGZF block). If the
            // terminal block isn't noted yet, wait for it.
            let end_offset = offset_in_block + record_len;
            let Some(end_vpos) =
                compute_end_vpos(&self.block_positions, block_number, block_start, end_offset)
            else {
                break;
            };

            let start_vpos = VirtualPosition::try_from((block_start, offset_in_block as u16))
                .unwrap_or(VirtualPosition::MIN);
            let chunk = Chunk::new(start_vpos, end_vpos);
            self.indexer.add_record(alignment_context, chunk).map_err(io::Error::other)?;
            self.entry_cache.pop_front();
        }

        // Prune block offsets no remaining or future record can reference. Every
        // pending entry has block_number >= the front's, and future records are
        // recorded with even larger block numbers, so blocks below the front are
        // dead. Advancing a watermark (block numbers are dense sequential
        // serials) makes this O(1) amortized and bounds `block_positions` to the
        // in-flight window.
        if let Some(front_block) = self.entry_cache.front().map(|e| e.block_number) {
            while self.pruned_below < front_block {
                self.block_positions.remove(&self.pruned_below);
                self.pruned_below += 1;
            }
        }

        Ok(())
    }

    /// Number of records still awaiting block-position resolution.
    #[must_use]
    pub(crate) fn pending(&self) -> usize {
        self.entry_cache.len()
    }

    /// Build the final BAI index.
    ///
    /// # Errors
    /// Returns an error if any record remains unresolved (its block position was
    /// never noted), which would indicate a lost or miscounted block.
    pub fn build(self, num_refs: usize) -> io::Result<bai::Index> {
        if !self.entry_cache.is_empty() {
            return Err(io::Error::other(format!(
                "Unflushed index entries remain: {} entries",
                self.entry_cache.len()
            )));
        }
        Ok(self.indexer.build(num_refs))
    }
}

impl Default for BaiBuilder {
    fn default() -> Self {
        Self::new()
    }
}

/// BAM writer that builds an index incrementally during write.
///
/// This writer generates a BAI (BAM index) file alongside the BAM output by tracking
/// virtual file positions and alignment information for each record. The index is
/// built incrementally as records are written, avoiding a second pass through the file.
///
/// Uses multi-threaded BGZF compression with htslib-style block caching for accurate
/// position tracking. Index entries are cached until their corresponding blocks are
/// written, then resolved with the actual compressed positions.
///
/// # Usage
///
/// ```ignore
/// let mut writer = create_indexing_bam_writer(output, &header, 6, 4)?;  // compression=6, threads=4
/// for record in records {
///     writer.write_raw_record(&record)?;
/// }
/// let index = writer.finish()?;  // Returns the BAI index
/// ```
pub struct IndexingBamWriter {
    inner: MultithreadedWriter<File>,
    block_info_rx: BlockInfoRx,
    bai: BaiBuilder,
    num_refs: usize,
    current_block_number: u64,
    current_block_offset: usize,
}

impl IndexingBamWriter {
    /// Create a new indexing BAM writer.
    fn new(inner: MultithreadedWriter<File>, num_refs: usize) -> Self {
        let block_info_rx = inner
            .block_info_receiver()
            .expect("block_info_receiver must be available for IndexingBamWriter")
            .clone();
        Self {
            current_block_number: inner.current_block_number(),
            current_block_offset: inner.buffer_offset(),
            inner,
            block_info_rx,
            bai: BaiBuilder::new(),
            num_refs,
        }
    }

    /// Write the BAM header.
    fn write_header(&mut self, header: &Header) -> io::Result<()> {
        write_bam_header(&mut self.inner, header)?;

        // Flush header to ensure it's in its own block(s)
        self.inner.flush()?;

        // Update position tracking after header
        self.current_block_number = self.inner.current_block_number();
        self.current_block_offset = self.inner.buffer_offset();

        Ok(())
    }

    /// Write a raw BAM record and update the index.
    ///
    /// Extracts alignment information from raw bytes for indexing:
    /// - Reference ID (tid)
    /// - Alignment start position
    /// - Alignment end position (calculated from CIGAR)
    /// - Mapped/unmapped status
    ///
    /// # Errors
    /// Returns an error if writing the record or flushing blocks fails.
    #[inline]
    #[allow(clippy::cast_possible_truncation)]
    pub fn write_raw_record(&mut self, record_bytes: &[u8]) -> io::Result<()> {
        // Capture position BEFORE write
        let block_number = self.current_block_number;
        let offset_in_block = self.current_block_offset;

        // Write record (4-byte size prefix + record bytes)
        let block_size = record_bytes.len() as u32;
        self.inner.write_all(&block_size.to_le_bytes())?;
        self.inner.write_all(record_bytes)?;

        let record_len = 4 + record_bytes.len();
        self.current_block_offset += record_len;

        // Check if a new block was started (buffer was flushed)
        let new_block_number = self.inner.current_block_number();
        if new_block_number > self.current_block_number {
            // Block boundary crossed - update tracking
            self.current_block_number = new_block_number;
            self.current_block_offset = self.inner.buffer_offset();
        }

        // Cache the record's position, then resolve any completed blocks.
        self.bai.record(block_number, offset_in_block, record_len, record_bytes);
        self.flush_completed_blocks()?;

        Ok(())
    }

    /// Drain block-position notifications and resolve any now-complete records.
    fn flush_completed_blocks(&mut self) -> io::Result<()> {
        while let Ok(info) = self.block_info_rx.try_recv() {
            self.bai.note_block(info.block_number, info.compressed_start);
        }
        self.bai.resolve()
    }

    /// Finish writing, flush the BGZF stream, and return the index.
    ///
    /// # Errors
    /// Returns an error if flushing, finalizing, or building the index fails.
    pub fn finish(mut self) -> io::Result<bai::Index> {
        // Flush any remaining data
        self.inner.flush()?;

        // Process all remaining blocks
        // Wait a bit for the writer thread to finish processing
        for _ in 0..100 {
            self.flush_completed_blocks()?;
            if self.bai.pending() == 0 {
                break;
            }
            std::thread::sleep(std::time::Duration::from_millis(10));
        }

        // Finish the writer (writes EOF marker)
        let _ = self.inner.finish().map_err(|e| io::Error::other(e.to_string()))?;

        // Final flush of any remaining entries
        self.flush_completed_blocks()?;

        // Build the final index (errors if any entry remained unresolved).
        let num_refs = self.num_refs;
        self.bai.build(num_refs)
    }
}

/// Create an indexing BAM writer that builds a BAI index during write.
///
/// This writer generates the BAM index incrementally as records are written,
/// avoiding a second pass through the file. The index can be written after
/// calling `finish()`.
///
/// Uses multi-threaded BGZF compression with htslib-style block caching for
/// accurate virtual position tracking. Index entries are cached until their
/// corresponding blocks are written, then resolved with actual positions.
///
/// # Arguments
/// * `path` - Path for the output BAM file
/// * `header` - SAM header (must have coordinate sort order)
/// * `compression_level` - BGZF compression level (0-12; 0 = uncompressed)
/// * `threads` - Number of compression threads (minimum 1)
///
/// # Returns
/// An `IndexingBamWriter` that will generate a BAI index on `finish()`.
///
/// # Errors
/// Returns an error if the file cannot be created or the header cannot be written.
///
/// # Panics
/// Panics if `NonZero::new` fails on the thread count (should not happen).
pub fn create_indexing_bam_writer<P: AsRef<Path>>(
    path: P,
    header: &Header,
    compression_level: u32,
    threads: usize,
) -> Result<IndexingBamWriter> {
    let path_ref = path.as_ref();
    if is_stdout_path(path_ref) {
        anyhow::bail!(
            "Cannot create an indexing BAM writer for stdout (indexing requires a seekable file)"
        );
    }
    let output_file = File::create(path_ref)
        .with_context(|| format!("Failed to create output BAM: {}", path_ref.display()))?;

    // Use multi-threaded writer with block position tracking
    let worker_count = NonZero::new(threads.max(1)).expect("threads.max(1) >= 1");
    let mut builder = MultithreadedWriterBuilder::default().set_worker_count(worker_count);

    #[allow(clippy::cast_possible_truncation)]
    if let Ok(cl) = CompressionLevel::new(compression_level as u8) {
        builder = builder.set_compression_level(cl);
    }

    let bgzf_writer = builder.build_from_writer(output_file);

    let num_refs = header.reference_sequences().len();
    let mut writer = IndexingBamWriter::new(bgzf_writer, num_refs);
    writer
        .write_header(header)
        .with_context(|| format!("Failed to write header to: {}", path_ref.display()))?;

    Ok(writer)
}

/// Append `.bai` to a BAM path to form its conventional sidecar index path.
///
/// This appends to the **full** path, which is the samtools convention, rather
/// than replacing the final extension:
///
/// - `foo.bam` → `foo.bam.bai`
/// - `foo` (no extension) → `foo.bai`
/// - `foo.sorted` → `foo.sorted.bai`
///
/// `Path::with_extension("bam.bai")` is not equivalent and must not be used
/// here: it replaces everything after the last `.`, so it is correct only for
/// paths ending in exactly `.bam`. For anything else it rewrites the basename
/// itself — `foo.sorted` and `foo.other` both collapse to `foo.bam.bai`, so the
/// index lands beside a file that does not exist and two sorts in one directory
/// silently overwrite each other's index.
#[must_use]
pub fn bai_sidecar_path<P: AsRef<Path>>(bam_path: P) -> PathBuf {
    let mut index_os = bam_path.as_ref().as_os_str().to_owned();
    index_os.push(".bai");
    PathBuf::from(index_os)
}

/// Build a BAI index for a finished coordinate-sorted BAM and write it to the
/// conventional `<bam_path>.bai` sidecar, returning the path written.
///
/// The sidecar path comes from [`bai_sidecar_path`], so a BAM whose path does
/// not end in `.bam` still gets a correctly named index.
///
/// **Invariant:** the BAM at `bam_path` must already be closed by its writer.
/// [`noodles::bam::fs::index`] re-reads the file from disk to compute virtual
/// offsets, so any buffered writes must have been flushed first.
///
/// Callers that already hold an in-memory index (for example an incrementally
/// indexing writer) should not use this — it would re-read the whole BAM. Pair
/// [`bai_sidecar_path`] with [`write_bai_index`] instead.
///
/// # Errors
/// Returns an error if indexing the BAM or writing the BAI sidecar fails.
pub fn write_bai_sidecar<P: AsRef<Path>>(bam_path: P) -> Result<PathBuf> {
    let bam_path = bam_path.as_ref();

    let index = noodles::bam::fs::index(bam_path)
        .with_context(|| format!("Failed to index {}", bam_path.display()))?;

    let index_path = bai_sidecar_path(bam_path);
    write_bai_index(&index_path, &index)
        .with_context(|| format!("Failed to write BAI to {}", index_path.display()))?;

    Ok(index_path)
}

/// Write a BAI index to a file.
///
/// # Arguments
/// * `path` - Path for the output BAI file; derive it with [`bai_sidecar_path`]
///   rather than by replacing the BAM path's extension
/// * `index` - The BAI index to write
///
/// # Errors
/// Returns an error if the file cannot be created or writing the index fails.
pub fn write_bai_index<P: AsRef<Path>>(path: P, index: &bai::Index) -> Result<()> {
    let path_ref = path.as_ref();
    // `bai::io::Writer` serializes the index as a great many tiny fields (each
    // bin and chunk writes 8-byte virtual offsets), so writing straight to the
    // file issues one `write()` syscall per field. For a WGS-scale index (tens
    // of MB) that unbuffered tail dominates the post-merge main-thread time.
    // Serialize into memory first, then write the whole sidecar in one call.
    let mut buf = Vec::new();
    bai::io::Writer::new(&mut buf)
        .write_index(index)
        .with_context(|| format!("Failed to serialize BAI index for: {}", path_ref.display()))?;
    std::fs::write(path_ref, &buf)
        .with_context(|| format!("Failed to write index to: {}", path_ref.display()))?;
    Ok(())
}

/// Create a BAM writer and write the header in one operation.
///
/// # Arguments
/// * `path` - Path for the output BAM file (use `-` or `/dev/stdout` for stdout)
/// * `header` - SAM header to write
/// * `threads` - Number of threads for BGZF compression (1 = single-threaded)
/// * `compression_level` - Compression level (0-12; 0 = uncompressed)
///
/// # Returns
/// A BAM writer with the header already written.
///
/// # Errors
/// Returns an error if the file cannot be created or the header cannot be written.
///
/// # Example
/// ```no_run
/// use fgumi_bam_io::writer::create_bam_writer;
/// use noodles::sam::Header;
/// use std::path::Path;
///
/// let header = Header::default();
/// // Single-threaded with fast compression
/// let mut writer = create_bam_writer(Path::new("output.bam"), &header, 1, 1).unwrap();
/// // Multi-threaded with 4 compression threads
/// let mut writer = create_bam_writer(Path::new("output.bam"), &header, 4, 1).unwrap();
/// ```
pub fn create_bam_writer<P: AsRef<Path>>(
    path: P,
    header: &Header,
    threads: usize,
    compression_level: u32,
) -> Result<BamWriter> {
    let path_ref = path.as_ref();
    let output = open_output_writer(path_ref)?;

    let mut bgzf_writer = make_bgzf_writer(output, threads, compression_level);

    // Write the header with the libdeflate-backed helper rather than the
    // noodles encoder: noodles' `Writer::write_header` corrupts the leading
    // BAM bytes when the underlying writer is `MultithreadedWriter`.
    write_bam_header(&mut bgzf_writer, header)
        .with_context(|| format!("Failed to write header to: {}", path_ref.display()))?;
    Ok(noodles::bam::io::Writer::from(bgzf_writer))
}

/// Create an optional BAM writer for rejects or filtered reads.
///
/// This is a convenience function for commands that optionally output rejected/filtered reads.
/// If the path is `None`, returns `None`. If the path is `Some`, creates a writer and writes
/// the header.
///
/// # Arguments
/// * `path` - Optional path for the output BAM file
/// * `header` - SAM header to write
/// * `threads` - Number of threads for BGZF compression (1 = single-threaded)
/// * `compression_level` - Compression level (0-12; 0 = uncompressed)
///
/// # Returns
/// An optional BAM writer, or None if path is None.
///
/// # Errors
/// Returns an error if the file cannot be created or the header cannot be written.
///
/// # Example
/// ```no_run
/// use fgumi_bam_io::writer::create_optional_bam_writer;
/// use noodles::sam::Header;
/// use std::path::PathBuf;
///
/// let header = Header::default();
/// let rejects = Some(PathBuf::from("rejects.bam"));
/// let mut writer = create_optional_bam_writer(rejects, &header, 1, 1).unwrap();
/// ```
pub fn create_optional_bam_writer<P: AsRef<Path>>(
    path: Option<P>,
    header: &Header,
    threads: usize,
    compression_level: u32,
) -> Result<Option<BamWriter>> {
    match path {
        Some(p) => Ok(Some(create_bam_writer(p, header, threads, compression_level)?)),
        None => Ok(None),
    }
}

/// Open an output writer for a path, supporting stdout via `-` or `/dev/stdout`.
pub(crate) fn open_output_writer<P: AsRef<Path>>(path: P) -> Result<Box<dyn Write + Send>> {
    let path_ref = path.as_ref();
    if is_stdout_path(path_ref) {
        Ok(Box::new(std::io::stdout()))
    } else {
        let file = File::create(path_ref)
            .with_context(|| format!("Failed to create output BAM: {}", path_ref.display()))?;
        Ok(Box::new(file))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::sam::header::record::value::{Map, map::ReferenceSequence};
    use std::num::NonZeroUsize;
    use tempfile::NamedTempFile;

    fn create_test_header() -> Header {
        let mut builder = Header::builder();
        let ref_seq = Map::<ReferenceSequence>::new(
            NonZeroUsize::new(100).expect("100 is non-zero constant"),
        );
        builder = builder.add_reference_sequence(b"chr1", ref_seq);
        builder.build()
    }

    #[test]
    fn test_create_bam_writer() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let header = create_test_header();

        let writer = create_bam_writer(temp_file.path(), &header, 1, 6);
        assert!(writer.is_ok());

        Ok(())
    }

    #[test]
    fn test_create_bam_writer_multithreaded() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let header = create_test_header();

        let writer = create_bam_writer(temp_file.path(), &header, 4, 6);
        assert!(writer.is_ok());

        Ok(())
    }

    #[test]
    fn test_create_bam_writer_invalid_path() {
        let header = create_test_header();
        let result = create_bam_writer("/invalid/path/output.bam", &header, 1, 6);
        assert!(result.is_err());
        if let Err(e) = result {
            let err_msg = e.to_string();
            assert!(err_msg.contains("Failed to create output BAM"));
        }
    }

    #[test]
    fn test_create_optional_bam_writer_some() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let header = create_test_header();

        let writer = create_optional_bam_writer(Some(temp_file.path()), &header, 1, 6)?;
        assert!(writer.is_some());

        Ok(())
    }

    #[test]
    fn test_create_optional_bam_writer_none() -> Result<()> {
        let header = create_test_header();

        let writer = create_optional_bam_writer(None::<&str>, &header, 1, 6)?;
        assert!(writer.is_none());

        Ok(())
    }

    #[test]
    fn test_create_optional_bam_writer_invalid_path() {
        let header = create_test_header();
        let result = create_optional_bam_writer(Some("/invalid/path/output.bam"), &header, 1, 6);
        assert!(result.is_err());
    }

    #[test]
    fn test_bgzf_writer_flush_single_threaded() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let output_file: Box<dyn Write + Send> = Box::new(File::create(temp_file.path())?);
        let mut writer = BgzfWriterEnum::single_threaded(BgzfWriter::new(output_file));

        // Write some data and flush
        writer.write_all(b"test data")?;
        let result = writer.flush();
        assert!(result.is_ok());

        Ok(())
    }

    #[test]
    fn test_bgzf_writer_flush_multithreaded() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let output_file: Box<dyn Write + Send> = Box::new(File::create(temp_file.path())?);
        let worker_count = NonZero::new(2).expect("2 is non-zero");
        let compression_level = CompressionLevel::new(6).expect("valid compression level");
        let mut writer = BgzfWriterEnum::multi_threaded(MultithreadedWriter::with_worker_count(
            worker_count,
            output_file,
            compression_level,
        ));

        // Write some data and flush
        writer.write_all(b"test data")?;
        let result = writer.flush();
        assert!(result.is_ok());

        Ok(())
    }

    #[test]
    fn test_bgzf_writer_finish_single_threaded() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let output_file: Box<dyn Write + Send> = Box::new(File::create(temp_file.path())?);
        let mut writer = BgzfWriterEnum::single_threaded(BgzfWriter::new(output_file));

        // Write some data
        writer.write_all(b"test data")?;

        // Finish the writer - this should flush and close properly
        let result = writer.finish();
        assert!(result.is_ok());

        Ok(())
    }

    #[test]
    fn test_bgzf_writer_finish_multithreaded() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let output_file: Box<dyn Write + Send> = Box::new(File::create(temp_file.path())?);
        let worker_count = NonZero::new(2).expect("2 is non-zero");
        let compression_level = CompressionLevel::new(6).expect("valid compression level");
        let mut writer = BgzfWriterEnum::multi_threaded(MultithreadedWriter::with_worker_count(
            worker_count,
            output_file,
            compression_level,
        ));

        // Write some data
        writer.write_all(b"test data")?;

        // Finish the writer - this ensures proper finalization
        let result = writer.finish();
        assert!(result.is_ok());

        Ok(())
    }

    #[test]
    fn test_bgzf_writer_write_single_threaded() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let output_file: Box<dyn Write + Send> = Box::new(File::create(temp_file.path())?);
        let mut writer = BgzfWriterEnum::single_threaded(BgzfWriter::new(output_file));

        // Test writing via the Write trait
        let bytes_written = writer.write(b"test")?;
        assert_eq!(bytes_written, 4);

        Ok(())
    }

    #[test]
    fn test_bgzf_writer_write_multithreaded() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let output_file: Box<dyn Write + Send> = Box::new(File::create(temp_file.path())?);
        let worker_count = NonZero::new(2).expect("2 is non-zero");
        let compression_level = CompressionLevel::new(6).expect("valid compression level");
        let mut writer = BgzfWriterEnum::multi_threaded(MultithreadedWriter::with_worker_count(
            worker_count,
            output_file,
            compression_level,
        ));

        // Test writing via the Write trait
        let bytes_written = writer.write(b"test")?;
        assert_eq!(bytes_written, 4);

        Ok(())
    }

    #[test]
    fn test_is_stdout_path() {
        use crate::paths::is_stdout_path;
        use std::path::Path;
        assert!(is_stdout_path("-"));
        assert!(is_stdout_path("/dev/stdout"));
        assert!(is_stdout_path(Path::new("-")));
        assert!(is_stdout_path(Path::new("/dev/stdout")));

        assert!(!is_stdout_path("output.bam"));
        assert!(!is_stdout_path("/path/to/file.bam"));
        assert!(!is_stdout_path(""));
        assert!(!is_stdout_path("/dev/null"));
    }

    #[test]
    fn test_indexing_writer_rejects_stdout() {
        let header = Header::default();
        let result = create_indexing_bam_writer("-", &header, 6, 2);
        assert!(result.is_err());
        let err = result.err().expect("result should be Err");
        assert!(err.to_string().contains("stdout"));
    }

    /// Create a minimal BAM record for testing.
    ///
    /// Creates a mapped record at the given reference and position with a simple CIGAR.
    #[allow(clippy::cast_possible_truncation)]
    fn create_test_bam_record(ref_id: i32, pos: i32, read_name: &[u8]) -> Vec<u8> {
        let name_with_null = read_name.len() + 1;
        let padding = (4 - (name_with_null % 4)) % 4;
        let l_read_name = (name_with_null + padding) as u8;

        let mapq: u8 = 60;
        let bin: u16 = 4681;
        let n_cigar_op: u16 = 1;
        let flag: u16 = 0;
        let l_seq: u32 = 10;
        let next_ref_id: i32 = -1;
        let next_pos: i32 = -1;
        let tlen: i32 = 0;

        let mut record = Vec::new();
        record.extend_from_slice(&ref_id.to_le_bytes());
        record.extend_from_slice(&pos.to_le_bytes());
        record.push(l_read_name);
        record.push(mapq);
        record.extend_from_slice(&bin.to_le_bytes());
        record.extend_from_slice(&n_cigar_op.to_le_bytes());
        record.extend_from_slice(&flag.to_le_bytes());
        record.extend_from_slice(&l_seq.to_le_bytes());
        record.extend_from_slice(&next_ref_id.to_le_bytes());
        record.extend_from_slice(&next_pos.to_le_bytes());
        record.extend_from_slice(&tlen.to_le_bytes());

        record.extend_from_slice(read_name);
        record.push(0);
        record.extend(std::iter::repeat_n(0u8, padding));

        let cigar_op: u32 = 10 << 4;
        record.extend_from_slice(&cigar_op.to_le_bytes());
        record.extend_from_slice(&[0x11, 0x11, 0x11, 0x11, 0x11]);
        record.extend_from_slice(&[30u8; 10]);

        record
    }

    #[test]
    fn test_create_indexing_bam_writer() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let header = create_test_header();

        let writer = create_indexing_bam_writer(temp_file.path(), &header, 6, 2);
        assert!(writer.is_ok());

        let writer = writer.expect("creating writer should succeed");
        let index = writer.finish();
        assert!(index.is_ok());

        Ok(())
    }

    #[test]
    fn test_indexing_bam_writer_with_records() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let header = create_test_header();

        let mut writer = create_indexing_bam_writer(temp_file.path(), &header, 6, 2)?;

        for i in 0..100 {
            let record = create_test_bam_record(0, i * 100, format!("read{i}").as_bytes());
            writer.write_raw_record(&record)?;
        }

        let index = writer.finish()?;
        assert!(!index.reference_sequences().is_empty());

        Ok(())
    }

    #[test]
    fn test_indexing_bam_writer_multithreaded() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let header = create_test_header();

        let mut writer = create_indexing_bam_writer(temp_file.path(), &header, 6, 4)?;

        for i in 0..1000 {
            let record = create_test_bam_record(0, i * 10, format!("read{i:04}").as_bytes());
            writer.write_raw_record(&record)?;
        }

        let index = writer.finish()?;
        assert!(!index.reference_sequences().is_empty());

        Ok(())
    }

    #[test]
    fn test_indexing_bam_writer_unmapped_records() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let header = create_test_header();

        let mut writer = create_indexing_bam_writer(temp_file.path(), &header, 6, 2)?;

        let record = create_test_bam_record(-1, -1, b"unmapped_read");
        writer.write_raw_record(&record)?;

        let record = create_test_bam_record(0, 100, b"mapped_read");
        writer.write_raw_record(&record)?;

        let index = writer.finish()?;
        assert!(!index.reference_sequences().is_empty());

        Ok(())
    }

    // ========================================================================
    // compute_end_vpos: BAI end-position resolution across BGZF blocks
    // ========================================================================

    /// Build a `block_positions` map mapping `[0..n)` to deterministic
    /// compressed offsets so virtual positions are easy to assert against.
    fn make_block_positions(n: u64) -> BlockPositions {
        // Use 0x10000-step starts so each block is at a distinct, recognizable
        // virtual offset.
        (0..n).map(|i| (i, i * 0x1_0000)).collect()
    }

    fn vpos(compressed: u64, uncompressed: u16) -> VirtualPosition {
        VirtualPosition::try_from((compressed, uncompressed))
            .expect("test inputs must form a valid VirtualPosition")
    }

    #[test]
    fn test_compute_end_vpos_in_block() {
        let positions = make_block_positions(1);
        let v = compute_end_vpos(&positions, 0, 0, 1024).expect("in-block must resolve");
        assert_eq!(v, vpos(0, 1024));
    }

    #[test]
    fn test_compute_end_vpos_at_block_boundary_stays_in_starting_block() {
        // end_offset == MAX_BLOCK_SIZE means the record's last byte fills the
        // start block exactly; the end_vpos sits at the end of that block.
        let positions = make_block_positions(1);
        let v =
            compute_end_vpos(&positions, 0, 0, MAX_BLOCK_SIZE).expect("boundary case must resolve");
        #[allow(clippy::cast_possible_truncation)]
        let expected = vpos(0, MAX_BLOCK_SIZE as u16);
        assert_eq!(v, expected);
    }

    #[test]
    fn test_compute_end_vpos_two_block_record() {
        let positions = make_block_positions(2);
        let v = compute_end_vpos(&positions, 0, 0, MAX_BLOCK_SIZE + 100)
            .expect("two-block case must resolve");
        assert_eq!(v, vpos(0x1_0000, 100));
    }

    #[test]
    fn test_compute_end_vpos_three_block_record() {
        // Long read / large aux payload spans three BGZF blocks.
        let positions = make_block_positions(3);
        let v = compute_end_vpos(&positions, 0, 0, 2 * MAX_BLOCK_SIZE + 50)
            .expect("three-block case must resolve");
        assert_eq!(
            v,
            vpos(0x2_0000, 50),
            "end_vpos should land in block 2 (skipping fully-consumed block 1)"
        );
    }

    #[test]
    fn test_compute_end_vpos_multi_block_lookup_uses_terminal_block() {
        // A record that spans blocks 5..7 (start block 5, ends in block 7).
        let mut positions = BlockPositions::default();
        positions.insert(5u64, 0x5_0000u64);
        // Block 6 missing on purpose — we only need the terminal block start.
        positions.insert(7u64, 0x7_0000u64);

        let v = compute_end_vpos(&positions, 5, 0x5_0000, 2 * MAX_BLOCK_SIZE + 7)
            .expect("terminal block start is sufficient");
        assert_eq!(v, vpos(0x7_0000, 7));
    }

    #[test]
    fn test_compute_end_vpos_returns_none_when_terminal_block_unknown() {
        // Record spans into block 1 but block 1's start hasn't been recorded.
        let positions = make_block_positions(1);
        let v = compute_end_vpos(&positions, 0, 0, MAX_BLOCK_SIZE + 100);
        assert!(v.is_none(), "must defer when the end block isn't ready");
    }

    #[test]
    fn test_compute_end_vpos_offset_within_starting_block() {
        // Start somewhere in the middle of block 3; record extends into block 4.
        let mut positions = BlockPositions::default();
        positions.insert(3u64, 0x3_0000u64);
        positions.insert(4u64, 0x4_0000u64);

        let start_offset = 60_000usize;
        let record_len = 10_000usize; // crosses one boundary (60000 + 10000 = 70000)
        let end_offset = start_offset + record_len;
        let v = compute_end_vpos(&positions, 3, 0x3_0000, end_offset)
            .expect("two-block case with non-zero start offset");
        // overflow = end_offset - MAX_BLOCK_SIZE = 70000 - 65280 = 4720
        #[allow(clippy::cast_possible_truncation)]
        let expected = vpos(0x4_0000, 4720u16);
        assert_eq!(v, expected);
    }

    #[test]
    fn test_write_raw_bytes() -> Result<()> {
        let temp = tempfile::NamedTempFile::new()?;
        let header = noodles::sam::Header::default();
        let mut writer = create_raw_bam_writer(temp.path(), &header, 1, 1)?;

        let record_bytes = vec![1, 2, 3, 4, 5];
        let mut formatted = Vec::new();
        #[allow(clippy::cast_possible_truncation)]
        let len = record_bytes.len() as u32;
        formatted.extend_from_slice(&len.to_le_bytes());
        formatted.extend_from_slice(&record_bytes);

        writer.write_raw_bytes(&formatted)?;
        writer.finish()?;

        assert!(temp.path().metadata()?.len() > 0);
        Ok(())
    }

    /// The samtools convention appends `.bai` to the whole BAM path rather than
    /// replacing the final extension. `Path::with_extension("bam.bai")` -- what
    /// these call sites used to do -- happens to agree only for paths ending in
    /// exactly `.bam`; for anything else it rewrites the basename, so distinct
    /// outputs collide on one index path.
    #[rstest::rstest]
    #[case::plain_bam("foo.bam", "foo.bam.bai")]
    #[case::no_extension("foo", "foo.bai")]
    #[case::non_bam_extension("foo.sorted", "foo.sorted.bai")]
    #[case::dotted_stem("sample.v2.bam", "sample.v2.bam.bai")]
    #[case::with_directory("a/b/out", "a/b/out.bai")]
    #[case::directory_with_dots("a.d/out.bam", "a.d/out.bam.bai")]
    fn test_bai_sidecar_path_appends_to_full_path(#[case] bam: &str, #[case] expected: &str) {
        assert_eq!(bai_sidecar_path(bam), PathBuf::from(expected));
    }

    /// The old expression and the new helper must agree exactly on `.bam` paths
    /// (so this change is a no-op for the common case) and must differ on the
    /// paths that were broken -- including two distinct outputs that previously
    /// collapsed onto the same index path and would silently clobber each other.
    #[test]
    fn test_bai_sidecar_path_differs_from_with_extension_only_off_dot_bam() {
        let with_ext = |p: &str| Path::new(p).with_extension("bam.bai");

        assert_eq!(bai_sidecar_path("foo.bam"), with_ext("foo.bam"), "no change for .bam paths");

        for broken in ["foo", "foo.sorted", "foo.other"] {
            assert_ne!(
                bai_sidecar_path(broken),
                with_ext(broken),
                "{broken} was mis-named by the old expression and must change",
            );
        }

        // The collision the old expression produced: distinct outputs, one index.
        assert_eq!(with_ext("foo.sorted"), with_ext("foo.other"));
        assert_ne!(bai_sidecar_path("foo.sorted"), bai_sidecar_path("foo.other"));
    }
}

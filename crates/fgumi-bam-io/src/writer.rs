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
use noodles_csi::binning_index::index::ReferenceSequence;
use noodles_csi::binning_index::index::reference_sequence::bin::{Bin, Chunk};
use noodles_csi::binning_index::index::reference_sequence::index::LinearIndex;
// The `ReferenceSequence` trait (aliased to avoid colliding with the concrete
// `index::ReferenceSequence` above) provides `metadata()`.
use noodles_csi::binning_index::ReferenceSequence as _;
use noodles_csi::binning_index::{BinningIndex, Indexer};
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
/// Returns `Some((ref_id, start, end, is_mapped))` for any read that has a
/// reference and a position — **including a placed-but-unmapped read** (e.g. the
/// unmapped mate of a mapped read, which carries its mate's `tid`/`pos` so it
/// sorts alongside it). Such reads are position-binned exactly as htslib/samtools
/// do, so region queries and `idxstats` see them; the `is_mapped` flag records
/// that they are unmapped without excluding them from the index. Only a truly
/// unplaced read (no reference or no position) returns `None`, which the indexer
/// counts toward the unplaced total.
#[inline]
#[allow(clippy::cast_sign_loss, clippy::cast_possible_wrap)]
pub(crate) fn extract_alignment_context(bam: &[u8]) -> Option<(usize, Position, Position, bool)> {
    let v = fgumi_raw_bam::RawRecordView::new(bam);
    let tid = v.ref_id();
    let pos = v.pos();
    let flags = v.flags();

    // Only a read with no reference / no position is unplaced. A read carrying a
    // valid tid+pos is binned by position regardless of its mapped flag — this is
    // what htslib does, and what makes placed-but-unmapped mates queryable.
    if tid < 0 || pos < 0 {
        return None;
    }

    let is_mapped = (flags & fgumi_raw_bam::flags::UNMAPPED) == 0;

    // Reference span from the CIGAR (byte-safe access; no 4-byte alignment
    // needed). A placed-but-unmapped read has no CIGAR, so this is 0; htslib
    // bins such a read over `[pos, pos + 1)`, so floor the span at one base
    // (matching `bam_endpos`, which returns `pos + max(rlen, 1)`).
    let ref_len = fgumi_raw_bam::reference_length_from_raw_bam(bam).max(1);

    // BAM positions are 0-based, Position is 1-based.
    let start = Position::try_from((pos + 1) as usize).ok()?;
    let end = Position::try_from((pos + ref_len) as usize).ok()?;

    Some((tid as usize, start, end, is_mapped))
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
    /// The in-progress run of consecutive records mapping to the same bin, whose
    /// chunks are coalesced into one (see [`Self::resolve`]). `None` before the
    /// first placed record and after any unplaced one.
    run: Option<ChunkRun>,
}

/// A run of consecutive same-bin records being coalesced into a single chunk.
///
/// noodles' `Indexer` merges two chunks only when they are exactly adjacent in
/// virtual-position space, which breaks at every BGZF block boundary and leaves
/// one chunk per block per bin. htslib instead emits one chunk per run of
/// same-bin records, spanning blocks. Anchoring every record's chunk at the
/// run's start (rather than the record's own start) makes noodles coalesce the
/// whole run — reproducing htslib's chunk shape. This mirrors the
/// `save_bin`/`save_off` state in htslib's `hts_idx_push` (`hts.c`).
struct ChunkRun {
    reference_sequence_id: usize,
    bin: usize,
    start: VirtualPosition,
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
            run: None,
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

            // Anchor the chunk at the current run's start so noodles coalesces
            // the whole run of same-bin records into one block-spanning chunk
            // (see [`ChunkRun`]). A new bin/reference — or an unplaced record —
            // starts a fresh run. The record's own context still drives bin
            // routing and the linear index inside `add_record`; only the chunk's
            // start position is widened.
            let chunk_start = if let Some((reference_sequence_id, start, end, _)) =
                alignment_context
            {
                let bin = reg2bin(start, end);
                // Read the current run out by value first, so reassigning
                // `self.run` below doesn't overlap the borrow.
                let current = self.run.as_ref().map(|r| (r.reference_sequence_id, r.bin, r.start));
                match current {
                    Some((run_ref, run_bin, run_start))
                        if run_ref == reference_sequence_id && run_bin == bin =>
                    {
                        run_start
                    }
                    _ => {
                        self.run = Some(ChunkRun { reference_sequence_id, bin, start: start_vpos });
                        start_vpos
                    }
                }
            } else {
                // Unplaced record: `add_record` ignores the chunk and counts it
                // as unplaced. End any run so the next placed record starts fresh.
                self.run = None;
                start_vpos
            };
            self.indexer
                .add_record(alignment_context, Chunk::new(chunk_start, end_vpos))
                .map_err(io::Error::other)?;
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
    /// The index is compacted (see `compact_bai_index`) so its on-disk size
    /// matches what htslib/samtools produce; noodles' `Indexer` alone emits a
    /// spec-valid but far larger index (see that function for why).
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
        let index = self.indexer.build(num_refs);
        // `reg2bin` (run-batching) hardcodes BAI's binning parameters, which must
        // match the ones the indexer actually routed records with — otherwise a
        // run's bin could disagree with the bin noodles binned into. They agree
        // today (`Indexer::default()` is 14/5); this catches a future drift.
        debug_assert_eq!(
            (u32::from(index.min_shift()), u32::from(index.depth())),
            (BAI_MIN_SHIFT, BAI_DEPTH),
            "reg2bin's BAI constants must match the indexer's binning parameters",
        );
        Ok(compact_bai_index(&index))
    }
}

impl Default for BaiBuilder {
    fn default() -> Self {
        Self::new()
    }
}

/// Minimum compressed-offset span, in bytes, below which htslib folds a bin
/// into its parent. Matches htslib's `HTS_MIN_MARKER_DIST` (`hts.c`); a 64 KiB
/// compressed span is roughly one BGZF block.
const HTS_MIN_MARKER_DIST: u64 = 0x1_0000;

/// Compact a freshly built BAI index the way htslib's `compress_binning`
/// (`hts.c`) does — a finalization step noodles' `Indexer` does not perform.
///
/// noodles builds a spec-valid index but, for dense short-read data, one that
/// is roughly an order of magnitude larger than samtools': it keeps a distinct
/// bin at every level down to the 16 KiB leaves (redundant with the linear
/// index it also stores) and merges only exactly-adjacent chunks. htslib folds
/// the small fine-grained bins into their parents and merges every chunk that
/// shares a BGZF block, then relies on the linear index for finer lookup. This
/// reproduces that pass on the built index.
///
/// Only the binning index changes: the linear index, metadata, and header are
/// carried over untouched (htslib leaves them alone here too), so region
/// queries return the same records — the index just shrinks.
fn compact_bai_index(index: &bai::Index) -> bai::Index {
    let depth = index.depth();

    let reference_sequences: Vec<ReferenceSequence<LinearIndex>> = index
        .reference_sequences()
        .iter()
        .map(|reference_sequence| {
            // Pull the bins into a plain map we can freely mutate, keyed by bin
            // number. Cloning here also lets us leave `indexmap` out of this
            // crate's direct dependencies (the rebuilt map's type is inferred
            // from `ReferenceSequence::new`).
            let bins: HashMap<usize, Vec<Chunk>> = reference_sequence
                .bins()
                .iter()
                .map(|(&bin_id, bin)| (bin_id, bin.chunks().to_vec()))
                .collect();
            // `compress_binning` already yields `(bin_id, Bin)`; collect into the
            // `IndexMap` `ReferenceSequence::new` wants (type inferred, so
            // `indexmap` stays out of this crate's direct dependencies).
            let bins = compress_binning(bins, depth).into_iter().collect();
            ReferenceSequence::new(
                bins,
                reference_sequence.index().clone(),
                reference_sequence.metadata().cloned(),
            )
        })
        .collect();

    let mut builder = bai::Index::builder()
        .set_min_shift(index.min_shift())
        .set_depth(depth)
        .set_reference_sequences(reference_sequences);
    if let Some(header) = index.header() {
        builder = builder.set_header(header.clone());
    }
    if let Some(count) = index.unplaced_unmapped_record_count() {
        builder = builder.set_unplaced_unmapped_record_count(count);
    }
    builder.build()
}

/// Port of htslib's `compress_binning` for one reference sequence's bins.
///
/// Returns the surviving bins (in ascending bin-number order for a
/// deterministic index) after two passes:
///
/// 1. **Fold small bins into their parent**, walking finest level to coarsest
///    so folding cascades: a bin whose chunks span less than
///    [`HTS_MIN_MARKER_DIST`] *compressed* bytes has its chunks appended to its
///    parent bin and is dropped — but only if the parent bin already exists
///    (htslib never creates one), matching htslib exactly. This removes the
///    per-16 KiB leaf bins that dominate an un-compacted index.
/// 2. **Merge chunks that share a BGZF block** within each surviving bin,
///    collapsing the many single-record chunks a coordinate sort produces.
///
/// `depth` is the index depth / number of non-root levels (5 for BAI).
fn compress_binning(mut bins: HashMap<usize, Vec<Chunk>>, depth: u8) -> Vec<(usize, Bin)> {
    // Pass 1: fold, finest level (`depth`) up to level 1. The root (bin 0) has
    // no parent and is never folded.
    for level in (1..=depth).rev() {
        let level_first = bin_level_first(level);
        // Snapshot the candidate bins: pass 1 removes folded bins as it goes.
        // Sorting keeps folding order (and thus the result) deterministic.
        let mut bin_ids: Vec<usize> =
            bins.keys().copied().filter(|&bin_id| bin_id >= level_first).collect();
        bin_ids.sort_unstable();

        for bin_id in bin_ids {
            let Some(chunks) = bins.get_mut(&bin_id) else { continue };
            if chunks.is_empty() {
                continue;
            }
            // htslib measures the span over start-sorted chunks: first chunk's
            // start to last chunk's end.
            chunks.sort_unstable_by_key(|chunk| u64::from(chunk.start()));
            let span = coffset(chunks.last().unwrap().end())
                .saturating_sub(coffset(chunks.first().unwrap().start()));
            if span >= HTS_MIN_MARKER_DIST {
                continue;
            }

            let parent = bin_parent(bin_id);
            if bins.contains_key(&parent) {
                let mut folded = bins.remove(&bin_id).unwrap();
                bins.get_mut(&parent).unwrap().append(&mut folded);
            }
            // Parent absent: keep the bin, as htslib does.
        }
    }

    // Pass 2: merge same-block chunks within each surviving bin, emitting bins
    // in ascending order.
    let mut bin_ids: Vec<usize> = bins.keys().copied().collect();
    bin_ids.sort_unstable();
    bin_ids
        .into_iter()
        .map(|bin_id| {
            let mut chunks = bins.remove(&bin_id).unwrap();
            chunks.sort_unstable_by_key(|chunk| u64::from(chunk.start()));
            (bin_id, Bin::new(merge_same_block_chunks(chunks)))
        })
        .collect()
}

/// Merge chunks that begin in the same BGZF block as the running chunk ends in
/// (comparing compressed offsets), matching htslib's post-fold chunk merge.
/// Input must be sorted by start virtual position.
fn merge_same_block_chunks(chunks: Vec<Chunk>) -> Vec<Chunk> {
    let mut merged: Vec<Chunk> = Vec::with_capacity(chunks.len());
    for chunk in chunks {
        if let Some(last) = merged.last_mut()
            && coffset(last.end()) >= coffset(chunk.start())
        {
            if u64::from(last.end()) < u64::from(chunk.end()) {
                *last = Chunk::new(last.start(), chunk.end());
            }
            continue;
        }
        merged.push(chunk);
    }
    merged
}

/// Compressed-file offset (top 48 bits) of a BGZF virtual position — htslib's
/// `voffset >> 16`.
fn coffset(virtual_position: VirtualPosition) -> u64 {
    u64::from(virtual_position) >> 16
}

/// First bin number at binning-tree `level` (0 = root). Matches htslib's
/// `hts_bin_first(l) = ((1 << 3l) - 1) / 7`.
fn bin_level_first(level: u8) -> usize {
    ((1usize << (3 * u32::from(level))) - 1) / 7
}

/// Parent of a non-root bin in the 8-ary binning tree. Matches htslib's
/// `hts_bin_parent(b) = (b - 1) >> 3`.
fn bin_parent(bin_id: usize) -> usize {
    (bin_id - 1) >> 3
}

/// BAI binning parameters, fixed by the format. These must match the values
/// noodles' `Indexer::default()` uses (14 / 5), because [`reg2bin`] here decides
/// where a run breaks and must agree with the bin noodles routes each record to.
const BAI_MIN_SHIFT: u32 = 14;
const BAI_DEPTH: u32 = 5;

/// Smallest bin fully containing the 1-based, closed interval `[start, end]`. Used
/// only to detect run boundaries; the on-disk bin is still assigned by noodles.
///
/// Delegates to [`fgumi_raw_bam::reg2bin`], the single implementation of the SAM
/// spec's binning scheme, rather than repeating it here — a second copy of the same
/// scheme is exactly what that one was factored out to prevent. Only the coordinate
/// convention is adapted: it takes a 0-based start and a *half-open* end, so a
/// 1-based inclusive `end` is already the half-open 0-based end.
///
/// Its `u16` return cannot truncate for any interval BAI can represent: the largest
/// BAI bin id is 37,448 (`((1 << 15) - 1) / 7 + (2^29 >> 14) - 1`), well inside
/// `u16::MAX`, and BAI cannot address beyond 2^29 at all.
fn reg2bin(start: Position, end: Position) -> usize {
    let beg = i32::try_from(usize::from(start) - 1).unwrap_or(i32::MAX);
    let end = i32::try_from(usize::from(end)).unwrap_or(i32::MAX);
    usize::from(fgumi_raw_bam::reg2bin(beg, end))
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

    /// Build a raw BAM record for a *placed-but-unmapped* read: a valid
    /// reference/position, the UNMAPPED flag set, and no CIGAR — as produced for
    /// the unmapped mate of a mapped read.
    #[allow(clippy::cast_possible_truncation)]
    fn create_placed_unmapped_bam_record(ref_id: i32, pos: i32, read_name: &[u8]) -> Vec<u8> {
        let name_with_null = read_name.len() + 1;
        let padding = (4 - (name_with_null % 4)) % 4;
        let l_read_name = (name_with_null + padding) as u8;

        let mut record = Vec::new();
        record.extend_from_slice(&ref_id.to_le_bytes());
        record.extend_from_slice(&pos.to_le_bytes());
        record.push(l_read_name);
        record.push(0); // mapq
        record.extend_from_slice(&0u16.to_le_bytes()); // bin
        record.extend_from_slice(&0u16.to_le_bytes()); // n_cigar_op = 0 (no CIGAR)
        record.extend_from_slice(&fgumi_raw_bam::flags::UNMAPPED.to_le_bytes()); // flag
        record.extend_from_slice(&0u32.to_le_bytes()); // l_seq = 0
        record.extend_from_slice(&(-1i32).to_le_bytes()); // next_ref_id
        record.extend_from_slice(&(-1i32).to_le_bytes()); // next_pos
        record.extend_from_slice(&0i32.to_le_bytes()); // tlen
        record.extend_from_slice(read_name);
        record.push(0);
        record.extend(std::iter::repeat_n(0u8, padding));
        record
    }

    #[test]
    fn extract_alignment_context_bins_placed_unmapped_reads() {
        // Placed-but-unmapped (valid ref+pos, UNMAPPED flag, no CIGAR) must be
        // binned over a 1-base span at its position, flagged unmapped — so region
        // queries and idxstats see it, matching htslib.
        let rec = create_placed_unmapped_bam_record(0, 100, b"unmapped_mate");
        let (ref_id, start, end, is_mapped) =
            extract_alignment_context(&rec).expect("placed-unmapped read must be binned");
        assert_eq!(ref_id, 0);
        assert_eq!(usize::from(start), 101, "0-based pos 100 -> 1-based 101");
        assert_eq!(usize::from(end), 101, "no CIGAR -> 1-base span [pos, pos+1)");
        assert!(!is_mapped, "flag must still record it as unmapped");

        // A mapped read is binned over its CIGAR span and flagged mapped.
        let mapped = create_test_bam_record(0, 100, b"mapped");
        let (_, m_start, m_end, m_mapped) =
            extract_alignment_context(&mapped).expect("mapped read must be binned");
        assert_eq!(usize::from(m_start), 101);
        assert_eq!(usize::from(m_end), 110, "10M CIGAR -> span of 10");
        assert!(m_mapped);

        // A truly unplaced read (no reference) is not binned.
        assert!(
            extract_alignment_context(&create_placed_unmapped_bam_record(-1, -1, b"unplaced"))
                .is_none(),
            "unplaced read -> None (counted as unplaced by the indexer)"
        );
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

    // ---- BAI index compaction (port of htslib's compress_binning) ----
    // (`vpos(compressed, uncompressed)` helper is defined above.)

    /// Collect the compacted bins into a plain map keyed by bin number.
    fn compacted(bins: HashMap<usize, Vec<Chunk>>) -> HashMap<usize, Vec<Chunk>> {
        compress_binning(bins, 5)
            .into_iter()
            .map(|(bin_id, bin)| (bin_id, bin.chunks().to_vec()))
            .collect()
    }

    #[test]
    fn compress_binning_folds_small_leaf_bin_into_existing_parent() {
        // Leaf bin 4681 (chunks within a tiny compressed span) and its L4
        // parent bin 585 = (4681 - 1) >> 3.
        let mut bins = HashMap::new();
        bins.insert(
            4681,
            vec![
                Chunk::new(vpos(1000, 0), vpos(1000, 500)),
                Chunk::new(vpos(1200, 0), vpos(1400, 10)),
            ],
        );
        bins.insert(585, vec![Chunk::new(vpos(900, 0), vpos(950, 0))]);

        let out = compacted(bins);

        assert!(!out.contains_key(&4681), "small leaf bin must fold into its parent");
        assert!(out.contains_key(&585), "parent bin (no grandparent present) survives");
        // The parent now covers both its own and the leaf's chunk range.
        let lo = out[&585].iter().map(|c| u64::from(c.start())).min().unwrap();
        let hi = out[&585].iter().map(|c| u64::from(c.end())).max().unwrap();
        assert!(lo <= u64::from(vpos(900, 0)) && hi >= u64::from(vpos(1400, 10)));
    }

    #[test]
    fn compress_binning_keeps_leaf_bin_when_parent_absent() {
        // Small leaf bin, but bin 585 is not present, so htslib (and this port)
        // keep it rather than synthesizing a parent.
        let mut bins = HashMap::new();
        bins.insert(4681, vec![Chunk::new(vpos(1000, 0), vpos(1000, 500))]);

        let out = compacted(bins);

        assert!(out.contains_key(&4681), "leaf with no existing parent bin must survive");
    }

    #[test]
    fn compress_binning_keeps_bin_spanning_at_least_one_block() {
        // Span exceeds HTS_MIN_MARKER_DIST (64 KiB compressed), so it stays put
        // even though its parent exists.
        let mut bins = HashMap::new();
        bins.insert(4681, vec![Chunk::new(vpos(0, 0), vpos(70_000, 0))]);
        bins.insert(585, vec![Chunk::new(vpos(0, 0), vpos(10, 0))]);

        let out = compacted(bins);

        assert!(out.contains_key(&4681), "leaf spanning >= 64 KiB compressed must not fold");
        assert!(out.contains_key(&585));
    }

    #[test]
    fn compress_binning_merges_chunks_in_same_bgzf_block() {
        // Bin 0 (the root, never folded) with three chunks that each begin in
        // the same compressed block the previous one ends in -> one chunk.
        let mut bins = HashMap::new();
        bins.insert(
            0,
            vec![
                Chunk::new(vpos(5000, 0), vpos(5000, 100)),
                Chunk::new(vpos(5000, 100), vpos(5000, 200)),
                Chunk::new(vpos(5000, 200), vpos(6000, 0)),
            ],
        );

        let out = compacted(bins);

        assert_eq!(out[&0].len(), 1, "chunks sharing a block must merge into one");
        assert_eq!(u64::from(out[&0][0].start()), u64::from(vpos(5000, 0)));
        assert_eq!(u64::from(out[&0][0].end()), u64::from(vpos(6000, 0)));
    }

    #[test]
    fn compress_binning_keeps_chunks_in_different_blocks() {
        let mut bins = HashMap::new();
        bins.insert(
            0,
            vec![
                Chunk::new(vpos(1000, 0), vpos(1000, 100)),
                Chunk::new(vpos(9000, 0), vpos(9000, 100)),
            ],
        );

        let out = compacted(bins);

        assert_eq!(out[&0].len(), 2, "chunks in different blocks must not merge");
    }

    #[test]
    fn compact_bai_index_folds_leaf_bin_and_preserves_reference_count() {
        // Ten mapped reads packed into the first 16 KiB window land in leaf bin
        // 4681; one read straddling the 16 KiB boundary lands in parent bin 585
        // (so the parent exists). All chunks sit in one compressed block, so the
        // leaf's span is tiny and must fold.
        let mut indexer = Indexer::<LinearIndex>::default();
        let context = |start: usize, end: usize| {
            Some((0usize, Position::new(start).unwrap(), Position::new(end).unwrap(), true))
        };
        for i in 0..10u16 {
            let start = vpos(100, i * 20);
            let end = vpos(100, i * 20 + 20);
            let pos = usize::from(i) * 30 + 1;
            indexer.add_record(context(pos, pos + 149), Chunk::new(start, end)).unwrap();
        }
        // Boundary-spanning read -> parent bin 585, later in the file.
        indexer
            .add_record(context(16_300, 16_449), Chunk::new(vpos(100, 250), vpos(100, 300)))
            .unwrap();

        let index = indexer.build(1);
        assert!(
            index.reference_sequences()[0].bins().contains_key(&4681),
            "sanity: the uncompacted noodles index keeps the 16 KiB leaf bin",
        );

        let index = compact_bai_index(&index);
        let reference_sequences = index.reference_sequences();
        assert_eq!(reference_sequences.len(), 1, "reference count is preserved");
        let bins = reference_sequences[0].bins();
        assert!(!bins.contains_key(&4681), "compaction folds the leaf bin away");
        assert!(bins.contains_key(&585), "the parent bin is retained");
    }

    #[test]
    fn reg2bin_matches_bai_bin_numbering() {
        let pos = |p: usize| Position::new(p).unwrap();
        // Reads within a single 16 KiB window -> consecutive leaf bins.
        assert_eq!(reg2bin(pos(1), pos(150)), 4681);
        assert_eq!(reg2bin(pos(16_385), pos(16_534)), 4682);
        // A read straddling a 16 KiB boundary but inside one 128 KiB window ->
        // the level-4 parent bin 585 (this is the bin leaf reads fold into).
        assert_eq!(reg2bin(pos(16_300), pos(16_449)), 585);
        // reg2bin(4681)'s parent must be 585, matching the fold target.
        assert_eq!(bin_parent(4681), 585);
    }
}

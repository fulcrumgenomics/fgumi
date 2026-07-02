//! Raw-bytes sorting implementation for BAM files.
//!
//! Items in this module are reachable only via the crate-root re-exports
//! ([`crate::RawExternalSorter`] etc.) and via siblings inside `fgumi-sort`.
//! Some helpers (legacy spill-chunk writer, alternative readers) are retained
//! for benchmarking and as fall-back paths.

#![allow(dead_code)]
// Pre-existing module body follows.
//! Instead of fully decoding each BAM record into `RecordBuf`, it uses noodles'
//! lazy `Record` type that stores raw bytes and only parses fields on-demand.
//!
//! # Performance Benefits
//!
//! - **3-4x lower memory usage**: Raw bytes are ~200-400 bytes vs ~800-1200 bytes decoded
//! - **No re-encoding overhead**: Records are written back as raw bytes
//! - **Lazy field access**: Only sort-key fields are parsed
//!
//! # Algorithm
//!
//! 1. Read BAM records as lazy `Record` objects (raw bytes)
//! 2. Extract only the fields needed for sort keys (tid, pos, flags, name, MI)
//! 3. Sort by keys while keeping raw records
//! 4. Write raw record bytes to output

use crate::inline::{
    CbKey32, InMemoryChunk, ProbeableBuffer, RecordBuffer, TEMPLATE_HEADER_SIZE, TemplateKey,
    TemplateKey24, TemplateKey40, TemplateLaneKey, TemplateRecordBuffer, TemplateRecordRef,
    TertKey32,
};
use crate::keys::{QuerynameComparator, RawSortKey, SortOrder};
use crate::memory_probe::{
    BufferProbeStats, ConsumerProbeStats, MergeProbe, SpillProbe, force_mi_collect, log_snapshot,
};
use crate::merge_slots::SortMergeSlot;
use crate::pooled_chunk_writer::PooledChunkWriter;
use crate::read_ahead::{PooledInputStream, RawReadAheadReader, RecordSource};
use crate::tmp_dir_alloc::TmpDirAllocator;
use crate::worker_pool::SortWorkerPool;
use anyhow::{Context as _, Result};
use crossbeam_channel::{Receiver, Sender, bounded};
use fgumi_bam_io::ProgressTracker;
use fgumi_bam_io::create_raw_bam_reader;
use fgumi_bam_io::{ChainedReader, TeeReader, is_stdin_path};
use fgumi_raw_bam::SamTag;
use log::{debug, info};
use noodles::sam::Header;
use noodles::sam::header::record::value::map::read_group::tag as rg_tag;
use noodles_bgzf::io::{
    MultithreadedWriter, Reader as BgzfReader, Writer as BgzfWriter, multithreaded_writer,
    writer::CompressionLevel,
};
use std::collections::{HashMap, HashSet};
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::num::NonZero;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::thread::{self, JoinHandle};
use std::time::{Duration, Instant};
use tempfile::TempDir;

// ============================================================================
// Per-Phase Timing for Sort Pipeline
// ============================================================================

/// Tracks wall-clock time spent in each phase of the sort pipeline.
///
/// Used to identify bottlenecks and validate thread architecture changes.
/// All times are cumulative (multiple spill cycles accumulate).
#[derive(Debug, Default)]
struct SortPhaseTimer {
    /// Time reading records from input BAM (includes BGZF decompression).
    read_secs: f64,
    /// Time sorting in-memory buffers (rayon parallel sort or single-threaded).
    sort_secs: f64,
    /// Time writing sorted chunks to temp files (BGZF compression).
    spill_write_secs: f64,
    /// Time consolidating temp files when limit exceeded.
    consolidate_secs: f64,
    /// Time in the final k-way merge phase (includes reader decompression + writer compression).
    merge_secs: f64,
    /// Time writing in-memory-only output (no merge needed).
    write_output_secs: f64,
    /// Number of spill cycles (sort + write).
    spill_count: usize,
    /// Number of consolidation merges.
    consolidate_count: usize,
    /// Total bytes written to spill files.
    total_spill_bytes: u64,
    /// Wall-clock start of the entire sort operation.
    overall_start: Option<Instant>,
    /// Tracks the start of the current read span (between spills).
    read_span_start: Option<Instant>,
}

impl SortPhaseTimer {
    fn new() -> Self {
        Self {
            overall_start: Some(Instant::now()),
            read_span_start: Some(Instant::now()),
            ..Default::default()
        }
    }

    /// End a read span (call before sort/spill). Returns elapsed read time.
    fn end_read_span(&mut self) -> Duration {
        if let Some(start) = self.read_span_start.take() {
            let elapsed = start.elapsed();
            self.read_secs += elapsed.as_secs_f64();
            elapsed
        } else {
            Duration::ZERO
        }
    }

    /// Start a new read span (call after spill write completes).
    fn begin_read_span(&mut self) {
        self.read_span_start = Some(Instant::now());
    }

    /// Time a closure and accumulate elapsed seconds into `field`.
    fn time<T>(field: &mut f64, f: impl FnOnce() -> T) -> T {
        let start = Instant::now();
        let result = f();
        *field += start.elapsed().as_secs_f64();
        result
    }

    /// Time a sort operation.
    fn time_sort<T>(&mut self, f: impl FnOnce() -> T) -> T {
        Self::time(&mut self.sort_secs, f)
    }

    /// Time a spill write operation.
    fn time_spill_write<T>(&mut self, f: impl FnOnce() -> Result<T>) -> Result<T> {
        let result = Self::time(&mut self.spill_write_secs, f);
        self.spill_count += 1;
        result
    }

    /// Record the size of a spill file.
    fn record_spill_size(&mut self, path: &Path) {
        if let Ok(meta) = std::fs::metadata(path) {
            self.total_spill_bytes += meta.len();
        }
    }

    /// Time a consolidation operation.
    fn time_consolidate(&mut self, f: impl FnOnce() -> Result<()>) -> Result<()> {
        let start = Instant::now();
        let result = f();
        let elapsed = start.elapsed().as_secs_f64();
        if elapsed > 0.001 {
            // Only count if consolidation actually happened
            self.consolidate_secs += elapsed;
            self.consolidate_count += 1;
        }
        result
    }

    /// Time the merge phase.
    fn time_merge<T>(&mut self, f: impl FnOnce() -> Result<T>) -> Result<T> {
        Self::time(&mut self.merge_secs, f)
    }

    /// Time writing in-memory-only output (no merge needed).
    fn time_write_output(&mut self, f: impl FnOnce() -> Result<()>) -> Result<()> {
        Self::time(&mut self.write_output_secs, f)
    }

    /// Log the final summary.
    #[allow(clippy::cast_precision_loss)]
    fn log_summary(&self, threads: usize) {
        let overall = self.overall_start.map_or(0.0, |s| s.elapsed().as_secs_f64());
        // Guard against division by zero when sort completes in negligible time.
        let overall_nonzero = if overall > 0.0 { overall } else { f64::EPSILON };
        let read_pct = 100.0 * self.read_secs / overall_nonzero;
        let sort_pct = 100.0 * self.sort_secs / overall_nonzero;
        let spill_pct = 100.0 * self.spill_write_secs / overall_nonzero;
        let spill_count = self.spill_count;
        let read_secs = self.read_secs;
        let sort_secs = self.sort_secs;
        let spill_secs = self.spill_write_secs;

        info!("=== Sort Phase Timing ===");
        info!("  Read + decompress: {read_secs:.1}s ({read_pct:.0}%)");
        info!("  In-memory sort:    {sort_secs:.1}s ({sort_pct:.0}%) [{spill_count} spills]");
        let spill_mb = self.total_spill_bytes as f64 / (1024.0 * 1024.0);
        info!(
            "  Spill write:       {spill_secs:.1}s ({spill_pct:.0}%) [{spill_count} writes, {spill_mb:.1} MB total]"
        );
        if self.consolidate_count > 0 {
            let cons_secs = self.consolidate_secs;
            let cons_pct = 100.0 * cons_secs / overall_nonzero;
            let cons_count = self.consolidate_count;
            info!("  Consolidation:     {cons_secs:.1}s ({cons_pct:.0}%) [{cons_count} merges]");
        }
        if self.merge_secs > 0.0 {
            let merge_secs = self.merge_secs;
            let merge_pct = 100.0 * merge_secs / overall_nonzero;
            info!("  K-way merge:       {merge_secs:.1}s ({merge_pct:.0}%)");
        }
        if self.write_output_secs > 0.0 {
            let write_secs = self.write_output_secs;
            let write_pct = 100.0 * write_secs / overall_nonzero;
            info!("  Write output:      {write_secs:.1}s ({write_pct:.0}%)");
        }
        info!("  Total wall clock:  {overall:.1}s");
        info!("  Threads: {threads}");
        info!("=========================");
    }
}

// ============================================================================
// Library Lookup for Template-Coordinate Sort
// ============================================================================

/// Deterministic hasher for cell barcode hashing in template-coordinate sort.
///
/// Uses arbitrary fixed seeds so that hash values are reproducible across runs.
#[must_use]
pub fn cb_hasher() -> ahash::RandomState {
    // Arbitrary fixed seeds — chosen for uniqueness, not cryptographic strength.
    ahash::RandomState::with_seeds(
        0xa1b2_c3d4_e5f6_0718,
        0x9182_7364_5546_3728,
        0xfede_dcba_0987_6543,
        0x0011_2233_4455_6677,
    )
}

/// Maps read group ID -> library ordinal for O(1) comparison.
///
/// Pre-computes ordinals by sorting library names alphabetically.
/// Empty/unknown library sorts first (ordinal 0).
#[derive(Clone)]
pub struct LibraryLookup {
    /// RG ID -> library ordinal
    rg_to_ordinal: HashMap<Vec<u8>, u32>,
    /// Deterministic hasher for read name hashing, constructed once for reuse.
    hasher: ahash::RandomState,
}

impl LibraryLookup {
    /// Build lookup from BAM header.
    #[must_use]
    #[allow(clippy::cast_possible_truncation)]
    pub fn from_header(header: &Header) -> Self {
        // Collect all unique library names from read groups
        let mut libraries: Vec<String> = header
            .read_groups()
            .iter()
            .filter_map(|(_, rg)| {
                rg.other_fields().get(&rg_tag::LIBRARY).map(std::string::ToString::to_string)
            })
            .collect();

        // Sort alphabetically and deduplicate
        libraries.sort();
        libraries.dedup();

        // Build library name -> ordinal mapping
        // Empty string gets ordinal 0, then libraries in sorted order
        let mut lib_to_ordinal: HashMap<String, u32> = HashMap::new();
        lib_to_ordinal.insert(String::new(), 0);
        for (i, lib) in libraries.iter().enumerate() {
            lib_to_ordinal.insert(lib.clone(), (i + 1) as u32);
        }

        // Build RG ID -> ordinal mapping
        let rg_to_ordinal: HashMap<Vec<u8>, u32> = header
            .read_groups()
            .iter()
            .map(|(id, rg)| {
                let lib = rg
                    .other_fields()
                    .get(&rg_tag::LIBRARY)
                    .map(std::string::ToString::to_string)
                    .unwrap_or_default();
                let ordinal = *lib_to_ordinal.get(&lib).unwrap_or(&0);
                (id.to_vec(), ordinal)
            })
            .collect();

        // Arbitrary fixed seeds — chosen for uniqueness, not cryptographic strength.
        let hasher = ahash::RandomState::with_seeds(
            0x517c_c1b7_2722_0a95,
            0x1234_5678_90ab_cdef,
            0xfedc_ba98_7654_3210,
            0x0123_4567_89ab_cdef,
        );

        Self { rg_to_ordinal, hasher }
    }

    /// Hash a read name deterministically.
    #[inline]
    #[must_use]
    pub fn hash_name(&self, name: &[u8]) -> u64 {
        self.hasher.hash_one(name)
    }

    /// Get library ordinal for a record (from RG tag in aux data).
    ///
    /// Only used in tests; production code uses [`ordinal_from_rg`](Self::ordinal_from_rg)
    /// with pre-extracted RG bytes from the single-pass aux scan.
    #[cfg(test)]
    #[must_use]
    pub fn get_ordinal(&self, bam: &[u8]) -> u32 {
        fgumi_raw_bam::RawRecordView::new(bam)
            .tags()
            .find_string(SamTag::RG)
            .and_then(|rg| self.rg_to_ordinal.get(rg))
            .copied()
            .unwrap_or(0)
    }

    /// Get library ordinal from pre-extracted RG tag bytes.
    #[inline]
    #[must_use]
    pub fn ordinal_from_rg(&self, rg: Option<&[u8]>) -> u32 {
        rg.and_then(|rg| self.rg_to_ordinal.get(rg)).copied().unwrap_or(0)
    }

    /// Number of distinct library ordinals realized by the header's read groups.
    ///
    /// Used to decide whether the template-key `tertiary` library lane can vary:
    /// `<= 1` means every header read group maps to the same library ordinal, so
    /// the lane is provisionally constant (decode-time verify still backstops the
    /// case where a record carries an RG absent from the header → ordinal 0).
    /// Returns 0 when the header declares no read groups.
    #[must_use]
    pub fn distinct_header_ordinals(&self) -> usize {
        self.rg_to_ordinal.values().copied().collect::<HashSet<u32>>().len()
    }
}

/// Number of records to prefetch per chunk during merge.
/// Larger buffer reduces I/O latency impact during merge.
const MERGE_PREFETCH_SIZE: usize = 1024;

/// Maximum number of temp files before consolidation (like samtools).
/// When this limit is reached, oldest files are merged to reduce file count.
const DEFAULT_MAX_TEMP_FILES: usize = 64;

/// Working estimate of the raw BAM bytes per template-coordinate record (excluding the
/// inline header and the `TemplateRecordRef<K>` index entry).  Used by the capacity
/// estimator in `sort_template_coordinate_impl` to split `effective_initial_capacity`
/// between the data arena and the ref index.
const EST_BAM_BYTES_PER_TEMPLATE_RECORD: usize = 250;

/// Counting semaphore for limiting concurrent chunk reader I/O.
/// Pre-filled with N tokens; readers acquire before decompressing, release after.
pub(crate) type ChunkReaderSemaphore = (Sender<()>, Receiver<()>);

/// Create a counting semaphore that allows `threads` concurrent readers.
pub(crate) fn make_reader_semaphore(threads: usize) -> Arc<ChunkReaderSemaphore> {
    let limit = threads.max(1);
    let (tx, rx) = bounded(limit);
    for _ in 0..limit {
        tx.send(()).expect("semaphore channel must not be disconnected during initialization");
    }
    Arc::new((tx, rx))
}

// ============================================================================
// Generic Keyed Temp File I/O (works with any RawSortKey)
// ============================================================================
//
// Stores pre-computed sort keys alongside each record for O(1) merge comparisons.
// Format: [key: serialized][len: 4 bytes][record: len bytes] per record

use std::marker::PhantomData;

/// Wrapper for temp chunk writers supporting both raw and compressed output.
enum ChunkWriterInner {
    /// Uncompressed raw output (fastest).
    Raw(BufWriter<std::fs::File>),
    /// Single-threaded BGZF-compressed output.
    SingleThreaded(BgzfWriter<BufWriter<std::fs::File>>),
    /// Multi-threaded BGZF-compressed output (faster for large chunks).
    MultiThreaded(MultithreadedWriter<BufWriter<std::fs::File>>),
}

impl Write for ChunkWriterInner {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            ChunkWriterInner::Raw(w) => w.write(buf),
            ChunkWriterInner::SingleThreaded(w) => w.write(buf),
            ChunkWriterInner::MultiThreaded(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            ChunkWriterInner::Raw(w) => w.flush(),
            ChunkWriterInner::SingleThreaded(w) => w.flush(),
            ChunkWriterInner::MultiThreaded(w) => w.flush(),
        }
    }
}

impl ChunkWriterInner {
    fn finish(self) -> Result<()> {
        match self {
            ChunkWriterInner::Raw(mut w) => {
                w.flush()?;
                Ok(())
            }
            ChunkWriterInner::SingleThreaded(w) => {
                w.finish()?;
                Ok(())
            }
            ChunkWriterInner::MultiThreaded(mut w) => {
                w.finish()?;
                Ok(())
            }
        }
    }
}

/// Generic writer for keyed temp chunks with pre-computed sort keys.
///
/// Works with any type implementing `RawSortKey`.
/// Supports optional BGZF compression for reduced disk usage.
pub(crate) struct GenericKeyedChunkWriter<K: RawSortKey> {
    writer: ChunkWriterInner,
    _marker: PhantomData<K>,
}

impl<K: RawSortKey> GenericKeyedChunkWriter<K> {
    /// Create a new keyed chunk writer with optional compression.
    ///
    /// - `compression_level` 0 = uncompressed (fastest, uses most disk).
    /// - `compression_level` > 0 = BGZF compression at specified level.
    /// - `threads` > 1 enables multi-threaded compression.
    ///
    /// # Errors
    ///
    /// Returns an error if the output file cannot be created.
    ///
    /// # Panics
    ///
    /// Panics if `threads` is greater than 1 but `NonZero::new` receives zero.
    pub fn create(path: &Path, compression_level: u32, threads: usize) -> Result<Self> {
        let file = std::fs::File::create(path)?;
        let buf = BufWriter::with_capacity(256 * 1024, file);

        let writer = if compression_level == 0 {
            ChunkWriterInner::Raw(buf)
        } else if threads > 1 {
            // Use multi-threaded BGZF for faster compression
            let worker_count = NonZero::new(threads).expect("threads > 1");
            let mut builder =
                multithreaded_writer::Builder::default().set_worker_count(worker_count);
            #[allow(clippy::cast_possible_truncation)]
            if let Some(level) = CompressionLevel::new(compression_level as u8) {
                builder = builder.set_compression_level(level);
            }
            ChunkWriterInner::MultiThreaded(builder.build_from_writer(buf))
        } else {
            // Single-threaded BGZF with specified compression level
            #[allow(clippy::cast_possible_truncation)]
            let level = CompressionLevel::new(compression_level as u8).unwrap_or_else(|| {
                CompressionLevel::new(6).expect("compression level 6 is always valid")
            });
            let writer = noodles_bgzf::io::writer::Builder::default()
                .set_compression_level(level)
                .build_from_writer(buf);
            ChunkWriterInner::SingleThreaded(writer)
        };

        Ok(Self { writer, _marker: PhantomData })
    }

    /// Write a keyed record.
    ///
    /// When `K::EMBEDDED_IN_RECORD` is true, the key is embedded in the record
    /// bytes so only the record is written. Otherwise, the key prefix is written
    /// followed by the record.
    ///
    /// # Errors
    ///
    /// Returns an error if writing to the underlying writer fails.
    #[inline]
    pub fn write_record(&mut self, key: &K, record: &[u8]) -> Result<()> {
        if !K::EMBEDDED_IN_RECORD {
            key.write_to(&mut self.writer)?;
        }
        // The merge reader reads this 4-byte length prefix back; a checked cast
        // fails loud before a truncated prefix could desync the spill stream.
        let record_len = u32::try_from(record.len())
            .map_err(|_| anyhow::anyhow!("BAM record too large ({} bytes)", record.len()))?;
        self.writer.write_all(&record_len.to_le_bytes())?;
        self.writer.write_all(record)?;
        Ok(())
    }

    /// Finish writing and flush.
    ///
    /// # Errors
    ///
    /// Returns an error if flushing the writer fails.
    pub fn finish(self) -> Result<()> {
        self.writer.finish()
    }
}

/// Result of reading a keyed record from a chunk: `Ok(Some(...))` for a record,
/// `Ok(None)` for EOF, or `Err(...)` for an I/O error.
type ChunkReadResult<K> = Result<Option<(K, Vec<u8>)>>;

/// Read exactly `buf.len()` bytes, distinguishing clean EOF from truncation.
///
/// Returns `Ok(true)` when `buf` is fully filled, `Ok(false)` when zero bytes are
/// available (clean EOF), and `Err` for any partial read or I/O error.
pub(crate) fn read_exact_or_eof<R: Read>(reader: &mut R, buf: &mut [u8]) -> std::io::Result<bool> {
    let mut offset = 0;
    while offset < buf.len() {
        match reader.read(&mut buf[offset..]) {
            Ok(0) => {
                return if offset == 0 {
                    Ok(false) // clean EOF
                } else {
                    Err(std::io::Error::new(
                        std::io::ErrorKind::UnexpectedEof,
                        format!("truncated chunk: read {} of {} bytes", offset, buf.len()),
                    ))
                };
            }
            Ok(n) => offset += n,
            Err(e) if e.kind() == std::io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    Ok(true)
}

/// Generic reader for keyed temp chunks with background prefetching.
///
/// Works with any type implementing `RawSortKey`.
/// Auto-detects BGZF compression via magic bytes.
pub(crate) struct GenericKeyedChunkReader<K: RawSortKey + 'static> {
    receiver: Receiver<ChunkReadResult<K>>,
    /// Return channel for empty buffers — the consumer sends its old buffer
    /// back so the producer can reuse the allocation instead of allocating.
    buf_return: Sender<Vec<u8>>,
    _handle: JoinHandle<()>,
}

impl<K: RawSortKey + 'static> GenericKeyedChunkReader<K> {
    /// Open a keyed chunk file for reading with background prefetching.
    /// Auto-detects BGZF/gzip compression via magic bytes (0x1f 0x8b).
    ///
    /// An optional `concurrency_limit` semaphore can be provided to cap the number
    /// of reader threads actively performing I/O. Readers acquire a token before
    /// reading a batch of records and release it after sending.
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be opened.
    #[allow(clippy::unnecessary_wraps)]
    pub fn open(path: &Path, concurrency_limit: Option<Arc<ChunkReaderSemaphore>>) -> Result<Self> {
        let (tx, rx) = bounded(MERGE_PREFETCH_SIZE);
        let (buf_tx, buf_rx) = bounded::<Vec<u8>>(MERGE_PREFETCH_SIZE);
        let path = path.to_path_buf();

        let handle = thread::spawn(move || {
            let file = match std::fs::File::open(&path) {
                Ok(f) => f,
                Err(e) => {
                    let _ = tx.send(Err(anyhow::anyhow!(
                        "Failed to open keyed chunk {}: {e}",
                        path.display()
                    )));
                    return;
                }
            };
            let mut buf_reader = BufReader::with_capacity(2 * 1024 * 1024, file);

            // Detect magic bytes at file start:
            //   BGZF/gzip  : 0x1f 0x8b
            //   ZSP1 spill : ASCII "ZSP1" followed by [u32 LE len][zstd frame]+
            //
            // Use `read_exact_or_eof` so a short `read()` (legal for `Read`)
            // can't truncate `ZSPILL_MAGIC` into a `None` and silently route a
            // zstd spill through the uncompressed path. Clean EOF (empty file)
            // is preserved as `read_n == 0` — same behavior as before.
            let mut magic = [0u8; 4];
            let read_n = match read_exact_or_eof(&mut buf_reader, &mut magic) {
                Ok(true) => 4,
                Ok(false) => 0,
                Err(e) => {
                    let _ = tx.send(Err(anyhow::anyhow!(
                        "Failed to read keyed chunk magic {}: {e}",
                        path.display()
                    )));
                    return;
                }
            };
            let codec = crate::codec::SpillCodec::from_magic(&magic[..read_n]);

            // Seek back to start so each decoder consumes the magic itself.
            if buf_reader.seek(SeekFrom::Start(0)).is_err() {
                let _ = tx
                    .send(Err(anyhow::anyhow!("Failed to seek in keyed chunk {}", path.display())));
                return;
            }

            match codec {
                Some(crate::codec::SpillCodec::Zstd) => {
                    match crate::zspill_stream::ZspillStreamReader::new(buf_reader) {
                        Ok(rdr) => Self::read_records(rdr, tx, buf_rx, concurrency_limit),
                        Err(e) => {
                            let _ =
                                tx.send(Err(anyhow::anyhow!("zstd spill reader open failed: {e}")));
                        }
                    }
                }
                Some(crate::codec::SpillCodec::Bgzf) => {
                    let bgzf_reader = BgzfReader::new(buf_reader);
                    Self::read_records(bgzf_reader, tx, buf_rx, concurrency_limit);
                }
                None => {
                    // Treat as uncompressed (existing test behavior).
                    Self::read_records(buf_reader, tx, buf_rx, concurrency_limit);
                }
            }
        });

        Ok(Self { receiver: rx, buf_return: buf_tx, _handle: handle })
    }

    /// Open a keyed chunk file with pool-based BGZF decompression.
    ///
    /// Instead of decompressing BGZF blocks on the reader thread, raw blocks are
    /// submitted to the shared `SortWorkerPool` for decompression. This matches
    /// Read records from a reader and send them through the channel.
    ///
    /// When a semaphore is provided, reads records in batches of 64: acquires
    /// a token, reads the batch (I/O + decompression), releases the token,
    /// then sends the batch through the channel. This prevents deadlock — the
    /// token is never held during a blocking `tx.send()`.
    #[allow(clippy::needless_pass_by_value)]
    fn read_records<R: Read>(
        mut reader: R,
        tx: crossbeam_channel::Sender<ChunkReadResult<K>>,
        buf_pool: crossbeam_channel::Receiver<Vec<u8>>,
        semaphore: Option<Arc<ChunkReaderSemaphore>>,
    ) {
        const BATCH_SIZE: usize = 64;

        loop {
            // Phase 1: Acquire token and read a batch of records from disk.
            if let Some(ref sem) = semaphore {
                let _ = sem.1.recv();
            }

            let mut batch: Vec<(K, Vec<u8>)> = Vec::with_capacity(BATCH_SIZE);
            let mut eof = false;
            let mut read_error: Option<String> = None;

            for _ in 0..BATCH_SIZE {
                if K::EMBEDDED_IN_RECORD {
                    // Keyless format: read record, extract key from BAM bytes.
                    let mut len_buf = [0u8; 4];
                    match read_exact_or_eof(&mut reader, &mut len_buf) {
                        Ok(true) => {}
                        Ok(false) => {
                            eof = true;
                            break;
                        }
                        Err(e) => {
                            read_error = Some(format!("Error reading chunk record length: {e}"));
                            break;
                        }
                    }
                    let len = u32::from_le_bytes(len_buf) as usize;
                    let mut record = buf_pool.try_recv().unwrap_or_default();
                    record.clear();
                    record.resize(len, 0);
                    if let Err(e) = reader.read_exact(&mut record) {
                        read_error = Some(format!("Error reading chunk record: {e}"));
                        break;
                    }
                    let key = K::extract_from_record(&record);
                    batch.push((key, record));
                } else {
                    // Keyed format: read key prefix, then record.
                    let key = match K::read_from(&mut reader) {
                        Ok(k) => k,
                        Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => {
                            eof = true;
                            break;
                        }
                        Err(e) => {
                            read_error = Some(format!("Error reading keyed chunk key: {e}"));
                            break;
                        }
                    };

                    let mut len_buf = [0u8; 4];
                    match reader.read_exact(&mut len_buf) {
                        Ok(()) => {}
                        Err(e) => {
                            read_error = Some(format!("Error reading keyed chunk length: {e}"));
                            break;
                        }
                    }
                    let len = u32::from_le_bytes(len_buf) as usize;

                    let mut record = buf_pool.try_recv().unwrap_or_default();
                    record.clear();
                    record.resize(len, 0);
                    if let Err(e) = reader.read_exact(&mut record) {
                        read_error = Some(format!("Error reading keyed chunk record: {e}"));
                        break;
                    }

                    batch.push((key, record));
                }
            }

            // Phase 2: Release token before any blocking channel sends.
            if let Some(ref sem) = semaphore {
                let _ = sem.0.send(());
            }

            // Phase 3: Send the batch through the channel (may block).
            for record in batch {
                if tx.send(Ok(Some(record))).is_err() {
                    return; // Receiver dropped
                }
            }

            if let Some(msg) = read_error {
                let _ = tx.send(Err(anyhow::anyhow!("{msg}")));
                break;
            }

            if eof {
                let _ = tx.send(Ok(None));
                break;
            }
        }
    }

    /// Read the next keyed record from the prefetch buffer into `buf`.
    ///
    /// On success the record bytes are swapped into `buf` and the sort key is
    /// returned. The old contents of `buf` are returned to the producer thread
    /// for reuse, avoiding per-record allocation on the disk path.
    ///
    /// # Errors
    ///
    /// Returns an error if the background reader encountered an I/O error.
    pub fn next_record(&mut self, buf: &mut Vec<u8>) -> Result<Option<K>> {
        match self.receiver.recv() {
            Ok(Ok(Some((key, mut data)))) => {
                std::mem::swap(buf, &mut data);
                // Return the old buffer to the producer for reuse.
                let _ = self.buf_return.try_send(data);
                Ok(Some(key))
            }
            Ok(Ok(None)) => Ok(None),
            Ok(Err(e)) => Err(e),
            // Channel disconnected — the producer thread panicked or was dropped
            // without sending an EOF sentinel. Treat as an error, not clean EOF.
            Err(_) => Err(anyhow::anyhow!("chunk reader thread terminated unexpectedly")),
        }
    }

    /// Try to read the next record without blocking.
    ///
    /// Returns `Some(Ok(Some(...)))` if a record is available, `Some(Ok(None))` if the
    /// stream is exhausted, `Some(Err(...))` on read error, or `None` if no record is
    /// currently available (channel empty).
    pub fn try_next_record(&mut self) -> Option<ChunkReadResult<K>> {
        match self.receiver.try_recv() {
            Ok(result) => Some(result),
            Err(crossbeam_channel::TryRecvError::Disconnected) => {
                Some(Err(anyhow::anyhow!("chunk reader thread terminated unexpectedly")))
            }
            Err(crossbeam_channel::TryRecvError::Empty) => None,
        }
    }
}

/// Container for the in-memory chunks passed into the merge.
///
/// Inline-buffer sorts (coordinate, template) produce
/// `Shared` chunks that share an `Arc<SegmentedBuf>` backing store —
/// zero per-record allocation, one memcpy per record at merge time.
///
/// The queryname streaming path accumulates records upstream as
/// individual `RawRecord`s and produces `Owned` chunks — the
/// per-record allocation is sunk cost, so we preserve the original
/// zero-copy `mem::swap` merge bridge.
pub enum MemorySources<K: RawSortKey + Default + 'static> {
    /// Inline-buffer (coordinate / template) residual chunks sharing an
    /// `Arc<SegmentedBuf>` — zero per-record allocation.
    Shared(Vec<InMemoryChunk<K>>),
    /// Queryname-style residual chunks (each record an owned `RawRecord`).
    Owned(Vec<Vec<(K, fgumi_raw_bam::RawRecord)>>),
}

impl<K: RawSortKey + Default + 'static> MemorySources<K> {
    fn num_non_empty(&self) -> usize {
        match self {
            Self::Shared(chunks) => chunks.iter().filter(|c| !c.is_empty()).count(),
            Self::Owned(chunks) => chunks.iter().filter(|c| !c.is_empty()).count(),
        }
    }

    fn push_into(self, sources: &mut Vec<ChunkSource<K>>) {
        match self {
            Self::Shared(chunks) => {
                for chunk in chunks {
                    if !chunk.is_empty() {
                        sources.push(ChunkSource::Memory { chunk, idx: usize::MAX });
                    }
                }
            }
            Self::Owned(chunks) => {
                for records in chunks {
                    if !records.is_empty() {
                        sources.push(ChunkSource::MemoryOwned { records, idx: usize::MAX });
                    }
                }
            }
        }
    }
}

/// Source for keyed chunks during merge (disk or in-memory).
///
/// Each variant owns its "current record" state internally so the merge loop
/// can borrow the current record's bytes via [`Self::current_bytes`] without
/// copying. Disk-backed variants own a `scratch: Vec<u8>` holding the most
/// recently read record's bytes; `Memory`/`MemoryOwned` borrow directly from
/// their backing store, eliminating the per-record memcpy on the in-memory path.
enum ChunkSource<K: RawSortKey + Default + 'static> {
    /// In-memory chunk from an inline-buffer sort (coordinate / template).
    /// All records borrow their bytes from a shared `Arc<SegmentedBuf>`;
    /// `current_bytes` borrows them directly (zero-copy).
    ///
    /// `idx` is the index of the CURRENT record. Sentinel `usize::MAX` means
    /// "no record loaded yet"; the first `advance` lands on 0.
    Memory { chunk: InMemoryChunk<K>, idx: usize },
    /// In-memory chunk from a queryname-style sort (each record an owned
    /// `RawRecord`). `current_bytes` borrows `records[idx].1` directly.
    ///
    /// `idx` uses the same `usize::MAX` sentinel as `Memory`.
    MemoryOwned { records: Vec<(K, fgumi_raw_bam::RawRecord)>, idx: usize },
    /// Pool-integrated disk source — workers read and decompress, main thread parses.
    PoolDisk {
        source_id: usize,
        /// Bytes for the current record (filled by the most recent `advance`).
        scratch: Vec<u8>,
    },
}

impl<K: RawSortKey + Default + 'static> ChunkSource<K> {
    /// Advance to the next record. Returns its key, or `None` at EOF.
    ///
    /// For `PoolDisk` sources, `consumer` must be `Some`. After this returns
    /// `Some(_)`, callers may invoke [`Self::current_bytes`] to borrow the
    /// record's bytes without copying.
    fn advance(&mut self, consumer: Option<&mut MainThreadChunkConsumer<K>>) -> Result<Option<K>> {
        match self {
            ChunkSource::PoolDisk { source_id, scratch } => {
                let c = consumer.ok_or_else(|| {
                    anyhow::anyhow!(
                        "PoolDisk source (id {source_id}) requires a MainThreadChunkConsumer \
                         but none was provided — this is a bug in the sort pipeline"
                    )
                })?;
                c.next_record(*source_id, scratch)
            }
            ChunkSource::Memory { chunk, idx } => {
                // Sentinel `usize::MAX` → first advance lands on 0.
                let next = if *idx == usize::MAX { 0 } else { *idx + 1 };
                if next < chunk.len() {
                    let key = chunk.take_key(next);
                    *idx = next;
                    Ok(Some(key))
                } else {
                    Ok(None)
                }
            }
            ChunkSource::MemoryOwned { records, idx } => {
                let next = if *idx == usize::MAX { 0 } else { *idx + 1 };
                if next < records.len() {
                    let key = std::mem::take(&mut records[next].0);
                    *idx = next;
                    Ok(Some(key))
                } else {
                    Ok(None)
                }
            }
        }
    }

    /// Bytes for the current record (the one whose key was returned by the
    /// most recent successful [`Self::advance`]).
    ///
    /// Caller MUST have received `Some(_)` from a prior `advance`. The merge
    /// loops uphold this: the init loop only retains a source whose first
    /// `advance` returned `Some`, and the main loop only calls this on the
    /// active tree winner — so the `PoolDisk` empty-scratch case is
    /// unreachable.
    fn current_bytes(&self) -> &[u8] {
        match self {
            ChunkSource::PoolDisk { scratch, .. } => scratch,
            ChunkSource::Memory { chunk, idx } => {
                debug_assert!(*idx != usize::MAX, "current_bytes called before advance");
                chunk.record_bytes(*idx)
            }
            ChunkSource::MemoryOwned { records, idx } => {
                debug_assert!(*idx != usize::MAX, "current_bytes called before advance");
                &records[*idx].1[..]
            }
        }
    }
}

// ============================================================================
// MainThreadChunkConsumer — pool-integrated merge reader
// ============================================================================
//
// The pool's worker threads read raw BGZF blocks from disk, decompress them,
// and insert the decompressed blocks into per-file `ReorderBuffer`s held in
// each `Phase2FileState`. The main thread (the merge loop) holds a snapshot of
// those file states and pulls one decompressed block at a time per source as
// the loser tree advances. There is no global queue and no per-source
// reorder/buffering layer in the consumer itself — the per-file
// `Phase2FileState.decompressed` reorder buffer IS the per-source buffer.
//
// Backpressure is enforced inside the worker pool: when a per-file reorder
// buffer reaches `PHASE2_DECOMP_CAP`, workers stop pulling new raw blocks for
// that file (with a deadlock-free escape hatch for the gap-filling serial).

/// Per-source byte-stream parser state owned by the main thread.
///
/// As blocks are pulled from `Phase2FileState.decompressed`, the bytes are
/// stashed in `current_buf` and consumed left-to-right by the record parser.
struct SourceParserState {
    /// Current decompressed block being consumed.
    current_buf: Vec<u8>,
    /// Read position within `current_buf`.
    current_pos: usize,
}

impl SourceParserState {
    fn new() -> Self {
        Self { current_buf: Vec::new(), current_pos: 0 }
    }

    fn remaining(&self) -> usize {
        self.current_buf.len() - self.current_pos
    }
}

/// Reads from per-file decompressed-block reorder buffers and presents records
/// to the merge loop.
///
/// The main thread drives all record consumption; no threads are spawned here.
/// Sort-pool workers do the disk reads and BGZF decompression in parallel and
/// publish results into per-file `Phase2FileState.decompressed` reorder
/// buffers, which this consumer drains in serial order.
///
/// # Type Parameter
///
/// `K` is the sort key type (`RawCoordinateKey`, `TemplateKey`, etc.).
pub(crate) struct MainThreadChunkConsumer<K: RawSortKey + 'static> {
    /// Snapshot of the pool's Phase 2 file vector. Indexed by `source_id`.
    files: Arc<Vec<crate::worker_pool::Phase2FileState>>,
    /// Per-source parser state.
    parser_state: Vec<SourceParserState>,
    /// Set by a pool worker when BGZF decompression of a chunk block fails.
    decompression_error: std::sync::Arc<std::sync::atomic::AtomicBool>,
    /// Set by a pool worker when a chunk file I/O read fails.
    chunk_read_error: std::sync::Arc<std::sync::atomic::AtomicBool>,
    /// Set by `do_shutdown` when a worker thread panicked unexpectedly.
    worker_panicked: std::sync::Arc<std::sync::atomic::AtomicBool>,
    _phantom: std::marker::PhantomData<K>,
}

impl<K: RawSortKey + 'static> MainThreadChunkConsumer<K> {
    /// Create a new consumer for the given pool file snapshot.
    #[must_use]
    pub(crate) fn new(
        files: Arc<Vec<crate::worker_pool::Phase2FileState>>,
        decompression_error: std::sync::Arc<std::sync::atomic::AtomicBool>,
        chunk_read_error: std::sync::Arc<std::sync::atomic::AtomicBool>,
        worker_panicked: std::sync::Arc<std::sync::atomic::AtomicBool>,
    ) -> Self {
        let parser_state = (0..files.len()).map(|_| SourceParserState::new()).collect();
        Self {
            files,
            parser_state,
            decompression_error,
            chunk_read_error,
            worker_panicked,
            _phantom: std::marker::PhantomData,
        }
    }

    /// Get the next record from a specific source.
    ///
    /// Pulls decompressed blocks from the source's per-file reorder buffer as
    /// needed. Parses the next record from the source's byte stream.
    ///
    /// Returns `Ok(Some(key))` with record bytes swapped into `buf`, `Ok(None)` at EOF.
    ///
    /// # Errors
    ///
    /// Returns an error if a record is truncated or if a worker reported an error.
    pub fn next_record(&mut self, source_id: usize, buf: &mut Vec<u8>) -> Result<Option<K>> {
        // Make sure we have data (or detect EOF) before parsing.
        if self.parser_state[source_id].remaining() == 0
            && !self.advance_to_next_block(source_id)?
        {
            return Ok(None);
        }
        self.parse_next_record(source_id, buf)
    }

    /// Pull the next decompressed block for `source_id` from the per-file
    /// reorder buffer, blocking the main thread (`std::thread::park`) until
    /// either a block becomes available, the source drains, or a worker error
    /// is reported.
    ///
    /// Returns `Ok(true)` if a new block was loaded into `current_buf`,
    /// `Ok(false)` if the source has produced all its data, or an error if a
    /// worker reported a fatal failure.
    fn advance_to_next_block(&mut self, source_id: usize) -> Result<bool> {
        let file = &self.files[source_id];
        loop {
            // Try to pop the next-in-order decompressed block.
            {
                let mut guard =
                    file.decompressed.lock().expect("phase2 decompressed mutex poisoned");
                if let Some(data) = guard.try_pop_next() {
                    drop(guard);
                    let st = &mut self.parser_state[source_id];
                    let exhausted = std::mem::replace(&mut st.current_buf, data);
                    st.current_pos = 0;
                    // Recycle the just-exhausted block back to the shared
                    // decompression buffer pool so zstd workers can reuse the
                    // allocation for the next frame (BGZF decompress returns its
                    // own buffers and does not check out, so skip it there). The
                    // initial `current_buf` is a zero-capacity `Vec` that was
                    // never checked out, so don't return it to the pool.
                    if matches!(file.codec, crate::codec::SpillCodec::Zstd)
                        && exhausted.capacity() != 0
                    {
                        file.buffer_pool.checkin(exhausted);
                    }
                    return Ok(true);
                }
            }

            // No block ready. Check error flags first — they take precedence
            // over EOF detection so a single-source sort that fails on its
            // last block surfaces the error rather than silently truncating.
            if self.decompression_error.load(std::sync::atomic::Ordering::Acquire) {
                return Err(anyhow::anyhow!(
                    "BGZF decompression error on chunk blocks (see log for details)"
                ));
            }
            if self.chunk_read_error.load(std::sync::atomic::Ordering::Acquire) {
                return Err(anyhow::anyhow!("I/O error reading chunk file (see log for details)"));
            }
            if self.worker_panicked.load(std::sync::atomic::Ordering::Acquire) {
                return Err(anyhow::anyhow!(
                    "a sort worker thread panicked unexpectedly (see log for details)"
                ));
            }

            // Source produced everything it ever will?
            if file.is_drained() {
                // Recycle the terminal zstd block: nothing replaces `current_buf`
                // after a source drains, so without this the last decompression
                // allocation for each source lingers in `parser_state` until the
                // whole merge ends — one retained buffer per spill file. (BGZF
                // returns its own buffers and never checks out; the initial
                // zero-capacity `Vec` was never checked out either.)
                let st = &mut self.parser_state[source_id];
                if matches!(file.codec, crate::codec::SpillCodec::Zstd)
                    && st.current_buf.capacity() != 0
                {
                    file.buffer_pool.checkin(std::mem::take(&mut st.current_buf));
                    st.current_pos = 0;
                }
                return Ok(false);
            }

            // Park until a worker unparks us. Workers unpark after pushing a
            // decompressed block, after setting reader.eof, and after error
            // flags. The loop re-checks all conditions on wake-up so spurious
            // wake-ups are harmless.
            std::thread::park();
        }
    }

    /// Parse the next record from a source's byte stream.
    ///
    /// Handles the format: for `EMBEDDED_IN_RECORD` keys, reads [len(4)][record(len)].
    /// For keyed format, reads [key][len(4)][record(len)].
    fn parse_next_record(&mut self, source_id: usize, buf: &mut Vec<u8>) -> Result<Option<K>> {
        let mut len_buf = [0u8; 4];

        if K::EMBEDDED_IN_RECORD {
            if !self.read_exact_from_source(source_id, &mut len_buf)? {
                return Ok(None);
            }
            let len = u32::from_le_bytes(len_buf) as usize;

            buf.clear();
            buf.resize(len, 0);
            if !self.read_exact_from_source(source_id, buf)? {
                return Err(anyhow::anyhow!("truncated record in chunk source {source_id}"));
            }
            let key = K::extract_from_record(buf);
            Ok(Some(key))
        } else {
            let Some(key) = self.read_key_from_source::<K>(source_id)? else {
                return Ok(None);
            };

            if !self.read_exact_from_source(source_id, &mut len_buf)? {
                return Err(anyhow::anyhow!("truncated record length in chunk source {source_id}"));
            }
            let len = u32::from_le_bytes(len_buf) as usize;

            buf.clear();
            buf.resize(len, 0);
            if !self.read_exact_from_source(source_id, buf)? {
                return Err(anyhow::anyhow!("truncated record in chunk source {source_id}"));
            }
            Ok(Some(key))
        }
    }

    /// Read exactly `out.len()` bytes from a source into `out`, pulling more
    /// blocks from the per-file reorder buffer as needed.
    ///
    /// Returns `Ok(false)` at clean EOF (zero bytes available), `Ok(true)` on success.
    fn read_exact_from_source(&mut self, source_id: usize, out: &mut [u8]) -> Result<bool> {
        let n = out.len();
        let mut filled = 0;

        while filled < n {
            if self.parser_state[source_id].remaining() == 0
                && !self.advance_to_next_block(source_id)?
            {
                if filled == 0 {
                    return Ok(false);
                }
                return Err(anyhow::anyhow!(
                    "truncated data in chunk source {source_id}: got {filled} of {n} bytes",
                ));
            }

            let st = &mut self.parser_state[source_id];
            let take = (n - filled).min(st.remaining());
            out[filled..filled + take]
                .copy_from_slice(&st.current_buf[st.current_pos..st.current_pos + take]);
            st.current_pos += take;
            filled += take;
        }

        Ok(true)
    }

    /// Read a sort key from a source's byte stream.
    ///
    /// Returns `Ok(None)` at clean EOF.
    fn read_key_from_source<KK: RawSortKey>(&mut self, source_id: usize) -> Result<Option<KK>> {
        let mut adapter = SourceReadAdapter { consumer: self, source_id, bytes_read: 0 };
        match KK::read_from(&mut adapter) {
            Ok(key) => Ok(Some(key)),
            Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => {
                if adapter.bytes_read == 0 {
                    Ok(None)
                } else {
                    Err(anyhow::anyhow!(
                        "truncated key in chunk source {source_id}: \
                         got {n} bytes then EOF",
                        n = adapter.bytes_read
                    ))
                }
            }
            Err(e) => Err(anyhow::anyhow!("error reading key from source {source_id}: {e}")),
        }
    }

    /// Gather probe statistics from the per-file Phase 2 state.
    fn probe_consumer_stats(&self) -> ConsumerProbeStats {
        let mut pending_blocks: u64 = 0;
        let mut pending_bytes: u64 = 0;
        let mut active_sources: u64 = 0;

        for file in self.files.iter() {
            let (blocks, bytes, active) = file.probe_stats();
            pending_blocks += blocks;
            pending_bytes += bytes;
            if active {
                active_sources += 1;
            }
        }

        ConsumerProbeStats {
            current_bytes: 0,
            current_capacity: 0,
            pending_blocks,
            pending_bytes,
            active_sources,
        }
    }
}

/// Adapter that implements `std::io::Read` over a `MainThreadChunkConsumer` source.
///
/// This allows `K::read_from(&mut reader)` to read from the pool-based byte stream.
/// `bytes_read` tracks how many bytes have been consumed so `read_key_from_source` can
/// distinguish a clean EOF (zero bytes seen) from a truncated read (some bytes then EOF).
struct SourceReadAdapter<'a, K: RawSortKey + 'static> {
    consumer: &'a mut MainThreadChunkConsumer<K>,
    source_id: usize,
    bytes_read: usize,
}

impl<K: RawSortKey + 'static> Read for SourceReadAdapter<'_, K> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        // Pull the next block if the current one is exhausted.
        if self.consumer.parser_state[self.source_id].remaining() == 0 {
            match self.consumer.advance_to_next_block(self.source_id) {
                Ok(true) => {}
                Ok(false) => return Ok(0), // clean EOF
                Err(e) => return Err(std::io::Error::other(e.to_string())),
            }
        }

        let st = &mut self.consumer.parser_state[self.source_id];
        let take = buf.len().min(st.remaining());
        buf[..take].copy_from_slice(&st.current_buf[st.current_pos..st.current_pos + take]);
        st.current_pos += take;
        self.bytes_read += take;
        Ok(take)
    }
}

/// Generates unique file paths for chunk and merged temp files.
///
/// Maintains monotonic counters for both chunk files (`chunk_0000.keyed`, ...)
/// and merged files (`merged_0000.keyed`, ...) to prevent naming collisions
/// after consolidation drains entries from the chunk file list.
///
/// When multiple temp directories are supplied via [`TmpDirAllocator`], chunk
/// and merged files are distributed across them in round-robin order.
struct ChunkNamer<'a> {
    alloc: &'a mut TmpDirAllocator,
    chunk_count: usize,
    merge_count: usize,
}

impl<'a> ChunkNamer<'a> {
    fn new(alloc: &'a mut TmpDirAllocator) -> Self {
        Self { alloc, chunk_count: 0, merge_count: 0 }
    }

    /// Returns the next unique chunk file path, drawing from the allocator.
    fn next_chunk_path(&mut self) -> Result<PathBuf> {
        let base = self.alloc.next()?;
        let path = base.join(format!("chunk_{:04}.keyed", self.chunk_count));
        self.chunk_count += 1;
        Ok(path)
    }

    /// Returns the next unique merged file path, drawing from the allocator.
    fn next_merged_path(&mut self) -> Result<PathBuf> {
        let base = self.alloc.next()?;
        let path = base.join(format!("merged_{:04}.keyed", self.merge_count));
        self.merge_count += 1;
        Ok(path)
    }
}

/// A spill write that is finishing in the background.
///
/// Used for pipelining: the I/O thread continues writing while the main thread
/// reads the next batch. The chunk path is stored alongside the handle so it
/// can be pushed to `chunk_files` only after the write completes.
struct PendingSpill {
    handle: crate::pooled_chunk_writer::SpillWriteHandle,
    chunk_path: PathBuf,
}

/// Build `BufferProbeStats` from any buffer implementing `ProbeableBuffer`.
#[allow(clippy::cast_possible_truncation)]
fn probe_stats(buf: &impl ProbeableBuffer) -> BufferProbeStats {
    BufferProbeStats {
        usage: buf.memory_usage() as u64,
        capacity: buf.allocated_capacity() as u64,
        records: buf.len() as u64,
        segments: buf.num_segments() as u64,
    }
}

/// Raw-bytes external sorter for BAM files.
///
/// This sorter uses lazy record parsing to minimize memory usage and avoid
/// re-encoding overhead. It's significantly faster than the RecordBuf-based
/// sorter for large files.
pub struct RawExternalSorter {
    /// Sort order to use.
    sort_order: SortOrder,
    /// Maximum memory to use for in-memory sorting.
    memory_limit: usize,
    /// Temporary directories for spill files.
    ///
    /// When empty, a single directory is created under the system default
    /// temp location. When one or more paths are given, spill files are
    /// distributed across them in free-space-aware round-robin order via
    /// [`TmpDirAllocator`].
    temp_dirs: Vec<PathBuf>,
    /// Number of threads for parallel operations — the base/fallback count used
    /// by both phases unless overridden by `sort_threads` / `merge_threads`.
    threads: usize,
    /// Optional override for the Phase-1 (accumulate/sort/spill) worker count.
    /// `None` ⇒ Phase 1 uses `threads`. Lower it to cede cores to an upstream
    /// producer (e.g. `bwa mem` in a pipeline) during accumulation.
    sort_threads: Option<usize>,
    /// Optional override for the Phase-2 (merge/write) worker count.
    /// `None` ⇒ Phase 2 uses `threads`. The merge/write phase scales with cores,
    /// so this is typically left at (or raised toward) the full core count.
    merge_threads: Option<usize>,
    /// Compression level for output.
    output_compression: u32,
    /// Compression level for temporary chunk files (0 = uncompressed).
    temp_compression: u32,
    /// Codec used for temporary chunk files (Bgzf or Zstd).
    spill_codec: crate::codec::SpillCodec,
    /// Program record info (version, `command_line`) for @PG header.
    pg_info: Option<(String, String)>,
    /// Maximum temp files before consolidation (0 = unlimited).
    max_temp_files: usize,
    /// Cell barcode tag for template-coordinate sort (e.g., `SamTag::CB`).
    /// When `Some`, CB hash is included in sort key for single-cell data.
    cell_tag: Option<SamTag>,
    /// Initial buffer capacity hint (bytes) for pre-allocation.
    ///
    /// Decoupled from `memory_limit` so that auto-detected limits can start with
    /// a modest allocation and let `Vec` grow on demand, while explicit limits
    /// pre-allocate the full budget upfront (preserving prior behavior).
    initial_capacity: Option<usize>,
    /// When true, wrap input in a `PrefetchReader` for async I/O.
    async_reader: bool,
    /// Which optional template-key lanes to retain (template-coordinate only).
    ///
    /// Defaults to [`KeyTypesSpec::Auto`], which provisions the narrowest key
    /// that fits the first record's optional lanes.
    key_types: KeyTypesSpec,
}

/// RAII guard that ensures Phase 2 teardown runs on every exit path between
/// `pool.set_phase(PHASE2)` and the explicit `deactivate()` in the merge loop.
/// Without this, any `?` early-return would leave the pool stuck in PHASE2 with
/// `phase2_files` still published. `deactivate()` drops the consumer (releasing
/// the Arc snapshot of the per-file vector), resets the phase to `LEGACY`, and
/// clears the pool's published file vector — in that order.
struct Phase2Guard<'a, K: RawSortKey + 'static> {
    pool: &'a Arc<SortWorkerPool>,
    consumer: Option<MainThreadChunkConsumer<K>>,
    active: bool,
}

impl<K: RawSortKey + 'static> Phase2Guard<'_, K> {
    fn consumer_mut(&mut self) -> Option<&mut MainThreadChunkConsumer<K>> {
        self.consumer.as_mut()
    }

    fn deactivate(&mut self) {
        if self.active {
            drop(self.consumer.take());
            self.pool.set_phase(crate::worker_pool::phase::LEGACY);
            self.pool.clear_phase2_files();
            self.active = false;
        }
    }
}

impl<K: RawSortKey + 'static> Drop for Phase2Guard<'_, K> {
    fn drop(&mut self) {
        self.deactivate();
    }
}

impl RawExternalSorter {
    /// The configured sort order. Read by the typed-step `SortAndSpill` adapter
    /// to pick the matching `*SortStream` variant.
    #[must_use]
    pub fn sort_order(&self) -> SortOrder {
        self.sort_order
    }

    /// The configured cell-barcode tag (template-coordinate CB hashing), if any.
    /// Read when building a record-input arena strategy so it matches the
    /// block-input template strategy exactly.
    #[must_use]
    pub fn cell_tag_value(&self) -> Option<SamTag> {
        self.cell_tag
    }

    /// The configured `--key-types` narrowing spec (template-coordinate only).
    #[must_use]
    pub fn key_types_spec(&self) -> KeyTypesSpec {
        self.key_types
    }

    /// The configured in-memory sort budget in bytes (per-run seal threshold).
    #[must_use]
    pub fn memory_limit_bytes(&self) -> usize {
        self.memory_limit
    }

    /// The configured sort thread count (per-chunk parallel sort width).
    #[must_use]
    pub fn num_threads(&self) -> usize {
        self.threads
    }

    /// Convert this configured sorter into a lean coordinate buffering sorter
    /// for the P6 `SortBuffer` step.
    ///
    /// Same buffer sizing / rayon pool as the retired `into_coordinate_stream`,
    /// but with **no** `SortWorkerPool`, temp dirs, or spill files — the
    /// `SortBuffer` step only ingests + par-sorts + materializes chunks; the
    /// `CompressSpill` step owns compression and disk I/O.
    ///
    /// # Errors
    ///
    /// Returns an error if `self.sort_order` is not `Coordinate`, the reference
    /// count exceeds `u32::MAX`, or the rayon sort pool cannot be built.
    pub fn into_coordinate_chunk_sorter(
        self,
        header: &Header,
    ) -> Result<crate::chunk_sorter::CoordinateChunkSorter> {
        anyhow::ensure!(
            matches!(self.sort_order, SortOrder::Coordinate),
            "into_coordinate_chunk_sorter requires SortOrder::Coordinate (got {:?})",
            self.sort_order,
        );
        let rayon_pool = self.build_sort_rayon_pool()?;
        let nref = u32::try_from(header.reference_sequences().len())
            .context("reference sequence count exceeds u32::MAX")?;
        let init_cap = self.effective_initial_capacity();
        // Per-record footprint estimate matches `into_coordinate_stream`.
        let estimated_records = init_cap / 240;
        // The data arena comes from the `ArenaPool` (acquired on the first
        // `ensure_arena`), so the buffer starts with NO pre-allocated data
        // segment (0 bytes) — only the refs Vec is pre-sized.
        let buffer = RecordBuffer::with_capacity(estimated_records, 0, nref);
        // Bounded reusable-arena pool (lever-1 RSS fix). Capacity 1 = legacy's
        // one-arena-at-a-time footprint: the next chunk cannot fill until the
        // prior chunk is spilled and its arena returns, bounding peak RSS to
        // ~one arena ≈ base. (Tunable via `--sort-arenas` in a later increment.)
        let arena_pool = crate::arena_pool::ArenaPool::new(1, crate::inline::SORT_SEGMENT_SIZE);
        Ok(crate::chunk_sorter::CoordinateChunkSorter::new(
            buffer,
            rayon_pool,
            self.memory_limit,
            arena_pool,
        ))
    }

    /// Convert this configured sorter into a lean template-coordinate buffering
    /// sorter for the P6 `SortBuffer` step (the template analogue of
    /// [`Self::into_coordinate_chunk_sorter`] — same sizing as the retired
    /// `into_template_coordinate_stream`, minus the pool/temp-dirs).
    ///
    /// # Errors
    ///
    /// Returns an error if `self.sort_order` is not `TemplateCoordinate` or the
    /// rayon sort pool cannot be built.
    pub fn into_template_chunk_sorter(
        self,
        header: &Header,
    ) -> Result<crate::chunk_sorter::TemplateChunkSorter> {
        anyhow::ensure!(
            matches!(self.sort_order, SortOrder::TemplateCoordinate),
            "into_template_chunk_sorter requires SortOrder::TemplateCoordinate (got {:?})",
            self.sort_order,
        );
        let rayon_pool = self.build_sort_rayon_pool()?;
        let lib_lookup = LibraryLookup::from_header(header);
        let header_library_varies = lib_lookup.distinct_header_ordinals() > 1;
        let cb_hasher = cb_hasher();
        // Capacity hints for the lazily-built narrowed-key buffer. The variant
        // (and thus the exact key width) is not known until the first record, so
        // a key-width-agnostic estimate is used (the per-K ref-size difference is
        // immaterial for a capacity hint).
        let init_cap = self.effective_initial_capacity();
        let bytes_per_record = 354;
        let estimated_records = init_cap / bytes_per_record;
        let estimated_data_bytes = init_cap * 86 / 100;
        Ok(crate::chunk_sorter::TemplateChunkSorter::new(
            rayon_pool,
            self.memory_limit,
            estimated_records,
            estimated_data_bytes,
            lib_lookup,
            cb_hasher,
            self.cell_tag,
            self.key_types,
            header_library_varies,
        ))
    }

    /// Create a new raw external sorter with the given sort order.
    #[must_use]
    pub fn new(sort_order: SortOrder) -> Self {
        Self {
            sort_order,
            memory_limit: 512 * 1024 * 1024, // 512 MB default
            temp_dirs: Vec::new(),
            threads: 1,
            sort_threads: None,
            merge_threads: None,
            output_compression: 6,
            temp_compression: 1, // Default: fast compression
            spill_codec: crate::codec::SpillCodec::default(),
            pg_info: None,
            max_temp_files: DEFAULT_MAX_TEMP_FILES,
            cell_tag: None,
            initial_capacity: None,
            async_reader: false,
            key_types: KeyTypesSpec::default(),
        }
    }

    /// Set the memory limit for in-memory sorting.
    #[must_use]
    pub fn memory_limit(mut self, limit: usize) -> Self {
        self.memory_limit = limit;
        self
    }

    /// Set a single temporary directory for spill files.
    ///
    /// Equivalent to calling [`Self::temp_dirs`] with a single-element vector.
    #[must_use]
    pub fn temp_dir(mut self, path: PathBuf) -> Self {
        self.temp_dirs = vec![path];
        self
    }

    /// Set multiple temporary directories for spill files.
    ///
    /// Spill files are distributed across the supplied directories in
    /// free-space-aware round-robin order. Passing an empty vector falls
    /// back to a single directory under the system temp location.
    #[must_use]
    pub fn temp_dirs(mut self, paths: Vec<PathBuf>) -> Self {
        self.temp_dirs = paths;
        self
    }

    /// Set the number of threads.
    #[must_use]
    pub fn threads(mut self, threads: usize) -> Self {
        self.threads = threads;
        self
    }

    /// Set the Phase-1 (accumulate/sort/spill) worker count. Defaults to
    /// [`threads`](Self::threads) when unset. Lower it to cede cores to an
    /// upstream producer in a pipeline.
    ///
    /// Honored by [`sort`](Self::sort) (via the worker-pool active cap) and by
    /// the streaming `into_*_stream` entry points (which size their Phase-1-only
    /// pool via `phase1_threads()`).
    #[must_use]
    pub fn sort_threads(mut self, n: usize) -> Self {
        self.sort_threads = Some(n);
        self
    }

    /// Set the Phase-2 (merge/write) worker count. Defaults to
    /// [`threads`](Self::threads) when unset.
    ///
    /// Honored by [`sort`](Self::sort) (via the worker-pool active cap). The
    /// streaming `into_*_stream` entry points do not read this value — their
    /// Phase 2 (decompress/merge) runs on the framework scheduler, so the chain
    /// builder caps streaming Phase-2 concurrency at the pipeline layer
    /// (`SortSpillDecompress::with_max_concurrency`) instead.
    #[must_use]
    pub fn merge_threads(mut self, n: usize) -> Self {
        self.merge_threads = Some(n);
        self
    }

    /// Effective Phase-1 worker count: `sort_threads` if set, else `threads`.
    ///
    /// Public so the streaming arena front (`SortBuffer::from_sorter`) sizes its
    /// per-chunk sort with the resolved `--sort-threads` override rather than the
    /// raw base `--threads` (`num_threads`), which silently ignores the override.
    #[must_use]
    pub fn phase1_threads(&self) -> usize {
        self.sort_threads.unwrap_or(self.threads).max(1)
    }

    /// Effective Phase-2 worker count: `merge_threads` if set, else `threads`.
    fn phase2_threads(&self) -> usize {
        self.merge_threads.unwrap_or(self.threads).max(1)
    }

    /// Set the output compression level.
    #[must_use]
    pub fn output_compression(mut self, level: u32) -> Self {
        self.output_compression = level;
        self
    }

    /// Set compression level for temporary chunk files.
    ///
    /// For [`SpillCodec::Bgzf`], level 0 disables compression (fastest, uses
    /// most disk space). For [`SpillCodec::Zstd`] there is no level-0 "stored"
    /// mode; pass a level >= 1 or [`Self::sort`] will reject the combination.
    /// Level 1 (default) provides fast compression with reasonable space savings.
    /// Higher levels provide better compression but are slower.
    ///
    /// [`SpillCodec::Bgzf`]: crate::codec::SpillCodec::Bgzf
    /// [`SpillCodec::Zstd`]: crate::codec::SpillCodec::Zstd
    #[must_use]
    pub fn temp_compression(mut self, level: u32) -> Self {
        self.temp_compression = level;
        self
    }

    /// Set the codec used for temporary chunk files.
    ///
    /// Defaults to [`SpillCodec::Zstd`](crate::codec::SpillCodec::Zstd) which is
    /// significantly faster than BGZF at comparable ratios for BAM-record data.
    #[must_use]
    pub fn spill_codec(mut self, codec: crate::codec::SpillCodec) -> Self {
        self.spill_codec = codec;
        self
    }

    /// Set program record info for @PG header entry.
    #[must_use]
    pub fn pg_info(mut self, version: String, command_line: String) -> Self {
        self.pg_info = Some((version, command_line));
        self
    }

    /// Set maximum temp files before consolidation.
    ///
    /// When the number of temp files exceeds this limit, the oldest files
    /// are merged together to reduce the count. Set to 0 for unlimited.
    /// Default is 64 (matching samtools).
    #[must_use]
    pub fn max_temp_files(mut self, max: usize) -> Self {
        self.max_temp_files = max;
        self
    }

    /// Set the cell barcode tag for template-coordinate sort.
    ///
    /// When set, the CB hash is included in the sort key so that templates
    /// from different cells at the same locus are not interleaved.
    #[must_use]
    pub fn cell_tag(mut self, tag: SamTag) -> Self {
        self.cell_tag = Some(tag);
        self
    }

    /// Set the initial buffer capacity hint (bytes).
    ///
    /// When set, buffer pre-allocation uses this value instead of `memory_limit`.
    /// This avoids huge upfront allocations when auto-detecting memory, while
    /// still allowing the buffer to grow up to `memory_limit` before spilling.
    #[must_use]
    pub fn initial_capacity(mut self, bytes: usize) -> Self {
        self.initial_capacity = Some(bytes);
        self
    }

    /// Set which optional template-key lanes to retain.
    ///
    /// Only affects the template-coordinate sort order. The default
    /// ([`KeyTypesSpec::Auto`]) provisions the narrowest key consistent with the
    /// first record's optional lanes; the decode-time verify guarantees the
    /// dropped lanes are constant across all records.
    #[must_use]
    pub fn key_types(mut self, spec: KeyTypesSpec) -> Self {
        self.key_types = spec;
        self
    }

    /// Enable/disable the async prefetch reader on input.
    ///
    /// When enabled, the input BAM is wrapped in a `PrefetchReader` before the
    /// BGZF layer, which overlaps block I/O with decompression.
    #[must_use]
    pub fn async_reader(mut self, enabled: bool) -> Self {
        self.async_reader = enabled;
        self
    }

    /// Returns the effective initial capacity for buffer pre-allocation.
    ///
    /// Uses `initial_capacity` if set, otherwise falls back to `memory_limit`.
    fn effective_initial_capacity(&self) -> usize {
        self.initial_capacity.unwrap_or(self.memory_limit).min(self.memory_limit)
    }

    /// Build a rayon thread pool sized to `self.threads`.
    ///
    /// The sort path uses `par_sort` and friends at several points. Rayon's
    /// global pool defaults to `num_cpus::get()`, which silently violates the
    /// user's `--threads` contract on machines where more physical cores are
    /// available. Every rayon call site is wrapped with `pool.install(...)`
    /// so that `rayon::current_num_threads()` returns `self.threads` and
    /// fan-out is bounded to the requested thread count.
    ///
    /// Oversubscription with the `SortWorkerPool` is not a concern because
    /// every call site is preceded by `drain_pending_spill`, which joins the
    /// prior chunk's I/O thread and therefore guarantees all sort workers are
    /// idle at the moment rayon fans out.
    fn build_sort_rayon_pool(&self) -> Result<rayon::ThreadPool> {
        rayon::ThreadPoolBuilder::new()
            .num_threads(self.phase1_threads())
            .thread_name(|i| format!("fgumi-sort-rayon-{i}"))
            .build()
            .map_err(|e| anyhow::anyhow!("failed to build rayon sort pool: {e}"))
    }

    /// Consolidate temp files if we've exceeded the limit.
    /// Wait for a pending spill to complete and, if one was present, run consolidation.
    ///
    /// Shared by all four sort functions: the in-loop "wait before next spill" drain and
    /// the post-loop "drain before merge" drain are identical aside from the key type.
    fn drain_pending_spill<K: RawSortKey + Default + 'static>(
        &self,
        pending: &mut Option<PendingSpill>,
        chunk_files: &mut Vec<PathBuf>,
        stats: &mut RawSortStats,
        timer: &mut SortPhaseTimer,
        namer: &mut ChunkNamer<'_>,
        pool: &std::sync::Arc<crate::worker_pool::SortWorkerPool>,
    ) -> Result<()> {
        if let Some(prev) = pending.take() {
            prev.handle.wait()?;
            timer.record_spill_size(&prev.chunk_path);
            chunk_files.push(prev.chunk_path);
            stats.chunks_written += 1;

            timer.time_consolidate(|| {
                self.maybe_consolidate_temp_files::<K>(chunk_files, namer, pool)
            })?;
        }
        Ok(())
    }

    ///
    /// Merges the oldest half of temp files into a single new file to reduce
    /// the total count while maintaining sort order.
    fn maybe_consolidate_temp_files<K: RawSortKey + Default + 'static>(
        &self,
        chunk_files: &mut Vec<PathBuf>,
        namer: &mut ChunkNamer<'_>,
        pool: &Arc<SortWorkerPool>,
    ) -> Result<()> {
        use crate::loser_tree::LoserTree;

        if self.max_temp_files == 0 || chunk_files.len() < self.max_temp_files {
            return Ok(());
        }

        // Need at least 2 files to consolidate meaningfully
        if self.max_temp_files < 2 {
            return Ok(());
        }

        // Merge oldest half of files into one (at least 2)
        let n_to_merge = (self.max_temp_files / 2).max(2).min(chunk_files.len());
        let files_to_merge: Vec<PathBuf> = chunk_files.drain(..n_to_merge).collect();

        debug!(
            "Consolidating {} temp files into 1 (total was {})...",
            n_to_merge,
            n_to_merge + chunk_files.len()
        );

        // Create merged output file
        let merged_path = namer.next_merged_path()?;

        // Open readers with semaphore to cap concurrent I/O. Consolidation runs
        // during Phase 1 (from `drain_pending_spill`, before the Phase-2 worker
        // switch), so bound the readers by the Phase-1 thread count to match the
        // rest of the Phase-1 sizing rather than the base thread count.
        let sem = make_reader_semaphore(self.phase1_threads());
        let mut readers: Vec<GenericKeyedChunkReader<K>> = files_to_merge
            .iter()
            .map(|p| GenericKeyedChunkReader::<K>::open(p, Some(Arc::clone(&sem))))
            .collect::<Result<Vec<_>>>()?;

        // Use pooled writer for parallel compression during consolidation.
        // Match the pool's spill codec so the merged file is the same format
        // as the per-chunk spill files it consolidates.
        let mut writer =
            PooledChunkWriter::<K>::new(Arc::clone(pool), &merged_path, pool.spill_codec())?;

        // Initialize loser tree with first record from each reader
        let mut initial_keys: Vec<K> = Vec::with_capacity(readers.len());
        let mut records: Vec<Vec<u8>> = Vec::with_capacity(readers.len());
        let mut source_map: Vec<usize> = Vec::with_capacity(readers.len());

        for (reader_idx, reader) in readers.iter_mut().enumerate() {
            let mut record = Vec::new();
            if let Some(key) = reader.next_record(&mut record)? {
                initial_keys.push(key);
                records.push(record);
                source_map.push(reader_idx);
            }
        }

        if initial_keys.is_empty() {
            writer.finish()?;
            // Insert at beginning to preserve stable order
            chunk_files.insert(0, merged_path);
            // Clean up old files
            for path in &files_to_merge {
                let _ = std::fs::remove_file(path);
            }
            return Ok(());
        }

        let mut tree = LoserTree::new(initial_keys);

        while tree.winner_is_active() {
            let winner = tree.winner();
            let reader_idx = source_map[winner];
            writer.write_record(tree.winner_key(), &records[winner])?;

            if let Some(next_key) = readers[reader_idx].next_record(&mut records[winner])? {
                tree.replace_winner(next_key);
            } else {
                tree.remove_winner();
            }
        }

        writer.finish()?;

        // Insert merged file at the beginning to preserve stable order for equal keys.
        // The merged file contains the oldest records, so it should be processed first.
        chunk_files.insert(0, merged_path);

        // Clean up old files
        for path in &files_to_merge {
            let _ = std::fs::remove_file(path);
        }

        debug!("Consolidation complete, {} temp files remain", chunk_files.len());

        Ok(())
    }

    /// Reject `temp_compression == 0` combined with `SpillCodec::Zstd`: zstd has
    /// no level-0 "stored" mode, so silently remapping to 1 would surprise
    /// callers who set `temp_compression(0)` expecting an uncompressed spill
    /// (which works for BGZF). Called by `sort()` and by every streaming entry
    /// point (`into_*_stream`) so neither the standalone nor the fused path can
    /// reach a misconfigured pool.
    ///
    /// # Errors
    ///
    /// Returns an error if `temp_compression == 0` and the spill codec is zstd.
    fn ensure_spill_codec_compat(&self) -> Result<()> {
        anyhow::ensure!(
            !(self.temp_compression == 0
                && matches!(self.spill_codec, crate::codec::SpillCodec::Zstd)),
            "temp_compression=0 is only supported with SpillCodec::Bgzf; \
             zstd does not have an uncompressed mode. Pass \
             spill_codec(SpillCodec::Bgzf) to keep level-0 spill, or pick a \
             zstd level >= 1."
        );
        Ok(())
    }

    /// Sort a BAM file using raw-bytes approach.
    ///
    /// # Errors
    ///
    /// Returns an error if reading, sorting, or writing the BAM file fails.
    pub fn sort(&self, input: &Path, output: &Path) -> Result<RawSortStats> {
        // Mirror the CLI guard so direct RawExternalSorter callers cannot reach
        // a misconfigured pool.
        self.ensure_spill_codec_compat()?;

        debug!("Starting raw-bytes sort with order: {:?}", self.sort_order);
        debug!("Memory limit: {} MB", self.memory_limit / (1024 * 1024));
        debug!("Threads: {}", self.threads);

        // Shared worker pool for parallel BGZF compress/decompress across all
        // phases. Sized to the larger of the two phase counts; the active-worker
        // cap (set below, then switched to the Phase-2 count at the merge/write
        // boundary — which may raise or lower it) governs how many run in each
        // phase.
        let pool = Arc::new(SortWorkerPool::new(
            self.phase1_threads().max(self.phase2_threads()),
            self.temp_compression,
            self.output_compression,
            self.spill_codec,
        ));
        // Phase 1 (accumulate/sort/spill) runs on the Phase-1 count; the
        // merge/write phase switches the cap to the Phase-2 count.
        pool.set_active_workers(self.phase1_threads());

        // Open input BAM and create record source
        // N+2 model: workers do ReadInputBlocks + DecompressInput,
        // main thread reads records directly from PooledInputStream.
        debug!(
            "Phase 1: Pool-integrated input reading ({} workers, N+2 model)",
            pool.num_workers()
        );
        let (record_source, header) = {
            let (reader, header) =
                create_raw_bam_reader_pool_integrated(input, &pool, self.async_reader)?;
            (RecordSource::direct(reader), header)
        };

        // Add @PG record if pg_info was provided
        let header = if let Some((ref version, ref command_line)) = self.pg_info {
            fgumi_bam_io::header::add_pg_record(header, version, command_line)?
        } else {
            header
        };

        // _temp_dirs: RAII handles; kept alive until sort returns.
        let (_temp_dirs, mut alloc) = self.create_temp_dirs()?;

        // Sort based on order
        match self.sort_order {
            SortOrder::Coordinate => {
                self.sort_coordinate_optimized(record_source, pool, &header, output, &mut alloc)
            }
            SortOrder::Queryname(comparator) => {
                self.sort_queryname(record_source, pool, &header, output, &mut alloc, comparator)
            }
            SortOrder::TemplateCoordinate => {
                self.sort_template_coordinate(record_source, pool, &header, output, &mut alloc)
            }
        }
    }

    /// Merge multiple pre-sorted BAM files into a single sorted BAM.
    ///
    /// Each input BAM must already be sorted in the order specified by
    /// `self.sort_order`. The output preserves the sort order.
    ///
    /// # Errors
    ///
    /// Returns an error if any input cannot be opened, or writing fails.
    pub fn merge_bams(&self, inputs: &[PathBuf], header: &Header, output: &Path) -> Result<u64> {
        use crate::inline::extract_coordinate_key_inline;
        use crate::keys::{
            QuerynameComparator, RawCoordinateKey, RawQuerynameKey, RawQuerynameLexKey, RawSortKey,
            SortContext,
        };

        debug!("Starting k-way merge of {} BAM files", inputs.len());

        let mut readers = Self::open_bam_prefetch_readers(inputs)?;
        let output_header = self.create_output_header(header);

        match self.sort_order {
            SortOrder::TemplateCoordinate => {
                let lib_lookup = LibraryLookup::from_header(header);
                let cell_tag = self.cell_tag;
                let hasher = cb_hasher();
                self.run_merge_loop(&mut readers, &output_header, output, |bam| {
                    extract_template_key_inline(bam, &lib_lookup, cell_tag, &hasher)
                })
            }
            SortOrder::Coordinate => {
                #[allow(clippy::cast_possible_truncation)]
                let nref = header.reference_sequences().len() as u32;
                self.run_merge_loop(&mut readers, &output_header, output, |bam| RawCoordinateKey {
                    sort_key: extract_coordinate_key_inline(bam, nref),
                })
            }
            SortOrder::Queryname(QuerynameComparator::Lexicographic) => {
                let ctx = SortContext::from_header(header);
                self.run_merge_loop(&mut readers, &output_header, output, |bam| {
                    RawQuerynameLexKey::extract(bam, &ctx)
                })
            }
            SortOrder::Queryname(QuerynameComparator::Natural) => {
                let ctx = SortContext::from_header(header);
                self.run_merge_loop(&mut readers, &output_header, output, |bam| {
                    RawQuerynameKey::extract(bam, &ctx)
                })
            }
        }
    }

    /// Open background prefetch readers for multiple BAM files.
    fn open_bam_prefetch_readers(inputs: &[PathBuf]) -> Result<Vec<RawReadAheadReader>> {
        inputs
            .iter()
            .map(|path| {
                let (reader, _header) = create_raw_bam_reader(path, 1)?;
                Ok(RawReadAheadReader::new(reader))
            })
            .collect()
    }

    /// K-way merge loop: extract keys on the merge thread, write to output.
    ///
    /// Uses reusable per-source record buffers to avoid per-record heap
    /// allocations during the merge.
    fn run_merge_loop<K: Ord>(
        &self,
        readers: &mut [RawReadAheadReader],
        output_header: &Header,
        output: &Path,
        extract_key: impl Fn(&[u8]) -> K,
    ) -> Result<u64> {
        use crate::loser_tree::LoserTree;

        // Initialize: collect first record + key from each reader using
        // reusable per-source buffers (one allocation per source, not per record).
        let mut initial_keys: Vec<K> = Vec::with_capacity(readers.len());
        let mut records: Vec<Vec<u8>> = Vec::with_capacity(readers.len());
        let mut source_map: Vec<usize> = Vec::with_capacity(readers.len());

        for (idx, reader) in readers.iter_mut().enumerate() {
            if let Some(raw_record) = reader.next() {
                let mut buf = Vec::with_capacity(raw_record.as_ref().len());
                buf.extend_from_slice(raw_record.as_ref());
                initial_keys.push(extract_key(&buf));
                records.push(buf);
                source_map.push(idx);
            }
        }

        if initial_keys.is_empty() {
            debug!("Merge complete: 0 records merged");
            let writer = fgumi_bam_io::create_raw_bam_writer(
                output,
                output_header,
                self.threads,
                self.output_compression,
            )?;
            writer.finish()?;
            return Ok(0);
        }

        let mut tree = LoserTree::new(initial_keys);

        let mut writer = fgumi_bam_io::create_raw_bam_writer(
            output,
            output_header,
            self.threads,
            self.output_compression,
        )?;

        let mut records_merged = 0u64;
        let merge_progress = ProgressTracker::new("Merged records").with_interval(1_000_000);

        while tree.winner_is_active() {
            let winner = tree.winner();

            writer.write_raw_record(&records[winner])?;
            records_merged += 1;
            merge_progress.log_if_needed(1);

            let reader_idx = source_map[winner];
            if let Some(raw_record) = readers[reader_idx].next() {
                let buf = &mut records[winner];
                buf.clear();
                buf.extend_from_slice(raw_record.as_ref());
                let new_key = extract_key(buf);
                tree.replace_winner(new_key);
            } else {
                tree.remove_winner();
            }
        }

        writer.finish()?;
        merge_progress.log_final();

        Ok(records_merged)
    }

    /// Optimized coordinate sort using inline buffer for reduced memory overhead.
    ///
    /// Uses `RecordBuffer` which stores records in a single contiguous allocation
    /// with pre-computed sort keys, eliminating per-record heap allocations.
    #[allow(clippy::cast_possible_truncation, clippy::too_many_lines)]
    fn sort_coordinate_optimized(
        &self,
        mut record_source: RecordSource,
        pool: Arc<SortWorkerPool>,
        header: &Header,
        output: &Path,
        alloc: &mut TmpDirAllocator,
    ) -> Result<RawSortStats> {
        use crate::keys::RawCoordinateKey;

        let mut stats = RawSortStats::default();
        let mut timer = SortPhaseTimer::new();

        // Get number of references (unmapped reads map to nref)
        let nref = header.reference_sequences().len() as u32;

        // Estimate capacity from initial_capacity (not memory_limit) to avoid
        // huge upfront allocations when auto-detecting memory.
        let init_cap = self.effective_initial_capacity();
        // Per-record footprint: ~200 bytes BAM + 8 header + 24 ref ≈ 232 bytes (rounded to 240 for headroom)
        let estimated_records = init_cap / 240;
        // Data bytes = init_cap minus ref overhead (24 bytes/record)
        let estimated_data_bytes = init_cap.saturating_sub(estimated_records * 24);

        let mut chunk_files: Vec<PathBuf> = Vec::new();
        let mut buffer = RecordBuffer::with_capacity(estimated_records, estimated_data_bytes, nref);
        let mut namer = ChunkNamer::new(alloc);
        let mut pending_spill: Option<PendingSpill> = None;
        let rayon_pool = self.build_sort_rayon_pool()?;

        let progress = ProgressTracker::new("Read records").with_interval(1_000_000);
        debug!("Phase 1: Reading and sorting chunks (inline buffer, keyed output)...");
        let mut probe = SpillProbe::new("phase1");

        // Borrow each record's bytes straight out of the decompressed block and
        // push them into the arena, skipping the intermediate `RawRecord` copy
        // (the borrowed slice is invalidated by the next `next_record_borrowed`
        // call, which is fine — `push_coordinate` copies the bytes into the buffer).
        while let Some(record) = record_source.next_record_borrowed()? {
            stats.total_records += 1;
            progress.log_if_needed(1);

            // Push directly to buffer - key extracted inline from raw bytes
            buffer.push_coordinate(record)?;

            if probe.should_sample_read(stats.total_records) {
                probe.log_mid_read(probe_stats(&buffer), Some(pool.phase1_queue_depths()));
            }

            // Check memory usage
            if buffer.memory_usage() >= self.memory_limit {
                timer.end_read_span();
                let bstats = probe_stats(&buffer);
                let depths = Some(pool.phase1_queue_depths());
                probe.pre_spill(bstats, depths);

                // Wait for any previous spill to complete before starting a new one
                self.drain_pending_spill::<RawCoordinateKey>(
                    &mut pending_spill,
                    &mut chunk_files,
                    &mut stats,
                    &mut timer,
                    &mut namer,
                    &pool,
                )?;
                probe.post_drain(probe_stats(&buffer), Some(pool.phase1_queue_depths()));

                let chunk_path = namer.next_chunk_path()?;

                timer.time_sort(|| {
                    rayon_pool.install(|| buffer.par_sort());
                });

                // Write keyed temp file with parallel BGZF compression via worker pool.
                // Use start_finish() for pipelining: I/O continues in background
                // while we read the next batch.
                let handle = timer.time_spill_write(|| {
                    let mut writer = PooledChunkWriter::<RawCoordinateKey>::new(
                        Arc::clone(&pool),
                        &chunk_path,
                        pool.spill_codec(),
                    )?;
                    for r in buffer.refs() {
                        let key = RawCoordinateKey { sort_key: r.sort_key };
                        let record_bytes = buffer.get_record(r);
                        writer.write_record(&key, record_bytes)?;
                    }
                    writer.start_finish()
                })?;

                pending_spill = Some(PendingSpill { handle, chunk_path });

                buffer.clear();
                force_mi_collect();
                probe.post_spill(Some(pool.phase1_queue_depths()));
                timer.begin_read_span();
            }
        }

        timer.end_read_span();
        progress.log_final();
        if let Some(err) = record_source.take_error() {
            return Err(anyhow::Error::from(err));
        }

        // Drain any pending spill before merge
        self.drain_pending_spill::<RawCoordinateKey>(
            &mut pending_spill,
            &mut chunk_files,
            &mut stats,
            &mut timer,
            &mut namer,
            &pool,
        )?;
        probe.phase1_end(buffer.memory_usage() as u64);

        // Ingest (accumulate/sort/spill) is done. Switch the worker cap to the
        // Phase-2 count for the writing phase (merge or in-memory write), which
        // is the part that scales and runs after any upstream producer.
        pool.set_active_workers(self.phase2_threads());

        if chunk_files.is_empty() {
            // All records fit in memory - no merge needed
            debug!("All records fit in memory, performing in-memory sort");

            timer.time_sort(|| {
                rayon_pool.install(|| buffer.par_sort());
            });

            timer.time_write_output(|| {
                use crate::pooled_bam_writer::PooledBamWriter;
                let output_header = self.create_output_header(header);
                let mut writer = PooledBamWriter::new(Arc::clone(&pool), output, &output_header)?;

                for record_bytes in buffer.iter_sorted() {
                    writer.write_raw_record(record_bytes)?;
                }
                writer.finish()?;
                Ok(())
            })?;
        } else {
            // Sort remaining records into separate sub-array chunks (avoids
            // intermediate merge back into a single sorted buffer); each
            // chunk becomes its own in-memory merge source.
            let memory_chunks: Vec<InMemoryChunk<RawCoordinateKey>> = if buffer.is_empty() {
                Vec::new()
            } else if self.phase1_threads() > 1 {
                timer.time_sort(|| {
                    rayon_pool.install(|| buffer.par_sort_into_chunks(self.phase1_threads()))
                })
            } else {
                timer.time_sort(|| {
                    rayon_pool.install(|| buffer.par_sort());
                });
                vec![buffer.drain_into_single_chunk()]
            };

            let memory_chunks = MemorySources::Shared(memory_chunks);
            let n_memory = memory_chunks.num_non_empty();
            debug!(
                "Phase 2: Merging {} chunks (keyed O(1) comparisons)...",
                chunk_files.len() + n_memory
            );

            // Merge disk chunks + in-memory chunks using O(1) key comparisons
            timer.time_merge(|| {
                self.merge_chunks_generic::<RawCoordinateKey>(
                    &chunk_files,
                    memory_chunks,
                    header,
                    output,
                    stats.total_records,
                    &pool,
                )
            })?;
        }

        stats.output_records = stats.total_records;
        if let Ok(pool) = Arc::try_unwrap(pool) {
            pool.shutdown();
        }
        timer.log_summary(self.threads);
        debug!("Sort complete: {} records processed", stats.total_records);

        Ok(stats)
    }

    /// Sort by queryname order using keyed temp files for O(1) merge comparisons.
    ///
    /// Dispatches to the appropriate key type based on the comparator:
    /// - `Lexicographic`: uses `RawQuerynameLexKey` (byte comparison)
    /// - `Natural`: uses `RawQuerynameKey` (natural numeric comparison)
    fn sort_queryname(
        &self,
        record_source: RecordSource,
        pool: Arc<SortWorkerPool>,
        header: &Header,
        output: &Path,
        alloc: &mut TmpDirAllocator,
        comparator: QuerynameComparator,
    ) -> Result<RawSortStats> {
        use crate::keys::{RawQuerynameKey, RawQuerynameLexKey};
        debug!("Using queryname sort with {comparator} comparator");
        match comparator {
            QuerynameComparator::Lexicographic => self.sort_queryname_keyed::<RawQuerynameLexKey>(
                record_source,
                pool,
                header,
                output,
                alloc,
            ),
            QuerynameComparator::Natural => self.sort_queryname_keyed::<RawQuerynameKey>(
                record_source,
                pool,
                header,
                output,
                alloc,
            ),
        }
    }

    /// Generic queryname sort using a specific key type.
    #[allow(clippy::too_many_lines)]
    fn sort_queryname_keyed<K: RawSortKey + Default + 'static>(
        &self,
        mut record_source: RecordSource,
        pool: Arc<SortWorkerPool>,
        header: &Header,
        output: &Path,
        alloc: &mut TmpDirAllocator,
    ) -> Result<RawSortStats> {
        use crate::keys::SortContext;

        let mut stats = RawSortStats::default();
        let mut timer = SortPhaseTimer::new();

        let ctx = SortContext::from_header(header);

        // Estimate capacity from initial_capacity to avoid huge upfront allocations.
        let init_cap = self.effective_initial_capacity();
        let estimated_records = init_cap / 300;

        let mut chunk_files: Vec<PathBuf> = Vec::new();
        let mut entries: Vec<(K, fgumi_raw_bam::RawRecord)> = Vec::with_capacity(estimated_records);
        let mut memory_used = 0usize;
        let mut namer = ChunkNamer::new(alloc);
        let mut pending_spill: Option<PendingSpill> = None;
        let rayon_pool = self.build_sort_rayon_pool()?;

        let progress = ProgressTracker::new("Read records").with_interval(1_000_000);
        debug!("Phase 1: Reading and sorting chunks (keyed output)...");
        let mut probe = SpillProbe::new("phase1");

        for record in record_source.by_ref() {
            stats.total_records += 1;
            progress.log_if_needed(1);

            // Extract key from raw bytes
            let key = K::extract(record.as_ref(), &ctx);

            // Estimate memory: record bytes + key overhead
            let record_size = record.as_ref().len() + 50; // approximate key size
            memory_used += record_size;

            entries.push((key, record));

            if probe.should_sample_read(stats.total_records) {
                let bstats = BufferProbeStats::simple(memory_used as u64, entries.len() as u64);
                probe.log_mid_read(bstats, Some(pool.phase1_queue_depths()));
            }

            // Check if we need to spill to disk
            if memory_used >= self.memory_limit {
                timer.end_read_span();
                let bstats = BufferProbeStats::simple(memory_used as u64, entries.len() as u64);
                let depths = Some(pool.phase1_queue_depths());
                probe.pre_spill(bstats, depths);

                // Wait for any previous spill to complete
                self.drain_pending_spill::<K>(
                    &mut pending_spill,
                    &mut chunk_files,
                    &mut stats,
                    &mut timer,
                    &mut namer,
                    &pool,
                )?;
                probe.post_drain(bstats, Some(pool.phase1_queue_depths()));

                let chunk_path = namer.next_chunk_path()?;

                timer.time_sort(|| {
                    use rayon::prelude::*;
                    rayon_pool.install(|| entries.par_sort_unstable_by(|a, b| a.0.cmp(&b.0)));
                });

                // Write keyed temp file with parallel BGZF compression via worker pool.
                let handle = timer.time_spill_write(|| {
                    let mut writer = PooledChunkWriter::<K>::new(
                        Arc::clone(&pool),
                        &chunk_path,
                        pool.spill_codec(),
                    )?;
                    for (key, record) in entries.drain(..) {
                        writer.write_record(&key, record.as_ref())?;
                    }
                    writer.start_finish()
                })?;

                pending_spill = Some(PendingSpill { handle, chunk_path });

                memory_used = 0;
                force_mi_collect();
                probe.post_spill(Some(pool.phase1_queue_depths()));
                timer.begin_read_span();
            }
        }

        timer.end_read_span();
        progress.log_final();
        if let Some(err) = record_source.take_error() {
            return Err(anyhow::Error::from(err));
        }

        // Drain any pending spill before merge
        self.drain_pending_spill::<K>(
            &mut pending_spill,
            &mut chunk_files,
            &mut stats,
            &mut timer,
            &mut namer,
            &pool,
        )?;
        probe.phase1_end(memory_used as u64);

        // Ingest done; switch the worker cap to the Phase-2 count for the write phase.
        pool.set_active_workers(self.phase2_threads());

        if chunk_files.is_empty() {
            // All records fit in memory
            debug!("All records fit in memory, performing in-memory sort");

            timer.time_sort(|| {
                use rayon::prelude::*;
                rayon_pool.install(|| entries.par_sort_unstable_by(|a, b| a.0.cmp(&b.0)));
            });

            timer.time_write_output(|| {
                use crate::pooled_bam_writer::PooledBamWriter;
                let output_header = self.create_output_header(header);
                let mut writer = PooledBamWriter::new(Arc::clone(&pool), output, &output_header)?;

                for (_key, record) in entries {
                    writer.write_raw_record(&record)?;
                }
                writer.finish()?;
                Ok(())
            })?;
        } else {
            // Sort remaining records into separate sub-array chunks (avoids
            // intermediate merge back into a single sorted buffer). The
            // queryname path accumulates records as individual `RawRecord`
            // allocations upstream (per-record alloc already paid; this is the
            // documented S3-015 tradeoff — see
            // docs/design/sort-queryname-arena-deferral.md), so we sort the
            // keyed `Vec`s in place and pack the owned `RawRecord`s into
            // `MemorySources::Owned` for the merge interface — unlike the
            // inline-buffer paths (coordinate, template), which produce
            // `InMemoryChunk`s sharing an `Arc<SegmentedBuf>`.
            let keyed_chunks: Vec<Vec<(K, fgumi_raw_bam::RawRecord)>> = if entries.is_empty() {
                Vec::new()
            } else if self.phase1_threads() > 1 {
                timer.time_sort(|| {
                    use rayon::prelude::*;
                    let chunk_size = entries.len().div_ceil(self.phase1_threads());
                    rayon_pool.install(|| {
                        entries.par_chunks_mut(chunk_size).for_each(|chunk| {
                            chunk.sort_unstable_by(|a, b| a.0.cmp(&b.0));
                        });
                    });
                    // Carve sub-chunks aligned with par_chunks_mut boundaries.
                    // par_chunks_mut produces [0..cs), [cs..2cs), ..., [n-tail..n).
                    // We peel the short tail first, then full chunks from the end,
                    // and reverse. Each split_off is O(chunk_size) → O(n) total.
                    let mut remaining = std::mem::take(&mut entries);
                    let num_chunks = remaining.len().div_ceil(chunk_size);
                    let mut chunks: Vec<Vec<(K, fgumi_raw_bam::RawRecord)>> =
                        Vec::with_capacity(num_chunks);
                    let tail_len = remaining.len() % chunk_size;
                    if tail_len != 0 {
                        let split_at = remaining.len() - tail_len;
                        chunks.push(remaining.split_off(split_at));
                    }
                    while !remaining.is_empty() {
                        let split_at = remaining.len().saturating_sub(chunk_size);
                        chunks.push(remaining.split_off(split_at));
                    }
                    chunks.reverse();
                    chunks
                })
            } else {
                timer.time_sort(|| {
                    entries.sort_unstable_by(|a, b| a.0.cmp(&b.0));
                });
                vec![entries]
            };
            let memory_chunks = MemorySources::Owned(keyed_chunks);

            let n_memory = memory_chunks.num_non_empty();
            debug!(
                "Phase 2: Merging {} chunks (keyed comparisons)...",
                chunk_files.len() + n_memory
            );

            // Merge disk chunks + in-memory records using keyed comparisons
            timer.time_merge(|| {
                self.merge_chunks_generic::<K>(
                    &chunk_files,
                    memory_chunks,
                    header,
                    output,
                    stats.total_records,
                    &pool,
                )
            })?;
        }

        stats.output_records = stats.total_records;
        if let Ok(pool) = Arc::try_unwrap(pool) {
            pool.shutdown();
        }
        timer.log_summary(self.threads);
        debug!("Sort complete: {} records processed", stats.total_records);

        Ok(stats)
    }

    /// Sort by template-coordinate order using inline buffer for reduced memory.
    ///
    /// Uses `TemplateRecordBuffer` which stores records in a single contiguous allocation
    /// with packed sort keys, eliminating per-record heap allocations for names.
    ///
    /// Writes keyed temp chunks that preserve pre-computed sort keys, enabling O(1)
    /// comparisons during merge (instead of expensive CIGAR/aux parsing).
    fn sort_template_coordinate(
        &self,
        mut record_source: RecordSource,
        pool: Arc<SortWorkerPool>,
        header: &Header,
        output: &Path,
        alloc: &mut TmpDirAllocator,
    ) -> Result<RawSortStats> {
        let lib_lookup = LibraryLookup::from_header(header);
        let cb_hasher = cb_hasher();

        // Peek the first record to provision the variant (extraction stays
        // full-width). The first record is handed (by reference) to the impl,
        // which processes it before draining the rest of `record_source`; the
        // wrapper only extracts `first_key` to choose `K` — it does not re-chain
        // the record into the iterator.
        let first_record = record_source.by_ref().next();
        let first_key = first_record.as_ref().map(|r| {
            extract_template_key_inline(r.as_ref(), &lib_lookup, self.cell_tag, &cb_hasher)
        });
        let header_library_varies = lib_lookup.distinct_header_ordinals() > 1;
        let variant =
            select_template_variant(first_key.as_ref(), self.key_types, header_library_varies);

        debug!(
            "template-coordinate variant: {} lanes (cb={}, tertiary={})",
            variant.lanes(),
            variant.cb,
            variant.tertiary
        );

        match (variant.cb, variant.tertiary) {
            (false, false) => self.sort_template_coordinate_impl::<TemplateKey24>(
                record_source,
                first_record.as_ref(),
                first_key,
                variant,
                pool,
                header,
                output,
                alloc,
                &lib_lookup,
                &cb_hasher,
            ),
            (true, false) => self.sort_template_coordinate_impl::<CbKey32>(
                record_source,
                first_record.as_ref(),
                first_key,
                variant,
                pool,
                header,
                output,
                alloc,
                &lib_lookup,
                &cb_hasher,
            ),
            (false, true) => self.sort_template_coordinate_impl::<TertKey32>(
                record_source,
                first_record.as_ref(),
                first_key,
                variant,
                pool,
                header,
                output,
                alloc,
                &lib_lookup,
                &cb_hasher,
            ),
            (true, true) => self.sort_template_coordinate_impl::<TemplateKey40>(
                record_source,
                first_record.as_ref(),
                first_key,
                variant,
                pool,
                header,
                output,
                alloc,
                &lib_lookup,
                &cb_hasher,
            ),
        }
    }

    /// Generic template-coordinate sort body, monomorphized over the chosen
    /// lane-key type `K`.
    ///
    /// The dispatch wrapper [`Self::sort_template_coordinate`] peeks the first
    /// record to choose `K`, then hands the captured `first_record` (already
    /// pulled out of `record_source`, borrowed here) plus the rest of the
    /// iterator. This impl processes `first_record` in a pre-loop block and
    /// drains the remainder in the read loop. `first_key` is the full key of the
    /// first record; the decode-time verify compares every subsequent record's
    /// dropped lanes against it.
    #[allow(clippy::too_many_lines)]
    // The generic body keeps the same parameter shape as the legacy function
    // plus the peeked-first-record handoff; the existing code already allows
    // this pattern for the merge helpers.
    #[allow(clippy::too_many_arguments)]
    fn sort_template_coordinate_impl<K: TemplateLaneKey + RawSortKey + Default + 'static>(
        &self,
        mut record_source: RecordSource,
        first_record: Option<&fgumi_raw_bam::RawRecord>,
        first_key: Option<TemplateKey>,
        variant: TemplateKeyVariant,
        pool: Arc<SortWorkerPool>,
        header: &Header,
        output: &Path,
        alloc: &mut TmpDirAllocator,
        lib_lookup: &LibraryLookup,
        cb_hasher: &ahash::RandomState,
    ) -> Result<RawSortStats> {
        let mut stats = RawSortStats::default();
        let mut timer = SortPhaseTimer::new();

        // The full key the decode-time verify compares dropped lanes against.
        // On empty input, `first_record` is `None` so the pre-loop push is
        // skipped and the read loop runs zero times — `first` is never read in
        // that case, so a default (all-zero) key is a harmless placeholder that
        // lets the function fall through to the post-loop empty-output path
        // (header-only BAM), matching the coordinate and queryname orders.
        let first = first_key.unwrap_or_default();

        // Per-record arena footprint (data side: inline header + BAM bytes; ref side:
        // the cached key + offset/len/pad). Ref size tracks the chosen key width K:
        //   K=TemplateKey24 → ref=40 B, bpr=298 B, data%≈86.6%
        //   K=TemplateKey32 → ref=48 B, bpr=306 B, data%≈84.3%
        //   K=TemplateKey40 → ref=56 B, bpr=314 B, data%≈82.2%
        let ref_bytes = std::mem::size_of::<TemplateRecordRef<K>>();
        let data_bytes_per_record = TEMPLATE_HEADER_SIZE + EST_BAM_BYTES_PER_TEMPLATE_RECORD;
        let bytes_per_record = data_bytes_per_record + ref_bytes;

        let init_cap = self.effective_initial_capacity();
        let estimated_records = (init_cap / bytes_per_record).max(1);
        let estimated_data_bytes = init_cap * data_bytes_per_record / bytes_per_record;

        let mut chunk_files: Vec<PathBuf> = Vec::new();
        let mut buffer =
            TemplateRecordBuffer::<K>::with_capacity(estimated_records, estimated_data_bytes);
        let mut namer = ChunkNamer::new(alloc);
        let mut pending_spill: Option<PendingSpill> = None;
        let rayon_pool = self.build_sort_rayon_pool()?;

        let progress = ProgressTracker::new("Read records").with_interval(1_000_000);
        debug!("Phase 1: Reading and sorting chunks (inline buffer)...");
        let mut probe = SpillProbe::new("phase1");

        // Process the captured first record before draining the rest: extract its
        // full key, verify the dropped lanes match `first` (trivially true for the
        // first record itself), and push the narrowed key. A single record cannot
        // exceed the memory limit, so no spill check is needed here.
        if let Some(record) = first_record {
            stats.total_records += 1;
            progress.log_if_needed(1);

            let bam_bytes = record.as_ref();
            let full = extract_template_key_inline(bam_bytes, lib_lookup, self.cell_tag, cb_hasher);
            if let Some(violation) = verify_dropped_lanes(&first, &full, variant) {
                let name = String::from_utf8_lossy(
                    fgumi_raw_bam::RawRecordView::new(bam_bytes).read_name(),
                )
                .into_owned();
                return Err(dropped_lane_error(&name, violation));
            }
            buffer.push(bam_bytes, K::from_full(&full))?;
        }

        // Borrow each record's bytes in place (see the coordinate ingest loop);
        // the key is extracted and the bytes copied into the buffer before the
        // borrow ends, so no owned `RawRecord` is needed here.
        while let Some(bam_bytes) = record_source.next_record_borrowed()? {
            stats.total_records += 1;
            progress.log_if_needed(1);

            // Extract the full template key, verify the lanes the chosen variant
            // dropped are constant relative to the first record, then push the
            // narrowed key.
            let full = extract_template_key_inline(bam_bytes, lib_lookup, self.cell_tag, cb_hasher);
            if let Some(violation) = verify_dropped_lanes(&first, &full, variant) {
                let name = String::from_utf8_lossy(
                    fgumi_raw_bam::RawRecordView::new(bam_bytes).read_name(),
                )
                .into_owned();
                return Err(dropped_lane_error(&name, violation));
            }
            buffer.push(bam_bytes, K::from_full(&full))?;

            if probe.should_sample_read(stats.total_records) {
                probe.log_mid_read(probe_stats(&buffer), Some(pool.phase1_queue_depths()));
            }

            // Check memory usage
            if buffer.memory_usage() >= self.memory_limit {
                timer.end_read_span();
                let bstats = probe_stats(&buffer);
                let depths = Some(pool.phase1_queue_depths());
                probe.pre_spill(bstats, depths);

                // Wait for any previous spill to complete
                self.drain_pending_spill::<K>(
                    &mut pending_spill,
                    &mut chunk_files,
                    &mut stats,
                    &mut timer,
                    &mut namer,
                    &pool,
                )?;
                probe.post_drain(probe_stats(&buffer), Some(pool.phase1_queue_depths()));

                let chunk_path = namer.next_chunk_path()?;

                timer.time_sort(|| {
                    rayon_pool.install(|| buffer.par_sort());
                });

                // Write keyed chunk with parallel BGZF compression via worker pool.
                let handle = timer.time_spill_write(|| {
                    let mut writer = PooledChunkWriter::<K>::new(
                        Arc::clone(&pool),
                        &chunk_path,
                        pool.spill_codec(),
                    )?;
                    for (key, record) in buffer.iter_sorted_keyed() {
                        writer.write_record(&key, record)?;
                    }
                    writer.start_finish()
                })?;

                pending_spill = Some(PendingSpill { handle, chunk_path });

                buffer.clear();
                force_mi_collect();
                probe.post_spill(Some(pool.phase1_queue_depths()));
                timer.begin_read_span();
            }
        }

        timer.end_read_span();
        progress.log_final();
        if let Some(err) = record_source.take_error() {
            return Err(anyhow::Error::from(err));
        }

        // Drain any pending spill before merge
        self.drain_pending_spill::<K>(
            &mut pending_spill,
            &mut chunk_files,
            &mut stats,
            &mut timer,
            &mut namer,
            &pool,
        )?;
        probe.phase1_end(buffer.memory_usage() as u64);

        // Ingest done; switch the worker cap to the Phase-2 count for the write phase.
        pool.set_active_workers(self.phase2_threads());

        if chunk_files.is_empty() {
            // All records fit in memory
            debug!("All records fit in memory, performing in-memory sort");

            timer.time_sort(|| {
                rayon_pool.install(|| buffer.par_sort());
            });

            timer.time_write_output(|| {
                use crate::pooled_bam_writer::PooledBamWriter;
                let output_header = self.create_output_header(header);
                let mut writer = PooledBamWriter::new(Arc::clone(&pool), output, &output_header)?;

                for record_bytes in buffer.iter_sorted() {
                    writer.write_raw_record(record_bytes)?;
                }
                writer.finish()?;
                Ok(())
            })?;
        } else {
            // Sort remaining records into separate sub-array chunks (avoids
            // intermediate merge back into a single sorted buffer)
            let memory_chunks: Vec<InMemoryChunk<K>> = if buffer.is_empty() {
                Vec::new()
            } else if self.phase1_threads() > 1 {
                timer.time_sort(|| {
                    rayon_pool.install(|| buffer.par_sort_into_chunks(self.phase1_threads()))
                })
            } else {
                timer.time_sort(|| {
                    rayon_pool.install(|| buffer.par_sort());
                });
                vec![buffer.drain_into_single_chunk()]
            };

            let memory_chunks = MemorySources::Shared(memory_chunks);
            let n_memory = memory_chunks.num_non_empty();
            debug!("Phase 2: Merging {} chunks...", chunk_files.len() + n_memory);

            // Merge using O(1) key comparisons
            timer.time_merge(|| {
                self.merge_chunks_generic::<K>(
                    &chunk_files,
                    memory_chunks,
                    header,
                    output,
                    stats.total_records,
                    &pool,
                )
            })?;
        }

        stats.output_records = stats.total_records;
        if let Ok(pool) = Arc::try_unwrap(pool) {
            pool.shutdown();
        }
        timer.log_summary(self.threads);
        debug!("Sort complete: {} records processed", stats.total_records);

        Ok(stats)
    }

    /// Build chunk sources from disk files and in-memory chunks.
    ///
    /// Disk files become `PoolDisk` sources: pool workers read and decompress
    /// them via work-stealing while the main thread parses records, so no
    /// per-source threads are spawned. In-memory chunks are appended as-is.
    /// (The index/verify path does not call this — it builds
    /// `GenericKeyedChunkReader` directly via `run_merge_loop`.)
    fn build_chunk_sources<K: RawSortKey + Default + 'static>(
        chunk_files: &[PathBuf],
        memory_chunks: MemorySources<K>,
        pool: &Arc<SortWorkerPool>,
    ) -> Result<Vec<ChunkSource<K>>> {
        let num_disk = chunk_files.len();
        let num_memory = memory_chunks.num_non_empty();
        let mut sources: Vec<ChunkSource<K>> = Vec::with_capacity(num_disk + num_memory);

        if !chunk_files.is_empty() {
            // Install per-file Phase 2 state on the pool, then create one
            // PoolDisk source per file. Workers cooperatively read+decompress
            // all files via work-stealing. (The legacy per-source-reader-thread
            // `ChunkSource::Disk` branch was removed in S3-003 slimming: it was
            // never reachable — `merge_chunks_generic` is the only caller and
            // always used the pool path; the index/verify path builds
            // `GenericKeyedChunkReader` directly via `run_merge_loop`.)
            pool.set_phase2_files(chunk_files)?;

            for source_id in 0..num_disk {
                sources.push(ChunkSource::PoolDisk { source_id, scratch: Vec::new() });
            }
        }

        memory_chunks.push_into(&mut sources);

        Ok(sources)
    }

    /// Generic merge for keyed chunks using `O(1)` key comparisons.
    ///
    /// This is the unified merge function that works with any `RawSortKey` type.
    /// It provides `O(1)` comparisons during merge for fixed-size keys (coordinate, template)
    /// and `O(name_len)` for variable-size keys (queryname).
    #[allow(clippy::too_many_lines)]
    fn merge_chunks_generic<K: RawSortKey + Default + 'static>(
        &self,
        chunk_files: &[PathBuf],
        memory_chunks: MemorySources<K>,
        header: &Header,
        output: &Path,
        total_records: u64,
        pool: &Arc<SortWorkerPool>,
    ) -> Result<u64> {
        use crate::loser_tree::LoserTree;
        use crate::pooled_bam_writer::PooledBamWriter;
        use crate::worker_pool::phase;

        let num_disk = chunk_files.len();

        if num_disk > 0 {
            debug!(
                "Pool-integrated merge: {} disk sources, {} pool workers (N+2 model)",
                num_disk,
                pool.num_workers()
            );
        }

        let mut sources = Self::build_chunk_sources::<K>(chunk_files, memory_chunks, pool)?;

        let num_sources = sources.len();
        debug!("Merging from {num_sources} sources...");

        // Create pool consumer for PoolDisk sources and activate Phase 2.
        // The consumer holds an Arc snapshot of the pool's per-file Phase 2
        // state — workers populate per-file reorder buffers, the consumer pops
        // from them.
        let mut guard: Phase2Guard<'_, K> = if num_disk > 0 {
            let files = pool.phase2_files();
            let consumer = MainThreadChunkConsumer::new(
                files,
                pool.decompress_error_flag(),
                pool.chunk_read_error_flag(),
                pool.worker_panicked_flag(),
            );
            pool.set_phase(phase::PHASE2);
            Phase2Guard { pool, consumer: Some(consumer), active: true }
        } else {
            Phase2Guard { pool, consumer: None, active: false }
        };

        let output_header = self.create_output_header(header);

        // Initialize loser tree with first record from each source.
        // Each source owns its current-record state internally; no
        // intermediate `records: Vec<Vec<u8>>` is needed.
        let mut initial_keys: Vec<K> = Vec::with_capacity(sources.len());
        let mut source_map: Vec<usize> = Vec::with_capacity(sources.len());

        for (idx, source) in sources.iter_mut().enumerate() {
            if let Some(key) = source.advance(guard.consumer_mut())? {
                initial_keys.push(key);
                source_map.push(idx);
            }
        }

        if initial_keys.is_empty() {
            debug!("Merge complete: 0 records merged");
            guard.deactivate();
            let writer = PooledBamWriter::new(Arc::clone(pool), output, &output_header)?;
            writer.finish()?;
            return Ok(0);
        }

        let mut tree = LoserTree::new(initial_keys);

        debug!("Merge thread budget: {} pool workers + 1 I/O + 1 main (N+2)", pool.num_workers());
        let mut writer = PooledBamWriter::new(Arc::clone(pool), output, &output_header)?;

        let mut records_merged = 0u64;
        let merge_progress = ProgressTracker::new("Merged records")
            .with_interval(1_000_000)
            .with_total(total_records);

        let mut merge_probe = MergeProbe::new();

        // Sub-phase timing: only paid when debug logging is enabled.
        let debug_timing = log::log_enabled!(log::Level::Debug);
        let merge_sample_interval: u64 = 1024;
        let mut merge_write_secs = 0.0f64;
        let mut merge_read_secs = 0.0f64;
        let mut merge_tree_secs = 0.0f64;
        let mut samples_taken: u64 = 0;
        let loop_start = Instant::now();

        while tree.winner_is_active() {
            let winner = tree.winner();
            let src_idx = source_map[winner];
            let sample_this = debug_timing && records_merged.is_multiple_of(merge_sample_interval);

            if sample_this {
                let t0 = Instant::now();
                writer.write_raw_record(sources[src_idx].current_bytes())?;
                merge_write_secs += t0.elapsed().as_secs_f64();
            } else {
                writer.write_raw_record(sources[src_idx].current_bytes())?;
            }

            records_merged += 1;
            merge_progress.log_if_needed(1);

            if merge_probe.should_sample(records_merged) {
                let depths = pool.phase1_queue_depths();
                let consumer_stats = guard.consumer_mut().map(|c| c.probe_consumer_stats());
                merge_probe.log_mid_with_depths(depths, consumer_stats);
            }

            if sample_this {
                let t0 = Instant::now();
                let next = sources[src_idx].advance(guard.consumer_mut())?;
                merge_read_secs += t0.elapsed().as_secs_f64();

                let t0 = Instant::now();
                if let Some(key) = next {
                    tree.replace_winner(key);
                } else {
                    tree.remove_winner();
                }
                merge_tree_secs += t0.elapsed().as_secs_f64();
                samples_taken += 1;
            } else {
                let next = sources[src_idx].advance(guard.consumer_mut())?;
                if let Some(key) = next {
                    tree.replace_winner(key);
                } else {
                    tree.remove_winner();
                }
            }
        }

        // Return to legacy mode before finishing writer (workers still need to compress).
        // The guard's deactivate() drops the consumer, resets phase, and clears the
        // pool's published file vector in the correct order.
        guard.deactivate();

        if debug_timing {
            let loop_total = loop_start.elapsed().as_secs_f64();
            #[allow(clippy::cast_precision_loss)]
            let scale =
                if samples_taken > 0 { records_merged as f64 / samples_taken as f64 } else { 1.0 };
            let est_write = merge_write_secs * scale;
            let est_read = merge_read_secs * scale;
            let est_tree = merge_tree_secs * scale;
            debug!(
                "Merge sub-phases (sampled {samples_taken}/{records_merged}, scale={scale:.1}x): write={est_write:.2}s read={est_read:.2}s tree={est_tree:.2}s total={loop_total:.2}s records={records_merged}"
            );
        }

        writer.finish()?;
        merge_progress.log_final();
        log_snapshot("phase2.end", 0);

        Ok(records_merged)
    }

    /// Create output header with appropriate sort order tags.
    fn create_output_header(&self, header: &Header) -> Header {
        super::create_output_header(self.sort_order, header)
    }

    /// Create the spill temp dirs + allocator for the P6 `CompressSpill` step.
    ///
    /// Public wrapper over [`Self::create_temp_dirs`]: honors `--temp-dir` /
    /// `FGUMI_TMP_DIRS` (the configured `temp_dirs`) exactly as the legacy
    /// streaming path does. The caller hands the [`Vec<TempDir>`] RAII handles to
    /// `CompressSpill` (held for the step's lifetime) and the [`TmpDirAllocator`]
    /// behind a shared mutex for per-file base-dir allocation.
    ///
    /// # Errors
    ///
    /// Returns an error if a temp directory cannot be created or has insufficient
    /// free space.
    pub fn create_spill_dirs(&self) -> Result<(Vec<TempDir>, TmpDirAllocator)> {
        self.create_temp_dirs()
    }

    /// Create per-base temp directories and an allocator over their subdirs.
    ///
    /// For each user-supplied base directory, a fresh sort-run subdirectory is
    /// created (via `tempfile::TempDir`). The returned [`Vec<TempDir>`] owns
    /// those handles so subdirs are removed on drop; the allocator hands out
    /// the corresponding subdir paths for chunk/merged file placement.
    ///
    /// When `temp_dirs` is empty, a single subdirectory is created under the
    /// system default temp location.
    fn create_temp_dirs(&self) -> Result<(Vec<TempDir>, TmpDirAllocator)> {
        use super::create_temp_dir;

        if self.temp_dirs.is_empty() {
            let td = create_temp_dir(None)?;
            let base = td.path().to_path_buf();
            let alloc = TmpDirAllocator::new(vec![base])?;
            return Ok((vec![td], alloc));
        }

        let mut handles = Vec::with_capacity(self.temp_dirs.len());
        let mut subdirs = Vec::with_capacity(self.temp_dirs.len());
        for base in &self.temp_dirs {
            let td = create_temp_dir(Some(base))?;
            subdirs.push(td.path().to_path_buf());
            handles.push(td);
        }
        let alloc = TmpDirAllocator::new(subdirs)?;
        Ok((handles, alloc))
    }
}

/// Extract a packed `TemplateKey` directly from BAM record bytes.
///
/// This function computes the template-coordinate sort key inline, avoiding
/// heap allocations for the read name by using a hash instead.
///
/// When `cell_tag` is `Some`, the CB (cellular barcode) tag value is hashed
/// and included in the sort key between neg2 and MI, matching fgbio's order.
#[must_use]
pub fn extract_template_key_inline(
    bam_bytes: &[u8],
    lib_lookup: &LibraryLookup,
    cell_tag: Option<SamTag>,
    cb_hasher: &ahash::RandomState,
) -> TemplateKey {
    use fgumi_raw_bam;
    use fgumi_raw_bam::{flags, mate_unclipped_5prime, unclipped_5prime_raw};

    // Single-pass extraction of all aux tags (MI, RG, cell barcode, MC)
    let aux = fgumi_raw_bam::extract_template_aux_tags(bam_bytes, cell_tag);
    let mi = aux.mi;
    let library = lib_lookup.ordinal_from_rg(aux.rg);
    let cb_hash = aux.cell.map_or(0u64, |cb_bytes| cb_hasher.hash_one(cb_bytes));

    // Extract fields from raw bytes
    let v = fgumi_raw_bam::RawRecordView::new(bam_bytes);
    let tid = v.ref_id();
    let pos = v.pos();
    let l_read_name = v.l_read_name() as usize;
    let flag = v.flags();
    let mate_tid = v.mate_ref_id();
    let mate_pos = v.mate_pos();

    // Extract flags
    let is_unmapped = (flag & flags::UNMAPPED) != 0;
    let mate_unmapped = (flag & flags::MATE_UNMAPPED) != 0;
    let is_reverse = (flag & flags::REVERSE) != 0;
    let mate_reverse = (flag & flags::MATE_REVERSE) != 0;
    let is_paired = (flag & flags::PAIRED) != 0;

    // Hash read name (exclude null terminator)
    let name_len = l_read_name.saturating_sub(1);
    let name = if name_len > 0 && 32 + name_len <= bam_bytes.len() {
        &bam_bytes[32..32 + name_len]
    } else {
        &[]
    };
    let name_hash = lib_lookup.hash_name(name);

    // Handle unmapped reads
    if is_unmapped {
        if is_paired && !mate_unmapped {
            // Unmapped read with mapped mate - use mate's position as primary key
            let mate_unclipped =
                aux.mc.map_or(mate_pos, |mc| mate_unclipped_5prime(mate_pos, mate_reverse, mc));

            return TemplateKey::new(
                mate_tid,
                mate_unclipped,
                mate_reverse,
                i32::MAX,
                i32::MAX,
                false,
                cb_hash,
                library,
                mi,
                name_hash,
                true, // Unmapped read is always "upper" relative to mapped mate
            );
        }

        // Completely unmapped - sort to end. Still carry the read's library and
        // MI lanes so a fully-unmapped read realizes the same tertiary as its
        // mapped, same-library peers (otherwise the dropped-lane verify treats a
        // single-library file as if the library varied; see #375).
        let is_read2 = (flag & 0x80) != 0; // is_last_segment flag
        return TemplateKey::unmapped(name_hash, cb_hash, library, mi, is_read2);
    }

    // Calculate unclipped 5' position for this read (zero-alloc: reads cigar directly)
    let this_pos = unclipped_5prime_raw(bam_bytes, pos, is_reverse);

    // Calculate mate's unclipped 5' position
    let mate_unclipped = if is_paired && !mate_unmapped {
        aux.mc.map_or(mate_pos, |mc| mate_unclipped_5prime(mate_pos, mate_reverse, mc))
    } else {
        mate_pos
    };

    // Determine canonical ordering
    let (tid1, tid2, pos1, pos2, neg1, neg2, is_upper) = if is_paired && !mate_unmapped {
        // Samtools logic: is_upper if pos > mate_pos, or (pos == mate_pos && this read is reverse)
        let is_upper = (tid, this_pos) > (mate_tid, mate_unclipped)
            || ((tid, this_pos) == (mate_tid, mate_unclipped) && is_reverse);

        if is_upper {
            // Swap: mate's position comes first
            (mate_tid, tid, mate_unclipped, this_pos, mate_reverse, is_reverse, true)
        } else {
            // No swap: this read's position comes first
            (tid, mate_tid, this_pos, mate_unclipped, is_reverse, mate_reverse, false)
        }
    } else {
        // Unpaired or mate unmapped - use MAX for tid2/pos2
        (tid, i32::MAX, this_pos, i32::MAX, is_reverse, false, false)
    };

    TemplateKey::new(tid1, pos1, neg1, tid2, pos2, neg2, cb_hash, library, mi, name_hash, is_upper)
}

// SortStats is defined in crate root (lib.rs); alias it here for internal use.
pub(crate) use crate::SortStats as RawSortStats;

// ============================================================================
// Pool-integrated BAM reader construction
//
// This function depends on sort-internal types (`SortWorkerPool`,
// `PooledInputStream`, `phase`), which is why it lives here in `fgumi-sort`
// rather than in `fgumi-bam-io`.
// ============================================================================

/// Create a raw BAM reader using the pool's Phase 1 integrated reading.
///
/// Workers in the pool do `ReadInputBlocks` + `DecompressInput`. The main
/// thread consumes decompressed bytes via `PooledInputStream`.
///
/// When `async_reader` is false, no extra threads are spawned: the pool's
/// block reader reads directly from the input file. When `async_reader` is
/// true, the input file is wrapped in a `PrefetchReader`, which spawns one
/// dedicated OS thread (`fgumi-prefetch`) that reads raw bytes ahead into a
/// bounded queue so the pool's block reader never blocks on disk I/O.
///
/// # Flow
///
/// 1. Parse header with noodles (first pass, single-threaded)
/// 2. Open file again, set as pool's input file
/// 3. Set pool to PHASE1 — workers start reading/decompressing
/// 4. Main thread skips header bytes from decompressed stream
/// 5. Returns `RawBamReader<PooledInputStream>` for direct record iteration
///
/// # Errors
///
/// Returns an error if the BAM file cannot be opened, the header cannot be
/// parsed, or header bytes cannot be skipped from the decompressed stream.
fn create_raw_bam_reader_pool_integrated<P: AsRef<Path>>(
    path: P,
    pool: &Arc<SortWorkerPool>,
    async_reader: bool,
) -> Result<(fgumi_raw_bam::RawBamReader<PooledInputStream>, Header)> {
    use crate::worker_pool::phase;
    use std::io;

    let path_ref = path.as_ref();

    let (header, reader): (Header, Box<dyn io::Read + Send>) = if is_stdin_path(path_ref) {
        let stdin = io::stdin();
        let tee = TeeReader::new(stdin);
        let bgzf = BgzfReader::new(tee);
        let mut noodles_reader = noodles::bam::io::Reader::from(bgzf);
        let header =
            noodles_reader.read_header().with_context(|| "Failed to read BAM header from stdin")?;

        let bgzf = noodles_reader.into_inner();
        let tee = bgzf.into_inner();
        let (buffered_bytes, stdin) = tee.into_parts();
        let chained = ChainedReader::new(buffered_bytes, stdin);
        (header, Box::new(io::BufReader::with_capacity(2 * 1024 * 1024, chained)))
    } else {
        let mut file = std::fs::File::open(path_ref)
            .with_context(|| format!("Failed to open input BAM: {}", path_ref.display()))?;

        let header = {
            let bgzf = BgzfReader::new(&mut file);
            let mut noodles_reader = noodles::bam::io::Reader::from(bgzf);
            noodles_reader
                .read_header()
                .with_context(|| format!("Failed to read header from: {}", path_ref.display()))?
        };

        file.seek(SeekFrom::Start(0))
            .with_context(|| format!("Failed to rewind input BAM: {}", path_ref.display()))?;

        let reader: Box<dyn io::Read + Send> = if async_reader {
            fgumi_bam_io::os_hints::advise_sequential(&file);
            log::debug!(
                "async sort reader enabled: spawning fgumi-prefetch thread for {}",
                path_ref.display()
            );
            Box::new(fgumi_bam_io::prefetch_reader::PrefetchReader::from_file(file))
        } else {
            Box::new(io::BufReader::with_capacity(2 * 1024 * 1024, file))
        };
        (header, reader)
    };

    pool.set_input_file(reader);
    pool.set_phase(phase::PHASE1);

    let mut pooled_input = PooledInputStream::new(
        pool.decompressed_input_queue(),
        pool.decompressed_input_done_flag(),
        pool.input_read_error_flag(),
        pool.decompress_error_flag(),
    );

    skip_bam_header(&mut pooled_input)
        .with_context(|| format!("Failed to skip header from: {}", path_ref.display()))?;

    let raw_reader = fgumi_raw_bam::RawBamReader::new(pooled_input);
    Ok((raw_reader, header))
}

/// Skip the BAM header from a reader positioned at the start of a BAM stream.
///
/// Reads and discards: magic (4 bytes), header text length + text,
/// `n_ref` + reference entries.
fn skip_bam_header<R: Read>(reader: &mut R) -> Result<()> {
    use std::io;
    let mut buf4 = [0u8; 4];

    reader.read_exact(&mut buf4)?;
    anyhow::ensure!(&buf4 == b"BAM\x01", "Not a BAM file (bad magic)");

    reader.read_exact(&mut buf4)?;
    let l_text = u32::from_le_bytes(buf4) as usize;
    let copied = io::copy(&mut reader.take(l_text as u64), &mut io::sink())?;
    anyhow::ensure!(
        copied == l_text as u64,
        "BAM header text truncated: expected {l_text} bytes, got {copied}"
    );

    reader.read_exact(&mut buf4)?;
    let n_ref = u32::from_le_bytes(buf4) as usize;

    for _ in 0..n_ref {
        reader.read_exact(&mut buf4)?;
        let l_name = u32::from_le_bytes(buf4) as usize;
        let copied = io::copy(&mut reader.take(l_name as u64), &mut io::sink())?;
        anyhow::ensure!(
            copied == l_name as u64,
            "BAM reference name truncated: expected {l_name} bytes, got {copied}"
        );
        reader.read_exact(&mut buf4)?;
    }

    Ok(())
}

// ============================================================================
// Template-key variant selection
// ============================================================================

/// Which optional lanes a chosen template key retains.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct TemplateKeyVariant {
    /// Retain the `cb_hash` lane.
    pub cb: bool,
    /// Retain the tertiary (library|mi) lane.
    pub tertiary: bool,
}

impl TemplateKeyVariant {
    /// Number of u64 lanes: 3 + cb + tertiary.
    #[must_use]
    pub fn lanes(self) -> usize {
        3 + usize::from(self.cb) + usize::from(self.tertiary)
    }
}

/// Parsed `--key-types` spec controlling which optional sort-key lanes are kept.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub enum KeyTypesSpec {
    /// Auto-detect from the first record (default).
    #[default]
    Auto,
    /// Force all lanes (N=5).
    Full,
    /// Force no optional lanes (N=3).
    None,
    /// Force a specific optional set.
    Explicit {
        /// Retain the `cb_hash` lane.
        cb: bool,
        /// Retain the tertiary (library|mi) lane.
        tertiary: bool,
    },
}

/// Mask for the MI component of the `tertiary` key lane (low 48 bits).
/// `tertiary` packs `library_ordinal << 48 | (mi_value << 1 | !mi_suffix)`, so the
/// high 16 bits are the library ordinal and the low 48 bits encode MI. A nonzero
/// low-48 region means the first record carries an MI tag.
const TERTIARY_MI_MASK: u64 = (1u64 << 48) - 1;

/// Choose the template-key variant for this sort.
///
/// `first_key` is the full key extracted from the first record (`None` for an
/// empty input — then the variant is irrelevant; default to lite). `spec` is the
/// parsed `--key-types` override. `header_library_varies` is `true` when the
/// header realizes more than one distinct library ordinal (see
/// [`LibraryLookup::distinct_header_ordinals`]).
///
/// Selection only *provisions* the variant; the decode-time verify (fill loop)
/// still *guarantees* the dropped lanes are constant. Under `Auto`, the tertiary
/// library lane is kept only when it can actually vary: the first record carries
/// an MI tag, or the header realizes more than one library ordinal. A single,
/// constant library ordinal in the high bits is droppable — the header informs
/// provisioning, while verify backstops the rare case of a record whose RG is
/// absent from the header (which realizes ordinal 0).
#[must_use]
pub fn select_template_variant(
    first_key: Option<&TemplateKey>,
    spec: KeyTypesSpec,
    header_library_varies: bool,
) -> TemplateKeyVariant {
    match spec {
        KeyTypesSpec::Full => TemplateKeyVariant { cb: true, tertiary: true },
        KeyTypesSpec::None => TemplateKeyVariant { cb: false, tertiary: false },
        KeyTypesSpec::Explicit { cb, tertiary } => TemplateKeyVariant { cb, tertiary },
        KeyTypesSpec::Auto => match first_key {
            None => TemplateKeyVariant { cb: false, tertiary: false },
            Some(k) => {
                // Keep the tertiary lane only if it can actually vary: MI present on
                // the first record (low-48 bits), or the header realizes more than
                // one library ordinal. A single (constant) library ordinal in the
                // high bits is droppable — decode-time verify backstops the rare
                // case of a record whose RG is absent from the header (ordinal 0).
                let mi_present = (k.tertiary & TERTIARY_MI_MASK) != 0;
                TemplateKeyVariant {
                    cb: k.cb_hash != 0,
                    tertiary: mi_present || header_library_varies,
                }
            }
        },
    }
}

// ============================================================================
// Dropped-lane violation detection
// ============================================================================

/// A dropped-lane violation discovered during decode-time verify.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DroppedLaneViolation {
    /// A record carried a CB tag absent from the first record.
    Cb,
    /// A record carried an MI value differing from the first record.
    Mi,
    /// A record's realized library ordinal differed from the first record.
    Library,
}

impl DroppedLaneViolation {
    /// The `--key-types` token that re-includes this lane.
    #[must_use]
    pub fn key_types_token(self) -> &'static str {
        match self {
            DroppedLaneViolation::Cb => "cb",
            DroppedLaneViolation::Mi => "mi",
            DroppedLaneViolation::Library => "library",
        }
    }
}

/// Verify that the lanes the chosen `variant` dropped are constant for `cur`
/// relative to the sort's `first` full key. Returns the first violation found,
/// or `None` if all dropped lanes match.
///
/// Tertiary sub-fields are decoded only to attribute the violation: the high 16
/// bits are the library ordinal; the low 48 bits are `(mi_value << 1) | !mi_suffix`
/// (`mi_value` occupies bits 47..1; `!mi_suffix` is bit 0).
#[must_use]
pub fn verify_dropped_lanes(
    first: &TemplateKey,
    cur: &TemplateKey,
    variant: TemplateKeyVariant,
) -> Option<DroppedLaneViolation> {
    if !variant.cb && cur.cb_hash != first.cb_hash {
        return Some(DroppedLaneViolation::Cb);
    }
    if !variant.tertiary && cur.tertiary != first.tertiary {
        // Attribute: library is the high 16 bits, mi the low 48.
        let lib_changed = (cur.tertiary >> 48) != (first.tertiary >> 48);
        return Some(if lib_changed {
            DroppedLaneViolation::Library
        } else {
            DroppedLaneViolation::Mi
        });
    }
    None
}

/// Build the actionable error message for a dropped-lane violation.
pub(crate) fn dropped_lane_error(name: &str, v: DroppedLaneViolation) -> anyhow::Error {
    let field = match v {
        DroppedLaneViolation::Cb => "CB",
        DroppedLaneViolation::Mi => "MI",
        DroppedLaneViolation::Library => "library",
    };
    anyhow::anyhow!(
        "record {name} carries a {field} value absent from the input's first record; \
         re-run with --key-types {}",
        v.key_types_token(),
    )
}

/// Open `path` as a single [`SortMergeSlot`] with the given `file_id`.
///
/// The spill codec is auto-detected from the file magic: for zstd (`ZSP1`) the
/// 4-byte magic is consumed here, leaving the reader at the first `[len][frame]`
/// record; for BGZF the `1f 8b` is part of the first block, so the reader is
/// rewound to byte 0. `SortSpillDecompress` reads `slot.codec` to decompress
/// accordingly.
///
/// `file_id` is the slot's stable merge-order identifier — `SortMerge` orders
/// sources by it so the `LoserTree` tie-break for equal sort keys matches the
/// legacy chunk-files order. The P6 `CompressSpill` step passes a chunk's
/// **logical spill index** here so the tie-break is independent of which worker
/// happened to write the file.
///
/// # Errors
///
/// Returns an error if the file cannot be opened, its magic cannot be read, or
/// (for BGZF) the rewind fails.
pub fn open_spill_slot(path: &std::path::Path, file_id: u32) -> Result<Arc<SortMergeSlot>> {
    let mut file = std::fs::File::open(path)
        .with_context(|| format!("failed to open spill chunk {}", path.display()))?;
    let mut magic = [0u8; 4];
    let filled = read_exact_or_eof(&mut file, &mut magic)
        .with_context(|| format!("failed to read spill chunk magic {}", path.display()))?;
    let codec = if filled {
        crate::codec::SpillCodec::from_magic(&magic).unwrap_or(crate::codec::SpillCodec::Bgzf)
    } else {
        crate::codec::SpillCodec::Bgzf
    };
    if matches!(codec, crate::codec::SpillCodec::Bgzf) {
        use std::io::Seek;
        file.seek(std::io::SeekFrom::Start(0))
            .with_context(|| format!("failed to rewind spill chunk {}", path.display()))?;
    }
    let reader = std::io::BufReader::new(file);
    Ok(Arc::new(SortMergeSlot::new(file_id, reader, codec)))
}

// ============================================================================
// CoordinateSortStream — streaming coordinate-sort handle
// ============================================================================

// ============================================================================
// QuerynameSortStream — streaming queryname-sort handle (lex or natural)
// ============================================================================

// ============================================================================
// TemplateCoordinateSortStream — streaming template-coordinate-sort handle
// ============================================================================

// ============================================================================
// Streaming-sort merge path (unified pipeline `SortMerge` step)
// ============================================================================
//
// Ported from the issue-#330 branch onto main's sort engine. This is the
// slot-backed, *non-blocking* merge driver used exclusively by the typed-step
// pipeline (`MergeDriver::from_slots`). It reads records out of
// `SortMergeSlot` queues, which the `SortSpillDecompress` step has *already*
// filled with decompressed `Vec<u8>` blocks — so this path is **codec-agnostic**
// (it never touches disk or the spill codec) and carries no chunk-corruption
// risk. The eager pool-integrated `MergeDriver::new` path from #330 is
// intentionally NOT ported: standalone `fgumi sort` uses main's own
// `merge_chunks_*` engine, and the pipeline only ever calls `from_slots`.

/// Result of a non-blocking source read. `WouldBlock` is only ever produced by
/// `Slot` sources whose decompressed queue is momentarily empty and not yet at
/// EOF; the in-memory source always returns `Ready`.
pub(crate) enum TryRead<K> {
    /// A record was read (`Some`) or the source is at clean EOF (`None`).
    Ready(Option<K>),
    /// The source has no bytes available right now but is not at EOF; the
    /// caller should yield and retry on a later dispatch.
    WouldBlock,
}

/// Per-slot byte-stream parser state for the non-blocking merge reader.
///
/// As decompressed blocks are popped from a `SortMergeSlot`, their bytes are
/// stashed in `current_buf` and consumed left-to-right by the record parser.
/// Named distinctly from main's own `SourceParserState` (the pool-merge path)
/// because this one additionally carries resumable cross-block framing state.
struct SlotParserState {
    /// Current decompressed block being consumed.
    current_buf: Vec<u8>,
    /// Read position within `current_buf`.
    current_pos: usize,
    /// Resumable framing state for [`slot_try_next_record`]. `Some` while a
    /// record spans a block boundary and the next block was not yet available,
    /// so the read returned `WouldBlock` mid-record; the next call resumes from
    /// here.
    pending: Option<PendingRecord>,
}

/// Partially-read record, retained across a `WouldBlock` so the non-blocking
/// slot reader can resume. The wire layout is `[key:ksize][len:4][body:len]`
/// for non-embedded keys (`ksize = K::SERIALIZED_SIZE`) and `[len:4][body:len]`
/// for embedded keys (`ksize = 0`, key extracted from the body). We collect the
/// key prefix (if any), then the length prefix, then the body.
#[derive(Default)]
struct PendingRecord {
    /// Expected serialized key size (`0` for embedded keys).
    key_size: usize,
    /// Key-prefix bytes collected so far (length grows to `key_size`).
    key_buf: Vec<u8>,
    /// Length-prefix bytes collected so far (`0..=4`).
    len_buf: [u8; 4],
    /// Count of `len_buf` bytes filled.
    len_got: usize,
    /// Total body length, known once `len_got == 4`. The body accumulates
    /// directly into the caller's record buffer; progress is `out.len()`.
    body_len: Option<usize>,
}

impl SlotParserState {
    fn new() -> Self {
        Self { current_buf: Vec::new(), current_pos: 0, pending: None }
    }

    fn remaining(&self) -> usize {
        self.current_buf.len() - self.current_pos
    }
}

/// Outcome of a non-blocking attempt to pull the next decompressed block from a
/// slot into `parser.current_buf`.
enum BlockLoad {
    /// A block was moved into `parser.current_buf`.
    Loaded,
    /// The slot is fully delivered (queue empty AND `queue_eof`).
    Eof,
    /// The queue is empty but the producer is still feeding this slot.
    WouldBlock,
}

/// Pop one decompressed block (FIFO) into `parser`, surfacing a decompression
/// error in preference to EOF, reporting clean EOF, or `WouldBlock` while the
/// producer is still active. Never parks.
fn slot_try_load_block(
    slot: &Arc<crate::merge_slots::SortMergeSlot>,
    parser: &mut SlotParserState,
) -> Result<BlockLoad> {
    use std::sync::atomic::Ordering;

    let mut guard = slot.decompressed.lock().expect("slot.decompressed mutex poisoned");
    if let Some(data) = guard.pop_front() {
        drop(guard);
        parser.current_buf = data;
        parser.current_pos = 0;
        return Ok(BlockLoad::Loaded);
    }
    // Queue empty. Both atomics are visible under the lock release-acquire
    // chain (producer stores them while holding this same mutex). Surface
    // error in preference to EOF.
    if slot.decomp_error.load(Ordering::Acquire) {
        drop(guard);
        anyhow::bail!("spill decompression error on slot {}", slot.file_id);
    }
    if slot.queue_eof.load(Ordering::Acquire) {
        return Ok(BlockLoad::Eof);
    }
    Ok(BlockLoad::WouldBlock)
}

/// Serialized key size on the wire for a slot source: `0` for embedded keys
/// (key extracted from the body), else `K::SERIALIZED_SIZE`. Variable-length
/// non-embedded keys are unsupported on the slot path (none exist in
/// production — every slot-path key is embedded).
fn slot_key_size<K: RawSortKey>() -> Result<usize> {
    if K::EMBEDDED_IN_RECORD {
        Ok(0)
    } else {
        K::SERIALIZED_SIZE.ok_or_else(|| {
            anyhow::anyhow!("non-embedded slot key must have a fixed serialized size")
        })
    }
}

/// Parse a sort key from `bytes`: extracted from the body for embedded keys, or
/// deserialized from the `key_size`-byte prefix for non-embedded keys.
fn slot_parse_key<K: RawSortKey + Default + 'static>(key_bytes: &[u8], body: &[u8]) -> Result<K> {
    if K::EMBEDDED_IN_RECORD {
        Ok(K::extract_from_record(body))
    } else {
        let mut r = key_bytes;
        K::read_from(&mut r).map_err(|e| anyhow::anyhow!("malformed slot key: {e}"))
    }
}

/// Non-blocking, resumable read of the next record from a slot.
///
/// Returns `Ready(Some(key))` with the record body in `out`, `Ready(None)` at
/// clean source EOF, or `WouldBlock` when a needed block is not yet available
/// (progress is saved in `parser.pending` so the next call resumes).
///
/// # Errors
///
/// Truncation (EOF mid-record), a malformed key, or a producer-side
/// decompression error.
fn slot_try_next_record<K: RawSortKey + Default + 'static>(
    slot: &Arc<crate::merge_slots::SortMergeSlot>,
    parser: &mut SlotParserState,
    out: &mut Vec<u8>,
) -> Result<TryRead<K>> {
    let key_size = slot_key_size::<K>()?;

    if parser.pending.is_none() {
        // Ensure at least one byte is available, or detect clean EOF at a
        // record boundary.
        if parser.remaining() == 0 {
            match slot_try_load_block(slot, parser)? {
                BlockLoad::Loaded => {}
                BlockLoad::Eof => return Ok(TryRead::Ready(None)),
                BlockLoad::WouldBlock => return Ok(TryRead::WouldBlock),
            }
        }
        // Fast path: the entire record lives in the current block → parse in
        // place with a single copy of the body into `out`.
        let header = key_size + 4;
        if parser.remaining() >= header {
            let p = parser.current_pos;
            let len = u32::from_le_bytes(
                parser.current_buf[p + key_size..p + header].try_into().expect("4-byte slice"),
            ) as usize;
            if parser.remaining() >= header + len {
                out.clear();
                out.extend_from_slice(&parser.current_buf[p + header..p + header + len]);
                let key = slot_parse_key::<K>(&parser.current_buf[p..p + key_size], out)?;
                parser.current_pos += header + len;
                return Ok(TryRead::Ready(Some(key)));
            }
        }
        // Slow path: record spans a block boundary (or fewer than `header`
        // bytes remain). Begin a resumable pending record.
        parser.pending = Some(PendingRecord { key_size, ..PendingRecord::default() });
    }

    slot_collect_pending::<K>(slot, parser, out)
}

/// Resumable continuation of [`slot_try_next_record`]'s slow path. Collects the
/// key prefix (if any), the 4-byte length prefix, then the body across as many
/// blocks as needed, returning `WouldBlock` (state retained in
/// `parser.pending`) whenever the next block isn't ready yet.
fn slot_collect_pending<K: RawSortKey + Default + 'static>(
    slot: &Arc<crate::merge_slots::SortMergeSlot>,
    parser: &mut SlotParserState,
    out: &mut Vec<u8>,
) -> Result<TryRead<K>> {
    // Stage 0: key prefix (non-embedded keys only; `key_size == 0` is a no-op).
    loop {
        let (key_size, key_have) = {
            let pending = parser.pending.as_ref().expect("pending set");
            (pending.key_size, pending.key_buf.len())
        };
        let need = key_size - key_have;
        if need == 0 {
            break;
        }
        if parser.current_pos == parser.current_buf.len() {
            match slot_try_load_block(slot, parser)? {
                BlockLoad::Loaded => {}
                BlockLoad::Eof => anyhow::bail!(
                    "truncated record key in slot {} ({key_have} of {key_size} key bytes at EOF)",
                    slot.file_id,
                ),
                BlockLoad::WouldBlock => return Ok(TryRead::WouldBlock),
            }
        }
        let avail = parser.current_buf.len() - parser.current_pos;
        let take = need.min(avail);
        let src = &parser.current_buf[parser.current_pos..parser.current_pos + take];
        parser.pending.as_mut().expect("pending set").key_buf.extend_from_slice(src);
        parser.current_pos += take;
    }

    // Stage 1: length prefix.
    loop {
        let len_got = parser.pending.as_ref().expect("pending set").len_got;
        if len_got == 4 {
            break;
        }
        if parser.current_pos == parser.current_buf.len() {
            match slot_try_load_block(slot, parser)? {
                BlockLoad::Loaded => {}
                BlockLoad::Eof => anyhow::bail!(
                    "truncated record length in slot {} ({len_got} of 4 length bytes at EOF)",
                    slot.file_id,
                ),
                BlockLoad::WouldBlock => return Ok(TryRead::WouldBlock),
            }
        }
        let avail = parser.current_buf.len() - parser.current_pos;
        let take = (4 - len_got).min(avail);
        let src = &parser.current_buf[parser.current_pos..parser.current_pos + take];
        let pending = parser.pending.as_mut().expect("pending set");
        pending.len_buf[len_got..len_got + take].copy_from_slice(src);
        pending.len_got += take;
        parser.current_pos += take;
    }

    // Compute body length once and reset `out` for accumulation.
    {
        let pending = parser.pending.as_mut().expect("pending set");
        if pending.body_len.is_none() {
            #[allow(clippy::cast_possible_truncation)]
            let len = u32::from_le_bytes(pending.len_buf) as usize;
            pending.body_len = Some(len);
            out.clear();
            out.reserve(len);
        }
    }

    // Stage 2: body.
    loop {
        let body_len = parser.pending.as_ref().expect("pending set").body_len.expect("len known");
        if out.len() >= body_len {
            break;
        }
        if parser.current_pos == parser.current_buf.len() {
            match slot_try_load_block(slot, parser)? {
                BlockLoad::Loaded => {}
                BlockLoad::Eof => anyhow::bail!(
                    "truncated record body in slot {} ({} of {body_len} bytes at EOF)",
                    slot.file_id,
                    out.len(),
                ),
                BlockLoad::WouldBlock => return Ok(TryRead::WouldBlock),
            }
        }
        let avail = parser.current_buf.len() - parser.current_pos;
        let take = (body_len - out.len()).min(avail);
        out.extend_from_slice(&parser.current_buf[parser.current_pos..parser.current_pos + take]);
        parser.current_pos += take;
    }

    let key = {
        let pending = parser.pending.as_ref().expect("pending set");
        slot_parse_key::<K>(&pending.key_buf, out)?
    };
    parser.pending = None;
    Ok(TryRead::Ready(Some(key)))
}

/// Merge source for the slot-backed `MergeDriver`: either a decompressing slot
/// (read non-blocking) or an already-sorted in-memory residual chunk.
enum SlotMergeSource<K: RawSortKey + Default + 'static> {
    /// Pipeline-integrated source — `SortSpillDecompress` fills
    /// `slot.decompressed` with decompressed blocks; this driver reads them
    /// in-order via the per-source [`SlotParserState`].
    Slot { slot: Arc<crate::merge_slots::SortMergeSlot>, parser: SlotParserState },
    /// In-memory sorted residual records from a queryname-style sort (each
    /// record an owned `RawRecord`). Read via a zero-copy `mem::swap` bridge.
    Memory { records: Vec<(K, fgumi_raw_bam::RawRecord)>, idx: usize },
    /// In-memory sorted residual from an inline-buffer (coordinate / template)
    /// sort, sharing an `Arc<SegmentedBuf>` — zero per-record allocation. Read
    /// copies each record's bytes into the caller's reused `buf`.
    MemoryShared { chunk: InMemoryChunk<K>, idx: usize },
}

impl<K: RawSortKey + Default + 'static> SlotMergeSource<K> {
    /// Non-blocking read of the next record into `buf`, returning the sort key.
    /// `Slot` sources may return `WouldBlock`; `Memory` always returns `Ready`.
    fn try_next_record(&mut self, buf: &mut Vec<u8>) -> Result<TryRead<K>> {
        match self {
            SlotMergeSource::Slot { slot, parser } => slot_try_next_record::<K>(slot, parser, buf),
            SlotMergeSource::Memory { records, idx } => {
                if *idx < records.len() {
                    let (ref mut key, ref mut data) = records[*idx];
                    // Bridge: RawRecord wraps Vec<u8>; swap via the inner vec to
                    // avoid re-allocating. The caller's buf is a plain Vec<u8>.
                    std::mem::swap(buf, data.as_mut_vec());
                    let key = std::mem::take(key);
                    *idx += 1;
                    Ok(TryRead::Ready(Some(key)))
                } else {
                    Ok(TryRead::Ready(None))
                }
            }
            SlotMergeSource::MemoryShared { chunk, idx } => {
                if *idx < chunk.len() {
                    // Bytes are borrowed from the shared `Arc<SegmentedBuf>`, so
                    // copy them into the caller's reused `buf` (one memcpy, no
                    // per-record allocation — the materialise-time copy is gone).
                    buf.clear();
                    buf.extend_from_slice(chunk.record_bytes(*idx));
                    let key = chunk.take_key(*idx);
                    *idx += 1;
                    Ok(TryRead::Ready(Some(key)))
                } else {
                    Ok(TryRead::Ready(None))
                }
            }
        }
    }
}

/// One step of a non-blocking k-way merge.
pub enum MergeStep<'a> {
    /// The next merged record. The bytes borrow the driver's internal winner
    /// buffer and are valid only until the next `try_step` call — the caller
    /// must copy them out before stepping again. The borrow checker enforces
    /// this, which is what makes the deferred-refill protocol sound.
    Produced(&'a [u8]),
    /// A source's decompressed queue is momentarily empty (and not at EOF), or
    /// a source hasn't been primed yet. No record is available right now; the
    /// caller should yield and retry on a later dispatch, by which point the
    /// producer will have refilled the slot.
    Stalled,
    /// The merge is exhausted.
    Done,
}

impl std::fmt::Debug for MergeStep<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Produced(bytes) => write!(f, "Produced({} bytes)", bytes.len()),
            Self::Stalled => f.write_str("Stalled"),
            Self::Done => f.write_str("Done"),
        }
    }
}

/// Object-safe view over `MergeDriver<K>` so the pipeline `SortMerge` step can
/// hold a single concrete type regardless of sort order.
pub trait MergeDriverDyn: Send {
    /// Produce the next merged record **without blocking**. Any step that needs
    /// a not-yet-ready slot block returns [`MergeStep::Stalled`].
    ///
    /// # Errors
    ///
    /// Propagates spill-file corruption or decompression errors from the
    /// underlying slot sources.
    fn try_step(&mut self) -> Result<MergeStep<'_>>;
    /// Total records emitted since the driver was constructed.
    fn records_merged(&self) -> u64;
}

/// Internal phase of the [`MergeDriver`] non-blocking state machine.
enum MergePhase<K> {
    /// Lazily priming the loser tree: each source's first record is pulled one
    /// at a time (non-blocking). `next` is the index of the next source to
    /// prime; `keys` accumulates the primed initial keys; `rec` is the partial
    /// record buffer for the source currently being primed (retained across a
    /// `Stalled` so a spanning first record resumes correctly).
    Priming { next: usize, keys: Vec<K>, rec: Vec<u8> },
    /// Active k-way merge. `pending_refill` is set after a record is emitted
    /// (`Produced`) and means the winner's source must be refilled before the
    /// next winner is produced.
    Merging { tree: crate::loser_tree::LoserTree<K>, pending_refill: bool },
    /// Exhausted.
    Done,
}

/// Pausable state machine for a non-blocking k-way merge over slot + in-memory
/// sources. Constructed via [`MergeDriver::from_slots`]; driven by the pipeline
/// `SortMerge` step via [`MergeDriverDyn::try_step`].
pub struct MergeDriver<K: RawSortKey + Default + Send + 'static> {
    sources: Vec<SlotMergeSource<K>>,
    /// Per-active-source current record buffer, indexed by loser-tree leaf
    /// (parallel to `source_map`). Grows during `Priming`.
    records: Vec<Vec<u8>>,
    /// Loser-tree leaf → `sources` index. Built during `Priming`.
    source_map: Vec<usize>,
    /// Non-blocking merge state machine.
    phase: MergePhase<K>,
    records_merged: u64,
    progress: ProgressTracker,
}

impl<K: RawSortKey + Default + Send + 'static> MergeDriver<K> {
    /// Construct a driver from `SortMergeSlot`s and already-sorted in-memory
    /// chunks. Used by the cooperative pipeline `SortMerge` step.
    ///
    /// Each slot becomes a slot source read **non-blocking** by
    /// [`MergeDriverDyn::try_step`]. Priming is **lazy**: the loser tree is
    /// built on the first `try_step` call(s) as each source yields its first
    /// record, so construction never blocks on a not-yet-decompressed slot.
    /// Construction is infallible (read errors surface from `try_step`), and
    /// the empty-input case is reported by the first `try_step` returning
    /// [`MergeStep::Done`].
    #[must_use]
    pub fn from_slots(
        slots: Vec<Arc<crate::merge_slots::SortMergeSlot>>,
        memory: MemorySources<K>,
        total_records: u64,
    ) -> Self {
        let num_slots = slots.len();
        let num_memory = memory.num_non_empty();
        let num_sources = num_slots + num_memory;

        if num_slots > 0 {
            info!(
                "Pipeline-integrated merge: {num_slots} slot sources + {num_memory} memory sources"
            );
        }

        let mut sources: Vec<SlotMergeSource<K>> = Vec::with_capacity(num_sources);
        for slot in slots {
            sources.push(SlotMergeSource::Slot { slot, parser: SlotParserState::new() });
        }
        match memory {
            MemorySources::Owned(chunks) => {
                for chunk in chunks {
                    if !chunk.is_empty() {
                        sources.push(SlotMergeSource::Memory { records: chunk, idx: 0 });
                    }
                }
            }
            MemorySources::Shared(chunks) => {
                for chunk in chunks {
                    if !chunk.is_empty() {
                        sources.push(SlotMergeSource::MemoryShared { chunk, idx: 0 });
                    }
                }
            }
        }

        let progress = ProgressTracker::new("Merged records")
            .with_interval(1_000_000)
            .with_total(total_records);

        Self {
            sources,
            records: Vec::with_capacity(num_sources),
            source_map: Vec::with_capacity(num_sources),
            // Lazy priming: tree built on first try_step as sources yield.
            phase: MergePhase::Priming {
                next: 0,
                keys: Vec::with_capacity(num_sources),
                rec: Vec::new(),
            },
            records_merged: 0,
            progress,
        }
    }
}

impl<K: RawSortKey + Default + Send + 'static> MergeDriverDyn for MergeDriver<K> {
    fn try_step(&mut self) -> Result<MergeStep<'_>> {
        // Move the phase out so the body works with owned `tree`/`keys`/`rec`
        // and can borrow `self.{sources,records,source_map}` freely without
        // aliasing the phase; the phase is written back before every return.
        // `Done` is the placeholder; any early `?` leaves the driver `Done`,
        // which is safe because the whole merge aborts on error.
        loop {
            match std::mem::replace(&mut self.phase, MergePhase::Done) {
                MergePhase::Priming { mut next, mut keys, mut rec } => {
                    let mut stalled = false;
                    while next < self.sources.len() {
                        match self.sources[next].try_next_record(&mut rec)? {
                            TryRead::WouldBlock => {
                                stalled = true;
                                break;
                            }
                            TryRead::Ready(Some(key)) => {
                                keys.push(key);
                                self.records.push(std::mem::take(&mut rec));
                                self.source_map.push(next);
                                next += 1;
                            }
                            TryRead::Ready(None) => {
                                rec.clear();
                                next += 1;
                            }
                        }
                    }
                    if stalled {
                        self.phase = MergePhase::Priming { next, keys, rec };
                        return Ok(MergeStep::Stalled);
                    }
                    if keys.is_empty() {
                        return Ok(MergeStep::Done);
                    }
                    info!("Merging from {} sources...", keys.len());
                    let tree = crate::loser_tree::LoserTree::new(keys);
                    self.phase = MergePhase::Merging { tree, pending_refill: false };
                    // Loop to emit the first winner.
                }
                MergePhase::Merging { mut tree, pending_refill } => {
                    // Resolve the refill deferred from the previous `Produced`
                    // (the read that would overwrite the just-emitted winner's
                    // buffer). On `WouldBlock`, yield with the refill still
                    // pending.
                    if pending_refill && tree.winner_is_active() {
                        let winner = tree.winner();
                        let src = self.source_map[winner];
                        match self.sources[src].try_next_record(&mut self.records[winner])? {
                            TryRead::WouldBlock => {
                                self.phase = MergePhase::Merging { tree, pending_refill: true };
                                return Ok(MergeStep::Stalled);
                            }
                            TryRead::Ready(Some(key)) => tree.replace_winner(key),
                            TryRead::Ready(None) => tree.remove_winner(),
                        }
                    }
                    if !tree.winner_is_active() {
                        return Ok(MergeStep::Done);
                    }
                    let winner = tree.winner();
                    self.records_merged += 1;
                    self.progress.log_if_needed(1);
                    // Defer this winner's refill to the next call so the
                    // returned borrow stays valid until the caller copies it.
                    self.phase = MergePhase::Merging { tree, pending_refill: true };
                    return Ok(MergeStep::Produced(&self.records[winner]));
                }
                MergePhase::Done => return Ok(MergeStep::Done),
            }
        }
    }

    fn records_merged(&self) -> u64 {
        self.records_merged
    }
}

impl<K: RawSortKey + Default + Send + 'static> Drop for MergeDriver<K> {
    fn drop(&mut self) {
        self.progress.log_final();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bstr::BString;
    use fgumi_sam::{PairBuilder, SamBuilder};
    use noodles::sam::header::record::value::Map;
    use noodles::sam::header::record::value::map::ReadGroup;

    // ========================================================================
    // LibraryLookup tests
    // ========================================================================

    #[test]
    fn test_library_lookup_empty_header() {
        let header = Header::builder().build();
        let lookup = LibraryLookup::from_header(&header);
        assert!(lookup.rg_to_ordinal.is_empty());
    }

    #[test]
    fn test_library_lookup_single_rg() {
        let rg = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("LibA"))
            .build()
            .expect("valid");
        let header = Header::builder().add_read_group(BString::from("rg1"), rg).build();

        let lookup = LibraryLookup::from_header(&header);
        assert_eq!(lookup.rg_to_ordinal.len(), 1);
        // LibA is the only library, so it gets ordinal 1 (0 is reserved for empty/unknown)
        assert_eq!(
            *lookup.rg_to_ordinal.get(b"rg1".as_slice()).expect("rg1 should be in ordinal map"),
            1
        );
    }

    #[test]
    fn test_library_lookup_multiple_libraries() {
        let rg_a = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("LibC"))
            .build()
            .expect("valid");
        let rg_b = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("LibA"))
            .build()
            .expect("valid");
        let rg_c = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("LibB"))
            .build()
            .expect("valid");

        let header = Header::builder()
            .add_read_group(BString::from("rg1"), rg_a)
            .add_read_group(BString::from("rg2"), rg_b)
            .add_read_group(BString::from("rg3"), rg_c)
            .build();

        let lookup = LibraryLookup::from_header(&header);
        assert_eq!(lookup.rg_to_ordinal.len(), 3);

        // Libraries sorted alphabetically: LibA=1, LibB=2, LibC=3
        let rg2 = *lookup.rg_to_ordinal.get(b"rg2".as_slice()).expect("rg2");
        let rg3 = *lookup.rg_to_ordinal.get(b"rg3".as_slice()).expect("rg3");
        let rg1 = *lookup.rg_to_ordinal.get(b"rg1".as_slice()).expect("rg1");
        assert_eq!(rg2, 1); // LibA
        assert_eq!(rg3, 2); // LibB
        assert_eq!(rg1, 3); // LibC
    }

    #[test]
    fn test_distinct_header_ordinals() {
        // (i) no read groups -> 0.
        let empty = LibraryLookup::from_header(&Header::builder().build());
        assert_eq!(empty.distinct_header_ordinals(), 0);

        // (ii) one RG with an LB -> 1.
        let rg_a = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("LibA"))
            .build()
            .expect("valid");
        let one = LibraryLookup::from_header(
            &Header::builder().add_read_group(BString::from("rg1"), rg_a).build(),
        );
        assert_eq!(one.distinct_header_ordinals(), 1);

        // (iii) two RGs with two different LBs -> 2.
        let rg_b = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("LibA"))
            .build()
            .expect("valid");
        let rg_c = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("LibB"))
            .build()
            .expect("valid");
        let two = LibraryLookup::from_header(
            &Header::builder()
                .add_read_group(BString::from("rg1"), rg_b)
                .add_read_group(BString::from("rg2"), rg_c)
                .build(),
        );
        assert_eq!(two.distinct_header_ordinals(), 2);

        // (iv) two RGs with the SAME LB -> 1.
        let rg_d = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("LibA"))
            .build()
            .expect("valid");
        let rg_e = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("LibA"))
            .build()
            .expect("valid");
        let same = LibraryLookup::from_header(
            &Header::builder()
                .add_read_group(BString::from("rg1"), rg_d)
                .add_read_group(BString::from("rg2"), rg_e)
                .build(),
        );
        assert_eq!(same.distinct_header_ordinals(), 1);

        // (v) one RG with an LB + one RG WITHOUT an LB -> 2 (ordinals {1, 0}).
        let rg_f = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("LibA"))
            .build()
            .expect("valid");
        let rg_no_lb = Map::<ReadGroup>::builder().build().expect("valid");
        let mixed = LibraryLookup::from_header(
            &Header::builder()
                .add_read_group(BString::from("rg1"), rg_f)
                .add_read_group(BString::from("rg2"), rg_no_lb)
                .build(),
        );
        assert_eq!(mixed.distinct_header_ordinals(), 2);
    }

    #[test]
    fn test_library_lookup_unknown_rg_returns_zero() {
        let rg = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("LibA"))
            .build()
            .expect("valid");
        let header = Header::builder().add_read_group(BString::from("rg1"), rg).build();

        let lookup = LibraryLookup::from_header(&header);
        // A BAM record with no RG tag or unknown RG should return ordinal 0
        // We can test get_ordinal with a minimal BAM record that has no aux data
        let mut bam = vec![0u8; 36];
        bam[8] = 4; // l_read_name = 4
        bam[32..36].copy_from_slice(b"rea\0");
        assert_eq!(lookup.get_ordinal(&bam), 0);
    }

    // ========================================================================
    // RawExternalSorter builder tests
    // ========================================================================

    #[test]
    fn test_raw_sorter_defaults() {
        let sorter = RawExternalSorter::new(SortOrder::Coordinate);
        assert_eq!(sorter.memory_limit, 512 * 1024 * 1024);
        assert!(sorter.temp_dirs.is_empty());
        assert_eq!(sorter.threads, 1);
        assert_eq!(sorter.output_compression, 6);
        assert_eq!(sorter.temp_compression, 1);
        assert!(sorter.pg_info.is_none());
        assert_eq!(sorter.max_temp_files, DEFAULT_MAX_TEMP_FILES);
    }

    #[test]
    fn test_raw_sorter_builder_chain() {
        let sorter = RawExternalSorter::new(SortOrder::Queryname(QuerynameComparator::default()))
            .memory_limit(1024)
            .temp_dir(PathBuf::from("/tmp/test"))
            .threads(8)
            .output_compression(9)
            .temp_compression(3)
            .pg_info("1.0".to_string(), "fgumi sort".to_string())
            .max_temp_files(128);

        assert_eq!(sorter.memory_limit, 1024);
        assert_eq!(sorter.temp_dirs, vec![PathBuf::from("/tmp/test")]);
        assert_eq!(sorter.threads, 8);
        assert_eq!(sorter.output_compression, 9);
        assert_eq!(sorter.temp_compression, 3);
        assert_eq!(sorter.pg_info, Some(("1.0".to_string(), "fgumi sort".to_string())));
        assert_eq!(sorter.max_temp_files, 128);
    }

    #[test]
    fn test_raw_sorter_phase_threads_builders() {
        let new = || RawExternalSorter::new(SortOrder::Coordinate).threads(8);
        // Default: both phases fall back to threads.
        let base = new();
        assert_eq!((base.sort_threads, base.merge_threads), (None, None));
        assert_eq!((base.phase1_threads(), base.phase2_threads()), (8, 8));
        // sort_threads overrides only Phase 1.
        let sort_only = new().sort_threads(2);
        assert_eq!((sort_only.phase1_threads(), sort_only.phase2_threads()), (2, 8));
        // The public `num_threads()` accessor returns the BASE count, unaffected
        // by the Phase-1 override — so the streaming front must read
        // `phase1_threads()`, not `num_threads()`, to honor `--sort-threads`
        // (the D1.1 regression was reading the base here). Pin the distinction.
        assert_eq!(sort_only.num_threads(), 8);
        assert_eq!(sort_only.phase1_threads(), 2);
        // merge_threads overrides only Phase 2.
        let merge_only = new().merge_threads(3);
        assert_eq!((merge_only.phase1_threads(), merge_only.phase2_threads()), (8, 3));
        // Both override independently.
        let both = new().sort_threads(2).merge_threads(16);
        assert_eq!((both.phase1_threads(), both.phase2_threads()), (2, 16));
        // Clamp to >= 1.
        let zero = new().sort_threads(0).merge_threads(0);
        assert_eq!((zero.phase1_threads(), zero.phase2_threads()), (1, 1));
    }

    /// Read every record's raw bytes from a BAM, in file order.
    fn records_bytes(path: &Path) -> Vec<Vec<u8>> {
        let (mut reader, _h) = create_raw_bam_reader(path, 1).expect("open bam");
        let mut out = Vec::new();
        let mut rec = fgumi_raw_bam::RawRecord::new();
        while reader.read_record(&mut rec).expect("read record") != 0 {
            out.push(rec.view().as_bytes().to_vec());
        }
        out
    }

    /// Splitting Phase-1 (`sort_threads`) and Phase-2 (`merge_threads`) from the
    /// base `threads` is a pure scheduling knob: output must be byte-identical to
    /// the unsplit sort, with spills forced. Matrixed over every order (`#[values]`)
    /// and every thread-split variant (`#[case]`) so a failure names the exact order
    /// and split. The cases mirror the prior inline matrix: `sort_low` (Phase-1 on
    /// 1), `merge_low` (Phase-2 on 1), `both` (Phase-1 on 1 / Phase-2 on 2), and
    /// `sort_high` (Phase-1 above base).
    #[rstest::rstest]
    #[case::sort_low(Some(1), None)]
    #[case::merge_low(None, Some(1))]
    #[case::both(Some(1), Some(2))]
    #[case::sort_high(Some(8), Some(2))]
    fn phase_thread_split_is_byte_identical(
        #[values(
            SortOrder::Coordinate,
            SortOrder::Queryname(QuerynameComparator::Natural),
            SortOrder::Queryname(QuerynameComparator::Lexicographic),
            SortOrder::TemplateCoordinate
        )]
        order: SortOrder,
        #[case] sort_t: Option<usize>,
        #[case] merge_t: Option<usize>,
    ) {
        use fgumi_sam::SamBuilder;
        let mut builder = SamBuilder::new();
        for i in 0..2000 {
            let _ = builder
                .add_pair()
                .name(&format!("read{i}"))
                .start1(i * 50 + 1)
                .start2(i * 50 + 101)
                .build();
        }
        let dir = tempfile::tempdir().expect("temp dir");
        let input = dir.path().join("in.bam");
        builder.write_bam(&input).expect("write bam");

        let mk = |out: &Path, sort_t: Option<usize>, merge_t: Option<usize>| {
            let mut s = RawExternalSorter::new(order)
                .threads(4)
                .memory_limit(32 * 1024) // force spills
                .spill_codec(crate::codec::SpillCodec::Bgzf)
                .temp_compression(0)
                .output_compression(0);
            if let Some(j) = sort_t {
                s = s.sort_threads(j);
            }
            if let Some(k) = merge_t {
                s = s.merge_threads(k);
            }
            s.sort(&input, out).expect("sort");
        };
        // Base: both phases on the unsplit `threads(4)`.
        let base = dir.path().join("base.bam");
        mk(&base, None, None);
        let golden = records_bytes(&base);
        // The split variant must reproduce the base output exactly.
        let out = dir.path().join("variant.bam");
        mk(&out, sort_t, merge_t);
        assert_eq!(
            records_bytes(&out),
            golden,
            "order {order:?} split sort_t={sort_t:?} merge_t={merge_t:?}: \
             thread split must not change output"
        );
    }

    /// Regression guard for the residual Phase-1 sort honoring `sort_threads`
    /// rather than the base `threads`. With `threads(1).sort_threads(4)` the
    /// residual sort/chunk decision must take the *parallel* `phase1_threads()`
    /// branch (previously it gated on `threads`, so it ran serially and ignored
    /// `sort_threads`). Spills are forced so the keyed-chunk residual path
    /// (`chunk_files` non-empty) is exercised; output must be byte-identical to
    /// the serial `threads(1)` sort for every order.
    #[rstest::rstest]
    #[case::coordinate(SortOrder::Coordinate)]
    #[case::queryname_natural(SortOrder::Queryname(QuerynameComparator::Natural))]
    #[case::queryname_lex(SortOrder::Queryname(QuerynameComparator::Lexicographic))]
    #[case::template_coordinate(SortOrder::TemplateCoordinate)]
    fn residual_phase1_honors_sort_threads_with_single_base_thread(#[case] order: SortOrder) {
        use fgumi_sam::SamBuilder;
        let mut builder = SamBuilder::new();
        for i in 0..2000 {
            let _ = builder
                .add_pair()
                .name(&format!("read{i}"))
                .start1(i * 50 + 1)
                .start2(i * 50 + 101)
                .build();
        }
        let dir = tempfile::tempdir().expect("temp dir");
        let input = dir.path().join("in.bam");
        builder.write_bam(&input).expect("write bam");

        let mk = |out: &Path, sort_t: Option<usize>| {
            let mut s = RawExternalSorter::new(order)
                .threads(1)
                .memory_limit(32 * 1024) // force spills -> exercise residual chunk path
                .spill_codec(crate::codec::SpillCodec::Bgzf)
                .temp_compression(0)
                .output_compression(0);
            if let Some(j) = sort_t {
                s = s.sort_threads(j);
            }
            s.sort(&input, out).expect("sort");
        };
        let serial = dir.path().join("serial.bam");
        mk(&serial, None); // threads(1): residual sort runs serially
        let golden = records_bytes(&serial);

        let parallel = dir.path().join("parallel.bam");
        mk(&parallel, Some(4)); // sort_threads(4): residual sort now parallel
        assert_eq!(
            records_bytes(&parallel),
            golden,
            "order {order:?}: threads(1).sort_threads(4) must match the serial sort output"
        );
    }

    #[test]
    fn test_raw_sorter_memory_limit() {
        let sorter = RawExternalSorter::new(SortOrder::Coordinate).memory_limit(256 * 1024 * 1024);
        assert_eq!(sorter.memory_limit, 256 * 1024 * 1024);
    }

    #[test]
    fn test_raw_sorter_temp_compression() {
        let sorter = RawExternalSorter::new(SortOrder::Coordinate).temp_compression(0);
        assert_eq!(sorter.temp_compression, 0);
    }

    /// `RawExternalSorter::sort` rejects `temp_compression=0` + `SpillCodec::Zstd`
    /// before doing any work, mirroring the CLI guard in `commands::sort`. zstd
    /// has no level-0 "stored" mode, and silently remapping to 1 would surprise
    /// API callers who pass 0 expecting uncompressed spills.
    #[test]
    fn test_raw_sorter_sort_rejects_temp_compression_zero_with_zstd() {
        use fgumi_sam::SamBuilder;

        let mut builder = SamBuilder::new();
        let _ = builder.add_pair().name("read0").start1(1).start2(101).build();

        let dir = tempfile::tempdir().expect("failed to create temp directory");
        let input = dir.path().join("input.bam");
        let output = dir.path().join("output.bam");
        builder.write_bam(&input).expect("failed to write BAM");

        let err = RawExternalSorter::new(SortOrder::Coordinate)
            .spill_codec(crate::codec::SpillCodec::Zstd)
            .temp_compression(0)
            .sort(&input, &output)
            .expect_err("sort should reject temp_compression=0 + zstd");
        let msg = err.to_string();
        assert!(
            msg.contains("temp_compression=0 is only supported with SpillCodec::Bgzf"),
            "unexpected error: {msg}"
        );
        assert!(!output.exists(), "no output should be produced on validation failure");
    }

    #[test]
    fn test_raw_sorter_max_temp_files() {
        let sorter = RawExternalSorter::new(SortOrder::Coordinate).max_temp_files(0);
        assert_eq!(sorter.max_temp_files, 0);
    }

    // ========================================================================
    // create_output_header tests
    // ========================================================================

    #[test]
    fn test_create_output_header_coordinate() {
        let sorter = RawExternalSorter::new(SortOrder::Coordinate);
        let header = Header::builder().build();
        let output_header = sorter.create_output_header(&header);

        let hd = output_header.header().expect("header should have HD record");
        let so = hd.other_fields().get(b"SO").expect("should have SO tag");
        assert_eq!(<_ as AsRef<[u8]>>::as_ref(so), b"coordinate");
    }

    #[test]
    fn test_create_output_header_queryname() {
        let sorter = RawExternalSorter::new(SortOrder::Queryname(QuerynameComparator::default()));
        let header = Header::builder().build();
        let output_header = sorter.create_output_header(&header);

        let hd = output_header.header().expect("header should have HD record");
        let so = hd.other_fields().get(b"SO").expect("should have SO tag");
        assert_eq!(<_ as AsRef<[u8]>>::as_ref(so), b"queryname");
    }

    #[test]
    fn test_create_output_header_template_coordinate() {
        let sorter = RawExternalSorter::new(SortOrder::TemplateCoordinate);
        let header = Header::builder().build();
        let output_header = sorter.create_output_header(&header);

        let hd = output_header.header().expect("header should have HD record");
        let fields = hd.other_fields();

        let so = fields.get(b"SO").expect("should have SO tag");
        assert_eq!(<_ as AsRef<[u8]>>::as_ref(so), b"unsorted");

        let go = fields.get(b"GO").expect("should have GO tag");
        assert_eq!(<_ as AsRef<[u8]>>::as_ref(go), b"query");

        let ss = fields.get(b"SS").expect("should have SS tag");
        assert_eq!(<_ as AsRef<[u8]>>::as_ref(ss), b"template-coordinate");
    }

    // ========================================================================
    // find_string_tag_in_record tests (via fgumi_raw_bam)
    // ========================================================================

    /// Helper to build minimal BAM bytes with aux data appended.
    /// Fixed 32-byte header + read name "rea\0" (`l_read_name=4`) + no cigar + no seq + aux.
    fn build_bam_with_aux(aux_data: &[u8]) -> Vec<u8> {
        let l_read_name: u8 = 4; // "rea" + null
        let mut bam = vec![0u8; 32];
        bam[8] = l_read_name;
        // n_cigar_op = 0 (bytes 12-13 already zero)
        // l_seq = 0 (bytes 16-19 already zero)
        // read name
        bam.extend_from_slice(b"rea\0");
        // no cigar, no seq, no qual
        // append aux data
        bam.extend_from_slice(aux_data);
        bam
    }

    #[rstest::rstest]
    #[case::present(b"RGZgroup1\0".as_slice(), Some(b"group1".as_slice()))]
    #[case::absent(b"".as_slice(), None)]
    fn test_find_rg_tag(#[case] aux_data: &[u8], #[case] expected: Option<&[u8]>) {
        let bam = build_bam_with_aux(aux_data);
        assert_eq!(
            fgumi_raw_bam::RawRecordView::new(&bam).tags().find_string(SamTag::RG),
            expected
        );
    }

    #[test]
    fn test_find_rg_tag_after_other_tags() {
        // XY:i:42 followed by RG:Z:mygroup — aux built dynamically for the integer tag
        let mut aux = Vec::new();
        aux.extend_from_slice(b"XYi");
        aux.extend_from_slice(&42i32.to_le_bytes());
        aux.extend_from_slice(b"RGZmygroup\0");
        let bam = build_bam_with_aux(&aux);
        assert_eq!(
            fgumi_raw_bam::RawRecordView::new(&bam).tags().find_string(SamTag::RG),
            Some(b"mygroup".as_slice())
        );
    }

    // ========================================================================
    // RawExternalSorter cell_tag builder tests
    // ========================================================================

    #[test]
    fn test_raw_sorter_cell_tag_default_is_none() {
        let sorter = RawExternalSorter::new(SortOrder::TemplateCoordinate);
        assert!(sorter.cell_tag.is_none());
    }

    #[test]
    fn test_raw_sorter_cell_tag_builder() {
        let sorter = RawExternalSorter::new(SortOrder::TemplateCoordinate).cell_tag(SamTag::CB);
        assert_eq!(sorter.cell_tag, Some(SamTag::CB));
    }

    /// Smoke test: an auto-detected lite (no-CB, no-tertiary) template-coordinate
    /// sort runs end-to-end and preserves the record count. With the default
    /// `KeyTypesSpec::Auto`, records lacking CB/MI/multi-library auxiliary data
    /// select the narrow `TemplateKey24` lane key; this exercises that path
    /// through Phase 1 and Phase 2 and confirms output completeness.
    #[test]
    fn template_sort_auto_lite_roundtrips_record_count() {
        let dir = tempfile::tempdir().unwrap();
        let (sorted, names) =
            create_sorted_bam(dir.path(), "lite", 50, 0, SortOrder::TemplateCoordinate);
        assert_eq!(collect_read_names(&sorted).len(), names.len() * 2);
    }

    /// Regression: a 0-record template-coordinate sort must still write a valid
    /// header-only output BAM (matching the coordinate and queryname orders),
    /// not skip the output entirely. The empty-input path skips the pre-loop
    /// first-record push, runs the read loop zero times, and falls through to
    /// the `chunk_files.is_empty()` fast path that writes the header-only BAM.
    #[test]
    fn template_sort_empty_input_writes_header_only_bam() {
        let dir = tempfile::tempdir().unwrap();
        let (sorted, names) =
            create_sorted_bam(dir.path(), "empty", 0, 0, SortOrder::TemplateCoordinate);
        assert!(names.is_empty());
        assert!(sorted.exists(), "empty-input sort must still produce an output file");
        assert_eq!(count_bam_records(&sorted), 0, "output must contain only the header");
    }

    // ========================================================================
    // extract_template_key_inline cell_tag tests
    // ========================================================================

    fn test_cb_hasher() -> ahash::RandomState {
        cb_hasher()
    }

    /// Build minimal BAM bytes with mapped read at (tid, pos) on the forward strand,
    /// with optional aux data appended.
    #[allow(clippy::cast_possible_truncation)]
    fn build_mapped_bam(tid: i32, pos: i32, name: &[u8], aux: &[u8]) -> Vec<u8> {
        let l_read_name = (name.len() + 1) as u8; // +1 for null terminator
        let mut bam = vec![0u8; 32];
        // ref_id at offset 0..4
        bam[0..4].copy_from_slice(&tid.to_le_bytes());
        // pos at offset 4..8
        bam[4..8].copy_from_slice(&pos.to_le_bytes());
        // l_read_name at offset 8
        bam[8] = l_read_name;
        // flags at offset 14..16: paired + proper pair = 0x03
        bam[14..16].copy_from_slice(&3u16.to_le_bytes());
        // mate_ref_id at offset 20..24: same tid
        bam[20..24].copy_from_slice(&tid.to_le_bytes());
        // mate_pos at offset 24..28: same pos
        bam[24..28].copy_from_slice(&pos.to_le_bytes());
        // read name
        bam.extend_from_slice(name);
        bam.push(0); // null terminator
        // no cigar, no seq, no qual
        // aux data
        bam.extend_from_slice(aux);
        bam
    }

    /// Build CB:Z: aux tag bytes.
    fn cb_aux(value: &[u8]) -> Vec<u8> {
        let mut aux = Vec::new();
        aux.extend_from_slice(b"CBZ");
        aux.extend_from_slice(value);
        aux.push(0); // null terminator
        aux
    }

    #[test]
    fn test_extract_template_key_cb_present_has_nonzero_hash() {
        let header = Header::builder().build();
        let lib_lookup = LibraryLookup::from_header(&header);
        let aux = cb_aux(b"ACGTACGT");
        let bam = build_mapped_bam(0, 100, b"read1", &aux);

        let key =
            extract_template_key_inline(&bam, &lib_lookup, Some(SamTag::CB), &test_cb_hasher());
        assert_ne!(key.cb_hash, 0, "CB present should produce non-zero cb_hash");
    }

    #[test]
    fn test_extract_template_key_cb_absent_has_zero_hash() {
        let header = Header::builder().build();
        let lib_lookup = LibraryLookup::from_header(&header);
        // No CB tag in aux data
        let bam = build_mapped_bam(0, 100, b"read1", &[]);

        let key =
            extract_template_key_inline(&bam, &lib_lookup, Some(SamTag::CB), &test_cb_hasher());
        assert_eq!(key.cb_hash, 0, "missing CB tag should produce cb_hash=0");
    }

    #[test]
    fn test_extract_template_key_cell_tag_none_has_zero_hash() {
        let header = Header::builder().build();
        let lib_lookup = LibraryLookup::from_header(&header);
        let aux = cb_aux(b"ACGTACGT");
        let bam = build_mapped_bam(0, 100, b"read1", &aux);

        let key = extract_template_key_inline(&bam, &lib_lookup, None, &test_cb_hasher());
        assert_eq!(key.cb_hash, 0, "cell_tag=None should produce cb_hash=0");
    }

    #[test]
    fn test_extract_template_key_different_cb_values_differ() {
        let header = Header::builder().build();
        let lib_lookup = LibraryLookup::from_header(&header);

        let aux1 = cb_aux(b"ACGTACGT");
        let bam1 = build_mapped_bam(0, 100, b"read1", &aux1);
        let key1 =
            extract_template_key_inline(&bam1, &lib_lookup, Some(SamTag::CB), &test_cb_hasher());

        let aux2 = cb_aux(b"TGCATGCA");
        let bam2 = build_mapped_bam(0, 100, b"read1", &aux2);
        let key2 =
            extract_template_key_inline(&bam2, &lib_lookup, Some(SamTag::CB), &test_cb_hasher());

        assert_ne!(
            key1.cb_hash, key2.cb_hash,
            "different CB values should produce different hashes"
        );
    }

    #[test]
    fn test_extract_template_key_cb_hash_is_deterministic() {
        let header = Header::builder().build();
        let lib_lookup = LibraryLookup::from_header(&header);
        let aux = cb_aux(b"ACGTACGT");
        let bam = build_mapped_bam(0, 100, b"read1", &aux);

        let key1 =
            extract_template_key_inline(&bam, &lib_lookup, Some(SamTag::CB), &test_cb_hasher());
        let key2 =
            extract_template_key_inline(&bam, &lib_lookup, Some(SamTag::CB), &test_cb_hasher());
        assert_eq!(key1.cb_hash, key2.cb_hash, "same input should produce same cb_hash");
    }

    #[test]
    fn test_extract_template_key_unmapped_with_cb() {
        let header = Header::builder().build();
        let lib_lookup = LibraryLookup::from_header(&header);
        let aux = cb_aux(b"ACGTACGT");

        // Build unmapped read (flag = PAIRED | UNMAPPED | MATE_UNMAPPED = 0x0D)
        let mut bam = vec![0u8; 32];
        bam[8] = 6; // l_read_name = 6 ("read1" + null)
        bam[14..16].copy_from_slice(&0x000Du16.to_le_bytes()); // flags
        // ref_id = -1 (unmapped)
        bam[0..4].copy_from_slice(&(-1i32).to_le_bytes());
        bam[4..8].copy_from_slice(&(-1i32).to_le_bytes()); // pos = -1
        bam[20..24].copy_from_slice(&(-1i32).to_le_bytes()); // mate_ref_id = -1
        bam[24..28].copy_from_slice(&(-1i32).to_le_bytes()); // mate_pos = -1
        bam.extend_from_slice(b"read1\0");
        bam.extend_from_slice(&aux);

        let key =
            extract_template_key_inline(&bam, &lib_lookup, Some(SamTag::CB), &test_cb_hasher());
        assert_ne!(key.cb_hash, 0, "unmapped read with CB should have non-zero cb_hash");
        assert_eq!(key.primary, u64::MAX, "unmapped both-mates should have MAX primary");
    }

    // ========================================================================
    // Consolidation chunk naming tests
    // ========================================================================

    /// Count records in a BAM file by reading with the raw BAM reader.
    fn count_bam_records(path: &std::path::Path) -> u64 {
        use crate::read_ahead::RawReadAheadReader;
        let (reader, _) = create_raw_bam_reader(path, 1).expect("failed to create raw BAM reader");
        RawReadAheadReader::new(reader).count() as u64
    }

    /// Read every record body of a BAM file into a `Vec<Vec<u8>>` for exact
    /// output comparison.
    fn read_all_bam_records(path: &std::path::Path) -> Vec<Vec<u8>> {
        use crate::read_ahead::RawReadAheadReader;
        let (reader, _) = create_raw_bam_reader(path, 1).expect("failed to create raw BAM reader");
        RawReadAheadReader::new(reader).map(|r| r.as_ref().to_vec()).collect()
    }

    /// Verifies that sort with consolidation preserves all records.
    ///
    /// Uses a tiny memory limit and low `max_temp_files` to force many chunks and
    /// consolidation. Before the fix (chunk files named by `chunk_files.len()`),
    /// consolidation would drain entries from the vector, shrinking its length,
    /// causing new chunks to collide with existing non-consolidated chunk files.
    ///
    /// Parameterized over both spill codecs so the consolidation merge writer
    /// path is exercised for each format. zstd uses `temp_compression(1)` since
    /// zstd has no level-0 "stored" mode (the CLI rejects 0+zstd; here we make
    /// the same choice explicitly so the test matches user-facing behavior).
    #[rstest::rstest]
    #[case::coordinate_bgzf(SortOrder::Coordinate, crate::codec::SpillCodec::Bgzf, 0)]
    #[case::queryname_bgzf(
        SortOrder::Queryname(QuerynameComparator::default()),
        crate::codec::SpillCodec::Bgzf,
        0
    )]
    #[case::queryname_natural_bgzf(
        SortOrder::Queryname(QuerynameComparator::Natural),
        crate::codec::SpillCodec::Bgzf,
        0
    )]
    #[case::template_coordinate_bgzf(
        SortOrder::TemplateCoordinate,
        crate::codec::SpillCodec::Bgzf,
        0
    )]
    #[case::coordinate_zstd(SortOrder::Coordinate, crate::codec::SpillCodec::Zstd, 1)]
    #[case::queryname_zstd(
        SortOrder::Queryname(QuerynameComparator::default()),
        crate::codec::SpillCodec::Zstd,
        1
    )]
    #[case::queryname_natural_zstd(
        SortOrder::Queryname(QuerynameComparator::Natural),
        crate::codec::SpillCodec::Zstd,
        1
    )]
    #[case::template_coordinate_zstd(
        SortOrder::TemplateCoordinate,
        crate::codec::SpillCodec::Zstd,
        1
    )]
    fn test_sort_with_consolidation_preserves_all_records(
        #[case] sort_order: SortOrder,
        #[case] spill_codec: crate::codec::SpillCodec,
        #[case] temp_compression: u32,
    ) {
        use fgumi_sam::SamBuilder;

        let num_pairs = 30;
        let mut builder = SamBuilder::new();
        for i in 0..num_pairs {
            let _ = builder
                .add_pair()
                .name(&format!("read{i}"))
                .start1(i * 200 + 1)
                .start2(i * 200 + 101)
                .build();
        }

        let dir = tempfile::tempdir().expect("failed to create temp directory");
        let input = dir.path().join("input.bam");
        let output = dir.path().join("output.bam");
        builder.write_bam(&input).expect("failed to write BAM");

        // Tiny memory limit forces many chunks; low max_temp_files triggers consolidation
        let stats = RawExternalSorter::new(sort_order)
            .memory_limit(1024) // 1 KB — each chunk holds very few records
            .max_temp_files(4)
            .spill_codec(spill_codec)
            .temp_compression(temp_compression)
            .output_compression(0)
            .sort(&input, &output)
            .expect("sort should succeed");

        assert!(
            stats.chunks_written >= 5,
            "expected at least 5 chunks to exercise post-consolidation naming, got {}",
            stats.chunks_written
        );

        // Count records in the output BAM to verify no data was lost
        let expected = (num_pairs * 2) as u64;
        let observed = count_bam_records(&output);
        assert_eq!(observed, expected, "chunk filename collision likely lost data");
    }

    /// Verifies that sort with many chunks exercises the pool-integrated
    /// merge readers during the final merge phase (not just consolidation).
    #[rstest::rstest]
    #[case::coordinate(SortOrder::Coordinate)]
    #[case::queryname(SortOrder::Queryname(QuerynameComparator::default()))]
    #[case::queryname_natural(SortOrder::Queryname(QuerynameComparator::Natural))]
    #[case::template_coordinate(SortOrder::TemplateCoordinate)]
    fn test_sort_many_chunks_with_semaphore(#[case] sort_order: SortOrder) {
        use fgumi_sam::SamBuilder;

        let num_pairs = 200;
        let mut builder = SamBuilder::new();
        for i in 0..num_pairs {
            let _ = builder
                .add_pair()
                .name(&format!("read{i}"))
                .start1(i * 200 + 1)
                .start2(i * 200 + 101)
                .build();
        }

        let dir = tempfile::tempdir().expect("failed to create temp directory");
        let input = dir.path().join("input.bam");
        let output = dir.path().join("output.bam");
        builder.write_bam(&input).expect("failed to write BAM");

        // Small memory limit + many records = guaranteed spill to multiple chunks
        // (no consolidation, so the semaphore must cap concurrent readers).
        // 32 KiB is small enough to force several spills across 400 records but
        // large enough to avoid creating hundreds of tiny chunks that saturate
        // OS threads on the CI runner and cause timeouts.
        let stats = RawExternalSorter::new(sort_order)
            .memory_limit(32 * 1024)
            .max_temp_files(0) // disable consolidation
            .threads(2) // semaphore allows 2 concurrent readers
            .spill_codec(crate::codec::SpillCodec::Bgzf) // level 0 = stored mode (zstd has no level 0)
            .temp_compression(0)
            .output_compression(0)
            .sort(&input, &output)
            .expect("sort should succeed");

        assert!(
            stats.chunks_written >= 2,
            "expected multiple chunks to exercise merge, got {}",
            stats.chunks_written
        );

        let expected = (num_pairs * 2) as u64;
        let observed = count_bam_records(&output);
        assert_eq!(observed, expected, "semaphore-capped merge lost data");
    }

    /// S3-008: the zstd Phase-2 decompress path recycles block buffers through
    /// the shared pool (worker checks out, merge consumer checks back in).
    /// Force many zstd spill chunks so the pool is cycled well past its
    /// capacity, and assert the recycled-buffer merge output is byte-identical
    /// to an in-memory sort of the same input (recycling must not corrupt,
    /// reorder, drop, or duplicate records).
    #[rstest::rstest]
    #[case::coordinate(SortOrder::Coordinate)]
    #[case::queryname_natural(SortOrder::Queryname(QuerynameComparator::Natural))]
    #[case::template_coordinate(SortOrder::TemplateCoordinate)]
    fn test_zstd_phase2_buffer_recycling_preserves_output(#[case] sort_order: SortOrder) {
        use fgumi_sam::SamBuilder;

        let num_pairs = 300;
        let mut builder = SamBuilder::new();
        for i in 0..num_pairs {
            let _ = builder
                .add_pair()
                .name(&format!("read{i}"))
                .start1((i * 7 % 4000) + 1)
                .start2((i * 7 % 4000) + 101)
                .build();
        }

        let dir = tempfile::tempdir().expect("failed to create temp directory");
        let input = dir.path().join("input.bam");
        builder.write_bam(&input).expect("failed to write BAM");

        // Spilled zstd sort: a tiny memory limit forces many chunks → many
        // Phase-2 zstd decompress blocks cycled through the recycling pool.
        let spilled = dir.path().join("spilled.bam");
        let stats = RawExternalSorter::new(sort_order)
            .memory_limit(16 * 1024)
            .threads(2)
            .spill_codec(crate::codec::SpillCodec::Zstd)
            .temp_compression(3)
            .output_compression(0)
            .sort(&input, &spilled)
            .expect("spilled zstd sort should succeed");
        assert!(
            stats.chunks_written >= 3,
            "expected several zstd spill chunks to cycle the buffer pool, got {}",
            stats.chunks_written
        );

        // In-memory reference sort (no spill, no Phase-2 decompression).
        let in_memory = dir.path().join("in_memory.bam");
        RawExternalSorter::new(sort_order)
            .memory_limit(256 * 1024 * 1024)
            .threads(2)
            .output_compression(0)
            .sort(&input, &in_memory)
            .expect("in-memory sort should succeed");

        assert_eq!(
            read_all_bam_records(&spilled),
            read_all_bam_records(&in_memory),
            "zstd-spilled recycled-buffer merge must match the in-memory sort exactly",
        );
    }

    // ========================================================================
    // Sub-array merge source tests
    // ========================================================================

    /// Verifies that multi-threaded sort produces the same output as single-threaded
    /// sort for each sort order.
    #[rstest::rstest]
    #[case::coordinate(SortOrder::Coordinate)]
    #[case::queryname(SortOrder::Queryname(QuerynameComparator::default()))]
    #[case::queryname_natural(SortOrder::Queryname(QuerynameComparator::Natural))]
    #[case::template_coordinate(SortOrder::TemplateCoordinate)]
    fn test_sort_sub_arrays_match_single_thread(#[case] sort_order: SortOrder) {
        use fgumi_sam::SamBuilder;

        let num_pairs = 50;
        let mut builder = SamBuilder::new();
        for i in 0..num_pairs {
            let _ = builder
                .add_pair()
                .name(&format!("read{i}"))
                .start1(i * 200 + 1)
                .start2(i * 200 + 101)
                .build();
        }

        let dir = tempfile::tempdir().expect("failed to create temp directory");
        let input = dir.path().join("input.bam");
        let output_st = dir.path().join("output_1t.bam");
        let output_mt = dir.path().join("output_2t.bam");
        builder.write_bam(&input).expect("failed to write BAM");

        // Sort single-threaded
        RawExternalSorter::new(sort_order)
            .memory_limit(16 * 1024) // force spill (50 pairs × ~300B ≈ 15KB) so merge path is exercised
            .threads(1)
            .spill_codec(crate::codec::SpillCodec::Bgzf) // level 0 = stored mode (zstd has no level 0)
            .temp_compression(0)
            .output_compression(0)
            .sort(&input, &output_st)
            .expect("sort should succeed");

        // Sort multi-threaded
        RawExternalSorter::new(sort_order)
            .memory_limit(16 * 1024)
            .threads(2)
            .spill_codec(crate::codec::SpillCodec::Bgzf)
            .temp_compression(0)
            .output_compression(0)
            .sort(&input, &output_mt)
            .expect("sort should succeed");

        let names_st = collect_read_names(&output_st);
        let names_mt = collect_read_names(&output_mt);

        let expected = num_pairs * 2;
        assert_eq!(names_st.len(), expected, "single-thread record count mismatch");
        assert_eq!(names_mt.len(), expected, "multi-thread record count mismatch");

        // Both outputs must have the same read names in the same order
        assert_eq!(names_st, names_mt, "multi-thread sort order differs from single-thread");
    }

    /// Verifies that the in-memory-only path (no spill to disk) works correctly.
    #[rstest::rstest]
    #[case::coordinate(SortOrder::Coordinate)]
    #[case::queryname(SortOrder::Queryname(QuerynameComparator::default()))]
    #[case::queryname_natural(SortOrder::Queryname(QuerynameComparator::Natural))]
    #[case::template_coordinate(SortOrder::TemplateCoordinate)]
    fn test_sort_sub_arrays_in_memory_only(#[case] sort_order: SortOrder) {
        use fgumi_sam::SamBuilder;

        let num_pairs = 20;
        let mut builder = SamBuilder::new();
        for i in 0..num_pairs {
            let _ = builder
                .add_pair()
                .name(&format!("read{i}"))
                .start1(i * 200 + 1)
                .start2(i * 200 + 101)
                .build();
        }

        let dir = tempfile::tempdir().expect("failed to create temp directory");
        let input = dir.path().join("input.bam");
        let output = dir.path().join("output.bam");
        builder.write_bam(&input).expect("failed to write BAM");

        // Large memory limit so everything stays in memory (no spill chunks).
        RawExternalSorter::new(sort_order)
            .memory_limit(10 * 1024 * 1024)
            .threads(2)
            .output_compression(0)
            .sort(&input, &output)
            .expect("sort should succeed");

        let expected = (num_pairs * 2) as u64;
        let observed = count_bam_records(&output);
        assert_eq!(observed, expected, "sort lost data");
    }

    // ========================================================================
    // merge_bams tests
    // ========================================================================

    /// Helper: create a `SamBuilder` with `num_pairs` pairs at non-overlapping positions,
    /// write an unsorted BAM, sort it with the given order, and return the sorted path.
    /// `start_offset` shifts all positions so different inputs have distinct read names/positions.
    fn create_sorted_bam(
        dir: &Path,
        prefix: &str,
        num_pairs: usize,
        start_offset: usize,
        sort_order: SortOrder,
    ) -> (PathBuf, Vec<String>) {
        let mut builder = SamBuilder::new();
        let mut names = Vec::with_capacity(num_pairs);
        for i in 0..num_pairs {
            let name = format!("{prefix}_read{i:04}");
            names.push(name.clone());
            let _ = builder
                .add_pair()
                .name(&name)
                .start1((start_offset + i * 200) + 1)
                .start2((start_offset + i * 200) + 101)
                .build();
        }
        let unsorted = dir.join(format!("{prefix}_unsorted.bam"));
        let sorted = dir.join(format!("{prefix}_sorted.bam"));
        builder.write_bam(&unsorted).expect("failed to write BAM");
        RawExternalSorter::new(sort_order)
            .output_compression(0)
            .sort(&unsorted, &sorted)
            .expect("sort should succeed");
        (sorted, names)
    }

    /// Collect all read names from a BAM file as strings.
    fn collect_read_names(path: &Path) -> Vec<String> {
        use crate::read_ahead::RawReadAheadReader;
        let (reader, _) = create_raw_bam_reader(path, 1).expect("failed to create raw BAM reader");
        RawReadAheadReader::new(reader)
            .map(|rec| {
                let name_bytes = fgumi_raw_bam::RawRecordView::new(rec.as_ref()).read_name();
                String::from_utf8(name_bytes.to_vec()).expect("read name should be valid UTF-8")
            })
            .collect()
    }

    /// Collect (`ref_id`, pos) tuples for every record in a BAM.
    fn collect_positions(path: &Path) -> Vec<(i32, i32)> {
        use crate::read_ahead::RawReadAheadReader;
        let (reader, _) = create_raw_bam_reader(path, 1).expect("failed to create raw BAM reader");
        RawReadAheadReader::new(reader)
            .map(|rec| {
                let bytes = rec.as_ref();
                {
                    let v = fgumi_raw_bam::RawRecordView::new(bytes);
                    (v.ref_id(), v.pos())
                }
            })
            .collect()
    }

    /// Helper to build a merge header from the `SamBuilder` default header.
    fn default_merge_header() -> Header {
        SamBuilder::new().header.clone()
    }

    #[test]
    fn test_merge_bams_coordinate_sort() {
        let dir = tempfile::tempdir().expect("failed to create temp directory");
        let (bam_a, _) = create_sorted_bam(dir.path(), "a", 10, 0, SortOrder::Coordinate);
        let (bam_b, _) = create_sorted_bam(dir.path(), "b", 10, 10_000, SortOrder::Coordinate);

        let merged = dir.path().join("merged.bam");
        let header = default_merge_header();
        let count = RawExternalSorter::new(SortOrder::Coordinate)
            .output_compression(0)
            .merge_bams(&[bam_a, bam_b], &header, &merged)
            .expect("sort should succeed");

        // 10 pairs * 2 records * 2 inputs = 40
        assert_eq!(count, 40);
        assert_eq!(count_bam_records(&merged), 40);

        // Verify coordinate sort order: (ref_id, pos) is non-decreasing
        let positions = collect_positions(&merged);
        for w in positions.windows(2) {
            assert!(w[0] <= w[1], "coordinate sort violated: {:?} > {:?}", w[0], w[1]);
        }
    }

    #[test]
    fn test_merge_bams_template_coordinate_sort() {
        let dir = tempfile::tempdir().expect("failed to create temp directory");
        let (bam_a, _) = create_sorted_bam(dir.path(), "a", 10, 0, SortOrder::TemplateCoordinate);
        let (bam_b, _) =
            create_sorted_bam(dir.path(), "b", 10, 10_000, SortOrder::TemplateCoordinate);

        let merged = dir.path().join("merged.bam");
        let header = default_merge_header();
        let count = RawExternalSorter::new(SortOrder::TemplateCoordinate)
            .output_compression(0)
            .merge_bams(&[bam_a, bam_b], &header, &merged)
            .expect("sort should succeed");

        assert_eq!(count, 40);
        assert_eq!(count_bam_records(&merged), 40);

        // Verify template-coordinate order/identity via an independent oracle
        // (mirrors `narrow_sort_output_passes_full_width_verify`): re-extract the
        // template key from every merged record and assert non-decreasing core
        // order. A count-only check would pass even if the merge reordered or
        // duplicated records.
        let (_, hdr) = create_raw_bam_reader(&merged, 1).expect("header");
        let lib = LibraryLookup::from_header(&hdr);
        let hasher = cb_hasher();
        let file = std::fs::File::open(&merged).expect("open merged");
        let mut raw_reader = crate::reader::RawBamRecordReader::new(file).expect("reader");
        raw_reader.skip_header().expect("skip header");
        let (total, violations, first) = crate::verify::verify_sort_order(
            raw_reader,
            |bam| extract_template_key_inline(bam, &lib, None, &hasher),
            |cur: &TemplateKey40, prev: &TemplateKey40| {
                cur.core_cmp(prev) == std::cmp::Ordering::Less
            },
        )
        .expect("verify runs");
        assert_eq!(total, 40, "verify must see every merged record");
        assert_eq!(
            violations, 0,
            "template-coordinate order violated after merge (first={first:?}, total={total})"
        );
    }

    #[test]
    fn test_merge_bams_queryname_sort() {
        let dir = tempfile::tempdir().expect("failed to create temp directory");
        let (bam_a, _) = create_sorted_bam(
            dir.path(),
            "a",
            10,
            0,
            SortOrder::Queryname(QuerynameComparator::default()),
        );
        let (bam_b, _) = create_sorted_bam(
            dir.path(),
            "b",
            10,
            10_000,
            SortOrder::Queryname(QuerynameComparator::default()),
        );

        let merged = dir.path().join("merged.bam");
        let header = default_merge_header();
        let count = RawExternalSorter::new(SortOrder::Queryname(QuerynameComparator::default()))
            .output_compression(0)
            .merge_bams(&[bam_a, bam_b], &header, &merged)
            .expect("sort should succeed");

        assert_eq!(count, 40);
        assert_eq!(count_bam_records(&merged), 40);

        // Verify queryname sort order: read names are non-decreasing
        let names = collect_read_names(&merged);
        for w in names.windows(2) {
            assert!(w[0] <= w[1], "queryname sort violated: {:?} > {:?}", w[0], w[1]);
        }
    }

    #[test]
    fn test_merge_bams_single_input() {
        let dir = tempfile::tempdir().expect("failed to create temp directory");
        let (bam_a, _) = create_sorted_bam(dir.path(), "a", 15, 0, SortOrder::Coordinate);

        let merged = dir.path().join("merged.bam");
        let header = default_merge_header();
        let count = RawExternalSorter::new(SortOrder::Coordinate)
            .output_compression(0)
            .merge_bams(&[bam_a], &header, &merged)
            .expect("sort should succeed");

        // 15 pairs * 2 = 30
        assert_eq!(count, 30);
        assert_eq!(count_bam_records(&merged), 30);
    }

    #[test]
    fn test_merge_bams_preserves_all_records() {
        let dir = tempfile::tempdir().expect("failed to create temp directory");
        let (bam_a, names_a) = create_sorted_bam(
            dir.path(),
            "a",
            5,
            0,
            SortOrder::Queryname(QuerynameComparator::default()),
        );
        let (bam_b, names_b) = create_sorted_bam(
            dir.path(),
            "b",
            5,
            10_000,
            SortOrder::Queryname(QuerynameComparator::default()),
        );

        let merged = dir.path().join("merged.bam");
        let header = default_merge_header();
        RawExternalSorter::new(SortOrder::Queryname(QuerynameComparator::default()))
            .output_compression(0)
            .merge_bams(&[bam_a, bam_b], &header, &merged)
            .expect("sort should succeed");

        let merged_names: std::collections::HashSet<String> =
            collect_read_names(&merged).into_iter().collect();

        // Every expected name (from both inputs) must appear in the merged output
        for name in names_a.iter().chain(names_b.iter()) {
            assert!(merged_names.contains(name), "read name {name:?} missing from merged output");
        }
    }

    #[test]
    fn test_merge_bams_many_inputs() {
        let dir = tempfile::tempdir().expect("failed to create temp directory");
        let k = 8;
        let pairs_per_input = 5;
        let mut inputs = Vec::with_capacity(k);

        for i in 0..k {
            let (bam, _) = create_sorted_bam(
                dir.path(),
                &format!("in{i}"),
                pairs_per_input,
                i * 50_000,
                SortOrder::Coordinate,
            );
            inputs.push(bam);
        }

        let merged = dir.path().join("merged.bam");
        let header = default_merge_header();
        let count = RawExternalSorter::new(SortOrder::Coordinate)
            .output_compression(0)
            .merge_bams(&inputs, &header, &merged)
            .expect("sort should succeed");

        let expected = (k * pairs_per_input * 2) as u64; // 8 * 5 * 2 = 80
        assert_eq!(count, expected);
        assert_eq!(count_bam_records(&merged), expected);

        // Verify coordinate sort order
        let positions = collect_positions(&merged);
        for w in positions.windows(2) {
            assert!(w[0] <= w[1], "coordinate sort violated with k={k}: {:?} > {:?}", w[0], w[1]);
        }
    }

    #[test]
    fn test_merge_bams_queryname_natural_sort() {
        let dir = tempfile::tempdir().expect("failed to create temp directory");
        let nat = SortOrder::Queryname(QuerynameComparator::Natural);
        let (bam_a, _) = create_sorted_bam(dir.path(), "a", 10, 0, nat);
        let (bam_b, _) = create_sorted_bam(dir.path(), "b", 10, 10_000, nat);

        let merged = dir.path().join("merged.bam");
        let header = default_merge_header();
        let count = RawExternalSorter::new(nat)
            .output_compression(0)
            .merge_bams(&[bam_a, bam_b], &header, &merged)
            .expect("merge should succeed");

        assert_eq!(count, 40);
        assert_eq!(count_bam_records(&merged), 40);

        // Verify natural queryname order: names are non-decreasing
        let names = collect_read_names(&merged);
        for w in names.windows(2) {
            assert!(w[0] <= w[1], "natural queryname sort violated: {:?} > {:?}", w[0], w[1]);
        }
    }

    // ========================================================================
    // SortPhaseTimer tests
    // ========================================================================

    #[test]
    fn test_sort_phase_timer_all_methods() {
        let mut timer = SortPhaseTimer::new();
        assert!(timer.overall_start.is_some());
        assert!(timer.read_span_start.is_some());

        // end_read_span accumulates and clears the span
        let elapsed = timer.end_read_span();
        assert!(timer.read_secs >= 0.0);
        assert!(timer.read_span_start.is_none());
        let _ = elapsed;

        // end_read_span with no active span → zero, no panic
        let elapsed2 = timer.end_read_span();
        assert_eq!(elapsed2, std::time::Duration::ZERO);
        assert!(timer.read_secs >= 0.0); // unchanged

        // begin_read_span restores the span
        timer.begin_read_span();
        assert!(timer.read_span_start.is_some());

        // time_sort measures elapsed
        timer.time_sort(|| {});
        assert!(timer.sort_secs >= 0.0);

        // time_spill_write: counts spills and returns the closure result
        let result = timer.time_spill_write(|| Ok::<u32, anyhow::Error>(42));
        assert_eq!(result.unwrap(), 42);
        assert_eq!(timer.spill_count, 1);
        assert!(timer.spill_write_secs >= 0.0);

        // record_spill_size: adds file size to total
        let dir = tempfile::tempdir().expect("tempdir");
        let path = dir.path().join("spill.bin");
        std::fs::write(&path, b"hello world").expect("write");
        timer.record_spill_size(&path);
        assert_eq!(timer.total_spill_bytes, 11);

        // record_spill_size with nonexistent path: no-op, no panic
        timer.record_spill_size(&dir.path().join("nonexistent.bin"));
        assert_eq!(timer.total_spill_bytes, 11);

        // time_consolidate: fast closure (<1ms) → consolidate_count stays 0
        timer.time_consolidate(|| Ok(())).expect("consolidate ok");
        assert_eq!(timer.consolidate_count, 0);

        // time_consolidate: slow closure (>1ms) → counted
        timer
            .time_consolidate(|| {
                std::thread::sleep(std::time::Duration::from_millis(10));
                Ok(())
            })
            .expect("consolidate ok");
        assert_eq!(timer.consolidate_count, 1);
        assert!(timer.consolidate_secs > 0.0);

        // time_merge
        timer.time_merge(|| Ok::<(), anyhow::Error>(())).expect("merge ok");
        assert!(timer.merge_secs >= 0.0);

        // time_write_output
        timer.time_write_output(|| Ok(())).expect("write ok");
        assert!(timer.write_output_secs >= 0.0);

        // log_summary must not panic (output goes to log sink)
        timer.log_summary(4);
    }

    // ========================================================================
    // Chunk boundary alignment tests
    // ========================================================================

    /// Regression test: `split_off` from the tail must align with `par_chunks_mut`
    /// boundaries so that every emitted chunk is individually sorted.
    ///
    /// The old code split fixed-size chunks from the tail, which crossed sorted
    /// chunk boundaries when `n % chunk_size != 0`.
    #[test]
    #[allow(clippy::cast_possible_truncation)]
    fn test_split_off_aligns_with_par_chunks_mut() {
        use rayon::prelude::*;

        let pool =
            rayon::ThreadPoolBuilder::new().num_threads(4).build().expect("build rayon pool");

        // Test a variety of (n, threads) combinations that exercise the tail case
        for &(n, threads) in &[(10, 3), (11, 4), (17, 3), (100, 7), (1000, 8), (13, 5)] {
            let mut entries: Vec<(u64, Vec<u8>)> =
                (0..n).rev().map(|i| (i as u64, vec![i as u8])).collect();

            let chunk_size = entries.len().div_ceil(threads);

            pool.install(|| {
                entries.par_chunks_mut(chunk_size).for_each(|chunk| {
                    chunk.sort_unstable_by(|a, b| a.0.cmp(&b.0));
                });
            });

            // Apply the same split logic as production code
            let mut remaining = std::mem::take(&mut entries);
            let num_chunks = remaining.len().div_ceil(chunk_size);
            let mut chunks: Vec<Vec<(u64, Vec<u8>)>> = Vec::with_capacity(num_chunks);
            let tail_len = remaining.len() % chunk_size;
            if tail_len != 0 {
                let split_at = remaining.len() - tail_len;
                chunks.push(remaining.split_off(split_at));
            }
            while !remaining.is_empty() {
                let split_at = remaining.len().saturating_sub(chunk_size);
                chunks.push(remaining.split_off(split_at));
            }
            chunks.reverse();

            // Every chunk must be individually sorted
            for (ci, chunk) in chunks.iter().enumerate() {
                for i in 1..chunk.len() {
                    assert!(
                        chunk[i - 1].0 <= chunk[i].0,
                        "n={n} threads={threads}: chunk {ci} not sorted at index {i} \
                         ({} > {})",
                        chunk[i - 1].0,
                        chunk[i].0,
                    );
                }
            }

            // Total record count must match
            let total: usize = chunks.iter().map(Vec::len).sum();
            assert_eq!(total, n, "n={n} threads={threads}: total mismatch");
        }
    }

    /// End-to-end sort across two temp directories. Forces multiple spill
    /// chunks and verifies the output is still correct. The round-robin
    /// distribution itself is covered by `TmpDirAllocator`'s unit tests; this
    /// test proves the plumbing from `temp_dirs(...)` through to the sort
    /// fns works and that multi-dir mode produces byte-identical output to
    /// single-dir mode.
    #[test]
    fn test_sort_with_two_temp_dirs_matches_single_dir() {
        use fgumi_sam::SamBuilder;

        let num_pairs = 200;
        let mut builder = SamBuilder::new();
        for i in 0..num_pairs {
            let _ = builder
                .add_pair()
                .name(&format!("read{i:05}"))
                .start1(i * 200 + 1)
                .start2(i * 200 + 101)
                .build();
        }

        let workdir = tempfile::tempdir().expect("workdir");
        let input = workdir.path().join("input.bam");
        let output_multi = workdir.path().join("output_multi.bam");
        let output_single = workdir.path().join("output_single.bam");
        builder.write_bam(&input).expect("write bam");

        let tmp_a = tempfile::tempdir().expect("tmp a");
        let tmp_b = tempfile::tempdir().expect("tmp b");

        // 8 KiB memory limit forces several spills across the dir rotation.
        let stats_multi = RawExternalSorter::new(SortOrder::Coordinate)
            .memory_limit(8 * 1024)
            .threads(1)
            .spill_codec(crate::codec::SpillCodec::Bgzf) // level 0 = stored mode (zstd has no level 0)
            .temp_compression(0)
            .output_compression(0)
            .temp_dirs(vec![tmp_a.path().to_path_buf(), tmp_b.path().to_path_buf()])
            .sort(&input, &output_multi)
            .expect("multi-dir sort should succeed");

        assert!(stats_multi.chunks_written >= 2, "expected multiple spill chunks");

        RawExternalSorter::new(SortOrder::Coordinate)
            .memory_limit(8 * 1024)
            .threads(1)
            .spill_codec(crate::codec::SpillCodec::Bgzf)
            .temp_compression(0)
            .output_compression(0)
            .sort(&input, &output_single)
            .expect("single-dir sort should succeed");

        let names_multi = collect_read_names(&output_multi);
        let names_single = collect_read_names(&output_single);
        assert_eq!(names_multi.len(), num_pairs * 2, "record count mismatch");
        assert_eq!(
            names_multi, names_single,
            "multi-dir and single-dir sort produced different record orders"
        );
    }

    // ========================================================================
    // read_exact_or_eof tests
    // ========================================================================

    /// A `Read` impl that returns at most one byte per `read()` call. This is
    /// legal for `Read` (which may return any `0 < n <= buf.len()`) and is what
    /// motivates using `read_exact_or_eof` over a single `read()` call for
    /// fixed-size magic-byte detection.
    struct OneBytePerRead<'a> {
        bytes: &'a [u8],
        pos: usize,
    }
    impl std::io::Read for OneBytePerRead<'_> {
        fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
            if buf.is_empty() || self.pos >= self.bytes.len() {
                return Ok(0);
            }
            buf[0] = self.bytes[self.pos];
            self.pos += 1;
            Ok(1)
        }
    }

    #[test]
    fn test_read_exact_or_eof_clean_eof() {
        // Empty stream -> Ok(false), buffer untouched.
        let mut src: &[u8] = &[];
        let mut buf = [0u8; 4];
        let filled = read_exact_or_eof(&mut src, &mut buf).expect("clean eof");
        assert!(!filled, "clean EOF should return Ok(false)");
        assert_eq!(buf, [0u8; 4]);
    }

    #[test]
    fn test_read_exact_or_eof_full_read() {
        // Source has exactly buf.len() bytes -> Ok(true), buffer filled.
        let mut src: &[u8] = b"ZSP1";
        let mut buf = [0u8; 4];
        let filled = read_exact_or_eof(&mut src, &mut buf).expect("full read");
        assert!(filled, "full read should return Ok(true)");
        assert_eq!(&buf, b"ZSP1");
    }

    #[test]
    fn test_read_exact_or_eof_short_reads_fill_buffer() {
        // Regression: a `Read` impl that returns 1 byte per call must still
        // fill the 4-byte buffer. A naive single `read()` would only get 1
        // byte and `SpillCodec::from_magic` would return None, silently
        // routing a real ZSP1 spill through the uncompressed fallback path.
        let mut src = OneBytePerRead { bytes: b"ZSP1", pos: 0 };
        let mut buf = [0u8; 4];
        let filled = read_exact_or_eof(&mut src, &mut buf).expect("short-read source");
        assert!(filled, "partial reads should still fill the buffer");
        assert_eq!(&buf, b"ZSP1");
        // And `from_magic` correctly identifies the codec from the filled
        // buffer — proving the original bug is fixed end-to-end.
        assert_eq!(
            crate::codec::SpillCodec::from_magic(&buf),
            Some(crate::codec::SpillCodec::Zstd)
        );
    }

    #[test]
    fn test_read_exact_or_eof_truncated_is_error() {
        // Source has fewer bytes than requested but is non-empty:
        // must error with UnexpectedEof (not silently report short fill).
        let mut src: &[u8] = b"ZSP"; // only 3 of 4 bytes
        let mut buf = [0u8; 4];
        let err = read_exact_or_eof(&mut src, &mut buf).expect_err("truncated read should error");
        assert_eq!(err.kind(), std::io::ErrorKind::UnexpectedEof);
    }

    // ========================================================================
    // verify_dropped_lanes / dropped_lane_error tests
    // ========================================================================

    fn full_key(cb: u64, tertiary: u64) -> TemplateKey {
        TemplateKey { primary: 0, secondary: 0, cb_hash: cb, tertiary, name_hash_upper: 0 }
    }

    const LITE: TemplateKeyVariant = TemplateKeyVariant { cb: false, tertiary: false };

    #[test]
    fn verify_detects_cb_appearing() {
        let v = verify_dropped_lanes(&full_key(0, 0), &full_key(123, 0), LITE);
        assert_eq!(v, Some(DroppedLaneViolation::Cb));
        assert_eq!(v.unwrap().key_types_token(), "cb");
    }

    #[test]
    fn verify_detects_mi_change() {
        let first = full_key(0, 0);
        let cur = full_key(0, 0b10); // mi_value bit set, library still 0
        assert_eq!(verify_dropped_lanes(&first, &cur, LITE), Some(DroppedLaneViolation::Mi));
    }

    #[test]
    fn verify_detects_library_change() {
        let first = full_key(0, 0);
        let cur = full_key(0, 1u64 << 48); // library ordinal 1
        assert_eq!(verify_dropped_lanes(&first, &cur, LITE), Some(DroppedLaneViolation::Library));
    }

    #[test]
    fn verify_passes_when_dropped_lanes_constant() {
        let variant = TemplateKeyVariant { cb: true, tertiary: false };
        assert_eq!(verify_dropped_lanes(&full_key(1, 5), &full_key(2, 5), variant), None);
    }

    #[test]
    fn verify_error_message_names_field_and_token() {
        let e = dropped_lane_error("read42", DroppedLaneViolation::Library);
        let msg = e.to_string();
        assert!(
            msg.contains("read42")
                && msg.contains("library")
                && msg.contains("--key-types library"),
            "{msg}"
        );
    }

    // ========================================================================
    // select_template_variant tests
    // ========================================================================

    fn key(cb: u64, tertiary: u64) -> TemplateKey {
        TemplateKey { primary: 1, secondary: 2, cb_hash: cb, tertiary, name_hash_upper: 3 }
    }

    #[test]
    fn auto_lite_when_no_cb_no_tertiary() {
        let v = select_template_variant(Some(&key(0, 0)), KeyTypesSpec::Auto, false);
        assert_eq!(v, TemplateKeyVariant { cb: false, tertiary: false });
        assert_eq!(v.lanes(), 3);
    }

    #[test]
    fn auto_cb_only_when_first_record_has_cb() {
        let v = select_template_variant(Some(&key(42, 0)), KeyTypesSpec::Auto, false);
        assert_eq!(v, TemplateKeyVariant { cb: true, tertiary: false });
        assert_eq!(v.lanes(), 4);
    }

    #[test]
    fn auto_tertiary_only_when_first_record_has_mi_or_library() {
        // 0xABCD is in the low-48 bits => MI present => keep tertiary, even though
        // the header library does not vary.
        let v = select_template_variant(Some(&key(0, 0xABCD)), KeyTypesSpec::Auto, false);
        assert_eq!(v, TemplateKeyVariant { cb: false, tertiary: true });
        assert_eq!(v.lanes(), 4);
    }

    #[test]
    fn auto_full_when_both_present() {
        let v = select_template_variant(Some(&key(7, 9)), KeyTypesSpec::Auto, false);
        assert_eq!(v, TemplateKeyVariant { cb: true, tertiary: true });
        assert_eq!(v.lanes(), 5);
    }

    // tertiary high bits only = a single library ordinal, no MI. With a non-varying
    // header, the constant library lane is dropped -> lite.
    #[test]
    fn auto_drops_constant_single_library_lane() {
        let k = key(0, 1u64 << 48); // library ordinal 1, no MI
        let v = select_template_variant(Some(&k), KeyTypesSpec::Auto, false);
        assert_eq!(v, TemplateKeyVariant { cb: false, tertiary: false });
        assert_eq!(v.lanes(), 3);
    }

    #[test]
    fn auto_keeps_tertiary_when_header_library_varies() {
        let k = key(0, 1u64 << 48); // library ordinal 1, no MI
        let v = select_template_variant(Some(&k), KeyTypesSpec::Auto, true); // >1 library
        assert_eq!(v, TemplateKeyVariant { cb: false, tertiary: true });
        assert_eq!(v.lanes(), 4);
    }

    #[test]
    fn auto_keeps_tertiary_for_mi_even_with_constant_library() {
        // MI present (low bits) -> keep tertiary even though header library is constant.
        let k = key(0, (1u64 << 48) | 0b10); // library ordinal 1 + MI value bit
        let v = select_template_variant(Some(&k), KeyTypesSpec::Auto, false);
        assert_eq!(v, TemplateKeyVariant { cb: false, tertiary: true });
    }

    #[test]
    fn none_forces_lite() {
        assert_eq!(
            select_template_variant(Some(&key(99, 99)), KeyTypesSpec::None, true).lanes(),
            3
        );
    }

    #[test]
    fn full_forces_all_lanes() {
        assert_eq!(select_template_variant(Some(&key(0, 0)), KeyTypesSpec::Full, false).lanes(), 5);
    }

    #[test]
    fn explicit_spec_passthrough() {
        assert_eq!(
            select_template_variant(
                None,
                KeyTypesSpec::Explicit { cb: true, tertiary: false },
                false
            )
            .lanes(),
            4
        );
    }

    #[test]
    fn empty_input_defaults_to_lite() {
        assert_eq!(select_template_variant(None, KeyTypesSpec::Auto, false).lanes(), 3);
    }

    // ========================================================================
    // Capacity estimator consistency tests (T8)
    // ========================================================================

    #[test]
    fn estimator_bytes_per_record_matches_ref_size() {
        use crate::inline::TemplateKey32;
        for (actual, expected) in [
            (std::mem::size_of::<TemplateRecordRef<TemplateKey24>>(), 40usize),
            (std::mem::size_of::<TemplateRecordRef<TemplateKey32>>(), 48),
            (std::mem::size_of::<TemplateRecordRef<TemplateKey40>>(), 56),
        ] {
            assert_eq!(actual, expected, "ref size must equal key + 16 B offset/len/pad");
        }
        let bpr_lite = (TEMPLATE_HEADER_SIZE + EST_BAM_BYTES_PER_TEMPLATE_RECORD)
            + std::mem::size_of::<TemplateRecordRef<TemplateKey24>>();
        assert_eq!(bpr_lite, 8 + 250 + 40);
    }

    // ========================================================================
    // T10: narrow-key byte-identity correctness proof
    // ========================================================================
    //
    // The whole point of `--key-types`: a narrow-key template-coordinate sort
    // must produce a record stream byte-identical to the full-key (`--key-types
    // full`) baseline on the same input. These tests build one unsorted BAM per
    // mode (with genuinely varying tags so the narrowed lane actually matters),
    // sort it under both the narrow spec and the full spec, and compare the
    // per-record BAM bodies in file order. They are non-vacuous: a sanity helper
    // asserts the constructed input produces the nonzero cb_hash / tertiary that
    // the mode claims, and the bonus `*_lite_sort_hard_errors` assertions prove a
    // mis-narrowed (lite) sort of the same input would refuse to run rather than
    // silently mis-sort.

    /// Number of distinct-position pairs per mode input — large enough that the
    /// spill-forced variant crosses a tiny memory limit and exercises the Phase-2
    /// k-way merge. Tied clusters (see `CLUSTER_SIZE` / `CLUSTER_COUNT`) are added
    /// on top of these.
    const T10_DISTINCT_PAIRS: usize = 300;

    /// Number of records in each coordinate-tied cluster. Within a cluster every
    /// record shares the same primary (tid, start1) AND secondary (start2, strands)
    /// and differs ONLY in the kept optional lane (`cb_hash` or tertiary), so that
    /// lane is the load-bearing tiebreaker above `name_hash_upper`. Made large so
    /// the chance that `name_hash` order coincidentally equals kept-lane order for
    /// every cluster (which would re-vacuate the test) is negligible.
    const CLUSTER_SIZE: usize = 8;

    /// Number of coordinate-tied clusters, each at a distinct shared coordinate.
    /// Several clusters further drive the coincidence probability to zero.
    const CLUSTER_COUNT: usize = 4;

    /// Total pairs (distinct + clustered) each mode input contains.
    const T10_PAIRS: usize = T10_DISTINCT_PAIRS + CLUSTER_SIZE * CLUSTER_COUNT;

    /// Collect each record's raw BAM bytes (excluding the header) in file order.
    ///
    /// Compares at the record-stream level rather than the file level: BGZF block
    /// batching is nondeterministic, so two semantically identical BAMs can differ
    /// byte-for-byte on disk while their decoded record bodies are identical.
    fn collect_record_bytes(path: &Path) -> Vec<Vec<u8>> {
        use crate::read_ahead::RawReadAheadReader;
        let (reader, _) = create_raw_bam_reader(path, 1).expect("reader");
        RawReadAheadReader::new(reader).map(|rec| rec.as_ref().to_vec()).collect()
    }

    /// Sort `input` with the given key-types spec and return the record-byte stream.
    ///
    /// `cell_tag(SamTag::CB)` is REQUIRED: `RawExternalSorter::new` defaults the
    /// cell tag to `None`, so without this the CB lane is never populated and any
    /// single-cell assertion would pass vacuously (`cb_hash` stays 0). This mirrors
    /// the CLI's `parse_cell_tag`.
    fn sort_and_collect(input: &Path, dir: &Path, tag: &str, spec: KeyTypesSpec) -> Vec<Vec<u8>> {
        let out = dir.join(format!("{tag}.bam"));
        RawExternalSorter::new(SortOrder::TemplateCoordinate)
            .output_compression(0)
            .cell_tag(SamTag::CB)
            .key_types(spec)
            .sort(input, &out)
            .expect("sort");
        collect_record_bytes(&out)
    }

    /// Insert a read group carrying an `LB:` (library) field into a builder header.
    ///
    /// Used by the multi-library inputs to realize two distinct library ordinals
    /// (sorted alphabetically, then 1-based) so the tertiary lane genuinely varies.
    fn add_rg_with_library(header: &mut Header, rg_id: &str, library: &str) {
        let rg = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from(library))
            .build()
            .expect("valid read group");
        header.read_groups_mut().insert(BString::from(rg_id), rg);
    }

    /// Extract the full-width template key for the first record of a BAM, using a
    /// fresh `LibraryLookup` from its header. Lets a builder sanity-check that the
    /// tags it wrote actually produce the nonzero `cb_hash` / tertiary it claims.
    fn first_record_full_key(path: &Path) -> TemplateKey40 {
        use crate::read_ahead::RawReadAheadReader;
        let (reader, header) = create_raw_bam_reader(path, 1).expect("reader");
        let lib = LibraryLookup::from_header(&header);
        let first = RawReadAheadReader::new(reader).next().expect("at least one record");
        extract_template_key_inline(first.as_ref(), &lib, Some(SamTag::CB), &cb_hasher())
    }

    /// Scramble a cluster member's slot index into a kept-lane value index so the
    /// kept-lane order is a NON-identity permutation of the name order within the
    /// cluster. `7` is coprime to `CLUSTER_SIZE` (8), so this is a bijection that
    /// fixes no element to its name position (the `+3` shift removes the only fixed
    /// point). This anti-correlates the kept lane with `name_hash` so a wrong
    /// narrowing that falls back to `name_hash` cannot coincidentally reproduce the
    /// correct (kept-lane) order.
    fn scrambled_lane_index(member: usize) -> usize {
        (member * 7 + 3) % CLUSTER_SIZE
    }

    /// Shared start1 for cluster `c`. Far above any distinct-position record
    /// (`i * 200` for `i < T10_DISTINCT_PAIRS`, i.e. below ~`60_000`) so clusters
    /// never collide with the general path and each cluster sits at its own
    /// coordinate. All members of a cluster share this exact coordinate.
    fn cluster_start1(c: usize) -> usize {
        10_000_000 + c * 10_000
    }

    /// Build an unsorted BAM with `T10_DISTINCT_PAIRS` uniquely-positioned pairs
    /// PLUS `CLUSTER_COUNT` coordinate-tied clusters of `CLUSTER_SIZE` pairs each.
    ///
    /// `tagger` decorates a uniquely-positioned pair (CB/MI/RG). `cluster_tagger`
    /// decorates a cluster member, receiving the scrambled lane index so the caller
    /// can set the kept optional lane (CB / MI / RG) to a value whose order is
    /// anti-correlated with the member's name. The uniquely-positioned records
    /// exercise the general sort/spill path; the clusters make the kept lane the
    /// load-bearing tiebreaker (shared primary+secondary, differ only in kept lane).
    fn build_mode_bam(
        dir: &Path,
        name: &str,
        mut header: Header,
        tagger: impl Fn(PairBuilder<'_>, usize) -> PairBuilder<'_>,
        cluster_tagger: impl Fn(PairBuilder<'_>, usize) -> PairBuilder<'_>,
    ) -> PathBuf {
        let mut builder = SamBuilder::new();
        builder.header = std::mem::take(&mut header);

        // Coordinate-tied clusters first: each member shares start1/start2/strands
        // with its cluster siblings and differs only in the kept optional lane.
        for c in 0..CLUSTER_COUNT {
            let start1 = cluster_start1(c);
            let start2 = start1 + 100;
            for member in 0..CLUSTER_SIZE {
                let pair = builder
                    .add_pair()
                    .name(&format!("read_c{c}_{member:02}"))
                    .start1(start1)
                    .start2(start2);
                let _ = cluster_tagger(pair, scrambled_lane_index(member)).build();
            }
        }

        // Uniquely-positioned pairs: interleave positions (even i ascend, odd i
        // descend) so file order is far from template-coordinate order.
        for i in 0..T10_DISTINCT_PAIRS {
            let pos = if i % 2 == 0 { i * 200 } else { (T10_DISTINCT_PAIRS - i) * 200 };
            let pair =
                builder.add_pair().name(&format!("read{i:05}")).start1(pos + 1).start2(pos + 101);
            let _ = tagger(pair, i).build();
        }

        let path = dir.join(format!("{name}.bam"));
        builder.write_bam(&path).expect("write bam");
        path
    }

    /// Distinct CB barcodes, one per cluster slot, so a tied cluster's members get
    /// distinct `cb_hash` values (the load-bearing single-cell tiebreaker).
    const CLUSTER_BARCODES: [&str; CLUSTER_SIZE] = [
        "AAAACCCC", "GGGGTTTT", "ACGTACGT", "TTTTGGGG", "CCCCAAAA", "TTTTACGT", "GACTGACT",
        "TGCATGCA",
    ];

    /// bulk (lite): no CB, no MI, single default read group. No load-bearing opt
    /// lane, so the clusters carry no distinguishing tag — bulk byte-identity is
    /// inherently about the general path, and the cluster members tie down to
    /// `name_hash` in both baseline and narrow (lite drops both lanes anyway).
    fn build_bulk_bam(dir: &Path) -> PathBuf {
        build_mode_bam(
            dir,
            "bulk",
            SamBuilder::new().header.clone(),
            |pair, _| pair,
            |pair, _| pair,
        )
    }

    /// single-cell: every pair carries `CB:Z:<barcode>`; no MI. Cluster members get
    /// DISTINCT barcodes (scrambled vs name) so `cb_hash` is the load-bearing lane.
    fn build_single_cell_bam(dir: &Path) -> PathBuf {
        const BARCODES: [&str; 4] = ["AAAACCCC", "GGGGTTTT", "ACGTACGT", "TTTTGGGG"];
        build_mode_bam(
            dir,
            "single_cell",
            SamBuilder::new().header.clone(),
            |pair, i| pair.attr("CB", BARCODES[i % BARCODES.len()]),
            |pair, lane| pair.attr("CB", CLUSTER_BARCODES[lane]),
        )
    }

    /// post-group: every pair carries `MI:Z:<id>`; no CB; single library. Cluster
    /// members get DISTINCT MI ids (scrambled vs name), so the tertiary lane
    /// (MI bits) is the load-bearing tiebreaker. Ids are `1..=CLUSTER_SIZE` so
    /// tertiary != 0 for every member.
    fn build_post_group_bam(dir: &Path) -> PathBuf {
        build_mode_bam(
            dir,
            "post_group",
            SamBuilder::new().header.clone(),
            |pair, i| pair.attr("MI", format!("{}", (i % 5) + 1)),
            |pair, lane| pair.attr("MI", format!("{}", lane + 1)),
        )
    }

    /// multi-lib: header with TWO `@RG` lines whose `LB:` differ; pairs split across
    /// the two RGs; no CB, no MI. The two libraries realize ordinals 1 and 2. Within
    /// a cluster, members alternate read groups by the SCRAMBLED lane index so the
    /// tertiary library bits (high-16) vary load-bearingly across the cluster and
    /// the library order is anti-correlated with the name order.
    fn build_multi_lib_bam(dir: &Path) -> PathBuf {
        let mut header = SamBuilder::new().header.clone();
        add_rg_with_library(&mut header, "rgA", "LibAlpha");
        add_rg_with_library(&mut header, "rgB", "LibBeta");
        build_mode_bam(
            dir,
            "multi_lib",
            header,
            |pair, i| pair.attr("RG", if i % 2 == 0 { "rgA" } else { "rgB" }),
            |pair, lane| pair.attr("RG", if lane % 2 == 0 { "rgA" } else { "rgB" }),
        )
    }

    /// full: CB + MI + multi-library together — all optional lanes vary. Cluster
    /// members get distinct CB and MI (scrambled vs name) plus alternating RG.
    /// Baseline (full) == narrow (full) so this mode is inherently trivial, but the
    /// clusters keep it meaningful (the kept lanes still vary load-bearingly).
    fn build_full_bam(dir: &Path) -> PathBuf {
        const BARCODES: [&str; 4] = ["AAAACCCC", "GGGGTTTT", "ACGTACGT", "TTTTGGGG"];
        let mut header = SamBuilder::new().header.clone();
        add_rg_with_library(&mut header, "rgA", "LibAlpha");
        add_rg_with_library(&mut header, "rgB", "LibBeta");
        build_mode_bam(
            dir,
            "full",
            header,
            |pair, i| {
                pair.attr("CB", BARCODES[i % BARCODES.len()])
                    .attr("MI", format!("{}", i % 5))
                    .attr("RG", if i % 2 == 0 { "rgA" } else { "rgB" })
            },
            |pair, lane| {
                pair.attr("CB", CLUSTER_BARCODES[lane])
                    .attr("MI", format!("{}", lane + 1))
                    .attr("RG", if lane % 2 == 0 { "rgA" } else { "rgB" })
            },
        )
    }

    /// For each mode, the narrow-key sort must produce a record stream IDENTICAL to
    /// the `--key-types full` baseline, and the `Auto`-detected variant must match
    /// it too. The narrow spec is the minimal one for the mode; an over-narrow
    /// (lite) sort of single-cell / post-group / multi-lib inputs would instead
    /// hard-error (see `*_lite_sort_hard_errors`), proving these are not vacuous.
    #[rstest::rstest]
    #[case::bulk(build_bulk_bam as fn(&Path) -> PathBuf, KeyTypesSpec::None)]
    #[case::single_cell(
        build_single_cell_bam as fn(&Path) -> PathBuf,
        KeyTypesSpec::Explicit { cb: true, tertiary: false }
    )]
    #[case::post_group(
        build_post_group_bam as fn(&Path) -> PathBuf,
        KeyTypesSpec::Explicit { cb: false, tertiary: true }
    )]
    #[case::multi_lib(
        build_multi_lib_bam as fn(&Path) -> PathBuf,
        KeyTypesSpec::Explicit { cb: false, tertiary: true }
    )]
    #[case::full(build_full_bam as fn(&Path) -> PathBuf, KeyTypesSpec::Full)]
    fn narrow_key_sort_byte_identical_to_full(
        #[case] build: fn(&Path) -> PathBuf,
        #[case] narrow: KeyTypesSpec,
    ) {
        let dir = tempfile::tempdir().unwrap();
        let input = build(dir.path());

        let baseline = sort_and_collect(&input, dir.path(), "full_baseline", KeyTypesSpec::Full);
        assert_eq!(baseline.len(), T10_PAIRS * 2, "baseline must retain every record");

        let narrowed = sort_and_collect(&input, dir.path(), "narrow", narrow);
        assert_eq!(narrowed, baseline, "narrow-key record stream must equal full-key stream");

        let auto = sort_and_collect(&input, dir.path(), "auto", KeyTypesSpec::Auto);
        assert_eq!(auto, baseline, "auto-detected variant must equal full-key stream");
    }

    /// single-library: one `@RG` carrying an `LB:`, every pair assigned to it; no CB,
    /// no MI. The lone library realizes ordinal 1, so the first record's tertiary is
    /// nonzero (high bits) yet the library lane is CONSTANT across every record — the
    /// WES/WGS scenario. `Auto` must therefore drop the tertiary lane to the 24-byte
    /// lite key while staying byte-identical to the full-key baseline.
    fn build_single_library_bam(dir: &Path) -> PathBuf {
        // Attach a library to the builder's single default read group ("A") so the
        // header realizes exactly ONE distinct library ordinal (1). Reads default to
        // RG "A", so every record maps to that lone library — the tertiary high bits
        // are a constant nonzero ordinal, with no MI and no CB.
        let mut header = SamBuilder::new().header.clone();
        let rg = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("LibSolo"))
            .build()
            .expect("valid read group");
        // The builder defaults to a single read group "A"; replace it in place with a
        // library-bearing one so there is exactly one read group (one ordinal).
        header.read_groups_mut().clear();
        header.read_groups_mut().insert(BString::from("A"), rg);
        build_mode_bam(dir, "single_library", header, |pair, _| pair, |pair, _| pair)
    }

    /// WES/WGS regression: a single-library input (one constant library ordinal, no
    /// CB, no MI) must auto-narrow to the 3-lane lite key AND sort byte-identically to
    /// the full-key baseline. Proves the constant library lane is provisionally
    /// droppable under `Auto` and that dropping it does not change the record stream.
    #[test]
    fn auto_drops_constant_library_lane_byte_identical() {
        let dir = tempfile::tempdir().unwrap();
        let input = build_single_library_bam(dir.path());

        // Production provisioning logic, reproduced: the first record carries a
        // nonzero (library) tertiary, but the header realizes a single ordinal, so
        // Auto must pick the 3-lane lite variant.
        let first = first_record_full_key(&input);
        let (_, header) = create_raw_bam_reader(&input, 1).expect("reader");
        let lib_lookup = LibraryLookup::from_header(&header);
        let header_library_varies = lib_lookup.distinct_header_ordinals() > 1;
        assert!(!header_library_varies, "single-library header must not vary");
        assert_ne!(first.tertiary, 0, "first record carries a (constant) library ordinal");

        let variant =
            select_template_variant(Some(&first), KeyTypesSpec::Auto, header_library_varies);
        assert_eq!(variant.lanes(), 3, "auto must narrow a single-library input to the lite key");

        // Byte-identity: auto (lite) == full-key baseline.
        let baseline = sort_and_collect(&input, dir.path(), "full_baseline", KeyTypesSpec::Full);
        assert_eq!(baseline.len(), T10_PAIRS * 2, "baseline must retain every record");
        let auto = sort_and_collect(&input, dir.path(), "auto", KeyTypesSpec::Auto);
        assert_eq!(auto, baseline, "auto (lite) record stream must equal full-key stream");

        // Cross-check: an explicit lite (`--key-types none`) sort is also identical,
        // proving the lite key is genuinely valid for this input.
        let lite = sort_and_collect(&input, dir.path(), "lite", KeyTypesSpec::None);
        assert_eq!(lite, baseline, "explicit lite record stream must equal full-key stream");
    }

    /// Spill-forced byte-identity: a tiny memory limit forces Phase 1 to spill and
    /// Phase 2 to k-way merge `PooledChunkWriter::<K>` chunks. This is the path
    /// where the serialized narrow-key width must match the Phase-1 width; a
    /// mismatch would corrupt the merge. Asserted on the bulk (lite) input.
    #[test]
    fn narrow_key_sort_byte_identical_when_spilling() {
        let dir = tempfile::tempdir().unwrap();
        let input = build_bulk_bam(dir.path());

        let sort_spilling = |tag: &str, spec: KeyTypesSpec| -> (Vec<Vec<u8>>, usize) {
            let out = dir.path().join(format!("{tag}.bam"));
            let stats = RawExternalSorter::new(SortOrder::TemplateCoordinate)
                .memory_limit(32 * 1024)
                .max_temp_files(0)
                .spill_codec(crate::codec::SpillCodec::Bgzf)
                .temp_compression(0)
                .output_compression(0)
                .cell_tag(SamTag::CB)
                .key_types(spec)
                .sort(&input, &out)
                .expect("sort");
            (collect_record_bytes(&out), stats.chunks_written)
        };

        let (baseline, base_chunks) = sort_spilling("spill_full", KeyTypesSpec::Full);
        let (narrowed, narrow_chunks) = sort_spilling("spill_narrow", KeyTypesSpec::None);

        assert!(base_chunks >= 2, "baseline must spill multiple chunks, got {base_chunks}");
        assert!(narrow_chunks >= 2, "narrow must spill multiple chunks, got {narrow_chunks}");
        assert_eq!(narrowed, baseline, "spilled narrow stream must equal spilled full stream");
    }

    /// Spill-forced byte-identity for ALL key widths: a tiny memory limit forces Phase 1
    /// to spill and Phase 2 to k-way merge `PooledChunkWriter::<K>` chunks across all four
    /// key widths (24-byte lite, 32-byte `CbKey32`, 32-byte `TertKey32`, 40-byte `TemplateKey40`).
    ///
    /// The existing `narrow_key_sort_byte_identical_when_spilling` test only covers the
    /// 24-byte (bulk/lite) path. A width/serialization regression in the 32- or 40-byte
    /// Phase-2 merge path would not be caught by that test alone. This test closes the gap
    /// by asserting both that each width's narrow-key spilled stream is byte-identical to
    /// the Full baseline spilled stream, AND that both runs actually spilled (>=2 chunks),
    /// so the Phase-2 merge path at each key width is genuinely exercised.
    #[rstest::rstest]
    #[case::bulk(build_bulk_bam as fn(&Path) -> PathBuf, KeyTypesSpec::None)]
    #[case::single_cell(
        build_single_cell_bam as fn(&Path) -> PathBuf,
        KeyTypesSpec::Explicit { cb: true, tertiary: false }
    )]
    #[case::post_group(
        build_post_group_bam as fn(&Path) -> PathBuf,
        KeyTypesSpec::Explicit { cb: false, tertiary: true }
    )]
    #[case::multi_lib(
        build_multi_lib_bam as fn(&Path) -> PathBuf,
        KeyTypesSpec::Explicit { cb: false, tertiary: true }
    )]
    #[case::full(build_full_bam as fn(&Path) -> PathBuf, KeyTypesSpec::Full)]
    fn narrow_key_spill_merge_byte_identical_all_widths(
        #[case] build: fn(&Path) -> PathBuf,
        #[case] narrow: KeyTypesSpec,
    ) {
        let dir = tempfile::tempdir().unwrap();
        let input = build(dir.path());

        let sort_spilling = |tag: &str, spec: KeyTypesSpec| -> (Vec<Vec<u8>>, usize) {
            let out = dir.path().join(format!("{tag}.bam"));
            let stats = RawExternalSorter::new(SortOrder::TemplateCoordinate)
                .memory_limit(32 * 1024)
                .max_temp_files(0)
                .spill_codec(crate::codec::SpillCodec::Bgzf)
                .temp_compression(0)
                .output_compression(0)
                .cell_tag(SamTag::CB)
                .key_types(spec)
                .sort(&input, &out)
                .expect("sort");
            (collect_record_bytes(&out), stats.chunks_written)
        };

        let (baseline, base_chunks) = sort_spilling("spill_full_baseline", KeyTypesSpec::Full);
        let (narrowed, narrow_chunks) = sort_spilling("spill_narrow", narrow);

        assert!(
            base_chunks >= 2,
            "Full-key baseline must spill multiple chunks (Phase-2 merge must run), \
             got {base_chunks}"
        );
        assert!(
            narrow_chunks >= 2,
            "Narrow-key sort must spill multiple chunks (Phase-2 merge must run at this \
             key width), got {narrow_chunks}"
        );
        assert_eq!(
            narrowed, baseline,
            "Narrow-key spilled record stream must be byte-identical to Full-key spilled \
             record stream — a mismatch indicates a width/serialization regression in the \
             Phase-2 merge path"
        );
    }

    /// Sanity (anti-vacuity): each mode's input must produce the nonzero `cb_hash` /
    /// tertiary it claims on its first record, otherwise the byte-identity tests
    /// would compare two identically-narrowed streams and prove nothing.
    #[test]
    fn mode_inputs_realize_claimed_lanes() {
        let dir = tempfile::tempdir().unwrap();

        let bulk = first_record_full_key(&build_bulk_bam(dir.path()));
        assert_eq!(bulk.cb_hash, 0, "bulk must have no CB");
        assert_eq!(bulk.tertiary, 0, "bulk must have no MI/library");

        let cell = first_record_full_key(&build_single_cell_bam(dir.path()));
        assert_ne!(cell.cb_hash, 0, "single-cell first record must have nonzero cb_hash");
        assert_eq!(cell.tertiary, 0, "single-cell must have no MI/library");

        let pg = first_record_full_key(&build_post_group_bam(dir.path()));
        assert_eq!(pg.cb_hash, 0, "post-group must have no CB");
        assert_ne!(pg.tertiary, 0, "post-group first record must have nonzero MI tertiary");
        assert_eq!(pg.tertiary >> 48, 0, "post-group must have no library (ordinal 0)");

        let multi = first_record_full_key(&build_multi_lib_bam(dir.path()));
        assert_eq!(multi.cb_hash, 0, "multi-lib must have no CB");
        assert_ne!(multi.tertiary >> 48, 0, "multi-lib first record must have nonzero library");

        let full = first_record_full_key(&build_full_bam(dir.path()));
        assert_ne!(full.cb_hash, 0, "full first record must have nonzero cb_hash");
        assert_ne!(full.tertiary, 0, "full first record must have nonzero tertiary");
        assert_ne!(full.tertiary >> 48, 0, "full first record must have nonzero library ordinal");
    }

    /// Guard the single-cell cluster against a silent `cb_hash` collision: every
    /// barcode in `CLUSTER_BARCODES` must hash to a DISTINCT `cb_hash` value under
    /// `cb_hasher()`, otherwise two cluster members would share the same sort key
    /// and the tiebreaker designed to make the narrow-key test non-vacuous would
    /// silently collapse.
    ///
    /// Uses `cb_hasher().hash_one(bc.as_bytes())` — exactly the expression in
    /// `extract_template_key_inline` (`cb_hasher.hash_one(cb_bytes)` where
    /// `cb_bytes` is the raw string value of the CB aux tag, i.e. ASCII barcode
    /// bytes without NUL terminator).
    #[test]
    fn cluster_barcodes_have_distinct_cb_hashes() {
        let hasher = cb_hasher();
        let mut hash_vals: Vec<u64> =
            CLUSTER_BARCODES.iter().map(|bc| hasher.hash_one(bc.as_bytes())).collect();
        let n = hash_vals.len();
        hash_vals.sort_unstable();
        hash_vals.dedup();
        assert_eq!(
            hash_vals.len(),
            n,
            "cluster barcodes must hash to distinct cb_hash values, else the single-cell \
             cluster is vacuous"
        );
    }

    /// A narrow (tertiary-only) sort output must re-verify correct at FULL width via
    /// `core_cmp`. This locks the consistency between the narrowed key used for
    /// sorting and the full key used by `fgumi sort --verify`: if narrowing dropped
    /// ordering information the merge depended on, the full-width re-check would
    /// find a violation.
    #[test]
    fn narrow_sort_output_passes_full_width_verify() {
        use std::cmp::Ordering;

        let dir = tempfile::tempdir().unwrap();
        let input = build_post_group_bam(dir.path());
        let out = dir.path().join("narrow.bam");
        RawExternalSorter::new(SortOrder::TemplateCoordinate)
            .output_compression(0)
            .cell_tag(SamTag::CB)
            .key_types(KeyTypesSpec::Explicit { cb: false, tertiary: true })
            .sort(&input, &out)
            .expect("sort");

        // Header only (for LibraryLookup); the auto reader is dropped immediately.
        let (_, header) = create_raw_bam_reader(&out, 1).expect("header");
        let lib = LibraryLookup::from_header(&header);
        let hasher = cb_hasher();

        // verify_sort_order wants RawBamRecordReader<File> — build it like execute_verify.
        let file = std::fs::File::open(&out).expect("open");
        let mut raw_reader = crate::reader::RawBamRecordReader::new(file).expect("reader");
        raw_reader.skip_header().expect("skip header");

        let (total, violations, first) = crate::verify::verify_sort_order(
            raw_reader,
            |bam| extract_template_key_inline(bam, &lib, Some(SamTag::CB), &hasher),
            |cur: &TemplateKey40, prev: &TemplateKey40| cur.core_cmp(prev) == Ordering::Less,
        )
        .expect("verify runs");

        assert_eq!(total, (T10_PAIRS * 2) as u64, "verify must see every record");
        assert_eq!(
            violations, 0,
            "narrow output must re-verify at full width (first={first:?}, total={total})"
        );
    }

    // ------------------------------------------------------------------------
    // T6 (deferred) hard-error integration tests
    //
    // Build a BAM whose FIRST pair lacks an optional field and a LATER pair
    // carries/changes it, then sort with the over-narrow spec that drops that
    // lane. The decode-time `verify_dropped_lanes` check must hard-error with a
    // message naming the `--key-types <token>` that re-includes the lane, rather
    // than silently mis-sorting. Each test forces the relevant lane to vary
    // across records while holding the dropped lanes constant; attribution
    // (cb / mi / library) follows `verify_dropped_lanes`.
    // ------------------------------------------------------------------------

    /// Build a two-pair BAM where the first pair omits a tag the second supplies.
    fn build_first_then_appears(
        dir: &Path,
        name: &str,
        header: Header,
        first: impl Fn(PairBuilder<'_>) -> PairBuilder<'_>,
        second: impl Fn(PairBuilder<'_>) -> PairBuilder<'_>,
    ) -> PathBuf {
        let mut builder = SamBuilder::new();
        builder.header = header;
        let _ = first(builder.add_pair().name("read00000").start1(1).start2(101)).build();
        let _ = second(builder.add_pair().name("read00001").start1(201).start2(301)).build();
        let path = dir.join(format!("{name}.bam"));
        builder.write_bam(&path).expect("write bam");
        path
    }

    /// CB appears after the first record under `--key-types none` → hard error
    /// instructing `--key-types cb`.
    #[test]
    fn lite_sort_hard_errors_when_cb_appears() {
        let dir = tempfile::tempdir().unwrap();
        let input = build_first_then_appears(
            dir.path(),
            "cb_appears",
            SamBuilder::new().header.clone(),
            |pair| pair,
            |pair| pair.attr("CB", "AAAACCCC"),
        );
        let out = dir.path().join("out.bam");
        let err = RawExternalSorter::new(SortOrder::TemplateCoordinate)
            .output_compression(0)
            .cell_tag(SamTag::CB)
            .key_types(KeyTypesSpec::None)
            .sort(&input, &out)
            .expect_err("CB appearing after a CB-free first record must hard-error under lite");
        let msg = err.to_string();
        assert!(msg.contains("--key-types cb"), "expected --key-types cb guidance, got: {msg}");
    }

    /// MI appears/changes after the first record under `--key-types none` → hard
    /// error instructing `--key-types mi`.
    #[test]
    fn lite_sort_hard_errors_when_mi_appears() {
        let dir = tempfile::tempdir().unwrap();
        let input = build_first_then_appears(
            dir.path(),
            "mi_appears",
            SamBuilder::new().header.clone(),
            |pair| pair,
            |pair| pair.attr("MI", "7"),
        );
        let out = dir.path().join("out.bam");
        let err = RawExternalSorter::new(SortOrder::TemplateCoordinate)
            .output_compression(0)
            .cell_tag(SamTag::CB)
            .key_types(KeyTypesSpec::None)
            .sort(&input, &out)
            .expect_err("MI appearing after an MI-free first record must hard-error under lite");
        let msg = err.to_string();
        assert!(msg.contains("--key-types mi"), "expected --key-types mi guidance, got: {msg}");
    }

    /// Library differs after the first record under `--key-types none` → hard error
    /// instructing `--key-types library`. The first pair carries a read group whose
    /// `LB:` realizes ordinal 1; the second pair carries NO RG tag, realizing
    /// ordinal 0 — exercising the "header undercounts realized libraries" case
    /// where a later record drops to an unseen library ordinal.
    #[test]
    fn lite_sort_hard_errors_when_library_differs() {
        let dir = tempfile::tempdir().unwrap();
        // Header has exactly one read group (with LB) → ordinal 1; the default
        // SamBuilder "A" group is replaced so the only library is LibAlpha.
        let mut header = Header::builder()
            .add_reference_sequence(
                BString::from("chr1"),
                Map::<noodles::sam::header::record::value::map::ReferenceSequence>::new(
                    std::num::NonZeroUsize::new(200_000_000).expect("nonzero"),
                ),
            )
            .build();
        add_rg_with_library(&mut header, "rgA", "LibAlpha");

        // First pair: RG=rgA (ordinal 1). Second pair: RG points at an id absent
        // from the header → realized ordinal 0, differing from the first.
        let input = build_first_then_appears(
            dir.path(),
            "library_differs",
            header,
            |pair| pair.attr("RG", "rgA"),
            |pair| pair.attr("RG", "rgUnknown"),
        );
        let out = dir.path().join("out.bam");
        let err = RawExternalSorter::new(SortOrder::TemplateCoordinate)
            .output_compression(0)
            .cell_tag(SamTag::CB)
            .key_types(KeyTypesSpec::None)
            .sort(&input, &out)
            .expect_err("a differing library ordinal must hard-error under lite");
        let msg = err.to_string();
        assert!(
            msg.contains("--key-types library"),
            "expected --key-types library guidance, got: {msg}"
        );
    }

    /// Regression (#375): a single-library header with MULTIPLE read groups (all
    /// sharing one `LB:`) must sort under auto `--key-types` WITHOUT a false
    /// dropped-lane "library" violation.
    ///
    /// Every read group realizes the SAME library ordinal, so
    /// `distinct_header_ordinals() == 1`, the tertiary library lane is genuinely
    /// constant, and auto correctly drops it. Decode-time verify must therefore
    /// see a constant library ordinal across records even though they carry
    /// DIFFERENT (but same-library) RG tags. Mirrors production multi-lane WGS
    /// BAMs (e.g. 12 `@RG` lines, one `LB`) that regressed to a spurious
    /// "carries a library value absent from the input's first record" error.
    #[test]
    fn auto_single_library_multi_readgroup_sorts_without_false_violation() {
        let dir = tempfile::tempdir().unwrap();
        // Build the header from scratch (no SamBuilder default read group, which
        // would add a second library) so the ONLY libraries are the ones we add.
        let mut header = Header::builder()
            .add_reference_sequence(
                BString::from("chr1"),
                Map::<noodles::sam::header::record::value::map::ReferenceSequence>::new(
                    std::num::NonZeroUsize::new(200_000_000).expect("nonzero"),
                ),
            )
            .build();
        // Two read groups, SAME library -> exactly one distinct library ordinal.
        add_rg_with_library(&mut header, "rgA", "LibAlpha");
        add_rg_with_library(&mut header, "rgB", "LibAlpha");
        assert_eq!(
            LibraryLookup::from_header(&header).distinct_header_ordinals(),
            1,
            "two read groups sharing one library must realize a single ordinal"
        );

        // Pairs split across the two read groups; both resolve to LibAlpha.
        let input = build_mode_bam(
            dir.path(),
            "single_lib_multi_rg",
            header,
            |pair, i| pair.attr("RG", if i % 2 == 0 { "rgA" } else { "rgB" }),
            |pair, lane| pair.attr("RG", if lane % 2 == 0 { "rgA" } else { "rgB" }),
        );
        let out = dir.path().join("out.bam");
        RawExternalSorter::new(SortOrder::TemplateCoordinate)
            .output_compression(0)
            .cell_tag(SamTag::CB)
            .key_types(KeyTypesSpec::Auto)
            .sort(&input, &out)
            .expect(
                "single-library input with multiple read groups must sort without a \
                 false dropped-lane library violation",
            );
    }

    /// Regression (#375): a single-library input that contains a COMPLETELY
    /// UNMAPPED read (both mates unmapped) must sort under auto `--key-types`
    /// WITHOUT a spurious dropped-lane "library" violation.
    ///
    /// `extract_template_key_inline` packs the read's library ordinal into the
    /// full key for mapped reads (and for unmapped-with-mapped-mate reads), but
    /// `TemplateKey::unmapped` zeroed the whole `tertiary` lane — so a
    /// fully-unmapped read carrying a valid RG realized library ordinal 0 while
    /// its mapped, same-library peers realized ordinal 1. With a single library
    /// the auto path drops the (provably-constant) library lane, and decode-time
    /// verify then saw 1 (the mapped first record) vs 0 (the unmapped read) and
    /// hard-errored with "carries a library value absent from the input's first
    /// record". Reproduces the 1kg WGS regression (unmapped pairs at the tail).
    #[test]
    fn auto_single_library_unmapped_read_sorts_without_false_violation() {
        let dir = tempfile::tempdir().unwrap();
        // Single read group with a library -> exactly one library ordinal, so the
        // auto path drops the tertiary library lane and relies on decode verify.
        let mut header = Header::builder()
            .add_reference_sequence(
                BString::from("chr1"),
                Map::<noodles::sam::header::record::value::map::ReferenceSequence>::new(
                    std::num::NonZeroUsize::new(200_000_000).expect("nonzero"),
                ),
            )
            .build();
        add_rg_with_library(&mut header, "rgA", "LibAlpha");
        assert_eq!(LibraryLookup::from_header(&header).distinct_header_ordinals(), 1);

        let mut builder = SamBuilder::new();
        builder.header = header;
        // A normal mapped pair (library ordinal 1) ...
        let _ =
            builder.add_pair().name("mapped").start1(1_000).start2(1_100).attr("RG", "rgA").build();
        // ... and a COMPLETELY unmapped pair carrying the same RG. Its realized
        // library ordinal is also 1, but the old unmapped key zeroed it.
        let _ =
            builder.add_pair().name("unmapped").unmapped1().unmapped2().attr("RG", "rgA").build();
        let input = dir.path().join("in.bam");
        builder.write_bam(&input).expect("write bam");

        let out = dir.path().join("out.bam");
        RawExternalSorter::new(SortOrder::TemplateCoordinate)
            .output_compression(0)
            .cell_tag(SamTag::CB)
            .key_types(KeyTypesSpec::Auto)
            .sort(&input, &out)
            .expect(
                "a completely-unmapped read in a single-library input must not trigger a \
                 false dropped-lane library violation",
            );
    }
}

#[cfg(test)]
mod from_slots_merge_tests {
    //! Direct unit tests for the slot-backed `MergeDriver::from_slots` merge
    //! driver: cross-block record parsing (incl. records spanning many blocks),
    //! non-blocking stall/resume on empty-non-EOF slots, the embedded and
    //! non-embedded sort-key arms, and memory-chunk mixing. Ported verbatim
    //! from the issue-#330 source branch (the slimmed `from_slots` API is
    //! identical at this layer). These exercise the parser/state-machine paths
    //! that the end-to-end `three_step_chain_*` tests do not isolate.
    use super::*;
    use std::io::{BufReader, Read, Write};
    use std::sync::Arc as StdArc;

    #[derive(Default, PartialEq, Eq, PartialOrd, Ord, Clone, Copy, Debug)]
    struct TestKey(u64);

    impl RawSortKey for TestKey {
        const SERIALIZED_SIZE: Option<usize> = Some(8);
        const EMBEDDED_IN_RECORD: bool = false;
        fn extract(_bam: &[u8], _ctx: &crate::keys::SortContext) -> Self {
            unimplemented!(
                "TestKey is only used in MergeDriver::from_slots tests, not for ingestion"
            )
        }
        fn write_to<W: Write>(&self, w: &mut W) -> std::io::Result<()> {
            w.write_all(&self.0.to_le_bytes())
        }
        fn read_from<R: Read>(r: &mut R) -> std::io::Result<Self> {
            let mut buf = [0u8; 8];
            r.read_exact(&mut buf)?;
            Ok(TestKey(u64::from_le_bytes(buf)))
        }
    }

    /// Serialize a sequence of `(TestKey, record_bytes)` pairs in the spill
    /// file's wire format: `[key(8)][len(4)][record(len)]` per entry.
    fn serialize_records(records: &[(TestKey, Vec<u8>)]) -> Vec<u8> {
        let mut out = Vec::new();
        for (key, rec) in records {
            key.write_to(&mut out).unwrap();
            #[allow(clippy::cast_possible_truncation)]
            let len = rec.len() as u32;
            out.write_all(&len.to_le_bytes()).unwrap();
            out.write_all(rec).unwrap();
        }
        out
    }

    /// Build a populated `SortMergeSlot` with `file_id` whose decompressed
    /// queue contains `records` partitioned across `n_blocks` "blocks"
    /// (so the slot parser exercises cross-block reads). Sets `queue_eof`
    /// so the slot reports drained once consumer empties it.
    fn populated_slot(
        file_id: u32,
        records: &[(TestKey, Vec<u8>)],
        n_blocks: usize,
    ) -> StdArc<SortMergeSlot> {
        assert!(n_blocks >= 1, "n_blocks must be >= 1");
        let bytes = serialize_records(records);
        let block_size = bytes.len().div_ceil(n_blocks);
        let slot = StdArc::new(SortMergeSlot::new(
            file_id,
            BufReader::new(tempfile::tempfile().expect("tempfile")),
            crate::codec::SpillCodec::Bgzf,
        ));
        {
            let mut dec = slot.decompressed.lock().unwrap();
            for chunk in bytes.chunks(block_size) {
                dec.push_back(chunk.to_vec());
            }
            slot.queue_eof.store(true, std::sync::atomic::Ordering::Release);
        }
        slot
    }

    /// Drive a `MergeDriver` to exhaustion via the non-blocking `try_step`,
    /// returning the emitted record bytes in emission order.
    ///
    /// Tests validate ordering by inspecting the human-readable record bytes
    /// (e.g. `b"AAA-1"`, `b"BBB-2"`) — the keys themselves are consumed by
    /// the parser before the record is exposed, so we have no direct access to
    /// them at this layer. These slots are pre-populated and EOF-marked before
    /// the merge, so `try_step` never returns `Stalled` here.
    fn drain_merge_driver_bytes<K: RawSortKey + Default + Send + 'static>(
        mut driver: MergeDriver<K>,
    ) -> Vec<Vec<u8>> {
        let mut out = Vec::new();
        loop {
            match driver.try_step().expect("try_step") {
                MergeStep::Produced(bytes) => out.push(bytes.to_vec()),
                MergeStep::Done => break,
                MergeStep::Stalled => panic!("unexpected Stalled on pre-populated EOF slots"),
            }
        }
        out
    }

    #[test]
    fn from_slots_merges_three_pre_populated_slots() {
        // Three slots, each with sorted records by TestKey. The merged output
        // should be globally sorted.
        let file0: Vec<(TestKey, Vec<u8>)> = vec![
            (TestKey(1), b"AAA-1".to_vec()),
            (TestKey(4), b"AAA-4".to_vec()),
            (TestKey(7), b"AAA-7".to_vec()),
        ];
        let file1: Vec<(TestKey, Vec<u8>)> = vec![
            (TestKey(2), b"BBB-2".to_vec()),
            (TestKey(5), b"BBB-5".to_vec()),
            (TestKey(8), b"BBB-8".to_vec()),
        ];
        let file2: Vec<(TestKey, Vec<u8>)> = vec![
            (TestKey(3), b"CCC-3".to_vec()),
            (TestKey(6), b"CCC-6".to_vec()),
            (TestKey(9), b"CCC-9".to_vec()),
        ];

        let slots = vec![
            populated_slot(0, &file0, 1),
            populated_slot(1, &file1, 1),
            populated_slot(2, &file2, 1),
        ];

        let driver = MergeDriver::<TestKey>::from_slots(slots, MemorySources::Owned(Vec::new()), 9);
        let emitted_bytes = drain_merge_driver_bytes(driver);
        let expected: Vec<Vec<u8>> = vec![
            b"AAA-1".to_vec(),
            b"BBB-2".to_vec(),
            b"CCC-3".to_vec(),
            b"AAA-4".to_vec(),
            b"BBB-5".to_vec(),
            b"CCC-6".to_vec(),
            b"AAA-7".to_vec(),
            b"BBB-8".to_vec(),
            b"CCC-9".to_vec(),
        ];
        assert_eq!(emitted_bytes, expected, "merged sequence not globally sorted by key");
    }

    #[test]
    fn from_slots_handles_records_spanning_block_boundaries() {
        // Same data, but each slot's bytes are split across 4 "decompressed
        // blocks" so the slot parser's read_exact + advance_to_next_block
        // path is exercised. The key+len header itself may straddle a block
        // boundary depending on byte counts.
        let records: Vec<(TestKey, Vec<u8>)> = (0u64..20)
            .map(|i| (TestKey(i), format!("record-{i:02}-with-some-padding").into_bytes()))
            .collect();

        let slots = vec![populated_slot(0, &records, 4)];
        let driver =
            MergeDriver::<TestKey>::from_slots(slots, MemorySources::Owned(Vec::new()), 20);
        let emitted_bytes = drain_merge_driver_bytes(driver);
        let expected: Vec<Vec<u8>> = records.into_iter().map(|(_, b)| b).collect();
        assert_eq!(emitted_bytes, expected, "cross-block parse corrupted record sequence");
    }

    #[test]
    fn from_slots_empty_sources_drain_to_zero_records() {
        let empty_slot = StdArc::new(SortMergeSlot::new(
            0,
            BufReader::new(tempfile::tempfile().expect("tempfile")),
            crate::codec::SpillCodec::Bgzf,
        ));
        // Mark drained without inserting any blocks.
        empty_slot.queue_eof.store(true, std::sync::atomic::Ordering::Release);
        // Lazy priming: `from_slots` always yields a driver; the empty
        // determination happens on the first `try_step`, which returns `Done`.
        let driver = MergeDriver::<TestKey>::from_slots(
            vec![empty_slot],
            MemorySources::Owned(Vec::new()),
            0,
        );
        let emitted = drain_merge_driver_bytes(driver);
        assert!(emitted.is_empty(), "all-empty merge must emit zero records");
    }

    /// The core BUG #3 regression: a slot whose decompressed queue is empty
    /// and `!queue_eof` must make `try_step` return `Stalled` (yield), NOT
    /// block the caller. The old blocking consumer (`block_ready.wait()`)
    /// would park here forever at a single worker. After the producer fills
    /// the slot and marks EOF, the merge resumes and drains to completion.
    #[test]
    fn from_slots_try_step_stalls_on_empty_non_eof_slot_then_resumes() {
        let slot = StdArc::new(SortMergeSlot::new(
            0,
            BufReader::new(tempfile::tempfile().expect("tempfile")),
            crate::codec::SpillCodec::Bgzf,
        ));
        // Empty queue, NOT eof — the producer is "still feeding".
        let mut driver = MergeDriver::<TestKey>::from_slots(
            vec![StdArc::clone(&slot)],
            MemorySources::Owned(Vec::new()),
            0,
        );

        // Priming cannot read the slot's first record → Stalled, not a hang.
        assert!(
            matches!(driver.try_step().expect("try_step"), MergeStep::Stalled),
            "empty non-eof slot must Stall, not block"
        );
        // Repeated calls keep stalling (idempotent yield).
        assert!(matches!(driver.try_step().expect("try_step"), MergeStep::Stalled));

        // Producer fills the slot and marks EOF.
        {
            let bytes = serialize_records(&[(TestKey(1), b"only-1".to_vec())]);
            let mut dec = slot.decompressed.lock().unwrap();
            dec.push_back(bytes);
            slot.queue_eof.store(true, std::sync::atomic::Ordering::Release);
        }

        // The merge now resumes: emits the record, then is Done.
        match driver.try_step().expect("try_step") {
            MergeStep::Produced(bytes) => assert_eq!(bytes, b"only-1"),
            MergeStep::Stalled => panic!("should have resumed, not Stalled"),
            MergeStep::Done => panic!("should have produced the record, not Done"),
        }
        assert!(matches!(driver.try_step().expect("try_step"), MergeStep::Done));
    }

    /// `try_step` must also stall (not block) when a slot drains to empty
    /// **mid-merge** (after priming) while still `!queue_eof`, then resume
    /// once the producer pushes the next block and marks EOF.
    #[test]
    fn from_slots_try_step_stalls_mid_merge_then_resumes() {
        let slot = StdArc::new(SortMergeSlot::new(
            0,
            BufReader::new(tempfile::tempfile().expect("tempfile")),
            crate::codec::SpillCodec::Bgzf,
        ));
        // One record present, NOT eof (more "coming").
        {
            let bytes = serialize_records(&[(TestKey(1), b"rec-A".to_vec())]);
            slot.decompressed.lock().unwrap().push_back(bytes);
        }
        let mut driver = MergeDriver::<TestKey>::from_slots(
            vec![StdArc::clone(&slot)],
            MemorySources::Owned(Vec::new()),
            0,
        );

        // Prime + emit the first record.
        match driver.try_step().expect("try_step") {
            MergeStep::Produced(bytes) => assert_eq!(bytes, b"rec-A"),
            other => panic!("expected Produced(rec-A), got {other:?}"),
        }
        // Deferred refill finds the slot drained + not eof → Stalled.
        assert!(
            matches!(driver.try_step().expect("try_step"), MergeStep::Stalled),
            "drained mid-merge non-eof slot must Stall"
        );

        // Producer pushes the next record and marks EOF.
        {
            let bytes = serialize_records(&[(TestKey(2), b"rec-B".to_vec())]);
            let mut dec = slot.decompressed.lock().unwrap();
            dec.push_back(bytes);
            slot.queue_eof.store(true, std::sync::atomic::Ordering::Release);
        }

        match driver.try_step().expect("try_step") {
            MergeStep::Produced(bytes) => assert_eq!(bytes, b"rec-B"),
            other => panic!("expected Produced(rec-B), got {other:?}"),
        }
        assert!(matches!(driver.try_step().expect("try_step"), MergeStep::Done));
    }

    /// Resumable-framer regression: feed a single **non-embedded** record one
    /// byte per block, asserting `try_step` `Stalled`s (retaining
    /// `parser.pending`) until the whole record is available, then reassembles
    /// the exact bytes. With 1-byte blocks the framer takes the slow path and
    /// resumes a `WouldBlock` across the key prefix, the length prefix, AND the
    /// body — the partial-record path no other test covers byte-for-byte.
    #[test]
    fn from_slots_try_step_resumes_record_fed_one_byte_at_a_time() {
        let rec_body = b"multi-stage-resumable-record-body".to_vec();
        let serialized = serialize_records(&[(TestKey(7), rec_body.clone())]);
        assert!(serialized.len() > 12, "must span the 8-byte key + 4-byte len header");

        let slot = StdArc::new(SortMergeSlot::new(
            0,
            BufReader::new(tempfile::tempfile().expect("tempfile")),
            crate::codec::SpillCodec::Bgzf,
        ));
        slot.decompressed.lock().unwrap().push_back(vec![serialized[0]]);

        let mut driver = MergeDriver::<TestKey>::from_slots(
            vec![StdArc::clone(&slot)],
            MemorySources::Owned(Vec::new()),
            1,
        );

        // Each step before the last byte consumes one byte into the pending
        // record and stalls (queue empty, not EOF).
        for (i, &b) in serialized.iter().enumerate().skip(1) {
            assert!(
                matches!(driver.try_step().expect("try_step"), MergeStep::Stalled),
                "expected Stalled with only {i} of {} bytes available",
                serialized.len(),
            );
            slot.decompressed.lock().unwrap().push_back(vec![b]);
        }

        // The step that consumes the final byte completes the record.
        match driver.try_step().expect("try_step") {
            MergeStep::Produced(bytes) => assert_eq!(bytes, rec_body.as_slice()),
            other => panic!("expected Produced(reassembled), got {other:?}"),
        }
        // No more bytes and not yet EOF → Stalled; EOF → Done.
        assert!(matches!(driver.try_step().expect("try_step"), MergeStep::Stalled));
        slot.queue_eof.store(true, std::sync::atomic::Ordering::Release);
        assert!(matches!(driver.try_step().expect("try_step"), MergeStep::Done));
    }

    /// Same resumable-framer property for an **embedded**-key record split mid
    /// body across two blocks (the key lives in the body, so this exercises the
    /// embedded `extract_from_record` path resuming after a `WouldBlock`).
    #[test]
    fn from_slots_try_step_resumes_embedded_record_split_mid_body() {
        // Embedded record: [len:4][body], key = first 8 body bytes.
        let record = embedded_record(123, b"-embedded-tail-bytes");
        let serialized = serialize_embedded(std::slice::from_ref(&record));
        // Cut mid-body: 4-byte len + 10 body bytes in the first block.
        let split = 4 + 10;
        assert!(split < serialized.len() && split > 4 + 8, "split mid-body, past the key");
        let (head, tail) = serialized.split_at(split);

        let slot = StdArc::new(SortMergeSlot::new(
            0,
            BufReader::new(tempfile::tempfile().expect("tempfile")),
            crate::codec::SpillCodec::Bgzf,
        ));
        slot.decompressed.lock().unwrap().push_back(head.to_vec());

        let mut driver = MergeDriver::<TestEmbeddedKey>::from_slots(
            vec![StdArc::clone(&slot)],
            MemorySources::Owned(Vec::new()),
            1,
        );

        // Body incomplete in the first block → Stalled mid-record.
        assert!(matches!(driver.try_step().expect("try_step"), MergeStep::Stalled));

        // Deliver the rest + EOF; the record reassembles exactly.
        {
            let mut dec = slot.decompressed.lock().unwrap();
            dec.push_back(tail.to_vec());
            slot.queue_eof.store(true, std::sync::atomic::Ordering::Release);
        }
        match driver.try_step().expect("try_step") {
            MergeStep::Produced(bytes) => assert_eq!(bytes, record.as_slice()),
            other => panic!("expected Produced(reassembled embedded), got {other:?}"),
        }
        assert!(matches!(driver.try_step().expect("try_step"), MergeStep::Done));
    }

    #[test]
    fn from_slots_merges_slots_with_memory_chunks() {
        let slot_records: Vec<(TestKey, Vec<u8>)> = vec![
            (TestKey(1), b"slot-1".to_vec()),
            (TestKey(3), b"slot-3".to_vec()),
            (TestKey(5), b"slot-5".to_vec()),
        ];
        let memory_records: Vec<(TestKey, fgumi_raw_bam::RawRecord)> = vec![
            (TestKey(2), fgumi_raw_bam::RawRecord::from(b"mem-2".to_vec())),
            (TestKey(4), fgumi_raw_bam::RawRecord::from(b"mem-4".to_vec())),
            (TestKey(6), fgumi_raw_bam::RawRecord::from(b"mem-6".to_vec())),
        ];

        let slots = vec![populated_slot(0, &slot_records, 2)];
        let driver = MergeDriver::<TestKey>::from_slots(
            slots,
            MemorySources::Owned(vec![memory_records]),
            6,
        );
        let emitted_bytes = drain_merge_driver_bytes(driver);
        let expected: Vec<Vec<u8>> = vec![
            b"slot-1".to_vec(),
            b"mem-2".to_vec(),
            b"slot-3".to_vec(),
            b"mem-4".to_vec(),
            b"slot-5".to_vec(),
            b"mem-6".to_vec(),
        ];
        assert_eq!(emitted_bytes, expected, "slot+memory merge not globally sorted");
    }

    // (v4: `from_slots_bails_on_unfilled_gap_in_decompressed` removed.
    //  v3.1's `decompressed: Mutex<ReorderBuffer<Vec<u8>>>` is now a
    //  `Mutex<VecDeque<Vec<u8>>>` (FIFO). There is no concept of
    //  "ordinals" or "gaps" — the producer pushes in order under the
    //  reader lock, so a missing-block-in-the-middle state is no
    //  longer representable.)

    // ------------------------------------------------------------------------
    // Coverage for the EMBEDDED_IN_RECORD = true parser arm.
    //
    // All three production sort keys (RawCoordinateKey, RawQuerynameKey,
    // TemplateKey) set EMBEDDED_IN_RECORD = true, but TestKey above uses
    // the non-embedded format. This test pins the embedded-format slot path
    // with a synthetic embedded key (TestEmbeddedKey) whose value lives at
    // bytes 0..8 of the record itself.

    /// Embedded-key test type: serialized format is `[len(4)][record(len)]`
    /// (no separate key prefix). The key is the first 8 LE bytes of the record.
    #[derive(Default, PartialEq, Eq, PartialOrd, Ord, Clone, Copy, Debug)]
    struct TestEmbeddedKey(u64);

    impl RawSortKey for TestEmbeddedKey {
        const SERIALIZED_SIZE: Option<usize> = Some(0);
        const EMBEDDED_IN_RECORD: bool = true;
        fn extract(_bam: &[u8], _ctx: &crate::keys::SortContext) -> Self {
            unimplemented!("TestEmbeddedKey is only used in MergeDriver::from_slots tests")
        }
        fn extract_from_record(bam: &[u8]) -> Self {
            assert!(bam.len() >= 8, "embedded-key test record must have ≥8 bytes");
            let mut buf = [0u8; 8];
            buf.copy_from_slice(&bam[..8]);
            TestEmbeddedKey(u64::from_le_bytes(buf))
        }
        fn write_to<W: Write>(&self, _w: &mut W) -> std::io::Result<()> {
            // No-op: embedded keys are not written separately; they live
            // inside the record bytes.
            Ok(())
        }
        fn read_from<R: Read>(_r: &mut R) -> std::io::Result<Self> {
            unreachable!(
                "EMBEDDED_IN_RECORD = true keys take the slot_try_next_record \
                 embedded arm, which never calls read_from"
            )
        }
    }

    /// Build a record whose first 8 bytes encode the embedded sort key.
    fn embedded_record(key: u64, tail: &[u8]) -> Vec<u8> {
        let mut v = Vec::with_capacity(8 + tail.len());
        v.extend_from_slice(&key.to_le_bytes());
        v.extend_from_slice(tail);
        v
    }

    /// Serialize records in the embedded-key spill format: `[len(4)][record(len)]`.
    fn serialize_embedded(records: &[Vec<u8>]) -> Vec<u8> {
        let mut out = Vec::new();
        for rec in records {
            #[allow(clippy::cast_possible_truncation)]
            let len = rec.len() as u32;
            out.write_all(&len.to_le_bytes()).unwrap();
            out.write_all(rec).unwrap();
        }
        out
    }

    fn embedded_populated_slot(
        file_id: u32,
        records: &[Vec<u8>],
        n_blocks: usize,
    ) -> StdArc<SortMergeSlot> {
        let bytes = serialize_embedded(records);
        let block_size = bytes.len().div_ceil(n_blocks.max(1));
        let slot = StdArc::new(SortMergeSlot::new(
            file_id,
            BufReader::new(tempfile::tempfile().expect("tempfile")),
            crate::codec::SpillCodec::Bgzf,
        ));
        {
            let mut dec = slot.decompressed.lock().unwrap();
            for chunk in bytes.chunks(block_size) {
                dec.push_back(chunk.to_vec());
            }
            slot.queue_eof.store(true, std::sync::atomic::Ordering::Release);
        }
        slot
    }

    #[test]
    fn from_slots_embedded_key_path_merges_correctly() {
        // Three slots, each with sorted records by embedded key. The merged
        // output must be globally sorted. Exercises `slot_try_next_record`'s
        // EMBEDDED_IN_RECORD arm + the `extract_from_record` key-from-bytes
        // path.
        let file0 =
            vec![embedded_record(1, b"-A"), embedded_record(4, b"-A"), embedded_record(7, b"-A")];
        let file1 =
            vec![embedded_record(2, b"-B"), embedded_record(5, b"-B"), embedded_record(8, b"-B")];
        let file2 =
            vec![embedded_record(3, b"-C"), embedded_record(6, b"-C"), embedded_record(9, b"-C")];
        let slots = vec![
            embedded_populated_slot(0, &file0, 1),
            embedded_populated_slot(1, &file1, 1),
            embedded_populated_slot(2, &file2, 1),
        ];

        let driver =
            MergeDriver::<TestEmbeddedKey>::from_slots(slots, MemorySources::Owned(Vec::new()), 9);
        let emitted = drain_merge_driver_bytes(driver);

        // Verify keys are emitted in sorted order. Tail bytes prove the
        // identity of each source.
        let expected = vec![
            embedded_record(1, b"-A"),
            embedded_record(2, b"-B"),
            embedded_record(3, b"-C"),
            embedded_record(4, b"-A"),
            embedded_record(5, b"-B"),
            embedded_record(6, b"-C"),
            embedded_record(7, b"-A"),
            embedded_record(8, b"-B"),
            embedded_record(9, b"-C"),
        ];
        assert_eq!(emitted, expected);
    }

    #[test]
    fn from_slots_skips_empty_slots_mixed_with_populated() {
        // Verify a slot whose decompressed queue is empty (and reader marked EOF)
        // is correctly elided from the LoserTree. Without correct handling
        // the priming step would either panic or produce a phantom source.
        let populated = vec![
            (TestKey(1), b"one".to_vec()),
            (TestKey(2), b"two".to_vec()),
            (TestKey(3), b"three".to_vec()),
        ];

        let empty_slot = StdArc::new(SortMergeSlot::new(
            99,
            BufReader::new(tempfile::tempfile().expect("tempfile")),
            crate::codec::SpillCodec::Bgzf,
        ));
        empty_slot.queue_eof.store(true, std::sync::atomic::Ordering::Release);
        // Insert empty-slot first, then a populated slot. The driver should
        // skip the empty one during priming and only emit records from the
        // populated slot.
        let slots = vec![empty_slot, populated_slot(0, &populated, 1)];

        let driver = MergeDriver::<TestKey>::from_slots(slots, MemorySources::Owned(Vec::new()), 3);
        let emitted = drain_merge_driver_bytes(driver);
        assert_eq!(emitted, vec![b"one".to_vec(), b"two".to_vec(), b"three".to_vec()]);
    }

    #[test]
    fn from_slots_records_merged_equals_emitted_total() {
        // `records_merged()` (surfaced by `SortMerge` at completion) must equal
        // the number of records actually emitted across all sources.
        let f0: Vec<(TestKey, Vec<u8>)> =
            (0u64..15).map(|i| (TestKey(i * 2), format!("a-{i:02}").into_bytes())).collect();
        let f1: Vec<(TestKey, Vec<u8>)> =
            (0u64..15).map(|i| (TestKey(i * 2 + 1), format!("b-{i:02}").into_bytes())).collect();
        let total = (f0.len() + f1.len()) as u64;

        let slots = vec![populated_slot(0, &f0, 3), populated_slot(1, &f1, 2)];
        let mut driver =
            MergeDriver::<TestKey>::from_slots(slots, MemorySources::Owned(Vec::new()), total);

        let mut emitted = 0u64;
        loop {
            match driver.try_step().expect("try_step") {
                MergeStep::Produced(_) => emitted += 1,
                MergeStep::Done => break,
                MergeStep::Stalled => panic!("unexpected Stalled on pre-populated EOF slots"),
            }
        }
        assert_eq!(emitted, total, "should emit every input record");
        assert_eq!(
            driver.records_merged(),
            total,
            "records_merged must equal the emitted record count"
        );
    }
}

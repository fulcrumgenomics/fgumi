//! Raw-bytes sorting implementation for BAM files.
//!
//! This module implements high-performance BAM sorting using lazy record parsing.
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

use crate::bam_io::create_raw_bam_reader;
use crate::progress::ProgressTracker;
#[cfg(test)]
use crate::sam::SamTag;
use crate::sort::inline_buffer::{
    ProbeableBuffer, RecordBuffer, TemplateKey, TemplateRecordBuffer,
};
use crate::sort::keys::{QuerynameComparator, RawSortKey, SortOrder};
use crate::sort::memory_probe::{
    BufferProbeStats, ConsumerProbeStats, MergeProbe, SpillProbe, force_mi_collect, log_snapshot,
};
use crate::sort::pooled_chunk_writer::PooledChunkWriter;
use crate::sort::read_ahead::{RawReadAheadReader, RecordSource};
use crate::sort::worker_pool::SortWorkerPool;
use anyhow::Result;
use crossbeam_channel::{Receiver, Sender, bounded};
use log::{debug, info};
use noodles::sam::Header;
use noodles::sam::header::record::value::map::read_group::tag as rg_tag;
use noodles_bgzf::io::{
    MultithreadedWriter, Reader as BgzfReader, Writer as BgzfWriter, multithreaded_writer,
    writer::CompressionLevel,
};
use std::collections::HashMap;
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
            .find_string(&SamTag::RG)
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
}

/// Number of records to prefetch per chunk during merge.
/// Larger buffer reduces I/O latency impact during merge.
const MERGE_PREFETCH_SIZE: usize = 1024;

/// Maximum number of temp files before consolidation (like samtools).
/// When this limit is reached, oldest files are merged to reduce file count.
const DEFAULT_MAX_TEMP_FILES: usize = 64;

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
pub struct GenericKeyedChunkWriter<K: RawSortKey> {
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
    #[allow(clippy::cast_possible_truncation)]
    pub fn write_record(&mut self, key: &K, record: &[u8]) -> Result<()> {
        if !K::EMBEDDED_IN_RECORD {
            key.write_to(&mut self.writer)?;
        }
        self.writer.write_all(&(record.len() as u32).to_le_bytes())?;
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
fn read_exact_or_eof<R: Read>(reader: &mut R, buf: &mut [u8]) -> std::io::Result<bool> {
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
pub struct GenericKeyedChunkReader<K: RawSortKey + 'static> {
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

            // Check for gzip/BGZF magic bytes: 0x1f 0x8b
            let mut magic = [0u8; 2];
            let is_compressed = if buf_reader.read_exact(&mut magic).is_ok() {
                magic == [0x1f, 0x8b]
            } else {
                false
            };

            // Seek back to start
            if buf_reader.seek(SeekFrom::Start(0)).is_err() {
                let _ = tx
                    .send(Err(anyhow::anyhow!("Failed to seek in keyed chunk {}", path.display())));
                return;
            }

            // Read using appropriate decoder
            if is_compressed {
                let bgzf_reader = BgzfReader::new(buf_reader);
                Self::read_records(bgzf_reader, tx, buf_rx, concurrency_limit);
            } else {
                Self::read_records(buf_reader, tx, buf_rx, concurrency_limit);
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

/// Source for keyed chunks during merge (disk or in-memory).
enum ChunkSource<K: RawSortKey + Default + 'static> {
    /// Disk-based chunk with prefetching reader (legacy path with per-source threads).
    Disk(GenericKeyedChunkReader<K>),
    /// In-memory sorted records from the inline buffer.
    Memory { records: Vec<(K, fgumi_raw_bam::RawRecord)>, idx: usize },
    /// Pool-integrated disk source — workers read and decompress, main thread parses.
    /// The `source_id` maps to the `MainThreadChunkConsumer`'s per-source buffer.
    PoolDisk { source_id: usize },
}

impl<K: RawSortKey + Default + 'static> ChunkSource<K> {
    /// Fill `buf` with the next record's bytes and return the sort key,
    /// or `None` at EOF.
    ///
    /// For `PoolDisk` sources, `consumer` must be `Some`.
    fn next_record(
        &mut self,
        buf: &mut Vec<u8>,
        consumer: Option<&mut MainThreadChunkConsumer<K>>,
    ) -> Result<Option<K>> {
        match self {
            ChunkSource::Disk(reader) => reader.next_record(buf),
            ChunkSource::Memory { records, idx } => {
                if *idx < records.len() {
                    let (ref mut key, ref mut data) = records[*idx];
                    // Bridge: RawRecord wraps Vec<u8>; swap via the inner vec to avoid
                    // re-allocating. The caller's buf is a plain Vec<u8> (merge scratch).
                    // TODO: change merge scratch to RawRecord to eliminate this bridge.
                    std::mem::swap(buf, data.as_mut_vec());
                    let key = std::mem::take(key);
                    *idx += 1;
                    Ok(Some(key))
                } else {
                    Ok(None)
                }
            }
            ChunkSource::PoolDisk { source_id } => consumer
                .ok_or_else(|| {
                    anyhow::anyhow!(
                        "PoolDisk source (id {source_id}) requires a MainThreadChunkConsumer \
                         but none was provided — this is a bug in the sort pipeline"
                    )
                })?
                .next_record(*source_id, buf),
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
    files: Arc<Vec<crate::sort::worker_pool::Phase2FileState>>,
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
        files: Arc<Vec<crate::sort::worker_pool::Phase2FileState>>,
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
                    st.current_buf = data;
                    st.current_pos = 0;
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
struct ChunkNamer<'a> {
    temp_path: &'a Path,
    chunk_count: usize,
    merge_count: usize,
}

impl<'a> ChunkNamer<'a> {
    fn new(temp_path: &'a Path) -> Self {
        Self { temp_path, chunk_count: 0, merge_count: 0 }
    }

    /// Returns the next unique chunk file path.
    fn next_chunk_path(&mut self) -> PathBuf {
        let path = self.temp_path.join(format!("chunk_{:04}.keyed", self.chunk_count));
        self.chunk_count += 1;
        path
    }

    /// Returns the next unique merged file path.
    fn next_merged_path(&mut self) -> PathBuf {
        let path = self.temp_path.join(format!("merged_{:04}.keyed", self.merge_count));
        self.merge_count += 1;
        path
    }
}

/// A spill write that is finishing in the background.
///
/// Used for pipelining: the I/O thread continues writing while the main thread
/// reads the next batch. The chunk path is stored alongside the handle so it
/// can be pushed to `chunk_files` only after the write completes.
struct PendingSpill {
    handle: crate::sort::pooled_chunk_writer::SpillWriteHandle,
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
    /// Temporary directory for spill files.
    temp_dir: Option<PathBuf>,
    /// Number of threads for parallel operations.
    threads: usize,
    /// Compression level for output.
    output_compression: u32,
    /// Compression level for temporary chunk files (0 = uncompressed).
    temp_compression: u32,
    /// Whether to write BAM index alongside output (coordinate sort only).
    write_index: bool,
    /// Program record info (version, `command_line`) for @PG header.
    pg_info: Option<(String, String)>,
    /// Maximum temp files before consolidation (0 = unlimited).
    max_temp_files: usize,
    /// Cell barcode tag for template-coordinate sort (e.g., `[b'C', b'B']`).
    /// When `Some`, CB hash is included in sort key for single-cell data.
    cell_tag: Option<[u8; 2]>,
    /// Initial buffer capacity hint (bytes) for pre-allocation.
    ///
    /// Decoupled from `memory_limit` so that auto-detected limits can start with
    /// a modest allocation and let `Vec` grow on demand, while explicit limits
    /// pre-allocate the full budget upfront (preserving prior behavior).
    initial_capacity: Option<usize>,
    /// When true, wrap input in a `PrefetchReader` for async I/O.
    async_reader: bool,
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
            self.pool.set_phase(crate::sort::worker_pool::phase::LEGACY);
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
    /// Create a new raw external sorter with the given sort order.
    #[must_use]
    pub fn new(sort_order: SortOrder) -> Self {
        Self {
            sort_order,
            memory_limit: 512 * 1024 * 1024, // 512 MB default
            temp_dir: None,
            threads: 1,
            output_compression: 6,
            temp_compression: 1, // Default: fast compression
            write_index: false,
            pg_info: None,
            max_temp_files: DEFAULT_MAX_TEMP_FILES,
            cell_tag: None,
            initial_capacity: None,
            async_reader: false,
        }
    }

    /// Set the memory limit for in-memory sorting.
    #[must_use]
    pub fn memory_limit(mut self, limit: usize) -> Self {
        self.memory_limit = limit;
        self
    }

    /// Set the temporary directory for spill files.
    #[must_use]
    pub fn temp_dir(mut self, path: PathBuf) -> Self {
        self.temp_dir = Some(path);
        self
    }

    /// Set the number of threads.
    #[must_use]
    pub fn threads(mut self, threads: usize) -> Self {
        self.threads = threads;
        self
    }

    /// Set the output compression level.
    #[must_use]
    pub fn output_compression(mut self, level: u32) -> Self {
        self.output_compression = level;
        self
    }

    /// Set compression level for temporary chunk files.
    ///
    /// Level 0 disables compression (fastest, uses most disk space).
    /// Level 1 (default) provides fast compression with reasonable space savings.
    /// Higher levels provide better compression but are slower.
    #[must_use]
    pub fn temp_compression(mut self, level: u32) -> Self {
        self.temp_compression = level;
        self
    }

    /// Enable writing BAM index alongside output.
    ///
    /// Only valid for coordinate sort. When enabled, writes `<output>.bai`
    /// alongside the output BAM file. Uses single-threaded compression
    /// for accurate virtual position tracking.
    #[must_use]
    pub fn write_index(mut self, enabled: bool) -> Self {
        self.write_index = enabled;
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
    pub fn cell_tag(mut self, tag: [u8; 2]) -> Self {
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

    /// Enable async prefetch reader for input I/O.
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
            .num_threads(self.threads.max(1))
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
        pool: &std::sync::Arc<crate::sort::worker_pool::SortWorkerPool>,
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
        use crate::sort::loser_tree::LoserTree;

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

        info!(
            "Consolidating {} temp files into 1 (total was {})...",
            n_to_merge,
            n_to_merge + chunk_files.len()
        );

        // Create merged output file
        let merged_path = namer.next_merged_path();

        // Open readers with semaphore to cap concurrent I/O.
        let sem = make_reader_semaphore(self.threads);
        let mut readers: Vec<GenericKeyedChunkReader<K>> = files_to_merge
            .iter()
            .map(|p| GenericKeyedChunkReader::<K>::open(p, Some(Arc::clone(&sem))))
            .collect::<Result<Vec<_>>>()?;

        // Use pooled writer for parallel compression during consolidation.
        let mut writer = PooledChunkWriter::<K>::new(Arc::clone(pool), &merged_path)?;

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

        info!("Consolidation complete, {} temp files remain", chunk_files.len());

        Ok(())
    }

    /// Sort a BAM file using raw-bytes approach.
    ///
    /// # Errors
    ///
    /// Returns an error if reading, sorting, or writing the BAM file fails.
    pub fn sort(&self, input: &Path, output: &Path) -> Result<RawSortStats> {
        info!("Starting raw-bytes sort with order: {:?}", self.sort_order);
        info!("Memory limit: {} MB", self.memory_limit / (1024 * 1024));
        info!("Threads: {}", self.threads);

        // Shared worker pool for parallel BGZF compress/decompress across all phases
        let pool = Arc::new(SortWorkerPool::new(
            self.threads.max(1),
            self.temp_compression,
            self.output_compression,
        ));

        // Open input BAM and create record source
        // N+2 model: workers do ReadInputBlocks + DecompressInput,
        // main thread reads records directly from PooledInputStream.
        info!("Phase 1: Pool-integrated input reading ({} workers, N+2 model)", pool.num_workers());
        let (record_source, header) = {
            let (reader, header) = crate::bam_io::create_raw_bam_reader_pool_integrated(
                input,
                &pool,
                self.async_reader,
            )?;
            (RecordSource::direct(reader), header)
        };

        // Add @PG record if pg_info was provided
        let header = if let Some((ref version, ref command_line)) = self.pg_info {
            crate::header::add_pg_record(header, version, command_line)?
        } else {
            header
        };

        // Create temp directory
        let temp_dir = self.create_temp_dir()?;
        let temp_path = temp_dir.path();

        // Sort based on order
        match self.sort_order {
            SortOrder::Coordinate => {
                self.sort_coordinate(record_source, pool, &header, output, temp_path)
            }
            SortOrder::Queryname(comparator) => {
                self.sort_queryname(record_source, pool, &header, output, temp_path, comparator)
            }
            SortOrder::TemplateCoordinate => {
                self.sort_template_coordinate(record_source, pool, &header, output, temp_path)
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
        use crate::sort::inline_buffer::extract_coordinate_key_inline;
        use crate::sort::keys::{
            QuerynameComparator, RawCoordinateKey, RawQuerynameKey, RawQuerynameLexKey, RawSortKey,
            SortContext,
        };

        info!("Starting k-way merge of {} BAM files", inputs.len());

        let mut readers = Self::open_bam_prefetch_readers(inputs)?;
        let output_header = self.create_output_header(header);

        match self.sort_order {
            SortOrder::TemplateCoordinate => {
                let lib_lookup = LibraryLookup::from_header(header);
                let cell_tag = self.cell_tag;
                let hasher = cb_hasher();
                self.run_merge_loop(&mut readers, &output_header, output, |bam| {
                    extract_template_key_inline(bam, &lib_lookup, cell_tag.as_ref(), &hasher)
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
        use crate::sort::loser_tree::LoserTree;

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
            info!("Merge complete: 0 records merged");
            let writer = crate::bam_io::create_raw_bam_writer(
                output,
                output_header,
                self.threads,
                self.output_compression,
            )?;
            writer.finish()?;
            return Ok(0);
        }

        let mut tree = LoserTree::new(initial_keys);

        let mut writer = crate::bam_io::create_raw_bam_writer(
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

    /// Sort by coordinate order using optimized radix sort for large arrays.
    fn sort_coordinate(
        &self,
        record_source: RecordSource,
        pool: Arc<SortWorkerPool>,
        header: &Header,
        output: &Path,
        temp_path: &Path,
    ) -> Result<RawSortStats> {
        if self.write_index {
            self.sort_coordinate_with_index(record_source, pool, header, output, temp_path)
        } else {
            self.sort_coordinate_optimized(record_source, pool, header, output, temp_path)
        }
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
        temp_path: &Path,
    ) -> Result<RawSortStats> {
        use crate::sort::keys::RawCoordinateKey;

        let mut stats = RawSortStats::default();
        let mut timer = SortPhaseTimer::new();

        // Get number of references (unmapped reads map to nref)
        let nref = header.reference_sequences().len() as u32;

        // Estimate capacity from initial_capacity (not memory_limit) to avoid
        // huge upfront allocations when auto-detecting memory.
        let init_cap = self.effective_initial_capacity();
        // Per-record footprint: ~200 bytes BAM + 16 header + 24 ref = ~240 bytes
        let estimated_records = init_cap / 240;
        // Data bytes = init_cap minus ref overhead (24 bytes/record)
        let estimated_data_bytes = init_cap.saturating_sub(estimated_records * 24);

        let mut chunk_files: Vec<PathBuf> = Vec::new();
        let mut buffer = RecordBuffer::with_capacity(estimated_records, estimated_data_bytes, nref);
        let mut namer = ChunkNamer::new(temp_path);
        let mut pending_spill: Option<PendingSpill> = None;
        let rayon_pool = self.build_sort_rayon_pool()?;

        let progress = ProgressTracker::new("Read records").with_interval(1_000_000);
        info!("Phase 1: Reading and sorting chunks (inline buffer, keyed output)...");
        let mut probe = SpillProbe::new("phase1");

        for record in record_source.by_ref() {
            stats.total_records += 1;
            progress.log_if_needed(1);

            // Push directly to buffer - key extracted inline from raw bytes
            buffer.push_coordinate(record.as_ref())?;

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

                let chunk_path = namer.next_chunk_path();

                timer.time_sort(|| {
                    rayon_pool.install(|| buffer.par_sort());
                });

                // Write keyed temp file with parallel BGZF compression via worker pool.
                // Use start_finish() for pipelining: I/O continues in background
                // while we read the next batch.
                let handle = timer.time_spill_write(|| {
                    let mut writer =
                        PooledChunkWriter::<RawCoordinateKey>::new(Arc::clone(&pool), &chunk_path)?;
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

        if chunk_files.is_empty() {
            // All records fit in memory - no merge needed
            info!("All records fit in memory, performing in-memory sort");

            timer.time_sort(|| {
                rayon_pool.install(|| buffer.par_sort());
            });

            timer.time_write_output(|| {
                use crate::sort::pooled_bam_writer::PooledBamWriter;
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
            let memory_chunks: Vec<Vec<(RawCoordinateKey, fgumi_raw_bam::RawRecord)>> = if buffer
                .is_empty()
            {
                Vec::new()
            } else if self.threads > 1 {
                timer.time_sort(|| rayon_pool.install(|| buffer.par_sort_into_chunks(self.threads)))
            } else {
                timer.time_sort(|| {
                    rayon_pool.install(|| buffer.par_sort());
                });
                let chunk = buffer
                    .refs()
                    .iter()
                    .map(|r| {
                        let key = RawCoordinateKey { sort_key: r.sort_key };
                        (key, fgumi_raw_bam::RawRecord::from(buffer.get_record(r).to_vec()))
                    })
                    .collect();
                vec![chunk]
            };

            let n_memory = memory_chunks.iter().filter(|c| !c.is_empty()).count();
            info!(
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
        info!("Sort complete: {} records processed", stats.total_records);

        Ok(stats)
    }

    /// Coordinate sort with BAM index generation.
    ///
    /// Similar to `sort_coordinate_optimized` but uses `IndexingBamWriter` to
    /// build the BAI index incrementally during write. Uses single-threaded
    /// compression for accurate virtual position tracking.
    #[allow(clippy::cast_possible_truncation, clippy::too_many_lines)]
    fn sort_coordinate_with_index(
        &self,
        mut record_source: RecordSource,
        pool: Arc<SortWorkerPool>,
        header: &Header,
        output: &Path,
        temp_path: &Path,
    ) -> Result<RawSortStats> {
        use crate::bam_io::{create_indexing_bam_writer, write_bai_index};
        use crate::sort::keys::RawCoordinateKey;

        info!("Indexing enabled: will write BAM index alongside output");

        let mut stats = RawSortStats::default();
        let mut timer = SortPhaseTimer::new();

        let nref = header.reference_sequences().len() as u32;
        let init_cap = self.effective_initial_capacity();
        // Per-record footprint: ~200 bytes BAM + 16 header + 24 ref = ~240 bytes
        let estimated_records = init_cap / 240;
        let estimated_data_bytes = init_cap.saturating_sub(estimated_records * 24);

        let mut chunk_files: Vec<PathBuf> = Vec::new();
        let mut buffer = RecordBuffer::with_capacity(estimated_records, estimated_data_bytes, nref);
        let mut namer = ChunkNamer::new(temp_path);
        let mut pending_spill: Option<PendingSpill> = None;
        let rayon_pool = self.build_sort_rayon_pool()?;

        info!("Phase 1: Reading and sorting chunks (inline buffer, keyed output)...");
        let mut probe = SpillProbe::new("phase1");

        for record in record_source.by_ref() {
            stats.total_records += 1;
            buffer.push_coordinate(record.as_ref())?;

            if probe.should_sample_read(stats.total_records) {
                probe.log_mid_read(probe_stats(&buffer), Some(pool.phase1_queue_depths()));
            }

            if buffer.memory_usage() >= self.memory_limit {
                timer.end_read_span();
                let bstats = probe_stats(&buffer);
                let depths = Some(pool.phase1_queue_depths());
                probe.pre_spill(bstats, depths);

                // Wait for any previous spill to complete
                self.drain_pending_spill::<RawCoordinateKey>(
                    &mut pending_spill,
                    &mut chunk_files,
                    &mut stats,
                    &mut timer,
                    &mut namer,
                    &pool,
                )?;
                probe.post_drain(probe_stats(&buffer), Some(pool.phase1_queue_depths()));

                let chunk_path = namer.next_chunk_path();

                timer.time_sort(|| {
                    rayon_pool.install(|| buffer.par_sort());
                });

                let handle = timer.time_spill_write(|| {
                    let mut writer =
                        PooledChunkWriter::<RawCoordinateKey>::new(Arc::clone(&pool), &chunk_path)?;
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
        info!("Read {} records total", stats.total_records);
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

        let output_header = self.create_output_header(header);

        if chunk_files.is_empty() {
            // All records fit in memory - no merge needed
            info!("All records fit in memory, performing in-memory sort");

            timer.time_sort(|| {
                rayon_pool.install(|| buffer.par_sort());
            });

            timer.time_write_output(|| {
                let mut writer = create_indexing_bam_writer(
                    output,
                    &output_header,
                    self.output_compression,
                    self.threads,
                )?;

                for record_bytes in buffer.iter_sorted() {
                    writer.write_raw_record(record_bytes)?;
                }

                let index = writer.finish()?;

                let index_path = output.with_extension("bam.bai");
                write_bai_index(&index_path, &index)?;
                info!("Wrote BAM index: {}", index_path.display());
                Ok(())
            })?;
        } else {
            // Sort remaining records into separate sub-array chunks
            let memory_chunks: Vec<Vec<(RawCoordinateKey, fgumi_raw_bam::RawRecord)>> = if buffer
                .is_empty()
            {
                Vec::new()
            } else if self.threads > 1 {
                timer.time_sort(|| rayon_pool.install(|| buffer.par_sort_into_chunks(self.threads)))
            } else {
                timer.time_sort(|| {
                    rayon_pool.install(|| buffer.par_sort());
                });
                let chunk = buffer
                    .refs()
                    .iter()
                    .map(|r| {
                        let key = RawCoordinateKey { sort_key: r.sort_key };
                        (key, fgumi_raw_bam::RawRecord::from(buffer.get_record(r).to_vec()))
                    })
                    .collect();
                vec![chunk]
            };

            let n_memory = memory_chunks.iter().filter(|c| !c.is_empty()).count();
            info!(
                "Phase 2: Merging {} chunks with index generation...",
                chunk_files.len() + n_memory
            );

            timer.time_merge(|| {
                let index = self.merge_chunks_with_index::<RawCoordinateKey>(
                    &chunk_files,
                    memory_chunks,
                    header,
                    output,
                    stats.total_records,
                    &pool,
                )?;

                let index_path = output.with_extension("bam.bai");
                write_bai_index(&index_path, &index)?;
                info!("Wrote BAM index: {}", index_path.display());
                Ok(())
            })?;
        }

        stats.output_records = stats.total_records;
        if let Ok(pool) = Arc::try_unwrap(pool) {
            pool.shutdown();
        }
        timer.log_summary(self.threads);
        info!("Sort complete: {} records processed", stats.total_records);

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
        temp_path: &Path,
        comparator: QuerynameComparator,
    ) -> Result<RawSortStats> {
        use crate::sort::keys::{RawQuerynameKey, RawQuerynameLexKey};
        info!("Using queryname sort with {comparator} comparator");
        match comparator {
            QuerynameComparator::Lexicographic => self.sort_queryname_keyed::<RawQuerynameLexKey>(
                record_source,
                pool,
                header,
                output,
                temp_path,
            ),
            QuerynameComparator::Natural => self.sort_queryname_keyed::<RawQuerynameKey>(
                record_source,
                pool,
                header,
                output,
                temp_path,
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
        temp_path: &Path,
    ) -> Result<RawSortStats> {
        use crate::sort::keys::SortContext;

        let mut stats = RawSortStats::default();
        let mut timer = SortPhaseTimer::new();

        let ctx = SortContext::from_header(header);

        // Estimate capacity from initial_capacity to avoid huge upfront allocations.
        let init_cap = self.effective_initial_capacity();
        let estimated_records = init_cap / 300;

        let mut chunk_files: Vec<PathBuf> = Vec::new();
        let mut entries: Vec<(K, fgumi_raw_bam::RawRecord)> = Vec::with_capacity(estimated_records);
        let mut memory_used = 0usize;
        let mut namer = ChunkNamer::new(temp_path);
        let mut pending_spill: Option<PendingSpill> = None;
        let rayon_pool = self.build_sort_rayon_pool()?;

        let progress = ProgressTracker::new("Read records").with_interval(1_000_000);
        info!("Phase 1: Reading and sorting chunks (keyed output)...");
        let mut probe = SpillProbe::new("phase1");

        for record in record_source.by_ref() {
            stats.total_records += 1;
            progress.log_if_needed(1);

            // Extract key from raw bytes
            let bam_bytes = record.as_ref();
            let key = K::extract(bam_bytes, &ctx);

            // Estimate memory: record bytes + key overhead
            let record_size = bam_bytes.len() + 50; // approximate key size
            memory_used += record_size;

            entries.push((key, fgumi_raw_bam::RawRecord::from(bam_bytes.to_vec())));

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

                let chunk_path = namer.next_chunk_path();

                timer.time_sort(|| {
                    use rayon::prelude::*;
                    rayon_pool.install(|| entries.par_sort_unstable_by(|a, b| a.0.cmp(&b.0)));
                });

                // Write keyed temp file with parallel BGZF compression via worker pool.
                let handle = timer.time_spill_write(|| {
                    let mut writer = PooledChunkWriter::<K>::new(Arc::clone(&pool), &chunk_path)?;
                    for (key, record) in entries.drain(..) {
                        writer.write_record(&key, &record)?;
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

        if chunk_files.is_empty() {
            // All records fit in memory
            info!("All records fit in memory, performing in-memory sort");

            timer.time_sort(|| {
                use rayon::prelude::*;
                rayon_pool.install(|| entries.par_sort_unstable_by(|a, b| a.0.cmp(&b.0)));
            });

            timer.time_write_output(|| {
                use crate::sort::pooled_bam_writer::PooledBamWriter;
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
            // intermediate merge back into a single sorted buffer)
            let memory_chunks: Vec<Vec<(K, fgumi_raw_bam::RawRecord)>> = if entries.is_empty() {
                Vec::new()
            } else if self.threads > 1 {
                timer.time_sort(|| {
                    use rayon::prelude::*;
                    let chunk_size = entries.len().div_ceil(self.threads.max(1));
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

            let n_memory = memory_chunks.iter().filter(|c| !c.is_empty()).count();
            info!(
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
        info!("Sort complete: {} records processed", stats.total_records);

        Ok(stats)
    }

    /// Sort by template-coordinate order using inline buffer for reduced memory.
    ///
    /// Uses `TemplateRecordBuffer` which stores records in a single contiguous allocation
    /// with packed sort keys, eliminating per-record heap allocations for names.
    ///
    /// Writes keyed temp chunks that preserve pre-computed sort keys, enabling O(1)
    /// comparisons during merge (instead of expensive CIGAR/aux parsing).
    #[allow(clippy::too_many_lines)]
    fn sort_template_coordinate(
        &self,
        mut record_source: RecordSource,
        pool: Arc<SortWorkerPool>,
        header: &Header,
        output: &Path,
        temp_path: &Path,
    ) -> Result<RawSortStats> {
        let mut stats = RawSortStats::default();
        let mut timer = SortPhaseTimer::new();

        // Build library lookup ONCE before sorting for O(1) ordinal lookups
        let lib_lookup = LibraryLookup::from_header(header);
        let cb_hasher = cb_hasher();

        // Estimate capacity from initial_capacity to avoid huge upfront allocations.
        // Memory layout per record:
        //   - data: 48 bytes (inline header) + ~250 bytes (BAM record) = 298 bytes
        //   - refs: 56 bytes (TemplateKey + u64 offset + u32 len + u32 pad)
        //   - Total: ~354 bytes per record
        let init_cap = self.effective_initial_capacity();
        let bytes_per_record = 354;
        let estimated_records = init_cap / bytes_per_record;
        // Allocate ~86% for data, ~14% for refs (48/338 ≈ 14%)
        let estimated_data_bytes = init_cap * 86 / 100;

        let mut chunk_files: Vec<PathBuf> = Vec::new();
        let mut buffer =
            TemplateRecordBuffer::with_capacity(estimated_records, estimated_data_bytes);
        let mut namer = ChunkNamer::new(temp_path);
        let mut pending_spill: Option<PendingSpill> = None;
        let rayon_pool = self.build_sort_rayon_pool()?;

        let progress = ProgressTracker::new("Read records").with_interval(1_000_000);
        info!("Phase 1: Reading and sorting chunks (inline buffer)...");
        let mut probe = SpillProbe::new("phase1");

        for record in record_source.by_ref() {
            stats.total_records += 1;
            progress.log_if_needed(1);

            // Extract template key and push to buffer
            let bam_bytes = record.as_ref();
            let key = extract_template_key_inline(
                bam_bytes,
                &lib_lookup,
                self.cell_tag.as_ref(),
                &cb_hasher,
            );
            buffer.push(bam_bytes, key)?;

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
                self.drain_pending_spill::<TemplateKey>(
                    &mut pending_spill,
                    &mut chunk_files,
                    &mut stats,
                    &mut timer,
                    &mut namer,
                    &pool,
                )?;
                probe.post_drain(probe_stats(&buffer), Some(pool.phase1_queue_depths()));

                let chunk_path = namer.next_chunk_path();

                timer.time_sort(|| {
                    rayon_pool.install(|| buffer.par_sort());
                });

                // Write keyed chunk with parallel BGZF compression via worker pool.
                let handle = timer.time_spill_write(|| {
                    let mut writer =
                        PooledChunkWriter::<TemplateKey>::new(Arc::clone(&pool), &chunk_path)?;
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
        self.drain_pending_spill::<TemplateKey>(
            &mut pending_spill,
            &mut chunk_files,
            &mut stats,
            &mut timer,
            &mut namer,
            &pool,
        )?;
        probe.phase1_end(buffer.memory_usage() as u64);

        if chunk_files.is_empty() {
            // All records fit in memory
            info!("All records fit in memory, performing in-memory sort");

            timer.time_sort(|| {
                rayon_pool.install(|| buffer.par_sort());
            });

            timer.time_write_output(|| {
                use crate::sort::pooled_bam_writer::PooledBamWriter;
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
            let memory_chunks: Vec<Vec<(TemplateKey, fgumi_raw_bam::RawRecord)>> = if buffer
                .is_empty()
            {
                Vec::new()
            } else if self.threads > 1 {
                timer.time_sort(|| rayon_pool.install(|| buffer.par_sort_into_chunks(self.threads)))
            } else {
                timer.time_sort(|| {
                    rayon_pool.install(|| buffer.par_sort());
                });
                let chunk = buffer
                    .iter_sorted_keyed()
                    .map(|(k, r)| (k, fgumi_raw_bam::RawRecord::from(r.to_vec())))
                    .collect();
                vec![chunk]
            };

            let n_memory = memory_chunks.iter().filter(|c| !c.is_empty()).count();
            info!("Phase 2: Merging {} chunks...", chunk_files.len() + n_memory);

            // Merge using O(1) key comparisons
            timer.time_merge(|| {
                self.merge_chunks_keyed(
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
        info!("Sort complete: {} records processed", stats.total_records);

        Ok(stats)
    }

    /// Build chunk sources from disk files and in-memory chunks.
    ///
    /// When `pool_decompress` is true, disk sources become `PoolDisk` variants —
    /// workers read and decompress in the background, main thread parses records.
    /// No per-source threads are spawned.
    ///
    /// When `pool_decompress` is false, disk sources use `GenericKeyedChunkReader`
    /// with its own per-source background thread. Used by the indexing path, which
    /// does not yet support pool-integrated decompression.
    fn build_chunk_sources<K: RawSortKey + Default + 'static>(
        chunk_files: &[PathBuf],
        memory_chunks: Vec<Vec<(K, fgumi_raw_bam::RawRecord)>>,
        reader_concurrency: usize,
        pool_decompress: bool,
        pool: &Arc<SortWorkerPool>,
    ) -> Result<Vec<ChunkSource<K>>> {
        let num_disk = chunk_files.len();
        let num_memory = memory_chunks.iter().filter(|c| !c.is_empty()).count();
        let mut sources: Vec<ChunkSource<K>> = Vec::with_capacity(num_disk + num_memory);

        if pool_decompress && !chunk_files.is_empty() {
            // Pool-integrated path: install per-file Phase 2 state on the
            // pool, then create one PoolDisk source per file. Workers
            // cooperatively read+decompress all files via work-stealing.
            pool.set_phase2_files(chunk_files)?;

            for source_id in 0..num_disk {
                sources.push(ChunkSource::PoolDisk { source_id });
            }
        } else {
            // Legacy path: per-source reader threads (used by indexing path)
            let sem = make_reader_semaphore(reader_concurrency);
            for path in chunk_files {
                sources.push(ChunkSource::Disk(GenericKeyedChunkReader::<K>::open(
                    path,
                    Some(Arc::clone(&sem)),
                )?));
            }
        }

        for chunk in memory_chunks {
            if !chunk.is_empty() {
                sources.push(ChunkSource::Memory { records: chunk, idx: 0 });
            }
        }

        Ok(sources)
    }

    /// Merge keyed chunks using O(1) key comparisons (delegates to generic merge).
    ///
    /// `memory_chunks` is a list of sorted in-memory sub-arrays, each of which
    /// becomes a separate merge source.
    fn merge_chunks_keyed(
        &self,
        chunk_files: &[PathBuf],
        memory_chunks: Vec<Vec<(TemplateKey, fgumi_raw_bam::RawRecord)>>,
        header: &Header,
        output: &Path,
        total_records: u64,
        pool: &Arc<SortWorkerPool>,
    ) -> Result<u64> {
        self.merge_chunks_generic::<TemplateKey>(
            chunk_files,
            memory_chunks,
            header,
            output,
            total_records,
            pool,
        )
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
        memory_chunks: Vec<Vec<(K, fgumi_raw_bam::RawRecord)>>,
        header: &Header,
        output: &Path,
        total_records: u64,
        pool: &Arc<SortWorkerPool>,
    ) -> Result<u64> {
        use crate::sort::loser_tree::LoserTree;
        use crate::sort::pooled_bam_writer::PooledBamWriter;
        use crate::sort::worker_pool::phase;

        let reader_concurrency: usize = 1;
        let num_disk = chunk_files.len();

        if num_disk > 0 {
            info!(
                "Pool-integrated merge: {} disk sources, {} pool workers (N+2 model)",
                num_disk,
                pool.num_workers()
            );
        }

        let mut sources = Self::build_chunk_sources::<K>(
            chunk_files,
            memory_chunks,
            reader_concurrency,
            true,
            pool,
        )?;

        let num_sources = sources.len();
        info!("Merging from {num_sources} sources...");

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

        // Initialize loser tree with first record from each source
        let mut initial_keys: Vec<K> = Vec::with_capacity(sources.len());
        let mut records: Vec<Vec<u8>> = Vec::with_capacity(sources.len());
        let mut source_map: Vec<usize> = Vec::with_capacity(sources.len());

        for (idx, source) in sources.iter_mut().enumerate() {
            let mut record = Vec::new();
            if let Some(key) = source.next_record(&mut record, guard.consumer_mut())? {
                initial_keys.push(key);
                records.push(record);
                source_map.push(idx);
            }
        }

        if initial_keys.is_empty() {
            info!("Merge complete: 0 records merged");
            guard.deactivate();
            let writer = PooledBamWriter::new(Arc::clone(pool), output, &output_header)?;
            writer.finish()?;
            return Ok(0);
        }

        let mut tree = LoserTree::new(initial_keys);

        info!("Merge thread budget: {} pool workers + 1 I/O + 1 main (N+2)", pool.num_workers());
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
            let sample_this = debug_timing && records_merged.is_multiple_of(merge_sample_interval);

            if sample_this {
                let t0 = Instant::now();
                writer.write_raw_record(&records[winner])?;
                merge_write_secs += t0.elapsed().as_secs_f64();
            } else {
                writer.write_raw_record(&records[winner])?;
            }

            records_merged += 1;
            merge_progress.log_if_needed(1);

            if merge_probe.should_sample(records_merged) {
                let depths = pool.phase1_queue_depths();
                let consumer_stats = guard.consumer_mut().map(|c| c.probe_consumer_stats());
                merge_probe.log_mid_with_depths(depths, consumer_stats);
            }

            let src_idx = source_map[winner];

            if sample_this {
                let t0 = Instant::now();
                let next =
                    sources[src_idx].next_record(&mut records[winner], guard.consumer_mut())?;
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
                let next =
                    sources[src_idx].next_record(&mut records[winner], guard.consumer_mut())?;
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

    /// Test helper: merge keyed chunk files into a BAM.
    ///
    /// Exposes `merge_chunks_generic` for use in tests where the spill pipeline
    /// produces chunk files that need merging.
    ///
    /// # Errors
    ///
    /// Returns an error if merging fails.
    #[cfg(test)]
    pub fn merge_chunks_for_test<K: RawSortKey + Default + 'static>(
        &self,
        chunk_files: &[PathBuf],
        memory_chunks: Vec<Vec<(K, fgumi_raw_bam::RawRecord)>>,
        header: &Header,
        output: &Path,
        pool: &Arc<SortWorkerPool>,
    ) -> Result<u64> {
        self.merge_chunks_generic::<K>(chunk_files, memory_chunks, header, output, 0, pool)
    }

    /// Merge keyed chunks with BAM index generation.
    ///
    /// Similar to `merge_chunks_generic` but uses `IndexingBamWriter` to build
    /// the BAI index incrementally during the merge. Returns the generated index.
    fn merge_chunks_with_index<K: RawSortKey + Default + 'static>(
        &self,
        chunk_files: &[PathBuf],
        memory_chunks: Vec<Vec<(K, fgumi_raw_bam::RawRecord)>>,
        header: &Header,
        output: &Path,
        total_records: u64,
        pool: &Arc<SortWorkerPool>,
    ) -> Result<noodles::bam::bai::Index> {
        use crate::bam_io::create_indexing_bam_writer;
        use crate::sort::loser_tree::LoserTree;

        // Thread budget: writer gets all threads (merge is output-compression-bound),
        // readers use semaphore concurrency of 1. Same split as merge_chunks_generic.
        // TODO: Replace create_indexing_bam_writer with PooledIndexingBamWriter
        let writer_threads = self.threads;
        let reader_concurrency: usize = 1;

        // The indexing path uses a non-pooled writer and has no pool consumer set up,
        // so use GenericKeyedChunkReader (pool_decompress=false).
        // TODO: Replace create_indexing_bam_writer with PooledIndexingBamWriter to
        // enable pool-integrated decompression for the indexing path.
        let mut sources = Self::build_chunk_sources::<K>(
            chunk_files,
            memory_chunks,
            reader_concurrency,
            false,
            pool,
        )?;

        let num_sources = sources.len();
        info!(
            "Merge thread budget (indexing): {writer_threads} writer + {reader_concurrency} reader + 1 main = {} total",
            writer_threads + reader_concurrency + 1
        );
        info!("Merging from {num_sources} sources (indexing)...");

        let output_header = self.create_output_header(header);

        // Initialize loser tree with first record from each source
        let mut initial_keys: Vec<K> = Vec::with_capacity(sources.len());
        let mut records: Vec<Vec<u8>> = Vec::with_capacity(sources.len());
        let mut source_map: Vec<usize> = Vec::with_capacity(sources.len());

        for (idx, source) in sources.iter_mut().enumerate() {
            let mut record = Vec::new();
            if let Some(key) = source.next_record(&mut record, None)? {
                initial_keys.push(key);
                records.push(record);
                source_map.push(idx);
            }
        }

        if initial_keys.is_empty() {
            let writer = create_indexing_bam_writer(
                output,
                &output_header,
                self.output_compression,
                writer_threads,
            )?;
            let index = writer.finish()?;
            info!("Merge complete: 0 records merged");
            return Ok(index);
        }

        let mut tree = LoserTree::new(initial_keys);

        let mut writer = create_indexing_bam_writer(
            output,
            &output_header,
            self.output_compression,
            writer_threads,
        )?;
        let merge_progress = ProgressTracker::new("Merged records")
            .with_interval(1_000_000)
            .with_total(total_records);

        while tree.winner_is_active() {
            let winner = tree.winner();
            writer.write_raw_record(&records[winner])?;
            merge_progress.log_if_needed(1);

            let src_idx = source_map[winner];
            if let Some(key) = sources[src_idx].next_record(&mut records[winner], None)? {
                tree.replace_winner(key);
            } else {
                tree.remove_winner();
            }
        }

        let index = writer.finish()?;
        merge_progress.log_final();
        Ok(index)
    }

    /// Create output header with appropriate sort order tags.
    fn create_output_header(&self, header: &Header) -> Header {
        super::create_output_header(self.sort_order, header)
    }

    /// Create temporary directory for spill files.
    fn create_temp_dir(&self) -> Result<TempDir> {
        super::create_temp_dir(self.temp_dir.as_deref())
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
    cell_tag: Option<&[u8; 2]>,
    cb_hasher: &ahash::RandomState,
) -> TemplateKey {
    use crate::sort::bam_fields;
    use bam_fields::{flags, mate_unclipped_5prime, unclipped_5prime_raw};

    // Single-pass extraction of all aux tags (MI, RG, cell barcode, MC)
    let aux = bam_fields::extract_template_aux_tags(bam_bytes, cell_tag);
    let mi = aux.mi;
    let library = lib_lookup.ordinal_from_rg(aux.rg);
    let cb_hash = aux.cell.map_or(0u64, |cb_bytes| cb_hasher.hash_one(cb_bytes));

    // Extract fields from raw bytes
    let v = bam_fields::RawRecordView::new(bam_bytes);
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

        // Completely unmapped - sort to end
        let is_read2 = (flag & 0x80) != 0; // is_last_segment flag
        return TemplateKey::unmapped(name_hash, cb_hash, is_read2);
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

// Use shared SortStats from the parent module.
pub use super::SortStats as RawSortStats;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sam::builder::SamBuilder;
    use bstr::BString;
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
        assert!(sorter.temp_dir.is_none());
        assert_eq!(sorter.threads, 1);
        assert_eq!(sorter.output_compression, 6);
        assert_eq!(sorter.temp_compression, 1);
        assert!(!sorter.write_index);
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
            .write_index(true)
            .pg_info("1.0".to_string(), "fgumi sort".to_string())
            .max_temp_files(128);

        assert_eq!(sorter.memory_limit, 1024);
        assert_eq!(sorter.temp_dir, Some(PathBuf::from("/tmp/test")));
        assert_eq!(sorter.threads, 8);
        assert_eq!(sorter.output_compression, 9);
        assert_eq!(sorter.temp_compression, 3);
        assert!(sorter.write_index);
        assert_eq!(sorter.pg_info, Some(("1.0".to_string(), "fgumi sort".to_string())));
        assert_eq!(sorter.max_temp_files, 128);
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
            fgumi_raw_bam::RawRecordView::new(&bam).tags().find_string(&SamTag::RG),
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
            fgumi_raw_bam::RawRecordView::new(&bam).tags().find_string(&SamTag::RG),
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
        let sorter = RawExternalSorter::new(SortOrder::TemplateCoordinate).cell_tag(*SamTag::CB);
        assert_eq!(sorter.cell_tag, Some(*SamTag::CB));
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

        let key = extract_template_key_inline(&bam, &lib_lookup, Some(b"CB"), &test_cb_hasher());
        assert_ne!(key.cb_hash, 0, "CB present should produce non-zero cb_hash");
    }

    #[test]
    fn test_extract_template_key_cb_absent_has_zero_hash() {
        let header = Header::builder().build();
        let lib_lookup = LibraryLookup::from_header(&header);
        // No CB tag in aux data
        let bam = build_mapped_bam(0, 100, b"read1", &[]);

        let key = extract_template_key_inline(&bam, &lib_lookup, Some(b"CB"), &test_cb_hasher());
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
        let key1 = extract_template_key_inline(&bam1, &lib_lookup, Some(b"CB"), &test_cb_hasher());

        let aux2 = cb_aux(b"TGCATGCA");
        let bam2 = build_mapped_bam(0, 100, b"read1", &aux2);
        let key2 = extract_template_key_inline(&bam2, &lib_lookup, Some(b"CB"), &test_cb_hasher());

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

        let key1 = extract_template_key_inline(&bam, &lib_lookup, Some(b"CB"), &test_cb_hasher());
        let key2 = extract_template_key_inline(&bam, &lib_lookup, Some(b"CB"), &test_cb_hasher());
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

        let key = extract_template_key_inline(&bam, &lib_lookup, Some(b"CB"), &test_cb_hasher());
        assert_ne!(key.cb_hash, 0, "unmapped read with CB should have non-zero cb_hash");
        assert_eq!(key.primary, u64::MAX, "unmapped both-mates should have MAX primary");
    }

    // ========================================================================
    // Consolidation chunk naming tests
    // ========================================================================

    /// Count records in a BAM file by reading with the raw BAM reader.
    fn count_bam_records(path: &std::path::Path) -> u64 {
        use crate::sort::read_ahead::RawReadAheadReader;
        let (reader, _) = create_raw_bam_reader(path, 1).expect("failed to create raw BAM reader");
        RawReadAheadReader::new(reader).count() as u64
    }

    /// Verifies that sort with consolidation preserves all records.
    ///
    /// Uses a tiny memory limit and low `max_temp_files` to force many chunks and
    /// consolidation. Before the fix (chunk files named by `chunk_files.len()`),
    /// consolidation would drain entries from the vector, shrinking its length,
    /// causing new chunks to collide with existing non-consolidated chunk files.
    #[rstest::rstest]
    #[case::coordinate(SortOrder::Coordinate, false)]
    #[case::coordinate_with_index(SortOrder::Coordinate, true)]
    #[case::queryname(SortOrder::Queryname(QuerynameComparator::default()), false)]
    #[case::queryname_natural(SortOrder::Queryname(QuerynameComparator::Natural), false)]
    #[case::template_coordinate(SortOrder::TemplateCoordinate, false)]
    fn test_sort_with_consolidation_preserves_all_records(
        #[case] sort_order: SortOrder,
        #[case] write_index: bool,
    ) {
        use crate::sam::builder::SamBuilder;

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
            .temp_compression(0)
            .output_compression(0)
            .write_index(write_index)
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
        use crate::sam::builder::SamBuilder;

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
        use crate::sam::builder::SamBuilder;

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
            .temp_compression(0)
            .output_compression(0)
            .sort(&input, &output_st)
            .expect("sort should succeed");

        // Sort multi-threaded
        RawExternalSorter::new(sort_order)
            .memory_limit(16 * 1024)
            .threads(2)
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
        use crate::sam::builder::SamBuilder;

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
        use crate::sort::read_ahead::RawReadAheadReader;
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
        use crate::sort::read_ahead::RawReadAheadReader;
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
}

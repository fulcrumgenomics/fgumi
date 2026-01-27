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

use crate::bam_io::create_bam_reader;
use crate::sort::inline_buffer::{RecordBuffer, TemplateKey, TemplateRecordBuffer};
use crate::sort::keys::SortOrder;
use crate::sort::radix::{heap_make, heap_sift_down};
use crate::sort::read_ahead::ReadAheadReader;
use anyhow::{Context, Result};
use bstr::BString;
use crossbeam_channel::{Receiver, bounded};
use log::info;
use noodles::sam::Header;
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::header::tag as header_tag;
use noodles::sam::header::record::value::map::read_group::tag as rg_tag;
use noodles_bgzf::io::{
    MultithreadedWriter, Reader as BgzfReader, Writer as BgzfWriter, multithreaded_writer,
    writer::CompressionLevel,
};
use std::num::NonZero;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::thread::{self, JoinHandle};
use tempfile::TempDir;

// ============================================================================
// Library Lookup for Template-Coordinate Sort
// ============================================================================

/// Maps read group ID -> library ordinal for O(1) comparison.
///
/// Pre-computes ordinals by sorting library names alphabetically.
/// Empty/unknown library sorts first (ordinal 0).
struct LibraryLookup {
    /// RG ID -> library ordinal
    rg_to_ordinal: HashMap<Vec<u8>, u32>,
}

impl LibraryLookup {
    /// Build lookup from BAM header.
    fn from_header(header: &Header) -> Self {
        // Collect all unique library names from read groups
        let mut libraries: Vec<String> = header
            .read_groups()
            .iter()
            .filter_map(|(_, rg)| rg.other_fields().get(&rg_tag::LIBRARY).map(|lb| lb.to_string()))
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
                    .map(|lb| lb.to_string())
                    .unwrap_or_default();
                let ordinal = *lib_to_ordinal.get(&lib).unwrap_or(&0);
                (id.to_vec(), ordinal)
            })
            .collect();

        Self { rg_to_ordinal }
    }

    /// Get library ordinal for a record (from RG tag in aux data).
    #[inline]
    fn get_ordinal(&self, bam: &[u8]) -> u32 {
        find_rg_tag_raw(bam).and_then(|rg| self.rg_to_ordinal.get(rg)).copied().unwrap_or(0)
    }
}

/// Find RG tag value in BAM auxiliary data.
fn find_rg_tag_raw(bam: &[u8]) -> Option<&[u8]> {
    let l_read_name = bam[8] as usize;
    let n_cigar_op = u16::from_le_bytes([bam[12], bam[13]]) as usize;
    let l_seq = u32::from_le_bytes([bam[16], bam[17], bam[18], bam[19]]) as usize;

    let aux_offset = 32 + l_read_name + n_cigar_op * 4 + l_seq.div_ceil(2) + l_seq;

    if aux_offset >= bam.len() {
        return None;
    }

    find_string_tag_in_aux(&bam[aux_offset..], *b"RG")
}

/// Find a string tag in auxiliary data, returning the value without type prefix.
fn find_string_tag_in_aux(aux: &[u8], tag: [u8; 2]) -> Option<&[u8]> {
    let mut pos = 0;
    while pos + 3 <= aux.len() {
        let tag_bytes = &aux[pos..pos + 2];
        let val_type = aux[pos + 2];

        if tag_bytes == tag && val_type == b'Z' {
            // Found the tag, find null terminator
            let start = pos + 3;
            let end = aux[start..].iter().position(|&b| b == 0)?;
            return Some(&aux[start..start + end]);
        }

        // Skip to next tag based on type
        pos += 3;
        match val_type {
            b'A' | b'c' | b'C' => pos += 1,
            b's' | b'S' => pos += 2,
            b'i' | b'I' | b'f' => pos += 4,
            b'Z' | b'H' => {
                // Null-terminated string
                while pos < aux.len() && aux[pos] != 0 {
                    pos += 1;
                }
                pos += 1; // Skip null
            }
            b'B' => {
                // Array: subtype (1 byte) + count (4 bytes) + data
                if pos + 5 > aux.len() {
                    return None;
                }
                let subtype = aux[pos];
                let count =
                    u32::from_le_bytes([aux[pos + 1], aux[pos + 2], aux[pos + 3], aux[pos + 4]])
                        as usize;
                pos += 5;
                let elem_size = match subtype {
                    b'c' | b'C' => 1,
                    b's' | b'S' => 2,
                    b'i' | b'I' | b'f' => 4,
                    _ => return None,
                };
                pos += count * elem_size;
            }
            _ => return None, // Unknown type
        }
    }
    None
}

/// Number of records to prefetch per chunk during merge.
/// Larger buffer reduces I/O latency impact during merge.
const MERGE_PREFETCH_SIZE: usize = 1024;

// ============================================================================
// Generic Keyed Temp File I/O (works with any RawSortKey)
// ============================================================================
//
// Stores pre-computed sort keys alongside each record for O(1) merge comparisons.
// Format: [key: serialized][len: 4 bytes][record: len bytes] per record

use crate::sort::keys::RawSortKey;
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
            // Single-threaded BGZF
            ChunkWriterInner::SingleThreaded(BgzfWriter::new(buf))
        };

        Ok(Self { writer, _marker: PhantomData })
    }

    /// Write a keyed record.
    #[inline]
    pub fn write_record(&mut self, key: &K, record: &[u8]) -> Result<()> {
        key.write_to(&mut self.writer)?;
        self.writer.write_all(&(record.len() as u32).to_le_bytes())?;
        self.writer.write_all(record)?;
        Ok(())
    }

    /// Finish writing and flush.
    pub fn finish(self) -> Result<()> {
        self.writer.finish()
    }
}

/// Generic reader for keyed temp chunks with background prefetching.
///
/// Works with any type implementing `RawSortKey`.
/// Auto-detects BGZF compression via magic bytes.
pub struct GenericKeyedChunkReader<K: RawSortKey + 'static> {
    receiver: Receiver<Option<(K, Vec<u8>)>>,
    _handle: JoinHandle<()>,
}

impl<K: RawSortKey + 'static> GenericKeyedChunkReader<K> {
    /// Open a keyed chunk file for reading with background prefetching.
    /// Auto-detects BGZF/gzip compression via magic bytes (0x1f 0x8b).
    pub fn open(path: &Path) -> Result<Self> {
        let (tx, rx) = bounded(MERGE_PREFETCH_SIZE);
        let path = path.to_path_buf();

        let handle = thread::spawn(move || {
            let file = match std::fs::File::open(&path) {
                Ok(f) => f,
                Err(e) => {
                    log::error!("Failed to open keyed chunk {}: {}", path.display(), e);
                    let _ = tx.send(None);
                    return;
                }
            };
            let mut buf_reader = BufReader::with_capacity(256 * 1024, file);

            // Check for gzip/BGZF magic bytes: 0x1f 0x8b
            let mut magic = [0u8; 2];
            let is_compressed = if buf_reader.read_exact(&mut magic).is_ok() {
                magic == [0x1f, 0x8b]
            } else {
                false
            };

            // Seek back to start
            if buf_reader.seek(SeekFrom::Start(0)).is_err() {
                log::error!("Failed to seek in keyed chunk {}", path.display());
                let _ = tx.send(None);
                return;
            }

            // Read using appropriate decoder
            if is_compressed {
                let bgzf_reader = BgzfReader::new(buf_reader);
                Self::read_records(bgzf_reader, tx);
            } else {
                Self::read_records(buf_reader, tx);
            }
        });

        Ok(Self { receiver: rx, _handle: handle })
    }

    /// Read records from a reader and send them through the channel.
    fn read_records<R: Read>(mut reader: R, tx: crossbeam_channel::Sender<Option<(K, Vec<u8>)>>) {
        loop {
            // Read key using the trait method
            let key = match K::read_from(&mut reader) {
                Ok(k) => k,
                Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => {
                    // Normal EOF
                    let _ = tx.send(None);
                    break;
                }
                Err(e) => {
                    log::error!("Error reading keyed chunk key: {}", e);
                    let _ = tx.send(None);
                    break;
                }
            };

            // Read record length
            let mut len_buf = [0u8; 4];
            if reader.read_exact(&mut len_buf).is_err() {
                log::error!("Error reading keyed chunk length");
                let _ = tx.send(None);
                break;
            }
            let len = u32::from_le_bytes(len_buf) as usize;

            // Read record data
            let mut record = vec![0u8; len];
            if reader.read_exact(&mut record).is_err() {
                log::error!("Error reading keyed chunk record");
                let _ = tx.send(None);
                break;
            }

            if tx.send(Some((key, record))).is_err() {
                break; // Receiver dropped
            }
        }
    }

    /// Read the next keyed record from the prefetch buffer.
    pub fn next_record(&mut self) -> Option<(K, Vec<u8>)> {
        match self.receiver.recv() {
            Ok(Some(entry)) => Some(entry),
            Ok(None) | Err(_) => None,
        }
    }
}

// ============================================================================
// Specialized Keyed Temp File I/O for Template-Coordinate Sort
// ============================================================================
//
// Stores pre-computed TemplateKey alongside each record so merge comparisons
// are O(1) integer comparisons instead of O(CIGAR_len + aux_scan) parsing.
//
// Format: [key: 32 bytes][len: 4 bytes][record: len bytes] per record

/// Writer for keyed temp chunks with pre-computed sort keys.
/// Supports optional BGZF compression for reduced disk usage.
struct KeyedChunkWriter {
    writer: ChunkWriterInner,
}

impl KeyedChunkWriter {
    /// Create a new keyed chunk writer with optional compression.
    ///
    /// - `compression_level` 0 = uncompressed (fastest, uses most disk).
    /// - `compression_level` > 0 = BGZF compression at specified level.
    /// - `threads` > 1 enables multi-threaded compression.
    fn create(path: &Path, compression_level: u32, threads: usize) -> Result<Self> {
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
            // Single-threaded BGZF
            ChunkWriterInner::SingleThreaded(BgzfWriter::new(buf))
        };

        Ok(Self { writer })
    }

    /// Write a keyed record.
    #[inline]
    fn write_record(&mut self, key: &TemplateKey, record: &[u8]) -> Result<()> {
        self.writer.write_all(&key.to_bytes())?;
        self.writer.write_all(&(record.len() as u32).to_le_bytes())?;
        self.writer.write_all(record)?;
        Ok(())
    }

    /// Finish writing and flush.
    fn finish(self) -> Result<()> {
        self.writer.finish()
    }
}

/// Reader for keyed temp chunks - returns pre-computed keys with records.
/// Auto-detects BGZF compression via magic bytes.
struct KeyedChunkReader {
    receiver: Receiver<Option<(TemplateKey, Vec<u8>)>>,
    _handle: JoinHandle<()>,
}

impl KeyedChunkReader {
    /// Open a keyed chunk file for reading with background prefetching.
    /// Auto-detects BGZF/gzip compression via magic bytes (0x1f 0x8b).
    fn open(path: &Path) -> Result<Self> {
        let (tx, rx) = bounded(MERGE_PREFETCH_SIZE);
        let path = path.to_path_buf();

        let handle = thread::spawn(move || {
            let file = match std::fs::File::open(&path) {
                Ok(f) => f,
                Err(e) => {
                    log::error!("Failed to open keyed chunk {}: {}", path.display(), e);
                    let _ = tx.send(None);
                    return;
                }
            };
            let mut buf_reader = BufReader::with_capacity(256 * 1024, file);

            // Check for gzip/BGZF magic bytes: 0x1f 0x8b
            let mut magic = [0u8; 2];
            let is_compressed = if buf_reader.read_exact(&mut magic).is_ok() {
                magic == [0x1f, 0x8b]
            } else {
                false
            };

            // Seek back to start
            if buf_reader.seek(SeekFrom::Start(0)).is_err() {
                log::error!("Failed to seek in keyed chunk {}", path.display());
                let _ = tx.send(None);
                return;
            }

            // Read using appropriate decoder
            if is_compressed {
                let bgzf_reader = BgzfReader::new(buf_reader);
                Self::read_records(bgzf_reader, tx);
            } else {
                Self::read_records(buf_reader, tx);
            }
        });

        Ok(Self { receiver: rx, _handle: handle })
    }

    /// Read records from a reader and send them through the channel.
    fn read_records<R: Read>(
        mut reader: R,
        tx: crossbeam_channel::Sender<Option<(TemplateKey, Vec<u8>)>>,
    ) {
        let mut key_buf = [0u8; 32];
        let mut len_buf = [0u8; 4];

        loop {
            // Read key
            match reader.read_exact(&mut key_buf) {
                Ok(()) => {}
                Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => {
                    // Normal EOF
                    let _ = tx.send(None);
                    break;
                }
                Err(e) => {
                    log::error!("Error reading keyed chunk key: {}", e);
                    let _ = tx.send(None);
                    break;
                }
            }

            // Read length
            if reader.read_exact(&mut len_buf).is_err() {
                log::error!("Error reading keyed chunk length");
                let _ = tx.send(None);
                break;
            }
            let len = u32::from_le_bytes(len_buf) as usize;

            // Read record
            let mut record = vec![0u8; len];
            if reader.read_exact(&mut record).is_err() {
                log::error!("Error reading keyed chunk record");
                let _ = tx.send(None);
                break;
            }

            let key = TemplateKey::from_bytes(&key_buf);
            if tx.send(Some((key, record))).is_err() {
                break; // Receiver dropped
            }
        }
    }

    /// Read the next keyed record from the prefetch buffer.
    fn next_record(&mut self) -> Option<(TemplateKey, Vec<u8>)> {
        match self.receiver.recv() {
            Ok(Some(entry)) => Some(entry),
            Ok(None) | Err(_) => None,
        }
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

    /// Sort a BAM file using raw-bytes approach.
    pub fn sort(&self, input: &Path, output: &Path) -> Result<RawSortStats> {
        info!("Starting raw-bytes sort with order: {:?}", self.sort_order);
        info!("Memory limit: {} MB", self.memory_limit / (1024 * 1024));
        info!("Threads: {}", self.threads);

        // Open input BAM with lazy record reading
        let (reader, header) = create_bam_reader(input, self.threads)?;

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
            SortOrder::Coordinate => self.sort_coordinate(reader, &header, output, temp_path),
            SortOrder::Queryname => self.sort_queryname(reader, &header, output, temp_path),
            SortOrder::TemplateCoordinate => {
                self.sort_template_coordinate(reader, &header, output, temp_path)
            }
        }
    }

    /// Sort by coordinate order using optimized radix sort for large arrays.
    fn sort_coordinate(
        &self,
        reader: crate::bam_io::BamReaderAuto,
        header: &Header,
        output: &Path,
        temp_path: &Path,
    ) -> Result<RawSortStats> {
        if self.write_index {
            self.sort_coordinate_with_index(reader, header, output, temp_path)
        } else {
            self.sort_coordinate_optimized(reader, header, output, temp_path)
        }
    }

    /// Optimized coordinate sort using inline buffer for reduced memory overhead.
    ///
    /// Uses `RecordBuffer` which stores records in a single contiguous allocation
    /// with pre-computed sort keys, eliminating per-record heap allocations.
    fn sort_coordinate_optimized(
        &self,
        reader: crate::bam_io::BamReaderAuto,
        header: &Header,
        output: &Path,
        temp_path: &Path,
    ) -> Result<RawSortStats> {
        use crate::sort::keys::RawCoordinateKey;

        let mut stats = RawSortStats::default();

        // Get number of references (unmapped reads map to nref)
        let nref = header.reference_sequences().len() as u32;

        // Estimate capacity: ~200 bytes per record average
        let estimated_records = self.memory_limit / 200;

        let mut chunk_files: Vec<PathBuf> = Vec::new();
        let mut buffer = RecordBuffer::with_capacity(estimated_records, self.memory_limit, nref);

        let read_ahead = ReadAheadReader::new(reader);

        info!("Phase 1: Reading and sorting chunks (inline buffer, keyed output)...");

        for record in read_ahead {
            stats.total_records += 1;

            // Push directly to buffer - key extracted inline from raw bytes
            buffer.push_coordinate(record.as_ref());

            // Check memory usage
            if buffer.memory_usage() >= self.memory_limit {
                let chunk_path = temp_path.join(format!("chunk_{:04}.keyed", chunk_files.len()));

                // Sort in place using parallel sort for large buffers
                if self.threads > 1 {
                    buffer.par_sort();
                } else {
                    buffer.sort();
                }

                // Write keyed temp file (stores sort key with each record for O(1) merge)
                let mut writer = GenericKeyedChunkWriter::<RawCoordinateKey>::create(
                    &chunk_path,
                    self.temp_compression,
                    self.threads,
                )?;
                for r in buffer.refs() {
                    let key = RawCoordinateKey { sort_key: r.sort_key, name_hash: r.name_hash };
                    let record_bytes = buffer.get_record(r);
                    writer.write_record(&key, record_bytes)?;
                }
                writer.finish()?;

                stats.chunks_written += 1;
                chunk_files.push(chunk_path);
                buffer.clear();
            }
        }

        info!("Read {} records total", stats.total_records);

        if chunk_files.is_empty() {
            // All records fit in memory - no merge needed
            info!("All records fit in memory, performing in-memory sort");

            if self.threads > 1 {
                buffer.par_sort();
            } else {
                buffer.sort();
            }

            let output_header = self.create_output_header(header);
            let mut writer = crate::bam_io::create_raw_bam_writer(
                output,
                &output_header,
                self.threads,
                self.output_compression,
            )?;

            for record_bytes in buffer.iter_sorted() {
                writer.write_raw_record(record_bytes)?;
            }
            writer.finish()?;
        } else {
            // Sort remaining records and prepare for merge
            if self.threads > 1 {
                buffer.par_sort();
            } else {
                buffer.sort();
            }

            // Extract keyed pairs from the in-memory buffer
            let memory_keyed: Vec<(RawCoordinateKey, Vec<u8>)> = if buffer.is_empty() {
                Vec::new()
            } else {
                buffer
                    .refs()
                    .iter()
                    .map(|r| {
                        let key = RawCoordinateKey { sort_key: r.sort_key, name_hash: r.name_hash };
                        let record_bytes = buffer.get_record(r).to_vec();
                        (key, record_bytes)
                    })
                    .collect()
            };

            info!(
                "Phase 2: Merging {} chunks (keyed O(1) comparisons)...",
                chunk_files.len() + usize::from(!memory_keyed.is_empty())
            );

            // Merge disk chunks + in-memory chunks using O(1) key comparisons
            self.merge_chunks_generic::<RawCoordinateKey>(
                &chunk_files,
                memory_keyed,
                header,
                output,
            )?;
        }

        stats.output_records = stats.total_records;
        info!("Sort complete: {} records processed", stats.total_records);

        Ok(stats)
    }

    /// Coordinate sort with BAM index generation.
    ///
    /// Similar to `sort_coordinate_optimized` but uses `IndexingBamWriter` to
    /// build the BAI index incrementally during write. Uses single-threaded
    /// compression for accurate virtual position tracking.
    fn sort_coordinate_with_index(
        &self,
        reader: crate::bam_io::BamReaderAuto,
        header: &Header,
        output: &Path,
        temp_path: &Path,
    ) -> Result<RawSortStats> {
        use crate::bam_io::{create_indexing_bam_writer, write_bai_index};
        use crate::sort::keys::RawCoordinateKey;

        info!("Indexing enabled: will write BAM index alongside output");

        let mut stats = RawSortStats::default();
        let nref = header.reference_sequences().len() as u32;
        let estimated_records = self.memory_limit / 200;

        let mut chunk_files: Vec<PathBuf> = Vec::new();
        let mut buffer = RecordBuffer::with_capacity(estimated_records, self.memory_limit, nref);
        let read_ahead = ReadAheadReader::new(reader);

        info!("Phase 1: Reading and sorting chunks (inline buffer, keyed output)...");

        for record in read_ahead {
            stats.total_records += 1;
            buffer.push_coordinate(record.as_ref());

            if buffer.memory_usage() >= self.memory_limit {
                let chunk_path = temp_path.join(format!("chunk_{:04}.keyed", chunk_files.len()));

                if self.threads > 1 {
                    buffer.par_sort();
                } else {
                    buffer.sort();
                }

                let mut writer = GenericKeyedChunkWriter::<RawCoordinateKey>::create(
                    &chunk_path,
                    self.temp_compression,
                    self.threads,
                )?;
                for r in buffer.refs() {
                    let key = RawCoordinateKey { sort_key: r.sort_key, name_hash: r.name_hash };
                    let record_bytes = buffer.get_record(r);
                    writer.write_record(&key, record_bytes)?;
                }
                writer.finish()?;

                stats.chunks_written += 1;
                chunk_files.push(chunk_path);
                buffer.clear();
            }
        }

        info!("Read {} records total", stats.total_records);

        let output_header = self.create_output_header(header);

        if chunk_files.is_empty() {
            // All records fit in memory - no merge needed
            info!("All records fit in memory, performing in-memory sort");

            if self.threads > 1 {
                buffer.par_sort();
            } else {
                buffer.sort();
            }

            // Use IndexingBamWriter for single-pass index generation
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

            // Write the index file
            let index_path = output.with_extension("bam.bai");
            write_bai_index(&index_path, &index)?;
            info!("Wrote BAM index: {}", index_path.display());
        } else {
            // Sort remaining records and prepare for merge
            if self.threads > 1 {
                buffer.par_sort();
            } else {
                buffer.sort();
            }

            let memory_keyed: Vec<(RawCoordinateKey, Vec<u8>)> = if buffer.is_empty() {
                Vec::new()
            } else {
                buffer
                    .refs()
                    .iter()
                    .map(|r| {
                        let key = RawCoordinateKey { sort_key: r.sort_key, name_hash: r.name_hash };
                        let record_bytes = buffer.get_record(r).to_vec();
                        (key, record_bytes)
                    })
                    .collect()
            };

            info!(
                "Phase 2: Merging {} chunks with index generation...",
                chunk_files.len() + usize::from(!memory_keyed.is_empty())
            );

            // Merge with index generation
            let index = self.merge_chunks_with_index::<RawCoordinateKey>(
                &chunk_files,
                memory_keyed,
                header,
                output,
            )?;

            // Write the index file
            let index_path = output.with_extension("bam.bai");
            write_bai_index(&index_path, &index)?;
            info!("Wrote BAM index: {}", index_path.display());
        }

        stats.output_records = stats.total_records;
        info!("Sort complete: {} records processed", stats.total_records);

        Ok(stats)
    }

    /// Sort by queryname order using keyed temp files for O(1) merge comparisons.
    ///
    /// Uses `RawQuerynameKey` which stores the full read name for natural ordering,
    /// enabling efficient merge without re-parsing BAM records.
    fn sort_queryname(
        &self,
        reader: crate::bam_io::BamReaderAuto,
        header: &Header,
        output: &Path,
        temp_path: &Path,
    ) -> Result<RawSortStats> {
        use crate::sort::keys::{RawQuerynameKey, SortContext};

        let mut stats = RawSortStats::default();
        let ctx = SortContext::from_header(header);

        // Estimate capacity: ~250 bytes per record + ~50 bytes for name + flags
        let estimated_records = self.memory_limit / 300;

        let mut chunk_files: Vec<PathBuf> = Vec::new();
        let mut entries: Vec<(RawQuerynameKey, Vec<u8>)> = Vec::with_capacity(estimated_records);
        let mut memory_used = 0usize;

        let read_ahead = ReadAheadReader::new(reader);

        info!("Phase 1: Reading and sorting chunks (keyed output)...");

        for record in read_ahead {
            stats.total_records += 1;

            // Extract key from raw bytes
            let bam_bytes = record.as_ref();
            let key = RawQuerynameKey::extract(bam_bytes, &ctx);

            // Estimate memory: key (name + 2 bytes flags) + record bytes
            let record_size = bam_bytes.len() + key.name.len() + 2;
            memory_used += record_size;

            entries.push((key, bam_bytes.to_vec()));

            // Check if we need to spill to disk
            if memory_used >= self.memory_limit {
                let chunk_path = temp_path.join(format!("chunk_{:04}.keyed", chunk_files.len()));

                // Sort the entries
                if self.threads > 1 {
                    use rayon::prelude::*;
                    entries.par_sort_unstable_by(|a, b| a.0.cmp(&b.0));
                } else {
                    entries.sort_unstable_by(|a, b| a.0.cmp(&b.0));
                }

                // Write keyed temp file
                let mut writer = GenericKeyedChunkWriter::<RawQuerynameKey>::create(
                    &chunk_path,
                    self.temp_compression,
                    self.threads,
                )?;
                for (key, record) in entries.drain(..) {
                    writer.write_record(&key, &record)?;
                }
                writer.finish()?;

                stats.chunks_written += 1;
                chunk_files.push(chunk_path);
                memory_used = 0;
            }
        }

        info!("Read {} records total", stats.total_records);

        if chunk_files.is_empty() {
            // All records fit in memory
            info!("All records fit in memory, performing in-memory sort");

            if self.threads > 1 {
                use rayon::prelude::*;
                entries.par_sort_unstable_by(|a, b| a.0.cmp(&b.0));
            } else {
                entries.sort_unstable_by(|a, b| a.0.cmp(&b.0));
            }

            let output_header = self.create_output_header(header);
            let mut writer = crate::bam_io::create_raw_bam_writer(
                output,
                &output_header,
                self.threads,
                self.output_compression,
            )?;

            for (_key, record) in entries {
                writer.write_raw_record(&record)?;
            }
            writer.finish()?;
        } else {
            // Sort remaining records
            if self.threads > 1 {
                use rayon::prelude::*;
                entries.par_sort_unstable_by(|a, b| a.0.cmp(&b.0));
            } else {
                entries.sort_unstable_by(|a, b| a.0.cmp(&b.0));
            }

            info!(
                "Phase 2: Merging {} chunks (keyed comparisons)...",
                chunk_files.len() + usize::from(!entries.is_empty())
            );

            // Merge disk chunks + in-memory records using keyed comparisons
            self.merge_chunks_generic::<RawQuerynameKey>(&chunk_files, entries, header, output)?;
        }

        stats.output_records = stats.total_records;
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
    fn sort_template_coordinate(
        &self,
        reader: crate::bam_io::BamReaderAuto,
        header: &Header,
        output: &Path,
        temp_path: &Path,
    ) -> Result<RawSortStats> {
        let mut stats = RawSortStats::default();

        // Build library lookup ONCE before sorting for O(1) ordinal lookups
        let lib_lookup = LibraryLookup::from_header(header);

        // Estimate capacity accounting for inline headers and cached-key refs
        // Memory layout per record:
        //   - data: 40 bytes (inline header) + ~250 bytes (BAM record) = 290 bytes
        //   - refs: 48 bytes (TemplateKey + u64 offset + u32 len + u32 pad)
        //   - Total: ~338 bytes per record
        // The larger refs trade memory for cache locality during sort
        let bytes_per_record = 338;
        let estimated_records = self.memory_limit / bytes_per_record;
        // Allocate ~86% for data, ~14% for refs (48/338 ≈ 14%)
        let estimated_data_bytes = self.memory_limit * 86 / 100;

        let mut chunk_files: Vec<PathBuf> = Vec::new();
        let mut buffer =
            TemplateRecordBuffer::with_capacity(estimated_records, estimated_data_bytes);

        let read_ahead = ReadAheadReader::new(reader);

        info!("Phase 1: Reading and sorting chunks (inline buffer)...");

        for record in read_ahead {
            stats.total_records += 1;

            // Extract template key and push to buffer
            let bam_bytes = record.as_ref();
            let key = extract_template_key_inline(bam_bytes, &lib_lookup);
            buffer.push(bam_bytes, key);

            // Check memory usage
            if buffer.memory_usage() >= self.memory_limit {
                // Use .keyed extension to distinguish from BAM chunks
                let chunk_path = temp_path.join(format!("chunk_{:04}.keyed", chunk_files.len()));

                // Sort in place using parallel sort for large buffers
                if self.threads > 1 {
                    buffer.par_sort();
                } else {
                    buffer.sort();
                }

                // Write keyed chunk preserving sort keys for O(1) merge comparisons
                let mut keyed_writer =
                    KeyedChunkWriter::create(&chunk_path, self.temp_compression, self.threads)?;
                for (key, record) in buffer.iter_sorted_keyed() {
                    keyed_writer.write_record(&key, record)?;
                }
                keyed_writer.finish()?;

                stats.chunks_written += 1;
                chunk_files.push(chunk_path);
                buffer.clear();
            }
        }

        info!("Read {} records total", stats.total_records);

        if chunk_files.is_empty() {
            // All records fit in memory
            info!("All records fit in memory, performing in-memory sort");

            if self.threads > 1 {
                buffer.par_sort();
            } else {
                buffer.sort();
            }

            let output_header = self.create_output_header(header);
            let mut writer = crate::bam_io::create_raw_bam_writer(
                output,
                &output_header,
                self.threads,
                self.output_compression,
            )?;

            for record_bytes in buffer.iter_sorted() {
                writer.write_raw_record(record_bytes)?;
            }
            writer.finish()?;
        } else {
            // Sort remaining records and prepare for merge
            if self.threads > 1 {
                buffer.par_sort();
            } else {
                buffer.sort();
            }

            // Collect keyed in-memory records for merge
            let memory_keyed: Vec<(TemplateKey, Vec<u8>)> = if buffer.is_empty() {
                Vec::new()
            } else {
                buffer.iter_sorted_keyed().map(|(k, r)| (k, r.to_vec())).collect()
            };

            info!(
                "Phase 2: Merging {} chunks...",
                chunk_files.len() + usize::from(!memory_keyed.is_empty())
            );

            // Merge using O(1) key comparisons
            self.merge_chunks_keyed(&chunk_files, memory_keyed, header, output)?;
        }

        stats.output_records = stats.total_records;
        info!("Sort complete: {} records processed", stats.total_records);

        Ok(stats)
    }

    /// Merge keyed chunks using O(1) key comparisons.
    ///
    /// This is the fast path for template-coordinate sorting. Instead of parsing
    /// CIGAR/aux data on every comparison (`O(CIGAR_len` + `aux_scan`)), we compare
    /// pre-computed 32-byte keys (O(1) - just 4 u64 comparisons).
    fn merge_chunks_keyed(
        &self,
        chunk_files: &[PathBuf],
        memory_keyed: Vec<(TemplateKey, Vec<u8>)>,
        header: &Header,
        output: &Path,
    ) -> Result<u64> {
        let num_disk = chunk_files.len();
        let has_memory = !memory_keyed.is_empty();
        info!("Merging from {} files and {} in-memory blocks...", num_disk, i32::from(has_memory));

        // Create output header with sort order tags
        let output_header = self.create_output_header(header);

        /// Source for keyed chunks during merge.
        enum KeyedChunkSource {
            /// Disk-based chunk with prefetching reader.
            Disk(KeyedChunkReader),
            /// In-memory sorted records.
            Memory { records: Vec<(TemplateKey, Vec<u8>)>, idx: usize },
        }

        impl KeyedChunkSource {
            fn next_record(&mut self) -> Option<(TemplateKey, Vec<u8>)> {
                match self {
                    KeyedChunkSource::Disk(reader) => reader.next_record(),
                    KeyedChunkSource::Memory { records, idx } => {
                        if *idx < records.len() {
                            // Use replace with zeroed key and empty vec to avoid Default requirement
                            let dummy = (TemplateKey::zeroed(), Vec::new());
                            let entry = std::mem::replace(&mut records[*idx], dummy);
                            *idx += 1;
                            Some(entry)
                        } else {
                            None
                        }
                    }
                }
            }
        }

        // Create unified chunk sources
        let mut sources: Vec<KeyedChunkSource> = Vec::with_capacity(num_disk + 1);

        // Add disk-based chunks with keyed readers
        for path in chunk_files {
            sources.push(KeyedChunkSource::Disk(KeyedChunkReader::open(path)?));
        }

        // Add in-memory keyed chunk
        if has_memory {
            sources.push(KeyedChunkSource::Memory { records: memory_keyed, idx: 0 });
        }

        /// Heap entry for keyed merge - stores pre-computed key for O(1) comparison.
        struct KeyedHeapEntry {
            key: TemplateKey,
            record: Vec<u8>,
            chunk_idx: usize,
        }

        // Initialize heap with first record from each chunk
        let mut heap: Vec<KeyedHeapEntry> = Vec::with_capacity(sources.len());
        for (chunk_idx, source) in sources.iter_mut().enumerate() {
            if let Some((key, record)) = source.next_record() {
                heap.push(KeyedHeapEntry { key, record, chunk_idx });
            }
        }

        let mut heap_size = heap.len();
        if heap_size == 0 {
            info!("Merge complete: 0 records merged");
            return Ok(0);
        }

        // O(1) comparison using pre-computed keys - THIS IS THE KEY OPTIMIZATION
        // No CIGAR parsing, no aux scanning, just 4 u64 comparisons
        let lt = |a: &KeyedHeapEntry, b: &KeyedHeapEntry| -> bool {
            a.key.cmp(&b.key) == Ordering::Greater
        };

        // Build initial heap
        heap_make(&mut heap, &lt);

        // Create raw output writer (bypasses noodles Record encoding for speed)
        let mut writer = crate::bam_io::create_raw_bam_writer(
            output,
            &output_header,
            self.threads,
            self.output_compression,
        )?;

        let mut records_merged = 0u64;

        // Output buffer to reduce write syscalls - stores raw bytes directly
        const OUTPUT_BUFFER_SIZE: usize = 2048;
        let mut output_buffer: Vec<Vec<u8>> = Vec::with_capacity(OUTPUT_BUFFER_SIZE);

        // Merge loop using fixed-array heap with O(1) comparisons
        while heap_size > 0 {
            // Get the minimum record (at heap[0])
            let chunk_idx = heap[0].chunk_idx;
            let record_bytes = std::mem::take(&mut heap[0].record);

            // Buffer the raw bytes directly (no Record conversion)
            output_buffer.push(record_bytes);
            records_merged += 1;

            // Flush buffer when full
            if output_buffer.len() >= OUTPUT_BUFFER_SIZE {
                for rec in output_buffer.drain(..) {
                    writer.write_raw_record(&rec)?;
                }
            }

            // Get next record from the same chunk
            if let Some((key, next_record)) = sources[chunk_idx].next_record() {
                // Replace top with new record and sift down
                heap[0].key = key;
                heap[0].record = next_record;
                heap_sift_down(&mut heap, 0, heap_size, &lt);
            } else {
                // Chunk exhausted - remove from heap
                heap_size -= 1;
                if heap_size > 0 {
                    // Move last element to top and sift down
                    heap.swap(0, heap_size);
                    heap_sift_down(&mut heap, 0, heap_size, &lt);
                }
            }
        }

        // Flush remaining buffered records
        for rec in output_buffer {
            writer.write_raw_record(&rec)?;
        }

        writer.finish()?;

        info!("Merge complete: {} records merged", records_merged);
        Ok(records_merged)
    }

    /// Generic merge for keyed chunks using `O(1)` key comparisons.
    ///
    /// This is the unified merge function that works with any `RawSortKey` type.
    /// It provides `O(1)` comparisons during merge for fixed-size keys (coordinate, template)
    /// and `O(name_len)` for variable-size keys (queryname).
    fn merge_chunks_generic<K: RawSortKey + Default + 'static>(
        &self,
        chunk_files: &[PathBuf],
        memory_keyed: Vec<(K, Vec<u8>)>,
        header: &Header,
        output: &Path,
    ) -> Result<u64> {
        let num_disk = chunk_files.len();
        let has_memory = !memory_keyed.is_empty();
        info!("Merging from {} files and {} in-memory blocks...", num_disk, i32::from(has_memory));

        // Create output header with sort order tags
        let output_header = self.create_output_header(header);

        /// Source for keyed chunks during merge.
        enum GenericKeyedChunkSource<K: RawSortKey + Default + 'static> {
            /// Disk-based chunk with prefetching reader.
            Disk(GenericKeyedChunkReader<K>),
            /// In-memory sorted records.
            Memory { records: Vec<(K, Vec<u8>)>, idx: usize },
        }

        impl<K: RawSortKey + Default + 'static> GenericKeyedChunkSource<K> {
            fn next_record(&mut self) -> Option<(K, Vec<u8>)> {
                match self {
                    GenericKeyedChunkSource::Disk(reader) => reader.next_record(),
                    GenericKeyedChunkSource::Memory { records, idx } => {
                        if *idx < records.len() {
                            // Take ownership of the record
                            let entry = std::mem::take(&mut records[*idx]);
                            *idx += 1;
                            Some(entry)
                        } else {
                            None
                        }
                    }
                }
            }
        }

        // Create unified chunk sources
        let mut sources: Vec<GenericKeyedChunkSource<K>> = Vec::with_capacity(num_disk + 1);

        // Add disk-based chunks with keyed readers
        for path in chunk_files {
            sources.push(GenericKeyedChunkSource::Disk(GenericKeyedChunkReader::<K>::open(path)?));
        }

        // Add in-memory keyed chunk
        if has_memory {
            sources.push(GenericKeyedChunkSource::Memory { records: memory_keyed, idx: 0 });
        }

        /// Heap entry for keyed merge - stores pre-computed key for O(1) comparison.
        struct GenericKeyedHeapEntry<K> {
            key: K,
            record: Vec<u8>,
            chunk_idx: usize,
        }

        // Initialize heap with first record from each chunk
        let mut heap: Vec<GenericKeyedHeapEntry<K>> = Vec::with_capacity(sources.len());
        for (chunk_idx, source) in sources.iter_mut().enumerate() {
            if let Some((key, record)) = source.next_record() {
                heap.push(GenericKeyedHeapEntry { key, record, chunk_idx });
            }
        }

        let mut heap_size = heap.len();
        if heap_size == 0 {
            info!("Merge complete: 0 records merged");
            return Ok(0);
        }

        // O(1) comparison using pre-computed keys
        let lt = |a: &GenericKeyedHeapEntry<K>, b: &GenericKeyedHeapEntry<K>| -> bool {
            a.key.cmp(&b.key) == Ordering::Greater
        };

        // Build initial heap
        heap_make(&mut heap, &lt);

        // Create raw output writer (bypasses noodles Record encoding for speed)
        let mut writer = crate::bam_io::create_raw_bam_writer(
            output,
            &output_header,
            self.threads,
            self.output_compression,
        )?;

        let mut records_merged = 0u64;

        // Output buffer to reduce write syscalls
        const OUTPUT_BUFFER_SIZE: usize = 2048;
        let mut output_buffer: Vec<Vec<u8>> = Vec::with_capacity(OUTPUT_BUFFER_SIZE);

        // Merge loop using fixed-array heap with key comparisons
        while heap_size > 0 {
            // Get the minimum record (at heap[0])
            let chunk_idx = heap[0].chunk_idx;
            let record_bytes = std::mem::take(&mut heap[0].record);

            // Buffer the raw bytes directly (no Record conversion)
            output_buffer.push(record_bytes);
            records_merged += 1;

            // Flush buffer when full
            if output_buffer.len() >= OUTPUT_BUFFER_SIZE {
                for rec in output_buffer.drain(..) {
                    writer.write_raw_record(&rec)?;
                }
            }

            // Get next record from the same chunk
            if let Some((key, next_record)) = sources[chunk_idx].next_record() {
                // Replace top with new record and sift down
                heap[0].key = key;
                heap[0].record = next_record;
                heap_sift_down(&mut heap, 0, heap_size, &lt);
            } else {
                // Chunk exhausted - remove from heap
                heap_size -= 1;
                if heap_size > 0 {
                    // Move last element to top and sift down
                    heap.swap(0, heap_size);
                    heap_sift_down(&mut heap, 0, heap_size, &lt);
                }
            }
        }

        // Flush remaining buffered records
        for rec in output_buffer {
            writer.write_raw_record(&rec)?;
        }

        writer.finish()?;

        info!("Merge complete: {} records merged", records_merged);
        Ok(records_merged)
    }

    /// Merge keyed chunks with BAM index generation.
    ///
    /// Similar to `merge_chunks_generic` but uses `IndexingBamWriter` to build
    /// the BAI index incrementally during the merge. Returns the generated index.
    fn merge_chunks_with_index<K: RawSortKey + Default + 'static>(
        &self,
        chunk_files: &[PathBuf],
        memory_keyed: Vec<(K, Vec<u8>)>,
        header: &Header,
        output: &Path,
    ) -> Result<noodles::bam::bai::Index> {
        use crate::bam_io::create_indexing_bam_writer;

        let num_disk = chunk_files.len();
        let has_memory = !memory_keyed.is_empty();
        info!("Merging from {} files and {} in-memory blocks...", num_disk, i32::from(has_memory));

        let output_header = self.create_output_header(header);

        // Source for keyed chunks during merge
        enum KeyedSource<K: RawSortKey + Default + 'static> {
            Disk(GenericKeyedChunkReader<K>),
            Memory { records: Vec<(K, Vec<u8>)>, idx: usize },
        }

        impl<K: RawSortKey + Default + 'static> KeyedSource<K> {
            fn next_record(&mut self) -> Option<(K, Vec<u8>)> {
                match self {
                    KeyedSource::Disk(reader) => reader.next_record(),
                    KeyedSource::Memory { records, idx } => {
                        if *idx < records.len() {
                            let entry = std::mem::take(&mut records[*idx]);
                            *idx += 1;
                            Some(entry)
                        } else {
                            None
                        }
                    }
                }
            }
        }

        // Create chunk sources
        let mut sources: Vec<KeyedSource<K>> = Vec::with_capacity(num_disk + 1);
        for path in chunk_files {
            sources.push(KeyedSource::Disk(GenericKeyedChunkReader::<K>::open(path)?));
        }
        if has_memory {
            sources.push(KeyedSource::Memory { records: memory_keyed, idx: 0 });
        }

        // Heap entry for merge
        struct HeapEntry<K> {
            key: K,
            record: Vec<u8>,
            chunk_idx: usize,
        }

        // Initialize heap
        let mut heap: Vec<HeapEntry<K>> = Vec::with_capacity(sources.len());
        for (chunk_idx, source) in sources.iter_mut().enumerate() {
            if let Some((key, record)) = source.next_record() {
                heap.push(HeapEntry { key, record, chunk_idx });
            }
        }

        let mut heap_size = heap.len();
        if heap_size == 0 {
            // Empty input - create empty indexed BAM
            let writer = create_indexing_bam_writer(
                output,
                &output_header,
                self.output_compression,
                self.threads,
            )?;
            let index = writer.finish()?;
            info!("Merge complete: 0 records merged");
            return Ok(index);
        }

        let lt =
            |a: &HeapEntry<K>, b: &HeapEntry<K>| -> bool { a.key.cmp(&b.key) == Ordering::Greater };
        heap_make(&mut heap, &lt);

        // Use IndexingBamWriter (single-threaded for accurate virtual positions)
        let mut writer = create_indexing_bam_writer(
            output,
            &output_header,
            self.output_compression,
            self.threads,
        )?;
        let mut records_merged = 0u64;

        // No buffering needed - IndexingBamWriter needs accurate virtual positions
        // for each record, so we write directly
        while heap_size > 0 {
            let chunk_idx = heap[0].chunk_idx;
            let record_bytes = std::mem::take(&mut heap[0].record);

            writer.write_raw_record(&record_bytes)?;
            records_merged += 1;

            if let Some((key, next_record)) = sources[chunk_idx].next_record() {
                heap[0].key = key;
                heap[0].record = next_record;
                heap_sift_down(&mut heap, 0, heap_size, &lt);
            } else {
                heap_size -= 1;
                if heap_size > 0 {
                    heap.swap(0, heap_size);
                    heap_sift_down(&mut heap, 0, heap_size, &lt);
                }
            }
        }

        let index = writer.finish()?;
        info!("Merge complete: {} records merged", records_merged);
        Ok(index)
    }

    /// Create output header with appropriate sort order tags.
    fn create_output_header(&self, header: &Header) -> Header {
        let mut builder = Header::builder();

        // Copy reference sequences
        for (name, seq) in header.reference_sequences() {
            builder = builder.add_reference_sequence(name.as_slice(), seq.clone());
        }

        // Copy read groups
        for (id, rg) in header.read_groups() {
            builder = builder.add_read_group(id.as_slice(), rg.clone());
        }

        // Copy programs
        for (id, pg) in header.programs().as_ref() {
            builder = builder.add_program(id.as_slice(), pg.clone());
        }

        // Copy comments
        for comment in header.comments() {
            builder = builder.add_comment(comment.clone());
        }

        // Set header record with sort order using insert API
        let hd = match self.sort_order {
            SortOrder::Coordinate => {
                Map::<noodles::sam::header::record::value::map::Header>::builder()
                    .insert(header_tag::SORT_ORDER, BString::from("coordinate"))
                    .build()
                    .expect("valid header")
            }
            SortOrder::Queryname => {
                Map::<noodles::sam::header::record::value::map::Header>::builder()
                    .insert(header_tag::SORT_ORDER, BString::from("queryname"))
                    .build()
                    .expect("valid header")
            }
            SortOrder::TemplateCoordinate => {
                Map::<noodles::sam::header::record::value::map::Header>::builder()
                    .insert(header_tag::SORT_ORDER, BString::from("unsorted"))
                    .insert(header_tag::GROUP_ORDER, BString::from("query"))
                    .insert(header_tag::SUBSORT_ORDER, BString::from("template-coordinate"))
                    .build()
                    .expect("valid header")
            }
        };

        builder = builder.set_header(hd);
        builder.build()
    }

    /// Create temporary directory for spill files.
    fn create_temp_dir(&self) -> Result<TempDir> {
        match &self.temp_dir {
            Some(base) => {
                std::fs::create_dir_all(base)?;
                TempDir::new_in(base).context("Failed to create temp directory")
            }
            None => TempDir::new().context("Failed to create temp directory"),
        }
    }
}

/// Extract a packed `TemplateKey` directly from BAM record bytes.
///
/// This function computes the template-coordinate sort key inline, avoiding
/// heap allocations for the read name by using a hash instead.
fn extract_template_key_inline(bam_bytes: &[u8], lib_lookup: &LibraryLookup) -> TemplateKey {
    use crate::sort::bam_fields::{
        find_mc_tag_in_record, find_mi_tag_in_record, flags, get_cigar_ops, mate_unclipped_5prime,
        unclipped_5prime,
    };
    use std::hash::{Hash, Hasher};

    // Extract MI tag using fast raw byte scan (bypasses noodles overhead)
    let mi = find_mi_tag_in_record(bam_bytes);

    // Get library ordinal
    let library = lib_lookup.get_ordinal(bam_bytes);

    // Extract fields from raw bytes
    let tid = i32::from_le_bytes([bam_bytes[0], bam_bytes[1], bam_bytes[2], bam_bytes[3]]);
    let pos = i32::from_le_bytes([bam_bytes[4], bam_bytes[5], bam_bytes[6], bam_bytes[7]]);
    let l_read_name = bam_bytes[8] as usize;
    let flag = u16::from_le_bytes([bam_bytes[14], bam_bytes[15]]);
    let mate_tid = i32::from_le_bytes([bam_bytes[20], bam_bytes[21], bam_bytes[22], bam_bytes[23]]);
    let mate_pos = i32::from_le_bytes([bam_bytes[24], bam_bytes[25], bam_bytes[26], bam_bytes[27]]);

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
    let name_hash = {
        let mut hasher = ahash::AHasher::default();
        name.hash(&mut hasher);
        hasher.finish()
    };

    // Handle unmapped reads
    if is_unmapped {
        if is_paired && !mate_unmapped {
            // Unmapped read with mapped mate - use mate's position as primary key
            let mate_unclipped = find_mc_tag_in_record(bam_bytes)
                .map(|mc| mate_unclipped_5prime(mate_pos, mate_reverse, mc))
                .unwrap_or(mate_pos);

            return TemplateKey::new(
                mate_tid,
                mate_unclipped,
                mate_reverse,
                i32::MAX,
                i32::MAX,
                false,
                library,
                mi,
                name_hash,
                true, // Unmapped read is always "upper" relative to mapped mate
            );
        }

        // Completely unmapped - sort to end
        let is_read2 = (flag & 0x80) != 0; // is_last_segment flag
        return TemplateKey::unmapped(name_hash, is_read2);
    }

    // Calculate unclipped 5' position for this read
    let cigar_ops = get_cigar_ops(bam_bytes);
    let this_pos = unclipped_5prime(pos, is_reverse, cigar_ops);

    // Calculate mate's unclipped 5' position
    let mate_unclipped = if is_paired && !mate_unmapped {
        find_mc_tag_in_record(bam_bytes)
            .map(|mc| mate_unclipped_5prime(mate_pos, mate_reverse, mc))
            .unwrap_or(mate_pos)
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

    TemplateKey::new(tid1, pos1, neg1, tid2, pos2, neg2, library, mi, name_hash, is_upper)
}

/// Statistics from a raw-bytes sort operation.
#[derive(Default, Debug)]
pub struct RawSortStats {
    /// Total records read from input.
    pub total_records: u64,
    /// Records written to output.
    pub output_records: u64,
    /// Number of temporary chunk files written.
    pub chunks_written: usize,
}

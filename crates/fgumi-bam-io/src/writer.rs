//! BAM writer factories and related types.
//!
//! This module provides writer factories for BAM output, including standard BAM writers,
//! raw-bytes BAM writers, and an incremental-indexing BAM writer.

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
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Write};
use std::num::NonZero;
use std::path::Path;

use crate::paths::is_stdout_path;
use crate::vendored::{BlockInfoRx, MultithreadedWriter, MultithreadedWriterBuilder};

/// Maximum uncompressed BGZF block size (64KB - 256 bytes for overhead).
const MAX_BLOCK_SIZE: usize = 65280;

/// Enum wrapping single-threaded and multi-threaded BGZF writers.
///
/// Uses `Box<dyn Write + Send>` to support both file and stdout output.
pub enum BgzfWriterEnum {
    /// Single-threaded BGZF writer
    SingleThreaded(BgzfWriter<Box<dyn Write + Send>>),
    /// Multi-threaded BGZF writer
    MultiThreaded(MultithreadedWriter<Box<dyn Write + Send>>),
}

impl Write for BgzfWriterEnum {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            BgzfWriterEnum::SingleThreaded(w) => w.write(buf),
            BgzfWriterEnum::MultiThreaded(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            BgzfWriterEnum::SingleThreaded(w) => w.flush(),
            BgzfWriterEnum::MultiThreaded(w) => w.flush(),
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
        match self {
            BgzfWriterEnum::SingleThreaded(mut w) => {
                w.flush()?;
                // Single-threaded writer writes EOF on drop
                Ok(())
            }
            BgzfWriterEnum::MultiThreaded(mut w) => {
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

        BgzfWriterEnum::MultiThreaded(builder.build_from_writer(output))
    } else {
        let level = noodles_bgzf::io::writer::CompressionLevel::new(compression_level as u8)
            .unwrap_or_default();
        let writer = noodles_bgzf::io::writer::Builder::default()
            .set_compression_level(level)
            .build_from_writer(output);
        BgzfWriterEnum::SingleThreaded(writer)
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
    indexer: Indexer<LinearIndex>,
    num_refs: usize,
    // Block-aware caching (htslib-style)
    entry_cache: Vec<CachedIndexEntry>,
    block_positions: HashMap<u64, u64>,
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
            indexer: Indexer::default(),
            num_refs,
            entry_cache: Vec::new(),
            block_positions: HashMap::new(),
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

        // Cache entry for later resolution
        let alignment_context = Self::extract_alignment_context(record_bytes);
        self.entry_cache.push(CachedIndexEntry {
            block_number,
            offset_in_block,
            record_len,
            alignment_context,
        });

        // Process any completed blocks
        self.flush_completed_blocks()?;

        Ok(())
    }

    /// Process block completion notifications and resolve cached index entries.
    #[allow(clippy::cast_possible_truncation)]
    fn flush_completed_blocks(&mut self) -> io::Result<()> {
        // Drain block info from receiver
        while let Ok(info) = self.block_info_rx.try_recv() {
            self.block_positions.insert(info.block_number, info.compressed_start);
        }

        // Flush cached entries for completed blocks
        let mut i = 0;
        while i < self.entry_cache.len() {
            let entry = &self.entry_cache[i];

            if let Some(&block_start) = self.block_positions.get(&entry.block_number) {
                // Block is complete - calculate virtual positions
                let start_vpos =
                    VirtualPosition::try_from((block_start, entry.offset_in_block as u16))
                        .unwrap_or(VirtualPosition::MIN);

                // End position: might be in same block or next
                let end_offset = entry.offset_in_block + entry.record_len;
                let end_vpos = if end_offset <= MAX_BLOCK_SIZE {
                    VirtualPosition::try_from((block_start, end_offset as u16))
                        .unwrap_or(VirtualPosition::MIN)
                } else {
                    // Record spans block boundary - need next block position
                    let next_block = entry.block_number + 1;
                    if let Some(&next_start) = self.block_positions.get(&next_block) {
                        let overflow = end_offset - MAX_BLOCK_SIZE;
                        VirtualPosition::try_from((next_start, overflow as u16))
                            .unwrap_or(VirtualPosition::MIN)
                    } else {
                        // Next block not ready yet - skip for now
                        i += 1;
                        continue;
                    }
                };

                let chunk = Chunk::new(start_vpos, end_vpos);
                if let Some(ctx) = entry.alignment_context {
                    self.indexer.add_record(Some(ctx), chunk).map_err(io::Error::other)?;
                } else {
                    // Unmapped record
                    self.indexer.add_record(None, chunk).map_err(io::Error::other)?;
                }

                self.entry_cache.remove(i);
            } else {
                i += 1;
            }
        }

        Ok(())
    }

    /// Extract alignment context from raw BAM bytes.
    ///
    /// Returns `Some((ref_id, start, end, is_mapped))` for mapped reads,
    /// or `None` for unmapped reads (needed for unmapped record counting).
    #[inline]
    #[allow(clippy::cast_sign_loss, clippy::cast_possible_wrap)]
    fn extract_alignment_context(bam: &[u8]) -> Option<(usize, Position, Position, bool)> {
        // Extract fields from BAM record
        let v = fgumi_raw_bam::RawRecordView::new(bam);
        let tid = v.ref_id();
        let pos = v.pos();
        let flags = v.flags();

        let is_unmapped = (flags & fgumi_raw_bam::flags::UNMAPPED) != 0;

        // Unmapped reads: return None (indexer handles them specially)
        if tid < 0 || is_unmapped {
            return None;
        }

        // Calculate alignment end from CIGAR using byte-safe access
        // (doesn't require 4-byte alignment like get_cigar_ops)
        let ref_len = fgumi_raw_bam::reference_length_from_raw_bam(bam);

        // BAM positions are 0-based, Position is 1-based
        let start = Position::try_from((pos + 1) as usize).ok()?;
        let end = Position::try_from((pos + ref_len) as usize).ok()?;

        Some((tid as usize, start, end, true))
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
            if self.entry_cache.is_empty() {
                break;
            }
            std::thread::sleep(std::time::Duration::from_millis(10));
        }

        // Finish the writer (writes EOF marker)
        let _ = self.inner.finish().map_err(|e| io::Error::other(e.to_string()))?;

        // Final flush of any remaining entries
        self.flush_completed_blocks()?;

        // All entries should be flushed now
        if !self.entry_cache.is_empty() {
            return Err(io::Error::other(format!(
                "Unflushed index entries remain: {} entries",
                self.entry_cache.len()
            )));
        }

        // Build the final index
        let index = self.indexer.build(self.num_refs);
        Ok(index)
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
/// * `compression_level` - BGZF compression level (1-12)
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

/// Write a BAI index to a file.
///
/// # Arguments
/// * `path` - Path for the output BAI file (typically `output.bam.bai`)
/// * `index` - The BAI index to write
///
/// # Errors
/// Returns an error if the file cannot be created or writing the index fails.
pub fn write_bai_index<P: AsRef<Path>>(path: P, index: &bai::Index) -> Result<()> {
    let path_ref = path.as_ref();
    let file = File::create(path_ref)
        .with_context(|| format!("Failed to create index file: {}", path_ref.display()))?;
    let mut writer = bai::io::Writer::new(file);
    writer
        .write_index(index)
        .with_context(|| format!("Failed to write index to: {}", path_ref.display()))?;
    Ok(())
}

/// Create a BAM writer and write the header in one operation.
///
/// # Arguments
/// * `path` - Path for the output BAM file (use `-` or `/dev/stdout` for stdout)
/// * `header` - SAM header to write
/// * `threads` - Number of threads for BGZF compression (1 = single-threaded)
/// * `compression_level` - Compression level (1-12)
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

    let bgzf_writer = make_bgzf_writer(output, threads, compression_level);

    let mut writer = noodles::bam::io::Writer::from(bgzf_writer);
    writer
        .write_header(header)
        .with_context(|| format!("Failed to write header to: {}", path_ref.display()))?;
    Ok(writer)
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
/// * `compression_level` - Compression level (1-12)
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
        let mut writer = BgzfWriterEnum::SingleThreaded(BgzfWriter::new(output_file));

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
        let mut writer = BgzfWriterEnum::MultiThreaded(MultithreadedWriter::with_worker_count(
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
        let mut writer = BgzfWriterEnum::SingleThreaded(BgzfWriter::new(output_file));

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
        let mut writer = BgzfWriterEnum::MultiThreaded(MultithreadedWriter::with_worker_count(
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
        let mut writer = BgzfWriterEnum::SingleThreaded(BgzfWriter::new(output_file));

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
        let mut writer = BgzfWriterEnum::MultiThreaded(MultithreadedWriter::with_worker_count(
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
}

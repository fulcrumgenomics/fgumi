//! BAM file I/O utilities.
//!
//! This module provides common utilities for creating BAM readers and writers with consistent
//! error handling and header management.
//!
//! # Threading Model
//!
//! BAM files use BGZF compression, which can be parallelized for both reading and writing:
//!
//! - **Single-threaded**: Use `threads=1` (lower overhead, good for small files)
//! - **Multi-threaded**: Use `threads>1` (higher throughput for large files)
//!
//! Multi-threaded I/O is beneficial for:
//! - Large BAM files where decompression/compression is a bottleneck
//! - Systems with fast storage (SSD/NVMe) where I/O isn't the limiting factor
//! - Pipelines that process many records in parallel

use anyhow::{Context, Result};
use noodles::bam::bai;
// Use old VirtualPosition for compatibility with noodles_csi::Chunk
use noodles::bgzf::VirtualPosition;
use noodles::core::Position;
use noodles::sam::Header;
// Use noodles_bgzf for standard BGZF types
use noodles_bgzf::io::{MultithreadedReader, Reader as BgzfReader, Writer as BgzfWriter};
// Use vendored MultithreadedWriter with position tracking (until upstream PR merges)
use crate::vendored::{BlockInfoRx, MultithreadedWriter, MultithreadedWriterBuilder};
// Use vendored RawBamReader for raw byte access (until noodles exposes Record bytes)
use crate::vendored::RawBamReader;
// Use bgzf crate for CompressionLevel
use bgzf::CompressionLevel;
use noodles_csi::binning_index::Indexer;
use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;
use noodles_csi::binning_index::index::reference_sequence::index::LinearIndex;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, Read, Write};
use std::num::NonZero;
use std::path::Path;

/// Maximum uncompressed BGZF block size (64KB - 256 bytes for overhead).
const MAX_BLOCK_SIZE: usize = 65280;

/// Enum wrapping single-threaded and multi-threaded BGZF readers.
///
/// This allows functions to accept either reader type through a unified interface.
pub enum BgzfReaderEnum {
    /// Single-threaded BGZF reader (lower overhead for small files)
    SingleThreaded(BgzfReader<File>),
    /// Multi-threaded BGZF reader (noodles built-in threading)
    MultiThreaded(MultithreadedReader<File>),
}

impl Read for BgzfReaderEnum {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match self {
            BgzfReaderEnum::SingleThreaded(r) => r.read(buf),
            BgzfReaderEnum::MultiThreaded(r) => r.read(buf),
        }
    }
}

impl BufRead for BgzfReaderEnum {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        match self {
            BgzfReaderEnum::SingleThreaded(r) => r.fill_buf(),
            BgzfReaderEnum::MultiThreaded(r) => r.fill_buf(),
        }
    }

    fn consume(&mut self, amt: usize) {
        match self {
            BgzfReaderEnum::SingleThreaded(r) => r.consume(amt),
            BgzfReaderEnum::MultiThreaded(r) => r.consume(amt),
        }
    }
}

/// Type alias for a BAM reader that supports both single and multi-threaded BGZF.
pub type BamReaderAuto = noodles::bam::io::Reader<BgzfReaderEnum>;

/// Enum wrapping single-threaded and multi-threaded BGZF writers
pub enum BgzfWriterEnum {
    /// Single-threaded BGZF writer
    SingleThreaded(BgzfWriter<File>),
    /// Multi-threaded BGZF writer
    MultiThreaded(MultithreadedWriter<File>),
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

/// Type alias for a BAM writer that supports both single and multi-threaded BGZF
pub type BamWriter = noodles::bam::io::Writer<BgzfWriterEnum>;

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
    #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
    pub fn write_header(&mut self, header: &Header) -> io::Result<()> {
        use std::io::Write;

        // BAM magic
        self.inner.write_all(b"BAM\x01")?;

        // Header text (SAM header serialized using noodles)
        let mut sam_writer = noodles::sam::io::Writer::new(Vec::new());
        sam_writer.write_header(header)?;
        let header_bytes = sam_writer.into_inner();
        let l_text = header_bytes.len() as i32;
        self.inner.write_all(&l_text.to_le_bytes())?;
        self.inner.write_all(&header_bytes)?;

        // Reference sequences
        let n_ref = header.reference_sequences().len() as i32;
        self.inner.write_all(&n_ref.to_le_bytes())?;

        for (name, map) in header.reference_sequences() {
            // l_name: length of name + null terminator
            let l_name = (name.len() + 1) as u32;
            self.inner.write_all(&l_name.to_le_bytes())?;
            self.inner.write_all(name)?;
            self.inner.write_all(&[0u8])?; // null terminator

            // l_ref: reference length
            let l_ref = map.length().get() as i32;
            self.inner.write_all(&l_ref.to_le_bytes())?;
        }

        Ok(())
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
    let output_file = File::create(path_ref)
        .with_context(|| format!("Failed to create output BAM: {}", path_ref.display()))?;

    let bgzf_writer = if threads > 1 {
        let worker_count = NonZero::new(threads).expect("threads > 1 checked above");
        let mut builder = MultithreadedWriterBuilder::default().set_worker_count(worker_count);

        #[allow(clippy::cast_possible_truncation)]
        if let Ok(cl) = CompressionLevel::new(compression_level as u8) {
            builder = builder.set_compression_level(cl);
        }

        BgzfWriterEnum::MultiThreaded(builder.build_from_writer(output_file))
    } else {
        BgzfWriterEnum::SingleThreaded(BgzfWriter::new(output_file))
    };

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
        let block_info_rx = inner.block_info_receiver().unwrap().clone();
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
    #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
    fn write_header(&mut self, header: &Header) -> io::Result<()> {
        // BAM magic
        self.inner.write_all(b"BAM\x01")?;

        // Header text (SAM header serialized using noodles)
        let mut sam_writer = noodles::sam::io::Writer::new(Vec::new());
        sam_writer.write_header(header)?;
        let header_bytes = sam_writer.into_inner();
        let l_text = header_bytes.len() as i32;
        self.inner.write_all(&l_text.to_le_bytes())?;
        self.inner.write_all(&header_bytes)?;

        // Reference sequences
        let n_ref = header.reference_sequences().len() as i32;
        self.inner.write_all(&n_ref.to_le_bytes())?;

        for (name, map) in header.reference_sequences() {
            let l_name = (name.len() + 1) as u32;
            self.inner.write_all(&l_name.to_le_bytes())?;
            self.inner.write_all(name)?;
            self.inner.write_all(&[0u8])?;

            let l_ref = map.length().get() as i32;
            self.inner.write_all(&l_ref.to_le_bytes())?;
        }

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

    /// Write a record buffer to the BAM file and update the index.
    ///
    /// This method encodes the record buffer to raw bytes using noodles' encoder,
    /// then writes it. Useful for the non-fast sorting path that uses decoded records.
    ///
    /// # Errors
    /// Returns an error if encoding or writing the record fails.
    pub fn write_record_buf(
        &mut self,
        header: &Header,
        record: &noodles::sam::alignment::record_buf::RecordBuf,
    ) -> io::Result<()> {
        use crate::vendored::bam_codec::encode_record_buf;

        // Encode record to raw bytes
        let mut buf = Vec::new();
        encode_record_buf(&mut buf, header, record)?;

        // Write using raw record method
        self.write_raw_record(&buf)
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
        let tid = i32::from_le_bytes([bam[0], bam[1], bam[2], bam[3]]);
        let pos = i32::from_le_bytes([bam[4], bam[5], bam[6], bam[7]]);
        let flags = u16::from_le_bytes([bam[14], bam[15]]);

        let is_unmapped = (flags & 0x4) != 0;

        // Unmapped reads: return None (indexer handles them specially)
        if tid < 0 || is_unmapped {
            return None;
        }

        // Calculate alignment end from CIGAR using byte-safe access
        // (doesn't require 4-byte alignment like get_cigar_ops)
        let ref_len = Self::calculate_reference_length(bam);

        // BAM positions are 0-based, Position is 1-based
        let start = Position::try_from((pos + 1) as usize).ok()?;
        let end = Position::try_from((pos + ref_len) as usize).ok()?;

        Some((tid as usize, start, end, true))
    }

    /// Calculate reference-consuming length from CIGAR using byte-safe access.
    ///
    /// This is similar to `reference_length_from_cigar` but reads CIGAR ops
    /// byte-by-byte to avoid alignment requirements.
    #[inline]
    #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
    fn calculate_reference_length(bam: &[u8]) -> i32 {
        let l_read_name = bam[8] as usize;
        let n_cigar_op = u16::from_le_bytes([bam[12], bam[13]]) as usize;

        if n_cigar_op == 0 {
            return 0;
        }

        let cigar_start = 32 + l_read_name;
        let cigar_end = cigar_start + n_cigar_op * 4;

        if cigar_end > bam.len() {
            return 0;
        }

        let mut ref_len = 0i32;

        // Read CIGAR ops byte-by-byte
        for i in 0..n_cigar_op {
            let offset = cigar_start + i * 4;
            let op = u32::from_le_bytes([
                bam[offset],
                bam[offset + 1],
                bam[offset + 2],
                bam[offset + 3],
            ]);

            let op_len = (op >> 4) as i32;
            let op_type = op & 0xF;

            // M (0), D (2), N (3), = (7), X (8) consume reference bases
            if matches!(op_type, 0 | 2 | 3 | 7 | 8) {
                ref_len += op_len;
            }
        }

        ref_len
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

/// Create a BAM reader and read its header.
///
/// # Arguments
/// * `path` - Path to the input BAM file
/// * `threads` - Number of decompression threads (1 = single-threaded)
///
/// # Threading
/// - `threads <= 1`: Single-threaded (lower overhead, good for small files)
/// - `threads > 1`: Multi-threaded (higher throughput for large files)
///
/// # Returns
/// A tuple of (BAM reader, header)
///
/// # Errors
/// Returns an error if the file cannot be opened or the header cannot be read
///
/// # Panics
/// Panics if `threads > 1` but `NonZero::new` fails (should not happen).
///
/// # Example
/// ```no_run
/// use fgumi_lib::bam_io::create_bam_reader;
/// use std::path::Path;
///
/// // Single-threaded
/// let (mut reader, header) = create_bam_reader(Path::new("input.bam"), 1).unwrap();
///
/// // Multi-threaded with 4 decompression threads
/// let (mut reader, header) = create_bam_reader(Path::new("input.bam"), 4).unwrap();
/// ```
pub fn create_bam_reader<P: AsRef<Path>>(
    path: P,
    threads: usize,
) -> Result<(BamReaderAuto, Header)> {
    let path_ref = path.as_ref();
    let file = File::open(path_ref)
        .with_context(|| format!("Failed to open input BAM: {}", path_ref.display()))?;

    let bgzf_reader = if threads > 1 {
        let worker_count = NonZero::new(threads).expect("threads > 1 checked above");
        BgzfReaderEnum::MultiThreaded(MultithreadedReader::with_worker_count(worker_count, file))
    } else {
        BgzfReaderEnum::SingleThreaded(BgzfReader::new(file))
    };

    let mut reader = noodles::bam::io::Reader::from(bgzf_reader);
    let header = reader
        .read_header()
        .with_context(|| format!("Failed to read header from: {}", path_ref.display()))?;

    Ok((reader, header))
}

/// Type alias for a raw BAM reader that supports both single and multi-threaded BGZF.
pub type RawBamReaderAuto = RawBamReader<BgzfReaderEnum>;

/// Create a raw BAM reader that yields raw bytes instead of noodles Record.
///
/// This is used by the raw sorting pipeline for high-performance byte-level access
/// without going through noodles' Record parsing.
///
/// # Arguments
/// * `path` - Path to the input BAM file
/// * `threads` - Number of threads for BGZF decompression (1 = single-threaded)
///
/// # Returns
/// A tuple of (`raw_reader`, `header`)
///
/// # Errors
/// Returns an error if the file cannot be opened or the header cannot be read
///
/// # Panics
/// Panics if `threads > 1` but `NonZero::new` fails (should not happen).
pub fn create_raw_bam_reader<P: AsRef<Path>>(
    path: P,
    threads: usize,
) -> Result<(RawBamReaderAuto, Header)> {
    let path_ref = path.as_ref();
    let file = File::open(path_ref)
        .with_context(|| format!("Failed to open input BAM: {}", path_ref.display()))?;

    let bgzf_reader = if threads > 1 {
        let worker_count = NonZero::new(threads).expect("threads > 1 checked above");
        BgzfReaderEnum::MultiThreaded(MultithreadedReader::with_worker_count(worker_count, file))
    } else {
        BgzfReaderEnum::SingleThreaded(BgzfReader::new(file))
    };

    // Use noodles to read the header, then extract the BGZF reader
    let mut noodles_reader = noodles::bam::io::Reader::from(bgzf_reader);
    let header = noodles_reader
        .read_header()
        .with_context(|| format!("Failed to read header from: {}", path_ref.display()))?;

    // Get back the BGZF reader (header has been consumed)
    let bgzf_reader = noodles_reader.into_inner();

    // Wrap in our raw reader
    let raw_reader = RawBamReader::new(bgzf_reader);

    Ok((raw_reader, header))
}

/// Create a BAM writer and write the header in one operation
///
/// # Arguments
/// * `path` - Path for the output BAM file
/// * `header` - SAM header to write
/// * `threads` - Number of threads for BGZF compression (1 = single-threaded)
/// * `compression_level` - Compression level (1-12)
///
/// # Returns
/// A BAM writer ready for writing records
///
/// # Errors
/// Returns an error if the file cannot be created or the header cannot be written
///
/// # Panics
/// Panics if `threads > 1` but `NonZero::new` fails (should not happen).
///
/// # Example
/// ```no_run
/// use fgumi_lib::bam_io::create_bam_writer;
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
    let output_file = File::create(path_ref)
        .with_context(|| format!("Failed to create output BAM: {}", path_ref.display()))?;

    let bgzf_writer = if threads > 1 {
        // Use multi-threaded BGZF writer with compression level
        let worker_count = NonZero::new(threads).expect("threads > 1 checked above");
        let mut builder = MultithreadedWriterBuilder::default().set_worker_count(worker_count);

        #[allow(clippy::cast_possible_truncation)]
        if let Ok(cl) = CompressionLevel::new(compression_level as u8) {
            builder = builder.set_compression_level(cl);
        }

        BgzfWriterEnum::MultiThreaded(builder.build_from_writer(output_file))
    } else {
        // Use single-threaded BGZF writer
        // Note: Single-threaded writer doesn't support compression level configuration
        // via builder pattern in the same way - would need noodles update
        BgzfWriterEnum::SingleThreaded(BgzfWriter::new(output_file))
    };

    let mut writer = noodles::bam::io::Writer::from(bgzf_writer);
    writer
        .write_header(header)
        .with_context(|| format!("Failed to write header to: {}", path_ref.display()))?;
    Ok(writer)
}

/// Create an optional BAM writer for rejects or filtered reads
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
/// An optional BAM writer, or None if path is None
///
/// # Errors
/// Returns an error if the file cannot be created or the header cannot be written
///
/// # Example
/// ```no_run
/// use fgumi_lib::bam_io::create_optional_bam_writer;
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

/// Check if a path refers to stdin.
///
/// Returns true if the path is "-" or "/dev/stdin".
///
/// # Example
/// ```
/// use fgumi_lib::bam_io::is_stdin_path;
/// use std::path::Path;
///
/// assert!(is_stdin_path(Path::new("-")));
/// assert!(is_stdin_path(Path::new("/dev/stdin")));
/// assert!(!is_stdin_path(Path::new("input.bam")));
/// ```
pub fn is_stdin_path<P: AsRef<Path>>(path: P) -> bool {
    let path_str = path.as_ref().to_string_lossy();
    path_str == "-" || path_str == "/dev/stdin"
}

/// A reader that buffers all bytes read through it.
///
/// This is used to capture raw bytes while parsing the BAM header,
/// so they can be replayed to the pipeline.
struct TeeReader<R> {
    inner: R,
    buffer: Vec<u8>,
}

impl<R: Read> TeeReader<R> {
    fn new(inner: R) -> Self {
        Self { inner, buffer: Vec::new() }
    }

    /// Consume the `TeeReader` and return the buffered bytes and inner reader.
    fn into_parts(self) -> (Vec<u8>, R) {
        (self.buffer, self.inner)
    }
}

impl<R: Read> Read for TeeReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let n = self.inner.read(buf)?;
        self.buffer.extend_from_slice(&buf[..n]);
        Ok(n)
    }
}

/// A reader that chains buffered data with a remaining stream.
///
/// First yields all bytes from the buffer, then reads from the inner reader.
pub struct ChainedReader<R> {
    buffer: io::Cursor<Vec<u8>>,
    inner: R,
    buffer_exhausted: bool,
}

impl<R: Read> ChainedReader<R> {
    /// Create a new chained reader from buffered data and an inner reader.
    pub fn new(buffer: Vec<u8>, inner: R) -> Self {
        Self { buffer: io::Cursor::new(buffer), inner, buffer_exhausted: false }
    }
}

impl<R: Read> Read for ChainedReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if !self.buffer_exhausted {
            let n = self.buffer.read(buf)?;
            if n > 0 {
                return Ok(n);
            }
            self.buffer_exhausted = true;
        }
        self.inner.read(buf)
    }
}

impl<R: Read + Send> ChainedReader<R> {
    /// Box this reader as a trait object for use with the pipeline.
    pub fn into_boxed(self) -> Box<dyn Read + Send>
    where
        R: 'static,
    {
        Box::new(self)
    }
}

/// Create a raw byte reader and header for pipeline use, supporting both files and stdin.
///
/// This function is designed for commands that need to pass a reader to the pipeline.
/// Unlike `create_bam_reader`, this returns a raw byte reader (not a BAM reader) that
/// can be passed directly to `run_bam_pipeline_*_from_reader` functions.
///
/// For files: Opens the file, reads the header, seeks back to start, returns the file.
/// For stdin: Buffers all bytes read while parsing header, returns a chained reader
///            that first yields the buffered bytes then continues from stdin.
///
/// # Arguments
/// * `path` - Path to the input BAM file, or "-" / "/dev/stdin" for stdin
///
/// # Returns
/// A tuple of (boxed reader, header). The reader is positioned at the start of the file
/// (including the header bytes) so the pipeline can parse the header again.
///
/// # Errors
/// Returns an error if the file cannot be opened or the header cannot be read
///
/// # Example
/// ```no_run
/// use fgumi_lib::bam_io::create_bam_reader_for_pipeline;
/// use std::path::Path;
///
/// // From file
/// let (reader, header) = create_bam_reader_for_pipeline(Path::new("input.bam")).unwrap();
///
/// // From stdin
/// let (reader, header) = create_bam_reader_for_pipeline(Path::new("-")).unwrap();
/// ```
pub fn create_bam_reader_for_pipeline<P: AsRef<Path>>(
    path: P,
) -> Result<(Box<dyn Read + Send>, Header)> {
    use std::io::{Seek, SeekFrom};

    let path_ref = path.as_ref();

    if is_stdin_path(path_ref) {
        // Read from stdin with buffering
        let stdin = io::stdin();
        let tee_reader = TeeReader::new(stdin);
        let bgzf_reader = BgzfReader::new(tee_reader);
        let mut bam_reader = noodles::bam::io::Reader::from(bgzf_reader);

        // Read header (this consumes and buffers the raw bytes)
        let header =
            bam_reader.read_header().with_context(|| "Failed to read header from stdin")?;

        // Get the buffered bytes and remaining stdin
        let bgzf_reader = bam_reader.into_inner();
        let tee_reader = bgzf_reader.into_inner();
        let (buffered_bytes, stdin) = tee_reader.into_parts();

        // Create a chained reader: buffered bytes first, then remaining stdin
        let chained = ChainedReader::new(buffered_bytes, stdin);

        Ok((Box::new(chained), header))
    } else {
        // Read from file - we can seek back
        let mut file = File::open(path_ref)
            .with_context(|| format!("Failed to open input BAM: {}", path_ref.display()))?;

        // Read header using noodles
        let bgzf_reader = BgzfReader::new(&file);
        let mut bam_reader = noodles::bam::io::Reader::from(bgzf_reader);
        let header = bam_reader
            .read_header()
            .with_context(|| format!("Failed to read header from: {}", path_ref.display()))?;

        // Seek file back to start
        file.seek(SeekFrom::Start(0))
            .with_context(|| format!("Failed to seek in file: {}", path_ref.display()))?;

        Ok((Box::new(file), header))
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
    fn test_create_bam_reader_nonexistent_file() {
        let result = create_bam_reader("/nonexistent/file.bam", 1);
        assert!(result.is_err());
        if let Err(e) = result {
            let err_msg = e.to_string();
            assert!(err_msg.contains("Failed to open input BAM"));
        }
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
    fn test_roundtrip_write_and_read() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let header = create_test_header();

        // Write (single-threaded)
        {
            let _writer = create_bam_writer(temp_file.path(), &header, 1, 6)?;
            // Writer is dropped, file is written
        }

        // Read (single-threaded)
        let (mut reader, read_header) = create_bam_reader(temp_file.path(), 1)?;

        // Verify header has our reference sequence
        assert_eq!(read_header.reference_sequences().len(), 1);

        // Verify we can iterate (even though there are no records)
        let records: Result<Vec<_>, _> = reader.records().collect();
        assert!(records.is_ok());

        Ok(())
    }

    #[test]
    fn test_roundtrip_write_and_read_multithreaded() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let header = create_test_header();

        // Write (multi-threaded)
        {
            let _writer = create_bam_writer(temp_file.path(), &header, 4, 6)?;
            // Writer is dropped, file is written
        }

        // Read (multi-threaded)
        let (mut reader, read_header) = create_bam_reader(temp_file.path(), 4)?;

        // Verify header has our reference sequence
        assert_eq!(read_header.reference_sequences().len(), 1);

        // Verify we can iterate (even though there are no records)
        let records: Result<Vec<_>, _> = reader.records().collect();
        assert!(records.is_ok());

        Ok(())
    }

    #[test]
    fn test_bgzf_writer_flush_single_threaded() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let output_file = File::create(temp_file.path())?;
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
        let output_file = File::create(temp_file.path())?;
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
        let output_file = File::create(temp_file.path())?;
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
        let output_file = File::create(temp_file.path())?;
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
        let output_file = File::create(temp_file.path())?;
        let mut writer = BgzfWriterEnum::SingleThreaded(BgzfWriter::new(output_file));

        // Test writing via the Write trait
        let bytes_written = writer.write(b"test")?;
        assert_eq!(bytes_written, 4);

        Ok(())
    }

    #[test]
    fn test_bgzf_writer_write_multithreaded() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let output_file = File::create(temp_file.path())?;
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
    fn test_is_stdin_path() {
        // Test stdin paths
        assert!(is_stdin_path("-"));
        assert!(is_stdin_path("/dev/stdin"));
        assert!(is_stdin_path(Path::new("-")));
        assert!(is_stdin_path(Path::new("/dev/stdin")));

        // Test non-stdin paths
        assert!(!is_stdin_path("input.bam"));
        assert!(!is_stdin_path("/path/to/file.bam"));
        assert!(!is_stdin_path(""));
        assert!(!is_stdin_path("/dev/null"));
    }

    #[test]
    fn test_create_bam_reader_for_pipeline_from_file() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let header = create_test_header();

        // Write a BAM file first
        {
            let _writer = create_bam_writer(temp_file.path(), &header, 1, 6)?;
        }

        // Read using create_bam_reader_for_pipeline
        let (mut reader, read_header) = create_bam_reader_for_pipeline(temp_file.path())?;
        assert_eq!(read_header.reference_sequences().len(), 1);

        // The reader returns raw bytes, verify we can read something
        let mut buf = [0u8; 16];
        let n = reader.read(&mut buf)?;
        assert!(n > 0, "Should read some bytes from the file");

        Ok(())
    }

    #[test]
    fn test_create_bam_reader_for_pipeline_nonexistent_file() {
        let result = create_bam_reader_for_pipeline("/nonexistent/file.bam");
        assert!(result.is_err());
        if let Err(e) = result {
            let err_msg = e.to_string();
            assert!(err_msg.contains("Failed to open input BAM"));
        }
    }

    #[test]
    fn test_chained_reader() {
        let buffer = vec![1, 2, 3, 4, 5];
        let remaining = io::Cursor::new(vec![6, 7, 8, 9, 10]);
        let mut chained = ChainedReader::new(buffer, remaining);

        let mut result = Vec::new();
        chained.read_to_end(&mut result).unwrap();

        assert_eq!(result, vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10]);
    }

    #[test]
    fn test_chained_reader_empty_buffer() {
        let buffer = vec![];
        let remaining = io::Cursor::new(vec![1, 2, 3]);
        let mut chained = ChainedReader::new(buffer, remaining);

        let mut result = Vec::new();
        chained.read_to_end(&mut result).unwrap();

        assert_eq!(result, vec![1, 2, 3]);
    }

    #[test]
    fn test_chained_reader_empty_remaining() {
        let buffer = vec![1, 2, 3];
        let remaining = io::Cursor::new(vec![]);
        let mut chained = ChainedReader::new(buffer, remaining);

        let mut result = Vec::new();
        chained.read_to_end(&mut result).unwrap();

        assert_eq!(result, vec![1, 2, 3]);
    }

    /// Create a minimal BAM record for testing.
    ///
    /// Creates a mapped record at the given reference and position with a simple CIGAR.
    /// The read name is padded to ensure CIGAR alignment (the read name length must make
    /// the CIGAR start at a 4-byte aligned offset from the record start).
    #[allow(clippy::cast_possible_truncation)]
    fn create_test_bam_record(ref_id: i32, pos: i32, read_name: &[u8]) -> Vec<u8> {
        // BAM record structure (without block_size prefix):
        // refID (4), pos (4), l_read_name (1), mapq (1), bin (2), n_cigar_op (2),
        // flag (2), l_seq (4), next_refID (4), next_pos (4), tlen (4),
        // read_name (l_read_name), cigar (n_cigar_op * 4), seq ((l_seq+1)/2), qual (l_seq)
        //
        // Fixed header is 32 bytes. CIGAR starts at 32 + l_read_name.
        // For CIGAR to be 4-byte aligned, l_read_name must be a multiple of 4.

        // Calculate padding needed to make l_read_name a multiple of 4
        let name_with_null = read_name.len() + 1;
        let padding = (4 - (name_with_null % 4)) % 4;
        let l_read_name = (name_with_null + padding) as u8;

        let mapq: u8 = 60;
        let bin: u16 = 4681; // bin for small coordinates
        let n_cigar_op: u16 = 1;
        let flag: u16 = 0; // mapped, forward strand
        let l_seq: u32 = 10;
        let next_ref_id: i32 = -1;
        let next_pos: i32 = -1;
        let tlen: i32 = 0;

        let mut record = Vec::new();

        // Core fields
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

        // Read name + null terminator + padding (null bytes)
        record.extend_from_slice(read_name);
        record.push(0); // null terminator
        record.extend(std::iter::repeat_n(0u8, padding));

        // CIGAR: 10M = (10 << 4) = 160
        let cigar_op: u32 = 10 << 4; // 10M
        record.extend_from_slice(&cigar_op.to_le_bytes());

        // Sequence: 10 bases = 5 bytes (4-bit encoded)
        // AAAAAAAAAA = 0x11 0x11 0x11 0x11 0x11
        record.extend_from_slice(&[0x11, 0x11, 0x11, 0x11, 0x11]);

        // Quality: 10 bytes of qual 30
        record.extend_from_slice(&[30u8; 10]);

        record
    }

    #[test]
    fn test_create_indexing_bam_writer() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let header = create_test_header();

        let writer = create_indexing_bam_writer(temp_file.path(), &header, 6, 2);
        assert!(writer.is_ok());

        // Just finish without writing records
        let writer = writer.unwrap();
        let index = writer.finish();
        assert!(index.is_ok());

        Ok(())
    }

    #[test]
    fn test_indexing_bam_writer_with_records() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let header = create_test_header();

        let mut writer = create_indexing_bam_writer(temp_file.path(), &header, 6, 2)?;

        // Write some test records
        for i in 0..100 {
            let record = create_test_bam_record(0, i * 100, format!("read{i}").as_bytes());
            writer.write_raw_record(&record)?;
        }

        let index = writer.finish()?;

        // Verify index has reference sequences
        assert!(!index.reference_sequences().is_empty());

        Ok(())
    }

    #[test]
    fn test_indexing_bam_writer_produces_valid_bam() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let header = create_test_header();

        // Write BAM with index
        {
            let mut writer = create_indexing_bam_writer(temp_file.path(), &header, 6, 2)?;

            for i in 0..10 {
                let record = create_test_bam_record(0, i * 100, format!("read{i}").as_bytes());
                writer.write_raw_record(&record)?;
            }

            let _index = writer.finish()?;
        }

        // Read back and verify
        let (mut reader, read_header) = create_bam_reader(temp_file.path(), 1)?;
        assert_eq!(read_header.reference_sequences().len(), 1);

        let records: Vec<_> = reader.records().collect();
        assert_eq!(records.len(), 10);

        Ok(())
    }

    #[test]
    fn test_indexing_bam_writer_multithreaded() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let header = create_test_header();

        let mut writer = create_indexing_bam_writer(temp_file.path(), &header, 6, 4)?;

        // Write enough records to potentially trigger multiple blocks
        for i in 0..1000 {
            let record = create_test_bam_record(0, i * 10, format!("read{i:04}").as_bytes());
            writer.write_raw_record(&record)?;
        }

        let index = writer.finish()?;

        // Verify index has reference sequences with bins
        assert!(!index.reference_sequences().is_empty());

        Ok(())
    }

    #[test]
    fn test_indexing_bam_writer_unmapped_records() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let header = create_test_header();

        let mut writer = create_indexing_bam_writer(temp_file.path(), &header, 6, 2)?;

        // Write unmapped record (ref_id = -1)
        let record = create_test_bam_record(-1, -1, b"unmapped_read");
        writer.write_raw_record(&record)?;

        // Write mapped record
        let record = create_test_bam_record(0, 100, b"mapped_read");
        writer.write_raw_record(&record)?;

        let index = writer.finish()?;
        assert!(!index.reference_sequences().is_empty());

        Ok(())
    }
}

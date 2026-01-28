//! External merge-sort implementation for BAM files.
//!
//! Implements a high-performance external sorting algorithm that handles BAM files
//! larger than available RAM by spilling sorted chunks to temporary files.
//!
//! # Algorithm
//!
//! 1. **Accumulate phase**: Read records into memory until limit reached
//! 2. **Sort phase**: Parallel sort using rayon
//! 3. **Spill phase**: Write sorted chunk to temp file with fast compression
//! 4. **Merge phase**: K-way merge using binary heap
//!
//! # Performance Features
//!
//! - **Lazy key extraction**: Only parse fields needed for sorting
//! - **Parallel sorting**: Uses rayon for in-memory parallel sort
//! - **Buffer recycling**: Reuses temp file buffers
//! - **Fast compression**: Level 1 for temp files, configurable for output

use crate::bam_io::{
    BamReaderAuto, create_bam_reader, create_bam_writer, create_indexing_bam_writer,
    write_bai_index,
};
use crate::sort::keys::{CoordinateKey, QuerynameKey, SortKey, SortOrder, TemplateCoordinateKey};
use anyhow::{Context, Result};
use bstr::BString;
use log::info;
use noodles::bam;
use noodles::sam::Header;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::header::tag as header_tag;
use rayon::prelude::*;
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};
use tempfile::TempDir;

/// Default memory limit per chunk (512 MB).
const DEFAULT_MEMORY_LIMIT: usize = 512 * 1024 * 1024;

/// Estimated bytes per BAM record for memory accounting.
/// This should be conservative to prevent OOM but not so high that we spill unnecessarily.
/// Actual `RecordBuf` + Key overhead is ~800-1200 bytes for typical BAM records.
#[allow(dead_code)]
const ESTIMATED_RECORD_SIZE: usize = 800;

/// Buffer size for reading temp files during merge.
const MERGE_BUFFER_SIZE: usize = 64 * 1024;

/// Compression level for temporary files (fast compression).
const TEMP_COMPRESSION_LEVEL: u32 = 1;

/// External sorter for BAM files.
pub struct ExternalSorter {
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
    /// Whether to generate a BAI index (coordinate sort only).
    write_index: bool,
    /// Program record info (version, `command_line`) for @PG header.
    pg_info: Option<(String, String)>,
}

impl ExternalSorter {
    /// Create a new external sorter with the given sort order.
    #[must_use]
    pub fn new(sort_order: SortOrder) -> Self {
        Self {
            sort_order,
            memory_limit: DEFAULT_MEMORY_LIMIT,
            temp_dir: None,
            threads: 1,
            output_compression: 6,
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

    /// Enable BAI index generation (coordinate sort only).
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

    /// Sort a BAM file.
    pub fn sort(&self, input: &Path, output: &Path) -> Result<SortStats> {
        info!("Starting sort with order: {:?}", self.sort_order);
        info!("Memory limit: {} MB", self.memory_limit / (1024 * 1024));
        info!("Threads: {}", self.threads);
        if self.write_index {
            info!("Write index: enabled");
        }

        // Open input BAM
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

        // Dispatch to sort order specific implementation
        match self.sort_order {
            SortOrder::Coordinate if self.write_index => {
                self.sort_coordinate_with_index(reader, &header, output, temp_path)
            }
            SortOrder::Coordinate => {
                self.sort_with_key::<CoordinateKey>(reader, &header, output, temp_path)
            }
            SortOrder::Queryname => {
                self.sort_with_key::<QuerynameKey>(reader, &header, output, temp_path)
            }
            SortOrder::TemplateCoordinate => {
                self.sort_with_key::<TemplateCoordinateKey>(reader, &header, output, temp_path)
            }
        }
    }

    /// Generic sort implementation for any key type.
    fn sort_with_key<K: SortKey + 'static>(
        &self,
        mut reader: BamReaderAuto,
        header: &Header,
        output: &Path,
        temp_path: &Path,
    ) -> Result<SortStats> {
        let max_records = self.memory_limit / ESTIMATED_RECORD_SIZE;
        let mut stats = SortStats::default();

        // Build context once (e.g., LibraryIndex for TemplateCoordinate)
        let ctx = K::build_context(header);

        // Phase 1: Read and sort chunks
        let mut chunk_files: Vec<PathBuf> = Vec::new();
        let mut records: Vec<(K, RecordBuf)> = Vec::with_capacity(max_records);
        let mut memory_used = 0usize;

        info!("Phase 1: Reading and sorting chunks...");

        for result in reader.record_bufs(header) {
            let record = result.context("Failed to read BAM record")?;
            stats.total_records += 1;

            // Extract sort key using pre-built context
            let key = K::from_record(&record, header, &ctx)?;

            // Estimate record memory usage
            let record_size = estimate_record_size(&record);
            memory_used += record_size;

            records.push((key, record));

            // Check if we need to spill to disk
            if memory_used >= self.memory_limit {
                let chunk_path = temp_path.join(format!("chunk_{:04}.bam", chunk_files.len()));
                self.sort_and_write_chunk(&mut records, header, &chunk_path)?;
                stats.chunks_written += 1;
                chunk_files.push(chunk_path);
                records.clear();
                memory_used = 0;
            }
        }

        info!("Read {} records total", stats.total_records);

        // Phase 2: Handle remaining records
        if chunk_files.is_empty() {
            // All records fit in memory - simple in-memory sort
            info!("All records fit in memory, performing in-memory sort");
            self.sort_and_write_final(&mut records, header, output)?;
        } else {
            // Need to merge chunks
            if !records.is_empty() {
                // Write remaining records as final chunk
                let chunk_path = temp_path.join(format!("chunk_{:04}.bam", chunk_files.len()));
                self.sort_and_write_chunk(&mut records, header, &chunk_path)?;
                stats.chunks_written += 1;
                chunk_files.push(chunk_path);
            }

            info!("Phase 2: Merging {} chunks...", chunk_files.len());

            // Phase 3: K-way merge
            self.merge_chunks::<K>(&chunk_files, header, output, &ctx)?;
        }

        stats.output_records = stats.total_records;
        info!("Sort complete: {} records processed", stats.total_records);

        Ok(stats)
    }

    /// Sort records in parallel and write to a chunk file.
    fn sort_and_write_chunk<K: SortKey + Send>(
        &self,
        records: &mut [(K, RecordBuf)],
        header: &Header,
        path: &Path,
    ) -> Result<()> {
        // Parallel sort using rayon
        if self.threads > 1 {
            records.par_sort_unstable_by(|(k1, _), (k2, _)| k1.cmp(k2));
        } else {
            records.sort_unstable_by(|(k1, _), (k2, _)| k1.cmp(k2));
        }

        // Write to temp file with fast compression
        let mut writer = create_bam_writer(path, header, 1, TEMP_COMPRESSION_LEVEL)?;

        for (_, record) in records.iter() {
            writer.write_alignment_record(header, record)?;
        }

        Ok(())
    }

    /// Sort records and write to final output.
    fn sort_and_write_final<K: SortKey + Send>(
        &self,
        records: &mut [(K, RecordBuf)],
        header: &Header,
        output: &Path,
    ) -> Result<()> {
        // Parallel sort using rayon
        if self.threads > 1 {
            records.par_sort_unstable_by(|(k1, _), (k2, _)| k1.cmp(k2));
        } else {
            records.sort_unstable_by(|(k1, _), (k2, _)| k1.cmp(k2));
        }

        // Modify header with sort order tags
        let output_header = self.create_output_header(header);

        // Write to output with configured compression
        let mut writer =
            create_bam_writer(output, &output_header, self.threads, self.output_compression)?;

        for (_, record) in records.iter() {
            writer.write_alignment_record(&output_header, record)?;
        }

        Ok(())
    }

    /// K-way merge of sorted chunk files.
    fn merge_chunks<K: SortKey + 'static>(
        &self,
        chunk_files: &[PathBuf],
        header: &Header,
        output: &Path,
        ctx: &K::Context,
    ) -> Result<()> {
        // Create output header with sort order tags
        let output_header = self.create_output_header(header);

        // Open chunk readers
        let mut chunk_readers: Vec<ChunkReader<K>> = chunk_files
            .iter()
            .enumerate()
            .map(|(idx, path)| ChunkReader::new(path, header, idx))
            .collect::<Result<Vec<_>>>()?;

        // Initialize heap with first record from each chunk
        let mut heap: BinaryHeap<Reverse<HeapEntry<K>>> =
            BinaryHeap::with_capacity(chunk_files.len());

        for reader in &mut chunk_readers {
            if let Some((key, record)) = reader.next(header, ctx)? {
                heap.push(Reverse(HeapEntry { key, record, chunk_idx: reader.idx }));
            }
        }

        // Create output writer
        let mut writer =
            create_bam_writer(output, &output_header, self.threads, self.output_compression)?;

        // Merge using heap
        while let Some(Reverse(entry)) = heap.pop() {
            // Write record to output
            writer.write_alignment_record(&output_header, &entry.record)?;

            // Get next record from the same chunk
            let reader = &mut chunk_readers[entry.chunk_idx];
            if let Some((key, record)) = reader.next(header, ctx)? {
                heap.push(Reverse(HeapEntry { key, record, chunk_idx: reader.idx }));
            }
        }

        Ok(())
    }

    /// Coordinate sort with BAI index generation.
    ///
    /// Similar to `sort_with_key` but uses `IndexingBamWriter` to build
    /// the BAI index incrementally during write.
    fn sort_coordinate_with_index(
        &self,
        mut reader: BamReaderAuto,
        header: &Header,
        output: &Path,
        temp_path: &Path,
    ) -> Result<SortStats> {
        let max_records = self.memory_limit / ESTIMATED_RECORD_SIZE;
        let mut stats = SortStats::default();

        // Build context (empty for CoordinateKey, but needed for generic consistency)
        #[allow(clippy::let_unit_value)]
        let ctx = CoordinateKey::build_context(header);

        info!("Indexing enabled: will write BAM index alongside output");

        // Phase 1: Read and sort chunks
        let mut chunk_files: Vec<PathBuf> = Vec::new();
        let mut records: Vec<(CoordinateKey, RecordBuf)> = Vec::with_capacity(max_records);
        let mut memory_used = 0usize;

        info!("Phase 1: Reading and sorting chunks...");

        for result in reader.record_bufs(header) {
            let record = result.context("Failed to read BAM record")?;
            stats.total_records += 1;

            let key = CoordinateKey::from_record(&record, header, &ctx)?;
            let record_size = estimate_record_size(&record);
            memory_used += record_size;

            records.push((key, record));

            if memory_used >= self.memory_limit {
                let chunk_path = temp_path.join(format!("chunk_{:04}.bam", chunk_files.len()));
                self.sort_and_write_chunk(&mut records, header, &chunk_path)?;
                stats.chunks_written += 1;
                chunk_files.push(chunk_path);
                records.clear();
                memory_used = 0;
            }
        }

        info!("Read {} records total", stats.total_records);

        let output_header = self.create_output_header(header);

        // Phase 2: Handle remaining records
        if chunk_files.is_empty() {
            // All records fit in memory
            info!("All records fit in memory, performing in-memory sort");
            self.sort_and_write_final_with_index(&mut records, &output_header, output)?;
        } else {
            // Need to merge chunks
            if !records.is_empty() {
                let chunk_path = temp_path.join(format!("chunk_{:04}.bam", chunk_files.len()));
                self.sort_and_write_chunk(&mut records, header, &chunk_path)?;
                stats.chunks_written += 1;
                chunk_files.push(chunk_path);
            }

            info!("Phase 2: Merging {} chunks with index generation...", chunk_files.len());
            self.merge_chunks_with_index(&chunk_files, header, &output_header, output, &ctx)?;
        }

        stats.output_records = stats.total_records;
        info!("Sort complete: {} records processed", stats.total_records);

        Ok(stats)
    }

    /// Sort records and write to final output with index generation.
    fn sort_and_write_final_with_index(
        &self,
        records: &mut [(CoordinateKey, RecordBuf)],
        output_header: &Header,
        output: &Path,
    ) -> Result<()> {
        // Parallel sort using rayon
        if self.threads > 1 {
            records.par_sort_unstable_by(|(k1, _), (k2, _)| k1.cmp(k2));
        } else {
            records.sort_unstable_by(|(k1, _), (k2, _)| k1.cmp(k2));
        }

        // Use IndexingBamWriter
        let mut writer = create_indexing_bam_writer(
            output,
            output_header,
            self.output_compression,
            self.threads,
        )?;

        for (_, record) in records.iter() {
            writer.write_record_buf(output_header, record)?;
        }

        let index = writer.finish()?;

        // Write the index file
        let index_path = output.with_extension("bam.bai");
        write_bai_index(&index_path, &index)?;
        info!("Wrote BAM index: {}", index_path.display());

        Ok(())
    }

    /// K-way merge with index generation.
    #[allow(clippy::trivially_copy_pass_by_ref)] // Context is generic, () for CoordinateKey
    fn merge_chunks_with_index(
        &self,
        chunk_files: &[PathBuf],
        header: &Header,
        output_header: &Header,
        output: &Path,
        ctx: &<CoordinateKey as SortKey>::Context,
    ) -> Result<()> {
        // Open chunk readers
        let mut chunk_readers: Vec<ChunkReader<CoordinateKey>> = chunk_files
            .iter()
            .enumerate()
            .map(|(idx, path)| ChunkReader::new(path, header, idx))
            .collect::<Result<Vec<_>>>()?;

        // Initialize heap
        let mut heap: BinaryHeap<Reverse<HeapEntry<CoordinateKey>>> =
            BinaryHeap::with_capacity(chunk_files.len());

        for reader in &mut chunk_readers {
            if let Some((key, record)) = reader.next(header, ctx)? {
                heap.push(Reverse(HeapEntry { key, record, chunk_idx: reader.idx }));
            }
        }

        // Create indexing writer
        let mut writer = create_indexing_bam_writer(
            output,
            output_header,
            self.output_compression,
            self.threads,
        )?;

        let mut records_merged = 0u64;

        // Merge using heap
        while let Some(Reverse(entry)) = heap.pop() {
            writer.write_record_buf(output_header, &entry.record)?;
            records_merged += 1;

            let reader = &mut chunk_readers[entry.chunk_idx];
            if let Some((key, record)) = reader.next(header, ctx)? {
                heap.push(Reverse(HeapEntry { key, record, chunk_idx: reader.idx }));
            }
        }

        let index = writer.finish()?;

        // Write the index file
        let index_path = output.with_extension("bam.bai");
        write_bai_index(&index_path, &index)?;
        info!("Wrote BAM index: {}", index_path.display());
        info!("Merge complete: {} records merged", records_merged);

        Ok(())
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
                // Template-coordinate uses: SO:unsorted, GO:query, SS:template-coordinate
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

/// Estimate memory usage of a BAM record plus its sort key.
fn estimate_record_size(record: &RecordBuf) -> usize {
    // Base RecordBuf struct size + variable-length fields
    let record_size = std::mem::size_of::<RecordBuf>()
        + record.name().map_or(0, |n| n.len())
        + record.sequence().len()
        + record.quality_scores().as_ref().len()
        + record.cigar().as_ref().len() * 4
        + 256; // Estimated tag overhead

    // Add estimated sort key size (varies by sort order, ~100-200 bytes typical)
    // TemplateCoordinateKey is the largest at ~150 bytes
    record_size + 200
}

/// Reader for a sorted chunk file.
struct ChunkReader<K> {
    reader: bam::io::Reader<noodles::bgzf::io::Reader<BufReader<File>>>,
    idx: usize,
    _phantom: std::marker::PhantomData<K>,
}

impl<K: SortKey> ChunkReader<K> {
    fn new(path: &Path, _header: &Header, idx: usize) -> Result<Self> {
        let file = File::open(path).context("Failed to open chunk file")?;
        let buf_reader = BufReader::with_capacity(MERGE_BUFFER_SIZE, file);
        let mut reader = bam::io::Reader::new(buf_reader);

        // Read and discard header
        reader.read_header()?;

        Ok(Self { reader, idx, _phantom: std::marker::PhantomData })
    }

    fn next(&mut self, header: &Header, ctx: &K::Context) -> Result<Option<(K, RecordBuf)>> {
        let mut record = RecordBuf::default();
        match self.reader.read_record_buf(header, &mut record) {
            Ok(0) => Ok(None), // EOF
            Ok(_) => {
                let key = K::from_record(&record, header, ctx)?;
                Ok(Some((key, record)))
            }
            Err(e) => Err(e.into()),
        }
    }
}

/// Entry in the merge heap.
struct HeapEntry<K> {
    key: K,
    record: RecordBuf,
    chunk_idx: usize,
}

impl<K: Ord> PartialEq for HeapEntry<K> {
    fn eq(&self, other: &Self) -> bool {
        self.key == other.key
    }
}

impl<K: Ord> Eq for HeapEntry<K> {}

impl<K: Ord> PartialOrd for HeapEntry<K> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<K: Ord> Ord for HeapEntry<K> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.key.cmp(&other.key)
    }
}

/// Statistics from a sort operation.
#[derive(Default, Debug)]
pub struct SortStats {
    /// Total records read from input.
    pub total_records: u64,
    /// Records written to output.
    pub output_records: u64,
    /// Number of temporary chunk files written.
    pub chunks_written: usize,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sam::builder::RecordBuilder;

    #[test]
    fn test_estimate_record_size() {
        let record = RecordBuilder::new().name("test").build();
        let size = estimate_record_size(&record);
        assert!(size > 0);
    }

    #[test]
    fn test_sorter_creation() {
        let sorter = ExternalSorter::new(SortOrder::Coordinate)
            .memory_limit(1024 * 1024)
            .threads(4)
            .output_compression(6);

        assert_eq!(sorter.memory_limit, 1024 * 1024);
        assert_eq!(sorter.threads, 4);
        assert_eq!(sorter.output_compression, 6);
    }
}

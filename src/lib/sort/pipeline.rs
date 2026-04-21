//! Pipelined merge with parallel I/O for BAM sorting.
//!
//! This module provides a multi-threaded merge implementation that overlaps
//! I/O with computation to maximize throughput during the merge phase of
//! external sort.
//!
//! # Architecture
//!
//! ```text
//! ┌─────────────┐    ┌─────────────┐    ┌─────────────┐
//! │ Reader Pool │───>│ Merge Heap  │───>│   Writer    │
//! │ (N threads) │    │ (1 thread)  │    │ (M threads) │
//! └─────────────┘    └─────────────┘    └─────────────┘
//!      │                   │                   │
//!      ▼                   ▼                   ▼
//!   Decompress          K-way              Compress
//!   in parallel         merge              in parallel
//! ```
//!
//! # Performance Benefits
//!
//! - **Parallel decompression**: Each chunk file is read by its own thread
//! - **Overlapped I/O**: Reading and decompression happen in background
//! - **Buffered prefetch**: Records are prefetched ahead of merge consumption
//! - **Multi-threaded output**: Writing uses parallel BGZF compression

use anyhow::{Context, Result, anyhow};
use crossbeam_channel::{Receiver, Sender, bounded};
use fgumi_raw_bam::{RawRecord, read_raw_record};
use noodles::bgzf;
use noodles::sam::Header;
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};
use std::thread::{self, JoinHandle};

use super::MERGE_BUFFER_SIZE;
use crate::bam_io::create_raw_bam_writer;

/// Number of records to prefetch per chunk reader.
const PREFETCH_BUFFER_SIZE: usize = 128;

/// A record with its sort key and source chunk index.
pub struct MergeEntry<K> {
    pub key: K,
    pub record: RawRecord,
    pub chunk_idx: usize,
}

impl<K: PartialEq> PartialEq for MergeEntry<K> {
    fn eq(&self, other: &Self) -> bool {
        self.key == other.key
    }
}

impl<K: Eq> Eq for MergeEntry<K> {}

impl<K: PartialOrd> PartialOrd for MergeEntry<K> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.key.partial_cmp(&other.key)
    }
}

impl<K: Ord> Ord for MergeEntry<K> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.key.cmp(&other.key)
    }
}

/// Configuration for parallel merge.
pub struct ParallelMergeConfig {
    /// Number of reader threads (typically `min(num_chunks, available_threads)`).
    pub reader_threads: usize,
    /// Number of writer threads for output compression.
    pub writer_threads: usize,
    /// Compression level for output.
    pub compression_level: u32,
}

impl Default for ParallelMergeConfig {
    fn default() -> Self {
        Self { reader_threads: 4, writer_threads: 4, compression_level: 6 }
    }
}

/// One message from a [`PrefetchingChunkReader`] — `Ok(Some(_))` carries a record,
/// `Ok(None)` signals clean EOF, and `Err(_)` carries a reader/I/O failure that
/// must abort the merge instead of being silently treated as EOF.
type ChunkReaderMsg = std::result::Result<Option<RawRecord>, String>;

/// Prefetching chunk reader that runs in a background thread.
///
/// This reader maintains a buffer of prefetched records, allowing the merge
/// thread to consume records without blocking on I/O.
struct PrefetchingChunkReader {
    /// Receiver for prefetched records.
    record_rx: Receiver<ChunkReaderMsg>,
    /// Handle to the reader thread.
    _handle: JoinHandle<()>,
    /// Chunk index for heap management.
    idx: usize,
}

impl PrefetchingChunkReader {
    /// Create a new prefetching reader for a chunk file.
    #[allow(clippy::unnecessary_wraps)]
    fn new(path: PathBuf, idx: usize) -> Result<Self> {
        // Channel for prefetched records
        let (record_tx, record_rx) = bounded(PREFETCH_BUFFER_SIZE);

        // Spawn reader thread. Reader-side failures are surfaced over the
        // channel so the merge loop can abort instead of truncating.
        let handle = thread::spawn(move || {
            Self::reader_thread(path, record_tx);
        });

        Ok(Self { record_rx, _handle: handle, idx })
    }

    /// Reader thread function.
    #[allow(clippy::needless_pass_by_value)]
    fn reader_thread(path: PathBuf, tx: Sender<ChunkReaderMsg>) {
        let file = match File::open(&path).context("Failed to open chunk file") {
            Ok(f) => f,
            Err(e) => {
                let _ = tx.send(Err(format!("{e:#}")));
                return;
            }
        };
        let buf_reader = BufReader::with_capacity(MERGE_BUFFER_SIZE, file);
        let bgzf_reader = bgzf::io::Reader::new(buf_reader);
        let mut noodles_reader = noodles::bam::io::Reader::from(bgzf_reader);

        // Read and discard header
        if let Err(e) = noodles_reader.read_header() {
            let _ = tx.send(Err(format!("Failed to read chunk header: {e}")));
            return;
        }

        // Extract the inner BGZF reader (header already consumed)
        let mut bgzf_reader = noodles_reader.into_inner();

        // Read raw records and send to channel
        let mut record = RawRecord::new();
        loop {
            match read_raw_record(&mut bgzf_reader, &mut record) {
                Ok(0) => {
                    // Clean EOF
                    let _ = tx.send(Ok(None));
                    break;
                }
                Ok(_) => {
                    // Take ownership of current record and replace with empty one
                    let owned_record = std::mem::take(&mut record);
                    if tx.send(Ok(Some(owned_record))).is_err() {
                        // Receiver dropped, exit
                        break;
                    }
                }
                Err(e) => {
                    let _ = tx.send(Err(format!("Error reading chunk: {e}")));
                    break;
                }
            }
        }
    }

    /// Get the next record from the prefetch buffer.
    ///
    /// Returns `Ok(Some(_))` for a record, `Ok(None)` for clean EOF, and
    /// `Err(_)` if the reader thread reported a failure or the channel was
    /// dropped unexpectedly.
    fn next(&self) -> Result<Option<RawRecord>> {
        match self.record_rx.recv() {
            Ok(Ok(record)) => Ok(record),
            Ok(Err(e)) => Err(anyhow!("{e}")),
            Err(e) => Err(anyhow!("chunk reader channel closed: {e}")),
        }
    }
}

/// Parallel merge implementation using prefetching readers.
///
/// # Errors
///
/// Returns an error if reading chunks, writing output, or merging fails.
#[allow(clippy::needless_pass_by_value)]
pub fn parallel_merge<K, F>(
    chunk_files: &[PathBuf],
    _header: &Header,
    output_header: &Header,
    output: &Path,
    extract_key: F,
    config: ParallelMergeConfig,
) -> Result<u64>
where
    K: Clone + Send + Sync + Ord,
    F: Fn(&RawRecord) -> K + Send + Sync,
{
    tracing::info!(
        "Starting parallel merge of {} chunks with {} reader threads",
        chunk_files.len(),
        config.reader_threads.min(chunk_files.len())
    );

    // Create prefetching readers for each chunk
    let chunk_readers: Vec<PrefetchingChunkReader> = chunk_files
        .iter()
        .enumerate()
        .map(|(idx, path)| PrefetchingChunkReader::new(path.clone(), idx))
        .collect::<Result<Vec<_>>>()?;

    // Initialize heap with first record from each chunk
    let mut heap: BinaryHeap<std::cmp::Reverse<MergeEntry<K>>> =
        BinaryHeap::with_capacity(chunk_files.len());

    for reader in &chunk_readers {
        if let Some(record) = reader.next()? {
            let key = extract_key(&record);
            heap.push(std::cmp::Reverse(MergeEntry { key, record, chunk_idx: reader.idx }));
        }
    }

    // Create output writer with multi-threaded compression
    let mut writer = create_raw_bam_writer(
        output,
        output_header,
        config.writer_threads,
        config.compression_level,
    )?;

    let mut records_merged = 0u64;

    // Merge loop
    while let Some(std::cmp::Reverse(entry)) = heap.pop() {
        // Write record to output
        writer.write_raw_record(&entry.record)?;
        records_merged += 1;

        // Get next record from the same chunk (non-blocking due to prefetch buffer)
        let reader = &chunk_readers[entry.chunk_idx];
        if let Some(record) = reader.next()? {
            let key = extract_key(&record);
            heap.push(std::cmp::Reverse(MergeEntry { key, record, chunk_idx: entry.chunk_idx }));
        }
    }

    tracing::info!("Parallel merge complete: {records_merged} records merged");

    Ok(records_merged)
}

/// Parallel merge with output buffering for even higher throughput.
///
/// This version adds an output buffer that accumulates records before
/// writing them to the output file, reducing the number of write calls.
///
/// # Errors
///
/// Returns an error if reading chunks, writing output, or merging fails.
#[allow(clippy::needless_pass_by_value)]
pub fn parallel_merge_buffered<K, F>(
    chunk_files: &[PathBuf],
    _header: &Header,
    output_header: &Header,
    output: &Path,
    extract_key: F,
    config: ParallelMergeConfig,
) -> Result<u64>
where
    K: Clone + Send + Sync + Ord,
    F: Fn(&RawRecord) -> K + Send + Sync,
{
    const OUTPUT_BUFFER_SIZE: usize = 1024;

    tracing::info!(
        "Starting buffered parallel merge of {} chunks with {} reader threads",
        chunk_files.len(),
        config.reader_threads.min(chunk_files.len())
    );

    // Create prefetching readers for each chunk
    let chunk_readers: Vec<PrefetchingChunkReader> = chunk_files
        .iter()
        .enumerate()
        .map(|(idx, path)| PrefetchingChunkReader::new(path.clone(), idx))
        .collect::<Result<Vec<_>>>()?;

    // Initialize heap with first record from each chunk
    let mut heap: BinaryHeap<std::cmp::Reverse<MergeEntry<K>>> =
        BinaryHeap::with_capacity(chunk_files.len());

    for reader in &chunk_readers {
        if let Some(record) = reader.next()? {
            let key = extract_key(&record);
            heap.push(std::cmp::Reverse(MergeEntry { key, record, chunk_idx: reader.idx }));
        }
    }

    // Create output writer with multi-threaded compression
    let mut writer = create_raw_bam_writer(
        output,
        output_header,
        config.writer_threads,
        config.compression_level,
    )?;

    let mut records_merged = 0u64;
    let mut output_buffer: Vec<RawRecord> = Vec::with_capacity(OUTPUT_BUFFER_SIZE);

    // Merge loop with output buffering
    while let Some(std::cmp::Reverse(entry)) = heap.pop() {
        output_buffer.push(entry.record);
        records_merged += 1;

        // Flush buffer if full
        if output_buffer.len() >= OUTPUT_BUFFER_SIZE {
            for record in output_buffer.drain(..) {
                writer.write_raw_record(&record)?;
            }
        }

        // Get next record from the same chunk
        let reader = &chunk_readers[entry.chunk_idx];
        if let Some(record) = reader.next()? {
            let key = extract_key(&record);
            heap.push(std::cmp::Reverse(MergeEntry { key, record, chunk_idx: entry.chunk_idx }));
        }
    }

    // Flush remaining buffered records
    for record in output_buffer {
        writer.write_raw_record(&record)?;
    }

    tracing::info!("Buffered parallel merge complete: {records_merged} records merged");

    Ok(records_merged)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_merge_entry_ordering() {
        let entry1 = MergeEntry { key: 1, record: RawRecord::new(), chunk_idx: 0 };
        let entry2 = MergeEntry { key: 2, record: RawRecord::new(), chunk_idx: 1 };

        assert!(entry1 < entry2);
    }

    #[test]
    fn test_config_default() {
        let config = ParallelMergeConfig::default();
        assert_eq!(config.reader_threads, 4);
        assert_eq!(config.writer_threads, 4);
        assert_eq!(config.compression_level, 6);
    }

    #[test]
    fn test_merge_entry_equal_keys() {
        let entry1 = MergeEntry { key: 5, record: RawRecord::new(), chunk_idx: 0 };
        let entry2 = MergeEntry { key: 5, record: RawRecord::new(), chunk_idx: 1 };

        assert_eq!(entry1.cmp(&entry2), Ordering::Equal);
    }

    #[test]
    fn test_merge_entry_greater_than() {
        let entry1 = MergeEntry { key: 2, record: RawRecord::new(), chunk_idx: 0 };
        let entry2 = MergeEntry { key: 1, record: RawRecord::new(), chunk_idx: 1 };

        assert!(entry1 > entry2);
    }

    #[test]
    fn test_merge_entry_ordering_ignores_chunk_idx() {
        let entry1 = MergeEntry { key: 42, record: RawRecord::new(), chunk_idx: 0 };
        let entry2 = MergeEntry { key: 42, record: RawRecord::new(), chunk_idx: 99 };

        assert_eq!(entry1.cmp(&entry2), Ordering::Equal);
    }

    #[test]
    fn test_merge_entry_partial_eq() {
        let entry1 = MergeEntry { key: 10, record: RawRecord::new(), chunk_idx: 0 };
        let entry2 = MergeEntry { key: 10, record: RawRecord::new(), chunk_idx: 3 };

        assert!(entry1 == entry2);
    }

    #[test]
    fn test_merge_entry_partial_eq_different() {
        let entry1 = MergeEntry { key: 10, record: RawRecord::new(), chunk_idx: 0 };
        let entry2 = MergeEntry { key: 20, record: RawRecord::new(), chunk_idx: 0 };

        assert!(entry1 != entry2);
    }

    #[test]
    fn test_merge_entry_string_keys() {
        let entry_a =
            MergeEntry { key: "apple".to_string(), record: RawRecord::new(), chunk_idx: 0 };
        let entry_b =
            MergeEntry { key: "banana".to_string(), record: RawRecord::new(), chunk_idx: 1 };
        let entry_c =
            MergeEntry { key: "cherry".to_string(), record: RawRecord::new(), chunk_idx: 2 };

        assert!(entry_a < entry_b);
        assert!(entry_b < entry_c);
        assert!(entry_a < entry_c);
    }

    #[test]
    fn test_merge_entry_in_binary_heap() {
        use std::cmp::Reverse;

        let mut heap = BinaryHeap::new();
        heap.push(Reverse(MergeEntry { key: 3, record: RawRecord::new(), chunk_idx: 0 }));
        heap.push(Reverse(MergeEntry { key: 1, record: RawRecord::new(), chunk_idx: 1 }));
        heap.push(Reverse(MergeEntry { key: 2, record: RawRecord::new(), chunk_idx: 2 }));

        // Should come out in ascending order: 1, 2, 3
        assert_eq!(heap.pop().expect("heap should have elements").0.key, 1);
        assert_eq!(heap.pop().expect("heap should have elements").0.key, 2);
        assert_eq!(heap.pop().expect("heap should have elements").0.key, 3);
        assert!(heap.is_empty());
    }

    #[test]
    fn test_config_custom_values() {
        let config =
            ParallelMergeConfig { reader_threads: 8, writer_threads: 16, compression_level: 9 };

        assert_eq!(config.reader_threads, 8);
        assert_eq!(config.writer_threads, 16);
        assert_eq!(config.compression_level, 9);
    }

    #[test]
    fn test_config_single_thread() {
        let config =
            ParallelMergeConfig { reader_threads: 1, writer_threads: 1, compression_level: 1 };

        assert_eq!(config.reader_threads, 1);
        assert_eq!(config.writer_threads, 1);
        assert_eq!(config.compression_level, 1);
    }

    #[test]
    fn test_merge_buffer_size() {
        assert_eq!(MERGE_BUFFER_SIZE, 65536);
    }

    #[test]
    fn test_prefetch_buffer_size() {
        assert_eq!(PREFETCH_BUFFER_SIZE, 128);
    }
}

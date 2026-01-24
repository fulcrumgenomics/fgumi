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

use anyhow::{Context, Result};
use crossbeam_channel::{Receiver, Sender, bounded};
use noodles::bam::{self, Record};
use noodles::bgzf;
use noodles::sam::Header;
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};
use std::thread::{self, JoinHandle};

use crate::bam_io::create_bam_writer;

/// Buffer size for reading temp files during merge.
const MERGE_BUFFER_SIZE: usize = 64 * 1024;

/// Number of records to prefetch per chunk reader.
const PREFETCH_BUFFER_SIZE: usize = 128;

/// A record with its sort key and source chunk index.
pub struct MergeEntry<K> {
    pub key: K,
    pub record: Record,
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

/// Prefetching chunk reader that runs in a background thread.
///
/// This reader maintains a buffer of prefetched records, allowing the merge
/// thread to consume records without blocking on I/O.
struct PrefetchingChunkReader {
    /// Receiver for prefetched records.
    record_rx: Receiver<Option<Record>>,
    /// Handle to the reader thread.
    _handle: JoinHandle<()>,
    /// Chunk index for heap management.
    idx: usize,
}

impl PrefetchingChunkReader {
    /// Create a new prefetching reader for a chunk file.
    fn new(path: PathBuf, idx: usize) -> Result<Self> {
        // Channel for prefetched records
        let (record_tx, record_rx) = bounded(PREFETCH_BUFFER_SIZE);

        // Spawn reader thread
        let handle = thread::spawn(move || {
            if let Err(e) = Self::reader_thread(path, record_tx) {
                log::error!("Chunk reader thread failed: {e}");
            }
        });

        Ok(Self { record_rx, _handle: handle, idx })
    }

    /// Reader thread function.
    fn reader_thread(path: PathBuf, tx: Sender<Option<Record>>) -> Result<()> {
        let file = File::open(&path).context("Failed to open chunk file")?;
        let buf_reader = BufReader::with_capacity(MERGE_BUFFER_SIZE, file);
        let bgzf_reader = bgzf::io::Reader::new(buf_reader);
        let mut reader = bam::io::Reader::from(bgzf_reader);

        // Read and discard header
        reader.read_header()?;

        // Read records and send to channel
        let mut record = Record::default();
        loop {
            match reader.read_record(&mut record) {
                Ok(0) => {
                    // EOF - send None to signal end
                    let _ = tx.send(None);
                    break;
                }
                Ok(_) => {
                    // Clone and send record (take ownership of current, reset for next)
                    let owned_record = std::mem::take(&mut record);
                    if tx.send(Some(owned_record)).is_err() {
                        // Receiver dropped, exit
                        break;
                    }
                }
                Err(e) => {
                    log::error!("Error reading chunk: {e}");
                    let _ = tx.send(None);
                    break;
                }
            }
        }

        Ok(())
    }

    /// Get the next record from the prefetch buffer.
    fn next(&self) -> Option<Record> {
        match self.record_rx.recv() {
            Ok(Some(record)) => Some(record),
            Ok(None) | Err(_) => None,
        }
    }
}

/// Parallel merge implementation using prefetching readers.
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
    F: Fn(&Record) -> K + Send + Sync,
{
    log::info!(
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
        if let Some(record) = reader.next() {
            let key = extract_key(&record);
            heap.push(std::cmp::Reverse(MergeEntry { key, record, chunk_idx: reader.idx }));
        }
    }

    // Create output writer with multi-threaded compression
    let mut writer =
        create_bam_writer(output, output_header, config.writer_threads, config.compression_level)?;

    let mut records_merged = 0u64;

    // Merge loop
    while let Some(std::cmp::Reverse(entry)) = heap.pop() {
        // Write record to output
        writer.write_record(output_header, &entry.record)?;
        records_merged += 1;

        // Get next record from the same chunk (non-blocking due to prefetch buffer)
        let reader = &chunk_readers[entry.chunk_idx];
        if let Some(record) = reader.next() {
            let key = extract_key(&record);
            heap.push(std::cmp::Reverse(MergeEntry { key, record, chunk_idx: entry.chunk_idx }));
        }
    }

    log::info!("Parallel merge complete: {} records merged", records_merged);

    Ok(records_merged)
}

/// Parallel merge with output buffering for even higher throughput.
///
/// This version adds an output buffer that accumulates records before
/// writing them to the output file, reducing the number of write calls.
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
    F: Fn(&Record) -> K + Send + Sync,
{
    const OUTPUT_BUFFER_SIZE: usize = 1024;

    log::info!(
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
        if let Some(record) = reader.next() {
            let key = extract_key(&record);
            heap.push(std::cmp::Reverse(MergeEntry { key, record, chunk_idx: reader.idx }));
        }
    }

    // Create output writer with multi-threaded compression
    let mut writer =
        create_bam_writer(output, output_header, config.writer_threads, config.compression_level)?;

    let mut records_merged = 0u64;
    let mut output_buffer: Vec<Record> = Vec::with_capacity(OUTPUT_BUFFER_SIZE);

    // Merge loop with output buffering
    while let Some(std::cmp::Reverse(entry)) = heap.pop() {
        output_buffer.push(entry.record);
        records_merged += 1;

        // Flush buffer if full
        if output_buffer.len() >= OUTPUT_BUFFER_SIZE {
            for record in output_buffer.drain(..) {
                writer.write_record(output_header, &record)?;
            }
        }

        // Get next record from the same chunk
        let reader = &chunk_readers[entry.chunk_idx];
        if let Some(record) = reader.next() {
            let key = extract_key(&record);
            heap.push(std::cmp::Reverse(MergeEntry { key, record, chunk_idx: entry.chunk_idx }));
        }
    }

    // Flush remaining buffered records
    for record in output_buffer {
        writer.write_record(output_header, &record)?;
    }

    log::info!("Buffered parallel merge complete: {} records merged", records_merged);

    Ok(records_merged)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_merge_entry_ordering() {
        let entry1 = MergeEntry { key: 1, record: Record::default(), chunk_idx: 0 };
        let entry2 = MergeEntry { key: 2, record: Record::default(), chunk_idx: 1 };

        assert!(entry1 < entry2);
    }

    #[test]
    fn test_config_default() {
        let config = ParallelMergeConfig::default();
        assert_eq!(config.reader_threads, 4);
        assert_eq!(config.writer_threads, 4);
        assert_eq!(config.compression_level, 6);
    }
}

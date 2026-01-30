//! Background read-ahead for BAM sorting.
//!
//! This module provides a `ReadAheadReader` that reads BAM records in a background
//! thread while the main thread processes records. This overlaps I/O with computation
//! for improved throughput.
//!
//! # Architecture
//!
//! ```text
//! ┌─────────────────┐    ┌─────────────────┐    ┌─────────────────┐
//! │  BGZF Reader    │───>│ Prefetch Buffer │───>│  Main Thread    │
//! │ (background)    │    │   (channel)     │    │  (processing)   │
//! └─────────────────┘    └─────────────────┘    └─────────────────┘
//! ```
//!
//! # Batched Reading
//!
//! Records are sent in batches (default 256 records per batch) to reduce
//! channel synchronization overhead by ~256x compared to per-record sends.

use crossbeam_channel::{Receiver, Sender, bounded};
use noodles::bam::Record;
use std::thread::{self, JoinHandle};

use crate::bam_io::BamReaderAuto;

/// Number of records per batch sent through the channel.
/// This reduces channel synchronization overhead by sending multiple records per send.
const BATCH_SIZE: usize = 256;

/// Number of batches to buffer in the channel.
/// Total prefetch = `BATCH_SIZE` * `CHANNEL_BUFFER_SIZE` = 256 * 16 = 4096 records.
const CHANNEL_BUFFER_SIZE: usize = 16;

/// Background reader that reads ahead while main thread processes.
///
/// This reader spawns a background thread that continuously reads records
/// from the BAM file and sends them through a channel in batches. The main
/// thread can then receive records without blocking on I/O.
///
/// Records are batched (256 per batch by default) to reduce channel
/// synchronization overhead by ~256x compared to per-record sends.
pub struct ReadAheadReader {
    /// Receiver for prefetched record batches (Option to allow closing before join).
    receiver: Option<Receiver<Vec<Record>>>,
    /// Handle to the reader thread.
    handle: Option<JoinHandle<()>>,
    /// Current batch being consumed.
    current_batch: Vec<Record>,
    /// Index into current batch.
    batch_index: usize,
}

impl ReadAheadReader {
    /// Create a new read-ahead reader with default buffer size.
    #[must_use]
    pub fn new(reader: BamReaderAuto) -> Self {
        Self::with_batch_size(reader, BATCH_SIZE, CHANNEL_BUFFER_SIZE)
    }

    /// Create a new read-ahead reader with specified batch and buffer sizes.
    #[must_use]
    pub fn with_batch_size(
        mut reader: BamReaderAuto,
        batch_size: usize,
        channel_buffer: usize,
    ) -> Self {
        let (tx, rx) = bounded(channel_buffer);

        let handle = thread::spawn(move || {
            Self::reader_thread(&mut reader, tx, batch_size);
        });

        Self { receiver: Some(rx), handle: Some(handle), current_batch: Vec::new(), batch_index: 0 }
    }

    /// Create with legacy per-record behavior (for compatibility).
    #[must_use]
    pub fn with_buffer_size(reader: BamReaderAuto, buffer_size: usize) -> Self {
        // Approximate the old behavior: buffer_size records total
        let buffer_size = buffer_size.max(1);
        let batch_size = BATCH_SIZE.min(buffer_size);
        let channel_buffer = buffer_size.div_ceil(batch_size);
        Self::with_batch_size(reader, batch_size, channel_buffer)
    }

    /// Reader thread function - reads records in batches.
    fn reader_thread(reader: &mut BamReaderAuto, tx: Sender<Vec<Record>>, batch_size: usize) {
        let mut record = Record::default();
        let mut batch = Vec::with_capacity(batch_size);

        loop {
            match reader.read_record(&mut record) {
                Ok(0) => {
                    // EOF - send any remaining records
                    if !batch.is_empty() {
                        let _ = tx.send(batch);
                    }
                    // Send empty batch to signal EOF
                    let _ = tx.send(Vec::new());
                    break;
                }
                Ok(_) => {
                    // Successfully read a record
                    batch.push(std::mem::take(&mut record));

                    // Send batch when full
                    if batch.len() >= batch_size {
                        if tx.send(batch).is_err() {
                            // Receiver dropped
                            break;
                        }
                        batch = Vec::with_capacity(batch_size);
                    }
                }
                Err(e) => {
                    log::error!("Error reading BAM record: {e}");
                    // Send any accumulated records before error
                    if !batch.is_empty() {
                        let _ = tx.send(batch);
                    }
                    let _ = tx.send(Vec::new());
                    break;
                }
            }
        }
    }

    /// Get the next record from the read-ahead buffer.
    ///
    /// Returns `Some(record)` if a record is available, `None` on EOF or error.
    #[inline]
    #[must_use]
    pub fn next_record(&mut self) -> Option<Record> {
        // Check if we have records in the current batch
        if self.batch_index < self.current_batch.len() {
            let record = std::mem::take(&mut self.current_batch[self.batch_index]);
            self.batch_index += 1;
            return Some(record);
        }

        // Need to fetch next batch
        let receiver = self.receiver.as_ref()?;
        match receiver.recv() {
            Ok(batch) if !batch.is_empty() => {
                self.current_batch = batch;
                self.batch_index = 1; // We'll return index 0
                Some(std::mem::take(&mut self.current_batch[0]))
            }
            Ok(_) | Err(_) => None, // Empty batch = EOF, or channel closed
        }
    }
}

impl Drop for ReadAheadReader {
    fn drop(&mut self) {
        // Close the receiver first to unblock any pending sends in the reader thread.
        // This prevents deadlock when the channel is full and the thread is blocked on send.
        drop(self.receiver.take());

        // Now wait for the reader thread to finish
        if let Some(handle) = self.handle.take() {
            let _ = handle.join();
        }
    }
}

/// Iterator adapter for `ReadAheadReader`.
impl Iterator for ReadAheadReader {
    type Item = Record;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_record()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_batch_config() {
        // Default: 256 records/batch × 16 batches = 4096 records total prefetch
        assert_eq!(BATCH_SIZE, 256);
        assert_eq!(CHANNEL_BUFFER_SIZE, 16);
        assert_eq!(BATCH_SIZE * CHANNEL_BUFFER_SIZE, 4096);
    }
}

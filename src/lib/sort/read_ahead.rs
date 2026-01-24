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

use crossbeam_channel::{Receiver, Sender, bounded};
use noodles::bam::Record;
use std::thread::{self, JoinHandle};

use crate::bam_io::BamReaderAuto;

/// Default number of records to prefetch.
const DEFAULT_PREFETCH_SIZE: usize = 4096;

/// Background reader that reads ahead while main thread processes.
///
/// This reader spawns a background thread that continuously reads records
/// from the BAM file and sends them through a channel. The main thread
/// can then receive records without blocking on I/O.
pub struct ReadAheadReader {
    /// Receiver for prefetched records.
    receiver: Receiver<Option<Record>>,
    /// Handle to the reader thread.
    handle: Option<JoinHandle<()>>,
}

impl ReadAheadReader {
    /// Create a new read-ahead reader with default buffer size.
    #[must_use]
    pub fn new(reader: BamReaderAuto) -> Self {
        Self::with_buffer_size(reader, DEFAULT_PREFETCH_SIZE)
    }

    /// Create a new read-ahead reader with specified buffer size.
    #[must_use]
    pub fn with_buffer_size(mut reader: BamReaderAuto, buffer_size: usize) -> Self {
        let (tx, rx) = bounded(buffer_size);

        let handle = thread::spawn(move || {
            Self::reader_thread(&mut reader, tx);
        });

        Self { receiver: rx, handle: Some(handle) }
    }

    /// Reader thread function.
    fn reader_thread(reader: &mut BamReaderAuto, tx: Sender<Option<Record>>) {
        let mut record = Record::default();

        loop {
            match reader.read_record(&mut record) {
                Ok(0) => {
                    // EOF
                    let _ = tx.send(None);
                    break;
                }
                Ok(_) => {
                    // Successfully read a record
                    let owned_record = std::mem::take(&mut record);
                    if tx.send(Some(owned_record)).is_err() {
                        // Receiver dropped
                        break;
                    }
                }
                Err(e) => {
                    log::error!("Error reading BAM record: {e}");
                    let _ = tx.send(None);
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
    pub fn next_record(&self) -> Option<Record> {
        match self.receiver.recv() {
            Ok(Some(record)) => Some(record),
            Ok(None) | Err(_) => None,
        }
    }
}

impl Drop for ReadAheadReader {
    fn drop(&mut self) {
        // Wait for the reader thread to finish
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
    fn test_default_prefetch_size() {
        assert_eq!(DEFAULT_PREFETCH_SIZE, 4096);
    }
}

//! Background read-ahead for BAM sorting.
//!
//! This module provides readers that prefetch BAM records in a background
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
//! # Readers
//!
//! - [`RawReadAheadReader`]: Yields [`RawRecord`] (raw bytes, no noodles Record),
//!   with I/O-error propagation via `take_error()`
//!
//! # Batched Reading
//!
//! Records are sent in batches (default 256 records per batch) to reduce
//! channel synchronization overhead by ~256x compared to per-record sends.

use crossbeam_channel::{Receiver, Sender, bounded};
use std::thread::{self, JoinHandle};

use fgumi_bam_io::RawBamReaderAuto;
use fgumi_raw_bam::{RawBamReader, RawRecord};
use std::io::Read as IoRead;

/// Number of records per batch sent through the channel.
/// This reduces channel synchronization overhead by sending multiple records per send.
const BATCH_SIZE: usize = 256;

/// Number of batches to buffer in the channel.
/// Total prefetch = `BATCH_SIZE` * `CHANNEL_BUFFER_SIZE` = 256 * 16 = 4096 records.
const CHANNEL_BUFFER_SIZE: usize = 16;

/// Background reader that reads raw BAM bytes ahead while main thread processes.
///
/// Yields [`RawRecord`] (raw bytes) instead of `noodles::bam::Record`. This enables:
/// - Zero-copy field extraction using fixed byte offsets
/// - Direct writes to output without re-encoding
/// - No dependency on noodles' internal Record representation
///
/// Use this for high-performance raw sorting pipelines.
pub struct RawReadAheadReader {
    /// Receiver for prefetched record batches (Option to allow closing before join).
    receiver: Option<Receiver<Vec<RawRecord>>>,
    /// Handle to the reader thread.
    handle: Option<JoinHandle<()>>,
    /// Current batch being consumed.
    current_batch: Vec<RawRecord>,
    /// Index into current batch.
    batch_index: usize,
    /// Captures any I/O error from the background reader thread.
    ///
    /// The background thread stores its error here before sending the EOF sentinel,
    /// so callers can retrieve it via `take_error()` after iteration.
    error_slot: std::sync::Arc<std::sync::Mutex<Option<std::io::Error>>>,
    /// Set once the reader thread's EOF sentinel has been received.
    ///
    /// Distinguishes the two ways `recv()` stops yielding records. The reader
    /// thread sends the sentinel on every exit path it controls — clean EOF and
    /// read error alike — so a sender dropped *without* one means the thread died
    /// (a panic), leaving the stream truncated at an arbitrary record. Without
    /// this flag both collapse to `None`, and a caller that trusts `take_error()`
    /// reads a partial input as the whole thing.
    saw_eof: bool,
}

impl RawReadAheadReader {
    /// Create a new raw read-ahead reader with default buffer size.
    #[must_use]
    pub fn new(reader: RawBamReaderAuto) -> Self {
        Self::from_reader(reader)
    }

    /// Create a raw read-ahead reader from any `RawBamReader<R>`.
    ///
    /// This accepts any `RawBamReader<R>` and spawns a background thread that
    /// reads records in batches.
    ///
    #[must_use]
    pub(crate) fn from_reader<R: IoRead + Send + 'static>(reader: RawBamReader<R>) -> Self {
        Self::from_reader_with_batch_size(reader, BATCH_SIZE, CHANNEL_BUFFER_SIZE)
    }

    /// Create a raw read-ahead reader from any `RawBamReader<R>` with specified
    /// batch and buffer sizes.
    #[must_use]
    pub(crate) fn from_reader_with_batch_size<R: IoRead + Send + 'static>(
        mut reader: RawBamReader<R>,
        batch_size: usize,
        channel_buffer: usize,
    ) -> Self {
        let error_slot = std::sync::Arc::new(std::sync::Mutex::new(None));
        let error_slot_thread = std::sync::Arc::clone(&error_slot);
        let (tx, rx) = bounded(channel_buffer);

        let handle = thread::spawn(move || {
            Self::reader_thread_generic(&mut reader, tx, batch_size, error_slot_thread);
        });

        Self {
            receiver: Some(rx),
            handle: Some(handle),
            current_batch: Vec::new(),
            batch_index: 0,
            error_slot,
            saw_eof: false,
        }
    }

    /// Reader thread function - reads raw records in batches from any reader type.
    #[allow(clippy::needless_pass_by_value)]
    fn reader_thread_generic<R: IoRead>(
        reader: &mut RawBamReader<R>,
        tx: Sender<Vec<RawRecord>>,
        batch_size: usize,
        error_slot: std::sync::Arc<std::sync::Mutex<Option<std::io::Error>>>,
    ) {
        let mut record = RawRecord::new();
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
                    // Successfully read a record - take ownership and replace with empty
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
                    log::error!("Error reading raw BAM record: {e}");
                    // Store error so the foreground thread can retrieve it via take_error()
                    if let Ok(mut slot) = error_slot.lock() {
                        *slot = Some(e);
                    }
                    // Send any accumulated records before the EOF sentinel
                    if !batch.is_empty() {
                        let _ = tx.send(batch);
                    }
                    let _ = tx.send(Vec::new());
                    break;
                }
            }
        }
    }

    /// Take any I/O error that occurred in the background reader thread.
    ///
    /// Returns `Some(err)` if the background thread stopped due to a read error
    /// rather than clean EOF. Call this after exhausting the iterator.
    #[must_use]
    pub fn take_error(&self) -> Option<std::io::Error> {
        self.error_slot.lock().ok()?.take()
    }

    /// Park an error reporting that the reader thread stopped without reaching EOF.
    ///
    /// A real read error is stored by the reader thread itself before it sends the
    /// sentinel, so an occupied slot always describes the underlying failure more
    /// precisely than this does and is left alone.
    fn park_reader_died_error(&self) {
        if let Ok(mut slot) = self.error_slot.lock()
            && slot.is_none()
        {
            *slot = Some(std::io::Error::new(
                std::io::ErrorKind::UnexpectedEof,
                "read-ahead reader thread stopped before the end of the input; \
                 the records read so far are only part of it",
            ));
        }
    }

    /// Get the next raw record from the read-ahead buffer.
    ///
    /// Returns `Some(record)` if a record is available, `None` on EOF or error.
    ///
    /// A `None` from a stream that never delivered its EOF sentinel parks an error
    /// for [`take_error`](Self::take_error), so a caller cannot mistake a dead
    /// reader thread for the end of the input.
    #[inline]
    #[must_use]
    pub fn next_record(&mut self) -> Option<RawRecord> {
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
            // Empty batch: the reader thread's EOF sentinel, sent on both its clean
            // and its failing exit path.
            Ok(_) => {
                self.saw_eof = true;
                None
            }
            // Channel closed. After the sentinel that is just a second poll of an
            // exhausted stream; before it, the reader thread died mid-input.
            Err(_) => {
                if !self.saw_eof {
                    self.park_reader_died_error();
                }
                None
            }
        }
    }
}

impl Drop for RawReadAheadReader {
    fn drop(&mut self) {
        // Close the receiver first to unblock any pending sends in the reader thread.
        drop(self.receiver.take());

        // Now wait for the reader thread to finish
        if let Some(handle) = self.handle.take() {
            let _ = handle.join();
        }
    }
}

/// Iterator adapter for `RawReadAheadReader`.
impl Iterator for RawReadAheadReader {
    type Item = RawRecord;

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

    #[test]
    fn test_prefetch_total_records() {
        // Verify the total prefetch capacity is exactly 4096 records.
        // This is the product of batch size and channel buffer size, and is
        // documented in the module-level comments as the expected total.
        let total = BATCH_SIZE * CHANNEL_BUFFER_SIZE;
        assert_eq!(total, 4096, "Total prefetch should be BATCH_SIZE * CHANNEL_BUFFER_SIZE = 4096");
    }

    use noodles::sam::Header;
    use noodles::sam::alignment::RecordBuf;
    use noodles::sam::alignment::io::Write as AlignmentWrite;
    use noodles::sam::alignment::record::Flags;
    use noodles::sam::header::record::value::Map;
    use noodles::sam::header::record::value::map::ReferenceSequence;
    use std::num::NonZeroUsize;
    use tempfile::NamedTempFile;

    use fgumi_bam_io::create_raw_bam_reader;

    /// Create a temporary BAM file with the given number of unmapped records.
    /// Returns the temp file handle (keeps the file alive) and the header used.
    fn create_test_bam_file(num_records: usize) -> (NamedTempFile, Header) {
        let ref_seq =
            Map::<ReferenceSequence>::new(NonZeroUsize::new(1000).expect("1000 is non-zero"));
        let header = Header::builder().add_reference_sequence("chr1", ref_seq).build();

        let tmp = NamedTempFile::new().expect("creating temp file/dir should succeed");
        let path = tmp.path().to_path_buf();

        {
            let file = std::fs::File::create(&path).expect("creating file should succeed");
            let mut writer = noodles::bam::io::Writer::new(file);
            writer.write_header(&header).expect("writing header should succeed");

            for i in 0..num_records {
                let name = format!("read{i}");
                let record =
                    RecordBuf::builder().set_name(&*name).set_flags(Flags::UNMAPPED).build();
                writer
                    .write_alignment_record(&header, &record)
                    .expect("writing alignment record should succeed");
            }
        }

        (tmp, header)
    }

    #[test]
    fn test_raw_read_ahead_empty() {
        let (tmp, _header) = create_test_bam_file(0);
        let (reader, _header) =
            create_raw_bam_reader(tmp.path(), 1).expect("creating BAM reader should succeed");

        let mut ra = RawReadAheadReader::new(reader);
        assert!(ra.next_record().is_none(), "Empty BAM should yield no raw records");
    }

    #[test]
    fn test_raw_read_ahead_multiple() {
        let num = 10;
        let (tmp, _header) = create_test_bam_file(num);
        let (reader, _header) =
            create_raw_bam_reader(tmp.path(), 1).expect("creating BAM reader should succeed");

        let ra = RawReadAheadReader::new(reader);
        let records: Vec<RawRecord> = ra.collect();
        assert_eq!(records.len(), num, "Raw read-ahead should yield exactly {num} records");
    }

    /// Build a `RawRecord` whose bytes are exactly `body`.
    fn raw_record_from(body: &[u8]) -> RawRecord {
        let mut rec = RawRecord::new();
        rec.as_mut_vec().extend_from_slice(body);
        rec
    }

    /// A sender dropped without the EOF sentinel is a dead reader thread, not the
    /// end of the input, and must not read as a clean EOF.
    ///
    /// The reader thread sends the sentinel on both exit paths it controls, so the
    /// only way the channel closes without one is a panic — after which the batches
    /// already sent are an arbitrary prefix of the input. `take_error()` returning
    /// `None` there would let `CheckedRecords` report clean EOF, and a comparison
    /// would call two inputs identical having read only part of one.
    ///
    /// The channel is driven directly rather than through a panicking reader: the
    /// state under test is the receiver's, and a real panic would only add noise to
    /// the test output.
    #[test]
    fn test_read_ahead_reports_a_reader_thread_that_died_before_eof() {
        let (tx, rx) = bounded::<Vec<RawRecord>>(4);
        tx.send(vec![raw_record_from(&[1u8; 8])]).expect("send one batch");
        drop(tx); // The thread died here: no EOF sentinel follows.

        let mut reader = RawReadAheadReader {
            receiver: Some(rx),
            handle: None,
            current_batch: Vec::new(),
            batch_index: 0,
            error_slot: std::sync::Arc::new(std::sync::Mutex::new(None)),
            saw_eof: false,
        };

        let seen: Vec<RawRecord> = std::iter::from_fn(|| reader.next()).collect();

        assert_eq!(seen.len(), 1, "the batch already sent must still be yielded");
        let err = reader.take_error().expect("a disconnect before EOF must park an error");
        assert_eq!(
            err.kind(),
            std::io::ErrorKind::UnexpectedEof,
            "a truncated stream is an unexpected EOF, got: {err}"
        );
    }

    /// The sentinel ends the stream cleanly, and polling an already-ended stream
    /// again must not manufacture the error above.
    ///
    /// Iterators are polled past their end routinely — `Iterator::count` and
    /// `for` both do — and the channel is closed by then, so the disconnect arm
    /// has to tell "closed after the sentinel" from "closed instead of it".
    #[test]
    fn test_read_ahead_polled_past_the_eof_sentinel_reports_no_error() {
        let (tx, rx) = bounded::<Vec<RawRecord>>(4);
        tx.send(vec![raw_record_from(&[1u8; 8])]).expect("send one batch");
        tx.send(Vec::new()).expect("send the EOF sentinel");
        drop(tx);

        let mut reader = RawReadAheadReader {
            receiver: Some(rx),
            handle: None,
            current_batch: Vec::new(),
            batch_index: 0,
            error_slot: std::sync::Arc::new(std::sync::Mutex::new(None)),
            saw_eof: false,
        };

        let seen: Vec<RawRecord> = std::iter::from_fn(|| reader.next()).collect();
        assert_eq!(seen.len(), 1, "one record precedes the sentinel");
        assert!(reader.next().is_none(), "an ended stream keeps returning None");

        assert!(
            reader.take_error().is_none(),
            "a stream that delivered its sentinel ended cleanly, however often it is polled"
        );
    }
}

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
//! - [`ReadAheadReader`]: Yields `noodles::bam::Record` (uses noodles parsing)
//! - [`RawReadAheadReader`]: Yields `RawRecord` (raw bytes, no noodles Record)
//!
//! # Batched Reading
//!
//! Records are sent in batches (default 256 records per batch) to reduce
//! channel synchronization overhead by ~256x compared to per-record sends.

use crossbeam_channel::{Receiver, Sender, bounded};
use noodles::bam::Record;
use std::thread::{self, JoinHandle};

use crate::bam_io::{BamReaderAuto, RawBamReaderAuto};
use fgumi_raw_bam::{RawBamReader, RawRecord};
use std::io::Read as IoRead;

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
    ///
    /// # Error handling
    ///
    /// Unlike [`RawReadAheadReader`], this reader does **not** propagate I/O errors
    /// to the foreground thread. When a read fails the error is logged and the thread
    /// sends an EOF sentinel, making the error indistinguishable from clean EOF.
    ///
    /// This is intentional: `ReadAheadReader` is used in paths (e.g. group/correct)
    /// where errors are surfaced by the noodles layer before they reach here, and
    /// adding an `error_slot` would complicate the public API for no practical gain.
    /// If you need error propagation, use [`RawReadAheadReader`] which provides
    /// `take_error()`.
    #[allow(clippy::needless_pass_by_value)]
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

// ============================================================================
// RawReadAheadReader - Raw bytes version (no noodles Record dependency)
// ============================================================================

/// Background reader that reads raw BAM bytes ahead while main thread processes.
///
/// This is similar to [`ReadAheadReader`] but yields [`RawRecord`] (raw bytes)
/// instead of `noodles::bam::Record`. This enables:
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
    /// # Note
    ///
    /// Do not pass a `RawBamReader<PooledInputStream>` here: the pool workers only
    /// unpark the pool's stored main-thread handle, which would cause a deadlock if
    /// the consumer thread is the background reader thread instead of the main thread.
    /// Use `RecordSource::direct` for pool-integrated reading.
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

    /// Get the next raw record from the read-ahead buffer.
    ///
    /// Returns `Some(record)` if a record is available, `None` on EOF or error.
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
            Ok(_) | Err(_) => None, // Empty batch = EOF, or channel closed
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

// ============================================================================
// PooledInputStream — pool-based decompressed input stream (no extra threads)
// ============================================================================

/// A `Read` implementation that consumes decompressed blocks from the pool's
/// `decompressed_input` `ArrayQueue` and presents a contiguous byte stream.
///
/// Workers do `ReadInputBlocks` + `DecompressInput`; this struct reassembles
/// the decompressed blocks in order for the main thread to parse records from.
///
/// Uses non-blocking `ArrayQueue::pop()` with `park()`/`unpark()` notification.
/// Workers call `unpark()` on the main thread after pushing blocks, so the main
/// thread sleeps at zero CPU when waiting and wakes instantly when data arrives.
/// EOF is detected via the shared `decompressed_input_done` `AtomicBool` flag.
pub struct PooledInputStream {
    /// `ArrayQueue` to pop decompressed blocks from pool workers.
    decompressed_input: std::sync::Arc<crossbeam_queue::ArrayQueue<(u64, Vec<u8>)>>,
    /// Shared flag: set when all input blocks have been decompressed.
    decompressed_input_done: std::sync::Arc<std::sync::atomic::AtomicBool>,
    /// Shared flag: set when a worker encountered an I/O error reading the input file.
    input_read_error: std::sync::Arc<std::sync::atomic::AtomicBool>,
    /// Shared flag: set when a worker failed to decompress a BGZF block.
    decompression_error: std::sync::Arc<std::sync::atomic::AtomicBool>,
    /// Reorder buffer: holds out-of-order decompressed blocks until their serial is next.
    reorder: crate::reorder_buffer::ReorderBuffer<Vec<u8>>,
    /// Current buffer being read from.
    current_buf: Vec<u8>,
    /// Read position within `current_buf`.
    current_pos: usize,
}

impl PooledInputStream {
    /// Create a new pooled input stream from the pool's shared state.
    #[must_use]
    pub fn new(
        decompressed_input: std::sync::Arc<crossbeam_queue::ArrayQueue<(u64, Vec<u8>)>>,
        decompressed_input_done: std::sync::Arc<std::sync::atomic::AtomicBool>,
        input_read_error: std::sync::Arc<std::sync::atomic::AtomicBool>,
        decompression_error: std::sync::Arc<std::sync::atomic::AtomicBool>,
    ) -> Self {
        Self {
            decompressed_input,
            decompressed_input_done,
            input_read_error,
            decompression_error,
            reorder: crate::reorder_buffer::ReorderBuffer::new(),
            current_buf: Vec::new(),
            current_pos: 0,
        }
    }

    /// Check if all input has been decompressed and the queue is drained.
    fn is_eof(&self) -> bool {
        self.decompressed_input_done.load(std::sync::atomic::Ordering::Acquire)
            && self.decompressed_input.is_empty()
    }

    /// Drain available blocks from the `ArrayQueue` into the reorder buffer.
    ///
    /// Applies a soft cap (`2 * queue.capacity()`) on the reorder buffer so
    /// a temporarily fast producer cannot build an unbounded backlog, while
    /// still draining unconditionally when the next expected sequence has not
    /// yet landed (otherwise the pop that would unblock the main thread could
    /// remain stuck in the `ArrayQueue` and deadlock the pipeline).
    fn drain_queue(&mut self) {
        // Deadlock-free soft cap on the reorder buffer to apply real
        // backpressure upstream. Without a cap, drain_queue unconditionally
        // empties the bounded ArrayQueue into the unbounded ReorderBuffer,
        // allowing it to grow to tens of GB when the main thread consumes
        // slower than workers produce.
        //
        // Invariant: if we *cannot* currently pop (next_seq is missing), we
        // must keep draining unconditionally — otherwise next_seq, which may
        // still be sitting in the ArrayQueue, can never be transferred, and
        // we deadlock.
        let reorder_cap = self.decompressed_input.capacity() * 2;
        loop {
            // Only enforce the cap when forward progress is possible.
            if self.reorder.can_pop() && self.reorder.buffer_len() >= reorder_cap {
                break;
            }
            match self.decompressed_input.pop() {
                Some((serial, data)) => self.reorder.insert(serial, data),
                None => break,
            }
        }
    }

    /// Get the next decompressed block in serial order.
    ///
    /// Drains the `ArrayQueue` into a reorder buffer, then checks if the next
    /// expected serial is available. If not, parks the thread until a worker
    /// calls `unpark()` after pushing new data.
    fn next_block(&mut self) -> Option<Vec<u8>> {
        loop {
            // Drain everything available into the reorder buffer
            self.drain_queue();

            // Check if the block we need is ready (O(1) via ReorderBuffer)
            if let Some(data) = self.reorder.try_pop_next() {
                return Some(data);
            }

            // Nothing available — check EOF
            if self.is_eof() {
                // Re-drain one more time before declaring EOF: a block could have
                // been pushed between the previous drain and the is_eof() check.
                // `decompressed_input_done` is only set after the block is in the
                // queue, so the risk is low, but a second drain makes this watertight.
                self.drain_queue();
                if let Some(data) = self.reorder.try_pop_next() {
                    return Some(data);
                }
                return None;
            }

            // Park until a worker pushes a block and calls unpark()
            std::thread::park();

            // After waking, check for errors before looping back to drain.
            // A worker may have set an error flag instead of pushing a block.
            if self.input_read_error.load(std::sync::atomic::Ordering::Acquire)
                || self.decompression_error.load(std::sync::atomic::Ordering::Acquire)
            {
                return None;
            }
        }
    }
}

impl IoRead for PooledInputStream {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        // Serve from current buffer first
        if self.current_pos < self.current_buf.len() {
            let available = &self.current_buf[self.current_pos..];
            let n = available.len().min(buf.len());
            buf[..n].copy_from_slice(&available[..n]);
            self.current_pos += n;
            return Ok(n);
        }

        // Check for I/O or decompression error before blocking on next_block()
        if self.input_read_error.load(std::sync::atomic::Ordering::Acquire) {
            return Err(std::io::Error::other(
                "I/O error reading input BAM blocks (see log for details)",
            ));
        }
        if self.decompression_error.load(std::sync::atomic::Ordering::Acquire) {
            return Err(std::io::Error::other(
                "BGZF decompression error on input blocks (see log for details)",
            ));
        }

        // Get next block
        if let Some(data) = self.next_block() {
            let n = data.len().min(buf.len());
            buf[..n].copy_from_slice(&data[..n]);
            if n < data.len() {
                self.current_buf = data;
                self.current_pos = n;
            } else {
                self.current_buf.clear();
                self.current_pos = 0;
            }
            Ok(n)
        } else {
            // Double-check error flags after draining (an error may have caused early EOF)
            if self.input_read_error.load(std::sync::atomic::Ordering::Acquire) {
                return Err(std::io::Error::other(
                    "I/O error reading input BAM blocks (see log for details)",
                ));
            }
            if self.decompression_error.load(std::sync::atomic::Ordering::Acquire) {
                return Err(std::io::Error::other(
                    "BGZF decompression error on input blocks (see log for details)",
                ));
            }
            Ok(0) // clean EOF
        }
    }
}

// ============================================================================
// RecordSource — unified iterator for pool and non-pool paths
// ============================================================================

/// Unified record source for Phase 1 reading.
///
/// Wraps either a `RawReadAheadReader` (non-pool path, has background thread)
/// or a direct `RawBamReader<PooledInputStream>` (pool path, no extra threads).
pub enum RecordSource {
    /// Legacy path: background thread prefetches records.
    ReadAhead(RawReadAheadReader),
    /// Pool path: main thread reads directly from pool's decompressed stream.
    ///
    /// The second field stores the first I/O error encountered during iteration,
    /// if any. Callers should call `take_error()` after the iteration loop to
    /// propagate errors rather than silently treating them as EOF.
    Direct(RawBamReader<PooledInputStream>, Option<std::io::Error>),
}

impl RecordSource {
    /// Wrap a pooled reader as the `Direct` variant.
    #[must_use]
    pub fn direct(reader: RawBamReader<PooledInputStream>) -> Self {
        Self::Direct(reader, None)
    }

    /// Take any I/O error that occurred during iteration.
    ///
    /// Returns `Some(err)` if the iterator stopped due to a read error rather than
    /// clean EOF. Call this after exhausting the iterator to detect truncated input.
    pub fn take_error(&mut self) -> Option<std::io::Error> {
        match self {
            Self::Direct(_, err) => err.take(),
            Self::ReadAhead(r) => r.take_error(),
        }
    }
}

impl Iterator for RecordSource {
    type Item = RawRecord;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            Self::ReadAhead(r) => r.next(),
            Self::Direct(reader, error_slot) => {
                let mut record = RawRecord::default();
                match reader.read_record(&mut record) {
                    Ok(0) => None, // EOF
                    Ok(_) => Some(record),
                    Err(e) => {
                        log::error!("Error reading raw BAM record: {e}");
                        // Preserve the first error; don't overwrite with later ones.
                        if error_slot.is_none() {
                            *error_slot = Some(e);
                        }
                        None
                    }
                }
            }
        }
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

    use crate::bam_io::{create_bam_reader, create_raw_bam_reader};

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
    fn test_read_ahead_reader_empty_bam() {
        let (tmp, _header) = create_test_bam_file(0);
        let (reader, _header) =
            create_bam_reader(tmp.path(), 1).expect("creating BAM reader should succeed");

        let mut ra = ReadAheadReader::new(reader);
        assert!(ra.next_record().is_none(), "Empty BAM should yield no records");
    }

    #[test]
    fn test_read_ahead_reader_single_record() {
        let (tmp, _header) = create_test_bam_file(1);
        let (reader, _header) =
            create_bam_reader(tmp.path(), 1).expect("creating BAM reader should succeed");

        let mut ra = ReadAheadReader::new(reader);
        assert!(ra.next_record().is_some(), "Should yield exactly one record");
        assert!(ra.next_record().is_none(), "Should yield no more records after the single record");
    }

    #[test]
    fn test_read_ahead_reader_multiple_records() {
        let num = 10;
        let (tmp, _header) = create_test_bam_file(num);
        let (reader, _header) =
            create_bam_reader(tmp.path(), 1).expect("creating BAM reader should succeed");

        let mut ra = ReadAheadReader::new(reader);
        let mut count = 0;
        while ra.next_record().is_some() {
            count += 1;
        }
        assert_eq!(count, num, "Should read exactly {num} records");
    }

    #[test]
    fn test_read_ahead_reader_iterator() {
        let num = 10;
        let (tmp, _header) = create_test_bam_file(num);
        let (reader, _header) =
            create_bam_reader(tmp.path(), 1).expect("creating BAM reader should succeed");

        let ra = ReadAheadReader::new(reader);
        let records: Vec<Record> = ra.collect();
        assert_eq!(records.len(), num, "Iterator should collect exactly {num} records");
    }

    #[test]
    fn test_read_ahead_reader_with_buffer_size() {
        let num = 10;
        let (tmp, _header) = create_test_bam_file(num);
        let (reader, _header) =
            create_bam_reader(tmp.path(), 1).expect("creating BAM reader should succeed");

        // Use a small buffer size (4) to exercise partial-batch logic
        let ra = ReadAheadReader::with_buffer_size(reader, 4);
        let records: Vec<Record> = ra.collect();
        assert_eq!(records.len(), num, "with_buffer_size(4) should still read all {num} records");
    }

    #[test]
    fn test_read_ahead_reader_drop_while_reading() {
        // Create many records so the background thread is likely still reading
        // when we drop the reader. This verifies the Drop impl does not deadlock.
        let num = 5000;
        let (tmp, _header) = create_test_bam_file(num);
        let (reader, _header) =
            create_bam_reader(tmp.path(), 1).expect("creating BAM reader should succeed");

        let mut ra = ReadAheadReader::new(reader);
        // Read only a few records
        for _ in 0..5 {
            let _ = ra.next_record();
        }
        // Drop the reader while there are still records pending in the channel
        // or still being read by the background thread.
        drop(ra);
        // If we reach here, the Drop impl did not deadlock.
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
}

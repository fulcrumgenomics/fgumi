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
    ///
    /// This only reliably reflects background errors when iteration ran to its
    /// natural EOF (`next_record` returned `None`). If the consumer stops early
    /// and drops the reader before EOF, `Drop` drops the channel receiver, which
    /// can make the background thread's `tx.send` fail so it breaks out of its
    /// loop *before* recording any subsequent read error — that error is then
    /// lost. Do not rely on `take_error()` after an early `break`.
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
    reorder: fgumi_bam_io::ReorderBuffer<Vec<u8>>,
    /// Current buffer being read from.
    current_buf: Vec<u8>,
    /// Read position within `current_buf`.
    current_pos: usize,
    /// Reusable scratch buffer for records (or their length prefixes) that
    /// straddle a decompressed-block boundary and therefore cannot be borrowed
    /// directly out of `current_buf`. See [`PooledInputStream::next_record_borrowed`].
    scratch: RawRecord,
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
            reorder: fgumi_bam_io::ReorderBuffer::new(),
            current_buf: Vec::new(),
            current_pos: 0,
            scratch: RawRecord::new(),
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

    /// Read the next raw BAM record, borrowing its bytes from the current
    /// decompressed block when possible.
    ///
    /// This is the borrow-in-place counterpart to [`read_raw_record`]: it removes
    /// the per-record `read_exact` copy into a `RawRecord` on the common path
    /// where the record body lies wholly within the current decompressed block.
    ///
    /// - **Fast path:** when the 4-byte `block_size` prefix and the full record
    ///   body are both contained in `current_buf`, returns a slice borrowed
    ///   directly out of `current_buf` (one copy: block → caller's arena).
    /// - **Slow path (block straddle):** when the prefix *or* the body spans a
    ///   decompressed-block boundary, the bytes are gathered into a reusable
    ///   internal scratch buffer via the byte-stream [`IoRead`] impl and a slice
    ///   into that scratch is returned.
    ///
    /// The returned slice borrows `self`; it is invalidated by the next call to
    /// any method on this stream. Returns `Ok(None)` at clean EOF.
    ///
    /// A `block_size` of 0 is treated as EOF, mirroring [`read_raw_record`].
    ///
    /// # Errors
    ///
    /// Returns an error if the underlying block stream errors (I/O or
    /// decompression), the prefix/body is truncated, or `block_size` overflows
    /// `usize`.
    pub fn next_record_borrowed(&mut self) -> std::io::Result<Option<&[u8]>> {
        // --- read the 4-byte block_size prefix ---
        let block_size = if self.current_buf.len() - self.current_pos >= 4 {
            // Fast path: the prefix is fully buffered — decode it in place.
            let p = self.current_pos;
            let n = u32::from_le_bytes([
                self.current_buf[p],
                self.current_buf[p + 1],
                self.current_buf[p + 2],
                self.current_buf[p + 3],
            ]);
            self.current_pos += 4;
            usize::try_from(n)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?
        } else {
            // Slow path: the prefix straddles a block boundary (or the buffer is
            // exhausted / at EOF). `read_block_size` reads across blocks via the
            // `IoRead` impl and returns 0 at clean EOF.
            fgumi_raw_bam::read_block_size(self)?
        };
        if block_size == 0 {
            return Ok(None); // EOF (or a zero-length prefix, treated as EOF)
        }

        // --- read the record body ---
        if self.current_buf.len() - self.current_pos >= block_size {
            // Fast path: the body lies wholly within the current block — borrow it.
            let start = self.current_pos;
            self.current_pos += block_size;
            Ok(Some(&self.current_buf[start..start + block_size]))
        } else {
            // Slow path: the body straddles a block boundary — gather it into the
            // reusable scratch buffer. Take the scratch out so `read_exact` can
            // borrow `self` mutably, then restore it (preserving its capacity).
            let mut scratch = std::mem::take(&mut self.scratch);
            let buf = scratch.as_mut_vec();
            buf.resize(block_size, 0);
            let res = self.read_exact(buf);
            self.scratch = scratch;
            res?;
            Ok(Some(self.scratch.as_ref()))
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
    ///
    /// The second field holds the most-recently-yielded owned record so that
    /// [`RecordSource::next_record_borrowed`] can lend a slice into it: this path
    /// delivers owned `RawRecord`s over a channel, so there is nothing in a
    /// shared block to borrow — the owned record is the thing we lend.
    #[allow(dead_code)] // retained for potential future use / benchmarking
    ReadAhead(RawReadAheadReader, Option<RawRecord>),
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
    ///
    /// Note: [`RecordSource::next_record_borrowed`] on the `Direct` path
    /// propagates errors directly via its `Result` rather than stashing them
    /// here, so for that path `take_error()` returns `None`. The `ReadAhead`
    /// path still surfaces background-thread errors here.
    pub fn take_error(&mut self) -> Option<std::io::Error> {
        match self {
            Self::Direct(_, err) => err.take(),
            Self::ReadAhead(r, _) => r.take_error(),
        }
    }

    /// Read the next record, borrowing its bytes in place where possible.
    ///
    /// On the `Direct` (pool) path this borrows straight out of the current
    /// decompressed block via [`PooledInputStream::next_record_borrowed`],
    /// avoiding the per-record copy into a `RawRecord`. On the `ReadAhead` path
    /// the background thread already owns the record, so we store it and lend a
    /// slice into it.
    ///
    /// The returned slice borrows `self` and is invalidated by the next call.
    /// Returns `Ok(None)` at clean EOF.
    ///
    /// # Errors
    ///
    /// Returns an error if the underlying `Direct` block stream errors or the
    /// input is truncated. The `ReadAhead` path defers errors to
    /// [`RecordSource::take_error`].
    pub fn next_record_borrowed(&mut self) -> std::io::Result<Option<&[u8]>> {
        match self {
            Self::Direct(reader, _error_slot) => reader.get_mut().next_record_borrowed(),
            Self::ReadAhead(r, held) => {
                *held = r.next_record();
                Ok(held.as_ref().map(RawRecord::as_ref))
            }
        }
    }
}

impl Iterator for RecordSource {
    type Item = RawRecord;

    // PERF NOTE (S3-015): the `Direct` arm below allocates a fresh `RawRecord`
    // per `next()`. This is INHERENT to owned iteration, not a missed reuse: the
    // sole hot consumer (`sort_queryname_keyed`) *stores* every yielded record
    // in its in-memory `entries` buffer until spill, so each live record needs
    // its own buffer — a scratch + `mem::take` would just hand out the scratch's
    // buffer and leave an empty one behind, reallocating on the very next read
    // (identical to `RawRecord::default()`). The principled fix is to put the
    // queryname path on the same byte arena the coordinate/template paths
    // already use (copy bytes once into a `SegmentedBuf`, store `(key, range)`
    // refs, pool the keys), which makes this owned `Iterator::next` dead. That
    // is a design + sort-bench effort tracked in
    // `docs/design/sort-queryname-arena-deferral.md`.
    fn next(&mut self) -> Option<Self::Item> {
        match self {
            Self::ReadAhead(r, _) => r.next(),
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

    // ── PooledInputStream::next_record_borrowed — block-straddle tests ───────
    //
    // These build a PooledInputStream whose decompressed blocks are a fixed
    // (often tiny) size, so records and their 4-byte length prefixes straddle
    // block boundaries. The borrow-in-place reader must reassemble straddling
    // records via its scratch buffer and produce byte-identical output to the
    // owned read_record path.

    use crossbeam_queue::ArrayQueue;
    use proptest::{prop_assert_eq, proptest};
    use rstest::rstest;
    use std::sync::Arc;
    use std::sync::atomic::AtomicBool;

    /// Frame a record body with its 4-byte little-endian `block_size` prefix.
    fn frame_record(body: &[u8]) -> Vec<u8> {
        let mut v = Vec::with_capacity(4 + body.len());
        v.extend_from_slice(&u32::try_from(body.len()).expect("body fits u32").to_le_bytes());
        v.extend_from_slice(body);
        v
    }

    /// Build a `PooledInputStream` whose decompressed blocks are exactly
    /// `block_len` bytes each (the last may be shorter). Small `block_len`
    /// values force prefixes and bodies to straddle block boundaries.
    fn pooled_stream_from(records: &[Vec<u8>], block_len: usize) -> PooledInputStream {
        assert!(block_len >= 1, "block_len must be positive");
        let mut stream = Vec::new();
        for body in records {
            stream.extend_from_slice(&frame_record(body));
        }
        let chunks: Vec<Vec<u8>> = if stream.is_empty() {
            Vec::new()
        } else {
            stream.chunks(block_len).map(<[u8]>::to_vec).collect()
        };
        let queue = Arc::new(ArrayQueue::new(chunks.len().max(1)));
        for (serial, chunk) in chunks.into_iter().enumerate() {
            let serial = u64::try_from(serial).expect("serial fits u64");
            queue.push((serial, chunk)).expect("queue has capacity for all chunks");
        }
        PooledInputStream::new(
            queue,
            Arc::new(AtomicBool::new(true)), // decompressed_input_done
            Arc::new(AtomicBool::new(false)), // input_read_error
            Arc::new(AtomicBool::new(false)), // decompression_error
        )
    }

    /// Drain all records via the borrowing API, copying each borrowed slice into
    /// an owned `Vec` for comparison.
    fn collect_borrowed(stream: &mut PooledInputStream) -> Vec<Vec<u8>> {
        let mut out = Vec::new();
        while let Some(rec) = stream.next_record_borrowed().expect("borrow read should succeed") {
            out.push(rec.to_vec());
        }
        out
    }

    /// Record bodies of assorted lengths with position-dependent content, so a
    /// misaligned read produces detectably wrong bytes.
    fn sample_record_bodies() -> Vec<Vec<u8>> {
        let lens = [1usize, 2, 3, 4, 5, 8, 13, 32, 33, 100, 255, 256, 257];
        lens.iter()
            .enumerate()
            .map(|(k, &len)| {
                let base = u8::try_from(k % 256).expect("k%256 fits u8");
                (0..len)
                    .map(|j| base.wrapping_add(u8::try_from(j % 256).expect("j%256 fits u8")))
                    .collect()
            })
            .collect()
    }

    #[rstest]
    fn test_next_record_borrowed_matches_input_across_block_sizes(
        // Tiny sizes force straddles; large sizes (> whole stream) take the
        // all-in-one-block fast path.
        #[values(1usize, 2, 3, 4, 5, 6, 7, 8, 16, 64, 1024, 65_535)] block_len: usize,
    ) {
        let records = sample_record_bodies();
        let mut stream = pooled_stream_from(&records, block_len);
        let got = collect_borrowed(&mut stream);
        assert_eq!(got, records, "record mismatch at block_len={block_len}");
    }

    #[test]
    fn test_next_record_borrowed_empty_stream() {
        let mut stream = pooled_stream_from(&[], 4);
        assert!(
            stream.next_record_borrowed().expect("empty stream read should succeed").is_none(),
            "empty stream should yield no records"
        );
    }

    #[test]
    fn test_next_record_borrowed_prefix_straddle() {
        // block_len = 2 guarantees every record's 4-byte prefix spans at least
        // two blocks, exercising the slow prefix path (fgumi_raw_bam::read_block_size).
        let records = vec![vec![0xAB; 10], vec![0xCD; 7], vec![0xEF; 1]];
        let mut stream = pooled_stream_from(&records, 2);
        assert_eq!(collect_borrowed(&mut stream), records);
    }

    #[test]
    fn test_next_record_borrowed_body_exact_boundary_and_straddle() {
        // With block_len = 10 the first framed record ([prefix(4)|body(6)] = 10
        // bytes) ends exactly on a block boundary; the next records straddle.
        let records = vec![vec![1u8; 6], vec![2u8; 9], vec![3u8; 3]];
        let mut stream = pooled_stream_from(&records, 10);
        assert_eq!(collect_borrowed(&mut stream), records);
    }

    #[rstest]
    fn test_next_record_borrowed_parity_with_read_record(
        #[values(1usize, 3, 7, 64)] block_len: usize,
    ) {
        // The borrow-in-place reader must produce byte-identical records to the
        // owned read_raw_record path over the same chunked stream.
        let records = sample_record_bodies();
        let mut borrowed_stream = pooled_stream_from(&records, block_len);
        let borrowed = collect_borrowed(&mut borrowed_stream);

        let owned_stream = pooled_stream_from(&records, block_len);
        let mut reader = fgumi_raw_bam::RawBamReader::new(owned_stream);
        let mut owned = Vec::new();
        let mut rec = fgumi_raw_bam::RawRecord::new();
        loop {
            let n = reader.read_record(&mut rec).expect("read_record should succeed");
            if n == 0 {
                break;
            }
            owned.push(rec.as_ref().to_vec());
        }

        assert_eq!(borrowed, owned, "borrowed vs owned mismatch at block_len={block_len}");
        assert_eq!(borrowed, records, "borrowed records must match input");
    }

    #[test]
    fn test_next_record_borrowed_truncated_body_errors() {
        // A framed record claiming a 20-byte body but only 5 bytes present must
        // surface an error (not silently return a short or wrong record).
        let mut stream = Vec::new();
        stream.extend_from_slice(&20u32.to_le_bytes());
        stream.extend_from_slice(&[7u8; 5]);
        let queue = Arc::new(ArrayQueue::new(1));
        queue.push((0u64, stream)).expect("push");
        let mut pooled = PooledInputStream::new(
            queue,
            Arc::new(AtomicBool::new(true)),
            Arc::new(AtomicBool::new(false)),
            Arc::new(AtomicBool::new(false)),
        );
        let err = pooled.next_record_borrowed().expect_err("truncated body should error");
        assert_eq!(err.kind(), std::io::ErrorKind::UnexpectedEof);
    }

    proptest! {
        /// Property: over randomized record bodies and decompressed-block sizes,
        /// the borrow-in-place reader yields records byte-identical to the input
        /// (and, transitively, to the owned `read_record` path). This widens the
        /// boundary coverage of the fixed straddle examples — any combination of
        /// record length and block length where a prefix or body crosses a block
        /// edge must still reassemble correctly via the scratch buffer.
        #[test]
        fn prop_next_record_borrowed_matches_input(
            // Up to 24 records, each 1..=300 bytes of arbitrary content.
            records in proptest::collection::vec(
                proptest::collection::vec(proptest::num::u8::ANY, 1..=300),
                0..=24,
            ),
            // Block sizes from 1 (every prefix straddles) up past the largest
            // record (whole-record fast path).
            block_len in 1usize..=512,
        ) {
            let mut stream = pooled_stream_from(&records, block_len);
            let got = collect_borrowed(&mut stream);
            prop_assert_eq!(got, records);
        }
    }
}

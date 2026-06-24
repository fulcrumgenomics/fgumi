//! Multithreaded BGZF writer with coarse progress counters.
//!
//! This module provides [`MultithreadedWriter`], a parallel BGZF compression writer
//! that maintains block ordering and exposes a few coarse progress counters.
//!
//! Some accessors and constructors are only exercised by tests; the
//! `dead_code` allow at module scope keeps them available for the test harness.

#![allow(dead_code)]
//!
//! # Progress Counters (not an exact block→offset index)
//!
//! The writer exposes a few coarse counters describing how far compression and
//! output have progressed:
//!
//! - [`current_block_number`](MultithreadedWriter::current_block_number): the block index
//!   currently being filled (incremented when a block is *enqueued* for compression)
//! - [`buffer_offset`](MultithreadedWriter::buffer_offset): bytes in the current block's staging buffer
//! - [`position`](MultithreadedWriter::position): total compressed bytes flushed by the writer thread
//! - [`blocks_written`](MultithreadedWriter::blocks_written): count of blocks flushed by the writer thread
//!
//! **These counters are NOT an exact block→compressed-offset map and must not be used
//! to build BAI/CSI virtual offsets.** `current_block_number` advances at enqueue time
//! while `position`/`blocks_written` reflect only what the writer thread has already
//! flushed, so with multiple blocks in flight there is no way to recover the unique
//! compressed start offset of a given block from these globals. Index construction
//! that needs exact offsets must obtain them another way (e.g. a real per-block
//! position API, or a second pass that re-reads the finished BGZF stream — the latter
//! is what fgumi's BAI generation does).

use std::io::{self, Write};
use std::mem;
use std::num::NonZero;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};
use std::thread::{self, JoinHandle};

use bytes::{Bytes, BytesMut};
use crossbeam_channel::{Receiver, Sender, bounded};

use bgzf::{BgzfError, CompressionLevel, Compressor};

// Re-define constants since they're not public in bgzf crate
/// The maximum uncompressed blocksize for BGZF compression.
const BGZF_BLOCK_SIZE: usize = 65280;

/// BGZF EOF marker block.
static BGZF_EOF: &[u8] = &[
    0x1f, 0x8b, // ID1, ID2
    0x08, // CM = DEFLATE
    0x04, // FLG = FEXTRA
    0x00, 0x00, 0x00, 0x00, // MTIME = 0
    0x00, // XFL = 0
    0xff, // OS = 255 (unknown)
    0x06, 0x00, // XLEN = 6
    0x42, 0x43, // SI1, SI2
    0x02, 0x00, // SLEN = 2
    0x1b, 0x00, // BSIZE = 27
    0x03, 0x00, // CDATA
    0x00, 0x00, 0x00, 0x00, // CRC32 = 0x00000000
    0x00, 0x00, 0x00, 0x00, // ISIZE = 0
];

type BgzfResult<T> = Result<T, BgzfError>;

// ============================================================================
// Types
// ============================================================================

// Internal channel types
type FrameParts = (Vec<u8>, u32, usize); // (compressed_data, crc32, uncompressed_size)
type BufferedTx = Sender<BgzfResult<FrameParts>>;
type BufferedRx = Receiver<BgzfResult<FrameParts>>;
type DeflateTx = Sender<(Bytes, BufferedTx)>;
type DeflateRx = Receiver<(Bytes, BufferedTx)>;
type WriteTx = Sender<(BufferedRx, u64)>;
type WriteRx = Receiver<(BufferedRx, u64)>;

// ============================================================================
// State
// ============================================================================

enum State<W> {
    Running {
        writer_handle: JoinHandle<BgzfResult<W>>,
        deflater_handles: Vec<JoinHandle<()>>,
        write_tx: WriteTx,
        deflate_tx: DeflateTx,
    },
    Done,
}

// ============================================================================
// MultithreadedWriter
// ============================================================================

/// A multithreaded BGZF writer with position tracking support.
///
/// Uses a thread pool to compress blocks in parallel while maintaining
/// block ordering in the output stream. Provides position tracking APIs
/// for building BAM indexes during concurrent writes.
///
/// # Thread Safety
///
/// The writer itself is not `Send` or `Sync`.
pub struct MultithreadedWriter<W>
where
    W: Write + Send + 'static,
{
    state: State<W>,
    buf: BytesMut,
    blocksize: usize,
    current_block_number: u64,
    position: Arc<AtomicU64>,
    blocks_written: Arc<AtomicU64>,
}

impl<W> MultithreadedWriter<W>
where
    W: Write + Send + 'static,
{
    /// Creates a new multithreaded writer with a single worker thread.
    ///
    /// For better parallelism, use [`with_worker_count`](Self::with_worker_count)
    /// or [`Builder`].
    pub fn new(inner: W, compression_level: CompressionLevel) -> Self {
        Self::with_worker_count(NonZero::<usize>::MIN, inner, compression_level)
    }

    /// Creates a new multithreaded writer with the specified number of worker threads.
    ///
    /// # Arguments
    ///
    /// * `worker_count` - Number of compression worker threads
    /// * `inner` - The underlying writer
    /// * `compression_level` - Compression level (1-12)
    pub fn with_worker_count(
        worker_count: NonZero<usize>,
        inner: W,
        compression_level: CompressionLevel,
    ) -> Self {
        Self::with_capacity(worker_count, inner, compression_level, BGZF_BLOCK_SIZE)
    }

    /// Creates a new multithreaded writer with custom block size.
    fn with_capacity(
        worker_count: NonZero<usize>,
        inner: W,
        compression_level: CompressionLevel,
        blocksize: usize,
    ) -> Self {
        let position = Arc::new(AtomicU64::new(0));
        let blocks_written = Arc::new(AtomicU64::new(0));

        // Create channels
        let (deflate_tx, deflate_rx) = bounded(worker_count.get());
        let (write_tx, write_rx) = bounded(worker_count.get());

        // Spawn deflater workers
        let deflater_handles = spawn_deflaters(compression_level, worker_count, deflate_rx);

        // Spawn writer thread
        let writer_handle =
            spawn_writer(inner, write_rx, Arc::clone(&position), Arc::clone(&blocks_written));

        Self {
            state: State::Running { writer_handle, deflater_handles, write_tx, deflate_tx },
            buf: BytesMut::with_capacity(blocksize),
            blocksize,
            current_block_number: 0,
            position,
            blocks_written,
        }
    }

    // === Progress Counters (coarse; not an exact block→offset index) ===

    /// Returns the index of the block currently being filled.
    ///
    /// This is incremented each time a block is *enqueued* for compression, so it
    /// runs ahead of the writer thread. It is a block *identifier*, not a position:
    /// it cannot be combined with [`position`](Self::position) to recover a block's
    /// exact compressed offset (see the module-level "Progress Counters" note).
    #[inline]
    #[must_use]
    pub fn current_block_number(&self) -> u64 {
        self.current_block_number
    }

    /// Returns the number of bytes in the staging buffer.
    ///
    /// This is the uncompressed offset within the current block.
    /// Resets to 0 after each flush.
    #[inline]
    #[must_use]
    pub fn buffer_offset(&self) -> usize {
        self.buf.len()
    }

    /// Returns the total number of compressed bytes written to the output.
    #[inline]
    #[must_use]
    pub fn position(&self) -> u64 {
        self.position.load(Ordering::Acquire)
    }

    /// Returns the count of blocks fully written to the output.
    ///
    /// This may lag behind [`current_block_number`](Self::current_block_number)
    /// during active compression.
    #[inline]
    #[must_use]
    pub fn blocks_written(&self) -> u64 {
        self.blocks_written.load(Ordering::Acquire)
    }

    // === Lifecycle ===

    /// Finishes the stream and returns the underlying writer.
    ///
    /// This method:
    /// 1. Flushes any remaining buffered data
    /// 2. Shuts down worker threads
    /// 3. Writes the BGZF EOF marker
    ///
    /// # Errors
    ///
    /// Returns an error if compression or I/O fails.
    pub fn finish(&mut self) -> BgzfResult<W> {
        // Flush remaining buffer
        if !self.buf.is_empty() {
            self.send()?;
        }

        // Swap state to Done
        let state = mem::replace(&mut self.state, State::Done);

        match state {
            State::Running { writer_handle, mut deflater_handles, write_tx, deflate_tx } => {
                // Drop channels to signal shutdown
                drop(deflate_tx);

                // Join deflater threads
                for handle in deflater_handles.drain(..) {
                    handle
                        .join()
                        .map_err(|_| BgzfError::Io(io::Error::other("Deflater thread panicked")))?;
                }

                // Drop write channel to signal writer to finish
                drop(write_tx);

                // Join writer thread (it appends EOF)
                writer_handle
                    .join()
                    .map_err(|_| BgzfError::Io(io::Error::other("Writer thread panicked")))?
            }
            State::Done => Err(BgzfError::Io(io::Error::other("finish() called twice"))),
        }
    }

    // === Internal ===

    #[inline]
    fn remaining(&self) -> usize {
        self.blocksize.saturating_sub(self.buf.len())
    }

    #[inline]
    fn has_remaining(&self) -> bool {
        self.remaining() > 0
    }

    fn send(&mut self) -> BgzfResult<()> {
        if self.buf.is_empty() {
            return Ok(());
        }

        let State::Running { write_tx, deflate_tx, .. } = &self.state else {
            return Err(BgzfError::Io(io::Error::other("Writer already finished")));
        };

        // Take buffer contents
        let data = self.buf.split().freeze();

        // Create per-block result channel
        let (buffered_tx, buffered_rx) = bounded(1);

        // Send to deflater
        deflate_tx
            .send((data, buffered_tx))
            .map_err(|_| BgzfError::Io(io::Error::other("Deflate channel closed")))?;

        // Send to writer (maintains ordering)
        write_tx
            .send((buffered_rx, self.current_block_number))
            .map_err(|_| BgzfError::Io(io::Error::other("Write channel closed")))?;

        self.current_block_number += 1;
        Ok(())
    }
}

impl<W> Write for MultithreadedWriter<W>
where
    W: Write + Send + 'static,
{
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let amt = self.remaining().min(buf.len());
        self.buf.extend_from_slice(&buf[..amt]);

        if !self.has_remaining() {
            self.send().map_err(|e| io::Error::other(e.to_string()))?;
        }

        Ok(amt)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.send().map_err(|e| io::Error::other(e.to_string()))
    }
}

impl<W> Drop for MultithreadedWriter<W>
where
    W: Write + Send + 'static,
{
    fn drop(&mut self) {
        if matches!(self.state, State::Running { .. }) {
            let _ = self.finish();
        }
    }
}

// ============================================================================
// Builder
// ============================================================================

/// Builder for [`MultithreadedWriter`] with configurable options.
#[derive(Debug, Clone)]
pub struct Builder {
    compression_level: CompressionLevel,
    worker_count: NonZero<usize>,
    blocksize: usize,
}

impl Default for Builder {
    fn default() -> Self {
        Self {
            compression_level: CompressionLevel::new(6)
                .expect("compression level 6 is always valid"),
            worker_count: NonZero::<usize>::MIN,
            blocksize: BGZF_BLOCK_SIZE,
        }
    }
}

impl Builder {
    /// Sets the compression level (1-12).
    #[must_use]
    pub fn set_compression_level(mut self, level: CompressionLevel) -> Self {
        self.compression_level = level;
        self
    }

    /// Sets the number of compression worker threads.
    #[must_use]
    pub fn set_worker_count(mut self, count: NonZero<usize>) -> Self {
        self.worker_count = count;
        self
    }

    /// Sets the maximum uncompressed block size.
    ///
    /// Must be <= `BGZF_BLOCK_SIZE` (65280).
    #[must_use]
    #[allow(dead_code)]
    pub fn set_blocksize(mut self, size: usize) -> Self {
        self.blocksize = size.min(BGZF_BLOCK_SIZE);
        self
    }

    /// Builds the [`MultithreadedWriter`] from the given writer.
    pub fn build_from_writer<W>(self, writer: W) -> MultithreadedWriter<W>
    where
        W: Write + Send + 'static,
    {
        MultithreadedWriter::with_capacity(
            self.worker_count,
            writer,
            self.compression_level,
            self.blocksize,
        )
    }
}

// ============================================================================
// Worker Functions
// ============================================================================

#[allow(clippy::needless_pass_by_value)]
fn spawn_deflaters(
    compression_level: CompressionLevel,
    worker_count: NonZero<usize>,
    deflate_rx: DeflateRx,
) -> Vec<JoinHandle<()>> {
    (0..worker_count.get())
        .map(|_| {
            let rx = deflate_rx.clone();
            let level = compression_level;
            thread::spawn(move || {
                let mut compressor = Compressor::new(level);
                let mut compress_buf = Vec::new();

                while let Ok((data, result_tx)) = rx.recv() {
                    let result = compress_block(&mut compressor, &data, &mut compress_buf);
                    let _ = result_tx.send(result);
                }
            })
        })
        .collect()
}

fn compress_block(
    compressor: &mut Compressor,
    input: &[u8],
    buffer: &mut Vec<u8>,
) -> BgzfResult<FrameParts> {
    buffer.clear();

    // Use the existing Compressor::compress method which handles header and footer
    compressor.compress(input, buffer)?;

    // Extract CRC32 from the footer (last 8 bytes: 4 bytes CRC32 + 4 bytes uncompressed size)
    let crc32 = u32::from_le_bytes([
        buffer[buffer.len() - 8],
        buffer[buffer.len() - 7],
        buffer[buffer.len() - 6],
        buffer[buffer.len() - 5],
    ]);

    // Move the compressed block out instead of cloning it (a full per-block
    // memcpy). Reinstall a fresh buffer pre-sized to the same capacity so the
    // next block compresses into it without regrowing from empty — `mem::take`
    // would leave a zero-capacity buffer and force that regrowth.
    let cap = buffer.capacity();
    let frame = std::mem::replace(buffer, Vec::with_capacity(cap));
    Ok((frame, crc32, input.len()))
}

fn spawn_writer<W>(
    mut writer: W,
    write_rx: WriteRx,
    position: Arc<AtomicU64>,
    blocks_written: Arc<AtomicU64>,
) -> JoinHandle<BgzfResult<W>>
where
    W: Write + Send + 'static,
{
    thread::spawn(move || {
        while let Ok((buffered_rx, _block_number)) = write_rx.recv() {
            // Wait for compression result
            let (compressed_data, _crc32, _uncompressed_size) = buffered_rx
                .recv()
                .map_err(|_| BgzfError::Io(io::Error::other("Compression channel closed")))??;

            let compressed_size = compressed_data.len();

            // Write compressed block
            writer.write_all(&compressed_data)?;

            // Update position
            position.fetch_add(compressed_size as u64, Ordering::Release);
            blocks_written.fetch_add(1, Ordering::Release);
        }

        // Write EOF marker
        writer.write_all(BGZF_EOF)?;
        position.fetch_add(BGZF_EOF.len() as u64, Ordering::Release);

        Ok(writer)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use bgzf::Reader;
    use std::io::Read;

    #[test]
    fn test_new_and_finish() {
        let mut writer = MultithreadedWriter::new(
            Vec::new(),
            CompressionLevel::new(6).expect("valid compression level 6"),
        );
        let result = writer.finish();
        assert!(result.is_ok());
    }

    #[test]
    fn test_roundtrip() {
        let mut writer = MultithreadedWriter::new(
            Vec::new(),
            CompressionLevel::new(6).expect("valid compression level 6"),
        );
        writer.write_all(b"hello world").expect("write_all should succeed");
        let data = writer.finish().expect("finish should succeed");

        let mut reader = Reader::new(&data[..]);
        let mut buf = Vec::new();
        reader.read_to_end(&mut buf).expect("read_to_end should succeed");
        assert_eq!(buf, b"hello world");
    }

    #[test]
    fn test_position_tracking_initial() {
        let writer = MultithreadedWriter::new(
            Vec::new(),
            CompressionLevel::new(6).expect("valid compression level 6"),
        );
        assert_eq!(writer.position(), 0);
        assert_eq!(writer.current_block_number(), 0);
        assert_eq!(writer.blocks_written(), 0);
        assert_eq!(writer.buffer_offset(), 0);
    }

    #[test]
    fn test_position_tracking_after_writes() {
        let mut writer = MultithreadedWriter::new(
            Vec::new(),
            CompressionLevel::new(6).expect("valid compression level 6"),
        );

        writer.write_all(b"hello").expect("write_all should succeed");
        assert_eq!(writer.buffer_offset(), 5);

        writer.flush().expect("flush should succeed");
        assert_eq!(writer.current_block_number(), 1);
        assert_eq!(writer.buffer_offset(), 0);

        writer.finish().expect("finish should succeed");
    }

    #[test]
    fn test_block_counters_track_flushed_blocks() {
        let mut writer = MultithreadedWriter::with_worker_count(
            NonZero::new(2).expect("non-zero value 2"),
            Vec::new(),
            CompressionLevel::new(6).expect("valid compression level 6"),
        );

        writer.write_all(b"block1").expect("write_all should succeed");
        writer.flush().expect("flush should succeed");
        writer.write_all(b"block2").expect("write_all should succeed");
        writer.flush().expect("flush should succeed");

        // Two blocks were sent for compression.
        assert_eq!(writer.current_block_number(), 2);

        writer.finish().expect("finish should succeed");

        // Both blocks reached the writer thread (position advanced past 0).
        assert_eq!(writer.blocks_written(), 2);
        assert!(writer.position() > 0);
    }

    #[test]
    fn test_multiple_workers() {
        for worker_count in [1, 2, 4] {
            let mut writer = MultithreadedWriter::with_worker_count(
                NonZero::new(worker_count).expect("worker_count is non-zero"),
                Vec::new(),
                CompressionLevel::new(6).expect("valid compression level 6"),
            );

            for i in 0..20 {
                writer.write_all(format!("block{i}").as_bytes()).expect("write_all should succeed");
                writer.flush().expect("flush should succeed");
            }

            let data = writer.finish().expect("finish should succeed");

            let mut reader = Reader::new(&data[..]);
            let mut buf = String::new();
            reader.read_to_string(&mut buf).expect("read_to_string should succeed");

            for i in 0..20 {
                assert!(buf.contains(&format!("block{i}")));
            }
        }
    }

    #[test]
    fn test_large_data() {
        let mut writer = MultithreadedWriter::with_worker_count(
            NonZero::new(4).expect("non-zero value 4"),
            Vec::new(),
            CompressionLevel::new(6).expect("valid compression level 6"),
        );

        let large_data = vec![b'A'; BGZF_BLOCK_SIZE * 5];
        writer.write_all(&large_data).expect("write_all should succeed");

        let compressed = writer.finish().expect("finish should succeed");

        let mut reader = Reader::new(&compressed[..]);
        let mut decompressed = Vec::new();
        reader.read_to_end(&mut decompressed).expect("read_to_end should succeed");

        assert_eq!(decompressed, large_data);
    }

    #[test]
    fn test_eof_marker() {
        let mut writer = MultithreadedWriter::new(
            Vec::new(),
            CompressionLevel::new(6).expect("valid compression level 6"),
        );
        writer.write_all(b"test").expect("write_all should succeed");
        let data = writer.finish().expect("finish should succeed");

        assert!(data.ends_with(BGZF_EOF));
    }

    #[test]
    fn test_builder() {
        let mut writer = Builder::default()
            .set_compression_level(CompressionLevel::new(9).expect("valid compression level 9"))
            .set_worker_count(NonZero::new(4).expect("non-zero value 4"))
            .set_blocksize(32768)
            .build_from_writer(Vec::new());

        writer.write_all(b"builder roundtrip").expect("write_all should succeed");
        let data = writer.finish().expect("finish should succeed");

        let mut reader = Reader::new(&data[..]);
        let mut buf = Vec::new();
        reader.read_to_end(&mut buf).expect("read_to_end should succeed");
        assert_eq!(buf, b"builder roundtrip");
    }
}

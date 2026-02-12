//! Multithreaded BGZF writer with position tracking.
//!
//! This module provides [`MultithreadedWriter`], a parallel BGZF compression writer
//! that maintains block ordering and supports position tracking for building indexes
//! during concurrent writes.
//!
//! # Position Tracking
//!
//! The writer provides APIs to track compressed file positions, enabling construction
//! of BAM indexes (BAI/CSI) without a second pass through the file:
//!
//! - [`current_block_number`](MultithreadedWriter::current_block_number): Block number being filled
//! - [`buffer_offset`](MultithreadedWriter::buffer_offset): Bytes in current block's staging buffer
//! - [`block_info_receiver`](MultithreadedWriter::block_info_receiver): Channel for block completion notifications

use std::io::{self, Write};
use std::mem;
use std::num::NonZero;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};
use std::thread::{self, JoinHandle};

use bytes::{Bytes, BytesMut};
use crossbeam_channel::{Receiver, Sender, bounded, unbounded};

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

/// Information about a completed BGZF block.
///
/// Sent through [`BlockInfoRx`] when a block is written to the output.
/// Use this to resolve virtual file positions for BAM index construction.
#[derive(Debug, Clone, Copy)]
pub struct BlockInfo {
    /// The block number (0-indexed, assigned when block is sent for compression).
    pub block_number: u64,
    /// The compressed file position where this block starts.
    pub compressed_start: u64,
    /// The size of the compressed block (including header and footer).
    pub compressed_size: usize,
    /// The size of the uncompressed data in this block.
    pub uncompressed_size: usize,
}

/// Receiver for [`BlockInfo`] notifications.
///
/// Clone this receiver and use `try_recv()` or `recv()` to receive
/// block completion notifications as they occur.
pub type BlockInfoRx = Receiver<BlockInfo>;

// Internal channel types
type FrameParts = (Vec<u8>, u32, usize); // (compressed_data, crc32, uncompressed_size)
type BufferedTx = Sender<BgzfResult<FrameParts>>;
type BufferedRx = Receiver<BgzfResult<FrameParts>>;
type DeflateTx = Sender<(Bytes, BufferedTx)>;
type DeflateRx = Receiver<(Bytes, BufferedTx)>;
type WriteTx = Sender<(BufferedRx, u64)>;
type WriteRx = Receiver<(BufferedRx, u64)>;
type BlockInfoTx = Sender<BlockInfo>;

// ============================================================================
// State
// ============================================================================

enum State<W> {
    Running {
        writer_handle: JoinHandle<BgzfResult<W>>,
        deflater_handles: Vec<JoinHandle<()>>,
        write_tx: WriteTx,
        deflate_tx: DeflateTx,
        block_info_rx: BlockInfoRx,
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
/// The writer itself is not `Send` or `Sync`, but the [`BlockInfoRx`] receiver
/// can be cloned and used from other threads.
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
        let (block_info_tx, block_info_rx) = unbounded();

        // Spawn deflater workers
        let deflater_handles = spawn_deflaters(compression_level, worker_count, deflate_rx);

        // Spawn writer thread
        let writer_handle = spawn_writer(
            inner,
            write_rx,
            Arc::clone(&position),
            Arc::clone(&blocks_written),
            block_info_tx,
        );

        Self {
            state: State::Running {
                writer_handle,
                deflater_handles,
                write_tx,
                deflate_tx,
                block_info_rx,
            },
            buf: BytesMut::with_capacity(blocksize),
            blocksize,
            current_block_number: 0,
            position,
            blocks_written,
        }
    }

    // === Position Tracking APIs ===

    /// Returns the receiver for block completion notifications.
    ///
    /// Returns `None` if the writer has already been finished.
    #[must_use]
    pub fn block_info_receiver(&self) -> Option<&BlockInfoRx> {
        match &self.state {
            State::Running { block_info_rx, .. } => Some(block_info_rx),
            State::Done => None,
        }
    }

    /// Returns the next block number to be assigned.
    ///
    /// This is incremented each time a block is sent for compression.
    /// Use this value when caching index entries to correlate with
    /// [`BlockInfo::block_number`] in completion notifications.
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
            State::Running {
                writer_handle,
                mut deflater_handles,
                write_tx,
                deflate_tx,
                block_info_rx: _,
            } => {
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
            compression_level: CompressionLevel::new(6).unwrap(),
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

    Ok((buffer.clone(), crc32, input.len()))
}

fn spawn_writer<W>(
    mut writer: W,
    write_rx: WriteRx,
    position: Arc<AtomicU64>,
    blocks_written: Arc<AtomicU64>,
    block_info_tx: BlockInfoTx,
) -> JoinHandle<BgzfResult<W>>
where
    W: Write + Send + 'static,
{
    thread::spawn(move || {
        while let Ok((buffered_rx, block_number)) = write_rx.recv() {
            // Wait for compression result
            let (compressed_data, _crc32, uncompressed_size) = buffered_rx
                .recv()
                .map_err(|_| BgzfError::Io(io::Error::other("Compression channel closed")))??;

            let compressed_start = position.load(Ordering::Acquire);
            let compressed_size = compressed_data.len();

            // Write compressed block
            writer.write_all(&compressed_data)?;

            // Update position
            position.fetch_add(compressed_size as u64, Ordering::Release);
            blocks_written.fetch_add(1, Ordering::Release);

            // Send block info notification
            let _ = block_info_tx.send(BlockInfo {
                block_number,
                compressed_start,
                compressed_size,
                uncompressed_size,
            });
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
        let mut writer = MultithreadedWriter::new(Vec::new(), CompressionLevel::new(6).unwrap());
        let result = writer.finish();
        assert!(result.is_ok());
    }

    #[test]
    fn test_roundtrip() {
        let mut writer = MultithreadedWriter::new(Vec::new(), CompressionLevel::new(6).unwrap());
        writer.write_all(b"hello world").unwrap();
        let data = writer.finish().unwrap();

        let mut reader = Reader::new(&data[..]);
        let mut buf = Vec::new();
        reader.read_to_end(&mut buf).unwrap();
        assert_eq!(buf, b"hello world");
    }

    #[test]
    fn test_position_tracking_initial() {
        let writer = MultithreadedWriter::new(Vec::new(), CompressionLevel::new(6).unwrap());
        assert_eq!(writer.position(), 0);
        assert_eq!(writer.current_block_number(), 0);
        assert_eq!(writer.blocks_written(), 0);
        assert_eq!(writer.buffer_offset(), 0);
    }

    #[test]
    fn test_position_tracking_after_writes() {
        let mut writer = MultithreadedWriter::new(Vec::new(), CompressionLevel::new(6).unwrap());

        writer.write_all(b"hello").unwrap();
        assert_eq!(writer.buffer_offset(), 5);

        writer.flush().unwrap();
        assert_eq!(writer.current_block_number(), 1);
        assert_eq!(writer.buffer_offset(), 0);

        writer.finish().unwrap();
    }

    #[test]
    fn test_block_info_notifications() {
        let mut writer = MultithreadedWriter::with_worker_count(
            NonZero::new(2).unwrap(),
            Vec::new(),
            CompressionLevel::new(6).unwrap(),
        );

        let rx = writer.block_info_receiver().unwrap().clone();

        writer.write_all(b"block1").unwrap();
        writer.flush().unwrap();
        writer.write_all(b"block2").unwrap();
        writer.flush().unwrap();

        writer.finish().unwrap();

        let infos: Vec<_> = rx.try_iter().collect();
        assert_eq!(infos.len(), 2);
        assert_eq!(infos[0].block_number, 0);
        assert_eq!(infos[1].block_number, 1);
        assert!(infos[1].compressed_start > 0);
    }

    #[test]
    fn test_multiple_workers() {
        for worker_count in [1, 2, 4] {
            let mut writer = MultithreadedWriter::with_worker_count(
                NonZero::new(worker_count).unwrap(),
                Vec::new(),
                CompressionLevel::new(6).unwrap(),
            );

            for i in 0..20 {
                writer.write_all(format!("block{i}").as_bytes()).unwrap();
                writer.flush().unwrap();
            }

            let data = writer.finish().unwrap();

            let mut reader = Reader::new(&data[..]);
            let mut buf = String::new();
            reader.read_to_string(&mut buf).unwrap();

            for i in 0..20 {
                assert!(buf.contains(&format!("block{i}")));
            }
        }
    }

    #[test]
    fn test_large_data() {
        let mut writer = MultithreadedWriter::with_worker_count(
            NonZero::new(4).unwrap(),
            Vec::new(),
            CompressionLevel::new(6).unwrap(),
        );

        let large_data = vec![b'A'; BGZF_BLOCK_SIZE * 5];
        writer.write_all(&large_data).unwrap();

        let compressed = writer.finish().unwrap();

        let mut reader = Reader::new(&compressed[..]);
        let mut decompressed = Vec::new();
        reader.read_to_end(&mut decompressed).unwrap();

        assert_eq!(decompressed, large_data);
    }

    #[test]
    fn test_eof_marker() {
        let mut writer = MultithreadedWriter::new(Vec::new(), CompressionLevel::new(6).unwrap());
        writer.write_all(b"test").unwrap();
        let data = writer.finish().unwrap();

        assert!(data.ends_with(BGZF_EOF));
    }

    #[test]
    fn test_builder() {
        let writer = Builder::default()
            .set_compression_level(CompressionLevel::new(9).unwrap())
            .set_worker_count(NonZero::new(4).unwrap())
            .set_blocksize(32768)
            .build_from_writer(Vec::new());

        assert!(writer.block_info_receiver().is_some());
    }
}

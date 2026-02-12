//! Parallel gzip writer using libdeflater.
//!
//! Uses the same architecture as `bgzf_writer`:
//! 1. Main thread buffers data into blocks
//! 2. Compression workers compress blocks in parallel
//! 3. I/O thread writes compressed blocks in order
//!
//! Output is valid gzip (concatenated gzip streams per RFC 1952).

use crossbeam_channel::{Receiver, Sender, bounded};
use libdeflater::{CompressionLvl, Compressor};
use std::collections::BTreeMap;
use std::io::{self, Write};
use std::sync::Arc;
use std::thread::{self, JoinHandle};

/// Default block size (64KB like BGZF).
const DEFAULT_BLOCK_SIZE: usize = 65536;

/// An uncompressed block ready for compression.
struct UncompressedBlock {
    serial: u64,
    data: Vec<u8>,
}

/// A compressed gzip block ready for writing, or a compression error.
enum CompressedBlock {
    Ok { serial: u64, data: Vec<u8> },
    Err { serial: u64, error: String },
}

/// Compress an uncompressed block into a gzip stream.
#[allow(clippy::needless_pass_by_value)]
fn compress_block(
    block: UncompressedBlock,
    compressor: &mut Compressor,
) -> io::Result<CompressedBlock> {
    let uncompressed = &block.data;

    // libdeflater gzip_compress_bound gives max compressed size
    let max_compressed = compressor.gzip_compress_bound(uncompressed.len());
    let mut compressed_data = vec![0u8; max_compressed];

    let compressed_len = compressor
        .gzip_compress(uncompressed, &mut compressed_data)
        .map_err(|e| io::Error::other(format!("Gzip compression failed: {e:?}")))?;

    compressed_data.truncate(compressed_len);

    Ok(CompressedBlock::Ok { serial: block.serial, data: compressed_data })
}

/// Configuration for parallel gzip writer.
#[derive(Debug, Clone)]
pub struct ParallelGzipConfig {
    /// Number of compression threads.
    pub compression_threads: usize,
    /// Compression level (1-12, default 6).
    pub compression_level: u8,
    /// Queue size for pending blocks (default: 2 * threads).
    pub queue_size: Option<usize>,
    /// Block size in bytes (default: 64KB).
    pub block_size: usize,
}

impl Default for ParallelGzipConfig {
    fn default() -> Self {
        Self {
            compression_threads: 4,
            compression_level: 6,
            queue_size: None,
            block_size: DEFAULT_BLOCK_SIZE,
        }
    }
}

impl ParallelGzipConfig {
    /// Create config with specified thread count.
    #[must_use]
    pub fn with_threads(threads: usize) -> Self {
        Self { compression_threads: threads.max(1), ..Default::default() }
    }

    fn effective_queue_size(&self) -> usize {
        self.queue_size.unwrap_or(self.compression_threads * 2)
    }
}

/// A parallel gzip writer that compresses blocks using multiple threads.
pub struct ParallelGzipWriter {
    /// Current block buffer being filled.
    block_buffer: Vec<u8>,
    /// Block size threshold.
    block_size: usize,
    /// Next serial number to assign.
    next_serial: u64,
    /// Channel to send uncompressed blocks to compression workers.
    compress_tx: Sender<UncompressedBlock>,
    /// Handles for compression worker threads.
    compression_handles: Vec<JoinHandle<()>>,
    /// Handle for the I/O writer thread.
    io_handle: Option<JoinHandle<io::Result<()>>>,
}

impl ParallelGzipWriter {
    /// Create a new parallel gzip writer.
    ///
    /// # Errors
    ///
    /// Returns an error if the compression level is invalid.
    pub fn new<W>(writer: W, config: &ParallelGzipConfig) -> io::Result<Self>
    where
        W: Write + Send + 'static,
    {
        let queue_size = config.effective_queue_size();
        let compression_level = CompressionLvl::new(i32::from(config.compression_level))
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, format!("{e:?}")))?;

        // Create channels
        let (compress_tx, compress_rx) = bounded::<UncompressedBlock>(queue_size);
        let (output_tx, output_rx) = bounded::<CompressedBlock>(queue_size);

        let compress_rx = Arc::new(compress_rx);
        let output_tx = Arc::new(output_tx);

        // Spawn compression workers
        let mut compression_handles = Vec::with_capacity(config.compression_threads);
        for _ in 0..config.compression_threads {
            let rx = Arc::clone(&compress_rx);
            let tx = Arc::clone(&output_tx);
            let level = compression_level;

            let handle = thread::spawn(move || {
                let mut compressor = Compressor::new(level);
                while let Ok(block) = rx.recv() {
                    let serial = block.serial;
                    let result = match compress_block(block, &mut compressor) {
                        Ok(compressed) => compressed,
                        Err(e) => CompressedBlock::Err { serial, error: e.to_string() },
                    };
                    if tx.send(result).is_err() {
                        break;
                    }
                }
            });
            compression_handles.push(handle);
        }

        drop(compress_rx);
        drop(output_tx);

        // Spawn I/O writer thread
        let io_handle =
            thread::spawn(move || -> io::Result<()> { Self::io_writer_loop(writer, output_rx) });

        Ok(Self {
            block_buffer: Vec::with_capacity(config.block_size),
            block_size: config.block_size,
            next_serial: 0,
            compress_tx,
            compression_handles,
            io_handle: Some(io_handle),
        })
    }

    /// The I/O writer thread main loop.
    #[allow(clippy::needless_pass_by_value)]
    fn io_writer_loop<W: Write>(
        mut writer: W,
        output_rx: Receiver<CompressedBlock>,
    ) -> io::Result<()> {
        let mut next_expected: u64 = 0;
        let mut pending: BTreeMap<u64, CompressedBlock> = BTreeMap::new();

        while let Ok(block) = output_rx.recv() {
            let serial = match &block {
                CompressedBlock::Ok { serial, .. } | CompressedBlock::Err { serial, .. } => *serial,
            };
            pending.insert(serial, block);

            // Write all consecutive blocks
            while let Some(block) = pending.remove(&next_expected) {
                match block {
                    CompressedBlock::Ok { data, .. } => writer.write_all(&data)?,
                    CompressedBlock::Err { error, .. } => {
                        return Err(io::Error::other(format!(
                            "compression failed for block {next_expected}: {error}"
                        )));
                    }
                }
                next_expected += 1;
            }
        }

        // Write any remaining pending blocks
        for (_, block) in pending {
            match block {
                CompressedBlock::Ok { data, .. } => writer.write_all(&data)?,
                CompressedBlock::Err { serial, error } => {
                    return Err(io::Error::other(format!(
                        "compression failed for block {serial}: {error}"
                    )));
                }
            }
        }

        writer.flush()?;
        Ok(())
    }

    /// Dispatch the current block buffer for compression.
    fn dispatch_block(&mut self) -> io::Result<()> {
        if self.block_buffer.is_empty() {
            return Ok(());
        }

        let block = UncompressedBlock {
            serial: self.next_serial,
            data: std::mem::replace(&mut self.block_buffer, Vec::with_capacity(self.block_size)),
        };

        self.next_serial += 1;

        self.compress_tx
            .send(block)
            .map_err(|_| io::Error::new(io::ErrorKind::BrokenPipe, "Compression channel closed"))
    }

    /// Finish writing and close all threads.
    ///
    /// # Errors
    ///
    /// Returns an error if flushing or joining worker threads fails.
    pub fn finish(mut self) -> io::Result<()> {
        // Dispatch any remaining data
        self.dispatch_block()?;

        // Close the compression channel
        drop(self.compress_tx);

        // Wait for compression workers (they return () so just check for panic)
        for handle in self.compression_handles {
            handle.join().map_err(|_| io::Error::other("Compression worker thread panicked"))?;
        }

        // Wait for I/O thread
        if let Some(handle) = self.io_handle.take() {
            handle.join().map_err(|_| io::Error::other("I/O thread panicked"))??;
        }

        Ok(())
    }
}

impl Write for ParallelGzipWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.block_buffer.extend_from_slice(buf);

        // Dispatch if buffer is full
        if self.block_buffer.len() >= self.block_size {
            self.dispatch_block()?;
        }

        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        // Dispatch current block if it has data
        self.dispatch_block()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::read::MultiGzDecoder;
    use std::fs::File;
    use std::io::Read;
    use tempfile::NamedTempFile;

    #[test]
    fn test_basic_compression() -> io::Result<()> {
        let temp = NamedTempFile::new()?;
        let path = temp.path().to_path_buf();

        {
            let file = File::create(&path)?;
            let config = ParallelGzipConfig::with_threads(2);
            let mut writer = ParallelGzipWriter::new(file, &config)?;

            writer.write_all(b"Hello, World!")?;
            writer.finish()?;
        }

        // Decompress and verify (use MultiGzDecoder for concatenated streams)
        let file = File::open(&path)?;
        let mut decoder = MultiGzDecoder::new(file);
        let mut decompressed = String::new();
        decoder.read_to_string(&mut decompressed)?;

        assert_eq!(decompressed, "Hello, World!");
        Ok(())
    }

    #[test]
    fn test_multi_block_compression() -> io::Result<()> {
        let temp = NamedTempFile::new()?;
        let path = temp.path().to_path_buf();
        let test_data = "ACGT".repeat(100_000); // ~400KB, multiple blocks

        {
            let file = File::create(&path)?;
            let config = ParallelGzipConfig {
                compression_threads: 4,
                block_size: 16384, // 16KB blocks for more parallelism
                ..Default::default()
            };
            let mut writer = ParallelGzipWriter::new(file, &config)?;

            writer.write_all(test_data.as_bytes())?;
            writer.finish()?;
        }

        // Decompress and verify (use MultiGzDecoder for concatenated streams)
        let file = File::open(&path)?;
        let mut decoder = MultiGzDecoder::new(file);
        let mut decompressed = String::new();
        decoder.read_to_string(&mut decompressed)?;

        assert_eq!(decompressed, test_data);
        Ok(())
    }

    #[test]
    fn test_single_thread() -> io::Result<()> {
        let temp = NamedTempFile::new()?;
        let path = temp.path().to_path_buf();

        {
            let file = File::create(&path)?;
            let config = ParallelGzipConfig::with_threads(1);
            let mut writer = ParallelGzipWriter::new(file, &config)?;

            writer.write_all(b"Single thread test")?;
            writer.finish()?;
        }

        // Decompress and verify (use MultiGzDecoder for concatenated streams)
        let file = File::open(&path)?;
        let mut decoder = MultiGzDecoder::new(file);
        let mut decompressed = String::new();
        decoder.read_to_string(&mut decompressed)?;

        assert_eq!(decompressed, "Single thread test");
        Ok(())
    }
}

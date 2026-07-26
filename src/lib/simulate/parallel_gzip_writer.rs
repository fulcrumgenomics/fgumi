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
    /// Compression level (0-12, default 6; 0 = uncompressed).
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
    ///
    /// `Option` so shutdown can close it — dropping the sender is what tells the
    /// workers to finish — while leaving the struct intact for [`Drop`] to run
    /// over a second time without repeating the work.
    compress_tx: Option<Sender<UncompressedBlock>>,
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
            compress_tx: Some(compress_tx),
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
            .as_ref()
            .ok_or_else(|| io::Error::new(io::ErrorKind::BrokenPipe, "Writer already finished"))?
            .send(block)
            .map_err(|_| io::Error::new(io::ErrorKind::BrokenPipe, "Compression channel closed"))
    }

    /// Finish writing and close all threads.
    ///
    /// # Errors
    ///
    /// Returns an error if flushing or joining worker threads fails. Every thread
    /// is joined regardless — see `shutdown`, which this delegates to.
    pub fn finish(mut self) -> io::Result<()> {
        self.shutdown()
    }

    /// Close the channel and join every thread, returning the first error.
    ///
    /// Joining is unconditional. It used to happen behind three `?`s — a failed
    /// final dispatch returned before any join, and a worker that panicked
    /// returned before the remaining workers and the I/O thread were joined — so
    /// the error surfaced while detached threads were still compressing and still
    /// writing to the output file. The caller would see the failure, and the file
    /// would keep growing behind it.
    ///
    /// The first error a join discovers wins, so the cause is reported rather than
    /// whatever the unwinding happened to hit last; later errors are dropped, not
    /// hidden behind an earlier `Ok`. The final dispatch's own error ranks below
    /// all of them because it is usually a symptom: an I/O thread that fails
    /// closes `output_rx`, the workers break out of their loop, and the dispatch
    /// then reports a generic "Compression channel closed" for a write that
    /// actually failed on a full disk.
    ///
    /// Idempotent: the sender and the handles are taken, so the [`Drop`] that runs
    /// after [`Self::finish`] finds nothing left to do.
    fn shutdown(&mut self) -> io::Result<()> {
        // Not `?`: a failed dispatch must not skip the joins below.
        let dispatch_result = self.dispatch_block();

        // Dropping the sender is what lets the workers see the channel close and
        // exit, so it has to happen before the joins or they would block forever.
        drop(self.compress_tx.take());

        let mut result: io::Result<()> = Ok(());
        for handle in self.compression_handles.drain(..) {
            if handle.join().is_err() && result.is_ok() {
                result = Err(io::Error::other("Compression worker thread panicked"));
            }
        }

        if let Some(handle) = self.io_handle.take() {
            let io_result = match handle.join() {
                Ok(io_result) => io_result,
                Err(_) => Err(io::Error::other("I/O thread panicked")),
            };
            if let Err(e) = io_result
                && result.is_ok()
            {
                result = Err(e);
            }
        }

        // Only if every thread shut down cleanly is the dispatch failure the cause.
        result.and(dispatch_result)
    }
}

impl Drop for ParallelGzipWriter {
    /// Join the worker threads even when [`Self::finish`] was never reached.
    ///
    /// Without this, any `?` between construction and `finish` — a failed
    /// `write_record`, or the first of two writers failing to finish — detached
    /// threads that were still compressing and writing. The caller unwound while
    /// its output file was still being written to, and whether those writes landed
    /// depended on how quickly the process exited.
    ///
    /// Errors cannot propagate from a drop, so they are logged rather than
    /// swallowed silently. A caller that wants them calls [`Self::finish`].
    fn drop(&mut self) {
        if let Err(e) = self.shutdown() {
            log::warn!("Error shutting down the parallel gzip writer: {e}");
        }
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

    /// A sink that fails every write, to drive the I/O thread's error path.
    struct FailingSink;

    impl Write for FailingSink {
        fn write(&mut self, _buf: &[u8]) -> io::Result<usize> {
            Err(io::Error::other("sink is broken"))
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    /// A failing I/O thread must surface its error from `finish`, not hang.
    ///
    /// The joins used to sit behind `?`, so an error on one of them returned with
    /// the remaining threads unjoined. This drives that path: the payload is small
    /// enough that nothing is dispatched until shutdown, so the final dispatch
    /// succeeds and the I/O thread's failure is the first error — which is the one
    /// the caller has to see.
    #[test]
    fn test_finish_reports_an_io_thread_failure() {
        let config = ParallelGzipConfig::with_threads(2);
        let mut writer =
            ParallelGzipWriter::new(FailingSink, &config).expect("writer construction");
        writer.write_all(b"ACGT").expect("buffered, not yet dispatched");

        let err = writer.finish().expect_err("a failing sink must surface as an error");
        assert!(
            err.to_string().contains("sink is broken"),
            "the I/O thread's own error must be reported, got: {err}"
        );
    }

    /// The I/O thread's error must outrank the channel closure it causes.
    ///
    /// With enough payload to dispatch mid-write, the failing sink kills the I/O
    /// thread first: `output_rx` closes, the compression workers break out of
    /// their loop, and the final dispatch in `shutdown` then fails with a generic
    /// "Compression channel closed". Seeding the shutdown result from that
    /// dispatch let the symptom win over the sink's real error, which is the one
    /// the caller has to see. Writes are expected to fail here for the same
    /// reason — the assertion is on what `finish` reports.
    #[test]
    fn test_finish_prefers_the_io_error_over_the_closed_channel() {
        let config = ParallelGzipConfig::with_threads(2);
        let mut writer =
            ParallelGzipWriter::new(FailingSink, &config).expect("writer construction");

        let block = vec![b'A'; DEFAULT_BLOCK_SIZE];
        for _ in 0..64 {
            if writer.write_all(&block).is_err() {
                break;
            }
        }
        // Leave a partial block buffered so shutdown has something left to dispatch
        // down the by-now-closed channel — that failed dispatch is the symptom.
        let _ = writer.write(b"ACGT");

        let err = writer.finish().expect_err("a failing sink must surface as an error");
        assert!(
            err.to_string().contains("sink is broken"),
            "the I/O thread's own error must outrank the channel-closed symptom, got: {err}"
        );
    }

    /// Dropping a writer whose sink fails must log, not panic.
    ///
    /// This is the path that runs when a caller is already unwinding: a write
    /// failed, the `?` propagated, and the writer is dropped on the way out. A
    /// panic from `Drop` during an unwind aborts the process, so the shutdown error
    /// has to be logged and swallowed here — the caller that wants it calls
    /// `finish`. The test passing without an abort is the assertion.
    #[test]
    fn test_drop_swallows_a_shutdown_error() {
        let config = ParallelGzipConfig::with_threads(2);
        let mut writer =
            ParallelGzipWriter::new(FailingSink, &config).expect("writer construction");
        writer.write_all(b"ACGT").expect("buffered, not yet dispatched");
        drop(writer);
    }

    /// Dropping a writer without calling `finish` must still complete the output.
    ///
    /// Any `?` between construction and `finish` drops the writer — a failed
    /// `write_record`, or the first of two writers failing to finish. Without a
    /// `Drop` that joins, the compression and I/O threads were detached and went
    /// on writing to the file while the caller unwound, so what ended up on disk
    /// depended on how quickly the process exited. Reading the file straight after
    /// the drop is the assertion: it can only be complete if the drop waited.
    #[test]
    fn test_drop_without_finish_completes_the_output() -> io::Result<()> {
        let temp = NamedTempFile::new()?;
        let path = temp.path().to_path_buf();
        let payload: Vec<u8> = (0..200_000u32).map(|i| b"ACGT"[(i % 4) as usize]).collect();

        {
            let config = ParallelGzipConfig::with_threads(4);
            let mut writer = ParallelGzipWriter::new(File::create(&path)?, &config)?;
            writer.write_all(&payload)?;
            // Deliberately no `finish()`: this is the path a `?` takes.
        }

        let mut decoded = Vec::new();
        MultiGzDecoder::new(File::open(&path)?).read_to_end(&mut decoded)?;
        assert_eq!(decoded, payload, "the dropped writer left an incomplete file");
        Ok(())
    }

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

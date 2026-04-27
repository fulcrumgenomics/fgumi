//! FASTQ file Source: multi-file, format-detecting source that emits
//! [`PerStreamChunk`] items.
//!
//! The source spawns one dedicated reader thread per input stream. Each thread
//! reads its file independently and pushes [`PerStreamChunk`] items — tagged
//! with a per-stream `batch_num` and a globally-unique pipeline sequence number —
//! into the shared output queue.
//!
//! # Formats
//!
//! - **BGZF**: raw compressed blocks are emitted without decompression
//!   (downstream stages parallel-decompress). `offsets` is `None`.
//! - **Gzip / Plain**: the stream is decompressed (if needed) into
//!   record-aligned byte chunks. Each chunk carries per-record byte offsets
//!   via `Some(offsets)`.

use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;

use crate::fastq_parse::read_n_fastq_records;
use crate::runall::engine::cancel::CancelToken;
use crate::runall::engine::fastq_types::PerStreamChunk;
use crate::runall::engine::source::{OutputQueue, Source};
use crate::runall::engine::stage::SequencedItem;

/// Buffer size for the initial `BufReader` wrapping the raw file handle.
const FILE_BUFFER_SIZE: usize = 1024 * 1024;

/// Default number of FASTQ records per emitted chunk on the gzip/plain path.
const DEFAULT_RECORDS_PER_BATCH: usize = 400;

/// Default number of raw BGZF blocks per emitted chunk on the BGZF path.
const DEFAULT_BLOCKS_PER_BATCH: usize = 4;

// ============================================================================
// Format detection
// ============================================================================

/// Compression / container format of a FASTQ input file.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FastqFormat {
    /// BGZF-compressed (blocked gzip with 'BC' extra subfield).
    Bgzf,
    /// Plain gzip (single stream, not BGZF).
    Gzip,
    /// Uncompressed plain text.
    Plain,
}

/// Detect the FASTQ container format by inspecting the first 18 bytes.
///
/// BGZF is distinguished from plain gzip by the presence of the `FEXTRA` flag
/// (bit 2 of byte 3) and the `BC` subfield identifier at offset 12.
///
/// # Errors
///
/// Returns an error if the file cannot be opened or read.
pub fn detect_fastq_format(path: &Path) -> Result<FastqFormat> {
    let mut f =
        File::open(path).with_context(|| format!("failed to open FASTQ {}", path.display()))?;
    let mut header = [0u8; 18];
    let n = f.read(&mut header)?;

    // Fewer than 2 bytes cannot be a gzip magic.
    if n < 2 {
        return Ok(FastqFormat::Plain);
    }

    // Gzip magic: 0x1f 0x8b
    if header[0] != 0x1f || header[1] != 0x8b {
        return Ok(FastqFormat::Plain);
    }

    // BGZF signature: FEXTRA flag set AND BC subfield at offset 12.
    if n >= 16 && (header[3] & 0x04) != 0 && header[12] == b'B' && header[13] == b'C' {
        return Ok(FastqFormat::Bgzf);
    }

    Ok(FastqFormat::Gzip)
}

// ============================================================================
// Source
// ============================================================================

/// Multi-file, format-detecting FASTQ source.
///
/// Spawns one reader thread per input path. Each thread pushes
/// [`PerStreamChunk`] items into the shared output queue with a monotonically
/// increasing per-stream `batch_num` and a globally-unique pipeline sequence
/// number.
pub struct FastqFileRead {
    paths: Vec<PathBuf>,
    formats: Vec<FastqFormat>,
    records_per_batch: usize,
    blocks_per_batch: usize,
}

impl FastqFileRead {
    /// Construct a new source by auto-detecting the format of each input path.
    ///
    /// # Errors
    ///
    /// Returns an error if any path cannot be opened for format detection.
    pub fn new(paths: Vec<PathBuf>) -> Result<Self> {
        let formats = paths.iter().map(|p| detect_fastq_format(p)).collect::<Result<Vec<_>>>()?;
        Ok(Self {
            paths,
            formats,
            records_per_batch: DEFAULT_RECORDS_PER_BATCH,
            blocks_per_batch: DEFAULT_BLOCKS_PER_BATCH,
        })
    }

    /// Override the records-per-batch for the gzip/plain path.
    #[must_use]
    pub fn with_records_per_batch(mut self, records_per_batch: usize) -> Self {
        self.records_per_batch = records_per_batch;
        self
    }

    /// Override the blocks-per-batch for the BGZF path.
    #[must_use]
    pub fn with_blocks_per_batch(mut self, blocks_per_batch: usize) -> Self {
        self.blocks_per_batch = blocks_per_batch;
        self
    }

    /// Paths of the input FASTQ files, in stream order.
    #[must_use]
    pub fn paths(&self) -> &[PathBuf] {
        &self.paths
    }

    /// Detected formats of the input FASTQ files, in stream order.
    #[must_use]
    pub fn formats(&self) -> &[FastqFormat] {
        &self.formats
    }
}

/// Per-thread stream reader loop. Emits [`PerStreamChunk`] items to `output`
/// with a globally-unique pipeline sequence number allocated from `next_seq`.
///
/// Records-read accounting: on the gzip/plain path the chunker knows the
/// exact record count (one per FASTQ record, counted from `offsets.len()`)
/// and emits `progress::records_read(n)` per emitted chunk. On the BGZF
/// path the raw blocks are not decompressed here, so the record count is
/// unknown at this layer and `records_read` is not emitted — the dashboard
/// header shows progress via `records_written` and ETA remains "unknown"
/// until the downstream `fastq_parse` stage finishes.
#[allow(clippy::too_many_arguments)]
fn run_stream(
    stream_idx: usize,
    path: &Path,
    fmt: FastqFormat,
    records_per_batch: usize,
    blocks_per_batch: usize,
    output: &Arc<dyn OutputQueue<PerStreamChunk>>,
    cancel: &CancelToken,
    next_seq: &Mutex<u64>,
    total_records: &AtomicU64,
) -> Result<()> {
    match fmt {
        FastqFormat::Bgzf => {
            run_stream_bgzf(stream_idx, path, blocks_per_batch, output, cancel, next_seq)
        }
        FastqFormat::Gzip | FastqFormat::Plain => run_stream_decompressed(
            stream_idx,
            path,
            matches!(fmt, FastqFormat::Gzip),
            records_per_batch,
            output,
            cancel,
            next_seq,
            total_records,
        ),
    }
}

/// BGZF path: read raw blocks and emit them verbatim.
fn run_stream_bgzf(
    stream_idx: usize,
    path: &Path,
    blocks_per_batch: usize,
    output: &Arc<dyn OutputQueue<PerStreamChunk>>,
    cancel: &CancelToken,
    next_seq: &Mutex<u64>,
) -> Result<()> {
    let file = File::open(path)
        .with_context(|| format!("failed to open BGZF FASTQ {}", path.display()))?;
    let mut reader = BufReader::with_capacity(FILE_BUFFER_SIZE, file);

    let mut batch_num: u64 = 0;
    loop {
        if cancel.is_cancelled() {
            return Ok(());
        }

        let blocks = fgumi_bgzf::read_raw_blocks(&mut reader, blocks_per_batch)
            .with_context(|| format!("failed to read BGZF blocks from {}", path.display()))?;
        let at_eof = blocks.is_empty();

        // Concatenate raw block bytes into a single buffer.
        let total_bytes: usize = blocks.iter().map(|b| b.data.len()).sum();
        let mut data: Vec<u8> = Vec::with_capacity(total_bytes);
        for block in blocks {
            data.extend_from_slice(&block.data);
        }

        let chunk = PerStreamChunk { stream_idx, batch_num, data, offsets: None, is_last: at_eof };
        let mem = chunk.data.capacity();

        let seq = {
            let mut guard = next_seq.lock().expect("next_seq mutex poisoned");
            let s = *guard;
            *guard += 1;
            s
        };
        let item = SequencedItem::new(seq, chunk, mem);
        if output.push_until_cancelled(item, cancel).is_err() {
            return Ok(());
        }

        batch_num += 1;
        if at_eof {
            return Ok(());
        }
    }
}

/// Open a decompressing reader for gzip (single-stream) or plain input.
///
/// Because we dispatch on detected file content rather than filename extension,
/// we can't use `fgoxide::io::Io::new_reader` (which keys off `.gz`).
fn open_decompressed_reader(path: &Path, is_gzip: bool) -> Result<Box<dyn BufRead + Send>> {
    let file =
        File::open(path).with_context(|| format!("failed to open FASTQ {}", path.display()))?;
    let buf = BufReader::with_capacity(FILE_BUFFER_SIZE, file);
    if is_gzip {
        Ok(Box::new(BufReader::with_capacity(FILE_BUFFER_SIZE, MultiGzDecoder::new(buf))))
    } else {
        Ok(Box::new(buf))
    }
}

/// Gzip/Plain path: decompress and emit record-aligned chunks.
#[allow(clippy::too_many_arguments)]
fn run_stream_decompressed(
    stream_idx: usize,
    path: &Path,
    is_gzip: bool,
    records_per_batch: usize,
    output: &Arc<dyn OutputQueue<PerStreamChunk>>,
    cancel: &CancelToken,
    next_seq: &Mutex<u64>,
    total_records: &AtomicU64,
) -> Result<()> {
    let mut reader = open_decompressed_reader(path, is_gzip)?;

    let mut batch_num: u64 = 0;
    loop {
        if cancel.is_cancelled() {
            return Ok(());
        }

        let (data, offsets, at_eof) = read_n_fastq_records(&mut reader, records_per_batch)
            .with_context(|| format!("failed to read FASTQ records from {}", path.display()))?;

        // `offsets` uses the sentinel convention: `offsets.len() == records + 1`.
        // A length of 1 (just the leading 0) means zero complete records read.
        if offsets.len() <= 1 && !at_eof {
            continue;
        }

        let record_count = (offsets.len() - 1) as u64;
        let mem = data.capacity() + offsets.capacity() * std::mem::size_of::<usize>();
        let chunk =
            PerStreamChunk { stream_idx, batch_num, data, offsets: Some(offsets), is_last: at_eof };

        let seq = {
            let mut guard = next_seq.lock().expect("next_seq mutex poisoned");
            let s = *guard;
            *guard += 1;
            s
        };
        let item = SequencedItem::new(seq, chunk, mem);
        if output.push_until_cancelled(item, cancel).is_err() {
            return Ok(());
        }

        if record_count > 0 {
            crate::progress::records_read(record_count);
            total_records.fetch_add(record_count, Ordering::Relaxed);
        }

        batch_num += 1;
        if at_eof {
            return Ok(());
        }
    }
}

impl Source for FastqFileRead {
    type Output = PerStreamChunk;

    fn run(
        self: Box<Self>,
        output: Box<dyn OutputQueue<Self::Output>>,
        cancel: CancelToken,
    ) -> Result<()> {
        let this = *self;

        // Share the output queue across per-stream reader threads.
        let output_arc: Arc<dyn OutputQueue<Self::Output>> = Arc::from(output);
        let cancel = Arc::new(cancel);
        let next_seq = Arc::new(Mutex::new(0u64));
        // Accumulates records emitted by the gzip/plain path across all
        // per-stream reader threads. BGZF contributes zero here (record
        // count is not known until downstream decompression).
        let total_records = Arc::new(AtomicU64::new(0));

        let records_per_batch = this.records_per_batch;
        let blocks_per_batch = this.blocks_per_batch;

        let mut handles = Vec::with_capacity(this.paths.len());
        for (stream_idx, (path, fmt)) in
            this.paths.into_iter().zip(this.formats.into_iter()).enumerate()
        {
            let output = output_arc.clone();
            let cancel = cancel.clone();
            let next_seq = next_seq.clone();
            let total_records = total_records.clone();
            let handle = thread::Builder::new()
                .name(format!("fastq_read_{stream_idx}"))
                .spawn(move || -> Result<()> {
                    run_stream(
                        stream_idx,
                        &path,
                        fmt,
                        records_per_batch,
                        blocks_per_batch,
                        &output,
                        &cancel,
                        &next_seq,
                        &total_records,
                    )
                })
                .with_context(|| format!("failed to spawn fastq_read_{stream_idx} thread"))?;
            handles.push(handle);
        }

        // Join all reader threads; on any failure, cancel siblings and capture
        // the first error to return.
        let mut first_err: Option<anyhow::Error> = None;
        for h in handles {
            match h.join() {
                Ok(Ok(())) => {}
                Ok(Err(e)) => {
                    cancel.cancel();
                    if first_err.is_none() {
                        first_err = Some(e);
                    }
                }
                Err(_) => {
                    cancel.cancel();
                    if first_err.is_none() {
                        first_err = Some(anyhow::anyhow!("fastq_read thread panicked"));
                    }
                }
            }
        }

        // Only signal reads_complete on a clean drain — a cancelled/errored
        // run would report a misleading total that doesn't match the number
        // of records the pipeline actually processed.
        if !cancel.is_cancelled() && first_err.is_none() {
            crate::progress::reads_complete(total_records.load(Ordering::Relaxed));
        }

        output_arc.close();
        match first_err {
            Some(e) => Err(e),
            None => Ok(()),
        }
    }

    fn name(&self) -> &'static str {
        "FastqFileRead"
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::runall::engine::memory::MemoryTracker;
    use crate::runall::engine::queue::StageQueue;
    use std::io::Write;

    #[test]
    fn test_detect_plain() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(tmp.path(), b"@r1\nACGT\n+\nIIII\n").unwrap();
        assert_eq!(detect_fastq_format(tmp.path()).unwrap(), FastqFormat::Plain);
    }

    #[test]
    fn test_detect_empty_file() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(tmp.path(), b"").unwrap();
        assert_eq!(detect_fastq_format(tmp.path()).unwrap(), FastqFormat::Plain);
    }

    #[test]
    fn test_detect_gzip() {
        use flate2::Compression;
        use flate2::write::GzEncoder;
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let f = std::fs::File::create(tmp.path()).unwrap();
        let mut enc = GzEncoder::new(f, Compression::default());
        enc.write_all(b"@r1\nACGT\n+\nIIII\n").unwrap();
        enc.finish().unwrap();
        assert_eq!(detect_fastq_format(tmp.path()).unwrap(), FastqFormat::Gzip);
    }

    #[test]
    fn test_read_n_fastq_records_returns_four_records() {
        let content: &[u8] =
            b"@r1\nACGT\n+\nIIII\n@r2\nCGTA\n+\nJJJJ\n@r3\nGATC\n+\nKKKK\n@r4\nTCAG\n+\nLLLL\n";
        let mut reader = std::io::Cursor::new(content);
        let (data, offsets, eof) = read_n_fastq_records(&mut reader, 10).unwrap();
        // Sentinel convention: len = records + 1 (leading 0 + trailing data.len()).
        assert_eq!(offsets.len(), 5);
        assert_eq!(&data, content);
        assert!(eof);
        assert_eq!(offsets, vec![0, 16, 32, 48, 64]);
    }

    #[test]
    fn test_read_n_fastq_records_respects_limit() {
        let content: &[u8] =
            b"@r1\nACGT\n+\nIIII\n@r2\nCGTA\n+\nJJJJ\n@r3\nGATC\n+\nKKKK\n@r4\nTCAG\n+\nLLLL\n";
        let mut reader = std::io::Cursor::new(content);
        let (data, offsets, eof) = read_n_fastq_records(&mut reader, 2).unwrap();
        assert_eq!(offsets.len(), 3);
        assert!(!eof);
        assert_eq!(data.len(), 32);
    }

    #[test]
    fn test_read_n_fastq_records_truncated_errors() {
        // 3 lines of a 4-line record -> truncated.
        let content: &[u8] = b"@r1\nACGT\n+\n";
        let mut reader = std::io::Cursor::new(content);
        let err = read_n_fastq_records(&mut reader, 10).unwrap_err();
        assert!(err.to_string().contains("truncated"));
    }

    #[test]
    fn test_read_n_fastq_records_empty_input() {
        let content: &[u8] = b"";
        let mut reader = std::io::Cursor::new(content);
        let (data, offsets, eof) = read_n_fastq_records(&mut reader, 4).unwrap();
        assert!(data.is_empty());
        // Empty input yields just the leading sentinel `[0]`.
        assert_eq!(offsets, vec![0]);
        assert!(eof);
    }

    /// End-to-end: run a source over a pair of plain FASTQ files and collect
    /// the emitted chunks from the output queue.
    #[test]
    fn test_fastq_file_read_plain_pair() {
        let tmp1 = tempfile::NamedTempFile::new().unwrap();
        let tmp2 = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(tmp1.path(), b"@r1\nACGT\n+\nIIII\n@r2\nCGTA\n+\nJJJJ\n").unwrap();
        std::fs::write(tmp2.path(), b"@r1\nTTTT\n+\nHHHH\n@r2\nGGGG\n+\nGGGG\n").unwrap();

        let paths = vec![tmp1.path().to_path_buf(), tmp2.path().to_path_buf()];
        let source = FastqFileRead::new(paths).unwrap().with_records_per_batch(10);
        assert_eq!(source.formats(), &[FastqFormat::Plain, FastqFormat::Plain]);

        let tracker = Arc::new(MemoryTracker::new(10_000_000));
        let queue: Arc<StageQueue<PerStreamChunk>> =
            Arc::new(StageQueue::new("fastq_out", 64, 1_000_000, tracker));

        let queue_for_source = queue.clone();
        let handle = std::thread::spawn(move || {
            Box::new(source).run(Box::new(queue_for_source), CancelToken::new()).unwrap();
        });
        handle.join().unwrap();

        assert!(queue.is_closed());
        let mut chunks = Vec::new();
        while let Some(item) = queue.pop() {
            chunks.push(item.item);
        }

        // Each stream emits one data chunk (2 records) and one terminal
        // `is_last` chunk. Depending on scheduling the terminal chunk may
        // arrive as a separate empty chunk OR combined with the data chunk
        // — with records_per_batch=10 and 2 records per file, we expect the
        // loop to read everything in one pass (returning at_eof=true), so
        // exactly one chunk per stream.
        assert_eq!(chunks.len(), 2, "expected one chunk per stream, got {chunks:?}");

        // Both streams must be represented with is_last=true.
        let mut by_stream = [false, false];
        for c in &chunks {
            assert!(c.is_last, "chunk should be terminal: {c:?}");
            // Sentinel convention: 2 records → offsets.len() == 3.
            assert_eq!(c.offsets.as_ref().map_or(0, std::vec::Vec::len), 3);
            by_stream[c.stream_idx] = true;
        }
        assert!(by_stream[0] && by_stream[1]);
    }

    #[test]
    fn test_fastq_file_read_gzip() {
        use flate2::Compression;
        use flate2::write::GzEncoder;
        let tmp = tempfile::NamedTempFile::new().unwrap();
        {
            let f = std::fs::File::create(tmp.path()).unwrap();
            let mut enc = GzEncoder::new(f, Compression::default());
            enc.write_all(b"@r1\nACGT\n+\nIIII\n@r2\nCGTA\n+\nJJJJ\n").unwrap();
            enc.finish().unwrap();
        }
        assert_eq!(detect_fastq_format(tmp.path()).unwrap(), FastqFormat::Gzip);

        let source =
            FastqFileRead::new(vec![tmp.path().to_path_buf()]).unwrap().with_records_per_batch(10);

        let tracker = Arc::new(MemoryTracker::new(10_000_000));
        let queue: Arc<StageQueue<PerStreamChunk>> =
            Arc::new(StageQueue::new("fastq_out", 64, 1_000_000, tracker));

        let queue_for_source = queue.clone();
        std::thread::spawn(move || {
            Box::new(source).run(Box::new(queue_for_source), CancelToken::new()).unwrap();
        })
        .join()
        .unwrap();

        assert!(queue.is_closed());
        let mut chunks = Vec::new();
        while let Some(item) = queue.pop() {
            chunks.push(item.item);
        }
        assert_eq!(chunks.len(), 1);
        let c = &chunks[0];
        assert!(c.is_last);
        assert_eq!(c.stream_idx, 0);
        // Sentinel convention: 2 records → offsets.len() == 3.
        assert_eq!(c.offsets.as_ref().unwrap().len(), 3);
    }

    #[test]
    fn test_source_name() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(tmp.path(), b"@r1\nA\n+\nI\n").unwrap();
        let source = FastqFileRead::new(vec![tmp.path().to_path_buf()]).unwrap();
        assert_eq!(Source::name(&source), "FastqFileRead");
    }
}

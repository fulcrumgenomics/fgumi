//! BAM reader factories and related types.
//!
//! This module provides factories for creating BAM readers with consistent
//! error handling, header management, and optional async prefetch.
//!
//! # Threading Model
//!
//! BAM files use BGZF compression, which can be parallelized for reading:
//!
//! - **Single-threaded**: Use `threads=1` (lower overhead, good for small files)
//! - **Multi-threaded**: Use `threads>1` (higher throughput for large files)
//!
//! Multi-threaded reading is beneficial for:
//! - Large BAM files where decompression is a bottleneck
//! - Systems with fast storage (SSD/NVMe) where I/O isn't the limiting factor
//! - Pipelines that process many records in parallel

use anyhow::{Context, Result};
use noodles::sam::Header;
// Use noodles_bgzf for standard BGZF types
use noodles_bgzf::io::{MultithreadedReader, Reader as BgzfReader};
// Use RawBamReader for raw byte access
// (until https://github.com/zaeleus/noodles/pull/373 merges)
use fgumi_raw_bam::RawBamReader;
use std::fs::File;
use std::io::{self, BufRead, Read};
use std::num::NonZero;
use std::path::Path;

use crate::paths::is_stdin_path;

/// Enum wrapping single-threaded and multi-threaded BGZF readers.
///
/// This allows functions to accept either reader type through a unified interface.
pub enum BgzfReaderEnum {
    /// Single-threaded BGZF reader (lower overhead for small files)
    SingleThreaded(BgzfReader<Box<dyn Read + Send>>),
    /// Multi-threaded BGZF reader (noodles built-in threading)
    MultiThreaded(MultithreadedReader<Box<dyn Read + Send>>),
}

impl Read for BgzfReaderEnum {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match self {
            BgzfReaderEnum::SingleThreaded(r) => r.read(buf),
            BgzfReaderEnum::MultiThreaded(r) => r.read(buf),
        }
    }
}

impl BufRead for BgzfReaderEnum {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        match self {
            BgzfReaderEnum::SingleThreaded(r) => r.fill_buf(),
            BgzfReaderEnum::MultiThreaded(r) => r.fill_buf(),
        }
    }

    fn consume(&mut self, amt: usize) {
        match self {
            BgzfReaderEnum::SingleThreaded(r) => r.consume(amt),
            BgzfReaderEnum::MultiThreaded(r) => r.consume(amt),
        }
    }
}

/// Type alias for a BAM reader that supports both single and multi-threaded BGZF.
pub type BamReaderAuto = noodles::bam::io::Reader<BgzfReaderEnum>;

/// Type alias for a raw BAM reader that supports both single and multi-threaded BGZF.
pub type RawBamReaderAuto = RawBamReader<BgzfReaderEnum>;

/// Options controlling how [`create_bam_reader_for_pipeline_with_opts`] opens
/// and wraps its input file.
#[derive(Debug, Clone, Copy, Default)]
pub struct PipelineReaderOpts {
    /// If true, wrap inputs in a `PrefetchReader` background thread
    /// so that disk reads happen on a dedicated I/O thread.
    pub async_reader: bool,
}

/// Create a [`BgzfReaderEnum`] from a reader, selecting single- or multi-threaded based on
/// `threads`.
fn make_bgzf_reader(reader: Box<dyn Read + Send>, threads: usize) -> BgzfReaderEnum {
    if threads > 1 {
        let worker_count = NonZero::new(threads).expect("threads > 1 checked above");
        BgzfReaderEnum::MultiThreaded(MultithreadedReader::with_worker_count(worker_count, reader))
    } else {
        BgzfReaderEnum::SingleThreaded(BgzfReader::new(reader))
    }
}

/// Create a BAM reader and read its header.
///
/// # Arguments
/// * `path` - Path to the input BAM file
/// * `threads` - Number of decompression threads (1 = single-threaded)
///
/// # Threading
/// - `threads <= 1`: Single-threaded (lower overhead, good for small files)
/// - `threads > 1`: Multi-threaded (higher throughput for large files)
///
/// # Returns
/// A tuple of (BAM reader, header)
///
/// # Errors
/// Returns an error if the file cannot be opened or the header cannot be read
///
/// # Panics
/// Panics if `threads > 1` but `NonZero::new` fails (should not happen).
///
/// # Example
/// ```no_run
/// use fgumi_bam_io::reader::create_bam_reader;
/// use std::path::Path;
///
/// // Single-threaded
/// let (mut reader, header) = create_bam_reader(Path::new("input.bam"), 1).unwrap();
///
/// // Multi-threaded with 4 decompression threads
/// let (mut reader, header) = create_bam_reader(Path::new("input.bam"), 4).unwrap();
/// ```
pub fn create_bam_reader<P: AsRef<Path>>(
    path: P,
    threads: usize,
) -> Result<(BamReaderAuto, Header)> {
    create_bam_reader_with_opts(path, threads, PipelineReaderOpts::default())
}

/// Variant of [`create_bam_reader`] that accepts [`PipelineReaderOpts`].
///
/// When `opts.async_reader` is true the file is wrapped in a `PrefetchReader`
/// before the BGZF decompression layer, overlapping disk I/O with BGZF decoding.
///
/// # Errors
///
/// Returns an error if the file cannot be opened or the BAM header cannot be parsed.
pub fn create_bam_reader_with_opts<P: AsRef<Path>>(
    path: P,
    threads: usize,
    opts: PipelineReaderOpts,
) -> Result<(BamReaderAuto, Header)> {
    let path_ref = path.as_ref();
    let bgzf_reader = open_bgzf_reader(path_ref, opts, threads, "BAM reader")?;

    let mut reader = noodles::bam::io::Reader::from(bgzf_reader);
    let header = reader
        .read_header()
        .with_context(|| format!("Failed to read header from: {}", path_ref.display()))?;

    Ok((reader, header))
}

/// Open `path` and return a BGZF reader configured per `opts`.
///
/// Centralizes the file-open + `posix_fadvise` hint + optional `PrefetchReader`
/// wrapping + `make_bgzf_reader` sequence shared by both [`create_bam_reader_with_opts`]
/// and [`create_raw_bam_reader_with_opts`]. `label` distinguishes the two
/// callers in the async-reader log line.
fn open_bgzf_reader(
    path: &Path,
    opts: PipelineReaderOpts,
    threads: usize,
    label: &str,
) -> Result<BgzfReaderEnum> {
    // stdin can't be re-opened: `File::open("-")` fails outright, and
    // `File::open("/dev/stdin")` re-reads fd 0 — losing any bytes a caller
    // already consumed. Read the process stdin directly so this BGZF reader is
    // the sole consumer (fixes the single-threaded double-open data loss).
    if is_stdin_path(path) {
        // Honor --async-reader for stdin too (parity with the file path below):
        // a pipe benefits from prefetch overlapping upstream production with us.
        let reader: Box<dyn Read + Send> = if opts.async_reader {
            log::info!("async {label} enabled: spawning fgumi-prefetch thread for stdin");
            Box::new(crate::prefetch_reader::PrefetchReader::new(io::stdin()))
        } else {
            Box::new(io::stdin())
        };
        let reader = crate::sam_input::normalize_to_bgzf(reader, path)?;
        return Ok(make_bgzf_reader(reader, threads));
    }

    let file = File::open(path)
        .with_context(|| format!("Failed to open input BAM: {}", path.display()))?;
    crate::os_hints::advise_sequential(&file);
    let reader: Box<dyn Read + Send> = if opts.async_reader {
        log::info!(
            "async {} enabled: spawning fgumi-prefetch thread for {}",
            label,
            path.display()
        );
        Box::new(crate::prefetch_reader::PrefetchReader::from_file(file))
    } else {
        Box::new(file)
    };
    let reader = crate::sam_input::normalize_to_bgzf(reader, path)?;
    Ok(make_bgzf_reader(reader, threads))
}

/// Create a raw BAM reader that yields raw bytes instead of noodles Record.
///
/// This is used by the raw sorting pipeline for high-performance byte-level access
/// without going through noodles' Record parsing.
///
/// # Arguments
/// * `path` - Path to the input BAM file
/// * `threads` - Number of threads for BGZF decompression (1 = single-threaded)
///
/// # Returns
/// A tuple of (`raw_reader`, `header`)
///
/// # Errors
/// Returns an error if the file cannot be opened or the header cannot be read
///
/// # Panics
/// Panics if `threads > 1` but `NonZero::new` fails (should not happen).
pub fn create_raw_bam_reader<P: AsRef<Path>>(
    path: P,
    threads: usize,
) -> Result<(RawBamReaderAuto, Header)> {
    create_raw_bam_reader_with_opts(path, threads, PipelineReaderOpts::default())
}

/// Variant of [`create_raw_bam_reader`] that accepts [`PipelineReaderOpts`].
///
/// When `opts.async_reader` is true the file is wrapped in a `PrefetchReader`
/// before the BGZF decompression layer, overlapping disk I/O with BGZF decoding.
///
/// # Errors
///
/// Returns an error if the file cannot be opened or the BAM header cannot be parsed.
pub fn create_raw_bam_reader_with_opts<P: AsRef<Path>>(
    path: P,
    threads: usize,
    opts: PipelineReaderOpts,
) -> Result<(RawBamReaderAuto, Header)> {
    let path_ref = path.as_ref();
    let bgzf_reader = open_bgzf_reader(path_ref, opts, threads, "raw BAM reader")?;

    // Use noodles to read the header, then extract the BGZF reader
    let mut noodles_reader = noodles::bam::io::Reader::from(bgzf_reader);
    let header = noodles_reader
        .read_header()
        .with_context(|| format!("Failed to read header from: {}", path_ref.display()))?;

    // Get back the BGZF reader (header has been consumed)
    let bgzf_reader = noodles_reader.into_inner();

    // Wrap in our raw reader
    let raw_reader = RawBamReader::new(bgzf_reader);

    Ok((raw_reader, header))
}

/// A reader that buffers all bytes read through it.
///
/// Wraps an inner reader and tees every byte read into an internal buffer.
/// Call [`into_parts`](Self::into_parts) to recover both the buffer and the
/// original inner reader so the consumed bytes can be replayed.  Used to
/// capture raw bytes while parsing a BAM header so they can be re-delivered
/// to a downstream reader via [`ChainedReader`].
pub struct TeeReader<R> {
    inner: R,
    buffer: Vec<u8>,
}

impl<R: Read> TeeReader<R> {
    /// Wrap `inner` in a tee reader.
    pub fn new(inner: R) -> Self {
        Self { inner, buffer: Vec::new() }
    }

    /// Consume the `TeeReader` and return the buffered bytes and inner reader.
    pub fn into_parts(self) -> (Vec<u8>, R) {
        (self.buffer, self.inner)
    }
}

impl<R: Read> Read for TeeReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let n = self.inner.read(buf)?;
        self.buffer.extend_from_slice(&buf[..n]);
        Ok(n)
    }
}

/// A reader that chains buffered data with a remaining stream.
///
/// First yields all bytes from the internal buffer, then transparently
/// delegates to the wrapped inner reader.  Typically constructed with bytes
/// captured by a [`TeeReader`] so that bytes consumed during header parsing
/// can be replayed to a downstream reader.
pub struct ChainedReader<R> {
    buffer: io::Cursor<Vec<u8>>,
    inner: R,
    buffer_exhausted: bool,
}

impl<R: Read> ChainedReader<R> {
    /// Create a new chained reader from buffered data and an inner reader.
    pub fn new(buffer: Vec<u8>, inner: R) -> Self {
        Self { buffer: io::Cursor::new(buffer), inner, buffer_exhausted: false }
    }
}

impl<R: Read> Read for ChainedReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if !self.buffer_exhausted {
            let n = self.buffer.read(buf)?;
            if n > 0 {
                return Ok(n);
            }
            self.buffer_exhausted = true;
        }
        self.inner.read(buf)
    }
}

/// Create a raw byte reader and header for pipeline use, supporting both files and stdin.
///
/// This function is designed for commands that need to pass a reader to the pipeline.
/// Unlike `create_bam_reader`, this returns a raw byte reader (not a BAM reader) that
/// can be passed directly to `run_bam_pipeline_*_from_reader` functions.
///
/// For files: Opens the file, reads the header, seeks back to start, returns the file.
/// For stdin: Buffers all bytes read while parsing header, returns a chained reader
///            that first yields the buffered bytes then continues from stdin.
///
/// This is a convenience wrapper that disables the async prefetch reader. Use
/// [`create_bam_reader_for_pipeline_with_opts`] to opt in.
///
/// # Arguments
/// * `path` - Path to the input BAM file, or "-" / "/dev/stdin" for stdin
///
/// # Returns
/// A tuple of (boxed reader, header). The reader is positioned at the start of the file
/// (including the header bytes) so the pipeline can parse the header again.
///
/// # Errors
/// Returns an error if the file cannot be opened or the header cannot be read
///
/// # Example
/// ```no_run
/// use fgumi_bam_io::reader::create_bam_reader_for_pipeline;
/// use std::path::Path;
///
/// // From file
/// let (reader, header) = create_bam_reader_for_pipeline(Path::new("input.bam")).unwrap();
///
/// // From stdin
/// let (reader, header) = create_bam_reader_for_pipeline(Path::new("-")).unwrap();
/// ```
pub fn create_bam_reader_for_pipeline<P: AsRef<Path>>(
    path: P,
) -> Result<(Box<dyn Read + Send>, Header)> {
    create_bam_reader_for_pipeline_with_opts(path, PipelineReaderOpts::default())
}

/// Variant of [`create_bam_reader_for_pipeline`] that accepts tuning options.
///
/// For regular files, `POSIX_FADV_SEQUENTIAL` is applied to the file descriptor
/// on Linux (a no-op elsewhere). If `opts.async_reader` is true the input is
/// wrapped in a `PrefetchReader` background thread — for regular files this also
/// applies kernel WILLNEED hints to the file descriptor before spawning the thread.
///
/// # Errors
///
/// Returns an error if the file cannot be opened or the header cannot be read.
pub fn create_bam_reader_for_pipeline_with_opts<P: AsRef<Path>>(
    path: P,
    opts: PipelineReaderOpts,
) -> Result<(Box<dyn Read + Send>, Header)> {
    let path_ref = path.as_ref();

    // Prefetch wraps the opened input, not the replaying reader: only the
    // `File` itself can carry the kernel WILLNEED hints that
    // `PrefetchReader::from_file` issues, and wrapping later would lose them.
    let reader: Box<dyn Read + Send> = if is_stdin_path(path_ref) {
        if opts.async_reader {
            log::info!("async BAM reader enabled: spawning fgumi-prefetch thread for stdin");
            Box::new(crate::prefetch_reader::PrefetchReader::new(io::stdin()))
        } else {
            Box::new(io::stdin())
        }
    } else {
        let file = File::open(path_ref)
            .with_context(|| format!("Failed to open input BAM: {}", path_ref.display()))?;

        // Tell the kernel to grow the per-fd read-ahead window. Best-effort;
        // failure is logged and ignored. On non-Linux targets this is a no-op.
        crate::os_hints::advise_sequential(&file);

        if opts.async_reader {
            log::info!(
                "async BAM reader enabled: spawning fgumi-prefetch thread for {}",
                path_ref.display()
            );
            Box::new(crate::prefetch_reader::PrefetchReader::from_file(file))
        } else {
            Box::new(file)
        }
    };

    // SAM text is transcoded to BGZF here so the rest of this function — and
    // every pipeline stage downstream of it — sees one stream format.
    let reader = crate::sam_input::normalize_to_bgzf(reader, path_ref)?;
    let (header, rest) = read_header_and_replay(reader, path_ref)?;

    Ok((rest, header))
}

/// Parse the BAM header from `reader`, then return a reader that replays the
/// whole stream from byte zero.
///
/// The bytes the header parse consumed are captured by a [`TeeReader`] and
/// chained back in front of the remainder, so callers get the header without
/// the input having to be seekable. That matters for stdin and FIFOs, and it
/// means the input is opened exactly once.
fn read_header_and_replay(
    reader: Box<dyn Read + Send>,
    path: &Path,
) -> Result<(Header, Box<dyn Read + Send>)> {
    let tee_reader = TeeReader::new(reader);
    let bgzf_reader = BgzfReader::new(tee_reader);
    let mut bam_reader = noodles::bam::io::Reader::from(bgzf_reader);

    let header = bam_reader
        .read_header()
        .with_context(|| format!("Failed to read header from: {}", path.display()))?;

    let (buffered_bytes, rest) = bam_reader.into_inner().into_inner().into_parts();

    Ok((header, Box::new(ChainedReader::new(buffered_bytes, rest))))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::writer::create_bam_writer;
    use noodles::sam::alignment::record::cigar::op::Kind as CigarKind;
    use noodles::sam::header::record::value::{Map, map::ReferenceSequence};
    use std::num::NonZeroUsize;
    use tempfile::NamedTempFile;

    fn create_test_header() -> Header {
        let mut builder = Header::builder();
        let ref_seq = Map::<ReferenceSequence>::new(
            NonZeroUsize::new(100).expect("100 is non-zero constant"),
        );
        builder = builder.add_reference_sequence(b"chr1", ref_seq);
        builder.build()
    }

    /// Minimal uncompressed SAM text: one reference, one aligned record.
    const SAM_TEXT: &str = "@HD\tVN:1.6\tSO:unsorted\n\
                            @SQ\tSN:chr1\tLN:100\n\
                            r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n";

    /// Writes [`SAM_TEXT`] to a temp file with the given suffix and returns it.
    fn write_sam_text(suffix: &str) -> NamedTempFile {
        let file =
            tempfile::Builder::new().suffix(suffix).tempfile().expect("create temp SAM file");
        std::fs::write(file.path(), SAM_TEXT).expect("write SAM text");
        file
    }

    /// Asserts a decoded header + records match [`SAM_TEXT`] field for field, so a
    /// transcode that silently corrupted a field cannot pass on record count alone.
    fn assert_matches_sam_text(header: &noodles::sam::Header, records: &[noodles::bam::Record]) {
        let refs = header.reference_sequences();
        assert_eq!(refs.len(), 1, "header @SQ lost");
        let (name, seq) = refs.get_index(0).expect("one @SQ");
        assert_eq!(AsRef::<[u8]>::as_ref(name), b"chr1", "@SQ SN");
        assert_eq!(usize::from(seq.length()), 100, "@SQ LN");

        assert_eq!(records.len(), 1, "record count");
        let r = &records[0];
        assert_eq!(r.name().map(std::convert::AsRef::as_ref), Some(b"r1" as &[u8]), "QNAME");
        assert_eq!(r.flags().bits(), 0, "FLAG");
        assert_eq!(r.alignment_start().expect("mapped").expect("pos").get(), 1, "POS (1-based)");
        let cigar: Vec<_> = r
            .cigar()
            .iter()
            .map(|op| op.map(|o| (o.kind(), o.len())))
            .collect::<io::Result<_>>()
            .expect("cigar");
        assert_eq!(cigar, [(CigarKind::Match, 4)], "CIGAR");
        assert_eq!(r.sequence().iter().collect::<Vec<u8>>(), b"ACGT", "SEQ");
        assert_eq!(r.quality_scores().as_ref(), &[40, 40, 40, 40], "QUAL ('IIII' = Q40)");
    }

    /// Every reader factory must accept uncompressed SAM, not just BGZF BAM:
    /// a command that takes BAM takes SAM. The typed reader is the surface
    /// used by whole-record consumers.
    #[test]
    fn test_create_bam_reader_accepts_uncompressed_sam() -> Result<()> {
        let input = write_sam_text(".sam");

        let (mut reader, header) = create_bam_reader(input.path(), 1)?;

        let records = reader.records().collect::<io::Result<Vec<_>>>()?;
        assert_matches_sam_text(&header, &records);
        Ok(())
    }

    /// The raw-byte reader backs the single-threaded fast paths, so it must
    /// accept SAM on the same terms as the typed reader.
    #[test]
    fn test_create_raw_bam_reader_accepts_uncompressed_sam() -> Result<()> {
        let input = write_sam_text(".sam");

        let (mut reader, header) = create_raw_bam_reader(input.path(), 1)?;

        assert_eq!(header.reference_sequences().len(), 1, "header @SQ lost");
        let mut record = fgumi_raw_bam::RawRecord::new();
        assert!(reader.read_record(&mut record)? > 0, "expected one record");
        // Identity, not just presence: a transcode that dropped a field would
        // still return a record, so pin the fields the SAM text specifies.
        let view = record.view();
        assert_eq!(view.read_name(), b"r1", "QNAME");
        assert_eq!(view.flags(), 0, "FLAG");
        assert_eq!(view.ref_id(), 0, "RNAME -> first @SQ");
        assert_eq!(view.pos(), 0, "POS (0-based in BAM)");
        assert_eq!(reader.read_record(&mut record)?, 0, "expected EOF");
        Ok(())
    }

    /// The pipeline reader hands the caller a byte stream that downstream
    /// stages BGZF-decode themselves, so the SAM path must present the same
    /// compressed-stream contract.
    #[test]
    fn test_create_bam_reader_for_pipeline_accepts_uncompressed_sam() -> Result<()> {
        let input = write_sam_text(".sam");

        let (stream, header) = create_bam_reader_for_pipeline(input.path())?;

        assert_eq!(header.reference_sequences().len(), 1, "header @SQ lost");
        let mut reader = noodles::bam::io::Reader::new(io::BufReader::new(stream));
        // The pipeline hands back a raw BGZF byte stream; decode it and assert the
        // replayed header + record match the SAM text, not merely that one arrived.
        let stream_header = reader.read_header()?;
        let records = reader.records().collect::<io::Result<Vec<_>>>()?;
        assert_matches_sam_text(&stream_header, &records);
        Ok(())
    }

    #[test]
    fn test_create_bam_reader_nonexistent_file() {
        let result = create_bam_reader("/nonexistent/file.bam", 1);
        assert!(result.is_err());
        if let Err(e) = result {
            let err_msg = e.to_string();
            assert!(err_msg.contains("Failed to open input BAM"));
        }
    }

    #[test]
    fn test_roundtrip_write_and_read() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let header = create_test_header();

        // Write (single-threaded)
        {
            let _writer = create_bam_writer(temp_file.path(), &header, 1, 6)?;
            // Writer is dropped, file is written
        }

        // Read (single-threaded)
        let (mut reader, read_header) = create_bam_reader(temp_file.path(), 1)?;

        // Verify header has our reference sequence
        assert_eq!(read_header.reference_sequences().len(), 1);

        // Verify we can iterate (even though there are no records)
        let records: Result<Vec<_>, _> = reader.records().collect();
        assert!(records.is_ok());

        Ok(())
    }

    #[test]
    fn test_roundtrip_write_and_read_multithreaded() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let header = create_test_header();

        // Write (multi-threaded)
        {
            let _writer = create_bam_writer(temp_file.path(), &header, 4, 6)?;
            // Writer is dropped, file is written
        }

        // Read (multi-threaded)
        let (mut reader, read_header) = create_bam_reader(temp_file.path(), 4)?;

        // Verify header has our reference sequence
        assert_eq!(read_header.reference_sequences().len(), 1);

        // Verify we can iterate (even though there are no records)
        let records: Result<Vec<_>, _> = reader.records().collect();
        assert!(records.is_ok());

        Ok(())
    }

    #[test]
    fn test_create_bam_reader_for_pipeline_from_file() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let header = create_test_header();

        // Write a BAM file first
        {
            let _writer = create_bam_writer(temp_file.path(), &header, 1, 6)?;
        }

        // Read using create_bam_reader_for_pipeline
        let (mut reader, read_header) = create_bam_reader_for_pipeline(temp_file.path())?;
        assert_eq!(read_header.reference_sequences().len(), 1);

        // The reader returns raw bytes, verify we can read something
        let mut buf = [0u8; 16];
        let n = reader.read(&mut buf)?;
        assert!(n > 0, "Should read some bytes from the file");

        Ok(())
    }

    #[test]
    fn test_create_bam_reader_for_pipeline_nonexistent_file() {
        let result = create_bam_reader_for_pipeline("/nonexistent/file.bam");
        assert!(result.is_err());
        if let Err(e) = result {
            let err_msg = e.to_string();
            assert!(err_msg.contains("Failed to open input BAM"));
        }
    }

    #[test]
    fn test_chained_reader() {
        let buffer = vec![1, 2, 3, 4, 5];
        let remaining = io::Cursor::new(vec![6, 7, 8, 9, 10]);
        let mut chained = ChainedReader::new(buffer, remaining);

        let mut result = Vec::new();
        chained.read_to_end(&mut result).expect("read_to_end should succeed");

        assert_eq!(result, vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10]);
    }

    #[test]
    fn test_chained_reader_empty_buffer() {
        let buffer = vec![];
        let remaining = io::Cursor::new(vec![1, 2, 3]);
        let mut chained = ChainedReader::new(buffer, remaining);

        let mut result = Vec::new();
        chained.read_to_end(&mut result).expect("read_to_end should succeed");

        assert_eq!(result, vec![1, 2, 3]);
    }

    #[test]
    fn test_chained_reader_empty_remaining() {
        let buffer = vec![1, 2, 3];
        let remaining = io::Cursor::new(vec![]);
        let mut chained = ChainedReader::new(buffer, remaining);

        let mut result = Vec::new();
        chained.read_to_end(&mut result).expect("read_to_end should succeed");

        assert_eq!(result, vec![1, 2, 3]);
    }

    #[test]
    fn test_create_bam_reader_for_pipeline_with_async_reader() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let header = create_test_header();

        // Write a BAM file first
        {
            let _writer = create_bam_writer(temp_file.path(), &header, 1, 6)?;
        }

        // Read using async reader opts — exercises the PrefetchReader branch
        let opts = PipelineReaderOpts { async_reader: true };
        let (mut reader, read_header) =
            create_bam_reader_for_pipeline_with_opts(temp_file.path(), opts)?;
        assert_eq!(read_header.reference_sequences().len(), 1);

        // The reader returns raw bytes; verify it is usable
        let mut buf = [0u8; 16];
        let n = reader.read(&mut buf)?;
        assert!(n > 0, "Should read some bytes from the async reader");

        // Ensure read_to_end works through the PrefetchReader
        let mut rest = Vec::new();
        reader.read_to_end(&mut rest)?;

        Ok(())
    }
}

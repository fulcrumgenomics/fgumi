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
use std::io::{self, BufRead, BufReader, Read};
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
#[derive(Debug, Clone, Copy)]
pub struct PipelineReaderOpts {
    /// If true, wrap inputs in a `PrefetchReader` background thread
    /// so that disk reads happen on a dedicated I/O thread.
    pub async_reader: bool,
    /// Whether the BGZF decode path should verify each block's CRC32
    /// checksum against its footer.
    ///
    /// This struct only opens the byte stream (or, for
    /// [`create_bam_reader_with_opts`]/[`create_raw_bam_reader_with_opts`],
    /// wraps it in noodles' own BGZF reader, which always verifies and has no
    /// skip knob) — it does not itself decode BGZF blocks, so **this field is
    /// not consumed today**. The multi-threaded unified pipeline gets its CRC
    /// policy from `PipelineConfig::verify_crc` (set in `build_pipeline_config`
    /// in `fgumi_lib`), not from here. This field becomes live once the
    /// raw-reader path is unified onto fgumi-bgzf (see #800), at which point
    /// [`create_raw_bam_reader_with_opts`] will honor it directly. It is kept
    /// now so callers already bundle the setting in one place. Defaults to
    /// `true` (verify) — the safe, pre-existing behavior.
    pub verify_crc: bool,
}

impl Default for PipelineReaderOpts {
    fn default() -> Self {
        Self { async_reader: false, verify_crc: true }
    }
}

/// Wrap `reader` in a BGZF decoder, multi-threaded when `threads > 1`.
///
/// Exposed so callers that must open their input as a stream — a pipe, a FIFO,
/// process substitution — can still get multi-threaded BGZF decode. Opening by
/// path is not a prerequisite for parallel decompression.
///
/// # Panics
///
/// Panics if `threads > 1` but `NonZero::new` fails, which cannot happen.
#[must_use]
pub fn make_bgzf_reader(reader: Box<dyn Read + Send>, threads: usize) -> BgzfReaderEnum {
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

/// Fill `buf` from `reader`, tolerating short reads and interruptions.
///
/// Returns how many bytes were read, which is less than `buf.len()` only at end
/// of input. Every format sniff in fgumi needs this: a single `read` on a pipe
/// may return fewer bytes than asked for even when more are coming (the upstream
/// producer has not flushed yet), and a signal during the read must not fail the
/// run. Classifying an input from a short first `read` is how BAM-on-stdin was
/// once misread as SAM.
///
/// The bytes are returned in `buf` rather than pushed back, so the caller can
/// replay them with a [`ChainedReader`] — the input may be a pipe, which cannot
/// be rewound.
///
/// # Errors
///
/// Returns any read error other than [`io::ErrorKind::Interrupted`].
pub fn read_prefix<R: Read>(reader: &mut R, buf: &mut [u8]) -> io::Result<usize> {
    let mut filled = 0;
    while filled < buf.len() {
        match reader.read(&mut buf[filled..]) {
            Ok(0) => break,
            Ok(n) => filled += n,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    Ok(filled)
}

/// Open `path` as a normalized BGZF byte stream — the one place a normalized
/// stream is minted.
///
/// Handles stdin (which cannot be re-opened: `File::open("-")` fails outright
/// and `File::open("/dev/stdin")` re-reads fd 0, losing bytes a caller already
/// consumed), applies the sequential-readahead hint to regular files, honors
/// `opts.async_reader`, and transcodes SAM text so callers only ever see BGZF.
///
/// Prefetch wraps the *opened* input rather than the normalized stream: only the
/// `File` itself can carry the kernel WILLNEED hints `PrefetchReader::from_file`
/// issues, and wrapping later would lose them.
///
/// `label` names this reader in the async-reader log line.
///
/// # Errors
///
/// Returns an error if the input cannot be opened, is empty, or is neither BAM
/// nor SAM.
fn open_normalized_with_opts(
    path: &Path,
    opts: PipelineReaderOpts,
    label: &str,
) -> Result<Box<dyn Read + Send>> {
    let opened: Box<dyn Read + Send> = if is_stdin_path(path) {
        if opts.async_reader {
            log::info!("async {label} enabled: spawning fgumi-prefetch thread for stdin");
            Box::new(crate::prefetch_reader::PrefetchReader::new(io::stdin()))
        } else {
            Box::new(io::stdin())
        }
    } else {
        let file = File::open(path)
            .with_context(|| format!("Failed to open input BAM: {}", path.display()))?;

        // Tell the kernel to grow the per-fd read-ahead window. Best-effort;
        // failure is logged and ignored. On non-Linux targets this is a no-op.
        crate::os_hints::advise_sequential(&file);

        if opts.async_reader {
            log::info!(
                "async {} enabled: spawning fgumi-prefetch thread for {}",
                label,
                path.display()
            );
            Box::new(crate::prefetch_reader::PrefetchReader::from_file(file))
        } else {
            Box::new(file)
        }
    };

    crate::sam_input::normalize_to_bgzf(opened, path)
}

/// Open `path` as a BGZF byte stream, transcoding it if it is uncompressed SAM.
///
/// The opener for consumers that BGZF-decode the bytes themselves — the raw
/// record readers in `fgumi-sort` and `compare` — rather than going through
/// [`create_bam_reader`] and friends.
///
/// # Errors
///
/// Returns an error if the input cannot be opened, is empty, or is neither BAM
/// nor SAM.
pub fn open_normalized_input<P: AsRef<Path>>(path: P) -> Result<Box<dyn Read + Send>> {
    open_normalized_with_opts(path.as_ref(), PipelineReaderOpts::default(), "BAM reader")
}

/// Buffer placed under the BGZF frame reader by [`open_bgzf_reader`], and under
/// the header parse in [`read_header_and_replay`].
///
/// Sized to hold several maximum-size (64 KiB) BGZF blocks, so a block read is
/// usually served from memory. Deliberately not the 2 MiB the sort engine gives
/// its own input: this buffer is per open reader, and `merge` holds one per input
/// file.
///
/// The two users differ in lifetime, not size: [`open_bgzf_reader`] leaves the
/// buffer resident under the frame reader, while [`read_header_and_replay`]
/// discards it once the header is parsed.
const BGZF_INPUT_BUFFER_SIZE: usize = 256 * 1024;

/// Open `path` and return a BGZF reader configured per `opts`.
///
/// Centralizes the file-open + `posix_fadvise` hint + optional `PrefetchReader`
/// wrapping + `make_bgzf_reader` sequence shared by both [`create_bam_reader_with_opts`]
/// and [`create_raw_bam_reader_with_opts`]. `label` distinguishes the two
/// callers in the async-reader log line.
///
/// The normalized stream is buffered before the BGZF reader sees it. The frame
/// reader asks for an 18-byte header and then a block body, so an unbuffered
/// source pays a syscall per read — and on stdin rather worse, since `io::Stdin`
/// takes a mutex and serves through an 8 KiB buffer on every call. `zipper` and
/// the sort engine's pool reader each already buffer for exactly this reason;
/// this is the same fix for the readers that go through here.
///
/// The buffer goes *under* the frame reader rather than around the raw file so it
/// is not paid twice: [`open_normalized_input`]'s consumers wrap the stream in
/// their own [`BufReader`] (`fgumi_sort::RawBamRecordReader`), and buffering
/// inside the shared opener would copy every byte through two buffers on that
/// path. A `PrefetchReader` already hands out large chunks, so `opts.async_reader`
/// gains nothing here and pays the same extra copy — it is left unwrapped.
fn open_bgzf_reader(
    path: &Path,
    opts: PipelineReaderOpts,
    threads: usize,
    label: &str,
) -> Result<BgzfReaderEnum> {
    let normalized = open_normalized_with_opts(path, opts, label)?;
    let buffered: Box<dyn Read + Send> = if opts.async_reader {
        normalized
    } else {
        Box::new(BufReader::with_capacity(BGZF_INPUT_BUFFER_SIZE, normalized))
    };
    Ok(make_bgzf_reader(buffered, threads))
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
///
/// This is [`std::io::Chain`] plus one property that matters here: the replay
/// buffer is **released** once drained. `Chain` holds its first reader for the
/// life of the stream, which on a many-contig BAM header is tens of megabytes
/// pinned for the whole run.
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
            // A zero-length `buf` yields zero bytes from any source, so it must
            // not be read as "the replay buffer is drained" — doing so would
            // abandon the buffered bytes and silently truncate the stream.
            if buf.is_empty() {
                return Ok(0);
            }

            let n = self.buffer.read(buf)?;
            if n > 0 {
                return Ok(n);
            }
            self.buffer_exhausted = true;
            // Replayed bytes are dead once delivered. Release the allocation
            // rather than holding a BAM header's worth of memory for the
            // lifetime of the stream.
            self.buffer = io::Cursor::new(Vec::new());
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

    // SAM text is transcoded inside the opener, so the rest of this function —
    // and every pipeline stage downstream of it — sees one stream format.
    let reader = open_normalized_with_opts(path_ref, opts, "BAM reader")?;
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
///
/// `reader` must already be a BGZF stream — call
/// [`normalize_to_bgzf`](crate::sam_input::normalize_to_bgzf) first if the
/// input may be uncompressed SAM.
///
/// The header parse itself is buffered, because the BGZF frame reader asks for
/// an 18-byte block header and then a block body: against an unbuffered source
/// every block costs a syscall, and on a many-contig header that is thousands
/// of them. Callers hand this the raw [`open_normalized_input`] output — a
/// plain `File` or `io::stdin()` — so the buffer belongs here rather than
/// repeated at each call site.
///
/// The buffer is *temporary*, and deliberately so. It is placed under the
/// [`TeeReader`] for the parse, then dismantled: its readahead is moved into
/// the replay buffer and the original stream is handed back unwrapped. Leaving
/// it in place would put a second buffer under consumers that already have one
/// — `fgumi_sort::RawBamRecordReader` and the sort worker pool both read BGZF
/// block by block, so a resident layer here would copy every byte of the input
/// a second time on the hot path. Draining it into the replay buffer costs one
/// memcpy of at most `BGZF_INPUT_BUFFER_SIZE`, which [`ChainedReader`] frees
/// as soon as it is delivered.
///
/// # Errors
///
/// Returns an error if the BAM header cannot be read from `reader`.
pub fn read_header_and_replay(
    reader: Box<dyn Read + Send>,
    path: &Path,
) -> Result<(Header, Box<dyn Read + Send>)> {
    let tee_reader = TeeReader::new(BufReader::with_capacity(BGZF_INPUT_BUFFER_SIZE, reader));
    let bgzf_reader = BgzfReader::new(tee_reader);
    let mut bam_reader = noodles::bam::io::Reader::from(bgzf_reader);

    let header = bam_reader
        .read_header()
        .with_context(|| format!("Failed to read header from: {}", path.display()))?;

    let (mut replay, buffered) = bam_reader.into_inner().into_inner().into_parts();

    // The buffer read past the header. Those bytes never went through the tee,
    // and on a pipe they cannot be re-read from the source, so they must be
    // appended to the replay buffer before the buffer is discarded — otherwise
    // the stream silently loses everything the last fill read ahead.
    replay.extend_from_slice(buffered.buffer());
    let rest = buffered.into_inner();

    Ok((header, Box::new(ChainedReader::new(replay, rest))))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::writer::create_bam_writer;
    use noodles::sam::alignment::record::cigar::op::Kind as CigarKind;
    use noodles::sam::header::record::value::{Map, map::ReferenceSequence};
    use std::num::NonZeroUsize;
    use std::sync::Arc;
    use std::sync::atomic::{AtomicUsize, Ordering};
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

    /// A zero-length read must not be mistaken for "the replay buffer is
    /// drained". Every fgumi input now passes through two nested
    /// `ChainedReader`s — one holding the sniffed magic bytes, one holding the
    /// BAM header — so swallowing the buffer here corrupts the stream for any
    /// consumer that happens to hand down an empty slice.
    #[test]
    fn test_chained_reader_zero_length_read_preserves_buffer() {
        let buffer = vec![1, 2, 3];
        let remaining = io::Cursor::new(vec![4, 5]);
        let mut chained = ChainedReader::new(buffer, remaining);

        assert_eq!(chained.read(&mut []).expect("empty read should succeed"), 0);

        let mut result = Vec::new();
        chained.read_to_end(&mut result).expect("read_to_end should succeed");

        assert_eq!(result, vec![1, 2, 3, 4, 5], "replay buffer lost after a zero-length read");
    }

    /// The replayed bytes are dead once delivered. Holding them costs a BAM
    /// header's worth of memory for the lifetime of the stream — tens of MB on
    /// a many-contig header — where the previous seek-based path buffered
    /// nothing.
    #[test]
    fn test_chained_reader_frees_buffer_once_exhausted() {
        let buffer = vec![7u8; 4096];
        let remaining = io::Cursor::new(vec![1, 2, 3]);
        let mut chained = ChainedReader::new(buffer, remaining);

        let mut result = Vec::new();
        chained.read_to_end(&mut result).expect("read_to_end should succeed");

        assert_eq!(result.len(), 4099);
        assert_eq!(
            chained.buffer.get_ref().capacity(),
            0,
            "replay buffer still resident after exhaustion"
        );
    }

    /// Transcodes [`SAM_TEXT`] into an in-memory BGZF stream — the input
    /// [`read_header_and_replay`] contractually expects.
    fn sam_text_as_bgzf() -> Vec<u8> {
        let source: Box<dyn Read + Send> = Box::new(io::Cursor::new(SAM_TEXT.as_bytes().to_vec()));
        let mut normalized = crate::sam_input::normalize_to_bgzf(source, Path::new("test.sam"))
            .expect("transcode SAM text to BGZF");
        let mut bytes = Vec::new();
        normalized.read_to_end(&mut bytes).expect("read transcoded BGZF");
        bytes
    }

    /// Records the largest `read` request that reaches the wrapped source, so a
    /// test can tell whether a buffer sits directly on top of it.
    struct ReadSizeSpy<R> {
        inner: R,
        max_request: Arc<AtomicUsize>,
    }

    impl<R: Read> Read for ReadSizeSpy<R> {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            self.max_request.fetch_max(buf.len(), Ordering::Relaxed);
            self.inner.read(buf)
        }
    }

    /// The BGZF frame reader asks for an 18-byte block header and then a block
    /// body, so an unbuffered source pays a syscall per block. `open_bgzf_reader`
    /// buffers for exactly this reason, and this is the *other* entry point onto a
    /// BGZF stream: its callers — the pipeline reader here and
    /// `fgumi_sort::open_raw_bam_record_reader_with_header` — hand it the raw
    /// `open_normalized_*` output, a plain `File` or `io::stdin()`.
    #[test]
    fn test_read_header_and_replay_buffers_the_source() -> Result<()> {
        let max_request = Arc::new(AtomicUsize::new(0));
        let source = Box::new(ReadSizeSpy {
            inner: io::Cursor::new(sam_text_as_bgzf()),
            max_request: Arc::clone(&max_request),
        });

        let (_header, _rest) = read_header_and_replay(source, Path::new("test.bam"))?;

        assert_eq!(
            max_request.load(Ordering::Relaxed),
            BGZF_INPUT_BUFFER_SIZE,
            "header parse reached the source unbuffered — every BGZF block costs a syscall"
        );
        Ok(())
    }

    /// The header parse's buffer is discarded before the stream is handed back,
    /// so its readahead must first be moved into the replay buffer: those bytes
    /// never went through the [`TeeReader`] and, on a pipe, cannot be re-read
    /// from the source. Drop them and the stream loses its first records with no
    /// error — the failure mode this asserts against, byte for byte.
    #[test]
    fn test_read_header_and_replay_replays_every_byte() -> Result<()> {
        let bgzf = sam_text_as_bgzf();
        let source = Box::new(io::Cursor::new(bgzf.clone()));

        let (header, mut rest) = read_header_and_replay(source, Path::new("test.bam"))?;
        assert_eq!(header.reference_sequences().len(), 1, "header @SQ lost");

        let mut replayed = Vec::new();
        rest.read_to_end(&mut replayed)?;
        assert_eq!(replayed, bgzf, "replayed stream is not byte-identical to the input");
        Ok(())
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
        let opts = PipelineReaderOpts { async_reader: true, ..PipelineReaderOpts::default() };
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

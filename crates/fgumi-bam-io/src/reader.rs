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
    /// Single-threaded BGZF reader backed by fgumi-bgzf, which honors the
    /// `verify_crc` policy (see [`FgumiBgzfReader`]). This is the arm used by
    /// the raw-reader path (`create_raw_bam_reader[_with_opts]`) for
    /// single-threaded input, so `--check-crc`/`--no-check-crc` take effect
    /// there.
    Fgumi(FgumiBgzfReader),
}

impl Read for BgzfReaderEnum {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match self {
            BgzfReaderEnum::SingleThreaded(r) => r.read(buf),
            BgzfReaderEnum::MultiThreaded(r) => r.read(buf),
            BgzfReaderEnum::Fgumi(r) => r.read(buf),
        }
    }
}

impl BufRead for BgzfReaderEnum {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        match self {
            BgzfReaderEnum::SingleThreaded(r) => r.fill_buf(),
            BgzfReaderEnum::MultiThreaded(r) => r.fill_buf(),
            BgzfReaderEnum::Fgumi(r) => r.fill_buf(),
        }
    }

    fn consume(&mut self, amt: usize) {
        match self {
            BgzfReaderEnum::SingleThreaded(r) => r.consume(amt),
            BgzfReaderEnum::MultiThreaded(r) => r.consume(amt),
            BgzfReaderEnum::Fgumi(r) => r.consume(amt),
        }
    }
}

/// Number of raw BGZF blocks the [`FgumiBgzfReader`] *frames* (reads and length
/// -validates) per refill of its frame queue.
///
/// Framing is decoupled from decoding: a batch this size is pulled from the
/// source at once, then decoded one block at a time on demand (see
/// [`FgumiBgzfReader::ensure_buffered`]). Two properties motivate the batch:
///
/// - **Correct EOF handling across concatenated streams.**
///   [`read_raw_blocks`](fgumi_bgzf::read_raw_blocks) skips BGZF EOF-marker
///   blocks and returns an empty vector only when it finds no real block within
///   the batch. A multithreaded writer (and `cat a.bam b.bam`) emits an EOF
///   marker between segments, so framing one block at a time would read that
///   lone marker, get an empty vector, and wrongly conclude end-of-stream before
///   the following data. A batch absorbs intermediate markers and returns real
///   blocks, so an empty vector reliably means true end of input.
/// - **Bounded memory.** Framed blocks are compressed (typically ~16 KiB), so a
///   batch of this size costs a few hundred KiB per reader — `merge` holds one
///   reader per input file, so this is kept modest rather than file-sized.
const FGUMI_FRAME_BATCH: usize = 16;

/// A single-threaded, safe streaming BGZF decoder built over fgumi-bgzf.
///
/// Unlike noodles' BGZF reader (which always verifies CRC32 and exposes no skip
/// knob), this decoder honors a `verify_crc` flag: with it `false`, block CRC32
/// verification is skipped (the decompressed-size check always runs), which is
/// the fast path for trusted input. It yields the **decompressed** BAM byte
/// stream, so it slots in wherever a noodles BGZF reader did — the header parse
/// and [`RawBamReader`] both see plain decompressed bytes.
///
/// Framing and decoding are decoupled. Blocks are framed in batches via
/// [`read_raw_blocks`](fgumi_bgzf::read_raw_blocks) (so intermediate BGZF EOF
/// markers are handled — see `FGUMI_FRAME_BATCH`) into a queue, and decoded
/// one at a time with
/// [`decompress_block_into_opts`](fgumi_bgzf::decompress_block_into_opts) into an
/// internal buffer served out across `read`/`fill_buf` calls with a cursor.
/// Decoding lazily (one queued block per refill) keeps decode work — and any
/// CRC/size error — tied to the bytes the caller actually reads, rather than
/// pulling a later block's error forward into, say, the header parse. Partial
/// reads and mid-block boundaries are handled by the cursor. No `unsafe` is used.
pub struct FgumiBgzfReader {
    /// The underlying compressed byte source (a normalized BGZF stream).
    inner: Box<dyn Read + Send>,
    /// Reusable libdeflater decompressor.
    decompressor: fgumi_bgzf::Decompressor,
    /// Whether to verify each block's CRC32 against its footer.
    verify_crc: bool,
    /// Framed-but-not-yet-decoded blocks, filled a [`FGUMI_FRAME_BATCH`] at a
    /// time and drained front to back.
    frames: std::collections::VecDeque<fgumi_bgzf::RawBgzfBlock>,
    /// Set once the source yields no further blocks, so framing stops.
    frames_done: bool,
    /// Decoded bytes of the current block not yet served to the caller.
    buffer: Vec<u8>,
    /// Cursor into `buffer`: `buffer[pos..]` is unread.
    pos: usize,
    /// Set once a block fails to decode, poisoning the reader. A decode error
    /// pops its frame and resets the cursor, so without this flag a later
    /// `read`/`fill_buf` would resume at the *next* frame and silently skip the
    /// corrupted block's records. Once set, every subsequent refill errors.
    failed: bool,
}

impl FgumiBgzfReader {
    /// Wrap `inner` in a streaming fgumi-bgzf decoder.
    ///
    /// `verify_crc` selects whether each block's CRC32 is checked; the
    /// decompressed-size check always runs regardless.
    #[must_use]
    pub fn new(inner: Box<dyn Read + Send>, verify_crc: bool) -> Self {
        Self {
            inner,
            decompressor: fgumi_bgzf::Decompressor::new(),
            verify_crc,
            frames: std::collections::VecDeque::new(),
            frames_done: false,
            buffer: Vec::new(),
            pos: 0,
            failed: false,
        }
    }

    /// Ensure `self.buffer[self.pos..]` holds at least one byte, framing and
    /// decoding more blocks if needed. Returns `Ok(false)` at end of input.
    ///
    /// Decodes at most one queued block per iteration and loops because a block
    /// can decode to zero bytes (e.g. `ISIZE == 0`) without being end of input;
    /// it stops only when the source yields no further blocks.
    fn ensure_buffered(&mut self) -> io::Result<bool> {
        if self.failed {
            return Err(io::Error::other(
                "fgumi-bgzf reader poisoned by an earlier BGZF block decode failure",
            ));
        }
        while self.pos >= self.buffer.len() {
            if self.frames.is_empty() {
                if self.frames_done {
                    return Ok(false);
                }
                let batch = fgumi_bgzf::read_raw_blocks(&mut self.inner, FGUMI_FRAME_BATCH)?;
                if batch.is_empty() {
                    // No real block in the batch => true end of input (EOF
                    // markers are skipped inside `read_raw_blocks`).
                    self.frames_done = true;
                    return Ok(false);
                }
                self.frames.extend(batch);
            }

            let block = self.frames.pop_front().expect("frames non-empty checked above");
            self.buffer.clear();
            self.pos = 0;
            // Poison the reader before propagating: this frame is already popped
            // and the cursor reset, so a retried refill would otherwise skip
            // straight to the next frame and hide the corrupted block.
            if let Err(e) = fgumi_bgzf::decompress_block_into_opts(
                &block,
                &mut self.decompressor,
                &mut self.buffer,
                self.verify_crc,
            ) {
                self.failed = true;
                return Err(e);
            }
        }
        Ok(true)
    }
}

impl Read for FgumiBgzfReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if buf.is_empty() {
            return Ok(0);
        }
        if !self.ensure_buffered()? {
            return Ok(0);
        }
        let available = &self.buffer[self.pos..];
        let n = available.len().min(buf.len());
        buf[..n].copy_from_slice(&available[..n]);
        self.pos += n;
        Ok(n)
    }
}

impl BufRead for FgumiBgzfReader {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        if !self.ensure_buffered()? {
            return Ok(&[]);
        }
        Ok(&self.buffer[self.pos..])
    }

    fn consume(&mut self, amt: usize) {
        self.pos = (self.pos + amt).min(self.buffer.len());
    }
}

/// Type alias for a BAM reader that supports both single and multi-threaded BGZF.
pub type BamReaderAuto = noodles::bam::io::Reader<BgzfReaderEnum>;

/// Type alias for a raw BAM reader that supports both single and multi-threaded BGZF.
pub type RawBamReaderAuto = RawBamReader<BgzfReaderEnum>;

/// Read-stream policy for a seekable BAM/SAM source: how many concurrent
/// positional reads to issue per fill window.
///
/// Ports the semantics of fgumi v0.7.0's `fgumi sort --read-streams` flag. The
/// mechanism (concurrent positional reads that raise the device's read queue
/// depth) lives in [`crate::scatter_reader`]; on a slow, deep-queue device
/// (e.g. EBS gp3) issuing several reads at once is markedly faster than the
/// single outstanding read a plain reader issues.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ReadStreams {
    /// Measure the device's single-stream throughput once, then pick a stream
    /// count from it. The default.
    #[default]
    Auto,
    /// Exactly this many concurrent streams; `1` is the plain sequential /
    /// async-prefetch reader (no scatter).
    Fixed(usize),
}

impl std::fmt::Display for ReadStreams {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Auto => f.write_str("auto"),
            Self::Fixed(n) => write!(f, "{n}"),
        }
    }
}

impl std::str::FromStr for ReadStreams {
    type Err = String;

    fn from_str(text: &str) -> Result<Self, Self::Err> {
        if text.eq_ignore_ascii_case("auto") {
            return Ok(Self::Auto);
        }
        match text.parse::<usize>() {
            Ok(0) => Err("--read-streams must be `auto` or at least 1".to_string()),
            Ok(n) => Ok(Self::Fixed(n)),
            Err(_) => Err(format!("expected `auto` or a positive number, got `{text}`")),
        }
    }
}

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
    /// This field is **honored** for the single-threaded readers built by
    /// [`create_bam_reader_with_opts`] and [`create_raw_bam_reader_with_opts`]:
    /// both route single-threaded input through fgumi-bgzf's decoder (see
    /// [`FgumiBgzfReader`]), which skips CRC32 verification when this is `false`
    /// (the decompressed-size check always runs). Multi-threaded input
    /// (`threads > 1`) still uses noodles' multithreaded reader, which always
    /// verifies; the multi-threaded unified pipeline instead gets its CRC policy
    /// from `PipelineConfig::verify_crc` (set in `build_pipeline_config` in
    /// `fgumi_lib`). Defaults to `true` (verify) — the safe, pre-existing
    /// behavior.
    pub verify_crc: bool,
    /// Concurrent positional-read policy for seekable regular-file inputs.
    ///
    /// `Fixed(1)` (the default) is the plain sequential / async-prefetch reader
    /// that every command except `fgumi sort` uses — only `sort` sets this from
    /// its own `--read-streams` flag. A higher count (or `Auto`) selects the
    /// [`crate::scatter_reader::ScatterReader`] for seekable files; non-seekable
    /// inputs (stdin, pipes) fall back to the sequential/async reader regardless.
    pub read_streams: ReadStreams,
}

impl Default for PipelineReaderOpts {
    fn default() -> Self {
        // NOTE: `Fixed(1)`, NOT `ReadStreams::default()` (which is `Auto`). Every
        // non-sort command relies on `..Default::default()` keeping today's
        // sequential/async behavior; defaulting to `Auto` here would silently
        // turn on scatter-read probing for the whole codebase.
        Self { async_reader: false, verify_crc: true, read_streams: ReadStreams::Fixed(1) }
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

        build_file_reader(file, path, &opts, label)?
    };

    crate::sam_input::normalize_to_bgzf(opened, path)
}

/// Choose the reader for an open regular file: the concurrent
/// [`crate::scatter_reader::ScatterReader`] when `--read-streams` asks for more
/// than one stream and the file is seekable (Unix only), otherwise the plain
/// sequential / async-prefetch reader. The scatter-vs-fallback decision (and its
/// "requested but unavailable" warning) is shared with `read_bam`'s entry point
/// via [`crate::scatter_reader::decide_reader`] so the two cannot drift.
fn build_file_reader(
    file: File,
    path: &Path,
    opts: &PipelineReaderOpts,
    label: &str,
) -> Result<Box<dyn Read + Send>> {
    let file = match crate::scatter_reader::decide_reader(file, opts.read_streams, path, label)
        .with_context(|| format!("open scatter reader for {}", path.display()))?
    {
        crate::scatter_reader::ScatterDecision::Scatter(scatter) => return Ok(scatter),
        crate::scatter_reader::ScatterDecision::Fallback(file) => file,
    };

    Ok(if opts.async_reader {
        log::info!("async {label} enabled: spawning fgumi-prefetch thread for {}", path.display());
        Box::new(crate::prefetch_reader::PrefetchReader::from_file(file))
    } else {
        Box::new(file)
    })
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
/// path. The async branch is buffered too: a `PrefetchReader` hands out large
/// chunks only when the caller asks for them, but the frame reader asks for 18
/// bytes then a block body, so an unbuffered async source pays the same per-read
/// overhead (plus a cross-thread handoff each time). The `PrefetchReader` stays
/// inside the buffer so its readahead still runs.
fn open_bgzf_reader(
    path: &Path,
    opts: PipelineReaderOpts,
    threads: usize,
    label: &str,
) -> Result<BgzfReaderEnum> {
    let normalized = open_normalized_with_opts(path, opts, label)?;
    let buffered: Box<dyn Read + Send> =
        Box::new(BufReader::with_capacity(BGZF_INPUT_BUFFER_SIZE, normalized));

    // Single-threaded input decodes through fgumi-bgzf, which honors
    // `opts.verify_crc` (noodles always verifies and has no skip knob). A
    // threaded fgumi-bgzf raw reader is out of scope, so `threads > 1` keeps
    // noodles' multithreaded decoder — its CRC policy is fixed on.
    if threads > 1 {
        Ok(make_bgzf_reader(buffered, threads))
    } else {
        Ok(BgzfReaderEnum::Fgumi(FgumiBgzfReader::new(buffered, opts.verify_crc)))
    }
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
    finish_raw_bam_reader(bgzf_reader, &path_ref.display().to_string())
}

/// Build a raw BAM reader over an already-opened, still-compressed BGZF byte
/// stream, decoding through fgumi-bgzf so `opts.verify_crc` is honored.
///
/// The single-threaded `correct` and `group` fast paths open their input exactly
/// once (required for stdin, which cannot be re-opened), parse the header for the
/// output, then must skip the header again to reach records. That second decode
/// is what this factory provides: routing it through [`FgumiBgzfReader`] — rather
/// than noodles' always-verify BGZF reader — makes `--no-check-crc` take effect
/// on those fast paths too (#800), matching the `--threads N` unified pipeline
/// (which already honors `PipelineConfig::verify_crc`).
///
/// The input is decoded single-threaded: this is the fast path taken precisely
/// when the caller passed no `--threads`, so no worker pool is warranted.
///
/// The returned header is the stream's own header, re-parsed here; callers that
/// have already augmented a header (e.g. synthesized @HD or added @PG) keep their
/// own copy and use this only to position the reader past the header bytes.
///
/// # Errors
///
/// Returns an error if the BAM header cannot be read from the stream.
pub fn create_raw_bam_reader_from_stream_with_opts(
    reader: Box<dyn Read + Send>,
    opts: PipelineReaderOpts,
) -> Result<(RawBamReaderAuto, Header)> {
    // Buffer before the frame reader for the same reason as `open_bgzf_reader`:
    // the BGZF frame reader asks for an 18-byte header and then a block body, so
    // an unbuffered source pays a syscall per block (and worse on stdin).
    let buffered: Box<dyn Read + Send> =
        Box::new(BufReader::with_capacity(BGZF_INPUT_BUFFER_SIZE, reader));
    let bgzf_reader = BgzfReaderEnum::Fgumi(FgumiBgzfReader::new(buffered, opts.verify_crc));
    finish_raw_bam_reader(bgzf_reader, "BAM stream")
}

/// Parse the BAM header from `bgzf_reader` and wrap the positioned decoder in a
/// [`RawBamReader`]. Shared tail of [`create_raw_bam_reader_with_opts`] and
/// [`create_raw_bam_reader_from_stream_with_opts`]: noodles is used only to
/// consume the header, after which the underlying BGZF decoder — already the
/// fgumi-bgzf or multithreaded arm chosen by the caller — is handed to the raw
/// reader for record access. `source` names the input in the header-read error.
fn finish_raw_bam_reader(
    bgzf_reader: BgzfReaderEnum,
    source: &str,
) -> Result<(RawBamReaderAuto, Header)> {
    let mut noodles_reader = noodles::bam::io::Reader::from(bgzf_reader);
    let header = noodles_reader
        .read_header()
        .with_context(|| format!("Failed to read header from: {source}"))?;

    // Get back the BGZF reader (header has been consumed) and wrap it.
    let raw_reader = RawBamReader::new(noodles_reader.into_inner());

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
/// The header block's CRC32 is **always** verified: the parse runs through
/// noodles' [`BgzfReader`], which has no CRC-skip knob. A pipeline's
/// `verify_crc = false` (`--no-check-crc`) therefore applies to the record body
/// only — a corrupted *header* block is rejected regardless. Threading the
/// policy here would mean swapping in the block-batching fgumi-bgzf reader,
/// whose readahead over-consumes the [`TeeReader`] and breaks the exact-byte
/// replay this function relies on, so the header parse deliberately stays on
/// noodles.
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
    use rstest::rstest;
    use std::num::NonZeroUsize;
    use std::str::FromStr;
    use std::sync::Arc;
    use std::sync::atomic::{AtomicUsize, Ordering};
    use tempfile::NamedTempFile;

    #[rstest]
    #[case::auto_lower("auto", ReadStreams::Auto)]
    #[case::auto_mixed_case("Auto", ReadStreams::Auto)]
    #[case::one("1", ReadStreams::Fixed(1))]
    #[case::four("4", ReadStreams::Fixed(4))]
    fn read_streams_from_str_accepts_valid(#[case] text: &str, #[case] expected: ReadStreams) {
        assert_eq!(ReadStreams::from_str(text).expect("valid"), expected);
    }

    #[rstest]
    #[case::zero("0")]
    #[case::negative("-1")]
    #[case::garbage("many")]
    fn read_streams_from_str_rejects_invalid(#[case] text: &str) {
        assert!(ReadStreams::from_str(text).is_err(), "'{text}' must be rejected");
    }

    #[rstest]
    #[case::auto(ReadStreams::Auto, "auto")]
    #[case::one(ReadStreams::Fixed(1), "1")]
    #[case::four(ReadStreams::Fixed(4), "4")]
    fn read_streams_display_round_trips(#[case] value: ReadStreams, #[case] expected: &str) {
        assert_eq!(value.to_string(), expected);
        assert_eq!(ReadStreams::from_str(expected).expect("round-trip"), value);
    }

    #[test]
    fn read_streams_default_is_auto() {
        assert_eq!(ReadStreams::default(), ReadStreams::Auto);
    }

    #[test]
    fn pipeline_reader_opts_default_read_streams_is_fixed_one() {
        // Load-bearing: the `PipelineReaderOpts` Default deliberately overrides
        // the `ReadStreams` enum default (`Auto`) with `Fixed(1)`, because every
        // non-sort command builds its opts via `..Default::default()`. A refactor
        // that let this fall back to `Auto` would silently enable scatter-read
        // probing across the whole codebase, so pin it.
        assert_eq!(PipelineReaderOpts::default().read_streams, ReadStreams::Fixed(1));
    }

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

    // ========================================================================
    // fgumi-bgzf raw-reader unify (#800)
    // ========================================================================

    /// Build SAM text with `n` aligned records over a single reference. Enough
    /// records to span more than one BGZF block once transcoded, so the raw
    /// reader's block-boundary refill logic is exercised.
    fn many_record_sam_text(n: usize) -> String {
        use std::fmt::Write as _;
        let mut s = String::from("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:100000\n");
        for i in 0..n {
            let pos = (i % 90_000) + 1;
            writeln!(s, "r{i}\t0\tchr1\t{pos}\t60\t4M\t*\t0\t0\tACGT\tIIII").expect("format");
        }
        s
    }

    /// Transcode SAM text to a BGZF BAM file on disk. For a large `sam` this
    /// spans several BGZF blocks, which the corrupted-CRC and roundtrip tests
    /// rely on.
    fn write_bgzf_bam(path: &Path, sam: &str) {
        let source: Box<dyn Read + Send> = Box::new(io::Cursor::new(sam.as_bytes().to_vec()));
        let mut normalized = crate::sam_input::normalize_to_bgzf(source, Path::new("test.sam"))
            .expect("transcode SAM to BGZF");
        let mut bytes = Vec::new();
        normalized.read_to_end(&mut bytes).expect("read transcoded BGZF");
        std::fs::write(path, bytes).expect("write BGZF BAM");
    }

    /// Flip a byte in the last BGZF block's CRC32 footer, adapted from PR2's
    /// `corrupt_last_block_crc` dedup integration helper. Requires the file to
    /// span at least two blocks so the corrupted block comes *after* reader
    /// construction succeeds: `FgumiBgzfReader` applies `verify_crc` uniformly to
    /// every block, so corrupting block 0 (the header's block) with
    /// `verify_crc: true` would fail inside reader construction (via `?`), before
    /// the intended `read_record`/count assertion can run.
    fn corrupt_last_block_crc(path: &Path) {
        let mut bytes = std::fs::read(path).expect("read bam for corruption");
        let mut cursor: &[u8] = &bytes;
        let blocks =
            fgumi_bgzf::read_raw_blocks(&mut cursor, 100_000).expect("read bgzf blocks from bam");
        assert!(
            blocks.len() >= 2,
            "test input must span >= 2 BGZF blocks so the corrupted block isn't the header's; \
             got {} -- generate more records",
            blocks.len()
        );
        let offset: usize =
            blocks[..blocks.len() - 1].iter().map(fgumi_bgzf::RawBgzfBlock::len).sum();
        let last = blocks.last().expect("checked len >= 2 above");
        // `read_raw_blocks` drops every BGZF EOF marker, so summing the returned
        // (real) block lengths yields the last block's on-disk offset only when no
        // marker sits *between* real blocks. Guard that: everything past the last
        // framed block must be whole trailing EOF markers (a writer may emit more
        // than one). An intermediate marker would leave real data here instead and
        // shift `crc_off` onto an unrelated byte, which the `>= 2` guard cannot
        // detect.
        let eof = &fgumi_bgzf::BGZF_EOF;
        let tail = &bytes[offset + last.len()..];
        assert!(
            tail.len().is_multiple_of(eof.len())
                && tail.chunks_exact(eof.len()).all(|chunk| chunk == &eof[..]),
            "bytes after the last framed block must be only trailing BGZF EOF markers; \
             an intermediate marker would invalidate the CRC offset"
        );
        let crc_off = offset + last.len() - fgumi_bgzf::BGZF_FOOTER_SIZE;
        bytes[crc_off] ^= 0x01;
        std::fs::write(path, bytes).expect("write corrupted bam");
    }

    /// Flip a byte in the *first* BGZF block's CRC32 footer — the block that
    /// carries the BAM header. The pipeline reader parses the header through
    /// noodles (which always verifies), so this corruption is rejected even with
    /// `verify_crc: false`.
    fn corrupt_first_block_crc(path: &Path) {
        let mut bytes = std::fs::read(path).expect("read bam for corruption");
        let mut cursor: &[u8] = &bytes;
        let blocks =
            fgumi_bgzf::read_raw_blocks(&mut cursor, 100_000).expect("read bgzf blocks from bam");
        let first_len = blocks.first().expect("at least one BGZF block").len();
        let crc_off = first_len - fgumi_bgzf::BGZF_FOOTER_SIZE;
        bytes[crc_off] ^= 0x01;
        std::fs::write(path, bytes).expect("write corrupted bam");
    }

    /// The pipeline reader parses the BAM header through noodles' BGZF reader,
    /// which always verifies CRC32 and has no skip knob. So `verify_crc: false`
    /// does **not** suppress a corrupted *header* block — it applies to the
    /// record body only. This pins that documented limitation (see
    /// [`read_header_and_replay`]).
    #[test]
    fn test_pipeline_reader_always_verifies_the_header_block() -> Result<()> {
        let file = tempfile::Builder::new().suffix(".bam").tempfile()?;
        write_bgzf_bam(file.path(), &many_record_sam_text(8000));
        corrupt_first_block_crc(file.path());

        let opts = PipelineReaderOpts { verify_crc: false, ..PipelineReaderOpts::default() };
        assert!(
            create_bam_reader_for_pipeline_with_opts(file.path(), opts).is_err(),
            "the header block is always CRC-verified, even with verify_crc: false"
        );
        Ok(())
    }

    /// Drain every record from a raw reader, returning the count.
    fn count_raw_records(reader: &mut RawBamReaderAuto) -> io::Result<usize> {
        let mut record = fgumi_raw_bam::RawRecord::new();
        let mut count = 0usize;
        while reader.read_record(&mut record)? > 0 {
            count += 1;
        }
        Ok(count)
    }

    /// The raw reader must honor `verify_crc`: with verification on, a corrupted
    /// block errors; with it off, the same file reads clean. Against the noodles
    /// path this could not be expressed — noodles always verifies — so the
    /// `verify_crc: false` case is the red test for the fgumi-bgzf unify (#800).
    #[test]
    fn test_raw_reader_honors_verify_crc_on_corrupted_block() -> Result<()> {
        let file = tempfile::Builder::new().suffix(".bam").tempfile()?;
        write_bgzf_bam(file.path(), &many_record_sam_text(8000));
        corrupt_last_block_crc(file.path());

        // verify_crc: true -> reading the corrupted block must error.
        let opts_verify = PipelineReaderOpts { verify_crc: true, ..PipelineReaderOpts::default() };
        let (mut reader, _header) = create_raw_bam_reader_with_opts(file.path(), 1, opts_verify)?;
        assert!(
            count_raw_records(&mut reader).is_err(),
            "verify_crc: true must reject a corrupted BGZF CRC32"
        );

        // verify_crc: false -> the same file reads clean (RED on the noodles path).
        let opts_skip = PipelineReaderOpts { verify_crc: false, ..PipelineReaderOpts::default() };
        let (mut reader, _header) = create_raw_bam_reader_with_opts(file.path(), 1, opts_skip)?;
        let count = count_raw_records(&mut reader)
            .expect("verify_crc: false must accept a corrupted BGZF CRC32");
        assert_eq!(count, 8000, "all records read with verify_crc: false");
        Ok(())
    }

    /// The from-stream raw reader must honor `verify_crc` exactly as the
    /// from-path reader does. This is the factory the single-threaded `correct`
    /// and `group` fast paths use on their already-opened stream (opened once so
    /// stdin — which cannot be re-opened — works), so `--no-check-crc` must take
    /// effect there too (#800). Against the noodles path these commands used
    /// before, the `verify_crc: false` case could not pass — noodles always
    /// verifies — so it is the red test for routing them through fgumi-bgzf.
    #[test]
    fn test_raw_reader_from_stream_honors_verify_crc() -> Result<()> {
        let file = tempfile::Builder::new().suffix(".bam").tempfile()?;
        write_bgzf_bam(file.path(), &many_record_sam_text(8000));
        corrupt_last_block_crc(file.path());
        let bytes = std::fs::read(file.path())?;

        // verify_crc: true -> reading the corrupted block must error.
        let opts_verify = PipelineReaderOpts { verify_crc: true, ..PipelineReaderOpts::default() };
        let stream: Box<dyn Read + Send> = Box::new(io::Cursor::new(bytes.clone()));
        let (mut reader, _header) =
            create_raw_bam_reader_from_stream_with_opts(stream, opts_verify)?;
        assert!(
            count_raw_records(&mut reader).is_err(),
            "verify_crc: true must reject a corrupted BGZF CRC32"
        );

        // verify_crc: false -> the same bytes read clean (RED on the noodles path).
        let opts_skip = PipelineReaderOpts { verify_crc: false, ..PipelineReaderOpts::default() };
        let stream: Box<dyn Read + Send> = Box::new(io::Cursor::new(bytes));
        let (mut reader, _header) = create_raw_bam_reader_from_stream_with_opts(stream, opts_skip)?;
        let count = count_raw_records(&mut reader)
            .expect("verify_crc: false must accept a corrupted BGZF CRC32");
        assert_eq!(count, 8000, "all records read with verify_crc: false");
        Ok(())
    }

    /// A file truncated mid-block — cut inside the final data block, including
    /// inside its 18-byte BGZF header — must surface as an error, never a clean
    /// `Ok(0)` end-of-stream. With `verify_crc` off the decompressed-size check
    /// is the only guard left, so both CRC modes must reject the truncation
    /// (see the path instruction: "a truncated or malformed block must surface
    /// as an error, never as a silent end-of-stream").
    #[rstest]
    #[case::inside_block_header(5)]
    #[case::mid_block_header(12)]
    #[case::header_only_body_cut(18)]
    #[case::inside_compressed_body(40)]
    fn test_raw_reader_rejects_truncated_final_block(
        #[case] kept_last_block_bytes: usize,
        #[values(true, false)] verify_crc: bool,
    ) -> Result<()> {
        // Build a multi-block BGZF BAM, then cut the file inside its final data
        // block at the requested offset (dropping any trailing EOF markers too).
        let bytes = {
            let file = tempfile::Builder::new().suffix(".bam").tempfile()?;
            write_bgzf_bam(file.path(), &many_record_sam_text(8000));
            std::fs::read(file.path())?
        };
        let blocks = {
            let mut cursor: &[u8] = &bytes;
            fgumi_bgzf::read_raw_blocks(&mut cursor, 100_000)?
        };
        assert!(blocks.len() >= 2, "need >= 2 blocks, got {}", blocks.len());
        let last_block_start: usize =
            blocks[..blocks.len() - 1].iter().map(fgumi_bgzf::RawBgzfBlock::len).sum();
        let last_block_len = blocks.last().expect("checked len >= 2 above").len();
        assert!(
            last_block_len > kept_last_block_bytes,
            "truncation offset {kept_last_block_bytes} must fall inside the last block \
             (len {last_block_len})"
        );

        let file = tempfile::Builder::new().suffix(".bam").tempfile()?;
        std::fs::write(file.path(), &bytes[..last_block_start + kept_last_block_bytes])?;

        // The header lives in the intact first block, so reader construction may
        // succeed; the truncation must then surface while draining records.
        // Either failure point counts — neither may yield a clean EOF.
        let opts = PipelineReaderOpts { verify_crc, ..PipelineReaderOpts::default() };
        let rejected = match create_raw_bam_reader_with_opts(file.path(), 1, opts) {
            Err(_) => true,
            Ok((mut reader, _header)) => count_raw_records(&mut reader).is_err(),
        };
        assert!(
            rejected,
            "truncated final block (kept {kept_last_block_bytes} bytes, verify_crc={verify_crc}) \
             must error, not return a clean EOF"
        );
        Ok(())
    }

    /// A clean multi-block file round-trips through the raw reader with
    /// verification on: every written record is read back.
    #[test]
    fn test_raw_reader_roundtrip_multiblock() -> Result<()> {
        let file = tempfile::Builder::new().suffix(".bam").tempfile()?;
        write_bgzf_bam(file.path(), &many_record_sam_text(8000));

        let opts = PipelineReaderOpts { verify_crc: true, ..PipelineReaderOpts::default() };
        let (mut reader, header) = create_raw_bam_reader_with_opts(file.path(), 1, opts)?;
        assert_eq!(header.reference_sequences().len(), 1, "header @SQ lost");
        assert_eq!(count_raw_records(&mut reader)?, 8000, "record count mismatch");
        Ok(())
    }

    /// A BGZF EOF marker appearing between data blocks (as a multithreaded
    /// writer emits between per-thread segments, or `cat a.bam b.bam` does) must
    /// not truncate the stream: the reader has to skip the marker and keep
    /// reading the blocks after it. This is the exact shape that a naive
    /// one-block-at-a-time refill mishandled.
    #[test]
    fn test_raw_reader_skips_intermediate_eof_marker() -> Result<()> {
        // Transcode to BGZF bytes, then splice a BGZF EOF marker after the first
        // block (not the header's block; it stays intact and parseable).
        let bytes = {
            let file = tempfile::Builder::new().suffix(".bam").tempfile()?;
            write_bgzf_bam(file.path(), &many_record_sam_text(8000));
            std::fs::read(file.path())?
        };
        let blocks = {
            let mut cursor: &[u8] = &bytes;
            fgumi_bgzf::read_raw_blocks(&mut cursor, 100_000)?
        };
        assert!(blocks.len() >= 2, "need >= 2 blocks, got {}", blocks.len());
        let first_block_len = blocks[0].len();

        let mut spliced = Vec::with_capacity(bytes.len() + fgumi_bgzf::BGZF_EOF.len());
        spliced.extend_from_slice(&bytes[..first_block_len]);
        spliced.extend_from_slice(&fgumi_bgzf::BGZF_EOF); // stray EOF marker mid-stream
        spliced.extend_from_slice(&bytes[first_block_len..]);

        let file = tempfile::Builder::new().suffix(".bam").tempfile()?;
        std::fs::write(file.path(), &spliced)?;

        let opts = PipelineReaderOpts { verify_crc: true, ..PipelineReaderOpts::default() };
        let (mut reader, _header) = create_raw_bam_reader_with_opts(file.path(), 1, opts)?;
        assert_eq!(
            count_raw_records(&mut reader)?,
            8000,
            "intermediate EOF marker truncated the record stream"
        );
        Ok(())
    }

    /// Byte-at-a-time reads must reconstruct exactly the same decompressed
    /// stream as one big read, across block boundaries. This pins the partial
    /// -read / cursor logic of the fgumi-bgzf streaming decoder.
    #[test]
    fn test_fgumi_bgzf_reader_partial_reads_match_full_decode() -> Result<()> {
        let bgzf = {
            let src: Box<dyn Read + Send> =
                Box::new(io::Cursor::new(many_record_sam_text(8000).into_bytes()));
            let mut normalized = crate::sam_input::normalize_to_bgzf(src, Path::new("test.sam"))?;
            let mut bytes = Vec::new();
            normalized.read_to_end(&mut bytes)?;
            bytes
        };

        // Full decode in one shot.
        let mut full = Vec::new();
        FgumiBgzfReader::new(Box::new(io::Cursor::new(bgzf.clone())), true)
            .read_to_end(&mut full)?;
        assert!(!full.is_empty(), "decoded stream should be non-empty");

        // Byte-at-a-time decode must match exactly.
        let mut reader = FgumiBgzfReader::new(Box::new(io::Cursor::new(bgzf.clone())), true);
        let mut one_byte_at_a_time = Vec::new();
        let mut byte = [0u8; 1];
        loop {
            let n = reader.read(&mut byte)?;
            if n == 0 {
                break;
            }
            one_byte_at_a_time.extend_from_slice(&byte[..n]);
        }
        assert_eq!(one_byte_at_a_time, full, "byte-at-a-time decode diverged from full decode");

        // A small odd-sized buffer (spanning block boundaries) must also match.
        let mut reader = FgumiBgzfReader::new(Box::new(io::Cursor::new(bgzf)), true);
        let mut small_chunks = Vec::new();
        let mut buf = [0u8; 7];
        loop {
            let n = reader.read(&mut buf)?;
            if n == 0 {
                break;
            }
            small_chunks.extend_from_slice(&buf[..n]);
        }
        assert_eq!(small_chunks, full, "7-byte-chunk decode diverged from full decode");
        Ok(())
    }

    /// A block decode failure must be **sticky**: once a block fails to decode,
    /// the reader must keep erroring rather than resume at the next frame and
    /// silently skip the corrupted block's records. Corrupt a *middle* block so
    /// valid blocks remain queued behind it — the exact shape where a
    /// non-sticky reader would hand back the following block's bytes.
    #[test]
    fn test_fgumi_bgzf_reader_poisons_after_decode_failure() -> Result<()> {
        let mut bgzf = {
            let src: Box<dyn Read + Send> =
                Box::new(io::Cursor::new(many_record_sam_text(8000).into_bytes()));
            let mut normalized = crate::sam_input::normalize_to_bgzf(src, Path::new("test.sam"))?;
            let mut bytes = Vec::new();
            normalized.read_to_end(&mut bytes)?;
            bytes
        };

        // Locate block boundaries and corrupt the CRC of block index 1 (a middle
        // block): blocks 0 and 2.. stay valid, so a non-sticky reader would skip
        // block 1 and continue serving block 2's bytes.
        let block_lens: Vec<usize> = {
            let mut cursor: &[u8] = &bgzf;
            fgumi_bgzf::read_raw_blocks(&mut cursor, 100_000)?
                .iter()
                .map(fgumi_bgzf::RawBgzfBlock::len)
                .collect()
        };
        assert!(block_lens.len() >= 3, "need >= 3 blocks, got {}", block_lens.len());
        let block1_offset = block_lens[0];
        let crc_off = block1_offset + block_lens[1] - fgumi_bgzf::BGZF_FOOTER_SIZE;
        bgzf[crc_off] ^= 0x01;

        let mut reader = FgumiBgzfReader::new(Box::new(io::Cursor::new(bgzf)), true);
        let mut buf = [0u8; 64];

        // Drain until the first error (must occur at the block 0 -> 1 boundary).
        loop {
            match reader.read(&mut buf) {
                Ok(0) => panic!("reader reached clean EOF without erroring on the corrupted block"),
                Ok(_) => {}
                Err(_) => break,
            }
        }

        // Sticky: every subsequent read must also error, never resume at block 2.
        for _ in 0..3 {
            assert!(
                reader.read(&mut buf).is_err(),
                "reader must stay poisoned after a decode failure, not skip the corrupted block"
            );
        }
        // `fill_buf` shares the poison too.
        assert!(reader.fill_buf().is_err(), "fill_buf must also honor the poison");
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

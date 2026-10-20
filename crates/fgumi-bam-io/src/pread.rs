//! Positional-read primitives shared by the sort's chain-native read paths.
//!
//! On a slow, deep-queue device (EBS gp3 is the motivating case) a single
//! outstanding `read()` leaves most of the device's bandwidth idle: gp3
//! sustains ~358 MB/s at queue-depth-1 but ~1177 MB/s with four reads in
//! flight. The sort raises the queue depth without any extra thread: a serial
//! planner step cuts each fill window into byte-range slices
//! ([`slice_ranges`], sized by [`ReadStreamsPolicy::slices_for`]), a `Parallel`
//! pool step reads each slice with one positional read
//! ([`read_at_exact`] over [`PositionalSource`]) into a recycled buffer from a
//! [`SliceBufferPool`], and the consumer frames the resulting [`SliceLease`]s.
//! Frames cut out of a slice borrow it ([`RawFrame::Borrowed`]); only a frame
//! that straddles two slices is copied ([`RawFrame::Owned`]).
//!
//! Positional reads take `&self` (POSIX `pread(2)` touches no shared file
//! offset), so concurrent reads on a shared `Arc<File>` need no
//! synchronization and no `unsafe`.
//!
//! **Unix only.** [`PositionalSource`] is implemented for [`std::fs::File`] on
//! Unix targets (fgumi supports Linux and macOS); the sort's spill reads depend
//! on it, so a non-Unix build fails to compile the sort rather than silently
//! taking another path.

use std::fmt;
use std::io;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::{Arc, Weak};
use std::time::Instant;

use crate::ReadStreams;

/// Bytes per fill window: one planner decision covers this much of a stream,
/// split into [`ReadStreamsPolicy::slices_for`] concurrent slices.
pub const FILL_BYTES: usize = 4 << 20;

/// Smallest slice a fill is split into. A window smaller than
/// `streams × MIN_SLICE_BYTES` uses fewer, larger slices: below ~512 KiB the
/// per-request overhead outweighs the queue-depth win.
const MIN_SLICE_BYTES: usize = 512 << 10;

/// Upper bound on the stream count, for both `Auto` and an explicit
/// `Fixed(n)`.
const MAX_STREAMS: usize = 8;

/// Fills per ratchet window: the ratchet decides once per this many fills
/// (32 MiB at [`FILL_BYTES`]).
const RATCHET_WINDOW_FILLS: u32 = 8;

/// The ratchet doubles the stream count when the device starved the framer
/// for more than this fraction of a window's wall time.
const RATCHET_STARVED_FRACTION: f64 = 0.25;

/// A byte-addressable source supporting concurrent positional reads.
///
/// Exists as a seam so the read steps can be tested against in-memory or
/// fault-injecting doubles (`test_sources`, behind the `test-utils` feature)
/// without a real file.
pub trait PositionalSource: Send + Sync {
    /// Read starting at `offset` into `buf`, returning the number of bytes read
    /// (may be short, as [`io::Read::read`]; `0` means EOF at `offset`).
    ///
    /// # Errors
    /// Propagates any I/O error from the underlying source.
    fn read_at(&self, buf: &mut [u8], offset: u64) -> io::Result<usize>;

    /// Total length of the source in bytes.
    ///
    /// # Errors
    /// Propagates any I/O error from querying the source's length.
    fn byte_len(&self) -> io::Result<u64>;

    /// Hint that `len` bytes from `offset` will be read soon, so the source
    /// can start fetching them before a later [`Self::read_at`] of the range.
    /// Advisory only: the read stays correct whatever the hint did, and waits
    /// on the device for any part not yet fetched. Allocation-free; the
    /// default does nothing.
    fn will_need(&self, offset: u64, len: u64) {
        let _ = (offset, len);
    }

    /// Whether a [`Self::read_at`] at `offset` would return without waiting
    /// for its source to produce bytes: a file always would (the device is
    /// slow, but it is not waiting on anyone), a stream only once it has bytes
    /// buffered, has ended or has failed ([`SequentialSource`]). A planner that
    /// issues a stream's reads to a shared pool asks this first, so no pool
    /// worker blocks on a producer that may itself be waiting on another of
    /// the pool's streams. The default is `true`.
    fn ready_at(&self, offset: u64) -> bool {
        let _ = offset;
        true
    }
}

#[cfg(unix)]
impl PositionalSource for std::fs::File {
    fn read_at(&self, buf: &mut [u8], offset: u64) -> io::Result<usize> {
        std::os::unix::fs::FileExt::read_at(self, buf, offset)
    }

    fn byte_len(&self) -> io::Result<u64> {
        Ok(self.metadata()?.len())
    }

    /// See [`crate::os_hints::advise_willneed_range`].
    fn will_need(&self, offset: u64, len: u64) {
        crate::os_hints::advise_willneed_range(self, offset, len);
    }
}

/// Read exactly `buf.len()` bytes at `offset`, looping over short reads and
/// retrying `Interrupted`. A `0`-byte read before the buffer is full means the
/// source ended earlier than its reported length promised (truncation or a
/// delete race), which fails closed with `UnexpectedEof`.
///
/// # Errors
/// `UnexpectedEof` on a short source; any other I/O error from the source.
pub fn read_at_exact<S: PositionalSource + ?Sized>(
    source: &S,
    buf: &mut [u8],
    offset: u64,
) -> io::Result<()> {
    let mut filled = 0usize;
    while filled < buf.len() {
        match source.read_at(&mut buf[filled..], offset + filled as u64) {
            Ok(0) => {
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    format!(
                        "positional read hit EOF at offset {} (wanted {} bytes)",
                        offset + filled as u64,
                        buf.len() - filled
                    ),
                ));
            }
            Ok(n) => filled += n,
            Err(ref e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    Ok(())
}

/// Bytes a [`SequentialSource`]'s pump buffers ahead of its reader before it
/// stops reading its stream. It reads once more only below this, so it holds
/// less than this plus one read ([`STREAM_PUMP_READ_BYTES`]).
///
/// The pump is what lets a coupled producer (`tee`, a demultiplexer: one
/// writer feeding several of a merge's inputs, blocking on whichever pipe is
/// full) make progress: every stream is drained into this buffer whatever the
/// merge reads next, so the producer is blocked only by a stream that is this
/// far, plus the merge's own read-ahead of it, ahead of the others in the
/// merge's order. 4 MiB is the merge's largest read of one stream (its hot
/// fill), so a read of an awaited stream that has caught up is one read.
pub(crate) const STREAM_PUMP_BYTES: usize = 4 << 20;

/// The most one `read` of a [`SequentialSource`]'s stream asks for, and the
/// capacity of each chunk the pump buffers: reads are appended to the last
/// chunk while it has room, so however small the producer's writes, the
/// buffer holds at most `STREAM_PUMP_BYTES / STREAM_PUMP_READ_BYTES + 2`
/// chunks. A pipe returns at most its buffer (64 KiB by default on Linux and
/// macOS) per read; a decompressing or transcoding stream returns what it has
/// decoded.
pub(crate) const STREAM_PUMP_READ_BYTES: usize = 256 << 10;

/// The most a [`SequentialSource`]'s pump allocates: its buffer of
/// [`STREAM_PUMP_BYTES`] plus one read, the unread part of a partly-read
/// front chunk and the unfilled part of the last, each under one chunk, and
/// its read buffer. 5 MiB.
#[cfg(test)]
pub(crate) const STREAM_PUMP_MAX_ALLOC_BYTES: usize =
    STREAM_PUMP_BYTES + 4 * STREAM_PUMP_READ_BYTES;

/// A non-seekable input (pipe, FIFO, process substitution) behind the
/// positional interface.
///
/// Its stream is read by a pump: a thread of its own (`fgumi-stream-pump`)
/// that reads the stream as fast as its producer writes it into a buffer of
/// about 4 MiB (it allocates at most 5 MiB), whatever the reader is doing. So
/// draining a pipe does not wait for the reader: a producer feeding several of
/// a merge's inputs in lockstep keeps writing to each while the merge reads
/// another, until one stream is that buffer ahead of the merge.
/// [`PositionalSource::ready_at`] reports whether the buffer has bytes (or the
/// stream has ended or failed), so a reader on a shared pool can avoid waiting
/// on it. The pump does not wake anyone when bytes arrive: such a reader
/// finds them when it next polls.
///
/// Reads must arrive at the stream's cursor (the read planner keeps one read
/// in flight per sequential slot), except reads below the header prefix, which
/// replay the bytes the header parse consumed through [`Self::header_reader`]
/// — so a reader of the source sees the stream from byte 0, as it would a
/// file. A read at the cursor returns what the buffer holds (up to the read's
/// length), waiting only while it holds nothing; `0` is the stream's EOF.
///
/// Dropping the source stops the pump once its current read returns: a pump
/// blocked reading a stream whose producer never writes again or closes keeps
/// the stream open until the process exits, and is never joined.
pub struct SequentialSource {
    path: std::path::PathBuf,
    state: parking_lot::Mutex<SeqState>,
    /// Lock-free mirror of `state.prefix.len()`, for [`PositionalSource::ready_at`],
    /// which must not wait behind a read holding `state`.
    prefix_len: AtomicU64,
    pump: Arc<Pump>,
}

/// The reader's side of a [`SequentialSource`].
struct SeqState {
    /// Stream bytes consumed so far (the header prefix included).
    cursor: u64,
    /// The bytes [`SequentialSource::header_reader`] returned, from offset 0.
    prefix: Vec<u8>,
}

/// The buffer a [`SequentialSource`]'s pump fills and its reader drains.
struct Pump {
    buf: parking_lot::Mutex<PumpBuf>,
    /// Signalled when bytes arrive or the stream ends (the reader waits on it).
    filled: parking_lot::Condvar,
    /// Signalled when the reader takes bytes or drops the source (the pump
    /// waits on it while the buffer is full).
    drained: parking_lot::Condvar,
}

/// How a pumped stream ended.
enum PumpEnd {
    /// Still being read.
    Open,
    /// At EOF.
    Eof,
    /// Failed; every later read reports this error (after the bytes before it).
    Failed(io::Error),
}

/// [`Pump`]'s state, under its lock.
struct PumpBuf {
    /// Bytes read from the stream and not yet taken, in order.
    chunks: std::collections::VecDeque<Vec<u8>>,
    /// Bytes of `chunks`' front chunk already taken.
    head: usize,
    /// Bytes in `chunks` not yet taken.
    len: usize,
    end: PumpEnd,
    /// The source was dropped: the pump stops.
    closed: bool,
}

impl Pump {
    /// Whether the reader can take bytes (or the stream's end) without waiting.
    fn ready(&self) -> bool {
        let b = self.buf.lock();
        b.len > 0 || !matches!(b.end, PumpEnd::Open)
    }

    /// Take up to `out.len()` buffered bytes, waiting only while none are
    /// buffered and the stream is open. `Ok(0)` is EOF.
    fn take(&self, out: &mut [u8]) -> io::Result<usize> {
        if out.is_empty() {
            return Ok(0);
        }
        let mut b = self.buf.lock();
        while b.len == 0 && matches!(b.end, PumpEnd::Open) {
            self.filled.wait(&mut b);
        }
        if b.len == 0 {
            return match &b.end {
                PumpEnd::Failed(e) => Err(io::Error::new(e.kind(), e.to_string())),
                PumpEnd::Eof | PumpEnd::Open => Ok(0),
            };
        }
        let mut n = 0;
        while n < out.len() {
            let head = b.head;
            let Some(front) = b.chunks.front() else { break };
            let k = (front.len() - head).min(out.len() - n);
            out[n..n + k].copy_from_slice(&front[head..head + k]);
            n += k;
            if head + k == front.len() {
                b.chunks.pop_front();
                b.head = 0;
            } else {
                b.head = head + k;
            }
        }
        b.len -= n;
        self.drained.notify_one();
        Ok(n)
    }

    /// Record how the stream ended and wake the reader.
    fn finish(&self, end: PumpEnd) {
        self.buf.lock().end = end;
        self.filled.notify_all();
    }

    /// The pump thread: open the stream (`open` runs here, so a stream that is
    /// classified or decoded as it is opened is not read on the caller's
    /// thread), then read it into the buffer until EOF, an error, or the
    /// source is dropped. A panic in the stream is the stream's error.
    fn run<F>(&self, open: F)
    where
        F: FnOnce() -> io::Result<Box<dyn io::Read + Send>>,
    {
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| self.pump(open)));
        let end = match result {
            Ok(Ok(true)) => PumpEnd::Eof,
            Ok(Ok(false)) => return,
            Ok(Err(e)) => PumpEnd::Failed(e),
            Err(payload) => PumpEnd::Failed(io::Error::other(crate::panic_message(
                payload.as_ref(),
                "the stream's reader panicked",
            ))),
        };
        self.finish(end);
    }

    /// [`Self::run`]'s loop: `Ok(true)` at EOF, `Ok(false)` once the source is
    /// dropped.
    fn pump<F>(&self, open: F) -> io::Result<bool>
    where
        F: FnOnce() -> io::Result<Box<dyn io::Read + Send>>,
    {
        let mut reader = open()?;
        let mut scratch = vec![0u8; STREAM_PUMP_READ_BYTES];
        loop {
            {
                let mut b = self.buf.lock();
                while b.len >= STREAM_PUMP_BYTES && !b.closed {
                    self.drained.wait(&mut b);
                }
                if b.closed {
                    return Ok(false);
                }
            }
            let n = match reader.read(&mut scratch) {
                Ok(0) => return Ok(true),
                Ok(n) => n,
                Err(e) if e.kind() == io::ErrorKind::Interrupted => continue,
                Err(e) => return Err(e),
            };
            let mut b = self.buf.lock();
            if b.closed {
                return Ok(false);
            }
            // Append to the last chunk while it has room, so small reads do
            // not each cost a chunk (see `STREAM_PUMP_READ_BYTES`).
            match b.chunks.back_mut() {
                Some(last) if last.capacity() - last.len() >= n => {
                    last.extend_from_slice(&scratch[..n]);
                }
                _ => {
                    let mut chunk = Vec::with_capacity(STREAM_PUMP_READ_BYTES);
                    chunk.extend_from_slice(&scratch[..n]);
                    b.chunks.push_back(chunk);
                }
            }
            b.len += n;
            drop(b);
            self.filled.notify_all();
        }
    }
}

impl SequentialSource {
    /// A source reading `reader`, a stream named `path` (for errors), through
    /// a pump thread started now.
    ///
    /// # Panics
    ///
    /// The pump thread cannot be spawned (see [`Self::open_normalized`], which
    /// returns that error).
    #[must_use]
    pub fn new(path: &std::path::Path, reader: Box<dyn io::Read + Send>) -> Self {
        Self::spawn(path, move || Ok(reader)).expect("failed to spawn a stream pump thread")
    }

    /// A source over the input at `path` (a pipe, FIFO, process substitution
    /// or stdin, or a SAM / plain-gzip file), normalized to BGZF as
    /// [`crate::open_normalized_input`] does. The input is opened, classified
    /// and decoded on the pump thread, which starts now, so this returns
    /// without waiting: a caller can open every stream of a producer that
    /// opens its FIFOs one at a time, in any order, and writes to none until
    /// all are open. An error opening or classifying the input (no such file,
    /// an empty or non-BAM/SAM stream) is the stream's error, which
    /// [`Self::stream_error`] reports and every read past the bytes before it
    /// returns.
    ///
    /// # Errors
    ///
    /// The pump thread cannot be spawned.
    pub fn open_normalized(path: &std::path::Path) -> io::Result<Self> {
        let name = path.to_path_buf();
        Self::spawn(path, move || {
            crate::reader::open_normalized_input(&name)
                .map_err(|e| io::Error::other(format!("{e:#}")))
        })
        .map_err(|e| {
            io::Error::new(
                e.kind(),
                format!("{}: cannot start a thread to read it: {e}", path.display()),
            )
        })
    }

    /// A source over the stream `open` returns, called on the pump thread,
    /// which starts now: whatever `open` reads (or waits for) to open,
    /// classify or decode the stream is read there, not on the caller's
    /// thread. An error from `open` is the stream's error.
    ///
    /// # Errors
    ///
    /// The pump thread cannot be spawned.
    pub(crate) fn spawn<F>(path: &std::path::Path, open: F) -> io::Result<Self>
    where
        F: FnOnce() -> io::Result<Box<dyn io::Read + Send>> + Send + 'static,
    {
        let pump = Arc::new(Pump {
            buf: parking_lot::Mutex::new(PumpBuf {
                chunks: std::collections::VecDeque::new(),
                head: 0,
                len: 0,
                end: PumpEnd::Open,
                closed: false,
            }),
            filled: parking_lot::Condvar::new(),
            drained: parking_lot::Condvar::new(),
        });
        let theirs = Arc::clone(&pump);
        std::thread::Builder::new()
            .name("fgumi-stream-pump".to_string())
            .spawn(move || theirs.run(open))?;
        Ok(Self {
            path: path.to_path_buf(),
            state: parking_lot::Mutex::new(SeqState { cursor: 0, prefix: Vec::new() }),
            prefix_len: AtomicU64::new(0),
            pump,
        })
    }

    /// Whether the pump's buffer is full and the stream still open: the pump
    /// has stopped reading, so the stream's producer blocks once the pipe's
    /// own buffer fills. For diagnosing a stall, not for flow control.
    #[must_use]
    pub fn is_backlogged(&self) -> bool {
        let b = self.pump.buf.lock();
        b.len >= STREAM_PUMP_BYTES && matches!(b.end, PumpEnd::Open)
    }

    /// Whether the pump has nothing for the reader and the stream is still
    /// open: a read at the cursor would wait for the producer.
    #[must_use]
    pub fn is_starved(&self) -> bool {
        !self.pump.ready()
    }

    /// The stream's error, once its pump has failed (opening, classifying,
    /// decoding or reading it, or a panic doing so): the error every read past
    /// the bytes before it returns. `None` while the stream is open or after a
    /// clean EOF. For a caller whose read failed inside a parser that wraps
    /// errors (the header parse), so it can report the stream's own error.
    #[must_use]
    pub fn stream_error(&self) -> Option<io::Error> {
        match &self.pump.buf.lock().end {
            PumpEnd::Failed(e) => Some(io::Error::new(e.kind(), e.to_string())),
            PumpEnd::Open | PumpEnd::Eof => None,
        }
    }

    /// A [`io::Read`] for parsing the stream's header: every byte it returns is
    /// kept in the prefix (replayed by [`PositionalSource::read_at`]) and
    /// advances the cursor.
    pub fn header_reader(&self) -> impl io::Read + '_ {
        HeaderReader(self)
    }

    /// Bytes the header parse consumed (the replayed prefix).
    #[must_use]
    pub fn prefix_len(&self) -> u64 {
        self.prefix_len.load(Ordering::Acquire)
    }
}

impl Drop for SequentialSource {
    /// Stop the pump: it exits once its current read returns (see the type's
    /// docs), dropping the stream.
    fn drop(&mut self) {
        self.pump.buf.lock().closed = true;
        self.pump.drained.notify_all();
    }
}

/// [`SequentialSource::header_reader`]'s reader.
struct HeaderReader<'a>(&'a SequentialSource);

impl io::Read for HeaderReader<'_> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let mut s = self.0.state.lock();
        let n = self.0.pump.take(buf)?;
        s.prefix.extend_from_slice(&buf[..n]);
        s.cursor += n as u64;
        self.0.prefix_len.store(s.prefix.len() as u64, Ordering::Release);
        Ok(n)
    }
}

impl PositionalSource for SequentialSource {
    /// `offset` below the prefix: a copy from it (possibly fewer bytes than
    /// `buf` holds); `offset` at the cursor: what the pump has buffered, up to
    /// `buf.len()`, waiting only while it has nothing (`0` is EOF); anything
    /// else: `InvalidInput` naming the path.
    fn read_at(&self, buf: &mut [u8], offset: u64) -> io::Result<usize> {
        let mut s = self.state.lock();
        let prefix_len = s.prefix.len() as u64;
        if offset < prefix_len {
            let at = usize::try_from(offset).expect("a prefix offset fits usize");
            let n = buf.len().min(s.prefix.len() - at);
            buf[..n].copy_from_slice(&s.prefix[at..at + n]);
            return Ok(n);
        }
        if offset != s.cursor {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "{}: out-of-order read at offset {offset}; the stream is at {}",
                    self.path.display(),
                    s.cursor
                ),
            ));
        }
        let n = self.pump.take(buf)?;
        s.cursor += n as u64;
        Ok(n)
    }

    /// `Unsupported`: a stream has no length (callers pass `len = None`).
    fn byte_len(&self) -> io::Result<u64> {
        Err(io::Error::new(
            io::ErrorKind::Unsupported,
            format!("{} is a stream and has no length", self.path.display()),
        ))
    }

    /// Below the prefix: always (a copy). Otherwise: whether the pump has
    /// bytes buffered, or the stream has ended or failed.
    fn ready_at(&self, offset: u64) -> bool {
        offset < self.prefix_len() || self.pump.ready()
    }
}

/// What [`open_seekable_bgzf`] found at a path.
pub enum SeekableInput {
    /// A regular BGZF file, opened (cursor at 0), and its length: readable
    /// positionally.
    Bgzf {
        /// The opened file.
        file: std::fs::File,
        /// Its length in bytes.
        len: u64,
    },
    /// An input that must be read sequentially, and why: `"is not a seekable
    /// regular file"` (stdin, a pipe, a FIFO, process substitution) or `"is not
    /// BGZF"` (SAM text, plain gzip).
    Sequential(&'static str),
}

/// Classify `path` as a positional (regular BGZF file) or sequential input:
/// the one rule every positional reader applies — `stat` first (never opening
/// a FIFO to sniff it), then a full format prefix read with
/// [`read_prefix`](crate::read_prefix) and [`classify_input`](crate::classify_input).
///
/// # Errors
///
/// Opening or reading a regular file failed.
pub fn open_seekable_bgzf(path: &std::path::Path) -> io::Result<SeekableInput> {
    use std::io::Seek as _;
    const NOT_SEEKABLE: &str = "is not a seekable regular file";
    if crate::is_stdin_path(path) || !std::fs::metadata(path)?.is_file() {
        return Ok(SeekableInput::Sequential(NOT_SEEKABLE));
    }
    let mut file = std::fs::File::open(path)?;
    let mut prefix = [0u8; crate::FORMAT_PREFIX_LEN];
    let n = crate::read_prefix(&mut file, &mut prefix)?;
    if crate::classify_input(&prefix[..n]) != crate::InputFormat::Bgzf {
        return Ok(SeekableInput::Sequential("is not BGZF"));
    }
    file.rewind()?;
    let len = file.metadata()?.len();
    Ok(SeekableInput::Bgzf { file, len })
}

/// Fill `buf` from `offset` by repeated `read_at`, stopping at a zero read
/// (EOF); `Interrupted` is retried. Returns the bytes read. The sequential twin
/// of [`read_at_exact`], which treats a short read as `UnexpectedEof`. A reader
/// sharing a thread pool with other streams should use [`read_once`] instead:
/// this waits until `buf` is full.
///
/// # Errors
/// Any I/O error from the source other than `Interrupted`.
#[deprecated(
    since = "0.8.1",
    note = "use `read_once`: a read that waits to fill its buffer can hold a pool worker \
            that a producer feeding several streams needs"
)]
pub fn read_until_eof<S: PositionalSource + ?Sized>(
    source: &S,
    buf: &mut [u8],
    offset: u64,
) -> io::Result<usize> {
    let mut filled = 0usize;
    while filled < buf.len() {
        match source.read_at(&mut buf[filled..], offset + filled as u64) {
            Ok(0) => break,
            Ok(n) => filled += n,
            Err(ref e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    Ok(filled)
}

/// One `read_at` of `buf` from `offset`, retrying only `Interrupted`. Returns
/// the bytes read, which may be fewer than `buf` holds — a stream returns what
/// it has buffered — and `0` only at EOF. The sequential twin of
/// [`read_at_exact`]: it does not wait to fill `buf`, because a reader waiting
/// for more of one stream holds its thread while the stream's producer may be
/// waiting on another stream (a producer feeding several inputs blocks on
/// whichever is full).
///
/// # Errors
/// Any I/O error from the source other than `Interrupted`.
pub fn read_once<S: PositionalSource + ?Sized>(
    source: &S,
    buf: &mut [u8],
    offset: u64,
) -> io::Result<usize> {
    loop {
        match source.read_at(buf, offset) {
            Err(ref e) if e.kind() == io::ErrorKind::Interrupted => {}
            r => return r,
        }
    }
}

/// One raise of the [`ReadStreamsPolicy`] ratchet, for `--sort-stats`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RatchetStep {
    /// Fills issued (over the whole run) when the raise happened.
    pub fill_index: u64,
    /// The stream count after the raise.
    pub streams: usize,
    /// Share of the deciding window's wall the framer was starved, in percent.
    pub starved_pct: u8,
}

/// Resolved `--read-streams`, shared by the input planner (which drives it)
/// and the spill planner (which only reads [`Self::streams`]).
///
/// `Auto` starts at one stream and doubles (to eight) at a window boundary
/// when the device starved the framer for more than a quarter of the
/// window's wall. It never decreases: too
/// few streams halves input throughput, while too many costs at most eight
/// syscalls per fill and no bytes, because the planner's ledger bounds
/// read-ahead. `Fixed(n)` pins the count and disables the ratchet.
pub struct ReadStreamsPolicy {
    streams: AtomicUsize,
    ratchet: Option<parking_lot::Mutex<Ratchet>>,
}

struct Ratchet {
    fills_in_window: u32,
    fills_total: u64,
    window_start: Instant,
    starved_ns: u64,
    last_visit: Option<Instant>,
    starved_at_last_visit: bool,
    history: Vec<RatchetStep>,
}

impl fmt::Debug for ReadStreamsPolicy {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("ReadStreamsPolicy")
            .field("streams", &self.streams())
            .field("auto", &self.is_auto())
            .finish_non_exhaustive()
    }
}

impl ReadStreamsPolicy {
    /// The `Auto` policy: one stream, ratchet armed.
    #[must_use]
    pub fn auto() -> Arc<Self> {
        Arc::new(Self {
            streams: AtomicUsize::new(1),
            ratchet: Some(parking_lot::Mutex::new(Ratchet {
                fills_in_window: 0,
                fills_total: 0,
                window_start: Instant::now(),
                starved_ns: 0,
                last_visit: None,
                starved_at_last_visit: false,
                history: Vec::new(),
            })),
        })
    }

    /// A pinned policy of `n` streams, clamped to `[1, 8]`.
    #[must_use]
    pub fn fixed(n: usize) -> Arc<Self> {
        Arc::new(Self { streams: AtomicUsize::new(n.clamp(1, MAX_STREAMS)), ratchet: None })
    }

    /// The policy the user's `--read-streams` value selects.
    #[must_use]
    pub fn from_flag(rs: ReadStreams) -> Arc<Self> {
        match rs {
            ReadStreams::Auto => Self::auto(),
            ReadStreams::Fixed(n) => Self::fixed(n),
        }
    }

    /// The current stream count.
    #[must_use]
    pub fn streams(&self) -> usize {
        self.streams.load(Ordering::Relaxed)
    }

    /// Whether the ratchet is armed (`Auto`).
    #[must_use]
    pub fn is_auto(&self) -> bool {
        self.ratchet.is_some()
    }

    /// Record one planner visit (see [`Self::observe_at`]) at the current time.
    pub fn observe(&self, starved_now: bool, fills_issued_this_visit: u32) {
        self.observe_at(Instant::now(), starved_now, fills_issued_this_visit);
    }

    /// Record one planner visit at `now`: integrate the starved wall since the
    /// previous visit (rectangle rule: the planner is visited every few µs, so
    /// quantisation is negligible), count the fills issued, and at a window
    /// boundary apply the doubling rule and start a new window. A no-op for a
    /// fixed policy.
    #[allow(
        clippy::cast_precision_loss,
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "nanosecond totals as f64 ratios; the percentage is clamped to 0..=100 before \
                  the u8 cast"
    )]
    pub fn observe_at(&self, now: Instant, starved_now: bool, fills_issued_this_visit: u32) {
        let Some(r) = &self.ratchet else { return };
        let mut r = r.lock();
        if let Some(last) = r.last_visit
            && r.starved_at_last_visit
        {
            r.starved_ns = r.starved_ns.saturating_add(
                u64::try_from(now.saturating_duration_since(last).as_nanos()).unwrap_or(u64::MAX),
            );
        }
        r.last_visit = Some(now);
        r.starved_at_last_visit = starved_now;
        r.fills_in_window += fills_issued_this_visit;
        r.fills_total += u64::from(fills_issued_this_visit);
        if r.fills_in_window < RATCHET_WINDOW_FILLS {
            return;
        }
        let wall = u64::try_from(now.saturating_duration_since(r.window_start).as_nanos())
            .unwrap_or(u64::MAX)
            .max(1);
        let share = r.starved_ns as f64 / wall as f64;
        let pct = (share * 100.0).clamp(0.0, 100.0) as u8;
        let cur = self.streams.load(Ordering::Relaxed);
        if share > RATCHET_STARVED_FRACTION && cur < MAX_STREAMS {
            let next = (cur * 2).min(MAX_STREAMS);
            self.streams.store(next, Ordering::Relaxed);
            let fill_index = r.fills_total;
            r.history.push(RatchetStep { fill_index, streams: next, starved_pct: pct });
            log::debug!(
                "read-streams ratchet: {cur} -> {next} at fill {fill_index} (starved {pct}%)"
            );
        }
        r.fills_in_window = 0;
        r.window_start = now;
        r.starved_ns = 0;
    }

    /// Slices to cut a `want`-byte fill into: the stream count, bounded by the
    /// `PreadSlices` clones that can run them at once (more slices than clones
    /// only queue) and by one slice per 512 KiB of `want`; at least one.
    #[must_use]
    pub fn slices_for(&self, want: usize, eligible_clones: usize) -> usize {
        self.streams().min(eligible_clones).min(want.div_ceil(MIN_SLICE_BYTES)).max(1)
    }

    /// Every raise so far, in order (empty for a fixed policy).
    #[must_use]
    pub fn history(&self) -> Vec<RatchetStep> {
        self.ratchet.as_ref().map(|r| r.lock().history.clone()).unwrap_or_default()
    }
}

/// The input planner's starvation signal: reads are outstanding up to (within
/// one fill of) the lookahead limit, so issuing more would exceed it, AND the
/// framer is waiting for the head slice. Out-of-order landings do not hide
/// starvation: `framer_waiting` reflects the head, not the byte count.
#[must_use]
pub fn starved_predicate(issued: u64, landed: u64, lookahead: u64, framer_waiting: bool) -> bool {
    issued > 0 && issued.saturating_sub(landed) + FILL_BYTES as u64 >= lookahead && framer_waiting
}

/// Split `[start, start + window)` into `slices` (at least one) contiguous
/// `(offset, len)` ranges in file order: equal shares, the last taking the
/// remainder.
///
/// # Panics
/// If one slice exceeds `u32::MAX` bytes (a fill is [`FILL_BYTES`]).
pub fn slice_ranges(start: u64, window: usize, slices: usize) -> impl Iterator<Item = (u64, u32)> {
    let slices = slices.max(1);
    let base = window / slices;
    (0..slices).map(move |i| {
        let off = start + (i * base) as u64;
        let len = if i + 1 == slices { window - i * base } else { base };
        (off, u32::try_from(len).expect("a read slice fits in u32"))
    })
}

/// Recycled slice buffers shared by the `PreadSlices` clones (which take
/// them) and every lease holder (whose last drop returns them).
///
/// A buffer is reused only for a slice of exactly its size, so every lease's
/// capacity equals its length: the bytes a reader charges for a slice (its
/// length) are the bytes the process holds for it, with no reuse slack for a
/// read-ahead bound to miss. Request sizes repeat (the fill and slice sizes
/// are fixed per path and per stream count), so exact reuse hits as often as a
/// looser fit would; a pool whose sizes change (the read-stream ratchet
/// shrinking slices) releases one mismatched buffer per miss and converges to
/// the sizes it is asked for. The released buffer is one of the rarest size in
/// the free list, so a one-off size (a spill run's short tail read) drains
/// first and the common size's buffers stay for the next hit. The gauges
/// charge capacity: `resident` is the capacity of live leases, `idle` the
/// capacity of the free list, and the peak is of their sum — the read path's
/// slice memory.
pub struct SliceBufferPool {
    free: parking_lot::Mutex<Vec<Vec<u8>>>,
    cap: usize,
    resident: AtomicU64,
    idle: AtomicU64,
    peak: AtomicU64,
}

impl fmt::Debug for SliceBufferPool {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("SliceBufferPool")
            .field("cap", &self.cap)
            .field("free", &self.free_len())
            .field("resident_bytes", &self.resident_bytes())
            .field("idle_bytes", &self.idle_bytes())
            .finish_non_exhaustive()
    }
}

/// Whether a pooled buffer of `capacity` may serve a `len`-byte slice: only
/// at exactly that size (see [`SliceBufferPool`]).
fn fits(capacity: usize, len: usize) -> bool {
    capacity == len
}

impl SliceBufferPool {
    /// A pool that keeps at most `cap` idle buffers.
    #[must_use]
    pub fn new(cap: usize) -> Arc<Self> {
        Arc::new(Self {
            free: parking_lot::Mutex::new(Vec::with_capacity(cap)),
            cap,
            resident: AtomicU64::new(0),
            idle: AtomicU64::new(0),
            peak: AtomicU64::new(0),
        })
    }

    /// A buffer of exactly `len` bytes and capacity: a pooled one of that
    /// capacity, else a fresh zeroed allocation (`calloc`, so the OS supplies zero pages
    /// rather than a memset). On a miss one pooled buffer of the rarest idle
    /// size is released, so a free list of the wrong sizes drains instead of
    /// pinning memory, while the common size keeps its buffers. Its
    /// contents are unspecified; a caller overwrites every byte before leasing
    /// it.
    #[must_use]
    pub fn take(&self, len: usize) -> Vec<u8> {
        let (hit, evicted) = {
            let mut free = self.free.lock();
            if let Some(i) = free.iter().position(|b| fits(b.capacity(), len)) {
                (Some(free.swap_remove(i)), None)
            } else {
                let count = |cap: usize| free.iter().filter(|b| b.capacity() == cap).count();
                let victim = (0..free.len()).min_by_key(|&i| count(free[i].capacity()));
                (None, victim.map(|i| free.swap_remove(i)))
            }
        };
        if let Some(e) = evicted {
            self.idle.fetch_sub(e.capacity() as u64, Ordering::Relaxed);
        }
        let Some(mut b) = hit else {
            return vec![0u8; len];
        };
        self.idle.fetch_sub(b.capacity() as u64, Ordering::Relaxed);
        // A buffer is returned at its full length; restore it if a holder
        // truncated it.
        b.resize(len, 0);
        b
    }

    /// Freeze a filled buffer into a shared read-only lease.
    #[must_use]
    pub fn lease(self: &Arc<Self>, bytes: Vec<u8>) -> SliceLease {
        let n = bytes.capacity() as u64;
        let now = self.resident.fetch_add(n, Ordering::Relaxed) + n;
        self.peak.fetch_max(now + self.idle.load(Ordering::Relaxed), Ordering::Relaxed);
        SliceLease(Arc::new(LeasedSlice { bytes, pool: Arc::downgrade(self) }))
    }

    fn checkin(&self, bytes: Vec<u8>) {
        let n = bytes.capacity() as u64;
        let mut free = self.free.lock();
        if free.len() < self.cap {
            // Idle before resident drops, so the sum never under-reads.
            self.idle.fetch_add(n, Ordering::Relaxed);
            free.push(bytes);
        }
        drop(free);
        self.resident.fetch_sub(n, Ordering::Relaxed);
    }

    /// Capacity of the buffers held by live leases.
    #[must_use]
    pub fn resident_bytes(&self) -> u64 {
        self.resident.load(Ordering::Relaxed)
    }

    /// Capacity of the idle buffers waiting for reuse.
    #[must_use]
    pub fn idle_bytes(&self) -> u64 {
        self.idle.load(Ordering::Relaxed)
    }

    /// The highest leased + idle capacity seen: the slice memory the pool held.
    #[must_use]
    pub fn peak_held_bytes(&self) -> u64 {
        self.peak.load(Ordering::Relaxed)
    }

    /// Idle buffers waiting for reuse.
    #[must_use]
    pub fn free_len(&self) -> usize {
        self.free.lock().len()
    }
}

struct LeasedSlice {
    bytes: Vec<u8>,
    pool: Weak<SliceBufferPool>,
}

impl Drop for LeasedSlice {
    fn drop(&mut self) {
        if let Some(p) = self.pool.upgrade() {
            p.checkin(std::mem::take(&mut self.bytes));
        }
    }
}

/// Refcounted read-only view of one read slice; `Clone` is an `Arc` clone and
/// the buffer returns to its pool when the last clone drops (on any thread).
/// It holds no back-pointer to a slot, so a stash → frame → lease chain cannot
/// form a cycle.
#[derive(Clone)]
pub struct SliceLease(Arc<LeasedSlice>);

impl SliceLease {
    /// Length of the slice in bytes.
    #[must_use]
    pub fn len(&self) -> usize {
        self.0.bytes.len()
    }

    /// Whether the slice is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.0.bytes.is_empty()
    }

    /// Capacity of the slice's buffer: what the lease holds and the pool's
    /// gauges charge. Equal to [`Self::len`] for a buffer from
    /// [`SliceBufferPool::take`] (exact-size reuse); a caller-built buffer
    /// passed to [`SliceBufferPool::lease`] may hold more.
    #[must_use]
    pub fn capacity(&self) -> usize {
        self.0.bytes.capacity()
    }
}

impl std::ops::Deref for SliceLease {
    type Target = [u8];
    fn deref(&self) -> &[u8] {
        &self.0.bytes
    }
}

impl fmt::Debug for SliceLease {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("SliceLease").field("len", &self.len()).finish()
    }
}

/// One raw frame (a BGZF block or a zstd frame body): borrowed from a leased
/// slice when wholly inside it, owned when it straddled slices.
///
/// `bytes::Bytes::from_owner(..).slice(..)` would give the same refcounted,
/// sub-sliceable buffer. This type stays explicit because the readers account
/// memory by what a frame pins: a borrowed frame holds its whole slice (which
/// the spill stash charges per slice, not per frame), an owned one only
/// itself, and the slice's return to its pool and the pool's gauges live on
/// the lease's drop. `Bytes` hides which of the two a frame is.
pub enum RawFrame {
    /// `lease[range]`; holding it keeps the whole slice out of its pool.
    Borrowed {
        /// The slice the frame lies in.
        lease: SliceLease,
        /// The frame's byte range within the slice.
        range: std::ops::Range<u32>,
    },
    /// An owned copy (a frame reassembled across slices, or any `Vec`).
    Owned(Vec<u8>),
}

impl RawFrame {
    /// The frame `lease[range]`, borrowing the slice.
    ///
    /// # Panics
    /// Panics if `range` is not within the slice. A slice's offsets fit `u32`:
    /// fills are cut into slices of at most `u32::MAX` bytes ([`slice_ranges`]).
    #[must_use]
    pub fn borrowed(lease: &SliceLease, range: std::ops::Range<usize>) -> Self {
        assert!(
            range.start <= range.end && range.end <= lease.len(),
            "frame range {range:?} outside a {}-byte slice",
            lease.len()
        );
        let to_u32 = |n: usize| u32::try_from(n).expect("a slice's offsets fit u32");
        Self::Borrowed { lease: lease.clone(), range: to_u32(range.start)..to_u32(range.end) }
    }

    /// The frame's bytes.
    #[inline]
    #[must_use]
    pub fn bytes(&self) -> &[u8] {
        match self {
            Self::Borrowed { lease, range } => &lease[range.start as usize..range.end as usize],
            Self::Owned(v) => v,
        }
    }

    /// Length of the frame in bytes.
    #[inline]
    #[must_use]
    pub fn len(&self) -> usize {
        self.bytes().len()
    }

    /// Whether the frame is empty.
    #[inline]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.bytes().is_empty()
    }
}

impl std::ops::Deref for RawFrame {
    type Target = [u8];
    fn deref(&self) -> &[u8] {
        self.bytes()
    }
}

impl From<Vec<u8>> for RawFrame {
    fn from(v: Vec<u8>) -> Self {
        RawFrame::Owned(v)
    }
}

impl fmt::Debug for RawFrame {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Borrowed { lease, range } => f
                .debug_struct("Borrowed")
                .field("slice_len", &lease.len())
                .field("range", range)
                .finish(),
            Self::Owned(v) => f.debug_struct("Owned").field("len", &v.len()).finish(),
        }
    }
}

/// In-memory and fault-injecting [`PositionalSource`] doubles for tests of
/// the read steps (here and in `fgumi-pipeline-io`).
#[cfg(any(test, feature = "test-utils"))]
pub mod test_sources {
    use super::PositionalSource;

    /// A reader returning at most 7 bytes per call, like a pipe with little
    /// buffered: partial reads are normal and must not be mistaken for EOF.
    pub struct Trickle(pub std::io::Cursor<Vec<u8>>);

    impl std::io::Read for Trickle {
        fn read(&mut self, b: &mut [u8]) -> io::Result<usize> {
            let n = b.len().min(7);
            self.0.read(&mut b[..n])
        }
    }
    use std::io;
    use std::time::Duration;

    fn copy_at(data: &[u8], buf: &mut [u8], offset: u64) -> usize {
        let Ok(off) = usize::try_from(offset) else { return 0 };
        if off >= data.len() {
            return 0;
        }
        let n = (data.len() - off).min(buf.len());
        buf[..n].copy_from_slice(&data[off..off + n]);
        n
    }

    /// In-memory source of known bytes.
    #[derive(Debug)]
    pub struct MemSource(Vec<u8>);

    impl MemSource {
        /// A source over `data`.
        #[must_use]
        pub fn new(data: Vec<u8>) -> Self {
            Self(data)
        }
    }

    impl PositionalSource for MemSource {
        fn read_at(&self, buf: &mut [u8], offset: u64) -> io::Result<usize> {
            Ok(copy_at(&self.0, buf, offset))
        }
        fn byte_len(&self) -> io::Result<u64> {
            Ok(self.0.len() as u64)
        }
    }

    /// A source that fails ("injected read failure") on the first read whose
    /// range covers `fail_at`.
    #[derive(Debug)]
    pub struct FaultySource {
        data: Vec<u8>,
        fail_at: u64,
    }

    impl FaultySource {
        /// A source over `data` that fails reads covering `fail_at`.
        #[must_use]
        pub fn new(data: Vec<u8>, fail_at: u64) -> Self {
            Self { data, fail_at }
        }
    }

    impl PositionalSource for FaultySource {
        fn read_at(&self, buf: &mut [u8], offset: u64) -> io::Result<usize> {
            let end = offset + buf.len() as u64;
            if offset <= self.fail_at && self.fail_at < end {
                return Err(io::Error::other("injected read failure"));
            }
            Ok(copy_at(&self.data, buf, offset))
        }
        fn byte_len(&self) -> io::Result<u64> {
            Ok(self.data.len() as u64)
        }
    }

    /// A source whose reported length (`claimed_len`) overstates its data, so
    /// a read past `data.len()` returns short — the truncation / delete-race
    /// case that must fail closed.
    #[derive(Debug)]
    pub struct ShortSource {
        data: Vec<u8>,
        claimed_len: u64,
    }

    impl ShortSource {
        /// A source over `data` that reports `claimed_len` bytes.
        #[must_use]
        pub fn new(data: Vec<u8>, claimed_len: u64) -> Self {
            Self { data, claimed_len }
        }
    }

    impl PositionalSource for ShortSource {
        fn read_at(&self, buf: &mut [u8], offset: u64) -> io::Result<usize> {
            Ok(copy_at(&self.data, buf, offset))
        }
        fn byte_len(&self) -> io::Result<u64> {
            Ok(self.claimed_len)
        }
    }

    /// In-memory source that sleeps `delay` per read (a slow device).
    #[derive(Debug)]
    pub struct TimedSource {
        data: Vec<u8>,
        delay: Duration,
    }

    impl TimedSource {
        /// A source over `data` whose every read takes at least `delay`.
        #[must_use]
        pub fn new(data: Vec<u8>, delay: Duration) -> Self {
            Self { data, delay }
        }
    }

    impl PositionalSource for TimedSource {
        fn read_at(&self, buf: &mut [u8], offset: u64) -> io::Result<usize> {
            if !self.delay.is_zero() {
                std::thread::sleep(self.delay);
            }
            Ok(copy_at(&self.data, buf, offset))
        }
        fn byte_len(&self) -> io::Result<u64> {
            Ok(self.data.len() as u64)
        }
    }
}

#[cfg(test)]
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "in-memory test doubles cast small, in-range offsets/lengths/bytes"
)]
mod tests {
    use super::test_sources::{FaultySource, MemSource, ShortSource};
    use super::*;
    use rstest::rstest;
    use std::io::Write;
    use std::time::Duration;
    use tempfile::NamedTempFile;

    // ---- slice_ranges ----

    #[rstest]
    #[case::one(0, 4 << 20, 1, vec![(0, 4 << 20)])]
    #[case::four_equal(
        100,
        4 << 20,
        4,
        vec![
            (100, 1 << 20),
            (100 + (1 << 20), 1 << 20),
            (100 + (2 << 20), 1 << 20),
            (100 + (3 << 20), 1 << 20),
        ]
    )]
    #[case::remainder_to_last(0, 10, 3, vec![(0, 3), (3, 3), (6, 4)])]
    fn slice_ranges_split_equally_last_takes_remainder(
        #[case] start: u64,
        #[case] window: usize,
        #[case] slices: usize,
        #[case] want: Vec<(u64, u32)>,
    ) {
        assert_eq!(slice_ranges(start, window, slices).collect::<Vec<_>>(), want);
    }

    // ---- slices_for ----

    #[rstest]
    #[case::fixed_four_clones_two(4, 4 << 20, 2, 2)]
    #[case::small_want_one_slice(8, 300 << 10, 16, 1)]
    #[case::two_min_slices(8, 1 << 20, 16, 2)]
    #[case::zero_clones_is_one(4, 4 << 20, 0, 1)]
    fn slices_never_exceed_eligible_clones(
        #[case] n: usize,
        #[case] want: usize,
        #[case] clones: usize,
        #[case] expect: usize,
    ) {
        assert_eq!(ReadStreamsPolicy::fixed(n).slices_for(want, clones), expect);
    }

    // ---- ratchet (deterministic clock via observe_at) ----

    /// Drive `windows` full windows; in each, the planner is visited every
    /// 100 µs and reports `starved` for the given share of visits. `now` is
    /// threaded through calls so consecutive drives continue one clock.
    fn drive(p: &ReadStreamsPolicy, now: &mut Instant, windows: u32, starved_share: f64) {
        let step = Duration::from_micros(100);
        let visits = 100u32;
        for _ in 0..windows {
            for v in 0..visits {
                let starved = f64::from(v) < starved_share * f64::from(visits);
                let fills = u32::from(v % (visits / RATCHET_WINDOW_FILLS) == 0);
                p.observe_at(*now, starved, fills);
                *now += step;
            }
        }
    }

    #[test]
    fn ratchet_stays_at_one_on_a_fast_source() {
        let p = ReadStreamsPolicy::auto();
        let mut now = Instant::now();
        drive(&p, &mut now, 6, 0.0);
        assert_eq!(p.streams(), 1);
    }

    #[test]
    fn ratchet_reaches_four_on_a_slow_source() {
        let p = ReadStreamsPolicy::auto();
        let mut now = Instant::now();
        drive(&p, &mut now, 2, 0.6);
        assert_eq!(p.streams(), 4, "two starved windows: 1 → 2 → 4");
        assert_eq!(p.history().len(), 2);
    }

    /// The discriminating case: the device keeps up but the framer is slow.
    /// The predicate — not the ratchet — must say "not starved", so the trace
    /// is fed through `starved_predicate` as the input planner does:
    /// outstanding reads are at the lookahead limit while the framer is busy
    /// (its input is non-empty, so `framer_waiting == false`).
    #[test]
    fn ratchet_stays_at_one_when_the_consumer_is_the_bottleneck() {
        let p = ReadStreamsPolicy::auto();
        let mut now = Instant::now();
        let lookahead = 8u64 << 20;
        for v in 0..800u32 {
            let issued = lookahead + (u64::from(v) << 10);
            let landed = issued - lookahead + FILL_BYTES as u64;
            let starved = starved_predicate(issued, landed, lookahead, false);
            p.observe_at(now, starved, u32::from(v % 12 == 0));
            now += Duration::from_micros(100);
        }
        assert_eq!(p.streams(), 1, "a slow consumer must not raise streams: {:?}", p.history());
    }

    #[test]
    fn ratchet_recovers_from_a_warm_head() {
        let p = ReadStreamsPolicy::auto();
        let mut now = Instant::now();
        drive(&p, &mut now, 1, 0.0); // the first 32 MiB at memcpy speed (page cache)
        drive(&p, &mut now, 3, 0.6); // the rest of the file is device-bound
        assert_eq!(p.streams(), MAX_STREAMS, "a monotone ratchet climbs after a warm head");
        assert!(
            p.history().iter().all(|h| h.starved_pct >= 25),
            "every raise came from the 25% rule: {:?}",
            p.history()
        );
    }

    #[test]
    fn ratchet_never_exceeds_max_or_decreases() {
        let p = ReadStreamsPolicy::auto();
        let mut now = Instant::now();
        drive(&p, &mut now, 10, 1.0);
        assert_eq!(p.streams(), MAX_STREAMS);
        drive(&p, &mut now, 10, 0.0);
        assert_eq!(p.streams(), MAX_STREAMS, "never decreases");
    }

    #[test]
    fn fixed_disables_the_ratchet() {
        let p = ReadStreamsPolicy::fixed(2);
        let mut now = Instant::now();
        drive(&p, &mut now, 10, 1.0);
        assert_eq!(p.streams(), 2);
        assert!(!p.is_auto());
        assert!(p.history().is_empty());
        assert_eq!(ReadStreamsPolicy::fixed(99).streams(), MAX_STREAMS);
        assert_eq!(ReadStreamsPolicy::fixed(0).streams(), 1);
    }

    #[rstest]
    #[case::auto(ReadStreams::Auto, 1, true)]
    #[case::fixed_one(ReadStreams::Fixed(1), 1, false)]
    #[case::fixed_four(ReadStreams::Fixed(4), 4, false)]
    fn from_flag_maps_the_cli_value(
        #[case] rs: ReadStreams,
        #[case] streams: usize,
        #[case] auto: bool,
    ) {
        let p = ReadStreamsPolicy::from_flag(rs);
        assert_eq!((p.streams(), p.is_auto()), (streams, auto));
    }

    // ---- the starved predicate (table, no clock) ----

    #[rstest]
    #[case::device_behind_framer_waiting(16 << 20, 8 << 20, 8 << 20, true, true)]
    #[case::device_behind_framer_busy(16 << 20, 8 << 20, 8 << 20, false, false)]
    #[case::device_keeping_up(16 << 20, 15 << 20, 8 << 20, true, false)]
    // A later slice landed but the head did not, so the framer waits.
    #[case::out_of_order_head_missing(16 << 20, 9 << 20, 8 << 20, true, true)]
    #[case::nothing_issued(0, 0, 8 << 20, true, false)]
    fn starved_predicate_matches_the_ledger_states(
        #[case] issued: u64,
        #[case] landed: u64,
        #[case] lookahead: u64,
        #[case] framer_waiting: bool,
        #[case] want: bool,
    ) {
        assert_eq!(starved_predicate(issued, landed, lookahead, framer_waiting), want);
    }

    // ---- pool + lease ----

    #[test]
    fn lease_returns_the_buffer_on_last_drop_from_any_thread() {
        let pool = SliceBufferPool::new(4);
        let mut b = pool.take(1 << 20);
        b[0] = 7;
        let lease = pool.lease(b);
        assert_eq!(pool.resident_bytes(), 1 << 20);
        let clone = lease.clone();
        drop(lease);
        assert_eq!(pool.free_len(), 0, "a live clone keeps it out of the pool");
        std::thread::spawn(move || drop(clone)).join().unwrap();
        assert_eq!(pool.free_len(), 1);
        assert_eq!(pool.resident_bytes(), 0);
        assert_eq!(pool.idle_bytes(), 1 << 20, "the returned buffer is idle, still held");
        assert_eq!(pool.peak_held_bytes(), 1 << 20);
    }

    /// Reuse is exact-size: a 1 MiB buffer serves another 1 MiB slice, and
    /// a request of any other size releases it rather than reusing it (so a
    /// changing request size drains the free list), so a lease's capacity is
    /// always its length.
    #[test]
    fn pool_reuses_only_an_exact_size_and_drops_mismatched_buffers() {
        let pool = SliceBufferPool::new(4);
        drop(pool.lease(pool.take(1 << 20)));
        assert_eq!((pool.resident_bytes(), pool.idle_bytes()), (0, 1 << 20));
        let reused = pool.lease(pool.take(1 << 20));
        assert_eq!(pool.free_len(), 0, "the 1 MiB buffer was reused");
        assert_eq!(pool.resident_bytes(), 1 << 20);
        drop(reused);
        let smaller = pool.lease(pool.take(600 << 10));
        assert_eq!((smaller.len(), smaller.capacity()), (600 << 10, 600 << 10));
        assert_eq!(pool.idle_bytes(), 0, "the mismatched idle buffer was released");
        assert_eq!(pool.resident_bytes(), 600 << 10);
        assert_eq!(pool.peak_held_bytes(), 1 << 20);
        drop(smaller);
        assert_eq!((pool.resident_bytes(), pool.idle_bytes()), (0, 600 << 10));
    }

    /// A miss releases a buffer of the rarest idle size, not whichever came
    /// back last: three full-size buffers and one tail-sized one, then a miss
    /// for a third size, keep the full-size buffers for the next full read.
    /// Catches: evicting by position (`pop`), which here drops a full-size
    /// buffer and forces the next full read to allocate.
    #[test]
    fn a_miss_releases_the_rarest_idle_size() {
        let pool = SliceBufferPool::new(8);
        let full: Vec<_> = (0..3).map(|_| pool.lease(pool.take(1 << 20))).collect();
        let tail = pool.lease(pool.take(5_000));
        drop(tail);
        drop(full); // the full-size buffers return last
        assert_eq!(pool.free_len(), 4);
        drop(pool.take(7_000)); // a miss
        assert_eq!(pool.free_len(), 3);
        assert_eq!(pool.idle_bytes(), 3 << 20, "the tail buffer was released");
    }

    #[test]
    fn pool_cap_drops_excess_and_reuses_the_allocation() {
        let pool = SliceBufferPool::new(1);
        let a = pool.lease(pool.take(64));
        let ptr = a.as_ptr();
        let b = pool.lease(pool.take(64));
        drop(a);
        drop(b);
        assert_eq!(pool.free_len(), 1, "the second return exceeds the cap and is dropped");
        let c = pool.take(64);
        assert_eq!((c.len(), c.capacity()), (64, 64));
        assert_eq!(c.as_ptr(), ptr, "reused the pooled allocation");
    }

    #[test]
    fn lease_outliving_its_pool_drops_cleanly() {
        let pool = SliceBufferPool::new(1);
        let lease = pool.lease(vec![1, 2, 3]);
        drop(pool);
        assert_eq!(&*lease, &[1, 2, 3]);
        drop(lease); // Weak::upgrade fails → plain drop, no panic
    }

    #[test]
    fn raw_frame_arms_expose_the_same_bytes() {
        let pool = SliceBufferPool::new(1);
        let lease = pool.lease(b"abcdef".to_vec());
        let b = RawFrame::borrowed(&lease, 2..5);
        let o = RawFrame::from(b"cde".to_vec());
        assert_eq!(b.bytes(), b"cde");
        assert_eq!(&*o, b"cde");
        assert_eq!(b.len(), 3);
        assert!(!b.is_empty());
        assert!(matches!(b, RawFrame::Borrowed { range, .. } if range == (2..5)));
    }

    #[test]
    #[should_panic(expected = "outside a 6-byte slice")]
    fn borrowed_frame_outside_its_slice_panics() {
        let pool = SliceBufferPool::new(1);
        let lease = pool.lease(b"abcdef".to_vec());
        let _ = RawFrame::borrowed(&lease, 4..7);
    }

    /// A leased buffer's capacity is what the lease charges, even when a
    /// caller hands the pool a buffer with spare capacity.
    #[test]
    fn lease_capacity_is_the_allocation() {
        let pool = SliceBufferPool::new(1);
        let mut v = Vec::with_capacity(64);
        v.extend_from_slice(&[0u8; 40]);
        let lease = pool.lease(v);
        assert_eq!((lease.len(), lease.capacity()), (40, 64));
        assert_eq!(pool.resident_bytes(), 64);
    }

    // ---- PositionalSource File impl ----

    #[cfg(unix)]
    #[test]
    fn file_positional_source_reads_offsets() {
        let mut f = NamedTempFile::new().expect("temp");
        f.write_all(b"0123456789").expect("write");
        f.flush().expect("flush");
        let file = f.reopen().expect("reopen");
        assert_eq!(PositionalSource::byte_len(&file).expect("len"), 10);
        let mut buf = [0u8; 3];
        let n = file.read_at(&mut buf, 4).expect("read_at");
        assert_eq!(&buf[..n], b"456");
    }

    // ---- read_at_exact ----

    #[test]
    fn read_at_exact_reads_the_whole_range() {
        let data: Vec<u8> = (0..200).map(|i| (i % 251) as u8).collect();
        let src = MemSource::new(data.clone());
        let mut buf = vec![0u8; 50];
        read_at_exact(&src, &mut buf, 100).unwrap();
        assert_eq!(buf, data[100..150]);
    }

    #[test]
    fn truncated_source_fails_closed_with_unexpected_eof() {
        // A source that claims more bytes than it holds (truncation / delete race)
        // must error, not silently return a short-but-clean slice.
        let src = ShortSource::new(vec![1u8; 50], 200);
        let mut buf = vec![0u8; 100];
        let err = read_at_exact(&src, &mut buf, 0).expect_err("truncation must fail closed");
        assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
        assert_eq!(err.to_string(), "positional read hit EOF at offset 50 (wanted 50 bytes)");
    }

    #[test]
    fn read_at_exact_surfaces_source_errors() {
        let src = FaultySource::new(vec![0u8; 100], 40);
        let mut buf = vec![0u8; 16];
        read_at_exact(&src, &mut buf, 0).expect("a range below the fault reads");
        let err = read_at_exact(&src, &mut buf, 32).expect_err("the fault must surface");
        assert_eq!(err.to_string(), "injected read failure");
    }

    /// A source that returns `Interrupted` once, then reads normally.
    struct InterruptOnce {
        inner: MemSource,
        fired: std::sync::atomic::AtomicBool,
    }

    impl PositionalSource for InterruptOnce {
        fn read_at(&self, buf: &mut [u8], offset: u64) -> io::Result<usize> {
            if !self.fired.swap(true, Ordering::Relaxed) {
                return Err(io::Error::from(io::ErrorKind::Interrupted));
            }
            self.inner.read_at(buf, offset)
        }
        fn byte_len(&self) -> io::Result<u64> {
            self.inner.byte_len()
        }
    }

    #[test]
    fn read_at_exact_retries_interrupted() {
        let src = InterruptOnce {
            inner: MemSource::new(b"abcdef".to_vec()),
            fired: std::sync::atomic::AtomicBool::new(false),
        };
        let mut buf = [0u8; 4];
        read_at_exact(&src, &mut buf, 1).unwrap();
        assert_eq!(&buf, b"bcde");
    }

    // ---- SequentialSource ----

    fn stream(n: usize) -> Vec<u8> {
        (0..n).map(|i| (i % 251) as u8).collect()
    }

    /// The header bytes `merge_headers` consumed are replayed from the prefix;
    /// reads past them continue the stream; EOF is a zero-byte read.
    #[test]
    fn sequential_source_serves_the_header_prefix_then_the_stream() {
        use super::test_sources::Trickle;
        use std::io::Read;
        let data = stream(10_000);
        let src = SequentialSource::new(
            std::path::Path::new("<pipe>"),
            Box::new(Trickle(std::io::Cursor::new(data.clone()))),
        );
        let mut head = vec![0u8; 700];
        src.header_reader().read_exact(&mut head).unwrap(); // what merge_headers does
        assert_eq!(head, data[..700]);
        let mut out = vec![0u8; data.len()];
        let mut at = 0;
        while at < out.len() {
            let n = src.read_at(&mut out[at..], at as u64).unwrap();
            assert!(n > 0, "a zero read before the end at {at}");
            at += n;
        }
        assert_eq!(out, data);
        assert_eq!(
            src.read_at(&mut [0u8; 16], data.len() as u64).unwrap(),
            0,
            "EOF is a zero read"
        );
    }

    /// R9: the error names the input; it is `InvalidInput`, never a wait.
    #[test]
    fn sequential_source_rejects_an_out_of_order_offset() {
        let src = SequentialSource::new(
            std::path::Path::new("/dev/fd/63"),
            Box::new(std::io::Cursor::new(stream(100))),
        );
        assert_eq!(src.read_at(&mut [0u8; 10], 0).unwrap(), 10);
        let err = src.read_at(&mut [0u8; 10], 50).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidInput);
        assert_eq!(
            err.to_string(),
            "/dev/fd/63: out-of-order read at offset 50; the stream is at 10"
        );
    }

    /// `read_until_eof`: loops partial reads, stops at a zero read, returns
    /// the bytes read.
    #[test]
    #[allow(deprecated)]
    fn read_until_eof_loops_partial_reads_and_stops_at_a_zero_read() {
        use super::test_sources::Trickle;
        let src = SequentialSource::new(
            std::path::Path::new("p"),
            Box::new(Trickle(std::io::Cursor::new(stream(100)))),
        );
        let mut buf = vec![0u8; 64];
        assert_eq!(read_until_eof(&src, &mut buf, 0).unwrap(), 64);
        assert_eq!(read_until_eof(&src, &mut buf, 64).unwrap(), 36);
        assert_eq!(read_until_eof(&src, &mut buf, 100).unwrap(), 0);
    }

    /// `read_once` (the sequential twin of `read_at_exact`): one read per
    /// call, which may be partial, and a zero read only at EOF; `Interrupted`
    /// is retried.
    #[test]
    fn read_once_returns_a_partial_read_and_zero_only_at_eof() {
        use super::test_sources::Trickle;
        let data = stream(100);
        let src = SequentialSource::new(
            std::path::Path::new("p"),
            Box::new(Trickle(std::io::Cursor::new(data.clone()))),
        );
        let mut buf = vec![0u8; 64];
        let mut got = Vec::new();
        while got.len() < data.len() {
            let n = read_once(&src, &mut buf, got.len() as u64).unwrap();
            assert!(n > 0, "a zero read before EOF at {}", got.len());
            got.extend_from_slice(&buf[..n]);
        }
        assert_eq!(got, data);
        assert_eq!(read_once(&src, &mut buf, 100).unwrap(), 0);
        let src = InterruptOnce {
            inner: MemSource::new(b"abcdef".to_vec()),
            fired: std::sync::atomic::AtomicBool::new(false),
        };
        assert_eq!(read_once(&src, &mut buf, 1).unwrap(), 5);
    }

    /// A reader that hands out `data` and then blocks until `release` is
    /// sent (or dropped), as a pipe whose producer has stalled does.
    struct Stalls {
        data: std::io::Cursor<Vec<u8>>,
        release: std::sync::mpsc::Receiver<()>,
    }

    impl std::io::Read for Stalls {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            let n = self.data.read(buf)?;
            if n == 0 {
                let _ = self.release.recv();
            }
            Ok(n)
        }
    }

    /// A reader that counts the bytes it hands out into its counter.
    struct Counting(std::io::Cursor<Vec<u8>>, Arc<AtomicU64>);

    impl std::io::Read for Counting {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            let n = self.0.read(buf)?;
            self.1.fetch_add(n as u64, Ordering::Relaxed);
            Ok(n)
        }
    }

    /// A reader that fails where its data ends.
    struct FailsAfter(std::io::Cursor<Vec<u8>>);

    impl std::io::Read for FailsAfter {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            match self.0.read(buf)? {
                0 => Err(io::Error::other("upstream failed")),
                n => Ok(n),
            }
        }
    }

    /// A reader that panics.
    struct Panics;

    impl std::io::Read for Panics {
        fn read(&mut self, _: &mut [u8]) -> io::Result<usize> {
            panic!("decoder bug")
        }
    }

    /// A reader that never ends, counting its reads; dropping it stores
    /// `u64::MAX` in the count.
    struct Endless(Arc<AtomicU64>);

    impl std::io::Read for Endless {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            self.0.fetch_add(1, Ordering::Relaxed);
            buf.fill(7);
            Ok(buf.len())
        }
    }

    impl Drop for Endless {
        fn drop(&mut self) {
            self.0.store(u64::MAX, Ordering::Relaxed);
        }
    }

    /// Spin (up to 10 s) until `cond` holds.
    fn eventually(cond: impl Fn() -> bool) {
        let deadline = std::time::Instant::now() + std::time::Duration::from_secs(10);
        while !cond() {
            assert!(std::time::Instant::now() < deadline, "condition never held");
            std::thread::sleep(std::time::Duration::from_millis(1));
        }
    }

    /// The pump drains the stream with no reader asking: a producer is never
    /// blocked by a reader that is busy elsewhere, until the buffer is full.
    /// It stops at `STREAM_PUMP_BYTES` (plus one read), and resumes as the
    /// reader takes bytes. Catches: a stream read only when its reader asks
    /// (the coupled-producer deadlock), and an unbounded pump.
    #[test]
    fn the_pump_drains_ahead_of_the_reader_up_to_its_bound() {
        let total = 3 * STREAM_PUMP_BYTES;
        let data = stream(total);
        let consumed = Arc::new(AtomicU64::new(0));
        let src = SequentialSource::new(
            std::path::Path::new("p"),
            Box::new(Counting(std::io::Cursor::new(data.clone()), Arc::clone(&consumed))),
        );
        let bound = (STREAM_PUMP_BYTES + STREAM_PUMP_READ_BYTES) as u64;
        // The pump stops once its buffer is full, and stays stopped.
        eventually(|| src.is_backlogged());
        let ahead = consumed.load(Ordering::Relaxed);
        assert!(ahead <= bound, "the pump read {ahead} bytes ahead, past its bound {bound}");
        for _ in 0..5 {
            std::thread::sleep(std::time::Duration::from_millis(10));
            assert_eq!(consumed.load(Ordering::Relaxed), ahead, "the pump read on while full");
        }
        let mut got = Vec::with_capacity(total);
        let mut buf = vec![0u8; 1 << 20];
        loop {
            let n = src.read_at(&mut buf, got.len() as u64).unwrap();
            if n == 0 {
                break;
            }
            got.extend_from_slice(&buf[..n]);
        }
        assert_eq!(got, data);
    }

    /// `ready_at`: true below the replayed prefix, and at the cursor once the
    /// pump holds bytes or the stream has ended, but not while the producer
    /// has written nothing new. A read at the cursor returns what is buffered
    /// without waiting for the rest of `buf`. Catches: a read issued to a pool
    /// worker that would block on a stalled producer.
    #[test]
    fn ready_at_reports_buffered_bytes_and_eof_and_a_read_takes_only_those() {
        use std::io::Read;
        let (release, rx) = std::sync::mpsc::channel();
        let src = SequentialSource::new(
            std::path::Path::new("p"),
            Box::new(Stalls { data: std::io::Cursor::new(stream(100)), release: rx }),
        );
        let mut head = [0u8; 10];
        src.header_reader().read_exact(&mut head).unwrap();
        assert!(src.ready_at(0), "the prefix is always ready");
        eventually(|| src.pump.buf.lock().len == 90);
        assert!(src.ready_at(10));
        let mut buf = vec![0u8; 4096];
        assert_eq!(src.read_at(&mut buf, 10).unwrap(), 90, "what is buffered, not a full read");
        assert_eq!(buf[..90], stream(100)[10..]);
        assert!(!src.ready_at(100), "nothing buffered and the stream still open");
        release.send(()).unwrap();
        eventually(|| src.ready_at(100));
        assert_eq!(src.read_at(&mut buf, 100).unwrap(), 0, "EOF");
    }

    /// A stream's error reaches the reader after the bytes before it, on every
    /// read from then on; an error opening the stream on the pump thread is
    /// the first read's; a panic in the stream is an error, not a hang.
    #[test]
    fn a_stream_error_follows_its_bytes_and_open_errors_and_panics_are_errors() {
        let src = SequentialSource::new(
            std::path::Path::new("p"),
            Box::new(FailsAfter(std::io::Cursor::new(stream(10)))),
        );
        let mut buf = [0u8; 64];
        let mut at = 0u64;
        let err = loop {
            match src.read_at(&mut buf, at) {
                Ok(n) => at += n as u64,
                Err(e) => break e,
            }
        };
        assert_eq!(at, 10, "every byte before the error");
        assert_eq!(err.to_string(), "upstream failed");
        assert_eq!(src.read_at(&mut buf, at).unwrap_err().to_string(), "upstream failed");
        assert!(src.ready_at(at), "a failed stream is ready (to report its error)");

        assert!(src.stream_error().is_some_and(|e| e.to_string() == "upstream failed"));

        let src = SequentialSource::spawn(std::path::Path::new("p"), || {
            Err(io::Error::new(io::ErrorKind::InvalidData, "not BAM"))
        })
        .unwrap();
        assert_eq!(src.read_at(&mut buf, 0).unwrap_err().kind(), io::ErrorKind::InvalidData);
        assert_eq!(src.stream_error().map(|e| e.kind()), Some(io::ErrorKind::InvalidData));
        let src = SequentialSource::new(std::path::Path::new("p"), Box::new(Panics));
        assert_eq!(src.read_at(&mut buf, 0).unwrap_err().to_string(), "decoder bug");
    }

    /// A producer writing a few bytes at a time does not cost the pump a chunk
    /// per write: reads are appended to the last chunk, so the buffer's
    /// allocation stays within `STREAM_PUMP_MAX_ALLOC_BYTES` once full.
    /// Catches: a byte bound that counts payload but not allocations (a
    /// million 1-byte `Vec`s for 1 MiB of tiny writes).
    #[test]
    fn small_writes_share_chunks_and_the_allocation_stays_bounded() {
        struct Tiny(u64);
        impl std::io::Read for Tiny {
            fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
                self.0 += 1;
                let n = buf.len().min(16);
                buf[..n].fill(3);
                Ok(n)
            }
        }
        let src = SequentialSource::new(std::path::Path::new("p"), Box::new(Tiny(0)));
        eventually(|| src.pump.buf.lock().len >= STREAM_PUMP_BYTES);
        let b = src.pump.buf.lock();
        assert!(
            b.chunks.len() <= STREAM_PUMP_BYTES / STREAM_PUMP_READ_BYTES + 2,
            "{} chunks for {} bytes of 16-byte reads",
            b.chunks.len(),
            b.len
        );
        let allocated: usize =
            b.chunks.iter().map(Vec::capacity).sum::<usize>() + STREAM_PUMP_READ_BYTES;
        assert!(allocated <= STREAM_PUMP_MAX_ALLOC_BYTES, "{allocated} bytes allocated");
    }

    /// `open_normalized` returns at once, whatever its input: the open (which
    /// for a FIFO waits for a writer) and the classification run on the pump
    /// thread, and their errors are the stream's. Catches: an open or format
    /// sniff on the caller's thread (a FIFO opened out of its producer's
    /// order deadlocks), and an open error lost or turned into a panic.
    #[cfg(unix)]
    #[test]
    fn open_normalized_opens_on_the_pump_thread_and_reports_open_errors() {
        let dir = tempfile::tempdir().unwrap();
        let fifo = dir.path().join("fifo");
        assert!(std::process::Command::new("mkfifo").arg(&fifo).status().unwrap().success());
        // No writer yet: opening it for reading would block this thread.
        let src = SequentialSource::open_normalized(&fifo).unwrap();
        assert!(!src.ready_at(0), "nothing written yet");
        let writer = std::thread::spawn({
            let fifo = fifo.clone();
            move || drop(std::fs::OpenOptions::new().write(true).open(fifo).unwrap())
        });
        writer.join().unwrap();
        eventually(|| src.ready_at(0));
        let err = src.read_at(&mut [0u8; 16], 0).unwrap_err();
        assert!(err.to_string().contains("No records read from"), "{err}");
        assert!(src.stream_error().is_some());

        let missing = dir.path().join("missing.bam");
        let src = SequentialSource::open_normalized(&missing).unwrap();
        eventually(|| src.ready_at(0));
        let err = src.read_at(&mut [0u8; 16], 0).unwrap_err();
        assert!(err.to_string().contains("missing.bam"), "{err}");
    }

    /// `is_backlogged` once the pump has filled its buffer from a stream still
    /// open, `is_starved` while it has nothing and the stream is open; neither
    /// once the stream has ended.
    #[test]
    fn backlogged_and_starved_report_the_pump_state() {
        let reads = Arc::new(AtomicU64::new(0));
        let full =
            SequentialSource::new(std::path::Path::new("p"), Box::new(Endless(Arc::clone(&reads))));
        eventually(|| full.is_backlogged());
        assert!(!full.is_starved());
        let (release, rx) = std::sync::mpsc::channel();
        let empty = SequentialSource::new(
            std::path::Path::new("q"),
            Box::new(Stalls { data: std::io::Cursor::new(Vec::new()), release: rx }),
        );
        assert!(empty.is_starved() && !empty.is_backlogged());
        drop(release);
        eventually(|| !empty.is_starved());
        assert!(!empty.is_backlogged(), "an ended stream is not backlogged");
    }

    /// Dropping the source stops a pump waiting for room. Catches: a pump
    /// thread (and its stream) kept alive by a full buffer forever.
    #[test]
    fn dropping_the_source_stops_a_full_pump() {
        let reads = Arc::new(AtomicU64::new(0));
        let src =
            SequentialSource::new(std::path::Path::new("p"), Box::new(Endless(Arc::clone(&reads))));
        eventually(|| src.pump.buf.lock().len >= STREAM_PUMP_BYTES);
        drop(src);
        eventually(|| reads.load(Ordering::Relaxed) == u64::MAX);
    }
}

#![deny(unsafe_code)]
//! Async userspace read-ahead adapter.
//!
//! [`PrefetchReader`] is a drop-in replacement for [`std::io::BufReader<File>`]
//! that performs asynchronous read-ahead on a dedicated OS thread. It exists to
//! decouple the blocking `read()` wait on a file descriptor from the pipeline
//! worker threads that consume the bytes.
//!
//! ## Motivation
//!
//! On Linux, the kernel's per-device read-ahead window (`read_ahead_kb`) is
//! 128 KB by default. A plain `BufReader<File>` is synchronous: when its
//! internal buffer drains, the next refill blocks in the kernel until pages
//! arrive from disk. During that stall the calling thread is parked and
//! cannot make progress on other work. In the fgumi unified pipeline, the
//! reader thread is also a pipeline worker, so a blocked read translates
//! directly into lost downstream throughput.
//!
//! `PrefetchReader` moves the blocking `read()` onto a dedicated producer
//! thread that pushes fixed-size chunks through a bounded
//! [`crossbeam_channel`]. Consumers see a normal [`std::io::Read`] interface;
//! internally, [`Read::read`] serves bytes out of the currently-held chunk and
//! only blocks when the producer has not yet delivered the next one. The
//! upshot is that stalls on the disk become independent of stalls in the
//! pipeline: the producer overlaps its disk wait with the consumer's CPU
//! work.
//!
//! This is essentially what the kernel's block-layer read-ahead does, but in
//! userspace — so it works without root, on any OS, and without having to
//! tune `/sys/block/*/queue/read_ahead_kb`.
//!
//! ## Lifecycle
//!
//! Constructing a `PrefetchReader` spawns exactly one OS thread, named
//! `fgumi-prefetch`. The thread owns the inner reader and exits when any of
//! the following happen:
//!
//! - The inner reader signals EOF (`Ok(0)` from `read`).
//! - The inner reader returns an error (the error is sent through the
//!   channel, then the thread exits).
//! - The [`PrefetchReader`] is dropped (the consumer-side receiver is
//!   destroyed, the producer's next `send` returns `Disconnected`, and the
//!   loop exits).
//!
//! `Drop` joins the producer thread, so leaks are impossible on well-behaved
//! inner readers. If the inner reader is currently parked in a long `read`
//! syscall, `Drop` will wait for it to return.

use std::fs::File;
use std::io::{self, Read};
use std::thread::{self, JoinHandle};

use crossbeam_channel::{Receiver, Sender, TryRecvError, bounded};

/// Default chunk size used by [`PrefetchReader::new`].
pub const DEFAULT_CHUNK_SIZE: usize = 4 * 1024 * 1024;

/// Default channel depth used by [`PrefetchReader::new`]. The producer will
/// keep up to this many filled chunks buffered ahead of the consumer.
pub const DEFAULT_PREFETCH_DEPTH: usize = 4;

/// How far ahead (in bytes) the producer thread asks the kernel to page in via
/// `posix_fadvise(POSIX_FADV_WILLNEED)` after each chunk fill. Only used when
/// the reader was constructed via [`PrefetchReader::from_file`], which captures
/// the raw fd needed for the syscall. 128 MiB is generous enough to cover
/// ~1 second of EBS gp3 baseline throughput (125 MiB/s), ensuring pages are
/// warm by the time the producer's `read()` reaches them.
///
/// Stored as `i64` to match `posix_fadvise`'s `off_t` signature directly.
const WILLNEED_LOOKAHEAD: i64 = 128 * 1024 * 1024;

/// An item handed from the producer thread to the consumer. Carries either a
/// filled chunk of bytes or a terminal I/O error.
type Item = io::Result<Vec<u8>>;

/// A `Read` adapter that performs asynchronous userspace prefetch on a
/// dedicated background thread.
///
/// See the [module docs](self) for the rationale and lifecycle model.
///
/// # Example
///
/// ```
/// use std::io::{Cursor, Read};
/// use fgumi_lib::prefetch_reader::PrefetchReader;
///
/// let data: Vec<u8> = (0..1024).map(|i| (i % 256) as u8).collect();
/// let mut reader = PrefetchReader::new(Cursor::new(data.clone()));
/// let mut out = Vec::new();
/// reader.read_to_end(&mut out).unwrap();
/// assert_eq!(out, data);
/// ```
#[derive(Debug)]
pub struct PrefetchReader {
    /// The chunk currently being served out to callers of `read`.
    ///
    /// Stored as `(data, position)`. `None` means we need to pull the next
    /// chunk from the channel on the next `read` call.
    current: Option<(Vec<u8>, usize)>,

    /// Receiving half of the prefetch channel. Held in `Option` so we can
    /// drop it explicitly in `Drop::drop` before joining the producer thread;
    /// dropping the receiver is what causes the producer's `send` to fail
    /// with `Disconnected`, which in turn terminates the producer loop.
    rx: Option<Receiver<Item>>,

    /// Join handle for the producer thread. `None` after join in `Drop`.
    handle: Option<JoinHandle<()>>,

    /// Number of times a consumer `read` call had to block on the channel
    /// because the producer had not yet delivered the next chunk. A high
    /// value relative to total reads suggests the prefetch depth is too
    /// shallow or the consumer is faster than the producer.
    consumer_stalls: u64,

    /// Total bytes handed back to callers of `read` so far.
    bytes_consumed: u64,
}

impl PrefetchReader {
    /// Construct a `PrefetchReader` from an arbitrary reader with default
    /// chunk size and prefetch depth.
    ///
    /// No `posix_fadvise(WILLNEED)` hints are issued because the inner reader
    /// may not be backed by a file descriptor. Use [`from_file`](Self::from_file)
    /// when wrapping a [`File`] to get kernel page-cache warming for free.
    #[must_use]
    pub fn new<R: Read + Send + 'static>(inner: R) -> Self {
        Self::with_config(inner, DEFAULT_CHUNK_SIZE, DEFAULT_PREFETCH_DEPTH)
    }

    /// Construct a `PrefetchReader` from an arbitrary reader with explicit
    /// `chunk_size` and `prefetch_depth`.
    ///
    /// Steady-state memory usage is bounded by `chunk_size * prefetch_depth`
    /// (plus one chunk being filled by the producer and one being consumed).
    ///
    /// # Panics
    ///
    /// Panics if `chunk_size == 0` or `prefetch_depth == 0`.
    #[must_use]
    pub fn with_config<R: Read + Send + 'static>(
        inner: R,
        chunk_size: usize,
        prefetch_depth: usize,
    ) -> Self {
        Self::build(inner, chunk_size, prefetch_depth, None)
    }

    /// Construct a `PrefetchReader` from a [`File`] with default chunk size
    /// and prefetch depth.
    ///
    /// On Linux the producer thread will call
    /// `posix_fadvise(POSIX_FADV_WILLNEED)` after each chunk fill to
    /// proactively page in the next [`WILLNEED_LOOKAHEAD`] bytes,
    /// making the reader independent of the kernel's default read-ahead
    /// window (`read_ahead_kb`). On non-Linux platforms this behaves
    /// identically to [`new`](Self::new).
    ///
    /// The file should be positioned at offset 0 when passed in — the
    /// WILLNEED hints assume reading starts from the beginning of the file.
    #[must_use]
    pub fn from_file(file: File) -> Self {
        Self::from_file_with_config(file, DEFAULT_CHUNK_SIZE, DEFAULT_PREFETCH_DEPTH)
    }

    /// Construct a `PrefetchReader` from a [`File`] with explicit `chunk_size`
    /// and `prefetch_depth`, plus kernel WILLNEED hints on Linux.
    ///
    /// # Panics
    ///
    /// Panics if `chunk_size == 0` or `prefetch_depth == 0`.
    #[must_use]
    pub fn from_file_with_config(file: File, chunk_size: usize, prefetch_depth: usize) -> Self {
        let hint_fd = crate::os_hints::hint_fd(&file);
        Self::build(file, chunk_size, prefetch_depth, hint_fd)
    }

    /// Shared construction logic.
    fn build<R: Read + Send + 'static>(
        inner: R,
        chunk_size: usize,
        prefetch_depth: usize,
        hint_fd: Option<i32>,
    ) -> Self {
        assert!(chunk_size > 0, "PrefetchReader chunk_size must be > 0");
        assert!(prefetch_depth > 0, "PrefetchReader prefetch_depth must be > 0");

        let (tx, rx) = bounded::<Item>(prefetch_depth);
        let handle = thread::Builder::new()
            .name("fgumi-prefetch".to_string())
            .spawn(move || producer_main(inner, chunk_size, hint_fd, &tx))
            .expect("failed to spawn fgumi-prefetch thread");

        Self {
            current: None,
            rx: Some(rx),
            handle: Some(handle),
            consumer_stalls: 0,
            bytes_consumed: 0,
        }
    }

    /// Total bytes served to callers of [`Read::read`] so far.
    #[must_use]
    pub fn bytes_consumed(&self) -> u64 {
        self.bytes_consumed
    }

    /// Number of times a `read` call had to block waiting for the producer
    /// to deliver the next chunk. Useful as a prototype-phase signal for
    /// whether [`DEFAULT_PREFETCH_DEPTH`] is large enough.
    #[must_use]
    pub fn consumer_stalls(&self) -> u64 {
        self.consumer_stalls
    }
}

/// Entry point for the producer thread. Wraps [`producer_loop`] so that a
/// panic inside the inner reader surfaces on the consumer side as an
/// [`io::Error`] instead of a silently-truncated stream.
fn producer_main<R: Read>(inner: R, chunk_size: usize, hint_fd: Option<i32>, tx: &Sender<Item>) {
    let tx_for_panic = tx.clone();
    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        producer_loop(inner, chunk_size, hint_fd, tx);
    }));
    if let Err(payload) = result {
        let msg = match payload.downcast_ref::<&'static str>() {
            Some(s) => (*s).to_string(),
            None => match payload.downcast_ref::<String>() {
                Some(s) => s.clone(),
                None => "fgumi-prefetch producer thread panicked".to_string(),
            },
        };
        let _ = tx_for_panic.send(Err(io::Error::other(msg)));
    }
}

/// The actual producer loop. Allocates a `chunk_size`-byte buffer, issues a
/// single `read()` call, and sends whatever bytes were returned through the
/// channel. Emitting after the first successful read (rather than looping to
/// fill the full chunk) ensures that short reads from pipes or throttled
/// readers are delivered promptly. For file-backed readers, the OS typically
/// returns a full page-aligned read in one call, so this rarely increases
/// channel traffic. Tolerates `Interrupted` errors. Exits on EOF, I/O error,
/// or when the consumer drops the receiver.
///
/// When `hint_fd` is `Some`, the loop calls
/// `posix_fadvise(POSIX_FADV_WILLNEED)` after each chunk to proactively page
/// in the next [`WILLNEED_LOOKAHEAD`] bytes. This removes the
/// dependence on the kernel's default `read_ahead_kb` setting.
fn producer_loop<R: Read>(
    mut inner: R,
    chunk_size: usize,
    hint_fd: Option<i32>,
    tx: &Sender<Item>,
) {
    let mut position: i64 = 0;

    loop {
        let mut buf = vec![0u8; chunk_size];
        let mut filled: usize = 0;
        let mut eof = false;

        loop {
            match inner.read(&mut buf[filled..]) {
                Ok(0) => {
                    eof = true;
                    break;
                }
                Ok(n) => {
                    filled += n;
                    break;
                }
                Err(e) if e.kind() == io::ErrorKind::Interrupted => (),
                Err(e) => {
                    // Flush any bytes we already managed to read *before*
                    // surfacing the error. Otherwise a mid-chunk failure
                    // would silently discard data that was successfully
                    // returned by the inner reader.
                    if filled > 0 {
                        buf.truncate(filled);
                        let _ = tx.send(Ok(buf));
                    }
                    let _ = tx.send(Err(e));
                    return;
                }
            }
        }

        if filled == 0 && eof {
            return;
        }

        position = position.saturating_add(i64::try_from(filled).unwrap_or(i64::MAX));

        // Ask the kernel to start paging in bytes ahead of our current
        // position. The call is non-blocking — the kernel initiates the I/O
        // and returns immediately.
        if let Some(fd) = hint_fd {
            crate::os_hints::advise_willneed_raw(fd, position, WILLNEED_LOOKAHEAD);
        }

        buf.truncate(filled);
        if tx.send(Ok(buf)).is_err() {
            // Consumer dropped the receiver; shutdown cleanly.
            return;
        }

        if eof {
            return;
        }
    }
}

impl Read for PrefetchReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if buf.is_empty() {
            return Ok(0);
        }

        loop {
            // Serve from the current chunk if any bytes remain.
            if let Some((data, pos)) = self.current.as_mut() {
                if *pos < data.len() {
                    let n = std::cmp::min(buf.len(), data.len() - *pos);
                    buf[..n].copy_from_slice(&data[*pos..*pos + n]);
                    *pos += n;
                    self.bytes_consumed += n as u64;
                    return Ok(n);
                }
                // Current chunk exhausted; drop it and pull the next one.
                self.current = None;
            }

            // No current chunk. Pull the next one from the producer.
            let Some(rx) = self.rx.as_ref() else {
                return Ok(0);
            };

            // Fast path: a chunk is already waiting.
            let item = match rx.try_recv() {
                Ok(item) => item,
                Err(TryRecvError::Disconnected) => return Ok(0),
                Err(TryRecvError::Empty) => {
                    // Slow path: producer hasn't delivered yet; we block.
                    self.consumer_stalls += 1;
                    match rx.recv() {
                        Ok(item) => item,
                        Err(_) => return Ok(0),
                    }
                }
            };

            match item {
                Ok(data) if !data.is_empty() => self.current = Some((data, 0)),
                Ok(_) => {} // defensive: skip empty chunks
                Err(e) => return Err(e),
            }
        }
    }
}

impl Drop for PrefetchReader {
    fn drop(&mut self) {
        // Drop the receiver first so the producer's next `send` returns
        // `Disconnected` and the loop exits. Then join the thread to
        // guarantee no leak.
        self.rx = None;
        self.current = None;
        if let Some(handle) = self.handle.take() {
            if handle.join().is_err() {
                log::debug!("fgumi-prefetch producer thread panicked during shutdown");
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;
    use std::io::Cursor;
    use std::sync::Arc;
    use std::sync::atomic::{AtomicBool, Ordering};

    /// Build a deterministic byte stream of the given length.
    fn sample_bytes(len: usize) -> Vec<u8> {
        (0..len).map(|i| u8::try_from(i % 251).expect("mod 251 fits in u8")).collect()
    }

    #[test]
    fn empty_input_returns_zero_immediately() {
        let mut reader = PrefetchReader::new(Cursor::new(Vec::<u8>::new()));
        let mut buf = [0u8; 16];
        assert_eq!(reader.read(&mut buf).unwrap(), 0);
        // Bytes counter unchanged.
        assert_eq!(reader.bytes_consumed(), 0);
    }

    #[test]
    fn from_file_reads_correctly() {
        use std::io::Write;
        let data = sample_bytes(50_000);
        let mut tmp = tempfile::NamedTempFile::new().expect("create temp file");
        tmp.write_all(&data).expect("write temp file");
        let file = File::open(tmp.path()).expect("reopen temp file");
        let mut reader = PrefetchReader::from_file(file);
        let mut out = Vec::new();
        reader.read_to_end(&mut out).unwrap();
        assert_eq!(out, data);
        assert_eq!(reader.bytes_consumed(), data.len() as u64);
    }

    #[test]
    fn read_to_end_small_matches_input() {
        let data = b"hello, fgumi prefetch".to_vec();
        let mut reader = PrefetchReader::new(Cursor::new(data.clone()));
        let mut out = Vec::new();
        reader.read_to_end(&mut out).unwrap();
        assert_eq!(out, data);
        assert_eq!(reader.bytes_consumed(), data.len() as u64);
    }

    #[test]
    fn read_to_end_large_matches_input() {
        // ~1 MiB — exercises multiple chunks through the channel.
        let data = sample_bytes(1_000_003);
        let mut reader = PrefetchReader::with_config(Cursor::new(data.clone()), 8 * 1024, 2);
        let mut out = Vec::new();
        reader.read_to_end(&mut out).unwrap();
        assert_eq!(out, data);
        assert_eq!(reader.bytes_consumed(), data.len() as u64);
    }

    #[test]
    fn tiny_chunk_size_and_many_small_reads() {
        let data = sample_bytes(5_000);
        let mut reader = PrefetchReader::with_config(Cursor::new(data.clone()), 17, 2);
        let mut out = Vec::new();
        let mut tmp = [0u8; 7];
        loop {
            let n = reader.read(&mut tmp).unwrap();
            if n == 0 {
                break;
            }
            out.extend_from_slice(&tmp[..n]);
        }
        assert_eq!(out, data);
    }

    #[test]
    fn repeated_read_after_eof_returns_zero_forever() {
        let mut reader = PrefetchReader::new(Cursor::new(b"abc".to_vec()));
        let mut out = Vec::new();
        reader.read_to_end(&mut out).unwrap();
        assert_eq!(out, b"abc");

        let mut tmp = [0u8; 8];
        for _ in 0..10 {
            assert_eq!(reader.read(&mut tmp).unwrap(), 0);
        }
    }

    #[test]
    fn drop_before_consuming_does_not_hang() {
        // 1 MB of data, tiny chunks, shallow channel: the producer will fill
        // the channel quickly and then block on send. When we drop the
        // reader, the send must unblock so `Drop::join` returns.
        let data = vec![0u8; 1_000_000];
        let _reader = PrefetchReader::with_config(Cursor::new(data), 4 * 1024, 2);
        // Drop happens at end of scope.
    }

    #[test]
    fn partial_read_then_drop_does_not_hang() {
        let data = sample_bytes(500_000);
        let mut reader = PrefetchReader::with_config(Cursor::new(data), 4 * 1024, 2);
        let mut tmp = [0u8; 32];
        // Read a bit, then drop without finishing.
        let n = reader.read(&mut tmp).unwrap();
        assert!(n > 0);
    }

    #[test]
    fn error_from_inner_reader_propagates_once() {
        struct AlwaysErr;
        impl Read for AlwaysErr {
            fn read(&mut self, _: &mut [u8]) -> io::Result<usize> {
                Err(io::Error::new(io::ErrorKind::PermissionDenied, "nope"))
            }
        }
        let mut reader = PrefetchReader::new(AlwaysErr);
        let mut buf = [0u8; 16];
        let err = reader.read(&mut buf).expect_err("first read should error");
        assert_eq!(err.kind(), io::ErrorKind::PermissionDenied);
        // After the error is propagated, subsequent reads see the producer
        // exit and return EOF.
        assert_eq!(reader.read(&mut buf).unwrap(), 0);
    }

    #[test]
    fn error_after_some_data_delivers_data_then_error() {
        // Reader that returns one chunk of data, then an error forever.
        struct DataThenErr {
            sent: bool,
        }
        impl Read for DataThenErr {
            fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
                if self.sent {
                    Err(io::Error::other("subsequent error"))
                } else {
                    self.sent = true;
                    let n = buf.len().min(128);
                    for (i, b) in buf.iter_mut().take(n).enumerate() {
                        *b = u8::try_from(i).expect("n <= 128 so i fits in u8");
                    }
                    Ok(n)
                }
            }
        }
        // chunk_size > 128 forces the producer to call read() again after
        // the first short read, which triggers the error.
        let mut reader = PrefetchReader::with_config(DataThenErr { sent: false }, 1024, 2);
        let mut out = Vec::new();
        let mut tmp = [0u8; 256];

        // First read should deliver the 128 bytes of data.
        let n = reader.read(&mut tmp).unwrap();
        assert_eq!(n, 128);
        out.extend_from_slice(&tmp[..n]);
        assert_eq!(out.len(), 128);

        // Next read should surface the error.
        let err = reader.read(&mut tmp).expect_err("second read should error");
        assert!(matches!(err.kind(), io::ErrorKind::Other | io::ErrorKind::UnexpectedEof));
    }

    #[test]
    fn interrupted_errors_are_retried_transparently() {
        struct FlakyThenEof {
            call: usize,
        }
        impl Read for FlakyThenEof {
            fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
                self.call += 1;
                if self.call <= 3 {
                    return Err(io::Error::new(io::ErrorKind::Interrupted, "try again"));
                }
                if self.call == 4 {
                    let n = buf.len().min(10);
                    for (i, b) in buf.iter_mut().take(n).enumerate() {
                        *b = u8::try_from(i + 1).expect("n <= 10 so i+1 fits in u8");
                    }
                    return Ok(n);
                }
                Ok(0)
            }
        }
        let mut reader = PrefetchReader::with_config(FlakyThenEof { call: 0 }, 64, 2);
        let mut out = Vec::new();
        reader.read_to_end(&mut out).unwrap();
        assert_eq!(out, (1..=10).collect::<Vec<u8>>());
    }

    #[test]
    fn drop_joins_producer_thread() {
        /// A reader that sets a flag inside its `Drop` impl so we can
        /// observe that the producer thread has actually been torn down.
        struct Tracked {
            flag: Arc<AtomicBool>,
            data: Vec<u8>,
            pos: usize,
        }
        impl Read for Tracked {
            fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
                if self.pos >= self.data.len() {
                    return Ok(0);
                }
                let n = buf.len().min(self.data.len() - self.pos);
                buf[..n].copy_from_slice(&self.data[self.pos..self.pos + n]);
                self.pos += n;
                Ok(n)
            }
        }
        impl Drop for Tracked {
            fn drop(&mut self) {
                self.flag.store(true, Ordering::SeqCst);
            }
        }

        let flag = Arc::new(AtomicBool::new(false));
        let inner = Tracked { flag: Arc::clone(&flag), data: sample_bytes(1024), pos: 0 };
        {
            let mut reader = PrefetchReader::with_config(inner, 64, 2);
            let mut out = Vec::new();
            reader.read_to_end(&mut out).unwrap();
            assert_eq!(out.len(), 1024);
            // reader is dropped here; `Tracked::drop` must fire on the
            // producer thread before `Drop::drop` returns.
        }
        assert!(
            flag.load(Ordering::SeqCst),
            "producer thread should have dropped the inner reader"
        );
    }

    proptest! {
        /// Property test: for any input bytes, any chunk size, any prefetch
        /// depth, and any consumer read size, `PrefetchReader` yields exactly
        /// the same byte sequence as a plain in-memory cursor.
        #[test]
        fn prop_byte_identical_to_cursor(
            data in prop::collection::vec(any::<u8>(), 0..8_192),
            chunk_size in 1usize..2_048,
            depth in 1usize..6,
            read_size in 1usize..256,
        ) {
            let expected = data.clone();
            let mut reader = PrefetchReader::with_config(
                Cursor::new(data),
                chunk_size,
                depth,
            );
            let mut out = Vec::with_capacity(expected.len());
            let mut tmp = vec![0u8; read_size];
            loop {
                let n = reader.read(&mut tmp).expect("read should not fail on Cursor");
                if n == 0 {
                    break;
                }
                out.extend_from_slice(&tmp[..n]);
            }
            prop_assert_eq!(out, expected);
        }
    }
}

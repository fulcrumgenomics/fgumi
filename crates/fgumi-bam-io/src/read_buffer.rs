#![deny(unsafe_code)]
//! A byte-bounded read-ahead buffer for streamed input.
//!
//! [`ReadBuffer`] is a [`Read`] adapter whose dedicated thread keeps reading
//! the inner reader into a buffer of up to a fixed number of bytes, whether or
//! not the consumer is reading. It exists for inputs that arrive through a
//! pipe: when the consumer stops reading for a few seconds (a sort sealing a
//! run: sorting it, compressing it and writing it to disk), a pipe's own
//! buffer fills in well under a second and blocks the upstream process. With a
//! [`ReadBuffer`] the upstream process keeps writing until this buffer, not
//! the pipe's, is full.
//!
//! It differs from [`PrefetchReader`](crate::prefetch_reader::PrefetchReader),
//! which hides disk latency with a few chunks of read-ahead, in three ways that
//! matter for a buffer of gigabytes fed by a pipe:
//!
//! - **Bounded by bytes.** A pipe returns at most its own capacity (64 KiB by
//!   default on Linux) per `read`. `PrefetchReader` sends each read as one
//!   chunk, so a deep chunk queue would hold mostly empty chunk capacity. Here
//!   the producer fills each chunk before sending it, except that it sends a
//!   partly filled chunk when the queue is empty right after a read returns, so
//!   a starved consumer gets those bytes at once. If the queue was not empty
//!   then and the consumer drains it while the producer is blocked in its next
//!   `read`, the bytes already read wait until that read returns (the upstream
//!   writes again), the input ends, or it fails. Bytes are never lost or
//!   reordered, and the end of input or an error delivers everything read
//!   before it.
//! - **Chunks are recycled.** Drained chunks go back to the producer instead of
//!   being freed, so a buffer that fills and drains at every stall does not
//!   page-fault a fresh allocation in each time.
//! - **It reports how full it got** (logged when it is dropped), which says
//!   whether the buffer absorbed the consumer's stalls or ran out.
//!
//! ## Memory
//!
//! At most `budget` bytes of chunk storage are allocated (for any budget of at
//! least 8 bytes): chunks are an eighth of the budget, up to
//! [`MAX_CHUNK_SIZE`], and the queue holds up to `budget / chunk - 2` of them,
//! plus one being filled and one being served. Chunks are allocated only as
//! the queue first grows, so a consumer that keeps up holds a few chunks, not
//! the budget. Once the consumer has read to the end of the input, the buffer
//! holds no chunks: the drained ones are freed rather than recycled.
//!
//! ## Lifecycle
//!
//! The thread (`fgumi-read-buffer`) owns the inner reader and exits at EOF,
//! after forwarding a read error (the bytes read before the error are
//! delivered first), or when the [`ReadBuffer`] is dropped. `Drop` joins it if
//! it has exited. Otherwise (blocked sending into a full buffer, or parked in a
//! `read` on a quiet upstream when the buffer is dropped early on an error
//! path) `Drop` detaches it rather than wait: the dropped receiver fails its
//! next send, so it exits once that send is attempted, or it ends with the
//! process. A read error is sticky: every `read` after it fails the same way.

use std::io::{self, Read};
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};
use std::thread::{self, JoinHandle};
use std::time::Instant;

use crossbeam_channel::{Receiver, Sender, TrySendError, bounded};

/// Largest chunk the producer fills before handing it to the consumer: the
/// most a consumer that has read to the end of its input can still hold.
pub const MAX_CHUNK_SIZE: usize = 4 * 1024 * 1024;

/// A filled chunk: its storage (always `chunk_size` long, so it can be reused
/// without zero-filling) and how many leading bytes are data.
struct Chunk {
    buf: Vec<u8>,
    len: usize,
}

/// What the producer hands the consumer: a chunk, or the error that ended the
/// stream.
type Item = io::Result<Chunk>;

/// Counters shared by the producer and the consumer. Read them through
/// [`ReadBuffer::stats`].
#[derive(Debug, Default)]
struct SharedStats {
    /// Bytes currently queued for the consumer.
    queued: AtomicU64,
    /// The most bytes ever queued at once.
    peak_queued: AtomicU64,
    /// Sends that found the queue full and waited for the consumer.
    full_waits: AtomicU64,
    /// Total time spent in those waits.
    full_wait_ns: AtomicU64,
    /// Bytes read from the inner reader.
    bytes_read: AtomicU64,
    /// Chunks allocated because none was waiting to be recycled.
    allocations: AtomicU64,
}

/// A snapshot of a [`ReadBuffer`]'s counters.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct ReadBufferStats {
    /// The byte budget the buffer was built with.
    pub(crate) budget: u64,
    /// The most bytes that were queued for the consumer at once.
    pub(crate) peak_queued: u64,
    /// How many times the buffer was full, so the producer stopped reading the
    /// inner reader until the consumer drained a chunk.
    pub(crate) full_waits: u64,
    /// Total time the producer spent waiting on a full buffer, in nanoseconds.
    /// While it waits, an upstream process writing into a pipe can block.
    pub(crate) full_wait_ns: u64,
    /// Bytes read from the inner reader so far.
    pub(crate) bytes_read: u64,
}

/// A [`Read`] adapter backed by a byte-bounded buffer that a dedicated thread
/// keeps filling. See the [module docs](self).
pub struct ReadBuffer {
    /// The chunk being served to `read` callers and the offset of its next
    /// unread byte.
    current: Option<(Chunk, usize)>,
    /// Filled chunks from the producer. Dropping it (with `self`) fails the
    /// producer's next send, which ends it.
    rx: Receiver<Item>,
    /// Drained chunks returned to the producer for reuse; `None` once the
    /// stream has ended, so drained chunks are freed instead (and the ones
    /// still waiting in the return queue are freed with it).
    recycle_tx: Option<Sender<Vec<u8>>>,
    /// The producer thread; `None` after it is joined.
    handle: Option<JoinHandle<()>>,
    /// Counters shared with the producer.
    stats: Arc<SharedStats>,
    /// The byte budget, for [`Self::stats`].
    budget: u64,
    /// The error that ended the stream, returned again by every later `read`
    /// so a caller that reads past it cannot mistake it for a clean end.
    failed: Option<(io::ErrorKind, String)>,
}

impl std::fmt::Debug for ReadBuffer {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ReadBuffer").field("stats", &self.stats()).finish_non_exhaustive()
    }
}

/// The chunk size and queue depth for a byte budget: chunks of an eighth of
/// the budget, at most [`MAX_CHUNK_SIZE`], and as many queued chunks as the
/// budget holds after the two in flight (at least one).
fn layout(budget: usize) -> (usize, usize) {
    let chunk_size = (budget / 8).clamp(1, MAX_CHUNK_SIZE);
    let depth = (budget / chunk_size).saturating_sub(2).max(1);
    (chunk_size, depth)
}

impl ReadBuffer {
    /// Start reading `inner` on a new thread into a buffer of at most about
    /// `budget` bytes (see the [module docs](self#memory) for the exact
    /// bound).
    ///
    /// # Errors
    ///
    /// Returns an `InvalidInput` error if `budget` is zero, and the spawn error
    /// if the thread cannot be started.
    pub fn try_new<R: Read + Send + 'static>(inner: R, budget: usize) -> io::Result<Self> {
        if budget == 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "a read buffer needs a budget of at least one byte",
            ));
        }
        let (chunk_size, depth) = layout(budget);
        Self::with_layout(inner, budget, chunk_size, depth)
    }

    /// [`Self::try_new`] with an explicit chunk size and queue depth, for
    /// tests.
    fn with_layout<R: Read + Send + 'static>(
        inner: R,
        budget: usize,
        chunk_size: usize,
        depth: usize,
    ) -> io::Result<Self> {
        let (tx, rx) = bounded::<Item>(depth);
        // Every chunk is in the queue, being filled, being served, or here, so
        // `depth + 2` slots never refuse a returned chunk.
        let (recycle_tx, recycle_rx) = bounded::<Vec<u8>>(depth + 2);
        let stats = Arc::new(SharedStats::default());
        let producer = Producer { chunk_size, tx, recycle_rx, stats: Arc::clone(&stats) };
        let handle = thread::Builder::new()
            .name("fgumi-read-buffer".to_string())
            .spawn(move || producer.run(inner))?;
        Ok(Self {
            current: None,
            rx,
            recycle_tx: Some(recycle_tx),
            handle: Some(handle),
            stats,
            budget: budget as u64,
            failed: None,
        })
    }

    /// A snapshot of the buffer's counters.
    #[must_use]
    pub(crate) fn stats(&self) -> ReadBufferStats {
        ReadBufferStats {
            budget: self.budget,
            peak_queued: self.stats.peak_queued.load(Ordering::Relaxed),
            full_waits: self.stats.full_waits.load(Ordering::Relaxed),
            full_wait_ns: self.stats.full_wait_ns.load(Ordering::Relaxed),
            bytes_read: self.stats.bytes_read.load(Ordering::Relaxed),
        }
    }

    /// Hand a drained chunk back to the producer. It is freed instead once the
    /// stream has ended, if the producer has exited, or (impossibly) if the
    /// return queue is full.
    fn recycle(&self, chunk: Chunk) {
        if let Some(recycle_tx) = &self.recycle_tx {
            let _ = recycle_tx.try_send(chunk.buf);
        }
    }

    /// The stream has ended (EOF or an error): stop recycling. Dropping the
    /// return queue's sender frees the chunks still in it once the producer,
    /// which no longer needs them, has dropped its receiver.
    fn release_recycled(&mut self) {
        self.recycle_tx = None;
    }
}

/// The producer thread's state.
struct Producer {
    chunk_size: usize,
    tx: Sender<Item>,
    recycle_rx: Receiver<Vec<u8>>,
    stats: Arc<SharedStats>,
}

impl Producer {
    /// Run the read loop; a panic in the inner reader reaches the consumer as
    /// an error rather than as a silently truncated stream.
    fn run<R: Read>(self, inner: R) {
        let result =
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| self.read_loop(inner)));
        if let Err(payload) = result {
            let msg = crate::panic_message(payload.as_ref(), "fgumi-read-buffer thread panicked");
            let _ = self.tx.send(Err(io::Error::other(msg)));
        }
    }

    /// Fill chunks from `inner` and send them until EOF, an error, or the
    /// consumer is gone.
    ///
    /// A chunk is sent when it is full, at EOF, before an error, or when the
    /// queue is empty right after a read returns — the consumer may be waiting,
    /// and holding those bytes back from it would only add latency. So at most
    /// one partly filled chunk is ever queued: one is sent only into an empty
    /// queue. The emptiness is checked only when a read returns, so if the
    /// queue was not empty then and the consumer drains it while this thread is
    /// blocked in the next `read`, the bytes already in the chunk wait until
    /// that read returns (the upstream writes again), EOF or an error. Nothing
    /// is lost or reordered, and EOF and errors send everything read before
    /// them.
    fn read_loop<R: Read>(&self, mut inner: R) {
        loop {
            let mut buf = self.recycle_rx.try_recv().unwrap_or_else(|_| {
                self.stats.allocations.fetch_add(1, Ordering::Relaxed);
                vec![0u8; self.chunk_size]
            });
            let mut len = 0;
            let mut end: Option<io::Error> = None;
            let mut eof = false;
            while len < buf.len() {
                match inner.read(&mut buf[len..]) {
                    Ok(0) => {
                        eof = true;
                        break;
                    }
                    Ok(n) => {
                        len += n;
                        self.stats.bytes_read.fetch_add(n as u64, Ordering::Relaxed);
                        if self.tx.is_empty() {
                            break;
                        }
                    }
                    Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
                    Err(e) => {
                        end = Some(e);
                        break;
                    }
                }
            }
            if len > 0 && !self.send(Chunk { buf, len }) {
                return;
            }
            if let Some(e) = end {
                let _ = self.tx.send(Err(e));
                return;
            }
            if eof {
                return;
            }
        }
    }

    /// Queue `chunk` for the consumer, waiting while the buffer is full.
    /// Returns `false` if the consumer is gone.
    fn send(&self, chunk: Chunk) -> bool {
        let len = chunk.len as u64;
        let queued = self.stats.queued.fetch_add(len, Ordering::Relaxed) + len;
        self.stats.peak_queued.fetch_max(queued, Ordering::Relaxed);
        let item = match self.tx.try_send(Ok(chunk)) {
            Ok(()) => return true,
            Err(TrySendError::Disconnected(_)) => return false,
            Err(TrySendError::Full(item)) => item,
        };
        // Counted when the wait starts, so a wait still in progress shows.
        self.stats.full_waits.fetch_add(1, Ordering::Relaxed);
        let started = Instant::now();
        let sent = self.tx.send(item).is_ok();
        let waited = u64::try_from(started.elapsed().as_nanos()).unwrap_or(u64::MAX);
        self.stats.full_wait_ns.fetch_add(waited, Ordering::Relaxed);
        sent
    }
}

impl Read for ReadBuffer {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        if out.is_empty() {
            return Ok(0);
        }
        if let Some((kind, msg)) = &self.failed {
            return Err(io::Error::new(*kind, msg.clone()));
        }
        loop {
            if let Some((chunk, pos)) = self.current.as_mut() {
                if *pos < chunk.len {
                    let n = out.len().min(chunk.len - *pos);
                    out[..n].copy_from_slice(&chunk.buf[*pos..*pos + n]);
                    *pos += n;
                    return Ok(n);
                }
                if let Some((chunk, _)) = self.current.take() {
                    self.recycle(chunk);
                }
            }
            // Disconnected: the producer exited after sending everything.
            let Ok(item) = self.rx.recv() else {
                self.release_recycled();
                return Ok(0);
            };
            match item {
                Ok(chunk) => {
                    self.stats.queued.fetch_sub(chunk.len as u64, Ordering::Relaxed);
                    self.current = Some((chunk, 0));
                }
                Err(e) => {
                    self.release_recycled();
                    self.failed = Some((e.kind(), e.to_string()));
                    return Err(e);
                }
            }
        }
    }
}

impl Drop for ReadBuffer {
    fn drop(&mut self) {
        // Join the producer only if it has already exited. One still running is
        // detached: one parked in a `read` on a quiet upstream (dropped early,
        // on an error path) would hold up the exit for as long as the upstream
        // stays quiet. `rx` is dropped with `self` right after this, which
        // fails the producer's next send, so it exits once it tries one (at
        // once if it is blocked on a full buffer), or with the process.
        self.current = None;
        if let Some(handle) = self.handle.take() {
            if !handle.is_finished() {
                log::debug!("fgumi-read-buffer thread still reading at drop; detached");
            } else if handle.join().is_err() {
                log::debug!("fgumi-read-buffer thread panicked during shutdown");
            }
        }
        log::info!("{}", self.stats().summary());
    }
}

impl ReadBufferStats {
    /// One line for the log: how much was read, how full the buffer got, and
    /// how long it was full (time during which a writer upstream could block).
    #[must_use]
    pub(crate) fn summary(&self) -> String {
        #[allow(clippy::cast_precision_loss)] // display only
        let mib = |bytes: u64| bytes as f64 / f64::from(1u32 << 20);
        #[allow(clippy::cast_precision_loss)] // display only
        let full_s = self.full_wait_ns as f64 / 1e9;
        format!(
            "Read buffer: read {:.1} MiB; held at most {:.1} of {:.1} MiB; full {} times for \
             {full_s:.1} s",
            mib(self.bytes_read),
            mib(self.peak_queued),
            mib(self.budget),
            self.full_waits,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;
    use rstest::rstest;
    use std::io::Cursor;
    use std::sync::atomic::AtomicUsize;
    use std::time::Duration;

    /// A convenient unit for test budgets.
    const KIB64: usize = 64 * 1024;

    /// Deterministic sample bytes.
    fn sample_bytes(len: usize) -> Vec<u8> {
        (0..len).map(|i| u8::try_from(i % 251).expect("i % 251 fits in u8")).collect()
    }

    /// An endless reader that hands out `step` bytes per call and counts them,
    /// so a test can see how far the producer got without the consumer.
    struct Counting {
        step: usize,
        read: Arc<AtomicUsize>,
    }

    impl Read for Counting {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            let n = buf.len().min(self.step);
            buf[..n].fill(7);
            self.read.fetch_add(n, Ordering::SeqCst);
            Ok(n)
        }
    }

    /// A [`ReadBuffer`] of `budget` bytes over `inner`.
    fn buffer<R: Read + Send + 'static>(inner: R, budget: usize) -> ReadBuffer {
        ReadBuffer::try_new(inner, budget).expect("start the read buffer")
    }

    /// Wait (up to 10 s) until the producer has found the buffer full more
    /// than `seen` times; returns the new count. `full_waits` is bumped only
    /// after a chunk is filled and the queue refused it, so the producer is
    /// then blocked with the buffer full and reads nothing more until the
    /// consumer drains a chunk.
    fn wait_until_full(reader: &ReadBuffer, seen: u64) -> u64 {
        let deadline = Instant::now() + Duration::from_secs(10);
        loop {
            let waits = reader.stats().full_waits;
            if waits > seen {
                return waits;
            }
            assert!(Instant::now() < deadline, "the buffer never filled: {:?}", reader.stats());
            thread::sleep(Duration::from_millis(1));
        }
    }

    /// Wraps a reader and signals `dropped` when it is dropped, which happens
    /// when the producer's read loop returns.
    struct DropSignal<R> {
        inner: R,
        dropped: std::sync::mpsc::Sender<()>,
    }

    impl<R: Read> Read for DropSignal<R> {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            self.inner.read(buf)
        }
    }

    impl<R> Drop for DropSignal<R> {
        fn drop(&mut self) {
            let _ = self.dropped.send(());
        }
    }

    #[rstest]
    #[case::tiny(64, 8, 6)]
    #[case::below_eight_bytes(4, 1, 2)]
    #[case::kib64(KIB64, 8 * 1024, 6)]
    #[case::small_budget(1024 * 1024, 128 * 1024, 6)]
    #[case::at_the_largest_chunk(32 << 20, MAX_CHUNK_SIZE, 6)]
    #[case::two_gib(2 << 30, MAX_CHUNK_SIZE, 510)]
    fn layout_follows_the_budget(
        #[case] budget: usize,
        #[case] chunk_size: usize,
        #[case] depth: usize,
    ) {
        assert_eq!(layout(budget), (chunk_size, depth));
    }

    #[rstest]
    #[case::empty(0)]
    #[case::one_byte(1)]
    #[case::many_chunks(3 * KIB64 + 17)]
    fn delivers_every_byte_in_order(#[case] len: usize) {
        let data = sample_bytes(len);
        let mut reader = buffer(Cursor::new(data.clone()), 4 * KIB64);
        let mut out = Vec::new();
        reader.read_to_end(&mut out).unwrap();
        assert_eq!(out, data);
        assert_eq!(reader.stats().bytes_read, len as u64);
        // EOF is sticky.
        assert_eq!(reader.read(&mut [0u8; 8]).unwrap(), 0);
        // At EOF the buffer stops recycling, so it holds no chunks.
        assert!(reader.recycle_tx.is_none() && reader.current.is_none());
    }

    #[test]
    fn a_zero_budget_is_an_error() {
        let err = ReadBuffer::try_new(Cursor::new(Vec::new()), 0).expect_err("zero budget");
        assert_eq!(err.kind(), io::ErrorKind::InvalidInput);
    }

    /// The point of the buffer: while the consumer reads nothing, the producer
    /// keeps reading the inner reader until the budget is queued, and then
    /// stops — it neither blocks early nor grows past the budget.
    #[test]
    fn keeps_reading_while_the_consumer_stalls_up_to_the_budget() {
        let budget = 32 * KIB64;
        let (chunk_size, depth) = layout(budget);
        let read = Arc::new(AtomicUsize::new(0));
        let reader = buffer(Counting { step: 4096, read: Arc::clone(&read) }, budget);
        wait_until_full(&reader, 0);
        let held = read.load(Ordering::SeqCst);
        // The queue is full (`depth` chunks, the first only partly filled: it
        // was sent into an empty queue) and one more chunk is filled and
        // waiting to be sent; the consumer holds none.
        assert!(held >= depth * chunk_size, "read only {held} bytes of a {budget}-byte buffer");
        assert!(held <= (depth + 1) * chunk_size, "read {held} bytes, past the {budget} budget");
        let stats = reader.stats();
        assert!(stats.peak_queued >= ((depth - 1) * chunk_size) as u64, "{stats:?}");
        assert!(stats.peak_queued <= budget as u64, "{stats:?}");
        assert!(stats.full_waits >= 1, "the producer never found the buffer full");
    }

    /// Draining the buffer lets the producer read further, and drained chunks
    /// are reused: across several fill/drain cycles the producer allocates no
    /// more chunks than the queue holds plus the two in flight.
    #[test]
    fn resumes_after_the_consumer_drains() {
        let budget = 16 * KIB64;
        let (_, depth) = layout(budget);
        let read = Arc::new(AtomicUsize::new(0));
        let mut reader = buffer(Counting { step: 65536, read: Arc::clone(&read) }, budget);
        let mut waits = wait_until_full(&reader, 0);
        let mut sink = vec![0u8; budget];
        for _ in 0..4 {
            let before = read.load(Ordering::SeqCst);
            reader.read_exact(&mut sink).unwrap();
            waits = wait_until_full(&reader, waits);
            let after = read.load(Ordering::SeqCst);
            assert!(after >= before + budget / 2, "producer did not resume: {before} -> {after}");
        }
        assert!(reader.stats().peak_queued <= budget as u64);
        let allocations = reader.stats.allocations.load(Ordering::Relaxed);
        assert!(allocations <= (depth + 2) as u64, "{allocations} chunks allocated, depth {depth}");
    }

    /// A pipe hands out small reads; the producer must still fill whole chunks
    /// while the consumer is away, or the buffer would hold mostly empty
    /// capacity.
    #[test]
    fn small_inner_reads_still_fill_the_buffer() {
        let budget = 32 * KIB64;
        let read = Arc::new(AtomicUsize::new(0));
        let reader = buffer(Counting { step: 512, read: Arc::clone(&read) }, budget);
        wait_until_full(&reader, 0);
        let (chunk_size, depth) = layout(budget);
        // At most one partly filled chunk is ever queued.
        assert!(reader.stats().peak_queued >= ((depth - 1) * chunk_size) as u64);
    }

    #[test]
    fn an_error_is_delivered_after_the_bytes_read_before_it() {
        struct DataThenErr {
            sent: bool,
        }
        impl Read for DataThenErr {
            fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
                if self.sent {
                    return Err(io::Error::new(io::ErrorKind::BrokenPipe, "upstream died"));
                }
                self.sent = true;
                let n = buf.len().min(1000);
                buf[..n].fill(9);
                Ok(n)
            }
        }
        let mut reader = buffer(DataThenErr { sent: false }, 4 * KIB64);
        let mut out = Vec::new();
        let err = reader.read_to_end(&mut out).expect_err("the error must reach the consumer");
        assert_eq!(err.kind(), io::ErrorKind::BrokenPipe);
        assert_eq!(out, vec![9u8; 1000]);
        // Sticky: reading past the error fails again, never a clean end.
        for _ in 0..2 {
            let again = reader.read(&mut [0u8; 8]).expect_err("the error must stay");
            assert_eq!(again.kind(), io::ErrorKind::BrokenPipe);
        }
    }

    #[test]
    fn a_panicking_reader_surfaces_as_an_error() {
        struct Panics;
        impl Read for Panics {
            fn read(&mut self, _: &mut [u8]) -> io::Result<usize> {
                panic!("inner reader blew up");
            }
        }
        let mut reader = buffer(Panics, KIB64);
        let err = reader.read(&mut [0u8; 8]).expect_err("panic must surface as an error");
        assert!(err.to_string().contains("inner reader blew up"), "{err}");
    }

    #[test]
    fn interrupted_reads_are_retried() {
        struct Flaky {
            calls: usize,
        }
        impl Read for Flaky {
            fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
                self.calls += 1;
                match self.calls {
                    1..=3 => Err(io::Error::new(io::ErrorKind::Interrupted, "again")),
                    4 => {
                        buf[..3].copy_from_slice(b"abc");
                        Ok(3)
                    }
                    _ => Ok(0),
                }
            }
        }
        let mut reader = buffer(Flaky { calls: 0 }, KIB64);
        let mut out = Vec::new();
        reader.read_to_end(&mut out).unwrap();
        assert_eq!(out, b"abc");
    }

    /// Dropping the reader with a full buffer and an endless inner reader must
    /// not hang, and the producer, blocked sending into the full buffer, must
    /// then exit: its send fails once the receiver is gone.
    #[test]
    fn drop_with_a_full_buffer_does_not_hang() {
        let read = Arc::new(AtomicUsize::new(0));
        let (dropped_tx, dropped) = std::sync::mpsc::channel();
        let counting = Counting { step: 4096, read: Arc::clone(&read) };
        let reader = buffer(DropSignal { inner: counting, dropped: dropped_tx }, 4 * KIB64);
        wait_until_full(&reader, 0);
        drop(reader);
        assert!(
            dropped.recv_timeout(Duration::from_secs(5)).is_ok(),
            "the producer did not exit after the reader was dropped"
        );
    }

    /// Dropping the reader while the producer is blocked in a `read` on a quiet
    /// upstream (as on an error path before the input ends) returns promptly,
    /// and the detached producer exits once that read returns.
    #[test]
    fn drop_while_the_producer_is_blocked_in_read_returns_promptly() {
        /// Hands out one byte, then blocks until the test ends.
        struct Quiet {
            gate: std::sync::mpsc::Receiver<()>,
            sent: bool,
        }
        impl Read for Quiet {
            fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
                if self.sent {
                    let _ = self.gate.recv();
                    return Ok(0);
                }
                self.sent = true;
                buf[0] = 1;
                Ok(1)
            }
        }
        let (gate_tx, gate) = std::sync::mpsc::channel::<()>();
        let (dropped_tx, dropped) = std::sync::mpsc::channel();
        let quiet = Quiet { gate, sent: false };
        let mut reader = buffer(DropSignal { inner: quiet, dropped: dropped_tx }, KIB64);
        assert_eq!(reader.read(&mut [0u8; 8]).unwrap(), 1);
        let (done_tx, done) = std::sync::mpsc::channel();
        thread::spawn(move || {
            drop(reader);
            let _ = done_tx.send(());
        });
        assert!(
            done.recv_timeout(Duration::from_secs(2)).is_ok(),
            "drop waited on a producer blocked in read"
        );
        drop(gate_tx);
        assert!(
            dropped.recv_timeout(Duration::from_secs(5)).is_ok(),
            "the detached producer did not exit after its read returned"
        );
    }

    proptest! {
        /// Any input, chunk size, depth and consumer read size yields exactly
        /// the input bytes.
        #[test]
        fn prop_byte_identical_to_the_input(
            data in prop::collection::vec(any::<u8>(), 0..8_192),
            chunk_size in 1usize..2_048,
            depth in 1usize..6,
            read_size in 1usize..256,
        ) {
            let mut reader = ReadBuffer::with_layout(Cursor::new(data.clone()), 1, chunk_size, depth)
                .expect("start the read buffer");
            let mut out = Vec::with_capacity(data.len());
            let mut tmp = vec![0u8; read_size];
            loop {
                let n = reader.read(&mut tmp).expect("a Cursor does not fail");
                if n == 0 {
                    break;
                }
                out.extend_from_slice(&tmp[..n]);
            }
            prop_assert_eq!(out, data);
        }
    }
}

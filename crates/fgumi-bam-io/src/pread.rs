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
}

#[cfg(unix)]
impl PositionalSource for std::fs::File {
    fn read_at(&self, buf: &mut [u8], offset: u64) -> io::Result<usize> {
        std::os::unix::fs::FileExt::read_at(self, buf, offset)
    }

    fn byte_len(&self) -> io::Result<u64> {
        Ok(self.metadata()?.len())
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
/// shrinking slices, or a file's short tail read) releases one mismatched
/// buffer per miss and converges to the sizes it is asked for. The gauges
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
    /// rather than a memset). On a miss one pooled buffer is released, so a
    /// free list of the wrong sizes drains instead of pinning memory. Its
    /// contents are unspecified; a caller overwrites every byte before leasing
    /// it.
    #[must_use]
    pub fn take(&self, len: usize) -> Vec<u8> {
        let (hit, evicted) = {
            let mut free = self.free.lock();
            match free.iter().position(|b| fits(b.capacity(), len)) {
                Some(i) => (Some(free.swap_remove(i)), None),
                None => (None, free.pop()),
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
}

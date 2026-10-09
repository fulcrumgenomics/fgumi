//! Concurrent positional-read source that raises a device's read queue depth.
//!
//! On a slow, deep-queue device (EBS gp3 is the motivating case) a single
//! outstanding `read()` leaves most of the device's bandwidth idle: gp3
//! sustains ~358 MB/s at queue-depth-1 but ~1177 MB/s with four reads in
//! flight. [`ScatterReader`] cuts each fill window into N byte-range slices,
//! fetches them with independent concurrent positional reads
//! (`pread`/[`std::os::unix::fs::FileExt::read_at`]), reassembles them **in file
//! order**, and presents the result as an ordinary forward [`Read`] stream —
//! so the downstream BGZF framer is untouched. This ports the mechanism of
//! fgumi v0.7.0's `fgumi sort --read-streams` flag (PRs #838/#846) onto the
//! current chain source.
//!
//! No `unsafe`: each slice is read into a `Vec<u8>` owned solely by the scoped
//! thread that fills it (via [`std::thread::scope`]), and positional reads take
//! `&self` (POSIX `pread(2)` touches no shared file offset), so concurrent
//! reads on a shared `Arc<File>` need no synchronization.

use std::io::{self, Read};
use std::path::Path;
use std::sync::Arc;
use std::time::{Duration, Instant};

/// Target aggregate throughput the [`ReadStreams::Auto`] stream count is chosen
/// to clear (`ceil(target / measured)`), from PR #846 (b358928f). Set from what
/// the consumer (decode + sort) can absorb, not the device — there is no point
/// reading faster than the pipeline can process. Validated on the arena chain:
/// a 38 GB coordinate sort on EBS gp3 (500 MB/s) ran ~14% faster at four streams
/// than at one, and `Auto` resolved to four there, as designed.
///
/// [`ReadStreams::Auto`]: crate::ReadStreams::Auto
pub(crate) const AUTO_TARGET_BYTES_PER_SEC: f64 = 1.2e9;

/// Upper bound on the chosen stream count (b358928f / 9a80d7fe). Applies to both
/// the `Auto` probe outcome and an explicit `Fixed(n)` (clamped in
/// [`ScatterReader::from_parts`]).
pub(crate) const MAX_STREAMS: usize = 8;

/// Number of fills read at a single stream before `Auto` commits to a count.
/// While probing, `active_streams` is 1, so each fill is a single positional
/// read — the same queue depth `Fixed(1)` would use, and the throughput the
/// probe is trying to measure. (It still travels `ScatterReader`'s own per-fill
/// buffer path, not the bare-`File` path `Fixed(1)` takes, so it is cheap but
/// not literally free.)
pub(crate) const AUTO_PROBE_FILLS: usize = 8;

/// Total bytes per fill window, split across the active stream count. Matches
/// `PrefetchReader`'s chunk size.
pub(crate) const DEFAULT_FILL_SIZE: usize = 4 * 1024 * 1024;

/// Smallest slice a fill is split into. A window smaller than
/// `active_streams * MIN_SLICE_SIZE` uses fewer, larger slices rather than many
/// tiny reads whose per-read overhead would swamp the queue-depth win.
pub(crate) const MIN_SLICE_SIZE: usize = 64 * 1024;

use crate::ReadStreams;

/// A byte-addressable source supporting concurrent positional reads.
///
/// Exists as a seam so [`ScatterReader`]'s reassembly / probe / fallback logic
/// can be unit-tested against an in-memory or fault-injecting double without a
/// real file.
pub trait PositionalSource: Send + Sync {
    /// Read starting at `offset` into `buf`, returning the number of bytes read
    /// (may be short, as [`Read::read`]; `0` means EOF at `offset`).
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

/// Whether `file` supports the scatter reader's concurrent positional reads.
///
/// Regular files only: a FIFO or socket opened by path fails `pread(2)` with
/// `ESPIPE`, so callers fall back to the sequential/async reader rather than
/// fail the run over a performance optimization. A `metadata()` failure (e.g. a
/// race with deletion) is treated as not-capable.
#[cfg(unix)]
#[must_use]
pub fn is_scatter_capable(file: &std::fs::File) -> bool {
    file.metadata().map(|m| m.is_file()).unwrap_or(false)
}

/// Non-Unix targets have no `read_at`; the scatter path is never selected.
#[cfg(not(unix))]
#[must_use]
pub fn is_scatter_capable(_file: &std::fs::File) -> bool {
    false
}

/// The reader [`decide_reader`] selected for an open regular `file`.
///
/// The scatter case is boxed as `dyn Read` rather than a concrete
/// `ScatterReader<File>` so the type is well-formed on non-Unix targets too
/// (where `File` is not a [`PositionalSource`] and the variant is never built).
pub enum ScatterDecision {
    /// The input is seekable and `--read-streams` asked for concurrency: read it
    /// through this scatter reader (a boxed [`ScatterReader`]).
    Scatter(Box<dyn Read + Send>),
    /// Scatter did not apply (`Fixed(1)`, a non-seekable input, or a non-Unix
    /// platform): the caller reads `file` with its own sequential/async reader.
    Fallback(std::fs::File),
}

/// Decide how to read an open regular `file` under the `--read-streams` policy.
///
/// Centralizes the scatter-vs-sequential decision (and its warning) so the two
/// production entry points — `open_normalized_with_opts` (the sort chain) and
/// [`read_bam`](crate) — cannot drift. Returns [`ScatterDecision::Scatter`] when
/// `read_streams` requests concurrency and `file` is a seekable regular file;
/// otherwise [`ScatterDecision::Fallback`], handing `file` back for the caller's
/// own sequential/async reader. Warns only when the user *explicitly* asked for
/// `Fixed(n>1)` that cannot be honored — `Auto` is a best-effort default and
/// falls back silently. `label` names the source in the scatter info log.
///
/// # Errors
/// Returns an I/O error only if the seekable file's length cannot be read.
pub fn decide_reader(
    file: std::fs::File,
    read_streams: ReadStreams,
    path: &Path,
    label: &str,
) -> io::Result<ScatterDecision> {
    // `Fixed(1)` is the plain reader: never scatter, never warn.
    if read_streams == ReadStreams::Fixed(1) {
        return Ok(ScatterDecision::Fallback(file));
    }
    #[cfg(unix)]
    {
        if is_scatter_capable(&file) {
            log::info!("read-streams={read_streams}: scatter {label} on {}", path.display());
            return Ok(ScatterDecision::Scatter(Box::new(ScatterReader::from_file(
                file,
                read_streams,
            )?)));
        }
        warn_read_streams_unavailable(
            read_streams,
            &path.display().to_string(),
            "is not a seekable regular file",
        );
    }
    #[cfg(not(unix))]
    {
        let _ = label;
        warn_read_streams_unavailable(
            read_streams,
            &path.display().to_string(),
            "is not supported on this platform",
        );
    }
    Ok(ScatterDecision::Fallback(file))
}

/// Warn that `--read-streams` cannot be honored for `subject` (`reason` says
/// why), but only when the user *explicitly* requested concurrency (`Fixed(n>1)`).
/// `Auto` is the default and falls back silently; `Fixed(1)` never reaches here.
/// Shared by [`decide_reader`] and the stdin paths so the wording is identical
/// everywhere a non-seekable input drops the knob.
pub fn warn_read_streams_unavailable(read_streams: ReadStreams, subject: &str, reason: &str) {
    if let ReadStreams::Fixed(n) = read_streams
        && n > 1
    {
        log::warn!(
            "--read-streams={n} requested but {subject} {reason}; \
             falling back to the sequential/async reader"
        );
    }
}

/// Pick a stream count from a measured single-stream throughput.
///
/// Pure and total: a non-positive measurement (which should not happen, but a
/// probe over an empty file could report it) degrades to `1` rather than
/// dividing by zero. The result is clamped to `[1, MAX_STREAMS]`.
#[must_use]
#[allow(
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "n is clamped into (1.0, MAX_STREAMS as f64) before the cast, so \
              the usize conversion is exact and non-negative"
)]
pub(crate) fn choose_stream_count(measured_bytes_per_sec: f64) -> usize {
    if measured_bytes_per_sec <= 0.0 {
        return 1;
    }
    let n = (AUTO_TARGET_BYTES_PER_SEC / measured_bytes_per_sec).ceil();
    if n <= 1.0 {
        1
    } else if n >= MAX_STREAMS as f64 {
        MAX_STREAMS
    } else {
        n as usize
    }
}

/// Read exactly `buf.len()` bytes at `offset`, looping over short reads and
/// retrying `Interrupted`. A `0`-byte read before the buffer is full means the
/// source ended earlier than its reported length promised.
fn read_at_exact<S: PositionalSource + ?Sized>(
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
                    "scatter reader hit EOF before filling a slice",
                ));
            }
            Ok(n) => filled += n,
            Err(ref e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    Ok(())
}

/// Split `[start, start + window)` into up to `max_slices` disjoint contiguous
/// ranges, each at least `min_slice` bytes (except possibly the final one when
/// `window` is not a multiple). Returns `(offset, len)` pairs in file order.
fn slice_ranges(
    start: u64,
    window: usize,
    max_slices: usize,
    min_slice: usize,
) -> Vec<(u64, usize)> {
    debug_assert!(window > 0 && max_slices >= 1 && min_slice >= 1);
    let by_min = window.div_ceil(min_slice).max(1);
    let n = max_slices.min(by_min);
    let base = window / n;
    let rem = window % n;
    let mut ranges = Vec::with_capacity(n);
    let mut off = start;
    for i in 0..n {
        let len = base + usize::from(i < rem);
        if len == 0 {
            continue;
        }
        ranges.push((off, len));
        off += len as u64;
    }
    ranges
}

/// A forward [`Read`] stream backed by concurrent positional reads.
///
/// See the module docs for the mechanism. Construct with [`Self::from_file`]
/// (the production entry point); the crate-internal `new` / `with_params`
/// constructors build one over any [`PositionalSource`] for tests.
pub struct ScatterReader<S: PositionalSource + 'static> {
    source: Arc<S>,
    file_len: u64,
    next_offset: u64,
    /// Resolved stream count. Under `Auto` it starts at 1 and jumps once, after
    /// [`AUTO_PROBE_FILLS`] probe fills.
    active_streams: usize,
    probing: bool,
    probe_fills_done: usize,
    probe_bytes: u64,
    probe_elapsed: Duration,
    /// The assembled current fill and how far it has been consumed.
    current: Option<(Vec<u8>, usize)>,
    /// An error observed while assembling a fill, surfaced only after the bytes
    /// that were successfully read before it have been served.
    pending_error: Option<io::Error>,
    fill_size: usize,
    min_slice_size: usize,
    /// Test hook: when set, the `Auto` probe uses this throughput instead of the
    /// wall-clock measurement, so probe-outcome tests are deterministic rather
    /// than dependent on real timing under a loaded test runner. Always `None`
    /// in production.
    forced_throughput: Option<f64>,
}

impl<S: PositionalSource + 'static> ScatterReader<S> {
    fn from_parts(
        source: S,
        file_len: u64,
        read_streams: ReadStreams,
        fill_size: usize,
        min_slice_size: usize,
    ) -> Self {
        let (active_streams, probing) = match read_streams {
            ReadStreams::Auto => (1, true),
            // Clamp an explicit count to the same `[1, MAX_STREAMS]` ceiling the
            // `Auto` path uses (`choose_stream_count`), so `--read-streams 1000`
            // cannot spawn far more concurrent reads per fill than the heuristic
            // ever would.
            ReadStreams::Fixed(n) => (n.clamp(1, MAX_STREAMS), false),
        };
        Self {
            source: Arc::new(source),
            file_len,
            next_offset: 0,
            active_streams,
            probing,
            probe_fills_done: 0,
            probe_bytes: 0,
            probe_elapsed: Duration::ZERO,
            current: None,
            pending_error: None,
            fill_size,
            min_slice_size,
            forced_throughput: None,
        }
    }

    /// Construct from a generic source of known length. Not part of the crate's
    /// public API — the production entry point is [`Self::from_file`], and tests
    /// use `with_params` for custom fill/slice sizing.
    #[must_use]
    pub(crate) fn new(source: S, file_len: u64, read_streams: ReadStreams) -> Self {
        Self::from_parts(source, file_len, read_streams, DEFAULT_FILL_SIZE, MIN_SLICE_SIZE)
    }

    /// Construct with explicit fill/slice sizing (used by tests to exercise
    /// multi-fill / multi-slice behavior on small inputs).
    #[cfg(test)]
    #[must_use]
    pub(crate) fn with_params(
        source: S,
        file_len: u64,
        read_streams: ReadStreams,
        fill_size: usize,
        min_slice_size: usize,
    ) -> Self {
        Self::from_parts(source, file_len, read_streams, fill_size.max(1), min_slice_size.max(1))
    }

    /// Currently-resolved stream count (for tests asserting the probe outcome).
    #[cfg(test)]
    pub(crate) fn active_streams(&self) -> usize {
        self.active_streams
    }

    /// Force the `Auto` probe to resolve against a fixed throughput, making
    /// probe-outcome tests deterministic (no wall-clock dependence).
    #[cfg(test)]
    pub(crate) fn force_probe_throughput(&mut self, bytes_per_sec: f64) {
        self.forced_throughput = Some(bytes_per_sec);
    }

    /// Assemble the next fill window into `self.current`, recording any I/O
    /// error in `self.pending_error` (after any bytes read before it). Never
    /// returns an error itself: a fill that fails still serves whatever prefix
    /// it read, then the stored error surfaces on the following `read`.
    fn fill_next(&mut self) {
        if self.next_offset >= self.file_len {
            return;
        }
        let remaining = self.file_len - self.next_offset;
        // `remaining` may exceed usize only on a 32-bit target with a file
        // larger than 4 GiB; clamp to `fill_size` either way (a fill is never
        // larger than `fill_size`), so the cast cannot truncate meaningfully.
        let window = usize::try_from(remaining).unwrap_or(self.fill_size).min(self.fill_size);
        let ranges =
            slice_ranges(self.next_offset, window, self.active_streams, self.min_slice_size);

        let start = if self.probing { Some(Instant::now()) } else { None };
        let assembled: Result<Vec<u8>, (Vec<u8>, io::Error)> = if ranges.len() == 1 {
            let (off, len) = ranges[0];
            let mut buf = vec![0u8; len];
            match read_at_exact(&*self.source, &mut buf, off) {
                // A single-slice read has no successfully-read prefix to preserve.
                Ok(()) => Ok(buf),
                Err(e) => Err((Vec::new(), e)),
            }
        } else {
            self.scatter_fill(&ranges)
        };
        if let Some(t) = start {
            self.probe_bytes += window as u64;
            self.probe_elapsed += t.elapsed();
        }

        match assembled {
            Ok(buf) => {
                self.current = Some((buf, 0));
                self.next_offset += window as u64;
            }
            Err((prefix, err)) => {
                // Serve whatever contiguous prefix was read, then surface the error.
                self.next_offset += prefix.len() as u64;
                if prefix.is_empty() {
                    self.pending_error = Some(err);
                } else {
                    self.current = Some((prefix, 0));
                    self.pending_error = Some(err);
                }
            }
        }

        if self.probing {
            self.probe_fills_done += 1;
            if self.probe_fills_done >= AUTO_PROBE_FILLS {
                let secs = self.probe_elapsed.as_secs_f64();
                // Precision loss is irrelevant: this is a throughput estimate
                // over a few fills, fed straight into a ceil-and-cap.
                #[allow(clippy::cast_precision_loss)]
                let measured = if secs > 0.0 { self.probe_bytes as f64 / secs } else { 0.0 };
                let throughput = self.forced_throughput.unwrap_or(measured);
                self.active_streams = choose_stream_count(throughput);
                self.probing = false;
            }
        }
    }

    /// Fetch every range concurrently (one scoped thread per range) and
    /// reassemble in file order. On success returns the full window; on error
    /// returns `(contiguous prefix successfully read, first error)`.
    fn scatter_fill(&self, ranges: &[(u64, usize)]) -> Result<Vec<u8>, (Vec<u8>, io::Error)> {
        let results: Vec<io::Result<Vec<u8>>> = std::thread::scope(|scope| {
            let handles: Vec<_> = ranges
                .iter()
                .map(|&(off, len)| {
                    let source = Arc::clone(&self.source);
                    scope.spawn(move || {
                        let mut buf = vec![0u8; len];
                        read_at_exact(&*source, &mut buf, off)?;
                        Ok(buf)
                    })
                })
                .collect();
            handles
                .into_iter()
                .map(|h| {
                    h.join().unwrap_or_else(|_| {
                        Err(io::Error::other("scatter reader slice thread panicked"))
                    })
                })
                .collect()
        });

        // Concatenate slots in order until the first error. Slots are
        // contiguous by construction, so the successes before the first failure
        // form a valid forward prefix.
        let mut assembled = Vec::new();
        for result in results {
            match result {
                Ok(buf) => assembled.extend_from_slice(&buf),
                Err(e) => return Err((assembled, e)),
            }
        }
        Ok(assembled)
    }
}

#[cfg(unix)]
impl ScatterReader<std::fs::File> {
    /// Construct from an open regular file (production entry point).
    ///
    /// # Errors
    /// Returns an I/O error if the file's length cannot be read.
    pub fn from_file(file: std::fs::File, read_streams: ReadStreams) -> io::Result<Self> {
        let file_len = PositionalSource::byte_len(&file)?;
        Ok(Self::new(file, file_len, read_streams))
    }
}

impl<S: PositionalSource + 'static> Read for ScatterReader<S> {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        if out.is_empty() {
            return Ok(0);
        }
        loop {
            if let Some((buf, pos)) = &mut self.current {
                if *pos < buf.len() {
                    let n = (buf.len() - *pos).min(out.len());
                    out[..n].copy_from_slice(&buf[*pos..*pos + n]);
                    *pos += n;
                    return Ok(n);
                }
                self.current = None;
            }
            if let Some(err) = self.pending_error.take() {
                return Err(err);
            }
            if self.next_offset >= self.file_len {
                return Ok(0);
            }
            self.fill_next();
            if self.current.is_none() && self.pending_error.is_none() {
                return Ok(0);
            }
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
    use super::*;
    use proptest::prelude::*;
    use rstest::rstest;
    use std::io::Write;
    use std::thread::ThreadId;
    use tempfile::NamedTempFile;

    /// In-memory source of known bytes.
    struct MemSource(Vec<u8>);

    impl PositionalSource for MemSource {
        fn read_at(&self, buf: &mut [u8], offset: u64) -> io::Result<usize> {
            let off = offset as usize;
            if off >= self.0.len() {
                return Ok(0);
            }
            let n = (self.0.len() - off).min(buf.len());
            buf[..n].copy_from_slice(&self.0[off..off + n]);
            Ok(n)
        }
        fn byte_len(&self) -> io::Result<u64> {
            Ok(self.0.len() as u64)
        }
    }

    /// A source that errors on any read whose range covers `fail_at`.
    struct FaultySource {
        data: Vec<u8>,
        fail_at: u64,
    }

    impl PositionalSource for FaultySource {
        fn read_at(&self, buf: &mut [u8], offset: u64) -> io::Result<usize> {
            let end = offset + buf.len() as u64;
            if offset <= self.fail_at && self.fail_at < end {
                return Err(io::Error::other("injected fault"));
            }
            let off = offset as usize;
            if off >= self.data.len() {
                return Ok(0);
            }
            let n = (self.data.len() - off).min(buf.len());
            buf[..n].copy_from_slice(&self.data[off..off + n]);
            Ok(n)
        }
        fn byte_len(&self) -> io::Result<u64> {
            Ok(self.data.len() as u64)
        }
    }

    /// A source that panics if read from any thread other than the constructing
    /// one, and reports a controllable per-fill latency for the probe.
    struct InstrumentedSource {
        data: Vec<u8>,
        owner: ThreadId,
        per_read_sleep: Duration,
    }

    impl InstrumentedSource {
        fn new(data: Vec<u8>, per_read_sleep: Duration) -> Self {
            Self { data, owner: std::thread::current().id(), per_read_sleep }
        }
    }

    impl PositionalSource for InstrumentedSource {
        fn read_at(&self, buf: &mut [u8], offset: u64) -> io::Result<usize> {
            assert_eq!(
                std::thread::current().id(),
                self.owner,
                "read issued off the owning thread",
            );
            if !self.per_read_sleep.is_zero() {
                std::thread::sleep(self.per_read_sleep);
            }
            let off = offset as usize;
            if off >= self.data.len() {
                return Ok(0);
            }
            let n = (self.data.len() - off).min(buf.len());
            buf[..n].copy_from_slice(&self.data[off..off + n]);
            Ok(n)
        }
        fn byte_len(&self) -> io::Result<u64> {
            Ok(self.data.len() as u64)
        }
    }

    /// Thread-safe in-memory source with a controllable per-read latency (used
    /// by the probe scale-up test, which legitimately fans out across threads).
    struct TimedSource {
        data: Vec<u8>,
        per_read_sleep: Duration,
    }

    impl PositionalSource for TimedSource {
        fn read_at(&self, buf: &mut [u8], offset: u64) -> io::Result<usize> {
            if !self.per_read_sleep.is_zero() {
                std::thread::sleep(self.per_read_sleep);
            }
            let off = offset as usize;
            if off >= self.data.len() {
                return Ok(0);
            }
            let n = (self.data.len() - off).min(buf.len());
            buf[..n].copy_from_slice(&self.data[off..off + n]);
            Ok(n)
        }
        fn byte_len(&self) -> io::Result<u64> {
            Ok(self.data.len() as u64)
        }
    }

    /// A source whose reported `byte_len` overstates its backing data, so a read
    /// that trusts the length hits EOF early — the truncation / delete-race case
    /// that must fail closed rather than silently return a short clean stream.
    struct ShortSource {
        data: Vec<u8>,
        claimed_len: u64,
    }

    impl PositionalSource for ShortSource {
        fn read_at(&self, buf: &mut [u8], offset: u64) -> io::Result<usize> {
            let off = offset as usize;
            if off >= self.data.len() {
                return Ok(0);
            }
            let n = (self.data.len() - off).min(buf.len());
            buf[..n].copy_from_slice(&self.data[off..off + n]);
            Ok(n)
        }
        fn byte_len(&self) -> io::Result<u64> {
            Ok(self.claimed_len)
        }
    }

    /// Distinct-thread tracker shared between a [`ThreadTrackingSource`] and the
    /// test, so the test can inspect it after the source has been moved into (and
    /// dropped by) the reader.
    type ThreadSet = Arc<std::sync::Mutex<std::collections::HashSet<ThreadId>>>;

    /// Thread-safe source that records the distinct threads it was read from, so
    /// a test can positively confirm the scatter path actually fanned reads out
    /// across threads (rather than silently reading everything sequentially).
    struct ThreadTrackingSource {
        data: Vec<u8>,
        threads: ThreadSet,
    }

    impl ThreadTrackingSource {
        fn new(data: Vec<u8>) -> (Self, ThreadSet) {
            let threads: ThreadSet =
                Arc::new(std::sync::Mutex::new(std::collections::HashSet::new()));
            (Self { data, threads: Arc::clone(&threads) }, threads)
        }
    }

    impl PositionalSource for ThreadTrackingSource {
        fn read_at(&self, buf: &mut [u8], offset: u64) -> io::Result<usize> {
            self.threads.lock().expect("lock").insert(std::thread::current().id());
            let off = offset as usize;
            if off >= self.data.len() {
                return Ok(0);
            }
            let n = (self.data.len() - off).min(buf.len());
            buf[..n].copy_from_slice(&self.data[off..off + n]);
            Ok(n)
        }
        fn byte_len(&self) -> io::Result<u64> {
            Ok(self.data.len() as u64)
        }
    }

    fn drain<S: PositionalSource + 'static>(mut r: ScatterReader<S>) -> io::Result<Vec<u8>> {
        let mut out = Vec::new();
        r.read_to_end(&mut out)?;
        Ok(out)
    }

    // ---- choose_stream_count ----

    #[rstest]
    #[case::ebs_gp3_from_f13d1093(358_000_000.0, 4)]
    #[case::instance_store_from_9a80d7fe(1_920_000_000.0, 1)]
    #[case::just_under_target(1_500_000_000.0, 1)]
    #[case::very_slow_device_caps_at_max(100_000_000.0, MAX_STREAMS)]
    #[case::zero_is_defensive(0.0, 1)]
    #[case::negative_is_defensive(-5.0, 1)]
    fn choose_stream_count_matches_blueprint(#[case] measured: f64, #[case] expected: usize) {
        assert_eq!(choose_stream_count(measured), expected);
    }

    #[rstest]
    #[case::one(1, 1)]
    #[case::within_ceiling(4, 4)]
    #[case::at_ceiling(MAX_STREAMS, MAX_STREAMS)]
    #[case::above_ceiling_clamps(1000, MAX_STREAMS)]
    fn fixed_count_clamps_to_max_streams(#[case] requested: usize, #[case] expected: usize) {
        // An explicit `--read-streams n` is bounded by the same ceiling the Auto
        // heuristic uses; a huge value cannot spawn an unbounded fan-out per fill.
        let reader = ScatterReader::with_params(
            MemSource(vec![0u8; 8]),
            8,
            ReadStreams::Fixed(requested),
            32,
            4,
        );
        assert_eq!(reader.active_streams(), expected);
    }

    // ---- is_scatter_capable ----

    #[cfg(unix)]
    #[test]
    fn regular_file_is_scatter_capable() {
        let f = NamedTempFile::new().expect("temp");
        assert!(is_scatter_capable(f.as_file()));
    }

    #[cfg(unix)]
    #[test]
    fn char_device_is_not_scatter_capable() {
        let dev_null = std::fs::File::open("/dev/null").expect("open /dev/null");
        assert!(!is_scatter_capable(&dev_null));
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

    // ---- byte-identity (the load-bearing correctness gate) ----

    proptest! {
        #[test]
        fn byte_identical_to_source(
            data in proptest::collection::vec(any::<u8>(), 0..5000),
            streams in prop::sample::select(vec![1usize, 2, 4, 8]),
            fill in 1usize..64,
            min_slice in 1usize..16,
        ) {
            let len = data.len() as u64;
            let reader = ScatterReader::with_params(
                MemSource(data.clone()),
                len,
                ReadStreams::Fixed(streams),
                fill,
                min_slice,
            );
            let out = drain(reader).expect("drain");
            prop_assert_eq!(out, data);
        }

        #[test]
        fn byte_identical_under_auto(
            data in proptest::collection::vec(any::<u8>(), 0..5000),
            fill in 1usize..64,
            min_slice in 1usize..16,
        ) {
            let len = data.len() as u64;
            let reader = ScatterReader::with_params(
                MemSource(data.clone()),
                len,
                ReadStreams::Auto,
                fill,
                min_slice,
            );
            let out = drain(reader).expect("drain");
            prop_assert_eq!(out, data);
        }
    }

    // ---- EOF / short-file case table ----

    #[rstest]
    #[case::empty(0)]
    #[case::one(1)]
    #[case::just_under_fill(31)]
    #[case::exactly_fill(32)]
    #[case::just_over_fill(33)]
    #[case::several_fills(32 * 3 + 17)]
    fn round_trips_across_sizes_and_streams(
        #[case] len: usize,
        #[values(ReadStreams::Fixed(1), ReadStreams::Fixed(4), ReadStreams::Auto)]
        streams: ReadStreams,
    ) {
        let data: Vec<u8> = (0..len).map(|i| (i % 251) as u8).collect();
        let mut reader =
            ScatterReader::with_params(MemSource(data.clone()), len as u64, streams, 32, 4);
        let mut out = Vec::new();
        reader.read_to_end(&mut out).expect("read");
        assert_eq!(out, data);
        // A second read after EOF returns Ok(0) forever.
        let mut extra = [0u8; 8];
        assert_eq!(reader.read(&mut extra).expect("post-eof read"), 0);
        assert_eq!(reader.read(&mut extra).expect("post-eof read"), 0);
    }

    // ---- error propagation ----

    #[test]
    fn serves_prefix_then_surfaces_error() {
        // 200 bytes, fill 32, 4 streams: the first fill covers [0,32); put the
        // fault at offset 40 so fill 1 succeeds and fill 2 ([32,64)) fails.
        let data: Vec<u8> = (0..200).map(|i| (i % 251) as u8).collect();
        let mut reader = ScatterReader::with_params(
            FaultySource { data: data.clone(), fail_at: 40 },
            200,
            ReadStreams::Fixed(4),
            32,
            4,
        );
        let mut out = Vec::new();
        let mut buf = [0u8; 16];
        let mut err = None;
        loop {
            match reader.read(&mut buf) {
                Ok(0) => break,
                Ok(n) => out.extend_from_slice(&buf[..n]),
                Err(e) => {
                    err = Some(e);
                    break;
                }
            }
        }
        // Fill 1 ([0,32)) plus fill 2's contiguous prefix up to the failing
        // slice ([32,40)) are delivered before the error surfaces; the slice
        // covering offset 40 ([40,48)) is where it fails.
        assert_eq!(out, data[..40]);
        assert!(err.is_some(), "the injected fault must surface");
    }

    #[test]
    fn single_slice_error_surfaces() {
        // Fixed(1) always takes fill_next's single-slice branch (ranges.len() == 1).
        // A fault in the first fill has no successfully-read prefix, so the error
        // must surface on the first read with no bytes served.
        let data: Vec<u8> = (0..200).map(|i| (i % 251) as u8).collect();
        let mut reader = ScatterReader::with_params(
            FaultySource { data, fail_at: 0 },
            200,
            ReadStreams::Fixed(1),
            32,
            4,
        );
        let mut buf = [0u8; 16];
        let err = reader.read(&mut buf).expect_err("the single-slice fault must surface");
        assert_eq!(err.to_string(), "injected fault");
    }

    #[test]
    fn truncated_source_fails_closed_with_unexpected_eof() {
        // A source that claims more bytes than it holds (truncation / delete race)
        // must error, not silently feed a short-but-clean stream to the framer.
        let data: Vec<u8> = (0..50).map(|i| (i % 251) as u8).collect();
        let mut reader = ScatterReader::with_params(
            ShortSource { data, claimed_len: 200 },
            200,
            ReadStreams::Fixed(1),
            32,
            4,
        );
        let mut out = Vec::new();
        let err = reader.read_to_end(&mut out).expect_err("truncation must fail closed");
        assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
    }

    // ---- scatter positively engages under Fixed(n>1) ----

    #[test]
    fn fixed_four_fans_reads_across_multiple_threads() {
        // The parity tests prove output is identical across stream counts, which by
        // itself cannot distinguish "scatter engaged" from "flag silently ignored".
        // This positively confirms Fixed(4) drives reads off the calling thread.
        let len = 4096usize * 8;
        let data: Vec<u8> = (0..len).map(|i| (i % 251) as u8).collect();
        let (source, threads) = ThreadTrackingSource::new(data.clone());
        // A window (fill 4096) split 4 ways with min_slice 64 yields 4 slices → 4
        // scoped threads per fill.
        let reader =
            ScatterReader::with_params(source, len as u64, ReadStreams::Fixed(4), 4096, 64);
        let out = drain(reader).expect("drain");
        assert_eq!(out, data, "scatter output must still be byte-identical");
        let distinct = threads.lock().expect("lock").len();
        assert!(distinct > 1, "Fixed(4) must fan reads across threads; saw {distinct} distinct");
    }

    // ---- Fixed(1) never spawns a thread ----

    #[test]
    fn fixed_one_never_leaves_the_calling_thread() {
        let data: Vec<u8> = (0..500).map(|i| (i % 251) as u8).collect();
        let reader = ScatterReader::with_params(
            InstrumentedSource::new(data.clone(), Duration::ZERO),
            500,
            ReadStreams::Fixed(1),
            32,
            4,
        );
        // Would panic inside read_at if any read ran off the owning thread.
        let out = drain(reader).expect("drain");
        assert_eq!(out, data);
    }

    // ---- Auto probe transition ----

    #[test]
    fn auto_stays_at_one_until_the_probe_completes_then_scales() {
        // Forced throughput removes wall-clock flakiness. 358 MB/s (the gp3
        // measurement) resolves to 4 streams once the probe window elapses.
        let fill = 4096usize;
        let len = fill * (AUTO_PROBE_FILLS + 4);
        let data: Vec<u8> = (0..len).map(|i| (i % 251) as u8).collect();
        // TimedSource is cross-thread safe (scaling up fans reads out); zero
        // sleep since timing no longer matters.
        let mut reader = ScatterReader::with_params(
            TimedSource { data: data.clone(), per_read_sleep: Duration::ZERO },
            len as u64,
            ReadStreams::Auto,
            fill,
            64,
        );
        reader.force_probe_throughput(358_000_000.0);
        assert_eq!(reader.active_streams(), 1, "must start at one stream");

        let mut out = Vec::new();
        let mut buf = [0u8; 1024];
        loop {
            let n = reader.read(&mut buf).expect("read");
            if n == 0 {
                break;
            }
            out.extend_from_slice(&buf[..n]);
        }
        assert_eq!(out, data, "output must be byte-identical regardless of probing");
        assert_eq!(reader.active_streams(), 4, "358 MB/s must resolve to four streams");
    }

    #[test]
    fn auto_stays_at_one_when_device_is_already_fast() {
        // Forced 1.92 GB/s (the instance-store measurement) resolves to 1 stream,
        // so the owner-checking source never sees an off-thread read.
        let fill = 4096usize;
        let len = fill * (AUTO_PROBE_FILLS + 2);
        let data: Vec<u8> = (0..len).map(|i| (i % 251) as u8).collect();
        let mut reader = ScatterReader::with_params(
            InstrumentedSource::new(data.clone(), Duration::ZERO),
            len as u64,
            ReadStreams::Auto,
            fill,
            64,
        );
        reader.force_probe_throughput(1_920_000_000.0);
        // InstrumentedSource panics off-thread; if Auto wrongly scaled up on a
        // fast device, the post-probe fills would spawn threads and panic.
        let mut out = Vec::new();
        reader.read_to_end(&mut out).expect("read");
        assert_eq!(out, data);
        assert_eq!(reader.active_streams(), 1, "a fast device stays at one stream");
    }

    // ---- slice_ranges cover the window exactly ----

    #[rstest]
    #[case::single(100, 1, 10)]
    #[case::even_split(100, 4, 1)]
    #[case::remainder(103, 4, 1)]
    #[case::min_slice_limits_count(100, 8, 40)]
    fn slice_ranges_cover_the_window(
        #[case] window: usize,
        #[case] max_slices: usize,
        #[case] min_slice: usize,
    ) {
        let ranges = slice_ranges(1000, window, max_slices, min_slice);
        assert!(!ranges.is_empty());
        assert!(ranges.len() <= max_slices);
        // Contiguous, starting at 1000, covering exactly `window` bytes.
        let mut expected_off = 1000u64;
        let mut total = 0usize;
        for &(off, len) in &ranges {
            assert_eq!(off, expected_off);
            expected_off += len as u64;
            total += len;
        }
        assert_eq!(total, window);
    }
}

//! Scattered reads for the sort's input and spill files.
//!
//! Both phases read a file through one blocking sequential `read()` at a time:
//! phase 1 behind the exclusive `ReadInputBlocks` step, phase 2 behind each
//! spill file's reader mutex. That is queue depth 1, and a device serves one
//! stream far worse than several -- measured on the production volume, 358 MB/s
//! against 1177 MB/s for four concurrent streams, with buffer size making no
//! difference at all (2, 8 and 32 MiB all landed within 0.1%). A bigger single
//! request is not more concurrency.
//!
//! Measured consequence on `1kg-wgs-HG00096` (t16, 44 spill runs): the merge
//! read 52.6 GB of spill back at **313 MB/s**, costing 168.3s of worker-busy
//! time, and phase 1's read span was 141.5s. Raising `read_ahead_kb` to 4096
//! fixes both, but it is root-only and system-wide, so a tool cannot rely on
//! it. `posix_fadvise` was measured as the per-file substitute and does not
//! work -- `SEQUENTIAL` 359 MB/s and `WILLNEED` 348-350 MB/s against 360 with
//! no hint, `strace` confirming the calls succeed and the kernel simply
//! declines. That leaves doing the concurrency ourselves.
//!
//! # Why bytes and not blocks
//!
//! BGZF blocks and zstd frames are variable-length, so block *N*'s offset is
//! unknowable without scanning from block 0 -- which is why framing is serial.
//! It can stay that way: framing is ~4% of the reader's time. So a fill fetches
//! fixed-size *byte* slices at known offsets and this type serves them in file
//! order through a plain [`Read`]. The framers above are generic over `Read`
//! and do not change.
//!
//! # Why the slices run on the pool
//!
//! The obvious implementation gives each reader its own threads. That was
//! measured, and it works -- and it is wrong, because `--threads N` then buys
//! the user N pool workers plus an ingest thread plus however many readers are
//! live. On a 16 vCPU host the 4-stream configuration ran 21 threads and
//! doubled involuntary context switches (1.39M to 2.81M). It won 31% of wall
//! clock only because that box had five idle cores to lend.
//!
//! So a fill is offered to the pool as [`FetchSlice`] work any worker may
//! steal, and the thread budget stays exactly `--threads` + the ingest thread.
//! The filling thread reads the first slice itself and then **reclaims** every
//! slice no worker has started, rather than waiting on it: submitting and
//! blocking would deadlock the moment every worker is a filler waiting for
//! slices nobody is left to run. Reclaim makes the worst case "exactly the
//! single sequential read we do today", which is a slowdown and never a wedge.

use std::collections::VecDeque;
use std::fs::File;
use std::io::{self, BufReader, Read};
use std::os::unix::fs::FileExt;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::{Arc, Condvar, Mutex};

use crossbeam_queue::ArrayQueue;

/// Bytes one fill fetches, across all its slices.
///
/// Twice the 2 MiB buffer the sequential path uses, so that at four streams a
/// slice is still 1 MiB: the read ladder measured four 1 MiB streams at
/// 975 MB/s against 1081 MB/s at 4 MiB, so slices much below a megabyte start
/// giving the gain back.
pub(crate) const FILL_BYTES: usize = 4 * 1024 * 1024;

/// Most streams the ramp will ever reach.
///
/// Eight already matched four on the storage that needs concurrency at all
/// (1080 MB/s against 1081 on EBS gp3), and past this each extra stream is
/// pool work that buys nothing.
const MAX_STREAMS: usize = 8;

/// Fills measured at one stream before choosing a stream count.
///
/// Long enough that one slow fill cannot move the decision, short enough that
/// the probe costs 32 MiB of reading at whatever the device does unaided.
const AUTO_PROBE_FILLS: usize = 8;

/// Single-stream throughput at or above which concurrency has nothing to buy.
///
/// Derived from the consumer, not from the device: the ingest thread reads
/// 43.1 GB in ~53s, so it can absorb about 810 MB/s, and a reader that already
/// beats that is not what the sort is waiting on. This leaves ~1.5x headroom
/// over that floor.
///
/// The two devices measured sit either side of it by a wide margin -- EBS gp3
/// sustains 358 MB/s on one stream and wants four; a local instance-store SSD
/// sustains 2214 MB/s and wants one, where forcing eight measured 1.8%
/// *slower*. Anything from roughly 900 MB/s to 2 GB/s separates them.
const AUTO_TARGET_BYTES_PER_SEC: f64 = 1_200_000_000.0;

/// Streams to use, from what one stream measured.
///
/// `ceil(target / measured)`: enough streams to reach a rate the consumer
/// cannot outrun, and no more. Measuring the *device* rather than the pipeline
/// is what makes this work -- an earlier version compared fetch time against
/// the reader's own elapsed time, which is a high fraction on every device
/// because a reader mostly reads, and it duly ramped both gp3 and a local SSD
/// to the cap.
fn streams_for_measured_rate(bytes: u64, nanos: u64) -> usize {
    if bytes == 0 || nanos == 0 {
        return 1;
    }
    #[allow(clippy::cast_precision_loss)]
    let rate = bytes as f64 * 1e9 / nanos as f64;
    let wanted = (AUTO_TARGET_BYTES_PER_SEC / rate).ceil();
    if !wanted.is_finite() || wanted <= 1.0 {
        return 1;
    }
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let wanted = wanted as usize;
    wanted.min(MAX_STREAMS)
}

/// Smallest slice worth offering to another worker.
///
/// Below this the per-request overhead and the fill's own barrier cost more
/// than the concurrency buys, so a small fill uses fewer slices than the stream
/// count allows rather than splitting into pinpricks.
const MIN_SLICE_BYTES: usize = 512 * 1024;

/// One slice of a fill: read a byte range into the buffer waiting for it.
///
/// Cloned into the pool's queue as an `Arc` and also kept by the filling
/// thread, so whichever gets to it first runs it and the other no-ops.
pub(crate) struct FetchSlice {
    file: Arc<File>,
    offset: u64,
    index: usize,
    state: Arc<FillState>,
}

/// The parts of a fill its slices report back into.
struct FillState {
    progress: Mutex<FillProgress>,
    completed: Condvar,
}

/// Slice buffers and how many are still outstanding.
struct FillProgress {
    slots: Vec<SlotState>,
    remaining: usize,
    error: Option<io::Error>,
}

/// A slice's buffer, and who has it.
///
/// The transition out of `Pending` under the mutex *is* the claim: exactly one
/// thread can take the buffer, so a worker and the filling thread racing for
/// the same slice cannot both read it, and no separate flag is needed.
enum SlotState {
    /// Sized and waiting for someone to read into it.
    Pending(Vec<u8>),
    /// Claimed; a thread is in `read_exact_at` right now.
    InFlight,
    /// Read, and holding its bytes.
    Done(Vec<u8>),
}

/// Fill slices offered to the pool, with the census of what it took.
pub(crate) struct FetchQueue {
    jobs: ArrayQueue<Arc<FetchSlice>>,
    offered: AtomicU64,
    taken: AtomicU64,
    /// The stream count one uncontended probe measured, or 0 while none has.
    ///
    /// Shared because the thing being measured is the *device*, and a merge
    /// reads K spill files at once: each reader would otherwise measure its own
    /// share of a contended device, conclude it needs more streams, and add
    /// contention. Measured on EBS gp3 -- phase 1 alone saw 335 MB/s and chose
    /// four, while the 44 spill readers each saw 224-229 MB/s and chose six or
    /// seven, which took the merge from 119.3s to 129.1s.
    chosen_streams: AtomicUsize,
}

/// A [`Read`] over a file whose buffer is filled by concurrent positional reads.
pub(crate) struct ScatterReader {
    file: Arc<File>,
    len: u64,
    offset: u64,
    ready: VecDeque<Vec<u8>>,
    /// How far the consumer has read into `ready`'s front slice. A cursor
    /// rather than draining the front: the framers read a few bytes at a time,
    /// and shifting a megabyte down by thirteen bytes per call would cost more
    /// than the reads this saves.
    front_pos: usize,
    free: Vec<Vec<u8>>,
    /// Where slices are offered. `None` reads every fill on this thread, which
    /// is what the tests and any pool-less caller get.
    fetch: Option<Arc<FetchQueue>>,
    streams: usize,
    /// Where to book fill time, when this reader is the one the phase-1
    /// reader report describes. A fill is exactly what `TimedReader` measures
    /// on the sequential arm -- the fetch, as opposed to the framing above it.
    stats: Option<Arc<crate::phase1_stats::ReaderStats>>,
    /// Fills whose slices have been offered but whose bytes are not collected,
    /// oldest first and contiguous in the file.
    pending: VecDeque<PendingFill>,
    /// How many fills to keep in flight ahead of the consumer, so the device
    /// works while the framer frames. Zero for the merge's spill readers: the
    /// merge already sits at its consumer-serial floor, so all lookahead would
    /// buy there is `K * depth * FILL_BYTES` of extra buffers.
    lookahead: usize,
    /// Whether to grow `streams` from what the fills measure. See
    /// [`ramped_streams`].
    auto: bool,
    /// Bytes and nanoseconds the probe has measured at one stream so far.
    probe_bytes: u64,
    probe_nanos: u64,
    probe_fills: usize,
}

/// A fill in flight: its slices, where they report, and where it began.
///
/// Its length is not stored because it is always `offset - start`: only the
/// most recently issued fill is ever pending, and `issue` advances `offset` by
/// exactly that fill's size.
struct PendingFill {
    slices: Vec<Arc<FetchSlice>>,
    state: Arc<FillState>,
    start: u64,
}

/// How one file's bytes are read.
pub(crate) enum SpillSource {
    /// One sequential stream -- what both phases did before, and what
    /// `--read-streams 1` still does.
    Sequential(BufReader<File>),
    /// Several positional reads at once, run by pool workers.
    Scattered(ScatterReader),
}

impl FetchQueue {
    /// A queue sized for `workers` workers.
    ///
    /// Capacity only bounds how many slices can be *waiting*; a push that finds
    /// it full is not an error, because the filling thread reclaims anything the
    /// pool did not take.
    pub(crate) fn new(workers: usize) -> Arc<Self> {
        Arc::new(Self {
            jobs: ArrayQueue::new(workers.max(1) * 8),
            offered: AtomicU64::new(0),
            taken: AtomicU64::new(0),
            chosen_streams: AtomicUsize::new(0),
        })
    }

    /// Whether any slice is waiting. Drives the step's eligibility.
    pub(crate) fn is_empty(&self) -> bool {
        self.jobs.is_empty()
    }

    /// Run one offered slice. `false` means there was nothing left to do --
    /// either the queue was empty, or the filling thread had already reclaimed
    /// the slice this worker popped.
    pub(crate) fn run_one(&self) -> bool {
        let Some(slice) = self.jobs.pop() else {
            return false;
        };
        if slice.run() {
            self.taken.fetch_add(1, Ordering::Relaxed);
            return true;
        }
        false
    }

    /// `(slices offered, slices a worker actually ran)`.
    ///
    /// The census for "did this engage at all". Zero taken means every fill was
    /// read by the thread that wanted it, which is the pre-change behaviour and
    /// which no output check would ever notice.
    pub(crate) fn census(&self) -> (u64, u64) {
        (self.offered.load(Ordering::Relaxed), self.taken.load(Ordering::Relaxed))
    }

    /// The stream count a probe has already settled on, if any.
    fn settled_streams(&self) -> Option<usize> {
        match self.chosen_streams.load(Ordering::Relaxed) {
            0 => None,
            n => Some(n),
        }
    }

    /// Publish what one uncontended probe measured, if nobody has yet.
    fn settle_streams(&self, streams: usize) {
        let _ =
            self.chosen_streams.compare_exchange(0, streams, Ordering::Relaxed, Ordering::Relaxed);
    }

    /// Offer a slice, ignoring a full queue -- the filler reclaims either way.
    fn offer(&self, slice: &Arc<FetchSlice>) {
        if self.jobs.push(Arc::clone(slice)).is_ok() {
            self.offered.fetch_add(1, Ordering::Relaxed);
        }
    }
}

impl FetchSlice {
    /// Read this slice, unless someone else already claimed it.
    ///
    /// Returns whether this call did the read, which is what lets a worker tell
    /// "I helped" from "the filler beat me to it".
    fn run(&self) -> bool {
        let mut buf = {
            let mut progress = self.state.progress.lock().expect("fill state poisoned");
            match std::mem::replace(&mut progress.slots[self.index], SlotState::InFlight) {
                SlotState::Pending(buf) => buf,
                claimed => {
                    progress.slots[self.index] = claimed;
                    return false;
                }
            }
        };
        let result = self.file.read_exact_at(&mut buf, self.offset);
        let mut progress = self.state.progress.lock().expect("fill state poisoned");
        match result {
            Ok(()) => progress.slots[self.index] = SlotState::Done(buf),
            Err(e) => {
                if progress.error.is_none() {
                    progress.error = Some(e);
                }
            }
        }
        progress.remaining -= 1;
        let last = progress.remaining == 0;
        drop(progress);
        if last {
            self.state.completed.notify_all();
        }
        true
    }
}

impl ScatterReader {
    /// Read `file` from `offset`, offering slices to `fetch`.
    ///
    /// # Errors
    ///
    /// Returns an error if the file's length cannot be determined.
    pub(crate) fn new(
        file: File,
        offset: u64,
        streams: usize,
        fetch: Option<Arc<FetchQueue>>,
    ) -> io::Result<Self> {
        let len = file.metadata()?.len();
        Ok(Self {
            file: Arc::new(file),
            len,
            offset,
            ready: VecDeque::new(),
            front_pos: 0,
            free: Vec::new(),
            fetch,
            streams: streams.max(1),
            stats: None,
            pending: VecDeque::new(),
            lookahead: 0,
            auto: false,
            probe_bytes: 0,
            probe_nanos: 0,
            probe_fills: 0,
        })
    }

    /// Grow the stream count from what the fills measure, starting at one.
    ///
    /// Costs nothing before it has measured anything -- one stream is exactly
    /// the pre-change path -- and only grows when the consumer is demonstrably
    /// waiting on fills. See [`ramped_streams`].
    #[must_use]
    pub(crate) fn auto_tuned(mut self) -> Self {
        self.auto = true;
        self.streams = 1;
        self
    }

    /// Keep `fills` in flight ahead of the consumer.
    ///
    /// Without this a fill is issued only once the previous one is exhausted,
    /// so the device idles for as long as framing and decompression take.
    /// Measured on the phase-1 read span: 69.4s at depth 0, 66.9s at depth 1,
    /// and 62.3s for a reader that prefetched continuously with its own
    /// threads -- so depth is worth something well past one.
    ///
    /// Costs `(fills + 1) * FILL_BYTES` of buffers for this reader.
    #[must_use]
    pub(crate) fn looking_ahead(mut self, fills: usize) -> Self {
        self.lookahead = fills;
        self
    }

    /// Book this reader's fills into `stats`, so framing and fetching stay
    /// separable in the phase-1 reader report.
    #[must_use]
    pub(crate) fn timed(mut self, stats: Arc<crate::phase1_stats::ReaderStats>) -> Self {
        self.stats = Some(stats);
        self
    }

    /// Build a reader for a [`crate::external::ReadStreams`] setting, tuning
    /// itself when asked.
    ///
    /// # Errors
    ///
    /// Returns an error if the file's length cannot be determined.
    pub(crate) fn for_streams(
        file: File,
        offset: u64,
        streams: crate::external::ReadStreams,
        fetch: Option<Arc<FetchQueue>>,
    ) -> io::Result<Self> {
        let reader = Self::new(file, offset, streams.initial(), fetch)?;
        Ok(if streams.is_auto() { reader.auto_tuned() } else { reader })
    }

    /// A buffer of exactly `len` bytes, recycled where possible.
    ///
    /// Sized rather than cleared and refilled: `read_exact_at` overwrites every
    /// byte, so zero-filling a recycled buffer of the right length is pure cost
    /// over the gigabytes a sort reads.
    fn take_buf(&mut self, len: usize) -> Vec<u8> {
        let mut buf = self.free.pop().unwrap_or_default();
        if buf.len() != len {
            buf.resize(len, 0);
        }
        buf
    }

    /// Offer the next [`FILL_BYTES`] as slices and advance the fetch offset.
    ///
    /// Returns `None` at EOF. Nothing is read here beyond what a worker happens
    /// to pick up: the bytes are claimed in [`Self::collect`], which is what
    /// lets a fill be in flight while the caller does something else.
    fn issue(&mut self) -> Option<PendingFill> {
        let remaining = self.len.saturating_sub(self.offset);
        if remaining == 0 {
            return None;
        }
        let want = usize::try_from(remaining.min(FILL_BYTES as u64)).expect("a fill fits usize");
        // Never split below `MIN_SLICE_BYTES`: a fill near EOF is small, and
        // cutting it into pinpricks costs more than the concurrency returns.
        let slices = self.streams.min(want.div_ceil(MIN_SLICE_BYTES)).max(1);
        let slice_len = want.div_ceil(slices);
        let start = self.offset;

        let bufs: Vec<Vec<u8>> = (0..slices)
            .map(|index| self.take_buf(slice_len.min(want - index * slice_len)))
            .collect();
        let state = Arc::new(FillState {
            progress: Mutex::new(FillProgress {
                slots: bufs.into_iter().map(SlotState::Pending).collect(),
                remaining: slices,
                error: None,
            }),
            completed: Condvar::new(),
        });
        let slices: Vec<Arc<FetchSlice>> = (0..slices)
            .map(|index| {
                Arc::new(FetchSlice {
                    file: Arc::clone(&self.file),
                    offset: start + (index * slice_len) as u64,
                    index,
                    state: Arc::clone(&state),
                })
            })
            .collect();

        // Slice 0 stays for the collecting thread, which has nothing better to
        // do than read it; the rest go where a worker can reach them.
        if let Some(queue) = &self.fetch {
            for slice in &slices[1..] {
                queue.offer(slice);
            }
        }
        self.offset += want as u64;
        Some(PendingFill { slices, state, start })
    }

    /// Claim whatever the pool has not taken, wait for the rest, and publish the
    /// fill's bytes in file order.
    fn collect(&mut self, pending: &PendingFill) -> io::Result<()> {
        // Walk every slice: the first is this thread's own share, and the rest
        // are reclaimed if no worker has started them. `run` no-ops on a slice
        // someone else took, so this is both "do my share" and "take back what
        // nobody wanted" in one pass -- and it is why a full queue, an empty
        // pool, or a pool of workers all busy filling cannot wedge here.
        for slice in &pending.slices {
            slice.run();
        }

        let mut progress = pending.state.progress.lock().expect("fill state poisoned");
        // Only slices a worker claimed before we got to them remain, and that
        // worker is inside `read_exact_at` rather than waiting on anything, so
        // this wait is bounded by one disk read.
        while progress.remaining > 0 {
            progress = pending.state.completed.wait(progress).expect("fill state poisoned");
        }
        if let Some(e) = progress.error.take() {
            drop(progress);
            // `issue` moved the fetch offset before the bytes were claimed, so
            // put it back: a caller that reads again must retry this range, not
            // skip it. Anything issued after it covers bytes past the failure
            // and is dropped for the same reason.
            self.offset = pending.start;
            self.pending.clear();
            return Err(e);
        }
        for slot in &mut progress.slots {
            match std::mem::replace(slot, SlotState::InFlight) {
                SlotState::Done(buf) => self.ready.push_back(buf),
                _ => unreachable!("a fill with no error has every slice done"),
            }
        }
        Ok(())
    }

    /// Fold one single-stream fill into the probe, and decide once it has seen
    /// enough. Decides exactly once; after that the reader stops measuring.
    fn probe_fill(&mut self, bytes: u64, nanos: u64) {
        // Self-guarding rather than relying on the caller to stop: deciding
        // twice is exactly how the previous design went wrong, and a probe that
        // cannot be re-armed cannot repeat it. Later fills also run at the
        // chosen stream count, so they no longer measure one stream and would
        // answer a different question.
        if !self.auto {
            return;
        }
        // Adopt an answer someone else already measured rather than measuring a
        // contended device. In a merge this is every spill reader but the first.
        if let Some(settled) = self.fetch.as_ref().and_then(|q| q.settled_streams()) {
            self.streams = settled;
            self.auto = false;
            return;
        }
        self.probe_bytes = self.probe_bytes.saturating_add(bytes);
        self.probe_nanos = self.probe_nanos.saturating_add(nanos);
        self.probe_fills += 1;
        if self.probe_fills < AUTO_PROBE_FILLS {
            return;
        }
        let chosen = streams_for_measured_rate(self.probe_bytes, self.probe_nanos);
        #[allow(clippy::cast_precision_loss)]
        let mbps = self.probe_bytes as f64 * 1e3 / self.probe_nanos.max(1) as f64;
        log::debug!("read streams: one stream measured {mbps:.0} MB/s, using {chosen}");
        if let Some(queue) = &self.fetch {
            queue.settle_streams(chosen);
        }
        self.streams = chosen;
        self.auto = false;
    }

    /// Make bytes available in `ready`, leaving it empty only at EOF.
    fn fill(&mut self) -> io::Result<()> {
        if self.pending.is_empty() {
            let Some(first) = self.issue() else {
                return Ok(());
            };
            self.pending.push_back(first);
        }
        let pending = self.pending.pop_front().expect("just ensured non-empty");
        let started = self.auto.then(std::time::Instant::now);
        let fill_start = pending.start;
        self.collect(&pending)?;
        if let Some(t0) = started {
            self.probe_fill(
                self.offset.saturating_sub(fill_start),
                u64::try_from(t0.elapsed().as_nanos()).unwrap_or(u64::MAX),
            );
        }
        // Topped up after collecting rather than before, so the fills issued
        // here overlap the framing of what just landed rather than the wait for
        // it. At depth 0 this loop does nothing and the reader is demand-driven.
        while self.pending.len() < self.lookahead {
            let Some(next) = self.issue() else {
                break;
            };
            self.pending.push_back(next);
        }
        Ok(())
    }
}

impl Read for ScatterReader {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        if out.is_empty() {
            return Ok(0);
        }
        if self.ready.is_empty() {
            let before = self.offset;
            let sched_before =
                self.stats.as_ref().and_then(|_| crate::phase1_stats::thread_schedstat());
            let started = std::time::Instant::now();
            let outcome = self.fill();
            if let Some(stats) = &self.stats {
                let elapsed = u64::try_from(started.elapsed().as_nanos()).unwrap_or(u64::MAX);
                if let (Some(a), Some(b)) = (sched_before, crate::phase1_stats::thread_schedstat())
                {
                    stats.record_refill_sched(
                        b.0.saturating_sub(a.0),
                        b.1.saturating_sub(a.1),
                        b.2.saturating_sub(a.2),
                    );
                }
                // A failed fill still cost time; crediting it zero bytes keeps
                // the throughput figure honest.
                let fetched = usize::try_from(self.offset - before).unwrap_or(0);
                stats.record_refill(elapsed, fetched);
            }
            outcome?;
            if self.ready.is_empty() {
                return Ok(0);
            }
        }
        let front = self.ready.front().expect("just filled");
        let n = out.len().min(front.len() - self.front_pos);
        out[..n].copy_from_slice(&front[self.front_pos..self.front_pos + n]);
        self.front_pos += n;
        if self.front_pos == front.len() {
            let done = self.ready.pop_front().expect("just borrowed");
            self.free.push(done);
            self.front_pos = 0;
        }
        Ok(n)
    }
}

#[cfg(test)]
impl ScatterReader {
    /// Offset of the next byte the consumer will read.
    fn position(&self) -> u64 {
        let unconsumed: usize =
            self.ready.iter().map(Vec::len).sum::<usize>().saturating_sub(self.front_pos);
        // Pending fills are contiguous and end at the fetch offset, so the
        // oldest one's start is where the in-flight region begins.
        let in_flight = self.pending.front().map_or(0, |p| self.offset - p.start);
        self.offset.saturating_sub(unconsumed as u64).saturating_sub(in_flight)
    }
}

#[cfg(test)]
impl SpillSource {
    /// Offset of the next byte the consumer will read.
    ///
    /// Test-only: production never asks. The codec-detection tests assert that
    /// a zstd spill starts past `ZSPILL_MAGIC` and a BGZF one at byte 0, and
    /// that has to hold on whichever arm the stream count selects.
    pub(crate) fn position(&mut self) -> io::Result<u64> {
        match self {
            Self::Sequential(reader) => std::io::Seek::stream_position(reader),
            Self::Scattered(reader) => Ok(reader.position()),
        }
    }
}

impl Read for SpillSource {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        match self {
            Self::Sequential(reader) => reader.read(out),
            Self::Scattered(reader) => reader.read(out),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    /// A temp file of `len` pseudo-random bytes, plus the bytes themselves.
    fn fixture(len: usize) -> (tempfile::NamedTempFile, Vec<u8>) {
        let mut state = 0x2545_F491_4F6C_DD1Du64;
        let bytes: Vec<u8> = (0..len)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                u8::try_from((state >> 24) & 0xff).expect("masked to a byte")
            })
            .collect();
        let mut f = tempfile::NamedTempFile::new().expect("temp file");
        f.write_all(&bytes).expect("write");
        f.flush().expect("flush");
        (f, bytes)
    }

    fn read_all(reader: &mut impl Read, buf_size: usize) -> Vec<u8> {
        let mut out = Vec::new();
        let mut buf = vec![0u8; buf_size];
        loop {
            let n = reader.read(&mut buf).expect("read");
            if n == 0 {
                break;
            }
            out.extend_from_slice(&buf[..n]);
        }
        out
    }

    /// A reader whose slices are offered to a queue **nobody drains**, which is
    /// the reclaim path: every test below therefore exercises it by default.
    fn scattered(path: &std::path::Path, offset: u64, streams: usize) -> ScatterReader {
        let file = File::open(path).expect("open");
        ScatterReader::new(file, offset, streams, Some(FetchQueue::new(4))).expect("reader")
    }

    #[test]
    fn test_read_streams_parses_what_a_user_would_type() {
        use crate::external::ReadStreams;
        use std::str::FromStr;
        assert_eq!(ReadStreams::from_str("auto"), Ok(ReadStreams::Auto));
        assert_eq!(ReadStreams::from_str("AUTO"), Ok(ReadStreams::Auto), "case-insensitive");
        assert_eq!(ReadStreams::from_str("1"), Ok(ReadStreams::Fixed(1)));
        assert_eq!(ReadStreams::from_str("4"), Ok(ReadStreams::Fixed(4)));
        // Zero streams cannot read anything, and silently promoting it to one
        // would hide a typo behind a plausible-looking run.
        assert!(ReadStreams::from_str("0").is_err(), "zero is rejected, not rounded up");
        assert!(ReadStreams::from_str("four").is_err());
        assert!(ReadStreams::from_str("-1").is_err());
        assert!(ReadStreams::from_str("").is_err());
    }

    #[test]
    fn test_read_streams_round_trips_through_its_display() {
        // clap prints the default in `--help` via `Display`, so a value it
        // shows has to be one it accepts.
        use crate::external::ReadStreams;
        use std::str::FromStr;
        for value in [ReadStreams::Auto, ReadStreams::Fixed(1), ReadStreams::Fixed(8)] {
            assert_eq!(ReadStreams::from_str(&value.to_string()), Ok(value));
        }
    }

    /// Bytes and nanos for a device sustaining `mbps` on one stream.
    fn probe_of(mbps: u64) -> (u64, u64) {
        (mbps * 1_000_000, 1_000_000_000)
    }

    #[test]
    fn test_the_probe_picks_the_stream_count_each_measured_device_wanted() {
        // The two devices actually measured, at the rate one stream sustained
        // on each. gp3 gained 28% from four streams; the instance-store SSD
        // lost 1.8% from being pushed past one.
        let (bytes, nanos) = probe_of(358);
        assert_eq!(streams_for_measured_rate(bytes, nanos), 4, "EBS gp3 at 358 MB/s wants four");
        let (bytes, nanos) = probe_of(2214);
        assert_eq!(
            streams_for_measured_rate(bytes, nanos),
            1,
            "a local SSD at 2214 MB/s wants one"
        );
    }

    #[test]
    fn test_the_probe_scales_between_those_two_points() {
        // "Smarter than 1 or 4" is the point: the count is computed from the
        // measurement, so storage between the two measured devices gets a
        // count between the two answers.
        assert_eq!(streams_for_measured_rate(probe_of(700).0, probe_of(700).1), 2);
        assert_eq!(streams_for_measured_rate(probe_of(1199).0, probe_of(1199).1), 2);
        assert_eq!(streams_for_measured_rate(probe_of(1200).0, probe_of(1200).1), 1);
        assert_eq!(streams_for_measured_rate(probe_of(5000).0, probe_of(5000).1), 1);
    }

    #[test]
    fn test_the_probe_is_capped_however_slow_the_device() {
        // A very slow mount would ask for dozens of streams; past the cap they
        // are pool work that buys nothing.
        assert_eq!(streams_for_measured_rate(probe_of(10).0, probe_of(10).1), MAX_STREAMS);
        assert_eq!(streams_for_measured_rate(1, u64::MAX), MAX_STREAMS);
    }

    #[test]
    fn test_the_probe_refuses_to_divide_by_a_degenerate_measurement() {
        // A zero-length or zero-byte probe would produce an infinity or a NaN,
        // and a NaN comparison is false, which would silently pin the reader at
        // whatever it happened to be.
        assert_eq!(streams_for_measured_rate(0, 1_000), 1);
        assert_eq!(streams_for_measured_rate(1_000, 0), 1);
        assert_eq!(streams_for_measured_rate(0, 0), 1);
    }

    #[test]
    fn test_the_probe_decides_once_and_then_stops_measuring() {
        // Deciding repeatedly is what the previous design did, and it ramped
        // every device to the cap. One decision on a device-rate measurement
        // has no such failure mode -- but only if it really does stop.
        let (file, _) = fixture(1024);
        let handle = File::open(file.path()).expect("open");
        let mut reader = ScatterReader::new(handle, 0, 1, Some(FetchQueue::new(4)))
            .expect("reader")
            .auto_tuned();
        for _ in 0..AUTO_PROBE_FILLS {
            reader.probe_fill(FILL_BYTES as u64, 12_000_000);
        }
        assert!(!reader.auto, "the probe must disarm itself");
        assert_eq!(reader.streams, 4, "4 MiB in 12ms is ~350 MB/s, which wants four");
        let settled = reader.streams;
        for _ in 0..AUTO_PROBE_FILLS * 4 {
            reader.probe_fill(FILL_BYTES as u64, 1);
        }
        assert_eq!(reader.streams, settled, "a disarmed probe must not revisit its answer");
    }

    #[test]
    fn test_readers_sharing_a_queue_share_one_measurement() {
        // What is being measured is the device, and a merge reads K spill files
        // at once. Measured on EBS gp3: phase 1 alone saw 335 MB/s and chose
        // four, while 44 concurrent spill readers each saw 224-229 MB/s of a
        // device they were contending for and chose six or seven -- which took
        // the merge from 119.3s to 129.1s. A second reader must adopt, not
        // re-measure.
        let (file, _) = fixture(1024);
        let queue = FetchQueue::new(4);
        let reader_with = || {
            ScatterReader::new(
                File::open(file.path()).expect("open"),
                0,
                1,
                Some(Arc::clone(&queue)),
            )
            .expect("reader")
            .auto_tuned()
        };

        let mut first = reader_with();
        for _ in 0..AUTO_PROBE_FILLS {
            first.probe_fill(FILL_BYTES as u64, 12_000_000);
        }
        assert_eq!(first.streams, 4, "the uncontended probe chooses four");

        // The second reader is handed a *slow* measurement, as a contended one
        // would be. It must ignore it and take the settled answer.
        let mut second = reader_with();
        second.probe_fill(FILL_BYTES as u64, 900_000_000);
        assert_eq!(second.streams, 4, "adopted, not re-measured");
        assert!(!second.auto, "and disarmed on the spot");
        assert_eq!(second.probe_fills, 0, "without accumulating a probe of its own");
    }

    #[test]
    fn test_an_incomplete_probe_decides_nothing() {
        // Deciding early would let one slow fill choose the count for the run.
        let (file, _) = fixture(1024);
        let handle = File::open(file.path()).expect("open");
        let mut reader = ScatterReader::new(handle, 0, 1, Some(FetchQueue::new(4)))
            .expect("reader")
            .auto_tuned();
        for _ in 0..AUTO_PROBE_FILLS - 1 {
            reader.probe_fill(FILL_BYTES as u64, 900_000_000);
        }
        assert_eq!(reader.streams, 1, "a probe one fill short must not decide");
        assert!(reader.auto, "and must still be armed");
    }

    #[test]
    fn test_delivers_the_file_verbatim_at_every_stream_count() {
        // The contract that matters: splitting a fill must be invisible. A slice
        // delivered out of order corrupts every BGZF block or zstd frame the
        // sort then tries to parse.
        let (file, expected) = fixture(FILL_BYTES * 2 + 98_765);
        for streams in [1usize, 2, 4, 8] {
            let got = read_all(&mut scattered(file.path(), 0, streams), 64 * 1024);
            assert_eq!(got.len(), expected.len(), "length at {streams} streams");
            assert!(got == expected, "bytes differ at {streams} streams");
        }
    }

    #[test]
    fn test_a_fill_completes_when_no_worker_ever_takes_a_slice() {
        // Deadlock-freedom, stated directly. Every worker could be a filler
        // waiting on slices nobody is left to run, so a fill that only ever
        // *waits* for the pool wedges the sort. It must reclaim instead.
        let (file, expected) = fixture(FILL_BYTES * 2);
        let queue = FetchQueue::new(4);
        let file_handle = File::open(file.path()).expect("open");
        let mut reader =
            ScatterReader::new(file_handle, 0, 4, Some(Arc::clone(&queue))).expect("reader");

        let (tx, rx) = std::sync::mpsc::channel();
        std::thread::spawn(move || drop(tx.send(read_all(&mut reader, 64 * 1024))));
        let got = rx
            .recv_timeout(std::time::Duration::from_secs(30))
            .expect("the fill waited on a worker that never came");
        assert!(got == expected, "reclaimed slices must still deliver the file");
        let (offered, taken) = queue.census();
        assert!(offered > 0, "slices were offered");
        assert_eq!(taken, 0, "no worker was running, so none can have been taken");
    }

    #[test]
    fn test_exactly_one_caller_ever_runs_a_slice() {
        // The claim is a state transition under the fill's mutex, and this is
        // what it has to guarantee: a worker and the reclaiming filler racing
        // for the same slice must not both read it (a wasted read, and two
        // buffers where one is expected) nor both skip it (a fill that never
        // completes). Deterministic -- it asserts a count, not who won.
        let (file, expected) = fixture(64 * 1024);
        let state = Arc::new(FillState {
            progress: Mutex::new(FillProgress {
                slots: vec![SlotState::Pending(vec![0u8; expected.len()])],
                remaining: 1,
                error: None,
            }),
            completed: Condvar::new(),
        });
        let slice = Arc::new(FetchSlice {
            file: Arc::new(File::open(file.path()).expect("open")),
            offset: 0,
            index: 0,
            state: Arc::clone(&state),
        });

        let racers: Vec<_> = (0..8)
            .map(|_| {
                let slice = Arc::clone(&slice);
                std::thread::spawn(move || usize::from(slice.run()))
            })
            .collect();
        let winners: usize = racers.into_iter().map(|h| h.join().expect("racer")).sum();
        assert_eq!(winners, 1, "exactly one caller may read a given slice");

        let mut progress = state.progress.lock().expect("state");
        assert_eq!(progress.remaining, 0, "the single read completed the fill");
        match std::mem::replace(&mut progress.slots[0], SlotState::InFlight) {
            SlotState::Done(buf) => assert!(buf == expected, "the winner read the right bytes"),
            _ => panic!("the slice should be done"),
        }
    }

    #[test]
    fn test_a_worker_draining_the_queue_never_corrupts_the_bytes() {
        // A live pool worker racing the filler for every slice of every fill.
        // Deliberately asserts bytes only: whether the worker wins any given
        // race depends on how long slice 0 takes, which is a page-cache hit
        // here and a disk read in production. `test_exactly_one_caller...`
        // covers the race itself; this covers the result.
        let (file, expected) = fixture(FILL_BYTES * 6);
        let queue = FetchQueue::new(4);
        let stop = Arc::new(std::sync::atomic::AtomicBool::new(false));
        let worker = {
            let (queue, stop) = (Arc::clone(&queue), Arc::clone(&stop));
            std::thread::spawn(move || {
                while !stop.load(Ordering::Acquire) {
                    if !queue.run_one() {
                        std::thread::yield_now();
                    }
                }
            })
        };
        let file_handle = File::open(file.path()).expect("open");
        let mut reader =
            ScatterReader::new(file_handle, 0, 4, Some(Arc::clone(&queue))).expect("reader");
        let got = read_all(&mut reader, 64 * 1024);
        stop.store(true, Ordering::Release);
        worker.join().expect("worker");

        assert!(got == expected, "bytes must match however the slices were shared");
    }

    /// A reader that keeps one fill in flight, like the phase-1 input reader.
    fn ahead(path: &std::path::Path, streams: usize) -> ScatterReader {
        let file = File::open(path).expect("open");
        ScatterReader::new(file, 0, streams, Some(FetchQueue::new(4)))
            .expect("reader")
            .looking_ahead(1)
    }

    #[test]
    fn test_lookahead_delivers_the_file_verbatim_at_every_stream_count() {
        // Prefetching the next fill must not change a byte, and it is the case
        // most likely to: two fills are alive at once, so a buffer recycled or
        // an offset advanced at the wrong moment corrupts the seam between them.
        let (file, expected) = fixture(FILL_BYTES * 3 + 12_345);
        for streams in [2usize, 4, 8] {
            let got = read_all(&mut ahead(file.path(), streams), 64 * 1024);
            assert_eq!(got.len(), expected.len(), "length at {streams} streams");
            assert!(got == expected, "bytes differ at {streams} streams");
        }
    }

    #[test]
    fn test_lookahead_keeps_the_next_fill_in_flight_while_the_consumer_reads() {
        // The whole point: the device should be working on the next fill while
        // the framer is still consuming this one. Without it the disk idles for
        // exactly as long as framing takes, which is what the demand-driven
        // reader measured as 18s of unrecovered read span.
        let (file, _) = fixture(FILL_BYTES * 3);
        let mut reader = ahead(file.path(), 4);
        assert!(reader.pending.is_empty(), "nothing is in flight before the first read");
        assert!(reader.read(&mut [0u8; 64]).expect("first read") > 0, "the fixture is not empty");
        assert_eq!(reader.pending.len(), 1, "the next fill should already be issued");
    }

    #[test]
    fn test_lookahead_keeps_the_depth_it_was_given_in_flight() {
        // Depth is the knob that decides how much of the framing time the fetch
        // can hide behind, so it has to actually take effect.
        let (file, _) = fixture(FILL_BYTES * 6);
        let handle = File::open(file.path()).expect("open");
        let mut reader = ScatterReader::new(handle, 0, 4, Some(FetchQueue::new(4)))
            .expect("reader")
            .looking_ahead(3);
        assert!(reader.read(&mut [0u8; 64]).expect("read") > 0, "the fixture is not empty");
        assert_eq!(reader.pending.len(), 3, "three fills should be in flight");
    }

    #[test]
    fn test_lookahead_does_not_run_ahead_of_the_end_of_the_file() {
        // The last fill has no successor. Issuing one anyway would read past EOF
        // and turn a clean finish into an error.
        let (file, expected) = fixture(FILL_BYTES + 16);
        let mut reader = ahead(file.path(), 4);
        assert!(read_all(&mut reader, 1024) == expected);
        assert!(reader.pending.is_empty(), "nothing may be in flight at EOF");
    }

    #[test]
    fn test_a_failed_fill_rewinds_so_the_next_read_retries_the_same_bytes() {
        // A fill is issued before it is collected, so the fetch offset moves
        // first. If a failure left it moved, the retry would silently skip the
        // range that failed -- losing records rather than reporting an error.
        let (file, _) = fixture(FILL_BYTES * 3);
        let mut reader = ahead(file.path(), 4);
        let before = reader.position();
        file.as_file().set_len(0).expect("truncate");
        assert!(reader.read(&mut [0u8; 4096]).is_err(), "the read reports the failure");
        assert_eq!(reader.position(), before, "a failed fill must not advance the reader");
        assert!(reader.read(&mut [0u8; 4096]).is_err(), "and it keeps reporting it");
    }

    #[test]
    fn test_position_ignores_bytes_that_are_only_in_flight() {
        // `position` is where the consumer is, not how far ahead the fetch has
        // run. Counting an in-flight fill would make it jump forward by a whole
        // `FILL_BYTES` the moment lookahead was enabled.
        let (file, _) = fixture(FILL_BYTES * 3);
        let mut reader = ahead(file.path(), 4);
        let consumed = reader.read(&mut [0u8; 100]).expect("read");
        assert_eq!(consumed, 100, "a 100-byte read from a full fill returns 100");
        assert_eq!(reader.position(), 100, "only consumed bytes count");
    }

    #[test]
    fn test_reading_starts_at_the_requested_offset() {
        // A zstd spill's reader is positioned past `ZSPILL_MAGIC` rather than at
        // byte 0. Ignoring the offset would feed the frame parser four bytes of
        // magic it does not expect.
        let (file, expected) = fixture(FILL_BYTES + 4096);
        let got = read_all(&mut scattered(file.path(), 4, 4), 32 * 1024);
        assert!(got == expected[4..], "did not start at byte 4");
    }

    #[test]
    fn test_a_failed_fill_keeps_failing_instead_of_hanging() {
        // Fills are demand-driven, so truncating the file after the reader has
        // measured it makes the next fill read past EOF deterministically. A
        // caller may call `read` again after an error; it must not park.
        let (file, _) = fixture(FILL_BYTES * 2);
        let mut reader = scattered(file.path(), 0, 4);
        file.as_file().set_len(0).expect("truncate");
        assert!(reader.read(&mut [0u8; 4096]).is_err(), "the first read reports the failure");
        assert!(reader.read(&mut [0u8; 4096]).is_err(), "and so does every read after it");
    }

    #[test]
    fn test_an_empty_file_reads_as_eof_rather_than_hanging() {
        let (file, _) = fixture(0);
        assert!(read_all(&mut scattered(file.path(), 0, 4), 4096).is_empty());
    }

    #[test]
    fn test_tiny_reads_reassemble_across_slice_boundaries() {
        // The framers read a few bytes at a time, so slice hand-offs land
        // mid-call-sequence. Getting the advance wrong duplicates or skips bytes
        // exactly at the seam between two slices of one fill.
        let (file, expected) = fixture(FILL_BYTES + 7);
        assert!(read_all(&mut scattered(file.path(), 0, 3), 13) == expected);
    }

    #[test]
    fn test_one_stream_offers_nothing() {
        // `--read-streams 1` must stay exactly as cheap as it is today: no
        // slicing, no offer, no barrier -- just the read the caller wanted.
        let (file, expected) = fixture(FILL_BYTES * 2);
        let queue = FetchQueue::new(4);
        let mut reader = ScatterReader::new(
            File::open(file.path()).expect("open"),
            0,
            1,
            Some(Arc::clone(&queue)),
        )
        .expect("reader");
        assert!(read_all(&mut reader, 64 * 1024) == expected);
        assert_eq!(queue.census().0, 0, "a single stream offered work to the pool");
    }

    #[test]
    fn test_buffers_recycled_through_many_fills_deliver_the_file_verbatim() {
        // Seven fills' worth, so every fill after the first reuses buffers the
        // consumer handed back. Reusing one at a stale length corrupts the seam
        // between one slice and the next, and only once recycling happens.
        let (file, expected) = fixture(FILL_BYTES * 6 + 3);
        assert!(read_all(&mut scattered(file.path(), 0, 2), 128 * 1024) == expected);
    }

    /// Threads in this process, from the kernel's own accounting.
    ///
    /// Linux-only, which is where CI and every benchmark host run. `nextest`
    /// gives each test its own process, so this is a clean count rather than a
    /// number shared with a test harness.
    #[cfg(target_os = "linux")]
    fn thread_count() -> usize {
        std::fs::read_dir("/proc/self/task").expect("procfs").count()
    }

    #[cfg(target_os = "linux")]
    #[test]
    fn test_scattered_reading_spawns_no_threads_of_its_own() {
        // The whole reason this reader offers slices to the pool instead of
        // owning threads: `--threads N` must buy N workers and the ingest
        // thread, full stop. The previous design spawned N more per reader, and
        // nothing in this suite noticed -- it took a 16 vCPU host running 21
        // threads and double the involuntary context switches to surface it.
        let (file, expected) = fixture(FILL_BYTES * 3);
        let before = thread_count();
        let mut reader = scattered(file.path(), 0, 8);
        assert!(read_all(&mut reader, 64 * 1024) == expected, "and it still reads correctly");
        assert_eq!(thread_count(), before, "the reader spawned threads outside the pool");
    }

    #[test]
    fn test_a_scattered_read_matches_the_sequential_one_it_replaces() {
        // The two `SpillSource` arms are chosen by a flag, so they have to be
        // interchangeable byte for byte -- otherwise the flag changes results,
        // not just speed.
        let (file, _) = fixture(FILL_BYTES * 2 + 1234);
        let mut sequential = SpillSource::Sequential(BufReader::with_capacity(
            2 * 1024 * 1024,
            File::open(file.path()).expect("open"),
        ));
        let mut scatter = SpillSource::Scattered(scattered(file.path(), 0, 4));
        assert!(read_all(&mut sequential, 8192) == read_all(&mut scatter, 8192));
    }
}

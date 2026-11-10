//! Fan-out primitives: results of work that pool threads do on behalf of a
//! serial step return to it through an [`OrderedInbox`].
//!
//! A [`FanoutWorker`] runs a [`Job`] one dispatch at a time and re-contends
//! its phase cap permit for every dispatch. By default a dispatch is one
//! slice, but a job may run to completion in one dispatch
//! ([`Job::run_dispatch`]): the sort's range merge does (its claim ticket runs
//! a claimed range or assembly to completion), because a started job the
//! merge waits on in order must not be put down behind other work. Such a job
//! holds its thread and its cap permit until it is done, so steps sharing the
//! cap wait for it; its only early exit is the pipeline's done signal, which
//! it checks between slices.
//!
//! # Lost-wake protocol
//!
//! Same argument as `fgumi-sort`'s `MergeDemand`: [`OrderedInbox::await_next`]
//! stores `waiting` and THEN re-checks the slot, both under the inbox mutex;
//! [`OrderedInbox::deliver`] inserts under that mutex and THEN reads
//! `waiting`. The mutex orders the two critical sections, so either the
//! owner's re-check sees the item or the deliverer sees `waiting` and wakes
//! the owner (after releasing the lock). The loom model `inbox_await_model`
//! in `tests/loom_wake.rs` drives this pair, so splitting `await_next`'s
//! check and record into two critical sections fails it.

use ahash::AHashMap;
use std::collections::BTreeMap;
use std::io;
use std::sync::Arc;

// The inbox's lock and gauge come from `loom` under `--cfg loom`, so
// `tests/loom_wake.rs` checks the REAL `await_next`/`deliver` pair (the same
// arrangement as `admission.rs` and `fgumi-sort`'s `merge_demand.rs`).
#[cfg(loom)]
use loom::sync::atomic::{AtomicU64, Ordering};
#[cfg(not(loom))]
use parking_lot::Mutex;
#[cfg(not(loom))]
use std::sync::atomic::{AtomicU64, Ordering};

use crate::admission::{PhaseCap, StepCap, StepCapRef, admit_input};
use crate::item::HeapSize;
use crate::runtime::contexts::StepCounters;
use crate::runtime::wake::{SideWake, StepWaker};
use crate::signal::PipelineSignal;
use crate::step::{CounterSpec, SideChannel, Step, StepCtx, StepKind, StepOutcome, StepProfile};

/// `loom`'s mutex behind `parking_lot`'s `lock()` signature, so the inbox
/// code is the same under `--cfg loom`.
#[cfg(loom)]
struct Mutex<T>(loom::sync::Mutex<T>);

#[cfg(loom)]
impl<T> Mutex<T> {
    fn new(value: T) -> Self {
        Self(loom::sync::Mutex::new(value))
    }

    fn lock(&self) -> loom::sync::MutexGuard<'_, T> {
        self.0.lock().expect("inbox lock poisoned")
    }
}

struct State<R> {
    /// Delivered results, by ordinal.
    pending: BTreeMap<u64, R>,
    /// Bytes reserved for ordinals not yet delivered.
    reserved: AHashMap<u64, u64>,
    /// Bytes held by each delivered result.
    held: AHashMap<u64, u64>,
    /// The next ordinal the owner takes in ordered mode.
    next: u64,
    /// The owner declared (via `await_next`) that it is parked.
    waiting: bool,
}

/// A fan-in channel: worker steps [`deliver`](Self::deliver) results
/// `(ordinal, R)` from pool threads; the owning serial step takes them,
/// in ordinal order when `ordered`, and parks when the next one is absent.
///
/// The inbox itself accepts every reservation and delivery. Its bound is the
/// owner's obligation: the owner keeps at most a window of ordinals posted
/// and not yet taken, and that window bounds what the inbox holds.
///
/// Each ordinal is reserved (optionally) and delivered at most once; a
/// duplicate delivery, a reservation of a delivered ordinal, an ordinal below
/// the next one taken (ordered mode), or the take method of the other mode
/// panics, release builds included.
pub struct OrderedInbox<R> {
    ordered: bool,
    state: Mutex<State<R>>,
    owner: Mutex<Option<StepWaker>>,
    /// Reserved plus delivered-and-untaken bytes; lock-free to read.
    bytes: AtomicU64,
    #[cfg(test)]
    wake_attempts: AtomicU64,
}

impl<R: HeapSize + Send + 'static> OrderedInbox<R> {
    /// A new inbox; `ordered` makes [`try_take_next`](Self::try_take_next) yield
    /// ordinals consecutively from 0.
    #[must_use]
    pub fn new(ordered: bool) -> Arc<Self> {
        Arc::new(Self {
            ordered,
            state: Mutex::new(State {
                pending: BTreeMap::new(),
                reserved: AHashMap::new(),
                held: AHashMap::new(),
                next: 0,
                waiting: false,
            }),
            owner: Mutex::new(None),
            bytes: AtomicU64::new(0),
            #[cfg(test)]
            wake_attempts: AtomicU64::new(0),
        })
    }

    /// Binds the waker of the owning step, used to wake it on delivery.
    pub fn bind_owner(&self, waker: StepWaker) {
        *self.owner.lock() = Some(waker);
    }

    /// Counts `estimate` bytes against the gauge for an ordinal about to be
    /// worked on. Reserving an ordinal again before its delivery replaces the
    /// estimate.
    ///
    /// # Panics
    ///
    /// If `ordinal` was already delivered and not taken, or, in ordered mode,
    /// is below the next ordinal to take.
    pub fn reserve(&self, ordinal: u64, estimate: u64) {
        let mut s = self.state.lock();
        assert!(
            !self.ordered || ordinal >= s.next,
            "OrderedInbox: reserve for ordinal {ordinal} below next {}",
            s.next
        );
        assert!(
            !s.pending.contains_key(&ordinal),
            "OrderedInbox: reserve for ordinal {ordinal}, which is already delivered"
        );
        if let Some(old) = s.reserved.insert(ordinal, estimate) {
            self.bytes.fetch_sub(old, Ordering::Relaxed);
        }
        self.bytes.fetch_add(estimate, Ordering::Relaxed);
    }

    /// Delivers the result for `ordinal`, replacing its reservation by the
    /// result's heap size, and wakes the owner if it is parked on this
    /// delivery. Returns the wake outcome, or [`SideWake::Unavailable`] when no
    /// wake was attempted or no owner is bound.
    ///
    /// # Panics
    ///
    /// If `ordinal` was already delivered and not taken, or, in ordered mode,
    /// is below the next ordinal to take.
    pub fn deliver(&self, ordinal: u64, item: R) -> SideWake {
        let wake = {
            let mut s = self.state.lock();
            assert!(
                !self.ordered || ordinal >= s.next,
                "OrderedInbox: delivery for ordinal {ordinal} below next {}",
                s.next
            );
            if let Some(est) = s.reserved.remove(&ordinal) {
                self.bytes.fetch_sub(est, Ordering::Relaxed);
            }
            let size = item.heap_size() as u64;
            let was_empty = s.pending.is_empty();
            let duplicate = s.pending.insert(ordinal, item).is_some();
            assert!(!duplicate, "OrderedInbox: duplicate delivery for ordinal {ordinal}");
            s.held.insert(ordinal, size);
            self.bytes.fetch_add(size, Ordering::Relaxed);
            let wake = s.waiting && if self.ordered { ordinal == s.next } else { was_empty };
            if wake {
                s.waiting = false;
            }
            wake
        };
        if !wake {
            return SideWake::Unavailable;
        }
        #[cfg(test)]
        self.wake_attempts.fetch_add(1, Ordering::Relaxed);
        let waker = self.owner.lock().clone();
        waker.map_or(SideWake::Unavailable, |w| w.wake())
    }

    fn take(&self, s: &mut State<R>, ordinal: u64) -> Option<(u64, R)> {
        let item = s.pending.remove(&ordinal)?;
        if let Some(b) = s.held.remove(&ordinal) {
            self.bytes.fetch_sub(b, Ordering::Relaxed);
        }
        s.waiting = false;
        Some((ordinal, item))
    }

    /// Takes the result with the next ordinal, if delivered. Ordered inbox
    /// only: an unordered inbox never advances the internal `next` counter this
    /// reads, so use [`Self::try_take_any`] there.
    ///
    /// # Panics
    ///
    /// If the inbox is unordered.
    pub fn try_take_next(&self) -> Option<(u64, R)> {
        assert!(
            self.ordered,
            "OrderedInbox: try_take_next called on an unordered inbox; use try_take_any"
        );
        let mut s = self.state.lock();
        let n = s.next;
        let r = self.take(&mut s, n)?;
        s.next = n + 1;
        Some(r)
    }

    /// Whether a result is ready to take (the next ordinal in ordered mode,
    /// any in unordered mode). When `false`, records that the owner is
    /// waiting, so the delivery it needs wakes it; the check and the record
    /// are atomic.
    pub fn await_next(&self) -> bool {
        let mut s = self.state.lock();
        let ready =
            if self.ordered { s.pending.contains_key(&s.next) } else { !s.pending.is_empty() };
        s.waiting = !ready;
        ready
    }

    /// [`Self::await_next`] with its check and its record in two critical
    /// sections, the check first: the negative control of the loom inbox
    /// model, which must lose a wake (a delivery between the two sees no
    /// waiter, and the record then parks the owner for good).
    #[cfg(loom)]
    #[doc(hidden)]
    pub fn await_next_check_before_record(&self) -> bool {
        let ready = {
            let s = self.state.lock();
            if self.ordered { s.pending.contains_key(&s.next) } else { !s.pending.is_empty() }
        };
        self.state.lock().waiting = !ready;
        ready
    }

    /// Takes the lowest-ordinal delivered result regardless of order.
    /// Unordered inbox only: on an ordered inbox it would remove an ordinal
    /// without advancing `next`, leaving the owner waiting on a result that is
    /// gone (and never woken), so use [`Self::try_take_next`] there.
    ///
    /// # Panics
    ///
    /// If the inbox is ordered.
    pub fn try_take_any(&self) -> Option<(u64, R)> {
        assert!(
            !self.ordered,
            "OrderedInbox: try_take_any called on an ordered inbox; use try_take_next"
        );
        let mut s = self.state.lock();
        let k = *s.pending.keys().next()?;
        self.take(&mut s, k)
    }

    /// Delivered results not yet taken.
    #[must_use]
    pub fn pending(&self) -> usize {
        self.state.lock().pending.len()
    }

    /// Bytes reserved or held (lock-free).
    #[must_use]
    pub fn in_flight_bytes(&self) -> u64 {
        self.bytes.load(Ordering::Relaxed)
    }
}

impl<R: HeapSize + Send + 'static> SideChannel for OrderedInbox<R> {
    fn in_flight_bytes(&self) -> u64 {
        Self::in_flight_bytes(self)
    }

    fn describe(&self) -> String {
        let s = self.state.lock();
        format!("inbox(next={}, pending={})", s.next, s.pending.len())
    }
}

/// How much work one [`Job::run_slice`] call may do before yielding the thread.
#[derive(Clone, Copy, Debug)]
pub struct SliceBudget {
    /// Soft cap on bytes processed in one slice.
    pub bytes: u64,
    /// Soft cap on records processed in one slice.
    pub records: u64,
}

/// The result of one [`Job::run_slice`].
#[derive(Debug)]
pub enum JobStep<R> {
    /// The job has more work; run another slice.
    Continue,
    /// The job is finished and produced this result.
    Done(R),
    /// The job ran and has nothing to deliver (a claim ticket that found its
    /// work already taken). The worker counts the slice and delivers nothing.
    Spent,
}

/// The result of one [`Job::run_dispatch`]: what its last slice returned and
/// how many slices it ran.
#[derive(Debug)]
pub struct Dispatched<R> {
    /// What the dispatch's last slice returned ([`JobStep::Continue`] for a
    /// job stopped early, which the worker keeps).
    pub step: JobStep<R>,
    /// Slices the dispatch ran: 1 by default; any number for a job that runs
    /// to completion; 0 for a claim ticket spent before it ran one.
    pub slices: u64,
}

impl<R> Dispatched<R> {
    /// One slice's result, as the default [`Job::run_dispatch`] returns it.
    #[must_use]
    pub fn one(step: JobStep<R>) -> Self {
        Self { step, slices: 1 }
    }
}

/// A unit of work a serial step posts for pool threads, run in slices.
///
/// A [`FanoutWorker`] calls [`Self::run_dispatch`] once per dispatch, under
/// one cap permit. By default that runs one slice, so no dispatch holds a
/// thread or a permit for long. A job that overrides it to run to completion
/// holds its thread and permit until it is done (see the module docs), and
/// must check its `stop` callback between slices.
///
/// # Obligations
///
/// The deadlock monitor sees fan-out work only through byte counts, so an
/// implementation must keep it visible:
///
/// - A job and its [`Self::Output`] must report a non-zero
///   [`HeapSize::heap_size`]. A queued job with a zero size is invisible to
///   the job edge's budget and to the monitor's in-flight probe, and so is a
///   delivered result to the inbox gauge; a stall on either then classifies as
///   `Starving` instead of `Stalled`/`Wedged`.
/// - A blocking read (a `pread`, a file or pipe read) inside
///   [`Self::run_slice`] must be wrapped in [`crate::liveness::external_wait`],
///   so a slow device reads as an external wait, not a wedge.
/// - The owner must [`OrderedInbox::reserve`] each ordinal before posting its
///   job. A job a worker has popped and is running across several slices is
///   on no edge; until it delivers, the reservation is the only charge the
///   monitor can see.
/// - A job that returns [`JobStep::Spent`] pops an input item like any other
///   and delivers nothing; its ordinal (if any) must be delivered by whoever
///   ran that work instead.
pub trait Job: Send + HeapSize + 'static {
    /// The result delivered to the owner's [`OrderedInbox`]. Must report a
    /// non-zero [`HeapSize::heap_size`] (see the trait's obligations).
    type Output: Send + HeapSize + 'static;

    /// The job's ordinal: the position of its result in the owner's order.
    /// The worker reads it after each [`Self::run_slice`] (to deliver a
    /// result or wrap an error), so a job whose work is chosen when it first
    /// runs (a claim ticket) reports the claimed ordinal from then on.
    fn ordinal(&self) -> u64;

    /// Runs one slice of work within `budget`.
    ///
    /// # Errors
    ///
    /// Any I/O or processing error; the runtime fails the pipeline with it.
    fn run_slice(&mut self, budget: SliceBudget) -> io::Result<JobStep<Self::Output>>;

    /// Runs one dispatch of work: by default one [`Self::run_slice`].
    ///
    /// A job may override this to run several slices (to completion) in one
    /// dispatch. It must then call `stop` between slices and, once `stop`
    /// returns `true` (the pipeline failed or was cancelled), return
    /// [`JobStep::Continue`] at once: the worker keeps the job and the runtime
    /// ends the run instead of dispatching it again. It reports the slices it
    /// ran in [`Dispatched::slices`] (the worker's `slices` counter).
    ///
    /// # Errors
    ///
    /// As [`Self::run_slice`].
    fn run_dispatch(
        &mut self,
        budget: SliceBudget,
        stop: &dyn Fn() -> bool,
    ) -> io::Result<Dispatched<Self::Output>> {
        let _ = stop;
        self.run_slice(budget).map(Dispatched::one)
    }
}

/// The error a [`FanoutWorker`] fails with when a [`Job::run_slice`] fails:
/// the job's ordinal around the job's own error. The wrapping [`io::Error`]
/// keeps the job error's [`io::ErrorKind`]. The display carries the job
/// error's text (`job N: <error>`), because the pipeline reports a step error
/// by its display alone; so that a reporter walking the chain does not print
/// that text twice, [`source`](std::error::Error::source) skips the job
/// error and returns the job error's own source. The job error itself stays
/// reachable through [`io::Error::get_ref`] and a downcast to `JobError`. The
/// step name is not repeated here: the driver already reports which step
/// failed.
#[derive(Debug)]
pub struct JobError {
    ordinal: u64,
    source: io::Error,
}

impl JobError {
    /// Wraps `source`, the error of the job with ordinal `ordinal`. Public so
    /// an owner that runs one of its jobs itself reports its error as a
    /// worker would.
    #[must_use]
    pub fn wrap(ordinal: u64, source: io::Error) -> io::Error {
        io::Error::new(source.kind(), Self { ordinal, source })
    }

    /// The job's own error.
    #[must_use]
    pub fn job_error(&self) -> &io::Error {
        &self.source
    }

    /// The ordinal of the job that failed.
    #[must_use]
    pub fn ordinal(&self) -> u64 {
        self.ordinal
    }
}

impl std::fmt::Display for JobError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "job {}: {}", self.ordinal, self.source)
    }
}

impl std::error::Error for JobError {
    /// The job error's own source: the display already includes the job
    /// error's text.
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        self.source.source()
    }
}

/// Step counters. `jobs` counts jobs started (popped), including any that later
/// fail; `slices` counts completed slices (as each dispatch reports them in
/// [`Dispatched::slices`], so a job run to completion in one dispatch counts
/// every slice it ran); `output_bytes` the heap bytes of delivered results.
const COUNTERS: &[CounterSpec] = &[
    CounterSpec::new("jobs", "jobs"),
    CounterSpec::new("slices", "slices"),
    CounterSpec::new("output_bytes", "bytes"),
];
const JOBS: usize = 0;
const SLICES: usize = 1;
const OUTPUT_BYTES: usize = 2;

/// A `Parallel` step that runs the [`Job`]s of its input edge, one
/// [`Job::run_dispatch`] (by default one slice) per dispatch under its phase
/// cap, and delivers each result to an
/// [`OrderedInbox`]. Each worker instance holds its own in-progress job
/// between dispatches; a copy made by [`Step::new_worker_copy`] starts with
/// none.
///
/// The worker relies on the [`Job`] obligations: non-zero job and output
/// heap sizes, blocking reads inside `run_slice` wrapped in
/// [`crate::liveness::external_wait`], and an owner that reserves each
/// ordinal on the inbox before posting it, since an in-progress job held by
/// a worker instance between slices is visible to the deadlock monitor only through
/// that reservation.
pub struct FanoutWorker<J: Job> {
    name: &'static str,
    inbox: Arc<OrderedInbox<J::Output>>,
    budget: SliceBudget,
    cap: Option<StepCap>,
    current: Option<J>,
}

impl<J: Job> FanoutWorker<J> {
    /// A worker named `name` delivering to `inbox`, running slices of `budget`.
    #[must_use]
    pub fn new(
        name: &'static str,
        inbox: Arc<OrderedInbox<J::Output>>,
        budget: SliceBudget,
    ) -> Self {
        Self { name, inbox, budget, cap: None, current: None }
    }

    /// Admits slices through `cap` (none by default).
    #[must_use]
    pub fn with_step_cap(mut self, cap: impl Into<Option<StepCap>>) -> Self {
        self.cap = cap.into();
        self
    }

    /// Runs one dispatch of `job`: keeps it for the next dispatch on
    /// `Continue`, delivers its result on `Done`, drops it on `Spent`. The
    /// ordinal is read after the dispatch (see [`Job::ordinal`]). `signal`
    /// is the pipeline's, polled by a job that runs several slices.
    fn advance(
        &mut self,
        mut job: J,
        counters: &StepCounters,
        signal: &PipelineSignal,
    ) -> io::Result<()> {
        let pipeline_done = || signal.is_done();
        let Dispatched { step, slices } = job
            .run_dispatch(self.budget, &pipeline_done)
            .map_err(|e| JobError::wrap(job.ordinal(), e))?;
        counters.add(SLICES, slices);
        match step {
            JobStep::Continue => self.current = Some(job),
            JobStep::Done(out) => {
                counters.add(OUTPUT_BYTES, out.heap_size() as u64);
                let _ = self.inbox.deliver(job.ordinal(), out);
            }
            JobStep::Spent => {}
        }
        Ok(())
    }

    fn admitted(
        &mut self,
        ctx: &mut StepCtx<'_, Self>,
        cap: Option<&PhaseCap>,
    ) -> io::Result<StepOutcome> {
        // A job in progress is work already claimed: take a permit directly.
        // Otherwise admit_input checks the edge before a permit is taken.
        let (_permit, job, first) = if let Some(job) = self.current.take() {
            let permit = match cap {
                Some(cap) => {
                    let Some(p) = cap.try_acquire() else {
                        self.current = Some(job);
                        return Ok(StepOutcome::Capped);
                    };
                    Some(p)
                }
                None => None,
            };
            (permit, job, false)
        } else {
            let permit = match admit_input(ctx.input, cap) {
                Ok(p) => p,
                Err(outcome) => return Ok(outcome),
            };
            match ctx.input.pop() {
                Some(job) => (permit, job, true),
                None => {
                    // Uncapped: admit_input did not look at the edge.
                    return Ok(if ctx.input.is_drained() {
                        StepOutcome::Finished
                    } else {
                        ctx.input.note_empty_poll();
                        StepOutcome::NoProgress
                    });
                }
            }
        };
        if first {
            ctx.counters.add(JOBS, 1);
        }
        self.advance(job, ctx.counters, ctx.signal)?;
        Ok(StepOutcome::Progress)
    }
}

impl<J: Job> Step for FanoutWorker<J> {
    type Input = J;
    type Outputs = ();

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: self.name,
            kind: StepKind::Parallel,
            sticky: false,
            output_queues: vec![],
            branch_ordering: vec![],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        let cap = self.cap.take();
        let outcome = self.admitted(ctx, cap.as_ref().and_then(StepCap::active));
        self.cap = cap;
        outcome
    }

    fn new_worker_copy(&self) -> Self {
        Self {
            name: self.name,
            inbox: Arc::clone(&self.inbox),
            budget: self.budget,
            cap: self.cap.clone(),
            current: None,
        }
    }

    fn counters(&self) -> &'static [CounterSpec] {
        COUNTERS
    }

    fn step_cap(&self) -> Option<StepCapRef<'_>> {
        self.cap.as_ref().map(StepCap::view)
    }

    fn side_channels(&self) -> Vec<Arc<dyn SideChannel>> {
        vec![Arc::clone(&self.inbox) as Arc<dyn SideChannel>]
    }
}

#[cfg(test)]
mod tests;

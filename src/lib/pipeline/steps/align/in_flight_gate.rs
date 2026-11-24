//! The subprocess align backend's in-flight gate: the byte budget bounding the
//! unmapped reads held while the aligner works, derived from the aligner's chunk
//! size, and the gate that enforces it (growing past it when the aligner needs
//! more input before its first output, up to a ceiling).

use std::io;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{Duration, Instant};

use bytesize::ByteSize;
use parking_lot::Mutex;

use super::IN_FLIGHT_BYTES_PER_BASE;

/// Multiplier over a single `-K` chunk for the in-flight budget. bwa's
/// `kt_pipeline` double-buffers (`p_nt=2`: one chunk aligning, one
/// queued), so the aligner can hold ~2 chunks before emitting; 4× adds
/// slack for the FASTQ-write buffer and pipeline jitter so the writer
/// rarely blocks before the aligner has a full chunk to emit (when it
/// does, the gate grows after a grace period — see `InFlightGate`).
const IN_FLIGHT_CHUNK_MULTIPLIER: u64 = 4;

/// Floor for the in-flight-unmapped budget. Keeps the budget workable
/// for small `-K` values (custom aligner commands, tests) where the
/// `-K`-derived value would be tiny, and covers aligners whose internal
/// buffering exceeds their nominal `-K`.
const IN_FLIGHT_MIN_BUDGET: u64 = 512 * 1024 * 1024;

/// Derive the in-flight-unmapped byte budget for the subprocess align path from
/// the aligner's `-K` chunk size (bases per batch). The budget bounds the
/// otherwise unbounded writer→reader token backlog: the writer blocks once the
/// in-flight unmapped bytes reach this, so a fast-draining aligner can't
/// accumulate the whole input's unmapped reads in RAM (issue #382).
///
/// Sized `≥` the aligner's real buffering so a **streaming** aligner
/// (one that emits output after ≤ `-K` bases — every real aligner)
/// always has a full chunk to emit before the writer blocks, so the
/// block is transient. An aligner that reads further ahead before its
/// first output (or a preset on short templates, whose records cost more
/// heap per base than [`IN_FLIGHT_BYTES_PER_BASE`]) overflows it; the gate
/// then grows up to the read-ahead ceiling rather than stall (see
/// [`InFlightGate`]). For the default `-K` of 150M bases this is ~1.8 GiB.
#[must_use]
pub(crate) fn in_flight_budget_for_chunk_size(chunk_size_bases: u64) -> u64 {
    chunk_size_bases
        .saturating_mul(IN_FLIGHT_BYTES_PER_BASE)
        .saturating_mul(IN_FLIGHT_CHUNK_MULTIPLIER)
        .max(IN_FLIGHT_MIN_BUDGET)
}

/// Hidden diagnostic: an environment variable that replaces the subprocess align
/// path's derived in-flight budget with an exact byte count, bypassing the
/// [`IN_FLIGHT_MIN_BUDGET`] floor. Tests use it to make a small input overflow
/// the budget; it is not a tuning knob (`--aligner::chunk-size` is).
pub(crate) const IN_FLIGHT_BUDGET_ENV: &str = "FGUMI_ALIGN_IN_FLIGHT_BUDGET";

/// Parse the value of [`IN_FLIGHT_BUDGET_ENV`]: a positive byte count overrides
/// the derived budget; unset means no override. Any other value is ignored
/// with a warning rather than failing the run, as for the other hidden
/// `FGUMI_*` diagnostics.
#[must_use]
pub(crate) fn in_flight_budget_override(value: Option<&str>) -> Option<u64> {
    let value = value?;
    match value.trim().parse::<u64>() {
        Ok(bytes) if bytes > 0 => Some(bytes),
        _ => {
            log::warn!(
                "{IN_FLIGHT_BUDGET_ENV}={value:?} is ignored: set it to a positive number of bytes"
            );
            None
        }
    }
}

/// How long the writer waits on a full in-flight budget, with no aligner output
/// in that time, before it concludes that the aligner needs more input before
/// it can emit anything and admits the batch anyway (see [`InFlightGate`]).
///
/// A wrong guess (an aligner that is merely busy) costs only memory: the input
/// admitted past the budget is what the aligner itself goes on to read, since
/// the stdin pipe still throttles the writer. So the grace can be short. 5 s is
/// also well inside the pipeline deadlock watchdog's default 10 s warning.
pub(crate) const IN_FLIGHT_STALL_GRACE: Duration = Duration::from_secs(5);

/// How long a batch waiting at the read-ahead ceiling may see no aligner output
/// consumed before the run fails (see [`InFlightGate`]). Far longer than the
/// grace: a failure here ends the run, so it must not fire for an aligner that
/// is merely slow per chunk; a real deadlock is still reported rather than
/// hanging forever.
pub(crate) const IN_FLIGHT_CEILING_TIMEOUT: Duration = Duration::from_mins(5);

/// Hidden diagnostic: an environment variable that replaces
/// [`IN_FLIGHT_CEILING_TIMEOUT`] with a number of seconds, so tests can reach the
/// failure at the ceiling without waiting five minutes. Not a tuning knob.
pub(crate) const IN_FLIGHT_CEILING_TIMEOUT_ENV: &str = "FGUMI_ALIGN_CEILING_TIMEOUT_SECS";

/// Parse the value of [`IN_FLIGHT_CEILING_TIMEOUT_ENV`]: a positive number of
/// seconds overrides the timeout; unset means no override. Any other value is
/// ignored with a warning, as for [`IN_FLIGHT_BUDGET_ENV`].
#[must_use]
pub(crate) fn ceiling_timeout_override(value: Option<&str>) -> Option<Duration> {
    let value = value?;
    match value.trim().parse::<u64>() {
        Ok(secs) if secs > 0 => Some(Duration::from_secs(secs)),
        _ => {
            log::warn!(
                "{IN_FLIGHT_CEILING_TIMEOUT_ENV}={value:?} is ignored: set it to a positive \
                 number of seconds"
            );
            None
        }
    }
}

/// The default read-ahead ceiling is this fraction (1 / N) of host memory. The
/// in-flight unmapped reads are resident in fgumi beside the aligner's own index
/// and buffers and the rest of the pipeline's queues and sort memory, so a
/// quarter of RAM leaves room for those while covering the read-ahead
/// measured for aligners that outgrow the derived budget (about 3× the budget
/// for 20 bp templates under a bwa preset; about 1 GiB of FASTQ for a
/// many-threaded aligner).
const READ_AHEAD_CEILING_HOST_MEMORY_DIVISOR: u64 = 4;

/// The hard ceiling the in-flight budget may grow to: `requested`
/// (`--aligner::max-read-ahead`) when set, else a quarter of `host_memory`
/// (see [`READ_AHEAD_CEILING_HOST_MEMORY_DIVISOR`]), and never below `budget`
/// by default, so the default never tightens the bound the budget already set.
/// An unknown host memory (`0`) leaves the ceiling at the budget: the gate then
/// fails a stall after the grace period instead of growing.
#[must_use]
pub(crate) fn read_ahead_ceiling(budget: u64, requested: Option<u64>, host_memory: u64) -> u64 {
    requested.unwrap_or_else(|| (host_memory / READ_AHEAD_CEILING_HOST_MEMORY_DIVISOR).max(budget))
}

/// Byte-budget gate bounding the subprocess align path's in-flight unmapped
/// reads (the writer→reader token backlog). The writer
/// [`acquire`](InFlightGate::acquire)s before feeding a batch to the aligner and
/// the reader [`release`](InFlightGate::release)s when it consumes the matching
/// token, so resident in-flight unmapped bytes stay near the budget.
///
/// The writer blocks only while in-flight is non-empty AND the batch would take
/// it past the budget — a single oversized batch always passes when the gate is
/// empty. A streaming aligner emits output within its `-K` buffering, which the
/// budget is sized to cover (see [`in_flight_budget_for_chunk_size`]), so the
/// reader drains and the block lifts.
///
/// An aligner may instead need more input than the budget before its first
/// output byte (a many-threaded aligner reading one buffer per thread, or a bwa
/// preset on short templates, whose records cost more heap per base than the
/// budget assumes). Blocking would then deadlock: the writer waits for output
/// and the aligner waits for input. So the gate never blocks indefinitely:
///
/// - **Grow.** When the writer has waited `stall_grace` on a full budget with no
///   reader progress in that time, the batch is admitted and the budget grows
///   to at least twice its size (capped at `ceiling`), with a warning the first
///   time. Until the next progress the gate is *stalled* and admits further
///   batches at once, growing as needed, so feeding the aligner costs one grace
///   period per stall, not one per batch. Progress ends the stall; the grown
///   budget stays. A wrong guess (an aligner that is merely slow) costs only
///   memory, so the grace is short.
/// - **Fail.** A batch that would take in-flight past `ceiling` waits for a
///   release. Only after `fail_after` without reader progress (far longer than
///   the grace, since failing a slow-but-working aligner is costly) does it fail,
///   with an error naming the ceiling and the options that move it.
///
/// Each consumed record signals [`record_progress`](InFlightGate::record_progress),
/// even within a large template or token; only a token release frees bytes.
///
/// No aligner output is consumed while the reader is blocked handing a batch
/// downstream, so that time does not count as quiet
/// ([`downstream_wait`](InFlightGate::downstream_wait)): backpressure from a
/// slow downstream stage is not mistaken for an aligner that needs input.
///
/// Growing only admits what the aligner reads: past the budget the writer is
/// still throttled by the aligner's stdin pipe. Pairing of aligner output to
/// unmapped templates is by token order, so admitting more batches does not
/// change the output.
///
/// If the reader exits (EOF or error), it latches `consumer_gone` so a blocked
/// writer wakes and bails rather than hanging.
pub(crate) struct InFlightGate {
    /// Hard limit on in-flight bytes once the gate holds more than one batch.
    ceiling: u64,
    /// How long a full gate may see no reader progress before it grows.
    stall_grace: Duration,
    /// How long a gate at its ceiling may see no reader progress before it fails.
    fail_after: Duration,
    /// Sampled by a waiting writer; record reads need neither the lock nor a wake.
    progress: AtomicU64,
    inner: Mutex<GateInner>,
    cond: parking_lot::Condvar,
}

// Each flag is an independent latch read under the one lock (reader position,
// stall, the two one-time warnings, consumer gone); a state enum would have to
// encode their product.
#[allow(clippy::struct_excessive_bools)]
struct GateInner {
    in_flight: u64,
    /// The current budget: the configured one, raised by stalls (≤ `ceiling`).
    budget: u64,
    /// Last progress sampled by the writer, including across acquire calls.
    progress_seen: u64,
    /// Whether the reader is blocked handing a batch downstream.
    downstream_blocked: bool,
    /// Whether the gate grew on a stall and has seen no reader progress since.
    stalled: bool,
    /// Whether the one-time growth warning was logged.
    warned_growth: bool,
    /// Whether the one-time waiting-at-the-ceiling warning was logged.
    warned_ceiling: bool,
    consumer_gone: bool,
}

/// What [`InFlightGate::acquire`] logs once its lock is released.
enum GateEvent {
    /// The budget grew on a stall.
    Grew { from: u64, to: u64, in_flight: u64, first: bool },
    /// A batch began waiting at the ceiling (logged once per gate).
    AtCeiling { in_flight: u64 },
}

impl InFlightGate {
    /// A gate starting at `budget` bytes that may grow to `ceiling` after a
    /// full gate sees no progress for `stall_grace`, and fails a batch past the
    /// ceiling after `fail_after` with no progress. A `budget` above `ceiling`
    /// is lowered to it.
    pub(crate) fn new(
        budget: u64,
        ceiling: u64,
        stall_grace: Duration,
        fail_after: Duration,
    ) -> Self {
        Self {
            ceiling,
            stall_grace,
            fail_after,
            progress: AtomicU64::new(0),
            inner: Mutex::new(GateInner {
                in_flight: 0,
                budget: budget.min(ceiling),
                progress_seen: 0,
                downstream_blocked: false,
                stalled: false,
                warned_growth: false,
                warned_ceiling: false,
                consumer_gone: false,
            }),
            cond: parking_lot::Condvar::new(),
        }
    }

    /// Reserve `n` in-flight bytes, blocking while the gate is full and
    /// non-empty, for at most about `stall_grace` without reader progress, or
    /// `fail_after` at the ceiling (see the type docs). Returns `Ok(false)` if
    /// the consumer (reader) has gone — the caller (writer) should then stop.
    ///
    /// # Errors
    ///
    /// Returns an error when admitting `n` would take in-flight bytes past the
    /// ceiling and the reader made no progress for `fail_after`. The batch is
    /// not admitted then.
    pub(crate) fn acquire(&self, n: u64) -> io::Result<bool> {
        let mut g = self.inner.lock();
        // Started on the first wait only, so the uncontended path reads no clock.
        let mut quiet_since: Option<Instant> = None;
        let event = loop {
            if g.consumer_gone {
                return Ok(false);
            }
            let wanted = g.in_flight.saturating_add(n);
            // The `in_flight == 0` arm lets a single oversized batch through
            // when the gate is empty: it can't be split, so holding it is
            // unavoidable, and blocking on it would deadlock.
            if g.in_flight == 0 || wanted <= g.budget {
                break None;
            }
            let now = Instant::now();
            let progress = self.progress.load(Ordering::Relaxed);
            let made_progress = progress != g.progress_seen;
            if made_progress {
                g.stalled = false;
                g.progress_seen = progress;
            }
            if made_progress || g.downstream_blocked || quiet_since.is_none() {
                // The reader made progress (or is blocked downstream, which is
                // not an aligner stall): restart the quiet clock.
                quiet_since = Some(now);
            }
            let quiet = now.duration_since(quiet_since.unwrap_or(now));
            if wanted <= self.ceiling {
                if g.stalled || quiet >= self.stall_grace {
                    let from = g.budget;
                    g.budget = wanted.max(from.saturating_mul(2)).min(self.ceiling);
                    g.stalled = true;
                    let first = !std::mem::replace(&mut g.warned_growth, true);
                    g.in_flight = wanted;
                    break Some(GateEvent::Grew {
                        from,
                        to: g.budget,
                        in_flight: wanted - n,
                        first,
                    });
                }
                let wake = quiet_since.unwrap_or(now) + self.stall_grace;
                let _ = self.cond.wait_until(&mut g, wake);
                continue;
            }
            // Past the ceiling: wait for the reader, failing only after
            // `fail_after` with no progress.
            if quiet >= self.fail_after {
                let in_flight = g.in_flight;
                drop(g);
                return Err(self.ceiling_error(in_flight, n));
            }
            if quiet >= self.stall_grace && !std::mem::replace(&mut g.warned_ceiling, true) {
                let in_flight = g.in_flight;
                drop(g);
                self.log_event(&GateEvent::AtCeiling { in_flight });
                g = self.inner.lock();
                continue;
            }
            let base = quiet_since.unwrap_or(now);
            let wake = if quiet < self.stall_grace {
                base + self.stall_grace
            } else {
                base + self.fail_after
            };
            // Record progress does not notify: sample it periodically rather
            // than adding up to another full ceiling timeout after output stops.
            let _ = self.cond.wait_until(&mut g, wake.min(now + Duration::from_secs(1)));
        };
        if event.is_none() {
            g.in_flight = g.in_flight.saturating_add(n);
        }
        drop(g);
        if let Some(event) = event {
            self.log_event(&event);
        }
        Ok(true)
    }

    /// Log a gate event (see [`GateEvent`]): the first growth and the first
    /// wait at the ceiling as warnings, later growths at debug level. The
    /// messages state what was observed (no aligner output was consumed) and
    /// its possible causes, not a single diagnosis.
    fn log_event(&self, event: &GateEvent) {
        match *event {
            GateEvent::Grew { from, to, in_flight, first: true } => log::warn!(
                "align: no aligner output was consumed for {grace:?} while {in_flight} of \
                 unmapped reads filled the in-flight budget; raising it from {from} to {to} (at \
                 most {ceiling}, --aligner::max-read-ahead). The aligner either reads more input \
                 than the budget before its first output or is slow per chunk. To avoid this \
                 pause, size --aligner::chunk-size to the aligner's read-ahead: with \
                 --aligner::command it sets only this budget, so raise it; with a preset it also \
                 sets -K, so lower it.",
                grace = self.stall_grace,
                in_flight = ByteSize(in_flight),
                from = ByteSize(from),
                to = ByteSize(to),
                ceiling = ByteSize(self.ceiling),
            ),
            GateEvent::Grew { from, to, .. } => log::debug!(
                "align: in-flight budget raised from {} to {} (at most {})",
                ByteSize(from),
                ByteSize(to),
                ByteSize(self.ceiling),
            ),
            GateEvent::AtCeiling { in_flight } => log::warn!(
                "align: {in_flight} of unmapped reads are in flight, the read-ahead ceiling \
                 ({ceiling}, --aligner::max-read-ahead), and no aligner output has been consumed \
                 for {grace:?}; waiting for the aligner, and failing if it writes nothing for \
                 {fail_after:?}",
                in_flight = ByteSize(in_flight),
                ceiling = ByteSize(self.ceiling),
                grace = self.stall_grace,
                fail_after = self.fail_after,
            ),
        }
    }

    /// The error for a batch of `n` bytes that would take `in_flight` past the
    /// ceiling with no aligner output consumed for `fail_after`.
    fn ceiling_error(&self, in_flight: u64, n: u64) -> io::Error {
        io::Error::other(format!(
            "the aligner read ahead more than {ceiling} ({ceiling_bytes} bytes) of unmapped reads \
             without its output being consumed: {in_flight} were in flight, the next batch needs \
             {n} more, and no output was consumed for {fail_after:?}. The aligner reads more \
             input before its first output than the ceiling allows, or has stopped. Raise \
             --aligner::max-read-ahead, or size --aligner::chunk-size to the aligner's \
             read-ahead (with --aligner::command it sets only the in-flight budget, so raise it; \
             with a preset it also sets -K, so lower it).",
            ceiling = ByteSize(self.ceiling),
            ceiling_bytes = self.ceiling,
            in_flight = ByteSize(in_flight),
            n = ByteSize(n),
            fail_after = self.fail_after,
        ))
    }

    /// Record consumed aligner output without releasing any unmapped bytes.
    /// A waiting writer samples this counter on its timed wakeups, so record
    /// reads do not contend on the byte-accounting lock or wake the writer.
    pub(crate) fn record_progress(&self) {
        // The counter carries no data; only a change in its value matters.
        self.progress.fetch_add(1, Ordering::Relaxed);
    }

    /// Release `n` in-flight bytes (reader consumed a token) and wake the
    /// writer if it is blocked. Ends a stall (see the type docs).
    pub(crate) fn release(&self, n: u64) {
        let mut g = self.inner.lock();
        g.in_flight = g.in_flight.saturating_sub(n);
        self.record_progress();
        g.stalled = false;
        drop(g);
        self.cond.notify_all();
    }

    /// Mark the reader blocked handing a batch downstream for the life of the
    /// returned guard: no aligner output is consumed then, but the aligner is
    /// not stalled, so the time does not count toward a growth or a failure.
    pub(crate) fn downstream_wait(&self) -> DownstreamWait<'_> {
        self.inner.lock().downstream_blocked = true;
        self.cond.notify_all();
        DownstreamWait(self)
    }

    /// Latch that the consumer (reader) has exited, so a blocked writer
    /// wakes and bails instead of hanging forever.
    pub(crate) fn mark_consumer_gone(&self) {
        let mut g = self.inner.lock();
        g.consumer_gone = true;
        drop(g);
        self.cond.notify_all();
    }

    /// The current budget, which stalls may have raised.
    #[cfg(test)]
    pub(crate) fn budget(&self) -> u64 {
        self.inner.lock().budget
    }

    /// The bytes currently reserved.
    #[cfg(test)]
    pub(crate) fn in_flight(&self) -> u64 {
        self.inner.lock().in_flight
    }
}

/// The reader's downstream wait (see [`InFlightGate::downstream_wait`]). Its
/// end counts as reader progress, restarting the writer's quiet clock.
pub(crate) struct DownstreamWait<'a>(&'a InFlightGate);

impl Drop for DownstreamWait<'_> {
    fn drop(&mut self) {
        let mut g = self.0.inner.lock();
        g.downstream_blocked = false;
        self.0.record_progress();
        drop(g);
        self.0.cond.notify_all();
    }
}

/// RAII guard that latches consumer-gone on drop — on normal return, on an
/// error return, AND on a panic unwinding out of `reader_loop`. Without the
/// guard the panic path would skip `mark_consumer_gone`, leaving a writer
/// parked in `gate.acquire` blocked forever and `Drop`'s `writer_thread.join()`
/// hung with it.
pub(crate) struct ConsumerGoneGuard<'a>(pub(crate) &'a InFlightGate);

impl Drop for ConsumerGoneGuard<'_> {
    fn drop(&mut self) {
        self.0.mark_consumer_gone();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    #[rstest::rstest]
    #[case::default_chunk(150_000_000, 1_800_000_000)]
    #[case::zero_chunk(0, 536_870_912)]
    #[case::small_chunk(1_000_000, 536_870_912)]
    fn in_flight_budget_derivation(#[case] chunk_size: u64, #[case] expected: u64) {
        assert_eq!(in_flight_budget_for_chunk_size(chunk_size), expected);
    }

    /// A gate whose stall grace is long enough that no test waits it out: it
    /// behaves as a plain blocking budget unless released.
    fn blocking_gate(budget: u64) -> InFlightGate {
        InFlightGate::new(budget, u64::MAX, Duration::from_hours(1), Duration::from_hours(1))
    }

    #[test]
    fn consumer_gone_latched_even_when_reader_scope_panics() {
        // A panic unwinding out of the reader scope must still latch
        // consumer-gone via `ConsumerGoneGuard`, so a writer parked in
        // `acquire` bails (returns false) instead of blocking forever. Without
        // the guard the panic would skip `mark_consumer_gone`, deadlocking the
        // writer and the `writer_thread.join()` in `Drop`.
        let gate = blocking_gate(1024);
        assert!(gate.acquire(10).unwrap(), "acquire succeeds before the consumer is gone");

        let panicked = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _consumer_gone = ConsumerGoneGuard(&gate);
            panic!("simulated reader_loop panic");
        }));
        assert!(panicked.is_err(), "the panic must propagate out of the guarded scope");

        // The guard's `Drop` latched consumer-gone during unwind, so a
        // subsequent writer reservation bails instead of blocking.
        assert!(
            !gate.acquire(10).unwrap(),
            "a writer must bail once the reader scope has panicked"
        );
    }

    /// A batch larger than the whole budget, and even the ceiling, must still
    /// pass when the gate is empty — it can't be split, so holding it is
    /// unavoidable, and blocking would deadlock.
    #[test]
    fn in_flight_gate_oversized_single_batch_passes_when_empty() {
        let gate = InFlightGate::new(100, 500, Duration::ZERO, Duration::ZERO);
        assert!(gate.acquire(1000).unwrap(), "oversized batch must pass when gate empty");
        assert_eq!(gate.budget(), 100, "admitting into an empty gate does not grow it");
    }

    #[test]
    fn in_flight_gate_blocks_until_release() {
        let gate = Arc::new(blocking_gate(100));
        assert!(gate.acquire(100).unwrap(), "first acquire fills the budget");

        let g2 = Arc::clone(&gate);
        let (tx, rx) = std::sync::mpsc::channel();
        let h = std::thread::spawn(move || {
            // in_flight is 100 (full, non-empty) → this blocks until release.
            assert!(g2.acquire(50).unwrap());
            tx.send(()).unwrap();
        });
        // Still blocked.
        assert!(
            rx.recv_timeout(Duration::from_millis(150)).is_err(),
            "second acquire must block while the gate is full"
        );
        gate.release(100);
        rx.recv_timeout(Duration::from_secs(5)).expect("acquire must unblock after release");
        h.join().unwrap();
        assert_eq!(gate.budget(), 100, "a release, not a stall, admitted the batch");
    }

    #[test]
    fn in_flight_gate_consumer_gone_unblocks_and_bails() {
        let gate = Arc::new(blocking_gate(100));
        assert!(gate.acquire(100).unwrap(), "fill the budget");

        let g2 = Arc::clone(&gate);
        let (tx, rx) = std::sync::mpsc::channel();
        let h = std::thread::spawn(move || {
            let ok = g2.acquire(50).unwrap(); // blocks (full) until consumer-gone
            tx.send(ok).unwrap();
        });
        assert!(rx.recv_timeout(Duration::from_millis(150)).is_err(), "must be blocked");
        gate.mark_consumer_gone();
        let ok = rx.recv_timeout(Duration::from_secs(5)).expect("must wake on consumer-gone");
        assert!(!ok, "acquire must return false once the consumer is gone");
        h.join().unwrap();
    }

    /// With no release during the grace period, a full gate admits the batch
    /// and grows: to twice the budget, or to what the batch needs if that is
    /// more, but never past the ceiling.
    #[rstest::rstest]
    #[case::doubles(100, 50, 1000, 200)]
    #[case::to_what_the_batch_needs(100, 250, 1000, 350)]
    #[case::capped_at_the_ceiling(100, 50, 160, 160)]
    #[case::to_exactly_the_ceiling(100, 60, 160, 160)]
    fn in_flight_gate_grows_on_a_stall(
        #[case] budget: u64,
        #[case] batch: u64,
        #[case] ceiling: u64,
        #[case] expected_budget: u64,
    ) {
        let gate = InFlightGate::new(budget, ceiling, Duration::ZERO, Duration::ZERO);
        assert!(gate.acquire(budget).unwrap(), "fill the budget");
        assert!(gate.acquire(batch).unwrap(), "a stalled gate admits the batch");
        assert_eq!(gate.budget(), expected_budget);
    }

    /// The writer waits the whole grace period on a full gate before it grows
    /// it, rather than growing on the first full acquire.
    #[test]
    fn in_flight_gate_waits_the_grace_before_growing() {
        let grace = Duration::from_millis(200);
        let gate = InFlightGate::new(100, 1000, grace, grace);
        assert!(gate.acquire(100).unwrap(), "fill the budget");
        let start = Instant::now();
        assert!(gate.acquire(50).unwrap(), "admitted after the grace period");
        assert!(start.elapsed() >= grace, "grew after {:?}, before the grace", start.elapsed());
        assert_eq!(gate.budget(), 200);
    }

    /// A release while the writer waits restarts the grace period: an aligner
    /// that is producing output is streaming, not stalled, so the gate does not
    /// grow until it has been quiet for a whole grace period.
    #[test]
    fn in_flight_gate_release_restarts_the_grace() {
        let grace = Duration::from_millis(300);
        let gate = Arc::new(InFlightGate::new(100, 1000, grace, grace));
        assert!(gate.acquire(100).unwrap(), "fill the budget");

        let g2 = Arc::clone(&gate);
        let start = Instant::now();
        let h = std::thread::spawn(move || {
            // Needs 100 + 50 - 10 = 140 > 100 even after the release below.
            assert!(g2.acquire(50).unwrap());
            start.elapsed()
        });
        std::thread::sleep(grace / 2);
        let released_at = start.elapsed();
        gate.release(10);
        let waited = h.join().unwrap();
        assert!(
            waited >= released_at + grace,
            "grew {waited:?} after the start, within a grace of the release at {released_at:?}"
        );
        assert_eq!(gate.budget(), 200);
    }

    /// Once a stall has grown the gate, further batches are admitted at once
    /// (still growing) until a release ends the stall; after that a full gate
    /// waits a whole grace period again.
    #[test]
    fn in_flight_gate_stall_lasts_until_a_release() {
        let grace = Duration::from_millis(200);
        let gate = InFlightGate::new(100, 10_000, grace, grace);
        assert!(gate.acquire(100).unwrap(), "fill the budget");
        assert!(gate.acquire(100).unwrap(), "first stall: grows to 200 after the grace");
        assert_eq!(gate.budget(), 200);

        // Still stalled (no release): 200 + 150 > 200 grows without waiting a
        // grace period. The bound is loose (half a grace), but a gate that
        // waited a whole grace per batch would fail it.
        let start = Instant::now();
        assert!(gate.acquire(150).unwrap());
        assert!(start.elapsed() < grace / 2, "a stalled gate admits at once");
        assert_eq!(gate.budget(), 400);

        // A release ends the stall; the grown budget stays.
        gate.release(350);
        assert_eq!(gate.budget(), 400);
        assert!(gate.acquire(350).unwrap(), "fits the grown budget: 0 + 350 <= 400");
        let start = Instant::now();
        assert!(gate.acquire(100).unwrap(), "full again: grows after a fresh grace");
        assert!(start.elapsed() >= grace, "a new stall must wait the grace again");
        assert_eq!(gate.budget(), 800);
    }

    /// Record progress between acquisitions ends a grown stall without freeing
    /// bytes, so another growth must wait a fresh grace period.
    #[test]
    fn in_flight_gate_record_progress_ends_a_grown_stall() {
        let grace = Duration::from_millis(200);
        let gate = InFlightGate::new(100, 1000, grace, grace);
        assert!(gate.acquire(100).unwrap());
        assert!(gate.acquire(100).unwrap());
        assert_eq!(gate.budget(), 200);

        gate.record_progress();
        assert_eq!(gate.in_flight(), 200, "record progress does not release bytes");
        let start = Instant::now();
        assert!(gate.acquire(50).unwrap());
        assert!(start.elapsed() >= grace, "record progress must restart the growth grace");
        assert_eq!(gate.budget(), 400);
        assert_eq!(gate.in_flight(), 250);
    }

    /// A batch that would take in-flight past the ceiling waits for a release
    /// and fails only once none came for the (longer) ceiling timeout, not
    /// after the growth grace; the error names the ceiling and both options
    /// that move it, and the batch is not admitted.
    #[rstest::rstest]
    #[case::first_full_acquire(100, None, 150)]
    #[case::after_growing_to_the_ceiling(100, Some(50), 160)]
    fn in_flight_gate_fails_past_the_ceiling(
        #[case] budget: u64,
        #[case] growing_batch: Option<u64>,
        #[case] ceiling: u64,
    ) {
        let grace = Duration::from_millis(50);
        let fail_after = Duration::from_millis(400);
        let gate = InFlightGate::new(budget, ceiling, grace, fail_after);
        assert!(gate.acquire(budget).unwrap(), "fill the budget");
        if let Some(batch) = growing_batch {
            assert!(gate.acquire(batch).unwrap(), "a stall grows the gate to the ceiling");
            assert_eq!(gate.budget(), ceiling);
        }
        let reserved = gate.in_flight();
        let start = Instant::now();
        let err = gate.acquire(1000).expect_err("past the ceiling must fail, not block");
        assert!(start.elapsed() >= fail_after, "fails only after the ceiling timeout");
        assert_eq!(gate.in_flight(), reserved, "the failing batch is not admitted");
        let msg = err.to_string();
        for needle in [
            "read ahead more than",
            &format!("({ceiling} bytes)"),
            "--aligner::max-read-ahead",
            "--aligner::chunk-size",
        ] {
            assert!(msg.contains(needle), "error must contain {needle:?}: {msg}");
        }
    }

    /// At the ceiling, a release before the timeout admits the batch instead of
    /// failing: a slow aligner that does produce output is not failed.
    #[test]
    fn in_flight_gate_at_the_ceiling_waits_for_a_slow_aligner() {
        let grace = Duration::from_millis(50);
        let gate = Arc::new(InFlightGate::new(100, 100, grace, Duration::from_secs(30)));
        assert!(gate.acquire(100).unwrap(), "fill the budget, which is the ceiling");
        let releaser = Arc::clone(&gate);
        let h = std::thread::spawn(move || {
            std::thread::sleep(grace * 4);
            releaser.release(100);
        });
        assert!(gate.acquire(50).unwrap(), "admitted once the aligner's output is consumed");
        h.join().unwrap();
        assert_eq!(gate.in_flight(), 50);
    }

    /// Time the reader spends blocked handing a batch downstream does not count
    /// as an aligner stall: a full gate grows only after a whole grace with no
    /// progress once that wait has ended, and at the ceiling it does not fail.
    #[rstest::rstest]
    #[case::grows_only_after_the_wait(1000)]
    #[case::does_not_fail_at_the_ceiling(100)]
    fn in_flight_gate_ignores_downstream_backpressure(#[case] ceiling: u64) {
        let grace = Duration::from_millis(100);
        let blocked_for = Duration::from_millis(500);
        let gate = Arc::new(InFlightGate::new(100, ceiling, grace, grace * 2));
        assert!(gate.acquire(100).unwrap(), "fill the budget");
        let reader = Arc::clone(&gate);
        let (blocked_tx, blocked_rx) = std::sync::mpsc::channel();
        let h = std::thread::spawn(move || {
            let downstream = reader.downstream_wait();
            blocked_tx.send(()).unwrap();
            std::thread::sleep(blocked_for);
            drop(downstream);
            if ceiling == 100 {
                // At the ceiling the batch can only be admitted by a release.
                std::thread::sleep(grace);
                reader.release(100);
            }
        });
        blocked_rx.recv().unwrap();
        let start = Instant::now();
        assert!(gate.acquire(50).unwrap(), "admitted, not failed");
        assert!(start.elapsed() >= blocked_for, "no growth while the reader is blocked downstream");
        h.join().unwrap();
    }

    /// The configured budget never exceeds the ceiling.
    #[test]
    fn in_flight_gate_budget_is_capped_at_the_ceiling() {
        assert_eq!(InFlightGate::new(1000, 100, Duration::ZERO, Duration::ZERO).budget(), 100);
    }

    /// `--aligner::max-read-ahead` wins when set; otherwise the ceiling is a
    /// quarter of host memory, but never below the budget.
    #[rstest::rstest]
    #[case::requested(100, Some(50), 64 << 30, 50)]
    #[case::quarter_of_host_memory(1 << 30, None, 64 << 30, 16 << 30)]
    #[case::at_least_the_budget(8 << 30, None, 16 << 30, 8 << 30)]
    #[case::unknown_host_memory(1 << 30, None, 0, 1 << 30)]
    fn read_ahead_ceiling_resolution(
        #[case] budget: u64,
        #[case] requested: Option<u64>,
        #[case] host_memory: u64,
        #[case] expected: u64,
    ) {
        assert_eq!(read_ahead_ceiling(budget, requested, host_memory), expected);
    }

    /// The hidden budget override takes a positive byte count; anything else
    /// is ignored.
    #[rstest::rstest]
    #[case::unset(None, None)]
    #[case::bytes(Some("65536"), Some(65_536))]
    #[case::padded(Some(" 42 "), Some(42))]
    #[case::zero(Some("0"), None)]
    #[case::not_a_number(Some("64KiB"), None)]
    fn in_flight_budget_override_parsing(
        #[case] value: Option<&str>,
        #[case] expected: Option<u64>,
    ) {
        assert_eq!(in_flight_budget_override(value), expected);
    }

    /// The hidden ceiling-timeout override takes a positive number of seconds;
    /// anything else is ignored.
    #[rstest::rstest]
    #[case::unset(None, None)]
    #[case::seconds(Some("2"), Some(Duration::from_secs(2)))]
    #[case::zero(Some("0"), None)]
    #[case::not_a_number(Some("2s"), None)]
    fn ceiling_timeout_override_parsing(
        #[case] value: Option<&str>,
        #[case] expected: Option<Duration>,
    ) {
        assert_eq!(ceiling_timeout_override(value), expected);
    }
}

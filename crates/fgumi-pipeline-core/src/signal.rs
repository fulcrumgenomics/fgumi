//! Pipeline-wide error broadcast and cancellation, folded into a single
//! shared atomic so the worker loop pays one relaxed load per iteration.

use std::io;
use std::sync::Arc;
use std::sync::OnceLock;
use std::sync::atomic::{AtomicU8, Ordering as AtomicOrdering};

const STATE_OK: u8 = 0;
const STATE_CANCELLED: u8 = 1;
const STATE_ERROR: u8 = 2;

/// Errors a pipeline run can return to its caller.
#[derive(Debug)]
pub enum PipelineError {
    /// A step's `try_run` returned `Err`. Carries the originating step's
    /// name and the underlying I/O error.
    Io { step: &'static str, source: io::Error },
    /// `cancel_handle().cancel()` was called from outside the pipeline.
    Cancelled,
    /// The chain has more `Exclusive` steps than the configured thread count.
    NotEnoughThreads { required: usize, available: usize },
    /// The deadlock monitor observed no global progress for `stalled_secs`
    /// while work was still stuck in queues/reorder buffers — a wedge. The
    /// pipeline is failed fast rather than left to hang forever.
    TimedOut { stalled_secs: u64 },
    /// The pipeline was built with the deadlock monitor armed
    /// (`deadlock_timeout_secs > 0`), but `step` declares a non-`ByteBounded`
    /// output transport (`spec`, e.g. `CountBounded`/`Unbounded`). The
    /// `in_flight_bytes` probe cannot see a wedge on such an edge, so fail-fast
    /// would be silently disabled there. The pipeline is rejected at startup
    /// rather than allowed to hang. This is a chain-construction error, not a
    /// runtime condition: production chains wire only `ByteBounded` transports.
    MonitorBlindTransport { step: &'static str, spec: String },
}

impl PipelineError {
    /// Reconstruct an owned copy of this error.
    ///
    /// `PipelineError` is not `Clone` because its `Io` variant holds an
    /// [`io::Error`], which is not `Clone`; this rebuilds that source from its
    /// `kind()` + display string. Used to turn the borrowed `signal.outcome()`
    /// into the owned error returned from the run drivers (`Pipeline::run` and
    /// the fused single-thread driver), which previously open-coded the same
    /// per-variant rematch in two places.
    pub(crate) fn reconstruct(&self) -> PipelineError {
        match self {
            Self::Cancelled => Self::Cancelled,
            Self::Io { step, source } => {
                Self::Io { step, source: io::Error::new(source.kind(), format!("{source}")) }
            }
            Self::NotEnoughThreads { required, available } => {
                Self::NotEnoughThreads { required: *required, available: *available }
            }
            Self::TimedOut { stalled_secs } => Self::TimedOut { stalled_secs: *stalled_secs },
            Self::MonitorBlindTransport { step, spec } => {
                Self::MonitorBlindTransport { step, spec: spec.clone() }
            }
        }
    }
}

impl std::fmt::Display for PipelineError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io { step, source } => write!(f, "step {step:?} failed: {source}"),
            Self::Cancelled => write!(f, "pipeline cancelled"),
            Self::NotEnoughThreads { required, available } => write!(
                f,
                "pipeline requires {required} threads for Exclusive steps, only {available} available"
            ),
            Self::TimedOut { stalled_secs } => write!(
                f,
                "pipeline deadlock detected: no progress for {stalled_secs}s with work still in flight"
            ),
            Self::MonitorBlindTransport { step, spec } => write!(
                f,
                "deadlock monitor is armed but step {step:?} declares a {spec} output transport, \
                 which is invisible to the in_flight_bytes probe — a wedge on that edge would \
                 silently disable fail-fast. Production transports must be QueueSpec::ByteBounded."
            ),
        }
    }
}

impl std::error::Error for PipelineError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io { source, .. } => Some(source),
            _ => None,
        }
    }
}

/// Single shared per-pipeline state covering both error broadcast and cancel.
/// Workers check `is_done()` once per iteration — one relaxed atomic load.
#[derive(Default)]
pub struct PipelineSignal {
    state: AtomicU8,
    payload: OnceLock<PipelineError>,
}

impl PipelineSignal {
    #[must_use]
    pub fn new() -> Arc<Self> {
        Arc::new(Self::default())
    }

    /// Hot-path: workers call this once per loop iteration.
    ///
    /// Uses `Relaxed` because the only thing a worker does on `true` is stop
    /// polling — workers never read the payload. The payload is read only by
    /// the run drivers via [`Self::to_result`], after every internal writer
    /// (workers, deadlock monitor) has been joined, so the join supplies the
    /// happens-before that makes a recorded error's payload visible.
    ///
    /// The one writer not covered by a join is an external
    /// [`CancelHandle::cancel`], which races the driver's terminal read with no
    /// synchronizing edge; [`Self::to_result`] derives `Cancelled` from the
    /// (coherence-visible) state rather than the payload to close that window.
    #[inline]
    #[must_use]
    pub fn is_done(&self) -> bool {
        self.state.load(AtomicOrdering::Relaxed) != STATE_OK
    }

    #[inline]
    #[must_use]
    pub fn is_cancelled(&self) -> bool {
        self.state.load(AtomicOrdering::Relaxed) == STATE_CANCELLED
    }

    /// First writer wins; later writers are silently dropped (the
    /// `compare_exchange` rejects the state transition and the
    /// `OnceLock::set` rejects the payload write).
    pub fn record_error(&self, err: PipelineError) {
        // Set the payload only if we won the state CAS, so the error surfaced by
        // `outcome()` always matches the writer that transitioned the state.
        // Otherwise writer A could win the state CAS while writer B wins the
        // `OnceLock::set`, leaving `state == ERROR` reporting B's error —
        // violating the documented "first writer wins" contract.
        if self
            .state
            .compare_exchange(
                STATE_OK,
                STATE_ERROR,
                AtomicOrdering::Release,
                AtomicOrdering::Relaxed,
            )
            .is_ok()
        {
            let _ = self.payload.set(err);
        }
    }

    /// First writer wins; later writers (including a `record_error` after
    /// `cancel`) are silently dropped.
    pub fn cancel(&self) {
        // Only publish the payload if THIS call won the state transition.
        // Setting it unconditionally races a concurrent `record_error`: if
        // that call wins the CAS (state → ERROR) but has not yet set its
        // payload, an unconditional `payload.set(Cancelled)` here can win the
        // `OnceLock` and leave `state == ERROR` while `outcome()` reports
        // `Cancelled` — an inconsistent state/payload pair. Guarding on the CAS
        // (as `record_error` does) keeps the two consistent: whoever wins the
        // state transition is the one that sets the payload.
        // (See the loom NOTE at the end of the tests module for why this
        // race-guard is not exercised by a unit test.)
        if self
            .state
            .compare_exchange(
                STATE_OK,
                STATE_CANCELLED,
                AtomicOrdering::Release,
                AtomicOrdering::Relaxed,
            )
            .is_ok()
        {
            let _ = self.payload.set(PipelineError::Cancelled);
        }
    }

    /// Read the recorded error payload, if any.
    ///
    /// The `OnceLock::get` acquire fence pairs with the `OnceLock::set` release
    /// on the writer side, so a reader synchronized with the writer (via the
    /// worker/monitor join in the run drivers) sees a recorded error's payload.
    ///
    /// Beware the gap this does NOT cover: `is_done() == true` does not imply
    /// `outcome().is_some()`. `cancel()`/`record_error()` publish the terminal
    /// `state` (CAS) before `payload.set()`, so a reader can observe the
    /// terminal state while the payload `OnceLock` is still empty — most
    /// reachably for an external `CancelHandle::cancel`, whose writer is never
    /// joined against the driver's read. Map an outcome to a run result through
    /// [`Self::to_result`], which handles that window, not by branching on
    /// `outcome()` directly.
    #[must_use]
    pub fn outcome(&self) -> Option<&PipelineError> {
        self.payload.get()
    }

    /// Map the recorded outcome to a run `Result`, as returned by the run
    /// drivers ([`crate::builder::Pipeline::run`] and the fused single-thread
    /// driver).
    ///
    /// A recorded error's payload is always visible here: every `record_error`
    /// writer (workers, the deadlock monitor, the fused driver itself) is
    /// joined before this read, so the join orders its `payload.set()` ahead of
    /// the read.
    ///
    /// An external [`CancelHandle::cancel`] is the exception — its writer is
    /// never joined against this read, so there is no happens-before edge to
    /// make its `OnceLock` payload visible, and it can leave `state ==
    /// CANCELLED` with `outcome() == None`. A naive `match outcome()` would then
    /// map a genuinely cancelled run (workers observed `is_done()` and stopped
    /// early) to `Ok(())`. `Cancelled` carries no payload data, so synthesize it
    /// from the coherence-visible `state` via [`Self::is_cancelled`] instead.
    pub(crate) fn to_result(&self) -> Result<(), PipelineError> {
        match self.outcome() {
            Some(err) => Err(err.reconstruct()),
            None if self.is_cancelled() => Err(PipelineError::Cancelled),
            None => Ok(()),
        }
    }
}

/// External handle for cancelling an in-flight pipeline.
#[derive(Clone)]
pub struct CancelHandle {
    signal: Arc<PipelineSignal>,
}

impl CancelHandle {
    /// Construct a `CancelHandle` from a shared signal. Used by
    /// `Pipeline::cancel_handle` (see `builder.rs`) to hand a handle to the
    /// caller of a built pipeline, before `Pipeline::run` kicks off worker
    /// threads.
    pub(crate) fn from_signal(signal: Arc<PipelineSignal>) -> Self {
        Self { signal }
    }

    pub fn cancel(&self) {
        self.signal.cancel();
    }

    #[must_use]
    pub fn is_cancelled(&self) -> bool {
        self.signal.is_cancelled()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fresh_signal_is_not_done() {
        let signal = PipelineSignal::new();
        assert!(!signal.is_done());
        assert!(!signal.is_cancelled());
    }

    #[test]
    fn record_error_marks_done() {
        let signal = PipelineSignal::new();
        signal.record_error(PipelineError::Io {
            step: "test_step",
            source: io::Error::other("boom"),
        });
        assert!(signal.is_done());
        assert!(!signal.is_cancelled());
        assert!(matches!(signal.outcome(), Some(PipelineError::Io { step: "test_step", .. })));
    }

    #[test]
    fn cancel_marks_done_and_cancelled() {
        let signal = PipelineSignal::new();
        signal.cancel();
        assert!(signal.is_done());
        assert!(signal.is_cancelled());
        assert!(matches!(signal.outcome(), Some(PipelineError::Cancelled)));
    }

    #[test]
    fn first_recorded_error_wins() {
        let signal = PipelineSignal::new();
        signal.record_error(PipelineError::Io { step: "first", source: io::Error::other("first") });
        signal
            .record_error(PipelineError::Io { step: "second", source: io::Error::other("second") });
        match signal.outcome().unwrap() {
            PipelineError::Io { step, .. } => assert_eq!(*step, "first"),
            other => panic!("expected Io, got {other:?}"),
        }
    }

    #[test]
    fn cancel_handle_propagates() {
        let signal = PipelineSignal::new();
        let handle = CancelHandle::from_signal(Arc::clone(&signal));
        assert!(!handle.is_cancelled());
        handle.cancel();
        assert!(handle.is_cancelled());
        assert!(signal.is_cancelled());
    }

    #[test]
    fn reconstruct_preserves_every_variant() {
        // `reconstruct` rebuilds an owned error from a borrow (PipelineError
        // isn't Clone). Each variant must round-trip; the `Io` source is
        // rebuilt from kind + display, the rest are copied verbatim.
        let io = PipelineError::Io { step: "s", source: io::Error::other("boom") }.reconstruct();
        match io {
            PipelineError::Io { step, source } => {
                assert_eq!(step, "s");
                assert_eq!(source.kind(), io::ErrorKind::Other);
                assert_eq!(source.to_string(), "boom");
            }
            other => panic!("expected Io, got {other:?}"),
        }

        assert!(matches!(PipelineError::Cancelled.reconstruct(), PipelineError::Cancelled));

        let net = PipelineError::NotEnoughThreads { required: 4, available: 2 }.reconstruct();
        assert!(matches!(net, PipelineError::NotEnoughThreads { required: 4, available: 2 }));

        let to = PipelineError::TimedOut { stalled_secs: 60 }.reconstruct();
        assert!(matches!(to, PipelineError::TimedOut { stalled_secs: 60 }));

        let blind = PipelineError::MonitorBlindTransport {
            step: "s",
            spec: "QueueSpec::Unbounded".to_string(),
        }
        .reconstruct();
        assert!(
            matches!(blind, PipelineError::MonitorBlindTransport { step: "s", spec } if spec == "QueueSpec::Unbounded")
        );
    }

    #[test]
    fn to_result_maps_clean_error_and_cancel() {
        let clean = PipelineSignal::new();
        assert!(clean.to_result().is_ok());

        let errored = PipelineSignal::new();
        errored.record_error(PipelineError::Io { step: "s", source: io::Error::other("boom") });
        assert!(matches!(errored.to_result(), Err(PipelineError::Io { step: "s", .. })));

        let cancelled = PipelineSignal::new();
        cancelled.cancel();
        assert!(matches!(cancelled.to_result(), Err(PipelineError::Cancelled)));
    }

    #[test]
    fn to_result_reports_cancel_when_payload_not_yet_published() {
        // Reproduce the external-cancel window: a `CancelHandle::cancel` on an
        // un-joined thread has published the terminal `state` (the CAS) but has
        // not yet run `payload.set(Cancelled)`. A driver reading the outcome in
        // that window sees `is_done() == true` but `outcome() == None`. Branching
        // on `outcome()` alone would map this genuinely cancelled run to `Ok(())`;
        // `to_result` must instead synthesize `Cancelled` from the state.
        let signal = PipelineSignal::new();
        signal.state.store(STATE_CANCELLED, AtomicOrdering::Release);
        assert!(signal.is_done(), "terminal state must read as done");
        assert!(signal.is_cancelled());
        assert!(signal.outcome().is_none(), "payload not set in this window");
        assert!(
            matches!(signal.to_result(), Err(PipelineError::Cancelled)),
            "a cancel observed before its payload is published must not map to Ok"
        );
    }

    // NOTE: the `cancel` vs `record_error` state/payload-consistency fix (guarding
    // `payload.set` behind the CAS) is deliberately NOT covered by a unit test.
    // The bug is a true concurrent interleaving — `cancel` losing the CAS in the
    // window between `record_error`'s CAS-win and its own `payload.set` — and
    // `OnceLock` masks the inconsistency in every *sequential* ordering, so no
    // single-threaded test can distinguish fixed from buggy. A spawned-thread
    // stress loop does not reliably hit the two-instruction window either (a
    // 2000-trial loop passed against the buggy code). Deterministically forcing
    // the interleaving would require `loom`; until that dependency is justified,
    // the fix rests on the inline argument at the `cancel` call site.
}

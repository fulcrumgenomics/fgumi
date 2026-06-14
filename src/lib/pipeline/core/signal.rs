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
    /// Uses `Relaxed` because the only thing the worker does on `true` is
    /// stop polling — it doesn't read the payload synchronously. Code that
    /// needs to read the payload (`outcome()`) goes through `OnceLock`,
    /// which provides its own Acquire/Release synchronization on the cell.
    /// So `Relaxed` here is sufficient: any worker that observes `true`
    /// and *then* calls `outcome()` will see the payload via the
    /// `OnceLock` acquire fence on `OnceLock::get`.
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
        let _ = self.state.compare_exchange(
            STATE_OK,
            STATE_CANCELLED,
            AtomicOrdering::Release,
            AtomicOrdering::Relaxed,
        );
        let _ = self.payload.set(PipelineError::Cancelled);
    }

    /// Read the recorded outcome. The `OnceLock::get` acquire fence pairs
    /// with the `OnceLock::set` release on the writer side, so any worker
    /// that observes `is_done() == true` and then calls `outcome()` sees
    /// the payload.
    #[must_use]
    pub fn outcome(&self) -> Option<&PipelineError> {
        self.payload.get()
    }
}

/// External handle for cancelling an in-flight pipeline.
#[derive(Clone)]
pub struct CancelHandle {
    signal: Arc<PipelineSignal>,
}

impl CancelHandle {
    /// Construct a `CancelHandle` from a shared signal. Used internally by
    /// `Pipeline::run` (Phase 2) to hand a handle to the caller before
    /// kicking off worker threads.
    #[allow(dead_code)] // wired up by Phase 2 Pipeline::run.
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
}

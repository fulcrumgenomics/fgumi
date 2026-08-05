//! One-shot lazy SAM/BAM header handle.
//!
//! `HeaderHandle` lets a downstream sink open its output file before the
//! full output header is known, and consume the header at first record
//! time once an upstream step has produced it. The primary motivating
//! consumer is the writer downstream of `AlignAndMergeStep`: the
//! aligner's `@PG` (and any `@RG`/`@CO` lines it adds) are runtime
//! contributions that aren't available at `Pipeline::build` time.
//!
//! ### Contract
//!
//! - `set` and `poison` are one-shot — calling either after the handle
//!   has been resolved returns `AlreadySetError` without altering state.
//! - `try_get` is non-blocking and dispatcher-friendly: a step that
//!   needs the header but observes `None` returns `StepOutcome::NoProgress`
//!   and is rescheduled.
//! - Cloning shares state (via `Arc<OnceLock<...>>`), so multiple
//!   readers see the same set-or-poison outcome. The cross-thread
//!   happens-before guarantee comes from `OnceLock`'s internal
//!   synchronization (`set` is Release, `get` is Acquire); the
//!   `Arc<OnceLock<...>>` wrapper exists only for shareability.
//!
//! ### Producer-side invariant
//!
//! Every consumer that polls `try_get` assumes the handle will
//! eventually resolve via `set` or `poison`. Producer steps owning a
//! `HeaderHandle` must poison their handle in their `Drop` impl as a
//! backstop against panics (or any exit path that bypasses normal
//! completion). The framework does not enforce this; it's a step-level
//! convention. A future revision may split the handle into typed
//! setter / reader halves so an orphaned setter is statically detectable.
//!
//! An orphaned setter — a producer that exits without resolving — does
//! **not** hang the run, provided the consumer follows the shape below.
//! Two framework backstops catch it, both pinned by tests in this module
//! (`orphaned_setter_*`):
//!
//! 1. **Producer left no items.** When a step reports `Finished` the
//!    framework closes its output edges, so the consumer's
//!    `is_input_drained()` arm fires and it can fail with a diagnostic.
//! 2. **Producer left items queued.** A header-blocked consumer pops
//!    nothing, so `is_drained()` (which is `drained && empty`) stays
//!    false and backstop 1 is unreachable — but those stranded items are
//!    in flight, so the deadlock monitor sees non-zero `in_flight_bytes`,
//!    classifies the run `Wedged` past the fatal timeout, and fails it
//!    with `PipelineError::TimedOut`.
//!
//! Backstop 2 depends on the stranded items carrying real heap bytes:
//! queue accounting is `HeapSize::heap_size()`-only, so a zero-heap item
//! type leaves `in_flight_bytes` at 0, which classifies as `Starving` and
//! resets the stall clock forever. That gap is general to the monitor
//! rather than specific to headers; see the `in_flight_bytes` docs in
//! `builder.rs`.
//!
//! ### Consumer-side contract
//!
//! A consumer must never poll `try_get` in isolation. Pair the probe with
//! `is_input_drained()` and treat "input drained before the header
//! resolved" as an error — that arm is backstop 1, and without it a
//! consumer whose producer left no items polls `None` forever.

use std::fmt;
use std::io;
use std::sync::{Arc, OnceLock};

use noodles::sam::Header;

/// Shared, one-shot header slot.
///
/// See module docs for the full contract.
#[derive(Clone, Default, Debug)]
pub struct HeaderHandle {
    inner: Arc<OnceLock<Result<Header, Arc<io::Error>>>>,
}

impl HeaderHandle {
    /// Construct an empty handle. `try_get` returns `None` until `set`
    /// or `poison` is called.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Construct a handle that already carries `header`. `try_get`
    /// returns `Some(Ok(&header))` immediately. Convenience for call
    /// sites that have an eager header today and want to opt into the
    /// handle-typed sink API without changing behavior.
    ///
    /// # Panics
    /// The internal `set` call is on a freshly-constructed cell and
    /// is therefore unreachable as a failure — the `expect` is a
    /// belt-and-braces assertion documenting the invariant.
    #[must_use]
    #[allow(clippy::missing_panics_doc)] // documented above; clippy
    // doesn't see the # Panics section across this constructor's
    // delegation to `set`.
    pub fn from_header(header: Header) -> Self {
        let handle = Self::new();
        handle.set(header).expect("fresh HeaderHandle accepts its first set");
        handle
    }

    /// Resolve the handle to `header`. Returns `AlreadySetError` if the
    /// handle was previously set or poisoned.
    ///
    /// # Errors
    /// Returns `AlreadySetError` on the second and subsequent calls.
    pub fn set(&self, header: Header) -> Result<(), AlreadySetError> {
        self.inner.set(Ok(header)).map_err(|_| AlreadySetError)
    }

    /// Resolve the handle to a failure. Subsequent `try_get` calls
    /// surface the error.
    ///
    /// # Errors
    /// Returns `AlreadySetError` on the second and subsequent calls.
    pub fn poison(&self, error: io::Error) -> Result<(), AlreadySetError> {
        self.inner.set(Err(Arc::new(error))).map_err(|_| AlreadySetError)
    }

    /// Non-blocking handle probe.
    ///
    /// Returns:
    /// - `None` if neither `set` nor `poison` has been called yet — the
    ///   caller should yield (e.g. return `StepOutcome::NoProgress`).
    /// - `Some(Ok(&header))` if `set(header)` was called.
    /// - `Some(Err(e))` if `poison(e)` was called. A fresh `io::Error`
    ///   is constructed on every call from the stored kind + message
    ///   (because `io::Error` is not `Clone`). **Note:** the
    ///   original error's `source()` chain and any structured
    ///   payload are not preserved — only kind + display string
    ///   round-trip. Producers that need to surface structured
    ///   diagnostic context should encode it into the display
    ///   string before poisoning the handle.
    #[must_use]
    pub fn try_get(&self) -> Option<io::Result<&Header>> {
        self.inner.get().map(|stored| match stored {
            Ok(header) => Ok(header),
            Err(err) => Err(io::Error::new(err.kind(), err.to_string())),
        })
    }

    /// `true` once `set` or `poison` has resolved this handle.
    #[must_use]
    pub fn is_set(&self) -> bool {
        self.inner.get().is_some()
    }
}

/// Returned by `set` / `poison` when the handle was previously
/// resolved.
#[derive(Debug, Clone, Copy)]
pub struct AlreadySetError;

impl fmt::Display for AlreadySetError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "HeaderHandle has already been set or poisoned")
    }
}

impl std::error::Error for AlreadySetError {}

#[cfg(test)]
mod tests {
    use super::*;
    use std::thread;

    #[test]
    fn new_handle_has_no_value() {
        let h = HeaderHandle::new();
        assert!(!h.is_set());
        assert!(h.try_get().is_none());
    }

    #[test]
    fn from_header_resolves_immediately() {
        let h = HeaderHandle::from_header(Header::default());
        assert!(h.is_set());
        let got = h.try_get().expect("set").expect("ok");
        assert_eq!(got, &Header::default());
    }

    #[test]
    fn set_then_try_get_returns_ok() {
        let h = HeaderHandle::new();
        h.set(Header::default()).expect("first set");
        assert!(h.is_set());
        let got = h.try_get().expect("set").expect("ok");
        assert_eq!(got, &Header::default());
    }

    #[test]
    fn poison_then_try_get_returns_err() {
        let h = HeaderHandle::new();
        h.poison(io::Error::other("aligner crashed")).expect("first poison");
        assert!(h.is_set());
        let err = h.try_get().expect("set").expect_err("poisoned");
        assert_eq!(err.kind(), io::ErrorKind::Other);
        assert_eq!(err.to_string(), "aligner crashed");
    }

    #[test]
    fn second_set_returns_already_set() {
        let h = HeaderHandle::new();
        h.set(Header::default()).expect("first");
        let err = h.set(Header::default()).expect_err("second set");
        let _ = err; // type asserted by the binding
    }

    #[test]
    fn set_then_poison_returns_already_set() {
        let h = HeaderHandle::new();
        h.set(Header::default()).expect("first");
        let _err: AlreadySetError =
            h.poison(io::Error::other("late")).expect_err("poison after set");
    }

    #[test]
    fn poison_then_set_returns_already_set() {
        let h = HeaderHandle::new();
        h.poison(io::Error::other("first")).expect("first");
        let _err: AlreadySetError = h.set(Header::default()).expect_err("set after poison");
    }

    #[test]
    fn clones_share_state() {
        let a = HeaderHandle::new();
        let b = a.clone();
        assert!(!a.is_set() && !b.is_set());
        a.set(Header::default()).expect("first");
        assert!(a.is_set() && b.is_set());
        assert!(b.try_get().expect("set").is_ok());
    }

    #[test]
    fn poison_preserves_error_kind_and_message_across_calls() {
        let h = HeaderHandle::new();
        h.poison(io::Error::new(io::ErrorKind::BrokenPipe, "stderr ring: foo")).unwrap();
        for _ in 0..3 {
            let e = h.try_get().unwrap().unwrap_err();
            assert_eq!(e.kind(), io::ErrorKind::BrokenPipe);
            assert_eq!(e.to_string(), "stderr ring: foo");
        }
    }

    #[test]
    fn set_from_another_thread_is_observable() {
        let h = HeaderHandle::new();
        let h2 = h.clone();
        let join = thread::spawn(move || {
            h2.set(Header::default()).expect("first set");
        });
        join.join().expect("thread join");
        assert!(h.is_set());
        assert!(h.try_get().unwrap().is_ok());
    }

    #[test]
    fn from_header_then_set_returns_already_set() {
        let h = HeaderHandle::from_header(Header::default());
        let _err: AlreadySetError =
            h.set(Header::default()).expect_err("from_header consumes the slot");
    }

    #[test]
    fn from_header_then_poison_returns_already_set() {
        let h = HeaderHandle::from_header(Header::default());
        let _err: AlreadySetError =
            h.poison(io::Error::other("late")).expect_err("from_header consumes the slot");
    }

    // ─────────────────────────────────────────────────────────────────────
    // Orphaned-setter backstops.
    //
    // These pin the two escapes described in the module's "Producer-side
    // invariant" section, so the claim there is tested rather than asserted.
    // ─────────────────────────────────────────────────────────────────────

    /// Shared fixtures for the orphaned-setter tests: a producer that finishes
    /// without ever resolving its handle, and a consumer shaped like the real
    /// BGZF writer (holds records until the header resolves; treats "input
    /// drained before the header arrived" as a hard error).
    mod orphan {
        use super::HeaderHandle;
        use std::io;

        use crate::item::HeapSize;
        use crate::outputs::Single;
        use crate::queues::QueueSpec;
        use crate::reorder::BranchOrdering;
        use crate::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};

        /// Item with a real heap payload, so a `ByteBounded` transport holding
        /// one reports non-zero `current_bytes` (the queue's accounting is
        /// `heap_size()`-only — `size_of::<T>()` is not counted).
        pub struct Block {
            pub payload: Vec<u8>,
        }
        impl HeapSize for Block {
            fn heap_size(&self) -> usize {
                self.payload.capacity()
            }
        }

        /// Producer that emits `remaining` blocks and then finishes, never
        /// calling `set`/`poison` on the handle it owns.
        #[derive(Clone)]
        pub struct OrphanedSetter {
            pub remaining: u32,
            pub block_bytes: usize,
            pub spec: QueueSpec,
            pub _handle: HeaderHandle,
        }
        impl Step for OrphanedSetter {
            type Input = ();
            type Outputs = Single<Block>;
            fn profile(&self) -> StepProfile {
                StepProfile {
                    name: "OrphanedSetter",
                    kind: StepKind::Exclusive,
                    sticky: false,
                    output_queues: vec![self.spec],
                    branch_ordering: vec![BranchOrdering::None],
                }
            }
            fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
                if self.remaining == 0 {
                    // Finishes WITHOUT resolving `_handle` — the orphan case.
                    return Ok(StepOutcome::Finished);
                }
                self.remaining -= 1;
                // Asserted, not discarded: backstop 2's whole premise is that
                // these items are stranded IN the queue. A silently dropped push
                // would leave `in_flight_bytes` at 0 and the test would pass for
                // the wrong reason. Budgets below always admit every block.
                ctx.outputs
                    .push(Block { payload: vec![0u8; self.block_bytes] })
                    .map_err(|_| io::Error::other("orphan test budget must admit every block"))?;
                Ok(StepOutcome::Progress)
            }
        }

        /// Consumer mirroring the real `WriteBgzfFile` sink.
        #[derive(Clone)]
        pub struct HeaderBlockedSink {
            pub handle: HeaderHandle,
        }
        impl Step for HeaderBlockedSink {
            type Input = Block;
            type Outputs = ();
            fn profile(&self) -> StepProfile {
                StepProfile {
                    name: "HeaderBlockedSink",
                    kind: StepKind::Exclusive,
                    sticky: false,
                    output_queues: vec![],
                    branch_ordering: vec![],
                }
            }
            fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
                // The header gates consumption: nothing may be written before
                // it resolves, so a blocked sink pops nothing.
                let header_ready = self.handle.try_get().transpose()?.is_some();
                if header_ready && ctx.input.pop().is_some() {
                    return Ok(StepOutcome::Progress);
                }
                if ctx.input.is_drained() {
                    if !header_ready {
                        return Err(io::Error::other(
                            "HeaderBlockedSink: input drained before HeaderHandle was resolved",
                        ));
                    }
                    return Ok(StepOutcome::Finished);
                }
                Ok(StepOutcome::NoProgress)
            }
        }
    }

    /// Backstop 1 — orphaned setter that emitted **nothing**: the consumer's
    /// `is_drained()` arm fires and the run fails fast with a diagnostic.
    ///
    /// `is_drained()` is `queue.is_drained() && queue.is_empty()`, so this arm
    /// is only reachable when the orphaned producer left no items behind. That
    /// precondition is what backstop 2 covers.
    #[test]
    fn orphaned_setter_with_empty_output_fails_the_run_fast() {
        use std::sync::mpsc;
        use std::time::Duration;

        use crate::queues::QueueSpec;
        use crate::{Pipeline, PipelineConfig};

        let handle = HeaderHandle::new();
        let for_run = handle.clone();
        let (tx, rx) = mpsc::channel();
        let join = thread::spawn(move || {
            let builder = Pipeline::builder();
            builder
                .chain(orphan::OrphanedSetter {
                    remaining: 0,
                    block_bytes: 0,
                    spec: QueueSpec::CountBounded { capacity: 8 },
                    _handle: for_run.clone(),
                })
                .chain(orphan::HeaderBlockedSink { handle: for_run })
                .into_sink_marker();
            let pipeline = builder.build().expect("pipeline build");
            let _ = tx.send(pipeline.run(PipelineConfig { threads: 2, ..Default::default() }));
        });

        let result = rx
            .recv_timeout(Duration::from_secs(10))
            .expect("an orphaned setter with an empty output must not hang the run");
        join.join().expect("run thread panicked");

        let err = result.expect_err("an unresolved HeaderHandle must fail the run");
        assert!(
            format!("{err}").contains("input drained before HeaderHandle was resolved"),
            "run must fail with the sink's drained-without-header error, got: {err}"
        );
        assert!(!handle.is_set(), "nobody ever resolved the handle");
    }

    /// Backstop 2 — orphaned setter that emitted items: those items sit unpopped
    /// (the blocked sink never consumes), so `is_drained()` stays false forever
    /// and the drained arm above can never fire. The deadlock monitor is what
    /// catches this: the stranded items are in flight on a `ByteBounded` edge,
    /// so `in_flight_bytes` is non-zero, `classify_stall` returns `Wedged` past
    /// the fatal timeout, and the run fails with `PipelineError::TimedOut`.
    ///
    /// The item must carry real heap bytes for that to hold — queue accounting
    /// is `heap_size()`-only, so a zero-heap item type would leave
    /// `in_flight_bytes` at 0, yielding `Starving`, which resets the stall clock
    /// on every poll and leaves the wedge uncatchable.
    #[test]
    fn orphaned_setter_with_inflight_bytes_is_caught_by_the_deadlock_monitor() {
        use std::sync::mpsc;
        use std::time::Duration;

        use crate::queues::QueueSpec;
        use crate::signal::PipelineError;
        use crate::{Pipeline, PipelineConfig};

        let handle = HeaderHandle::new();
        let for_run = handle.clone();
        let (tx, rx) = mpsc::channel();
        let join = thread::spawn(move || {
            let builder = Pipeline::builder();
            builder
                .chain(orphan::OrphanedSetter {
                    remaining: 4,
                    block_bytes: 4096,
                    spec: QueueSpec::ByteBounded { limit_bytes: 1 << 20 },
                    _handle: for_run.clone(),
                })
                .chain(orphan::HeaderBlockedSink { handle: for_run })
                .into_sink_marker();
            let pipeline = builder.build().expect("pipeline build");
            let stats = pipeline.stats();
            let _ = tx.send(pipeline.run(PipelineConfig {
                threads: 2,
                stats: Some(stats),
                // warn at 1s, fatal at 1s * DEADLOCK_FATAL_MULTIPLE.
                deadlock_timeout_secs: 1,
                ..Default::default()
            }));
        });

        let result = rx
            .recv_timeout(Duration::from_secs(60))
            .expect("the deadlock monitor must fail a header-wedged run, not let it hang");
        join.join().expect("run thread panicked");

        let err = result.expect_err("a header-wedged run must fail");
        assert!(
            matches!(err, PipelineError::TimedOut { .. }),
            "the wedge must surface as TimedOut, got: {err:?}"
        );
        assert!(!handle.is_set(), "nobody ever resolved the handle");
    }

    #[test]
    fn poison_visible_to_all_clones() {
        let a = HeaderHandle::new();
        let b = a.clone();
        let c = a.clone();
        a.poison(io::Error::new(io::ErrorKind::BrokenPipe, "boom")).unwrap();
        for clone in [&a, &b, &c] {
            let err = clone.try_get().expect("set").expect_err("poisoned");
            assert_eq!(err.kind(), io::ErrorKind::BrokenPipe);
            assert_eq!(err.to_string(), "boom");
        }
    }
}

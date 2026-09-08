//! A `StepKind::Serial` step that consumes one ordered stream of items and
//! performs a single finishing action when the stream is exhausted. Framework-
//! generic: this crate has zero dependency on the root `fgumi` package (only
//! ahash/anyhow/crossbeam-queue/log/noodles/parking_lot), so `Item`/`Self` are
//! free to be application-specific types defined outside this crate — see
//! `src/lib/inline_metrics_collector.rs::CoordinateGroupCollector` for the
//! concrete instance this was built for.
//!
//! **Serial-step trait shape**, confirmed by reading `GroupByMi`
//! (`src/lib/pipeline/steps/group/mi.rs`) and the no-output sink
//! `WriteRawFile` (`crates/fgumi-pipeline-io/src/sink/write_raw.rs`): a
//! `Serial` step is a plain [`Step`] impl — the framework itself holds the
//! single shared instance behind a `Mutex` and serializes `try_run` calls, so
//! the step body needs no internal locking of its own (unlike `WriteRawFile`,
//! which additionally wraps its writer in a `Mutex` because it is also
//! `sticky` — that is orthogonal to `Serial` and not needed here). A step
//! with no output branch (this one — it is always the terminal step on its
//! branch, never wired into a downstream queue) declares `type Outputs = ()`
//! with empty `output_queues`/`branch_ordering` vectors in its
//! [`StepProfile`], exactly as `WriteRawFile` does.

use anyhow::Result;

use crate::item::HeapSize;
use crate::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};

/// Consumes items one at a time via `record`, in the order they arrive, then
/// performs one finishing action via `finish` when the upstream input is
/// exhausted. Implementations typically buffer partial state across `record`
/// calls (e.g. group-by-key accumulation) and flush it in `finish`.
pub trait MetricsReducer: Send + 'static {
    /// The item type popped off the ordered input stream. Bounded by
    /// [`HeapSize`] (not just `Send + 'static`) because [`MetricsCollectorStep`]
    /// wires it in as a [`Step::Input`], and every step's input type must
    /// satisfy that bound to flow through the framework's queues.
    type Item: Send + HeapSize + 'static;

    /// Record one item. Called once per item, in the order the upstream
    /// producer emitted them.
    ///
    /// # Errors
    ///
    /// Returns any error raised while folding `item` into this reducer's
    /// buffered state.
    fn record(&mut self, item: Self::Item) -> Result<()>;

    /// Consume `self` and perform the one finishing action (e.g. flush any
    /// still-buffered state, then hand the result to a caller-supplied
    /// callback). Called exactly once, after the upstream input is fully
    /// drained.
    ///
    /// # Errors
    ///
    /// Returns any error raised while finishing.
    fn finish(self) -> Result<()>;
}

/// Wraps a [`MetricsReducer`] as a `StepKind::Serial` pipeline step with no
/// output branch — it is always the terminal step on whichever branch feeds
/// it.
///
/// The reducer is held as `Option<M>` (rather than a bare `M`, which the
/// design sketch this was built from used) because [`MetricsReducer::finish`]
/// consumes `self` by value: `try_run` takes `&mut self`, so the only way to
/// hand the reducer to `finish` is to take it out of the `Option` on the
/// drained path, leaving `None` behind as the step's permanently-finished
/// state.
pub struct MetricsCollectorStep<M: MetricsReducer> {
    reducer: Option<M>,
}

impl<M: MetricsReducer> MetricsCollectorStep<M> {
    #[must_use]
    pub fn new(reducer: M) -> Self {
        Self { reducer: Some(reducer) }
    }
}

impl<M: MetricsReducer> Step for MetricsCollectorStep<M> {
    type Input = M::Item;
    type Outputs = ();

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "MetricsCollectorStep",
            kind: StepKind::Serial,
            sticky: false,
            output_queues: vec![],
            branch_ordering: vec![],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> std::io::Result<StepOutcome> {
        // Already finished (a prior call took `reducer` and called `finish`).
        // Only reachable if the framework calls `try_run` again after
        // `Finished` — defensive, matching `WriteRawFile`'s own `guard.take()`
        // short-circuit.
        let Some(reducer) = self.reducer.as_mut() else {
            return Ok(StepOutcome::Finished);
        };

        if let Some(item) = ctx.input.pop() {
            reducer.record(item).map_err(std::io::Error::other)?;
            return Ok(StepOutcome::Progress);
        }

        if ctx.input.is_drained() {
            let reducer = self.reducer.take().expect("checked Some above");
            reducer.finish().map_err(std::io::Error::other)?;
            return Ok(StepOutcome::Finished);
        }

        Ok(StepOutcome::NoProgress)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct SumReducer {
        sum: u64,
        on_finish: Box<dyn FnOnce(u64) + Send>,
    }

    impl MetricsReducer for SumReducer {
        type Item = u64;
        fn record(&mut self, item: u64) -> anyhow::Result<()> {
            self.sum += item;
            Ok(())
        }
        fn finish(self) -> anyhow::Result<()> {
            (self.on_finish)(self.sum);
            Ok(())
        }
    }

    #[test]
    fn record_accumulates_and_finish_fires_exactly_once_with_the_final_value() {
        let observed = std::sync::Arc::new(std::sync::Mutex::new(None));
        let observed_clone = std::sync::Arc::clone(&observed);
        let mut reducer = SumReducer {
            sum: 0,
            on_finish: Box::new(move |sum| *observed_clone.lock().unwrap() = Some(sum)),
        };
        reducer.record(3).unwrap();
        reducer.record(4).unwrap();
        reducer.record(5).unwrap();
        reducer.finish().unwrap();
        assert_eq!(*observed.lock().unwrap(), Some(12));
    }

    /// `MetricsCollectorStep::new` wraps the reducer so it can later be taken
    /// by value for `finish`; this pins that construction alone doesn't lose
    /// or double-wrap it.
    #[test]
    fn new_wraps_the_reducer_as_present() {
        let step = MetricsCollectorStep::new(SumReducer { sum: 0, on_finish: Box::new(|_| {}) });
        assert!(step.reducer.is_some());
    }

    #[test]
    fn profile_advertises_serial_with_no_output_branch() {
        let step = MetricsCollectorStep::new(SumReducer { sum: 0, on_finish: Box::new(|_| {}) });
        let p = step.profile();
        assert_eq!(p.name, "MetricsCollectorStep");
        assert_eq!(p.kind, StepKind::Serial);
        assert!(!p.sticky);
        assert!(p.output_queues.is_empty());
        assert!(p.branch_ordering.is_empty());
    }
}

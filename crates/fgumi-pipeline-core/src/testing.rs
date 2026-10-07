//! Drive one [`Step`] by hand, outside a pipeline.
//!
//! For unit tests of a step's `try_run` contract — admission, held output, the
//! drain hand-off — where building and running a whole pipeline would hide
//! which call returned what. Not used by the runtime.

use std::io;
use std::sync::Arc;

use crate::builder::InstrumentationLevel;
use crate::handles::{BranchInputHandle, OutputQueueSet};
use crate::item::HeapSize;
use crate::outputs::StepOutputs;
use crate::queues::{ItemQueue, UnboundedQueue};
use crate::runtime::contexts::StepCounters;
use crate::step::{InputHandle, OutputHandles, Step, StepCtx, StepOutcome};

/// One step's input queue and output queues, wired the way the builder would
/// wire them (output queues from the step's own [`crate::StepProfile`]), with
/// no scheduler: the test calls [`StepProbe::try_run`] itself.
pub struct StepProbe<S: Step> {
    input_queue: Arc<UnboundedQueue<S::Input>>,
    input: BranchInputHandle<S::Input>,
    outputs: OutputHandles<S::Outputs>,
    output_queues: OutputQueueSet,
    counters: StepCounters,
}

impl<S: Step> StepProbe<S> {
    /// Build the queues for `step` (its input is unbounded; its outputs follow
    /// its profile).
    #[must_use]
    pub fn new(step: &S) -> Self {
        let profile = step.profile();
        let (output_queues, view) = <S::Outputs as StepOutputs>::build_queues(
            &profile.output_queues,
            &profile.branch_ordering,
            InstrumentationLevel::Off,
        );
        let input_queue = Arc::new(UnboundedQueue::new());
        let input =
            BranchInputHandle::direct(Arc::clone(&input_queue) as Arc<dyn ItemQueue<S::Input>>);
        Self {
            input_queue,
            input,
            outputs: OutputHandles::new(view),
            output_queues,
            counters: StepCounters::disabled(),
        }
    }

    /// Queue one input item.
    ///
    /// # Panics
    ///
    /// Panics if the input was already closed with [`Self::close_input`].
    pub fn push_input(&self, item: S::Input) {
        assert!(self.input_queue.try_push(item).is_ok(), "unbounded input never rejects");
    }

    /// Mark the input drained (upstream finished).
    pub fn close_input(&self) {
        self.input_queue.mark_drained();
    }

    /// Whether the step's input queue is empty (no item left to pop).
    #[must_use]
    pub fn input_is_empty(&self) -> bool {
        self.input.is_empty()
    }

    /// The consumer end of output branch `branch`, typed as `T`.
    ///
    /// # Panics
    ///
    /// Panics if `T` is not that branch's item type or the branch was taken.
    pub fn take_output<T: Send + HeapSize + 'static>(
        &mut self,
        branch: usize,
    ) -> BranchInputHandle<T> {
        self.output_queues.take_typed_input(branch)
    }

    /// Call `step.try_run` once against these queues.
    ///
    /// # Errors
    ///
    /// Whatever the step's `try_run` returns.
    pub fn try_run(&self, step: &mut S) -> io::Result<StepOutcome> {
        let mut ctx =
            StepCtx { input: &self.input, outputs: &self.outputs, counters: &self.counters };
        step.try_run(&mut ctx)
    }
}

/// Pin the admission contract every capped pool step owes (see
/// [`crate::admit_input`]), against `make`'s step and one input `item`, under a
/// one-permit cap named `cap_name`:
///
/// - a full cap refuses (`Capped`, counted) without popping the item;
/// - once the permit is released, the step pops and processes it, and releases
///   its own permit when `try_run` returns;
/// - an empty poll, open or drained, takes no permit and is not a refusal;
/// - a worker copy shares the same counter (the cap belongs to the phase);
/// - `make(None)` is uncapped.
///
/// # Panics
///
/// Panics if the step breaks any part of the contract.
pub fn assert_admission_contract<S: Step>(
    cap_name: &'static str,
    make: impl Fn(Option<Arc<crate::PhaseCap>>) -> S,
    item: S::Input,
) {
    let cap = crate::PhaseCap::new(cap_name, 1);
    let mut step = make(Some(Arc::clone(&cap)));
    let probe = StepProbe::new(&step);
    probe.push_input(item);

    let held = cap.try_acquire().expect("the only permit");
    assert_eq!(probe.try_run(&mut step).expect("try_run"), StepOutcome::Capped);
    assert!(!probe.input_is_empty(), "a refused clone must not pop the item");
    assert_eq!(cap.refused(), 1, "the refusal is counted");

    drop(held);
    assert_eq!(probe.try_run(&mut step).expect("try_run"), StepOutcome::Progress);
    assert!(probe.input_is_empty(), "the admitted clone pops the item");
    assert_eq!(cap.active(), 0, "the permit is released when try_run returns");

    // Empty polls take no permit: with the only permit held elsewhere, an open
    // empty input reports NoProgress (not Capped) and a drained one Finished
    // (an admit-first step would report Capped for both), and neither is a
    // refusal.
    let held = cap.try_acquire().expect("the only permit");
    assert_eq!(probe.try_run(&mut step).expect("try_run"), StepOutcome::NoProgress);
    probe.close_input();
    assert_eq!(probe.try_run(&mut step).expect("try_run"), StepOutcome::Finished);
    assert_eq!(cap.refused(), 1, "empty polls are not refusals");

    let copy = step.new_worker_copy();
    assert!(
        std::ptr::eq(copy.phase_cap().expect("capped copy"), Arc::as_ptr(&cap)),
        "a worker copy shares the step's counter"
    );
    drop(held);
    assert!(make(None).phase_cap().is_none(), "uncapped by default");
}

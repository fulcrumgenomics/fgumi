//! Refill hints: how a stage asks the pool to refill its input first.
//!
//! A pipeline-level type, so both the steps that raise a hint (`steps`) and the
//! chain builder that installs the scheduler for it (`chains`) can name it
//! without `steps` importing from `chains`.

use std::sync::Arc;
use std::sync::atomic::AtomicBool;

use crate::pipeline::core::runtime::RefillSource;
use crate::pipeline::core::topology::{BranchIdx, StepIdx};

/// A stage's refill hint: the scheduler's [`RefillSource`] (while its signal is
/// up and the feeding producer's output queue holds less than its cap, the
/// pool walks the steps up to and including that producer upstream-first),
/// plus when the hint applies.
///
/// A chain may carry several hints, one per stage (the chain builder rejects a
/// stage adding a second). `requires_drain_first` hints apply only when the
/// chain chose drain-first scheduling (the in-process aligner's input refill);
/// the others apply under either base walk direction (the sort merge's
/// starvation).
#[derive(Clone)]
pub(crate) struct RefillHint {
    /// What the scheduler is handed when the hint applies.
    pub(crate) source: RefillSource,
    /// Whether the hint applies only when the chain chose drain-first dispatch.
    pub(crate) requires_drain_first: bool,
}

impl RefillHint {
    /// A hint raised on `signal` for the stage fed by `feed` (a producer step
    /// and output branch), capping the read-ahead on that queue at `cap_bytes`
    /// (`u64::MAX` = uncapped).
    pub(crate) fn new(
        signal: Arc<AtomicBool>,
        feed: (StepIdx, BranchIdx),
        cap_bytes: u64,
        requires_drain_first: bool,
    ) -> Self {
        Self {
            source: RefillSource {
                signal,
                refill_through: feed.0,
                refill_branch: feed.1,
                cap_bytes,
            },
            requires_drain_first,
        }
    }
}

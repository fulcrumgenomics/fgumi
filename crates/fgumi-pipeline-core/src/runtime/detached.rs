//! Detached-step extraction and the dedicated **driver thread** that runs them
//! off the work-stealing pool.
//!
//! A [`StepKind::Detached`] step is excluded
//! from the pool ([`build_worker_storage`](crate::runtime::build_worker_storage)
//! gives every pool worker a `Skip` entry for it) and instead runs on a dedicated
//! OS thread spawned at run start alongside the deadlock-monitor /
//! queue-rebalancer and joined after the workers. This is the legacy sort's
//! "N + 2" threading: N pool workers do the parallel (compression-bound) work
//! while off-pool driver threads do the serial coordination + I/O, so neither
//! steals a pool worker slot.
//!
//! ## The driver IS a 1-thread pool (`run_detached_driver`)
//!
//! A driver thread does **not** have a bespoke drive loop. It runs the *same*
//! `run_worker_loop` the pool uses — over a
//! purpose-built storage row where its group's steps are `Owned` and every other
//! step is `Skip` ([`build_driver_storage`]) — with a
//! [`WorkerCore::driver`](crate::runtime::WorkerCore) (Park backoff — a
//! 10 µs → 2 ms `park_timeout` ramp, ended early by the `unpark` every producer
//! of its live steps delivers through the wake plan — off-pool stats
//! attribution) and the [`DrainFirstScheduler`]. So a driver is literally a
//! 1-thread (or, for a group, still 1-thread over several `Owned` steps) instance
//! of the pool loop. Several detached steps sharing a
//! [`DetachedGroup::Shared`](crate::step::DetachedGroup) label are driven by ONE
//! thread that round-robins them; [`DetachedGroup::PerStep`] (the default) keeps
//! one thread per step. A driver group may also host one `Parallel` step's
//! single clone (`PoolPlacement::ExcludeReader` at one worker, on its Detached
//! producer's driver): [`extract_hosted_parallel_steps`] inserts it in chain
//! order, and its drain counter is 1.
//!
//! Because the group's steps and their phases are temporally disjoint (e.g. the
//! sort's phase-1 admit/sort/frame finish and leave the live set before the
//! phase-2 merge runs), one driver thread covers a whole phase's coordination
//! without oversubscription — exactly like main's single main thread.
//!
//! ## Lock-ordering acyclicity (non-negotiable, the deadlock proof)
//!
//! A driver thread must NEVER hold one queue's internal lock while parking on
//! another. `run_worker_loop` upholds that by construction:
//!
//!   1. Each `try_run_erased` is a single non-blocking call. The step body pops
//!      from its input transport (a lock-free `crossbeam ArrayQueue` — `try_pop`,
//!      no lock held across the call) and pushes to its outputs (`try_push`,
//!      likewise); a full output / empty input is reported back as `NoProgress` /
//!      `Contention`. It never *blocks* inside `try_run`.
//!   2. The only blocking a driver does is the loop's timer park
//!      (`park_timeout` of `WorkerCore::backoff_deadline`, the Park policy's
//!      ramp), which holds NO queue lock.
//!   3. The loop tries EVERY live step in a pass before it parks (round-robin,
//!      park only after a full no-progress pass). So when one grouped step is
//!      blocked, a sibling on the same driver still runs — a park-on-first-idle
//!      loop would wedge (see `driver_round_robins_all_live_before_parking`).
//!
//! So there is no cycle: a driver parks only *between* `try_run` calls, never
//! while holding a transport lock, so a two-sided step (consuming from the pool
//! AND producing to it, both bounded) cannot deadlock. The
//! `detached_two_sided_no_deadlock` test pins this with a wall-clock watchdog.

use std::any::Any;
use std::collections::HashMap;
use std::sync::Arc;

use crate::erased::{ErasedStep, ErasedStepCtx};
use crate::runtime::contexts::ChainContexts;
use crate::runtime::drain::StepDrainCounter;
use crate::runtime::driver::run_worker_loop;
use crate::runtime::placement::ParallelHosts;
use crate::runtime::scheduler::DrainFirstScheduler;
use crate::runtime::stats::PipelineStats;
use crate::runtime::storage::WorkerStepEntry;
use crate::runtime::worker_core::WorkerCore;
use crate::signal::PipelineSignal;
use crate::step::{
    Affinity, DetachedGroup, OutputsViewAny, PoolPlacement, StepKind, StepOutcome, StepProfile,
};
use crate::topology::StepIdx;

/// Sentinel left in the chain's `steps` vec in place of a step whose real
/// instance moved onto a driver thread: a `Detached` step
/// ([`extract_detached_steps`]) or a `Parallel` step's single hosted clone
/// ([`extract_hosted_parallel_steps`]). Keeping a same-position placeholder
/// preserves the `step_idx`-aligned indexing that `build_worker_storage` and
/// `ChainContexts` rely on. It reports the moved step's cached `name`, `kind`,
/// `detached_group` and `pool_placement` — `build_worker_storage` reads only
/// `kind()` (`Detached`, or a `Parallel` step with no host worker: every worker
/// gets `Skip`) and then drops the box — and every other method panics, to
/// catch a framework bug if one is ever invoked.
struct ExtractedPlaceholder {
    name: &'static str,
    kind: StepKind,
    /// The real step's group, carried through the swap.
    ///
    /// Not consulted on the run path — `extract_detached_steps` reads the real
    /// step's group *before* installing this placeholder. Preserved because
    /// `extract_detached_steps` is `pub` and leaves these placeholders in the
    /// caller's slice: an external caller regrouping from that slice would read
    /// `PerStep` for a step that declared `Shared(..)` and split one shared
    /// group across separate driver threads. `profile()` panics for the same
    /// class of reason — returning fabricated metadata misinforms that caller.
    group: DetachedGroup,
    /// The real step's placement, carried through the swap like `group`.
    placement: PoolPlacement,
}

impl ExtractedPlaceholder {
    /// The placeholder of `real`, which is about to move onto a driver thread.
    fn of(real: &dyn ErasedStep) -> Self {
        Self {
            name: real.name(),
            kind: real.kind(),
            group: real.detached_group(),
            placement: real.pool_placement(),
        }
    }

    fn moved(&self, method: &str) -> ! {
        panic!(
            "ExtractedPlaceholder::{method} invoked for {:?} ({:?}) — the real instance runs on \
             its driver thread",
            self.name, self.kind
        );
    }
}

impl ErasedStep for ExtractedPlaceholder {
    fn profile(&self) -> StepProfile {
        // Panics like every other placeholder method rather than returning empty
        // `output_queues` / `branch_ordering`. `Pipeline::run` computes the drain
        // counters AFTER the steps are extracted, from the cached `kind()`
        // accessor, and everything else after extraction reads the cached
        // accessors too — so nothing on the run path calls this. But
        // `extract_detached_steps` is `pub` and leaves placeholders in the
        // caller's slice: silently reporting "no outputs" for a step that
        // declares them would misinform any later caller.
        panic!(
            "ExtractedPlaceholder::profile invoked for {:?} — the real instance was extracted; \
             read the cached kind()/name() accessors instead",
            self.name
        );
    }
    fn name(&self) -> &'static str {
        self.name
    }
    fn kind(&self) -> StepKind {
        self.kind
    }
    fn sticky(&self) -> bool {
        false
    }
    fn affinity(&self) -> Affinity {
        Affinity::None
    }
    fn detached_group(&self) -> DetachedGroup {
        // See the field doc: off the run path, but must not report `PerStep` for a
        // step that declared `Shared(..)`.
        self.group
    }
    fn pool_placement(&self) -> PoolPlacement {
        self.placement
    }
    fn try_run_erased(&mut self, _ctx: &mut ErasedStepCtx<'_>) -> std::io::Result<StepOutcome> {
        self.moved("try_run_erased")
    }
    fn clone_boxed(&self) -> Box<dyn ErasedStep> {
        self.moved("clone_boxed")
    }
    fn build_input_handle(
        &self,
        _producer_set: &mut crate::handles::OutputQueueSet,
        _branch_idx: usize,
    ) -> Box<dyn Any + Send + Sync> {
        self.moved("build_input_handle")
    }
    fn build_output_set(
        &self,
        _level: crate::builder::InstrumentationLevel,
    ) -> (crate::handles::OutputQueueSet, OutputsViewAny) {
        self.moved("build_output_set")
    }
    fn build_fused_output_set(
        &self,
        _level: crate::builder::InstrumentationLevel,
    ) -> (crate::handles::OutputQueueSet, OutputsViewAny) {
        self.moved("build_fused_output_set")
    }
    fn wrap_outputs_view(&self, _view: OutputsViewAny) -> Box<dyn Any + Send + Sync> {
        self.moved("wrap_outputs_view")
    }
    fn mark_outputs_drained(&self, _outputs: &(dyn Any + Send + Sync)) {
        self.moved("mark_outputs_drained")
    }
    fn is_source(&self) -> bool {
        false
    }
}

/// One dedicated driver thread's worth of extracted steps, in chain
/// (`StepIdx`) order. The caller spawns one OS thread per group and drives it
/// with `run_detached_driver`.
pub struct DetachedDriverGroup {
    /// Members in chain order: the group's `Detached` steps, plus any
    /// `Parallel` step whose single clone is hosted on this driver (added by
    /// [`Self::insert_hosted`] after [`Self::new`]). Never empty, and the first
    /// member is always `Detached` (a hosted step comes after its Detached
    /// producer). Private so those invariants cannot be bypassed by an external
    /// caller building the struct directly (which could otherwise trigger
    /// `steps[0]` panics in `primary_step` / `label`, or name the driver after a
    /// hosted step).
    steps: Vec<(StepIdx, Box<dyn ErasedStep>)>,
}

impl DetachedDriverGroup {
    /// Wrap a driver thread's extracted Detached steps, enforcing the invariants
    /// every consumer relies on: the group is **non-empty** (`primary_step` /
    /// `label` index `steps[0]`) and **every step is [`StepKind::Detached`]** (a
    /// non-detached step must not run off-pool on a dedicated driver thread). A
    /// hosted `Parallel` clone joins later, through [`Self::insert_hosted`].
    ///
    /// # Panics
    ///
    /// Panics if `steps` is empty or contains a non-`Detached` step.
    #[must_use]
    fn new(steps: Vec<(StepIdx, Box<dyn ErasedStep>)>) -> Self {
        assert!(!steps.is_empty(), "detached driver group must be non-empty");
        assert!(
            steps.iter().all(|(_, step)| step.kind() == StepKind::Detached),
            "detached driver group may only contain Detached steps"
        );
        Self { steps }
    }

    /// Add a `Parallel` step whose single clone `plan_parallel_hosts` placed on
    /// this driver, at its chain position.
    ///
    /// # Panics
    ///
    /// Panics if `step` is not [`StepKind::Parallel`], or if it would precede
    /// every member (the first member must stay the group's `Detached` step: it
    /// names the driver thread and keys its stats; `PipelineBuilder` registers
    /// a consumer after its producer, so a hosted step always follows it).
    pub fn insert_hosted(&mut self, idx: StepIdx, step: Box<dyn ErasedStep>) {
        assert_eq!(
            step.kind(),
            StepKind::Parallel,
            "only a Parallel step's clone is hosted on a driver (`{}`)",
            step.name()
        );
        let at = self.steps.partition_point(|(i, _)| i.0 < idx.0);
        assert!(
            at > 0,
            "hosted Parallel step `{}` (step {}) must follow its Detached producer in the \
             driver group",
            step.name(),
            idx.0
        );
        self.steps.insert(at, (idx, step));
    }

    /// The group's step indices in chain order, hosted members included.
    #[cfg(test)]
    pub fn step_indices(&self) -> impl Iterator<Item = StepIdx> + '_ {
        self.steps.iter().map(|(idx, _)| *idx)
    }

    /// The group's representative step (first in chain order). Used as the
    /// driver thread's off-pool stats key (`WorkerCore::driver`) and as a stable
    /// label. Non-empty by construction.
    #[must_use]
    pub fn primary_step(&self) -> StepIdx {
        self.steps[0].0
    }

    /// The name of the group's representative step (first in chain order), used
    /// to label a `PerStep` driver thread. Non-empty by construction.
    #[must_use]
    pub fn primary_name(&self) -> &'static str {
        self.steps[0].1.name()
    }

    /// The group label — derived from the steps (every step in the group reports
    /// the same [`DetachedGroup`]), used to name the driver thread. Non-empty by
    /// construction.
    #[must_use]
    pub fn label(&self) -> DetachedGroup {
        self.steps[0].1.detached_group()
    }

    /// Consume the group, yielding its members for [`build_driver_storage`].
    #[must_use]
    fn into_steps(self) -> Vec<(StepIdx, Box<dyn ErasedStep>)> {
        self.steps
    }
}

/// The label a driver is known by: its [`DetachedGroup::Shared`] label, or for
/// a [`DetachedGroup::PerStep`] driver the name of its one Detached step.
/// `group` and `primary_name` are the group's first step's (a hosted step never
/// comes first). `Pipeline::run` names the driver thread from it
/// ([`driver_thread_name`]) and `dag_at` renders it as `host=driver(<label>)`,
/// so the two cannot disagree.
#[must_use]
pub(crate) fn driver_label(group: DetachedGroup, primary_name: &'static str) -> &'static str {
    match group {
        DetachedGroup::Shared(label) => label,
        DetachedGroup::PerStep => primary_name,
    }
}

/// The driver thread's OS name: `fgumi-driver-<label>` for a shared group,
/// `fgumi-detached-<label>` for a per-step driver ([`driver_label`]).
#[must_use]
pub(crate) fn driver_thread_name(group: DetachedGroup, primary_name: &'static str) -> String {
    let label = driver_label(group, primary_name);
    match group {
        DetachedGroup::Shared(_) => format!("fgumi-driver-{label}"),
        DetachedGroup::PerStep => format!("fgumi-detached-{label}"),
    }
}

/// Move each Parallel step whose hosts name a driver into that driver's group
/// ([`DetachedDriverGroup::insert_hosted`]), leaving an `ExtractedPlaceholder`
/// so the surviving slots keep their `step_idx` positions. Called by
/// `Pipeline::run` after [`extract_detached_steps`] and `plan_parallel_hosts`,
/// before `build_worker_storage` consumes `steps`.
///
/// # Panics
///
/// Panics if `hosts` does not cover every step, if a host names a driver index
/// with no group, or as [`DetachedDriverGroup::insert_hosted`] does (the step
/// is not `Parallel`, or would precede its group's Detached steps).
pub fn extract_hosted_parallel_steps(
    steps: &mut [Box<dyn ErasedStep>],
    hosts: &[ParallelHosts],
    groups: &mut [DetachedDriverGroup],
) {
    assert_eq!(hosts.len(), steps.len(), "hosts must cover every step");
    for (idx, h) in hosts.iter().enumerate() {
        let Some(d) = h.driver else { continue };
        let placeholder: Box<dyn ErasedStep> =
            Box::new(ExtractedPlaceholder::of(steps[idx].as_ref()));
        let real = std::mem::replace(&mut steps[idx], placeholder);
        groups
            .get_mut(d.0)
            .expect("host driver index has a group")
            .insert_hosted(StepIdx(idx), real);
    }
}

/// Remove every [`StepKind::Detached`] step's
/// real instance from `steps`, replacing each in place with an
/// `ExtractedPlaceholder` so the surviving slots keep their `step_idx`
/// positions (which `build_worker_storage` and `ChainContexts` index by).
///
/// Groups the extracted steps by their [`DetachedGroup`]: every
/// [`DetachedGroup::Shared`] label collects onto ONE group (one driver thread);
/// each [`DetachedGroup::PerStep`] step becomes its own singleton group (the
/// legacy one-thread-per-step behavior — the default, so non-sort chains are
/// unchanged). Within a group and across groups, order follows chain order
/// (first appearance). Called by `Pipeline::run` **before** `build_worker_storage`
/// consumes `steps`, while the (read-only) `ChainContexts` have already been
/// built from `&steps`.
#[must_use]
pub fn extract_detached_steps(steps: &mut [Box<dyn ErasedStep>]) -> Vec<DetachedDriverGroup> {
    // Accumulate each driver thread's steps as a raw vec, then wrap through
    // `DetachedDriverGroup::new` so the non-empty / all-Detached invariants are
    // enforced in one place rather than trusting each construction site. The
    // grouping itself is `driver_index_of`'s, so the wake plan's routing and
    // the driver threads cannot disagree.
    let driver_of = driver_index_of(steps);
    let n_groups = driver_of.iter().flatten().map(|d| d.0 + 1).max().unwrap_or(0);
    let mut group_steps: Vec<Vec<(StepIdx, Box<dyn ErasedStep>)>> =
        (0..n_groups).map(|_| Vec::new()).collect();
    for (idx, (slot, driver)) in steps.iter_mut().zip(driver_of).enumerate() {
        let Some(driver) = driver else { continue };
        let placeholder: Box<dyn ErasedStep> = Box::new(ExtractedPlaceholder::of(slot.as_ref()));
        let real = std::mem::replace(slot, placeholder);
        group_steps[driver.0].push((StepIdx(idx), real));
    }
    group_steps.into_iter().map(DetachedDriverGroup::new).collect()
}

/// The driver thread each step runs on: `Some(d)` for a [`StepKind::Detached`]
/// step, `None` otherwise. Every [`DetachedGroup::Shared`] label maps to ONE
/// driver; each [`DetachedGroup::PerStep`] step gets its own. Drivers are
/// numbered by first appearance in chain order — the order
/// [`extract_detached_steps`] returns its groups in.
#[must_use]
pub(crate) fn driver_index_of(
    steps: &[Box<dyn ErasedStep>],
) -> Vec<Option<crate::runtime::wake::DriverIdx>> {
    use crate::runtime::wake::DriverIdx;
    let mut shared_index: HashMap<&'static str, usize> = HashMap::new();
    let mut n_drivers = 0usize;
    steps
        .iter()
        .map(|step| {
            if step.kind() != StepKind::Detached {
                return None;
            }
            let d = match step.detached_group() {
                DetachedGroup::PerStep => {
                    n_drivers += 1;
                    n_drivers - 1
                }
                DetachedGroup::Shared(label) => *shared_index.entry(label).or_insert_with(|| {
                    n_drivers += 1;
                    n_drivers - 1
                }),
            };
            Some(DriverIdx(d))
        })
        .collect()
}

/// Build a driver thread's storage row: a full-length `Vec<WorkerStepEntry>`
/// (one entry per step of `drain_counters`, indexed by global `step_idx` like
/// every other row) where each of `group_steps` is `Owned` and every other slot
/// is `Skip`. The driver thread runs `run_worker_loop` over this row exactly as
/// a pool worker runs over its own row.
///
/// # Panics
///
/// - if a group step's `StepDrainCounter` is not 1. The driver runs exactly one
///   instance of each member, so only a counter of 1 reaches 0 and closes the
///   member's output: a Detached step's always is, and a `Parallel` step's is
///   only when `plan_parallel_hosts` gave its single clone to this driver. A
///   Parallel step that also has pool clones counts them, and a driver running
///   one more instance would never take that counter to 0, hanging the
///   downstream consumer. Checked before the driver runs anything, while no
///   other thread can have touched the counter (every pool worker `Skip`s a
///   driver's members);
/// - if two group steps map to the same `step_idx` (a double registration), or a
///   step index is out of range.
#[must_use]
pub fn build_driver_storage(
    group_steps: Vec<(StepIdx, Box<dyn ErasedStep>)>,
    drain_counters: &[Arc<StepDrainCounter>],
) -> Vec<WorkerStepEntry> {
    let n_total_steps = drain_counters.len();
    let mut row: Vec<WorkerStepEntry> = (0..n_total_steps).map(|_| WorkerStepEntry::Skip).collect();
    for (idx, step) in group_steps {
        assert!(
            idx.0 < n_total_steps,
            "driver group step index {} out of range (chain has {n_total_steps} steps)",
            idx.0
        );
        let remaining = drain_counters[idx.0].remaining();
        assert_eq!(
            remaining,
            1,
            "driver group step `{}` ({:?}) has drain counter {remaining}, not 1; a driver runs \
             one instance of it, so its shared output would never close",
            step.name(),
            step.kind()
        );
        assert!(
            matches!(row[idx.0], WorkerStepEntry::Skip),
            "driver group step index {} registered twice (dual registration)",
            idx.0
        );
        row[idx.0] = WorkerStepEntry::Owned { step };
    }
    row
}

/// Drive one [`DetachedDriverGroup`] to completion on the calling (dedicated)
/// thread — the unified "1-thread pool". Builds the group's `Owned`/`Skip`
/// storage row ([`build_driver_storage`]) and runs the *same*
/// [`run_worker_loop`] the pool uses, with a [`WorkerCore::driver`] (Park
/// backoff, off-pool stats) and the [`DrainFirstScheduler`] (drain/seal
/// downstream before producing more — frees the sort's capacity-1 arena fastest).
///
/// `drain_counters` is the full per-step slice (init 1 for each detached step and
/// each hosted Parallel clone, so the single finisher closes its output edges —
/// the downstream consumer's end-of-stream signal). On a cancel before
/// `Finished`, outputs are NOT closed (the run is tearing down; the recorded
/// error/cancel is what propagates) — `run_worker_loop`'s top-of-loop `is_done`
/// break upholds this.
///
/// `board` / `state_slot` — the per-thread state board for scheduling telemetry
/// and this driver thread's slot in it; `None` (telemetry off) keeps stamping a
/// no-op. Detached drivers take slots after the pool workers (`threads..`), so
/// their slot never collides with a pool worker's.
#[allow(clippy::too_many_arguments)] // one driver's worth of shared state plus the telemetry board/slot
pub(crate) fn run_detached_driver(
    group: DetachedDriverGroup,
    contexts: &Arc<ChainContexts>,
    drain_counters: &[Arc<StepDrainCounter>],
    signal: &Arc<PipelineSignal>,
    stats: Option<&Arc<PipelineStats>>,
    liveness: &crate::liveness::LivenessCounter,
    board: Option<&crate::runtime::worker_state::WorkerStateBoard>,
    state_slot: usize,
    // The chain's wake plan; a driver is a notifier/unparker and is woken by
    // `unpark`, never an event-count waiter. Passed with `pinned = true` below
    // so the driver keeps its `Park` backoff and never deep-parks on the
    // shared event-count.
    wake: &crate::runtime::wake::WakePlan,
    // A test override of this driver's idle-backoff bounds
    // (`PipelineConfig::test_backoff`); `None` in every real run.
    backoff_override: Option<u64>,
) {
    let primary = group.primary_step();
    debug_assert_eq!(drain_counters.len(), contexts.inputs.len(), "one counter per step");
    let mut row = build_driver_storage(group.into_steps(), drain_counters);
    let mut worker = WorkerCore::driver(primary).with_backoff_override(backoff_override);
    run_worker_loop(
        &mut worker,
        &mut row,
        contexts,
        drain_counters,
        signal,
        stats,
        liveness,
        &DrainFirstScheduler,
        board,
        state_slot,
        wake,
        // `pinned`: a detached driver is not a pool waiter. This routes its idle
        // to the `Park` backoff rather than the shared event-count wait, while
        // the plan above still routes its wakes.
        true,
    );
}

#[cfg(test)]
mod tests {
    use std::io;
    use std::sync::atomic::{AtomicU32, Ordering};
    use std::time::Duration;

    use super::*;
    use crate::builder::InstrumentationLevel;
    use crate::erased::TypedStep;
    use crate::outputs::Single;
    use crate::queues::QueueSpec;
    use crate::reorder::BranchOrdering;
    use crate::step::{InputHandle, OutputHandles, Step, StepCtx, StepKind, StepProfile};

    /// `() -> u32` source stub used only so `build_output_set` constructs the
    /// transport that becomes the Detached step's input edge. Never run.
    #[derive(Clone)]
    struct SrcStub {
        capacity: usize,
    }
    impl Step for SrcStub {
        type Input = ();
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "Src",
                kind: StepKind::Exclusive,
                sticky: false,
                output_queues: vec![QueueSpec::CountBounded { capacity: self.capacity }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            Ok(StepOutcome::NoProgress)
        }
    }

    /// `u32 -> u32` pass-through Detached step: pop one item per `try_run`,
    /// push it on (holding it on output-full backpressure), report `Finished`
    /// once input drains and nothing is held.
    #[derive(Clone)]
    struct PassThroughDetached {
        held: Option<u32>,
    }
    impl Step for PassThroughDetached {
        type Input = u32;
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "PassThroughDetached",
                kind: StepKind::Detached,
                sticky: false,
                output_queues: vec![QueueSpec::CountBounded { capacity: 4 }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            // A rejected push reports `NoProgress`, not `Contention`: the driver
            // treats them identically, but `Contention` means "a Serial step's
            // mutex was held by another worker" and feeds `contention_count`,
            // which the bottleneck verdict turns into its SPIN ratio. Using it
            // for ordinary output backpressure invents contention that never
            // happened.
            if let Some(v) = self.held.take() {
                if ctx.outputs.push(v).is_err() {
                    self.held = Some(v);
                    return Ok(StepOutcome::NoProgress);
                }
                return Ok(StepOutcome::Progress);
            }
            match ctx.input.pop() {
                Some(v) => match ctx.outputs.push(v) {
                    Ok(()) => Ok(StepOutcome::Progress),
                    Err(unpushed) => {
                        self.held = Some(unpushed.into_item());
                        Ok(StepOutcome::NoProgress)
                    }
                },
                None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
                None => Ok(StepOutcome::NoProgress),
            }
        }
    }

    /// Assemble a `Source -> PassThroughDetached -> Sink` shaped context by
    /// hand: build the producer's output set (= the Detached input edge), wire
    /// the Detached step's input from it, build the Detached step's own output
    /// set (= the downstream consumer's input edge), and return the pieces.
    #[allow(clippy::type_complexity)]
    fn build_one_detached(
        src_capacity: usize,
    ) -> (
        Box<dyn ErasedStep>,                       // the detached step
        Arc<ChainContexts>,                        // contexts for step_idx 1
        Arc<Box<dyn std::any::Any + Send + Sync>>, // producer outputs (push side)
        crate::handles::BranchInputHandle<u32>,    // downstream consumer (pop side)
    ) {
        let producer: Box<dyn ErasedStep> =
            Box::new(TypedStep::new(SrcStub { capacity: src_capacity }));
        let (mut producer_set, producer_view) =
            producer.build_output_set(InstrumentationLevel::Off);
        let producer_outputs_any = producer.wrap_outputs_view(producer_view);

        let det: Box<dyn ErasedStep> = Box::new(TypedStep::new(PassThroughDetached { held: None }));
        let det_input = det.build_input_handle(&mut producer_set, 0);
        let (mut det_set, det_view) = det.build_output_set(InstrumentationLevel::Off);
        let det_outputs_any = det.wrap_outputs_view(det_view);
        let det_output_consumer = det_set.take_typed_input::<u32>(0);

        let contexts = Arc::new(ChainContexts {
            inputs: vec![Box::new(()), det_input, Box::new(())],
            outputs: vec![Box::new(()), det_outputs_any, Box::new(())],
            bounded_queues: vec![],
            holder_only_queues: vec![],
            edges: vec![],
            step_counters: (0..3).map(|_| crate::runtime::StepCounters::disabled()).collect(),
        });
        (det, contexts, Arc::new(producer_outputs_any), det_output_consumer)
    }

    /// Drive a single detached `step` at `step_idx` through the unified driver —
    /// a `PerStep` group of one — on the calling thread, with a full-length
    /// `drain_counters` slice (init 1 each). Mirrors what `builder.rs` step 4d
    /// does for a one-step group.
    fn drive_single(
        step: Box<dyn ErasedStep>,
        step_idx: StepIdx,
        contexts: &Arc<ChainContexts>,
        signal: &Arc<PipelineSignal>,
    ) {
        let drain_counters: Vec<Arc<StepDrainCounter>> =
            (0..contexts.inputs.len()).map(|_| StepDrainCounter::new(1)).collect();
        let group = DetachedDriverGroup::new(vec![(step_idx, step)]);
        run_detached_driver(
            group,
            contexts,
            &drain_counters,
            signal,
            None,
            &crate::liveness::LivenessCounter::new(1),
            None,
            0,
            &crate::runtime::wake::WakePlan::legacy(None),
            None,
        );
    }

    /// Items pushed onto a Detached step's input flow through to its output;
    /// the step finishes once its input is drained and closes the output edge.
    #[test]
    fn detached_step_flows_items_and_finishes() {
        let (det, contexts, producer_outputs_any, consumer) = build_one_detached(64);
        let producer_outputs =
            producer_outputs_any.downcast_ref::<OutputHandles<Single<u32>>>().unwrap();
        producer_outputs.push(10).unwrap();
        producer_outputs.push(20).unwrap();
        producer_outputs.push(30).unwrap();
        producer_outputs.mark_all_drained();

        let signal = PipelineSignal::new();
        drive_single(det, StepIdx(1), &contexts, &signal);

        let mut got = Vec::new();
        while let Some(v) = consumer.pop() {
            got.push(v);
        }
        assert_eq!(got, vec![10, 20, 30]);
        assert!(InputHandle::is_drained(&consumer), "output closed on Finished");
        assert!(!signal.is_done(), "clean completion, no error");
    }

    /// Zero items in (input drained from the start): the Detached step cleanly
    /// drains its output and returns.
    #[test]
    fn detached_step_zero_items_clean_drain() {
        let (det, contexts, producer_outputs_any, consumer) = build_one_detached(64);
        let producer_outputs =
            producer_outputs_any.downcast_ref::<OutputHandles<Single<u32>>>().unwrap();
        producer_outputs.mark_all_drained(); // no items

        let signal = PipelineSignal::new();
        drive_single(det, StepIdx(1), &contexts, &signal);

        assert!(consumer.pop().is_none(), "no items produced");
        assert!(InputHandle::is_drained(&consumer), "output closed on clean empty drain");
    }

    /// A two-sided Detached step — consumer of a bounded pool-fed input AND
    /// producer to a bounded pool-drained output, both tiny — completes without
    /// deadlock. A wall-clock watchdog aborts (fails the test) if it wedges.
    /// This is the sort merge's topology.
    #[test]
    fn detached_two_sided_no_deadlock() {
        const N: u32 = 5_000;

        // cap-2 input AND cap-4 output both force interleaved backpressure.
        let (det, contexts, producer_outputs_any, consumer) = build_one_detached(2);
        let signal = PipelineSignal::new();

        // Watchdog: a deadlock parks forever; abort so the test FAILS loudly.
        let done = Arc::new(std::sync::atomic::AtomicBool::new(false));
        {
            let done = Arc::clone(&done);
            std::thread::spawn(move || {
                for _ in 0..200 {
                    std::thread::sleep(Duration::from_millis(50));
                    if done.load(Ordering::SeqCst) {
                        return;
                    }
                }
                eprintln!("detached_two_sided_no_deadlock: WEDGED (deadlock)");
                std::process::abort();
            });
        }

        // Producer thread: push N items into the cap-2 input (backpressure),
        // then close it so the Detached step's input drains.
        let pushed = Arc::new(AtomicU32::new(0));
        let producer_handle = {
            let producer_outputs_any = Arc::clone(&producer_outputs_any);
            let pushed = Arc::clone(&pushed);
            std::thread::spawn(move || {
                let outputs =
                    producer_outputs_any.downcast_ref::<OutputHandles<Single<u32>>>().unwrap();
                let mut held: Option<u32> = None;
                let mut next = 0u32;
                loop {
                    if let Some(v) = held.take() {
                        match outputs.push(v) {
                            Ok(()) => {
                                pushed.fetch_add(1, Ordering::Relaxed);
                            }
                            Err(unpushed) => {
                                held = Some(unpushed.into_item());
                                std::thread::yield_now();
                            }
                        }
                        continue;
                    }
                    if next >= N {
                        break;
                    }
                    match outputs.push(next) {
                        Ok(()) => {
                            pushed.fetch_add(1, Ordering::Relaxed);
                            next += 1;
                        }
                        Err(unpushed) => {
                            // The value at `next` is now held for retry; advance
                            // `next` so the fresh-push branch doesn't re-emit it
                            // after `held` flushes (which would double-count).
                            held = Some(unpushed.into_item());
                            next += 1;
                            std::thread::yield_now();
                        }
                    }
                }
                outputs.mark_all_drained();
            })
        };

        // Consumer thread: pop everything the Detached step produces. Keep the
        // VALUES, not just a count — a count-only assertion passes for a step that
        // emits N copies of one item, or that duplicates the held value while
        // dropping a popped one, which is exactly the loss this test exists to
        // catch.
        let received = Arc::new(parking_lot::Mutex::new(Vec::<u32>::new()));
        let consumer_handle = {
            let received = Arc::clone(&received);
            std::thread::spawn(move || {
                loop {
                    if let Some(v) = consumer.pop() {
                        received.lock().push(v);
                    } else if InputHandle::is_drained(&consumer) {
                        break;
                    } else {
                        std::thread::yield_now();
                    }
                }
            })
        };

        // Drive the Detached step (as a one-step group) on this thread to
        // completion via the unified driver.
        drive_single(det, StepIdx(1), &contexts, &signal);

        producer_handle.join().unwrap();
        consumer_handle.join().unwrap();
        done.store(true, Ordering::SeqCst);

        assert_eq!(pushed.load(Ordering::Relaxed), N, "all items pushed");
        // Sorted multiset, not the sequence: the producer's backpressure branch
        // holds `next` and advances, so a retried item can arrive after a later
        // one. Order is not the invariant here; every distinct item arriving
        // exactly once is.
        let mut got = received.lock().clone();
        got.sort_unstable();
        assert_eq!(
            got,
            (0..N).collect::<Vec<u32>>(),
            "every distinct item must flow through the two-sided Detached step exactly once"
        );
    }

    /// A `Detached` step declaring a `Shared` group label. Used only to exercise
    /// `extract_detached_steps` grouping — never actually run.
    #[derive(Clone)]
    struct SharedDetached {
        label: &'static str,
    }
    impl Step for SharedDetached {
        type Input = u32;
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "SharedDetached",
                kind: StepKind::Detached,
                sticky: false,
                output_queues: vec![QueueSpec::CountBounded { capacity: 4 }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn detached_group(&self) -> crate::step::DetachedGroup {
            crate::step::DetachedGroup::Shared(self.label)
        }
        fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            Ok(StepOutcome::Finished)
        }
    }

    /// A `Parallel` step — on a driver only as a clone hosted there alone.
    #[derive(Clone)]
    struct ParallelStub;
    impl Step for ParallelStub {
        type Input = u32;
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "ParallelStub",
                kind: StepKind::Parallel,
                sticky: false,
                output_queues: vec![QueueSpec::CountBounded { capacity: 4 }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            Ok(StepOutcome::NoProgress)
        }
        fn new_worker_copy(&self) -> Self {
            self.clone()
        }
    }

    /// The wake plan's driver map agrees with the driver groups
    /// `extract_detached_steps` builds.
    #[test]
    fn driver_index_of_agrees_with_extracted_groups() {
        use crate::runtime::wake::DriverIdx;
        let mut steps: Vec<Box<dyn ErasedStep>> = vec![
            Box::new(TypedStep::new(SrcStub { capacity: 4 })),
            Box::new(TypedStep::new(SharedDetached { label: "coord" })),
            Box::new(TypedStep::new(PassThroughDetached { held: None })),
            Box::new(TypedStep::new(SharedDetached { label: "coord" })),
            Box::new(TypedStep::new(SharedDetached { label: "io" })),
        ];
        let driver_of = driver_index_of(&steps);
        assert_eq!(
            driver_of,
            vec![
                None,
                Some(DriverIdx(0)),
                Some(DriverIdx(1)),
                Some(DriverIdx(0)),
                Some(DriverIdx(2))
            ]
        );
        let groups = extract_detached_steps(&mut steps);
        for (d, g) in groups.iter().enumerate() {
            for idx in g.step_indices() {
                assert_eq!(driver_of[idx.0], Some(DriverIdx(d)), "step {}", idx.0);
            }
        }
    }

    /// `extract_detached_steps` collects every `Shared(label)` onto one group,
    /// keeps each `PerStep` (default) step as its own singleton, preserves chain
    /// order within and across groups, and leaves non-detached steps in place.
    #[test]
    fn extract_groups_shared_together_and_perstep_alone() {
        use crate::step::DetachedGroup;
        // idx0 non-detached; idx1/idx3 Shared("coord"); idx2 PerStep; idx4 Shared("io").
        let mut steps: Vec<Box<dyn ErasedStep>> = vec![
            Box::new(TypedStep::new(SrcStub { capacity: 4 })),
            Box::new(TypedStep::new(SharedDetached { label: "coord" })),
            Box::new(TypedStep::new(PassThroughDetached { held: None })),
            Box::new(TypedStep::new(SharedDetached { label: "coord" })),
            Box::new(TypedStep::new(SharedDetached { label: "io" })),
        ];
        let groups = extract_detached_steps(&mut steps);

        assert_eq!(groups.len(), 3, "coord{{1,3}}, perstep{{2}}, io{{4}}");
        // Order is first-appearance in chain order.
        assert_eq!(groups[0].label(), DetachedGroup::Shared("coord"));
        assert_eq!(
            groups[0].steps.iter().map(|(i, _)| i.0).collect::<Vec<_>>(),
            vec![1, 3],
            "shared group keeps both steps in chain order"
        );
        assert_eq!(groups[0].primary_step(), StepIdx(1));
        assert_eq!(groups[1].label(), DetachedGroup::PerStep);
        assert_eq!(groups[1].steps.iter().map(|(i, _)| i.0).collect::<Vec<_>>(), vec![2]);
        assert_eq!(groups[2].label(), DetachedGroup::Shared("io"));

        // Non-detached step survives; detached slots became placeholders that
        // still report `Detached` (so `build_worker_storage` Skips them on pool).
        assert_eq!(steps[0].kind(), StepKind::Exclusive);
        assert_eq!(steps[1].kind(), StepKind::Detached);
        assert_eq!(steps[2].kind(), StepKind::Detached);

        // Each placeholder reports the group its real step declared, not a blanket
        // `PerStep`. `extract_detached_steps` is `pub` and hands this slice back,
        // so a caller regrouping from it would otherwise split the `coord` group
        // across separate driver threads.
        assert_eq!(
            steps[1].detached_group(),
            DetachedGroup::Shared("coord"),
            "placeholder must not downgrade a Shared group to PerStep"
        );
        assert_eq!(steps[2].detached_group(), DetachedGroup::PerStep);
    }

    #[test]
    #[should_panic(expected = "must be non-empty")]
    fn driver_group_rejects_empty() {
        // An empty group would panic later in `primary_step`/`label` (`steps[0]`);
        // the constructor rejects it up front.
        let _ = DetachedDriverGroup::new(vec![]);
    }

    #[test]
    #[should_panic(expected = "only contain Detached steps")]
    fn driver_group_rejects_non_detached_step() {
        // `SrcStub` is Exclusive, not Detached — running it off-pool on a driver
        // thread is a bug, so the constructor rejects the group.
        let step: Box<dyn ErasedStep> = Box::new(TypedStep::new(SrcStub { capacity: 1 }));
        let _ = DetachedDriverGroup::new(vec![(StepIdx(0), step)]);
    }

    /// `n` drain counters, each initialized to `init`.
    fn counters(n: usize, init: usize) -> Vec<Arc<StepDrainCounter>> {
        (0..n).map(|_| StepDrainCounter::new(init)).collect()
    }

    /// `build_driver_storage` makes the group's steps `Owned` and every other
    /// slot `Skip`, at the correct global indices.
    #[test]
    fn build_driver_storage_owns_group_skips_rest() {
        let det: Box<dyn ErasedStep> = Box::new(TypedStep::new(PassThroughDetached { held: None }));
        let row = build_driver_storage(vec![(StepIdx(2), det)], &counters(5, 1));
        assert_eq!(row.len(), 5);
        assert!(matches!(row[2], WorkerStepEntry::Owned { .. }), "group step is Owned");
        for i in [0usize, 1, 3, 4] {
            assert!(matches!(row[i], WorkerStepEntry::Skip), "non-group slot {i} is Skip");
        }
    }

    /// G3: a member whose drain counter counts pool clones (here a `Parallel`
    /// step at 2) is rejected: one driver instance could never take it to 0.
    #[test]
    #[should_panic(expected = "has drain counter 2, not 1")]
    fn build_driver_storage_rejects_parallel_group_step() {
        let par: Box<dyn ErasedStep> = Box::new(TypedStep::new(ParallelStub));
        let _ = build_driver_storage(vec![(StepIdx(0), par)], &counters(2, 2));
    }

    /// A `Parallel` member whose counter is 1 — a clone `plan_parallel_hosts`
    /// hosted on this driver alone — is accepted; the same step at 2 is the
    /// panic above.
    #[test]
    fn build_driver_storage_accepts_a_hosted_parallel_member() {
        let det: Box<dyn ErasedStep> = Box::new(TypedStep::new(PassThroughDetached { held: None }));
        let par: Box<dyn ErasedStep> = Box::new(TypedStep::new(ParallelStub));
        let row = build_driver_storage(vec![(StepIdx(0), det), (StepIdx(1), par)], &counters(2, 1));
        assert!(matches!(row[1], WorkerStepEntry::Owned { .. }));
    }

    /// A hosted clone joins its driver's group at its chain position and leaves
    /// a `Parallel`-reporting placeholder in `steps`; a Parallel step with worker
    /// hosts stays where it is.
    #[test]
    fn extract_hosted_parallel_steps_inserts_in_chain_order() {
        use crate::runtime::wake::DriverIdx;
        // One shared driver over steps 1 and 3, so the hosted step 2 must land
        // BETWEEN them (an append would put it last).
        let mut steps: Vec<Box<dyn ErasedStep>> = vec![
            Box::new(TypedStep::new(SrcStub { capacity: 1 })),
            Box::new(TypedStep::new(SharedDetached { label: "g" })),
            Box::new(TypedStep::new(ParallelStub)),
            Box::new(TypedStep::new(SharedDetached { label: "g" })),
            Box::new(TypedStep::new(ParallelStub)),
        ];
        let mut groups = extract_detached_steps(&mut steps);
        assert_eq!(groups.len(), 1, "one Shared label, one driver");
        let hosts = vec![
            ParallelHosts::default(),
            ParallelHosts::default(),
            ParallelHosts { workers: vec![], driver: Some(DriverIdx(0)) },
            ParallelHosts::default(),
            ParallelHosts { workers: vec![0], driver: None },
        ];
        extract_hosted_parallel_steps(&mut steps, &hosts, &mut groups);
        let members: Vec<(usize, StepKind)> =
            groups[0].steps.iter().map(|(i, s)| (i.0, s.kind())).collect();
        assert_eq!(
            members,
            vec![(1, StepKind::Detached), (2, StepKind::Parallel), (3, StepKind::Detached)]
        );
        assert_eq!(groups[0].primary_step(), StepIdx(1), "the first member stays Detached");
        assert_eq!(steps[2].kind(), StepKind::Parallel, "placeholder reports Parallel");
        assert_eq!(steps[2].name(), "ParallelStub", "placeholder keeps the name");
        assert_eq!(steps[2].pool_placement(), steps[4].pool_placement(), "and the placement");
        // The pool-placed step is still the real instance (cloning it works).
        let _ = steps[4].clone_boxed();
    }

    /// A hosted step may not become a group's first member: the first member
    /// names the driver thread and keys its stats, and must stay Detached.
    #[test]
    #[should_panic(expected = "must follow its Detached producer")]
    fn insert_hosted_rejects_a_step_before_the_groups_detached_steps() {
        let det: Box<dyn ErasedStep> = Box::new(TypedStep::new(PassThroughDetached { held: None }));
        let mut group = DetachedDriverGroup::new(vec![(StepIdx(3), det)]);
        group.insert_hosted(StepIdx(1), Box::new(TypedStep::new(ParallelStub)));
    }

    /// Only a `Parallel` step's clone is hosted.
    #[test]
    #[should_panic(expected = "only a Parallel step's clone is hosted")]
    fn insert_hosted_rejects_a_non_parallel_step() {
        let det: Box<dyn ErasedStep> = Box::new(TypedStep::new(PassThroughDetached { held: None }));
        let mut group = DetachedDriverGroup::new(vec![(StepIdx(0), det)]);
        group.insert_hosted(StepIdx(1), Box::new(TypedStep::new(SrcStub { capacity: 1 })));
    }

    /// G3: two group steps at the same index is a dual registration — rejected.
    #[test]
    #[should_panic(expected = "registered twice")]
    fn build_driver_storage_rejects_dual_registration() {
        let a: Box<dyn ErasedStep> = Box::new(TypedStep::new(PassThroughDetached { held: None }));
        let b: Box<dyn ErasedStep> = Box::new(TypedStep::new(PassThroughDetached { held: None }));
        let _ = build_driver_storage(vec![(StepIdx(1), a), (StepIdx(1), b)], &counters(3, 1));
    }

    /// G3: a group step index past the step count is the documented out-of-range
    /// panic, matching the guard the sibling `build_worker_storage` applies —
    /// without it the bare `row[idx.0]` index would panic with an opaque message.
    #[test]
    #[should_panic(expected = "out of range")]
    fn build_driver_storage_rejects_out_of_range_index() {
        let det: Box<dyn ErasedStep> = Box::new(TypedStep::new(PassThroughDetached { held: None }));
        let _ = build_driver_storage(vec![(StepIdx(5), det)], &counters(3, 1));
    }
}

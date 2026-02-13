//! Scheduler trait and common types for pipeline step selection.
//!
//! This module provides pluggable scheduling strategies for the unified pipeline.
//! Different schedulers make different tradeoffs between simplicity, contention,
//! and adaptability.

use clap::ValueEnum;

use super::base::{ActiveSteps, PipelineStep};

#[doc(hidden)]
pub mod backpressure_proportional;
#[doc(hidden)]
pub mod balanced_chase;
#[doc(hidden)]
pub mod balanced_chase_drain;
#[doc(hidden)]
pub mod chase_bottleneck;
#[doc(hidden)]
pub mod epsilon_greedy;
#[doc(hidden)]
pub mod fixed_priority;
#[doc(hidden)]
pub mod hybrid_adaptive;
#[doc(hidden)]
pub mod learned_affinity;
#[doc(hidden)]
pub mod optimized_chase;
#[doc(hidden)]
pub mod sticky_work_stealing;
#[doc(hidden)]
pub mod thompson_sampling;
#[doc(hidden)]
pub mod thompson_with_priors;
#[doc(hidden)]
pub mod two_phase;
#[doc(hidden)]
pub mod ucb;

#[doc(hidden)]
pub use backpressure_proportional::BackpressureProportionalScheduler;
#[doc(hidden)]
pub use balanced_chase::BalancedChaseScheduler;
#[doc(hidden)]
pub use balanced_chase_drain::BalancedChaseDrainScheduler;
#[doc(hidden)]
pub use chase_bottleneck::ChaseBottleneckScheduler;
#[doc(hidden)]
pub use epsilon_greedy::EpsilonGreedyScheduler;
#[doc(hidden)]
pub use fixed_priority::FixedPriorityScheduler;
#[doc(hidden)]
pub use hybrid_adaptive::HybridAdaptiveScheduler;
#[doc(hidden)]
pub use learned_affinity::LearnedAffinityScheduler;
#[doc(hidden)]
pub use optimized_chase::OptimizedChaseScheduler;
#[doc(hidden)]
pub use sticky_work_stealing::StickyWorkStealingScheduler;
#[doc(hidden)]
pub use thompson_sampling::ThompsonSamplingScheduler;
#[doc(hidden)]
pub use thompson_with_priors::ThompsonWithPriorsScheduler;
#[doc(hidden)]
pub use two_phase::TwoPhaseScheduler;
#[doc(hidden)]
pub use ucb::UCBScheduler;

/// Scheduler strategy for pipeline execution.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, ValueEnum)]
pub enum SchedulerStrategy {
    /// Fixed priority scheduling based on thread role.
    ///
    /// Thread 0 prioritizes reading, Thread N-1 prioritizes writing,
    /// middle threads rotate among parallel steps. Includes backpressure
    /// override when output queue fills.
    #[value(name = "fixed-priority")]
    FixedPriority,

    /// Chase-bottleneck scheduling with dynamic adaptation.
    ///
    /// Threads follow work: move downstream when output blocked,
    /// move upstream when input empty, stay sticky on success.
    /// Automatically rebalances as pipeline stages progress.
    /// Shows ~10% improvement at medium thread counts (4 threads).
    #[value(name = "chase-bottleneck")]
    ChaseBottleneck,

    /// Thompson Sampling with Beta distributions.
    ///
    /// Uses Bayesian inference to balance exploration and exploitation.
    /// Each step maintains a Beta(α, β) distribution updated on success/failure.
    #[value(name = "thompson-sampling")]
    ThompsonSampling,

    /// Upper Confidence Bound algorithm.
    ///
    /// Prioritizes steps with high success rate plus uncertainty bonus.
    /// Naturally explores under-tried steps while exploiting successful ones.
    #[value(name = "ucb")]
    UCB,

    /// Epsilon-Greedy exploration/exploitation.
    ///
    /// With probability ε (10%), explores randomly.
    /// Otherwise exploits the step with highest observed success rate.
    #[value(name = "epsilon-greedy")]
    EpsilonGreedy,

    /// Thompson Sampling with thread-specific priors.
    ///
    /// Like Thompson Sampling, but initializes with biased priors based on
    /// thread role (reader→Read, writer→Write, etc.).
    #[value(name = "thompson-with-priors")]
    ThompsonWithPriors,

    /// Hybrid: switches between fixed-priority and chase-bottleneck.
    ///
    /// Starts with fixed-priority for efficiency. After consecutive failures,
    /// switches to chase-bottleneck for adaptability. Returns when stable.
    #[value(name = "hybrid-adaptive")]
    HybridAdaptive,

    /// Backpressure-proportional with EMA weights.
    ///
    /// Dynamically adjusts step weights based on queue depths.
    /// Downstream steps get higher priority when output backs up.
    #[value(name = "backpressure-proportional")]
    BackpressureProportional,

    /// Two-phase: startup/steady-state/drain optimization.
    ///
    /// Uses chase-bottleneck during startup (fill pipeline) and drain (empty pipeline).
    /// Uses fixed-priority during steady-state for efficiency.
    #[value(name = "two-phase")]
    TwoPhase,

    /// Sticky work-stealing with home steps.
    ///
    /// Each thread has a "home" step based on role. When home step has no work,
    /// steals from adjacent steps first, then any step. Periodically returns home.
    #[value(name = "sticky-work-stealing")]
    StickyWorkStealing,

    /// Learned affinity with decaying exploration.
    ///
    /// Tracks success rates per step and builds learned priority order.
    /// Exploration rate decays over time to converge on optimal strategy.
    #[value(name = "learned-affinity")]
    LearnedAffinity,

    /// Optimized chase-bottleneck with profiling-based improvements.
    ///
    /// Enhanced chase-bottleneck with:
    /// - Compress-biased prioritization (Compress is the bottleneck)
    /// - Exclusive step avoidance (non-specialists deprioritize exclusive steps)
    /// - Bottleneck stickiness (stay on Compress/Serialize longer)
    /// - Contention backoff (avoid exclusive steps after contention)
    #[value(name = "optimized-chase")]
    OptimizedChase,

    /// Balanced chase scheduler focused on even work distribution.
    ///
    /// Key insight: exclusive specialists (T0=Read, T7=Write) should help
    /// with bottleneck steps instead of staying sticky. After completing
    /// exclusive work, immediately pivot to Compress/Serialize.
    #[value(name = "balanced-chase")]
    BalancedChase,

    /// Balanced chase with drain mode for output backpressure (default).
    ///
    /// Like balanced-chase, but when Serialize fails due to Q6 being full,
    /// enters drain mode: prioritize Compress until backpressure clears.
    /// This is backpressure-driven rather than using a fixed iteration count.
    #[default]
    #[value(name = "balanced-chase-drain")]
    BalancedChaseDrain,
}

/// Direction for step traversal in chase-bottleneck scheduler.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Direction {
    /// Move towards output (drain downstream)
    Forward,
    /// Move towards input (fill upstream)
    Backward,
}

/// Backpressure state from queue depths and memory pressure.
#[derive(Debug, Clone, Copy)]
#[allow(clippy::struct_excessive_bools)]
pub struct BackpressureState {
    /// True when output queue exceeds high-water mark.
    pub output_high: bool,
    /// True when input queue is below low-water mark.
    pub input_low: bool,
    /// True when reading is complete.
    pub read_done: bool,
    /// True when memory tracker is at or above limit.
    /// When set, schedulers should prioritize steps that drain memory (e.g., Process).
    pub memory_high: bool,
    /// True when memory has drained below the low-water mark (50% of limit).
    /// Used with hysteresis: enter drain mode when `memory_high`, exit when `memory_drained`.
    pub memory_drained: bool,
}

impl Default for BackpressureState {
    fn default() -> Self {
        Self {
            output_high: false,
            input_low: false,
            read_done: false,
            memory_high: false,
            // Default to true: when there's no memory tracking/pressure,
            // we consider memory "drained" (fine) to allow drain mode to exit.
            memory_drained: true,
        }
    }
}

/// Trait for pipeline schedulers.
///
/// Schedulers determine which pipeline step a thread should attempt next.
/// Different strategies make different tradeoffs between simplicity,
/// contention, and adaptability.
pub trait Scheduler: Send {
    /// Get the priority order of steps to try.
    ///
    /// Returns a slice of steps in priority order. The worker loop
    /// tries each step in order until one succeeds.
    fn get_priorities(&mut self, backpressure: BackpressureState) -> &[PipelineStep];

    /// Record the outcome of a step attempt.
    ///
    /// Used by adaptive schedulers to learn from successes and failures.
    /// The `was_contention` flag indicates lock contention vs queue empty/full.
    fn record_outcome(&mut self, step: PipelineStep, success: bool, was_contention: bool);

    /// Get the thread ID for this scheduler instance.
    fn thread_id(&self) -> usize;

    /// Get the number of threads in the pipeline.
    fn num_threads(&self) -> usize;

    /// Get the active steps for this pipeline.
    fn active_steps(&self) -> &ActiveSteps;

    /// Returns the exclusive step this thread owns, if any.
    ///
    /// Maps thread IDs to active exclusive steps from both ends:
    /// - T0 → first active exclusive step
    /// - T\[N-1\] → last active exclusive step
    /// - T1 → second exclusive step (from front)
    /// - T\[N-2\] → second-to-last exclusive step (from back)
    /// - Others: None (parallel-only workers)
    ///
    /// This preserves backward compatibility: with all 4 exclusive steps
    /// (`Read`, `FindBoundaries`, `Group`, `Write`) and 8 threads, the mapping is
    /// T0→Read, T1→FindBoundaries, T6→Group, T7→Write — same as the
    /// original hardcoded mapping.
    ///
    /// For `num_threads < 4`: Returns `None` for all threads (relaxed ownership).
    fn exclusive_step_owned(&self) -> Option<PipelineStep> {
        let num_threads = self.num_threads();
        let thread_id = self.thread_id();

        if num_threads < 4 {
            return None;
        }

        let exclusive = self.active_steps().exclusive_steps();
        if exclusive.is_empty() {
            return None;
        }

        // T0 → first exclusive step
        if thread_id == 0 {
            return Some(exclusive[0]);
        }

        // TN-1 → last exclusive step (if more than one)
        if thread_id == num_threads - 1 && exclusive.len() > 1 {
            return Some(*exclusive.last().unwrap());
        }

        // Interior exclusive steps: assign from both ends
        let num_interior = exclusive.len().saturating_sub(2);
        if num_interior == 0 {
            return None;
        }

        // Front half: T1→exclusive[1], T2→exclusive[2], ...
        let front_count = num_interior.div_ceil(2);
        if thread_id >= 1 && thread_id <= front_count {
            return Some(exclusive[thread_id]);
        }

        // Back half: TN-2→exclusive[len-2], TN-3→exclusive[len-3], ...
        let back_count = num_interior - front_count;
        let back_thread_offset = num_threads - 1 - thread_id;
        if back_thread_offset >= 1 && back_thread_offset <= back_count {
            return Some(exclusive[exclusive.len() - 1 - back_thread_offset]);
        }

        None
    }

    /// Returns true if this thread should attempt the given step.
    ///
    /// For exclusive steps: only the owner should attempt.
    /// For parallel steps: all threads can attempt.
    ///
    /// Note: This is the default behavior. In drain phase (`read_done && group_done`),
    /// the worker loop may override this to allow all threads to help with exclusive steps.
    fn should_attempt_step(&self, step: PipelineStep) -> bool {
        if !step.is_exclusive() {
            return true;
        }

        let num_threads = self.num_threads();
        if num_threads < 4 {
            return true;
        }

        match self.exclusive_step_owned() {
            Some(owned) => owned == step,
            None => false,
        }
    }

    /// Returns true if this thread should attempt the given step, considering drain mode.
    ///
    /// In drain mode (`read_done && group_done`), all threads can help with exclusive steps
    /// to avoid starvation. Otherwise, only the owner can attempt exclusive steps.
    fn should_attempt_step_with_drain(&self, step: PipelineStep, drain_mode: bool) -> bool {
        if !step.is_exclusive() {
            return true;
        }

        // In drain mode, all threads can help with exclusive steps
        if drain_mode {
            return true;
        }

        let num_threads = self.num_threads();
        if num_threads < 4 {
            return true;
        }

        match self.exclusive_step_owned() {
            Some(owned) => owned == step,
            None => false,
        }
    }
}

/// Create a scheduler based on strategy.
#[must_use]
pub fn create_scheduler(
    strategy: SchedulerStrategy,
    thread_id: usize,
    num_threads: usize,
    active_steps: ActiveSteps,
) -> Box<dyn Scheduler> {
    match strategy {
        SchedulerStrategy::FixedPriority => {
            Box::new(FixedPriorityScheduler::new(thread_id, num_threads, active_steps))
        }
        SchedulerStrategy::ChaseBottleneck => {
            Box::new(ChaseBottleneckScheduler::new(thread_id, num_threads, active_steps))
        }
        SchedulerStrategy::ThompsonSampling => {
            Box::new(ThompsonSamplingScheduler::new(thread_id, num_threads, active_steps))
        }
        SchedulerStrategy::UCB => Box::new(UCBScheduler::new(thread_id, num_threads, active_steps)),
        SchedulerStrategy::EpsilonGreedy => {
            Box::new(EpsilonGreedyScheduler::new(thread_id, num_threads, active_steps))
        }
        SchedulerStrategy::ThompsonWithPriors => {
            Box::new(ThompsonWithPriorsScheduler::new(thread_id, num_threads, active_steps))
        }
        SchedulerStrategy::HybridAdaptive => {
            Box::new(HybridAdaptiveScheduler::new(thread_id, num_threads, active_steps))
        }
        SchedulerStrategy::BackpressureProportional => {
            Box::new(BackpressureProportionalScheduler::new(thread_id, num_threads, active_steps))
        }
        SchedulerStrategy::TwoPhase => {
            Box::new(TwoPhaseScheduler::new(thread_id, num_threads, active_steps))
        }
        SchedulerStrategy::StickyWorkStealing => {
            Box::new(StickyWorkStealingScheduler::new(thread_id, num_threads, active_steps))
        }
        SchedulerStrategy::LearnedAffinity => {
            Box::new(LearnedAffinityScheduler::new(thread_id, num_threads, active_steps))
        }
        SchedulerStrategy::OptimizedChase => {
            Box::new(OptimizedChaseScheduler::new(thread_id, num_threads, active_steps))
        }
        SchedulerStrategy::BalancedChase => {
            Box::new(BalancedChaseScheduler::new(thread_id, num_threads, active_steps))
        }
        SchedulerStrategy::BalancedChaseDrain => {
            Box::new(BalancedChaseDrainScheduler::new(thread_id, num_threads, active_steps))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    fn all() -> ActiveSteps {
        ActiveSteps::all()
    }

    #[test]
    fn test_scheduler_strategy_default() {
        assert_eq!(SchedulerStrategy::default(), SchedulerStrategy::BalancedChaseDrain);
    }

    #[rstest]
    #[case::fixed_priority(SchedulerStrategy::FixedPriority)]
    #[case::chase_bottleneck(SchedulerStrategy::ChaseBottleneck)]
    #[case::thompson_sampling(SchedulerStrategy::ThompsonSampling)]
    #[case::ucb(SchedulerStrategy::UCB)]
    #[case::epsilon_greedy(SchedulerStrategy::EpsilonGreedy)]
    #[case::thompson_with_priors(SchedulerStrategy::ThompsonWithPriors)]
    #[case::hybrid_adaptive(SchedulerStrategy::HybridAdaptive)]
    #[case::backpressure_proportional(SchedulerStrategy::BackpressureProportional)]
    #[case::two_phase(SchedulerStrategy::TwoPhase)]
    #[case::sticky_work_stealing(SchedulerStrategy::StickyWorkStealing)]
    #[case::learned_affinity(SchedulerStrategy::LearnedAffinity)]
    #[case::optimized_chase(SchedulerStrategy::OptimizedChase)]
    #[case::balanced_chase(SchedulerStrategy::BalancedChase)]
    #[case::balanced_chase_drain(SchedulerStrategy::BalancedChaseDrain)]
    fn test_create_scheduler(#[case] strategy: SchedulerStrategy) {
        let scheduler = create_scheduler(strategy, 0, 8, all());
        assert_eq!(scheduler.thread_id(), 0);
    }

    #[test]
    fn test_exclusive_step_ownership_small_thread_counts() {
        // 2 threads: relaxed ownership (all can attempt all)
        let s0 = create_scheduler(SchedulerStrategy::BalancedChase, 0, 2, all());
        let s1 = create_scheduler(SchedulerStrategy::BalancedChase, 1, 2, all());

        assert!(s0.exclusive_step_owned().is_none());
        assert!(s1.exclusive_step_owned().is_none());
        assert!(s0.should_attempt_step(PipelineStep::Group));
        assert!(s1.should_attempt_step(PipelineStep::Group));

        // 3 threads: also relaxed
        let s0 = create_scheduler(SchedulerStrategy::BalancedChase, 0, 3, all());
        let s1 = create_scheduler(SchedulerStrategy::BalancedChase, 1, 3, all());
        let s2 = create_scheduler(SchedulerStrategy::BalancedChase, 2, 3, all());

        assert!(s0.should_attempt_step(PipelineStep::FindBoundaries));
        assert!(s1.should_attempt_step(PipelineStep::Group));
        assert!(s2.should_attempt_step(PipelineStep::Write));
    }

    #[test]
    fn test_exclusive_step_ownership_eight_threads() {
        // 8 threads: strict ownership
        for thread_id in 0..8 {
            let scheduler = create_scheduler(SchedulerStrategy::BalancedChase, thread_id, 8, all());

            let expected_ownership = match thread_id {
                0 => Some(PipelineStep::Read),
                1 => Some(PipelineStep::FindBoundaries),
                6 => Some(PipelineStep::Group), // N-2 = 6
                7 => Some(PipelineStep::Write), // N-1 = 7
                _ => None,
            };

            assert_eq!(
                scheduler.exclusive_step_owned(),
                expected_ownership,
                "Thread {thread_id} ownership mismatch"
            );
        }
    }

    #[test]
    fn test_should_attempt_step_parallel_always_allowed() {
        let scheduler = create_scheduler(SchedulerStrategy::BalancedChase, 3, 8, all());

        // Parallel steps: all threads can attempt
        assert!(scheduler.should_attempt_step(PipelineStep::Decompress));
        assert!(scheduler.should_attempt_step(PipelineStep::Decode));
        assert!(scheduler.should_attempt_step(PipelineStep::Process));
        assert!(scheduler.should_attempt_step(PipelineStep::Serialize));
        assert!(scheduler.should_attempt_step(PipelineStep::Compress));
    }

    #[test]
    fn test_should_attempt_step_exclusive_only_owner() {
        // T3 (middle worker) should not attempt any exclusive steps
        let t3 = create_scheduler(SchedulerStrategy::BalancedChase, 3, 8, all());
        assert!(!t3.should_attempt_step(PipelineStep::Read));
        assert!(!t3.should_attempt_step(PipelineStep::FindBoundaries));
        assert!(!t3.should_attempt_step(PipelineStep::Group));
        assert!(!t3.should_attempt_step(PipelineStep::Write));

        // T6 (Group owner) should only attempt Group among exclusive steps
        let t6 = create_scheduler(SchedulerStrategy::BalancedChase, 6, 8, all());
        assert!(!t6.should_attempt_step(PipelineStep::Read));
        assert!(!t6.should_attempt_step(PipelineStep::FindBoundaries));
        assert!(t6.should_attempt_step(PipelineStep::Group)); // Owner!
        assert!(!t6.should_attempt_step(PipelineStep::Write));
    }

    #[test]
    fn test_exclusive_ownership_all_strategies() {
        let strategies = [
            SchedulerStrategy::FixedPriority,
            SchedulerStrategy::ChaseBottleneck,
            SchedulerStrategy::BalancedChase,
            SchedulerStrategy::OptimizedChase,
            SchedulerStrategy::ThompsonSampling,
            SchedulerStrategy::UCB,
            SchedulerStrategy::EpsilonGreedy,
            SchedulerStrategy::ThompsonWithPriors,
            SchedulerStrategy::HybridAdaptive,
            SchedulerStrategy::BackpressureProportional,
            SchedulerStrategy::TwoPhase,
            SchedulerStrategy::StickyWorkStealing,
            SchedulerStrategy::LearnedAffinity,
            SchedulerStrategy::BalancedChaseDrain,
        ];

        for strategy in strategies {
            let t6 = create_scheduler(strategy, 6, 8, all());
            assert_eq!(
                t6.exclusive_step_owned(),
                Some(PipelineStep::Group),
                "{strategy:?} T6 should own Group"
            );

            let t3 = create_scheduler(strategy, 3, 8, all());
            assert_eq!(t3.exclusive_step_owned(), None, "{strategy:?} T3 should own nothing");
        }
    }

    #[test]
    fn test_four_thread_edge_case() {
        // 4 threads: minimal strict ownership
        // T0 = Read, T1 = FindBoundaries, T2 = Group (N-2), T3 = Write (N-1)
        let t0 = create_scheduler(SchedulerStrategy::BalancedChase, 0, 4, all());
        let t1 = create_scheduler(SchedulerStrategy::BalancedChase, 1, 4, all());
        let t2 = create_scheduler(SchedulerStrategy::BalancedChase, 2, 4, all());
        let t3 = create_scheduler(SchedulerStrategy::BalancedChase, 3, 4, all());

        assert_eq!(t0.exclusive_step_owned(), Some(PipelineStep::Read));
        assert_eq!(t1.exclusive_step_owned(), Some(PipelineStep::FindBoundaries));
        assert_eq!(t2.exclusive_step_owned(), Some(PipelineStep::Group));
        assert_eq!(t3.exclusive_step_owned(), Some(PipelineStep::Write));

        // All can do parallel steps
        assert!(t0.should_attempt_step(PipelineStep::Process));
        assert!(t1.should_attempt_step(PipelineStep::Process));
        assert!(t2.should_attempt_step(PipelineStep::Process));
        assert!(t3.should_attempt_step(PipelineStep::Process));
    }
}

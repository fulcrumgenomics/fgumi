//! Scheduler strategy enum shared by the scheduler dispatch and the CLI.
//!
//! Relocated out of `unified_pipeline::scheduler` (C5/R6a) to a neutral root
//! module rather than `commands::common`: it is consumed both by the
//! scheduler dispatch (`unified_pipeline::scheduler::create_scheduler`, a
//! pipeline-layer concern) and by the CLI's `SchedulerOptions`
//! (`commands::common`, the CLI layer). Homing it under `commands` would have
//! made the pipeline layer depend upward on the CLI layer through
//! `unified_pipeline::scheduler`'s shim; this module is downward-reachable
//! from both. It also survives C6, since the enum outlives
//! `unified_pipeline`.

/// Scheduler strategy for pipeline execution.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, clap::ValueEnum)]
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

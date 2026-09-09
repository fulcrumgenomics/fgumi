//! Scheduler strategy enum for the legacy `--scheduler` CLI flag.
//!
//! The typed-step chain engine does not use a pluggable scheduler strategy, so
//! this enum no longer drives any dispatch: it backs the hidden `--scheduler`
//! flag on `SchedulerOptions` (`commands::common`), and
//! `warn_unwired_pipeline_flags` warns when it is set to a non-default value so
//! a developer sees the flag is inert rather than silently ignored. It was
//! relocated to this neutral root module in C5/R6a (out of the legacy
//! multi-thread engine's `scheduler` module, removed in R6/C6).

/// Scheduler strategy that once selected a dispatch policy for the removed
/// legacy multi-thread engine.
///
/// The typed-step chain engine ignores this setting (see the module docs), so
/// none of the variants below drive any current dispatch. Their descriptions
/// are retained only as a historical record of what each `--scheduler` value
/// meant to the legacy engine.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, clap::ValueEnum)]
pub enum SchedulerStrategy {
    /// Historical: fixed-priority scheduling based on thread role.
    ///
    /// Thread 0 prioritized reading, thread N-1 prioritized writing,
    /// middle threads rotated among parallel steps, with a backpressure
    /// override when the output queue filled.
    #[value(name = "fixed-priority")]
    FixedPriority,

    /// Historical: chase-bottleneck scheduling with dynamic adaptation.
    ///
    /// Threads followed work: downstream when output blocked, upstream when
    /// input empty, sticky on success, rebalancing as pipeline stages
    /// progressed (~10% improvement at medium thread counts).
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

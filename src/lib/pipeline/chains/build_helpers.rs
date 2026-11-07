//! Shared helpers for per-command chain builders.

use std::sync::Arc;

use anyhow::Result;

use crate::commands::common::{QueueMemoryOptions, SchedulerOptions};
use crate::pipeline::core::builder::{InstrumentationLevel, Pipeline, PipelineConfig};
use crate::pipeline::core::runtime::stats::{PipelineStats, TimingResolution};

/// Truthiness of a boolean-ish environment variable: `true` iff `name` is set
/// and its value is neither `"0"` nor `"false"`. Shared by the
/// `FGUMI_PIPELINE_STATS` diagnostic checks here and in
/// [`ChainBuilder::build`](super::builder::ChainBuilder::build) so the two
/// cannot drift if the parsing rules ever change.
pub(crate) fn env_flag_enabled(name: &str) -> bool {
    flag_os_value_enabled(std::env::var_os(name).as_deref())
}

/// Pure predicate behind [`env_flag_enabled`]: `true` iff `value` is present
/// and is neither `"0"` nor `"false"` (so any other non-empty — or empty —
/// value counts as enabled). Split out from the env read so the parsing rule
/// is unit-testable without mutating the process environment (mutating it would
/// require `unsafe` under edition 2024, which this crate forbids).
fn flag_os_value_enabled(value: Option<&std::ffi::OsStr>) -> bool {
    value.is_some_and(|v| v != "0" && v != "false")
}

/// Build the [`PipelineConfig`] for a chain run and optionally attach
/// [`PipelineStats`]. Centralizes the threads + deadlock-timeout +
/// queue-memory + stats wiring that every per-command builder needs
/// identically.
///
/// Returns the constructed config and the `Arc<PipelineStats>` if stats
/// were attached. Stats are attached (`Some`) when any of:
/// - `--pipeline-stats` is set (`scheduler.collect_stats()`), or
/// - `--deadlock-timeout > 0`, since the deadlock monitor reads stats
///   counters to detect stalls, or
/// - `FGUMI_PIPELINE_STATS` is set (the standalone diagnostic knob), or
/// - tracing is on (`config.instrumentation.is_on()`), whose end-of-run
///   edge table and wall clock render through the stats snapshot.
///
/// Stats attached only for the deadlock monitor time dispatches on the coarse
/// clock (see [`stats_timing`]); every other reason times them precisely.
///
/// The caller decides whether to push a [`PipelineStatsFinalizeHook`]
/// based on `user_wants_stats` (which is separate from the deadlock
/// monitor's need for stats counters — the monitor wants the counters
/// attached but doesn't need the post-run log).
///
/// [`PipelineStatsFinalizeHook`]: crate::pipeline::chains::PipelineStatsFinalizeHook
///
/// # Errors
///
/// Returns an error if `queue_memory.calculate_memory_limit` fails
/// (e.g. invalid memory limit specification).
pub(crate) fn build_pipeline_config_for_chain(
    pipeline: &Pipeline,
    num_threads: usize,
    scheduler: &SchedulerOptions,
    queue_memory: &QueueMemoryOptions,
) -> Result<(PipelineConfig, Option<Arc<PipelineStats>>)> {
    let mut config = PipelineConfig { threads: num_threads, ..Default::default() };
    config.deadlock_timeout_secs = scheduler.deadlock_timeout_secs();
    config.queue_memory_total = Some(queue_memory.calculate_memory_limit(num_threads)?);
    // Per-edge instrumentation (`--pipeline-trace` / `FGUMI_PIPELINE_TRACE`),
    // threaded through the command-layer `SchedulerOptions`. Left at `Off` / `None`
    // when the flag/env are unset, which keeps `PipelineConfig::default()`'s
    // zero-overhead path (no sampler thread, non-instrumented queues).
    config.instrumentation = scheduler.instrumentation_level();
    config.trace_path = scheduler.trace_path();
    // A `--pipeline-trace-out` path is only consulted at `timeline`/`deep` (the
    // levels that write the per-tick TSV). If a path was given but the *resolved*
    // level (after the `FGUMI_PIPELINE_TRACE` override) does not write one, the
    // flag is a silent no-op — warn rather than drop it without a trace. Checking
    // `config.instrumentation` here, not the raw flag, keeps the env override in
    // scope, so an env-enabled `timeline` run with a path does not false-warn.
    if config.trace_path.is_some() && !config.instrumentation.timeline() {
        log::warn!(
            "--pipeline-trace-out was set but the resolved trace level ({:?}) does not write a \
             timeline TSV (only `timeline` or `deep` do); the path is ignored",
            config.instrumentation
        );
    }

    // Per-tick thread/bandwidth telemetry (`--pipeline-telemetry-out` /
    // `FGUMI_PIPELINE_TELEMETRY_OUT`). Telemetry needs registered edges (edges
    // register only when instrumentation is on), so bump `Off` to `Summary` when
    // telemetry is requested; an explicit `--pipeline-trace` level set above wins.
    config.telemetry = scheduler.telemetry_config();
    if config.telemetry.is_some() && config.instrumentation == InstrumentationLevel::Off {
        config.instrumentation = InstrumentationLevel::Summary;
    }
    if config.telemetry.is_some() {
        config.rss_probe = Some(std::sync::Arc::new(fgumi_sort::process_rss_bytes));
    }

    let user_wants_stats = scheduler.collect_stats();
    let monitor_needs_stats = config.deadlock_timeout_secs > 0;
    // Tracing implies stats: the end-of-run edge table + verdict render through
    // the stats snapshot, and the wall clock comes from PipelineStats.
    let trace_wants_stats = config.instrumentation.is_on();
    // DIAGNOSTIC: `FGUMI_PIPELINE_STATS=1` force-creates the stats Arc even on
    // paths (standalone `fgumi sort`) that neither flatten `--pipeline-stats`
    // nor enable the deadlock monitor, so the per-step timing report can be dumped.
    let env_wants_stats = env_flag_enabled("FGUMI_PIPELINE_STATS");
    let stats = if user_wants_stats || monitor_needs_stats || env_wants_stats || trace_wants_stats {
        let timing = stats_timing(user_wants_stats, env_wants_stats, trace_wants_stats);
        let s = pipeline.stats_with_timing(timing);
        config.stats = Some(Arc::clone(&s));
        Some(s)
    } else {
        None
    };

    Ok((config, stats))
}

/// The clock to time a chain's stats on. Stats asked for explicitly
/// (`--pipeline-stats`, `FGUMI_PIPELINE_STATS`, tracing) are timed precisely.
/// Stats that exist only for the deadlock monitor's stall snapshot, and the
/// end-of-run table that comes with them, are on by default, so they use the
/// coarse clock: a precise read on every dispatch is a system call on a host
/// whose clock source user space cannot read (e.g. `hpet`).
///
/// Takes each reason separately (rather than their union) so a test can pin
/// every route to precise timing without touching the process environment.
fn stats_timing(
    user_wants_stats: bool,
    env_wants_stats: bool,
    trace_wants_stats: bool,
) -> TimingResolution {
    if user_wants_stats || env_wants_stats || trace_wants_stats {
        TimingResolution::Precise
    } else {
        TimingResolution::Coarse
    }
}

#[cfg(test)]
mod tests {
    use std::ffi::OsStr;
    use std::io;
    use std::sync::Arc;

    use rstest::rstest;

    use super::{
        TimingResolution, build_pipeline_config_for_chain, env_flag_enabled, flag_os_value_enabled,
        stats_timing,
    };
    use crate::commands::common::{
        DEFAULT_DEADLOCK_TIMEOUT_SECS, DeadlockTimeoutOptions, QueueMemoryOptions, SchedulerOptions,
    };
    use crate::pipeline::core::builder::{InstrumentationLevel, Pipeline};
    use crate::pipeline::core::outputs::Single;
    use crate::pipeline::core::queues::QueueSpec;
    use crate::pipeline::core::reorder::BranchOrdering;
    use crate::pipeline::core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};

    /// Any one reason to want stats (`--pipeline-stats`, `FGUMI_PIPELINE_STATS`,
    /// tracing or telemetry) makes them precise; none (the monitor alone) leaves
    /// them coarse.
    #[rstest]
    fn stats_timing_is_precise_only_when_requested(
        #[values(false, true)] user: bool,
        #[values(false, true)] env: bool,
        #[values(false, true)] trace: bool,
    ) {
        let expected =
            if user || env || trace { TimingResolution::Precise } else { TimingResolution::Coarse };
        assert_eq!(stats_timing(user, env, trace), expected);
    }

    /// A source that finishes at once, so the two-step chain below builds.
    struct StubSource;
    impl Step for StubSource {
        type Input = ();
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "Source",
                kind: StepKind::Exclusive,
                sticky: true,
                output_queues: vec![QueueSpec::CountBounded { capacity: 64 }],
                branch_ordering: vec![BranchOrdering::ByOrdinal],
            }
        }
        fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            Ok(StepOutcome::Finished)
        }
    }

    struct StubSink;
    impl Step for StubSink {
        type Input = u32;
        type Outputs = ();
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "Sink",
                kind: StepKind::Exclusive,
                sticky: false,
                output_queues: vec![],
                branch_ordering: vec![],
            }
        }
        fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            Ok(StepOutcome::NoProgress)
        }
    }

    /// The clock `build_pipeline_config_for_chain` actually attaches: coarse
    /// when the stats exist only for the deadlock monitor (the default), precise
    /// when they are asked for with `--pipeline-stats` or by tracing.
    #[rstest]
    #[case::monitor_only_is_coarse(
        false,
        InstrumentationLevel::Off,
        false,
        TimingResolution::Coarse
    )]
    #[case::pipeline_stats_is_precise(
        true,
        InstrumentationLevel::Off,
        false,
        TimingResolution::Precise
    )]
    #[case::tracing_is_precise(
        false,
        InstrumentationLevel::Summary,
        false,
        TimingResolution::Precise
    )]
    #[case::telemetry_is_precise(false, InstrumentationLevel::Off, true, TimingResolution::Precise)]
    fn chain_config_times_stats_on_the_expected_clock(
        #[case] pipeline_stats: bool,
        #[case] pipeline_trace: InstrumentationLevel,
        #[case] telemetry: bool,
        #[case] expected: TimingResolution,
    ) {
        // These variables override the options below (any set
        // `FGUMI_PIPELINE_TRACE` value replaces `--pipeline-trace`; either
        // telemetry variable turns telemetry, and so tracing, on), and the
        // crate cannot clear them (`remove_var` is unsafe). This test checks
        // the option-driven choice, so it skips, saying so, when one is set.
        let overriding: Vec<&str> = [
            "FGUMI_PIPELINE_TRACE",
            "FGUMI_PIPELINE_TELEMETRY_OUT",
            "FGUMI_PIPELINE_TELEMETRY_INTERVAL",
        ]
        .into_iter()
        .filter(|v| std::env::var_os(v).is_some())
        .chain(env_flag_enabled("FGUMI_PIPELINE_STATS").then_some("FGUMI_PIPELINE_STATS"))
        .collect();
        if !overriding.is_empty() {
            eprintln!("skipping: {} set in the environment", overriding.join(", "));
            return;
        }
        let builder = Pipeline::builder();
        builder.chain(StubSource).chain(StubSink).into_sink_marker();
        let pipeline = builder.build().expect("chain builds");
        let scheduler = SchedulerOptions {
            pipeline_stats,
            pipeline_trace,
            // Only enables telemetry; nothing is written until a run.
            pipeline_telemetry_out: telemetry
                .then(|| std::env::temp_dir().join("telemetry-unused")),
            deadlock: DeadlockTimeoutOptions { deadlock_timeout: DEFAULT_DEADLOCK_TIMEOUT_SECS },
            ..SchedulerOptions::default()
        };
        let (config, stats) = build_pipeline_config_for_chain(
            &pipeline,
            2,
            &scheduler,
            &QueueMemoryOptions::default(),
        )
        .expect("config builds");
        let stats = stats.expect("the deadlock monitor attaches stats");
        assert_eq!(stats.timing(), expected);
        let attached = config.stats.expect("the config carries the same stats");
        assert!(Arc::ptr_eq(&stats, &attached));
    }

    #[test]
    fn flag_value_unset_is_disabled() {
        assert!(!flag_os_value_enabled(None));
    }

    #[test]
    fn flag_value_zero_and_false_are_disabled() {
        assert!(!flag_os_value_enabled(Some(OsStr::new("0"))));
        assert!(!flag_os_value_enabled(Some(OsStr::new("false"))));
    }

    #[test]
    fn flag_value_other_values_are_enabled() {
        // Any present value that is not exactly "0"/"false" enables — including
        // an empty string, per the documented contract.
        for v in ["1", "true", "yes", "on", "", "FALSE", "00"] {
            assert!(flag_os_value_enabled(Some(OsStr::new(v))), "value {v:?} should enable");
        }
    }
}

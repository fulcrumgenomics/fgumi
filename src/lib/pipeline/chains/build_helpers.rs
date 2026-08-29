//! Shared helpers for per-command chain builders.

use std::sync::Arc;

use anyhow::Result;

use crate::commands::common::{QueueMemoryOptions, SchedulerOptions};
use crate::pipeline::core::builder::{Pipeline, PipelineConfig};
use crate::pipeline::core::runtime::stats::PipelineStats;

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
    // `main-runall`'s command-layer `SchedulerOptions` does not carry per-edge
    // instrumentation (`--pipeline-trace`); `config.instrumentation` / `trace_path`
    // stay at their `Off` / `None` defaults from `PipelineConfig::default()`.

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
        let s = pipeline.stats();
        config.stats = Some(Arc::clone(&s));
        Some(s)
    } else {
        None
    };

    Ok((config, stats))
}

#[cfg(test)]
mod tests {
    use std::ffi::OsStr;

    use super::flag_os_value_enabled;

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

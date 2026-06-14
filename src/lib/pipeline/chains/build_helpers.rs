//! Shared helpers for per-command chain builders.

use std::sync::Arc;

use anyhow::Result;

use crate::commands::common::{QueueMemoryOptions, SchedulerOptions};
use crate::pipeline::core::builder::{Pipeline, PipelineConfig};
use crate::pipeline::core::runtime::stats::PipelineStats;

/// Build the [`PipelineConfig`] for a chain run and optionally attach
/// [`PipelineStats`]. Centralizes the threads + deadlock-timeout +
/// queue-memory + stats wiring that every per-command builder needs
/// identically.
///
/// Returns the constructed config and the `Arc<PipelineStats>` if stats
/// were attached (`Some` when either `--pipeline-stats` is set OR
/// `--deadlock-timeout > 0`, since the deadlock monitor reads stats
/// counters to detect stalls).
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
pub fn build_pipeline_config_for_chain(
    pipeline: &Pipeline,
    num_threads: usize,
    scheduler: &SchedulerOptions,
    queue_memory: &QueueMemoryOptions,
) -> Result<(PipelineConfig, Option<Arc<PipelineStats>>)> {
    let mut config = PipelineConfig { threads: num_threads, ..Default::default() };
    config.deadlock_timeout_secs = scheduler.deadlock_timeout_secs();
    config.queue_memory_total = Some(queue_memory.calculate_memory_limit(num_threads)?);

    let user_wants_stats = scheduler.collect_stats();
    let monitor_needs_stats = config.deadlock_timeout_secs > 0;
    let stats = if user_wants_stats || monitor_needs_stats {
        let s = pipeline.stats();
        config.stats = Some(Arc::clone(&s));
        Some(s)
    } else {
        None
    };

    Ok((config, stats))
}

//! Post-pipeline cleanup hooks.

use std::sync::Arc;

use anyhow::Result;

use crate::pipeline::core::builder::{Pipeline, PipelineConfig};
use crate::pipeline::core::runtime::stats::PipelineStats;

/// Post-pipeline cleanup. Each stage with metrics, summary logging,
/// `--min-corrected`-style gates, or rejects-file finalization
/// registers one or more hooks during chain build. The caller iterates
/// `built.finalize.into_iter().try_for_each(FinalizeHook::finalize)`
/// after `pipeline.run()` returns.
///
/// Trait object so heterogeneous hooks (correct's metrics drain,
/// consensus's summary log, AAM's records-aligned counter, etc.) can
/// share a `Vec<Box<dyn FinalizeHook>>`.
pub trait FinalizeHook: Send {
    /// Run the post-pipeline action. Called exactly once, after
    /// `Pipeline::run` returns and before the caller exits.
    ///
    /// # Errors
    ///
    /// Returns the underlying I/O / aggregation / gate-check error.
    /// Callers typically chain via `try_for_each`.
    fn finalize(self: Box<Self>) -> Result<()>;
}

/// Logs [`PipelineStats`] after the pipeline finishes, when
/// `--pipeline-stats` was requested. Builders construct this and push
/// it onto [`BuiltPipeline::finalize`] as the first hook (so it logs
/// before command-specific finalize hooks).
pub struct PipelineStatsFinalizeHook {
    pub stats: Arc<PipelineStats>,
}

impl FinalizeHook for PipelineStatsFinalizeHook {
    fn finalize(self: Box<Self>) -> Result<()> {
        crate::commands::common::log_new_pipeline_stats(Some(self.stats));
        Ok(())
    }
}

/// Logs wall-time elapsed since the hook was constructed. One instance is
/// registered per chain — by [`ChainBuilder::build`], not by individual
/// `add_<stage>` methods — so for a fused multi-stage chain like
/// `[Group, Simplex]` a single line is logged rather than two lines each
/// reporting the full pipeline duration.
///
/// The `chain_label` identifies the chain shape in the log line, e.g.
/// `"dedup"` for a single-stage dedup chain or `"group→simplex"` for a
/// fused group+simplex chain.
pub struct StageTimingFinalizeHook {
    pub chain_label: String,
    pub start: std::time::Instant,
}

impl StageTimingFinalizeHook {
    /// Construct the hook with a static stage name and `start` captured now.
    ///
    /// Kept for unit tests and any call site that still uses static names;
    /// production code should prefer [`Self::new_with_label`].
    #[must_use]
    pub fn new(stage_name: &'static str) -> Self {
        Self { chain_label: stage_name.to_owned(), start: std::time::Instant::now() }
    }

    /// Construct the hook with an owned `chain_label` string and `start`
    /// captured at call time. [`ChainBuilder::build`] uses this form,
    /// generating the label from `spec.stages` (e.g. `"dedup"` or
    /// `"group→simplex"`).
    #[must_use]
    pub fn new_with_label(chain_label: String) -> Self {
        Self { chain_label, start: std::time::Instant::now() }
    }
}

impl FinalizeHook for StageTimingFinalizeHook {
    fn finalize(self: Box<Self>) -> Result<()> {
        let elapsed = self.start.elapsed();
        log::info!("Pipeline `{}` ran in {:.2}s", self.chain_label, elapsed.as_secs_f64());
        Ok(())
    }
}

/// Output of [`crate::pipeline::chains::build_for`]. The
/// caller calls [`BuiltPipeline::run`] which executes the pipeline and
/// drains the finalize hooks in order.
pub struct BuiltPipeline {
    pub pipeline: Pipeline,
    pub config: PipelineConfig,
    pub finalize: Vec<Box<dyn FinalizeHook>>,
}

impl BuiltPipeline {
    /// Convenience: run the pipeline and drain finalize hooks in order.
    ///
    /// Finalize hooks are *always* drained in a finally-style path, even
    /// when `Pipeline::run` fails: the hooks own the metrics drains,
    /// summary logging, `--min-corrected`-style gates, and rejects-file
    /// finalization that must happen regardless of how the run ended.
    /// Every hook is run even if an earlier hook errors, so one failing
    /// hook cannot prevent later cleanup. The first error encountered
    /// (pipeline run first, then hooks in registration order) is
    /// returned.
    ///
    /// # Errors
    ///
    /// Returns the first pipeline-run error or, if the pipeline
    /// succeeded, the first finalize-hook error.
    pub fn run(self) -> Result<()> {
        let BuiltPipeline { pipeline, config, finalize } = self;

        // Capture the run result but do not short-circuit: the finalize
        // hooks must still drain on the error path.
        let run_result = pipeline.run(config).map_err(|e| anyhow::anyhow!("Pipeline::run: {e:?}"));

        drain_finalize(run_result, finalize)
    }
}

/// Drains every finalize hook in registration order, then combines the
/// pipeline-run result with the first hook error.
///
/// Hooks are always drained — including on the `run_result == Err` path —
/// and a failing hook does not prevent later hooks from running. The
/// pipeline-run error takes precedence over any hook error.
fn drain_finalize(run_result: Result<()>, finalize: Vec<Box<dyn FinalizeHook>>) -> Result<()> {
    // Drain every hook, keeping the first hook error.
    let mut first_hook_error: Option<anyhow::Error> = None;
    for hook in finalize {
        if let Err(e) = hook.finalize() {
            if first_hook_error.is_none() {
                first_hook_error = Some(e);
            }
        }
    }

    // The pipeline-run error takes precedence over hook errors.
    run_result.and_then(|()| first_hook_error.map_or(Ok(()), Err))
}

#[cfg(test)]
mod tests {
    use std::sync::atomic::{AtomicUsize, Ordering};

    use super::*;

    /// Test hook that bumps a shared counter on finalize and optionally
    /// returns an error, so tests can assert every hook ran.
    struct CountingHook {
        ran: Arc<AtomicUsize>,
        fail: bool,
    }

    impl FinalizeHook for CountingHook {
        fn finalize(self: Box<Self>) -> Result<()> {
            self.ran.fetch_add(1, Ordering::SeqCst);
            if self.fail { Err(anyhow::anyhow!("hook failed")) } else { Ok(()) }
        }
    }

    fn counting_hook(ran: &Arc<AtomicUsize>, fail: bool) -> Box<dyn FinalizeHook> {
        Box::new(CountingHook { ran: Arc::clone(ran), fail })
    }

    #[test]
    fn drain_finalize_runs_all_hooks_when_pipeline_fails() {
        // Even though the pipeline "failed", all hooks must still run so
        // metrics drains / rejects finalization are not skipped.
        let ran = Arc::new(AtomicUsize::new(0));
        let hooks = vec![counting_hook(&ran, false), counting_hook(&ran, false)];
        let result = drain_finalize(Err(anyhow::anyhow!("pipeline boom")), hooks);
        assert_eq!(ran.load(Ordering::SeqCst), 2, "all hooks should run on the error path");
        // The pipeline error takes precedence and is propagated.
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("pipeline boom"));
    }

    #[test]
    fn drain_finalize_runs_later_hooks_after_a_failing_hook() {
        // A failing hook must not prevent later cleanup hooks from running.
        let ran = Arc::new(AtomicUsize::new(0));
        let hooks =
            vec![counting_hook(&ran, true), counting_hook(&ran, false), counting_hook(&ran, false)];
        let result = drain_finalize(Ok(()), hooks);
        assert_eq!(ran.load(Ordering::SeqCst), 3, "all hooks should run despite an early failure");
        // With a successful run, the first hook error is surfaced.
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("hook failed"));
    }

    #[test]
    fn drain_finalize_ok_when_run_and_hooks_succeed() {
        let ran = Arc::new(AtomicUsize::new(0));
        let hooks = vec![counting_hook(&ran, false)];
        let result = drain_finalize(Ok(()), hooks);
        assert_eq!(ran.load(Ordering::SeqCst), 1);
        assert!(result.is_ok());
    }

    #[test]
    fn stage_timing_hook_finalize_succeeds() {
        // Construct the hook with a test stage name and immediately
        // finalize. Verifies the hook is wired correctly and produces
        // no error. The "Pipeline `X` ran in Ys" log line is emitted via
        // log::info!; the call path is exercised by this test.
        let hook = StageTimingFinalizeHook::new("test_stage");
        assert_eq!(hook.chain_label, "test_stage");
        let result = Box::new(hook).finalize();
        assert!(result.is_ok(), "StageTimingFinalizeHook::finalize should succeed: {result:?}");
    }

    #[test]
    fn stage_timing_hook_new_with_label() {
        // Verify that new_with_label stores the label correctly.
        let hook = StageTimingFinalizeHook::new_with_label("group→simplex".to_owned());
        assert_eq!(hook.chain_label, "group→simplex");
        let result = Box::new(hook).finalize();
        assert!(
            result.is_ok(),
            "StageTimingFinalizeHook::finalize (new_with_label) should succeed: {result:?}"
        );
    }

    #[test]
    fn pipeline_stats_hook_finalize_succeeds() {
        // Construct a PipelineStatsFinalizeHook with minimal stats
        // and finalize it. Verifies the hook wiring and that
        // finalize() succeeds without error.
        let stats: Arc<PipelineStats> = Arc::new(PipelineStats::new(vec!["test_step"]));
        let hook = PipelineStatsFinalizeHook { stats };
        let result = Box::new(hook).finalize();
        assert!(result.is_ok(), "PipelineStatsFinalizeHook::finalize should succeed: {result:?}");
    }
}

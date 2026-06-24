//! Post-pipeline cleanup hooks.

use std::sync::Arc;

use anyhow::Result;

use crate::pipeline::core::builder::{Pipeline, PipelineConfig};
use crate::pipeline::core::runtime::stats::PipelineStats;

// The bare `FinalizeHook` trait lives in `fgumi-pipeline-core` so CLI crates
// built directly on the pipeline core (e.g. `fgumi-sort-cli`) can implement it
// without redefining a parallel trait (X1-005). It is re-exported here so the
// canonical `crate::pipeline::chains::FinalizeHook` path — used by every
// in-crate `impl FinalizeHook` and by `BuiltPipeline` below — is unchanged.
pub use crate::pipeline::core::FinalizeHook;

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
    /// Hooks that run **only after a fully successful run** — the pipeline
    /// completed and every always-drain `finalize` hook succeeded.
    ///
    /// Use this for actions that publish a derived artifact from the output
    /// the run produced (e.g. writing a `.bai` sidecar by re-reading the
    /// finished BAM). On the error path the output is incomplete, so running
    /// these would publish a stale/partial artifact — exactly the
    /// `IndexBamFinalizeHook` footgun the standalone `fgumi sort` flow guards
    /// against by gating the index write behind `run_result?`.
    pub finalize_on_success: Vec<Box<dyn FinalizeHook>>,
}

impl BuiltPipeline {
    /// Convenience: run the pipeline and drain finalize hooks in order.
    ///
    /// `finalize` hooks are *always* drained in a finally-style path, even
    /// when `Pipeline::run` fails: the hooks own the metrics drains,
    /// summary logging, `--min-corrected`-style gates, and rejects-file
    /// finalization that must happen regardless of how the run ended.
    /// Every hook is run even if an earlier hook errors, so one failing
    /// hook cannot prevent later cleanup. The first error encountered
    /// (pipeline run first, then hooks in registration order) is
    /// returned.
    ///
    /// `finalize_on_success` hooks run only once the run and every
    /// always-drain hook have succeeded — see the field docs.
    ///
    /// # Errors
    ///
    /// Returns the first pipeline-run error or, if the pipeline
    /// succeeded, the first finalize-hook error.
    pub fn run(self) -> Result<()> {
        let BuiltPipeline { pipeline, config, finalize, finalize_on_success } = self;

        // Capture the run result but do not short-circuit: the finalize
        // hooks must still drain on the error path.
        let run_result = pipeline.run(config).map_err(|e| anyhow::anyhow!("Pipeline::run: {e:?}"));

        drain_finalize(run_result, finalize, finalize_on_success)
    }
}

/// Drains every always-run finalize hook in registration order, then — only
/// if the run and all of those hooks succeeded — runs the success-gated hooks.
///
/// Always-run hooks are drained even on the `run_result == Err` path, and a
/// failing hook does not prevent later always-run hooks from running. The
/// pipeline-run error takes precedence over any hook error. Success-gated
/// hooks are skipped entirely unless everything above succeeded; once running,
/// the first of them to error short-circuits the rest (they publish derived
/// artifacts, so there is no "drain everything" obligation on their path).
fn drain_finalize(
    run_result: Result<()>,
    finalize: Vec<Box<dyn FinalizeHook>>,
    finalize_on_success: Vec<Box<dyn FinalizeHook>>,
) -> Result<()> {
    // Drain every always-run hook, keeping the first hook error.
    let mut first_hook_error: Option<anyhow::Error> = None;
    for hook in finalize {
        if let Err(e) = hook.finalize() {
            if first_hook_error.is_none() {
                first_hook_error = Some(e);
            }
        }
    }

    // The pipeline-run error takes precedence over hook errors. Resolve the
    // combined result of the run plus the always-run hooks first.
    run_result.and_then(|()| first_hook_error.map_or(Ok(()), Err))?;

    // Reached only on full success: now safe to publish derived artifacts.
    for hook in finalize_on_success {
        hook.finalize()?;
    }
    Ok(())
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
        let result = drain_finalize(Err(anyhow::anyhow!("pipeline boom")), hooks, vec![]);
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
        let result = drain_finalize(Ok(()), hooks, vec![]);
        assert_eq!(ran.load(Ordering::SeqCst), 3, "all hooks should run despite an early failure");
        // With a successful run, the first hook error is surfaced.
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("hook failed"));
    }

    #[test]
    fn drain_finalize_ok_when_run_and_hooks_succeed() {
        let ran = Arc::new(AtomicUsize::new(0));
        let hooks = vec![counting_hook(&ran, false)];
        let result = drain_finalize(Ok(()), hooks, vec![]);
        assert_eq!(ran.load(Ordering::SeqCst), 1);
        assert!(result.is_ok());
    }

    #[test]
    fn drain_finalize_runs_success_hooks_only_on_full_success() {
        // Success-gated hooks run once the pipeline and all always-run hooks
        // succeed.
        let always = Arc::new(AtomicUsize::new(0));
        let on_success = Arc::new(AtomicUsize::new(0));
        let result = drain_finalize(
            Ok(()),
            vec![counting_hook(&always, false)],
            vec![counting_hook(&on_success, false)],
        );
        assert!(result.is_ok());
        assert_eq!(always.load(Ordering::SeqCst), 1);
        assert_eq!(on_success.load(Ordering::SeqCst), 1, "success hook should run on full success");
    }

    #[test]
    fn drain_finalize_skips_success_hooks_when_pipeline_fails() {
        // A failed run must not run success-gated hooks — the output is
        // incomplete, so e.g. a `.bai` sidecar would be stale. The always-run
        // hooks still drain.
        let always = Arc::new(AtomicUsize::new(0));
        let on_success = Arc::new(AtomicUsize::new(0));
        let result = drain_finalize(
            Err(anyhow::anyhow!("pipeline boom")),
            vec![counting_hook(&always, false)],
            vec![counting_hook(&on_success, false)],
        );
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("pipeline boom"));
        assert_eq!(always.load(Ordering::SeqCst), 1, "always-run hook still drains on failure");
        assert_eq!(on_success.load(Ordering::SeqCst), 0, "success hook must not run on failure");
    }

    #[test]
    fn drain_finalize_skips_success_hooks_when_an_always_hook_fails() {
        // If an always-run hook fails (even with a successful pipeline run),
        // the output is not trustworthy, so success-gated hooks are skipped.
        let always = Arc::new(AtomicUsize::new(0));
        let on_success = Arc::new(AtomicUsize::new(0));
        let result = drain_finalize(
            Ok(()),
            vec![counting_hook(&always, true)],
            vec![counting_hook(&on_success, false)],
        );
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("hook failed"));
        assert_eq!(
            on_success.load(Ordering::SeqCst),
            0,
            "success hook must not run after a failure"
        );
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

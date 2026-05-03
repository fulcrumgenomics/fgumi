//! Chain-runner dispatch trait and registry for `fgumi runall`.
//!
//! `fgumi runall` runs many `(start_from, stop_after)` chain shapes,
//! each of which may want a bespoke pipeline implementation (for
//! example, `(Extract, Extract)` runs purely on
//! [`crate::unified_pipeline::fastq`] while `(Correct, Correct)` runs
//! purely on [`crate::unified_pipeline::bam`]). Rather than encoding
//! the dispatch logic in a fast-growing `if/else` ladder inside
//! `execute()`, runall uses a [`ChainRegistry`] of trait objects:
//!
//! ```text
//!                ┌──────────────┐
//!  CLI args ──>  │ ChainRegistry │  ──>  ChainRunner::run(ctx)
//!                └──────────────┘
//!                  ↑          │
//!                  └─ first runner whose
//!                     ChainRunner::supports(...) returns true
//! ```
//!
//! PR 1 wires a single fallback runner ([`crate::commands::runall::chains::NotImplementedRunner`])
//! that matches every shape and bails with the documented
//! "not yet implemented" error. PR 2/3 push real runners onto the front
//! of the registry; the fallback handles every shape they don't.

use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::AtomicBool;

use anyhow::{Result, bail};

use super::options::{StartFrom, StopAfter};

/// Cancel token shared between the runall execute loop and the running
/// chain. The signal handler installed by
/// [`super::infra::signal::install_signal_handler_once`] flips the
/// underlying flag on SIGINT/SIGTERM; chain runners must poll it
/// during long-running work and return cleanly when it goes high.
///
/// PR 1 only registers the [`crate::commands::runall::chains::NotImplementedRunner`]
/// fallback, which never polls the token; the helpers below are
/// therefore unused in PR 1 and `#[allow(dead_code)]`-gated until a
/// real chain runner consumes them.
#[derive(Clone, Debug)]
pub(crate) struct CancelToken {
    flag: Arc<AtomicBool>,
}

impl CancelToken {
    /// Wrap an existing shared flag (used by `signal::install_signal_handler_once`).
    pub(crate) fn from_flag(flag: Arc<AtomicBool>) -> Self {
        Self { flag }
    }

    /// Has cancellation been requested?
    #[allow(dead_code)]
    pub(crate) fn is_cancelled(&self) -> bool {
        self.flag.load(std::sync::atomic::Ordering::SeqCst)
    }

    /// Borrow the shared flag (used by `execute()` for the post-run check).
    #[allow(dead_code)]
    pub(crate) fn flag(&self) -> &Arc<AtomicBool> {
        &self.flag
    }
}

/// CLI-level inputs that every runner uses to decide whether it
/// supports the current invocation. Borrowed (not owned) so we don't
/// have to clone the `Runall` struct just to dispatch.
///
/// `metrics_path` is consumed by chain runners that need to gate
/// themselves off `--metrics` (since metrics support arrives chain
/// by chain). PR 1's [`crate::commands::runall::chains::NotImplementedRunner`]
/// ignores it; `#[allow(dead_code)]` keeps the field warning-free
/// until a real runner reads it.
#[derive(Debug)]
pub(crate) struct DispatchContext<'a> {
    /// Where in the pipeline graph the chain starts.
    pub(crate) start_from: StartFrom,
    /// Where in the pipeline graph the chain stops (inclusive).
    pub(crate) stop_after: StopAfter,
    /// Whether the user passed `--metrics`. Many chain runners gate
    /// themselves off when metrics are requested but not yet wired up,
    /// falling through to the next runner (or the
    /// `NotImplementedRunner` fallback).
    #[allow(dead_code)]
    pub(crate) metrics_path: Option<&'a Path>,
}

/// Per-invocation context handed to a chain runner's `run()` method.
///
/// Owns the staging-output path (the `.tmp` file under the user's
/// `--output`) and the cancel token. Borrows everything CLI-shaped
/// from [`DispatchContext`] so a runner has the same view of the
/// invocation that the registry used to choose it.
///
/// PR 1's [`crate::commands::runall::chains::NotImplementedRunner`]
/// only consults `dispatch`; the other fields are
/// `#[allow(dead_code)]`-gated until a real runner consumes them.
#[derive(Debug)]
pub(crate) struct ChainContext<'a> {
    /// What the registry used to pick this runner.
    pub(crate) dispatch: DispatchContext<'a>,
    /// Path the runner must write its output to (`.bam.tmp`).
    #[allow(dead_code)]
    pub(crate) tmp_output: PathBuf,
    /// Cancel token: poll periodically and return `Ok(())` on cancel.
    #[allow(dead_code)]
    pub(crate) cancel: CancelToken,
    /// Captured `argv` for `@PG` records and explain output.
    #[allow(dead_code)]
    pub(crate) command_line: &'a str,
}

/// A chain-runner implementation.
///
/// Runners are registered in a [`ChainRegistry`] and consulted in
/// insertion order — the first runner whose `supports(...)` returns
/// `true` is asked to `run(...)` the chain.
pub(crate) trait ChainRunner: Send + Sync {
    /// Stable short name for `--explain` output and error messages.
    fn name(&self) -> &'static str;

    /// Returns `true` if this runner handles the given dispatch
    /// context. Runners that gate on flag combinations (e.g.
    /// "I don't yet support `--metrics`") return `false` for the
    /// gated cases so the registry falls through to the next runner.
    fn supports(&self, ctx: &DispatchContext<'_>) -> bool;

    /// Run the chain end-to-end against the prepared `ChainContext`.
    /// On cancellation, return `Ok(())`; the caller (in `execute()`)
    /// inspects the cancel flag separately and converts a cancelled
    /// `Ok(())` into an `Err` so the atomic-output guard discards
    /// the temp file.
    fn run(&self, ctx: ChainContext<'_>) -> Result<()>;
}

/// Insertion-ordered registry of chain runners.
///
/// Runners are consulted in registration order; the first one whose
/// `supports(...)` returns `true` handles the dispatch. The PR 1
/// fallback ([`crate::commands::runall::chains::NotImplementedRunner`])
/// always sits last and matches every shape it sees, so the registry
/// always finds *some* runner — the only way `dispatch()` returns the
/// "not yet implemented" `bail!` is if the caller deliberately built
/// a registry without the fallback (used in
/// `tests::registry_falls_back_when_no_runner_matches`).
pub(crate) struct ChainRegistry {
    runners: Vec<Box<dyn ChainRunner>>,
}

impl ChainRegistry {
    /// Build an empty registry.
    pub(crate) fn new() -> Self {
        Self { runners: Vec::new() }
    }

    /// Push a runner onto the back of the registry. Earlier-registered
    /// runners are consulted first.
    pub(crate) fn register(mut self, runner: Box<dyn ChainRunner>) -> Self {
        self.runners.push(runner);
        self
    }

    /// Find the first runner whose `supports(...)` is true and run it.
    /// Returns the documented "not yet implemented" error if no runner
    /// matches.
    pub(crate) fn dispatch(&self, ctx: ChainContext<'_>) -> Result<()> {
        for runner in &self.runners {
            if runner.supports(&ctx.dispatch) {
                log::debug!(
                    "runall dispatching {:?}->{:?} to {}",
                    ctx.dispatch.start_from,
                    ctx.dispatch.stop_after,
                    runner.name(),
                );
                return runner.run(ctx);
            }
        }
        bail!(
            "the chain --start-from {:?} --stop-after {:?} is not yet implemented",
            ctx.dispatch.start_from,
            ctx.dispatch.stop_after,
        )
    }

    /// Find the runner that *would* handle the given dispatch context,
    /// without invoking it. Used by `--explain` to print a chain-aware
    /// runner name.
    pub(crate) fn select(&self, ctx: &DispatchContext<'_>) -> Option<&dyn ChainRunner> {
        self.runners.iter().find(|r| r.supports(ctx)).map(std::convert::AsRef::as_ref)
    }
}

#[cfg(test)]
pub(crate) mod tests {
    //! Unit tests for the dispatch/registry plumbing, plus a
    //! [`TestChainRunner`] that integration tests can register to
    //! exercise dispatch + atomic-output + signal-handling end-to-end
    //! without needing a real chain.

    use super::*;
    use std::sync::Mutex;
    use std::sync::atomic::{AtomicUsize, Ordering};

    /// Test-only chain runner that:
    ///
    /// * matches a configured `(StartFrom, StopAfter)` pair, and
    /// * on `run()`, increments a call counter, writes a one-byte BAM
    ///   stand-in to the temp output path, then returns the configured
    ///   result.
    ///
    /// Lives inside `dispatch::tests` (not in `chains/`) so that
    /// integration tests can pull it in via `#[cfg(test)]` without
    /// polluting the production chain registry.
    pub(crate) struct TestChainRunner {
        name: &'static str,
        start: StartFrom,
        stop: StopAfter,
        rejects_metrics: bool,
        call_count: AtomicUsize,
        result: Mutex<Option<anyhow::Error>>,
    }

    impl TestChainRunner {
        pub(crate) fn new(name: &'static str, start: StartFrom, stop: StopAfter) -> Self {
            Self {
                name,
                start,
                stop,
                rejects_metrics: false,
                call_count: AtomicUsize::new(0),
                result: Mutex::new(None),
            }
        }

        /// Borrow the call counter (lets a test verify the runner
        /// fired exactly once).
        pub(crate) fn call_count(&self) -> usize {
            self.call_count.load(Ordering::SeqCst)
        }
    }

    impl ChainRunner for TestChainRunner {
        fn name(&self) -> &'static str {
            self.name
        }

        fn supports(&self, ctx: &DispatchContext<'_>) -> bool {
            if self.rejects_metrics && ctx.metrics_path.is_some() {
                return false;
            }
            ctx.start_from == self.start && ctx.stop_after == self.stop
        }

        fn run(&self, ctx: ChainContext<'_>) -> Result<()> {
            self.call_count.fetch_add(1, Ordering::SeqCst);
            std::fs::write(&ctx.tmp_output, b"BAM-stub").expect("write tmp output");
            let mut slot = self.result.lock().unwrap();
            if let Some(err) = slot.take() {
                return Err(err);
            }
            Ok(())
        }
    }

    fn dispatch_ctx(start: StartFrom, stop: StopAfter) -> DispatchContext<'static> {
        DispatchContext { start_from: start, stop_after: stop, metrics_path: None }
    }

    fn chain_ctx(start: StartFrom, stop: StopAfter, tmp_output: PathBuf) -> ChainContext<'static> {
        let cancel = CancelToken::from_flag(Arc::new(AtomicBool::new(false)));
        ChainContext { dispatch: dispatch_ctx(start, stop), tmp_output, cancel, command_line: "" }
    }

    #[test]
    fn registry_dispatches_to_first_supporting_runner() {
        let runner_a = TestChainRunner::new("A", StartFrom::Extract, StopAfter::Extract);
        let runner_b = TestChainRunner::new("B", StartFrom::Correct, StopAfter::Correct);

        let registry =
            ChainRegistry::new().register(Box::new(runner_a)).register(Box::new(runner_b));

        // Both runners reject metrics by default; verify selection.
        let ctx_a = dispatch_ctx(StartFrom::Extract, StopAfter::Extract);
        assert_eq!(registry.select(&ctx_a).map(ChainRunner::name), Some("A"));

        let ctx_b = dispatch_ctx(StartFrom::Correct, StopAfter::Correct);
        assert_eq!(registry.select(&ctx_b).map(ChainRunner::name), Some("B"));

        let ctx_unknown = dispatch_ctx(StartFrom::Group, StopAfter::Group);
        assert!(registry.select(&ctx_unknown).is_none());
    }

    #[test]
    fn registry_falls_back_when_no_runner_matches() {
        let registry = ChainRegistry::new(); // empty
        let dir = tempfile::TempDir::new().unwrap();
        let tmp = dir.path().join("out.bam.tmp");
        let err = registry
            .dispatch(chain_ctx(StartFrom::Sort, StopAfter::Filter, tmp.clone()))
            .expect_err("empty registry must bail");
        let msg = format!("{err:#}");
        assert!(msg.contains("not yet implemented"), "got: {msg}");
        assert!(msg.contains("Sort"), "error must include start_from");
        assert!(msg.contains("Filter"), "error must include stop_after");
    }

    #[test]
    fn registry_writes_tmp_output_when_runner_succeeds() {
        let dir = tempfile::TempDir::new().unwrap();
        let tmp = dir.path().join("out.bam.tmp");
        let registry = ChainRegistry::new().register(Box::new(TestChainRunner::new(
            "T",
            StartFrom::Extract,
            StopAfter::Extract,
        )));

        registry
            .dispatch(chain_ctx(StartFrom::Extract, StopAfter::Extract, tmp.clone()))
            .expect("dispatch ok");

        assert!(tmp.exists());
        assert_eq!(std::fs::read(&tmp).unwrap(), b"BAM-stub");
    }
}

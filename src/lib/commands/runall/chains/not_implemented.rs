//! Fallback chain runner that bails with a "not yet implemented" error.
//!
//! In PR 1 there are no real chain runners: every `(start_from,
//! stop_after)` pair flows through this runner, which always
//! `supports(...)` and always returns the documented error from
//! `run(...)`. PRs 2 and 3 push real runners onto the front of the
//! registry; this fallback keeps handling the (strictly larger) set
//! of shapes that haven't been wired up yet.

use anyhow::{Result, bail};

use crate::commands::runall::dispatch::{ChainContext, ChainRunner, DispatchContext};
use crate::commands::runall::infra::explain::NOT_IMPLEMENTED_RUNNER_NAME;

/// Fallback chain runner. See module docs.
pub(crate) struct NotImplementedRunner;

impl ChainRunner for NotImplementedRunner {
    fn name(&self) -> &'static str {
        NOT_IMPLEMENTED_RUNNER_NAME
    }

    /// Always supports — this is the registry's fallback.
    fn supports(&self, _ctx: &DispatchContext<'_>) -> bool {
        true
    }

    /// Always bails with the documented error.
    fn run(&self, ctx: ChainContext<'_>) -> Result<()> {
        bail!(
            "the chain --start-from {:?} --stop-after {:?} is not yet implemented",
            ctx.dispatch.start_from,
            ctx.dispatch.stop_after,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::commands::runall::dispatch::CancelToken;
    use crate::commands::runall::options::{StartFrom, StopAfter};
    use std::path::PathBuf;
    use std::sync::Arc;
    use std::sync::atomic::AtomicBool;

    fn ctx(start: StartFrom, stop: StopAfter) -> ChainContext<'static> {
        ChainContext {
            dispatch: DispatchContext { start_from: start, stop_after: stop, metrics_path: None },
            tmp_output: PathBuf::from("/tmp/unused.bam.tmp"),
            cancel: CancelToken::from_flag(Arc::new(AtomicBool::new(false))),
            command_line: "",
        }
    }

    #[test]
    fn supports_every_chain_shape() {
        let runner = NotImplementedRunner;
        for start in [
            StartFrom::Extract,
            StartFrom::Correct,
            StartFrom::Fastq,
            StartFrom::Sort,
            StartFrom::Group,
        ] {
            for stop in [
                StopAfter::Extract,
                StopAfter::Correct,
                StopAfter::Fastq,
                StopAfter::Align,
                StopAfter::Zipper,
                StopAfter::Sort,
                StopAfter::Group,
                StopAfter::Consensus,
                StopAfter::Filter,
            ] {
                let dctx =
                    DispatchContext { start_from: start, stop_after: stop, metrics_path: None };
                assert!(runner.supports(&dctx));
            }
        }
    }

    #[test]
    fn run_bails_with_chain_in_message() {
        let runner = NotImplementedRunner;
        let err = runner.run(ctx(StartFrom::Extract, StopAfter::Filter)).unwrap_err();
        let msg = format!("{err:#}");
        assert!(msg.contains("not yet implemented"), "got: {msg}");
        assert!(msg.contains("Extract"));
        assert!(msg.contains("Filter"));
    }

    #[test]
    fn name_matches_explain_constant() {
        assert_eq!(NotImplementedRunner.name(), NOT_IMPLEMENTED_RUNNER_NAME);
    }
}

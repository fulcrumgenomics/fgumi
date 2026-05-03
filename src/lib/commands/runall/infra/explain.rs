//! `--explain` text renderer for `fgumi runall`.
//!
//! `--explain` short-circuits the command before validation, file I/O,
//! and signal-handler setup. It prints a concise human-readable
//! description of the chain that *would* run for the given
//! `(start_from, stop_after)` pair, plus the runner that handles it.
//! For chain shapes that are not yet implemented, the explainer still
//! prints the shape and clearly marks the chain as unimplemented so
//! users can see what runall plans to wire up next.

use crate::commands::runall::options::{StartFrom, StopAfter};

/// Render a multi-line explanation of a chain shape.
///
/// `runner_name` is the short name of the chain runner that handles
/// this `(start, stop)` pair (e.g. `"ExtractChainRunner"`), or
/// `"NotImplementedRunner"` for the PR 1 fallback.
pub(crate) fn render_explain(
    start: StartFrom,
    stop: StopAfter,
    runner_name: &str,
    command_line: &str,
) -> String {
    use std::fmt::Write as _;

    let mut out = String::new();
    out.push_str("fgumi runall — chain plan\n");
    let _ = writeln!(out, "  command-line     : {command_line}");
    let _ = writeln!(out, "  --start-from     : {start:?}");
    let _ = writeln!(out, "  --stop-after     : {stop:?}");
    let _ = writeln!(out, "  chain runner     : {runner_name}");
    let _ = writeln!(out, "  stage chain      : {}", describe_stage_chain(start, stop));

    if runner_name == NOT_IMPLEMENTED_RUNNER_NAME {
        out.push_str(
            "\nThis chain is not yet implemented. The runall command \
             accepts the flags but the underlying chain runner has not \
             yet been wired up; the command will exit with a clean \
             error explaining the same.\n",
        );
    }
    out
}

/// Stable name for the PR 1 fallback runner. Kept here (rather than in
/// `chains/not_implemented.rs`) so the explain renderer never has to
/// import the chain module just to identify it.
pub(crate) const NOT_IMPLEMENTED_RUNNER_NAME: &str = "NotImplementedRunner";

/// Describe the canonical stage chain for a `(start, stop)` pair.
///
/// This is purely descriptive — no flag is "decided" by this string —
/// so it stays in sync with the doc-comments on `StartFrom`/`StopAfter`
/// rather than driving any real dispatch.
fn describe_stage_chain(start: StartFrom, stop: StopAfter) -> String {
    let stages = canonical_stages(start, stop);
    stages.join(" -> ")
}

/// Return the stage names a chain is expected to traverse for a given
/// `(start, stop)` shape.
fn canonical_stages(start: StartFrom, stop: StopAfter) -> Vec<&'static str> {
    let all_stages: &[(StopAfter, &'static str)] = &[
        (StopAfter::Extract, "extract"),
        (StopAfter::Correct, "correct"),
        (StopAfter::Fastq, "fastq"),
        (StopAfter::Align, "align"),
        (StopAfter::Zipper, "zipper"),
        (StopAfter::Sort, "sort"),
        (StopAfter::Group, "group"),
        (StopAfter::Consensus, "consensus"),
        (StopAfter::Filter, "filter"),
    ];

    // Find the canonical "first stage" the chain runs based on `start`.
    let first_stop = match start {
        StartFrom::Extract => StopAfter::Extract,
        StartFrom::Correct => StopAfter::Correct,
        StartFrom::Fastq => StopAfter::Fastq,
        StartFrom::Sort => StopAfter::Sort,
        StartFrom::Group => StopAfter::Group,
    };

    all_stages
        .iter()
        .filter(|(s, _)| *s >= first_stop && *s <= stop)
        .map(|(_, name)| *name)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn explain_renders_chain_shape_with_runner_name() {
        let out = render_explain(
            StartFrom::Extract,
            StopAfter::Filter,
            "TestRunner",
            "fgumi runall --start-from extract",
        );
        assert!(out.contains("--start-from     : Extract"));
        assert!(out.contains("--stop-after     : Filter"));
        assert!(out.contains("chain runner     : TestRunner"));
        assert!(out.contains("extract -> correct -> fastq"));
        assert!(out.contains("filter"));
    }

    #[test]
    fn explain_marks_unimplemented_chains() {
        let out = render_explain(
            StartFrom::Sort,
            StopAfter::Filter,
            NOT_IMPLEMENTED_RUNNER_NAME,
            "fgumi runall",
        );
        assert!(out.contains(NOT_IMPLEMENTED_RUNNER_NAME));
        assert!(out.contains("not yet implemented"));
    }

    #[test]
    fn implemented_chain_does_not_print_unimplemented_warning() {
        let out = render_explain(
            StartFrom::Extract,
            StopAfter::Extract,
            "ExtractChainRunner",
            "fgumi runall",
        );
        assert!(!out.contains("not yet implemented"));
    }

    #[test]
    fn canonical_stages_for_extract_to_extract() {
        assert_eq!(canonical_stages(StartFrom::Extract, StopAfter::Extract), vec!["extract"]);
    }

    #[test]
    fn canonical_stages_for_extract_to_correct() {
        assert_eq!(
            canonical_stages(StartFrom::Extract, StopAfter::Correct),
            vec!["extract", "correct"]
        );
    }

    #[test]
    fn canonical_stages_for_sort_to_filter() {
        assert_eq!(
            canonical_stages(StartFrom::Sort, StopAfter::Filter),
            vec!["sort", "group", "consensus", "filter"]
        );
    }

    #[test]
    fn canonical_stages_for_group_to_group() {
        assert_eq!(canonical_stages(StartFrom::Group, StopAfter::Group), vec!["group"]);
    }
}

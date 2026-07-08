//! `fgumi sort` command wiring.
//!
//! The CLI surface (`Sort`, `SortOptions`, `SortOrderArg`, …) is defined in the
//! framework-light `fgumi-sort-cli` crate so its option types can be reused by
//! `runall` (via `MultiSortOptions`). The *execution* path lives here in the
//! umbrella crate so standalone `fgumi sort` runs through the same `build_for`
//! chain builder — and therefore the same work-stealing pool — as every other
//! command. `fgumi-sort-cli` no longer carries a sort execution path of its
//! own — it provides only the CLI surface, the `--verify` read-and-check path,
//! and the `IndexBamFinalizeHook` BAI indexer.

use anyhow::{Result, bail};
use log::warn;

pub use fgumi_sort_cli::sort::{
    MultiSortOptions, Sort, SortOptions, SortOrderArg, TMP_DIRS_ENV, parse_cell_tag,
    resolve_tmp_dirs,
};

use crate::commands::common::{
    DEFAULT_DEADLOCK_TIMEOUT_SECS, QueueMemoryOptions, SchedulerOptions, ThreadingOptions,
};
use crate::pipeline::chains::{ChainSpec, SinkSpec, SourceSpec, Stage, StageOptionsBag, build_for};

/// Run the `fgumi sort` command through the chain builder.
///
/// Validates the CLI options, then either runs the standalone `--verify`
/// read-and-check path (delegated to `fgumi-sort-cli`, which is pre-pipeline
/// and engine-independent) or builds a sole-`[Stage::Sort]` [`ChainSpec`] and
/// executes it via [`build_for`] on the work-stealing pool.
///
/// # Errors
///
/// Returns CLI validation errors or any pipeline build/run error.
pub fn execute_sort_command(cmd: &Sort, command_line: &str) -> Result<()> {
    if cmd.verify && cmd.output.is_some() {
        bail!("--verify cannot be used with --output");
    }
    if cmd.verify && cmd.write_index {
        bail!("--write-index cannot be used with --verify");
    }

    // Validate inputs. Exempt stdin paths: the streaming source reads stdin
    // once. `--verify` re-scans the input, which a non-seekable stdin can't
    // satisfy, so reject stdin there up front.
    if fgumi_bam_io::is_stdin_path(&cmd.input) {
        if cmd.verify {
            bail!(
                "fgumi sort --verify cannot read from stdin (it re-scans the input); \
                 provide a file path instead"
            );
        }
    } else {
        crate::validation::validate_file_exists(&cmd.input, "Input BAM")?;
    }

    if !cmd.verify && cmd.output.is_none() {
        bail!("Either --output or --verify must be specified");
    }

    if cmd.verify {
        // Pre-pipeline read-and-check path — engine-independent, unchanged.
        return cmd.execute_verify();
    }

    if cmd.write_index && !matches!(cmd.order, SortOrderArg::Coordinate) {
        bail!("--write-index is only valid for coordinate sort");
    }

    // Order-independent sort-option semantics (block-batch range,
    // temp-compression/codec compatibility) — the single shared validator that
    // `fgumi runall` also runs, so the two CLI surfaces can't drift.
    cmd.options.validate()?;

    // Order-*dependent* rule, specific to the standalone command's order surface.
    if cmd.options.key_types.is_some() && !matches!(cmd.order, SortOrderArg::TemplateCoordinate) {
        warn!("--key-types is ignored for --order {:?}", cmd.order);
    }

    build_for(build_sort_spec(cmd, command_line))?.run()
}

/// Build the sole-`[Stage::Sort]` [`ChainSpec`] for a standalone `fgumi sort`
/// run. Assumes `cmd` has already passed validation and carries an `--output`.
///
/// Split out from [`execute_sort_command`] so the scheduler / instrumentation
/// wiring is unit-testable without executing a pipeline.
fn build_sort_spec(cmd: &Sort, command_line: &str) -> ChainSpec {
    let output = cmd.output.as_ref().expect("output required for sort mode");

    // Copy order into SortOptions so the chain builder's `add_sort` reads it.
    let mut sort_opts = cmd.options.clone();
    sort_opts.order = cmd.order;

    let sink = if cmd.write_index {
        SinkSpec::BamWithIndex(output.clone())
    } else {
        SinkSpec::Bam(output.clone())
    };

    ChainSpec {
        stages: vec![Stage::Sort],
        source: SourceSpec::Bam(cmd.input.clone()),
        sink,
        stage_opts: StageOptionsBag { sort: Some(sort_opts), ..Default::default() },
        threading: ThreadingOptions { threads: Some(cmd.threads) },
        compression: cmd.compression.clone(),
        // Carry the sort command's diagnostic instrumentation flags (sort is a
        // typed-step pipeline now, so --pipeline-stats / --pipeline-trace apply
        // as they do for the other typed-step commands).
        //
        // `deadlock_timeout` is set explicitly to the shared CLI default:
        // `SchedulerOptions::default()` zeroes it, which *disables* deadlock
        // detection, and `fgumi sort` exposes no `--deadlock-timeout` flag of its
        // own — so without this the fail-fast/diagnostic safety net the field doc
        // calls out for large sort merges would be silently off. Referencing
        // `DEFAULT_DEADLOCK_TIMEOUT_SECS` (the same source of truth as the
        // `--deadlock-timeout` clap default) keeps the two from drifting.
        scheduler: SchedulerOptions {
            pipeline_stats: cmd.pipeline_stats,
            pipeline_trace: cmd.pipeline_trace.clone(),
            pipeline_trace_out: cmd.pipeline_trace_out.clone(),
            deadlock_timeout: DEFAULT_DEADLOCK_TIMEOUT_SECS,
            ..Default::default()
        },
        queue_memory: QueueMemoryOptions::default(),
        async_reader: cmd.async_reader,
        command_line: command_line.to_string(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;

    #[test]
    fn sort_spec_restores_deadlock_timeout_default() {
        // `SchedulerOptions::default()` zeroes `deadlock_timeout` (detection
        // disabled); standalone sort must restore the CLI default so its
        // merge/spill stages keep the fail-fast safety net every other
        // typed-step command has.
        let cmd = Sort::try_parse_from(["fgumi-sort", "-i", "in.bam", "-o", "out.bam"])
            .expect("parse sort args");
        let spec = build_sort_spec(&cmd, "fgumi sort -i in.bam -o out.bam");
        assert_eq!(
            spec.scheduler.deadlock_timeout_secs(),
            10,
            "standalone sort must keep deadlock detection enabled (10s default)"
        );
    }

    #[test]
    fn sort_spec_forwards_instrumentation_flags() {
        let cmd = Sort::try_parse_from([
            "fgumi-sort",
            "-i",
            "in.bam",
            "-o",
            "out.bam",
            "--pipeline-stats",
            "true",
            "--pipeline-trace",
            "summary",
        ])
        .expect("parse sort args");
        let spec = build_sort_spec(&cmd, "cmdline");
        assert!(spec.scheduler.collect_stats(), "--pipeline-stats must flow through");
        assert_eq!(spec.scheduler.pipeline_trace, "summary");
    }
}

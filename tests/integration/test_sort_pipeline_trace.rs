//! Integration tests for the `--pipeline-trace` / `--pipeline-trace-out` /
//! `FGUMI_PIPELINE_TRACE` instrumentation surface on `fgumi sort`.
//!
//! These pin the *honoring* path end-to-end: that the CLI flag (and the env
//! var) actually reach `PipelineConfig::instrumentation` and cause the
//! end-of-run per-edge table + bottleneck verdict to render (and the timeline
//! TSV to be written), rather than being parsed and silently dropped — the
//! characteristic CLI-layer failure this whole change exists to avoid. They run
//! the built binary as a subprocess with its own environment, which is the only
//! way to exercise the `FGUMI_PIPELINE_TRACE` precedence without mutating the
//! test process's environment (`std::env::set_var` is `unsafe` under edition
//! 2024, which this crate forbids).

use std::path::{Path, PathBuf};
use std::process::Command;

use rstest::rstest;
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, create_umi_family, write_bam};

/// Stable substring of the end-of-run per-edge instrumentation report
/// (`fgumi_pipeline_core::builder`, `snapshot_with_edges`) — only ever logged
/// when the resolved instrumentation level is above `Off`.
const EDGE_TABLE_SUBSTRING: &str = "Pipeline edges";

/// Writes a small unsorted UMI BAM. A handful of families is enough: the sort
/// chain has the same per-edge topology regardless of size, so the edge table
/// renders whether or not the sort spills.
fn write_fixture(path: &Path) {
    let header = create_minimal_header("chr1", 100_000);
    let records: Vec<_> = (0..40)
        .flat_map(|i| create_umi_family("ACGT", 2, &format!("fam_{i:06}"), "ACGTACGTAC", 35))
        .collect();
    write_bam(path, &header, &records);
}

/// Runs `fgumi sort` at `info` verbosity over a fresh fixture, threading the
/// optional `--pipeline-trace <flag>`, `--pipeline-trace-out <path>`, and
/// `FGUMI_PIPELINE_TRACE=<env>`. Returns the captured stderr; asserts success.
fn run_sort(flag: Option<&str>, env_trace: Option<&str>, trace_out: Option<&Path>) -> String {
    let tmp = TempDir::new().expect("tempdir");
    let input: PathBuf = tmp.path().join("unsorted.bam");
    write_fixture(&input);
    let output = tmp.path().join("sorted.bam");

    let mut cmd = Command::new(env!("CARGO_BIN_EXE_fgumi"));
    cmd.env("RUST_LOG", "info");
    cmd.env_remove("FGUMI_PIPELINE_TRACE");
    if let Some(value) = env_trace {
        cmd.env("FGUMI_PIPELINE_TRACE", value);
    }
    cmd.args(["sort", "-i"]).arg(&input).arg("-o").arg(&output).args(["--order", "coordinate"]);
    if let Some(level) = flag {
        cmd.args(["--pipeline-trace", level]);
    }
    if let Some(path) = trace_out {
        cmd.arg("--pipeline-trace-out").arg(path);
    }

    let result = cmd.output().expect("run fgumi sort");
    let stderr = String::from_utf8_lossy(&result.stderr).into_owned();
    assert!(result.status.success(), "fgumi sort failed:\n{stderr}");
    stderr
}

/// The `--pipeline-trace` flag and the `FGUMI_PIPELINE_TRACE` env var each
/// reach `PipelineConfig` and gate the end-of-run edge table; when both are set
/// the env value wins (in either direction). This is the honoring path that no
/// in-process test can reach.
#[rstest]
// flag alone → the level it names decides.
#[case::flag_deep(Some("deep"), None, true)]
#[case::default_off(None, None, false)]
// env alone → enables tracing without the flag (the standalone-enable path).
#[case::env_summary(None, Some("summary"), true)]
// env set → precedence over the flag, in both directions.
#[case::env_off_overrides_flag_deep(Some("deep"), Some("off"), false)]
#[case::env_deep_overrides_flag_off(Some("off"), Some("deep"), true)]
fn pipeline_trace_gates_the_edge_table(
    #[case] flag: Option<&str>,
    #[case] env_trace: Option<&str>,
    #[case] expect_table: bool,
) {
    let stderr = run_sort(flag, env_trace, None);
    assert_eq!(
        stderr.contains(EDGE_TABLE_SUBSTRING),
        expect_table,
        "expected edge-table presence={expect_table} for flag={flag:?} env={env_trace:?}; \
         stderr:\n{stderr}"
    );
}

/// The valid env-precedence cases above prove tracing turns *on or off*, but not
/// *which* level the env value resolves to: a resolver that mapped `summary` to
/// `timeline`/`deep`, or fell back to `deep` on an unrecognized value, would slip
/// through an on/off check. Pinning the resolved level needs a level-sensitive
/// probe — pair each env value with a conflicting `--pipeline-trace` flag (to
/// prove the env wins) *and* a `--pipeline-trace-out` path, then read the
/// TSV-write behavior, which differs by resolved level: `timeline`/`deep` write
/// the TSV, while `off`/`summary` skip it and warn the path is ignored.
#[rstest]
// env `summary` overrides CLI `deep`: resolves to `summary` (edge table on) but,
// being below `timeline`, skips the TSV and warns the path is ignored — so the
// env did *not* resolve to `deep`/`timeline`, which would have written it.
#[case::env_summary_overrides_flag_deep(Some("deep"), "summary", true, false, true)]
// env `deep` overrides CLI `off`: resolves to `deep` — edge table on and the TSV
// is written, so the env did *not* resolve to `off`/`summary`.
#[case::env_deep_overrides_flag_off(Some("off"), "deep", true, true, false)]
// unrecognized env overrides CLI `deep`: resolves to `off` (edge table gone)
// despite the flag; the path, now below `timeline`, is skipped with the warning.
#[case::env_invalid_overrides_flag_deep(Some("deep"), "summry", false, false, true)]
fn env_level_overrides_flag_and_resolves_specific_level(
    #[case] flag: Option<&str>,
    #[case] env_trace: &str,
    #[case] expect_table: bool,
    #[case] expect_tsv: bool,
    #[case] expect_path_ignored_warning: bool,
) {
    let tmp = TempDir::new().expect("tempdir");
    let trace_tsv = tmp.path().join("trace.tsv");
    let stderr = run_sort(flag, Some(env_trace), Some(&trace_tsv));

    assert_eq!(
        stderr.contains(EDGE_TABLE_SUBSTRING),
        expect_table,
        "expected edge-table presence={expect_table} for flag={flag:?} env={env_trace:?}; \
         stderr:\n{stderr}"
    );
    assert_eq!(
        trace_tsv.exists(),
        expect_tsv,
        "expected trace-out TSV presence={expect_tsv} for flag={flag:?} env={env_trace:?}; \
         stderr:\n{stderr}"
    );
    assert_eq!(
        stderr.contains("the path is ignored"),
        expect_path_ignored_warning,
        "expected path-ignored warning={expect_path_ignored_warning} for flag={flag:?} \
         env={env_trace:?}; stderr:\n{stderr}"
    );
}

/// `--pipeline-trace <level> --pipeline-trace-out <path>` writes the per-tick
/// timeline TSV to the requested path — proving `--pipeline-trace-out` is
/// threaded into `PipelineConfig::trace_path`, not just parsed. Both levels that
/// write the TSV are covered: `timeline` (the lowest such level, whose whole
/// purpose is the TSV) and `deep`. A defect where `timeline` fails to initialize
/// or write `trace_path` would otherwise slip through a `deep`-only test.
#[rstest]
#[case::timeline("timeline")]
#[case::deep("deep")]
fn pipeline_trace_out_writes_the_timeline_tsv(#[case] level: &str) {
    let tmp = TempDir::new().expect("tempdir");
    let trace_tsv = tmp.path().join("trace.tsv");
    let stderr = run_sort(Some(level), None, Some(&trace_tsv));

    assert!(stderr.contains(EDGE_TABLE_SUBSTRING), "{level} trace did not render the edge table");
    assert!(trace_tsv.exists(), "--pipeline-trace-out path was not written:\n{stderr}");
    let tsv = std::fs::read_to_string(&trace_tsv).expect("read trace tsv");
    let mut lines = tsv.lines();
    assert!(
        lines.next().is_some_and(|header| header.starts_with("t_ms")),
        "timeline TSV missing its `t_ms` header row:\n{tsv:.200}"
    );
    // The sampler writes a row per tick and a guaranteed final row for a run too
    // short to tick, so a `timeline`/`deep` TSV always carries at least one data
    // row — a header-only file means the per-tick writer never ran.
    assert!(
        lines.next().is_some(),
        "timeline TSV has only its `t_ms` header, no per-tick data row:\n{tsv:.400}"
    );
}

/// A `--pipeline-trace-out` path is only consulted at the levels that write the
/// per-tick TSV (`timeline`/`deep`). When the resolved level does *not* write one
/// (`off` or `summary`), `build_pipeline_config_for_chain` warns that the path is
/// ignored rather than dropping it silently, and no TSV is created. This pins the
/// warn-and-skip branch, complementing the `timeline`/`deep` honoring cases above;
/// without it, a regression that silently swallowed the path at these levels would
/// go unnoticed.
#[rstest]
#[case::off("off")]
#[case::summary("summary")]
fn pipeline_trace_out_ignored_below_timeline_warns(#[case] level: &str) {
    let tmp = TempDir::new().expect("tempdir");
    let trace_tsv = tmp.path().join("trace.tsv");
    let stderr = run_sort(Some(level), None, Some(&trace_tsv));

    assert!(
        stderr.contains("the path is ignored"),
        "expected a warning that the trace-out path is ignored at `{level}`:\n{stderr}"
    );
    assert!(
        !trace_tsv.exists(),
        "no timeline TSV must be written at `{level}`, but the path was created:\n{stderr}"
    );
}

/// An unrecognized `FGUMI_PIPELINE_TRACE` value falls back to `Off` (tracing
/// disabled) but, unlike the clap-validated flag, would otherwise do so
/// silently — so a warning naming the bad value is emitted.
#[test]
fn unrecognized_env_value_warns_and_disables_tracing() {
    let stderr = run_sort(None, Some("summry"), None);
    assert!(
        !stderr.contains(EDGE_TABLE_SUBSTRING),
        "an unrecognized FGUMI_PIPELINE_TRACE must not enable tracing:\n{stderr}"
    );
    assert!(
        stderr.contains("is not a recognized level"),
        "expected a warning naming the unrecognized FGUMI_PIPELINE_TRACE value:\n{stderr}"
    );
    assert!(
        stderr.contains("summry"),
        "the warning must name the offending value (`summry`), not just report a bad level:\n{stderr}"
    );
}

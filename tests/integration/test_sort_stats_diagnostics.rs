//! Integration tests for `fgumi sort --sort-stats`.
//!
//! `--sort-stats` is meant to gate the sort's performance diagnostics: on a
//! spilling sort, `SortMerge`'s "Sort merge diag:" line (merge-loop stalls /
//! contention / backpressure counters); on a sort that fits entirely in
//! memory, the "Sort fast-path diag:" note (there is no k-way merge to
//! diagnose there). Before this fix the chain-based `fgumi sort` ignored the
//! flag entirely on the spilling path (the diag line was emitted
//! unconditionally, and `--sort-stats` itself only produced a stale "is
//! ignored by the sort chain" warning) and produced nothing at all on the
//! fast path (the flag silently doing nothing there too). These tests pin
//! both paths to something real.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::{Mutex, Once, OnceLock};

use rstest::rstest;
use tempfile::TempDir;

use fgumi_lib::commands::common::{
    CompressionOptions, MaxTempFiles, MemoryLimit, MemoryReserve, QueueMemoryOptions,
    SchedulerOptions, ThreadingOptions,
};
use fgumi_lib::commands::group::GroupOptions;
use fgumi_lib::commands::sort::{SortOptions, SortOrderArg};
use fgumi_lib::pipeline::chains::{
    ChainSpec, SinkSpec, SourceSpec, Stage, StageOptionsBag, build_for,
};
use fgumi_lib::sam::SamTag;

use crate::helpers::assertions::{assert_bam_sorted, string_tag};
use crate::helpers::bam_generator::{
    create_minimal_header, create_umi_family, create_umi_family_at_pos, shuffled_umi_families,
    write_bam,
};

/// The stable substring pinning the `SortMerge` k-way-merge diagnostic line
/// (`crates/fgumi-pipeline-io/src/sort/merge.rs`, `emit_batches_cooperative`).
/// Only ever logged when a real merge ran (a spilling sort).
const MERGE_DIAG_SUBSTRING: &str = "Sort merge diag:";

/// The stable substring pinning the in-memory fast-path diagnostic note
/// (`crates/fgumi-pipeline-io/src/sort/merge.rs`, `emit_fast_batches`). Only
/// ever logged when the sort fit entirely in memory (no k-way merge ran).
/// Deliberately a different prefix from [`MERGE_DIAG_SUBSTRING`] so the two
/// are never mistaken for one another.
const FAST_PATH_DIAG_SUBSTRING: &str = "Sort fast-path diag:";

/// Shared buffer the capture logger appends every formatted `log` record to.
/// A single global buffer (not thread-local) so a diagnostic emitted on a sort
/// worker thread is still seen from the test thread.
static CAPTURED_LOGS: OnceLock<Mutex<Vec<String>>> = OnceLock::new();

/// A minimal `log::Log` that records each message into [`CAPTURED_LOGS`]. Used
/// by the in-process chain tests below, which run `build_for(..).run()`
/// directly and so cannot capture a subprocess's stderr the way the CLI tests
/// do.
struct CaptureLogger;

impl log::Log for CaptureLogger {
    fn enabled(&self, _: &log::Metadata<'_>) -> bool {
        true
    }

    fn log(&self, record: &log::Record<'_>) {
        if let Some(buf) = CAPTURED_LOGS.get() {
            buf.lock().expect("capture-log mutex poisoned").push(format!("{}", record.args()));
        }
    }

    fn flush(&self) {}
}

/// Installs [`CaptureLogger`] as the global logger (idempotent) and returns the
/// shared buffer, cleared for the caller. `log::set_boxed_logger` succeeds only
/// once per process; a later `Err` (already set) is fine because the buffer is
/// the same. Nextest runs each test in its own process, so no other test in
/// this binary competes for the global logger.
fn capture_logs() -> &'static Mutex<Vec<String>> {
    static INIT: Once = Once::new();
    let buf = CAPTURED_LOGS.get_or_init(|| Mutex::new(Vec::new()));
    INIT.call_once(|| {
        let _ = log::set_boxed_logger(Box::new(CaptureLogger));
        log::set_max_level(log::LevelFilter::Info);
    });
    buf.lock().expect("capture-log mutex poisoned").clear();
    buf
}

/// Writes a BAM with `families` three-read UMI families, all at the same
/// position.
fn write_bam_fixture(path: &Path, families: usize) {
    let header = create_minimal_header("chr1", 100_000);
    let records: Vec<_> = (0..families)
        .flat_map(|i| create_umi_family("ACGT", 3, &format!("fam_{i:06}"), "ACGTACGTAC", 35))
        .collect();
    write_bam(path, &header, &records);
}

/// Runs `fgumi sort` at info verbosity and returns its stderr.
fn sort_and_capture_logs(families: usize, max_memory: &str, extra_args: &[&str]) -> String {
    let tmp = TempDir::new().expect("tempdir");
    let input: PathBuf = tmp.path().join("unsorted.bam");
    write_bam_fixture(&input, families);
    let output = tmp.path().join("sorted.bam");

    let result = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .env("RUST_LOG", "info")
        .args(["sort", "-i"])
        .arg(&input)
        .arg("-o")
        .arg(&output)
        .args(["--order", "coordinate", "-m", max_memory])
        .args(extra_args)
        .output()
        .expect("run fgumi sort");

    let stderr = String::from_utf8_lossy(&result.stderr).into_owned();
    assert!(result.status.success(), "fgumi sort failed:\n{stderr}");
    stderr
}

/// The stable header pinning the restored per-phase timing roll-up
/// (`src/lib/pipeline/chains/commands/sort.rs`,
/// `SortPhaseTimingFinalizeHook`). Emitted only when a `PipelineStats`
/// collector is attached (`--pipeline-stats` / `FGUMI_PIPELINE_STATS=1`).
const PHASE_TIMING_HEADER: &str = "=== Sort Phase Timing ===";

/// Parse the percentage from a `=== Sort Phase Timing ===` phase line
/// (`  <label>   <secs>s   <pct>%`), or `None` if no line contains `label`.
/// The percentage is the last whitespace-separated token, e.g. `9.8%`.
fn phase_timing_pct(stderr: &str, label: &str) -> Option<f64> {
    stderr
        .lines()
        .find(|line| line.contains(label))
        .and_then(|line| line.split_whitespace().last())
        .and_then(|token| token.strip_suffix('%'))
        .and_then(|pct| pct.trim().parse::<f64>().ok())
}

/// Runs `fgumi sort` at info verbosity with `FGUMI_PIPELINE_STATS=1` set (which
/// force-attaches the pipeline stats collector even though standalone
/// `fgumi sort` does not flatten `--pipeline-stats`), and returns its stderr.
fn sort_and_capture_logs_with_pipeline_stats(
    families: usize,
    max_memory: &str,
    extra_args: &[&str],
) -> String {
    let tmp = TempDir::new().expect("tempdir");
    let input: PathBuf = tmp.path().join("unsorted.bam");
    write_bam_fixture(&input, families);
    let output = tmp.path().join("sorted.bam");

    let result = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .env("RUST_LOG", "info")
        .env("FGUMI_PIPELINE_STATS", "1")
        .args(["sort", "-i"])
        .arg(&input)
        .arg("-o")
        .arg(&output)
        .args(["--order", "coordinate", "-m", max_memory])
        .args(extra_args)
        .output()
        .expect("run fgumi sort");

    let stderr = String::from_utf8_lossy(&result.stderr).into_owned();
    assert!(result.status.success(), "fgumi sort failed:\n{stderr}");
    stderr
}

/// Runs `fgumi sort` with a small `--max-memory` (forcing spills + a real
/// k-way merge) and returns its stderr.
///
/// Sized so a small `--max-memory` forces multiple spill runs and therefore a
/// real k-way merge (not the single-chunk fast path, which never logs the
/// merge diagnostic at all): mirrors the fixture size
/// `test_streaming_output::sort_streams_to_stdout`'s `spill_and_merge` case
/// already validates as spilling (`test_streaming_output::the_spilling_fixture_spills`).
fn sort_spilling(extra_args: &[&str]) -> String {
    sort_and_capture_logs(20_000, "1M", extra_args)
}

/// Asserts the given stderr proves a real spill + k-way merge happened,
/// independent of `--sort-stats`: the standalone-sort summary's "Spill runs:"
/// line (`SortSummaryFinalizeHook`, unconditional) is only logged when
/// `runs_written > 0`, which is exactly the condition (`!slots.is_empty()`)
/// that rules out `SortMerge`'s single-chunk fast path. Guards the
/// merge-diagnostic assertions below against a memory-accounting change
/// silently making [`sort_spilling`]'s fixture fit in memory (which would
/// make "the diagnostic is present" vacuously true only because the fast-path
/// note happens to share no text with it, rather than because a real merge
/// diagnostic fired).
fn assert_really_spilled(stderr: &str) {
    assert!(
        stderr.contains("Spill runs:"),
        "fixture did not spill (no k-way merge ran) -- the merge-diagnostic assertions in this \
         test are vacuous without a real merge; stderr:\n{stderr}"
    );
}

/// Strip `env_logger`'s `[<time> <LEVEL> <module>] ` prefix: every line `fgumi`
/// logs through `log::info!` carries it (default `env_logger` format,
/// `src/main.rs`), so column parsing must start after the first `"] "`.
fn log_body(line: &str) -> &str {
    line.split_once("] ").map_or(line, |(_, body)| body)
}

/// Writes a BAM of `families` three-read UMI families at seeded pseudo-random
/// positions, in generation order, so the file is NOT coordinate-sorted. Every
/// tenth family reuses the previous family's position, so equal coordinate keys
/// exist across families. Names are hashed so name order is not index order either.
fn write_shuffled_bam_fixture(path: &Path, families: usize, seed: u64) {
    let header = create_minimal_header("chr1", 100_000);
    write_bam(path, &header, &shuffled_umi_families(families, seed));
}

/// Runs `fgumi sort --order <order>` on the shuffled fixture at info verbosity
/// (and `FGUMI_PIPELINE_STATS=1` when `with_stats`); returns stderr.
fn sort_shuffled(
    order: &str,
    families: usize,
    max_memory: &str,
    extra_args: &[&str],
    with_stats: bool,
) -> String {
    let tmp = TempDir::new().expect("tempdir");
    let input: PathBuf = tmp.path().join("unsorted.bam");
    write_shuffled_bam_fixture(&input, families, 0x5EED);
    let output = tmp.path().join("sorted.bam");
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_fgumi"));
    cmd.env("RUST_LOG", "info");
    if with_stats {
        cmd.env("FGUMI_PIPELINE_STATS", "1");
    }
    let result = cmd
        .args(["sort", "-i"])
        .arg(&input)
        .arg("-o")
        .arg(&output)
        .args(["--order", order, "-m", max_memory])
        .args(extra_args)
        .output()
        .expect("run fgumi sort");
    let stderr = String::from_utf8_lossy(&result.stderr).into_owned();
    assert!(result.status.success(), "fgumi sort failed:\n{stderr}");
    stderr
}

/// `fgumi sort --order <order>` on the shuffled multi-source fixture (unsorted
/// in every order, with equal-key ties) under a small `--max-memory`, at info
/// verbosity, with `extra` args; returns stderr. Callers guard on
/// `Merge sources: N` with `N >= 2`.
fn sort_spilling_with(order: &str, extra: &[&str]) -> String {
    sort_shuffled(order, 20_000, "1M", extra, false)
}

/// `N` from the standalone-sort summary line `Merge sources: N`.
fn merge_sources(stderr: &str) -> Option<u64> {
    stderr
        .lines()
        .map(log_body)
        .find_map(|b| b.trim().strip_prefix("Merge sources: "))
        .and_then(|n| n.trim().parse().ok())
}

/// Parse `(N, U, T)` from the first
/// ``detached `<step>`: ... ; N parks (unparked U, timed out T, avg ...)`` line.
fn detached_parks(stderr: &str, step: &str) -> Option<(u64, u64, u64)> {
    let needle = format!("detached `{step}`:");
    let line = stderr.lines().map(log_body).find(|l| l.contains(&needle))?;
    let after = line.split("; ").nth(1)?; // "N parks (unparked U, timed out T, avg ..."
    let n = after.split_whitespace().next()?.parse().ok()?;
    let field = |label: &str| -> Option<u64> {
        let rest = after.split_once(label)?.1;
        rest.split(|c: char| !c.is_ascii_digit()).next()?.parse().ok()
    };
    Some((n, field("(unparked ")?, field(", timed out ")?))
}

/// On a sort that merges several runs, the coordination driver's waits are
/// booked to the step that last made progress on it, not to the first step of
/// its group (`ReadBlocks`). Before the attribution fix every park of the
/// shared `sort-coord` driver showed up under `ReadBlocks` and its other steps
/// (`FindBoundariesAndSort`, `SpillGather`, `SortMerge`) all printed `0 parks`.
/// Which of them parks depends on scheduling (the merge may never wait on a
/// loaded host), so the check is on the group, not on one step.
#[test]
fn sort_stats_attribute_coord_driver_parks_to_the_progressing_step() {
    let stderr = sort_shuffled("coordinate", 20_000, "512K", &["--threads", "4"], true);
    let sources =
        merge_sources(&stderr).unwrap_or_else(|| panic!("no `Merge sources:` line:\n{stderr}"));
    assert!(sources > 1, "the fixture must merge more than one source (got {sources}):\n{stderr}");
    let parks = |step: &str| {
        detached_parks(&stderr, step)
            .unwrap_or_else(|| panic!("no `detached \\`{step}\\`` line:\n{stderr}"))
    };
    // Every park ended one of two ways, so each line's split must add up.
    for step in ["ReadBlocks", "FindBoundariesAndSort", "SpillGather", "SortMerge"] {
        let (n, unparked, timed_out) = parks(step);
        assert_eq!(
            unparked + timed_out,
            n,
            "`{step}`: unparked + timed out must equal its parks:\n{stderr}"
        );
    }
    let booked_past_the_primary: u64 =
        ["FindBoundariesAndSort", "SpillGather", "SortMerge"].iter().map(|s| parks(s).0).sum();
    assert!(
        booked_past_the_primary > 0,
        "coord-driver parks must be booked to the step that last progressed, not all to \
         ReadBlocks:\n{stderr}"
    );
}

/// Parse the integer in column `column` (by its header name, e.g. `gated_off`)
/// of the first `wake` table's row for `step`. The stats table is logged twice
/// (by `fgumi_pipeline_core::builder` and by `commands::common`); both copies
/// are identical, and the first is used.
fn wake_column(stderr: &str, step: &str, column: &str) -> Option<u64> {
    let bodies: Vec<&str> = stderr.lines().map(log_body).collect();
    let header_at = bodies
        .iter()
        .position(|b| b.trim_start().starts_with("wake ") && b.contains("gated_off"))?;
    let column = bodies[header_at].split_whitespace().position(|h| h == column)?;
    bodies[header_at + 1..]
        .iter()
        .take_while(|b| !b.trim_start().starts_with("worker"))
        .find(|b| b.split_whitespace().next() == Some(step))
        .and_then(|b| b.split_whitespace().nth(column))
        .and_then(|t| t.parse().ok())
}

/// On a multi-source sort the Detached `SpillWrite` producer is push-gated
/// (most of its `Progress` dispatches push nothing to the decompress edge), the
/// pool producers unpark the drivers, and `SortMerge` parks carry the
/// unparked/timed-out split. (Standalone `fgumi sort` arms the deadlock monitor,
/// which attaches the stats collector, so the stats tables — this one included
/// — are logged on every run; there is no stats-off CLI path to contrast with.)
#[test]
fn sort_stats_render_the_wake_table_under_pipeline_stats() {
    let stderr = sort_shuffled("coordinate", 20_000, "512K", &["--threads", "4"], true);
    let sources =
        merge_sources(&stderr).unwrap_or_else(|| panic!("no `Merge sources:` line:\n{stderr}"));
    assert!(sources > 1, "the fixture must merge more than one source (got {sources}):\n{stderr}");
    let has_table = stderr
        .lines()
        .map(log_body)
        .any(|b| b.trim_start().starts_with("wake ") && b.contains("gated_off"));
    assert!(has_table, "the wake table must be rendered:\n{stderr}");
    let gated_off = wake_column(&stderr, "SpillWrite", "gated_off").expect("SpillWrite gated_off");
    let notify = wake_column(&stderr, "SpillWrite", "notify").expect("SpillWrite notify");
    assert!(
        gated_off > 0 && notify < gated_off,
        "SpillWrite must be push-gated: notify={notify} gated_off={gated_off}\n{stderr}"
    );
    let unparks = wake_column(&stderr, "BgzfCompress", "unpark").expect("BgzfCompress unpark");
    assert!(unparks > 0, "BgzfCompress → WriteBgzfFile is a Driver target\n{stderr}");
    let (parks, unparked, timed_out) =
        detached_parks(&stderr, "SortMerge").expect("SortMerge park split");
    assert_eq!(unparked + timed_out, parks, "every SortMerge park ended one way:\n{stderr}");
    let group_parks: u64 = ["ReadBlocks", "FindBoundariesAndSort", "SpillGather", "SortMerge"]
        .iter()
        .map(|s| detached_parks(&stderr, s).map_or(0, |p| p.0))
        .sum();
    assert!(group_parks > 0, "the coordination driver parked at least once:\n{stderr}");
}

/// `--sort-stats` gates the `SortMerge` k-way-merge diagnostic line on a
/// sort that merges at least two sources: absent by default, present when
/// passed.
#[rstest]
#[case::without_flag(&[], false)]
#[case::with_flag(&["--sort-stats"], true)]
fn sort_stats_gates_merge_diagnostics(#[case] extra_args: &[&str], #[case] expect_present: bool) {
    let stderr = sort_spilling_with("coordinate", extra_args);
    assert_really_spilled(&stderr);
    let n = merge_sources(&stderr).unwrap_or_else(|| panic!("no `Merge sources:` line:\n{stderr}"));
    assert!(n >= 2, "the fixture must merge at least 2 sources (got {n}):\n{stderr}");

    let present = stderr.contains(MERGE_DIAG_SUBSTRING);
    assert_eq!(
        present, expect_present,
        "expected `{MERGE_DIAG_SUBSTRING}` presence={expect_present}, got:\n{stderr}"
    );
    // The flag must not resurrect the stale "ignored" warning.
    assert!(
        !stderr.contains("is ignored by the sort chain"),
        "stale --sort-stats warning still present:\n{stderr}"
    );
}

/// The merge-demand lines are gated on `--sort-stats`, on a sort that merges
/// at least two sources, in every order.
#[rstest]
fn sort_stats_gates_merge_demand_lines(
    #[values("coordinate", "queryname", "template-coordinate")] order: &str,
    #[values(false, true)] with_flag: bool,
) {
    let extra: &[&str] = if with_flag { &["--sort-stats"] } else { &[] };
    let stderr = sort_spilling_with(order, extra);
    let n = merge_sources(&stderr).unwrap_or_else(|| panic!("no `Merge sources:` line:\n{stderr}"));
    assert!(n >= 2, "{order}: the fixture must merge at least 2 sources (got {n}):\n{stderr}");
    for needle in [
        "Merge demand:",
        "Awaited slot at stall:",
        "Pool at stall:",
        "Consumer served itself:",
        "Merge output:",
    ] {
        assert_eq!(stderr.contains(needle), with_flag, "{order}: `{needle}` presence:\n{stderr}");
    }
}

/// `Byte fetch (input):` is gated on `--sort-stats` on the native input path
/// (the fixture is a regular BGZF file at the default `--read-streams auto`),
/// and its request histogram is populated, not merely rendered.
#[rstest]
fn sort_stats_gates_the_input_byte_fetch_line(#[values(false, true)] with_flag: bool) {
    let extra: &[&str] = if with_flag { &["--sort-stats"] } else { &[] };
    let stderr = sort_spilling_with("coordinate", extra);
    let line =
        stderr.lines().map(log_body).find_map(|l| l.trim().strip_prefix("Byte fetch (input): "));
    assert_eq!(line.is_some(), with_flag, "`Byte fetch (input):` presence:\n{stderr}");
    if let Some(rest) = line {
        let slices: u64 = rest.split_whitespace().next().and_then(|n| n.parse().ok()).unwrap();
        assert!(slices > 0, "the input histogram recorded no reads: {rest}");
    }
}

/// The merge supply's lines — the read-ahead budget, the spill byte fetch,
/// the measured spill block size, the disk reads, and the stash and claims —
/// are gated on `--sort-stats`, on a sort that merges at least two sources, in
/// every order.
#[rstest]
fn sort_stats_gates_spill_supply_lines(
    #[values("coordinate", "queryname", "template-coordinate")] order: &str,
    #[values(false, true)] with_flag: bool,
) {
    let extra: &[&str] = if with_flag { &["--sort-stats"] } else { &[] };
    let stderr = sort_spilling_with(order, extra);
    let n = merge_sources(&stderr).unwrap_or_else(|| panic!("no `Merge sources:` line:\n{stderr}"));
    assert!(n >= 2, "{order}: the fixture must merge at least 2 sources (got {n}):\n{stderr}");
    let bodies: Vec<&str> = stderr.lines().map(log_body).collect();
    for prefix in [
        "Spill supply: read-ahead budget R=",
        "Byte fetch: ",
        "Spill blocks: ",
        "Spill disk read: ",
        "Spill supply: frames ",
    ] {
        let present = bodies.iter().any(|b| b.starts_with(prefix));
        assert_eq!(present, with_flag, "{order}: `{prefix}` presence:\n{stderr}");
    }
}

/// On a multi-source merge, `--sort-stats` reports how often the merge served
/// itself from a stalled slot's stash. Whether a stall finds stashed blocks is
/// scheduling (the unit test `stalled_merge_self_serves_from_the_stash` forces
/// it), but the line must agree with itself and with the supply line whatever
/// the schedule: a park avoided served at least one block, and every block the
/// merge decompressed is a claim booked to the consumer.
#[rstest]
fn sort_stats_reports_self_serve(
    #[values("coordinate", "queryname", "template-coordinate")] order: &str,
) {
    let stderr = sort_spilling_with(order, &["--sort-stats"]);
    let n = merge_sources(&stderr).unwrap_or_else(|| panic!("no `Merge sources:` line:\n{stderr}"));
    assert!(n >= 2, "{order}: the fixture must merge at least 2 sources (got {n}):\n{stderr}");
    let body = |prefix: &str| -> &str {
        stderr
            .lines()
            .map(log_body)
            .find_map(|b| b.trim().strip_prefix(prefix))
            .unwrap_or_else(|| panic!("{order}: no `{prefix}` line:\n{stderr}"))
    };
    // "{E} parks avoided by decompressing already-read blocks inline ({B} blocks ..."
    let line = body("Consumer served itself: ");
    let episodes: u64 = line.split_whitespace().next().and_then(|n| n.parse().ok()).unwrap();
    let blocks: u64 =
        line.split_once('(').and_then(|(_, r)| r.split_whitespace().next()?.parse().ok()).unwrap();
    assert!(blocks >= episodes, "{order}: {line}");
    // "... claims {W} by workers + {C} by the consumer; ..."
    let supply = body("Spill supply: frames ");
    let consumer: u64 = supply
        .split_once(" by workers + ")
        .and_then(|(_, r)| r.split_whitespace().next()?.parse().ok())
        .unwrap();
    assert_eq!(consumer, blocks, "{order}: every merge-decompressed block is a consumer claim");
}

/// On a sort that fits entirely in memory (the single-chunk fast path, no
/// k-way merge), `--sort-stats` gates the fast-path note: absent by default,
/// present when the flag is passed. It must not be a silent no-op with the
/// flag, and — the disabled case pins that the note is truly gated, not
/// emitted unconditionally — it must stay silent without it. The merge
/// diagnostic never appears on this path regardless of the flag, because no
/// k-way merge runs.
#[rstest]
#[case::without_flag(&[], false)]
#[case::with_flag(&["--sort-stats"], true)]
fn sort_stats_gates_fast_path_note(#[case] extra_args: &[&str], #[case] expect_present: bool) {
    // A handful of records under a large `--max-memory`: fits in one
    // in-memory chunk, so `SortMerge` takes `FastPath` and never reaches
    // `emit_batches_cooperative` (where `MERGE_DIAG_SUBSTRING` is logged).
    let stderr = sort_and_capture_logs(4, "768M", extra_args);

    assert!(
        !stderr.contains("Spill runs:"),
        "fixture unexpectedly spilled; this case must exercise the fast path:\n{stderr}"
    );

    let present = stderr.contains(FAST_PATH_DIAG_SUBSTRING);
    assert_eq!(
        present, expect_present,
        "expected `{FAST_PATH_DIAG_SUBSTRING}` presence={expect_present} on a non-spilling sort, \
         got:\n{stderr}"
    );
    // No k-way merge runs on the fast path, so the merge diagnostic must never
    // appear regardless of the flag.
    assert!(
        !stderr.contains(MERGE_DIAG_SUBSTRING),
        "the k-way-merge diagnostic must not appear when no merge ran:\n{stderr}"
    );
}

// ─────────────────────────────────────────────────────────────────────────
// Restored `=== Sort Phase Timing ===` per-phase diagnostic
// ─────────────────────────────────────────────────────────────────────────
//
// The chain cutover dropped the owned engine's unconditional phase-timing
// block (the retired `RawExternalSorter::sort` was its only caller). It is
// restored as a roll-up derived from the pipeline's per-step timings, gated on
// the stats collector being attached so the default sort path is unaffected.

/// The `=== Sort Phase Timing ===` block is gated on the pipeline stats
/// collector being attached: emitted on a spilling sort under
/// `FGUMI_PIPELINE_STATS=1`, and absent on the default sort path so profiling
/// instrumentation stays opt-in and the default path keeps its zero-overhead
/// behavior. Mirrors this file's `sort_stats_gates_*` present/absent tables.
#[rstest]
#[case::enabled(true, true)]
#[case::disabled(false, false)]
fn sort_phase_timing_is_gated_on_pipeline_stats(
    #[case] with_pipeline_stats: bool,
    #[case] expect_present: bool,
) {
    let stderr = if with_pipeline_stats {
        sort_and_capture_logs_with_pipeline_stats(20_000, "1M", &[])
    } else {
        sort_spilling(&[])
    };
    assert_really_spilled(&stderr);

    assert_eq!(
        stderr.contains(PHASE_TIMING_HEADER),
        expect_present,
        "expected `{PHASE_TIMING_HEADER}` presence={expect_present} \
         (with_pipeline_stats={with_pipeline_stats}); stderr:\n{stderr}"
    );
    // When present, assert the phases that must dominate a spilling coordinate
    // sort report a NONZERO percentage — not merely that the label is printed.
    // The block emits every label unconditionally (a zero-time phase still
    // prints `0.000s 0.0%`), so a label-presence check would pass even if a
    // classifier bug routed a phase's steps to the wrong bucket; parsing the
    // percentage catches that. (Spill write is confirmed separately by
    // `assert_really_spilled`, and can round to 0.0% on a fast run, so it is not
    // asserted nonzero here.)
    if expect_present {
        for phase in ["in-memory sort", "k-way merge", "write output"] {
            let pct = phase_timing_pct(&stderr, phase).unwrap_or_else(|| {
                panic!("phase-timing block missing a parseable `{phase}` line; stderr:\n{stderr}")
            });
            assert!(
                pct > 0.0,
                "phase `{phase}` reported {pct}% on a spilling sort — a mis-bucketed classifier \
                 arm would zero a phase like this; stderr:\n{stderr}"
            );
        }
    }
}

/// `--sort-stats` alone (its own merge-diagnostic flag) must NOT trigger the
/// phase-timing block — that block is gated on the pipeline stats collector,
/// not on `--sort-stats`. Pins the two diagnostics as independently gated so a
/// future change cannot silently couple them.
#[test]
fn sort_stats_flag_alone_does_not_emit_phase_timing() {
    let stderr = sort_spilling(&["--sort-stats"]);
    assert_really_spilled(&stderr);
    assert!(
        stderr.contains(MERGE_DIAG_SUBSTRING),
        "the --sort-stats merge diagnostic should still fire; stderr:\n{stderr}"
    );
    assert!(
        !stderr.contains(PHASE_TIMING_HEADER),
        "--sort-stats alone must not emit the phase-timing block; stderr:\n{stderr}"
    );
}

// ─────────────────────────────────────────────────────────────────────────
// Coverage: the INTERMEDIATE `add_sort` branch's `.with_sort_stats(...)`
// ─────────────────────────────────────────────────────────────────────────
//
// The tests above exercise `--sort-stats` only through the standalone
// `fgumi sort` CLI, which always builds a `[Stage::Sort]`-terminal chain
// (`Sort::to_chain_spec`, `src/lib/commands/sort.rs`). That leaves the
// INTERMEDIATE branch of `add_sort` (`src/lib/pipeline/chains/builder.rs`,
// the arm reached when `position != StagePosition::Terminal` — the fused
// sort→group path) completely unexercised by any existing test: it is
// `SortMerge<RecordBatchOutput>::new(..).with_sort_stats(sort.sort_stats)`,
// followed by `DecodeFromRecords` and the `current_tail`/`chain_tail_kind`
// bookkeeping for a non-terminal sort.
//
// No CLI command builds a multi-stage, `Stage::Sort`-non-terminal chain yet
// (`runall` has not landed on this branch), but the chain builder itself has
// no such restriction, and `ChainSpec`/`build_for` are a deliberate,
// existing seam for exactly this: `validate.rs`'s Rule 6 comment notes "No
// CLI builds this chain today ... this protects a directly-constructed
// programmatic chain", and `test_chain_bam_with_index.rs` already
// establishes the pattern of building a `ChainSpec` directly and running it
// via `build_for(spec)?.run()` with no command-level driver involved. This
// test follows that same pattern for `[Stage::Sort, Stage::Group]` — the
// fused sort→group path the intermediate branch's own doc comment names.
#[test]
fn fused_sort_then_group_chain_exercises_intermediate_add_sort_branch() {
    let dir = TempDir::new().expect("create temp dir");
    let input_bam = dir.path().join("in.bam");
    let output_bam = dir.path().join("out.bam");

    // Three UMI families at distinct, deliberately out-of-order positions, so
    // the intermediate sort has real reordering work to do before group sees
    // template-coordinate-ordered input.
    let header = create_minimal_header("chr1", 100_000);
    let records: Vec<_> = create_umi_family_at_pos("AAAAAAAA", 3, "fam_a", "ACGTACGTAC", 35, 3_000)
        .into_iter()
        .chain(create_umi_family_at_pos("CCCCCCCC", 3, "fam_c", "ACGTACGTAC", 35, 1_000))
        .chain(create_umi_family_at_pos("GGGGGGGG", 3, "fam_g", "ACGTACGTAC", 35, 2_000))
        .collect();
    write_bam(&input_bam, &header, &records);

    // Capture the in-process chain's logs so we can prove `sort_stats: true`
    // actually reaches `SortMerge::with_sort_stats` on the intermediate branch
    // -- not just that the BAM came out sorted (which a regression dropping the
    // flag before `with_sort_stats` would leave unchanged). This fixture is
    // tiny and non-spilling, so the sort takes the single-chunk fast path and,
    // with sort-stats enabled, emits the fast-path diagnostic note.
    let logs = capture_logs();

    let sort_options = SortOptions {
        order: SortOrderArg::TemplateCoordinate,
        key_types: None,
        max_memory: MemoryLimit::Fixed(64 * 1024 * 1024),
        memory_reserve: MemoryReserve::Auto,
        memory_per_thread: true,
        tmp_dirs: Vec::new(),
        sort_threads: None,
        merge_threads: None,
        temp_compression: 1,
        temp_codec: fgumi_sort::SpillCodec::default(),
        max_temp_files: MaxTempFiles::Auto,
        // The point of this test: exercise `with_sort_stats(true)` on the
        // intermediate branch too, not just the flag being threaded through
        // without effect.
        sort_stats: true,
    };

    let spec = ChainSpec {
        stages: vec![Stage::Sort, Stage::Group],
        source: SourceSpec::Bam(input_bam.clone()),
        sink: SinkSpec::Bam(output_bam.clone()),
        stage_opts: StageOptionsBag {
            sort: Some(sort_options),
            group: Some(GroupOptions::default()),
            ..Default::default()
        },
        threading: ThreadingOptions { threads: None },
        compression: CompressionOptions::default(),
        scheduler: SchedulerOptions::default(),
        queue_memory: QueueMemoryOptions::default(),
        async_reader: false,
        read_streams: fgumi_bam_io::ReadStreams::Fixed(1),
        verify_crc: true,
        command_line: "fgumi <fused sort+group test>".to_string(),
    };

    build_for(spec)
        .expect("build_for should accept a fused [Stage::Sort, Stage::Group] chain")
        .run()
        .expect("running the fused sort->group chain should succeed with --sort-stats enabled");

    assert!(output_bam.exists(), "fused sort->group output BAM was not written");

    // The intermediate sort ran with `sort_stats: true` on a non-spilling
    // single-chunk fixture, so `SortMerge`'s fast-path note must have fired.
    // This is the assertion that pins propagation through `with_sort_stats`:
    // drop the flag before that call and this line fails while every BAM
    // assertion below still passes.
    let captured = logs.lock().expect("capture-log mutex poisoned");
    assert!(
        !captured.iter().any(|line| line.contains(MERGE_DIAG_SUBSTRING)),
        "the k-way-merge diagnostic must not appear on the non-spilling fast path; logs:\n{}",
        captured.join("\n")
    );
    assert!(
        captured.iter().any(|line| line.contains(FAST_PATH_DIAG_SUBSTRING)),
        "expected the fast-path note `{FAST_PATH_DIAG_SUBSTRING}` from the intermediate sort with \
         `sort_stats: true`; logs:\n{}",
        captured.join("\n")
    );
    drop(captured);

    // `Group`'s own runtime precondition (`require_group_input_ordering`) only
    // accepts template-coordinate-sorted input, so a correct, non-erroring
    // output here proves the intermediate sort really ran (`position !=
    // Terminal`) and its header transform + `SortMerge<RecordBatchOutput>` +
    // `DecodeFromRecords` handoff to group all worked end to end -- not just
    // that the call compiled.
    assert_bam_sorted(&output_bam, "template-coordinate", Some("mi"));

    let mut reader =
        noodles::bam::io::Reader::new(std::fs::File::open(&output_bam).expect("open output BAM"));
    reader.read_header().expect("read output BAM header");
    let mi_values: Vec<String> =
        reader.records().map(|r| string_tag(&r.expect("read output record"), SamTag::MI)).collect();
    assert_eq!(mi_values.len(), 9, "expected all 9 input records in the fused output");

    // Group MI values by their molecule id, stripping any trailing "/A" or
    // "/B" strand suffix (fgumi's own convention for pairing both strands of
    // a molecule under one id) so this assertion does not depend on whether
    // that suffix applies to single-end reads.
    let mut counts: HashMap<String, usize> = HashMap::new();
    for mi in &mi_values {
        let base = mi.split('/').next().unwrap_or(mi).to_string();
        *counts.entry(base).or_insert(0) += 1;
    }
    assert_eq!(
        counts.len(),
        3,
        "expected 3 distinct molecule groups (one per UMI family), got {counts:?}"
    );
    assert!(
        counts.values().all(|&n| n == 3),
        "expected each molecule group to have all 3 of its family's reads, got {counts:?}"
    );
}

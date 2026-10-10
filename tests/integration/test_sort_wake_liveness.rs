//! Liveness of sort chains under push-gated wakes: every side-channel reader in
//! the chain (merge demand, spill-run handoff, refill hint) and every gated
//! producer has a wake other than its park timer.
//!
//! Each case runs with every worker and driver on a 30 s park timer under a
//! 10 s watchdog, so a single wake that only a timer delivers fails the case. A
//! healthy run takes well under a second.

use std::ffi::OsStr;
use std::path::{Path, PathBuf};
use std::sync::{Mutex, Once, OnceLock};
use std::thread::ThreadId;
use std::time::Duration;

use fgumi_bam_io::ReadStreams;
use fgumi_lib::commands::common::MemoryLimit;
use fgumi_lib::commands::sort::SortOrderArg;
use fgumi_lib::pipeline::chains::{Stage, build_for};
use fgumi_pipeline_core::runtime::{TestBackoff, TestBackoffTarget};
use rstest::rstest;
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, write_bam};
use crate::helpers::chain_spec::{logged_count, sort_options, spec, uniquely_named_records};
use crate::helpers::cli::run_args_and_capture_logs;
use crate::helpers::cutover::decompressed_records_without_pg;
use crate::helpers::{read_bam_output, run_with_timeout};

/// Every worker's and driver's park timer in the slow run.
const PARK_TIMER_US: u64 = 30_000_000;
/// Shorter than one park timer, so one timer-only wake fails the case.
const WATCHDOG: Duration = Duration::from_secs(10);
/// Input records; with [`BUDGET`] a standalone sort spills well over 3 runs.
const N_RECORDS: usize = 30_000;
/// Total sort memory budget (not per thread).
const BUDGET: usize = 256 * 1024;

/// A spilling sort at t2, every thread on a 30 s park timer, completes inside a
/// 10 s watchdog, and its output is checked against the input independently:
/// `fgumi sort --verify` for the order and the read-name multiset (names are
/// unique) for completeness. The spill precondition (≥ 3 merge sources) is read
/// from the summary the same in-process run logs, so it vouches for the chain
/// under test. `Auto` is the CLI's default `--read-streams` and builds the
/// native pread front end (`PlanInputReads → PreadInputSlices →
/// FrameBgzfBlocks`); `Fixed(1)` builds the sequential `ReadBgzfBlocks` front.
#[rstest]
#[case::template_coordinate_auto(SortOrderArg::TemplateCoordinate, ReadStreams::Auto)]
#[case::template_coordinate_fixed1(SortOrderArg::TemplateCoordinate, ReadStreams::Fixed(1))]
#[case::coordinate_auto(SortOrderArg::Coordinate, ReadStreams::Auto)]
fn sort_chain_completes_under_long_timers_at_t2(
    #[case] order: SortOrderArg,
    #[case] read_streams: ReadStreams,
) {
    let dir = TempDir::new().unwrap();
    let input = dir.path().join("in.bam");
    write_bam(
        &input,
        &create_minimal_header("chr1", 1_000_000),
        &uniquely_named_records(N_RECORDS),
    );
    let output = dir.path().join("slow.bam");

    let logs = capture_logs();
    let (input2, output2) = (input.clone(), output.clone());
    let run_thread = run_with_timeout("wedged: a side channel has no wake", WATCHDOG, move || {
        run_chain_at_t2(vec![Stage::Sort], order, read_streams, &input2, &output2, true);
        std::thread::current().id()
    });

    let log = logs_of(logs, run_thread);
    assert!(
        logged_count(&log, "Merge sources:").unwrap_or(0) >= 3,
        "fixture must spill ≥ 3 runs on the chain under test:\n{log}"
    );
    let order_arg = match order {
        SortOrderArg::TemplateCoordinate => "template-coordinate",
        SortOrderArg::Coordinate => "coordinate",
        other => panic!("no CLI name mapped for {other:?}"),
    };
    run_args_and_capture_logs(&[
        OsStr::new("sort"),
        OsStr::new("--verify"),
        OsStr::new("-i"),
        output.as_os_str(),
        OsStr::new("--order"),
        OsStr::new(order_arg),
    ]);
    assert_eq!(read_names(&output), read_names(&input), "every input record, exactly once");
}

/// `[Sort, Group]` (template-coordinate, native input) at t2 under the same
/// timers and watchdog: the group stage's producers are gated too. Group
/// output has no independent oracle here, so it is checked against a
/// normal-timer run of the same chain and for the input's mapped read names
/// (group drops unmapped reads; every `RX` is unique, so nothing else merges
/// or drops). This chain logs no sort
/// summary, so the spill is not asserted; the budget is the one the standalone
/// case shows spilling.
#[test]
fn sort_group_chain_completes_under_long_timers_at_t2() {
    let dir = TempDir::new().unwrap();
    let input = dir.path().join("in.bam");
    write_bam(
        &input,
        &create_minimal_header("chr1", 1_000_000),
        &uniquely_named_records(N_RECORDS),
    );
    let normal = dir.path().join("normal.bam");
    let slow = dir.path().join("slow.bam");

    let stages = || vec![Stage::Sort, Stage::Group];
    let tc = SortOrderArg::TemplateCoordinate;
    let (i, o) = (input.clone(), normal.clone());
    run_with_timeout("normal-timer sort+group run", Duration::from_secs(60), move || {
        run_chain_at_t2(stages(), tc, ReadStreams::Auto, &i, &o, false);
    });
    let (i, o) = (input.clone(), slow.clone());
    run_with_timeout("wedged: a gated producer has no wake", WATCHDOG, move || {
        run_chain_at_t2(stages(), tc, ReadStreams::Auto, &i, &o, true);
    });

    let (slow_names, mapped_names) = (read_names(&slow), read_names_of_mapped(&input));
    assert!(
        slow_names == mapped_names,
        "group output must hold every mapped input record exactly once: {} records, want {}",
        slow_names.len(),
        mapped_names.len()
    );
    let (slow_bytes, normal_bytes) =
        (decompressed_records_without_pg(&slow), decompressed_records_without_pg(&normal));
    assert!(
        slow_bytes == normal_bytes,
        "slow-timer output {} ({} bytes) differs from normal-timer output {} ({} bytes)",
        slow.display(),
        slow_bytes.len(),
        normal.display(),
        normal_bytes.len()
    );
}

/// The chain `fgumi sort`/`runall` would build for `stages`, in-process (so
/// `test_backoff` can be set): 2 threads, a [`BUDGET`] total sort budget, the
/// given sort order and `--read-streams` policy; `slow` puts every worker and
/// driver on a [`PARK_TIMER_US`] park timer.
fn run_chain_at_t2(
    stages: Vec<Stage>,
    order: SortOrderArg,
    read_streams: ReadStreams,
    input: &Path,
    output: &Path,
    slow: bool,
) -> PathBuf {
    let mut sort = sort_options(None, None);
    sort.order = order;
    sort.max_memory = MemoryLimit::Fixed(BUDGET);
    sort.memory_per_thread = false;
    let mut chain = spec(stages, input, output, 2, sort);
    chain.read_streams = read_streams;
    let mut built = build_for(chain).expect("builds");
    if slow {
        built.config.test_backoff = vec![
            TestBackoff { target: TestBackoffTarget::AllWorkers, us: PARK_TIMER_US },
            TestBackoff { target: TestBackoffTarget::AllDrivers, us: PARK_TIMER_US },
        ];
    }
    built.run().expect("runs");
    output.to_path_buf()
}

/// The sorted read names of a BAM's records.
fn read_names(path: &Path) -> Vec<Vec<u8>> {
    names_where(path, |_| true)
}

/// The sorted read names of a BAM's mapped records.
fn read_names_of_mapped(path: &Path) -> Vec<Vec<u8>> {
    names_where(path, |r| !r.flags().is_unmapped())
}

fn names_where(
    path: &Path,
    keep: impl Fn(&noodles::sam::alignment::RecordBuf) -> bool,
) -> Vec<Vec<u8>> {
    let (_, records) = read_bam_output(path);
    let mut names: Vec<Vec<u8>> = records
        .iter()
        .filter(|r| keep(r))
        .map(|r| r.name().map(|n| n.to_vec()).unwrap_or_default())
        .collect();
    names.sort_unstable();
    names
}

/// Log lines captured in-process, tagged with the thread that logged them.
static CAPTURED: OnceLock<Mutex<Vec<(ThreadId, String)>>> = OnceLock::new();

struct CaptureLogger;

impl log::Log for CaptureLogger {
    fn enabled(&self, metadata: &log::Metadata<'_>) -> bool {
        metadata.level() <= log::Level::Info
    }

    fn log(&self, record: &log::Record<'_>) {
        if let Some(buf) = CAPTURED.get() {
            let line = (std::thread::current().id(), record.args().to_string());
            buf.lock().expect("capture mutex").push(line);
        }
    }

    fn flush(&self) {}
}

/// Installs [`CaptureLogger`] as the global logger (once per process; a later
/// install attempt is a no-op) and returns the shared buffer.
fn capture_logs() -> &'static Mutex<Vec<(ThreadId, String)>> {
    static INIT: Once = Once::new();
    let buf = CAPTURED.get_or_init(|| Mutex::new(Vec::new()));
    INIT.call_once(|| {
        let _ = log::set_boxed_logger(Box::new(CaptureLogger));
        log::set_max_level(log::LevelFilter::Info);
    });
    buf
}

/// The lines `thread` logged, joined. The sort summary is logged by the
/// thread that calls `run` (a success-gated finalize hook), so filtering by
/// it keeps another test's run in the same process out of the precondition.
fn logs_of(buf: &Mutex<Vec<(ThreadId, String)>>, thread: ThreadId) -> String {
    let buf = buf.lock().expect("capture mutex");
    buf.iter().filter(|(t, _)| *t == thread).map(|(_, l)| l.as_str()).collect::<Vec<_>>().join("\n")
}

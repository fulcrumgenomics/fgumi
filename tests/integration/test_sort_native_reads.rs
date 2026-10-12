//! The sort's chain-native input path (`PlanInputReads → PreadInputSlices →
//! FrameBgzfBlocks`): when it is selected, and that it changes nothing about
//! the output.

use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{Duration, Instant};

use rstest::rstest;
use tempfile::TempDir;

use crate::helpers::bam_generator::{
    create_minimal_header, shuffled_umi_families, transcode_bam_to_sam, write_bam,
};
use crate::helpers::read_bam_output;

/// `families` three-read UMI families at seeded pseudo-random positions with
/// hashed names, in generation order: unsorted in every order, and every tenth
/// family reuses the previous position, so equal keys exist.
fn write_shuffled_bam(path: &Path, families: usize) {
    write_bam(path, &create_minimal_header("chr1", 100_000), &shuffled_records(families));
}

/// [`write_shuffled_bam`]'s records in uncompressed (stored) BGZF blocks, so
/// a modest record count spans several 4 MiB read fills.
fn write_shuffled_stored_bam(path: &Path, families: usize) {
    let header = create_minimal_header("chr1", 100_000);
    let mut writer = fgumi_bam_io::create_raw_bam_writer(path, &header, 1, 0).unwrap();
    for record in shuffled_records(families) {
        writer.write_raw_record(record.as_ref()).unwrap();
    }
    writer.finish().unwrap();
}

/// The records of [`write_shuffled_bam`]: the shared shuffled fixture at seed
/// `0x5EED`.
fn shuffled_records(families: usize) -> Vec<fgumi_raw_bam::RawRecord> {
    shuffled_umi_families(families, 0x5EED)
}

/// The input shapes the native-path selection must distinguish.
#[derive(Clone, Copy, Debug)]
enum Input {
    /// A regular BGZF BAM.
    RegularBam,
    /// A BAM streamed through a FIFO.
    Fifo,
    /// A BAM under 64 KiB streamed through a FIFO: the writer finishes and
    /// closes before the reader is done, so a second open would block forever.
    SmallFifo,
    /// SAM text gzip-compressed (not BGZF) under a `.bam` name.
    PlainGzipBam,
    /// SAM text.
    Sam,
    /// A BAM on stdin (`-i -`).
    Stdin,
}

impl Input {
    /// Whether the input reaches the sort as BAM (and so through a raw-block
    /// source).
    fn is_bam(self) -> bool {
        !matches!(self, Input::Sam)
    }
}

/// Run `fgumi sort` on `input` with `FGUMI_PIPELINE_STATS=1` (so the stats
/// table names every step) under a 120 s watchdog; returns stderr.
fn run_sort_with_pipeline_stats(input: Input, extra: &[&str]) -> String {
    let tmp = TempDir::new().unwrap();
    let bam = tmp.path().join("src.bam");
    let families = if matches!(input, Input::SmallFifo) { 100 } else { 2_000 };
    write_shuffled_bam(&bam, families);
    let mut writer: Option<std::thread::JoinHandle<()>> = None;
    let path: PathBuf = match input {
        Input::RegularBam => bam.clone(),
        Input::Stdin => PathBuf::from("-"),
        Input::Sam => {
            let sam = tmp.path().join("in.sam");
            transcode_bam_to_sam(&bam, &sam);
            sam
        }
        Input::PlainGzipBam => {
            let sam = tmp.path().join("in.sam");
            transcode_bam_to_sam(&bam, &sam);
            let gz = tmp.path().join("in.bam");
            let mut enc = flate2::write::GzEncoder::new(
                std::fs::File::create(&gz).unwrap(),
                flate2::Compression::default(),
            );
            enc.write_all(&std::fs::read(&sam).unwrap()).unwrap();
            enc.finish().unwrap();
            gz
        }
        Input::Fifo | Input::SmallFifo => {
            let bytes = std::fs::read(&bam).unwrap();
            if matches!(input, Input::SmallFifo) {
                assert!(bytes.len() < 64 << 10, "small FIFO payload is {} bytes", bytes.len());
            }
            let fifo = tmp.path().join("in.bam");
            let status = Command::new("mkfifo").arg(&fifo).status().expect("run mkfifo");
            assert!(status.success());
            let f2 = fifo.clone();
            writer = Some(std::thread::spawn(move || {
                let mut w = std::fs::OpenOptions::new().write(true).open(&f2).unwrap();
                w.write_all(&bytes).unwrap();
            }));
            fifo
        }
    };
    let mut child = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .env("RUST_LOG", "info")
        .env("FGUMI_PIPELINE_STATS", "1")
        .args(["sort", "-i"])
        .arg(&path)
        .arg("-o")
        .arg(tmp.path().join("out.bam"))
        .args(["--order", "coordinate"])
        .args(extra)
        .stdin(if matches!(input, Input::Stdin) {
            std::process::Stdio::from(std::fs::File::open(&bam).unwrap())
        } else {
            std::process::Stdio::null()
        })
        .stderr(std::process::Stdio::piped())
        .spawn()
        .expect("run fgumi sort");
    let deadline = Instant::now() + Duration::from_secs(120);
    loop {
        if child.try_wait().unwrap().is_some() {
            break;
        }
        if Instant::now() > deadline {
            child.kill().unwrap();
            panic!("fgumi sort on {input:?} hung");
        }
        std::thread::sleep(Duration::from_millis(50));
    }
    let out = child.wait_with_output().unwrap();
    if let Some(w) = writer {
        w.join().unwrap();
    }
    let stderr = String::from_utf8_lossy(&out.stderr).into_owned();
    assert!(out.status.success(), "fgumi sort on {input:?} failed:\n{stderr}");
    stderr
}

/// The native path is taken only for a seekable BGZF file
/// with `--read-streams` other than `1`, and never by opening a FIFO to find
/// out (a small FIFO would hang).
#[rstest]
#[case::regular_bgzf_auto(Input::RegularBam, "auto", true)]
#[case::regular_bgzf_fixed1(Input::RegularBam, "1", false)]
#[case::fifo(Input::Fifo, "4", false)]
#[case::fifo_small_payload_auto(Input::SmallFifo, "auto", false)]
#[case::plain_gzip_named_bam(Input::PlainGzipBam, "auto", false)]
#[case::sam(Input::Sam, "auto", false)]
#[case::stdin(Input::Stdin, "auto", false)]
fn native_path_is_chosen_only_for_seekable_bgzf(
    #[case] input: Input,
    #[case] streams: &str,
    #[case] native: bool,
) {
    let stderr =
        run_sort_with_pipeline_stats(input, &["--read-streams", streams, "--threads", "4"]);
    assert_eq!(stderr.contains("PlanInputReads"), native, "{stderr}");
    assert_eq!(stderr.contains("ReadBgzfBlocks"), !native && input.is_bam(), "{stderr}");
}

/// An explicit `--read-streams N > 1` the input cannot honour says so.
#[rstest]
#[case::plain_gzip(Input::PlainGzipBam)]
#[case::sam(Input::Sam)]
#[case::fifo(Input::Fifo)]
#[case::stdin(Input::Stdin)]
fn explicit_streams_on_a_non_native_input_warn(#[case] input: Input) {
    let stderr = run_sort_with_pipeline_stats(input, &["--read-streams", "4", "--threads", "4"]);
    assert!(stderr.contains("--read-streams=4 applies only to seekable BGZF input"), "{stderr}");
}

/// A sorted output and the `Merge sources: N` it reported.
struct Sorted {
    records: Vec<noodles::sam::alignment::RecordBuf>,
    merge_sources: u64,
}

/// Sort the shuffled spilling fixture with `args`.
fn sort_to_records(order: &str, args: &[&str]) -> Sorted {
    let tmp = TempDir::new().unwrap();
    let input = tmp.path().join("in.bam");
    write_shuffled_bam(&input, 20_000);
    let output = tmp.path().join("out.bam");
    let out = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .env("RUST_LOG", "info")
        .args(["sort", "-i"])
        .arg(&input)
        .arg("-o")
        .arg(&output)
        .args(["--order", order])
        .args(args)
        .output()
        .expect("run fgumi sort");
    let stderr = String::from_utf8_lossy(&out.stderr).into_owned();
    assert!(out.status.success(), "fgumi sort failed:\n{stderr}");
    let merge_sources = stderr
        .lines()
        .find_map(|l| l.split_once("Merge sources: ").map(|(_, n)| n.trim().to_owned()))
        .and_then(|n| n.parse().ok())
        .unwrap_or(0);
    Sorted { records: read_bam_output(&output).1, merge_sources }
}

/// Output identity across read-streams values and thread counts, in every
/// order the arena front serves, on a fixture unsorted in all of them with
/// equal-key ties; both runs must merge ≥ 2 sources.
#[rstest]
fn native_reads_are_byte_identical(
    #[values(
        "template-coordinate",
        "coordinate",
        "queryname::natural",
        "queryname::lexicographic"
    )]
    order: &str,
    #[values("1", "auto", "4")] streams: &str,
    #[values("1", "4")] threads: &str,
) {
    // `--memory-per-thread false`: `-m` is the total at every thread count, so
    // both runs spill.
    let fixed = ["-m", "1M", "--memory-per-thread", "false"];
    let a =
        sort_to_records(order, &[&["--read-streams", "1", "--threads", "1"][..], &fixed].concat());
    let b = sort_to_records(
        order,
        &[&["--read-streams", streams, "--threads", threads][..], &fixed].concat(),
    );
    assert!(a.merge_sources >= 2 && b.merge_sources >= 2, "{order}: must merge ≥ 2 sources");
    assert_eq!(a.records.len(), 60_000);
    assert!(a.records == b.records, "{order} streams={streams} threads={threads}: records differ");
}

/// This process's OS thread count.
fn os_thread_count() -> usize {
    #[cfg(target_os = "linux")]
    {
        std::fs::read_dir("/proc/self/task").expect("procfs").count()
    }
    #[cfg(not(target_os = "linux"))]
    {
        // `ps -M` prints a header line and one line per thread.
        let out = Command::new("ps")
            .args(["-M", "-p", &std::process::id().to_string()])
            .output()
            .expect("run ps");
        String::from_utf8_lossy(&out.stdout).lines().count() - 1
    }
}

/// No thread is created for reads (a budget large enough that no run is sorted
/// before the input is read, so nothing else starts a thread in between). The
/// reference is the thread count at the planner's last fill, when the pipeline
/// has long finished starting: no fill may see more threads than that (a read
/// thread started at any fill, early ones included, raises it), and once a
/// fill sees the full count every later fill sees exactly it. Early fills may
/// see fewer while the pool and drivers are still being spawned.
#[cfg(unix)]
#[test]
fn native_reads_spawn_no_threads() {
    use clap::Parser as _;
    use fgumi_lib::commands::command::Command as _;
    use fgumi_pipeline_io::source::plan_input_reads::test_hooks;
    use std::sync::{Arc, Mutex};

    let tmp = TempDir::new().unwrap();
    let input = tmp.path().join("in.bam");
    // Stored blocks keep the fixture over several 4 MiB fills.
    write_shuffled_stored_bam(&input, 150_000);
    assert!(std::fs::metadata(&input).unwrap().len() > 24 << 20, "several fills");
    let counts = Arc::new(Mutex::new(Vec::new()));
    let c2 = Arc::clone(&counts);
    test_hooks::set(move |_last| c2.lock().unwrap().push(os_thread_count()));
    let sort = fgumi_lib::commands::sort::Sort::try_parse_from([
        "sort",
        "-i",
        input.to_str().unwrap(),
        "-o",
        tmp.path().join("out.bam").to_str().unwrap(),
        "--order",
        "coordinate",
        "--threads",
        "4",
        "--read-streams",
        "4",
        "-m",
        "2G",
    ])
    .unwrap();
    let result = sort.execute("fgumi sort test");
    test_hooks::clear();
    result.unwrap();
    let counts = counts.lock().unwrap();
    assert!(counts.len() >= 5, "several fills observed: {counts:?}");
    let full = *counts.last().expect("fills");
    assert!(counts.iter().all(|&c| c <= full), "a thread started while reading: {counts:?}");
    let started = counts.iter().position(|&c| c == full).expect("the last fill is full");
    assert!(
        counts[started..].iter().all(|&c| c == full),
        "the thread count moved after start-up: {counts:?}"
    );
    assert!(started + 3 <= counts.len(), "start-up spanned most of the read: {counts:?}");
}

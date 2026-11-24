//! Liveness of `runall`'s subprocess align stage when the aligner reads further
//! ahead than the in-flight budget before writing any output.
//!
//! The stand-in aligner here is an `awk` program that reads its whole stdin
//! before it prints anything (header included), then reports every read
//! unmapped. The hidden `FGUMI_ALIGN_IN_FLIGHT_BUDGET` override shrinks the
//! budget to 64 KiB, so a few MiB of input overflow it, and the hidden
//! `FGUMI_ALIGN_CEILING_TIMEOUT_SECS` shortens the wait before a run at the
//! ceiling fails. A gate that blocked on the full budget would wait forever for
//! output the aligner writes only at end of input, so each run is bounded by a
//! watchdog that kills it and fails the test instead of hanging CI.

use std::path::Path;
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

use fgumi_raw_bam::{RawRecord, SamBuilder, flags};
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, create_test_reference, write_bam};
use crate::helpers::cli::utf8_path as p;
use crate::helpers::read_bam_output;

/// A stand-in aligner that reads all of its interleaved FASTQ input before
/// writing anything, then reports every read unmapped (flag 77 for R1, 141 for
/// R2) under the `@SQ` of [`create_test_reference`].
const READ_ALL_THEN_EMIT_ALIGNER: &str = "awk 'NR % 4 == 1 { name[++n] = substr($1, 2) } \
    NR % 4 == 2 { seq[n] = $0 } NR % 4 == 0 { qual[n] = $0 } \
    END { OFS = \"\\t\"; print \"@HD\", \"VN:1.6\"; print \"@SQ\", \"SN:chr1\", \"LN:10000\"; \
    for (i = 1; i <= n; i++) print name[i], (i % 2 ? 77 : 141), \"*\", 0, 0, \"*\", \"*\", 0, 0, \
    seq[i], qual[i] }' # {ref}";

/// The in-flight budget the runs start from, far below the input's size.
const SMALL_BUDGET_BYTES: &str = "65536";

/// Template pairs in the input: about 16 MiB of unmapped records in memory,
/// fed to the aligner in several batches.
const PAIRS: usize = 50_000;

/// How long a run may take before it is declared hung: a run waits one or two
/// 5 s stall grace periods, plus the shortened ceiling timeout when it fails.
const RUN_TIMEOUT: Duration = Duration::from_secs(45);

/// Writes `PAIRS` unmapped 20 bp read pairs, named in queryname order.
fn write_unmapped_pairs(path: &Path) {
    let header = create_minimal_header("chr1", 10_000);
    let records: Vec<RawRecord> = (0..PAIRS)
        .flat_map(|i| {
            let name = format!("q{i:06}");
            [flags::FIRST_SEGMENT, flags::LAST_SEGMENT].map(|segment| {
                let mut b = SamBuilder::new();
                b.read_name(name.as_bytes())
                    .flags(flags::PAIRED | flags::UNMAPPED | flags::MATE_UNMAPPED | segment)
                    .sequence(b"ACGTACGTACGTACGTACGT")
                    .qualities(&[30; 20]);
                b.build()
            })
        })
        .collect();
    write_bam(path, &header, &records);
}

/// The outcome of a bounded `fgumi` run.
struct BoundedRun {
    success: bool,
    stderr: String,
}

/// Runs `runall --start-from align --stop-after zipper` over `input` with
/// [`READ_ALL_THEN_EMIT_ALIGNER`] and a [`SMALL_BUDGET_BYTES`] budget, killing
/// it and panicking if it has not exited within [`RUN_TIMEOUT`]. Its stderr goes
/// to a file, so a chatty run cannot block on a full pipe.
fn run_runall_bounded(dir: &Path, input: &Path, output: &Path, extra: &[&str]) -> BoundedRun {
    let reference = create_test_reference(dir);
    let stderr_path = dir.join("stderr.txt");
    let stderr_file = std::fs::File::create(&stderr_path).expect("create stderr file");
    let mut child = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args(["runall", "--start-from", "align", "--stop-after", "zipper"])
        .args(["-i", p(input), "--ref", p(&reference), "-o", p(output)])
        .args(["--aligner::command", READ_ALL_THEN_EMIT_ALIGNER])
        .args(extra)
        .env("FGUMI_ALIGN_IN_FLIGHT_BUDGET", SMALL_BUDGET_BYTES)
        .env("FGUMI_ALIGN_CEILING_TIMEOUT_SECS", "2")
        .env("RUST_LOG", "info")
        .stdin(Stdio::null())
        .stdout(Stdio::null())
        .stderr(stderr_file)
        .spawn()
        .expect("spawn fgumi runall");

    let deadline = Instant::now() + RUN_TIMEOUT;
    let status = loop {
        if let Some(status) = child.try_wait().expect("poll fgumi runall") {
            break status;
        }
        if Instant::now() >= deadline {
            let _ = child.kill();
            let _ = child.wait();
            let stderr = std::fs::read_to_string(&stderr_path).unwrap_or_default();
            panic!(
                "runall did not exit within {RUN_TIMEOUT:?}: the align stage is stalled on a \
                 full in-flight budget. stderr:\n{stderr}"
            );
        }
        std::thread::sleep(Duration::from_millis(50));
    };
    let stderr = std::fs::read_to_string(&stderr_path).expect("read stderr");
    BoundedRun { success: status.success(), stderr }
}

/// An aligner that reads all of its input before its first output byte
/// completes: the gate admits input past the full budget after the stall grace
/// period, warns once, and every read comes back paired with its own unmapped
/// half.
#[test]
fn aligner_reading_past_the_budget_before_output_completes() {
    let tmp = TempDir::new().expect("tempdir");
    let input = tmp.path().join("unmapped.bam");
    let output = tmp.path().join("out.bam");
    write_unmapped_pairs(&input);

    let run = run_runall_bounded(tmp.path(), &input, &output, &[]);
    assert!(run.success, "runall must complete:\n{}", run.stderr);
    for needle in ["raising it from 64.0 KiB", "--aligner::chunk-size"] {
        assert!(
            run.stderr.contains(needle),
            "the growth warning must be logged and point at --aligner::chunk-size \
             (missing {needle:?}):\n{}",
            run.stderr
        );
    }
    assert_eq!(
        run.stderr.matches("no aligner output was consumed for").count(),
        1,
        "the growth warning is logged once, however often the gate grows:\n{}",
        run.stderr
    );

    let (_, records) = read_bam_output(&output);
    assert_eq!(records.len(), 2 * PAIRS, "every read must be written");
    for (i, record) in records.iter().enumerate() {
        let expected = format!("q{:06}", i / 2);
        assert_eq!(record.name().map(ToString::to_string), Some(expected), "record {i}");
        assert!(record.flags().is_unmapped(), "record {i} must be unmapped");
        assert_eq!(record.flags().is_first_segment(), i % 2 == 0, "record {i} segment");
    }
}

/// Past the read-ahead ceiling the run fails with a diagnostic naming both
/// remedies, instead of hanging: the gate first grows to the 8 MiB ceiling,
/// then the next batch would pass it.
#[test]
fn aligner_reading_past_the_ceiling_fails_with_a_diagnostic() {
    let tmp = TempDir::new().expect("tempdir");
    let input = tmp.path().join("unmapped.bam");
    let output = tmp.path().join("out.bam");
    write_unmapped_pairs(&input);

    let run = run_runall_bounded(tmp.path(), &input, &output, &["--aligner::max-read-ahead", "8"]);
    assert!(!run.success, "runall must fail past the ceiling:\n{}", run.stderr);
    for needle in [
        "raising it from 64.0 KiB",
        "read ahead more than 8.0 MiB",
        "--aligner::chunk-size",
        "--aligner::max-read-ahead",
    ] {
        assert!(run.stderr.contains(needle), "stderr must contain {needle:?}:\n{}", run.stderr);
    }
}

//! `runall --metrics` writes a per-stage TSV file.
//!
//! Exercises the `--metrics` flag end-to-end via the `fgumi` binary using the
//! small grouped-BAM fixture (`--start-from group`), which does not require
//! external tools (bwa, samtools) and runs in well under a second.
//!
//! The tests verify:
//!
//! - The TSV is created when the flag is set, with the exact advertised
//!   header, at least one data row, and parseable typed columns.
//! - Omitting the flag does not create a stray file and does not fail.
//! - Pointing the flag at an unwritable path surfaces an error to the caller
//!   rather than silently succeeding.

use std::io::BufRead;
use std::path::Path;
use std::process::Command;

use crate::helpers::grouped_bam_fixture::build_small_grouped_bam;

fn fgumi_bin() -> &'static str {
    env!("CARGO_BIN_EXE_fgumi")
}

use fgumi_lib::runall::engine::stats::METRICS_TSV_HEADER as EXPECTED_HEADER;

fn base_args(input: &Path, output: &Path) -> Vec<std::ffi::OsString> {
    vec![
        "runall".into(),
        "--start-from".into(),
        "group".into(),
        "--threads".into(),
        "1".into(),
        "--input".into(),
        input.into(),
        "--output".into(),
        output.into(),
        // runall's default --stop-after is filter, which now requires min-reads
        // to match standalone `fgumi filter` semantics.
        "--filter::min-reads".into(),
        "1".into(),
    ]
}

fn run_fgumi(args: &[&std::ffi::OsStr]) -> std::process::Output {
    Command::new(fgumi_bin()).args(args).output().expect("spawn fgumi runall")
}

#[test]
fn runall_writes_metrics_tsv_when_flag_set() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input.bam");
    let output = tmp.path().join("out.bam");
    let metrics = tmp.path().join("metrics.tsv");
    build_small_grouped_bam(&input);

    let mut args = base_args(&input, &output);
    args.push("--metrics".into());
    args.push(metrics.clone().into());
    let arg_refs: Vec<&std::ffi::OsStr> = args.iter().map(AsRef::as_ref).collect();

    let out = run_fgumi(&arg_refs);
    assert!(
        out.status.success(),
        "fgumi runall failed: status={:?}\nstdout:\n{}\nstderr:\n{}",
        out.status,
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr),
    );

    assert!(metrics.exists(), "metrics TSV was not created at {}", metrics.display());

    let file = std::fs::File::open(&metrics).expect("open metrics TSV");
    let lines: Vec<String> =
        std::io::BufReader::new(file).lines().map(std::result::Result::unwrap).collect();

    // 1. Header matches exactly.
    assert_eq!(lines[0], EXPECTED_HEADER, "header mismatch");

    // 2. At least one data row.
    assert!(lines.len() >= 2, "expected at least one data row, got {} lines", lines.len());

    // 3-5. Every data row has 4 tab-separated fields; records_in / records_out
    //      parse as u64; wall_time_secs parses as f64.
    for (i, line) in lines.iter().enumerate().skip(1) {
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields.len(), 4, "row {i} has {} fields, expected 4: {line:?}", fields.len());
        let _: f64 = fields[1]
            .parse()
            .unwrap_or_else(|e| panic!("row {i} wall_time_secs not a float: {e} (line={line:?})"));
        let _: u64 = fields[2]
            .parse()
            .unwrap_or_else(|e| panic!("row {i} records_in not an integer: {e} (line={line:?})"));
        let _: u64 = fields[3]
            .parse()
            .unwrap_or_else(|e| panic!("row {i} records_out not an integer: {e} (line={line:?})"));
    }
}

#[test]
fn runall_without_metrics_flag_does_not_create_file_or_fail() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input.bam");
    let output = tmp.path().join("out.bam");
    let would_be_metrics = tmp.path().join("metrics.tsv");
    build_small_grouped_bam(&input);

    let args = base_args(&input, &output);
    let arg_refs: Vec<&std::ffi::OsStr> = args.iter().map(AsRef::as_ref).collect();

    let out = run_fgumi(&arg_refs);
    assert!(
        out.status.success(),
        "fgumi runall without --metrics should succeed: status={:?}\nstderr:\n{}",
        out.status,
        String::from_utf8_lossy(&out.stderr),
    );
    assert!(output.exists(), "expected output BAM at {}", output.display());
    assert!(
        !would_be_metrics.exists(),
        "metrics TSV was written to {} without --metrics flag",
        would_be_metrics.display()
    );
}

#[test]
fn runall_with_unwritable_metrics_path_surfaces_error() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input.bam");
    let output = tmp.path().join("out.bam");
    // A path under a non-existent parent directory -> File::create will fail.
    let bad_metrics = tmp.path().join("nonexistent-subdir").join("metrics.tsv");
    build_small_grouped_bam(&input);

    let mut args = base_args(&input, &output);
    args.push("--metrics".into());
    args.push(bad_metrics.clone().into());
    let arg_refs: Vec<&std::ffi::OsStr> = args.iter().map(AsRef::as_ref).collect();

    let out = run_fgumi(&arg_refs);
    assert!(
        !out.status.success(),
        "fgumi runall should fail when --metrics path is unwritable but succeeded: \
         stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr),
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("metrics") || stderr.contains("TSV") || stderr.contains("tsv"),
        "expected error message to reference metrics TSV, got stderr:\n{stderr}"
    );
}

/// `--filter::stats <path>` writes a 4-row TSV (`total_reads`, `passed_reads`,
/// `failed_reads`, `pass_rate`) matching the schema of standalone
/// `fgumi filter --stats`. Guards against runall silently dropping filter
/// metrics that users rely on.
#[test]
fn runall_writes_filter_stats_tsv_when_flag_set() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input.bam");
    let output = tmp.path().join("out.bam");
    let stats = tmp.path().join("filter_stats.tsv");
    build_small_grouped_bam(&input);

    let mut args = base_args(&input, &output);
    args.push("--filter::stats".into());
    args.push(stats.clone().into());
    let arg_refs: Vec<&std::ffi::OsStr> = args.iter().map(AsRef::as_ref).collect();

    let out = run_fgumi(&arg_refs);
    assert!(
        out.status.success(),
        "fgumi runall failed: stderr:\n{}",
        String::from_utf8_lossy(&out.stderr),
    );

    assert!(stats.exists(), "filter stats TSV was not created at {}", stats.display());
    let contents = std::fs::read_to_string(&stats).expect("read filter stats TSV");
    let lines: Vec<&str> = contents.lines().collect();
    assert_eq!(lines.len(), 4, "expected 4 rows, got:\n{contents}");

    // Schema: `<key>\t<value>` in a known order. Matches
    // `commands::filter::write_filter_stats`.
    let parts: Vec<Vec<&str>> = lines.iter().map(|l| l.split('\t').collect()).collect();
    assert_eq!(parts[0][0], "total_reads");
    assert_eq!(parts[1][0], "passed_reads");
    assert_eq!(parts[2][0], "failed_reads");
    assert_eq!(parts[3][0], "pass_rate");

    // Values parse.
    let _: u64 = parts[0][1].parse().expect("total_reads numeric");
    let _: u64 = parts[1][1].parse().expect("passed_reads numeric");
    let _: u64 = parts[2][1].parse().expect("failed_reads numeric");
    let _: f64 = parts[3][1].parse().expect("pass_rate numeric");
}

//! `sort --read-buffer`: a piped input is read ahead into a buffer taken out of
//! `--max-memory`, which must never change the output.
//!
//! The buffer changes how much of the budget is left for the sort (so where its
//! runs end and how many it spills), never which records come out or in what
//! order. The command line lands in the output's `@PG` header, so raw output
//! bytes (BAM and `.bai`) can only be compared between runs with the same
//! arguments: stdin from a pipe (buffered) against stdin redirected from the
//! file (a regular file, so not buffered). Runs with different `--read-buffer`
//! values are compared on everything but their `@PG` lines.

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

use rstest::rstest;
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, transcode_bam_to_sam, write_bam};
use crate::helpers::chain_spec::records;
use crate::helpers::cutover::decompressed_records_without_pg;

/// Records in the fixture: ~19 MB of records, so the 4 MiB budget below spills
/// several runs.
const RECORDS: usize = 200_000;

/// A fixed total budget small enough to spill several runs from the fixture,
/// and large enough that a quarter of it (the read-buffer cap, 1 MiB) holds a
/// few chunks.
const MAX_MEMORY: &str = "4MiB";

/// Write the unsorted fixture: records that tie on position and name, each
/// with a unique `RX`, so a change in the order of tied records shows.
fn write_fixture(dir: &Path) -> PathBuf {
    write_fixture_of(dir, RECORDS)
}

/// [`write_fixture`] with `n` records.
fn write_fixture_of(dir: &Path, n: usize) -> PathBuf {
    let input = dir.join("unsorted.bam");
    let records = records(n, |i| format!("q{}", (i * 31) % (n / 4)));
    write_bam(&input, &create_minimal_header("chr1", 1_000_000), &records);
    input
}

/// How the sort's input reaches it.
#[derive(Clone, Copy, Debug)]
enum Feed {
    /// `-i -`, fed through a pipe by a thread of this test.
    Stdin,
    /// `-i - < <fixture>`: stdin is the regular file itself.
    StdinFromFile,
    /// `-i <fixture>`: a regular file.
    File,
}

/// Run `fgumi sort -i - / <file> -o <output> --order coordinate --write-index
/// --read-buffer <read_buffer>` under a fixed 4 MiB budget, with `stdin_bytes`
/// (the fixture's bytes, or a prefix of them) piped in for [`Feed::Stdin`] and
/// the fixture file itself as stdin for [`Feed::StdinFromFile`].
fn run_sort(
    feed: Feed,
    input: &Path,
    stdin_bytes: Option<Vec<u8>>,
    output: &Path,
    read_buffer: &str,
) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_fgumi"));
    command.env("RUST_LOG", "info").arg("sort").arg("-i");
    match feed {
        Feed::Stdin | Feed::StdinFromFile => command.arg("-"),
        Feed::File => command.arg(input),
    };
    command
        .arg("-o")
        .arg(output)
        .args(["--order", "coordinate", "--write-index", "--threads", "2"])
        .args(["--max-memory", MAX_MEMORY, "--memory-per-thread", "false"])
        .args(["--read-buffer", read_buffer])
        .stdin(match feed {
            Feed::Stdin => Stdio::piped(),
            Feed::StdinFromFile => Stdio::from(fs::File::open(input).expect("open fixture")),
            Feed::File => Stdio::null(),
        })
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());
    let mut child = command.spawn().expect("spawn fgumi sort");
    // Write stdin from a thread: the sort reads it while this thread waits.
    let writer = stdin_bytes.map(|bytes| {
        let mut stdin = child.stdin.take().expect("stdin was piped");
        std::thread::spawn(move || {
            // The sort may exit early (e.g. on a truncated input) and close
            // the pipe; that is not this writer's failure to report.
            let _ = stdin.write_all(&bytes);
        })
    });
    let result = child.wait_with_output().expect("wait for fgumi sort");
    if let Some(writer) = writer {
        writer.join().expect("stdin writer thread");
    }
    result
}

/// Run [`run_sort`] and assert it succeeded; returns its stderr.
fn sort_ok(feed: Feed, input: &Path, output: &Path, read_buffer: &str) -> String {
    let stdin_bytes = matches!(feed, Feed::Stdin).then(|| fs::read(input).expect("read fixture"));
    let result = run_sort(feed, input, stdin_bytes, output, read_buffer);
    let stderr = String::from_utf8_lossy(&result.stderr).into_owned();
    assert!(result.status.success(), "fgumi sort --read-buffer {read_buffer} failed:\n{stderr}");
    stderr
}

/// The line the read buffer logs when it is closed, naming how much it read.
const BUFFER_SUMMARY: &str = "Read buffer: read ";

/// The line the startup log prints when the input is read without a buffer.
const NOT_A_STREAM: &str = "Read buffer: off (the input is not a pipe, FIFO or other stream)";

/// Assert the input was (`used`) or was not read through the read buffer, by
/// the summary line the buffer logs when it is closed.
fn assert_buffer_used(stderr: &str, used: bool) {
    assert_eq!(
        stderr.contains(BUFFER_SUMMARY),
        used,
        "expected the read buffer to be {}used:\n{stderr}",
        if used { "" } else { "un" }
    );
}

/// The BAM and `.bai` bytes written to `output`.
fn output_bytes(output: &Path) -> (Vec<u8>, Vec<u8>) {
    let bai = PathBuf::from(format!("{}.bai", output.display()));
    (fs::read(output).expect("read BAM"), fs::read(bai).expect("read BAI"))
}

/// The same command, once reading a pipe through the default read buffer and
/// once reading stdin redirected from the file (not buffered, so the sort keeps
/// the whole budget and its runs end elsewhere), writes byte-identical BAM and
/// index files while spilling.
#[test]
fn buffered_and_unbuffered_stdin_write_identical_bytes() {
    let dir = TempDir::new().expect("tempdir");
    let input = write_fixture(dir.path());
    let output = dir.path().join("sorted.bam");

    let unbuffered = sort_ok(Feed::StdinFromFile, &input, &output, "2048");
    assert!(unbuffered.contains(NOT_A_STREAM), "{unbuffered}");
    assert_buffer_used(&unbuffered, false);
    let unbuffered_bytes = output_bytes(&output);

    let buffered = sort_ok(Feed::Stdin, &input, &output, "2048");
    assert!(buffered.contains("Read buffer: 1.0 MiB"), "{buffered}");
    assert_buffer_used(&buffered, true);
    assert!(spill_runs(&unbuffered) > 1, "the fixture must spill:\n{unbuffered}");
    // The buffer comes out of the sort's budget, so the buffered sort's runs
    // are smaller and more of them spill.
    assert!(
        spill_runs(&buffered) > spill_runs(&unbuffered),
        "the read buffer was not taken out of the sort's budget:\n{buffered}\n{unbuffered}"
    );
    let buffered_bytes = output_bytes(&output);

    assert!(buffered_bytes.0 == unbuffered_bytes.0, "BAM bytes differ with the read buffer");
    assert!(buffered_bytes.1 == unbuffered_bytes.1, "BAI bytes differ with the read buffer");
}

/// No buffer, a buffer below the cap and a capped buffer sort a pipe to the
/// same records (everything but the `@PG` line, which carries the argument).
#[test]
fn read_buffer_size_does_not_change_the_records() {
    let dir = TempDir::new().expect("tempdir");
    let input = write_fixture(dir.path());
    let outputs: Vec<_> = ["0", "128K", "2048"]
        .into_iter()
        .map(|read_buffer| {
            let output = dir.path().join(format!("sorted-{read_buffer}.bam"));
            let stderr = sort_ok(Feed::Stdin, &input, &output, read_buffer);
            assert_buffer_used(&stderr, read_buffer != "0");
            (read_buffer, decompressed_records_without_pg(&output))
        })
        .collect();
    for (read_buffer, records) in &outputs[1..] {
        assert!(*records == outputs[0].1, "records differ with --read-buffer {read_buffer}");
    }
}

/// The number of spill runs the sort reported.
fn spill_runs(stderr: &str) -> u64 {
    crate::helpers::chain_spec::logged_count(stderr, "Spill runs:").unwrap_or(0)
}

/// The same records come out of a piped input as out of the file itself.
#[test]
fn piped_input_with_a_read_buffer_matches_the_file_input() {
    let dir = TempDir::new().expect("tempdir");
    let input = write_fixture(dir.path());
    let from_file = dir.path().join("from_file.bam");
    let from_pipe = dir.path().join("from_pipe.bam");
    assert_buffer_used(&sort_ok(Feed::File, &input, &from_file, "2048"), false);
    assert_buffer_used(&sort_ok(Feed::Stdin, &input, &from_pipe, "2048"), true);
    assert!(
        decompressed_records_without_pg(&from_file) == decompressed_records_without_pg(&from_pipe),
        "piped and file input sorted to different records"
    );
}

/// SAM text piped in through the default read buffer sorts to the same records
/// as the BAM file it was transcoded from. A quarter of the fixture keeps the
/// test quick: SAM text is slow to write and parse in a debug build.
#[test]
fn piped_sam_with_a_read_buffer_matches_the_file_input() {
    let dir = TempDir::new().expect("tempdir");
    let input = write_fixture_of(dir.path(), RECORDS / 4);
    let sam = dir.path().join("unsorted.sam");
    transcode_bam_to_sam(&input, &sam);
    let from_file = dir.path().join("from_file.bam");
    let from_pipe = dir.path().join("from_sam_pipe.bam");
    sort_ok(Feed::File, &input, &from_file, "2048");
    let stderr = sort_ok(Feed::Stdin, &sam, &from_pipe, "2GiB");
    assert!(stderr.contains("Read buffer: 1.0 MiB"), "{stderr}");
    assert_buffer_used(&stderr, true);
    assert!(
        decompressed_records_without_pg(&from_file) == decompressed_records_without_pg(&from_pipe),
        "piped SAM and file input sorted to different records"
    );
}

/// The startup line states the buffer in force: off when asked, off for a
/// regular file (which keeps its whole budget), and capped at a quarter of the
/// budget for a pipe; the input is read through the buffer only when it is on.
#[rstest]
#[case::off_for_a_pipe(Feed::Stdin, "0000", "Read buffer: off (--read-buffer 0)", false)]
#[case::off_for_a_file(Feed::File, "2048", NOT_A_STREAM, false)]
#[case::capped_for_a_pipe(
    Feed::Stdin,
    "2048",
    "capped at a quarter of the 4.0 MiB memory budget",
    true
)]
#[case::fits_for_a_pipe(
    Feed::Stdin,
    "128K",
    "Read buffer: 125.0 KiB (taken from the 4.0 MiB",
    true
)]
#[case::off_for_redirected_stdin(Feed::StdinFromFile, "2048", NOT_A_STREAM, false)]
fn startup_log_names_the_read_buffer_in_force(
    #[case] feed: Feed,
    #[case] read_buffer: &str,
    #[case] expected: &str,
    #[case] buffered: bool,
) {
    let dir = TempDir::new().expect("tempdir");
    let input = write_fixture(dir.path());
    let stderr = sort_ok(feed, &input, &dir.path().join("sorted.bam"), read_buffer);
    assert!(stderr.contains(expected), "expected {expected:?} in:\n{stderr}");
    assert_buffer_used(&stderr, buffered);
}

/// A pipe that ends mid-stream must fail the sort with the same error as it
/// does without the buffer: the buffer forwards end of input, it does not
/// paper over it.
#[test]
fn a_truncated_pipe_fails_the_sort() {
    let dir = TempDir::new().expect("tempdir");
    let input = write_fixture(dir.path());
    let mut bytes = fs::read(&input).expect("read fixture");
    bytes.truncate(bytes.len() * 2 / 3);
    let errors: Vec<Vec<String>> = ["0000", "2048"]
        .into_iter()
        .map(|read_buffer| {
            let output = dir.path().join(format!("sorted-{read_buffer}.bam"));
            let result = run_sort(Feed::Stdin, &input, Some(bytes.clone()), &output, read_buffer);
            let stderr = String::from_utf8_lossy(&result.stderr);
            assert!(!result.status.success(), "a truncated input must fail the sort:\n{stderr}");
            assert!(!output.exists(), "a failed sort must not leave an output");
            let errors: Vec<String> =
                stderr.lines().filter(|l| l.starts_with("Error:")).map(str::to_owned).collect();
            assert!(
                !errors.is_empty(),
                "no Error: line with --read-buffer {read_buffer}:\n{stderr}"
            );
            errors
        })
        .collect();
    assert_eq!(errors[0], errors[1], "the buffer changed the error");
    assert!(
        errors[0].iter().any(|line| line.contains("truncated BGZF block")),
        "the error does not name the truncation: {errors:?}"
    );
}

/// A FIFO — even one named `*.bam` — is a stream: it is read through the
/// buffer and sorts to the same records as the file.
#[test]
fn a_fifo_input_is_buffered() {
    let dir = TempDir::new().expect("tempdir");
    let input = write_fixture(dir.path());
    let fifo = dir.path().join("fifo.bam");
    let status = Command::new("mkfifo").arg(&fifo).status().expect("run mkfifo");
    assert!(status.success(), "mkfifo failed");

    let from_file = dir.path().join("from_file.bam");
    assert_buffer_used(&sort_ok(Feed::File, &input, &from_file, "2048"), false);

    // The writer blocks opening the FIFO until the sort opens it to read.
    let writer = {
        let (fifo, bytes) = (fifo.clone(), fs::read(&input).expect("read fixture"));
        std::thread::spawn(move || {
            let mut out = fs::OpenOptions::new().write(true).open(&fifo).expect("open FIFO");
            let _ = out.write_all(&bytes);
        })
    };
    let from_fifo = dir.path().join("from_fifo.bam");
    let stderr = sort_ok(Feed::File, &fifo, &from_fifo, "2048");
    writer.join().expect("FIFO writer thread");
    assert!(stderr.contains("Read buffer: 1.0 MiB"), "{stderr}");
    assert_buffer_used(&stderr, true);
    assert!(
        decompressed_records_without_pg(&from_file) == decompressed_records_without_pg(&from_fifo),
        "FIFO and file input sorted to different records"
    );
}

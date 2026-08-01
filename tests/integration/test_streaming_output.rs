//! Integration tests for streaming *output* support (`-o -`, `-o /dev/stdout`).
//!
//! [`test_input_source_matrix`](crate::test_input_source_matrix) holds every
//! command to the stdout axis once, through whatever code path its default
//! arguments select. That is the breadth; this file is the depth, for the two
//! commands whose writer is chosen by something other than `--output`:
//!
//! - `sort` picks a different writer per phase — an in-memory write when
//!   everything fits under `--max-memory`, a k-way merge when it does not — and
//!   `--write-index` selects a third that needs a seekable file and so cannot be
//!   satisfied by a pipe at all.
//! - `extract` writes through the unified pipeline when `--threads` makes it
//!   parallel and through a plain raw-BAM writer when it does not. Only one of
//!   those two used to reach stdout, which made the bug depend on a flag that has
//!   nothing to do with where the output goes.
//!
//! The failure mode throughout is silent: a writer that opens `-` as a filename
//! creates a regular file called `-`, writes the output into it, and exits zero
//! with an empty pipe. So every test here runs the command in its own working
//! directory and asserts nothing was left behind there, rather than trusting the
//! exit status.

use std::path::Path;
use std::process::{Command, Output};

use rstest::rstest;
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, create_umi_family, write_bam};
use crate::helpers::read_bam_output;

/// Runs `fgumi` with `args` from inside `cwd`, and returns the completed output.
///
/// The working directory matters: it is where a command that mistakes `-` for a
/// filename creates its file, and giving each run a fresh one is what lets
/// [`assert_streamed`] tell "streamed to the pipe" from "wrote a file called
/// `-`" without depending on the state of the crate root.
fn run_in<S: AsRef<std::ffi::OsStr>>(cwd: &Path, args: &[S]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args(args)
        .current_dir(cwd)
        .output()
        .expect("failed to run fgumi")
}

/// Asserts a run succeeded, streamed something, and left no file behind.
///
/// Returns the streamed bytes, so a caller can go on to check they are the same
/// output the command writes to a file. `label` names the case in any failure,
/// since every caller runs the same assertions over several spellings.
fn assert_streamed(run: &Output, cwd: &Path, label: &str) -> Vec<u8> {
    assert!(
        run.status.success(),
        "{label}: fgumi failed: {}",
        String::from_utf8_lossy(&run.stderr)
    );
    assert!(
        !cwd.join("-").exists(),
        "{label}: created a regular file named `-` instead of writing to stdout"
    );
    assert!(!run.stdout.is_empty(), "{label}: wrote nothing to stdout");
    run.stdout.clone()
}

/// Writes `stdout` to a file so it can be read back with the shared BAM helper,
/// which normalizes the `@PG` `CL` line the two runs legitimately differ on.
fn bam_from_stdout(dir: &Path, name: &str, stdout: &[u8]) -> std::path::PathBuf {
    let path = dir.join(name);
    std::fs::write(&path, stdout).expect("write captured stdout");
    path
}

/// A BAM of `families` UMI families, three reads each, at staggered positions.
///
/// Sized by caller: the sort tests need one input that fits in memory and one
/// that cannot, and the only difference between those two runs should be the
/// number of records.
fn write_test_bam(path: &Path, families: usize) {
    let header = create_minimal_header("chr1", 100_000);
    let records: Vec<_> = (0..families)
        .flat_map(|i| create_umi_family("ACGTACGT", 3, &format!("fam_{i:06}"), "ACGTACGTAC", 35))
        .collect();
    write_bam(path, &header, &records);
}

/// `sort` reaches stdout through every one of its write paths.
///
/// The phase that writes the output depends on whether the records fit under
/// `--max-memory`: under it, a single in-memory write; over it, a k-way merge of
/// spilled chunks. Those are two different writers, and only exercising the
/// small case would leave the merge path — the one every large sort takes —
/// untested. `--threads 2` additionally moves the write onto the pool.
///
/// Each case asserts the streamed bytes are the same BAM the file run produced,
/// so a path that streams *something* but drops or reorders records still fails.
#[rstest]
#[case::in_memory(64, "768M", "1")]
#[case::spill_and_merge(20_000, "1M", "1")]
#[case::spill_and_merge_threaded(20_000, "1M", "2")]
fn sort_streams_to_stdout(
    #[case] families: usize,
    #[case] max_memory: &str,
    #[case] threads: &str,
) {
    let dir = TempDir::new().expect("create temp dir");
    let input = dir.path().join("input.bam");
    write_test_bam(&input, families);

    let expected = dir.path().join("expected.bam");
    let file_args = [
        "sort",
        "-i",
        input.to_str().unwrap(),
        "-o",
        expected.to_str().unwrap(),
        "--order",
        "coordinate",
        "-m",
        max_memory,
        "--threads",
        threads,
    ];
    let file_run = run_in(dir.path(), &file_args);
    assert!(
        file_run.status.success(),
        "the file baseline failed, so the stdout comparison would prove nothing: {}",
        String::from_utf8_lossy(&file_run.stderr)
    );

    // A fresh working directory per streamed run, so a stray `-` is attributable.
    for spelling in ["-", "/dev/stdout"] {
        let cwd = dir.path().join(format!("cwd{spelling}").replace('/', "_"));
        std::fs::create_dir_all(&cwd).expect("create run dir");

        let mut args = file_args;
        args[4] = spelling;
        let run = run_in(&cwd, &args);
        let label = format!("sort -o {spelling} (m={max_memory}, threads={threads})");
        let stdout = assert_streamed(&run, &cwd, &label);

        let streamed = bam_from_stdout(dir.path(), "streamed.bam", &stdout);
        assert_eq!(
            read_bam_output(&streamed),
            read_bam_output(&expected),
            "{label}: streamed a different BAM than the same sort wrote to a file"
        );
    }
}

/// A merge that spilled really did spill.
///
/// [`sort_streams_to_stdout`]'s `spill_and_merge` cases are only meaningful if
/// the record count and memory limit actually push the sort past an in-memory
/// write. That is a property of the fixture, not of the code under test, so it
/// is asserted separately rather than assumed: if a future change to the buffer
/// sizing makes 20k records fit in 1M, this fails and says so instead of quietly
/// turning two of the three cases above into duplicates of the first.
#[test]
fn the_spilling_fixture_spills() {
    let dir = TempDir::new().expect("create temp dir");
    let input = dir.path().join("input.bam");
    write_test_bam(&input, 20_000);

    let run = run_in(
        dir.path(),
        &[
            "sort",
            "-i",
            input.to_str().unwrap(),
            "-o",
            dir.path().join("out.bam").to_str().unwrap(),
            "--order",
            "coordinate",
            "-m",
            "1M",
        ],
    );
    let log = String::from_utf8_lossy(&run.stderr);
    assert!(
        log.contains("Temporary chunks:"),
        "20k records under `-m 1M` no longer spill, so the merge path is not being covered:\n{log}"
    );
}

/// `--write-index` and stdout are mutually exclusive, and saying so beats
/// half-doing it.
///
/// A BAI has to be written to a sidecar path next to a seekable file, which a
/// pipe has neither of. The rejection is the same one `create_indexing_bam_writer`
/// already makes; the point of asserting it here is that it must happen *before*
/// anything is written, rather than the sort streaming a BAM and then failing to
/// index it — or, worse, falling back to opening `-` as a file.
#[rstest]
#[case::dash("-")]
#[case::dev_stdout("/dev/stdout")]
fn sort_rejects_write_index_with_stdout(#[case] spelling: &str) {
    let dir = TempDir::new().expect("create temp dir");
    let input = dir.path().join("input.bam");
    write_test_bam(&input, 8);

    let run = run_in(
        dir.path(),
        &[
            "sort",
            "-i",
            input.to_str().unwrap(),
            "-o",
            spelling,
            "--order",
            "coordinate",
            "--write-index",
        ],
    );

    assert!(!run.status.success(), "`--write-index -o {spelling}` must be rejected, not attempted");
    let stderr = String::from_utf8_lossy(&run.stderr);
    assert!(
        stderr.contains("--write-index"),
        "the error must name the flag that cannot be satisfied, got: {stderr}"
    );
    assert!(!dir.path().join("-").exists(), "the rejected run still created a file named `-`");
    assert!(!dir.path().join("-.bai").exists(), "the rejected run still created a `-.bai`");
}

/// `extract` reaches stdout whether or not `--threads` makes it parallel.
///
/// The two cases select different writers — a plain raw-BAM writer when the flag
/// is absent, the unified pipeline when it is present — so `-o -` streamed for a
/// bare `fgumi extract` and silently wrote a file named `-` once `--threads` was
/// given. A user has no reason to expect where their output lands to depend on
/// how many threads they asked for.
///
/// The first case omits `--threads` rather than passing `1`: threading mode is
/// `SingleThreaded` only when the flag is absent, and `--threads 1` is
/// `Threads(1)`, which takes the parallel path with a single worker. Passing `1`
/// would therefore have run the pipeline writer twice and never covered the raw
/// one.
#[rstest]
#[case::flag_omitted(None)]
#[case::parallel(Some("4"))]
fn extract_streams_to_stdout(#[case] threads: Option<&str>) {
    let dir = TempDir::new().expect("create temp dir");
    let input = dir.path().join("input.fq");
    let fastq = (0..8).fold(String::new(), |mut text, i| {
        use std::fmt::Write as _;
        writeln!(text, "@read_{i}\nACGTACGTACGTACGTAC\n+\nIIIIIIIIIIIIIIIIII")
            .expect("writing to a String cannot fail");
        text
    });
    std::fs::write(&input, fastq).expect("write FASTQ");

    let expected = dir.path().join("expected.bam");
    let args = |output: &str| {
        let mut args = vec![
            "extract".to_string(),
            "-i".to_string(),
            input.to_str().unwrap().to_string(),
            "-o".to_string(),
            output.to_string(),
            "-r".to_string(),
            "8M+T".to_string(),
            "--sample".to_string(),
            "S1".to_string(),
            "--library".to_string(),
            "L1".to_string(),
        ];
        if let Some(threads) = threads {
            args.extend(["--threads".to_string(), threads.to_string()]);
        }
        args
    };
    let file_args = args(expected.to_str().unwrap());
    let file_run = run_in(dir.path(), &file_args);
    assert!(
        file_run.status.success(),
        "the file baseline failed, so the stdout comparison would prove nothing: {}",
        String::from_utf8_lossy(&file_run.stderr)
    );

    let cwd = dir.path().join("cwd");
    std::fs::create_dir_all(&cwd).expect("create run dir");
    let dash_args = args("-");
    let run = run_in(&cwd, &dash_args);
    let label = format!("extract -o - (threads={})", threads.unwrap_or("<omitted>"));
    let stdout = assert_streamed(&run, &cwd, &label);

    let streamed = bam_from_stdout(dir.path(), "streamed.bam", &stdout);
    assert_eq!(
        read_bam_output(&streamed),
        read_bam_output(&expected),
        "{label}: streamed a different BAM than it wrote to a file"
    );
}

/// `fastq` accepts `-o -` for the stdout it already supports by omission.
///
/// Omitting `--output` entirely has always written plain FASTQ to stdout, so the
/// stream itself is not new — but `-o -` wrote a file named `-`, which makes the
/// command the odd one out among a CLI where `-` means stdout everywhere else.
/// Both spellings must produce exactly what the file run produced.
#[rstest]
#[case::dash(Some("-"))]
#[case::dev_stdout(Some("/dev/stdout"))]
#[case::omitted(None)]
fn fastq_streams_to_stdout(#[case] output: Option<&str>) {
    let dir = TempDir::new().expect("create temp dir");
    let input = dir.path().join("input.bam");
    write_test_bam(&input, 4);

    let expected = dir.path().join("expected.fq");
    let file_run = run_in(
        dir.path(),
        &["fastq", "-i", input.to_str().unwrap(), "-o", expected.to_str().unwrap()],
    );
    assert!(
        file_run.status.success(),
        "the file baseline failed, so the stdout comparison would prove nothing: {}",
        String::from_utf8_lossy(&file_run.stderr)
    );

    let cwd = dir.path().join("cwd");
    std::fs::create_dir_all(&cwd).expect("create run dir");
    let mut args = vec!["fastq", "-i", input.to_str().unwrap()];
    if let Some(spelling) = output {
        args.extend(["-o", spelling]);
    }
    let run = run_in(&cwd, &args);
    let label = format!("fastq -o {}", output.unwrap_or("<omitted>"));
    let stdout = assert_streamed(&run, &cwd, &label);

    assert_eq!(
        stdout,
        std::fs::read(&expected).expect("read the file baseline"),
        "{label}: streamed different FASTQ than the same command wrote to a file"
    );
}

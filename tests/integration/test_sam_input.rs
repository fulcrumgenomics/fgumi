//! Every command that accepts BAM must accept uncompressed SAM.
//!
//! fgumi opens input through several different readers — the single-threaded
//! raw-byte reader, the multi-threaded pipeline reader, and the sort engine's
//! pool-integrated reader — and a command picks between them based on flags
//! like `--threads` that have nothing to do with the input format. That makes
//! "does this command accept SAM?" easy to answer differently per command, and
//! even per flag combination within one command.
//!
//! These tests pin the invariant instead: for each command, the same data as
//! `.bam` and as `.sam` must both succeed and produce the same output — header
//! and records alike.

use bstr::BString;
use noodles::bam;
use noodles::sam;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::header::record::value::map::program::tag;
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};
use std::process::Command;

use rstest::rstest;
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, create_umi_family, write_bam};

/// Writes a small grouped-input BAM and the byte-for-byte equivalent SAM.
///
/// Both files hold the same records, so any difference in a command's output
/// between the two is a difference in how the input was read.
fn write_input_pair(dir: &Path) -> (PathBuf, PathBuf) {
    let bam_path = dir.join("input.bam");
    let sam_path = dir.join("input.sam");

    // `create_minimal_header` already advertises SO:unsorted GO:query
    // SS:template-coordinate, which is what the grouping commands require.
    let header = create_minimal_header("chr1", 10_000);

    let mut records = create_umi_family("ACGTACGT", 3, "fam_a", "ACGTACGTAC", 35);
    records.extend(create_umi_family("TTTTGGGG", 2, "fam_b", "ACGTACGTAC", 35));
    write_bam(&bam_path, &header, &records);

    transcode_bam_to_sam(&bam_path, &sam_path);

    (bam_path, sam_path)
}

/// Rewrites a BAM as uncompressed SAM text.
///
/// Done in-process with noodles rather than by shelling out to `samtools`, so
/// the test does not depend on an external tool being installed.
fn transcode_bam_to_sam(bam_path: &Path, sam_path: &Path) {
    let mut reader =
        bam::io::Reader::new(BufReader::new(File::open(bam_path).expect("open BAM to transcode")));
    let header = reader.read_header().expect("read BAM header");

    let mut writer = sam::io::Writer::new(File::create(sam_path).expect("create SAM"));
    writer.write_header(&header).expect("write SAM header");
    for result in reader.record_bufs(&header) {
        let record = result.expect("read BAM record");
        writer.write_alignment_record(&header, &record).expect("write SAM record");
    }
}

/// Reads a BAM's entire header and records, for comparing two runs' outputs.
///
/// Comparing the header as well as the records is what makes a dropped `@SQ`,
/// `@RG` or sort-order line a test failure rather than something only the
/// records would have to expose. The one part that legitimately differs between
/// runs is normalized away: fgumi stamps a `@PG` line whose `CL` field is the
/// command line, which names the input path.
fn read_output(path: &Path) -> (sam::Header, Vec<RecordBuf>) {
    let mut reader =
        bam::io::Reader::new(BufReader::new(File::open(path).expect("open output BAM")));
    let mut header = reader.read_header().expect("read output header");
    let records = reader
        .record_bufs(&header)
        .collect::<std::io::Result<Vec<_>>>()
        .expect("read output records");

    for (_, program) in header.programs_mut().as_mut() {
        if let Some(command_line) = program.other_fields_mut().get_mut(&tag::COMMAND_LINE) {
            *command_line = BString::from("<normalized>");
        }
    }

    (header, records)
}

/// Runs a command with `input` substituted for the `{input}` placeholder.
fn run(args: &[&str], input: &Path, output: &Path) -> std::process::Output {
    let input = input.to_str().expect("input path is UTF-8");
    let output = output.to_str().expect("output path is UTF-8");
    let args: Vec<&str> = args
        .iter()
        .map(|arg| match *arg {
            "{input}" => input,
            "{output}" => output,
            other => other,
        })
        .collect();

    Command::new(env!("CARGO_BIN_EXE_fgumi")).args(&args).output().expect("failed to run fgumi")
}

/// Runs a command with `input` fed in on stdin rather than named as a path.
fn run_with_stdin(args: &[&str], input: &Path, output: &Path) -> std::process::Output {
    let output = output.to_str().expect("output path is UTF-8");
    let args: Vec<&str> =
        args.iter().map(|arg| if *arg == "{output}" { output } else { *arg }).collect();

    let stdin = File::open(input).expect("open input to pipe");
    Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args(&args)
        .stdin(std::process::Stdio::from(stdin))
        .output()
        .expect("failed to run fgumi")
}

/// A command reading BAM must read the same data as SAM, and must produce the
/// same records from both.
///
/// `--threads` is varied per case because it selects between entirely separate
/// orchestrations (single-threaded fast path vs. pipeline); input format
/// support must not depend on which one runs.
#[rstest]
#[case::sort_single_threaded(&["sort", "-i", "{input}", "-o", "{output}", "--order", "queryname"])]
#[case::sort_multi_threaded(&["sort", "-i", "{input}", "-o", "{output}", "--order", "queryname", "--threads", "2"])]
#[case::group_single_threaded(&["group", "-i", "{input}", "-o", "{output}", "--strategy", "identity", "--edits", "0"])]
#[case::group_multi_threaded(&["group", "-i", "{input}", "-o", "{output}", "--strategy", "identity", "--edits", "0", "--threads", "2"])]
fn command_accepts_sam_input_identically_to_bam(#[case] args: &[&str]) {
    let dir = TempDir::new().expect("create temp dir");
    let (bam_input, sam_input) = write_input_pair(dir.path());

    let bam_output = dir.path().join("bam_output.bam");
    let bam_run = run(args, &bam_input, &bam_output);
    assert!(
        bam_run.status.success(),
        "BAM input failed: {}",
        String::from_utf8_lossy(&bam_run.stderr)
    );

    let sam_output = dir.path().join("sam_output.bam");
    let sam_run = run(args, &sam_input, &sam_output);
    assert!(
        sam_run.status.success(),
        "SAM input failed where BAM succeeded: {}",
        String::from_utf8_lossy(&sam_run.stderr)
    );

    assert_eq!(
        read_output(&sam_output),
        read_output(&bam_output),
        "SAM and BAM input produced different output",
    );
}

/// Stdin is the one input that cannot be rewound or reopened, which is why the
/// header is parsed through a tee and replayed. Piping SAM in is therefore the
/// case that a regression to seek-based header parsing would break first.
#[rstest]
#[case::sam_on_stdin("input.sam")]
#[case::bam_on_stdin("input.bam")]
fn sort_accepts_piped_input_identically_to_a_file(#[case] piped: &str) {
    let dir = TempDir::new().expect("create temp dir");
    let (bam_input, _) = write_input_pair(dir.path());

    let piped_output = dir.path().join("piped_output.bam");
    let piped_run = run_with_stdin(
        &["sort", "-i", "-", "-o", "{output}", "--order", "queryname"],
        &dir.path().join(piped),
        &piped_output,
    );
    assert!(
        piped_run.status.success(),
        "piped input failed: {}",
        String::from_utf8_lossy(&piped_run.stderr)
    );

    let file_output = dir.path().join("file_output.bam");
    let file_run = run(
        &["sort", "-i", "{input}", "-o", "{output}", "--order", "queryname"],
        &bam_input,
        &file_output,
    );
    assert!(
        file_run.status.success(),
        "BAM file input failed: {}",
        String::from_utf8_lossy(&file_run.stderr)
    );

    assert_eq!(
        read_output(&piped_output),
        read_output(&file_output),
        "piped input produced different output than the BAM file run",
    );
}

/// Format detection reads the input's magic bytes, so a file's extension
/// cannot make a command misread it. A `.bam` holding SAM text and a `.sam`
/// holding BGZF both have to work.
#[rstest]
#[case::sam_text_named_bam("misnamed.bam", false)]
#[case::bgzf_named_sam("misnamed.sam", true)]
fn input_format_follows_content_not_extension(#[case] name: &str, #[case] write_bgzf: bool) {
    let dir = TempDir::new().expect("create temp dir");
    let (bam_input, sam_input) = write_input_pair(dir.path());
    let args = &["sort", "-i", "{input}", "-o", "{output}", "--order", "queryname"];

    let misnamed = dir.path().join(name);
    let source = if write_bgzf { &bam_input } else { &sam_input };
    std::fs::copy(source, &misnamed).expect("copy to misnamed path");

    let output = dir.path().join("out.bam");
    let result = run(args, &misnamed, &output);
    assert!(
        result.status.success(),
        "misnamed input rejected: {}",
        String::from_utf8_lossy(&result.stderr)
    );

    // Succeeding is not enough: the run must produce what the correctly named
    // input produces, so a misread that silently drops records fails here.
    let canonical_output = dir.path().join("canonical.bam");
    let canonical = run(args, &bam_input, &canonical_output);
    assert!(
        canonical.status.success(),
        "canonical BAM input rejected: {}",
        String::from_utf8_lossy(&canonical.stderr)
    );

    assert_eq!(
        read_output(&output),
        read_output(&canonical_output),
        "misnamed input produced different output than the canonical BAM run",
    );
}

/// `fastq` reads through the single-threaded raw-byte reader and writes FASTQ,
/// so it exercises a reader surface the BAM-output cases do not and its output
/// can be compared byte for byte.
#[test]
fn fastq_command_accepts_sam_input_identically_to_bam() {
    let dir = TempDir::new().expect("create temp dir");
    let (bam_input, sam_input) = write_input_pair(dir.path());

    let bam_output = dir.path().join("bam_output.fq");
    let bam_run = run(&["fastq", "-i", "{input}", "-o", "{output}"], &bam_input, &bam_output);
    assert!(
        bam_run.status.success(),
        "BAM input failed: {}",
        String::from_utf8_lossy(&bam_run.stderr)
    );

    let sam_output = dir.path().join("sam_output.fq");
    let sam_run = run(&["fastq", "-i", "{input}", "-o", "{output}"], &sam_input, &sam_output);
    assert!(
        sam_run.status.success(),
        "SAM input failed where BAM succeeded: {}",
        String::from_utf8_lossy(&sam_run.stderr)
    );

    assert_eq!(
        std::fs::read(&sam_output).expect("read SAM-derived FASTQ"),
        std::fs::read(&bam_output).expect("read BAM-derived FASTQ"),
        "SAM and BAM input produced different FASTQ",
    );
}

//! Integration tests for streaming input support (stdin/pipes).
//!
//! These tests spawn cat processes whose stdout is piped to fgumi commands.
//! The child processes are properly cleaned up when their stdout is consumed.
#![allow(clippy::zombie_processes)]

use noodles::bam;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use std::fs::{self, File};
use std::io::{BufReader, Read};
use std::path::PathBuf;
use std::process::{Command, Stdio};

use rstest::rstest;
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, create_umi_family, to_record_buf};

/// Test that the group command works correctly with piped input.
#[test]
fn test_group_command_with_piped_input() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_from_file = temp_dir.path().join("output_file.bam");
    let output_from_pipe = temp_dir.path().join("output_pipe.bam");

    // Create test BAM file
    create_test_input_bam(&input_bam);

    // Run group command with file input (baseline)
    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "group",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_from_file.to_str().unwrap(),
            "--strategy",
            "identity",
            "--edits",
            "0",
            "--compression-level",
            "1",
            "--threads",
            "2",
        ])
        .status()
        .expect("Failed to run group command with file input");
    assert!(status.success(), "Group command with file input failed");
    assert!(output_from_file.exists(), "Output BAM from file not created");

    // Run group command with piped input (cat file | fgumi group)
    let cat_child = Command::new("cat")
        .arg(input_bam.to_str().unwrap())
        .stdout(Stdio::piped())
        .spawn()
        .expect("Failed to spawn cat");

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "group",
            "--input",
            "-", // stdin
            "--output",
            output_from_pipe.to_str().unwrap(),
            "--strategy",
            "identity",
            "--edits",
            "0",
            "--compression-level",
            "1",
            "--threads",
            "2",
        ])
        .stdin(cat_child.stdout.unwrap())
        .status()
        .expect("Failed to run group command with piped input");
    assert!(status.success(), "Group command with piped input failed");
    assert!(output_from_pipe.exists(), "Output BAM from pipe not created");

    // Compare outputs - records should be identical (headers may differ due to @PG command line)
    compare_bam_records(&output_from_file, &output_from_pipe);
}

/// Test that /dev/stdin path also works.
#[test]
fn test_group_command_with_dev_stdin_path() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    // Create test BAM file
    create_test_input_bam(&input_bam);

    // Run group command using /dev/stdin
    let cat_child = Command::new("cat")
        .arg(input_bam.to_str().unwrap())
        .stdout(Stdio::piped())
        .spawn()
        .expect("Failed to spawn cat");

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "group",
            "--input",
            "/dev/stdin",
            "--output",
            output_bam.to_str().unwrap(),
            "--strategy",
            "identity",
            "--edits",
            "0",
            "--compression-level",
            "1",
            "--threads",
            "2",
        ])
        .stdin(cat_child.stdout.unwrap())
        .status()
        .expect("Failed to run group command with /dev/stdin");

    assert!(status.success(), "Group command with /dev/stdin failed");
    assert!(output_bam.exists(), "Output BAM not created");
}

/// Test simplex command with piped input.
#[test]
fn test_simplex_command_with_piped_input() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_from_file = temp_dir.path().join("output_file.bam");
    let output_from_pipe = temp_dir.path().join("output_pipe.bam");

    // Create test BAM file with grouped reads (MI tag)
    create_grouped_test_bam(&input_bam);

    // Run simplex with file input
    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "simplex",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_from_file.to_str().unwrap(),
            "--min-reads",
            "1",
            "--min-consensus-base-quality",
            "0",
            "--compression-level",
            "1",
            "--threads",
            "2",
        ])
        .status()
        .expect("Failed to run simplex command with file input");
    assert!(status.success(), "Simplex command with file input failed");
    assert!(output_from_file.exists(), "Output BAM from file not created");

    // Run simplex with piped input
    let cat_child = Command::new("cat")
        .arg(input_bam.to_str().unwrap())
        .stdout(Stdio::piped())
        .spawn()
        .expect("Failed to spawn cat");

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "simplex",
            "--input",
            "-",
            "--output",
            output_from_pipe.to_str().unwrap(),
            "--min-reads",
            "1",
            "--min-consensus-base-quality",
            "0",
            "--compression-level",
            "1",
            "--threads",
            "2",
        ])
        .stdin(cat_child.stdout.unwrap())
        .status()
        .expect("Failed to run simplex command with piped input");
    assert!(status.success(), "Simplex command with piped input failed");
    assert!(output_from_pipe.exists(), "Output BAM from pipe not created");

    // Guard against a vacuous pass: simplex must actually emit consensus records,
    // otherwise the file-vs-pipe comparison below is trivially true (0 == 0) even
    // if stdin were dropped entirely.
    assert!(
        count_bam_records(&output_from_file) > 0,
        "simplex produced no consensus records — parity check would be vacuous"
    );

    // Compare outputs - records should be identical (headers may differ due to @PG command line)
    compare_bam_records(&output_from_file, &output_from_pipe);
}

/// Convert a BAM file to a parallel SAM text file containing the same
/// records. Lets `fgumi group --input foo.sam` reuse the same record
/// set as the BAM-input baseline for parity tests.
fn convert_bam_to_sam(bam_path: &PathBuf, sam_path: &PathBuf) {
    use noodles::sam::alignment::io::Write as _;
    let mut bam_reader = bam::io::reader::Builder.build_from_path(bam_path).expect("open bam");
    let header = bam_reader.read_header().expect("read bam header");
    let mut sam_writer =
        noodles::sam::io::Writer::new(fs::File::create(sam_path).expect("create sam"));
    sam_writer.write_header(&header).expect("write sam header");
    for record in bam_reader.records() {
        let record = record.expect("bam record");
        sam_writer.write_alignment_record(&header, &record).expect("write sam record");
    }
}

/// SAM-input parity for the new typed-step pipeline. Generates a SAM file
/// from the same record set as the BAM baseline, runs `fgumi group
/// --use-new-pipeline --input foo.sam`, and compares the output BAMs.
/// Exercises the `ReadSamChunks` + `ParseSamChunk` parallel-parse path.
#[test]
fn test_group_command_with_sam_input_new_pipeline_matches_bam_baseline() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let bam_path = temp_dir.path().join("input.bam");
    let sam_path = temp_dir.path().join("input.sam");
    let out_bam_baseline = temp_dir.path().join("output_bam.bam");
    let out_sam = temp_dir.path().join("output_sam.bam");

    create_test_input_bam(&bam_path);
    convert_bam_to_sam(&bam_path, &sam_path);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "group",
            "--input",
            bam_path.to_str().unwrap(),
            "--output",
            out_bam_baseline.to_str().unwrap(),
            "--strategy",
            "identity",
            "--edits",
            "0",
            "--compression-level",
            "1",
            "--threads",
            "2",
        ])
        .status()
        .expect("Failed to run group --use-new-pipeline with BAM input");
    assert!(status.success(), "BAM-input baseline failed");

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "group",
            "--input",
            sam_path.to_str().unwrap(),
            "--output",
            out_sam.to_str().unwrap(),
            "--strategy",
            "identity",
            "--edits",
            "0",
            "--compression-level",
            "1",
            "--threads",
            "2",
        ])
        .status()
        .expect("Failed to run group --use-new-pipeline with SAM input");
    assert!(status.success(), "SAM-input run failed");

    compare_bam_records(&out_bam_baseline, &out_sam);
}

/// SAM-input parity for simplex. Generates a SAM from the same grouped
/// records as the BAM baseline (with MI tags), runs simplex on both,
/// compares the consensus outputs.
#[test]
fn test_simplex_command_with_sam_input_new_pipeline_matches_bam_baseline() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let bam_path = temp_dir.path().join("input.bam");
    let sam_path = temp_dir.path().join("input.sam");
    let out_bam_baseline = temp_dir.path().join("output_bam.bam");
    let out_sam = temp_dir.path().join("output_sam.bam");

    create_grouped_test_bam(&bam_path);
    convert_bam_to_sam(&bam_path, &sam_path);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "simplex",
            "--input",
            bam_path.to_str().unwrap(),
            "--output",
            out_bam_baseline.to_str().unwrap(),
            "--min-reads",
            "1",
            "--min-consensus-base-quality",
            "0",
            "--compression-level",
            "1",
            "--threads",
            "2",
        ])
        .status()
        .expect("BAM baseline run");
    assert!(status.success(), "simplex BAM baseline failed");

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "simplex",
            "--input",
            sam_path.to_str().unwrap(),
            "--output",
            out_sam.to_str().unwrap(),
            "--min-reads",
            "1",
            "--min-consensus-base-quality",
            "0",
            "--compression-level",
            "1",
            "--threads",
            "2",
        ])
        .status()
        .expect("SAM run");
    assert!(status.success(), "simplex SAM run failed");

    compare_bam_records(&out_bam_baseline, &out_sam);
}

/// SAM-input parity for `correct` — exercises the `(Sam, true)` arm of
/// correct.rs's 4-way match (the `into_multi()` shape that emits both
/// kept records and rejects on parallel branches).
#[test]
fn test_correct_command_with_sam_input_new_pipeline_with_rejects_matches_bam_baseline() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let bam_path = temp_dir.path().join("input.bam");
    let sam_path = temp_dir.path().join("input.sam");
    let umi_file = temp_dir.path().join("umis.txt");
    let out_bam_baseline = temp_dir.path().join("output_bam.bam");
    let out_bam_rejects_baseline = temp_dir.path().join("rejects_bam.bam");
    let out_sam = temp_dir.path().join("output_sam.bam");
    let out_sam_rejects = temp_dir.path().join("rejects_sam.bam");

    create_test_input_bam(&bam_path);
    convert_bam_to_sam(&bam_path, &sam_path);
    // The UMIs used by `create_umi_family` are AAAAAAAA and CCCCCCCC.
    std::fs::write(&umi_file, b"AAAAAAAA\nCCCCCCCC\n").expect("write umi file");

    for (input, output, rejects, label) in [
        (&bam_path, &out_bam_baseline, &out_bam_rejects_baseline, "BAM"),
        (&sam_path, &out_sam, &out_sam_rejects, "SAM"),
    ] {
        let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
            .args([
                "correct",
                "--input",
                input.to_str().unwrap(),
                "--output",
                output.to_str().unwrap(),
                "--rejects",
                rejects.to_str().unwrap(),
                "--umi-files",
                umi_file.to_str().unwrap(),
                "--max-mismatches",
                "1",
                "--min-distance",
                "1",
                "--compression-level",
                "1",
                "--threads",
                "2",
            ])
            .status()
            .unwrap_or_else(|_| panic!("Failed to run correct {label}"));
        assert!(status.success(), "correct {label} failed");
    }
    compare_bam_records(&out_bam_baseline, &out_sam);
    compare_bam_records(&out_bam_rejects_baseline, &out_sam_rejects);
}

/// Build an unmapped-consensus BAM with depth tags so filter accepts it
/// without `--ref` (mapped reads would require ref for NM/UQ/MD tag
/// regeneration). Mirrors the fixture in `test_filter_command.rs`.
fn create_unmapped_consensus_bam(path: &PathBuf) {
    use fgumi_lib::sam::SamTag;
    use fgumi_raw_bam::{SamBuilder as RawSamBuilder, flags};

    let header = create_minimal_header("chr1", 10000);
    let records: Vec<_> = ["good1", "good2", "low_depth"]
        .iter()
        .enumerate()
        .map(|(i, name)| {
            let depth = if *name == "low_depth" { 1u16 } else { 10u16 };
            let mut b = RawSamBuilder::new();
            b.read_name(name.as_bytes())
                .flags(flags::UNMAPPED)
                .sequence(b"ACGTACGT")
                .qualities(&[35; 8]);
            b.add_array_u16(SamTag::CD_BASES, &[depth; 8]).add_array_u16(SamTag::CE_BASES, &[0; 8]);
            let _ = i; // silence unused
            b.build()
        })
        .collect();

    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");
    for raw in records {
        let rec =
            fgumi_raw_bam::raw_record_to_record_buf(&raw, &header).expect("raw -> record_buf");
        writer.write_alignment_record(&header, &rec).expect("Failed to write record");
    }
    writer.try_finish().expect("Failed to finish BAM");
}

/// SAM-input parity for `filter --rejects` — exercises the `(Sam, ...)`
/// arm of `execute_single_read_with_rejects_new` (the `into_multi()`
/// fan-out shape that emits kept + rejected on parallel branches).
#[test]
fn test_filter_command_with_rejects_sam_input_matches_bam_baseline() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let bam_path = temp_dir.path().join("input.bam");
    let sam_path = temp_dir.path().join("input.sam");
    let out_bam_baseline = temp_dir.path().join("output_bam.bam");
    let out_bam_rejects_baseline = temp_dir.path().join("rejects_bam.bam");
    let out_sam = temp_dir.path().join("output_sam.bam");
    let out_sam_rejects = temp_dir.path().join("rejects_sam.bam");

    create_unmapped_consensus_bam(&bam_path);
    convert_bam_to_sam(&bam_path, &sam_path);

    for (input, output, rejects, label) in [
        (&bam_path, &out_bam_baseline, &out_bam_rejects_baseline, "BAM"),
        (&sam_path, &out_sam, &out_sam_rejects, "SAM"),
    ] {
        let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
            .args([
                "filter",
                "--input",
                input.to_str().unwrap(),
                "--output",
                output.to_str().unwrap(),
                "--rejects",
                rejects.to_str().unwrap(),
                "--min-mean-base-quality",
                "0",
                "--min-reads",
                "3",
                "--compression-level",
                "1",
                "--threads",
                "2",
            ])
            .status()
            .unwrap_or_else(|_| panic!("Failed to run filter {label}"));
        assert!(status.success(), "filter {label} failed");
    }
    compare_bam_records(&out_bam_baseline, &out_sam);
    compare_bam_records(&out_bam_rejects_baseline, &out_sam_rejects);
}

/// SAM-on-stdin parity (`cat foo.sam | fgumi group --input -`). The
/// `InputSource::open` magic-byte peek picks the SAM-chunks source for
/// stdin (no `\x1f` BGZF magic, falls through to text parsing).
#[test]
fn test_group_command_with_sam_stdin_new_pipeline_matches_bam_baseline() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let bam_path = temp_dir.path().join("input.bam");
    let sam_path = temp_dir.path().join("input.sam");
    let out_bam_baseline = temp_dir.path().join("output_bam.bam");
    let out_sam_stdin = temp_dir.path().join("output_sam_stdin.bam");

    create_test_input_bam(&bam_path);
    convert_bam_to_sam(&bam_path, &sam_path);

    // BAM baseline (file input).
    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "group",
            "--input",
            bam_path.to_str().unwrap(),
            "--output",
            out_bam_baseline.to_str().unwrap(),
            "--strategy",
            "identity",
            "--edits",
            "0",
            "--compression-level",
            "1",
            "--threads",
            "2",
        ])
        .status()
        .expect("BAM baseline run");
    assert!(status.success(), "BAM-input baseline failed");

    // SAM piped through stdin.
    let cat_child = Command::new("cat")
        .arg(sam_path.to_str().unwrap())
        .stdout(Stdio::piped())
        .spawn()
        .expect("spawn cat");
    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "group",
            "--input",
            "-",
            "--output",
            out_sam_stdin.to_str().unwrap(),
            "--strategy",
            "identity",
            "--edits",
            "0",
            "--compression-level",
            "1",
            "--threads",
            "2",
        ])
        .stdin(cat_child.stdout.unwrap())
        .status()
        .expect("SAM stdin run");
    assert!(status.success(), "SAM-stdin run failed");

    compare_bam_records(&out_bam_baseline, &out_sam_stdin);
}

/// Same shape as `test_group_command_with_piped_input` but routed through
/// the typed-step pipeline (`--use-new-pipeline`). Exercises the
/// `read_bam_auto` dispatcher introduced for issue #330 Phase 1 T1.3.4:
/// when the input path is `-`, the source step is `read_bam_stdin`
/// (which buffers the BAM header via `TeeReader` and replays it ahead
/// of the remaining stdin stream).
#[test]
fn test_group_command_with_piped_input_new_pipeline() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_from_file = temp_dir.path().join("output_file.bam");
    let output_from_pipe = temp_dir.path().join("output_pipe.bam");

    create_test_input_bam(&input_bam);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "group",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_from_file.to_str().unwrap(),
            "--strategy",
            "identity",
            "--edits",
            "0",
            "--compression-level",
            "1",
            "--threads",
            "2",
        ])
        .status()
        .expect("Failed to run group --use-new-pipeline with file input");
    assert!(status.success(), "group --use-new-pipeline (file) failed");
    assert!(output_from_file.exists(), "file-input output not created");

    let cat_child = Command::new("cat")
        .arg(input_bam.to_str().unwrap())
        .stdout(Stdio::piped())
        .spawn()
        .expect("Failed to spawn cat");

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "group",
            "--input",
            "-",
            "--output",
            output_from_pipe.to_str().unwrap(),
            "--strategy",
            "identity",
            "--edits",
            "0",
            "--compression-level",
            "1",
            "--threads",
            "2",
        ])
        .stdin(cat_child.stdout.unwrap())
        .status()
        .expect("Failed to run group --use-new-pipeline with piped input");
    assert!(status.success(), "group --use-new-pipeline (stdin) failed");
    assert!(output_from_pipe.exists(), "stdin output not created");

    compare_bam_records(&output_from_file, &output_from_pipe);
}

/// Helper to read file contents for comparison.
#[allow(dead_code)]
fn read_file_contents(path: &PathBuf) -> Vec<u8> {
    let file = File::open(path).expect("Failed to open file");
    let mut reader = BufReader::new(file);
    let mut contents = Vec::new();
    reader.read_to_end(&mut contents).expect("Failed to read file");
    contents
}

/// Helper to compare BAM records (ignoring header differences like @PG command line).
///
/// Decodes to eager `RecordBuf`s and compares them for full equality — every
/// field including bases, qualities, flags, CIGAR, positions, and aux tags — so
/// a record-level corruption from mishandled input (dropped or garbled bytes)
/// fails here, not just a count/name mismatch.
fn compare_bam_records(path1: &PathBuf, path2: &PathBuf) {
    use noodles::sam::alignment::RecordBuf;

    let mut reader1 =
        bam::io::reader::Builder.build_from_path(path1).expect("Failed to open BAM 1");
    let header1 = reader1.read_header().expect("Failed to read header 1");

    let mut reader2 =
        bam::io::reader::Builder.build_from_path(path2).expect("Failed to open BAM 2");
    let header2 = reader2.read_header().expect("Failed to read header 2");

    let records1: Vec<RecordBuf> =
        reader1.record_bufs(&header1).map(|r| r.expect("Failed to read record 1")).collect();
    let records2: Vec<RecordBuf> =
        reader2.record_bufs(&header2).map(|r| r.expect("Failed to read record 2")).collect();

    assert_eq!(records1.len(), records2.len(), "BAM files have different record counts");

    for (i, (r1, r2)) in records1.iter().zip(records2.iter()).enumerate() {
        assert_eq!(r1, r2, "Record {i} differs (full record identity)");
    }
}

/// Create a test BAM with UMI families for grouping.
fn create_test_input_bam(path: &PathBuf) {
    let header = create_minimal_header("chr1", 10000);

    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");

    // Create multiple UMI families
    let family1 = create_umi_family("AAAAAAAA", 5, "family1", "ACGTACGT", 30);
    let family2 = create_umi_family("CCCCCCCC", 3, "family2", "TGCATGCA", 30);

    for record in family1.iter().chain(family2.iter()) {
        writer
            .write_alignment_record(&header, &to_record_buf(record))
            .expect("Failed to write record");
    }

    writer.try_finish().expect("Failed to finish BAM");
}

/// Create a test BAM with already-grouped reads (MI tag set).
fn create_grouped_test_bam(path: &PathBuf) {
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::data::field::Value;

    let header = create_minimal_header("chr1", 10000);

    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");

    // Create grouped reads with MI tag.
    //
    // MI MUST be a STRING tag: `fgumi group` writes the MoleculeId as a string,
    // and the consensus callers group reads by reading MI as a string (see
    // `MiGrouper::get_mi_tag` via `find_string_tag_in_record`). An integer MI is
    // invisible to the grouper — every read is dropped, so simplex/duplex emit
    // zero consensus records (which would silently make the parity tests vacuous).
    let mi_tag = Tag::from(fgumi_lib::sam::SamTag::MI);

    let records = create_umi_family("AAAAAAAA", 5, "mol1", "ACGTACGT", 30);
    for raw in &records {
        let mut rec = to_record_buf(raw);
        rec.data_mut().insert(mi_tag, Value::String("1".into()));
        writer.write_alignment_record(&header, &rec).expect("Failed to write record");
    }

    // Second family with different MI
    let records2 = create_umi_family("CCCCCCCC", 3, "mol2", "TGCATGCA", 30);
    for raw in &records2 {
        let mut rec = to_record_buf(raw);
        rec.data_mut().insert(mi_tag, Value::String("2".into()));
        writer.write_alignment_record(&header, &rec).expect("Failed to write record");
    }

    writer.try_finish().expect("Failed to finish BAM");
}

/// Run `fgumi` with `args`; if `stdin_bam` is `Some`, pipe that BAM to its
/// stdin via `cat`. Asserts the process exits successfully.
fn run_fgumi_checked(args: &[String], stdin_bam: Option<&std::path::Path>) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_fgumi"));
    command.args(args);
    let status = if let Some(bam) = stdin_bam {
        let cat = Command::new("cat")
            .arg(bam.to_str().unwrap())
            .stdout(Stdio::piped())
            .spawn()
            .expect("Failed to spawn cat");
        command.stdin(cat.stdout.unwrap()).status()
    } else {
        command.status()
    };
    let status = status.unwrap_or_else(|e| panic!("Failed to run fgumi {args:?}: {e}"));
    assert!(status.success(), "fgumi {args:?} failed (stdin_piped={})", stdin_bam.is_some());
}

/// Assert a command emits byte-identical records reading `input_bam` from a
/// file vs piped stdin (`-`) vs `/dev/stdin`. `build_args(input, output)`
/// returns the full argv for one invocation; `label` namespaces temp outputs.
fn assert_stdin_input_parity(
    input_bam: &std::path::Path,
    temp: &std::path::Path,
    label: &str,
    build_args: impl Fn(&str, &str) -> Vec<String>,
) {
    let out_file = temp.join(format!("{label}_file.bam"));
    let out_dash = temp.join(format!("{label}_dash.bam"));
    let out_dev = temp.join(format!("{label}_dev.bam"));

    run_fgumi_checked(&build_args(input_bam.to_str().unwrap(), out_file.to_str().unwrap()), None);
    run_fgumi_checked(&build_args("-", out_dash.to_str().unwrap()), Some(input_bam));
    run_fgumi_checked(&build_args("/dev/stdin", out_dev.to_str().unwrap()), Some(input_bam));

    // Guard against a vacuous pass: if the command emitted no records, a
    // file-vs-stdin parity check would be trivially true even if stdin were
    // dropped entirely. Require the file baseline to produce real output.
    assert!(
        count_bam_records(&out_file) > 0,
        "{label}: file baseline produced no records — parity check would be vacuous"
    );

    compare_bam_records(&out_file, &out_dash);
    compare_bam_records(&out_file, &out_dev);
}

/// Count records in a BAM (eager decode). Used to reject vacuous parity tests.
fn count_bam_records(path: &std::path::Path) -> usize {
    let mut reader = bam::io::reader::Builder.build_from_path(path).expect("open BAM for count");
    let header = reader.read_header().expect("read header for count");
    let mut count = 0;
    for record in reader.record_bufs(&header) {
        record.expect("read record for count");
        count += 1;
    }
    count
}

fn s(v: &str) -> String {
    v.to_string()
}

#[rstest]
#[case("1")]
#[case("2")]
fn test_sort_reads_stdin_once(#[case] threads: &str) {
    let dir = TempDir::new().unwrap();
    let input = dir.path().join("input.bam");
    create_test_input_bam(&input);
    assert_stdin_input_parity(&input, dir.path(), "sort", |i, o| {
        vec![
            s("sort"),
            s("--input"),
            s(i),
            s("--output"),
            s(o),
            s("--order"),
            s("coordinate"),
            s("--threads"),
            s(threads),
            s("--compression-level"),
            s("1"),
        ]
    });
}

/// `sort --verify` reads its input twice (header probe + a fresh record
/// re-scan via `File::open`), which a non-seekable stdin stream can't satisfy,
/// so it must be rejected up front with a clear message rather than failing
/// deep in a re-open. (Plain `sort` streams stdin once and IS supported.)
///
/// The rejection keys off `is_stdin_path`, so both `-` and `/dev/stdin` must be
/// rejected — `/dev/stdin` is a real path (`Path::exists` is `true`), so a gate
/// keyed on the literal `-` would let it slip into the double-read path.
#[rstest]
#[case("-")]
#[case("/dev/stdin")]
fn test_sort_verify_rejects_stdin_with_clear_message(#[case] input_arg: &str) {
    let dir = TempDir::new().unwrap();
    let input = dir.path().join("input.bam");
    create_test_input_bam(&input);
    let cat = Command::new("cat")
        .arg(input.to_str().unwrap())
        .stdout(Stdio::piped())
        .spawn()
        .expect("spawn cat");
    let output = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args(["sort", "--verify", "--input", input_arg, "--order", "coordinate"])
        .stdin(cat.stdout.unwrap())
        .output()
        .unwrap_or_else(|e| panic!("run sort --verify -i {input_arg}: {e}"));
    assert!(!output.status.success(), "sort --verify from stdin ({input_arg}) must be rejected");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("stdin"),
        "sort --verify stdin ({input_arg}) rejection should mention stdin; stderr: {stderr}"
    );
}

/// Single-threaded `simplex -i -` exercises the shared `open_bgzf_reader` stdin
/// path (the multi-threaded path does not — which is how the single-threaded
/// `File::open("-")` double-open regression shipped). With `--async-reader`, the
/// stdin branch must also spawn the prefetch thread and still read correctly.
/// `assert_stdin_input_parity` covers both `-` and `/dev/stdin` and guards
/// against a vacuous pass (it requires the file baseline to emit records).
#[rstest]
#[case::sync(false)]
#[case::async_reader(true)]
fn test_simplex_single_threaded_reads_stdin_once(#[case] async_reader: bool) {
    let dir = TempDir::new().expect("Failed to create temp dir");
    let input = dir.path().join("input.bam");
    create_grouped_test_bam(&input);
    // No `--threads` → single-threaded fast path.
    assert_stdin_input_parity(&input, dir.path(), "simplex", move |i, o| {
        let mut args = vec![
            s("simplex"),
            s("--input"),
            s(i),
            s("--output"),
            s(o),
            s("--min-reads"),
            s("1"),
            s("--min-consensus-base-quality"),
            s("0"),
            s("--compression-level"),
            s("1"),
        ];
        if async_reader {
            args.push(s("--async-reader"));
        }
        args
    });
}

fn push_threads(args: &mut Vec<String>, threads: Option<&str>) {
    if let Some(t) = threads {
        args.push(s("--threads"));
        args.push(s(t));
    }
}

#[rstest]
#[case("1")]
#[case("2")]
fn test_group_reads_stdin_once(#[case] threads: &str) {
    let dir = TempDir::new().unwrap();
    let input = dir.path().join("input.bam");
    create_test_input_bam(&input);
    assert_stdin_input_parity(&input, dir.path(), "group", |i, o| {
        vec![
            s("group"),
            s("--input"),
            s(i),
            s("--output"),
            s(o),
            s("--strategy"),
            s("identity"),
            s("--edits"),
            s("0"),
            s("--threads"),
            s(threads),
            s("--compression-level"),
            s("1"),
        ]
    });
}

#[rstest]
#[case("1")]
#[case("2")]
fn test_dedup_reads_stdin_once(#[case] threads: &str) {
    let dir = TempDir::new().unwrap();
    let input = dir.path().join("input.bam");
    create_test_input_bam(&input);
    assert_stdin_input_parity(&input, dir.path(), "dedup", |i, o| {
        vec![
            s("dedup"),
            s("--input"),
            s(i),
            s("--output"),
            s(o),
            s("--strategy"),
            s("identity"),
            s("--edits"),
            s("0"),
            s("--threads"),
            s(threads),
            s("--compression-level"),
            s("1"),
        ]
    });
}

/// correct has a single-threaded fast path (no `--threads`) plus the
/// multi-threaded path; cover both (single-threaded exercises the reader reuse).
#[rstest]
#[case::single_threaded(None)]
#[case::multi_threaded(Some("2"))]
fn test_correct_reads_stdin_once(#[case] threads: Option<&str>) {
    let dir = TempDir::new().unwrap();
    let input = dir.path().join("input.bam");
    create_test_input_bam(&input);
    // create_test_input_bam tags reads with RX = AAAAAAAA / CCCCCCCC.
    assert_stdin_input_parity(&input, dir.path(), "correct", move |i, o| {
        let mut args = vec![
            s("correct"),
            s("--input"),
            s(i),
            s("--output"),
            s(o),
            s("--umis"),
            s("AAAAAAAA"),
            s("--umis"),
            s("CCCCCCCC"),
            s("--min-distance"),
            s("1"),
            s("--compression-level"),
            s("1"),
        ];
        push_threads(&mut args, threads);
        args
    });
}

/// codec has a single-threaded fast path (no `--threads`) plus the
/// multi-threaded path; cover both.
#[rstest]
#[case::single_threaded(None)]
#[case::multi_threaded(Some("2"))]
fn test_codec_reads_stdin_once(#[case] threads: Option<&str>) {
    let dir = TempDir::new().unwrap();
    let input = dir.path().join("codec.bam");
    let pairs: Vec<_> = (0..3)
        .map(|i| {
            crate::test_codec_command::create_codec_read_pair(
                &format!("read{i}"),
                b"ACGTACGT",
                b"ACGTACGT",
                &[30; 8],
                &[30; 8],
                100,
                "UMI001",
                None,
            )
        })
        .collect();
    crate::test_codec_command::create_codec_test_bam(&input, pairs);
    assert_stdin_input_parity(&input, dir.path(), "codec", move |i, o| {
        let mut args = vec![
            s("codec"),
            s("--input"),
            s(i),
            s("--output"),
            s(o),
            s("--min-reads"),
            s("1"),
            s("--min-duplex-length"),
            s("1"),
            s("--compression-level"),
            s("1"),
        ];
        push_threads(&mut args, threads);
        args
    });
}

/// duplex has a single-threaded fast path (no `--threads`) plus the
/// multi-threaded path; cover both.
#[rstest]
#[case::single_threaded(None)]
#[case::multi_threaded(Some("2"))]
fn test_duplex_reads_stdin_once(#[case] threads: Option<&str>) {
    let dir = TempDir::new().unwrap();
    let input = dir.path().join("duplex.bam");
    let molecule =
        crate::test_duplex_command::create_duplex_molecule("MI001", "ACGTACGT", 30, 100, 3);
    crate::test_duplex_command::create_duplex_bam(&input, vec![molecule]);
    assert_stdin_input_parity(&input, dir.path(), "duplex", move |i, o| {
        let mut args = vec![
            s("duplex"),
            s("--input"),
            s(i),
            s("--output"),
            s(o),
            s("--min-reads"),
            s("1"),
            s("--compression-level"),
            s("1"),
        ];
        push_threads(&mut args, threads);
        args
    });
}

#[rstest]
#[case("1")]
#[case("2")]
fn test_filter_reads_stdin_once(#[case] threads: &str) {
    use crate::helpers::bam_generator::create_test_reference;
    let dir = TempDir::new().unwrap();
    let reference = create_test_reference(dir.path()).to_str().unwrap().to_string();
    let input = dir.path().join("consensus.bam");
    crate::test_filter_command::write_filter_consensus_bam(&input);
    assert_stdin_input_parity(&input, dir.path(), "filter", move |i, o| {
        // --ref is required when filtering mapped reads (NM/UQ/MD).
        vec![
            s("filter"),
            s("--input"),
            s(i),
            s("--output"),
            s(o),
            s("--ref"),
            reference.clone(),
            s("--min-reads"),
            s("1"),
            s("--max-no-call-fraction"),
            s("1.0"),
            s("--threads"),
            s(threads),
            s("--compression-level"),
            s("1"),
        ]
    });
}

/// clip has a single-threaded fast path (no `--threads`) plus the
/// multi-threaded path; cover both.
#[rstest]
#[case::single_threaded(None)]
#[case::multi_threaded(Some("2"))]
fn test_clip_reads_stdin_once(#[case] threads: Option<&str>) {
    use crate::helpers::bam_generator::{create_paired_umi_family, create_test_reference};
    let dir = TempDir::new().unwrap();
    let reference = create_test_reference(dir.path()).to_str().unwrap().to_string();
    let input = dir.path().join("paired.bam");
    let header = create_minimal_header("chr1", 10000);
    let mut writer = bam::io::Writer::new(fs::File::create(&input).expect("create paired BAM"));
    writer.write_header(&header).expect("write header");
    for raw in &create_paired_umi_family("ACGTACGT", 3, "clip", "ACGTACGTAC", "TGCATGCATG", 30) {
        writer.write_alignment_record(&header, &to_record_buf(raw)).expect("write paired record");
    }
    writer.try_finish().expect("finish BAM");
    assert_stdin_input_parity(&input, dir.path(), "clip", move |i, o| {
        let mut args = vec![
            s("clip"),
            s("--input"),
            s(i),
            s("--output"),
            s(o),
            s("--reference"),
            reference.clone(),
            s("--clip-overlapping-reads"),
            s("--compression-level"),
            s("1"),
        ];
        push_threads(&mut args, threads);
        args
    });
}

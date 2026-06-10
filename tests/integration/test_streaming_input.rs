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
fn compare_bam_records(path1: &PathBuf, path2: &PathBuf) {
    let mut reader1 =
        bam::io::reader::Builder.build_from_path(path1).expect("Failed to open BAM 1");
    let _header1 = reader1.read_header().expect("Failed to read header 1");

    let mut reader2 =
        bam::io::reader::Builder.build_from_path(path2).expect("Failed to open BAM 2");
    let _header2 = reader2.read_header().expect("Failed to read header 2");

    // Compare record counts and content
    let records1: Vec<_> = reader1.records().map(|r| r.expect("Failed to read record")).collect();
    let records2: Vec<_> = reader2.records().map(|r| r.expect("Failed to read record")).collect();

    assert_eq!(records1.len(), records2.len(), "BAM files have different record counts");

    for (i, (r1, r2)) in records1.iter().zip(records2.iter()).enumerate() {
        // Compare key fields
        assert_eq!(r1.name(), r2.name(), "Record {i} has different name");
        assert_eq!(
            r1.sequence().len(),
            r2.sequence().len(),
            "Record {i} has different sequence length"
        );
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

    // Create grouped reads with MI tag
    let mi_tag = Tag::from(fgumi_lib::sam::SamTag::MI);
    let records = create_umi_family("AAAAAAAA", 5, "mol1", "ACGTACGT", 30);

    // Add MI tag to each record (convert to RecordBuf first to enable tag mutation)
    let mi_value = 1;
    for raw in &records {
        let mut rec = to_record_buf(raw);
        rec.data_mut().insert(mi_tag, Value::from(mi_value));
        writer.write_alignment_record(&header, &rec).expect("Failed to write record");
    }

    // Second family with different MI
    let mi_value2 = 2;
    let records2 = create_umi_family("CCCCCCCC", 3, "mol2", "TGCATGCA", 30);
    for raw in &records2 {
        let mut rec = to_record_buf(raw);
        rec.data_mut().insert(mi_tag, Value::from(mi_value2));
        writer.write_alignment_record(&header, &rec).expect("Failed to write record");
    }

    writer.try_finish().expect("Failed to finish BAM");
}

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

use crate::helpers::bam_generator::{create_minimal_header, create_umi_family};

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
            "--raw-tag",
            "RX",
            "--assign-tag",
            "MI",
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
            "--raw-tag",
            "RX",
            "--assign-tag",
            "MI",
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
            "--raw-tag",
            "RX",
            "--assign-tag",
            "MI",
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
        writer.write_alignment_record(&header, record).expect("Failed to write record");
    }

    writer.finish(&header).expect("Failed to finish BAM");
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
    let mi_tag = Tag::new(b'M', b'I');
    let records = create_umi_family("AAAAAAAA", 5, "mol1", "ACGTACGT", 30);

    // Add MI tag to each record
    let mi_value = 1;
    for record in records {
        let mut rec = record.clone();
        // Get existing data and add MI tag
        let data = rec.data_mut();
        data.insert(mi_tag, Value::from(mi_value));
        writer.write_alignment_record(&header, &rec).expect("Failed to write record");
    }

    // Second family with different MI
    let mi_value2 = 2;
    let records2 = create_umi_family("CCCCCCCC", 3, "mol2", "TGCATGCA", 30);
    for record in records2 {
        let mut rec = record.clone();
        let data = rec.data_mut();
        data.insert(mi_tag, Value::from(mi_value2));
        writer.write_alignment_record(&header, &rec).expect("Failed to write record");
    }

    writer.finish(&header).expect("Failed to finish BAM");
}

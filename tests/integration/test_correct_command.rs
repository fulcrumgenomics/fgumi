//! End-to-end CLI tests for the correct command.
//!
//! These tests run the actual `fgumi correct` binary and validate:
//! 1. Basic UMI correction against a whitelist
//! 2. Metrics output
//! 3. Rejected reads output

use noodles::bam;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use std::fs;
use std::path::PathBuf;
use std::process::Command;
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, create_umi_family};

/// Write a BAM with UMI-tagged reads.
fn create_umi_bam(
    path: &PathBuf,
    families: Vec<Vec<noodles::sam::alignment::record_buf::RecordBuf>>,
) {
    let header = create_minimal_header("chr1", 10000);
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");

    for family in families {
        for record in family {
            writer.write_alignment_record(&header, &record).expect("Failed to write record");
        }
    }
    writer.try_finish().expect("Failed to finish BAM");
}

/// Write a UMI whitelist file.
fn create_whitelist(path: &PathBuf, umis: &[&str]) {
    fs::write(path, umis.join("\n")).expect("Failed to write whitelist");
}

/// Test basic UMI correction.
#[test]
fn test_correct_command_basic() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let whitelist = temp_dir.path().join("whitelist.txt");

    // Create reads: 5 with correct UMI "ACGTACGT" and 2 with 1bp error "ACGTACGA"
    let correct_reads = create_umi_family("ACGTACGT", 5, "correct", "AAAAGGGG", 30);
    let error_reads = create_umi_family("ACGTACGA", 2, "error", "AAAAGGGG", 30);
    create_umi_bam(&input_bam, vec![correct_reads, error_reads]);
    create_whitelist(&whitelist, &["ACGTACGT"]);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "correct",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--umi-files",
            whitelist.to_str().unwrap(),
            "--max-mismatches",
            "1",
            "--min-distance",
            "1",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run correct command");

    assert!(status.success(), "Correct command failed");
    assert!(output_bam.exists(), "Output BAM not created");

    // Verify output has records
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let count = reader.records().count();
    assert_eq!(count, 7, "Should have all 7 reads in output");
}

/// Test correct command with metrics output.
#[test]
fn test_correct_command_with_metrics() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let whitelist = temp_dir.path().join("whitelist.txt");
    let metrics = temp_dir.path().join("metrics.tsv");

    let reads = create_umi_family("ACGTACGT", 3, "read", "AAAAGGGG", 30);
    create_umi_bam(&input_bam, vec![reads]);
    create_whitelist(&whitelist, &["ACGTACGT"]);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "correct",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--umi-files",
            whitelist.to_str().unwrap(),
            "--max-mismatches",
            "1",
            "--min-distance",
            "1",
            "--metrics",
            metrics.to_str().unwrap(),
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run correct command");

    assert!(status.success(), "Correct command with metrics failed");
    assert!(metrics.exists(), "Metrics file not created");

    let content = fs::read_to_string(&metrics).unwrap();
    assert!(!content.is_empty(), "Metrics file should not be empty");
}

/// Test correct command with rejects output.
#[test]
fn test_correct_command_with_rejects() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let rejects_bam = temp_dir.path().join("rejects.bam");
    let whitelist = temp_dir.path().join("whitelist.txt");

    // UMI "TTTTTTTT" has edit distance 8 from "ACGTACGT" — won't correct
    let uncorrectable = create_umi_family("TTTTTTTT", 2, "far", "AAAAGGGG", 30);
    let correctable = create_umi_family("ACGTACGT", 3, "exact", "AAAAGGGG", 30);
    create_umi_bam(&input_bam, vec![correctable, uncorrectable]);
    create_whitelist(&whitelist, &["ACGTACGT"]);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "correct",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--umi-files",
            whitelist.to_str().unwrap(),
            "--max-mismatches",
            "1",
            "--min-distance",
            "1",
            "--rejects",
            rejects_bam.to_str().unwrap(),
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run correct command");

    assert!(status.success(), "Correct command with rejects failed");
    assert!(rejects_bam.exists(), "Rejects BAM not created");
}

/// Test correct command with --allow-c-to-t flag: C→T converted UMIs should match whitelist.
#[test]
fn test_correct_command_allow_c_to_t() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let rejects_bam = temp_dir.path().join("rejects.bam");
    let whitelist = temp_dir.path().join("whitelist.txt");

    // Whitelist UMI has C at positions 1,3: "ACGCACGC"
    // Observed UMI has C→T at those positions: "ATGTATGT"
    // Without --allow-c-to-t this would be 4 mismatches; with --allow-c-to-t it should be 0
    let exact_reads = create_umi_family("ACGCACGC", 3, "exact", "AAAAGGGG", 30);
    let converted_reads = create_umi_family("ATGTATGT", 2, "converted", "AAAAGGGG", 30);
    create_umi_bam(&input_bam, vec![exact_reads, converted_reads]);
    create_whitelist(&whitelist, &["ACGCACGC"]);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "correct",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--umi-files",
            whitelist.to_str().unwrap(),
            "--max-mismatches",
            "0",
            "--min-distance",
            "1",
            "--allow-c-to-t",
            "--rejects",
            rejects_bam.to_str().unwrap(),
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run correct command");

    assert!(status.success(), "Correct command with --allow-c-to-t failed");
    assert!(output_bam.exists(), "Output BAM not created");

    // All 5 reads should be in output (exact match + C→T match both pass at max_mismatches=0)
    // and all reads should carry the canonical whitelist RX tag
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let header = reader.read_header().unwrap();
    let rx_tag = [b'R', b'X'];
    let mut output_count = 0;
    for result in reader.record_bufs(&header) {
        let record = result.expect("Failed to read record");
        output_count += 1;
        let rx = record.data().get(&rx_tag).expect("Record should have RX tag");
        let rx_str = match rx {
            noodles::sam::alignment::record_buf::data::field::Value::String(s) => s.to_string(),
            _ => panic!("RX tag should be a string"),
        };
        assert_eq!(
            rx_str, "ACGCACGC",
            "All corrected reads should carry the canonical whitelist RX"
        );
    }
    assert_eq!(output_count, 5, "Should have all 5 reads in output with --allow-c-to-t");

    // Rejects should be empty
    let mut reject_reader = bam::io::Reader::new(fs::File::open(&rejects_bam).unwrap());
    let _header = reject_reader.read_header().unwrap();
    let reject_count = reject_reader.records().count();
    assert_eq!(reject_count, 0, "Should have 0 rejects with --allow-c-to-t");
}

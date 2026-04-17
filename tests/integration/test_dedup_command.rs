//! End-to-end CLI tests for the dedup command.
//!
//! These tests run the actual `fgumi dedup` binary and validate:
//! 1. Basic duplicate marking
//! 2. Metrics output
//! 3. Remove duplicates mode

use fgumi_raw_bam::{RawRecord, SamBuilder, flags};
use noodles::bam;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use std::fs;
use std::path::PathBuf;
use std::process::Command;
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, to_record_buf};

/// Create a template-coordinate sorted BAM with UMI-tagged reads.
///
/// Template-coordinate sort groups reads by position, then by name within each position.
/// The header must have SO:unsorted GO:query SS:template-coordinate tags.
fn create_sorted_bam(path: &PathBuf, records: Vec<RawRecord>) {
    let header = create_minimal_header("chr1", 10000);
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");

    for record in records {
        writer
            .write_alignment_record(&header, &to_record_buf(&record))
            .expect("Failed to write record");
    }
    writer.try_finish().expect("Failed to finish BAM");
}

/// Create a group of paired-end reads at the same position with the same UMI
/// (simulating PCR duplicates).
fn create_duplicate_group(base_name: &str, umi: &str, count: usize, start: i32) -> Vec<RawRecord> {
    let mut records = Vec::new();
    for i in 0..count {
        let name = format!("{base_name}_{i}");

        let r1 = {
            let mut b = SamBuilder::new();
            b.read_name(name.as_bytes())
                .sequence(b"ACGTACGT")
                .qualities(&[30; 8])
                .flags(flags::PAIRED | flags::FIRST_SEGMENT)
                .ref_id(0)
                .pos(start - 1)
                .mapq(60)
                .cigar_ops(&[8 << 4]) // 8M
                .mate_ref_id(0)
                .mate_pos(start + 99)
                .template_length(108)
                .add_string_tag(b"RX", umi.as_bytes())
                .add_string_tag(b"MC", b"8M");
            b.build()
        };

        let r2 = {
            let mut b = SamBuilder::new();
            b.read_name(name.as_bytes())
                .sequence(b"ACGTACGT")
                .qualities(&[30; 8])
                .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE)
                .ref_id(0)
                .pos(start + 99)
                .mapq(60)
                .cigar_ops(&[8 << 4]) // 8M
                .mate_ref_id(0)
                .mate_pos(start - 1)
                .template_length(-108)
                .add_string_tag(b"RX", umi.as_bytes())
                .add_string_tag(b"MC", b"8M");
            b.build()
        };

        records.push(r1);
        records.push(r2);
    }
    records
}

/// Test basic dedup command (mark duplicates).
#[test]
fn test_dedup_command_basic() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    // Create 3 duplicate pairs at position 100 with same UMI, and 2 at position 500
    let mut records = create_duplicate_group("dup1", "ACGTACGT", 3, 100);
    records.extend(create_duplicate_group("dup2", "TGCATGCA", 2, 500));
    create_sorted_bam(&input_bam, records);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "dedup",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--strategy",
            "identity",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run dedup command");

    assert!(status.success(), "Dedup command failed");
    assert!(output_bam.exists(), "Output BAM not created");

    // All reads should be present (duplicates are marked, not removed)
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let count = reader.records().count();
    assert_eq!(count, 10, "All reads should be in output (marked, not removed)");
}

/// Test dedup command with metrics output.
#[test]
fn test_dedup_command_with_metrics() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let metrics_path = temp_dir.path().join("metrics.txt");

    let records = create_duplicate_group("dup1", "ACGTACGT", 3, 100);
    create_sorted_bam(&input_bam, records);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "dedup",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--strategy",
            "identity",
            "--metrics",
            metrics_path.to_str().unwrap(),
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run dedup command");

    assert!(status.success(), "Dedup command with metrics failed");
    assert!(metrics_path.exists(), "Metrics file not created");
}

/// Test dedup command with remove-duplicates flag.
#[test]
fn test_dedup_command_remove_duplicates() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    // 3 duplicate pairs → should keep 1 pair (2 records), remove 2 pairs (4 records)
    let records = create_duplicate_group("dup1", "ACGTACGT", 3, 100);
    create_sorted_bam(&input_bam, records);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "dedup",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--strategy",
            "identity",
            "--remove-duplicates",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run dedup command");

    assert!(status.success(), "Dedup command with --remove-duplicates failed");
    assert!(output_bam.exists(), "Output BAM not created");

    // With remove-duplicates, only the best pair should remain
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let count = reader.records().count();
    assert!(count < 6, "Remove-duplicates should produce fewer reads than input");
    assert!(count >= 2, "Should keep at least one pair");
}

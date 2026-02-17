//! End-to-end CLI tests for the filter command.
//!
//! These tests run the actual `fgumi filter` binary and validate:
//! 1. Basic consensus read filtering
//! 2. Rejected reads output
//! 3. Statistics output

use fgumi_lib::sam::builder::{ConsensusTagsBuilder, RecordBuilder};
use noodles::bam;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record_buf::RecordBuf;
use std::fs;
use std::path::Path;
use std::process::Command;
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, create_test_reference};

/// Create a consensus BAM with records that have consensus tags.
fn create_consensus_bam(path: &Path, records: Vec<RecordBuf>) {
    let header = create_minimal_header("chr1", 10000);
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");

    for record in records {
        writer.write_alignment_record(&header, &record).expect("Failed to write record");
    }
    writer.finish(&header).expect("Failed to finish BAM");
}

/// Test basic filter command with passing reads.
#[test]
fn test_filter_command_basic() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let ref_path = create_test_reference(temp_dir.path());

    // Create consensus reads with good quality and per-base tags (cd/ce).
    let r1 = RecordBuilder::new()
        .name("cons1")
        .sequence("ACGTACGT")
        .qualities(&[35; 8])
        .reference_sequence_id(0)
        .alignment_start(100)
        .mapping_quality(60)
        .cigar("8M")
        .consensus_tags(
            ConsensusTagsBuilder::new().per_base_depths(&[10; 8]).per_base_errors(&[0; 8]),
        )
        .build();

    let r2 = RecordBuilder::new()
        .name("cons2")
        .sequence("ACGTACGT")
        .qualities(&[35; 8])
        .reference_sequence_id(0)
        .alignment_start(200)
        .mapping_quality(60)
        .cigar("8M")
        .consensus_tags(
            ConsensusTagsBuilder::new().per_base_depths(&[5; 8]).per_base_errors(&[0; 8]),
        )
        .build();

    create_consensus_bam(&input_bam, vec![r1, r2]);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "filter",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--ref",
            ref_path.to_str().unwrap(),
            "--min-reads",
            "1",
            "--max-no-call-fraction",
            "1.0",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run filter command");

    assert!(status.success(), "Filter command failed");
    assert!(output_bam.exists(), "Output BAM not created");

    // Verify output has records (reads may have masked bases but should still be present)
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let count = reader.records().count();
    assert_eq!(count, 2, "Both reads should be present in output");
}

/// Test filter command with reads that fail due to low per-base depth.
///
/// The filter masks bases where the per-base depth (cd tag) is below --min-reads.
/// If enough bases are masked, the no-call fraction exceeds the threshold and the
/// read is rejected.
#[test]
fn test_filter_command_rejects_low_depth() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let rejects_bam = temp_dir.path().join("rejects.bam");
    let ref_path = create_test_reference(temp_dir.path());

    // Good read: per-base depth 10 (all above min-reads=3), no bases masked
    let good = RecordBuilder::new()
        .name("good")
        .sequence("ACGTACGT")
        .qualities(&[35; 8])
        .reference_sequence_id(0)
        .alignment_start(100)
        .mapping_quality(60)
        .cigar("8M")
        .consensus_tags(
            ConsensusTagsBuilder::new().per_base_depths(&[10; 8]).per_base_errors(&[0; 8]),
        )
        .build();

    // Low-depth read: per-base depth 1 (all below min-reads=3), all bases masked
    let low_depth = RecordBuilder::new()
        .name("low_depth")
        .sequence("ACGTACGT")
        .qualities(&[35; 8])
        .reference_sequence_id(0)
        .alignment_start(200)
        .mapping_quality(60)
        .cigar("8M")
        .consensus_tags(
            ConsensusTagsBuilder::new().per_base_depths(&[1; 8]).per_base_errors(&[0; 8]),
        )
        .build();

    create_consensus_bam(&input_bam, vec![good, low_depth]);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "filter",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--ref",
            ref_path.to_str().unwrap(),
            "--min-reads",
            "3",
            "--rejects",
            rejects_bam.to_str().unwrap(),
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run filter command");

    assert!(status.success(), "Filter command with rejects failed");
    assert!(output_bam.exists(), "Output BAM not created");
    assert!(rejects_bam.exists(), "Rejects BAM not created");

    // Verify the good read passed and the low-depth read was rejected
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let output_count = reader.records().count();
    assert_eq!(output_count, 1, "Only the good read should pass filtering");

    let mut reject_reader = bam::io::Reader::new(fs::File::open(&rejects_bam).unwrap());
    let _header = reject_reader.read_header().unwrap();
    let reject_count = reject_reader.records().count();
    assert_eq!(reject_count, 1, "The low-depth read should be rejected");
}

/// Test filter command with statistics output.
#[test]
fn test_filter_command_with_stats() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let stats_path = temp_dir.path().join("stats.txt");
    let ref_path = create_test_reference(temp_dir.path());

    let record = RecordBuilder::new()
        .name("cons1")
        .sequence("ACGTACGT")
        .qualities(&[35; 8])
        .reference_sequence_id(0)
        .alignment_start(100)
        .mapping_quality(60)
        .cigar("8M")
        .consensus_tags(ConsensusTagsBuilder::new().depth_max(10).depth_min(8).error_rate(0.005))
        .build();

    create_consensus_bam(&input_bam, vec![record]);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "filter",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--ref",
            ref_path.to_str().unwrap(),
            "--min-reads",
            "1",
            "--stats",
            stats_path.to_str().unwrap(),
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run filter command");

    assert!(status.success(), "Filter command with stats failed");
    assert!(stats_path.exists(), "Stats file not created");
    let content = fs::read_to_string(&stats_path).expect("Failed to read stats file");
    assert!(!content.trim().is_empty(), "Stats file should not be empty");
    assert!(content.contains("total_reads"), "Stats should contain total_reads");
}

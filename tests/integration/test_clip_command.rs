//! End-to-end CLI tests for the clip command.
//!
//! These tests run the actual `fgumi clip` binary and validate:
//! 1. Basic read clipping
//! 2. Fixed clipping (5' and 3' ends)
//! 3. Metrics output

use fgumi_lib::sam::builder::RecordBuilder;
use noodles::bam;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record_buf::RecordBuf;
use std::fs;
use std::path::Path;
use std::process::Command;
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, create_test_reference};

/// Create a BAM with paired reads.
fn create_paired_bam(path: &Path, pairs: Vec<(RecordBuf, RecordBuf)>) {
    let header = create_minimal_header("chr1", 10000);
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");

    for (r1, r2) in pairs {
        writer.write_alignment_record(&header, &r1).expect("Failed to write R1");
        writer.write_alignment_record(&header, &r2).expect("Failed to write R2");
    }
    writer.finish(&header).expect("Failed to finish BAM");
}

/// Test basic clip command with fixed-end clipping.
#[test]
fn test_clip_command_basic() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let ref_path = create_test_reference(temp_dir.path());

    // Create a paired-end read
    let r1 = RecordBuilder::new()
        .name("read1")
        .sequence("ACGTACGT")
        .qualities(&[30; 8])
        .paired(true)
        .first_segment(true)
        .reference_sequence_id(0)
        .alignment_start(100)
        .mapping_quality(60)
        .cigar("8M")
        .mate_reference_sequence_id(0)
        .mate_alignment_start(104)
        .template_length(12)
        .build();

    let r2 = RecordBuilder::new()
        .name("read1")
        .sequence("ACGTACGT")
        .qualities(&[30; 8])
        .paired(true)
        .first_segment(false)
        .reverse_complement(true)
        .reference_sequence_id(0)
        .alignment_start(104)
        .mapping_quality(60)
        .cigar("8M")
        .mate_reference_sequence_id(0)
        .mate_alignment_start(100)
        .template_length(-12)
        .build();

    create_paired_bam(&input_bam, vec![(r1, r2)]);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "clip",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--reference",
            ref_path.to_str().unwrap(),
            "--read-one-five-prime",
            "1",
            "--read-one-three-prime",
            "1",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run clip command");

    assert!(status.success(), "Clip command failed");
    assert!(output_bam.exists(), "Output BAM not created");

    // Verify output has the reads
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let count = reader.records().count();
    assert_eq!(count, 2, "Should have both reads in output");
}

/// Test clip command with metrics output.
#[test]
fn test_clip_command_with_metrics() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let metrics_path = temp_dir.path().join("metrics.txt");
    let ref_path = create_test_reference(temp_dir.path());

    let r1 = RecordBuilder::new()
        .name("read1")
        .sequence("ACGTACGT")
        .qualities(&[30; 8])
        .paired(true)
        .first_segment(true)
        .reference_sequence_id(0)
        .alignment_start(100)
        .mapping_quality(60)
        .cigar("8M")
        .mate_reference_sequence_id(0)
        .mate_alignment_start(200)
        .template_length(108)
        .build();

    let r2 = RecordBuilder::new()
        .name("read1")
        .sequence("ACGTACGT")
        .qualities(&[30; 8])
        .paired(true)
        .first_segment(false)
        .reverse_complement(true)
        .reference_sequence_id(0)
        .alignment_start(200)
        .mapping_quality(60)
        .cigar("8M")
        .mate_reference_sequence_id(0)
        .mate_alignment_start(100)
        .template_length(-108)
        .build();

    create_paired_bam(&input_bam, vec![(r1, r2)]);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "clip",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--reference",
            ref_path.to_str().unwrap(),
            "--read-one-five-prime",
            "2",
            "--metrics",
            metrics_path.to_str().unwrap(),
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run clip command");

    assert!(status.success(), "Clip command with metrics failed");
    assert!(metrics_path.exists(), "Metrics file not created");
}

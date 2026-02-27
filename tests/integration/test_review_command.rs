//! End-to-end CLI tests for the review command.
//!
//! These tests run the actual `fgumi review` binary and validate:
//! 1. Missing required arguments produce a clear error
//! 2. Basic execution with minimal valid inputs

use fgumi_lib::sam::builder::RecordBuilder;
use noodles::bam;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use std::fs;
use std::io::Write;
use std::path::Path;
use std::process::Command;
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_coordinate_sorted_header, create_test_reference};

/// Create a coordinate-sorted, indexed BAM file.
fn create_indexed_bam(path: &Path, records: &[noodles::sam::alignment::record_buf::RecordBuf]) {
    let header = create_coordinate_sorted_header("chr1", 10000);

    // Write BAM and drop to flush
    {
        let mut writer =
            bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
        writer.write_header(&header).expect("Failed to write header");
        for record in records {
            writer.write_alignment_record(&header, record).expect("Failed to write record");
        }
    }

    // Create BAI index
    let index = bam::fs::index(path).expect("Failed to create BAM index");
    let index_path = path.with_extension("bam.bai");
    let mut index_writer =
        noodles::bam::bai::io::Writer::new(fs::File::create(&index_path).unwrap());
    index_writer.write_index(&index).unwrap();
}

/// Create a minimal VCF file with one variant.
fn create_test_vcf(path: &Path) {
    let mut vcf = fs::File::create(path).unwrap();
    writeln!(vcf, "##fileformat=VCFv4.3").unwrap();
    writeln!(vcf, "##contig=<ID=chr1,length=10000>").unwrap();
    writeln!(vcf, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
    writeln!(vcf, "chr1\t100\t.\tA\tC\t.\tPASS\t.").unwrap();
    vcf.flush().unwrap();
}

/// Test that missing required arguments produce a non-zero exit code.
#[test]
fn test_review_missing_required_args() {
    // No arguments at all — should fail with clap error
    let output = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args(["review"])
        .output()
        .expect("Failed to run review command");

    assert!(!output.status.success(), "Review should fail without required args");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("required") || stderr.contains("Usage") || stderr.contains("error"),
        "Stderr should mention missing required arguments, got: {stderr}"
    );
}

/// Test that review fails clearly when required file inputs are missing.
#[test]
fn test_review_missing_input_files() {
    let temp_dir = TempDir::new().unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "review",
            "--input",
            "/nonexistent/variants.vcf",
            "--consensus-bam",
            "/nonexistent/consensus.bam",
            "--grouped-bam",
            "/nonexistent/grouped.bam",
            "--ref",
            "/nonexistent/ref.fa",
            "--output",
            temp_dir.path().join("output").to_str().unwrap(),
        ])
        .output()
        .expect("Failed to run review command");

    assert!(!output.status.success(), "Review should fail for nonexistent input files");
}

/// Test basic review execution with valid inputs.
///
/// Creates minimal VCF, consensus BAM, grouped BAM, and reference files.
/// Verifies the command completes and produces output files.
#[test]
fn test_review_basic_execution() {
    let temp_dir = TempDir::new().unwrap();
    let vcf_path = temp_dir.path().join("variants.vcf");
    let consensus_bam = temp_dir.path().join("consensus.bam");
    let grouped_bam = temp_dir.path().join("grouped.bam");
    let ref_path = create_test_reference(temp_dir.path());
    let output_prefix = temp_dir.path().join("review_out");

    // Create VCF with one variant at chr1:100
    create_test_vcf(&vcf_path);

    // Create a consensus BAM with a read spanning the variant position,
    // with MI tag for molecule identifier
    let consensus_records = vec![
        RecordBuilder::new()
            .name("cons1")
            .sequence("ACGTACGT")
            .qualities(&[30; 8])
            .reference_sequence_id(0)
            .alignment_start(97)
            .mapping_quality(60)
            .cigar("8M")
            .tag("MI", "1")
            .build(),
    ];
    create_indexed_bam(&consensus_bam, &consensus_records);

    // Create a grouped BAM with raw reads (MI tag matches consensus)
    let grouped_records = vec![
        RecordBuilder::new()
            .name("raw1")
            .sequence("ACGTACGT")
            .qualities(&[30; 8])
            .reference_sequence_id(0)
            .alignment_start(97)
            .mapping_quality(60)
            .cigar("8M")
            .tag("MI", "1")
            .build(),
        RecordBuilder::new()
            .name("raw2")
            .sequence("ACGTACGT")
            .qualities(&[30; 8])
            .reference_sequence_id(0)
            .alignment_start(97)
            .mapping_quality(60)
            .cigar("8M")
            .tag("MI", "1")
            .build(),
    ];
    create_indexed_bam(&grouped_bam, &grouped_records);

    let output = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "review",
            "--input",
            vcf_path.to_str().unwrap(),
            "--consensus-bam",
            consensus_bam.to_str().unwrap(),
            "--grouped-bam",
            grouped_bam.to_str().unwrap(),
            "--ref",
            ref_path.to_str().unwrap(),
            "--output",
            output_prefix.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to run review command");

    assert!(
        output.status.success(),
        "Review command failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Verify output files were created
    let consensus_out = output_prefix.with_extension("consensus.bam");
    let grouped_out = output_prefix.with_extension("grouped.bam");
    assert!(consensus_out.exists(), "Consensus output BAM should exist");
    assert!(grouped_out.exists(), "Grouped output BAM should exist");
}

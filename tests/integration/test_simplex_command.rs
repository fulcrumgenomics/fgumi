//! End-to-end CLI tests for the simplex command.
//!
//! These tests run the actual `fgumi simplex` binary and validate:
//! 1. Basic simplex consensus calling from grouped reads
//! 2. Statistics output
//! 3. Rejected reads output

use noodles::bam;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record_buf::RecordBuf;
use std::fs;
use std::path::Path;
use std::process::Command;
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, create_umi_family};

/// Write grouped BAM file (reads grouped by MI tag).
fn create_grouped_bam(path: &Path, families: Vec<(&str, Vec<RecordBuf>)>) {
    let header = create_minimal_header("chr1", 10000);
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");

    for (mi, records) in families {
        for mut record in records {
            // Add MI tag for molecule identity (simplex expects grouped reads with MI)
            use noodles::sam::alignment::record::data::field::Tag;
            use noodles::sam::alignment::record_buf::data::field::Value;
            let mi_tag = Tag::new(b'M', b'I');
            record.data_mut().insert(mi_tag, Value::from(mi));
            writer.write_alignment_record(&header, &record).expect("Failed to write record");
        }
    }
    writer.finish(&header).expect("Failed to finish BAM");
}

/// Test basic simplex consensus calling.
#[test]
fn test_simplex_command_basic_consensus() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    // Create two families: 5 reads each
    let family1 = create_umi_family("ACGT", 5, "fam1", "ACGTACGT", 30);
    let family2 = create_umi_family("TGCA", 5, "fam2", "TTTTAAAA", 30);
    create_grouped_bam(&input_bam, vec![("1", family1), ("2", family2)]);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "simplex",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--min-reads",
            "2",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run simplex command");

    assert!(status.success(), "Simplex command failed");
    assert!(output_bam.exists(), "Output BAM not created");

    // Read output and verify consensus reads were produced
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let mut count = 0;
    for result in reader.records() {
        let record = result.expect("Failed to read record");
        // Verify consensus tags exist
        let cd_tag = [b'c', b'D'];
        assert!(record.data().get(&cd_tag).is_some(), "Consensus should have cD tag");
        count += 1;
    }
    assert!(count > 0, "Should have produced consensus reads");
}

/// Test simplex command with statistics output.
#[test]
fn test_simplex_command_with_stats() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let stats_path = temp_dir.path().join("stats.txt");

    let family = create_umi_family("ACGT", 5, "fam1", "ACGTACGT", 30);
    create_grouped_bam(&input_bam, vec![("1", family)]);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "simplex",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--min-reads",
            "1",
            "--stats",
            stats_path.to_str().unwrap(),
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run simplex command");

    assert!(status.success(), "Simplex command with stats failed");
    assert!(stats_path.exists(), "Stats file not created");
}

/// Test simplex command with rejects output.
#[test]
fn test_simplex_command_with_rejects() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let rejects_bam = temp_dir.path().join("rejects.bam");

    // Create one family with 5 reads (passes min-reads=2) and one with 1 read (fails)
    let family1 = create_umi_family("ACGT", 5, "fam1", "ACGTACGT", 30);
    let family2 = create_umi_family("TGCA", 1, "fam2", "TTTTAAAA", 30);
    create_grouped_bam(&input_bam, vec![("1", family1), ("2", family2)]);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "simplex",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--min-reads",
            "2",
            "--rejects",
            rejects_bam.to_str().unwrap(),
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run simplex command");

    assert!(status.success(), "Simplex command with rejects failed");
    assert!(rejects_bam.exists(), "Rejects BAM not created");
}

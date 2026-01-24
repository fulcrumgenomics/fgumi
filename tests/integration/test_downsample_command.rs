//! Integration tests for the downsample command.

use fgumi_lib::sam::builder::RecordBuilder;
use noodles::bam;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value;
use std::fs;
use std::path::PathBuf;
use std::process::Command;
use tempfile::TempDir;

use crate::helpers::bam_generator::create_minimal_header;

/// MI tag constant
const MI_TAG: Tag = Tag::new(b'M', b'I');

/// Create a grouped BAM file with MI tags (simulating output from group).
fn create_grouped_bam(path: &PathBuf, families: Vec<(&str, usize)>) {
    let header = create_minimal_header("chr1", 10000);

    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));

    writer.write_header(&header).expect("Failed to write header");

    // Create records grouped by MI tag
    let mut read_idx = 0;
    for (mi, count) in families {
        for _ in 0..count {
            let record = RecordBuilder::new()
                .name(&format!("read_{read_idx}"))
                .sequence("ACGT") // Minimal sequence
                .reference_sequence_id(0)
                .alignment_start(100)
                .mapping_quality(60)
                .tag("MI", mi)
                .build();

            writer.write_alignment_record(&header, &record).expect("Failed to write record");
            read_idx += 1;
        }
    }

    writer.finish(&header).expect("Failed to finish BAM");
}

/// Read records from a BAM file.
fn read_bam_records(path: &PathBuf) -> Vec<noodles::sam::alignment::RecordBuf> {
    let mut reader = bam::io::reader::Builder.build_from_path(path).expect("Failed to open BAM");
    let header = reader.read_header().expect("Failed to read header");

    reader.record_bufs(&header).map(|r| r.expect("Failed to read record")).collect()
}

/// Count the number of unique MI values in a BAM file.
fn count_unique_mis(path: &PathBuf) -> usize {
    let records = read_bam_records(path);
    let mis: std::collections::HashSet<String> = records
        .iter()
        .filter_map(|r| {
            r.data().get(&MI_TAG).map(|v| {
                if let Value::String(s) = v {
                    s.to_string()
                } else {
                    panic!("MI tag is not a string")
                }
            })
        })
        .collect();
    mis.len()
}

/// Test basic downsampling functionality.
#[test]
fn test_downsample_basic() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    // Create input BAM with 10 families of 5 reads each
    let families: Vec<(&str, usize)> =
        (0..10).map(|i| (Box::leak(format!("{i}").into_boxed_str()) as &str, 5)).collect();
    create_grouped_bam(&input_bam, families);

    // Run downsample with fraction=0.5 and seed for reproducibility
    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "downsample",
            "-i",
            input_bam.to_str().unwrap(),
            "-o",
            output_bam.to_str().unwrap(),
            "-f",
            "0.5",
            "--seed",
            "42",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run downsample command");

    assert!(status.success(), "Downsample command failed");
    assert!(output_bam.exists(), "Output BAM not created");

    // Verify that some families were kept (with seed 42 and fraction 0.5, should be around 5)
    let output_records = read_bam_records(&output_bam);
    assert!(!output_records.is_empty(), "Output should have some records");
    assert!(output_records.len() < 50, "Output should have fewer records than input (50)");

    let output_families = count_unique_mis(&output_bam);
    assert!(output_families > 0, "Should have kept some families");
    assert!(output_families < 10, "Should have fewer than 10 families");
}

/// Test downsampling with rejects output.
#[test]
fn test_downsample_with_rejects() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let rejects_bam = temp_dir.path().join("rejects.bam");

    // Create input BAM with 10 families
    let families: Vec<(&str, usize)> =
        (0..10).map(|i| (Box::leak(format!("{i}").into_boxed_str()) as &str, 5)).collect();
    create_grouped_bam(&input_bam, families);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "downsample",
            "-i",
            input_bam.to_str().unwrap(),
            "-o",
            output_bam.to_str().unwrap(),
            "-f",
            "0.5",
            "--seed",
            "42",
            "--rejects",
            rejects_bam.to_str().unwrap(),
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run downsample command");

    assert!(status.success(), "Downsample command failed");
    assert!(output_bam.exists(), "Output BAM not created");
    assert!(rejects_bam.exists(), "Rejects BAM not created");

    // Kept + rejected should equal input
    let output_count = read_bam_records(&output_bam).len();
    let rejects_count = read_bam_records(&rejects_bam).len();

    assert_eq!(output_count + rejects_count, 50, "Total records should be preserved");
}

/// Test determinism with same seed.
#[test]
fn test_downsample_deterministic() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output1_bam = temp_dir.path().join("output1.bam");
    let output2_bam = temp_dir.path().join("output2.bam");

    // Create input BAM
    let families: Vec<(&str, usize)> =
        (0..20).map(|i| (Box::leak(format!("{i}").into_boxed_str()) as &str, 3)).collect();
    create_grouped_bam(&input_bam, families);

    // Run twice with same seed
    for output in [&output1_bam, &output2_bam] {
        let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
            .args([
                "downsample",
                "-i",
                input_bam.to_str().unwrap(),
                "-o",
                output.to_str().unwrap(),
                "-f",
                "0.3",
                "--seed",
                "12345",
                "--compression-level",
                "1",
            ])
            .status()
            .expect("Failed to run downsample command");

        assert!(status.success(), "Downsample command failed");
    }

    // Both outputs should be identical
    let records1 = read_bam_records(&output1_bam);
    let records2 = read_bam_records(&output2_bam);

    assert_eq!(records1.len(), records2.len(), "Same seed should produce same count");

    // Check that the same families were selected
    let mi_set1: std::collections::HashSet<_> = records1
        .iter()
        .filter_map(|r| {
            r.data().get(&MI_TAG).map(|v| {
                if let Value::String(s) = v {
                    s.to_string()
                } else {
                    panic!("MI tag is not a string")
                }
            })
        })
        .collect();

    let mi_set2: std::collections::HashSet<_> = records2
        .iter()
        .filter_map(|r| {
            r.data().get(&MI_TAG).map(|v| {
                if let Value::String(s) = v {
                    s.to_string()
                } else {
                    panic!("MI tag is not a string")
                }
            })
        })
        .collect();

    assert_eq!(mi_set1, mi_set2, "Same seed should select same families");
}

/// Test histogram output.
#[test]
fn test_downsample_histograms() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let hist_kept = temp_dir.path().join("hist_kept.txt");
    let hist_rejected = temp_dir.path().join("hist_rejected.txt");

    // Create input BAM with varied family sizes
    let families = vec![("0", 1), ("1", 2), ("2", 3), ("3", 4), ("4", 5)];
    create_grouped_bam(&input_bam, families);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "downsample",
            "-i",
            input_bam.to_str().unwrap(),
            "-o",
            output_bam.to_str().unwrap(),
            "-f",
            "0.5",
            "--seed",
            "42",
            "--histogram-kept",
            hist_kept.to_str().unwrap(),
            "--histogram-rejected",
            hist_rejected.to_str().unwrap(),
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run downsample command");

    assert!(status.success(), "Downsample command failed");
    assert!(hist_kept.exists(), "Kept histogram not created");
    assert!(hist_rejected.exists(), "Rejected histogram not created");

    // Check histogram format
    let kept_contents = fs::read_to_string(&hist_kept).unwrap();
    assert!(kept_contents.contains("family_size\tcount"), "Histogram should have header");

    let rejected_contents = fs::read_to_string(&hist_rejected).unwrap();
    assert!(rejected_contents.contains("family_size\tcount"), "Histogram should have header");
}

/// Test with fraction=1.0 (keep all).
#[test]
fn test_downsample_keep_all() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    let families = vec![("0", 5), ("1", 3), ("2", 7)];
    create_grouped_bam(&input_bam, families);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "downsample",
            "-i",
            input_bam.to_str().unwrap(),
            "-o",
            output_bam.to_str().unwrap(),
            "-f",
            "1.0",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run downsample command");

    assert!(status.success(), "Downsample command failed");

    // All records should be kept
    let output_records = read_bam_records(&output_bam);
    assert_eq!(output_records.len(), 15, "All 15 records should be kept with fraction=1.0");
}

/// Test error handling for invalid fraction.
#[test]
fn test_downsample_invalid_fraction() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    create_grouped_bam(&input_bam, vec![("0", 5)]);

    // Test fraction = 0.0
    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "downsample",
            "-i",
            input_bam.to_str().unwrap(),
            "-o",
            output_bam.to_str().unwrap(),
            "-f",
            "0.0",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run downsample command");

    assert!(!status.success(), "Fraction=0.0 should fail");

    // Test fraction > 1.0
    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "downsample",
            "-i",
            input_bam.to_str().unwrap(),
            "-o",
            output_bam.to_str().unwrap(),
            "-f",
            "1.5",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run downsample command");

    assert!(!status.success(), "Fraction>1.0 should fail");
}

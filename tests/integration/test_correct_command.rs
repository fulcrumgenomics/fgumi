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

use crate::helpers::bam_generator::{create_minimal_header, create_umi_family, to_record_buf};
use fgumi_raw_bam::RawRecord;

/// Write a BAM with UMI-tagged reads.
fn create_umi_bam(path: &PathBuf, families: Vec<Vec<RawRecord>>) {
    let header = create_minimal_header("chr1", 10000);
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");

    for family in families {
        for record in &family {
            writer
                .write_alignment_record(&header, &to_record_buf(record))
                .expect("Failed to write record");
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

/// Exercises the multi-threaded rejects-streaming path end-to-end and asserts:
/// 1. The rejects BAM is a valid BGZF stream with the terminating EOF block.
/// 2. The `@HD` header reports `SO:unsorted` because rejects are emitted in
///    mutex-acquisition order rather than input order.
/// 3. Every uncorrectable input record appears exactly once in the rejects
///    BAM — the writer does not drop records under worker contention and does
///    not emit duplicates.
#[test]
fn test_correct_command_rejects_streaming_threaded_integrity() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let rejects_bam = temp_dir.path().join("rejects.bam");
    let whitelist = temp_dir.path().join("whitelist.txt");

    // Each record has a unique QNAME, so each record is its own template.
    // `CorrectUmis` uses a 1000-template batch, so size the uncorrectable
    // families so that `far_a + far_b + far_c` exceeds that boundary and
    // multiple batches flow through the 4-thread pool concurrently, letting
    // more than one worker race to flush rejects.
    let far_family_size: u32 = 400;
    let corr_small = create_umi_family("ACGTACGT", 3, "c1_exact", "AAAAGGGG", 30);
    let corr_big = create_umi_family("ACGTACGT", 30, "c2_exact", "AAAAGGGG", 30);
    let far_a = create_umi_family("TTTTTTTT", far_family_size as usize, "far_a", "AAAAGGGG", 30);
    let far_b = create_umi_family("GGGGGGGG", far_family_size as usize, "far_b", "AAAAGGGG", 30);
    let far_c = create_umi_family("CCCCCCCC", far_family_size as usize, "far_c", "AAAAGGGG", 30);

    // Expected reject names mirror `create_umi_family`'s "{base_name}_{i}"
    // convention for the three uncorrectable families.
    let mut expected_names: std::collections::HashSet<String> = std::collections::HashSet::new();
    for base in ["far_a", "far_b", "far_c"] {
        for i in 0..far_family_size {
            expected_names.insert(format!("{base}_{i}"));
        }
    }
    let expected_rejects = expected_names.len();

    create_umi_bam(&input_bam, vec![corr_small, corr_big, far_a, far_b, far_c]);
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
            "--threads",
            "4",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run correct command");

    assert!(status.success(), "Correct command with threaded rejects failed");
    assert!(rejects_bam.exists(), "Rejects BAM not created");

    crate::helpers::assertions::assert_has_bgzf_eof(&rejects_bam);
    crate::helpers::assertions::assert_header_unsorted(&rejects_bam);

    let mut reader = bam::io::Reader::new(fs::File::open(&rejects_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let mut seen: std::collections::HashMap<String, u32> = std::collections::HashMap::new();
    let mut total: usize = 0;
    for result in reader.records() {
        let record = result.expect("Failed to read reject record");
        let name = record.name().expect("reject record missing read name").to_string();
        *seen.entry(name).or_insert(0) += 1;
        total += 1;
    }

    assert_eq!(
        total, expected_rejects,
        "rejects BAM should contain one record per uncorrectable input read",
    );
    assert_eq!(seen.len(), expected_rejects, "unexpected reject-name set size");
    for name in &expected_names {
        assert_eq!(
            seen.remove(name),
            Some(1),
            "expected uncorrectable record {name} exactly once in the rejects BAM",
        );
    }
    assert!(seen.is_empty(), "rejects BAM contained unexpected records: {seen:?}");
}

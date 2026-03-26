//! Integration tests for the group command with new metrics infrastructure.

use fgoxide::io::DelimFile;
use fgumi_lib::metrics::UmiGroupingMetrics;
use noodles::bam;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record::data::field::Tag;
use std::collections::HashSet;
use std::fs;
use std::path::PathBuf;
use std::process::Command;
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, create_umi_family};

/// Collect MI tag values from all records in a BAM file.
fn collect_mi_tags(bam_path: &std::path::Path) -> Vec<String> {
    let mi_tag = Tag::new(b'M', b'I');
    let mut reader = bam::io::Reader::new(std::fs::File::open(bam_path).unwrap());
    let header = reader.read_header().unwrap();
    reader
        .record_bufs(&header)
        .filter_map(std::result::Result::ok)
        .filter_map(|r| {
            r.data().get(&mi_tag).and_then(|v| {
                if let noodles::sam::alignment::record_buf::data::field::Value::String(s) = v {
                    Some(String::from_utf8_lossy(s.as_ref()).into_owned())
                } else {
                    None
                }
            })
        })
        .collect()
}

/// Helper to write a BAM with UMI-tagged reads for group tests.
fn create_group_input_bam(
    path: &PathBuf,
    families: Vec<Vec<noodles::sam::alignment::record_buf::RecordBuf>>,
) {
    let header = create_minimal_header("chr1", 10000);
    let mut writer =
        bam::io::Writer::new(std::fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");
    for family in families {
        for record in family {
            writer.write_alignment_record(&header, &record).expect("Failed to write record");
        }
    }
    writer.try_finish().expect("Failed to finish BAM");
}

/// Test that the group command properly writes metrics in the new format.
#[test]
fn test_group_command_writes_new_metrics() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let metrics_file = temp_dir.path().join("metrics.txt");

    // Create test BAM file with multiple UMI families
    create_test_input_bam(&input_bam);

    // Run the group command
    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "group",
            "--input",
            input_bam.to_str().unwrap(),
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
            "--grouping-metrics",
            metrics_file.to_str().unwrap(),
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run group command");

    assert!(status.success(), "Group command failed");
    assert!(output_bam.exists(), "Output BAM not created");
    assert!(metrics_file.exists(), "Metrics file not created");

    // Read and validate metrics
    let metrics: Vec<UmiGroupingMetrics> =
        DelimFile::default().read_tsv(&metrics_file).expect("Failed to read metrics file");

    assert_eq!(metrics.len(), 1, "Expected exactly one metrics record");

    let metric = &metrics[0];

    // Verify basic fields are populated
    assert!(metric.total_records > 0, "total_records should be positive");
    assert!(metric.accepted_records > 0, "accepted_records should be positive");
    assert!(metric.unique_molecule_ids > 0, "unique_molecule_ids should be positive");

    // Rejection fields are validated by their presence in the struct (no runtime check needed
    // since they are unsigned integers that are always >= 0)
    let _ = metric.discarded_non_pf;
    let _ = metric.discarded_poor_alignment;
    let _ = metric.discarded_ns_in_umi;
    let _ = metric.discarded_umi_too_short;
}

/// Test that the group command handles UMIs with N bases correctly.
#[test]
fn test_group_command_rejects_n_bases_in_umi() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let metrics_file = temp_dir.path().join("metrics.txt");

    // Create BAM with UMIs containing N bases
    create_test_bam_with_n_umis(&input_bam);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "group",
            "--input",
            input_bam.to_str().unwrap(),
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
            "--grouping-metrics",
            metrics_file.to_str().unwrap(),
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run group command");

    assert!(status.success());

    // Read metrics and verify N rejection tracking
    let metrics: Vec<UmiGroupingMetrics> =
        DelimFile::default().read_tsv(&metrics_file).expect("Failed to read metrics");

    let metric = &metrics[0];
    assert!(metric.discarded_ns_in_umi > 0, "Should have discarded reads with N in UMI");
}

/// Helper function to create a test BAM file with multiple UMI families.
fn create_test_input_bam(path: &PathBuf) {
    let header = create_minimal_header("chr1", 10000);

    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));

    writer.write_header(&header).expect("Failed to write header");

    // Create 3 UMI families
    let family1 = create_umi_family("AAAAAAAA", 10, "family1", "ACGTACGT", 30);
    let family2 = create_umi_family("CCCCCCCC", 5, "family2", "TGCATGCA", 30);
    let family3 = create_umi_family("GGGGGGGG", 15, "family3", "ATCGATCG", 30);

    for record in family1.iter().chain(family2.iter()).chain(family3.iter()) {
        writer.write_alignment_record(&header, record).expect("Failed to write record");
    }

    writer.try_finish().expect("Failed to finish BAM");
}

/// Helper function to create a BAM file with UMIs containing N bases.
fn create_test_bam_with_n_umis(path: &PathBuf) {
    let header = create_minimal_header("chr1", 10000);

    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));

    writer.write_header(&header).expect("Failed to write header");

    // Create families with good and bad UMIs
    let good_family = create_umi_family("AAAAAAAA", 5, "good", "ACGTACGT", 30);
    let bad_family = create_umi_family("NNNNNNNN", 3, "bad", "TGCATGCA", 30);

    for record in good_family.iter().chain(bad_family.iter()) {
        writer.write_alignment_record(&header, record).expect("Failed to write record");
    }

    writer.try_finish().expect("Failed to finish BAM");
}

/// Without --truncate, reads with 5-bp and 6-bp UMIs should receive different MI tags
/// (they are in separate length-shape groups even at the same position).
#[test]
fn test_group_variable_length_umis_separated_without_truncate() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    let five_mer = create_umi_family("AAAAA", 3, "five", "AAAAGGGG", 30);
    let six_mer = create_umi_family("AAAAAA", 3, "six", "AAAAGGGG", 30);
    create_group_input_bam(&input_bam, vec![five_mer, six_mer]);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "group",
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
        .unwrap();
    assert!(status.success());

    let mi_tags = collect_mi_tags(&output_bam);
    let unique_mi: HashSet<_> = mi_tags.iter().collect();
    assert!(
        unique_mi.len() >= 2,
        "Expected >=2 MI groups for different-length UMIs, got: {unique_mi:?}"
    );
}

/// With --truncate --min-umi-length=5, reads with 5-bp and 6-bp UMIs that truncate to the
/// same sequence should be grouped together.
#[test]
fn test_group_truncate_merges_variable_length_umis() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    let five_mer = create_umi_family("AAAAA", 3, "five", "AAAAGGGG", 30);
    let six_mer = create_umi_family("AAAAAA", 3, "six", "AAAAGGGG", 30);
    create_group_input_bam(&input_bam, vec![five_mer, six_mer]);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "group",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--strategy",
            "identity",
            "--min-umi-length",
            "5",
            "--truncate",
            "--compression-level",
            "1",
        ])
        .status()
        .unwrap();
    assert!(status.success());

    let mi_tags = collect_mi_tags(&output_bam);
    let unique_mi: HashSet<_> = mi_tags.iter().collect();
    assert_eq!(
        unique_mi.len(),
        1,
        "Truncated UMIs should all share one MI group, got: {unique_mi:?}"
    );
}

/// --truncate without --min-umi-length should fail with a non-zero exit code.
#[test]
fn test_group_truncate_without_min_umi_length_fails() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    // Empty BAM is fine; the error is caught in validation before I/O.
    create_group_input_bam(&input_bam, vec![]);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "group",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--strategy",
            "identity",
            "--truncate",
        ])
        .status()
        .unwrap();
    assert!(!status.success(), "Should fail when --truncate given without --min-umi-length");
}

//! Integration tests for the group command with new metrics infrastructure.

use fgoxide::io::DelimFile;
use fgumi_lib::metrics::UmiGroupingMetrics;
use fgumi_lib::sam::builder::RecordBuilder;
use noodles::bam;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::sam::alignment::record_buf::data::field::value::Value as DataValue;
use std::collections::HashSet;
use std::fs;
use std::path::PathBuf;
use std::process::Command;
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, create_umi_family, to_record_buf};

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
        writer
            .write_alignment_record(&header, &to_record_buf(record))
            .expect("Failed to write record");
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
        writer
            .write_alignment_record(&header, &to_record_buf(record))
            .expect("Failed to write record");
    }

    writer.try_finish().expect("Failed to finish BAM");
}

/// Create a group of single-end reads at a specific position with the same UMI.
///
/// All reads share the same reference, position, and UMI, so they form a single
/// UMI family within their position group.
fn create_umi_family_at(base_name: &str, umi: &str, count: usize, start: usize) -> Vec<RecordBuf> {
    (0..count)
        .map(|i| {
            RecordBuilder::new()
                .name(&format!("{base_name}_{i}"))
                .sequence("ACGTACGT")
                .qualities(&[30; 8])
                .reference_sequence_id(0)
                .alignment_start(start)
                .mapping_quality(60)
                .cigar("8M")
                .tag("RX", umi)
                .build()
        })
        .collect()
}

/// Regression test for issue #269: the `group` command should emit consecutive
/// MI integer values in the range `[0, N-1]`, where `N` is the number of distinct
/// UMI families. fgbio's `GroupReadsByUmi` has this property; fgumi should too.
///
/// The input is constructed so that several position groups each contain more
/// templates than distinct UMI families, forcing the old per-thread block
/// allocator (which advanced by `templates.len()`) to leave gaps in the
/// emitted MI integer space.
#[test]
fn test_group_command_assigns_consecutive_mi_values() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    // Build the input BAM: three position groups, each with templates.len() > distinct UMIs.
    //
    // Group 1 (position 100):  4 reads, 1 UMI       -> 1 distinct MI, 3 wasted slots
    // Group 2 (position 500):  5 reads, 2 UMIs      -> 2 distinct MIs, 3 wasted slots
    // Group 3 (position 1000): 3 reads, 1 UMI       -> 1 distinct MI, 2 wasted slots
    // Total distinct MIs expected: 4
    let mut records: Vec<RecordBuf> = Vec::new();
    records.extend(create_umi_family_at("g1", "AAAAAAAA", 4, 100));
    records.extend(create_umi_family_at("g2a", "CCCCCCCC", 3, 500));
    records.extend(create_umi_family_at("g2b", "GGGGGGGG", 2, 500));
    records.extend(create_umi_family_at("g3", "TTTTTTTT", 3, 1000));

    let header = create_minimal_header("chr1", 10000);
    let mut writer =
        bam::io::Writer::new(fs::File::create(&input_bam).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");
    for record in &records {
        writer.write_alignment_record(&header, record).expect("Failed to write record");
    }
    writer.try_finish().expect("Failed to finish BAM");

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "group",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--strategy",
            "identity",
            "--edits",
            "0",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run group command");
    assert!(status.success(), "Group command failed");
    assert!(output_bam.exists(), "Output BAM not created");

    // Collect the set of MI tag values written to the output.
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let out_header = reader.read_header().expect("read output header");
    let mi_tag = Tag::from(*b"MI");

    let mut mi_values: HashSet<u64> = HashSet::new();
    for result in reader.record_bufs(&out_header) {
        let record = result.expect("read record");
        if let Some(DataValue::String(mi)) = record.data().get(&mi_tag) {
            let s = std::str::from_utf8(mi).expect("MI tag is utf8");
            let value: u64 = s.parse().expect("MI tag is a non-negative integer");
            mi_values.insert(value);
        }
    }

    assert_eq!(mi_values.len(), 4, "expected 4 distinct MI values (one per UMI family)");

    let max_mi = *mi_values.iter().max().expect("at least one MI value");
    let expected_max = (mi_values.len() as u64) - 1;
    assert_eq!(
        max_mi, expected_max,
        "MI values should be consecutive 0..N-1; got {mi_values:?} (max {max_mi}, expected {expected_max})"
    );
}

/// Regression test for issue #269 covering the multi-threaded pipeline path.
///
/// Same invariant (consecutive MI integers 0..N-1) but forces the unified BAM
/// pipeline by running with multiple threads. The counter in that path is an
/// `AtomicU64` shared across parallel serialize workers, so this test guards
/// against the per-worker block-allocation form of the bug as well.
#[test]
fn test_group_command_assigns_consecutive_mi_values_multi_threaded() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    // Build a larger input with many position groups so multiple worker threads
    // can reserve MI blocks in parallel. Each position group has more templates
    // than distinct UMI families so the old per-template block size would leak
    // gaps into the emitted MI space.
    let mut records: Vec<RecordBuf> = Vec::new();
    let num_groups = 20;
    for g in 0..num_groups {
        let start = 100 + (g * 500);
        // Two UMI families per position group, 4 and 3 reads respectively.
        records.extend(create_umi_family_at(&format!("g{g}_a"), "AAAAAAAA", 4, start));
        records.extend(create_umi_family_at(&format!("g{g}_b"), "CCCCCCCC", 3, start));
    }

    let header = create_minimal_header("chr1", 1_000_000);
    let mut writer =
        bam::io::Writer::new(fs::File::create(&input_bam).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");
    for record in &records {
        writer.write_alignment_record(&header, record).expect("Failed to write record");
    }
    writer.try_finish().expect("Failed to finish BAM");

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "group",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--strategy",
            "identity",
            "--edits",
            "0",
            "--threads",
            "4",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run group command");
    assert!(status.success(), "Group command failed");

    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let out_header = reader.read_header().expect("read output header");
    let mi_tag = Tag::from(*b"MI");

    let mut mi_values: HashSet<u64> = HashSet::new();
    for result in reader.record_bufs(&out_header) {
        let record = result.expect("read record");
        if let Some(DataValue::String(mi)) = record.data().get(&mi_tag) {
            let s = std::str::from_utf8(mi).expect("MI tag is utf8");
            let value: u64 = s.parse().expect("MI tag is a non-negative integer");
            mi_values.insert(value);
        }
    }

    let expected_distinct = num_groups * 2;
    assert_eq!(
        mi_values.len(),
        expected_distinct,
        "expected {expected_distinct} distinct MI values (two per position group)"
    );
    let max_mi = *mi_values.iter().max().expect("at least one MI value");
    let expected_max = (mi_values.len() as u64) - 1;
    assert_eq!(
        max_mi, expected_max,
        "MI values should be consecutive 0..N-1 under multi-threaded pipeline; got max {max_mi}, expected {expected_max}"
    );
}

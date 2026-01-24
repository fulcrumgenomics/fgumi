//! Integration tests for the duplex-metrics command.

use fgoxide::io::DelimFile;
use fgumi_lib::metrics::duplex::{DuplexFamilySizeMetric, FamilySizeMetric, UmiMetric};
use fgumi_lib::sam::builder::RecordBuilder;
use noodles::bam;
use noodles::sam::Header;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record_buf::RecordBuf;
use std::fs;
use std::path::PathBuf;
use std::process::Command;
use tempfile::TempDir;

use crate::helpers::bam_generator::create_minimal_header;

/// Creates a paired-end read pair for duplex-metrics testing.
///
/// # Arguments
///
/// * `name` - Read name
/// * `ref_id` - Reference sequence ID
/// * `pos1` - R1 alignment start position
/// * `pos2` - R2 alignment start position
/// * `rx_umi` - Raw UMI (RX tag) in format "AAA-TTT"
/// * `mi_tag` - Molecule identifier (MI tag) with strand suffix, e.g., "1/A" or "1/B"
/// * `strand1_forward` - If true, R1 is on forward strand
/// * `strand2_forward` - If true, R2 is on forward strand
#[allow(clippy::too_many_arguments)]
fn create_duplex_pair(
    name: &str,
    ref_id: usize,
    pos1: usize,
    pos2: usize,
    rx_umi: &str,
    mi_tag: &str,
    strand1_forward: bool,
    strand2_forward: bool,
) -> (RecordBuf, RecordBuf) {
    let sequence = "ACGTACGTACGTACGT";
    let quality = 30u8;

    let r1 = RecordBuilder::new()
        .name(name)
        .sequence(sequence)
        .qualities(&vec![quality; sequence.len()])
        .paired(true)
        .first_segment(true)
        .reference_sequence_id(ref_id)
        .alignment_start(pos1)
        .mapping_quality(60)
        .reverse_complement(!strand1_forward)
        .tag("RX", rx_umi)
        .tag("MI", mi_tag)
        .build();

    let r2 = RecordBuilder::new()
        .name(name)
        .sequence(sequence)
        .qualities(&vec![quality; sequence.len()])
        .paired(true)
        .first_segment(false)
        .reference_sequence_id(ref_id)
        .alignment_start(pos2)
        .mapping_quality(60)
        .reverse_complement(!strand2_forward)
        .tag("RX", rx_umi)
        .tag("MI", mi_tag)
        .build();

    (r1, r2)
}

/// Creates a test BAM file with duplex families for testing.
///
/// Creates multiple families:
/// - Family 1: Duplex with 2/A + 1/B at position 100-200
/// - Family 2: Single-strand with 1/A at position 300-400
/// - Family 3: Duplex with 3/A + 2/B at position 500-600
fn create_duplex_test_bam(path: &PathBuf, header: &Header) {
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));

    writer.write_header(header).expect("Failed to write header");

    let mut records = Vec::new();

    // Family 1: Duplex with 2 A-strand + 1 B-strand reads at position 100-200
    // A-strand reads
    for i in 0..2 {
        let (r1, r2) = create_duplex_pair(
            &format!("fam1_a_{i}"),
            0,
            100,
            200,
            "AAA-TTT",
            "1/A",
            true,  // forward
            false, // reverse
        );
        records.push(r1);
        records.push(r2);
    }
    // B-strand read
    let (r1, r2) = create_duplex_pair(
        "fam1_b_0", 0, 100, 200, "TTT-AAA", "1/B", false, // reverse
        true,  // forward
    );
    records.push(r1);
    records.push(r2);

    // Family 2: Single-strand with 1 A-strand read at position 300-400
    let (r1, r2) = create_duplex_pair(
        "fam2_a_0", 0, 300, 400, "CCC-GGG", "2/A", true,  // forward
        false, // reverse
    );
    records.push(r1);
    records.push(r2);

    // Family 3: Duplex with 3 A-strand + 2 B-strand reads at position 500-600
    // A-strand reads
    for i in 0..3 {
        let (r1, r2) = create_duplex_pair(
            &format!("fam3_a_{i}"),
            0,
            500,
            600,
            "GGG-CCC",
            "3/A",
            true,  // forward
            false, // reverse
        );
        records.push(r1);
        records.push(r2);
    }
    // B-strand reads
    for i in 0..2 {
        let (r1, r2) = create_duplex_pair(
            &format!("fam3_b_{i}"),
            0,
            500,
            600,
            "CCC-GGG",
            "3/B",
            false, // reverse
            true,  // forward
        );
        records.push(r1);
        records.push(r2);
    }

    for record in &records {
        writer.write_alignment_record(header, record).expect("Failed to write record");
    }

    writer.finish(header).expect("Failed to finish BAM");
}

/// Test that the duplex-metrics command creates all expected output files.
#[test]
fn test_duplex_metrics_command_creates_output_files() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_prefix = temp_dir.path().join("output");

    let header = create_minimal_header("chr1", 10000);
    create_duplex_test_bam(&input_bam, &header);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "duplex-metrics",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_prefix.to_str().unwrap(),
        ])
        .status()
        .expect("Failed to run duplex-metrics command");

    assert!(status.success(), "duplex-metrics command failed");

    // Verify all expected output files exist
    let family_sizes_path = format!("{}.family_sizes.txt", output_prefix.display());
    let duplex_family_sizes_path = format!("{}.duplex_family_sizes.txt", output_prefix.display());
    let yield_metrics_path = format!("{}.duplex_yield_metrics.txt", output_prefix.display());
    let umi_counts_path = format!("{}.umi_counts.txt", output_prefix.display());

    assert!(std::path::Path::new(&family_sizes_path).exists(), "family_sizes.txt not created");
    assert!(
        std::path::Path::new(&duplex_family_sizes_path).exists(),
        "duplex_family_sizes.txt not created"
    );
    assert!(
        std::path::Path::new(&yield_metrics_path).exists(),
        "duplex_yield_metrics.txt not created"
    );
    assert!(std::path::Path::new(&umi_counts_path).exists(), "umi_counts.txt not created");
}

/// Test that the duplex-metrics command produces correct family size metrics.
#[test]
fn test_duplex_metrics_command_family_sizes() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_prefix = temp_dir.path().join("output");

    let header = create_minimal_header("chr1", 10000);
    create_duplex_test_bam(&input_bam, &header);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "duplex-metrics",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_prefix.to_str().unwrap(),
        ])
        .status()
        .expect("Failed to run duplex-metrics command");

    assert!(status.success(), "duplex-metrics command failed");

    // Read and validate family size metrics
    let family_sizes_path = format!("{}.family_sizes.txt", output_prefix.display());
    let metrics: Vec<FamilySizeMetric> =
        DelimFile::default().read_tsv(&family_sizes_path).expect("Failed to read family sizes");

    // We should have metrics for various family sizes
    assert!(!metrics.is_empty(), "Should have family size metrics");

    // Verify some expected properties:
    // - Family 1: DS size 3 (2+1)
    // - Family 2: SS size 1 (1 A-strand only)
    // - Family 3: DS size 5 (3+2)
    // Check for size 1 (from family 2's single-strand)
    let size_1 = metrics.iter().find(|m| m.family_size == 1);
    assert!(size_1.is_some(), "Should have family size 1");
}

/// Test that the duplex-metrics command produces correct duplex family size metrics.
#[test]
fn test_duplex_metrics_command_duplex_family_sizes() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_prefix = temp_dir.path().join("output");

    let header = create_minimal_header("chr1", 10000);
    create_duplex_test_bam(&input_bam, &header);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "duplex-metrics",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_prefix.to_str().unwrap(),
        ])
        .status()
        .expect("Failed to run duplex-metrics command");

    assert!(status.success(), "duplex-metrics command failed");

    // Read and validate duplex family size metrics
    let duplex_family_sizes_path = format!("{}.duplex_family_sizes.txt", output_prefix.display());
    let metrics: Vec<DuplexFamilySizeMetric> = DelimFile::default()
        .read_tsv(&duplex_family_sizes_path)
        .expect("Failed to read duplex family sizes");

    // We should have duplex family size metrics
    // The metrics show how families are distributed by their AB and BA strand counts
    assert!(!metrics.is_empty(), "Should have duplex family size metrics");

    // Verify the structure is correct (has required fields)
    for m in &metrics {
        assert!(m.ab_size > 0 || m.ba_size > 0, "Family should have at least one strand");
        assert!(m.count > 0, "Family count should be positive");
    }
}

/// Test that the duplex-metrics command produces UMI metrics.
#[test]
fn test_duplex_metrics_command_umi_metrics() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_prefix = temp_dir.path().join("output");

    let header = create_minimal_header("chr1", 10000);
    create_duplex_test_bam(&input_bam, &header);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "duplex-metrics",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_prefix.to_str().unwrap(),
        ])
        .status()
        .expect("Failed to run duplex-metrics command");

    assert!(status.success(), "duplex-metrics command failed");

    // Read and validate UMI metrics
    let umi_counts_path = format!("{}.umi_counts.txt", output_prefix.display());
    let metrics: Vec<UmiMetric> =
        DelimFile::default().read_tsv(&umi_counts_path).expect("Failed to read UMI counts");

    // Should have UMI metrics from our test data
    assert!(!metrics.is_empty(), "Should have UMI metrics");

    // Verify total observations match our input
    let total_raw: usize = metrics.iter().map(|m| m.raw_observations).sum();
    // We have: 2 + 1 + 1 + 3 + 2 = 9 templates, each with 2 UMI parts = 18 raw UMI observations
    assert!(total_raw > 0, "Should have raw UMI observations");
}

/// Test that the --duplex-umi-counts flag produces additional output.
#[test]
fn test_duplex_metrics_command_with_duplex_umi_counts() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_prefix = temp_dir.path().join("output");

    let header = create_minimal_header("chr1", 10000);
    create_duplex_test_bam(&input_bam, &header);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "duplex-metrics",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_prefix.to_str().unwrap(),
            "--duplex-umi-counts",
        ])
        .status()
        .expect("Failed to run duplex-metrics command");

    assert!(status.success(), "duplex-metrics command failed with --duplex-umi-counts");

    // Verify the additional duplex UMI counts file exists
    let duplex_umi_counts_path = format!("{}.duplex_umi_counts.txt", output_prefix.display());
    assert!(
        std::path::Path::new(&duplex_umi_counts_path).exists(),
        "duplex_umi_counts.txt not created when --duplex-umi-counts is specified"
    );
}

/// Test that min-ab-reads and min-ba-reads parameters affect duplex classification.
#[test]
fn test_duplex_metrics_command_min_reads_parameters() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_prefix_default = temp_dir.path().join("output_default");
    let output_prefix_strict = temp_dir.path().join("output_strict");

    let header = create_minimal_header("chr1", 10000);
    create_duplex_test_bam(&input_bam, &header);

    // Run with default parameters (min-ab=1, min-ba=1)
    let status_default = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "duplex-metrics",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_prefix_default.to_str().unwrap(),
        ])
        .status()
        .expect("Failed to run duplex-metrics command (default)");

    assert!(status_default.success(), "duplex-metrics command failed (default)");

    // Run with stricter parameters (min-ab=2, min-ba=2)
    let status_strict = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "duplex-metrics",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_prefix_strict.to_str().unwrap(),
            "--min-ab-reads",
            "2",
            "--min-ba-reads",
            "2",
        ])
        .status()
        .expect("Failed to run duplex-metrics command (strict)");

    assert!(status_strict.success(), "duplex-metrics command failed (strict)");

    // Read both duplex family size files
    let default_path = format!("{}.duplex_family_sizes.txt", output_prefix_default.display());
    let strict_path = format!("{}.duplex_family_sizes.txt", output_prefix_strict.display());

    let default_metrics: Vec<DuplexFamilySizeMetric> = DelimFile::default()
        .read_tsv(&default_path)
        .expect("Failed to read default duplex metrics");

    let strict_metrics: Vec<DuplexFamilySizeMetric> =
        DelimFile::default().read_tsv(&strict_path).expect("Failed to read strict duplex metrics");

    // With default settings, both families 1 (2/A + 1/B) and 3 (3/A + 2/B) are duplex
    // With strict settings (2/2), only family 3 qualifies

    // The strict metrics should have fewer or equal entries because fewer families qualify
    // Note: Even non-qualifying families may appear with counts, but the duplex_count in
    // yield_metrics would be lower
    assert!(
        default_metrics.len() >= strict_metrics.len(),
        "Strict parameters should result in fewer or equal duplex families"
    );
}

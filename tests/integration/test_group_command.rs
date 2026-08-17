//! Integration tests for the group command with new metrics infrastructure.
//!
//! These tests invoke `GroupReadsByUmi::execute()` in-process.

use clap::Parser;
use fgoxide::io::DelimFile;
use fgumi_lib::commands::command::Command;
use fgumi_lib::commands::group::GroupReadsByUmi;
use fgumi_lib::metrics::UmiGroupingMetrics;
use noodles::bam;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use std::fs;
use std::path::PathBuf;
use tempfile::TempDir;

use crate::helpers::assertions::assert_bam_sorted;
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
    let cmd = GroupReadsByUmi::try_parse_from([
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
    .expect("failed to parse group args");
    cmd.execute("fgumi group").expect("Failed to run group command");
    assert!(output_bam.exists(), "Output BAM not created");
    assert!(metrics_file.exists(), "Metrics file not created");

    // group's `@HD` advertises SS:template-coordinate ("Output is always written in
    // template-coordinate order"); gate that claim against the bytes it wrote.
    assert_bam_sorted(&output_bam, "template-coordinate", Some("mi"));

    // Read and validate metrics
    let metrics: Vec<UmiGroupingMetrics> =
        DelimFile::default().read_tsv(&metrics_file).expect("Failed to read metrics file");

    assert_eq!(metrics.len(), 1, "Expected exactly one metrics record");

    let metric = &metrics[0];

    // Assert the exact serialized fgbio columns against `create_test_input_bam`'s known
    // structure: three families of 10 + 5 + 15 single-end reads, all PF, all mapq 60, all
    // carrying a valid 8bp UMI with no Ns — so every record is accepted and nothing is
    // discarded. The `.grouping_metrics.txt` matches fgbio's 5-column `UmiGroupingMetric`;
    // `total_records`/`unique_molecule_ids` are fgumi-internal (`#[serde(skip)]`) and read
    // back as zero, so they are not asserted here.
    assert_eq!(metric.accepted_records, 30, "all 30 input records should be accepted");
    assert_eq!(metric.discarded_non_pf, 0);
    assert_eq!(metric.discarded_poor_alignment, 0);
    assert_eq!(metric.discarded_ns_in_umi, 0);
    assert_eq!(metric.discarded_umi_too_short, 0);
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

    let cmd = GroupReadsByUmi::try_parse_from([
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
    .expect("failed to parse group args");
    cmd.execute("fgumi group").expect("Failed to run group command");

    // Read metrics and verify N rejection tracking
    let metrics: Vec<UmiGroupingMetrics> =
        DelimFile::default().read_tsv(&metrics_file).expect("Failed to read metrics");

    let metric = &metrics[0];
    assert!(metric.discarded_ns_in_umi > 0, "Should have discarded reads with N in UMI");
}

/// `group`'s metrics file is fgbio's `UmiGroupingMetric` and must stay exactly
/// five columns with fgbio's spellings, whatever fgumi tracks internally. The
/// struct-level test in `fgumi-metrics` pins the serialization; this pins what the
/// command actually writes to disk, so an internal refactor of the filter
/// accounting cannot change the file fgbio has to be able to read.
#[test]
fn test_grouping_metrics_header_is_fgbio_five_columns() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let metrics_file = temp_dir.path().join("metrics.txt");

    create_test_input_bam(&input_bam);

    let cmd = GroupReadsByUmi::try_parse_from([
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
    .expect("failed to parse group args");
    cmd.execute("fgumi group").expect("Failed to run group command");

    let content = fs::read_to_string(&metrics_file).expect("Failed to read metrics file");
    assert_eq!(
        content.lines().next(),
        Some(
            "accepted_sam_records\tdiscarded_non_pf\tdiscarded_poor_alignment\t\
             discarded_ns_in_umi\tdiscarded_umis_to_short"
        )
    );
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

// ============================================================================
// --check-crc / --no-check-crc on the single-threaded fast path (#800)
// ============================================================================

/// Flip a byte in the last BGZF block's CRC32 footer, so decoding that block
/// fails only when CRC verification is on. Requires the file to span at least
/// two blocks so the corrupted block is not the header's block (the header
/// parse always verifies).
fn corrupt_last_block_crc(path: &std::path::Path) {
    let mut bytes = fs::read(path).expect("read bam for corruption");
    let mut cursor: &[u8] = &bytes;
    let blocks = fgumi_lib::bgzf_reader::read_raw_blocks(&mut cursor, 100_000)
        .expect("read bgzf blocks from test bam");
    assert!(
        blocks.len() >= 2,
        "test input must span >= 2 BGZF blocks so the corrupted block isn't the header's; \
         got {} -- generate more records",
        blocks.len()
    );
    let offset: usize =
        blocks[..blocks.len() - 1].iter().map(fgumi_lib::bgzf_reader::RawBgzfBlock::len).sum();
    let last = blocks.last().expect("checked len >= 2 above");
    let crc_off = offset + last.len() - fgumi_lib::bgzf_reader::BGZF_FOOTER_SIZE;
    bytes[crc_off] ^= 0x01;
    fs::write(path, bytes).expect("write corrupted bam");
}

/// Write a single UMI family large enough to span more than one BGZF block (all
/// reads sit at one position on reference 0, so the input is template-coordinate
/// sorted), so the corrupted last block is a record block, not the header's.
fn create_multiblock_group_bam(path: &PathBuf) {
    let header = create_minimal_header("chr1", 10000);
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");
    for record in &create_umi_family("AAAAAAAA", 3000, "read", "ACGTACGT", 30) {
        writer
            .write_alignment_record(&header, &to_record_buf(record))
            .expect("Failed to write record");
    }
    writer.try_finish().expect("Failed to finish BAM");
}

/// Build group's argv for the single-threaded fast path (no `--threads`),
/// appending any extra flags (e.g. `--no-check-crc`).
fn group_args<'a>(input: &'a str, output: &'a str, extra: &[&'a str]) -> Vec<&'a str> {
    let mut args = vec![
        "group",
        "--input",
        input,
        "--output",
        output,
        "--strategy",
        "identity",
        "--edits",
        "0",
        "--compression-level",
        "1",
    ];
    args.extend_from_slice(extra);
    args
}

/// `--no-check-crc` must let group's single-threaded reader accept a corrupted
/// BGZF CRC32: it decodes through fgumi-bgzf, honoring the flag (#800). Against
/// the noodles reader this path used before, this could not pass.
#[test]
fn test_group_no_check_crc_accepts_corrupted_crc() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    create_multiblock_group_bam(&input_bam);
    corrupt_last_block_crc(&input_bam);

    let cmd = GroupReadsByUmi::try_parse_from(group_args(
        input_bam.to_str().unwrap(),
        output_bam.to_str().unwrap(),
        &["--no-check-crc"],
    ))
    .expect("failed to parse group args");
    cmd.execute("fgumi group")
        .expect("--no-check-crc must accept a corrupted BGZF CRC32 and complete");
    assert!(output_bam.exists(), "Output BAM not created");
}

/// Default (verify-on for file input) must reject the same corrupted CRC32.
#[test]
fn test_group_rejects_corrupted_crc_by_default() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    create_multiblock_group_bam(&input_bam);
    corrupt_last_block_crc(&input_bam);

    let cmd = GroupReadsByUmi::try_parse_from(group_args(
        input_bam.to_str().unwrap(),
        output_bam.to_str().unwrap(),
        &[],
    ))
    .expect("failed to parse group args");
    let err = cmd
        .execute("fgumi group")
        .expect_err("default (verify-on for file input) must reject a corrupted BGZF CRC32");
    let message = format!("{err:#}");
    assert!(message.to_uppercase().contains("CRC32"), "error should mention CRC32: {message}");
}

/// `--check-crc` must also reject the corrupted CRC32 (forces verification on).
#[test]
fn test_group_check_crc_rejects_corrupted_crc() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    create_multiblock_group_bam(&input_bam);
    corrupt_last_block_crc(&input_bam);

    let cmd = GroupReadsByUmi::try_parse_from(group_args(
        input_bam.to_str().unwrap(),
        output_bam.to_str().unwrap(),
        &["--check-crc"],
    ))
    .expect("failed to parse group args");
    let err =
        cmd.execute("fgumi group").expect_err("--check-crc must reject a corrupted BGZF CRC32");
    let message = format!("{err:#}");
    assert!(message.to_uppercase().contains("CRC32"), "error should mention CRC32: {message}");
}

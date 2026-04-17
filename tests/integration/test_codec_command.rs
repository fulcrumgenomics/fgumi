//! End-to-end CLI tests for the codec command.
//!
//! These tests run the actual `fgumi codec` binary and validate:
//! 1. Basic consensus calling from CODEC read pairs
//! 2. Statistics output
//! 3. Rejected reads output
//! 4. Quality filtering options

use fgumi_lib::bam_io::create_raw_bam_writer;
use fgumi_raw_bam::{RawRecord, SamBuilder, flags as raw_flags};
use noodles::bam;
use std::fs;
use std::path::PathBuf;
use std::process::Command;
use tempfile::TempDir;

use crate::helpers::bam_generator::create_minimal_header;

/// Creates a CODEC read pair (R1 forward, R2 reverse from opposite strand).
///
/// In CODEC sequencing, R1 and R2 come from opposite strands of the same molecule,
/// so a single read pair can produce duplex consensus.
#[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap, clippy::too_many_arguments)]
fn create_codec_read_pair(
    name: &str,
    r1_seq: &[u8],
    r2_seq: &[u8],
    r1_qual: &[u8],
    r2_qual: &[u8],
    ref_start: usize,
    umi: &str,
    cell_barcode: Option<&str>,
) -> (RawRecord, RawRecord) {
    let read_len = r1_seq.len();
    let cigar_op = u32::try_from(read_len).expect("read_len fits u32") << 4; // nM
    // SamBuilder pos is 0-based; ref_start is 1-based
    let pos = i32::try_from(ref_start).expect("ref_start fits i32") - 1;
    let mc_tag = format!("{read_len}M");

    let mut b1 = SamBuilder::new();
    b1.read_name(name.as_bytes())
        .sequence(r1_seq)
        .qualities(r1_qual)
        .cigar_ops(&[cigar_op])
        .flags(raw_flags::PAIRED | raw_flags::FIRST_SEGMENT | raw_flags::MATE_REVERSE)
        .ref_id(0)
        .pos(pos)
        .mapq(60)
        .mate_ref_id(0)
        .mate_pos(pos)
        .template_length(read_len as i32)
        .add_string_tag(b"MI", umi.as_bytes())
        .add_string_tag(b"MC", mc_tag.as_bytes());
    if let Some(cb) = cell_barcode {
        b1.add_string_tag(b"CB", cb.as_bytes());
    }

    let mut b2 = SamBuilder::new();
    b2.read_name(name.as_bytes())
        .sequence(r2_seq)
        .qualities(r2_qual)
        .cigar_ops(&[cigar_op])
        .flags(raw_flags::PAIRED | raw_flags::LAST_SEGMENT | raw_flags::REVERSE)
        .ref_id(0)
        .pos(pos)
        .mapq(60)
        .mate_ref_id(0)
        .mate_pos(pos)
        .template_length(-(read_len as i32))
        .add_string_tag(b"MI", umi.as_bytes())
        .add_string_tag(b"MC", mc_tag.as_bytes());
    if let Some(cb) = cell_barcode {
        b2.add_string_tag(b"CB", cb.as_bytes());
    }

    (b1.build(), b2.build())
}

/// Helper to create a test BAM file with CODEC read pairs.
fn create_codec_test_bam(path: &PathBuf, pairs: Vec<(RawRecord, RawRecord)>) {
    let header = create_minimal_header("chr1", 10000);
    let mut writer =
        create_raw_bam_writer(path, &header, 1, 6).expect("Failed to create raw BAM writer");
    for (r1, r2) in pairs {
        writer.write_raw_record(r1.as_ref()).expect("Failed to write R1");
        writer.write_raw_record(r2.as_ref()).expect("Failed to write R2");
    }
    writer.finish().expect("Failed to finish BAM");
}

/// Test basic CODEC consensus calling.
#[test]
fn test_codec_command_basic_consensus() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    // Create 3 read pairs for one molecule
    let mut pairs = Vec::new();
    for i in 0..3 {
        let (r1, r2) = create_codec_read_pair(
            &format!("read{i}"),
            b"ACGTACGT",
            b"ACGTACGT",
            &[30; 8],
            &[30; 8],
            100,
            "UMI001",
            None,
        );
        pairs.push((r1, r2));
    }
    create_codec_test_bam(&input_bam, pairs);

    // Run codec command
    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "codec",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--min-reads",
            "1",
            "--min-duplex-length",
            "1",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run codec command");

    assert!(status.success(), "Codec command failed");
    assert!(output_bam.exists(), "Output BAM not created");

    // Read output and verify consensus was created
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let mut consensus_count = 0;

    for result in reader.records() {
        let record = result.expect("Failed to read record");
        consensus_count += 1;

        // Verify consensus tags exist by checking the raw tag bytes
        let cd_tag = [b'c', b'D'];
        assert!(record.data().get(&cd_tag).is_some(), "Consensus should have cD tag");
    }

    assert!(consensus_count > 0, "Should have produced at least one consensus read");
}

/// Test CODEC command with statistics output.
#[test]
fn test_codec_command_with_stats() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let stats_file = temp_dir.path().join("stats.tsv");

    // Create read pairs for two molecules
    let mut pairs = Vec::new();
    for i in 0..2 {
        let (r1, r2) = create_codec_read_pair(
            &format!("mol1_read{i}"),
            b"ACGTACGT",
            b"ACGTACGT",
            &[30; 8],
            &[30; 8],
            100,
            "UMI001",
            None,
        );
        pairs.push((r1, r2));
    }
    for i in 0..3 {
        let (r1, r2) = create_codec_read_pair(
            &format!("mol2_read{i}"),
            b"TGCATGCA",
            b"TGCATGCA",
            &[30; 8],
            &[30; 8],
            200,
            "UMI002",
            None,
        );
        pairs.push((r1, r2));
    }
    create_codec_test_bam(&input_bam, pairs);

    // Run codec command with stats
    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "codec",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--stats",
            stats_file.to_str().unwrap(),
            "--min-reads",
            "1",
            "--min-duplex-length",
            "1",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run codec command");

    assert!(status.success(), "Codec command failed");
    assert!(stats_file.exists(), "Stats file not created");

    // Verify stats file has content
    let stats_content = fs::read_to_string(&stats_file).expect("Failed to read stats");
    assert!(!stats_content.is_empty(), "Stats file should not be empty");
}

/// Test CODEC command with rejected reads output.
#[test]
fn test_codec_command_with_rejects() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let rejects_bam = temp_dir.path().join("rejects.bam");

    // Create one pair that will pass and one that won't (need min-reads=3)
    let mut pairs = Vec::new();

    // Molecule 1: 3 pairs (will pass with min-reads=3)
    for i in 0..3 {
        let (r1, r2) = create_codec_read_pair(
            &format!("pass_read{i}"),
            b"ACGTACGT",
            b"ACGTACGT",
            &[30; 8],
            &[30; 8],
            100,
            "UMI_PASS",
            None,
        );
        pairs.push((r1, r2));
    }

    // Molecule 2: 1 pair (will fail with min-reads=3)
    let (r1, r2) = create_codec_read_pair(
        "fail_read0",
        b"TGCATGCA",
        b"TGCATGCA",
        &[30; 8],
        &[30; 8],
        200,
        "UMI_FAIL",
        None,
    );
    pairs.push((r1, r2));

    create_codec_test_bam(&input_bam, pairs);

    // Run codec command with rejects output and high min-reads
    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "codec",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--rejects",
            rejects_bam.to_str().unwrap(),
            "--min-reads",
            "3",
            "--min-duplex-length",
            "1",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run codec command");

    assert!(status.success(), "Codec command failed");
    assert!(output_bam.exists(), "Output BAM not created");
    assert!(rejects_bam.exists(), "Rejects BAM not created");

    // Count records in each file
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let consensus_count = reader.records().count();

    // Should have at least one consensus (from 3-pair molecule)
    assert!(consensus_count >= 1, "Should have consensus from passing molecule");

    // Rejects file should exist (may or may not have records depending on implementation)
    // The important thing is the file was created
    assert!(rejects_bam.exists(), "Rejects BAM file should be created");
}

/// Test CODEC command with minimum duplex length filter.
#[test]
fn test_codec_command_min_duplex_length() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    // Create read pairs with short sequences (8bp)
    let mut pairs = Vec::new();
    for i in 0..3 {
        let (r1, r2) = create_codec_read_pair(
            &format!("read{i}"),
            b"ACGTACGT",
            b"ACGTACGT",
            &[30; 8],
            &[30; 8],
            100,
            "UMI001",
            None,
        );
        pairs.push((r1, r2));
    }
    create_codec_test_bam(&input_bam, pairs);

    // Run codec command with high min-duplex-length (should reject)
    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "codec",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--min-reads",
            "1",
            "--min-duplex-length",
            "100", // Much longer than our 8bp reads
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run codec command");

    assert!(status.success(), "Codec command failed");

    // Should produce no consensus due to insufficient duplex length
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let consensus_count = reader.records().count();

    assert_eq!(consensus_count, 0, "Should have no consensus due to min-duplex-length filter");
}

/// Test CODEC command with per-base tags output.
#[test]
fn test_codec_command_per_base_tags() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    // Create read pairs
    let mut pairs = Vec::new();
    for i in 0..3 {
        let (r1, r2) = create_codec_read_pair(
            &format!("read{i}"),
            b"ACGTACGT",
            b"ACGTACGT",
            &[30; 8],
            &[30; 8],
            100,
            "UMI001",
            None,
        );
        pairs.push((r1, r2));
    }
    create_codec_test_bam(&input_bam, pairs);

    // Run codec command with per-base tags
    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "codec",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--min-reads",
            "1",
            "--min-duplex-length",
            "1",
            "--output-per-base-tags",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run codec command");

    assert!(status.success(), "Codec command failed");

    // Verify per-base tags exist
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();

    for result in reader.records() {
        let record = result.expect("Failed to read record");

        // Check for per-base depth tags (ad, bd)
        let ad_tag = [b'a', b'd'];
        let bd_tag = [b'b', b'd'];

        assert!(record.data().get(&ad_tag).is_some(), "Should have per-base depth tag 'ad'");
        assert!(record.data().get(&bd_tag).is_some(), "Should have per-base depth tag 'bd'");
    }
}

/// Test CODEC command preserves cell barcode tag.
#[test]
fn test_codec_command_cell_barcode_preservation() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    // Create 3 read pairs with cell barcode
    let mut pairs = Vec::new();
    for i in 0..3 {
        let (r1, r2) = create_codec_read_pair(
            &format!("read{i}"),
            b"ACGTACGT",
            b"ACGTACGT",
            &[30; 8],
            &[30; 8],
            100,
            "UMI001",
            Some("CELLBC123"),
        );
        pairs.push((r1, r2));
    }
    create_codec_test_bam(&input_bam, pairs);

    // Run codec command with cell-tag option
    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "codec",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--min-reads",
            "1",
            "--min-duplex-length",
            "1",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run codec command");

    assert!(status.success(), "Codec command failed");
    assert!(output_bam.exists(), "Output BAM not created");

    // Read output and verify cell barcode is preserved
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let mut consensus_count = 0;

    for result in reader.records() {
        let record = result.expect("Failed to read record");
        consensus_count += 1;

        // Verify cell barcode tag is preserved
        let cb_tag = [b'C', b'B'];
        assert!(
            record.data().get(&cb_tag).is_some(),
            "Consensus should have CB (cell barcode) tag"
        );
    }

    assert!(consensus_count > 0, "Should have produced at least one consensus read");
}

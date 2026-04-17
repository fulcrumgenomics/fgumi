//! Integration tests for the simplex-metrics command.

use fgoxide::io::DelimFile;
use fgumi_lib::metrics::shared::UmiMetric;
use fgumi_lib::metrics::simplex::{SimplexFamilySizeMetric, SimplexYieldMetric};
use fgumi_raw_bam::{RawRecord, SamBuilder, flags};
use noodles::bam;
use noodles::sam::Header;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use std::fs;
use std::path::PathBuf;
use std::process::Command;
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, to_record_buf};

/// Creates a paired-end read pair for simplex-metrics testing.
///
/// `pos1`/`pos2` are SAM-style 1-based positions; both must be `>= 1`.
fn create_simplex_pair(
    name: &str,
    ref_id: i32,
    pos1: i32,
    pos2: i32,
    rx_umi: &str,
    mi_tag: &str,
) -> (RawRecord, RawRecord) {
    const SEQUENCE: &[u8; 16] = b"ACGTACGTACGTACGT";
    const QUALS: [u8; 16] = [30; 16];

    assert!(pos1 >= 1 && pos2 >= 1, "pos1 and pos2 must be 1-based SAM positions (>= 1)");
    let cigar_op = u32::try_from(SEQUENCE.len()).expect("sequence.len() fits u32") << 4;
    let tlen = pos2 - pos1 + i32::try_from(SEQUENCE.len()).expect("sequence.len() fits i32");

    let r1 = {
        let mut b = SamBuilder::new();
        b.read_name(name.as_bytes())
            .sequence(SEQUENCE)
            .qualities(&QUALS)
            .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE)
            .ref_id(ref_id)
            .pos(pos1 - 1)
            .mapq(60)
            .cigar_ops(&[cigar_op])
            .mate_ref_id(ref_id)
            .mate_pos(pos2 - 1)
            .template_length(tlen)
            .add_string_tag(b"RX", rx_umi.as_bytes())
            .add_string_tag(b"MI", mi_tag.as_bytes());
        b.build()
    };

    let r2 = {
        let mut b = SamBuilder::new();
        b.read_name(name.as_bytes())
            .sequence(SEQUENCE)
            .qualities(&QUALS)
            .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE)
            .ref_id(ref_id)
            .pos(pos2 - 1)
            .mapq(60)
            .cigar_ops(&[cigar_op])
            .mate_ref_id(ref_id)
            .mate_pos(pos1 - 1)
            .template_length(-tlen)
            .add_string_tag(b"RX", rx_umi.as_bytes())
            .add_string_tag(b"MI", mi_tag.as_bytes());
        b.build()
    };

    (r1, r2)
}

/// Creates a test BAM file with simplex families for testing.
///
/// Creates multiple families:
/// - Family 1: 3 reads with MI 1/A at position 100-200 (SS family size 3)
/// - Family 2: 1 read with MI 2/A at position 100-200 (SS family size 1, same CS as family 1)
/// - Family 3: 2 reads with MI 3/A at position 300-400 (SS family size 2)
fn create_simplex_test_bam(path: &PathBuf, header: &Header) {
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(header).expect("Failed to write header");

    let mut records = Vec::new();

    // Family 1: 3 reads at pos 100-200 with MI 1/A
    for i in 0..3 {
        let (r1, r2) = create_simplex_pair(&format!("fam1_{i}"), 0, 100, 200, "AAA", "1/A");
        records.push(r1);
        records.push(r2);
    }

    // Family 2: 1 read at pos 100-200 with MI 2/A (same coordinate group as family 1)
    let (r1, r2) = create_simplex_pair("fam2_0", 0, 100, 200, "CCC", "2/A");
    records.push(r1);
    records.push(r2);

    // Family 3: 2 reads at pos 300-400 with MI 3/A
    for i in 0..2 {
        let (r1, r2) = create_simplex_pair(&format!("fam3_{i}"), 0, 300, 400, "GGG", "3/A");
        records.push(r1);
        records.push(r2);
    }

    for record in &records {
        writer
            .write_alignment_record(header, &to_record_buf(record))
            .expect("Failed to write record");
    }
    writer.finish(header).expect("Failed to finish BAM");
}

#[test]
fn test_simplex_metrics_command_creates_output_files() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_prefix = temp_dir.path().join("output");

    let header = create_minimal_header("chr1", 10000);
    create_simplex_test_bam(&input_bam, &header);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "simplex-metrics",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_prefix.to_str().unwrap(),
        ])
        .status()
        .expect("Failed to run simplex-metrics command");

    assert!(status.success(), "simplex-metrics command failed");

    let family_sizes_path = format!("{}.family_sizes.txt", output_prefix.display());
    let yield_metrics_path = format!("{}.simplex_yield_metrics.txt", output_prefix.display());
    let umi_counts_path = format!("{}.umi_counts.txt", output_prefix.display());

    assert!(std::path::Path::new(&family_sizes_path).exists(), "family_sizes.txt not created");
    assert!(
        std::path::Path::new(&yield_metrics_path).exists(),
        "simplex_yield_metrics.txt not created"
    );
    assert!(std::path::Path::new(&umi_counts_path).exists(), "umi_counts.txt not created");
}

#[test]
fn test_simplex_metrics_command_family_sizes() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_prefix = temp_dir.path().join("output");

    let header = create_minimal_header("chr1", 10000);
    create_simplex_test_bam(&input_bam, &header);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "simplex-metrics",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_prefix.to_str().unwrap(),
        ])
        .status()
        .expect("Failed to run simplex-metrics command");

    assert!(status.success());

    let family_sizes_path = format!("{}.family_sizes.txt", output_prefix.display());
    let metrics: Vec<SimplexFamilySizeMetric> =
        DelimFile::default().read_tsv(&family_sizes_path).expect("Failed to read family sizes");

    assert!(!metrics.is_empty(), "Should have family size metrics");

    // Family 2 has SS size 1; Family 3 has SS size 2; Family 1 has SS size 3
    let size_1 = metrics.iter().find(|m| m.family_size == 1).unwrap();
    assert_eq!(size_1.ss_count, 1, "One SS family of size 1 (family 2)");

    let size_2 = metrics.iter().find(|m| m.family_size == 2).unwrap();
    assert_eq!(size_2.ss_count, 1, "One SS family of size 2 (family 3)");

    let size_3 = metrics.iter().find(|m| m.family_size == 3).unwrap();
    assert_eq!(size_3.ss_count, 1, "One SS family of size 3 (family 1)");
}

#[test]
fn test_simplex_metrics_command_yield_metrics() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_prefix = temp_dir.path().join("output");

    let header = create_minimal_header("chr1", 10000);
    create_simplex_test_bam(&input_bam, &header);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "simplex-metrics",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_prefix.to_str().unwrap(),
        ])
        .status()
        .expect("Failed to run simplex-metrics command");

    assert!(status.success());

    let yield_path = format!("{}.simplex_yield_metrics.txt", output_prefix.display());
    let metrics: Vec<SimplexYieldMetric> =
        DelimFile::default().read_tsv(&yield_path).expect("Failed to read yield metrics");

    assert_eq!(metrics.len(), 20, "Should have 20 downsampling fractions");

    // 100% fraction should have all families
    let full = metrics.last().expect("Should have at least one yield metric");
    assert!((full.fraction - 1.0).abs() < 0.01);
    assert_eq!(full.ss_families, 3, "Should have 3 SS families at 100%");
    assert_eq!(full.cs_families, 2, "Should have 2 CS families at 100%");
}

#[test]
fn test_simplex_metrics_command_umi_metrics() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_prefix = temp_dir.path().join("output");

    let header = create_minimal_header("chr1", 10000);
    create_simplex_test_bam(&input_bam, &header);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "simplex-metrics",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_prefix.to_str().unwrap(),
        ])
        .status()
        .expect("Failed to run simplex-metrics command");

    assert!(status.success());

    let umi_counts_path = format!("{}.umi_counts.txt", output_prefix.display());
    let metrics: Vec<UmiMetric> =
        DelimFile::default().read_tsv(&umi_counts_path).expect("Failed to read UMI counts");

    assert!(!metrics.is_empty(), "Should have UMI metrics");
    let total_raw: usize = metrics.iter().map(|m| m.raw_observations).sum();
    assert!(total_raw > 0, "Should have raw UMI observations");
}

#[test]
fn test_simplex_metrics_min_reads_parameter() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_default = temp_dir.path().join("output_default");
    let output_strict = temp_dir.path().join("output_strict");

    let header = create_minimal_header("chr1", 10000);
    create_simplex_test_bam(&input_bam, &header);

    // Run with default min-reads=1
    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "simplex-metrics",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_default.to_str().unwrap(),
        ])
        .status()
        .expect("Failed to run simplex-metrics (default)");
    assert!(status.success());

    // Run with min-reads=3
    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "simplex-metrics",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_strict.to_str().unwrap(),
            "--min-reads",
            "3",
        ])
        .status()
        .expect("Failed to run simplex-metrics (strict)");
    assert!(status.success());

    let default_yield_path = format!("{}.simplex_yield_metrics.txt", output_default.display());
    let strict_yield_path = format!("{}.simplex_yield_metrics.txt", output_strict.display());

    let default_metrics: Vec<SimplexYieldMetric> =
        DelimFile::default().read_tsv(&default_yield_path).unwrap();
    let strict_metrics: Vec<SimplexYieldMetric> =
        DelimFile::default().read_tsv(&strict_yield_path).unwrap();

    // At 100% fraction
    let default_full = default_metrics.last().expect("Should have at least one yield metric");
    let strict_full = strict_metrics.last().expect("Should have at least one yield metric");

    // With min-reads=1, all 3 SS families are consensus-eligible
    // With min-reads=3, only family 1 (size 3) is consensus-eligible
    assert_eq!(default_full.ss_consensus_families, 3);
    assert_eq!(strict_full.ss_consensus_families, 1);
}

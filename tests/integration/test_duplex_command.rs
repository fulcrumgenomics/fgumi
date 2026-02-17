//! End-to-end CLI tests for the duplex command.
//!
//! These tests run the actual `fgumi duplex` binary and validate:
//! 1. Basic duplex consensus calling from paired-UMI grouped reads
//! 2. Statistics output
//! 3. Rejected reads output

use fgumi_lib::sam::builder::RecordBuilder;
use noodles::bam;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::sam::alignment::record_buf::data::field::Value;
use std::fs;
use std::path::Path;
use std::process::Command;
use tempfile::TempDir;

use crate::helpers::bam_generator::create_minimal_header;

/// Create a paired-end read pair for duplex consensus testing.
///
/// Duplex consensus requires reads grouped by MI tag with /A and /B strand suffixes.
/// For /A strand: R1 forward at pos, R2 reverse at pos+100 (FR orientation)
/// For /B strand: R1 reverse at pos+100, R2 forward at pos (RF orientation)
fn create_duplex_read_pair(
    name: &str,
    mi_tag: &str,
    sequence: &str,
    quality: u8,
    ref_start: usize,
    is_b_strand: bool,
) -> (RecordBuf, RecordBuf) {
    let read_len = sequence.len();
    let cigar = format!("{read_len}M");

    let (r1_start, r2_start, r1_rev, r2_rev) = if is_b_strand {
        // B strand: RF orientation — R1 reverse at far position, R2 forward at near
        (ref_start + 100, ref_start, true, false)
    } else {
        // A strand: FR orientation — R1 forward at near position, R2 reverse at far
        (ref_start, ref_start + 100, false, true)
    };

    #[expect(
        clippy::cast_possible_truncation,
        clippy::cast_possible_wrap,
        reason = "test data with known small values"
    )]
    let tlen: i32 = if is_b_strand { -((read_len + 100) as i32) } else { (read_len + 100) as i32 };

    let mut r1 = RecordBuilder::new()
        .name(name)
        .sequence(sequence)
        .qualities(&vec![quality; read_len])
        .paired(true)
        .first_segment(true)
        .reverse_complement(r1_rev)
        .mate_reverse_complement(r2_rev)
        .reference_sequence_id(0)
        .alignment_start(r1_start)
        .mapping_quality(60)
        .cigar(&cigar)
        .mate_reference_sequence_id(0)
        .mate_alignment_start(r2_start)
        .template_length(tlen)
        .build();

    let mut r2 = RecordBuilder::new()
        .name(name)
        .sequence(sequence)
        .qualities(&vec![quality; read_len])
        .paired(true)
        .first_segment(false)
        .reverse_complement(r2_rev)
        .mate_reverse_complement(r1_rev)
        .reference_sequence_id(0)
        .alignment_start(r2_start)
        .mapping_quality(60)
        .cigar(&cigar)
        .mate_reference_sequence_id(0)
        .mate_alignment_start(r1_start)
        .template_length(-tlen)
        .build();

    // Add MI tag
    let mi = Tag::new(b'M', b'I');
    r1.data_mut().insert(mi, Value::from(mi_tag));
    r2.data_mut().insert(mi, Value::from(mi_tag));

    (r1, r2)
}

/// Create a BAM with duplex-grouped reads (MI tags with /A and /B strand suffixes).
fn create_duplex_bam(path: &Path, molecules: Vec<Vec<(RecordBuf, RecordBuf)>>) {
    let header = create_minimal_header("chr1", 10000);
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");

    for pairs in molecules {
        for (r1, r2) in pairs {
            writer.write_alignment_record(&header, &r1).expect("Failed to write R1");
            writer.write_alignment_record(&header, &r2).expect("Failed to write R2");
        }
    }
    writer.finish(&header).expect("Failed to finish BAM");
}

/// Create a single duplex molecule with `depth` read pairs on each strand.
fn create_duplex_molecule(
    mi_id: &str,
    sequence: &str,
    quality: u8,
    ref_start: usize,
    depth: usize,
) -> Vec<(RecordBuf, RecordBuf)> {
    let mut molecule = Vec::new();
    for i in 0..depth {
        let (r1, r2) = create_duplex_read_pair(
            &format!("ab_{i}"),
            &format!("{mi_id}/A"),
            sequence,
            quality,
            ref_start,
            false,
        );
        molecule.push((r1, r2));
    }
    for i in 0..depth {
        let (r1, r2) = create_duplex_read_pair(
            &format!("ba_{i}"),
            &format!("{mi_id}/B"),
            sequence,
            quality,
            ref_start,
            true,
        );
        molecule.push((r1, r2));
    }
    molecule
}

/// Test basic duplex consensus calling.
#[test]
fn test_duplex_command_basic() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    let molecule = create_duplex_molecule("1", "ACGTACGT", 30, 100, 3);
    create_duplex_bam(&input_bam, vec![molecule]);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "duplex",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--min-reads",
            "1",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run duplex command");

    assert!(status.success(), "Duplex command failed");
    assert!(output_bam.exists(), "Output BAM not created");

    // Verify consensus reads were produced
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let mut count = 0;
    for result in reader.records() {
        let record = result.expect("Failed to read record");
        let cd_tag = [b'c', b'D'];
        assert!(record.data().get(&cd_tag).is_some(), "Duplex consensus should have cD tag");
        count += 1;
    }
    assert!(count > 0, "Should have produced duplex consensus reads");
}

/// Test duplex command with statistics output.
#[test]
fn test_duplex_command_with_stats() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let stats_path = temp_dir.path().join("stats.txt");

    let molecule = create_duplex_molecule("1", "ACGTACGT", 30, 100, 3);
    create_duplex_bam(&input_bam, vec![molecule]);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "duplex",
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
        .expect("Failed to run duplex command");

    assert!(status.success(), "Duplex command with stats failed");
    assert!(stats_path.exists(), "Stats file not created");
}

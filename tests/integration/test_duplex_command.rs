//! End-to-end CLI tests for the duplex command.
//!
//! These tests run the actual `fgumi duplex` binary and validate:
//! 1. Basic duplex consensus calling from paired-UMI grouped reads
//! 2. Statistics output
//! 3. Rejected reads output

use fgumi_lib::sam::SamTag;
use fgumi_raw_bam::{RawRecord, SamBuilder, flags};
use noodles::bam;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use std::fs;
use std::path::Path;
use std::process::Command;
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, to_record_buf};

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
    ref_start: i32,
    is_b_strand: bool,
) -> (RawRecord, RawRecord) {
    let seq = sequence.as_bytes();
    let read_len = seq.len();
    let cigar_op = u32::try_from(read_len).expect("read_len fits u32") << 4;

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

    let r1_flags = flags::PAIRED
        | flags::FIRST_SEGMENT
        | if r1_rev { flags::REVERSE } else { 0 }
        | if r2_rev { flags::MATE_REVERSE } else { 0 };

    let r2_flags = flags::PAIRED
        | flags::LAST_SEGMENT
        | if r2_rev { flags::REVERSE } else { 0 }
        | if r1_rev { flags::MATE_REVERSE } else { 0 };

    let r1 = {
        let mut b = SamBuilder::new();
        b.read_name(name.as_bytes())
            .sequence(seq)
            .qualities(&vec![quality; read_len])
            .flags(r1_flags)
            .ref_id(0)
            .pos(r1_start - 1)
            .mapq(60)
            .cigar_ops(&[cigar_op])
            .mate_ref_id(0)
            .mate_pos(r2_start - 1)
            .template_length(tlen)
            .add_string_tag(SamTag::MI, mi_tag.as_bytes());
        b.build()
    };

    let r2 = {
        let mut b = SamBuilder::new();
        b.read_name(name.as_bytes())
            .sequence(seq)
            .qualities(&vec![quality; read_len])
            .flags(r2_flags)
            .ref_id(0)
            .pos(r2_start - 1)
            .mapq(60)
            .cigar_ops(&[cigar_op])
            .mate_ref_id(0)
            .mate_pos(r1_start - 1)
            .template_length(-tlen)
            .add_string_tag(SamTag::MI, mi_tag.as_bytes());
        b.build()
    };

    (r1, r2)
}

/// Create a BAM with duplex-grouped reads (MI tags with /A and /B strand suffixes).
fn create_duplex_bam(path: &Path, molecules: Vec<Vec<(RawRecord, RawRecord)>>) {
    let header = create_minimal_header("chr1", 10000);
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");

    for pairs in molecules {
        for (r1, r2) in pairs {
            writer
                .write_alignment_record(&header, &to_record_buf(&r1))
                .expect("Failed to write R1");
            writer
                .write_alignment_record(&header, &to_record_buf(&r2))
                .expect("Failed to write R2");
        }
    }
    writer.try_finish().expect("Failed to finish BAM");
}

/// Create a single duplex molecule with `depth` read pairs on each strand.
fn create_duplex_molecule(
    mi_id: &str,
    sequence: &str,
    quality: u8,
    ref_start: i32,
    depth: usize,
) -> Vec<(RawRecord, RawRecord)> {
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
        let cd_tag = SamTag::CD.to_noodles_tag();
        assert!(record.data().get(&cd_tag).is_some(), "Duplex consensus should have cD tag");
        count += 1;
    }
    assert!(count > 0, "Should have produced duplex consensus reads");
}

/// Test duplex command with rejects output.
///
/// Runs the multi-threaded pipeline with `--min-reads 2` and a singleton /A
/// template that fails the single-strand `min-reads` check. Verifies that:
///
/// 1. The rejects BAM is created and has its `@HD` header rewritten to
///    `SO:unsorted` because reject writes are emitted in mutex-acquisition
///    order rather than input order.
/// 2. The singleton template's raw records are actually streamed to the
///    rejects BAM rather than being silently dropped.
#[test]
fn test_duplex_command_with_rejects() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let rejects_bam = temp_dir.path().join("rejects.bam");

    // Kept molecule: 3 /A pairs + 3 /B pairs
    let kept = create_duplex_molecule("1", "ACGTACGT", 30, 100, 3);
    // Singleton /A pair (fails the /A strand's min-reads=2 check)
    let singleton = vec![create_duplex_read_pair("solo", "2/A", "ACGTACGT", 30, 500, false)];
    create_duplex_bam(&input_bam, vec![kept, singleton]);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "duplex",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--rejects",
            rejects_bam.to_str().unwrap(),
            "--min-reads",
            "2",
            "--threads",
            "2",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run duplex command");

    assert!(status.success(), "Duplex command with rejects failed");
    assert!(rejects_bam.exists(), "Rejects BAM not created");

    crate::helpers::assertions::assert_header_unsorted(&rejects_bam);

    // The singleton /A template (R1 + R2) should be streamed to the rejects BAM.
    let mut reader = bam::io::Reader::new(fs::File::open(&rejects_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let mut reject_count = 0;
    for result in reader.records() {
        result.expect("Failed to read reject record");
        reject_count += 1;
    }
    assert_eq!(
        reject_count, 2,
        "Singleton paired-end /A template should be streamed to rejects BAM (R1 + R2)"
    );
}

/// Regression test: each rejected input record must appear in the rejects BAM
/// exactly once.
///
/// The duplex caller accumulates rejects from two sources: records dropped at
/// the single-strand (ss) layer, and records dropped at the duplex layer. When
/// the duplex layer rejects a group it returns every raw input record — which
/// already includes anything the ss layer would have rejected. Appending both
/// sources naively would emit the overlapping records twice.
///
/// This test drives several independent rejection paths in a single pipeline
/// run and asserts that the rejects BAM contains each `(read_name, flags)`
/// tuple at most once, with the total record count matching the expected
/// input-record count.
#[test]
fn test_duplex_command_rejects_contain_no_duplicates() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let rejects_bam = temp_dir.path().join("rejects.bam");

    // Kept molecule: produces a consensus, contributes zero rejects.
    let kept = create_duplex_molecule("1", "ACGTACGT", 30, 100, 3);
    // Four singleton molecules (2 records each) — two on each strand — that
    // reject via the insufficient-reads path. Using molecule IDs in a single
    // `singletons` vec avoids the similar_names clippy lint.
    let singletons: Vec<Vec<(RawRecord, RawRecord)>> = vec![
        vec![create_duplex_read_pair("soloA1", "2/A", "ACGTACGT", 30, 500, false)],
        vec![create_duplex_read_pair("soloA2", "3/A", "ACGTACGT", 30, 800, false)],
        vec![create_duplex_read_pair("soloB1", "4/B", "ACGTACGT", 30, 1100, true)],
        vec![create_duplex_read_pair("soloB2", "5/B", "ACGTACGT", 30, 1400, true)],
    ];

    // Expected rejects = 4 singletons × 2 records = 8. The kept molecule
    // contributes 0 rejects.
    let expected_rejects = u32::try_from(singletons.len() * 2).expect("rejects count fits in u32");

    let mut molecules = vec![kept];
    molecules.extend(singletons);
    create_duplex_bam(&input_bam, molecules);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "duplex",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--rejects",
            rejects_bam.to_str().unwrap(),
            "--min-reads",
            "2",
            "--threads",
            "2",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run duplex command");

    assert!(status.success(), "Duplex command with rejects failed");
    assert!(rejects_bam.exists(), "Rejects BAM not created");

    crate::helpers::assertions::assert_has_bgzf_eof(&rejects_bam);
    crate::helpers::assertions::assert_header_unsorted(&rejects_bam);

    let mut reader = bam::io::Reader::new(fs::File::open(&rejects_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let mut seen: std::collections::HashMap<(String, u16), u32> = std::collections::HashMap::new();
    let mut total = 0u32;
    for result in reader.records() {
        let record = result.expect("Failed to read reject record");
        let name = record.name().expect("reject record missing read name").to_string();
        let flags = u16::from(record.flags());
        *seen.entry((name, flags)).or_insert(0) += 1;
        total += 1;
    }

    assert_eq!(
        total, expected_rejects,
        "rejects BAM record count should match the expected rejected-input count",
    );
    for ((name, flags), count) in &seen {
        assert_eq!(
            *count, 1,
            "record ({name}, flags={flags}) appears {count} times in the rejects BAM",
        );
    }
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

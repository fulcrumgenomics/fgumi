//! Integration tests for the fastq command.

use fgumi_raw_bam::{SamBuilder, flags};
use noodles::bam;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use std::fs;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::process::Command;
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, to_record_buf};

/// Create a BAM file with paired-end reads for testing.
fn create_paired_bam(path: &PathBuf, read_pairs: Vec<(&str, &str, &str, &str, &str, bool)>) {
    let header = create_minimal_header("chr1", 10000);

    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));

    writer.write_header(&header).expect("Failed to write header");

    // Create paired-end records
    for (name, seq1, qual1, seq2, qual2, r2_reverse) in read_pairs {
        let q1: Vec<u8> = qual1.bytes().map(|b| b - 33).collect();
        let q2: Vec<u8> = qual2.bytes().map(|b| b - 33).collect();

        // R1 (forward strand)
        let r1 = {
            let mut b = SamBuilder::new();
            b.read_name(name.as_bytes())
                .sequence(seq1.as_bytes())
                .qualities(&q1)
                .flags(flags::PAIRED | flags::FIRST_SEGMENT)
                .ref_id(0)
                .pos(99)
                .mapq(60);
            b.build()
        };

        writer.write_alignment_record(&header, &to_record_buf(&r1)).expect("Failed to write R1");

        // R2 (optionally reverse complemented)
        let r2_flags =
            flags::PAIRED | flags::LAST_SEGMENT | if r2_reverse { flags::REVERSE } else { 0 };
        let r2 = {
            let mut b = SamBuilder::new();
            b.read_name(name.as_bytes())
                .sequence(seq2.as_bytes())
                .qualities(&q2)
                .flags(r2_flags)
                .ref_id(0)
                .pos(199)
                .mapq(60);
            b.build()
        };

        writer.write_alignment_record(&header, &to_record_buf(&r2)).expect("Failed to write R2");
    }

    writer.try_finish().expect("Failed to finish BAM");
}

/// Parse FASTQ records from a file.
fn parse_fastq_records(path: &PathBuf) -> Vec<(String, String, String)> {
    let file = fs::File::open(path).expect("Failed to open FASTQ");
    let reader = BufReader::new(file);
    let lines: Vec<String> = reader.lines().map(|l| l.unwrap()).collect();

    let mut records = Vec::new();
    for chunk in lines.chunks(4) {
        if chunk.len() == 4 {
            let name = chunk[0].trim_start_matches('@').to_string();
            let seq = chunk[1].clone();
            let qual = chunk[3].clone();
            records.push((name, seq, qual));
        }
    }
    records
}

/// Test basic fastq conversion.
#[test]
fn test_fastq_basic() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_fq = temp_dir.path().join("output.fq");

    // Create input BAM with 2 read pairs
    create_paired_bam(
        &input_bam,
        vec![
            ("read1", "ACGTACGT", "IIIIIIII", "TGCATGCA", "IIIIIIII", false),
            ("read2", "AAAACCCC", "IIIIIIII", "GGGGTTTT", "IIIIIIII", false),
        ],
    );

    // Run fastq command with output redirected to file
    let output = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args(["fastq", "-i", input_bam.to_str().unwrap()])
        .output()
        .expect("Failed to run fastq command");

    assert!(
        output.status.success(),
        "Fastq command failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Write stdout to file for parsing
    fs::write(&output_fq, &output.stdout).expect("Failed to write output");

    // Verify output
    let records = parse_fastq_records(&output_fq);
    assert_eq!(records.len(), 4, "Should have 4 FASTQ records (2 pairs)");

    // Check read names have correct suffixes
    assert_eq!(records[0].0, "read1/1");
    assert_eq!(records[1].0, "read1/2");
    assert_eq!(records[2].0, "read2/1");
    assert_eq!(records[3].0, "read2/2");

    // Check sequences
    assert_eq!(records[0].1, "ACGTACGT");
    assert_eq!(records[1].1, "TGCATGCA");
    assert_eq!(records[2].1, "AAAACCCC");
    assert_eq!(records[3].1, "GGGGTTTT");
}

/// Test fastq with reverse complemented reads.
#[test]
fn test_fastq_reverse_complement() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_fq = temp_dir.path().join("output.fq");

    // Create input BAM with R2 reverse complemented
    // The stored sequence is "AAAA", but since it's reverse complemented,
    // the output should be reverse complement: "TTTT"
    create_paired_bam(
        &input_bam,
        vec![
            // R2 is stored as "AAAA" but marked as reverse complemented
            // Output should be reverse complement: "TTTT"
            ("read1", "ACGTACGT", "IIIIIIII", "AAAA", "IIII", true),
        ],
    );

    let output = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args(["fastq", "-i", input_bam.to_str().unwrap()])
        .output()
        .expect("Failed to run fastq command");

    assert!(
        output.status.success(),
        "Fastq command failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    fs::write(&output_fq, &output.stdout).expect("Failed to write output");

    let records = parse_fastq_records(&output_fq);
    assert_eq!(records.len(), 2);

    // R1 should be unchanged
    assert_eq!(records[0].1, "ACGTACGT");

    // R2 should be reverse complemented: AAAA -> TTTT
    assert_eq!(records[1].1, "TTTT", "R2 should be reverse complemented from AAAA to TTTT");
}

/// Test fastq with no-suffix option.
#[test]
fn test_fastq_no_suffix() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_fq = temp_dir.path().join("output.fq");

    create_paired_bam(&input_bam, vec![("read1", "ACGT", "IIII", "TGCA", "IIII", false)]);

    let output = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "fastq",
            "-i",
            input_bam.to_str().unwrap(),
            "-n", // no suffix
        ])
        .output()
        .expect("Failed to run fastq command");

    assert!(output.status.success());

    fs::write(&output_fq, &output.stdout).expect("Failed to write output");

    let records = parse_fastq_records(&output_fq);
    assert_eq!(records.len(), 2);

    // Read names should NOT have /1 and /2 suffixes
    assert_eq!(records[0].0, "read1");
    assert_eq!(records[1].0, "read1");
}

/// Test quality score encoding (Phred+33).
#[test]
fn test_fastq_quality_encoding() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_fq = temp_dir.path().join("output.fq");

    // Create BAM with specific quality scores
    // ASCII 'I' (73) - 33 = quality 40
    // ASCII '!' (33) - 33 = quality 0
    // ASCII '~' (126) - 33 = quality 93 (max)
    create_paired_bam(&input_bam, vec![("read1", "ACGT", "!I?~", "TGCA", "IIII", false)]);

    let output = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args(["fastq", "-i", input_bam.to_str().unwrap()])
        .output()
        .expect("Failed to run fastq command");

    assert!(output.status.success());

    fs::write(&output_fq, &output.stdout).expect("Failed to write output");

    let records = parse_fastq_records(&output_fq);
    assert_eq!(records.len(), 2);

    // Quality should be preserved as Phred+33 ASCII
    assert_eq!(records[0].2, "!I?~", "Quality scores should be preserved");
}

/// Helper to create a BAM with secondary/supplementary reads for flag filtering tests.
fn create_bam_with_flags(path: &PathBuf) {
    let header = create_minimal_header("chr1", 10000);

    let file = fs::File::create(path).expect("Failed to create BAM file");
    let mut writer = bam::io::Writer::new(file);
    writer.write_header(&header).expect("Failed to write header");

    // Create a primary read
    let primary = {
        let mut b = SamBuilder::new();
        b.read_name(b"primary")
            .sequence(b"ACGT")
            .qualities(&[30, 30, 30, 30])
            .ref_id(0)
            .pos(99)
            .mapq(60);
        b.build()
    };
    writer
        .write_alignment_record(&header, &to_record_buf(&primary))
        .expect("Failed to write primary");

    // Create a secondary read (flag 0x100)
    let secondary = {
        let mut b = SamBuilder::new();
        b.read_name(b"secondary")
            .sequence(b"TGCA")
            .qualities(&[30, 30, 30, 30])
            .flags(flags::SECONDARY)
            .ref_id(0)
            .pos(99)
            .mapq(60);
        b.build()
    };
    writer
        .write_alignment_record(&header, &to_record_buf(&secondary))
        .expect("Failed to write secondary");

    // Create a supplementary read (flag 0x800)
    let supplementary = {
        let mut b = SamBuilder::new();
        b.read_name(b"supplementary")
            .sequence(b"GGGG")
            .qualities(&[30, 30, 30, 30])
            .flags(flags::SUPPLEMENTARY)
            .ref_id(0)
            .pos(99)
            .mapq(60);
        b.build()
    };
    writer
        .write_alignment_record(&header, &to_record_buf(&supplementary))
        .expect("Failed to write supplementary");

    writer.try_finish().expect("Failed to finish BAM");
}

/// Test flag filtering with -F option.
#[test]
fn test_fastq_exclude_flags() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_fq = temp_dir.path().join("output.fq");

    create_bam_with_flags(&input_bam);

    // Run with default flags (excludes secondary and supplementary)
    let output = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args(["fastq", "-i", input_bam.to_str().unwrap()])
        .output()
        .expect("Failed to run fastq command");

    assert!(
        output.status.success(),
        "Fastq command failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    fs::write(&output_fq, &output.stdout).expect("Failed to write output");

    let records = parse_fastq_records(&output_fq);
    assert_eq!(
        records.len(),
        1,
        "Should only have primary read (secondary and supplementary excluded)"
    );
    assert!(records[0].0.starts_with("primary"));
}

/// Test with multiple threads.
#[test]
fn test_fastq_multithreaded() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_fq = temp_dir.path().join("output.fq");

    // Create a larger set of reads
    let read_pairs: Vec<(&str, &str, &str, &str, &str, bool)> = (0..10)
        .map(|i| {
            let name: &'static str = Box::leak(format!("read{i}").into_boxed_str());
            (name, "ACGTACGT", "IIIIIIII", "TGCATGCA", "IIIIIIII", false)
        })
        .collect();

    create_paired_bam(&input_bam, read_pairs);

    let output = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args(["fastq", "-i", input_bam.to_str().unwrap(), "-@", "4"])
        .output()
        .expect("Failed to run fastq command");

    assert!(
        output.status.success(),
        "Multithreaded fastq failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    fs::write(&output_fq, &output.stdout).expect("Failed to write output");

    let records = parse_fastq_records(&output_fq);
    assert_eq!(records.len(), 20, "Should have 20 FASTQ records (10 pairs)");
}

/// Test hex flag parsing.
#[test]
fn test_fastq_hex_flags() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_fq = temp_dir.path().join("output.fq");

    create_paired_bam(&input_bam, vec![("read1", "ACGT", "IIII", "TGCA", "IIII", false)]);

    // Use hex notation for flags
    let output = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "fastq",
            "-i",
            input_bam.to_str().unwrap(),
            "-F",
            "0x900", // hex notation
        ])
        .output()
        .expect("Failed to run fastq command");

    assert!(output.status.success());

    fs::write(&output_fq, &output.stdout).expect("Failed to write output");

    let records = parse_fastq_records(&output_fq);
    assert_eq!(records.len(), 2);
}

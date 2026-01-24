//! Integration tests for the extract command.
//!
//! Tests end-to-end extract workflows with different compression formats.

use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

use flate2::Compression;
use flate2::write::GzEncoder;
use noodles::bam;
use noodles::sam::alignment::RecordBuf;
use noodles_bgzf::io::Writer as BgzfWriter;
use tempfile::TempDir;

/// Type alias for FASTQ test records (name, sequence, quality).
type FastqRecords = Vec<(&'static str, &'static str, &'static str)>;

// ============================================================================
// Helper Functions
// ============================================================================

/// Create a plain (uncompressed) FASTQ file.
fn create_plain_fastq(dir: &TempDir, name: &str, records: &[(&str, &str, &str)]) -> PathBuf {
    let path = dir.path().join(name);
    let mut file = File::create(&path).unwrap();
    for (name, seq, qual) in records {
        writeln!(file, "@{name}").unwrap();
        writeln!(file, "{seq}").unwrap();
        writeln!(file, "+").unwrap();
        writeln!(file, "{qual}").unwrap();
    }
    path
}

/// Create a gzip-compressed FASTQ file.
fn create_gzip_fastq(dir: &TempDir, name: &str, records: &[(&str, &str, &str)]) -> PathBuf {
    let path = dir.path().join(name);
    let file = File::create(&path).unwrap();
    let mut encoder = GzEncoder::new(file, Compression::default());
    for (name, seq, qual) in records {
        writeln!(encoder, "@{name}").unwrap();
        writeln!(encoder, "{seq}").unwrap();
        writeln!(encoder, "+").unwrap();
        writeln!(encoder, "{qual}").unwrap();
    }
    encoder.finish().unwrap();
    path
}

/// Create a BGZF-compressed FASTQ file.
fn create_bgzf_fastq(dir: &TempDir, name: &str, records: &[(&str, &str, &str)]) -> PathBuf {
    let path = dir.path().join(name);
    let file = File::create(&path).unwrap();
    let mut writer = BgzfWriter::new(file);
    for (name, seq, qual) in records {
        writeln!(writer, "@{name}").unwrap();
        writeln!(writer, "{seq}").unwrap();
        writeln!(writer, "+").unwrap();
        writeln!(writer, "{qual}").unwrap();
    }
    writer.finish().unwrap();
    path
}

/// Read BAM records from a file.
fn read_bam_records(path: &PathBuf) -> Vec<RecordBuf> {
    let mut reader = File::open(path).map(bam::io::Reader::new).unwrap();
    let header = reader.read_header().unwrap();
    reader.record_bufs(&header).map(|r| r.expect("Failed to read BAM record")).collect()
}

/// Standard test records for paired-end FASTQs.
fn paired_end_records() -> (FastqRecords, FastqRecords) {
    let r1 = vec![
        ("read1", "AAAAACGTACGTAAAA", "IIIIIIIIIIIIIIII"),
        ("read2", "TTTTTTTTTTTTTTTT", "IIIIIIIIIIIIIIII"),
        ("read3", "CCCCGGGGAAAATTTT", "IIIIIIIIIIIIIIII"),
    ];
    let r2 = vec![
        ("read1", "GGGGCGTACGTACCCC", "IIIIIIIIIIIIIIII"),
        ("read2", "AAAAAAAAAAAAGGGG", "IIIIIIIIIIIIIIII"),
        ("read3", "TTTTCCCCGGGGAAAA", "IIIIIIIIIIIIIIII"),
    ];
    (r1, r2)
}

// ============================================================================
// Gzip Compression Tests
// ============================================================================

#[test]
fn test_extract_gzip_single_threaded() {
    let tmp = TempDir::new().unwrap();
    let (r1_records, r2_records) = paired_end_records();

    let r1 = create_gzip_fastq(&tmp, "r1.fq.gz", &r1_records);
    let r2 = create_gzip_fastq(&tmp, "r2.fq.gz", &r2_records);
    let output = tmp.path().join("output.bam");

    // Run extract command (single-threaded)
    let status = std::process::Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "extract",
            "--inputs",
            r1.to_str().unwrap(),
            r2.to_str().unwrap(),
            "--output",
            output.to_str().unwrap(),
            "--read-structures",
            "5M+T",
            "5M+T",
            "--sample",
            "test_sample",
            "--library",
            "test_library",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to execute extract command");

    assert!(status.success(), "Extract command failed");

    // Verify output
    let records = read_bam_records(&output);
    assert_eq!(records.len(), 6, "Should have 6 records (3 pairs × 2 reads)");
}

#[test]
fn test_extract_gzip_multithreaded() {
    let tmp = TempDir::new().unwrap();
    let (r1_records, r2_records) = paired_end_records();

    let r1 = create_gzip_fastq(&tmp, "r1.fq.gz", &r1_records);
    let r2 = create_gzip_fastq(&tmp, "r2.fq.gz", &r2_records);
    let output = tmp.path().join("output.bam");

    // Run extract command with --threads 4
    let status = std::process::Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "extract",
            "--inputs",
            r1.to_str().unwrap(),
            r2.to_str().unwrap(),
            "--output",
            output.to_str().unwrap(),
            "--read-structures",
            "5M+T",
            "5M+T",
            "--sample",
            "test_sample",
            "--library",
            "test_library",
            "--threads",
            "4",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to execute extract command");

    assert!(status.success(), "Extract command with --threads failed");

    // Verify output
    let records = read_bam_records(&output);
    assert_eq!(records.len(), 6, "Should have 6 records (3 pairs × 2 reads)");
}

#[test]
fn test_extract_gzip_threads_mode() {
    let tmp = TempDir::new().unwrap();
    let (r1_records, r2_records) = paired_end_records();

    let r1 = create_gzip_fastq(&tmp, "r1.fq.gz", &r1_records);
    let r2 = create_gzip_fastq(&tmp, "r2.fq.gz", &r2_records);
    let output = tmp.path().join("output.bam");

    // Run extract command with --threads 4
    let status = std::process::Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "extract",
            "--inputs",
            r1.to_str().unwrap(),
            r2.to_str().unwrap(),
            "--output",
            output.to_str().unwrap(),
            "--read-structures",
            "5M+T",
            "5M+T",
            "--sample",
            "test_sample",
            "--library",
            "test_library",
            "--threads",
            "4",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to execute extract command");

    assert!(status.success(), "Extract command with --threads failed");

    // Verify output
    let records = read_bam_records(&output);
    assert_eq!(records.len(), 6, "Should have 6 records (3 pairs × 2 reads)");
}

// ============================================================================
// BGZF Compression Tests
// ============================================================================

#[test]
fn test_extract_bgzf_single_threaded() {
    let tmp = TempDir::new().unwrap();
    let (r1_records, r2_records) = paired_end_records();

    let r1 = create_bgzf_fastq(&tmp, "r1.fq.bgz", &r1_records);
    let r2 = create_bgzf_fastq(&tmp, "r2.fq.bgz", &r2_records);
    let output = tmp.path().join("output.bam");

    // Run extract command (single-threaded)
    let status = std::process::Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "extract",
            "--inputs",
            r1.to_str().unwrap(),
            r2.to_str().unwrap(),
            "--output",
            output.to_str().unwrap(),
            "--read-structures",
            "5M+T",
            "5M+T",
            "--sample",
            "test_sample",
            "--library",
            "test_library",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to execute extract command");

    assert!(status.success(), "Extract command with BGZF input failed");

    // Verify output
    let records = read_bam_records(&output);
    assert_eq!(records.len(), 6, "Should have 6 records (3 pairs × 2 reads)");
}

#[test]
fn test_extract_bgzf_multithreaded() {
    let tmp = TempDir::new().unwrap();
    let (r1_records, r2_records) = paired_end_records();

    let r1 = create_bgzf_fastq(&tmp, "r1.fq.bgz", &r1_records);
    let r2 = create_bgzf_fastq(&tmp, "r2.fq.bgz", &r2_records);
    let output = tmp.path().join("output.bam");

    // Run extract command with --threads 4
    let status = std::process::Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "extract",
            "--inputs",
            r1.to_str().unwrap(),
            r2.to_str().unwrap(),
            "--output",
            output.to_str().unwrap(),
            "--read-structures",
            "5M+T",
            "5M+T",
            "--sample",
            "test_sample",
            "--library",
            "test_library",
            "--threads",
            "4",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to execute extract command");

    assert!(status.success(), "Extract command with BGZF + --threads failed");

    // Verify output
    let records = read_bam_records(&output);
    assert_eq!(records.len(), 6, "Should have 6 records (3 pairs × 2 reads)");
}

#[test]
fn test_extract_bgzf_threads_mode() {
    let tmp = TempDir::new().unwrap();
    let (r1_records, r2_records) = paired_end_records();

    let r1 = create_bgzf_fastq(&tmp, "r1.fq.bgz", &r1_records);
    let r2 = create_bgzf_fastq(&tmp, "r2.fq.bgz", &r2_records);
    let output = tmp.path().join("output.bam");

    // Run extract command with --threads 4
    let status = std::process::Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "extract",
            "--inputs",
            r1.to_str().unwrap(),
            r2.to_str().unwrap(),
            "--output",
            output.to_str().unwrap(),
            "--read-structures",
            "5M+T",
            "5M+T",
            "--sample",
            "test_sample",
            "--library",
            "test_library",
            "--threads",
            "4",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to execute extract command");

    assert!(status.success(), "Extract command with BGZF + --threads failed");

    // Verify output
    let records = read_bam_records(&output);
    assert_eq!(records.len(), 6, "Should have 6 records (3 pairs × 2 reads)");
}

// ============================================================================
// Content Verification Tests
// ============================================================================

#[test]
fn test_extract_gzip_verifies_umi_extraction() {
    use noodles::sam::alignment::record::data::field::Tag;

    let tmp = TempDir::new().unwrap();

    // Create FASTQs with known UMI sequences (first 5 bases)
    // UMI = first 5 bases, Template = remaining bases
    let r1 = create_gzip_fastq(
        &tmp,
        "r1.fq.gz",
        &[
            ("read1", "ACGTATTTTTTT", "IIIIIIIIIIII"), // UMI: ACGTA, Template: TTTTTTT
        ],
    );
    let r2 = create_gzip_fastq(
        &tmp,
        "r2.fq.gz",
        &[
            ("read1", "TGCATAAAAAAA", "IIIIIIIIIIII"), // UMI: TGCAT, Template: AAAAAAA
        ],
    );
    let output = tmp.path().join("output.bam");

    let status = std::process::Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "extract",
            "--inputs",
            r1.to_str().unwrap(),
            r2.to_str().unwrap(),
            "--output",
            output.to_str().unwrap(),
            "--read-structures",
            "5M+T",
            "5M+T",
            "--sample",
            "test_sample",
            "--library",
            "test_library",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to execute extract command");

    assert!(status.success());

    let records = read_bam_records(&output);
    assert_eq!(records.len(), 2);

    // Verify UMI was extracted and stored in RX tag
    for record in &records {
        let rx_tag = record.data().get(&Tag::from(*b"RX"));
        assert!(rx_tag.is_some(), "RX tag should be present");
    }

    // Verify template bases (remaining after UMI extraction)
    // R1: ACGTA removed, leaving TTTTTTT
    // R2: TGCAT removed, leaving AAAAAAA
    let r1_record = &records[0];
    let r2_record = &records[1];

    assert_eq!(r1_record.sequence().as_ref(), b"TTTTTTT", "R1 template should be TTTTTTT");
    assert_eq!(r2_record.sequence().as_ref(), b"AAAAAAA", "R2 template should be AAAAAAA");
}

#[test]
fn test_extract_bgzf_verifies_umi_extraction() {
    let tmp = TempDir::new().unwrap();

    // Create FASTQs with known UMI sequences (first 5 bases)
    let r1 = create_bgzf_fastq(
        &tmp,
        "r1.fq.bgz",
        &[
            ("read1", "ACGTATTTTTTT", "IIIIIIIIIIII"), // UMI: ACGTA, Template: TTTTTTT
        ],
    );
    let r2 = create_bgzf_fastq(
        &tmp,
        "r2.fq.bgz",
        &[
            ("read1", "TGCATAAAAAAA", "IIIIIIIIIIII"), // UMI: TGCAT, Template: AAAAAAA
        ],
    );
    let output = tmp.path().join("output.bam");

    let status = std::process::Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "extract",
            "--inputs",
            r1.to_str().unwrap(),
            r2.to_str().unwrap(),
            "--output",
            output.to_str().unwrap(),
            "--read-structures",
            "5M+T",
            "5M+T",
            "--sample",
            "test_sample",
            "--library",
            "test_library",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to execute extract command");

    assert!(status.success());

    let records = read_bam_records(&output);
    assert_eq!(records.len(), 2);

    // Verify template bases
    let r1_record = &records[0];
    let r2_record = &records[1];

    assert_eq!(r1_record.sequence().as_ref(), b"TTTTTTT", "R1 template should be TTTTTTT");
    assert_eq!(r2_record.sequence().as_ref(), b"AAAAAAA", "R2 template should be AAAAAAA");
}

// ============================================================================
// Mixed Compression Tests
// ============================================================================

#[test]
fn test_extract_mixed_gzip_and_bgzf() {
    let tmp = TempDir::new().unwrap();
    let (r1_records, r2_records) = paired_end_records();

    // R1 is gzip, R2 is BGZF
    let r1 = create_gzip_fastq(&tmp, "r1.fq.gz", &r1_records);
    let r2 = create_bgzf_fastq(&tmp, "r2.fq.bgz", &r2_records);
    let output = tmp.path().join("output.bam");

    let status = std::process::Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "extract",
            "--inputs",
            r1.to_str().unwrap(),
            r2.to_str().unwrap(),
            "--output",
            output.to_str().unwrap(),
            "--read-structures",
            "5M+T",
            "5M+T",
            "--sample",
            "test_sample",
            "--library",
            "test_library",
            "--threads",
            "4",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to execute extract command");

    assert!(status.success(), "Extract with mixed gzip/BGZF should succeed");

    let records = read_bam_records(&output);
    assert_eq!(records.len(), 6, "Should have 6 records");
}

#[test]
fn test_extract_plain_and_compressed() {
    let tmp = TempDir::new().unwrap();
    let (r1_records, r2_records) = paired_end_records();

    // R1 is plain, R2 is gzip
    let r1 = create_plain_fastq(&tmp, "r1.fq", &r1_records);
    let r2 = create_gzip_fastq(&tmp, "r2.fq.gz", &r2_records);
    let output = tmp.path().join("output.bam");

    let status = std::process::Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "extract",
            "--inputs",
            r1.to_str().unwrap(),
            r2.to_str().unwrap(),
            "--output",
            output.to_str().unwrap(),
            "--read-structures",
            "5M+T",
            "5M+T",
            "--sample",
            "test_sample",
            "--library",
            "test_library",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to execute extract command");

    assert!(status.success(), "Extract with plain/gzip mix should succeed");

    let records = read_bam_records(&output);
    assert_eq!(records.len(), 6, "Should have 6 records");
}

// ============================================================================
// Parallel Parse Output Equivalence Tests
// ============================================================================

/// Test that parallel parse produces identical output across different thread counts.
///
/// This verifies that the parallel parse optimization (the key t8 scaling fix)
/// produces correct, deterministic output regardless of thread count.
#[test]
fn test_parallel_parse_output_equivalence_across_threads() {
    let tmp = TempDir::new().unwrap();

    // Create larger test data to exercise parallel parsing
    let records: Vec<(&str, &str, &str)> = (0..100)
        .map(|i| {
            // Leak strings to get static references (fine for tests)
            let name: &'static str = Box::leak(format!("read{i}").into_boxed_str());
            let seq: &'static str = Box::leak(format!("ACGT{i:016}").into_boxed_str());
            let qual: &'static str = Box::leak("I".repeat(20).into_boxed_str());
            (name, seq, qual)
        })
        .collect();

    let r1 = create_gzip_fastq(&tmp, "r1.fq.gz", &records);
    let r2 = create_gzip_fastq(&tmp, "r2.fq.gz", &records);

    // Run with different thread counts and compare outputs
    let mut outputs: Vec<Vec<RecordBuf>> = Vec::new();

    for threads in [1, 2, 4, 8] {
        let output = tmp.path().join(format!("output_t{threads}.bam"));

        let status = std::process::Command::new(env!("CARGO_BIN_EXE_fgumi"))
            .args([
                "extract",
                "--inputs",
                r1.to_str().unwrap(),
                r2.to_str().unwrap(),
                "--output",
                output.to_str().unwrap(),
                "--read-structures",
                "5M+T",
                "5M+T",
                "--sample",
                "test",
                "--library",
                "test",
                "--threads",
                &threads.to_string(),
                "--compression-level",
                "1",
            ])
            .status()
            .expect("Failed to execute extract command");

        assert!(status.success(), "Extract with {threads} threads failed");

        let records = read_bam_records(&output);
        outputs.push(records);
    }

    // Verify all outputs have the same number of records
    let expected_count = outputs[0].len();
    for (i, output) in outputs.iter().enumerate() {
        assert_eq!(
            output.len(),
            expected_count,
            "Thread count {} produced different record count",
            [1, 2, 4, 8][i]
        );
    }

    // Verify record sequences match across all thread counts
    for record_idx in 0..expected_count {
        let expected_seq = outputs[0][record_idx].sequence().as_ref();
        for (thread_idx, output) in outputs.iter().enumerate().skip(1) {
            let actual_seq = output[record_idx].sequence().as_ref();
            assert_eq!(
                expected_seq,
                actual_seq,
                "Record {record_idx} sequence mismatch between t1 and t{}",
                [1, 2, 4, 8][thread_idx]
            );
        }
    }
}

/// Test that running the same input multiple times produces identical output.
///
/// This verifies determinism of the parallel parse pipeline.
#[test]
fn test_parallel_parse_determinism() {
    let tmp = TempDir::new().unwrap();
    let (r1_records, r2_records) = paired_end_records();

    let r1 = create_gzip_fastq(&tmp, "r1.fq.gz", &r1_records);
    let r2 = create_gzip_fastq(&tmp, "r2.fq.gz", &r2_records);

    let mut outputs: Vec<Vec<RecordBuf>> = Vec::new();

    // Run 3 times with 4 threads
    for run in 0..3 {
        let output = tmp.path().join(format!("output_run{run}.bam"));

        let status = std::process::Command::new(env!("CARGO_BIN_EXE_fgumi"))
            .args([
                "extract",
                "--inputs",
                r1.to_str().unwrap(),
                r2.to_str().unwrap(),
                "--output",
                output.to_str().unwrap(),
                "--read-structures",
                "5M+T",
                "5M+T",
                "--sample",
                "test",
                "--library",
                "test",
                "--threads",
                "4",
                "--compression-level",
                "1",
            ])
            .status()
            .expect("Failed to execute extract command");

        assert!(status.success(), "Run {run} failed");
        outputs.push(read_bam_records(&output));
    }

    // Verify all runs produce identical output
    for run in 1..3 {
        assert_eq!(
            outputs[0].len(),
            outputs[run].len(),
            "Run {run} produced different record count"
        );

        for (i, (r0, r_run)) in outputs[0].iter().zip(outputs[run].iter()).enumerate() {
            assert_eq!(
                r0.sequence().as_ref(),
                r_run.sequence().as_ref(),
                "Record {i} sequence differs between run 0 and run {run}"
            );
            assert_eq!(
                r0.quality_scores().as_ref(),
                r_run.quality_scores().as_ref(),
                "Record {i} quality differs between run 0 and run {run}"
            );
        }
    }
}

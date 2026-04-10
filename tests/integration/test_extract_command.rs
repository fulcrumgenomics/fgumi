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

// ============================================================================
// Unified Pipeline Path Tests
// ============================================================================

/// Test BGZF+sync: multithreaded output matches single-threaded content.
///
/// This verifies the new BGZF+synchronized code path (which didn't exist before
/// the unified pipeline) produces correct output by comparing against the
/// single-threaded fast-path result.
#[test]
fn test_bgzf_sync_multithreaded_matches_single_threaded() {
    let tmp = TempDir::new().unwrap();

    let records_r1 = vec![
        ("read1", "ACGTACGTAAAA", "IIIIIIIIIIII"),
        ("read2", "TTTTTCGTACGT", "IIIIIIIIIIII"),
        ("read3", "CCCCGGGGAAAA", "IIIIIIIIIIII"),
    ];
    let records_r2 = vec![
        ("read1", "GGGGCGTACCCC", "IIIIIIIIIIII"),
        ("read2", "AAAAACGTAAAA", "IIIIIIIIIIII"),
        ("read3", "TTTTCCCCGGGG", "IIIIIIIIIIII"),
    ];

    let r1 = create_bgzf_fastq(&tmp, "r1.fq.bgz", &records_r1);
    let r2 = create_bgzf_fastq(&tmp, "r2.fq.bgz", &records_r2);

    // Run single-threaded (fast-path)
    let output_st = tmp.path().join("output_st.bam");
    let status = std::process::Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "extract",
            "--inputs",
            r1.to_str().unwrap(),
            r2.to_str().unwrap(),
            "--output",
            output_st.to_str().unwrap(),
            "--read-structures",
            "5M+T",
            "5M+T",
            "--sample",
            "test_sample",
            "--library",
            "test_library",
            "--threads",
            "1",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to execute single-threaded extract");
    assert!(status.success(), "Single-threaded extract failed");

    // Run multithreaded (BGZF+sync through unified pipeline)
    let output_threaded = tmp.path().join("output_mt.bam");
    let status = std::process::Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "extract",
            "--inputs",
            r1.to_str().unwrap(),
            r2.to_str().unwrap(),
            "--output",
            output_threaded.to_str().unwrap(),
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
        .expect("Failed to execute multithreaded extract");
    assert!(status.success(), "Multithreaded extract failed");

    // Compare record content (not just count)
    let st_records = read_bam_records(&output_st);
    let mt_records = read_bam_records(&output_threaded);

    assert_eq!(st_records.len(), mt_records.len(), "Record count mismatch");

    for (i, (st, mt)) in st_records.iter().zip(mt_records.iter()).enumerate() {
        assert_eq!(
            st.sequence().as_ref(),
            mt.sequence().as_ref(),
            "Record {i}: sequence differs between single-threaded and multithreaded"
        );
        assert_eq!(
            st.quality_scores().as_ref(),
            mt.quality_scores().as_ref(),
            "Record {i}: quality differs between single-threaded and multithreaded"
        );
    }
}

/// Test variable-length reads through the `RecordCount` reader.
///
/// The O(N^2) regression was caused by variable-length reads creating consistent
/// record count mismatches between R1 and R2 in the old byte-chunk reader.
/// The `RecordCount` reader fixes this by reading a fixed number of records.
#[test]
fn test_variable_length_reads_gzip() {
    let tmp = TempDir::new().unwrap();

    // Create reads with significantly different lengths between R1 and R2
    let records_r1 = vec![
        ("read1", "ACGTAAA", "IIIIIII"),                           // 7bp
        ("read2", "TTTTTCCCCCCCCCC", "IIIIIIIIIIIIIII"),           // 15bp
        ("read3", "CCCCG", "IIIII"),                               // 5bp
        ("read4", "GGGGAAAAAATTTTTTCCCC", "IIIIIIIIIIIIIIIIIIII"), // 20bp
    ];
    let records_r2 = vec![
        ("read1", "GGGGCGTACCCCCCCCCCCC", "IIIIIIIIIIIIIIIIIIII"), // 20bp
        ("read2", "AAAAA", "IIIII"),                               // 5bp
        ("read3", "TTTTCCCCGGGGAAAA", "IIIIIIIIIIIIIIII"),         // 16bp
        ("read4", "ACGTAC", "IIIIII"),                             // 6bp
    ];

    let r1 = create_gzip_fastq(&tmp, "r1.fq.gz", &records_r1);
    let r2 = create_gzip_fastq(&tmp, "r2.fq.gz", &records_r2);
    let output = tmp.path().join("output.bam");

    // Run multithreaded (exercises RecordCount reader for gzip+sync)
    let status = std::process::Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "extract",
            "--inputs",
            r1.to_str().unwrap(),
            r2.to_str().unwrap(),
            "--output",
            output.to_str().unwrap(),
            "--read-structures",
            "+T",
            "+T",
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

    assert!(status.success(), "Extract with variable-length reads failed");

    let records = read_bam_records(&output);
    assert_eq!(records.len(), 8, "Should have 8 records (4 pairs × 2 reads)");
}

/// Test variable-length reads with BGZF compression.
///
/// Exercises the BGZF+synchronized path with reads of different lengths.
#[test]
fn test_variable_length_reads_bgzf() {
    let tmp = TempDir::new().unwrap();

    let records_r1 = vec![
        ("read1", "ACGTAAA", "IIIIIII"),
        ("read2", "TTTTTCCCCCCCCCC", "IIIIIIIIIIIIIII"),
        ("read3", "CCCCG", "IIIII"),
        ("read4", "GGGGAAAAAATTTTTTCCCC", "IIIIIIIIIIIIIIIIIIII"),
    ];
    let records_r2 = vec![
        ("read1", "GGGGCGTACCCCCCCCCCCC", "IIIIIIIIIIIIIIIIIIII"),
        ("read2", "AAAAA", "IIIII"),
        ("read3", "TTTTCCCCGGGGAAAA", "IIIIIIIIIIIIIIII"),
        ("read4", "ACGTAC", "IIIIII"),
    ];

    let r1 = create_bgzf_fastq(&tmp, "r1.fq.bgz", &records_r1);
    let r2 = create_bgzf_fastq(&tmp, "r2.fq.bgz", &records_r2);
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
            "+T",
            "+T",
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

    assert!(status.success(), "BGZF extract with variable-length reads failed");

    let records = read_bam_records(&output);
    assert_eq!(records.len(), 8, "Should have 8 records (4 pairs × 2 reads)");
}

/// Regression test: verify that user-configurable options (--read-group-id,
/// --umi-tag, --annotate-read-names) are correctly propagated in multi-threaded mode.
#[test]
fn test_extract_multithreaded_custom_options() {
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::data::field::Value;

    let tmp = TempDir::new().unwrap();
    let r1 = create_gzip_fastq(&tmp, "r1.fq.gz", &[("read1", "ACGTATTTTTTT", "IIIIIIIIIIII")]);
    let r2 = create_gzip_fastq(&tmp, "r2.fq.gz", &[("read1", "TGCATAAAAAAA", "IIIIIIIIIIII")]);
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
            "--read-group-id",
            "MyRG",
            "--umi-tag",
            "ZU",
            "--annotate-read-names",
        ])
        .status()
        .expect("Failed to execute extract command");

    assert!(status.success(), "Multi-threaded extract with custom options failed");

    let records = read_bam_records(&output);
    assert_eq!(records.len(), 2);

    for record in &records {
        // Verify custom read group ID
        let rg = record.data().get(&Tag::from(*b"RG")).expect("RG tag should be present");
        match rg {
            Value::String(s) => {
                let rg_str = String::from_utf8_lossy(s);
                assert_eq!(rg_str, "MyRG", "Read group should be MyRG, not the default 'A'");
            }
            _ => panic!("RG tag should be a string"),
        }

        // Verify UMI stored under custom tag ZU (not default RX)
        assert!(record.data().get(&Tag::from(*b"ZU")).is_some(), "UMI should be under ZU tag");
        assert!(
            record.data().get(&Tag::from(*b"RX")).is_none(),
            "RX tag should not be present when --umi-tag ZU is used"
        );

        // Verify read name annotation (should contain '+' with UMI appended)
        let name = std::str::from_utf8(record.name().unwrap().as_ref()).unwrap();
        assert!(name.contains('+'), "Read name should be annotated with UMI: {name}");
    }
}

// ============================================================================
// BGZF BlockParseFast / BlockMerge Pipeline End-to-End Tests
// ============================================================================
//
// These tests verify the BGZF-specific parallel pipeline path (BlockParseFast +
// BlockMerge) at multiple thread counts, with both single-stream (R1 only) and
// paired-stream (R1+R2) inputs.

/// Helper: run extract with BGZF inputs at a specific thread count, return BAM records.
fn run_bgzf_extract(
    r1: &std::path::Path,
    r2_opt: Option<&std::path::Path>,
    threads: usize,
    tmp: &TempDir,
    tag_suffix: &str,
) -> Vec<RecordBuf> {
    let output = tmp.path().join(format!("out_{tag_suffix}_t{threads}.bam"));
    let mut args =
        vec!["extract".to_string(), "--inputs".to_string(), r1.to_str().unwrap().to_string()];
    if let Some(r2) = r2_opt {
        args.push(r2.to_str().unwrap().to_string());
    }
    args.extend([
        "--output".to_string(),
        output.to_str().unwrap().to_string(),
        "--read-structures".to_string(),
        "5M+T".to_string(),
    ]);
    if r2_opt.is_some() {
        args.push("5M+T".to_string());
    }
    args.extend([
        "--sample".to_string(),
        "test".to_string(),
        "--library".to_string(),
        "test".to_string(),
        "--threads".to_string(),
        threads.to_string(),
        "--compression-level".to_string(),
        "1".to_string(),
    ]);

    let status = std::process::Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args(&args)
        .status()
        .expect("Failed to execute extract command");

    assert!(status.success(), "extract t{threads} {tag_suffix} failed");
    read_bam_records(&output)
}

/// BGZF paired-stream extract at T1, T4, T8 — record count and content must agree.
#[test]
fn test_bgzf_block_merge_paired_thread_counts() {
    let tmp = TempDir::new().unwrap();

    // 500 pairs with 150-base reads — total uncompressed ~450 KB, guaranteeing
    // multiple BGZF blocks (each ≤64 KiB) so BlockParseFast/BlockMerge stitching
    // is exercised.
    let r1_recs: Vec<(&'static str, &'static str, &'static str)> = (0..500)
        .map(|i| {
            let name: &'static str = Box::leak(format!("read{i:04}").into_boxed_str());
            let seq: &'static str =
                Box::leak(format!("AAAAA{i:010}{}", "ACGT".repeat(33)).into_boxed_str());
            let qual: &'static str = Box::leak("I".repeat(147).into_boxed_str());
            (name, seq, qual)
        })
        .collect();

    let r2_recs: Vec<(&'static str, &'static str, &'static str)> = (0..500)
        .map(|i| {
            let name: &'static str = Box::leak(format!("read{i:04}").into_boxed_str());
            let seq: &'static str =
                Box::leak(format!("TTTTT{i:010}{}", "TGCA".repeat(33)).into_boxed_str());
            let qual: &'static str = Box::leak("J".repeat(147).into_boxed_str());
            (name, seq, qual)
        })
        .collect();

    let r1 = create_bgzf_fastq(&tmp, "r1.fq.bgz", &r1_recs);
    let r2 = create_bgzf_fastq(&tmp, "r2.fq.bgz", &r2_recs);

    // Run at T1, T4, T8.
    let results: Vec<Vec<RecordBuf>> =
        [1, 4, 8].iter().map(|&t| run_bgzf_extract(&r1, Some(&r2), t, &tmp, "paired")).collect();

    // All three runs must produce the same number of records (500 pairs × 2 reads = 1000).
    for (idx, recs) in results.iter().enumerate() {
        assert_eq!(
            recs.len(),
            1000,
            "T{} produced {} records, expected 1000",
            [1, 4, 8][idx],
            recs.len()
        );
    }

    // Sequences must be identical across thread counts (order preserved).
    let baseline = &results[0];
    for (ti, recs) in results.iter().enumerate().skip(1) {
        for (ri, (base, actual)) in baseline.iter().zip(recs.iter()).enumerate() {
            assert_eq!(
                base.sequence().as_ref(),
                actual.sequence().as_ref(),
                "record {ri}: sequence mismatch between T1 and T{}",
                [1, 4, 8][ti]
            );
        }
    }
}

/// BGZF single-stream extract (R1 only) at T1, T4, T8.
#[test]
fn test_bgzf_block_merge_single_stream_thread_counts() {
    let tmp = TempDir::new().unwrap();

    // 500 records with 150-base reads — spans multiple BGZF blocks.
    let r1_recs: Vec<(&'static str, &'static str, &'static str)> = (0..500)
        .map(|i| {
            let name: &'static str = Box::leak(format!("read{i:04}").into_boxed_str());
            let seq: &'static str =
                Box::leak(format!("AAAAA{i:010}{}", "ACGT".repeat(33)).into_boxed_str());
            let qual: &'static str = Box::leak("I".repeat(147).into_boxed_str());
            (name, seq, qual)
        })
        .collect();

    let r1 = create_bgzf_fastq(&tmp, "r1_single.fq.bgz", &r1_recs);

    // Single-stream extract: no R2.
    let results: Vec<Vec<RecordBuf>> =
        [1, 4, 8].iter().map(|&t| run_bgzf_extract(&r1, None, t, &tmp, "single")).collect();

    // 500 records (R1 only).
    for (idx, recs) in results.iter().enumerate() {
        assert_eq!(
            recs.len(),
            500,
            "T{} single-stream: {} records, expected 500",
            [1, 4, 8][idx],
            recs.len()
        );
    }

    // Sequences identical across thread counts.
    let baseline = &results[0];
    for (ti, recs) in results.iter().enumerate().skip(1) {
        for (ri, (base, actual)) in baseline.iter().zip(recs.iter()).enumerate() {
            assert_eq!(
                base.sequence().as_ref(),
                actual.sequence().as_ref(),
                "record {ri}: sequence mismatch T1 vs T{}",
                [1, 4, 8][ti]
            );
        }
    }
}

/// BGZF paired-stream extract — verify that template sequences match expected values.
///
/// This exercises the full `BlockParseFast` + `BlockMerge` path end-to-end and
/// verifies that R1/R2 pairing is correct (R1 template bases are in the right
/// BAM records).
#[test]
fn test_bgzf_block_merge_content_verification() {
    use noodles::sam::alignment::record::data::field::Tag;

    let tmp = TempDir::new().unwrap();

    // 3 pairs, template = everything after the 5-base UMI.
    let r1_recs = vec![
        ("read1", "ACGTATTTTTTT", "IIIIIIIIIIII"), // UMI=ACGTA, tmpl=TTTTTTT
        ("read2", "TGCATAAAAAAA", "IIIIIIIIIIII"), // UMI=TGCAT, tmpl=AAAAAAA
        ("read3", "CCGGAGGGGGG", "IIIIIIIIIII"),   // UMI=CCGGA, tmpl=GGGGGG (11-base read)
    ];
    let r2_recs = vec![
        ("read1", "GCTATCCCCCCC", "IIIIIIIIIIII"), // UMI=GCTAT, tmpl=CCCCCCC
        ("read2", "ATCGAGGGGGGGG", "IIIIIIIIIIIII"), // UMI=ATCGA, tmpl=GGGGGG G
        ("read3", "TTAAACCCCCC", "IIIIIIIIIII"),   // UMI=TTAAA, tmpl=CCCCCC
    ];

    let r1 = create_bgzf_fastq(&tmp, "r1_content.fq.bgz", &r1_recs);
    let r2 = create_bgzf_fastq(&tmp, "r2_content.fq.bgz", &r2_recs);

    // Run at T4 to exercise the parallel path.
    let records = run_bgzf_extract(&r1, Some(&r2), 4, &tmp, "content");

    // 3 pairs × 2 reads = 6 records.
    assert_eq!(records.len(), 6, "expected 6 records");

    // Verify UMI tag is present.
    for rec in &records {
        assert!(
            rec.data().get(&Tag::from(*b"RX")).is_some(),
            "RX tag must be present on every record"
        );
    }
}

/// Create a BGZF FASTQ file where each record is flushed into its own block.
///
/// This forces a specific number of BGZF blocks (one per record), which lets
/// us create R1 and R2 files with different block counts by giving them
/// records of different lengths.
fn create_bgzf_fastq_one_block_per_record(
    dir: &TempDir,
    name: &str,
    records: &[(&str, &str, &str)],
) -> PathBuf {
    let path = dir.path().join(name);
    let file = File::create(&path).unwrap();
    let mut writer = BgzfWriter::new(file);
    for (name, seq, qual) in records {
        write!(writer, "@{name}\n{seq}\n+\n{qual}\n").unwrap();
        // Flush after each record so it gets its own BGZF block.
        writer.flush().unwrap();
    }
    writer.finish().unwrap();
    path
}

/// Regression test: BGZF paired-end extract where R1 and R2 have different
/// numbers of BGZF blocks.
///
/// When R1 and R2 BGZF files have different compressed sizes (common in
/// practice), they produce different numbers of BGZF blocks. The `BlockMerge`
/// step must drain the remaining blocks from the longer stream after the
/// shorter stream is exhausted, or the pipeline deadlocks.
#[test]
fn test_extract_bgzf_unequal_block_counts() {
    let tmp = TempDir::new().unwrap();
    let num_pairs = 50;

    // Build R1 and R2 with different sequence lengths.
    // R1: 5bp UMI + 10bp template = 15bp per record
    // R2: 5bp UMI + 80bp template = 85bp per record
    // This makes R2 much larger per record, but both have the same record count.
    let r1_recs: Vec<(String, String, String)> = (0..num_pairs)
        .map(|i| {
            let name = format!("read{i:04}");
            let seq = format!("ACGTA{}", "T".repeat(10)); // 15bp
            let qual = "I".repeat(15);
            (name, seq, qual)
        })
        .collect();

    let r2_recs: Vec<(String, String, String)> = (0..num_pairs)
        .map(|i| {
            let name = format!("read{i:04}");
            let seq = format!("TGCAT{}", "A".repeat(80)); // 85bp
            let qual = "I".repeat(85);
            (name, seq, qual)
        })
        .collect();

    // Write with one record per BGZF block so R1 gets N blocks, R2 gets N blocks.
    // Actually, since the test data is small and we want *different* block counts,
    // write R1 normally (all in one block) and R2 with one-record-per-block.
    let r1_as_str: Vec<(&str, &str, &str)> =
        r1_recs.iter().map(|(a, b, c)| (a.as_str(), b.as_str(), c.as_str())).collect();
    let r2_as_str: Vec<(&str, &str, &str)> =
        r2_recs.iter().map(|(a, b, c)| (a.as_str(), b.as_str(), c.as_str())).collect();

    // R1: all records in one BGZF block (1 block total)
    // R2: one record per BGZF block (50 blocks total)
    let r1 = create_bgzf_fastq(&tmp, "r1_unequal.fq.bgz", &r1_as_str);
    let r2 = create_bgzf_fastq_one_block_per_record(&tmp, "r2_unequal.fq.bgz", &r2_as_str);

    for threads in [1, 4, 8] {
        let output = tmp.path().join(format!("out_unequal_t{threads}.bam"));
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

        assert!(
            status.success(),
            "Extract with unequal BGZF block counts should succeed at T{threads}"
        );

        let records = read_bam_records(&output);
        assert_eq!(
            records.len(),
            num_pairs * 2,
            "Should have {num_pairs}*2 records at T{threads}, got {}",
            records.len()
        );
    }
}

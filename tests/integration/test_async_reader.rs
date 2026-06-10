//! Integration tests for the FASTQ-side `--async-reader` flag on `fgumi extract`.
//!
//! Each test runs `fgumi extract` twice — once without `--async-reader` and once
//! with it — then asserts the outputs are identical record-by-record. This ensures
//! the prefetch reader is transparent across BGZF, gzip, and plain FASTQ inputs
//! (single- and multi-threaded).
//!
//! BAM-side `--async-reader` was deleted in Phase 1 (D2 of issue #330); the
//! corresponding tests were removed alongside the flag.

use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

use flate2::Compression;
use flate2::write::GzEncoder;
use noodles::bam;
use noodles::sam::alignment::RecordBuf;
use noodles_bgzf::io::Writer as BgzfWriter;
use tempfile::TempDir;

// ============================================================================
// Helpers
// ============================================================================

/// Type alias for FASTQ test records (name, sequence, quality).
type FastqRecords = Vec<(&'static str, &'static str, &'static str)>;

/// Standard paired-end test records with a 5-base UMI prefix on each read.
fn paired_end_records() -> (FastqRecords, FastqRecords) {
    let r1 = vec![
        ("read1", "AAAAACGTACGTAAAA", "IIIIIIIIIIIIIIII"),
        ("read2", "TTTTTCCCCCCCCCCC", "IIIIIIIIIIIIIIII"),
        ("read3", "CCCCCTTTTTAAAAAG", "IIIIIIIIIIIIIIII"),
    ];
    let r2 = vec![
        ("read1", "GGGGGCGTACGTCCCC", "IIIIIIIIIIIIIIII"),
        ("read2", "AAAAATTTTTTTTGGG", "IIIIIIIIIIIIIIII"),
        ("read3", "TTTTTCCCCGGGGAAA", "IIIIIIIIIIIIIIII"),
    ];
    (r1, r2)
}

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
fn read_bam_records(path: &Path) -> Vec<RecordBuf> {
    let mut reader = File::open(path).map(bam::io::Reader::new).unwrap();
    let header = reader.read_header().unwrap();
    reader.record_bufs(&header).map(|r| r.expect("Failed to read BAM record")).collect()
}

/// Read the BAM header from a file.
fn read_bam_header(path: &Path) -> noodles::sam::Header {
    let mut reader = File::open(path).map(bam::io::Reader::new).unwrap();
    reader.read_header().unwrap()
}

/// Assert that two BAM files contain identical records and identical headers
/// (ignoring `@PG` entries, which differ by command line).
fn assert_bam_records_equal(path1: &Path, path2: &Path) {
    // Compare headers with @PG entries stripped: those differ by command line
    // (--async-reader is in one and not the other) but everything else —
    // @HD, @SQ, @RG, @CO, sort order — must match.
    let mut header1 = read_bam_header(path1);
    let mut header2 = read_bam_header(path2);
    header1.programs_mut().as_mut().clear();
    header2.programs_mut().as_mut().clear();
    assert_eq!(header1, header2, "BAM headers differ (excluding @PG)");

    let records1 = read_bam_records(path1);
    let records2 = read_bam_records(path2);

    assert_eq!(
        records1.len(),
        records2.len(),
        "Record count differs: {} (without --async-reader) vs {} (with --async-reader)",
        records1.len(),
        records2.len(),
    );

    for (i, (r1, r2)) in records1.iter().zip(records2.iter()).enumerate() {
        assert_eq!(r1.name(), r2.name(), "Record {i} name mismatch");
        assert_eq!(r1.flags(), r2.flags(), "Record {i} flags mismatch");
        assert_eq!(
            r1.alignment_start(),
            r2.alignment_start(),
            "Record {i} alignment_start mismatch"
        );
        assert_eq!(
            r1.mapping_quality(),
            r2.mapping_quality(),
            "Record {i} mapping_quality mismatch"
        );
        assert_eq!(r1.cigar(), r2.cigar(), "Record {i} cigar mismatch");
        assert_eq!(r1.sequence(), r2.sequence(), "Record {i} sequence mismatch");
        assert_eq!(r1.quality_scores(), r2.quality_scores(), "Record {i} quality mismatch");
        assert_eq!(r1.data(), r2.data(), "Record {i} aux tags mismatch");
    }
}

/// Run `fgumi extract` with the given args, return the output BAM path.
fn run_extract(
    tmp: &TempDir,
    r1: &Path,
    r2: &Path,
    output_name: &str,
    threads: Option<usize>,
    async_reader: bool,
) -> PathBuf {
    let output = tmp.path().join(output_name);
    let mut args = vec![
        "extract".to_string(),
        "--inputs".into(),
        r1.to_str().unwrap().into(),
        r2.to_str().unwrap().into(),
        "--output".into(),
        output.to_str().unwrap().into(),
        "--read-structures".into(),
        "5M+T".into(),
        "5M+T".into(),
        "--sample".into(),
        "test_sample".into(),
        "--library".into(),
        "test_library".into(),
        "--compression-level".into(),
        "1".into(),
    ];
    if let Some(t) = threads {
        args.push("--threads".into());
        args.push(t.to_string());
    }
    if async_reader {
        args.push("--async-reader".into());
    }

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args(&args)
        .status()
        .expect("Failed to execute extract command");
    assert!(status.success(), "Extract command failed (async_reader={async_reader})");
    output
}
// ============================================================================
// Extract: BGZF FASTQ inputs
// ============================================================================

#[test]
fn test_extract_bgzf_async_reader_single_threaded() {
    let tmp = TempDir::new().unwrap();
    let (r1_recs, r2_recs) = paired_end_records();
    let r1 = create_bgzf_fastq(&tmp, "r1.fq.bgz", &r1_recs);
    let r2 = create_bgzf_fastq(&tmp, "r2.fq.bgz", &r2_recs);

    let baseline = run_extract(&tmp, &r1, &r2, "baseline.bam", None, false);
    let async_out = run_extract(&tmp, &r1, &r2, "async.bam", None, true);
    assert_bam_records_equal(&baseline, &async_out);
}

#[test]
fn test_extract_bgzf_async_reader_multithreaded() {
    let tmp = TempDir::new().unwrap();
    let (r1_recs, r2_recs) = paired_end_records();
    let r1 = create_bgzf_fastq(&tmp, "r1.fq.bgz", &r1_recs);
    let r2 = create_bgzf_fastq(&tmp, "r2.fq.bgz", &r2_recs);

    let baseline = run_extract(&tmp, &r1, &r2, "baseline.bam", Some(4), false);
    let async_out = run_extract(&tmp, &r1, &r2, "async.bam", Some(4), true);
    assert_bam_records_equal(&baseline, &async_out);
}

// ============================================================================
// Extract: gzip FASTQ inputs
// ============================================================================

#[test]
fn test_extract_gzip_async_reader_single_threaded() {
    let tmp = TempDir::new().unwrap();
    let (r1_recs, r2_recs) = paired_end_records();
    let r1 = create_gzip_fastq(&tmp, "r1.fq.gz", &r1_recs);
    let r2 = create_gzip_fastq(&tmp, "r2.fq.gz", &r2_recs);

    let baseline = run_extract(&tmp, &r1, &r2, "baseline.bam", None, false);
    let async_out = run_extract(&tmp, &r1, &r2, "async.bam", None, true);
    assert_bam_records_equal(&baseline, &async_out);
}

#[test]
fn test_extract_gzip_async_reader_multithreaded() {
    let tmp = TempDir::new().unwrap();
    let (r1_recs, r2_recs) = paired_end_records();
    let r1 = create_gzip_fastq(&tmp, "r1.fq.gz", &r1_recs);
    let r2 = create_gzip_fastq(&tmp, "r2.fq.gz", &r2_recs);

    let baseline = run_extract(&tmp, &r1, &r2, "baseline.bam", Some(4), false);
    let async_out = run_extract(&tmp, &r1, &r2, "async.bam", Some(4), true);
    assert_bam_records_equal(&baseline, &async_out);
}

// ============================================================================
// Extract: plain (uncompressed) FASTQ inputs
// ============================================================================

#[test]
fn test_extract_plain_async_reader_single_threaded() {
    let tmp = TempDir::new().unwrap();
    let (r1_recs, r2_recs) = paired_end_records();
    let r1 = create_plain_fastq(&tmp, "r1.fq", &r1_recs);
    let r2 = create_plain_fastq(&tmp, "r2.fq", &r2_recs);

    let baseline = run_extract(&tmp, &r1, &r2, "baseline.bam", None, false);
    let async_out = run_extract(&tmp, &r1, &r2, "async.bam", None, true);
    assert_bam_records_equal(&baseline, &async_out);
}

#[test]
fn test_extract_plain_async_reader_multithreaded() {
    let tmp = TempDir::new().unwrap();
    let (r1_recs, r2_recs) = paired_end_records();
    let r1 = create_plain_fastq(&tmp, "r1.fq", &r1_recs);
    let r2 = create_plain_fastq(&tmp, "r2.fq", &r2_recs);

    let baseline = run_extract(&tmp, &r1, &r2, "baseline.bam", Some(4), false);
    let async_out = run_extract(&tmp, &r1, &r2, "async.bam", Some(4), true);
    assert_bam_records_equal(&baseline, &async_out);
}

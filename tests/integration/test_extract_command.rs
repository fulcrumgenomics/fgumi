//! Integration tests for the extract command.
//!
//! Tests end-to-end extract workflows with different compression formats.

use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

use clap::Parser;
use fgumi_lib::commands::command::Command;
use fgumi_lib::commands::extract::Extract;
use flate2::Compression;
use flate2::write::GzEncoder;
use noodles::bam;
use noodles::sam::alignment::RecordBuf;
use noodles_bgzf::io::Writer as BgzfWriter;
use rstest::rstest;
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
fn read_bam_records(path: &Path) -> Vec<RecordBuf> {
    let mut reader = File::open(path).map(bam::io::Reader::new).unwrap();
    let header = reader.read_header().unwrap();
    reader.record_bufs(&header).map(|r| r.expect("Failed to read BAM record")).collect()
}

/// The number of records in `path` if it can be read as a BAM, or `None` if the
/// file is absent or does not parse as one. Used to inspect the partial output a
/// rejected run may leave behind without panicking on a truncated file.
fn try_bam_record_count(path: &Path) -> Option<usize> {
    let mut reader = File::open(path).map(bam::io::Reader::new).ok()?;
    let header = reader.read_header().ok()?;
    let mut count = 0;
    for record in reader.record_bufs(&header) {
        record.ok()?;
        count += 1;
    }
    Some(count)
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
    let cmd = Extract::try_parse_from([
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
    .expect("failed to parse extract args");
    cmd.execute("fgumi extract").expect("Extract command failed");

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
    let cmd = Extract::try_parse_from([
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
    .expect("failed to parse extract args");
    cmd.execute("fgumi extract").expect("Extract command with --threads failed");

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
    let cmd = Extract::try_parse_from([
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
    .expect("failed to parse extract args");
    cmd.execute("fgumi extract").expect("Extract command with --threads failed");

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
    let cmd = Extract::try_parse_from([
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
    .expect("failed to parse extract args");
    cmd.execute("fgumi extract").expect("Extract command with BGZF input failed");

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
    let cmd = Extract::try_parse_from([
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
    .expect("failed to parse extract args");
    cmd.execute("fgumi extract").expect("Extract command with BGZF + --threads failed");

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
    let cmd = Extract::try_parse_from([
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
    .expect("failed to parse extract args");
    cmd.execute("fgumi extract").expect("Extract command with BGZF + --threads failed");

    // Verify output
    let records = read_bam_records(&output);
    assert_eq!(records.len(), 6, "Should have 6 records (3 pairs × 2 reads)");
}

// ============================================================================
// Content Verification Tests
// ============================================================================

#[test]
fn test_extract_gzip_verifies_umi_extraction() {
    use fgumi_lib::sam::SamTag;
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

    let cmd = Extract::try_parse_from([
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
    .expect("failed to parse extract args");
    cmd.execute("fgumi extract").expect("Failed to execute extract command");

    let records = read_bam_records(&output);
    assert_eq!(records.len(), 2);

    // Verify UMI was extracted and stored in RX tag
    for record in &records {
        let rx_tag = record.data().get(&Tag::from(SamTag::RX));
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

    let cmd = Extract::try_parse_from([
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
    .expect("failed to parse extract args");
    cmd.execute("fgumi extract").expect("Failed to execute extract command");

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

    let cmd = Extract::try_parse_from([
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
    .expect("failed to parse extract args");
    cmd.execute("fgumi extract").expect("Extract with mixed gzip/BGZF should succeed");

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

    let cmd = Extract::try_parse_from([
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
    .expect("failed to parse extract args");
    cmd.execute("fgumi extract").expect("Extract with plain/gzip mix should succeed");

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

        let cmd = Extract::try_parse_from([
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
        .expect("failed to parse extract args");
        cmd.execute("fgumi extract")
            .unwrap_or_else(|e| panic!("Extract with {threads} threads failed: {e}"));

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

/// Extract's output must be independent of the worker count: a no-`--threads`
/// run (the chain at a single worker) must equal a `--threads N` run on decoded
/// records AND normalized header, with `--threads 1` explicitly covered. Since
/// the legacy serial path was retired, both invocations route through the
/// declarative chain builder; this pins single-worker output against
/// multi-worker output (the earlier across-threads test only compares
/// `--threads` values against each other, never the no-`--threads` default).
///
/// Fixed-value anchors (record count, each mate's template sequence, the combined
/// `RX`, and the `RG`) are pinned in addition to the cross-worker-count equality,
/// so a regression in the shared record-building code is caught even when both
/// runs still agree with each other. (Byte-parity against the pre-removal serial
/// binary lives in `test_extract_cutover_parity.rs`.)
#[rstest]
#[case::threads_one(1)]
#[case::threads_four(4)]
fn no_threads_matches_threaded(#[case] threads: usize) {
    use fgumi_lib::sam::SamTag;
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::data::field::Value;

    let tmp = TempDir::new().unwrap();
    // 5M+T over a 10 bp read: a 5 bp UMI then a 5 bp template. One paired template.
    let r1 = create_plain_fastq(&tmp, "r1.fq", &[("pair0", "AAAAACCCCC", "IIIIIIIIII")]);
    let r2 = create_plain_fastq(&tmp, "r2.fq", &[("pair0", "GGGGGTTTTT", "IIIIIIIIII")]);

    let run = |threads: Option<usize>, out: &Path| {
        let mut args: Vec<String> = vec![
            "extract".into(),
            "--inputs".into(),
            r1.to_str().unwrap().into(),
            r2.to_str().unwrap().into(),
            "--output".into(),
            out.to_str().unwrap().into(),
            "--read-structures".into(),
            "5M+T".into(),
            "5M+T".into(),
            "--sample".into(),
            "psample".into(),
            "--library".into(),
            "plib".into(),
            "--compression-level".into(),
            "1".into(),
        ];
        if let Some(t) = threads {
            args.push("--threads".into());
            args.push(t.to_string());
        }
        Extract::try_parse_from(args)
            .expect("failed to parse extract args")
            .execute("fgumi extract")
            .expect("extract failed");
    };

    let oracle_out = tmp.path().join("oracle.bam");
    let chain_out = tmp.path().join(format!("chain_t{threads}.bam"));
    run(None, &oracle_out); // chain at a single worker (no --threads)
    run(Some(threads), &chain_out); // chain (--threads N, incl. N == 1)

    let (oracle_hdr, oracle_recs) = crate::helpers::read_bam_output(&oracle_out);
    let (chain_hdr, chain_recs) = crate::helpers::read_bam_output(&chain_out);
    assert_eq!(chain_hdr, oracle_hdr, "normalized header parity (threads={threads})");
    assert_eq!(chain_recs, oracle_recs, "record parity (threads={threads})");

    // Fixed-value anchors on the single-worker output (which the threaded run now equals): one
    // record per mate, the template halves after the 5 bp UMI, the combined
    // R1-R2 UMI in RX, and the default read group on every record.
    assert_eq!(oracle_recs.len(), 2, "one BAM record per mate");
    assert_eq!(oracle_recs[0].sequence().as_ref(), b"CCCCC", "R1 template after 5M UMI");
    assert_eq!(oracle_recs[1].sequence().as_ref(), b"TTTTT", "R2 template after 5M UMI");
    let string_tag = |rec: &RecordBuf, tag: SamTag| -> Vec<u8> {
        match rec.data().get(&Tag::from(tag)) {
            Some(Value::String(s)) => s.to_vec(),
            other => panic!("expected a string {tag:?} tag, got {other:?}"),
        }
    };
    for rec in &oracle_recs {
        assert_eq!(string_tag(rec, SamTag::RX), b"AAAAA-GGGGG", "combined R1-R2 UMI in RX");
        assert_eq!(string_tag(rec, SamTag::RG), b"A", "default read group id");
    }
}

/// Cross-worker-count parity across the barcode / single-tag / annotate branches
/// of the record builder. `no_threads_matches_threaded` covers only RX/RG, and
/// every threading-parameterized tag test pins `store_cell_quals` / `single_tag`
/// / `annotate_read_names` / `store_sample_barcode_qualities` off, so nothing else
/// exercises those branches. A single-end `3C3B3M+T` read with every tag flag on
/// drives CB/CY, BC/QT, RX/QX, the `--single-tag` copy, and
/// `--annotate-read-names`; the single-worker (no-`--threads`) and multi-worker
/// (`--threads N`) runs must agree byte-for-byte, and the emitted tag values are
/// pinned so a shared-code regression is caught even when the two runs still agree.
#[rstest]
#[case::threads_one(1)]
#[case::threads_four(4)]
fn no_threads_matches_threaded_barcode_and_annotate_tags(#[case] threads: usize) {
    use fgumi_lib::sam::SamTag;
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::data::field::Value;

    let tmp = TempDir::new().unwrap();
    // 3C 3B 3M then template over a 13 bp read "AAACCCGGGTTTT": cell=AAA,
    // sample=CCC, umi=GGG, template=TTTT. (3 bp segments, so the pinned barcode/UMI
    // payloads are not 2-char tag-shaped literals the tag-literal lint would flag.)
    // The quality string is per-position distinct AND unambiguously Phred+33 (every
    // char < ASCII 64), so each segment's quality slice is pinnable to a known value
    // and cannot be shifted by Phred+64 encoding detection: cell qual "+,-",
    // sample qual "./0", umi qual "123", template qual "4567".
    let input = create_plain_fastq(&tmp, "r1.fq", &[("readX", "AAACCCGGGTTTT", "+,-./01234567")]);

    let run = |threads: Option<usize>, out: &Path| {
        let mut args: Vec<String> = vec![
            "extract".into(),
            "--inputs".into(),
            input.to_str().unwrap().into(),
            "--output".into(),
            out.to_str().unwrap().into(),
            "--read-structures".into(),
            "3C3B3M+T".into(),
            "--sample".into(),
            "s".into(),
            "--library".into(),
            "l".into(),
            "--store-umi-quals".into(),
            "--store-cell-quals".into(),
            "--store-sample-barcode-qualities".into(),
            "--annotate-read-names".into(),
            "--single-tag".into(),
            "MI".into(),
            "--compression-level".into(),
            "1".into(),
        ];
        if let Some(t) = threads {
            args.push("--threads".into());
            args.push(t.to_string());
        }
        Extract::try_parse_from(args)
            .expect("failed to parse extract args")
            .execute("fgumi extract")
            .expect("extract failed");
    };

    let oracle_out = tmp.path().join("oracle.bam");
    let chain_out = tmp.path().join(format!("chain_t{threads}.bam"));
    run(None, &oracle_out);
    run(Some(threads), &chain_out);

    let (oracle_hdr, oracle_recs) = crate::helpers::read_bam_output(&oracle_out);
    let (chain_hdr, chain_recs) = crate::helpers::read_bam_output(&chain_out);
    assert_eq!(chain_hdr, oracle_hdr, "header parity (threads={threads})");
    assert_eq!(chain_recs, oracle_recs, "record parity (threads={threads})");

    assert_eq!(oracle_recs.len(), 1, "one record for the single-end read");
    let rec = &oracle_recs[0];
    let string_tag = |rec: &RecordBuf, tag: SamTag| -> Vec<u8> {
        match rec.data().get(&Tag::from(tag)) {
            Some(Value::String(s)) => s.to_vec(),
            other => panic!("expected a string {tag:?} tag, got {other:?}"),
        }
    };
    assert_eq!(rec.sequence().as_ref(), b"TTTT", "template after 3C3B3M");
    assert_eq!(string_tag(rec, SamTag::CB), b"AAA", "cell barcode");
    assert_eq!(string_tag(rec, SamTag::BC), b"CCC", "sample barcode");
    assert_eq!(string_tag(rec, SamTag::RX), b"GGG", "UMI");
    assert_eq!(string_tag(rec, SamTag::MI), b"GGG", "--single-tag copy of the UMI");
    // The store_* flags are on, so the quality tags carry the per-segment quality
    // slices verbatim (the input is unambiguously Phred+33, so no encoding shift).
    assert_eq!(string_tag(rec, SamTag::CY), b"+,-", "cell-barcode quals (CY)");
    assert_eq!(string_tag(rec, SamTag::QT), b"./0", "sample-barcode quals (QT)");
    assert_eq!(string_tag(rec, SamTag::QX), b"123", "UMI quals (QX)");
    // --annotate-read-names appends `+<UMI>` to the read name.
    let name: &[u8] = rec.name().expect("read name must be present").as_ref();
    assert_eq!(name, &b"readX+GGG"[..], "annotated read name");
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

        let cmd = Extract::try_parse_from([
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
        .expect("failed to parse extract args");
        cmd.execute("fgumi extract").unwrap_or_else(|e| panic!("Run {run} failed: {e}"));
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

/// Test BGZF+sync: multi-worker output matches single-worker content.
///
/// This verifies the BGZF+synchronized code path produces correct output by
/// comparing a `--threads 1` (single-worker) chain run against a multi-worker
/// chain run over BGZF-compressed input.
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

    // Run at a single worker (--threads 1)
    let output_st = tmp.path().join("output_st.bam");
    let cmd = Extract::try_parse_from([
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
    .expect("failed to parse extract args");
    cmd.execute("fgumi extract").expect("Failed to execute single-worker extract");

    // Run multithreaded (BGZF+sync through unified pipeline)
    let output_threaded = tmp.path().join("output_mt.bam");
    let cmd = Extract::try_parse_from([
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
    .expect("failed to parse extract args");
    cmd.execute("fgumi extract").expect("Failed to execute multithreaded extract");

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
    let cmd = Extract::try_parse_from([
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
    .expect("failed to parse extract args");
    cmd.execute("fgumi extract").expect("Extract with variable-length reads failed");

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

    let cmd = Extract::try_parse_from([
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
    .expect("failed to parse extract args");
    cmd.execute("fgumi extract").expect("BGZF extract with variable-length reads failed");

    let records = read_bam_records(&output);
    assert_eq!(records.len(), 8, "Should have 8 records (4 pairs × 2 reads)");
}

/// Regression test: verify that user-configurable options (--read-group-id,
/// --annotate-read-names) are correctly propagated in multi-threaded mode.
#[test]
fn test_extract_multithreaded_custom_options() {
    use fgumi_lib::sam::SamTag;
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::data::field::Value;

    let tmp = TempDir::new().unwrap();
    let r1 = create_gzip_fastq(&tmp, "r1.fq.gz", &[("read1", "ACGTATTTTTTT", "IIIIIIIIIIII")]);
    let r2 = create_gzip_fastq(&tmp, "r2.fq.gz", &[("read1", "TGCATAAAAAAA", "IIIIIIIIIIII")]);
    let output = tmp.path().join("output.bam");

    let cmd = Extract::try_parse_from([
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
        "--annotate-read-names",
    ])
    .expect("failed to parse extract args");
    cmd.execute("fgumi extract").expect("Multi-threaded extract with custom options failed");

    let records = read_bam_records(&output);
    assert_eq!(records.len(), 2);

    for record in &records {
        // Verify custom read group ID
        let rg = record.data().get(&Tag::from(SamTag::RG)).expect("RG tag should be present");
        match rg {
            Value::String(s) => {
                let rg_str = String::from_utf8_lossy(s);
                assert_eq!(rg_str, "MyRG", "Read group should be MyRG, not the default 'A'");
            }
            _ => panic!("RG tag should be a string"),
        }

        // Verify UMI stored under standard RX tag
        assert!(record.data().get(&Tag::from(SamTag::RX)).is_some(), "UMI should be under RX tag");

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

    let cmd = Extract::try_parse_from(args.iter().map(String::as_str))
        .expect("failed to parse extract args");
    cmd.execute("fgumi extract")
        .unwrap_or_else(|e| panic!("extract t{threads} {tag_suffix} failed: {e}"));
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
    use fgumi_lib::sam::SamTag;
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
            rec.data().get(&Tag::from(SamTag::RX)).is_some(),
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
        let cmd = Extract::try_parse_from([
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
        .expect("failed to parse extract args");
        cmd.execute("fgumi extract").unwrap_or_else(|e| {
            panic!("Extract with unequal BGZF block counts should succeed at T{threads}: {e}")
        });

        let records = read_bam_records(&output);
        assert_eq!(
            records.len(),
            num_pairs * 2,
            "Should have {num_pairs}*2 records at T{threads}, got {}",
            records.len()
        );
    }
}

// ============================================================================
// Compression Level Boundary Tests (issue #360)
// ============================================================================

/// `--compression-level 0` must produce a valid, readable BAM through both the
/// single- and multi-threaded unified-pipeline paths. Byte-level correctness
/// (stored vs DEFLATE blocks) is asserted by the `InlineBgzfCompressor` unit
/// tests in `crates/fgumi-bgzf`.
#[rstest]
#[case::single_threaded(1)]
#[case::multi_threaded(4)]
fn test_extract_compression_level_zero(#[case] threads: usize) {
    let tmp = TempDir::new().unwrap();
    let (r1_records, r2_records) = paired_end_records();
    let r1 = create_gzip_fastq(&tmp, "r1.fq.gz", &r1_records);
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
            "--threads",
            &threads.to_string(),
            "--compression-level",
            "0",
        ])
        .status()
        .expect("Failed to execute extract command");

    assert!(status.success(), "Extract with --compression-level 0 failed at T{threads}");

    let records = read_bam_records(&output);
    assert_eq!(records.len(), 6, "Expected 6 records at T{threads}");
}

/// `--compression-level` values outside `0..=12` must be rejected by clap
/// before any work is done, with a range-validation error on stderr.
#[test]
fn test_extract_compression_level_out_of_range_rejected() {
    let tmp = TempDir::new().unwrap();
    let (r1_records, r2_records) = paired_end_records();
    let r1 = create_gzip_fastq(&tmp, "r1.fq.gz", &r1_records);
    let r2 = create_gzip_fastq(&tmp, "r2.fq.gz", &r2_records);
    let output = tmp.path().join("output.bam");

    let result = std::process::Command::new(env!("CARGO_BIN_EXE_fgumi"))
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
            "13",
        ])
        .output()
        .expect("Failed to execute extract command");

    assert!(!result.status.success(), "Extract with --compression-level 13 should have failed");
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(
        stderr.contains("0..=12") || stderr.contains("invalid value"),
        "Expected clap range-validation error mentioning `0..=12` or `invalid value`, \
         got stderr: {stderr}"
    );
}

// ============================================================================
// Mismatched-Length FASTQ Pair (issue #773)
// ============================================================================

/// Build `count` synthetic paired-end FASTQ records named `read00000..`.
///
/// Sequences are constant; only the name distinguishes records, which is what
/// the identity assertions below key on.
fn numbered_records(count: usize) -> Vec<(String, String, String)> {
    (0..count)
        .map(|i| {
            (format!("read{i:06}"), "ACGTACGTACGTACGT".to_string(), "IIIIIIIIIIIIIIII".to_string())
        })
        .collect()
}

/// Borrow a `Vec<(String, String, String)>` as the `&str` triples the FASTQ
/// writers take.
fn as_str_records(records: &[(String, String, String)]) -> Vec<(&str, &str, &str)> {
    records.iter().map(|(a, b, c)| (a.as_str(), b.as_str(), c.as_str())).collect()
}

/// The FASTQ compression formats `extract` dispatches on. Each selects a
/// different pipeline: BGZF uses `BlockParseFast`/`BlockMerge`, gzip and plain
/// use `FindBoundaries`/`Decode`.
#[derive(Clone, Copy, Debug)]
enum FastqFlavor {
    Plain,
    Gzip,
    Bgzf,
}

impl FastqFlavor {
    fn write(self, dir: &TempDir, name: &str, records: &[(&str, &str, &str)]) -> PathBuf {
        match self {
            Self::Plain => create_plain_fastq(dir, name, records),
            Self::Gzip => create_gzip_fastq(dir, name, records),
            Self::Bgzf => create_bgzf_fastq(dir, name, records),
        }
    }

    /// Write raw FASTQ `bytes` through this flavor's compressor. Unlike
    /// [`Self::write`], the caller controls the exact bytes, so a file that does
    /// not end in a newline can be produced.
    fn write_bytes(self, dir: &TempDir, name: &str, bytes: &[u8]) -> PathBuf {
        let path = dir.path().join(name);
        let file = File::create(&path).unwrap();
        match self {
            Self::Plain => {
                let mut file = file;
                file.write_all(bytes).unwrap();
            }
            Self::Gzip => {
                let mut encoder = GzEncoder::new(file, Compression::default());
                encoder.write_all(bytes).unwrap();
                encoder.finish().unwrap();
            }
            Self::Bgzf => {
                let mut writer = BgzfWriter::new(file);
                writer.write_all(bytes).unwrap();
                writer.finish().unwrap();
            }
        }
        path
    }
}

/// Serialize `records` as FASTQ text, optionally omitting the newline after the
/// final record. A file that does not end in a newline leaves its last record in
/// `suffix_bytes` on the BGZF path, which is the shape the no-trailing-newline
/// regression test exercises.
fn fastq_bytes(records: &[(&str, &str, &str)], trailing_newline: bool) -> Vec<u8> {
    let mut bytes = Vec::new();
    for (name, seq, qual) in records {
        bytes.extend_from_slice(format!("@{name}\n{seq}\n+\n{qual}\n").as_bytes());
    }
    if !trailing_newline && bytes.last() == Some(&b'\n') {
        bytes.pop();
    }
    bytes
}

/// Run `fgumi extract` over a paired FASTQ input, returning the command result.
fn run_extract_pair(r1: &Path, r2: &Path, output: &Path, threads: usize) -> anyhow::Result<()> {
    let cmd = Extract::try_parse_from([
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
        &threads.to_string(),
        "--compression-level",
        "1",
    ])
    .expect("failed to parse extract args");
    cmd.execute("fgumi extract")
}

/// Run `fgumi extract` over a paired FASTQ input with no `--threads` flag (the
/// chain at a single worker), returning the command result.
fn run_extract_pair_single_threaded(r1: &Path, r2: &Path, output: &Path) -> anyhow::Result<()> {
    let cmd = Extract::try_parse_from([
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
        "--compression-level",
        "1",
    ])
    .expect("failed to parse extract args");
    cmd.execute("fgumi extract")
}

/// Every read name in `path`, in file order (one entry per BAM record, so a
/// properly paired template contributes its name twice).
fn bam_read_names(path: &Path) -> Vec<String> {
    read_bam_records(path)
        .iter()
        .map(|record| {
            let name: &[u8] = record.name().expect("BAM record must have a name").as_ref();
            String::from_utf8(name.to_vec()).expect("read name must be UTF-8")
        })
        .collect()
}

/// Describe, by record identity, how a BAM written from a mismatched FASTQ pair
/// differs from the long input stream.
///
/// Used only to build the failure message when `extract` wrongly accepts such a
/// pair, so the failure names the records that were lost rather than reporting a
/// bare count mismatch.
fn describe_mismatched_output(output: &Path, long_stream: &[(String, String, String)]) -> String {
    let emitted = bam_read_names(output);
    let present: std::collections::HashSet<&str> = emitted.iter().map(String::as_str).collect();
    let missing: Vec<&str> = long_stream
        .iter()
        .map(|(name, _, _)| name.as_str())
        .filter(|name| !present.contains(name))
        .collect();
    format!(
        "{} of {} input records reached the output, {} were dropped (first missing: {:?}), \
         and {} records were emitted unpaired",
        present.len(),
        long_stream.len(),
        missing.len(),
        missing.first(),
        unpaired_record_count(output),
    )
}

/// Number of BAM records in `path` that are not flagged as paired (`0x1`).
fn unpaired_record_count(path: &Path) -> usize {
    read_bam_records(path).iter().filter(|record| !record.flags().is_segmented()).count()
}

/// A matched pair must survive extract with every template intact — asserted by
/// exact read-name identity, not by a count.
///
/// This is the control for [`test_extract_rejects_mismatched_fastq_pair`]: the
/// rejection added for issue #773 must not reject a well-formed input, and the
/// record set it emits must be exactly the input record set. `RECORD_COUNT`
/// straddles the per-batch record limit (200 at 1 thread, 800 at 4+) so the
/// final batch is partial — the shape in which the old truncation silently
/// dropped records.
#[rstest]
#[case::plain_t1(FastqFlavor::Plain, 1)]
#[case::plain_t4(FastqFlavor::Plain, 4)]
#[case::gzip_t1(FastqFlavor::Gzip, 1)]
#[case::gzip_t4(FastqFlavor::Gzip, 4)]
#[case::bgzf_t1(FastqFlavor::Bgzf, 1)]
#[case::bgzf_t4(FastqFlavor::Bgzf, 4)]
fn test_extract_matched_pair_emits_every_template(
    #[case] flavor: FastqFlavor,
    #[case] threads: usize,
) {
    const RECORD_COUNT: usize = 950;

    let tmp = TempDir::new().unwrap();
    let records = numbered_records(RECORD_COUNT);
    let as_str = as_str_records(&records);
    let r1 = flavor.write(&tmp, "r1.fq", &as_str);
    let r2 = flavor.write(&tmp, "r2.fq", &as_str);
    let output = tmp.path().join("matched.bam");

    run_extract_pair(&r1, &r2, &output, threads).unwrap_or_else(|e| {
        panic!("matched pair must extract cleanly ({flavor:?}, t{threads}): {e}")
    });

    // Identity: exactly two records per input template, one per mate, and the
    // emitted name multiset equals the input name multiset.
    let mut emitted = bam_read_names(&output);
    emitted.sort();
    let mut expected: Vec<String> =
        records.iter().flat_map(|(name, _, _)| [name.clone(), name.clone()]).collect();
    expected.sort();
    assert_eq!(
        emitted, expected,
        "emitted read names must equal the input read names ({flavor:?}, t{threads})"
    );
    assert_eq!(
        unpaired_record_count(&output),
        0,
        "every emitted record must be flagged paired ({flavor:?}, t{threads})"
    );
}

/// A well-formed matched pair whose FASTQ files do not end in a trailing newline
/// must still extract cleanly (#773).
///
/// `detect_suffix_start` treats the final record of a file without a trailing
/// newline as an incomplete fragment and parks it in `suffix_bytes`. On the BGZF
/// path the last block has no successor to stitch that fragment with, so before
/// the EOF flush a valid matched pair simply hung — completion required the
/// suffix buffers to be empty, and they never were. `RECORD_COUNT` spans several
/// BGZF blocks so both cross-block stitching and the terminal flush are covered.
/// gzip and plain (which never park a suffix) are included so the guarantee reads
/// as uniform across flavors.
#[rstest]
#[case::plain_t1(FastqFlavor::Plain, 1)]
#[case::plain_t4(FastqFlavor::Plain, 4)]
#[case::gzip_t1(FastqFlavor::Gzip, 1)]
#[case::gzip_t4(FastqFlavor::Gzip, 4)]
#[case::bgzf_t1(FastqFlavor::Bgzf, 1)]
#[case::bgzf_t4(FastqFlavor::Bgzf, 4)]
fn test_extract_matched_pair_without_trailing_newline(
    #[case] flavor: FastqFlavor,
    #[case] threads: usize,
) {
    const RECORD_COUNT: usize = 3_000;

    let tmp = TempDir::new().unwrap();
    let records = numbered_records(RECORD_COUNT);
    let bytes = fastq_bytes(&as_str_records(&records), false);
    let r1 = flavor.write_bytes(&tmp, "r1.fq", &bytes);
    let r2 = flavor.write_bytes(&tmp, "r2.fq", &bytes);
    let output = tmp.path().join("matched.bam");

    run_extract_pair(&r1, &r2, &output, threads).unwrap_or_else(|e| {
        panic!("matched pair without a trailing newline must extract cleanly ({flavor:?}, t{threads}): {e}")
    });

    // Identity: the final record (the one parked in `suffix_bytes`) must be
    // present, so the emitted name multiset equals the input name multiset.
    let mut emitted = bam_read_names(&output);
    emitted.sort();
    let mut expected: Vec<String> =
        records.iter().flat_map(|(name, _, _)| [name.clone(), name.clone()]).collect();
    expected.sort();
    assert_eq!(
        emitted, expected,
        "emitted read names must equal the input read names ({flavor:?}, t{threads})"
    );
    assert_eq!(
        unpaired_record_count(&output),
        0,
        "every emitted record must be flagged paired ({flavor:?}, t{threads})"
    );
}

/// A FASTQ pair whose two streams hold different numbers of records must be
/// rejected, not silently truncated or silently emitted as fragments (#773).
///
/// The two shapes reach the rejection through different code, so both are run.
/// In `within_batch` the streams first disagree inside a batch index they both
/// reached, which the old `align_stream_records` truncated to the shorter stream,
/// dropping the surplus outright. In `whole_batches` the short stream ends on an
/// exact batch boundary, so the first divergent batch index carries *only* the
/// long stream — which the old code emitted as single-end fragments (gzip and
/// plain) or parked in the block merger forever (BGZF).
///
/// The short length is a multiple of the pipeline's per-batch record count at
/// every thread count under test (200 at 1 thread, 800 at 4+), which is what puts
/// `whole_batches` on the missing-stream path instead of the truncation path.
///
/// Assertions are by record identity: when the run wrongly succeeds, the panic
/// names the specific records lost and the specific records emitted unpaired.
#[rstest]
#[case::plain_t1(FastqFlavor::Plain, 1)]
#[case::plain_t4(FastqFlavor::Plain, 4)]
#[case::gzip_t1(FastqFlavor::Gzip, 1)]
#[case::gzip_t4(FastqFlavor::Gzip, 4)]
#[case::bgzf_t1(FastqFlavor::Bgzf, 1)]
#[case::bgzf_t4(FastqFlavor::Bgzf, 4)]
fn test_extract_rejects_mismatched_fastq_pair(#[case] flavor: FastqFlavor, #[case] threads: usize) {
    /// `(shape, short_count, long_count)`. 800 is a whole number of batches at
    /// both 200 and 800 records per batch; 850 leaves a partial one.
    const SHAPES: [(&str, usize, usize); 2] =
        [("within_batch", 850, 950), ("whole_batches", 800, 2_600)];

    for (shape, short_count, long_count) in SHAPES {
        // Run both orientations: the surplus in R1 and the surplus in R2.
        for surplus_in_r1 in [true, false] {
            let tmp = TempDir::new().unwrap();
            let long = numbered_records(long_count);
            let short = numbered_records(short_count);
            let (r1_recs, r2_recs) = if surplus_in_r1 { (&long, &short) } else { (&short, &long) };
            let r1 = flavor.write(&tmp, "r1.fq", &as_str_records(r1_recs));
            let r2 = flavor.write(&tmp, "r2.fq", &as_str_records(r2_recs));
            let output = tmp.path().join("mismatched.bam");
            let label = format!("{flavor:?}, t{threads}, {shape}, surplus_in_r1={surplus_in_r1}");

            let Err(error) = run_extract_pair(&r1, &r2, &output, threads) else {
                panic!(
                    "extract accepted a mismatched FASTQ pair ({label}): {}",
                    describe_mismatched_output(&output, &long)
                );
            };

            // The message must identify which input ran out, so the operator
            // knows which file to re-fetch — not merely that something is wrong.
            let message = format!("{error:#}");
            let expected_direction =
                if surplus_in_r1 { "R2 ended before R1" } else { "R1 ended before R2" };
            assert!(
                message.contains("out of sync"),
                "rejection must say the sources are out of sync ({label}): {message}"
            );
            assert!(
                message.contains(expected_direction),
                "rejection must name the stream that ended first ({label}): {message}"
            );

            // `extract` streams its output, so by the time it rejects the pair a
            // BAM header and some record blocks may already be on disk — a
            // rejected run does NOT guarantee the output is absent. What it must
            // guarantee is that the leftover is never a *complete-looking* BAM: if
            // it parses at all, it holds fewer records than the long stream alone
            // would, so a downstream tool cannot mistake it for a full extraction
            // (the truncated-but-valid-output failure class of #773).
            if let Some(emitted) = try_bam_record_count(&output) {
                assert!(
                    emitted < 2 * long_count,
                    "a rejected run must not leave a complete-looking BAM ({label}): \
                     found {emitted} record(s), long stream has {long_count}"
                );
            }
        }
    }
}

/// A no-`--threads` run (the chain at a single worker) must reject a mismatched
/// pair with a message that names which stream ended first, matching the
/// multi-worker pipeline and the help text's promise (#773). This pins the
/// directional out-of-sync wording on the no-`--threads` invocation specifically
/// — the other rejection test always passes `--threads N`.
#[rstest]
#[case::plain(FastqFlavor::Plain)]
#[case::gzip(FastqFlavor::Gzip)]
#[case::bgzf(FastqFlavor::Bgzf)]
fn test_extract_single_threaded_rejects_mismatched_fastq_pair(#[case] flavor: FastqFlavor) {
    for surplus_in_r1 in [true, false] {
        let tmp = TempDir::new().unwrap();
        let long = numbered_records(950);
        let short = numbered_records(850);
        let (r1_recs, r2_recs) = if surplus_in_r1 { (&long, &short) } else { (&short, &long) };
        let r1 = flavor.write(&tmp, "r1.fq", &as_str_records(r1_recs));
        let r2 = flavor.write(&tmp, "r2.fq", &as_str_records(r2_recs));
        let output = tmp.path().join("mismatched.bam");
        let label = format!("{flavor:?}, single-threaded, surplus_in_r1={surplus_in_r1}");

        let Err(error) = run_extract_pair_single_threaded(&r1, &r2, &output) else {
            panic!("single-threaded extract accepted a mismatched FASTQ pair ({label})");
        };

        let message = format!("{error:#}");
        let expected_direction =
            if surplus_in_r1 { "R2 ended before R1" } else { "R1 ended before R2" };
        assert!(
            message.contains("out of sync"),
            "rejection must say the sources are out of sync ({label}): {message}"
        );
        assert!(
            message.contains(expected_direction),
            "rejection must name the stream that ended first ({label}): {message}"
        );
    }
}

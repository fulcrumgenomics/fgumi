//! Integration tests for the fastq command.
//!
//! These tests invoke the fastq command in-process via `Command::execute()`.

use clap::Parser;
use fgumi_lib::commands::command::Command;
use fgumi_lib::commands::fastq::Fastq;
use fgumi_lib::sam::SamTag;
use fgumi_raw_bam::{SamBuilder, flags};
use noodles::bam;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use std::fs;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
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

    // Run fastq command with output written directly to file
    let cmd = Fastq::try_parse_from([
        "fastq",
        "-i",
        input_bam.to_str().unwrap(),
        "-o",
        output_fq.to_str().unwrap(),
    ])
    .expect("failed to parse fastq args");
    cmd.execute("fgumi fastq").expect("fastq command failed");

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

    let cmd = Fastq::try_parse_from([
        "fastq",
        "-i",
        input_bam.to_str().unwrap(),
        "-o",
        output_fq.to_str().unwrap(),
    ])
    .expect("failed to parse fastq args");
    cmd.execute("fgumi fastq").expect("fastq command failed");

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

    let cmd = Fastq::try_parse_from([
        "fastq",
        "-i",
        input_bam.to_str().unwrap(),
        "-n",
        "-o",
        output_fq.to_str().unwrap(),
    ])
    .expect("failed to parse fastq args");
    cmd.execute("fgumi fastq").expect("fastq command failed");

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

    let cmd = Fastq::try_parse_from([
        "fastq",
        "-i",
        input_bam.to_str().unwrap(),
        "-o",
        output_fq.to_str().unwrap(),
    ])
    .expect("failed to parse fastq args");
    cmd.execute("fgumi fastq").expect("fastq command failed");

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
    let cmd = Fastq::try_parse_from([
        "fastq",
        "-i",
        input_bam.to_str().unwrap(),
        "-o",
        output_fq.to_str().unwrap(),
    ])
    .expect("failed to parse fastq args");
    cmd.execute("fgumi fastq").expect("fastq command failed");

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

    let cmd = Fastq::try_parse_from([
        "fastq",
        "-i",
        input_bam.to_str().unwrap(),
        "-@",
        "4",
        "-o",
        output_fq.to_str().unwrap(),
    ])
    .expect("failed to parse fastq args");
    cmd.execute("fgumi fastq").expect("fastq command failed");

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
    let cmd = Fastq::try_parse_from([
        "fastq",
        "-i",
        input_bam.to_str().unwrap(),
        "-F",
        "0x900", // hex notation
        "-o",
        output_fq.to_str().unwrap(),
    ])
    .expect("failed to parse fastq args");
    cmd.execute("fgumi fastq").expect("fastq command failed");

    let records = parse_fastq_records(&output_fq);
    assert_eq!(records.len(), 2);
}

/// `--output` pointing at the same file as `--input` must fail before any
/// truncation, since `File::create` would otherwise clobber the BAM data.
#[test]
fn test_fastq_output_same_as_input_rejected() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");

    create_paired_bam(&input_bam, vec![("read1", "ACGT", "IIII", "TGCA", "IIII", false)]);
    let input_size_before = std::fs::metadata(&input_bam).expect("stat input").len();

    let cmd = Fastq::try_parse_from([
        "fastq",
        "-i",
        input_bam.to_str().unwrap(),
        "-o",
        input_bam.to_str().unwrap(),
    ])
    .expect("failed to parse fastq args");
    let err =
        cmd.execute("fgumi fastq").expect_err("execute must reject identical --input/--output");
    assert!(
        err.to_string().contains("is the same file as --input"),
        "unexpected error message: {err}"
    );

    // Most importantly: the input BAM must not have been truncated.
    let input_size_after = std::fs::metadata(&input_bam).expect("stat input").len();
    assert_eq!(
        input_size_before, input_size_after,
        "input BAM was truncated/clobbered by --output=--input"
    );
}

/// `--output` pointing at a symlink that resolves to `--input` must also be
/// rejected. Lexical `PathBuf` comparison misses this case, so the validator
/// also canonicalises both sides when the output already exists.
#[cfg(unix)]
#[test]
fn test_fastq_output_symlink_to_input_rejected() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_link = temp_dir.path().join("output.bam");

    create_paired_bam(&input_bam, vec![("read1", "ACGT", "IIII", "TGCA", "IIII", false)]);
    std::os::unix::fs::symlink(&input_bam, &output_link).expect("create symlink");
    let input_size_before = std::fs::metadata(&input_bam).expect("stat input").len();

    let cmd = Fastq::try_parse_from([
        "fastq",
        "-i",
        input_bam.to_str().unwrap(),
        "-o",
        output_link.to_str().unwrap(),
    ])
    .expect("failed to parse fastq args");
    let err =
        cmd.execute("fgumi fastq").expect_err("execute must reject --output symlinked to --input");
    assert!(
        err.to_string().contains("is the same file as --input"),
        "unexpected error message: {err}"
    );

    // The input BAM must not have been truncated through the symlink.
    let input_size_after = std::fs::metadata(&input_bam).expect("stat input").len();
    assert_eq!(
        input_size_before, input_size_after,
        "input BAM was truncated/clobbered through symlinked --output"
    );
}

/// Create a single-end BAM whose reads carry an `RX` UMI tag.
fn create_bam_with_umis(path: &PathBuf, reads: &[(&str, &str)]) {
    let header = create_minimal_header("chr1", 10000);
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");

    for (name, umi) in reads {
        let record = {
            let mut b = SamBuilder::new();
            b.read_name(name.as_bytes())
                .sequence(b"ACGTACGT")
                .qualities(&[30u8; 8])
                .flags(0)
                .ref_id(0)
                .pos(99)
                .mapq(60);
            b.add_string_tag(fgumi_lib::sam::SamTag::RX, umi.as_bytes());
            b.build()
        };
        writer
            .write_alignment_record(&header, &to_record_buf(&record))
            .expect("Failed to write record");
    }
    writer.try_finish().expect("Failed to finish BAM");
}

/// `-a` appends the record's UMI to the read name, rewriting fgumi's stored `-`
/// duplex separator to `+` — reproducing `samtools fastq -U` / the DRAGEN layout.
#[test]
fn test_fastq_annotates_read_names_with_umi() {
    let dir = TempDir::new().expect("temp dir");
    let input_bam = dir.path().join("umi.bam");
    let output_fq = dir.path().join("out.fq");
    create_bam_with_umis(&input_bam, &[("readA", "ACGT"), ("readB", "ACGT-TTTT")]);

    let cmd = Fastq::try_parse_from([
        "fastq",
        "-i",
        input_bam.to_str().unwrap(),
        "-o",
        output_fq.to_str().unwrap(),
        "-a",
        "true",
    ])
    .expect("failed to parse fastq args");
    cmd.execute("fgumi fastq").expect("fastq command failed");

    let records = parse_fastq_records(&output_fq);
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].0, "readA:ACGT", "simplex UMI appended after ':'");
    assert_eq!(records[1].0, "readB:ACGT+TTTT", "duplex '-' rewritten to '+'");

    // Without -a the names are untouched.
    let plain_fq = dir.path().join("plain.fq");
    let cmd = Fastq::try_parse_from([
        "fastq",
        "-i",
        input_bam.to_str().unwrap(),
        "-o",
        plain_fq.to_str().unwrap(),
    ])
    .expect("failed to parse fastq args");
    cmd.execute("fgumi fastq").expect("fastq command failed");
    assert_eq!(parse_fastq_records(&plain_fq)[0].0, "readA");
}

/// A `.gz` output path must be real BGZF, not plain text under a compressed name.
#[test]
fn test_fastq_gz_output_is_real_bgzf() {
    let dir = TempDir::new().expect("temp dir");
    let input_bam = dir.path().join("umi.bam");
    let output_gz = dir.path().join("out.fq.gz");
    let output_fq = dir.path().join("out.fq");
    create_bam_with_umis(&input_bam, &[("readA", "ACGT"), ("readB", "TTTT")]);

    for output in [&output_gz, &output_fq] {
        let cmd = Fastq::try_parse_from([
            "fastq",
            "-i",
            input_bam.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .expect("failed to parse fastq args");
        cmd.execute("fgumi fastq").expect("fastq command failed");
    }

    let bytes = fs::read(&output_gz).expect("read gz output");
    assert_eq!(&bytes[..2], &[0x1f, 0x8b], "output must carry gzip magic, not plain text");
    assert!(bytes.ends_with(&fgumi_bgzf::BGZF_EOF), "BGZF stream must end with the EOF marker");

    // Decompresses byte-for-byte to the FASTQ the uncompressed path produces.
    let mut cursor = std::io::Cursor::new(&bytes);
    let blocks = fgumi_bgzf::read_raw_blocks(&mut cursor, 64).expect("read blocks");
    let mut decompressor = libdeflater::Decompressor::new();
    let mut decoded = Vec::new();
    for block in &blocks {
        decoded.extend_from_slice(
            &fgumi_bgzf::decompress_block(block, &mut decompressor).expect("decompress"),
        );
    }
    let text = String::from_utf8(decoded).expect("utf8");
    let expected = fs::read_to_string(&output_fq).expect("read plain fq output");
    // Guard the oracle itself, so a shared defect in both paths cannot pass vacuously.
    assert_eq!(expected.lines().count(), 8, "plain-path oracle: two records, four lines each");
    assert!(text.starts_with("@readA\n"), "decompressed FASTQ should start with readA");
    assert_eq!(text, expected, "BGZF output must decompress to the plain-path FASTQ byte-for-byte");
}

// ============================================================================
// Paired split output (`-1`/`-2`/`-0`) — the BAM→FASTQ paired path on the chain.
// ============================================================================

/// Build a BAM from explicit records `(name, flags, seq, qual_ascii, rx)`,
/// written in the given order. `rx` adds an `RX` tag (UMI) when `Some`.
fn create_bam_with_records(path: &PathBuf, records: &[(&str, u16, &str, &str, Option<&str>)]) {
    let header = create_minimal_header("chr1", 10000);
    let mut writer = bam::io::Writer::new(fs::File::create(path).expect("create BAM"));
    writer.write_header(&header).expect("write header");
    for (name, flag_bits, seq, qual, rx) in records {
        let q: Vec<u8> = qual.bytes().map(|b| b - 33).collect();
        let mut b = SamBuilder::new();
        b.read_name(name.as_bytes()).sequence(seq.as_bytes()).qualities(&q).flags(*flag_bits);
        if let Some(rx) = rx {
            b.add_string_tag(SamTag::RX, rx.as_bytes());
        }
        let rec = b.build();
        writer.write_alignment_record(&header, &to_record_buf(&rec)).expect("write record");
    }
    writer.try_finish().expect("finish BAM");
}

/// Names in the order they appear (record names, suffixes included if any).
fn names(records: &[(String, String, String)]) -> Vec<String> {
    records.iter().map(|(n, _, _)| n.clone()).collect()
}

/// Paired split routes R1 → `--out1`, R2 → `--out2`, and single-end/ambiguous
/// reads → `--out0`, mirroring `samtools fastq -1/-2/-0`, and omits the `/1` `/2`
/// suffix so R1 and R2 read names match.
#[test]
fn test_fastq_paired_split_routes_and_omits_suffix() {
    let dir = TempDir::new().expect("temp dir");
    let input = dir.path().join("in.bam");
    let (r1, r2, r0) =
        (dir.path().join("r1.fq"), dir.path().join("r2.fq"), dir.path().join("r0.fq"));
    create_bam_with_records(
        &input,
        &[
            ("p1", flags::PAIRED | flags::FIRST_SEGMENT, "ACGTACGT", "IIIIIIII", None),
            ("p1", flags::PAIRED | flags::LAST_SEGMENT, "TTTTGGGG", "JJJJJJJJ", None),
            ("p2", flags::PAIRED | flags::FIRST_SEGMENT, "AAAACCCC", "KKKKKKKK", None),
            ("p2", flags::PAIRED | flags::LAST_SEGMENT, "GGGGAAAA", "LLLLLLLL", None),
            // Single-end (neither FIRST nor LAST) → "other".
            ("solo", 0, "CCCCTTTT", "MMMMMMMM", None),
            // Ambiguous (BOTH segment bits set) → "other" via a separate
            // routing branch than the single-end case above.
            (
                "ambig",
                flags::PAIRED | flags::FIRST_SEGMENT | flags::LAST_SEGMENT,
                "GATTACAG",
                "NNNNNNNN",
                None,
            ),
        ],
    );

    Fastq::try_parse_from([
        "fastq",
        "-i",
        input.to_str().unwrap(),
        "-1",
        r1.to_str().unwrap(),
        "-2",
        r2.to_str().unwrap(),
        "-0",
        r0.to_str().unwrap(),
    ])
    .expect("parse")
    .execute("fgumi fastq")
    .expect("paired fastq");

    let (rec1, rec2, rec0) =
        (parse_fastq_records(&r1), parse_fastq_records(&r2), parse_fastq_records(&r0));
    // Assert full record identity (name, sequence, quality), not just routing:
    // a transcode that corrupted a quality string or a later record's sequence
    // would pass a names-only or `[0]`-only check. Split mode omits the `/1` `/2`
    // suffix, so R1 and R2 names match.
    assert_eq!(
        rec1,
        vec![
            ("p1".into(), "ACGTACGT".into(), "IIIIIIII".into()),
            ("p2".into(), "AAAACCCC".into(), "KKKKKKKK".into()),
        ],
        "R1 file: FIRST_SEGMENT reads, in order, no /1 suffix"
    );
    assert_eq!(
        rec2,
        vec![
            ("p1".into(), "TTTTGGGG".into(), "JJJJJJJJ".into()),
            ("p2".into(), "GGGGAAAA".into(), "LLLLLLLL".into()),
        ],
        "R2 file: LAST_SEGMENT reads, in order, no /2 suffix"
    );
    assert_eq!(
        rec0,
        vec![
            ("solo".into(), "CCCCTTTT".into(), "MMMMMMMM".into()),
            ("ambig".into(), "GATTACAG".into(), "NNNNNNNN".into()),
        ],
        "out0 file: the single-end read and the both-segment-bits ambiguous read"
    );
}

/// Degenerate record: an empty `SEQ` with zero-length quality still fans out to
/// the correct branch and round-trips as the four-line record `@name\n\n+\n\n`.
/// Here `extract_sequence_into` and `encode_quality_into` both produce empty
/// buffers — a path no other paired-split case exercises (every other builder
/// caller passes a non-empty seq/qual).
#[test]
fn test_fastq_paired_split_emits_degenerate_empty_record() {
    let dir = TempDir::new().expect("temp dir");
    let input = dir.path().join("in.bam");
    let (r1, r2, r0) =
        (dir.path().join("r1.fq"), dir.path().join("r2.fq"), dir.path().join("r0.fq"));
    create_bam_with_records(
        &input,
        &[("empty", flags::PAIRED | flags::FIRST_SEGMENT, "", "", None)],
    );

    Fastq::try_parse_from([
        "fastq",
        "-i",
        input.to_str().unwrap(),
        "-1",
        r1.to_str().unwrap(),
        "-2",
        r2.to_str().unwrap(),
        "-0",
        r0.to_str().unwrap(),
    ])
    .expect("parse")
    .execute("fgumi fastq")
    .expect("degenerate paired fastq");

    // FIRST_SEGMENT with empty seq/qual → R1 as `@empty\n\n+\n\n`.
    assert_eq!(
        fs::read_to_string(&r1).expect("read r1"),
        "@empty\n\n+\n\n",
        "R1 must carry the empty-seq record round-tripped as @empty\\n\\n+\\n\\n"
    );
    assert!(fs::read_to_string(&r2).expect("read r2").is_empty(), "R2 file must be byte-empty");
    assert!(fs::read_to_string(&r0).expect("read r0").is_empty(), "out0 file must be byte-empty");
}

/// An input with only R1 records must still complete (not hang): the encode
/// step emits an empty block on the R2 and "other" branches for every batch, so
/// their reorder stages are never starved. The R2/other files come out empty.
#[test]
fn test_fastq_paired_split_only_r1_completes_with_empty_r2() {
    let dir = TempDir::new().expect("temp dir");
    let input = dir.path().join("in.bam");
    let (r1, r2, r0) =
        (dir.path().join("r1.fq"), dir.path().join("r2.fq"), dir.path().join("r0.fq"));
    create_bam_with_records(
        &input,
        &[
            ("a", flags::PAIRED | flags::FIRST_SEGMENT, "ACGT", "IIII", None),
            ("b", flags::PAIRED | flags::FIRST_SEGMENT, "TTTT", "JJJJ", None),
        ],
    );

    Fastq::try_parse_from([
        "fastq",
        "-i",
        input.to_str().unwrap(),
        "--threads",
        "4",
        "-1",
        r1.to_str().unwrap(),
        "-2",
        r2.to_str().unwrap(),
        "-0",
        r0.to_str().unwrap(),
    ])
    .expect("parse")
    .execute("fgumi fastq")
    .expect("paired fastq with only R1s");

    assert_eq!(names(&parse_fastq_records(&r1)), vec!["a", "b"], "all R1s present");
    assert!(fs::read_to_string(&r2).expect("read r2").is_empty(), "R2 file is empty");
    assert!(fs::read_to_string(&r0).expect("read r0").is_empty(), "out0 file is empty");
}

/// Symmetric to the R1-only case (sibling parity): an input with only R2
/// (`LAST_SEGMENT`) records must still complete (not hang) — the encode step emits
/// an empty block on the R1 and "other" branches for every batch, so their
/// reorder stages are never starved. The R1/other files come out empty.
#[test]
fn test_fastq_paired_split_only_r2_completes_with_empty_r1() {
    let dir = TempDir::new().expect("temp dir");
    let input = dir.path().join("in.bam");
    let (r1, r2, r0) =
        (dir.path().join("r1.fq"), dir.path().join("r2.fq"), dir.path().join("r0.fq"));
    create_bam_with_records(
        &input,
        &[
            ("a", flags::PAIRED | flags::LAST_SEGMENT, "ACGT", "IIII", None),
            ("b", flags::PAIRED | flags::LAST_SEGMENT, "TTTT", "JJJJ", None),
        ],
    );

    Fastq::try_parse_from([
        "fastq",
        "-i",
        input.to_str().unwrap(),
        "--threads",
        "4",
        "-1",
        r1.to_str().unwrap(),
        "-2",
        r2.to_str().unwrap(),
        "-0",
        r0.to_str().unwrap(),
    ])
    .expect("parse")
    .execute("fgumi fastq")
    .expect("paired fastq with only R2s");

    assert_eq!(names(&parse_fastq_records(&r2)), vec!["a", "b"], "all R2s present");
    assert!(fs::read_to_string(&r1).expect("read r1").is_empty(), "R1 file is empty");
    assert!(fs::read_to_string(&r0).expect("read r0").is_empty(), "out0 file is empty");
}

/// A `.gz` paired output is real BGZF and decompresses byte-for-byte to the
/// plain-path FASTQ.
#[test]
fn test_fastq_paired_split_gz_is_valid_bgzf() {
    let dir = TempDir::new().expect("temp dir");
    let input = dir.path().join("in.bam");
    create_bam_with_records(
        &input,
        &[
            ("p1", flags::PAIRED | flags::FIRST_SEGMENT, "ACGTACGT", "IIIIIIII", None),
            ("p1", flags::PAIRED | flags::LAST_SEGMENT, "TTTTGGGG", "JJJJJJJJ", None),
        ],
    );
    let plain_r1 = dir.path().join("plain_r1.fq");
    let plain_r2 = dir.path().join("plain_r2.fq");
    Fastq::try_parse_from([
        "fastq",
        "-i",
        input.to_str().unwrap(),
        "-1",
        plain_r1.to_str().unwrap(),
        "-2",
        plain_r2.to_str().unwrap(),
    ])
    .expect("parse")
    .execute("fgumi fastq")
    .expect("plain paired");

    let gz_r1 = dir.path().join("r1.fq.gz");
    let gz_r2 = dir.path().join("r2.fq.gz");
    Fastq::try_parse_from([
        "fastq",
        "-i",
        input.to_str().unwrap(),
        "-1",
        gz_r1.to_str().unwrap(),
        "-2",
        gz_r2.to_str().unwrap(),
    ])
    .expect("parse")
    .execute("fgumi fastq")
    .expect("gz paired");

    let bytes = fs::read(&gz_r1).expect("read gz r1");
    assert_eq!(&bytes[..2], &[0x1f, 0x8b], "R1.gz must carry gzip magic");
    assert!(bytes.ends_with(&fgumi_bgzf::BGZF_EOF), "R1.gz must end with the BGZF EOF marker");

    // Sibling-sink parity: the R2 `.gz` sink must be BGZF too, else a per-sink
    // compression bug (plain text under a `.gz` name) would slip past the
    // R1-only checks above.
    let r2_bytes = fs::read(&gz_r2).expect("read gz r2");
    assert_eq!(&r2_bytes[..2], &[0x1f, 0x8b], "R2.gz must carry gzip magic");
    assert!(r2_bytes.ends_with(&fgumi_bgzf::BGZF_EOF), "R2.gz must end with the BGZF EOF marker");

    // Sibling-sink parity: format checks alone pass on an empty or wrong R2
    // stream, so decompress R2 and pin it to the exact R2 record too.
    let mut r2_cursor = std::io::Cursor::new(&r2_bytes);
    let r2_blocks = fgumi_bgzf::read_raw_blocks(&mut r2_cursor, 64).expect("read r2 blocks");
    let mut r2_decompressor = libdeflater::Decompressor::new();
    let mut r2_decoded = Vec::new();
    for block in &r2_blocks {
        r2_decoded.extend_from_slice(
            &fgumi_bgzf::decompress_block(block, &mut r2_decompressor).expect("decompress r2"),
        );
    }
    let expected_r2 = fs::read(&plain_r2).expect("read plain r2");
    assert!(!expected_r2.is_empty(), "plain-path R2 oracle must be non-empty");
    assert_eq!(
        expected_r2,
        b"@p1\nTTTTGGGG\n+\nJJJJJJJJ\n".to_vec(),
        "oracle: the single R2 record"
    );
    assert_eq!(r2_decoded, expected_r2, "R2.gz must decompress to the plain-path R2 byte-for-byte");

    let mut cursor = std::io::Cursor::new(&bytes);
    let blocks = fgumi_bgzf::read_raw_blocks(&mut cursor, 64).expect("read blocks");
    let mut decompressor = libdeflater::Decompressor::new();
    let mut decoded = Vec::new();
    for block in &blocks {
        decoded.extend_from_slice(
            &fgumi_bgzf::decompress_block(block, &mut decompressor).expect("decompress"),
        );
    }
    let expected = fs::read(&plain_r1).expect("read plain r1");
    // Guard the oracle itself: if a shared defect made both paths emit nothing,
    // an empty-vs-empty comparison would pass vacuously. Pin the actual R1 record.
    assert!(!expected.is_empty(), "plain-path oracle must be non-empty");
    assert_eq!(expected, b"@p1\nACGTACGT\n+\nIIIIIIII\n".to_vec(), "oracle: the single R1 record");
    assert_eq!(decoded, expected, "R1.gz must decompress to the plain-path R1 byte-for-byte");
}

/// Assert `path` is a real BGZF stream (gzip magic + BGZF EOF marker) whose
/// decompressed bytes equal `expected_plain` exactly. Shared by the mixed-case
/// suffix test so each sink is pinned to content, not just format.
fn assert_bgzf_decompresses_to(path: &Path, expected_plain: &[u8], label: &str) {
    let bytes = fs::read(path).unwrap_or_else(|e| panic!("read {label}: {e}"));
    assert_eq!(&bytes[..2], &[0x1f, 0x8b], "{label} must carry gzip magic");
    assert!(bytes.ends_with(&fgumi_bgzf::BGZF_EOF), "{label} must end with the BGZF EOF marker");
    let mut cursor = std::io::Cursor::new(&bytes);
    let blocks = fgumi_bgzf::read_raw_blocks(&mut cursor, 64).expect("read blocks");
    let mut decompressor = libdeflater::Decompressor::new();
    let mut decoded = Vec::new();
    for block in &blocks {
        decoded.extend_from_slice(
            &fgumi_bgzf::decompress_block(block, &mut decompressor).expect("decompress"),
        );
    }
    assert!(!expected_plain.is_empty(), "{label} oracle must be non-empty");
    assert_eq!(decoded, expected_plain, "{label} must decompress to the plain-path bytes");
}

/// Compressed-suffix detection is case-insensitive and spans the `bgzf`/`bgz`
/// arms, not just lowercase `.gz`: an R1 `.bgzf` sink and an uppercase R2 `.BGZ`
/// sink must both emit BGZF (not plain text under a compressed name). Mirrors
/// `test_fastq_paired_split_gz_is_valid_bgzf` for the other suffix-table arms.
#[test]
fn test_fastq_paired_split_mixed_case_bgzf_suffixes() {
    let dir = TempDir::new().expect("temp dir");
    let input = dir.path().join("in.bam");
    create_bam_with_records(
        &input,
        &[
            ("p1", flags::PAIRED | flags::FIRST_SEGMENT, "ACGTACGT", "IIIIIIII", None),
            ("p1", flags::PAIRED | flags::LAST_SEGMENT, "TTTTGGGG", "JJJJJJJJ", None),
        ],
    );

    // Plain-path oracles for byte-for-byte content comparison.
    let plain_r1 = dir.path().join("plain_r1.fq");
    let plain_r2 = dir.path().join("plain_r2.fq");
    Fastq::try_parse_from([
        "fastq",
        "-i",
        input.to_str().unwrap(),
        "-1",
        plain_r1.to_str().unwrap(),
        "-2",
        plain_r2.to_str().unwrap(),
    ])
    .expect("parse")
    .execute("fgumi fastq")
    .expect("plain paired");
    let expected_r1 = fs::read(&plain_r1).expect("read plain r1");
    let expected_r2 = fs::read(&plain_r2).expect("read plain r2");

    // R1 exercises the lowercase `.bgzf` arm; R2 the uppercase `.BGZ` arm.
    let bgzf_r1 = dir.path().join("r1.bgzf");
    let bgz_r2 = dir.path().join("r2.BGZ");
    Fastq::try_parse_from([
        "fastq",
        "-i",
        input.to_str().unwrap(),
        "-1",
        bgzf_r1.to_str().unwrap(),
        "-2",
        bgz_r2.to_str().unwrap(),
    ])
    .expect("parse")
    .execute("fgumi fastq")
    .expect("mixed-case bgzf paired");

    assert_bgzf_decompresses_to(&bgzf_r1, &expected_r1, "R1.bgzf");
    assert_bgzf_decompresses_to(&bgz_r2, &expected_r2, "R2.BGZ");
}

/// `--annotate-read-names` appends the UMI in paired-split mode too: a duplex
/// UMI `AAAA-CCCC` becomes `name:AAAA+CCCC`, a simplex UMI `name:GGGG`.
#[test]
fn test_fastq_paired_split_annotates_umi() {
    let dir = TempDir::new().expect("temp dir");
    let input = dir.path().join("in.bam");
    let (r1, r2) = (dir.path().join("r1.fq"), dir.path().join("r2.fq"));
    create_bam_with_records(
        &input,
        &[
            ("dup", flags::PAIRED | flags::FIRST_SEGMENT, "ACGT", "IIII", Some("AAAA-CCCC")),
            ("dup", flags::PAIRED | flags::LAST_SEGMENT, "TTTT", "JJJJ", Some("AAAA-CCCC")),
            ("smp", flags::PAIRED | flags::FIRST_SEGMENT, "GGGG", "KKKK", Some("GGGG")),
            ("smp", flags::PAIRED | flags::LAST_SEGMENT, "CCCC", "LLLL", Some("GGGG")),
        ],
    );

    Fastq::try_parse_from([
        "fastq",
        "-i",
        input.to_str().unwrap(),
        "--annotate-read-names",
        "-1",
        r1.to_str().unwrap(),
        "-2",
        r2.to_str().unwrap(),
    ])
    .expect("parse")
    .execute("fgumi fastq")
    .expect("annotated paired fastq");

    // Duplex `-` rewritten to `+`; simplex passed through. R1 and R2 names match.
    assert_eq!(names(&parse_fastq_records(&r1)), vec!["dup:AAAA+CCCC", "smp:GGGG"]);
    assert_eq!(names(&parse_fastq_records(&r2)), vec!["dup:AAAA+CCCC", "smp:GGGG"]);
}

/// `--out1` without `--out2` is a clap error (they are required together).
/// Assert the specific `MissingRequiredArgument` kind, not just `is_err()`: a
/// parse failure for any unrelated reason (e.g. `-1`/`-2` renamed or removed —
/// the "flag silently does nothing" failure this test exists to catch) would
/// satisfy a bare `is_err()`.
#[test]
fn test_fastq_paired_requires_out2() {
    let err = Fastq::try_parse_from(["fastq", "-i", "in.bam", "-1", "r1.fq"])
        .expect_err("--out1 without --out2 must be rejected by clap");
    assert_eq!(
        err.kind(),
        clap::error::ErrorKind::MissingRequiredArgument,
        "--out1 without --out2 must fail as a missing-required-argument error, got: {err}"
    );
}

/// `--out1`/`--out2` conflict with the interleaved `--output`. Assert the
/// specific `ArgumentConflict` kind (see `test_fastq_paired_requires_out2` for
/// why `is_err()` alone is too weak).
#[test]
fn test_fastq_paired_conflicts_with_output() {
    let err = Fastq::try_parse_from([
        "fastq", "-i", "in.bam", "-o", "out.fq", "-1", "r1.fq", "-2", "r2.fq",
    ])
    .expect_err("--output with --out1/--out2 must be rejected by clap");
    assert_eq!(
        err.kind(),
        clap::error::ErrorKind::ArgumentConflict,
        "--output with --out1/--out2 must fail as an argument-conflict error, got: {err}"
    );
}

/// Paired `--out1 -` (R1 → stdout) with `--out2` a file and `--out0` omitted must
/// be rejected: with no `--out0`, "other" reads default to stdout too, so R1 and
/// the "other" stream would interleave on one stdout. The implicit stdout target
/// must reach the collision guard even though it is not a user-specified path.
#[test]
fn test_fastq_paired_out1_stdout_collides_with_default_other_stdout() {
    let dir = TempDir::new().expect("temp dir");
    let input = dir.path().join("in.bam");
    create_bam_with_records(
        &input,
        &[("p1", flags::PAIRED | flags::FIRST_SEGMENT, "ACGT", "IIII", None)],
    );

    let r2 = dir.path().join("r2.fq");
    let err = Fastq::try_parse_from([
        "fastq",
        "-i",
        input.to_str().unwrap(),
        "-1",
        "-", // R1 → stdout
        "-2",
        r2.to_str().unwrap(),
        // no --out0: "other" also defaults to stdout, colliding with -1 -
    ])
    .expect("parse")
    .execute("fgumi fastq")
    .expect_err("--out1 - with a default-stdout --out0 must be rejected");
    assert!(err.to_string().contains("stdout"), "got: {err}");
    assert!(!r2.exists(), "collision rejection must occur before opening the R2 sink");
}

/// Two paired outputs naming the same file are rejected before any write —
/// including when the two spellings differ only by a leading `./` and the file
/// does not yet exist, the alias case a naive lexical or exists()-gated check
/// misses (two sink threads would truncate one physical file).
#[test]
fn test_fastq_paired_duplicate_output_rejected() {
    let dir = TempDir::new().expect("temp dir");
    let input = dir.path().join("in.bam");
    create_bam_with_records(
        &input,
        &[("p1", flags::PAIRED | flags::FIRST_SEGMENT, "ACGT", "IIII", None)],
    );

    // Identical spellings.
    let same = dir.path().join("same.fq");
    let err = Fastq::try_parse_from([
        "fastq",
        "-i",
        input.to_str().unwrap(),
        "-1",
        same.to_str().unwrap(),
        "-2",
        same.to_str().unwrap(),
    ])
    .expect("parse")
    .execute("fgumi fastq")
    .expect_err("identical -1/-2 paths must be rejected");
    assert!(err.to_string().contains("write to"), "got: {err}");
    assert!(!same.exists(), "rejection must happen before any sink is opened");

    // Aliased spellings of a not-yet-existing file: `<dir>/out.fq` vs
    // `<dir>/./out.fq`. Both are absolute (no cwd mutation, so this is safe under
    // thread-parallel `cargo test`), neither exists yet, and both canonicalize to
    // one path — the alias case an exists()-gated or lexical check would miss.
    let plain = dir.path().join("out.fq");
    let dotted = dir.path().join(".").join("out.fq");
    let err = Fastq::try_parse_from([
        "fastq",
        "-i",
        input.to_str().unwrap(),
        "-1",
        plain.to_str().unwrap(),
        "-2",
        dotted.to_str().unwrap(),
    ])
    .expect("parse")
    .execute("fgumi fastq")
    .expect_err("`out.fq` and `./out.fq` are the same file and must be rejected");
    assert!(err.to_string().contains("write to"), "got: {err}");
    assert!(!plain.exists(), "rejection must happen before any sink is opened");

    // `--out0` participates in the same collision guard, but the cases above
    // only pair `-1` with `-2`. A guard built from the R1/R2 pair alone would
    // let `-1 p -0 p` open one physical file from two sink threads. Pair
    // `--out0` with `--out1` to cover the third slot.
    let shared = dir.path().join("shared.fq");
    let distinct = dir.path().join("distinct.fq");
    let err = Fastq::try_parse_from([
        "fastq",
        "-i",
        input.to_str().unwrap(),
        "-1",
        shared.to_str().unwrap(),
        "-2",
        distinct.to_str().unwrap(),
        "-0",
        shared.to_str().unwrap(),
    ])
    .expect("parse")
    .execute("fgumi fastq")
    .expect_err("--out0 colliding with --out1 must be rejected");
    assert!(err.to_string().contains("write to"), "got: {err}");
    assert!(!shared.exists(), "rejection must happen before any sink is opened");
}

/// Run `fgumi fastq` in paired-split mode with the given `--out1`/`--out2`/`--out0`
/// paths and assert that it is rejected because one of them is `--input`, before
/// any sink is opened — the input BAM must be left byte-for-byte intact. Shared by
/// the direct-path and symlink alias tests below.
fn assert_paired_output_alias_rejected(input: &Path, out1: &Path, out2: &Path, out0: &Path) {
    let input_size_before = std::fs::metadata(input).expect("stat input").len();
    let err = Fastq::try_parse_from([
        "fastq",
        "-i",
        input.to_str().unwrap(),
        "-1",
        out1.to_str().unwrap(),
        "-2",
        out2.to_str().unwrap(),
        "-0",
        out0.to_str().unwrap(),
    ])
    .expect("parse")
    .execute("fgumi fastq")
    .expect_err("a paired output aliasing --input must be rejected");
    assert!(err.to_string().contains("is the same file as --input"), "got: {err}");
    // The guard runs before any sink is opened, so the input BAM is untouched.
    let input_size_after = std::fs::metadata(input).expect("stat input").len();
    assert_eq!(
        input_size_before, input_size_after,
        "input BAM was truncated/clobbered by a paired output aliasing --input"
    );
}

/// Each paired output (`--out1`, `--out2`, `--out0`) pointing directly at `--input`
/// must be rejected before any write — the paired-path analogue of the interleaved
/// `test_fastq_output_same_as_input_rejected`. All four flags route through the same
/// `reject_write_aliasing_input` over the `(path, flag)` output slice, but only the
/// single `--output` case was covered until now.
#[test]
fn test_fastq_paired_output_same_as_input_rejected() {
    let dir = TempDir::new().expect("temp dir");
    let input = dir.path().join("in.bam");
    create_bam_with_records(
        &input,
        &[("p1", flags::PAIRED | flags::FIRST_SEGMENT, "ACGT", "IIII", None)],
    );
    let r1 = dir.path().join("r1.fq");
    let r2 = dir.path().join("r2.fq");
    let other = dir.path().join("other.fq");

    // Exactly one output aliases --input in each case; the siblings are distinct
    // valid paths so clap's `requires` chain is satisfied and the only collision
    // under test is <aliased> == --input.
    assert_paired_output_alias_rejected(&input, &input, &r2, &other); // --out1
    assert_paired_output_alias_rejected(&input, &r1, &input, &other); // --out2
    assert_paired_output_alias_rejected(&input, &r1, &r2, &input); // --out0
}

/// The same guard must also catch a paired output reached through a symlink that
/// resolves to `--input` — a lexical `PathBuf` compare misses this, so the
/// validator canonicalises both sides. Paired-path analogue of the interleaved
/// `test_fastq_output_symlink_to_input_rejected`.
#[cfg(unix)]
#[test]
fn test_fastq_paired_output_symlink_to_input_rejected() {
    let dir = TempDir::new().expect("temp dir");
    let input = dir.path().join("in.bam");
    create_bam_with_records(
        &input,
        &[("p1", flags::PAIRED | flags::FIRST_SEGMENT, "ACGT", "IIII", None)],
    );
    let link = dir.path().join("alias.link");
    std::os::unix::fs::symlink(&input, &link).expect("create symlink");
    let r1 = dir.path().join("r1.fq");
    let r2 = dir.path().join("r2.fq");
    let other = dir.path().join("other.fq");

    // The symlink resolves to --input, so it must be rejected in each output slot.
    assert_paired_output_alias_rejected(&input, &link, &r2, &other); // --out1
    assert_paired_output_alias_rejected(&input, &r1, &link, &other); // --out2
    assert_paired_output_alias_rejected(&input, &r1, &r2, &link); // --out0
}

/// Paired split honors `-F`/`-f` flag filters the same way the interleaved path
/// does (sibling-path parity): a read excluded by `--exclude-flags` never
/// reaches any of `--out1`/`--out2`/`--out0`, and `--require-flags` keeps only
/// matching reads.
#[test]
fn test_fastq_paired_split_applies_flag_filters() {
    let dir = TempDir::new().expect("temp dir");
    let input = dir.path().join("in.bam");
    let (r1, r2, r0) =
        (dir.path().join("r1.fq"), dir.path().join("r2.fq"), dir.path().join("r0.fq"));
    // p1 is a normal pair; dup is a duplicate-flagged (0x400) pair that
    // --exclude-flags must drop from every output.
    create_bam_with_records(
        &input,
        &[
            ("p1", flags::PAIRED | flags::FIRST_SEGMENT, "ACGTACGT", "IIIIIIII", None),
            ("p1", flags::PAIRED | flags::LAST_SEGMENT, "TTTTGGGG", "JJJJJJJJ", None),
            ("dup", flags::PAIRED | flags::FIRST_SEGMENT | flags::DUPLICATE, "AAAA", "IIII", None),
            ("dup", flags::PAIRED | flags::LAST_SEGMENT | flags::DUPLICATE, "CCCC", "JJJJ", None),
        ],
    );

    Fastq::try_parse_from([
        "fastq",
        "-i",
        input.to_str().unwrap(),
        "-F",
        "0x400", // exclude duplicates
        "-1",
        r1.to_str().unwrap(),
        "-2",
        r2.to_str().unwrap(),
        "-0",
        r0.to_str().unwrap(),
    ])
    .expect("parse")
    .execute("fgumi fastq")
    .expect("filtered paired fastq");

    // Only p1 survives; the duplicate-flagged pair is dropped from all outputs.
    assert_eq!(names(&parse_fastq_records(&r1)), vec!["p1"], "R1: duplicate excluded");
    assert_eq!(names(&parse_fastq_records(&r2)), vec!["p1"], "R2: duplicate excluded");
    assert!(fs::read_to_string(&r0).expect("read r0").is_empty(), "out0: nothing routed here");

    // Require-flags is the sibling direction: `-f`/`--require-flags` keeps only
    // reads whose flags include every required bit, so requiring 0x400 keeps the
    // duplicate-flagged pair and drops the normal one — the inverse of the
    // exclude case above.
    let (rf1, rf2, rf0) =
        (dir.path().join("rf1.fq"), dir.path().join("rf2.fq"), dir.path().join("rf0.fq"));
    Fastq::try_parse_from([
        "fastq",
        "-i",
        input.to_str().unwrap(),
        "-f",
        "0x400", // require duplicates
        "-1",
        rf1.to_str().unwrap(),
        "-2",
        rf2.to_str().unwrap(),
        "-0",
        rf0.to_str().unwrap(),
    ])
    .expect("parse")
    .execute("fgumi fastq")
    .expect("require-flags paired fastq");

    // Only the duplicate-flagged pair matches the required flag.
    assert_eq!(names(&parse_fastq_records(&rf1)), vec!["dup"], "R1: only required kept");
    assert_eq!(names(&parse_fastq_records(&rf2)), vec!["dup"], "R2: only required kept");
    assert!(fs::read_to_string(&rf0).expect("read rf0").is_empty(), "out0: nothing routed here");
}

/// `-K`/`--bwa-chunk-size` is deprecated and ignored: a run passing it still
/// succeeds and produces the same output as one that omits it (the value no
/// longer affects batching, which the pipeline now manages).
#[test]
fn test_fastq_bwa_chunk_size_deprecated_but_accepted() {
    let dir = TempDir::new().expect("temp dir");
    let input = dir.path().join("in.bam");
    create_paired_bam(
        &input,
        vec![("read1", "ACGTACGT", "IIIIIIII", "TGCATGCA", "IIIIIIII", false)],
    );

    let with_k = dir.path().join("with_k.fq");
    Fastq::try_parse_from([
        "fastq",
        "-i",
        input.to_str().unwrap(),
        "-K",
        "12345", // arbitrary non-default value; must be accepted and ignored
        "-o",
        with_k.to_str().unwrap(),
    ])
    .expect("`-K <value>` must still parse")
    .execute("fgumi fastq")
    .expect("`-K` is deprecated but must not fail the run");

    let without_k = dir.path().join("without_k.fq");
    Fastq::try_parse_from([
        "fastq",
        "-i",
        input.to_str().unwrap(),
        "-o",
        without_k.to_str().unwrap(),
    ])
    .expect("parse")
    .execute("fgumi fastq")
    .expect("baseline run");

    // The ignored knob must not change the output. Guard the oracle: if a shared
    // defect made both runs emit nothing, a bare `read(&with_k) == read(&without_k)`
    // would pass vacuously, so assert the output is non-empty first.
    let with_k_bytes = fs::read(&with_k).expect("read with_k");
    assert!(!with_k_bytes.is_empty(), "oracle must be non-empty, not both-empty");
    assert_eq!(
        with_k_bytes,
        fs::read(&without_k).expect("read without_k"),
        "--bwa-chunk-size must not affect output (it is ignored)"
    );
}

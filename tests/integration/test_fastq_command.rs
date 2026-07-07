//! Integration tests for the fastq command.
//!
//! These tests invoke the fastq command in-process via `Command::execute()`.

use clap::Parser;
use fgumi_lib::commands::command::Command;
use fgumi_lib::commands::fastq::Fastq;
use fgumi_raw_bam::{SamBuilder, flags};
use fgumi_tag::SamTag;
use noodles::bam;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use std::fs;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, to_record_buf};

/// Create a BAM with a single read pair carrying a UMI in the given tag.
///
/// `umi_tag` is `None` to write no UMI tag at all (passthrough test).
fn create_pair_with_umi_tag(path: &Path, name: &str, umi_tag: Option<(SamTag, &str)>) {
    let header = create_minimal_header("chr1", 10000);
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");

    for (seq, segment_flag) in
        [("ACGTACGT", flags::FIRST_SEGMENT), ("TGCATGCA", flags::LAST_SEGMENT)]
    {
        let mut b = SamBuilder::new();
        b.read_name(name.as_bytes())
            .sequence(seq.as_bytes())
            .qualities(&[40, 40, 40, 40, 40, 40, 40, 40])
            .flags(flags::PAIRED | segment_flag)
            .ref_id(0)
            .pos(99)
            .mapq(60);
        if let Some((tag, value)) = umi_tag {
            b.add_string_tag(tag, value.as_bytes());
        }
        writer
            .write_alignment_record(&header, &to_record_buf(&b.build()))
            .expect("Failed to write record");
    }
    writer.try_finish().expect("Failed to finish BAM");
}

/// Run `fgumi fastq` with the given extra args and return parsed FASTQ records.
fn run_fastq_with_args(input_bam: &Path, extra: &[&str]) -> Vec<(String, String, String)> {
    let temp = input_bam.parent().expect("parent");
    let output_fq = temp.join("umi_out.fq");
    let mut args =
        vec!["fastq", "-i", input_bam.to_str().unwrap(), "-o", output_fq.to_str().unwrap()];
    args.extend_from_slice(extra);
    let cmd = Fastq::try_parse_from(args).expect("failed to parse fastq args");
    cmd.execute("fgumi fastq").expect("fastq command failed");
    parse_fastq_records(&output_fq)
}

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

/// `--umi-in-header` with defaults reproduces `samtools fastq -U`: the duplex
/// UMI stored as `AAAAAAAA-CCCCCCCC` in `RX` becomes `:AAAAAAAA+CCCCCCCC`
/// appended before the /1 /2 suffix.
#[test]
fn test_fastq_umi_in_header_dragen_format() {
    let temp_dir = TempDir::new().expect("temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    create_pair_with_umi_tag(&input_bam, "read1", Some((SamTag::RX, "AAAAAAAA-CCCCCCCC")));

    let records = run_fastq_with_args(&input_bam, &["--annotate-read-names"]);
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].0, "read1:AAAAAAAA+CCCCCCCC/1");
    assert_eq!(records[1].0, "read1:AAAAAAAA+CCCCCCCC/2");
}

/// `-n` drops the /1 /2 suffix but keeps the UMI annotation.
#[test]
fn test_fastq_umi_in_header_no_suffix() {
    let temp_dir = TempDir::new().expect("temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    create_pair_with_umi_tag(&input_bam, "read1", Some((SamTag::RX, "AAAAAAAA-CCCCCCCC")));

    let records = run_fastq_with_args(&input_bam, &["--annotate-read-names", "-n"]);
    assert_eq!(records[0].0, "read1:AAAAAAAA+CCCCCCCC");
    assert_eq!(records[1].0, "read1:AAAAAAAA+CCCCCCCC");
}

/// Custom delimiters: underscore name separator, empty UMI separator.
#[test]
fn test_fastq_umi_in_header_custom_delims() {
    let temp_dir = TempDir::new().expect("temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    create_pair_with_umi_tag(&input_bam, "read1", Some((SamTag::RX, "AAAA-GGGG")));

    let records = run_fastq_with_args(
        &input_bam,
        &["--annotate-read-names", "-n", "--umi-name-delim", "_", "--umi-sep", ""],
    );
    assert_eq!(records[0].0, "read1_AAAAGGGG");
}

/// When `RX` is absent the fallback `OX` tag is used (default tag list `RX,OX`).
#[test]
fn test_fastq_umi_in_header_ox_fallback() {
    let temp_dir = TempDir::new().expect("temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    create_pair_with_umi_tag(&input_bam, "read1", Some((SamTag::OX, "TTTT-GGGG")));

    let records = run_fastq_with_args(&input_bam, &["--annotate-read-names", "-n"]);
    assert_eq!(records[0].0, "read1:TTTT+GGGG");
}

/// With the default tag list `RX,OX`, `RX` takes priority when both are present.
#[test]
fn test_fastq_umi_in_header_rx_preferred_over_ox() {
    let temp_dir = TempDir::new().expect("temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    // Write both RX and OX; RX must win.
    let header = create_minimal_header("chr1", 10000);
    let mut writer = bam::io::Writer::new(fs::File::create(&input_bam).expect("create BAM"));
    writer.write_header(&header).expect("write header");
    let mut b = SamBuilder::new();
    b.read_name(b"read1")
        .sequence(b"ACGT")
        .qualities(&[40, 40, 40, 40])
        .flags(flags::PAIRED | flags::FIRST_SEGMENT)
        .ref_id(0)
        .pos(99)
        .mapq(60);
    b.add_string_tag(SamTag::RX, b"RRRR");
    b.add_string_tag(SamTag::OX, b"OOOO");
    writer.write_alignment_record(&header, &to_record_buf(&b.build())).expect("write");
    writer.try_finish().expect("finish");

    let records = run_fastq_with_args(&input_bam, &["--annotate-read-names", "-n"]);
    assert_eq!(records[0].0, "read1:RRRR");
}

/// A record with no UMI tag passes through with the name unchanged.
#[test]
fn test_fastq_umi_in_header_missing_tag_passthrough() {
    let temp_dir = TempDir::new().expect("temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    create_pair_with_umi_tag(&input_bam, "read1", None);

    let records = run_fastq_with_args(&input_bam, &["--annotate-read-names", "-n"]);
    assert_eq!(records[0].0, "read1");
}

/// Without `--umi-in-header`, the tag is ignored and the name is unchanged.
#[test]
fn test_fastq_umi_not_in_header_by_default() {
    let temp_dir = TempDir::new().expect("temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    create_pair_with_umi_tag(&input_bam, "read1", Some((SamTag::RX, "AAAA-GGGG")));

    let records = run_fastq_with_args(&input_bam, &["-n"]);
    assert_eq!(records[0].0, "read1");
}

/// Read a FASTQ file, transparently decompressing BGZF/gzip (`.gz`/`.bgz`) via
/// multi-member gzip decoding (BGZF is a series of gzip blocks). Both extensions
/// trigger BGZF writing on the output side, so the reader honors both.
fn read_fastq_maybe_gz(path: &Path) -> Vec<(String, String, String)> {
    let file = fs::File::open(path).expect("open fastq");
    let is_gz = matches!(path.extension().and_then(|e| e.to_str()), Some("gz" | "bgz"));
    let reader: Box<dyn BufRead> = if is_gz {
        Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    let lines: Vec<String> = reader.lines().map(|l| l.unwrap()).collect();
    let mut records = Vec::new();
    for chunk in lines.chunks(4) {
        if chunk.len() == 4 {
            records.push((
                chunk[0].trim_start_matches('@').to_string(),
                chunk[1].clone(),
                chunk[3].clone(),
            ));
        }
    }
    records
}

/// Paired `-1`/`-2` output routes R1 to out1 and R2 to out2, with the UMI
/// annotation applied to both.
#[test]
fn test_fastq_paired_split_output() {
    let temp_dir = TempDir::new().expect("temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    create_pair_with_umi_tag(&input_bam, "read1", Some((SamTag::RX, "AAAA-GGGG")));
    let r1 = temp_dir.path().join("R1.fq");
    let r2 = temp_dir.path().join("R2.fq");

    let cmd = Fastq::try_parse_from([
        "fastq",
        "-i",
        input_bam.to_str().unwrap(),
        "--annotate-read-names",
        "-1",
        r1.to_str().unwrap(),
        "-2",
        r2.to_str().unwrap(),
    ])
    .expect("parse");
    cmd.execute("fgumi fastq").expect("run");

    let recs1 = read_fastq_maybe_gz(&r1);
    let recs2 = read_fastq_maybe_gz(&r2);
    assert_eq!(recs1.len(), 1, "one R1 read");
    assert_eq!(recs2.len(), 1, "one R2 read");
    // Split mode omits the /1 /2 suffix (like samtools -1/-2); R1 and R2 names match.
    assert_eq!(recs1[0].0, "read1:AAAA+GGGG");
    assert_eq!(recs2[0].0, "read1:AAAA+GGGG");
    // R1 sequence is the FIRST_SEGMENT read; R2 the LAST_SEGMENT read.
    assert_eq!(recs1[0].1, "ACGTACGT");
    assert_eq!(recs2[0].1, "TGCATGCA");
}

/// A `.gz` output path is written as BGZF (gzip magic bytes + BGZF `BC` subfield)
/// and round-trips to the expected records.
#[test]
fn test_fastq_paired_output_is_bgzf() {
    let temp_dir = TempDir::new().expect("temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    create_pair_with_umi_tag(&input_bam, "read1", None);
    let r1 = temp_dir.path().join("R1.fq.gz");
    let r2 = temp_dir.path().join("R2.fq.gz");

    let cmd = Fastq::try_parse_from([
        "fastq",
        "-i",
        input_bam.to_str().unwrap(),
        "-1",
        r1.to_str().unwrap(),
        "-2",
        r2.to_str().unwrap(),
    ])
    .expect("parse");
    cmd.execute("fgumi fastq").expect("run");

    // Verify BGZF framing: gzip magic 1f 8b, and the BGZF `BC` extra subfield.
    let raw = fs::read(&r1).expect("read gz");
    assert!(raw.len() > 18, "non-empty bgzf");
    assert_eq!(&raw[0..2], &[0x1f, 0x8b], "gzip magic");
    assert_eq!(raw[3] & 0x04, 0x04, "FEXTRA flag set (BGZF has an extra field)");
    assert_eq!([raw[12], raw[13]], [b'B', b'C'], "BGZF BC subfield identifier");

    // And it decodes to the expected read (split mode omits the /1 suffix).
    let recs1 = read_fastq_maybe_gz(&r1);
    assert_eq!(recs1.len(), 1);
    assert_eq!(recs1[0].0, "read1");
}

/// The `.bgz` extension triggers BGZF writing on the same footing as `.gz`,
/// so a `.bgz`-suffixed paired output is BGZF-framed and round-trips. Guards
/// against a writer bug specific to `.bgz` extension detection.
#[test]
fn test_fastq_paired_output_bgz_extension_is_bgzf() {
    let temp_dir = TempDir::new().expect("temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    create_pair_with_umi_tag(&input_bam, "read1", None);
    let r1 = temp_dir.path().join("R1.fq.bgz");
    let r2 = temp_dir.path().join("R2.fq.bgz");

    let cmd = Fastq::try_parse_from([
        "fastq",
        "-i",
        input_bam.to_str().unwrap(),
        "-1",
        r1.to_str().unwrap(),
        "-2",
        r2.to_str().unwrap(),
    ])
    .expect("parse");
    cmd.execute("fgumi fastq").expect("run");

    // Verify BGZF framing on the `.bgz` output: gzip magic + BGZF `BC` subfield.
    let raw = fs::read(&r1).expect("read bgz");
    assert!(raw.len() > 18, "non-empty bgzf");
    assert_eq!(&raw[0..2], &[0x1f, 0x8b], "gzip magic");
    assert_eq!(raw[3] & 0x04, 0x04, "FEXTRA flag set (BGZF has an extra field)");
    assert_eq!([raw[12], raw[13]], [b'B', b'C'], "BGZF BC subfield identifier");

    // And both `.bgz` mates decode to the expected read (split mode omits /1 /2).
    let recs1 = read_fastq_maybe_gz(&r1);
    let recs2 = read_fastq_maybe_gz(&r2);
    assert_eq!(recs1.len(), 1);
    assert_eq!(recs2.len(), 1);
    assert_eq!(recs1[0].0, "read1");
    assert_eq!(recs2[0].0, "read1");
}

/// An interleaved `--output` ending in `.gz` is also written as BGZF.
#[test]
fn test_fastq_interleaved_gz_is_bgzf() {
    let temp_dir = TempDir::new().expect("temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    create_pair_with_umi_tag(&input_bam, "read1", None);
    let out = temp_dir.path().join("out.fq.gz");

    let cmd = Fastq::try_parse_from([
        "fastq",
        "-i",
        input_bam.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ])
    .expect("parse");
    cmd.execute("fgumi fastq").expect("run");

    let raw = fs::read(&out).expect("read gz");
    assert_eq!(&raw[0..2], &[0x1f, 0x8b], "gzip magic");
    assert_eq!([raw[12], raw[13]], [b'B', b'C'], "BGZF BC subfield");
    let recs = read_fastq_maybe_gz(&out);
    assert_eq!(recs.len(), 2, "both mates interleaved");
    // Assert record identity, not just count: R1 then R2, with the /1 /2 suffix.
    assert_eq!(recs[0].0, "read1/1");
    assert_eq!(recs[1].0, "read1/2");
    assert_eq!(recs[0].1, "ACGTACGT");
    assert_eq!(recs[1].1, "TGCATGCA");
}

/// `--out1` without `--out2` is rejected by clap (`requires`).
#[test]
fn test_fastq_out1_requires_out2() {
    let temp_dir = TempDir::new().expect("temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let r1 = temp_dir.path().join("R1.fq");
    let parsed = Fastq::try_parse_from([
        "fastq",
        "-i",
        input_bam.to_str().unwrap(),
        "-1",
        r1.to_str().unwrap(),
    ]);
    assert!(parsed.is_err(), "--out1 without --out2 must be rejected");
}

/// `-U` remains a working alias for `--annotate-read-names`.
#[test]
fn test_fastq_annotate_u_alias() {
    let temp_dir = TempDir::new().expect("temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    create_pair_with_umi_tag(&input_bam, "read1", Some((SamTag::RX, "AAAA-GGGG")));
    let records = run_fastq_with_args(&input_bam, &["-U", "-n"]);
    assert_eq!(records[0].0, "read1:AAAA+GGGG");
}

/// Two outputs pointing at the same file must be rejected before any writing,
/// so two writer streams cannot silently truncate/corrupt one file.
#[test]
fn test_fastq_duplicate_output_paths_rejected() {
    let temp_dir = TempDir::new().expect("temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    create_pair_with_umi_tag(&input_bam, "read1", None);
    let same = temp_dir.path().join("both.fq");

    let cmd = Fastq::try_parse_from([
        "fastq",
        "-i",
        input_bam.to_str().unwrap(),
        "-1",
        same.to_str().unwrap(),
        "-2",
        same.to_str().unwrap(),
    ])
    .expect("parse");
    let err = cmd.execute("fgumi fastq").expect_err("must reject --out1 == --out2");
    assert!(err.to_string().contains("used more than once"), "unexpected error: {err}");
    assert!(!same.exists(), "no output file should be created when paths collide");
}

/// Two lexically-different but equivalent output paths that don't exist yet —
/// here the same file reached through a real directory and a symlink to it —
/// must still be rejected. The dedup canonicalizes the (existing) parent
/// directory, so the aliases collide before two sink threads can truncate the
/// same not-yet-existing target. Without parent-canonicalization the raw
/// lexical paths differ and the collision slips through.
#[cfg(unix)]
#[test]
fn test_fastq_aliased_nonexistent_output_paths_rejected() {
    let temp_dir = TempDir::new().expect("temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    create_pair_with_umi_tag(&input_bam, "read1", None);

    // `real/` exists; `link/` is a symlink to it. Same not-yet-existing file
    // `r.fq` reached two lexically-distinct ways that resolve to one inode.
    let real_dir = temp_dir.path().join("real");
    fs::create_dir(&real_dir).expect("mkdir real");
    let link_dir = temp_dir.path().join("link");
    std::os::unix::fs::symlink(&real_dir, &link_dir).expect("symlink");
    let out1 = real_dir.join("r.fq");
    let out2 = link_dir.join("r.fq");
    assert!(!out1.exists(), "target must not exist yet for this to test the non-existent path");
    assert_ne!(out1, out2, "the two paths must be lexically distinct");

    let cmd = Fastq::try_parse_from([
        "fastq",
        "-i",
        input_bam.to_str().unwrap(),
        "-1",
        out1.to_str().unwrap(),
        "-2",
        out2.to_str().unwrap(),
    ])
    .expect("parse");
    let err = cmd.execute("fgumi fastq").expect_err("must reject aliased non-existent outputs");
    assert!(err.to_string().contains("used more than once"), "unexpected error: {err}");
    assert!(!out1.exists(), "no output file should be created when paths collide");
}

/// `/dev/null` is exempt from the duplicate-path check: routing several streams
/// to it (e.g. `--out0 /dev/null`) is intentional, as with `samtools fastq`.
#[cfg(unix)]
#[test]
fn test_fastq_dev_null_out0_allowed() {
    let temp_dir = TempDir::new().expect("temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    create_pair_with_umi_tag(&input_bam, "read1", None);
    let r1 = temp_dir.path().join("R1.fq");
    let r2 = temp_dir.path().join("R2.fq");

    let cmd = Fastq::try_parse_from([
        "fastq",
        "-i",
        input_bam.to_str().unwrap(),
        "-1",
        r1.to_str().unwrap(),
        "-2",
        r2.to_str().unwrap(),
        "-0",
        "/dev/null",
    ])
    .expect("parse");
    cmd.execute("fgumi fastq").expect("run with --out0 /dev/null");
    assert_eq!(read_fastq_maybe_gz(&r1).len(), 1);
    assert_eq!(read_fastq_maybe_gz(&r2).len(), 1);
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
    assert!(err.to_string().contains("must differ"), "unexpected error message: {err}");

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
    assert!(err.to_string().contains("must differ"), "unexpected error message: {err}");

    // The input BAM must not have been truncated through the symlink.
    let input_size_after = std::fs::metadata(&input_bam).expect("stat input").len();
    assert_eq!(
        input_size_before, input_size_after,
        "input BAM was truncated/clobbered through symlinked --output"
    );
}

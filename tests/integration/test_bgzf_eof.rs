//! Integration tests verifying that BAM pipeline output files contain
//! the standard 28-byte BGZF EOF marker block.
//!
//! GitHub issue #125: BAM files written by pipeline-based commands were missing
//! the BGZF EOF block, causing `samtools quickcheck` to fail.

use fgumi_raw_bam::{RawRecord, SamBuilder, flags};
use noodles::bam;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value;
use std::fs;
use std::path::Path;
use std::process::{Command, Stdio};
use tempfile::TempDir;

use crate::helpers::assertions::assert_has_bgzf_eof;
use crate::helpers::bam_generator::{
    create_minimal_header, create_test_reference, create_umi_family, to_record_buf,
};

/// Write a BAM file with the given records using the standard test header.
fn write_test_bam(path: &Path, records: &[RawRecord]) {
    let header = create_minimal_header("chr1", 10000);
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");
    for record in records {
        writer
            .write_alignment_record(&header, &to_record_buf(record))
            .expect("Failed to write record");
    }
    writer.try_finish().expect("Failed to finish BAM");
}

/// Write a BAM file with MI-tagged records (grouped input for consensus callers).
fn write_grouped_bam(path: &Path, families: &[(&str, &[RawRecord])]) {
    let header = create_minimal_header("chr1", 10000);
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");

    let mi_tag = Tag::from(fgumi_lib::sam::SamTag::MI);
    for &(mi, records) in families {
        for record in records {
            let mut buf = to_record_buf(record);
            buf.data_mut().insert(mi_tag, Value::from(mi));
            writer.write_alignment_record(&header, &buf).expect("Failed to write record");
        }
    }
    writer.try_finish().expect("Failed to finish BAM");
}

/// Create paired-end duplicate reads at a given position with MC tags.
fn create_dedup_reads(name: &[u8], umi: &str, start: i32) -> Vec<RawRecord> {
    let r1 = {
        let mut b = SamBuilder::new();
        b.read_name(name)
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE)
            .ref_id(0)
            .pos(start - 1)
            .mapq(60)
            .cigar_ops(&[8 << 4]) // 8M
            .mate_ref_id(0)
            .mate_pos(start + 99)
            .template_length(108)
            .add_string_tag(b"RX", umi.as_bytes())
            .add_string_tag(b"MC", b"8M");
        b.build()
    };

    let r2 = {
        let mut b = SamBuilder::new();
        b.read_name(name)
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE)
            .ref_id(0)
            .pos(start + 99)
            .mapq(60)
            .cigar_ops(&[8 << 4]) // 8M
            .mate_ref_id(0)
            .mate_pos(start - 1)
            .template_length(-108)
            .add_string_tag(b"RX", umi.as_bytes())
            .add_string_tag(b"MC", b"8M");
        b.build()
    };

    vec![r1, r2]
}

/// Create a duplex read pair with /A or /B strand MI suffix.
fn create_duplex_pair(
    name: &[u8],
    mi_tag: &str,
    sequence: &[u8],
    ref_start: i32,
    is_b_strand: bool,
) -> (RawRecord, RawRecord) {
    let read_len = sequence.len();
    let cigar_op = u32::try_from(read_len).expect("read_len fits u32") << 4; // NM
    let (r1_start, r2_start, r1_rev, r2_rev) = if is_b_strand {
        (ref_start + 100, ref_start, true, false)
    } else {
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
        b.read_name(name)
            .sequence(sequence)
            .qualities(&vec![30; read_len])
            .flags(r1_flags)
            .ref_id(0)
            .pos(r1_start - 1)
            .mapq(60)
            .cigar_ops(&[cigar_op])
            .mate_ref_id(0)
            .mate_pos(r2_start - 1)
            .template_length(tlen)
            .add_string_tag(b"MI", mi_tag.as_bytes());
        b.build()
    };

    let r2 = {
        let mut b = SamBuilder::new();
        b.read_name(name)
            .sequence(sequence)
            .qualities(&vec![30; read_len])
            .flags(r2_flags)
            .ref_id(0)
            .pos(r2_start - 1)
            .mapq(60)
            .cigar_ops(&[cigar_op])
            .mate_ref_id(0)
            .mate_pos(r1_start - 1)
            .template_length(-tlen)
            .add_string_tag(b"MI", mi_tag.as_bytes());
        b.build()
    };

    (r1, r2)
}

/// Create a single-template family (will be rejected by --min-reads > 1).
fn create_rejected_family(name: &[u8], umi: &str, start: i32) -> Vec<RawRecord> {
    let r1 = {
        let mut b = SamBuilder::new();
        b.read_name(name)
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE)
            .ref_id(0)
            .pos(start - 1)
            .mapq(60)
            .cigar_ops(&[8 << 4]) // 8M
            .mate_ref_id(0)
            .mate_pos(start + 99)
            .template_length(108)
            .add_string_tag(b"RX", umi.as_bytes())
            .add_string_tag(b"MC", b"8M");
        b.build()
    };

    let r2 = {
        let mut b = SamBuilder::new();
        b.read_name(name)
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE)
            .ref_id(0)
            .pos(start + 99)
            .mapq(60)
            .cigar_ops(&[8 << 4]) // 8M
            .mate_ref_id(0)
            .mate_pos(start - 1)
            .template_length(-108)
            .add_string_tag(b"RX", umi.as_bytes())
            .add_string_tag(b"MC", b"8M");
        b.build()
    };

    vec![r1, r2]
}

/// Create a CODEC read pair (FR orientation, same position).
fn create_codec_pair(name: &[u8], umi: &str, ref_start: i32) -> (RawRecord, RawRecord) {
    let r1 = {
        let mut b = SamBuilder::new();
        b.read_name(name)
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE)
            .ref_id(0)
            .pos(ref_start - 1)
            .mapq(60)
            .cigar_ops(&[8 << 4]) // 8M
            .mate_ref_id(0)
            .mate_pos(ref_start - 1)
            .template_length(8)
            .add_string_tag(b"MI", umi.as_bytes())
            .add_string_tag(b"MC", b"8M");
        b.build()
    };

    let r2 = {
        let mut b = SamBuilder::new();
        b.read_name(name)
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE)
            .ref_id(0)
            .pos(ref_start - 1)
            .mapq(60)
            .cigar_ops(&[8 << 4]) // 8M
            .mate_ref_id(0)
            .mate_pos(ref_start - 1)
            .template_length(-8)
            .add_string_tag(b"MI", umi.as_bytes())
            .add_string_tag(b"MC", b"8M");
        b.build()
    };

    (r1, r2)
}

// ============================================================================
// Tests
// ============================================================================

#[test]
fn test_group_output_has_bgzf_eof() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    let family = create_umi_family("AAAAAAAA", 5, "fam", "ACGTACGT", 30);
    write_test_bam(&input_bam, &family);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "group",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--strategy",
            "identity",
            "--edits",
            "0",
            "--threads",
            "2",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run group command");

    assert!(status.success(), "group command failed");
    assert_has_bgzf_eof(&output_bam);
}

#[test]
fn test_dedup_output_has_bgzf_eof() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    let mut records = create_dedup_reads(b"dup_0", "ACGTACGT", 100);
    records.extend(create_dedup_reads(b"dup_1", "ACGTACGT", 100));
    write_test_bam(&input_bam, &records);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "dedup",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--strategy",
            "identity",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run dedup command");

    assert!(status.success(), "dedup command failed");
    assert_has_bgzf_eof(&output_bam);
}

#[test]
fn test_simplex_output_has_bgzf_eof() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    let family = create_umi_family("ACGT", 5, "fam", "ACGTACGT", 30);
    write_grouped_bam(&input_bam, &[("1", &family)]);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "simplex",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--min-reads",
            "2",
            "--threads",
            "2",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run simplex command");

    assert!(status.success(), "simplex command failed");
    assert_has_bgzf_eof(&output_bam);
}

#[test]
fn test_duplex_output_has_bgzf_eof() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    // Create a duplex molecule: 3 /A pairs and 3 /B pairs
    let mut records = Vec::new();
    for i in 0..3 {
        let (r1, r2) =
            create_duplex_pair(format!("a_{i}").as_bytes(), "1/A", b"ACGTACGT", 100, false);
        records.push(r1);
        records.push(r2);
    }
    for i in 0..3 {
        let (r1, r2) =
            create_duplex_pair(format!("b_{i}").as_bytes(), "1/B", b"ACGTACGT", 100, true);
        records.push(r1);
        records.push(r2);
    }
    write_test_bam(&input_bam, &records);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "duplex",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--min-reads",
            "1",
            "--threads",
            "2",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run duplex command");

    assert!(status.success(), "duplex command failed");
    assert_has_bgzf_eof(&output_bam);
}

#[test]
fn test_codec_output_has_bgzf_eof() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    let mut records = Vec::new();
    for i in 0..3 {
        let (r1, r2) = create_codec_pair(format!("read{i}").as_bytes(), "1", 100);
        records.push(r1);
        records.push(r2);
    }
    write_test_bam(&input_bam, &records);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "codec",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--min-reads",
            "1",
            "--threads",
            "2",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run codec command");

    assert!(status.success(), "codec command failed");
    assert_has_bgzf_eof(&output_bam);
}

#[test]
fn test_filter_output_has_bgzf_eof() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let ref_path = create_test_reference(temp_dir.path());

    let record = {
        let mut b = SamBuilder::new();
        b.read_name(b"cons1")
            .sequence(b"ACGTACGT")
            .qualities(&[35; 8])
            .ref_id(0)
            .pos(99)
            .mapq(60)
            .cigar_ops(&[8 << 4]) // 8M
            .add_array_u16(b"cd", &[10; 8])
            .add_array_u16(b"ce", &[0; 8]);
        b.build()
    };

    write_test_bam(&input_bam, &[record]);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "filter",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--ref",
            ref_path.to_str().unwrap(),
            "--min-reads",
            "1",
            "--max-no-call-fraction",
            "1.0",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run filter command");

    assert!(status.success(), "filter command failed");
    assert_has_bgzf_eof(&output_bam);
}

#[test]
fn test_clip_output_has_bgzf_eof() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let ref_path = create_test_reference(temp_dir.path());

    let r1 = {
        let mut b = SamBuilder::new();
        b.read_name(b"read1")
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE)
            .ref_id(0)
            .pos(99)
            .mapq(60)
            .cigar_ops(&[8 << 4]) // 8M
            .mate_ref_id(0)
            .mate_pos(103)
            .template_length(12);
        b.build()
    };

    let r2 = {
        let mut b = SamBuilder::new();
        b.read_name(b"read1")
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE)
            .ref_id(0)
            .pos(103)
            .mapq(60)
            .cigar_ops(&[8 << 4]) // 8M
            .mate_ref_id(0)
            .mate_pos(99)
            .template_length(-12);
        b.build()
    };

    write_test_bam(&input_bam, &[r1, r2]);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "clip",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--reference",
            ref_path.to_str().unwrap(),
            "--read-one-five-prime",
            "1",
            "--read-one-three-prime",
            "1",
            "--threads",
            "2",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run clip command");

    assert!(status.success(), "clip command failed");
    assert_has_bgzf_eof(&output_bam);
}

#[test]
fn test_correct_output_has_bgzf_eof() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let whitelist = temp_dir.path().join("whitelist.txt");

    let reads = create_umi_family("ACGTACGT", 3, "read", "AAAAGGGG", 30);
    write_test_bam(&input_bam, &reads);
    fs::write(&whitelist, "ACGTACGT").expect("Failed to write whitelist");

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "correct",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--umi-files",
            whitelist.to_str().unwrap(),
            "--max-mismatches",
            "1",
            "--min-distance",
            "1",
            "--threads",
            "2",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run correct command");

    assert!(status.success(), "correct command failed");
    assert_has_bgzf_eof(&output_bam);
}

/// Regression test for the streaming pipe scenario from issue #125:
/// `fgumi simplex ... -o /dev/stdout | fgumi filter -i /dev/stdin ...`
///
/// Verifies that piped BAM output contains the BGZF EOF marker when the
/// downstream command writes to a file.
#[test]
#[allow(clippy::zombie_processes)]
fn test_piped_simplex_to_filter_has_bgzf_eof() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let ref_path = create_test_reference(temp_dir.path());

    // Create grouped input with consensus-worthy families
    let family = create_umi_family("ACGT", 5, "fam", "ACGTACGT", 30);
    write_grouped_bam(&input_bam, &[("1", &family)]);

    // Spawn simplex writing to stdout
    let mut simplex = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "simplex",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            "/dev/stdout",
            "--min-reads",
            "2",
            "--threads",
            "2",
            "--compression-level",
            "1",
        ])
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .expect("Failed to spawn simplex");

    // Spawn filter reading from simplex stdout
    let filter = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "filter",
            "--input",
            "/dev/stdin",
            "--output",
            output_bam.to_str().unwrap(),
            "--ref",
            ref_path.to_str().unwrap(),
            "--min-reads",
            "1",
            "--max-no-call-fraction",
            "1.0",
            "--compression-level",
            "1",
        ])
        .stdin(simplex.stdout.take().unwrap())
        .stderr(Stdio::null())
        .output()
        .expect("Failed to run filter");

    let simplex_status = simplex.wait().expect("Failed to wait for simplex");
    assert!(simplex_status.success(), "simplex command failed in pipe");
    assert!(filter.status.success(), "filter command failed in pipe");
    assert_has_bgzf_eof(&output_bam);
}

// ============================================================================
// Rejects writer BGZF EOF tests
// ============================================================================

#[test]
fn test_simplex_rejects_has_bgzf_eof() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let rejects_bam = temp_dir.path().join("rejects.bam");

    // One family with 3 reads (kept) and one singleton family (rejected by --min-reads 2)
    let kept = create_umi_family("ACGT", 3, "kept", "ACGTACGT", 30);
    let rejected = create_rejected_family(b"reject", "TTTT", 500);
    write_grouped_bam(&input_bam, &[("1", &kept), ("2", &rejected)]);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "simplex",
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
        .expect("Failed to run simplex command");

    assert!(status.success(), "simplex command with rejects failed");
    assert_has_bgzf_eof(&output_bam);
    assert_has_bgzf_eof(&rejects_bam);
}

#[test]
fn test_duplex_rejects_has_bgzf_eof() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let rejects_bam = temp_dir.path().join("rejects.bam");

    // 3 /A pairs and 3 /B pairs (kept), plus a singleton /A pair (rejected by --min-reads 2)
    let mut records = Vec::new();
    for i in 0..3 {
        let (r1, r2) =
            create_duplex_pair(format!("a_{i}").as_bytes(), "1/A", b"ACGTACGT", 100, false);
        records.push(r1);
        records.push(r2);
    }
    for i in 0..3 {
        let (r1, r2) =
            create_duplex_pair(format!("b_{i}").as_bytes(), "1/B", b"ACGTACGT", 100, true);
        records.push(r1);
        records.push(r2);
    }
    // Singleton family that will be rejected
    let (r1, r2) = create_duplex_pair(b"reject_a", "2/A", b"ACGTACGT", 500, false);
    records.push(r1);
    records.push(r2);
    write_test_bam(&input_bam, &records);

    // Use single-threaded mode: duplex multi-threaded pipeline doesn't collect
    // raw reject bytes, so rejects writing only happens in the single-threaded path.
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
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run duplex command");

    assert!(status.success(), "duplex command with rejects failed");
    assert_has_bgzf_eof(&output_bam);
    assert_has_bgzf_eof(&rejects_bam);
}

#[test]
fn test_codec_rejects_has_bgzf_eof() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let rejects_bam = temp_dir.path().join("rejects.bam");

    // One family with 3 reads (kept) and one singleton (rejected by --min-reads 2)
    let mut records = Vec::new();
    for i in 0..3 {
        let (r1, r2) = create_codec_pair(format!("read{i}").as_bytes(), "1", 100);
        records.push(r1);
        records.push(r2);
    }
    let (r1, r2) = create_codec_pair(b"reject", "2", 500);
    records.push(r1);
    records.push(r2);
    write_test_bam(&input_bam, &records);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "codec",
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
        .expect("Failed to run codec command");

    assert!(status.success(), "codec command with rejects failed");
    assert_has_bgzf_eof(&output_bam);
    assert_has_bgzf_eof(&rejects_bam);
}

#[test]
fn test_correct_rejects_has_bgzf_eof() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let rejects_bam = temp_dir.path().join("rejects.bam");
    let whitelist = temp_dir.path().join("whitelist.txt");

    // Reads with matching UMI (kept) and non-matching UMI (rejected)
    let mut records = create_umi_family("ACGTACGT", 3, "kept", "AAAAGGGG", 30);
    records.extend(create_umi_family("NNNNNNNN", 2, "reject", "AAAAGGGG", 30));
    write_test_bam(&input_bam, &records);
    fs::write(&whitelist, "ACGTACGT").expect("Failed to write whitelist");

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "correct",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--rejects",
            rejects_bam.to_str().unwrap(),
            "--umi-files",
            whitelist.to_str().unwrap(),
            "--max-mismatches",
            "1",
            "--min-distance",
            "1",
            "--threads",
            "2",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run correct command");

    assert!(status.success(), "correct command with rejects failed");
    assert_has_bgzf_eof(&output_bam);
    assert_has_bgzf_eof(&rejects_bam);
}

#[test]
fn test_correct_single_threaded_rejects_has_bgzf_eof() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let rejects_bam = temp_dir.path().join("rejects.bam");
    let whitelist = temp_dir.path().join("whitelist.txt");

    let mut records = create_umi_family("ACGTACGT", 3, "kept", "AAAAGGGG", 30);
    records.extend(create_umi_family("NNNNNNNN", 2, "reject", "AAAAGGGG", 30));
    write_test_bam(&input_bam, &records);
    fs::write(&whitelist, "ACGTACGT").expect("Failed to write whitelist");

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "correct",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--rejects",
            rejects_bam.to_str().unwrap(),
            "--umi-files",
            whitelist.to_str().unwrap(),
            "--max-mismatches",
            "1",
            "--min-distance",
            "1",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run correct command");

    assert!(status.success(), "correct single-threaded with rejects failed");
    assert_has_bgzf_eof(&output_bam);
    assert_has_bgzf_eof(&rejects_bam);
}

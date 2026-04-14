//! End-to-end CLI tests for the zipper command.
//!
//! These tests run the actual `fgumi zipper` binary and validate:
//! 1. Basic merge of unmapped + mapped BAMs with tag transfer
//! 2. Tag removal via `--tags-to-remove`
//! 3. Error on missing input files

use fgumi_lib::sam::SamTag;
use fgumi_lib::sam::builder::RecordBuilder;
use noodles::bam;
use noodles::sam;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::RecordBuf;
use std::fs;
use std::io::Write;
use std::path::Path;
use std::process::{Command, Stdio};
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, create_test_reference};

/// Create a queryname-sorted unmapped BAM with UMI tags (RX/QX).
fn create_unmapped_bam(path: &Path, records: &[RecordBuf]) {
    let header = sam::Header::default();
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create unmapped BAM"));
    writer.write_header(&header).expect("Failed to write header");
    for record in records {
        writer.write_alignment_record(&header, record).expect("Failed to write record");
    }
    writer.finish(&header).expect("Failed to finish BAM");
}

/// Create a mapped SAM file with aligned reads (same read names as unmapped).
fn create_mapped_sam(path: &Path, header: &sam::Header, records: &[RecordBuf]) {
    let file = fs::File::create(path).expect("Failed to create mapped SAM");
    let mut writer = sam::io::Writer::new(file);
    writer.write_header(header).expect("Failed to write header");
    for record in records {
        writer.write_alignment_record(header, record).expect("Failed to write record");
    }
}

/// Create a mapped BAM file with aligned reads (same read names as unmapped).
fn create_mapped_bam(path: &Path, header: &sam::Header, records: &[RecordBuf]) {
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create mapped BAM"));
    writer.write_header(header).expect("Failed to write header");
    for record in records {
        writer.write_alignment_record(header, record).expect("Failed to write record");
    }
    writer.finish(header).expect("Failed to finish BAM");
}

/// Test basic zipper merge — unmapped tags are transferred to mapped reads.
#[test]
fn test_zipper_basic_merge() {
    let temp_dir = TempDir::new().unwrap();
    let unmapped_bam = temp_dir.path().join("unmapped.bam");
    let mapped_sam = temp_dir.path().join("mapped.sam");
    let output_bam = temp_dir.path().join("output.bam");
    let ref_path = create_test_reference(temp_dir.path());

    // Unmapped reads with RX/QX tags
    let unmapped_records = vec![
        RecordBuilder::new()
            .name("read1")
            .sequence("ACGTACGT")
            .qualities(&[30; 8])
            .unmapped(true)
            .tag("RX", "AACCGGTT")
            .tag("QX", "IIIIIIII")
            .build(),
        RecordBuilder::new()
            .name("read2")
            .sequence("TGCATGCA")
            .qualities(&[30; 8])
            .unmapped(true)
            .tag("RX", "GGTTCCAA")
            .tag("QX", "IIIIIIII")
            .build(),
    ];
    create_unmapped_bam(&unmapped_bam, &unmapped_records);

    // Mapped reads (same names, aligned to chr1)
    let mapped_header = create_minimal_header("chr1", 10000);
    let mapped_records = vec![
        RecordBuilder::new()
            .name("read1")
            .sequence("ACGTACGT")
            .qualities(&[30; 8])
            .reference_sequence_id(0)
            .alignment_start(100)
            .mapping_quality(60)
            .cigar("8M")
            .build(),
        RecordBuilder::new()
            .name("read2")
            .sequence("TGCATGCA")
            .qualities(&[30; 8])
            .reference_sequence_id(0)
            .alignment_start(200)
            .mapping_quality(60)
            .cigar("8M")
            .build(),
    ];
    create_mapped_sam(&mapped_sam, &mapped_header, &mapped_records);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "zipper",
            "--input",
            mapped_sam.to_str().unwrap(),
            "--unmapped",
            unmapped_bam.to_str().unwrap(),
            "--reference",
            ref_path.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run zipper command");

    assert!(status.success(), "Zipper command failed");
    assert!(output_bam.exists(), "Output BAM not created");

    // Verify output records have UMI tags transferred
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let header = reader.read_header().unwrap();
    let records: Vec<RecordBuf> = reader.record_bufs(&header).map(|r| r.unwrap()).collect();
    assert_eq!(records.len(), 2, "Should have 2 records in output");

    let rx_tag = Tag::from(SamTag::RX);
    for record in &records {
        assert!(record.data().get(&rx_tag).is_some(), "Output record should have RX tag");
    }
}

/// Test zipper with `--tags-to-remove` strips specified tags.
#[test]
fn test_zipper_tag_removal() {
    let temp_dir = TempDir::new().unwrap();
    let unmapped_bam = temp_dir.path().join("unmapped.bam");
    let mapped_sam = temp_dir.path().join("mapped.sam");
    let output_bam = temp_dir.path().join("output.bam");
    let ref_path = create_test_reference(temp_dir.path());

    // Unmapped read with RX and a custom tag XY
    let unmapped_records = vec![
        RecordBuilder::new()
            .name("read1")
            .sequence("ACGTACGT")
            .qualities(&[30; 8])
            .unmapped(true)
            .tag("RX", "AACCGGTT")
            .tag("XY", "REMOVE_ME")
            .build(),
    ];
    create_unmapped_bam(&unmapped_bam, &unmapped_records);

    let mapped_header = create_minimal_header("chr1", 10000);
    let mapped_records = vec![
        RecordBuilder::new()
            .name("read1")
            .sequence("ACGTACGT")
            .qualities(&[30; 8])
            .reference_sequence_id(0)
            .alignment_start(100)
            .mapping_quality(60)
            .cigar("8M")
            .build(),
    ];
    create_mapped_sam(&mapped_sam, &mapped_header, &mapped_records);

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "zipper",
            "--input",
            mapped_sam.to_str().unwrap(),
            "--unmapped",
            unmapped_bam.to_str().unwrap(),
            "--reference",
            ref_path.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--tags-to-remove",
            "XY",
            "--compression-level",
            "1",
        ])
        .status()
        .expect("Failed to run zipper command");

    assert!(status.success(), "Zipper command with --tags-to-remove failed");

    // Verify XY tag was removed but RX tag was kept
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let header = reader.read_header().unwrap();
    let records: Vec<RecordBuf> = reader.record_bufs(&header).map(|r| r.unwrap()).collect();
    assert_eq!(records.len(), 1);

    let rx_tag = Tag::from(SamTag::RX);
    let xy_tag = Tag::from(SamTag::new(b'X', b'Y'));
    assert!(records[0].data().get(&rx_tag).is_some(), "RX tag should be present");
    assert!(records[0].data().get(&xy_tag).is_none(), "XY tag should have been removed");
}

/// Test zipper errors with non-zero exit code for missing input file.
#[test]
fn test_zipper_missing_input() {
    let temp_dir = TempDir::new().unwrap();
    let output_bam = temp_dir.path().join("output.bam");
    let ref_path = create_test_reference(temp_dir.path());

    let missing_mapped = temp_dir.path().join("missing.mapped.sam");
    let missing_unmapped = temp_dir.path().join("missing.unmapped.bam");

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "zipper",
            "--input",
            missing_mapped.to_str().unwrap(),
            "--unmapped",
            missing_unmapped.to_str().unwrap(),
            "--reference",
            ref_path.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
        ])
        .status()
        .expect("Failed to run zipper command");

    assert!(!status.success(), "Zipper should fail for nonexistent input");
}

/// Test zipper accepts BAM as mapped input (--input).
#[test]
fn test_zipper_bam_mapped_input() {
    let temp_dir = TempDir::new().unwrap();
    let unmapped_bam = temp_dir.path().join("unmapped.bam");
    let mapped_bam = temp_dir.path().join("mapped.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let ref_path = create_test_reference(temp_dir.path());

    let unmapped_records = vec![
        RecordBuilder::new()
            .name("read1")
            .sequence("ACGTACGT")
            .qualities(&[30; 8])
            .unmapped(true)
            .tag("RX", "AACCGGTT")
            .tag("QX", "IIIIIIII")
            .build(),
        RecordBuilder::new()
            .name("read2")
            .sequence("TGCATGCA")
            .qualities(&[30; 8])
            .unmapped(true)
            .tag("RX", "GGTTCCAA")
            .tag("QX", "IIIIIIII")
            .build(),
    ];
    create_unmapped_bam(&unmapped_bam, &unmapped_records);

    let mapped_header = create_minimal_header("chr1", 10000);
    let mapped_records = vec![
        RecordBuilder::new()
            .name("read1")
            .sequence("ACGTACGT")
            .qualities(&[30; 8])
            .reference_sequence_id(0)
            .alignment_start(100)
            .mapping_quality(60)
            .cigar("8M")
            .build(),
        RecordBuilder::new()
            .name("read2")
            .sequence("TGCATGCA")
            .qualities(&[30; 8])
            .reference_sequence_id(0)
            .alignment_start(200)
            .mapping_quality(60)
            .cigar("8M")
            .build(),
    ];
    create_mapped_bam(&mapped_bam, &mapped_header, &mapped_records);

    let output = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "zipper",
            "--input",
            mapped_bam.to_str().unwrap(),
            "--unmapped",
            unmapped_bam.to_str().unwrap(),
            "--reference",
            ref_path.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--compression-level",
            "1",
        ])
        .output()
        .expect("Failed to run zipper command");

    assert!(
        output.status.success(),
        "Zipper command failed with BAM input: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output_bam.exists(), "Output BAM not created");

    // Verify output records have UMI tags transferred
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let header = reader.read_header().unwrap();
    let records: Vec<RecordBuf> = reader.record_bufs(&header).map(|r| r.unwrap()).collect();
    assert_eq!(records.len(), 2, "Should have 2 records in output");

    let rx_tag = Tag::from(SamTag::RX);
    for record in &records {
        assert!(record.data().get(&rx_tag).is_some(), "Output record should have RX tag");
    }

    // Verify warning about BAM input was emitted
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("BAM input detected"), "Should warn about BAM input. stderr: {stderr}");
}

/// Test zipper accepts BAM piped to stdin (auto-detected via BGZF magic bytes).
///
/// Common in nf-core / Galaxy pipelines that produce BAM from intermediate steps
/// (e.g. `bwameth.py | samtools view -b | fgumi zipper`). Without auto-detection,
/// fgumi treats stdin as SAM text and crashes with a confusing
/// "invalid flags / lexical parse error" error from the SAM parser misreading
/// BGZF binary bytes.
#[test]
fn test_zipper_bam_stdin_input() {
    let temp_dir = TempDir::new().unwrap();
    let unmapped_bam = temp_dir.path().join("unmapped.bam");
    let mapped_bam = temp_dir.path().join("mapped.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let ref_path = create_test_reference(temp_dir.path());

    let unmapped_records = vec![
        RecordBuilder::new()
            .name("read1")
            .sequence("ACGTACGT")
            .qualities(&[30; 8])
            .unmapped(true)
            .tag("RX", "AACCGGTT")
            .tag("QX", "IIIIIIII")
            .build(),
        RecordBuilder::new()
            .name("read2")
            .sequence("TGCATGCA")
            .qualities(&[30; 8])
            .unmapped(true)
            .tag("RX", "GGTTCCAA")
            .tag("QX", "IIIIIIII")
            .build(),
    ];
    create_unmapped_bam(&unmapped_bam, &unmapped_records);

    let mapped_header = create_minimal_header("chr1", 10000);
    let mapped_records = vec![
        RecordBuilder::new()
            .name("read1")
            .sequence("ACGTACGT")
            .qualities(&[30; 8])
            .reference_sequence_id(0)
            .alignment_start(100)
            .mapping_quality(60)
            .cigar("8M")
            .build(),
        RecordBuilder::new()
            .name("read2")
            .sequence("TGCATGCA")
            .qualities(&[30; 8])
            .reference_sequence_id(0)
            .alignment_start(200)
            .mapping_quality(60)
            .cigar("8M")
            .build(),
    ];
    create_mapped_bam(&mapped_bam, &mapped_header, &mapped_records);

    // Pipe the BAM bytes into the child's stdin; do NOT pass --input so that
    // the default ("-") triggers the stdin code path.
    let bam_bytes = fs::read(&mapped_bam).expect("read mapped BAM bytes");
    let mut child = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "zipper",
            "--unmapped",
            unmapped_bam.to_str().unwrap(),
            "--reference",
            ref_path.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--compression-level",
            "1",
        ])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("Failed to spawn zipper command");

    child
        .stdin
        .as_mut()
        .expect("Failed to open child stdin")
        .write_all(&bam_bytes)
        .expect("Failed to write BAM bytes to stdin");

    let output = child.wait_with_output().expect("Failed to wait for zipper");

    assert!(
        output.status.success(),
        "Zipper command failed with BAM on stdin: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output_bam.exists(), "Output BAM not created");

    // Verify output records have UMI tags transferred
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let header = reader.read_header().unwrap();
    let records: Vec<RecordBuf> = reader.record_bufs(&header).map(|r| r.unwrap()).collect();
    assert_eq!(records.len(), 2, "Should have 2 records in output");

    let rx_tag = Tag::from(SamTag::RX);
    for record in &records {
        assert!(record.data().get(&rx_tag).is_some(), "Output record should have RX tag");
    }
}

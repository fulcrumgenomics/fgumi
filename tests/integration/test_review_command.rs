//! End-to-end CLI tests for the review command.
//!
//! These tests run the actual `fgumi review` binary and validate:
//! 1. Missing required arguments produce a clear error
//! 2. Basic execution with minimal valid inputs

use fgumi_lib::sam::builder::RecordBuilder;
use noodles::bam;
use noodles::sam;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;
use tempfile::TempDir;

/// Create a reference FASTA + .fai + .dict that matches the test header (chr1, 10000bp).
fn create_test_reference(dir: &Path) -> PathBuf {
    let ref_path = dir.join("ref.fa");
    let fai_path = dir.join("ref.fa.fai");
    let dict_path = dir.join("ref.dict");

    let ref_seq = "ACGTACGT".repeat(1250); // 10000 bases
    let mut fasta = fs::File::create(&ref_path).unwrap();
    writeln!(fasta, ">chr1").unwrap();
    writeln!(fasta, "{ref_seq}").unwrap();
    fasta.flush().unwrap();

    // FASTA index (.fai): name\tlength\toffset\tlinebases\tlinewidth
    let fai_content = "chr1\t10000\t6\t10000\t10001\n";
    fs::write(&fai_path, fai_content).unwrap();

    // Sequence dictionary
    let mut dict = fs::File::create(&dict_path).unwrap();
    writeln!(dict, "@HD\tVN:1.6\tSO:unsorted").unwrap();
    writeln!(dict, "@SQ\tSN:chr1\tLN:10000").unwrap();
    dict.flush().unwrap();

    ref_path
}

/// Create a coordinate-sorted header for BAM indexing.
fn create_coordinate_sorted_header() -> sam::Header {
    use bstr::BString;
    use noodles::sam::header::record::value::map::Header as HeaderRecord;
    use noodles::sam::header::record::value::map::Map as HeaderRecordMap;
    use noodles::sam::header::record::value::map::header::tag::Tag as HeaderTag;
    use noodles::sam::header::record::value::{Map, map::ReferenceSequence};
    use std::num::NonZeroUsize;

    let HeaderTag::Other(sort_order_tag) = HeaderTag::from([b'S', b'O']) else { unreachable!() };

    let header_map = HeaderRecordMap::<HeaderRecord>::builder()
        .insert(sort_order_tag, "coordinate")
        .build()
        .expect("valid header map");

    let reference_sequence =
        Map::<ReferenceSequence>::new(NonZeroUsize::new(10000).expect("non-zero"));

    sam::Header::builder()
        .set_header(header_map)
        .add_reference_sequence(BString::from("chr1"), reference_sequence)
        .build()
}

/// Create a coordinate-sorted, indexed BAM file.
fn create_indexed_bam(path: &PathBuf, records: &[noodles::sam::alignment::record_buf::RecordBuf]) {
    let header = create_coordinate_sorted_header();

    // Write BAM and drop to flush
    {
        let mut writer =
            bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
        writer.write_header(&header).expect("Failed to write header");
        for record in records {
            writer.write_alignment_record(&header, record).expect("Failed to write record");
        }
    }

    // Create BAI index
    let index = bam::fs::index(path).expect("Failed to create BAM index");
    let index_path = path.with_extension("bam.bai");
    let mut index_writer =
        noodles::bam::bai::io::Writer::new(fs::File::create(&index_path).unwrap());
    index_writer.write_index(&index).unwrap();
}

/// Create a minimal VCF file with one variant.
fn create_test_vcf(path: &PathBuf) {
    let mut vcf = fs::File::create(path).unwrap();
    writeln!(vcf, "##fileformat=VCFv4.3").unwrap();
    writeln!(vcf, "##contig=<ID=chr1,length=10000>").unwrap();
    writeln!(vcf, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
    writeln!(vcf, "chr1\t100\t.\tA\tC\t.\tPASS\t.").unwrap();
    vcf.flush().unwrap();
}

/// Test that missing required arguments produce a non-zero exit code.
#[test]
fn test_review_missing_required_args() {
    // No arguments at all — should fail with clap error
    let output = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args(["review"])
        .output()
        .expect("Failed to run review command");

    assert!(!output.status.success(), "Review should fail without required args");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("required") || stderr.contains("Usage") || stderr.contains("error"),
        "Stderr should mention missing required arguments, got: {stderr}"
    );
}

/// Test that review fails clearly when required file inputs are missing.
#[test]
fn test_review_missing_input_files() {
    let temp_dir = TempDir::new().unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "review",
            "--input",
            "/nonexistent/variants.vcf",
            "--consensus-bam",
            "/nonexistent/consensus.bam",
            "--grouped-bam",
            "/nonexistent/grouped.bam",
            "--ref",
            "/nonexistent/ref.fa",
            "--output",
            temp_dir.path().join("output").to_str().unwrap(),
        ])
        .output()
        .expect("Failed to run review command");

    assert!(!output.status.success(), "Review should fail for nonexistent input files");
}

/// Test basic review execution with valid inputs.
///
/// Creates minimal VCF, consensus BAM, grouped BAM, and reference files.
/// Verifies the command completes and produces output files.
#[test]
fn test_review_basic_execution() {
    let temp_dir = TempDir::new().unwrap();
    let vcf_path = temp_dir.path().join("variants.vcf");
    let consensus_bam = temp_dir.path().join("consensus.bam");
    let grouped_bam = temp_dir.path().join("grouped.bam");
    let ref_path = create_test_reference(temp_dir.path());
    let output_prefix = temp_dir.path().join("review_out");

    // Create VCF with one variant at chr1:100
    create_test_vcf(&vcf_path);

    // Create a consensus BAM with a read spanning the variant position,
    // with MI tag for molecule identifier
    let consensus_records = vec![
        RecordBuilder::new()
            .name("cons1")
            .sequence("ACGTACGT")
            .qualities(&[30; 8])
            .reference_sequence_id(0)
            .alignment_start(97)
            .mapping_quality(60)
            .cigar("8M")
            .tag("MI", "1")
            .build(),
    ];
    create_indexed_bam(&consensus_bam, &consensus_records);

    // Create a grouped BAM with raw reads (MI tag matches consensus)
    let grouped_records = vec![
        RecordBuilder::new()
            .name("raw1")
            .sequence("ACGTACGT")
            .qualities(&[30; 8])
            .reference_sequence_id(0)
            .alignment_start(97)
            .mapping_quality(60)
            .cigar("8M")
            .tag("MI", "1")
            .build(),
        RecordBuilder::new()
            .name("raw2")
            .sequence("ACGTACGT")
            .qualities(&[30; 8])
            .reference_sequence_id(0)
            .alignment_start(97)
            .mapping_quality(60)
            .cigar("8M")
            .tag("MI", "1")
            .build(),
    ];
    create_indexed_bam(&grouped_bam, &grouped_records);

    let output = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args([
            "review",
            "--input",
            vcf_path.to_str().unwrap(),
            "--consensus-bam",
            consensus_bam.to_str().unwrap(),
            "--grouped-bam",
            grouped_bam.to_str().unwrap(),
            "--ref",
            ref_path.to_str().unwrap(),
            "--output",
            output_prefix.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to run review command");

    assert!(
        output.status.success(),
        "Review command failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Verify output files were created
    let consensus_out = output_prefix.with_extension("consensus.bam");
    let grouped_out = output_prefix.with_extension("grouped.bam");
    assert!(consensus_out.exists(), "Consensus output BAM should exist");
    assert!(grouped_out.exists(), "Grouped output BAM should exist");
}

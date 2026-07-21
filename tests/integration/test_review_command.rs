//! In-process tests for the review command.
//!
//! These tests invoke the Review command directly via in-process calls and validate:
//! 1. Missing required arguments produce a clap parse error
//! 2. Missing input files produce a validation error
//! 3. Basic execution with minimal valid inputs

use clap::Parser;
use fgumi_lib::commands::command::Command;
use fgumi_lib::commands::review::Review;
use fgumi_lib::sam::SamTag;
use fgumi_raw_bam::{RawRecord, SamBuilder};
use noodles::bam;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use rstest::rstest;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use tempfile::TempDir;

use crate::helpers::bam_generator::{
    create_coordinate_sorted_header, create_test_reference, to_record_buf,
};

/// Append `.bai` to the full BAM path, spelled out here (not via the production
/// helper under test) so the test owns an independent expectation of the samtools
/// sidecar convention.
fn appended_bai_path(bam_path: &Path) -> PathBuf {
    let mut index_os = bam_path.as_os_str().to_owned();
    index_os.push(".bai");
    PathBuf::from(index_os)
}

/// Create a coordinate-sorted, indexed BAM file and return the path of the BAI written.
///
/// The index path is derived here *independently* of the production helper under
/// test — by literally appending `.bai` to the full BAM path (see
/// [`appended_bai_path`]) — so this fixture is a faithful oracle. If
/// `bai_sidecar_path` / `Review::read_bam_index` regressed to the old
/// extension-replacing convention (`Path::with_extension("bam.bai")`), the index
/// would land at a path the command no longer looks in for any BAM not named
/// `*.bam`, and the test would fail. Using the production helper on both sides
/// would instead let a naming regression pass whenever writer and reader drift
/// together.
fn create_indexed_bam(path: &Path, records: &[RawRecord]) -> PathBuf {
    let header = create_coordinate_sorted_header("chr1", 10000);

    // Write BAM and drop to flush
    {
        let mut writer =
            bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
        writer.write_header(&header).expect("Failed to write header");
        for record in records {
            writer
                .write_alignment_record(&header, &to_record_buf(record))
                .expect("Failed to write record");
        }
    }

    let index_path = appended_bai_path(path);
    let index = bam::fs::index(path).expect("Failed to create BAM index");
    let mut index_writer =
        noodles::bam::bai::io::Writer::new(fs::File::create(&index_path).unwrap());
    index_writer.write_index(&index).expect("Failed to write BAM index");
    index_path
}

/// Create a minimal VCF file with one variant.
fn create_test_vcf(path: &Path) {
    let mut vcf = fs::File::create(path).unwrap();
    writeln!(vcf, "##fileformat=VCFv4.3").unwrap();
    writeln!(vcf, "##contig=<ID=chr1,length=10000>").unwrap();
    writeln!(vcf, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
    writeln!(vcf, "chr1\t100\t.\tA\tC\t.\tPASS\t.").unwrap();
    vcf.flush().unwrap();
}

/// Test that missing required arguments produce a clap parse error.
#[test]
fn test_review_missing_required_args() {
    // No arguments at all — should fail with clap error
    let err =
        Review::try_parse_from(["review"]).expect_err("Review should fail without required args");
    assert_eq!(
        err.kind(),
        clap::error::ErrorKind::MissingRequiredArgument,
        "Expected MissingRequiredArgument, got {:?}: {err}",
        err.kind()
    );
}

/// Test that review fails clearly when required file inputs are missing.
#[test]
fn test_review_missing_input_files() {
    let temp_dir = TempDir::new().unwrap();

    let cmd = Review::try_parse_from([
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
    .expect("failed to parse review args");
    let err =
        cmd.execute("fgumi review").expect_err("Review should fail for nonexistent input files");
    let err_msg = format!("{err:#}");
    assert!(
        err_msg.contains("does not exist"),
        "Validation error should mention missing file, got: {err_msg}"
    );
}

/// Test basic review execution with valid inputs.
///
/// Creates minimal VCF, consensus BAM, grouped BAM, and reference files.
/// Verifies the command completes and produces output files.
///
/// Parameterized over the input BAM naming so the index-lookup contract is pinned
/// for the cases that distinguish the samtools sidecar convention from the old
/// extension-replacing one: a plain `*.bam` (where the two agree), a non-`.bam`
/// suffix, and an extensionless name. The command reads its indexes through
/// `Review::read_bam_index`, so it only succeeds if it finds the BAI at the
/// appended path this fixture wrote.
#[rstest]
#[case::dot_bam("consensus.bam", "grouped.bam")]
#[case::non_bam_suffix("consensus.sorted", "grouped.sorted")]
#[case::extensionless("consensus", "grouped")]
fn test_review_basic_execution(#[case] consensus_name: &str, #[case] grouped_name: &str) {
    let temp_dir = TempDir::new().unwrap();
    let vcf_path = temp_dir.path().join("variants.vcf");
    let consensus_bam = temp_dir.path().join(consensus_name);
    let grouped_bam = temp_dir.path().join(grouped_name);
    let ref_path = create_test_reference(temp_dir.path());
    let output_prefix = temp_dir.path().join("review_out");

    // Create VCF with one variant at chr1:100
    create_test_vcf(&vcf_path);

    // Create a consensus BAM with a read spanning the variant position,
    // with MI tag for molecule identifier
    let consensus_records = vec![{
        let mut b = SamBuilder::new();
        b.read_name(b"cons1")
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .ref_id(0)
            .pos(96) // 0-based for 1-based pos 97
            .mapq(60)
            .cigar_ops(&[8 << 4]) // 8M
            .add_string_tag(SamTag::MI, b"1");
        b.build()
    }];
    let consensus_index = create_indexed_bam(&consensus_bam, &consensus_records);

    // Create a grouped BAM with raw reads (MI tag matches consensus)
    let grouped_records = vec![
        {
            let mut b = SamBuilder::new();
            b.read_name(b"raw1")
                .sequence(b"ACGTACGT")
                .qualities(&[30; 8])
                .ref_id(0)
                .pos(96)
                .mapq(60)
                .cigar_ops(&[8 << 4]) // 8M
                .add_string_tag(SamTag::MI, b"1");
            b.build()
        },
        {
            let mut b = SamBuilder::new();
            b.read_name(b"raw2")
                .sequence(b"ACGTACGT")
                .qualities(&[30; 8])
                .ref_id(0)
                .pos(96)
                .mapq(60)
                .cigar_ops(&[8 << 4]) // 8M
                .add_string_tag(SamTag::MI, b"1");
            b.build()
        },
    ];
    let grouped_index = create_indexed_bam(&grouped_bam, &grouped_records);

    // The fixture wrote each index at the appended sidecar path (`<bam>.bai`). Assert
    // it is there, and that nothing was written at the extension-replaced path the old
    // buggy convention produced — so a regression that reads from the wrong path cannot
    // be masked by an index happening to sit there.
    for (bam, index) in [(&consensus_bam, &consensus_index), (&grouped_bam, &grouped_index)] {
        assert!(index.exists(), "fixture must write the BAI sidecar at {}", index.display());
        let extension_replaced = bam.with_extension("bam.bai");
        if extension_replaced != *index {
            assert!(
                !extension_replaced.exists(),
                "no index should exist at the extension-replaced path {}",
                extension_replaced.display(),
            );
        }
    }

    let cmd = Review::try_parse_from([
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
    .expect("failed to parse review args");
    cmd.execute("fgumi review").unwrap_or_else(|e| panic!("Review command failed: {e:#}"));

    // Verify output files were created
    let consensus_out = output_prefix.with_extension("consensus.bam");
    let grouped_out = output_prefix.with_extension("grouped.bam");
    assert!(consensus_out.exists(), "Consensus output BAM should exist");
    assert!(grouped_out.exists(), "Grouped output BAM should exist");
}

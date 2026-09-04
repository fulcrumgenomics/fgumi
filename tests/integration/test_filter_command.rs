//! End-to-end tests for the filter command.
//!
//! These tests invoke `Filter::execute()` in-process and validate:
//! 1. Basic consensus read filtering
//! 2. Rejected reads output
//! 3. Statistics output

use clap::Parser;
use fgumi_lib::commands::command::Command;
use fgumi_lib::commands::filter::Filter;
use fgumi_lib::sam::SamTag;
use fgumi_raw_bam::{RawRecord, SamBuilder as RawSamBuilder, flags};
use noodles::bam;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use rstest::rstest;
use std::fs;
use std::path::Path;
use tempfile::TempDir;

use crate::helpers::bam_generator::{
    create_coordinate_sorted_header, create_minimal_header, create_test_reference,
};

/// Convert a raw BAM record to a `RecordBuf` using the given header.
///
/// Used to bridge the raw-byte builder API with the noodles BAM writer for test
/// file creation. Mapped records reference the header's sequence dictionary by
/// index, so callers must pass the same header they use to write the BAM.
fn to_record_buf(
    raw: &RawRecord,
    header: &noodles::sam::Header,
) -> noodles::sam::alignment::RecordBuf {
    fgumi_raw_bam::raw_record_to_record_buf(raw, header)
        .expect("raw_record_to_record_buf should succeed in test")
}

/// Create a consensus BAM with the given header and consensus-tagged records.
fn create_consensus_bam_with_header(
    path: &Path,
    header: &noodles::sam::Header,
    records: Vec<RawRecord>,
) {
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(header).expect("Failed to write header");

    for record in records {
        writer
            .write_alignment_record(header, &to_record_buf(&record, header))
            .expect("Failed to write record");
    }
    writer.try_finish().expect("Failed to finish BAM");
}

/// Create a consensus BAM with a query-grouped (template-coordinate) header.
fn create_consensus_bam(path: &Path, records: Vec<RawRecord>) {
    let header = create_minimal_header("chr1", 10000);
    create_consensus_bam_with_header(path, &header, records);
}

/// Two mapped-consensus records with CD/CE per-base tags that pass `filter`
/// (`cons1` depth 10, `cons2` depth 5). Shared by the query-grouped writer and
/// the coordinate-sorted guard test.
fn filter_consensus_records() -> Vec<RawRecord> {
    let r1 = {
        let mut b = RawSamBuilder::new();
        b.read_name(b"cons1")
            .ref_id(0)
            .pos(99)
            .mapq(60)
            .cigar_ops(&[8 << 4]) // 8M
            .sequence(b"ACGTACGT")
            .qualities(&[35; 8]);
        b.add_int_tag(SamTag::CD, 10).add_float_tag(SamTag::CE, 0.0_f32);
        b.add_array_u16(SamTag::CD_BASES, &[10; 8]).add_array_u16(SamTag::CE_BASES, &[0; 8]);
        b.build()
    };
    let r2 = {
        let mut b = RawSamBuilder::new();
        b.read_name(b"cons2")
            .ref_id(0)
            .pos(199)
            .mapq(60)
            .cigar_ops(&[8 << 4]) // 8M
            .sequence(b"ACGTACGT")
            .qualities(&[35; 8]);
        b.add_int_tag(SamTag::CD, 5).add_float_tag(SamTag::CE, 0.0_f32);
        b.add_array_u16(SamTag::CD_BASES, &[5; 8]).add_array_u16(SamTag::CE_BASES, &[0; 8]);
        b.build()
    };
    vec![r1, r2]
}

/// Write a small mapped-consensus BAM (two reads with CD/CE per-base tags that
/// pass `filter`). Shared by `test_filter_command` and the stdin parity test.
pub(crate) fn write_filter_consensus_bam(path: &Path) {
    create_consensus_bam(path, filter_consensus_records());
}

/// Write the same passing consensus reads, but with a header carrying only `@SQ`
/// (no `@HD` line) -- mirrors the header-less inputs `filter` must normalize.
fn write_headerless_consensus_bam(path: &Path) {
    use noodles::sam::header::record::value::Map;
    use noodles::sam::header::record::value::map::ReferenceSequence;

    let header = noodles::sam::Header::builder()
        .add_reference_sequence(
            "chr1",
            Map::<ReferenceSequence>::new(
                std::num::NonZero::new(10000).expect("non-zero reference length"),
            ),
        )
        .build();
    assert!(header.header().is_none(), "precondition: input must lack @HD");

    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");
    for record in filter_consensus_records() {
        writer
            .write_alignment_record(&header, &to_record_buf(&record, &header))
            .expect("Failed to write record");
    }
    writer.try_finish().expect("Failed to finish BAM");
}

/// Test basic filter command with passing reads.
#[test]
fn test_filter_command_basic() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let ref_path = create_test_reference(temp_dir.path());

    write_filter_consensus_bam(&input_bam);

    let cmd = Filter::try_parse_from([
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
    .expect("failed to parse filter args");
    cmd.execute("fgumi filter").expect("Filter command failed");
    assert!(output_bam.exists(), "Output BAM not created");

    // Verify output has records (reads may have masked bases but should still be present)
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let count = reader.records().count();
    assert_eq!(count, 2, "Both reads should be present in output");
}

/// FILT3-02: filter must reject coordinate-sorted (non-query-grouped) input,
/// matching fgbio's `Bams.requireQueryGrouped`. On coordinate-sorted input a
/// template's mates scatter, every read degrades to its own "template", and the
/// both-primaries-pass logic silently corrupts the result — so we hard-fail.
#[test]
fn test_filter_rejects_coordinate_sorted_input() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let ref_path = create_test_reference(temp_dir.path());

    // Same consensus records as write_filter_consensus_bam, but a coordinate header.
    let header = create_coordinate_sorted_header("chr1", 10000);
    create_consensus_bam_with_header(&input_bam, &header, filter_consensus_records());

    let cmd = Filter::try_parse_from([
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
    .expect("failed to parse filter args");

    let err = cmd.execute("fgumi filter").expect_err("must reject coordinate-sorted input");
    let msg = err.to_string();
    assert!(msg.contains("queryname sorted or query grouped"), "unexpected error message: {msg}");
}

/// FILT3-02: filter must reject header-less input. A header-less BAM synthesizes
/// `@HD VN:1.6 SO:unsorted` (via `ensure_hd_record`), which is neither queryname
/// sorted nor query grouped, so `require_query_grouped` rejects it — matching
/// fgbio's `Bams.requireQueryGrouped`. The `@HD` synthesis itself is still exercised
/// end-to-end by `correct`/`review` (which have no query-grouped guard).
#[test]
fn test_filter_rejects_headerless_input() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let ref_path = create_test_reference(temp_dir.path());

    write_headerless_consensus_bam(&input_bam);

    let cmd = Filter::try_parse_from([
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
    .expect("failed to parse filter args");

    let err = cmd.execute("fgumi filter").expect_err("must reject header-less input");
    let msg = err.to_string();
    assert!(msg.contains("queryname sorted or query grouped"), "unexpected error message: {msg}");
}

/// Test filter command with reads that fail due to low per-base depth.
///
/// The filter masks bases where the per-base depth (cd tag) is below --min-reads.
/// If enough bases are masked, the no-call fraction exceeds the threshold and the
/// read is rejected.
#[test]
fn test_filter_command_rejects_low_depth() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let rejects_bam = temp_dir.path().join("rejects.bam");
    let ref_path = create_test_reference(temp_dir.path());

    // Good read: per-base depth 10 (all above min-reads=3), no bases masked
    let good = {
        let mut b = RawSamBuilder::new();
        b.read_name(b"good")
            .ref_id(0)
            .pos(99)
            .mapq(60)
            .cigar_ops(&[8 << 4]) // 8M
            .sequence(b"ACGTACGT")
            .qualities(&[35; 8]);
        b.add_int_tag(SamTag::CD, 10).add_float_tag(SamTag::CE, 0.0_f32);
        b.add_array_u16(SamTag::CD_BASES, &[10; 8]).add_array_u16(SamTag::CE_BASES, &[0; 8]);
        b.build()
    };

    // Low-depth read: per-base depth 1 (all below min-reads=3), all bases masked
    let low_depth = {
        let mut b = RawSamBuilder::new();
        b.read_name(b"low_depth")
            .ref_id(0)
            .pos(199)
            .mapq(60)
            .cigar_ops(&[8 << 4]) // 8M
            .sequence(b"ACGTACGT")
            .qualities(&[35; 8]);
        b.add_int_tag(SamTag::CD, 1).add_float_tag(SamTag::CE, 0.0_f32);
        b.add_array_u16(SamTag::CD_BASES, &[1; 8]).add_array_u16(SamTag::CE_BASES, &[0; 8]);
        b.build()
    };

    create_consensus_bam(&input_bam, vec![good, low_depth]);

    let cmd = Filter::try_parse_from([
        "filter",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--ref",
        ref_path.to_str().unwrap(),
        "--min-reads",
        "3",
        "--rejects",
        rejects_bam.to_str().unwrap(),
        "--compression-level",
        "1",
    ])
    .expect("failed to parse filter args");
    cmd.execute("fgumi filter").expect("Filter command with rejects failed");
    assert!(output_bam.exists(), "Output BAM not created");
    assert!(rejects_bam.exists(), "Rejects BAM not created");

    // Verify the good read passed and the low-depth read was rejected
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let output_count = reader.records().count();
    assert_eq!(output_count, 1, "Only the good read should pass filtering");

    let mut reject_reader = bam::io::Reader::new(fs::File::open(&rejects_bam).unwrap());
    let _header = reject_reader.read_header().unwrap();
    let reject_count = reject_reader.records().count();
    assert_eq!(reject_count, 1, "The low-depth read should be rejected");
}

/// Test filter command with statistics output.
#[test]
fn test_filter_command_with_stats() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let stats_path = temp_dir.path().join("stats.txt");
    let ref_path = create_test_reference(temp_dir.path());

    let record = {
        let mut b = RawSamBuilder::new();
        b.read_name(b"cons1")
            .ref_id(0)
            .pos(99)
            .mapq(60)
            .cigar_ops(&[8 << 4]) // 8M
            .sequence(b"ACGTACGT")
            .qualities(&[35; 8]);
        b.add_int_tag(SamTag::CD, 10)
            .add_int_tag(SamTag::CM, 8)
            .add_float_tag(SamTag::CE, 0.005_f32);
        b.build()
    };

    create_consensus_bam(&input_bam, vec![record]);

    let cmd = Filter::try_parse_from([
        "filter",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--ref",
        ref_path.to_str().unwrap(),
        "--min-reads",
        "1",
        "--stats",
        stats_path.to_str().unwrap(),
        "--compression-level",
        "1",
    ])
    .expect("failed to parse filter args");
    cmd.execute("fgumi filter").expect("Filter command with stats failed");
    assert!(stats_path.exists(), "Stats file not created");
    let content = fs::read_to_string(&stats_path).expect("Failed to read stats file");
    assert!(!content.trim().is_empty(), "Stats file should not be empty");
    assert!(content.contains("total_reads"), "Stats should contain total_reads");
}

/// Test filter command without --ref on unmapped consensus reads.
///
/// When all reads are unmapped (e.g. pre-alignment consensus pipeline), the reference
/// is never consulted for tag regeneration, so --ref should be optional.
#[test]
fn test_filter_command_no_ref_unmapped_reads() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    // Create unmapped consensus reads with good per-base tags
    let r1 = {
        let mut b = RawSamBuilder::new();
        b.read_name(b"cons1").flags(flags::UNMAPPED).sequence(b"ACGTACGT").qualities(&[35; 8]);
        b.add_int_tag(SamTag::CD, 10).add_float_tag(SamTag::CE, 0.0_f32);
        b.add_array_u16(SamTag::CD_BASES, &[10; 8]).add_array_u16(SamTag::CE_BASES, &[0; 8]);
        b.build()
    };

    let r2 = {
        let mut b = RawSamBuilder::new();
        b.read_name(b"cons2").flags(flags::UNMAPPED).sequence(b"ACGTACGT").qualities(&[35; 8]);
        b.add_int_tag(SamTag::CD, 5).add_float_tag(SamTag::CE, 0.0_f32);
        b.add_array_u16(SamTag::CD_BASES, &[5; 8]).add_array_u16(SamTag::CE_BASES, &[0; 8]);
        b.build()
    };

    create_consensus_bam(&input_bam, vec![r1, r2]);

    // Run filter WITHOUT --ref
    let cmd = Filter::try_parse_from([
        "filter",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--min-reads",
        "1",
        "--max-no-call-fraction",
        "1.0",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse filter args");
    cmd.execute("fgumi filter").expect("Filter command without --ref failed");
    assert!(output_bam.exists(), "Output BAM not created");

    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let count = reader.records().count();
    assert_eq!(count, 2, "Both unmapped reads should pass filtering without --ref");
}

/// Test filter command without --ref, with --rejects output for unmapped reads.
///
/// One good unmapped read and one low-depth unmapped read. The good read should
/// pass and the low-depth read should be rejected.
#[test]
fn test_filter_command_no_ref_with_rejects() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let rejects_bam = temp_dir.path().join("rejects.bam");

    // Good unmapped read: per-base depth 10 (above min-reads=3)
    let good = {
        let mut b = RawSamBuilder::new();
        b.read_name(b"good").flags(flags::UNMAPPED).sequence(b"ACGTACGT").qualities(&[35; 8]);
        b.add_int_tag(SamTag::CD, 10).add_float_tag(SamTag::CE, 0.0_f32);
        b.add_array_u16(SamTag::CD_BASES, &[10; 8]).add_array_u16(SamTag::CE_BASES, &[0; 8]);
        b.build()
    };

    // Low-depth unmapped read: per-base depth 1 (below min-reads=3), all bases masked
    let low_depth = {
        let mut b = RawSamBuilder::new();
        b.read_name(b"low_depth").flags(flags::UNMAPPED).sequence(b"ACGTACGT").qualities(&[35; 8]);
        b.add_int_tag(SamTag::CD, 1).add_float_tag(SamTag::CE, 0.0_f32);
        b.add_array_u16(SamTag::CD_BASES, &[1; 8]).add_array_u16(SamTag::CE_BASES, &[0; 8]);
        b.build()
    };

    create_consensus_bam(&input_bam, vec![good, low_depth]);

    // Run filter WITHOUT --ref, with --rejects
    let cmd = Filter::try_parse_from([
        "filter",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--min-reads",
        "3",
        "--rejects",
        rejects_bam.to_str().unwrap(),
        "--compression-level",
        "1",
    ])
    .expect("failed to parse filter args");
    cmd.execute("fgumi filter").expect("Filter command without --ref with rejects failed");
    assert!(output_bam.exists(), "Output BAM not created");
    assert!(rejects_bam.exists(), "Rejects BAM not created");

    // Verify the good read passed
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let output_count = reader.records().count();
    assert_eq!(output_count, 1, "Only the good read should pass filtering");

    // Verify the low-depth read was rejected
    let mut reject_reader = bam::io::Reader::new(fs::File::open(&rejects_bam).unwrap());
    let _header = reject_reader.read_header().unwrap();
    let reject_count = reject_reader.records().count();
    assert_eq!(reject_count, 1, "The low-depth read should be rejected");
}

/// Test that filter command without --ref fails when given mapped reads.
///
/// Mapped reads require the reference for NM/UQ/MD tag regeneration after base masking.
/// The command should fail with a clear error rather than emit stale tags.
#[test]
fn test_filter_command_no_ref_mapped_reads_fails() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    // Create a mapped consensus read (has ref_id, alignment_start, cigar)
    let mapped = {
        let mut b = RawSamBuilder::new();
        b.read_name(b"mapped_read")
            .ref_id(0)
            .pos(99)
            .mapq(60)
            .cigar_ops(&[8 << 4]) // 8M
            .sequence(b"ACGTACGT")
            .qualities(&[35; 8]);
        b.add_int_tag(SamTag::CD, 10).add_float_tag(SamTag::CE, 0.0_f32);
        b.add_array_u16(SamTag::CD_BASES, &[10; 8]).add_array_u16(SamTag::CE_BASES, &[0; 8]);
        b.build()
    };

    create_consensus_bam(&input_bam, vec![mapped]);

    // Run filter WITHOUT --ref on mapped reads — should fail
    let cmd = Filter::try_parse_from([
        "filter",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--min-reads",
        "1",
        "--max-no-call-fraction",
        "1.0",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse filter args");
    let result = cmd.execute("fgumi filter");
    assert!(result.is_err(), "Filter command should fail for mapped reads without --ref");
    let err_msg = format!("{:#}", result.unwrap_err());
    assert!(
        err_msg.contains("--ref is required"),
        "Error should mention --ref requirement, got: {err_msg}"
    );
}

//////////////////////////////////////////////////////////////////////////////
// R3.1 chain-cutover parity tests: `--threads N` (declarative chain builder)
// vs. no-`--threads` (legacy unified-pipeline oracle).
//////////////////////////////////////////////////////////////////////////////

/// Run `filter` on `input` writing `output`, with `extra` args appended
/// (e.g. `--threads`, `--rejects`, `--stats`, `--filter-by-template`). Always
/// passes `--ref` (mapped fixtures need it for tag regeneration) and
/// `--min-reads 3` -- deliberately *not* overriding `--max-no-call-fraction`
/// (default 0.2), so a fully-masked low-depth read is actually rejected rather
/// than passed through. Asserts the run succeeds.
fn filter_run(input: &Path, output: &Path, ref_path: &Path, extra: &[&str]) {
    let mut args = vec![
        "filter",
        "--input",
        input.to_str().unwrap(),
        "--output",
        output.to_str().unwrap(),
        "--ref",
        ref_path.to_str().unwrap(),
        "--min-reads",
        "3",
        "--compression-level",
        "1",
    ];
    args.extend_from_slice(extra);
    Filter::try_parse_from(args)
        .expect("failed to parse filter args")
        .execute("fgumi filter")
        .unwrap_or_else(|e| panic!("filter run failed: {e:#}"));
}

/// Read a BAM's records back as decoded `RecordBuf`s, for record-for-record
/// comparison (order-sensitive, catches drops/dupes/reorders that a bare count
/// would miss).
fn read_filter_records(path: &Path) -> Vec<noodles::sam::alignment::RecordBuf> {
    let mut reader = bam::io::Reader::new(fs::File::open(path).unwrap());
    let header = reader.read_header().unwrap();
    reader.record_bufs(&header).collect::<std::io::Result<Vec<_>>>().expect("read filter records")
}

/// Build `n` two-mate templates (R1/R2 sharing a read name, both mapped
/// primaries). Even-indexed templates have both mates at passing depth
/// (`cD 10 >= --min-reads 3`); odd-indexed templates have a passing R1 but a
/// failing (`cD 1`, fully masked) R2. This is what makes both filter step
/// factories do genuinely non-trivial, non-vacuous work:
/// - filter-by-template mode drops the *whole* odd template (fgbio's
///   "all primaries must pass" rule), so both mates of an odd template are
///   rejected even though R1 alone would have passed;
/// - single-read mode keeps every R1 and drops only the failing R2s.
fn build_mixed_depth_templates(n: usize) -> Vec<RawRecord> {
    let mut records = Vec::with_capacity(n * 2);
    for i in 0..i32::try_from(n).expect("n fits in i32") {
        let name = format!("t{i}");
        let pos = 100 + i * 20;
        let r2_depth: u16 = if i % 2 == 0 { 10 } else { 1 };

        let mut r1 = RawSamBuilder::new();
        r1.read_name(name.as_bytes())
            .flags(flags::PAIRED | flags::FIRST_SEGMENT)
            .ref_id(0)
            .pos(pos)
            .mapq(60)
            .cigar_ops(&[8 << 4]) // 8M
            .sequence(b"ACGTACGT")
            .qualities(&[35; 8]);
        r1.add_int_tag(SamTag::CD, 10).add_float_tag(SamTag::CE, 0.0_f32);
        r1.add_array_u16(SamTag::CD_BASES, &[10; 8]).add_array_u16(SamTag::CE_BASES, &[0; 8]);
        records.push(r1.build());

        let mut r2 = RawSamBuilder::new();
        r2.read_name(name.as_bytes())
            .flags(flags::PAIRED | flags::LAST_SEGMENT)
            .ref_id(0)
            .pos(pos + 8)
            .mapq(60)
            .cigar_ops(&[8 << 4]) // 8M
            .sequence(b"ACGTACGT")
            .qualities(&[35; 8]);
        r2.add_int_tag(SamTag::CD, i32::from(r2_depth)).add_float_tag(SamTag::CE, 0.0_f32);
        r2.add_array_u16(SamTag::CD_BASES, &[r2_depth; 8]).add_array_u16(SamTag::CE_BASES, &[0; 8]);
        records.push(r2.build());
    }
    records
}

/// The chain (`--threads N`) path produces output identical to the non-chain
/// (no-`--threads`) path, for both `--threads 1` (the minimal chain engine)
/// and `--threads 4` (genuinely parallel), and for both `filter-by-template`
/// step factories (template mode groups via `GroupByQueryname` before
/// filtering; single-read mode filters each record independently -- distinct
/// step factories in `chains/commands/filter.rs`). `build_mixed_depth_templates`
/// guarantees a non-vacuous mix of passing and rejected reads in every case.
#[rstest]
#[case::threads1_template(1, true)]
#[case::threads4_template(4, true)]
#[case::threads1_single_read(1, false)]
#[case::threads4_single_read(4, false)]
fn test_filter_chain_matches_single_threaded(
    #[case] threads: usize,
    #[case] filter_by_template: bool,
) {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let ref_path = create_test_reference(temp_dir.path());
    create_consensus_bam(&input_bam, build_mixed_depth_templates(8));

    let template_flag = if filter_by_template { "true" } else { "false" };
    let threads_str = threads.to_string();

    let oracle_out = temp_dir.path().join("oracle.bam");
    filter_run(&input_bam, &oracle_out, &ref_path, &["--filter-by-template", template_flag]);

    let chain_out = temp_dir.path().join("chain.bam");
    filter_run(
        &input_bam,
        &chain_out,
        &ref_path,
        &["--filter-by-template", template_flag, "--threads", &threads_str],
    );

    // Compare the COMPLETE normalized header (not just `@HD`): `read_bam_output`
    // normalizes the `@PG` command-line field that legitimately differs by
    // `--threads`, so a chain regression that changed `@RG`/`@SQ`/`@CO` or any
    // other header record — not only sort order — is caught too.
    let (oracle_header, expected) = crate::helpers::read_bam_output(&oracle_out);
    let (chain_header, actual) = crate::helpers::read_bam_output(&chain_out);
    // Non-vacuous: filtering must actually drop records, or a pass-through
    // regression on BOTH paths would satisfy `actual == expected` silently.
    // `build_mixed_depth_templates(8)` = 16 records; template mode drops the 4
    // odd templates whole (keeps 8), single-read mode drops only the 4 failing
    // R2s (keeps 12).
    let expected_kept = if filter_by_template { 8 } else { 12 };
    assert_eq!(
        expected.len(),
        expected_kept,
        "oracle must keep exactly {expected_kept} of 16 records \
         (filter_by_template={filter_by_template})"
    );
    assert_eq!(
        actual, expected,
        "chain (threads={threads}, filter_by_template={filter_by_template}) output must match \
         the non-chain path record-for-record"
    );
    // Guard against a vacuous header comparison: if a regression dropped `@HD` on
    // both paths, the header equality could still pass on two equally-broken
    // headers, so assert the oracle actually declares `@HD` first.
    assert!(oracle_header.header().is_some(), "oracle output must declare @HD");
    assert_eq!(
        chain_header, oracle_header,
        "chain (threads={threads}, filter_by_template={filter_by_template}) output header must \
         match the non-chain path (complete normalized header, not just @HD)"
    );
}

/// Header carrying two `@RG` lines with distinct `LB` values, for
/// `test_filter_chain_matches_single_threaded_with_rg_and_cb_variation`. The
/// `source_group_key_config` perf fix routes the chain filter path's first
/// stage to a `name_hash_only` `GroupKeyConfig`, which never resolves a
/// per-record library index or walks the CIGAR for position during decode —
/// unlike the filter-by-template oracle's own decode config
/// (`Filter::build_filter_pipeline_config`'s `new_raw_no_cell`), which
/// resolves the library index and the position but, like `name_hash_only`,
/// never extracts `CB` either way. So RG (library index) and position are the
/// fields that genuinely differ between the two decode configs, and are what
/// this parity test needs to vary to be non-vacuous; the `CB` variation below
/// is incidental fixture realism, not a discriminating input — the oracle
/// never reads `CB` regardless of this change, so it cannot expose a
/// regression here.
fn create_consensus_header_with_read_groups(
    ref_name: &str,
    ref_len: usize,
) -> noodles::sam::Header {
    use bstr::BString;
    use noodles::sam::header::record::value::Map;
    use noodles::sam::header::record::value::map::Map as HeaderRecordMap;
    use noodles::sam::header::record::value::map::header::tag::Tag as HeaderTag;
    use noodles::sam::header::record::value::map::read_group::tag as rg_tag;
    use noodles::sam::header::record::value::map::{
        Header as HeaderRecord, ReadGroup, ReferenceSequence,
    };
    use std::num::NonZeroUsize;

    let mut header_builder = HeaderRecordMap::<HeaderRecord>::builder();
    for &(tag_bytes, value) in
        &[(*b"SO", "unsorted"), (*b"GO", "query"), (*b"SS", "template-coordinate")]
    {
        let HeaderTag::Other(tag) = HeaderTag::from(tag_bytes) else { unreachable!() };
        header_builder = header_builder.insert(tag, value);
    }
    let header_map = header_builder.build().expect("valid header map");

    let reference_sequence = Map::<ReferenceSequence>::new(
        NonZeroUsize::new(ref_len).expect("reference length must be non-zero"),
    );

    let rg_a = Map::<ReadGroup>::builder()
        .insert(rg_tag::LIBRARY, String::from("libA"))
        .build()
        .expect("building read group RG1 should succeed");
    let rg_b = Map::<ReadGroup>::builder()
        .insert(rg_tag::LIBRARY, String::from("libB"))
        .build()
        .expect("building read group RG2 should succeed");

    noodles::sam::Header::builder()
        .set_header(header_map)
        .add_reference_sequence(BString::from(ref_name), reference_sequence)
        .add_read_group(BString::from("RG1"), rg_a)
        .add_read_group(BString::from("RG2"), rg_b)
        .build()
}

/// Like [`build_mixed_depth_templates`], but each mate additionally carries an
/// `RG` tag (alternating `RG1`/`RG2`, so both `@RG` lines in
/// [`create_consensus_header_with_read_groups`] are exercised) and a `CB` tag
/// that varies per template (`CB{i}`). Same mixed pass/fail depth pattern, so
/// filtering still does non-trivial work.
fn build_mixed_depth_templates_with_rg_and_cb(n: usize) -> Vec<RawRecord> {
    let mut records = Vec::with_capacity(n * 2);
    for i in 0..i32::try_from(n).expect("n fits in i32") {
        let name = format!("t{i}");
        let pos = 100 + i * 20;
        let r2_depth: u16 = if i % 2 == 0 { 10 } else { 1 };
        let rg_id: &[u8] = if i % 2 == 0 { b"RG1" } else { b"RG2" };
        let cb = format!("CB{i}");

        let mut r1 = RawSamBuilder::new();
        r1.read_name(name.as_bytes())
            .flags(flags::PAIRED | flags::FIRST_SEGMENT)
            .ref_id(0)
            .pos(pos)
            .mapq(60)
            .cigar_ops(&[8 << 4]) // 8M
            .sequence(b"ACGTACGT")
            .qualities(&[35; 8]);
        r1.add_int_tag(SamTag::CD, 10).add_float_tag(SamTag::CE, 0.0_f32);
        r1.add_array_u16(SamTag::CD_BASES, &[10; 8]).add_array_u16(SamTag::CE_BASES, &[0; 8]);
        r1.add_string_tag(SamTag::RG, rg_id).add_string_tag(SamTag::CB, cb.as_bytes());
        records.push(r1.build());

        let mut r2 = RawSamBuilder::new();
        r2.read_name(name.as_bytes())
            .flags(flags::PAIRED | flags::LAST_SEGMENT)
            .ref_id(0)
            .pos(pos + 8)
            .mapq(60)
            .cigar_ops(&[8 << 4]) // 8M
            .sequence(b"ACGTACGT")
            .qualities(&[35; 8]);
        r2.add_int_tag(SamTag::CD, i32::from(r2_depth)).add_float_tag(SamTag::CE, 0.0_f32);
        r2.add_array_u16(SamTag::CD_BASES, &[r2_depth; 8]).add_array_u16(SamTag::CE_BASES, &[0; 8]);
        r2.add_string_tag(SamTag::RG, rg_id).add_string_tag(SamTag::CB, cb.as_bytes());
        records.push(r2.build());
    }
    records
}

/// The chain (`--threads N`) path matches the non-chain oracle record-for-record
/// even when the discarded position/RG portion of the group key genuinely
/// varies per record (distinct `@RG`/`LB` per template and distinct positions
/// — the fields the oracle's decode config resolves but `name_hash_only`
/// does not, per `create_consensus_header_with_read_groups`'s doc comment).
/// `test_filter_chain_matches_single_threaded` already covers the
/// `threads`/`filter_by_template` matrix on RG-free, single-position records;
/// this test's unique contribution is exercising the exact fields
/// `source_group_key_config`'s `name_hash_only` routing stops computing for
/// the chain filter path, confirming the skip is genuinely invisible in the
/// output. Records also carry a per-template `CB` tag for fixture realism,
/// but it is not discriminating: the oracle's own decode config never
/// extracts `CB` either, so `CB` variation cannot expose a regression here.
#[test]
fn test_filter_chain_matches_single_threaded_with_rg_and_cb_variation() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let ref_path = create_test_reference(temp_dir.path());
    let header = create_consensus_header_with_read_groups("chr1", 10000);
    create_consensus_bam_with_header(
        &input_bam,
        &header,
        build_mixed_depth_templates_with_rg_and_cb(8),
    );

    let oracle_out = temp_dir.path().join("oracle.bam");
    filter_run(&input_bam, &oracle_out, &ref_path, &["--filter-by-template", "true"]);

    let chain_out = temp_dir.path().join("chain.bam");
    filter_run(
        &input_bam,
        &chain_out,
        &ref_path,
        &["--filter-by-template", "true", "--threads", "4"],
    );

    let (oracle_header, expected) = crate::helpers::read_bam_output(&oracle_out);
    let (chain_header, actual) = crate::helpers::read_bam_output(&chain_out);
    // Non-vacuous: filter-by-template drops the 4 odd templates whole (keeps 8
    // of 16), same as the RG/CB-free matrix test.
    assert_eq!(
        expected.len(),
        8,
        "oracle must keep exactly 8 of 16 records with RG/CB variation present"
    );
    assert_eq!(
        actual, expected,
        "chain output must match the non-chain path record-for-record with RG/CB variation \
         present, proving source_group_key_config's name_hash_only skip for Stage::Filter \
         changes no output"
    );
    assert_eq!(
        chain_header, oracle_header,
        "chain output header must match the non-chain path with RG/CB variation present"
    );
}

/// The `--rejects` output BAM matches record-for-record between the chain and
/// oracle paths, not just the kept output. Template mode over
/// `build_mixed_depth_templates` rejects both mates of every odd-indexed
/// template (8 records across 4 templates), so the rejects comparison is
/// non-vacuous.
#[test]
fn test_filter_chain_rejects_bam_matches_single_threaded() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let ref_path = create_test_reference(temp_dir.path());
    create_consensus_bam(&input_bam, build_mixed_depth_templates(8));

    let oracle_out = temp_dir.path().join("oracle.bam");
    let oracle_rejects = temp_dir.path().join("oracle.rejects.bam");
    filter_run(
        &input_bam,
        &oracle_out,
        &ref_path,
        &["--rejects", oracle_rejects.to_str().unwrap()],
    );

    let chain_out = temp_dir.path().join("chain.bam");
    let chain_rejects = temp_dir.path().join("chain.rejects.bam");
    filter_run(
        &input_bam,
        &chain_out,
        &ref_path,
        &["--rejects", chain_rejects.to_str().unwrap(), "--threads", "4"],
    );

    let expected_kept = read_filter_records(&oracle_out);
    let actual_kept = read_filter_records(&chain_out);
    assert!(!expected_kept.is_empty(), "oracle kept output must be non-empty");
    assert_eq!(actual_kept, expected_kept, "kept output must match record-for-record");

    let expected_rejects = read_filter_records(&oracle_rejects);
    let actual_rejects = read_filter_records(&chain_rejects);
    assert!(
        !expected_rejects.is_empty(),
        "oracle rejects output must be non-empty (guard against a vacuous pass)"
    );
    assert_eq!(
        actual_rejects, expected_rejects,
        "rejects output must match record-for-record between chain and oracle"
    );
}

/// The `--stats` file is byte-identical between the chain and oracle paths.
#[test]
fn test_filter_chain_stats_file_matches_single_threaded() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let ref_path = create_test_reference(temp_dir.path());
    create_consensus_bam(&input_bam, build_mixed_depth_templates(8));

    let oracle_out = temp_dir.path().join("oracle.bam");
    let oracle_stats = temp_dir.path().join("oracle.stats.txt");
    filter_run(&input_bam, &oracle_out, &ref_path, &["--stats", oracle_stats.to_str().unwrap()]);

    let chain_out = temp_dir.path().join("chain.bam");
    let chain_stats = temp_dir.path().join("chain.stats.txt");
    filter_run(
        &input_bam,
        &chain_out,
        &ref_path,
        &["--stats", chain_stats.to_str().unwrap(), "--threads", "4"],
    );

    let oracle_content = fs::read_to_string(&oracle_stats).expect("read oracle stats");
    let chain_content = fs::read_to_string(&chain_stats).expect("read chain stats");
    assert!(!oracle_content.trim().is_empty(), "oracle stats file should not be empty");
    assert!(oracle_content.contains("total_reads"), "stats should contain total_reads");
    // Non-vacuous: the stats must record that filtering actually rejected reads,
    // or a pass-through regression on both paths would still byte-match.
    let failed = oracle_content
        .lines()
        .find_map(|l| l.strip_prefix("failed_reads\t"))
        .and_then(|v| v.trim().parse::<u64>().ok())
        .expect("stats file must report failed_reads");
    assert!(failed > 0, "filtering must reject at least one read; stats:\n{oracle_content}");
    assert_eq!(
        oracle_content, chain_content,
        "stats file must be byte-identical across chain and oracle modes"
    );
}

/// Build a query-grouped BAM with `n` **paired** templates: each template `t{i}`
/// is two mates (R1/R2) sharing a query name, both mapped primaries, both
/// passing consensus reads (`cD 10 >= --min-reads 3`). Two properties make this
/// a genuine multi-batch / carry-over exercise (a one-record-per-name fixture
/// triggers neither):
/// - The `2 * n` records span multiple byte-aligned input record-batches
///   (`DecompressedBlock`s are record-aligned but NOT template-aligned), so some
///   template's two mates land in consecutive input batches — exercising
///   `GroupByQueryname`'s cross-input-batch `current_template` carry-over
///   (`pipeline/steps/group/queryname.rs`), which a single-record-per-name
///   fixture never triggers (every group flushes within one batch).
/// - The `n` templates also push past `GroupByQueryname`'s 1000-template
///   `DEFAULT_TARGET_BATCH_COUNT`, so its output-emit boundary is crossed too.
fn create_multi_batch_filter_input(path: &Path, n: usize) {
    let mut records = Vec::with_capacity(n * 2);
    for i in 0..i32::try_from(n).expect("n fits in i32") {
        let name = format!("t{i}");
        let pos = 100 + i;
        for (flag, mate_pos) in [
            (flags::PAIRED | flags::FIRST_SEGMENT, pos),
            (flags::PAIRED | flags::LAST_SEGMENT, pos + 8),
        ] {
            let mut b = RawSamBuilder::new();
            b.read_name(name.as_bytes())
                .flags(flag)
                .ref_id(0)
                .pos(mate_pos)
                .mapq(60)
                .cigar_ops(&[8 << 4]) // 8M
                .sequence(b"ACGTACGT")
                .qualities(&[35; 8]);
            b.add_int_tag(SamTag::CD, 10).add_float_tag(SamTag::CE, 0.0_f32);
            b.add_array_u16(SamTag::CD_BASES, &[10; 8]).add_array_u16(SamTag::CE_BASES, &[0; 8]);
            records.push(b.build());
        }
    }
    create_consensus_bam(path, records);
}

/// The chain's multi-worker path (`--threads 4`) must match the single-threaded
/// oracle across `GroupByQueryname`'s batch machinery: both its cross-input-batch
/// `current_template` carry-over (a paired template whose two mates land in
/// consecutive input record-batches) AND its 1000-template output-emit boundary.
/// `create_multi_batch_filter_input(1200)` emits 1200 paired templates (2400
/// records) to exercise both -- the small-fixture parity tests above (well under
/// one batch, single-record templates) exercise neither. Compare `--threads 4`
/// output record-for-record against the oracle.
#[test]
fn test_filter_chain_threads4_matches_single_threaded_multi_batch() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let ref_path = create_test_reference(temp_dir.path());
    create_multi_batch_filter_input(&input_bam, 1200);

    let oracle_out = temp_dir.path().join("oracle.bam");
    filter_run(&input_bam, &oracle_out, &ref_path, &[]);

    let chain_out = temp_dir.path().join("chain.bam");
    filter_run(&input_bam, &chain_out, &ref_path, &["--threads", "4"]);

    let expected = read_filter_records(&oracle_out);
    let actual = read_filter_records(&chain_out);
    // 1200 paired templates, all passing -> all 2400 records kept on both paths.
    assert_eq!(
        expected.len(),
        2400,
        "expected all 2400 records (1200 paired passing templates) in the oracle output"
    );
    assert_eq!(
        actual, expected,
        "chain (--threads 4) output must match the single-threaded oracle record-for-record \
         across GroupByQueryname input-batch carry-over and output-emit boundaries"
    );
}

/// FILT3-02 on the chain path: mirrors `test_filter_rejects_coordinate_sorted_input`
/// but with `--threads 4`, proving the `require_query_grouped` guard added to
/// `add_filter` (builder.rs) rejects coordinate-sorted input on the chain path
/// exactly as the legacy path does -- not just when the legacy oracle is used.
#[test]
fn test_filter_chain_rejects_coordinate_sorted_input() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let ref_path = create_test_reference(temp_dir.path());

    let header = create_coordinate_sorted_header("chr1", 10000);
    create_consensus_bam_with_header(&input_bam, &header, filter_consensus_records());

    let cmd = Filter::try_parse_from([
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
        "--threads",
        "4",
    ])
    .expect("failed to parse filter args");

    let err = cmd
        .execute("fgumi filter")
        .expect_err("chain path must also reject coordinate-sorted input");
    let msg = err.to_string();
    assert!(msg.contains("queryname sorted or query grouped"), "unexpected error message: {msg}");
}

/// FILT3-02 on the chain path, header-less variant: mirrors
/// `test_filter_rejects_headerless_input` but with `--threads 4`. Both paths now
/// synthesize `@HD VN:1.6 SO:unsorted` for header-less input — the legacy path
/// via `ensure_hd_record` in `execute()`, the chain path via `ensure_hd_record`
/// in `ChainBuilder::new` (before `add_pg_record`) — and `require_query_grouped`
/// rejects `SO:unsorted`, so the chain path rejects header-less input with the
/// same message. Pins chain-path header-less rejection so a future
/// header-handling refactor cannot silently break it.
#[test]
fn test_filter_chain_rejects_headerless_input() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let ref_path = create_test_reference(temp_dir.path());

    write_headerless_consensus_bam(&input_bam);

    let cmd = Filter::try_parse_from([
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
        "--threads",
        "4",
    ])
    .expect("failed to parse filter args");

    let err =
        cmd.execute("fgumi filter").expect_err("chain path must also reject header-less input");
    let msg = err.to_string();
    assert!(msg.contains("queryname sorted or query grouped"), "unexpected error message: {msg}");
}

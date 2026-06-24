//! End-to-end tests for the dedup command.
//!
//! These tests invoke `MarkDuplicates::execute()` in-process and validate:
//! 1. Basic duplicate marking
//! 2. Metrics output
//! 3. Remove duplicates mode

use clap::Parser;
use fgumi_lib::commands::command::Command;
use fgumi_lib::commands::dedup::MarkDuplicates;
use fgumi_lib::sam::SamTag;
use fgumi_raw_bam::{RawRecord, SamBuilder, flags};
use noodles::bam;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use std::fs;
use std::path::PathBuf;
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, to_record_buf};

/// Create a template-coordinate sorted BAM with UMI-tagged reads.
///
/// Template-coordinate sort groups reads by position, then by name within each position.
/// The header must have SO:unsorted GO:query SS:template-coordinate tags.
fn create_sorted_bam(path: &PathBuf, records: Vec<RawRecord>) {
    let header = create_minimal_header("chr1", 10000);
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");

    for record in records {
        writer
            .write_alignment_record(&header, &to_record_buf(&record))
            .expect("Failed to write record");
    }
    writer.try_finish().expect("Failed to finish BAM");
}

/// Create a group of paired-end reads at the same position with the same UMI
/// (simulating PCR duplicates).
fn create_duplicate_group(base_name: &str, umi: &str, count: usize, start: i32) -> Vec<RawRecord> {
    let mut records = Vec::new();
    for i in 0..count {
        let name = format!("{base_name}_{i}");

        let r1 = {
            let mut b = SamBuilder::new();
            b.read_name(name.as_bytes())
                .sequence(b"ACGTACGT")
                .qualities(&[30; 8])
                .flags(flags::PAIRED | flags::FIRST_SEGMENT)
                .ref_id(0)
                .pos(start - 1)
                .mapq(60)
                .cigar_ops(&[8 << 4]) // 8M
                .mate_ref_id(0)
                .mate_pos(start + 99)
                .template_length(108)
                .add_string_tag(SamTag::RX, umi.as_bytes())
                .add_string_tag(SamTag::MC, b"8M");
            b.build()
        };

        let r2 = {
            let mut b = SamBuilder::new();
            b.read_name(name.as_bytes())
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
                .add_string_tag(SamTag::RX, umi.as_bytes())
                .add_string_tag(SamTag::MC, b"8M");
            b.build()
        };

        records.push(r1);
        records.push(r2);
    }
    records
}

/// Test basic dedup command (mark duplicates).
#[test]
fn test_dedup_command_basic() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    // Create 3 duplicate pairs at position 100 with same UMI, and 2 at position 500
    let mut records = create_duplicate_group("dup1", "ACGTACGT", 3, 100);
    records.extend(create_duplicate_group("dup2", "TGCATGCA", 2, 500));
    create_sorted_bam(&input_bam, records);

    let cmd = MarkDuplicates::try_parse_from([
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
    .expect("failed to parse dedup args");
    cmd.execute("fgumi dedup").expect("Dedup command failed");
    assert!(output_bam.exists(), "Output BAM not created");

    // All reads should be present (duplicates are marked, not removed)
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let count = reader.records().count();
    assert_eq!(count, 10, "All reads should be in output (marked, not removed)");
}

/// Test dedup command with metrics output.
#[test]
fn test_dedup_command_with_metrics() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let metrics_path = temp_dir.path().join("metrics.txt");

    let records = create_duplicate_group("dup1", "ACGTACGT", 3, 100);
    create_sorted_bam(&input_bam, records);

    let cmd = MarkDuplicates::try_parse_from([
        "dedup",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--strategy",
        "identity",
        "--metrics",
        metrics_path.to_str().unwrap(),
        "--compression-level",
        "1",
    ])
    .expect("failed to parse dedup args");
    cmd.execute("fgumi dedup").expect("Dedup command with metrics failed");
    assert!(metrics_path.exists(), "Metrics file not created");
}

/// Test dedup command with remove-duplicates flag.
#[test]
fn test_dedup_command_remove_duplicates() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    // 3 duplicate pairs → should keep 1 pair (2 records), remove 2 pairs (4 records)
    let records = create_duplicate_group("dup1", "ACGTACGT", 3, 100);
    create_sorted_bam(&input_bam, records);

    let cmd = MarkDuplicates::try_parse_from([
        "dedup",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--strategy",
        "identity",
        "--remove-duplicates",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse dedup args");
    cmd.execute("fgumi dedup").expect("Dedup command with --remove-duplicates failed");
    assert!(output_bam.exists(), "Output BAM not created");

    // With --remove-duplicates, exactly one pair (2 records) survives: the 3
    // duplicate pairs (6 records) collapse to the single best representative
    // pair, so the other 2 pairs (4 records) are dropped. Assert the exact
    // expected count rather than loose bounds (S9b-008).
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let count = reader.records().count();
    assert_eq!(count, 2, "remove-duplicates over 3 duplicate pairs must keep exactly one pair");
}

/// SAM-input parity: dedup's typed-step path accepts both BAM and SAM via
/// `InputSource::open_with_opts`. Writes the same record set to both `.bam`
/// and `.sam`, runs `fgumi dedup` on each, and confirms the output BAMs
/// agree on record name + flags (duplicate flag 0x400 is the load-bearing
/// bit dedup mutates).
#[test]
fn test_dedup_command_with_sam_input_matches_bam_baseline() {
    let temp_dir = TempDir::new().unwrap();
    let bam_input = temp_dir.path().join("input.bam");
    let sam_input = temp_dir.path().join("input.sam");
    let out_bam_baseline = temp_dir.path().join("output_bam.bam");
    let out_sam = temp_dir.path().join("output_sam.bam");

    // Mirror `test_dedup_command_basic` so we exercise the same shape:
    // two distinct UMI families at two positions.
    let mut records = create_duplicate_group("dup1", "ACGTACGT", 3, 100);
    records.extend(create_duplicate_group("dup2", "TGCATGCA", 2, 500));
    create_sorted_bam(&bam_input, records);

    // Convert BAM → SAM by replaying records through the noodles SAM writer.
    {
        use noodles::sam::alignment::io::Write as _;
        let mut bam_reader =
            bam::io::reader::Builder.build_from_path(&bam_input).expect("open bam");
        let header = bam_reader.read_header().expect("read bam header");
        let mut sam_writer =
            noodles::sam::io::Writer::new(fs::File::create(&sam_input).expect("create sam"));
        sam_writer.write_header(&header).expect("write sam header");
        for record in bam_reader.records() {
            let record = record.expect("bam record");
            sam_writer.write_alignment_record(&header, &record).expect("write sam record");
        }
    }

    let run_dedup = |input: &PathBuf, output: &PathBuf| {
        let status = std::process::Command::new(env!("CARGO_BIN_EXE_fgumi"))
            .args([
                "dedup",
                "--input",
                input.to_str().unwrap(),
                "--output",
                output.to_str().unwrap(),
                "--strategy",
                "identity",
                "--compression-level",
                "1",
                "--threads",
                "2",
            ])
            .status()
            .expect("Failed to run dedup");
        assert!(status.success(), "dedup failed on {}", input.display());
    };

    run_dedup(&bam_input, &out_bam_baseline);
    run_dedup(&sam_input, &out_sam);

    // Compare record name + flags between the two outputs. The duplicate
    // flag (0x400) is set by dedup's mark step, so identical flags across
    // both runs proves the SAM source preamble feeds the same records
    // into the chain as the BAM preamble.
    let read_records = |path: &PathBuf| {
        let mut reader = bam::io::reader::Builder.build_from_path(path).expect("open bam");
        let _header = reader.read_header().expect("read header");
        reader.records().map(|r| r.expect("read record")).collect::<Vec<_>>()
    };
    let bam_records = read_records(&out_bam_baseline);
    let sam_records = read_records(&out_sam);

    assert_eq!(
        bam_records.len(),
        sam_records.len(),
        "BAM and SAM inputs produced different record counts ({} vs {})",
        bam_records.len(),
        sam_records.len(),
    );
    for (i, (bam_rec, sam_rec)) in bam_records.iter().zip(sam_records.iter()).enumerate() {
        assert_eq!(bam_rec.name(), sam_rec.name(), "Record {i} name mismatch");
        assert_eq!(
            bam_rec.flags(),
            sam_rec.flags(),
            "Record {i} flag mismatch (duplicate-flag 0x400 disagrees between BAM and SAM input)",
        );
    }
}

/// Regression test for OOM with large position groups in `--no-umi` mode.
///
/// WES data can have extreme depth pileups at capture targets, creating positions
/// with thousands of reads. With `--no-umi` (identity strategy), ALL reads at the
/// same position form ONE group. This test exercises a 5,000-template position
/// group to verify the pipeline completes without unbounded memory growth.
#[test]
#[allow(clippy::too_many_lines)]
fn test_dedup_no_umi_large_position_group() {
    use bstr::BString;
    use fgumi_raw_bam::raw_record_to_record_buf;
    use fgumi_raw_bam::{SamBuilder as RawSamBuilder, flags, testutil::encode_op};
    use noodles::sam::Header;
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::RecordBuf;
    use noodles::sam::alignment::record_buf::data::field::value::Value as DataValue;
    use noodles::sam::header::record::value::map::Map as HeaderRecordMap;
    use noodles::sam::header::record::value::map::header::tag::Tag as HeaderTag;
    use noodles::sam::header::record::value::{
        Map, map::Header as HeaderRecord, map::ReferenceSequence,
    };
    use std::collections::HashSet;
    use std::num::NonZeroUsize;

    const NUM_TEMPLATES: usize = 5_000;
    // PROPERLY_PAIRED flag = 0x2 (not exposed as a named constant in flags module).
    const PROPERLY_PAIRED: u16 = 0x2;
    // MATE_REVERSE flag = 0x20 (not exposed as a named constant in flags module).
    // R1 is forward / R2 is reverse; we set MATE_REVERSE on R1 so the pair is FR
    // self-consistent (mate reverse flag matches mate's REVERSE).
    const MATE_REVERSE: u16 = 0x20;

    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    // Build a template-coordinate sorted header with one reference.
    let header = {
        let mut builder = HeaderRecordMap::<HeaderRecord>::builder();
        for (tag_bytes, value) in
            [(*b"SO", "unsorted"), (*b"GO", "query"), (*b"SS", "template-coordinate")]
        {
            let HeaderTag::Other(tag) = HeaderTag::from(tag_bytes) else { unreachable!() };
            builder = builder.insert(tag, value);
        }
        Header::builder()
            .set_header(builder.build().expect("valid header map"))
            .add_reference_sequence(
                BString::from("chr1"),
                Map::<ReferenceSequence>::new(NonZeroUsize::new(10_000).expect("non-zero length")),
            )
            .build()
    };

    // Write 5,000 paired-end templates all at position 100 with no UMI tag.
    // Template-coordinate sort groups by position then name, so we write
    // R1 and R2 together for each template in name-sorted order.
    {
        let mut writer = bam::io::Writer::new(fs::File::create(&input_bam).expect("create BAM"));
        writer.write_header(&header).expect("write header");

        for i in 0..NUM_TEMPLATES {
            let name = format!("read_{i:05}");

            let r1 = {
                let mut b = RawSamBuilder::new();
                b.read_name(name.as_bytes())
                    .sequence(b"ACGTACGT")
                    .qualities(&[30; 8])
                    .flags(flags::PAIRED | PROPERLY_PAIRED | flags::FIRST_SEGMENT | MATE_REVERSE)
                    .ref_id(0)
                    .pos(99)
                    .mapq(60)
                    .cigar_ops(&[encode_op(0, 8)])
                    .mate_ref_id(0)
                    .mate_pos(199)
                    .template_length(108);
                b.add_string_tag(SamTag::MC, b"8M");
                b.build()
            };

            let r2 = {
                let mut b = RawSamBuilder::new();
                b.read_name(name.as_bytes())
                    .sequence(b"ACGTACGT")
                    .qualities(&[30; 8])
                    .flags(flags::PAIRED | PROPERLY_PAIRED | flags::REVERSE | flags::LAST_SEGMENT)
                    .ref_id(0)
                    .pos(199)
                    .mapq(60)
                    .cigar_ops(&[encode_op(0, 8)])
                    .mate_ref_id(0)
                    .mate_pos(99)
                    .template_length(-108);
                b.add_string_tag(SamTag::MC, b"8M");
                b.build()
            };

            let r1_buf = raw_record_to_record_buf(&r1, &header).expect("convert R1");
            let r2_buf = raw_record_to_record_buf(&r2, &header).expect("convert R2");
            writer.write_alignment_record(&header, &r1_buf).expect("write R1");
            writer.write_alignment_record(&header, &r2_buf).expect("write R2");
        }
        writer.try_finish().expect("finish BAM");
    }

    // Run dedup with --no-umi via MarkDuplicates::execute().
    let cmd = MarkDuplicates::try_parse_from([
        "dedup",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--no-umi",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse dedup args");
    cmd.execute("fgumi dedup").expect("dedup --no-umi should succeed with large position group");

    assert!(output_bam.exists(), "output BAM should be created");

    // Read back the output and verify duplicate marking and MI tags.
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).expect("open output BAM"));
    let out_header = reader.read_header().expect("read output header");

    let mut total_records = 0usize;
    let mut duplicate_records = 0usize;
    let mut non_duplicate_records = 0usize;
    let mut mi_values = HashSet::new();
    let mut mi_count = 0usize;
    let mut non_dup_names = HashSet::new();

    let mi_tag = Tag::from(SamTag::MI);

    for result in reader.record_bufs(&out_header) {
        let record: RecordBuf = result.expect("read record");
        total_records += 1;

        let is_dup = record.flags().is_duplicate();
        if is_dup {
            duplicate_records += 1;
        } else {
            non_duplicate_records += 1;
            if let Some(name) = record.name() {
                non_dup_names.insert(name.to_owned());
            }
        }

        if let Some(DataValue::String(mi)) = record.data().get(&mi_tag) {
            mi_values.insert(mi.to_owned());
            mi_count += 1;
        }
    }

    // All 10,000 records (5,000 pairs) should be present.
    assert_eq!(total_records, NUM_TEMPLATES * 2, "all records should be present in output");

    // Exactly one template (2 records) should NOT be marked as duplicate.
    assert_eq!(
        non_duplicate_records, 2,
        "exactly one pair should be non-duplicate (the best-scoring template)"
    );
    assert_eq!(
        duplicate_records,
        (NUM_TEMPLATES - 1) * 2,
        "all other pairs should be marked as duplicates"
    );

    // The two non-duplicate records should share the same read name.
    assert_eq!(non_dup_names.len(), 1, "non-duplicate records should be from one template");

    // Every record should have an MI tag.
    assert_eq!(mi_count, NUM_TEMPLATES * 2, "all records should have an MI tag");

    // All templates should share the same MI value (one group in identity strategy).
    assert_eq!(mi_values.len(), 1, "all records should share a single MI tag value");
}

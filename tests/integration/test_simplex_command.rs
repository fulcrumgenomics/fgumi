//! End-to-end CLI tests for the simplex command.
//!
//! These tests invoke `Simplex::execute()` in-process and validate:
//! 1. Basic simplex consensus calling from grouped reads
//! 2. Statistics output
//! 3. Rejected reads output

use bstr::BString;
use clap::Parser;
use fgumi_lib::commands::command::Command;
use fgumi_lib::commands::simplex::Simplex;
use fgumi_lib::sam::SamTag;
use noodles::bam;
use noodles::sam::Header;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::ReadGroup;
use noodles::sam::header::record::value::map::read_group::tag as rg_tag;
use std::fs;
use std::num::NonZeroUsize;
use std::path::Path;
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, create_umi_family, to_record_buf};

/// Write grouped BAM file (reads grouped by MI tag).
fn create_grouped_bam(path: &Path, families: Vec<(&str, Vec<fgumi_raw_bam::RawRecord>)>) {
    let header = create_minimal_header("chr1", 10000);
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");

    for (mi, records) in families {
        for raw in &records {
            // Convert to RecordBuf, add MI tag, write
            use noodles::sam::alignment::record::data::field::Tag;
            use noodles::sam::alignment::record_buf::data::field::Value;
            let mut record = to_record_buf(raw);
            let mi_tag = Tag::from(fgumi_lib::sam::SamTag::MI);
            record.data_mut().insert(mi_tag, Value::from(mi));
            writer.write_alignment_record(&header, &record).expect("Failed to write record");
        }
    }
    writer.try_finish().expect("Failed to finish BAM");
}

/// Test basic simplex consensus calling.
#[test]
fn test_simplex_command_basic_consensus() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    // Create two families: 5 reads each
    let family1 = create_umi_family("ACGT", 5, "fam1", "ACGTACGT", 30);
    let family2 = create_umi_family("TGCA", 5, "fam2", "TTTTAAAA", 30);
    create_grouped_bam(&input_bam, vec![("1", family1), ("2", family2)]);

    let cmd = Simplex::try_parse_from([
        "simplex",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--min-reads",
        "2",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse simplex args");
    cmd.execute("test").expect("Simplex command failed");
    assert!(output_bam.exists(), "Output BAM not created");

    // Read output and verify consensus reads were produced
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let mut count = 0;
    for result in reader.records() {
        let record = result.expect("Failed to read record");
        // Verify consensus tags exist
        let cd_tag = SamTag::CD.to_noodles_tag();
        assert!(record.data().get(&cd_tag).is_some(), "Consensus should have cD tag");
        count += 1;
    }
    assert!(count > 0, "Should have produced consensus reads");
}

/// Test simplex command with statistics output.
#[test]
fn test_simplex_command_with_stats() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let stats_path = temp_dir.path().join("stats.txt");

    let family = create_umi_family("ACGT", 5, "fam1", "ACGTACGT", 30);
    create_grouped_bam(&input_bam, vec![("1", family)]);

    let cmd = Simplex::try_parse_from([
        "simplex",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--min-reads",
        "1",
        "--stats",
        stats_path.to_str().unwrap(),
        "--compression-level",
        "1",
    ])
    .expect("failed to parse simplex args");
    cmd.execute("test").expect("Simplex command with stats failed");
    assert!(stats_path.exists(), "Stats file not created");
}

/// Test simplex command with rejects output.
#[test]
fn test_simplex_command_with_rejects() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let rejects_bam = temp_dir.path().join("rejects.bam");

    // Create one family with 5 reads (passes min-reads=2) and one with 1 read (fails)
    let family1 = create_umi_family("ACGT", 5, "fam1", "ACGTACGT", 30);
    let family2 = create_umi_family("TGCA", 1, "fam2", "TTTTAAAA", 30);
    create_grouped_bam(&input_bam, vec![("1", family1), ("2", family2)]);

    let cmd = Simplex::try_parse_from([
        "simplex",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--min-reads",
        "2",
        "--rejects",
        rejects_bam.to_str().unwrap(),
        "--compression-level",
        "1",
    ])
    .expect("failed to parse simplex args");
    cmd.execute("test").expect("Simplex command with rejects failed");
    assert!(rejects_bam.exists(), "Rejects BAM not created");
}

/// Creates a header with multiple read groups, each with SM, LB, PL tags.
fn create_header_with_read_groups(ref_name: &str, ref_len: usize) -> Header {
    use noodles::sam::header::record::value::map::Header as HeaderRecord;
    use noodles::sam::header::record::value::map::Map as HeaderRecordMap;
    use noodles::sam::header::record::value::map::ReferenceSequence;
    use noodles::sam::header::record::value::map::header::tag::Tag as HeaderTag;

    let HeaderTag::Other(sort_order_tag) = HeaderTag::from([b'S', b'O']) else { unreachable!() };
    let HeaderTag::Other(group_order_tag) = HeaderTag::from([b'G', b'O']) else { unreachable!() };
    let HeaderTag::Other(sub_sort_tag) = HeaderTag::from([b'S', b'S']) else { unreachable!() };

    let header_map = HeaderRecordMap::<HeaderRecord>::builder()
        .insert(sort_order_tag, "unsorted")
        .insert(group_order_tag, "query")
        .insert(sub_sort_tag, "template-coordinate")
        .build()
        .expect("valid header map");

    let reference_sequence = Map::<ReferenceSequence>::new(
        NonZeroUsize::new(ref_len).expect("reference length must be non-zero"),
    );

    let rg1 = Map::<ReadGroup>::builder()
        .insert(rg_tag::SAMPLE, String::from("SampleA"))
        .insert(rg_tag::LIBRARY, String::from("LibA"))
        .insert(rg_tag::PLATFORM, String::from("illumina"))
        .build()
        .expect("valid RG1");

    let rg2 = Map::<ReadGroup>::builder()
        .insert(rg_tag::SAMPLE, String::from("SampleA"))
        .insert(rg_tag::LIBRARY, String::from("LibB"))
        .insert(rg_tag::PLATFORM, String::from("ILLUMINA"))
        .build()
        .expect("valid RG2");

    Header::builder()
        .set_header(header_map)
        .add_reference_sequence(BString::from(ref_name), reference_sequence)
        .add_read_group(BString::from("RG1"), rg1)
        .add_read_group(BString::from("RG2"), rg2)
        .build()
}

/// Write grouped BAM file with a custom header that has multiple read groups.
fn create_grouped_bam_with_header(
    path: &Path,
    header: &Header,
    families: Vec<(&str, Vec<fgumi_raw_bam::RawRecord>)>,
) {
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(header).expect("Failed to write header");

    for (mi, records) in families {
        for raw in &records {
            use noodles::sam::alignment::record::data::field::Tag;
            use noodles::sam::alignment::record_buf::data::field::Value;
            let mut record = to_record_buf(raw);
            let mi_tag = Tag::from(fgumi_lib::sam::SamTag::MI);
            record.data_mut().insert(mi_tag, Value::from(mi));
            writer.write_alignment_record(header, &record).expect("Failed to write record");
        }
    }
    writer.finish(header).expect("Failed to finish BAM");
}

/// Test that simplex output header collapses read group attributes from input.
#[test]
fn test_simplex_command_collapses_read_group_attributes() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    // Create input BAM with two read groups having different LB but same SM
    let header = create_header_with_read_groups("chr1", 10000);
    let family = create_umi_family("ACGT", 5, "fam1", "ACGTACGT", 30);
    create_grouped_bam_with_header(&input_bam, &header, vec![("1", family)]);

    let cmd = Simplex::try_parse_from([
        "simplex",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--min-reads",
        "2",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse simplex args");
    cmd.execute("test").expect("Simplex command failed");

    // Read the output header and verify collapsed read group attributes
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let output_header = reader.read_header().unwrap();
    let read_groups = output_header.read_groups();

    assert_eq!(read_groups.len(), 1, "Should have exactly one output read group");

    let rg = read_groups.get(&BString::from("A")).expect("Output RG 'A' not found");

    // SM: both input RGs have "SampleA" -> deduplicated to single value
    assert_eq!(
        rg.other_fields().get(&rg_tag::SAMPLE).map(std::string::ToString::to_string),
        Some("SampleA".to_string()),
        "SM tag should be 'SampleA'",
    );

    // LB: "LibA" and "LibB" -> comma-joined
    assert_eq!(
        rg.other_fields().get(&rg_tag::LIBRARY).map(std::string::ToString::to_string),
        Some("LibA,LibB".to_string()),
        "LB tag should be comma-joined distinct values",
    );

    // PL: "illumina" and "ILLUMINA" -> uppercased and deduplicated
    assert_eq!(
        rg.other_fields().get(&rg_tag::PLATFORM).map(std::string::ToString::to_string),
        Some("ILLUMINA".to_string()),
        "PL tag should be uppercased and deduplicated",
    );
}

/// Verifies that the rejects BAM advertises the **input** header (RGs included),
/// while the primary output BAM advertises the **consensus** header (RGs collapsed).
///
/// This is the end-to-end check for the `secondary_output_header` parameter on
/// [`run_bam_pipeline_from_reader_with_secondary`][rbprws]. A regression where a
/// caller accidentally passes the primary `output_header` (consensus header) for
/// the secondary would collapse rejects' RGs to "A" — this test catches that.
///
/// [rbprws]: fgumi_lib::unified_pipeline::run_bam_pipeline_from_reader_with_secondary
#[test]
fn test_simplex_command_rejects_inherits_input_read_groups() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let rejects_bam = temp_dir.path().join("rejects.bam");

    // Build an input with two distinct read groups (RG1, RG2) that the
    // consensus header would collapse to a single RG "A".
    let header = create_header_with_read_groups("chr1", 10000);
    // One family that passes consensus (kept), one singleton that doesn't (rejected).
    let kept = create_umi_family("ACGT", 5, "kept", "ACGTACGT", 30);
    let rejected = create_umi_family("TGCA", 1, "reject_singleton", "TTTTAAAA", 30);
    create_grouped_bam_with_header(&input_bam, &header, vec![("1", kept), ("2", rejected)]);

    let cmd = Simplex::try_parse_from([
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
    .expect("failed to parse simplex args");
    cmd.execute("test").expect("Simplex command with rejects failed");

    // Primary output BAM: consensus header collapses RG1+RG2 -> single "A".
    let primary_header =
        bam::io::Reader::new(fs::File::open(&output_bam).unwrap()).read_header().unwrap();
    let primary_rgs: Vec<BString> = primary_header.read_groups().keys().cloned().collect();
    assert_eq!(
        primary_rgs,
        vec![BString::from("A")],
        "primary output should have the collapsed consensus RG 'A'",
    );

    // Rejects BAM: input header preserved verbatim — both RG1 and RG2 present.
    let mut reader = bam::io::Reader::new(fs::File::open(&rejects_bam).unwrap());
    let rejects_header = reader.read_header().unwrap();
    let mut rejects_rgs: Vec<BString> = rejects_header.read_groups().keys().cloned().collect();
    rejects_rgs.sort();
    assert_eq!(
        rejects_rgs,
        vec![BString::from("RG1"), BString::from("RG2")],
        "rejects BAM should inherit the input header's read groups verbatim",
    );

    // Sanity-check that the rejection path actually fired: the singleton
    // input record (TGCA family, depth=1) must land in the rejects BAM.
    // Assert the read *name* (not just the count) so a regression that
    // rejects the wrong family (e.g. the kept "ACGT" family instead of the
    // singleton "TGCA" family) is caught — a count-only check would still
    // pass on a swap.
    let reject_names: Vec<String> = reader
        .records()
        .map(|result| {
            result
                .expect("Failed to read rejects record")
                .name()
                .expect("reject record missing read name")
                .to_string()
        })
        .collect();
    assert_eq!(
        reject_names.len(),
        1,
        "rejects BAM should contain the singleton record dropped by --min-reads",
    );
    assert!(
        reject_names[0].starts_with("reject_singleton"),
        "expected the singleton family in rejects, got {reject_names:?}",
    );
}

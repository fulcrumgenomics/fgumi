//! End-to-end CLI tests for the simplex command.
//!
//! These tests invoke `Simplex::execute()` in-process and validate:
//! 1. Basic simplex consensus calling from grouped reads
//! 2. Statistics output
//! 3. Rejected reads output

use bstr::BString;
use clap::Parser;
use fgumi_dna::reverse_complement;
use fgumi_lib::commands::command::Command;
use fgumi_lib::commands::simplex::Simplex;
use fgumi_lib::sam::SamTag;
use fgumi_raw_bam::flags as raw_flags;
use noodles::bam;
use noodles::sam::Header;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::ReadGroup;
use noodles::sam::header::record::value::map::read_group::tag as rg_tag;
use rstest::rstest;
use std::fs;
use std::num::NonZeroUsize;
use std::path::Path;
use tempfile::TempDir;

use crate::helpers::assertions::{int_tag, string_tag};
use crate::helpers::bam_generator::{
    create_coordinate_sorted_header, create_minimal_header, create_paired_umi_family,
    create_umi_family, to_record_buf,
};
use crate::helpers::read_bam_output;

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
    cmd.execute("fgumi simplex").expect("Simplex command failed");
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
    cmd.execute("fgumi simplex").expect("Simplex command with stats failed");
    assert!(stats_path.exists(), "Stats file not created");
}

/// A simplex consensus record reduced to the fields a caller cannot get right by accident:
/// its name, flags, sequence, molecule identifier, UMI, and the raw-read depth it claims.
#[derive(Debug, PartialEq, Eq)]
struct SimplexConsensusIdentity {
    name: String,
    flags: u16,
    sequence: String,
    mi: String,
    rx: String,
    depth: i64,
}

impl SimplexConsensusIdentity {
    /// Project an emitted BAM record onto the fields this test pins.
    fn from_record(record: &bam::Record) -> Self {
        Self {
            name: String::from_utf8_lossy(
                record.name().expect("consensus record must be named").as_ref(),
            )
            .into_owned(),
            flags: record.flags().bits(),
            sequence: record.sequence().iter().map(char::from).collect(),
            mi: string_tag(record, SamTag::MI),
            rx: string_tag(record, SamTag::RX),
            depth: int_tag(record, SamTag::CD),
        }
    }
}

/// The reported `consensus_reads_emitted` must equal the number of consensus records
/// actually written, in both the single-threaded and pipeline (`--threads N`) paths.
///
/// `ConsensusCallingStats::consensus_reads` counts emitted reads, so a caller that counts
/// molecules or templates instead reports a number that does not match its own output.
#[rstest]
#[case::single_threaded(None)]
#[case::threaded(Some("2"))]
fn test_simplex_stats_consensus_reads_matches_records_written(#[case] threads: Option<&str>) {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let stats_path = temp_dir.path().join("stats.txt");

    // Paired families, so the caller's R1/R2 branch -- the one that must record two
    // consensus reads per template -- is the branch under test.
    let families = (1..=3)
        .map(|i| {
            create_paired_umi_family("ACGT", 3, &format!("fam{i}"), "ACGTACGT", "TTTTAAAA", 30)
        })
        .collect::<Vec<_>>();
    create_grouped_bam(
        &input_bam,
        vec![("1", families[0].clone()), ("2", families[1].clone()), ("3", families[2].clone())],
    );

    let mut args = vec![
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
    ];
    if let Some(threads) = threads {
        args.extend_from_slice(&["--threads", threads]);
    }
    Simplex::try_parse_from(args)
        .expect("failed to parse simplex args")
        .execute("fgumi simplex")
        .expect("Simplex command failed");

    // Identity, not just cardinality: every family is built from three identical read
    // pairs, so each consensus pair is fully determined -- an R1/R2 pair per family,
    // named and MI-tagged by the family's molecule id, carrying that family's two
    // sequences, its UMI, and its three raw reads per strand.
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let records =
        reader.records().collect::<Result<Vec<_>, _>>().expect("failed to read the consensus BAM");
    let observed = records.iter().map(SimplexConsensusIdentity::from_record).collect::<Vec<_>>();
    let unmapped_pair = raw_flags::PAIRED | raw_flags::UNMAPPED | raw_flags::MATE_UNMAPPED;
    let expected = (1..=3)
        .flat_map(|family| {
            [(raw_flags::FIRST_SEGMENT, "ACGTACGT"), (raw_flags::LAST_SEGMENT, "TTTTAAAA")].map(
                |(segment, sequence)| SimplexConsensusIdentity {
                    name: format!(":{family}"),
                    flags: unmapped_pair | segment,
                    sequence: sequence.to_string(),
                    mi: family.to_string(),
                    rx: "ACGT".to_string(),
                    depth: 3,
                },
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(observed, expected, "simplex must emit one fully determined R1/R2 pair per family");

    let records_written = records.len();

    let stats = fs::read_to_string(&stats_path).expect("read stats");
    let emitted: usize = stats
        .lines()
        .find_map(|line| line.strip_prefix("consensus_reads_emitted\t"))
        .and_then(|rest| rest.split('\t').next())
        .expect("stats must contain a consensus_reads_emitted row")
        .parse()
        .expect("consensus_reads_emitted must be an integer");

    assert_eq!(
        emitted, records_written,
        "consensus_reads_emitted must count consensus reads written to the output BAM"
    );
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
    cmd.execute("fgumi simplex").expect("Simplex command with rejects failed");
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
    cmd.execute("fgumi simplex").expect("Simplex command failed");

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
    cmd.execute("fgumi simplex").expect("Simplex command with rejects failed");

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

/// Builds a tag family of `mapped_pairs` fully mapped pairs plus one pair whose R2 is unmapped.
///
/// The mapped R2s carry `r2_mapped_seq` and the unmapped R2 carries `r2_unmapped_seq`, so the
/// emitted R2 consensus shows which reads contributed to it. All records share one position so
/// they land in a single tag family.
fn half_mapped_family(
    mapped_pairs: usize,
    r1_seq: &str,
    r2_mapped_seq: &str,
    r2_unmapped_seq: &str,
) -> Vec<fgumi_raw_bam::RawRecord> {
    use fgumi_raw_bam::{SamBuilder, flags};

    let r1_cigar = u32::try_from(r1_seq.len()).expect("fits u32") << 4;
    let r2_cigar = u32::try_from(r2_mapped_seq.len()).expect("fits u32") << 4;
    let template_len = i32::try_from(100 + r2_mapped_seq.len()).expect("fits i32");
    let mut records = Vec::new();

    let mut push_pair = |name: &str, r2_seq: &str, r2_unmapped: bool| {
        let mut b1 = SamBuilder::new();
        b1.read_name(name.as_bytes())
            .ref_id(0)
            .pos(99)
            .mapq(60)
            .flags(flags::PAIRED | flags::FIRST_SEGMENT)
            .mate_ref_id(0)
            .mate_pos(199)
            .template_length(template_len)
            .cigar_ops(&[r1_cigar])
            .sequence(r1_seq.as_bytes())
            .qualities(&vec![30u8; r1_seq.len()])
            .add_string_tag(SamTag::RX, b"ACGT");
        records.push(b1.build());

        let mut b2 = SamBuilder::new();
        let r2_flags = if r2_unmapped {
            flags::PAIRED | flags::LAST_SEGMENT | flags::UNMAPPED
        } else {
            flags::PAIRED | flags::LAST_SEGMENT
        };
        b2.read_name(name.as_bytes())
            .ref_id(0)
            .pos(199)
            .mapq(if r2_unmapped { 0 } else { 60 })
            .flags(r2_flags)
            .mate_ref_id(0)
            .mate_pos(99)
            .template_length(-template_len)
            .sequence(r2_seq.as_bytes())
            .qualities(&vec![30u8; r2_seq.len()])
            .add_string_tag(SamTag::RX, b"ACGT");
        if !r2_unmapped {
            b2.cigar_ops(&[r2_cigar]);
        }
        records.push(b2.build());
    };

    for i in 0..mapped_pairs {
        push_pair(&format!("mapped_{i}"), r2_mapped_seq, false);
    }
    push_pair("halfmapped", r2_unmapped_seq, true);

    records
}

/// The unmapped end of a half-mapped pair must not contribute to a consensus when mapped reads
/// are present in the same family. fgumi always prefers mapped reads.
#[test]
fn test_simplex_ignores_unmapped_end_when_mapped_reads_present() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    let family = half_mapped_family(2, "ACGTACGTAC", "GGGGGGGGGG", "TTTTTTTTTT");
    create_grouped_bam(&input_bam, vec![("1", family)]);

    let cmd = Simplex::try_parse_from([
        "simplex",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--min-reads",
        "1",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse simplex args");
    cmd.execute("fgumi simplex").expect("Simplex command failed");

    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let records: Vec<bam::Record> =
        reader.records().map(|r| r.expect("failed to read record")).collect();
    assert_eq!(records.len(), 2, "expected one consensus read pair");

    let r1 = records
        .iter()
        .find(|r| r.flags().bits() & raw_flags::FIRST_SEGMENT != 0)
        .expect("missing R1 consensus");
    let r2 = records
        .iter()
        .find(|r| r.flags().bits() & raw_flags::LAST_SEGMENT != 0)
        .expect("missing R2 consensus");

    // All three R1s are mapped, so all three contribute.
    assert_eq!(int_tag(r1, SamTag::CD), 3, "R1 should use all three mapped reads");

    // Only the two mapped R2s may contribute; the unmapped R2's bases must not reach the consensus.
    assert_eq!(int_tag(r2, SamTag::CD), 2, "R2 must use only the two mapped reads");
    let r2_bases: String = r2.sequence().iter().map(char::from).collect();
    assert_eq!(r2_bases, "GGGGGGGGGG", "R2 consensus must come from the mapped reads");
}

/// Under `--allow-unmapped`, unmapped reads are admitted to grouping, but must still not displace
/// mapped reads within a set being consensus called.
#[test]
fn test_simplex_prefers_mapped_reads_even_when_unmapped_are_allowed() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    let family = half_mapped_family(2, "ACGTACGTAC", "GGGGGGGGGG", "TTTTTTTTTT");
    create_grouped_bam(&input_bam, vec![("1", family)]);

    let cmd = Simplex::try_parse_from([
        "simplex",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--min-reads",
        "1",
        "--allow-unmapped",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse simplex args");
    cmd.execute("fgumi simplex").expect("Simplex command failed");

    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let records: Vec<bam::Record> =
        reader.records().map(|r| r.expect("failed to read record")).collect();
    let r2 = records
        .iter()
        .find(|r| r.flags().bits() & raw_flags::LAST_SEGMENT != 0)
        .expect("missing R2 consensus");

    assert_eq!(int_tag(r2, SamTag::CD), 2, "mapped reads must win even with --allow-unmapped");
    let r2_bases: String = r2.sequence().iter().map(char::from).collect();
    assert_eq!(r2_bases, "GGGGGGGGGG", "R2 consensus must come from the mapped reads");
}

/// `--allow-unmapped` must still produce a consensus for a wholly unmapped family; this is the
/// capability (e.g. ribosome display) that keeps fgumi from simply requiring mapped reads.
#[test]
fn test_simplex_allow_unmapped_still_consenses_a_wholly_unmapped_family() {
    use fgumi_raw_bam::{SamBuilder, flags};

    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    let mut records = Vec::new();
    for i in 0..3 {
        for (segment, seq) in
            [(flags::FIRST_SEGMENT, "ACGTACGTAC"), (flags::LAST_SEGMENT, "GGGGGGGGGG")]
        {
            let mut b = SamBuilder::new();
            b.read_name(format!("unmapped_{i}").as_bytes())
                .ref_id(-1)
                .pos(-1)
                .mapq(0)
                .flags(flags::PAIRED | segment | flags::UNMAPPED | flags::MATE_UNMAPPED)
                .mate_ref_id(-1)
                .mate_pos(-1)
                .sequence(seq.as_bytes())
                .qualities(&vec![30u8; seq.len()])
                .add_string_tag(SamTag::RX, b"ACGT");
            records.push(b.build());
        }
    }
    create_grouped_bam(&input_bam, vec![("1", records)]);

    let cmd = Simplex::try_parse_from([
        "simplex",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--min-reads",
        "1",
        "--allow-unmapped",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse simplex args");
    cmd.execute("fgumi simplex").expect("Simplex command failed");

    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let records: Vec<bam::Record> =
        reader.records().map(|r| r.expect("failed to read record")).collect();
    assert_eq!(records.len(), 2, "a wholly unmapped family must still yield a consensus pair");
    for record in &records {
        assert_eq!(int_tag(record, SamTag::CD), 3, "all three unmapped reads should contribute");
    }
}

/// Read length of an [`indel_readthrough_family`] read: `2S124M1D3M` and `3S124M2S` both consume
/// 129 query bases.
const INDEL_READTHROUGH_READ_LEN: usize = 129;

/// R1's 1-based alignment start in [`indel_readthrough_family`]; R2 starts one base to its left.
const INDEL_READTHROUGH_R1_START: usize = 200;

/// A deterministic synthetic reference base at a 1-based position.
///
/// Deliberately aperiodic over the span the fixture uses, so a clip of the right size taken from
/// the wrong end, or one base off, changes the consensus rather than landing on an identical
/// repeat.
fn synthetic_ref_base(pos_1based: usize) -> u8 {
    b"ACGT"[(pos_1based * 7 + pos_1based / 5 + pos_1based / 11) % 4]
}

/// The reference-orientation bases of an [`indel_readthrough_family`] read pair, `(r1, r2)`.
///
/// Both reads take their aligned bases from [`synthetic_ref_base`], so they agree everywhere
/// their alignments overlap. R1's `1D` skips reference position 324, which R2 never reaches.
///
/// The soft-clipped adapter tails are a fixed called base rather than `N`: consensus calling
/// drops trailing no-call positions, which would shorten a correctly clipped read and hide the
/// very thing this fixture measures.
fn indel_readthrough_read_bases() -> (Vec<u8>, Vec<u8>) {
    const ADAPTER: u8 = b'A';
    let r1_start = INDEL_READTHROUGH_R1_START;

    // R1 = 2S124M1D3M: 2 adapter bases, 124 aligned from 200, a skipped 324, 3 aligned from 325.
    let mut r1 = vec![ADAPTER; 2];
    r1.extend((r1_start..r1_start + 124).map(synthetic_ref_base));
    r1.extend((r1_start + 125..r1_start + 128).map(synthetic_ref_base));

    // R2 = 3S124M2S at 199: 3 adapter bases, 124 aligned from 199, 2 adapter bases.
    let mut r2 = vec![ADAPTER; 3];
    r2.extend((r1_start - 1..r1_start + 123).map(synthetic_ref_base));
    r2.extend([ADAPTER; 2]);

    assert_eq!(r1.len(), INDEL_READTHROUGH_READ_LEN);
    assert_eq!(r2.len(), INDEL_READTHROUGH_READ_LEN);
    (r1, r2)
}

/// Builds `depth` read pairs that sequence through the insert and into the adapter, with a
/// 1-base deletion three bases from the positive read's 3' end — the fulcrumgenomics/fgumi#752
/// shape, in the end-to-end setting where its cost is visible.
///
/// - R1 (positive): `2S124M1D3M` at 1-based 200 (reference span 200-327)
/// - R2 (negative): `3S124M2S`   at 1-based 199 (reference span 199-322)
///
/// R2's soft-only unclipped *reference* end is its alignment end plus its 2-base trailing soft
/// clip — a query distance added to a reference coordinate — which lands inside R1's `1D`.
/// Measuring the clip in query space from the last reference position the two reads share (322)
/// gives R1 4 query bases past it against R2's 2, so R1 gives up 2 and keeps 127.
fn indel_readthrough_family(depth: usize) -> Vec<fgumi_raw_bam::RawRecord> {
    use fgumi_raw_bam::SamBuilder;

    // BAM stores SEQ on the forward reference strand for both mates, so the negative strand's
    // record carries its reference-orientation bases as-is. Getting this backwards matters here:
    // `simplex` consensus-calls overlapping bases by default, so two strands that disagree at a
    // reference position are both masked to `N` and the sequence assertions go blind.
    let (r1_seq, r2_seq) = indel_readthrough_read_bases();
    // 2S124M1D3M and 3S124M2S, packed as (len << 4) | op_code.
    let r1_cigar = [(2u32 << 4) | 4, (124u32 << 4), (1u32 << 4) | 2, 3u32 << 4];
    let r2_cigar = [(3u32 << 4) | 4, (124u32 << 4), (2u32 << 4) | 4];
    let quals = [35u8; INDEL_READTHROUGH_READ_LEN];
    // 0-based positions: R1 at 1-based 200, R2 one base to its left.
    let r1_pos = i32::try_from(INDEL_READTHROUGH_R1_START).expect("start fits i32") - 1;
    let r2_pos = r1_pos - 1;

    let mut records = Vec::new();
    for i in 0..depth {
        let read_name = format!("rt_{i}");

        let mut b1 = SamBuilder::new();
        b1.read_name(read_name.as_bytes())
            .ref_id(0)
            .pos(r1_pos)
            .mapq(60)
            .flags(raw_flags::PAIRED | raw_flags::FIRST_SEGMENT | raw_flags::MATE_REVERSE)
            .mate_ref_id(0)
            .mate_pos(r2_pos)
            .template_length(124)
            .cigar_ops(&r1_cigar)
            .sequence(&r1_seq)
            .qualities(&quals)
            .add_string_tag(SamTag::RX, b"ACGT")
            .add_string_tag(SamTag::MC, b"3S124M2S");
        records.push(b1.build());

        let mut b2 = SamBuilder::new();
        b2.read_name(read_name.as_bytes())
            .ref_id(0)
            .pos(r2_pos)
            .mapq(60)
            .flags(raw_flags::PAIRED | raw_flags::LAST_SEGMENT | raw_flags::REVERSE)
            .mate_ref_id(0)
            .mate_pos(r1_pos)
            .template_length(-124)
            .cigar_ops(&r2_cigar)
            .sequence(&r2_seq)
            .qualities(&quals)
            .add_string_tag(SamTag::RX, b"ACGT")
            .add_string_tag(SamTag::MC, b"2S124M1D3M");
        records.push(b2.build());
    }
    records
}

/// A read-through family with an indel at the overlap boundary must survive to a consensus,
/// with the positive strand clipped by the 2 bases it genuinely overhangs.
///
/// The overlap clip is computed before anything else the caller does, and the superseded
/// reference-space arithmetic answered with the positive read's *entire* 129 bases here: the
/// mate boundary landed inside its `1D`, where the read has no position at all, and the whole
/// read length was clipped. Every R1 in the family was trimmed to nothing, banked under
/// `zero_length_after_trimming`, and no R1 consensus came out — silently, since the R2 half
/// still emitted normally. Measuring the clip in query space keeps 127 of the 129 bases.
#[test]
fn test_simplex_indel_at_overlap_boundary_still_calls_a_consensus() {
    const DEPTH: usize = 3;

    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    create_grouped_bam(&input_bam, vec![("1", indel_readthrough_family(DEPTH))]);

    Simplex::try_parse_from([
        "simplex",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--min-reads",
        "3",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse simplex args")
    .execute("fgumi simplex")
    .expect("Simplex command failed");

    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let records: Vec<bam::Record> =
        reader.records().map(|r| r.expect("failed to read record")).collect();

    assert_eq!(records.len(), 2, "both strands must emit a consensus, R1 included");
    for record in &records {
        assert_eq!(
            int_tag(record, SamTag::CD),
            i64::try_from(DEPTH).expect("depth fits in i64"),
            "every raw read must contribute"
        );
    }

    // Each strand overhangs its mate by exactly 2 query bases, and those 2 come off the
    // sequencing 3' end. Pinning the full consensus sequence — not just its length — is what
    // says *which* bases were dropped: a clip of the right size taken from the wrong end, or
    // taken in the mate's coordinate frame, produces a 127-base read with different content.
    // Simplex emits consensus reads in read orientation, so the negative strand's expectation
    // is built from the reverse complement of its stored sequence.
    let (r1_bases, r2_bases) = indel_readthrough_read_bases();
    let expected_r1: Vec<u8> = r1_bases[..127].to_vec();
    // The negative strand is emitted in sequencing order, which is the reverse complement of the
    // reference-orientation bases the BAM stores; its clip comes off that sequence's 3' end.
    let expected_r2: Vec<u8> = reverse_complement(&r2_bases)[..127].to_vec();
    for record in &records {
        let bases: Vec<u8> = record.sequence().iter().collect();
        let is_r1 = record.flags().is_first_segment();
        let expected = if is_r1 { &expected_r1 } else { &expected_r2 };
        assert_eq!(
            &bases,
            expected,
            "the {} consensus must be its read minus the 2 bases it overhangs, 3'-end first",
            if is_r1 { "R1" } else { "R2" }
        );
    }
}

//////////////////////////////////////////////////////////////////////////////
// Chain-path worker-count determinism tests
//
// `simplex` always routes through the declarative chain builder (the legacy
// single-threaded fast path is retired on a `consensus` build). These tests
// diff a no-`--threads` run (the chain at a single worker) against a
// `--threads N` run (the chain at N workers): the "non-chain path" / "oracle"
// naming below now means the single-worker chain, not the retired serial loop.
// Byte-parity against the pre-removal serial binary lives in
// `test_consensus_cutover_parity.rs`.
//////////////////////////////////////////////////////////////////////////////

/// Read a BAM's records back as decoded `RecordBuf`s, for record-for-record
/// comparison (order-sensitive, catches drops/dupes/reorders that a bare count
/// would miss), via the shared `read_bam_output` helper. (The full-header half
/// of `read_bam_output` is used directly in the main parity test.)
fn read_consensus_records(path: &Path) -> Vec<noodles::sam::alignment::RecordBuf> {
    read_bam_output(path).1
}

/// Run `simplex` on `input` writing `output`, with `extra` args appended
/// (e.g. `--threads`, `--min-reads`, output flags). Asserts the run succeeds.
fn simplex_run(input: &Path, output: &Path, extra: &[&str]) {
    let mut args =
        vec!["simplex", "--input", input.to_str().unwrap(), "--output", output.to_str().unwrap()];
    args.extend_from_slice(extra);
    Simplex::try_parse_from(args)
        .expect("failed to parse simplex args")
        .execute("fgumi simplex")
        .expect("simplex run failed");
}

/// A mix of single-end and paired-end MI families, several distinct MI groups
/// each deep enough to survive `--min-reads 2` -- the default parity fixture
/// for axes that don't need a more specialized shape.
fn create_parity_families(count: usize) -> Vec<(String, Vec<fgumi_raw_bam::RawRecord>)> {
    (0..count)
        .map(|i| {
            let mi = i.to_string();
            let depth = 3 + (i % 3);
            let family = if i % 2 == 0 {
                create_umi_family("ACGT", depth, &format!("fam{i}"), "ACGTACGTAC", 30)
            } else {
                create_paired_umi_family(
                    "TGCA",
                    depth,
                    &format!("fam{i}"),
                    "ACGTACGTAC",
                    "TTTTAAAAGG",
                    30,
                )
            };
            (mi, family)
        })
        .collect()
}

/// Writes [`create_parity_families`]-shaped families as the `(&str, Vec<RawRecord>)`
/// pairs [`create_grouped_bam`] expects.
fn write_parity_bam(path: &Path, families: &[(String, Vec<fgumi_raw_bam::RawRecord>)]) {
    let refs: Vec<(&str, Vec<fgumi_raw_bam::RawRecord>)> =
        families.iter().map(|(mi, records)| (mi.as_str(), records.clone())).collect();
    create_grouped_bam(path, refs);
}

/// Axis 1 + 7: the chain (`--threads N`) path produces output records
/// record-for-record identical to the non-chain (no-`--threads`) path, at both
/// `--threads 1` (the minimal chain engine) and `--threads 4` (genuinely
/// parallel, with enough distinct MI families to cross `GroupByMi` batch
/// boundaries).
#[rstest]
#[case::threads_1(1)]
#[case::threads_4(4)]
fn test_simplex_chain_matches_single_threaded(#[case] threads: usize) {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");

    // 40 distinct MI families is enough to span multiple in-flight batches at
    // `--threads 4` without slowing the test down.
    let families = create_parity_families(40);
    write_parity_bam(&input_bam, &families);

    let oracle_out = temp_dir.path().join("oracle.bam");
    simplex_run(&input_bam, &oracle_out, &["--min-reads", "2"]);

    let chain_out = temp_dir.path().join("chain.bam");
    let threads_str = threads.to_string();
    simplex_run(&input_bam, &chain_out, &["--min-reads", "2", "--threads", &threads_str]);

    let (oracle_header, expected) = read_bam_output(&oracle_out);
    let (chain_header, actual) = read_bam_output(&chain_out);
    assert!(!expected.is_empty(), "oracle output must be non-empty (guard against a vacuous pass)");
    assert_eq!(
        actual, expected,
        "chain --threads {threads} output must match the non-chain path record-for-record"
    );
    // Whole-header parity (read_bam_output normalizes the @PG CL that differs by
    // --threads): also catches a dropped @SQ/@RG/@HD/@CO.
    assert_eq!(
        chain_header, oracle_header,
        "chain and non-chain output headers must match (with @PG CL normalized)"
    );
}

/// The `--ref requires --methylation-mode` validation (relocated to run before
/// the `--threads` chain dispatch) must reject on BOTH paths -- passing `--ref`
/// without `--methylation-mode` errors regardless of `--threads`. The bail sits
/// in the reader-free pre-flight, before the reference is opened, so a
/// nonexistent `--ref` path still hits it.
#[rstest]
#[case::single_threaded(&[])]
#[case::threaded(&["--threads", "2"])]
fn test_simplex_rejects_ref_without_methylation_mode(#[case] extra: &[&str]) {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let families = create_parity_families(4);
    write_parity_bam(&input_bam, &families);
    let out = temp_dir.path().join("out.bam");
    let fake_ref = temp_dir.path().join("ref.fa");

    let mut args = vec![
        "simplex",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
        "--min-reads",
        "1",
        "--ref",
        fake_ref.to_str().unwrap(),
    ];
    args.extend_from_slice(extra);
    let err = Simplex::try_parse_from(args)
        .expect("failed to parse simplex args")
        .execute("fgumi simplex")
        .expect_err("simplex must reject --ref without --methylation-mode");
    assert!(
        err.to_string().contains("--ref requires --methylation-mode"),
        "unexpected error: {err:#}",
    );
}

/// Axis 2: `--rejects` output from the chain path matches the non-chain path
/// record-for-record, and is non-vacuous (the singleton family must actually be
/// rejected on both paths). Also pins rejects **header** parity: the rejects
/// BAM carries the raw input header verbatim on both paths (the PR #332
/// contract), so the chain path must not leak its own `@PG` into rejects.
#[test]
fn test_simplex_chain_rejects_matches_single_threaded() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");

    let mut families = create_parity_families(10);
    // A singleton family that fails --min-reads 2 on both paths.
    families
        .push(("singleton".to_string(), create_umi_family("GATC", 1, "reject", "GGGGCCCCTT", 30)));
    write_parity_bam(&input_bam, &families);

    let oracle_out = temp_dir.path().join("oracle.bam");
    let oracle_rejects = temp_dir.path().join("oracle.rejects.bam");
    simplex_run(
        &input_bam,
        &oracle_out,
        &["--min-reads", "2", "--rejects", oracle_rejects.to_str().unwrap()],
    );

    let chain_out = temp_dir.path().join("chain.bam");
    let chain_rejects = temp_dir.path().join("chain.rejects.bam");
    simplex_run(
        &input_bam,
        &chain_out,
        &["--min-reads", "2", "--rejects", chain_rejects.to_str().unwrap(), "--threads", "4"],
    );

    let expected_out = read_consensus_records(&oracle_out);
    let actual_out = read_consensus_records(&chain_out);
    assert_eq!(actual_out, expected_out, "chain primary output must match the non-chain path");

    let (oracle_rejects_header, expected_rejects) = read_bam_output(&oracle_rejects);
    let (chain_rejects_header, actual_rejects) = read_bam_output(&chain_rejects);
    assert!(
        !expected_rejects.is_empty(),
        "oracle rejects must be non-empty (guard against a vacuous pass)"
    );
    assert_eq!(
        actual_rejects, expected_rejects,
        "chain --rejects output must match the non-chain path record-for-record"
    );

    // Rejects-header provenance: rejects are raw-input records, so both paths
    // must write them under the raw input header verbatim. `read_bam_output`
    // normalizes the `@PG` CL, so the input header (which carries no fgumi `@PG`)
    // must equal both rejects headers — a chain path that injected its own `@PG`
    // into the rejects sink would fail this.
    let (input_header, _) = read_bam_output(&input_bam);
    assert_eq!(
        chain_rejects_header, oracle_rejects_header,
        "chain and non-chain rejects headers must match"
    );
    assert_eq!(
        chain_rejects_header, input_header,
        "rejects header must be the raw input header verbatim (no injected @PG)"
    );
}

/// Axis 3: `--stats` output from the chain path is byte-identical to the
/// non-chain path (recipe: byte-compare, not just field-by-field).
#[test]
fn test_simplex_chain_stats_matches_single_threaded() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");

    let families = create_parity_families(20);
    write_parity_bam(&input_bam, &families);

    let oracle_out = temp_dir.path().join("oracle.bam");
    let oracle_stats = temp_dir.path().join("oracle.stats.txt");
    simplex_run(
        &input_bam,
        &oracle_out,
        &["--min-reads", "2", "--stats", oracle_stats.to_str().unwrap()],
    );

    let chain_out = temp_dir.path().join("chain.bam");
    let chain_stats = temp_dir.path().join("chain.stats.txt");
    simplex_run(
        &input_bam,
        &chain_out,
        &["--min-reads", "2", "--stats", chain_stats.to_str().unwrap(), "--threads", "4"],
    );

    let expected = fs::read(&oracle_stats).expect("read oracle stats");
    let actual = fs::read(&chain_stats).expect("read chain stats");
    assert!(!expected.is_empty(), "oracle stats must be non-empty (guard against a vacuous pass)");
    assert_eq!(
        actual, expected,
        "chain --stats output must be byte-identical to the non-chain path"
    );
}

/// Axis 4: methylation mode (`--methylation-mode em-seq` + `--ref`) output
/// parity between the chain and non-chain paths. Builds a minimal grouped
/// input with mapped, fully-methylated reads against a synthetic all-C
/// reference, since no dedicated methylation test infra exists in this file.
#[test]
fn test_simplex_chain_methylation_matches_single_threaded() {
    use fgumi_raw_bam::SamBuilder;
    use fgumi_sam::builder::create_test_fasta;

    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");

    // An all-C reference so every read shows a fully-methylated (unconverted) signal.
    let ref_seq = "C".repeat(200);
    let ref_fasta = create_test_fasta(&[("chr1", &ref_seq)]).expect("build test fasta");

    let header = create_minimal_header("chr1", 200);
    let seq = "CCCCCCCCCC";
    let cigar = u32::try_from(seq.len()).expect("fits u32") << 4;
    let mut families = Vec::new();
    for family_idx in 0..3 {
        let records: Vec<fgumi_raw_bam::RawRecord> = (0..3)
            .map(|read_idx| {
                let mut b = SamBuilder::new();
                b.read_name(format!("f{family_idx}_r{read_idx}").as_bytes())
                    .ref_id(0)
                    .pos(0)
                    .mapq(60)
                    .flags(0)
                    .cigar_ops(&[cigar])
                    .sequence(seq.as_bytes())
                    .qualities(&vec![30u8; seq.len()])
                    .add_string_tag(SamTag::RX, b"ACGT");
                b.build()
            })
            .collect();
        families.push((family_idx.to_string(), records));
    }
    let refs: Vec<(&str, Vec<fgumi_raw_bam::RawRecord>)> =
        families.iter().map(|(mi, records)| (mi.as_str(), records.clone())).collect();
    create_grouped_bam_with_header(&input_bam, &header, refs);

    let methylation_args = &[
        "--min-reads",
        "2",
        "--methylation-mode",
        "em-seq",
        "--ref",
        ref_fasta.path().to_str().unwrap(),
    ];

    let oracle_out = temp_dir.path().join("oracle.bam");
    simplex_run(&input_bam, &oracle_out, methylation_args);

    let mut chain_args = methylation_args.to_vec();
    chain_args.extend_from_slice(&["--threads", "4"]);
    let chain_out = temp_dir.path().join("chain.bam");
    simplex_run(&input_bam, &chain_out, &chain_args);

    let expected = read_consensus_records(&oracle_out);
    let actual = read_consensus_records(&chain_out);
    assert!(!expected.is_empty(), "oracle output must be non-empty (guard against a vacuous pass)");
    assert_eq!(
        actual, expected,
        "chain methylation-mode output must match the non-chain path record-for-record"
    );

    // Non-vacuity: the methylation tags this mode adds must actually be present.
    let cu_tag = noodles::sam::alignment::record::data::field::Tag::from([b'c', b'u']);
    assert!(
        expected.iter().all(|r| r.data().get(&cu_tag).is_some()),
        "methylation-mode consensus reads must carry the cu tag"
    );
}

/// Axis 5: `--consensus-call-overlapping-bases false` output parity between the
/// chain and non-chain paths. Uses the overlapping read-through fixture (real
/// mate overlap), so this exercises the *disabled* overlap path -- the default
/// (enabled) case is already covered by [`test_simplex_chain_matches_single_threaded`].
#[test]
fn test_simplex_chain_overlapping_disabled_matches_single_threaded() {
    const DEPTH: usize = 3;

    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    create_grouped_bam(&input_bam, vec![("1", indel_readthrough_family(DEPTH))]);

    let args = &["--min-reads", "3", "--consensus-call-overlapping-bases", "false"];

    let oracle_out = temp_dir.path().join("oracle.bam");
    simplex_run(&input_bam, &oracle_out, args);

    let mut chain_args = args.to_vec();
    chain_args.extend_from_slice(&["--threads", "4"]);
    let chain_out = temp_dir.path().join("chain.bam");
    simplex_run(&input_bam, &chain_out, &chain_args);

    let expected = read_consensus_records(&oracle_out);
    let actual = read_consensus_records(&chain_out);
    assert!(!expected.is_empty(), "oracle output must be non-empty (guard against a vacuous pass)");
    assert_eq!(
        actual, expected,
        "chain --consensus-call-overlapping-bases false output must match the non-chain path"
    );
}

/// Axis 6: `--allow-unmapped` output parity between the chain and non-chain
/// paths, using a wholly-unmapped family (the case that most exercises the
/// flag -- see `test_simplex_allow_unmapped_still_consenses_a_wholly_unmapped_family`).
#[test]
fn test_simplex_chain_allow_unmapped_matches_single_threaded() {
    use fgumi_raw_bam::{SamBuilder, flags};

    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");

    let mut records = Vec::new();
    for i in 0..3 {
        for (segment, seq) in
            [(flags::FIRST_SEGMENT, "ACGTACGTAC"), (flags::LAST_SEGMENT, "GGGGGGGGGG")]
        {
            let mut b = SamBuilder::new();
            b.read_name(format!("unmapped_{i}").as_bytes())
                .ref_id(-1)
                .pos(-1)
                .mapq(0)
                .flags(flags::PAIRED | segment | flags::UNMAPPED | flags::MATE_UNMAPPED)
                .mate_ref_id(-1)
                .mate_pos(-1)
                .sequence(seq.as_bytes())
                .qualities(&vec![30u8; seq.len()])
                .add_string_tag(SamTag::RX, b"ACGT");
            records.push(b.build());
        }
    }
    create_grouped_bam(&input_bam, vec![("1", records)]);

    let args = &["--min-reads", "1", "--allow-unmapped"];

    let oracle_out = temp_dir.path().join("oracle.bam");
    simplex_run(&input_bam, &oracle_out, args);

    let mut chain_args = args.to_vec();
    chain_args.extend_from_slice(&["--threads", "4"]);
    let chain_out = temp_dir.path().join("chain.bam");
    simplex_run(&input_bam, &chain_out, &chain_args);

    let expected = read_consensus_records(&oracle_out);
    let actual = read_consensus_records(&chain_out);
    assert!(!expected.is_empty(), "oracle output must be non-empty (guard against a vacuous pass)");
    assert_eq!(
        actual, expected,
        "chain --allow-unmapped output must match the non-chain path record-for-record"
    );
}

/// Axis 8 (proves Task 2A): a coordinate-sorted input is REJECTED on the chain
/// (`--threads 4`) path with the same "not sorted correctly for consensus
/// calling" error the non-chain path already raises. `add_simplex` has no
/// out-of-order/duplicate detection of its own (`GroupByMi` doesn't either), so
/// without the `check_consensus_sort_order` guard a mis-sorted input would be
/// silently mis-grouped instead of rejected.
#[test]
fn test_simplex_chain_rejects_coordinate_sorted_input() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");

    let header = create_coordinate_sorted_header("chr1", 10000);
    let family = create_umi_family("ACGT", 3, "fam1", "ACGTACGTAC", 30);
    create_grouped_bam_with_header(&input_bam, &header, vec![("1", family)]);

    let output_bam = temp_dir.path().join("output.bam");

    // Non-chain path (no --threads): must already reject.
    let oracle_err = Simplex::try_parse_from([
        "simplex",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--min-reads",
        "2",
    ])
    .expect("failed to parse simplex args")
    .execute("fgumi simplex")
    .expect_err("non-chain path must reject coordinate-sorted input");
    let oracle_message = format!("{oracle_err:#}");
    assert!(
        oracle_message.contains("not sorted correctly for consensus calling"),
        "non-chain error should mention sort order: {oracle_message}"
    );

    // Chain path (--threads 4): must reject with the same message, proving the
    // Task 2A `check_consensus_sort_order` guard in `add_simplex`.
    let chain_out = temp_dir.path().join("chain_output.bam");
    let chain_err = Simplex::try_parse_from([
        "simplex",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        chain_out.to_str().unwrap(),
        "--min-reads",
        "2",
        "--threads",
        "4",
    ])
    .expect("failed to parse simplex args")
    .execute("fgumi simplex")
    .expect_err("chain path must reject coordinate-sorted input");
    let chain_message = format!("{chain_err:#}");
    assert!(
        chain_message.contains("not sorted correctly for consensus calling"),
        "chain error should mention sort order: {chain_message}"
    );
    assert!(!chain_out.exists(), "chain path must not create an output file on rejection");
}

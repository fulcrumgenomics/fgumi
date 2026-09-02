//! End-to-end in-process tests for the codec command.
//!
//! These tests invoke the `Codec` command's `execute()` method in-process and validate:
//! 1. Basic consensus calling from CODEC read pairs
//! 2. Statistics output
//! 3. Rejected reads output
//! 4. Quality filtering options

use clap::Parser;
use fgumi_bam_io::create_raw_bam_writer;
use fgumi_dna::reverse_complement;
use fgumi_lib::commands::codec::{Codec, CodecOptions};
use fgumi_lib::commands::command::Command;
use fgumi_lib::pipeline::chains::{
    ChainSpec, SingleStageContext, Stage, StageOptionsBag, build_for,
};
use fgumi_lib::sam::SamTag;
use fgumi_raw_bam::{RawRecord, RawRecordView, SamBuilder, flags as raw_flags};
use noodles::bam;
use rstest::rstest;
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use tempfile::TempDir;

use crate::helpers::assertions::{int_tag, string_tag};
use crate::helpers::bam_generator::{create_coordinate_sorted_header, create_minimal_header};
use crate::helpers::read_bam_output;

/// Creates a CODEC read pair (R1 forward, R2 reverse from opposite strand).
///
/// In CODEC sequencing, R1 and R2 come from opposite strands of the same molecule,
/// so a single read pair can produce duplex consensus.
#[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap, clippy::too_many_arguments)]
pub(crate) fn create_codec_read_pair(
    name: &str,
    r1_seq: &[u8],
    r2_seq: &[u8],
    r1_qual: &[u8],
    r2_qual: &[u8],
    ref_start: usize,
    umi: &str,
    cell_barcode: Option<&str>,
) -> (RawRecord, RawRecord) {
    let r1_len = r1_seq.len();
    let r2_len = r2_seq.len();
    let r1_cigar_op = u32::try_from(r1_len).expect("r1_len fits u32") << 4; // nM
    let r2_cigar_op = u32::try_from(r2_len).expect("r2_len fits u32") << 4; // nM
    // SamBuilder pos is 0-based; ref_start is 1-based
    let pos = i32::try_from(ref_start).expect("ref_start fits i32") - 1;
    // MC is the mate's CIGAR, so b1 carries R2's tag and vice versa.
    let r1_mc_tag = format!("{r1_len}M");
    let r2_mc_tag = format!("{r2_len}M");

    let mut b1 = SamBuilder::new();
    b1.read_name(name.as_bytes())
        .sequence(r1_seq)
        .qualities(r1_qual)
        .cigar_ops(&[r1_cigar_op])
        .flags(raw_flags::PAIRED | raw_flags::FIRST_SEGMENT | raw_flags::MATE_REVERSE)
        .ref_id(0)
        .pos(pos)
        .mapq(60)
        .mate_ref_id(0)
        .mate_pos(pos)
        .template_length(r1_len as i32)
        .add_string_tag(SamTag::MI, umi.as_bytes())
        .add_string_tag(SamTag::MC, r2_mc_tag.as_bytes());
    if let Some(cb) = cell_barcode {
        b1.add_string_tag(SamTag::CB, cb.as_bytes());
    }

    let mut b2 = SamBuilder::new();
    b2.read_name(name.as_bytes())
        .sequence(r2_seq)
        .qualities(r2_qual)
        .cigar_ops(&[r2_cigar_op])
        .flags(raw_flags::PAIRED | raw_flags::LAST_SEGMENT | raw_flags::REVERSE)
        .ref_id(0)
        .pos(pos)
        .mapq(60)
        .mate_ref_id(0)
        .mate_pos(pos)
        .template_length(-(r2_len as i32))
        .add_string_tag(SamTag::MI, umi.as_bytes())
        .add_string_tag(SamTag::MC, r1_mc_tag.as_bytes());
    if let Some(cb) = cell_barcode {
        b2.add_string_tag(SamTag::CB, cb.as_bytes());
    }

    (b1.build(), b2.build())
}

/// Helper to create a test BAM file with CODEC read pairs.
pub(crate) fn create_codec_test_bam(path: &PathBuf, pairs: Vec<(RawRecord, RawRecord)>) {
    let header = create_minimal_header("chr1", 10000);
    let mut writer =
        create_raw_bam_writer(path, &header, 1, 6).expect("Failed to create raw BAM writer");
    for (r1, r2) in pairs {
        writer.write_raw_record(r1.as_ref()).expect("Failed to write R1");
        writer.write_raw_record(r2.as_ref()).expect("Failed to write R2");
    }
    writer.finish().expect("Failed to finish BAM");
}

/// Test basic CODEC consensus calling.
#[test]
fn test_codec_command_basic_consensus() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    // Create 3 read pairs for one molecule
    let mut pairs = Vec::new();
    for i in 0..3 {
        let (r1, r2) = create_codec_read_pair(
            &format!("read{i}"),
            b"ACGTACGT",
            b"ACGTACGT",
            &[30; 8],
            &[30; 8],
            100,
            "UMI001",
            None,
        );
        pairs.push((r1, r2));
    }
    create_codec_test_bam(&input_bam, pairs);

    // Run codec command
    let cmd = Codec::try_parse_from([
        "codec",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--min-reads",
        "1",
        "--min-duplex-length",
        "1",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse codec args");
    cmd.execute("fgumi codec").expect("Failed to run codec command");
    assert!(output_bam.exists(), "Output BAM not created");

    // Read output and verify consensus was created
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let mut consensus_count = 0;

    for result in reader.records() {
        let record = result.expect("Failed to read record");
        consensus_count += 1;

        // Verify consensus tags exist by checking the raw tag bytes
        let cd_tag = SamTag::CD.to_noodles_tag();
        assert!(record.data().get(&cd_tag).is_some(), "Consensus should have cD tag");
    }

    assert!(consensus_count > 0, "Should have produced at least one consensus read");
}

/// CODEC3-05 (end-to-end): `--single-strand-qual` must mask the single-strand tails, which are
/// padded with lowercase `n`. The offset molecule spans ref[0..40) with 10 single-strand
/// positions on each end. Asserts the *complete* emitted output — record count, all consensus
/// bases, and the full per-base quality vector — for both option states against fixed
/// expectations (pinned to fgumi's output, externally verified byte-identical to fgbio 4.1.0),
/// so wrong bases, drifted interior qualities, or an extra/dropped record all fail the test.
/// Before the fix the option had no effect on the tails (the mask check missed lowercase `n`).
#[test]
fn test_codec_command_masks_single_strand_tails() {
    // The two same-UMI molecules collapse to a single 40bp consensus. The overlap window
    // (ref[10..30]) disagrees between the two synthetic strands, so its bases are called `N`; the
    // single-strand tails (ref[0..10] from R1, ref[30..40] from R2) carry each read's base. Only
    // the per-base qualities change under masking. All values are pinned to what fgumi emits
    // (externally verified byte-identical to fgbio 4.1.0), so the test fails on any drift in
    // bases, interior qualities, or record count.
    const EXPECTED_BASES: &[u8; 40] = b"AAAAAAAAAANNNNNNNNNNNNNNNNNNNNGGGGGGGGGG";
    // ref[10..30] is the 20bp overlap (single-strand-disagreement `N` calls); tails (ref[0..10],
    // ref[30..40]) are the lowercase-'n'-padded single-strand positions the option masks.
    const OVERLAP_QUALS: [u8; 20] = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];
    // Without --single-strand-qual the single-strand tails keep their unmasked Q45.
    let expected_baseline_quals: Vec<u8> = {
        let mut q = vec![45u8; 40];
        q[10..30].copy_from_slice(&OVERLAP_QUALS);
        q
    };
    // With --single-strand-qual 2 the single-strand tails are masked to Q2; the overlap
    // (already Q2 from the disagreement calls) is untouched.
    let expected_masked_quals: Vec<u8> = {
        let mut q = vec![2u8; 40];
        q[10..30].copy_from_slice(&OVERLAP_QUALS);
        q
    };

    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    create_codec_test_bam(
        &input_bam,
        vec![create_offset_codec_pair("m1", "UMI1"), create_offset_codec_pair("m2", "UMI1")],
    );

    // Run the codec command with the given single-strand-qual (None = option omitted) and return
    // every emitted consensus record's (bases, qualities). Reading *all* records (not just the
    // first) guards against extra/dropped records; returning full base+quality vectors guards
    // against wrong bases or interior qualities that a tail-only threshold check would miss.
    let run_codec = |ssq: Option<&str>| -> Vec<(Vec<u8>, Vec<u8>)> {
        let output_bam = temp_dir.path().join(format!("out_{}.bam", ssq.unwrap_or("none")));
        let mut args = vec![
            "codec",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_bam.to_str().unwrap(),
            "--min-reads",
            "1",
            "--min-duplex-length",
            "1",
            "--max-duplex-disagreements",
            "1000",
            "--max-duplex-disagreement-rate",
            "1.0",
            "--compression-level",
            "1",
        ];
        if let Some(q) = ssq {
            args.push("--single-strand-qual");
            args.push(q);
        }
        Codec::try_parse_from(args)
            .expect("failed to parse codec args")
            .execute("fgumi codec")
            .expect("Failed to run codec command");
        let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
        let _header = reader.read_header().unwrap();
        reader
            .records()
            .map(|r| {
                let record = r.expect("read record");
                let bases = record.sequence().iter().collect::<Vec<u8>>();
                let quals = record.quality_scores().as_ref().to_vec();
                (bases, quals)
            })
            .collect()
    };

    let baseline = run_codec(None);
    let masked = run_codec(Some("2"));

    assert_eq!(baseline.len(), 1, "exactly one consensus record without masking");
    assert_eq!(masked.len(), 1, "exactly one consensus record with masking");

    let (baseline_bases, baseline_quals) = &baseline[0];
    let (masked_bases, masked_quals) = &masked[0];

    // Masking must not perturb the called bases.
    assert_eq!(baseline_bases, EXPECTED_BASES, "unmasked consensus bases");
    assert_eq!(masked_bases, EXPECTED_BASES, "masked consensus bases");

    // Full quality-vector identity for both states — before the fix the tails stayed at Q45 even
    // with --single-strand-qual set (the mask check missed the lowercase-'n' padding).
    assert_eq!(baseline_quals, &expected_baseline_quals, "unmasked consensus qualities");
    assert_eq!(masked_quals, &expected_masked_quals, "masked consensus qualities");
}

/// Test CODEC command with statistics output.
#[test]
fn test_codec_command_with_stats() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let stats_file = temp_dir.path().join("stats.tsv");

    // Create read pairs for two molecules
    let mut pairs = Vec::new();
    for i in 0..2 {
        let (r1, r2) = create_codec_read_pair(
            &format!("mol1_read{i}"),
            b"ACGTACGT",
            b"ACGTACGT",
            &[30; 8],
            &[30; 8],
            100,
            "UMI001",
            None,
        );
        pairs.push((r1, r2));
    }
    for i in 0..3 {
        let (r1, r2) = create_codec_read_pair(
            &format!("mol2_read{i}"),
            b"TGCATGCA",
            b"TGCATGCA",
            &[30; 8],
            &[30; 8],
            200,
            "UMI002",
            None,
        );
        pairs.push((r1, r2));
    }
    create_codec_test_bam(&input_bam, pairs);

    // Run codec command with stats
    let cmd = Codec::try_parse_from([
        "codec",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--stats",
        stats_file.to_str().unwrap(),
        "--min-reads",
        "1",
        "--min-duplex-length",
        "1",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse codec args");
    cmd.execute("fgumi codec").expect("Failed to run codec command");
    assert!(stats_file.exists(), "Stats file not created");

    // Verify stats file has content
    let stats_content = fs::read_to_string(&stats_file).expect("Failed to read stats");
    assert!(!stats_content.is_empty(), "Stats file should not be empty");
}

/// A codec consensus record reduced to the fields a caller cannot get right by accident:
/// its name, flags, sequence, molecule identifier, and the raw-read depths it claims.
#[derive(Debug, PartialEq, Eq)]
struct CodecConsensusIdentity {
    name: String,
    flags: u16,
    sequence: String,
    mi: String,
    total_depth: i64,
    a_depth: i64,
    b_depth: i64,
}

impl CodecConsensusIdentity {
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
            total_depth: int_tag(record, SamTag::CD),
            a_depth: int_tag(record, SamTag::AD),
            b_depth: int_tag(record, SamTag::BD),
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
fn test_codec_stats_consensus_reads_matches_records_written(#[case] threads: Option<&str>) {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let stats_file = temp_dir.path().join("stats.tsv");

    let mut pairs = Vec::new();
    for molecule in 0..3 {
        for read in 0..3 {
            let (r1, r2) = create_codec_read_pair(
                &format!("mol{molecule}_read{read}"),
                b"ACGTACGT",
                b"ACGTACGT",
                &[30; 8],
                &[30; 8],
                100 + molecule * 100,
                &format!("UMI00{molecule}"),
                None,
            );
            pairs.push((r1, r2));
        }
    }
    create_codec_test_bam(&input_bam, pairs);

    let mut args = vec![
        "codec",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--stats",
        stats_file.to_str().unwrap(),
        "--min-reads",
        "1",
        "--min-duplex-length",
        "1",
        "--compression-level",
        "1",
    ];
    if let Some(threads) = threads {
        args.extend_from_slice(&["--threads", threads]);
    }
    Codec::try_parse_from(args)
        .expect("failed to parse codec args")
        .execute("fgumi codec")
        .expect("Failed to run codec command");

    // Identity, not just cardinality: every molecule is built from three identical read
    // pairs, so each consensus record is fully determined -- one unmapped record per
    // molecule, named and MI-tagged by its UMI, carrying the molecule's sequence and its
    // six raw reads split evenly across the two strands.
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let records =
        reader.records().collect::<Result<Vec<_>, _>>().expect("failed to read the consensus BAM");
    let observed = records.iter().map(CodecConsensusIdentity::from_record).collect::<Vec<_>>();
    let expected = (0..3)
        .map(|molecule| CodecConsensusIdentity {
            name: format!(":UMI00{molecule}"),
            flags: raw_flags::UNMAPPED,
            sequence: "ACGTACGT".to_string(),
            mi: format!("UMI00{molecule}"),
            total_depth: 6,
            a_depth: 3,
            b_depth: 3,
        })
        .collect::<Vec<_>>();
    assert_eq!(observed, expected, "codec must emit one fully determined consensus per molecule");

    let records_written = records.len();

    let stats = fs::read_to_string(&stats_file).expect("Failed to read stats");
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

/// #748: `codec --stats` must emit the five codec-only rows fgbio's
/// `CallCodecConsensusReads --stats` writes, with fgbio's exact key names and
/// description strings, so the two tools' stats files stay diffable. All five were
/// absent, which made a diff report missing keys rather than value differences.
///
/// The counts are checked against the output BAM rather than hard-coded, so the
/// rows are pinned to what the run actually emitted.
///
/// The `with_rejected_molecule` case adds a molecule whose two strands disagree at
/// every base, rejected by `--max-duplex-disagreements 0`. It pins the "emitted"
/// half of the contract: a rejected molecule bumps `consensus_reads_rejected_hdd`
/// and contributes to *none* of the three base counters. Counting it would inflate
/// `duplex_disagreement_rate` by exactly the molecules discarded for disagreeing.
#[rstest]
#[case::all_molecules_emitted(false)]
#[case::with_rejected_molecule(true)]
fn test_codec_stats_emits_fgbio_codec_only_rows(#[case] add_disagreeing_molecule: bool) {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let stats_file = temp_dir.path().join("stats.tsv");

    let mut pairs = Vec::new();
    for molecule in 0..3 {
        for read in 0..3 {
            let (r1, r2) = create_codec_read_pair(
                &format!("mol{molecule}_read{read}"),
                b"ACGTACGT",
                b"ACGTACGT",
                &[30; 8],
                &[30; 8],
                100 + molecule * 100,
                &format!("UMI00{molecule}"),
                None,
            );
            pairs.push((r1, r2));
        }
    }
    if add_disagreeing_molecule {
        // Same shape as the molecules above, but the second strand carries a
        // different base at every position, so the molecule is rejected.
        for read in 0..3 {
            let (r1, r2) = create_codec_read_pair(
                &format!("hdd_read{read}"),
                b"ACGTACGT",
                b"TTTTTTTT",
                &[30; 8],
                &[30; 8],
                400,
                "UMI_HDD",
                None,
            );
            pairs.push((r1, r2));
        }
    }
    create_codec_test_bam(&input_bam, pairs);

    Codec::try_parse_from(vec![
        "codec",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--stats",
        stats_file.to_str().unwrap(),
        "--min-reads",
        "1",
        "--min-duplex-length",
        "1",
        "--max-duplex-disagreements",
        "0",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse codec args")
    .execute("fgumi codec")
    .expect("Failed to run codec command");

    // The three agreeing molecules are always emitted; the disagreeing one never is,
    // so the emitted-base totals are identical across both cases.
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let records =
        reader.records().collect::<Result<Vec<_>, _>>().expect("failed to read the consensus BAM");
    assert_eq!(records.len(), 3, "only the three agreeing molecules are consensus-called");
    let bases_written: usize = records.iter().map(|record| record.sequence().len()).sum();
    assert!(bases_written > 0, "the run must emit at least one consensus base");

    let stats = fs::read_to_string(&stats_file).expect("Failed to read stats");
    let row = |key: &str| -> (String, String) {
        let mut fields = stats
            .lines()
            .find(|line| line.starts_with(&format!("{key}\t")))
            .unwrap_or_else(|| panic!("stats must contain a `{key}` row"))
            .split('\t');
        let _key = fields.next();
        let value = fields.next().expect("row must have a value column").to_string();
        let description = fields.next().expect("row must have a description column").to_string();
        (value, description)
    };

    // fgbio's descriptions, verbatim (`CodecConsensusCaller.statistics`).
    let expected_rejected = usize::from(add_disagreeing_molecule);
    assert_eq!(
        row("consensus_reads_rejected_hdd"),
        (
            expected_rejected.to_string(),
            "Consensus Reads Rejected: High Duplex Disagreement".to_string()
        ),
        "one molecule per disagreeing template is counted, once"
    );
    assert_eq!(
        row("consensus_bases_emitted"),
        (bases_written.to_string(), "Total consensus bases emitted in consensus reads".to_string()),
        "consensus_bases_emitted must total the bases written to the output BAM"
    );
    assert_eq!(
        row("consensus_duplex_bases_emitted"),
        (
            bases_written.to_string(),
            "Consensus bases emitted with support from both strands of the duplex".to_string()
        ),
        "R1 and R2 fully overlap in the emitted molecules, and the rejected one emits nothing"
    );
    assert_eq!(
        row("duplex_disagreement_base_count"),
        (
            "0".to_string(),
            "Number of consensus bases at which the top and bottom strands disagreed".to_string()
        ),
        "the emitted molecules' strands agree; the rejected molecule's disagreements do not count"
    );
    assert_eq!(
        row("duplex_disagreement_rate"),
        (
            "0.000000".to_string(),
            "Rate of top/bottom strand disagreement within duplex regions of consensus reads"
                .to_string()
        ),
        "the rate covers emitted output only, so a rejected molecule cannot inflate it"
    );
}

/// Test CODEC command with rejected reads output.
#[test]
fn test_codec_command_with_rejects() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let rejects_bam = temp_dir.path().join("rejects.bam");

    // Create one pair that will pass and one that won't (need min-reads=3)
    let mut pairs = Vec::new();

    // Molecule 1: 3 pairs (will pass with min-reads=3)
    for i in 0..3 {
        let (r1, r2) = create_codec_read_pair(
            &format!("pass_read{i}"),
            b"ACGTACGT",
            b"ACGTACGT",
            &[30; 8],
            &[30; 8],
            100,
            "UMI_PASS",
            None,
        );
        pairs.push((r1, r2));
    }

    // Molecule 2: 1 pair (will fail with min-reads=3)
    let (r1, r2) = create_codec_read_pair(
        "fail_read0",
        b"TGCATGCA",
        b"TGCATGCA",
        &[30; 8],
        &[30; 8],
        200,
        "UMI_FAIL",
        None,
    );
    pairs.push((r1, r2));

    create_codec_test_bam(&input_bam, pairs);

    // Run codec command with rejects output and high min-reads
    let cmd = Codec::try_parse_from([
        "codec",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--rejects",
        rejects_bam.to_str().unwrap(),
        "--min-reads",
        "3",
        "--min-duplex-length",
        "1",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse codec args");
    cmd.execute("fgumi codec").expect("Failed to run codec command");
    assert!(output_bam.exists(), "Output BAM not created");
    assert!(rejects_bam.exists(), "Rejects BAM not created");

    // Count records in each file
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let consensus_count = reader.records().count();

    // Should have at least one consensus (from 3-pair molecule)
    assert!(consensus_count >= 1, "Should have consensus from passing molecule");

    // The single-pair UMI_FAIL molecule (1 < min-reads=3) should be written to
    // the rejects BAM as 2 records (R1 + R2). Match the assertion shape used by
    // test_simplex_command_with_rejects / test_duplex_command_with_rejects.
    let mut rejects_reader = bam::io::Reader::new(fs::File::open(&rejects_bam).unwrap());
    let _rejects_header = rejects_reader.read_header().unwrap();
    let rejects_count = rejects_reader.records().count();
    assert_eq!(rejects_count, 2, "Rejects BAM should contain both reads of the failing molecule");
}

/// Builds a CODEC FR pair carrying an explicit CIGAR on both reads.
///
/// [`create_codec_read_pair`] derives an all-`M` CIGAR from the sequence length, so it
/// cannot express the differing indel patterns the minority-alignment filter keys on.
#[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
fn create_codec_read_pair_with_cigar(
    name: &str,
    seq: &[u8],
    ref_start: usize,
    umi: &str,
    cigar_ops: &[u32],
) -> (RawRecord, RawRecord) {
    let pos = i32::try_from(ref_start).expect("ref_start fits i32") - 1;
    let qual = vec![30u8; seq.len()];

    let mut b1 = SamBuilder::new();
    b1.read_name(name.as_bytes())
        .sequence(seq)
        .qualities(&qual)
        .cigar_ops(cigar_ops)
        .flags(raw_flags::PAIRED | raw_flags::FIRST_SEGMENT | raw_flags::MATE_REVERSE)
        .ref_id(0)
        .pos(pos)
        .mapq(60)
        .mate_ref_id(0)
        .mate_pos(pos)
        .template_length(seq.len() as i32)
        .add_string_tag(SamTag::MI, umi.as_bytes());

    let mut b2 = SamBuilder::new();
    b2.read_name(name.as_bytes())
        .sequence(seq)
        .qualities(&qual)
        .cigar_ops(cigar_ops)
        .flags(raw_flags::PAIRED | raw_flags::LAST_SEGMENT | raw_flags::REVERSE)
        .ref_id(0)
        .pos(pos)
        .mapq(60)
        .mate_ref_id(0)
        .mate_pos(pos)
        .template_length(-(seq.len() as i32))
        .add_string_tag(SamTag::MI, umi.as_bytes());

    (b1.build(), b2.build())
}

/// #751: the rejects BAM must reconcile with `--stats`.
///
/// The molecule here emits a consensus *and* drops reads: `minority` carries a minority
/// indel pattern (`raw_reads_rejected_for_minority_alignment`) and `fragment` is unpaired
/// (`raw_reads_rejected_for_non_paired_reads`). Both rejections happen inside a group
/// that goes on to produce a consensus, which is exactly the case that used to be counted
/// in `--stats` and written nowhere — `raw_reads_rejected` was non-zero against a rejects
/// BAM of zero records. Asserted in both the single-threaded and pipeline paths, since
/// they drain the caller's rejects through different code.
#[rstest]
#[case::single_threaded(None)]
#[case::threaded(Some("2"))]
fn test_codec_rejects_bam_reconciles_with_stats(#[case] threads: Option<&str>) {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let rejects_bam = temp_dir.path().join("rejects.bam");
    let stats_file = temp_dir.path().join("stats.tsv");

    // Three majority (30M) templates and one minority (10M2D20M) template, all in one
    // molecule, so the majority reads still consense while the minority pair is dropped.
    let ungapped = [30u32 << 4];
    let gapped = [10u32 << 4, (2u32 << 4) | 2, 20u32 << 4];
    let seq = [b'A'; 30];
    let mut pairs = Vec::new();
    for i in 0..3 {
        pairs.push(create_codec_read_pair_with_cigar(
            &format!("majority{i}"),
            &seq,
            1,
            "UMI_MIX",
            &ungapped,
        ));
    }
    pairs.push(create_codec_read_pair_with_cigar("minority", &seq, 1, "UMI_MIX", &gapped));

    // An unpaired read in the same molecule: rejected on its own while the FR pairs consense.
    let mut fragment = SamBuilder::new();
    fragment
        .read_name(b"fragment")
        .sequence(&seq)
        .qualities(&[30u8; 30])
        .cigar_ops(&ungapped)
        .flags(0)
        .ref_id(0)
        .pos(0)
        .mapq(60)
        .add_string_tag(SamTag::MI, b"UMI_MIX");

    let mut input_records: Vec<RawRecord> = Vec::new();
    for (r1, r2) in pairs {
        input_records.push(r1);
        input_records.push(r2);
    }
    input_records.push(fragment.build());

    // The --rejects contract is byte-for-byte preservation of the input record, so key the
    // three expected rejects by (name, flags) and compare whole serialized records below.
    let expected_rejects: HashMap<(Vec<u8>, u16), Vec<u8>> = input_records
        .iter()
        .filter(|record| {
            matches!(RawRecordView::new(record).read_name(), b"minority" | b"fragment")
        })
        .map(|record| {
            let view = RawRecordView::new(record);
            ((view.read_name().to_vec(), view.flags()), record.as_ref().to_vec())
        })
        .collect();
    assert_eq!(expected_rejects.len(), 3, "the minority pair and the fragment must be rejected");

    let header = create_minimal_header("chr1", 10000);
    let mut writer =
        create_raw_bam_writer(&input_bam, &header, 1, 6).expect("Failed to create raw BAM writer");
    for record in &input_records {
        writer.write_raw_record(record.as_ref()).expect("Failed to write record");
    }
    writer.finish().expect("Failed to finish BAM");

    let mut args = vec![
        "codec",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--rejects",
        rejects_bam.to_str().unwrap(),
        "--stats",
        stats_file.to_str().unwrap(),
        "--min-reads",
        "1",
        "--min-duplex-length",
        "1",
        "--compression-level",
        "1",
    ];
    if let Some(threads) = threads {
        args.extend_from_slice(&["--threads", threads]);
    }
    Codec::try_parse_from(args)
        .expect("failed to parse codec args")
        .execute("fgumi codec")
        .expect("Failed to run codec command");

    let mut output_reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _output_header = output_reader.read_header().unwrap();
    assert_eq!(
        output_reader.records().count(),
        1,
        "the majority-alignment reads must still emit a consensus"
    );

    let stats = fs::read_to_string(&stats_file).expect("Failed to read stats");
    let stat = |key: &str| -> usize {
        stats
            .lines()
            .find_map(|line| line.strip_prefix(&format!("{key}\t")))
            .and_then(|rest| rest.split('\t').next())
            .unwrap_or_else(|| panic!("stats must contain a {key} row"))
            .parse()
            .unwrap_or_else(|_| panic!("{key} must be an integer"))
    };
    assert_eq!(
        stat("raw_reads_rejected_for_minority_alignment"),
        2,
        "the minority template's R1 and R2 must be rejected"
    );
    assert_eq!(
        stat("raw_reads_rejected_for_non_paired_reads"),
        1,
        "the fragment read must be rejected"
    );

    // Read the rejects back as raw records so the whole serialized record is compared, not
    // just the fields noodles decodes ergonomically: a rejects BAM that dropped tags or
    // rewrote a CIGAR would satisfy a name-and-count assertion.
    let (mut rejects_reader, _rejects_header) =
        fgumi_bam_io::create_raw_bam_reader(&rejects_bam, 1).expect("open rejects BAM");
    let mut seen: Vec<(Vec<u8>, u16)> = Vec::new();
    let mut record = RawRecord::new();
    while rejects_reader.read_record(&mut record).expect("read reject record") != 0 {
        let key = (RawRecordView::new(&record).read_name().to_vec(), record.flags());
        let expected = expected_rejects.get(&key).unwrap_or_else(|| {
            panic!("unexpected reject record {}", String::from_utf8_lossy(&key.0))
        });
        assert_eq!(
            record.as_ref(),
            expected.as_slice(),
            "each reject must be byte-for-byte identical to its input record — every field \
             and tag preserved on the --rejects path"
        );
        assert!(!seen.contains(&key), "a reject record was written more than once");
        seen.push(key);
    }

    assert_eq!(
        seen.len(),
        expected_rejects.len(),
        "the rejects BAM must hold exactly the reads counted as rejected"
    );
    assert_eq!(
        seen.len(),
        stat("raw_reads_rejected"),
        "the rejects BAM record count must equal raw_reads_rejected"
    );

    // The --rejects path preserves input order: `drain_marked_rejects` walks each group's
    // records in input order (see `CodecConsensusCaller::drain_marked_rejects`), and this
    // molecule is a single group, so both the single-threaded and pipeline paths emit the
    // rejects in the order they appeared in the input. Assert the exact ordered sequence —
    // including both the minority-alignment pair and the fragment reject — so a reordering
    // regression on the mixed-rejection path cannot pass on membership and count alone.
    let expected_reject_order: Vec<(Vec<u8>, u16)> = input_records
        .iter()
        .filter(|record| {
            matches!(RawRecordView::new(record).read_name(), b"minority" | b"fragment")
        })
        .map(|record| {
            let view = RawRecordView::new(record);
            (view.read_name().to_vec(), view.flags())
        })
        .collect();
    assert_eq!(
        seen, expected_reject_order,
        "the rejects BAM must emit the rejected reads in their original input order"
    );
}

/// Test CODEC command with minimum duplex length filter.
#[test]
fn test_codec_command_min_duplex_length() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    // Create read pairs with short sequences (8bp)
    let mut pairs = Vec::new();
    for i in 0..3 {
        let (r1, r2) = create_codec_read_pair(
            &format!("read{i}"),
            b"ACGTACGT",
            b"ACGTACGT",
            &[30; 8],
            &[30; 8],
            100,
            "UMI001",
            None,
        );
        pairs.push((r1, r2));
    }
    create_codec_test_bam(&input_bam, pairs);

    // Run codec command with high min-duplex-length (should reject)
    let cmd = Codec::try_parse_from([
        "codec",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--min-reads",
        "1",
        "--min-duplex-length",
        "100", // Much longer than our 8bp reads
        "--compression-level",
        "1",
    ])
    .expect("failed to parse codec args");
    cmd.execute("fgumi codec").expect("Failed to run codec command");

    // Should produce no consensus due to insufficient duplex length
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let consensus_count = reader.records().count();

    assert_eq!(consensus_count, 0, "Should have no consensus due to min-duplex-length filter");
}

/// Test CODEC command with per-base tags output.
#[test]
fn test_codec_command_per_base_tags() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    // Create read pairs
    let mut pairs = Vec::new();
    for i in 0..3 {
        let (r1, r2) = create_codec_read_pair(
            &format!("read{i}"),
            b"ACGTACGT",
            b"ACGTACGT",
            &[30; 8],
            &[30; 8],
            100,
            "UMI001",
            None,
        );
        pairs.push((r1, r2));
    }
    create_codec_test_bam(&input_bam, pairs);

    // Run codec command with per-base tags
    let cmd = Codec::try_parse_from([
        "codec",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--min-reads",
        "1",
        "--min-duplex-length",
        "1",
        "--output-per-base-tags",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse codec args");
    cmd.execute("fgumi codec").expect("Failed to run codec command");

    // Verify per-base tags exist
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();

    for result in reader.records() {
        let record = result.expect("Failed to read record");

        // Check for per-base depth tags (ad, bd)
        let ad_tag = SamTag::AD_BASES.to_noodles_tag();
        let bd_tag = SamTag::BD_BASES.to_noodles_tag();

        assert!(record.data().get(&ad_tag).is_some(), "Should have per-base depth tag 'ad'");
        assert!(record.data().get(&bd_tag).is_some(), "Should have per-base depth tag 'bd'");
    }
}

/// Test CODEC command preserves cell barcode tag.
#[test]
fn test_codec_command_cell_barcode_preservation() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    // Create 3 read pairs with cell barcode
    let mut pairs = Vec::new();
    for i in 0..3 {
        let (r1, r2) = create_codec_read_pair(
            &format!("read{i}"),
            b"ACGTACGT",
            b"ACGTACGT",
            &[30; 8],
            &[30; 8],
            100,
            "UMI001",
            Some("CELLBC123"),
        );
        pairs.push((r1, r2));
    }
    create_codec_test_bam(&input_bam, pairs);

    // Run codec command with cell-tag option
    let cmd = Codec::try_parse_from([
        "codec",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--min-reads",
        "1",
        "--min-duplex-length",
        "1",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse codec args");
    cmd.execute("fgumi codec").expect("Failed to run codec command");
    assert!(output_bam.exists(), "Output BAM not created");

    // Read output and verify cell barcode is preserved
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let mut consensus_count = 0;

    for result in reader.records() {
        let record = result.expect("Failed to read record");
        consensus_count += 1;

        // Verify cell barcode tag is preserved
        let cb_tag = SamTag::CB.to_noodles_tag();
        assert!(
            record.data().get(&cb_tag).is_some(),
            "Consensus should have CB (cell barcode) tag"
        );
    }

    assert!(consensus_count > 0, "Should have produced at least one consensus read");
}

/// Builds an FR CODEC pair with R1 and R2 covering an offset window so that
/// `start2 - start1` positions on each side fall outside the duplex overlap and
/// register as duplex disagreements (matching the unit-level `create_fr_pair`
/// in `crates/fgumi-consensus/src/codec_caller.rs`). Used by the
/// duplex-disagreement reject tests to drive `Codec::run` end-to-end through
/// the typed `CodecConsensusError` recovery path (issue #338).
#[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
fn create_offset_codec_pair(name: &str, umi: &str) -> (RawRecord, RawRecord) {
    // 40bp synthetic reference; R1 covers ref[0..30] (positions 1..30),
    // R2 covers ref[10..40] (positions 11..40). Overlap is ref[10..30] = 20bp;
    // 10 single-strand positions on each side count as duplex disagreements.
    const REF: &[u8; 40] = b"AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT";
    const READ_LEN: usize = 30;
    const QUAL: u8 = 35;

    let r1_seq = REF[..READ_LEN].to_vec();
    let r2_seq_fwd = REF[10..10 + READ_LEN].to_vec();
    // R2 is the reverse strand; BAM stores its sequence as reverse-complement
    // of the reference window so that single-strand consensus matches REF.
    let r2_seq = reverse_complement(&r2_seq_fwd);
    let cigar_op = (READ_LEN as u32) << 4; // 30M

    let r1_pos = 0_i32; // 0-based pos for ref position 1
    let r2_pos = 10_i32; // 0-based pos for ref position 11
    let mc_tag = format!("{READ_LEN}M");
    // Insert size from leftmost (R1) to rightmost (R2 + read length).
    let template_len = (10 + READ_LEN) as i32;

    let mut b1 = SamBuilder::new();
    b1.read_name(name.as_bytes())
        .sequence(&r1_seq)
        .qualities(&[QUAL; READ_LEN])
        .cigar_ops(&[cigar_op])
        .flags(raw_flags::PAIRED | raw_flags::FIRST_SEGMENT | raw_flags::MATE_REVERSE)
        .ref_id(0)
        .pos(r1_pos)
        .mapq(60)
        .mate_ref_id(0)
        .mate_pos(r2_pos)
        .template_length(template_len)
        .add_string_tag(SamTag::MI, umi.as_bytes())
        .add_string_tag(SamTag::MC, mc_tag.as_bytes());

    let mut b2 = SamBuilder::new();
    b2.read_name(name.as_bytes())
        .sequence(&r2_seq)
        .qualities(&[QUAL; READ_LEN])
        .cigar_ops(&[cigar_op])
        .flags(raw_flags::PAIRED | raw_flags::LAST_SEGMENT | raw_flags::REVERSE)
        .ref_id(0)
        .pos(r2_pos)
        .mapq(60)
        .mate_ref_id(0)
        .mate_pos(r1_pos)
        .template_length(-template_len)
        .add_string_tag(SamTag::MI, umi.as_bytes())
        .add_string_tag(SamTag::MC, mc_tag.as_bytes());

    (b1.build(), b2.build())
}

/// Verifies that the single-threaded `Codec::run` path recovers from a
/// duplex-disagreement molecule via the typed `CodecConsensusError` instead
/// of bailing — covers `src/lib/commands/codec.rs:383-397` (issue #338).
#[test]
fn test_codec_command_recovers_from_duplex_disagreement() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    // Single offset pair → 20 single-strand positions = 20 disagreements,
    // which trips the count threshold below.
    let (r1, r2) = create_offset_codec_pair("disagree_read", "UMI_DISAGREE");
    create_codec_test_bam(&input_bam, vec![(r1, r2)]);

    let cmd = Codec::try_parse_from([
        "codec",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--min-reads",
        "1",
        "--min-duplex-length",
        "1",
        // Strict: any disagreement rejects the molecule. Forces
        // `is_duplex_disagreement()` to return true so the loop
        // continues instead of returning the wrapped error.
        "--max-duplex-disagreements",
        "0",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse codec args");
    cmd.execute("fgumi codec")
        .expect("Codec command must succeed when only failure is a recoverable reject");

    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let consensus_count = reader.records().count();
    assert_eq!(consensus_count, 0, "Disagreeing molecule should produce no consensus");
}

/// Verifies that the parallel `--threads` path recovers from a
/// duplex-disagreement molecule via the typed `CodecConsensusError` — covers
/// what the single-threaded test above does not exercise.
///
/// Since the R3.7 codec-command cutover, `default = ["consensus"]` means
/// `--threads 2` here routes through the declarative chain builder
/// (`Codec::execute_chain` / `add_codec`), not the legacy
/// `Codec::execute_threads_mode`. This integration binary is built with the
/// workspace default features, so that chain path is what this test exercises —
/// a de-facto chain-path disagreement-recovery regression test.
#[test]
fn test_codec_command_recovers_from_duplex_disagreement_threaded() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let rejects_bam = temp_dir.path().join("rejects.bam");

    let (r1, r2) = create_offset_codec_pair("disagree_read_mt", "UMI_DISAGREE_MT");
    // Capture each read's *complete* raw BAM bytes (every field and tag), keyed by
    // flags, before the records are moved into the input BAM. The `--rejects`
    // contract is byte-for-byte preservation of the original records, so asserting
    // full-record identity — not just name + sequence — is what catches a pipeline
    // that silently drops or rewrites a field (quals, CIGAR, pos, mate info,
    // template length, MI/MC tags) on the reject path. R1 and R2 share a read name,
    // so flags is the distinguishing key.
    let expected_by_flags: std::collections::HashMap<u16, Vec<u8>> =
        [&r1, &r2].into_iter().map(|read| (read.flags(), read.as_ref().to_vec())).collect();
    create_codec_test_bam(&input_bam, vec![(r1, r2)]);

    let cmd = Codec::try_parse_from([
        "codec",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--rejects",
        rejects_bam.to_str().unwrap(),
        "--threads",
        "2",
        "--min-reads",
        "1",
        "--min-duplex-length",
        "1",
        "--max-duplex-disagreements",
        "0",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse codec args");
    cmd.execute("fgumi codec")
        .expect("Threaded codec command must succeed when only failure is a recoverable reject");

    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let consensus_count = reader.records().count();
    assert_eq!(consensus_count, 0, "Disagreeing molecule should produce no consensus");

    // Exercising --rejects forces the parallel-mode reject-collection path
    // through the typed-disagreement arm (`is_duplex_disagreement()`) — on a
    // consensus build that's the chain's `build_codec_consensus_step_with_rejects`
    // step (`add_codec`); on a non-consensus build it's the legacy
    // `execute_threads_mode`'s `flush_byte_records` call. CODEC3-08:
    // `consensus_reads_typed` preserves the disagreeing molecule's raw records
    // for the --rejects output before returning the recoverable error (fgbio
    // routes these to its rejectsWriter), so the two input reads land in the
    // rejects BAM on either path.
    assert!(rejects_bam.exists(), "Rejects BAM file should be created");

    // Read the rejects back as raw BAM records so we can compare the *entire*
    // serialized record (all fields + tags), not just the handful noodles decodes
    // ergonomically. Each reject must be byte-for-byte identical to one generated
    // read (keyed by flags), with no read matched twice.
    let (mut rejects_reader, _rejects_header) =
        fgumi_bam_io::create_raw_bam_reader(&rejects_bam, 1).expect("open rejects BAM");
    let mut matched_flags: std::collections::HashSet<u16> = std::collections::HashSet::new();
    let mut reject_count = 0;
    let mut record = RawRecord::new();
    while rejects_reader.read_record(&mut record).expect("read reject record") != 0 {
        reject_count += 1;

        let flags = record.flags();
        let expected_bytes = expected_by_flags
            .get(&flags)
            .unwrap_or_else(|| panic!("reject record has unexpected flags {flags}"));
        assert_eq!(
            record.as_ref(),
            expected_bytes.as_slice(),
            "reject record (flags={flags}) must be byte-for-byte identical to the generated \
             read — every field and tag preserved on the --rejects path"
        );
        assert!(
            matched_flags.insert(flags),
            "reject record with flags {flags} appears more than once"
        );
    }

    assert_eq!(
        reject_count, 2,
        "The disagreeing molecule's R1 and R2 should be written to the rejects BAM"
    );
    assert_eq!(
        matched_flags.len(),
        expected_by_flags.len(),
        "both generated reads (R1 and R2) should be present in the rejects BAM"
    );
}

/// Builds a CODEC pair whose R2 is unmapped, matching the geometry of
/// [`create_codec_read_pair`] so it lands in the same molecule.
#[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
fn create_half_mapped_codec_pair(
    name: &str,
    r1_seq: &[u8],
    r2_seq: &[u8],
    ref_start: usize,
    umi: &str,
) -> (RawRecord, RawRecord) {
    let r1_len = r1_seq.len();
    let r2_len = r2_seq.len();
    let r1_cigar_op = u32::try_from(r1_len).expect("r1_len fits u32") << 4;
    let pos = i32::try_from(ref_start).expect("ref_start fits i32") - 1;

    let mut b1 = SamBuilder::new();
    b1.read_name(name.as_bytes())
        .sequence(r1_seq)
        .qualities(&vec![30u8; r1_len])
        .cigar_ops(&[r1_cigar_op])
        .flags(raw_flags::PAIRED | raw_flags::FIRST_SEGMENT | raw_flags::MATE_REVERSE)
        .ref_id(0)
        .pos(pos)
        .mapq(60)
        .mate_ref_id(0)
        .mate_pos(pos)
        .template_length(r1_len as i32)
        .add_string_tag(SamTag::MI, umi.as_bytes());

    let mut b2 = SamBuilder::new();
    b2.read_name(name.as_bytes())
        .sequence(r2_seq)
        .qualities(&vec![30u8; r2_len])
        .flags(
            raw_flags::PAIRED | raw_flags::LAST_SEGMENT | raw_flags::REVERSE | raw_flags::UNMAPPED,
        )
        .ref_id(0)
        .pos(pos)
        .mapq(0)
        .mate_ref_id(0)
        .mate_pos(pos)
        .template_length(-(r2_len as i32))
        .add_string_tag(SamTag::MI, umi.as_bytes());

    (b1.build(), b2.build())
}

/// Neither end of a half-mapped pair may contribute to a CODEC consensus when mapped reads are
/// present in the same molecule. fgumi always prefers mapped reads.
#[test]
fn test_codec_ignores_half_mapped_pair_when_mapped_reads_present() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    let mut pairs = Vec::new();
    for i in 0..2 {
        pairs.push(create_codec_read_pair(
            &format!("mapped{i}"),
            b"ACGTACGT",
            b"ACGTACGT",
            &[30; 8],
            &[30; 8],
            100,
            "UMI001",
            None,
        ));
    }
    pairs.push(create_half_mapped_codec_pair(
        "halfmapped",
        b"ACGTACGT",
        b"TTTTTTTT",
        100,
        "UMI001",
    ));
    create_codec_test_bam(&input_bam, pairs);

    let cmd = Codec::try_parse_from([
        "codec",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--min-reads",
        "1",
        "--min-duplex-length",
        "1",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse codec args");
    cmd.execute("fgumi codec").expect("Failed to run codec command");

    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let records: Vec<bam::Record> =
        reader.records().map(|r| r.expect("failed to read record")).collect();
    assert_eq!(records.len(), 1, "expected a single CODEC consensus read");

    // Only the two fully mapped pairs contribute to each single-strand consensus.
    assert_eq!(int_tag(&records[0], SamTag::AD), 2, "R1 strand must use only the mapped reads");
    assert_eq!(int_tag(&records[0], SamTag::BD), 2, "R2 strand must use only the mapped reads");
}

// ─────────────────────────────────────────────────────────────────────────────
// R3.7 codec-command cutover: chain-vs-oracle parity tests
//
// The no-`--threads` single-threaded fast path (`Codec::execute`) is the
// in-process parity oracle; `--threads N` now routes onto the declarative
// chain builder (`Codec::execute_chain` -> `add_codec`). These tests prove
// byte/record identity between the two paths across the codec knobs, plus the
// Task 2 (sort-order guard) and Task 2A (validation mirroring) hardening.
// ─────────────────────────────────────────────────────────────────────────────

/// Parses and runs `codec` with the given input/output/extra args, panicking
/// on failure. Shared by the parity tests below.
fn codec_run(input: &Path, output: &Path, extra: &[&str]) {
    let mut args: Vec<&str> =
        vec!["codec", "--input", input.to_str().unwrap(), "--output", output.to_str().unwrap()];
    args.extend_from_slice(extra);
    Codec::try_parse_from(args)
        .expect("failed to parse codec args")
        .execute("fgumi codec")
        .expect("codec run failed");
}

/// Reads only the records of a codec output BAM, delegating to the shared
/// [`read_bam_output`] helper (full-header comparisons use that directly).
fn read_consensus_records(path: &Path) -> Vec<noodles::sam::alignment::RecordBuf> {
    read_bam_output(path).1
}

/// Builds `count` distinct-MI, full-overlap CODEC read pairs (one full duplex
/// molecule per pair) spread across increasing reference positions, using
/// [`create_codec_read_pair`]. Each pair is written to the BAM back-to-back,
/// which keeps `create_codec_test_bam`'s template-coordinate-sorted header
/// truthful (MI groups are contiguous in input order).
///
/// A large `count` pushes the chain's `GroupByMi` batching machinery to close
/// an internal batch mid-molecule-set, exercising cross-input-batch MI
/// bookkeeping that a handful of paired templates (or any single-record-per-
/// name fixture) never reaches.
fn create_paired_codec_templates(count: usize) -> Vec<(RawRecord, RawRecord)> {
    (0..count)
        .map(|i| {
            let name = format!("read{i:04}");
            let umi = format!("UMI{i:04}");
            let pos = 100 + i * 20;
            create_codec_read_pair(
                &name,
                b"ACGTACGT",
                b"ACGTACGT",
                &[30; 8],
                &[30; 8],
                pos,
                &umi,
                None,
            )
        })
        .collect()
}

/// Task 3 item 1: plain chain (`--threads 4`) vs oracle, full-header + record
/// parity.
#[test]
fn test_codec_chain_matches_single_threaded() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    create_codec_test_bam(&input_bam, create_paired_codec_templates(10));

    let oracle_out = temp_dir.path().join("oracle.bam");
    codec_run(&input_bam, &oracle_out, &["--min-reads", "1"]);

    let chain_out = temp_dir.path().join("chain.bam");
    codec_run(&input_bam, &chain_out, &["--min-reads", "1", "--threads", "4"]);

    let (oracle_header, expected) = read_bam_output(&oracle_out);
    let (chain_header, actual) = read_bam_output(&chain_out);
    assert!(!expected.is_empty(), "oracle output must be non-empty (guard against a vacuous pass)");
    assert_eq!(
        actual, expected,
        "chain --threads 4 output must match the non-chain path record-for-record"
    );
    assert_eq!(
        chain_header, oracle_header,
        "chain and non-chain output headers must match (with @PG CL normalized)"
    );
}

/// Knob-parity: set several non-default content-affecting codec knobs and assert
/// the chain (`--threads 4`) output matches the single-threaded oracle
/// record-for-record. Guards against `add_codec`'s hand-built
/// `CodecConsensusOptions` dropping or mis-wiring a knob that the oracle's
/// `to_codec_options()` honors — the default-only parity tests would pass green
/// on such a drop, since both paths would then agree on the (identical) default.
#[test]
fn test_codec_chain_matches_single_threaded_with_nondefault_knobs() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    create_codec_test_bam(&input_bam, create_paired_codec_templates(10));

    // Non-default values for every content-affecting knob the chain builds by
    // hand in `add_codec` (quality overrides, outer-base trimming, duplex-length
    // gate). Both paths get the SAME flags, so any divergence means the chain
    // dropped or mis-wired one of them.
    let knobs: &[&str] = &[
        "--min-reads",
        "1",
        "--min-duplex-length",
        "2",
        "--outer-bases-length",
        "3",
        "--outer-bases-qual",
        "10",
        "--single-strand-qual",
        "20",
    ];
    let mut chain_args: Vec<&str> = knobs.to_vec();
    chain_args.extend_from_slice(&["--threads", "4"]);

    let oracle_out = temp_dir.path().join("oracle.bam");
    codec_run(&input_bam, &oracle_out, knobs);
    let chain_out = temp_dir.path().join("chain.bam");
    codec_run(&input_bam, &chain_out, &chain_args);

    let (oracle_header, expected) = read_bam_output(&oracle_out);
    let (chain_header, actual) = read_bam_output(&chain_out);
    assert!(!expected.is_empty(), "oracle output must be non-empty (guard against a vacuous pass)");
    assert_eq!(
        actual, expected,
        "chain --threads 4 output with non-default knobs must match the oracle record-for-record"
    );
    assert_eq!(
        chain_header, oracle_header,
        "chain and non-chain output headers must match (with @PG CL normalized)"
    );
}

/// Task 3 item 4: multi-batch parity — many distinct-MI paired templates so
/// cross-input-batch MI carry-over inside the chain's `GroupByMi` is actually
/// exercised (a handful of templates all land in one batch and would never
/// reach this code path).
#[test]
fn test_codec_chain_multi_batch_matches_single_threaded() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    create_codec_test_bam(&input_bam, create_paired_codec_templates(200));

    let oracle_out = temp_dir.path().join("oracle.bam");
    codec_run(&input_bam, &oracle_out, &["--min-reads", "1"]);

    let chain_out = temp_dir.path().join("chain.bam");
    codec_run(&input_bam, &chain_out, &["--min-reads", "1", "--threads", "4"]);

    let expected = read_consensus_records(&oracle_out);
    let actual = read_consensus_records(&chain_out);
    assert!(!expected.is_empty(), "oracle output must be non-empty (guard against a vacuous pass)");
    assert_eq!(
        actual.len(),
        200,
        "every one of the 200 distinct-MI molecules must produce a consensus read"
    );
    assert_eq!(
        actual, expected,
        "chain multi-batch output must match the non-chain path record-for-record"
    );
}

/// Task 3 item 2: `--rejects` parity — record parity **and** header parity.
/// The rejects BAM must carry the raw input header verbatim on both paths (the
/// PR #332 contract; Task 2B), so a chain path that leaked its own `@PG` into
/// rejects would fail this.
#[test]
fn test_codec_chain_rejects_matches_single_threaded() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");

    let mut pairs = create_paired_codec_templates(5);
    // A duplex-disagreement molecule that `--max-duplex-disagreements 0` rejects
    // on both paths.
    pairs.push(create_offset_codec_pair("disagree_read", "UMI_DISAGREE"));
    create_codec_test_bam(&input_bam, pairs);

    let base_args = &["--min-reads", "1", "--max-duplex-disagreements", "0"];

    let oracle_out = temp_dir.path().join("oracle.bam");
    let oracle_rejects = temp_dir.path().join("oracle.rejects.bam");
    let mut oracle_args = base_args.to_vec();
    oracle_args.extend_from_slice(&["--rejects", oracle_rejects.to_str().unwrap()]);
    codec_run(&input_bam, &oracle_out, &oracle_args);

    let chain_out = temp_dir.path().join("chain.bam");
    let chain_rejects = temp_dir.path().join("chain.rejects.bam");
    let mut chain_args = base_args.to_vec();
    chain_args.extend_from_slice(&["--rejects", chain_rejects.to_str().unwrap(), "--threads", "4"]);
    codec_run(&input_bam, &chain_out, &chain_args);

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

    // Rejects-header provenance (Task 2B): rejects are raw-input records, so
    // both paths must write them under the raw input header verbatim.
    // `read_bam_output` normalizes the `@PG` CL, so the input header (which
    // carries no fgumi `@PG`) must equal both rejects headers.
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

/// Task 3 item 3: `--stats` output from the chain path is byte-identical to
/// the non-chain path (byte-compare, not field-by-field).
#[test]
fn test_codec_chain_stats_matches_single_threaded() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    create_codec_test_bam(&input_bam, create_paired_codec_templates(15));

    let oracle_out = temp_dir.path().join("oracle.bam");
    let oracle_stats = temp_dir.path().join("oracle.stats.txt");
    codec_run(
        &input_bam,
        &oracle_out,
        &["--min-reads", "1", "--stats", oracle_stats.to_str().unwrap()],
    );

    let chain_out = temp_dir.path().join("chain.bam");
    let chain_stats = temp_dir.path().join("chain.stats.txt");
    codec_run(
        &input_bam,
        &chain_out,
        &["--min-reads", "1", "--stats", chain_stats.to_str().unwrap(), "--threads", "4"],
    );

    let expected = fs::read(&oracle_stats).expect("read oracle stats");
    let actual = fs::read(&chain_stats).expect("read chain stats");
    assert!(!expected.is_empty(), "oracle stats must be non-empty (guard against a vacuous pass)");
    assert_eq!(
        actual, expected,
        "chain --stats output must be byte-identical to the non-chain path"
    );
}

/// Task 3 item 5 (proves Task 2): a coordinate-sorted input is REJECTED on
/// both the non-chain and chain (`--threads 4`) paths with the same
/// "not sorted correctly for consensus calling" error. `add_codec` has no
/// out-of-order/duplicate detection of its own (`GroupByMi` doesn't either),
/// so without the `check_consensus_sort_order` guard a mis-sorted input would
/// be silently mis-grouped instead of rejected.
#[test]
fn test_codec_chain_rejects_coordinate_sorted_input() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");

    let header = create_coordinate_sorted_header("chr1", 10000);
    let pair = create_codec_read_pair(
        "read1",
        b"ACGTACGT",
        b"ACGTACGT",
        &[30; 8],
        &[30; 8],
        100,
        "UMI1",
        None,
    );
    let mut writer =
        create_raw_bam_writer(&input_bam, &header, 1, 6).expect("failed to create raw BAM writer");
    writer.write_raw_record(pair.0.as_ref()).expect("failed to write R1");
    writer.write_raw_record(pair.1.as_ref()).expect("failed to write R2");
    writer.finish().expect("failed to finish BAM");

    let output_bam = temp_dir.path().join("output.bam");

    let oracle_err = Codec::try_parse_from([
        "codec",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--min-reads",
        "1",
    ])
    .expect("failed to parse codec args")
    .execute("fgumi codec")
    .expect_err("codec must reject non-template-coordinate input");

    let chain_err = Codec::try_parse_from([
        "codec",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--min-reads",
        "1",
        "--threads",
        "4",
    ])
    .expect("failed to parse codec args")
    .execute("fgumi codec")
    .expect_err("codec chain path must reject non-template-coordinate input");

    assert!(
        oracle_err.to_string().contains("template-coordinate"),
        "oracle error must suggest template-coordinate sorting: {oracle_err}"
    );
    assert_eq!(
        oracle_err.to_string(),
        chain_err.to_string(),
        "both paths must reject a coordinate-sorted input with the same error message"
    );
}

/// Task 3 item 6a: an out-of-range `--max-duplex-disagreement-rate` is
/// rejected identically on both `--threads`/no-`--threads` paths with the same
/// error. This passes with or without Task 2A's mirrored checks in
/// `add_codec`, because `Codec::execute()`'s top-level `validate()` rejects it
/// before either branch dispatches — it is a parity/regression guard, not an
/// isolation test for `add_codec`'s new checks (see the builder-level test
/// below for that).
#[rstest]
#[case::single_threaded(&[])]
#[case::threaded(&["--threads", "4"])]
fn test_codec_rejects_out_of_range_disagreement_rate_end_to_end(#[case] extra: &[&str]) {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    create_codec_test_bam(&input_bam, create_paired_codec_templates(2));
    let output_bam = temp_dir.path().join("output.bam");

    let mut args = vec![
        "codec",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--min-reads",
        "1",
        "--max-duplex-disagreement-rate",
        "1.5",
    ];
    args.extend_from_slice(extra);

    let err = Codec::try_parse_from(args)
        .expect("failed to parse codec args")
        .execute("fgumi codec")
        .expect_err("codec must reject an out-of-range --max-duplex-disagreement-rate");
    assert!(
        err.to_string().contains("max-duplex-disagreement-rate must be between 0.0 and 1.0"),
        "unexpected error: {err:#}"
    );
}

/// Task 3 item 6b (the test that actually pins Task 2A): drives `add_codec`
/// directly via `ChainSpec::single_stage(Stage::Codec, ..)` with an
/// out-of-range `CodecOptions`, bypassing `Codec::execute`'s top-level
/// `validate()` entirely. This is the test that fails if Task 2A's mirrored
/// checks were absent from `add_codec` — the end-to-end test above cannot
/// distinguish "the builder validates" from "the CLI's `validate()` already
/// caught it first".
#[test]
fn test_add_codec_rejects_out_of_range_disagreement_rate_bypassing_cli_validate() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    create_codec_test_bam(&input_bam, create_paired_codec_templates(2));
    let output_bam = temp_dir.path().join("output.bam");

    let cmd = Codec::try_parse_from([
        "codec",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--min-reads",
        "1",
    ])
    .expect("failed to parse codec args");

    // Bypass `Codec::execute`'s top-level `validate()` entirely by going
    // straight through `to_codec_options()` into a hand-built `ChainSpec` --
    // never calling `Command::execute`/`Codec::validate` at all.
    let mut codec_opts = cmd.to_codec_options();
    codec_opts.max_duplex_disagreement_rate = 1.5;

    let ctx = SingleStageContext {
        io: &cmd.io,
        threading: &cmd.threading,
        compression: &cmd.compression,
        scheduler: &cmd.scheduler_opts,
        queue_memory: &cmd.queue_memory,
        command_line: "test",
    };
    let stage_opts = StageOptionsBag { codec: Some(codec_opts), ..Default::default() };
    let spec = ChainSpec::single_stage(Stage::Codec, stage_opts, &ctx);

    let err = build_for(spec).and_then(fgumi_lib::pipeline::chains::BuiltPipeline::run).expect_err(
        "add_codec must reject an out-of-range disagreement rate itself, even when \
             Codec::execute's top-level validate() is bypassed entirely",
    );
    assert!(
        err.to_string().contains("max-duplex-disagreement-rate must be between 0.0 and 1.0"),
        "unexpected error: {err:#}"
    );
}

/// Task 2A (exhaustive): pin every bound `CodecOptions::validate()` enforces.
/// `validate()` is the single source of truth that BOTH `Codec::validate()` (the
/// CLI pre-flight) and `add_codec` (the chain builder) delegate to, so a guard
/// dropped from it silently weakens both paths at once — a green end-to-end
/// parity test cannot catch that, because both paths would drop the guard
/// together. Driving `validate()` directly (a pure function; no BAM/I/O needed)
/// is the isolation test that would fail if any one check were removed.
#[rstest]
#[case::error_rate_pre_umi_zero(
    |o: &mut CodecOptions| o.error_rate_pre_umi = 0,
    "error-rate-pre-umi must be > 0"
)]
#[case::error_rate_post_umi_zero(
    |o: &mut CodecOptions| o.error_rate_post_umi = 0,
    "error-rate-post-umi must be > 0"
)]
#[case::min_reads_zero(
    |o: &mut CodecOptions| o.min_reads = 0,
    "min-reads must be >= 1"
)]
#[case::max_reads_below_min_reads(
    |o: &mut CodecOptions| { o.min_reads = 3; o.max_reads = Some(2); },
    "max-reads (2) must be >= min-reads (3)"
)]
#[case::min_duplex_length_zero(
    |o: &mut CodecOptions| o.min_duplex_length = 0,
    "min-duplex-length must be >= 1"
)]
#[case::disagreement_rate_nan(
    |o: &mut CodecOptions| o.max_duplex_disagreement_rate = f64::NAN,
    "max-duplex-disagreement-rate must be a finite number"
)]
#[case::disagreement_rate_above_one(
    |o: &mut CodecOptions| o.max_duplex_disagreement_rate = 1.5,
    "max-duplex-disagreement-rate must be between 0.0 and 1.0"
)]
#[case::disagreement_rate_negative(
    |o: &mut CodecOptions| o.max_duplex_disagreement_rate = -0.1,
    "max-duplex-disagreement-rate must be between 0.0 and 1.0"
)]
#[case::single_strand_qual_above_max_phred(
    |o: &mut CodecOptions| o.single_strand_qual = Some(200),
    "single-strand-qual (200) exceeds maximum Phred score"
)]
#[case::outer_bases_qual_above_max_phred(
    |o: &mut CodecOptions| o.outer_bases_qual = Some(200),
    "outer-bases-qual (200) exceeds maximum Phred score"
)]
fn test_codec_options_validate_rejects_invalid(
    #[case] mutate: fn(&mut CodecOptions),
    #[case] expected_msg: &str,
) {
    // A valid baseline projected from a normally-parsed command. `validate()` is
    // pure, so the referenced paths are never opened.
    let cmd = Codec::try_parse_from([
        "codec",
        "--input",
        "in.bam",
        "--output",
        "out.bam",
        "--min-reads",
        "1",
    ])
    .expect("failed to parse codec args");
    let mut opts = cmd.to_codec_options();
    opts.validate().expect("baseline CodecOptions must be valid (or the case proves nothing)");

    mutate(&mut opts);
    let err = opts.validate().expect_err("validate() must reject the mutated options");
    assert!(
        err.to_string().contains(expected_msg),
        "expected error containing {expected_msg:?}, got: {err:#}"
    );
}

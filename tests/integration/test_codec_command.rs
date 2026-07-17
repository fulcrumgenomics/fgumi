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
use fgumi_lib::commands::codec::Codec;
use fgumi_lib::commands::command::Command;
use fgumi_lib::sam::SamTag;
use fgumi_raw_bam::{RawRecord, SamBuilder, flags as raw_flags};
use noodles::bam;
use std::fs;
use std::path::PathBuf;
use tempfile::TempDir;

use crate::helpers::bam_generator::create_minimal_header;

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

/// Verifies that the parallel `Codec::execute_threads_mode` path recovers
/// from a duplex-disagreement molecule via the typed `CodecConsensusError`
/// — covers `src/lib/commands/codec.rs:617-630` (issue #338) which the
/// single-threaded test above does not exercise.
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

    // Exercising --rejects forces the parallel-mode `flush_byte_records`
    // call inside the typed-disagreement arm (`is_duplex_disagreement()`
    // branch in `execute_threads_mode`). CODEC3-08: `consensus_reads_typed`
    // now preserves the disagreeing molecule's raw records for the --rejects
    // output before returning the recoverable error (fgbio routes these to its
    // rejectsWriter), so the two input reads land in the rejects BAM.
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

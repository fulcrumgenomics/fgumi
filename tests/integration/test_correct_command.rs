//! End-to-end tests for the correct command.
//!
//! These tests invoke `CorrectUmis::execute()` in-process and validate:
//! 1. Basic UMI correction against a whitelist
//! 2. Metrics output
//! 3. Rejected reads output

use clap::Parser;
use fgumi_lib::commands::command::Command;
use fgumi_lib::commands::correct::{CorrectUmis, Target};
use fgumi_lib::sam::SamTag;
use fgumi_raw_bam::{RawRecord, SamBuilder, flags};
use noodles::bam;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value;
use rstest::rstest;
use std::collections::BTreeMap;
use std::fs;
use std::path::PathBuf;
use tempfile::TempDir;

use crate::helpers::bam_generator::{
    create_family_with_tag, create_minimal_header, create_umi_family, to_record_buf,
};
use crate::helpers::read_bam_output;

/// Write a BAM with UMI-tagged reads.
fn create_umi_bam(path: &PathBuf, families: Vec<Vec<RawRecord>>) {
    let header = create_minimal_header("chr1", 10000);
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");

    for family in families {
        for record in &family {
            writer
                .write_alignment_record(&header, &to_record_buf(record))
                .expect("Failed to write record");
        }
    }
    writer.try_finish().expect("Failed to finish BAM");
}

/// Write a UMI whitelist file.
fn create_whitelist(path: &PathBuf, umis: &[&str]) {
    fs::write(path, umis.join("\n")).expect("Failed to write whitelist");
}

/// Read the output BAM into an independent per-record oracle: a map from read
/// name to that record's values for `tags` (aligned to `tags`, `None` when a
/// tag is absent). Enables order-independent, per-record identity assertions
/// (`{QNAME -> (RX, BC, OX, ob)}`) instead of aggregate counts. Panics on a
/// duplicate read name so the oracle stays unambiguous.
fn records_by_name(path: &PathBuf, tags: &[SamTag]) -> BTreeMap<String, Vec<Option<String>>> {
    let noodles_tags: Vec<Tag> = tags.iter().map(|t| t.to_noodles_tag()).collect();
    let mut reader = bam::io::Reader::new(fs::File::open(path).expect("open bam"));
    let header = reader.read_header().expect("read header");
    let mut out: BTreeMap<String, Vec<Option<String>>> = BTreeMap::new();
    for result in reader.record_bufs(&header) {
        let record = result.expect("read record");
        let name = record
            .name()
            .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
            .expect("record missing read name");
        let values = noodles_tags
            .iter()
            .zip(tags)
            .map(|(noodles_tag, tag)| {
                record.data().get(noodles_tag).map(|value| match value {
                    Value::String(s) => s.to_string(),
                    other => panic!("{tag} tag is not a string: {other:?}"),
                })
            })
            .collect();
        assert!(out.insert(name.clone(), values).is_none(), "duplicate read name {name}");
    }
    out
}

/// Build the expected per-record oracle for `records_by_name(.., &[seq_tag,
/// original_tag])` over a fixture of `exact_{0..exact}` reads (correct as-is,
/// no original stored) and `corr_{0..corrected}` reads (one-mismatch, corrected
/// with the original stored). `store_original` reflects whether the run kept
/// the original value (`false` under `--dont-store-original`).
fn expected_target_records(
    exact: usize,
    corrected: usize,
    store_original: bool,
) -> BTreeMap<String, Vec<Option<String>>> {
    let seq = Some("ACGTACGT".to_string());
    let original = store_original.then(|| "ACGTACGA".to_string());
    let mut expected: BTreeMap<String, Vec<Option<String>>> = BTreeMap::new();
    for i in 0..exact {
        expected.insert(format!("exact_{i}"), vec![seq.clone(), None]);
    }
    for i in 0..corrected {
        expected.insert(format!("corr_{i}"), vec![seq.clone(), original.clone()]);
    }
    expected
}

/// Test basic UMI correction.
#[test]
fn test_correct_command_basic() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let whitelist = temp_dir.path().join("whitelist.txt");

    // Create reads: 5 with correct UMI "ACGTACGT" and 2 with 1bp error "ACGTACGA"
    let correct_reads = create_umi_family("ACGTACGT", 5, "correct", "AAAAGGGG", 30);
    let error_reads = create_umi_family("ACGTACGA", 2, "error", "AAAAGGGG", 30);
    create_umi_bam(&input_bam, vec![correct_reads, error_reads]);
    create_whitelist(&whitelist, &["ACGTACGT"]);

    let cmd = CorrectUmis::try_parse_from([
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
        "--compression-level",
        "1",
    ])
    .expect("failed to parse correct args");
    cmd.execute("fgumi correct").expect("Correct command failed");
    assert!(output_bam.exists(), "Output BAM not created");

    // Verify output has records
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let count = reader.records().count();
    assert_eq!(count, 7, "Should have all 7 reads in output");
}

/// Test correct command with metrics output.
#[test]
fn test_correct_command_with_metrics() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let whitelist = temp_dir.path().join("whitelist.txt");
    let metrics = temp_dir.path().join("metrics.tsv");

    let reads = create_umi_family("ACGTACGT", 3, "read", "AAAAGGGG", 30);
    create_umi_bam(&input_bam, vec![reads]);
    create_whitelist(&whitelist, &["ACGTACGT"]);

    let cmd = CorrectUmis::try_parse_from([
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
        "--metrics",
        metrics.to_str().unwrap(),
        "--compression-level",
        "1",
    ])
    .expect("failed to parse correct args");
    cmd.execute("fgumi correct").expect("Correct command with metrics failed");
    assert!(metrics.exists(), "Metrics file not created");

    let content = fs::read_to_string(&metrics).unwrap();
    assert!(!content.is_empty(), "Metrics file should not be empty");
}

/// Test correct command with rejects output.
#[test]
fn test_correct_command_with_rejects() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let rejects_bam = temp_dir.path().join("rejects.bam");
    let whitelist = temp_dir.path().join("whitelist.txt");

    // UMI "TTTTTTTT" has edit distance 8 from "ACGTACGT" — won't correct
    let uncorrectable = create_umi_family("TTTTTTTT", 2, "far", "AAAAGGGG", 30);
    let correctable = create_umi_family("ACGTACGT", 3, "exact", "AAAAGGGG", 30);
    create_umi_bam(&input_bam, vec![correctable, uncorrectable]);
    create_whitelist(&whitelist, &["ACGTACGT"]);

    let cmd = CorrectUmis::try_parse_from([
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
        "--rejects",
        rejects_bam.to_str().unwrap(),
        "--compression-level",
        "1",
    ])
    .expect("failed to parse correct args");
    cmd.execute("fgumi correct").expect("Correct command with rejects failed");
    assert!(rejects_bam.exists(), "Rejects BAM not created");
}

/// Exercises the multi-threaded rejects-streaming path end-to-end and asserts:
/// 1. The rejects BAM is a valid BGZF stream with the terminating EOF block.
/// 2. The `@HD` sort fields (`SO`/`GO`/`SS`) match the input BAM, because
///    rejects flow through the unified pipeline's secondary output in
///    batch/input order (a subset of an SO-X stream is still SO-X).
/// 3. Every uncorrectable input record appears exactly once in the rejects
///    BAM — the writer does not drop records under worker contention and does
///    not emit duplicates.
#[test]
fn test_correct_command_rejects_streaming_threaded_integrity() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let rejects_bam = temp_dir.path().join("rejects.bam");
    let whitelist = temp_dir.path().join("whitelist.txt");

    // Each record has a unique QNAME, so each record is its own template.
    // `CorrectUmis` uses a 1000-template batch, so size the uncorrectable
    // families so that `far_a + far_b + far_c` exceeds that boundary and
    // multiple batches flow through the 4-thread pool concurrently, letting
    // more than one worker race to flush rejects.
    let far_family_size: u32 = 400;
    let corr_small = create_umi_family("ACGTACGT", 3, "c1_exact", "AAAAGGGG", 30);
    let corr_big = create_umi_family("ACGTACGT", 30, "c2_exact", "AAAAGGGG", 30);
    let far_a = create_umi_family("TTTTTTTT", far_family_size as usize, "far_a", "AAAAGGGG", 30);
    let far_b = create_umi_family("GGGGGGGG", far_family_size as usize, "far_b", "AAAAGGGG", 30);
    let far_c = create_umi_family("CCCCCCCC", far_family_size as usize, "far_c", "AAAAGGGG", 30);

    // Expected reject names mirror `create_umi_family`'s "{base_name}_{i}"
    // convention for the three uncorrectable families. The vector encodes the
    // batch-input order the pipeline now guarantees for rejects (far_a, then
    // far_b, then far_c — each family's records emitted in their original
    // 0..far_family_size order).
    let expected_order: Vec<String> = ["far_a", "far_b", "far_c"]
        .into_iter()
        .flat_map(|base| (0..far_family_size).map(move |i| format!("{base}_{i}")))
        .collect();
    let expected_names: std::collections::HashSet<String> =
        expected_order.iter().cloned().collect();
    let expected_rejects = expected_names.len();

    create_umi_bam(&input_bam, vec![corr_small, corr_big, far_a, far_b, far_c]);
    create_whitelist(&whitelist, &["ACGTACGT"]);

    let cmd = CorrectUmis::try_parse_from([
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
        "--rejects",
        rejects_bam.to_str().unwrap(),
        "--threads",
        "4",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse correct args");
    cmd.execute("fgumi correct").expect("Correct command with threaded rejects failed");
    assert!(rejects_bam.exists(), "Rejects BAM not created");

    crate::helpers::assertions::assert_has_bgzf_eof(&rejects_bam);
    crate::helpers::assertions::assert_rejects_header_matches_input(&rejects_bam, &input_bam);

    let mut reader = bam::io::Reader::new(fs::File::open(&rejects_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let mut observed_order: Vec<String> = Vec::with_capacity(expected_rejects);
    for result in reader.records() {
        let record = result.expect("Failed to read reject record");
        let name = record.name().expect("reject record missing read name").to_string();
        observed_order.push(name);
    }

    assert_eq!(
        observed_order.len(),
        expected_rejects,
        "rejects BAM should contain one record per uncorrectable input read",
    );
    // Strict equality on the full sequence catches both drops/dupes (covered
    // by the count check) and a regression to mutex-acquisition (or any other
    // non-input) ordering.
    assert_eq!(observed_order, expected_order, "rejects should be emitted in batch-input order",);
    let observed_names: std::collections::HashSet<String> =
        observed_order.iter().cloned().collect();
    assert_eq!(observed_names, expected_names, "unexpected reject-name set");
}

/// Corrected value lands in the target's sequence tag (`RX` for UMI, `BC` for
/// barcode) and the pre-correction original lands in the target's original
/// tag (`OX` for UMI, `ob` for barcode); exact-match reads carry no original.
///
/// Asserts per-record identity (`{QNAME -> (seq_tag, original_tag)}`) against
/// an independently constructed expected map, so a broken tag mapping — or one
/// that drops/duplicates records — fails rather than passing an aggregate
/// count. `seq_tag`/`original_tag` come from the case table (literal tags),
/// not from `Target::sequence_tag()`, so the oracle does not derive its
/// expectations from the code under test. Both the single-thread (`None`) and
/// multi-threaded (`Some(2)`) pipeline paths are covered per target.
#[rstest]
#[case::umi_serial(Target::Umi, SamTag::RX, SamTag::OX, None)]
#[case::umi_threaded(Target::Umi, SamTag::RX, SamTag::OX, Some(2))]
#[case::barcode_serial(Target::Barcode, SamTag::BC, SamTag::OB, None)]
#[case::barcode_threaded(Target::Barcode, SamTag::BC, SamTag::OB, Some(2))]
fn correct_writes_corrected_value_and_original_by_target(
    #[case] target: Target,
    #[case] seq_tag: SamTag,
    #[case] original_tag: SamTag,
    #[case] threads: Option<usize>,
) {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let whitelist = temp_dir.path().join("whitelist.txt");

    // 3 exact-match + 2 one-mismatch (correctable) reads on the target tag.
    let exact = create_family_with_tag(seq_tag, "ACGTACGT", 3, "exact", "AAAAGGGG", 30);
    let corrected = create_family_with_tag(seq_tag, "ACGTACGA", 2, "corr", "AAAAGGGG", 30);
    create_umi_bam(&input_bam, vec![exact, corrected]);
    create_whitelist(&whitelist, &["ACGTACGT"]);

    let target_str = match target {
        Target::Umi => "umi",
        Target::Barcode => "barcode",
    };
    let mut args = vec![
        "correct".to_string(),
        "-i".to_string(),
        input_bam.to_str().unwrap().to_string(),
        "-o".to_string(),
        output_bam.to_str().unwrap().to_string(),
        "-U".to_string(),
        whitelist.to_str().unwrap().to_string(),
        "--max-mismatches".to_string(),
        "1".to_string(),
        "--min-distance".to_string(),
        "1".to_string(),
        "--target".to_string(),
        target_str.to_string(),
    ];
    if let Some(n) = threads {
        args.push("--threads".to_string());
        args.push(n.to_string());
    }
    let cmd = CorrectUmis::try_parse_from(args).expect("parse");
    cmd.execute("test").expect("correct runs");

    // Independent per-record oracle: each named record must carry exactly the
    // expected sequence-tag and original-tag values.
    let actual = records_by_name(&output_bam, &[seq_tag, original_tag]);
    let expected = expected_target_records(3, 2, true);
    assert_eq!(actual, expected);
}

/// `--dont-store-original` suppresses the original-tag write for both targets,
/// on both pipeline paths. The corrected reads are still kept and carry the
/// corrected value in the sequence tag — asserted via the per-record oracle so
/// the suppression check cannot pass vacuously on an empty (all-rejected)
/// output. `seq_tag`/`original_tag` are literal case-table tags, decoupled from
/// `Target::sequence_tag()`.
#[rstest]
#[case::umi_serial(Target::Umi, SamTag::RX, SamTag::OX, None)]
#[case::umi_threaded(Target::Umi, SamTag::RX, SamTag::OX, Some(2))]
#[case::barcode_serial(Target::Barcode, SamTag::BC, SamTag::OB, None)]
#[case::barcode_threaded(Target::Barcode, SamTag::BC, SamTag::OB, Some(2))]
fn correct_dont_store_original_suppresses_original_by_target(
    #[case] target: Target,
    #[case] seq_tag: SamTag,
    #[case] original_tag: SamTag,
    #[case] threads: Option<usize>,
) {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let whitelist = temp_dir.path().join("whitelist.txt");

    // 3 exact-match + 2 one-mismatch (correctable) reads on the target tag, so
    // the oracle can confirm kept records exist (guarding a vacuous pass).
    let exact = create_family_with_tag(seq_tag, "ACGTACGT", 3, "exact", "AAAAGGGG", 30);
    let corrected = create_family_with_tag(seq_tag, "ACGTACGA", 2, "corr", "AAAAGGGG", 30);
    create_umi_bam(&input_bam, vec![exact, corrected]);
    create_whitelist(&whitelist, &["ACGTACGT"]);

    let target_str = match target {
        Target::Umi => "umi",
        Target::Barcode => "barcode",
    };
    let mut args = vec![
        "correct".to_string(),
        "-i".to_string(),
        input_bam.to_str().unwrap().to_string(),
        "-o".to_string(),
        output_bam.to_str().unwrap().to_string(),
        "-U".to_string(),
        whitelist.to_str().unwrap().to_string(),
        "--max-mismatches".to_string(),
        "1".to_string(),
        "--min-distance".to_string(),
        "1".to_string(),
        "--target".to_string(),
        target_str.to_string(),
        "--dont-store-original".to_string(),
    ];
    if let Some(n) = threads {
        args.push("--threads".to_string());
        args.push(n.to_string());
    }
    let cmd = CorrectUmis::try_parse_from(args).expect("parse");
    cmd.execute("test").expect("correct runs");

    // Corrected reads are still present and corrected, but NO original was
    // stored (store_original = false) — even for the corrected reads.
    let actual = records_by_name(&output_bam, &[seq_tag, original_tag]);
    let expected = expected_target_records(3, 2, false);
    assert_eq!(actual, expected);
}

/// Records carrying BOTH an `RX` (UMI) tag and a `BC` (barcode) tag, when run
/// through `--target barcode`, should have only `BC`/`ob` touched: `BC` is
/// corrected as expected, while `RX` is left byte-for-byte as the input wrote
/// it and no `OX` tag is written. This guards against a correction pass
/// bleeding across the two tag pairs when both are present on a record.
#[test]
fn correct_barcode_leaves_umi_tags_untouched() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let whitelist = temp_dir.path().join("whitelist.txt");

    // Build families via the RX helper, then add a second BC tag to every
    // record so each carries both tags: an untouched, unmatched-to-any-fixed
    // UMI on RX and a barcode on BC that we want corrected.
    let build_dual_tag_family = |bc: &str, depth: usize, base_name: &str| -> Vec<RawRecord> {
        create_family_with_tag(SamTag::RX, "TTTTTTTT", depth, base_name, "AAAAGGGG", 30)
            .into_iter()
            .map(|mut record| {
                fgumi_raw_bam::update_string_tag(record.as_mut_vec(), SamTag::BC, bc.as_bytes());
                record
            })
            .collect()
    };

    let exact = build_dual_tag_family("ACGTACGT", 3, "exact");
    let corrected = build_dual_tag_family("ACGTACGA", 2, "corr");
    create_umi_bam(&input_bam, vec![exact, corrected]);
    create_whitelist(&whitelist, &["ACGTACGT"]);

    let cmd = CorrectUmis::try_parse_from([
        "correct",
        "-i",
        input_bam.to_str().unwrap(),
        "-o",
        output_bam.to_str().unwrap(),
        "-U",
        whitelist.to_str().unwrap(),
        "--max-mismatches",
        "1",
        "--min-distance",
        "1",
        "--target",
        "barcode",
    ])
    .expect("parse");
    cmd.execute("test").expect("correct runs");

    // Per-record identity over both tag pairs (`RX`, `BC`, `OX`, `ob`): `BC` is
    // corrected to the fixed barcode with the original stashed in `ob` for the
    // corrected reads only, while the UMI pair is left exactly as the input
    // wrote it — `RX` still the uncorrected "TTTTTTTT" and no `OX` on any
    // record. Asserting full per-record identity (not aggregates) catches any
    // cross-contamination between the two tag pairs.
    let actual = records_by_name(&output_bam, &[SamTag::RX, SamTag::BC, SamTag::OX, SamTag::OB]);
    let rx = || Some("TTTTTTTT".to_string());
    let bc = || Some("ACGTACGT".to_string());
    let mut expected: BTreeMap<String, Vec<Option<String>>> = BTreeMap::new();
    for i in 0..3 {
        // exact-match barcode reads: BC unchanged in value, no `ob`.
        expected.insert(format!("exact_{i}"), vec![rx(), bc(), None, None]);
    }
    for i in 0..2 {
        // corrected barcode reads: BC corrected, original barcode in `ob`.
        expected.insert(format!("corr_{i}"), vec![rx(), bc(), None, Some("ACGTACGA".to_string())]);
    }
    assert_eq!(actual, expected);
}

// ============================================================================
// --check-crc / --no-check-crc on the single-threaded fast path (#800)
// ============================================================================

/// Flip a byte in the last BGZF block's CRC32 footer, so decoding that block
/// fails only when CRC verification is on. Requires the file to span at least
/// two blocks so the corrupted block is not the header's block (the header
/// parse always verifies).
fn corrupt_last_block_crc(path: &std::path::Path) {
    let mut bytes = fs::read(path).expect("read bam for corruption");
    let mut cursor: &[u8] = &bytes;
    let blocks = fgumi_lib::bgzf_reader::read_raw_blocks(&mut cursor, 100_000)
        .expect("read bgzf blocks from test bam");
    assert!(
        blocks.len() >= 2,
        "test input must span >= 2 BGZF blocks so the corrupted block isn't the header's; \
         got {} -- generate more records",
        blocks.len()
    );
    let offset: usize =
        blocks[..blocks.len() - 1].iter().map(fgumi_lib::bgzf_reader::RawBgzfBlock::len).sum();
    let last = blocks.last().expect("checked len >= 2 above");
    let crc_off = offset + last.len() - fgumi_lib::bgzf_reader::BGZF_FOOTER_SIZE;
    bytes[crc_off] ^= 0x01;
    fs::write(path, bytes).expect("write corrupted bam");
}

/// Write a single UMI family large enough to span more than one BGZF block, so
/// the corrupted last block is a record block rather than the header's.
fn create_multiblock_umi_bam(path: &PathBuf, umi: &str) {
    create_umi_bam(path, vec![create_umi_family(umi, 3000, "read", "ACGTACGT", 30)]);
}

/// Build correct's argv for the single-threaded fast path (no `--threads`),
/// appending any extra flags (e.g. `--no-check-crc`).
fn correct_args<'a>(input: &'a str, output: &'a str, extra: &[&'a str]) -> Vec<&'a str> {
    let mut args = vec![
        "correct",
        "--input",
        input,
        "--output",
        output,
        "--umis",
        "ACGTACGT",
        "--max-mismatches",
        "1",
        "--min-distance",
        "1",
        "--compression-level",
        "1",
    ];
    args.extend_from_slice(extra);
    args
}

/// `--no-check-crc` must let correct's single-threaded reader accept a corrupted
/// BGZF CRC32: it decodes through fgumi-bgzf, honoring the flag (#800). Against
/// the noodles reader this path used before, this could not pass.
#[test]
fn test_correct_no_check_crc_accepts_corrupted_crc() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    create_multiblock_umi_bam(&input_bam, "ACGTACGT");
    corrupt_last_block_crc(&input_bam);

    let cmd = CorrectUmis::try_parse_from(correct_args(
        input_bam.to_str().unwrap(),
        output_bam.to_str().unwrap(),
        &["--no-check-crc"],
    ))
    .expect("failed to parse correct args");
    cmd.execute("fgumi correct")
        .expect("--no-check-crc must accept a corrupted BGZF CRC32 and complete");
    assert!(output_bam.exists(), "Output BAM not created");
}

/// Default (verify-on for file input) must reject the same corrupted CRC32.
#[test]
fn test_correct_rejects_corrupted_crc_by_default() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    create_multiblock_umi_bam(&input_bam, "ACGTACGT");
    corrupt_last_block_crc(&input_bam);

    let cmd = CorrectUmis::try_parse_from(correct_args(
        input_bam.to_str().unwrap(),
        output_bam.to_str().unwrap(),
        &[],
    ))
    .expect("failed to parse correct args");
    let err = cmd
        .execute("fgumi correct")
        .expect_err("default (verify-on for file input) must reject a corrupted BGZF CRC32");
    let message = format!("{err:#}");
    assert!(message.to_uppercase().contains("CRC32"), "error should mention CRC32: {message}");
}

/// `--check-crc` must also reject the corrupted CRC32 (forces verification on).
#[test]
fn test_correct_check_crc_rejects_corrupted_crc() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    create_multiblock_umi_bam(&input_bam, "ACGTACGT");
    corrupt_last_block_crc(&input_bam);

    let cmd = CorrectUmis::try_parse_from(correct_args(
        input_bam.to_str().unwrap(),
        output_bam.to_str().unwrap(),
        &["--check-crc"],
    ))
    .expect("failed to parse correct args");
    let err =
        cmd.execute("fgumi correct").expect_err("--check-crc must reject a corrupted BGZF CRC32");
    let message = format!("{err:#}");
    assert!(message.to_uppercase().contains("CRC32"), "error should mention CRC32: {message}");
}

// ============================================================================
// Chain (`--threads N`) vs single-threaded oracle parity (R3.2 cutover)
// ============================================================================

/// Read the complete `RecordBuf` for every output record, in output order, via
/// the shared `read_bam_output` helper. Positional `Vec` equality across two
/// runs catches re-ordering, drops, duplicates, and any per-record field
/// difference. (The full-header half of `read_bam_output` is used directly in
/// the main parity test for whole-header comparison.)
fn read_correct_records(path: &std::path::Path) -> Vec<noodles::sam::alignment::RecordBuf> {
    read_bam_output(path).1
}

/// Runs `correct --umi-files <whitelist> --max-mismatches 1
/// --compression-level 1` plus `extra` (e.g. `--threads`, `--rejects`,
/// `--metrics`, `--min-distance`), panicking with the command's error on
/// failure. `--min-distance` is left to `extra` (rather than baked in here)
/// so a case that needs a non-default value doesn't pass the flag twice --
/// clap rejects a `--min-distance` argument repeated across the base args
/// and `extra`.
fn run_correct(
    input: &std::path::Path,
    output: &std::path::Path,
    whitelist: &std::path::Path,
    extra: &[&str],
    tag: &str,
) {
    let mut args = vec![
        "correct",
        "--input",
        input.to_str().unwrap(),
        "--output",
        output.to_str().unwrap(),
        "--umi-files",
        whitelist.to_str().unwrap(),
        "--max-mismatches",
        "1",
        "--compression-level",
        "1",
    ];
    args.extend_from_slice(extra);
    let cmd = CorrectUmis::try_parse_from(args).expect("parse correct args");
    cmd.execute("fgumi correct").unwrap_or_else(|e| panic!("correct ({tag}) failed: {e:#}"));
}

/// `--threads N` (chain) vs no-`--threads` (single-thread oracle) on a
/// non-vacuous fixture: one exact-match family (kept as-is), one correctable
/// (1-mismatch) family (corrected), and one uncorrectable (far) family
/// (dropped from the main output, since no `--rejects` is given). Exercises
/// the keep/rewrite/reject paths together so the parity check cannot pass on
/// a trivial pass-through.
#[rstest]
#[case::threads1(1)]
#[case::threads4(4)]
fn test_correct_chain_matches_single_threaded(#[case] threads: usize) {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let whitelist = temp_dir.path().join("whitelist.txt");

    let exact = create_umi_family("ACGTACGT", 3, "exact", "AAAAGGGG", 30);
    let corr = create_umi_family("ACGTACGA", 2, "corr", "AAAAGGGG", 30);
    let far = create_umi_family("TTTTTTTT", 2, "far", "AAAAGGGG", 30);
    create_umi_bam(&input_bam, vec![exact, corr, far]);
    create_whitelist(&whitelist, &["ACGTACGT"]);

    let oracle_out = temp_dir.path().join("oracle.bam");
    let chain_out = temp_dir.path().join("chain.bam");
    run_correct(&input_bam, &oracle_out, &whitelist, &["--min-distance", "1"], "oracle");
    let threads_arg = threads.to_string();
    run_correct(
        &input_bam,
        &chain_out,
        &whitelist,
        &["--min-distance", "1", "--threads", &threads_arg],
        "chain",
    );

    // 3 exact + 2 corrected are kept; the 2 far-family records are dropped
    // (no --rejects). This count is computed from the fixture, not derived
    // from the code under test, so the parity check below cannot pass
    // vacuously on two empty/degenerate outputs.
    let expected_kept = 5;
    let (oracle_header, oracle_records) = read_bam_output(&oracle_out);
    let (chain_header, chain_records) = read_bam_output(&chain_out);
    assert_eq!(oracle_records.len(), expected_kept, "oracle output should keep exactly 5 records");
    assert_eq!(
        oracle_records, chain_records,
        "chain (--threads {threads}) output must match the single-threaded oracle record-for-record",
    );
    // Whole-header parity (read_bam_output normalizes the @PG CL, which
    // legitimately differs by --threads): also catches a dropped @SQ/@RG/@HD/@CO.
    assert_eq!(
        oracle_header, chain_header,
        "chain and oracle output headers must match (with @PG CL normalized)",
    );
}

/// The `--rejects` output must also match record-for-record between the
/// chain and oracle paths -- rejects flow through the chain's own
/// `correct_step_with_rejects` branch, a code path the main-output-only
/// parity test above does not exercise.
#[test]
fn test_correct_chain_matches_single_threaded_rejects() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let whitelist = temp_dir.path().join("whitelist.txt");

    let exact = create_umi_family("ACGTACGT", 3, "exact", "AAAAGGGG", 30);
    let far = create_umi_family("TTTTTTTT", 4, "far", "AAAAGGGG", 30);
    create_umi_bam(&input_bam, vec![exact, far]);
    create_whitelist(&whitelist, &["ACGTACGT"]);

    let oracle_out = temp_dir.path().join("oracle.bam");
    let oracle_rejects = temp_dir.path().join("oracle.rejects.bam");
    let chain_out = temp_dir.path().join("chain.bam");
    let chain_rejects = temp_dir.path().join("chain.rejects.bam");

    run_correct(
        &input_bam,
        &oracle_out,
        &whitelist,
        &["--min-distance", "1", "--rejects", oracle_rejects.to_str().unwrap()],
        "oracle",
    );
    run_correct(
        &input_bam,
        &chain_out,
        &whitelist,
        &["--min-distance", "1", "--rejects", chain_rejects.to_str().unwrap(), "--threads", "4"],
        "chain",
    );

    let oracle_rejects_records = read_correct_records(&oracle_rejects);
    assert!(!oracle_rejects_records.is_empty(), "oracle rejects must be non-empty");
    assert_eq!(oracle_rejects_records.len(), 4, "the far family (4 reads) is entirely rejected");
    assert_eq!(
        oracle_rejects_records,
        read_correct_records(&chain_rejects),
        "chain (--threads 4) rejects output must match the oracle record-for-record",
    );
    assert_eq!(
        read_correct_records(&oracle_out),
        read_correct_records(&chain_out),
        "main output must also match between chain and oracle",
    );
}

/// Byte-compares the `--metrics` TSV between the chain and oracle paths.
/// Metrics counting runs through a per-thread accumulator on the chain path
/// vs. a single accumulator on the oracle, even though both eventually call
/// the same `UmiCorrectionMetrics::write_metrics` writer -- so this axis is
/// load-bearing: it is the one test that would catch a per-thread counting
/// divergence that a record-parity check alone would miss.
#[test]
fn test_correct_metrics_files_match_across_modes() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let whitelist = temp_dir.path().join("whitelist.txt");

    let exact = create_umi_family("ACGTACGT", 5, "exact", "AAAAGGGG", 30);
    let corr = create_umi_family("ACGTACGA", 3, "corr", "AAAAGGGG", 30);
    create_umi_bam(&input_bam, vec![exact, corr]);
    create_whitelist(&whitelist, &["ACGTACGT"]);

    let oracle_metrics = temp_dir.path().join("oracle.metrics.tsv");
    let chain_metrics = temp_dir.path().join("chain.metrics.tsv");
    run_correct(
        &input_bam,
        &temp_dir.path().join("oracle.bam"),
        &whitelist,
        &["--min-distance", "1", "--metrics", oracle_metrics.to_str().unwrap()],
        "oracle",
    );
    run_correct(
        &input_bam,
        &temp_dir.path().join("chain.bam"),
        &whitelist,
        &["--min-distance", "1", "--metrics", chain_metrics.to_str().unwrap(), "--threads", "4"],
        "chain",
    );

    // Non-vacuous: the corrected family's one-mismatch match must show up as
    // a real, non-zero count -- not an all-zero table.
    let parsed = fgumi_lib::metrics::UmiCorrectionMetrics::read_metrics(&oracle_metrics)
        .expect("parse oracle metrics");
    assert!(
        parsed.iter().any(|m| m.one_mismatch_matches > 0),
        "metrics must record a non-zero one-mismatch match count: {parsed:?}",
    );

    let oracle_content = fs::read_to_string(&oracle_metrics).expect("read oracle metrics");
    let chain_content = fs::read_to_string(&chain_metrics).expect("read chain metrics");
    assert_eq!(oracle_content, chain_content, "metrics TSV must be byte-identical across modes");
}

/// Build a **paired** template: two mates (R1/R2) sharing query name `name`,
/// both carrying the same `RX` UMI, both mapped primaries. Unlike
/// `create_umi_family` (which emits distinct-named single reads), these two
/// records share one query name, so `GroupByQueryname` groups them as one
/// multi-record template. When such a pair straddles a byte-aligned input
/// record-batch boundary (`DecompressedBlock`s are record-aligned but not
/// template-aligned), it exercises the cross-batch `current_template`
/// carry-over that a single-record-per-name fixture never triggers.
fn create_paired_umi_template(name: &str, umi: &str, pos: i32) -> Vec<RawRecord> {
    [(flags::PAIRED | flags::FIRST_SEGMENT, pos), (flags::PAIRED | flags::LAST_SEGMENT, pos + 8)]
        .into_iter()
        .map(|(flag, mate_pos)| {
            let mut b = SamBuilder::new();
            b.read_name(name.as_bytes())
                .ref_id(0)
                .pos(mate_pos)
                .mapq(60)
                .flags(flag)
                .cigar_ops(&[8 << 4]) // 8M
                .sequence(b"ACGTACGT")
                .qualities(&[30; 8]);
            b.add_string_tag(SamTag::RX, umi.as_bytes());
            b.build()
        })
        .collect()
}

/// The chain's multi-worker path reorders across batches and accumulates
/// per-batch state -- including `GroupByQueryname`'s cross-input-batch
/// `current_template` carry-over, which only fires when a query-name group
/// (a multi-record template) straddles an input record-batch boundary.
/// `create_paired_umi_template` emits `N` two-mate templates (`2 * N`
/// records) so the input spans multiple byte-aligned input record-batches
/// (exercising the carry-over) and comfortably exceeds `GroupByQueryname`'s
/// default template-batch size (exercising the output-emit boundary). A
/// single-record-per-name fixture exercises neither. Compare `--threads 4`
/// against the oracle record-for-record.
#[test]
fn test_correct_chain_matches_single_threaded_multi_batch() {
    const N: usize = 1200;

    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let whitelist = temp_dir.path().join("whitelist.txt");

    let families: Vec<Vec<RawRecord>> = (0..N)
        .map(|i| {
            create_paired_umi_template(
                &format!("fam{i}"),
                "ACGTACGT",
                100 + i32::try_from(i).expect("i fits in i32"),
            )
        })
        .collect();
    create_umi_bam(&input_bam, families);
    create_whitelist(&whitelist, &["ACGTACGT"]);

    let oracle_out = temp_dir.path().join("oracle.bam");
    let chain_out = temp_dir.path().join("chain.bam");
    run_correct(&input_bam, &oracle_out, &whitelist, &["--min-distance", "1"], "oracle");
    run_correct(
        &input_bam,
        &chain_out,
        &whitelist,
        &["--min-distance", "1", "--threads", "4"],
        "chain",
    );

    let oracle_records = read_correct_records(&oracle_out);
    assert_eq!(
        oracle_records.len(),
        N * 2,
        "all {} records ({N} paired templates) should be kept",
        N * 2
    );
    assert_eq!(
        oracle_records,
        read_correct_records(&chain_out),
        "chain (--threads 4) output must match the single-threaded oracle record-for-record across batches",
    );
}

/// `--min-distance 0` parity (proves Task 2A): the chain path routes through
/// `add_correct`'s `CorrectOptions::validate`, which used to bail on
/// `min_distance_diff == 0` before the relaxation. Both paths must succeed
/// and match record-for-record; a panic/error here would mean the "no
/// unguarded `- 1` subtraction" safety analysis backing the relaxation was
/// wrong.
#[test]
fn test_correct_chain_matches_single_threaded_min_distance_zero() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let whitelist = temp_dir.path().join("whitelist.txt");

    let exact = create_umi_family("ACGTACGT", 3, "exact", "AAAAGGGG", 30);
    let corr = create_umi_family("ACGTACGA", 2, "corr", "AAAAGGGG", 30);
    create_umi_bam(&input_bam, vec![exact, corr]);
    create_whitelist(&whitelist, &["ACGTACGT"]);

    let oracle_out = temp_dir.path().join("oracle.bam");
    let chain_out = temp_dir.path().join("chain.bam");
    run_correct(&input_bam, &oracle_out, &whitelist, &["--min-distance", "0"], "oracle");
    run_correct(
        &input_bam,
        &chain_out,
        &whitelist,
        &["--min-distance", "0", "--threads", "4"],
        "chain",
    );

    let oracle_records = read_correct_records(&oracle_out);
    assert_eq!(oracle_records.len(), 5, "both families should be kept under --min-distance 0");
    assert_eq!(
        oracle_records,
        read_correct_records(&chain_out),
        "--min-distance 0 output must match between the chain and single-threaded oracle",
    );
}

/// The chain (`--threads`) path rejects empty input under a positive
/// `--min-corrected`, where the legacy single-threaded oracle silently passes
/// (its `0 / 0 = NaN` kept-ratio compares false against the floor -- a latent
/// bug). This is an intentional correctness improvement of the chain path, not
/// a regression: an empty run cannot meet a positive minimum, so the chain's
/// explicit error is the right behavior. Pin it so it is not lost; the legacy
/// engine's NaN-pass is left to the follow-up PR that removes that engine.
#[test]
fn test_correct_chain_rejects_empty_input_with_min_corrected() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let whitelist = temp_dir.path().join("whitelist.txt");

    // Header-only BAM: no records at all.
    create_umi_bam(&input_bam, vec![]);
    create_whitelist(&whitelist, &["ACGTACGT"]);

    let cmd = CorrectUmis::try_parse_from([
        "correct",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        temp_dir.path().join("chain.bam").to_str().unwrap(),
        "--umi-files",
        whitelist.to_str().unwrap(),
        "--max-mismatches",
        "1",
        "--min-distance",
        "1",
        "--min-corrected",
        "0.5",
        "--compression-level",
        "1",
        "--threads",
        "4",
    ])
    .expect("parse correct args");
    let err = cmd
        .execute("fgumi correct")
        .expect_err("chain path must reject empty input under a positive --min-corrected");
    assert!(err.to_string().contains("No reads were processed"), "unexpected error: {err:#}");
}

/// Write a header-less UMI-tagged BAM (no `@HD`) -- the shape
/// `ChainBuilder::new`'s `ensure_hd_record` call (Task 2B) must synthesize
/// `@HD` for. Mirrors `write_headerless_consensus_bam` in
/// `test_filter_command.rs`.
fn write_headerless_umi_bam(path: &std::path::Path) {
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
    for record in create_umi_family("ACGTACGT", 3, "read", "AAAAGGGG", 30) {
        writer
            .write_alignment_record(&header, &to_record_buf(&record))
            .expect("Failed to write record");
    }
    writer.try_finish().expect("Failed to finish BAM");
}

/// `@HD` synthesis parity (proves Task 2B): correct enforces no
/// query-grouped guard, so header-less input flows through to completion on
/// both paths. Before the `ChainBuilder::new` fix, the chain path would have
/// emitted a header without `@HD` while the oracle synthesized one -- a real
/// output-byte divergence.
#[rstest]
#[case::threads1(1)]
#[case::threads4(4)]
fn test_correct_chain_synthesizes_hd_for_headerless_input(#[case] threads: usize) {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let whitelist = temp_dir.path().join("whitelist.txt");
    write_headerless_umi_bam(&input_bam);
    create_whitelist(&whitelist, &["ACGTACGT"]);

    let oracle_out = temp_dir.path().join("oracle.bam");
    let chain_out = temp_dir.path().join("chain.bam");
    run_correct(&input_bam, &oracle_out, &whitelist, &["--min-distance", "1"], "oracle");
    let threads_arg = threads.to_string();
    run_correct(
        &input_bam,
        &chain_out,
        &whitelist,
        &["--min-distance", "1", "--threads", &threads_arg],
        "chain",
    );

    let (oracle_header, _) = read_bam_output(&oracle_out);
    let (chain_header, _) = read_bam_output(&chain_out);
    assert!(
        oracle_header.header().is_some(),
        "oracle output must synthesize @HD for headerless input",
    );
    assert!(
        chain_header.header().is_some(),
        "chain output must synthesize @HD for headerless input",
    );
    assert_eq!(
        chain_header.header(),
        oracle_header.header(),
        "chain (--threads {threads}) must synthesize the same @HD as the oracle for headerless input",
    );
}

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
use fgumi_raw_bam::RawRecord;
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

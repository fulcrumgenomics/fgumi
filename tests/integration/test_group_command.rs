//! Integration tests for the group command with new metrics infrastructure.
//!
//! These tests invoke `GroupReadsByUmi::execute()` in-process.

use clap::Parser;
use fgoxide::io::DelimFile;
use fgumi_lib::commands::command::Command;
use fgumi_lib::commands::group::GroupReadsByUmi;
use fgumi_lib::metrics::UmiGroupingMetrics;
use noodles::bam;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use rstest::rstest;
use std::fs;
use std::path::PathBuf;
use tempfile::TempDir;

use crate::helpers::assertions::assert_bam_sorted;
use crate::helpers::bam_generator::{
    create_minimal_header, create_query_grouped_header, create_umi_family,
    create_umi_family_at_pos, to_record_buf,
};
use fgumi_lib::sam::SamTag;
use fgumi_raw_bam::{RawRecord, SamBuilder, flags};

/// Test that the group command properly writes metrics in the new format.
#[test]
fn test_group_command_writes_new_metrics() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let metrics_file = temp_dir.path().join("metrics.txt");

    // Create test BAM file with multiple UMI families
    create_test_input_bam(&input_bam);

    // Run the group command
    let cmd = GroupReadsByUmi::try_parse_from([
        "group",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--strategy",
        "identity",
        "--edits",
        "0",
        "--grouping-metrics",
        metrics_file.to_str().unwrap(),
        "--compression-level",
        "1",
    ])
    .expect("failed to parse group args");
    cmd.execute("fgumi group").expect("Failed to run group command");
    assert!(output_bam.exists(), "Output BAM not created");
    assert!(metrics_file.exists(), "Metrics file not created");

    // group's `@HD` advertises SS:template-coordinate ("Output is always written in
    // template-coordinate order"); gate that claim against the bytes it wrote.
    assert_bam_sorted(&output_bam, "template-coordinate", Some("mi"));

    // Read and validate metrics
    let metrics: Vec<UmiGroupingMetrics> =
        DelimFile::default().read_tsv(&metrics_file).expect("Failed to read metrics file");

    assert_eq!(metrics.len(), 1, "Expected exactly one metrics record");

    let metric = &metrics[0];

    // Assert the exact serialized fgbio columns against `create_test_input_bam`'s known
    // structure: three families of 10 + 5 + 15 single-end reads, all PF, all mapq 60, all
    // carrying a valid 8bp UMI with no Ns — so every record is accepted and nothing is
    // discarded. The `.grouping_metrics.txt` matches fgbio's 5-column `UmiGroupingMetric`;
    // `total_records`/`unique_molecule_ids` are fgumi-internal (`#[serde(skip)]`) and read
    // back as zero, so they are not asserted here.
    assert_eq!(metric.accepted_records, 30, "all 30 input records should be accepted");
    assert_eq!(metric.discarded_non_pf, 0);
    assert_eq!(metric.discarded_poor_alignment, 0);
    assert_eq!(metric.discarded_ns_in_umi, 0);
    assert_eq!(metric.discarded_umi_too_short, 0);
}

/// Test that the group command handles UMIs with N bases correctly.
#[test]
fn test_group_command_rejects_n_bases_in_umi() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let metrics_file = temp_dir.path().join("metrics.txt");

    // Create BAM with UMIs containing N bases
    create_test_bam_with_n_umis(&input_bam);

    let cmd = GroupReadsByUmi::try_parse_from([
        "group",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--strategy",
        "identity",
        "--edits",
        "0",
        "--grouping-metrics",
        metrics_file.to_str().unwrap(),
        "--compression-level",
        "1",
    ])
    .expect("failed to parse group args");
    cmd.execute("fgumi group").expect("Failed to run group command");

    // Read metrics and verify N rejection tracking
    let metrics: Vec<UmiGroupingMetrics> =
        DelimFile::default().read_tsv(&metrics_file).expect("Failed to read metrics");

    let metric = &metrics[0];
    assert!(metric.discarded_ns_in_umi > 0, "Should have discarded reads with N in UMI");
}

/// `group`'s metrics file is fgbio's `UmiGroupingMetric` and must stay exactly
/// five columns with fgbio's spellings, whatever fgumi tracks internally. The
/// struct-level test in `fgumi-metrics` pins the serialization; this pins what the
/// command actually writes to disk, so an internal refactor of the filter
/// accounting cannot change the file fgbio has to be able to read.
#[test]
fn test_grouping_metrics_header_is_fgbio_five_columns() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let metrics_file = temp_dir.path().join("metrics.txt");

    create_test_input_bam(&input_bam);

    let cmd = GroupReadsByUmi::try_parse_from([
        "group",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--strategy",
        "identity",
        "--edits",
        "0",
        "--grouping-metrics",
        metrics_file.to_str().unwrap(),
        "--compression-level",
        "1",
    ])
    .expect("failed to parse group args");
    cmd.execute("fgumi group").expect("Failed to run group command");

    let content = fs::read_to_string(&metrics_file).expect("Failed to read metrics file");
    assert_eq!(
        content.lines().next(),
        Some(
            "accepted_sam_records\tdiscarded_non_pf\tdiscarded_poor_alignment\t\
             discarded_ns_in_umi\tdiscarded_umis_to_short"
        )
    );
}

/// Helper function to create a test BAM file with multiple UMI families.
fn create_test_input_bam(path: &PathBuf) {
    let header = create_minimal_header("chr1", 10000);

    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));

    writer.write_header(&header).expect("Failed to write header");

    // Create 3 UMI families
    let family1 = create_umi_family("AAAAAAAA", 10, "family1", "ACGTACGT", 30);
    let family2 = create_umi_family("CCCCCCCC", 5, "family2", "TGCATGCA", 30);
    let family3 = create_umi_family("GGGGGGGG", 15, "family3", "ATCGATCG", 30);

    for record in family1.iter().chain(family2.iter()).chain(family3.iter()) {
        writer
            .write_alignment_record(&header, &to_record_buf(record))
            .expect("Failed to write record");
    }

    writer.try_finish().expect("Failed to finish BAM");
}

/// Helper function to create a BAM file with UMIs containing N bases.
fn create_test_bam_with_n_umis(path: &PathBuf) {
    let header = create_minimal_header("chr1", 10000);

    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));

    writer.write_header(&header).expect("Failed to write header");

    // Create families with good and bad UMIs
    let good_family = create_umi_family("AAAAAAAA", 5, "good", "ACGTACGT", 30);
    let bad_family = create_umi_family("NNNNNNNN", 3, "bad", "TGCATGCA", 30);

    for record in good_family.iter().chain(bad_family.iter()) {
        writer
            .write_alignment_record(&header, &to_record_buf(record))
            .expect("Failed to write record");
    }

    writer.try_finish().expect("Failed to finish BAM");
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

/// Write a single UMI family large enough to span more than one BGZF block (all
/// reads sit at one position on reference 0, so the input is template-coordinate
/// sorted), so the corrupted last block is a record block, not the header's.
fn create_multiblock_group_bam(path: &PathBuf) {
    let header = create_minimal_header("chr1", 10000);
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");
    for record in &create_umi_family("AAAAAAAA", 3000, "read", "ACGTACGT", 30) {
        writer
            .write_alignment_record(&header, &to_record_buf(record))
            .expect("Failed to write record");
    }
    writer.try_finish().expect("Failed to finish BAM");
}

/// Build group's argv for the single-threaded fast path (no `--threads`),
/// appending any extra flags (e.g. `--no-check-crc`).
fn group_args<'a>(input: &'a str, output: &'a str, extra: &[&'a str]) -> Vec<&'a str> {
    let mut args = vec![
        "group",
        "--input",
        input,
        "--output",
        output,
        "--strategy",
        "identity",
        "--edits",
        "0",
        "--compression-level",
        "1",
    ];
    args.extend_from_slice(extra);
    args
}

/// `--no-check-crc` must let group's single-threaded reader accept a corrupted
/// BGZF CRC32: it decodes through fgumi-bgzf, honoring the flag (#800). Against
/// the noodles reader this path used before, this could not pass.
#[test]
fn test_group_no_check_crc_accepts_corrupted_crc() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    create_multiblock_group_bam(&input_bam);
    // Baseline: group the intact file to pin the expected records.
    let baseline = run_group_records(temp_dir.path(), "baseline", &input_bam, &[]);
    assert!(!baseline.is_empty(), "baseline run must produce grouped records");

    // Corrupt only the last block's CRC32 footer, then re-group with
    // --no-check-crc: it must accept the file AND produce records identical to
    // the intact run (a bit-flipped CRC footer leaves the payload intact, so a
    // correct skip-CRC decode is byte-for-byte identical — not merely non-empty).
    corrupt_last_block_crc(&input_bam);
    let accepted = run_group_records(temp_dir.path(), "nocrc", &input_bam, &["--no-check-crc"]);
    assert_eq!(
        baseline, accepted,
        "--no-check-crc must accept the corrupted CRC and yield records identical to the intact run",
    );
}

/// Default (verify-on for file input) must reject the same corrupted CRC32.
#[test]
fn test_group_rejects_corrupted_crc_by_default() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    create_multiblock_group_bam(&input_bam);
    corrupt_last_block_crc(&input_bam);

    let cmd = GroupReadsByUmi::try_parse_from(group_args(
        input_bam.to_str().unwrap(),
        output_bam.to_str().unwrap(),
        &[],
    ))
    .expect("failed to parse group args");
    let err = cmd
        .execute("fgumi group")
        .expect_err("default (verify-on for file input) must reject a corrupted BGZF CRC32");
    let message = format!("{err:#}");
    assert!(message.to_uppercase().contains("CRC32"), "error should mention CRC32: {message}");
}

/// `--check-crc` must also reject the corrupted CRC32 (forces verification on).
#[test]
fn test_group_check_crc_rejects_corrupted_crc() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    create_multiblock_group_bam(&input_bam);
    corrupt_last_block_crc(&input_bam);

    let cmd = GroupReadsByUmi::try_parse_from(group_args(
        input_bam.to_str().unwrap(),
        output_bam.to_str().unwrap(),
        &["--check-crc"],
    ))
    .expect("failed to parse group args");
    let err =
        cmd.execute("fgumi group").expect_err("--check-crc must reject a corrupted BGZF CRC32");
    let message = format!("{err:#}");
    assert!(message.to_uppercase().contains("CRC32"), "error should mention CRC32: {message}");
}

// ============================================================================
// Chain-backed `--threads N` path: parity with the single-threaded fast path
// (the `group` pilot rewire). The `--threads` path now runs on the declarative
// chain builder; these tests pin the behaviors the single-threaded oracle
// cannot see: the chain's CRC policy and its record-ordering acceptance, plus
// cross-mode output parity.
// ============================================================================

/// Read the complete `RecordBuf` for every output record, in output order,
/// asserting an `MI:Z` tag is present on each.
///
/// Positional `Vec` equality across two runs catches both re-ordering and any
/// per-record difference (MI re-numbering, dropped/added records).
fn read_group_records(path: &std::path::Path) -> Vec<noodles::sam::alignment::RecordBuf> {
    use noodles::sam::alignment::record_buf::data::field::Value;
    let mut reader = bam::io::Reader::new(fs::File::open(path).expect("open output bam"));
    let header = reader.read_header().expect("read output header");
    let mi = fgumi_lib::sam::SamTag::MI.to_noodles_tag();
    reader
        .record_bufs(&header)
        .map(|r| {
            let r = r.expect("read output record");
            // Group must assign an MI:Z to every accepted output record; assert
            // that invariant here, then return the COMPLETE record so the parity
            // checks compare every BAM field, not only QNAME/flags/MI (a field
            // that changed outside those three would otherwise slip through).
            let qname = r.name().expect("output record has a name").to_string();
            match r.data().get(&mi) {
                Some(Value::String(_)) => {}
                other => panic!("expected MI:Z string on output record {qname}, got {other:?}"),
            }
            r
        })
        .collect()
}

/// Run `group` on `input`, returning the output records. `extra` is appended to
/// the base argv (e.g. `["--threads", "1"]`, `["--allow-unmapped"]`).
fn run_group_records(
    dir: &std::path::Path,
    tag: &str,
    input: &std::path::Path,
    extra: &[&str],
) -> Vec<noodles::sam::alignment::RecordBuf> {
    let output = dir.join(format!("{tag}.bam"));
    let cmd = GroupReadsByUmi::try_parse_from(group_args(
        input.to_str().unwrap(),
        output.to_str().unwrap(),
        extra,
    ))
    .expect("failed to parse group args");
    cmd.execute("fgumi group").unwrap_or_else(|e| panic!("group ({tag}) failed: {e:#}"));
    read_group_records(&output)
}

/// The chain-backed `--threads N` path must honor the same CRC-verification
/// policy as the single-threaded fast path — the plumbing that threads
/// `effective_check_crc()` into the chain's BGZF decode. Before that plumbing
/// the chain hardcoded verify-on, so `--no-check-crc --threads N` would have
/// wrongly rejected a corrupted-CRC file. The single-threaded CRC tests above
/// cannot see this path.
#[rstest]
#[case::no_check_crc_accepts(&["--no-check-crc"], true)]
#[case::default_rejects(&[], false)]
#[case::check_crc_rejects(&["--check-crc"], false)]
fn test_group_threaded_crc_policy(#[case] extra: &[&str], #[case] expect_ok: bool) {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");

    create_multiblock_group_bam(&input_bam);
    // Baseline: group the intact file through the same chain (--threads 1) so the
    // accept case can assert record identity, not merely non-emptiness — a chain-path
    // regression could silently drop the corrupted last block's records yet still
    // write a non-empty file.
    let baseline = expect_ok
        .then(|| run_group_records(temp_dir.path(), "baseline", &input_bam, &["--threads", "1"]));
    corrupt_last_block_crc(&input_bam);

    let mut args = extra.to_vec();
    args.extend_from_slice(&["--threads", "1"]);
    let cmd = GroupReadsByUmi::try_parse_from(group_args(
        input_bam.to_str().unwrap(),
        output_bam.to_str().unwrap(),
        &args,
    ))
    .expect("failed to parse group args");
    let result = cmd.execute("fgumi group");

    if expect_ok {
        result.expect("--no-check-crc --threads must accept a corrupted BGZF CRC32 and complete");
        let records = read_group_records(&output_bam);
        assert_eq!(
            baseline.expect("baseline is captured for the accept case"),
            records,
            "--no-check-crc --threads must accept the corrupted CRC and yield records \
             identical to the intact run, not merely a non-empty file",
        );
    } else {
        let err = result.expect_err("chain path must reject a corrupted BGZF CRC32");
        let message = format!("{err:#}");
        assert!(message.to_uppercase().contains("CRC32"), "error should mention CRC32: {message}");
    }
}

/// The chain-backed `--threads 1` output must match the single-threaded fast
/// path record-for-record on a normal template-coordinate input. This is the
/// cross-mode parity oracle the within-mode determinism tests cannot provide.
#[test]
fn test_group_chain_matches_single_threaded() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    create_test_input_bam(&input_bam);

    let single = run_group_records(temp_dir.path(), "single", &input_bam, &[]);
    let chained = run_group_records(temp_dir.path(), "chained", &input_bam, &["--threads", "1"]);

    assert!(
        !single.is_empty(),
        "sanity: expected grouped records so the parity check is not vacuous"
    );
    assert_eq!(
        single, chained,
        "chain (--threads 1) output must match the single-threaded fast path record-for-record",
    );
}

/// A query-grouped (GO:query, not template-coordinate) input under
/// `--allow-unmapped` must be accepted by the chain and produce output
/// identical to the single-threaded path. The chain builder previously required
/// `SO:queryname` and would have rejected this input; sharing the classifier
/// (`require_group_input_ordering`) fixed that.
#[test]
fn test_group_allow_unmapped_query_grouped_chain_parity() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");

    // Query-grouped (GO:query) input, not template-coordinate sorted.
    let header = create_query_grouped_header("chr1", 10000);
    let mut writer =
        bam::io::Writer::new(fs::File::create(&input_bam).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");
    let family1 = create_umi_family("AAAAAAAA", 10, "family1", "ACGTACGT", 30);
    let family2 = create_umi_family("CCCCCCCC", 5, "family2", "TGCATGCA", 30);
    for record in family1.iter().chain(family2.iter()) {
        writer
            .write_alignment_record(&header, &to_record_buf(record))
            .expect("Failed to write record");
    }
    writer.try_finish().expect("Failed to finish BAM");

    let single = run_group_records(temp_dir.path(), "single", &input_bam, &["--allow-unmapped"]);
    let chained = run_group_records(
        temp_dir.path(),
        "chained",
        &input_bam,
        &["--allow-unmapped", "--threads", "1"],
    );

    // This fixture/flag combo is content-pinned by no other test, so guard
    // against a vacuous pass where both paths regressed to empty/all-dropped
    // output: 10 + 5 input reads across two UMI families, all accepted.
    assert_eq!(single.len(), 15, "expected all 15 query-grouped reads in the output");
    assert_eq!(
        single, chained,
        "chain and single-threaded --allow-unmapped output must match on a query-grouped input",
    );
}

/// The chain's finalize hook writes metrics via `write_metrics_for_chain`, a
/// separate path from the single-threaded `write_all_metrics`. This pins that
/// the two produce identical metrics files (grouping-metrics and the
/// family-size histogram) on the same input, so the two orchestrations cannot
/// drift on secondary output. Metrics files are deterministic TSV count tables
/// with no embedded paths/timestamps, so a byte comparison is stable.
#[test]
fn test_group_metrics_files_match_across_modes() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    create_test_input_bam(&input_bam);

    // Run `group` for one mode, returning (grouping-metrics, family-size-histogram) contents.
    let run = |tag: &str, extra: &[&str]| -> (String, String) {
        let output = temp_dir.path().join(format!("{tag}.bam"));
        let grouping = temp_dir.path().join(format!("{tag}.grouping.txt"));
        let famsize = temp_dir.path().join(format!("{tag}.famsize.txt"));
        let mut args: Vec<&str> = vec![
            "--grouping-metrics",
            grouping.to_str().unwrap(),
            "--family-size-histogram",
            famsize.to_str().unwrap(),
        ];
        args.extend_from_slice(extra);
        let cmd = GroupReadsByUmi::try_parse_from(group_args(
            input_bam.to_str().unwrap(),
            output.to_str().unwrap(),
            &args,
        ))
        .expect("failed to parse group args");
        cmd.execute("fgumi group").unwrap_or_else(|e| panic!("group ({tag}) failed: {e:#}"));
        (
            fs::read_to_string(&grouping).expect("read grouping metrics"),
            fs::read_to_string(&famsize).expect("read family-size histogram"),
        )
    };

    let single = run("single", &[]);
    let chained = run("chained", &["--threads", "1"]);

    // Pin a concrete value so the equality checks below cannot pass on two empty
    // or degenerate files: create_test_input_bam has 30 accepted records, and
    // `accepted_sam_records` is the first grouping-metrics column.
    assert!(
        single.0.lines().nth(1).is_some_and(|row| row.starts_with("30\t")),
        "grouping-metrics should report 30 accepted records; got:\n{}",
        single.0
    );

    assert_eq!(
        single.0, chained.0,
        "grouping-metrics file must be identical across single-threaded and chain modes",
    );
    assert_eq!(
        single.1, chained.1,
        "family-size histogram must be identical across single-threaded and chain modes",
    );
}

/// Build a template-coordinate-sorted BAM with `n` single-read families at `n`
/// distinct ascending positions, so the input spans `n` position groups and —
/// for `n` past the pipeline's 500-template batch size — multiple batches.
fn create_multi_batch_input_bam(path: &std::path::Path, n: usize) {
    let header = create_minimal_header("chr1", 10_000_000);
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");
    for i in 0..n {
        for record in
            &create_umi_family_at_pos("AAAAAAAA", 1, &format!("fam{i}"), "ACGTACGT", 30, i + 1)
        {
            writer
                .write_alignment_record(&header, &to_record_buf(record))
                .expect("Failed to write record");
        }
    }
    writer.try_finish().expect("Failed to finish BAM");
}

/// The chain's multi-worker path (`--threads N>1`) reorders across batches and
/// accumulates per-batch `MoleculeId` offsets — machinery the `--threads 1` parity
/// tests do not exercise. Compare `--threads 4` output record-for-record against
/// the single-threaded oracle on an input spanning many position groups and
/// multiple pipeline batches (> the 500-template batch size).
#[test]
fn test_group_chain_threads4_matches_single_threaded_multi_batch() {
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    create_multi_batch_input_bam(&input_bam, 600);

    let single = run_group_records(temp_dir.path(), "single", &input_bam, &[]);
    let threads4 = run_group_records(temp_dir.path(), "threads4", &input_bam, &["--threads", "4"]);

    // 600 positions × 1 read → 600 output records spanning multiple batches.
    assert_eq!(single.len(), 600, "expected 600 grouped records (one per position group)");
    assert_eq!(
        single, threads4,
        "chain (--threads 4) output must match the single-threaded oracle record-for-record across batches",
    );
}

/// End-to-end regression for issue #901: a secondary/supplementary read carrying
/// a `tc` tag (stamped by `fgumi zipper` for template-coordinate ordering) is
/// documented as filtered from grouping, so it must not change the molecule
/// count. Before the fix, such a read received a resolved position key and, being
/// interleaved between two same-UMI reads of one molecule, split that molecule
/// across position groups — inflating the reported molecule count and diverging
/// from fgbio. This pins the reported symptom (distinct MI count), not just the
/// internal position-group shape the unit tests cover.
///
/// Records are written in a deliberately split-inducing order (the tc-keyed
/// supplementary sits between the two primaries); `group` trusts the declared
/// template-coordinate order and does not re-sort, so this reproduces the
/// interleaving that real template-coordinate sort can produce.
fn create_supplementary_split_bam(path: &PathBuf) {
    let header = create_minimal_header("chr1", 20000);
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");

    // Two single-read "molecules" sharing one UMI at the same position — one
    // molecule that must group together.
    let mut records = Vec::new();
    records.extend(create_umi_family_at_pos("AAAAAAAA", 1, "b1", "ACGTACGT", 30, 200));

    // A supplementary read carrying the same UMI and a `tc` tag whose resolved
    // template coordinate (pos 1) differs from the primaries' (pos 200).
    let mut supp = SamBuilder::new();
    supp.read_name(b"supp")
        .ref_id(0)
        .pos(4999)
        .mapq(60)
        .flags(flags::SUPPLEMENTARY)
        .cigar_ops(&[8u32 << 4]) // 8M
        .sequence(b"ACGTACGT")
        .qualities(&[30u8; 8]);
    supp.add_string_tag(SamTag::RX, b"AAAAAAAA");
    supp.add_array_i32(SamTag::TC, &[0, 1, 0, 0, 1, 0]);
    records.push(supp.build());

    records.extend(create_umi_family_at_pos("AAAAAAAA", 1, "b2", "ACGTACGT", 30, 200));

    for record in &records {
        writer
            .write_alignment_record(&header, &to_record_buf(record))
            .expect("Failed to write record");
    }
    writer.try_finish().expect("Failed to finish BAM");
}

/// Build the paired-end analogue of [`create_supplementary_split_bam`]: two
/// paired templates sharing one UMI at the same 5' positions (one molecule), with
/// a `tc`-keyed supplementary interleaved between them. This mirrors real client
/// data, where the reads are paired and carry `MC` tags.
fn create_paired_supplementary_split_bam(path: &PathBuf) {
    let header = create_minimal_header("chr1", 20000);
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");

    // One paired template (R1 fwd @200, R2 rev @400) with MC tags, sharing a UMI.
    let paired = |name: &[u8]| -> Vec<RawRecord> {
        let mut r1 = SamBuilder::new();
        r1.read_name(name)
            .ref_id(0)
            .pos(199)
            .mapq(60)
            .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE)
            .mate_ref_id(0)
            .mate_pos(399)
            .cigar_ops(&[8u32 << 4])
            .sequence(b"ACGTACGT")
            .qualities(&[30u8; 8])
            .add_string_tag(SamTag::RX, b"AAAAAAAA")
            .add_string_tag(SamTag::MC, b"8M");
        let mut r2 = SamBuilder::new();
        r2.read_name(name)
            .ref_id(0)
            .pos(399)
            .mapq(60)
            .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE)
            .mate_ref_id(0)
            .mate_pos(199)
            .cigar_ops(&[8u32 << 4])
            .sequence(b"ACGTACGT")
            .qualities(&[30u8; 8])
            .add_string_tag(SamTag::RX, b"AAAAAAAA")
            .add_string_tag(SamTag::MC, b"8M");
        vec![r1.build(), r2.build()]
    };

    // A supplementary carrying the same UMI and a `tc` tag whose resolved
    // template coordinate differs from the primaries', written between the two
    // templates so it would break contiguity if not filtered.
    let mut supp = SamBuilder::new();
    supp.read_name(b"supp")
        .ref_id(0)
        .pos(4999)
        .mapq(60)
        .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::SUPPLEMENTARY)
        .mate_ref_id(0)
        .mate_pos(399)
        .cigar_ops(&[8u32 << 4])
        .sequence(b"ACGTACGT")
        .qualities(&[30u8; 8])
        .add_string_tag(SamTag::RX, b"AAAAAAAA")
        .add_array_i32(SamTag::TC, &[0, 1, 0, 0, 1, 0]);

    let mut records = paired(b"pair1");
    records.push(supp.build());
    records.extend(paired(b"pair2"));

    for record in &records {
        writer
            .write_alignment_record(&header, &to_record_buf(record))
            .expect("Failed to write record");
    }
    writer.try_finish().expect("Failed to finish BAM");
}

/// Shared oracle for the #901 regression: grouped output must contain exactly one
/// molecule (one distinct MI) and no secondary/supplementary records.
fn assert_single_molecule_no_secondary(
    records: &[noodles::sam::alignment::RecordBuf],
    expected_qnames: &[&str],
) {
    use noodles::sam::alignment::record_buf::data::field::Value;
    let mi = SamTag::MI.to_noodles_tag();
    let distinct_mi: std::collections::BTreeSet<String> = records
        .iter()
        .filter_map(|r| match r.data().get(&mi) {
            Some(Value::String(s)) => Some(s.to_string()),
            _ => None,
        })
        .collect();
    assert_eq!(
        distinct_mi.len(),
        1,
        "the tc-keyed supplementary must not split the single molecule (got MIs {distinct_mi:?})"
    );
    // Assert record identity, not just the count: a regression that dropped one
    // of the two primaries would still leave a single MI and no flagged record,
    // so pin the exact set of output QNAMEs.
    let output_qnames: std::collections::BTreeSet<String> =
        records.iter().filter_map(|r| r.name().map(ToString::to_string)).collect();
    let expected: std::collections::BTreeSet<String> =
        expected_qnames.iter().map(|s| (*s).to_string()).collect();
    assert_eq!(
        output_qnames, expected,
        "both primaries must survive: expected QNAMEs {expected:?}, got {output_qnames:?}"
    );
    for r in records {
        let flag = r.flags();
        assert!(
            !flag.is_secondary() && !flag.is_supplementary(),
            "no secondary/supplementary record should appear in grouped output"
        );
    }
}

/// Regression for issue #901 across both the single-threaded fast path and the
/// `--threads N` chain path, for single-end and paired input. A `tc`-keyed
/// supplementary interleaved between two same-UMI reads of one molecule must not
/// split it or appear in the output, on every path.
#[rstest]
#[case::single_end(false)]
#[case::paired(true)]
fn test_group_supplementary_does_not_inflate_molecule_count(
    #[case] paired: bool,
    #[values(&[] as &[&str], &["--threads", "2"])] extra: &[&str],
) {
    let temp_dir = TempDir::new().expect("create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    if paired {
        create_paired_supplementary_split_bam(&input_bam);
    } else {
        create_supplementary_split_bam(&input_bam);
    }

    let records = run_group_records(temp_dir.path(), "supp_split", &input_bam, extra);
    // Single-end input names its two primaries `b1_0`/`b2_0`; the paired input
    // names its two templates `pair1`/`pair2` (R1 and R2 share the QNAME).
    let expected_qnames: &[&str] = if paired { &["pair1", "pair2"] } else { &["b1_0", "b2_0"] };
    assert_single_molecule_no_secondary(&records, expected_qnames);
}

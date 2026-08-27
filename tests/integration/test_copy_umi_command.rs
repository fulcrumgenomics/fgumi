//! End-to-end tests for the `copy-umi` command.
//!
//! These invoke `CopyUmi::execute()` in-process and validate that the UMI in the
//! read name is copied to the `RX` tag with fgbio `CopyUmiFromReadName` parity:
//! last-field location, `+`→`-` joining, `r`-prefix reverse-complementation (and
//! its strip-only override), optional name trimming, tag overwrite/fail policy,
//! and fail-fast on malformed names.

use clap::Parser;
use fgumi_lib::commands::command::Command;
use fgumi_lib::commands::copy_umi::CopyUmi;
use fgumi_raw_bam::{RawRecord, SamBuilder};
use rstest::rstest;
use std::collections::HashMap;
use std::path::Path;
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, transcode_bam_to_sam, write_bam};

/// A minimal mapped record carrying the given read name (sequence/quals are fixed
/// filler — `copy-umi` never touches them).
fn record_named(name: &str) -> RawRecord {
    let mut b = SamBuilder::new();
    b.read_name(name.as_bytes())
        .sequence(b"ACGT")
        .qualities(&[30; 4])
        .flags(0)
        .ref_id(0)
        .pos(99)
        .mapq(60)
        .cigar_ops(&[4 << 4]);
    b.build()
}

/// As [`record_named`], but pre-populated with an `RX` tag to exercise the
/// overwrite / fail-if-present policy.
fn record_named_with_rx(name: &str, rx: &str) -> RawRecord {
    let mut b = SamBuilder::new();
    b.read_name(name.as_bytes())
        .sequence(b"ACGT")
        .qualities(&[30; 4])
        .flags(0)
        .ref_id(0)
        .pos(99)
        .mapq(60)
        .cigar_ops(&[4 << 4])
        .add_string_tag(fgumi_lib::sam::SamTag::RX, rx.as_bytes());
    b.build()
}

/// Write `records` to a BAM in `dir` and return its path.
fn write_input(dir: &Path, records: &[RawRecord]) -> std::path::PathBuf {
    let path = dir.join("in.bam");
    let header = create_minimal_header("chr1", 10_000);
    write_bam(&path, &header, records);
    path
}

/// Parse + run `copy-umi -i input -o output [extra...]`.
fn run_copy_umi(input: &Path, output: &Path, extra: &[&str]) -> anyhow::Result<()> {
    let mut args: Vec<String> = vec![
        "copy-umi".into(),
        "-i".into(),
        input.display().to_string(),
        "-o".into(),
        output.display().to_string(),
    ];
    args.extend(extra.iter().map(|s| (*s).to_string()));
    let cmd = CopyUmi::try_parse_from(args).expect("failed to parse copy-umi args");
    cmd.execute("copy-umi test")
}

/// Read (`name`, `RX`) for every record in a BAM by transcoding to SAM text.
fn read_name_and_rx(dir: &Path, bam: &Path) -> Vec<(String, Option<String>)> {
    let sam = dir.join("out.sam");
    transcode_bam_to_sam(bam, &sam);
    std::fs::read_to_string(&sam)
        .expect("read sam")
        .lines()
        .filter(|l| !l.starts_with('@'))
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            let name = fields[0].to_string();
            let rx = fields.iter().find_map(|f| f.strip_prefix("RX:Z:").map(str::to_string));
            (name, rx)
        })
        .collect()
}

// ============================================================================
// fgbio-parity: read name -> RX, reverse-complement ON (default)
// ============================================================================

#[rstest]
#[case::single("inst:1:FC:1:1101:5:7:ACGT", "ACGT")]
#[case::short_name("1:2:3:GGGG", "GGGG")]
#[case::dual_plus_to_hyphen("blah:AAAA+CCCC", "AAAA-CCCC")]
#[case::already_hyphen("blah:AAAA-CCCC", "AAAA-CCCC")]
#[case::r_prefixed_revcomped("blah:rAAAA", "TTTT")]
#[case::r_only_prefixed_segment("blah:rAAAA+CCCC", "TTTT-CCCC")]
fn copies_umi_to_rx(#[case] name: &str, #[case] expected_rx: &str) {
    let dir = TempDir::new().unwrap();
    let input = write_input(dir.path(), &[record_named(name)]);
    let output = dir.path().join("out.bam");
    run_copy_umi(&input, &output, &[]).expect("copy-umi failed");

    let rows = read_name_and_rx(dir.path(), &output);
    assert_eq!(rows.len(), 1);
    assert_eq!(rows[0].0, name, "read name should be unchanged without --remove-umi");
    assert_eq!(rows[0].1.as_deref(), Some(expected_rx));
}

// ============================================================================
// Reverse-complement OFF (strip-only): fgpyo / UMI-tools behavior
// ============================================================================

#[rstest]
#[case::r_stripped_not_revcomped("blah:rAAAA", "AAAA")]
#[case::r_stripped_dual("blah:rAAAA+CCCC", "AAAA-CCCC")]
fn strip_only_does_not_reverse_complement(#[case] name: &str, #[case] expected_rx: &str) {
    let dir = TempDir::new().unwrap();
    let input = write_input(dir.path(), &[record_named(name)]);
    let output = dir.path().join("out.bam");
    run_copy_umi(&input, &output, &["--reverse-complement-r-umis", "false"])
        .expect("copy-umi failed");

    let rows = read_name_and_rx(dir.path(), &output);
    assert_eq!(rows.len(), 1, "one input record must yield exactly one output record");
    assert_eq!(rows[0].1.as_deref(), Some(expected_rx));
}

// ============================================================================
// --remove-umi trims the UMI (and its delimiter) off the read name
// ============================================================================

#[test]
fn remove_umi_trims_read_name() {
    let dir = TempDir::new().unwrap();
    let input = write_input(dir.path(), &[record_named("inst:1:FC:1:1101:5:7:ACGT+TTGG")]);
    let output = dir.path().join("out.bam");
    run_copy_umi(&input, &output, &["--remove-umi"]).expect("copy-umi failed");

    let rows = read_name_and_rx(dir.path(), &output);
    assert_eq!(rows.len(), 1, "one input record must yield exactly one output record");
    assert_eq!(rows[0].0, "inst:1:FC:1:1101:5:7");
    assert_eq!(rows[0].1.as_deref(), Some("ACGT-TTGG"));
}

// ============================================================================
// Pre-existing RX: overwrite by default, error with --fail-if-tag-present
// ============================================================================

#[test]
fn overwrites_existing_rx_by_default() {
    let dir = TempDir::new().unwrap();
    let input = write_input(dir.path(), &[record_named_with_rx("blah:ACGT", "STALEVAL")]);
    let output = dir.path().join("out.bam");
    run_copy_umi(&input, &output, &[]).expect("copy-umi failed");

    let rows = read_name_and_rx(dir.path(), &output);
    assert_eq!(rows.len(), 1, "one input record must yield exactly one output record");
    assert_eq!(rows[0].1.as_deref(), Some("ACGT"), "stale RX should be replaced");
}

#[test]
fn fail_if_tag_present_errors_on_existing_rx() {
    let dir = TempDir::new().unwrap();
    let input = write_input(dir.path(), &[record_named_with_rx("blah:ACGT", "STALEVAL")]);
    let output = dir.path().join("out.bam");
    let err = run_copy_umi(&input, &output, &["--fail-if-tag-present"])
        .expect_err("should error when RX already present");
    assert!(
        format!("{err:#}").contains("already has an RX tag"),
        "error should name the pre-existing RX tag, got: {err:#}"
    );
}

// ============================================================================
// Fail-fast on malformed read names — each error must name its actual cause,
// so an unrelated failure cannot pass the test.
// ============================================================================

#[rstest]
#[case::single_field("NAME", "no ':'-delimited UMI field")]
#[case::empty_last_field("1:2:3:", "empty UMI field")]
#[case::illegal_char("NAME:CCKC", "illegal character")]
#[case::coordinate_not_a_umi("inst:1:FC:1:1101:5:7:10799", "illegal character")]
#[case::empty_after_normalization("blah:r", "normalizes to an empty UMI")]
fn errors_on_malformed_name(#[case] name: &str, #[case] expected_msg: &str) {
    let dir = TempDir::new().unwrap();
    let input = write_input(dir.path(), &[record_named(name)]);
    let output = dir.path().join("out.bam");
    let err =
        run_copy_umi(&input, &output, &[]).expect_err(&format!("expected failure for '{name}'"));
    assert!(
        format!("{err:#}").contains(expected_msg),
        "error for '{name}' should contain '{expected_msg}', got: {err:#}"
    );
}

#[test]
fn remove_umi_errors_when_it_would_empty_the_read_name() {
    // A leading-delimiter name (`:ACGT`) has its UMI as the whole tail; trimming
    // it would leave an empty QNAME, which must be rejected rather than written.
    let dir = TempDir::new().unwrap();
    let input = write_input(dir.path(), &[record_named(":ACGT")]);
    let output = dir.path().join("out.bam");
    let err = run_copy_umi(&input, &output, &["--remove-umi"])
        .expect_err("should refuse to empty the read name");
    assert!(
        format!("{err:#}").contains("would leave an empty read name"),
        "error should explain the empty-name refusal, got: {err:#}"
    );
}

// ============================================================================
// Idempotency: running twice (no trim) yields identical RX
// ============================================================================

#[test]
fn idempotent_without_trim() {
    let dir = TempDir::new().unwrap();
    let input = write_input(dir.path(), &[record_named("blah:ACGT+TTGG")]);
    let once = dir.path().join("once.bam");
    let twice = dir.path().join("twice.bam");
    run_copy_umi(&input, &once, &[]).expect("first pass");
    run_copy_umi(&once, &twice, &[]).expect("second pass");

    let a = read_name_and_rx(dir.path(), &once);
    let dir2 = TempDir::new().unwrap();
    let b = read_name_and_rx(dir2.path(), &twice);
    assert_eq!(a.len(), 1, "first pass must yield exactly one output record");
    assert_eq!(b.len(), 1, "second pass must yield exactly one output record");
    assert_eq!(a, b, "copy-umi should be idempotent when the name still carries the UMI");
}

// ============================================================================
// Metrics
// ============================================================================

/// Read the single metrics row as a `column -> value` map.
fn read_metrics_row(path: &Path) -> HashMap<String, String> {
    let text = std::fs::read_to_string(path).expect("read metrics");
    let mut lines = text.lines();
    let header: Vec<String> =
        lines.next().expect("header line").split('\t').map(str::to_string).collect();
    let values: Vec<String> =
        lines.next().expect("value line").split('\t').map(str::to_string).collect();
    assert_eq!(
        header.len(),
        values.len(),
        "metrics header and value rows must have equal column counts"
    );
    assert!(lines.next().is_none(), "metrics file must contain exactly one data row");
    header.into_iter().zip(values).collect()
}

#[test]
fn writes_metrics_row() {
    let dir = TempDir::new().unwrap();
    let input = write_input(
        dir.path(),
        &[
            record_named("blah:ACGT"),
            record_named_with_rx("blah:CCGG", "STALE"),
            record_named("blah:TTAA"),
        ],
    );
    let output = dir.path().join("out.bam");
    let metrics = dir.path().join("metrics.tsv");
    run_copy_umi(&input, &output, &["-M", &metrics.display().to_string()])
        .expect("copy-umi failed");

    let row = read_metrics_row(&metrics);
    assert_eq!(row["total_records"], "3");
    assert_eq!(row["rx_written"], "3");
    assert_eq!(row["rx_overwritten"], "1");
    assert_eq!(row["names_trimmed"], "0");
}

#[test]
fn metrics_count_trimmed_names_under_remove_umi() {
    let dir = TempDir::new().unwrap();
    let input = write_input(dir.path(), &[record_named("blah:ACGT"), record_named("blah:TTAA")]);
    let output = dir.path().join("out.bam");
    let metrics = dir.path().join("metrics.tsv");
    run_copy_umi(&input, &output, &["--remove-umi", "-M", &metrics.display().to_string()])
        .expect("copy-umi failed");

    let row = read_metrics_row(&metrics);
    assert_eq!(row["total_records"], "2");
    assert_eq!(row["names_trimmed"], "2", "every record's name was trimmed");
}

// ============================================================================
// Non-destructive contract: SEQ and unrelated tags survive; only RX changes
// ============================================================================

/// Every whitespace-split field of every non-header SAM line.
fn read_sam_fields(dir: &Path, bam: &Path) -> Vec<Vec<String>> {
    let sam = dir.join("fields.sam");
    transcode_bam_to_sam(bam, &sam);
    std::fs::read_to_string(&sam)
        .expect("read sam")
        .lines()
        .filter(|l| !l.starts_with('@'))
        .map(|line| line.split('\t').map(str::to_string).collect())
        .collect()
}

#[test]
fn preserves_sequence_and_unrelated_tags() {
    // A record with a stale RX, an unrelated BX tag, and a known SEQ/POS.
    let mut b = SamBuilder::new();
    b.read_name(b"blah:ACGT")
        .sequence(b"ACGTACGT")
        .qualities(&[30; 8])
        .flags(0)
        .ref_id(0)
        .pos(99)
        .mapq(60)
        .cigar_ops(&[8 << 4])
        .add_string_tag(fgumi_lib::sam::SamTag::RX, b"STALE")
        .add_string_tag(fgumi_lib::sam::SamTag::OX, b"KEEPME");
    let dir = TempDir::new().unwrap();
    let input = write_input(dir.path(), &[b.build()]);
    let output = dir.path().join("out.bam");
    run_copy_umi(&input, &output, &[]).expect("copy-umi failed");

    // Exact-preservation oracle: transcode both input and output and compare the
    // whole record. Only RX may change; everything else must be byte-for-byte equal.
    let in_rows = read_sam_fields(dir.path(), &input);
    let out_rows = read_sam_fields(dir.path(), &output);
    assert_eq!(in_rows.len(), 1, "one input record must yield exactly one output record");
    assert_eq!(out_rows.len(), 1, "one input record must yield exactly one output record");
    let (in_rec, out_rec) = (&in_rows[0], &out_rows[0]);

    // The 11 mandatory SAM fields (QNAME..QUAL) pass through untouched: QNAME,
    // FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL.
    assert_eq!(
        in_rec[..11],
        out_rec[..11],
        "only RX may change; the 11 mandatory fields must be identical"
    );

    // The optional-tag multiset must match exactly once the expected RX rewrite is
    // applied — catching any dropped, added, duplicated, or mutated tag, and a
    // second RX. Sort both sides so tag order is not asserted.
    let mut expected_tags: Vec<String> = in_rec[11..]
        .iter()
        .map(|t| if t.starts_with("RX:") { "RX:Z:ACGT".to_string() } else { t.clone() })
        .collect();
    expected_tags.sort();
    let mut got_tags: Vec<String> = out_rec[11..].to_vec();
    got_tags.sort();
    assert_eq!(
        got_tags, expected_tags,
        "only RX:Z:STALE -> RX:Z:ACGT; every other tag preserved with no add/drop/dup"
    );
}

// ============================================================================
// Write-target guards: a write path that aliases the input is rejected
// ============================================================================

#[test]
fn rejects_output_aliasing_input() {
    let dir = TempDir::new().unwrap();
    let input = write_input(dir.path(), &[record_named("blah:ACGT")]);
    // Snapshot the whole input BAM so we can prove the rejection had no write
    // side effect — an error message alone can be produced *after* a writer has
    // already truncated the aliased file.
    let before = std::fs::read(&input).expect("read input before");
    // -o == -i would clobber the BAM being read; the guard must reject it.
    let err = run_copy_umi(&input, &input, &[]).expect_err("should refuse to clobber the input");
    assert!(
        format!("{err:#}").contains("same file as --input"),
        "error should name the input-aliasing conflict, got: {err:#}"
    );
    // The guard fires before any writer opens, so the input must be byte-for-byte
    // intact — not truncated or partially rewritten.
    let after = std::fs::read(&input).expect("read input after");
    assert_eq!(after, before, "the input BAM must be untouched when the alias is rejected");
}

#[test]
fn rejects_metrics_aliasing_input() {
    // The metrics path is a second write target and gets the same clobber guard.
    let dir = TempDir::new().unwrap();
    let input = write_input(dir.path(), &[record_named("blah:ACGT")]);
    let output = dir.path().join("out.bam");
    let before = std::fs::read(&input).expect("read input before");
    let err = run_copy_umi(&input, &output, &["-M", &input.display().to_string()])
        .expect_err("should refuse to write metrics over the input");
    assert!(
        format!("{err:#}").contains("same file as --input"),
        "error should name the input-aliasing conflict, got: {err:#}"
    );
    // Rejection happens before any writer opens: the input is untouched and no
    // output BAM was created.
    let after = std::fs::read(&input).expect("read input after");
    assert_eq!(after, before, "the input BAM must be untouched when the metrics alias is rejected");
    assert!(!output.exists(), "no output BAM should be created when the alias is rejected");
}

#[test]
fn rejects_output_and_metrics_colliding() {
    // Two write targets resolving to one path is rejected before any writer opens.
    let dir = TempDir::new().unwrap();
    let input = write_input(dir.path(), &[record_named("blah:ACGT")]);
    let collide = dir.path().join("out.bam");
    let before = std::fs::read(&input).expect("read input before");
    let err = run_copy_umi(&input, &collide, &["-M", &collide.display().to_string()])
        .expect_err("should refuse two writers to one path");
    let msg = format!("{err:#}");
    // The guard must fire specifically on the two-writers-one-path collision
    // (naming both flags), not merely surface some unrelated parse/IO error.
    assert!(
        msg.contains("both write to") && msg.contains("--output") && msg.contains("--metrics"),
        "error should name the --output/--metrics collision, got: {msg}"
    );
    // Rejection happens before either writer opens, so nothing was created and
    // the input is untouched.
    assert!(!collide.exists(), "no output should be written when the collision is rejected");
    let after = std::fs::read(&input).expect("read input after");
    assert_eq!(after, before, "the input BAM must be untouched when the collision is rejected");
}

// ============================================================================
// --field-delimiter selects a non-default separator
// ============================================================================

#[test]
fn respects_non_default_field_delimiter() {
    // With `_` as the delimiter, the UMI is the last `_` field; the `:` in the
    // name is ordinary content, not a delimiter.
    let dir = TempDir::new().unwrap();
    let input = write_input(dir.path(), &[record_named("a:b_c_ACGT")]);
    let output = dir.path().join("out.bam");
    run_copy_umi(&input, &output, &["--field-delimiter", "_"]).expect("copy-umi failed");

    let rows = read_name_and_rx(dir.path(), &output);
    assert_eq!(rows.len(), 1, "one input record must yield exactly one output record");
    assert_eq!(rows[0].1.as_deref(), Some("ACGT"));
}

#[test]
fn rejects_non_ascii_field_delimiter() {
    // A single non-ASCII char parses as a `char` but is rejected at validation.
    let dir = TempDir::new().unwrap();
    let input = write_input(dir.path(), &[record_named("blah:ACGT")]);
    let output = dir.path().join("out.bam");
    let err = run_copy_umi(&input, &output, &["--field-delimiter", "€"])
        .expect_err("non-ASCII delimiter must be rejected");
    assert!(
        format!("{err:#}").contains("single ASCII character"),
        "error should explain the ASCII requirement, got: {err:#}"
    );
}

// ============================================================================
// The parallel pipeline path (--threads) produces the same result
// ============================================================================

#[test]
fn runs_on_multiple_threads() {
    let dir = TempDir::new().unwrap();
    let names = ["blah:AAAA", "blah:CCCC", "blah:GGGG", "blah:TTTT"];
    let input = write_input(dir.path(), &names.iter().map(|n| record_named(n)).collect::<Vec<_>>());
    let output = dir.path().join("out.bam");
    run_copy_umi(&input, &output, &["--threads", "4"]).expect("copy-umi failed");

    // Compare the complete (QNAME, RX) rows in input order, not just the RX
    // column: a threaded-path defect that reordered records or swapped read
    // names while still emitting the expected RX set would pass an RX-only check.
    let rows = read_name_and_rx(dir.path(), &output);
    let want: Vec<(String, Option<String>)> = names
        .iter()
        .map(|n| {
            // copy-umi copies the name's UMI into RX and leaves the name intact
            // (no --remove-umi); `blah:<UMI>` → RX `<UMI>`.
            let umi = n.strip_prefix("blah:").expect("test fixture name shape is blah:<UMI>");
            ((*n).to_string(), Some(umi.to_string()))
        })
        .collect();
    assert_eq!(
        rows, want,
        "the parallel path must preserve every record's QNAME and its copied RX, in input order"
    );
}

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
#[case::r_stripped_not_revcomped("blah:rAAAA", "blah:rAAAA", "AAAA")]
#[case::r_stripped_dual("blah:rAAAA+CCCC", "blah:rAAAA+CCCC", "AAAA-CCCC")]
fn strip_only_does_not_reverse_complement(
    #[case] name: &str,
    #[case] expected_name: &str,
    #[case] expected_rx: &str,
) {
    let dir = TempDir::new().unwrap();
    let input = write_input(dir.path(), &[record_named(name)]);
    let output = dir.path().join("out.bam");
    run_copy_umi(&input, &output, &["--reverse-complement-r-umis", "false"])
        .expect("copy-umi failed");

    let rows = read_name_and_rx(dir.path(), &output);
    assert_eq!(rows.len(), 1, "one input record must yield exactly one output record");
    assert_eq!(rows[0].0, expected_name, "read name should be unchanged");
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

/// `--fail-if-tag-present` must abort on both the serial oracle and the
/// `--threads` chain path — proving the flag is wired on the parallel path too,
/// not just the serial one.
#[rstest]
#[case::serial(&[][..])]
#[case::chain(&["--threads", "2"][..])]
fn fail_if_tag_present_errors_on_existing_rx(#[case] thread_flags: &[&str]) {
    let dir = TempDir::new().unwrap();
    let input = write_input(dir.path(), &[record_named_with_rx("blah:ACGT", "STALEVAL")]);
    let output = dir.path().join("out.bam");
    let mut flags = vec!["--fail-if-tag-present"];
    flags.extend_from_slice(thread_flags);
    let err =
        run_copy_umi(&input, &output, &flags).expect_err("should error when RX already present");
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
    assert_eq!(rows[0].0, "a:b_c_ACGT", "read name should be unchanged");
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

// ============================================================================
// Chain-vs-oracle parity (the R3 cutover contract)
//
// The no-`--threads` serial path (`run_single_threaded`) is the in-process
// oracle; `--threads N` routes through the chain builder. `--threads 1` is a
// distinct chain-at-single-worker path from both, so the matrix includes it.
// ============================================================================

/// Build a copy-umi input of `n` records with `:`-delimited UMI read names.
/// `n=500` crosses pipeline batch boundaries, exercising the cross-slot metrics
/// reduction and the multi-batch ordering.
fn build_parity_input(dir: &Path, n: usize) -> std::path::PathBuf {
    let bases = [b'A', b'C', b'G', b'T'];
    let records: Vec<RawRecord> = (0..n)
        .map(|i| {
            let umi: String = (0..8).map(|j| bases[(i + j) % 4] as char).collect();
            record_named(&format!("inst:1:FC:1:1101:{i}:{i}:{umi}"))
        })
        .collect();
    write_input(dir, &records)
}

#[rstest]
#[case::t1(1)]
#[case::t2(2)]
#[case::t4(4)]
fn chain_matches_single_threaded(#[case] threads: u8) {
    let dir = TempDir::new().unwrap();
    let input = build_parity_input(dir.path(), 500);

    let oracle_out = dir.path().join("oracle.bam");
    run_copy_umi(&input, &oracle_out, &[]).expect("serial oracle run");

    let chain_out = dir.path().join("chain.bam");
    run_copy_umi(&input, &chain_out, &["--threads", &threads.to_string()]).expect("chain run");

    let (oracle_hdr, oracle_recs) = crate::helpers::read_bam_output(&oracle_out);
    let (chain_hdr, chain_recs) = crate::helpers::read_bam_output(&chain_out);
    assert_eq!(chain_hdr, oracle_hdr, "normalized header parity (threads={threads})");
    assert_eq!(chain_recs, oracle_recs, "record parity (threads={threads})");
    // Non-vacuous: every record survived and carries the copied RX.
    assert_eq!(chain_recs.len(), 500, "all records emitted");
}

/// Every non-default, output-producing knob must be wired identically on the
/// chain (`--threads`) path and the serial oracle. The default-flag parity test
/// above exercises none of them, so a knob dropped or mis-wired on the parallel
/// path would go unnoticed. The input carries `r`-prefixed and dual UMIs so both
/// `--remove-umi` (name trim) and `--reverse-complement-r-umis false` (strip vs
/// reverse-complement) change the output — if the chain dropped a flag, its
/// records would diverge from the oracle's.
/// The `expected` case column is the absolute `(read_name, RX)` the oracle must
/// produce for each of the three input records under the case's flags. Pinning
/// these fixed values — not just chain-vs-oracle parity — is what catches a bug
/// present on *both* paths (e.g. a flag mis-parsed in shared code): parity alone
/// would still pass. Verified against the command's observed output; note
/// `--reverse-complement-r-umis false` still maps the `+` dual-UMI delimiter to
/// `-` (it only disables the `r`-prefix reverse-complement), so the strip case
/// yields `AAAA-CCCC`, not `AAAA+CCCC`.
#[rstest]
#[case::remove_umi(
    &["--remove-umi"][..],
    &[
        ("inst:1:FC:1:1101:5:7", "TTTT-CCCC"),
        ("inst:1:FC:1:1101:5:8", "CCCC"),
        ("inst:1:FC:1:1101:5:9", "TTTT"),
    ],
)]
#[case::strip_only(
    &["--reverse-complement-r-umis", "false"][..],
    &[
        ("inst:1:FC:1:1101:5:7:rAAAA+CCCC", "AAAA-CCCC"),
        ("inst:1:FC:1:1101:5:8:rGGGG", "GGGG"),
        ("inst:1:FC:1:1101:5:9:TTTT", "TTTT"),
    ],
)]
#[case::remove_and_strip(
    &["--remove-umi", "--reverse-complement-r-umis", "false"][..],
    &[
        ("inst:1:FC:1:1101:5:7", "AAAA-CCCC"),
        ("inst:1:FC:1:1101:5:8", "GGGG"),
        ("inst:1:FC:1:1101:5:9", "TTTT"),
    ],
)]
fn chain_matches_oracle_with_nondefault_options(
    #[case] option_flags: &[&str],
    #[case] expected: &[(&str, &str)],
) {
    let dir = TempDir::new().unwrap();
    let records = [
        record_named("inst:1:FC:1:1101:5:7:rAAAA+CCCC"),
        record_named("inst:1:FC:1:1101:5:8:rGGGG"),
        record_named("inst:1:FC:1:1101:5:9:TTTT"),
    ];
    let input = write_input(dir.path(), &records);

    let oracle_out = dir.path().join("oracle.bam");
    run_copy_umi(&input, &oracle_out, option_flags).expect("serial oracle run");

    let mut chain_flags = option_flags.to_vec();
    chain_flags.extend_from_slice(&["--threads", "4"]);
    let chain_out = dir.path().join("chain.bam");
    run_copy_umi(&input, &chain_out, &chain_flags).expect("chain run");

    // Absolute anchor: the oracle must produce exactly these names+RX. A flag
    // dropped in shared code would change these, even though chain==oracle held.
    let oracle_name_rx = read_name_and_rx(dir.path(), &oracle_out);
    let expected_name_rx: Vec<(String, Option<String>)> =
        expected.iter().map(|(name, rx)| ((*name).to_string(), Some((*rx).to_string()))).collect();
    assert_eq!(
        oracle_name_rx, expected_name_rx,
        "oracle (name, RX) must match the pinned expected values (flags={option_flags:?})",
    );

    let (oracle_hdr, oracle_recs) = crate::helpers::read_bam_output(&oracle_out);
    let (chain_hdr, chain_recs) = crate::helpers::read_bam_output(&chain_out);
    assert_eq!(chain_hdr, oracle_hdr, "normalized header parity (flags={option_flags:?})");
    assert_eq!(chain_recs, oracle_recs, "record parity (flags={option_flags:?})");
    assert_eq!(chain_recs.len(), 3, "all records emitted (flags={option_flags:?})");
}

/// `--metrics` TSV must be byte-identical between the serial oracle and the
/// chain. `fgumi_metrics::write_metrics` embeds no timestamp, so a raw byte
/// compare is stable.
#[test]
fn chain_metrics_match_single_threaded() {
    let dir = TempDir::new().unwrap();
    // Include a record with a pre-existing RX so `rx_overwritten` is non-zero.
    let mut records: Vec<RawRecord> =
        (0..50).map(|i| record_named(&format!("q:{i}:ACGT"))).collect();
    records.push(record_named_with_rx("q:99:TTTT", "AAAA"));
    let input = write_input(dir.path(), &records);

    let oracle_out = dir.path().join("oracle.bam");
    let oracle_metrics = dir.path().join("oracle.tsv");
    run_copy_umi(&input, &oracle_out, &["-M", &oracle_metrics.display().to_string()])
        .expect("serial run");

    let chain_out = dir.path().join("chain.bam");
    let chain_metrics = dir.path().join("chain.tsv");
    run_copy_umi(
        &input,
        &chain_out,
        &["--threads", "4", "-M", &chain_metrics.display().to_string()],
    )
    .expect("chain run");

    let oracle_tsv = std::fs::read(&oracle_metrics).expect("read oracle metrics");
    let chain_tsv = std::fs::read(&chain_metrics).expect("read chain metrics");
    assert_eq!(chain_tsv, oracle_tsv, "the --metrics TSV must be byte-identical across paths");

    // Absolute anchor: pin the metrics row to fixed expected counts, not just
    // chain==oracle parity. 50 fresh names + 1 with a pre-existing RX = 51
    // records; every record gets an RX (rx_written 51), the one carrying a stale
    // RX is overwritten (rx_overwritten 1), and no `--remove-umi` means no name
    // trims (names_trimmed 0). Shared miscounting would slip past a byte compare.
    let row = read_metrics_row(&oracle_metrics);
    assert_eq!(row.get("total_records").map(String::as_str), Some("51"), "row: {row:?}");
    assert_eq!(row.get("rx_written").map(String::as_str), Some("51"), "row: {row:?}");
    assert_eq!(row.get("rx_overwritten").map(String::as_str), Some("1"), "row: {row:?}");
    assert_eq!(row.get("names_trimmed").map(String::as_str), Some("0"), "row: {row:?}");
}

/// Fail-fast under `--threads`: one record with an empty last field aborts the
/// run, naming the read name. A single bad record keeps the surfaced
/// error deterministic under the parallel pipeline.
#[test]
fn fail_fast_under_threads_single_bad_record() {
    let dir = TempDir::new().unwrap();
    let bad = "read-with-empty-umi:";
    let input = write_input(dir.path(), &[record_named("q:1:ACGT"), record_named(bad)]);
    let output = dir.path().join("out.bam");
    let err = run_copy_umi(&input, &output, &["--threads", "2"])
        .expect_err("empty UMI field must abort the --threads run");
    let msg = format!("{err:#}");
    assert!(msg.contains(bad), "chain error names the read name; got: {msg}");
}

/// Two-level fail-fast: an illegal-character UMI is a
/// `with_context`-wrapped error. Both paths must abort naming the read name and
/// the outer "extracting UMI from read name" context. The serial oracle also
/// preserves the inner reason (full `anyhow` chain); the chain path's
/// `reconstruct` flattens to Display-only, so the inner reason is not asserted on
/// the chain path — this documents the framework flattening rather than asserting
/// a parity that does not hold.
#[test]
fn two_level_fail_fast_names_read_name_on_both_paths() {
    let dir = TempDir::new().unwrap();
    let bad = "q:1:ZZZZ"; // Z is outside ACGTN-, so normalization fails.
    let input = write_input(dir.path(), &[record_named(bad)]);

    let serial_out = dir.path().join("serial.bam");
    let serial_err = run_copy_umi(&input, &serial_out, &[])
        .expect_err("illegal UMI char must abort the serial run");
    let serial_msg = format!("{serial_err:#}");
    assert!(serial_msg.contains(bad), "serial names the read name; got: {serial_msg}");
    assert!(
        serial_msg.contains("extracting UMI from read name"),
        "serial carries the outer context; got: {serial_msg}"
    );

    let chain_out = dir.path().join("chain.bam");
    let chain_err = run_copy_umi(&input, &chain_out, &["--threads", "2"])
        .expect_err("illegal UMI char must abort the chain run");
    let chain_msg = format!("{chain_err:#}");
    assert!(chain_msg.contains(bad), "chain names the read name; got: {chain_msg}");
    assert!(
        chain_msg.contains("extracting UMI from read name"),
        "chain carries the outer context; got: {chain_msg}"
    );
}

/// Header-less input: both the chain and the serial oracle synthesize
/// `@HD VN:1.6 SO:unsorted` and stay in header parity (`ChainBuilder::new`
/// calls `ensure_hd_record`, so the chain path
/// synthesizes @HD just like the oracle).
#[test]
fn chain_and_oracle_synthesize_hd_for_headerless_input() {
    use noodles::sam::Header as SamHeader;
    use noodles::sam::header::record::value::{Map, map::ReferenceSequence};
    use std::num::NonZeroUsize;

    let dir = TempDir::new().unwrap();
    // A header with a reference sequence but NO @HD line.
    let header = SamHeader::builder()
        .add_reference_sequence(
            bstr::BString::from("chr1"),
            Map::<ReferenceSequence>::new(NonZeroUsize::new(10_000).unwrap()),
        )
        .build();
    assert!(header.header().is_none(), "fixture must be header-less");
    let input = dir.path().join("headerless.bam");
    write_bam(&input, &header, &[record_named("q:1:ACGT")]);

    let oracle_out = dir.path().join("oracle.bam");
    run_copy_umi(&input, &oracle_out, &[]).expect("serial run");
    let chain_out = dir.path().join("chain.bam");
    run_copy_umi(&input, &chain_out, &["--threads", "2"]).expect("chain run");

    let (oracle_hdr, _) = crate::helpers::read_bam_output(&oracle_out);
    let (chain_hdr, _) = crate::helpers::read_bam_output(&chain_out);
    assert!(oracle_hdr.header().is_some(), "serial synthesizes @HD");
    assert!(chain_hdr.header().is_some(), "chain synthesizes @HD");
    assert_eq!(chain_hdr, oracle_hdr, "header-less-input @HD parity");

    // Absolute anchor: pin the exact synthesized @HD contract on BOTH paths, not
    // just their equality — otherwise both could synthesize the same *wrong* @HD
    // and pass. `ChainBuilder::new`/the oracle synthesize `@HD VN:1.6 SO:unsorted`.
    let hd_line = |bam: &Path, tag: &str| -> String {
        let sam = dir.path().join(format!("{tag}.sam"));
        transcode_bam_to_sam(bam, &sam);
        std::fs::read_to_string(&sam)
            .expect("read sam")
            .lines()
            .find(|l| l.starts_with("@HD"))
            .unwrap_or_else(|| panic!("{tag} output must carry an @HD line"))
            .to_string()
    };
    for (tag, bam) in [("oracle", &oracle_out), ("chain", &chain_out)] {
        let hd = hd_line(bam, tag);
        assert!(hd.contains("VN:1.6"), "{tag} @HD must declare VN:1.6, got: {hd}");
        assert!(hd.contains("SO:unsorted"), "{tag} @HD must declare SO:unsorted, got: {hd}");
    }
}

// ============================================================================
// CRC-policy parity: the serial oracle (no `--threads`) and the
// chain (`--threads N`) derive `verify_crc` from independent sources — the
// serial path from `self.io.pipeline_reader_opts()`, the chain from
// `spec.verify_crc` — so pin them equal: both must honor `--check-crc` /
// `--no-check-crc` identically (default verify-on for file input rejects a
// corrupted block; `--no-check-crc` accepts it).
// ============================================================================

/// Flip a byte in the last BGZF block's CRC32 footer, so decoding that block
/// fails only when CRC verification is on. The input must span >= 2 BGZF blocks
/// so the corrupted block is not the header's block (which always verifies during
/// the header parse). Mirrors the sibling command tests (clip/group/dedup/…).
fn corrupt_last_block_crc(path: &Path) {
    use fgumi_lib::bgzf_reader::{BGZF_FOOTER_SIZE, RawBgzfBlock, read_raw_blocks};
    let mut bytes = std::fs::read(path).expect("read bam for corruption");
    let mut cursor: &[u8] = &bytes;
    let blocks = read_raw_blocks(&mut cursor, 100_000).expect("read bgzf blocks from test bam");
    assert!(
        blocks.len() >= 2,
        "test input must span >= 2 BGZF blocks so the corrupted block isn't the header's; got {}",
        blocks.len()
    );
    let offset: usize = blocks[..blocks.len() - 1].iter().map(RawBgzfBlock::len).sum();
    let last = blocks.last().expect("checked len >= 2 above");
    let crc_off = offset + last.len() - BGZF_FOOTER_SIZE;
    bytes[crc_off] ^= 0x01;
    std::fs::write(path, bytes).expect("write corrupted bam");
}

/// Both the serial oracle and the `--threads` chain must honor the CRC policy
/// identically: the default (verify-on for file input) rejects a corrupted BGZF
/// block, while `--no-check-crc` accepts it and completes.
#[rstest]
#[case::serial(&[][..])]
#[case::chain(&["--threads", "4"][..])]
fn crc_policy_parity_serial_and_chain(#[case] path_flags: &[&str]) {
    let dir = TempDir::new().unwrap();
    // 6000 records span several BGZF blocks, so the corrupted last block is a
    // record block, not the header's.
    let input = build_parity_input(dir.path(), 6000);
    corrupt_last_block_crc(&input);

    // Default policy: both paths reject the corrupted block.
    let out_reject = dir.path().join("reject.bam");
    assert!(
        run_copy_umi(&input, &out_reject, path_flags).is_err(),
        "default CRC policy must reject a corrupted block (flags={path_flags:?})"
    );

    // `--no-check-crc`: both paths accept the corrupted block and complete.
    let out_accept = dir.path().join("accept.bam");
    let mut accept_flags = path_flags.to_vec();
    accept_flags.push("--no-check-crc");
    run_copy_umi(&input, &out_accept, &accept_flags).unwrap_or_else(|e| {
        panic!("--no-check-crc must accept a corrupted block (flags={path_flags:?}): {e}")
    });
    // Assert the accepted output is record-for-record identical to an
    // *uncorrupted* run of the same input — a bare count (6000) would pass even
    // if `--no-check-crc` silently dropped, duplicated, or reordered the
    // corrupted block's records. Corrupting the CRC leaves the payload intact, so
    // the two outputs must match exactly (order + content).
    let (_hdr, recs) = crate::helpers::read_bam_output(&out_accept);
    assert_eq!(
        recs.len(),
        6000,
        "--no-check-crc must accept AND retain all records (flags={path_flags:?})"
    );
    let clean_dir = TempDir::new().unwrap();
    let clean_input = build_parity_input(clean_dir.path(), 6000);
    let clean_out = clean_dir.path().join("clean.bam");
    run_copy_umi(&clean_input, &clean_out, path_flags).expect("uncorrupted reference run");
    let (_clean_hdr, clean_recs) = crate::helpers::read_bam_output(&clean_out);
    assert_eq!(
        recs, clean_recs,
        "--no-check-crc output must be record-for-record identical to the uncorrupted run \
         (flags={path_flags:?})"
    );
}

// ============================================================================
// Malformed-record safety: a framing-consistent-but-malformed record (an
// `l_read_name` that runs past the record end) must fail with a CLEAN error on
// BOTH the serial oracle and the `--threads` chain — never panic. Before the
// serial path revalidated framing, `copy_umi_into_record`'s `read_name()` slice
// would panic (index out of bounds) on such a record, a regression from the
// pre-cutover pipeline, which decoded (and thus validated) every record.
// ============================================================================

/// Write a one-record BAM whose record has a corrupted `l_read_name` (claims a
/// 200-byte read name the short record cannot hold), via the RAW writer so the
/// corruption survives (`write_bam` would re-normalize it through noodles).
fn write_malformed_l_read_name_bam(dir: &Path) -> std::path::PathBuf {
    let path = dir.join("malformed.bam");
    let header = create_minimal_header("chr1", 10_000);
    let mut record = record_named("inst:1:FC:1:1101:5:7:ACGT");
    // Byte 8 of the record body is `l_read_name`; 200 runs far past this short
    // record's end, so an unguarded `read_name()` would slice out of bounds.
    record.as_mut_vec()[8] = 200;
    let mut writer =
        fgumi_bam_io::create_raw_bam_writer(&path, &header, 1, 6).expect("create raw BAM writer");
    writer.write_raw_record(record.as_ref()).expect("write malformed record");
    writer.finish().expect("finish malformed BAM");
    path
}

#[rstest]
#[case::serial(&[][..])]
#[case::chain(&["--threads", "2"][..])]
fn malformed_record_errors_cleanly_on_both_paths(#[case] thread_flags: &[&str]) {
    let dir = TempDir::new().unwrap();
    let input = write_malformed_l_read_name_bam(dir.path());
    let output = dir.path().join("out.bam");
    let err = run_copy_umi(&input, &output, thread_flags)
        .expect_err("a malformed l_read_name record must error, not panic");
    assert!(
        format!("{err:#}").contains("read-name region runs past record end"),
        "expected the clean framing error on both paths (flags={thread_flags:?}), got: {err:#}"
    );
}

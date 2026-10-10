//! End-to-end tests for the merge command's input sort-order guard (MERGE3-01).
//!
//! Two mechanisms are exercised through `Merge::execute`:
//! 1. the fast header-declared check (rejects a declared-order conflict before
//!    reading records, e.g. a coordinate BAM into the template-coordinate default), and
//! 2. the streaming monotonicity verify (catches actual disorder in a
//!    bare/undeclared-header input, and leaves no partial output).

use fgumi_lib::commands::command::Command;
use fgumi_lib::commands::merge::Merge;
use fgumi_lib::commands::sort::SortOrderArg;
use fgumi_lib::sam::SamTag;
use fgumi_raw_bam::{RawRecord, SamBuilder, flags};
use noodles::bam;
use noodles::sam::Header;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use std::fs;
use std::path::{Path, PathBuf};
use tempfile::TempDir;

use crate::helpers::assertions::assert_bam_sorted;
use crate::helpers::bam_generator::{
    create_coordinate_sorted_header, create_minimal_header, to_record_buf,
};

/// Path to the built `fgumi` binary, for the end-to-end merge-output test.
const FGUMI: &str = env!("CARGO_BIN_EXE_fgumi");

/// A mapped 20M record on chr1 (ref 0) at the given 1-based position.
fn mapped_record(name: &[u8], pos: i32) -> RawRecord {
    let mut b = SamBuilder::new();
    b.read_name(name)
        .ref_id(0)
        .pos(pos - 1) // 1-based → 0-based
        .mapq(60)
        .cigar_ops(&[20 << 4]) // 20M
        .sequence(b"ACGTACGTACGTACGTACGT")
        .qualities(&[30; 20]);
    b.build()
}

/// Write `records` to `path` under `header`.
fn write_bam(path: &Path, header: &Header, records: &[RawRecord]) {
    let mut writer = bam::io::Writer::new(fs::File::create(path).expect("create bam"));
    writer.write_header(header).expect("write header");
    for r in records {
        writer.write_alignment_record(header, &to_record_buf(r)).expect("write record");
    }
    writer.try_finish().expect("finish bam");
}

/// Assert no `.fgumi-merge-*.tmp` staging sibling remains in `dir` (the atomic
/// temp+persist must clean up its temp on both the success and rejection paths).
fn assert_no_leftover_temp(dir: &Path) {
    let leftover_temps: Vec<_> = fs::read_dir(dir)
        .unwrap()
        .filter_map(Result::ok)
        .filter(|e| e.file_name().to_string_lossy().contains("fgumi-merge-"))
        .map(|e| e.file_name())
        .collect();
    assert!(leftover_temps.is_empty(), "merge temp files must be cleaned up: {leftover_temps:?}");
}

/// Read a merged BAM into `(read name, 1-based alignment start)` pairs in file
/// order — an identity-aware oracle for asserting merge interleaving.
fn read_name_pos(path: &Path) -> Vec<(String, usize)> {
    let mut reader = bam::io::Reader::new(fs::File::open(path).unwrap());
    let _ = reader.read_header().unwrap();
    reader
        .records()
        .map(|r| {
            let rec = r.unwrap();
            let name = String::from_utf8(rec.name().unwrap().to_vec()).unwrap();
            let pos = rec.alignment_start().unwrap().unwrap().get();
            (name, pos)
        })
        .collect()
}

/// A header with a `chr1` reference and **no** `@HD` sort-order tag — its declared
/// order classifies as "none", so it passes the header check and is verified by
/// the streaming monotonicity check during the merge.
fn bare_header() -> Header {
    use noodles::sam::header::record::value::Map;
    use noodles::sam::header::record::value::map::ReferenceSequence;
    use std::num::NonZeroUsize;
    Header::builder()
        .add_reference_sequence(
            "chr1",
            Map::<ReferenceSequence>::new(NonZeroUsize::new(10000).expect("nonzero")),
        )
        .build()
}

fn merge(
    output: &Path,
    inputs: Vec<std::path::PathBuf>,
    order: SortOrderArg,
) -> anyhow::Result<()> {
    merge_with_command_line(output, inputs, order, "fgumi merge")
}

/// [`merge`], recording `command_line` as the merge's command line.
fn merge_with_command_line(
    output: &Path,
    inputs: Vec<std::path::PathBuf>,
    order: SortOrderArg,
    command_line: &str,
) -> anyhow::Result<()> {
    Merge {
        output: output.to_path_buf(),
        inputs,
        input_list: None,
        order,
        threads: 1,
        compression_level: 1,
    }
    .execute(command_line)
}

/// MERGE3-01 header check: a coordinate-sorted input fed to the default
/// `--order template-coordinate` merge is rejected before any records are read.
#[test]
fn test_merge_rejects_coordinate_input_into_template_coordinate_default() {
    let dir = TempDir::new().unwrap();
    let header = create_coordinate_sorted_header("chr1", 10000);
    let a = dir.path().join("a.bam");
    let b = dir.path().join("b.bam");
    write_bam(&a, &header, &[mapped_record(b"reada0", 100), mapped_record(b"reada1", 200)]);
    write_bam(&b, &header, &[mapped_record(b"readb0", 150), mapped_record(b"readb1", 250)]);
    let out = dir.path().join("merged.bam");

    // Default order is template-coordinate; the coordinate inputs conflict.
    let err = merge(&out, vec![a, b], SortOrderArg::TemplateCoordinate)
        .expect_err("must reject coordinate input into template-coordinate merge");
    let msg = err.to_string();
    assert!(msg.contains("sorted by coordinate"), "unexpected error: {msg}");
    assert!(msg.contains("a.bam"), "error should name the offending input: {msg}");
    assert!(!out.exists(), "no output should be written on a declared-order conflict");
}

/// The CLI default for `--order` is template-coordinate — the fgumi-specific
/// footgun this guard defends (`fgumi merge a.bam b.bam` on plain
/// coordinate-sorted BAMs would otherwise silently corrupt). The reject test
/// above passes `TemplateCoordinate` through the `merge()` helper, bypassing
/// clap; parse without `--order` here so the default the guard relies on is
/// itself pinned.
#[test]
fn test_merge_order_defaults_to_template_coordinate() {
    use clap::Parser;
    let parsed = Merge::try_parse_from(["fgumi-merge", "-o", "out.bam", "in.bam"])
        .expect("merge parses without --order");
    assert_eq!(parsed.order, SortOrderArg::TemplateCoordinate);
}

/// A valid merge (coordinate inputs, `--order coordinate`) succeeds and produces
/// globally coordinate-sorted output.
#[test]
fn test_merge_valid_coordinate_succeeds() {
    let dir = TempDir::new().unwrap();
    let header = create_coordinate_sorted_header("chr1", 10000);
    let a = dir.path().join("a.bam");
    let b = dir.path().join("b.bam");
    write_bam(&a, &header, &[mapped_record(b"reada0", 100), mapped_record(b"reada1", 300)]);
    write_bam(&b, &header, &[mapped_record(b"readb0", 200), mapped_record(b"readb1", 400)]);
    let out = dir.path().join("merged.bam");

    merge(&out, vec![a, b], SortOrderArg::Coordinate).expect("valid coordinate merge");
    assert!(out.exists());

    // The atomic temp+rename must leave no `.fgumi-merge-*.tmp` sibling on the
    // success path either — same cleanup contract the failure test checks.
    assert_no_leftover_temp(dir.path());

    // Assert record *identity* at each position, not just that positions are
    // sorted — a name-aware check fails on interleaving bugs that leave the
    // coordinate column sorted but swap records between inputs.
    let got = read_name_pos(&out);
    assert_eq!(
        got,
        vec![
            ("reada0".to_string(), 100),
            ("readb0".to_string(), 200),
            ("reada1".to_string(), 300),
            ("readb1".to_string(), 400),
        ],
        "merged output must interleave the two inputs in coordinate order"
    );
}

/// MERGE3-01 streaming verify: a bare-header input whose records are actually
/// out of coordinate order passes the header check but is caught record-by-record
/// during the merge — and no partial/temp output is left behind (temp + rename).
#[test]
fn test_merge_streaming_verify_rejects_missorted_bare_header_no_partial_output() {
    let dir = TempDir::new().unwrap();
    let coord = create_coordinate_sorted_header("chr1", 10000);
    let good = dir.path().join("good.bam");
    write_bam(&good, &coord, &[mapped_record(b"readg0", 100)]);

    // Bare header (declared order "none"), records descending in position.
    let bad = dir.path().join("bad.bam");
    write_bam(&bad, &bare_header(), &[mapped_record(b"readx", 300), mapped_record(b"ready", 100)]);

    let out = dir.path().join("merged.bam");
    let err = merge(&out, vec![good, bad.clone()], SortOrderArg::Coordinate)
        .expect_err("must reject an actually-mis-sorted input");
    let msg = err.to_string();
    assert!(msg.contains("not sorted in coordinate order"), "unexpected error: {msg}");
    assert!(msg.contains("bad.bam"), "error should name the offending input: {msg}");

    // No partial output and no leftover temp file. The merge temp is a uniquely
    // named `.fgumi-merge-*.tmp` sibling, so scan the directory for that prefix
    // rather than a fixed name (a RAII NamedTempFile removes it on any error).
    assert!(!out.exists(), "a rejected merge must leave no output file");
    assert_no_leftover_temp(dir.path());
}

/// MERGE3-01 streaming verify, destination preservation: a streaming-order
/// failure must not touch a *pre-existing* destination. The atomic temp+persist
/// writes to a sibling temp and only renames on success, so a rejected merge must
/// leave an already-present output byte-for-byte unchanged (the sibling
/// absent-output test above only proves an absent output stays absent).
#[test]
fn test_merge_streaming_verify_preserves_existing_destination() {
    let dir = TempDir::new().unwrap();
    let coord = create_coordinate_sorted_header("chr1", 10000);
    let good = dir.path().join("good.bam");
    write_bam(&good, &coord, &[mapped_record(b"readg0", 100)]);

    // Bare header (declared order "none"), records descending in position — caught
    // by the streaming monotonicity check, not the fast header-declared check.
    let bad = dir.path().join("bad.bam");
    write_bam(&bad, &bare_header(), &[mapped_record(b"readx", 300), mapped_record(b"ready", 100)]);

    // Pre-create the destination with sentinel bytes that are *not* a valid BAM,
    // so any accidental write (partial output, clobber) is detectable.
    let out = dir.path().join("merged.bam");
    let sentinel: &[u8] = b"pre-existing sentinel output that must survive a rejected merge";
    fs::write(&out, sentinel).unwrap();

    let err = merge(&out, vec![good, bad.clone()], SortOrderArg::Coordinate)
        .expect_err("must reject an actually-mis-sorted input");
    assert!(err.to_string().contains("not sorted in coordinate order"), "unexpected error: {err}");

    // The destination must still exist with exactly its original bytes...
    assert!(out.exists(), "a rejected merge must not remove a pre-existing destination");
    assert_eq!(
        fs::read(&out).unwrap(),
        sentinel,
        "a rejected merge must leave a pre-existing destination byte-for-byte unchanged"
    );
    // ...and no leftover temp sibling from the discarded staged output.
    assert_no_leftover_temp(dir.path());
}

/// MERGE3-01 streaming verify boundary: two consecutive records at the *same*
/// coordinate are a tie, not a decrease, so the monotonicity check (`<`, not
/// `<=`) must accept them. Guards the `>` vs `>=` off-by-one in the comparator
/// that the strictly-decreasing rejection test above cannot surface.
#[test]
fn test_merge_streaming_verify_accepts_equal_adjacent_positions() {
    let dir = TempDir::new().unwrap();
    // Bare header (declared order "none") routes records through the streaming
    // monotonicity check rather than the fast header-declared check.
    let input = dir.path().join("ties.bam");
    write_bam(
        &input,
        &bare_header(),
        &[mapped_record(b"read0", 100), mapped_record(b"read1", 100)],
    );

    let out = dir.path().join("merged.bam");
    merge(&out, vec![input], SortOrderArg::Coordinate)
        .expect("equal adjacent positions must be accepted, not rejected as mis-sorted");
    assert!(out.exists());

    // The tied records must pass through in their original order, unaltered.
    let got = read_name_pos(&out);
    assert_eq!(
        got,
        vec![("read0".to_string(), 100), ("read1".to_string(), 100)],
        "tied-position records must be accepted and preserved in order"
    );
}

/// The atomic temp+persist output must carry normal file permissions, not the
/// `0600` a `NamedTempFile` is created with. Overwriting an existing destination
/// keeps that destination's mode (matching the pre-atomic-temp `File::create`
/// path), which is deterministic regardless of the test process's umask.
#[cfg(unix)]
#[test]
fn test_merge_output_preserves_existing_destination_mode() {
    use std::os::unix::fs::PermissionsExt;

    let dir = TempDir::new().unwrap();
    let header = create_coordinate_sorted_header("chr1", 10000);
    let a = dir.path().join("a.bam");
    let b = dir.path().join("b.bam");
    write_bam(&a, &header, &[mapped_record(b"reada0", 100)]);
    write_bam(&b, &header, &[mapped_record(b"readb0", 200)]);

    // Pre-create the destination with a group/other-readable mode (0644); the
    // temp NamedTempFile is 0600, so without the mode fix the merge would leave
    // the output owner-only.
    let out = dir.path().join("merged.bam");
    fs::write(&out, b"placeholder").unwrap();
    fs::set_permissions(&out, fs::Permissions::from_mode(0o644)).unwrap();

    merge(&out, vec![a, b], SortOrderArg::Coordinate).expect("valid coordinate merge");

    let mode = fs::metadata(&out).unwrap().permissions().mode() & 0o777;
    assert_eq!(mode, 0o644, "merge must preserve the destination's existing mode, got {mode:#o}");
}

/// A newly created merge output (destination absent) must carry the mode a plain
/// `File::create` produces (`0o666 & !umask`), not the staging `NamedTempFile`'s
/// private `0600`. Compare against a control created via `File::create` in the
/// same process so the expectation tracks the test's umask deterministically —
/// an independent oracle for the new-file branch of `target_file_mode`.
#[cfg(unix)]
#[test]
fn test_merge_new_output_matches_file_create_mode() {
    use std::os::unix::fs::PermissionsExt;

    let dir = TempDir::new().unwrap();
    let header = create_coordinate_sorted_header("chr1", 10000);
    let a = dir.path().join("a.bam");
    let b = dir.path().join("b.bam");
    write_bam(&a, &header, &[mapped_record(b"reada0", 100)]);
    write_bam(&b, &header, &[mapped_record(b"readb0", 200)]);

    // Oracle: the mode File::create yields under this process's umask.
    let control = dir.path().join("control");
    fs::File::create(&control).unwrap();
    let expected = fs::metadata(&control).unwrap().permissions().mode() & 0o777;

    // Merge into a previously absent destination.
    let out = dir.path().join("merged.bam");
    assert!(!out.exists(), "destination must not exist before the merge");
    merge(&out, vec![a, b], SortOrderArg::Coordinate).expect("valid coordinate merge");

    let mode = fs::metadata(&out).unwrap().permissions().mode() & 0o777;
    assert_eq!(
        mode, expected,
        "new merge output must match File::create's mode ({expected:#o}), got {mode:#o}"
    );
}

/// The atomic temp+persist replaces its destination with a regular file, so a
/// non-regular output (FIFO, device, socket, directory) must be rejected up front
/// rather than clobbered. A directory is the portable stand-in for the whole
/// non-regular class — it exercises the same `is_file()` guard without needing a
/// platform-specific `mkfifo`.
#[test]
fn test_merge_rejects_non_regular_output_destination() {
    let dir = TempDir::new().unwrap();
    let header = create_coordinate_sorted_header("chr1", 10000);
    let a = dir.path().join("a.bam");
    let b = dir.path().join("b.bam");
    write_bam(&a, &header, &[mapped_record(b"reada0", 100)]);
    write_bam(&b, &header, &[mapped_record(b"readb0", 200)]);

    // A directory sitting where the output path points is a non-regular file.
    let out = dir.path().join("out_is_a_dir");
    fs::create_dir(&out).unwrap();

    let err = merge(&out, vec![a, b], SortOrderArg::Coordinate)
        .expect_err("must reject a non-regular output destination");
    assert!(
        err.to_string().contains("must be a regular file or stdout"),
        "unexpected error: {err}"
    );
    // The destination directory must be left untouched, and no temp staged.
    assert!(out.is_dir(), "the non-regular destination must be left in place");
    assert_no_leftover_temp(dir.path());
}

/// The atomic temp+persist must follow a symlinked output to its real target
/// (matching the pre-atomic-temp `File::create`, which follows symlinks) rather
/// than replacing the link with a regular file. A `latest.bam -> run.bam`
/// symlink must survive the merge, with the merged output landing in `run.bam`.
#[cfg(unix)]
#[test]
fn test_merge_output_follows_symlink_destination() {
    let dir = TempDir::new().unwrap();
    let header = create_coordinate_sorted_header("chr1", 10000);
    let a = dir.path().join("a.bam");
    let b = dir.path().join("b.bam");
    write_bam(&a, &header, &[mapped_record(b"reada0", 100)]);
    write_bam(&b, &header, &[mapped_record(b"readb0", 200)]);

    // A pre-existing real target with a relative symlink pointing at it.
    let real = dir.path().join("run.bam");
    write_bam(&real, &header, &[mapped_record(b"stale", 1)]);
    let link = dir.path().join("latest.bam");
    std::os::unix::fs::symlink("run.bam", &link).unwrap();

    merge(&link, vec![a, b], SortOrderArg::Coordinate).expect("valid coordinate merge");

    // The symlink itself must be preserved, not clobbered by a regular file...
    assert!(
        fs::symlink_metadata(&link).unwrap().file_type().is_symlink(),
        "symlink output must be preserved, not replaced by a regular file"
    );
    // ...and the merged records must land in the real target it points to.
    let mut reader = bam::io::Reader::new(fs::File::open(&real).unwrap());
    let _ = reader.read_header().unwrap();
    let names: Vec<String> = reader
        .records()
        .map(|r| String::from_utf8(r.unwrap().name().unwrap().to_vec()).unwrap())
        .collect();
    assert_eq!(
        names,
        vec!["reada0".to_string(), "readb0".to_string()],
        "merged records must land in the symlink's target file"
    );
}

/// MERGE3-01 streaming verify across every `--order`, not just coordinate: a
/// bare-header input whose records are out of order *for the requested order* is
/// caught record-by-record during the merge, the offending input is named, and no
/// partial output or temp is left behind. The mis-ordered pair is order-specific
/// (descending position for coordinate/template-coordinate, a descending read
/// name for the two queryname orders) so each case actually exercises that order's
/// key comparator rather than re-testing the coordinate path.
#[rstest::rstest]
#[case::coordinate(SortOrderArg::Coordinate, "coordinate", "first", 300, "second", 100)]
#[case::queryname(SortOrderArg::Queryname, "queryname", "nameB", 100, "nameA", 100)]
#[case::queryname_natural(
    SortOrderArg::QuerynameNatural,
    "queryname::natural",
    "read2",
    100,
    "read1",
    100
)]
#[case::template_coordinate(
    SortOrderArg::TemplateCoordinate,
    "template-coordinate",
    "first",
    300,
    "second",
    100
)]
fn test_merge_streaming_verify_rejects_missorted_for_each_order(
    #[case] order: SortOrderArg,
    #[case] order_label: &str,
    #[case] name0: &str,
    #[case] pos0: i32,
    #[case] name1: &str,
    #[case] pos1: i32,
) {
    let dir = TempDir::new().unwrap();
    // Bare header (declared order "none") routes records through the streaming
    // monotonicity check rather than the fast header-declared check.
    let bad = dir.path().join("bad.bam");
    write_bam(
        &bad,
        &bare_header(),
        &[mapped_record(name0.as_bytes(), pos0), mapped_record(name1.as_bytes(), pos1)],
    );
    let out = dir.path().join("merged.bam");

    let err = merge(&out, vec![bad], order)
        .expect_err("must reject an input mis-sorted for the requested order");
    let msg = err.to_string();
    assert!(
        msg.contains(&format!("not sorted in {order_label} order")),
        "error must name the {order_label} order: {msg}"
    );
    assert!(msg.contains("bad.bam"), "error should name the offending input: {msg}");
    assert!(!out.exists(), "a rejected merge must leave no output file");
    assert_no_leftover_temp(dir.path());
}

/// The valid-order counterpart to the per-order rejection table: a correctly
/// sorted bare-header input merges for every `--order`, and the records pass
/// through with their identities intact. Asserting the `(name, position)` oracle
/// — not just a record count — catches an order-specific key bug that reorders or
/// drops records while leaving the count unchanged.
#[rstest::rstest]
#[case::coordinate(SortOrderArg::Coordinate, "low", 100, "high", 300)]
#[case::queryname(SortOrderArg::Queryname, "nameA", 100, "nameB", 100)]
#[case::queryname_natural(SortOrderArg::QuerynameNatural, "read1", 100, "read2", 100)]
#[case::template_coordinate(SortOrderArg::TemplateCoordinate, "low", 100, "high", 300)]
fn test_merge_streaming_verify_accepts_sorted_for_each_order(
    #[case] order: SortOrderArg,
    #[case] name0: &str,
    #[case] pos0: i32,
    #[case] name1: &str,
    #[case] pos1: i32,
) {
    let dir = TempDir::new().unwrap();
    let input = dir.path().join("sorted.bam");
    write_bam(
        &input,
        &bare_header(),
        &[mapped_record(name0.as_bytes(), pos0), mapped_record(name1.as_bytes(), pos1)],
    );
    let out = dir.path().join("merged.bam");

    merge(&out, vec![input], order).expect("a correctly-ordered input must merge");
    assert!(out.exists());

    let got = read_name_pos(&out);
    assert_eq!(
        got,
        vec![
            (name0.to_string(), usize::try_from(pos0).unwrap()),
            (name1.to_string(), usize::try_from(pos1).unwrap()),
        ],
        "merged records must retain their identity and input order"
    );
}

// -----------------------------------------------------------------------------
// End-to-end output-order gate (drives the built `fgumi` binary): `fgumi merge`
// advertises `SS:template-coordinate` on its output, so build two genuinely
// template-coordinate sorted inputs, merge them, and assert the merged bytes
// verify as template-coordinate sorted via `fgumi sort --verify`.
// -----------------------------------------------------------------------------

/// A mapped F1R2 pair at `start` (0-based) on ref 0 with the given name/UMI.
fn mapped_pair(name: &str, umi: &str, start: i32) -> Vec<RawRecord> {
    let mut r1 = SamBuilder::new();
    r1.read_name(name.as_bytes())
        .sequence(b"ACGTACGT")
        .qualities(&[30; 8])
        .flags(flags::PAIRED | flags::FIRST_SEGMENT)
        .ref_id(0)
        .pos(start)
        .mapq(60)
        .cigar_ops(&[8 << 4])
        .mate_ref_id(0)
        .mate_pos(start + 100)
        .template_length(108)
        .add_string_tag(SamTag::RX, umi.as_bytes())
        .add_string_tag(SamTag::MC, b"8M");
    let mut r2 = SamBuilder::new();
    r2.read_name(name.as_bytes())
        .sequence(b"ACGTACGT")
        .qualities(&[30; 8])
        .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE)
        .ref_id(0)
        .pos(start + 100)
        .mapq(60)
        .cigar_ops(&[8 << 4])
        .mate_ref_id(0)
        .mate_pos(start)
        .template_length(-108)
        .add_string_tag(SamTag::RX, umi.as_bytes())
        .add_string_tag(SamTag::MC, b"8M");
    vec![r1.build(), r2.build()]
}

/// Write `records` to a temp BAM, then sort it into `path` in genuine
/// template-coordinate order with fgumi's own sorter.
fn write_template_coordinate_bam(dir: &Path, name: &str, records: Vec<RawRecord>) -> PathBuf {
    let unsorted = dir.join(format!("{name}.unsorted.bam"));
    let sorted = dir.join(format!("{name}.bam"));
    let header = create_minimal_header("chr1", 100_000);
    let mut writer = bam::io::Writer::new(fs::File::create(&unsorted).unwrap());
    writer.write_header(&header).unwrap();
    for record in records {
        writer.write_alignment_record(&header, &to_record_buf(&record)).unwrap();
    }
    writer.try_finish().unwrap();

    let status = std::process::Command::new(FGUMI)
        .args([
            "sort",
            "-i",
            unsorted.to_str().unwrap(),
            "-o",
            sorted.to_str().unwrap(),
            "--order",
            "template-coordinate",
            "--key-types",
            "mi",
        ])
        .status()
        .expect("spawn fgumi sort");
    assert!(status.success(), "failed to sort merge input {name}");
    sorted
}

#[test]
fn test_merge_output_is_template_coordinate_sorted() {
    let dir = TempDir::new().unwrap();
    // Two inputs whose positions interleave, so a correct merge must reorder across
    // files (not just concatenate).
    let in1 = write_template_coordinate_bam(
        dir.path(),
        "in1",
        [mapped_pair("a", "ACGTACGT", 200), mapped_pair("c", "ACGTACGT", 5_000)]
            .into_iter()
            .flatten()
            .collect(),
    );
    let in2 = write_template_coordinate_bam(
        dir.path(),
        "in2",
        [mapped_pair("b", "TGCATGCA", 1_000), mapped_pair("d", "TGCATGCA", 8_000)]
            .into_iter()
            .flatten()
            .collect(),
    );

    let merged = dir.path().join("merged.bam");
    let status = std::process::Command::new(FGUMI)
        .args([
            "merge",
            "-o",
            merged.to_str().unwrap(),
            "--order",
            "template-coordinate",
            in1.to_str().unwrap(),
            in2.to_str().unwrap(),
        ])
        .status()
        .expect("spawn fgumi merge");
    assert!(status.success(), "fgumi merge failed");
    assert!(merged.exists(), "merged BAM not created");

    // The whole point of merge is the output claim: verify it holds.
    assert_bam_sorted(&merged, "template-coordinate", Some("mi"));
}

// -----------------------------------------------------------------------------
// Inputs that reuse an @RG/@PG ID: identical records are combined, and a record
// whose ID is already used for different content gets a fresh ID, which that
// input's read tags follow.
// -----------------------------------------------------------------------------

/// A `chr1` header with no `@HD` sort order (so inputs are verified record by
/// record during the merge) plus the tab-separated @RG/@PG `lines`.
fn header_with_lines(lines: &[&str]) -> Header {
    let mut text = String::from("@SQ\tSN:chr1\tLN:10000\n");
    for line in lines {
        text.push_str(line);
        text.push('\n');
    }
    text.parse().expect("valid SAM header text")
}

/// A mapped 20M record on chr1 at the 1-based `pos`, with the given `RG`/`PG` tags.
fn tagged_record(name: &[u8], pos: i32, rg: Option<&str>, pg: Option<&str>) -> RawRecord {
    let mut b = SamBuilder::new();
    b.read_name(name)
        .ref_id(0)
        .pos(pos - 1)
        .mapq(60)
        .cigar_ops(&[20 << 4])
        .sequence(b"ACGTACGTACGTACGTACGT")
        .qualities(&[30; 20]);
    if let Some(rg) = rg {
        b.add_string_tag(SamTag::RG, rg.as_bytes());
    }
    if let Some(pg) = pg {
        b.add_string_tag(SamTag::PG, pg.as_bytes());
    }
    b.build()
}

/// `(read name, RG, PG)` for every record of `path`, in file order.
fn read_name_rg_pg(path: &Path) -> Vec<(String, Option<String>, Option<String>)> {
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::data::field::Value;
    let mut reader = bam::io::Reader::new(fs::File::open(path).unwrap());
    let header = reader.read_header().unwrap();
    let string_tag = |rec: &noodles::sam::alignment::RecordBuf, tag: SamTag| match rec
        .data()
        .get(&Tag::from(tag))
    {
        Some(Value::String(value)) => Some(value.to_string()),
        _ => None,
    };
    reader
        .record_bufs(&header)
        .map(|r| {
            let rec = r.unwrap();
            let name = String::from_utf8(rec.name().unwrap().to_vec()).unwrap();
            (name, string_tag(&rec, SamTag::RG), string_tag(&rec, SamTag::PG))
        })
        .collect()
}

/// The @PG line `fgumi merge` (run by [`merge`]) adds, with ID `id`, chained to `previous`.
fn merge_pg_line(id: &str, previous: Option<&str>) -> String {
    merge_pg_line_with_command_line(id, "fgumi merge", previous)
}

/// [`merge_pg_line`] with `command_line` as the recorded `CL`.
fn merge_pg_line_with_command_line(id: &str, command_line: &str, previous: Option<&str>) -> String {
    let line = format!(
        "@PG\tID:{id}\tPN:fgumi\tVN:{}\tCL:{command_line}",
        fgumi_lib::version::VERSION.as_str()
    );
    match previous {
        Some(pp) => format!("{line}\tPP:{pp}"),
        None => line,
    }
}

/// The @RG and @PG lines of the header of `path`, in header order, as SAM text.
fn read_group_and_program_lines(path: &Path) -> Vec<String> {
    let mut reader = bam::io::Reader::new(fs::File::open(path).unwrap());
    let header = reader.read_header().unwrap();
    let mut writer = noodles::sam::io::Writer::new(Vec::new());
    writer.write_header(&header).unwrap();
    String::from_utf8(writer.get_ref().clone())
        .unwrap()
        .lines()
        .filter(|line| line.starts_with("@RG") || line.starts_with("@PG"))
        .map(str::to_string)
        .collect()
}

/// Two inputs reuse `RG:A` and `PG:bwa` for different read groups and programs: the
/// second input's are renamed `A.1` and `bwa.1`, and its reads' tags follow, in every
/// merge order. A tag naming an ID its input's header does not declare is kept.
#[rstest::rstest]
#[case::coordinate(SortOrderArg::Coordinate)]
#[case::queryname(SortOrderArg::Queryname)]
#[case::queryname_natural(SortOrderArg::QuerynameNatural)]
#[case::template_coordinate(SortOrderArg::TemplateCoordinate)]
fn test_merge_rewrites_read_tags_of_conflicting_ids(#[case] order: SortOrderArg) {
    let dir = TempDir::new().unwrap();
    let in1 = dir.path().join("in1.bam");
    let in2 = dir.path().join("in2.bam");
    write_bam(
        &in1,
        &header_with_lines(&["@RG\tID:A\tLB:libX", "@PG\tID:bwa\tCL:bwa ref1.fa"]),
        &[tagged_record(b"r1", 100, Some("A"), Some("bwa"))],
    );
    write_bam(
        &in2,
        &header_with_lines(&["@RG\tID:A\tLB:libY", "@PG\tID:bwa\tCL:bwa ref2.fa"]),
        &[
            tagged_record(b"r2", 200, Some("A"), Some("bwa")),
            tagged_record(b"r3", 300, Some("undeclared"), None),
        ],
    );

    let merged = dir.path().join("merged.bam");
    merge(&merged, vec![in1, in2], order).expect("merge should succeed");

    assert_eq!(
        read_group_and_program_lines(&merged),
        [
            "@RG\tID:A\tLB:libX".to_string(),
            "@RG\tID:A.1\tLB:libY".to_string(),
            "@PG\tID:bwa\tCL:bwa ref1.fa".to_string(),
            "@PG\tID:bwa.1\tCL:bwa ref2.fa".to_string(),
            merge_pg_line("fgumi", Some("bwa.1")),
        ]
    );
    let mut tags = read_name_rg_pg(&merged);
    tags.sort();
    let owned = |name: &str, rg: Option<&str>, pg: Option<&str>| {
        (name.to_string(), rg.map(str::to_string), pg.map(str::to_string))
    };
    assert_eq!(
        tags,
        [
            owned("r1", Some("A"), Some("bwa")),
            owned("r2", Some("A.1"), Some("bwa.1")),
            owned("r3", Some("undeclared"), None),
        ]
    );
}

/// Shards of one sample share identical @RG/@PG records: each is written once and the
/// reads' tags are unchanged (`samtools merge` would rename them without `-c`/`-p`).
#[test]
fn test_merge_combines_identical_shard_headers() {
    let dir = TempDir::new().unwrap();
    let lines = ["@RG\tID:A\tLB:libX", "@PG\tID:bwa\tCL:bwa ref.fa"];
    let in1 = dir.path().join("in1.bam");
    let in2 = dir.path().join("in2.bam");
    write_bam(
        &in1,
        &header_with_lines(&lines),
        &[tagged_record(b"r1", 100, Some("A"), Some("bwa"))],
    );
    write_bam(
        &in2,
        &header_with_lines(&lines),
        &[tagged_record(b"r2", 200, Some("A"), Some("bwa"))],
    );

    let merged = dir.path().join("merged.bam");
    merge(&merged, vec![in1, in2], SortOrderArg::Coordinate).expect("merge should succeed");

    let mut expected_lines: Vec<String> = lines.iter().map(|line| (*line).to_string()).collect();
    expected_lines.push(merge_pg_line("fgumi", Some("bwa")));
    assert_eq!(read_group_and_program_lines(&merged), expected_lines);
    let tags = read_name_rg_pg(&merged);
    let expected: Vec<_> = ["r1", "r2"]
        .iter()
        .map(|name| (name.to_string(), Some("A".to_string()), Some("bwa".to_string())))
        .collect();
    assert_eq!(tags, expected);
}

/// The template-coordinate key orders templates at the same position by library. Both
/// inputs call their read group `A`, with libraries `libA` and `libZ`; once the second
/// input's reads are retagged `A.1` (`libZ`), its read `b` sorts after `a` (`libA`) at
/// position 300. Were `b` keyed by the first input's `A` (`libA`), the two would tie on
/// library and the name hash, which puts `b` first, would decide. The `refill` case puts
/// another record before each tie, so `a` and `b` are each read after their input's
/// first record rather than as it.
#[rstest::rstest]
#[case::first_record(&[], &[], &["a", "b"])]
#[case::refill(&[("c", 100)], &[("d", 50)], &["d", "c", "a", "b"])]
fn test_merge_template_coordinate_orders_by_renamed_read_group_library(
    #[case] in1_before: &[(&str, i32)],
    #[case] in2_before: &[(&str, i32)],
    #[case] expected: &[&str],
) {
    let dir = TempDir::new().unwrap();
    let records = |before: &[(&str, i32)], tied: &str| -> Vec<RawRecord> {
        before
            .iter()
            .map(|(name, pos)| tagged_record(name.as_bytes(), *pos, Some("A"), None))
            .chain(std::iter::once(tagged_record(tied.as_bytes(), 300, Some("A"), None)))
            .collect()
    };
    let in1 = dir.path().join("in1.bam");
    let in2 = dir.path().join("in2.bam");
    write_bam(&in1, &header_with_lines(&["@RG\tID:A\tLB:libA"]), &records(in1_before, "a"));
    write_bam(&in2, &header_with_lines(&["@RG\tID:A\tLB:libZ"]), &records(in2_before, "b"));

    let merged = dir.path().join("merged.bam");
    merge(&merged, vec![in1, in2], SortOrderArg::TemplateCoordinate).expect("merge should succeed");

    let names: Vec<String> = read_name_rg_pg(&merged).into_iter().map(|(name, ..)| name).collect();
    assert_eq!(names, expected);
}

/// `fgumi merge` adds one @PG of its own, with its command line, chained to the newest
/// program in the merged header (the last chain end in header order), as every other
/// fgumi command does, and under a fresh ID when an input already has `ID:fgumi`. A
/// single input is recorded too.
#[rstest::rstest]
#[case::no_input_programs(vec![vec![], vec![]], vec![merge_pg_line("fgumi", None)])]
#[case::single_input(
    vec![vec!["@PG\tID:bwa\tCL:bwa ref.fa"]],
    vec!["@PG\tID:bwa\tCL:bwa ref.fa".to_string(), merge_pg_line("fgumi", Some("bwa"))],
)]
#[case::chains_to_newest_program(
    vec![
        vec!["@PG\tID:bwa\tCL:bwa ref.fa", "@PG\tID:samtools\tPP:bwa\tCL:samtools sort"],
        vec!["@PG\tID:bwa\tCL:bwa ref.fa", "@PG\tID:samtools\tPP:bwa\tCL:samtools sort"],
    ],
    vec![
        "@PG\tID:bwa\tCL:bwa ref.fa".to_string(),
        "@PG\tID:samtools\tPP:bwa\tCL:samtools sort".to_string(),
        merge_pg_line("fgumi", Some("samtools")),
    ],
)]
#[case::fresh_id_after_fgumi_programs(
    vec![
        vec!["@PG\tID:bwa\tCL:bwa ref.fa", "@PG\tID:fgumi\tPP:bwa\tCL:fgumi zipper"],
        vec!["@PG\tID:bwa\tCL:bwa ref.fa", "@PG\tID:fgumi\tPP:bwa\tCL:fgumi zipper"],
    ],
    vec![
        "@PG\tID:bwa\tCL:bwa ref.fa".to_string(),
        "@PG\tID:fgumi\tPP:bwa\tCL:fgumi zipper".to_string(),
        merge_pg_line("fgumi.1", Some("fgumi")),
    ],
)]
// Sharded fgumi runs: the inputs' fgumi @PGs differ, so the second is renamed `fgumi.1`
// and the merge's own @PG takes the next free ID and chains to it.
#[case::fresh_id_after_conflicting_fgumi_programs(
    vec![vec!["@PG\tID:fgumi\tCL:fgumi zipper a"], vec!["@PG\tID:fgumi\tCL:fgumi zipper b"]],
    vec![
        "@PG\tID:fgumi\tCL:fgumi zipper a".to_string(),
        "@PG\tID:fgumi.1\tCL:fgumi zipper b".to_string(),
        merge_pg_line("fgumi.2", Some("fgumi.1")),
    ],
)]
fn test_merge_adds_its_own_program_record(
    #[case] inputs: Vec<Vec<&str>>,
    #[case] expected_programs: Vec<String>,
) {
    let dir = TempDir::new().unwrap();
    let paths: Vec<PathBuf> = inputs
        .iter()
        .enumerate()
        .map(|(i, lines)| {
            let path = dir.path().join(format!("in{i}.bam"));
            let pos = 100 * (i32::try_from(i).unwrap() + 1);
            let name = format!("r{i}");
            write_bam(
                &path,
                &header_with_lines(lines),
                &[tagged_record(name.as_bytes(), pos, None, None)],
            );
            path
        })
        .collect();

    let merged = dir.path().join("merged.bam");
    merge(&merged, paths, SortOrderArg::Coordinate).expect("merge should succeed");

    assert_eq!(read_group_and_program_lines(&merged), expected_programs);
}

/// The merge's @PG records the command line it was run with, not a fixed string, with
/// the tab a SAM header cannot hold turned into a space rather than failing the output.
#[test]
fn test_merge_program_record_records_command_line() {
    let dir = TempDir::new().unwrap();
    let paths: Vec<PathBuf> = (0..2)
        .map(|i| {
            let path = dir.path().join(format!("in{i}.bam"));
            let name = format!("r{i}");
            write_bam(
                &path,
                &header_with_lines(&[]),
                &[tagged_record(name.as_bytes(), 100 * (i + 1), None, None)],
            );
            path
        })
        .collect();

    let merged = dir.path().join("merged.bam");
    merge_with_command_line(
        &merged,
        paths,
        SortOrderArg::Coordinate,
        "fgumi merge --order coordinate\t-o merged.bam in0.bam in1.bam",
    )
    .expect("merge should succeed");

    assert_eq!(
        read_group_and_program_lines(&merged),
        [merge_pg_line_with_command_line(
            "fgumi",
            "fgumi merge --order coordinate -o merged.bam in0.bam in1.bam",
            None,
        )]
    );
}

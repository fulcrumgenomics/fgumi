//! Fused-mode `--group::*` metrics parity: runall must produce the same
//! family-size-histogram/grouping-metrics output as standalone `fgumi group`
//! with the same flags, once the anti-goal block is removed.

use crate::helpers::bam_generator::{create_minimal_header, write_bam};
use fgumi_raw_bam::{RawRecord, SamBuilder, flags};
use rstest::rstest;
use std::process::Command;

/// Runs the compiled `fgumi` binary as a subprocess, matching this repo's
/// real integration-test convention (`tests/integration/helpers/cli.rs`).
fn run_fgumi(args: &[&str]) {
    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args(args)
        .status()
        .unwrap_or_else(|e| panic!("failed to spawn fgumi {args:?}: {e}"));
    assert!(status.success(), "fgumi {args:?} failed with {status}");
}

/// A single-end, mapped read carrying an `RX` UMI tag, suitable as `group`
/// input under the `adjacency` strategy.
fn grouped_record(name: &str, pos: i32, umi: &str) -> RawRecord {
    let mut b = SamBuilder::new();
    b.read_name(name.as_bytes())
        .sequence(b"ACGTACGTAC")
        .qualities(&[30; 10])
        .flags(0)
        .ref_id(0)
        .pos(pos)
        .mapq(60)
        .cigar_ops(&[10 << 4]);
    b.add_string_tag(fgumi_raw_bam::SamTag::RX, umi.as_bytes());
    b.build()
}

#[rstest]
#[case::single_thread(1)]
#[case::multi_thread(4)]
fn runall_group_family_size_histogram_matches_standalone_group(#[case] threads: usize) {
    let dir = tempfile::tempdir().expect("tempdir");
    let header = create_minimal_header("chr1", 10_000);
    let bam_path = dir.path().join("in.bam");
    write_bam(&bam_path, &header, &[grouped_record("r1", 99, "AAAAAAAA")]);

    let fused_histogram = dir.path().join("fused.family_sizes.txt");
    run_fgumi(&[
        "runall",
        "-i",
        bam_path.to_str().unwrap(),
        "-o",
        dir.path().join("fused_out.bam").to_str().unwrap(),
        "--start-from",
        "group",
        "--stop-after",
        "group",
        "--group::strategy",
        "adjacency",
        "--group::family-size-histogram",
        fused_histogram.to_str().unwrap(),
        "--threads",
        &threads.to_string(),
    ]);

    let standalone_histogram = dir.path().join("standalone.family_sizes.txt");
    run_fgumi(&[
        "group",
        "-i",
        bam_path.to_str().unwrap(),
        "-o",
        dir.path().join("standalone_out.bam").to_str().unwrap(),
        "-s",
        "adjacency",
        "--family-size-histogram",
        standalone_histogram.to_str().unwrap(),
    ]);

    assert_files_eq(
        &fused_histogram,
        &standalone_histogram,
        "group family-size histogram: fused vs standalone",
    );
}

#[test]
fn runall_group_metrics_and_consensus_metrics_coexist_in_one_fused_run() {
    let dir = tempfile::tempdir().expect("tempdir");
    let header = create_minimal_header("chr1", 10_000);
    let bam_path = dir.path().join("in.bam");
    write_bam(&bam_path, &header, &[grouped_record("r1", 99, "AAAAAAAA")]);

    let group_histogram = dir.path().join("group.family_sizes.txt");
    let simplex_prefix = dir.path().join("simplex");
    run_fgumi(&[
        "runall",
        "-i",
        bam_path.to_str().unwrap(),
        "-o",
        dir.path().join("out.bam").to_str().unwrap(),
        "--start-from",
        "group",
        "--stop-after",
        "consensus",
        "--consensus",
        "simplex",
        "--group::strategy",
        "adjacency",
        "--group::family-size-histogram",
        group_histogram.to_str().unwrap(),
        "--simplex::min-reads",
        "1",
        "--simplex::metrics",
        simplex_prefix.to_str().unwrap(),
    ]);

    assert!(group_histogram.is_file(), "group's own metrics file must exist");
    assert!(
        dir.path().join("simplex.family_sizes.txt").is_file(),
        "simplex's inline metrics file must exist independently of group's"
    );
}

/// `--all-metrics` must never overwrite a metrics path the user set
/// explicitly on a per-stage flag — it only fills in options still `None`.
///
/// Uses `--group::metrics` (the `metrics_prefix` field) rather than
/// `--group::family-size-histogram`: `--all-metrics` only ever derives
/// group's `metrics_prefix` (see the corrected fill-in and
/// `all_metrics_fills_in_every_applicable_stage_present_in_the_chain`'s doc
/// comment below), so `--group::family-size-histogram` is never a path
/// `--all-metrics` would derive in the first place — asserting explicit-wins
/// against it would be vacuous.
#[test]
fn all_metrics_never_overwrites_an_explicit_per_stage_flag() {
    let dir = tempfile::tempdir().expect("tempdir");
    let header = create_minimal_header("chr1", 10_000);
    let bam_path = dir.path().join("in.bam");
    write_bam(&bam_path, &header, &[grouped_record("r1", 100, "AAAAAAAA")]);

    let explicit_prefix = dir.path().join("explicit");
    let all_metrics_prefix = dir.path().join("all");

    run_fgumi(&[
        "runall",
        "-i",
        bam_path.to_str().unwrap(),
        "-o",
        dir.path().join("out.bam").to_str().unwrap(),
        "--start-from",
        "group",
        "--stop-after",
        "group",
        "--group::strategy",
        "adjacency",
        "--group::metrics",
        explicit_prefix.to_str().unwrap(),
        "--all-metrics",
        all_metrics_prefix.to_str().unwrap(),
    ]);

    assert!(
        dir.path().join("explicit.family_sizes.txt").is_file(),
        "explicit flag's prefix must be used"
    );
    assert!(
        !dir.path().join("all.group.family_sizes.txt").exists(),
        "the derived prefix must NOT also be written — explicit wins outright"
    );
}

/// `--all-metrics` fills in every applicable stage's metrics option(s) for
/// every stage present in the chain, deriving each path/prefix from the
/// single `--all-metrics` prefix, and touches nothing for stages absent
/// from this chain (correct/filter/duplex/codec here).
///
/// Uses `--consensus simplex` over single-end input — a shape already
/// proven safe by
/// `runall_group_metrics_and_consensus_metrics_coexist_in_one_fused_run`
/// above — purely to exercise the same `derived_metrics_prefix`/
/// `derived_metrics_path` fill-in code for every stage type; it is not
/// standing in for anything paired-end-specific. (An earlier version of
/// this comment noted that paired-end input combined with a fused
/// per-stage consensus `--metrics` flag hit a T1 inline-metrics-tap bug —
/// `coordinate_group_from_processed_position` read the `MI` SAM tag off
/// the raw record before `MiAssignGroups` had written it. That bug is
/// fixed; see `runall_group_to_codec_metrics_succeeds_on_paired_end_input`
/// below for the paired-end regression coverage.)
#[test]
fn all_metrics_fills_in_every_applicable_stage_present_in_the_chain() {
    let dir = tempfile::tempdir().expect("tempdir");
    let header = create_minimal_header("chr1", 10_000);
    let bam_path = dir.path().join("in.bam");
    write_bam(&bam_path, &header, &[grouped_record("r1", 99, "AAAAAAAA")]);
    let prefix = dir.path().join("all");

    run_fgumi(&[
        "runall",
        "-i",
        bam_path.to_str().unwrap(),
        "-o",
        dir.path().join("out.bam").to_str().unwrap(),
        "--start-from",
        "group",
        "--stop-after",
        "consensus",
        "--consensus",
        "simplex",
        "--group::strategy",
        "adjacency",
        "--simplex::min-reads",
        "1",
        "--all-metrics",
        prefix.to_str().unwrap(),
    ]);

    // `--all-metrics` derives ONLY group's `metrics_prefix` (`all.group`);
    // its fan-out then writes the complete canonical group metric set —
    // Task 1's three pinned filenames — with no collision.
    assert!(
        dir.path().join("all.group.family_sizes.txt").is_file(),
        "group's own --metrics prefix output"
    );
    assert!(dir.path().join("all.group.grouping_metrics.txt").is_file());
    assert!(dir.path().join("all.group.position_group_sizes.txt").is_file());
    assert!(
        !dir.path().join("all.group.family_size_histogram.txt").exists(),
        "--all-metrics must NOT separately derive --group::family-size-histogram: it would be a \
         redundant near-duplicate of the metrics_prefix fan-out's family_sizes.txt"
    );
    assert!(dir.path().join("all.simplex.family_sizes.txt").is_file());
    assert!(dir.path().join("all.simplex.umi_counts.txt").is_file());
    // No correct/filter/duplex/codec files — those stages are not in this chain.
    assert!(!dir.path().join("all.correct.metrics.txt").exists());
    assert!(!dir.path().join("all.filter.stats.txt").exists());
    assert!(!dir.path().join("all.duplex.family_sizes.txt").exists());
    assert!(!dir.path().join("all.codec.family_sizes.txt").exists());
}

/// A paired-end, mapped R1/R2 pair carrying a dual (`-`-delimited) `RX` UMI
/// tag and no `MI` tag — the shape `group` consumes as ungrouped input.
/// `MI` is deliberately absent: it is `group`'s job to assign it, and the
/// whole point of the regression test below is that the fused pipeline must
/// not need it on the record before that point.
fn paired_record(name: &str, pos: i32, umi: &str) -> (RawRecord, RawRecord) {
    let seq = b"ACGTACGTAC";
    let quals = [30u8; 10];

    let mut b1 = SamBuilder::new();
    b1.read_name(name.as_bytes())
        .sequence(seq)
        .qualities(&quals)
        .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE)
        .ref_id(0)
        .pos(pos)
        .mapq(60)
        .cigar_ops(&[10 << 4]) // 10M
        .mate_ref_id(0)
        .mate_pos(pos + 10)
        .template_length(20);
    b1.add_string_tag(fgumi_raw_bam::SamTag::RX, umi.as_bytes());
    b1.add_string_tag(fgumi_raw_bam::SamTag::MC, b"10M");
    let r1 = b1.build();

    let mut b2 = SamBuilder::new();
    b2.read_name(name.as_bytes())
        .sequence(seq)
        .qualities(&quals)
        .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE)
        .ref_id(0)
        .pos(pos + 10)
        .mapq(60)
        .cigar_ops(&[10 << 4]) // 10M
        .mate_ref_id(0)
        .mate_pos(pos)
        .template_length(-20);
    b2.add_string_tag(fgumi_raw_bam::SamTag::RX, umi.as_bytes());
    b2.add_string_tag(fgumi_raw_bam::SamTag::MC, b"10M");
    let r2 = b2.build();

    (r1, r2)
}

/// Regression test for the verified T1 bug: the fused inline-metrics tap
/// (`coordinate_group_from_processed_position`, called from
/// `build_group_process_step` in `chains/commands/group.rs`) used to read
/// `mi` from the record's `MI` SAM aux tag via `build_template_info`. At
/// that point in a fused `runall` pipeline, `group` has already assigned
/// the in-memory `Template.mi` FIELD, but the `MI` tag is not written onto
/// the record until BAM serialization runs later — so any fused
/// `runall --start-from group --consensus <mode> --<mode>::metrics=<prefix>`
/// run over PAIRED-END (multi-record-template) input failed with "missing
/// the required MI tag" (repro: plain `--codec::metrics`). This is exactly
/// that scenario: paired-end input, `--codec::metrics` set, fused
/// `group -> codec`.
///
/// Against the pre-fix adapter (which called `build_template_info`,
/// unconditionally reading the `MI` aux tag) this `run_fgumi` call panics —
/// `fgumi runall ... failed with exit status: 1` — because the process
/// exits non-zero. With the fix (`build_template_info_with_mi`, sourcing
/// `mi` from `Template.mi`) the run succeeds and produces codec's usual
/// inline-metrics file set (codec's finalize hook uses
/// `MetricsThresholds::Duplex`, so it writes the same file set
/// `duplex-metrics` writes).
#[test]
fn runall_group_to_codec_metrics_succeeds_on_paired_end_input() {
    let dir = tempfile::tempdir().expect("tempdir");
    let header = create_minimal_header("chr1", 10_000);
    let bam_path = dir.path().join("in.bam");

    let (r1, r2) = paired_record("r1", 99, "AAAAAAAA-TTTTTTTT");
    write_bam(&bam_path, &header, &[r1, r2]);

    let codec_prefix = dir.path().join("codec");
    run_fgumi(&[
        "runall",
        "-i",
        bam_path.to_str().unwrap(),
        "-o",
        dir.path().join("out.bam").to_str().unwrap(),
        "--start-from",
        "group",
        "--stop-after",
        "consensus",
        "--consensus",
        "codec",
        "--group::strategy",
        "adjacency",
        "--codec::min-reads",
        "1",
        "--codec::metrics",
        codec_prefix.to_str().unwrap(),
    ]);

    assert!(dir.path().join("codec.family_sizes.txt").is_file(), "family_sizes.txt must exist");
    assert!(
        dir.path().join("codec.duplex_family_sizes.txt").is_file(),
        "duplex_family_sizes.txt must exist"
    );
    assert!(dir.path().join("codec.umi_counts.txt").is_file(), "umi_counts.txt must exist");
    assert!(
        dir.path().join("codec.duplex_yield_metrics.txt").is_file(),
        "duplex_yield_metrics.txt must exist"
    );
}

/// Assert two metrics files exist and are byte-for-byte identical, and that the
/// comparison is non-vacuous.
///
/// Both sides of these comparisons come from the same `fgumi` binary, so an
/// empty-vs-empty match would pass without pinning anything. Every metrics
/// artifact compared through this helper is a line-oriented TSV whose first line
/// is a header and whose data begins on the second line, so a file with fewer
/// than two lines is header-only — the signature of a silently-broken input
/// (e.g. single-end reads feeding the simplex metrics path, which counts paired
/// templates only, leaving `simplex_yield_metrics` / `umi_counts` empty). Require
/// at least one data row on *every* compared file so no degenerate fixture can
/// make a comparison vacuous. (A bare "file is non-empty" check would be useless:
/// every metrics writer emits a header regardless.)
fn assert_files_eq(actual: &std::path::Path, expected: &std::path::Path, ctx: &str) {
    let a = std::fs::read_to_string(actual)
        .unwrap_or_else(|e| panic!("{ctx}: missing actual {}: {e}", actual.display()));
    let e = std::fs::read_to_string(expected)
        .unwrap_or_else(|e| panic!("{ctx}: missing expected {}: {e}", expected.display()));
    assert!(
        e.lines().count() >= 2,
        "{ctx}: {} has no data rows — the comparison would be vacuous",
        expected.display(),
    );
    assert_eq!(a, e, "{ctx}: {} diverges from {}", actual.display(), expected.display());
}

/// A small multi-family, multi-position input, so the group and consensus
/// metrics files carry non-trivial content (family sizes 3, 2, 1 across two
/// coordinate groups) rather than a single degenerate family. Uses real
/// PAIRED templates: group metrics count families fine from single-end reads,
/// but the SIMPLEX consensus metrics count paired templates only, so a
/// single-end input would leave `<prefix>.simplex.*` empty and make the
/// consensus half of the content comparisons vacuous.
fn multi_family_records() -> Vec<RawRecord> {
    let mut records = Vec::new();
    for (name, pos, umi) in [
        ("a1", 100, "AAAAAAAA"),
        ("a2", 100, "AAAAAAAA"),
        ("a3", 100, "AAAAAAAA"),
        ("b1", 100, "CCCCCCCC"),
        ("b2", 100, "CCCCCCCC"),
        ("c1", 500, "GGGGGGGG"),
    ] {
        let (r1, r2) = paired_record(name, pos, umi);
        records.push(r1);
        records.push(r2);
    }
    records
}

/// `--all-metrics` must produce byte-identical metric files to setting the
/// equivalent per-stage flags explicitly — the aggregator only fills the same
/// options those flags set, so their *content* must match, not merely their
/// existence. (The existing `all_metrics_fills_in_*` test asserts creation;
/// this closes the content gap.) Covers group's full three-file set and
/// simplex's three files including `simplex_yield_metrics.txt`.
#[test]
fn all_metrics_content_matches_explicit_per_stage_flags() {
    let dir = tempfile::tempdir().expect("tempdir");
    let header = create_minimal_header("chr1", 10_000);
    let bam = dir.path().join("in.bam");
    write_bam(&bam, &header, &multi_family_records());

    // Run 1: one --all-metrics prefix fills every present stage's metrics.
    let all = dir.path().join("all");
    run_fgumi(&[
        "runall",
        "-i",
        bam.to_str().unwrap(),
        "-o",
        dir.path().join("o1.bam").to_str().unwrap(),
        "--start-from",
        "group",
        "--stop-after",
        "consensus",
        "--consensus",
        "simplex",
        "--group::strategy",
        "adjacency",
        "--simplex::min-reads",
        "1",
        "--all-metrics",
        all.to_str().unwrap(),
    ]);

    // Run 2: the same pipeline, metrics requested via explicit per-stage flags.
    let g = dir.path().join("g");
    let s = dir.path().join("s");
    run_fgumi(&[
        "runall",
        "-i",
        bam.to_str().unwrap(),
        "-o",
        dir.path().join("o2.bam").to_str().unwrap(),
        "--start-from",
        "group",
        "--stop-after",
        "consensus",
        "--consensus",
        "simplex",
        "--group::strategy",
        "adjacency",
        "--simplex::min-reads",
        "1",
        "--group::metrics",
        g.to_str().unwrap(),
        "--simplex::metrics",
        s.to_str().unwrap(),
    ]);

    for suffix in ["family_sizes.txt", "grouping_metrics.txt", "position_group_sizes.txt"] {
        assert_files_eq(
            &dir.path().join(format!("all.group.{suffix}")),
            &dir.path().join(format!("g.{suffix}")),
            "--all-metrics group output vs explicit --group::metrics",
        );
    }
    for suffix in ["family_sizes.txt", "simplex_yield_metrics.txt", "umi_counts.txt"] {
        assert_files_eq(
            &dir.path().join(format!("all.simplex.{suffix}")),
            &dir.path().join(format!("s.{suffix}")),
            "--all-metrics simplex output vs explicit --simplex::metrics",
        );
    }
}

/// Fused `runall --group::metrics <prefix>` must produce the same three-file
/// group metric set (`family_sizes` / `grouping_metrics` /
/// `position_group_sizes`) as standalone `fgumi group --metrics <prefix>`.
/// The existing fused-vs-standalone parity test only compares the
/// `--family-size-histogram` file; this closes the `--metrics`-prefix content
/// gap, single- and multi-threaded.
#[rstest]
#[case::single_thread(1)]
#[case::multi_thread(4)]
fn runall_group_metrics_prefix_matches_standalone_group(#[case] threads: usize) {
    let dir = tempfile::tempdir().expect("tempdir");
    let header = create_minimal_header("chr1", 10_000);
    let bam = dir.path().join("in.bam");
    write_bam(&bam, &header, &multi_family_records());

    let fused = dir.path().join("fused");
    run_fgumi(&[
        "runall",
        "-i",
        bam.to_str().unwrap(),
        "-o",
        dir.path().join("fused_out.bam").to_str().unwrap(),
        "--start-from",
        "group",
        "--stop-after",
        "group",
        "--group::strategy",
        "adjacency",
        "--group::metrics",
        fused.to_str().unwrap(),
        "--threads",
        &threads.to_string(),
    ]);

    let standalone = dir.path().join("standalone");
    run_fgumi(&[
        "group",
        "-i",
        bam.to_str().unwrap(),
        "-o",
        dir.path().join("standalone_out.bam").to_str().unwrap(),
        "-s",
        "adjacency",
        "--metrics",
        standalone.to_str().unwrap(),
    ]);

    for suffix in ["family_sizes.txt", "grouping_metrics.txt", "position_group_sizes.txt"] {
        assert_files_eq(
            &dir.path().join(format!("fused.{suffix}")),
            &dir.path().join(format!("standalone.{suffix}")),
            "fused --group::metrics vs standalone group --metrics",
        );
    }
}

/// An unmapped single-end read carrying an `RX` UMI — the shape the `correct`
/// stage consumes (`--start-from correct`). The correction itself is
/// incidental here; the point is that `--all-metrics` fills in and writes
/// correct's metrics file when a correct stage is in the chain.
fn unmapped_umi_record(name: &str, umi: &str) -> RawRecord {
    let mut b = SamBuilder::new();
    b.read_name(name.as_bytes())
        .sequence(b"ACGTACGTAC")
        .qualities(&[30; 10])
        .flags(flags::UNMAPPED);
    b.add_string_tag(fgumi_raw_bam::SamTag::RX, umi.as_bytes());
    b.build()
}

/// `--all-metrics` must fill in and write the CORRECT stage's metrics when a
/// correct stage is present, with content identical to an explicit
/// `--correct::metrics`. Existing coverage only asserts correct's file is
/// ABSENT when the stage is out of the chain (see
/// `all_metrics_fills_in_every_applicable_stage_present_in_the_chain`); this
/// is the positive case.
#[test]
fn all_metrics_fills_correct_metrics_when_correct_stage_present() {
    let dir = tempfile::tempdir().expect("tempdir");
    let header = create_minimal_header("chr1", 10_000);
    let bam = dir.path().join("unmapped.bam");
    write_bam(
        &bam,
        &header,
        &[
            unmapped_umi_record("r1", "AAAAAAAA"),
            unmapped_umi_record("r2", "AAAAAAAT"),
            unmapped_umi_record("r3", "CCCCCCCC"),
        ],
    );
    let whitelist = dir.path().join("umis.txt");
    std::fs::write(&whitelist, "AAAAAAAA\nCCCCCCCC\n").expect("write whitelist");

    let all = dir.path().join("all");
    run_fgumi(&[
        "runall",
        "-i",
        bam.to_str().unwrap(),
        "-o",
        dir.path().join("o1.bam").to_str().unwrap(),
        "--start-from",
        "correct",
        "--stop-after",
        "correct",
        "--correct::umi-files",
        whitelist.to_str().unwrap(),
        "--correct::max-mismatches",
        "1",
        "--correct::min-distance",
        "1",
        "--all-metrics",
        all.to_str().unwrap(),
    ]);
    let all_correct = dir.path().join("all.correct.metrics.txt");
    assert!(all_correct.is_file(), "--all-metrics must write correct's metrics file");

    let explicit = dir.path().join("explicit.correct.txt");
    run_fgumi(&[
        "runall",
        "-i",
        bam.to_str().unwrap(),
        "-o",
        dir.path().join("o2.bam").to_str().unwrap(),
        "--start-from",
        "correct",
        "--stop-after",
        "correct",
        "--correct::umi-files",
        whitelist.to_str().unwrap(),
        "--correct::max-mismatches",
        "1",
        "--correct::min-distance",
        "1",
        "--correct::metrics",
        explicit.to_str().unwrap(),
    ]);
    assert_files_eq(
        &all_correct,
        &explicit,
        "--all-metrics correct vs explicit --correct::metrics",
    );
}

/// `--all-metrics` must fill in and write the FILTER stage's stats when a
/// filter stage is present, with content identical to an explicit
/// `--filter::stats`. The consensus BAM the filter stage needs is produced by
/// a preliminary fused group -> simplex run.
#[test]
fn all_metrics_fills_filter_stats_when_filter_stage_present() {
    let dir = tempfile::tempdir().expect("tempdir");
    let header = create_minimal_header("chr1", 10_000);
    let bam = dir.path().join("in.bam");
    write_bam(&bam, &header, &multi_family_records());

    // A consensus BAM to feed the filter-only chain.
    let consensus = dir.path().join("consensus.bam");
    run_fgumi(&[
        "runall",
        "-i",
        bam.to_str().unwrap(),
        "-o",
        consensus.to_str().unwrap(),
        "--start-from",
        "group",
        "--stop-after",
        "consensus",
        "--consensus",
        "simplex",
        "--group::strategy",
        "adjacency",
        "--simplex::min-reads",
        "1",
    ]);

    let all = dir.path().join("all");
    run_fgumi(&[
        "runall",
        "-i",
        consensus.to_str().unwrap(),
        "-o",
        dir.path().join("f1.bam").to_str().unwrap(),
        "--start-from",
        "filter",
        "--stop-after",
        "filter",
        "--filter::min-reads",
        "1",
        "--all-metrics",
        all.to_str().unwrap(),
    ]);
    let all_stats = dir.path().join("all.filter.stats.txt");
    assert!(all_stats.is_file(), "--all-metrics must write filter's stats file");

    let explicit = dir.path().join("explicit.filter.txt");
    run_fgumi(&[
        "runall",
        "-i",
        consensus.to_str().unwrap(),
        "-o",
        dir.path().join("f2.bam").to_str().unwrap(),
        "--start-from",
        "filter",
        "--stop-after",
        "filter",
        "--filter::min-reads",
        "1",
        "--filter::stats",
        explicit.to_str().unwrap(),
    ]);
    assert_files_eq(&all_stats, &explicit, "--all-metrics filter vs explicit --filter::stats");
}

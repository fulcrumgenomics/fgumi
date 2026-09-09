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

    assert_eq!(
        std::fs::read_to_string(&fused_histogram).expect("fused histogram exists"),
        std::fs::read_to_string(&standalone_histogram).expect("standalone histogram exists"),
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
#[test]
fn all_metrics_never_overwrites_an_explicit_per_stage_flag() {
    let dir = tempfile::tempdir().expect("tempdir");
    let header = create_minimal_header("chr1", 10_000);
    let bam_path = dir.path().join("in.bam");
    write_bam(&bam_path, &header, &[grouped_record("r1", 100, "AAAAAAAA")]);

    let explicit_path = dir.path().join("explicit.family_sizes.txt");
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
        "--group::family-size-histogram",
        explicit_path.to_str().unwrap(),
        "--all-metrics",
        all_metrics_prefix.to_str().unwrap(),
    ]);

    assert!(explicit_path.is_file(), "explicit flag's path must be used");
    assert!(
        !dir.path().join("all.group.family_size_histogram.txt").exists(),
        "the derived path must NOT also be written — explicit wins outright"
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

    assert!(dir.path().join("all.group.family_size_histogram.txt").is_file());
    assert!(dir.path().join("all.group.grouping_metrics.txt").is_file());
    assert!(
        dir.path().join("all.group.family_sizes.txt").is_file(),
        "group's own --metrics prefix output"
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

//! Fused-mode `--group::*` metrics parity: runall must produce the same
//! family-size-histogram/grouping-metrics output as standalone `fgumi group`
//! with the same flags, once the anti-goal block is removed.

use crate::helpers::bam_generator::{create_minimal_header, write_bam};
use fgumi_raw_bam::{RawRecord, SamBuilder};
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

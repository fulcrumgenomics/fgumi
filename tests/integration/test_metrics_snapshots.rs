//! Snapshot the `--metrics` TSV from runall.
//!
//! Any unintended change to per-stage record counts becomes a reviewable
//! diff. Update snapshots via `cargo insta review` after intentional
//! metric changes.

use std::process::Command;

use crate::helpers::grouped_bam_fixture::build_small_grouped_bam;

fn fgumi_bin() -> &'static str {
    env!("CARGO_BIN_EXE_fgumi")
}

#[test]
fn runall_group_to_filter_metrics_tsv() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input.bam");
    build_small_grouped_bam(&input);
    let output = tmp.path().join("output.bam");
    let metrics = tmp.path().join("metrics.tsv");

    let status = Command::new(fgumi_bin())
        .args(["runall", "--start-from", "group", "--threads", "1"])
        .args(["--filter::min-reads", "1"])
        .arg("--input")
        .arg(&input)
        .arg("--output")
        .arg(&output)
        .arg("--metrics")
        .arg(&metrics)
        .status()
        .expect("spawn fgumi");
    assert!(status.success(), "runall failed");

    let tsv = std::fs::read_to_string(&metrics).expect("read metrics TSV");
    let filtered = drop_wall_time_column(&tsv);

    insta::assert_snapshot!(filtered);
}

/// Drop the `wall_time_secs` column from every row so snapshots are stable.
fn drop_wall_time_column(tsv: &str) -> String {
    let mut lines = tsv.lines();
    let header = lines.next().expect("at least a header");
    let header_fields: Vec<&str> = header.split('\t').collect();
    let wall_col = header_fields
        .iter()
        .position(|&f| f == "wall_time_secs")
        .expect("header contains wall_time_secs");

    let drop = |line: &str| -> String {
        line.split('\t')
            .enumerate()
            .filter(|(i, _)| *i != wall_col)
            .map(|(_, f)| f)
            .collect::<Vec<_>>()
            .join("\t")
    };

    let mut out = String::new();
    out.push_str(&drop(header));
    out.push('\n');
    for line in lines {
        out.push_str(&drop(line));
        out.push('\n');
    }
    out
}

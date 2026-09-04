//! End-to-end determinism tests for the `retag` command on the chain builder.
//!
//! retag always routes through the declarative chain builder (the legacy serial
//! path is retired). These assert the output is independent of the worker count:
//! a no-`--threads` run (the chain at a single worker) produces identical output
//! records + normalized header and an identical `--metrics` TSV to a
//! `--threads N` run. `--threads 1` is exercised explicitly as a distinct
//! single-worker configuration. (Byte-parity against the pre-removal serial
//! binary lives in `test_retag_cutover_parity.rs`.)

use clap::Parser;
use fgumi_lib::commands::command::Command;
use fgumi_lib::commands::retag::Retag;
use fgumi_lib::sam::SamTag;
use fgumi_raw_bam::{RawRecord, SamBuilder};
use rstest::rstest;
use std::path::{Path, PathBuf};
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, write_bam};
use crate::helpers::read_bam_output;

/// A minimal mapped record with the given name and optional RX/BX tags. Sequence
/// and quals are fixed filler — retag only rewrites tags.
fn rec(name: &str, rx: Option<&str>, bx: Option<&str>) -> RawRecord {
    let mut b = SamBuilder::new();
    b.read_name(name.as_bytes())
        .sequence(b"ACGT")
        .qualities(&[30; 4])
        .flags(0)
        .ref_id(0)
        .pos(99)
        .mapq(60)
        .cigar_ops(&[4 << 4]);
    if let Some(v) = rx {
        b.add_string_tag(SamTag::RX, v.as_bytes());
    }
    if let Some(v) = bx {
        b.add_string_tag("BX".parse::<SamTag>().unwrap(), v.as_bytes());
    }
    b.build()
}

/// Write `records` to a BAM in `dir/<name>` and return its path.
fn write_input(dir: &Path, name: &str, records: &[RawRecord]) -> PathBuf {
    let path = dir.join(name);
    let header = create_minimal_header("chr1", 10_000);
    write_bam(&path, &header, records);
    path
}

/// A mixed input: records with RX, one with a pre-existing BX (to exercise
/// overwrite), and one with no RX (to exercise `src_missing`).
fn mixed_input(dir: &Path, count: usize) -> PathBuf {
    let mut records = Vec::with_capacity(count);
    for i in 0..count {
        match i % 4 {
            0 => records.push(rec(&format!("r{i}"), Some("ACGT"), None)),
            1 => records.push(rec(&format!("r{i}"), Some("TTTT"), Some("OLD"))),
            2 => records.push(rec(&format!("r{i}"), None, None)),
            _ => records.push(rec(&format!("r{i}"), Some("GGCC"), None)),
        }
    }
    write_input(dir, "in.bam", &records)
}

/// Parse + run `retag`. `threads = None` omits `--threads` (the chain at a
/// single worker); `Some(n)` passes `--threads n`. Both take the chain.
fn run_retag(
    input: &Path,
    output: &Path,
    threads: Option<usize>,
    ops: &[&str],
    metrics: Option<&Path>,
) -> anyhow::Result<()> {
    let mut args: Vec<String> = vec![
        "retag".into(),
        "-i".into(),
        input.display().to_string(),
        "-o".into(),
        output.display().to_string(),
    ];
    if let Some(t) = threads {
        args.push("--threads".into());
        args.push(t.to_string());
    }
    if let Some(m) = metrics {
        args.push("-M".into());
        args.push(m.display().to_string());
    }
    args.extend(ops.iter().map(|s| (*s).to_string()));
    let cmd = Retag::try_parse_from(args).expect("failed to parse retag args");
    cmd.execute("retag test")
}

// ── worker-count invariance: record + header ───────────────────────────────

#[rstest]
#[case::threads_1(1)]
#[case::threads_2(2)]
#[case::threads_4(4)]
fn output_matches_across_worker_counts(#[case] threads: usize) {
    let dir = TempDir::new().unwrap();
    let input = mixed_input(dir.path(), 40);
    // Fan RX out to BX and CB, then drop RX — exercises applied/overwritten/missing
    // across three positional ops.
    let ops = ["RX::copy::BX", "RX::copy::CB", "RX::delete"];

    let single_worker_out = dir.path().join("single_worker.bam");
    run_retag(&input, &single_worker_out, None, &ops, None).expect("single-worker run");

    let multi_worker_out = dir.path().join("multi_worker.bam");
    run_retag(&input, &multi_worker_out, Some(threads), &ops, None).expect("multi-worker run");

    let single_worker = read_bam_output(&single_worker_out);
    let multi_worker = read_bam_output(&multi_worker_out);
    assert_eq!(single_worker.1.len(), 40, "all records emitted");
    assert_eq!(
        multi_worker, single_worker,
        "--threads {threads} output must match the single-worker run"
    );
}

// ── --metrics TSV byte-parity ──────────────────────────────────────────────

#[rstest]
#[case::threads_1(1)]
#[case::threads_2(2)]
#[case::threads_4(4)]
fn metrics_match_across_worker_counts(#[case] threads: usize) {
    let dir = TempDir::new().unwrap();
    let input = mixed_input(dir.path(), 40);
    let ops = ["RX::copy::BX", "RX::move::CB"];

    let single_worker_out = dir.path().join("single_worker.bam");
    let single_worker_m = dir.path().join("single_worker.tsv");
    run_retag(&input, &single_worker_out, None, &ops, Some(&single_worker_m))
        .expect("single-worker run");

    let multi_worker_out = dir.path().join("multi_worker.bam");
    let multi_worker_m = dir.path().join("multi_worker.tsv");
    run_retag(&input, &multi_worker_out, Some(threads), &ops, Some(&multi_worker_m))
        .expect("multi-worker run");

    let single_worker_tsv =
        std::fs::read_to_string(&single_worker_m).expect("read single-worker tsv");
    let multi_worker_tsv = std::fs::read_to_string(&multi_worker_m).expect("read multi-worker tsv");
    assert_eq!(
        multi_worker_tsv, single_worker_tsv,
        "--metrics TSV must be byte-identical across worker counts"
    );
    // Non-vacuous: at least one op actually applied to some records.
    assert!(
        single_worker_tsv.lines().any(|l| {
            let cols: Vec<&str> = l.split('\t').collect();
            cols.len() >= 3 && cols[2].parse::<u64>().map(|n| n > 0).unwrap_or(false)
        }),
        "expected a non-zero records_applied somewhere in the metrics:\n{single_worker_tsv}"
    );
}

// ── warn-on-zero-match: a never-matching op reports records_applied == 0 ────

#[rstest]
#[case::threads_1(1)]
#[case::threads_2(2)]
fn zero_match_op_consistent_across_worker_counts(#[case] threads: usize) {
    let dir = TempDir::new().unwrap();
    let input = mixed_input(dir.path(), 20);
    // `ZZ` is present on no record, so this op matches zero records at any worker count.
    let ops = ["ZZ::delete"];

    let single_worker_out = dir.path().join("single_worker.bam");
    let single_worker_m = dir.path().join("single_worker.tsv");
    run_retag(&input, &single_worker_out, None, &ops, Some(&single_worker_m))
        .expect("single-worker run");

    let multi_worker_out = dir.path().join("multi_worker.bam");
    let multi_worker_m = dir.path().join("multi_worker.tsv");
    run_retag(&input, &multi_worker_out, Some(threads), &ops, Some(&multi_worker_m))
        .expect("multi-worker run");

    let single_worker_tsv = std::fs::read_to_string(&single_worker_m).unwrap();
    let multi_worker_tsv = std::fs::read_to_string(&multi_worker_m).unwrap();
    assert_eq!(
        multi_worker_tsv, single_worker_tsv,
        "zero-match metrics must match across worker counts"
    );
    // The durable observable of the warn: records_applied == 0 for the op.
    // TSV columns: operation, kind, records_applied, dst_overwritten, src_missing.
    let zz_row = single_worker_tsv
        .lines()
        .find(|l| l.starts_with("ZZ::delete\t"))
        .unwrap_or_else(|| panic!("metrics should carry the ZZ::delete row:\n{single_worker_tsv}"));
    let fields: Vec<&str> = zz_row.split('\t').collect();
    assert_eq!(
        fields.get(2).copied(),
        Some("0"),
        "records_applied must be 0 for a zero-match op, got row: {zz_row}"
    );
    // Output records are unchanged (no ZZ to delete), and the two BAMs agree.
    assert_eq!(read_bam_output(&multi_worker_out), read_bam_output(&single_worker_out));
}

// ── empty input: zero-record parity (degenerate boundary for the hooks) ─────

#[rstest]
#[case::threads_1(1)]
#[case::threads_4(4)]
fn empty_input_consistent_across_worker_counts(#[case] threads: usize) {
    let dir = TempDir::new().unwrap();
    // Header-only BAM: exercises the finalize hooks' zero-record path.
    let input = write_input(dir.path(), "empty.bam", &[]);
    let ops = ["RX::copy::BX", "ZZ::delete"];

    let single_worker_out = dir.path().join("single_worker.bam");
    let single_worker_m = dir.path().join("single_worker.tsv");
    run_retag(&input, &single_worker_out, None, &ops, Some(&single_worker_m))
        .expect("single-worker run");

    let multi_worker_out = dir.path().join("multi_worker.bam");
    let multi_worker_m = dir.path().join("multi_worker.tsv");
    run_retag(&input, &multi_worker_out, Some(threads), &ops, Some(&multi_worker_m))
        .expect("multi-worker run");

    // Header-only output, identical across worker counts.
    assert_eq!(read_bam_output(&multi_worker_out), read_bam_output(&single_worker_out));
    // Metrics TSVs byte-identical (all counts zero at any worker count).
    assert_eq!(
        std::fs::read_to_string(&multi_worker_m).unwrap(),
        std::fs::read_to_string(&single_worker_m).unwrap(),
        "empty-input metrics must match across worker counts"
    );
}

// ── multi-batch: enough records to cross pipeline batch boundaries ──────────

#[test]
fn multi_batch_consistent_across_worker_counts() {
    let dir = TempDir::new().unwrap();
    let input = mixed_input(dir.path(), 5_000);
    let ops = ["RX::copy::BX", "RX::delete"];

    let single_worker_out = dir.path().join("single_worker.bam");
    let single_worker_m = dir.path().join("single_worker.tsv");
    run_retag(&input, &single_worker_out, None, &ops, Some(&single_worker_m))
        .expect("single-worker run");

    let multi_worker_out = dir.path().join("multi_worker.bam");
    let multi_worker_m = dir.path().join("multi_worker.tsv");
    run_retag(&input, &multi_worker_out, Some(4), &ops, Some(&multi_worker_m))
        .expect("multi-worker run");

    assert_eq!(
        read_bam_output(&multi_worker_out),
        read_bam_output(&single_worker_out),
        "multi-batch multi-worker output must match the single-worker run"
    );
    assert_eq!(
        std::fs::read_to_string(&multi_worker_m).unwrap(),
        std::fs::read_to_string(&single_worker_m).unwrap(),
        "multi-batch summed metrics must match across worker counts"
    );
}

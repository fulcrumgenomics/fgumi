//! Determinism regression tests for the `group` and `dedup` commands.
//!
//! These tests guard against a class of non-determinism caused by iterating
//! `AHashMap`-keyed orientation subgroups (`(bool, bool) -> Vec<usize>`)
//! without sorting before assigning `MoleculeId`s. With Rust's per-process
//! random hasher, two runs of the same command on the same input would walk
//! orientation subgroups in different orders and produce identical molecule
//! groupings but mismatched `MI:Z` numbering across runs — which then
//! propagates into downstream consensus QNAMEs.
//!
//! The failing case requires multiple distinct orientations so that
//! `subgroups` has more than one entry; otherwise iteration order is trivially
//! stable.
//!
//! Both tests build a small paired-end BAM containing several UMI families
//! distributed across both FR (R1 forward, R2 reverse) and RF (R1 reverse,
//! R2 forward) orientations, run the command twice on the same input, and
//! assert that the per-record `MI:Z` tags are identical across runs.

use fgumi_lib::sam::SamTag;
use noodles::bam;
use noodles::sam::alignment::record_buf::data::field::Value;
use rstest::rstest;
use std::fs;
use std::path::Path;
use std::process::Command;
use tempfile::TempDir;

use crate::helpers::mi_determinism_fixture::build_mixed_orientation_bam;

/// Read all `(QNAME, flags, MI:Z)` tuples from a BAM in input order.
///
/// Returning the tuple lets the test compare both per-output-position (which
/// catches sort-order differences induced by re-numbered MIs) and
/// per-(QNAME, flags) (which catches MI re-numbering on the same record).
///
/// The MI tag is extracted via the typed `Value::String` accessor rather
/// than the `Debug` representation so the assertion compares the actual MI
/// string rather than coupling to `noodles`' `Value` Debug surface.
fn read_qname_mi(path: &Path) -> Vec<(String, u16, String)> {
    let mut reader = bam::io::Reader::new(fs::File::open(path).expect("open output BAM"));
    let header = reader.read_header().expect("read header");
    let mi = SamTag::MI.to_noodles_tag();
    reader
        .record_bufs(&header)
        .map(|r| {
            let r = r.expect("read record");
            let qname = r.name().expect("missing name").to_string();
            let flags = u16::from(r.flags());
            let mi_val = match r.data().get(&mi) {
                Some(Value::String(s)) => s.to_string(),
                Some(other) => panic!("expected MI tag to be a string, got {other:?}"),
                None => panic!("missing MI tag on output record {qname}"),
            };
            (qname, flags, mi_val)
        })
        .collect()
}

/// Run a fgumi subcommand `n_runs` times on the same input and assert that
/// the per-record MI tags are byte-identical across every pair of runs.
///
/// Each run is a fresh process invocation so `AHashMap`'s per-process random
/// hasher seed varies between runs; this is exactly the property that exposes
/// the orientation-subgroup iteration bug. Running more than two times and
/// comparing pairwise raises the probability of catching the variance even
/// when only a handful of subgroup keys are present.
///
/// `threads = None` omits the `--threads` flag entirely. For `fgumi group`
/// this triggers the dedicated `execute_single_threaded` fast path
/// (`src/lib/commands/group.rs`'s `if self.threading.threads.is_none()` at
/// the top of `execute`); `Some("1")` runs the parallel pipeline with one
/// worker, which is a different code path. Both must be deterministic.
fn assert_mi_deterministic(subcommand: &str, threads: Option<&str>, n_runs: usize) {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    build_mixed_orientation_bam(&input_bam);

    let mut all_runs: Vec<Vec<(String, u16, String)>> = Vec::with_capacity(n_runs);
    for run in 0..n_runs {
        let out = temp_dir.path().join(format!("out_{run}.bam"));
        let mut args: Vec<&str> = vec![
            subcommand,
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            out.to_str().unwrap(),
            "--strategy",
            "paired",
            "--compression-level",
            "1",
        ];
        if let Some(t) = threads {
            args.extend(["--threads", t]);
        }
        let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
            .args(&args)
            .status()
            .expect("Failed to run command");
        assert!(status.success(), "{subcommand} run {run} failed");
        all_runs.push(read_qname_mi(&out));
    }

    let threads_label = threads.unwrap_or("<none>");
    let baseline = &all_runs[0];
    let baseline_by_key: std::collections::HashMap<(String, u16), &str> =
        baseline.iter().map(|(q, f, m)| ((q.clone(), *f), m.as_str())).collect();

    for (run, run_records) in all_runs.iter().enumerate().skip(1) {
        assert_eq!(
            baseline.len(),
            run_records.len(),
            "{subcommand} produced different record counts across runs",
        );

        // Per-(QNAME, flags) comparison: catches MI tag re-numbering on the
        // same logical record, even if the output is re-sorted by MI.
        let mut keyed_mismatches = 0usize;
        let mut sample_diffs: Vec<String> = Vec::new();
        for (qname, flags, mi_run) in run_records {
            let mi_base = baseline_by_key.get(&(qname.clone(), *flags));
            let Some(mi_base) = mi_base else {
                continue;
            };
            if *mi_base != mi_run.as_str() {
                keyed_mismatches += 1;
                if sample_diffs.len() < 5 {
                    sample_diffs
                        .push(format!("{qname} flags={flags}: run0={mi_base} run{run}={mi_run}",));
                }
            }
        }

        // Positional comparison: catches re-sorting that occurs because
        // templates are sorted by MI on the way out.
        let pos_mismatches =
            baseline.iter().zip(run_records.iter()).filter(|(x, y)| x != y).count();

        assert_eq!(
            pos_mismatches, 0,
            "{subcommand} --threads {threads_label}: output order differs between run 0 and run {run} ({pos_mismatches} positional mismatches)",
        );

        assert_eq!(
            keyed_mismatches,
            0,
            "{subcommand} --threads {threads_label}: per-(QNAME,flags) MI tags differ between run 0 and run {run} ({keyed_mismatches}/{} records). Positional mismatches: {pos_mismatches}. Examples: {sample_diffs:?}",
            baseline.len(),
        );
    }
}

/// Repeated runs of `fgumi group` and `fgumi dedup` with `--strategy paired`
/// must assign identical `MI:Z` tags to every record across runs and across
/// every supported threading mode.
///
/// The three threads cases exercise three distinct code paths:
///
/// * `None` — `group` only: the dedicated `execute_single_threaded` fast
///   path (taken when `--threads` is omitted), which assigns `MoleculeId`s
///   inline without going through the unified pipeline scheduler.
/// * `Some("1")` — the unified pipeline with a single worker.
/// * `Some("4")` — the unified pipeline with four workers, which was the
///   originally failing case before the serial-ordered `MI Assign` hook.
#[rstest]
#[case::group_no_threads("group", None)]
#[case::group_threads_1("group", Some("1"))]
#[case::group_threads_4("group", Some("4"))]
#[case::dedup_threads_1("dedup", Some("1"))]
#[case::dedup_threads_4("dedup", Some("4"))]
fn test_paired_mi_deterministic(#[case] subcommand: &str, #[case] threads: Option<&str>) {
    assert_mi_deterministic(subcommand, threads, 6);
}

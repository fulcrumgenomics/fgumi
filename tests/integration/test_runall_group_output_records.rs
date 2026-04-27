//! Regression: `runall --start-from group --threads 1` must produce a
//! non-empty output BAM whose record count matches the simplex-consensus
//! result from a sequential `fgumi simplex` on the same input.
//!
//! Filed in response to a real-data report at HEAD `4aebe2f` where the
//! pipeline-internal Consensus stage emitted records but the writer
//! produced zero records on disk. The pre-existing
//! `test_metrics_snapshots::runall_group_to_filter_metrics_tsv` test
//! caught the symptom (Filter `0/0` in the metrics TSV) but only asserted
//! on the metrics file, not the output BAM, so the regression slipped
//! through.
//!
//! Today's failure mode: pipeline produces 0 records. Expected: 3
//! consensus reads, one per UMI family in `build_small_grouped_bam`.

use std::process::Command;

use crate::helpers::grouped_bam_fixture::build_small_grouped_bam;
use crate::helpers::mi_determinism_fixture::{build_duplex_pair_bam, build_mixed_orientation_bam};

fn fgumi_bin() -> &'static str {
    env!("CARGO_BIN_EXE_fgumi")
}

fn count_bam_records(path: &std::path::Path) -> usize {
    use noodles::bam;
    use std::fs;
    let mut reader = bam::io::Reader::new(fs::File::open(path).expect("open output bam"));
    let _ = reader.read_header().expect("header");
    reader.records().count()
}

#[test]
fn runall_group_to_filter_writes_consensus_records_threads_1() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input.bam");
    build_small_grouped_bam(&input);
    let output = tmp.path().join("output.bam");

    let status = Command::new(fgumi_bin())
        .args(["runall", "--start-from", "group", "--threads", "1"])
        .args(["--filter::min-reads", "1"])
        .arg("--input")
        .arg(&input)
        .arg("--output")
        .arg(&output)
        .status()
        .expect("spawn fgumi runall");
    assert!(status.success(), "fgumi runall did not exit cleanly");

    let n = count_bam_records(&output);
    assert_eq!(
        n, 3,
        "runall --start-from group on the 3-UMI-family fixture should emit 3 \
         simplex consensus records; got {n}",
    );
}

/// Same as the threads=4 duplex test but with `--threads 1` to isolate
/// any parallelism-related interactions with `MiAssignStage`'s segment
/// boundary.
#[test]
fn runall_group_to_filter_writes_consensus_records_threads_1_duplex() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input.bam");
    build_duplex_pair_bam(&input);
    let output = tmp.path().join("output.bam");

    let status = Command::new(fgumi_bin())
        .args([
            "runall",
            "--start-from",
            "group",
            "--threads",
            "1",
            "--consensus",
            "duplex",
            "--group::strategy",
            "paired",
            "--group::max-edits",
            "1",
        ])
        .args(["--filter::min-reads", "1"])
        .arg("--input")
        .arg(&input)
        .arg("--output")
        .arg(&output)
        .status()
        .expect("spawn fgumi runall");
    assert!(status.success(), "fgumi runall did not exit cleanly");

    let n = count_bam_records(&output);
    assert!(n > 0, "single-threaded duplex runall must emit records (got 0)",);
}

/// Diagnostic: same as the failing duplex test but with maximally permissive
/// filter thresholds. If THIS still produces zero records, the bug is in
/// the consensus→filter wiring (records never reach filter) and not in the
/// filter's predicate.
#[test]
fn runall_group_to_filter_writes_consensus_records_threads_4_duplex_permissive() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input.bam");
    build_duplex_pair_bam(&input);
    let output = tmp.path().join("output.bam");

    let status = Command::new(fgumi_bin())
        .args([
            "runall",
            "--start-from",
            "group",
            "--threads",
            "4",
            "--consensus",
            "duplex",
            "--group::strategy",
            "paired",
            "--group::max-edits",
            "1",
            "--filter::min-reads",
            "1",
            "--filter::max-read-error-rate",
            "1.0",
            "--filter::max-base-error-rate",
            "1.0",
            "--filter::max-no-call-fraction",
            "1.0",
            "--filter::min-base-quality",
            "0",
            "--filter::filter-by-template",
            "false",
        ])
        .arg("--input")
        .arg(&input)
        .arg("--output")
        .arg(&output)
        .status()
        .expect("spawn fgumi runall");
    assert!(status.success(), "fgumi runall did not exit cleanly");

    let n = count_bam_records(&output);
    assert!(
        n > 0,
        "permissive duplex runall must emit records (got 0 — wiring bug, not filter predicate)",
    );
}

/// Mirrors the user-reported failure: `--threads 4 --consensus duplex` on a
/// paired-UMI fixture. The pre-existing `test_runall_mi_determinism` tests
/// use `--stop-after group`, which never reaches Consensus or Filter; this
/// test pushes records through the full `Consensus -> Filter -> sink` chain.
#[test]
fn runall_group_to_filter_writes_consensus_records_threads_4_duplex() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input.bam");
    build_duplex_pair_bam(&input);
    let output = tmp.path().join("output.bam");

    let status = Command::new(fgumi_bin())
        .args([
            "runall",
            "--start-from",
            "group",
            "--threads",
            "4",
            "--consensus",
            "duplex",
            "--group::strategy",
            "paired",
            "--group::max-edits",
            "1",
        ])
        .args(["--filter::min-reads", "1"])
        .arg("--input")
        .arg(&input)
        .arg("--output")
        .arg(&output)
        .status()
        .expect("spawn fgumi runall");
    assert!(status.success(), "fgumi runall did not exit cleanly");

    let n = count_bam_records(&output);
    assert!(
        n > 0,
        "runall --start-from group --consensus duplex must emit at least one \
         duplex consensus record; got 0 (regression — pipeline silently \
         dropped every record)",
    );
}

/// Same as the duplex case but for simplex consensus.
#[test]
fn runall_group_to_filter_writes_consensus_records_threads_4() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input.bam");
    build_mixed_orientation_bam(&input);
    let output = tmp.path().join("output.bam");

    let status = Command::new(fgumi_bin())
        .args([
            "runall",
            "--start-from",
            "group",
            "--threads",
            "4",
            "--group::strategy",
            "paired",
            "--group::max-edits",
            "1",
        ])
        .args(["--filter::min-reads", "1"])
        .arg("--input")
        .arg(&input)
        .arg("--output")
        .arg(&output)
        .status()
        .expect("spawn fgumi runall");
    assert!(status.success(), "fgumi runall did not exit cleanly");

    let n = count_bam_records(&output);
    // The fixture has 40 position groups × 256 paired UMIs × 2 orientations =
    // 20480 templates spread evenly across 40 positions. Simplex consensus
    // collapses each MI group to one consensus read; with `--filter::min-reads 1`
    // the filter accepts every consensus read. Exact count depends on how the
    // paired-UMI assigner clusters orientations per position; what this test
    // asserts is the regression-class invariant — runall must not silently emit
    // an empty output BAM.
    assert!(
        n > 0,
        "runall --start-from group on a 40-position-group fixture must emit \
         at least one consensus record; got 0 (regression — pipeline silently \
         dropped every record)",
    );
}

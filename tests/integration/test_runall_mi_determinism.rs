//! Determinism regression test for the `fgumi runall` pipeline's group stage.
//!
//! Exercises the same non-determinism bug PR #319 addressed for standalone
//! `group`/`dedup`, but for runall's own engine in
//! `src/lib/runall/engine/stages/group_assign.rs`. The bug: workers call a
//! shared `Arc<dyn UmiAssigner>` whose internal `AtomicU64` counter increments
//! in arrival order, so cross-run `MI:Z` integers vary even though groupings
//! are stable.
//!
//! The test uses `--start-from sort --stop-after group` to capture the raw
//! grouped BAM (with integer MI:Z tags) before consensus collapses them to UMI
//! strings. Without `--stop-after group`, the non-determinism in MI integers is
//! invisible in the consensus output because consensus QNAME and MI tags are
//! derived from the UMI string, not the integer.

use std::process::Command;
use tempfile::TempDir;

use crate::helpers::mi_determinism_fixture::{
    assert_mi_deterministic_across_runs, build_mixed_orientation_bam, read_qname_mi,
};

/// Run `fgumi runall --start-from sort --stop-after group` `n_runs` times on
/// the same input BAM and assert that the per-`(QNAME, flags)` `MI:Z` integers
/// are byte-identical across every pair of runs.
///
/// Each run is a fresh process invocation so `AHashMap`'s per-process random
/// hasher seed varies between runs; this is exactly the property that exposes
/// the orientation-subgroup iteration bug in `GroupAssignStage`.
fn assert_runall_group_mi_deterministic(threads: &str, n_runs: usize) {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    build_mixed_orientation_bam(&input_bam);

    let mut all_runs: Vec<Vec<(String, u16, String)>> = Vec::with_capacity(n_runs);
    for run in 0..n_runs {
        let out = temp_dir.path().join(format!("out_{run}.bam"));
        let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
            .args([
                "runall",
                "--start-from",
                "sort",
                "--stop-after",
                "group",
                "--input",
                input_bam.to_str().unwrap(),
                "--output",
                out.to_str().unwrap(),
                "--group::strategy",
                "paired",
                "--group::max-edits",
                "1",
                "--group::min-mapq",
                "0",
                "--threads",
                threads,
                "--compression-level",
                "1",
            ])
            .status()
            .expect("Failed to run runall");
        assert!(status.success(), "runall run {run} failed");
        all_runs.push(read_qname_mi(&out));
    }

    let label = format!("runall group --threads {threads}");
    assert_mi_deterministic_across_runs(&all_runs, n_runs, &label);
}

/// Repeated runs of `fgumi runall --start-from sort --stop-after group
/// --group::strategy paired --threads 4` must assign identical `MI:Z` integer
/// tags to every record.
///
/// This test is expected to **fail** on code with the `GroupAssignStage`
/// orientation-subgroup non-determinism bug (pre-fix).
#[test]
fn test_runall_group_paired_mi_deterministic_threads_4() {
    assert_runall_group_mi_deterministic("4", 6);
}

/// Single-threaded path is deterministic by construction — there is no
/// inter-worker counter race, so MI integers are assigned in a single
/// sequential pass. This test pins that property as a regression guard.
#[test]
fn test_runall_group_paired_mi_deterministic_threads_1() {
    assert_runall_group_mi_deterministic("1", 6);
}

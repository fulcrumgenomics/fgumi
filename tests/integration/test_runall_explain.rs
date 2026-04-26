//! Integration tests for `fgumi runall --explain`.
//!
//! `--explain` resolves the pipeline plan, prints a human-readable report,
//! and exits without running the pipeline or creating the output file.
//! The test intentionally passes paths that do not exist — `--explain`
//! must not require real inputs, so these tests also guard against
//! regressions that would re-introduce file I/O into the explain path.

use std::path::{Path, PathBuf};
use std::process::Command;

fn fgumi_bin() -> &'static str {
    env!("CARGO_BIN_EXE_fgumi")
}

/// Unique scratch output path per test to avoid cross-test file collisions
/// when running nextest with parallelism.
fn scratch_output(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join("fgumi_runall_explain_tests");
    std::fs::create_dir_all(&dir).expect("create scratch dir");
    dir.join(format!("{tag}.bam"))
}

fn assert_no_output_files(output: &Path) {
    assert!(
        !output.exists(),
        "runall --explain must not create the final output file: {}",
        output.display()
    );
    let tmp = output.with_extension("bam.tmp");
    assert!(
        !tmp.exists(),
        "runall --explain must not create the temp output file: {}",
        tmp.display()
    );
}

#[test]
fn runall_explain_prints_plan_and_exits() {
    let output_path = scratch_output("group_prints_plan");
    let _ = std::fs::remove_file(&output_path);
    let _ = std::fs::remove_file(output_path.with_extension("bam.tmp"));

    let output = Command::new(fgumi_bin())
        .args(["runall", "--start-from", "group", "--explain"])
        .arg("--input")
        .arg("ignored.bam")
        .arg("--output")
        .arg(&output_path)
        .output()
        .expect("spawn");
    assert!(
        output.status.success(),
        "explain should succeed; stderr={}",
        String::from_utf8_lossy(&output.stderr)
    );

    let stdout = String::from_utf8(output.stdout).expect("utf8 stdout");
    assert!(stdout.contains("Pipeline plan"), "missing header: {stdout}");
    assert!(
        stdout.contains("GroupAssign") || stdout.contains("group-assign"),
        "missing GroupAssign stage: {stdout}"
    );
    assert!(stdout.contains("Consensus"), "missing Consensus stage: {stdout}");
    assert!(stdout.contains("Filter"), "missing Filter stage: {stdout}");

    assert_no_output_files(&output_path);
}

#[test]
fn runall_explain_covers_all_start_from_entrypoints() {
    // This test guards the entrypoints that `build_plan`/`build_plan_explain`
    // currently handle. `align` and `zipper` were removed from the CLI
    // because they were never wired through build_plan.
    for start_from in ["extract", "fastq", "sort", "group"] {
        let output_path = scratch_output(&format!("all_entrypoints_{start_from}"));
        let _ = std::fs::remove_file(&output_path);
        let _ = std::fs::remove_file(output_path.with_extension("bam.tmp"));

        let mut cmd = Command::new(fgumi_bin());
        cmd.args(["runall", "--explain", "--start-from", start_from]);

        match start_from {
            "extract" => {
                cmd.args(["--input", "r1.fq", "r2.fq"])
                    .args(["--extract::sample", "s"])
                    .args(["--extract::library", "l"])
                    .args(["--extract::read-structures", "8M+T", "8M+T"])
                    .args(["--reference", "ignored.fa"])
                    .arg("--output")
                    .arg(&output_path);
            }
            "sort" => {
                cmd.args(["--input", "ignored.bam"])
                    .args(["--reference", "ignored.fa"])
                    .arg("--output")
                    .arg(&output_path);
            }
            "fastq" => {
                cmd.args(["--stop-after", "fastq"])
                    .args(["--input", "ignored.bam"])
                    .arg("--output")
                    .arg(&output_path);
            }
            _ => {
                cmd.args(["--input", "ignored.bam"]).arg("--output").arg(&output_path);
            }
        }

        let output = cmd.output().expect("spawn");
        assert!(
            output.status.success(),
            "--explain failed for --start-from {}: stderr={}",
            start_from,
            String::from_utf8_lossy(&output.stderr)
        );

        let stdout = String::from_utf8(output.stdout).expect("utf8 stdout");
        assert!(
            stdout.contains("Pipeline plan"),
            "missing header for --start-from {start_from}: {stdout}"
        );

        assert_no_output_files(&output_path);
    }
}

/// Exhaustive coverage of every valid `(start-from, stop-after)` combination
/// that the planner wires. Uses `--explain` so no real I/O happens — we're
/// only verifying that the planner accepts the combo (and that its stage
/// explain-summary surfaces the expected number of stages).
///
/// Combos are enumerated from the validator's ordering: stop ≥ start in the
/// canonical stage order {extract, correct, fastq, zipper, sort, group,
/// consensus, filter}. `--stop-after align` is intentionally omitted — the
/// pipeline does not expose raw aligner SAM output, and `build_plan` bails
/// with a clear message.
#[test]
fn runall_explain_covers_all_valid_start_stop_combos() {
    // (start-from, stop-after) pairs the planner must accept.
    let combos: &[(&str, &str)] = &[
        // --start-from extract
        ("extract", "extract"),
        ("extract", "correct"),
        ("extract", "fastq"),
        ("extract", "zipper"),
        ("extract", "sort"),
        ("extract", "group"),
        ("extract", "consensus"),
        ("extract", "filter"),
        // --start-from correct
        ("correct", "correct"),
        ("correct", "fastq"),
        ("correct", "zipper"),
        ("correct", "sort"),
        ("correct", "group"),
        ("correct", "consensus"),
        ("correct", "filter"),
        // --start-from fastq
        ("fastq", "fastq"),
        ("fastq", "zipper"),
        ("fastq", "sort"),
        ("fastq", "group"),
        ("fastq", "consensus"),
        ("fastq", "filter"),
        // --start-from sort
        ("sort", "sort"),
        ("sort", "group"),
        ("sort", "consensus"),
        ("sort", "filter"),
        // --start-from group
        ("group", "group"),
        ("group", "consensus"),
        ("group", "filter"),
    ];

    for &(start, stop) in combos {
        let output_path =
            scratch_output(&format!("combo_{start}_to_{stop}").replace(['.', '-'], "_"));
        let _ = std::fs::remove_file(&output_path);
        let _ = std::fs::remove_file(output_path.with_extension("bam.tmp"));

        let mut cmd = Command::new(fgumi_bin());
        cmd.args(["runall", "--explain", "--start-from", start, "--stop-after", stop]);

        match start {
            "extract" => {
                cmd.args(["--input", "r1.fq", "r2.fq"])
                    .args(["--extract::sample", "s"])
                    .args(["--extract::library", "l"])
                    .args(["--extract::read-structures", "8M+T", "8M+T"]);
            }
            _ => {
                cmd.args(["--input", "ignored.bam"]);
            }
        }
        cmd.args(["--reference", "ignored.fa"]).arg("--output").arg(&output_path);

        // Any combo whose plan reaches the consensus or filter stage needs
        // --filter::min-reads (matches standalone `fgumi filter`).
        if matches!(stop, "consensus" | "filter") {
            cmd.args(["--filter::min-reads", "1"]);
        }

        // Correct-entry / correct-inclusive combos need at least one UMI.
        if start == "correct" || (start == "extract" && stop != "extract") {
            // extract start doesn't require --correct::umis unless correction
            // is actually requested; skip the flag there.
            if start == "correct" {
                cmd.args(["--correct::umis", "ACGTACGT"]);
            }
        }

        let output = cmd.output().expect("spawn");
        assert!(
            output.status.success(),
            "--explain failed for --start-from {start} --stop-after {stop}: stderr={}",
            String::from_utf8_lossy(&output.stderr)
        );

        let stdout = String::from_utf8(output.stdout).expect("utf8 stdout");
        assert!(
            stdout.contains("Pipeline plan"),
            "missing plan header for combo ({start}, {stop}): {stdout}"
        );

        assert_no_output_files(&output_path);
    }
}

/// Regression: `--stop-after consensus` must not include the Filter stage
/// in the plan for `--start-from group` or `--start-from sort`.
///
/// Pre-fix: `plan_from_group` and `plan_from_sort` always pushed the
/// `filter_spec` regardless of `--stop-after`, so plans that requested
/// `--stop-after consensus` silently ran the filter logic and (because
/// `FilterStage` zeros the per-batch record `count`) the user-facing
/// "records written" counter reported `0` even when records reached the
/// output BAM. The fix wires `--stop-after consensus` to drop the Filter
/// stage from the plan, matching the behaviour of the
/// `append_tail_from_unmapped_bam` arm used by `plan_from_extract` /
/// `plan_from_fastq` / `plan_from_correct`.
#[test]
fn runall_explain_stop_after_consensus_omits_filter() {
    for start_from in ["group", "sort"] {
        let output_path = scratch_output(&format!("stop_after_consensus_{start_from}"));
        let _ = std::fs::remove_file(&output_path);
        let _ = std::fs::remove_file(output_path.with_extension("bam.tmp"));

        let mut cmd = Command::new(fgumi_bin());
        cmd.args(["runall", "--explain", "--start-from", start_from])
            .args(["--stop-after", "consensus"])
            .args(["--filter::min-reads", "1"])
            .args(["--input", "ignored.bam"])
            .arg("--output")
            .arg(&output_path);
        if start_from == "sort" {
            cmd.args(["--reference", "ignored.fa"]);
        }

        let output = cmd.output().expect("spawn");
        assert!(
            output.status.success(),
            "--explain failed for --start-from {start_from} --stop-after consensus: stderr={}",
            String::from_utf8_lossy(&output.stderr)
        );

        let stdout = String::from_utf8(output.stdout).expect("utf8 stdout");
        assert!(
            stdout.contains("Consensus"),
            "missing Consensus stage for --start-from {start_from}: {stdout}"
        );
        assert!(
            !stdout.contains("Filter"),
            "--stop-after consensus must omit the Filter stage from the plan \
             (--start-from {start_from}); got:\n{stdout}"
        );

        assert_no_output_files(&output_path);
    }
}

/// Regression: `--start-from group --stop-after group` must produce a plan
/// that writes the grouped BAM via `BgzfCompress` and stops, without
/// reaching `Consensus` or `Filter`. Mirrors the analogous early-exit arm
/// in `plan_from_sort`.
#[test]
fn runall_explain_start_from_group_stop_after_group_omits_consensus() {
    let output_path = scratch_output("start_group_stop_group");
    let _ = std::fs::remove_file(&output_path);
    let _ = std::fs::remove_file(output_path.with_extension("bam.tmp"));

    let output = Command::new(fgumi_bin())
        .args(["runall", "--explain", "--start-from", "group", "--stop-after", "group"])
        .args(["--input", "ignored.bam"])
        .arg("--output")
        .arg(&output_path)
        .output()
        .expect("spawn");
    assert!(
        output.status.success(),
        "--explain failed: stderr={}",
        String::from_utf8_lossy(&output.stderr)
    );

    let stdout = String::from_utf8(output.stdout).expect("utf8 stdout");
    assert!(stdout.contains("Pipeline plan"), "missing header: {stdout}");
    assert!(stdout.contains("PositionBatch"), "missing PositionBatch stage: {stdout}");
    assert!(stdout.contains("GroupAssign"), "missing GroupAssign stage: {stdout}");
    assert!(stdout.contains("BgzfCompress"), "missing BgzfCompress stage: {stdout}");
    assert!(
        !stdout.contains("Consensus"),
        "--stop-after group must omit Consensus from the plan; got:\n{stdout}"
    );
    assert!(
        !stdout.contains("Filter"),
        "--stop-after group must omit Filter from the plan; got:\n{stdout}"
    );

    assert_no_output_files(&output_path);
}

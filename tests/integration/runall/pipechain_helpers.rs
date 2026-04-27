//! Shared helpers for the runall pipechain combinator tests in
//! `test_pipechain_combinations.rs`.
//!
//! Each test in that file picks a `(start_from, stop_after)` pair, builds the
//! corresponding standalone chain (v1) and the equivalent
//! `runall --start-from S --stop-after T` invocation (v2), and compares the
//! outputs via `fgumi compare bams --command <preset>` (or byte-equality for
//! FASTQ output).
//!
//! Fixtures are produced via `fgumi simulate` subcommands (gated behind the
//! `simulate` feature). Comparison uses `fgumi compare bams` (gated behind the
//! `compare` feature). The whole module is therefore feature-gated behind
//! both via `runall/mod.rs`.

#![allow(dead_code, reason = "helpers used selectively across test functions")]

use std::path::{Path, PathBuf};
use std::process::Command;

/// Path to the `fgumi` binary built by the integration-test harness.
pub(crate) fn fgumi_bin() -> &'static str {
    env!("CARGO_BIN_EXE_fgumi")
}

// ---------------------------------------------------------------------------
// fgumi simulate wrappers
// ---------------------------------------------------------------------------

/// Run `fgumi simulate fastq-reads` to produce paired-end gzipped FASTQs.
///
/// Returns `(R1, R2)` paths.  Uses a fixed seed and small molecule count to
/// keep the fixture deterministic and the tests fast.
pub(crate) fn simulate_fastq_reads(
    out_dir: &Path,
    num_molecules: usize,
    reference: &Path,
) -> (PathBuf, PathBuf) {
    let r1 = out_dir.join("sim_R1.fastq.gz");
    let r2 = out_dir.join("sim_R2.fastq.gz");
    let truth = out_dir.join("sim_truth.tsv");
    let output = Command::new(fgumi_bin())
        .args([
            "simulate",
            "fastq-reads",
            "--r1",
            r1.to_str().unwrap(),
            "--r2",
            r2.to_str().unwrap(),
            "--truth",
            truth.to_str().unwrap(),
            "--reference",
            reference.to_str().unwrap(),
            "--num-molecules",
            &num_molecules.to_string(),
            "--seed",
            "42",
            "--threads",
            "1",
        ])
        .output()
        .expect("spawn fgumi simulate fastq-reads");
    assert!(
        output.status.success(),
        "simulate fastq-reads failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    (r1, r2)
}

/// Run `fgumi simulate correct-reads` to produce an unmapped BAM with `RX`
/// tags plus a UMI includelist file.
pub(crate) fn simulate_correct_reads(out_dir: &Path, num_reads: usize) -> (PathBuf, PathBuf) {
    let bam = out_dir.join("sim_correct.bam");
    let includelist = out_dir.join("sim_correct_umis.txt");
    let truth = out_dir.join("sim_correct_truth.tsv");
    let output = Command::new(fgumi_bin())
        .args([
            "simulate",
            "correct-reads",
            "--output",
            bam.to_str().unwrap(),
            "--includelist",
            includelist.to_str().unwrap(),
            "--truth",
            truth.to_str().unwrap(),
            "--num-reads",
            &num_reads.to_string(),
            "--seed",
            "42",
            "--threads",
            "1",
            "--compression-level",
            "1",
        ])
        .output()
        .expect("spawn fgumi simulate correct-reads");
    assert!(
        output.status.success(),
        "simulate correct-reads failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    (bam, includelist)
}

/// Run `fgumi simulate mapped-reads` to produce a template-coordinate-sorted
/// mapped BAM with `RX` tags.
pub(crate) fn simulate_mapped_reads(
    out_dir: &Path,
    num_molecules: usize,
    reference: &Path,
) -> PathBuf {
    let bam = out_dir.join("sim_mapped.bam");
    let truth = out_dir.join("sim_mapped_truth.tsv");
    let output = Command::new(fgumi_bin())
        .args([
            "simulate",
            "mapped-reads",
            "--output",
            bam.to_str().unwrap(),
            "--truth",
            truth.to_str().unwrap(),
            "--reference",
            reference.to_str().unwrap(),
            "--num-molecules",
            &num_molecules.to_string(),
            "--seed",
            "42",
            "--threads",
            "1",
            "--compression-level",
            "1",
        ])
        .output()
        .expect("spawn fgumi simulate mapped-reads");
    assert!(
        output.status.success(),
        "simulate mapped-reads failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    bam
}

/// Run `fgumi simulate grouped-reads` (simplex / non-duplex) to produce a
/// template-coordinate-sorted BAM with `MI` tags assigned.
pub(crate) fn simulate_grouped_reads(
    out_dir: &Path,
    num_molecules: usize,
    reference: &Path,
) -> PathBuf {
    let bam = out_dir.join("sim_grouped.bam");
    let truth = out_dir.join("sim_grouped_truth.tsv");
    let output = Command::new(fgumi_bin())
        .args([
            "simulate",
            "grouped-reads",
            "--output",
            bam.to_str().unwrap(),
            "--truth",
            truth.to_str().unwrap(),
            "--reference",
            reference.to_str().unwrap(),
            "--num-molecules",
            &num_molecules.to_string(),
            "--seed",
            "42",
            "--threads",
            "1",
            "--compression-level",
            "1",
        ])
        .output()
        .expect("spawn fgumi simulate grouped-reads");
    assert!(
        output.status.success(),
        "simulate grouped-reads failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    bam
}

// ---------------------------------------------------------------------------
// fgumi standalone-stage wrappers
// ---------------------------------------------------------------------------

/// Run `fgumi extract` on paired FASTQ inputs to produce an unmapped BAM with
/// UMIs in the `RX` tag.
pub(crate) fn run_extract(r1: &Path, r2: &Path, out: &Path) {
    let output = Command::new(fgumi_bin())
        .args([
            "extract",
            "--inputs",
            r1.to_str().unwrap(),
            r2.to_str().unwrap(),
            "--output",
            out.to_str().unwrap(),
            "--sample",
            "S",
            "--library",
            "L",
            "--read-structures",
            "8M+T",
            "+T",
            "--compression-level",
            "1",
        ])
        .output()
        .expect("spawn fgumi extract");
    assert!(
        output.status.success(),
        "fgumi extract failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

/// Run `fgumi correct` on an unmapped BAM and includelist.  No rejects/metrics
/// outputs (kept simple for parity tests).
pub(crate) fn run_correct(input: &Path, includelist: &Path, out: &Path) {
    let output = Command::new(fgumi_bin())
        .args([
            "correct",
            "--input",
            input.to_str().unwrap(),
            "--output",
            out.to_str().unwrap(),
            "--umi-files",
            includelist.to_str().unwrap(),
            "--max-mismatches",
            "2",
            "--min-distance",
            "2",
            "--compression-level",
            "1",
            "--threads",
            "1",
        ])
        .output()
        .expect("spawn fgumi correct");
    assert!(
        output.status.success(),
        "fgumi correct failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

/// Run `fgumi fastq` on an unmapped BAM, capturing stdout to `out`.
pub(crate) fn run_fastq(input: &Path, out: &Path) {
    let output = Command::new(fgumi_bin())
        .args(["fastq", "--input", input.to_str().unwrap()])
        .output()
        .expect("spawn fgumi fastq");
    assert!(
        output.status.success(),
        "fgumi fastq failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    std::fs::write(out, &output.stdout).expect("write fgumi fastq stdout");
}

/// Run `bwa mem -p -K 150000000 -t 1` on an interleaved FASTQ and capture
/// stdout to `out_sam`.
///
/// Uses `-K 150000000` to match runall's default `--aligner::chunk-size` and
/// `-t 1` for determinism.
pub(crate) fn run_bwa_mem(reference: &Path, fastq: &Path, out_sam: &Path) {
    let output = Command::new("bwa")
        .args([
            "mem",
            "-p",
            "-K",
            "150000000",
            "-t",
            "1",
            reference.to_str().unwrap(),
            fastq.to_str().unwrap(),
        ])
        .output()
        .expect("spawn bwa mem");
    assert!(
        output.status.success(),
        "bwa mem failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    std::fs::write(out_sam, &output.stdout).expect("write bwa mem stdout");
}

/// Run `fgumi zipper` to merge an unmapped BAM with a mapped SAM.
pub(crate) fn run_zipper(unmapped: &Path, mapped_sam: &Path, reference: &Path, out: &Path) {
    let output = Command::new(fgumi_bin())
        .args([
            "zipper",
            "--unmapped",
            unmapped.to_str().unwrap(),
            "--input",
            mapped_sam.to_str().unwrap(),
            "--output",
            out.to_str().unwrap(),
            "--reference",
            reference.to_str().unwrap(),
            "--compression-level",
            "1",
        ])
        .output()
        .expect("spawn fgumi zipper");
    assert!(
        output.status.success(),
        "fgumi zipper failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

/// Run `fgumi sort --order template-coordinate` on a BAM.
pub(crate) fn run_sort(input: &Path, out: &Path) {
    let output = Command::new(fgumi_bin())
        .args([
            "sort",
            "--input",
            input.to_str().unwrap(),
            "--output",
            out.to_str().unwrap(),
            "--order",
            "template-coordinate",
            "--threads",
            "1",
        ])
        .output()
        .expect("spawn fgumi sort");
    assert!(
        output.status.success(),
        "fgumi sort failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

/// Run `fgumi group --strategy <strategy>` on a template-coordinate-sorted BAM.
pub(crate) fn run_group(input: &Path, strategy: &str, out: &Path) {
    let output = Command::new(fgumi_bin())
        .args([
            "group",
            "--input",
            input.to_str().unwrap(),
            "--output",
            out.to_str().unwrap(),
            "--strategy",
            strategy,
        ])
        .output()
        .expect("spawn fgumi group");
    assert!(
        output.status.success(),
        "fgumi group failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

/// Run `fgumi simplex` on a grouped BAM with `--min-reads`.
pub(crate) fn run_simplex_consensus(input: &Path, out: &Path, min_reads: usize) {
    let output = Command::new(fgumi_bin())
        .args([
            "simplex",
            "--input",
            input.to_str().unwrap(),
            "--output",
            out.to_str().unwrap(),
            "--min-reads",
            &min_reads.to_string(),
            "--compression-level",
            "1",
        ])
        .output()
        .expect("spawn fgumi simplex");
    assert!(
        output.status.success(),
        "fgumi simplex failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

/// Run `fgumi filter` with the canonical simplex-stable parameter set used
/// across these pipechain tests.
pub(crate) fn run_filter(input: &Path, out: &Path) {
    let output = Command::new(fgumi_bin())
        .args([
            "filter",
            "--input",
            input.to_str().unwrap(),
            "--output",
            out.to_str().unwrap(),
            "--min-base-quality",
            "40",
            "--max-read-error-rate",
            "0.025",
            "--max-base-error-rate",
            "0.1",
            "--max-no-call-fraction",
            "0.2",
            "--min-reads",
            "1",
            "--compression-level",
            "1",
        ])
        .output()
        .expect("spawn fgumi filter");
    assert!(
        output.status.success(),
        "fgumi filter failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

// ---------------------------------------------------------------------------
// runall + comparison wrappers
// ---------------------------------------------------------------------------

/// Run `fgumi runall` with the given argument list. The caller is responsible
/// for spelling out every flag (start/stop, inputs, output, stage flags, etc.)
/// because each pipechain test has its own argument shape.
pub(crate) fn run_runall(args: &[&str]) {
    let output =
        Command::new(fgumi_bin()).arg("runall").args(args).output().expect("spawn fgumi runall");
    assert!(
        output.status.success(),
        "fgumi runall {} failed (exit {}):\nstdout:\n{}\nstderr:\n{}",
        args.join(" "),
        output.status,
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );
}

/// Compare two BAM files using `fgumi compare bams --command <preset>`.
///
/// Returns `true` if the comparator declared the files equivalent; on failure
/// dumps the comparator's stdout/stderr so the test failure message is useful.
pub(crate) fn compare_bams(a: &Path, b: &Path, command: &str) -> bool {
    let output = Command::new(fgumi_bin())
        .args(["compare", "bams", a.to_str().unwrap(), b.to_str().unwrap(), "--command", command])
        .output()
        .expect("spawn fgumi compare bams");
    if !output.status.success() {
        eprintln!(
            "fgumi compare bams --command {command} stdout:\n{}",
            String::from_utf8_lossy(&output.stdout)
        );
        eprintln!(
            "fgumi compare bams --command {command} stderr:\n{}",
            String::from_utf8_lossy(&output.stderr)
        );
    }
    output.status.success()
}

/// Compare two FASTQ files (typically the output of `--stop-after fastq`) by
/// raw byte equality.  The helper panics with size info on mismatch so the
/// failure message is useful.
pub(crate) fn assert_fastq_bytes_equal(a: &Path, b: &Path) {
    let av = std::fs::read(a).expect("read v1 fastq");
    let bv = std::fs::read(b).expect("read v2 fastq");
    assert_eq!(av.len(), bv.len(), "FASTQ lengths differ: v1={} v2={}", av.len(), bv.len());
    assert_eq!(av, bv, "FASTQ bytes differ");
}

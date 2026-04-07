//! End-to-end regression tests using simulate and compare commands.
//!
//! These tests validate that full pipelines produce deterministic, consistent output
//! by generating synthetic data with `simulate`, running pipeline commands, and
//! verifying outputs with `compare`.  No golden files are checked in — expected
//! output is generated fresh each run.

use std::ffi::OsString;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use tempfile::TempDir;

// ---------------------------------------------------------------------------
// Helpers: command invocation
// ---------------------------------------------------------------------------

/// Run a fgumi subcommand and return the full output.
fn fgumi(args: &[OsString]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_fgumi")).args(args).output().expect("failed to execute fgumi")
}

/// Run a fgumi subcommand, assert it succeeded, and return stdout.
fn fgumi_ok(args: &[OsString]) -> String {
    let output = fgumi(args);
    assert!(
        output.status.success(),
        "fgumi {args:?} failed:\nstdout: {}\nstderr: {}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );
    String::from_utf8_lossy(&output.stdout).to_string()
}

/// Build a command argument list from mixed string and path arguments.
macro_rules! args {
    ($($arg:expr),+ $(,)?) => {
        &[$( OsString::from($arg) ),+]
    };
}

// ---------------------------------------------------------------------------
// Helpers: simulate
// ---------------------------------------------------------------------------

/// Generate grouped-reads BAM using simulate with deterministic seed.
fn simulate_grouped_reads(output: &Path, truth: &Path, seed: u32, num_molecules: u32) {
    fgumi_ok(args![
        "simulate",
        "grouped-reads",
        "-o",
        output,
        "--truth",
        truth,
        "--num-molecules",
        &num_molecules.to_string(),
        "--seed",
        &seed.to_string(),
        "--read-length",
        "100",
        "--umi-length",
        "6",
        "--min-family-size",
        "2",
    ]);
}

/// Generate FASTQ reads using simulate with deterministic seed.
fn simulate_fastq_reads(r1: &Path, r2: &Path, truth: &Path, seed: u32, num_molecules: u32) {
    fgumi_ok(args![
        "simulate",
        "fastq-reads",
        "-1",
        r1,
        "-2",
        r2,
        "--truth",
        truth,
        "--num-molecules",
        &num_molecules.to_string(),
        "--seed",
        &seed.to_string(),
        "--read-length",
        "100",
        "--umi-length",
        "6",
        "--read-structure-r1",
        "6M94T",
        "--read-structure-r2",
        "100T",
        "--min-family-size",
        "2",
    ]);
}

// ---------------------------------------------------------------------------
// Helpers: pipeline steps
// ---------------------------------------------------------------------------

/// Run simplex consensus calling with single-threaded deterministic execution.
fn run_simplex(input: &Path, output: &Path, min_reads: u32) {
    fgumi_ok(args![
        "simplex",
        "-i",
        input,
        "-o",
        output,
        "--threads",
        "1",
        "--min-reads",
        &min_reads.to_string(),
    ]);
}

/// Run filter on a consensus BAM.
fn run_filter(input: &Path, output: &Path, min_reads: u32, min_base_quality: u32) {
    fgumi_ok(args![
        "filter",
        "-i",
        input,
        "-o",
        output,
        "--min-reads",
        &min_reads.to_string(),
        "--min-base-quality",
        &min_base_quality.to_string(),
    ]);
}

/// Run dedup on a grouped BAM.
fn run_dedup(input: &Path, output: &Path) {
    fgumi_ok(args!["dedup", "--input", input, "--output", output]);
}

// ---------------------------------------------------------------------------
// Helpers: compare
// ---------------------------------------------------------------------------

/// Compare two BAM files using the given mode, returning the full output.
fn compare_bams(bam1: &Path, bam2: &Path, mode: &str) -> Output {
    fgumi(args!["compare", "bams", bam1, bam2, "--mode", mode])
}

/// Assert that two BAM files are identical according to the given compare mode.
fn assert_bams_identical(bam1: &Path, bam2: &Path, mode: &str, context: &str) {
    let output = compare_bams(bam1, bam2, mode);
    assert!(
        output.status.success(),
        "{context}:\nstdout: {}\nstderr: {}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );
}

// ---------------------------------------------------------------------------
// Helpers: test setup
// ---------------------------------------------------------------------------

/// Create a `TempDir` and simulate grouped reads, returning (tmpdir, grouped bam path).
fn setup_grouped_reads(seed: u32, num_molecules: u32) -> (TempDir, PathBuf) {
    let tmp = TempDir::new().expect("failed to create temp dir");
    let grouped = tmp.path().join("grouped.bam");
    let truth = tmp.path().join("truth.tsv");
    simulate_grouped_reads(&grouped, &truth, seed, num_molecules);
    (tmp, grouped)
}

// ---------------------------------------------------------------------------
// Determinism: simulate produces identical output with same seed
// ---------------------------------------------------------------------------

#[test]
fn test_simulate_grouped_reads_deterministic() {
    let tmp = TempDir::new().expect("failed to create temp dir");
    let bam1 = tmp.path().join("grouped1.bam");
    let bam2 = tmp.path().join("grouped2.bam");
    let truth1 = tmp.path().join("truth1.tsv");
    let truth2 = tmp.path().join("truth2.tsv");

    simulate_grouped_reads(&bam1, &truth1, 42, 100);
    simulate_grouped_reads(&bam2, &truth2, 42, 100);

    assert_bams_identical(
        &bam1,
        &bam2,
        "full",
        "Two runs with same seed should produce identical BAMs",
    );
}

// ---------------------------------------------------------------------------
// Simplex pipeline: grouped-reads -> simplex -> deterministic output
// ---------------------------------------------------------------------------

#[test]
fn test_simplex_pipeline_deterministic() {
    let (tmp, grouped) = setup_grouped_reads(42, 200);

    let simplex1 = tmp.path().join("simplex1.bam");
    let simplex2 = tmp.path().join("simplex2.bam");
    run_simplex(&grouped, &simplex1, 1);
    run_simplex(&grouped, &simplex2, 1);

    assert_bams_identical(
        &simplex1,
        &simplex2,
        "content",
        "Two simplex runs should produce identical BAMs",
    );
}

// ---------------------------------------------------------------------------
// Filter pipeline: grouped-reads -> simplex -> filter -> deterministic output
// ---------------------------------------------------------------------------

#[test]
fn test_simplex_filter_pipeline_deterministic() {
    let (tmp, grouped) = setup_grouped_reads(99, 200);

    let simplex = tmp.path().join("simplex.bam");
    run_simplex(&grouped, &simplex, 1);

    let filtered1 = tmp.path().join("filtered1.bam");
    let filtered2 = tmp.path().join("filtered2.bam");
    run_filter(&simplex, &filtered1, 2, 10);
    run_filter(&simplex, &filtered2, 2, 10);

    assert_bams_identical(
        &filtered1,
        &filtered2,
        "content",
        "Two filter runs should produce identical BAMs",
    );
}

// ---------------------------------------------------------------------------
// Full pipeline: fastq -> extract -> group -> simplex -> filter
// ---------------------------------------------------------------------------

#[test]
fn test_full_pipeline_extract_to_filter() {
    let tmp = TempDir::new().expect("failed to create temp dir");

    // Generate synthetic FASTQ
    let r1 = tmp.path().join("r1.fq.gz");
    let r2 = tmp.path().join("r2.fq.gz");
    let truth = tmp.path().join("truth.tsv");
    simulate_fastq_reads(&r1, &r2, &truth, 42, 200);

    // Run the full pipeline twice to verify determinism
    for suffix in ["a", "b"] {
        let extracted = tmp.path().join(format!("extracted_{suffix}.bam"));
        fgumi_ok(args![
            "extract",
            "--inputs",
            &r1,
            &r2,
            "--output",
            &extracted,
            "--read-structures",
            "6M94T",
            "100T",
            "--sample",
            "test_sample",
            "--library",
            "test_lib",
        ]);

        let grouped = tmp.path().join(format!("grouped_{suffix}.bam"));
        fgumi_ok(args![
            "group",
            "--input",
            &extracted,
            "--output",
            &grouped,
            "--strategy",
            "identity",
            "--edits",
            "0",
        ]);

        let simplex = tmp.path().join(format!("simplex_{suffix}.bam"));
        run_simplex(&grouped, &simplex, 1);

        let filtered = tmp.path().join(format!("filtered_{suffix}.bam"));
        run_filter(&simplex, &filtered, 2, 10);
    }

    // Compare the two independent pipeline runs
    assert_bams_identical(
        &tmp.path().join("filtered_a.bam"),
        &tmp.path().join("filtered_b.bam"),
        "content",
        "Two full pipeline runs should produce identical BAMs",
    );
}

// ---------------------------------------------------------------------------
// Dedup pipeline: grouped-reads -> dedup -> deterministic
// ---------------------------------------------------------------------------

#[test]
fn test_dedup_pipeline_deterministic() {
    let (tmp, grouped) = setup_grouped_reads(77, 200);

    let dedup1 = tmp.path().join("dedup1.bam");
    let dedup2 = tmp.path().join("dedup2.bam");
    run_dedup(&grouped, &dedup1);
    run_dedup(&grouped, &dedup2);

    assert_bams_identical(
        &dedup1,
        &dedup2,
        "content",
        "Two dedup runs should produce identical BAMs",
    );
}

// ---------------------------------------------------------------------------
// Different seeds produce different output (sanity check)
// ---------------------------------------------------------------------------

#[test]
fn test_different_seeds_produce_different_output() {
    let tmp = TempDir::new().expect("failed to create temp dir");
    let bam1 = tmp.path().join("seed1.bam");
    let bam2 = tmp.path().join("seed2.bam");
    let truth1 = tmp.path().join("truth1.tsv");
    let truth2 = tmp.path().join("truth2.tsv");

    simulate_grouped_reads(&bam1, &truth1, 42, 100);
    simulate_grouped_reads(&bam2, &truth2, 99, 100);

    let output = compare_bams(&bam1, &bam2, "content");
    assert_eq!(
        output.status.code(),
        Some(1),
        "Expected compare to report content differences (exit 1), got {:?}\nstdout: {}\nstderr: {}",
        output.status.code(),
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );
}

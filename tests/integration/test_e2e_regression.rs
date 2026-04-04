//! End-to-end regression tests using simulate and compare commands.
//!
//! These tests validate that full pipelines produce deterministic, consistent output
//! by generating synthetic data with `simulate`, running pipeline commands, and
//! verifying outputs with `compare`.  No golden files are checked in — expected
//! output is generated fresh each run.

use std::path::Path;
use std::process::{Command, Output};
use tempfile::TempDir;

/// Run a fgumi subcommand and return the full output.
fn fgumi(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_fgumi")).args(args).output().expect("failed to execute fgumi")
}

/// Run a fgumi subcommand, assert it succeeded, and return stdout.
fn fgumi_ok(args: &[&str]) -> String {
    let output = fgumi(args);
    assert!(
        output.status.success(),
        "fgumi {} failed:\nstdout: {}\nstderr: {}",
        args.join(" "),
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );
    String::from_utf8_lossy(&output.stdout).to_string()
}

/// Helper to convert a Path to &str, panicking with a message if not valid UTF-8.
fn path_str(p: &Path) -> &str {
    p.to_str().expect("path is valid UTF-8")
}

/// Generate grouped-reads BAM using simulate with deterministic seed.
fn simulate_grouped_reads(output: &Path, truth: &Path, seed: u32, num_molecules: u32) {
    fgumi_ok(&[
        "simulate",
        "grouped-reads",
        "-o",
        path_str(output),
        "--truth",
        path_str(truth),
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
    fgumi_ok(&[
        "simulate",
        "fastq-reads",
        "-1",
        path_str(r1),
        "-2",
        path_str(r2),
        "--truth",
        path_str(truth),
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

/// Compare two BAM files using the given mode, returning (success, stdout).
fn compare_bams(bam1: &Path, bam2: &Path, mode: &str) -> (bool, String) {
    let output = fgumi(&["compare", "bams", path_str(bam1), path_str(bam2), "--mode", mode]);
    let stdout = String::from_utf8_lossy(&output.stdout).to_string();
    (output.status.success(), stdout)
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

    let (success, stdout) = compare_bams(&bam1, &bam2, "full");
    assert!(success, "Two runs with same seed should produce identical BAMs:\n{stdout}");
}

// ---------------------------------------------------------------------------
// Simplex pipeline: grouped-reads → simplex → deterministic output
// ---------------------------------------------------------------------------

#[test]
fn test_simplex_pipeline_deterministic() {
    let tmp = TempDir::new().expect("failed to create temp dir");

    // Generate grouped input
    let grouped = tmp.path().join("grouped.bam");
    let truth = tmp.path().join("truth.tsv");
    simulate_grouped_reads(&grouped, &truth, 42, 200);

    // Run simplex twice with single thread (deterministic)
    let simplex1 = tmp.path().join("simplex1.bam");
    let simplex2 = tmp.path().join("simplex2.bam");

    fgumi_ok(&[
        "simplex",
        "-i",
        path_str(&grouped),
        "-o",
        path_str(&simplex1),
        "--threads",
        "1",
        "--min-reads",
        "1",
    ]);

    fgumi_ok(&[
        "simplex",
        "-i",
        path_str(&grouped),
        "-o",
        path_str(&simplex2),
        "--threads",
        "1",
        "--min-reads",
        "1",
    ]);

    let (success, stdout) = compare_bams(&simplex1, &simplex2, "content");
    assert!(success, "Two simplex runs should produce identical BAMs:\n{stdout}");
}

// ---------------------------------------------------------------------------
// Filter pipeline: grouped-reads → simplex → filter → deterministic output
// ---------------------------------------------------------------------------

#[test]
fn test_simplex_filter_pipeline_deterministic() {
    let tmp = TempDir::new().expect("failed to create temp dir");

    // Generate grouped input
    let grouped = tmp.path().join("grouped.bam");
    let truth = tmp.path().join("truth.tsv");
    simulate_grouped_reads(&grouped, &truth, 99, 200);

    // Simplex
    let simplex = tmp.path().join("simplex.bam");
    fgumi_ok(&[
        "simplex",
        "-i",
        path_str(&grouped),
        "-o",
        path_str(&simplex),
        "--threads",
        "1",
        "--min-reads",
        "1",
    ]);

    // Filter twice
    let filtered1 = tmp.path().join("filtered1.bam");
    let filtered2 = tmp.path().join("filtered2.bam");

    for filtered in [&filtered1, &filtered2] {
        fgumi_ok(&[
            "filter",
            "-i",
            path_str(&simplex),
            "-o",
            path_str(filtered),
            "--min-reads",
            "2",
            "--min-base-quality",
            "10",
        ]);
    }

    let (success, stdout) = compare_bams(&filtered1, &filtered2, "content");
    assert!(success, "Two filter runs should produce identical BAMs:\n{stdout}");
}

// ---------------------------------------------------------------------------
// Full pipeline: fastq → extract → group → simplex → filter
// ---------------------------------------------------------------------------

#[test]
fn test_full_pipeline_extract_to_filter() {
    let tmp = TempDir::new().expect("failed to create temp dir");

    // Step 1: Generate synthetic FASTQ
    let r1 = tmp.path().join("r1.fq.gz");
    let r2 = tmp.path().join("r2.fq.gz");
    let truth = tmp.path().join("truth.tsv");
    simulate_fastq_reads(&r1, &r2, &truth, 42, 200);
    assert!(r1.exists(), "R1 FASTQ should exist");
    assert!(r2.exists(), "R2 FASTQ should exist");

    // Step 2: Extract UMIs
    let extracted = tmp.path().join("extracted.bam");
    fgumi_ok(&[
        "extract",
        "--inputs",
        path_str(&r1),
        path_str(&r2),
        "--output",
        path_str(&extracted),
        "--read-structures",
        "6M94T",
        "100T",
        "--sample",
        "test_sample",
        "--library",
        "test_lib",
    ]);
    assert!(extracted.exists(), "Extracted BAM should exist");

    // Step 3: Group by UMI
    let grouped = tmp.path().join("grouped.bam");
    fgumi_ok(&[
        "group",
        "--input",
        path_str(&extracted),
        "--output",
        path_str(&grouped),
        "--raw-tag",
        "RX",
        "--assign-tag",
        "MI",
        "--strategy",
        "identity",
        "--edits",
        "0",
    ]);
    assert!(grouped.exists(), "Grouped BAM should exist");

    // Step 4: Simplex consensus calling
    let simplex = tmp.path().join("simplex.bam");
    fgumi_ok(&[
        "simplex",
        "-i",
        path_str(&grouped),
        "-o",
        path_str(&simplex),
        "--threads",
        "1",
        "--min-reads",
        "1",
    ]);
    assert!(simplex.exists(), "Simplex BAM should exist");

    // Step 5: Filter
    let filtered = tmp.path().join("filtered.bam");
    fgumi_ok(&[
        "filter",
        "-i",
        path_str(&simplex),
        "-o",
        path_str(&filtered),
        "--min-reads",
        "2",
        "--min-base-quality",
        "10",
    ]);
    assert!(filtered.exists(), "Filtered BAM should exist");

    // Verify the filtered BAM has records (pipeline didn't silently drop everything)
    let output = fgumi_ok(&[
        "compare",
        "bams",
        path_str(&filtered),
        path_str(&filtered),
        "--mode",
        "content",
    ]);
    assert!(output.contains("IDENTICAL"), "Self-comparison should report IDENTICAL:\n{output}");
}

// ---------------------------------------------------------------------------
// Dedup pipeline: grouped-reads → dedup → deterministic
// ---------------------------------------------------------------------------

#[test]
fn test_dedup_pipeline_deterministic() {
    let tmp = TempDir::new().expect("failed to create temp dir");

    let grouped = tmp.path().join("grouped.bam");
    let truth = tmp.path().join("truth.tsv");
    simulate_grouped_reads(&grouped, &truth, 77, 200);

    let dedup1 = tmp.path().join("dedup1.bam");
    let dedup2 = tmp.path().join("dedup2.bam");

    for dedup in [&dedup1, &dedup2] {
        fgumi_ok(&["dedup", "--input", path_str(&grouped), "--output", path_str(dedup)]);
    }

    let (success, stdout) = compare_bams(&dedup1, &dedup2, "content");
    assert!(success, "Two dedup runs should produce identical BAMs:\n{stdout}");
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

    let (success, _stdout) = compare_bams(&bam1, &bam2, "content");
    assert!(!success, "Different seeds should produce different BAMs");
}

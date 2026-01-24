//! Integration tests for simulate command template-coordinate sorting.
//!
//! These tests verify that the simulated BAM files are properly sorted
//! in template-coordinate order, matching `samtools sort --template-coordinate`.

use std::process::Command;
use tempfile::TempDir;

/// Check if samtools is available in PATH.
fn samtools_available() -> bool {
    Command::new("samtools").arg("--version").output().map(|o| o.status.success()).unwrap_or(false)
}

/// Check if fgumi binary is built with simulate feature.
fn fgumi_simulate_available() -> bool {
    let output = Command::new("cargo")
        .args(["build", "--features", "simulate", "--message-format=short"])
        .output();
    output.map(|o| o.status.success()).unwrap_or(false)
}

/// Run samtools sort --template-coordinate and compare with original.
/// Returns true if the files are identical (original is already sorted).
fn verify_template_coordinate_sorted(bam_path: &std::path::Path) -> bool {
    let temp_dir = TempDir::new().unwrap();
    let sorted_path = temp_dir.path().join("sorted.bam");

    // Sort with samtools
    let status = Command::new("samtools")
        .args([
            "sort",
            "--template-coordinate",
            "-o",
            sorted_path.to_str().unwrap(),
            bam_path.to_str().unwrap(),
        ])
        .status()
        .expect("Failed to run samtools sort");

    if !status.success() {
        return false;
    }

    // Extract read names and positions from both files
    let original = Command::new("samtools")
        .args(["view", bam_path.to_str().unwrap()])
        .output()
        .expect("Failed to read original BAM");

    let sorted = Command::new("samtools")
        .args(["view", sorted_path.to_str().unwrap()])
        .output()
        .expect("Failed to read sorted BAM");

    // Compare outputs
    original.stdout == sorted.stdout
}

/// Test that mapped-reads produces template-coordinate sorted output.
#[test]
#[ignore = "requires samtools and fgumi with simulate feature"]
fn test_mapped_reads_matches_samtools() {
    if !samtools_available() {
        eprintln!("Skipping: samtools not available");
        return;
    }
    if !fgumi_simulate_available() {
        eprintln!("Skipping: fgumi simulate feature not available");
        return;
    }

    let temp_dir = TempDir::new().unwrap();
    let bam_path = temp_dir.path().join("mapped.bam");
    let truth_path = temp_dir.path().join("truth.tsv");

    // Generate mapped reads
    let status = Command::new("cargo")
        .args([
            "run",
            "--features",
            "simulate",
            "--release",
            "--",
            "simulate",
            "mapped-reads",
            "-o",
            bam_path.to_str().unwrap(),
            "--truth",
            truth_path.to_str().unwrap(),
            "--num-molecules",
            "1000",
            "--seed",
            "42",
        ])
        .status()
        .expect("Failed to run fgumi simulate mapped-reads");

    assert!(status.success(), "fgumi simulate mapped-reads failed");
    assert!(bam_path.exists(), "BAM file not created");

    // Verify sort order matches samtools
    assert!(
        verify_template_coordinate_sorted(&bam_path),
        "mapped-reads output does not match samtools sort --template-coordinate"
    );
}

/// Test that grouped-reads (simplex) produces template-coordinate sorted output.
#[test]
#[ignore = "requires samtools and fgumi with simulate feature"]
fn test_grouped_reads_simplex_matches_samtools() {
    if !samtools_available() {
        eprintln!("Skipping: samtools not available");
        return;
    }
    if !fgumi_simulate_available() {
        eprintln!("Skipping: fgumi simulate feature not available");
        return;
    }

    let temp_dir = TempDir::new().unwrap();
    let bam_path = temp_dir.path().join("grouped.bam");
    let truth_path = temp_dir.path().join("truth.tsv");

    // Generate grouped reads (simplex)
    let status = Command::new("cargo")
        .args([
            "run",
            "--features",
            "simulate",
            "--release",
            "--",
            "simulate",
            "grouped-reads",
            "-o",
            bam_path.to_str().unwrap(),
            "--truth",
            truth_path.to_str().unwrap(),
            "--num-molecules",
            "1000",
            "--seed",
            "42",
        ])
        .status()
        .expect("Failed to run fgumi simulate grouped-reads");

    assert!(status.success(), "fgumi simulate grouped-reads failed");
    assert!(bam_path.exists(), "BAM file not created");

    // Verify sort order matches samtools
    assert!(
        verify_template_coordinate_sorted(&bam_path),
        "grouped-reads simplex output does not match samtools sort --template-coordinate"
    );
}

/// Test that grouped-reads (duplex) produces template-coordinate sorted output.
#[test]
#[ignore = "requires samtools and fgumi with simulate feature"]
fn test_grouped_reads_duplex_matches_samtools() {
    if !samtools_available() {
        eprintln!("Skipping: samtools not available");
        return;
    }
    if !fgumi_simulate_available() {
        eprintln!("Skipping: fgumi simulate feature not available");
        return;
    }

    let temp_dir = TempDir::new().unwrap();
    let bam_path = temp_dir.path().join("grouped_duplex.bam");
    let truth_path = temp_dir.path().join("truth.tsv");

    // Generate grouped reads (duplex)
    let status = Command::new("cargo")
        .args([
            "run",
            "--features",
            "simulate",
            "--release",
            "--",
            "simulate",
            "grouped-reads",
            "-o",
            bam_path.to_str().unwrap(),
            "--truth",
            truth_path.to_str().unwrap(),
            "--num-molecules",
            "1000",
            "--duplex",
            "--seed",
            "42",
        ])
        .status()
        .expect("Failed to run fgumi simulate grouped-reads --duplex");

    assert!(status.success(), "fgumi simulate grouped-reads --duplex failed");
    assert!(bam_path.exists(), "BAM file not created");

    // Verify sort order matches samtools
    assert!(
        verify_template_coordinate_sorted(&bam_path),
        "grouped-reads duplex output does not match samtools sort --template-coordinate"
    );
}

/// Test with larger dataset to catch edge cases.
#[test]
#[ignore = "requires samtools and fgumi with simulate feature; slow"]
fn test_large_dataset_matches_samtools() {
    if !samtools_available() {
        eprintln!("Skipping: samtools not available");
        return;
    }
    if !fgumi_simulate_available() {
        eprintln!("Skipping: fgumi simulate feature not available");
        return;
    }

    let temp_dir = TempDir::new().unwrap();
    let bam_path = temp_dir.path().join("large.bam");
    let truth_path = temp_dir.path().join("truth.tsv");

    // Generate larger dataset to catch edge cases
    let status = Command::new("cargo")
        .args([
            "run",
            "--features",
            "simulate",
            "--release",
            "--",
            "simulate",
            "grouped-reads",
            "-o",
            bam_path.to_str().unwrap(),
            "--truth",
            truth_path.to_str().unwrap(),
            "--num-molecules",
            "10000",
            "--duplex",
            "--seed",
            "12345",
        ])
        .status()
        .expect("Failed to run fgumi simulate grouped-reads");

    assert!(status.success());
    assert!(verify_template_coordinate_sorted(&bam_path));
}

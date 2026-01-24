//! Integration tests for sort --write-index functionality.
//!
//! These tests verify that the BAM index generation works correctly
//! for both fast (raw-bytes) and standard (record buffer) sorting modes.

use std::path::Path;
use std::process::Command;
use tempfile::TempDir;

/// Check if samtools is available in PATH.
fn samtools_available() -> bool {
    Command::new("samtools").arg("--version").output().map(|o| o.status.success()).unwrap_or(false)
}

/// Check if the fgumi binary exists.
fn fgumi_binary_path() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_BIN_EXE_fgumi"))
}

/// Create a small test BAM file using samtools.
fn create_test_bam(dir: &Path) -> std::path::PathBuf {
    let bam_path = dir.join("test_input.bam");

    // Create a simple SAM file and convert to BAM
    let sam_content = r"@HD	VN:1.6	SO:unsorted
@SQ	SN:chr1	LN:10000
@SQ	SN:chr2	LN:10000
read1	0	chr1	100	60	10M	*	0	0	ACGTACGTAC	IIIIIIIIII
read2	0	chr1	200	60	10M	*	0	0	ACGTACGTAC	IIIIIIIIII
read3	0	chr1	300	60	10M	*	0	0	ACGTACGTAC	IIIIIIIIII
read4	0	chr2	100	60	10M	*	0	0	ACGTACGTAC	IIIIIIIIII
read5	0	chr2	200	60	10M	*	0	0	ACGTACGTAC	IIIIIIIIII
read6	4	*	0	0	*	*	0	0	ACGTACGTAC	IIIIIIIIII
";

    let sam_path = dir.join("test_input.sam");
    std::fs::write(&sam_path, sam_content).expect("Failed to write SAM file");

    // Convert SAM to BAM using samtools
    let status = Command::new("samtools")
        .args(["view", "-b", "-o", bam_path.to_str().unwrap(), sam_path.to_str().unwrap()])
        .status()
        .expect("Failed to run samtools view");

    assert!(status.success(), "samtools view failed");
    bam_path
}

/// Verify that a BAM index allows region queries.
fn verify_index_works(bam_path: &Path) -> bool {
    // Try to query a region - this will fail if index is invalid
    let output = Command::new("samtools")
        .args(["view", "-c", bam_path.to_str().unwrap(), "chr1:100-300"])
        .output()
        .expect("Failed to run samtools view");

    if !output.status.success() {
        eprintln!("samtools view failed: {}", String::from_utf8_lossy(&output.stderr));
        return false;
    }

    // Should find at least one read
    let count: i32 = String::from_utf8_lossy(&output.stdout).trim().parse().unwrap_or(0);
    count > 0
}

/// Test sort --fast --write-index creates a valid index.
#[test]
#[ignore = "requires samtools"]
fn test_sort_fast_write_index() {
    if !samtools_available() {
        eprintln!("Skipping: samtools not available");
        return;
    }

    let temp_dir = TempDir::new().unwrap();
    let input_bam = create_test_bam(temp_dir.path());
    let sorted_bam_path = temp_dir.path().join("sorted.bam");
    let index_path = temp_dir.path().join("sorted.bam.bai");

    // Run fgumi sort --fast --write-index
    let status = Command::new(fgumi_binary_path())
        .args([
            "sort",
            "--fast",
            "-i",
            input_bam.to_str().unwrap(),
            "-o",
            sorted_bam_path.to_str().unwrap(),
            "--order",
            "coordinate",
            "--write-index",
            "-m",
            "10M",
        ])
        .status()
        .expect("Failed to run fgumi sort");

    assert!(status.success(), "fgumi sort --fast --write-index failed");
    assert!(sorted_bam_path.exists(), "Output BAM not created");
    assert!(index_path.exists(), "Output BAI not created");

    // Verify index works
    assert!(verify_index_works(&sorted_bam_path), "Index does not work for region queries");
}

/// Test sort (non-fast) --write-index creates a valid index.
#[test]
#[ignore = "requires samtools"]
fn test_sort_standard_write_index() {
    if !samtools_available() {
        eprintln!("Skipping: samtools not available");
        return;
    }

    let temp_dir = TempDir::new().unwrap();
    let input_bam = create_test_bam(temp_dir.path());
    let sorted_bam_path = temp_dir.path().join("sorted.bam");
    let index_path = temp_dir.path().join("sorted.bam.bai");

    // Run fgumi sort (non-fast) --write-index
    let status = Command::new(fgumi_binary_path())
        .args([
            "sort",
            "-i",
            input_bam.to_str().unwrap(),
            "-o",
            sorted_bam_path.to_str().unwrap(),
            "--order",
            "coordinate",
            "--write-index",
            "-m",
            "10M",
        ])
        .status()
        .expect("Failed to run fgumi sort");

    assert!(status.success(), "fgumi sort --write-index failed");
    assert!(sorted_bam_path.exists(), "Output BAM not created");
    assert!(index_path.exists(), "Output BAI not created");

    // Verify index works
    assert!(verify_index_works(&sorted_bam_path), "Index does not work for region queries");
}

/// Test that --write-index with multi-threading still produces valid index.
#[test]
#[ignore = "requires samtools"]
fn test_sort_write_index_multithreaded() {
    if !samtools_available() {
        eprintln!("Skipping: samtools not available");
        return;
    }

    let temp_dir = TempDir::new().unwrap();
    let input_bam = create_test_bam(temp_dir.path());
    let sorted_bam_path = temp_dir.path().join("sorted.bam");
    let index_path = temp_dir.path().join("sorted.bam.bai");

    // Run fgumi sort --fast --write-index with multiple threads
    let status = Command::new(fgumi_binary_path())
        .args([
            "sort",
            "--fast",
            "-i",
            input_bam.to_str().unwrap(),
            "-o",
            sorted_bam_path.to_str().unwrap(),
            "--order",
            "coordinate",
            "--write-index",
            "-m",
            "10M",
            "-@",
            "4",
        ])
        .status()
        .expect("Failed to run fgumi sort");

    assert!(status.success(), "fgumi sort --fast --write-index -@ 4 failed");
    assert!(sorted_bam_path.exists(), "Output BAM not created");
    assert!(index_path.exists(), "Output BAI not created");

    // Verify index works
    assert!(verify_index_works(&sorted_bam_path), "Index does not work for region queries");
}

/// Test that --write-index fails for non-coordinate sort orders.
#[test]
fn test_write_index_requires_coordinate_sort() {
    let temp_dir = TempDir::new().unwrap();
    let sorted_bam_path = temp_dir.path().join("sorted.bam");

    // Run fgumi sort with queryname order and --write-index (should fail)
    let output = Command::new(fgumi_binary_path())
        .args([
            "sort",
            "-i",
            "/dev/null", // Won't get to reading input
            "-o",
            sorted_bam_path.to_str().unwrap(),
            "--order",
            "queryname",
            "--write-index",
        ])
        .output()
        .expect("Failed to run fgumi sort");

    assert!(!output.status.success(), "fgumi sort --order queryname --write-index should fail");

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("--write-index is only valid for coordinate sort"),
        "Expected error about coordinate sort, got: {stderr}"
    );
}

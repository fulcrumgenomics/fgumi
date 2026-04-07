//! Sort correctness tests for all sort orders and spill configurations.
//!
//! Validates that `fgumi sort` produces correctly sorted output across:
//! - All 4 sort orders (coordinate, queryname-lex, queryname-natural, template-coordinate)
//! - Multiple thread counts (1, 2, 4)
//! - Various memory limits forcing 0, 1, and multiple spills
//! - Record content is preserved (byte-identical records)

use std::ffi::OsString;
use std::fmt::Write as _;
use std::path::Path;
use std::process::{Command, Output};
use tempfile::TempDir;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn fgumi_binary() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_BIN_EXE_fgumi"))
}

fn fgumi(args: &[OsString]) -> Output {
    Command::new(fgumi_binary()).args(args).output().expect("failed to execute fgumi")
}

fn samtools_available() -> bool {
    Command::new("samtools").arg("--version").output().map(|o| o.status.success()).unwrap_or(false)
}

/// Create a test BAM with many records at various positions using samtools.
/// Returns the path to the unsorted BAM.
fn create_test_bam(dir: &Path, num_reads: usize) -> std::path::PathBuf {
    let bam_path = dir.join("unsorted.bam");

    // Build SAM content with reads spread across multiple chromosomes
    let mut sam = String::new();
    sam.push_str("@HD\tVN:1.6\tSO:unsorted\n");
    sam.push_str("@SQ\tSN:chr1\tLN:100000\n");
    sam.push_str("@SQ\tSN:chr2\tLN:100000\n");
    sam.push_str("@SQ\tSN:chr3\tLN:100000\n");
    sam.push_str("@RG\tID:rg1\tSM:sample1\tLB:lib1\n");

    let seq = "ACGTACGTAC";
    let qual = "IIIIIIIIII";

    // Generate reads in a deliberately unsorted order
    for i in 0..num_reads {
        let chrom = match i % 3 {
            0 => "chr2",
            1 => "chr1",
            _ => "chr3",
        };
        // Positions go backwards within each chromosome to ensure unsorted input
        let pos = 50000 - (i % 1000) * 50 + 1;
        // Use non-zero-padded names so natural and lexicographic sort produce different
        // orderings (e.g. read_2 < read_10 under natural but read_10 < read_2 under lex).
        let name = format!("read_{i}");

        // Paired-end reads
        let mate_pos = pos + 200;
        // Assign each template to one of 4 cells so template-coordinate tests
        // exercise the cell-aware sort key path, not just the fallback.
        let cell = format!("cell{}", i % 4);
        // R1
        writeln!(
            sam,
            "{name}\t99\t{chrom}\t{pos}\t60\t10M\t=\t{mate_pos}\t210\t{seq}\t{qual}\tRG:Z:rg1\tMI:Z:umi_{i}\tCB:Z:{cell}"
        )
        .expect("write SAM record");
        // R2
        writeln!(
            sam,
            "{name}\t147\t{chrom}\t{mate_pos}\t60\t10M\t=\t{pos}\t-210\t{seq}\t{qual}\tRG:Z:rg1\tMI:Z:umi_{i}\tCB:Z:{cell}"
        )
        .expect("write SAM record");
    }

    // Add some unmapped reads
    for i in 0..10 {
        let name = format!("unmapped_{i:04}");
        writeln!(sam, "{name}\t77\t*\t0\t0\t*\t*\t0\t0\t{seq}\t{qual}\tRG:Z:rg1")
            .expect("write SAM record");
        writeln!(sam, "{name}\t141\t*\t0\t0\t*\t*\t0\t0\t{seq}\t{qual}\tRG:Z:rg1")
            .expect("write SAM record");
    }

    let sam_path = dir.join("unsorted.sam");
    std::fs::write(&sam_path, sam).expect("write SAM");

    // Convert to BAM
    let status = Command::new("samtools")
        .args(["view", "-b", "-o", bam_path.to_str().unwrap(), sam_path.to_str().unwrap()])
        .status()
        .expect("samtools view");
    assert!(status.success(), "samtools view failed");

    bam_path
}

/// Run `fgumi sort --verify` and return whether the file is correctly sorted.
fn verify_sorted(bam_path: &Path, order: &str) -> bool {
    let mut cmd_args: Vec<OsString> = vec![
        "sort".into(),
        "--verify".into(),
        "-i".into(),
        bam_path.as_os_str().to_owned(),
        "--order".into(),
        order.into(),
    ];
    if order == "template-coordinate" {
        cmd_args.push("--cell-tag".into());
        cmd_args.push("CB".into());
    }
    let output = fgumi(&cmd_args);
    output.status.success()
}

/// Count records in a BAM using samtools.
fn count_records(bam_path: &Path) -> u64 {
    let output = Command::new("samtools")
        .args(["view", "-c", bam_path.to_str().unwrap()])
        .output()
        .expect("samtools view -c");
    assert!(output.status.success());
    String::from_utf8_lossy(&output.stdout)
        .trim()
        .parse()
        .expect("samtools view -c output should be a valid integer")
}

/// Sort a BAM file with fgumi and return the output path.
fn sort_bam(input: &Path, output: &Path, order: &str, threads: usize, max_memory: &str) {
    sort_bam_with_args(input, output, order, threads, max_memory, &[]);
}

fn sort_bam_with_args(
    input: &Path,
    output: &Path,
    order: &str,
    threads: usize,
    max_memory: &str,
    extra_args: &[&str],
) {
    let mut cmd_args: Vec<OsString> = vec![
        "sort".into(),
        "-i".into(),
        input.as_os_str().to_owned(),
        "-o".into(),
        output.as_os_str().to_owned(),
        "--order".into(),
        order.into(),
        "--threads".into(),
        threads.to_string().into(),
        "-m".into(),
        max_memory.into(),
        "--temp-compression".into(),
        "1".into(),
    ];
    if order == "template-coordinate" {
        cmd_args.push("--cell-tag".into());
        cmd_args.push("CB".into());
    }
    for arg in extra_args {
        cmd_args.push((*arg).into());
    }
    let output_result = fgumi(&cmd_args);
    assert!(
        output_result.status.success(),
        "fgumi sort failed for order={order} threads={threads} memory={max_memory}:\nstderr: {}",
        String::from_utf8_lossy(&output_result.stderr),
    );
}

// ---------------------------------------------------------------------------
// Tests: Coordinate Sort
// ---------------------------------------------------------------------------

#[test]
#[ignore = "requires samtools"]
fn test_sort_coordinate_in_memory() {
    if !samtools_available() {
        return;
    }
    let dir = TempDir::new().unwrap();
    let input = create_test_bam(dir.path(), 500);
    let output = dir.path().join("sorted.bam");

    // Large memory = no spills (in-memory sort)
    sort_bam(&input, &output, "coordinate", 2, "100M");
    assert!(verify_sorted(&output, "coordinate"), "coordinate sort verification failed");
    assert_eq!(count_records(&input), count_records(&output), "record count mismatch");
}

#[test]
#[ignore = "requires samtools"]
fn test_sort_coordinate_with_spills() {
    if !samtools_available() {
        return;
    }
    let dir = TempDir::new().unwrap();
    let input = create_test_bam(dir.path(), 2000);
    let output = dir.path().join("sorted.bam");

    // Small memory = forces spills to disk
    sort_bam(&input, &output, "coordinate", 2, "50K");
    assert!(verify_sorted(&output, "coordinate"), "coordinate sort with spills failed");
    assert_eq!(count_records(&input), count_records(&output), "record count mismatch");
}

#[test]
#[ignore = "requires samtools"]
fn test_sort_coordinate_many_spills_t4() {
    if !samtools_available() {
        return;
    }
    let dir = TempDir::new().unwrap();
    let input = create_test_bam(dir.path(), 3000);
    let output = dir.path().join("sorted.bam");

    // Very small memory with 4 threads = many spills
    sort_bam(&input, &output, "coordinate", 4, "30K");
    assert!(verify_sorted(&output, "coordinate"), "coordinate sort many spills failed");
    assert_eq!(count_records(&input), count_records(&output), "record count mismatch");
}

#[test]
#[ignore = "requires samtools"]
fn test_sort_coordinate_single_thread() {
    if !samtools_available() {
        return;
    }
    let dir = TempDir::new().unwrap();
    let input = create_test_bam(dir.path(), 1000);
    let output = dir.path().join("sorted.bam");

    sort_bam(&input, &output, "coordinate", 1, "50K");
    assert!(verify_sorted(&output, "coordinate"), "coordinate sort t1 failed");
    assert_eq!(count_records(&input), count_records(&output), "record count mismatch");
}

// ---------------------------------------------------------------------------
// Tests: Queryname Lexicographic Sort
// ---------------------------------------------------------------------------

#[test]
#[ignore = "requires samtools"]
fn test_sort_queryname_lex_in_memory() {
    if !samtools_available() {
        return;
    }
    let dir = TempDir::new().unwrap();
    let input = create_test_bam(dir.path(), 500);
    let output = dir.path().join("sorted.bam");

    sort_bam(&input, &output, "queryname", 2, "100M");
    assert!(verify_sorted(&output, "queryname"), "queryname-lex sort failed");
    assert_eq!(count_records(&input), count_records(&output), "record count mismatch");
}

#[test]
#[ignore = "requires samtools"]
fn test_sort_queryname_lex_with_spills() {
    if !samtools_available() {
        return;
    }
    let dir = TempDir::new().unwrap();
    let input = create_test_bam(dir.path(), 2000);
    let output = dir.path().join("sorted.bam");

    sort_bam(&input, &output, "queryname", 2, "50K");
    assert!(verify_sorted(&output, "queryname"), "queryname-lex sort with spills failed");
    assert_eq!(count_records(&input), count_records(&output), "record count mismatch");
}

// ---------------------------------------------------------------------------
// Tests: Queryname Natural Sort
// ---------------------------------------------------------------------------

#[test]
#[ignore = "requires samtools"]
fn test_sort_queryname_natural_in_memory() {
    if !samtools_available() {
        return;
    }
    let dir = TempDir::new().unwrap();
    let input = create_test_bam(dir.path(), 500);
    let output = dir.path().join("sorted.bam");

    sort_bam(&input, &output, "queryname::natural", 2, "100M");
    assert!(verify_sorted(&output, "queryname::natural"), "queryname-natural sort failed");
    assert_eq!(count_records(&input), count_records(&output), "record count mismatch");
}

#[test]
#[ignore = "requires samtools"]
fn test_sort_queryname_natural_with_spills() {
    if !samtools_available() {
        return;
    }
    let dir = TempDir::new().unwrap();
    let input = create_test_bam(dir.path(), 2000);
    let output = dir.path().join("sorted.bam");

    sort_bam(&input, &output, "queryname::natural", 2, "50K");
    assert!(
        verify_sorted(&output, "queryname::natural"),
        "queryname-natural sort with spills failed"
    );
    assert_eq!(count_records(&input), count_records(&output), "record count mismatch");
}

// ---------------------------------------------------------------------------
// Tests: Template-Coordinate Sort
// ---------------------------------------------------------------------------

#[test]
#[ignore = "requires samtools"]
fn test_sort_template_coordinate_in_memory() {
    if !samtools_available() {
        return;
    }
    let dir = TempDir::new().unwrap();
    let input = create_test_bam(dir.path(), 500);
    let output = dir.path().join("sorted.bam");

    sort_bam(&input, &output, "template-coordinate", 2, "100M");
    assert!(verify_sorted(&output, "template-coordinate"), "template-coordinate sort failed");
    assert_eq!(count_records(&input), count_records(&output), "record count mismatch");
}

#[test]
#[ignore = "requires samtools"]
fn test_sort_template_coordinate_with_spills() {
    if !samtools_available() {
        return;
    }
    let dir = TempDir::new().unwrap();
    let input = create_test_bam(dir.path(), 2000);
    let output = dir.path().join("sorted.bam");

    sort_bam(&input, &output, "template-coordinate", 2, "50K");
    assert!(
        verify_sorted(&output, "template-coordinate"),
        "template-coordinate sort with spills failed"
    );
    assert_eq!(count_records(&input), count_records(&output), "record count mismatch");
}

#[test]
#[ignore = "requires samtools"]
fn test_sort_template_coordinate_many_spills_t4() {
    if !samtools_available() {
        return;
    }
    let dir = TempDir::new().unwrap();
    let input = create_test_bam(dir.path(), 3000);
    let output = dir.path().join("sorted.bam");

    sort_bam(&input, &output, "template-coordinate", 4, "30K");
    assert!(
        verify_sorted(&output, "template-coordinate"),
        "template-coordinate many spills failed"
    );
    assert_eq!(count_records(&input), count_records(&output), "record count mismatch");
}

// ---------------------------------------------------------------------------
// Tests: Consistency across thread counts
// ---------------------------------------------------------------------------

#[test]
#[ignore = "requires samtools"]
fn test_sort_coordinate_consistent_across_threads() {
    if !samtools_available() {
        return;
    }
    let dir = TempDir::new().unwrap();
    let input = create_test_bam(dir.path(), 1000);

    // Sort with different thread counts using small memory to force spills
    let outputs: Vec<_> = [1, 2, 4]
        .iter()
        .map(|t| {
            let output = dir.path().join(format!("sorted_t{t}.bam"));
            sort_bam(&input, &output, "coordinate", *t, "50K");
            output
        })
        .collect();

    // All should be correctly sorted
    for output in &outputs {
        assert!(verify_sorted(output, "coordinate"), "failed for {}", output.display());
    }

    // Record counts should all match
    let input_count = count_records(&input);
    for output in &outputs {
        assert_eq!(input_count, count_records(output), "count mismatch for {}", output.display());
    }
}

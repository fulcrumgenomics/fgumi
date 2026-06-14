//! Integration tests for sort --write-index functionality.
//!
//! Verifies that BAM index generation produces a valid `.bai` that
//! samtools can query, across the coordinate sort code path.
//!
//! fgumi sort is invoked in-process via `Sort::execute()`; samtools
//! is invoked as a subprocess (external tool).

use clap::Parser;
use fgumi_lib::commands::command::Command as FgumiCommand;
use fgumi_lib::commands::sort::Sort;
use rstest::rstest;
use std::ffi::OsStr;
use std::path::Path;
use std::process::Command;
use tempfile::TempDir;

/// Check if samtools is available in PATH.
fn samtools_available() -> bool {
    Command::new("samtools").arg("--version").output().map(|o| o.status.success()).unwrap_or(false)
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

/// Sort `input` into `output` (coordinate order) in-process via `Sort::execute`.
///
/// `write_index` toggles `--write-index`; `threads` sets `-@`. Paths are passed
/// as `&OsStr` so `try_parse_from` does not UTF-8 round-trip non-UTF-8 temp dirs.
fn run_coordinate_sort(input: &Path, output: &Path, write_index: bool, threads: usize) {
    let threads_str = threads.to_string();
    let mut args: Vec<&OsStr> = vec![
        OsStr::new("sort"),
        OsStr::new("-i"),
        input.as_os_str(),
        OsStr::new("-o"),
        output.as_os_str(),
        OsStr::new("--order"),
        OsStr::new("coordinate"),
        OsStr::new("-m"),
        OsStr::new("10M"),
        OsStr::new("-@"),
        OsStr::new(&threads_str),
    ];
    if write_index {
        args.push(OsStr::new("--write-index"));
    }

    let cmd = Sort::try_parse_from(args).expect("failed to parse sort args");
    cmd.execute("fgumi sort").expect("fgumi sort failed");
}

/// Return the sorted read names returned by a samtools region query against a
/// (BAI-indexed) coordinate-sorted BAM.
fn region_query_names(bam_path: &Path, region: &str) -> Vec<String> {
    let output = Command::new("samtools")
        .args(["view", bam_path.to_str().unwrap(), region])
        .output()
        .expect("Failed to run samtools view");
    assert!(
        output.status.success(),
        "samtools region query failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let mut names: Vec<String> = String::from_utf8_lossy(&output.stdout)
        .lines()
        .filter_map(|line| line.split('\t').next().map(str::to_string))
        .collect();
    names.sort();
    names
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

/// Test sort --write-index creates a valid index.
#[test]
#[ignore = "requires samtools"]
fn test_sort_write_index() {
    if !samtools_available() {
        eprintln!("Skipping: samtools not available");
        return;
    }

    let temp_dir = TempDir::new().unwrap();
    let input_bam = create_test_bam(temp_dir.path());
    let sorted_bam_path = temp_dir.path().join("sorted.bam");
    let index_path = temp_dir.path().join("sorted.bam.bai");

    // Run fgumi sort (non-fast) --write-index. Pass `&OsStr` directly so
    // `try_parse_from` doesn't UTF-8 round-trip the temp paths (would panic
    // on a non-UTF-8 temp dir).
    let cmd = Sort::try_parse_from([
        OsStr::new("sort"),
        OsStr::new("-i"),
        input_bam.as_os_str(),
        OsStr::new("-o"),
        sorted_bam_path.as_os_str(),
        OsStr::new("--order"),
        OsStr::new("coordinate"),
        OsStr::new("--write-index"),
        OsStr::new("-m"),
        OsStr::new("10M"),
    ])
    .expect("failed to parse sort args");

    cmd.execute("fgumi sort").expect("fgumi sort --write-index failed");
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

    // Run fgumi sort --write-index with multiple threads
    let cmd = Sort::try_parse_from([
        OsStr::new("sort"),
        OsStr::new("-i"),
        input_bam.as_os_str(),
        OsStr::new("-o"),
        sorted_bam_path.as_os_str(),
        OsStr::new("--order"),
        OsStr::new("coordinate"),
        OsStr::new("--write-index"),
        OsStr::new("-m"),
        OsStr::new("10M"),
        OsStr::new("-@"),
        OsStr::new("4"),
    ])
    .expect("failed to parse sort args");

    cmd.execute("fgumi sort").expect("fgumi sort --write-index -@ 4 failed");
    assert!(sorted_bam_path.exists(), "Output BAM not created");
    assert!(index_path.exists(), "Output BAI not created");

    // Verify index works
    assert!(verify_index_works(&sorted_bam_path), "Index does not work for region queries");
}

/// Regression guard for issue #330 Item 3: the `--write-index`-on and -off
/// coordinate-sort runs must produce **byte-identical** BAM output.
///
/// Before the in-sort indexer was removed, `--write-index` wrote the BAM
/// through a different writer backend (noodles' `MultithreadedWriter`) than
/// the off path (`PooledBamWriter`), so the two BGZF streams diverged
/// byte-for-byte. Now that BAI generation is a decoupled post-write pass
/// (`IndexBamFinalizeHook`), both runs go through the identical writer and the
/// BAM bytes must match exactly. This is checked across thread counts because
/// the writer backend is the only axis that ever differed.
#[rstest]
#[case::threads_1(1)]
#[case::threads_4(4)]
#[ignore = "requires samtools"]
fn test_sort_write_index_byte_identical_to_off(#[case] threads: usize) {
    if !samtools_available() {
        eprintln!("Skipping: samtools not available");
        return;
    }

    let temp_dir = TempDir::new().unwrap();
    let input_bam = create_test_bam(temp_dir.path());

    let with_index = temp_dir.path().join("with_index.bam");
    let without_index = temp_dir.path().join("without_index.bam");

    run_coordinate_sort(&input_bam, &with_index, true, threads);
    run_coordinate_sort(&input_bam, &without_index, false, threads);

    assert!(
        with_index.with_extension("bam.bai").exists(),
        "BAI not produced for --write-index run (threads={threads})"
    );

    let with_bytes = std::fs::read(&with_index).expect("read with-index BAM");
    let without_bytes = std::fs::read(&without_index).expect("read without-index BAM");
    assert_eq!(
        with_bytes, without_bytes,
        "--write-index BAM bytes diverge from --write-index-off bytes (threads={threads}); \
         the two runs must share the same writer backend"
    );
}

/// Cross-check that fgumi's BAI returns exactly the same region-query record set
/// as a samtools-generated index over the same coordinate-sorted BAM, across
/// thread counts. Guards the post-write `IndexBamFinalizeHook` against subtle
/// virtual-offset / binning errors that `verify_index_works` (count > 0) misses.
#[rstest]
#[case::t1_chr1_100_300(1, "chr1:100-300")]
#[case::t1_chr1_1_10000(1, "chr1:1-10000")]
#[case::t1_chr2(1, "chr2")]
#[case::t1_chr2_200_200(1, "chr2:200-200")]
#[case::t4_chr1_100_300(4, "chr1:100-300")]
#[case::t4_chr1_1_10000(4, "chr1:1-10000")]
#[case::t4_chr2(4, "chr2")]
#[case::t4_chr2_200_200(4, "chr2:200-200")]
#[ignore = "requires samtools"]
fn test_sort_write_index_matches_samtools_bai(#[case] threads: usize, #[case] region: &str) {
    if !samtools_available() {
        eprintln!("Skipping: samtools not available");
        return;
    }

    let temp_dir = TempDir::new().unwrap();
    let input_bam = create_test_bam(temp_dir.path());

    // fgumi: sort + index in one pass.
    let fgumi_bam = temp_dir.path().join("fgumi_sorted.bam");
    run_coordinate_sort(&input_bam, &fgumi_bam, true, threads);

    // samtools reference: sort the same input, then index it.
    let sam_bam = temp_dir.path().join("samtools_sorted.bam");
    let status = Command::new("samtools")
        .args(["sort", "-o", sam_bam.to_str().unwrap(), input_bam.to_str().unwrap()])
        .status()
        .expect("Failed to run samtools sort");
    assert!(status.success(), "samtools sort failed");
    let status = Command::new("samtools")
        .args(["index", sam_bam.to_str().unwrap()])
        .status()
        .expect("Failed to run samtools index");
    assert!(status.success(), "samtools index failed");

    let fgumi_names = region_query_names(&fgumi_bam, region);
    let samtools_names = region_query_names(&sam_bam, region);
    assert_eq!(
        fgumi_names, samtools_names,
        "fgumi BAI region query disagrees with samtools BAI for {region} (threads={threads})"
    );
}

/// Test that --write-index fails for non-coordinate sort orders.
#[test]
fn test_write_index_requires_coordinate_sort() {
    let temp_dir = TempDir::new().unwrap();
    let sorted_bam_path = temp_dir.path().join("sorted.bam");

    // Run fgumi sort with queryname order and --write-index (should fail)
    let cmd = Sort::try_parse_from([
        OsStr::new("sort"),
        OsStr::new("-i"),
        OsStr::new("/dev/null"), // Won't get to reading input
        OsStr::new("-o"),
        sorted_bam_path.as_os_str(),
        OsStr::new("--order"),
        OsStr::new("queryname"),
        OsStr::new("--write-index"),
    ])
    .expect("failed to parse sort args");

    let result = cmd.execute("fgumi sort");
    assert!(result.is_err(), "fgumi sort --order queryname --write-index should fail");

    let err_msg = format!("{:#}", result.unwrap_err());
    assert!(
        err_msg.contains("--write-index is only valid for coordinate sort"),
        "Expected error about coordinate sort, got: {err_msg}"
    );
}

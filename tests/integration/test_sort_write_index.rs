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
use std::ffi::OsStr;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};
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

/// Build a larger unsorted BAM: `mapped_per_contig` mapped reads on each of
/// chr1/chr2 at scattered positions, plus some unmapped reads. Large enough to
/// span many BGZF blocks so cross-block virtual offsets are exercised.
fn create_large_test_bam(dir: &Path, mapped_per_contig: usize) -> PathBuf {
    let mut sam = String::from(
        "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000000\n@SQ\tSN:chr2\tLN:1000000\n",
    );
    // Interleave contigs and scatter positions (unsorted input); fgumi sorts it.
    for i in 0..mapped_per_contig {
        let p1 = 1 + (i * 37) % 999_000;
        let p2 = 1 + (i * 53 + 11) % 999_000;
        let _ = writeln!(sam, "r1_{i}\t0\tchr1\t{p1}\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII");
        let _ = writeln!(sam, "r2_{i}\t0\tchr2\t{p2}\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII");
    }
    for i in 0..100 {
        let _ = writeln!(sam, "u_{i}\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII");
    }
    // Placed-but-unmapped reads: FUNMAP (0x4) set, but with a reference + position
    // (e.g. the unmapped mate of a mapped read). samtools position-bins these so
    // they appear in region queries; fgumi's indexer classifies them as unmapped
    // and does not position-bin them (pre-existing behavior, unchanged by the
    // pooled-writer work). The test therefore compares MAPPED-only retrieval.
    for i in 0..2000 {
        let p1 = 1 + (i * 101) % 999_000;
        let p2 = 1 + (i * 149 + 7) % 999_000;
        let _ = writeln!(sam, "pu1_{i}\t4\tchr1\t{p1}\t0\t*\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII");
        let _ = writeln!(sam, "pu2_{i}\t4\tchr2\t{p2}\t0\t*\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII");
    }

    let sam_path = dir.join("large_input.sam");
    std::fs::write(&sam_path, sam).expect("write SAM");
    let bam_path = dir.join("large_input.bam");
    let status = Command::new("samtools")
        .args(["view", "-b", "-o", bam_path.to_str().unwrap(), sam_path.to_str().unwrap()])
        .status()
        .expect("run samtools view");
    assert!(status.success(), "samtools view failed");
    bam_path
}

/// `samtools view -c [extra] <bam> <region>` — record count in a region via the
/// `.bai`. `extra` carries filter flags (e.g. `-F 4` to count mapped reads only).
fn region_count(bam: &Path, region: &str, extra: &[&str]) -> i64 {
    let out = Command::new("samtools")
        .arg("view")
        .arg("-c")
        .args(extra)
        .args([bam.to_str().unwrap(), region])
        .output()
        .expect("run samtools view -c");
    assert!(
        out.status.success(),
        "samtools view -c {region} failed on {}: {}",
        bam.display(),
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8_lossy(&out.stdout).trim().parse().expect("parse count")
}

/// Per-reference MAPPED counts (`ref\tlen\tmapped`) from idxstats. Excludes the
/// unmapped column, which legitimately differs for placed-but-unmapped reads.
fn idxstats_mapped(bam: &Path) -> String {
    idxstats(bam)
        .lines()
        .map(|line| {
            let mut f = line.split('\t');
            format!(
                "{}\t{}\t{}",
                f.next().unwrap_or_default(),
                f.next().unwrap_or_default(),
                f.next().unwrap_or_default()
            )
        })
        .collect::<Vec<_>>()
        .join("\n")
}

/// `samtools idxstats <bam>` — per-reference mapped/unmapped counts from the `.bai`.
fn idxstats(bam: &Path) -> String {
    let out = Command::new("samtools")
        .args(["idxstats", bam.to_str().unwrap()])
        .output()
        .expect("run samtools idxstats");
    assert!(
        out.status.success(),
        "samtools idxstats failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8_lossy(&out.stdout).into_owned()
}

/// The `.bai` fgumi writes must be equivalent to one samtools builds for the
/// *same* output bytes: identical region counts and identical idxstats.
///
/// This is the strong correctness check for the pooled indexing writer — it
/// validates the recovered virtual offsets against samtools (the reference
/// implementation) across a spill-forcing, multi-block, multi-threaded sort.
#[test]
#[ignore = "requires samtools"]
fn test_sort_write_index_matches_samtools_index() {
    if !samtools_available() {
        eprintln!("Skipping: samtools not available");
        return;
    }

    let temp_dir = TempDir::new().unwrap();
    // 15k reads/contig + unmapped → many BGZF blocks; -m 64K forces many spills
    // and a consolidation pass, so the fixed pooled indexed merge is exercised.
    let input_bam = create_large_test_bam(temp_dir.path(), 15_000);
    let sorted = temp_dir.path().join("sorted.bam");
    let fgumi_bai = temp_dir.path().join("sorted.bam.bai");

    let cmd = Sort::try_parse_from([
        OsStr::new("sort"),
        OsStr::new("-i"),
        input_bam.as_os_str(),
        OsStr::new("-o"),
        sorted.as_os_str(),
        OsStr::new("--order"),
        OsStr::new("coordinate"),
        OsStr::new("--write-index"),
        OsStr::new("-m"),
        OsStr::new("64K"),
        OsStr::new("-@"),
        OsStr::new("4"),
    ])
    .expect("parse sort args");
    cmd.execute("fgumi sort").expect("fgumi sort --write-index failed");
    assert!(sorted.exists() && fgumi_bai.exists(), "sorted BAM + fgumi .bai must exist");

    // Reference index: let samtools index the identical output bytes.
    let reference = temp_dir.path().join("reference.bam");
    std::fs::copy(&sorted, &reference).expect("copy sorted bam");
    let status = Command::new("samtools")
        .args(["index", reference.to_str().unwrap()])
        .status()
        .expect("run samtools index");
    assert!(status.success(), "samtools index failed");

    // Region counts via fgumi's .bai must equal those via samtools' .bai.
    let regions = [
        "chr1",
        "chr2",
        "chr1:1-100000",
        "chr1:500000-600000",
        "chr2:1-1000000",
        "chr2:250000-250500",
        "chr1:999500-1000000",
        "chr1:1-1000000",
    ];
    // MAPPED-read retrieval (`-F 4`) via fgumi's .bai must exactly match samtools'
    // .bai over the same bytes — this validates the recovered virtual offsets.
    // (Placed-but-unmapped reads are position-binned by samtools but not by fgumi;
    // see create_large_test_bam. That is a pre-existing indexing choice, not a
    // property of the pooled writer, so it is deliberately excluded here.)
    let mut total = 0i64;
    for r in regions {
        let via_fgumi = region_count(&sorted, r, &["-F", "4"]);
        let via_samtools = region_count(&reference, r, &["-F", "4"]);
        assert_eq!(
            via_fgumi, via_samtools,
            "region {r}: fgumi .bai mapped count {via_fgumi} != samtools .bai {via_samtools}"
        );
        total += via_fgumi;
    }
    assert!(total > 0, "expected non-empty mapped region queries");

    // Per-reference mapped counts must match exactly.
    assert_eq!(
        idxstats_mapped(&sorted),
        idxstats_mapped(&reference),
        "per-reference mapped counts via fgumi .bai must match samtools"
    );
}

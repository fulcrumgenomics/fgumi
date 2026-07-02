//! Integration tests for sort --write-index functionality.
//!
//! Verifies that BAM index generation produces a valid `.bai` that
//! samtools can query, across the coordinate sort code path.
//!
//! The input fixture is built in-process with `SamBuilder` + the noodles BAM
//! writer (no samtools), so the fgumi-only tests run by default. The samtools
//! cross-check tests are NO LONGER `#[ignore]`'d: they run by default and are
//! gated only by a runtime `samtools_available()` check, so they execute in CI
//! (samtools is installed in the workflow) and skip gracefully on a local dev
//! box without samtools.

use clap::Parser;
use fgumi_lib::commands::sort::Sort;
use fgumi_raw_bam::{RawRecord, SamBuilder, flags};
use noodles::bam;
use noodles::sam::Header;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use rstest::rstest;
use std::ffi::OsStr;
use std::path::Path;
use std::process::Command;
use tempfile::TempDir;

use crate::helpers::bam_generator::to_record_buf;

/// Check if samtools is available in PATH.
fn samtools_available() -> bool {
    which::which("samtools").is_ok()
}

/// Create a small unsorted multi-chromosome test BAM in-process (no samtools).
///
/// Six single-end records: five mapped across chr1/chr2 plus one unmapped, in a
/// deliberately unsorted emission order so the coordinate sort has work to do.
fn create_test_bam(dir: &Path) -> std::path::PathBuf {
    use bstr::BString;
    use noodles::sam::header::record::value::map::header::tag::Tag as HeaderTag;
    use noodles::sam::header::record::value::map::{Map, ReferenceSequence};
    use std::num::NonZeroUsize;

    let bam_path = dir.join("test_input.bam");

    let HeaderTag::Other(so_tag) = HeaderTag::from(*b"SO") else { unreachable!() };
    let header_map = Map::<noodles::sam::header::record::value::map::Header>::builder()
        .insert(so_tag, "unsorted")
        .build()
        .expect("valid header map");
    let mut hb = Header::builder().set_header(header_map);
    for chrom in ["chr1", "chr2"] {
        hb = hb.add_reference_sequence(
            BString::from(chrom),
            Map::<ReferenceSequence>::new(NonZeroUsize::new(10_000).expect("non-zero")),
        );
    }
    let header = hb.build();

    let seq = b"ACGTACGTAC";
    let qual = vec![40u8; seq.len()];
    let cigar = u32::try_from(seq.len()).expect("seq len fits u32") << 4; // 10M

    let mapped = |name: &str, ref_id: i32, pos: i32| -> RawRecord {
        let mut b = SamBuilder::new();
        b.read_name(name.as_bytes())
            .ref_id(ref_id)
            .pos(pos)
            .mapq(60)
            .flags(0)
            .cigar_ops(&[cigar])
            .sequence(seq)
            .qualities(&qual);
        b.build()
    };

    // Emit out of coordinate order on purpose.
    let mut records = vec![
        mapped("read4", 1, 99),  // chr2:100
        mapped("read2", 0, 199), // chr1:200
        mapped("read5", 1, 199), // chr2:200
        mapped("read1", 0, 99),  // chr1:100
        mapped("read3", 0, 299), // chr1:300
    ];
    // One unmapped read.
    let mut un = SamBuilder::new();
    un.read_name(b"read6").flags(flags::UNMAPPED).sequence(seq).qualities(&qual);
    records.push(un.build());

    let mut writer =
        bam::io::Writer::new(std::fs::File::create(&bam_path).expect("create input BAM"));
    writer.write_header(&header).expect("write header");
    for r in &records {
        writer.write_alignment_record(&header, &to_record_buf(r)).expect("write record");
    }
    writer.try_finish().expect("finish BAM");

    bam_path
}

/// Sort `input` into `output` (coordinate order) in-process via
/// `execute_sort_command` (the production standalone-sort entry point).
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
    fgumi_lib::commands::sort::execute_sort_command(&cmd, "fgumi sort").expect("fgumi sort failed");
}

/// Return the sorted read names returned by a samtools region query against a
/// (BAI-indexed) coordinate-sorted BAM.
fn region_query_names(bam_path: &Path, region: &str) -> Vec<String> {
    let output = Command::new("samtools")
        .arg("view")
        .arg(bam_path)
        .arg(region)
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

/// Verify that a BAM index allows region queries AND returns the correct slice.
/// Asserting record identities (not just a non-zero count) catches a BAI that
/// returns the wrong reads for the interval, which `samtools view -c > 0` would
/// miss. The `chr1:100-300` interval covers read1/read2/read3 in
/// `create_test_bam` (read4/read5 are on chr2).
fn verify_index_works(bam_path: &Path) -> bool {
    region_query_names(bam_path, "chr1:100-300")
        == vec!["read1".to_string(), "read2".to_string(), "read3".to_string()]
}

/// Test sort --write-index creates a valid index.
#[test]
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

    fgumi_lib::commands::sort::execute_sort_command(&cmd, "fgumi sort")
        .expect("fgumi sort --write-index failed");
    assert!(sorted_bam_path.exists(), "Output BAM not created");
    assert!(index_path.exists(), "Output BAI not created");

    // Verify index works
    assert!(verify_index_works(&sorted_bam_path), "Index does not work for region queries");
}

/// Test that --write-index with multi-threading still produces valid index.
#[test]
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

    fgumi_lib::commands::sort::execute_sort_command(&cmd, "fgumi sort")
        .expect("fgumi sort --write-index -@ 4 failed");
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
fn test_sort_write_index_byte_identical_to_off(#[case] threads: usize) {
    // No `samtools_available()` gate: this is a pure-fgumi regression (both runs
    // go through `run_coordinate_sort`), so it must always run — including on
    // minimal CI jobs without samtools.
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
        .arg("sort")
        .arg("-o")
        .arg(&sam_bam)
        .arg(&input_bam)
        .status()
        .expect("Failed to run samtools sort");
    assert!(status.success(), "samtools sort failed");
    let status = Command::new("samtools")
        .arg("index")
        .arg(&sam_bam)
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

    let result = fgumi_lib::commands::sort::execute_sort_command(&cmd, "fgumi sort");
    assert!(result.is_err(), "fgumi sort --order queryname --write-index should fail");

    let err_msg = format!("{:#}", result.unwrap_err());
    assert!(
        err_msg.contains("--write-index is only valid for coordinate sort"),
        "Expected error about coordinate sort, got: {err_msg}"
    );
}

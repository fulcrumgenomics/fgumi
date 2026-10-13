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
use std::process::Command;
use tempfile::TempDir;

use clap::Parser;
use fgumi_lib::commands::command::Command as FgumiCommand;
use fgumi_lib::commands::sort::Sort;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn fgumi_sort_in_process(args: &[OsString]) -> anyhow::Result<()> {
    // `clap::Parser::try_parse_from` accepts any iterator of `Into<OsString>`
    // values, so feeding `OsString` through avoids a UTF-8 round-trip that
    // would panic on a non-UTF-8 temp-dir path.
    let cmd = Sort::try_parse_from(args.iter().cloned()).expect("failed to parse sort args");
    cmd.execute("fgumi sort")
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
    let cmd_args: Vec<OsString> = vec![
        "sort".into(),
        "--verify".into(),
        "-i".into(),
        bam_path.as_os_str().to_owned(),
        "--order".into(),
        order.into(),
    ];
    fgumi_sort_in_process(&cmd_args).is_ok()
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

/// Collect a sorted multiset of per-record identities (`QNAME`, `FLAG`, `RNAME`, `POS`).
///
/// Sorting reorders records, so identity has to be compared order-independently.
/// A record count alone cannot distinguish "all records present" from "one record
/// dropped and another duplicated" — which is exactly the failure mode the sort
/// engine's spill/merge paths can produce under memory pressure. `RNAME` and `POS`
/// are included so a record whose alignment was corrupted in the spill/merge round
/// trip is caught too, not just one that vanished; fields beyond `POS` (`SEQ`,
/// `QUAL`, tags) are not compared.
fn record_identities(bam_path: &Path) -> Vec<String> {
    let output = Command::new("samtools")
        .args(["view", bam_path.to_str().unwrap()])
        .output()
        .expect("samtools view");
    assert!(output.status.success(), "samtools view failed for {}", bam_path.display());
    let mut identities: Vec<String> = String::from_utf8_lossy(&output.stdout)
        .lines()
        .map(|line| {
            let mut fields = line.split('\t');
            let qname = fields.next().unwrap_or_default();
            let flag = fields.next().unwrap_or_default();
            let rname = fields.next().unwrap_or_default();
            let pos = fields.next().unwrap_or_default();
            format!("{qname}\t{flag}\t{rname}\t{pos}")
        })
        .collect();
    identities.sort_unstable();
    identities
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
    for arg in extra_args {
        cmd_args.push((*arg).into());
    }
    fgumi_sort_in_process(&cmd_args).unwrap_or_else(|e| {
        panic!("fgumi sort failed for order={order} threads={threads} memory={max_memory}: {e:#}")
    });
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

/// At `--threads 1` the template split front hosts key extraction on the
/// `sort-coord` driver; a spilling sort there must emit every record, sorted,
/// in exactly the order the four-thread sort emits.
#[test]
#[ignore = "requires samtools"]
fn test_sort_template_coordinate_single_thread_with_spills() {
    if !samtools_available() {
        return;
    }
    let dir = TempDir::new().unwrap();
    let input = create_test_bam(dir.path(), 3000);
    let t1 = dir.path().join("sorted.t1.bam");
    let t4 = dir.path().join("sorted.t4.bam");

    sort_bam(&input, &t1, "template-coordinate", 1, "30K");
    sort_bam(&input, &t4, "template-coordinate", 4, "30K");
    assert!(verify_sorted(&t1, "template-coordinate"), "template-coordinate t1 failed");
    assert_eq!(record_identities(&input), record_identities(&t1), "t1 record multiset differs");
    let view = |path: &Path| {
        let out = Command::new("samtools").args(["view"]).arg(path).output().expect("samtools");
        assert!(out.status.success());
        out.stdout
    };
    assert!(view(&t1) == view(&t4), "t1 and t4 template-coordinate outputs differ");
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

/// Phase-1 input-serial reservation must not drop a worker's trailing blocks.
///
/// `try_read_input_blocks` reserves its batch's serial range while holding the
/// input-file lock. If serials were instead assigned after releasing the lock, a
/// worker preempted in that window would leave its already-read blocks
/// unaccounted for in `input_read_serial` while another worker read the empty EOF
/// batch and set `input_eof`; the done-check would fire early and those blocks
/// would be pushed after the merge finalized — a short output BAM with exit 0.
///
/// The race is timing-dependent, so this is a stress guard rather than a
/// deterministic reproducer: many input blocks and a high worker count maximize
/// contention on the input lock, and any lost block shows up as a count mismatch.
#[test]
#[ignore = "requires samtools"]
fn test_sort_multithreaded_does_not_truncate_output() {
    if !samtools_available() {
        return;
    }
    let dir = TempDir::new().unwrap();
    // Enough reads to span many BGZF input blocks, so multiple workers contend
    // on the input-file lock rather than one worker draining the whole file.
    let input = create_test_bam(dir.path(), 50_000);
    let input_count = count_records(&input);
    // Independent oracle: the exact multiset of records that must survive the sort.
    // Count alone would pass a regression that dropped one record and duplicated
    // another.
    let input_identities = record_identities(&input);

    for iteration in 0..3 {
        let output = dir.path().join(format!("sorted_stress_{iteration}.bam"));
        // Small memory forces spills; 8 threads maximizes lock contention.
        sort_bam(&input, &output, "coordinate", 8, "1M");

        assert_eq!(
            input_count,
            count_records(&output),
            "sorted output lost records on iteration {iteration} — Phase-1 serial \
             reservation regression (silent truncation)"
        );
        assert_eq!(
            input_identities,
            record_identities(&output),
            "sorted output is not the same multiset of records as the input on \
             iteration {iteration} (records dropped, duplicated, or their \
             QNAME/FLAG/RNAME/POS altered)"
        );
        assert!(verify_sorted(&output, "coordinate"), "iteration {iteration} not sorted");
    }
}

// ---------------------------------------------------------------------------
// Chain-level twins of the retired engine's tests (fgumi-sort `external.rs`).
// These need no samtools: their oracle is another `fgumi sort` run (or the
// error text), not an external tool.
// ---------------------------------------------------------------------------

use crate::helpers::cutover::decompressed_records_without_pg;
use fgumi_sam::SamBuilder;
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::{ReadGroup, ReferenceSequence};
use rstest::rstest;

/// Count a BAM's records natively, for the tests that run in every CI job
/// (`count_records` shells out to samtools, which only the samtools job has).
fn count_records_native(bam_path: &Path) -> u64 {
    let mut reader = noodles::bam::io::reader::Builder.build_from_path(bam_path).unwrap();
    let _ = reader.read_header().unwrap();
    let mut n = 0u64;
    for record in reader.records() {
        let _ = record.unwrap();
        n += 1;
    }
    n
}

/// `fgumi sort` in process, returning its error rather than panicking.
fn try_sort(
    input: &Path,
    output: &Path,
    order: &str,
    threads: usize,
    extra: &[&str],
) -> anyhow::Result<()> {
    let mut args: Vec<OsString> = vec![
        "sort".into(),
        "-i".into(),
        input.as_os_str().to_owned(),
        "-o".into(),
        output.as_os_str().to_owned(),
        "--order".into(),
        order.into(),
        "--threads".into(),
        threads.to_string().into(),
    ];
    args.extend(extra.iter().map(OsString::from));
    fgumi_sort_in_process(&args)
}

/// A header with one 200 Mb contig and the given `(id, library)` read groups.
fn header_with_read_groups(groups: &[(&str, &str)]) -> noodles::sam::Header {
    let mut builder = noodles::sam::Header::builder().add_reference_sequence(
        "chr1",
        Map::<ReferenceSequence>::new(std::num::NonZeroUsize::new(200_000_000).expect("nonzero")),
    );
    for (id, library) in groups {
        let rg = Map::<ReadGroup>::builder()
            .insert(noodles::sam::header::record::value::map::read_group::tag::LIBRARY, *library)
            .build()
            .expect("read group");
        builder = builder.add_read_group(*id, rg);
    }
    builder.build()
}

/// Twin of `lite_sort_hard_errors_when_{cb_appears,mi_appears,library_differs}`:
/// under `--key-types none` the first record fixes the dropped lanes, and a
/// later record with a different value fails the sort with the token that
/// re-includes the lane. The library case's second record names a read group
/// absent from the header (library ordinal 0 against the first's 1).
#[rstest]
#[case::cb("CB", "AAAACCCC", "--key-types cb")]
#[case::mi("MI", "7", "--key-types mi")]
#[case::library("RG", "rgUnknown", "--key-types library")]
fn a_dropped_lane_that_varies_fails_the_sort(
    #[case] tag: &str,
    #[case] value: &str,
    #[case] token: &str,
    #[values(1, 4)] threads: usize,
) {
    let dir = TempDir::new().unwrap();
    let mut builder = SamBuilder::new();
    builder.header = header_with_read_groups(&[("rgA", "LibAlpha")]);
    let first = builder.add_pair().name("read00000").start1(1).start2(101);
    let _ = if tag == "RG" { first.attr("RG", "rgA") } else { first }.build();
    let _ = builder.add_pair().name("read00001").start1(201).start2(301).attr(tag, value).build();
    let input = dir.path().join("in.bam");
    builder.write_bam(&input).unwrap();
    let err = try_sort(
        &input,
        &dir.path().join("out.bam"),
        "template-coordinate",
        threads,
        &["--key-types", "none"],
    )
    .expect_err("a varying dropped lane must fail the sort");
    assert!(format!("{err:#}").contains(token), "expected {token:?} guidance: {err:#}");
}

/// Twin of `auto_single_library_{multi_readgroup,unmapped_read}_sorts_without_false_violation`
/// (#375): one library under several read groups, and a fully unmapped pair in
/// that library, both sort under auto `--key-types` with no false library
/// violation — in memory and spilling.
#[rstest]
fn a_single_library_input_sorts_under_auto_without_a_false_violation(
    #[values("768M", "64K")] max_memory: &str,
    #[values(1, 4)] threads: usize,
) {
    let dir = TempDir::new().unwrap();
    let mut builder = SamBuilder::new();
    builder.header = header_with_read_groups(&[("rgA", "LibAlpha"), ("rgB", "LibAlpha")]);
    for i in 0..2_000usize {
        let rg = if i % 2 == 0 { "rgA" } else { "rgB" };
        let start = 1 + (i * 7919) % 100_000;
        let _ = builder
            .add_pair()
            .name(&format!("read{i:05}"))
            .start1(start)
            .start2(start + 100)
            .attr("RG", rg)
            .build();
    }
    let _ = builder.add_pair().name("unmapped").unmapped1().unmapped2().attr("RG", "rgA").build();
    let input = dir.path().join("in.bam");
    builder.write_bam(&input).unwrap();
    try_sort(
        &input,
        &dir.path().join("out.bam"),
        "template-coordinate",
        threads,
        &["-m", max_memory, "--memory-per-thread", "false"],
    )
    .expect("a single-library input must sort under auto without a false library violation");
}

/// Twin of `owned_template_sort_reports_a_violation_whose_read_name_overruns_its_body`:
/// a violating record whose `l_read_name` overruns its body fails the sort with
/// the dropped-lane violation, the record named by a placeholder, rather than
/// a panic.
#[test]
fn a_violation_on_a_record_with_an_overrunning_read_name_is_an_error() {
    use std::io::Write as _;
    let dir = TempDir::new().unwrap();
    let header = header_with_read_groups(&[]);
    let mut cb = b"CBZ".to_vec();
    cb.extend_from_slice(b"AAAA\0");
    let first = fgumi_raw_bam::testutil::make_bam_bytes(0, 100, 0x3, b"r0a", &[], 20, 0, 150, &cb);
    let mut bad =
        fgumi_raw_bam::testutil::make_bam_bytes(0, 101, 0x3, b"r1a", &[], 20, 0, 151, &[]);
    bad[8] = 200;
    let input = dir.path().join("in.bam");
    {
        let mut w = noodles::bgzf::io::Writer::new(std::fs::File::create(&input).unwrap());
        fgumi_bam_io::write_bam_header(&mut w, &header).unwrap();
        for r in [&first, &bad] {
            w.write_all(&u32::try_from(r.len()).unwrap().to_le_bytes()).unwrap();
            w.write_all(r).unwrap();
        }
        w.finish().unwrap();
    }
    for threads in [1, 4] {
        let result = std::panic::catch_unwind(|| {
            try_sort(
                &input,
                &dir.path().join(format!("out{threads}.bam")),
                "template-coordinate",
                threads,
                &["--key-types", "none"],
            )
        });
        let err = result.expect("the sort must not panic").expect_err("the record is invalid");
        let msg = format!("{err:#}");
        assert!(
            msg.contains("<malformed read name>") && msg.contains("--key-types cb"),
            "t{threads}: the dropped-lane violation, naming the malformed record: {msg}"
        );
    }
}

/// Template records whose optional lanes are constant: one CB on every record,
/// no MI, one library.
fn constant_lane_bam(dir: &Path, n: usize) -> std::path::PathBuf {
    let mut builder = SamBuilder::new();
    builder.header = header_with_read_groups(&[("rgA", "LibAlpha")]);
    for i in 0..n {
        let start = 1 + (i * 7919) % 50_000;
        let _ = builder
            .add_pair()
            .name(&format!("q{i:05}"))
            .start1(start)
            .start2(start + 150)
            .attr("RG", "rgA")
            .attr("CB", "AAAACCCC")
            .build();
    }
    let input = dir.join("constant_lanes.bam");
    builder.write_bam(&input).unwrap();
    input
}

/// Twin of `narrow_key_sort_byte_identical_to_full`, `..._when_spilling`,
/// `narrow_key_spill_merge_byte_identical_all_widths` and
/// `auto_drops_constant_library_lane_byte_identical`: when the dropped lanes
/// are constant, every `--key-types` width writes the same records as the full
/// key, in memory and across a spill and merge.
#[rstest]
fn every_key_width_matches_the_full_key_when_dropped_lanes_are_constant(
    #[values("none", "cb", "mi", "cb,mi", "auto")] key_types: &str,
    #[values("768M", "64K")] max_memory: &str,
    #[values(1, 4)] threads: usize,
) {
    let dir = TempDir::new().unwrap();
    let input = constant_lane_bam(dir.path(), 3_000);
    let full = dir.path().join("full.bam");
    try_sort(&input, &full, "template-coordinate", 1, &["--key-types", "full"]).unwrap();
    let narrow = dir.path().join("narrow.bam");
    let mut extra = vec!["-m", max_memory, "--memory-per-thread", "false"];
    if key_types != "auto" {
        extra.extend(["--key-types", key_types]);
    }
    try_sort(&input, &narrow, "template-coordinate", threads, &extra).unwrap();
    assert!(
        decompressed_records_without_pg(&narrow) == decompressed_records_without_pg(&full),
        "--key-types {key_types} at {max_memory} t{threads} differs from the full key"
    );
}

/// Twin of `test_sort_with_two_temp_dirs_matches_single_dir`: spilling across
/// two temp dirs writes what spilling into one does.
#[rstest]
fn two_temp_dirs_match_a_single_dir(
    #[values("coordinate", "queryname", "template-coordinate")] order: &str,
    #[values(1, 4)] threads: usize,
) {
    let dir = TempDir::new().unwrap();
    let (a, b) = (dir.path().join("a"), dir.path().join("b"));
    std::fs::create_dir_all(&a).unwrap();
    std::fs::create_dir_all(&b).unwrap();
    let input = constant_lane_bam(dir.path(), 3_000);
    let one = dir.path().join("one.bam");
    let two = dir.path().join("two.bam");
    let (a, b) = (a.to_str().unwrap(), b.to_str().unwrap());
    let base = ["-m", "64K", "--memory-per-thread", "false"];
    try_sort(&input, &one, order, threads, &[&base[..], &["-T", a]].concat()).unwrap();
    try_sort(&input, &two, order, threads, &[&base[..], &["-T", a, "-T", b]].concat()).unwrap();
    assert!(
        decompressed_records_without_pg(&one) == decompressed_records_without_pg(&two),
        "{order} t{threads}: two temp dirs changed the output"
    );
}

/// Twin of `template_sort_empty_input_writes_header_only_bam` and
/// `test_merge_chunks_*_with_only_empty_sources`: an input with no records
/// writes a valid, header-only BAM in every order.
#[rstest]
fn an_empty_input_writes_a_header_only_bam(
    #[values("coordinate", "queryname", "queryname::natural", "template-coordinate")] order: &str,
    #[values(1, 4)] threads: usize,
) {
    let dir = TempDir::new().unwrap();
    let input = dir.path().join("empty.bam");
    SamBuilder::new().write_bam(&input).unwrap();
    let output = dir.path().join("out.bam");
    try_sort(&input, &output, order, threads, &[]).unwrap();
    assert_eq!(count_records_native(&output), 0, "{order} t{threads}");
}

/// Twin of `test_sort_coordinate_with_index_in_memory_writes_sidecar`,
/// `test_sort_coordinate_with_index_spilled_preserves_records`,
/// `test_appending_run_with_write_index` and
/// `test_write_index_oversized_record_builds_loadable_bai`: `--write-index`
/// writes a loadable BAI and does not change the records, in memory and
/// spilling, with one record wider than a BGZF block.
#[rstest]
fn write_index_writes_a_loadable_bai_and_the_same_records(
    #[values("768M", "64K")] max_memory: &str,
    #[values(1, 4)] threads: usize,
) {
    let dir = TempDir::new().unwrap();
    let input = oversized_record_bam(dir.path());
    let plain = dir.path().join("plain.bam");
    let indexed = dir.path().join("indexed.bam");
    let base = ["-m", max_memory, "--memory-per-thread", "false"];
    try_sort(&input, &plain, "coordinate", threads, &base).unwrap();
    try_sort(&input, &indexed, "coordinate", threads, &[&base[..], &["--write-index"]].concat())
        .unwrap();
    let bai = indexed.with_extension("bam.bai");
    let index =
        noodles::bam::bai::fs::read(&bai).unwrap_or_else(|e| panic!("{}: {e}", bai.display()));
    assert_eq!(index.reference_sequences().len(), 1);
    assert!(decompressed_records_without_pg(&plain) == decompressed_records_without_pg(&indexed));
}

/// 2,000 pairs plus one pair whose reads are 70 kb, so each of its records is
/// wider than a BGZF block.
fn oversized_record_bam(dir: &Path) -> std::path::PathBuf {
    let mut builder = SamBuilder::new();
    builder.header = header_with_read_groups(&[]);
    for i in 0..2_000usize {
        let start = 1 + (i * 7919) % 150_000;
        let _ =
            builder.add_pair().name(&format!("q{i:05}")).start1(start).start2(start + 150).build();
    }
    let long = "ACGT".repeat(17_500);
    let _ = builder
        .add_pair()
        .name("wide")
        .start1(75_000)
        .start2(80_000)
        .bases1(&long)
        .bases2(&long)
        .build();
    let input = dir.join("oversized.bam");
    builder.write_bam(&input).unwrap();
    input
}

/// Twin of `test_oversized_record_survives_spill_and_pool_merge`: a record wider
/// than a BGZF block survives a spill and merge, so the spilled output matches
/// the in-memory one.
#[rstest]
fn a_record_wider_than_a_bgzf_block_survives_a_spill_and_merge(
    #[values("coordinate", "queryname", "template-coordinate")] order: &str,
    #[values(1, 4)] threads: usize,
) {
    let dir = TempDir::new().unwrap();
    let input = oversized_record_bam(dir.path());
    let in_memory = dir.path().join("mem.bam");
    let spilled = dir.path().join("spill.bam");
    try_sort(&input, &in_memory, order, threads, &["-m", "768M"]).unwrap();
    try_sort(&input, &spilled, order, threads, &["-m", "64K", "--memory-per-thread", "false"])
        .unwrap();
    assert_eq!(count_records_native(&spilled), 4_002);
    assert!(
        decompressed_records_without_pg(&in_memory) == decompressed_records_without_pg(&spilled),
        "{order} t{threads}: the spilled output differs from the in-memory one"
    );
}

/// Twin of `test_create_output_header_{coordinate,queryname,template_coordinate}`:
/// the output `@HD` declares the order the sort wrote.
#[rstest]
#[case::coordinate("coordinate", &[("SO", "coordinate")])]
#[case::queryname("queryname", &[("SO", "queryname")])]
#[case::template_coordinate(
    "template-coordinate",
    &[("SO", "unsorted"), ("GO", "query"), ("SS", "unsorted:template-coordinate")]
)]
fn the_output_header_declares_the_order(#[case] order: &str, #[case] want: &[(&str, &str)]) {
    let dir = TempDir::new().unwrap();
    let input = constant_lane_bam(dir.path(), 10);
    let output = dir.path().join("out.bam");
    try_sort(&input, &output, order, 1, &[]).unwrap();
    let mut reader = noodles::bam::io::reader::Builder.build_from_path(&output).unwrap();
    let header = reader.read_header().unwrap();
    let hd = header.header().expect("an @HD line");
    for (tag, value) in want {
        let key: [u8; 2] = tag.as_bytes().try_into().unwrap();
        let got = hd.other_fields().get(&key).map(std::string::ToString::to_string);
        assert_eq!(got.as_deref(), Some(*value), "{order}: @HD {tag}");
    }
}

/// Twin of `fgumi-sort`'s `test_output_compression.rs`: `--compression-level`
/// reaches the written BAM, and `--temp-compression` never does, in memory and
/// spilling, in every order.
#[rstest]
fn the_output_level_reaches_the_bam_and_the_temp_level_does_not(
    #[values("coordinate", "queryname", "template-coordinate")] order: &str,
    #[values("768M", "64K")] max_memory: &str,
    #[values(1, 4)] threads: usize,
) {
    let dir = TempDir::new().unwrap();
    let input = constant_lane_bam(dir.path(), 5_000);
    let run = |name: &str, level: &str, temp: &str| {
        let out = dir.path().join(format!("{name}.bam"));
        try_sort(
            &input,
            &out,
            order,
            threads,
            &[
                "-m",
                max_memory,
                "--memory-per-thread",
                "false",
                "--compression-level",
                level,
                "--temp-compression",
                temp,
                "--temp-codec",
                "bgzf",
            ],
        )
        .unwrap();
        (std::fs::metadata(&out).unwrap().len(), decompressed_records_without_pg(&out))
    };
    let (stored, stored_records) = run("l0_t1", "0", "1");
    let (packed, packed_records) = run("l9_t1", "9", "1");
    let (stored_hot_temp, _) = run("l0_t9", "0", "9");
    assert!(stored_records == packed_records, "the level changed the records");
    assert!(packed < stored / 2, "level 9 ({packed}) did not shrink level 0 ({stored})");
    assert_eq!(stored, stored_hot_temp, "the temp level reached the output");
}

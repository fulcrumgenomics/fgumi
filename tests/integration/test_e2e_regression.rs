//! End-to-end regression tests using simulate and compare commands.
//!
//! These tests validate that full pipelines produce deterministic, consistent output
//! by generating synthetic data with `simulate`, running pipeline commands, and
//! verifying outputs with `compare`.  No golden files are checked in — expected
//! output is generated fresh each run.

use clap::Parser;
use fgumi_lib::commands::command::Command as FgumiCommand;
use fgumi_lib::commands::compare::{CompareBams, CompareMismatch};
use fgumi_lib::commands::dedup::MarkDuplicates;
use fgumi_lib::commands::extract::Extract;
use fgumi_lib::commands::filter::Filter;
use fgumi_lib::commands::group::GroupReadsByUmi;
use fgumi_lib::commands::simplex::Simplex;
use fgumi_lib::commands::simulate::fastq_reads::FastqReads;
use fgumi_lib::commands::simulate::grouped_reads::GroupedReads;
use fgumi_lib::commands::sort::Sort;
use fgumi_lib::sam::SamTag;
use noodles::bam;
use std::ffi::{OsStr, OsString};
use std::fs::File;
use std::io::Write;
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

/// Build a command argument list from mixed string and path arguments.
macro_rules! args {
    ($($arg:expr),+ $(,)?) => {
        &[$( OsString::from($arg) ),+]
    };
}

// ---------------------------------------------------------------------------
// Helpers: simulate
// ---------------------------------------------------------------------------

/// Create a test reference FASTA file with a single chromosome.
fn create_test_reference(dir: &Path) -> PathBuf {
    let ref_path = dir.join("ref.fa");
    let mut f = std::fs::File::create(&ref_path).expect("failed to create ref FASTA");
    writeln!(f, ">chr1").unwrap();
    // 10 kb of repeating ACGT — large enough for any simulated insert size
    f.write_all(&b"ACGT".repeat(2500)).unwrap();
    writeln!(f).unwrap();
    f.flush().unwrap();
    ref_path
}

/// Generate grouped-reads BAM using simulate with deterministic seed.
fn simulate_grouped_reads(
    output: &Path,
    truth: &Path,
    reference: &Path,
    seed: u32,
    num_molecules: u32,
) {
    let num_molecules_str = num_molecules.to_string();
    let seed_str = seed.to_string();
    let cmd = GroupedReads::try_parse_from([
        OsStr::new("grouped-reads"),
        OsStr::new("-o"),
        output.as_os_str(),
        OsStr::new("--truth"),
        truth.as_os_str(),
        OsStr::new("--reference"),
        reference.as_os_str(),
        OsStr::new("--num-molecules"),
        OsStr::new(&num_molecules_str),
        OsStr::new("--seed"),
        OsStr::new(&seed_str),
        OsStr::new("--read-length"),
        OsStr::new("100"),
        OsStr::new("--umi-length"),
        OsStr::new("6"),
        OsStr::new("--min-family-size"),
        OsStr::new("2"),
    ])
    .expect("failed to parse grouped-reads args");
    cmd.execute("fgumi simulate grouped-reads").expect("simulate grouped-reads failed");
}

/// Generate FASTQ reads using simulate with deterministic seed.
fn simulate_fastq_reads(
    r1: &Path,
    r2: &Path,
    truth: &Path,
    reference: &Path,
    seed: u32,
    num_molecules: u32,
) {
    let num_molecules_str = num_molecules.to_string();
    let seed_str = seed.to_string();
    let cmd = FastqReads::try_parse_from([
        OsStr::new("fastq-reads"),
        OsStr::new("-1"),
        r1.as_os_str(),
        OsStr::new("-2"),
        r2.as_os_str(),
        OsStr::new("--truth"),
        truth.as_os_str(),
        OsStr::new("--reference"),
        reference.as_os_str(),
        OsStr::new("--num-molecules"),
        OsStr::new(&num_molecules_str),
        OsStr::new("--seed"),
        OsStr::new(&seed_str),
        OsStr::new("--read-length"),
        OsStr::new("100"),
        OsStr::new("--umi-length"),
        OsStr::new("6"),
        OsStr::new("--read-structure-r1"),
        OsStr::new("6M94T"),
        OsStr::new("--read-structure-r2"),
        OsStr::new("100T"),
        OsStr::new("--min-family-size"),
        OsStr::new("2"),
    ])
    .expect("failed to parse fastq-reads args");
    cmd.execute("fgumi simulate fastq-reads").expect("simulate fastq-reads failed");
}

// ---------------------------------------------------------------------------
// Helpers: pipeline steps
// ---------------------------------------------------------------------------

/// Run simplex consensus calling with single-threaded deterministic execution.
fn run_simplex(input: &Path, output: &Path, min_reads: u32) {
    let min_reads_str = min_reads.to_string();
    let cmd = Simplex::try_parse_from([
        OsStr::new("simplex"),
        OsStr::new("-i"),
        input.as_os_str(),
        OsStr::new("-o"),
        output.as_os_str(),
        OsStr::new("--threads"),
        OsStr::new("1"),
        OsStr::new("--min-reads"),
        OsStr::new(&min_reads_str),
    ])
    .expect("failed to parse simplex args");
    cmd.execute("fgumi simplex").expect("simplex failed");
}

/// Run filter on a consensus BAM.
fn run_filter(input: &Path, output: &Path, min_reads: u32, min_base_quality: u32) {
    let min_reads_str = min_reads.to_string();
    let min_base_quality_str = min_base_quality.to_string();
    let cmd = Filter::try_parse_from([
        OsStr::new("filter"),
        OsStr::new("-i"),
        input.as_os_str(),
        OsStr::new("-o"),
        output.as_os_str(),
        OsStr::new("--min-reads"),
        OsStr::new(&min_reads_str),
        OsStr::new("--min-base-quality"),
        OsStr::new(&min_base_quality_str),
    ])
    .expect("failed to parse filter args");
    cmd.execute("fgumi filter").expect("filter failed");
}

/// Run dedup on a grouped BAM.
fn run_dedup(input: &Path, output: &Path) {
    let cmd = MarkDuplicates::try_parse_from([
        OsStr::new("dedup"),
        OsStr::new("--input"),
        input.as_os_str(),
        OsStr::new("--output"),
        output.as_os_str(),
    ])
    .expect("failed to parse dedup args");
    cmd.execute("fgumi dedup").expect("dedup failed");
}

// ---------------------------------------------------------------------------
// Helpers: compare
// ---------------------------------------------------------------------------

/// Run `CompareBams::execute()` in-process. Returns true on match
/// (IDENTICAL or EQUIVALENT), false on `CompareMismatch` (DIFFER); panics on
/// any other anyhow error.
fn compare_bams_in_process(bam1: &Path, bam2: &Path, mode: &str) -> bool {
    let cmd = CompareBams::try_parse_from([
        OsStr::new("bams"),
        bam1.as_os_str(),
        bam2.as_os_str(),
        OsStr::new("--mode"),
        OsStr::new(mode),
    ])
    .expect("failed to parse compare bams args");
    match cmd.execute("fgumi compare bams") {
        Ok(()) => true,
        Err(e) if e.is::<CompareMismatch>() => false,
        Err(e) => panic!("compare bams hit unexpected error: {e:#}"),
    }
}

/// Assert that two BAM files are identical (or equivalent) according to the
/// given compare mode.
fn assert_bams_identical(bam1: &Path, bam2: &Path, mode: &str, context: &str) {
    assert!(compare_bams_in_process(bam1, bam2, mode), "{context}: compare bams reported DIFFER");
}

// ---------------------------------------------------------------------------
// Helpers: test setup
// ---------------------------------------------------------------------------

/// Create a `TempDir` and simulate grouped reads, returning (tmpdir, grouped bam path).
fn setup_grouped_reads(seed: u32, num_molecules: u32) -> (TempDir, PathBuf) {
    let tmp = TempDir::new().expect("failed to create temp dir");
    let reference = create_test_reference(tmp.path());
    let grouped = tmp.path().join("grouped.bam");
    let truth = tmp.path().join("truth.tsv");
    simulate_grouped_reads(&grouped, &truth, &reference, seed, num_molecules);
    (tmp, grouped)
}

// ---------------------------------------------------------------------------
// Determinism: simulate produces identical output with same seed
// ---------------------------------------------------------------------------

#[test]
fn test_simulate_grouped_reads_deterministic() {
    let tmp = TempDir::new().expect("failed to create temp dir");
    let reference = create_test_reference(tmp.path());
    let bam1 = tmp.path().join("grouped1.bam");
    let bam2 = tmp.path().join("grouped2.bam");
    let truth1 = tmp.path().join("truth1.tsv");
    let truth2 = tmp.path().join("truth2.tsv");

    simulate_grouped_reads(&bam1, &truth1, &reference, 42, 100);
    simulate_grouped_reads(&bam2, &truth2, &reference, 42, 100);

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
    let reference = create_test_reference(tmp.path());

    // Generate synthetic FASTQ
    let r1 = tmp.path().join("r1.fq.gz");
    let r2 = tmp.path().join("r2.fq.gz");
    let truth = tmp.path().join("truth.tsv");
    simulate_fastq_reads(&r1, &r2, &truth, &reference, 42, 200);

    // Run the full pipeline twice to verify determinism
    for suffix in ["a", "b"] {
        let extracted = tmp.path().join(format!("extracted_{suffix}.bam"));
        Extract::try_parse_from([
            OsStr::new("extract"),
            OsStr::new("--inputs"),
            r1.as_os_str(),
            r2.as_os_str(),
            OsStr::new("--output"),
            extracted.as_os_str(),
            OsStr::new("--read-structures"),
            OsStr::new("6M94T"),
            OsStr::new("100T"),
            OsStr::new("--sample"),
            OsStr::new("test_sample"),
            OsStr::new("--library"),
            OsStr::new("test_lib"),
        ])
        .expect("failed to parse extract args")
        .execute("fgumi extract")
        .expect("extract failed");

        // Extract emits SO:unsorted GO:query (no SS); `fgumi group` requires
        // template-coordinate sorted input, so put the sort step between them.
        let sorted = tmp.path().join(format!("sorted_{suffix}.bam"));
        Sort::try_parse_from([
            OsStr::new("sort"),
            OsStr::new("--input"),
            extracted.as_os_str(),
            OsStr::new("--output"),
            sorted.as_os_str(),
            OsStr::new("--order"),
            OsStr::new("template-coordinate"),
        ])
        .expect("failed to parse sort args")
        .execute("fgumi sort")
        .expect("sort failed");

        let grouped = tmp.path().join(format!("grouped_{suffix}.bam"));
        GroupReadsByUmi::try_parse_from([
            OsStr::new("group"),
            OsStr::new("--input"),
            sorted.as_os_str(),
            OsStr::new("--output"),
            grouped.as_os_str(),
            OsStr::new("--strategy"),
            OsStr::new("identity"),
            OsStr::new("--edits"),
            OsStr::new("0"),
        ])
        .expect("failed to parse group args")
        .execute("fgumi group")
        .expect("group failed");

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
// Strict template-coordinate input: extract → group must fail without sort
// ---------------------------------------------------------------------------

/// CLI-level regression test for the strict template-coordinate sort
/// requirement. `fgumi extract` writes `SO:unsorted GO:query` with no `SS`
/// tag; piping that directly into `fgumi group` must fail with an actionable
/// error message naming `SS:template-coordinate` and pointing at `fgumi sort`.
///
/// This test runs the actual `fgumi` binary (via `CARGO_BIN_EXE_fgumi`)
/// instead of calling `cmd.execute(...)` in-process so that it pins the
/// user-visible behavior — exit status and stderr — that downstream tooling
/// and CI pipelines depend on.
#[test]
fn test_group_rejects_extract_output_without_sort_step() {
    let tmp = TempDir::new().expect("failed to create temp dir");
    let reference = create_test_reference(tmp.path());

    let r1 = tmp.path().join("r1.fq.gz");
    let r2 = tmp.path().join("r2.fq.gz");
    let truth = tmp.path().join("truth.tsv");
    simulate_fastq_reads(&r1, &r2, &truth, &reference, 42, 50);

    let extracted = tmp.path().join("extracted.bam");
    Extract::try_parse_from([
        OsStr::new("extract"),
        OsStr::new("--inputs"),
        r1.as_os_str(),
        r2.as_os_str(),
        OsStr::new("--output"),
        extracted.as_os_str(),
        OsStr::new("--read-structures"),
        OsStr::new("6M94T"),
        OsStr::new("100T"),
        OsStr::new("--sample"),
        OsStr::new("test_sample"),
        OsStr::new("--library"),
        OsStr::new("test_lib"),
    ])
    .expect("failed to parse extract args")
    .execute("fgumi extract")
    .expect("extract failed");

    // Run `fgumi group` on the extract output directly — must fail.
    let grouped = tmp.path().join("grouped.bam");
    let output = fgumi(args![
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

    assert!(
        !output.status.success(),
        "fgumi group must reject extract output (non-template-coordinate sorted); \
         got exit status {:?}\nstdout: {}\nstderr: {}",
        output.status.code(),
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("template-coordinate"),
        "stderr should mention template-coordinate: {stderr}",
    );
    assert!(
        stderr.contains("SS:template-coordinate"),
        "stderr should name the SS:template-coordinate tag specifically: {stderr}",
    );
    assert!(
        stderr.contains("fgumi sort"),
        "stderr should point at `fgumi sort` as the remediation: {stderr}",
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
    let reference = create_test_reference(tmp.path());
    let bam1 = tmp.path().join("seed1.bam");
    let bam2 = tmp.path().join("seed2.bam");
    let truth1 = tmp.path().join("truth1.tsv");
    let truth2 = tmp.path().join("truth2.tsv");

    simulate_grouped_reads(&bam1, &truth1, &reference, 42, 100);
    simulate_grouped_reads(&bam2, &truth2, &reference, 99, 100);

    assert!(
        !compare_bams_in_process(&bam1, &bam2, "content"),
        "Expected compare bams to report content differences between distinct seeds",
    );
}

// ---------------------------------------------------------------------------
// Methylation pipeline: grouped-reads --methylation-mode em-seq ->
//                       simplex --methylation-mode em-seq --ref
// ---------------------------------------------------------------------------

#[test]
fn test_e2e_methylation_pipeline() {
    let tmp = TempDir::new().expect("failed to create temp dir");

    // Build a reference with CpG sites (ACGTCCGG pattern has CpG at positions 2-3 and 5-6)
    let ref_path = tmp.path().join("ref.fa");
    let mut f = std::fs::File::create(&ref_path).expect("failed to create ref FASTA");
    writeln!(f, ">chr1").unwrap();
    f.write_all(&b"ACGTCCGG".repeat(1250)).unwrap(); // 10 kb
    writeln!(f).unwrap();
    f.flush().unwrap();

    // Generate grouped reads with EM-Seq methylation
    let grouped = tmp.path().join("grouped.bam");
    let truth = tmp.path().join("truth.tsv");
    GroupedReads::try_parse_from([
        OsStr::new("grouped-reads"),
        OsStr::new("-o"),
        grouped.as_os_str(),
        OsStr::new("--truth"),
        truth.as_os_str(),
        OsStr::new("--reference"),
        ref_path.as_os_str(),
        OsStr::new("--num-molecules"),
        OsStr::new("50"),
        OsStr::new("--seed"),
        OsStr::new("42"),
        OsStr::new("--read-length"),
        OsStr::new("100"),
        OsStr::new("--umi-length"),
        OsStr::new("6"),
        OsStr::new("--min-family-size"),
        OsStr::new("3"),
        OsStr::new("--methylation-mode"),
        OsStr::new("em-seq"),
    ])
    .expect("failed to parse grouped-reads args")
    .execute("fgumi simulate grouped-reads")
    .expect("simulate grouped-reads failed");

    // Run simplex consensus with methylation mode
    let simplex = tmp.path().join("simplex.bam");
    Simplex::try_parse_from([
        OsStr::new("simplex"),
        OsStr::new("-i"),
        grouped.as_os_str(),
        OsStr::new("-o"),
        simplex.as_os_str(),
        OsStr::new("--threads"),
        OsStr::new("1"),
        OsStr::new("--min-reads"),
        OsStr::new("1"),
        OsStr::new("--methylation-mode"),
        OsStr::new("em-seq"),
        OsStr::new("--ref"),
        ref_path.as_os_str(),
    ])
    .expect("failed to parse simplex args")
    .execute("fgumi simplex")
    .expect("simplex failed");

    // Verify the simplex output actually carries methylation emission.
    // The simplex BAM should contain at least one consensus record with the
    // MM/ML methylation tag pair and the cu/ct TAPS-specific count arrays.
    assert_simplex_has_methylation_tags(&simplex);
}

/// Assert that a simplex BAM contains non-empty MM/ML and cu/ct tags on at
/// least one record.
///
/// MM is a string tag (hex 'MM'), ML is a byte-array tag (hex 'ML').
/// cu and ct are i16-array tags added by the EM-Seq/TAPS methylation caller.
#[allow(clippy::similar_names)] // MM/ML and cu/ct are the natural tag names.
fn assert_simplex_has_methylation_tags(bam_path: &Path) {
    let file = File::open(bam_path).expect("failed to open simplex BAM");
    let mut reader = bam::io::Reader::new(file);
    let _header = reader.read_header().expect("failed to read BAM header");

    let mm_tag = SamTag::MM.to_noodles_tag();
    let ml_tag = SamTag::ML.to_noodles_tag();
    let cu_tag = SamTag::CU.to_noodles_tag();
    let ct_tag = SamTag::CT.to_noodles_tag();

    let mut total = 0usize;
    let mut with_mm = 0usize;
    let mut with_ml = 0usize;
    let mut with_cu = 0usize;
    let mut with_ct = 0usize;
    for result in reader.records() {
        let record = result.expect("failed to read simplex record");
        let data = record.data();
        total += 1;
        if data.get(&mm_tag).is_some() {
            with_mm += 1;
        }
        if data.get(&ml_tag).is_some() {
            with_ml += 1;
        }
        if data.get(&cu_tag).is_some() {
            with_cu += 1;
        }
        if data.get(&ct_tag).is_some() {
            with_ct += 1;
        }
    }

    assert!(total > 0, "Simplex BAM should contain at least one record");
    assert!(with_mm > 0, "Expected at least one record with an MM tag; got {total} records");
    assert!(with_ml > 0, "Expected at least one record with an ML tag; got {total} records");
    assert!(with_cu > 0, "Expected at least one record with a cu tag; got {total} records");
    assert!(with_ct > 0, "Expected at least one record with a ct tag; got {total} records");
}

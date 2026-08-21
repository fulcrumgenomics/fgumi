//! End-to-end regression tests using simulate and compare commands.
//!
//! These tests validate that full pipelines produce deterministic, consistent output
//! by generating synthetic data with `simulate`, running pipeline commands, and
//! verifying outputs with `compare`.  No golden files are checked in — expected
//! output is generated fresh each run.

use clap::Parser;
use fgumi_lib::commands::command::Command as FgumiCommand;
use fgumi_lib::commands::compare::{CompareBams, CompareMismatch, compare_headers};
use fgumi_lib::commands::dedup::MarkDuplicates;
use fgumi_lib::commands::extract::Extract;
use fgumi_lib::commands::filter::Filter;
use fgumi_lib::commands::group::GroupReadsByUmi;
use fgumi_lib::commands::simplex::Simplex;
use fgumi_lib::commands::simulate::consensus_reads::ConsensusReads;
use fgumi_lib::commands::simulate::correct_reads::CorrectReads;
use fgumi_lib::commands::simulate::fastq_reads::FastqReads;
use fgumi_lib::commands::simulate::grouped_reads::GroupedReads;
use fgumi_lib::commands::sort::Sort;
use fgumi_lib::sam::SamTag;
use noodles::bam;
use rstest::rstest;
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

/// Sort a BAM into genuine template-coordinate order.
///
/// `simulate grouped-reads` declares `SS:template-coordinate` in its output header but does
/// not itself guarantee the emitted records are actually in that order (see the `NOTE` in
/// `grouped_reads.rs`'s molecule-sort-key comment — this is a documented, tracked gap fixed
/// wholesale by #576, not something this test suite should paper over by weakening the
/// `@HD` claim). Real pipelines always sort before a template-coordinate-order-dependent
/// step, so tests that content-compare `simulate grouped-reads` output (or anything
/// downstream that preserves record order, like `dedup`) must do the same or `fgumi compare
/// bams`'s universal sort-order verification precondition correctly rejects the merely
/// declared-but-not-actual order before a single record is paired.
fn sort_template_coordinate(input: &Path, output: &Path) {
    Sort::try_parse_from([
        OsStr::new("sort"),
        OsStr::new("--input"),
        input.as_os_str(),
        OsStr::new("--output"),
        output.as_os_str(),
        OsStr::new("--order"),
        OsStr::new("template-coordinate"),
    ])
    .expect("failed to parse sort args")
    .execute("fgumi sort")
    .expect("sort failed");
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

// ---------------------------------------------------------------------------
// Helpers: record-byte-exact comparison
// ---------------------------------------------------------------------------
//
// `compare_bams_in_process(..., "content")` (used directly by
// `test_different_seeds_produce_different_output` below) normalizes tag order, integer tag
// width, the `tc` tag, and `@PG`/`@CO` (see `raw_compare::content_key_exact`) -- by design,
// so genuinely equivalent output from different tools/writers still matches. But that same
// normalization means a *same-seed, same-binary* determinism test using "content" mode
// cannot detect writer nondeterminism confined to exactly those normalized dimensions (e.g.
// two runs emitting the same tags in a different order, or a different integer tag width for
// the same value). The helpers below instead compare every record's raw on-disk bytes
// exactly, so a same-seed run is held to a strictly stronger bar: not just semantically
// equivalent, but byte-identical -- which is why every same-seed determinism test in this
// file asserts via [`assert_bams_record_byte_identical`] alone, not the weaker
// content-mode compare.

/// Reads every record from `path`, in on-disk order, as its full raw BAM bytes via the same
/// non-noodles reader the sort/compare engines use internally. Nothing in the record body is
/// masked — including the `bin` field: for a *same-seed, same-binary* determinism check `bin`
/// is a deterministic pure function of the record's already-compared position + CIGAR span
/// (`fgumi group`/`simulate` compute it via `region_to_bin`; pass-through commands preserve
/// the source bytes), so it can never flag a false difference, and keeping it holds this
/// helper to its "byte-identical" contract — catching any future writer that emits `bin`
/// nondeterministically. (`compare`'s `content` mode deliberately ignores `bin` as a non-SAM
/// index artifact; this stricter byte-exact backstop is a different, complementary check.)
///
/// The header is intentionally never read past (`skip_header`): `@PG` carries this process's
/// own invocation (e.g. absolute temp-dir paths), which differs run-to-run even for
/// genuinely deterministic record output, so a byte-exact header comparison would produce
/// false positives unrelated to the writer determinism this helper exists to check.
fn read_all_record_bytes(path: &Path) -> Vec<Vec<u8>> {
    let file = File::open(path).expect("open BAM for record-byte-exact read");
    let mut reader = fgumi_sort::RawBamRecordReader::new(file).expect("open raw BAM reader");
    reader.skip_header().expect("skip header");
    let mut records = Vec::new();
    while let Some(record) = reader.next_record().expect("read raw BAM record") {
        records.push(record.as_ref().to_vec());
    }
    records
}

/// The parsed SAM header of a BAM file, for the header half of
/// [`assert_bams_record_byte_identical`]'s comparison.
fn read_bam_header(path: &Path) -> noodles::sam::Header {
    let file = File::open(path).expect("open BAM for header read");
    let mut reader = bam::io::Reader::new(file);
    reader.read_header().expect("read BAM header")
}

/// The raw `@HD SS` sub-sort tag, if present. `compare_headers` deliberately *normalizes* `SS`
/// (bare-vs-prefixed) and compares `SO`/`GO` only semantically, so a changed queryname/
/// template-coordinate subsort would slip past it — but for a same-seed determinism check the
/// `@HD` must be byte-stable, so this helper lets the caller compare `SS` exactly.
fn hd_ss(header: &noodles::sam::Header) -> Option<Vec<u8>> {
    header.header().and_then(|hd| hd.other_fields().get(b"SS").map(|ss| ss.to_vec()))
}

/// Assert that two BAM files carry exactly the same records, in the same on-disk order,
/// down to the raw byte (raw *header bytes* excluded — see [`read_all_record_bytes`] — but the
/// behaviorally relevant header fields `@SQ`/`@RG`/`SO`/`GO` are still compared via
/// [`compare_headers`], which normalizes the per-invocation `@PG`/`@CO`, plus an exact `@HD SS`
/// sub-sort check that `compare_headers` would otherwise normalize away). Unlike a
/// `"content"`-mode compare (see `compare_bams_in_process`), the record comparison does NOT
/// normalize tag order, integer tag width, or the `tc` tag, so it catches same-seed writer
/// nondeterminism confined to exactly those dimensions.
fn assert_bams_record_byte_identical(bam1: &Path, bam2: &Path, context: &str) {
    // The raw *header bytes* are excluded from the byte-exact record comparison (see
    // `read_all_record_bytes`: `@PG` carries run-specific temp paths that differ run-to-run
    // even for deterministic output). But skipping the header entirely would let a divergent
    // `@SQ`/`@RG`/`SO`/`GO` slip past these "identical BAM" tests, so compare the behaviorally
    // relevant header fields separately via `compare_headers` (which normalizes `@PG`/`@CO`).
    let (header1, header2) = (read_bam_header(bam1), read_bam_header(bam2));
    if let Some(diffs) = compare_headers(&header1, &header2) {
        panic!(
            "{context}: BAM headers differ on behaviorally-relevant fields (@SQ/@RG/SO/GO): {}",
            diffs.join("; ")
        );
    }
    // `compare_headers` normalizes `@HD SS`, so compare it exactly here — a same-seed run must
    // emit the identical sub-sort, and a changed subsort semantics must not read as "identical".
    assert_eq!(
        hd_ss(&header1),
        hd_ss(&header2),
        "{context}: @HD SS sub-sort tag differs between same-seed outputs"
    );

    let records1 = read_all_record_bytes(bam1);
    let records2 = read_all_record_bytes(bam2);
    assert_eq!(
        records1.len(),
        records2.len(),
        "{context}: record count differs byte-exact ({} vs {})",
        records1.len(),
        records2.len()
    );
    for (i, (r1, r2)) in records1.iter().zip(records2.iter()).enumerate() {
        assert_eq!(
            r1, r2,
            "{context}: record {i} is not byte-identical -- \
             same-seed run produced different bytes for the same logical record"
        );
    }
}

// ---------------------------------------------------------------------------
// Helpers: test setup
// ---------------------------------------------------------------------------

/// Create a `TempDir` and simulate grouped reads, returning (tmpdir, grouped bam path).
///
/// The returned BAM is sorted into genuine template-coordinate order (see
/// [`sort_template_coordinate`]) so it is safe to feed directly into template-coordinate-
/// order-dependent commands (`dedup`) and to content-compare downstream output against.
fn setup_grouped_reads(seed: u32, num_molecules: u32) -> (TempDir, PathBuf) {
    let tmp = TempDir::new().expect("failed to create temp dir");
    let reference = create_test_reference(tmp.path());
    let unsorted = tmp.path().join("grouped_unsorted.bam");
    let truth = tmp.path().join("truth.tsv");
    simulate_grouped_reads(&unsorted, &truth, &reference, seed, num_molecules);
    let grouped = tmp.path().join("grouped.bam");
    sort_template_coordinate(&unsorted, &grouped);
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

    // Compare the raw simulate output directly, WITHOUT an intervening sort. This is a
    // determinism test, and `assert_bams_record_byte_identical` reads records in raw on-disk
    // order with no sort-order precondition — so sorting first would only canonicalize record
    // order and thereby mask any run-to-run *ordering* nondeterminism in `simulate`. The
    // single-threaded (`--threads 1`, the helper's default) sorted-molecule-id streaming
    // generation emits records in a deterministic order for a given seed, so the same-seed
    // outputs are byte-identical as-is; holding them to that stricter bar is the point.
    assert_bams_record_byte_identical(
        &bam1,
        &bam2,
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

    assert_bams_record_byte_identical(
        &simplex1,
        &simplex2,
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

    assert_bams_record_byte_identical(
        &filtered1,
        &filtered2,
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
    assert_bams_record_byte_identical(
        &tmp.path().join("filtered_a.bam"),
        &tmp.path().join("filtered_b.bam"),
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

    assert_bams_record_byte_identical(
        &dedup1,
        &dedup2,
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

    // See `sort_template_coordinate`'s doc comment: sort before content-comparing so the
    // universal sort-order verification precondition sees genuinely ordered input rather
    // than rejecting the pair before content is ever compared.
    let sorted1 = tmp.path().join("seed1_sorted.bam");
    let sorted2 = tmp.path().join("seed2_sorted.bam");
    sort_template_coordinate(&bam1, &sorted1);
    sort_template_coordinate(&bam2, &sorted2);

    assert!(
        !compare_bams_in_process(&sorted1, &sorted2, "content"),
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

/// A writer that dies must report *its* error, not the send failure it causes.
///
/// The writer thread owns the output files; the generator only learns something
/// went wrong because its channel closed. The command used to return on that send
/// failure without joining the writer, so it reported "Failed to send record to
/// writer: sending on a disconnected channel" — the symptom — and discarded the
/// writer's own error along with the unjoined handle. It also left the writer
/// running over its output while the caller unwound.
///
/// A directory where R1 should go makes the writer fail at its first `File::create`,
/// which is the cheapest way to reach that path deterministically.
#[test]
fn simulate_fastq_reads_reports_the_writer_error_not_the_send_failure() {
    let dir = TempDir::new().expect("create temp dir");
    let reference = create_test_reference(dir.path());

    // R1's path is a directory, so the writer cannot create it.
    let r1 = dir.path().join("r1.fastq.gz");
    std::fs::create_dir(&r1).expect("create the blocking directory");

    let cmd = FastqReads::try_parse_from([
        OsStr::new("fastq-reads"),
        OsStr::new("-1"),
        r1.as_os_str(),
        OsStr::new("-2"),
        dir.path().join("r2.fastq.gz").as_os_str(),
        OsStr::new("--truth"),
        dir.path().join("truth.txt").as_os_str(),
        OsStr::new("--reference"),
        reference.as_os_str(),
        OsStr::new("--num-molecules"),
        OsStr::new("200"),
    ])
    .expect("failed to parse fastq-reads args");

    let err = cmd.execute("fgumi simulate fastq-reads").expect_err("the writer cannot succeed");
    let rendered = format!("{err:#}");
    assert!(
        rendered.contains("Failed to create"),
        "the writer's own error must be reported, got: {rendered}"
    );
    assert!(
        !rendered.contains("Failed to send record to writer"),
        "the send failure is the symptom, not the cause: {rendered}"
    );
}

/// The BAM-writing siblings of the test above must report the writer's error too.
///
/// `consensus_reads`, `correct_reads` and `fastq_reads` each spawn a writer thread
/// and each had the same early return, so fixing one and testing only that one is
/// how the other two would quietly regress. A directory where the output BAM
/// should go makes the writer fail at `create_raw_bam_writer`, which reaches the
/// same path as the FASTQ case.
#[rstest]
#[case::consensus_reads(false)]
#[case::correct_reads(true)]
fn simulate_bam_commands_report_the_writer_error_not_the_send_failure(#[case] correct: bool) {
    let dir = TempDir::new().expect("create temp dir");
    let reference = create_test_reference(dir.path());

    // The output BAM's path is a directory, so the writer cannot create it.
    let output = dir.path().join("out.bam");
    std::fs::create_dir(&output).expect("create the blocking directory");

    let err = if correct {
        let includelist = dir.path().join("umis.txt");
        std::fs::write(&includelist, "ACGTACGT\nTTTTGGGG\n").expect("write include list");
        let cmd = CorrectReads::try_parse_from([
            OsStr::new("correct-reads"),
            OsStr::new("-o"),
            output.as_os_str(),
            OsStr::new("--includelist"),
            includelist.as_os_str(),
            OsStr::new("--truth"),
            dir.path().join("truth.txt").as_os_str(),
            OsStr::new("--num-reads"),
            OsStr::new("200"),
        ])
        .expect("failed to parse correct-reads args");
        cmd.execute("fgumi simulate correct-reads").expect_err("the writer cannot succeed")
    } else {
        let cmd = ConsensusReads::try_parse_from([
            OsStr::new("consensus-reads"),
            OsStr::new("-o"),
            output.as_os_str(),
            OsStr::new("--reference"),
            reference.as_os_str(),
            OsStr::new("--num-reads"),
            OsStr::new("200"),
        ])
        .expect("failed to parse consensus-reads args");
        cmd.execute("fgumi simulate consensus-reads").expect_err("the writer cannot succeed")
    };

    let rendered = format!("{err:#}");
    assert!(
        rendered.contains("Failed to create output BAM"),
        "the writer's own error must be reported, got: {rendered}"
    );
    assert!(
        !rendered.contains("Failed to send record to writer"),
        "the send failure is the symptom, not the cause: {rendered}"
    );
}

// ---------------------------------------------------------------------------
// Regression (#687): `simulate grouped-reads` output must be accepted by the
// grouped-input engine it exists to feed.
// ---------------------------------------------------------------------------

/// Generate a `simulate grouped-reads` BAM, optionally in duplex mode.
///
/// Mirrors [`simulate_grouped_reads`] but adds `--duplex` when requested, so the
/// regression test can cover both the simplex (`<id>`) and duplex (`<id>/A`,
/// `<id>/B`) MI-tag shapes.
fn simulate_grouped_reads_mode(
    output: &Path,
    truth: &Path,
    reference: &Path,
    seed: u32,
    num_molecules: u32,
    duplex: bool,
) {
    let num_molecules_str = num_molecules.to_string();
    let seed_str = seed.to_string();
    let mut argv: Vec<&OsStr> = vec![
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
    ];
    if duplex {
        argv.push(OsStr::new("--duplex"));
    }
    let cmd = GroupedReads::try_parse_from(argv).expect("failed to parse grouped-reads args");
    cmd.execute("fgumi simulate grouped-reads").expect("simulate grouped-reads failed");
}

/// Run `fgumi compare bams <a> <b> --command group` in-process, returning the raw
/// `execute` result so the caller can distinguish a clean match (`Ok`) from the
/// grouped-input precondition rejection (`Err`) that #687 is about.
fn compare_bams_group_command(bam1: &Path, bam2: &Path) -> anyhow::Result<()> {
    let cmd = CompareBams::try_parse_from([
        OsStr::new("bams"),
        bam1.as_os_str(),
        bam2.as_os_str(),
        OsStr::new("--command"),
        OsStr::new("group"),
    ])
    .expect("failed to parse compare bams args");
    cmd.execute("fgumi compare bams")
}

/// `simulate grouped-reads` exists to produce grouping test input, so its output
/// must round-trip through the engine that consumes grouped BAMs: comparing the
/// file against a byte-identical copy of itself must report a match, not fail the
/// grouped-input precondition.
///
/// Regression test for #687: MI ids were minted as a pre-sort counter and emerged
/// non-monotonic in the template-coordinate-sorted file, so
/// `compare bams --command group` rejected `simulate`'s own output with
/// "input is not grouped: MI base N is not greater than a previously seen base M".
#[rstest]
#[case::simplex(false)]
#[case::duplex(true)]
fn simulate_grouped_reads_roundtrips_through_compare_group(#[case] duplex: bool) {
    let tmp = TempDir::new().expect("failed to create temp dir");
    let reference = create_test_reference(tmp.path());
    let bam = tmp.path().join("grouped.bam");
    let truth = tmp.path().join("truth.tsv");
    simulate_grouped_reads_mode(&bam, &truth, &reference, 42, 500, duplex);

    // A byte-identical copy: comparing a file against itself cannot legitimately
    // DIFFER, so any failure here is the grouped-input precondition rejecting the
    // simulated input, not a comparison result.
    let copy = tmp.path().join("grouped.copy.bam");
    std::fs::copy(&bam, &copy).expect("copy grouped bam");

    compare_bams_group_command(&bam, &copy).unwrap_or_else(|e| {
        panic!(
            "`compare bams --command group` rejected simulate grouped-reads output \
             (duplex={duplex}): {e:#}"
        )
    });
}

/// The truth TSV must stay consistent with the renumbered BAM. After #687's
/// post-sort MI renumber, a consumer that joins truth rows to BAM records by read
/// name must still see the same MI: this asserts every BAM record's MI equals its
/// truth row's `mi_tag`, the truth `molecule_id` equals the MI base, and the
/// renumbered bases are exactly `1..=N` in file order (monotonic and contiguous).
#[rstest]
#[case::simplex(false)]
#[case::duplex(true)]
fn simulate_grouped_reads_truth_matches_renumbered_bam(#[case] duplex: bool) {
    let tmp = TempDir::new().expect("failed to create temp dir");
    let reference = create_test_reference(tmp.path());
    let bam = tmp.path().join("grouped.bam");
    let truth = tmp.path().join("truth.tsv");
    simulate_grouped_reads_mode(&bam, &truth, &reference, 7, 300, duplex);

    // truth: read_name -> (molecule_id, mi_tag)
    let truth_txt = std::fs::read_to_string(&truth).expect("read truth");
    let mut truth_by_name: std::collections::HashMap<String, (String, String)> =
        std::collections::HashMap::new();
    let mut lines = truth_txt.lines();
    lines.next().expect("truth header"); // skip header
    for line in lines {
        let cols: Vec<&str> = line.split('\t').collect();
        assert!(cols.len() > 3, "malformed truth row: {line:?}");
        // Reject duplicate read names: each pair has a unique name, and a
        // silently-overwritten row would let a stale/duplicate entry pass the
        // relation check below.
        let prior =
            truth_by_name.insert(cols[0].to_string(), (cols[2].to_string(), cols[3].to_string()));
        assert!(prior.is_none(), "duplicate truth read name: {}", cols[0]);
    }

    let mut prev_base: Option<u64> = None;
    let mut distinct_bases: Vec<u64> = Vec::new();
    let mut bam_names: std::collections::HashSet<String> = std::collections::HashSet::new();
    for bytes in read_all_record_bytes(&bam) {
        let view = fgumi_raw_bam::RawRecordView::new(&bytes);
        let qname: String =
            String::from_utf8(view.read_name().iter().copied().take_while(|&b| b != 0).collect())
                .expect("read name is UTF-8");
        let mi_bytes = fgumi_raw_bam::find_string_tag_in_record(&bytes, SamTag::MI)
            .expect("record has MI tag");
        let mi = std::str::from_utf8(mi_bytes).expect("MI is UTF-8").to_string();

        let (truth_mol_id, truth_mi) = truth_by_name
            .get(&qname)
            .unwrap_or_else(|| panic!("BAM read {qname} absent from truth"));
        assert_eq!(&mi, truth_mi, "BAM MI != truth mi_tag for {qname} (duplex={duplex})");

        let base: u64 = mi.split('/').next().unwrap().parse().expect("MI base is an integer");
        assert_eq!(
            &base.to_string(),
            truth_mol_id,
            "truth molecule_id != MI base for {qname} (duplex={duplex})"
        );

        if prev_base != Some(base) {
            if let Some(p) = prev_base {
                assert!(base > p, "MI base not strictly increasing in file order: {p} then {base}");
            }
            distinct_bases.push(base);
            prev_base = Some(base);
        }
        bam_names.insert(qname);
    }

    let n = distinct_bases.len() as u64;
    assert_eq!(
        distinct_bases,
        (1..=n).collect::<Vec<_>>(),
        "renumbered bases must be exactly 1..=N in file order (duplex={duplex})"
    );

    // Check the relation in both directions: every truth row is realized in the
    // BAM and vice versa, so a stale/extra truth row (or a dropped BAM read) fails
    // rather than slipping through the per-record lookup above.
    let truth_names: std::collections::HashSet<String> = truth_by_name.into_keys().collect();
    assert_eq!(
        bam_names, truth_names,
        "BAM read-name set != truth read-name set (duplex={duplex})"
    );
}

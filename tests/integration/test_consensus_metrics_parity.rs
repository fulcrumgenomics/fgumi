//! Correctness anchor: inline consensus metrics (`runall --simplex::metrics`
//! / `--duplex::metrics` / `--codec::metrics`, and the standalone unified
//! consensus command's own `--metrics`) must match the separate-pass
//! `simplex-metrics`/`duplex-metrics` commands' numeric output exactly.
//!
//! Three independently-computed paths are compared for every case below:
//!
//! - **Ground truth**: `fgumi group` then `fgumi simplex-metrics` /
//!   `fgumi duplex-metrics` on the grouped BAM (a fully separate pass).
//! - **T1 (fused)**: `fgumi runall --start-from group --consensus <mode>
//!   --<mode>::metrics=<prefix>` on the same raw, ungrouped input — trusts
//!   `GroupByPosition`'s pre-formed coordinate-group boundaries directly.
//! - **T2 (standalone)**: the standalone unified consensus command
//!   (`simplex`/`duplex`/`codec`) with `--metrics=<prefix>`, fed the *same*
//!   grouped BAM ground truth used above — re-derives coordinate-group
//!   boundaries itself via `ReadInfoKey`-equality (`split_into_runs`/
//!   `classify_batch_runs`/`reassemble_boundary`), independently of
//!   `GroupByPosition`.
//!
//! T1 and T2 share the underlying accumulation math (Tasks 5/7's reducer
//! functions), so any divergence from ground truth — or between T1 and T2 —
//! traces to a boundary-detection bug in one of the two independently
//! implemented mechanisms, not to the arithmetic itself. Comparison is exact
//! string equality of each output file's full text: every path drives the
//! *same* `write_metrics_auto` TSV writer over `f64`/`usize` fields computed
//! from the *same* input, so a genuine parity holds bit-for-bit, and any
//! difference (formatting or value) is real evidence of a bug, not a
//! tolerance artifact — matching this repo's existing metrics-parity
//! convention (`test_group_metrics_fused_parity.rs`).
//!
//! `--group::edits 0` / `-e 0` is used throughout the extended matrix (not
//! Step 1's original single-family case, which mirrors the task brief
//! verbatim) so UMI-to-MI assignment is exact-match only: many-molecule
//! shapes use programmatically distinct UMIs, and `edits 0` removes any risk
//! of the adjacency/paired assigner's default edit-distance-1 clustering
//! merging two that were meant to stay separate MI families.

use crate::helpers::bam_generator::{create_minimal_header, write_bam};
use fgumi_raw_bam::{RawRecord, SamBuilder, SamTag, flags};
use rstest::rstest;
use std::path::{Path, PathBuf};
use std::process::Command;

/// Runs the compiled `fgumi` binary as a subprocess, matching this repo's
/// real integration-test convention (`tests/integration/helpers/cli.rs`).
fn run_fgumi(args: &[&str]) {
    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args(args)
        .status()
        .unwrap_or_else(|e| panic!("failed to spawn fgumi {args:?}: {e}"));
    assert!(status.success(), "fgumi {args:?} failed with {status}");
}

/// Appends `.{suffix}` to a metrics prefix path (mirrors
/// `commands::group::with_extension`, which integration tests cannot import
/// as it is `pub(crate)`).
fn suffixed(prefix: &Path, suffix: &str) -> PathBuf {
    PathBuf::from(format!("{}.{suffix}", prefix.display()))
}

/// Asserts that the `.{suffix}` metrics file under `actual_prefix` is byte-
/// identical to the one under `expected_prefix`, panicking with both paths
/// named on either a missing file or a content mismatch.
fn assert_metrics_file_eq(
    actual_prefix: &Path,
    expected_prefix: &Path,
    suffix: &str,
    context: &str,
) {
    let actual = std::fs::read_to_string(suffixed(actual_prefix, suffix)).unwrap_or_else(|e| {
        panic!("{context}: missing {suffix} at {}: {e}", actual_prefix.display())
    });
    let expected = std::fs::read_to_string(suffixed(expected_prefix, suffix)).unwrap_or_else(|e| {
        panic!("{context}: missing {suffix} at {}: {e}", expected_prefix.display())
    });
    // Non-vacuity guard: a `family_sizes` file with only a header (no data
    // rows) means the input produced no counted families — which happens
    // silently if the test data is single-end, since every metrics path drops
    // non-PAIRED records. An empty-vs-empty comparison would pass without
    // exercising the metric at all, so require at least one data row on the
    // family-size files (the ones that are empty exactly when nothing was
    // counted).
    if suffix.ends_with("family_sizes.txt") {
        assert!(
            expected.lines().count() >= 2,
            "{context}: ground-truth {suffix} has no data rows — the parity assertion would be \
             vacuous (is the test input PAIRED? simplex metrics count paired templates only)",
        );
    }
    assert_eq!(actual, expected, "{context}: {suffix} diverges from ground truth");
}

// ============================================================================
// Step 1 (task brief, verbatim modulo the `with_extension_appended` fix the
// brief itself calls out): simplex, fused vs. standalone vs. separate-pass,
// single-threaded, one family per case.
// ============================================================================

/// Builds one simplex read PAIR (R1 + R2, one single-strand UMI) sharing
/// `name`, `umi`, and optionally an `mi` tag. R1 sits at `r1_pos` spanning
/// `r1_len`; R2 is stacked immediately after it (reverse strand), forming a
/// valid FR template.
///
/// Simplex metrics count PAIRED templates only — every metrics path
/// (`pair_records_by_read_name`, and the separate-pass `simplex-metrics`)
/// drops non-`PAIRED` records — so these builders MUST emit real pairs. With
/// single-end records (`flags(0)`) both the inline output and the ground
/// truth come out empty and the parity assertions pass vacuously.
fn simplex_pair(
    name: &str,
    umi: &str,
    mi: Option<&str>,
    r1_pos: i32,
    r1_len: usize,
) -> (RawRecord, RawRecord) {
    let r2_pos = r1_pos + i32::try_from(r1_len).expect("r1_len fits i32");
    let r1_cigar = u32::try_from(r1_len).expect("r1_len fits u32") << 4;
    let r2_cigar = 10u32 << 4;
    let tlen = (r2_pos + 10) - r1_pos;

    let mut b1 = SamBuilder::new();
    b1.read_name(name.as_bytes())
        .sequence(&vec![b'A'; r1_len])
        .qualities(&vec![30u8; r1_len])
        .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE)
        .ref_id(0)
        .pos(r1_pos)
        .mapq(60)
        .cigar_ops(&[r1_cigar])
        .mate_ref_id(0)
        .mate_pos(r2_pos)
        .template_length(tlen);
    b1.add_string_tag(SamTag::RX, umi.as_bytes());
    b1.add_string_tag(SamTag::MC, b"10M");
    if let Some(mi) = mi {
        b1.add_string_tag(SamTag::MI, mi.as_bytes());
    }
    let r1 = b1.build();

    let mut b2 = SamBuilder::new();
    b2.read_name(name.as_bytes())
        .sequence(b"ACGTACGTAC")
        .qualities(&[30u8; 10])
        .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE)
        .ref_id(0)
        .pos(r2_pos)
        .mapq(60)
        .cigar_ops(&[r2_cigar])
        .mate_ref_id(0)
        .mate_pos(r1_pos)
        .template_length(-tlen);
    b2.add_string_tag(SamTag::RX, umi.as_bytes());
    b2.add_string_tag(SamTag::MC, format!("{r1_len}M").as_bytes());
    if let Some(mi) = mi {
        b2.add_string_tag(SamTag::MI, mi.as_bytes());
    }
    let r2 = b2.build();

    (r1, r2)
}

/// Builds `family_size` read pairs, all sharing one MI (one UMI family) at
/// one position, so the whole group forms exactly one CS/SS family.
fn one_family_records(family_size: usize) -> Vec<RawRecord> {
    (0..family_size)
        .flat_map(|i| {
            let (r1, r2) = simplex_pair(&format!("read-{i}"), "ACGT-TGCA", Some("0"), 100, 10);
            [r1, r2]
        })
        .collect()
}

#[rstest]
#[case::singleton_family(1)]
#[case::small_family(5)]
#[case::large_family(200)]
fn simplex_inline_fused_matches_separate_pass_single_threaded(#[case] family_size: usize) {
    let dir = tempfile::tempdir().expect("tempdir");
    let header = create_minimal_header("chr1", 10_000);
    let bam_path = dir.path().join("in.bam");
    write_bam(&bam_path, &header, &one_family_records(family_size));

    // (1) Fused: runall --simplex::metrics
    let fused_prefix = dir.path().join("fused");
    run_fgumi(&[
        "runall",
        "-i",
        bam_path.to_str().unwrap(),
        "-o",
        dir.path().join("fused_out.bam").to_str().unwrap(),
        "--start-from",
        "group",
        "--stop-after",
        "consensus",
        "--consensus",
        "simplex",
        "--group::strategy",
        "adjacency",
        "--simplex::min-reads",
        "1",
        "--simplex::metrics",
        fused_prefix.to_str().unwrap(),
        "--threads",
        "1",
    ]);

    // (2) Ground truth: fgumi group, then separate-pass simplex-metrics
    let grouped_bam = dir.path().join("grouped.bam");
    run_fgumi(&[
        "group",
        "-i",
        bam_path.to_str().unwrap(),
        "-o",
        grouped_bam.to_str().unwrap(),
        "-s",
        "adjacency",
    ]);

    let ground_truth_prefix = dir.path().join("ground_truth");
    run_fgumi(&[
        "simplex-metrics",
        "-i",
        grouped_bam.to_str().unwrap(),
        "-o",
        ground_truth_prefix.to_str().unwrap(),
        "--min-reads",
        "1",
    ]);

    // (3) Standalone consensus command's own --metrics, fed the grouped BAM (T2)
    let standalone_prefix = dir.path().join("standalone");
    run_fgumi(&[
        "simplex",
        "-i",
        grouped_bam.to_str().unwrap(),
        "-o",
        dir.path().join("standalone_out.bam").to_str().unwrap(),
        "--min-reads",
        "1",
        "--metrics",
        standalone_prefix.to_str().unwrap(),
    ]);

    for suffix in ["family_sizes.txt", "umi_counts.txt", "simplex_yield_metrics.txt"] {
        assert_metrics_file_eq(
            &fused_prefix,
            &ground_truth_prefix,
            suffix,
            "fused vs ground truth",
        );
        assert_metrics_file_eq(
            &standalone_prefix,
            &ground_truth_prefix,
            suffix,
            "standalone vs ground truth",
        );
    }
}

// ============================================================================
// Extended matrix (Step 3): {simplex, duplex, codec} x family shapes x
// thread counts, plus non-default min-reads, intervals, and the dedicated
// multi-batch/multi-worker T2 hazard case.
// ============================================================================

/// Maps an index to an 8-base, exact-match-distinguishable pseudo-UMI.
///
/// Base-4 encoding of `i` over `{A,C,G,T}` — unique for any `i < 4^8 =
/// 65536`, comfortably above every family count used below. Uniqueness (not
/// edit-distance separation) is all that is required: every `group`
/// invocation in this file passes `--edits 0` / `-e 0`, so UMI assignment is
/// exact-match only and clustering cannot merge two distinct indices anyway.
fn indexed_umi(i: usize) -> String {
    const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
    let mut n = i;
    let mut out = [0u8; 8];
    for slot in out.iter_mut().rev() {
        *slot = BASES[n % 4];
        n /= 4;
    }
    String::from_utf8(out.to_vec()).expect("ACGT bytes are valid utf8")
}

/// A family-size/MI-count shape exercised across every mode and thread
/// count in the matrix below.
#[derive(Clone, Copy)]
enum Shape {
    /// One MI, one read/pair.
    Singleton,
    /// One MI, a handful of reads/pairs.
    Small,
    /// One MI, many reads/pairs.
    Large,
    /// One coordinate group (position group) split across many distinct
    /// MIs, well under `GroupByMi::DEFAULT_TARGET_BATCH_COUNT` (50).
    ManyMis,
    /// One coordinate group whose MI count straddles the 50-MI batch
    /// boundary `GroupByMi` packages completed groups into (exercising the
    /// standalone T2 path's `split_into_runs`/`classify_batch_runs`
    /// batch-boundary handling).
    StraddleBatch,
}

impl Shape {
    fn n_families(self) -> usize {
        match self {
            Shape::Singleton | Shape::Small | Shape::Large => 1,
            Shape::ManyMis => 15,
            Shape::StraddleBatch => 55,
        }
    }

    fn depth(self) -> usize {
        match self {
            Shape::Singleton | Shape::StraddleBatch => 1,
            Shape::Small => 5,
            Shape::Large => 200,
            Shape::ManyMis => 2,
        }
    }
}

/// Builds single-end, ungrouped (no `MI` tag) simplex-mode input for `shape`:
/// `n_families` distinct single-index UMIs at one position, `depth` reads
/// each — so `group -s identity -e 0`/`-s adjacency -e 0` assigns exactly
/// `n_families` distinct MI values, all in one coordinate group.
fn simplex_shape_records(shape: Shape) -> Vec<RawRecord> {
    let mut records = Vec::new();
    for f in 0..shape.n_families() {
        let umi = indexed_umi(f);
        for i in 0..shape.depth() {
            let (r1, r2) = simplex_pair(&format!("f{f}_r{i}"), &umi, None, 100, 10);
            records.push(r1);
            records.push(r2);
        }
    }
    records
}

/// One FR read pair at `(r1_pos, r2_pos)` (0-based) sharing `umi` (dual,
/// `-`-delimited), with R1/R2 lengths independently controllable so two
/// families that share a coordinate-group key can still have distinct
/// genomic spans (used by the interval-boundary cases below). No `MI` tag:
/// `group -s paired` assigns it.
#[allow(clippy::too_many_arguments)]
fn duplex_pair(
    name: &str,
    umi: &str,
    r1_pos: i32,
    r1_len: usize,
    r2_pos: i32,
    r2_len: usize,
) -> (RawRecord, RawRecord) {
    let r1_seq = vec![b'A'; r1_len];
    let r2_seq = vec![b'A'; r2_len];
    let r1_cigar = u32::try_from(r1_len).expect("r1_len fits u32") << 4;
    let r2_cigar = u32::try_from(r2_len).expect("r2_len fits u32") << 4;
    let r1_end = r1_pos + i32::try_from(r1_len).expect("r1_len fits i32");
    let r2_end = r2_pos + i32::try_from(r2_len).expect("r2_len fits i32");
    let tlen = r2_end.max(r1_end) - r1_pos.min(r2_pos);

    let mut b1 = SamBuilder::new();
    b1.read_name(name.as_bytes())
        .sequence(&r1_seq)
        .qualities(&vec![30u8; r1_len])
        .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE)
        .ref_id(0)
        .pos(r1_pos)
        .mapq(60)
        .cigar_ops(&[r1_cigar])
        .mate_ref_id(0)
        .mate_pos(r2_pos)
        .template_length(tlen);
    b1.add_string_tag(SamTag::RX, umi.as_bytes());
    b1.add_string_tag(SamTag::MC, format!("{r2_len}M").as_bytes());
    let r1 = b1.build();

    let mut b2 = SamBuilder::new();
    b2.read_name(name.as_bytes())
        .sequence(&r2_seq)
        .qualities(&vec![30u8; r2_len])
        .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE)
        .ref_id(0)
        .pos(r2_pos)
        .mapq(60)
        .cigar_ops(&[r2_cigar])
        .mate_ref_id(0)
        .mate_pos(r1_pos)
        .template_length(-tlen);
    b2.add_string_tag(SamTag::RX, umi.as_bytes());
    b2.add_string_tag(SamTag::MC, format!("{r1_len}M").as_bytes());
    let r2 = b2.build();

    (r1, r2)
}

/// Builds ungrouped duplex/codec-mode input for `shape`: `n_families`
/// distinct dual-UMI, single-strand molecules at one position, `depth` read
/// pairs each. `group -s paired -e 0` assigns each its own `/A`-suffixed MI
/// (single strand is sufficient — a molecule need not have a reciprocal `/B`
/// partner to form its own `MiGroup`; see `test_paired_mode_with_absent_umi_*`
/// in `group.rs`), so this yields exactly `n_families` distinct MI groups
/// in one coordinate group, matching [`simplex_shape_records`]'s shape.
fn duplex_shape_records(shape: Shape) -> Vec<RawRecord> {
    let mut records = Vec::new();
    for f in 0..shape.n_families() {
        let umi = format!("{}-{}", indexed_umi(f), indexed_umi(f + 10_000));
        for i in 0..shape.depth() {
            let (r1, r2) = duplex_pair(&format!("f{f}_r{i}"), &umi, 100, 10, 200, 10);
            records.push(r1);
            records.push(r2);
        }
    }
    records
}

/// Builds one full duplex molecule (both `/A` and `/B` strands present, with
/// independently-controllable per-strand depth) for the non-default
/// `min-reads` sub-case. `-s paired` unifies the two strands into one base
/// MI: strand B swaps both the UMI halves and the R1/R2 positions relative
/// to strand A, mirroring the `test_paired_mode_with_absent_umi_on_*`
/// convention in `group.rs`.
fn full_duplex_molecule(
    umi_a: &str,
    umi_b: &str,
    ab_depth: usize,
    ba_depth: usize,
) -> Vec<RawRecord> {
    let mut records = Vec::new();
    for i in 0..ab_depth {
        let (r1, r2) = duplex_pair(&format!("ab{i}"), umi_a, 100, 10, 200, 10);
        records.push(r1);
        records.push(r2);
    }
    for i in 0..ba_depth {
        let (r1, r2) = duplex_pair(&format!("ba{i}"), umi_b, 200, 10, 100, 10);
        records.push(r1);
        records.push(r2);
    }
    records
}

/// Two SS families sharing one `ReadInfoKey` coordinate-group (identical R1
/// start/strand and R2 start/strand) but different genomic spans (R1 length
/// varies), used by the interval-boundary cases: a `near` family whose span
/// stays close to the shared start, and a `far` family whose R1 extends well
/// downstream. An interval positioned past `near`'s end but inside `far`'s
/// span excludes exactly one member of an otherwise-single coordinate group.
fn split_interval_duplex_family() -> Vec<RawRecord> {
    let mut records = Vec::new();
    let (r1, r2) = duplex_pair("near", "AAAAAAAA-CCCCCCCC", 100, 10, 200, 10);
    records.push(r1);
    records.push(r2);
    let (r1, r2) = duplex_pair("far", "GGGGGGGG-TTTTTTTT", 100, 600, 200, 10);
    records.push(r1);
    records.push(r2);
    records
}

/// Single-end analog of [`split_interval_duplex_family`] for simplex mode:
/// two single-end reads sharing one position/strand key but differing read
/// length, so their spans differ while remaining one coordinate group.
fn split_interval_simplex_family() -> Vec<RawRecord> {
    let mut records = Vec::new();
    // "near" (R1 spans 100..110) falls outside the 500..600 interval; "far"
    // (R1 spans 100..700) overlaps it — so interval filtering must drop one
    // template and keep the other. Both must be PAIRED or they are dropped
    // before the interval filter ever runs.
    for (name, umi, len) in [("near", "AAAAAAAAAA", 10usize), ("far", "CCCCCCCCCC", 600usize)] {
        let (r1, r2) = simplex_pair(name, umi, None, 100, len);
        records.push(r1);
        records.push(r2);
    }
    records
}

/// Writes a BED-format interval file (0-based, half-open).
fn write_bed_interval(dir: &Path, chrom: &str, start: i32, end: i32) -> PathBuf {
    let path = dir.join("intervals.bed");
    std::fs::write(&path, format!("{chrom}\t{start}\t{end}\n")).expect("write bed interval");
    path
}

/// Writes a Picard-interval-list-format interval file (1-based, closed;
/// `parse_intervals` only requires a line starting with `@` to switch modes,
/// so the header need not be a fully valid SAM header).
fn write_picard_interval(dir: &Path, chrom: &str, start: i32, end: i32) -> PathBuf {
    let path = dir.join("intervals.interval_list");
    std::fs::write(
        &path,
        format!("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:{chrom}\tLN:10000\n{chrom}\t{start}\t{end}\t+\tregion\n"),
    )
    .expect("write picard interval list");
    path
}

// ---------------------------------------------------------------------------
// Simplex matrix
// ---------------------------------------------------------------------------

fn simplex_min_reads_arg(min_reads: usize) -> String {
    min_reads.to_string()
}

fn run_simplex_fused(
    dir: &Path,
    input: &Path,
    tag: &str,
    threads: usize,
    min_reads: usize,
    intervals: Option<&Path>,
) -> PathBuf {
    let prefix = dir.join(format!("fused_{tag}"));
    let out = dir.join(format!("fused_{tag}_out.bam"));
    let min_reads_arg = simplex_min_reads_arg(min_reads);
    let threads_arg = threads.to_string();
    let mut args = vec![
        "runall",
        "-i",
        input.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--start-from",
        "group",
        "--stop-after",
        "consensus",
        "--consensus",
        "simplex",
        "--group::strategy",
        "adjacency",
        "--group::edits",
        "0",
        "--simplex::min-reads",
        &min_reads_arg,
        "--simplex::metrics",
        prefix.to_str().unwrap(),
        "--threads",
        &threads_arg,
    ];
    if let Some(iv) = intervals {
        args.push("--simplex::intervals");
        args.push(iv.to_str().unwrap());
    }
    run_fgumi(&args);
    prefix
}

/// Builds the grouped-BAM + ground-truth simplex-metrics pair for `input`.
/// Returns `(grouped_bam, ground_truth_prefix)`.
fn simplex_ground_truth(
    dir: &Path,
    input: &Path,
    min_reads: usize,
    intervals: Option<&Path>,
) -> (PathBuf, PathBuf) {
    let grouped = dir.join("grouped.bam");
    run_fgumi(&[
        "group",
        "-i",
        input.to_str().unwrap(),
        "-o",
        grouped.to_str().unwrap(),
        "-s",
        "adjacency",
        "-e",
        "0",
    ]);

    let prefix = dir.join("ground_truth");
    let min_reads_arg = simplex_min_reads_arg(min_reads);
    let mut args = vec![
        "simplex-metrics",
        "-i",
        grouped.to_str().unwrap(),
        "-o",
        prefix.to_str().unwrap(),
        "--min-reads",
        &min_reads_arg,
    ];
    if let Some(iv) = intervals {
        args.push("--intervals");
        args.push(iv.to_str().unwrap());
    }
    run_fgumi(&args);
    (grouped, prefix)
}

fn run_simplex_standalone(
    dir: &Path,
    grouped: &Path,
    tag: &str,
    threads: Option<usize>,
    min_reads: usize,
    intervals: Option<&Path>,
) -> PathBuf {
    let prefix = dir.join(format!("standalone_{tag}"));
    let out = dir.join(format!("standalone_{tag}_out.bam"));
    let min_reads_arg = simplex_min_reads_arg(min_reads);
    let threads_arg = threads.map(|t| t.to_string());
    let mut args = vec![
        "simplex",
        "-i",
        grouped.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--min-reads",
        &min_reads_arg,
        "--metrics",
        prefix.to_str().unwrap(),
    ];
    if let Some(t) = &threads_arg {
        args.push("--threads");
        args.push(t);
    }
    if let Some(iv) = intervals {
        args.push("--intervals");
        args.push(iv.to_str().unwrap());
    }
    run_fgumi(&args);
    prefix
}

const SIMPLEX_SUFFIXES: [&str; 3] =
    ["family_sizes.txt", "umi_counts.txt", "simplex_yield_metrics.txt"];
const DUPLEX_SUFFIXES: [&str; 4] =
    ["family_sizes.txt", "duplex_family_sizes.txt", "umi_counts.txt", "duplex_yield_metrics.txt"];

#[rstest]
#[case::singleton(Shape::Singleton)]
#[case::small(Shape::Small)]
#[case::large(Shape::Large)]
#[case::many_mis(Shape::ManyMis)]
#[case::straddle_batch(Shape::StraddleBatch)]
fn simplex_matrix(#[case] shape: Shape, #[values(1, 2, 8)] threads: usize) {
    let dir = tempfile::tempdir().expect("tempdir");
    let header = create_minimal_header("chr1", 10_000);
    let bam_path = dir.path().join("in.bam");
    write_bam(&bam_path, &header, &simplex_shape_records(shape));

    let fused_prefix = run_simplex_fused(dir.path(), &bam_path, "matrix", threads, 1, None);
    let (grouped, ground_truth_prefix) = simplex_ground_truth(dir.path(), &bam_path, 1, None);
    let standalone_prefix =
        run_simplex_standalone(dir.path(), &grouped, "matrix", Some(threads), 1, None);

    for suffix in SIMPLEX_SUFFIXES {
        assert_metrics_file_eq(
            &fused_prefix,
            &ground_truth_prefix,
            suffix,
            "simplex fused vs truth",
        );
        assert_metrics_file_eq(
            &standalone_prefix,
            &ground_truth_prefix,
            suffix,
            "simplex standalone vs truth",
        );
    }
}

#[rstest]
#[case::bed(false)]
#[case::picard(true)]
fn simplex_intervals(#[case] picard: bool) {
    let dir = tempfile::tempdir().expect("tempdir");
    let header = create_minimal_header("chr1", 10_000);
    let bam_path = dir.path().join("in.bam");
    write_bam(&bam_path, &header, &split_interval_simplex_family());

    let intervals = if picard {
        write_picard_interval(dir.path(), "chr1", 501, 600)
    } else {
        write_bed_interval(dir.path(), "chr1", 500, 600)
    };

    let fused_prefix = run_simplex_fused(dir.path(), &bam_path, "iv", 1, 1, Some(&intervals));
    let (grouped, ground_truth_prefix) =
        simplex_ground_truth(dir.path(), &bam_path, 1, Some(&intervals));
    let standalone_prefix =
        run_simplex_standalone(dir.path(), &grouped, "iv", None, 1, Some(&intervals));

    for suffix in SIMPLEX_SUFFIXES {
        assert_metrics_file_eq(
            &fused_prefix,
            &ground_truth_prefix,
            suffix,
            "simplex fused vs truth (intervals)",
        );
        assert_metrics_file_eq(
            &standalone_prefix,
            &ground_truth_prefix,
            suffix,
            "simplex standalone vs truth (intervals)",
        );
    }
}

// ---------------------------------------------------------------------------
// Duplex matrix
// ---------------------------------------------------------------------------

fn duplex_min_reads_arg(min_reads: &[usize]) -> String {
    min_reads.iter().map(usize::to_string).collect::<Vec<_>>().join(",")
}

/// fgbio's `padTo(3, last)` (also `DuplexConsensusCaller::min_yx_reads_for`):
/// `min_ab = min_reads.get(1).unwrap_or(last)`, `min_ba =
/// min_reads.get(2).or(last)`. Mirrors `builder.rs`'s inline-metrics wiring
/// exactly so the separate-pass ground truth uses the same thresholds the
/// fused/standalone paths derive internally.
fn duplex_ab_ba(min_reads: &[usize]) -> (usize, usize) {
    let last = *min_reads.last().expect("non-empty min_reads");
    let ab = min_reads.get(1).copied().unwrap_or(last);
    let ba = min_reads.get(2).copied().unwrap_or(last);
    (ab, ba)
}

fn run_duplex_fused(
    dir: &Path,
    input: &Path,
    tag: &str,
    threads: usize,
    min_reads: &[usize],
    intervals: Option<&Path>,
) -> PathBuf {
    let prefix = dir.join(format!("fused_{tag}"));
    let out = dir.join(format!("fused_{tag}_out.bam"));
    let min_reads_arg = duplex_min_reads_arg(min_reads);
    let threads_arg = threads.to_string();
    let mut args = vec![
        "runall",
        "-i",
        input.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--start-from",
        "group",
        "--stop-after",
        "consensus",
        "--consensus",
        "duplex",
        "--group::strategy",
        "paired",
        "--group::edits",
        "0",
        "--duplex::min-reads",
        &min_reads_arg,
        "--duplex::metrics",
        prefix.to_str().unwrap(),
        "--threads",
        &threads_arg,
    ];
    if let Some(iv) = intervals {
        args.push("--duplex::intervals");
        args.push(iv.to_str().unwrap());
    }
    run_fgumi(&args);
    prefix
}

/// Builds the grouped-BAM + ground-truth duplex-metrics pair for `input`,
/// via the `min-ab-reads`/`min-ba-reads` fgbio `padTo(3, last)` projection of
/// `min_reads`. Returns `(grouped_bam, ground_truth_prefix)`.
fn duplex_ground_truth(
    dir: &Path,
    input: &Path,
    min_reads: &[usize],
    intervals: Option<&Path>,
) -> (PathBuf, PathBuf) {
    let grouped = dir.join("grouped.bam");
    run_fgumi(&[
        "group",
        "-i",
        input.to_str().unwrap(),
        "-o",
        grouped.to_str().unwrap(),
        "-s",
        "paired",
        "-e",
        "0",
    ]);

    let (ab, ba) = duplex_ab_ba(min_reads);
    let ab_arg = ab.to_string();
    let ba_arg = ba.to_string();
    let prefix = dir.join("ground_truth");
    let mut args = vec![
        "duplex-metrics",
        "-i",
        grouped.to_str().unwrap(),
        "-o",
        prefix.to_str().unwrap(),
        "--min-ab-reads",
        &ab_arg,
        "--min-ba-reads",
        &ba_arg,
    ];
    if let Some(iv) = intervals {
        args.push("--intervals");
        args.push(iv.to_str().unwrap());
    }
    run_fgumi(&args);
    (grouped, prefix)
}

fn run_duplex_standalone(
    dir: &Path,
    grouped: &Path,
    tag: &str,
    threads: Option<usize>,
    min_reads: &[usize],
    intervals: Option<&Path>,
) -> PathBuf {
    let prefix = dir.join(format!("standalone_{tag}"));
    let out = dir.join(format!("standalone_{tag}_out.bam"));
    let min_reads_arg = duplex_min_reads_arg(min_reads);
    let threads_arg = threads.map(|t| t.to_string());
    let mut args = vec![
        "duplex",
        "-i",
        grouped.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--min-reads",
        &min_reads_arg,
        "--metrics",
        prefix.to_str().unwrap(),
    ];
    if let Some(t) = &threads_arg {
        args.push("--threads");
        args.push(t);
    }
    if let Some(iv) = intervals {
        args.push("--intervals");
        args.push(iv.to_str().unwrap());
    }
    run_fgumi(&args);
    prefix
}

#[rstest]
#[case::singleton(Shape::Singleton)]
#[case::small(Shape::Small)]
#[case::large(Shape::Large)]
#[case::many_mis(Shape::ManyMis)]
#[case::straddle_batch(Shape::StraddleBatch)]
fn duplex_matrix(#[case] shape: Shape, #[values(1, 2, 8)] threads: usize) {
    let dir = tempfile::tempdir().expect("tempdir");
    let header = create_minimal_header("chr1", 10_000);
    let bam_path = dir.path().join("in.bam");
    write_bam(&bam_path, &header, &duplex_shape_records(shape));

    let min_reads = [1usize];
    let fused_prefix = run_duplex_fused(dir.path(), &bam_path, "matrix", threads, &min_reads, None);
    let (grouped, ground_truth_prefix) =
        duplex_ground_truth(dir.path(), &bam_path, &min_reads, None);
    let standalone_prefix =
        run_duplex_standalone(dir.path(), &grouped, "matrix", Some(threads), &min_reads, None);

    for suffix in DUPLEX_SUFFIXES {
        assert_metrics_file_eq(
            &fused_prefix,
            &ground_truth_prefix,
            suffix,
            "duplex fused vs truth",
        );
        assert_metrics_file_eq(
            &standalone_prefix,
            &ground_truth_prefix,
            suffix,
            "duplex standalone vs truth",
        );
    }
}

/// Non-default `--duplex::min-reads` sub-case. `validate_min_reads` requires
/// each value no larger than the one before it (`total >= XY >= YX`), so
/// `3,3,1` (not the task brief's `2,3,1`, which violates that ordering and is
/// rejected by the CLI — a documentation slip in the brief, fixed here) pads
/// to `(total=3, xy=3, yx=1)`, projecting to `--min-ab-reads=3
/// --min-ba-reads=1` on the separate-pass ground truth, matching the value
/// pair the brief names. Uses a real two-strand molecule (`/A` depth 4, `/B`
/// depth 1) so the AB/BA split is meaningful rather than trivially symmetric.
#[test]
fn duplex_non_default_min_reads() {
    let dir = tempfile::tempdir().expect("tempdir");
    let header = create_minimal_header("chr1", 10_000);
    let bam_path = dir.path().join("in.bam");
    let records = full_duplex_molecule("AAAAAAAA-CCCCCCCC", "CCCCCCCC-AAAAAAAA", 4, 1);
    write_bam(&bam_path, &header, &records);

    let min_reads = [3usize, 3, 1];
    let fused_prefix = run_duplex_fused(dir.path(), &bam_path, "mr", 2, &min_reads, None);
    let (grouped, ground_truth_prefix) =
        duplex_ground_truth(dir.path(), &bam_path, &min_reads, None);
    let standalone_prefix =
        run_duplex_standalone(dir.path(), &grouped, "mr", Some(2), &min_reads, None);

    for suffix in DUPLEX_SUFFIXES {
        assert_metrics_file_eq(
            &fused_prefix,
            &ground_truth_prefix,
            suffix,
            "duplex fused vs truth (non-default min-reads)",
        );
        assert_metrics_file_eq(
            &standalone_prefix,
            &ground_truth_prefix,
            suffix,
            "duplex standalone vs truth (non-default min-reads)",
        );
    }
}

#[rstest]
#[case::bed(false)]
#[case::picard(true)]
fn duplex_intervals(#[case] picard: bool) {
    let dir = tempfile::tempdir().expect("tempdir");
    let header = create_minimal_header("chr1", 10_000);
    let bam_path = dir.path().join("in.bam");
    write_bam(&bam_path, &header, &split_interval_duplex_family());

    let intervals = if picard {
        write_picard_interval(dir.path(), "chr1", 501, 600)
    } else {
        write_bed_interval(dir.path(), "chr1", 500, 600)
    };

    let min_reads = [1usize];
    let fused_prefix =
        run_duplex_fused(dir.path(), &bam_path, "iv", 1, &min_reads, Some(&intervals));
    let (grouped, ground_truth_prefix) =
        duplex_ground_truth(dir.path(), &bam_path, &min_reads, Some(&intervals));
    let standalone_prefix =
        run_duplex_standalone(dir.path(), &grouped, "iv", None, &min_reads, Some(&intervals));

    for suffix in DUPLEX_SUFFIXES {
        assert_metrics_file_eq(
            &fused_prefix,
            &ground_truth_prefix,
            suffix,
            "duplex fused vs truth (intervals)",
        );
        assert_metrics_file_eq(
            &standalone_prefix,
            &ground_truth_prefix,
            suffix,
            "duplex standalone vs truth (intervals)",
        );
    }
}

// ---------------------------------------------------------------------------
// Codec matrix — codec's inline finalize hook uses `MetricsThresholds::Duplex`
// with `min_ab_reads == min_ba_reads == codec.min_reads` (a single symmetric
// scalar threshold; builder.rs's `add_codec`), so ground truth is
// `duplex-metrics --min-ab-reads=N --min-ba-reads=N` on the grouped BAM, and
// the output file set is duplex-metrics' four files, matching
// `runall_group_to_codec_metrics_succeeds_on_paired_end_input` in
// `test_group_metrics_fused_parity.rs`.
// ---------------------------------------------------------------------------

fn run_codec_fused(
    dir: &Path,
    input: &Path,
    tag: &str,
    threads: usize,
    min_reads: usize,
    intervals: Option<&Path>,
) -> PathBuf {
    let prefix = dir.join(format!("fused_{tag}"));
    let out = dir.join(format!("fused_{tag}_out.bam"));
    let min_reads_arg = min_reads.to_string();
    let threads_arg = threads.to_string();
    let mut args = vec![
        "runall",
        "-i",
        input.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--start-from",
        "group",
        "--stop-after",
        "consensus",
        "--consensus",
        "codec",
        "--group::strategy",
        "adjacency",
        "--group::edits",
        "0",
        "--codec::min-reads",
        &min_reads_arg,
        "--codec::metrics",
        prefix.to_str().unwrap(),
        "--threads",
        &threads_arg,
    ];
    if let Some(iv) = intervals {
        args.push("--codec::intervals");
        args.push(iv.to_str().unwrap());
    }
    run_fgumi(&args);
    prefix
}

fn codec_ground_truth(
    dir: &Path,
    input: &Path,
    min_reads: usize,
    intervals: Option<&Path>,
) -> (PathBuf, PathBuf) {
    let grouped = dir.join("grouped.bam");
    run_fgumi(&[
        "group",
        "-i",
        input.to_str().unwrap(),
        "-o",
        grouped.to_str().unwrap(),
        "-s",
        "adjacency",
        "-e",
        "0",
    ]);

    let min_reads_arg = min_reads.to_string();
    let prefix = dir.join("ground_truth");
    let mut args = vec![
        "duplex-metrics",
        "-i",
        grouped.to_str().unwrap(),
        "-o",
        prefix.to_str().unwrap(),
        "--min-ab-reads",
        &min_reads_arg,
        "--min-ba-reads",
        &min_reads_arg,
    ];
    if let Some(iv) = intervals {
        args.push("--intervals");
        args.push(iv.to_str().unwrap());
    }
    run_fgumi(&args);
    (grouped, prefix)
}

fn run_codec_standalone(
    dir: &Path,
    grouped: &Path,
    tag: &str,
    threads: Option<usize>,
    min_reads: usize,
    intervals: Option<&Path>,
) -> PathBuf {
    let prefix = dir.join(format!("standalone_{tag}"));
    let out = dir.join(format!("standalone_{tag}_out.bam"));
    let min_reads_arg = min_reads.to_string();
    let threads_arg = threads.map(|t| t.to_string());
    let mut args = vec![
        "codec",
        "-i",
        grouped.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--min-reads",
        &min_reads_arg,
        "--metrics",
        prefix.to_str().unwrap(),
    ];
    if let Some(t) = &threads_arg {
        args.push("--threads");
        args.push(t);
    }
    if let Some(iv) = intervals {
        args.push("--intervals");
        args.push(iv.to_str().unwrap());
    }
    run_fgumi(&args);
    prefix
}

#[rstest]
#[case::singleton(Shape::Singleton)]
#[case::small(Shape::Small)]
#[case::large(Shape::Large)]
#[case::many_mis(Shape::ManyMis)]
#[case::straddle_batch(Shape::StraddleBatch)]
fn codec_matrix(#[case] shape: Shape, #[values(1, 2, 8)] threads: usize) {
    let dir = tempfile::tempdir().expect("tempdir");
    let header = create_minimal_header("chr1", 10_000);
    let bam_path = dir.path().join("in.bam");
    write_bam(&bam_path, &header, &duplex_shape_records(shape));

    let fused_prefix = run_codec_fused(dir.path(), &bam_path, "matrix", threads, 1, None);
    let (grouped, ground_truth_prefix) = codec_ground_truth(dir.path(), &bam_path, 1, None);
    let standalone_prefix =
        run_codec_standalone(dir.path(), &grouped, "matrix", Some(threads), 1, None);

    for suffix in DUPLEX_SUFFIXES {
        assert_metrics_file_eq(&fused_prefix, &ground_truth_prefix, suffix, "codec fused vs truth");
        assert_metrics_file_eq(
            &standalone_prefix,
            &ground_truth_prefix,
            suffix,
            "codec standalone vs truth",
        );
    }
}

/// Non-default `--codec::min-reads` sub-case (codec's threshold is a single
/// symmetric scalar, unlike duplex's up-to-3-value vector).
#[test]
fn codec_non_default_min_reads() {
    let dir = tempfile::tempdir().expect("tempdir");
    let header = create_minimal_header("chr1", 10_000);
    let bam_path = dir.path().join("in.bam");
    let records = full_duplex_molecule("AAAAAAAA-CCCCCCCC", "CCCCCCCC-AAAAAAAA", 4, 3);
    write_bam(&bam_path, &header, &records);

    let fused_prefix = run_codec_fused(dir.path(), &bam_path, "mr", 2, 3, None);
    let (grouped, ground_truth_prefix) = codec_ground_truth(dir.path(), &bam_path, 3, None);
    let standalone_prefix = run_codec_standalone(dir.path(), &grouped, "mr", Some(2), 3, None);

    for suffix in DUPLEX_SUFFIXES {
        assert_metrics_file_eq(
            &fused_prefix,
            &ground_truth_prefix,
            suffix,
            "codec fused vs truth (non-default min-reads)",
        );
        assert_metrics_file_eq(
            &standalone_prefix,
            &ground_truth_prefix,
            suffix,
            "codec standalone vs truth (non-default min-reads)",
        );
    }
}

#[rstest]
#[case::bed(false)]
#[case::picard(true)]
fn codec_intervals(#[case] picard: bool) {
    let dir = tempfile::tempdir().expect("tempdir");
    let header = create_minimal_header("chr1", 10_000);
    let bam_path = dir.path().join("in.bam");
    write_bam(&bam_path, &header, &split_interval_duplex_family());

    let intervals = if picard {
        write_picard_interval(dir.path(), "chr1", 501, 600)
    } else {
        write_bed_interval(dir.path(), "chr1", 500, 600)
    };

    let fused_prefix = run_codec_fused(dir.path(), &bam_path, "iv", 1, 1, Some(&intervals));
    let (grouped, ground_truth_prefix) =
        codec_ground_truth(dir.path(), &bam_path, 1, Some(&intervals));
    let standalone_prefix =
        run_codec_standalone(dir.path(), &grouped, "iv", None, 1, Some(&intervals));

    for suffix in DUPLEX_SUFFIXES {
        assert_metrics_file_eq(
            &fused_prefix,
            &ground_truth_prefix,
            suffix,
            "codec fused vs truth (intervals)",
        );
        assert_metrics_file_eq(
            &standalone_prefix,
            &ground_truth_prefix,
            suffix,
            "codec standalone vs truth (intervals)",
        );
    }
}

// ============================================================================
// The 3-batch/3-worker case: this task's own correctness anchor for the
// standalone T2 path's cross-batch boundary-detection mechanism specifically
// (`split_into_runs`/`classify_batch_runs`/`reassemble_boundary`). One
// coordinate group with 130 distinct MIs — spanning 3 of `GroupByMi`'s
// 50-MI-group batches (`DEFAULT_TARGET_BATCH_COUNT = 50`: 50 + 50 + 30) —
// run at `--threads 8`. Only the *standalone* T2 path is exposed to this
// hazard (T1 trusts `GroupByPosition`'s pre-formed groups directly and never
// re-derives boundaries from a batched MI stream), so this compares T2
// against ground truth only, per the brief.
// ============================================================================

#[test]
fn three_batch_multi_worker_t2_matches_ground_truth() {
    let dir = tempfile::tempdir().expect("tempdir");
    let header = create_minimal_header("chr1", 10_000);
    let bam_path = dir.path().join("in.bam");

    let mut records = Vec::new();
    for f in 0..130 {
        let umi = indexed_umi(f);
        let (r1, r2) = simplex_pair(&format!("f{f}"), &umi, None, 100, 10);
        records.push(r1);
        records.push(r2);
    }
    write_bam(&bam_path, &header, &records);

    let (grouped, ground_truth_prefix) = simplex_ground_truth(dir.path(), &bam_path, 1, None);
    let standalone_prefix =
        run_simplex_standalone(dir.path(), &grouped, "3batch", Some(8), 1, None);

    for suffix in SIMPLEX_SUFFIXES {
        assert_metrics_file_eq(
            &standalone_prefix,
            &ground_truth_prefix,
            suffix,
            "T2 (3-batch, 8 workers) vs ground truth",
        );
    }
}

// ============================================================================
// Multi-key batch-straddle case (parallel T2 producer, Task 3's own
// correctness anchor): 3 DISTINCT coordinate keys, 40 MI families each (120
// MI groups total). `GroupByMi`'s 50-MI-group batches (`target_batch_count`)
// then split this stream into 3 batches (50 + 50 + 20) whose boundaries do
// NOT align with the coordinate-key boundaries (every 40 groups) — so within
// a single batch, `split_into_runs`/`classify_batch_runs` sees a Head run
// (the tail of one key), one or more interior (complete) runs, and a Tail run
// (the start of the next key), exercising all three run kinds within the same
// batch as well as across batches. Only the standalone T2 path is exposed to
// this hazard (see the 3-batch case above), so this compares T2 against
// ground truth only.
// ============================================================================

const MULTI_KEY_FAMILIES_PER_KEY: usize = 40;

#[test]
fn multi_key_batch_straddle_t2_matches_ground_truth() {
    let dir = tempfile::tempdir().expect("tempdir");
    let header = create_minimal_header("chr1", 10_000);
    let bam_path = dir.path().join("in.bam");

    let mut records = Vec::new();
    for k in 0..3 {
        // Well-separated positions so each `k` forms its own coordinate group.
        let r1_pos = 100 + i32::try_from(k).expect("k fits i32") * 2000;
        for f in 0..MULTI_KEY_FAMILIES_PER_KEY {
            let umi = indexed_umi(k * MULTI_KEY_FAMILIES_PER_KEY + f);
            let (r1, r2) = simplex_pair(&format!("k{k}f{f}"), &umi, None, r1_pos, 10);
            records.push(r1);
            records.push(r2);
        }
    }
    write_bam(&bam_path, &header, &records);

    let (grouped, ground_truth_prefix) = simplex_ground_truth(dir.path(), &bam_path, 1, None);
    let standalone_prefix =
        run_simplex_standalone(dir.path(), &grouped, "multikey", Some(8), 1, None);

    for suffix in SIMPLEX_SUFFIXES {
        assert_metrics_file_eq(
            &standalone_prefix,
            &ground_truth_prefix,
            suffix,
            "T2 (multi-key batch straddle, 8 workers) vs ground truth",
        );
    }
}

// ============================================================================
// Duplex analogs of the two parallel-T2 correctness anchors above (Task 4):
// many single-strand `duplex_pair` families at one coordinate key spanning
// three `GroupByMi` batches, and three distinct coordinate keys whose MI
// counts straddle those batch boundaries. Same rationale as the simplex
// cases: only the standalone T2 path re-derives coordinate-group boundaries
// from a batched MI stream, so these compare T2 against the separate-pass
// `duplex-metrics` ground truth only.
// ============================================================================

#[test]
fn duplex_three_batch_multi_worker_t2_matches_ground_truth() {
    let dir = tempfile::tempdir().expect("tempdir");
    let header = create_minimal_header("chr1", 10_000);
    let bam_path = dir.path().join("in.bam");

    let mut records = Vec::new();
    for f in 0..130 {
        let umi = format!("{}-{}", indexed_umi(f), indexed_umi(f + 10_000));
        let (r1, r2) = duplex_pair(&format!("f{f}"), &umi, 100, 10, 200, 10);
        records.push(r1);
        records.push(r2);
    }
    write_bam(&bam_path, &header, &records);

    let min_reads = [1usize];
    let (grouped, ground_truth_prefix) =
        duplex_ground_truth(dir.path(), &bam_path, &min_reads, None);
    let standalone_prefix =
        run_duplex_standalone(dir.path(), &grouped, "3batch", Some(8), &min_reads, None);

    for suffix in DUPLEX_SUFFIXES {
        assert_metrics_file_eq(
            &standalone_prefix,
            &ground_truth_prefix,
            suffix,
            "duplex T2 (3-batch, 8 workers) vs ground truth",
        );
    }
}

#[test]
fn duplex_multi_key_batch_straddle_t2_matches_ground_truth() {
    let dir = tempfile::tempdir().expect("tempdir");
    let header = create_minimal_header("chr1", 10_000);
    let bam_path = dir.path().join("in.bam");

    let mut records = Vec::new();
    for k in 0..3 {
        // Well-separated positions so each `k` forms its own coordinate group.
        let r1_pos = 100 + i32::try_from(k).expect("k fits i32") * 2000;
        let r2_pos = r1_pos + 100;
        for f in 0..MULTI_KEY_FAMILIES_PER_KEY {
            let idx = k * MULTI_KEY_FAMILIES_PER_KEY + f;
            let umi = format!("{}-{}", indexed_umi(idx), indexed_umi(idx + 10_000));
            let (r1, r2) = duplex_pair(&format!("k{k}f{f}"), &umi, r1_pos, 10, r2_pos, 10);
            records.push(r1);
            records.push(r2);
        }
    }
    write_bam(&bam_path, &header, &records);

    let min_reads = [1usize];
    let (grouped, ground_truth_prefix) =
        duplex_ground_truth(dir.path(), &bam_path, &min_reads, None);
    let standalone_prefix =
        run_duplex_standalone(dir.path(), &grouped, "multikey", Some(8), &min_reads, None);

    for suffix in DUPLEX_SUFFIXES {
        assert_metrics_file_eq(
            &standalone_prefix,
            &ground_truth_prefix,
            suffix,
            "duplex T2 (multi-key batch straddle, 8 workers) vs ground truth",
        );
    }
}

// ============================================================================
// Codec analogs of the two parallel-T2 correctness anchors above (Task 5):
// many single-strand `duplex_pair` families at one coordinate key spanning
// three `GroupByMi` batches, and three distinct coordinate keys whose MI
// counts straddle those batch boundaries. Same rationale as the simplex/duplex
// cases: only the standalone T2 path re-derives coordinate-group boundaries
// from a batched MI stream, so these compare T2 against the separate-pass
// `duplex-metrics` ground truth only (codec shares the duplex-shaped metrics
// files, see `codec_ground_truth`/`DUPLEX_SUFFIXES` above).
// ============================================================================

#[test]
fn codec_three_batch_multi_worker_t2_matches_ground_truth() {
    let dir = tempfile::tempdir().expect("tempdir");
    let header = create_minimal_header("chr1", 10_000);
    let bam_path = dir.path().join("in.bam");

    let mut records = Vec::new();
    for f in 0..130 {
        let umi = format!("{}-{}", indexed_umi(f), indexed_umi(f + 10_000));
        let (r1, r2) = duplex_pair(&format!("f{f}"), &umi, 100, 10, 200, 10);
        records.push(r1);
        records.push(r2);
    }
    write_bam(&bam_path, &header, &records);

    let (grouped, ground_truth_prefix) = codec_ground_truth(dir.path(), &bam_path, 1, None);
    let standalone_prefix = run_codec_standalone(dir.path(), &grouped, "3batch", Some(8), 1, None);

    for suffix in DUPLEX_SUFFIXES {
        assert_metrics_file_eq(
            &standalone_prefix,
            &ground_truth_prefix,
            suffix,
            "codec T2 (3-batch, 8 workers) vs ground truth",
        );
    }
}

#[test]
fn codec_multi_key_batch_straddle_t2_matches_ground_truth() {
    let dir = tempfile::tempdir().expect("tempdir");
    let header = create_minimal_header("chr1", 10_000);
    let bam_path = dir.path().join("in.bam");

    let mut records = Vec::new();
    for k in 0..3 {
        // Well-separated positions so each `k` forms its own coordinate group.
        let r1_pos = 100 + i32::try_from(k).expect("k fits i32") * 2000;
        let r2_pos = r1_pos + 100;
        for f in 0..MULTI_KEY_FAMILIES_PER_KEY {
            let idx = k * MULTI_KEY_FAMILIES_PER_KEY + f;
            let umi = format!("{}-{}", indexed_umi(idx), indexed_umi(idx + 10_000));
            let (r1, r2) = duplex_pair(&format!("k{k}f{f}"), &umi, r1_pos, 10, r2_pos, 10);
            records.push(r1);
            records.push(r2);
        }
    }
    write_bam(&bam_path, &header, &records);

    let (grouped, ground_truth_prefix) = codec_ground_truth(dir.path(), &bam_path, 1, None);
    let standalone_prefix =
        run_codec_standalone(dir.path(), &grouped, "multikey", Some(8), 1, None);

    for suffix in DUPLEX_SUFFIXES {
        assert_metrics_file_eq(
            &standalone_prefix,
            &ground_truth_prefix,
            suffix,
            "codec T2 (multi-key batch straddle, 8 workers) vs ground truth",
        );
    }
}

//! Pipechain combinator tests for the 20 `(start_from, stop_after)` pairs that
//! aren't already covered by the dedicated per-stop test files.
//!
//! Each test exercises one pair, building two outputs:
//!
//! - **v1**: a chain of standalone `fgumi <stage>` invocations covering the
//!   same span.
//! - **v2**: `fgumi runall --start-from S --stop-after T` over the same input.
//!
//! Outputs are compared via `fgumi compare bams --command <preset>` (or
//! byte-equality for `--stop-after fastq`).
//!
//! Fixtures are produced via `fgumi simulate` (single seed, small molecule
//! count) so each test stays fast and deterministic. All v1 / v2 invocations
//! pin `--threads 1` and `--compression-level 1` to keep byte/record output
//! deterministic.
//!
//! Pairs covered here (rows = `--start-from`, columns = `--stop-after`):
//!
//! | from \ to | extract | correct | fastq | zipper | sort | group | consensus | filter |
//! |-----------|---------|---------|-------|--------|------|-------|-----------|--------|
//! | extract   | (1)     | (1)     | YES   | (1)    | (1)  | (1)   | YES       | (1)    |
//! | correct   |         | (1)     | YES   | YES    | YES  | YES   | YES       | YES    |
//! | fastq     |         |         | (1)   | YES    | YES  | YES   | YES       | YES    |
//! | sort      |         |         |       |        | YES  | YES   | YES       | YES    |
//! | group     |         |         |       |        |      | YES   | YES       | YES    |
//!
//! `(1)` = covered elsewhere in `tests/integration/runall/test_*.rs`.
//! `YES` = covered here.

use super::pipechain_helpers::{
    assert_fastq_bytes_equal, compare_bams, run_bwa_mem, run_correct, run_extract, run_fastq,
    run_filter, run_group, run_runall, run_simplex_consensus, run_sort, run_zipper,
    simulate_correct_reads, simulate_fastq_reads, simulate_grouped_reads, simulate_mapped_reads,
};
use crate::helpers::references::check_alignment_prerequisites;

const NUM_MOLECULES: usize = 30;
const NUM_CORRECT_READS: usize = 50;

// ===========================================================================
// from extract
// ===========================================================================

#[test]
fn test_extract_to_fastq() {
    // No alignment needed for extract -> fastq, but we still need a reference
    // for `simulate fastq-reads`.  Use the tracked PhiX (skips locally if
    // missing).
    let Some(reference) = check_alignment_prerequisites() else {
        return;
    };
    let tmp = tempfile::tempdir().unwrap();
    let (r1, r2) = simulate_fastq_reads(tmp.path(), NUM_MOLECULES, &reference);
    let v1_extract = tmp.path().join("v1_extract.bam");
    let v1_fastq = tmp.path().join("v1.fastq");
    let v2_fastq = tmp.path().join("v2.fastq");

    // v1: extract -> fastq (stdout)
    run_extract(&r1, &r2, &v1_extract);
    run_fastq(&v1_extract, &v1_fastq);

    // v2: runall extract..fastq
    run_runall(&[
        "--start-from",
        "extract",
        "--stop-after",
        "fastq",
        "--input",
        r1.to_str().unwrap(),
        r2.to_str().unwrap(),
        "--output",
        v2_fastq.to_str().unwrap(),
        "--extract::sample",
        "S",
        "--extract::library",
        "L",
        "--extract::read-structures",
        "8M+T",
        "+T",
        "--threads",
        "1",
        "--compression-level",
        "1",
    ]);

    assert_fastq_bytes_equal(&v1_fastq, &v2_fastq);
}

#[test]
fn test_extract_to_consensus() {
    let Some(reference) = check_alignment_prerequisites() else {
        return;
    };
    let tmp = tempfile::tempdir().unwrap();
    let (r1, r2) = simulate_fastq_reads(tmp.path(), NUM_MOLECULES, &reference);
    let v1_out = tmp.path().join("v1.bam");
    let v2_out = tmp.path().join("v2.bam");

    // v1: extract -> fastq -> bwa -> zipper -> sort -> group -> simplex
    let unmapped = tmp.path().join("v1_unmapped.bam");
    let fastq = tmp.path().join("v1_interleaved.fastq");
    let mapped = tmp.path().join("v1_mapped.sam");
    let zipped = tmp.path().join("v1_zipped.bam");
    let sorted = tmp.path().join("v1_sorted.bam");
    let grouped = tmp.path().join("v1_grouped.bam");
    run_extract(&r1, &r2, &unmapped);
    run_fastq(&unmapped, &fastq);
    run_bwa_mem(&reference, &fastq, &mapped);
    run_zipper(&unmapped, &mapped, &reference, &zipped);
    run_sort(&zipped, &sorted);
    run_group(&sorted, "adjacency", &grouped);
    run_simplex_consensus(&grouped, &v1_out, 1);

    // v2: runall extract..consensus
    run_runall(&[
        "--start-from",
        "extract",
        "--stop-after",
        "consensus",
        "--input",
        r1.to_str().unwrap(),
        r2.to_str().unwrap(),
        "--output",
        v2_out.to_str().unwrap(),
        "--extract::sample",
        "S",
        "--extract::library",
        "L",
        "--extract::read-structures",
        "8M+T",
        "+T",
        "--reference",
        reference.to_str().unwrap(),
        "--aligner::preset",
        "bwa-mem",
        "--aligner::threads",
        "1",
        "--consensus",
        "simplex",
        "--group::strategy",
        "adjacency",
        "--filter::min-reads",
        "1",
        "--threads",
        "1",
        "--compression-level",
        "1",
    ]);

    assert!(compare_bams(&v1_out, &v2_out, "simplex"), "extract->consensus mismatch");
}

// ===========================================================================
// from correct
// ===========================================================================

#[test]
fn test_correct_to_fastq() {
    // simulate correct-reads doesn't need a reference; --stop-after fastq
    // doesn't either.
    let tmp = tempfile::tempdir().unwrap();
    let (input_bam, includelist) = simulate_correct_reads(tmp.path(), NUM_CORRECT_READS);
    let v1_corrected = tmp.path().join("v1_corrected.bam");
    let v1_fastq = tmp.path().join("v1.fastq");
    let v2_fastq = tmp.path().join("v2.fastq");

    // v1: correct -> fastq
    run_correct(&input_bam, &includelist, &v1_corrected);
    run_fastq(&v1_corrected, &v1_fastq);

    // v2: runall correct..fastq
    run_runall(&[
        "--start-from",
        "correct",
        "--stop-after",
        "fastq",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        v2_fastq.to_str().unwrap(),
        "--correct::umi-files",
        includelist.to_str().unwrap(),
        "--correct::max-mismatches",
        "2",
        "--correct::min-distance",
        "2",
        "--threads",
        "1",
        "--compression-level",
        "1",
    ]);

    assert_fastq_bytes_equal(&v1_fastq, &v2_fastq);
}

#[test]
fn test_correct_to_zipper() {
    let Some(reference) = check_alignment_prerequisites() else {
        return;
    };
    let tmp = tempfile::tempdir().unwrap();
    let (input_bam, includelist) = simulate_correct_reads(tmp.path(), NUM_CORRECT_READS);
    let v1_out = tmp.path().join("v1.bam");
    let v2_out = tmp.path().join("v2.bam");

    // v1: correct -> fastq -> bwa -> zipper
    let corrected = tmp.path().join("v1_corrected.bam");
    let fastq = tmp.path().join("v1_interleaved.fastq");
    let mapped = tmp.path().join("v1_mapped.sam");
    run_correct(&input_bam, &includelist, &corrected);
    run_fastq(&corrected, &fastq);
    run_bwa_mem(&reference, &fastq, &mapped);
    run_zipper(&corrected, &mapped, &reference, &v1_out);

    // v2: runall correct..zipper
    run_runall(&[
        "--start-from",
        "correct",
        "--stop-after",
        "zipper",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        v2_out.to_str().unwrap(),
        "--correct::umi-files",
        includelist.to_str().unwrap(),
        "--correct::max-mismatches",
        "2",
        "--correct::min-distance",
        "2",
        "--reference",
        reference.to_str().unwrap(),
        "--aligner::preset",
        "bwa-mem",
        "--aligner::threads",
        "1",
        "--threads",
        "1",
        "--compression-level",
        "1",
    ]);

    assert!(compare_bams(&v1_out, &v2_out, "zipper"), "correct->zipper mismatch");
}

#[test]
fn test_correct_to_sort() {
    let Some(reference) = check_alignment_prerequisites() else {
        return;
    };
    let tmp = tempfile::tempdir().unwrap();
    let (input_bam, includelist) = simulate_correct_reads(tmp.path(), NUM_CORRECT_READS);
    let v1_out = tmp.path().join("v1.bam");
    let v2_out = tmp.path().join("v2.bam");

    // v1: correct -> fastq -> bwa -> zipper -> sort
    let corrected = tmp.path().join("v1_corrected.bam");
    let fastq = tmp.path().join("v1_interleaved.fastq");
    let mapped = tmp.path().join("v1_mapped.sam");
    let zipped = tmp.path().join("v1_zipped.bam");
    run_correct(&input_bam, &includelist, &corrected);
    run_fastq(&corrected, &fastq);
    run_bwa_mem(&reference, &fastq, &mapped);
    run_zipper(&corrected, &mapped, &reference, &zipped);
    run_sort(&zipped, &v1_out);

    // v2: runall correct..sort
    run_runall(&[
        "--start-from",
        "correct",
        "--stop-after",
        "sort",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        v2_out.to_str().unwrap(),
        "--correct::umi-files",
        includelist.to_str().unwrap(),
        "--correct::max-mismatches",
        "2",
        "--correct::min-distance",
        "2",
        "--reference",
        reference.to_str().unwrap(),
        "--aligner::preset",
        "bwa-mem",
        "--aligner::threads",
        "1",
        "--threads",
        "1",
        "--compression-level",
        "1",
    ]);

    assert!(compare_bams(&v1_out, &v2_out, "sort"), "correct->sort mismatch");
}

#[test]
fn test_correct_to_group() {
    let Some(reference) = check_alignment_prerequisites() else {
        return;
    };
    let tmp = tempfile::tempdir().unwrap();
    let (input_bam, includelist) = simulate_correct_reads(tmp.path(), NUM_CORRECT_READS);
    let v1_out = tmp.path().join("v1.bam");
    let v2_out = tmp.path().join("v2.bam");

    // v1: correct -> fastq -> bwa -> zipper -> sort -> group
    let corrected = tmp.path().join("v1_corrected.bam");
    let fastq = tmp.path().join("v1_interleaved.fastq");
    let mapped = tmp.path().join("v1_mapped.sam");
    let zipped = tmp.path().join("v1_zipped.bam");
    let sorted = tmp.path().join("v1_sorted.bam");
    run_correct(&input_bam, &includelist, &corrected);
    run_fastq(&corrected, &fastq);
    run_bwa_mem(&reference, &fastq, &mapped);
    run_zipper(&corrected, &mapped, &reference, &zipped);
    run_sort(&zipped, &sorted);
    run_group(&sorted, "adjacency", &v1_out);

    // v2: runall correct..group
    run_runall(&[
        "--start-from",
        "correct",
        "--stop-after",
        "group",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        v2_out.to_str().unwrap(),
        "--correct::umi-files",
        includelist.to_str().unwrap(),
        "--correct::max-mismatches",
        "2",
        "--correct::min-distance",
        "2",
        "--reference",
        reference.to_str().unwrap(),
        "--aligner::preset",
        "bwa-mem",
        "--aligner::threads",
        "1",
        "--group::strategy",
        "adjacency",
        "--threads",
        "1",
        "--compression-level",
        "1",
    ]);

    assert!(compare_bams(&v1_out, &v2_out, "group"), "correct->group mismatch");
}

#[test]
fn test_correct_to_consensus() {
    let Some(reference) = check_alignment_prerequisites() else {
        return;
    };
    let tmp = tempfile::tempdir().unwrap();
    let (input_bam, includelist) = simulate_correct_reads(tmp.path(), NUM_CORRECT_READS);
    let v1_out = tmp.path().join("v1.bam");
    let v2_out = tmp.path().join("v2.bam");

    // v1: correct -> fastq -> bwa -> zipper -> sort -> group -> simplex
    let corrected = tmp.path().join("v1_corrected.bam");
    let fastq = tmp.path().join("v1_interleaved.fastq");
    let mapped = tmp.path().join("v1_mapped.sam");
    let zipped = tmp.path().join("v1_zipped.bam");
    let sorted = tmp.path().join("v1_sorted.bam");
    let grouped = tmp.path().join("v1_grouped.bam");
    run_correct(&input_bam, &includelist, &corrected);
    run_fastq(&corrected, &fastq);
    run_bwa_mem(&reference, &fastq, &mapped);
    run_zipper(&corrected, &mapped, &reference, &zipped);
    run_sort(&zipped, &sorted);
    run_group(&sorted, "adjacency", &grouped);
    run_simplex_consensus(&grouped, &v1_out, 1);

    // v2: runall correct..consensus
    run_runall(&[
        "--start-from",
        "correct",
        "--stop-after",
        "consensus",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        v2_out.to_str().unwrap(),
        "--correct::umi-files",
        includelist.to_str().unwrap(),
        "--correct::max-mismatches",
        "2",
        "--correct::min-distance",
        "2",
        "--reference",
        reference.to_str().unwrap(),
        "--aligner::preset",
        "bwa-mem",
        "--aligner::threads",
        "1",
        "--consensus",
        "simplex",
        "--group::strategy",
        "adjacency",
        "--filter::min-reads",
        "1",
        "--threads",
        "1",
        "--compression-level",
        "1",
    ]);

    assert!(compare_bams(&v1_out, &v2_out, "simplex"), "correct->consensus mismatch");
}

#[test]
fn test_correct_to_filter() {
    let Some(reference) = check_alignment_prerequisites() else {
        return;
    };
    let tmp = tempfile::tempdir().unwrap();
    let (input_bam, includelist) = simulate_correct_reads(tmp.path(), NUM_CORRECT_READS);
    let v1_out = tmp.path().join("v1.bam");
    let v2_out = tmp.path().join("v2.bam");

    // v1: correct -> fastq -> bwa -> zipper -> sort -> group -> simplex -> filter
    let corrected = tmp.path().join("v1_corrected.bam");
    let fastq = tmp.path().join("v1_interleaved.fastq");
    let mapped = tmp.path().join("v1_mapped.sam");
    let zipped = tmp.path().join("v1_zipped.bam");
    let sorted = tmp.path().join("v1_sorted.bam");
    let grouped = tmp.path().join("v1_grouped.bam");
    let consensus = tmp.path().join("v1_consensus.bam");
    run_correct(&input_bam, &includelist, &corrected);
    run_fastq(&corrected, &fastq);
    run_bwa_mem(&reference, &fastq, &mapped);
    run_zipper(&corrected, &mapped, &reference, &zipped);
    run_sort(&zipped, &sorted);
    run_group(&sorted, "adjacency", &grouped);
    run_simplex_consensus(&grouped, &consensus, 1);
    run_filter(&consensus, &v1_out);

    // v2: runall correct..filter
    run_runall(&[
        "--start-from",
        "correct",
        "--stop-after",
        "filter",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        v2_out.to_str().unwrap(),
        "--correct::umi-files",
        includelist.to_str().unwrap(),
        "--correct::max-mismatches",
        "2",
        "--correct::min-distance",
        "2",
        "--reference",
        reference.to_str().unwrap(),
        "--aligner::preset",
        "bwa-mem",
        "--aligner::threads",
        "1",
        "--consensus",
        "simplex",
        "--group::strategy",
        "adjacency",
        "--filter::min-base-quality",
        "40",
        "--filter::max-read-error-rate",
        "0.025",
        "--filter::max-base-error-rate",
        "0.1",
        "--filter::max-no-call-fraction",
        "0.2",
        "--filter::min-reads",
        "1",
        "--threads",
        "1",
        "--compression-level",
        "1",
    ]);

    assert!(compare_bams(&v1_out, &v2_out, "filter"), "correct->filter mismatch");
}

// ===========================================================================
// from fastq
//
// `runall --start-from fastq` consumes an unmapped BAM (extract output), not
// raw FASTQ.  Each test below pre-runs `fgumi extract` once to build that
// fixture, then feeds the same BAM into both v1 and v2.
// ===========================================================================

/// Build the unmapped-BAM fixture used by every `from fastq` test.
fn make_extracted_fixture(tmp_dir: &std::path::Path, reference: &std::path::Path) -> PathBuf {
    let (r1, r2) = simulate_fastq_reads(tmp_dir, NUM_MOLECULES, reference);
    let extracted = tmp_dir.join("fixture_extracted.bam");
    run_extract(&r1, &r2, &extracted);
    extracted
}

use std::path::PathBuf;

#[test]
fn test_fastq_to_zipper() {
    let Some(reference) = check_alignment_prerequisites() else {
        return;
    };
    let tmp = tempfile::tempdir().unwrap();
    let extracted = make_extracted_fixture(tmp.path(), &reference);
    let v1_out = tmp.path().join("v1.bam");
    let v2_out = tmp.path().join("v2.bam");

    // v1: fastq -> bwa -> zipper
    let fastq = tmp.path().join("v1_interleaved.fastq");
    let mapped = tmp.path().join("v1_mapped.sam");
    run_fastq(&extracted, &fastq);
    run_bwa_mem(&reference, &fastq, &mapped);
    run_zipper(&extracted, &mapped, &reference, &v1_out);

    // v2: runall fastq..zipper
    run_runall(&[
        "--start-from",
        "fastq",
        "--stop-after",
        "zipper",
        "--input",
        extracted.to_str().unwrap(),
        "--output",
        v2_out.to_str().unwrap(),
        "--reference",
        reference.to_str().unwrap(),
        "--aligner::preset",
        "bwa-mem",
        "--aligner::threads",
        "1",
        "--threads",
        "1",
        "--compression-level",
        "1",
    ]);

    assert!(compare_bams(&v1_out, &v2_out, "zipper"), "fastq->zipper mismatch");
}

#[test]
fn test_fastq_to_sort() {
    let Some(reference) = check_alignment_prerequisites() else {
        return;
    };
    let tmp = tempfile::tempdir().unwrap();
    let extracted = make_extracted_fixture(tmp.path(), &reference);
    let v1_out = tmp.path().join("v1.bam");
    let v2_out = tmp.path().join("v2.bam");

    let fastq = tmp.path().join("v1_interleaved.fastq");
    let mapped = tmp.path().join("v1_mapped.sam");
    let zipped = tmp.path().join("v1_zipped.bam");
    run_fastq(&extracted, &fastq);
    run_bwa_mem(&reference, &fastq, &mapped);
    run_zipper(&extracted, &mapped, &reference, &zipped);
    run_sort(&zipped, &v1_out);

    run_runall(&[
        "--start-from",
        "fastq",
        "--stop-after",
        "sort",
        "--input",
        extracted.to_str().unwrap(),
        "--output",
        v2_out.to_str().unwrap(),
        "--reference",
        reference.to_str().unwrap(),
        "--aligner::preset",
        "bwa-mem",
        "--aligner::threads",
        "1",
        "--threads",
        "1",
        "--compression-level",
        "1",
    ]);

    assert!(compare_bams(&v1_out, &v2_out, "sort"), "fastq->sort mismatch");
}

#[test]
fn test_fastq_to_group() {
    let Some(reference) = check_alignment_prerequisites() else {
        return;
    };
    let tmp = tempfile::tempdir().unwrap();
    let extracted = make_extracted_fixture(tmp.path(), &reference);
    let v1_out = tmp.path().join("v1.bam");
    let v2_out = tmp.path().join("v2.bam");

    let fastq = tmp.path().join("v1_interleaved.fastq");
    let mapped = tmp.path().join("v1_mapped.sam");
    let zipped = tmp.path().join("v1_zipped.bam");
    let sorted = tmp.path().join("v1_sorted.bam");
    run_fastq(&extracted, &fastq);
    run_bwa_mem(&reference, &fastq, &mapped);
    run_zipper(&extracted, &mapped, &reference, &zipped);
    run_sort(&zipped, &sorted);
    run_group(&sorted, "adjacency", &v1_out);

    run_runall(&[
        "--start-from",
        "fastq",
        "--stop-after",
        "group",
        "--input",
        extracted.to_str().unwrap(),
        "--output",
        v2_out.to_str().unwrap(),
        "--reference",
        reference.to_str().unwrap(),
        "--aligner::preset",
        "bwa-mem",
        "--aligner::threads",
        "1",
        "--group::strategy",
        "adjacency",
        "--threads",
        "1",
        "--compression-level",
        "1",
    ]);

    assert!(compare_bams(&v1_out, &v2_out, "group"), "fastq->group mismatch");
}

#[test]
fn test_fastq_to_consensus() {
    let Some(reference) = check_alignment_prerequisites() else {
        return;
    };
    let tmp = tempfile::tempdir().unwrap();
    let extracted = make_extracted_fixture(tmp.path(), &reference);
    let v1_out = tmp.path().join("v1.bam");
    let v2_out = tmp.path().join("v2.bam");

    let fastq = tmp.path().join("v1_interleaved.fastq");
    let mapped = tmp.path().join("v1_mapped.sam");
    let zipped = tmp.path().join("v1_zipped.bam");
    let sorted = tmp.path().join("v1_sorted.bam");
    let grouped = tmp.path().join("v1_grouped.bam");
    run_fastq(&extracted, &fastq);
    run_bwa_mem(&reference, &fastq, &mapped);
    run_zipper(&extracted, &mapped, &reference, &zipped);
    run_sort(&zipped, &sorted);
    run_group(&sorted, "adjacency", &grouped);
    run_simplex_consensus(&grouped, &v1_out, 1);

    run_runall(&[
        "--start-from",
        "fastq",
        "--stop-after",
        "consensus",
        "--input",
        extracted.to_str().unwrap(),
        "--output",
        v2_out.to_str().unwrap(),
        "--reference",
        reference.to_str().unwrap(),
        "--aligner::preset",
        "bwa-mem",
        "--aligner::threads",
        "1",
        "--consensus",
        "simplex",
        "--group::strategy",
        "adjacency",
        "--filter::min-reads",
        "1",
        "--threads",
        "1",
        "--compression-level",
        "1",
    ]);

    assert!(compare_bams(&v1_out, &v2_out, "simplex"), "fastq->consensus mismatch");
}

#[test]
fn test_fastq_to_filter() {
    let Some(reference) = check_alignment_prerequisites() else {
        return;
    };
    let tmp = tempfile::tempdir().unwrap();
    let extracted = make_extracted_fixture(tmp.path(), &reference);
    let v1_out = tmp.path().join("v1.bam");
    let v2_out = tmp.path().join("v2.bam");

    let fastq = tmp.path().join("v1_interleaved.fastq");
    let mapped = tmp.path().join("v1_mapped.sam");
    let zipped = tmp.path().join("v1_zipped.bam");
    let sorted = tmp.path().join("v1_sorted.bam");
    let grouped = tmp.path().join("v1_grouped.bam");
    let consensus = tmp.path().join("v1_consensus.bam");
    run_fastq(&extracted, &fastq);
    run_bwa_mem(&reference, &fastq, &mapped);
    run_zipper(&extracted, &mapped, &reference, &zipped);
    run_sort(&zipped, &sorted);
    run_group(&sorted, "adjacency", &grouped);
    run_simplex_consensus(&grouped, &consensus, 1);
    run_filter(&consensus, &v1_out);

    run_runall(&[
        "--start-from",
        "fastq",
        "--stop-after",
        "filter",
        "--input",
        extracted.to_str().unwrap(),
        "--output",
        v2_out.to_str().unwrap(),
        "--reference",
        reference.to_str().unwrap(),
        "--aligner::preset",
        "bwa-mem",
        "--aligner::threads",
        "1",
        "--consensus",
        "simplex",
        "--group::strategy",
        "adjacency",
        "--filter::min-base-quality",
        "40",
        "--filter::max-read-error-rate",
        "0.025",
        "--filter::max-base-error-rate",
        "0.1",
        "--filter::max-no-call-fraction",
        "0.2",
        "--filter::min-reads",
        "1",
        "--threads",
        "1",
        "--compression-level",
        "1",
    ]);

    assert!(compare_bams(&v1_out, &v2_out, "filter"), "fastq->filter mismatch");
}

// ===========================================================================
// from sort
//
// Input is a coordinate / template-coordinate sorted mapped BAM with `RX`
// tags.  We use `simulate mapped-reads` (which emits template-coord-sorted)
// for the fixture; runall's `sort` stage (re-)sorts and the v1 standalone
// chain does likewise for the `sort` stop, so behaviours match.
// ===========================================================================

#[test]
fn test_sort_to_sort() {
    let Some(reference) = check_alignment_prerequisites() else {
        return;
    };
    let tmp = tempfile::tempdir().unwrap();
    let input_bam = simulate_mapped_reads(tmp.path(), NUM_MOLECULES, &reference);
    let v1_out = tmp.path().join("v1.bam");
    let v2_out = tmp.path().join("v2.bam");

    // v1: just sort the input
    run_sort(&input_bam, &v1_out);

    // v2: runall sort..sort
    run_runall(&[
        "--start-from",
        "sort",
        "--stop-after",
        "sort",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        v2_out.to_str().unwrap(),
        "--reference",
        reference.to_str().unwrap(),
        "--threads",
        "1",
        "--compression-level",
        "1",
    ]);

    assert!(compare_bams(&v1_out, &v2_out, "sort"), "sort->sort mismatch");
}

#[test]
fn test_sort_to_group() {
    let Some(reference) = check_alignment_prerequisites() else {
        return;
    };
    let tmp = tempfile::tempdir().unwrap();
    let input_bam = simulate_mapped_reads(tmp.path(), NUM_MOLECULES, &reference);
    let v1_out = tmp.path().join("v1.bam");
    let v2_out = tmp.path().join("v2.bam");

    let sorted = tmp.path().join("v1_sorted.bam");
    run_sort(&input_bam, &sorted);
    run_group(&sorted, "adjacency", &v1_out);

    run_runall(&[
        "--start-from",
        "sort",
        "--stop-after",
        "group",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        v2_out.to_str().unwrap(),
        "--reference",
        reference.to_str().unwrap(),
        "--group::strategy",
        "adjacency",
        "--threads",
        "1",
        "--compression-level",
        "1",
    ]);

    assert!(compare_bams(&v1_out, &v2_out, "group"), "sort->group mismatch");
}

#[test]
fn test_sort_to_consensus() {
    let Some(reference) = check_alignment_prerequisites() else {
        return;
    };
    let tmp = tempfile::tempdir().unwrap();
    let input_bam = simulate_mapped_reads(tmp.path(), NUM_MOLECULES, &reference);
    let v1_out = tmp.path().join("v1.bam");
    let v2_out = tmp.path().join("v2.bam");

    let sorted = tmp.path().join("v1_sorted.bam");
    let grouped = tmp.path().join("v1_grouped.bam");
    run_sort(&input_bam, &sorted);
    run_group(&sorted, "adjacency", &grouped);
    run_simplex_consensus(&grouped, &v1_out, 1);

    run_runall(&[
        "--start-from",
        "sort",
        "--stop-after",
        "consensus",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        v2_out.to_str().unwrap(),
        "--reference",
        reference.to_str().unwrap(),
        "--consensus",
        "simplex",
        "--group::strategy",
        "adjacency",
        "--filter::min-reads",
        "1",
        "--threads",
        "1",
        "--compression-level",
        "1",
    ]);

    assert!(compare_bams(&v1_out, &v2_out, "simplex"), "sort->consensus mismatch");
}

#[test]
fn test_sort_to_filter() {
    let Some(reference) = check_alignment_prerequisites() else {
        return;
    };
    let tmp = tempfile::tempdir().unwrap();
    let input_bam = simulate_mapped_reads(tmp.path(), NUM_MOLECULES, &reference);
    let v1_out = tmp.path().join("v1.bam");
    let v2_out = tmp.path().join("v2.bam");

    let sorted = tmp.path().join("v1_sorted.bam");
    let grouped = tmp.path().join("v1_grouped.bam");
    let consensus = tmp.path().join("v1_consensus.bam");
    run_sort(&input_bam, &sorted);
    run_group(&sorted, "adjacency", &grouped);
    run_simplex_consensus(&grouped, &consensus, 1);
    run_filter(&consensus, &v1_out);

    run_runall(&[
        "--start-from",
        "sort",
        "--stop-after",
        "filter",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        v2_out.to_str().unwrap(),
        "--reference",
        reference.to_str().unwrap(),
        "--consensus",
        "simplex",
        "--group::strategy",
        "adjacency",
        "--filter::min-base-quality",
        "40",
        "--filter::max-read-error-rate",
        "0.025",
        "--filter::max-base-error-rate",
        "0.1",
        "--filter::max-no-call-fraction",
        "0.2",
        "--filter::min-reads",
        "1",
        "--threads",
        "1",
        "--compression-level",
        "1",
    ]);

    assert!(compare_bams(&v1_out, &v2_out, "filter"), "sort->filter mismatch");
}

// ===========================================================================
// from group
//
// Input is a template-coordinate-sorted mapped BAM with MI tags already
// assigned.  Per the runall CLI doc, `--start-from group` overwrites pre-
// existing MI tags, so the standalone `--stop-after group` chain is just
// `fgumi group --strategy adjacency` on the same input (which also re-assigns
// MI tags).
// ===========================================================================

#[test]
fn test_group_to_group() {
    let Some(reference) = check_alignment_prerequisites() else {
        return;
    };
    let tmp = tempfile::tempdir().unwrap();
    let input_bam = simulate_grouped_reads(tmp.path(), NUM_MOLECULES, &reference);
    let v1_out = tmp.path().join("v1.bam");
    let v2_out = tmp.path().join("v2.bam");

    // v1: re-group the input
    run_group(&input_bam, "adjacency", &v1_out);

    // v2: runall group..group (also re-groups)
    run_runall(&[
        "--start-from",
        "group",
        "--stop-after",
        "group",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        v2_out.to_str().unwrap(),
        "--reference",
        reference.to_str().unwrap(),
        "--group::strategy",
        "adjacency",
        "--threads",
        "1",
        "--compression-level",
        "1",
    ]);

    assert!(compare_bams(&v1_out, &v2_out, "group"), "group->group mismatch");
}

#[test]
fn test_group_to_consensus() {
    let Some(reference) = check_alignment_prerequisites() else {
        return;
    };
    let tmp = tempfile::tempdir().unwrap();
    let input_bam = simulate_grouped_reads(tmp.path(), NUM_MOLECULES, &reference);
    let v1_out = tmp.path().join("v1.bam");
    let v2_out = tmp.path().join("v2.bam");

    let grouped = tmp.path().join("v1_grouped.bam");
    run_group(&input_bam, "adjacency", &grouped);
    run_simplex_consensus(&grouped, &v1_out, 1);

    run_runall(&[
        "--start-from",
        "group",
        "--stop-after",
        "consensus",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        v2_out.to_str().unwrap(),
        "--reference",
        reference.to_str().unwrap(),
        "--consensus",
        "simplex",
        "--group::strategy",
        "adjacency",
        "--filter::min-reads",
        "1",
        "--threads",
        "1",
        "--compression-level",
        "1",
    ]);

    assert!(compare_bams(&v1_out, &v2_out, "simplex"), "group->consensus mismatch");
}

#[test]
fn test_group_to_filter() {
    let Some(reference) = check_alignment_prerequisites() else {
        return;
    };
    let tmp = tempfile::tempdir().unwrap();
    let input_bam = simulate_grouped_reads(tmp.path(), NUM_MOLECULES, &reference);
    let v1_out = tmp.path().join("v1.bam");
    let v2_out = tmp.path().join("v2.bam");

    let grouped = tmp.path().join("v1_grouped.bam");
    let consensus = tmp.path().join("v1_consensus.bam");
    run_group(&input_bam, "adjacency", &grouped);
    run_simplex_consensus(&grouped, &consensus, 1);
    run_filter(&consensus, &v1_out);

    run_runall(&[
        "--start-from",
        "group",
        "--stop-after",
        "filter",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        v2_out.to_str().unwrap(),
        "--reference",
        reference.to_str().unwrap(),
        "--consensus",
        "simplex",
        "--group::strategy",
        "adjacency",
        "--filter::min-base-quality",
        "40",
        "--filter::max-read-error-rate",
        "0.025",
        "--filter::max-base-error-rate",
        "0.1",
        "--filter::max-no-call-fraction",
        "0.2",
        "--filter::min-reads",
        "1",
        "--threads",
        "1",
        "--compression-level",
        "1",
    ]);

    assert!(compare_bams(&v1_out, &v2_out, "filter"), "group->filter mismatch");
}

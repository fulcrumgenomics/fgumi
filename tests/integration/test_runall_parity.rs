//! Runall parity tests (#33).
//!
//! Verifies that each accepted `(--start-from, --stop-after)` pair
//! produces output equivalent to the corresponding non-fused
//! invocation:
//!
//!   * **Class A — single-stage parity (5 tests).** Runs
//!     `fgumi runall --start-from X --stop-after X` and `fgumi X` on
//!     the same input, asserts the output BAMs contain the same
//!     record stream. Pins that runall's per-stage execution path
//!     matches the standalone command's execution path. All five run
//!     by default.
//!
//!   * **Class B — composition / multi-stage parity (7 tests).** Runs
//!     `fgumi runall --start-from S --stop-after T` and the equivalent
//!     staged chain (`fgumi S → tmp.bam → fgumi ... → fgumi T`),
//!     asserts the output BAMs contain the same record stream. Pins
//!     that runall's FUSED chain matches the standalone STAGED chain.
//!     All seven (`{Sort,Group}→{Simplex,Duplex,Codec}` plus
//!     `Sort→Group`) run by default.
//!
//! Total: 12 parity tests covering every valid `(start, stop)` pair
//! the validator accepts (5 diagonal + 7 off-diagonal = 12; the
//! remaining 13 cells of the 5×5 matrix are structurally invalid per
//! `RunAllStage::validate_with` and have CLI-level rejection tests at
//! the bottom of this file).

#![allow(clippy::needless_pass_by_value)]
// The non-consensus parity oracles (sort/group/zipper/correct/align) run in
// every build; the consensus-targeting cases (`*_simplex/duplex/codec*`) are
// gated per-test with `#[cfg(feature = "consensus")]` below. In a reduced
// (`--no-default-features`) build those gated tests compile out, leaving the
// consensus-only fixtures and helper imports they used unused — silence the
// resulting warnings in that build only (they are all live under `consensus`).
#![cfg_attr(not(feature = "consensus"), allow(dead_code, unused_imports))]

use std::ffi::OsString;
use std::fs;
use std::path::{Path, PathBuf};

use noodles::bam;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use rstest::rstest;
use tempfile::TempDir;

use crate::helpers::bam_generator::{
    create_both_unmapped_pair, create_minimal_header, create_paired_umi_family,
    create_paired_umi_family_at, create_umi_family, create_umi_family_at, to_record_buf,
};
use crate::helpers::cli_runner::{
    ParityArgs, Stage, fgumi, fgumi_binary, run_runall, run_runall_consensus_to_filter,
    run_staged_chain, run_standalone, run_standalone_filter,
};
use crate::helpers::parity::{
    assert_bam_headers_equivalent_ignoring_pg, assert_bams_record_equivalent,
    assert_bams_record_equivalent_nonempty, read_bam_records,
};

// ────────────────────────── Fixtures ──────────────────────────
//
// Three fixture families cover the 12 combos:
//
//   * **Simplex** (`unsorted_simplex_fixture`, `sorted_simplex_fixture`,
//     `grouped_simplex_fixture`) — single-end aligned BAM with single
//     RX tag per family. Used by all Simplex-targeting tests and by
//     Sort→Sort / Group→Group.
//   * **Duplex** (`unsorted_duplex_fixture`, `sorted_duplex_fixture`,
//     `grouped_duplex_fixture`) — paired-end aligned BAM with paired
//     UMI (`AAAA-CCCC`) and MC tags. Used by all Duplex-targeting
//     tests (`--strategy paired`).
//   * **Codec** (`unsorted_codec_fixture`, `sorted_codec_fixture`,
//     `grouped_codec_fixture`) — paired-end aligned BAM with a single
//     RX tag per family + MC tags, grouped with `--strategy adjacency`.
//     CODEC's protocol is "R1 reads one strand, R2 reads the other",
//     so the right input is *paired-end, single-UMI* — not the
//     simplex single-end fixture. Used by all Codec-targeting tests.
//
// Each family is built by chaining the standalone commands at fixture-
// build time (sorted = unsorted → sort, grouped = sorted → group). The
// chain is small and deterministic, so re-running it per test costs
// well under a second.

fn unsorted_simplex_fixture(dir: &Path) -> PathBuf {
    let path = dir.join("unsorted_simplex.bam");
    let header = create_minimal_header("chr1", 10000);
    let mut writer = bam::io::Writer::new(fs::File::create(&path).expect("create fixture BAM"));
    writer.write_header(&header).expect("write fixture header");
    // 4 families × 4 reads each at distinct positions on chr1; the
    // file is written in family-interleaved order so it is not
    // sorted by template-coordinate.
    let families = [
        ("ACGTACGT", "fam_a", "ACGTACGTACGTACGT"),
        ("TGCATGCA", "fam_b", "TGCATGCATGCATGCA"),
        ("CCAATTGG", "fam_c", "CCAATTGGCCAATTGG"),
        ("GGTTAACC", "fam_d", "GGTTAACCGGTTAACC"),
    ];
    let mut all_records = Vec::new();
    for (umi, name, seq) in families {
        for r in create_umi_family(umi, 4, name, seq, 30) {
            all_records.push(r);
        }
    }
    // Shuffle (deterministic): reverse so the file order doesn't
    // match template-coordinate.
    all_records.reverse();
    for raw in &all_records {
        writer.write_alignment_record(&header, &to_record_buf(raw)).expect("write fixture record");
    }
    writer.try_finish().expect("finish fixture BAM");
    path
}

fn sorted_simplex_fixture(dir: &Path) -> PathBuf {
    let unsorted = unsorted_simplex_fixture(dir);
    let sorted = dir.join("sorted_simplex.bam");
    let args = ParityArgs::for_simplex_like();
    let output = run_standalone(Stage::Sort, &unsorted, &sorted, &args);
    assert!(
        output.status.success(),
        "sorted_simplex_fixture: sort failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    sorted
}

fn grouped_simplex_fixture(dir: &Path) -> PathBuf {
    let sorted = sorted_simplex_fixture(dir);
    let grouped = dir.join("grouped_simplex.bam");
    let args = ParityArgs::for_simplex_like();
    let output = run_standalone(Stage::Group, &sorted, &grouped, &args);
    assert!(
        output.status.success(),
        "grouped_simplex_fixture: group failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    grouped
}

fn unsorted_duplex_fixture(dir: &Path) -> PathBuf {
    let path = dir.join("unsorted_duplex.bam");
    let header = create_minimal_header("chr1", 10000);
    let mut writer = bam::io::Writer::new(fs::File::create(&path).expect("create fixture BAM"));
    writer.write_header(&header).expect("write fixture header");
    // 4 paired-UMI families × 4 read-pairs each. Paired UMI format
    // ("AAAA-CCCC") is what the `paired` strategy expects.
    let families = [
        ("AAAA-CCCC", "fam_a", "ACGTACGT", "TTTTAAAA"),
        ("GGGG-TTTT", "fam_b", "TGCATGCA", "CCCCGGGG"),
        ("CCAA-TTGG", "fam_c", "CCAATTGG", "GGTTAACC"),
        ("ATCG-GCTA", "fam_d", "GGTTAACC", "CCAATTGG"),
    ];
    let mut all_records = Vec::new();
    for (umi, name, r1_seq, r2_seq) in families {
        for r in create_paired_umi_family(umi, 4, name, r1_seq, r2_seq, 30) {
            all_records.push(r);
        }
    }
    all_records.reverse();
    for raw in &all_records {
        writer.write_alignment_record(&header, &to_record_buf(raw)).expect("write fixture record");
    }
    writer.try_finish().expect("finish fixture BAM");
    path
}

fn sorted_duplex_fixture(dir: &Path) -> PathBuf {
    let unsorted = unsorted_duplex_fixture(dir);
    let sorted = dir.join("sorted_duplex.bam");
    let args = ParityArgs::for_duplex();
    let output = run_standalone(Stage::Sort, &unsorted, &sorted, &args);
    assert!(
        output.status.success(),
        "sorted_duplex_fixture: sort failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    sorted
}

fn grouped_duplex_fixture(dir: &Path) -> PathBuf {
    let sorted = sorted_duplex_fixture(dir);
    let grouped = dir.join("grouped_duplex.bam");
    let args = ParityArgs::for_duplex();
    let output = run_standalone(Stage::Group, &sorted, &grouped, &args);
    assert!(
        output.status.success(),
        "grouped_duplex_fixture: group failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    grouped
}

/// Template-coordinate-sorted duplex fixture that ALSO contains one
/// both-unmapped read pair (a valid paired UMI, both mates `UNMAPPED`).
///
/// Built so the both-unmapped pair survives `group --allow-unmapped` but is
/// then dropped by the duplex `record_filter` on the non-fused path. Used by
/// [`new_002_group_to_duplex_allow_unmapped_parity`] to pin that the fused
/// group→duplex bridge applies the same filter (NEW-002).
fn sorted_duplex_with_unmapped_fixture(dir: &Path) -> PathBuf {
    let unsorted = dir.join("unsorted_duplex_with_unmapped.bam");
    let header = create_minimal_header("chr1", 10000);
    let mut writer = bam::io::Writer::new(fs::File::create(&unsorted).expect("create fixture BAM"));
    writer.write_header(&header).expect("write fixture header");

    // The standard 4 fully-mapped paired-UMI families.
    let families = [
        ("AAAA-CCCC", "fam_a", "ACGTACGT", "TTTTAAAA"),
        ("GGGG-TTTT", "fam_b", "TGCATGCA", "CCCCGGGG"),
        ("CCAA-TTGG", "fam_c", "CCAATTGG", "GGTTAACC"),
        ("ATCG-GCTA", "fam_d", "GGTTAACC", "CCAATTGG"),
    ];
    let mut all_records = Vec::new();
    for (umi, name, r1_seq, r2_seq) in families {
        for r in create_paired_umi_family(umi, 4, name, r1_seq, r2_seq, 30) {
            all_records.push(r);
        }
    }
    // One both-unmapped pair with its own paired UMI.
    for r in create_both_unmapped_pair("TTAA-GGCC", "fam_unmapped", "ACGTACGT", "TTTTAAAA", 30) {
        all_records.push(r);
    }
    all_records.reverse();
    for raw in &all_records {
        writer.write_alignment_record(&header, &to_record_buf(raw)).expect("write fixture record");
    }
    writer.try_finish().expect("finish fixture BAM");

    let sorted = dir.join("sorted_duplex_with_unmapped.bam");
    let args = ParityArgs::for_duplex();
    let output = run_standalone(Stage::Sort, &unsorted, &sorted, &args);
    assert!(
        output.status.success(),
        "sorted_duplex_with_unmapped_fixture: sort failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    sorted
}

fn unsorted_codec_fixture(dir: &Path) -> PathBuf {
    let path = dir.join("unsorted_codec.bam");
    let header = create_minimal_header("chr1", 10000);
    let mut writer = bam::io::Writer::new(fs::File::create(&path).expect("create fixture BAM"));
    writer.write_header(&header).expect("write fixture header");
    // 4 paired-end families × 4 read-pairs. Each pair shares one
    // single-strand UMI in the RX tag; CODEC pairs the two strands at
    // consensus time, not at grouping time. `--strategy adjacency`
    // (not `paired`) is the documented grouping for CODEC.
    let families = [
        ("ACGTACGT", "fam_a", "ACGTACGT", "TTTTAAAA"),
        ("TGCATGCA", "fam_b", "TGCATGCA", "CCCCGGGG"),
        ("CCAATTGG", "fam_c", "CCAATTGG", "GGTTAACC"),
        ("GGTTAACC", "fam_d", "GGTTAACC", "CCAATTGG"),
    ];
    let mut all_records = Vec::new();
    for (umi, name, r1_seq, r2_seq) in families {
        for r in create_paired_umi_family(umi, 4, name, r1_seq, r2_seq, 30) {
            all_records.push(r);
        }
    }
    all_records.reverse();
    for raw in &all_records {
        writer.write_alignment_record(&header, &to_record_buf(raw)).expect("write fixture record");
    }
    writer.try_finish().expect("finish fixture BAM");
    path
}

fn sorted_codec_fixture(dir: &Path) -> PathBuf {
    let unsorted = unsorted_codec_fixture(dir);
    let sorted = dir.join("sorted_codec.bam");
    let args = ParityArgs::for_simplex_like();
    let output = run_standalone(Stage::Sort, &unsorted, &sorted, &args);
    assert!(
        output.status.success(),
        "sorted_codec_fixture: sort failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    sorted
}

fn grouped_codec_fixture(dir: &Path) -> PathBuf {
    let sorted = sorted_codec_fixture(dir);
    let grouped = dir.join("grouped_codec.bam");
    let args = ParityArgs::for_simplex_like();
    let output = run_standalone(Stage::Group, &sorted, &grouped, &args);
    assert!(
        output.status.success(),
        "grouped_codec_fixture: group failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    grouped
}

// ────────────────────────── Class A — single-stage parity ──────────────────────────
//
// Each Class A test runs `fgumi runall --start-from X --stop-after X`
// against `fgumi X` on the same input.

/// `Sort → Sort` parity vs standalone `fgumi sort`. Runall delegates to
/// `Sort::execute`, so the two are record- and header-equivalent — identical
/// records and `@HD`/`@SQ`/`@RG`, ignoring only the `@PG` command-line
/// provenance (which legitimately differs between runall and standalone).
#[test]
fn parity_a_sort_to_sort() {
    let tmp = TempDir::new().unwrap();
    let fixture = unsorted_simplex_fixture(tmp.path());
    let runall_out = tmp.path().join("runall.bam");
    let standalone_out = tmp.path().join("standalone.bam");
    let args = ParityArgs::for_simplex_like();

    let r = run_runall(Stage::Sort, Stage::Sort, &fixture, &runall_out, &args);
    assert!(r.status.success(), "runall: {}", String::from_utf8_lossy(&r.stderr));

    let r = run_standalone(Stage::Sort, &fixture, &standalone_out, &args);
    assert!(r.status.success(), "standalone: {}", String::from_utf8_lossy(&r.stderr));

    assert_bams_record_equivalent_nonempty(&runall_out, &standalone_out);
    // S9a-005: the sort-order metadata (`@HD SO/GO/SS`), `@SQ` dictionary, and
    // `@RG` set must match too — a fused runall that emitted the wrong sort
    // order would be invisible to the record-only equivalence check.
    assert_bam_headers_equivalent_ignoring_pg(&runall_out, &standalone_out);
}

/// `Group → Group` parity vs standalone `fgumi group`. Runall delegates to
/// `GroupReadsByUmi::execute`, so the two are record- and header-equivalent —
/// identical records and `@HD`/`@SQ`/`@RG`, ignoring only the `@PG`
/// command-line provenance (which legitimately differs).
#[test]
fn parity_a_group_to_group() {
    let tmp = TempDir::new().unwrap();
    let fixture = sorted_simplex_fixture(tmp.path());
    let runall_out = tmp.path().join("runall.bam");
    let standalone_out = tmp.path().join("standalone.bam");
    let args = ParityArgs::for_simplex_like();

    let r = run_runall(Stage::Group, Stage::Group, &fixture, &runall_out, &args);
    assert!(r.status.success(), "runall: {}", String::from_utf8_lossy(&r.stderr));

    let r = run_standalone(Stage::Group, &fixture, &standalone_out, &args);
    assert!(r.status.success(), "standalone: {}", String::from_utf8_lossy(&r.stderr));

    assert_bams_record_equivalent_nonempty(&runall_out, &standalone_out);
    // S9a-005: group preserves the sort metadata + dictionary + read groups.
    assert_bam_headers_equivalent_ignoring_pg(&runall_out, &standalone_out);
}

/// `Simplex → Simplex` parity vs standalone `fgumi simplex`.
/// Uses the consensus-only delegation fast path
/// (`execute_consensus_only`) —
/// runall constructs a `Simplex` struct from its own flags and calls
/// `Simplex::execute` directly, skipping the group step entirely.
/// Parity therefore holds by construction (same standalone command,
/// same input); this test pins the delegation glue (input forwarding,
/// option-struct construction).
#[cfg(feature = "consensus")]
#[test]
fn parity_a_simplex_to_simplex() {
    let tmp = TempDir::new().unwrap();
    let fixture = grouped_simplex_fixture(tmp.path());
    let runall_out = tmp.path().join("runall.bam");
    let standalone_out = tmp.path().join("standalone.bam");
    let args = ParityArgs::for_simplex_like();

    let r = run_runall(Stage::Simplex, Stage::Simplex, &fixture, &runall_out, &args);
    assert!(r.status.success(), "runall: {}", String::from_utf8_lossy(&r.stderr));

    let r = run_standalone(Stage::Simplex, &fixture, &standalone_out, &args);
    assert!(r.status.success(), "standalone: {}", String::from_utf8_lossy(&r.stderr));

    // S9a-001: simplex on the grouped fixture emits one fragment consensus per
    // MI group, so the stream is non-empty — pin that floor so a regression
    // that dropped all consensus records can't pass vacuously.
    assert_bams_record_equivalent_nonempty(&runall_out, &standalone_out);
    // S9a-005: the consensus header collapses read groups; pin that the fused
    // and standalone paths agree on the @HD/@SQ/@RG metadata.
    assert_bam_headers_equivalent_ignoring_pg(&runall_out, &standalone_out);
}

/// `Duplex → Duplex` parity vs standalone `fgumi duplex`. Runs through
/// the consensus-only delegation fast path. See
/// `parity_a_simplex_to_simplex` for the delegation contract this pins.
#[cfg(feature = "consensus")]
#[test]
fn parity_a_duplex_to_duplex() {
    let tmp = TempDir::new().unwrap();
    let fixture = grouped_duplex_fixture(tmp.path());
    let runall_out = tmp.path().join("runall.bam");
    let standalone_out = tmp.path().join("standalone.bam");
    let args = ParityArgs::for_duplex();

    let r = run_runall(Stage::Duplex, Stage::Duplex, &fixture, &runall_out, &args);
    assert!(r.status.success(), "runall: {}", String::from_utf8_lossy(&r.stderr));

    let r = run_standalone(Stage::Duplex, &fixture, &standalone_out, &args);
    assert!(r.status.success(), "standalone: {}", String::from_utf8_lossy(&r.stderr));

    // S9a-001: the duplex chain emits duplex consensus records — non-empty.
    assert_bams_record_equivalent_nonempty(&runall_out, &standalone_out);
    assert_bam_headers_equivalent_ignoring_pg(&runall_out, &standalone_out);
}

/// `Codec → Codec` parity vs standalone `fgumi codec`. Runs through the
/// consensus-only delegation fast path. Uses `grouped_codec_fixture`
/// (paired-end, single-UMI, grouped with `--strategy adjacency`).
///
/// NOTE (S9a-001): this fixture lacks the FR-overlap flag shape CODEC
/// consensus requires (R1 `PAIRED|FIRST|MATE_REVERSE` / R2 `PAIRED|LAST|REVERSE`
/// at the same position), so every pair is rejected and BOTH the fused and
/// standalone paths emit **0 consensus records**. The parity therefore holds
/// only because both sides are empty. The explicit `== 0` assertion below makes
/// that deliberate zero VISIBLE: if a future change gives the codec fixture
/// FR-overlap flags (making it non-empty), this test fails loudly and forces an
/// upgrade to a real non-empty parity check (the `zipper_*_codec` smoke tests
/// already cover the non-empty FR-overlap shape via `zipper_codec_fixture`).
#[cfg(feature = "consensus")]
#[test]
fn parity_a_codec_to_codec() {
    let tmp = TempDir::new().unwrap();
    let fixture = grouped_codec_fixture(tmp.path());
    let runall_out = tmp.path().join("runall.bam");
    let standalone_out = tmp.path().join("standalone.bam");
    let args = ParityArgs::for_simplex_like();

    let r = run_runall(Stage::Codec, Stage::Codec, &fixture, &runall_out, &args);
    assert!(r.status.success(), "runall: {}", String::from_utf8_lossy(&r.stderr));

    let r = run_standalone(Stage::Codec, &fixture, &standalone_out, &args);
    assert!(r.status.success(), "standalone: {}", String::from_utf8_lossy(&r.stderr));

    assert_bams_record_equivalent(&runall_out, &standalone_out);
    // Header parity matters most in this empty-stream case: a wrong fused
    // @HD/@SQ/@RG would otherwise stay green behind the zero-record assertion.
    assert_bam_headers_equivalent_ignoring_pg(&runall_out, &standalone_out);
    assert_eq!(
        read_bam_records(&runall_out).len(),
        0,
        "codec fixture lacks the FR-overlap shape → 0 consensus by design (S9a-001)"
    );
}

// ────────────────────────── Class B — multi-stage parity ──────────────────────────
//
// Each Class B test runs `fgumi runall --start-from S --stop-after T`
// against the equivalent staged chain. All seven Class B combos are
// accepted by the validator and run by default.

/// `Sort → Group` parity vs staged `fgumi sort | fgumi group`.
/// Runall delegates to `Sort::execute` then `GroupReadsByUmi::execute`
/// (via a tempfile), so the output is record- and header-equivalent
/// (ignoring `@PG` command-line provenance) by construction.
#[test]
fn parity_b_sort_to_group() {
    let tmp = TempDir::new().unwrap();
    let fixture = unsorted_simplex_fixture(tmp.path());
    let runall_out = tmp.path().join("runall.bam");
    let staged_out = tmp.path().join("staged.bam");
    let args = ParityArgs::for_simplex_like();

    let r = run_runall(Stage::Sort, Stage::Group, &fixture, &runall_out, &args);
    assert!(r.status.success(), "runall: {}", String::from_utf8_lossy(&r.stderr));

    let r =
        run_staged_chain(&[Stage::Sort, Stage::Group], &fixture, &staged_out, tmp.path(), &args);
    assert!(r.status.success(), "staged: {}", String::from_utf8_lossy(&r.stderr));

    assert_bams_record_equivalent_nonempty(&runall_out, &staged_out);
    // S9a-005: the stack changes header propagation/update paths, so also pin
    // @HD/@SQ/@RG parity — a fused-chain header regression (e.g. SO:unsorted
    // after sort, a dropped @SQ, uncollapsed @RG) is invisible to the
    // record-only comparison above.
    assert_bam_headers_equivalent_ignoring_pg(&runall_out, &staged_out);
}

/// `Sort → Simplex` parity vs staged `fgumi sort | fgumi group | fgumi simplex`.
#[cfg(feature = "consensus")]
#[test]
fn parity_b_sort_to_simplex() {
    let tmp = TempDir::new().unwrap();
    let fixture = unsorted_simplex_fixture(tmp.path());
    let runall_out = tmp.path().join("runall.bam");
    let staged_out = tmp.path().join("staged.bam");
    let args = ParityArgs::for_simplex_like();

    let r = run_runall(Stage::Sort, Stage::Simplex, &fixture, &runall_out, &args);
    assert!(r.status.success(), "runall: {}", String::from_utf8_lossy(&r.stderr));

    let r = run_staged_chain(
        &[Stage::Sort, Stage::Group, Stage::Simplex],
        &fixture,
        &staged_out,
        tmp.path(),
        &args,
    );
    assert!(r.status.success(), "staged: {}", String::from_utf8_lossy(&r.stderr));

    assert_bams_record_equivalent_nonempty(&runall_out, &staged_out);
    // S9a-005: the stack changes header propagation/update paths, so also pin
    // @HD/@SQ/@RG parity — a fused-chain header regression (e.g. SO:unsorted
    // after sort, a dropped @SQ, uncollapsed @RG) is invisible to the
    // record-only comparison above.
    assert_bam_headers_equivalent_ignoring_pg(&runall_out, &staged_out);
}

/// `Sort → Duplex` parity vs staged `fgumi sort | fgumi group | fgumi duplex`.
/// `create_paired_umi_family` sets MC tags programmatically so the
/// fixture flows through `GroupByPosition` (which requires MC on
/// paired-end inputs).
#[cfg(feature = "consensus")]
#[test]
fn parity_b_sort_to_duplex() {
    let tmp = TempDir::new().unwrap();
    let fixture = unsorted_duplex_fixture(tmp.path());
    let runall_out = tmp.path().join("runall.bam");
    let staged_out = tmp.path().join("staged.bam");
    let args = ParityArgs::for_duplex();

    let r = run_runall(Stage::Sort, Stage::Duplex, &fixture, &runall_out, &args);
    assert!(r.status.success(), "runall: {}", String::from_utf8_lossy(&r.stderr));

    let r = run_staged_chain(
        &[Stage::Sort, Stage::Group, Stage::Duplex],
        &fixture,
        &staged_out,
        tmp.path(),
        &args,
    );
    assert!(r.status.success(), "staged: {}", String::from_utf8_lossy(&r.stderr));

    assert_bams_record_equivalent_nonempty(&runall_out, &staged_out);
    // S9a-005: the stack changes header propagation/update paths, so also pin
    // @HD/@SQ/@RG parity — a fused-chain header regression (e.g. SO:unsorted
    // after sort, a dropped @SQ, uncollapsed @RG) is invisible to the
    // record-only comparison above.
    assert_bam_headers_equivalent_ignoring_pg(&runall_out, &staged_out);
}

/// `Sort → Codec` parity vs staged `fgumi sort | fgumi group | fgumi codec`.
/// Uses the codec fixture family (paired-end, single-strand UMI,
/// `--strategy adjacency`) so the consensus stage actually exercises
/// CODEC's duplex logic.
#[cfg(feature = "consensus")]
#[test]
fn parity_b_sort_to_codec() {
    let tmp = TempDir::new().unwrap();
    let fixture = unsorted_codec_fixture(tmp.path());
    let runall_out = tmp.path().join("runall.bam");
    let staged_out = tmp.path().join("staged.bam");
    let args = ParityArgs::for_simplex_like();

    let r = run_runall(Stage::Sort, Stage::Codec, &fixture, &runall_out, &args);
    assert!(r.status.success(), "runall: {}", String::from_utf8_lossy(&r.stderr));

    let r = run_staged_chain(
        &[Stage::Sort, Stage::Group, Stage::Codec],
        &fixture,
        &staged_out,
        tmp.path(),
        &args,
    );
    assert!(r.status.success(), "staged: {}", String::from_utf8_lossy(&r.stderr));

    assert_bams_record_equivalent(&runall_out, &staged_out);
    // S9a-005: header parity holds even for the zero-record codec branch (both
    // sides still emit @HD/@SQ/@RG) — assert it before the zero-count pin below.
    assert_bam_headers_equivalent_ignoring_pg(&runall_out, &staged_out);
    // S9a-001: the codec fixture lacks the FR-overlap shape → 0 consensus by
    // design; pin the deliberate zero so a future non-empty fixture forces a
    // real parity check.
    assert_eq!(read_bam_records(&runall_out).len(), 0, "codec fixture → 0 consensus by design");
}

/// `Group → Simplex` parity vs staged `fgumi group | fgumi simplex`.
#[cfg(feature = "consensus")]
#[test]
fn parity_b_group_to_simplex() {
    let tmp = TempDir::new().unwrap();
    let fixture = sorted_simplex_fixture(tmp.path());
    let runall_out = tmp.path().join("runall.bam");
    let staged_out = tmp.path().join("staged.bam");
    let args = ParityArgs::for_simplex_like();

    let r = run_runall(Stage::Group, Stage::Simplex, &fixture, &runall_out, &args);
    assert!(r.status.success(), "runall: {}", String::from_utf8_lossy(&r.stderr));

    let r =
        run_staged_chain(&[Stage::Group, Stage::Simplex], &fixture, &staged_out, tmp.path(), &args);
    assert!(r.status.success(), "staged: {}", String::from_utf8_lossy(&r.stderr));

    assert_bams_record_equivalent_nonempty(&runall_out, &staged_out);
    // S9a-005: the stack changes header propagation/update paths, so also pin
    // @HD/@SQ/@RG parity — a fused-chain header regression (e.g. SO:unsorted
    // after sort, a dropped @SQ, uncollapsed @RG) is invisible to the
    // record-only comparison above.
    assert_bam_headers_equivalent_ignoring_pg(&runall_out, &staged_out);
}

/// `Group → Duplex` parity vs staged `fgumi group | fgumi duplex`.
/// See [`parity_b_sort_to_duplex`] for the MC-tag fixture note.
#[cfg(feature = "consensus")]
#[test]
fn parity_b_group_to_duplex() {
    let tmp = TempDir::new().unwrap();
    let fixture = sorted_duplex_fixture(tmp.path());
    let runall_out = tmp.path().join("runall.bam");
    let staged_out = tmp.path().join("staged.bam");
    let args = ParityArgs::for_duplex();

    let r = run_runall(Stage::Group, Stage::Duplex, &fixture, &runall_out, &args);
    assert!(r.status.success(), "runall: {}", String::from_utf8_lossy(&r.stderr));

    let r =
        run_staged_chain(&[Stage::Group, Stage::Duplex], &fixture, &staged_out, tmp.path(), &args);
    assert!(r.status.success(), "staged: {}", String::from_utf8_lossy(&r.stderr));

    assert_bams_record_equivalent_nonempty(&runall_out, &staged_out);
    // S9a-005: the stack changes header propagation/update paths, so also pin
    // @HD/@SQ/@RG parity — a fused-chain header regression (e.g. SO:unsorted
    // after sort, a dropped @SQ, uncollapsed @RG) is invisible to the
    // record-only comparison above.
    assert_bam_headers_equivalent_ignoring_pg(&runall_out, &staged_out);
}

/// NEW-002: under `--group::allow-unmapped`, a both-unmapped read pair survives
/// `group` but must be dropped by the duplex `record_filter` on BOTH the fused
/// `runall group→duplex` path and the staged `fgumi group | fgumi duplex` path.
///
/// Assert that a `group --allow-unmapped` output BAM retains the fixture's sole
/// both-unmapped template (`fam_unmapped`, UMI `TTAA-GGCC`; see
/// `sorted_duplex_with_unmapped_fixture`) by IDENTITY, not just by count. A
/// count-only check would still pass if the group stage retained two of some
/// *other* record while dropping a `fam_unmapped` mate.
fn assert_group_retains_both_unmapped_pair(group_bam: &Path) {
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::data::field::Value;
    let mi_tag = Tag::from(fgumi_lib::sam::SamTag::MI);
    let rx_tag = Tag::from(fgumi_lib::sam::SamTag::RX);
    let group_records = read_bam_records(group_bam);
    let retained: Vec<_> = group_records
        .iter()
        .filter(|rec| rec.flags().is_unmapped() && rec.data().get(&mi_tag).is_some())
        .collect();
    assert_eq!(
        retained.len(),
        2,
        "group --allow-unmapped output must retain BOTH mates of the one both-unmapped \
         template (= 2 MI-tagged unmapped records); a count other than 2 means a mate was \
         dropped or duplicated. Found {} such records in {}",
        retained.len(),
        group_bam.display()
    );
    // Both retained records must be the `fam_unmapped` pair, carrying its UMI.
    for rec in &retained {
        let qname: &[u8] = rec.name().expect("retained record has a read name").as_ref();
        assert_eq!(
            qname, b"fam_unmapped",
            "retained MI-tagged unmapped record has an unexpected qname (not fam_unmapped)"
        );
        let Some(Value::String(rx)) = rec.data().get(&rx_tag) else {
            panic!("retained unmapped record is missing its RX (UMI) tag");
        };
        let rx_bytes: &[u8] = rx.as_ref();
        assert_eq!(
            rx_bytes, b"TTAA-GGCC",
            "retained unmapped record carries the wrong UMI (not the fam_unmapped UMI)"
        );
    }
    // Both retained mates must share ONE molecule id: presence alone (the filter
    // above) would still pass if a regression split `fam_unmapped` into two
    // different MI values, changing the grouped-template identity.
    let mi_values: std::collections::BTreeSet<Vec<u8>> = retained
        .iter()
        .map(|rec| {
            let Some(Value::String(mi)) = rec.data().get(&mi_tag) else {
                panic!("retained unmapped record is missing its MI tag");
            };
            let mi_bytes: &[u8] = mi.as_ref();
            mi_bytes.to_vec()
        })
        .collect();
    assert_eq!(
        mi_values.len(),
        1,
        "the retained both-unmapped pair must share one MI value; got {mi_values:?}"
    );
    // ...and they are exactly one R1 + one R2 (no mate dropped or duplicated).
    let first = retained.iter().filter(|r| r.flags().is_first_segment()).count();
    let last = retained.iter().filter(|r| r.flags().is_last_segment()).count();
    assert_eq!(
        (first, last),
        (1, 1),
        "the retained both-unmapped pair must be exactly one R1 + one R2, got \
         ({first} first, {last} last) segments"
    );
}

/// Before the fix the fused `templates_to_mi_step` bridge did not apply the
/// duplex filter, so the both-unmapped pair leaked into the consensus caller on
/// the fused path only — diverging from the staged path. This test fails on the
/// unfixed bridge and passes once the bridge re-applies `duplex_record_filter`.
#[test]
fn new_002_group_to_duplex_allow_unmapped_parity() {
    let tmp = TempDir::new().unwrap();
    let fixture = sorted_duplex_with_unmapped_fixture(tmp.path());
    let runall_out = tmp.path().join("runall.bam");
    let group_bam = tmp.path().join("staged_group.bam");
    let staged_out = tmp.path().join("staged.bam");
    let args = ParityArgs::for_duplex();

    // Fused: runall group→duplex with allow-unmapped on the group stage.
    let fused_args: Vec<OsString> = vec![
        "runall".into(),
        "--start-from".into(),
        Stage::Group.runall_stage_value().into(),
        "--stop-after".into(),
        Stage::Duplex.runall_stage_value().into(),
        "--consensus".into(),
        "duplex".into(),
        "--input".into(),
        fixture.as_os_str().to_owned(),
        "--output".into(),
        runall_out.as_os_str().to_owned(),
        "--group::strategy".into(),
        args.strategy.into(),
        "--group::edits".into(),
        args.edits.to_string().into(),
        "--group::allow-unmapped".into(),
        "true".into(),
        "--duplex::min-reads".into(),
        args.min_reads.into(),
        "--threads".into(),
        args.threads.to_string().into(),
    ];
    let r = fgumi(&fused_args);
    assert!(r.status.success(), "fused runall: {}", String::from_utf8_lossy(&r.stderr));

    // Staged: standalone group --allow-unmapped, then standalone duplex.
    let group_args: Vec<OsString> = vec![
        "group".into(),
        "--input".into(),
        fixture.as_os_str().to_owned(),
        "--output".into(),
        group_bam.as_os_str().to_owned(),
        "--strategy".into(),
        args.strategy.into(),
        "--edits".into(),
        args.edits.to_string().into(),
        "--allow-unmapped".into(),
        "--threads".into(),
        args.threads.to_string().into(),
    ];
    let r = fgumi(&group_args);
    assert!(r.status.success(), "staged group: {}", String::from_utf8_lossy(&r.stderr));

    // Independent intermediate check: the group stage with --allow-unmapped must
    // retain the both-unmapped template (MI-tagged, by identity) BEFORE duplex
    // runs. This pins the bridge/filter contract directly — without it, a
    // regression that drops the unmapped template at the group stage could still
    // pass the final fused-vs-staged comparison if both paths dropped it
    // identically.
    assert_group_retains_both_unmapped_pair(&group_bam);

    let duplex_args: Vec<OsString> = vec![
        "duplex".into(),
        "--input".into(),
        group_bam.as_os_str().to_owned(),
        "--output".into(),
        staged_out.as_os_str().to_owned(),
        "--min-reads".into(),
        args.min_reads.into(),
        "--threads".into(),
        args.threads.to_string().into(),
    ];
    let r = fgumi(&duplex_args);
    assert!(r.status.success(), "staged duplex: {}", String::from_utf8_lossy(&r.stderr));

    assert_bams_record_equivalent_nonempty(&runall_out, &staged_out);
    // S9a-005: the stack changes header propagation/update paths, so also pin
    // @HD/@SQ/@RG parity — a fused-chain header regression (e.g. SO:unsorted
    // after sort, a dropped @SQ, uncollapsed @RG) is invisible to the
    // record-only comparison above.
    assert_bam_headers_equivalent_ignoring_pg(&runall_out, &staged_out);

    // Negative oracle: path-equivalence alone would still pass if BOTH paths
    // wrongly emitted output derived from the both-unmapped family. The contract
    // is that the duplex filter DROPS `fam_unmapped` before consensus, so no
    // output consensus read may carry its UMI. The four mapped families use
    // distinct UMIs (AAAA-CCCC / GGGG-TTTT / CCAA-TTGG / ATCG-GCTA), so the
    // `fam_unmapped` UMI (TTAA-GGCC) appearing in any output `RX` would mean it
    // leaked into consensus. Assert it is absent from both paths' output.
    let rx_tag = noodles::sam::alignment::record::data::field::Tag::from(SamTag::RX);
    for (label, path) in [("runall", &runall_out), ("staged", &staged_out)] {
        let mut reader = noodles::bam::io::Reader::new(std::fs::File::open(path).unwrap());
        let header = reader.read_header().unwrap();
        for result in reader.record_bufs(&header) {
            let rec = result.expect("read output record");
            if let Some(noodles::sam::alignment::record_buf::data::field::Value::String(rx)) =
                rec.data().get(&rx_tag)
            {
                let rx_bytes: &[u8] = rx.as_ref();
                assert_ne!(
                    rx_bytes, b"TTAA-GGCC",
                    "{label} output has a consensus read with the both-unmapped family's UMI — \
                     fam_unmapped was not dropped by the duplex filter",
                );
            }
        }
    }
}

/// `Group → Codec` parity vs staged `fgumi group | fgumi codec`.
/// Uses the codec fixture family; see [`parity_b_sort_to_codec`] for
/// the rationale.
#[cfg(feature = "consensus")]
#[test]
fn parity_b_group_to_codec() {
    let tmp = TempDir::new().unwrap();
    let fixture = sorted_codec_fixture(tmp.path());
    let runall_out = tmp.path().join("runall.bam");
    let staged_out = tmp.path().join("staged.bam");
    let args = ParityArgs::for_simplex_like();

    let r = run_runall(Stage::Group, Stage::Codec, &fixture, &runall_out, &args);
    assert!(r.status.success(), "runall: {}", String::from_utf8_lossy(&r.stderr));

    let r =
        run_staged_chain(&[Stage::Group, Stage::Codec], &fixture, &staged_out, tmp.path(), &args);
    assert!(r.status.success(), "staged: {}", String::from_utf8_lossy(&r.stderr));

    assert_bams_record_equivalent(&runall_out, &staged_out);
    // S9a-005: header parity holds even for the zero-record codec branch.
    assert_bam_headers_equivalent_ignoring_pg(&runall_out, &staged_out);
    // S9a-001: codec fixture lacks FR-overlap → 0 consensus by design.
    assert_eq!(read_bam_records(&runall_out).len(), 0, "codec fixture → 0 consensus by design");
}

// ──────────────── Fused consensus → filter parity (issue #330 option a) ────────────────
//
// `runall --start-from group --consensus <mode> --stop-after filter` fuses
// the consensus caller and the filter stage into one in-memory pipeline (no
// intermediate consensus BAM). Each test compares that against the staged
// equivalent: produce the consensus BAM with `runall ... --stop-after
// consensus`, then run standalone `fgumi filter` on it. Both sides apply the
// same default filter parameters, so the filtered records must match.

/// `Group → Simplex → Filter` (fused) vs staged consensus BAM + `fgumi filter`.
#[cfg(feature = "consensus")]
#[test]
fn parity_group_to_simplex_to_filter() {
    let tmp = TempDir::new().unwrap();
    let fixture = sorted_simplex_fixture(tmp.path());
    let fused_out = tmp.path().join("fused.bam");
    let consensus_bam = tmp.path().join("consensus.bam");
    let staged_out = tmp.path().join("staged.bam");
    let args = ParityArgs::for_simplex_like();

    // Fused: group → simplex → filter in one runall invocation.
    let r =
        run_runall_consensus_to_filter(Stage::Group, Stage::Simplex, &fixture, &fused_out, &args);
    assert!(r.status.success(), "fused: {}", String::from_utf8_lossy(&r.stderr));

    // Staged: produce the consensus BAM, then filter it standalone.
    let r = run_runall(Stage::Group, Stage::Simplex, &fixture, &consensus_bam, &args);
    assert!(r.status.success(), "staged consensus: {}", String::from_utf8_lossy(&r.stderr));
    let r = run_standalone_filter(&consensus_bam, &staged_out, &args);
    assert!(r.status.success(), "staged filter: {}", String::from_utf8_lossy(&r.stderr));

    // S9a-001 / S5b2-004: simplex consensus at --min-reads 1 survives the
    // filter, so the fused consensus→filter stream is non-empty.
    assert_bams_record_equivalent_nonempty(&fused_out, &staged_out);
    // S9a-005: the fused consensus→filter path also updates the header; pin
    // @HD/@SQ/@RG parity against the staged consensus + standalone-filter chain.
    assert_bam_headers_equivalent_ignoring_pg(&fused_out, &staged_out);
}

/// `Group → Duplex → Filter` (fused) vs staged consensus BAM + `fgumi filter`.
#[cfg(feature = "consensus")]
#[test]
fn parity_group_to_duplex_to_filter() {
    let tmp = TempDir::new().unwrap();
    let fixture = sorted_duplex_fixture(tmp.path());
    let fused_out = tmp.path().join("fused.bam");
    let consensus_bam = tmp.path().join("consensus.bam");
    let staged_out = tmp.path().join("staged.bam");
    let args = ParityArgs::for_duplex();

    let r =
        run_runall_consensus_to_filter(Stage::Group, Stage::Duplex, &fixture, &fused_out, &args);
    assert!(r.status.success(), "fused: {}", String::from_utf8_lossy(&r.stderr));

    let r = run_runall(Stage::Group, Stage::Duplex, &fixture, &consensus_bam, &args);
    assert!(r.status.success(), "staged consensus: {}", String::from_utf8_lossy(&r.stderr));
    let r = run_standalone_filter(&consensus_bam, &staged_out, &args);
    assert!(r.status.success(), "staged filter: {}", String::from_utf8_lossy(&r.stderr));

    // S9a-001 / S5b2-004: this pins the fused consensus→filter topology (the
    // duplex strand-strip bridge) — fused MUST equal staged. The duplex fixture
    // is single-strand (all reads of each family carry the same paired UMI →
    // AB depth 4, BA depth 0), so `fgumi filter --min-reads 1` drops every
    // duplex consensus read for failing the per-strand depth floor (BA = 0 < 1).
    // Both paths therefore yield 0 — assert that deliberate zero explicitly
    // (rather than letting the equivalence pass vacuously). `parity_b_group_to_
    // duplex` already pins the NON-empty duplex consensus stream before filtering.
    assert_bams_record_equivalent(&fused_out, &staged_out);
    // S9a-005: header parity holds even for the zero-record duplex→filter branch.
    assert_bam_headers_equivalent_ignoring_pg(&fused_out, &staged_out);
    assert_eq!(
        read_bam_records(&fused_out).len(),
        0,
        "single-strand duplex fixture → 0 reads survive `filter --min-reads 1` (BA depth 0)"
    );
}

/// `Group → Codec → Filter` (fused) vs staged consensus BAM + `fgumi filter`.
#[cfg(feature = "consensus")]
#[test]
fn parity_group_to_codec_to_filter() {
    let tmp = TempDir::new().unwrap();
    let fixture = sorted_codec_fixture(tmp.path());
    let fused_out = tmp.path().join("fused.bam");
    let consensus_bam = tmp.path().join("consensus.bam");
    let staged_out = tmp.path().join("staged.bam");
    let args = ParityArgs::for_simplex_like();

    let r = run_runall_consensus_to_filter(Stage::Group, Stage::Codec, &fixture, &fused_out, &args);
    assert!(r.status.success(), "fused: {}", String::from_utf8_lossy(&r.stderr));

    let r = run_runall(Stage::Group, Stage::Codec, &fixture, &consensus_bam, &args);
    assert!(r.status.success(), "staged consensus: {}", String::from_utf8_lossy(&r.stderr));
    let r = run_standalone_filter(&consensus_bam, &staged_out, &args);
    assert!(r.status.success(), "staged filter: {}", String::from_utf8_lossy(&r.stderr));

    assert_bams_record_equivalent(&fused_out, &staged_out);
    // S9a-005: header parity holds even for the zero-record codec→filter branch.
    assert_bam_headers_equivalent_ignoring_pg(&fused_out, &staged_out);
    // S9a-001: codec fixture → 0 consensus → 0 after filter, by design.
    assert_eq!(read_bam_records(&fused_out).len(), 0, "codec fixture → 0 consensus by design");
}

// ───────────────────── Run-to-run / thread-count determinism ─────────────────────
//
// The fused pipeline must emit a deterministic record stream regardless of
// thread count: the same input produces the same records in the same order
// whether run at --threads 1 or --threads 8, and run-to-run. This guards the
// `ByItemOrdinal` output-ordering contract against a regression that would let
// a step emit in thread-completion order (which the benchmark fleet would see
// as run-vs-run / runall-vs-standalone divergence). Headers (`@PG`) are
// intentionally excluded — they record the literal command line.

/// Build a moderate sorted duplex fixture: many paired-UMI duplex molecules
/// (top + bottom strand each) spread across distinct template-coordinate
/// positions, so there are enough distinct templates to span multiple
/// consensus batches and parallel workers (where thread-completion order
/// would diverge from input order if ordering were broken).
fn deterministic_duplex_fixture(dir: &Path, positions: usize) -> PathBuf {
    let unsorted = dir.join("det_unsorted.bam");
    let header = create_minimal_header("chr1", positions * 200 + 1000);
    let mut writer = bam::io::Writer::new(fs::File::create(&unsorted).expect("create det fixture"));
    writer.write_header(&header).expect("write det header");
    let mut all = Vec::new();
    for p in 0..positions {
        let r1_pos = 99 + p * 200;
        // Top and bottom strand of one duplex molecule: the `paired` strategy
        // canonicalizes the swapped UMI ("A-B" / "B-A") to one MI (/A + /B).
        for r in create_paired_umi_family_at(
            "AAAA-CCCC",
            3,
            &format!("f{p}t"),
            "ACGTACGT",
            "TTTTAAAA",
            30,
            r1_pos,
        ) {
            all.push(r);
        }
        for r in create_paired_umi_family_at(
            "CCCC-AAAA",
            3,
            &format!("f{p}b"),
            "ACGTACGT",
            "TTTTAAAA",
            30,
            r1_pos,
        ) {
            all.push(r);
        }
    }
    all.reverse();
    for raw in &all {
        writer.write_alignment_record(&header, &to_record_buf(raw)).expect("write det record");
    }
    writer.try_finish().expect("finish det fixture");

    let sorted = dir.join("det_sorted.bam");
    let args = ParityArgs::for_duplex();
    let out = run_standalone(Stage::Sort, &unsorted, &sorted, &args);
    assert!(out.status.success(), "det sort: {}", String::from_utf8_lossy(&out.stderr));
    sorted
}

/// The fused `group → duplex` pipeline must produce a byte-identical record
/// stream run-to-run at `--threads 8`, and the same record stream at
/// `--threads 1` as at `--threads 8` (thread-count-invariant output).
#[cfg(feature = "consensus")]
#[test]
fn runall_duplex_record_stream_is_deterministic() {
    let tmp = TempDir::new().unwrap();
    let fixture = deterministic_duplex_fixture(tmp.path(), 200);

    let run = |out: &Path, threads: u32| {
        let mut args = ParityArgs::for_duplex();
        args.threads = threads;
        let r = run_runall(Stage::Group, Stage::Duplex, &fixture, out, &args);
        assert!(
            r.status.success(),
            "runall duplex (threads={threads}): {}",
            String::from_utf8_lossy(&r.stderr)
        );
    };

    let out_a = tmp.path().join("run_t8_a.bam");
    let out_b = tmp.path().join("run_t8_b.bam");
    let out_t1 = tmp.path().join("run_t1.bam");
    run(&out_a, 8);
    run(&out_b, 8);
    run(&out_t1, 1);

    // Run-to-run determinism at high thread count.
    assert_bams_record_equivalent(&out_a, &out_b);
    // Thread-count invariance: t8 output == t1 output.
    assert_bams_record_equivalent(&out_a, &out_t1);
}

/// Build a many-position single-end simplex fixture spanning enough distinct
/// template-coordinate positions that the sort/group/consensus stages run
/// across multiple parallel workers — the regime where a broken
/// `ByItemOrdinal` contract would let output diverge by thread-completion
/// order. Returns a SORTED BAM (so `group` / `simplex` can consume it).
fn deterministic_simplex_fixture(dir: &Path, positions: usize) -> PathBuf {
    let unsorted = dir.join("det_simplex_unsorted.bam");
    let header = create_minimal_header("chr1", positions * 200 + 1000);
    let mut writer =
        bam::io::Writer::new(fs::File::create(&unsorted).expect("create det simplex fixture"));
    writer.write_header(&header).expect("write det simplex header");
    let mut all = Vec::new();
    // A small set of valid-DNA UMIs (non-ACGT chars would be rejected by the
    // UMI grouping stage); cycle through them so distinct families form across
    // positions. Each family is mapped to its own template-coordinate position
    // (100 bp apart) via `create_umi_family_at`, so the families genuinely span
    // many coordinates — exercising the sort/group/consensus path across
    // multiple parallel workers rather than collapsing onto one position.
    let umis = ["ACGTACGT", "TGCATGCA", "CCAATTGG", "GGTTAACC", "AACCGGTT", "TTGGCCAA"];
    for p in 0..positions {
        let umi = umis[p % umis.len()];
        let pos = i32::try_from(100 + p * 100).expect("family position fits in i32");
        for r in create_umi_family_at(pos, umi, 3, &format!("ds{p}_{umi}"), "ACGTACGTACGT", 30) {
            all.push(r);
        }
    }
    all.reverse();
    for raw in &all {
        writer
            .write_alignment_record(&header, &to_record_buf(raw))
            .expect("write det simplex record");
    }
    writer.try_finish().expect("finish det simplex fixture");

    let sorted = dir.join("det_simplex_sorted.bam");
    let args = ParityArgs::for_simplex_like();
    let out = run_standalone(Stage::Sort, &unsorted, &sorted, &args);
    assert!(out.status.success(), "det simplex sort: {}", String::from_utf8_lossy(&out.stderr));
    sorted
}

/// S9a-008: extend the determinism guard beyond `group → duplex` to the
/// streaming `sort` and the `group → simplex` parallel paths — the other
/// stages whose output most plausibly emits in thread-completion order if the
/// `ByItemOrdinal` contract regressed. Each asserts run-to-run equivalence at
/// `--threads 8` AND thread-count invariance (`t8 == t1`).
#[test]
fn runall_sort_record_stream_is_deterministic() {
    let tmp = TempDir::new().unwrap();
    // The sort engine itself is the parallel stage under test — feed it the
    // UNSORTED, many-position duplex fixture and sort it via runall.
    let unsorted = tmp.path().join("sort_det_unsorted.bam");
    {
        let positions = 200usize;
        let header = create_minimal_header("chr1", positions * 200 + 1000);
        let mut writer =
            bam::io::Writer::new(fs::File::create(&unsorted).expect("create sort det fixture"));
        writer.write_header(&header).expect("write sort det header");
        let mut all = Vec::new();
        for p in 0..positions {
            for r in create_paired_umi_family_at(
                "AAAA-CCCC",
                3,
                &format!("sd{p}"),
                "ACGTACGT",
                "TTTTAAAA",
                30,
                99 + p * 200,
            ) {
                all.push(r);
            }
        }
        all.reverse();
        for raw in &all {
            writer
                .write_alignment_record(&header, &to_record_buf(raw))
                .expect("write sort det record");
        }
        writer.try_finish().expect("finish sort det fixture");
    }

    let run = |out: &Path, threads: u32| {
        let mut args = ParityArgs::for_simplex_like();
        args.threads = threads;
        let r = run_runall(Stage::Sort, Stage::Sort, &unsorted, out, &args);
        assert!(
            r.status.success(),
            "runall sort (threads={threads}): {}",
            String::from_utf8_lossy(&r.stderr)
        );
    };
    let out_a = tmp.path().join("sort_t8_a.bam");
    let out_b = tmp.path().join("sort_t8_b.bam");
    let out_t1 = tmp.path().join("sort_t1.bam");
    run(&out_a, 8);
    run(&out_b, 8);
    run(&out_t1, 1);
    assert_bams_record_equivalent_nonempty(&out_a, &out_b);
    assert_bams_record_equivalent_nonempty(&out_a, &out_t1);
}

/// S9a-008: `group → simplex` determinism (run-to-run + thread-count invariant).
#[test]
fn runall_simplex_record_stream_is_deterministic() {
    let tmp = TempDir::new().unwrap();
    let fixture = deterministic_simplex_fixture(tmp.path(), 200);

    let run = |out: &Path, threads: u32| {
        let mut args = ParityArgs::for_simplex_like();
        args.threads = threads;
        let r = run_runall(Stage::Group, Stage::Simplex, &fixture, out, &args);
        assert!(
            r.status.success(),
            "runall simplex (threads={threads}): {}",
            String::from_utf8_lossy(&r.stderr)
        );
    };
    let out_a = tmp.path().join("simplex_t8_a.bam");
    let out_b = tmp.path().join("simplex_t8_b.bam");
    let out_t1 = tmp.path().join("simplex_t1.bam");
    run(&out_a, 8);
    run(&out_b, 8);
    run(&out_t1, 1);
    assert_bams_record_equivalent_nonempty(&out_a, &out_b);
    assert_bams_record_equivalent_nonempty(&out_a, &out_t1);
}

// ───────────────────────── stdin + FASTQ-source coverage ─────────────────────────
//
// X3-001/X3-002: runall is the one command that FUSES the stdin-aware sources,
// yet it had no stdin test and no FASTQ-entry-point (`--start-from extract`)
// end-to-end test. These close both gaps: the "reads the source exactly once"
// invariant (which fusion is most likely to break by re-opening / double-pulling
// a non-seekable `-`) and the FASTQ→consensus fusion path.

/// X3-001: `fgumi runall --start-from group --stop-after simplex --input -`
/// reading a BAM piped to stdin must produce the SAME (non-empty) record stream
/// as the identical runall reading the BAM from a file path. A regression that
/// made the fused chain pull stdin twice (or seek it) would read an empty stream
/// the second time and diverge; every BAM-source-from-file test would still pass.
#[test]
fn runall_reads_bam_from_stdin_once() {
    use std::process::{Command, Stdio};

    let tmp = TempDir::new().unwrap();
    let fixture = sorted_simplex_fixture(tmp.path());
    let from_file = tmp.path().join("from_file.bam");
    let from_stdin = tmp.path().join("from_stdin.bam");
    let args = ParityArgs::for_simplex_like();

    // Baseline: file input.
    let r = run_runall(Stage::Group, Stage::Simplex, &fixture, &from_file, &args);
    assert!(r.status.success(), "runall (file input): {}", String::from_utf8_lossy(&r.stderr));

    // Piped: `cat fixture.bam | fgumi runall ... --input -`.
    let mut cat = Command::new("cat")
        .arg(&fixture)
        .stdout(Stdio::piped())
        .spawn()
        .expect("spawn cat for stdin pipe");
    let cat_stdout = cat.stdout.take().expect("cat stdout");
    let stdin_out = Command::new(fgumi_binary())
        .args([
            "runall",
            "--start-from",
            "group",
            "--stop-after",
            "consensus",
            "--consensus",
            "simplex",
            "--input",
            "-",
            "--output",
            from_stdin.to_str().unwrap(),
            "--group::strategy",
            args.strategy,
            "--group::edits",
            "1",
            "--simplex::min-reads",
            args.min_reads,
            "--threads",
            "1",
        ])
        .stdin(cat_stdout)
        .output()
        .expect("run runall with stdin");
    // Reap the `cat` producer so it isn't left unwaited.
    let _ = cat.wait().expect("wait on cat");
    assert!(
        stdin_out.status.success(),
        "runall (stdin input): {}",
        String::from_utf8_lossy(&stdin_out.stderr)
    );

    // Same record stream both ways, and non-empty (so the stdin source was read
    // exactly once — a double-pull would yield 0 records on the second read).
    assert_bams_record_equivalent_nonempty(&from_file, &from_stdin);
}

/// CG-5 (audit E2): every runall failure test rejects *before* the pipeline runs
/// (invalid stage pair, missing `--unmapped`/`--ref`/`--umis`, bad flag combo).
/// None checks that the fused BAM-in pipeline fails *cleanly mid-run* when a
/// downstream stage hits malformed input. Feed `runall --start-from group
/// --stop-after consensus --consensus simplex` a BAM truncated mid-stream and
/// assert a clean non-zero exit — no panic / `unreachable`, and (via a
/// wall-clock watchdog) no hang.
#[test]
fn runall_group_simplex_fails_cleanly_on_truncated_bam() {
    use std::io::Read as _;
    use std::process::{Command, Stdio};
    use std::time::{Duration, Instant};

    let tmp = TempDir::new().unwrap();
    let fixture = sorted_simplex_fixture(tmp.path());

    // Truncate the valid grouped-input BAM: keep the header + some BGZF blocks,
    // then cut, so the reader errors partway through a real record stream rather
    // than failing an up-front header parse.
    let bytes = fs::read(&fixture).expect("read fixture");
    assert!(bytes.len() > 200, "fixture unexpectedly tiny ({} bytes)", bytes.len());
    let truncated = tmp.path().join("truncated.bam");
    let keep = bytes.len() * 3 / 5; // 60% — past the header, into the block stream
    fs::write(&truncated, &bytes[..keep]).expect("write truncated bam");

    let out_path = tmp.path().join("out.bam");
    let args = ParityArgs::for_simplex_like();
    let mut child = Command::new(fgumi_binary())
        .args([
            "runall",
            "--start-from",
            "group",
            "--stop-after",
            "consensus",
            "--consensus",
            "simplex",
            "--input",
            truncated.to_str().unwrap(),
            "--output",
            out_path.to_str().unwrap(),
            "--group::strategy",
            args.strategy,
            "--group::edits",
            "1",
            "--simplex::min-reads",
            args.min_reads,
            "--threads",
            "2",
        ])
        .stdout(Stdio::null())
        .stderr(Stdio::piped())
        .spawn()
        .expect("spawn fgumi runall");

    // Drain stderr on a worker so a hung child (which never closes the pipe)
    // cannot block the watchdog loop below.
    let mut stderr_pipe = child.stderr.take().expect("child stderr");
    let reader = std::thread::spawn(move || {
        let mut s = String::new();
        let _ = stderr_pipe.read_to_string(&mut s);
        s
    });

    // Wall-clock watchdog: a mid-run corruption must fail fast, never hang.
    let deadline = Instant::now() + Duration::from_secs(60);
    let status = loop {
        if let Some(status) = child.try_wait().expect("try_wait") {
            break status;
        }
        if Instant::now() >= deadline {
            let _ = child.kill();
            panic!("runall did not exit within 60s on truncated input (hang?)");
        }
        std::thread::sleep(Duration::from_millis(20));
    };
    let stderr = reader.join().expect("join stderr reader");

    assert!(!status.success(), "runall must fail on a truncated BAM; stderr:\n{stderr}");
    assert!(
        !stderr.contains("panicked") && !stderr.contains("unreachable"),
        "runall must fail cleanly (no panic/unreachable) on truncated input; stderr:\n{stderr}"
    );
    // A clean failure must not leave a half-written output masquerading as success.
    // (Absent or present-but-the-command-failed is fine; we only forbid a
    // success-looking exit, asserted above.)
}

/// Write a gzip-compressed FASTQ from `(name, seq, qual)` records.
fn write_gzip_fastq(path: &Path, records: &[(&str, &str, &str)]) {
    use std::io::Write as _;
    let file = fs::File::create(path).expect("create gzip fastq");
    let mut encoder = flate2::write::GzEncoder::new(file, flate2::Compression::default());
    for (name, seq, qual) in records {
        writeln!(encoder, "@{name}\n{seq}\n+\n{qual}").expect("write fastq record");
    }
    encoder.finish().expect("finish gzip fastq");
}

/// Build a small paired gzip FASTQ (`r1.fq.gz`, `r2.fq.gz`) for the X3-002
/// fusion test: 3 families × 3 read-pairs, BOTH mates carrying a 6 bp UMI
/// prefix (read structure `6M16T` each) so extract builds a paired UMI
/// (`<r1umi>-<r2umi>`) — the shape `--group::strategy paired` + duplex
/// consensus consume. Within a family all pairs share the same UMI + template
/// bases so they collapse into one MI group. Returns `(r1, r2)`.
fn write_paired_umi_fastq(dir: &Path) -> (PathBuf, PathBuf) {
    let r1 = dir.join("r1.fq.gz");
    let r2 = dir.join("r2.fq.gz");
    // (r1 umi, r2 umi, r1 template (16bp), r2 template (16bp)).
    let families = [
        ("AAAAAA", "CCCCCC", "ACGTACGTACGTACGT", "GGTTAACCGGTTAACC"),
        ("GGGGGG", "TTTTTT", "TGCATGCATGCATGCA", "AACCGGTTAACCGGTT"),
        ("ACACAC", "TGTGTG", "TTGGCCAATTGGCCAA", "CCAATTGGCCAATTGG"),
    ];
    let qual = "I".repeat(22); // 6 (UMI) + 16 (template), both mates
    let mut r1_owned: Vec<(String, String)> = Vec::new();
    let mut r2_owned: Vec<(String, String)> = Vec::new();
    for (fi, (r1_umi, r2_umi, r1_tmpl, r2_tmpl)) in families.iter().enumerate() {
        for i in 0..3 {
            let name = format!("fam{fi}_{i}");
            r1_owned.push((name.clone(), format!("{r1_umi}{r1_tmpl}"))); // 6 + 16 = 22
            r2_owned.push((name, format!("{r2_umi}{r2_tmpl}"))); // 6 + 16 = 22
        }
    }
    let r1_slices: Vec<(&str, &str, &str)> =
        r1_owned.iter().map(|(n, s)| (n.as_str(), s.as_str(), qual.as_str())).collect();
    let r2_slices: Vec<(&str, &str, &str)> =
        r2_owned.iter().map(|(n, s)| (n.as_str(), s.as_str(), qual.as_str())).collect();
    write_gzip_fastq(&r1, &r1_slices);
    write_gzip_fastq(&r2, &r2_slices);
    (r1, r2)
}

/// X3-002: the headline FASTQ→consensus fusion. `fgumi runall --start-from
/// extract` over a paired gzip FASTQ must produce the SAME (non-empty) record
/// stream as the equivalent chain where extract is run STANDALONE and its BAM
/// is fed to `runall --start-from align`. Every other runall fixture is a
/// synthesized BAM; this gives the fused extract source (`SourceSpec::Fastqs`)
/// and the paired-FASTQ pull machinery its only end-to-end coverage.
///
/// NB: runall AUTO-INSERTS `Stage::Align` after extract (a FASTQ source must be
/// aligned before sort/group), so `--start-from extract` requires `--ref` + an
/// aligner. We use the deterministic mock-aligner shell script (as the AAM
/// tests do) so the test needs no real bwa. The staged side feeds the SAME
/// extracted BAM into `--start-from align` with the SAME mock aligner, isolating
/// the one difference under test: whether extract is fused into the chain.
#[test]
fn runall_extract_to_duplex_fastq_fusion_matches_staged() {
    let tmp = TempDir::new().unwrap();
    let reference = crate::helpers::bam_generator::create_test_reference(tmp.path());
    // PAIRED mock aligner: extract emits paired R1/R2 sharing a queryname, so a
    // single-end mock (flag 0 on both) would collide on the queryname. The
    // paired mock alternates R1(99)/R2(147) as proper pairs with MC tags, the
    // shape `GroupByPosition` + duplex consensus need.
    let mock = write_paired_mock_aligner_script(tmp.path());
    let aligner_cmd = format!("{} {{ref}}", mock.display());

    let (r1, r2) = write_paired_umi_fastq(tmp.path());

    // ── Fused: runall --start-from extract (auto-inserts align) → simplex. ──
    let fused_out = tmp.path().join("fused_extract_simplex.bam");
    let fused = fgumi(&[
        "runall".into(),
        "--start-from".into(),
        "extract".into(),
        "--stop-after".into(),
        "consensus".into(),
        "--consensus".into(),
        "duplex".into(),
        "--ref".into(),
        reference.as_os_str().to_owned(),
        "--aligner::command".into(),
        aligner_cmd.clone().into(),
        "--extract::inputs".into(),
        r1.as_os_str().to_owned(),
        r2.as_os_str().to_owned(),
        "--extract::read-structures".into(),
        "6M16T".into(),
        "6M16T".into(),
        "--extract::sample".into(),
        "s".into(),
        "--extract::library".into(),
        "l".into(),
        "--group::strategy".into(),
        "paired".into(),
        "--group::edits".into(),
        "0".into(),
        "--duplex::min-reads".into(),
        "1,1,0".into(),
        "--output".into(),
        fused_out.as_os_str().to_owned(),
        "--threads".into(),
        "1".into(),
    ]);
    assert!(
        fused.status.success(),
        "fused extract→duplex: {}",
        String::from_utf8_lossy(&fused.stderr)
    );

    // ── Staged: standalone extract, then runall --start-from align on the
    // extracted BAM with the SAME mock aligner. ──
    let extracted = tmp.path().join("extracted.bam");
    let e = fgumi(&[
        "extract".into(),
        "--inputs".into(),
        r1.as_os_str().to_owned(),
        r2.as_os_str().to_owned(),
        "--read-structures".into(),
        "6M16T".into(),
        "6M16T".into(),
        "--sample".into(),
        "s".into(),
        "--library".into(),
        "l".into(),
        "--output".into(),
        extracted.as_os_str().to_owned(),
    ]);
    assert!(e.status.success(), "staged extract: {}", String::from_utf8_lossy(&e.stderr));

    let staged_out = tmp.path().join("staged_align_duplex.bam");
    let staged = fgumi(&[
        "runall".into(),
        "--start-from".into(),
        "align".into(),
        "--stop-after".into(),
        "consensus".into(),
        "--consensus".into(),
        "duplex".into(),
        "--input".into(),
        extracted.as_os_str().to_owned(),
        "--ref".into(),
        reference.as_os_str().to_owned(),
        "--aligner::command".into(),
        aligner_cmd.into(),
        "--group::strategy".into(),
        "paired".into(),
        "--group::edits".into(),
        "0".into(),
        "--duplex::min-reads".into(),
        "1,1,0".into(),
        "--output".into(),
        staged_out.as_os_str().to_owned(),
        "--threads".into(),
        "1".into(),
    ]);
    assert!(
        staged.status.success(),
        "staged align→duplex: {}",
        String::from_utf8_lossy(&staged.stderr)
    );

    // The fused FASTQ→consensus chain must match the standalone-extract +
    // align-onward chain, and be non-empty (the FASTQ entry path actually
    // produced consensus records through the full extract→align→…→simplex path).
    assert_bams_record_equivalent_nonempty(&fused_out, &staged_out);
    // S9a-005: the FASTQ-entrypoint chain also builds the output header from
    // scratch (extract→align→…), so pin @HD/@SQ/@RG parity against the staged chain.
    assert_bam_headers_equivalent_ignoring_pg(&fused_out, &staged_out);
}

// ────────────────────────── CLI-level rejection tests ──────────────────────────
//
// These pin the validator's structural rules at the binary entry
// point (clap parsing + `RunAll::execute`'s `validate_stages` call),
// catching regressions that the unit-level `validate_with_*` tests in
// `commands::runall::tests` wouldn't surface (e.g. a clap-level flag
// rename that silently drops `--start-from`).
//
// We test three STRUCTURALLY DISTINCT rejection causes, each pinned by its
// OWN stderr fragment (S9a-010 — the shared `"--start-from"` substring did not
// distinguish them; that flag name appears in essentially any error that
// echoes it, so it could not tell a correct rejection from an unrelated one):
//   * Invalid start value (`--start-from simplex` / `duplex`) — the
//     per-algorithm names are NOT valid stage values (the algorithm moved to
//     `--consensus`), so CLAP rejects them with `invalid value '<name>'`
//     BEFORE the validator runs.
//   * Backwards order (`group → sort`) — `group.ord() > sort.ord()`, rejected
//     by `RunAllStage::validate_with` with the `comes after` ordinal message.

fn assert_runall_rejects(start: Stage, stop: Stage, expected_stderr_fragment: &str) {
    let tmp = TempDir::new().unwrap();
    let fixture = unsorted_simplex_fixture(tmp.path());
    let out = tmp.path().join("out.bam");
    let _ = fgumi_binary();
    let args: Vec<OsString> = vec![
        "runall".into(),
        "--start-from".into(),
        start.cli_value().into(),
        "--stop-after".into(),
        stop.cli_value().into(),
        "--input".into(),
        fixture.as_os_str().to_owned(),
        "--output".into(),
        out.as_os_str().to_owned(),
        "--group::strategy".into(),
        "adjacency".into(),
        "--group::edits".into(),
        "1".into(),
        "--threads".into(),
        "1".into(),
    ];
    let output = fgumi(&args);
    assert!(
        !output.status.success(),
        "expected runall {start:?} → {stop:?} to FAIL but it succeeded; stderr={}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains(expected_stderr_fragment),
        "stderr did not contain {expected_stderr_fragment:?}; got: {stderr}"
    );
}

/// `--start-from simplex` is not a valid stage value — the per-algorithm
/// names moved to `--consensus`. Clap rejects it with `invalid value
/// 'simplex'` before any validator rule runs. (S9a-010: pin the clap
/// value-rejection message, not the bare `--start-from` flag name.)
#[cfg(feature = "consensus")]
#[test]
fn rejects_cross_flavor_consensus_pair() {
    assert_runall_rejects(Stage::Simplex, Stage::Duplex, "invalid value 'simplex'");
}

/// `--start-from duplex` is likewise an invalid stage value → clap
/// `invalid value 'duplex'`.
#[cfg(feature = "consensus")]
#[test]
fn rejects_reverse_order_consensus_to_consensus() {
    assert_runall_rejects(Stage::Duplex, Stage::Simplex, "invalid value 'duplex'");
}

/// Backwards order: `group → sort` has `group.ord() > sort.ord()`, rejected
/// by the validator's ordinal rule. Pin the `comes after` message so this
/// test verifies the ORDER rule fired, not merely that the flag was echoed.
#[test]
fn rejects_backwards_group_to_sort() {
    assert_runall_rejects(
        Stage::Group,
        Stage::Sort,
        "--start-from group comes after --stop-after sort",
    );
}

/// S5c2-003: runall `--group::*` flag combos that standalone `fgumi group`
/// rejects must also be rejected on the fused path. This uses standalone
/// `fgumi group` (with the same flags, `--group::` prefix stripped) as the
/// independent ORACLE: it asserts standalone rejects the combo with
/// `expected_fragment`, then asserts runall rejects it the same way. Comparing
/// against the live standalone failure (not just a hard-coded substring) keeps
/// the two code paths aligned if the standalone error text ever changes.
fn assert_runall_group_combo_rejects(extra_group_flags: &[&str], expected_fragment: &str) {
    let tmp = TempDir::new().unwrap();
    let fixture = sorted_duplex_fixture(tmp.path());

    // ── Oracle: standalone `fgumi group` with the un-prefixed flags. ──
    let standalone_out = tmp.path().join("standalone_out.bam");
    let mut standalone_args: Vec<OsString> = vec![
        "group".into(),
        "--input".into(),
        fixture.as_os_str().to_owned(),
        "--output".into(),
        standalone_out.as_os_str().to_owned(),
        "--threads".into(),
        "1".into(),
    ];
    for f in extra_group_flags {
        // Translate the runall-prefixed `--group::<flag>` into the standalone
        // `--<flag>`; non-prefixed tokens (values) pass through unchanged.
        let translated = f
            .strip_prefix("--group::")
            .map_or_else(|| (*f).to_string(), |bare| format!("--{bare}"));
        standalone_args.push(translated.into());
    }
    let standalone = fgumi(&standalone_args);
    let standalone_stderr = String::from_utf8_lossy(&standalone.stderr);
    assert!(
        !standalone.status.success(),
        "oracle: expected standalone group combo {extra_group_flags:?} to FAIL but it \
         succeeded; stderr={standalone_stderr}"
    );
    assert!(
        standalone_stderr.contains(expected_fragment),
        "oracle: standalone stderr did not contain {expected_fragment:?}; got: {standalone_stderr}"
    );

    // ── Runall must reject the same combo with the same error fragment. ──
    let out = tmp.path().join("out.bam");
    let mut args: Vec<OsString> = vec![
        "runall".into(),
        "--start-from".into(),
        Stage::Group.runall_stage_value().into(),
        "--stop-after".into(),
        Stage::Group.runall_stage_value().into(),
        "--input".into(),
        fixture.as_os_str().to_owned(),
        "--output".into(),
        out.as_os_str().to_owned(),
        "--threads".into(),
        "1".into(),
    ];
    for f in extra_group_flags {
        args.push((*f).into());
    }
    let output = fgumi(&args);
    assert!(
        !output.status.success(),
        "expected runall group combo {extra_group_flags:?} to FAIL but it succeeded; stderr={}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains(expected_fragment),
        "runall stderr did not contain {expected_fragment:?} (the standalone oracle rejected \
         with it); got: {stderr}"
    );
}

#[rstest]
#[case(
    &["--group::strategy", "paired", "--group::min-umi-length", "4"],
    "Paired strategy cannot be used with --min-umi-length"
)]
#[case(
    &["--group::strategy", "paired", "--group::no-umi", "true"],
    "--no-umi cannot be used with --strategy paired"
)]
fn rejects_invalid_group_combo(
    #[case] extra_group_flags: &[&str],
    #[case] expected_fragment: &str,
) {
    assert_runall_group_combo_rejects(extra_group_flags, expected_fragment);
}

// ────────────────────────── Zipper fixture + tests ──────────────────────────
//
// Zipper takes a (mapped, unmapped, reference) triple. The runall
// dispatch for `--start-from zipper` builds a Zipper, runs it, then
// chains downstream stages through tempfiles via
// `execute_zipper_then_downstream`. Class A parity here pins the
// delegation glue; the chained downstream cases get smoke coverage
// (success + non-empty BAM) — a full per-stage parity sweep would
// duplicate the work the upstream parity tests already cover.

use crate::helpers::bam_generator::create_test_reference;
use crate::helpers::cli_runner::{ZipperInputs, run_runall_zipper, run_standalone_zipper};
use fgumi_lib::sam::SamTag;
use fgumi_raw_bam::{RawRecord, SamBuilder, flags};
use std::path::Path as StdPath;

/// Build a queryname-sorted zipper fixture: an unmapped BAM with
/// RX/QX tags + a mapped BAM with the same read names aligned to
/// chr1 + a reference FASTA (with `.dict`).
fn zipper_fixture(dir: &StdPath) -> ZipperInputs {
    let unmapped_path = dir.join("zipper_unmapped.bam");
    let mapped_path = dir.join("zipper_mapped.bam");
    let reference = create_test_reference(dir);

    let names = ["read1", "read2", "read3", "read4"];

    // Unmapped: each read has an RX/QX tag, no alignment.
    let unmapped_records: Vec<RawRecord> = names
        .iter()
        .enumerate()
        .map(|(i, name)| {
            let mut b = SamBuilder::new();
            b.read_name(name.as_bytes())
                .sequence(b"ACGTACGT")
                .qualities(&[30; 8])
                .flags(flags::UNMAPPED)
                .add_string_tag(SamTag::RX, format!("AACC{i}").as_bytes())
                .add_string_tag(SamTag::QX, b"IIIIIIII");
            b.build()
        })
        .collect();
    write_qname_sorted_bam(&unmapped_path, &noodles::sam::Header::default(), &unmapped_records);

    // Mapped: same read names, aligned at distinct positions on chr1.
    let mapped_header = create_minimal_header("chr1", 10000);
    let mapped_records: Vec<RawRecord> = names
        .iter()
        .enumerate()
        .map(|(i, name)| {
            let mut b = SamBuilder::new();
            b.read_name(name.as_bytes())
                .ref_id(0)
                .pos(i32::try_from(100 + 50 * i).unwrap_or(0))
                .mapq(60)
                .flags(0)
                .cigar_ops(&[8u32 << 4])
                .sequence(b"ACGTACGT")
                .qualities(&[30; 8]);
            b.build()
        })
        .collect();
    write_qname_sorted_bam(&mapped_path, &mapped_header, &mapped_records);

    ZipperInputs { mapped: mapped_path, unmapped: unmapped_path, reference }
}

fn write_qname_sorted_bam(path: &StdPath, header: &noodles::sam::Header, records: &[RawRecord]) {
    let mut writer = bam::io::Writer::new(fs::File::create(path).expect("create zipper fixture"));
    writer.write_header(header).expect("write fixture header");
    for raw in records {
        writer.write_alignment_record(header, &to_record_buf(raw)).expect("write fixture record");
    }
    writer.try_finish().expect("finish zipper fixture");
}

/// Build a zipper fixture that supports a downstream **duplex**
/// consensus call.
///
/// The default `zipper_fixture` (4 unpaired single-end records with
/// distinct RX tags) is unsuitable for duplex consensus: even a
/// correctly-wired chain produces 0 consensus reads from it because
/// paired-strategy grouping ends up with `0/A`, `1/A`, `2/A`, `3/A`
/// (no `/B` partners). That made the wire-up bug in
/// `execute_duplex_pipeline` invisible to
/// `zipper_to_duplex_fused_pipeline` (status passed; header-only BAM
/// satisfied `meta.len() > 0`) until a 4-stage real-data bench
/// surfaced it.
///
/// This fixture creates one duplex molecule with two strands × two
/// paired-end read pairs per strand (8 records per BAM, 4 per
/// strand). The paired-UMI RX values `"AAAA-CCCC"` (strand A) and
/// `"CCCC-AAAA"` (strand B) are the format `fgumi group --strategy
/// paired` expects: after grouping the records carry MI tags `0/A`
/// and `0/B`, and `fgumi duplex --min-reads 1,1,0` produces a single
/// consensus read.
///
/// Both BAMs are queryname-sorted (alphabetical: `mol_a_0`,
/// `mol_a_1`, `mol_b_0`, `mol_b_1`) and have matching read names so
/// zipper merge can pair them. The mapped BAM carries
/// `MC` (mate-cigar) tags — `GroupByPosition` requires these on
/// paired-end inputs so template span can be computed without a
/// second pass.
fn zipper_duplex_fixture(dir: &StdPath) -> ZipperInputs {
    use bstr::ByteSlice;
    let unmapped_path = dir.join("zipper_duplex_unmapped.bam");
    let mapped_path = dir.join("zipper_duplex_mapped.bam");
    let reference = create_test_reference(dir);

    // One molecule with two strands. Paired UMIs reverse each other —
    // that's what `paired` strategy uses to recognise A/B duplex
    // partners.
    let pairs = [
        ("mol_a_0", "AAAA-CCCC"),
        ("mol_a_1", "AAAA-CCCC"),
        ("mol_b_0", "CCCC-AAAA"),
        ("mol_b_1", "CCCC-AAAA"),
    ];

    // R1 covers chr1:100-107 (8 bp), R2 covers chr1:200-207 (8 bp).
    // Template length: from R1 start (100) to R2 end (208) = 108.
    let r1_seq: &[u8] = b"ACGTACGT";
    let r2_seq: &[u8] = b"TTTTAAAA";
    let r1_cigar: u32 = 8u32 << 4; // 8M
    let r2_cigar: u32 = 8u32 << 4;
    let r1_mc = b"8M";
    let r2_mc = b"8M";
    let tlen: i32 = 108;

    let mut unmapped_records: Vec<RawRecord> = Vec::new();
    let mut mapped_records: Vec<RawRecord> = Vec::new();

    for (name, umi) in &pairs {
        // R1 — unmapped variant: PAIRED | UNMAPPED | FIRST_SEGMENT,
        // RX/QX tags present, no alignment fields.
        let mut u1 = SamBuilder::new();
        u1.read_name(name.as_bytes())
            .sequence(r1_seq)
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::UNMAPPED | flags::MATE_UNMAPPED | flags::FIRST_SEGMENT)
            .add_string_tag(SamTag::RX, umi.as_bytes())
            .add_string_tag(SamTag::QX, b"IIIIIIII");
        unmapped_records.push(u1.build());

        let mut u2 = SamBuilder::new();
        u2.read_name(name.as_bytes())
            .sequence(r2_seq)
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::UNMAPPED | flags::MATE_UNMAPPED | flags::LAST_SEGMENT)
            .add_string_tag(SamTag::RX, umi.as_bytes())
            .add_string_tag(SamTag::QX, b"IIIIIIII");
        unmapped_records.push(u2.build());

        // R1 — mapped variant: PAIRED | FIRST_SEGMENT, aligned at
        // chr1:99 (0-based 99 = 1-based 100), mate at 199. RX
        // intentionally absent — zipper merges it in from the
        // unmapped side.
        let mut m1 = SamBuilder::new();
        m1.read_name(name.as_bytes())
            .ref_id(0)
            .pos(99)
            .mapq(60)
            .flags(flags::PAIRED | flags::FIRST_SEGMENT)
            .mate_ref_id(0)
            .mate_pos(199)
            .template_length(tlen)
            .cigar_ops(&[r1_cigar])
            .sequence(r1_seq)
            .qualities(&[30; 8])
            .add_string_tag(SamTag::MC, r1_mc.as_bstr());
        mapped_records.push(m1.build());

        let mut m2 = SamBuilder::new();
        m2.read_name(name.as_bytes())
            .ref_id(0)
            .pos(199)
            .mapq(60)
            .flags(flags::PAIRED | flags::LAST_SEGMENT)
            .mate_ref_id(0)
            .mate_pos(99)
            .template_length(-tlen)
            .cigar_ops(&[r2_cigar])
            .sequence(r2_seq)
            .qualities(&[30; 8])
            .add_string_tag(SamTag::MC, r2_mc.as_bstr());
        mapped_records.push(m2.build());
    }

    write_qname_sorted_bam(&unmapped_path, &noodles::sam::Header::default(), &unmapped_records);
    let mapped_header = create_minimal_header("chr1", 10000);
    write_qname_sorted_bam(&mapped_path, &mapped_header, &mapped_records);

    ZipperInputs { mapped: mapped_path, unmapped: unmapped_path, reference }
}

/// Build a zipper fixture that supports a downstream **CODEC**
/// consensus call.
///
/// CODEC recognises duplex pairs by **FR-overlap flag shape** (R1
/// `PAIRED | FIRST_SEGMENT | MATE_REVERSE`, R2 `PAIRED |
/// LAST_SEGMENT | REVERSE`, both at the SAME position), NOT by
/// paired-UMI grouping. Without those flags, CODEC consensus
/// rejects every pair and emits 0 records — which is what
/// `parity_b_sort_to_codec`'s `unsorted_codec_fixture` actually
/// does in practice (the existing test passes only because both
/// fused and staged produce 0 records, so they're trivially
/// equivalent).
///
/// This fixture mirrors the `create_codec_read_pair` helper from
/// `tests/integration/test_codec_command.rs`, which IS known to
/// produce non-zero CODEC output: same FR-overlap flags, same
/// position for R1 + R2, single-strand UMI in RX, 3 pairs per
/// family for sufficient depth.
///
/// After `zipper → sort → group(adjacency) → codec`, the chain
/// produces ≥ 1 consensus.
fn zipper_codec_fixture(dir: &StdPath) -> ZipperInputs {
    use bstr::ByteSlice;
    let unmapped_path = dir.join("zipper_codec_unmapped.bam");
    let mapped_path = dir.join("zipper_codec_mapped.bam");
    let reference = create_test_reference(dir);

    let r1_seq: &[u8] = b"ACGTACGT";
    let r2_seq: &[u8] = b"ACGTACGT"; // matching R1 (FR overlap pair)
    let r1_cigar: u32 = 8u32 << 4;
    let r2_cigar: u32 = 8u32 << 4;
    let r1_mc = b"8M";
    let r2_mc = b"8M";
    let pos: i32 = 99;
    let umi = "ACGTACGT";

    let mut unmapped_records: Vec<RawRecord> = Vec::new();
    let mut mapped_records: Vec<RawRecord> = Vec::new();

    for i in 0..3 {
        let name = format!("read_{i}");

        let mut u1 = SamBuilder::new();
        u1.read_name(name.as_bytes())
            .sequence(r1_seq)
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::UNMAPPED | flags::MATE_UNMAPPED | flags::FIRST_SEGMENT)
            .add_string_tag(SamTag::RX, umi.as_bytes())
            .add_string_tag(SamTag::QX, b"IIIIIIII");
        unmapped_records.push(u1.build());

        let mut u2 = SamBuilder::new();
        u2.read_name(name.as_bytes())
            .sequence(r2_seq)
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::UNMAPPED | flags::MATE_UNMAPPED | flags::LAST_SEGMENT)
            .add_string_tag(SamTag::RX, umi.as_bytes())
            .add_string_tag(SamTag::QX, b"IIIIIIII");
        unmapped_records.push(u2.build());

        // R1: PAIRED | FIRST_SEGMENT | MATE_REVERSE — the
        // `MATE_REVERSE` flag tells codec the R2 partner is on
        // the reverse strand (FR-overlap shape).
        let mut m1 = SamBuilder::new();
        m1.read_name(name.as_bytes())
            .ref_id(0)
            .pos(pos)
            .mapq(60)
            .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE)
            .mate_ref_id(0)
            .mate_pos(pos)
            .template_length(i32::try_from(r1_seq.len()).expect("r1_seq length fits i32"))
            .cigar_ops(&[r1_cigar])
            .sequence(r1_seq)
            .qualities(&[30; 8])
            .add_string_tag(SamTag::MC, r1_mc.as_bstr());
        mapped_records.push(m1.build());

        // R2: PAIRED | LAST_SEGMENT | REVERSE — same pos as R1.
        let mut m2 = SamBuilder::new();
        m2.read_name(name.as_bytes())
            .ref_id(0)
            .pos(pos)
            .mapq(60)
            .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE)
            .mate_ref_id(0)
            .mate_pos(pos)
            .template_length(-i32::try_from(r2_seq.len()).expect("r2_seq length fits i32"))
            .cigar_ops(&[r2_cigar])
            .sequence(r2_seq)
            .qualities(&[30; 8])
            .add_string_tag(SamTag::MC, r2_mc.as_bstr());
        mapped_records.push(m2.build());
    }

    write_qname_sorted_bam(&unmapped_path, &noodles::sam::Header::default(), &unmapped_records);
    let mapped_header = create_minimal_header("chr1", 10000);
    write_qname_sorted_bam(&mapped_path, &mapped_header, &mapped_records);

    ZipperInputs { mapped: mapped_path, unmapped: unmapped_path, reference }
}

/// Class A — `zipper → zipper` parity vs standalone `fgumi zipper`.
/// Pins the runall consensus-only-style delegation: runall constructs
/// a `Zipper` struct from its own flags and calls `Zipper::execute`.
///
/// Both runall and standalone now route through the typed-step
/// pipeline framework (the legacy iterator path was deleted in
/// Phase 1 T1.10), so this end-to-end test verifies the runall
/// dispatcher wires Zipper correctly on the integration fixture's
/// single-end queryname-sorted shape.
#[test]
fn parity_a_zipper_to_zipper() {
    let tmp = TempDir::new().unwrap();
    let inputs = zipper_fixture(tmp.path());
    let runall_out = tmp.path().join("runall.bam");
    let standalone_out = tmp.path().join("standalone.bam");
    let args = ParityArgs::for_simplex_like();

    let r = run_runall_zipper(Stage::Zipper, &inputs, &runall_out, &args);
    assert!(r.status.success(), "runall: {}", String::from_utf8_lossy(&r.stderr));

    let r = run_standalone_zipper(&inputs, &standalone_out, &args);
    assert!(r.status.success(), "standalone: {}", String::from_utf8_lossy(&r.stderr));

    // S9a-001: zipper merges the mapped + unmapped reads → a non-empty stream.
    assert_bams_record_equivalent_nonempty(&runall_out, &standalone_out);
    assert_bam_headers_equivalent_ignoring_pg(&runall_out, &standalone_out);
}

/// Smoke check — `fgumi runall --start-from zipper --stop-after zipper`
/// reaches `Zipper::execute_new_pipeline`. Pinned via the distinctive
/// log line at the typed-step framework entry point
/// (`zipper::NEW_PIPELINE_START_LOG`). The legacy iterator path was
/// deleted in Phase 1 T1.10 so this test now verifies the surviving
/// path runs to completion via the runall fan-in. `RUST_LOG=info` is
/// pinned on the spawned process so an ambient stricter env
/// (e.g. `RUST_LOG=warn`) on the developer's shell cannot silently
/// mask the log line — the binary's default level (set in `main.rs`)
/// is info, but env overrides win.
#[test]
fn runall_zipper_self_pair_uses_new_pipeline() {
    use crate::helpers::cli_runner::{Stage, fgumi_binary};
    use fgumi_lib::commands::zipper::NEW_PIPELINE_START_LOG;
    use std::ffi::OsString;
    use std::process::Command;

    let tmp = TempDir::new().unwrap();
    let inputs = zipper_fixture(tmp.path());
    let out = tmp.path().join("out.bam");
    let args = ParityArgs::for_simplex_like();

    let cli_args: Vec<OsString> = vec![
        "runall".into(),
        "--start-from".into(),
        Stage::Zipper.cli_value().into(),
        "--stop-after".into(),
        Stage::Zipper.cli_value().into(),
        "--input".into(),
        inputs.mapped.as_os_str().to_owned(),
        "--unmapped".into(),
        inputs.unmapped.as_os_str().to_owned(),
        "--ref".into(),
        inputs.reference.as_os_str().to_owned(),
        "--output".into(),
        out.as_os_str().to_owned(),
        "--group::strategy".into(),
        args.strategy.into(),
        "--group::edits".into(),
        args.edits.to_string().into(),
        "--threads".into(),
        args.threads.to_string().into(),
    ];
    let r = Command::new(fgumi_binary())
        .env("RUST_LOG", "info")
        .args(&cli_args)
        .output()
        .expect("failed to execute fgumi");
    assert!(r.status.success(), "runall: {}", String::from_utf8_lossy(&r.stderr));

    let stderr = String::from_utf8_lossy(&r.stderr);
    assert!(
        stderr.contains(NEW_PIPELINE_START_LOG),
        "runall --start-from zipper --stop-after zipper should run the typed-step \
         framework path, but the new-pipeline log line ({NEW_PIPELINE_START_LOG:?}) is \
         absent. Stderr:\n{stderr}"
    );
}

/// Smoke test — `--start-from zipper --stop-after sort` produces a
/// non-empty BAM. This exercises the chained zipper→sort path inside
/// `execute_zipper_then_downstream`, which the per-stage parity
/// tests do not reach.
#[test]
fn smoke_zipper_to_sort_chain() {
    let tmp = TempDir::new().unwrap();
    let inputs = zipper_fixture(tmp.path());
    let runall_out = tmp.path().join("runall.bam");
    let args = ParityArgs::for_simplex_like();

    let r = run_runall_zipper(Stage::Sort, &inputs, &runall_out, &args);
    assert!(r.status.success(), "runall: {}", String::from_utf8_lossy(&r.stderr));
    // Sanity check: BAM has at least one record (zipper transferred
    // the unmapped tags onto the mapped reads).
    let mut reader = bam::io::reader::Builder.build_from_path(&runall_out).expect("open bam");
    let _ = reader.read_header().expect("read header");
    let records = reader
        .record_bufs(&noodles::sam::Header::default())
        .collect::<Result<Vec<_>, _>>()
        .expect("decode records");
    let n = records.len();
    // The zipper fixture has 4 reads and sort doesn't drop records, so the
    // zipper→sort output must carry exactly 4 — pin it rather than `>= 1`.
    assert_eq!(n, 4, "expected exactly 4 records in zipper→sort output, got {n}");
}

/// CLI rejection — `--start-from zipper` without `--unmapped` should
/// fail with a clear message.
#[test]
fn rejects_zipper_without_unmapped() {
    let tmp = TempDir::new().unwrap();
    let inputs = zipper_fixture(tmp.path());
    let out = tmp.path().join("out.bam");
    let args: Vec<std::ffi::OsString> = vec![
        "runall".into(),
        "--start-from".into(),
        "zipper".into(),
        "--stop-after".into(),
        "zipper".into(),
        "--input".into(),
        inputs.mapped.as_os_str().to_owned(),
        // (intentionally omit --unmapped)
        "--ref".into(),
        inputs.reference.as_os_str().to_owned(),
        "--output".into(),
        out.as_os_str().to_owned(),
        "--threads".into(),
        "1".into(),
    ];
    let r = fgumi(&args);
    assert!(!r.status.success(), "expected zipper-without-unmapped to fail");
    let stderr = String::from_utf8_lossy(&r.stderr);
    assert!(stderr.contains("--unmapped"), "stderr should mention --unmapped, got: {stderr}");
}

// ────────────────────────── AlignAndMerge validation tests ──────────────────────────
//
// C3 lands the CLI surface for --start-from align-and-merge but no
// execution path (those go in C4). These tests pin that the CLI flag
// validation fires correctly: missing required flags produce clear
// errors; mutually-exclusive flag combos are rejected; the
// "not-yet-implemented" message is shown when validation passes but
// execution is reached.

fn fgumi_with_args(args: &[std::ffi::OsString]) -> std::process::Output {
    fgumi(args)
}

/// Common base flags for --start-from align-and-merge tests. The
/// caller appends preset/command-specific flags.
fn aam_base_args(input: &Path, ref_path: &Path, output: &Path) -> Vec<std::ffi::OsString> {
    vec![
        "runall".into(),
        "--start-from".into(),
        "align".into(),
        "--stop-after".into(),
        "zipper".into(),
        "--input".into(),
        input.as_os_str().to_owned(),
        "--output".into(),
        output.as_os_str().to_owned(),
        "--ref".into(),
        ref_path.as_os_str().to_owned(),
        "--threads".into(),
        "1".into(),
    ]
}

/// Without `--aligner::preset` or `--aligner::command`, validation
/// errors with a clear "requires one of" message.
#[test]
fn aam_rejects_when_no_preset_or_command() {
    let tmp = TempDir::new().unwrap();
    let inputs = zipper_fixture(tmp.path());
    let out = tmp.path().join("out.bam");
    let args = aam_base_args(&inputs.unmapped, &inputs.reference, &out);
    let r = fgumi_with_args(&args);
    assert!(!r.status.success());
    let stderr = String::from_utf8_lossy(&r.stderr);
    assert!(
        stderr.contains("requires one of `--aligner::preset`"),
        "expected preset-or-command requirement, got: {stderr}"
    );
}

/// `--aligner::preset` and `--aligner::command` together is rejected.
#[test]
fn aam_rejects_preset_and_command_together() {
    let tmp = TempDir::new().unwrap();
    let inputs = zipper_fixture(tmp.path());
    let out = tmp.path().join("out.bam");
    let mut args = aam_base_args(&inputs.unmapped, &inputs.reference, &out);
    args.push("--aligner::preset".into());
    args.push("bwa-mem3".into());
    args.push("--aligner::command".into());
    args.push("bwa-mem3 mem {ref} /dev/stdin".into());
    let r = fgumi_with_args(&args);
    assert!(!r.status.success());
    let stderr = String::from_utf8_lossy(&r.stderr);
    assert!(stderr.contains("mutually exclusive"), "expected exclusivity error, got: {stderr}");
}

/// `--aligner::command` without `{ref}` placeholder is rejected.
#[test]
fn aam_rejects_command_missing_ref_placeholder() {
    let tmp = TempDir::new().unwrap();
    let inputs = zipper_fixture(tmp.path());
    let out = tmp.path().join("out.bam");
    let mut args = aam_base_args(&inputs.unmapped, &inputs.reference, &out);
    args.push("--aligner::command".into());
    args.push("bwa-mem3 mem -t 8 /dev/stdin".into());
    let r = fgumi_with_args(&args);
    assert!(!r.status.success());
    let stderr = String::from_utf8_lossy(&r.stderr);
    assert!(stderr.contains("{ref}"), "expected {{ref}} placeholder requirement, got: {stderr}");
}

/// `--aligner-bin` with `--aligner::command` is rejected
/// (preset-mode-only knob).
#[test]
fn aam_rejects_aligner_bin_with_command_mode() {
    let tmp = TempDir::new().unwrap();
    let inputs = zipper_fixture(tmp.path());
    let out = tmp.path().join("out.bam");
    let mut args = aam_base_args(&inputs.unmapped, &inputs.reference, &out);
    args.push("--aligner-bin".into());
    args.push("/usr/local/bin/bwa-mem3".into());
    args.push("--aligner::command".into());
    args.push("bwa-mem3 mem {ref} /dev/stdin".into());
    let r = fgumi_with_args(&args);
    assert!(!r.status.success());
    let stderr = String::from_utf8_lossy(&r.stderr);
    assert!(
        stderr.contains("--aligner-bin is only valid with --aligner::preset"),
        "expected preset-only error, got: {stderr}"
    );
}

/// `--start-from align-and-merge` without `--ref` is rejected.
#[test]
fn aam_rejects_when_ref_missing() {
    let tmp = TempDir::new().unwrap();
    let inputs = zipper_fixture(tmp.path());
    let out = tmp.path().join("out.bam");
    let args: Vec<std::ffi::OsString> = vec![
        "runall".into(),
        "--start-from".into(),
        "align".into(),
        "--stop-after".into(),
        "zipper".into(),
        "--input".into(),
        inputs.mapped.as_os_str().to_owned(),
        "--output".into(),
        out.as_os_str().to_owned(),
        "--threads".into(),
        "1".into(),
        "--aligner::command".into(),
        "bwa-mem3 mem {ref} /dev/stdin".into(),
    ];
    let r = fgumi_with_args(&args);
    assert!(!r.status.success());
    let stderr = String::from_utf8_lossy(&r.stderr);
    assert!(stderr.contains("--ref"), "expected --ref requirement, got: {stderr}");
}

/// Build a mock-aligner shell script that reads FASTQ from stdin and
/// emits a valid SAM file to stdout. Each input record produces one
/// mapped SAM record on `chr1` at position 1000 (or 1500 for /2),
/// with a CIGAR matching the sequence length. Suitable for end-to-end
/// AAM tests that don't have a real aligner available.
/// Like `write_mock_aligner_script`, but assumes the input FASTQ is
/// strict R1/R2 alternation and emits paired-end SAM records: R1 with
/// flag 99 (`PAIRED|PROPER_PAIR|MATE_REVERSE|FIRST_IN_PAIR`), R2 with
/// flag 147 (`PAIRED|PROPER_PAIR|REVERSE|SECOND_IN_PAIR`). Use with
/// fixtures whose records are paired-end (`zipper_duplex_fixture`,
/// `zipper_codec_fixture`). Single-end fixtures should keep the
/// default `write_mock_aligner_script` — pairing them arbitrarily
/// here would inject paired-end flags onto records that don't
/// actually pair, confusing downstream group/consensus steps.
fn write_paired_mock_aligner_script(dir: &Path) -> PathBuf {
    let path = dir.join("mock-aligner-pe.sh");
    let script = r#"#!/bin/bash
set -euo pipefail
awk '
BEGIN {
    print "@HD\tVN:1.6"
    print "@SQ\tSN:chr1\tLN:10000"
    nread = 0
}
NR%4==1 { name = substr($1, 2); nread++ }
NR%4==2 {
    seq = $1
    len = length(seq)
    qual = ""
    for (i = 0; i < len; i++) qual = qual "I"
    pos = 1000
    cigar = len "M"
    # Alternate R1 (flag 99) and R2 (flag 147), same name, mate at
    # same chrom/pos, tlen=+/-len. PROPER_PAIR + MC tag so downstream
    # `GroupByPosition` accepts the pair (template-span computation
    # requires MC on paired-end inputs).
    if (nread % 2 == 1) {
        printf "%s\t99\tchr1\t%d\t60\t%s\t=\t%d\t%d\t%s\t%s\tMC:Z:%s\n", name, pos, cigar, pos, len, seq, qual, cigar
    } else {
        printf "%s\t147\tchr1\t%d\t60\t%s\t=\t%d\t-%d\t%s\t%s\tMC:Z:%s\n", name, pos, cigar, pos, len, seq, qual, cigar
    }
}
'
"#;
    fs::write(&path, script).expect("write paired mock aligner script");
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        let mut perms = fs::metadata(&path).unwrap().permissions();
        perms.set_mode(0o755);
        fs::set_permissions(&path, perms).unwrap();
    }
    path
}

fn write_mock_aligner_script(dir: &Path) -> PathBuf {
    let path = dir.join("mock-aligner.sh");
    // Mock aligner: reads FASTQ from stdin (with `no_suffix = true`
    // — fgumi's AAM FASTQ writer no longer appends `/1` / `/2`
    // suffixes) and emits a valid SAM record per read. Output flag
    // is 0 (mapped, primary, unpaired) so queryname pairing
    // matches the original unmapped BAM byte-for-byte.
    //
    // Uses portable POSIX `awk` (mawk / gawk / BSD awk all handle
    // these features identically); test fixtures avoid `mawk` as a
    // hard dependency because CI images vary.
    let script = r#"#!/bin/bash
set -euo pipefail
awk '
BEGIN { print "@HD\tVN:1.6"; print "@SQ\tSN:chr1\tLN:10000" }
NR%4==1 { name = substr($1, 2) }
NR%4==2 {
    seq = $1
    len = length(seq)
    pos = 1000
    # Print SAM: QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL
    # Mapped, primary, unpaired (flag = 0). Quality = "I" * len.
    qual = ""
    for (i = 0; i < len; i++) qual = qual "I"
    printf "%s\t0\tchr1\t%d\t60\t%dM\t*\t0\t0\t%s\t%s\n", name, pos, len, seq, qual
}
'
"#;
    fs::write(&path, script).expect("write mock aligner script");
    // chmod +x so /bin/bash -c can run it.
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        let mut perms = fs::metadata(&path).unwrap().permissions();
        perms.set_mode(0o755);
        fs::set_permissions(&path, perms).unwrap();
    }
    path
}

/// End-to-end smoke test for `--start-from align-and-merge --stop-after
/// zipper`. Uses a hand-rolled shell script as the "aligner" so we
/// can exercise the full AAM chain (BAM → FASTQ → aligner → zipper)
/// without depending on a real `bwa-mem3` / `bwa` binary in the
/// developer's environment. Real-aligner parity tests come in C5
/// when CI installs the bwa-mem3 / bwa binaries from bioconda.
#[test]
fn aam_end_to_end_with_mock_aligner() {
    let tmp = TempDir::new().unwrap();
    let inputs = zipper_fixture(tmp.path());
    let out = tmp.path().join("aam-merged.bam");
    let mock = write_mock_aligner_script(tmp.path());

    let mut args = aam_base_args(&inputs.unmapped, &inputs.reference, &out);
    // Run the mock aligner script as the aligner. `cat` discards the
    // path positional arg — our mock reads FASTQ from stdin via the
    // awk pipe, and the trailing `{ref}` is just present to satisfy
    // the placeholder requirement.
    args.push("--aligner::command".into());
    args.push(format!("{} {{ref}}", mock.display()).into());
    let r = fgumi_with_args(&args);
    assert!(
        r.status.success(),
        "AAM end-to-end (mock aligner) failed: stdout={} stderr={}",
        String::from_utf8_lossy(&r.stdout),
        String::from_utf8_lossy(&r.stderr)
    );
    // Output BAM must exist and be non-empty (the merge produced at
    // least the header + some records).
    let meta = fs::metadata(&out).expect("merged BAM must exist");
    assert!(meta.len() > 0, "merged BAM was empty");
    // Sanity check: at least one record made it through merge.
    let mut reader = bam::io::reader::Builder.build_from_path(&out).expect("open merged BAM");
    let _ = reader.read_header().expect("read merged header");
    let n_records = reader.record_bufs(&noodles::sam::Header::default()).count();
    assert!(n_records >= 1, "expected at least 1 record in merged BAM, got {n_records}");
}

/// Companion to `aam_end_to_end_with_mock_aligner` exercising the
/// new fused `--start-from align --stop-after sort` path. This is
/// the regression gate for the "no `aam-merged.bam` tempfile
/// bridge" change: a runall invocation that previously wrote two
/// BAM files (aam-merged.bam + sorted-output) now writes one in a
/// single `Pipeline::run`.
#[test]
fn aam_to_sort_fused_pipeline_with_mock_aligner() {
    let tmp = TempDir::new().unwrap();
    let inputs = zipper_fixture(tmp.path());
    let out = tmp.path().join("aam-sorted.bam");
    let mock = write_mock_aligner_script(tmp.path());

    let mut args = aam_base_args(&inputs.unmapped, &inputs.reference, &out);
    // Override the default `--stop-after zipper` baked into
    // `aam_base_args` to exercise the fused AAM→Sort path.
    let stop_idx = args
        .iter()
        .position(|a| a == std::ffi::OsStr::new("--stop-after"))
        .expect("--stop-after flag in aam_base_args");
    args[stop_idx + 1] = "sort".into();
    args.push("--sort::max-memory".into());
    args.push("256M".into());
    args.push("--aligner::command".into());
    args.push(format!("{} {{ref}}", mock.display()).into());
    let r = fgumi_with_args(&args);
    assert!(
        r.status.success(),
        "AAM→Sort fused (mock aligner) failed: stdout={} stderr={}",
        String::from_utf8_lossy(&r.stdout),
        String::from_utf8_lossy(&r.stderr)
    );
    // Output BAM must exist and have at least one record.
    let meta = fs::metadata(&out).expect("sorted BAM must exist");
    assert!(meta.len() > 0, "sorted BAM was empty");
    let mut reader = bam::io::reader::Builder.build_from_path(&out).expect("open sorted BAM");
    let _ = reader.read_header().expect("read sorted header");
    let n_records = reader.record_bufs(&noodles::sam::Header::default()).count();
    assert!(n_records >= 1, "expected at least 1 record in sorted BAM, got {n_records}");

    // The fused path's promise: no AAM/sort tempfile bridge BAM is left
    // behind. S9a-004: assert the EXACT set of top-level `.bam` files —
    // the two `zipper_fixture` inputs plus the user's `--output` — and
    // that NO tempfile-bridge name appears. The old `len() <= 3` bound was
    // self-defeating: 3 legitimate BAMs + 1 leaked `aam-merged.bam` = 4,
    // but a leak of a name already counted (or any off-by-one) could slip
    // past a loose count. The exact allow-set + per-name leak check (the
    // pattern the simplex companion already uses) catches a single leak.
    let bam_names: std::collections::BTreeSet<String> = fs::read_dir(tmp.path())
        .unwrap()
        .filter_map(Result::ok)
        .filter(|e| e.path().extension().is_some_and(|x| x == "bam"))
        .filter_map(|e| e.file_name().to_str().map(str::to_owned))
        .collect();
    let expected: std::collections::BTreeSet<String> =
        ["zipper_mapped.bam", "zipper_unmapped.bam", "aam-sorted.bam"]
            .into_iter()
            .map(str::to_owned)
            .collect();
    assert_eq!(
        bam_names, expected,
        "unexpected top-level BAM set (tempfile leak?): got {bam_names:?}, expected {expected:?}"
    );
    for tempfile_name in ["aam-merged.bam", "sorted.bam", "after_sort.bam", "after_group.bam"] {
        assert!(
            !bam_names.contains(tempfile_name),
            "fused-pipeline tempfile {tempfile_name} leaked into the user-visible dir: {bam_names:?}"
        );
    }
}

/// Companion to `aam_to_sort_fused_pipeline_with_mock_aligner`:
/// exercises the fully-fused `--start-from align --stop-after
/// simplex` chain. This is the regression gate for the "no
/// `sorted.bam` / `after_group.bam` tempfile bridges" change for
/// the consensus path: AAM → Sort → Group → Consensus runs in a
/// single `Pipeline::run` with no intermediate disk hops.
///
/// Mock aligner is deterministic enough to produce stable
/// alignments; we don't assert byte-equivalence against the
/// tempfile-bridge path because the merged BAM's record order
/// depends on the in-pipeline vs across-tempfile sort cadence —
/// instead we assert the consensus BAM is non-empty.
#[cfg(feature = "consensus")]
#[test]
fn aam_to_simplex_fused_pipeline_with_mock_aligner() {
    let tmp = TempDir::new().unwrap();
    let inputs = zipper_fixture(tmp.path());
    let out = tmp.path().join("aam-simplex.bam");
    let mock = write_mock_aligner_script(tmp.path());

    let mut args = aam_base_args(&inputs.unmapped, &inputs.reference, &out);
    let stop_idx = args
        .iter()
        .position(|a| a == std::ffi::OsStr::new("--stop-after"))
        .expect("--stop-after flag in aam_base_args");
    args[stop_idx + 1] = "consensus".into();
    args.push("--consensus".into());
    args.push("simplex".into());
    args.push("--sort::max-memory".into());
    args.push("256M".into());
    args.push("--group::strategy".into());
    args.push("adjacency".into());
    args.push("--simplex::min-reads".into());
    args.push("1".into());
    args.push("--aligner::command".into());
    args.push(format!("{} {{ref}}", mock.display()).into());

    let r = fgumi_with_args(&args);
    assert!(
        r.status.success(),
        "AAM→Simplex fused (mock aligner) failed: stdout={} stderr={}",
        String::from_utf8_lossy(&r.stdout),
        String::from_utf8_lossy(&r.stderr)
    );
    let meta = fs::metadata(&out).expect("simplex BAM must exist");
    assert!(meta.len() > 0, "simplex BAM was empty");
    // Smoke-test the chain end-to-end: the output BAM must open and
    // parse a header. Record-count assertion is left to the
    // paired-end `aam_to_duplex_fused_pipeline_with_mock_aligner_pe`
    // and the higher-fidelity `parity_b_*` tests — the default
    // single-end mock aligner emits flag=0 records, and `fgumi
    // simplex` requires both R1 and R2 of a paired template to emit
    // consensus, so this chain produces 0 records by construction.
    let mut reader = bam::io::reader::Builder.build_from_path(&out).expect("open simplex BAM");
    let _ = reader.read_header().expect("read simplex header");

    // Verify no `sorted.bam` / `after_group.bam` leftover in the
    // user-visible working dir (would indicate the fused path
    // didn't kick in — tempfiles should live in
    // `tempfile::tempdir()` inside the pipeline, never in
    // the user-visible test dir).
    let leftover_names: Vec<String> = fs::read_dir(tmp.path())
        .unwrap()
        .filter_map(Result::ok)
        .filter_map(|e| e.file_name().to_str().map(str::to_owned))
        .collect();
    for tempfile_name in ["sorted.bam", "after_sort.bam", "after_group.bam", "aam-merged.bam"] {
        assert!(
            !leftover_names.iter().any(|n| n == tempfile_name),
            "fused-pipeline tempfile {tempfile_name} leaked into the user-visible dir: {leftover_names:?}"
        );
    }
}

#[cfg(feature = "consensus")]
#[test]
fn aam_to_duplex_fused_pipeline_with_mock_aligner() {
    let tmp = TempDir::new().unwrap();
    // Use the paired-UMI strand A/B fixture and the paired-end mock
    // aligner so the consensus path actually has well-formed input
    // (proper-pair flags, both R1 and R2 with matching queryname).
    // The default single-end mock aligner would emit `flag=0` for
    // every record, which collapses the duplex grouping shape and
    // makes the test header-only — exactly the wire-up-bug gap that
    // motivated 4f7a149.
    //
    // Pass `inputs.unmapped` (not `inputs.mapped`) because AAM uses
    // the input BAM both as the FASTQ source AND as the unmapped
    // half for `merge_raw`. The duplex fixture's mapped half has
    // `RX` deliberately omitted (zipper-merge-from-unmapped is the
    // production shape); the unmapped half carries `RX`. Passing
    // the mapped BAM here produces records with no RX after merge,
    // which collapses paired-UMI grouping to zero groups.
    let inputs = zipper_duplex_fixture(tmp.path());
    let out = tmp.path().join("aam-duplex.bam");
    let mock = write_paired_mock_aligner_script(tmp.path());

    let mut args = aam_base_args(&inputs.unmapped, &inputs.reference, &out);
    let stop_idx = args
        .iter()
        .position(|a| a == std::ffi::OsStr::new("--stop-after"))
        .expect("--stop-after flag in aam_base_args");
    args[stop_idx + 1] = "consensus".into();
    args.push("--consensus".into());
    args.push("duplex".into());
    args.push("--sort::max-memory".into());
    args.push("256M".into());
    args.push("--group::strategy".into());
    args.push("paired".into());
    args.push("--duplex::min-reads".into());
    args.push("1,1,0".into());
    args.push("--aligner::command".into());
    args.push(format!("{} {{ref}}", mock.display()).into());

    let r = fgumi_with_args(&args);
    assert!(
        r.status.success(),
        "AAM→Duplex fused (mock aligner) failed: stdout={} stderr={}",
        String::from_utf8_lossy(&r.stdout),
        String::from_utf8_lossy(&r.stderr)
    );
    let meta = fs::metadata(&out).expect("duplex BAM must exist");
    assert!(meta.len() > 0, "duplex BAM was empty");
    let mut reader = bam::io::reader::Builder.build_from_path(&out).expect("open duplex BAM");
    let _ = reader.read_header().expect("read duplex header");
    let n_records = reader.record_bufs(&noodles::sam::Header::default()).count();
    assert!(n_records >= 1, "expected >= 1 record in duplex BAM, got {n_records}");
}

#[cfg(feature = "consensus")]
#[test]
fn aam_to_codec_fused_pipeline_with_mock_aligner() {
    let tmp = TempDir::new().unwrap();
    // CODEC consensus requires the FR-overlap flag shape
    // (R1 `PAIRED | FIRST_SEGMENT | MATE_REVERSE`, R2 `PAIRED |
    // LAST_SEGMENT | REVERSE`, both at the SAME position). The
    // paired mock aligner here emits the proper-pair shape
    // (`MATE_REVERSE` on R1, `REVERSE` on R2, same chrom/pos), which
    // satisfies the FR-overlap check. Combined with
    // `zipper_codec_fixture`'s 3-pair-per-family input, the chain
    // produces at least one consensus. Pass `inputs.unmapped` for
    // the same reason as the duplex test above — AAM consumes its
    // input both as FASTQ source and as the unmapped half for
    // `merge_raw`; RX lives only on the unmapped side.
    let inputs = zipper_codec_fixture(tmp.path());
    let out = tmp.path().join("aam-codec.bam");
    let mock = write_paired_mock_aligner_script(tmp.path());

    let mut args = aam_base_args(&inputs.unmapped, &inputs.reference, &out);
    let stop_idx = args
        .iter()
        .position(|a| a == std::ffi::OsStr::new("--stop-after"))
        .expect("--stop-after flag in aam_base_args");
    args[stop_idx + 1] = "consensus".into();
    args.push("--consensus".into());
    args.push("codec".into());
    args.push("--sort::max-memory".into());
    args.push("256M".into());
    args.push("--group::strategy".into());
    args.push("adjacency".into());
    args.push("--codec::min-reads".into());
    args.push("1".into());
    args.push("--aligner::command".into());
    args.push(format!("{} {{ref}}", mock.display()).into());

    let r = fgumi_with_args(&args);
    assert!(
        r.status.success(),
        "AAM→Codec fused (mock aligner) failed: stdout={} stderr={}",
        String::from_utf8_lossy(&r.stdout),
        String::from_utf8_lossy(&r.stderr)
    );
    let meta = fs::metadata(&out).expect("codec BAM must exist");
    assert!(meta.len() > 0, "codec BAM was empty");
    let mut reader = bam::io::reader::Builder.build_from_path(&out).expect("open codec BAM");
    let _ = reader.read_header().expect("read codec header");
    let n_records = reader.record_bufs(&noodles::sam::Header::default()).count();
    assert!(n_records >= 1, "expected >= 1 record in codec BAM, got {n_records}");
}

/// Regression gate for the fused `--start-from zipper --stop-after
/// simplex` chain: Zipper → Sort → Group → Consensus runs in a
/// single `Pipeline::run` with no intermediate tempfiles.
///
/// Reuses `zipper_codec_fixture` (FR-overlap paired-end records,
/// single-strand UMI, 3 pairs per family — same shape that
/// `tests/integration/test_codec_command.rs::create_codec_read_pair`
/// is known to produce non-zero consensus on). Simplex `--min-reads
/// 1` emits ≥ 1 consensus per MI group from this fixture.
///
/// The original `zipper_fixture` (4 unpaired single-end records)
/// turns out to produce 0 simplex consensus reads on a correctly-
/// wired chain — so a record-count assertion on it can't tell the
/// wire-up bug shape apart from the "fixture too thin" shape.
#[cfg(feature = "consensus")]
#[test]
fn zipper_to_simplex_fused_pipeline() {
    let tmp = TempDir::new().unwrap();
    let inputs = zipper_codec_fixture(tmp.path());
    let out = tmp.path().join("zipper-simplex.bam");
    let args = ParityArgs::for_simplex_like();
    let r = run_runall_zipper(Stage::Simplex, &inputs, &out, &args);
    assert!(
        r.status.success(),
        "Zipper→Simplex fused failed: stdout={} stderr={}",
        String::from_utf8_lossy(&r.stdout),
        String::from_utf8_lossy(&r.stderr)
    );
    assert!(fs::metadata(&out).expect("simplex BAM").len() > 0);

    let mut reader = bam::io::reader::Builder.build_from_path(&out).expect("open simplex BAM");
    let _ = reader.read_header().expect("read simplex header");
    let n_records = reader.record_bufs(&noodles::sam::Header::default()).count();
    assert!(
        n_records >= 1,
        "zipper→simplex fused chain produced 0 consensus reads — likely a wire-up \
         regression (`execute_simplex_pipeline`'s `else if let Some(zipper) = zipper` \
         arm is missing or unreachable)"
    );

    // No intermediate tempfiles in the user-visible dir.
    let leftover_names: Vec<String> = fs::read_dir(tmp.path())
        .unwrap()
        .filter_map(Result::ok)
        .filter_map(|e| e.file_name().to_str().map(str::to_owned))
        .collect();
    for tempfile_name in ["zipped.bam", "after_sort.bam", "after_group.bam"] {
        assert!(
            !leftover_names.iter().any(|n| n == tempfile_name),
            "fused-pipeline tempfile {tempfile_name} leaked: {leftover_names:?}"
        );
    }
}

/// Regression gate for `--start-from zipper --stop-after duplex`.
///
/// Uses `zipper_duplex_fixture` (paired-UMI strand-A/strand-B records
/// in both mapped + unmapped BAMs) so a correctly-wired fused chain
/// produces ≥ 1 consensus read. The previous fixture
/// (`zipper_fixture`, 4 unpaired single-end records) couldn't surface
/// the wire-up bug because both buggy and correct chains produced 0
/// records on it (no `/A`+`/B` partners after paired grouping).
///
/// The wire-up bug being guarded: `execute_duplex_pipeline` was
/// missing the `else if let Some(zipper) = zipper` arm so
/// `Some(zipper)` silently fell through to the no-sort baseline
/// branch, processing the mapped BAM alone (no RX tags →
/// everything filtered → header-only BAM, which still satisfies
/// `meta.len() > 0`). The record-count assertion below makes that
/// regression visible as a test failure rather than a silent zero.
#[cfg(feature = "consensus")]
#[test]
fn zipper_to_duplex_fused_pipeline() {
    let tmp = TempDir::new().unwrap();
    let inputs = zipper_duplex_fixture(tmp.path());
    let out = tmp.path().join("zipper-duplex.bam");
    let args = ParityArgs::for_duplex();
    let r = run_runall_zipper(Stage::Duplex, &inputs, &out, &args);
    assert!(
        r.status.success(),
        "Zipper→Duplex fused failed: stdout={} stderr={}",
        String::from_utf8_lossy(&r.stdout),
        String::from_utf8_lossy(&r.stderr)
    );
    assert!(fs::metadata(&out).expect("duplex BAM").len() > 0);

    let mut reader = bam::io::reader::Builder.build_from_path(&out).expect("open duplex BAM");
    let _ = reader.read_header().expect("read duplex header");
    let n_records = reader.record_bufs(&noodles::sam::Header::default()).count();
    assert!(
        n_records >= 1,
        "zipper→duplex fused chain produced 0 consensus reads — likely a wire-up \
         regression (`execute_duplex_pipeline`'s `else if let Some(zipper) = zipper` \
         arm is missing or unreachable, so the chain fell through to the no-sort \
         baseline and processed the mapped BAM without RX tags)"
    );
}

/// Regression gate for `--start-from zipper --stop-after codec`.
///
/// Uses `zipper_codec_fixture` which mirrors the FR-overlap shape
/// `tests/integration/test_codec_command.rs::create_codec_read_pair`
/// uses — same flags (`MATE_REVERSE`/`REVERSE`), same R1+R2
/// position, single-strand UMI, 3 pairs per family. Standalone
/// codec is known to emit ≥ 1 consensus on this shape, so the
/// fused chain should too.
///
/// Same wire-up bug as `zipper_to_duplex_fused_pipeline` guards —
/// `execute_codec_pipeline` was previously missing its
/// `else if let Some(zipper) = zipper` arm.
#[cfg(feature = "consensus")]
#[test]
fn zipper_to_codec_fused_pipeline() {
    let tmp = TempDir::new().unwrap();
    let inputs = zipper_codec_fixture(tmp.path());
    let out = tmp.path().join("zipper-codec.bam");
    let args = ParityArgs::for_simplex_like();
    let r = run_runall_zipper(Stage::Codec, &inputs, &out, &args);
    assert!(
        r.status.success(),
        "Zipper→Codec fused failed: stdout={} stderr={}",
        String::from_utf8_lossy(&r.stdout),
        String::from_utf8_lossy(&r.stderr)
    );
    assert!(fs::metadata(&out).expect("codec BAM").len() > 0);

    let mut reader = bam::io::reader::Builder.build_from_path(&out).expect("open codec BAM");
    let _ = reader.read_header().expect("read codec header");
    let n_records = reader.record_bufs(&noodles::sam::Header::default()).count();
    assert!(
        n_records >= 1,
        "zipper→codec fused chain produced 0 consensus reads — likely a wire-up \
         regression (`execute_codec_pipeline`'s `else if let Some(zipper) = zipper` \
         arm is missing or unreachable)"
    );
}

/// Regression gate for the fused `--start-from align --stop-after
/// group` chain: AAM → Sort → Group → Write runs in a single
/// `Pipeline::run`. The new tap point — `MiGroupsToRecordBatch`
/// after `TemplatesToMiGroups` — is what makes group output land
/// directly in the BGZF compress step (no `aam-merged.bam` /
/// `sorted.bam` bridges).
///
/// Note: there's no explicit assertion that the v1 bridge files
/// don't exist anywhere — the v1 code created them in
/// `tempfile::tempdir()` (system tmpdir, not the test's
/// `TempDir`), so they would never appear in `tmp.path()` even on
/// the v1 code. The fusion guarantee comes from chain assembly
/// (no intermediate `Write` step exists in the typed-step
/// builder), enforced by `Pipeline::build()`.
#[test]
fn aam_to_group_fused_pipeline_with_mock_aligner() {
    let tmp = TempDir::new().unwrap();
    let inputs = zipper_fixture(tmp.path());
    let out = tmp.path().join("aam-group.bam");
    let mock = write_mock_aligner_script(tmp.path());

    let mut args = aam_base_args(&inputs.unmapped, &inputs.reference, &out);
    let stop_idx = args
        .iter()
        .position(|a| a == std::ffi::OsStr::new("--stop-after"))
        .expect("--stop-after flag in aam_base_args");
    args[stop_idx + 1] = "group".into();
    args.push("--sort::max-memory".into());
    args.push("256M".into());
    args.push("--group::strategy".into());
    args.push("adjacency".into());
    args.push("--aligner::command".into());
    args.push(format!("{} {{ref}}", mock.display()).into());

    let r = fgumi_with_args(&args);
    assert!(
        r.status.success(),
        "AAM→Group fused (mock aligner) failed: stdout={} stderr={}",
        String::from_utf8_lossy(&r.stdout),
        String::from_utf8_lossy(&r.stderr)
    );
    let meta = fs::metadata(&out).expect("group BAM must exist");
    assert!(meta.len() > 0, "group BAM was empty");
    // Smoke-test the chain end-to-end: the output BAM must open
    // and parse a header without erroring out. Record counts vary
    // by mock-aligner / fixture pairing — single-end records from
    // the mock aligner don't form valid templates for the default
    // group filter, so the record count assertion is left to
    // higher-fidelity tests in `parity_b_sort_to_group` (which
    // uses a paired-end fixture).
    let mut reader = bam::io::reader::Builder.build_from_path(&out).expect("open group BAM");
    let _ = reader.read_header().expect("read group header");
}

/// Regression gate for the fused `--start-from zipper --stop-after
/// sort` chain: Zipper → Sort → Write runs in a single
/// `Pipeline::run` (no `zipped.bam` bridge).
#[test]
fn zipper_to_sort_fused_pipeline() {
    let tmp = TempDir::new().unwrap();
    let inputs = zipper_fixture(tmp.path());
    let out = tmp.path().join("zipper-sorted.bam");
    let args = ParityArgs::for_simplex_like();
    let r = run_runall_zipper(Stage::Sort, &inputs, &out, &args);
    assert!(
        r.status.success(),
        "Zipper→Sort fused failed: stdout={} stderr={}",
        String::from_utf8_lossy(&r.stdout),
        String::from_utf8_lossy(&r.stderr)
    );
    assert!(fs::metadata(&out).expect("sorted BAM").len() > 0);
    // Sort doesn't filter records, so 4 reads in → 4 reads out;
    // pin that as a strong sanity check (vs. just header-only
    // BAM survival).
    let mut reader = bam::io::reader::Builder.build_from_path(&out).expect("open sorted BAM");
    let _ = reader.read_header().expect("read sorted header");
    let n_records = reader.record_bufs(&noodles::sam::Header::default()).count();
    assert_eq!(n_records, 4, "expected exactly 4 records in sorted BAM, got {n_records}");
}

/// Regression gate for the fused `--start-from zipper --stop-after
/// group` chain: Zipper → Sort → Group → Write runs in a single
/// `Pipeline::run` (no `zipped.bam` / `after_sort.bam` bridges).
#[test]
fn zipper_to_group_fused_pipeline() {
    let tmp = TempDir::new().unwrap();
    let inputs = zipper_fixture(tmp.path());
    let out = tmp.path().join("zipper-group.bam");
    let args = ParityArgs::for_simplex_like();
    let r = run_runall_zipper(Stage::Group, &inputs, &out, &args);
    assert!(
        r.status.success(),
        "Zipper→Group fused failed: stdout={} stderr={}",
        String::from_utf8_lossy(&r.stdout),
        String::from_utf8_lossy(&r.stderr)
    );
    assert!(fs::metadata(&out).expect("group BAM").len() > 0);
    // See `aam_to_group_fused_pipeline_with_mock_aligner` for why
    // the record count is not strict-checked here (the zipper
    // fixture is also single-end mocked).
    let mut reader = bam::io::reader::Builder.build_from_path(&out).expect("open group BAM");
    let _ = reader.read_header().expect("read group header");
}

/// `--rejects` is consensus-only; on non-consensus
/// `--stop-after` (`group` here) the flag must be ignored without
/// truncate-creating a 0-byte BAM at the path. Regression gate
/// for the silently-ignore behavior `build_pipeline_context`
/// gates via `track_rejects = is_enabled() && stop_after.is_consensus()`.
#[test]
fn rejects_out_with_stop_after_group_does_not_create_file() {
    let tmp = TempDir::new().unwrap();
    let fixture = unsorted_simplex_fixture(tmp.path());
    let rejects = tmp.path().join("rejects.bam");
    let out = tmp.path().join("group.bam");

    let cli_args: Vec<std::ffi::OsString> = vec![
        "runall".into(),
        "--start-from".into(),
        "sort".into(),
        "--stop-after".into(),
        "group".into(),
        "--input".into(),
        fixture.as_os_str().to_owned(),
        "--output".into(),
        out.as_os_str().to_owned(),
        "--rejects".into(),
        rejects.as_os_str().to_owned(),
        "--group::strategy".into(),
        "adjacency".into(),
        "--threads".into(),
        "2".into(),
    ];
    let r = fgumi_with_args(&cli_args);
    assert!(
        r.status.success(),
        "stop-after group with --rejects should succeed: stdout={} stderr={}",
        String::from_utf8_lossy(&r.stdout),
        String::from_utf8_lossy(&r.stderr)
    );
    assert!(fs::metadata(&out).expect("group BAM").len() > 0);
    // Critical: the rejects file must NOT exist — the gate
    // suppresses `build_rejects_writer` for non-consensus stops.
    // A 0-byte file at `rejects` would mean the writer was
    // opened but never finalized (the regression this test
    // guards against).
    assert!(
        !rejects.exists(),
        "rejects file was created for non-consensus --stop-after group; \
         the rejects writer should have been suppressed via the \
         is_consensus() gate"
    );
}

/// When the aligner subprocess exits 0 with no output at all
/// (e.g. user passes `--aligner::command "true"`), AAM detects the
/// empty stdout up-front and bails with a clear "aligner produced
/// no output" message rather than letting Zipper fail later with
/// a cryptic "invalid SAM header".
#[test]
fn aam_rejects_empty_aligner_output() {
    let tmp = TempDir::new().unwrap();
    let inputs = zipper_fixture(tmp.path());
    let out = tmp.path().join("out.bam");
    let mut args = aam_base_args(&inputs.unmapped, &inputs.reference, &out);
    // `cat > /dev/null` reads stdin to EOF and writes nothing —
    // simulating an aligner that gracefully consumed input but
    // emitted zero bytes.
    args.push("--aligner::command".into());
    args.push("cat > /dev/null # {ref}".into());
    let r = fgumi_with_args(&args);
    assert!(!r.status.success(), "empty aligner output should fail");
    let stderr = String::from_utf8_lossy(&r.stderr);
    // The AAM-as-Step pipeline surfaces the empty-output case via
    // the stable `align_and_merge::ERR_ALIGNER_EXITED_BEFORE_OUTPUT`
    // marker. The legacy orchestrator phrased it as "aligner produced
    // no output"; we accept either so the test isn't churned by
    // rewordings of the user-facing message.
    assert!(
        stderr.contains("exited before emitting any output")
            || stderr.contains("aligner produced no output"),
        "expected explicit empty-output error, got: {stderr}"
    );
}

/// When the aligner subprocess produces nonsense output (here:
/// `echo` instead of a real aligner), the AAM chain fails cleanly
/// at the zipper-merge step. Pins that errors propagate as a clean
/// non-zero exit (NOT a panic / unreachable).
#[test]
fn aam_fails_cleanly_on_bogus_aligner_output() {
    let tmp = TempDir::new().unwrap();
    let inputs = zipper_fixture(tmp.path());
    let out = tmp.path().join("out.bam");
    let mut args = aam_base_args(&inputs.unmapped, &inputs.reference, &out);
    // `cat > /dev/null` drains the writer's stdin to EOF *before* the bogus
    // line is emitted (same guard the sibling `aam_rejects_empty_aligner_output`
    // uses). A bare `echo dummy {ref}` exits immediately and closes its stdin,
    // which races the AAM writer's flush: under slower runs (e.g. llvm-cov
    // coverage) the flush loses the race and surfaces a "Broken pipe" I/O error
    // instead of the deterministic `@SQ count` mismatch this test pins. Draining
    // stdin first removes the race while leaving stdout (the bogus, non-SAM
    // `dummy {ref}` line ⇒ 0 `@SQ` records) unchanged.
    args.push("--aligner::command".into());
    args.push("cat > /dev/null; echo dummy {ref}".into());
    let r = fgumi_with_args(&args);
    assert!(!r.status.success(), "AAM should fail on bogus aligner output");
    let stderr = String::from_utf8_lossy(&r.stderr);
    assert!(!stderr.contains("panicked"), "unexpected panic: {stderr}");
    assert!(!stderr.contains("unreachable"), "unexpected unreachable: {stderr}");
    // S9a-006: a clean non-zero exit is necessary but not sufficient — pin
    // that the failure is the EXPECTED AAM-step error, not an incidental clap
    // / I/O failure before the aligner is reached. `echo dummy {ref}` emits a
    // non-SAM line, so the aligner's parsed header carries 0 `@SQ` records,
    // which the AAM `@SQ`-consistency check rejects against the reference
    // dict's 1 `@SQ`. The `"AAM"` step marker plus the stable `"@SQ count"`
    // fragment confirm the failure happens at the aligner-merge step.
    assert!(
        stderr.contains("AAM") && stderr.contains("@SQ count"),
        "expected the stable AAM @SQ-mismatch error at the merge step, got: {stderr}"
    );
}

/// `--methylation-mode em-seq` paired with `--start-from align-and-merge`
/// is rejected (methylation auto-routing through presets is deferred
/// to a follow-up PR; EM-seq users go through `--aligner::command`).
#[test]
fn aam_rejects_methylation_mode_paired_with_aam() {
    let tmp = TempDir::new().unwrap();
    let inputs = zipper_fixture(tmp.path());
    let out = tmp.path().join("out.bam");
    let mut args = aam_base_args(&inputs.unmapped, &inputs.reference, &out);
    args.push("--methylation-mode".into());
    args.push("em-seq".into());
    args.push("--aligner::command".into());
    args.push("echo dummy {ref}".into());
    let r = fgumi_with_args(&args);
    assert!(!r.status.success(), "methylation-mode + AAM should fail");
    let stderr = String::from_utf8_lossy(&r.stderr);
    assert!(
        stderr.contains("--methylation-mode is not yet supported"),
        "expected methylation rejection, got: {stderr}"
    );
}

/// `--stop-after align` is rejected at the clap parse layer: `align`
/// is absent from the `StopAfter` value-enum (no useful stop point
/// between alignment and zipper-merge — raw aligner output without
/// merge loses every original tag), so clap reports an invalid value
/// before any runtime validation runs.
#[test]
fn aam_rejects_stop_after_align_and_merge() {
    let tmp = TempDir::new().unwrap();
    let inputs = zipper_fixture(tmp.path());
    let out = tmp.path().join("out.bam");
    let args: Vec<std::ffi::OsString> = vec![
        "runall".into(),
        "--start-from".into(),
        "align".into(),
        "--stop-after".into(),
        "align".into(),
        "--input".into(),
        inputs.mapped.as_os_str().to_owned(),
        "--output".into(),
        out.as_os_str().to_owned(),
        "--ref".into(),
        inputs.reference.as_os_str().to_owned(),
        "--threads".into(),
        "1".into(),
        "--aligner::command".into(),
        "echo dummy {ref}".into(),
    ];
    let r = fgumi_with_args(&args);
    assert!(!r.status.success(), "--stop-after align should fail");
    let stderr = String::from_utf8_lossy(&r.stderr);
    assert!(
        stderr.contains("invalid value 'align'") && stderr.contains("--stop-after"),
        "expected clap parse-layer rejection of --stop-after align, got: {stderr}"
    );
}

// ────────────────────────── AAM real-aligner parity tests ──────────────────────────
//
// These exercise the AAM chain end-to-end against a real aligner
// binary (bwa-mem3 / bwa). They are gated by `which::which()` so they
// skip cleanly when the aligner isn't installed (most local dev
// environments); CI installs both binaries from bioconda before
// running the test suite (see `.github/workflows/check.yml`).
//
// The parity assertion is: runall AAM via `--aligner::preset <X>` must
// produce the same record stream as `fgumi fastq | <X> mem -p ... |
// fgumi zipper`. Same input, same aligner argv, same merge function —
// equivalence is by construction.

use crate::helpers::bam_generator::build_aligner_index;

/// Build a small reference + run the supplied aligner's `index` command on it.
/// Returns `(ref_path, _temp_dir)`; `_temp_dir` extends the lifetime of the
/// indexed files.
///
/// S9a-002: this is called ONLY after the caller has confirmed the binary is
/// on `PATH`, so EVERY failure here (tempdir creation, `<binary> index`
/// non-zero exit) is a real setup regression and `panic!`s loudly rather than
/// being laundered into a silent skip. A skip that is indistinguishable from a
/// pass would mask an index-build regression; the only legitimate skip
/// condition — a missing binary — is handled by the caller's `which` gate
/// before this runs.
fn build_aam_reference_for(binary_name: &str) -> (PathBuf, TempDir) {
    let tmp = TempDir::new().expect("AAM parity: create reference tempdir");
    let ref_path = crate::helpers::bam_generator::create_test_reference(tmp.path());
    build_aligner_index(&ref_path, binary_name).unwrap_or_else(|e| {
        panic!(
            "AAM parity: `{binary_name} index` failed with the binary present (real setup \
             regression, not a skip): {e}"
        )
    });
    (ref_path, tmp)
}

/// Run an AAM parity test: invoke runall with the supplied
/// `aligner_args` and the equivalent shell-pipeline staged chain,
/// then assert the merged BAMs are record-equivalent.
///
/// Three CRITICAL contracts the staged pipeline mirrors so parity
/// holds:
///   * `fgumi fastq -n` — the `-n` flag suppresses `/1` / `/2`
///     suffixes on FASTQ names. Runall AAM's internal writer uses
///     `no_suffix = true`; we MUST match that on the shell side
///     or the SAM QNAME-comparison in Zipper will diverge (and
///     bwa's silent suffix-stripping behaviour is
///     version-dependent — don't rely on it).
///   * `-K {DEFAULT_ALIGNER_CHUNK_SIZE}` — sourced from the same
///     constant runall uses, so a future tuning of the default
///     won't drift the test.
///   * `-t 1` — pin single-threaded aligner output for stable
///     record ordering.
fn run_aam_parity_test(aligner_binary: &str, aligner_args: &[(&str, &str)]) {
    // S9a-002: a MISSING aligner binary is the ONLY condition that skips —
    // legitimate for a binary-less dev laptop running with `--include-ignored`.
    // Any OTHER setup failure (tempdir, index build with the binary present) is
    // a real regression and panics in `build_aam_reference_for` rather than
    // green-passing with zero assertions. (CI installs both binaries and
    // verifies them on PATH before these `#[ignore]`'d tests run, so this skip
    // is unreachable in CI.)
    if which::which(aligner_binary).is_err() {
        eprintln!(
            "skipping AAM parity test: `{aligner_binary}` not found on PATH \
             (install via bioconda for CI)"
        );
        return;
    }
    let (ref_path, _ref_temp) = build_aam_reference_for(aligner_binary);

    let tmp = TempDir::new().unwrap();
    // The fixture uses single-end UNMAPPED reads with RX tags
    // (no PAIRED/FIRST/LAST flags). bwa `-p` will pair
    // consecutive records positionally regardless of flags;
    // determinism holds for a fixed `-t 1` bwa binary even
    // though MAPQ on the ambiguous repeated-reference fixture
    // is 0. A more realistic AAM fixture is a follow-up
    // (#33 task: dedicated AAM input fixture with paired flags
    // + diverse reference); the parity test today is
    // structural ("runall reproduces the shell pipeline") not
    // qualitative.
    let inputs = zipper_fixture(tmp.path());
    let runall_out = tmp.path().join("runall-aam.bam");
    let staged_out = tmp.path().join("staged-aam.bam");

    // Side 1: runall AAM.
    let mut args = aam_base_args(&inputs.unmapped, &ref_path, &runall_out);
    for (flag, value) in aligner_args {
        args.push((*flag).into());
        args.push((*value).into());
    }
    let r = fgumi_with_args(&args);
    assert!(
        r.status.success(),
        "runall AAM ({aligner_binary}) failed: stderr={}",
        String::from_utf8_lossy(&r.stderr)
    );

    // Side 2: shell-pipeline staged chain. `-n` on `fgumi fastq`
    // matches runall AAM's no-suffix FASTQ writer. `-K` is
    // sourced from the same const runall uses so future drift
    // can't desync the two sides.
    let staged_status = std::process::Command::new("/bin/bash")
        .args([
            "-c",
            &format!(
                "{bin} fastq -n -i {input} | \
                 {aligner_binary} mem -p -K {chunk} -t 1 {ref} /dev/stdin | \
                 {bin} zipper -i - -u {input} -r {ref} -o {out}",
                bin = fgumi_binary().display(),
                input = inputs.unmapped.display(),
                ref = ref_path.display(),
                out = staged_out.display(),
                chunk = fgumi_lib::aligner::DEFAULT_ALIGNER_CHUNK_SIZE,
            ),
        ])
        .status()
        .expect("spawn staged shell pipeline");
    assert!(staged_status.success(), "staged shell pipeline failed");

    assert_bams_record_equivalent(&runall_out, &staged_out);
}

/// Parity test: `runall --start-from align-and-merge --aligner::preset
/// bwa-mem3` vs the equivalent shell-pipeline staged chain.
/// `#[ignore]`'d because it needs bwa-mem3 on PATH (CI installs it
/// from bioconda; run locally with `cargo nextest run --include-ignored
/// aam_parity` once `bwa-mem3` is on PATH).
#[test]
#[ignore = "requires bwa-mem3 on PATH (install via bioconda)"]
fn aam_parity_preset_bwa_mem3() {
    run_aam_parity_test(
        "bwa-mem3",
        &[("--aligner::preset", "bwa-mem3"), ("--aligner::threads", "1")],
    );
}

/// Parity test: `runall --start-from align-and-merge --aligner::preset bwa`
/// vs the equivalent shell-pipeline staged chain. `#[ignore]`'d
/// because it needs bwa on PATH.
#[test]
#[ignore = "requires bwa on PATH (install via bioconda)"]
fn aam_parity_preset_bwa() {
    run_aam_parity_test("bwa", &[("--aligner::preset", "bwa"), ("--aligner::threads", "1")]);
}

/// Parity test: `runall --start-from align-and-merge --aligner::command
/// "bwa-mem3 mem ..."` (command mode) vs the equivalent shell-pipeline.
/// Pins that command mode is aligner-agnostic at the user surface.
/// `#[ignore]`'d because it needs bwa-mem3 on PATH.
#[test]
#[ignore = "requires bwa-mem3 on PATH (install via bioconda)"]
fn aam_parity_command_mode_bwa_mem3() {
    // Build the command template using the same chunk_size const
    // the preset mode uses, so command-mode parity isn't driven
    // by a different `-K` value than preset mode.
    let chunk = fgumi_lib::aligner::DEFAULT_ALIGNER_CHUNK_SIZE;
    let cmd = format!("bwa-mem3 mem -p -K {chunk} -t 1 {{ref}} /dev/stdin");
    run_aam_parity_test("bwa-mem3", &[("--aligner::command", &cmd)]);
}

// ────────────────────────── Correct parity tests (EC-C3) ──────────────────────────

/// Build an unsorted-UMI fixture specifically for the correct
/// parity test: four UMI families whose UMIs are deliberately
/// off-whitelist by 0, 1, 2, and 3 mismatches respectively against
/// the canonical whitelist UMI `ACGTACGT`. With
/// `--max-mismatches 2 --min-distance 1` the corrector accepts the
/// 0/1/2-mismatch families (overwriting RX, stashing the original
/// in OX) and rejects the 3-mismatch one. This exercises every
/// branch of the correction code path (no-op, single rewrite,
/// rejection) under both invocations, so parity must hold across
/// `RX` rewrite, `OX` stash, rejected-read emission, and metrics
/// counter increments — not just trivial pass-through.
fn correct_parity_fixture(dir: &Path) -> PathBuf {
    let path = dir.join("correct_parity_fixture.bam");
    let header = create_minimal_header("chr1", 10_000);
    let mut writer =
        bam::io::Writer::new(fs::File::create(&path).expect("create correct-parity fixture"));
    writer.write_header(&header).expect("write correct-parity header");

    // (umi, family_name, sequence) — same length as ACGTACGT (8) for
    // all variants. ACGTACGT exact match; ACGTACGC = 1 mismatch;
    // ACGTACCC = 2 mismatches; ACGTCCCC = 3 mismatches (exceeds
    // --max-mismatches 2, rejected).
    let families = [
        ("ACGTACGT", "fam_exact", "ACGTACGTACGTACGT"),
        ("ACGTACGC", "fam_1mm", "TTTTAAAACCCCGGGG"),
        ("ACGTACCC", "fam_2mm", "GGGGTTTTAAAACCCC"),
        ("ACGTCCCC", "fam_4mm_rejected", "CCCCAAAATTTTGGGG"),
    ];
    let mut all_records = Vec::new();
    for (umi, name, seq) in families {
        for r in create_umi_family(umi, 2, name, seq, 30) {
            all_records.push(r);
        }
    }
    for raw in &all_records {
        writer
            .write_alignment_record(&header, &to_record_buf(raw))
            .expect("write correct-parity record");
    }
    writer.try_finish().expect("finish correct-parity fixture");
    path
}

/// Class A — `correct → correct` parity vs standalone `fgumi correct`.
/// Pins the runall delegation: runall constructs a `CorrectUmis`
/// struct from its own flags and calls `CorrectUmis::execute`, so the
/// output is record- and header-equivalent (ignoring `@PG` command-line
/// provenance) by construction. The fixture
/// (`correct_parity_fixture`) is built specifically to exercise
/// every branch of the correction code path — no-op, single-
/// mismatch rewrite (RX overwritten + OX stashed), and rejection
/// (4-mismatch family fails) — under both invocations, plus
/// `--metrics` so the rejects-BAM and metrics-TSV outputs are also
/// parity-checked. This would catch a delegation bug that dropped
/// any of `rejects_opts`, `--correct::max-mismatches`,
/// `--correct::min-distance`, or the `OX`-stash behavior on a
/// successful correction.
#[test]
fn parity_a_correct_to_correct() {
    let tmp = TempDir::new().unwrap();
    let fixture = correct_parity_fixture(tmp.path());
    let whitelist = tmp.path().join("whitelist.txt");
    fs::write(&whitelist, "ACGTACGT\n").expect("write whitelist");

    let runall_out = tmp.path().join("runall.bam");
    let runall_rejects = tmp.path().join("runall.rejects.bam");
    let runall_metrics = tmp.path().join("runall.metrics.tsv");
    let standalone_out = tmp.path().join("standalone.bam");
    let standalone_rejects = tmp.path().join("standalone.rejects.bam");
    let standalone_metrics = tmp.path().join("standalone.metrics.tsv");

    // runall
    let runall_args: Vec<OsString> = vec![
        "runall".into(),
        "--start-from".into(),
        "correct".into(),
        "--stop-after".into(),
        "correct".into(),
        "--input".into(),
        fixture.as_os_str().to_owned(),
        "--output".into(),
        runall_out.as_os_str().to_owned(),
        "--rejects".into(),
        runall_rejects.as_os_str().to_owned(),
        "--correct::umi-files".into(),
        whitelist.as_os_str().to_owned(),
        "--correct::min-distance".into(),
        "1".into(),
        "--correct::max-mismatches".into(),
        "2".into(),
        "--correct::metrics".into(),
        runall_metrics.as_os_str().to_owned(),
        "--threads".into(),
        "1".into(),
    ];
    let r = fgumi(&runall_args);
    assert!(r.status.success(), "runall: {}", String::from_utf8_lossy(&r.stderr));

    // standalone
    let standalone_args: Vec<OsString> = vec![
        "correct".into(),
        "--input".into(),
        fixture.as_os_str().to_owned(),
        "--output".into(),
        standalone_out.as_os_str().to_owned(),
        "--rejects".into(),
        standalone_rejects.as_os_str().to_owned(),
        "--umi-files".into(),
        whitelist.as_os_str().to_owned(),
        "--min-distance".into(),
        "1".into(),
        "--max-mismatches".into(),
        "2".into(),
        "--metrics".into(),
        standalone_metrics.as_os_str().to_owned(),
        "--threads".into(),
        "1".into(),
    ];
    let r = fgumi(&standalone_args);
    assert!(r.status.success(), "standalone: {}", String::from_utf8_lossy(&r.stderr));

    // The fixture is designed to emit both kept reads and rejects, so pin both
    // streams non-empty — otherwise a shared-empty regression would silently
    // satisfy the record/header equivalence checks on both paths.
    assert_bams_record_equivalent_nonempty(&runall_out, &standalone_out);
    assert_bams_record_equivalent_nonempty(&runall_rejects, &standalone_rejects);
    // S9a-005: both the accepted-output and rejects BAMs must carry equivalent
    // @HD/@SQ/@RG headers between the fused runall correct stage and standalone.
    assert_bam_headers_equivalent_ignoring_pg(&runall_out, &standalone_out);
    assert_bam_headers_equivalent_ignoring_pg(&runall_rejects, &standalone_rejects);
    // Metrics TSVs are plain text; compare exactly (both are produced
    // by the same `UmiCorrectionMetrics::write_metrics` writer, so any
    // forwarding bug would surface as different counts on one side).
    let runall_metrics_text = fs::read_to_string(&runall_metrics).expect("read runall metrics tsv");
    let standalone_metrics_text =
        fs::read_to_string(&standalone_metrics).expect("read standalone metrics tsv");
    assert_eq!(runall_metrics_text, standalone_metrics_text);
}

/// CLI rejection — `--start-from correct` without either
/// `--correct::umis` or `--correct::umi-files` should fail with a
/// clear message naming both flags.
#[test]
fn rejects_correct_without_umis() {
    let tmp = TempDir::new().unwrap();
    let fixture = unsorted_simplex_fixture(tmp.path());
    let out = tmp.path().join("out.bam");
    let args: Vec<OsString> = vec![
        "runall".into(),
        "--start-from".into(),
        "correct".into(),
        "--stop-after".into(),
        "correct".into(),
        "--input".into(),
        fixture.as_os_str().to_owned(),
        "--output".into(),
        out.as_os_str().to_owned(),
        "--correct::min-distance".into(),
        "1".into(),
        // (intentionally omit --correct::umis / --correct::umi-files)
        "--threads".into(),
        "1".into(),
    ];
    let r = fgumi(&args);
    assert!(!r.status.success(), "expected correct-without-umis to fail");
    let stderr = String::from_utf8_lossy(&r.stderr);
    assert!(
        stderr.contains("UMI or UMI file"),
        "stderr should mention UMI requirement, got: {stderr}"
    );
}

/// CLI rejection — `--start-from correct --stop-after sort` (any
/// cross-stage stop past `correct`) requires `--ref` for the AAM
/// aligner. The dispatcher validates aligner config up-front (before
/// loading UMI files / running correct), so a missing `--ref` should
/// surface a clear aligner-reference error.
#[test]
fn rejects_correct_to_sort_chain_without_ref() {
    let tmp = TempDir::new().unwrap();
    let fixture = unsorted_simplex_fixture(tmp.path());
    let whitelist = tmp.path().join("whitelist.txt");
    fs::write(&whitelist, "ACGTACGT\n").expect("write whitelist");
    let out = tmp.path().join("out.bam");
    let args: Vec<OsString> = vec![
        "runall".into(),
        "--start-from".into(),
        "correct".into(),
        "--stop-after".into(),
        "sort".into(),
        "--input".into(),
        fixture.as_os_str().to_owned(),
        "--output".into(),
        out.as_os_str().to_owned(),
        "--correct::umi-files".into(),
        whitelist.as_os_str().to_owned(),
        "--correct::min-distance".into(),
        "1".into(),
        "--group::strategy".into(),
        "identity".into(),
        "--threads".into(),
        "1".into(),
    ];
    let r = fgumi(&args);
    assert!(!r.status.success(), "expected correct→sort without --ref to fail");
    let stderr = String::from_utf8_lossy(&r.stderr);
    assert!(stderr.contains("--ref"), "stderr should mention --ref requirement, got: {stderr}");
}

// ───────────────── correct → AAM → ... fused smoke tests ─────────────────
//
// These tests exercise the `--start-from correct --stop-after X`
// dispatcher (T2.3) end-to-end for every non-`correct` `X`. Each test
// builds a UMI-tagged BAM + a whitelist matching the fixture's UMIs,
// writes a mock aligner script, and asserts the chain succeeds and
// emits ≥ 1 record. They mirror the existing
// `aam_to_*_fused_pipeline_with_mock_aligner` family in shape.
//
// **Smoke-only, not byte-equivalence parity.** Full parity against the
// staged equivalent (`fgumi correct` → `fgumi runall --start-from
// align → ...`) would require three separate fgumi processes per test;
// per-stage parity is pinned by the existing AAM parity tests
// (`aam_to_*_fused_pipeline_with_mock_aligner` and the `parity_b_*`
// family). The per-command `parity_*_new_vs_legacy` unit tests were
// deleted in Phase 1 once each command's legacy path was removed.
// What these tests add is the dispatcher-level guarantee that the
// correct → AAM chain wires up at all and produces non-empty output
// for each consensus / non-consensus stop-after value.

/// Build a (umi-tagged BAM, reference, whitelist) triple for
/// single-end correct→AAM tests. Reuses
/// `unsorted_simplex_fixture`'s 4-family mapped BAM (UMIs
/// `ACGTACGT`, `TGCATGCA`, `CCAATTGG`, `GGTTAACC`) and writes a
/// matching whitelist so every record passes correct's filter.
///
/// **Input shape caveat**: this reuses `unsorted_simplex_fixture`,
/// which writes MAPPED records (pos, mapq, flags). Production
/// `--start-from correct` expects an UNMAPPED UMI-tagged BAM — correct's
/// output then flows to the aligner as a FASTQ source. The tests work
/// with the mapped fixture because (a) correct doesn't enforce input
/// shape and (b) the mock aligner script is flag-agnostic (reads
/// sequence/quality, ignores flags + pos). For production input shapes,
/// the AAM step would re-align these records, which is structurally
/// fine but redundant. The shape mismatch is acceptable for smoke
/// testing the dispatcher's wire-up; per-stage byte-equivalence is
/// covered by `aam_to_*_fused_pipeline_with_mock_aligner` (the
/// per-stage `parity_*_new_vs_legacy` unit tests were deleted in
/// Phase 1 once each command's legacy path was removed).
fn correct_to_aam_fixture(dir: &Path) -> (PathBuf, PathBuf, PathBuf) {
    let input = unsorted_simplex_fixture(dir);
    let reference = crate::helpers::bam_generator::create_test_reference(dir);
    let whitelist = dir.join("whitelist.txt");
    fs::write(&whitelist, "ACGTACGT\nTGCATGCA\nCCAATTGG\nGGTTAACC\n").expect("write whitelist");
    (input, reference, whitelist)
}

/// Build a paired-end (umi-tagged BAM, reference, whitelist) triple
/// for duplex / codec correct→AAM tests. Reuses
/// `unsorted_duplex_fixture`'s 4-family paired-UMI BAM and
/// whitelists both halves of each paired UMI (`AAAA`/`CCCC`,
/// `GGGG`/`TTTT`, `CCAA`/`TTGG`, `ATCG`/`GCTA`) so every record
/// passes correct's filter — paired-UMI correct keys off each half
/// independently.
///
/// **Input shape caveat**: this reuses `unsorted_duplex_fixture`,
/// which writes MAPPED records (pos, mapq, flags). Production
/// `--start-from correct` expects an UNMAPPED UMI-tagged BAM — correct's
/// output then flows to the aligner as a FASTQ source. The tests work
/// with the mapped fixture because (a) correct doesn't enforce input
/// shape and (b) the mock aligner script is flag-agnostic (reads
/// sequence/quality, ignores flags + pos). For production input shapes,
/// the AAM step would re-align these records, which is structurally
/// fine but redundant. The shape mismatch is acceptable for smoke
/// testing the dispatcher's wire-up; per-stage byte-equivalence is
/// covered by `aam_to_*_fused_pipeline_with_mock_aligner` (the
/// per-stage `parity_*_new_vs_legacy` unit tests were deleted in
/// Phase 1 once each command's legacy path was removed).
fn correct_to_aam_paired_fixture(dir: &Path) -> (PathBuf, PathBuf, PathBuf) {
    let input = unsorted_duplex_fixture(dir);
    let reference = crate::helpers::bam_generator::create_test_reference(dir);
    let whitelist = dir.join("whitelist.txt");
    fs::write(&whitelist, "AAAA\nCCCC\nGGGG\nTTTT\nCCAA\nTTGG\nATCG\nGCTA\n")
        .expect("write whitelist");
    (input, reference, whitelist)
}

/// Build the common args used by every `--start-from correct
/// --stop-after X` smoke test. Caller appends stage-specific flags
/// (`--sort::max-memory`, `--group::strategy`, `--min-reads`, etc.)
/// before invoking `fgumi_with_args`.
fn correct_to_aam_base_args(
    input: &Path,
    ref_path: &Path,
    whitelist: &Path,
    output: &Path,
    stop_after: &str,
    mock_aligner_path: &Path,
) -> Vec<OsString> {
    // The three consensus algorithms collapse to the `consensus` stage;
    // the algorithm is passed separately via `--consensus`.
    let is_consensus = matches!(stop_after, "simplex" | "duplex" | "codec");
    let stop_value = if is_consensus { "consensus" } else { stop_after };
    let mut args: Vec<OsString> = vec![
        "runall".into(),
        "--start-from".into(),
        "correct".into(),
        "--stop-after".into(),
        stop_value.into(),
        "--input".into(),
        input.as_os_str().to_owned(),
        "--output".into(),
        output.as_os_str().to_owned(),
        "--ref".into(),
        ref_path.as_os_str().to_owned(),
        "--correct::umi-files".into(),
        whitelist.as_os_str().to_owned(),
        "--correct::min-distance".into(),
        "1".into(),
        "--aligner::command".into(),
        format!("{} {{ref}}", mock_aligner_path.display()).into(),
        "--threads".into(),
        "2".into(),
    ];
    if is_consensus {
        args.push("--consensus".into());
        args.push(stop_after.into());
    }
    args
}

/// Smoke: `--start-from correct --stop-after zipper`. Confirms the
/// fused correct → AAM → write chain runs to completion and emits a
/// non-empty merged BAM.
#[test]
fn correct_to_zipper_fused_pipeline() {
    let tmp = TempDir::new().unwrap();
    let (input, reference, whitelist) = correct_to_aam_fixture(tmp.path());
    let out = tmp.path().join("correct-to-zipper.bam");
    let mock = write_mock_aligner_script(tmp.path());

    let args = correct_to_aam_base_args(&input, &reference, &whitelist, &out, "zipper", &mock);

    let r = fgumi_with_args(&args);
    assert!(
        r.status.success(),
        "correct→Zipper fused failed: stdout={} stderr={}",
        String::from_utf8_lossy(&r.stdout),
        String::from_utf8_lossy(&r.stderr)
    );
    let meta = fs::metadata(&out).expect("merged BAM must exist");
    assert!(meta.len() > 0, "merged BAM was empty");
    let mut reader = bam::io::reader::Builder.build_from_path(&out).expect("open merged BAM");
    let _ = reader.read_header().expect("read merged header");
    let n_records = reader.record_bufs(&noodles::sam::Header::default()).count();
    assert!(n_records >= 1, "expected ≥ 1 record in merged BAM, got {n_records}");
}

/// Smoke: `--start-from correct --stop-after sort`. Confirms the
/// fused correct → AAM → Sort → write chain runs to completion and
/// emits a non-empty sorted BAM.
#[test]
fn correct_to_sort_fused_pipeline() {
    let tmp = TempDir::new().unwrap();
    let (input, reference, whitelist) = correct_to_aam_fixture(tmp.path());
    let out = tmp.path().join("correct-to-sort.bam");
    let mock = write_mock_aligner_script(tmp.path());

    let mut args = correct_to_aam_base_args(&input, &reference, &whitelist, &out, "sort", &mock);
    args.push("--sort::max-memory".into());
    args.push("256M".into());

    let r = fgumi_with_args(&args);
    assert!(
        r.status.success(),
        "correct→Sort fused failed: stdout={} stderr={}",
        String::from_utf8_lossy(&r.stdout),
        String::from_utf8_lossy(&r.stderr)
    );
    let meta = fs::metadata(&out).expect("sorted BAM must exist");
    assert!(meta.len() > 0, "sorted BAM was empty");
    let mut reader = bam::io::reader::Builder.build_from_path(&out).expect("open sorted BAM");
    let _ = reader.read_header().expect("read sorted header");
    let n_records = reader.record_bufs(&noodles::sam::Header::default()).count();
    assert!(n_records >= 1, "expected ≥ 1 record in sorted BAM, got {n_records}");
}

/// Smoke: `--start-from correct --stop-after group`. Confirms the
/// fused correct → AAM → Sort → Group → write chain runs to
/// completion. Record-count assertion is header-only because
/// single-end records from the default mock aligner don't form valid
/// templates for the default group filter — same shape as
/// `aam_to_group_fused_pipeline_with_mock_aligner`.
#[test]
fn correct_to_group_fused_pipeline() {
    let tmp = TempDir::new().unwrap();
    let (input, reference, whitelist) = correct_to_aam_fixture(tmp.path());
    let out = tmp.path().join("correct-to-group.bam");
    let mock = write_mock_aligner_script(tmp.path());

    let mut args = correct_to_aam_base_args(&input, &reference, &whitelist, &out, "group", &mock);
    args.push("--sort::max-memory".into());
    args.push("256M".into());
    args.push("--group::strategy".into());
    args.push("adjacency".into());

    let r = fgumi_with_args(&args);
    assert!(
        r.status.success(),
        "correct→Group fused failed: stdout={} stderr={}",
        String::from_utf8_lossy(&r.stdout),
        String::from_utf8_lossy(&r.stderr)
    );
    let meta = fs::metadata(&out).expect("group BAM must exist");
    assert!(meta.len() > 0, "group BAM was empty");
    let mut reader = bam::io::reader::Builder.build_from_path(&out).expect("open group BAM");
    let _ = reader.read_header().expect("read group header");
}

/// Smoke: `--start-from correct --stop-after simplex`. Confirms the
/// fused correct → AAM → Sort → Group → Consensus → write chain
/// runs to completion. Same record-count caveat as
/// `aam_to_simplex_fused_pipeline_with_mock_aligner`: the default
/// single-end mock aligner emits `flag=0` records, and simplex
/// requires both R1 + R2 of a pair to emit consensus, so this
/// produces a header-only BAM by construction. Smoke-test asserts
/// the chain succeeds + the output BAM opens cleanly.
#[cfg(feature = "consensus")]
#[test]
fn correct_to_simplex_fused_pipeline() {
    let tmp = TempDir::new().unwrap();
    let (input, reference, whitelist) = correct_to_aam_fixture(tmp.path());
    let out = tmp.path().join("correct-to-simplex.bam");
    let mock = write_mock_aligner_script(tmp.path());

    let mut args = correct_to_aam_base_args(&input, &reference, &whitelist, &out, "simplex", &mock);
    args.push("--sort::max-memory".into());
    args.push("256M".into());
    args.push("--group::strategy".into());
    args.push("adjacency".into());
    args.push("--simplex::min-reads".into());
    args.push("1".into());

    let r = fgumi_with_args(&args);
    assert!(
        r.status.success(),
        "correct→Simplex fused failed: stdout={} stderr={}",
        String::from_utf8_lossy(&r.stdout),
        String::from_utf8_lossy(&r.stderr)
    );
    let meta = fs::metadata(&out).expect("simplex BAM must exist");
    assert!(meta.len() > 0, "simplex BAM was empty");
    let mut reader = bam::io::reader::Builder.build_from_path(&out).expect("open simplex BAM");
    let _ = reader.read_header().expect("read simplex header");
}

/// Smoke: `--start-from correct --stop-after duplex`. Uses the
/// paired-end UMI fixture + paired mock aligner so the chain has
/// well-formed input (proper-pair flags, both R1 and R2 with
/// matching qnames). Asserts ≥ 1 consensus record — same shape as
/// `aam_to_duplex_fused_pipeline_with_mock_aligner`.
#[cfg(feature = "consensus")]
#[test]
fn correct_to_duplex_fused_pipeline() {
    let tmp = TempDir::new().unwrap();
    let (input, reference, whitelist) = correct_to_aam_paired_fixture(tmp.path());
    let out = tmp.path().join("correct-to-duplex.bam");
    let mock = write_paired_mock_aligner_script(tmp.path());

    let mut args = correct_to_aam_base_args(&input, &reference, &whitelist, &out, "duplex", &mock);
    args.push("--sort::max-memory".into());
    args.push("256M".into());
    args.push("--group::strategy".into());
    args.push("paired".into());
    args.push("--duplex::min-reads".into());
    args.push("1,1,0".into());

    let r = fgumi_with_args(&args);
    assert!(
        r.status.success(),
        "correct→Duplex fused failed: stdout={} stderr={}",
        String::from_utf8_lossy(&r.stdout),
        String::from_utf8_lossy(&r.stderr)
    );
    let meta = fs::metadata(&out).expect("duplex BAM must exist");
    assert!(meta.len() > 0, "duplex BAM was empty");
    let mut reader = bam::io::reader::Builder.build_from_path(&out).expect("open duplex BAM");
    let _ = reader.read_header().expect("read duplex header");
    let n_records = reader.record_bufs(&noodles::sam::Header::default()).count();
    assert!(n_records >= 1, "expected ≥ 1 record in duplex BAM, got {n_records}");
}

/// Smoke: `--start-from correct --stop-after codec`. Uses the
/// paired-end UMI fixture + paired mock aligner (which emits the
/// FR-overlap shape codec requires). Asserts ≥ 1 consensus record —
/// same shape as `aam_to_codec_fused_pipeline_with_mock_aligner`.
#[cfg(feature = "consensus")]
#[test]
fn correct_to_codec_fused_pipeline() {
    let tmp = TempDir::new().unwrap();
    let (input, reference, whitelist) = correct_to_aam_paired_fixture(tmp.path());
    let out = tmp.path().join("correct-to-codec.bam");
    let mock = write_paired_mock_aligner_script(tmp.path());

    let mut args = correct_to_aam_base_args(&input, &reference, &whitelist, &out, "codec", &mock);
    args.push("--sort::max-memory".into());
    args.push("256M".into());
    args.push("--group::strategy".into());
    args.push("adjacency".into());
    args.push("--codec::min-reads".into());
    args.push("1".into());

    let r = fgumi_with_args(&args);
    assert!(
        r.status.success(),
        "correct→Codec fused failed: stdout={} stderr={}",
        String::from_utf8_lossy(&r.stdout),
        String::from_utf8_lossy(&r.stderr)
    );
    let meta = fs::metadata(&out).expect("codec BAM must exist");
    assert!(meta.len() > 0, "codec BAM was empty");
    let mut reader = bam::io::reader::Builder.build_from_path(&out).expect("open codec BAM");
    let _ = reader.read_header().expect("read codec header");
    let n_records = reader.record_bufs(&noodles::sam::Header::default()).count();
    assert!(n_records >= 1, "expected ≥ 1 record in codec BAM, got {n_records}");
}

/// X3-006: combining `--rejects` with a cross-stage `--start-from correct
/// --stop-after sort` chain does NOT capture correct's UMI rejects — the fused
/// chain uses `correct_step_kept_only`, which has no UMI-rejects branch, so the
/// rejected reads are silently DROPPED (a near-data-loss the dispatcher guards
/// only with a warn log).
///
/// The old version of this test used a whitelist containing EVERY family's UMI,
/// so 0 reads were rejected — it could not distinguish "rejects correctly
/// discarded" from "no rejects existed at all" (the vacuous-pass pattern
/// S9a-001 calls out). This version uses a PARTIAL whitelist that genuinely
/// rejects one family, then PROVES the discard is real and bounded:
///   (a) the warn log fires,
///   (b) the fused `--rejects` BAM is EMPTY (correct's rejects were discarded,
///       not routed here), and
///   (c) a standalone `fgumi correct --rejects` on the SAME input produces a
///       NON-empty rejects BAM — i.e. real reads ARE lost, justifying the warning.
#[test]
fn correct_to_sort_with_rejects_emits_warning() {
    let tmp = TempDir::new().unwrap();
    // `unsorted_simplex_fixture` carries families ACGTACGT/TGCATGCA/CCAATTGG/
    // GGTTAACC. Drop GGTTAACC from the whitelist so its 4 reads are rejected by
    // correct (the remaining whitelist UMIs are all far from GGTTAACC, so it
    // cannot be corrected to any of them within `--correct::min-distance 1`).
    let input = unsorted_simplex_fixture(tmp.path());
    let reference = crate::helpers::bam_generator::create_test_reference(tmp.path());
    let whitelist = tmp.path().join("partial_whitelist.txt");
    fs::write(&whitelist, "ACGTACGT\nTGCATGCA\nCCAATTGG\n").expect("write partial whitelist");

    let out = tmp.path().join("correct-to-sort-rejects.bam");
    let rejects = tmp.path().join("rejects.bam");
    let mock = write_mock_aligner_script(tmp.path());

    let mut args = correct_to_aam_base_args(&input, &reference, &whitelist, &out, "sort", &mock);
    args.push("--rejects".into());
    args.push(rejects.as_os_str().to_owned());
    args.push("--sort::max-memory".into());
    args.push("256M".into());

    // Pin `RUST_LOG=info` on the child so an ambient `RUST_LOG=warn`
    // (or stricter) on the developer's shell doesn't silently mask
    // the warn-level log line — the binary's default level (set in
    // `main.rs`) is info, but env overrides win.
    let r = std::process::Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .env("RUST_LOG", "info")
        .args(&args)
        .output()
        .expect("failed to execute fgumi");
    assert!(
        r.status.success(),
        "correct→sort with --rejects failed: stdout={} stderr={}",
        String::from_utf8_lossy(&r.stdout),
        String::from_utf8_lossy(&r.stderr)
    );

    // (a) The warning fires.
    let stderr = String::from_utf8_lossy(&r.stderr);
    assert!(
        stderr.contains("discards correct's UMI rejects"),
        "expected UMI-rejects-discarded warning, got stderr:\n{stderr}"
    );

    // (b) The fused `--rejects` capture is EMPTY — correct's rejects were
    // discarded by `correct_step_kept_only`, and sort produces no rejects of its
    // own. The file is either absent (sort emitted nothing to it) or a
    // header-only BAM; either way, ZERO reject reads were captured.
    let fused_reject_count = if rejects.exists() { read_bam_records(&rejects).len() } else { 0 };
    assert_eq!(
        fused_reject_count, 0,
        "fused correct→sort --rejects must capture 0 reads (correct's rejects are discarded), got {fused_reject_count}"
    );

    // (c) Standalone `fgumi correct --rejects` on the SAME input + whitelist
    // captures the rejected family — proving the discard above lost REAL data,
    // exactly what the warning claims. The standalone correct path uses the
    // 2-output `correct_step_with_rejects` branch.
    let standalone_kept = tmp.path().join("standalone-correct.bam");
    let standalone_rejects = tmp.path().join("standalone-correct-rejects.bam");
    let correct_args: Vec<OsString> = vec![
        "correct".into(),
        "--input".into(),
        input.as_os_str().to_owned(),
        "--output".into(),
        standalone_kept.as_os_str().to_owned(),
        "--umi-files".into(),
        whitelist.as_os_str().to_owned(),
        "--min-distance".into(),
        "1".into(),
        "--rejects".into(),
        standalone_rejects.as_os_str().to_owned(),
        "--threads".into(),
        "1".into(),
    ];
    let r = fgumi(&correct_args);
    assert!(
        r.status.success(),
        "standalone correct --rejects failed: {}",
        String::from_utf8_lossy(&r.stderr)
    );
    let standalone_reject_records = read_bam_records(&standalone_rejects);
    assert_eq!(
        standalone_reject_records.len(),
        4,
        "standalone correct --rejects must capture the ENTIRE off-whitelist family that the \
         fused chain discarded — expected exactly the 4 GGTTAACC reads, got {}. A non-4 count \
         means only part of the family was rejected (or the wrong reads were captured).",
        standalone_reject_records.len()
    );
    // Pin the captured family by IDENTITY, not just count: every rejected read
    // must be a `fam_d` read carrying the off-whitelist UMI `GGTTAACC`, and the
    // four must be exactly fam_d_0..fam_d_3 — proving the discarded data was the
    // specific off-whitelist family, not some unrelated four reads.
    {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value;
        let rx_tag = Tag::from(fgumi_lib::sam::SamTag::RX);
        let mut reject_qnames: Vec<Vec<u8>> = Vec::new();
        for rec in &standalone_reject_records {
            let Some(Value::String(rx)) = rec.data().get(&rx_tag) else {
                panic!("rejected record is missing its RX (UMI) tag");
            };
            let rx_bytes: &[u8] = rx.as_ref();
            assert_eq!(
                rx_bytes, b"GGTTAACC",
                "rejected record carries the wrong UMI — only the off-whitelist GGTTAACC \
                 family should be rejected"
            );
            let qname: &[u8] = rec.name().expect("rejected record has a read name").as_ref();
            reject_qnames.push(qname.to_vec());
        }
        reject_qnames.sort();
        let expected_qnames: Vec<Vec<u8>> =
            (0..4).map(|i| format!("fam_d_{i}").into_bytes()).collect();
        assert_eq!(
            reject_qnames, expected_qnames,
            "rejected reads must be exactly the four GGTTAACC family members fam_d_0..fam_d_3"
        );
    }
}

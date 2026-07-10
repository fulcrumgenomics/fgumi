//! Mutation / metamorphic soundness harness for `fgumi compare`.
//!
//! Every other `compare` test proves the comparator recognizes something is right;
//! this file proves it can't be fooled into saying so when something is wrong. The
//! technique (metamorphic testing) is: build a pair of inputs that are already
//! known to `MATCH`, apply exactly **one** targeted mutation to one side, and
//! assert the *expected* verdict changes (or, for an accepted-divergence case,
//! deliberately does not). A single mutation that should flip the verdict to
//! `DIFFER` but doesn't is a proven soundness hole in the oracle CI relies on to
//! catch real regressions — so this file is deliberately organized as a flat
//! catalog of small, single-purpose cases rather than a few broad ones, and each
//! case is tagged in [`MUTATION_CATALOG`] (see [`mutation_catalog_covers_every_compared_dimension`]) so a
//! future compared field/tag/row-dimension that never got a mutation case is
//! caught by a coverage assertion, not just missed silently.
//!
//! Two classes of mutation are exercised per the compare-hardening design spec's
//! §"Soundness":
//!
//! - **Must DIFFER** — a real divergence the oracle must catch (a false `MATCH`
//!   here is a soundness hole in the CI gate).
//! - **Must NOT DIFFER** — an accepted/tolerated divergence (a false `DIFFER`
//!   here is a regression: a spurious CI failure on legitimate output).
//!
//! Every engine landed by the compare-hardening effort is exercised directly
//! (not only through the `fgumi compare bams`/`metrics` CLI, though metrics goes
//! through its `Command` entry point since that's its only public surface):
//! [`positional_compare`] under all three [`ContentPredicate`] variants,
//! [`keyjoin_compare`] (`group`'s merge-join + MI bijection), [`sort_verify_compare`]
//! (`sort`'s independent-baseline order check), [`compare_headers`], and
//! `fgumi compare metrics`'s row key-join (via [`CompareMetrics`]).
//!
//! Where a fixture in `harden-compare-corpus/fixtures` already exercises an
//! equivalent adversarial case (e.g. the accepted `r2-ucc-02` per-cell-consensus
//! partitioning divergence), this harness does not re-derive it from real fgbio
//! output — those fixtures are evidence for a *different* effort (the BS4
//! adversarial corpus / fgbio-parity burn-down); this harness's job is proving
//! the comparator's own logic is sound on deliberately minimal, programmatic
//! inputs, so every case here is built from scratch via `SamBuilder`.

use clap::Parser;
use fgumi_lib::commands::command::Command as FgumiCommand;
use fgumi_lib::commands::compare::{
    CompareMetrics, CompareMismatch, ContentPredicate, compare_headers, keyjoin_compare,
    positional_compare, sort_verify_compare,
};
use fgumi_lib::sam::SamTag;
use fgumi_raw_bam::{RawRecord, SamBuilder, flags};
use noodles::sam::Header;
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::header::tag as hd_tag;
use noodles::sam::header::record::value::map::program::tag as pg_tag;
use noodles::sam::header::record::value::map::read_group::tag as rg_tag;
use noodles::sam::header::record::value::map::{
    Header as HeaderRecord, Program, ReadGroup, ReferenceSequence,
};
use rstest::rstest;
use std::num::NonZeroUsize;
use std::path::Path;
use tempfile::TempDir;

use crate::helpers::bam_generator::{
    create_coordinate_sorted_header, create_minimal_header, keyjoin_cfg, mi_record, write_bam,
};
use crate::helpers::write_tsv;

// ---------------------------------------------------------------------------
// Shared helpers
// ---------------------------------------------------------------------------

/// Packs a `(length, op)` pair into the BAM CIGAR op-code encoding
/// (`length << 4 | op`; `op` 0 = `M`, 4 = `S`).
fn cigar_op(len: u32, op: u32) -> u32 {
    (len << 4) | op
}

/// A two-`@SQ` header (`chr1`, `chr2`), needed for the RNEXT mutation case (which sets
/// a mate reference id that must resolve to a *different* declared reference).
fn two_ref_header() -> Header {
    let seq1 = Map::<ReferenceSequence>::new(NonZeroUsize::new(10000).expect("non-zero"));
    let seq2 = Map::<ReferenceSequence>::new(NonZeroUsize::new(10000).expect("non-zero"));
    Header::builder()
        .add_reference_sequence(bstr::BString::from("chr1"), seq1)
        .add_reference_sequence(bstr::BString::from("chr2"), seq2)
        .build()
}

// =============================================================================
// Section 1 — Positional engine (`positional_compare`), `ContentPredicate::Exact`
// =============================================================================
//
// `positional_compare` pairs record `i` of bam1 with record `i` of bam2 and checks
// `RecordKey` equality *before* content — see `engines/positional.rs`'s doc comment.
// This section proves every core SAM field and every compared tag class flips the
// verdict to DIFFER on a single mutation, that tag order does not, and that record
// add/drop/swap are all caught (the last one being the "swap ≡ drop+add" masking
// test the no-resync design exists to defeat).

/// A fully-populated mapped, paired baseline record: every core SAM field the
/// comparator inspects (per `FIELD_NAMES`: QNAME, FLAG, RNAME, POS, MAPQ, CIGAR,
/// RNEXT, PNEXT, TLEN, SEQ, QUAL) has a non-default value, so every core-field
/// mutation case below changes exactly one field relative to this baseline.
/// `customize` is applied last, so a mutation case only needs to override the one
/// field it targets.
fn core_record_variant(customize: impl FnOnce(&mut SamBuilder)) -> RawRecord {
    let mut b = SamBuilder::new();
    b.read_name(b"read1")
        .flags(flags::PAIRED | flags::FIRST_SEGMENT)
        .ref_id(0)
        .pos(99)
        .mapq(60)
        .cigar_ops(&[cigar_op(10, 0)]) // 10M
        .sequence(b"ACGTACGTAC")
        .qualities(&[30; 10])
        .mate_ref_id(0)
        .mate_pos(199)
        .template_length(150)
        .add_int_tag(SamTag::NM, 1)
        .add_string_tag(SamTag::RX, b"AAAA");
    customize(&mut b);
    b.build()
}

fn mutate_seq(b: &mut SamBuilder) {
    b.sequence(b"ACGTACGTAG").qualities(&[30; 10]);
}
fn mutate_qual(b: &mut SamBuilder) {
    b.qualities(&[30, 30, 30, 30, 30, 30, 30, 30, 30, 10]);
}
fn mutate_pos(b: &mut SamBuilder) {
    b.pos(199);
}
fn mutate_flag(b: &mut SamBuilder) {
    // Toggle DUPLICATE (0x400): must not touch FIRST/LAST/SECONDARY/SUPPLEMENTARY,
    // or the mutation would flip `RecordKey` and be caught as a key mismatch
    // instead of the intended content diff.
    b.flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::DUPLICATE);
}
fn mutate_cigar(b: &mut SamBuilder) {
    b.cigar_ops(&[cigar_op(9, 0), cigar_op(1, 4)]); // 9M1S (same query length as 10M)
}
fn mutate_mapq(b: &mut SamBuilder) {
    b.mapq(30);
}
fn mutate_rnext(b: &mut SamBuilder) {
    b.mate_ref_id(1);
}
fn mutate_pnext(b: &mut SamBuilder) {
    b.mate_pos(250);
}
fn mutate_tlen(b: &mut SamBuilder) {
    b.template_length(-150);
}
fn mutate_rname(b: &mut SamBuilder) {
    // The record's own reference id (RNAME), as opposed to `mutate_rnext`'s mate
    // reference id (RNEXT) above. `core_record_variant`'s baseline sets `ref_id(0)`
    // (`chr1` under `two_ref_header`); this moves it to `chr2`. RNAME is excluded
    // from `RecordKey` for primary alignments (see `record_key.rs`'s module docs —
    // "a primary record whose POS legitimately differs... must still pair with its
    // counterpart"; the same reasoning applies to RNAME), so this must surface as a
    // content diff, not a key mismatch.
    b.ref_id(1);
}

#[test]
fn positional_baseline_identical_records_match() {
    // Sanity check that the shared baseline really is a MATCH before any mutation
    // case relies on "baseline vs baseline+mutation" as its MATCH/DIFFER pair.
    let tmp = TempDir::new().unwrap();
    let header = two_ref_header();
    let records = vec![core_record_variant(|_| {})];
    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records);
    write_bam(&bam2, &header, &records);
    let outcome = positional_compare(&bam1, &bam2, 1, 64, 10, ContentPredicate::Exact)
        .expect("positional_compare should succeed");
    assert!(outcome.is_match(), "unmutated baseline must MATCH: {outcome:?}");
}

#[rstest]
#[case::seq(mutate_seq as fn(&mut SamBuilder), "SEQ")]
#[case::qual(mutate_qual as fn(&mut SamBuilder), "QUAL")]
#[case::pos(mutate_pos as fn(&mut SamBuilder), "POS")]
#[case::flag(mutate_flag as fn(&mut SamBuilder), "FLAG")]
#[case::cigar(mutate_cigar as fn(&mut SamBuilder), "CIGAR")]
#[case::mapq(mutate_mapq as fn(&mut SamBuilder), "MAPQ")]
#[case::rnext(mutate_rnext as fn(&mut SamBuilder), "RNEXT")]
#[case::pnext(mutate_pnext as fn(&mut SamBuilder), "PNEXT")]
#[case::tlen(mutate_tlen as fn(&mut SamBuilder), "TLEN")]
#[case::rname(mutate_rname as fn(&mut SamBuilder), "RNAME")]
fn positional_exact_core_field_mutation_differs(
    #[case] mutate: fn(&mut SamBuilder),
    #[case] field_name: &str,
) {
    let tmp = TempDir::new().unwrap();
    let header = two_ref_header();
    let rec_a = core_record_variant(|_| {});
    let rec_b = core_record_variant(mutate);

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &[rec_a]);
    write_bam(&bam2, &header, &[rec_b]);

    let outcome = positional_compare(&bam1, &bam2, 1, 64, 10, ContentPredicate::Exact)
        .expect("positional_compare should succeed");
    assert_eq!(
        outcome.key_mismatch_at, None,
        "a core-field-only mutation must not change RecordKey: {outcome:?}"
    );
    assert_eq!(outcome.content_diffs, 1, "exactly one content diff expected: {outcome:?}");
    assert!(!outcome.is_match(), "a core-field mutation must DIFFER: {outcome:?}");
    assert!(
        outcome.diff_details.iter().any(|d| d.contains(&format!("{field_name}:"))),
        "diff must name the {field_name} field, got: {outcome:?}"
    );
}

#[test]
fn positional_exact_qname_change_is_caught_as_key_mismatch() {
    // QNAME is a compared core field (see `FIELD_NAMES`), but unlike the other ten
    // core fields exercised by the rstest above, it cannot appear as a *content*
    // diff: `RecordKey` includes the read name for identity pairing (see
    // `record_key.rs`), so a QNAME-only mutation flips the key itself. The
    // positional engine's no-resync design must catch this as a key mismatch and
    // stop pairing, not silently re-pair by index and report zero diffs.
    let tmp = TempDir::new().unwrap();
    let header = two_ref_header();
    let rec_a = core_record_variant(|_| {});
    let rec_b = core_record_variant(|b| {
        b.read_name(b"read2");
    });

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &[rec_a]);
    write_bam(&bam2, &header, &[rec_b]);

    let outcome = positional_compare(&bam1, &bam2, 1, 64, 10, ContentPredicate::Exact)
        .expect("positional_compare should succeed");
    assert_eq!(
        outcome.key_mismatch_at,
        Some(0),
        "a QNAME-only mutation must flip RecordKey: {outcome:?}"
    );
    assert_eq!(
        outcome.content_diffs, 0,
        "pairing must stop at the key mismatch before any content diff is recorded: {outcome:?}"
    );
    assert!(!outcome.is_match(), "a QNAME change must DIFFER: {outcome:?}");
}

// ---- tags -------------------------------------------------------------

fn record_with_optional_rx(has_rx: bool) -> RawRecord {
    let mut b = SamBuilder::new();
    b.read_name(b"read1").sequence(b"ACGT").qualities(&[30; 4]).add_int_tag(SamTag::NM, 1);
    if has_rx {
        b.add_string_tag(SamTag::RX, b"AAAA");
    }
    b.build()
}

fn record_with_nm(nm: i32) -> RawRecord {
    let mut b = SamBuilder::new();
    b.read_name(b"read1").sequence(b"ACGT").qualities(&[30; 4]).add_int_tag(SamTag::NM, nm);
    b.build()
}

fn record_with_array_tag(values: &[u8]) -> RawRecord {
    // `ML` (base modification probability array, `B:C`) is a real SAM-spec array tag
    // unrelated to any other mutation case in this file -- a convenient stand-in for
    // "some `B:C`-typed array tag" rather than inventing a synthetic non-tag literal.
    let mut b = SamBuilder::new();
    b.read_name(b"read1").sequence(b"ACGT").qualities(&[30; 4]).add_array_u8(SamTag::ML, values);
    b.build()
}

/// Standalone positional/Exact scenarios that reduce cleanly to a single
/// match/DIFFER verdict via the shared [`positional_differs`] helper — tag drop,
/// tag value change, array-tag element change, tag reorder (the one case that
/// must NOT flip the verdict — `Exact` compares tags order-independently, see
/// `content.rs`'s `tag_order_match`), and a dropped trailing record and MI value
/// change (both plain content diffs under `Exact`). Cases that additionally pin
/// down *how* the engine reported a divergence (an exact `key_mismatch_at`
/// index, or exact `bam1_count`/`bam2_count`) stay as their own `#[test]`s below
/// since `positional_differs` only surfaces the final match/DIFFER bool.
#[rstest]
#[case::tag_drop(
    vec![record_with_optional_rx(true)],
    vec![record_with_optional_rx(false)],
    ContentPredicate::Exact,
    false
)]
#[case::tag_value_change(
    vec![record_with_nm(1)],
    vec![record_with_nm(2)],
    ContentPredicate::Exact,
    false
)]
#[case::array_tag_element_change(
    vec![record_with_array_tag(&[1, 2, 3, 4])],
    vec![record_with_array_tag(&[1, 2, 9, 4])],
    ContentPredicate::Exact,
    false
)]
#[case::tag_reorder_is_not_a_diff(
    vec![{
        let mut b = SamBuilder::new();
        b.read_name(b"read1")
            .sequence(b"ACGT")
            .qualities(&[30; 4])
            .add_int_tag(SamTag::NM, 1)
            .add_string_tag(SamTag::RX, b"AAAA");
        b.build()
    }],
    vec![{
        let mut b = SamBuilder::new();
        b.read_name(b"read1")
            .sequence(b"ACGT")
            .qualities(&[30; 4])
            .add_string_tag(SamTag::RX, b"AAAA")
            .add_int_tag(SamTag::NM, 1);
        b.build()
    }],
    ContentPredicate::Exact,
    true
)]
#[case::dropped_trailing_record(
    vec![
        keyed_record(b"read1", flags::FIRST_SEGMENT, 100),
        keyed_record(b"read2", flags::FIRST_SEGMENT, 200),
    ],
    vec![keyed_record(b"read1", flags::FIRST_SEGMENT, 100)], // read2 dropped
    ContentPredicate::Exact,
    false
)]
#[case::mi_value_change(
    vec![{
        let mut b = SamBuilder::new();
        b.read_name(b"read1").sequence(b"ACGT").qualities(&[30; 4]).add_int_tag(SamTag::MI, 5);
        b.build()
    }],
    vec![{
        let mut b = SamBuilder::new();
        b.read_name(b"read1").sequence(b"ACGT").qualities(&[30; 4]).add_int_tag(SamTag::MI, 99);
        b.build()
    }],
    ContentPredicate::Exact,
    false
)]
fn positional_exact_standalone_scenarios(
    #[case] recs_a: Vec<RawRecord>,
    #[case] recs_b: Vec<RawRecord>,
    #[case] pred: ContentPredicate,
    #[case] expect_match: bool,
) {
    let header = create_minimal_header("chr1", 10000);
    let differed = positional_differs(&header, &recs_a, &recs_b, pred);
    assert_eq!(
        !differed, expect_match,
        "recs_a vs recs_b under {pred:?}: expected match={expect_match}, got differed={differed}"
    );
}

// ---- records / order ----------------------------------------------------

/// A simple mapped, paired-flag record with a given name/segment/position — distinct
/// names give distinct `RecordKey`s, which is what the record/order mutations below need.
fn keyed_record(name: &[u8], segment_flag: u16, pos: i32) -> RawRecord {
    let mut b = SamBuilder::new();
    b.read_name(name)
        .flags(flags::PAIRED | segment_flag)
        .ref_id(0)
        .pos(pos - 1)
        .mapq(60)
        .sequence(b"ACGT")
        .qualities(&[30; 4]);
    b.build()
}

#[test]
fn positional_exact_duplicated_record_differs() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    let r1 = keyed_record(b"read1", flags::FIRST_SEGMENT, 100);
    let r2 = keyed_record(b"read2", flags::FIRST_SEGMENT, 200);
    let records1 = vec![r1.clone(), r2.clone()];
    let records2 = vec![r1, r2.clone(), r2]; // read2 duplicated
    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let outcome = positional_compare(&bam1, &bam2, 1, 64, 10, ContentPredicate::Exact)
        .expect("positional_compare should succeed");
    assert_eq!(outcome.bam1_count, 2);
    assert_eq!(outcome.bam2_count, 3);
    assert!(!outcome.is_match(), "a duplicated record must DIFFER: {outcome:?}");
}

#[test]
fn positional_exact_swap_of_distinct_key_records_differs() {
    // The swap ≡ drop+add masking test: read2 and read3 have distinct RecordKeys
    // (different QNAMEs) and are adjacent-swapped between the two files. A
    // resyncing comparator could match read2(file1) against read2(file2) a
    // position later and see no diff at all; the no-resync positional engine must
    // instead catch the desync at the exact index where the keys first disagree.
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    let r1 = keyed_record(b"read1", flags::FIRST_SEGMENT, 100);
    let r2 = keyed_record(b"read2", flags::FIRST_SEGMENT, 200);
    let r3 = keyed_record(b"read3", flags::FIRST_SEGMENT, 300);
    let records1 = vec![r1.clone(), r2.clone(), r3.clone()];
    let records2 = vec![r1, r3, r2]; // read2/read3 swapped
    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let outcome = positional_compare(&bam1, &bam2, 1, 64, 10, ContentPredicate::Exact)
        .expect("positional_compare should succeed");
    assert_eq!(
        outcome.key_mismatch_at,
        Some(1),
        "the swap must be caught as a key mismatch at the first disagreeing index: {outcome:?}"
    );
    assert!(!outcome.is_match(), "a key-distinct swap must DIFFER: {outcome:?}");
}

// ---- ExactMinusMi via positional_compare --------------------------------
//
// `ExactMinusMi` is already unit-tested directly against `content_diffs` in
// `engines/content.rs`, but until now it was only exercised end-to-end via
// `keyjoin_compare` (Section 4) — never through `positional_compare`, even
// though `ExactMinusMi` is routed through the positional engine in production
// (the `group` comparison preset). These two cases close that gap.

fn record_with_mi_and_nm(mi: i32, nm: i32) -> RawRecord {
    let mut b = SamBuilder::new();
    b.read_name(b"read1")
        .sequence(b"ACGT")
        .qualities(&[30; 4])
        .add_int_tag(SamTag::MI, mi)
        .add_int_tag(SamTag::NM, nm);
    b.build()
}

#[test]
fn positional_exact_minus_mi_tolerates_mi_only_value_change() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &[record_with_mi_and_nm(5, 1)]);
    write_bam(&bam2, &header, &[record_with_mi_and_nm(99, 1)]);

    let outcome = positional_compare(&bam1, &bam2, 1, 64, 10, ContentPredicate::ExactMinusMi)
        .expect("positional_compare should succeed");
    assert_eq!(outcome.content_diffs, 0, "{outcome:?}");
    assert!(
        outcome.is_match(),
        "an MI-only value change must NOT DIFFER under ExactMinusMi via positional_compare: {outcome:?}"
    );
}

#[test]
fn positional_exact_minus_mi_still_flags_non_mi_tag_change() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &[record_with_mi_and_nm(5, 1)]);
    write_bam(&bam2, &header, &[record_with_mi_and_nm(5, 2)]);

    let outcome = positional_compare(&bam1, &bam2, 1, 64, 10, ContentPredicate::ExactMinusMi)
        .expect("positional_compare should succeed");
    assert_eq!(outcome.content_diffs, 1, "{outcome:?}");
    assert!(
        !outcome.is_match(),
        "a non-MI tag change must still DIFFER under ExactMinusMi via positional_compare: {outcome:?}"
    );
}

// =============================================================================
// Section 2 — `ContentPredicate::ExactConsensus` (consensus saturation carve-out)
// =============================================================================

fn consensus_record_cd(cd: i32) -> RawRecord {
    let mut b = SamBuilder::new();
    b.read_name(b"read1").sequence(b"ACGT").qualities(&[30; 4]).add_int_tag(SamTag::CD, cd);
    b.build()
}

#[rstest]
#[case::identical_matches_under_exact(100, 100, ContentPredicate::Exact, true)]
#[case::saturated_pair_tolerated_under_consensus(
    32767,
    59920,
    ContentPredicate::ExactConsensus,
    true
)]
#[case::saturated_pair_still_differs_under_plain_exact(
    32767,
    59920,
    ContentPredicate::Exact,
    false
)]
#[case::non_saturated_pair_still_differs_under_consensus(
    100,
    200,
    ContentPredicate::ExactConsensus,
    false
)]
fn positional_consensus_cd_saturation_carveout(
    #[case] cd_a: i32,
    #[case] cd_b: i32,
    #[case] pred: ContentPredicate,
    #[case] expect_match: bool,
) {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &[consensus_record_cd(cd_a)]);
    write_bam(&bam2, &header, &[consensus_record_cd(cd_b)]);

    let outcome = positional_compare(&bam1, &bam2, 1, 64, 10, pred)
        .expect("positional_compare should succeed");
    assert_eq!(outcome.is_match(), expect_match, "cD {cd_a} vs {cd_b} under {pred:?}: {outcome:?}");
}

#[test]
fn positional_consensus_cm_and_ce_tolerated_alongside_saturated_cd() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    let rec = |cd: i32, cm: i32, ce: f32| {
        let mut b = SamBuilder::new();
        b.read_name(b"read1")
            .sequence(b"ACGT")
            .qualities(&[30; 4])
            .add_int_tag(SamTag::CD, cd)
            .add_int_tag(SamTag::CM, cm)
            .add_float_tag(SamTag::CE, ce);
        b.build()
    };
    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &[rec(32767, 32767, 0.01)]);
    write_bam(&bam2, &header, &[rec(59920, 40000, 0.02)]);

    let outcome = positional_compare(&bam1, &bam2, 1, 64, 10, ContentPredicate::ExactConsensus)
        .expect("positional_compare should succeed");
    assert!(outcome.is_match(), "cM/cE must be tolerated alongside a saturated cD: {outcome:?}");
}

// ---- Duplex depth-saturation carve-out (aD/bD at 32767, cD at 65534) ----
//
// fgbio's *duplex* consensus additionally saturates the per-strand `aD`/`bD` at
// `i16::MAX` (32767) and the *combined* `cD` at `2 * i16::MAX` = 65534 (the sum of two
// independently-saturated `Short`s). Real Phase-7 duplex data measured `aD 59034 vs
// 32767` and `cD 118954 vs 65534`.

fn consensus_record_with_tag(tag: SamTag, value: i32) -> RawRecord {
    let mut b = SamBuilder::new();
    b.read_name(b"read1").sequence(b"ACGT").qualities(&[30; 4]).add_int_tag(tag, value);
    b.build()
}

#[rstest]
#[case::ad_saturated_32767(SamTag::AD, 32767, 59034, ContentPredicate::ExactConsensus, true)]
#[case::bd_saturated_32767(SamTag::BD, 32767, 61000, ContentPredicate::ExactConsensus, true)]
#[case::cd_saturated_65534_duplex(
    SamTag::CD,
    65534,
    118_954,
    ContentPredicate::ExactConsensus,
    true
)]
#[case::ad_non_saturated_still_differs(
    SamTag::AD,
    100,
    200,
    ContentPredicate::ExactConsensus,
    false
)]
#[case::cd_saturated_65534_still_differs_under_plain_exact(
    SamTag::CD,
    65534,
    118_954,
    ContentPredicate::Exact,
    false
)]
fn positional_consensus_duplex_depth_saturation_carveout(
    #[case] tag: SamTag,
    #[case] value_a: i32,
    #[case] value_b: i32,
    #[case] pred: ContentPredicate,
    #[case] expect_match: bool,
) {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &[consensus_record_with_tag(tag, value_a)]);
    write_bam(&bam2, &header, &[consensus_record_with_tag(tag, value_b)]);

    let outcome = positional_compare(&bam1, &bam2, 1, 64, 10, pred)
        .expect("positional_compare should succeed");
    assert_eq!(
        outcome.is_match(),
        expect_match,
        "{tag:?} {value_a} vs {value_b} under {pred:?}: {outcome:?}"
    );
}

#[test]
fn positional_consensus_ad_undercount_relative_to_saturation_still_differs() {
    // An UNDERCOUNT relative to a saturated value (fgbio saturated at 32767, fgumi
    // reports something *smaller*) must never be tolerated -- this is not the "true
    // depth is at least as large" shape the carve-out permits, so it must still DIFFER
    // even under ExactConsensus. This is the soundness proof for finding B: a real
    // depth-undercount regression on a saturated family must still be caught.
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &[consensus_record_with_tag(SamTag::AD, 32767)]);
    write_bam(&bam2, &header, &[consensus_record_with_tag(SamTag::AD, 5)]);

    let outcome = positional_compare(&bam1, &bam2, 1, 64, 10, ContentPredicate::ExactConsensus)
        .expect("positional_compare should succeed");
    assert!(
        !outcome.is_match(),
        "an aD undercount relative to a saturated value must still DIFFER: {outcome:?}"
    );
}

fn per_base_consensus_record(cd_bases: &[i16]) -> RawRecord {
    let mut b = SamBuilder::new();
    b.read_name(b"read1")
        .sequence(b"ACGT")
        .qualities(&[30; 4])
        .add_array_i16(SamTag::CD_BASES, cd_bases);
    b.build()
}

#[rstest]
#[case::exact(ContentPredicate::Exact)]
#[case::exact_consensus(ContentPredicate::ExactConsensus)]
fn positional_per_base_consensus_tag_mutation_always_differs(#[case] pred: ContentPredicate) {
    // The per-base `cd`/`ce` (`B:s`) arrays saturate identically in both tools and
    // are NOT part of the per-read `cD`/`cM`/`cE` saturation carve-out — a change
    // here must DIFFER under every predicate, including `ExactConsensus`.
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &[per_base_consensus_record(&[10, 10, 10, 10])]);
    write_bam(&bam2, &header, &[per_base_consensus_record(&[10, 10, 9, 10])]);

    let outcome = positional_compare(&bam1, &bam2, 1, 64, 10, pred)
        .expect("positional_compare should succeed");
    assert!(
        !outcome.is_match(),
        "a per-base consensus tag change must DIFFER under {pred:?}: {outcome:?}"
    );
}

// =============================================================================
// Section 3 — `compare_headers` (`@HD`/`@SQ`/`@RG` compared, `@PG`/`@CO` normalized)
// =============================================================================

fn header_with_sq(name: &str, len: usize) -> Header {
    let seq = Map::<ReferenceSequence>::new(NonZeroUsize::new(len).expect("non-zero"));
    Header::builder().add_reference_sequence(bstr::BString::from(name), seq).build()
}

fn header_with_two_sq(first: (&str, usize), second: (&str, usize)) -> Header {
    let s1 = Map::<ReferenceSequence>::new(NonZeroUsize::new(first.1).expect("non-zero"));
    let s2 = Map::<ReferenceSequence>::new(NonZeroUsize::new(second.1).expect("non-zero"));
    Header::builder()
        .add_reference_sequence(bstr::BString::from(first.0), s1)
        .add_reference_sequence(bstr::BString::from(second.0), s2)
        .build()
}

fn header_with_hd(so: Option<&str>, go: Option<&str>, ss: Option<&str>) -> Header {
    let mut hd = Map::<HeaderRecord>::default();
    if let Some(so) = so {
        hd.other_fields_mut().insert(hd_tag::SORT_ORDER, bstr::BString::from(so));
    }
    if let Some(go) = go {
        hd.other_fields_mut().insert(hd_tag::GROUP_ORDER, bstr::BString::from(go));
    }
    if let Some(ss) = ss {
        hd.other_fields_mut().insert(hd_tag::SUBSORT_ORDER, bstr::BString::from(ss));
    }
    Header::builder().set_header(hd).build()
}

fn header_with_rg_sample(id: &str, sample: &str) -> Header {
    let mut rg = Map::<ReadGroup>::default();
    rg.other_fields_mut().insert(rg_tag::SAMPLE, bstr::BString::from(sample));
    Header::builder().add_read_group(bstr::BString::from(id), rg).build()
}

fn header_with_pg_version(id: &str, version: &str) -> Header {
    let mut pg = Map::<Program>::default();
    pg.other_fields_mut().insert(pg_tag::VERSION, bstr::BString::from(version));
    Header::builder().add_program(bstr::BString::from(id), pg).build()
}

fn header_with_comment(text: &str) -> Header {
    Header::builder().add_comment(bstr::BString::from(text)).build()
}

#[test]
fn compare_headers_sq_name_change_differs() {
    let h1 = header_with_sq("chr1", 1000);
    let h2 = header_with_sq("chr2", 1000);
    assert!(compare_headers(&h1, &h1).is_none(), "baseline must MATCH itself");
    let diffs = compare_headers(&h1, &h2).expect("@SQ name divergence must DIFFER");
    assert!(diffs.iter().any(|d| d.starts_with("@SQ")), "{diffs:?}");
}

#[test]
fn compare_headers_sq_length_change_differs() {
    let h1 = header_with_sq("chr1", 1000);
    let h2 = header_with_sq("chr1", 2000);
    let diffs = compare_headers(&h1, &h2).expect("@SQ length divergence must DIFFER");
    assert!(diffs.iter().any(|d| d.starts_with("@SQ")), "{diffs:?}");
}

#[test]
fn compare_headers_sq_order_change_differs() {
    let h1 = header_with_two_sq(("chr1", 1000), ("chr2", 2000));
    let h2 = header_with_two_sq(("chr2", 2000), ("chr1", 1000));
    let diffs = compare_headers(&h1, &h2).expect("@SQ order divergence must DIFFER");
    assert!(diffs.iter().any(|d| d.starts_with("@SQ")), "{diffs:?}");
}

#[test]
fn compare_headers_hd_ss_change_differs() {
    let h1 = header_with_hd(Some("queryname"), None, Some("queryname:natural"));
    let h2 = header_with_hd(Some("queryname"), None, Some("queryname:lexicographical"));
    assert!(compare_headers(&h1, &h1).is_none(), "baseline must MATCH itself");
    let diffs = compare_headers(&h1, &h2).expect("@HD SS divergence must DIFFER");
    assert!(diffs.iter().any(|d| d.starts_with("@HD")), "{diffs:?}");
}

#[test]
fn compare_headers_rg_sample_change_differs() {
    let h1 = header_with_rg_sample("rg1", "sampleA");
    let h2 = header_with_rg_sample("rg1", "sampleB");
    let diffs = compare_headers(&h1, &h2).expect("@RG SM divergence must DIFFER");
    assert!(diffs.iter().any(|d| d.starts_with("@RG")), "{diffs:?}");
}

#[test]
fn compare_headers_pg_version_change_is_not_a_diff() {
    // @PG is tool-invocation metadata (name/version/command line) that legitimately
    // differs between fgumi and fgbio even on functionally identical output --
    // it is normalized away entirely, never compared.
    let h1 = header_with_pg_version("fgumi", "1.0.0");
    let h2 = header_with_pg_version("fgbio", "4.1.0");
    assert!(compare_headers(&h1, &h2).is_none(), "@PG must be normalized, not compared");
}

#[test]
fn compare_headers_co_change_is_not_a_diff() {
    let h1 = header_with_comment("built by fgumi");
    let h2 = header_with_comment("built by fgbio");
    assert!(compare_headers(&h1, &h2).is_none(), "@CO must be normalized, not compared");
}

// =============================================================================
// Section 4 — `keyjoin_compare` (`group`: merge-join + MI bijection)
// =============================================================================
//
// `keyjoin_cfg` and `mi_record` are shared with `test_compare_bams.rs` and now
// live in `helpers::bam_generator` (see the top-of-file import).

#[test]
fn keyjoin_mi_renumber_preserving_partition_is_not_a_diff() {
    // The accepted divergence `group` key-join exists for: fgumi and fgbio number
    // molecules independently, so a pure MI relabel that preserves which reads
    // share a molecule must NOT DIFFER.
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    let records1 = vec![
        mi_record(b"read1", 100, "1"),
        mi_record(b"read2", 200, "1"),
        mi_record(b"read3", 300, "2"),
    ];
    let records2 = vec![
        mi_record(b"read1", 100, "5"),
        mi_record(b"read2", 200, "5"),
        mi_record(b"read3", 300, "9"),
    ];
    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let outcome =
        keyjoin_compare(&bam1, &bam2, &keyjoin_cfg(&tmp)).expect("keyjoin_compare should succeed");
    assert_eq!(outcome.mi_bijection_mismatches, 0, "{outcome:?}");
    assert!(outcome.is_match(), "MI renumbering alone must NOT DIFFER: {outcome:?}");
}

#[test]
fn keyjoin_mi_bijection_violation_differs() {
    // Same reads as the tolerated case above, but read2/read3 are split into two
    // molecules on one side: a genuine grouping (partition) change, which the MI
    // bijection tracker must catch even though content and record presence are
    // otherwise untouched.
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    let records1 = vec![mi_record(b"read1", 100, "1"), mi_record(b"read2", 200, "1")];
    let records2 = vec![mi_record(b"read1", 100, "1"), mi_record(b"read2", 200, "2")];
    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let outcome =
        keyjoin_compare(&bam1, &bam2, &keyjoin_cfg(&tmp)).expect("keyjoin_compare should succeed");
    assert!(outcome.mi_bijection_mismatches > 0, "{outcome:?}");
    assert!(!outcome.is_match(), "a grouping split must DIFFER: {outcome:?}");
}

#[test]
fn keyjoin_reordered_records_same_grouping_is_not_a_diff() {
    // `group`'s output order legitimately differs between tools (they number
    // molecules independently); the key-join canonicalizes both sides to
    // queryname order internally, so a pure reorder with an unchanged partition
    // must NOT DIFFER.
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    let r1 = mi_record(b"read1", 100, "1");
    let r2 = mi_record(b"read2", 200, "1");
    let r3 = mi_record(b"read3", 300, "2");
    let records1 = vec![r1.clone(), r2.clone(), r3.clone()];
    let records2 = vec![r3, r1, r2];
    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let outcome =
        keyjoin_compare(&bam1, &bam2, &keyjoin_cfg(&tmp)).expect("keyjoin_compare should succeed");
    assert!(
        outcome.is_match(),
        "reordering with an unchanged partition must NOT DIFFER: {outcome:?}"
    );
}

#[test]
fn keyjoin_content_diff_with_intact_grouping_differs() {
    // BS1 regression proof: a non-MI content bug on one matched pair, with the MI
    // partition completely untouched, must still DIFFER -- the key-join's job is
    // not only MI equivalence but content equality too.
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    let good = |name: &[u8], pos: i32, seq: &[u8]| {
        let mut b = SamBuilder::new();
        b.read_name(name)
            .sequence(seq)
            .qualities(&vec![30; seq.len()])
            .ref_id(0)
            .pos(pos - 1)
            .mapq(60)
            .add_string_tag(SamTag::MI, b"1");
        b.build()
    };
    let records1 = vec![good(b"read1", 100, b"ACGTACGT"), good(b"read2", 200, b"ACGTACGT")];
    let records2 = vec![good(b"read1", 100, b"ACGTACGT"), good(b"read2", 200, b"ACGTACGA")];
    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let outcome =
        keyjoin_compare(&bam1, &bam2, &keyjoin_cfg(&tmp)).expect("keyjoin_compare should succeed");
    assert_eq!(outcome.mi_bijection_mismatches, 0, "the grouping itself is untouched: {outcome:?}");
    assert_eq!(outcome.content_diffs, 1, "{outcome:?}");
    assert!(
        !outcome.is_match(),
        "a content diff must DIFFER even with intact grouping: {outcome:?}"
    );
}

#[test]
fn keyjoin_dropped_record_is_a_presence_diff() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    let records1 = vec![mi_record(b"read1", 100, "1"), mi_record(b"read2", 200, "1")];
    let records2 = vec![mi_record(b"read1", 100, "1")]; // read2 dropped
    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let outcome =
        keyjoin_compare(&bam1, &bam2, &keyjoin_cfg(&tmp)).expect("keyjoin_compare should succeed");
    assert_eq!(outcome.only_in_bam1, 1, "{outcome:?}");
    assert!(!outcome.is_match(), "a dropped record must DIFFER: {outcome:?}");
}

// =============================================================================
// Section 5 — `sort_verify_compare` (`sort`: independent-baseline order check)
// =============================================================================

fn coord_record(name: &[u8], pos: i32) -> RawRecord {
    let mut b = SamBuilder::new();
    b.read_name(name).ref_id(0).pos(pos - 1).mapq(60).sequence(b"ACGT").qualities(&[30; 4]);
    b.build()
}

#[test]
fn sort_verify_baseline_identical_correctly_sorted_bams_match() {
    let tmp = TempDir::new().unwrap();
    let header = create_coordinate_sorted_header("chr1", 10000);
    let records = vec![coord_record(b"r1", 50), coord_record(b"r2", 100), coord_record(b"r3", 150)];
    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records);
    write_bam(&bam2, &header, &records);

    let outcome =
        sort_verify_compare(&bam1, &bam2, 10).expect("sort_verify_compare should succeed");
    assert!(outcome.is_match(), "an unmutated, correctly-sorted baseline must MATCH: {outcome:?}");
}

#[test]
fn sort_verify_order_violation_differs() {
    let tmp = TempDir::new().unwrap();
    let header = create_coordinate_sorted_header("chr1", 10000);
    let records1 =
        vec![coord_record(b"r1", 50), coord_record(b"r2", 100), coord_record(b"r3", 150)];
    // r3 now sorts before r2: bam2 is no longer non-decreasing under (tid, pos).
    let records2 = vec![coord_record(b"r1", 50), coord_record(b"r2", 100), coord_record(b"r3", 75)];
    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let outcome =
        sort_verify_compare(&bam1, &bam2, 10).expect("sort_verify_compare should succeed");
    assert!(outcome.bam2_violations > 0, "{outcome:?}");
    assert!(!outcome.is_match(), "a sort-order violation must DIFFER: {outcome:?}");
}

#[test]
fn sort_verify_swap_within_equal_key_run_is_not_a_diff() {
    // Two records tied on the same (tid, pos): coordinate sort does not tie-break
    // on read name, so their relative order is unconstrained (the SORT-01 tie
    // residue) -- swapping them must NOT DIFFER.
    let tmp = TempDir::new().unwrap();
    let header = create_coordinate_sorted_header("chr1", 10000);
    let a = coord_record(b"readA", 100);
    let b_rec = coord_record(b"readB", 100);
    let records1 = vec![a.clone(), b_rec.clone()];
    let records2 = vec![b_rec, a];
    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let outcome =
        sort_verify_compare(&bam1, &bam2, 10).expect("sort_verify_compare should succeed");
    assert_eq!(outcome.run_mismatches, 0, "{outcome:?}");
    assert!(outcome.is_match(), "a tie-run swap must NOT DIFFER: {outcome:?}");
}

#[test]
fn sort_verify_dropped_record_within_equal_key_run_differs() {
    let tmp = TempDir::new().unwrap();
    let header = create_coordinate_sorted_header("chr1", 10000);
    let a = coord_record(b"readA", 100);
    let b_rec = coord_record(b"readB", 100);
    let records1 = vec![a.clone(), b_rec];
    let records2 = vec![a]; // readB removed from the tied run
    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let outcome =
        sort_verify_compare(&bam1, &bam2, 10).expect("sort_verify_compare should succeed");
    assert!(outcome.run_mismatches > 0, "{outcome:?}");
    assert!(!outcome.is_match(), "removing a record from a tied run must DIFFER: {outcome:?}");
}

// =============================================================================
// Section 6 — `fgumi compare metrics` (row key-join)
// =============================================================================

/// Runs `CompareMetrics::execute()` in-process; `true` on match, `false` on
/// `CompareMismatch` (DIFFER). Panics on any other error.
fn run_metrics_compare(file1: &Path, file2: &Path) -> bool {
    let cmd = CompareMetrics::try_parse_from([
        "metrics",
        file1.to_str().unwrap(),
        file2.to_str().unwrap(),
    ])
    .expect("failed to parse compare metrics args");
    match cmd.execute("fgumi compare metrics") {
        Ok(()) => true,
        Err(e) if e.is::<CompareMismatch>() => false,
        Err(e) => panic!("compare metrics hit unexpected error: {e:#}"),
    }
}

#[test]
fn metrics_baseline_identical_files_match() {
    let tmp = TempDir::new().unwrap();
    let content = "umi\tcount\nAAA\t3\nCCC\t5\n";
    let f1 = write_tsv(tmp.path(), "a.txt", content);
    let f2 = write_tsv(tmp.path(), "b.txt", content);
    assert!(run_metrics_compare(&f1, &f2), "an unmutated baseline must MATCH");
}

#[test]
fn metrics_integer_count_change_differs() {
    let tmp = TempDir::new().unwrap();
    let f1 = write_tsv(tmp.path(), "a.txt", "umi\tcount\nAAA\t3\n");
    let f2 = write_tsv(tmp.path(), "b.txt", "umi\tcount\nAAA\t4\n");
    assert!(!run_metrics_compare(&f1, &f2), "an integer count change must DIFFER");
}

#[test]
fn metrics_float_beyond_epsilon_differs() {
    let tmp = TempDir::new().unwrap();
    let f1 = write_tsv(tmp.path(), "a.txt", "key\tvalue\nA\t1.0\n");
    let f2 = write_tsv(tmp.path(), "b.txt", "key\tvalue\nA\t1.5\n");
    assert!(!run_metrics_compare(&f1, &f2), "a float difference beyond tolerance must DIFFER");
}

#[test]
fn metrics_float_within_epsilon_representation_only_is_not_a_diff() {
    // fgbio's rounded `DecimalFormat("0.######")` output vs fgumi's full-precision
    // float -- the one accepted metrics divergence.
    let tmp = TempDir::new().unwrap();
    let f1 = write_tsv(tmp.path(), "a.txt", "key\tvalue\nA\t0.333333\n");
    let f2 = write_tsv(tmp.path(), "b.txt", "key\tvalue\nA\t0.3333333333333333\n");
    assert!(run_metrics_compare(&f1, &f2), "float representation-only must NOT DIFFER");
}

#[test]
fn metrics_dropped_row_differs() {
    let tmp = TempDir::new().unwrap();
    let f1 = write_tsv(tmp.path(), "a.txt", "key\tvalue\nA\t1\nB\t2\n");
    let f2 = write_tsv(tmp.path(), "b.txt", "key\tvalue\nA\t1\n");
    assert!(!run_metrics_compare(&f1, &f2), "a dropped row must DIFFER");
}

#[test]
fn metrics_added_row_differs() {
    let tmp = TempDir::new().unwrap();
    let f1 = write_tsv(tmp.path(), "a.txt", "key\tvalue\nA\t1\n");
    let f2 = write_tsv(tmp.path(), "b.txt", "key\tvalue\nA\t1\nB\t2\n");
    assert!(!run_metrics_compare(&f1, &f2), "an added row must DIFFER");
}

#[test]
fn metrics_permuted_row_order_is_not_a_diff() {
    let tmp = TempDir::new().unwrap();
    let f1 = write_tsv(tmp.path(), "a.txt", "key\tvalue\nA\t1\nB\t2\nC\t3\n");
    let f2 = write_tsv(tmp.path(), "b.txt", "key\tvalue\nC\t3\nA\t1\nB\t2\n");
    assert!(run_metrics_compare(&f1, &f2), "row order must NOT be a diff under the key-join");
}

#[test]
fn metrics_different_column_set_differs() {
    let tmp = TempDir::new().unwrap();
    let f1 = write_tsv(tmp.path(), "a.txt", "umi\tcount\tfraction\nAAA\t3\t0.5\n");
    let f2 = write_tsv(tmp.path(), "b.txt", "umi\tcount\nAAA\t3\n");
    assert!(!run_metrics_compare(&f1, &f2), "a differing column set must DIFFER");
}

// =============================================================================
// Coverage guard
// =============================================================================
//
// The catalog below used to be checked only against `REQUIRED_SUBSTRINGS`, a second
// hand-maintained list of strings -- two arrays checked against each other, neither
// tied to this file's actual `#[test]`/`#[case]` functions or to the real comparator.
// Deleting (or never writing) the test for a compared dimension left both arrays
// self-consistent and the guard green. Each `CatalogEntry` now also carries a
// `verify: fn() -> bool` closure that independently rebuilds a minimal fixture for
// that exact mutation and re-runs it through the real engine (`positional_compare`,
// `compare_headers`, `keyjoin_compare`, `sort_verify_compare`, or
// `CompareMetrics::execute`), so `mutation_catalog_covers_every_compared_dimension`
// re-derives every row's expected verdict from production code -- not merely from
// this array's own `must_differ` field -- and would fail if a compared
// field/tag/row-dimension actually stopped being compared, independent of whether
// any *other* test in this file still happens to cover it.

/// One row of the mutation catalog this harness implements: which engine/primitive
/// it exercises, a short name for the specific mutation, whether it must `DIFFER`
/// (`true`) or must NOT (`false`, an accepted/tolerated divergence), and a `verify`
/// closure that independently re-derives the actual verdict from the real
/// comparator (see the module-level comment above).
struct CatalogEntry {
    engine: &'static str,
    mutation: &'static str,
    must_differ: bool,
    /// Rebuilds a minimal fixture for this exact mutation and runs it through the
    /// real comparator, returning `true` if the comparator reported DIFFER. This is
    /// what makes the guard load-bearing: it does not consult `must_differ`, any
    /// other field on this struct, or any other test in the file.
    verify: fn() -> bool,
}

// ---- Shared fixture-and-compare helpers used by `verify` closures below -------
//
// Each helper writes a tiny BAM (or TSV) pair to its own `TempDir` and returns
// whether the real engine reported DIFFER, mirroring the section-specific tests
// above but reduced to a single bool so `MUTATION_CATALOG`'s `verify` closures stay
// one-liners.

/// Re-derives the "single core-field mutation must DIFFER" property (Section 1's
/// `positional_exact_core_field_mutation_differs` rstest body) via the real
/// `positional_compare` engine, for a given field-mutation closure.
fn verify_core_field_mutation_differs(mutate: fn(&mut SamBuilder)) -> bool {
    let tmp = TempDir::new().unwrap();
    let header = two_ref_header();
    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &[core_record_variant(|_| {})]);
    write_bam(&bam2, &header, &[core_record_variant(mutate)]);
    let outcome = positional_compare(&bam1, &bam2, 1, 64, 10, ContentPredicate::Exact)
        .expect("positional_compare should succeed");
    outcome.key_mismatch_at.is_none() && outcome.content_diffs == 1 && !outcome.is_match()
}

/// Runs `positional_compare` under `pred` over two single-record (or short) BAMs
/// built from `recs_a`/`recs_b`, returning `true` if it reported DIFFER.
fn positional_differs(
    header: &Header,
    recs_a: &[RawRecord],
    recs_b: &[RawRecord],
    pred: ContentPredicate,
) -> bool {
    let tmp = TempDir::new().unwrap();
    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, header, recs_a);
    write_bam(&bam2, header, recs_b);
    !positional_compare(&bam1, &bam2, 1, 64, 10, pred)
        .expect("positional_compare should succeed")
        .is_match()
}

/// Returns `true` if `compare_headers` reports a divergence between `h1` and `h2`.
fn headers_differ(h1: &Header, h2: &Header) -> bool {
    compare_headers(h1, h2).is_some()
}

/// Runs `keyjoin_compare` over two short BAMs built from `recs_a`/`recs_b`,
/// returning `true` if it reported DIFFER.
fn keyjoin_differs(recs_a: &[RawRecord], recs_b: &[RawRecord]) -> bool {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, recs_a);
    write_bam(&bam2, &header, recs_b);
    !keyjoin_compare(&bam1, &bam2, &keyjoin_cfg(&tmp))
        .expect("keyjoin_compare should succeed")
        .is_match()
}

/// Runs `sort_verify_compare` over two coordinate-sorted-header BAMs built from
/// `recs_a`/`recs_b`, returning `true` if it reported DIFFER.
fn sort_verify_differs(recs_a: &[RawRecord], recs_b: &[RawRecord]) -> bool {
    let tmp = TempDir::new().unwrap();
    let header = create_coordinate_sorted_header("chr1", 10000);
    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, recs_a);
    write_bam(&bam2, &header, recs_b);
    !sort_verify_compare(&bam1, &bam2, 10).expect("sort_verify_compare should succeed").is_match()
}

/// Runs `fgumi compare metrics` over two ad hoc TSV contents, returning `true` if
/// it reported DIFFER.
fn metrics_content_differs(content_a: &str, content_b: &str) -> bool {
    let tmp = TempDir::new().unwrap();
    let f1 = write_tsv(tmp.path(), "a.txt", content_a);
    let f2 = write_tsv(tmp.path(), "b.txt", content_b);
    !run_metrics_compare(&f1, &f2)
}

/// The full catalog of mutation cases implemented in this file. Each entry's
/// `verify` closure independently reruns the real comparator on a minimal fixture
/// -- see the section comment above for why that (not just this array's shape)
/// is what makes the guard test meaningful.
const MUTATION_CATALOG: &[CatalogEntry] = &[
    // Section 1: positional / Exact
    CatalogEntry {
        engine: "positional/Exact",
        mutation: "SEQ base flip",
        must_differ: true,
        verify: || verify_core_field_mutation_differs(mutate_seq),
    },
    CatalogEntry {
        engine: "positional/Exact",
        mutation: "QUAL byte change",
        must_differ: true,
        verify: || verify_core_field_mutation_differs(mutate_qual),
    },
    CatalogEntry {
        engine: "positional/Exact",
        mutation: "POS shift",
        must_differ: true,
        verify: || verify_core_field_mutation_differs(mutate_pos),
    },
    CatalogEntry {
        engine: "positional/Exact",
        mutation: "FLAG bit toggle",
        must_differ: true,
        verify: || verify_core_field_mutation_differs(mutate_flag),
    },
    CatalogEntry {
        engine: "positional/Exact",
        mutation: "CIGAR change",
        must_differ: true,
        verify: || verify_core_field_mutation_differs(mutate_cigar),
    },
    CatalogEntry {
        engine: "positional/Exact",
        mutation: "MAPQ change",
        must_differ: true,
        verify: || verify_core_field_mutation_differs(mutate_mapq),
    },
    CatalogEntry {
        engine: "positional/Exact",
        mutation: "RNEXT change",
        must_differ: true,
        verify: || verify_core_field_mutation_differs(mutate_rnext),
    },
    CatalogEntry {
        engine: "positional/Exact",
        mutation: "PNEXT change",
        must_differ: true,
        verify: || verify_core_field_mutation_differs(mutate_pnext),
    },
    CatalogEntry {
        engine: "positional/Exact",
        mutation: "TLEN change",
        must_differ: true,
        verify: || verify_core_field_mutation_differs(mutate_tlen),
    },
    CatalogEntry {
        engine: "positional/Exact",
        mutation: "RNAME change",
        must_differ: true,
        verify: || verify_core_field_mutation_differs(mutate_rname),
    },
    CatalogEntry {
        engine: "positional/Exact",
        mutation: "QNAME change (caught as key mismatch, not a content diff)",
        must_differ: true,
        verify: || {
            positional_differs(
                &two_ref_header(),
                &[core_record_variant(|_| {})],
                &[core_record_variant(|b| {
                    b.read_name(b"read2");
                })],
                ContentPredicate::Exact,
            )
        },
    },
    CatalogEntry {
        engine: "positional/Exact",
        mutation: "tag drop",
        must_differ: true,
        verify: || {
            positional_differs(
                &create_minimal_header("chr1", 10000),
                &[record_with_optional_rx(true)],
                &[record_with_optional_rx(false)],
                ContentPredicate::Exact,
            )
        },
    },
    CatalogEntry {
        engine: "positional/Exact",
        mutation: "tag value change",
        must_differ: true,
        verify: || {
            positional_differs(
                &create_minimal_header("chr1", 10000),
                &[record_with_nm(1)],
                &[record_with_nm(2)],
                ContentPredicate::Exact,
            )
        },
    },
    CatalogEntry {
        engine: "positional/Exact",
        mutation: "array-tag element change",
        must_differ: true,
        verify: || {
            positional_differs(
                &create_minimal_header("chr1", 10000),
                &[record_with_array_tag(&[1, 2, 3, 4])],
                &[record_with_array_tag(&[1, 2, 9, 4])],
                ContentPredicate::Exact,
            )
        },
    },
    CatalogEntry {
        engine: "positional/Exact",
        mutation: "tag reorder",
        must_differ: false,
        verify: || {
            let rec_a = {
                let mut b = SamBuilder::new();
                b.read_name(b"read1")
                    .sequence(b"ACGT")
                    .qualities(&[30; 4])
                    .add_int_tag(SamTag::NM, 1)
                    .add_string_tag(SamTag::RX, b"AAAA");
                b.build()
            };
            let rec_b = {
                let mut b = SamBuilder::new();
                b.read_name(b"read1")
                    .sequence(b"ACGT")
                    .qualities(&[30; 4])
                    .add_string_tag(SamTag::RX, b"AAAA")
                    .add_int_tag(SamTag::NM, 1);
                b.build()
            };
            positional_differs(
                &create_minimal_header("chr1", 10000),
                &[rec_a],
                &[rec_b],
                ContentPredicate::Exact,
            )
        },
    },
    CatalogEntry {
        engine: "positional/Exact",
        mutation: "drop a record",
        must_differ: true,
        verify: || {
            positional_differs(
                &create_minimal_header("chr1", 10000),
                &[
                    keyed_record(b"read1", flags::FIRST_SEGMENT, 100),
                    keyed_record(b"read2", flags::FIRST_SEGMENT, 200),
                ],
                &[keyed_record(b"read1", flags::FIRST_SEGMENT, 100)],
                ContentPredicate::Exact,
            )
        },
    },
    CatalogEntry {
        engine: "positional/Exact",
        mutation: "duplicate a record",
        must_differ: true,
        verify: || {
            let r1 = keyed_record(b"read1", flags::FIRST_SEGMENT, 100);
            let r2 = keyed_record(b"read2", flags::FIRST_SEGMENT, 200);
            positional_differs(
                &create_minimal_header("chr1", 10000),
                &[r1.clone(), r2.clone()],
                &[r1, r2.clone(), r2],
                ContentPredicate::Exact,
            )
        },
    },
    CatalogEntry {
        engine: "positional/Exact",
        mutation: "swap two distinct-key records",
        must_differ: true,
        verify: || {
            let r1 = keyed_record(b"read1", flags::FIRST_SEGMENT, 100);
            let r2 = keyed_record(b"read2", flags::FIRST_SEGMENT, 200);
            let r3 = keyed_record(b"read3", flags::FIRST_SEGMENT, 300);
            positional_differs(
                &create_minimal_header("chr1", 10000),
                &[r1.clone(), r2.clone(), r3.clone()],
                &[r1, r3, r2],
                ContentPredicate::Exact,
            )
        },
    },
    CatalogEntry {
        engine: "positional/Exact",
        mutation: "MI value change",
        must_differ: true,
        verify: || {
            let rec = |mi: i32| {
                let mut b = SamBuilder::new();
                b.read_name(b"read1")
                    .sequence(b"ACGT")
                    .qualities(&[30; 4])
                    .add_int_tag(SamTag::MI, mi);
                b.build()
            };
            positional_differs(
                &create_minimal_header("chr1", 10000),
                &[rec(5)],
                &[rec(99)],
                ContentPredicate::Exact,
            )
        },
    },
    CatalogEntry {
        engine: "positional/ExactMinusMi",
        mutation: "ExactMinusMi tolerates an MI-only value change",
        must_differ: false,
        verify: || {
            positional_differs(
                &create_minimal_header("chr1", 10000),
                &[record_with_mi_and_nm(5, 1)],
                &[record_with_mi_and_nm(99, 1)],
                ContentPredicate::ExactMinusMi,
            )
        },
    },
    CatalogEntry {
        engine: "positional/ExactMinusMi",
        mutation: "ExactMinusMi still flags a non-MI tag change",
        must_differ: true,
        verify: || {
            positional_differs(
                &create_minimal_header("chr1", 10000),
                &[record_with_mi_and_nm(5, 1)],
                &[record_with_mi_and_nm(5, 2)],
                ContentPredicate::ExactMinusMi,
            )
        },
    },
    // Section 2: ExactConsensus
    CatalogEntry {
        engine: "positional/ExactConsensus",
        mutation: "cD saturated pair tolerated",
        must_differ: false,
        verify: || {
            positional_differs(
                &create_minimal_header("chr1", 10000),
                &[consensus_record_cd(32767)],
                &[consensus_record_cd(59920)],
                ContentPredicate::ExactConsensus,
            )
        },
    },
    CatalogEntry {
        engine: "positional/Exact",
        mutation: "cD saturated pair still differs under plain Exact",
        must_differ: true,
        verify: || {
            positional_differs(
                &create_minimal_header("chr1", 10000),
                &[consensus_record_cd(32767)],
                &[consensus_record_cd(59920)],
                ContentPredicate::Exact,
            )
        },
    },
    CatalogEntry {
        engine: "positional/ExactConsensus",
        mutation: "cD non-saturated pair still differs",
        must_differ: true,
        verify: || {
            positional_differs(
                &create_minimal_header("chr1", 10000),
                &[consensus_record_cd(100)],
                &[consensus_record_cd(200)],
                ContentPredicate::ExactConsensus,
            )
        },
    },
    CatalogEntry {
        engine: "positional/ExactConsensus",
        mutation: "cM/cE tolerated alongside saturated cD",
        must_differ: false,
        verify: || {
            let rec = |cd: i32, cm: i32, ce: f32| {
                let mut b = SamBuilder::new();
                b.read_name(b"read1")
                    .sequence(b"ACGT")
                    .qualities(&[30; 4])
                    .add_int_tag(SamTag::CD, cd)
                    .add_int_tag(SamTag::CM, cm)
                    .add_float_tag(SamTag::CE, ce);
                b.build()
            };
            positional_differs(
                &create_minimal_header("chr1", 10000),
                &[rec(32767, 32767, 0.01)],
                &[rec(59920, 40000, 0.02)],
                ContentPredicate::ExactConsensus,
            )
        },
    },
    CatalogEntry {
        engine: "positional/ExactConsensus",
        mutation: "per-base consensus tag change",
        must_differ: true,
        verify: || {
            positional_differs(
                &create_minimal_header("chr1", 10000),
                &[per_base_consensus_record(&[10, 10, 10, 10])],
                &[per_base_consensus_record(&[10, 10, 9, 10])],
                ContentPredicate::ExactConsensus,
            )
        },
    },
    // Section 2 (duplex): the ExactConsensus saturation carve-out extended to the
    // per-strand aD/bD (32767) and combined cD (65534 = 2 * i16::MAX) duplex depth
    // tags -- Phase-7 finding B.
    CatalogEntry {
        engine: "positional/ExactConsensus",
        mutation: "duplex aD saturated pair tolerated",
        must_differ: false,
        verify: || {
            positional_differs(
                &create_minimal_header("chr1", 10000),
                &[consensus_record_with_tag(SamTag::AD, 32767)],
                &[consensus_record_with_tag(SamTag::AD, 59034)],
                ContentPredicate::ExactConsensus,
            )
        },
    },
    CatalogEntry {
        engine: "positional/ExactConsensus",
        mutation: "duplex bD saturated pair tolerated",
        must_differ: false,
        verify: || {
            positional_differs(
                &create_minimal_header("chr1", 10000),
                &[consensus_record_with_tag(SamTag::BD, 32767)],
                &[consensus_record_with_tag(SamTag::BD, 61000)],
                ContentPredicate::ExactConsensus,
            )
        },
    },
    CatalogEntry {
        engine: "positional/ExactConsensus",
        mutation: "duplex combined cD saturated at 65534 tolerated",
        must_differ: false,
        verify: || {
            positional_differs(
                &create_minimal_header("chr1", 10000),
                &[consensus_record_with_tag(SamTag::CD, 65534)],
                &[consensus_record_with_tag(SamTag::CD, 118_954)],
                ContentPredicate::ExactConsensus,
            )
        },
    },
    CatalogEntry {
        engine: "positional/Exact",
        mutation: "duplex combined cD saturated at 65534 still differs under plain Exact",
        must_differ: true,
        verify: || {
            positional_differs(
                &create_minimal_header("chr1", 10000),
                &[consensus_record_with_tag(SamTag::CD, 65534)],
                &[consensus_record_with_tag(SamTag::CD, 118_954)],
                ContentPredicate::Exact,
            )
        },
    },
    CatalogEntry {
        engine: "positional/ExactConsensus",
        mutation: "duplex aD non-saturated pair still differs",
        must_differ: true,
        verify: || {
            positional_differs(
                &create_minimal_header("chr1", 10000),
                &[consensus_record_with_tag(SamTag::AD, 100)],
                &[consensus_record_with_tag(SamTag::AD, 200)],
                ContentPredicate::ExactConsensus,
            )
        },
    },
    CatalogEntry {
        engine: "positional/ExactConsensus",
        mutation: "duplex aD undercount relative to saturation still differs",
        must_differ: true,
        verify: || {
            positional_differs(
                &create_minimal_header("chr1", 10000),
                &[consensus_record_with_tag(SamTag::AD, 32767)],
                &[consensus_record_with_tag(SamTag::AD, 5)],
                ContentPredicate::ExactConsensus,
            )
        },
    },
    // Section 3: compare_headers
    CatalogEntry {
        engine: "compare_headers",
        mutation: "@SQ name change",
        must_differ: true,
        verify: || headers_differ(&header_with_sq("chr1", 1000), &header_with_sq("chr2", 1000)),
    },
    CatalogEntry {
        engine: "compare_headers",
        mutation: "@SQ length change",
        must_differ: true,
        verify: || headers_differ(&header_with_sq("chr1", 1000), &header_with_sq("chr1", 2000)),
    },
    CatalogEntry {
        engine: "compare_headers",
        mutation: "@SQ order change",
        must_differ: true,
        verify: || {
            headers_differ(
                &header_with_two_sq(("chr1", 1000), ("chr2", 2000)),
                &header_with_two_sq(("chr2", 2000), ("chr1", 1000)),
            )
        },
    },
    CatalogEntry {
        engine: "compare_headers",
        mutation: "@HD SS change",
        must_differ: true,
        verify: || {
            headers_differ(
                &header_with_hd(Some("queryname"), None, Some("queryname:natural")),
                &header_with_hd(Some("queryname"), None, Some("queryname:lexicographical")),
            )
        },
    },
    CatalogEntry {
        engine: "compare_headers",
        mutation: "@RG SM change",
        must_differ: true,
        verify: || {
            headers_differ(
                &header_with_rg_sample("rg1", "sampleA"),
                &header_with_rg_sample("rg1", "sampleB"),
            )
        },
    },
    CatalogEntry {
        engine: "compare_headers",
        mutation: "@PG version change",
        must_differ: false,
        verify: || {
            headers_differ(
                &header_with_pg_version("fgumi", "1.0.0"),
                &header_with_pg_version("fgbio", "4.1.0"),
            )
        },
    },
    CatalogEntry {
        engine: "compare_headers",
        mutation: "@CO change",
        must_differ: false,
        verify: || {
            headers_differ(
                &header_with_comment("built by fgumi"),
                &header_with_comment("built by fgbio"),
            )
        },
    },
    // Section 4: keyjoin (group)
    CatalogEntry {
        engine: "keyjoin",
        mutation: "MI renumber preserving partition",
        must_differ: false,
        verify: || {
            keyjoin_differs(
                &[
                    mi_record(b"read1", 100, "1"),
                    mi_record(b"read2", 200, "1"),
                    mi_record(b"read3", 300, "2"),
                ],
                &[
                    mi_record(b"read1", 100, "5"),
                    mi_record(b"read2", 200, "5"),
                    mi_record(b"read3", 300, "9"),
                ],
            )
        },
    },
    CatalogEntry {
        engine: "keyjoin",
        mutation: "MI bijection violation",
        must_differ: true,
        verify: || {
            keyjoin_differs(
                &[mi_record(b"read1", 100, "1"), mi_record(b"read2", 200, "1")],
                &[mi_record(b"read1", 100, "1"), mi_record(b"read2", 200, "2")],
            )
        },
    },
    CatalogEntry {
        engine: "keyjoin",
        mutation: "reordered records same grouping",
        must_differ: false,
        verify: || {
            let r1 = mi_record(b"read1", 100, "1");
            let r2 = mi_record(b"read2", 200, "1");
            let r3 = mi_record(b"read3", 300, "2");
            keyjoin_differs(&[r1.clone(), r2.clone(), r3.clone()], &[r3, r1, r2])
        },
    },
    CatalogEntry {
        engine: "keyjoin",
        mutation: "content diff with intact grouping",
        must_differ: true,
        verify: || {
            let good = |name: &[u8], pos: i32, seq: &[u8]| {
                let mut b = SamBuilder::new();
                b.read_name(name)
                    .sequence(seq)
                    .qualities(&vec![30; seq.len()])
                    .ref_id(0)
                    .pos(pos - 1)
                    .mapq(60)
                    .add_string_tag(SamTag::MI, b"1");
                b.build()
            };
            keyjoin_differs(
                &[good(b"read1", 100, b"ACGTACGT"), good(b"read2", 200, b"ACGTACGT")],
                &[good(b"read1", 100, b"ACGTACGT"), good(b"read2", 200, b"ACGTACGA")],
            )
        },
    },
    CatalogEntry {
        engine: "keyjoin",
        mutation: "dropped record presence diff",
        must_differ: true,
        verify: || {
            keyjoin_differs(
                &[mi_record(b"read1", 100, "1"), mi_record(b"read2", 200, "1")],
                &[mi_record(b"read1", 100, "1")],
            )
        },
    },
    // Section 5: sort_verify_compare
    CatalogEntry {
        engine: "sort_verify",
        mutation: "sort order violation",
        must_differ: true,
        verify: || {
            sort_verify_differs(
                &[coord_record(b"r1", 50), coord_record(b"r2", 100), coord_record(b"r3", 150)],
                &[coord_record(b"r1", 50), coord_record(b"r2", 100), coord_record(b"r3", 75)],
            )
        },
    },
    CatalogEntry {
        engine: "sort_verify",
        mutation: "swap within equal sort-key run",
        must_differ: false,
        verify: || {
            let a = coord_record(b"readA", 100);
            let b_rec = coord_record(b"readB", 100);
            sort_verify_differs(&[a.clone(), b_rec.clone()], &[b_rec, a])
        },
    },
    CatalogEntry {
        engine: "sort_verify",
        mutation: "dropped record within equal sort-key run",
        must_differ: true,
        verify: || {
            let a = coord_record(b"readA", 100);
            let b_rec = coord_record(b"readB", 100);
            sort_verify_differs(&[a.clone(), b_rec], &[a])
        },
    },
    // Section 6: metrics
    CatalogEntry {
        engine: "metrics",
        mutation: "integer count change",
        must_differ: true,
        verify: || metrics_content_differs("umi\tcount\nAAA\t3\n", "umi\tcount\nAAA\t4\n"),
    },
    CatalogEntry {
        engine: "metrics",
        mutation: "float beyond epsilon change",
        must_differ: true,
        verify: || metrics_content_differs("key\tvalue\nA\t1.0\n", "key\tvalue\nA\t1.5\n"),
    },
    CatalogEntry {
        engine: "metrics",
        mutation: "float within epsilon representation only",
        must_differ: false,
        verify: || {
            metrics_content_differs(
                "key\tvalue\nA\t0.333333\n",
                "key\tvalue\nA\t0.3333333333333333\n",
            )
        },
    },
    CatalogEntry {
        engine: "metrics",
        mutation: "dropped row change",
        must_differ: true,
        verify: || metrics_content_differs("key\tvalue\nA\t1\nB\t2\n", "key\tvalue\nA\t1\n"),
    },
    CatalogEntry {
        engine: "metrics",
        mutation: "added row change",
        must_differ: true,
        verify: || metrics_content_differs("key\tvalue\nA\t1\n", "key\tvalue\nA\t1\nB\t2\n"),
    },
    CatalogEntry {
        engine: "metrics",
        mutation: "permuted row order",
        must_differ: false,
        verify: || {
            metrics_content_differs(
                "key\tvalue\nA\t1\nB\t2\nC\t3\n",
                "key\tvalue\nC\t3\nA\t1\nB\t2\n",
            )
        },
    },
    CatalogEntry {
        engine: "metrics",
        mutation: "different column set",
        must_differ: true,
        verify: || {
            metrics_content_differs("umi\tcount\tfraction\nAAA\t3\t0.5\n", "umi\tcount\nAAA\t3\n")
        },
    },
];

#[test]
fn mutation_catalog_covers_every_compared_dimension() {
    // (1) Every "dimension" the compare-hardening design spec's mutation catalog
    // calls out (core SAM fields, tag classes, header record types, metrics row
    // dimensions, sort-key cases) must be named by at least one `MUTATION_CATALOG`
    // entry's `mutation` string -- this catches a future compared field/tag/
    // row-dimension added without any catalog entry at all.
    const REQUIRED_SUBSTRINGS: &[&str] = &[
        // core SAM fields
        "SEQ base",
        "QUAL byte",
        "POS shift",
        "FLAG bit",
        "CIGAR change",
        "MAPQ change",
        "RNEXT change",
        "PNEXT change",
        "TLEN change",
        "RNAME change",
        "QNAME",
        // tag classes
        "tag drop",
        "tag value change",
        "array-tag element",
        "tag reorder",
        "per-base consensus tag",
        "MI value change",
        "MI renumber",
        "MI bijection",
        // records / order
        "drop a record",
        "duplicate a record",
        "swap two distinct-key records",
        // consensus saturation carve-out
        "cD saturated pair tolerated",
        "cD non-saturated pair",
        "duplex aD saturated pair tolerated",
        "duplex bD saturated pair tolerated",
        "duplex combined cD saturated at 65534 tolerated",
        "duplex aD non-saturated pair",
        "duplex aD undercount relative to saturation",
        // header record types
        "@SQ name",
        "@SQ length",
        "@SQ order",
        "@HD SS",
        "@RG SM",
        "@PG version",
        "@CO change",
        // metrics row dimensions
        "integer count change",
        "float beyond epsilon",
        "float within epsilon",
        "dropped row change",
        "added row change",
        "permuted row order",
        "different column set",
        // sort-key cases
        "sort order violation",
        "swap within equal sort-key run",
        "dropped record within equal sort-key run",
    ];

    for needle in REQUIRED_SUBSTRINGS {
        assert!(
            MUTATION_CATALOG.iter().any(|entry| entry.mutation.contains(needle)),
            "no mutation catalog entry covers required dimension {needle:?} -- a compared \
             field/tag/row-dimension may be missing its mutation test"
        );
    }

    // (2) The load-bearing check: re-run every entry's `verify` closure against the
    // real comparator and confirm it agrees with the entry's claimed `must_differ`.
    // Unlike (1), this does not trust the catalog's own text -- it trusts only the
    // production `positional_compare`/`compare_headers`/`keyjoin_compare`/
    // `sort_verify_compare`/`CompareMetrics` entry points. A catalog row can no
    // longer "silently pass" merely by existing: if the compared dimension it
    // names stops actually being compared (or never was), this loop fails.
    for entry in MUTATION_CATALOG {
        let differed = (entry.verify)();
        assert_eq!(
            differed, entry.must_differ,
            "catalog entry {:?}/{:?} claims must_differ={}, but re-running it against the real \
             comparator produced DIFFER={differed}",
            entry.engine, entry.mutation, entry.must_differ
        );
    }

    let must_differ_count = MUTATION_CATALOG.iter().filter(|e| e.must_differ).count();
    let tolerated_count = MUTATION_CATALOG.len() - must_differ_count;
    assert!(
        MUTATION_CATALOG.len() >= 40,
        "mutation catalog looks too small ({} entries) for the design spec's catalog",
        MUTATION_CATALOG.len()
    );
    assert!(
        must_differ_count >= 25,
        "expected a majority of must-DIFFER cases: {must_differ_count}"
    );
    assert!(tolerated_count >= 8, "expected several accepted-divergence cases: {tolerated_count}");

    // Every catalog entry names a non-empty engine -- guards against a copy-paste
    // row with a blank/placeholder engine name slipping in unnoticed.
    assert!(MUTATION_CATALOG.iter().all(|e| !e.engine.is_empty()));
}

//! Content-equality predicates for a single pair of BAM records.
//!
//! A [`ContentPredicate`] decides whether two BAM records are "the same" for comparison
//! purposes; [`content_diffs`] evaluates a predicate against a pair of raw BAM records
//! and, when they differ, returns human-readable diff strings describing exactly which
//! core field(s) or tag(s) disagree.
//!
//! This module intentionally has no notion of record *position* within a stream — that is
//! the job of the positional engine (`super::positional`). Keeping the two concerns
//! separate is what lets the positional engine treat any key/position/length mismatch as
//! an immediate DIFFER without ever re-interpreting content equality as license to resync.

use std::collections::BTreeMap;

use fgumi_raw_bam::{RawRecordView, TagValue};
use noodles::sam::Header;

use crate::sam::SamTag;

use super::super::bams::{FIELD_NAMES, get_core_fields_raw};
use super::super::raw_compare::{
    raw_compare_structured, raw_core_fields_equal, raw_records_byte_equal, raw_tags_byte_equal,
};

/// A predicate for deciding whether two BAM records are content-equal.
///
/// Content predicates ignore each record's position within its stream (see
/// [`super::positional`] for the ordering-aware engine) and compare only the byte-level
/// content of a single pair of records.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ContentPredicate {
    /// Every core SAM field must match exactly (via [`raw_compare_structured`]'s
    /// `core_match`), and every aux tag must match under an order-independent, semantic
    /// comparison of tag values (`core_match`'s sibling `tag_order_match`, computed from
    /// [`super::super::raw_compare::raw_tags_equal_order_independent`]). Tag order and
    /// on-disk integer width are *not* significant — only the decoded values are.
    ///
    /// One accepted divergence applies to *every* predicate (not just `Exact`): aux tags
    /// that fgumi emits but fgbio never persists — currently just `tc`, the
    /// template-coordinate sort key `fgumi zipper` writes and `fgumi dedup` consumes — are
    /// ignored entirely (see [`super::super::raw_compare::is_ignored_fgumi_only_tag`] and
    /// the design spec §"Accepted divergences"). The carve-out is narrow: only that tag is
    /// dropped, so a real difference in any other tag still reports a DIFFER.
    Exact,
    /// Consensus BAM output (`simplex`/`duplex`/`codec`, and `filter` output run over
    /// consensus reads). This compares **exactly like [`Self::Exact`]** — including the
    /// per-read depth/disagreement tags `cD`/`cM`/`aD`/`aM`/`bD`/`bM` and their `cE`/`aE`/`bE`
    /// error rates.
    ///
    /// Historically this predicate carried a saturation-aware carve-out that tolerated a
    /// depth divergence: fgbio stores per-base consensus depth as a `Short` `Array` capped
    /// at `i16::MAX` (32767) and derives the scalar depth tags from those capped arrays,
    /// while fgumi's codec/duplex callers emitted uncapped depths. That divergence has been
    /// fixed at the source — fgumi now clamps every scalar depth tag to fgbio's `Short`
    /// ceiling (per-strand `aD`/`bD` at 32767; combined `cD` caps each strand per base
    /// *before* summing, matching fgbio's `totalDepths`) — so the tags are bit-identical and
    /// the tolerance is no longer needed. Comparing them exactly is now both sound and
    /// complete; the old carve-out could not be made sound for the combined `cD` anyway
    /// (fgbio caps-then-sums, which is unrecoverable from the scalar).
    ///
    /// Retained as a distinct predicate for the consensus/`filter` presets; now that the
    /// carve-out is gone it is behaviourally equivalent to [`Self::Exact`] and could be
    /// consolidated with it in a follow-up.
    ExactConsensus,
    /// [`Self::Exact`], but excludes the `MI` tag entirely from tag comparison — neither its
    /// *value* nor its *presence* is compared. `group` legitimately renumbers MI values across
    /// tools (compare-hardening design spec, §"Accepted divergences", "`group` MI numbering"),
    /// and the key-join engine (`super::keyjoin`) verifies MI equivalence separately via a
    /// bijection over the joined records' MI values, plus its own "missing MI" accounting — so
    /// this predicate must not also react to an MI value or presence difference, which would
    /// double up (and could disagree) with those separate checks. Every other tag and every
    /// core SAM field stay exact.
    ExactMinusMi,
}

/// Compares two raw BAM records under `pred`.
///
/// Returns `None` if the records are content-equal under `pred`, or
/// `Some(diffs)` with one human-readable string per mismatched core field or aux tag
/// otherwise.
///
/// `header_a` and `header_b` are used *only* to resolve RNAME/RNEXT reference names for
/// `a` and `b` respectively when rendering diff strings (mirroring
/// [`get_core_fields_raw`]'s per-side header usage); the equality decision itself
/// operates on raw `ref_id` bytes and never consults either header.
#[must_use]
pub fn content_diffs(
    a: &[u8],
    b: &[u8],
    pred: ContentPredicate,
    header_a: &Header,
    header_b: &Header,
) -> Option<Vec<String>> {
    // Byte-identical records are content-equal under every `ContentPredicate`: `Exact`
    // compares them exactly, `ExactConsensus` is now equivalent to `Exact`, and
    // `ExactMinusMi` only *drops* the MI check (never adds strictness), so a byte-for-byte
    // match can never be a `Some` under any of the three. Checking this first skips building
    // any tag maps on the common (matching) case.
    if raw_records_byte_equal(a, b) {
        return None;
    }

    match pred {
        ContentPredicate::Exact => exact_diffs(a, b, header_a, header_b),
        ContentPredicate::ExactConsensus => exact_consensus_diffs(a, b, header_a, header_b),
        ContentPredicate::ExactMinusMi => exact_minus_mi_diffs(a, b, header_a, header_b),
    }
}

/// Shared "collect diffs, then wrap" shape used by every [`ContentPredicate`] implementation:
/// given the two boolean match outcomes, return `None` if both matched, otherwise lazily
/// invoke whichever of `core`/`tags` corresponds to a `false` outcome, concatenate their diff
/// strings, and fall back to a single "could not be localized" note if that concatenation is
/// still empty (e.g. malformed data that only the underlying byte-level comparison detects) —
/// a real DIFFER must never silently report zero diff lines.
///
/// `core` and `tags` are `FnOnce` closures (not plain `Vec<String>` arguments) so the caller
/// only pays for building a diff-string `Vec` on the side that actually mismatched.
fn finish_diffs(
    core_match: bool,
    tags_match: bool,
    core: impl FnOnce() -> Vec<String>,
    tags: impl FnOnce() -> Vec<String>,
) -> Option<Vec<String>> {
    if core_match && tags_match {
        return None;
    }

    let mut diffs = Vec::new();
    if !core_match {
        diffs.extend(core());
    }
    if !tags_match {
        diffs.extend(tags());
    }
    if diffs.is_empty() {
        diffs.push("records differ but the mismatch could not be localized".to_string());
    }
    Some(diffs)
}

/// Implements [`ContentPredicate::Exact`].
fn exact_diffs(a: &[u8], b: &[u8], header_a: &Header, header_b: &Header) -> Option<Vec<String>> {
    let result = raw_compare_structured(a, b);
    finish_diffs(
        result.core_match,
        result.tag_order_match,
        || core_field_diffs(a, b, header_a, header_b),
        || tag_diffs(a, b),
    )
}

/// Implements [`ContentPredicate::ExactConsensus`] (see its doc-comment for the exact
/// carve-out rules and the accepted-divergence citation).
fn exact_consensus_diffs(
    a: &[u8],
    b: &[u8],
    header_a: &Header,
    header_b: &Header,
) -> Option<Vec<String>> {
    let core_match = raw_core_fields_equal(a, b);
    let tags_match = consensus_tags_match(a, b);
    finish_diffs(
        core_match,
        tags_match,
        || core_field_diffs(a, b, header_a, header_b),
        || consensus_tag_diffs(a, b),
    )
}

/// Implements [`ContentPredicate::ExactMinusMi`] (see its doc-comment).
fn exact_minus_mi_diffs(
    a: &[u8],
    b: &[u8],
    header_a: &Header,
    header_b: &Header,
) -> Option<Vec<String>> {
    let core_match = raw_core_fields_equal(a, b);
    let tags_match = tags_match_excluding_mi(a, b);
    finish_diffs(
        core_match,
        tags_match,
        || core_field_diffs(a, b, header_a, header_b),
        || tag_diffs_excluding_mi(a, b),
    )
}

/// [`tag_typed_map`], with the `MI` tag removed (if present).
fn tag_typed_map_excluding_mi(raw: &[u8]) -> BTreeMap<[u8; 2], TagValue<'_>> {
    let mut tags = tag_typed_map(raw);
    tags.remove(&*SamTag::MI);
    tags
}

/// Returns `true` if `a` and `b` have the same tag set and every tag value is equal,
/// ignoring `MI` entirely (on either side, in any combination of present/absent/differing).
fn tags_match_excluding_mi(a: &[u8], b: &[u8]) -> bool {
    // Byte-identical tag regions are trivially equal under every comparison this function
    // could otherwise perform, so this skips building either side's tag map on the common
    // (matching) case — mirroring `ContentPredicate::Exact`'s fast path in `raw_compare_structured`.
    if raw_tags_byte_equal(a, b) {
        return true;
    }
    let tags_a = tag_typed_map_excluding_mi(a);
    let tags_b = tag_typed_map_excluding_mi(b);
    tags_a.len() == tags_b.len() && tags_a.iter().all(|(tag, va)| tags_b.get(tag) == Some(va))
}

/// Shared "union tag names, sort+dedup, `filter_map` a `{tag}: {va:?} vs {vb:?}` line for
/// unequal pairs" shape used by every per-predicate tag-diff builder ([`tag_diffs_excluding_mi`],
/// [`consensus_tag_diffs`], [`tag_diffs`]) — these differ only in how a single tag's pair of
/// (possibly absent) values is judged equal.
///
/// `values_equal` receives the tag name plus each side's value (`None` if that tag is absent
/// on that side, which only happens when the tag is present on the other side, since `tag`
/// always comes from the union of both maps' keys) and decides whether the tag counts as equal.
fn tag_diffs_with(
    tags_a: &BTreeMap<[u8; 2], TagValue<'_>>,
    tags_b: &BTreeMap<[u8; 2], TagValue<'_>>,
    values_equal: impl Fn([u8; 2], Option<&TagValue<'_>>, Option<&TagValue<'_>>) -> bool,
) -> Vec<String> {
    let mut all_tags: Vec<&[u8; 2]> = tags_a.keys().chain(tags_b.keys()).collect();
    all_tags.sort();
    all_tags.dedup();
    all_tags
        .into_iter()
        .filter_map(|tag| {
            let (va, vb) = (tags_a.get(tag), tags_b.get(tag));
            (!values_equal(*tag, va, vb))
                .then(|| format!("{}: {va:?} vs {vb:?}", String::from_utf8_lossy(tag)))
        })
        .collect()
}

/// Builds one diff string per aux tag (excluding `MI`) present on only one side, or whose
/// value differs, between `a` and `b`.
fn tag_diffs_excluding_mi(a: &[u8], b: &[u8]) -> Vec<String> {
    let tags_a = tag_typed_map_excluding_mi(a);
    let tags_b = tag_typed_map_excluding_mi(b);
    tag_diffs_with(&tags_a, &tags_b, |_tag, va, vb| va == vb)
}

/// Extracts a record's aux tags as a map of tag name -> zero-copy typed value.
///
/// Values are kept as [`TagValue`] (not stringified): `TagValue` derives `PartialEq` and
/// widens integer tags to `i64` regardless of on-disk width (`c`/`C`/`s`/`S`/`i`/`I`), so
/// comparing `TagValue`s directly already gives the same width-insensitive equality that
/// [`super::super::raw_compare::raw_tags_equal_order_independent`] uses to decide
/// `tag_order_match` — no separate stringified map is needed for [`tag_diffs`] or
/// [`consensus_tag_diffs`].
fn tag_typed_map(raw: &[u8]) -> BTreeMap<[u8; 2], TagValue<'_>> {
    let mut tags: BTreeMap<[u8; 2], TagValue<'_>> =
        RawRecordView::new(raw).tags().iter_typed().collect();
    // Drop fgumi-only tags (e.g. `tc`) that fgbio never persists, so they are ignored — both
    // presence and value — by every predicate that builds its comparison off this map
    // (`ExactConsensus`, `ExactMinusMi`, and all diff-string rendering). `Exact`'s match
    // decision applies the same filter in `raw_tags_equal_order_independent`.
    tags.retain(|tag, _| !super::super::raw_compare::is_ignored_fgumi_only_tag(*tag));
    tags
}

/// Returns `true` if `a` and `b` have the same tag set and every tag value is equal
/// (exactly). Consensus depth/error tags are compared exactly like every other tag now
/// that fgumi clamps them to fgbio's `Short` ceiling (see [`ContentPredicate::ExactConsensus`]).
fn consensus_tags_match(a: &[u8], b: &[u8]) -> bool {
    // Byte-identical tag regions are trivially equal under every comparison this function
    // could otherwise perform, so this skips building either side's tag map on the common
    // (matching) case — mirroring `ContentPredicate::Exact`'s fast path in `raw_compare_structured`.
    if raw_tags_byte_equal(a, b) {
        return true;
    }
    let tags_a = tag_typed_map(a);
    let tags_b = tag_typed_map(b);
    if tags_a.len() != tags_b.len() {
        return false;
    }
    tags_a.iter().all(|(tag, va)| tags_b.get(tag).is_some_and(|vb| va == vb))
}

/// Builds one diff string per aux tag that is present on only one side, or whose value
/// differs (exactly), between `a` and `b`.
fn consensus_tag_diffs(a: &[u8], b: &[u8]) -> Vec<String> {
    let tags_a = tag_typed_map(a);
    let tags_b = tag_typed_map(b);
    tag_diffs_with(&tags_a, &tags_b, |_tag, va, vb| va == vb)
}

/// Builds one diff string per core SAM field (see [`FIELD_NAMES`]) whose rendered value
/// differs between `a` and `b`.
fn core_field_diffs(a: &[u8], b: &[u8], header_a: &Header, header_b: &Header) -> Vec<String> {
    let fields_a = get_core_fields_raw(a, header_a);
    let fields_b = get_core_fields_raw(b, header_b);
    fields_a
        .iter()
        .zip(fields_b.iter())
        .zip(FIELD_NAMES.iter())
        .filter(|((va, vb), _)| va != vb)
        .map(|((va, vb), name)| format!("{name}: {va:?} vs {vb:?}"))
        .collect()
}

/// Builds one diff string per aux tag that is present on only one side, or whose value
/// differs, between `a` and `b`. Used by [`ContentPredicate::Exact`] (via [`exact_diffs`]),
/// which requires every tag value to match exactly — unlike [`consensus_tag_diffs`], no
/// carve-out is applied.
fn tag_diffs(a: &[u8], b: &[u8]) -> Vec<String> {
    let tags_a = tag_typed_map(a);
    let tags_b = tag_typed_map(b);
    tag_diffs_with(&tags_a, &tags_b, |_tag, va, vb| va == vb)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fgumi_raw_bam::{RawRecord, SamBuilder};
    use noodles::sam::Header;
    use rstest::rstest;

    fn hdr() -> Header {
        Header::default()
    }

    #[test]
    fn identical_records_have_no_diffs() {
        let a = SamBuilder::new().read_name(b"r").sequence(b"ACGT").build();
        let b = SamBuilder::new().read_name(b"r").sequence(b"ACGT").build();
        assert!(
            content_diffs(a.as_ref(), b.as_ref(), ContentPredicate::Exact, &hdr(), &hdr())
                .is_none()
        );
    }

    #[test]
    fn seq_difference_is_reported() {
        let a = SamBuilder::new().read_name(b"r").sequence(b"ACGT").build();
        let b = SamBuilder::new().read_name(b"r").sequence(b"ACGA").build();
        let diffs = content_diffs(a.as_ref(), b.as_ref(), ContentPredicate::Exact, &hdr(), &hdr())
            .expect("a single SEQ base difference must be reported");
        assert!(
            diffs.iter().any(|d| d.starts_with("SEQ:")),
            "diff must name the SEQ field, got: {diffs:?}"
        );
    }

    #[test]
    fn qual_difference_is_reported() {
        let a = SamBuilder::new()
            .read_name(b"r")
            .sequence(b"ACGT")
            .qualities(&[30, 30, 30, 30])
            .build();
        let b = SamBuilder::new()
            .read_name(b"r")
            .sequence(b"ACGT")
            .qualities(&[30, 30, 30, 10])
            .build();
        let diffs = content_diffs(a.as_ref(), b.as_ref(), ContentPredicate::Exact, &hdr(), &hdr())
            .expect("a single QUAL base difference must be reported");
        assert!(
            diffs.iter().any(|d| d.starts_with("QUAL:")),
            "diff must name the QUAL field, got: {diffs:?}"
        );
    }

    #[test]
    fn tag_difference_is_reported() {
        let a =
            SamBuilder::new().read_name(b"r").sequence(b"ACGT").add_int_tag(SamTag::NM, 1).build();
        let b =
            SamBuilder::new().read_name(b"r").sequence(b"ACGT").add_int_tag(SamTag::NM, 2).build();
        let diffs = content_diffs(a.as_ref(), b.as_ref(), ContentPredicate::Exact, &hdr(), &hdr())
            .expect("a differing tag value must be reported");
        assert!(
            diffs.iter().any(|d| d.starts_with("NM:")),
            "diff must name the differing NM tag, got: {diffs:?}"
        );
    }

    #[test]
    fn tag_order_is_not_a_difference() {
        // Same two tags, different on-record order: `Exact` compares tags
        // order-independently, so this must NOT be reported as a diff.
        let a = SamBuilder::new()
            .read_name(b"r")
            .sequence(b"ACGT")
            .add_int_tag(SamTag::NM, 1)
            .add_string_tag(SamTag::RG, b"grp")
            .build();
        let b = SamBuilder::new()
            .read_name(b"r")
            .sequence(b"ACGT")
            .add_string_tag(SamTag::RG, b"grp")
            .add_int_tag(SamTag::NM, 1)
            .build();
        assert!(
            content_diffs(a.as_ref(), b.as_ref(), ContentPredicate::Exact, &hdr(), &hdr())
                .is_none(),
            "tag order must not affect Exact equality"
        );
    }

    // ==================== ExactConsensus: exact depth/error tags ====================
    //
    // fgumi clamps every scalar consensus depth tag to fgbio's `Short` ceiling (see
    // `nh/fix-consensus-depth-clamp-parity`), so `cD`/`cM`/`aD`/`aM`/`bD`/`bM` and their
    // `cE`/`aE`/`bE` error rates are bit-identical to fgbio. `ExactConsensus` therefore
    // compares them EXACTLY, with no saturation tolerance: any difference — including a
    // pre-fix "saturated" pair like `32767` vs `59920`, or the duplex `65534` sum — is a
    // real divergence and must DIFFER.

    fn depth_tag_record(tag: SamTag, value: i32) -> RawRecord {
        let mut b = SamBuilder::new();
        b.read_name(b"r").sequence(b"ACGT").add_int_tag(tag, value);
        b.build()
    }

    /// Every consensus depth tag is compared exactly under both `ExactConsensus` and plain
    /// `Exact`: a former "saturation" pair now DIFFERs, and an identical pair matches.
    #[rstest]
    #[case::cd(SamTag::CD)]
    #[case::cm(SamTag::CM)]
    #[case::ad(SamTag::AD)]
    #[case::am(SamTag::AM)]
    #[case::bd(SamTag::BD)]
    #[case::bm(SamTag::BM)]
    fn exact_consensus_compares_depth_tags_exactly(#[case] tag: SamTag) {
        let saturated = depth_tag_record(tag, 32767);
        let uncapped = depth_tag_record(tag, 59920);
        for pred in [ContentPredicate::ExactConsensus, ContentPredicate::Exact] {
            assert!(
                content_diffs(saturated.as_ref(), uncapped.as_ref(), pred, &hdr(), &hdr())
                    .is_some_and(|d| !d.is_empty()),
                "{tag:?} 32767 vs 59920 must DIFFER under {pred:?} (no saturation tolerance)"
            );
        }
        let same = depth_tag_record(tag, 32767);
        assert!(
            content_diffs(
                saturated.as_ref(),
                same.as_ref(),
                ContentPredicate::ExactConsensus,
                &hdr(),
                &hdr()
            )
            .is_none(),
            "identical {tag:?} must match under ExactConsensus"
        );
    }

    /// Consensus error tags (`cE`/`aE`/`bE`) are also compared exactly: any float
    /// difference DIFFERs, with no gating on depth saturation.
    #[rstest]
    #[case::ce(SamTag::CE)]
    #[case::ae(SamTag::AE)]
    #[case::be(SamTag::BE)]
    fn exact_consensus_compares_error_tags_exactly(#[case] tag: SamTag) {
        let a =
            SamBuilder::new().read_name(b"r").sequence(b"ACGT").add_float_tag(tag, 0.01).build();
        let b =
            SamBuilder::new().read_name(b"r").sequence(b"ACGT").add_float_tag(tag, 0.02).build();
        assert!(
            content_diffs(a.as_ref(), b.as_ref(), ContentPredicate::ExactConsensus, &hdr(), &hdr())
                .is_some_and(|d| !d.is_empty()),
            "{tag:?} error-rate difference must DIFFER"
        );
        let same =
            SamBuilder::new().read_name(b"r").sequence(b"ACGT").add_float_tag(tag, 0.01).build();
        assert!(
            content_diffs(
                a.as_ref(),
                same.as_ref(),
                ContentPredicate::ExactConsensus,
                &hdr(),
                &hdr()
            )
            .is_none(),
            "identical {tag:?} must match"
        );
    }

    #[test]
    fn exact_consensus_still_exact_on_unrelated_tags() {
        // A non-consensus tag (NM) must still be compared exactly under ExactConsensus.
        let a =
            SamBuilder::new().read_name(b"r").sequence(b"ACGT").add_int_tag(SamTag::NM, 1).build();
        let b =
            SamBuilder::new().read_name(b"r").sequence(b"ACGT").add_int_tag(SamTag::NM, 2).build();
        let diffs =
            content_diffs(a.as_ref(), b.as_ref(), ContentPredicate::ExactConsensus, &hdr(), &hdr())
                .expect("a differing NM tag must be reported under ExactConsensus");
        assert!(
            diffs.iter().any(|d| d.starts_with("NM:")),
            "diff must name the differing NM tag, got: {diffs:?}"
        );
    }

    // ==================== ExactMinusMi carve-out tests ====================
    //
    // `ExactMinusMi` excludes the `MI` tag entirely (value and presence) from the
    // comparison, since `group` legitimately renumbers MI and the key-join engine
    // verifies MI equivalence separately via its own bijection/missing-MI accounting
    // (see `ContentPredicate::ExactMinusMi`'s doc-comment).

    #[rstest]
    #[case::mi_value_differs_only(5, 99, 1, 1, false)] // MI differs, NM same -> not a diff
    #[case::mi_and_nm_both_differ(5, 99, 1, 2, true)] // MI differs, NM also differs -> diff (NM only)
    #[case::identical(5, 5, 1, 1, false)]
    fn exact_minus_mi_ignores_mi_value(
        #[case] mi_a: i32,
        #[case] mi_b: i32,
        #[case] nm_a: i32,
        #[case] nm_b: i32,
        #[case] expect_diff: bool,
    ) {
        let a = SamBuilder::new()
            .read_name(b"r")
            .sequence(b"ACGT")
            .add_int_tag(SamTag::MI, mi_a)
            .add_int_tag(SamTag::NM, nm_a)
            .build();
        let b = SamBuilder::new()
            .read_name(b"r")
            .sequence(b"ACGT")
            .add_int_tag(SamTag::MI, mi_b)
            .add_int_tag(SamTag::NM, nm_b)
            .build();
        let diffs =
            content_diffs(a.as_ref(), b.as_ref(), ContentPredicate::ExactMinusMi, &hdr(), &hdr());
        assert_eq!(diffs.is_some(), expect_diff, "diffs={diffs:?}");
        if let Some(diffs) = diffs {
            assert!(
                diffs.iter().all(|d| !d.starts_with("MI:")),
                "MI must never be reported under ExactMinusMi, got: {diffs:?}"
            );
        }
    }

    #[test]
    fn exact_minus_mi_ignores_mi_presence_difference() {
        // MI present on `a`, absent on `b` — still not a diff under ExactMinusMi (the
        // key-join engine's separate "missing MI" accounting owns this case).
        let a =
            SamBuilder::new().read_name(b"r").sequence(b"ACGT").add_int_tag(SamTag::MI, 5).build();
        let b = SamBuilder::new().read_name(b"r").sequence(b"ACGT").build();
        assert!(
            content_diffs(a.as_ref(), b.as_ref(), ContentPredicate::ExactMinusMi, &hdr(), &hdr())
                .is_none(),
            "MI presence-only difference must not be reported under ExactMinusMi"
        );
    }

    #[test]
    fn exact_minus_mi_still_flags_core_field_differences() {
        let a = SamBuilder::new().read_name(b"r").sequence(b"ACGT").build();
        let b = SamBuilder::new().read_name(b"r").sequence(b"ACGA").build();
        let diffs =
            content_diffs(a.as_ref(), b.as_ref(), ContentPredicate::ExactMinusMi, &hdr(), &hdr())
                .expect("a SEQ difference must still be reported under ExactMinusMi");
        assert!(diffs.iter().any(|d| d.starts_with("SEQ:")));
    }

    #[test]
    fn exact_minus_mi_still_flags_non_mi_tag_differences() {
        let a =
            SamBuilder::new().read_name(b"r").sequence(b"ACGT").add_int_tag(SamTag::NM, 1).build();
        let b =
            SamBuilder::new().read_name(b"r").sequence(b"ACGT").add_int_tag(SamTag::NM, 2).build();
        let diffs =
            content_diffs(a.as_ref(), b.as_ref(), ContentPredicate::ExactMinusMi, &hdr(), &hdr())
                .expect("a differing non-MI tag must still be reported under ExactMinusMi");
        assert!(
            diffs.iter().any(|d| d.starts_with("NM:")),
            "diff must name the differing NM tag, got: {diffs:?}"
        );
    }

    #[test]
    fn exact_predicate_still_flags_mi_value_difference() {
        // Sanity check that plain Exact is unaffected by ExactMinusMi's carve-out: an
        // MI-only difference IS a diff under Exact.
        let a =
            SamBuilder::new().read_name(b"r").sequence(b"ACGT").add_int_tag(SamTag::MI, 5).build();
        let b =
            SamBuilder::new().read_name(b"r").sequence(b"ACGT").add_int_tag(SamTag::MI, 99).build();
        let diffs = content_diffs(a.as_ref(), b.as_ref(), ContentPredicate::Exact, &hdr(), &hdr())
            .expect("an MI value difference must be reported under plain Exact");
        assert!(
            diffs.iter().any(|d| d.starts_with("MI:")),
            "diff must name the differing MI tag, got: {diffs:?}"
        );
    }

    #[rstest]
    #[case::exact(ContentPredicate::Exact)]
    #[case::exact_consensus(ContentPredicate::ExactConsensus)]
    #[case::exact_minus_mi(ContentPredicate::ExactMinusMi)]
    fn tc_tag_only_on_fgumi_side_is_tolerated(#[case] pred: ContentPredicate) {
        // fgumi carries the fgumi-only `tc` tag; the fgbio-style record does not. Every
        // predicate must treat this as a MATCH (fgbio never persists `tc`).
        let fgumi = SamBuilder::new()
            .read_name(b"r")
            .sequence(b"ACGT")
            .add_int_tag(SamTag::NM, 1)
            .add_int_tag(SamTag::TC, 99)
            .build();
        let fgbio =
            SamBuilder::new().read_name(b"r").sequence(b"ACGT").add_int_tag(SamTag::NM, 1).build();
        assert!(
            content_diffs(fgumi.as_ref(), fgbio.as_ref(), pred, &hdr(), &hdr()).is_none(),
            "an fgumi-only `tc` tag must be ignored under {pred:?}"
        );
    }

    #[rstest]
    #[case::exact(ContentPredicate::Exact)]
    #[case::exact_consensus(ContentPredicate::ExactConsensus)]
    #[case::exact_minus_mi(ContentPredicate::ExactMinusMi)]
    fn tc_carveout_is_narrow_other_tag_still_differs(#[case] pred: ContentPredicate) {
        // The carve-out drops only `tc`: a genuine difference in another tag (NM) still
        // DIFFERs even when an ignored `tc` is present on the fgumi side.
        let fgumi = SamBuilder::new()
            .read_name(b"r")
            .sequence(b"ACGT")
            .add_int_tag(SamTag::NM, 1)
            .add_int_tag(SamTag::TC, 99)
            .build();
        let fgbio =
            SamBuilder::new().read_name(b"r").sequence(b"ACGT").add_int_tag(SamTag::NM, 2).build();
        let diffs = content_diffs(fgumi.as_ref(), fgbio.as_ref(), pred, &hdr(), &hdr())
            .expect("a real NM difference must still be reported despite the ignored tc");
        assert!(
            diffs.iter().any(|d| d.starts_with("NM:")),
            "diff must name the differing NM tag, got: {diffs:?}"
        );
    }
}

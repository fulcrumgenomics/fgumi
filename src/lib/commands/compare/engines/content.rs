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
    /// [`Self::Exact`], plus a single accepted divergence for consensus BAM output
    /// (`simplex`/`duplex`/`codec`, and `filter` output run over consensus reads, since
    /// `filter` only drops reads and passes these tags through unchanged): fgbio stores
    /// per-read consensus depth and disagreement-count tags as a signed 16-bit `Short`
    /// and saturates at `i16::MAX` (32767) for families deeper than that, while fgumi
    /// reports the true (larger) depth. fgumi's uncapped depth is the more correct value,
    /// so this divergence is accepted rather than fixed (see the compare-hardening design
    /// spec, §"Accepted divergences", "Consensus per-read depth saturation"; decision-b,
    /// 2026-07-08).
    ///
    /// This applies to both simplex/codec (`cD`/`cM`/`cE` only) and duplex output,
    /// where fgbio *also* saturates the per-strand `aD`/`aM` and `bD`/`bM` at 32767
    /// each, and the *combined* `cD`/`cM` at `2 * i16::MAX` = 65534 (the sum of two
    /// independently-saturated `Short`s) — e.g. measured Phase-7 duplex data: `aD 59034
    /// vs 32767`, `cD 118954 vs 65534`.
    ///
    /// The carve-out is deliberately narrow (a "sound" tolerance, not a blanket ignore):
    /// - Int depth tags `cD`/`cM`/`aD`/`aM`/`bD`/`bM` are each accepted as equal if
    ///   `a == b`, OR (`SAT.contains(a) && b >= a`), OR (`SAT.contains(b) && a >= b`),
    ///   where `SAT = {32767, 65534}` (`i16::MAX` and `2 * i16::MAX`) — i.e. one side hit
    ///   a known saturation point and the other side's true depth is at least as large.
    ///   An UNDERCOUNT relative to a saturated value (e.g. fgbio `32767` vs fgumi `5`) is
    ///   never tolerated — `5 < 32767`, so this still `DIFFER`s.
    /// - The corresponding float error tag (`cE` for `cD`, `aE` for `aD`, `bE` for `bD`)
    ///   is accepted as differing *only* on a record whose own D-tag mismatch was itself
    ///   accepted via the saturation rule above (`cE` gated on `cD`'s outcome, `aE` on
    ///   `aD`'s, `bE` on `bD`'s — independently per strand); otherwise it must match
    ///   exactly.
    /// - Every other tag (including the per-base `cd`/`ce`/`ad`/`ae`/`bd`/`be` `B:s`
    ///   arrays, which saturate identically in both tools) and every core SAM field stay
    ///   exact, so a real depth regression on a normal (non-saturated) family still
    ///   reports a DIFFER.
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
    // trivially agrees, and both `ExactConsensus` and `ExactMinusMi` only ever *add*
    // tolerance on top of `Exact` (the former tolerates saturated-depth divergence, the
    // latter simply drops the MI check — neither ever adds *strictness*), so a byte-for-byte
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

/// The saturation points fgbio's `Short`-typed depth tags clamp to: `i16::MAX` (32767)
/// for the per-strand `aD`/`aM`/`bD`/`bM` (and for simplex/codec `cD`/`cM`), and
/// `2 * i16::MAX` (65534) for *duplex* `cD`/`cM`, which sum two independently-saturated
/// per-strand `Short`s.
const SAT: [i64; 2] = [i16::MAX as i64, 2 * i16::MAX as i64];

/// Returns `true` if two depth-tag values (`cD`/`cM`/`aD`/`aM`/`bD`/`bM`) are equal
/// under the saturation-aware carve-out: either they agree exactly, or one side has hit
/// a known saturation point in [`SAT`] and the other side's true (uncapped) value is at
/// least as large. An undercount relative to a saturated value (e.g. `32767` vs `5`) is
/// deliberately excluded — `5 < 32767`, so neither disjunct fires and this still counts
/// as a mismatch.
fn depth_tag_values_agree(a: i64, b: i64) -> bool {
    a == b || (SAT.contains(&a) && b >= a) || (SAT.contains(&b) && a >= b)
}

/// Returns `true` if `a` and `b` both carry an int-valued `d_tag` whose values disagree
/// but are accepted as equal purely via the saturation branch of
/// [`depth_tag_values_agree`] (i.e. not already exactly equal). This is the trigger
/// condition that additionally permits the D-tag's corresponding error-rate tag (`cE`
/// for `cD`, `aE` for `aD`, `bE` for `bD`) to differ on the same record.
fn d_tag_saturation_accepted(
    tags_a: &BTreeMap<[u8; 2], TagValue<'_>>,
    tags_b: &BTreeMap<[u8; 2], TagValue<'_>>,
    d_tag: SamTag,
) -> bool {
    match (tags_a.get(&*d_tag), tags_b.get(&*d_tag)) {
        (Some(TagValue::Int(va)), Some(TagValue::Int(vb))) => {
            va != vb && depth_tag_values_agree(*va, *vb)
        }
        _ => false,
    }
}

/// Precomputed per-record-pair saturation outcomes for the three depth families this
/// carve-out covers (`cD` combined, `aD`/`bD` per-strand), used to gate whether the
/// corresponding error-rate tag (`cE`/`aE`/`bE`) may differ on that same record. Each
/// flag is independent: a duplex record can have `aD` saturate while `bD` does not (or
/// vice versa), and each strand's error tag must be gated on its own D-tag only.
struct DepthSaturation {
    cd: bool,
    ad: bool,
    bd: bool,
}

impl DepthSaturation {
    /// Computes all three saturation flags once for a record pair's tag maps.
    fn compute(
        tags_a: &BTreeMap<[u8; 2], TagValue<'_>>,
        tags_b: &BTreeMap<[u8; 2], TagValue<'_>>,
    ) -> Self {
        Self {
            cd: d_tag_saturation_accepted(tags_a, tags_b, SamTag::CD),
            ad: d_tag_saturation_accepted(tags_a, tags_b, SamTag::AD),
            bd: d_tag_saturation_accepted(tags_a, tags_b, SamTag::BD),
        }
    }
}

/// Returns `true` if the values of a single named tag are equal under
/// [`ContentPredicate::ExactConsensus`]'s carve-out.
///
/// `saturation` must be precomputed once per record pair by [`DepthSaturation::compute`]
/// and threaded through, since each error-rate tag's carve-out depends on its own D-tag's
/// saturation outcome, not on the error-rate tag's own value.
fn consensus_tag_value_equal(
    tag: [u8; 2],
    va: &TagValue<'_>,
    vb: &TagValue<'_>,
    saturation: &DepthSaturation,
) -> bool {
    if [SamTag::CD, SamTag::CM, SamTag::AD, SamTag::AM, SamTag::BD, SamTag::BM]
        .iter()
        .any(|t| tag == **t)
    {
        match (va, vb) {
            (TagValue::Int(a), TagValue::Int(b)) => depth_tag_values_agree(*a, *b),
            _ => va == vb,
        }
    } else if tag == *SamTag::CE {
        saturation.cd || va == vb
    } else if tag == *SamTag::AE {
        saturation.ad || va == vb
    } else if tag == *SamTag::BE {
        saturation.bd || va == vb
    } else {
        va == vb
    }
}

/// Extracts a record's aux tags as a map of tag name -> zero-copy typed value.
///
/// Values are kept as [`TagValue`] (not stringified): `TagValue` derives `PartialEq` and
/// widens integer tags to `i64` regardless of on-disk width (`c`/`C`/`s`/`S`/`i`/`I`), so
/// comparing `TagValue`s directly already gives the same width-insensitive equality that
/// [`super::super::raw_compare::raw_tags_equal_order_independent`] uses to decide
/// `tag_order_match` — no separate stringified map is needed for [`tag_diffs`], and
/// [`consensus_tag_value_equal`] can pattern-match on the `Int`/`Float` variants for the
/// `cD`/`cM`/`cE`/`aD`/`aM`/`aE`/`bD`/`bM`/`bE` carve-out.
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
/// under [`ContentPredicate::ExactConsensus`] (see [`consensus_tag_value_equal`]).
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
    let saturation = DepthSaturation::compute(&tags_a, &tags_b);
    tags_a.iter().all(|(tag, va)| {
        tags_b.get(tag).is_some_and(|vb| consensus_tag_value_equal(*tag, va, vb, &saturation))
    })
}

/// Builds one diff string per aux tag that is present on only one side, or whose value
/// differs under [`ContentPredicate::ExactConsensus`], between `a` and `b`.
fn consensus_tag_diffs(a: &[u8], b: &[u8]) -> Vec<String> {
    let tags_a = tag_typed_map(a);
    let tags_b = tag_typed_map(b);
    let saturation = DepthSaturation::compute(&tags_a, &tags_b);
    tag_diffs_with(&tags_a, &tags_b, |tag, va, vb| match (va, vb) {
        (Some(a_val), Some(b_val)) => consensus_tag_value_equal(tag, a_val, b_val, &saturation),
        _ => false,
    })
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

    // ==================== ExactConsensus carve-out tests ====================
    //
    // These are mutation-style tests: each one flips a single value (a saturated `cD`,
    // a non-saturated `cD`, or a `cE` difference) and asserts the DIFFER/no-DIFFER
    // outcome the accepted-divergence carve-out promises (see
    // `ContentPredicate::ExactConsensus`'s doc-comment and the compare-hardening
    // design spec's §"Accepted divergences").

    #[test]
    fn exact_consensus_tolerates_saturated_cd_difference() {
        // fgbio saturates cD at i16::MAX (32767); fgumi reports the true depth
        // (59920). Under ExactConsensus this must NOT be a diff.
        let a = SamBuilder::new()
            .read_name(b"r")
            .sequence(b"ACGT")
            .add_int_tag(SamTag::CD, 32767)
            .build();
        let b = SamBuilder::new()
            .read_name(b"r")
            .sequence(b"ACGT")
            .add_int_tag(SamTag::CD, 59920)
            .build();
        assert!(
            content_diffs(a.as_ref(), b.as_ref(), ContentPredicate::ExactConsensus, &hdr(), &hdr())
                .is_none(),
            "a cD difference where one side saturated at 32767 must be tolerated"
        );
    }

    #[test]
    fn exact_consensus_still_flags_non_saturated_cd_difference() {
        // Both values are below the saturation point (32767): this is a real depth
        // discrepancy on a normal-depth family and the carve-out must NOT blind it.
        let a = SamBuilder::new()
            .read_name(b"r")
            .sequence(b"ACGT")
            .add_int_tag(SamTag::CD, 100)
            .build();
        let b = SamBuilder::new()
            .read_name(b"r")
            .sequence(b"ACGT")
            .add_int_tag(SamTag::CD, 200)
            .build();
        let diffs =
            content_diffs(a.as_ref(), b.as_ref(), ContentPredicate::ExactConsensus, &hdr(), &hdr())
                .expect("a non-saturated cD difference must still be reported as a diff");
        assert!(
            diffs.iter().any(|d| d.starts_with("cD:")),
            "diff must name the differing cD tag, got: {diffs:?}"
        );
    }

    #[test]
    fn exact_predicate_still_flags_saturated_cd_difference() {
        // The saturation carve-out is exclusive to ExactConsensus: plain Exact must
        // still DIFFER on the same 32767-vs-59920 pair used by the tolerated case above.
        let a = SamBuilder::new()
            .read_name(b"r")
            .sequence(b"ACGT")
            .add_int_tag(SamTag::CD, 32767)
            .build();
        let b = SamBuilder::new()
            .read_name(b"r")
            .sequence(b"ACGT")
            .add_int_tag(SamTag::CD, 59920)
            .build();
        let diffs = content_diffs(a.as_ref(), b.as_ref(), ContentPredicate::Exact, &hdr(), &hdr())
            .expect("plain Exact must not apply the consensus saturation carve-out");
        assert!(
            diffs.iter().any(|d| d.starts_with("cD:")),
            "diff must name the differing cD tag, got: {diffs:?}"
        );
    }

    #[test]
    fn exact_consensus_tolerates_cd_saturated_cm_and_ce_differences() {
        // cM gets its own independent saturation carve-out, and cE (derived from the
        // saturated depth) is tolerated on a record whose cD accepted-via-saturation.
        let a = SamBuilder::new()
            .read_name(b"r")
            .sequence(b"ACGT")
            .add_int_tag(SamTag::CD, 32767)
            .add_int_tag(SamTag::CM, 32767)
            .add_float_tag(SamTag::CE, 0.01)
            .build();
        let b = SamBuilder::new()
            .read_name(b"r")
            .sequence(b"ACGT")
            .add_int_tag(SamTag::CD, 59920)
            .add_int_tag(SamTag::CM, 40000)
            .add_float_tag(SamTag::CE, 0.02)
            .build();
        assert!(
            content_diffs(a.as_ref(), b.as_ref(), ContentPredicate::ExactConsensus, &hdr(), &hdr())
                .is_none(),
            "cM and cE differences must be tolerated alongside a saturated cD"
        );
    }

    #[test]
    fn exact_consensus_flags_ce_difference_when_cd_not_saturated() {
        // cD matches exactly (no saturation triggered), so cE must still compare
        // exact — the cE carve-out is keyed off cD's saturation outcome, not off cE's
        // own value.
        let a = SamBuilder::new()
            .read_name(b"r")
            .sequence(b"ACGT")
            .add_int_tag(SamTag::CD, 100)
            .add_float_tag(SamTag::CE, 0.01)
            .build();
        let b = SamBuilder::new()
            .read_name(b"r")
            .sequence(b"ACGT")
            .add_int_tag(SamTag::CD, 100)
            .add_float_tag(SamTag::CE, 0.02)
            .build();
        let diffs =
            content_diffs(a.as_ref(), b.as_ref(), ContentPredicate::ExactConsensus, &hdr(), &hdr())
                .expect("a cE difference without a saturated cD must still be reported");
        assert!(
            diffs.iter().any(|d| d.starts_with("cE:")),
            "diff must name the differing cE tag, got: {diffs:?}"
        );
    }

    // ---- Duplex depth-saturation carve-out (aD/aM/bD/bM at 32767, cD/cM at 65534) ----
    //
    // fgbio's *duplex* consensus also saturates the per-strand `aD`/`aM`/`bD`/`bM` at
    // `i16::MAX` (32767, the same single-strand saturation point as simplex `cD`/`cM`),
    // and the *combined* `cD`/`cM` at `2 * i16::MAX` = 65534 (the sum of two saturated
    // `Short`s). Real Phase-7 data: `aD 59034 vs 32767`, `cD 118954 vs 65534`. This
    // table exercises every depth tag at its own saturation point, confirms a
    // non-saturated difference on the same tags still DIFFERs, and proves an
    // UNDERCOUNT (fgbio saturated, fgumi *smaller*) is never tolerated.

    fn depth_tag_record(tag: SamTag, value: i32) -> RawRecord {
        let mut b = SamBuilder::new();
        b.read_name(b"r").sequence(b"ACGT").add_int_tag(tag, value);
        b.build()
    }

    #[rstest]
    // -- saturated at 32767 (per-strand aD/aM/bD/bM, and simplex-style cD/cM) --
    #[case::ad_saturated_32767(SamTag::AD, 32767, 59034, true)]
    #[case::am_saturated_32767(SamTag::AM, 32767, 40000, true)]
    #[case::bd_saturated_32767(SamTag::BD, 32767, 61000, true)]
    #[case::bm_saturated_32767(SamTag::BM, 32767, 41000, true)]
    // -- duplex combined cD/cM saturate at 65534 = 2 * i16::MAX --
    #[case::cd_saturated_65534_duplex(SamTag::CD, 65534, 118_954, true)]
    #[case::cm_saturated_65534_duplex(SamTag::CM, 65534, 90000, true)]
    // -- non-saturated differences on the same tags must still DIFFER --
    #[case::ad_non_saturated_still_differs(SamTag::AD, 100, 200, false)]
    #[case::bd_non_saturated_still_differs(SamTag::BD, 100, 200, false)]
    #[case::cd_non_saturated_still_differs_duplex(SamTag::CD, 40000, 50000, false)]
    // -- an undercount relative to a saturated value is never tolerated: fgbio's
    //    saturated depth is a floor, and 5 < 32767 (or < 65534) is a real regression,
    //    not the "true depth is at least as large" carve-out shape --
    #[case::ad_undercount_not_tolerated(SamTag::AD, 32767, 5, false)]
    #[case::cd_undercount_not_tolerated_duplex(SamTag::CD, 65534, 5, false)]
    fn exact_consensus_depth_tag_saturation_carveout(
        #[case] tag: SamTag,
        #[case] value_a: i32,
        #[case] value_b: i32,
        #[case] expect_tolerated: bool,
    ) {
        let a = depth_tag_record(tag, value_a);
        let b = depth_tag_record(tag, value_b);
        let diffs =
            content_diffs(a.as_ref(), b.as_ref(), ContentPredicate::ExactConsensus, &hdr(), &hdr());
        assert_eq!(
            diffs.is_none(),
            expect_tolerated,
            "{tag:?}: {value_a} vs {value_b} (expect_tolerated={expect_tolerated}): diffs={diffs:?}"
        );
    }

    #[rstest]
    #[case::ad_saturated_32767(SamTag::AD, 32767, 59034)]
    #[case::bd_saturated_32767(SamTag::BD, 32767, 61000)]
    #[case::cd_saturated_65534_duplex(SamTag::CD, 65534, 118_954)]
    fn exact_predicate_still_flags_saturated_duplex_depth_difference(
        #[case] tag: SamTag,
        #[case] value_a: i32,
        #[case] value_b: i32,
    ) {
        // The saturation carve-out is exclusive to ExactConsensus: plain Exact must
        // still DIFFER on the same saturated-vs-true-depth pairs used above, for every
        // depth tag (not just cD).
        let a = depth_tag_record(tag, value_a);
        let b = depth_tag_record(tag, value_b);
        let diffs = content_diffs(a.as_ref(), b.as_ref(), ContentPredicate::Exact, &hdr(), &hdr())
            .expect("plain Exact must not apply the consensus saturation carve-out");
        assert!(
            !diffs.is_empty(),
            "plain Exact must still report a diff for {tag:?}: {value_a} vs {value_b}"
        );
    }

    #[test]
    fn exact_consensus_tolerates_ae_when_ad_saturated() {
        // aE (the per-strand-A error rate) is tolerated only on a record whose aD
        // itself saturated -- mirroring cE's carve-out, but gated on aD specifically.
        let a = SamBuilder::new()
            .read_name(b"r")
            .sequence(b"ACGT")
            .add_int_tag(SamTag::AD, 32767)
            .add_float_tag(SamTag::AE, 0.01)
            .build();
        let b = SamBuilder::new()
            .read_name(b"r")
            .sequence(b"ACGT")
            .add_int_tag(SamTag::AD, 59034)
            .add_float_tag(SamTag::AE, 0.05)
            .build();
        assert!(
            content_diffs(a.as_ref(), b.as_ref(), ContentPredicate::ExactConsensus, &hdr(), &hdr())
                .is_none(),
            "aE must be tolerated alongside a saturated aD"
        );
    }

    #[test]
    fn exact_consensus_tolerates_be_when_bd_saturated() {
        let a = SamBuilder::new()
            .read_name(b"r")
            .sequence(b"ACGT")
            .add_int_tag(SamTag::BD, 32767)
            .add_float_tag(SamTag::BE, 0.01)
            .build();
        let b = SamBuilder::new()
            .read_name(b"r")
            .sequence(b"ACGT")
            .add_int_tag(SamTag::BD, 61000)
            .add_float_tag(SamTag::BE, 0.05)
            .build();
        assert!(
            content_diffs(a.as_ref(), b.as_ref(), ContentPredicate::ExactConsensus, &hdr(), &hdr())
                .is_none(),
            "bE must be tolerated alongside a saturated bD"
        );
    }

    #[test]
    fn exact_consensus_flags_ae_difference_when_ad_not_saturated() {
        // aD matches exactly (no saturation triggered), so aE must still compare
        // exact -- the aE carve-out is keyed off aD's saturation outcome, not aE's
        // own value.
        let a = SamBuilder::new()
            .read_name(b"r")
            .sequence(b"ACGT")
            .add_int_tag(SamTag::AD, 100)
            .add_float_tag(SamTag::AE, 0.01)
            .build();
        let b = SamBuilder::new()
            .read_name(b"r")
            .sequence(b"ACGT")
            .add_int_tag(SamTag::AD, 100)
            .add_float_tag(SamTag::AE, 0.02)
            .build();
        let diffs =
            content_diffs(a.as_ref(), b.as_ref(), ContentPredicate::ExactConsensus, &hdr(), &hdr())
                .expect("an aE difference without a saturated aD must still be reported");
        assert!(
            diffs.iter().any(|d| d.starts_with("aE:")),
            "diff must name the differing aE tag, got: {diffs:?}"
        );
    }

    #[test]
    fn exact_consensus_gates_be_independently_of_ad_saturation() {
        // aD saturates (so aE is tolerated) but bD does NOT saturate on the same
        // record pair: bE must still be reported exactly, proving each strand's
        // error-rate carve-out is gated on *its own* depth tag, not any depth tag.
        let a = SamBuilder::new()
            .read_name(b"r")
            .sequence(b"ACGT")
            .add_int_tag(SamTag::AD, 32767)
            .add_float_tag(SamTag::AE, 0.01)
            .add_int_tag(SamTag::BD, 100)
            .add_float_tag(SamTag::BE, 0.01)
            .build();
        let b = SamBuilder::new()
            .read_name(b"r")
            .sequence(b"ACGT")
            .add_int_tag(SamTag::AD, 59034)
            .add_float_tag(SamTag::AE, 0.05)
            .add_int_tag(SamTag::BD, 100)
            .add_float_tag(SamTag::BE, 0.02)
            .build();
        let diffs =
            content_diffs(a.as_ref(), b.as_ref(), ContentPredicate::ExactConsensus, &hdr(), &hdr())
                .expect("bE must still be reported even though aE is tolerated");
        assert!(
            diffs.iter().any(|d| d.starts_with("bE:")),
            "diff must name the differing bE tag, got: {diffs:?}"
        );
        assert!(
            !diffs.iter().any(|d| d.starts_with("aD:") || d.starts_with("aE:")),
            "aD/aE must not be reported -- aD saturated, so aE is tolerated: {diffs:?}"
        );
    }

    #[test]
    fn exact_consensus_still_exact_on_unrelated_tags() {
        // A tag unrelated to the consensus carve-out (NM) must still be compared
        // exactly under ExactConsensus.
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

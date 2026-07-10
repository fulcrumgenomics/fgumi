//! A collision-resistant identity key for pairing records across two BAM files.
//!
//! Records that share a QNAME — the two ends of a pair, and the primary plus
//! any secondary/supplementary alignments of one end — must remain distinct
//! when the comparator pairs records by identity. The previous key hashed only
//! `(name, is_read1)`, so a primary and a secondary alignment of the same end
//! collapsed to one bucket and a content difference on multi-mapped reads could
//! be masked. `RecordKey` keeps segment, secondary, and supplementary explicit,
//! plus — only for secondary/supplementary alignments — a `(ref_id, pos)`
//! discriminator.
//!
//! Locus is deliberately *not* part of the key for primary alignments: a
//! primary record whose POS legitimately differs between the two BAMs under
//! comparison (e.g. fgumi vs fgbio) must still pair with its counterpart so
//! the POS difference surfaces as a content diff, not as "record present in
//! one file only". Only secondary/supplementary alignments — which can
//! legitimately have several equally-valid loci for the same segment — use
//! locus as a discriminator.
//!
//! This makes `RecordKey` collision-*resistant*, not collision-free: two
//! secondary (or two supplementary) alignments of the same segment at the
//! same locus but with a different CIGAR/strand still collide. That residual
//! is accepted for this phase; a later phase may add a CIGAR/strand
//! discriminator.
//!
//! In the paths this key currently rewires (group/consensus grouping
//! comparisons), every record is a primary — `group` emits primary R1/R2 and
//! consensus reads are unmapped (`pos == -1`) — so `multimap_locus` is always
//! `None` and the key reduces to `(name, segment)`: the same join semantics as
//! the old `(name, is_read1)` key, plus the segment/secondary/supplementary
//! fix. That is why the existing integration tests stay green.

use fgumi_raw_bam::RawRecord;

/// Which read of a template a record belongs to.
///
/// All four FIRST/LAST bit combinations map to a distinct variant so that
/// same-QNAME records with different segment flags never collide into one key —
/// including the malformed `FIRST | LAST` state, which must not be conflated with
/// a plain `Fragment` (neither flag set).
///
/// **The variant declaration order is load-bearing** (it defines the derived
/// `Ord` used by [`RecordKey`]). `FirstAndLast` is deliberately ordered *before*
/// `First`/`Last`, i.e. it does **not** match `fgumi_sort`'s
/// `queryname_flag_order` bit-packing, which sorts the `FIRST | LAST` flag word
/// (`0xc0`) *last*. That intentional disagreement is what keeps the key-join F1
/// soundness guard ([`keyjoin_compare`](super::engines::keyjoin::keyjoin_compare))
/// able to detect — and hard-fail on — a canonicalized stream containing such a
/// malformed record, rather than silently joining it. The
/// `test_keyjoin_compare_hard_errors_on_record_key_order_violation` integration
/// test locks this coupling: reordering `FirstAndLast` to match the flag packing
/// would break it.
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub enum Segment {
    /// Neither FIRST nor LAST segment set (unpaired / fragment read).
    Fragment,
    /// Both `FIRST_OF_PAIR` and `LAST_OF_PAIR` set — a malformed/ambiguous flag
    /// state that is kept distinct from `Fragment` so the two never share an
    /// identity key. Ordered here (before `First`/`Last`) on purpose — see the
    /// enum-level note above.
    FirstAndLast,
    /// `FIRST_OF_PAIR` (flag 0x40) only.
    First,
    /// `LAST_OF_PAIR` (flag 0x80) only.
    Last,
}

/// Collision-resistant identity of a single BAM record for cross-file pairing.
///
/// Primary alignments key on `(name, segment)` alone — see the module docs
/// for why locus must not affect primary identity. Secondary and
/// supplementary alignments additionally key on `multimap_locus` to keep
/// multiple alignments of the same segment distinct.
#[derive(Clone, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub struct RecordKey {
    /// QNAME bytes.
    pub name: Vec<u8>,
    /// R1 / R2 / fragment.
    pub segment: Segment,
    /// SECONDARY alignment (flag 0x100).
    pub secondary: bool,
    /// SUPPLEMENTARY alignment (flag 0x800).
    pub supplementary: bool,
    /// `(ref_id, pos)` of this alignment, but only for secondary/supplementary
    /// records — `None` for primary alignments, which must pair regardless of
    /// POS. Disambiguates multiple secondary/supplementary alignments of the
    /// same segment.
    pub multimap_locus: Option<(i32, i32)>,
}

/// Derive a record's [`Segment`] from its FIRST/LAST-segment flags, without allocating.
///
/// Shared by [`record_key`] and the [`record_keys_match`] hot path so the two never
/// disagree on what counts as `First`/`Last`/`Fragment`.
#[inline]
fn segment_of(raw: &RawRecord) -> Segment {
    match (raw.is_first_segment(), raw.is_last_segment()) {
        (false, false) => Segment::Fragment,
        (true, false) => Segment::First,
        (false, true) => Segment::Last,
        (true, true) => Segment::FirstAndLast,
    }
}

/// Build the collision-resistant identity key for a raw record.
pub fn record_key(raw: &RawRecord) -> RecordKey {
    let segment = segment_of(raw);
    let secondary = raw.is_secondary();
    let supplementary = raw.is_supplementary();
    // Only multi-mapped (secondary/supplementary) alignments use locus as a
    // discriminator; a primary's POS may legitimately differ between the two
    // BAMs under comparison and must not affect whether it pairs with its
    // counterpart.
    let multimap_locus = (secondary || supplementary).then(|| (raw.ref_id(), raw.pos()));
    RecordKey { name: raw.read_name().to_vec(), segment, secondary, supplementary, multimap_locus }
}

/// Returns `true` iff `a` and `b` have the same collision-resistant identity —
/// i.e. `record_key(a) == record_key(b)` — without allocating an owned
/// [`RecordKey`] (in particular, without copying either read name).
///
/// This is the hot-path equivalent of comparing two owned `record_key()` results:
/// it derives segment/secondary/supplementary directly from each record's flag accessors
/// (`is_first_segment`/`is_last_segment`/`is_secondary`/`is_supplementary`), compares those
/// plus the read-name bytes (a borrowed slice comparison), and only consults `(ref_id, pos)`
/// when at least one side is secondary/supplementary — the same short-circuiting a
/// `RecordKey`'s `PartialEq` would perform, minus the allocation.
#[must_use]
pub fn record_keys_match(a: &RawRecord, b: &RawRecord) -> bool {
    if segment_of(a) != segment_of(b) {
        return false;
    }

    let secondary_a = a.is_secondary();
    let secondary_b = b.is_secondary();
    if secondary_a != secondary_b {
        return false;
    }

    let supplementary_a = a.is_supplementary();
    let supplementary_b = b.is_supplementary();
    if supplementary_a != supplementary_b {
        return false;
    }

    if a.read_name() != b.read_name() {
        return false;
    }

    // Only multi-mapped (secondary/supplementary) alignments use locus as a
    // discriminator, matching `record_key`'s `multimap_locus` construction.
    if secondary_a || supplementary_a {
        return a.ref_id() == b.ref_id() && a.pos() == b.pos();
    }

    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use fgumi_raw_bam::SamBuilder;
    use proptest::prelude::*;
    use rstest::rstest;
    // Reuse `fgumi_raw_bam`'s own flag-bit constants (rather than redefining them here) for
    // building raw test records under the old short names these tests already use throughout.
    use fgumi_raw_bam::flags::{
        FIRST_SEGMENT as FLAG_FIRST, LAST_SEGMENT as FLAG_LAST, SECONDARY as FLAG_SECONDARY,
        SUPPLEMENTARY as FLAG_SUPPLEMENTARY,
    };

    /// Build one record with the given name, raw flag word, ref id, and pos.
    fn rec_at(name: &[u8], flags: u16, ref_id: i32, pos: i32) -> RawRecord {
        SamBuilder::new().read_name(name).flags(flags).ref_id(ref_id).pos(pos).build()
    }

    /// Build one unmapped record with the given name and raw flag word.
    fn rec(name: &[u8], flags: u16) -> RawRecord {
        rec_at(name, flags, -1, -1)
    }

    #[test]
    fn distinguishes_r1_from_r2_with_same_name() {
        let r1 = record_key(&rec(b"read", FLAG_FIRST));
        let r2 = record_key(&rec(b"read", FLAG_LAST));
        assert_ne!(r1, r2, "R1 and R2 of the same template must differ");
        assert_eq!(r1.segment, Segment::First);
        assert_eq!(r2.segment, Segment::Last);
    }

    #[test]
    fn distinguishes_primary_secondary_supplementary_same_segment() {
        let primary = record_key(&rec(b"read", FLAG_FIRST));
        let secondary = record_key(&rec(b"read", FLAG_FIRST | FLAG_SECONDARY));
        let supp = record_key(&rec(b"read", FLAG_FIRST | FLAG_SUPPLEMENTARY));
        assert_ne!(primary, secondary, "primary vs secondary must differ");
        assert_ne!(primary, supp, "primary vs supplementary must differ");
        assert_ne!(secondary, supp, "secondary vs supplementary must differ");
    }

    #[test]
    fn distinguishes_two_supplementaries_by_ref_and_pos() {
        // Two supplementary alignments of the same segment differ only by locus.
        let a = record_key(&rec_at(b"read", FLAG_FIRST | FLAG_SUPPLEMENTARY, 1, 100));
        let b = record_key(&rec_at(b"read", FLAG_FIRST | FLAG_SUPPLEMENTARY, 1, 250));
        assert_ne!(a, b, "same-segment supplementaries at different loci must differ");
        assert_eq!(a.multimap_locus, Some((1, 100)));
        assert_eq!(b.multimap_locus, Some((1, 250)));
    }

    #[test]
    fn primaries_pair_regardless_of_locus() {
        // A primary's POS may legitimately differ between the two BAMs under
        // comparison (e.g. fgumi vs fgbio); it must still pair with its
        // counterpart so the POS difference surfaces as a content diff rather
        // than as "record present in one file only".
        let a = record_key(&rec_at(b"read", FLAG_FIRST, 1, 100));
        let b = record_key(&rec_at(b"read", FLAG_FIRST, 2, 999));
        assert_eq!(a, b, "primary records must pair by (name, segment) regardless of locus");
        assert_eq!(a.multimap_locus, None);
    }

    #[test]
    fn equal_for_identical_records() {
        let a = record_key(&rec(b"read", FLAG_FIRST));
        let b = record_key(&rec(b"read", FLAG_FIRST));
        assert_eq!(a, b, "identical records must produce equal keys");
    }

    #[test]
    fn fragment_when_neither_first_nor_last() {
        assert_eq!(record_key(&rec(b"frag", 0)).segment, Segment::Fragment);
    }

    #[test]
    fn first_and_last_set_does_not_collapse_into_fragment() {
        // A record with both FIRST and LAST set is a distinct (malformed but real) flag
        // state; it must not share a key with a plain fragment (neither flag set) of the
        // same name, or the two would collide and pair arbitrarily.
        let both = record_key(&rec(b"read", FLAG_FIRST | FLAG_LAST));
        let fragment = record_key(&rec(b"read", 0));
        assert_eq!(both.segment, Segment::FirstAndLast);
        assert_ne!(
            both, fragment,
            "FIRST|LAST and a plain fragment of the same name must not collide"
        );
    }

    /// `record_keys_match` is a non-allocating equivalent of `record_key(a) == record_key(b)`,
    /// not a different notion of identity — it must agree with the owned comparison on every
    /// case. Each scenario is its own `#[case]` so a failure names the scenario, not an index.
    #[rstest]
    #[case::identical(rec(b"read", FLAG_FIRST), rec(b"read", FLAG_FIRST))]
    #[case::different_names(rec(b"read_a", FLAG_FIRST), rec(b"read_b", FLAG_FIRST))]
    #[case::r1_vs_r2(rec(b"read", FLAG_FIRST), rec(b"read", FLAG_LAST))]
    #[case::primary_vs_secondary(rec(b"read", FLAG_FIRST), rec(b"read", FLAG_FIRST | FLAG_SECONDARY))]
    #[case::primary_vs_supplementary(
        rec(b"read", FLAG_FIRST),
        rec(b"read", FLAG_FIRST | FLAG_SUPPLEMENTARY)
    )]
    #[case::two_secondaries_different_locus(
        rec_at(b"read", FLAG_FIRST | FLAG_SECONDARY, 1, 100),
        rec_at(b"read", FLAG_FIRST | FLAG_SECONDARY, 1, 250)
    )]
    #[case::two_secondaries_same_locus(
        rec_at(b"read", FLAG_FIRST | FLAG_SECONDARY, 1, 100),
        rec_at(b"read", FLAG_FIRST | FLAG_SECONDARY, 1, 100)
    )]
    #[case::primaries_different_locus_still_pair(
        rec_at(b"read", FLAG_FIRST, 1, 100),
        rec_at(b"read", FLAG_FIRST, 2, 999)
    )]
    #[case::fragments(rec(b"frag", 0), rec(b"frag", 0))]
    fn record_keys_match_agrees_with_owned_record_key_equality(
        #[case] a: RawRecord,
        #[case] b: RawRecord,
    ) {
        assert_eq!(
            record_keys_match(&a, &b),
            record_key(&a) == record_key(&b),
            "record_keys_match must agree with record_key(a) == record_key(b) for {a:?} vs {b:?}"
        );
    }

    proptest! {
        /// The fixed case table above can only spot-check the `record_keys_match` ==
        /// `record_key(a) == record_key(b)` invariant; this proves it over arbitrary
        /// records. Names are drawn from a tiny alphabet and full-random flags cover every
        /// FIRST/LAST/SECONDARY/SUPPLEMENTARY combination, so equal names — and thus the
        /// segment/secondary/supplementary/locus discriminators — are exercised frequently,
        /// not just the trivially-unequal-name path. A single divergence here would be a
        /// silent soundness hole in the positional engine's record pairing.
        #[test]
        fn record_keys_match_agrees_with_owned_equality_over_arbitrary_records(
            name_a in "[a-c]{1,3}",
            flags_a in any::<u16>(),
            ref_a in -1i32..4,
            pos_a in -1i32..500,
            name_b in "[a-c]{1,3}",
            flags_b in any::<u16>(),
            ref_b in -1i32..4,
            pos_b in -1i32..500,
        ) {
            let a = rec_at(name_a.as_bytes(), flags_a, ref_a, pos_a);
            let b = rec_at(name_b.as_bytes(), flags_b, ref_b, pos_b);
            prop_assert_eq!(record_keys_match(&a, &b), record_key(&a) == record_key(&b));
        }
    }
}

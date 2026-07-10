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
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub enum Segment {
    /// Neither FIRST nor LAST segment set (unpaired / fragment read).
    Fragment,
    /// `FIRST_OF_PAIR` (flag 0x40).
    First,
    /// `LAST_OF_PAIR` (flag 0x80).
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
        (true, false) => Segment::First,
        (false, true) => Segment::Last,
        _ => Segment::Fragment,
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

    /// `record_keys_match` must agree with `record_key(a) == record_key(b)` on every
    /// case exercised above — it is a non-allocating equivalent of that comparison,
    /// not a different notion of identity.
    #[test]
    fn record_keys_match_agrees_with_owned_record_key_equality() {
        let cases: Vec<(RawRecord, RawRecord)> = vec![
            // Identical records.
            (rec(b"read", FLAG_FIRST), rec(b"read", FLAG_FIRST)),
            // Different read names.
            (rec(b"read_a", FLAG_FIRST), rec(b"read_b", FLAG_FIRST)),
            // Same name, different segment (R1 vs R2).
            (rec(b"read", FLAG_FIRST), rec(b"read", FLAG_LAST)),
            // Same name/segment, primary vs secondary.
            (rec(b"read", FLAG_FIRST), rec(b"read", FLAG_FIRST | FLAG_SECONDARY)),
            // Same name/segment, primary vs supplementary.
            (rec(b"read", FLAG_FIRST), rec(b"read", FLAG_FIRST | FLAG_SUPPLEMENTARY)),
            // Two secondary alignments of the same segment at different loci.
            (
                rec_at(b"read", FLAG_FIRST | FLAG_SECONDARY, 1, 100),
                rec_at(b"read", FLAG_FIRST | FLAG_SECONDARY, 1, 250),
            ),
            // Two secondary alignments of the same segment at the same locus.
            (
                rec_at(b"read", FLAG_FIRST | FLAG_SECONDARY, 1, 100),
                rec_at(b"read", FLAG_FIRST | FLAG_SECONDARY, 1, 100),
            ),
            // Primaries at different loci must still pair (locus ignored for primaries).
            (rec_at(b"read", FLAG_FIRST, 1, 100), rec_at(b"read", FLAG_FIRST, 2, 999)),
            // Fragment reads (neither FIRST nor LAST).
            (rec(b"frag", 0), rec(b"frag", 0)),
        ];

        for (a, b) in &cases {
            let owned_equal = record_key(a) == record_key(b);
            assert_eq!(
                record_keys_match(a, b),
                owned_equal,
                "record_keys_match must agree with record_key(a) == record_key(b) for {a:?} vs {b:?}"
            );
        }
    }
}

//! Per-order run boundary key, used by the arena spill run-former to decide when
//! a freshly-sealed sorted chunk *extends* the currently-open run.
//!
//! # Run extension
//!
//! The owned external sorter coalesces contiguous sorted spill chunks into a
//! single run: when a new chunk's minimum key is `>=` the open run's maximum
//! key, every record in the run precedes every record in the new chunk, so
//! concatenating the new chunk's blocks after the run keeps the run globally
//! sorted. On an already-sorted input, all chunks are contiguous and collapse to
//! one run — a trivial final merge. The arena spill path reproduces this in
//! `SpillGather`; [`RunBound`] is the owned, type-erased key it compares.
//!
//! # Why an owned enum
//!
//! Each sort order carries a different key type ([`RawCoordinateKey`],
//! [`RawQuerynameLexKey`], [`RawQuerynameKey`], and the four narrowed template
//! lanes). The run-former sees chunks type-erased through
//! [`TemplateMemChunk`](crate::TemplateMemChunk) / `MemoryChunkErased`, so it
//! needs one owned type that can hold whichever order's boundary key and compare
//! two bounds of the *same* order. A bound outlives the chunk it came from (the
//! open run's max is retained while later chunks are ingested and the source
//! chunk is spilled and dropped), so it must own its key rather than borrow it.
//!
//! # Queryname tie-break is parity-safe
//!
//! Queryname keys ([`RawQuerynameKey`], [`RawQuerynameLexKey`]) carry an
//! intra-chunk ingest `pos` as their final `cmp` tie-break. Within a chunk `pos`
//! increases with ingest order, so the open run's max (`key_at(len-1)`, a large
//! `pos`) never precedes-or-equals the next chunk's min (`key_at(0)`, `pos`
//! reset near 0) on an exact name+flags tie — extension is conservatively
//! declined. This is safe: `pos` is *not* serialized (the on-disk key is
//! `[name_len][name][flags]`), so an extended run and a non-extended run produce
//! byte-identical spill output. Coordinate and template keys are content-based
//! (no `pos`), so their ties *do* extend — also byte-identical either way.

use crate::inline::{CbKey32, TemplateKey, TemplateKey24, TertKey32};
use crate::keys::{RawCoordinateKey, RawQuerynameKey, RawQuerynameLexKey};

/// The minimum or maximum sort key of an in-memory sorted chunk, owned and
/// type-erased over the sort order so the arena spill run-former can compare the
/// open run's maximum against an incoming chunk's minimum.
///
/// One variant per concrete sort key: the three whole-record orders plus the
/// four narrowed template-coordinate lanes. Two bounds are only ever compared
/// when they come from the same sort run (hence the same variant); see
/// [`precedes_or_equal`](Self::precedes_or_equal).
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum RunBound {
    /// Coordinate order, `K = RawCoordinateKey`.
    Coordinate(RawCoordinateKey),
    /// Queryname order, lexicographic comparator, `K = RawQuerynameLexKey`.
    QuerynameLex(RawQuerynameLexKey),
    /// Queryname order, natural comparator, `K = RawQuerynameKey`.
    QuerynameNatural(RawQuerynameKey),
    /// Template-coordinate, 24-byte core-only lane.
    TemplateK24(TemplateKey24),
    /// Template-coordinate, 32-byte `cb_hash` lane.
    TemplateCb32(CbKey32),
    /// Template-coordinate, 32-byte tertiary lane.
    TemplateTert32(TertKey32),
    /// Template-coordinate, full 40-byte key (all lanes).
    TemplateK40(TemplateKey),
}

impl RunBound {
    /// True iff a chunk whose minimum bound is `next_min` extends the open run
    /// whose maximum bound is `self` — i.e. `self <= next_min`, so the run stays
    /// globally sorted when the incoming chunk's blocks are appended after it.
    ///
    /// Compares within a single sort order. Bounds of differing variants can
    /// only arise from a programming error (a sort run is one order); such a
    /// mismatch conservatively returns `false` (declines to extend) rather than
    /// risk building an out-of-order run.
    #[must_use]
    pub fn precedes_or_equal(&self, next_min: &RunBound) -> bool {
        match (self, next_min) {
            (Self::Coordinate(a), Self::Coordinate(b)) => a <= b,
            (Self::QuerynameLex(a), Self::QuerynameLex(b)) => a <= b,
            (Self::QuerynameNatural(a), Self::QuerynameNatural(b)) => a <= b,
            (Self::TemplateK24(a), Self::TemplateK24(b)) => a <= b,
            (Self::TemplateCb32(a), Self::TemplateCb32(b)) => a <= b,
            (Self::TemplateTert32(a), Self::TemplateTert32(b)) => a <= b,
            (Self::TemplateK40(a), Self::TemplateK40(b)) => a <= b,
            _ => false,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::keys::RawSortKey;

    /// Coordinate keys are content-based (no `pos` tie-break): an open run whose
    /// max is `<=` an incoming chunk's min extends; a `>` boundary does not.
    #[test]
    fn coordinate_extends_when_max_precedes_or_equals_min() {
        let nref = 100;
        let lo = RunBound::Coordinate(RawCoordinateKey::new(1, 1000, false, nref));
        let hi = RunBound::Coordinate(RawCoordinateKey::new(1, 2000, false, nref));

        // open run max (lo) precedes incoming min (hi) → extend.
        assert!(lo.precedes_or_equal(&hi));
        // reverse order → do not extend.
        assert!(!hi.precedes_or_equal(&lo));
    }

    /// An exact boundary tie on a content-based (coordinate) key *does* extend:
    /// equal keys satisfy `<=`, and the spill bytes are identical either way.
    #[test]
    fn coordinate_exact_tie_extends() {
        let k = RawCoordinateKey::new(3, 500, true, 100);
        let a = RunBound::Coordinate(k);
        let b = RunBound::Coordinate(k);
        assert!(a.precedes_or_equal(&b));
    }

    /// Queryname keys carry `pos`, so an exact name+flags boundary tie does NOT
    /// extend: the open run's max (`pos` large) is not `<=` the incoming min
    /// (`pos` reset). Parity-safe because `pos` is not serialized. Uses the lex
    /// comparator; the natural comparator shares the same `pos` tie-break.
    #[test]
    fn queryname_exact_name_tie_does_not_extend_due_to_pos() {
        // Same name + flags, differing only in the intra-chunk ingest position:
        // the open run's max record was ingested later (larger pos) than the
        // incoming chunk's first record (pos reset near 0).
        let mut open_max = RawQuerynameLexKey::new(b"readAAA".to_vec(), 0);
        open_max.set_position(9);
        let mut incoming_min = RawQuerynameLexKey::new(b"readAAA".to_vec(), 0);
        incoming_min.set_position(0);

        let open = RunBound::QuerynameLex(open_max);
        let incoming = RunBound::QuerynameLex(incoming_min);
        assert!(
            !open.precedes_or_equal(&incoming),
            "an exact name tie must NOT extend: pos makes max > min"
        );
    }

    /// A strictly-smaller queryname name in the open run's max still extends
    /// regardless of `pos` — the name comparison decides before `pos`.
    #[test]
    fn queryname_strictly_smaller_name_extends() {
        let mut open_max = RawQuerynameLexKey::new(b"readAAA".to_vec(), 0);
        open_max.set_position(9);
        let mut incoming_min = RawQuerynameLexKey::new(b"readBBB".to_vec(), 0);
        incoming_min.set_position(0);

        let open = RunBound::QuerynameLex(open_max);
        let incoming = RunBound::QuerynameLex(incoming_min);
        assert!(open.precedes_or_equal(&incoming));
    }

    /// Mismatched variants can only arise from a bug; the compare declines to
    /// extend rather than risk an out-of-order run.
    #[test]
    fn variant_mismatch_never_extends() {
        let coord = RunBound::Coordinate(RawCoordinateKey::new(0, 0, false, 100));
        let name = RunBound::QuerynameLex(RawQuerynameLexKey::new(b"r".to_vec(), 0));
        assert!(!coord.precedes_or_equal(&name));
        assert!(!name.precedes_or_equal(&coord));
    }
}

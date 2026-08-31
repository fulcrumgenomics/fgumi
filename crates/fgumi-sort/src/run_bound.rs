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
//! # Tie behaviour is parity-safe (all orders)
//!
//! An exact boundary tie (the open run's max key equals the incoming chunk's min
//! key) *extends* the run in the arena path, for every order. That is safe
//! whether a tie extends or splits, because equal-key records keep the same
//! relative order either way: within one physical run they are concatenated in
//! (chunk seal order, intra-chunk order); across two runs the merge combines the
//! slots in the same (`file_id`, intra-chunk) order — the byte output is identical.
//!
//! Note the `pos` field on the queryname keys ([`RawQuerynameKey`],
//! [`RawQuerynameLexKey`]) does **not** change this. `pos` is a per-chunk ingest
//! index and is the keys' final `cmp` tie-break, but it is set (via
//! [`RawSortKey::set_position`](crate::RawSortKey::set_position)) **only** in the
//! owned engine's `RunFormer`, which never uses `RunBound`. The arena ingest that
//! builds `RunBound` breaks equal-key ties by arena offset in the sort comparator
//! and never calls `set_position`, so every arena queryname key carries `pos == 0`
//! — an exact name+flags tie therefore compares `Equal` and the run extends. And
//! even if `pos` did differ, it is not serialized (the on-disk key is
//! `[name_len][name][flags]`), so extend-vs-split is byte-identical regardless.

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

/// The run-formation summary line both sort engines log, so the wording and the
/// `extended = chunks_spilled - runs` arithmetic have a single source of truth:
/// the owned engine's `log_run_formation` and the arena's
/// `SpillGather::run_formation_summary` both format through this.
///
/// Every spilled chunk either starts a run or extends an existing one, so the
/// number that extended is exactly `chunks_spilled - runs`.
#[must_use]
pub fn format_run_formation(runs: usize, chunks_spilled: usize) -> String {
    format!(
        "Spill runs: {runs} from {chunks_spilled} chunks ({} extended an existing run)",
        chunks_spilled.saturating_sub(runs)
    )
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

    /// The ARENA path never sets `pos` (it breaks ties by arena offset in the
    /// sort comparator, not via `set_position`), so every arena queryname key has
    /// `pos == 0`. An exact name+flags boundary tie therefore compares `Equal` and
    /// the run EXTENDS — this is the real production behaviour, and it is
    /// parity-safe (see the module docs: extend-vs-split is byte-identical).
    #[test]
    fn queryname_exact_name_tie_with_default_pos_extends() {
        // As built by the arena ingest: pos defaults to 0 on both sides.
        let open_max = RawQuerynameLexKey::new(b"readAAA".to_vec(), 0);
        let incoming_min = RawQuerynameLexKey::new(b"readAAA".to_vec(), 0);
        assert_eq!(open_max.name(), incoming_min.name());

        let open = RunBound::QuerynameLex(open_max);
        let incoming = RunBound::QuerynameLex(incoming_min);
        assert!(
            open.precedes_or_equal(&incoming),
            "with the default pos==0 the arena builds, an exact tie extends"
        );
    }

    /// `RunBound` compares the *full* key, so its comparator honours the `pos`
    /// tie-break when the two keys carry different `pos` values. This pins the
    /// comparator itself; note the differing-`pos` state does NOT arise in the
    /// arena path (see `queryname_exact_name_tie_with_default_pos_extends`) — it
    /// is only reachable through the owned engine's `set_position`.
    #[test]
    fn run_bound_queryname_compare_honours_pos_tiebreak() {
        let mut later = RawQuerynameLexKey::new(b"readAAA".to_vec(), 0);
        later.set_position(9);
        let mut earlier = RawQuerynameLexKey::new(b"readAAA".to_vec(), 0);
        earlier.set_position(0);

        // Same name+flags; larger pos does not precede-or-equal the smaller pos.
        assert!(!RunBound::QuerynameLex(later).precedes_or_equal(&RunBound::QuerynameLex(earlier)));
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

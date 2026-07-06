#![allow(dead_code)] // Wired into the sort ingest in increment 3; until then exercised only by tests.

//! Arena layout planner for the parallel-inflate sort ingest.
//!
//! Given the ISIZE (uncompressed size) of each incoming BGZF block, [`ArenaLayout`]
//! computes the gap-aware byte offset at which that block's decompressed bytes will
//! land in a [`crate::segmented_buf::SegmentedBuf`] arena, and decides when to seal
//! the current arena so no arena exceeds the `--max-memory` budget. It is a pure
//! arithmetic planner: it owns no arena and performs no I/O. The offsets it predicts
//! are byte-identical to [`crate::segmented_buf::SegmentedBuf::reserve_contiguous`]
//! (and therefore to `grow_uninit`), so increment 3 can drive a real arena from these
//! placements. Seal points mirror the existing spill trigger
//! (`memory_usage() >= memory_limit`).

/// Where one BGZF block's decompressed bytes land in the arena sequence, plus
/// whether placing it seals the current arena (it reached the memory budget).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct BlockPlacement {
    /// Which arena (monotonically increasing run index, == spill `file_id` order).
    pub(crate) arena_id: u64,
    /// Gap-aware byte offset of the block's first byte within `arena_id`.
    pub(crate) start: usize,
    /// Length of the block's slot in bytes (== the block's ISIZE).
    pub(crate) len: usize,
    /// True iff this block is the last in `arena_id` (the arena reached the budget);
    /// the next `plan` call begins a fresh arena at offset 0.
    pub(crate) seals_arena: bool,
}

/// Plans gap-aware arena offsets and seal points for a stream of BGZF blocks.
///
/// Pure arithmetic — owns no arena. `segment_size` and `budget` are bytes. The
/// offsets reproduce [`SegmentedBuf::reserve_contiguous`](crate::segmented_buf::SegmentedBuf::reserve_contiguous);
/// seal points reproduce the existing `memory_usage() >= memory_limit` spill trigger.
pub(crate) struct ArenaLayout {
    segment_size: usize,
    budget: usize,
    arena_id: u64,
    /// Gap-inclusive bytes placed in the current arena so far (== what
    /// `SegmentedBuf::len()` would report for this arena).
    arena_len: usize,
}

impl ArenaLayout {
    /// Create a planner for the given `segment_size` and byte `budget`.
    ///
    /// # Note
    ///
    /// `budget` should be `> 0`. With `budget == 0`, [`plan`](Self::plan) seals
    /// every block into its own single-block arena (even a zero-length one, since
    /// `0 >= 0`) — safe and bounded, but degenerate. Production always resolves
    /// `--max-memory` to a positive budget; callers that subdivide it (the
    /// `budget / N` arena pool) must keep each share positive.
    pub(crate) fn new(segment_size: usize, budget: usize) -> Self {
        Self { segment_size: segment_size.max(1), budget, arena_id: 0, arena_len: 0 }
    }

    /// Place a block of `isize` uncompressed bytes; return its [`BlockPlacement`].
    ///
    /// # Panics
    ///
    /// Panics if `isize > segment_size` (a slot must fit within one segment).
    pub(crate) fn plan(&mut self, isize: usize) -> BlockPlacement {
        assert!(
            isize <= self.segment_size,
            "block of {} bytes exceeds segment size {}",
            isize,
            self.segment_size,
        );

        // Gap-pad to the next segment boundary, mirroring SegmentedBuf::reserve_contiguous.
        let seg_used = self.arena_len % self.segment_size;
        if seg_used + isize > self.segment_size && seg_used > 0 {
            self.arena_len += self.segment_size - seg_used;
        }

        let start = self.arena_len;
        self.arena_len += isize;

        // Seal AFTER the block that reaches the budget (inclusive >=), matching the
        // existing spill trigger. The exact seal point is parity-neutral.
        let seals_arena = self.arena_len >= self.budget;
        let placement = BlockPlacement { arena_id: self.arena_id, start, len: isize, seals_arena };
        if seals_arena {
            self.arena_id += 1;
            self.arena_len = 0;
        }
        placement
    }

    /// The arena id the next [`plan`](Self::plan) call will place into.
    pub(crate) fn current_arena_id(&self) -> u64 {
        self.arena_id
    }

    /// True when the current arena holds at least one unsealed block — the EOF
    /// "seal the residual run" signal.
    pub(crate) fn has_open_arena(&self) -> bool {
        self.arena_len > 0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::segmented_buf::SegmentedBuf;
    use proptest::prelude::*;

    proptest! {
        // For any sequence of block sizes (incl. zero-length), with a budget that
        // never seals, the planner's offsets equal reserve_contiguous's offsets.
        #[test]
        fn prop_offsets_match_reserve_contiguous(isizes in prop::collection::vec(0usize..=100, 0..50)) {
            let seg = 100usize;
            let mut layout = ArenaLayout::new(seg, usize::MAX);
            let mut oracle = SegmentedBuf::with_capacity(0, seg);
            for &n in &isizes {
                let p = layout.plan(n);
                let want = oracle.reserve_contiguous(n);
                oracle.extend_in_place(&vec![0u8; n]);
                prop_assert_eq!(p.start, want);
                prop_assert_eq!(p.len, n);
                prop_assert!(!p.seals_arena);
                prop_assert_eq!(p.arena_id, 0);
            }
        }

        // For any block sizes and any positive budget: arena_id increments by exactly
        // one and only right after a seal; within an arena offsets are the prefix sum;
        // a sealed arena's cumulative reached the budget while its predecessor did not.
        #[test]
        fn prop_seal_invariants(
            isizes in prop::collection::vec(1usize..=50, 1..40),
            budget in 1usize..=300,
        ) {
            let mut layout = ArenaLayout::new(1_000_000, budget); // no intra-arena gaps for ≤50B blocks
            let mut prev_arena: u64 = 0;
            let mut running: usize = 0;     // cumulative within the current arena
            let mut prev_sealed = false;    // did the previous block seal?
            for &n in &isizes {
                let p = layout.plan(n);
                if prev_sealed {
                    prop_assert_eq!(p.arena_id, prev_arena + 1, "arena_id must increment by 1 after a seal");
                    running = 0;
                } else {
                    prop_assert_eq!(p.arena_id, prev_arena, "arena_id must not change without a seal");
                }
                prop_assert_eq!(p.start, running, "offset must be the within-arena prefix sum");
                running += n;
                prop_assert_eq!(p.seals_arena, running >= budget, "seal iff cumulative reached budget");
                prev_arena = p.arena_id;
                prev_sealed = p.seals_arena;
            }
        }
    }

    #[allow(unsafe_code)]
    #[test]
    fn planned_offsets_drive_grow_uninit_round_trip() {
        // End-to-end bridge to increment 1: the planner predicts exactly where
        // grow_uninit lands each block in a real arena (single arena, MAX budget),
        // and the bytes round-trip through slice_mut/slice.
        let seg = 100usize;
        // Includes a zero-length block (index 2) so the unsafe path exercises
        // grow_uninit(0): its predicted offset must still match grow_uninit's.
        let isizes = [60usize, 60, 0, 30, 40, 50];
        let mut layout = ArenaLayout::new(seg, usize::MAX);
        let mut buf = SegmentedBuf::with_capacity(0, seg);

        for (i, &n) in isizes.iter().enumerate() {
            let p = layout.plan(n);
            // SAFETY: the slot is fully written via slice_mut before any read.
            let off = unsafe { buf.grow_uninit(n) };
            assert_eq!(off, p.start, "grow_uninit landed where the planner predicted");
            unsafe { buf.slice_mut(off, n) }.fill(u8::try_from(i).unwrap());
        }

        let mut replay = ArenaLayout::new(seg, usize::MAX);
        for (i, &n) in isizes.iter().enumerate() {
            let p = replay.plan(n);
            assert_eq!(buf.slice(p.start, n), &vec![u8::try_from(i).unwrap(); n][..]);
        }
    }

    #[test]
    fn plan_offsets_match_reserve_contiguous_across_segment_gaps() {
        // segment_size 100 forces gaps; budget MAX so nothing seals (single arena).
        let seg = 100usize;
        let isizes = [60usize, 60, 30, 40, 50];
        let mut layout = ArenaLayout::new(seg, usize::MAX);
        let mut oracle = SegmentedBuf::with_capacity(0, seg);

        for &n in &isizes {
            let p = layout.plan(n);
            let want = oracle.reserve_contiguous(n);
            oracle.extend_in_place(&vec![0u8; n]); // advance the oracle like a real write
            assert_eq!(
                p.start, want,
                "planner offset diverged from reserve_contiguous for isize={n}"
            );
            assert_eq!(p.len, n);
            assert_eq!(p.arena_id, 0);
            assert!(!p.seals_arena, "MAX budget must never seal");
        }
        // Explicit expected offsets: 0, gap→100, 160, gap→200, 240.
        let mut l2 = ArenaLayout::new(seg, usize::MAX);
        let starts: Vec<usize> = isizes.iter().map(|&n| l2.plan(n).start).collect();
        assert_eq!(starts, vec![0, 100, 160, 200, 240]);
    }

    #[test]
    fn seal_fires_after_block_that_reaches_budget_then_resets() {
        // segment_size huge so no intra-arena gaps; budget 100; 40-byte blocks.
        let mut layout = ArenaLayout::new(1_000_000, 100);
        let p0 = layout.plan(40); // arena_len 40
        let p1 = layout.plan(40); // arena_len 80
        let p2 = layout.plan(40); // arena_len 120 >= 100 → seals
        let p3 = layout.plan(40); // fresh arena

        assert_eq!((p0.arena_id, p0.start, p0.seals_arena), (0, 0, false));
        assert_eq!((p1.arena_id, p1.start, p1.seals_arena), (0, 40, false));
        assert_eq!((p2.arena_id, p2.start, p2.seals_arena), (0, 80, true));
        assert_eq!((p3.arena_id, p3.start, p3.seals_arena), (1, 0, false));
        assert!(layout.has_open_arena(), "arena 1 holds an unsealed 40-byte block");
        assert_eq!(layout.current_arena_id(), 1);
    }

    #[test]
    fn single_block_at_or_over_budget_seals_immediately() {
        let mut layout = ArenaLayout::new(1_000_000, 30);
        let p = layout.plan(50); // 50 >= 30
        assert_eq!((p.arena_id, p.start, p.len, p.seals_arena), (0, 0, 50, true));
        assert!(!layout.has_open_arena(), "the sealing block left no open arena");
    }

    #[test]
    fn has_open_arena_false_before_any_block_and_after_clean_seal() {
        let mut layout = ArenaLayout::new(1_000_000, 100);
        assert!(!layout.has_open_arena(), "no blocks placed yet");
        layout.plan(100); // exactly reaches budget → seals, arena_len resets to 0
        assert!(!layout.has_open_arena(), "exact-budget block sealed cleanly");
    }

    #[test]
    fn zero_length_block_is_placed_without_advancing() {
        let mut layout = ArenaLayout::new(100, usize::MAX);
        let p0 = layout.plan(0);
        let p1 = layout.plan(10);
        assert_eq!((p0.start, p0.len), (0, 0));
        assert_eq!((p1.start, p1.len), (0, 10)); // zero-length block did not move the cursor
    }

    #[test]
    #[should_panic(expected = "exceeds segment size")]
    fn plan_isize_over_segment_size_panics() {
        let mut layout = ArenaLayout::new(100, usize::MAX);
        let _ = layout.plan(101);
    }
}

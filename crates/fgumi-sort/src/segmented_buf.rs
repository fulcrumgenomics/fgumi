#![deny(unsafe_code)]
//! Segmented byte buffer that grows without copying existing data.
//!
//! Unlike `Vec<u8>`, which must reallocate and copy all existing data when it
//! outgrows its capacity, `SegmentedBuf` appends new fixed-size segments.
//! This avoids the O(n) memcpy cost on each doubling and the transient peak
//! memory of holding both old and new buffers simultaneously.
//!
//! Designed as a drop-in replacement for `Vec<u8>` in the sort buffer, where
//! records are appended sequentially and later accessed by `(offset, len)`.

/// Default segment size: 256 MiB.
const DEFAULT_SEGMENT_SIZE: usize = 256 * 1024 * 1024;

/// A growable byte buffer backed by fixed-size segments.
///
/// Each segment is an independent heap allocation, so growing the buffer never
/// copies data already written to *earlier* segments — when the current segment
/// is full, a new one is allocated rather than the whole buffer being doubled.
/// Within the current segment the backing `Vec` may still reallocate and move
/// its own bytes unless it was pre-sized (see
/// [`reserve_full_capacity`](Self::reserve_full_capacity)). Pre-sizing fixes the
/// *inner* buffer only; it does NOT make it sound to run a growth while a
/// `slice_mut` borrow is live, because the growth may also push to the OUTER
/// segment vector. See the `# Safety` note on
/// [`grow_uninit`](Self::grow_uninit) for the supported shape.
///
/// Records written via [`extend_from_slice`](Self::extend_from_slice) are
/// guaranteed not to span segment boundaries: if the current segment lacks
/// room, a new segment is started first.
pub struct SegmentedBuf {
    /// The fixed capacity of each segment.
    segment_size: usize,
    /// Backing storage: one `Vec<u8>` per segment, each with capacity `segment_size`.
    segments: Vec<Vec<u8>>,
    /// Total bytes stored across all segments.
    total_len: usize,
    /// Index of the segment currently being filled (the write cursor). Writes
    /// land in `segments[cur]`; on overflow the cursor advances and **reuses**
    /// `segments[cur+1]` if it already exists (retained by `reset_for_reuse`),
    /// else a new segment is pushed. The overflow path pads `total_len` to the
    /// segment boundary, so `total_len / segment_size == cur` holds at every
    /// segment *transition* — which is what keeps `locate()` exact. It is not a
    /// standing invariant: mid-segment `total_len` is `cur * segment_size +
    /// (bytes written into segment cur)`, so the quotient equals `cur` only
    /// because that remainder is always < `segment_size`. For an actively-growing buffer `cur` is the last
    /// index, so behaviour is identical to the pre-cursor `segments.last()` path.
    cur: usize,
}

#[allow(dead_code)] // Some methods are only exercised from tests for now.
impl SegmentedBuf {
    /// Create a new buffer with the given initial capacity hint and segment size.
    ///
    /// The first segment is pre-allocated immediately.  Additional segments are
    /// allocated on demand.
    #[must_use]
    pub fn with_capacity(capacity: usize, segment_size: usize) -> Self {
        let segment_size = segment_size.max(1);
        // Pre-allocate the first segment up to segment_size.  The capacity hint
        // tells the outer `segments` Vec how many segments to expect, avoiding
        // re-allocations as the buffer grows.
        let first_cap = capacity.min(segment_size);
        let estimated_segments = (capacity / segment_size).max(1);
        let mut segments = Vec::with_capacity(estimated_segments);
        segments.push(Vec::with_capacity(first_cap));
        Self { segment_size, segments, total_len: 0, cur: 0 }
    }

    /// Create a new buffer with a default segment size of 256 MiB.
    #[must_use]
    pub fn new() -> Self {
        Self::with_capacity(0, DEFAULT_SEGMENT_SIZE)
    }

    /// Ensure `needed` bytes fit contiguously in the current segment, starting a
    /// fresh (gap-padded) segment if they don't.
    ///
    /// This is the single source of truth for the "pad `total_len` to the next
    /// segment boundary, then [`advance_segment`](Self::advance_segment)" offset
    /// arithmetic shared by [`extend_from_slice`](Self::extend_from_slice),
    /// [`reserve_contiguous`](Self::reserve_contiguous), and
    /// [`grow_uninit`](Self::grow_uninit). Keeping it in one place is a
    /// correctness requirement, not just DRY: the three call sites must agree
    /// byte-for-byte on where a record lands, and a divergence between them
    /// would silently corrupt sort output identity. Anything that needs to
    /// predict these offsets ahead of the write must drive this buffer rather
    /// than re-derive the arithmetic. Callers read `self.total_len` as the write
    /// offset *after* this returns.
    #[inline]
    fn make_room(&mut self, needed: usize) {
        let seg = &self.segments[self.cur];
        if seg.len() + needed > self.segment_size {
            // Pad total_len to the next segment boundary so locate() works.
            let remainder = self.total_len % self.segment_size;
            if remainder > 0 {
                self.total_len += self.segment_size - remainder;
            }
            self.advance_segment();
        }
    }

    /// Append bytes to the buffer, returning the global offset of the write.
    ///
    /// If the current segment does not have enough remaining capacity for the
    /// entire slice, a new segment is allocated first.  This guarantees the
    /// written bytes are contiguous within a single segment, so callers can
    /// later retrieve them with a single `(offset, len)` reference.
    ///
    /// **Important:** always use the returned offset — do not use [`len`](Self::len)
    /// as a pre-write offset, because gap padding may shift it.
    ///
    /// # Panics
    ///
    /// Panics if `data.len() > segment_size` — a single write must fit in one
    /// segment.
    pub fn extend_from_slice(&mut self, data: &[u8]) -> usize {
        assert!(
            data.len() <= self.segment_size,
            "write of {} bytes exceeds segment size {}",
            data.len(),
            self.segment_size,
        );

        self.make_room(data.len());

        let offset = self.total_len;
        self.segments[self.cur].extend_from_slice(data);
        self.total_len += data.len();
        offset
    }

    /// Advance the write cursor to the next segment, **reusing** a retained
    /// (already-cleared) segment if one exists, else allocating a fresh one.
    /// Reuse is what makes `reset_for_reuse` + refill allocation-free.
    #[inline]
    fn advance_segment(&mut self) {
        self.cur += 1;
        if self.cur == self.segments.len() {
            self.segments.push(Vec::with_capacity(self.segment_size));
        }
        // Always on, release included: `cur` and `total_len` are kept in step by
        // the gap padding in `make_room`, and a non-empty segment here breaks
        // that silently — `locate` then resolves every later offset into the
        // wrong segment and the sort emits wrong bytes with no error. One
        // `is_empty()` per segment transition (one per `segment_size` bytes, so
        // per 256 MiB on the sort path) is not a measurable cost next to that.
        assert!(
            self.segments[self.cur].is_empty(),
            "advance_segment landed on a non-empty segment (reset_for_reuse must clear all)"
        );
        // Uphold the pointer-stability invariant `grow_uninit`'s SAFETY note
        // relies on: the segment the cursor lands in always holds `segment_size`
        // capacity, so a later `grow_uninit` inside it cannot reallocate and move
        // a slot `slice_mut` already handed out. `reserve_full_capacity` primes
        // only `segments[cur]`, and this transition happens inside `make_room`
        // where no caller can observe it and re-prime. A no-op on every path
        // today — segments are pushed `with_capacity(segment_size)` and retained
        // ones were built the same way — but nothing else enforces it.
        let segment_size = self.segment_size;
        let seg = &mut self.segments[self.cur];
        if seg.capacity() < segment_size {
            seg.reserve_exact(segment_size);
        }
    }

    /// Ensure the current segment's backing `Vec` has capacity == `segment_size`.
    ///
    /// After this, any sequence of [`grow_uninit`](Self::grow_uninit) /
    /// [`extend_from_slice`](Self::extend_from_slice) calls that stay within this
    /// segment never reallocate it — so a `&mut [u8]` previously handed out by
    /// [`slice_mut`](Self::slice_mut) for an earlier slot cannot be invalidated by
    /// a later grow. This is what makes the parallel-inflate arena sound: the
    /// serial admit path calls this once per fresh segment before handing any slot
    /// to an inflate worker. No-op if the segment already holds `segment_size`.
    pub fn reserve_full_capacity(&mut self) {
        let seg = &mut self.segments[self.cur];
        let have = seg.capacity();
        if have < self.segment_size {
            seg.reserve_exact(self.segment_size - seg.len());
        }
    }

    /// Ensure at least `additional` bytes fit in the current segment.
    ///
    /// If the remaining capacity in the current segment is less than
    /// `additional`, a new segment is started (with gap padding).
    /// Returns the global offset where the next write will land.
    ///
    /// Use this before a sequence of `extend_from_slice` calls that must
    /// stay contiguous (e.g. header + record body).
    ///
    /// # Panics
    ///
    /// Panics if `additional > segment_size`.
    pub fn reserve_contiguous(&mut self, additional: usize) -> usize {
        assert!(
            additional <= self.segment_size,
            "reserve of {} bytes exceeds segment size {}",
            additional,
            self.segment_size,
        );

        self.make_room(additional);
        self.total_len
    }

    /// Append bytes to the current segment without checking capacity.
    ///
    /// Caller must ensure room via [`reserve_contiguous`](Self::reserve_contiguous)
    /// first. This is the fast path for multi-part writes (header + body).
    ///
    /// # Panics
    ///
    /// Panics (in debug builds) if the write would overflow the current
    /// segment — call [`reserve_contiguous`](Self::reserve_contiguous) first to
    /// start a fresh one. In release the `debug_assert` is compiled out and the
    /// segment simply grows past `segment_size`, which breaks the invariant that
    /// a record never spans a segment boundary.
    #[inline]
    pub fn extend_in_place(&mut self, data: &[u8]) {
        debug_assert!(
            self.segments[self.cur].len() + data.len() <= self.segment_size,
            "extend_in_place exceeds segment capacity; use reserve_contiguous first"
        );
        self.segments[self.cur].extend_from_slice(data);
        self.total_len += data.len();
    }

    /// Reserve a contiguous slot of `additional` **uninitialized-but-live** bytes,
    /// returning the global offset of its first byte. The slot is guaranteed to lie
    /// within a single segment (gap-padding is applied exactly as
    /// [`reserve_contiguous`](Self::reserve_contiguous)), so it can later be written
    /// through [`slice_mut`](Self::slice_mut) and read through [`slice`](Self::slice).
    ///
    /// Unlike `reserve_contiguous` + `extend_in_place`, this grows the segment's
    /// length WITHOUT copying or zero-filling — it is the admit-path primitive for
    /// the parallel-inflate arena, where an inflate worker writes every byte of the
    /// slot exactly once. The returned offsets reproduce `reserve_contiguous`'s
    /// accounting byte-for-byte.
    ///
    /// # Safety
    ///
    /// On return the `additional` bytes at the returned offset are LIVE but
    /// UNINITIALIZED. The caller MUST fully initialize the entire slot (e.g. via
    /// [`slice_mut`](Self::slice_mut)) before reading any of those bytes through
    /// [`slice`](Self::slice) or [`slice_mut`](Self::slice_mut). Reading the slot before it is
    /// fully written is undefined behavior.
    ///
    /// **Pointer stability, and why it is not sufficient for concurrency:**
    /// `grow_uninit` calls `Vec::reserve`, which *can* reallocate and move a
    /// segment's storage — invalidating any `&mut [u8]` previously handed out by
    /// [`slice_mut`](Self::slice_mut) for an earlier slot in the same segment.
    /// [`reserve_full_capacity`](Self::reserve_full_capacity), called once per
    /// fresh segment, fixes that segment at `segment_size` so no subsequent
    /// `grow_uninit` within it reallocates.
    ///
    /// That addresses the *inner* buffer only, and it does NOT make it sound to
    /// run `grow_uninit` while a `slice_mut` borrow is live on another thread.
    /// Two further hazards remain, and neither is fixed by pre-sizing:
    ///
    /// - `grow_uninit` may reach `advance_segment`, which pushes to
    ///   `self.segments`. That push reallocates and frees the OUTER
    ///   `Vec<Vec<u8>>` whenever it is at capacity, leaving a concurrent
    ///   `slice_mut` indexing freed memory — a use-after-free, not merely a
    ///   stale read. It is capacity-dependent, so it fires intermittently.
    /// - `grow_uninit`'s `set_len` writes the inner `Vec`'s length field while a
    ///   concurrent `slice_mut` reads it for its bounds assert — a data race,
    ///   on every call rather than only at a capacity boundary.
    ///
    /// Neither is closed by pre-sizing, and neither can be closed at all through
    /// this API: `grow_uninit` takes `&mut self` and `slice_mut` takes `&self`,
    /// so overlapping them is an aliasing violation whatever the caller does.
    /// That is why the overlap is forbidden outright by `slice_mut`'s third
    /// precondition rather than made safe — Miri reports UB on exactly this
    /// pattern. The supported shape is therefore to reserve ALL slots
    /// for a segment first, and only then hand them to workers — see
    /// `slice_mut_concurrent_disjoint_writes_are_sound`, which does that. In the
    /// single-threaded case the borrow checker already forbids the overlap.
    ///
    /// # Panics
    ///
    /// Panics if `additional > segment_size`.
    #[allow(unsafe_code)]
    #[allow(clippy::uninit_vec)] // intentional: slot is live-but-uninitialized; contract requires caller to fully write it
    pub unsafe fn grow_uninit(&mut self, additional: usize) -> usize {
        assert!(
            additional <= self.segment_size,
            "grow of {} bytes exceeds segment size {}",
            additional,
            self.segment_size,
        );

        self.make_room(additional);

        let offset = self.total_len;
        let seg = &mut self.segments[self.cur];
        let new_len = seg.len() + additional;
        seg.reserve(additional);
        // SAFETY: after `reserve(additional)`, `capacity() >= new_len`, so
        // `set_len(new_len)` only extends `len` over already-allocated bytes. The
        // grown region `old_len..new_len` is intentionally left uninitialized; the
        // documented `# Safety` contract requires the caller to fully write it via
        // `slice_mut` before any read. `u8` has no drop glue and no validity
        // invariant, so growing `len` over uninitialized bytes is not itself UB —
        // only a read-before-write would be, which the contract forbids.
        //
        // NB: this `reserve` MAY reallocate and move the segment's storage; that
        // is safe here only because no `slice_mut` borrow is live across this call
        // — either the borrow checker forbids it (`&mut self`), or the caller
        // pre-sized the segment via `reserve_full_capacity` so this `reserve` is a
        // no-op (see the `# Safety` "Pointer stability" note).
        unsafe {
            seg.set_len(new_len);
        }
        self.total_len += additional;
        offset
    }

    /// Obtain a mutable view of a previously-reserved slot by global `(offset, len)`.
    ///
    /// Takes `&self` (not `&mut self`) so that multiple **disjoint** slots can be
    /// written concurrently from different threads, each holding its own
    /// `&mut [u8]` into the shared buffer — the parallel-inflate use case.
    ///
    /// # Safety
    ///
    /// The caller MUST guarantee that:
    /// 1. `(offset, len)` was produced by [`grow_uninit`](Self::grow_uninit) (or
    ///    `reserve_contiguous` + a matching in-place fill), so it lies within a
    ///    single segment's live region and is in bounds; and
    /// 2. no other `slice`/`slice_mut` borrow of ANY overlapping range is live for
    ///    the duration of the returned reference — callers partition the buffer
    ///    into non-overlapping slots and write each exactly once; and
    /// 3. no `&mut self` method runs for the duration of the returned reference.
    ///    That is every mutating method on this type —
    ///    [`grow_uninit`](Self::grow_uninit),
    ///    [`extend_from_slice`](Self::extend_from_slice),
    ///    [`reserve_contiguous`](Self::reserve_contiguous),
    ///    [`extend_in_place`](Self::extend_in_place),
    ///    [`reserve_full_capacity`](Self::reserve_full_capacity),
    ///    [`reset_for_reuse`](Self::reset_for_reuse), and [`clear`](Self::clear).
    ///
    /// Violating (1) or (2) is undefined behavior (an aliasing `&mut`, or an
    /// out-of-bounds reference). Violating (3) is a data race and can be a
    /// use-after-free: this method reads `self.segments[..]` and the segment's
    /// `len`, both of which a concurrent `grow_uninit` may reallocate or write.
    /// Pre-sizing the segment does not rescue it — see `grow_uninit`'s pointer-
    /// stability note. The supported concurrent shape is to reserve every slot
    /// for a segment *before* handing any of them out, which is what the
    /// parallel-inflate path does.
    ///
    /// # Panics
    ///
    /// Panics if the range spans a segment boundary or exceeds the segment's
    /// live length.
    #[must_use]
    #[allow(unsafe_code)]
    #[allow(clippy::mut_from_ref)] // intentional: &self lets disjoint slots be written concurrently; caller's disjointness contract is the safety invariant
    pub unsafe fn slice_mut(&self, offset: usize, len: usize) -> &mut [u8] {
        // A zero-length slot never needs to be located: `offset` may legitimately
        // sit on a segment boundary (e.g. after `grow_uninit(0)` when the current
        // segment is exactly full), where `locate` would index past the last
        // segment and panic. Return an empty slice instead.
        if len == 0 {
            // Narrow to the case this exists for: `grow_uninit(0)` at an exactly
            // full segment, where `offset` legitimately sits one past the last
            // byte and `locate` would index past the last segment. A zero-length
            // read at a genuinely bogus offset must still fail loudly rather than
            // silently return empty.
            assert!(
                offset <= self.total_len,
                "slice_mut offset {offset} is past the end of the buffer ({})",
                self.total_len,
            );
            return &mut [];
        }
        let (seg_idx, seg_offset) = self.locate(offset);
        let seg = &self.segments[seg_idx];
        assert!(
            seg_offset + len <= seg.len(),
            "slice_mut ({offset}..{}) spans segment boundary (seg {seg_idx}, seg_offset {seg_offset}, seg_len {})",
            offset + len,
            seg.len(),
        );
        // SAFETY: `seg_offset + len <= seg.len()` (asserted) keeps the range inside
        // the segment's live region, so the pointer + length are in bounds. The
        // `&mut [u8]` synthesized from a shared `&self` is sound only under the
        // caller's contract (point 2) that this range is disjoint from every other
        // concurrently-borrowed range, so the produced `&mut` aliases no other
        // `&`/`&mut`. The bytes are live (slot reserved via grow_uninit).
        unsafe {
            let base = seg.as_ptr().add(seg_offset).cast_mut();
            std::slice::from_raw_parts_mut(base, len)
        }
    }

    /// Total bytes stored (including gap padding at segment boundaries).
    #[inline]
    #[must_use]
    pub fn len(&self) -> usize {
        self.total_len
    }

    /// Whether the buffer is empty.
    #[inline]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.total_len == 0
    }

    /// Retrieve a contiguous slice by global `(offset, len)`.
    ///
    /// Because [`extend_from_slice`](Self::extend_from_slice) never splits a
    /// write across segments, any range that was written in a single call is
    /// guaranteed to be within one segment.
    ///
    /// # Panics
    ///
    /// Panics if the range is out of bounds or spans a segment boundary.
    #[inline]
    #[must_use]
    pub fn slice(&self, offset: usize, len: usize) -> &[u8] {
        // A zero-length range never needs to be located: `offset` may legitimately
        // sit on a segment boundary (e.g. after a zero-length reservation when the
        // current segment is exactly full), where `locate` would index past the
        // last segment and panic. Return an empty slice instead.
        if len == 0 {
            // Same narrowing as `slice_mut`: permit the one-past-the-end offset a
            // zero-length reservation produces, reject anything beyond it.
            assert!(
                offset <= self.total_len,
                "slice offset {offset} is past the end of the buffer ({})",
                self.total_len,
            );
            return &[];
        }
        let (seg_idx, seg_offset) = self.locate(offset);
        let seg = &self.segments[seg_idx];
        assert!(
            seg_offset + len <= seg.len(),
            "slice ({offset}..{}) spans segment boundary (seg {seg_idx}, seg_offset {seg_offset}, seg_len {})",
            offset + len,
            seg.len(),
        );
        &seg[seg_offset..seg_offset + len]
    }

    /// Total allocated capacity in bytes across all segments.
    #[must_use]
    pub fn allocated_capacity(&self) -> usize {
        self.segments.iter().map(Vec::capacity).sum()
    }

    /// Number of allocated segments.
    #[must_use]
    pub fn num_segments(&self) -> usize {
        self.segments.len()
    }

    /// Clear all data, retaining only the first segment's allocation (drops
    /// segments `1..`). Used by callers that want memory released between
    /// unrelated sorts (legacy `.sort()`, the streaming path). For
    /// allocation-free arena reuse across spill chunks use
    /// [`reset_for_reuse`](Self::reset_for_reuse) instead.
    pub fn clear(&mut self) {
        self.segments.truncate(1);
        self.segments[0].clear();
        self.total_len = 0;
        self.cur = 0;
    }

    /// Reset for reuse **retaining every segment allocation** (capacity kept).
    /// The next fill reuses the retained segments via the write cursor, so a
    /// pooled arena cycles through fill→spill→reset without reallocating its
    /// segments. Unlike [`clear`](Self::clear), `num_segments()` and
    /// `allocated_capacity()` are unchanged (they hold the peak across reuses).
    pub fn reset_for_reuse(&mut self) {
        for seg in &mut self.segments {
            seg.clear();
        }
        self.total_len = 0;
        self.cur = 0;
    }

    /// Translate a global byte offset to `(segment_index, offset_within_segment)`.
    ///
    /// # Invariant
    ///
    /// All segments except the last are exactly `segment_size` bytes long (full).
    /// A write that doesn't fit in the current segment is always redirected to a
    /// fresh segment with gap-padding applied to `total_len`, so the arithmetic
    /// `offset / segment_size` and `offset % segment_size` are exact.
    #[inline]
    fn locate(&self, offset: usize) -> (usize, usize) {
        let seg_idx = offset / self.segment_size;
        let seg_offset = offset % self.segment_size;
        debug_assert!(
            seg_idx < self.segments.len(),
            "locate({offset}): seg_idx {seg_idx} out of bounds (len {})",
            self.segments.len()
        );
        (seg_idx, seg_offset)
    }

    /// The segment size.
    #[must_use]
    pub fn segment_size(&self) -> usize {
        self.segment_size
    }
}

impl Default for SegmentedBuf {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[test]
    fn test_new_is_empty() {
        let buf = SegmentedBuf::new();
        assert!(buf.is_empty());
        assert_eq!(buf.len(), 0);
        assert_eq!(buf.num_segments(), 1); // one pre-allocated segment
    }

    #[test]
    fn reset_for_reuse_retains_segments_and_reuses_storage() {
        // segment_size 16; ten-byte writes each start a fresh segment (10+10 > 16).
        let mut buf = SegmentedBuf::with_capacity(0, 16);
        let mut offs = Vec::new();
        for i in 0..4u8 {
            offs.push(buf.extend_from_slice(&[i; 10]));
        }
        assert_eq!(buf.num_segments(), 4, "four 10B writes into 16B segments → 4 segments");
        for (i, &o) in offs.iter().enumerate() {
            assert_eq!(buf.slice(o, 10), &[u8::try_from(i).unwrap(); 10][..]);
        }
        let cap_before = buf.allocated_capacity();

        // Reset for reuse: data cleared, but the 4 segment allocations are RETAINED.
        buf.reset_for_reuse();
        assert!(buf.is_empty());
        assert_eq!(buf.len(), 0);
        assert_eq!(buf.num_segments(), 4, "reset retains the segment allocations");
        assert_eq!(buf.allocated_capacity(), cap_before, "no realloc on reset");

        // Refill: must REUSE the retained segments (count stays 4, capacity unchanged),
        // and locate()/slice() must remain correct over the reused segments.
        let mut offs2 = Vec::new();
        for i in 0..4u8 {
            offs2.push(buf.extend_from_slice(&[100 + i; 10]));
        }
        assert_eq!(buf.num_segments(), 4, "refill reuses retained segments (no growth)");
        assert_eq!(buf.allocated_capacity(), cap_before, "refill reused storage (no new alloc)");
        assert_eq!(offs2, offs, "offsets reproduce identically after reuse");
        for (i, &o) in offs2.iter().enumerate() {
            assert_eq!(buf.slice(o, 10), &[100 + u8::try_from(i).unwrap(); 10][..]);
        }
    }

    #[test]
    fn test_extend_and_len() {
        let mut buf = SegmentedBuf::with_capacity(0, 1024);
        let o1 = buf.extend_from_slice(b"hello");
        assert_eq!(o1, 0);
        assert_eq!(buf.len(), 5);
        assert!(!buf.is_empty());

        let o2 = buf.extend_from_slice(b" world");
        assert_eq!(o2, 5);
        assert_eq!(buf.len(), 11);
    }

    #[test]
    fn test_slice_retrieval() {
        let mut buf = SegmentedBuf::with_capacity(0, 1024);
        let o1 = buf.extend_from_slice(b"hello");
        let o2 = buf.extend_from_slice(b" world");

        assert_eq!(buf.slice(o1, 5), b"hello");
        assert_eq!(buf.slice(o2, 6), b" world");
        // Both in same segment, so contiguous range works
        assert_eq!(buf.slice(o1, 11), b"hello world");
    }

    #[test]
    fn test_segment_boundary() {
        // Segment size of 10 bytes
        let mut buf = SegmentedBuf::with_capacity(0, 10);

        // Fill first segment exactly
        let o1 = buf.extend_from_slice(b"0123456789");
        assert_eq!(buf.num_segments(), 1);
        assert_eq!(buf.len(), 10);

        // Next write starts a new segment
        let o2 = buf.extend_from_slice(b"abcde");
        assert_eq!(buf.num_segments(), 2);
        assert_eq!(buf.len(), 15);

        // Data is retrievable across segments
        assert_eq!(buf.slice(o1, 10), b"0123456789");
        assert_eq!(buf.slice(o2, 5), b"abcde");
    }

    #[test]
    fn test_spill_to_new_segment_when_not_enough_room() {
        let mut buf = SegmentedBuf::with_capacity(0, 10);

        let o1 = buf.extend_from_slice(b"1234567"); // 7 bytes in seg 0
        assert_eq!(o1, 0);
        assert_eq!(buf.num_segments(), 1);

        // 5 bytes won't fit in remaining 3 bytes, so spills to seg 1
        // total_len jumps from 7 → 10 (gap) → 15
        let o2 = buf.extend_from_slice(b"abcde");
        assert_eq!(o2, 10); // starts at segment boundary
        assert_eq!(buf.num_segments(), 2);
        assert_eq!(buf.len(), 15); // 10 (seg 0 padded) + 5

        assert_eq!(buf.slice(o1, 7), b"1234567");
        assert_eq!(buf.slice(o2, 5), b"abcde");
    }

    #[test]
    fn test_offset_accounting_with_gaps() {
        // When a write spills to a new segment, total_len must include the
        // gap at the end of the previous segment so that locate() works.
        let mut buf = SegmentedBuf::with_capacity(0, 10);

        let o0 = buf.extend_from_slice(b"aaa"); // seg 0
        assert_eq!(o0, 0);

        let o1 = buf.extend_from_slice(b"bbb"); // seg 0
        assert_eq!(o1, 3);

        // 6 bytes used in seg 0, 5-byte write won't fit → spills to seg 1
        let o2 = buf.extend_from_slice(b"ccccc");
        assert_eq!(o2, 10); // gap of 4 bytes at end of seg 0

        assert_eq!(buf.len(), 15); // 10 (seg 0 padded) + 5

        // All three slices are retrievable
        assert_eq!(buf.slice(o0, 3), b"aaa");
        assert_eq!(buf.slice(o1, 3), b"bbb");
        assert_eq!(buf.slice(o2, 5), b"ccccc");
    }

    #[test]
    #[should_panic(expected = "exceeds segment size")]
    fn test_write_exceeding_segment_panics() {
        let mut buf = SegmentedBuf::with_capacity(0, 10);
        buf.extend_from_slice(b"12345678901"); // 11 bytes > 10
    }

    #[test]
    fn test_clear_resets() {
        let mut buf = SegmentedBuf::with_capacity(0, 1024);
        buf.extend_from_slice(b"some data");
        buf.clear();

        assert!(buf.is_empty());
        assert_eq!(buf.len(), 0);
        assert_eq!(buf.num_segments(), 1);
    }

    /// `clear` truncates to one segment, so it MUST also rewind the write cursor
    /// — otherwise `cur` still points at a segment that no longer exists and the
    /// next write indexes out of bounds. A single-segment buffer cannot catch
    /// this: `cur` is already 0 there, so dropping the reset passes.
    #[test]
    fn clear_rewinds_the_write_cursor_after_multiple_segments() {
        let mut buf = SegmentedBuf::with_capacity(0, 16);
        for _ in 0..4 {
            buf.extend_from_slice(&[7u8; 10]); // 10 > 16/2, so each forces a new segment
        }
        assert!(buf.num_segments() > 1, "precondition: the buffer spans several segments");

        buf.clear();
        assert_eq!(buf.num_segments(), 1, "clear drops the extra segments");

        // The write that would panic if `cur` still pointed past segment 0.
        let off = buf.extend_from_slice(b"after clear");
        assert_eq!(off, 0, "the first post-clear write lands at offset 0");
        assert_eq!(buf.slice(off, 11), b"after clear");
    }

    #[allow(unsafe_code)]
    #[test]
    fn grow_uninit_then_slice_mut_round_trips() {
        // segment_size 100; first slot fits, second forces a gap to a new segment.
        let mut buf = SegmentedBuf::with_capacity(0, 100);

        // SAFETY: each slot is fully written via slice_mut before any slice() read.
        let o0 = unsafe { buf.grow_uninit(60) };
        assert_eq!(o0, 0);
        assert_eq!(buf.len(), 60);
        assert_eq!(buf.num_segments(), 1);
        unsafe { buf.slice_mut(o0, 60) }.fill(0xAB);

        // 60 used, a 60-byte slot won't fit in the remaining 40 → gap to seg 1.
        let o1 = unsafe { buf.grow_uninit(60) };
        assert_eq!(o1, 100, "gap-padded to the segment boundary (matches reserve_contiguous)");
        assert_eq!(buf.num_segments(), 2);
        assert_eq!(buf.len(), 160);
        unsafe { buf.slice_mut(o1, 60) }.fill(0xCD);

        // Both slots read back exactly through the safe slice() path.
        assert_eq!(buf.slice(o0, 60), &[0xAB; 60][..]);
        assert_eq!(buf.slice(o1, 60), &[0xCD; 60][..]);
    }

    #[allow(unsafe_code)]
    #[test]
    fn grow_uninit_offsets_match_reserve_contiguous() {
        // grow_uninit must reproduce reserve_contiguous's offset/gap accounting so
        // the ISIZE prefix-sum planner and the existing reserve path agree.
        let sizes = [3usize, 3, 5, 90, 7];
        let mut a = SegmentedBuf::with_capacity(0, 100);
        let mut b = SegmentedBuf::with_capacity(0, 100);
        for &n in &sizes {
            let ra = a.reserve_contiguous(n);
            a.extend_in_place(&vec![0u8; n]);
            // SAFETY: slot fully written immediately below.
            let rb = unsafe { b.grow_uninit(n) };
            unsafe { b.slice_mut(rb, n) }.fill(0);
            assert_eq!(ra, rb, "grow_uninit offset diverged from reserve_contiguous for n={n}");
        }
        assert_eq!(a.len(), b.len());
        assert_eq!(a.num_segments(), b.num_segments());
    }

    #[allow(unsafe_code)]
    #[test]
    #[should_panic(expected = "exceeds segment size")]
    fn grow_uninit_oversize_panics() {
        let mut buf = SegmentedBuf::with_capacity(0, 10);
        // SAFETY: call panics before any slot is reserved.
        let _ = unsafe { buf.grow_uninit(11) };
    }

    /// Regression: a zero-length reservation whose offset lands exactly on a
    /// segment boundary (the current segment is full) must not make `slice` /
    /// `slice_mut` index past the last segment and panic. Both accessors return
    /// an empty slice for a zero-length range.
    #[allow(unsafe_code)]
    #[test]
    fn zero_length_slice_at_full_segment_boundary_does_not_panic() {
        let mut buf = SegmentedBuf::with_capacity(0, 10);

        // Fill the (only) segment exactly full.
        buf.extend_from_slice(b"0123456789");
        assert_eq!(buf.num_segments(), 1);
        assert_eq!(buf.len(), 10);

        // A zero-length grow returns an offset on the boundary (== total_len),
        // which `locate` would map to a not-yet-allocated segment.
        // SAFETY: the returned slot has length 0, so no byte is ever read or
        // written through it.
        let off = unsafe { buf.grow_uninit(0) };
        assert_eq!(off, 10, "zero-length grow lands on the segment boundary");
        assert_eq!(buf.num_segments(), 1, "zero-length grow allocates nothing");

        // Neither accessor may panic; both yield an empty slice.
        assert!(buf.slice(off, 0).is_empty());
        // SAFETY: zero-length slot — no bytes are aliased, read, or written.
        assert!(unsafe { buf.slice_mut(off, 0) }.is_empty());
    }

    #[test]
    fn test_many_segments() {
        let mut buf = SegmentedBuf::with_capacity(0, 100);

        // Write 50 records of 80 bytes each. First record fits in seg 0 (80 bytes),
        // second doesn't fit (80+80=160 > 100), so each record gets its own segment.
        let mut offsets = Vec::new();
        for i in 0u8..50 {
            let record = vec![i; 80];
            let offset = buf.extend_from_slice(&record);
            offsets.push(offset);
        }

        // Each record in its own 100-byte segment
        assert_eq!(buf.num_segments(), 50);

        // Verify all records readable
        #[allow(clippy::cast_possible_truncation)]
        for (i, &offset) in offsets.iter().enumerate() {
            let data = buf.slice(offset, 80);
            assert_eq!(data[0], i as u8);
        }
    }

    #[test]
    fn test_realistic_record_pattern() {
        // Simulate RecordBuffer usage: header (16 bytes) + record (~200-400 bytes)
        // written as two separate extend_from_slice calls at the same logical offset.
        let mut buf = SegmentedBuf::with_capacity(0, 1024);

        let header = [0u8; 16];
        let record = [42u8; 300];

        let offset = buf.extend_from_slice(&header);
        let record_offset = buf.extend_from_slice(&record);

        // Both in same segment, so header starts at offset, record right after
        assert_eq!(record_offset, offset + 16);
        assert_eq!(buf.slice(offset, 16), &header);
        assert_eq!(buf.slice(record_offset, 300), &record);
    }

    #[test]
    fn test_memory_usage_includes_gaps() {
        // memory_usage should reflect the virtual address space used (including gaps),
        // not just the bytes of actual data written.
        let mut buf = SegmentedBuf::with_capacity(0, 100);

        buf.extend_from_slice(&[0u8; 60]); // seg 0: 60 bytes used
        buf.extend_from_slice(&[0u8; 60]); // spills to seg 1: gap of 40 at end of seg 0

        // len() includes the gap: 100 + 60 = 160
        assert_eq!(buf.len(), 160);
    }

    #[test]
    fn test_reserve_contiguous_then_multi_part_write() {
        // Simulate RecordBuffer: reserve space for header+record, then write parts.
        let mut buf = SegmentedBuf::with_capacity(0, 100);

        // First record: 16-byte header + 60-byte record = 76 bytes
        let offset = buf.reserve_contiguous(76);
        assert_eq!(offset, 0);
        buf.extend_in_place(&[0xAA; 16]); // header
        buf.extend_in_place(&[0xBB; 60]); // record

        // Second record: 16+60 = 76 bytes, only 24 left in seg 0 → spill
        let offset2 = buf.reserve_contiguous(76);
        assert_eq!(offset2, 100); // new segment
        buf.extend_in_place(&[0xCC; 16]);
        buf.extend_in_place(&[0xDD; 60]);

        assert_eq!(buf.num_segments(), 2);
        assert_eq!(buf.slice(0, 16), &[0xAA; 16]);
        assert_eq!(buf.slice(16, 60), &[0xBB; 60]);
        assert_eq!(buf.slice(100, 16), &[0xCC; 16]);
        assert_eq!(buf.slice(116, 60), &[0xDD; 60]);
    }

    #[test]
    fn test_consecutive_writes_same_segment() {
        // Multiple small writes that all fit in one segment should be contiguous.
        let mut buf = SegmentedBuf::with_capacity(0, 1024);

        let o1 = buf.extend_from_slice(b"aaaa");
        let o2 = buf.extend_from_slice(b"bbbb");
        let o3 = buf.extend_from_slice(b"cccc");

        assert_eq!(o1, 0);
        assert_eq!(o2, 4);
        assert_eq!(o3, 8);
        assert_eq!(buf.num_segments(), 1);

        // Can read the entire contiguous range
        assert_eq!(buf.slice(0, 12), b"aaaabbbbcccc");
    }

    #[allow(unsafe_code)]
    #[test]
    fn slice_mut_concurrent_disjoint_writes_are_sound() {
        use std::thread;

        // Reserve N disjoint slots on the serial admit path, then write them all
        // concurrently — each thread owns exactly one slot. This is the
        // parallel-inflate access pattern.
        //
        // What this does NOT prove: `miri.yml` runs only `-p fgumi-raw-bam sort`
        // and `-p fgumi-pipeline-core erased`, so nothing checks this crate under
        // miri. The test exercises the pattern and would catch a slot-arithmetic
        // bug via the content assertions below, but it is not evidence of
        // race-freedom — a benign-looking interleaving that violates `slice_mut`'s
        // disjointness contract would still pass here.
        let n_slots = 8usize;
        let slot_len = 40usize; // 40-byte slots in a 100-byte segment → 2 slots/segment
        let mut buf = SegmentedBuf::with_capacity(0, 100);
        let mut slots = Vec::with_capacity(n_slots);
        for _ in 0..n_slots {
            // SAFETY: every reserved slot is fully written exactly once below before
            // any read; offsets are distinct so the ranges are non-overlapping.
            slots.push(unsafe { buf.grow_uninit(slot_len) });
        }

        let buf_ref = &buf;
        thread::scope(|scope| {
            for (i, &offset) in slots.iter().enumerate() {
                let byte = u8::try_from(i).unwrap();
                scope.spawn(move || {
                    // SAFETY: each thread writes a distinct slot (distinct offset,
                    // fixed len); the ISIZE-prefix-sum analogue here is the distinct
                    // grow_uninit offsets, so the ranges are pairwise disjoint and no
                    // two &mut alias.
                    let dst = unsafe { buf_ref.slice_mut(offset, slot_len) };
                    dst.fill(byte);
                });
            }
        });

        // After the scope joins (write→read happens-before), every slot reads back
        // its writer's byte through the safe slice() path.
        for (i, &offset) in slots.iter().enumerate() {
            let byte = u8::try_from(i).unwrap();
            assert_eq!(buf.slice(offset, slot_len), &vec![byte; slot_len][..]);
        }
    }

    /// `is_err()` alone is satisfied by ANY panic, including one from unrelated
    /// bounds arithmetic, so it cannot tell you the guard under test fired. Both
    /// cases therefore assert on the message. The second case is the one the name
    /// promises: a range that starts in one segment and runs into the next — the
    /// original test only over-read within a single segment, so the true
    /// boundary-spanning path was never exercised.
    #[rstest]
    // One 40-byte slot in a 100-byte segment; asking for 80 runs past the live
    // region of the only segment.
    #[case::past_the_live_region(&[40], 80)]
    // 60 + 60 > 100, so the second slot lands in segment 1 and segment 0 stays
    // 60 bytes live. Asking for 61 from offset 0 crosses out of segment 0 —
    // which is what the name promises and what the old single-slot test could
    // not reach.
    #[case::across_a_segment_boundary(&[60, 60], 61)]
    fn slice_mut_outside_one_segments_live_region_panics(
        #[case] slots: &[usize],
        #[case] ask: usize,
    ) {
        let msg = slice_mut_panic_message(slots, ask);
        assert!(
            msg.contains("spans segment boundary"),
            "must fail slice_mut's own bounds guard, not something incidental; got: {msg}",
        );
    }

    /// Reserve two slots, then ask for `ask` bytes from the first and return the
    /// panic message. Split out because `#[allow(unsafe_code)]` does not reach
    /// the per-case functions `#[rstest]` generates.
    #[allow(unsafe_code)]
    fn slice_mut_panic_message(slots: &[usize], ask: usize) -> String {
        let mut buf = SegmentedBuf::with_capacity(0, 100);
        let mut first = 0;
        for (i, &len) in slots.iter().enumerate() {
            // SAFETY: every reserved slot is fully written immediately below,
            // before anything reads it.
            let off = unsafe { buf.grow_uninit(len) };
            unsafe { buf.slice_mut(off, len) }.fill(u8::try_from(i).expect("few slots"));
            if i == 0 {
                first = off;
            }
        }
        let a = first;

        let err = std::panic::catch_unwind(|| {
            // SAFETY: the call asserts and panics before producing a reference.
            let _ = unsafe { buf.slice_mut(a, ask) };
        })
        .expect_err("slice_mut outside the segment's live region must panic");
        err.downcast_ref::<String>()
            .map_or_else(|| (*err.downcast_ref::<&str>().unwrap_or(&"")).to_string(), Clone::clone)
    }

    #[test]
    fn reserve_full_capacity_makes_segment_zero_full() {
        // A pool-fresh arena has a capacity-0 segment 0; reserve_full_capacity must
        // bring it to segment_size so the first grow_uninit cannot reallocate it.
        let seg = 4096usize;
        let mut buf = SegmentedBuf::with_capacity(0, seg);
        buf.reserve_full_capacity();
        assert!(buf.num_segments() == 1);
        assert!(buf.allocated_capacity() >= seg, "segment 0 must hold >= segment_size capacity");
    }

    /// `reserve_full_capacity` reserves `segment_size - seg.len()`, and the
    /// subtraction matters: `Vec::reserve_exact` takes an amount ADDITIONAL to
    /// the current length, so reserving the full `segment_size` on a partly
    /// filled segment over-allocates by whatever is already written — up to 2x
    /// on a nearly-full segment. Testing only the empty case cannot see that.
    #[test]
    fn reserve_full_capacity_does_not_over_allocate_a_partly_filled_segment() {
        let seg = 4096usize;
        let mut buf = SegmentedBuf::with_capacity(0, seg);
        buf.extend_from_slice(&[0u8; 3000]);
        buf.reserve_full_capacity();
        assert!(
            buf.allocated_capacity() >= seg,
            "must still reach segment_size so later grows cannot realloc",
        );
        assert!(
            buf.allocated_capacity() < seg + 3000,
            "must reserve segment_size - len, not segment_size on top of len; got {}",
            buf.allocated_capacity(),
        );
    }

    #[allow(unsafe_code)]
    #[test]
    fn grow_uninit_is_realloc_free_after_reserve_full_capacity() {
        // After reserve_full_capacity, growing many blocks that sum to <= segment_size
        // must (a) stay in ONE segment (no gap/advance) and (b) never move the segment's
        // backing buffer — so a slice_mut handed out for an early block stays valid.
        let seg = 4096usize;
        let mut buf = SegmentedBuf::with_capacity(0, seg);
        buf.reserve_full_capacity();

        // SAFETY: each slot is fully written before any read; offsets are distinct.
        let off0 = unsafe { buf.grow_uninit(100) };
        let ptr0 = unsafe { buf.slice_mut(off0, 100) }.as_ptr() as usize;
        unsafe { buf.slice_mut(off0, 100) }.fill(0xA1);

        // Grow several more blocks; total stays under segment_size → one segment.
        let mut prev_end = off0 + 100;
        for i in 0..20u8 {
            let n = 100usize;
            let off = unsafe { buf.grow_uninit(n) };
            assert_eq!(off, prev_end, "no gap: offsets are contiguous within one segment");
            unsafe { buf.slice_mut(off, n) }.fill(0xB0 + i);
            prev_end = off + n;
        }
        assert_eq!(buf.num_segments(), 1, "all growth stayed in the pre-sized segment 0");

        // The early block's backing pointer is unchanged → no realloc moved it.
        let ptr0_after = buf.slice(off0, 100).as_ptr() as usize;
        assert_eq!(ptr0, ptr0_after, "segment 0 buffer must not have been reallocated");
        assert_eq!(buf.slice(off0, 100), &[0xA1; 100][..]);
    }

    #[allow(unsafe_code)]
    #[test]
    fn grow_uninit_reuses_arena_after_reset_without_stale_reads() {
        // Increment-5 pattern: a pooled arena is filled, reset_for_reuse'd (len → 0,
        // segment capacity RETAINED), then refilled. The post-reset grow_uninit takes
        // the seg.reserve path on an already-grown segment; confirm offsets reproduce
        // and reads return the fresh bytes, never the retained bytes from the first fill.
        let mut buf = SegmentedBuf::with_capacity(0, 100);

        // First fill: two 40-byte slots (both in segment 0; 40 + 40 <= 100).
        // SAFETY: each slot is fully written via slice_mut before any read.
        let a0 = unsafe { buf.grow_uninit(40) };
        unsafe { buf.slice_mut(a0, 40) }.fill(0x11);
        let a1 = unsafe { buf.grow_uninit(40) };
        unsafe { buf.slice_mut(a1, 40) }.fill(0x22);
        assert_eq!(buf.slice(a0, 40), &[0x11; 40][..]);
        assert_eq!(buf.slice(a1, 40), &[0x22; 40][..]);

        // Reset for reuse: data cleared, segment allocations retained.
        buf.reset_for_reuse();
        assert!(buf.is_empty());

        // Refill: offsets must reproduce the first fill, and reads must return the
        // freshly-written bytes — never the retained 0x11/0x22 still physically present
        // in the reused segment's capacity.
        // SAFETY: each slot is fully written before any read.
        let b0 = unsafe { buf.grow_uninit(40) };
        unsafe { buf.slice_mut(b0, 40) }.fill(0x33);
        let b1 = unsafe { buf.grow_uninit(40) };
        unsafe { buf.slice_mut(b1, 40) }.fill(0x44);
        assert_eq!(b0, a0, "post-reset grow_uninit reproduces the first-fill offset");
        assert_eq!(b1, a1, "post-reset grow_uninit reproduces the first-fill offset");
        assert_eq!(buf.slice(b0, 40), &[0x33; 40][..]);
        assert_eq!(buf.slice(b1, 40), &[0x44; 40][..]);
    }
}

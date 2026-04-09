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
/// Each segment is an independent heap allocation. Appending data never moves
/// existing bytes — when the current segment is full, a new one is allocated.
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
        Self { segment_size, segments, total_len: 0 }
    }

    /// Create a new buffer with a default segment size of 256 MiB.
    #[must_use]
    pub fn new() -> Self {
        Self::with_capacity(0, DEFAULT_SEGMENT_SIZE)
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

        let seg = self.segments.last().expect("segments is never empty");
        if seg.len() + data.len() > self.segment_size {
            // Pad total_len to the next segment boundary so locate() works.
            let remainder = self.total_len % self.segment_size;
            if remainder > 0 {
                self.total_len += self.segment_size - remainder;
            }
            self.segments.push(Vec::with_capacity(self.segment_size));
        }

        let offset = self.total_len;
        self.segments.last_mut().expect("segments is never empty").extend_from_slice(data);
        self.total_len += data.len();
        offset
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

        let seg = self.segments.last().expect("segments is never empty");
        if seg.len() + additional > self.segment_size {
            let remainder = self.total_len % self.segment_size;
            if remainder > 0 {
                self.total_len += self.segment_size - remainder;
            }
            self.segments.push(Vec::with_capacity(self.segment_size));
        }
        self.total_len
    }

    /// Append bytes to the current segment without checking capacity.
    ///
    /// Caller must ensure room via [`reserve_contiguous`](Self::reserve_contiguous)
    /// first. This is the fast path for multi-part writes (header + body).
    ///
    /// # Panics
    ///
    /// Panics (in debug builds) if the segments vec is empty.
    #[inline]
    pub fn extend_in_place(&mut self, data: &[u8]) {
        debug_assert!(
            self.segments.last().is_some_and(|s| s.len() + data.len() <= self.segment_size),
            "extend_in_place exceeds segment capacity; use reserve_contiguous first"
        );
        self.segments.last_mut().expect("segments is never empty").extend_from_slice(data);
        self.total_len += data.len();
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

    /// Clear all data, retaining the first segment's allocation.
    pub fn clear(&mut self) {
        self.segments.truncate(1);
        self.segments[0].clear();
        self.total_len = 0;
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
    use super::*;

    #[test]
    fn test_new_is_empty() {
        let buf = SegmentedBuf::new();
        assert!(buf.is_empty());
        assert_eq!(buf.len(), 0);
        assert_eq!(buf.num_segments(), 1); // one pre-allocated segment
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
}

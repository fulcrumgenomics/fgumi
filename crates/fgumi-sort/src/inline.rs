//! Inline record buffer for samtools-style memory layout.
//!
//! This module provides a contiguous buffer that stores BAM records inline
//! with pre-computed sort keys, eliminating per-record heap allocations.
//!
//! Key benefits over `Vec<(Key, Record)>`:
//! - Single large allocation instead of millions of small ones
//! - Better cache locality - sequential memory access
//! - Sort by index array, not by moving actual record data
//! - ~24 bytes overhead per record vs ~110 bytes
//!
//! Some helpers (alternative buffer types, reference-array variants) are
//! retained for benchmarking and are reachable only inside the crate.

#![allow(dead_code)]

use crate::keys::{RawCoordinateKey, RawSortKey, SortContext};
use crate::radix::bytes_needed_u64;
use crate::segmented_buf::SegmentedBuf;
use fgumi_raw_bam::RawRecordView;
use std::cmp::Ordering;
use std::io::{Read, Write};
use std::sync::Arc;

// ============================================================================
// In-memory sorted chunk (shared-buffer)
// ============================================================================

/// A sorted in-memory chunk produced by `par_sort_into_chunks`.
///
/// All chunks produced by a single sort share the original
/// `RecordBuffer`'s `SegmentedBuf` via `Arc`, instead of each record
/// owning its own `Vec<u8>`. The previous design called
/// `RawRecord::from(buffer.get_record(r).to_vec())` once per record,
/// paying one `mi_malloc` + memcpy per record (~16 M of each for the
/// twist-umi sort). The shared-buffer design pays zero per-record
/// allocations and reads record bytes directly out of the shared
/// segments at merge time.
///
/// The merge consumer still copies record bytes into its own scratch
/// `Vec<u8>` via `extend_from_slice` — the original `mem::swap`
/// zero-copy bridge no longer applies because the bytes are now
/// borrowed from the shared `SegmentedBuf`, not owned per record.
/// In exchange, materialization is allocation-free and the peak
/// memory of the sort drops by ~1× (the materialization no longer
/// transiently doubles the buffer).
pub(crate) struct InMemoryChunk<K> {
    /// Shared backing store for record bytes (the original sort
    /// buffer's `SegmentedBuf`). All sibling chunks from one
    /// `par_sort_into_chunks` call share this Arc, so the segments
    /// are freed only when the last chunk drops. The buffer stays
    /// resident for the entire merge phase, unlike the pre-PR
    /// behavior where per-record `Vec<u8>`s could free incrementally
    /// as the merge consumed them — steady-state RSS during the
    /// merge is slightly higher, but the materialization-peak
    /// memory drops by ~1× (pre-PR transiently held both the buffer
    /// and a full copy in per-record `Vec<u8>`s).
    data: Arc<SegmentedBuf>,
    /// `(sort_key, byte_offset_in_data, len)` per record, in sorted
    /// order. `offset` is `u64` because a single `SegmentedBuf` can
    /// exceed 4 GiB at high `--max-memory` × `--threads`; `len` is
    /// `u32` because individual BAM records are < 4 GiB.
    records: Vec<(K, u64, u32)>,
}

impl<K> InMemoryChunk<K> {
    /// Construct an empty chunk backed by an empty shared buffer.
    #[must_use]
    pub(crate) fn empty() -> Self {
        Self { data: Arc::new(SegmentedBuf::default()), records: Vec::new() }
    }

    /// Construct a chunk holding the given records, all referencing
    /// the same shared data buffer.
    #[must_use]
    pub(crate) fn from_parts(data: Arc<SegmentedBuf>, records: Vec<(K, u64, u32)>) -> Self {
        Self { data, records }
    }

    /// Number of records in the chunk.
    #[must_use]
    pub(crate) fn len(&self) -> usize {
        self.records.len()
    }

    /// Whether the chunk is empty.
    #[must_use]
    pub(crate) fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    /// Borrow the `i`th record's bytes from the shared data buffer.
    #[must_use]
    #[allow(clippy::cast_possible_truncation)] // offset/len fit in usize on all supported targets
    pub(crate) fn record_bytes(&self, i: usize) -> &[u8] {
        let (_, offset, len) = &self.records[i];
        self.data.slice(*offset as usize, *len as usize)
    }

    /// Borrow the `i`th record's sort key.
    #[must_use]
    pub(crate) fn key_at(&self, i: usize) -> &K {
        &self.records[i].0
    }

    /// Move out the `i`th record's key (replaces it with
    /// `K::default()`). Used by the merge loop after the
    /// loser-tree consumes the key.
    pub(crate) fn take_key(&mut self, i: usize) -> K
    where
        K: Default,
    {
        std::mem::take(&mut self.records[i].0)
    }
}

impl<K> Default for InMemoryChunk<K> {
    fn default() -> Self {
        Self::empty()
    }
}

// ============================================================================
// Buffer Probe Trait
// ============================================================================

/// Common metrics shared by `RecordBuffer` and `TemplateRecordBuffer` for
/// memory-probe instrumentation.
pub trait ProbeableBuffer {
    /// Logical bytes stored (data + refs).
    fn memory_usage(&self) -> usize;
    /// Total allocated capacity (segments + refs Vec).
    fn allocated_capacity(&self) -> usize;
    /// Number of records in the buffer.
    fn len(&self) -> usize;
    /// Whether the buffer is empty.
    fn is_empty(&self) -> bool {
        self.len() == 0
    }
    /// Number of data segments.
    fn num_segments(&self) -> usize;
}

// ============================================================================
// Packed Sort Keys
// ============================================================================

/// Packed sort key for coordinate ordering.
///
/// Format: `(tid << 34) | ((pos+1) << 1) | reverse`
///
/// This allows single u64 comparison for most records. The +1 on pos
/// ensures pos=0 doesn't collide with unmapped (which uses MAX).
#[repr(transparent)]
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Debug)]
pub struct PackedCoordinateKey(pub u64);

impl PackedCoordinateKey {
    /// Create a packed coordinate key.
    ///
    /// # Arguments
    /// * `tid` - Reference sequence ID (-1 for unmapped)
    /// * `pos` - 0-based alignment position
    /// * `reverse` - True if reverse complemented
    /// * `nref` - Number of reference sequences (for unmapped handling)
    #[inline]
    #[must_use]
    #[allow(clippy::cast_sign_loss)]
    pub fn new(tid: i32, pos: i32, reverse: bool, _nref: u32) -> Self {
        if tid < 0 {
            return Self::unmapped();
        }
        let tid = tid as u32;
        // Pack: tid in high bits, (pos+1) in middle, reverse in LSB
        // Using pos+1 so that pos=0 doesn't become 0 in the key
        #[allow(clippy::cast_lossless)] // Explicit bit packing requires precise control
        let key = (u64::from(tid) << 34)
            | ((i64::from(pos) as u64).wrapping_add(1) << 1)
            | u64::from(reverse);
        Self(key)
    }

    /// Create a key for unmapped records (sorts after all mapped).
    #[inline]
    #[must_use]
    pub fn unmapped() -> Self {
        Self(u64::MAX)
    }
}

// ============================================================================
// Record References (Index for Sorting)
// ============================================================================

/// Reference to a record in the buffer (used for sorting).
///
/// This is a lightweight handle that can be sorted efficiently.
/// The actual record data stays in place in the buffer.
///
/// Note: No read name tie-breaking is used, matching samtools behavior.
/// Equal records maintain their original input order (stable sort).
#[repr(C)]
#[derive(Copy, Clone, Debug)]
pub struct RecordRef {
    /// Packed primary sort key for fast comparison.
    pub sort_key: u64,
    /// Offset into `RecordBuffer` where record header starts.
    pub offset: u64,
    /// Length of raw BAM data (excluding inline header).
    pub len: u32,
    /// Padding for 8-byte alignment.
    padding: u32,
}

impl RecordRef {
    /// Construct a `RecordRef` from its public fields (padding is internal).
    /// Primarily for tests and benchmarks that build sort indices directly.
    #[must_use]
    pub fn new(sort_key: u64, offset: u64, len: u32) -> Self {
        Self { sort_key, offset, len, padding: 0 }
    }
}

impl PartialEq for RecordRef {
    fn eq(&self, other: &Self) -> bool {
        self.sort_key == other.sort_key
    }
}

impl Eq for RecordRef {}

impl PartialOrd for RecordRef {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for RecordRef {
    fn cmp(&self, other: &Self) -> Ordering {
        self.sort_key.cmp(&other.sort_key)
    }
}

// ============================================================================
// Inline Header (stored before each record in buffer)
// ============================================================================

/// Inline record header stored in buffer before raw BAM data.
///
/// This header is written once when the record is added and records the
/// length of the BAM bytes that follow it. The sort key is *not* stored
/// here: it is cached in the parallel `RecordRef` index, which is the only
/// place the sort reads it from. Record bytes are accessed via
/// `offset + HEADER_SIZE`, so the header itself is never read back during
/// sorting or merge — only its 8-byte footprint matters.
#[repr(C)]
#[derive(Copy, Clone, Debug)]
struct InlineHeader {
    /// Length of following raw BAM data.
    record_len: u32,
    /// Padding for 8-byte alignment.
    padding: u32,
}

const HEADER_SIZE: usize = std::mem::size_of::<InlineHeader>(); // 8 bytes
const _: () = assert!(
    std::mem::size_of::<InlineHeader>() == HEADER_SIZE,
    "HEADER_SIZE must match size_of::<InlineHeader>()"
);

// ============================================================================
// RecordBuffer - Main Data Structure
// ============================================================================

/// Contiguous buffer for inline record storage.
///
/// Records are stored sequentially in a single allocation:
/// ```text
/// [Header0][BAM0][Header1][BAM1][Header2][BAM2]...
/// ```
///
/// An index array of `RecordRef` is maintained for sorting.
/// Sorting only reorders the index; records stay in place.
pub struct RecordBuffer {
    /// Segmented byte storage for all records (headers + BAM data).
    data: SegmentedBuf,
    /// Index of record references for sorting.
    refs: Vec<RecordRef>,
    /// Number of reference sequences (for unmapped handling).
    nref: u32,
}

/// Segment size for in-memory sort buffers: 256 MiB.
///
/// Both `RecordBuffer` and `TemplateRecordBuffer` use this segment size so that
/// a single BAM record (≤128 MiB in practice) always fits within one segment.
const SORT_SEGMENT_SIZE: usize = 256 * 1024 * 1024;

/// Shared implementation for `par_sort_into_chunks` on both buffer types.
///
/// Sorts refs in place (parallel radix sort, partitioned by `chunk_size`),
/// then drains `self.data` into a single `Arc<SegmentedBuf>` shared
/// across all produced chunks. Each chunk's `records` Vec is built from
/// its refs' `(key, offset + header_size, len)` triples — no per-record
/// allocation or memcpy.
///
/// `$header_size` is the byte distance from each ref's `offset` field
/// to the start of its actual record bytes inside `data`
/// (`HEADER_SIZE` for `RecordBuffer`, `TEMPLATE_HEADER_SIZE` for
/// `TemplateRecordBuffer`).
macro_rules! par_sort_into_chunks_impl {
    ($self:expr, $threads:expr, $sort_fn:expr, $header_size:expr, $key_fn:expr) => {{
        use rayon::prelude::*;
        use std::sync::Arc;

        let n = $self.refs.len();
        let header_size: u64 = $header_size as u64;

        if $threads <= 1 || n < RADIX_THRESHOLD * 2 || n <= 10_000 {
            $sort_fn(&mut $self.refs);
            let data = Arc::new(std::mem::take(&mut $self.data));
            let records: Vec<_> =
                $self.refs.iter().map(|r| ($key_fn(r), r.offset + header_size, r.len)).collect();
            // Clear refs to leave the buffer in a consistent drained
            // state (both `data` and `refs` empty). Callers must drop
            // the buffer after this — `iter_sorted()` would yield zero
            // records and `get_record()` with a stale `RecordRef` would
            // index-panic against the empty `SegmentedBuf`.
            $self.refs.clear();
            return vec![InMemoryChunk::from_parts(data, records)];
        }

        let chunk_size = n.div_ceil($threads);

        $self.refs.par_chunks_mut(chunk_size).for_each(|chunk| {
            $sort_fn(chunk);
        });

        let data = Arc::new(std::mem::take(&mut $self.data));

        let chunks: Vec<_> = $self
            .refs
            .chunks(chunk_size)
            .map(|chunk_refs| {
                let records: Vec<_> = chunk_refs
                    .iter()
                    .map(|r| ($key_fn(r), r.offset + header_size, r.len))
                    .collect();
                InMemoryChunk::from_parts(Arc::clone(&data), records)
            })
            .collect();
        $self.refs.clear();
        chunks
    }};
}

impl RecordBuffer {
    /// Create a new buffer with estimated capacity.
    ///
    /// # Arguments
    /// * `estimated_records` - Expected number of records
    /// * `estimated_bytes` - Expected total bytes of BAM data
    /// * `nref` - Number of reference sequences in header
    #[must_use]
    pub fn with_capacity(estimated_records: usize, estimated_bytes: usize, nref: u32) -> Self {
        Self {
            data: SegmentedBuf::with_capacity(
                estimated_bytes + estimated_records * HEADER_SIZE,
                SORT_SEGMENT_SIZE,
            ),
            refs: Vec::with_capacity(estimated_records),
            nref,
        }
    }

    /// Push a record for coordinate sorting.
    ///
    /// Extracts the sort key inline from raw BAM bytes (zero-copy).
    ///
    /// # Errors
    ///
    /// Returns an error if the record (plus header) exceeds the segment size
    /// (256 MiB) or if the record length exceeds `u32::MAX`.
    #[inline]
    pub fn push_coordinate(&mut self, record: &[u8]) -> anyhow::Result<()> {
        // BAM fixed-length block is 32 bytes; coordinate fields end at offset 15.
        const MIN_BAM_RECORD_LEN: usize = 16;
        anyhow::ensure!(
            record.len() >= MIN_BAM_RECORD_LEN,
            "BAM record is truncated: need at least {} bytes to extract coordinate fields, got {}",
            MIN_BAM_RECORD_LEN,
            record.len(),
        );

        let total_bytes = HEADER_SIZE + record.len();
        anyhow::ensure!(
            total_bytes <= SORT_SEGMENT_SIZE,
            "BAM record of {} bytes (+ {} byte header) exceeds segment size of {} bytes; \
             this is likely a malformed BAM file",
            record.len(),
            HEADER_SIZE,
            SORT_SEGMENT_SIZE,
        );
        let len = u32::try_from(record.len())
            .map_err(|_| anyhow::anyhow!("record length {} exceeds u32::MAX", record.len()))?;

        // Extract sort key from raw BAM bytes. The key is cached only in the
        // `RecordRef` index below (the live reader); it is intentionally not
        // stored in the inline header, which now records only the length.
        let sort_key = extract_coordinate_key_inline(record, self.nref);

        // Reserve contiguous space for header + record
        let offset = self.data.reserve_contiguous(total_bytes) as u64;

        // Write inline header (8 bytes: record_len + padding)
        let header = InlineHeader { record_len: len, padding: 0 };
        self.data.extend_in_place(&header.record_len.to_le_bytes());
        self.data.extend_in_place(&header.padding.to_le_bytes());

        // Write raw BAM data
        self.data.extend_in_place(record);

        // Add to index
        self.refs.push(RecordRef { sort_key, offset, len, padding: 0 });
        Ok(())
    }

    /// Sort the index by key (records stay in place).
    ///
    /// Uses radix sort for O(n×k) performance instead of O(n log n) comparison sort.
    /// Falls back to insertion sort for small arrays.
    pub fn sort(&mut self) {
        radix_sort_record_refs(&mut self.refs);
    }

    /// Sort using parallel radix sort (for large arrays).
    ///
    /// Divides data into chunks, sorts each with radix sort, then merges.
    /// For very large arrays this is faster than single-threaded radix sort.
    pub fn par_sort(&mut self) {
        // For parallel sort, we use chunked radix sort
        // Radix sort is already O(n×k), so parallelizing chunks provides
        // linear speedup without the merge overhead of comparison-based parallel sorts
        parallel_radix_sort_record_refs(&mut self.refs);
    }

    /// Sort in parallel and return each sub-array as a separate sorted chunk.
    ///
    /// Instead of merging the parallel sort sub-arrays back into one sorted
    /// buffer (as `par_sort` does), this returns each sub-array as its own
    /// `InMemoryChunk<RawCoordinateKey>` so they can be passed as separate
    /// merge sources to the k-way merge, avoiding the intermediate merge step.
    ///
    /// When `threads <= 1`, returns a single chunk.
    ///
    /// After this returns, `self.data` and `self.refs` are both empty —
    /// the produced chunks share the original data via `Arc<SegmentedBuf>`.
    /// The caller must drop the buffer. `iter_sorted()` would yield zero
    /// records and `get_record()` called with a stale `RecordRef` would
    /// index-panic against the empty `SegmentedBuf`.
    pub(crate) fn par_sort_into_chunks(
        &mut self,
        threads: usize,
    ) -> Vec<InMemoryChunk<RawCoordinateKey>> {
        // Each chunk derives its own bound inside its byte-0 counting pass, so
        // there is no whole-buffer scan to share and no bound to pass around. A
        // per-chunk bound is also never wider than a shared one, so chunks whose
        // keys happen to be narrow run fewer passes.
        let sort_chunk = |refs: &mut [RecordRef]| radix_sort_record_refs(refs);
        par_sort_into_chunks_impl!(self, threads, sort_chunk, HEADER_SIZE, |r: &RecordRef| {
            RawCoordinateKey { sort_key: r.sort_key }
        })
    }

    /// Drain the (already-sorted) buffer into a single in-memory chunk
    /// whose records share an `Arc<SegmentedBuf>` backing store.
    ///
    /// Used by sort entry points that call `par_sort` directly (e.g.
    /// the single-thread non-rayon path) and then need to hand the
    /// sorted records to the merge as one chunk. After this returns,
    /// both `self.data` and `self.refs` are empty; the caller must
    /// drop the buffer. `iter_sorted()` would yield zero records and
    /// `get_record()` called with a stale `RecordRef` would
    /// index-panic against the empty `SegmentedBuf`.
    #[must_use]
    pub(crate) fn drain_into_single_chunk(&mut self) -> InMemoryChunk<RawCoordinateKey> {
        let data = Arc::new(std::mem::take(&mut self.data));
        let records: Vec<_> = self
            .refs
            .iter()
            .map(|r| {
                let key = RawCoordinateKey { sort_key: r.sort_key };
                (key, r.offset + HEADER_SIZE as u64, r.len)
            })
            .collect();
        self.refs.clear();
        InMemoryChunk::from_parts(data, records)
    }

    /// Get record bytes by reference.
    #[inline]
    #[must_use]
    #[allow(clippy::cast_possible_truncation)]
    pub fn get_record(&self, r: &RecordRef) -> &[u8] {
        self.data.slice(r.offset as usize + HEADER_SIZE, r.len as usize)
    }

    /// Iterate over sorted records.
    pub fn iter_sorted(&self) -> impl Iterator<Item = &[u8]> {
        self.refs.iter().map(|r| self.get_record(r))
    }

    /// Get the record references.
    #[must_use]
    pub fn refs(&self) -> &[RecordRef] {
        &self.refs
    }

    /// Memory usage in bytes (actual data stored, not capacity).
    #[must_use]
    pub fn memory_usage(&self) -> usize {
        self.data.len() + self.refs.len() * std::mem::size_of::<RecordRef>()
    }

    /// Total allocated capacity in bytes (data segments + refs Vec).
    #[must_use]
    pub fn allocated_capacity(&self) -> usize {
        self.data.allocated_capacity() + self.refs.capacity() * std::mem::size_of::<RecordRef>()
    }

    /// Number of data segments in the underlying buffer.
    #[must_use]
    pub fn num_segments(&self) -> usize {
        self.data.num_segments()
    }

    /// Number of records.
    #[must_use]
    pub fn len(&self) -> usize {
        self.refs.len()
    }

    /// Check if buffer is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.refs.is_empty()
    }

    /// Clear the buffer for reuse.
    pub fn clear(&mut self) {
        self.data.clear();
        self.refs.clear();
    }

    /// Get the number of reference sequences.
    #[must_use]
    pub fn nref(&self) -> u32 {
        self.nref
    }
}

impl ProbeableBuffer for RecordBuffer {
    fn memory_usage(&self) -> usize {
        self.memory_usage()
    }

    fn allocated_capacity(&self) -> usize {
        self.allocated_capacity()
    }

    fn len(&self) -> usize {
        self.len()
    }

    fn num_segments(&self) -> usize {
        self.num_segments()
    }
}

// ============================================================================
// Key Extraction Functions
// ============================================================================

/// Extract coordinate sort key directly from raw BAM bytes.
///
/// This is a zero-copy operation that reads fields at fixed offsets.
///
/// For coordinate sorting (following samtools behavior):
/// - Reads with valid tid (>= 0) are sorted by (tid, pos, reverse)
/// - Reads with tid = -1 (no reference) sort at the end
/// - Unmapped reads with a valid tid use that tid for sorting (typically mate's position)
///
/// Note: No read name tie-breaking is used, matching samtools behavior.
/// Equal records maintain their original input order (stable sort).
#[inline]
#[must_use]
pub fn extract_coordinate_key_inline(bam: &[u8], nref: u32) -> u64 {
    let tid = fgumi_raw_bam::ref_id(bam);
    let pos = fgumi_raw_bam::pos(bam);
    let reverse = RawRecordView::new(bam).flags() & fgumi_raw_bam::flags::REVERSE != 0;

    // Pack key based on tid (samtools behavior):
    // - tid >= 0: sort by (tid, pos, reverse) even if unmapped flag is set
    // - tid < 0: unmapped with no reference, sort at end
    if tid < 0 {
        PackedCoordinateKey::unmapped().0
    } else {
        PackedCoordinateKey::new(tid, pos, reverse, nref).0
    }
}

// ============================================================================
// Template-Coordinate Support
// ============================================================================

/// Extended key for template-coordinate sorting.
///
/// Template-coordinate requires comparing multiple fields:
/// tid1, tid2, pos1, pos2, neg1, neg2, CB, library, MI, name, `is_upper`
///
/// We pack these into 5 u64 values for efficient comparison.
/// The `name_hash_upper` field packs both `name_hash` and `is_upper`:
/// - Upper 63 bits: name hash (groups same names together)
/// - Lowest bit: `is_upper` (false=0, true=1)
///
/// This ensures reads from the same template stay together (same hash),
/// with `is_upper=false` sorting before `is_upper=true`.
#[repr(C)]
#[derive(Copy, Clone, Debug, Eq, PartialEq, bytemuck::Pod, bytemuck::Zeroable)]
pub struct TemplateKey {
    /// Packed: (tid1 << 48) | (tid2 << 32) | pos1
    /// Comparison order matches samtools: tid1, tid2, pos1
    pub primary: u64,
    /// Packed: (pos2 << 32) | (!neg1 << 1) | !neg2
    /// neg flags inverted so reverse (neg=true) sorts before forward (neg=false)
    pub secondary: u64,
    /// Hash of the CB (cellular barcode) tag value.
    /// Inserted between neg2 and MI to match fgbio's sort order.
    /// 0 when no cell tag is present or `cell_tag` is disabled.
    pub cb_hash: u64,
    /// Packed: (library << 48) | (`mi_value` << 1) | `mi_suffix`
    pub tertiary: u64,
    /// Packed: (`name_hash` << 1) | `is_upper`
    /// This ensures same-name records group together, with `is_upper` as tie-breaker.
    pub name_hash_upper: u64,
}

impl TemplateKey {
    /// Create a new template key from extracted fields.
    #[allow(clippy::too_many_arguments, clippy::cast_sign_loss)]
    #[must_use]
    pub fn new(
        tid1: i32,
        pos1: i32,
        neg1: bool,
        tid2: i32,
        pos2: i32,
        neg2: bool,
        cb_hash: u64,
        library: u32,
        mi: (u64, bool),
        name_hash: u64,
        is_upper: bool,
    ) -> Self {
        // Handle i32::MAX specially (indicates unmapped mate)
        let tid1_packed = if tid1 == i32::MAX { 0xFFFF_u64 } else { (tid1.max(0) as u64) & 0xFFFF };
        let tid2_packed = if tid2 == i32::MAX { 0xFFFF_u64 } else { (tid2.max(0) as u64) & 0xFFFF };
        // Convert signed positions to unsigned preserving sort order (XOR with sign bit).
        // This ensures negative positions sort correctly before positive ones.
        // i32::MAX maps to 0xFFFFFFFF which sorts last (used for unmapped mate sentinel).
        let pos1_packed = u64::from((pos1 as u32) ^ 0x8000_0000) & 0xFFFF_FFFF;
        let pos2_packed = u64::from((pos2 as u32) ^ 0x8000_0000) & 0xFFFF_FFFF;

        // Pack primary: tid1 (bits 63-48), tid2 (bits 47-32), pos1 (bits 31-0)
        // This ensures comparison order matches samtools: tid1, tid2, pos1
        let p1 = (tid1_packed << 48) | (tid2_packed << 32) | pos1_packed;

        // Pack secondary: pos2 (bits 63-32), neg1 (bit 1), neg2 (bit 0)
        // Invert neg flags: samtools sorts reverse (neg=true) BEFORE forward (neg=false)
        // By storing !neg, we get: !true=0 < !false=1, so reverse sorts first
        let p2 = (pos2_packed << 32) | (u64::from(!neg1) << 1) | u64::from(!neg2);

        // Pack tertiary: library (high 16), mi_value (middle), mi_suffix (bit 0)
        // Note: /B suffix should sort after /A, so we use !is_a as the bit
        let p3 = ((u64::from(library) & 0xFFFF) << 48)
            | ((mi.0 & 0xFFFF_FFFF_FFFF) << 1)
            | u64::from(!mi.1);

        // Pack name_hash and is_upper: hash in upper 63 bits, is_upper in bit 0
        // This ensures same-name records group together, with is_upper=false before is_upper=true
        let p4 = (name_hash << 1) | u64::from(is_upper);

        Self { primary: p1, secondary: p2, cb_hash, tertiary: p3, name_hash_upper: p4 }
    }

    /// Create a key for completely unmapped records.
    ///
    /// `library` and `mi` are packed into `tertiary` identically to
    /// [`TemplateKey::new`] so a fully-unmapped read realizes the same library
    /// and MI lanes as its mapped, same-template peers. Zeroing `tertiary` here
    /// (the prior behavior) made an unmapped read carrying a valid RG realize
    /// library ordinal 0 while its mapped same-library peers realized ordinal N,
    /// which the auto `--key-types` dropped-lane verify then flagged as a
    /// spurious "library" violation on single-library inputs (#375).
    #[must_use]
    pub fn unmapped(
        name_hash: u64,
        cb_hash: u64,
        library: u32,
        mi: (u64, bool),
        is_read2: bool,
    ) -> Self {
        // Pack tertiary identically to `new`: library (high 16), mi_value
        // (middle), mi_suffix (bit 0).
        let tertiary = ((u64::from(library) & 0xFFFF) << 48)
            | ((mi.0 & 0xFFFF_FFFF_FFFF) << 1)
            | u64::from(!mi.1);
        Self {
            primary: u64::MAX,
            secondary: u64::MAX,
            cb_hash,
            tertiary,
            name_hash_upper: (name_hash << 1) | u64::from(is_read2),
        }
    }

    /// Create a zeroed key (used as dummy for memory operations).
    #[inline]
    #[must_use]
    pub fn zeroed() -> Self {
        Self { primary: 0, secondary: 0, cb_hash: 0, tertiary: 0, name_hash_upper: 0 }
    }
}

impl Default for TemplateKey {
    fn default() -> Self {
        Self::zeroed()
    }
}

impl TemplateKey {
    /// Serialize to bytes for storage in keyed temp files.
    #[inline]
    #[must_use]
    pub fn to_bytes(&self) -> [u8; 40] {
        let mut buf = [0u8; 40];
        buf[0..8].copy_from_slice(&self.primary.to_le_bytes());
        buf[8..16].copy_from_slice(&self.secondary.to_le_bytes());
        buf[16..24].copy_from_slice(&self.cb_hash.to_le_bytes());
        buf[24..32].copy_from_slice(&self.tertiary.to_le_bytes());
        buf[32..40].copy_from_slice(&self.name_hash_upper.to_le_bytes());
        buf
    }

    /// Deserialize from bytes read from keyed temp files.
    #[inline]
    #[must_use]
    pub fn from_bytes(buf: &[u8; 40]) -> Self {
        Self {
            primary: u64::from_le_bytes([
                buf[0], buf[1], buf[2], buf[3], buf[4], buf[5], buf[6], buf[7],
            ]),
            secondary: u64::from_le_bytes([
                buf[8], buf[9], buf[10], buf[11], buf[12], buf[13], buf[14], buf[15],
            ]),
            cb_hash: u64::from_le_bytes([
                buf[16], buf[17], buf[18], buf[19], buf[20], buf[21], buf[22], buf[23],
            ]),
            tertiary: u64::from_le_bytes([
                buf[24], buf[25], buf[26], buf[27], buf[28], buf[29], buf[30], buf[31],
            ]),
            name_hash_upper: u64::from_le_bytes([
                buf[32], buf[33], buf[34], buf[35], buf[36], buf[37], buf[38], buf[39],
            ]),
        }
    }
}

impl PartialOrd for TemplateKey {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for TemplateKey {
    fn cmp(&self, other: &Self) -> Ordering {
        self.primary
            .cmp(&other.primary)
            .then_with(|| self.secondary.cmp(&other.secondary))
            .then_with(|| self.cb_hash.cmp(&other.cb_hash))
            .then_with(|| self.tertiary.cmp(&other.tertiary))
            // name_hash_upper comparison handles both name grouping AND is_upper ordering
            .then_with(|| self.name_hash_upper.cmp(&other.name_hash_upper))
    }
}

impl TemplateKey {
    /// Compare only core fields (tid1, tid2, pos1, pos2, neg1, neg2, CB, library, MI).
    ///
    /// This ignores the `name_hash` tie-breaker, allowing verification to accept
    /// both fgumi and samtools sorted files (which differ only in tie-breaking).
    #[inline]
    #[must_use]
    pub fn core_cmp(&self, other: &Self) -> Ordering {
        self.primary
            .cmp(&other.primary)
            .then_with(|| self.secondary.cmp(&other.secondary))
            .then_with(|| self.cb_hash.cmp(&other.cb_hash))
            .then_with(|| self.tertiary.cmp(&other.tertiary))
    }
}

/// Full 5-lane template-coordinate key (alias of [`TemplateKey`]).
///
/// Used as the canonical full-width key from which narrow keys are derived and
/// against which decode-time verify compares. Identical layout / `Ord` / radix.
pub type TemplateKey40 = TemplateKey;

impl RawSortKey for TemplateKey {
    const SERIALIZED_SIZE: Option<usize> = Some(40);

    /// # Panics
    ///
    /// Always panics. `TemplateKey` extraction requires a [`LibraryLookup`](crate::external::LibraryLookup)
    /// context not available through the `RawSortKey` trait interface. All
    /// callers must use `extract_template_key_inline()` instead.
    fn extract(_bam: &[u8], _ctx: &SortContext) -> Self {
        unreachable!(
            "TemplateKey::extract() should not be called directly. \
             Use extract_template_key_inline() with LibraryLookup instead."
        )
    }

    #[inline]
    fn write_to<W: Write>(&self, writer: &mut W) -> std::io::Result<()> {
        writer.write_all(&self.to_bytes())
    }

    #[inline]
    fn read_from<R: Read>(reader: &mut R) -> std::io::Result<Self> {
        let mut buf = [0u8; 40];
        reader.read_exact(&mut buf)?;
        Ok(Self::from_bytes(&buf))
    }
}

/// Lite template-coordinate key: only the always-present lanes
/// (`primary`, `secondary`, `name_hash_upper`). 24 bytes.
///
/// Valid only when neither optional word (`cb_hash`, `tertiary`) varies across
/// records — verified at decode time. Lanes are in canonical sort precedence.
#[repr(C)]
#[derive(Copy, Clone, Debug, Eq, PartialEq, Default, bytemuck::Pod, bytemuck::Zeroable)]
pub struct TemplateKey24 {
    /// Packed coordinate primary lane (tid1<<48 | tid2<<32 | pos1).
    pub primary: u64,
    /// Packed secondary lane (pos2<<32 | !neg1<<1 | !neg2).
    pub secondary: u64,
    /// Packed (`name_hash`<<1 | `is_upper`).
    pub name_hash_upper: u64,
}

impl TemplateKey24 {
    /// Serialize to bytes for storage in keyed temp files.
    #[inline]
    #[must_use]
    pub fn to_bytes(self) -> [u8; 24] {
        let mut buf = [0u8; 24];
        buf[0..8].copy_from_slice(&self.primary.to_le_bytes());
        buf[8..16].copy_from_slice(&self.secondary.to_le_bytes());
        buf[16..24].copy_from_slice(&self.name_hash_upper.to_le_bytes());
        buf
    }

    /// Deserialize from bytes read from keyed temp files.
    #[inline]
    #[must_use]
    pub fn from_bytes(buf: &[u8; 24]) -> Self {
        Self {
            primary: u64::from_le_bytes([
                buf[0], buf[1], buf[2], buf[3], buf[4], buf[5], buf[6], buf[7],
            ]),
            secondary: u64::from_le_bytes([
                buf[8], buf[9], buf[10], buf[11], buf[12], buf[13], buf[14], buf[15],
            ]),
            name_hash_upper: u64::from_le_bytes([
                buf[16], buf[17], buf[18], buf[19], buf[20], buf[21], buf[22], buf[23],
            ]),
        }
    }

    /// Build a lite key from a full extracted key, dropping the optional lanes.
    /// Valid only when `cb_hash` and `tertiary` are constant across all records.
    #[inline]
    #[must_use]
    pub fn from_full(full: &TemplateKey40) -> Self {
        Self {
            primary: full.primary,
            secondary: full.secondary,
            name_hash_upper: full.name_hash_upper,
        }
    }

    /// Compare only the non-tie-breaker lanes (everything except `name_hash_upper`).
    #[inline]
    #[must_use]
    pub fn core_cmp(&self, other: &Self) -> Ordering {
        self.primary.cmp(&other.primary).then_with(|| self.secondary.cmp(&other.secondary))
    }
}

impl PartialOrd for TemplateKey24 {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for TemplateKey24 {
    fn cmp(&self, other: &Self) -> Ordering {
        self.primary
            .cmp(&other.primary)
            .then_with(|| self.secondary.cmp(&other.secondary))
            .then_with(|| self.name_hash_upper.cmp(&other.name_hash_upper))
    }
}

/// Single-optional template key: `primary`, `secondary`, one optional lane
/// (`opt` = `cb_hash` OR tertiary), `name_hash_upper`. 32 bytes.
///
/// Used when exactly one optional word is non-constant. The `opt` lane source
/// (`cb_hash` vs tertiary) is fixed by the constructor chosen at sort start; the
/// layout and `Ord` are identical for both — see `from_full_cb` / `from_full_tertiary` (T3).
#[repr(C)]
#[derive(Copy, Clone, Debug, Eq, PartialEq, Default, bytemuck::Pod, bytemuck::Zeroable)]
pub struct TemplateKey32 {
    /// Packed coordinate primary lane (tid1<<48 | tid2<<32 | pos1).
    pub primary: u64,
    /// Packed secondary lane (pos2<<32 | !neg1<<1 | !neg2).
    pub secondary: u64,
    /// The single retained optional lane (`cb_hash` or tertiary).
    pub opt: u64,
    /// Packed (`name_hash`<<1 | `is_upper`).
    pub name_hash_upper: u64,
}

impl TemplateKey32 {
    /// Serialize to bytes for storage in keyed temp files.
    #[inline]
    #[must_use]
    pub fn to_bytes(self) -> [u8; 32] {
        let mut buf = [0u8; 32];
        buf[0..8].copy_from_slice(&self.primary.to_le_bytes());
        buf[8..16].copy_from_slice(&self.secondary.to_le_bytes());
        buf[16..24].copy_from_slice(&self.opt.to_le_bytes());
        buf[24..32].copy_from_slice(&self.name_hash_upper.to_le_bytes());
        buf
    }

    /// Deserialize from bytes read from keyed temp files.
    #[inline]
    #[must_use]
    pub fn from_bytes(buf: &[u8; 32]) -> Self {
        Self {
            primary: u64::from_le_bytes([
                buf[0], buf[1], buf[2], buf[3], buf[4], buf[5], buf[6], buf[7],
            ]),
            secondary: u64::from_le_bytes([
                buf[8], buf[9], buf[10], buf[11], buf[12], buf[13], buf[14], buf[15],
            ]),
            opt: u64::from_le_bytes([
                buf[16], buf[17], buf[18], buf[19], buf[20], buf[21], buf[22], buf[23],
            ]),
            name_hash_upper: u64::from_le_bytes([
                buf[24], buf[25], buf[26], buf[27], buf[28], buf[29], buf[30], buf[31],
            ]),
        }
    }

    /// Build a single-optional key keeping `cb_hash` in the `opt` lane
    /// (single-cell, pre-group). Valid only when `tertiary` is constant.
    #[inline]
    #[must_use]
    pub fn from_full_cb(full: &TemplateKey40) -> Self {
        Self {
            primary: full.primary,
            secondary: full.secondary,
            opt: full.cb_hash,
            name_hash_upper: full.name_hash_upper,
        }
    }

    /// Build a single-optional key keeping `tertiary` (library<<48 | mi) in the
    /// `opt` lane (post-group / multi-library). Valid only when `cb_hash` is constant.
    #[inline]
    #[must_use]
    pub fn from_full_tertiary(full: &TemplateKey40) -> Self {
        Self {
            primary: full.primary,
            secondary: full.secondary,
            opt: full.tertiary,
            name_hash_upper: full.name_hash_upper,
        }
    }

    /// Compare only the non-tie-breaker lanes (everything except `name_hash_upper`).
    ///
    /// Includes `opt` because it is a core sort lane (either `cb_hash` or tertiary,
    /// both core fields in `TemplateKey40::core_cmp`).
    #[inline]
    #[must_use]
    pub fn core_cmp(&self, other: &Self) -> Ordering {
        self.primary
            .cmp(&other.primary)
            .then_with(|| self.secondary.cmp(&other.secondary))
            .then_with(|| self.opt.cmp(&other.opt))
    }
}

impl PartialOrd for TemplateKey32 {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for TemplateKey32 {
    fn cmp(&self, other: &Self) -> Ordering {
        self.primary
            .cmp(&other.primary)
            .then_with(|| self.secondary.cmp(&other.secondary))
            .then_with(|| self.opt.cmp(&other.opt))
            .then_with(|| self.name_hash_upper.cmp(&other.name_hash_upper))
    }
}

impl RawSortKey for TemplateKey24 {
    const SERIALIZED_SIZE: Option<usize> = Some(24);

    /// # Panics
    ///
    /// Always panics. Narrow template keys are derived from a full extracted
    /// key via `from_full`, never extracted directly; extraction is full-width
    /// only (see [`TemplateKey::extract`]).
    fn extract(_bam: &[u8], _ctx: &SortContext) -> Self {
        unreachable!(
            "TemplateKey24::extract() should not be called directly. \
             Narrow keys are derived from a full TemplateKey40 via from_full()."
        )
    }

    #[inline]
    fn write_to<W: Write>(&self, writer: &mut W) -> std::io::Result<()> {
        writer.write_all(&self.to_bytes())
    }

    #[inline]
    fn read_from<R: Read>(reader: &mut R) -> std::io::Result<Self> {
        let mut buf = [0u8; 24];
        reader.read_exact(&mut buf)?;
        Ok(Self::from_bytes(&buf))
    }
}

impl RawSortKey for TemplateKey32 {
    const SERIALIZED_SIZE: Option<usize> = Some(32);

    /// # Panics
    ///
    /// Always panics. Narrow template keys are derived from a full extracted
    /// key via `from_full`, never extracted directly; extraction is full-width
    /// only (see [`TemplateKey::extract`]).
    fn extract(_bam: &[u8], _ctx: &SortContext) -> Self {
        unreachable!(
            "TemplateKey32::extract() should not be called directly. \
             Narrow keys are derived from a full TemplateKey40 via from_full()."
        )
    }

    #[inline]
    fn write_to<W: Write>(&self, writer: &mut W) -> std::io::Result<()> {
        writer.write_all(&self.to_bytes())
    }

    #[inline]
    fn read_from<R: Read>(reader: &mut R) -> std::io::Result<Self> {
        let mut buf = [0u8; 32];
        reader.read_exact(&mut buf)?;
        Ok(Self::from_bytes(&buf))
    }
}

/// Lane-addressable template key usable by the generic radix and buffer.
///
/// `primary` is always lane 0 (radix phase 1). `remaining_fields` lists the
/// post-primary lanes in sort precedence, of arity `LANES - 1`.
pub trait TemplateLaneKey:
    Copy + Clone + Ord + Default + bytemuck::Pod + Send + Sync + 'static
{
    /// Number of u64 lanes (3, 4, or 5).
    const LANES: usize;
    /// The primary (first) lane — radix phase 1 always sorts on this.
    fn primary(&self) -> u64;
    /// Post-primary lanes, in sort precedence (length `LANES - 1`).
    fn remaining_fields() -> &'static [fn(&Self) -> u64];
    /// Build from a full extracted key (drops constant lanes).
    fn from_full(full: &TemplateKey40) -> Self;
}

impl TemplateLaneKey for TemplateKey24 {
    const LANES: usize = 3;

    #[inline]
    fn primary(&self) -> u64 {
        self.primary
    }

    fn remaining_fields() -> &'static [fn(&Self) -> u64] {
        const F: [fn(&TemplateKey24) -> u64; 2] = [|k| k.secondary, |k| k.name_hash_upper];
        &F
    }

    #[inline]
    fn from_full(full: &TemplateKey40) -> Self {
        TemplateKey24::from_full(full)
    }
}

impl TemplateLaneKey for TemplateKey32 {
    const LANES: usize = 4;

    #[inline]
    fn primary(&self) -> u64 {
        self.primary
    }

    fn remaining_fields() -> &'static [fn(&Self) -> u64] {
        const F: [fn(&TemplateKey32) -> u64; 3] =
            [|k| k.secondary, |k| k.opt, |k| k.name_hash_upper];
        &F
    }

    /// `TemplateKey32` is retained as the byte/`Ord`/radix substrate (and is
    /// exercised directly by the radix-sort tests). The production cb-vs-tertiary
    /// dispatch uses the `CbKey32` / `TertKey32` newtypes instead, so this
    /// `from_full` is not on the production sort path; it sources `opt` from
    /// `cb_hash` (`from_full_cb`) only so the trait is total. The radix/`Ord` are
    /// independent of which optional field `opt` carries.
    #[inline]
    fn from_full(full: &TemplateKey40) -> Self {
        TemplateKey32::from_full_cb(full)
    }
}

impl TemplateLaneKey for TemplateKey40 {
    const LANES: usize = 5;

    #[inline]
    fn primary(&self) -> u64 {
        self.primary
    }

    fn remaining_fields() -> &'static [fn(&Self) -> u64] {
        const F: [fn(&TemplateKey40) -> u64; 4] =
            [|k| k.secondary, |k| k.cb_hash, |k| k.tertiary, |k| k.name_hash_upper];
        &F
    }

    #[inline]
    fn from_full(full: &TemplateKey40) -> Self {
        *full
    }
}

// ============================================================================
// Single-optional newtype monomorphizations (CbKey32 / TertKey32)
//
// Both wrap `TemplateKey32` transparently and have identical byte layout, `Ord`,
// and radix behavior. They differ only in which optional lane the `opt` word is
// sourced from at extraction: `CbKey32::from_full` uses `cb_hash`, while
// `TertKey32::from_full` uses the tertiary (library<<48 | mi) word. Two distinct
// types let the dispatch in `sort_template_coordinate` pick the right source via
// monomorphization, with zero runtime cost (the inner key is `#[repr(C)]` Pod and
// the wrapper is `#[repr(transparent)]`, so `bytemuck::Pod` is sound to derive).
// ============================================================================

/// 32-byte single-optional template key whose `opt` lane carries `cb_hash`.
///
/// Used when exactly the `cb_hash` lane is non-constant (single-cell, pre-group).
#[repr(transparent)]
#[derive(
    Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Default, bytemuck::Pod, bytemuck::Zeroable,
)]
pub struct CbKey32(pub TemplateKey32);

/// 32-byte single-optional template key whose `opt` lane carries the tertiary
/// (library<<48 | mi) word.
///
/// Used when exactly the tertiary lane is non-constant (post-group / multi-library).
#[repr(transparent)]
#[derive(
    Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Default, bytemuck::Pod, bytemuck::Zeroable,
)]
pub struct TertKey32(pub TemplateKey32);

impl TemplateLaneKey for CbKey32 {
    const LANES: usize = 4;

    #[inline]
    fn primary(&self) -> u64 {
        self.0.primary
    }

    fn remaining_fields() -> &'static [fn(&Self) -> u64] {
        const F: [fn(&CbKey32) -> u64; 3] =
            [|k| k.0.secondary, |k| k.0.opt, |k| k.0.name_hash_upper];
        &F
    }

    #[inline]
    fn from_full(full: &TemplateKey40) -> Self {
        CbKey32(TemplateKey32::from_full_cb(full))
    }
}

impl TemplateLaneKey for TertKey32 {
    const LANES: usize = 4;

    #[inline]
    fn primary(&self) -> u64 {
        self.0.primary
    }

    fn remaining_fields() -> &'static [fn(&Self) -> u64] {
        const F: [fn(&TertKey32) -> u64; 3] =
            [|k| k.0.secondary, |k| k.0.opt, |k| k.0.name_hash_upper];
        &F
    }

    #[inline]
    fn from_full(full: &TemplateKey40) -> Self {
        TertKey32(TemplateKey32::from_full_tertiary(full))
    }
}

impl RawSortKey for CbKey32 {
    const SERIALIZED_SIZE: Option<usize> = Some(32);

    /// # Panics
    ///
    /// Always panics. Narrow template keys are derived from a full extracted
    /// key via `from_full`, never extracted directly.
    fn extract(_bam: &[u8], _ctx: &SortContext) -> Self {
        unreachable!(
            "CbKey32::extract() should not be called directly. \
             Narrow keys are derived from a full TemplateKey40 via from_full()."
        )
    }

    #[inline]
    fn write_to<W: Write>(&self, writer: &mut W) -> std::io::Result<()> {
        writer.write_all(&self.0.to_bytes())
    }

    #[inline]
    fn read_from<R: Read>(reader: &mut R) -> std::io::Result<Self> {
        let mut buf = [0u8; 32];
        reader.read_exact(&mut buf)?;
        Ok(CbKey32(TemplateKey32::from_bytes(&buf)))
    }
}

impl RawSortKey for TertKey32 {
    const SERIALIZED_SIZE: Option<usize> = Some(32);

    /// # Panics
    ///
    /// Always panics. Narrow template keys are derived from a full extracted
    /// key via `from_full`, never extracted directly.
    fn extract(_bam: &[u8], _ctx: &SortContext) -> Self {
        unreachable!(
            "TertKey32::extract() should not be called directly. \
             Narrow keys are derived from a full TemplateKey40 via from_full()."
        )
    }

    #[inline]
    fn write_to<W: Write>(&self, writer: &mut W) -> std::io::Result<()> {
        writer.write_all(&self.0.to_bytes())
    }

    #[inline]
    fn read_from<R: Read>(reader: &mut R) -> std::io::Result<Self> {
        let mut buf = [0u8; 32];
        reader.read_exact(&mut buf)?;
        Ok(TertKey32(TemplateKey32::from_bytes(&buf)))
    }
}

/// Inline header stored before each record in the data buffer.
///
/// The header records only the length of the BAM bytes that follow it. The
/// full sort key is cached in the parallel `TemplateRecordRef` index, which
/// is the only place the sort reads it from. Record bytes are accessed via
/// `offset + TEMPLATE_HEADER_SIZE`, so the header itself is never read back
/// during sorting or merge — only its 8-byte footprint matters.
///
/// Layout (8 bytes total):
/// - `record_len`: 4 bytes (u32)
/// - padding: 4 bytes (u32)
#[repr(C)]
#[derive(Copy, Clone, Debug)]
pub struct TemplateInlineHeader {
    /// Length of following raw BAM data.
    pub record_len: u32,
    /// Padding for 8-byte alignment.
    pub padding: u32,
}

/// Size of `TemplateInlineHeader` in bytes.
pub const TEMPLATE_HEADER_SIZE: usize = 8; // 4 (record_len) + 4 (padding)
const _: () = assert!(
    std::mem::size_of::<TemplateInlineHeader>() == TEMPLATE_HEADER_SIZE,
    "TEMPLATE_HEADER_SIZE must match size_of::<TemplateInlineHeader>()"
);

impl TemplateInlineHeader {
    /// Serialize header to a byte array.
    #[inline]
    #[must_use]
    #[allow(clippy::wrong_self_convention)]
    pub fn to_bytes(self) -> [u8; TEMPLATE_HEADER_SIZE] {
        let mut buf = [0u8; TEMPLATE_HEADER_SIZE];
        buf[0..4].copy_from_slice(&self.record_len.to_le_bytes());
        buf[4..8].copy_from_slice(&self.padding.to_le_bytes());
        buf
    }

    /// Read header from a byte slice.
    #[inline]
    #[must_use]
    pub fn read_from(data: &[u8]) -> Self {
        let record_len = u32::from_le_bytes([data[0], data[1], data[2], data[3]]);
        Self { record_len, padding: 0 }
    }
}

/// Record reference for template-coordinate sorting with cached key.
///
/// Generic over `K: TemplateLaneKey`; caches the sort key `K` inline for
/// O(1) comparisons during sort. The cached-key width varies with the chosen
/// key type: 24 bytes for `TemplateKey24`, 32 for `TemplateKey32`, 40 for
/// `TemplateKey40`. All comparison data is in the ref itself, avoiding random
/// access to the multi-GB data buffer during sorting.
#[repr(C)]
#[derive(Copy, Clone, Debug)]
pub struct TemplateRecordRef<K: TemplateLaneKey> {
    /// Cached sort key for O(1) comparisons without data buffer access.
    pub key: K,
    /// Offset to inline header in data buffer.
    pub offset: u64,
    /// Length of raw BAM data (excluding inline header).
    pub len: u32,
    /// Padding for alignment.
    pub padding: u32,
}

impl<K: TemplateLaneKey> PartialEq for TemplateRecordRef<K> {
    fn eq(&self, other: &Self) -> bool {
        self.offset == other.offset
    }
}

impl<K: TemplateLaneKey> Eq for TemplateRecordRef<K> {}

impl<K: TemplateLaneKey> PartialOrd for TemplateRecordRef<K> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<K: TemplateLaneKey> Ord for TemplateRecordRef<K> {
    fn cmp(&self, other: &Self) -> Ordering {
        // Note: This Ord impl is NOT used for sorting - we use sort_by with data access
        self.offset.cmp(&other.offset)
    }
}

/// Template-coordinate record buffer with inline headers, generic over the lane-key type `K`.
///
/// Inline headers (8 bytes each) store only record lengths; sort keys are cached
/// in the parallel `TemplateRecordRef<K>` index for O(1) comparison. Ref size
/// varies with `K`: 40 bytes for `TemplateKey24`, 48 for `TemplateKey32`,
/// 56 for `TemplateKey40`.
///
/// Memory layout:
/// ```text
/// data: [Header0][Record0][Header1][Record1]...
/// refs: [Ref0][Ref1]...
/// ```
///
/// During sorting, we use a custom comparator that:
/// 1. Compares sort keys cached in refs (fast, O(1))
/// 2. On ties, uses the full `K: Ord` ordering
pub struct TemplateRecordBuffer<K: TemplateLaneKey> {
    /// Segmented byte storage: inline headers + record data.
    data: SegmentedBuf,
    /// Minimal index for sorting.
    refs: Vec<TemplateRecordRef<K>>,
}

//
impl<K: TemplateLaneKey> TemplateRecordBuffer<K> {
    /// Create a new buffer with estimated capacity.
    #[must_use]
    pub fn with_capacity(estimated_records: usize, estimated_bytes: usize) -> Self {
        let header_bytes = estimated_records * TEMPLATE_HEADER_SIZE;
        Self {
            data: SegmentedBuf::with_capacity(estimated_bytes + header_bytes, SORT_SEGMENT_SIZE),
            refs: Vec::with_capacity(estimated_records),
        }
    }

    /// Push a record with a pre-computed template key.
    ///
    /// # Errors
    ///
    /// Returns an error if the record (plus header) exceeds the segment size
    /// (256 MiB) or if the record length exceeds `u32::MAX`.
    #[inline]
    pub fn push(&mut self, record: &[u8], key: K) -> anyhow::Result<()> {
        let total_bytes = TEMPLATE_HEADER_SIZE + record.len();
        anyhow::ensure!(
            total_bytes <= SORT_SEGMENT_SIZE,
            "BAM record of {} bytes (+ {} byte header) exceeds segment size of {} bytes; \
             this is likely a malformed BAM file",
            record.len(),
            TEMPLATE_HEADER_SIZE,
            SORT_SEGMENT_SIZE,
        );
        let record_len = u32::try_from(record.len())
            .map_err(|_| anyhow::anyhow!("record length {} exceeds u32::MAX", record.len()))?;

        // Reserve contiguous space for header + record
        let offset = self.data.reserve_contiguous(total_bytes) as u64;

        // Write inline header (8 bytes: record_len + padding). The full sort
        // key is cached only in the `TemplateRecordRef` index below (the live
        // reader); it is intentionally not stored in the inline header.
        let header = TemplateInlineHeader { record_len, padding: 0 };
        self.data.extend_in_place(&header.to_bytes());

        // Write raw BAM data
        self.data.extend_in_place(record);

        // Add ref with cached key for O(1) sort comparisons
        self.refs.push(TemplateRecordRef { key, offset, len: record_len, padding: 0 });
        Ok(())
    }

    /// Sort the index by cached key using stable LSD radix sort.
    ///
    /// Uses multi-field radix sort which is stable (preserves relative order
    /// of records with equal keys). This ensures deterministic output that
    /// matches samtools when records have identical sort keys.
    pub fn sort(&mut self) {
        radix_sort_template_refs(&mut self.refs);
    }

    /// Sort using parallel radix sort with stable k-way merge.
    ///
    /// Each chunk is sorted with stable radix sort, then merged with a
    /// heap that uses `chunk_idx` as tie-breaker to preserve input order.
    pub fn par_sort(&mut self) {
        parallel_radix_sort_template_refs(&mut self.refs);
    }

    /// Get record bytes by reference.
    #[inline]
    #[must_use]
    #[allow(clippy::cast_possible_truncation)]
    pub fn get_record(&self, r: &TemplateRecordRef<K>) -> &[u8] {
        self.data.slice(r.offset as usize + TEMPLATE_HEADER_SIZE, r.len as usize)
    }

    /// Iterate over sorted records.
    pub fn iter_sorted(&self) -> impl Iterator<Item = &[u8]> {
        self.refs.iter().map(|r| self.get_record(r))
    }

    /// Get the record references.
    #[must_use]
    pub fn refs(&self) -> &[TemplateRecordRef<K>] {
        &self.refs
    }

    /// Memory usage in bytes (actual data stored, not capacity).
    #[must_use]
    pub fn memory_usage(&self) -> usize {
        self.data.len() + self.refs.len() * std::mem::size_of::<TemplateRecordRef<K>>()
    }

    /// Total allocated capacity in bytes (data segments + refs Vec).
    #[must_use]
    pub fn allocated_capacity(&self) -> usize {
        self.data.allocated_capacity()
            + self.refs.capacity() * std::mem::size_of::<TemplateRecordRef<K>>()
    }

    /// Number of data segments in the underlying buffer.
    #[must_use]
    pub fn num_segments(&self) -> usize {
        self.data.num_segments()
    }

    /// Number of records.
    #[must_use]
    pub fn len(&self) -> usize {
        self.refs.len()
    }

    /// Check if buffer is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.refs.is_empty()
    }

    /// Get the key for a record reference (returns cached key from ref).
    #[inline]
    #[must_use]
    #[allow(clippy::unused_self)]
    pub fn get_key(&self, r: &TemplateRecordRef<K>) -> K {
        r.key
    }

    /// Iterate over sorted (key, record) pairs.
    /// Used for writing keyed temp chunks that preserve sort keys.
    pub fn iter_sorted_keyed(&self) -> impl Iterator<Item = (K, &[u8])> {
        self.refs.iter().map(|r| (self.get_key(r), self.get_record(r)))
    }

    /// Clear the buffer for reuse.
    pub fn clear(&mut self) {
        self.data.clear();
        self.refs.clear();
    }

    /// Sort in parallel and return each sub-array as a separate sorted chunk.
    ///
    /// Instead of merging the parallel sort sub-arrays back into one sorted
    /// buffer (as `par_sort` does), this returns each sub-array as its own
    /// `InMemoryChunk<TemplateKey>` so they can be passed as separate merge
    /// sources to the k-way merge, avoiding the intermediate merge step.
    ///
    /// When `threads <= 1`, returns a single chunk.
    ///
    /// After this returns, `self.data` and `self.refs` are both empty —
    /// the produced chunks share the original data via `Arc<SegmentedBuf>`.
    /// The caller must drop the buffer. `iter_sorted_keyed()` would yield
    /// zero records and `get_record()` called with a stale
    /// `TemplateRecordRef` would index-panic against the empty
    /// `SegmentedBuf`.
    pub(crate) fn par_sort_into_chunks(&mut self, threads: usize) -> Vec<InMemoryChunk<K>> {
        par_sort_into_chunks_impl!(
            self,
            threads,
            radix_sort_template_refs,
            TEMPLATE_HEADER_SIZE,
            |r: &TemplateRecordRef<K>| r.key
        )
    }

    /// Drain the (already-sorted) buffer into a single in-memory chunk
    /// whose records share an `Arc<SegmentedBuf>` backing store.
    ///
    /// See [`RecordBuffer::drain_into_single_chunk`] for usage and
    /// lifetime notes — both `data` and `refs` are cleared after this
    /// call. `iter_sorted_keyed()` would yield zero records and
    /// `get_record()` called with a stale `TemplateRecordRef` would
    /// index-panic against the empty `SegmentedBuf`.
    #[must_use]
    pub(crate) fn drain_into_single_chunk(&mut self) -> InMemoryChunk<K> {
        let data = Arc::new(std::mem::take(&mut self.data));
        let records: Vec<_> = self
            .refs
            .iter()
            .map(|r| (r.key, r.offset + TEMPLATE_HEADER_SIZE as u64, r.len))
            .collect();
        self.refs.clear();
        InMemoryChunk::from_parts(data, records)
    }
}

impl<K: TemplateLaneKey> ProbeableBuffer for TemplateRecordBuffer<K> {
    fn memory_usage(&self) -> usize {
        self.memory_usage()
    }

    fn allocated_capacity(&self) -> usize {
        self.allocated_capacity()
    }

    fn len(&self) -> usize {
        self.len()
    }

    fn num_segments(&self) -> usize {
        self.num_segments()
    }
}

// ============================================================================
// Radix Sort for RecordRef
// ============================================================================

// NOTE: QuerynameRecordBuffer was removed — it had zero callers and proved to
// regress when wired in (reconstruct_record allocates per record on write-out).

/// Threshold below which we use insertion sort instead of radix sort.
const RADIX_THRESHOLD: usize = 256;

/// Radix sort for `RecordRef` arrays using LSD (Least Significant Digit) approach.
///
/// Sorts by the `sort_key` field using 8-bit radix (256 buckets).
/// This is O(n×k) where k is the number of bytes to sort (typically 5-8).
///
/// Derives the maximum non-sentinel key and the byte-0 histogram in a single
/// fused counting pass, then dispatches to `radix_sort_sized` directly. See
/// [`radix_sort_record_refs_with_max`] for the (non-production) entry point
/// that accepts a caller-supplied bound instead.
///
/// That counting pass deliberately ignores `u64::MAX` keys (the unmapped
/// coordinate sentinel), which is what keeps `bytes_needed` at the ~5–6 bytes a mapped
/// coordinate key occupies instead of widening to the full 8 the moment a
/// single unmapped read appears — and coordinate-sorted BAMs essentially always
/// carry an unmapped tail. Ignoring it needs no contract from the caller:
/// `u64::MAX` is the largest `u64` *and* truncates to the all-`0xFF` maximum at
/// every radix width, so it sorts after every other key and stays stable among
/// its peers no matter how many passes run.
///
/// # Stability
/// Radix sort is inherently stable - records with equal keys maintain their
/// relative input order, matching samtools behavior.
pub fn radix_sort_record_refs(refs: &mut [RecordRef]) {
    let n = refs.len();
    if n < RADIX_THRESHOLD {
        insertion_sort_refs(refs);
        return;
    }

    // Derive the bound inside the byte-0 counting pass rather than in a separate
    // traversal. That pass already loads every `sort_key` to build its
    // histogram, so tracking the max costs a compare per record instead of a
    // second pass over the array -- and because it happens before any element is
    // scattered, the width is known in time to size the run.
    //
    // This is why `RecordBuffer` does not maintain a running bound across
    // pushes: doing so would buy the same saved traversal at the cost of a field
    // to keep in sync, and a bound that is wrong in release silently mis-sorts.
    // Here the bound is recomputed from the data every time and cannot drift.
    let mut counts = [0usize; 256];
    let mut max_key = 0u64;
    for r in refs.iter() {
        counts[(r.sort_key & 0xFF) as usize] += 1;
        if r.sort_key != u64::MAX {
            max_key = max_key.max(r.sort_key);
        }
    }

    radix_sort_sized(refs, bytes_needed_for(max_key), Some(&counts));
}

/// Number of radix byte-passes needed to order keys up to `max_key`.
///
/// Never zero: a bound of 0 does not mean the keys are all equal, since a zero
/// key can coexist with `u64::MAX` sentinels and those do not compare equal. One
/// pass separates them and is a stable no-op when the keys really are identical.
///
/// Widens by one byte when `max_key` exactly fills its byte count: a mapped
/// `max_key` that is all-`0xFF` at `bytes_needed` bytes (e.g. `0xFFFF` from
/// `tid = 0, pos = 32766, reverse`) truncates to the same value the `u64::MAX`
/// sentinel does at that width, so the two would tie across every pass and the
/// sentinel would be left in input order rather than sorting to the tail. Only
/// the max can collide — any smaller key would need all its low bytes set, which
/// would make it the max — so one extra byte forces a pass where the sentinel
/// (`0xFF`) and the key (`0x00`) differ. `bytes_needed` is at most 8 already (a
/// non-sentinel key is `< u64::MAX`), and the guard only fires below 8, so it
/// never overflows.
#[inline]
fn bytes_needed_for(max_key: u64) -> usize {
    let mut bytes_needed =
        if max_key == 0 { 1 } else { ((64 - max_key.leading_zeros()) as usize).div_ceil(8) };
    if bytes_needed < 8 && max_key == u64::MAX >> ((8 - bytes_needed) * 8) {
        bytes_needed += 1;
    }
    bytes_needed
}

/// Largest key in `refs` that is not the unmapped sentinel (`u64::MAX`), or 0
/// if there is no such key. See [`radix_sort_record_refs`] for why excluding
/// the sentinel preserves the sort order.
///
/// The sort itself derives this inside its byte-0 counting pass rather than in a
/// separate traversal; this standalone form remains for tests and benchmarks
/// that need the bound without running a sort.
#[inline]
fn max_non_sentinel_key(refs: &[RecordRef]) -> u64 {
    refs.iter().map(|r| r.sort_key).filter(|&k| k != u64::MAX).max().unwrap_or(0)
}

/// Radix sort for `RecordRef` arrays using a caller-supplied maximum key.
///
/// No production path uses this: [`radix_sort_record_refs`] derives the bound
/// from the data inside its first counting pass, which costs no extra traversal
/// and cannot disagree with the keys it is sorting. This entry point exists for
/// benchmarks and tests that need to pin a specific radix width -- passing
/// `u64::MAX`, for instance, forces the full 8 passes to measure what narrowing
/// saves. Prefer [`radix_sort_record_refs`] everywhere else: its bound carries
/// no precondition, whereas this one silently mis-sorts if `max_key` is wrong in
/// a release build.
///
/// `max_key` sizes the number of radix byte-passes (`bytes_needed`). Each key
/// must be **either `<= max_key` or exactly `u64::MAX`** (the unmapped
/// coordinate sentinel). Keys `<= max_key` fit within `bytes_needed` bytes and
/// are ordered correctly; `u64::MAX` keys truncate to all-`0xFF` under any
/// width, so they sort *after* every other key and remain stable among
/// themselves. This lets a coordinate sort size the passes from the maximum
/// *mapped* key (~5–6 bytes) instead of the full 8 even when unmapped reads are
/// present. A key strictly between `max_key` and `u64::MAX` would under-size the
/// byte count and mis-order records — debug builds assert against this.
#[allow(clippy::uninit_vec, unsafe_code)]
pub fn radix_sort_record_refs_with_max(refs: &mut [RecordRef], max_key: u64) {
    let n = refs.len();
    if n < RADIX_THRESHOLD {
        // Use insertion sort for small arrays
        insertion_sort_refs(refs);
        return;
    }

    debug_assert!(
        refs.iter().all(|r| r.sort_key <= max_key || r.sort_key == u64::MAX),
        "radix_sort_record_refs_with_max: a key exceeds max_key ({max_key}) without being the \
         unmapped sentinel (u64::MAX); bytes_needed would be too small and records would mis-sort",
    );

    radix_sort_sized(refs, bytes_needed_for(max_key), None);
}

/// Run `bytes_needed` LSD radix passes over `refs`.
///
/// `first_counts`, when supplied, is the byte-0 histogram the caller already
/// built (see [`radix_sort_record_refs`], which builds it while deriving the
/// bound). Reusing it keeps the fused bound derivation free of an extra
/// traversal; passing `None` just counts byte 0 here as every later pass does.
#[allow(clippy::uninit_vec, unsafe_code)]
fn radix_sort_sized(
    refs: &mut [RecordRef],
    bytes_needed: usize,
    first_counts: Option<&[usize; 256]>,
) {
    let n = refs.len();
    debug_assert!(n >= RADIX_THRESHOLD, "callers handle short slices with insertion sort");
    debug_assert!(bytes_needed >= 1, "a zero-width sort would leave zero keys behind sentinels");

    // Allocate auxiliary buffer
    let mut aux: Vec<RecordRef> = Vec::with_capacity(n);
    // SAFETY: `aux` is written exactly once per radix pass (the scatter loop
    // below) before any element is read; `RecordRef` is `Copy`/`Pod`, so leaving
    // it uninitialized until the first scatter is sound. See CLAUDE.md
    // "Approved hot-path unsafe".
    unsafe {
        aux.set_len(n);
    }

    // SAFETY: `src`/`dst` always point at the disjoint, properly-aligned
    // `[RecordRef]` storage of `refs`/`aux` (same length `n`); each pass writes
    // every `dst` slot exactly once (the scatter loop) before that buffer is read
    // as `src` on the next pass. See CLAUDE.md "Approved hot-path unsafe".
    let mut src = refs as *mut [RecordRef];
    let mut dst = aux.as_mut_slice() as *mut [RecordRef];
    let mut precomputed = first_counts.copied();

    // LSD radix sort - byte by byte from least significant
    for byte_idx in 0..bytes_needed {
        let src_slice = unsafe { &*src };
        let dst_slice = unsafe { &mut *dst };

        // Count occurrences of each byte value, unless the caller already did.
        let mut counts = precomputed.take().unwrap_or_else(|| {
            let mut counts = [0usize; 256];
            for r in src_slice {
                let byte = ((r.sort_key >> (byte_idx * 8)) & 0xFF) as usize;
                counts[byte] += 1;
            }
            counts
        });

        // Convert to cumulative offsets
        let mut total = 0;
        for count in &mut counts {
            let c = *count;
            *count = total;
            total += c;
        }

        // Scatter elements to destination (stable - preserves order within buckets)
        for r in src_slice {
            let byte = ((r.sort_key >> (byte_idx * 8)) & 0xFF) as usize;
            let dest_idx = counts[byte];
            counts[byte] += 1;
            dst_slice[dest_idx] = *r;
        }

        // Swap src and dst
        std::mem::swap(&mut src, &mut dst);
    }

    // If odd number of passes, copy back to original buffer
    if bytes_needed % 2 == 1 {
        let src_slice = unsafe { &*src };
        refs.copy_from_slice(src_slice);
    }
}

/// Parallel radix sort for `RecordRef` arrays.
///
/// Divides the array into chunks, sorts each chunk with radix sort,
/// then performs k-way merge. This provides near-linear speedup.
pub fn parallel_radix_sort_record_refs(refs: &mut [RecordRef]) {
    let max_key = max_non_sentinel_key(refs);
    parallel_radix_sort_record_refs_with_max(refs, max_key);
}

/// Parallel radix sort using a precomputed maximum key (see
/// [`radix_sort_record_refs_with_max`] for the `max_key` contract). The same
/// upper bound is handed to every chunk's sort, so no chunk re-scans for its
/// own maximum.
pub fn parallel_radix_sort_record_refs_with_max(refs: &mut [RecordRef], max_key: u64) {
    use rayon::prelude::*;

    let n = refs.len();
    if n < RADIX_THRESHOLD * 2 {
        // Small array - just use single-threaded radix sort
        radix_sort_record_refs_with_max(refs, max_key);
        return;
    }

    // Get number of threads from rayon
    let n_threads = rayon::current_num_threads();

    // For very large arrays, parallel chunked sort + merge is faster
    if n_threads > 1 && n > 10_000 {
        let chunk_size = n.div_ceil(n_threads);

        // Sort each chunk in parallel using radix sort. Every chunk's max key is
        // bounded by the whole-slice `max_key`, so the shared bound is safe.
        refs.par_chunks_mut(chunk_size).for_each(|chunk| {
            radix_sort_record_refs_with_max(chunk, max_key);
        });

        // K-way merge the sorted chunks
        // For simplicity, we use a heap-based merge into auxiliary storage
        let chunk_boundaries: Vec<_> = refs
            .chunks(chunk_size)
            .scan(0, |pos, chunk| {
                let start = *pos;
                *pos += chunk.len();
                Some(start..*pos)
            })
            .collect();

        merge_sorted_chunks(refs, &chunk_boundaries);
    } else {
        // Single-threaded radix sort
        radix_sort_record_refs_with_max(refs, max_key);
    }
}

/// Merge k sorted chunks in place using auxiliary storage.
fn merge_sorted_chunks(refs: &mut [RecordRef], chunk_ranges: &[std::ops::Range<usize>]) {
    use crate::radix::{heap_make, heap_sift_down};

    // Heap entry: (sort_key, chunk_idx, position_in_chunk)
    struct HeapEntry {
        key: u64,
        chunk_idx: usize,
        pos: usize,
    }

    if chunk_ranges.len() <= 1 {
        return;
    }

    let n = refs.len();
    let mut result: Vec<RecordRef> = Vec::with_capacity(n);

    // Initialize heap with first element from each chunk
    let mut heap: Vec<HeapEntry> = Vec::with_capacity(chunk_ranges.len());
    for (chunk_idx, range) in chunk_ranges.iter().enumerate() {
        if !range.is_empty() {
            heap.push(HeapEntry { key: refs[range.start].sort_key, chunk_idx, pos: range.start });
        }
    }

    if heap.is_empty() {
        return;
    }

    // Min-heap: smaller keys should be at the top
    // Use `chunk_idx` as tie-breaker for stability
    let lt = |a: &HeapEntry, b: &HeapEntry| -> bool { (a.key, a.chunk_idx) > (b.key, b.chunk_idx) };

    heap_make(&mut heap, &lt);
    let mut heap_size = heap.len();

    // Merge
    while heap_size > 0 {
        let entry = &heap[0];
        result.push(refs[entry.pos]);

        let chunk_idx = entry.chunk_idx;
        let next_pos = entry.pos + 1;
        let range = &chunk_ranges[chunk_idx];

        if next_pos < range.end {
            // More elements in this chunk
            heap[0] = HeapEntry { key: refs[next_pos].sort_key, chunk_idx, pos: next_pos };
            heap_sift_down(&mut heap, 0, heap_size, &lt);
        } else {
            // Chunk exhausted
            heap_size -= 1;
            if heap_size > 0 {
                heap.swap(0, heap_size);
                heap_sift_down(&mut heap, 0, heap_size, &lt);
            }
        }
    }

    // Copy result back
    refs.copy_from_slice(&result);
}

/// Binary insertion sort for small arrays of `RecordRef`.
fn insertion_sort_refs(refs: &mut [RecordRef]) {
    for i in 1..refs.len() {
        let key = refs[i].sort_key;
        let insert_pos = refs[..i].partition_point(|r| r.sort_key <= key);
        if insert_pos < i {
            refs[insert_pos..=i].rotate_right(1);
        }
    }
}

// ============================================================================
// Radix Sort for TemplateRecordRef (multi-field 4×u64 key)
// ============================================================================

/// MSD (Most Significant Digit) hybrid radix sort for `TemplateRecordRef` arrays.
///
/// Exploits the fact that the `primary` field (tid1/tid2/pos1) is highly discriminating:
/// most records have unique primary values in typical sequencing data.
///
/// Strategy:
/// 1. Radix sort by `primary` field — O(n × k₁) where k₁ = adaptive bytes needed
/// 2. Find runs of equal `primary` values (single O(n) scan)
/// 3. Sub-sort only the runs that need it:
///    - Size 1: already sorted (majority of records)
///    - Size 2–64: insertion sort by remaining fields
///    - Size >64: LSD radix sort remaining fields
///
/// This avoids full O(n) passes for fields 2–5 when most records are already
/// resolved by the primary field alone.
#[allow(clippy::uninit_vec, unsafe_code)]
pub fn radix_sort_template_refs<K: TemplateLaneKey>(refs: &mut [TemplateRecordRef<K>]) {
    let n = refs.len();
    if n < RADIX_THRESHOLD {
        insertion_sort_template_refs(refs);
        return;
    }

    // Allocate auxiliary buffer (reused across phases).
    let mut aux: Vec<TemplateRecordRef<K>> = Vec::with_capacity(n);
    // SAFETY: `aux` is written exactly once per radix pass (the scatter loop in
    // `radix_sort_template_field`) before any element is read; `TemplateRecordRef<K>`
    // is `Copy`/`Pod`, so leaving it uninitialized until the first scatter is sound.
    // (Unchanged invariant from the non-generic version — see CLAUDE.md
    // "Approved hot-path unsafe".)
    unsafe {
        aux.set_len(n);
    }

    // Phase 1: Radix sort by primary field. The lane accessor is a `fn(&K) -> u64`
    // pointer (non-capturing closure coercion), keeping the field-getter type fixed
    // across the recursive `sub_sort_runs` calls below — a capturing closure would
    // grow a fresh type per recursion level and blow the monomorphization recursion
    // limit.
    let primary: fn(&K) -> u64 = |k| k.primary();
    let max_primary = refs.iter().map(|r| primary(&r.key)).max().unwrap_or(0);
    let bytes_needed = bytes_needed_u64(max_primary);
    if bytes_needed > 0 {
        radix_sort_template_field(refs, &mut aux, primary, bytes_needed);
    }

    // Phase 2: Find equal-primary runs and sub-sort by remaining fields
    sub_sort_runs(refs, &mut aux, primary, K::remaining_fields());
}

/// Threshold for sub-sort runs: below this, use insertion sort.
const SUB_SORT_INSERTION_THRESHOLD: usize = 64;

/// Find runs of equal values for `run_field` and sub-sort each run by `remaining_fields`.
fn sub_sort_runs<K: TemplateLaneKey>(
    refs: &mut [TemplateRecordRef<K>],
    aux: &mut [TemplateRecordRef<K>],
    run_field: fn(&K) -> u64,
    remaining_fields: &[fn(&K) -> u64],
) {
    if remaining_fields.is_empty() {
        return;
    }

    let n = refs.len();
    let mut start = 0;
    while start < n {
        let val = run_field(&refs[start].key);
        let mut end = start + 1;
        while end < n && run_field(&refs[end].key) == val {
            end += 1;
        }

        let run = &mut refs[start..end];
        let run_len = run.len();
        if run_len > 1 {
            if run_len <= SUB_SORT_INSERTION_THRESHOLD {
                insertion_sort_template_refs(run);
            } else {
                let next_field = remaining_fields[0];
                let max_val = run.iter().map(|r| next_field(&r.key)).max().unwrap_or(0);
                let bytes_needed = bytes_needed_u64(max_val);
                if bytes_needed > 0 {
                    let run_aux = &mut aux[start..end];
                    radix_sort_template_field(run, run_aux, next_field, bytes_needed);
                }

                if remaining_fields.len() > 1 {
                    sub_sort_runs(run, &mut aux[start..end], next_field, &remaining_fields[1..]);
                }
            }
        }
        start = end;
    }
}

/// Radix sort a single u64 field of `TemplateRecordRef` using raw pointers.
#[allow(clippy::uninit_vec, unsafe_code)]
fn radix_sort_template_field<K: TemplateLaneKey>(
    refs: &mut [TemplateRecordRef<K>],
    aux: &mut [TemplateRecordRef<K>],
    get_field: fn(&K) -> u64,
    bytes_needed: usize,
) {
    let n = refs.len();

    // Use raw pointers to avoid borrow checker issues with swapping.
    // SAFETY: `src`/`dst` always point at the disjoint, properly-aligned
    // `[TemplateRecordRef<K>]` storage of `refs`/`aux` (same length `n`); each
    // pass writes every `dst` slot exactly once (the scatter loop) before that
    // buffer is read as `src` on the next pass. Unchanged invariant from the
    // non-generic version — see CLAUDE.md "Approved hot-path unsafe".
    let mut src = refs as *mut [TemplateRecordRef<K>];
    let mut dst = aux as *mut [TemplateRecordRef<K>];

    for byte_idx in 0..bytes_needed {
        let src_slice = unsafe { &*src };
        let dst_slice = unsafe { &mut *dst };

        // Count occurrences of each byte value
        let mut counts = [0usize; 256];
        for r in src_slice {
            let byte = ((get_field(&r.key) >> (byte_idx * 8)) & 0xFF) as usize;
            counts[byte] += 1;
        }

        // Convert to cumulative offsets
        let mut total = 0;
        for count in &mut counts {
            let c = *count;
            *count = total;
            total += c;
        }

        // Scatter elements to destination
        for item in src_slice.iter().take(n) {
            let byte = ((get_field(&item.key) >> (byte_idx * 8)) & 0xFF) as usize;
            let dest_idx = counts[byte];
            counts[byte] += 1;
            dst_slice[dest_idx] = *item;
        }

        // Swap src and dst
        std::mem::swap(&mut src, &mut dst);
    }

    // If odd number of passes, copy back to original buffer
    if bytes_needed % 2 == 1 {
        let src_slice = unsafe { &*src };
        refs.copy_from_slice(src_slice);
    }
}

/// Parallel radix sort for `TemplateRecordRef` arrays.
///
/// Uses stable radix sort per chunk with stable k-way merge (`chunk_idx` tie-breaker).
pub fn parallel_radix_sort_template_refs<K: TemplateLaneKey>(refs: &mut [TemplateRecordRef<K>]) {
    use rayon::prelude::*;

    let n = refs.len();
    if n < RADIX_THRESHOLD * 2 {
        radix_sort_template_refs(refs);
        return;
    }

    let n_threads = rayon::current_num_threads();

    if n_threads > 1 && n > 10_000 {
        let chunk_size = n.div_ceil(n_threads);

        // Sort each chunk in parallel
        refs.par_chunks_mut(chunk_size).for_each(|chunk| {
            radix_sort_template_refs(chunk);
        });

        // K-way merge the sorted chunks
        let chunk_boundaries: Vec<_> = refs
            .chunks(chunk_size)
            .scan(0, |pos, chunk| {
                let start = *pos;
                *pos += chunk.len();
                Some(start..*pos)
            })
            .collect();

        merge_sorted_template_chunks(refs, &chunk_boundaries);
    } else {
        radix_sort_template_refs(refs);
    }
}

/// Merge k sorted chunks of `TemplateRecordRef` in place.
///
/// Uses `chunk_idx` as tie-breaker for stable merge (preserves input order for equal keys).
fn merge_sorted_template_chunks<K: TemplateLaneKey>(
    refs: &mut [TemplateRecordRef<K>],
    chunk_ranges: &[std::ops::Range<usize>],
) {
    use crate::radix::{heap_make, heap_sift_down};

    struct HeapEntry<K> {
        key: K,
        chunk_idx: usize,
        pos: usize,
    }

    if chunk_ranges.len() <= 1 {
        return;
    }

    let n = refs.len();
    let mut result: Vec<TemplateRecordRef<K>> = Vec::with_capacity(n);

    let mut heap: Vec<HeapEntry<K>> = Vec::with_capacity(chunk_ranges.len());
    for (chunk_idx, range) in chunk_ranges.iter().enumerate() {
        if !range.is_empty() {
            heap.push(HeapEntry { key: refs[range.start].key, chunk_idx, pos: range.start });
        }
    }

    if heap.is_empty() {
        return;
    }

    // Min-heap with chunk_idx tie-breaker for stability
    let lt = |a: &HeapEntry<K>, b: &HeapEntry<K>| -> bool {
        match a.key.cmp(&b.key) {
            std::cmp::Ordering::Greater => true,
            std::cmp::Ordering::Less => false,
            std::cmp::Ordering::Equal => a.chunk_idx > b.chunk_idx,
        }
    };

    heap_make(&mut heap, &lt);
    let mut heap_size = heap.len();

    while heap_size > 0 {
        let entry = &heap[0];
        result.push(refs[entry.pos]);

        let chunk_idx = entry.chunk_idx;
        let next_pos = entry.pos + 1;
        let range = &chunk_ranges[chunk_idx];

        if next_pos < range.end {
            heap[0] = HeapEntry { key: refs[next_pos].key, chunk_idx, pos: next_pos };
            heap_sift_down(&mut heap, 0, heap_size, &lt);
        } else {
            heap_size -= 1;
            if heap_size > 0 {
                heap.swap(0, heap_size);
                heap_sift_down(&mut heap, 0, heap_size, &lt);
            }
        }
    }

    refs.copy_from_slice(&result);
}

/// Binary insertion sort for small arrays of `TemplateRecordRef`.
fn insertion_sort_template_refs<K: TemplateLaneKey>(refs: &mut [TemplateRecordRef<K>]) {
    for i in 1..refs.len() {
        let key = &refs[i].key;
        let insert_pos = refs[..i].partition_point(|r| r.key <= *key);
        if insert_pos < i {
            refs[insert_pos..=i].rotate_right(1);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_packed_coord_key_ordering() {
        // Lower tid should come first
        assert!(
            PackedCoordinateKey::new(0, 100, false, 10)
                < PackedCoordinateKey::new(1, 100, false, 10)
        );

        // Lower pos should come first
        assert!(
            PackedCoordinateKey::new(0, 100, false, 10)
                < PackedCoordinateKey::new(0, 200, false, 10)
        );

        // Forward should come before reverse (false < true)
        assert!(
            PackedCoordinateKey::new(0, 100, false, 10)
                < PackedCoordinateKey::new(0, 100, true, 10)
        );

        // Unmapped should come last
        assert!(PackedCoordinateKey::new(9, 1_000_000, true, 10) < PackedCoordinateKey::unmapped());
    }

    #[test]
    fn test_template_key_ordering() {
        let k1 = TemplateKey::new(0, 100, false, 0, 200, false, 0, 0, (1, true), 0, false);
        let k2 = TemplateKey::new(0, 100, false, 0, 200, false, 0, 0, (2, true), 0, false);
        assert!(k1 < k2);

        // /A suffix should come before /B
        let ka = TemplateKey::new(0, 100, false, 0, 200, false, 0, 0, (1, true), 0, false);
        let kb = TemplateKey::new(0, 100, false, 0, 200, false, 0, 0, (1, false), 0, false);
        assert!(ka < kb);

        // Same name hash: is_upper=false should come before is_upper=true
        let lower = TemplateKey::new(0, 100, false, 0, 200, false, 0, 0, (1, true), 12345, false);
        let upper = TemplateKey::new(0, 100, false, 0, 200, false, 0, 0, (1, true), 12345, true);
        assert!(lower < upper, "is_upper=false should sort before is_upper=true");

        // Different name hashes should group separately
        let first_hash_lo =
            TemplateKey::new(0, 100, false, 0, 200, false, 0, 0, (1, true), 100, false);
        let first_hash_hi =
            TemplateKey::new(0, 100, false, 0, 200, false, 0, 0, (1, true), 100, true);
        let second_hash =
            TemplateKey::new(0, 100, false, 0, 200, false, 0, 0, (1, true), 200, false);
        // first_hash records should come before second_hash records
        assert!(first_hash_lo < second_hash);
        assert!(first_hash_hi < second_hash);
    }

    #[test]
    fn test_radix_sort_record_refs() {
        // Create test refs with various keys
        let mut refs = vec![
            RecordRef { sort_key: 100, offset: 0, len: 10, padding: 0 },
            RecordRef { sort_key: 50, offset: 100, len: 10, padding: 0 },
            RecordRef { sort_key: 200, offset: 200, len: 10, padding: 0 },
            RecordRef { sort_key: 50, offset: 300, len: 10, padding: 0 }, // Duplicate key
            RecordRef { sort_key: 1, offset: 400, len: 10, padding: 0 },
        ];

        radix_sort_record_refs(&mut refs);

        // Verify sorted order
        assert_eq!(refs[0].sort_key, 1);
        assert_eq!(refs[1].sort_key, 50);
        assert_eq!(refs[2].sort_key, 50);
        assert_eq!(refs[3].sort_key, 100);
        assert_eq!(refs[4].sort_key, 200);

        // Verify stability: duplicate keys should maintain original order
        // offset=100 was before offset=300 in input, should remain so
        assert_eq!(refs[1].offset, 100);
        assert_eq!(refs[2].offset, 300);
    }

    #[test]
    #[allow(clippy::cast_sign_loss)]
    fn test_radix_sort_large() {
        // Test with larger array to trigger radix sort (> RADIX_THRESHOLD)
        let mut refs: Vec<RecordRef> = (0..1000)
            .map(|i| RecordRef {
                sort_key: (999 - i) as u64, // Reverse order
                offset: i as u64 * 100,
                len: 10,
                padding: 0,
            })
            .collect();

        radix_sort_record_refs(&mut refs);

        // Verify sorted
        for (i, r) in refs.iter().enumerate() {
            assert_eq!(r.sort_key, i as u64, "Expected sort_key {i} at index {i}");
        }
    }

    #[test]
    fn test_radix_sort_empty() {
        let mut refs: Vec<RecordRef> = Vec::new();
        radix_sort_record_refs(&mut refs);
        assert!(refs.is_empty());
    }

    #[test]
    fn test_radix_sort_single() {
        let mut refs = vec![RecordRef { sort_key: 42, offset: 0, len: 10, padding: 0 }];
        radix_sort_record_refs(&mut refs);
        assert_eq!(refs.len(), 1);
        assert_eq!(refs[0].sort_key, 42);
    }

    #[test]
    #[allow(clippy::cast_sign_loss)]
    fn test_radix_sort_all_same_keys() {
        // All same keys should maintain original order (stability)
        let mut refs: Vec<RecordRef> = (0..100)
            .map(|i| RecordRef { sort_key: 42, offset: i as u64 * 100, len: 10, padding: 0 })
            .collect();

        radix_sort_record_refs(&mut refs);

        // Verify all keys are 42 and order is preserved
        for (i, r) in refs.iter().enumerate() {
            assert_eq!(r.sort_key, 42);
            assert_eq!(r.offset, i as u64 * 100, "Stability violated at index {i}");
        }
    }

    #[test]
    #[allow(clippy::cast_sign_loss)]
    fn test_radix_sort_all_zero_keys() {
        let mut refs: Vec<RecordRef> = (0..50)
            .map(|i| RecordRef { sort_key: 0, offset: i as u64 * 100, len: 10, padding: 0 })
            .collect();

        radix_sort_record_refs(&mut refs);

        // Should be stable for all-zero keys
        for (i, r) in refs.iter().enumerate() {
            assert_eq!(r.sort_key, 0);
            assert_eq!(r.offset, i as u64 * 100);
        }
    }

    #[test]
    fn test_radix_sort_max_key() {
        // Test with maximum u64 values
        let mut refs = vec![
            RecordRef { sort_key: u64::MAX, offset: 0, len: 10, padding: 0 },
            RecordRef { sort_key: 0, offset: 100, len: 10, padding: 0 },
            RecordRef { sort_key: u64::MAX / 2, offset: 200, len: 10, padding: 0 },
        ];

        radix_sort_record_refs(&mut refs);

        assert_eq!(refs[0].sort_key, 0);
        assert_eq!(refs[1].sort_key, u64::MAX / 2);
        assert_eq!(refs[2].sort_key, u64::MAX);
    }

    #[test]
    #[allow(clippy::cast_sign_loss)]
    fn test_parallel_radix_sort() {
        // Test parallel sort with enough elements to trigger parallelism
        let mut refs: Vec<RecordRef> = (0..50_000)
            .map(|i| RecordRef {
                sort_key: (49_999 - i) as u64, // Reverse order
                offset: i as u64 * 100,
                len: 10,
                padding: 0,
            })
            .collect();

        parallel_radix_sort_record_refs(&mut refs);

        // Verify sorted
        for (i, r) in refs.iter().enumerate() {
            assert_eq!(r.sort_key, i as u64, "Expected sort_key {i} at index {i}");
        }
    }

    #[test]
    #[allow(clippy::cast_sign_loss)]
    fn test_parallel_radix_sort_stability() {
        // Test that parallel sort maintains stability for equal keys
        let mut refs: Vec<RecordRef> = (0..20_000)
            .map(|i| RecordRef {
                sort_key: (i / 100) as u64, // Groups of 100 with same key
                offset: i as u64,           // Use offset to track original order
                len: 10,
                padding: 0,
            })
            .collect();

        parallel_radix_sort_record_refs(&mut refs);

        // Verify sorted and stable within groups
        for i in 1..refs.len() {
            assert!(refs[i - 1].sort_key <= refs[i].sort_key, "Not sorted at index {i}");
            // Within same key group, offsets should be in order
            if refs[i - 1].sort_key == refs[i].sort_key {
                assert!(refs[i - 1].offset < refs[i].offset, "Stability violated at index {i}");
            }
        }
    }

    /// Build a `TemplateRecordBuffer<K>` of the given keys, sort it, and assert
    /// the resulting key order equals the keys sorted by `K: Ord`. Each record's
    /// payload is its insertion index so the buffer accepts a distinct push.
    fn assert_buffer_sorts<K: TemplateLaneKey + std::fmt::Debug>(keys: &[K]) {
        let mut buf = TemplateRecordBuffer::<K>::with_capacity(keys.len(), keys.len() * 8);
        for (i, k) in keys.iter().enumerate() {
            // Payload is the insertion index; its value is irrelevant to the sort
            // (only the cached key drives ordering) — it just makes each push distinct.
            let rec = (i as u64).to_le_bytes();
            buf.push(&rec, *k).expect("push");
        }
        buf.sort();
        let got: Vec<K> = buf.refs().iter().map(|r| r.key).collect();
        let mut expect = keys.to_vec();
        expect.sort();
        assert_eq!(got, expect, "radix sort must match Ord");
    }

    #[test]
    fn radix_sorts_template_key24() {
        // Large n (>= RADIX_THRESHOLD) exercises the radix path.
        // `name_hash_upper` uses a large odd constant so the lane spans the full
        // u64 width and drives an 8-byte radix pass (exercises `bytes_needed_u64`
        // for wide values on the narrow key type).
        let keys: Vec<TemplateKey24> = (0..5000u64)
            .map(|i| TemplateKey24 {
                primary: i.wrapping_mul(2_654_435_761) % 97,
                secondary: i % 13,
                name_hash_upper: i.wrapping_mul(0x9E37_79B9_7F4A_7C15),
            })
            .collect();
        assert_buffer_sorts(&keys);

        // Small n (< RADIX_THRESHOLD) exercises the insertion-sort path.
        let small: Vec<TemplateKey24> = (0..16u64)
            .map(|i| TemplateKey24 { primary: (17 - i) % 5, secondary: i % 3, name_hash_upper: i })
            .collect();
        assert_buffer_sorts(&small);
    }

    #[test]
    fn radix_sorts_template_key32() {
        // `name_hash_upper` uses a large odd constant so the lane spans the full
        // u64 width and drives an 8-byte radix pass (exercises `bytes_needed_u64`
        // for wide values on the narrow key type).
        let keys: Vec<TemplateKey32> = (0..5000u64)
            .map(|i| TemplateKey32 {
                primary: i.wrapping_mul(2_654_435_761) % 97,
                secondary: i % 13,
                opt: i.wrapping_mul(40503) % 31,
                name_hash_upper: i.wrapping_mul(0x9E37_79B9_7F4A_7C15),
            })
            .collect();
        assert_buffer_sorts(&keys);

        let small: Vec<TemplateKey32> = (0..16u64)
            .map(|i| TemplateKey32 {
                primary: (17 - i) % 5,
                secondary: i % 3,
                opt: (7 - i % 7),
                name_hash_upper: i,
            })
            .collect();
        assert_buffer_sorts(&small);
    }

    #[test]
    fn radix_sorts_template_key40() {
        // `name_hash_upper` uses a large odd constant so the lane spans the full
        // u64 width and drives an 8-byte radix pass.
        let keys: Vec<TemplateKey40> = (0..5000u64)
            .map(|i| TemplateKey40 {
                primary: i.wrapping_mul(2_654_435_761) % 97,
                secondary: i % 13,
                cb_hash: i.wrapping_mul(6_364_136_223_846_793_005) % 53,
                tertiary: i.wrapping_mul(40503) % 31,
                name_hash_upper: i.wrapping_mul(0x9E37_79B9_7F4A_7C15),
            })
            .collect();
        assert_buffer_sorts(&keys);

        let small: Vec<TemplateKey40> = (0..16u64)
            .map(|i| TemplateKey40 {
                primary: (17 - i) % 5,
                secondary: i % 3,
                cb_hash: (11 - i % 11),
                tertiary: (7 - i % 7),
                name_hash_upper: i,
            })
            .collect();
        assert_buffer_sorts(&small);
    }

    #[test]
    #[allow(clippy::cast_sign_loss)]
    fn test_radix_sort_template_refs_stability() {
        // Test that template radix sort is stable for equal keys
        // Create refs with identical TemplateKey but different offsets to track order
        let key = TemplateKey::new(0, 100, false, 0, 200, false, 0, 0, (1, true), 12345, false);

        let mut refs: Vec<TemplateRecordRef<TemplateKey>> = (0..500)
            .map(|i| TemplateRecordRef {
                key,              // All same key
                offset: i as u64, // Use offset to track original order
                len: 10,
                padding: 0,
            })
            .collect();

        radix_sort_template_refs(&mut refs);

        // All keys equal, so offsets should maintain original order (stability)
        for i in 1..refs.len() {
            assert!(
                refs[i - 1].offset < refs[i].offset,
                "Template radix sort stability violated at index {}: offset {} should be < {}",
                i,
                refs[i - 1].offset,
                refs[i].offset
            );
        }
    }

    #[test]
    #[allow(clippy::cast_sign_loss)]
    fn test_parallel_radix_sort_template_refs_stability() {
        // Test that parallel template sort maintains stability for equal keys
        // Use groups of records with same key to verify stability within groups
        let mut refs: Vec<TemplateRecordRef<TemplateKey>> = (0..20_000)
            .map(|i| {
                // Create groups of 100 with same key (different lib_name_hash for each group)
                let group = i / 100;
                let key = TemplateKey::new(
                    0,
                    100,
                    false,
                    0,
                    200,
                    false,
                    0,
                    0,
                    (1, true),
                    group as u64, // Different hash per group
                    false,
                );
                TemplateRecordRef {
                    key,
                    offset: i as u64, // Use offset to track original order
                    len: 10,
                    padding: 0,
                }
            })
            .collect();

        parallel_radix_sort_template_refs(&mut refs);

        // Verify sorted by key and stable within groups
        for i in 1..refs.len() {
            let prev_key = &refs[i - 1].key;
            let curr_key = &refs[i].key;
            assert!(prev_key <= curr_key, "Not sorted at index {i}");
            // Within same key group, offsets should be in order (stability)
            if prev_key == curr_key {
                assert!(
                    refs[i - 1].offset < refs[i].offset,
                    "Parallel template sort stability violated at index {}: offset {} should be < {}",
                    i,
                    refs[i - 1].offset,
                    refs[i].offset
                );
            }
        }
    }

    #[test]
    fn test_template_record_buffer_sort_stability() {
        // Test that TemplateRecordBuffer::sort() is stable
        // This is the method actually used during template-coordinate sorting
        let mut buffer = TemplateRecordBuffer::with_capacity(100, 10000);

        // Add records with identical keys - they should maintain insertion order after sort
        let key = TemplateKey::new(0, 100, false, 0, 200, false, 0, 0, (1, true), 12345, false);

        // Create distinct records (different sequence bytes) with same sort key
        for i in 0..100u8 {
            // Minimal valid BAM record: 4 bytes block_size prefix not included in our data
            // ref_id (4) + pos (4) + name_len (1) + mapq (1) + bin (2) + n_cigar (2) + flag (2) + seq_len (4)
            // + mate_ref_id (4) + mate_pos (4) + tlen (4) + name + cigar + seq + qual
            let mut record = vec![
                0, 0, 0, 0, // ref_id = 0
                100, 0, 0, 0, // pos = 100
                2, // name_len = 2 (including null)
                0, // mapq = 0
                0, 0, // bin
                0, 0, // n_cigar_op = 0
                99, 0, // flag = 99 (paired, proper, mate reverse, first)
                1, 0, 0, 0, // seq_len = 1
                0, 0, 0, 0, // mate_ref_id = 0
                200, 0, 0, 0, // mate_pos = 200
                0, 0, 0, 0, // tlen = 0
                b'A', 0, // read name "A\0"
                // no cigar
                i,    // seq (1 byte for 2 bases) - use i to make records distinguishable
                0xFF, // qual
            ];
            // Pad to make it valid
            while record.len() < 40 {
                record.push(0);
            }
            buffer.push(&record, key).expect("push should succeed in tests");
        }

        buffer.sort();

        // Verify records maintain insertion order (tracked by sequence byte value)
        let mut prev_seq_byte = None;
        for rec in buffer.iter_sorted() {
            // Extract the distinguishing byte (seq field at offset 34 after header parsing)
            let seq_byte = rec.get(34).copied().unwrap_or(0);
            if let Some(prev) = prev_seq_byte {
                assert!(
                    prev < seq_byte,
                    "TemplateRecordBuffer::sort() stability violated: {prev} should be < {seq_byte}"
                );
            }
            prev_seq_byte = Some(seq_byte);
        }
    }

    // ========================================================================
    // TemplateKey cb_hash tests
    // ========================================================================

    #[test]
    fn test_template_key_cb_hash_ordering() {
        // Same position but different cb_hash — lower hash sorts first
        let k1 = TemplateKey::new(0, 100, false, 0, 200, false, 10, 0, (1, true), 0, false);
        let k2 = TemplateKey::new(0, 100, false, 0, 200, false, 20, 0, (1, true), 0, false);
        assert!(k1 < k2, "lower cb_hash should sort before higher cb_hash");
    }

    #[test]
    fn test_template_key_cb_hash_zero_sorts_first() {
        // cb_hash=0 (no CB) should sort before non-zero cb_hash
        let k_no_cb = TemplateKey::new(0, 100, false, 0, 200, false, 0, 0, (1, true), 0, false);
        let k_cb = TemplateKey::new(0, 100, false, 0, 200, false, 42, 0, (1, true), 0, false);
        assert!(k_no_cb < k_cb, "cb_hash=0 should sort before non-zero cb_hash");
    }

    #[test]
    fn test_template_key_cb_hash_between_secondary_and_tertiary() {
        // Verify cb_hash is compared AFTER secondary but BEFORE tertiary.
        // Same primary+secondary, different cb_hash and tertiary values
        let k1 = TemplateKey::new(0, 100, false, 0, 200, false, 10, 0, (1, true), 0, false);
        let k2 = TemplateKey::new(0, 100, false, 0, 200, false, 20, 0, (0, true), 0, false);
        // k1 has lower cb_hash, so it should sort first regardless of tertiary
        assert!(k1 < k2, "cb_hash should be compared before tertiary (library/MI)");
    }

    #[test]
    fn test_template_key_serialization_with_cb_hash() {
        let key =
            TemplateKey::new(1, 500, true, 2, 600, false, 0xDEAD_BEEF, 3, (7, true), 999, false);
        let bytes = key.to_bytes();
        assert_eq!(bytes.len(), 40);
        let roundtrip = TemplateKey::from_bytes(&bytes);
        assert_eq!(key, roundtrip, "serialization roundtrip should preserve cb_hash");
        assert_eq!(roundtrip.cb_hash, 0xDEAD_BEEF);
    }

    #[test]
    fn test_template_key_unmapped_with_cb_hash() {
        let key = TemplateKey::unmapped(12345, 0xCAFE, 0, (0, true), false);
        assert_eq!(key.primary, u64::MAX);
        assert_eq!(key.secondary, u64::MAX);
        assert_eq!(key.cb_hash, 0xCAFE, "unmapped should preserve cb_hash");
        assert_eq!(key.tertiary, 0, "no library/MI -> zero tertiary");
    }

    /// `TemplateKey::unmapped` must pack the library ordinal (and MI) into
    /// `tertiary` exactly as `TemplateKey::new` does, so a fully-unmapped read
    /// realizes the same lanes as its mapped, same-library peers (#375).
    #[test]
    fn test_template_key_unmapped_packs_library_and_mi() {
        let unmapped = TemplateKey::unmapped(1, 0, 3, (7, true), false);
        assert_eq!(unmapped.tertiary >> 48, 3, "library ordinal occupies tertiary high-16");
        // Same library + MI fed to `new` must produce the same tertiary.
        let mapped = TemplateKey::new(0, 100, false, 0, 200, false, 0, 3, (7, true), 1, false);
        assert_eq!(
            unmapped.tertiary, mapped.tertiary,
            "unmapped must pack library/MI identically to TemplateKey::new"
        );
    }

    #[test]
    fn test_template_key_core_cmp_includes_cb_hash() {
        let k1 = TemplateKey::new(0, 100, false, 0, 200, false, 10, 0, (1, true), 0, false);
        let k2 = TemplateKey::new(0, 100, false, 0, 200, false, 20, 0, (1, true), 0, false);
        assert_eq!(k1.core_cmp(&k2), std::cmp::Ordering::Less, "core_cmp should include cb_hash");

        // Same cb_hash, different name_hash — core_cmp should be Equal
        let k3 = TemplateKey::new(0, 100, false, 0, 200, false, 10, 0, (1, true), 100, false);
        let k4 = TemplateKey::new(0, 100, false, 0, 200, false, 10, 0, (1, true), 200, false);
        assert_eq!(k3.core_cmp(&k4), std::cmp::Ordering::Equal, "core_cmp should ignore name_hash");
    }

    #[test]
    fn test_template_key_zeroed_has_zero_cb_hash() {
        let key = TemplateKey::zeroed();
        assert_eq!(key.cb_hash, 0);
    }

    #[test]
    fn test_template_key_default_has_zero_cb_hash() {
        let key = TemplateKey::default();
        assert_eq!(key.cb_hash, 0);
    }

    /// Helper: build a minimal raw BAM record byte array with a distinguishing index.
    fn make_bam_record(index: u16) -> Vec<u8> {
        let mut record = vec![
            0, 0, 0, 0, // ref_id = 0
            100, 0, 0, 0, // pos = 100
            2, // name_len = 2 (including null)
            0, // mapq = 0
            0, 0, // bin
            0, 0, // n_cigar_op = 0
            99, 0, // flag = 99 (paired, proper, mate reverse, first)
            1, 0, 0, 0, // seq_len = 1
            0, 0, 0, 0, // mate_ref_id = 0
            200, 0, 0, 0, // mate_pos = 200
            0, 0, 0, 0, // tlen = 0
            b'A', 0, // read name "A\0"
        ];
        // Encode the index into two bytes so each record is distinguishable
        record.push((index & 0xFF) as u8);
        record.push((index >> 8) as u8);
        record.push(0xFF); // qual
        // Pad to at least 40 bytes
        while record.len() < 40 {
            record.push(0);
        }
        record
    }

    /// Helper: build a minimal raw BAM record for coordinate sort with a specific position.
    fn make_coordinate_bam_record(tid: i32, pos: i32) -> Vec<u8> {
        let mut record = Vec::with_capacity(40);
        record.extend_from_slice(&tid.to_le_bytes()); // ref_id
        record.extend_from_slice(&pos.to_le_bytes()); // pos
        record.push(2); // name_len = 2 (including null)
        record.push(0); // mapq = 0
        record.extend_from_slice(&0u16.to_le_bytes()); // bin
        record.extend_from_slice(&0u16.to_le_bytes()); // n_cigar_op = 0
        record.extend_from_slice(&0u16.to_le_bytes()); // flags = 0 (forward)
        record.extend_from_slice(&1u32.to_le_bytes()); // seq_len = 1
        record.extend_from_slice(&(-1i32).to_le_bytes()); // mate_ref_id
        record.extend_from_slice(&(-1i32).to_le_bytes()); // mate_pos
        record.extend_from_slice(&0i32.to_le_bytes()); // tlen
        record.push(b'A'); // read name
        record.push(0); // null terminator
        // Pad to at least 34 bytes (fgumi_raw_bam::flags reads at offset 14-15)
        while record.len() < 40 {
            record.push(0);
        }
        record
    }

    #[test]
    fn test_record_buffer_excludes_unmapped_sentinel_and_sorts_correctly() {
        // Unmapped reads (tid < 0) carry the u64::MAX sentinel key. They must NOT
        // inflate the bound the radix derives (which would force it to the full
        // 8-byte width), yet must still sort to the tail. Exercise a mix above the
        // radix threshold so the real radix path (not insertion sort) runs.
        let nref = 5;
        let mut buffer = RecordBuffer::with_capacity(2048, 2048 * 64, nref);

        // 1500 mapped records with small keys (tid 0-3), interleaved with 500
        // unmapped (tid = -1 -> u64::MAX). Push order is intentionally scrambled.
        let mut expected_mapped_max = 0u64;
        for i in 0..2000u32 {
            if i % 4 == 3 {
                buffer.push_coordinate(&make_coordinate_bam_record(-1, -1)).expect("push unmapped");
            } else {
                let tid = i32::try_from(i % 3).unwrap();
                let pos = i32::try_from((i * 7) % 100_000).unwrap();
                buffer.push_coordinate(&make_coordinate_bam_record(tid, pos)).expect("push mapped");
                expected_mapped_max =
                    expected_mapped_max.max(PackedCoordinateKey::new(tid, pos, false, nref).0);
            }
        }

        // The derived bound excludes the sentinel: it equals the max mapped key,
        // far below u64::MAX (so bytes_needed stays ~3-4, not 8).
        let bound = max_non_sentinel_key(buffer.refs());
        assert_eq!(bound, expected_mapped_max);
        assert_ne!(bound, u64::MAX);
        assert!(bound < (1u64 << 48), "mapped max should be well under 8 bytes");

        // Sorting with the mapped-only max must still produce a fully ordered
        // result: mapped keys ascending, all unmapped (u64::MAX) at the tail.
        buffer.sort();
        let keys: Vec<u64> = buffer.refs().iter().map(|r| r.sort_key).collect();
        assert!(keys.windows(2).all(|w| w[0] <= w[1]), "result is not sorted by sort_key");
        let unmapped = keys.iter().filter(|&&k| k == u64::MAX).count();
        assert_eq!(unmapped, 500, "all unmapped reads should be present");
        assert!(
            keys[keys.len() - unmapped..].iter().all(|&k| k == u64::MAX),
            "unmapped reads must sort to the tail",
        );
    }

    /// `radix_sort_record_refs` and its parallel sibling scan for their own
    /// bound, and that scan must ignore the unmapped sentinel — otherwise a
    /// single unmapped read widens `bytes_needed` to the full 8 and the whole
    /// optimization is lost on exactly the inputs it targets (real coordinate
    /// BAMs, which always carry an unmapped tail). Sorting must stay correct
    /// either way, so assert both the ordering and the narrowed width.
    #[rstest::rstest]
    #[case::serial(false)]
    #[case::parallel(true)]
    fn test_scanning_radix_entry_points_ignore_unmapped_sentinel(#[case] parallel: bool) {
        // Well above RADIX_THRESHOLD * 2 so the real radix (and, for the
        // parallel case, the chunk-and-merge path) runs rather than insertion sort.
        let n = 20_000usize;
        let mut refs: Vec<RecordRef> = Vec::with_capacity(n);
        let mut expected_mapped_max = 0u64;
        for i in 0..n {
            // Every 4th record is unmapped; the rest get small, scrambled keys.
            let key = if i % 4 == 3 {
                u64::MAX
            } else {
                let k = ((i as u64 * 7919) % 100_000) + 1;
                expected_mapped_max = expected_mapped_max.max(k);
                k
            };
            refs.push(RecordRef::new(key, i as u64, 1));
        }
        let unmapped_count = refs.iter().filter(|r| r.sort_key == u64::MAX).count();

        // The bound the entry points derive must be the mapped max, not u64::MAX.
        assert_eq!(max_non_sentinel_key(&refs), expected_mapped_max);
        assert!(
            expected_mapped_max < (1u64 << 24),
            "mapped keys should need ~3 bytes, so the sentinel-inclusive max would cost 5 extra passes",
        );

        if parallel {
            parallel_radix_sort_record_refs(&mut refs);
        } else {
            radix_sort_record_refs(&mut refs);
        }

        let keys: Vec<u64> = refs.iter().map(|r| r.sort_key).collect();
        assert!(keys.windows(2).all(|w| w[0] <= w[1]), "result is not sorted by sort_key");
        assert_eq!(
            keys.iter().filter(|&&k| k == u64::MAX).count(),
            unmapped_count,
            "no unmapped record may be lost or duplicated",
        );
        assert!(
            keys[keys.len() - unmapped_count..].iter().all(|&k| k == u64::MAX),
            "unmapped records must sort to the tail",
        );
    }

    /// A buffer holding nothing but unmapped reads derives a bound of 0, which
    /// drives `bytes_needed` to 0 and takes the early return. That is a real
    /// input shape (an unmapped-only BAM), not a degenerate one, so it must come
    /// back with every record intact rather than panicking or dropping the tail.
    #[test]
    fn test_radix_sort_all_unmapped_records() {
        let n = 1_000usize; // above RADIX_THRESHOLD so the radix path is taken
        let mut refs: Vec<RecordRef> =
            (0..n).map(|i| RecordRef::new(u64::MAX, i as u64, 1)).collect();

        assert_eq!(max_non_sentinel_key(&refs), 0, "an all-unmapped buffer has no mapped bound");
        radix_sort_record_refs(&mut refs);

        assert_eq!(refs.len(), n, "no record may be dropped");
        assert!(refs.iter().all(|r| r.sort_key == u64::MAX), "all keys stay the sentinel");
        // Equal keys means the sort is a no-op, so the stable order is the input
        // order -- offsets must still read 0..n.
        assert!(
            refs.iter().enumerate().all(|(i, r)| r.offset == i as u64),
            "equal keys must preserve input order",
        );
    }

    /// A zero key mixed with unmapped sentinels is the one case where deriving
    /// the bound from mapped keys alone yields `max_key == 0` while the keys are
    /// *not* all equal. Skipping the sort there would leave the sentinels ahead
    /// of the zero keys. A zero key is reachable, not hypothetical:
    /// `PackedCoordinateKey::new` packs `tid = 0, pos = -1, reverse = false` to
    /// exactly 0, which a malformed record can carry.
    #[test]
    fn test_radix_sort_zero_keys_mixed_with_unmapped_sentinel() {
        // Confirm the packing really does produce a zero key, so this test keeps
        // tracking the reachable case rather than a synthetic one.
        assert_eq!(
            PackedCoordinateKey::new(0, -1, false, 25).0,
            0,
            "tid=0, pos=-1 should pack to a zero key",
        );

        // Above RADIX_THRESHOLD so the radix path runs rather than insertion sort.
        let n = 1_000usize;
        let mut refs: Vec<RecordRef> = (0..n)
            .map(|i| {
                // Alternate sentinel and zero keys, sentinels first, so an
                // unsorted result is immediately visible.
                let key = if i % 2 == 0 { u64::MAX } else { 0 };
                RecordRef::new(key, i as u64, 1)
            })
            .collect();

        assert_eq!(max_non_sentinel_key(&refs), 0, "the mapped bound here really is 0");
        radix_sort_record_refs(&mut refs);

        let keys: Vec<u64> = refs.iter().map(|r| r.sort_key).collect();
        assert!(
            keys.windows(2).all(|w| w[0] <= w[1]),
            "zero keys must sort ahead of the unmapped sentinels, not stay in input order",
        );
        assert_eq!(keys.iter().filter(|&&k| k == 0).count(), n / 2, "no zero-key record lost");
        assert_eq!(keys.iter().filter(|&&k| k == u64::MAX).count(), n / 2, "no sentinel lost");
    }

    /// `parallel_radix_sort_record_refs` falls back to the single-threaded sort
    /// when the slice is above the insertion-sort threshold but below the
    /// parallel cutoff. That branch has its own call into the bounded sort, so
    /// exercise it at a size that lands between the two.
    #[test]
    fn test_parallel_radix_sort_single_threaded_fallback_range() {
        // RADIX_THRESHOLD * 2 = 512 < n < 10_000 -> parallel cutoff not met.
        let n = 5_000usize;
        let mut refs: Vec<RecordRef> =
            (0..n).map(|i| RecordRef::new(((n - i) as u64) + 1, i as u64, 1)).collect();

        parallel_radix_sort_record_refs(&mut refs);

        let keys: Vec<u64> = refs.iter().map(|r| r.sort_key).collect();
        assert!(keys.windows(2).all(|w| w[0] <= w[1]), "fallback path must still sort");
        assert_eq!(keys.len(), n, "no record may be dropped");
    }

    /// The bound is now derived inside the byte-0 counting pass and its
    /// histogram is reused for that pass's scatter. Reusing a histogram across
    /// the width decision is the part that could go wrong, so this pins that the
    /// fused entry point agrees exactly -- ordering and stable tie-break -- with
    /// the same sort driven by an externally supplied bound.
    #[rstest::rstest]
    #[case::mixed_with_sentinels(5_000, 5, false)]
    #[case::no_sentinels(5_000, 0, false)]
    #[case::mostly_sentinels(5_000, 90, false)]
    #[case::just_over_threshold(RADIX_THRESHOLD + 1, 20, false)]
    #[case::wide_keys(20_000, 10, true)]
    fn test_fused_bound_matches_externally_supplied_bound(
        #[case] n: usize,
        #[case] sentinel_pct: usize,
        #[case] wide: bool,
    ) {
        let mut refs: Vec<RecordRef> = Vec::with_capacity(n);
        for i in 0..n {
            // Deterministic scramble so equal keys occur and stability is observable.
            let key = if sentinel_pct > 0 && (i * 37) % 100 < sentinel_pct {
                u64::MAX
            } else if wide {
                // Spread non-sentinel keys across the low 7 bytes so the derived
                // bound needs 4-8 radix passes, exercising `bytes_needed_for`'s
                // wide-width path (and its all-`0xFF` collision widening) rather
                // than the sub-2-byte range the other cases stay within.
                (i as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15) & 0x00FF_FFFF_FFFF_FFFF
            } else {
                ((i as u64 * 7919) % 50_000) + 1
            };
            refs.push(RecordRef::new(key, i as u64, 1));
        }

        let mut fused = refs.clone();
        radix_sort_record_refs(&mut fused);

        let mut supplied = refs.clone();
        radix_sort_record_refs_with_max(&mut supplied, max_non_sentinel_key(&refs));

        let as_pairs =
            |v: &[RecordRef]| v.iter().map(|r| (r.sort_key, r.offset)).collect::<Vec<_>>();
        assert_eq!(
            as_pairs(&fused),
            as_pairs(&supplied),
            "fused bound derivation changed the output for n={n}, sentinels={sentinel_pct}%",
        );

        // And it really is sorted, not merely consistent with the other path.
        assert!(fused.windows(2).all(|w| w[0].sort_key <= w[1].sort_key), "result is not sorted");
    }

    /// Sorting with a bound derived from the mapped keys must be
    /// indistinguishable from sorting at full 8-byte width. This is the
    /// output-identity check behind the pass-narrowing: the sentinel's
    /// truncation to all-`0xFF` has to leave the ordering *and* the stable
    /// tie-break order untouched, not merely sorted-looking.
    #[test]
    fn test_narrowed_radix_matches_full_width_output_exactly() {
        let n = 5_000usize;
        let mut refs: Vec<RecordRef> = Vec::with_capacity(n);
        for i in 0..n {
            // Deliberate duplicate keys so stability is observable: distinct
            // `offset` values act as the tie-break witness.
            let key = if i % 5 == 0 { u64::MAX } else { ((i as u64 * 37) % 500) + 1 };
            refs.push(RecordRef::new(key, i as u64, 1));
        }

        let mut narrowed = refs.clone();
        radix_sort_record_refs(&mut narrowed);

        let mut full_width = refs.clone();
        radix_sort_record_refs_with_max(&mut full_width, u64::MAX);

        let as_pairs =
            |v: &[RecordRef]| v.iter().map(|r| (r.sort_key, r.offset)).collect::<Vec<_>>();
        assert_eq!(
            as_pairs(&narrowed),
            as_pairs(&full_width),
            "narrowing the radix width changed the output; the sentinel must sort identically \
             at every width, including the stable order among equal keys",
        );
    }

    /// The one mapped key that collides with the `u64::MAX` sentinel under a
    /// narrowed radix is an all-`0xFF` key (`0xFF`, `0xFFFF`, `0xFFFFFF`, ...):
    /// it fills every bit of `bytes_needed` bytes, so it truncates to the same
    /// all-`0xFF` value the sentinel does at that width, and the two become
    /// indistinguishable across every pass. Such a key is ordinary — `0xFFFF`
    /// packs from `tid = 0, pos = 32766, reverse = true` — so a sentinel pushed
    /// ahead of it in input order would otherwise sort ahead of a mapped read.
    /// Placing the sentinels first makes any such mis-order visible, and the
    /// full-8-byte sort is an independent oracle for the correct output.
    #[rstest::rstest]
    #[case::one_byte(0xFFu64)]
    #[case::two_bytes(0xFFFFu64)]
    #[case::three_bytes(0xFF_FFFFu64)]
    fn test_narrowed_radix_separates_all_ones_key_from_sentinel(#[case] all_ones_key: u64) {
        // `0xFFFF` is a real coordinate key: tid=0, pos=32766, reverse=true.
        if all_ones_key == 0xFFFF {
            assert_eq!(
                PackedCoordinateKey::new(0, 32766, true, 25).0,
                0xFFFF,
                "0xFFFF must be a reachable mapped key, not a synthetic one",
            );
        }

        // Above RADIX_THRESHOLD so the real radix path runs, not insertion sort.
        let n = 1_000usize;
        let refs: Vec<RecordRef> = (0..n)
            .map(|i| {
                // Sentinels on even indices, the all-`0xFF` mapped key on odd:
                // every sentinel sits immediately before a colliding mapped key,
                // so leaving them in input order is unmistakably unsorted.
                let key = if i % 2 == 0 { u64::MAX } else { all_ones_key };
                RecordRef::new(key, i as u64, 1)
            })
            .collect();
        let sentinel_count = refs.iter().filter(|r| r.sort_key == u64::MAX).count();

        // The derived bound is exactly the all-`0xFF` mapped key, which is what
        // triggers the collision the fix guards against.
        assert_eq!(max_non_sentinel_key(&refs), all_ones_key);

        let mut narrowed = refs.clone();
        radix_sort_record_refs(&mut narrowed);

        // Full-width sort as an independent oracle.
        let mut full_width = refs.clone();
        radix_sort_record_refs_with_max(&mut full_width, u64::MAX);

        let as_pairs =
            |v: &[RecordRef]| v.iter().map(|r| (r.sort_key, r.offset)).collect::<Vec<_>>();
        assert_eq!(
            as_pairs(&narrowed),
            as_pairs(&full_width),
            "an all-0xFF mapped key ({all_ones_key:#x}) tied with the sentinel under the \
             narrowed radix; the two must be separated at the same output as a full-width sort",
        );

        let keys: Vec<u64> = narrowed.iter().map(|r| r.sort_key).collect();
        assert!(
            keys.windows(2).all(|w| w[0] <= w[1]),
            "result is not sorted: the all-0xFF mapped key and the sentinel were left tied",
        );
        assert!(
            keys[keys.len() - sentinel_count..].iter().all(|&k| k == u64::MAX),
            "every unmapped sentinel must sort to the tail, after the all-0xFF mapped key",
        );
    }

    /// Assert that `chunks` are each individually sorted and that the total
    /// record count across all chunks equals `expected_total`.
    fn assert_chunks_sorted_and_complete<K: RawSortKey + Default + Ord + 'static>(
        chunks: &[InMemoryChunk<K>],
        expected_total: usize,
    ) {
        for (ci, chunk) in chunks.iter().enumerate() {
            for i in 1..chunk.len() {
                assert!(
                    chunk.key_at(i - 1) <= chunk.key_at(i),
                    "chunk {ci} not sorted at index {i}"
                );
            }
        }
        let total: usize = chunks.iter().map(InMemoryChunk::len).sum();
        assert_eq!(total, expected_total, "total records across chunks should equal input count");
    }

    // `par_sort_into_chunks` returns exactly 1 chunk on the single-threaded
    // fallback (`threads == 1` or `n` below `RADIX_THRESHOLD * 2` / 10_000),
    // and > 1 chunks on the parallel path. We cover both ends across both
    // buffer types via `rstest`. The parallel cases must run inside a rayon
    // pool so `rayon::current_num_threads() > 1` inside the sort.
    //
    // (n, threads): n=100 stays under the parallel-path floor so we exercise
    // the fallback; n=10_500 clears both the > 10_000 and > RADIX_THRESHOLD*2
    // (512) gates.

    #[rstest::rstest]
    #[case::template_single_threaded(100, 1)]
    #[case::template_parallel(10_500, 4)]
    fn test_par_sort_into_chunks_template(#[case] n: usize, #[case] threads: usize) {
        let run = || {
            let mut buffer = TemplateRecordBuffer::with_capacity(n, n * 50);
            for i in 0..n {
                // Reverse-vary the primary sort field so sorting is non-trivial.
                let pos = i32::try_from(n - i).expect("test n fits in i32");
                let idx = u16::try_from(i).expect("test n fits in u16");
                let key = TemplateKey::new(0, pos, false, 0, 200, false, 0, 0, (1, true), 0, false);
                buffer.push(&make_bam_record(idx), key).expect("push should succeed in tests");
            }
            buffer.par_sort_into_chunks(threads)
        };

        let chunks = if threads > 1 {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .expect("failed to build rayon thread pool");
            pool.install(run)
        } else {
            run()
        };

        if threads == 1 {
            assert_eq!(chunks.len(), 1, "single-threaded should produce exactly 1 chunk");
            assert_eq!(chunks[0].len(), n, "the single chunk should contain all records");
        } else {
            assert!(
                chunks.len() > 1,
                "expected multiple chunks from parallel path, got {}",
                chunks.len()
            );
        }
        assert_chunks_sorted_and_complete(&chunks, n);
    }

    #[rstest::rstest]
    #[case::coordinate_single_threaded(100, 1)]
    #[case::coordinate_parallel(10_500, 4)]
    fn test_par_sort_into_chunks_coordinate(#[case] n: usize, #[case] threads: usize) {
        let nref = 10u32;
        let run = || {
            let mut buffer = RecordBuffer::with_capacity(n, n * 50, nref);
            for i in 0..n {
                // Reverse-vary the position so sorting is non-trivial.
                let pos = i32::try_from(n - i).expect("test n fits in i32");
                buffer
                    .push_coordinate(&make_coordinate_bam_record(0, pos))
                    .expect("push_coordinate should succeed in tests");
            }
            buffer.par_sort_into_chunks(threads)
        };

        let chunks = if threads > 1 {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .expect("failed to build rayon thread pool");
            pool.install(run)
        } else {
            run()
        };

        if threads == 1 {
            assert_eq!(chunks.len(), 1, "single-threaded should produce exactly 1 chunk");
            assert_eq!(chunks[0].len(), n, "the single chunk should contain all records");
        } else {
            assert!(
                chunks.len() > 1,
                "expected multiple chunks from parallel path, got {}",
                chunks.len()
            );
        }
        assert_chunks_sorted_and_complete(&chunks, n);
    }

    /// `par_sort_into_chunks` has each chunk derive its own bound inside its
    /// byte-0 counting pass, so the sentinel exclusion has to hold per chunk:
    /// each chunk must come back sorted and no unmapped record may be lost, on
    /// both the single-threaded early-return path and the multi-threaded path
    /// (which drain `refs` through different branches of the macro).
    #[rstest::rstest]
    #[case::single_threaded(100, 1)]
    #[case::parallel(10_500, 4)]
    fn test_par_sort_into_chunks_handles_unmapped_sentinel(
        #[case] n: usize,
        #[case] threads: usize,
    ) {
        let nref = 10u32;
        let mut buffer = RecordBuffer::with_capacity(n, n * 50, nref);
        let mut expected_unmapped = 0usize;
        for i in 0..n {
            // Every 4th record unmapped, so each chunk's bound must exclude the
            // sentinel or the chunks would all sort at full width.
            if i % 4 == 3 {
                buffer
                    .push_coordinate(&make_coordinate_bam_record(-1, -1))
                    .expect("push_coordinate should succeed in tests");
                expected_unmapped += 1;
            } else {
                let pos = i32::try_from(n - i).expect("test n fits in i32");
                buffer
                    .push_coordinate(&make_coordinate_bam_record(0, pos))
                    .expect("push_coordinate should succeed in tests");
            }
        }

        let chunks = if threads > 1 {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .expect("failed to build rayon thread pool");
            pool.install(|| buffer.par_sort_into_chunks(threads))
        } else {
            buffer.par_sort_into_chunks(threads)
        };

        assert_chunks_sorted_and_complete(&chunks, n);
        let unmapped: usize = chunks
            .iter()
            .map(|c| (0..c.len()).filter(|&i| c.key_at(i).sort_key == u64::MAX).count())
            .sum();
        assert_eq!(unmapped, expected_unmapped, "no unmapped record may be lost across chunks");
    }

    /// Helper: build a `SegmentedBuf` containing the given byte chunks
    /// concatenated, plus a vec of `(offset, len)` for each chunk.
    fn build_segmented_buf_with_records(records: &[&[u8]]) -> (SegmentedBuf, Vec<(u64, u32)>) {
        let total_bytes: usize = records.iter().map(|r| r.len()).sum();
        let mut buf = SegmentedBuf::with_capacity(total_bytes.max(1), total_bytes.max(1));
        let mut offsets = Vec::with_capacity(records.len());
        let mut cursor: u64 = 0;
        for r in records {
            buf.extend_in_place(r);
            offsets.push((cursor, u32::try_from(r.len()).unwrap()));
            cursor += r.len() as u64;
        }
        (buf, offsets)
    }

    #[test]
    fn test_in_memory_chunk_from_parts_reads_records() {
        let (buf, offsets) =
            build_segmented_buf_with_records(&[b"first", b"second-record", b"third"]);
        let data = Arc::new(buf);
        let records = vec![
            (10u32, offsets[0].0, offsets[0].1),
            (20u32, offsets[1].0, offsets[1].1),
            (30u32, offsets[2].0, offsets[2].1),
        ];
        let chunk = InMemoryChunk::<u32>::from_parts(data, records);

        assert_eq!(chunk.len(), 3);
        assert!(!chunk.is_empty());
        assert_eq!(chunk.record_bytes(0), b"first");
        assert_eq!(chunk.record_bytes(1), b"second-record");
        assert_eq!(chunk.record_bytes(2), b"third");
        assert_eq!(chunk.key_at(0), &10);
        assert_eq!(chunk.key_at(1), &20);
        assert_eq!(chunk.key_at(2), &30);
    }

    #[test]
    fn test_in_memory_chunk_take_key_replaces_with_default() {
        let (buf, offsets) = build_segmented_buf_with_records(&[b"payload"]);
        let chunk = InMemoryChunk::<u32>::from_parts(
            Arc::new(buf),
            vec![(42u32, offsets[0].0, offsets[0].1)],
        );
        let mut chunk = chunk;
        let k = chunk.take_key(0);
        assert_eq!(k, 42);
        // Subsequent take returns the default, which is what the merge
        // loop relies on once a record has been consumed.
        assert_eq!(chunk.take_key(0), 0);
    }

    #[test]
    fn test_in_memory_chunk_shared_arc_across_sibling_chunks() {
        // The whole point of the shared-buffer design: sibling chunks
        // from a single par_sort_into_chunks call share one
        // Arc<SegmentedBuf>, so the data isn't duplicated and is freed
        // only when the last chunk drops.
        let (buf, offsets) = build_segmented_buf_with_records(&[b"alpha", b"beta", b"gamma"]);
        let data = Arc::new(buf);
        let chunk_a = InMemoryChunk::<u32>::from_parts(
            Arc::clone(&data),
            vec![(1, offsets[0].0, offsets[0].1)],
        );
        let chunk_b = InMemoryChunk::<u32>::from_parts(
            Arc::clone(&data),
            vec![(2, offsets[1].0, offsets[1].1), (3, offsets[2].0, offsets[2].1)],
        );
        assert_eq!(Arc::strong_count(&data), 3, "data + 2 chunks share the Arc");
        assert_eq!(chunk_a.record_bytes(0), b"alpha");
        assert_eq!(chunk_b.record_bytes(0), b"beta");
        assert_eq!(chunk_b.record_bytes(1), b"gamma");

        drop(chunk_a);
        assert_eq!(Arc::strong_count(&data), 2);
        drop(chunk_b);
        assert_eq!(Arc::strong_count(&data), 1);
    }

    #[test]
    fn test_in_memory_chunk_empty_default() {
        let chunk = InMemoryChunk::<u32>::empty();
        assert_eq!(chunk.len(), 0);
        assert!(chunk.is_empty());
        // Default produces an empty chunk too.
        let chunk_default: InMemoryChunk<u32> = InMemoryChunk::default();
        assert!(chunk_default.is_empty());
    }

    #[test]
    #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
    fn test_drain_into_single_chunk_round_trips_record_bytes() {
        // End-to-end: insert records into a real RecordBuffer, sort,
        // drain into a single chunk, and verify that record_bytes(i)
        // returns the same bytes the original input held (offset +
        // HEADER_SIZE arithmetic must be correct).
        let nref = 4u32;
        let n: usize = 32;
        let mut buffer = RecordBuffer::with_capacity(n, n * 64, nref);

        let mut originals: Vec<Vec<u8>> = Vec::with_capacity(n);
        for i in 0..n {
            // Reverse position so the sort is non-trivial.
            let rec = make_coordinate_bam_record(0, (n - i) as i32);
            originals.push(rec.clone());
            buffer.push_coordinate(&rec).expect("push_coordinate should succeed");
        }

        buffer.par_sort();
        let chunk = buffer.drain_into_single_chunk();

        // Buffer's refs and data must both be empty after the drain.
        assert!(buffer.refs().is_empty(), "refs cleared after drain");
        assert!(buffer.is_empty());

        assert_eq!(chunk.len(), n);

        // The chunk contains the same byte multisets as the input. Use
        // a sorted compare since the input order was reversed.
        let mut chunk_bytes: Vec<Vec<u8>> =
            (0..chunk.len()).map(|i| chunk.record_bytes(i).to_vec()).collect();
        chunk_bytes.sort();
        originals.sort();
        assert_eq!(chunk_bytes, originals, "chunk records round-trip input bytes");
    }

    #[test]
    #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
    fn test_par_sort_into_chunks_shares_arc_across_siblings() {
        // Multi-thread path: verify all returned chunks share the same
        // Arc<SegmentedBuf>, i.e. the macro really does avoid duplicating
        // the source data into per-chunk owned buffers.
        let nref = 4u32;
        let n: usize = 10_500;
        let pool =
            rayon::ThreadPoolBuilder::new().num_threads(4).build().expect("rayon pool build");

        let chunks = pool.install(|| {
            let mut buffer = RecordBuffer::with_capacity(n, n * 64, nref);
            for i in 0..n {
                let rec = make_coordinate_bam_record(0, (n - i) as i32);
                buffer.push_coordinate(&rec).expect("push_coordinate");
            }
            buffer.par_sort_into_chunks(4)
        });

        assert!(chunks.len() > 1, "expected multiple chunks");

        // All chunks must point at the same Arc. ptr_eq compares the
        // Arc allocation address.
        let first_data = Arc::as_ptr(&chunks[0].data);
        for (i, c) in chunks.iter().enumerate().skip(1) {
            assert!(
                std::ptr::eq(Arc::as_ptr(&c.data), first_data),
                "chunk {i} should share the same Arc as chunk 0"
            );
        }
        // Strong count = number of chunks (after the macro's local `data`
        // binding has been dropped).
        assert_eq!(Arc::strong_count(&chunks[0].data), chunks.len());
    }

    mod template_key32 {
        use super::*;

        #[test]
        fn template_key32_ord_is_lexicographic_over_lanes() {
            let a = TemplateKey32 { primary: 1, secondary: 1, opt: 1, name_hash_upper: 1 };
            let b = TemplateKey32 { primary: 1, secondary: 1, opt: 1, name_hash_upper: 2 };
            let c = TemplateKey32 { primary: 1, secondary: 1, opt: 2, name_hash_upper: 0 };
            let d = TemplateKey32 { primary: 1, secondary: 2, opt: 0, name_hash_upper: 0 };
            assert!(a < b && b < c && c < d);
        }

        #[test]
        fn template_key32_to_from_bytes_roundtrip() {
            let k = TemplateKey32 {
                primary: 0xDEAD_BEEF_CAFE_1234,
                secondary: 0xABCD_0001_0002_0003,
                opt: 0x0102_0304_0506_0708,
                name_hash_upper: 0xFFFF_DEAD_BEEF_0099,
            };
            assert_eq!(TemplateKey32::from_bytes(&k.to_bytes()), k);
        }

        #[test]
        fn template_key32_core_cmp_includes_opt_excludes_name_hash() {
            let a = TemplateKey32 { primary: 1, secondary: 1, opt: 1, name_hash_upper: 5 };
            let b = TemplateKey32 { primary: 1, secondary: 1, opt: 2, name_hash_upper: 0 };
            let c = TemplateKey32 { primary: 1, secondary: 1, opt: 1, name_hash_upper: 9 };
            assert_eq!(a.core_cmp(&b), std::cmp::Ordering::Less, "opt is a core lane");
            assert_eq!(a.core_cmp(&c), std::cmp::Ordering::Equal, "name_hash ignored in core_cmp");
        }
    }

    mod template_key24 {
        use super::*;

        #[test]
        fn template_key24_ord_is_lexicographic_over_lanes() {
            let a = TemplateKey24 { primary: 1, secondary: 5, name_hash_upper: 9 };
            let b = TemplateKey24 { primary: 1, secondary: 5, name_hash_upper: 10 };
            let c = TemplateKey24 { primary: 1, secondary: 6, name_hash_upper: 0 };
            let d = TemplateKey24 { primary: 2, secondary: 0, name_hash_upper: 0 };
            assert!(a < b, "name_hash breaks ties");
            assert!(b < c, "secondary dominates name_hash");
            assert!(c < d, "primary dominates all");
            assert_eq!(a.cmp(&a), std::cmp::Ordering::Equal);
        }

        #[test]
        fn template_key24_to_from_bytes_roundtrip() {
            let k = TemplateKey24 {
                primary: 0xDEAD_BEEF_CAFE_1234,
                secondary: 0xABCD_0001_0002_0003,
                name_hash_upper: 0xFFFF_DEAD_BEEF_0099,
            };
            assert_eq!(TemplateKey24::from_bytes(&k.to_bytes()), k);
        }

        #[test]
        fn template_key24_core_cmp_ignores_name_hash() {
            let a = TemplateKey24 { primary: 1, secondary: 2, name_hash_upper: 100 };
            let b = TemplateKey24 { primary: 1, secondary: 2, name_hash_upper: 200 };
            assert_eq!(a.core_cmp(&b), std::cmp::Ordering::Equal);
        }
    }

    mod from_full {
        use super::*;

        #[test]
        fn template_key24_from_full_drops_optional_lanes() {
            let full = TemplateKey40 {
                primary: 11,
                secondary: 22,
                cb_hash: 99,
                tertiary: 88,
                name_hash_upper: 33,
            };
            let lite = TemplateKey24::from_full(&full);
            assert_eq!(lite, TemplateKey24 { primary: 11, secondary: 22, name_hash_upper: 33 });
        }

        #[test]
        fn template_key32_from_full_cb_uses_cb_hash_as_opt() {
            let full = TemplateKey40 {
                primary: 1,
                secondary: 2,
                cb_hash: 7,
                tertiary: 0,
                name_hash_upper: 4,
            };
            let k = TemplateKey32::from_full_cb(&full);
            assert_eq!(k, TemplateKey32 { primary: 1, secondary: 2, opt: 7, name_hash_upper: 4 });
        }

        #[test]
        fn template_key32_from_full_tertiary_uses_tertiary_as_opt() {
            let full = TemplateKey40 {
                primary: 1,
                secondary: 2,
                cb_hash: 0,
                tertiary: 55,
                name_hash_upper: 4,
            };
            let k = TemplateKey32::from_full_tertiary(&full);
            assert_eq!(k, TemplateKey32 { primary: 1, secondary: 2, opt: 55, name_hash_upper: 4 });
        }

        #[test]
        fn narrow_ord_matches_full_when_dropped_lanes_constant() {
            let recs = [
                TemplateKey40 {
                    primary: 2,
                    secondary: 0,
                    cb_hash: 5,
                    tertiary: 5,
                    name_hash_upper: 0,
                },
                TemplateKey40 {
                    primary: 1,
                    secondary: 9,
                    cb_hash: 5,
                    tertiary: 5,
                    name_hash_upper: 1,
                },
                TemplateKey40 {
                    primary: 1,
                    secondary: 9,
                    cb_hash: 5,
                    tertiary: 5,
                    name_hash_upper: 0,
                },
                TemplateKey40 {
                    primary: 1,
                    secondary: 1,
                    cb_hash: 5,
                    tertiary: 5,
                    name_hash_upper: 9,
                },
            ];
            for a in &recs {
                for b in &recs {
                    assert_eq!(
                        TemplateKey24::from_full(a).cmp(&TemplateKey24::from_full(b)),
                        a.cmp(b),
                        "lite ordering must match full when cb/tertiary constant",
                    );
                }
            }
        }
    }

    mod proptest_msd {
        use super::*;
        use proptest::{prop_assert_eq, proptest};
        use std::collections::hash_map::DefaultHasher;
        use std::hash::{Hash, Hasher};

        fn make_ref(key: TemplateKey, offset: u64) -> TemplateRecordRef<TemplateKey> {
            TemplateRecordRef { key, offset, len: 10, padding: 0 }
        }

        fn hash_pair(a: u64, b: u64) -> u64 {
            let mut h = DefaultHasher::new();
            (a, b).hash(&mut h);
            h.finish()
        }

        proptest! {
            /// Oracle test: MSD hybrid sort must produce the same ordering as a
            /// reference sort-by-key on randomized inputs with heavy primary collisions.
            #[test]
            fn msd_sort_matches_reference(
                n_primaries in 1_usize..=4,
                seed in proptest::num::u64::ANY,
            ) {
                let primaries: Vec<u64> = (0..n_primaries)
                    .map(|i| hash_pair(seed, i as u64))
                    .collect();

                // Build 300+ refs (above RADIX_THRESHOLD) with random keys sharing primaries
                let n = 300;
                let mut refs: Vec<TemplateRecordRef<TemplateKey>> = Vec::with_capacity(n);
                for i in 0..n {
                    let h = hash_pair(seed, (i + n_primaries) as u64);
                    let primary = primaries[i % n_primaries];
                    let key = TemplateKey {
                        primary,
                        secondary: h,
                        cb_hash: h.wrapping_mul(2_654_435_761),
                        tertiary: h.wrapping_mul(40503),
                        name_hash_upper: h.rotate_left(17),
                    };
                    refs.push(make_ref(key, i as u64));
                }

                let mut expected = refs.clone();
                expected.sort_by(|a, b| a.key.cmp(&b.key));

                radix_sort_template_refs(&mut refs);

                for i in 0..n {
                    prop_assert_eq!(
                        refs[i].key, expected[i].key,
                        "Mismatch at index {}: MSD key {:?} != reference {:?}",
                        i, refs[i].key, expected[i].key
                    );
                }
            }

            /// Oracle test with fully random keys (no forced primary collisions).
            #[test]
            fn msd_sort_matches_reference_random_keys(
                seed in proptest::num::u64::ANY,
            ) {
                let n = 500;
                let mut refs: Vec<TemplateRecordRef<TemplateKey>> = Vec::with_capacity(n);
                for i in 0..n {
                    let h = hash_pair(seed, i as u64);
                    let key = TemplateKey {
                        primary: h,
                        secondary: h.wrapping_mul(6_364_136_223_846_793_005),
                        cb_hash: h.wrapping_mul(2_654_435_761),
                        tertiary: h.wrapping_mul(40503),
                        name_hash_upper: h.rotate_left(17),
                    };
                    refs.push(make_ref(key, i as u64));
                }

                let mut expected = refs.clone();
                expected.sort_by(|a, b| a.key.cmp(&b.key));

                radix_sort_template_refs(&mut refs);

                for i in 0..n {
                    prop_assert_eq!(
                        refs[i].key, expected[i].key,
                        "Mismatch at index {}: MSD key {:?} != reference {:?}",
                        i, refs[i].key, expected[i].key
                    );
                }
            }
        }
    }
}

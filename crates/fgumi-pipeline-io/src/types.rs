//! Concrete data types that flow through the BAM step library.
//!
//! Every flowing type carries an explicit `batch_serial: u64` field and
//! impls both [`HeapSize`] (so byte-bounded queues can budget memory) and
//! [`Ordered`] (so `BranchOrdering::ByItemOrdinal` reorder stages preserve
//! global ordering across multi-step Parallel transforms).

use fgumi_pipeline_core::{HeapSize, Ordered};
use fgumi_raw_bam::RawRecord;

// ─────────────────────────────────────────────────────────────────────────────
// BamIndexManifest / RecordIndexEntry — inline BAI sidecar metadata.
// ─────────────────────────────────────────────────────────────────────────────

/// One record's virtual-offset bookkeeping within a compressed [`BgzfBlock`],
/// as needed to build a BAI index inline with compression.
#[derive(Debug, Clone, Copy)]
pub struct RecordIndexEntry {
    /// Offset of the record's 4-byte `block_size` length prefix, measured
    /// from the start of the *uncompressed* bytes of the batch that produced
    /// this block (i.e. the BGZF virtual offset's uncompressed-offset half,
    /// before the physical block offset is added in).
    pub uoffset: u32,
    /// Total record length in bytes: `4` (the length prefix) plus the body
    /// length.
    pub len: u32,
    /// The record's alignment context (reference id, start/end, mapped
    /// flag), or `None` for a record with no reference/position (fully
    /// unplaced unmapped reads are not indexable by coordinate).
    pub ctx: Option<fgumi_bam_io::AlignmentContext>,
}

/// Sidecar metadata riding alongside a compressed [`BgzfBlock`], carrying
/// everything an inline BAI indexer needs to attribute records to virtual
/// offsets without re-parsing the compressed bytes.
///
/// # Invariants
///
/// - `phys_comp_len.iter().map(|&n| n as u64).sum::<u64>() == bytes.len() as
///   u64`, where `bytes` is the sibling [`BgzfBlock::bytes`] this manifest
///   describes: the physical lengths of the constituent BGZF blocks sum to
///   the total compressed byte count.
/// - `phys_comp_len.len() == (uncompressed_size as
///   usize).div_ceil(fgumi_bgzf::BGZF_MAX_BLOCK_SIZE)`, where
///   `uncompressed_size` is the sibling [`BgzfBlock::uncompressed_size`]:
///   one physical block per `BGZF_MAX_BLOCK_SIZE`-sized (or smaller, for the
///   final one) chunk of uncompressed input.
/// - `records` is in batch order, which is final coordinate order.
/// - Each entry's `uoffset` is the offset of the record's 4-byte length
///   prefix within the batch's uncompressed bytes; `len` is `4 +` the body
///   length.
#[derive(Debug, Clone)]
pub struct BamIndexManifest {
    /// Physical (compressed) length of each constituent BGZF block, in the
    /// order those blocks appear in the sibling [`BgzfBlock::bytes`].
    pub phys_comp_len: Vec<u32>,
    /// Per-record virtual-offset bookkeeping, in batch (= final coordinate)
    /// order.
    pub records: Vec<RecordIndexEntry>,
}

// ─────────────────────────────────────────────────────────────────────────────
// BgzfBlock — raw compressed BGZF block + read-order serial.
// ─────────────────────────────────────────────────────────────────────────────

/// Raw compressed BGZF block + parsed metadata. Carries read-order serial.
///
/// Sentinel/EOF blocks have `bytes` containing the 28-byte BGZF EOF marker
/// and `uncompressed_size = 0`.
#[derive(Debug)]
pub struct BgzfBlock {
    /// Read-order serial. Set by `ReadBgzfBlocks` based on block read index.
    pub batch_serial: u64,
    pub bytes: Vec<u8>,
    /// Decompressed size, parsed from the BGZF block header.
    pub uncompressed_size: u32,
    /// Sidecar index metadata for this block's records, populated only on the
    /// inline-BAI-indexed compress path (`None` for every other producer,
    /// including sentinel/EOF blocks).
    pub index: Option<Box<BamIndexManifest>>,
}

impl HeapSize for BgzfBlock {
    fn heap_size(&self) -> usize {
        // Byte-bounded queues budget on resident heap, so account for the full
        // allocation (`capacity`), not just the populated prefix (`len`).
        let manifest = self.index.as_ref().map_or(0, |m| {
            std::mem::size_of::<BamIndexManifest>()
                + m.phys_comp_len.capacity() * std::mem::size_of::<u32>()
                + m.records.capacity() * std::mem::size_of::<RecordIndexEntry>()
        });
        self.bytes.capacity() + manifest
    }
}

impl Ordered for BgzfBlock {
    fn ordinal(&self) -> u64 {
        self.batch_serial
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// DecompressedBlock — raw bytes from a BGZF decompression, record-aligned
// or not.
// ─────────────────────────────────────────────────────────────────────────────

/// Decompressed bytes from one or more BGZF blocks. Carries a serial for
/// ordering purposes; record alignment is the consumer's responsibility.
#[derive(Debug)]
pub struct DecompressedBlock {
    pub batch_serial: u64,
    pub bytes: Vec<u8>,
}

impl HeapSize for DecompressedBlock {
    fn heap_size(&self) -> usize {
        // Account for the full allocation (`capacity`), not just `len` — see
        // the `BgzfBlock` impl above.
        self.bytes.capacity()
    }
}

impl Ordered for DecompressedBlock {
    fn ordinal(&self) -> u64 {
        self.batch_serial
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// RecordBatch — parsed BAM records grouped into a batch.
// ─────────────────────────────────────────────────────────────────────────────

/// A batch of parsed BAM records, stored as a flat backing buffer + per-record
/// `(start, end)` ranges.
#[derive(Debug)]
pub struct RecordBatch {
    batch_serial: u64,
    /// All record bodies concatenated, in batch order.
    backing: Vec<u8>,
    /// (start, end) byte ranges into `backing`, one per record.
    ranges: Vec<(u32, u32)>,
}

impl RecordBatch {
    /// Construct a batch from a pre-parsed backing buffer and `(start, end)` ranges.
    ///
    /// Each range must satisfy `start <= end <= backing.len()`; callers own that
    /// invariant (the boundary scans that produce these ranges already validate
    /// it). Violating it panics later in
    /// [`iter_record_bytes`](Self::iter_record_bytes) when the range is sliced,
    /// far from the site that produced it — the `debug_assert!` pins the failure
    /// at the constructor boundary instead. Note that `total_bytes` reports
    /// `backing.len()`, which differs from the sum of the record bodies if the
    /// ranges do not cover the whole buffer.
    #[must_use]
    pub fn from_parsed(batch_serial: u64, backing: Vec<u8>, ranges: Vec<(u32, u32)>) -> Self {
        debug_assert!(
            ranges.iter().all(|&(s, e)| s <= e && e as usize <= backing.len()),
            "RecordBatch::from_parsed: range outside backing buffer"
        );
        Self { batch_serial, backing, ranges }
    }

    /// Convenience constructor: serializes a slice of `RawRecord`s into the
    /// flat representation.
    ///
    /// # Panics
    ///
    /// Panics if the concatenated record bodies exceed `u32::MAX` bytes.
    #[must_use]
    pub fn new(batch_serial: u64, records: &[RawRecord]) -> Self {
        let total: usize = records.iter().map(RawRecord::len).sum();
        let mut backing = Vec::with_capacity(total);
        let mut ranges = Vec::with_capacity(records.len());
        for rec in records {
            let start = u32::try_from(backing.len()).expect("backing fits in u32");
            backing.extend_from_slice(rec.as_ref());
            let end = u32::try_from(backing.len()).expect("backing fits in u32");
            ranges.push((start, end));
        }
        Self { batch_serial, backing, ranges }
    }

    /// Self-managed ordinal.
    #[must_use]
    pub fn batch_serial(&self) -> u64 {
        self.batch_serial
    }

    /// Number of records in the batch.
    #[must_use]
    pub fn len(&self) -> usize {
        self.ranges.len()
    }

    /// `true` iff the batch contains zero records.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.ranges.is_empty()
    }

    /// Total bytes across all record bodies.
    #[must_use]
    pub fn total_bytes(&self) -> usize {
        self.backing.len()
    }

    /// Iterate the record bodies as borrowed byte slices into the backing buffer.
    pub fn iter_record_bytes(&self) -> impl Iterator<Item = &[u8]> + '_ {
        let backing = &self.backing[..];
        self.ranges.iter().map(move |&(s, e)| &backing[s as usize..e as usize])
    }
}

/// Builder for emit-side `RecordBatch` construction.
#[derive(Debug)]
pub struct RecordBatchBuilder {
    batch_serial: u64,
    backing: Vec<u8>,
    ranges: Vec<(u32, u32)>,
}

impl RecordBatchBuilder {
    /// Build an empty builder with reserved capacity.
    #[must_use]
    pub fn with_capacity(batch_serial: u64, bytes_cap: usize, records_cap: usize) -> Self {
        Self {
            batch_serial,
            backing: Vec::with_capacity(bytes_cap),
            ranges: Vec::with_capacity(records_cap),
        }
    }

    /// Append one record's body bytes.
    ///
    /// # Panics
    ///
    /// Panics if accumulated bytes would exceed `u32::MAX`.
    pub fn push_record_bytes(&mut self, bytes: &[u8]) {
        let start = u32::try_from(self.backing.len()).expect("backing fits in u32");
        self.backing.extend_from_slice(bytes);
        let end = u32::try_from(self.backing.len()).expect("backing fits in u32");
        self.ranges.push((start, end));
    }

    /// Number of records appended so far.
    #[must_use]
    pub fn len(&self) -> usize {
        self.ranges.len()
    }

    /// `true` iff no records have been appended.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.ranges.is_empty()
    }

    /// Total record bytes appended so far.
    #[must_use]
    pub fn total_bytes(&self) -> usize {
        self.backing.len()
    }

    /// Finalize and produce the `RecordBatch`. Consumes the builder.
    #[must_use]
    pub fn build(self) -> RecordBatch {
        RecordBatch { batch_serial: self.batch_serial, backing: self.backing, ranges: self.ranges }
    }
}

impl HeapSize for RecordBatch {
    fn heap_size(&self) -> usize {
        // Account for the full allocation (`capacity`) of both buffers, not
        // just their populated prefixes — see the `BgzfBlock` impl above.
        self.backing.capacity() + self.ranges.capacity() * std::mem::size_of::<(u32, u32)>()
    }
}

impl Ordered for RecordBatch {
    fn ordinal(&self) -> u64 {
        self.batch_serial
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bgzf_block_heap_size_matches_bytes_capacity() {
        let b = BgzfBlock {
            batch_serial: 0,
            bytes: vec![0u8; 1024],
            uncompressed_size: 4096,
            index: None,
        };
        assert_eq!(b.heap_size(), 1024);
        assert_eq!(b.ordinal(), 0);
    }

    #[test]
    fn bgzf_block_index_defaults_none_and_heap_counts_manifest() {
        let plain =
            BgzfBlock { batch_serial: 0, bytes: vec![0u8; 100], uncompressed_size: 0, index: None };
        assert!(plain.index.is_none());
        assert_eq!(plain.heap_size(), plain.bytes.capacity());

        let manifest = Box::new(BamIndexManifest {
            phys_comp_len: vec![10, 20],
            records: vec![RecordIndexEntry { uoffset: 0, len: 50, ctx: None }],
        });
        let indexed = BgzfBlock {
            batch_serial: 1,
            bytes: vec![0u8; 30],
            uncompressed_size: 0,
            index: Some(manifest),
        };
        let m = indexed.index.as_ref().expect("manifest present");
        assert_eq!(
            indexed.heap_size(),
            indexed.bytes.capacity()
                + std::mem::size_of::<BamIndexManifest>()
                + m.phys_comp_len.capacity() * std::mem::size_of::<u32>()
                + m.records.capacity() * std::mem::size_of::<RecordIndexEntry>(),
            "heap_size must account for the manifest allocation exactly"
        );
    }

    #[test]
    fn decompressed_block_heap_size_matches_bytes_capacity() {
        let b = DecompressedBlock { batch_serial: 7, bytes: vec![0u8; 4096] };
        assert_eq!(b.heap_size(), 4096);
        assert_eq!(b.ordinal(), 7);
    }

    #[test]
    fn heap_size_counts_allocated_capacity_not_logical_len() {
        // A buffer with spare capacity (e.g. after `with_capacity`) holds more
        // resident heap than its `len` — byte-bounded queues must budget the
        // full allocation or they undercount and bypass configured limits.
        let mut bytes = Vec::with_capacity(8192);
        bytes.extend_from_slice(&[0u8; 100]);
        assert!(bytes.capacity() >= 8192 && bytes.len() == 100);
        let cap = bytes.capacity();

        let block = BgzfBlock { batch_serial: 0, bytes, uncompressed_size: 0, index: None };
        assert_eq!(block.heap_size(), cap);

        let mut backing = Vec::with_capacity(4096);
        backing.extend_from_slice(&[0u8; 10]);
        let mut ranges = Vec::with_capacity(64);
        ranges.push((0u32, 10u32));
        let backing_cap = backing.capacity();
        let ranges_cap = ranges.capacity();
        let batch = RecordBatch::from_parsed(0, backing, ranges);
        assert_eq!(batch.heap_size(), backing_cap + ranges_cap * std::mem::size_of::<(u32, u32)>());
    }

    #[test]
    fn record_batch_total_bytes_sums_record_lengths() {
        let r1: RawRecord = vec![0u8; 100].into();
        let r2: RawRecord = vec![0u8; 200].into();
        let batch = RecordBatch::new(3, &[r1, r2]);
        assert_eq!(batch.len(), 2);
        assert_eq!(batch.total_bytes(), 300);
        // `Vec::with_capacity` may over-allocate, so assert against the actual
        // allocated capacities rather than the requested sizes (see the sibling
        // `heap_size_counts_allocated_capacity_not_logical_len` test).
        assert_eq!(
            batch.heap_size(),
            batch.backing.capacity() + batch.ranges.capacity() * std::mem::size_of::<(u32, u32)>()
        );
        assert_eq!(batch.ordinal(), 3);
    }

    #[test]
    fn record_batch_from_parsed_round_trips_ranges() {
        let backing = b"AAABBBBCC".to_vec();
        let ranges = vec![(0u32, 3u32), (3u32, 7u32), (7u32, 9u32)];
        let batch = RecordBatch::from_parsed(11, backing, ranges);
        let got: Vec<&[u8]> = batch.iter_record_bytes().collect();
        assert_eq!(got, vec![&b"AAA"[..], &b"BBBB"[..], &b"CC"[..]]);
        assert_eq!(batch.len(), 3);
        assert_eq!(batch.total_bytes(), 9);
        assert_eq!(batch.batch_serial(), 11);
    }

    #[test]
    fn record_batch_builder_collects_records() {
        let mut b = RecordBatchBuilder::with_capacity(0, 64, 4);
        b.push_record_bytes(&[0u8; 50]);
        b.push_record_bytes(&[0u8; 75]);
        assert_eq!(b.len(), 2);
        assert!(!b.is_empty());
        assert_eq!(b.total_bytes(), 125);
        let batch = b.build();
        let got: Vec<usize> = batch.iter_record_bytes().map(<[u8]>::len).collect();
        assert_eq!(got, vec![50, 75]);
        // `heap_size` budgets allocated capacity (not logical length): the
        // builder was seeded with a 64-byte backing buffer but 125 bytes were
        // pushed, so `backing` reallocated and its capacity now exceeds 125.
        assert_eq!(
            batch.heap_size(),
            batch.backing.capacity() + batch.ranges.capacity() * std::mem::size_of::<(u32, u32)>()
        );
        assert!(batch.heap_size() >= 125 + 2 * std::mem::size_of::<(u32, u32)>());
    }
}

//! Concrete data types that flow through the BAM step library.
//!
//! Every flowing type carries an explicit `batch_serial: u64` field and
//! impls both [`HeapSize`] (so byte-bounded queues can budget memory) and
//! [`Ordered`] (so `BranchOrdering::ByItemOrdinal` reorder stages preserve
//! global ordering across multi-step Parallel transforms).

use fgumi_pipeline_core::{HeapSize, Ordered};
use fgumi_raw_bam::RawRecord;

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
}

impl HeapSize for BgzfBlock {
    fn heap_size(&self) -> usize {
        // Byte-bounded queues budget on resident heap, so account for the full
        // allocation (`capacity`), not just the populated prefix (`len`).
        self.bytes.capacity()
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
    #[must_use]
    pub fn from_parsed(batch_serial: u64, backing: Vec<u8>, ranges: Vec<(u32, u32)>) -> Self {
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
        let b = BgzfBlock { batch_serial: 0, bytes: vec![0u8; 1024], uncompressed_size: 4096 };
        assert_eq!(b.heap_size(), 1024);
        assert_eq!(b.ordinal(), 0);
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

        let block = BgzfBlock { batch_serial: 0, bytes, uncompressed_size: 0 };
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

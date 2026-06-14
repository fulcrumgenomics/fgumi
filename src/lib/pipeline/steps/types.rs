//! Concrete data types that flow through the BAM step library.
//!
//! Every flowing type carries an explicit `batch_serial: u64` field and
//! impls both [`HeapSize`] (so byte-bounded queues can budget memory) and
//! [`Ordered`] (so `BranchOrdering::ByItemOrdinal` reorder stages preserve
//! global ordering across multi-step Parallel transforms).
//!
//! Serial propagation: every transform copies the input's `batch_serial`
//! onto its output items. For Serial steps that span batch boundaries
//! (`FindBamBoundaries`, `GroupBam`), the output's serial is the
//! serial of the *last* contributing input — this preserves monotonicity
//! since the output ordinal of batch N is always ≥ the output ordinal of
//! batch N-1.
//!
//! `BamTemplateBatch` is named to disambiguate from
//! `crate::template::TemplateBatch` (a legacy type alias for
//! `Vec<Template>`); the new framework's batch carries metadata alongside
//! the templates.

use fgumi_raw_bam::RawRecord;

use crate::pipeline::core::item::{HeapSize, Ordered};
use crate::template::Template;
use fgumi_bam_io::DecodedRecord;

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
        self.bytes.len()
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
/// ordering purposes; record alignment is the consumer's responsibility
/// (see `FindBamBoundaries`).
#[derive(Debug)]
pub struct DecompressedBlock {
    pub batch_serial: u64,
    pub bytes: Vec<u8>,
}

impl HeapSize for DecompressedBlock {
    fn heap_size(&self) -> usize {
        self.bytes.len()
    }
}

impl Ordered for DecompressedBlock {
    fn ordinal(&self) -> u64 {
        self.batch_serial
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// RecordBatch — parsed BAM records grouped into a batch.
//
// Internal representation: a single flat `backing` buffer holding all record
// bodies concatenated, plus per-record (start, end) byte ranges into it. This
// avoids the per-record `Vec<u8>` allocation that an
// earlier `Vec<RawRecord>` shape required during high-throughput parse and
// merge-emit paths (ParseBamRecords parses N records in a block; Sort's
// `merge_into` emits N records per batch — both can run zero-alloc per record
// at the inner loop, with one `Vec` allocation per batch instead).
//
// Public API:
//   - `iter_record_bytes()`           — borrowing iterator of `&[u8]` views
//   - `len()`, `is_empty()`           — record count
//   - `total_bytes()`, `batch_serial()` — accessors
//   - [`RecordBatchBuilder`]          — emit-side builder (`push_record_bytes`)
//   - `RecordBatch::from_parsed(...)` — ParseBam-side constructor
//   - `RecordBatch::new(...)`         — convenience for tests; serializes a
//                                       `Vec<RawRecord>` into the flat repr.
// ─────────────────────────────────────────────────────────────────────────────

/// A batch of parsed BAM records, stored as a flat backing buffer + per-record
/// `(start, end)` ranges. See module-level note for rationale.
#[derive(Debug)]
pub struct RecordBatch {
    batch_serial: u64,
    /// All record bodies concatenated, in batch order.
    backing: Vec<u8>,
    /// (start, end) byte ranges into `backing`, one per record. `u32` is
    /// sufficient since BAM record body size is u32-bounded and we never
    /// concat more than `u32::MAX` bytes per batch in practice
    /// (`tuning.per_step_byte_limit` defaults to 4 MiB).
    ranges: Vec<(u32, u32)>,
}

impl RecordBatch {
    /// Construct a batch from a pre-parsed backing buffer and `(start, end)`
    /// ranges. Used by `ParseBamRecords` to wrap a `DecompressedBlock`'s
    /// bytes without re-allocating per record.
    #[must_use]
    pub fn from_parsed(batch_serial: u64, backing: Vec<u8>, ranges: Vec<(u32, u32)>) -> Self {
        Self { batch_serial, backing, ranges }
    }

    /// Convenience constructor: serializes a slice of `RawRecord`s into the
    /// flat representation. Intended for tests; production emitters should
    /// use [`RecordBatchBuilder`] to avoid the intermediate `Vec<RawRecord>`
    /// allocation entirely.
    ///
    /// # Panics
    ///
    /// Panics if the concatenated record bodies exceed `u32::MAX` bytes
    /// (the `(start, end)` range type for each record is `(u32, u32)`).
    /// In practice this is unreachable for typical BAM workloads —
    /// `tuning.per_step_byte_limit` defaults to 4 MiB.
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

    /// Iterate the record bodies as borrowed byte slices into the backing
    /// buffer. Zero-allocation per record.
    pub fn iter_record_bytes(&self) -> impl Iterator<Item = &[u8]> + '_ {
        let backing = &self.backing[..];
        self.ranges.iter().map(move |&(s, e)| &backing[s as usize..e as usize])
    }
}

/// Builder for emit-side `RecordBatch` construction. Used by `Sort`'s
/// `merge_into` closure: append record bodies one at a time into a shared
/// backing buffer; finalize with [`Self::build`]. One backing allocation
/// per emitted batch instead of one per record.
#[derive(Debug)]
pub struct RecordBatchBuilder {
    batch_serial: u64,
    backing: Vec<u8>,
    ranges: Vec<(u32, u32)>,
}

impl RecordBatchBuilder {
    /// Build an empty builder with reserved capacity for `bytes_cap` bytes of
    /// record body and `records_cap` records.
    #[must_use]
    pub fn with_capacity(batch_serial: u64, bytes_cap: usize, records_cap: usize) -> Self {
        Self {
            batch_serial,
            backing: Vec::with_capacity(bytes_cap),
            ranges: Vec::with_capacity(records_cap),
        }
    }

    /// Append one record's body bytes. Memcpys into the backing buffer; no
    /// per-record allocation.
    ///
    /// # Panics
    ///
    /// Panics if accumulated bytes would exceed `u32::MAX` (each
    /// record's `(start, end)` range is `(u32, u32)`). Unreachable in
    /// practice — `tuning.per_step_byte_limit` defaults to 4 MiB.
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
        // Backing bytes plus the ranges Vec's heap. (Record count × 8 bytes
        // since each range is `(u32, u32)`.)
        self.backing.len() + self.ranges.len() * std::mem::size_of::<(u32, u32)>()
    }
}

impl Ordered for RecordBatch {
    fn ordinal(&self) -> u64 {
        self.batch_serial
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// DecodedRecordBatch — parsed records with pre-computed `GroupKey`s.
// Output of the parallel `DecodeRecords` step; input to the serial
// `GroupBam` step. Mirrors legacy's split where key computation runs in
// the parallel Decode step (`pipeline/bam.rs:329-403`) so that
// the serial Group step only does fast hash-and-name comparisons.
// ─────────────────────────────────────────────────────────────────────────────

/// A batch of `DecodedRecord`s — raw record bytes paired with a
/// pre-computed `GroupKey`. Carries `total_bytes` for O(1) `HeapSize`.
///
/// `total_bytes` is cached at construction (`new`) and assumes the
/// `records` vector is **not** mutated in a size-changing way afterwards.
/// The fields are `pub` for ergonomics, but in-place edits that change a
/// record's `raw_bytes().len()` would desync `heap_size()` from reality —
/// rebuild via [`Self::new`] instead of mutating in place.
#[derive(Debug)]
pub struct DecodedRecordBatch {
    pub batch_serial: u64,
    pub records: Vec<DecodedRecord>,
    pub total_bytes: usize,
}

impl DecodedRecordBatch {
    #[must_use]
    pub fn new(batch_serial: u64, records: Vec<DecodedRecord>) -> Self {
        let total_bytes = records.iter().map(|d| d.raw_bytes().len()).sum();
        Self { batch_serial, records, total_bytes }
    }
}

impl HeapSize for DecodedRecordBatch {
    fn heap_size(&self) -> usize {
        self.total_bytes
    }
}

impl Ordered for DecodedRecordBatch {
    fn ordinal(&self) -> u64 {
        self.batch_serial
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BamTemplateBatch — grouped templates carrying batch metadata.
// Disambiguated from `crate::template::TemplateBatch` (= `Vec<Template>`).
// ─────────────────────────────────────────────────────────────────────────────

/// A batch of templates (same-queryname record groups), with batch serial
/// for ordering and `total_bytes` for memory accounting.
#[derive(Debug)]
pub struct BamTemplateBatch {
    pub batch_serial: u64,
    pub templates: Vec<Template>,
    pub total_bytes: usize,
}

impl BamTemplateBatch {
    #[must_use]
    pub fn new(batch_serial: u64, templates: Vec<Template>) -> Self {
        let total_bytes = templates.iter().map(Template::heap_size).sum();
        Self { batch_serial, templates, total_bytes }
    }

    /// Construct from a pre-computed `total_bytes`, skipping the
    /// per-template `heap_size` walk inside `new`. Use when the
    /// caller already summed `heap_size` in a prior pass (e.g.,
    /// the AAM merge loop, which folds heap accounting into the
    /// same loop that builds the `Vec<Template>`).
    #[must_use]
    pub fn from_parts(batch_serial: u64, templates: Vec<Template>, total_bytes: usize) -> Self {
        Self { batch_serial, templates, total_bytes }
    }
}

impl HeapSize for BamTemplateBatch {
    fn heap_size(&self) -> usize {
        self.total_bytes
    }
}

impl Ordered for BamTemplateBatch {
    fn ordinal(&self) -> u64 {
        self.batch_serial
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SamChunk — a slab of SAM text bytes ending on a newline boundary, plus an
// inline line-offset table. Emitted by `ReadSamChunks` (Serial+Reader) and
// consumed by `ParseSamChunk` (Parallel). The parallel parser slices each
// line via `bytes[offsets[i]..offsets[i+1]]` — no cross-chunk dependencies
// because the Read step never splits a record across two chunks.
// ─────────────────────────────────────────────────────────────────────────────

/// A chunk of SAM text plus its sentinel-form line-offset table.
///
/// `line_offsets` has `N+1` entries describing `N` complete records. Line
/// `i` is `bytes[line_offsets[i] as usize..line_offsets[i+1] as usize]`
/// and includes its trailing `\n`. The chunk always ends on a newline
/// boundary (the Read step holds any trailing partial line as carryover).
#[derive(Debug)]
pub struct SamChunk {
    pub batch_serial: u64,
    pub bytes: Vec<u8>,
    pub line_offsets: Vec<u32>,
}

impl SamChunk {
    /// Number of complete records carried by this chunk.
    #[must_use]
    pub fn record_count(&self) -> usize {
        self.line_offsets.len().saturating_sub(1)
    }
}

impl HeapSize for SamChunk {
    fn heap_size(&self) -> usize {
        self.bytes.len() + self.line_offsets.len() * std::mem::size_of::<u32>()
    }
}

impl Ordered for SamChunk {
    fn ordinal(&self) -> u64 {
        self.batch_serial
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bgzf_block_heap_size_matches_bytes_len() {
        let b = BgzfBlock { batch_serial: 0, bytes: vec![0u8; 1024], uncompressed_size: 4096 };
        assert_eq!(b.heap_size(), 1024);
        assert_eq!(b.ordinal(), 0);
    }

    #[test]
    fn decompressed_block_heap_size_matches_bytes_len() {
        let b = DecompressedBlock { batch_serial: 7, bytes: vec![0u8; 4096] };
        assert_eq!(b.heap_size(), 4096);
        assert_eq!(b.ordinal(), 7);
    }

    #[test]
    fn record_batch_total_bytes_sums_record_lengths() {
        // For the size-accounting test we just need records with known
        // byte lengths; the bytes need not be valid BAM records.
        let r1: RawRecord = vec![0u8; 100].into();
        let r2: RawRecord = vec![0u8; 200].into();
        let batch = RecordBatch::new(3, &[r1, r2]);
        assert_eq!(batch.len(), 2);
        assert_eq!(batch.total_bytes(), 300);
        // `heap_size` includes the ranges Vec (8 bytes per entry).
        assert_eq!(batch.heap_size(), 300 + 2 * std::mem::size_of::<(u32, u32)>());
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
        assert_eq!(batch.heap_size(), 125 + 2 * std::mem::size_of::<(u32, u32)>());
    }

    #[test]
    fn template_batch_carries_serial() {
        let batch = BamTemplateBatch::new(42, vec![]);
        assert_eq!(batch.heap_size(), 0);
        assert_eq!(batch.ordinal(), 42);
    }

    #[test]
    fn template_batch_new_total_bytes_sums_heap_sizes() {
        // Template::new(name) → heap_size() == name.len() (no records).
        // Two templates with distinct known names; no RawRecord builders needed.
        let t1 = crate::template::Template::new(b"read1".to_vec()); // 5 bytes
        let t2 = crate::template::Template::new(b"read_two".to_vec()); // 8 bytes
        let expected = t1.heap_size() + t2.heap_size(); // 13
        let batch = BamTemplateBatch::new(7, vec![t1, t2]);
        assert_eq!(batch.total_bytes, expected);
        assert_eq!(batch.heap_size(), expected); // HeapSize delegates to total_bytes
        assert_eq!(batch.ordinal(), 7);
    }

    #[test]
    fn template_batch_from_parts_trusts_caller_total_bytes() {
        // from_parts must NOT recompute total_bytes — it trusts the caller
        // (e.g. AAM merge loop that folds heap accounting into its own loop).
        let batch = BamTemplateBatch::from_parts(99, vec![], 12345);
        assert_eq!(batch.heap_size(), 12345);
        assert_eq!(batch.ordinal(), 99);
    }

    #[test]
    fn sam_chunk_record_count_is_offsets_minus_one() {
        // Sentinel-form line_offsets: N+1 entries describe N records.
        // `bytes` and `line_offsets` come from `split_complete_lines`.
        let bytes = b"abc\nde\n".to_vec();
        let line_offsets = vec![0u32, 4, 7];
        let chunk = SamChunk { batch_serial: 5, bytes, line_offsets };
        assert_eq!(chunk.record_count(), 2);
    }

    #[test]
    fn sam_chunk_record_count_is_zero_for_empty_offsets() {
        let chunk = SamChunk { batch_serial: 0, bytes: Vec::new(), line_offsets: Vec::new() };
        assert_eq!(chunk.record_count(), 0);
    }

    #[test]
    fn sam_chunk_heap_size_sums_bytes_and_offsets() {
        let chunk =
            SamChunk { batch_serial: 1, bytes: vec![0u8; 500], line_offsets: vec![0u32; 10] };
        assert_eq!(chunk.heap_size(), 500 + 10 * std::mem::size_of::<u32>());
        assert_eq!(chunk.ordinal(), 1);
    }
}

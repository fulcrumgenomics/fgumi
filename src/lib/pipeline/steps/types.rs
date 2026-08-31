//! Concrete data types that flow through the BAM step library.
//!
//! Every flowing type carries an explicit `batch_serial: u64` field and
//! impls both [`HeapSize`] (so byte-bounded queues can budget memory) and
//! [`Ordered`] (so `BranchOrdering::ByItemOrdinal` reorder stages preserve
//! global ordering across multi-step Parallel transforms).
//!
//! Serial propagation: every transform copies the input's `batch_serial`
//! onto its output items when the mapping is 1:1. Steps that aggregate or
//! split across input batches (`GroupBam`, `GroupByPosition`,
//! `GroupByQueryname`) instead mint their own monotonic counter from 0, one
//! value per emitted batch: `ByItemOrdinal` requires a *contiguous* sequence,
//! and inheriting an input serial would skip ordinals and stall the
//! downstream `ReorderStage`.
//!
//! `BamTemplateBatch` is named to disambiguate from
//! `crate::template::TemplateBatch` (a legacy type alias for
//! `Vec<Template>`); the new framework's batch carries metadata alongside
//! the templates.

use crate::pipeline::core::item::{HeapSize, Ordered};
use crate::template::Template;
use fgumi_bam_io::{DecodedRecord, MemoryEstimate};

// ─────────────────────────────────────────────────────────────────────────────
// Light flowing types (BgzfBlock, DecompressedBlock, RecordBatch,
// RecordBatchBuilder) moved to the `fgumi-pipeline-io` crate. Re-export them so
// `crate::pipeline::steps::types::*` paths keep resolving and there is ONE
// canonical definition of each. The heavy types below (BamTemplateBatch,
// DecodedRecordBatch, SamChunk) stay in the umbrella because they depend on
// `crate::template::Template`.
// ─────────────────────────────────────────────────────────────────────────────

pub use fgumi_pipeline_io::types::{
    BamIndexManifest, BgzfBlock, DecompressedBlock, RecordBatch, RecordBatchBuilder,
    RecordIndexEntry,
};

/// Bytes held by a `Vec`'s backing store, counted by **allocated capacity**.
///
/// The batch types below cache a `total_bytes` that feeds byte-bounded
/// admission control, so it must describe the memory the queue is really
/// holding. Summing only the elements' own heap misses the container itself: a
/// 1,000-element `Vec<DecodedRecord>` is tens of KiB before a single record
/// body is counted, and a `Vec` grown by repeated `push` typically carries up
/// to 2x its logical length in reserved slots. Counting `capacity()` rather
/// than `len()` is the same rule the per-record accounting already follows.
///
/// Takes the `Vec`'s `capacity()` rather than a borrow of the `Vec` itself, so
/// the reserved-slot count is explicit at every call site and the helper cannot
/// silently degrade to `len()` through a slice deref.
fn container_bytes<T>(capacity: usize) -> usize {
    capacity * std::mem::size_of::<T>()
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
/// `total_bytes` is cached at construction (`new`) and sums each record's
/// **allocated capacity** ([`MemoryEstimate::estimate_heap_size`]), not its
/// logical `raw_bytes().len()`, plus the `Vec`'s own backing store (also by
/// allocated capacity). This feeds the pipeline's byte-backpressure budget,
/// which must track the true heap footprint: accounting by `len` under-counts
/// an over-allocated record, and ignoring the container under-counts the
/// reserved slots — either lets the pipeline buffer far past its budget.
///
/// The fields are **private** so the cached `total_bytes` can never drift out
/// of sync with `records`: the only constructor is [`Self::new`] (which
/// recomputes the cache), reads go through [`Self::records`] /
/// [`Self::total_bytes`], and taking ownership of the records consumes the
/// batch via [`Self::into_records`]. There is no way to mutate `records` in
/// place, so byte-bounded admission control always sees the true footprint.
#[derive(Debug)]
pub struct DecodedRecordBatch {
    batch_serial: u64,
    records: Vec<DecodedRecord>,
    total_bytes: usize,
}

impl DecodedRecordBatch {
    #[must_use]
    pub fn new(batch_serial: u64, records: Vec<DecodedRecord>) -> Self {
        let total_bytes = records.iter().map(MemoryEstimate::estimate_heap_size).sum::<usize>()
            + container_bytes::<DecodedRecord>(records.capacity());
        Self { batch_serial, records, total_bytes }
    }

    /// The batch's ordering serial.
    #[must_use]
    pub fn batch_serial(&self) -> u64 {
        self.batch_serial
    }

    /// Immutable view of the batch's decoded records.
    #[must_use]
    pub fn records(&self) -> &[DecodedRecord] {
        &self.records
    }

    /// Cached heap footprint: each record's allocated capacity plus the
    /// `Vec`'s own backing store.
    #[must_use]
    pub fn total_bytes(&self) -> usize {
        self.total_bytes
    }

    /// Consume the batch and return ownership of its records. The cached
    /// `total_bytes` is dropped with the batch, so callers move the records
    /// out without any risk of leaving a stale byte count behind.
    #[must_use]
    pub fn into_records(self) -> Vec<DecodedRecord> {
        self.records
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
///
/// The fields are **private** so the cached `total_bytes` can never drift out
/// of sync with `templates`: construct via [`Self::new`] (which recomputes the
/// cache) or [`Self::from_parts`] (which trusts a caller that already summed
/// `heap_size`), read through the accessors, and take ownership of the
/// templates by consuming the batch via [`Self::into_parts`]. There is no way
/// to mutate `templates` in place, so byte-bounded admission control always
/// sees the true footprint.
#[derive(Debug)]
pub struct BamTemplateBatch {
    batch_serial: u64,
    templates: Vec<Template>,
    total_bytes: usize,
}

impl BamTemplateBatch {
    #[must_use]
    pub fn new(batch_serial: u64, templates: Vec<Template>) -> Self {
        let total_bytes = templates.iter().map(Template::heap_size).sum::<usize>()
            + container_bytes::<Template>(templates.capacity());
        Self { batch_serial, templates, total_bytes }
    }

    /// Construct from a pre-computed `total_bytes`, skipping the
    /// per-template `heap_size` walk inside `new`. Use when the
    /// caller already summed `heap_size` in a prior pass (e.g.,
    /// the AAM merge loop, which folds heap accounting into the
    /// same loop that builds the `Vec<Template>`).
    #[must_use]
    pub fn from_parts(batch_serial: u64, templates: Vec<Template>, total_bytes: usize) -> Self {
        // The caller supplies the summed per-template heap; the container's
        // reserved capacity is added here because the caller cannot know it.
        let total_bytes = total_bytes + container_bytes::<Template>(templates.capacity());
        Self { batch_serial, templates, total_bytes }
    }

    /// The batch's ordering serial.
    #[must_use]
    pub fn batch_serial(&self) -> u64 {
        self.batch_serial
    }

    /// Immutable view of the batch's templates.
    #[must_use]
    pub fn templates(&self) -> &[Template] {
        &self.templates
    }

    /// Cached heap footprint: summed `Template::heap_size` plus the `Vec`'s
    /// own backing store.
    #[must_use]
    pub fn total_bytes(&self) -> usize {
        self.total_bytes
    }

    /// Consume the batch and return its serial and templates. The cached
    /// `total_bytes` is dropped with the batch, so callers move the templates
    /// out (and re-wrap via [`Self::new`]) without risking a stale byte count.
    #[must_use]
    pub fn into_parts(self) -> (u64, Vec<Template>) {
        (self.batch_serial, self.templates)
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
        self.bytes.capacity() + container_bytes::<u32>(self.line_offsets.capacity())
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

    // Tests for the light types (BgzfBlock, DecompressedBlock, RecordBatch,
    // RecordBatchBuilder) moved with those types to the `fgumi-pipeline-io`
    // crate. Only the heavy-type tests remain below.

    #[test]
    fn template_batch_carries_serial() {
        let batch = BamTemplateBatch::new(42, vec![]);
        assert_eq!(batch.heap_size(), 0);
        assert_eq!(batch.ordinal(), 42);
    }

    #[test]
    fn decoded_record_batch_heap_size_counts_capacity_not_len() {
        // `heap_size` feeds the pipeline's byte-backpressure budget, so it must
        // report the records' true heap footprint (allocated capacity), not
        // their logical length. Accounting by `len` under-counts an
        // over-allocated record — exactly the failure mode where a SAM-parsed
        // record with a 4 KiB buffer holding ~60 bytes fooled backpressure and
        // let the pipeline buffer far past its memory budget.
        use fgumi_bam_io::{GroupKey, MemoryEstimate};
        let mut raw = Vec::with_capacity(4096);
        raw.extend_from_slice(&[0u8; 64]); // len 64, capacity 4096
        let rec = DecodedRecord::from_raw_bytes(raw, GroupKey::default());
        let cap = rec.estimate_heap_size();
        assert!(cap >= 4096, "sanity: over-allocated record capacity, got {cap}");
        let batch = DecodedRecordBatch::new(0, vec![rec]);
        let expected = cap + std::mem::size_of::<DecodedRecord>(); // one container slot
        assert_eq!(
            batch.heap_size(),
            expected,
            "heap_size must reflect allocated capacity ({cap}) plus the container, \
             not logical len — else byte-backpressure under-counts real memory",
        );
    }

    #[test]
    fn decoded_record_batch_accessors_and_into_records_round_trip() {
        // The fields are private; callers read via the accessors and take
        // ownership by consuming the batch. `into_records` must hand back the
        // exact records `new` was given, and `batch_serial` must round-trip.
        use fgumi_bam_io::GroupKey;
        let rec = DecodedRecord::from_raw_bytes(vec![7u8; 32], GroupKey::default());
        let batch = DecodedRecordBatch::new(3, vec![rec]);
        assert_eq!(batch.batch_serial(), 3);
        assert_eq!(batch.records().len(), 1);
        assert_eq!(batch.total_bytes(), batch.heap_size());
        let records = batch.into_records();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].raw_bytes().len(), 32);
    }

    #[test]
    fn bam_template_batch_accessors_and_into_parts_round_trip() {
        // Mirror of the DecodedRecordBatch contract: private fields, read via
        // accessors, consume via `into_parts` to move the templates out and
        // re-wrap through `new` (which recomputes the byte cache).
        let t1 = crate::template::Template::new(b"read1".to_vec());
        let expected_bytes = t1.heap_size() + std::mem::size_of::<crate::template::Template>();
        let batch = BamTemplateBatch::new(9, vec![t1]);
        assert_eq!(batch.batch_serial(), 9);
        assert_eq!(batch.templates().len(), 1);
        assert_eq!(batch.total_bytes(), expected_bytes);
        let (serial, templates) = batch.into_parts();
        assert_eq!(serial, 9);
        // Re-wrapping the moved-out templates recomputes an identical cache.
        let rewrapped = BamTemplateBatch::new(serial, templates);
        assert_eq!(rewrapped.total_bytes(), expected_bytes);
    }

    #[test]
    fn template_batch_new_total_bytes_sums_heap_sizes() {
        // Template::new(name) → heap_size() == name.len() (no records).
        // Two templates with distinct known names; no RawRecord builders needed.
        let t1 = crate::template::Template::new(b"read1".to_vec()); // 5 bytes
        let t2 = crate::template::Template::new(b"read_two".to_vec()); // 8 bytes
        let expected =
            t1.heap_size() + t2.heap_size() + 2 * std::mem::size_of::<crate::template::Template>(); // 13 + container
        let batch = BamTemplateBatch::new(7, vec![t1, t2]);
        assert_eq!(batch.total_bytes, expected);
        assert_eq!(batch.heap_size(), expected); // HeapSize delegates to total_bytes
        assert_eq!(batch.ordinal(), 7);
    }

    #[test]
    fn template_batch_from_parts_trusts_caller_total_bytes() {
        // from_parts must NOT re-walk the templates — it trusts the caller's
        // summed `heap_size` (e.g. the AAM merge loop, which folds heap
        // accounting into its own loop). It still adds the container term,
        // which the caller has no way to know.
        let batch = BamTemplateBatch::from_parts(99, vec![], 12345);
        assert_eq!(batch.heap_size(), 12345, "an empty Vec has no container allocation");
        assert_eq!(batch.ordinal(), 99);
    }

    // ========================================================================
    // Container capacity in the byte budget
    //
    // `heap_size` feeds byte-bounded admission control, so it must report the
    // batch's TRUE footprint. Each of these types owns a `Vec` whose backing
    // store is real memory the queue is holding: a 1000-element
    // `Vec<DecodedRecord>` is tens of KiB of container before a single record
    // body is counted. Omitting it lets a full queue exceed its configured
    // limit — the same class of under-count the per-record capacity accounting
    // above already guards against.
    // ========================================================================

    /// Reserved-but-unused slots are allocated memory, so the batch's byte
    /// count must grow with `capacity()`, not with `len()`.
    #[test]
    fn decoded_record_batch_counts_reserved_container_capacity() {
        use fgumi_bam_io::GroupKey;
        let record = || DecodedRecord::from_raw_bytes(vec![7u8; 32], GroupKey::default());

        let tight = DecodedRecordBatch::new(0, vec![record()]);

        let mut padded = Vec::with_capacity(256);
        padded.push(record());
        let padded = DecodedRecordBatch::new(0, padded);

        assert!(
            padded.heap_size() > tight.heap_size(),
            "a batch holding 256 reserved slots must out-weigh a tight one \
             ({} vs {})",
            padded.heap_size(),
            tight.heap_size(),
        );
        assert_eq!(
            padded.heap_size() - tight.heap_size(),
            (256 - 1) * std::mem::size_of::<DecodedRecord>(),
            "the difference must be exactly the extra reserved slots",
        );
    }

    /// Same contract for templates, through both constructors.
    #[test]
    fn template_batch_counts_reserved_container_capacity() {
        let template = || crate::template::Template::new(b"read1".to_vec());

        let mut padded = Vec::with_capacity(64);
        padded.push(template());
        let heap_only: usize = padded.iter().map(crate::template::Template::heap_size).sum();

        let via_new = BamTemplateBatch::new(0, padded);
        assert_eq!(
            via_new.total_bytes(),
            heap_only + 64 * std::mem::size_of::<crate::template::Template>(),
            "new() must add the container's reserved capacity to the summed heap",
        );

        // `from_parts` receives the caller's summed heap and must add the same
        // container term — otherwise the two constructors disagree on the
        // footprint of identical batches.
        let (_, templates) = via_new.into_parts();
        let via_parts = BamTemplateBatch::from_parts(0, templates, heap_only);
        assert_eq!(
            via_parts.total_bytes(),
            via_new_total_bytes_for_comparison(heap_only),
            "from_parts must agree with new() on the same templates",
        );
    }

    /// The expected `total_bytes` for the 64-slot fixture above, expressed once
    /// so both constructors are compared against the same definition.
    fn via_new_total_bytes_for_comparison(heap_only: usize) -> usize {
        heap_only + 64 * std::mem::size_of::<crate::template::Template>()
    }

    /// A `SamChunk`'s `bytes` come from a reader buffer that routinely carries
    /// spare capacity, so accounting by `len()` under-counts the slab the queue
    /// is actually holding.
    #[test]
    fn sam_chunk_heap_size_counts_capacity_not_len() {
        let mut bytes = Vec::with_capacity(8192);
        bytes.extend_from_slice(&[0u8; 500]);
        let mut line_offsets = Vec::with_capacity(64);
        line_offsets.extend_from_slice(&[0u32; 10]);

        let chunk = SamChunk { batch_serial: 1, bytes, line_offsets };
        assert_eq!(
            chunk.heap_size(),
            8192 + 64 * std::mem::size_of::<u32>(),
            "heap_size must report allocated capacity, not logical length",
        );
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

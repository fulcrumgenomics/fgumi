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

use crate::pipeline::core::item::{HeapSize, Ordered};
use crate::template::Template;
use fgumi_bam_io::DecodedRecord;

// ─────────────────────────────────────────────────────────────────────────────
// Light flowing types (BgzfBlock, DecompressedBlock, RecordBatch,
// RecordBatchBuilder) moved to the `fgumi-pipeline-io` crate. Re-export them so
// `crate::pipeline::steps::types::*` paths keep resolving and there is ONE
// canonical definition of each. The heavy types below (BamTemplateBatch,
// DecodedRecordBatch, SamChunk) stay in the umbrella because they depend on
// `crate::template::Template`.
// ─────────────────────────────────────────────────────────────────────────────

pub use fgumi_pipeline_io::types::{BgzfBlock, DecompressedBlock, RecordBatch, RecordBatchBuilder};

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

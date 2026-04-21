//! Item types for the decomposed BAM output chain.
//!
//! The extract (and future processing-path) pipelines serialize BAM records
//! in one stage, compress them in a parallel stage, and write them to a file
//! in a sequential sink. Each of those boundaries has an item type:
//!
//! - [`SerializedBatch`] — output of the serialize boundary (`ExtractStage` in
//!   plan 2c-3; `BamSerialize` in a future plan).
//! - [`CompressedBatch`] — output of the compress boundary (`BgzfCompress`).
//!
//! Both carry a monotonic `ordinal` stamped at the last Sequential coordination
//! stage upstream (`FastqPair` / `FastqBlockMerge` today). The ordinal is
//! preserved unchanged through the Parallel section and consumed by the
//! Sequential [`crate::runall::engine::sink::bam_file_write::BamFileWrite`] sink's
//! [`crate::runall::engine::reorder::ReorderBuffer`] to produce deterministic
//! output order.
//!
//! # Naming
//!
//! [`RawBytes`] and [`CompressedBytes`] are named for the state of their
//! `data` field. The "these are BAM bytes" context comes from the module
//! path. Keeping them as distinct types (rather than one type plus a phantom
//! marker) gives type-level certainty that a function expecting compressed
//! bytes cannot silently receive raw bytes.

/// Uncompressed, length-prefixed BAM record bytes for one output stream.
///
/// `data` is the concatenation of zero or more BAM records, each preceded by
/// its 4-byte little-endian `block_size`. The bytes are NOT BGZF-compressed.
#[derive(Debug, Clone)]
pub struct RawBytes {
    /// Length-prefixed BAM records, concatenated. Not BGZF-compressed.
    pub data: Vec<u8>,
    /// Number of BAM records encoded in `data`.
    pub record_count: u64,
}

/// BGZF-compressed bytes for one output stream.
///
/// `data` is a contiguous concatenation of complete BGZF blocks (each block's
/// header + compressed payload + footer). Ready to append to an output file
/// that has already had its BAM header written. `record_count` is preserved
/// from the upstream [`RawBytes`] — `data.len()` is byte count, not records.
#[derive(Debug, Clone)]
pub struct CompressedBytes {
    /// Concatenated BGZF blocks. Ready to append to an output file.
    pub data: Vec<u8>,
    /// Number of BAM records encoded (before compression).
    pub record_count: u64,
}

/// Sentinel ordinal used to mark batches that should be skipped by downstream
/// sinks (e.g. empty `TemplateBatch` items that carry no records and must not
/// enter the `ReorderBuffer`).
/// Output of the serialize boundary (`ExtractStage` in plan 2c-3).
///
/// `secondary` is `Some` only for dual-output pipelines (e.g. filter with a
/// rejects BAM). For single-output pipelines (extract), it is always `None`.
#[derive(Debug, Clone)]
pub struct SerializedBatch {
    pub primary: RawBytes,
    pub secondary: Option<RawBytes>,
    /// Monotonic batch ordinal stamped at the last Sequential coordination
    /// stage upstream. Preserved unchanged through `BgzfCompress`; consumed by
    /// `BamFileWrite`'s `ReorderBuffer`.
    pub ordinal: u64,
}

/// Output of the compress boundary (`BgzfCompress`).
#[derive(Debug, Clone)]
pub struct CompressedBatch {
    pub primary: CompressedBytes,
    pub secondary: Option<CompressedBytes>,
    pub ordinal: u64,
}

/// Plain (uncompressed) FASTQ bytes for one output batch.
///
/// Emitted by `ToFastq` (plan 2c-5a); consumed by `FastqFileSink` for
/// plain-FASTQ file output. A future `SubprocessStdinSink` (plan 2c-5b)
/// will consume the same type to pipe into an aligner's stdin.
///
/// `data` is zero or more FASTQ records back-to-back. Paired-end output
/// is interleaved: R1 record, then R2 record, then the next pair, etc.
/// Every record ends with a `\n` (LF), including the last — the buffer
/// is written directly to the output file or pipe.
#[derive(Debug, Clone)]
pub struct SerializedFastqBatch {
    /// FASTQ records (interleaved for paired-end).
    pub data: Vec<u8>,
    /// Number of FASTQ records in `data` (counts individual reads, not
    /// pairs — a paired batch with 10 pairs has `record_count == 20`).
    pub record_count: u64,
    /// Monotonic batch ordinal carried from the upstream `SerializedBatch`
    /// unchanged so the sink reorders correctly.
    pub ordinal: u64,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_raw_bytes_clone_preserves_fields() {
        let rb = RawBytes { data: vec![1, 2, 3, 4], record_count: 1 };
        let rb2 = rb.clone();
        assert_eq!(rb.data, rb2.data);
        assert_eq!(rb.record_count, rb2.record_count);
    }

    #[test]
    fn test_compressed_bytes_clone_preserves_fields() {
        let cb = CompressedBytes { data: vec![0x1f, 0x8b, 0x08], record_count: 7 };
        let cb2 = cb.clone();
        assert_eq!(cb.data, cb2.data);
        assert_eq!(cb.record_count, cb2.record_count);
    }

    #[test]
    fn test_serialized_batch_single_output_has_no_secondary() {
        let sb = SerializedBatch {
            primary: RawBytes { data: vec![42; 100], record_count: 3 },
            secondary: None,
            ordinal: 7,
        };
        assert!(sb.secondary.is_none());
        assert_eq!(sb.ordinal, 7);
        assert_eq!(sb.primary.record_count, 3);
    }

    #[test]
    fn test_serialized_batch_dual_output_preserves_both() {
        let sb = SerializedBatch {
            primary: RawBytes { data: vec![1; 10], record_count: 2 },
            secondary: Some(RawBytes { data: vec![2; 20], record_count: 1 }),
            ordinal: 99,
        };
        assert_eq!(sb.primary.record_count, 2);
        assert_eq!(sb.secondary.as_ref().unwrap().record_count, 1);
        assert_eq!(sb.ordinal, 99);
    }

    #[test]
    fn test_compressed_batch_clone_preserves_ordinal() {
        let cb = CompressedBatch {
            primary: CompressedBytes { data: vec![0; 50], record_count: 5 },
            secondary: None,
            ordinal: 123,
        };
        let cb2 = cb.clone();
        assert_eq!(cb.ordinal, cb2.ordinal);
        assert_eq!(cb.primary.record_count, cb2.primary.record_count);
    }

    #[test]
    fn test_serialized_fastq_batch_clone_preserves_fields() {
        let fb = SerializedFastqBatch {
            data: b"@r1\nACGT\n+\nIIII\n".to_vec(),
            record_count: 1,
            ordinal: 42,
        };
        let fb2 = fb.clone();
        assert_eq!(fb.data, fb2.data);
        assert_eq!(fb.record_count, fb2.record_count);
        assert_eq!(fb.ordinal, fb2.ordinal);
    }

    #[test]
    fn test_serialized_fastq_batch_empty() {
        let fb = SerializedFastqBatch { data: Vec::new(), record_count: 0, ordinal: 0 };
        assert!(fb.data.is_empty());
        assert_eq!(fb.record_count, 0);
    }
}

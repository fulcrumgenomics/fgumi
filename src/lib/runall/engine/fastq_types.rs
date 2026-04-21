//! Shared item types flowing through the FASTQ input pipeline stages.
//!
//! These mirror `PerStreamChunk` and `BlockParsed` from the unified pipeline
//! on main (`src/lib/unified_pipeline/fastq.rs`).

/// A per-stream chunk emitted by `FastqFileRead`.
///
/// For BGZF input: `data` is raw undecompressed BGZF blocks, `offsets` is None.
/// For gzip/plain input: `data` is decompressed record-aligned bytes, `offsets`
/// is Some with per-record byte offsets (the last offset == `data.len()`).
#[derive(Debug)]
pub struct PerStreamChunk {
    /// Which input stream this chunk came from (0=R1, 1=R2, etc.).
    pub stream_idx: usize,
    /// Monotonically increasing batch number per stream. Used by pair/merge
    /// stages to synchronize across streams.
    pub batch_num: u64,
    /// Raw BGZF or decompressed record-aligned bytes.
    pub data: Vec<u8>,
    /// Record-boundary offsets for gzip path (None for BGZF).
    pub offsets: Option<Vec<usize>>,
    /// Whether this chunk's stream has reached EOF after this chunk.
    pub is_last: bool,
}

/// Result of parallel per-chunk FASTQ parsing in the BGZF path.
///
/// `prefix_bytes` and `suffix_bytes` carry incomplete records that span
/// block boundaries; downstream `FastqBlockMerge` stitches them together.
#[derive(Debug)]
pub struct BlockParsed {
    /// Block index = per-stream batch number. Used for ordering at merge.
    pub block_idx: u64,
    /// Which input stream this block came from.
    pub stream_idx: usize,
    /// Fully-parsed complete records (zero-copy or owned).
    pub records: Vec<OwnedFastqRecord>,
    /// Leading incomplete record fragment (belongs to previous block).
    pub prefix_bytes: Vec<u8>,
    /// Trailing incomplete record fragment (belongs to next block).
    pub suffix_bytes: Vec<u8>,
    /// Whether this is the last block for this stream.
    pub is_last: bool,
}

/// Result of parallel parsing in the gzip path (no prefix/suffix because
/// `FastqFileRead` returns record-aligned data).
#[derive(Debug)]
pub struct FastqParsedStream {
    pub batch_num: u64,
    pub stream_idx: usize,
    pub records: Vec<OwnedFastqRecord>,
    pub is_last: bool,
}

/// Owned FASTQ record. Parser wraps `fgumi_simd_fastq::FastqRecord<'a>` and
/// owns the bytes so it can cross queue boundaries.
#[derive(Debug, Clone)]
pub struct OwnedFastqRecord {
    pub name: Vec<u8>,
    pub sequence: Vec<u8>,
    pub quality: Vec<u8>,
}

/// Multi-stream bundle: one record per input stream, zipped by position.
/// Mirrors `crate::grouper::FastqTemplate`.
#[derive(Debug)]
pub struct FastqTemplateV2 {
    /// One record per input stream, ordered by `stream_idx`.
    pub records: Vec<OwnedFastqRecord>,
    /// Common read name (shared across all streams).
    pub name: Vec<u8>,
}

/// A batch of [`FastqTemplateV2`] templates with a monotonic ordinal stamped
/// by the last Sequential coordination stage upstream (`FastqPair` /
/// `FastqBlockMerge`). The ordinal is preserved unchanged through the
/// Parallel section (`ExtractStage`, `BgzfCompress`) and consumed by
/// `BamFileWrite`'s `ReorderBuffer` to produce deterministic output order.
///
/// Downstream stages that do not care about ordering may ignore the field.
#[derive(Debug)]
pub struct TemplateBatch {
    pub templates: Vec<FastqTemplateV2>,
    pub ordinal: u64,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_template_batch_construct() {
        let t = FastqTemplateV2 { records: vec![], name: b"r1".to_vec() };
        let batch = TemplateBatch { templates: vec![t], ordinal: 42 };
        assert_eq!(batch.ordinal, 42);
        assert_eq!(batch.templates.len(), 1);
        assert_eq!(&batch.templates[0].name, b"r1");
    }
}

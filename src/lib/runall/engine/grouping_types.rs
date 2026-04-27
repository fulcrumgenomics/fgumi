//! Pipeline domain types for the group/consensus/filter chain.
//!
//! All byte-payload types in `pipeline` use the same length-prefixed
//! concat-byte encoding (see
//! [`crate::runall::engine::output_types::SerializedBatch`]). These domain
//! types add small typed metadata (position key, MI) on top of the same
//! encoding so stages can pattern-match on what a batch represents without
//! re-parsing the bytes.
//!
//! Record iteration uses the shared [`iter_length_prefixed`] helper which
//! walks a `&[u8]` as `[block_size u32 LE][record body][block_size …]`.
//!
//! The [`CoalesceInput`] trait is implemented by any typed batch that can be
//! flattened into a byte stream. `Coalesce` is generic over this trait so
//! that the same reorder-and-concatenate engine handles both byte-level
//! batches (`SerializedBatch`, emitted by Extract/Correct/AlignAndMerge/Sort)
//! and structured batches (`MiGroupBatch`, emitted by `GroupAssignStage`).
//!
//! ## MI tag injection
//!
//! `GroupAssignStage` does **not** inject the MI tag into the raw BAM record
//! bytes. It carries the local [`MoleculeId`] on each [`MiGroup`] instead.
//! A downstream serial-ordered stage applies the cumulative cross-batch
//! offset and writes the final MI tag bytes into the records.

use fgumi_umi::MoleculeId;

use crate::runall::engine::output_types::{RawBytes, SerializedBatch};

/// A batch of BAM records sharing the same genomic position key.
///
/// Emitted by
/// [`crate::runall::engine::stages::position_batch::PositionBatchStage`]; consumed
/// by [`crate::runall::engine::stages::group_assign::GroupAssignStage`].
#[derive(Debug, Clone)]
pub struct PositionGroupBatch {
    /// Length-prefixed concat records (body bytes only; each prefixed by a
    /// 4-byte little-endian `block_size`).
    pub data: Vec<u8>,
    /// Number of records encoded in [`data`](Self::data).
    pub record_count: u64,
    /// Monotonic batch ordinal assigned by `PositionBatchStage`.
    pub ordinal: u64,
    /// Template-coordinate position key: `(primary, secondary)`.
    pub position_key: (u64, u64),
}

/// A group of BAM records assigned to the same molecule ID.
///
/// Nested inside a [`MiGroupBatch`] — `GroupAssignStage` emits one
/// `MiGroupBatch` per input `PositionGroupBatch`, with each inner `MiGroup`
/// holding the records for a single MI.
#[derive(Debug, Clone)]
pub struct MiGroup {
    /// Length-prefixed concat records. **The MI tag is NOT yet injected**
    /// when emitted by `GroupAssignStage` — a downstream serial-ordered
    /// stage writes it after applying the cumulative cross-batch offset.
    /// Records flowing through the runall pipeline downstream of that stage
    /// carry the MI tag.
    pub data: Vec<u8>,
    /// Number of records encoded in [`data`](Self::data).
    pub record_count: u64,
    /// Numeric MI for ordering. Local (`0..N-1` per batch) when emitted by
    /// `GroupAssignStage`; replaced by the global value by the downstream
    /// serial-ordered MI assign stage.
    pub mi: u64,
    /// Local `MoleculeId` (variant tag preserved). The downstream MI assign
    /// stage calls `with_offset(base)` to convert to a global ID before
    /// writing the MI tag bytes. `None` after that stage runs (the offset
    /// has been applied and the value lives in `mi` and in the record bytes).
    pub local_mi: Option<MoleculeId>,
}

/// A batch of MI groups coming from one position group.
///
/// Emitted by `GroupAssignStage` carrying *local* MI integers; rewritten in
/// place by a downstream serial-ordered MI assign stage with global integers
/// and BAM-record MI tag bytes; consumed by `ConsensusStage` or
/// `Coalesce<MiGroupBatch>` on the `--stop-after group` path. The MI groups
/// are ordered by MI ascending so downstream consumers can assume a stable
/// within-batch ordering for byte-compare determinism.
#[derive(Debug, Clone)]
pub struct MiGroupBatch {
    /// The MI groups for this position, ordered by MI ascending.
    pub groups: Vec<MiGroup>,
    /// Monotonic batch ordinal, inherited from the incoming
    /// `PositionGroupBatch.ordinal` so downstream stages can reorder by
    /// pipeline-serial order.
    pub ordinal: u64,
    /// Template-coordinate position key: `(primary, secondary)`.
    pub position_key: (u64, u64),
    /// Two-byte BAM tag for the assigned molecule ID (e.g., `b"MI"`).
    /// Threaded from `GroupAssignStage` to the downstream MI assign stage
    /// so the latter knows which tag to write into each record.
    pub assign_tag: [u8; 2],
}

/// A typed pipeline batch that can be coalesced into bytes for BGZF
/// compression. Implemented by `SerializedBatch` trivially and by
/// `MiGroupBatch` by concatenating each inner `MiGroup.data` in order.
pub trait CoalesceInput: Send + 'static {
    /// The batch's ordinal — used by `Coalesce` to reorder arrivals before
    /// concatenation.
    fn ordinal(&self) -> u64;

    /// Consume the batch and return its byte representation:
    /// `(primary, optional secondary)`. The primary stream is the one that
    /// flows through `BgzfCompress → BamFileWrite` into the primary BAM
    /// output; the optional secondary stream is used for rejects BAMs
    /// (e.g., Correct's rejected reads).
    fn into_parts(self) -> (RawBytes, Option<RawBytes>);
}

impl CoalesceInput for SerializedBatch {
    fn ordinal(&self) -> u64 {
        self.ordinal
    }

    fn into_parts(self) -> (RawBytes, Option<RawBytes>) {
        (self.primary, self.secondary)
    }
}

impl CoalesceInput for MiGroupBatch {
    fn ordinal(&self) -> u64 {
        self.ordinal
    }

    fn into_parts(self) -> (RawBytes, Option<RawBytes>) {
        let total: usize = self.groups.iter().map(|g| g.data.len()).sum();
        let mut data = Vec::with_capacity(total);
        let mut record_count: u64 = 0;
        for g in self.groups {
            data.extend_from_slice(&g.data);
            record_count += g.record_count;
        }
        (RawBytes { data, record_count }, None)
    }
}

/// Iterator over length-prefixed records in a concat-byte buffer.
///
/// Yields the record body (without the `block_size` prefix) for each record,
/// or `Err` if the buffer is truncated mid-record.
#[must_use]
pub fn iter_length_prefixed(data: &[u8]) -> LengthPrefixedIter<'_> {
    LengthPrefixedIter { data, cursor: 0 }
}

/// Iterator produced by [`iter_length_prefixed`].
pub struct LengthPrefixedIter<'a> {
    data: &'a [u8],
    cursor: usize,
}

impl<'a> Iterator for LengthPrefixedIter<'a> {
    type Item = anyhow::Result<&'a [u8]>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.cursor >= self.data.len() {
            return None;
        }
        if self.data.len() < self.cursor + 4 {
            let err = anyhow::anyhow!(
                "iter_length_prefixed: truncated block_size at offset {}",
                self.cursor
            );
            self.cursor = self.data.len(); // stop further iteration
            return Some(Err(err));
        }
        let block_size =
            u32::from_le_bytes(self.data[self.cursor..self.cursor + 4].try_into().unwrap())
                as usize;
        let rec_start = self.cursor + 4;
        let rec_end = rec_start + block_size;
        if self.data.len() < rec_end {
            let err = anyhow::anyhow!(
                "iter_length_prefixed: truncated record at offset {} (block_size={block_size})",
                self.cursor
            );
            self.cursor = self.data.len();
            return Some(Err(err));
        }
        let slice = &self.data[rec_start..rec_end];
        self.cursor = rec_end;
        Some(Ok(slice))
    }
}

/// Serialize a slice of per-record byte buffers into the concat length-prefixed
/// format used by [`PositionGroupBatch`] and [`MiGroup`].
///
/// # Errors
///
/// Returns an error if any record exceeds `u32::MAX` bytes.
pub fn serialize_records(records: &[Vec<u8>]) -> anyhow::Result<(Vec<u8>, u64)> {
    let total: usize = records.iter().map(|r| 4 + r.len()).sum();
    let mut data = Vec::with_capacity(total);
    let mut count: u64 = 0;
    for rec in records {
        let block_size = u32::try_from(rec.len()).map_err(|_| {
            anyhow::anyhow!("serialize_records: record exceeds u32::MAX bytes ({})", rec.len())
        })?;
        data.extend_from_slice(&block_size.to_le_bytes());
        data.extend_from_slice(rec);
        count += 1;
    }
    Ok((data, count))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_iter_empty() {
        let it = iter_length_prefixed(&[]);
        let v: Vec<_> = it.collect();
        assert!(v.is_empty());
    }

    #[test]
    fn test_iter_single_record() {
        let mut buf = Vec::new();
        let rec: &[u8] = b"HELLO";
        buf.extend_from_slice(&u32::try_from(rec.len()).unwrap().to_le_bytes());
        buf.extend_from_slice(rec);
        let v: Vec<&[u8]> = iter_length_prefixed(&buf).map(Result::unwrap).collect();
        assert_eq!(v, vec![rec]);
    }

    #[test]
    fn test_iter_multiple_records() {
        let recs: [&[u8]; 3] = [b"A", b"BB", b"CCC"];
        let (data, count) =
            serialize_records(&recs.iter().map(|r| r.to_vec()).collect::<Vec<_>>()).unwrap();
        assert_eq!(count, 3);
        let parsed: Vec<&[u8]> = iter_length_prefixed(&data).map(Result::unwrap).collect();
        assert_eq!(parsed, recs);
    }

    #[test]
    fn test_iter_truncated_prefix() {
        let buf = vec![0x01, 0x00]; // only 2 bytes, need 4 for prefix
        let mut it = iter_length_prefixed(&buf);
        assert!(matches!(it.next(), Some(Err(_))));
        assert!(it.next().is_none(), "iterator stops after error");
    }

    #[test]
    fn test_iter_truncated_record() {
        let mut buf = Vec::new();
        buf.extend_from_slice(&10u32.to_le_bytes()); // claims 10 bytes
        buf.extend_from_slice(b"ABC"); // only 3 present
        let mut it = iter_length_prefixed(&buf);
        assert!(matches!(it.next(), Some(Err(_))));
    }

    #[test]
    fn test_round_trip() {
        let records: Vec<Vec<u8>> = vec![b"alpha".to_vec(), b"beta".to_vec(), b"gamma".to_vec()];
        let (data, count) = serialize_records(&records).unwrap();
        assert_eq!(count, records.len() as u64);
        let out: Vec<Vec<u8>> = iter_length_prefixed(&data).map(|r| r.unwrap().to_vec()).collect();
        assert_eq!(out, records);
    }
}

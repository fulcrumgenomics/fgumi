//! Molecular Identifier (MI) group utilities for streaming BAM processing.
//!
//! This module provides an iterator for grouping consecutive BAM records by their MI tag.
//! This is useful for streaming consensus calling where input is already sorted by MI groups
//! (e.g., output from `group`).

use anyhow::Result;
use fgumi_raw_bam::RawRecord;
use std::io;

use crate::unified_pipeline::{BatchWeight, DecodedRecord, MemoryEstimate};

use crate::unified_pipeline::Grouper;
use std::collections::VecDeque;

// ============================================================================
// Raw-byte MI grouping for consensus callers
// ============================================================================

/// A single MI group holding raw-byte BAM records.
#[derive(Debug, Clone)]
pub struct MiGroup {
    /// The MI tag value (e.g., "0", "1/A", "1/B")
    pub mi: String,
    /// Raw BAM records sharing this MI value
    pub records: Vec<RawRecord>,
}

impl MiGroup {
    /// Creates a new raw MI group.
    #[must_use]
    pub fn new(mi: String, records: Vec<RawRecord>) -> Self {
        Self { mi, records }
    }
}

impl BatchWeight for MiGroup {
    fn batch_weight(&self) -> usize {
        self.records.len()
    }
}

impl MemoryEstimate for MiGroup {
    fn estimate_heap_size(&self) -> usize {
        let mi_size = self.mi.capacity();
        let records_size: usize = self.records.iter().map(RawRecord::capacity).sum();
        let records_vec_overhead = self.records.capacity() * std::mem::size_of::<RawRecord>();
        mi_size + records_size + records_vec_overhead
    }
}

/// A batch of raw MI groups for parallel processing.
#[derive(Default)]
pub struct MiGroupBatch {
    /// The raw MI groups in this batch
    pub groups: Vec<MiGroup>,
}

impl MiGroupBatch {
    /// Creates a new empty raw MI group batch.
    #[must_use]
    pub fn new() -> Self {
        Self { groups: Vec::new() }
    }

    /// Creates a new raw MI group batch with pre-allocated capacity.
    #[must_use]
    pub fn with_capacity(capacity: usize) -> Self {
        Self { groups: Vec::with_capacity(capacity) }
    }

    /// Returns the number of groups in the batch.
    #[must_use]
    pub fn len(&self) -> usize {
        self.groups.len()
    }

    /// Returns true if the batch contains no groups.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.groups.is_empty()
    }

    /// Clears all groups from the batch.
    pub fn clear(&mut self) {
        self.groups.clear();
    }
}

impl BatchWeight for MiGroupBatch {
    fn batch_weight(&self) -> usize {
        self.groups.iter().map(|g| g.records.len()).sum()
    }
}

impl MemoryEstimate for MiGroupBatch {
    fn estimate_heap_size(&self) -> usize {
        let groups_size: usize = self.groups.iter().map(MemoryEstimate::estimate_heap_size).sum();
        let groups_vec_overhead = self.groups.capacity() * std::mem::size_of::<MiGroup>();
        groups_size + groups_vec_overhead
    }
}

/// Type alias for raw-byte MI tag transformation function.
type MiTransformFn = Box<dyn Fn(&[u8]) -> String + Send + Sync>;

/// Type alias for raw-byte record filter function.
type RecordFilterFn = Box<dyn Fn(&[u8]) -> bool + Send + Sync>;

/// A Grouper that groups raw-byte BAM records by MI tag.
///
/// Records arrive as raw BAM bytes (see [`DecodedRecord::from_raw_bytes`]). MI tags are
/// extracted directly from the raw bytes using `bam_fields::find_string_tag_in_record()`
/// without parsing into `RecordBuf`.
pub struct MiGrouper {
    /// The MI tag bytes to search for
    tag: [u8; 2],
    /// Number of MI groups per batch
    batch_size: usize,
    /// Current MI value being accumulated
    current_mi: Option<String>,
    /// Records in current MI group (raw bytes)
    current_records: Vec<RawRecord>,
    /// Completed groups waiting to be batched
    pending_groups: VecDeque<MiGroup>,
    /// Whether `finish()` has been called
    finished: bool,
    /// Optional MI tag transformation function
    mi_transform: Option<MiTransformFn>,
    /// Optional record filter
    record_filter: Option<RecordFilterFn>,
    /// Optional cell barcode tag for composite grouping (MI + cell barcode)
    cell_tag: Option<[u8; 2]>,
}

impl MiGrouper {
    /// Create a new `MiGrouper`.
    ///
    /// # Panics
    ///
    /// Panics if `tag_name` is not exactly 2 characters.
    #[must_use]
    pub fn new(tag_name: &str, batch_size: usize) -> Self {
        assert!(tag_name.len() == 2, "Tag name must be exactly 2 characters");
        let tag_bytes = tag_name.as_bytes();

        Self {
            tag: [tag_bytes[0], tag_bytes[1]],
            batch_size: batch_size.max(1),
            current_mi: None,
            current_records: Vec::new(),
            pending_groups: VecDeque::new(),
            finished: false,
            mi_transform: None,
            record_filter: None,
            cell_tag: None,
        }
    }

    /// Create a `MiGrouper` with record filtering and MI tag transformation.
    ///
    /// # Panics
    ///
    /// Panics if `tag_name` is not exactly 2 characters.
    pub fn with_filter_and_transform<F, T>(
        tag_name: &str,
        batch_size: usize,
        record_filter: F,
        mi_transform: T,
    ) -> Self
    where
        F: Fn(&[u8]) -> bool + Send + Sync + 'static,
        T: Fn(&[u8]) -> String + Send + Sync + 'static,
    {
        assert!(tag_name.len() == 2, "Tag name must be exactly 2 characters");
        let tag_bytes = tag_name.as_bytes();

        Self {
            tag: [tag_bytes[0], tag_bytes[1]],
            batch_size: batch_size.max(1),
            current_mi: None,
            current_records: Vec::new(),
            pending_groups: VecDeque::new(),
            finished: false,
            mi_transform: Some(Box::new(mi_transform)),
            record_filter: Some(Box::new(record_filter)),
            cell_tag: None,
        }
    }

    /// Sets an optional cell barcode tag for composite grouping.
    ///
    /// When set, the group key becomes `"MI_VALUE\tCELL_VALUE"` so that reads from
    /// different cells with the same MI tag are placed in separate groups.
    #[must_use]
    pub fn with_cell_tag(mut self, cell_tag: Option<[u8; 2]>) -> Self {
        self.cell_tag = cell_tag;
        self
    }

    /// Get MI tag from raw BAM bytes, optionally applying transformation.
    ///
    /// When `cell_tag` is set, returns a composite key of `"MI\tCELL"`.
    fn get_mi_tag(&self, bam: &[u8]) -> Option<String> {
        use crate::sort::bam_fields;
        let value = bam_fields::find_string_tag_in_record(bam, self.tag)?;
        let mut key = if let Some(ref transform) = self.mi_transform {
            transform(value)
        } else {
            String::from_utf8_lossy(value).into_owned()
        };
        if let Some(ct) = &self.cell_tag {
            key.push('\t');
            if let Some(cell_value) = bam_fields::find_string_tag_in_record(bam, ct) {
                key.push_str(&String::from_utf8_lossy(cell_value));
            }
        }
        Some(key)
    }

    /// Check if a raw record passes the filter.
    fn should_keep(&self, bam: &[u8]) -> bool {
        match &self.record_filter {
            Some(filter) => filter(bam),
            None => true,
        }
    }

    /// Flush current MI group to pending.
    fn flush_current_group(&mut self) {
        if let Some(mi) = self.current_mi.take() {
            if !self.current_records.is_empty() {
                let records = std::mem::take(&mut self.current_records);
                self.pending_groups.push_back(MiGroup::new(mi, records));
            }
        }
    }

    /// Try to form complete batches from pending groups.
    fn drain_batches(&mut self) -> Vec<MiGroupBatch> {
        let mut batches = Vec::new();
        while self.pending_groups.len() >= self.batch_size {
            let groups: Vec<MiGroup> = self.pending_groups.drain(..self.batch_size).collect();
            batches.push(MiGroupBatch { groups });
        }
        batches
    }
}

impl Grouper for MiGrouper {
    type Group = MiGroupBatch;

    fn add_records(&mut self, records: Vec<DecodedRecord>) -> io::Result<Vec<Self::Group>> {
        for decoded in records {
            let raw = decoded.into_raw_bytes();

            // Apply record filter if configured
            if !self.should_keep(&raw) {
                continue;
            }

            // Skip records without an MI tag: the consensus caller requires the
            // MI tag and would fail when processing them. `MiGroupIterator`
            // applies the same skip policy for parity across paths.
            let Some(mi) = self.get_mi_tag(&raw) else {
                continue;
            };

            // Check if this starts a new MI group
            match &self.current_mi {
                Some(current) if current == &mi => {
                    self.current_records.push(raw);
                }
                Some(_) => {
                    self.flush_current_group();
                    self.current_mi = Some(mi);
                    self.current_records.push(raw);
                }
                None => {
                    self.current_mi = Some(mi);
                    self.current_records.push(raw);
                }
            }
        }

        Ok(self.drain_batches())
    }

    fn finish(&mut self) -> io::Result<Option<Self::Group>> {
        if self.finished {
            return Ok(None);
        }
        self.finished = true;

        self.flush_current_group();

        if self.pending_groups.is_empty() {
            Ok(None)
        } else {
            let groups: Vec<MiGroup> = self.pending_groups.drain(..).collect();
            Ok(Some(MiGroupBatch { groups }))
        }
    }

    fn has_pending(&self) -> bool {
        !self.pending_groups.is_empty() || self.current_mi.is_some()
    }
}

/// An iterator that groups consecutive raw BAM records by their MI tag, yielding
/// one group per distinct MI value from the input stream.
#[allow(clippy::type_complexity)]
pub struct MiGroupIterator<I>
where
    I: Iterator<Item = Result<RawRecord>>,
{
    record_iter: I,
    tag: [u8; 2],
    /// Optional cell barcode tag for composite grouping (MI + cell barcode)
    cell_tag: Option<[u8; 2]>,
    current_mi: Option<String>,
    current_group: Vec<RawRecord>,
    done: bool,
    /// Pending error to return after flushing a group
    pending_error: Option<anyhow::Error>,
    /// Optional MI tag transformation function
    mi_transform: Option<Box<dyn Fn(&[u8]) -> String>>,
}

impl<I> MiGroupIterator<I>
where
    I: Iterator<Item = Result<RawRecord>>,
{
    /// Creates a new `MiGroupIterator`.
    ///
    /// # Panics
    ///
    /// Panics if `tag_name` is not exactly 2 characters.
    pub fn new(record_iter: I, tag_name: &str) -> Self {
        assert!(tag_name.len() == 2, "Tag name must be exactly 2 characters");
        let tag_bytes = tag_name.as_bytes();
        Self {
            record_iter,
            tag: [tag_bytes[0], tag_bytes[1]],
            cell_tag: None,
            current_mi: None,
            current_group: Vec::new(),
            done: false,
            pending_error: None,
            mi_transform: None,
        }
    }

    /// Creates a new `MiGroupIterator` with MI tag transformation.
    ///
    /// # Panics
    ///
    /// Panics if `tag_name` is not exactly 2 characters.
    pub fn with_transform<F>(record_iter: I, tag_name: &str, mi_transform: F) -> Self
    where
        F: Fn(&[u8]) -> String + 'static,
    {
        assert!(tag_name.len() == 2, "Tag name must be exactly 2 characters");
        let tag_bytes = tag_name.as_bytes();
        Self {
            record_iter,
            tag: [tag_bytes[0], tag_bytes[1]],
            cell_tag: None,
            current_mi: None,
            current_group: Vec::new(),
            done: false,
            pending_error: None,
            mi_transform: Some(Box::new(mi_transform)),
        }
    }

    /// Sets an optional cell barcode tag for composite grouping.
    ///
    /// When set, the group key becomes `"MI_VALUE\tCELL_VALUE"` so that reads from
    /// different cells with the same MI tag are placed in separate groups.
    #[must_use]
    pub fn with_cell_tag(mut self, cell_tag: Option<[u8; 2]>) -> Self {
        self.cell_tag = cell_tag;
        self
    }

    /// Extracts the group key from raw BAM bytes.
    ///
    /// When `cell_tag` is set, returns a composite key of `"MI\tCELL"`.
    /// Otherwise returns just the MI tag value.
    fn get_mi(&self, bam: &[u8]) -> Option<String> {
        use crate::sort::bam_fields;
        let value = bam_fields::find_string_tag_in_record(bam, self.tag)?;
        let mut key = if let Some(ref transform) = self.mi_transform {
            transform(value)
        } else {
            String::from_utf8_lossy(value).into_owned()
        };
        if let Some(ct) = &self.cell_tag {
            key.push('\t');
            if let Some(cell_value) = bam_fields::find_string_tag_in_record(bam, ct) {
                key.push_str(&String::from_utf8_lossy(cell_value));
            }
        }
        Some(key)
    }
}

impl<I> Iterator for MiGroupIterator<I>
where
    I: Iterator<Item = Result<RawRecord>>,
{
    type Item = Result<(String, Vec<RawRecord>)>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }

        // Check for pending error (return after we've flushed any pending group)
        if let Some(e) = self.pending_error.take() {
            self.done = true;
            return Some(Err(e));
        }

        loop {
            match self.record_iter.next() {
                None => {
                    self.done = true;
                    if self.current_group.is_empty() {
                        return None;
                    }
                    let mi = self.current_mi.take().unwrap_or_default();
                    let group = std::mem::take(&mut self.current_group);
                    return Some(Ok((mi, group)));
                }
                Some(Err(e)) => {
                    if !self.current_group.is_empty() {
                        self.pending_error = Some(e);
                        let mi = self.current_mi.take().unwrap_or_default();
                        let group = std::mem::take(&mut self.current_group);
                        return Some(Ok((mi, group)));
                    }
                    self.done = true;
                    return Some(Err(e));
                }
                Some(Ok(raw)) => {
                    // Skip records without an MI tag; the consensus caller requires
                    // an MI tag and would fail when processing them. `MiGrouper`
                    // applies the same skip policy for parity across paths.
                    let Some(mi) = self.get_mi(&raw) else {
                        continue;
                    };

                    if self.current_group.is_empty() {
                        self.current_mi = Some(mi);
                        self.current_group.push(raw);
                    } else if self.current_mi.as_ref() == Some(&mi) {
                        self.current_group.push(raw);
                    } else {
                        let old_mi = self.current_mi.take().unwrap_or_default();
                        let group = std::mem::take(&mut self.current_group);
                        self.current_mi = Some(mi);
                        self.current_group.push(raw);
                        return Some(Ok((old_mi, group)));
                    }
                }
            }
        }
    }
}

#[cfg(test)]
#[allow(clippy::similar_names)]
mod tests {
    use super::*;
    use crate::sam::SamTag;
    use crate::umi::extract_mi_base;

    // ========================================================================
    // Helpers for raw BAM byte construction
    // ========================================================================

    // Build a minimal raw BAM record with the given string tag.
    //
    // The BAM binary format (after the 4-byte block_size prefix, which we omit):
    //   bytes 0..4   : refID  (i32 LE)
    //   bytes 4..8   : pos    (i32 LE)
    //   byte  8      : l_read_name (u8, includes NUL)
    //   byte  9      : mapq
    //   bytes 10..12 : bin    (u16 LE)
    //   bytes 12..14 : n_cigar_op (u16 LE)
    //   bytes 14..16 : flag   (u16 LE)
    //   bytes 16..20 : l_seq  (u32 LE)
    //   bytes 20..24 : next_refID (i32 LE)
    //   bytes 24..28 : next_pos   (i32 LE)
    //   bytes 28..32 : tlen       (i32 LE)
    //   then: read_name (l_read_name bytes, NUL-terminated)
    //   then: cigar    (n_cigar_op * 4 bytes)
    //   then: seq      (ceil(l_seq/2) bytes, 4-bit encoded)
    //   then: qual     (l_seq bytes)
    //   then: aux data
    #[allow(clippy::cast_possible_truncation)]
    fn make_raw_bam_with_tag(tag_name: &str, tag_value: &str) -> RawRecord {
        let name = b"read";
        let l_read_name: u8 = (name.len() + 1) as u8; // +1 for NUL
        let seq_len: u32 = 4; // 4 bases (ACGT)
        let seq_bytes = seq_len.div_ceil(2) as usize;

        // Aux data: tag[0] tag[1] 'Z' value NUL
        let tag_bytes = tag_name.as_bytes();
        let aux: Vec<u8> =
            [&[tag_bytes[0], tag_bytes[1], b'Z'], tag_value.as_bytes(), &[0u8]].concat();

        let total = 32 + l_read_name as usize + seq_bytes + seq_len as usize + aux.len();
        let mut buf = vec![0u8; total];

        // refID = -1 (unmapped)
        buf[0..4].copy_from_slice(&(-1i32).to_le_bytes());
        // pos = -1 (unmapped)
        buf[4..8].copy_from_slice(&(-1i32).to_le_bytes());
        // l_read_name
        buf[8] = l_read_name;
        // n_cigar_op = 0
        buf[12..14].copy_from_slice(&0u16.to_le_bytes());
        // l_seq
        buf[16..20].copy_from_slice(&seq_len.to_le_bytes());
        // next_refID = -1
        buf[20..24].copy_from_slice(&(-1i32).to_le_bytes());
        // next_pos = -1
        buf[24..28].copy_from_slice(&(-1i32).to_le_bytes());

        // read_name + NUL
        let name_start = 32;
        buf[name_start..name_start + name.len()].copy_from_slice(name);
        buf[name_start + name.len()] = 0;

        // seq and qual are left as zeros

        // aux data
        let aux_start = 32 + l_read_name as usize + seq_bytes + seq_len as usize;
        buf[aux_start..aux_start + aux.len()].copy_from_slice(&aux);

        RawRecord::from(buf)
    }

    /// Build a minimal raw BAM record with no aux tags.
    #[allow(clippy::cast_possible_truncation)]
    fn make_raw_bam_without_tag() -> RawRecord {
        let name = b"read";
        let l_read_name: u8 = (name.len() + 1) as u8;
        let seq_len: u32 = 4;
        let seq_bytes = seq_len.div_ceil(2) as usize;

        let total = 32 + l_read_name as usize + seq_bytes + seq_len as usize;
        let mut buf = vec![0u8; total];

        buf[0..4].copy_from_slice(&(-1i32).to_le_bytes());
        buf[4..8].copy_from_slice(&(-1i32).to_le_bytes());
        buf[8] = l_read_name;
        buf[12..14].copy_from_slice(&0u16.to_le_bytes());
        buf[16..20].copy_from_slice(&seq_len.to_le_bytes());
        buf[20..24].copy_from_slice(&(-1i32).to_le_bytes());
        buf[24..28].copy_from_slice(&(-1i32).to_le_bytes());

        let name_start = 32;
        buf[name_start..name_start + name.len()].copy_from_slice(name);
        buf[name_start + name.len()] = 0;

        RawRecord::from(buf)
    }

    // ========================================================================
    // MiGroupIterator tests
    // ========================================================================

    #[test]
    fn test_empty_iterator() {
        let records: Vec<Result<RawRecord>> = vec![];
        let mut iter = MiGroupIterator::new(records.into_iter(), "MI");
        assert!(iter.next().is_none());
    }

    #[test]
    fn test_single_group() {
        let records: Vec<Result<RawRecord>> = vec![
            Ok(make_raw_bam_with_tag("MI", "0")),
            Ok(make_raw_bam_with_tag("MI", "0")),
            Ok(make_raw_bam_with_tag("MI", "0")),
        ];
        let mut iter = MiGroupIterator::new(records.into_iter(), "MI");

        let result = iter.next().expect("iterator should yield item").expect("item should be Ok");
        assert_eq!(result.0, "0");
        assert_eq!(result.1.len(), 3);

        assert!(iter.next().is_none());
    }

    #[test]
    fn test_multiple_groups() {
        let records: Vec<Result<RawRecord>> = vec![
            Ok(make_raw_bam_with_tag("MI", "0")),
            Ok(make_raw_bam_with_tag("MI", "0")),
            Ok(make_raw_bam_with_tag("MI", "1")),
            Ok(make_raw_bam_with_tag("MI", "1")),
            Ok(make_raw_bam_with_tag("MI", "1")),
            Ok(make_raw_bam_with_tag("MI", "2")),
        ];
        let mut iter = MiGroupIterator::new(records.into_iter(), "MI");

        let result = iter.next().expect("iterator should yield item").expect("item should be Ok");
        assert_eq!(result.0, "0");
        assert_eq!(result.1.len(), 2);

        let result = iter.next().expect("iterator should yield item").expect("item should be Ok");
        assert_eq!(result.0, "1");
        assert_eq!(result.1.len(), 3);

        let result = iter.next().expect("iterator should yield item").expect("item should be Ok");
        assert_eq!(result.0, "2");
        assert_eq!(result.1.len(), 1);

        assert!(iter.next().is_none());
    }

    #[test]
    fn test_skips_records_without_mi_tag() {
        let records: Vec<Result<RawRecord>> = vec![
            Ok(make_raw_bam_with_tag("MI", "0")),
            Ok(make_raw_bam_without_tag()),
            Ok(make_raw_bam_with_tag("MI", "0")),
            Ok(make_raw_bam_without_tag()),
            Ok(make_raw_bam_with_tag("MI", "1")),
        ];
        let mut iter = MiGroupIterator::new(records.into_iter(), "MI");

        let result = iter.next().expect("iterator should yield item").expect("item should be Ok");
        assert_eq!(result.0, "0");
        assert_eq!(result.1.len(), 2); // Skipped records without MI

        let result = iter.next().expect("iterator should yield item").expect("item should be Ok");
        assert_eq!(result.0, "1");
        assert_eq!(result.1.len(), 1);

        assert!(iter.next().is_none());
    }

    #[test]
    fn test_error_propagation() {
        let records: Vec<Result<RawRecord>> = vec![
            Ok(make_raw_bam_with_tag("MI", "0")),
            Err(anyhow::anyhow!("test error")),
            Ok(make_raw_bam_with_tag("MI", "1")),
        ];
        let mut iter = MiGroupIterator::new(records.into_iter(), "MI");

        // First group before error
        let result = iter.next().expect("iterator should yield item").expect("item should be Ok");
        assert_eq!(result.0, "0");
        assert_eq!(result.1.len(), 1);

        // Error should be returned after flushing the pending group
        let err = iter.next().expect("iterator should yield item");
        assert!(err.is_err());

        assert!(iter.next().is_none());
    }

    #[test]
    fn test_error_with_no_pending_group() {
        let records: Vec<Result<RawRecord>> =
            vec![Err(anyhow::anyhow!("immediate error")), Ok(make_raw_bam_with_tag("MI", "0"))];
        let mut iter = MiGroupIterator::new(records.into_iter(), "MI");

        // Error should be returned directly
        let result = iter.next().expect("iterator should yield item");
        assert!(result.is_err());

        // Iterator should be done after error
        assert!(iter.next().is_none());
    }

    #[test]
    fn test_custom_tag() {
        let records: Vec<Result<RawRecord>> =
            vec![Ok(make_raw_bam_with_tag("RX", "ACGT")), Ok(make_raw_bam_with_tag("RX", "ACGT"))];
        let mut iter = MiGroupIterator::new(records.into_iter(), "RX");

        let result = iter.next().expect("iterator should yield item").expect("item should be Ok");
        assert_eq!(result.0, "ACGT");
        assert_eq!(result.1.len(), 2);

        assert!(iter.next().is_none());
    }

    #[test]
    #[should_panic(expected = "Tag name must be exactly 2 characters")]
    fn test_invalid_tag_length() {
        let records: Vec<Result<RawRecord>> = vec![];
        let _ = MiGroupIterator::new(records.into_iter(), "M");
    }

    #[test]
    fn test_with_transform() {
        // Simulate duplex reads: 1/A, 1/B should group under "1"
        let records: Vec<Result<RawRecord>> = vec![
            Ok(make_raw_bam_with_tag("MI", "1/A")),
            Ok(make_raw_bam_with_tag("MI", "1/A")),
            Ok(make_raw_bam_with_tag("MI", "1/B")),
            Ok(make_raw_bam_with_tag("MI", "1/B")),
            Ok(make_raw_bam_with_tag("MI", "2/A")),
            Ok(make_raw_bam_with_tag("MI", "2/B")),
        ];
        let mut iter = MiGroupIterator::with_transform(records.into_iter(), "MI", |raw| {
            let s = String::from_utf8_lossy(raw);
            extract_mi_base(&s).to_string()
        });

        // First group: base MI "1" with 4 reads
        let result = iter.next().expect("iterator should yield item").expect("item should be Ok");
        assert_eq!(result.0, "1");
        assert_eq!(result.1.len(), 4);

        // Second group: base MI "2" with 2 reads
        let result = iter.next().expect("iterator should yield item").expect("item should be Ok");
        assert_eq!(result.0, "2");
        assert_eq!(result.1.len(), 2);

        assert!(iter.next().is_none());
    }

    #[test]
    fn test_with_transform_empty() {
        let records: Vec<Result<RawRecord>> = vec![];
        let mut iter = MiGroupIterator::with_transform(records.into_iter(), "MI", |raw| {
            String::from_utf8_lossy(raw).to_uppercase()
        });
        assert!(iter.next().is_none());
    }

    #[test]
    #[should_panic(expected = "Tag name must be exactly 2 characters")]
    fn test_with_transform_invalid_tag_length() {
        let records: Vec<Result<RawRecord>> = vec![];
        let _ = MiGroupIterator::with_transform(records.into_iter(), "ABC", |raw| {
            String::from_utf8_lossy(raw).into_owned()
        });
    }

    #[test]
    fn test_get_mi_without_transform() {
        let iter = MiGroupIterator::new(std::iter::empty::<Result<RawRecord>>(), "MI");
        let bam = make_raw_bam_with_tag("MI", "42");
        assert_eq!(iter.get_mi(&bam), Some("42".to_string()));
    }

    #[test]
    fn test_get_mi_with_transform() {
        let iter =
            MiGroupIterator::with_transform(std::iter::empty::<Result<RawRecord>>(), "MI", |raw| {
                let s = String::from_utf8_lossy(raw);
                s.to_uppercase()
            });
        let bam = make_raw_bam_with_tag("MI", "abc");
        assert_eq!(iter.get_mi(&bam), Some("ABC".to_string()));
    }

    #[test]
    fn test_get_mi_missing_tag() {
        let iter = MiGroupIterator::new(std::iter::empty::<Result<RawRecord>>(), "MI");
        let bam = make_raw_bam_without_tag();
        assert_eq!(iter.get_mi(&bam), None);
    }

    #[test]
    fn test_get_mi_wrong_tag() {
        let iter = MiGroupIterator::new(std::iter::empty::<Result<RawRecord>>(), "MI");
        let bam = make_raw_bam_with_tag("RX", "ACGT");
        assert_eq!(iter.get_mi(&bam), None);
    }

    /// Build a minimal raw BAM record with two string tags.
    #[allow(clippy::cast_possible_truncation)]
    fn make_raw_bam_with_two_tags(tag1: &str, val1: &str, tag2: &str, val2: &str) -> RawRecord {
        let name = b"read";
        let l_read_name: u8 = (name.len() + 1) as u8;
        let seq_len: u32 = 4;
        let seq_bytes = seq_len.div_ceil(2) as usize;

        let t1 = tag1.as_bytes();
        let t2 = tag2.as_bytes();
        let aux: Vec<u8> = [
            &[t1[0], t1[1], b'Z'],
            val1.as_bytes(),
            &[0u8],
            &[t2[0], t2[1], b'Z'],
            val2.as_bytes(),
            &[0u8],
        ]
        .concat();

        let total = 32 + l_read_name as usize + seq_bytes + seq_len as usize + aux.len();
        let mut buf = vec![0u8; total];

        buf[0..4].copy_from_slice(&(-1i32).to_le_bytes());
        buf[4..8].copy_from_slice(&(-1i32).to_le_bytes());
        buf[8] = l_read_name;
        buf[12..14].copy_from_slice(&0u16.to_le_bytes());
        buf[16..20].copy_from_slice(&seq_len.to_le_bytes());
        buf[20..24].copy_from_slice(&(-1i32).to_le_bytes());
        buf[24..28].copy_from_slice(&(-1i32).to_le_bytes());

        let name_start = 32;
        buf[name_start..name_start + name.len()].copy_from_slice(name);
        buf[name_start + name.len()] = 0;

        let aux_start = 32 + l_read_name as usize + seq_bytes + seq_len as usize;
        buf[aux_start..aux_start + aux.len()].copy_from_slice(&aux);

        RawRecord::from(buf)
    }

    #[test]
    fn test_cell_tag_composite_grouping() {
        // Same MI but different cell barcodes should form separate groups
        let records: Vec<Result<RawRecord>> = vec![
            Ok(make_raw_bam_with_two_tags("MI", "1", "CB", "ACGT")),
            Ok(make_raw_bam_with_two_tags("MI", "1", "CB", "ACGT")),
            Ok(make_raw_bam_with_two_tags("MI", "1", "CB", "TGCA")),
            Ok(make_raw_bam_with_two_tags("MI", "1", "CB", "TGCA")),
        ];
        let mut iter =
            MiGroupIterator::new(records.into_iter(), "MI").with_cell_tag(Some(*SamTag::CB));

        // First group: MI=1, CB=ACGT
        let result = iter.next().expect("iterator should yield item").expect("item should be Ok");
        assert_eq!(result.0, "1\tACGT");
        assert_eq!(result.1.len(), 2);

        // Second group: MI=1, CB=TGCA
        let result = iter.next().expect("iterator should yield item").expect("item should be Ok");
        assert_eq!(result.0, "1\tTGCA");
        assert_eq!(result.1.len(), 2);

        assert!(iter.next().is_none());
    }

    #[test]
    fn test_cell_tag_none_groups_by_mi_only() {
        // Without cell_tag, same MI records group together regardless of other tags
        let records: Vec<Result<RawRecord>> = vec![
            Ok(make_raw_bam_with_two_tags("MI", "1", "CB", "ACGT")),
            Ok(make_raw_bam_with_two_tags("MI", "1", "CB", "TGCA")),
        ];
        let mut iter = MiGroupIterator::new(records.into_iter(), "MI").with_cell_tag(None);

        let result = iter.next().expect("iterator should yield item").expect("item should be Ok");
        assert_eq!(result.0, "1");
        assert_eq!(result.1.len(), 2); // Both records in same group

        assert!(iter.next().is_none());
    }

    #[test]
    fn test_cell_tag_missing_cell_value() {
        // Records with MI but no cell tag get grouped under "MI\t" (empty cell)
        let records: Vec<Result<RawRecord>> = vec![
            Ok(make_raw_bam_with_tag("MI", "1")),
            Ok(make_raw_bam_with_tag("MI", "1")),
            Ok(make_raw_bam_with_two_tags("MI", "1", "CB", "ACGT")),
        ];
        let mut iter =
            MiGroupIterator::new(records.into_iter(), "MI").with_cell_tag(Some(*SamTag::CB));

        // First group: MI=1, no CB (key = "1\t")
        let result = iter.next().expect("iterator should yield item").expect("item should be Ok");
        assert_eq!(result.0, "1\t");
        assert_eq!(result.1.len(), 2);

        // Second group: MI=1, CB=ACGT (key = "1\tACGT")
        let result = iter.next().expect("iterator should yield item").expect("item should be Ok");
        assert_eq!(result.0, "1\tACGT");
        assert_eq!(result.1.len(), 1);

        assert!(iter.next().is_none());
    }

    #[test]
    fn test_cell_tag_with_transform() {
        // Cell tag should work with MI transform (duplex-style /A /B stripping)
        let records: Vec<Result<RawRecord>> = vec![
            Ok(make_raw_bam_with_two_tags("MI", "1/A", "CB", "ACGT")),
            Ok(make_raw_bam_with_two_tags("MI", "1/B", "CB", "ACGT")),
            Ok(make_raw_bam_with_two_tags("MI", "1/A", "CB", "TGCA")),
        ];
        let mut iter = MiGroupIterator::with_transform(records.into_iter(), "MI", |raw| {
            let s = String::from_utf8_lossy(raw);
            extract_mi_base(&s).to_string()
        })
        .with_cell_tag(Some(*SamTag::CB));

        // All three have base MI "1", but different cells
        // First group: base MI=1, CB=ACGT (2 records: 1/A and 1/B)
        let result = iter.next().expect("iterator should yield item").expect("item should be Ok");
        assert_eq!(result.0, "1\tACGT");
        assert_eq!(result.1.len(), 2);

        // Second group: base MI=1, CB=TGCA
        let result = iter.next().expect("iterator should yield item").expect("item should be Ok");
        assert_eq!(result.0, "1\tTGCA");
        assert_eq!(result.1.len(), 1);

        assert!(iter.next().is_none());
    }

    // ========================================================================
    // MiGrouper tests (Grouper trait impl)
    // ========================================================================

    /// Helper: create a `DecodedRecord` from raw bytes with a dummy group key.
    fn make_raw_decoded_record(tag_name: &str, tag_value: &str) -> DecodedRecord {
        let raw = make_raw_bam_with_tag(tag_name, tag_value);
        let key = crate::unified_pipeline::GroupKey::single(0, 0, 0, 0, 0, 0);
        DecodedRecord::from_raw_bytes(raw, key)
    }

    /// Helper: create a `DecodedRecord` from raw bytes with no tags.
    fn make_raw_decoded_record_no_tag() -> DecodedRecord {
        let raw = make_raw_bam_without_tag();
        let key = crate::unified_pipeline::GroupKey::single(0, 0, 0, 0, 0, 0);
        DecodedRecord::from_raw_bytes(raw, key)
    }

    #[test]
    fn test_grouper_single_mi_group() {
        let mut grouper = MiGrouper::new("MI", 10);

        let records = vec![
            make_raw_decoded_record("MI", "0"),
            make_raw_decoded_record("MI", "0"),
            make_raw_decoded_record("MI", "0"),
        ];

        let batches = grouper.add_records(records).expect("add_records should succeed");
        // Only 1 MI group so far, batch_size=10, no complete batch yet
        assert!(batches.is_empty());
        assert!(grouper.has_pending());

        let final_batch =
            grouper.finish().expect("finish should succeed").expect("should return final batch");
        assert_eq!(final_batch.groups.len(), 1);
        assert_eq!(final_batch.groups[0].mi, "0");
        assert_eq!(final_batch.groups[0].records.len(), 3);
    }

    #[test]
    fn test_grouper_multiple_mi_groups() {
        let mut grouper = MiGrouper::new("MI", 10);

        let records = vec![
            make_raw_decoded_record("MI", "0"),
            make_raw_decoded_record("MI", "0"),
            make_raw_decoded_record("MI", "1"),
            make_raw_decoded_record("MI", "1"),
            make_raw_decoded_record("MI", "1"),
            make_raw_decoded_record("MI", "2"),
        ];

        let batches = grouper.add_records(records).expect("add_records should succeed");
        assert!(batches.is_empty()); // 3 groups < batch_size 10

        let final_batch =
            grouper.finish().expect("finish should succeed").expect("should return final batch");
        assert_eq!(final_batch.groups.len(), 3);
        assert_eq!(final_batch.groups[0].mi, "0");
        assert_eq!(final_batch.groups[0].records.len(), 2);
        assert_eq!(final_batch.groups[1].mi, "1");
        assert_eq!(final_batch.groups[1].records.len(), 3);
        assert_eq!(final_batch.groups[2].mi, "2");
        assert_eq!(final_batch.groups[2].records.len(), 1);
    }

    #[test]
    fn test_grouper_batch_size_triggers() {
        let mut grouper = MiGrouper::new("MI", 2);

        let records = vec![
            make_raw_decoded_record("MI", "0"),
            make_raw_decoded_record("MI", "1"),
            make_raw_decoded_record("MI", "2"),
            make_raw_decoded_record("MI", "3"),
            make_raw_decoded_record("MI", "4"),
        ];

        let batches = grouper.add_records(records).expect("add_records should succeed");
        // 5 MI values => 4 completed groups (0,1,2,3) while "4" is still current
        // batch_size=2 => 2 batches of 2 groups each
        assert_eq!(batches.len(), 2);
        assert_eq!(batches[0].groups.len(), 2);
        assert_eq!(batches[0].groups[0].mi, "0");
        assert_eq!(batches[0].groups[1].mi, "1");
        assert_eq!(batches[1].groups.len(), 2);
        assert_eq!(batches[1].groups[0].mi, "2");
        assert_eq!(batches[1].groups[1].mi, "3");

        // "4" is still pending
        assert!(grouper.has_pending());

        let final_batch =
            grouper.finish().expect("finish should succeed").expect("should return final batch");
        assert_eq!(final_batch.groups.len(), 1);
        assert_eq!(final_batch.groups[0].mi, "4");
    }

    #[test]
    fn test_grouper_skips_records_without_mi_tag() {
        // Matches `MiGroupIterator`: records without an MI tag are skipped so
        // the consensus caller never sees them (it requires the MI tag).
        let mut grouper = MiGrouper::new("MI", 10);

        let records = vec![
            make_raw_decoded_record("MI", "0"),
            make_raw_decoded_record_no_tag(),
            make_raw_decoded_record("MI", "0"),
        ];

        let batches = grouper.add_records(records).expect("add_records should succeed");
        assert!(batches.is_empty());

        let final_batch =
            grouper.finish().expect("finish should succeed").expect("should return final batch");
        // The no-MI record is dropped; the two MI=0 records coalesce into one group.
        assert_eq!(final_batch.groups.len(), 1);
        assert_eq!(final_batch.groups[0].mi, "0");
        assert_eq!(final_batch.groups[0].records.len(), 2);
    }

    #[test]
    fn test_grouper_finish_empty() {
        let mut grouper = MiGrouper::new("MI", 10);
        assert!(!grouper.has_pending());

        let final_batch = grouper.finish().expect("finish should succeed");
        assert!(final_batch.is_none());
    }

    #[test]
    fn test_grouper_finish_idempotent() {
        let mut grouper = MiGrouper::new("MI", 10);

        let records = vec![make_raw_decoded_record("MI", "0")];
        grouper.add_records(records).expect("add_records should succeed");

        let batch1 = grouper.finish().expect("finish should succeed");
        assert!(batch1.is_some());

        let batch2 = grouper.finish().expect("finish should succeed");
        assert!(batch2.is_none());
    }

    #[test]
    fn test_grouper_has_pending_states() {
        let mut grouper = MiGrouper::new("MI", 10);

        // Initially nothing pending
        assert!(!grouper.has_pending());

        // After adding records, current_mi is set
        let records = vec![make_raw_decoded_record("MI", "0")];
        grouper.add_records(records).expect("add_records should succeed");
        assert!(grouper.has_pending());

        // After adding a different MI, first group moves to pending_groups
        let records = vec![make_raw_decoded_record("MI", "1")];
        grouper.add_records(records).expect("add_records should succeed");
        assert!(grouper.has_pending());
    }

    #[test]
    fn test_grouper_with_filter_and_transform() {
        // Build records with flag field set: use byte 14..16 for flag.
        // We create a filter that checks the flag field for secondary alignment (0x100).
        let mut grouper = MiGrouper::with_filter_and_transform(
            "MI",
            10,
            |bam: &[u8]| {
                // Check flag field
                let flag = fgumi_raw_bam::RawRecordView::new(bam).flags();
                flag & fgumi_raw_bam::flags::SECONDARY == 0 // keep if NOT secondary
            },
            |raw: &[u8]| {
                let s = String::from_utf8_lossy(raw);
                extract_mi_base(&s).to_string()
            },
        );

        // Create records: 1/A (primary), 1/A (secondary flag), 1/B (primary)
        let rec_primary = make_raw_bam_with_tag("MI", "1/A");
        // rec_primary flag is 0 (primary) - default

        let mut rec_secondary = make_raw_bam_with_tag("MI", "1/A");
        // Set secondary flag (0x100) at bytes 14..16
        rec_secondary[14..16].copy_from_slice(&0x100u16.to_le_bytes());

        let rec_b = make_raw_bam_with_tag("MI", "1/B");

        let key = crate::unified_pipeline::GroupKey::single(0, 0, 0, 0, 0, 0);
        let records = vec![
            DecodedRecord::from_raw_bytes(rec_primary, key),
            DecodedRecord::from_raw_bytes(rec_secondary, key),
            DecodedRecord::from_raw_bytes(rec_b, key),
        ];

        let batches = grouper.add_records(records).expect("add_records should succeed");
        assert!(batches.is_empty());

        let final_batch =
            grouper.finish().expect("finish should succeed").expect("should return final batch");
        // Secondary was filtered out, and 1/A and 1/B transform to "1"
        assert_eq!(final_batch.groups.len(), 1);
        assert_eq!(final_batch.groups[0].mi, "1");
        assert_eq!(final_batch.groups[0].records.len(), 2); // primary 1/A + 1/B
    }

    #[test]
    fn test_grouper_cell_tag_composite_grouping() {
        // Records with the same MI but different cell barcodes must form separate groups.
        // The composite key (MI\tCB) is used internally for grouping only; it is stored in
        // `MiGroup::mi` for logging but never written back to the output records.
        fn make_two_tag_decoded(mi: &str, cb: &str) -> DecodedRecord {
            let raw = make_raw_bam_with_two_tags("MI", mi, "CB", cb);
            let key = crate::unified_pipeline::GroupKey::single(0, 0, 0, 0, 0, 0);
            DecodedRecord::from_raw_bytes(raw, key)
        }

        let mut grouper = MiGrouper::new("MI", 10).with_cell_tag(Some(*SamTag::CB));

        // Two records per (MI, CB) combination; MI=1,CB=ACGT and MI=1,CB=TGCA must be split.
        let records = vec![
            make_two_tag_decoded("1", "ACGT"),
            make_two_tag_decoded("1", "ACGT"),
            make_two_tag_decoded("1", "TGCA"),
            make_two_tag_decoded("1", "TGCA"),
        ];

        let batches = grouper.add_records(records).expect("add_records should succeed");
        assert!(batches.is_empty()); // 2 groups < batch_size 10

        let final_batch =
            grouper.finish().expect("finish should succeed").expect("should return final batch");
        assert_eq!(final_batch.groups.len(), 2);
        // Composite key is stored in the mi field
        assert_eq!(final_batch.groups[0].mi, "1\tACGT");
        assert_eq!(final_batch.groups[0].records.len(), 2);
        assert_eq!(final_batch.groups[1].mi, "1\tTGCA");
        assert_eq!(final_batch.groups[1].records.len(), 2);

        // Output records must be unchanged (composite key not written back as a tag).
        for group in &final_batch.groups {
            for raw in &group.records {
                use crate::sort::bam_fields;
                // MI tag still holds the original value, not the composite key
                let mi_val = bam_fields::find_string_tag_in_record(raw, SamTag::MI);
                assert_eq!(mi_val, Some(b"1".as_ref()), "MI tag must retain original value");
            }
        }
    }

    // ========================================================================
    // Batch type tests
    // ========================================================================

    #[test]
    fn test_mi_group_new() {
        let raw = make_raw_bam_with_tag("MI", "42");
        let group = MiGroup::new("42".to_string(), vec![raw]);
        assert_eq!(group.mi, "42");
        assert_eq!(group.records.len(), 1);
    }

    #[test]
    fn test_mi_group_batch_weight() {
        let group = MiGroup::new(
            "0".to_string(),
            vec![make_raw_bam_with_tag("MI", "0"), make_raw_bam_with_tag("MI", "0")],
        );
        assert_eq!(group.batch_weight(), 2);
    }

    #[test]
    fn test_mi_group_memory_estimate() {
        let group = MiGroup::new("0".to_string(), vec![make_raw_bam_with_tag("MI", "0")]);
        let size = group.estimate_heap_size();
        assert!(size > 0);
    }

    #[test]
    fn test_mi_group_batch_new() {
        let batch = MiGroupBatch::new();
        assert!(batch.groups.is_empty());
    }

    #[test]
    fn test_mi_group_batch_default() {
        let batch = MiGroupBatch::default();
        assert!(batch.groups.is_empty());
    }

    #[test]
    fn test_mi_group_batch_weight_method() {
        let mut batch = MiGroupBatch::new();
        assert_eq!(batch.batch_weight(), 0);

        batch.groups.push(MiGroup::new(
            "0".to_string(),
            vec![make_raw_bam_with_tag("MI", "0"), make_raw_bam_with_tag("MI", "0")],
        ));
        assert_eq!(batch.batch_weight(), 2);

        batch.groups.push(MiGroup::new("1".to_string(), vec![make_raw_bam_with_tag("MI", "1")]));
        assert_eq!(batch.batch_weight(), 3);
    }

    #[test]
    fn test_mi_group_batch_memory_estimate() {
        let batch = MiGroupBatch::new();
        let _ = batch.estimate_heap_size(); // Should not panic

        let mut batch = MiGroupBatch::new();
        batch.groups.push(MiGroup::new("0".to_string(), vec![make_raw_bam_with_tag("MI", "0")]));
        let size = batch.estimate_heap_size();
        assert!(size > 0);
    }

    #[test]
    fn test_mi_group_batch_with_capacity() {
        let batch = MiGroupBatch::with_capacity(16);
        assert!(batch.is_empty());
        assert_eq!(batch.len(), 0);
    }

    #[test]
    fn test_mi_group_batch_len_and_is_empty() {
        let mut batch = MiGroupBatch::new();
        assert!(batch.is_empty());

        batch.groups.push(MiGroup::new("0".to_string(), vec![make_raw_bam_with_tag("MI", "0")]));
        assert!(!batch.is_empty());
        assert_eq!(batch.len(), 1);

        batch.groups.push(MiGroup::new(
            "1".to_string(),
            vec![make_raw_bam_with_tag("MI", "1"), make_raw_bam_with_tag("MI", "1")],
        ));
        assert_eq!(batch.len(), 2);
    }

    #[test]
    fn test_mi_group_batch_clear() {
        let mut batch = MiGroupBatch::new();
        batch.groups.push(MiGroup::new("0".to_string(), vec![make_raw_bam_with_tag("MI", "0")]));
        assert!(!batch.is_empty());

        batch.clear();
        assert!(batch.is_empty());
        assert_eq!(batch.len(), 0);
    }

    #[test]
    fn test_mi_group_batch_weight_single() {
        let group = MiGroup::new(
            "0".to_string(),
            vec![
                make_raw_bam_with_tag("MI", "0"),
                make_raw_bam_with_tag("MI", "0"),
                make_raw_bam_with_tag("MI", "0"),
            ],
        );
        assert_eq!(group.batch_weight(), 3);
    }

    // ========================================================================
    // MiGroupIterator::with_transform edge case tests
    // ========================================================================

    #[test]
    fn test_with_transform_single_group() {
        // 0/A, 0/B, 0/A all transform to "0" → one group of 3
        let records: Vec<Result<RawRecord>> = vec![
            Ok(make_raw_bam_with_tag("MI", "0/A")),
            Ok(make_raw_bam_with_tag("MI", "0/B")),
            Ok(make_raw_bam_with_tag("MI", "0/A")),
        ];
        let mut iter = MiGroupIterator::with_transform(records.into_iter(), "MI", |raw| {
            let s = String::from_utf8_lossy(raw);
            extract_mi_base(&s).to_string()
        });

        let result = iter.next().expect("iterator should yield item").expect("item should be Ok");
        assert_eq!(result.0, "0");
        assert_eq!(result.1.len(), 3);

        assert!(iter.next().is_none());
    }

    #[test]
    fn test_with_transform_error_propagation() {
        let records: Vec<Result<RawRecord>> = vec![
            Ok(make_raw_bam_with_tag("MI", "0/A")),
            Err(anyhow::anyhow!("test error")),
            Ok(make_raw_bam_with_tag("MI", "1/B")),
        ];
        let mut iter = MiGroupIterator::with_transform(records.into_iter(), "MI", |raw| {
            let s = String::from_utf8_lossy(raw);
            extract_mi_base(&s).to_string()
        });

        // First group before the error
        let result = iter.next().expect("iterator should yield item").expect("item should be Ok");
        assert_eq!(result.0, "0");
        assert_eq!(result.1.len(), 1);

        // Error is returned next
        let result = iter.next().expect("iterator should yield item");
        assert!(result.is_err());

        // Iterator is done after the error
        assert!(iter.next().is_none());
    }

    #[test]
    fn test_with_transform_custom_function() {
        // Custom transform: uppercase. "abc", "ABC", "Abc" all → "ABC" → single group.
        let records: Vec<Result<RawRecord>> = vec![
            Ok(make_raw_bam_with_tag("MI", "abc")),
            Ok(make_raw_bam_with_tag("MI", "ABC")),
            Ok(make_raw_bam_with_tag("MI", "Abc")),
        ];
        let mut iter = MiGroupIterator::with_transform(records.into_iter(), "MI", |raw| {
            String::from_utf8_lossy(raw).to_uppercase()
        });

        let result = iter.next().expect("iterator should yield item").expect("item should be Ok");
        assert_eq!(result.0, "ABC");
        assert_eq!(result.1.len(), 3);

        assert!(iter.next().is_none());
    }
}

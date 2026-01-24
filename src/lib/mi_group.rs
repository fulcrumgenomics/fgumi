//! Molecular Identifier (MI) group utilities for streaming BAM processing.
//!
//! This module provides an iterator for grouping consecutive BAM records by their MI tag.
//! This is useful for streaming consensus calling where input is already sorted by MI groups
//! (e.g., output from `group`).

use anyhow::{Result, bail};
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::RecordBuf;
use std::io;

use crate::unified_pipeline::{BatchWeight, DecodedRecord, MemoryEstimate};

/// An iterator that groups consecutive BAM records by their MI (Molecular Identifier) tag.
///
/// This iterator assumes records with the same MI value are consecutive in the input stream,
/// which is the output format from `group`. Each call to `next()` returns all records
/// sharing the same MI value as a single group.
///
/// # Example
///
/// ```ignore
/// use fgumi_lib::mi_group::MiGroupIterator;
///
/// let record_iter = reader.record_bufs(&header).map(|r| r.map_err(Into::into));
/// let mi_groups = MiGroupIterator::new(record_iter, "MI");
///
/// for result in mi_groups {
///     let (mi_value, records) = result?;
///     // Process all records with this MI value
/// }
/// ```
pub struct MiGroupIterator<I>
where
    I: Iterator<Item = Result<RecordBuf>>,
{
    record_iter: I,
    tag: Tag,
    current_mi: Option<String>,
    current_group: Vec<RecordBuf>,
    pending_error: Option<anyhow::Error>,
    done: bool,
}

impl<I> MiGroupIterator<I>
where
    I: Iterator<Item = Result<RecordBuf>>,
{
    /// Creates a new `MiGroupIterator`.
    ///
    /// # Arguments
    ///
    /// * `record_iter` - Iterator over records to group by MI tag
    /// * `tag_name` - The two-character tag name (e.g., "MI")
    ///
    /// # Panics
    ///
    /// Panics if `tag_name` is not exactly 2 characters.
    pub fn new(record_iter: I, tag_name: &str) -> Self {
        assert!(tag_name.len() == 2, "Tag name must be exactly 2 characters");
        let tag_bytes = tag_name.as_bytes();
        let tag = Tag::from([tag_bytes[0], tag_bytes[1]]);

        MiGroupIterator {
            record_iter,
            tag,
            current_mi: None,
            current_group: Vec::new(),
            pending_error: None,
            done: false,
        }
    }

    /// Extracts the MI tag value from a record.
    fn get_mi(&self, record: &RecordBuf) -> Result<Option<String>> {
        if let Some(tag_value) = record.data().get(&self.tag) {
            match tag_value {
                noodles::sam::alignment::record_buf::data::field::Value::String(s) => {
                    Ok(Some(s.to_string()))
                }
                _ => {
                    bail!("MI tag must be a string value");
                }
            }
        } else {
            Ok(None)
        }
    }
}

impl<I> Iterator for MiGroupIterator<I>
where
    I: Iterator<Item = Result<RecordBuf>>,
{
    /// Each item is a Result containing (`MI_value`, `Vec<RecordBuf>`) or an error.
    type Item = Result<(String, Vec<RecordBuf>)>;

    /// Returns the next MI group from the record iterator.
    ///
    /// Reads records until the MI tag changes, then returns the completed group.
    /// On EOF, returns any remaining group, then None.
    ///
    /// Records without an MI tag are skipped.
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
                    // EOF - return any remaining group
                    self.done = true;
                    if self.current_group.is_empty() {
                        return None;
                    }
                    let mi = self.current_mi.take().unwrap_or_default();
                    let group = std::mem::take(&mut self.current_group);
                    return Some(Ok((mi, group)));
                }
                Some(Err(e)) => {
                    // If we have a pending group, return it first and save the error
                    if !self.current_group.is_empty() {
                        self.pending_error = Some(e);
                        let mi = self.current_mi.take().unwrap_or_default();
                        let group = std::mem::take(&mut self.current_group);
                        return Some(Ok((mi, group)));
                    }
                    self.done = true;
                    return Some(Err(e));
                }
                Some(Ok(record)) => {
                    // Extract MI tag
                    let mi = match self.get_mi(&record) {
                        Ok(Some(mi)) => mi,
                        Ok(None) => {
                            // Skip records without MI tag
                            continue;
                        }
                        Err(e) => {
                            // If we have a pending group, return it first and save the error
                            if !self.current_group.is_empty() {
                                self.pending_error = Some(e);
                                let mi = self.current_mi.take().unwrap_or_default();
                                let group = std::mem::take(&mut self.current_group);
                                return Some(Ok((mi, group)));
                            }
                            self.done = true;
                            return Some(Err(e));
                        }
                    };

                    if self.current_group.is_empty() {
                        // First record or first record after returning a group
                        self.current_mi = Some(mi);
                        self.current_group.push(record);
                    } else if self.current_mi.as_ref() == Some(&mi) {
                        // Same MI, add to current group
                        self.current_group.push(record);
                    } else {
                        // Different MI - return current group and start new one
                        let old_mi = self.current_mi.take().unwrap_or_default();
                        let group = std::mem::take(&mut self.current_group);
                        self.current_mi = Some(mi);
                        self.current_group.push(record);
                        return Some(Ok((old_mi, group)));
                    }
                }
            }
        }
    }
}

/// An iterator that groups consecutive BAM records by a transformed key.
///
/// Similar to `MiGroupIterator` but applies a key transformation function to the MI tag
/// before grouping. This is useful for duplex consensus calling where reads with
/// "1/A" and "1/B" should be grouped together under key "1".
///
/// # Example
///
/// ```ignore
/// use fgumi_lib::mi_group::MiGroupIteratorWithTransform;
/// use fgumi_lib::umi::extract_mi_base;
///
/// let record_iter = reader.record_bufs(&header).map(|r| r.map_err(Into::into));
/// let mi_groups = MiGroupIteratorWithTransform::new(
///     record_iter,
///     "MI",
///     |mi| extract_mi_base(mi).to_string(),
/// );
///
/// for result in mi_groups {
///     let (base_mi, records) = result?;
///     // records contains all reads with MI tags like "base_mi/A" and "base_mi/B"
/// }
/// ```
pub struct MiGroupIteratorWithTransform<I, F>
where
    I: Iterator<Item = Result<RecordBuf>>,
    F: Fn(&str) -> String,
{
    record_iter: I,
    tag: Tag,
    key_transform: F,
    current_key: Option<String>,
    current_group: Vec<RecordBuf>,
    pending_error: Option<anyhow::Error>,
    done: bool,
}

impl<I, F> MiGroupIteratorWithTransform<I, F>
where
    I: Iterator<Item = Result<RecordBuf>>,
    F: Fn(&str) -> String,
{
    /// Creates a new `MiGroupIteratorWithTransform`.
    ///
    /// # Arguments
    ///
    /// * `record_iter` - Iterator over records to group
    /// * `tag_name` - The two-character tag name (e.g., "MI")
    /// * `key_transform` - Function to transform the tag value into a grouping key
    ///
    /// # Panics
    ///
    /// Panics if `tag_name` is not exactly 2 characters.
    pub fn new(record_iter: I, tag_name: &str, key_transform: F) -> Self {
        assert!(tag_name.len() == 2, "Tag name must be exactly 2 characters");
        let tag_bytes = tag_name.as_bytes();
        let tag = Tag::from([tag_bytes[0], tag_bytes[1]]);

        MiGroupIteratorWithTransform {
            record_iter,
            tag,
            key_transform,
            current_key: None,
            current_group: Vec::new(),
            pending_error: None,
            done: false,
        }
    }

    /// Extracts the tag value from a record and transforms it using the key function.
    fn get_key(&self, record: &RecordBuf) -> Result<Option<String>> {
        if let Some(tag_value) = record.data().get(&self.tag) {
            match tag_value {
                noodles::sam::alignment::record_buf::data::field::Value::String(s) => {
                    let raw = s.to_string();
                    Ok(Some((self.key_transform)(&raw)))
                }
                _ => {
                    bail!("Tag must be a string value");
                }
            }
        } else {
            Ok(None)
        }
    }
}

impl<I, F> Iterator for MiGroupIteratorWithTransform<I, F>
where
    I: Iterator<Item = Result<RecordBuf>>,
    F: Fn(&str) -> String,
{
    /// Each item is a Result containing (key, `Vec<RecordBuf>`) or an error.
    type Item = Result<(String, Vec<RecordBuf>)>;

    /// Returns the next group from the record iterator.
    ///
    /// Reads records until the transformed key changes, then returns the completed group.
    /// On EOF, returns any remaining group, then None.
    ///
    /// Records without the tag are skipped.
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
                    // EOF - return any remaining group
                    self.done = true;
                    if self.current_group.is_empty() {
                        return None;
                    }
                    let key = self.current_key.take().unwrap_or_default();
                    let group = std::mem::take(&mut self.current_group);
                    return Some(Ok((key, group)));
                }
                Some(Err(e)) => {
                    // If we have a pending group, return it first and save the error
                    if !self.current_group.is_empty() {
                        self.pending_error = Some(e);
                        let key = self.current_key.take().unwrap_or_default();
                        let group = std::mem::take(&mut self.current_group);
                        return Some(Ok((key, group)));
                    }
                    self.done = true;
                    return Some(Err(e));
                }
                Some(Ok(record)) => {
                    // Extract and transform tag
                    let key = match self.get_key(&record) {
                        Ok(Some(key)) => key,
                        Ok(None) => {
                            // Skip records without tag
                            continue;
                        }
                        Err(e) => {
                            // If we have a pending group, return it first and save the error
                            if !self.current_group.is_empty() {
                                self.pending_error = Some(e);
                                let key = self.current_key.take().unwrap_or_default();
                                let group = std::mem::take(&mut self.current_group);
                                return Some(Ok((key, group)));
                            }
                            self.done = true;
                            return Some(Err(e));
                        }
                    };

                    if self.current_group.is_empty() {
                        // First record or first record after returning a group
                        self.current_key = Some(key);
                        self.current_group.push(record);
                    } else if self.current_key.as_ref() == Some(&key) {
                        // Same key, add to current group
                        self.current_group.push(record);
                    } else {
                        // Different key - return current group and start new one
                        let old_key = self.current_key.take().unwrap_or_default();
                        let group = std::mem::take(&mut self.current_group);
                        self.current_key = Some(key);
                        self.current_group.push(record);
                        return Some(Ok((old_key, group)));
                    }
                }
            }
        }
    }
}

// ============================================================================
// Parallel processing support for MI groups
// ============================================================================

/// A single MI group: the MI tag value and all records with that MI.
///
/// This represents records that share the same Molecular Identifier (MI) tag value,
/// used by consensus calling commands to group reads from the same source molecule.
#[derive(Debug, Clone)]
pub struct MiGroup {
    /// The MI tag value (e.g., "0", "1/A", "1/B")
    pub mi: String,
    /// All records sharing this MI value
    pub records: Vec<RecordBuf>,
}

impl MiGroup {
    /// Creates a new MI group.
    #[must_use]
    pub fn new(mi: String, records: Vec<RecordBuf>) -> Self {
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
        // mi: String (heap allocated)
        let mi_size = self.mi.capacity();

        // records: Vec<RecordBuf>
        let records_size: usize = self.records.iter().map(MemoryEstimate::estimate_heap_size).sum();
        let records_vec_overhead = self.records.capacity() * std::mem::size_of::<RecordBuf>();

        mi_size + records_size + records_vec_overhead
    }
}

impl MemoryEstimate for MiGroupBatch {
    fn estimate_heap_size(&self) -> usize {
        let groups_size: usize = self.groups.iter().map(MemoryEstimate::estimate_heap_size).sum();
        let groups_vec_overhead = self.groups.capacity() * std::mem::size_of::<MiGroup>();
        groups_size + groups_vec_overhead
    }
}

/// A batch of MI groups for parallel processing.
///
/// This structure holds a collection of MI groups that will be processed together
/// in parallel.
///
/// # Performance
///
/// The batch is designed to be reused across iterations to minimize allocations.
#[derive(Default)]
pub struct MiGroupBatch {
    /// The MI groups in this batch
    pub groups: Vec<MiGroup>,
}

impl MiGroupBatch {
    /// Creates a new empty MI group batch.
    #[must_use]
    pub fn new() -> Self {
        Self { groups: Vec::new() }
    }

    /// Creates a new MI group batch with pre-allocated capacity.
    #[must_use]
    pub fn with_capacity(capacity: usize) -> Self {
        Self { groups: Vec::with_capacity(capacity) }
    }

    /// Returns the number of groups in the batch.
    #[must_use]
    pub fn len(&self) -> usize {
        self.groups.len()
    }

    /// Returns true if the batch is empty.
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

// ============================================================================
// MiGrouper for 7-step unified pipeline
// ============================================================================

use crate::unified_pipeline::Grouper;
use std::collections::VecDeque;

/// Type alias for record filter function used in [`MiGrouper`].
type RecordFilterFn = Box<dyn Fn(&RecordBuf) -> bool + Send + Sync>;

/// Type alias for MI tag transformation function used in [`MiGrouper`].
type MiTransformFn = Box<dyn Fn(&str) -> String + Send + Sync>;

/// A Grouper that reads pre-grouped BAM records by MI tag for the 7-step pipeline.
///
/// Input BAM is expected to be sorted/grouped by MI tag (output from group command).
/// This grouper reads consecutive records with the same MI tag and yields them as groups.
///
/// Unlike `OwnedMiGroupReader` which implements the `Reader` trait for the parallel
/// framework, this implements the `Grouper` trait for the 7-step unified pipeline.
/// The key difference is that `Grouper::add_bytes` receives raw decompressed bytes
/// and must deserialize BAM records inline.
///
/// # Example
///
/// ```ignore
/// use fgumi_lib::mi_group::MiGrouper;
/// use fgumi_lib::unified_pipeline::Grouper;
///
/// let grouper = MiGrouper::new("MI", 100);
/// // Use with run_bam_pipeline_with_grouper...
/// ```
///
/// For duplex consensus calling with record filtering and MI transformation:
///
/// ```ignore
/// let grouper = MiGrouper::with_filter_and_transform(
///     "MI",
///     100,
///     |r| !r.flags().is_secondary(), // filter function
///     |mi| extract_mi_base(mi).to_string(),   // transform function
/// );
/// ```
pub struct MiGrouper {
    /// The MI tag to group by (e.g., "MI")
    tag: Tag,
    /// Number of MI groups per batch
    batch_size: usize,
    /// Current MI value being accumulated
    current_mi: Option<String>,
    /// Records in current MI group
    current_records: Vec<RecordBuf>,
    /// Completed groups waiting to be batched
    pending_groups: VecDeque<MiGroup>,
    /// Whether `finish()` has been called
    finished: bool,
    /// Optional record filter (returns true to keep record)
    record_filter: Option<RecordFilterFn>,
    /// Optional MI tag transformation function
    mi_transform: Option<MiTransformFn>,
}

impl MiGrouper {
    /// Create a new `MiGrouper`.
    ///
    /// # Arguments
    /// * `tag_name` - The MI tag name (e.g., "MI")
    /// * `batch_size` - Number of MI groups per batch (100 is typical)
    ///
    /// # Panics
    ///
    /// Panics if `tag_name` is not exactly 2 characters.
    #[must_use]
    pub fn new(tag_name: &str, batch_size: usize) -> Self {
        assert!(tag_name.len() == 2, "Tag name must be exactly 2 characters");
        let tag_bytes = tag_name.as_bytes();
        let tag = Tag::from([tag_bytes[0], tag_bytes[1]]);

        Self {
            tag,
            batch_size: batch_size.max(1),
            current_mi: None,
            current_records: Vec::new(),
            pending_groups: VecDeque::new(),
            finished: false,
            record_filter: None,
            mi_transform: None,
        }
    }

    /// Create a `MiGrouper` with record filtering and MI tag transformation.
    ///
    /// This is useful for duplex consensus calling where we need to:
    /// - Filter out secondary/supplementary reads
    /// - Transform MI tags by stripping /A and /B suffixes
    ///
    /// # Arguments
    /// * `tag_name` - The MI tag name (e.g., "MI")
    /// * `batch_size` - Number of MI groups per batch (100 is typical)
    /// * `record_filter` - Function that returns true to keep a record
    /// * `mi_transform` - Function to transform MI tag value (e.g., strip /A /B suffix)
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
        F: Fn(&RecordBuf) -> bool + Send + Sync + 'static,
        T: Fn(&str) -> String + Send + Sync + 'static,
    {
        assert!(tag_name.len() == 2, "Tag name must be exactly 2 characters");
        let tag_bytes = tag_name.as_bytes();
        let tag = Tag::from([tag_bytes[0], tag_bytes[1]]);

        Self {
            tag,
            batch_size: batch_size.max(1),
            current_mi: None,
            current_records: Vec::new(),
            pending_groups: VecDeque::new(),
            finished: false,
            record_filter: Some(Box::new(record_filter)),
            mi_transform: Some(Box::new(mi_transform)),
        }
    }

    /// Get the MI tag value from a record, optionally applying transformation.
    fn get_mi_tag(&self, record: &RecordBuf) -> Option<String> {
        record.data().get(&self.tag).and_then(|v| match v {
            noodles::sam::alignment::record_buf::data::field::Value::String(s) => {
                let raw = s.to_string();
                // Apply transform if configured
                if let Some(ref transform) = self.mi_transform {
                    Some(transform(&raw))
                } else {
                    Some(raw)
                }
            }
            _ => None,
        })
    }

    /// Check if a record passes the filter (or return true if no filter).
    fn should_keep(&self, record: &RecordBuf) -> bool {
        match &self.record_filter {
            Some(filter) => filter(record),
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
            let record = decoded.record;
            // Apply record filter if configured - skip records that don't pass
            if !self.should_keep(&record) {
                continue;
            }

            // Get MI tag from record (with optional transformation)
            let mi = self.get_mi_tag(&record).unwrap_or_default();

            // Check if this starts a new MI group
            match &self.current_mi {
                Some(current) if current == &mi => {
                    // Same MI, add to current group
                    self.current_records.push(record);
                }
                Some(_) => {
                    // Different MI, flush current and start new
                    self.flush_current_group();
                    self.current_mi = Some(mi);
                    self.current_records.push(record);
                }
                None => {
                    // First record
                    self.current_mi = Some(mi);
                    self.current_records.push(record);
                }
            }
        }

        // Return complete batches
        Ok(self.drain_batches())
    }

    fn finish(&mut self) -> io::Result<Option<Self::Group>> {
        if self.finished {
            return Ok(None);
        }
        self.finished = true;

        // Flush any remaining current group
        self.flush_current_group();

        // Return any remaining groups as final batch
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

#[cfg(test)]
#[allow(clippy::similar_names)]
mod tests {
    use super::*;
    use crate::sam::builder::RecordBuilder;
    use crate::umi::extract_mi_base;

    fn create_record_with_mi(mi: &str) -> RecordBuf {
        RecordBuilder::new()
            .sequence("ACGT") // Minimal sequence
            .tag("MI", mi)
            .build()
    }

    fn create_record_without_mi() -> RecordBuf {
        RecordBuilder::new()
            .sequence("ACGT") // Minimal sequence
            .build()
    }

    #[test]
    fn test_empty_iterator() {
        let records: Vec<Result<RecordBuf>> = vec![];
        let mut iter = MiGroupIterator::new(records.into_iter(), "MI");
        assert!(iter.next().is_none());
    }

    #[test]
    fn test_single_group() {
        let records: Vec<Result<RecordBuf>> = vec![
            Ok(create_record_with_mi("0")),
            Ok(create_record_with_mi("0")),
            Ok(create_record_with_mi("0")),
        ];
        let mut iter = MiGroupIterator::new(records.into_iter(), "MI");

        let result = iter.next().unwrap().unwrap();
        assert_eq!(result.0, "0");
        assert_eq!(result.1.len(), 3);

        assert!(iter.next().is_none());
    }

    #[test]
    fn test_multiple_groups() {
        let records: Vec<Result<RecordBuf>> = vec![
            Ok(create_record_with_mi("0")),
            Ok(create_record_with_mi("0")),
            Ok(create_record_with_mi("1")),
            Ok(create_record_with_mi("1")),
            Ok(create_record_with_mi("1")),
            Ok(create_record_with_mi("2")),
        ];
        let mut iter = MiGroupIterator::new(records.into_iter(), "MI");

        let result = iter.next().unwrap().unwrap();
        assert_eq!(result.0, "0");
        assert_eq!(result.1.len(), 2);

        let result = iter.next().unwrap().unwrap();
        assert_eq!(result.0, "1");
        assert_eq!(result.1.len(), 3);

        let result = iter.next().unwrap().unwrap();
        assert_eq!(result.0, "2");
        assert_eq!(result.1.len(), 1);

        assert!(iter.next().is_none());
    }

    #[test]
    fn test_skips_records_without_mi_tag() {
        let records: Vec<Result<RecordBuf>> = vec![
            Ok(create_record_with_mi("0")),
            Ok(create_record_without_mi()),
            Ok(create_record_with_mi("0")),
            Ok(create_record_without_mi()),
            Ok(create_record_with_mi("1")),
        ];
        let mut iter = MiGroupIterator::new(records.into_iter(), "MI");

        let result = iter.next().unwrap().unwrap();
        assert_eq!(result.0, "0");
        assert_eq!(result.1.len(), 2); // Skipped record without MI

        let result = iter.next().unwrap().unwrap();
        assert_eq!(result.0, "1");
        assert_eq!(result.1.len(), 1);

        assert!(iter.next().is_none());
    }

    #[test]
    fn test_error_propagation() {
        let records: Vec<Result<RecordBuf>> = vec![
            Ok(create_record_with_mi("0")),
            Err(anyhow::anyhow!("test error")),
            Ok(create_record_with_mi("1")),
        ];
        let mut iter = MiGroupIterator::new(records.into_iter(), "MI");

        // First group before error
        let result = iter.next().unwrap().unwrap();
        assert_eq!(result.0, "0");
        assert_eq!(result.1.len(), 1);

        // Error
        let result = iter.next().unwrap();
        assert!(result.is_err());

        // Iterator should be done after error
        assert!(iter.next().is_none());
    }

    #[test]
    fn test_custom_tag() {
        let record1 = RecordBuilder::new().sequence("ACGT").tag("RX", "ACGT").build();

        let record2 = RecordBuilder::new().sequence("ACGT").tag("RX", "ACGT").build();

        let records: Vec<Result<RecordBuf>> = vec![Ok(record1), Ok(record2)];
        let mut iter = MiGroupIterator::new(records.into_iter(), "RX");

        let result = iter.next().unwrap().unwrap();
        assert_eq!(result.0, "ACGT");
        assert_eq!(result.1.len(), 2);

        assert!(iter.next().is_none());
    }

    #[test]
    #[should_panic(expected = "Tag name must be exactly 2 characters")]
    fn test_invalid_tag_length() {
        let records: Vec<Result<RecordBuf>> = vec![];
        let _ = MiGroupIterator::new(records.into_iter(), "M");
    }

    // Tests for MiGroupIteratorWithTransform
    #[test]
    fn test_transform_groups_by_base_mi() {
        // Simulate duplex reads: 1/A, 1/A, 1/B, 1/B, 2/A, 2/B
        let records: Vec<Result<RecordBuf>> = vec![
            Ok(create_record_with_mi("1/A")),
            Ok(create_record_with_mi("1/A")),
            Ok(create_record_with_mi("1/B")),
            Ok(create_record_with_mi("1/B")),
            Ok(create_record_with_mi("2/A")),
            Ok(create_record_with_mi("2/B")),
        ];
        let mut iter = MiGroupIteratorWithTransform::new(records.into_iter(), "MI", |mi| {
            extract_mi_base(mi).to_string()
        });

        // First group: base MI "1" with 4 reads (2 /A + 2 /B)
        let result = iter.next().unwrap().unwrap();
        assert_eq!(result.0, "1");
        assert_eq!(result.1.len(), 4);

        // Second group: base MI "2" with 2 reads (1 /A + 1 /B)
        let result = iter.next().unwrap().unwrap();
        assert_eq!(result.0, "2");
        assert_eq!(result.1.len(), 2);

        assert!(iter.next().is_none());
    }

    #[test]
    fn test_transform_empty_iterator() {
        let records: Vec<Result<RecordBuf>> = vec![];
        let mut iter = MiGroupIteratorWithTransform::new(records.into_iter(), "MI", |mi| {
            extract_mi_base(mi).to_string()
        });
        assert!(iter.next().is_none());
    }

    #[test]
    fn test_transform_single_group() {
        let records: Vec<Result<RecordBuf>> = vec![
            Ok(create_record_with_mi("0/A")),
            Ok(create_record_with_mi("0/B")),
            Ok(create_record_with_mi("0/A")),
        ];
        let mut iter = MiGroupIteratorWithTransform::new(records.into_iter(), "MI", |mi| {
            extract_mi_base(mi).to_string()
        });

        let result = iter.next().unwrap().unwrap();
        assert_eq!(result.0, "0");
        assert_eq!(result.1.len(), 3);

        assert!(iter.next().is_none());
    }

    #[test]
    fn test_transform_error_propagation() {
        let records: Vec<Result<RecordBuf>> = vec![
            Ok(create_record_with_mi("0/A")),
            Err(anyhow::anyhow!("test error")),
            Ok(create_record_with_mi("1/B")),
        ];
        let mut iter = MiGroupIteratorWithTransform::new(records.into_iter(), "MI", |mi| {
            extract_mi_base(mi).to_string()
        });

        // First group before error
        let result = iter.next().unwrap().unwrap();
        assert_eq!(result.0, "0");
        assert_eq!(result.1.len(), 1);

        // Error
        let result = iter.next().unwrap();
        assert!(result.is_err());

        // Iterator should be done after error
        assert!(iter.next().is_none());
    }

    #[test]
    fn test_transform_custom_function() {
        // Test with a custom transformation function (uppercase)
        let records: Vec<Result<RecordBuf>> = vec![
            Ok(create_record_with_mi("abc")),
            Ok(create_record_with_mi("ABC")),
            Ok(create_record_with_mi("Abc")),
        ];
        let mut iter =
            MiGroupIteratorWithTransform::new(records.into_iter(), "MI", str::to_uppercase);

        // All should be grouped together since they uppercase to "ABC"
        let result = iter.next().unwrap().unwrap();
        assert_eq!(result.0, "ABC");
        assert_eq!(result.1.len(), 3);

        assert!(iter.next().is_none());
    }
}

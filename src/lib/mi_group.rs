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
            let record = decoded.into_record().ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "MiGrouper requires parsed records")
            })?;
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

// ============================================================================
// Raw-byte MI grouping for consensus callers
// ============================================================================

/// A single MI group holding raw-byte BAM records.
#[derive(Debug, Clone)]
pub struct RawMiGroup {
    /// The MI tag value (e.g., "0", "1/A", "1/B")
    pub mi: String,
    /// Raw BAM records sharing this MI value
    pub records: Vec<Vec<u8>>,
}

impl RawMiGroup {
    /// Creates a new raw MI group.
    #[must_use]
    pub fn new(mi: String, records: Vec<Vec<u8>>) -> Self {
        Self { mi, records }
    }
}

impl BatchWeight for RawMiGroup {
    fn batch_weight(&self) -> usize {
        self.records.len()
    }
}

impl MemoryEstimate for RawMiGroup {
    fn estimate_heap_size(&self) -> usize {
        let mi_size = self.mi.capacity();
        let records_size: usize = self.records.iter().map(|r| r.capacity()).sum();
        let records_vec_overhead = self.records.capacity() * std::mem::size_of::<Vec<u8>>();
        mi_size + records_size + records_vec_overhead
    }
}

/// A batch of raw MI groups for parallel processing.
#[derive(Default)]
pub struct RawMiGroupBatch {
    /// The raw MI groups in this batch
    pub groups: Vec<RawMiGroup>,
}

impl RawMiGroupBatch {
    /// Creates a new empty raw MI group batch.
    #[must_use]
    pub fn new() -> Self {
        Self { groups: Vec::new() }
    }
}

impl BatchWeight for RawMiGroupBatch {
    fn batch_weight(&self) -> usize {
        self.groups.iter().map(|g| g.records.len()).sum()
    }
}

impl MemoryEstimate for RawMiGroupBatch {
    fn estimate_heap_size(&self) -> usize {
        let groups_size: usize = self.groups.iter().map(MemoryEstimate::estimate_heap_size).sum();
        let groups_vec_overhead = self.groups.capacity() * std::mem::size_of::<RawMiGroup>();
        groups_size + groups_vec_overhead
    }
}

/// Type alias for raw-byte MI tag transformation function.
type RawMiTransformFn = Box<dyn Fn(&[u8]) -> String + Send + Sync>;

/// Type alias for raw-byte record filter function.
type RawRecordFilterFn = Box<dyn Fn(&[u8]) -> bool + Send + Sync>;

/// A Grouper that groups raw-byte BAM records by MI tag.
///
/// This is the raw-byte equivalent of `MiGrouper`. Instead of parsing records
/// into `RecordBuf`, it extracts MI tags directly from raw BAM bytes using
/// `bam_fields::find_string_tag_in_record()`.
pub struct RawMiGrouper {
    /// The MI tag bytes to search for
    tag: [u8; 2],
    /// Number of MI groups per batch
    batch_size: usize,
    /// Current MI value being accumulated
    current_mi: Option<String>,
    /// Records in current MI group (raw bytes)
    current_records: Vec<Vec<u8>>,
    /// Completed groups waiting to be batched
    pending_groups: VecDeque<RawMiGroup>,
    /// Whether `finish()` has been called
    finished: bool,
    /// Optional MI tag transformation function
    mi_transform: Option<RawMiTransformFn>,
    /// Optional record filter
    record_filter: Option<RawRecordFilterFn>,
}

impl RawMiGrouper {
    /// Create a new `RawMiGrouper`.
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
        }
    }

    /// Create a `RawMiGrouper` with record filtering and MI tag transformation.
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
        }
    }

    /// Get MI tag from raw BAM bytes, optionally applying transformation.
    fn get_mi_tag(&self, bam: &[u8]) -> Option<String> {
        use crate::sort::bam_fields;
        let value = bam_fields::find_string_tag_in_record(bam, &self.tag)?;
        if let Some(ref transform) = self.mi_transform {
            Some(transform(value))
        } else {
            Some(String::from_utf8_lossy(value).into_owned())
        }
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
                self.pending_groups.push_back(RawMiGroup::new(mi, records));
            }
        }
    }

    /// Try to form complete batches from pending groups.
    fn drain_batches(&mut self) -> Vec<RawMiGroupBatch> {
        let mut batches = Vec::new();
        while self.pending_groups.len() >= self.batch_size {
            let groups: Vec<RawMiGroup> = self.pending_groups.drain(..self.batch_size).collect();
            batches.push(RawMiGroupBatch { groups });
        }
        batches
    }
}

impl Grouper for RawMiGrouper {
    type Group = RawMiGroupBatch;

    fn add_records(&mut self, records: Vec<DecodedRecord>) -> io::Result<Vec<Self::Group>> {
        for decoded in records {
            let raw = decoded.into_raw_bytes().ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "RawMiGrouper requires raw byte records")
            })?;

            // Apply record filter if configured
            if !self.should_keep(&raw) {
                continue;
            }

            // Get MI tag from raw bytes (use empty string if absent, matching MiGrouper behavior)
            let mi = self.get_mi_tag(&raw).unwrap_or_default();

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
            let groups: Vec<RawMiGroup> = self.pending_groups.drain(..).collect();
            Ok(Some(RawMiGroupBatch { groups }))
        }
    }

    fn has_pending(&self) -> bool {
        !self.pending_groups.is_empty() || self.current_mi.is_some()
    }
}

/// An iterator that groups consecutive raw BAM records by their MI tag.
///
/// This is the raw-byte equivalent of `MiGroupIterator` for the single-threaded path.
#[allow(clippy::type_complexity)]
pub struct RawMiGroupIterator<I>
where
    I: Iterator<Item = Result<Vec<u8>>>,
{
    record_iter: I,
    tag: [u8; 2],
    /// Optional cell barcode tag for composite grouping (MI + cell barcode)
    cell_tag: Option<[u8; 2]>,
    current_mi: Option<String>,
    current_group: Vec<Vec<u8>>,
    done: bool,
    /// Pending error to return after flushing a group
    pending_error: Option<anyhow::Error>,
    /// Optional MI tag transformation function
    mi_transform: Option<Box<dyn Fn(&[u8]) -> String>>,
}

impl<I> RawMiGroupIterator<I>
where
    I: Iterator<Item = Result<Vec<u8>>>,
{
    /// Creates a new `RawMiGroupIterator`.
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

    /// Creates a new `RawMiGroupIterator` with MI tag transformation.
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
        let value = bam_fields::find_string_tag_in_record(bam, &self.tag)?;
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

impl<I> Iterator for RawMiGroupIterator<I>
where
    I: Iterator<Item = Result<Vec<u8>>>,
{
    type Item = Result<(String, Vec<Vec<u8>>)>;

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
    fn make_raw_bam_with_tag(tag_name: &str, tag_value: &str) -> Vec<u8> {
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

        buf
    }

    /// Build a minimal raw BAM record with no aux tags.
    fn make_raw_bam_without_tag() -> Vec<u8> {
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

        buf
    }

    // ========================================================================
    // RawMiGroupIterator tests
    // ========================================================================

    #[test]
    fn test_raw_empty_iterator() {
        let records: Vec<Result<Vec<u8>>> = vec![];
        let mut iter = RawMiGroupIterator::new(records.into_iter(), "MI");
        assert!(iter.next().is_none());
    }

    #[test]
    fn test_raw_single_group() {
        let records: Vec<Result<Vec<u8>>> = vec![
            Ok(make_raw_bam_with_tag("MI", "0")),
            Ok(make_raw_bam_with_tag("MI", "0")),
            Ok(make_raw_bam_with_tag("MI", "0")),
        ];
        let mut iter = RawMiGroupIterator::new(records.into_iter(), "MI");

        let result = iter.next().unwrap().unwrap();
        assert_eq!(result.0, "0");
        assert_eq!(result.1.len(), 3);

        assert!(iter.next().is_none());
    }

    #[test]
    fn test_raw_multiple_groups() {
        let records: Vec<Result<Vec<u8>>> = vec![
            Ok(make_raw_bam_with_tag("MI", "0")),
            Ok(make_raw_bam_with_tag("MI", "0")),
            Ok(make_raw_bam_with_tag("MI", "1")),
            Ok(make_raw_bam_with_tag("MI", "1")),
            Ok(make_raw_bam_with_tag("MI", "1")),
            Ok(make_raw_bam_with_tag("MI", "2")),
        ];
        let mut iter = RawMiGroupIterator::new(records.into_iter(), "MI");

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
    fn test_raw_skips_records_without_mi_tag() {
        let records: Vec<Result<Vec<u8>>> = vec![
            Ok(make_raw_bam_with_tag("MI", "0")),
            Ok(make_raw_bam_without_tag()),
            Ok(make_raw_bam_with_tag("MI", "0")),
            Ok(make_raw_bam_without_tag()),
            Ok(make_raw_bam_with_tag("MI", "1")),
        ];
        let mut iter = RawMiGroupIterator::new(records.into_iter(), "MI");

        let result = iter.next().unwrap().unwrap();
        assert_eq!(result.0, "0");
        assert_eq!(result.1.len(), 2); // Skipped records without MI

        let result = iter.next().unwrap().unwrap();
        assert_eq!(result.0, "1");
        assert_eq!(result.1.len(), 1);

        assert!(iter.next().is_none());
    }

    #[test]
    fn test_raw_error_propagation() {
        let records: Vec<Result<Vec<u8>>> = vec![
            Ok(make_raw_bam_with_tag("MI", "0")),
            Err(anyhow::anyhow!("test error")),
            Ok(make_raw_bam_with_tag("MI", "1")),
        ];
        let mut iter = RawMiGroupIterator::new(records.into_iter(), "MI");

        // First group before error
        let result = iter.next().unwrap().unwrap();
        assert_eq!(result.0, "0");
        assert_eq!(result.1.len(), 1);

        // Error should be returned after flushing the pending group
        let err = iter.next().unwrap();
        assert!(err.is_err());

        assert!(iter.next().is_none());
    }

    #[test]
    fn test_raw_error_with_no_pending_group() {
        let records: Vec<Result<Vec<u8>>> =
            vec![Err(anyhow::anyhow!("immediate error")), Ok(make_raw_bam_with_tag("MI", "0"))];
        let mut iter = RawMiGroupIterator::new(records.into_iter(), "MI");

        // Error should be returned directly
        let result = iter.next().unwrap();
        assert!(result.is_err());

        // Iterator should be done after error
        assert!(iter.next().is_none());
    }

    #[test]
    fn test_raw_custom_tag() {
        let records: Vec<Result<Vec<u8>>> =
            vec![Ok(make_raw_bam_with_tag("RX", "ACGT")), Ok(make_raw_bam_with_tag("RX", "ACGT"))];
        let mut iter = RawMiGroupIterator::new(records.into_iter(), "RX");

        let result = iter.next().unwrap().unwrap();
        assert_eq!(result.0, "ACGT");
        assert_eq!(result.1.len(), 2);

        assert!(iter.next().is_none());
    }

    #[test]
    #[should_panic(expected = "Tag name must be exactly 2 characters")]
    fn test_raw_invalid_tag_length() {
        let records: Vec<Result<Vec<u8>>> = vec![];
        let _ = RawMiGroupIterator::new(records.into_iter(), "M");
    }

    #[test]
    fn test_raw_with_transform() {
        // Simulate duplex reads: 1/A, 1/B should group under "1"
        let records: Vec<Result<Vec<u8>>> = vec![
            Ok(make_raw_bam_with_tag("MI", "1/A")),
            Ok(make_raw_bam_with_tag("MI", "1/A")),
            Ok(make_raw_bam_with_tag("MI", "1/B")),
            Ok(make_raw_bam_with_tag("MI", "1/B")),
            Ok(make_raw_bam_with_tag("MI", "2/A")),
            Ok(make_raw_bam_with_tag("MI", "2/B")),
        ];
        let mut iter = RawMiGroupIterator::with_transform(records.into_iter(), "MI", |raw| {
            let s = String::from_utf8_lossy(raw);
            extract_mi_base(&s).to_string()
        });

        // First group: base MI "1" with 4 reads
        let result = iter.next().unwrap().unwrap();
        assert_eq!(result.0, "1");
        assert_eq!(result.1.len(), 4);

        // Second group: base MI "2" with 2 reads
        let result = iter.next().unwrap().unwrap();
        assert_eq!(result.0, "2");
        assert_eq!(result.1.len(), 2);

        assert!(iter.next().is_none());
    }

    #[test]
    fn test_raw_with_transform_empty() {
        let records: Vec<Result<Vec<u8>>> = vec![];
        let mut iter = RawMiGroupIterator::with_transform(records.into_iter(), "MI", |raw| {
            String::from_utf8_lossy(raw).to_uppercase()
        });
        assert!(iter.next().is_none());
    }

    #[test]
    #[should_panic(expected = "Tag name must be exactly 2 characters")]
    fn test_raw_with_transform_invalid_tag_length() {
        let records: Vec<Result<Vec<u8>>> = vec![];
        let _ = RawMiGroupIterator::with_transform(records.into_iter(), "ABC", |raw| {
            String::from_utf8_lossy(raw).into_owned()
        });
    }

    #[test]
    fn test_raw_get_mi_without_transform() {
        let iter = RawMiGroupIterator::new(std::iter::empty::<Result<Vec<u8>>>(), "MI");
        let bam = make_raw_bam_with_tag("MI", "42");
        assert_eq!(iter.get_mi(&bam), Some("42".to_string()));
    }

    #[test]
    fn test_raw_get_mi_with_transform() {
        let iter = RawMiGroupIterator::with_transform(
            std::iter::empty::<Result<Vec<u8>>>(),
            "MI",
            |raw| {
                let s = String::from_utf8_lossy(raw);
                s.to_uppercase()
            },
        );
        let bam = make_raw_bam_with_tag("MI", "abc");
        assert_eq!(iter.get_mi(&bam), Some("ABC".to_string()));
    }

    #[test]
    fn test_raw_get_mi_missing_tag() {
        let iter = RawMiGroupIterator::new(std::iter::empty::<Result<Vec<u8>>>(), "MI");
        let bam = make_raw_bam_without_tag();
        assert_eq!(iter.get_mi(&bam), None);
    }

    #[test]
    fn test_raw_get_mi_wrong_tag() {
        let iter = RawMiGroupIterator::new(std::iter::empty::<Result<Vec<u8>>>(), "MI");
        let bam = make_raw_bam_with_tag("RX", "ACGT");
        assert_eq!(iter.get_mi(&bam), None);
    }

    /// Build a minimal raw BAM record with two string tags.
    fn make_raw_bam_with_two_tags(tag1: &str, val1: &str, tag2: &str, val2: &str) -> Vec<u8> {
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

        buf
    }

    #[test]
    fn test_raw_cell_tag_composite_grouping() {
        // Same MI but different cell barcodes should form separate groups
        let records: Vec<Result<Vec<u8>>> = vec![
            Ok(make_raw_bam_with_two_tags("MI", "1", "CB", "ACGT")),
            Ok(make_raw_bam_with_two_tags("MI", "1", "CB", "ACGT")),
            Ok(make_raw_bam_with_two_tags("MI", "1", "CB", "TGCA")),
            Ok(make_raw_bam_with_two_tags("MI", "1", "CB", "TGCA")),
        ];
        let mut iter =
            RawMiGroupIterator::new(records.into_iter(), "MI").with_cell_tag(Some([b'C', b'B']));

        // First group: MI=1, CB=ACGT
        let result = iter.next().unwrap().unwrap();
        assert_eq!(result.0, "1\tACGT");
        assert_eq!(result.1.len(), 2);

        // Second group: MI=1, CB=TGCA
        let result = iter.next().unwrap().unwrap();
        assert_eq!(result.0, "1\tTGCA");
        assert_eq!(result.1.len(), 2);

        assert!(iter.next().is_none());
    }

    #[test]
    fn test_raw_cell_tag_none_groups_by_mi_only() {
        // Without cell_tag, same MI records group together regardless of other tags
        let records: Vec<Result<Vec<u8>>> = vec![
            Ok(make_raw_bam_with_two_tags("MI", "1", "CB", "ACGT")),
            Ok(make_raw_bam_with_two_tags("MI", "1", "CB", "TGCA")),
        ];
        let mut iter = RawMiGroupIterator::new(records.into_iter(), "MI").with_cell_tag(None);

        let result = iter.next().unwrap().unwrap();
        assert_eq!(result.0, "1");
        assert_eq!(result.1.len(), 2); // Both records in same group

        assert!(iter.next().is_none());
    }

    #[test]
    fn test_raw_cell_tag_missing_cell_value() {
        // Records with MI but no cell tag get grouped under "MI\t" (empty cell)
        let records: Vec<Result<Vec<u8>>> = vec![
            Ok(make_raw_bam_with_tag("MI", "1")),
            Ok(make_raw_bam_with_tag("MI", "1")),
            Ok(make_raw_bam_with_two_tags("MI", "1", "CB", "ACGT")),
        ];
        let mut iter =
            RawMiGroupIterator::new(records.into_iter(), "MI").with_cell_tag(Some([b'C', b'B']));

        // First group: MI=1, no CB (key = "1\t")
        let result = iter.next().unwrap().unwrap();
        assert_eq!(result.0, "1\t");
        assert_eq!(result.1.len(), 2);

        // Second group: MI=1, CB=ACGT (key = "1\tACGT")
        let result = iter.next().unwrap().unwrap();
        assert_eq!(result.0, "1\tACGT");
        assert_eq!(result.1.len(), 1);

        assert!(iter.next().is_none());
    }

    #[test]
    fn test_raw_cell_tag_with_transform() {
        // Cell tag should work with MI transform (duplex-style /A /B stripping)
        let records: Vec<Result<Vec<u8>>> = vec![
            Ok(make_raw_bam_with_two_tags("MI", "1/A", "CB", "ACGT")),
            Ok(make_raw_bam_with_two_tags("MI", "1/B", "CB", "ACGT")),
            Ok(make_raw_bam_with_two_tags("MI", "1/A", "CB", "TGCA")),
        ];
        let mut iter = RawMiGroupIterator::with_transform(records.into_iter(), "MI", |raw| {
            let s = String::from_utf8_lossy(raw);
            extract_mi_base(&s).to_string()
        })
        .with_cell_tag(Some([b'C', b'B']));

        // All three have base MI "1", but different cells
        // First group: base MI=1, CB=ACGT (2 records: 1/A and 1/B)
        let result = iter.next().unwrap().unwrap();
        assert_eq!(result.0, "1\tACGT");
        assert_eq!(result.1.len(), 2);

        // Second group: base MI=1, CB=TGCA
        let result = iter.next().unwrap().unwrap();
        assert_eq!(result.0, "1\tTGCA");
        assert_eq!(result.1.len(), 1);

        assert!(iter.next().is_none());
    }

    // ========================================================================
    // RawMiGrouper tests (Grouper trait impl)
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

    /// Helper: create a parsed `DecodedRecord` (for testing error paths).
    fn make_parsed_decoded_record(mi: &str) -> DecodedRecord {
        let record = create_record_with_mi(mi);
        let key = crate::unified_pipeline::GroupKey::single(0, 0, 0, 0, 0, 0);
        DecodedRecord::new(record, key)
    }

    #[test]
    fn test_raw_grouper_single_mi_group() {
        let mut grouper = RawMiGrouper::new("MI", 10);

        let records = vec![
            make_raw_decoded_record("MI", "0"),
            make_raw_decoded_record("MI", "0"),
            make_raw_decoded_record("MI", "0"),
        ];

        let batches = grouper.add_records(records).unwrap();
        // Only 1 MI group so far, batch_size=10, no complete batch yet
        assert!(batches.is_empty());
        assert!(grouper.has_pending());

        let final_batch = grouper.finish().unwrap().unwrap();
        assert_eq!(final_batch.groups.len(), 1);
        assert_eq!(final_batch.groups[0].mi, "0");
        assert_eq!(final_batch.groups[0].records.len(), 3);
    }

    #[test]
    fn test_raw_grouper_multiple_mi_groups() {
        let mut grouper = RawMiGrouper::new("MI", 10);

        let records = vec![
            make_raw_decoded_record("MI", "0"),
            make_raw_decoded_record("MI", "0"),
            make_raw_decoded_record("MI", "1"),
            make_raw_decoded_record("MI", "1"),
            make_raw_decoded_record("MI", "1"),
            make_raw_decoded_record("MI", "2"),
        ];

        let batches = grouper.add_records(records).unwrap();
        assert!(batches.is_empty()); // 3 groups < batch_size 10

        let final_batch = grouper.finish().unwrap().unwrap();
        assert_eq!(final_batch.groups.len(), 3);
        assert_eq!(final_batch.groups[0].mi, "0");
        assert_eq!(final_batch.groups[0].records.len(), 2);
        assert_eq!(final_batch.groups[1].mi, "1");
        assert_eq!(final_batch.groups[1].records.len(), 3);
        assert_eq!(final_batch.groups[2].mi, "2");
        assert_eq!(final_batch.groups[2].records.len(), 1);
    }

    #[test]
    fn test_raw_grouper_batch_size_triggers() {
        let mut grouper = RawMiGrouper::new("MI", 2);

        let records = vec![
            make_raw_decoded_record("MI", "0"),
            make_raw_decoded_record("MI", "1"),
            make_raw_decoded_record("MI", "2"),
            make_raw_decoded_record("MI", "3"),
            make_raw_decoded_record("MI", "4"),
        ];

        let batches = grouper.add_records(records).unwrap();
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

        let final_batch = grouper.finish().unwrap().unwrap();
        assert_eq!(final_batch.groups.len(), 1);
        assert_eq!(final_batch.groups[0].mi, "4");
    }

    #[test]
    fn test_raw_grouper_groups_records_without_mi_under_empty_key() {
        let mut grouper = RawMiGrouper::new("MI", 10);

        let records = vec![
            make_raw_decoded_record("MI", "0"),
            make_raw_decoded_record_no_tag(),
            make_raw_decoded_record("MI", "0"),
        ];

        let batches = grouper.add_records(records).unwrap();
        assert!(batches.is_empty());

        // Records without MI are now grouped under "" (matching MiGrouper behavior)
        let final_batch = grouper.finish().unwrap().unwrap();
        assert_eq!(final_batch.groups.len(), 3);
        assert_eq!(final_batch.groups[0].mi, "0");
        assert_eq!(final_batch.groups[0].records.len(), 1);
        assert_eq!(final_batch.groups[1].mi, "");
        assert_eq!(final_batch.groups[1].records.len(), 1);
        assert_eq!(final_batch.groups[2].mi, "0");
        assert_eq!(final_batch.groups[2].records.len(), 1);
    }

    #[test]
    fn test_raw_grouper_rejects_parsed_records() {
        let mut grouper = RawMiGrouper::new("MI", 10);

        let records = vec![make_parsed_decoded_record("0")];
        let result = grouper.add_records(records);
        assert!(result.is_err());
    }

    #[test]
    fn test_raw_grouper_finish_empty() {
        let mut grouper = RawMiGrouper::new("MI", 10);
        assert!(!grouper.has_pending());

        let final_batch = grouper.finish().unwrap();
        assert!(final_batch.is_none());
    }

    #[test]
    fn test_raw_grouper_finish_idempotent() {
        let mut grouper = RawMiGrouper::new("MI", 10);

        let records = vec![make_raw_decoded_record("MI", "0")];
        grouper.add_records(records).unwrap();

        let batch1 = grouper.finish().unwrap();
        assert!(batch1.is_some());

        let batch2 = grouper.finish().unwrap();
        assert!(batch2.is_none());
    }

    #[test]
    fn test_raw_grouper_has_pending_states() {
        let mut grouper = RawMiGrouper::new("MI", 10);

        // Initially nothing pending
        assert!(!grouper.has_pending());

        // After adding records, current_mi is set
        let records = vec![make_raw_decoded_record("MI", "0")];
        grouper.add_records(records).unwrap();
        assert!(grouper.has_pending());

        // After adding a different MI, first group moves to pending_groups
        let records = vec![make_raw_decoded_record("MI", "1")];
        grouper.add_records(records).unwrap();
        assert!(grouper.has_pending());
    }

    #[test]
    fn test_raw_grouper_with_filter_and_transform() {
        // Build records with flag field set: use byte 14..16 for flag.
        // We create a filter that checks the flag field for secondary alignment (0x100).
        let mut grouper = RawMiGrouper::with_filter_and_transform(
            "MI",
            10,
            |bam: &[u8]| {
                // Check flag field at offset 14..16
                let flag = u16::from_le_bytes([bam[14], bam[15]]);
                flag & 0x100 == 0 // keep if NOT secondary
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

        let batches = grouper.add_records(records).unwrap();
        assert!(batches.is_empty());

        let final_batch = grouper.finish().unwrap().unwrap();
        // Secondary was filtered out, and 1/A and 1/B transform to "1"
        assert_eq!(final_batch.groups.len(), 1);
        assert_eq!(final_batch.groups[0].mi, "1");
        assert_eq!(final_batch.groups[0].records.len(), 2); // primary 1/A + 1/B
    }

    // ========================================================================
    // MiGrouper tests (Grouper trait impl)
    // ========================================================================

    // Helper: create a parsed `DecodedRecord` for `MiGrouper` tests.
    fn make_decoded_record_for_mi_grouper(mi: &str) -> DecodedRecord {
        let record = create_record_with_mi(mi);
        let key = crate::unified_pipeline::GroupKey::single(0, 0, 0, 0, 0, 0);
        DecodedRecord::new(record, key)
    }

    #[test]
    fn test_mi_grouper_single_group() {
        let mut grouper = MiGrouper::new("MI", 10);

        let records = vec![
            make_decoded_record_for_mi_grouper("0"),
            make_decoded_record_for_mi_grouper("0"),
            make_decoded_record_for_mi_grouper("0"),
        ];

        let batches = grouper.add_records(records).unwrap();
        assert!(batches.is_empty());
        assert!(grouper.has_pending());

        let final_batch = grouper.finish().unwrap().unwrap();
        assert_eq!(final_batch.groups.len(), 1);
        assert_eq!(final_batch.groups[0].mi, "0");
        assert_eq!(final_batch.groups[0].records.len(), 3);
    }

    #[test]
    fn test_mi_grouper_multiple_groups() {
        let mut grouper = MiGrouper::new("MI", 10);

        let records = vec![
            make_decoded_record_for_mi_grouper("0"),
            make_decoded_record_for_mi_grouper("0"),
            make_decoded_record_for_mi_grouper("1"),
            make_decoded_record_for_mi_grouper("1"),
            make_decoded_record_for_mi_grouper("1"),
            make_decoded_record_for_mi_grouper("2"),
        ];

        let batches = grouper.add_records(records).unwrap();
        assert!(batches.is_empty());

        let final_batch = grouper.finish().unwrap().unwrap();
        assert_eq!(final_batch.groups.len(), 3);
        assert_eq!(final_batch.groups[0].mi, "0");
        assert_eq!(final_batch.groups[0].records.len(), 2);
        assert_eq!(final_batch.groups[1].mi, "1");
        assert_eq!(final_batch.groups[1].records.len(), 3);
        assert_eq!(final_batch.groups[2].mi, "2");
        assert_eq!(final_batch.groups[2].records.len(), 1);
    }

    #[test]
    fn test_mi_grouper_batch_size_triggers() {
        let mut grouper = MiGrouper::new("MI", 2);

        let records = vec![
            make_decoded_record_for_mi_grouper("0"),
            make_decoded_record_for_mi_grouper("1"),
            make_decoded_record_for_mi_grouper("2"),
            make_decoded_record_for_mi_grouper("3"),
            make_decoded_record_for_mi_grouper("4"),
        ];

        let batches = grouper.add_records(records).unwrap();
        assert_eq!(batches.len(), 2);
        assert_eq!(batches[0].groups.len(), 2);
        assert_eq!(batches[0].groups[0].mi, "0");
        assert_eq!(batches[0].groups[1].mi, "1");
        assert_eq!(batches[1].groups.len(), 2);
        assert_eq!(batches[1].groups[0].mi, "2");
        assert_eq!(batches[1].groups[1].mi, "3");

        assert!(grouper.has_pending());

        let final_batch = grouper.finish().unwrap().unwrap();
        assert_eq!(final_batch.groups.len(), 1);
        assert_eq!(final_batch.groups[0].mi, "4");
    }

    #[test]
    fn test_mi_grouper_rejects_raw_records() {
        let mut grouper = MiGrouper::new("MI", 10);

        let records = vec![make_raw_decoded_record("MI", "0")];
        let result = grouper.add_records(records);
        assert!(result.is_err());
    }

    #[test]
    fn test_mi_grouper_finish_empty() {
        let mut grouper = MiGrouper::new("MI", 10);
        assert!(!grouper.has_pending());

        let final_batch = grouper.finish().unwrap();
        assert!(final_batch.is_none());
    }

    #[test]
    fn test_mi_grouper_finish_idempotent() {
        let mut grouper = MiGrouper::new("MI", 10);

        let records = vec![make_decoded_record_for_mi_grouper("0")];
        grouper.add_records(records).unwrap();

        let batch1 = grouper.finish().unwrap();
        assert!(batch1.is_some());

        let batch2 = grouper.finish().unwrap();
        assert!(batch2.is_none());
    }

    #[test]
    fn test_mi_grouper_with_filter_and_transform() {
        let mut grouper = MiGrouper::with_filter_and_transform(
            "MI",
            10,
            |_r| true, // keep all records
            |mi| extract_mi_base(mi).to_string(),
        );

        let records = vec![
            make_decoded_record_for_mi_grouper("1/A"),
            make_decoded_record_for_mi_grouper("1/B"),
            make_decoded_record_for_mi_grouper("2/A"),
        ];

        let batches = grouper.add_records(records).unwrap();
        assert!(batches.is_empty());

        let final_batch = grouper.finish().unwrap().unwrap();
        assert_eq!(final_batch.groups.len(), 2);
        assert_eq!(final_batch.groups[0].mi, "1");
        assert_eq!(final_batch.groups[0].records.len(), 2);
        assert_eq!(final_batch.groups[1].mi, "2");
        assert_eq!(final_batch.groups[1].records.len(), 1);
    }

    #[test]
    fn test_mi_grouper_has_pending_states() {
        let mut grouper = MiGrouper::new("MI", 10);

        assert!(!grouper.has_pending());

        let records = vec![make_decoded_record_for_mi_grouper("0")];
        grouper.add_records(records).unwrap();
        assert!(grouper.has_pending());

        let records = vec![make_decoded_record_for_mi_grouper("1")];
        grouper.add_records(records).unwrap();
        assert!(grouper.has_pending());
    }

    #[test]
    fn test_mi_grouper_record_without_mi_gets_empty_string() {
        let mut grouper = MiGrouper::new("MI", 10);

        let record = create_record_without_mi();
        let key = crate::unified_pipeline::GroupKey::single(0, 0, 0, 0, 0, 0);
        let records = vec![DecodedRecord::new(record, key)];

        let batches = grouper.add_records(records).unwrap();
        assert!(batches.is_empty());

        // Record without MI tag gets default empty string as MI value
        let final_batch = grouper.finish().unwrap().unwrap();
        assert_eq!(final_batch.groups.len(), 1);
        assert_eq!(final_batch.groups[0].mi, "");
        assert_eq!(final_batch.groups[0].records.len(), 1);
    }

    // ========================================================================
    // Batch type tests
    // ========================================================================

    #[test]
    fn test_mi_group_batch_new() {
        let batch = MiGroupBatch::new();
        assert!(batch.is_empty());
        assert_eq!(batch.len(), 0);
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

        batch.groups.push(MiGroup::new("0".to_string(), vec![create_record_with_mi("0")]));
        assert!(!batch.is_empty());
        assert_eq!(batch.len(), 1);

        batch.groups.push(MiGroup::new(
            "1".to_string(),
            vec![create_record_with_mi("1"), create_record_with_mi("1")],
        ));
        assert_eq!(batch.len(), 2);
    }

    #[test]
    fn test_mi_group_batch_clear() {
        let mut batch = MiGroupBatch::new();
        batch.groups.push(MiGroup::new("0".to_string(), vec![create_record_with_mi("0")]));
        assert!(!batch.is_empty());

        batch.clear();
        assert!(batch.is_empty());
        assert_eq!(batch.len(), 0);
    }

    #[test]
    fn test_mi_group_batch_weight() {
        let mut batch = MiGroupBatch::new();
        assert_eq!(batch.batch_weight(), 0);

        batch.groups.push(MiGroup::new(
            "0".to_string(),
            vec![create_record_with_mi("0"), create_record_with_mi("0")],
        ));
        assert_eq!(batch.batch_weight(), 2);

        batch.groups.push(MiGroup::new("1".to_string(), vec![create_record_with_mi("1")]));
        assert_eq!(batch.batch_weight(), 3);
    }

    #[test]
    fn test_mi_group_batch_default() {
        let batch = MiGroupBatch::default();
        assert!(batch.is_empty());
    }

    #[test]
    fn test_mi_group_batch_memory_estimate() {
        let batch = MiGroupBatch::new();
        // Empty batch should have minimal heap size (just vec overhead)
        let _ = batch.estimate_heap_size(); // Should not panic

        let mut batch = MiGroupBatch::new();
        batch.groups.push(MiGroup::new("0".to_string(), vec![create_record_with_mi("0")]));
        let size = batch.estimate_heap_size();
        assert!(size > 0);
    }

    #[test]
    fn test_mi_group_new() {
        let group = MiGroup::new("42".to_string(), vec![create_record_with_mi("42")]);
        assert_eq!(group.mi, "42");
        assert_eq!(group.records.len(), 1);
    }

    #[test]
    fn test_mi_group_batch_weight_single() {
        let group = MiGroup::new(
            "0".to_string(),
            vec![
                create_record_with_mi("0"),
                create_record_with_mi("0"),
                create_record_with_mi("0"),
            ],
        );
        assert_eq!(group.batch_weight(), 3);
    }

    #[test]
    fn test_mi_group_memory_estimate() {
        let group = MiGroup::new("0".to_string(), vec![create_record_with_mi("0")]);
        let size = group.estimate_heap_size();
        assert!(size > 0);
    }

    // ========================================================================
    // Raw batch type tests
    // ========================================================================

    #[test]
    fn test_raw_mi_group_new() {
        let raw = make_raw_bam_with_tag("MI", "42");
        let group = RawMiGroup::new("42".to_string(), vec![raw]);
        assert_eq!(group.mi, "42");
        assert_eq!(group.records.len(), 1);
    }

    #[test]
    fn test_raw_mi_group_batch_weight() {
        let group = RawMiGroup::new(
            "0".to_string(),
            vec![make_raw_bam_with_tag("MI", "0"), make_raw_bam_with_tag("MI", "0")],
        );
        assert_eq!(group.batch_weight(), 2);
    }

    #[test]
    fn test_raw_mi_group_memory_estimate() {
        let group = RawMiGroup::new("0".to_string(), vec![make_raw_bam_with_tag("MI", "0")]);
        let size = group.estimate_heap_size();
        assert!(size > 0);
    }

    #[test]
    fn test_raw_mi_group_batch_new() {
        let batch = RawMiGroupBatch::new();
        assert!(batch.groups.is_empty());
    }

    #[test]
    fn test_raw_mi_group_batch_default() {
        let batch = RawMiGroupBatch::default();
        assert!(batch.groups.is_empty());
    }

    #[test]
    fn test_raw_mi_group_batch_weight_method() {
        let mut batch = RawMiGroupBatch::new();
        assert_eq!(batch.batch_weight(), 0);

        batch.groups.push(RawMiGroup::new(
            "0".to_string(),
            vec![make_raw_bam_with_tag("MI", "0"), make_raw_bam_with_tag("MI", "0")],
        ));
        assert_eq!(batch.batch_weight(), 2);

        batch.groups.push(RawMiGroup::new("1".to_string(), vec![make_raw_bam_with_tag("MI", "1")]));
        assert_eq!(batch.batch_weight(), 3);
    }

    #[test]
    fn test_raw_mi_group_batch_memory_estimate() {
        let batch = RawMiGroupBatch::new();
        let _ = batch.estimate_heap_size(); // Should not panic

        let mut batch = RawMiGroupBatch::new();
        batch.groups.push(RawMiGroup::new("0".to_string(), vec![make_raw_bam_with_tag("MI", "0")]));
        let size = batch.estimate_heap_size();
        assert!(size > 0);
    }
}

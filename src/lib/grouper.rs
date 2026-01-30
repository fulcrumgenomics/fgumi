//! Grouper implementations for the 9-step pipeline.
//!
//! Each grouper is responsible for grouping decoded BAM records according
//! to command-specific rules and emitting complete groups for processing.
//!
//! Note: BAM record parsing is now handled by the Decode step in the pipeline.
//! Groupers receive pre-decoded `RecordBuf` vectors.

use std::io;

use noodles::sam::alignment::RecordBuf;

use crate::template::{Template, TemplateBatch};
use crate::unified_pipeline::{BatchWeight, DecodedRecord, Grouper, MemoryEstimate};

// ============================================================================
// BatchWeight Implementations
// ============================================================================

/// Single records have weight 1 (for batch_size-like behavior).
impl BatchWeight for RecordBuf {
    fn batch_weight(&self) -> usize {
        1
    }
}

/// Template batches have weight equal to the number of templates.
impl BatchWeight for TemplateBatch {
    fn batch_weight(&self) -> usize {
        self.len()
    }
}

// ============================================================================
// SingleRecordGrouper
// ============================================================================

/// A grouper that emits each record as its own "group".
///
/// Used by commands that process records independently:
/// - filter: Apply filters to each record
/// - clip: Clip each record
/// - correct: Correct UMIs in each record
#[derive(Default)]
pub struct SingleRecordGrouper;

impl SingleRecordGrouper {
    /// Create a new single-record grouper.
    #[must_use]
    pub fn new() -> Self {
        Self
    }
}

impl Grouper for SingleRecordGrouper {
    type Group = RecordBuf;

    fn add_records(&mut self, records: Vec<DecodedRecord>) -> io::Result<Vec<Self::Group>> {
        // Each record is its own group - extract and return them directly
        Ok(records.into_iter().map(|d| d.record).collect())
    }

    fn finish(&mut self) -> io::Result<Option<Self::Group>> {
        // No state to flush - records are passed through immediately
        Ok(None)
    }

    fn has_pending(&self) -> bool {
        // Never has pending state
        false
    }
}

// ============================================================================
// TemplateGrouper
// ============================================================================

use std::collections::VecDeque;

/// A Grouper that batches BAM records into Templates by QNAME.
///
/// Input BAM must be query-name sorted or grouped by QNAME.
/// Records with the same QNAME are grouped into a [`Template`], then
/// templates are batched for efficient parallel processing.
///
/// # Example
///
/// ```ignore
/// use fgumi_lib::grouper::TemplateGrouper;
/// use fgumi_lib::unified_pipeline::Grouper;
///
/// let grouper = TemplateGrouper::new(1000); // 1000 templates per batch
/// // Use with run_bam_pipeline_with_grouper...
/// ```
pub struct TemplateGrouper {
    /// Number of templates per batch.
    batch_size: usize,
    /// Current template name being accumulated.
    current_name: Option<Vec<u8>>,
    /// Hash of current template name (for fast comparison).
    current_name_hash: Option<u64>,
    /// Records for current template.
    current_records: Vec<RecordBuf>,
    /// Completed templates waiting to be batched.
    pending_templates: VecDeque<Template>,
}

impl TemplateGrouper {
    /// Create a new [`TemplateGrouper`].
    ///
    /// # Arguments
    /// * `batch_size` - Number of templates per batch (e.g., 1000)
    #[must_use]
    pub fn new(batch_size: usize) -> Self {
        Self {
            batch_size: batch_size.max(1),
            current_name: None,
            current_name_hash: None,
            current_records: Vec::new(),
            pending_templates: VecDeque::new(),
        }
    }

    /// Flush current template to pending queue if non-empty.
    fn flush_current_template(&mut self) -> io::Result<()> {
        if !self.current_records.is_empty() {
            let template = Template::from_records(std::mem::take(&mut self.current_records))
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            self.pending_templates.push_back(template);
            self.current_name = None;
            self.current_name_hash = None;
        }
        Ok(())
    }

    /// Collect completed batches from pending templates.
    fn collect_batches(&mut self) -> Vec<TemplateBatch> {
        let mut batches = Vec::new();

        while self.pending_templates.len() >= self.batch_size {
            let mut batch = Vec::with_capacity(self.batch_size);
            for _ in 0..self.batch_size {
                if let Some(template) = self.pending_templates.pop_front() {
                    batch.push(template);
                }
            }
            batches.push(batch);
        }

        batches
    }
}

impl Grouper for TemplateGrouper {
    type Group = TemplateBatch;

    fn add_records(&mut self, records: Vec<DecodedRecord>) -> io::Result<Vec<Self::Group>> {
        // Group records by QNAME (using pre-computed name_hash for fast comparison)
        for decoded in records {
            let record = decoded.record;
            let name_hash = decoded.key.name_hash;
            let name = record.name().map(|n| Vec::from(<_ as AsRef<[u8]>>::as_ref(n)));

            match (self.current_name_hash, name_hash) {
                (Some(current_hash), new_hash) if current_hash == new_hash => {
                    // Same template (based on hash) - add record
                    self.current_records.push(record);
                }
                (Some(_), _) => {
                    // Different template - flush current and start new
                    self.flush_current_template()?;
                    self.current_name = name;
                    self.current_name_hash = Some(name_hash);
                    self.current_records.push(record);
                }
                (None, _) => {
                    // First record
                    self.current_name = name;
                    self.current_name_hash = Some(name_hash);
                    self.current_records.push(record);
                }
            }
        }

        // Return any completed batches
        Ok(self.collect_batches())
    }

    fn finish(&mut self) -> io::Result<Option<Self::Group>> {
        // Flush any remaining template
        self.flush_current_template()?;

        // Return remaining templates as a final batch (may be < batch_size)
        if self.pending_templates.is_empty() {
            return Ok(None);
        }

        let mut batch = Vec::with_capacity(self.pending_templates.len());
        while let Some(template) = self.pending_templates.pop_front() {
            batch.push(template);
        }
        Ok(Some(batch))
    }

    fn has_pending(&self) -> bool {
        self.current_name.is_some() || !self.pending_templates.is_empty()
    }
}

// ============================================================================
// PositionGrouper
// ============================================================================

use noodles::sam::alignment::record::data::field::Tag;

use crate::read_info::LibraryLookup;
use crate::unified_pipeline::GroupKey;

/// A position group containing all templates at the same genomic position.
#[derive(Debug)]
pub struct PositionGroup {
    /// The position key for this group (pre-computed from first record).
    pub group_key: GroupKey,
    /// All templates at this position.
    pub templates: Vec<Template>,
    /// Pre-assigned MI base offset for global uniqueness.
    /// Each position group gets a unique range of MIs: `[base_mi, base_mi + num_umis)`.
    pub base_mi: u64,
}

impl BatchWeight for PositionGroup {
    /// Returns the number of templates in this group.
    fn batch_weight(&self) -> usize {
        self.templates.len()
    }
}

impl MemoryEstimate for PositionGroup {
    fn estimate_heap_size(&self) -> usize {
        // group_key: GroupKey is a Copy type (no heap allocation)
        // templates: Vec<Template>
        let templates_size: usize =
            self.templates.iter().map(MemoryEstimate::estimate_heap_size).sum();
        let templates_vec_overhead = self.templates.capacity() * std::mem::size_of::<Template>();
        // base_mi: u64 is inline
        templates_size + templates_vec_overhead
    }
}

impl MemoryEstimate for ProcessedPositionGroup {
    fn estimate_heap_size(&self) -> usize {
        // templates: Vec<Template>
        let templates_size: usize =
            self.templates.iter().map(MemoryEstimate::estimate_heap_size).sum();
        let templates_vec_overhead = self.templates.capacity() * std::mem::size_of::<Template>();

        // family_sizes: AHashMap<usize, u64>
        // Each entry is ~24 bytes (key + value + overhead)
        let family_sizes_size = self.family_sizes.len() * 24;

        // filter_metrics: FilterMetrics is mostly inline (u64 fields)
        // Just a small struct overhead

        templates_size + templates_vec_overhead + family_sizes_size
    }
}

/// Metrics tracking what was filtered during processing.
#[derive(Default, Clone, Debug)]
pub struct FilterMetrics {
    /// Total templates seen before filtering.
    pub total_templates: u64,
    /// Templates accepted after filtering.
    pub accepted_templates: u64,
    /// Templates discarded because they weren't passing filter.
    pub discarded_non_pf: u64,
    /// Templates discarded due to poor alignment.
    pub discarded_poor_alignment: u64,
    /// Templates discarded due to Ns in UMI.
    pub discarded_ns_in_umi: u64,
    /// Templates discarded due to UMI being too short.
    pub discarded_umi_too_short: u64,
}

impl FilterMetrics {
    /// Create new empty metrics.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Merge another `FilterMetrics` into this one.
    pub fn merge(&mut self, other: &FilterMetrics) {
        self.total_templates += other.total_templates;
        self.accepted_templates += other.accepted_templates;
        self.discarded_non_pf += other.discarded_non_pf;
        self.discarded_poor_alignment += other.discarded_poor_alignment;
        self.discarded_ns_in_umi += other.discarded_ns_in_umi;
        self.discarded_umi_too_short += other.discarded_umi_too_short;
    }
}

/// Result of processing a position group through UMI assignment.
#[derive(Debug)]
pub struct ProcessedPositionGroup {
    /// Templates with MI tags assigned, sorted by MI then name.
    /// `Template.mi` contains local IDs (0, 1, 2, ...) - add `base_mi` for global IDs.
    pub templates: Vec<Template>,
    /// Family size counts for this position group.
    pub family_sizes: ahash::AHashMap<usize, u64>,
    /// Filter metrics for this position group (thread-local, merged later).
    pub filter_metrics: FilterMetrics,
    /// Total input records processed (for progress tracking).
    pub input_record_count: u64,
    /// Base MI offset for this position group (for global uniqueness).
    pub base_mi: u64,
}

/// Configuration for position grouping.
#[derive(Clone)]
pub struct PositionGrouperConfig {
    /// Cell barcode tag (e.g., "CB").
    pub cell_tag: Tag,
    /// Library lookup from read group to library.
    pub library_lookup: LibraryLookup,
    /// Whether to include secondary and supplementary reads in templates.
    /// When true, secondary/supplementary reads are included in the template
    /// that contains their corresponding primary read. Requires the input to be
    /// template-coordinate sorted with `pa` tags on secondary/supplementary reads.
    /// Default: false (secondary/supplementary reads are skipped).
    pub include_secondary_supplementary: bool,
}

impl PositionGrouperConfig {
    /// Creates a new configuration.
    #[must_use]
    pub fn new(cell_tag: Tag, library_lookup: LibraryLookup) -> Self {
        Self { cell_tag, library_lookup, include_secondary_supplementary: false }
    }

    /// Creates a new configuration that includes secondary/supplementary reads.
    #[must_use]
    pub fn with_secondary_supplementary(cell_tag: Tag, library_lookup: LibraryLookup) -> Self {
        Self { cell_tag, library_lookup, include_secondary_supplementary: true }
    }
}

/// A grouper that groups templates by genomic position.
///
/// Used by the `group` command to collect all templates at the same position
/// for UMI assignment. Expects input sorted by template-coordinate (mates adjacent,
/// sorted by position).
///
/// Groups records by query name into Templates, then groups Templates by
/// position (using pre-computed `GroupKey`) into position groups.
pub struct PositionGrouper {
    /// Configuration for grouping.
    config: PositionGrouperConfig,
    /// Current template name hash being accumulated.
    current_name_hash: Option<u64>,
    /// Records for current template.
    current_template_records: Vec<RecordBuf>,
    /// `GroupKeys` for all records in the current template.
    current_template_keys: Vec<GroupKey>,
    /// Current position group key being accumulated.
    current_position_key: Option<GroupKey>,
    /// Templates at current position.
    current_templates: Vec<Template>,
    /// Next available MI base offset for global uniqueness.
    /// Incremented by template count when emitting each position group.
    next_mi_base: u64,
}

impl PositionGrouper {
    /// Create a new position grouper with the specified configuration.
    #[must_use]
    pub fn new(config: PositionGrouperConfig) -> Self {
        Self {
            config,
            current_name_hash: None,
            current_template_records: Vec::new(),
            current_template_keys: Vec::new(),
            current_position_key: None,
            current_templates: Vec::new(),
            next_mi_base: 0,
        }
    }

    /// Combine multiple `GroupKeys` from a template's records into a single position key.
    ///
    /// For paired-end reads without MC tag, each record only has its own position.
    /// We combine them to get the full position info.
    ///
    /// Secondary and supplementary reads have default UNKNOWN position values in their keys
    /// (because `compute_group_key` skips them). We filter these out before combining
    /// to avoid corrupting the template's position key.
    ///
    /// This implementation avoids allocation by finding valid keys through iteration
    /// rather than collecting into a Vec.
    fn combine_keys(keys: &[GroupKey]) -> GroupKey {
        if keys.is_empty() {
            return GroupKey::default();
        }

        // Fast path: single key (very common case)
        if keys.len() == 1 {
            return if keys[0].ref_id1 == GroupKey::UNKNOWN_REF {
                // Preserve library/cell metadata to avoid merging across libraries
                let k = keys[0];
                GroupKey {
                    library_idx: k.library_idx,
                    cell_hash: k.cell_hash,
                    name_hash: k.name_hash,
                    ..GroupKey::default()
                }
            } else {
                keys[0]
            };
        }

        // Find first two valid keys without allocation.
        // Valid keys have ref_id1 != UNKNOWN_REF (i.e., not secondary/supplementary).
        let mut first: Option<&GroupKey> = None;
        let mut second: Option<&GroupKey> = None;

        for key in keys {
            if key.ref_id1 != GroupKey::UNKNOWN_REF {
                if first.is_none() {
                    // Check if this key has valid mate info - if so, use it directly
                    if key.ref_id2 != GroupKey::UNKNOWN_REF {
                        return *key;
                    }
                    first = Some(key);
                } else if second.is_none() {
                    second = Some(key);
                    break; // We have two valid keys, no need to continue
                }
            }
        }

        // Handle cases based on how many valid keys we found
        match (first, second) {
            (None, _) => {
                // No valid keys - preserve library/cell/name metadata for template grouping
                let k = keys[0];
                GroupKey {
                    library_idx: k.library_idx,
                    cell_hash: k.cell_hash,
                    name_hash: k.name_hash,
                    ..GroupKey::default()
                }
            }
            (Some(f), None) => {
                // Single valid key
                *f
            }
            (Some(f), Some(s)) => {
                // Two valid keys without mate info - combine their positions
                GroupKey::paired(
                    f.ref_id1,
                    f.pos1,
                    f.strand1,
                    s.ref_id1,
                    s.pos1,
                    s.strand1,
                    f.library_idx,
                    f.cell_hash,
                    f.name_hash,
                )
            }
        }
    }

    /// Build a Template from accumulated records.
    fn build_current_template(&mut self) -> io::Result<Option<(Template, GroupKey)>> {
        if self.current_template_records.is_empty() {
            return Ok(None);
        }

        let records = std::mem::take(&mut self.current_template_records);
        let keys = std::mem::take(&mut self.current_template_keys);
        self.current_name_hash = None;

        let combined_key = Self::combine_keys(&keys);

        Template::from_records(records)
            .map(|t| Some((t, combined_key)))
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }

    /// Process a single record, potentially emitting completed groups.
    fn process_record(&mut self, decoded: DecodedRecord) -> io::Result<Vec<PositionGroup>> {
        let record = decoded.record;
        let key = decoded.key;

        // Skip secondary and supplementary reads unless configured to include them
        let flags = record.flags();
        if !self.config.include_secondary_supplementary
            && (flags.is_secondary() || flags.is_supplementary())
        {
            return Ok(Vec::new());
        }

        let mut completed_groups = Vec::new();
        let name_hash = key.name_hash;

        // Check if this record belongs to the current template (using name hash)
        let same_template = self.current_name_hash.is_none_or(|h| h == name_hash);

        if same_template || self.current_template_records.is_empty() {
            // Add to current template
            if self.current_name_hash.is_none() {
                self.current_name_hash = Some(name_hash);
            }
            self.current_template_keys.push(key);
            self.current_template_records.push(record);
        } else {
            // New template - finish current one
            if let Some((template, template_key)) = self.build_current_template()? {
                // Compare position keys (ignoring name_hash)
                let template_pos_key = template_key.position_key();

                if let Some(ref current_key) = self.current_position_key {
                    if current_key.position_key() == template_pos_key {
                        // Same position - add to current group
                        self.current_templates.push(template);
                    } else {
                        // New position - emit current group
                        let finished_templates = std::mem::take(&mut self.current_templates);
                        let finished_key =
                            self.current_position_key.take().expect("inside Some branch");

                        // Assign MI base and reserve range for this group
                        let base_mi = self.next_mi_base;
                        self.next_mi_base += finished_templates.len() as u64;

                        completed_groups.push(PositionGroup {
                            group_key: finished_key,
                            templates: finished_templates,
                            base_mi,
                        });

                        // Start new position group
                        self.current_position_key = Some(template_key);
                        self.current_templates.push(template);
                    }
                } else {
                    // First template - start new position group
                    self.current_position_key = Some(template_key);
                    self.current_templates.push(template);
                }
            }

            // Start new template
            self.current_name_hash = Some(name_hash);
            self.current_template_keys.push(key);
            self.current_template_records.push(record);
        }

        Ok(completed_groups)
    }
}

impl Grouper for PositionGrouper {
    type Group = PositionGroup;

    fn add_records(&mut self, records: Vec<DecodedRecord>) -> io::Result<Vec<Self::Group>> {
        // Process each record using pre-computed GroupKey for fast comparison
        let mut completed_groups = Vec::new();
        for decoded in records {
            let groups = self.process_record(decoded)?;
            completed_groups.extend(groups);
        }

        Ok(completed_groups)
    }

    fn finish(&mut self) -> io::Result<Option<Self::Group>> {
        // Build any remaining template
        if let Some((template, template_key)) = self.build_current_template()? {
            let template_pos_key = template_key.position_key();

            if let Some(ref current_key) = self.current_position_key {
                if current_key.position_key() == template_pos_key {
                    // Same position - add to current group
                    self.current_templates.push(template);
                } else {
                    // Different position - this shouldn't happen at EOF, but handle it
                    // Just add to current group as a fallback
                    self.current_templates.push(template);
                }
            } else {
                // First template
                self.current_position_key = Some(template_key);
                self.current_templates.push(template);
            }
        }

        // Emit final position group
        if !self.current_templates.is_empty() {
            if let Some(key) = self.current_position_key.take() {
                let templates = std::mem::take(&mut self.current_templates);

                // Assign MI base and reserve range for this group
                let base_mi = self.next_mi_base;
                self.next_mi_base += templates.len() as u64;

                return Ok(Some(PositionGroup { group_key: key, templates, base_mi }));
            }
        }

        Ok(None)
    }

    fn has_pending(&self) -> bool {
        !self.current_template_records.is_empty() || !self.current_templates.is_empty()
    }
}

// ============================================================================
// PositionMiGrouper - Groups by position then by MI for duplex metrics
// ============================================================================

use crate::umi::extract_mi_base;
use noodles::sam::alignment::record::Cigar;
use noodles::sam::alignment::record::cigar::op::Kind;
use std::collections::HashMap;

/// Coordinate grouping key matching fgbio's `ReadInfo` structure.
/// Uses unclipped 5' positions of both R1 and R2, normalized so lower position comes first.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct ReadInfoKey {
    /// Reference index for position 1 (lower).
    pub ref_index1: Option<usize>,
    /// Unclipped 5' position for position 1.
    pub start1: i32,
    /// Strand for position 1 (true = reverse).
    pub strand1: bool,
    /// Reference index for position 2 (higher).
    pub ref_index2: Option<usize>,
    /// Unclipped 5' position for position 2.
    pub start2: i32,
    /// Strand for position 2 (true = reverse).
    pub strand2: bool,
}

impl ReadInfoKey {
    /// Create an "unknown" key for unpaired/unmapped reads (excluded from metrics).
    #[must_use]
    pub fn unknown() -> Self {
        Self {
            ref_index1: None,
            start1: 0,
            strand1: false,
            ref_index2: None,
            start2: 0,
            strand2: false,
        }
    }

    /// Returns true if this key represents a valid paired read.
    #[must_use]
    pub fn is_valid(&self) -> bool {
        self.ref_index1.is_some() && self.ref_index2.is_some()
    }

    /// Compute `ReadInfoKey` from template records (R1 + R2).
    ///
    /// Uses same filtering criteria as duplex-metrics:
    /// - R1: segmented, mapped, mate mapped, first segment, primary
    /// - R2: segmented, mapped, mate mapped, last segment, primary
    ///
    /// Returns `ReadInfoKey::unknown()` if R1/R2 cannot be identified or are on different references.
    #[must_use]
    pub fn from_template_records(records: &[RecordBuf]) -> Self {
        Self::compute_from_iter(records.iter())
    }

    /// Compute `ReadInfoKey` from template record references (R1 + R2).
    ///
    /// Same as `from_template_records` but takes references to avoid cloning.
    #[must_use]
    pub fn from_template_refs(records: &[&RecordBuf]) -> Self {
        Self::compute_from_iter(records.iter().copied())
    }

    /// Compute `ReadInfoKey` directly from R1 and R2 records.
    ///
    /// This is the most efficient method when R1 and R2 are already identified.
    /// Returns `ReadInfoKey::unknown()` if records are on different references.
    #[must_use]
    pub fn from_r1_r2(r1: &RecordBuf, r2: &RecordBuf) -> Self {
        let r1_ref = r1.reference_sequence_id();
        let r2_ref = r2.reference_sequence_id();

        // Only process if same reference
        if r1_ref != r2_ref || r1_ref.is_none() {
            return Self::unknown();
        }

        let r1_5prime = unclipped_five_prime_position(r1);
        let r2_5prime = unclipped_five_prime_position(r2);
        let r1_strand = r1.flags().is_reverse_complemented();
        let r2_strand = r2.flags().is_reverse_complemented();

        match (r1_5prime, r2_5prime) {
            (Some(s1), Some(s2)) => {
                // Normalize: lower position comes first (matching fgbio)
                if (r1_ref, s1) <= (r2_ref, s2) {
                    Self {
                        ref_index1: r1_ref,
                        start1: s1,
                        strand1: r1_strand,
                        ref_index2: r2_ref,
                        start2: s2,
                        strand2: r2_strand,
                    }
                } else {
                    Self {
                        ref_index1: r2_ref,
                        start1: s2,
                        strand1: r2_strand,
                        ref_index2: r1_ref,
                        start2: s1,
                        strand2: r1_strand,
                    }
                }
            }
            _ => Self::unknown(),
        }
    }

    /// Internal helper to compute ReadInfoKey from an iterator of record references.
    fn compute_from_iter<'a>(records: impl Iterator<Item = &'a RecordBuf>) -> Self {
        let records: Vec<&RecordBuf> = records.collect();

        // Find R1: segmented, mapped, mate mapped, first segment, primary
        let r1 = records.iter().find(|r| {
            let f = r.flags();
            f.is_segmented()
                && !f.is_unmapped()
                && !f.is_mate_unmapped()
                && f.is_first_segment()
                && !f.is_secondary()
                && !f.is_supplementary()
        });

        // Find R2: segmented, mapped, mate mapped, last segment, primary
        let r2 = records.iter().find(|r| {
            let f = r.flags();
            f.is_segmented()
                && !f.is_unmapped()
                && !f.is_mate_unmapped()
                && f.is_last_segment()
                && !f.is_secondary()
                && !f.is_supplementary()
        });

        match (r1, r2) {
            (Some(r1), Some(r2)) => {
                let r1_ref = r1.reference_sequence_id();
                let r2_ref = r2.reference_sequence_id();

                // Only process if same reference
                if r1_ref != r2_ref || r1_ref.is_none() {
                    return Self::unknown();
                }

                let r1_5prime = unclipped_five_prime_position(r1);
                let r2_5prime = unclipped_five_prime_position(r2);
                let r1_strand = r1.flags().is_reverse_complemented();
                let r2_strand = r2.flags().is_reverse_complemented();

                match (r1_5prime, r2_5prime) {
                    (Some(s1), Some(s2)) => {
                        // Normalize: lower position comes first (matching fgbio)
                        if (r1_ref, s1) <= (r2_ref, s2) {
                            Self {
                                ref_index1: r1_ref,
                                start1: s1,
                                strand1: r1_strand,
                                ref_index2: r2_ref,
                                start2: s2,
                                strand2: r2_strand,
                            }
                        } else {
                            Self {
                                ref_index1: r2_ref,
                                start1: s2,
                                strand1: r2_strand,
                                ref_index2: r1_ref,
                                start2: s1,
                                strand2: r1_strand,
                            }
                        }
                    }
                    _ => Self::unknown(),
                }
            }
            _ => Self::unknown(),
        }
    }
}

/// Computes the unclipped 5' position for a read, matching fgbio's `positionOf`.
///
/// For forward strand reads: `unclipped_start = alignment_start - leading_soft_clips`
/// For reverse strand reads: `unclipped_end = alignment_end + trailing_soft_clips`
fn unclipped_five_prime_position(record: &RecordBuf) -> Option<i32> {
    let is_reverse = record.flags().is_reverse_complemented();
    let cigar = record.cigar();

    if is_reverse {
        // For reverse strand, 5' is at the end
        // unclippedEnd = alignmentEnd + trailing soft clips
        let alignment_end = record.alignment_end().map(|p| usize::from(p) as i32)?;

        // Count trailing soft clips (last element if it's S)
        let trailing_clips: i32 = cigar
            .iter()
            .filter_map(std::result::Result::ok)
            .last()
            .filter(|op| op.kind() == Kind::SoftClip)
            .map_or(0, |op| op.len() as i32);

        Some(alignment_end + trailing_clips)
    } else {
        // For forward strand, 5' is at the start
        // unclippedStart = alignmentStart - leading soft clips
        let alignment_start = record.alignment_start().map(|p| usize::from(p) as i32)?;

        // Count leading soft clips (first element if it's S)
        let leading_clips: i32 = cigar
            .iter()
            .find_map(std::result::Result::ok)
            .filter(|op| op.kind() == Kind::SoftClip)
            .map_or(0, |op| op.len() as i32);

        Some(alignment_start - leading_clips)
    }
}

/// A position-MI group containing records at the same position, organized by MI.
///
/// This structure is used for collecting duplex metrics during consensus calling.
/// It groups records first by genomic position (using `ReadInfoKey`), then by
/// molecular identifier (MI).
///
/// Family size definitions:
/// - **CS** (Coordinate & Strand): Total templates at this position (`total_templates`)
/// - **SS** (Single Strand): Templates with same position AND full MI tag (including /A or /B)
/// - **DS** (Double Strand): Templates with same position AND base MI (without /A, /B)
#[derive(Debug)]
pub struct PositionMiGroup {
    /// The coordinate key for this group (computed from R1+R2 unclipped 5' positions).
    pub read_info_key: ReadInfoKey,
    /// Records organized by base MI (without /A, /B suffix).
    /// Each entry maps base_mi -> records with that MI.
    pub mi_groups: HashMap<String, Vec<RecordBuf>>,
    /// Total number of templates at this position (CS family size).
    pub total_templates: usize,
    /// Total number of records at this position.
    pub total_records: usize,
}

impl BatchWeight for PositionMiGroup {
    /// Returns the number of templates in this group.
    fn batch_weight(&self) -> usize {
        self.total_templates
    }
}

impl MemoryEstimate for PositionMiGroup {
    fn estimate_heap_size(&self) -> usize {
        // read_info_key: ReadInfoKey has minimal heap allocation
        // mi_groups: HashMap<String, Vec<RecordBuf>>
        let mi_groups_size: usize = self
            .mi_groups
            .iter()
            .map(|(k, v)| {
                let key_size = k.capacity();
                let vec_size: usize = v.iter().map(MemoryEstimate::estimate_heap_size).sum();
                let vec_overhead = v.capacity() * std::mem::size_of::<RecordBuf>();
                key_size + vec_size + vec_overhead
            })
            .sum();
        let hashmap_overhead = self.mi_groups.capacity() * 48; // rough estimate
        mi_groups_size + hashmap_overhead
    }
}

/// A batch of position-MI groups for parallel processing.
#[derive(Default)]
pub struct PositionMiGroupBatch {
    /// The position-MI groups in this batch.
    pub groups: Vec<PositionMiGroup>,
}

impl PositionMiGroupBatch {
    /// Creates a new empty batch.
    #[must_use]
    pub fn new() -> Self {
        Self { groups: Vec::new() }
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
}

impl BatchWeight for PositionMiGroupBatch {
    fn batch_weight(&self) -> usize {
        self.groups.iter().map(|g| g.total_templates).sum()
    }
}

impl MemoryEstimate for PositionMiGroupBatch {
    fn estimate_heap_size(&self) -> usize {
        let groups_size: usize = self.groups.iter().map(MemoryEstimate::estimate_heap_size).sum();
        let groups_vec_overhead = self.groups.capacity() * std::mem::size_of::<PositionMiGroup>();
        groups_size + groups_vec_overhead
    }
}

/// A grouper that groups records by position, then by MI.
///
/// Used by the `duplex` command with `--metrics-output` to collect duplex metrics
/// during consensus calling. Groups records first by genomic position using
/// `ReadInfoKey` (matching duplex-metrics exactly), then organizes them by base MI.
///
/// Expects input sorted by template-coordinate with MI tags already assigned
/// (output from `group` command).
pub struct PositionMiGrouper {
    /// MI tag name (e.g., "MI").
    mi_tag: Tag,
    /// Number of position groups per batch.
    batch_size: usize,

    // Template accumulation (group records by QNAME first)
    /// Current template name hash being accumulated.
    current_name_hash: Option<u64>,
    /// Records for the current template.
    current_template_records: Vec<RecordBuf>,

    // Position group (using ReadInfoKey computed from template)
    /// Current ReadInfoKey for the position group.
    current_read_info_key: Option<ReadInfoKey>,
    /// Records at current position, organized by base MI.
    current_mi_groups: HashMap<String, Vec<RecordBuf>>,
    /// Number of templates at current position (for CS family count).
    current_template_count: usize,
    /// Total records at current position.
    current_record_count: usize,

    /// Completed groups waiting to be batched.
    pending_groups: VecDeque<PositionMiGroup>,
    /// Whether `finish()` has been called.
    finished: bool,
}

impl PositionMiGrouper {
    /// Create a new `PositionMiGrouper`.
    ///
    /// # Arguments
    /// * `mi_tag_name` - The MI tag name (e.g., "MI")
    /// * `batch_size` - Number of position groups per batch
    ///
    /// # Panics
    /// Panics if `mi_tag_name` is not exactly 2 characters.
    #[must_use]
    pub fn new(mi_tag_name: &str, batch_size: usize) -> Self {
        assert!(mi_tag_name.len() == 2, "Tag name must be exactly 2 characters");
        let tag_bytes = mi_tag_name.as_bytes();
        let mi_tag = Tag::from([tag_bytes[0], tag_bytes[1]]);

        Self {
            mi_tag,
            batch_size: batch_size.max(1),
            current_name_hash: None,
            current_template_records: Vec::new(),
            current_read_info_key: None,
            current_mi_groups: HashMap::new(),
            current_template_count: 0,
            current_record_count: 0,
            pending_groups: VecDeque::new(),
            finished: false,
        }
    }

    /// Extract the MI tag value from a record.
    fn get_record_mi(&self, record: &RecordBuf) -> Option<String> {
        if let Some(val) = record.data().get(&self.mi_tag) {
            use noodles::sam::alignment::record_buf::data::field::Value;
            if let Value::String(s) = val {
                return Some(s.to_string());
            }
        }
        None
    }


    /// Process accumulated template: compute `ReadInfoKey` and add to position group.
    fn flush_template(&mut self) {
        if self.current_template_records.is_empty() {
            return;
        }

        let records = std::mem::take(&mut self.current_template_records);
        let read_info_key = ReadInfoKey::from_template_records(&records);
        self.current_name_hash = None;

        // Skip templates with invalid ReadInfoKey (unpaired, unmapped, etc.)
        if !read_info_key.is_valid() {
            return;
        }

        // Check if position changed
        if self.current_read_info_key.as_ref() != Some(&read_info_key) {
            self.flush_current_position_group();
            self.current_read_info_key = Some(read_info_key);
        }

        // Extract MI from first record with MI tag
        let full_mi = records.iter().find_map(|r| self.get_record_mi(r)).unwrap_or_default();
        let base_mi = extract_mi_base(&full_mi).to_string();

        // Add records to MI group
        let record_count = records.len();
        self.current_mi_groups.entry(base_mi).or_default().extend(records);
        self.current_template_count += 1;
        self.current_record_count += record_count;
    }

    /// Flush current position group to pending queue.
    fn flush_current_position_group(&mut self) {
        if let Some(key) = self.current_read_info_key.take() {
            if !self.current_mi_groups.is_empty() || self.current_record_count > 0 {
                let mi_groups = std::mem::take(&mut self.current_mi_groups);
                let total_templates = self.current_template_count;
                let total_records = self.current_record_count;
                self.current_template_count = 0;
                self.current_record_count = 0;

                self.pending_groups.push_back(PositionMiGroup {
                    read_info_key: key,
                    mi_groups,
                    total_templates,
                    total_records,
                });
            }
        }
    }

    /// Try to form complete batches from pending groups.
    fn drain_batches(&mut self) -> Vec<PositionMiGroupBatch> {
        let mut batches = Vec::new();
        while self.pending_groups.len() >= self.batch_size {
            let groups: Vec<PositionMiGroup> =
                self.pending_groups.drain(..self.batch_size).collect();
            batches.push(PositionMiGroupBatch { groups });
        }
        batches
    }
}

impl Grouper for PositionMiGrouper {
    type Group = PositionMiGroupBatch;

    fn add_records(&mut self, records: Vec<DecodedRecord>) -> io::Result<Vec<Self::Group>> {
        for decoded in records {
            let record = decoded.record;
            let key = decoded.key;

            // Skip secondary and supplementary reads
            let flags = record.flags();
            if flags.is_secondary() || flags.is_supplementary() {
                continue;
            }

            // Check if same template (by name hash)
            if self.current_name_hash != Some(key.name_hash) {
                // New template - flush previous
                self.flush_template();
                self.current_name_hash = Some(key.name_hash);
            }

            self.current_template_records.push(record);
        }

        Ok(self.drain_batches())
    }

    fn finish(&mut self) -> io::Result<Option<Self::Group>> {
        if self.finished {
            return Ok(None);
        }
        self.finished = true;

        // Flush any remaining template
        self.flush_template();

        // Flush any remaining position group
        self.flush_current_position_group();

        // Return remaining groups as final batch
        if self.pending_groups.is_empty() {
            Ok(None)
        } else {
            let groups: Vec<PositionMiGroup> = self.pending_groups.drain(..).collect();
            Ok(Some(PositionMiGroupBatch { groups }))
        }
    }

    fn has_pending(&self) -> bool {
        !self.current_template_records.is_empty()
            || !self.current_mi_groups.is_empty()
            || !self.pending_groups.is_empty()
    }
}

// ============================================================================
// FASTQ Record Parsing Utilities
// ============================================================================

/// A parsed FASTQ record.
#[derive(Debug, Clone)]
pub struct FastqRecord {
    /// Read name (without @ prefix).
    pub name: Vec<u8>,
    /// Sequence bases.
    pub sequence: Vec<u8>,
    /// Quality scores (Phred+33 encoded).
    pub quality: Vec<u8>,
}

/// Result of parsing a FASTQ record.
#[derive(Debug)]
enum FastqParseResult {
    /// Record incomplete - need more data.
    Incomplete,
    /// Parse error.
    Error(io::Error),
}

/// Parse FASTQ records from a byte buffer.
///
/// FASTQ format:
/// ```text
/// @read_name
/// SEQUENCE
/// +
/// QUALITY
/// ```
///
/// Returns parsed records and leftover bytes for incomplete records.
pub fn parse_fastq_records(data: &[u8]) -> io::Result<(Vec<FastqRecord>, Vec<u8>)> {
    let mut records = Vec::new();
    let mut pos = 0;

    while pos < data.len() {
        // Find the start of a record (@ character at start of line)
        if data[pos] != b'@' {
            // Skip until we find @ at start of line or run out of data
            while pos < data.len() && data[pos] != b'@' {
                pos += 1;
            }
            if pos >= data.len() {
                return Ok((records, Vec::new()));
            }
        }

        // Try to parse a complete record
        match parse_single_fastq_record(&data[pos..]) {
            Ok((record, consumed)) => {
                records.push(record);
                pos += consumed;
            }
            Err(FastqParseResult::Incomplete) => {
                // Not enough data for a complete record
                return Ok((records, data[pos..].to_vec()));
            }
            Err(FastqParseResult::Error(e)) => {
                return Err(e);
            }
        }
    }

    Ok((records, Vec::new()))
}

/// Parse a single FASTQ record from the beginning of a buffer.
/// Returns (record, `bytes_consumed`) or an error.
fn parse_single_fastq_record(data: &[u8]) -> Result<(FastqRecord, usize), FastqParseResult> {
    let mut pos = 0;

    // Line 1: @name
    if data.is_empty() || data[0] != b'@' {
        return Err(FastqParseResult::Error(io::Error::new(
            io::ErrorKind::InvalidData,
            "FASTQ record must start with @",
        )));
    }
    pos += 1; // Skip @

    let name_end = find_newline(&data[pos..]).ok_or(FastqParseResult::Incomplete)?;
    let name = data[pos..pos + name_end].to_vec();
    pos += name_end + 1; // +1 for newline

    // Line 2: sequence
    if pos >= data.len() {
        return Err(FastqParseResult::Incomplete);
    }
    let seq_end = find_newline(&data[pos..]).ok_or(FastqParseResult::Incomplete)?;
    let sequence = data[pos..pos + seq_end].to_vec();
    pos += seq_end + 1;

    // Line 3: + (separator)
    if pos >= data.len() {
        return Err(FastqParseResult::Incomplete);
    }
    let plus_end = find_newline(&data[pos..]).ok_or(FastqParseResult::Incomplete)?;
    if data[pos] != b'+' {
        return Err(FastqParseResult::Error(io::Error::new(
            io::ErrorKind::InvalidData,
            "FASTQ separator line must start with +",
        )));
    }
    pos += plus_end + 1;

    // Line 4: quality
    if pos >= data.len() {
        return Err(FastqParseResult::Incomplete);
    }
    let qual_end = find_newline(&data[pos..]).ok_or(FastqParseResult::Incomplete)?;
    let quality = data[pos..pos + qual_end].to_vec();
    pos += qual_end + 1;

    // Validate lengths match
    if sequence.len() != quality.len() {
        return Err(FastqParseResult::Error(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Sequence length ({}) != quality length ({})", sequence.len(), quality.len()),
        )));
    }

    Ok((FastqRecord { name, sequence, quality }, pos))
}

/// Find the position of the next newline character.
fn find_newline(data: &[u8]) -> Option<usize> {
    data.iter().position(|&b| b == b'\n')
}

/// Strip /1 or /2 suffix from read name for comparison.
#[must_use]
pub fn strip_read_suffix(name: &[u8]) -> &[u8] {
    if name.len() >= 2 {
        let suffix = &name[name.len() - 2..];
        if suffix == b"/1" || suffix == b"/2" {
            return &name[..name.len() - 2];
        }
    }
    name
}

// ============================================================================
// FastqGrouper - Groups FASTQ records from multiple input streams
// ============================================================================

/// A template containing FASTQ records from all input files.
/// For paired-end: R1 + R2 records with matching names.
#[derive(Debug)]
pub struct FastqTemplate {
    /// Records from each input file (index matches input order).
    pub records: Vec<FastqRecord>,
    /// The common read name (should match across all records).
    pub name: Vec<u8>,
}

impl MemoryEstimate for FastqRecord {
    fn estimate_heap_size(&self) -> usize {
        self.name.capacity() + self.sequence.capacity() + self.quality.capacity()
    }
}

impl MemoryEstimate for FastqTemplate {
    fn estimate_heap_size(&self) -> usize {
        let records_heap: usize = self.records.iter().map(MemoryEstimate::estimate_heap_size).sum();
        let records_vec_overhead = self.records.capacity() * std::mem::size_of::<FastqRecord>();
        self.name.capacity() + records_heap + records_vec_overhead
    }
}

/// Groups FASTQ records from multiple synchronized input streams into templates.
///
/// This grouper expects decompressed bytes from multiple FASTQ files,
/// provided via `add_bytes_for_stream`. It parses records and groups them by
/// read name into templates.
pub struct FastqGrouper {
    /// Number of input streams.
    num_inputs: usize,
    /// Leftover bytes for each input stream.
    leftovers: Vec<Vec<u8>>,
    /// Parsed but not yet grouped records for each stream.
    pending_records: Vec<VecDeque<FastqRecord>>,
}

impl FastqGrouper {
    /// Create a new FASTQ grouper for the specified number of input streams.
    #[must_use]
    pub fn new(num_inputs: usize) -> Self {
        log::debug!("FastqGrouper::new: num_inputs={num_inputs}");
        Self {
            num_inputs,
            leftovers: vec![Vec::new(); num_inputs],
            pending_records: (0..num_inputs).map(|_| VecDeque::new()).collect(),
        }
    }

    /// Add decompressed bytes from a specific input stream.
    ///
    /// **Note:** This method parses the bytes into records under the caller's lock.
    /// For better parallel scaling, use `add_records_for_stream` with pre-parsed records.
    pub fn add_bytes_for_stream(&mut self, stream_idx: usize, data: &[u8]) -> io::Result<()> {
        log::trace!(
            "FastqGrouper::add_bytes_for_stream: stream_idx={}, num_inputs={}, data_len={}",
            stream_idx,
            self.num_inputs,
            data.len()
        );
        if stream_idx >= self.num_inputs {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("Stream index {} out of range (max {})", stream_idx, self.num_inputs - 1),
            ));
        }

        // Combine with leftover
        let combined = if self.leftovers[stream_idx].is_empty() {
            data.to_vec()
        } else {
            let mut combined = std::mem::take(&mut self.leftovers[stream_idx]);
            combined.extend_from_slice(data);
            combined
        };

        // Parse FASTQ records
        let (records, leftover) = parse_fastq_records(&combined)?;
        self.leftovers[stream_idx] = leftover;

        // Add to pending
        self.pending_records[stream_idx].extend(records);

        Ok(())
    }

    /// Add pre-parsed records directly to a specific input stream.
    ///
    /// This method enables the parallel Parse optimization by accepting already-parsed
    /// records instead of raw bytes. The parsing happens in a parallel step, while
    /// grouping (which requires sequential access) only does the grouping work.
    ///
    /// **This is the key method for fixing the t8 scaling bottleneck.**
    pub fn add_records_for_stream(
        &mut self,
        stream_idx: usize,
        records: Vec<FastqRecord>,
    ) -> io::Result<()> {
        log::trace!(
            "FastqGrouper::add_records_for_stream: stream_idx={}, num_inputs={}, records_len={}",
            stream_idx,
            self.num_inputs,
            records.len()
        );
        if stream_idx >= self.num_inputs {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("Stream index {} out of range (max {})", stream_idx, self.num_inputs - 1),
            ));
        }

        // Add directly to pending - no parsing needed!
        self.pending_records[stream_idx].extend(records);

        Ok(())
    }

    /// Check if there are leftover bytes for any stream.
    ///
    /// When using the parallel Parse optimization, this should always return false
    /// because boundary finding handles leftovers separately.
    #[must_use]
    pub fn has_leftover_bytes(&self) -> bool {
        self.leftovers.iter().any(|l| !l.is_empty())
    }

    /// Try to emit complete templates (records from all streams with matching names).
    pub fn drain_complete_templates(&mut self) -> io::Result<Vec<FastqTemplate>> {
        let mut templates = Vec::new();

        // Keep emitting while we have at least one record in each stream
        while self.pending_records.iter().all(|q| !q.is_empty()) {
            // Validate names match and get base_name (in block so names is dropped before pop)
            let base_name = {
                // Peek at the first record from each stream
                let names: Vec<_> =
                    self.pending_records.iter().map(|q| &q.front().unwrap().name).collect();

                // Copy base_name immediately
                let base_name = strip_read_suffix(names[0]).to_vec();

                // Validate all names match (strip /1, /2 suffixes for comparison)
                for (i, name) in names.iter().enumerate().skip(1) {
                    let other_base = strip_read_suffix(name);
                    if base_name != other_base {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!(
                                "FASTQ files out of sync: stream 0 has '{}', stream {} has '{}'",
                                String::from_utf8_lossy(&base_name),
                                i,
                                String::from_utf8_lossy(other_base),
                            ),
                        ));
                    }
                }

                base_name // names dropped here
            };

            // Pop records from all streams
            let records: Vec<_> =
                self.pending_records.iter_mut().map(|q| q.pop_front().unwrap()).collect();

            templates.push(FastqTemplate { name: base_name, records });
        }

        Ok(templates)
    }

    /// Check if there are pending records or leftover bytes.
    #[must_use]
    pub fn has_pending(&self) -> bool {
        self.leftovers.iter().any(|l| !l.is_empty())
            || self.pending_records.iter().any(|q| !q.is_empty())
    }

    /// Finish processing and return any remaining template.
    pub fn finish(&mut self) -> io::Result<Option<FastqTemplate>> {
        // Check for any remaining complete templates
        let templates = self.drain_complete_templates()?;
        if templates.len() == 1 {
            return Ok(templates.into_iter().next());
        } else if templates.len() > 1 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Multiple templates remaining at finish",
            ));
        }

        // Check for incomplete data
        if self.leftovers.iter().any(|l| !l.is_empty()) {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "Incomplete FASTQ record at EOF",
            ));
        }

        if self.pending_records.iter().any(|q| !q.is_empty()) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Unmatched FASTQ records at EOF - files out of sync",
            ));
        }

        Ok(None)
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_single_record_grouper_empty() {
        let mut grouper = SingleRecordGrouper::new();
        assert!(!grouper.has_pending());

        let result = grouper.finish().unwrap();
        assert!(result.is_none());
    }

    // FASTQ parsing tests
    #[test]
    fn test_parse_single_fastq_record() {
        let data = b"@read1\nACGT\n+\nIIII\n";
        let (record, consumed) = parse_single_fastq_record(data).unwrap();
        assert_eq!(record.name, b"read1");
        assert_eq!(record.sequence, b"ACGT");
        assert_eq!(record.quality, b"IIII");
        assert_eq!(consumed, data.len());
    }

    #[test]
    fn test_parse_fastq_records_multiple() {
        let data = b"@read1\nACGT\n+\nIIII\n@read2\nTGCA\n+\nJJJJ\n";
        let (records, leftover) = parse_fastq_records(data).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name, b"read1");
        assert_eq!(records[1].name, b"read2");
        assert!(leftover.is_empty());
    }

    #[test]
    fn test_parse_fastq_incomplete_record() {
        let data = b"@read1\nACGT\n+\n";
        let (records, leftover) = parse_fastq_records(data).unwrap();
        assert!(records.is_empty());
        assert_eq!(leftover, data);
    }

    #[test]
    fn test_parse_fastq_with_leftover() {
        let data = b"@read1\nACGT\n+\nIIII\n@read2\nTG";
        let (records, leftover) = parse_fastq_records(data).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].name, b"read1");
        assert_eq!(leftover, b"@read2\nTG");
    }

    #[test]
    fn test_strip_read_suffix() {
        assert_eq!(strip_read_suffix(b"read1/1"), b"read1");
        assert_eq!(strip_read_suffix(b"read1/2"), b"read1");
        assert_eq!(strip_read_suffix(b"read1"), b"read1");
        assert_eq!(strip_read_suffix(b"a"), b"a");
        assert_eq!(strip_read_suffix(b""), b"" as &[u8]);
    }

    // FastqGrouper tests
    #[test]
    fn test_fastq_grouper_paired() {
        let mut grouper = FastqGrouper::new(2);

        // Add R1 record
        grouper.add_bytes_for_stream(0, b"@read1/1\nACGT\n+\nIIII\n").unwrap();
        // Add R2 record
        grouper.add_bytes_for_stream(1, b"@read1/2\nTGCA\n+\nJJJJ\n").unwrap();

        let templates = grouper.drain_complete_templates().unwrap();
        assert_eq!(templates.len(), 1);
        assert_eq!(templates[0].name, b"read1");
        assert_eq!(templates[0].records.len(), 2);
        assert_eq!(templates[0].records[0].sequence, b"ACGT");
        assert_eq!(templates[0].records[1].sequence, b"TGCA");
    }

    #[test]
    fn test_fastq_grouper_multiple_templates() {
        let mut grouper = FastqGrouper::new(2);

        // Add multiple R1 records
        grouper
            .add_bytes_for_stream(0, b"@read1/1\nACGT\n+\nIIII\n@read2/1\nAAAA\n+\nIIII\n")
            .unwrap();
        // Add multiple R2 records
        grouper
            .add_bytes_for_stream(1, b"@read1/2\nTGCA\n+\nJJJJ\n@read2/2\nTTTT\n+\nJJJJ\n")
            .unwrap();

        let templates = grouper.drain_complete_templates().unwrap();
        assert_eq!(templates.len(), 2);
        assert_eq!(templates[0].name, b"read1");
        assert_eq!(templates[1].name, b"read2");
    }

    #[test]
    fn test_fastq_grouper_partial_then_complete() {
        let mut grouper = FastqGrouper::new(2);

        // Add R1 record only
        grouper.add_bytes_for_stream(0, b"@read1/1\nACGT\n+\nIIII\n").unwrap();

        // No complete templates yet
        let templates = grouper.drain_complete_templates().unwrap();
        assert!(templates.is_empty());
        assert!(grouper.has_pending());

        // Now add R2 record
        grouper.add_bytes_for_stream(1, b"@read1/2\nTGCA\n+\nJJJJ\n").unwrap();

        let templates = grouper.drain_complete_templates().unwrap();
        assert_eq!(templates.len(), 1);
    }

    #[test]
    fn test_fastq_grouper_finish_empty() {
        let mut grouper = FastqGrouper::new(2);
        let result = grouper.finish().unwrap();
        assert!(result.is_none());
    }

    #[test]
    fn test_fastq_grouper_out_of_sync() {
        let mut grouper = FastqGrouper::new(2);

        // Add mismatched records
        grouper.add_bytes_for_stream(0, b"@read1/1\nACGT\n+\nIIII\n").unwrap();
        grouper.add_bytes_for_stream(1, b"@read2/2\nTGCA\n+\nJJJJ\n").unwrap();

        let result = grouper.drain_complete_templates();
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("out of sync"));
    }

    // ========================================================================
    // combine_keys tests
    // ========================================================================

    #[test]
    fn test_combine_keys_empty() {
        let keys: Vec<GroupKey> = vec![];
        let result = PositionGrouper::combine_keys(&keys);
        assert_eq!(result.ref_id1, GroupKey::UNKNOWN_REF);
    }

    #[test]
    fn test_combine_keys_single_valid() {
        let key = GroupKey::single(5, 1000, 0, 0, 0, 12345);
        let keys = vec![key];
        let result = PositionGrouper::combine_keys(&keys);
        assert_eq!(result.ref_id1, 5);
        assert_eq!(result.pos1, 1000);
        assert_eq!(result.name_hash, 12345);
    }

    #[test]
    fn test_combine_keys_single_unknown() {
        // Secondary/supplementary reads have UNKNOWN ref_id1
        let key = GroupKey { name_hash: 99999, ..GroupKey::default() };
        let keys = vec![key];
        let result = PositionGrouper::combine_keys(&keys);
        // Should return default with preserved name_hash
        assert_eq!(result.ref_id1, GroupKey::UNKNOWN_REF);
        assert_eq!(result.name_hash, 99999);
    }

    #[test]
    fn test_combine_keys_two_valid_without_mate_info() {
        // R1 and R2 each have their own position but no mate info
        let r1_key = GroupKey::single(0, 100, 0, 0, 0, 12345);
        let r2_key = GroupKey::single(0, 300, 1, 0, 0, 12345);
        let keys = vec![r1_key, r2_key];
        let result = PositionGrouper::combine_keys(&keys);

        // Should combine: r1's position as pos1, r2's position as pos2
        assert_eq!(result.ref_id1, 0);
        assert_eq!(result.pos1, 100);
        assert_eq!(result.strand1, 0);
        assert_eq!(result.ref_id2, 0);
        assert_eq!(result.pos2, 300);
        assert_eq!(result.strand2, 1);
    }

    #[test]
    fn test_combine_keys_key_with_mate_info() {
        // Key that already has mate info should be returned directly
        let key = GroupKey::paired(0, 100, 0, 0, 300, 1, 0, 0, 12345);
        let keys = vec![key];
        let result = PositionGrouper::combine_keys(&keys);

        assert_eq!(result.ref_id1, 0);
        assert_eq!(result.pos1, 100);
        assert_eq!(result.ref_id2, 0);
        assert_eq!(result.pos2, 300);
    }

    #[test]
    fn test_combine_keys_mixed_valid_and_unknown() {
        // Mix of valid key (primary) and UNKNOWN key (secondary/supplementary)
        let valid_key = GroupKey::single(5, 500, 0, 0, 0, 12345);
        let unknown_key = GroupKey { name_hash: 12345, ..GroupKey::default() };
        let keys = vec![unknown_key, valid_key];
        let result = PositionGrouper::combine_keys(&keys);

        // Should skip unknown and use valid key
        assert_eq!(result.ref_id1, 5);
        assert_eq!(result.pos1, 500);
    }

    #[test]
    fn test_combine_keys_all_unknown() {
        // All secondary/supplementary - no valid primary
        let unknown_key1 = GroupKey { name_hash: 12345, ..GroupKey::default() };
        let unknown_key2 = GroupKey { name_hash: 12345, ..GroupKey::default() };
        let keys = vec![unknown_key1, unknown_key2];
        let result = PositionGrouper::combine_keys(&keys);

        // Should return default with preserved name_hash
        assert_eq!(result.ref_id1, GroupKey::UNKNOWN_REF);
        assert_eq!(result.name_hash, 12345);
    }

    #[test]
    fn test_combine_keys_first_has_mate_info() {
        // First valid key has mate info - should return it directly
        let key_with_mate = GroupKey::paired(0, 100, 0, 0, 200, 1, 0, 0, 12345);
        let key_without_mate = GroupKey::single(0, 100, 0, 0, 0, 12345);
        let keys = vec![key_with_mate, key_without_mate];
        let result = PositionGrouper::combine_keys(&keys);

        assert_eq!(result.ref_id2, 0);
        assert_eq!(result.pos2, 200);
    }
}

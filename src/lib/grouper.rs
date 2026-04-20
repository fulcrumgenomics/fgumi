//! Grouper implementations for the 9-step pipeline.
//!
//! Each grouper is responsible for grouping decoded BAM records according
//! to command-specific rules and emitting complete groups for processing.
//!
//! Note: BAM record parsing is handled by the Decode step in the pipeline.
//! Groupers receive pre-decoded raw-byte records.

use std::io;

use noodles::sam::alignment::RecordBuf;

use crate::sort::bam_fields;
use crate::template::{Template, TemplateBatch};
use crate::unified_pipeline::{BatchWeight, DecodedRecord, Grouper, MemoryEstimate};
use fgumi_raw_bam::{RawRecord, raw_record_to_record_buf};

// ============================================================================
// BatchWeight Implementations
// ============================================================================

/// Single records have weight 1 (for batch_size-like behavior).
impl BatchWeight for RecordBuf {
    fn batch_weight(&self) -> usize {
        1
    }
}

/// Raw byte records have weight 1 (for batch_size-like behavior).
impl BatchWeight for Vec<u8> {
    fn batch_weight(&self) -> usize {
        1
    }
}

/// [`RawRecord`] has weight 1 (same semantics as the `Vec<u8>` impl above).
impl BatchWeight for RawRecord {
    fn batch_weight(&self) -> usize {
        1
    }
}

/// Heap size of a [`RawRecord`] is its buffer capacity.
impl MemoryEstimate for RawRecord {
    fn estimate_heap_size(&self) -> usize {
        self.capacity()
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

/// A grouper that emits each record as its own "group" (decoded to `RecordBuf`).
///
/// Used by pass-through pipeline tests. Production filter/clip/correct commands
/// use [`SingleRawRecordGrouper`] instead to avoid noodles decode/encode.
#[derive(Default)]
pub struct SingleRecordGrouper {
    header: noodles::sam::Header,
}

impl SingleRecordGrouper {
    /// Create a new single-record grouper with an empty header. Only safe for
    /// unmapped records; mapped records require [`Self::with_header`] so
    /// reference-sequence IDs can be resolved during decoding.
    #[must_use]
    pub fn new() -> Self {
        Self { header: noodles::sam::Header::default() }
    }

    /// Create a single-record grouper that decodes against `header`. Required
    /// when records are mapped so reference-sequence IDs resolve correctly.
    #[must_use]
    pub fn with_header(header: noodles::sam::Header) -> Self {
        Self { header }
    }
}

impl Grouper for SingleRecordGrouper {
    type Group = RecordBuf;

    fn add_records(&mut self, records: Vec<DecodedRecord>) -> io::Result<Vec<Self::Group>> {
        records
            .into_iter()
            .map(|d| {
                raw_record_to_record_buf(&d.into_raw_bytes(), &self.header)
                    .map_err(io::Error::other)
            })
            .collect()
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
// SingleRawRecordGrouper
// ============================================================================

/// A grouper that emits each raw-byte record as its own "group".
///
/// Used by commands that process records independently using raw bytes
/// (e.g., filter with raw-byte pipeline).
#[derive(Default)]
pub struct SingleRawRecordGrouper;

impl SingleRawRecordGrouper {
    /// Create a new single raw record grouper.
    #[must_use]
    pub fn new() -> Self {
        Self
    }
}

impl Grouper for SingleRawRecordGrouper {
    type Group = RawRecord;

    fn add_records(&mut self, records: Vec<DecodedRecord>) -> io::Result<Vec<Self::Group>> {
        Ok(records.into_iter().map(DecodedRecord::into_raw_bytes).collect())
    }

    fn finish(&mut self) -> io::Result<Option<Self::Group>> {
        Ok(None)
    }

    fn has_pending(&self) -> bool {
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
    /// Raw records for the current template being built.
    current_raw_records: Vec<RawRecord>,
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
            current_raw_records: Vec::new(),
            pending_templates: VecDeque::new(),
        }
    }

    /// Flush current template to pending queue if non-empty.
    fn flush_current_template(&mut self) -> io::Result<()> {
        if !self.current_raw_records.is_empty() {
            let template = Template::from_records(std::mem::take(&mut self.current_raw_records))
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
        // Group records by QNAME: use name_hash as a fast pre-check, then confirm
        // the QNAME bytes match to defend against hash collisions merging unrelated
        // records into the same template.
        for decoded in records {
            let name_hash = decoded.key.name_hash;
            let raw = decoded.data;
            let read_name = bam_fields::read_name(&raw);
            let same_template = match (self.current_name_hash, self.current_name.as_deref()) {
                (Some(h), Some(name)) => h == name_hash && name == read_name,
                _ => false,
            };
            if same_template {
                self.current_raw_records.push(raw);
            } else {
                self.flush_current_template()?;
                self.current_name = Some(read_name.to_vec());
                self.current_name_hash = Some(name_hash);
                self.current_raw_records.push(raw);
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
        self.current_name.is_some()
            || !self.pending_templates.is_empty()
            || !self.current_raw_records.is_empty()
    }
}

use crate::unified_pipeline::GroupKey;

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
    /// `Template.mi` contains local IDs (0, 1, 2, ...) - the serialize step adds a global offset.
    pub templates: Vec<Template>,
    /// Family size counts for this position group.
    pub family_sizes: ahash::AHashMap<usize, u64>,
    /// Filter metrics for this position group (thread-local, merged later).
    pub filter_metrics: FilterMetrics,
    /// Total input records processed (for progress tracking).
    pub input_record_count: u64,
    /// Number of distinct numeric molecule IDs assigned in this group (i.e., the
    /// size of the block the serialize step must reserve from the global MI
    /// counter). Equals `max(MoleculeId::id()) + 1` across assigned templates,
    /// or `0` when no templates have an assigned MI. This can be strictly less
    /// than `templates.len()` because templates in the same UMI family share
    /// a `MoleculeId`, and because `PairedA(id)` / `PairedB(id)` share the same
    /// numeric `id`. Using this value for the block size (rather than the
    /// template count) ensures the emitted MI tag integers are consecutive
    /// `0..N-1`, matching fgbio's `GroupReadsByUmi`.
    pub distinct_mi_count: u64,
}

// ============================================================================
// RawPositionGroup + RecordPositionGrouper
// ============================================================================

/// Position key tuple returned by [`GroupKey::position_key`].
type PositionKeyTuple = (i32, i32, u8, i32, i32, u8, u16, u64);

/// A position group containing decoded records at the same genomic position.
///
/// Contains raw [`DecodedRecord`]s rather than built [`Template`]s. Template
/// construction is deferred to the parallel Process step via
/// [`build_templates_from_records`].
///
/// Used by [`RecordPositionGrouper`] which requires MC tags so that each record's
/// pre-computed [`GroupKey`] contains the complete paired position.
#[derive(Debug)]
pub struct RawPositionGroup {
    /// The position key for this group (from the first record's `GroupKey`).
    pub group_key: GroupKey,
    /// Raw decoded records at this position (not yet grouped into Templates).
    pub records: Vec<DecodedRecord>,
}

impl BatchWeight for RawPositionGroup {
    fn batch_weight(&self) -> usize {
        // Estimate ~2 records per template for paired-end data
        self.records.len().div_ceil(2)
    }
}

impl MemoryEstimate for RawPositionGroup {
    fn estimate_heap_size(&self) -> usize {
        let records_size: usize = self.records.iter().map(MemoryEstimate::estimate_heap_size).sum();
        let records_vec_overhead = self.records.capacity() * std::mem::size_of::<DecodedRecord>();
        records_size + records_vec_overhead
    }
}

/// A lightweight position grouper that compares per-record [`GroupKey::position_key`]
/// values to detect group boundaries.
///
/// Does **not** build Templates or combine keys — it accumulates raw
/// [`DecodedRecord`]s and emits [`RawPositionGroup`]s. Template construction
/// is deferred to the parallel Process step via [`build_templates_from_records`].
///
/// **Requirement:** Paired-end reads must have MC tags so that [`compute_group_key`]
/// produces complete [`GroupKey::paired`] values. Without MC tags, R1 and R2 would
/// get different `position_key()` values and be incorrectly split.
///
/// By default, secondary/supplementary reads are skipped (they have UNKNOWN
/// position keys). Use [`with_secondary_supplementary`](Self::with_secondary_supplementary)
/// to include them — they are coalesced by `name_hash` into the group of their
/// adjacent primary read (requires template-coordinate sorted input).
///
/// [`compute_group_key`]: crate::read_info::compute_group_key
pub struct RecordPositionGrouper {
    /// Current position key being accumulated (tuple for fast comparison).
    current_position_key: Option<PositionKeyTuple>,
    /// Full `GroupKey` for emitting with the group.
    current_group_key: Option<GroupKey>,
    /// Records at the current position.
    current_records: Vec<DecodedRecord>,
    /// Whether MC tags have been validated on paired records.
    mc_validated: bool,
    /// Whether to include secondary and supplementary reads in groups.
    /// When true, secondary/supplementary reads (which have UNKNOWN position keys)
    /// are kept and coalesced by `name_hash` into the group of their primary read.
    /// When false (default), they are skipped.
    include_secondary_supplementary: bool,
}

impl RecordPositionGrouper {
    /// Create a new record-level position grouper.
    ///
    /// Secondary/supplementary reads are skipped by default.
    #[must_use]
    pub fn new() -> Self {
        Self {
            current_position_key: None,
            current_group_key: None,
            current_records: Vec::new(),
            mc_validated: false,
            include_secondary_supplementary: false,
        }
    }

    /// Create a record-level position grouper that includes secondary/supplementary reads.
    ///
    /// Secondary/supplementary reads have UNKNOWN position keys and are coalesced
    /// by `name_hash` into the group of their adjacent primary read.
    #[must_use]
    pub fn with_secondary_supplementary() -> Self {
        Self { include_secondary_supplementary: true, ..Self::new() }
    }

    /// Validate that a paired primary record has an MC tag.
    ///
    /// Skips validation for records that are unmapped or whose mates are unmapped,
    /// since unmapped reads have no CIGAR to report in an MC tag. Every eligible
    /// paired primary record is checked — `mc_validated` is kept only as a
    /// "seen at least one valid MC" marker, not a short-circuit.
    fn validate_mc_tag(&mut self, decoded: &DecodedRecord) -> io::Result<()> {
        use fgumi_raw_bam::RawRecordView;

        let raw = &decoded.data;
        let flg = RawRecordView::new(raw).flags();
        let is_paired = (flg & bam_fields::flags::PAIRED) != 0;
        let is_secondary = (flg & bam_fields::flags::SECONDARY) != 0;
        let is_supplementary = (flg & bam_fields::flags::SUPPLEMENTARY) != 0;
        let is_unmapped = (flg & bam_fields::flags::UNMAPPED) != 0;
        let is_mate_unmapped = (flg & bam_fields::flags::MATE_UNMAPPED) != 0;

        if is_paired && !is_secondary && !is_supplementary && !is_unmapped && !is_mate_unmapped {
            if bam_fields::find_mc_tag_in_record(raw).is_none() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "RecordPositionGrouper requires MC tags on paired-end reads. \
                     Run `fgumi zipper` to add MC tags before `fgumi group`.",
                ));
            }
            self.mc_validated = true;
        }

        Ok(())
    }

    /// Feed a single record to the grouper, returning a completed group if the record
    /// starts a new position. Equivalent to `add_records(vec![decoded])` but avoids the
    /// per-call allocation — useful in streaming hot loops.
    ///
    /// # Errors
    /// Returns an error if MC-tag validation fails on a paired primary record.
    pub fn add_record(&mut self, decoded: DecodedRecord) -> io::Result<Option<RawPositionGroup>> {
        self.process_record(decoded)
    }

    /// Process a single decoded record, potentially emitting a completed group.
    fn process_record(&mut self, decoded: DecodedRecord) -> io::Result<Option<RawPositionGroup>> {
        // Skip secondary and supplementary reads (they have UNKNOWN ref_id1)
        // unless configured to include them (for dedup, which needs them in templates).
        if decoded.key.ref_id1 == GroupKey::UNKNOWN_REF && !self.include_secondary_supplementary {
            return Ok(None);
        }

        // Validate MC tag on first paired primary record
        self.validate_mc_tag(&decoded)?;

        let record_pos_key = decoded.key.position_key();

        match self.current_position_key {
            Some(current_pos) if current_pos == record_pos_key => {
                // Same position — accumulate
                self.current_records.push(decoded);
                Ok(None)
            }
            Some(_)
                if self.current_records.last().is_some_and(|last| {
                    // Hash match is a fast pre-check; confirm QNAME bytes to guard
                    // against hash collisions merging unrelated templates.
                    last.key.name_hash == decoded.key.name_hash
                        && bam_fields::read_name(&last.data) == bam_fields::read_name(&decoded.data)
                }) =>
            {
                // Different position but same template (name_hash + QNAME match with
                // previous record). This happens for paired reads with unmapped mates
                // in template-coordinate sorted input: R1 is mapped at some position
                // while R2 is unmapped (position -1:0), but they're adjacent by QNAME.
                // Keep them in the same group so they form a complete template.
                self.current_records.push(decoded);
                Ok(None)
            }
            Some(_) => {
                // Different position — emit current group, start new one
                let finished_records = std::mem::take(&mut self.current_records);
                let finished_key = self
                    .current_group_key
                    .take()
                    .expect("current_group_key set when current_position_key is set");

                let group = RawPositionGroup { group_key: finished_key, records: finished_records };

                // Start new position group
                self.current_position_key = Some(record_pos_key);
                self.current_group_key = Some(decoded.key);
                self.current_records.push(decoded);

                Ok(Some(group))
            }
            None => {
                // First record
                self.current_position_key = Some(record_pos_key);
                self.current_group_key = Some(decoded.key);
                self.current_records.push(decoded);
                Ok(None)
            }
        }
    }
}

impl Default for RecordPositionGrouper {
    fn default() -> Self {
        Self::new()
    }
}

impl Grouper for RecordPositionGrouper {
    type Group = RawPositionGroup;

    fn add_records(&mut self, records: Vec<DecodedRecord>) -> io::Result<Vec<Self::Group>> {
        let mut completed_groups = Vec::new();
        for decoded in records {
            if let Some(group) = self.process_record(decoded)? {
                completed_groups.push(group);
            }
        }
        Ok(completed_groups)
    }

    fn finish(&mut self) -> io::Result<Option<Self::Group>> {
        if !self.current_records.is_empty() {
            debug_assert!(
                self.current_group_key.is_some(),
                "RecordPositionGrouper has {} buffered records but no group key",
                self.current_records.len()
            );
            if let Some(key) = self.current_group_key.take() {
                let records = std::mem::take(&mut self.current_records);
                self.current_position_key = None;
                return Ok(Some(RawPositionGroup { group_key: key, records }));
            }
        }
        Ok(None)
    }

    fn has_pending(&self) -> bool {
        !self.current_records.is_empty()
    }
}

/// Build [`Template`]s from raw decoded records.
///
/// Groups decoded records by `name_hash` and builds a [`Template`] from each group.
///
/// This is a generic helper: `extract` pulls the per-record payload out of each
/// [`DecodedRecord`], and `build` converts a batch of payloads into a [`Template`].
fn group_by_name_and_build<T>(
    records: Vec<DecodedRecord>,
    extract: impl Fn(DecodedRecord) -> io::Result<T>,
    build: impl Fn(Vec<T>) -> anyhow::Result<Template>,
) -> io::Result<Vec<Template>> {
    let mut templates = Vec::new();
    let mut current_name_hash: Option<u64> = None;
    let mut current_name: Option<Vec<u8>> = None;
    let mut current_items: Vec<T> = Vec::new();

    for decoded in records {
        let name_hash = decoded.key.name_hash;
        let read_name = bam_fields::read_name(&decoded.data).to_vec();
        let item = extract(decoded)?;

        // name_hash is a fast pre-check; confirm QNAME bytes to guard against
        // hash collisions merging unrelated templates.
        let same = current_name_hash == Some(name_hash)
            && current_name.as_deref() == Some(read_name.as_slice());
        if same {
            current_items.push(item);
        } else if current_name_hash.is_some() {
            let template = build(std::mem::take(&mut current_items))
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            templates.push(template);
            current_name_hash = Some(name_hash);
            current_name = Some(read_name);
            current_items.push(item);
        } else {
            current_name_hash = Some(name_hash);
            current_name = Some(read_name);
            current_items.push(item);
        }
    }

    if !current_items.is_empty() {
        let template =
            build(current_items).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        templates.push(template);
    }

    Ok(templates)
}

/// Groups records by `name_hash` (QNAME) and builds a [`Template`] from each group.
/// Records must be grouped by QNAME (guaranteed by template-coordinate sort order).
///
/// Designed to run in the parallel Process step after [`RecordPositionGrouper`]
/// emits [`RawPositionGroup`]s.
///
/// # Errors
///
/// Returns an error if template construction from records fails.
pub fn build_templates_from_records(records: Vec<DecodedRecord>) -> io::Result<Vec<Template>> {
    group_by_name_and_build(records, |d| Ok(d.data), Template::from_records)
}

use crate::fastq_parse::{FastqRecord, parse_fastq_records, strip_read_suffix};

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
    ///
    /// # Errors
    ///
    /// Returns an error if the stream index is out of range or parsing fails.
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
    ///
    /// # Errors
    ///
    /// Returns an error if the stream index is out of range.
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
    ///
    /// # Errors
    ///
    /// Returns an error if FASTQ files are out of sync (mismatched read names).
    ///
    /// # Panics
    ///
    /// Panics if a pending record queue is unexpectedly empty after checking non-emptiness.
    pub fn drain_complete_templates(&mut self) -> io::Result<Vec<FastqTemplate>> {
        let mut templates = Vec::new();

        // Keep emitting while we have at least one record in each stream
        while self.pending_records.iter().all(|q| !q.is_empty()) {
            // Validate names match and get base_name (in block so names is dropped before pop)
            let base_name = {
                // Peek at the first record from each stream
                let names: Vec<&[u8]> = self
                    .pending_records
                    .iter()
                    .map(|q| {
                        q.front()
                            .expect("pending queue must be non-empty inside all-non-empty loop")
                            .name()
                    })
                    .collect();

                // Copy base_name immediately
                let base_name = strip_read_suffix(names[0]).to_vec();

                // Validate all names match (strip /1, /2 suffixes for comparison)
                for (i, &name) in names.iter().enumerate().skip(1) {
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
            let records: Vec<_> = self
                .pending_records
                .iter_mut()
                .map(|q| {
                    q.pop_front()
                        .expect("pending queue must be non-empty inside all-non-empty loop")
                })
                .collect();

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
    ///
    /// # Errors
    ///
    /// Returns an error if there is incomplete or unmatched data at EOF.
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

        let result = grouper.finish().expect("finish should succeed");
        assert!(result.is_none());
    }

    #[test]
    fn test_single_raw_record_grouper_empty() {
        let mut grouper = SingleRawRecordGrouper::new();
        assert!(!grouper.has_pending());

        let result = grouper.finish().expect("finish should succeed");
        assert!(result.is_none());
    }

    #[test]
    fn test_single_raw_record_grouper_emits_each_record() {
        use crate::unified_pipeline::{DecodedRecord, GroupKey};

        let mut grouper = SingleRawRecordGrouper::new();
        let raw1 = RawRecord::from(vec![1u8; 36]);
        let raw2 = RawRecord::from(vec![2u8; 36]);
        let records = vec![
            DecodedRecord::from_raw_bytes(raw1.clone(), GroupKey::default()),
            DecodedRecord::from_raw_bytes(raw2.clone(), GroupKey::default()),
        ];

        let groups = grouper.add_records(records).expect("add_records should succeed");
        assert_eq!(groups.len(), 2);
        assert_eq!(groups[0], raw1);
        assert_eq!(groups[1], raw2);
    }

    // FastqGrouper tests
    #[test]
    fn test_fastq_grouper_paired() {
        let mut grouper = FastqGrouper::new(2);

        // Add R1 record
        grouper
            .add_bytes_for_stream(0, b"@read1/1\nACGT\n+\nIIII\n")
            .expect("add_bytes_for_stream failed");
        // Add R2 record
        grouper
            .add_bytes_for_stream(1, b"@read1/2\nTGCA\n+\nJJJJ\n")
            .expect("add_bytes_for_stream failed");

        let templates =
            grouper.drain_complete_templates().expect("drain_complete_templates failed");
        assert_eq!(templates.len(), 1);
        assert_eq!(templates[0].name, b"read1");
        assert_eq!(templates[0].records.len(), 2);
        assert_eq!(templates[0].records[0].sequence(), b"ACGT");
        assert_eq!(templates[0].records[1].sequence(), b"TGCA");
    }

    #[test]
    fn test_fastq_grouper_multiple_templates() {
        let mut grouper = FastqGrouper::new(2);

        // Add multiple R1 records
        grouper
            .add_bytes_for_stream(0, b"@read1/1\nACGT\n+\nIIII\n@read2/1\nAAAA\n+\nIIII\n")
            .expect("add_bytes_for_stream should succeed for R1");
        // Add multiple R2 records
        grouper
            .add_bytes_for_stream(1, b"@read1/2\nTGCA\n+\nJJJJ\n@read2/2\nTTTT\n+\nJJJJ\n")
            .expect("add_bytes_for_stream should succeed for R2");

        let templates =
            grouper.drain_complete_templates().expect("drain_complete_templates failed");
        assert_eq!(templates.len(), 2);
        assert_eq!(templates[0].name, b"read1");
        assert_eq!(templates[1].name, b"read2");
    }

    #[test]
    fn test_fastq_grouper_partial_then_complete() {
        let mut grouper = FastqGrouper::new(2);

        // Add R1 record only
        grouper
            .add_bytes_for_stream(0, b"@read1/1\nACGT\n+\nIIII\n")
            .expect("add_bytes_for_stream failed");

        // No complete templates yet
        let templates =
            grouper.drain_complete_templates().expect("drain_complete_templates failed");
        assert!(templates.is_empty());
        assert!(grouper.has_pending());

        // Now add R2 record
        grouper
            .add_bytes_for_stream(1, b"@read1/2\nTGCA\n+\nJJJJ\n")
            .expect("add_bytes_for_stream failed");

        let templates =
            grouper.drain_complete_templates().expect("drain_complete_templates failed");
        assert_eq!(templates.len(), 1);
    }

    #[test]
    fn test_fastq_grouper_finish_empty() {
        let mut grouper = FastqGrouper::new(2);
        let result = grouper.finish().expect("finish should succeed");
        assert!(result.is_none());
    }

    #[test]
    fn test_fastq_grouper_out_of_sync() {
        let mut grouper = FastqGrouper::new(2);

        // Add mismatched records
        grouper
            .add_bytes_for_stream(0, b"@read1/1\nACGT\n+\nIIII\n")
            .expect("add_bytes_for_stream failed");
        grouper
            .add_bytes_for_stream(1, b"@read2/2\nTGCA\n+\nJJJJ\n")
            .expect("add_bytes_for_stream failed");

        let result = grouper.drain_complete_templates();
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("out of sync"));
    }

    // ========================================================================
    // RecordPositionGrouper tests
    // ========================================================================

    use fgumi_raw_bam::{SamBuilder as RawSamBuilder, flags as raw_flags};

    /// Helper: create a `DecodedRecord` with the given `GroupKey` and flags/tags.
    fn make_decoded(
        key: GroupKey,
        paired: bool,
        first_segment: bool,
        mc: Option<&str>,
    ) -> DecodedRecord {
        use fgumi_raw_bam::testutil::encode_op;
        let mut flags: u16 = 0;
        if paired {
            flags |= raw_flags::PAIRED;
        }
        if first_segment {
            flags |= raw_flags::FIRST_SEGMENT;
        }
        let mut b = RawSamBuilder::new();
        b.read_name(b"read1")
            .sequence(b"ACGT")
            .qualities(&[30; 4])
            .flags(flags)
            .ref_id(0)
            .pos(99) // alignment_start=100 in 1-based => 0-based pos=99
            .cigar_ops(&[encode_op(0, 4)]);
        if let Some(mc_val) = mc {
            b.add_string_tag(b"MC", mc_val.as_bytes());
        }
        DecodedRecord::from_raw_bytes(b.build(), key)
    }

    /// Helper: create a secondary/supplementary `DecodedRecord` with UNKNOWN key.
    fn make_secondary_decoded(name_hash: u64) -> DecodedRecord {
        let key = GroupKey { name_hash, ..GroupKey::default() };
        let mut b = RawSamBuilder::new();
        b.read_name(b"read1").sequence(b"ACGT").qualities(&[30; 4]).flags(raw_flags::SECONDARY);
        DecodedRecord::from_raw_bytes(b.build(), key)
    }

    #[test]
    fn test_record_position_grouper_empty() {
        let mut grouper = RecordPositionGrouper::new();
        assert!(!grouper.has_pending());
        let result = grouper.finish().expect("finish should succeed");
        assert!(result.is_none());
    }

    #[test]
    fn test_record_position_grouper_single_unpaired_record() {
        let mut grouper = RecordPositionGrouper::new();
        let key = GroupKey::single(0, 100, 0, 0, 0, 12345);
        let decoded = make_decoded(key, false, false, None);

        let groups = grouper.add_records(vec![decoded]).expect("add_records should succeed");
        assert!(groups.is_empty()); // Not emitted yet
        assert!(grouper.has_pending());

        let final_group =
            grouper.finish().expect("finish should succeed").expect("should emit final group");
        assert_eq!(final_group.records.len(), 1);
        assert_eq!(final_group.group_key.ref_id1, 0);
        assert_eq!(final_group.group_key.pos1, 100);
    }

    #[test]
    fn test_record_position_grouper_same_position_multiple_records() {
        let mut grouper = RecordPositionGrouper::new();
        let key1 = GroupKey::paired(0, 100, 0, 0, 200, 1, 0, 0, 11111);
        let key2 = GroupKey::paired(0, 100, 0, 0, 200, 1, 0, 0, 22222);
        let key3 = GroupKey::paired(0, 100, 0, 0, 200, 1, 0, 0, 33333);

        let records = vec![
            make_decoded(key1, true, true, Some("4M")),
            make_decoded(key2, true, true, Some("4M")),
            make_decoded(key3, true, true, Some("4M")),
        ];

        let groups = grouper.add_records(records).expect("add_records should succeed");
        assert!(groups.is_empty()); // All same position — not emitted yet

        let final_group =
            grouper.finish().expect("finish should succeed").expect("should emit final group");
        assert_eq!(final_group.records.len(), 3);
    }

    #[test]
    fn test_record_position_grouper_different_positions() {
        let mut grouper = RecordPositionGrouper::new();
        let key_pos1 = GroupKey::single(0, 100, 0, 0, 0, 11111);
        let key_pos2 = GroupKey::single(0, 200, 0, 0, 0, 22222);
        let key_pos3 = GroupKey::single(0, 300, 0, 0, 0, 33333);

        let records = vec![
            make_decoded(key_pos1, false, false, None),
            make_decoded(key_pos2, false, false, None),
            make_decoded(key_pos3, false, false, None),
        ];

        let groups = grouper.add_records(records).expect("add_records should succeed");
        // First two positions emitted when boundary detected
        assert_eq!(groups.len(), 2);
        assert_eq!(groups[0].group_key.pos1, 100);
        assert_eq!(groups[0].records.len(), 1);
        assert_eq!(groups[1].group_key.pos1, 200);
        assert_eq!(groups[1].records.len(), 1);

        // Third position still pending
        let final_group =
            grouper.finish().expect("finish should succeed").expect("should emit final group");
        assert_eq!(final_group.group_key.pos1, 300);
    }

    #[test]
    fn test_record_position_grouper_skips_secondary() {
        let mut grouper = RecordPositionGrouper::new();
        let primary_key = GroupKey::single(0, 100, 0, 0, 0, 11111);

        let records = vec![
            make_decoded(primary_key, false, false, None),
            make_secondary_decoded(11111), // Should be skipped
        ];

        let groups = grouper.add_records(records).expect("add_records should succeed");
        assert!(groups.is_empty());

        let final_group =
            grouper.finish().expect("finish should succeed").expect("should emit final group");
        assert_eq!(final_group.records.len(), 1); // Only primary kept
    }

    #[test]
    fn test_record_position_grouper_paired_same_position_key() {
        // R1 and R2 with MC tags produce identical position_key after normalization
        let mut grouper = RecordPositionGrouper::new();
        // Both R1 and R2 normalize to the same position key
        let r1_key = GroupKey::paired(0, 100, 0, 0, 200, 1, 0, 0, 12345);
        let r2_key = GroupKey::paired(0, 100, 0, 0, 200, 1, 0, 0, 12345);
        // position_key() should be identical for r1 and r2
        assert_eq!(r1_key.position_key(), r2_key.position_key());

        let records = vec![
            make_decoded(r1_key, true, true, Some("4M")),
            make_decoded(r2_key, true, false, Some("4M")),
        ];

        let groups = grouper.add_records(records).expect("add_records should succeed");
        assert!(groups.is_empty()); // Same position — all in one group

        let final_group =
            grouper.finish().expect("finish should succeed").expect("should emit final group");
        assert_eq!(final_group.records.len(), 2);
    }

    #[test]
    fn test_record_position_grouper_groups_records_by_position() {
        let mut grouper = RecordPositionGrouper::new();
        // Group 1: 3 records at position 100
        let key1 = GroupKey::single(0, 100, 0, 0, 0, 11111);
        let key2 = GroupKey::single(0, 100, 0, 0, 0, 22222);
        let key3 = GroupKey::single(0, 100, 0, 0, 0, 33333);
        // Group 2: 2 records at position 200
        let key4 = GroupKey::single(0, 200, 0, 0, 0, 44444);
        let key5 = GroupKey::single(0, 200, 0, 0, 0, 55555);
        // Group 3: 1 record at position 300
        let key6 = GroupKey::single(0, 300, 0, 0, 0, 66666);

        let records = vec![
            make_decoded(key1, false, false, None),
            make_decoded(key2, false, false, None),
            make_decoded(key3, false, false, None),
            make_decoded(key4, false, false, None),
            make_decoded(key5, false, false, None),
            make_decoded(key6, false, false, None),
        ];

        let groups = grouper.add_records(records).expect("add_records should succeed");
        assert_eq!(groups.len(), 2);
        assert_eq!(groups[0].records.len(), 3);
        assert_eq!(groups[1].records.len(), 2);

        let final_group =
            grouper.finish().expect("finish should succeed").expect("should emit final group");
        assert_eq!(final_group.records.len(), 1);
    }

    #[test]
    fn test_record_position_grouper_coalesces_unmapped_mate_by_name_hash() {
        // In template-coordinate sort, a mapped R1 and its unmapped R2 are adjacent.
        // R1 gets GroupKey::single(5, 100, ...) and R2 gets GroupKey::single(-1, 0, ...).
        // They should stay in the same group because they share name_hash.
        let mut grouper = RecordPositionGrouper::new();
        let name_hash = 12345_u64;

        // R1: mapped at chr5:100, mate unmapped — no MC tag
        let r1_key = GroupKey::single(5, 100, 0, 0, 0, name_hash);
        let r1 = {
            use fgumi_raw_bam::testutil::encode_op;
            let mut b = RawSamBuilder::new();
            b.read_name(b"read1")
                .sequence(b"ACGT")
                .qualities(&[30; 4])
                .flags(raw_flags::PAIRED | raw_flags::FIRST_SEGMENT | raw_flags::MATE_UNMAPPED)
                .ref_id(5)
                .pos(99)
                .cigar_ops(&[encode_op(0, 4)]);
            DecodedRecord::from_raw_bytes(b.build(), r1_key)
        };

        // R2: unmapped, mate mapped — different position key
        let r2_key = GroupKey::single(-1, 0, 0, 0, 0, name_hash);
        let r2 = {
            let mut b = RawSamBuilder::new();
            b.read_name(b"read1")
                .sequence(b"TGCA")
                .qualities(&[30; 4])
                .flags(raw_flags::PAIRED | raw_flags::LAST_SEGMENT | raw_flags::UNMAPPED);
            DecodedRecord::from_raw_bytes(b.build(), r2_key)
        };

        // Verify position keys differ (this is the bug scenario)
        assert_ne!(r1_key.position_key(), r2_key.position_key());

        let groups = grouper.add_records(vec![r1, r2]).expect("add_records should succeed");
        // Both should be in the same group — no emission yet
        assert!(groups.is_empty());

        let final_group =
            grouper.finish().expect("finish should succeed").expect("should emit final group");
        assert_eq!(final_group.records.len(), 2, "R1 and R2 should be in the same group");
        assert_eq!(final_group.group_key.ref_id1, 5, "Group key should use R1's position");
    }

    #[test]
    fn test_record_position_grouper_does_not_coalesce_different_name_hash() {
        // Records with different name_hashes at different positions should NOT coalesce.
        let mut grouper = RecordPositionGrouper::new();

        let r1_key = GroupKey::single(0, 100, 0, 0, 0, 11111);
        let r2_key = GroupKey::single(0, 200, 0, 0, 0, 22222);

        let records = vec![
            make_decoded(r1_key, false, false, None),
            make_decoded(r2_key, false, false, None),
        ];

        let groups = grouper.add_records(records).expect("add_records should succeed");
        // Position boundary with different name_hash → group emitted
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].records.len(), 1);
    }

    #[test]
    fn test_record_position_grouper_mc_validation_skips_unmapped_mate() {
        // Paired records with unmapped mates have no MC tag — validation should skip them.
        use fgumi_raw_bam::testutil::encode_op;
        let mut grouper = RecordPositionGrouper::new();
        let key = GroupKey::single(0, 100, 0, 0, 0, 12345);
        let mut b = RawSamBuilder::new();
        b.read_name(b"read1")
            .sequence(b"ACGT")
            .qualities(&[30; 4])
            .flags(raw_flags::PAIRED | raw_flags::FIRST_SEGMENT | raw_flags::MATE_UNMAPPED)
            .ref_id(0)
            .pos(99)
            .cigar_ops(&[encode_op(0, 4)]);
        let decoded = DecodedRecord::from_raw_bytes(b.build(), key);

        // Should NOT error even though there's no MC tag
        let result = grouper.add_records(vec![decoded]);
        assert!(result.is_ok());
        assert!(!grouper.mc_validated); // Not validated — skipped due to unmapped mate
    }

    #[test]
    fn test_record_position_grouper_mc_validation_fails_without_mc() {
        let mut grouper = RecordPositionGrouper::new();
        // Paired record without MC tag
        let key = GroupKey::paired(0, 100, 0, 0, 200, 1, 0, 0, 12345);
        let decoded = make_decoded(key, true, true, None); // No MC tag

        let result = grouper.add_records(vec![decoded]);
        assert!(result.is_err());
        let err_msg = result.unwrap_err().to_string();
        assert!(err_msg.contains("MC tags"), "Error should mention MC tags: {err_msg}");
    }

    #[test]
    fn test_record_position_grouper_mc_validation_passes_with_mc() {
        let mut grouper = RecordPositionGrouper::new();
        let key = GroupKey::paired(0, 100, 0, 0, 200, 1, 0, 0, 12345);
        let decoded = make_decoded(key, true, true, Some("4M"));

        let result = grouper.add_records(vec![decoded]);
        assert!(result.is_ok());
    }

    #[test]
    fn test_record_position_grouper_mc_validation_catches_later_missing_mc() {
        // Regression: after a record with a valid MC tag, a later paired primary
        // missing MC must still fail validation (no short-circuit).
        let mut grouper = RecordPositionGrouper::new();
        let key_ok = GroupKey::paired(0, 100, 0, 0, 200, 1, 0, 0, 11111);
        let ok_record = make_decoded(key_ok, true, true, Some("4M"));
        grouper.add_records(vec![ok_record]).expect("first record with MC should validate");
        assert!(grouper.mc_validated);

        let key_bad = GroupKey::paired(0, 300, 0, 0, 400, 1, 0, 0, 22222);
        let bad_record = make_decoded(key_bad, true, true, None);
        let result = grouper.add_records(vec![bad_record]);
        assert!(result.is_err(), "later paired primary missing MC must fail");
        let err_msg = result.unwrap_err().to_string();
        assert!(err_msg.contains("MC tags"), "Error should mention MC tags: {err_msg}");
    }

    #[test]
    fn test_record_position_grouper_mc_validation_skips_unpaired() {
        // Unpaired records should not trigger MC validation
        let mut grouper = RecordPositionGrouper::new();
        let key = GroupKey::single(0, 100, 0, 0, 0, 12345);
        let decoded = make_decoded(key, false, false, None); // Unpaired, no MC tag

        let result = grouper.add_records(vec![decoded]);
        assert!(result.is_ok());
        assert!(!grouper.mc_validated); // Should not be validated for unpaired
    }

    #[test]
    fn test_record_position_grouper_default_impl() {
        let grouper = RecordPositionGrouper::default();
        assert!(!grouper.has_pending());
    }

    // ========================================================================
    // build_templates_from_records tests
    // ========================================================================

    #[test]
    fn test_build_templates_empty() {
        let result = build_templates_from_records(vec![]).expect("build templates from records");
        assert!(result.is_empty());
    }

    #[test]
    fn test_build_templates_single_record() {
        use fgumi_raw_bam::testutil::encode_op;
        let key = GroupKey::single(0, 100, 0, 0, 0, 12345);
        let mut b = RawSamBuilder::new();
        b.read_name(b"read1")
            .sequence(b"ACGT")
            .qualities(&[30; 4])
            .flags(raw_flags::PAIRED | raw_flags::FIRST_SEGMENT)
            .ref_id(0)
            .pos(99)
            .mapq(60)
            .cigar_ops(&[encode_op(0, 4)]);
        let decoded = DecodedRecord::from_raw_bytes(b.build(), key);

        let templates =
            build_templates_from_records(vec![decoded]).expect("build templates from records");
        assert_eq!(templates.len(), 1);
    }

    #[test]
    fn test_build_templates_paired_same_name_hash() {
        // R1 and R2 with same name_hash should produce one template
        use fgumi_raw_bam::testutil::encode_op;
        let r1_key = GroupKey::paired(0, 100, 0, 0, 200, 1, 0, 0, 12345);
        let r2_key = GroupKey::paired(0, 100, 0, 0, 200, 1, 0, 0, 12345);

        let r1 = {
            let mut b = RawSamBuilder::new();
            b.read_name(b"read1")
                .sequence(b"ACGT")
                .qualities(&[30; 4])
                .flags(raw_flags::PAIRED | raw_flags::FIRST_SEGMENT)
                .ref_id(0)
                .pos(99)
                .cigar_ops(&[encode_op(0, 4)]);
            DecodedRecord::from_raw_bytes(b.build(), r1_key)
        };
        let r2 = {
            let mut b = RawSamBuilder::new();
            b.read_name(b"read1")
                .sequence(b"TGCA")
                .qualities(&[30; 4])
                .flags(raw_flags::PAIRED | raw_flags::LAST_SEGMENT | raw_flags::REVERSE)
                .ref_id(0)
                .pos(199)
                .cigar_ops(&[encode_op(0, 4)]);
            DecodedRecord::from_raw_bytes(b.build(), r2_key)
        };

        let templates =
            build_templates_from_records(vec![r1, r2]).expect("build templates from records");
        assert_eq!(templates.len(), 1);
    }

    #[test]
    fn test_build_templates_multiple_qnames() {
        // Records with different name_hashes should produce separate templates
        use fgumi_raw_bam::testutil::encode_op;
        let key1 = GroupKey::single(0, 100, 0, 0, 0, 11111);
        let key2 = GroupKey::single(0, 100, 0, 0, 0, 22222);
        let key3 = GroupKey::single(0, 100, 0, 0, 0, 33333);

        let make_rec = |name: &[u8]| {
            let mut b = RawSamBuilder::new();
            b.read_name(name)
                .sequence(b"ACGT")
                .qualities(&[30; 4])
                .flags(0)
                .ref_id(0)
                .pos(99)
                .cigar_ops(&[encode_op(0, 4)]);
            b.build()
        };

        let records: Vec<DecodedRecord> = vec![
            DecodedRecord::from_raw_bytes(make_rec(b"readA"), key1),
            DecodedRecord::from_raw_bytes(make_rec(b"readB"), key2),
            DecodedRecord::from_raw_bytes(make_rec(b"readC"), key3),
        ];

        let templates =
            build_templates_from_records(records).expect("build templates from records");
        assert_eq!(templates.len(), 3);
    }

    #[test]
    fn test_build_templates_from_raw_bytes() {
        use crate::sort::bam_fields;
        use crate::unified_pipeline::{DecodedRecord, GroupKey};

        let key = GroupKey::single(0, 100, 0, 0, 0, 12345);
        let raw = bam_fields::make_bam_bytes(
            0,                                                            // tid
            100,                                                          // pos
            bam_fields::flags::PAIRED | bam_fields::flags::FIRST_SEGMENT, // flags
            b"read1",                                                     // name
            &[bam_fields::encode_op(0, 4)],                               // 4M cigar
            4,                                                            // seq_len
            0,                                                            // mate_tid
            200,                                                          // mate_pos
            &[],                                                          // aux
        );
        let decoded = DecodedRecord::from_raw_bytes(raw, key);

        let templates =
            build_templates_from_records(vec![decoded]).expect("build templates from records");
        assert_eq!(templates.len(), 1);
    }

    #[test]
    fn test_build_templates_from_raw_bytes_paired() {
        use crate::sort::bam_fields;
        use crate::unified_pipeline::{DecodedRecord, GroupKey};

        let key1 = GroupKey::paired(0, 100, 0, 0, 200, 1, 0, 0, 12345);
        let key2 = GroupKey::paired(0, 100, 0, 0, 200, 1, 0, 0, 12345);

        let r1 = bam_fields::make_bam_bytes(
            0,
            100,
            bam_fields::flags::PAIRED | bam_fields::flags::FIRST_SEGMENT,
            b"read1",
            &[bam_fields::encode_op(0, 4)],
            4,
            0,
            200,
            &[],
        );
        let r2 = bam_fields::make_bam_bytes(
            0,
            200,
            bam_fields::flags::PAIRED
                | bam_fields::flags::LAST_SEGMENT
                | bam_fields::flags::REVERSE,
            b"read1",
            &[bam_fields::encode_op(0, 4)],
            4,
            0,
            100,
            &[],
        );

        let records =
            vec![DecodedRecord::from_raw_bytes(r1, key1), DecodedRecord::from_raw_bytes(r2, key2)];

        let templates =
            build_templates_from_records(records).expect("build templates from records");
        assert_eq!(templates.len(), 1);
    }

    // ========================================================================
    // RawPositionGroup trait impl tests
    // ========================================================================

    #[test]
    fn test_raw_position_group_batch_weight() {
        let key = GroupKey::single(0, 100, 0, 0, 0, 12345);
        let records = vec![
            make_decoded(GroupKey::single(0, 100, 0, 0, 0, 11111), false, false, None),
            make_decoded(GroupKey::single(0, 100, 0, 0, 0, 22222), false, false, None),
            make_decoded(GroupKey::single(0, 100, 0, 0, 0, 33333), false, false, None),
            make_decoded(GroupKey::single(0, 100, 0, 0, 0, 44444), false, false, None),
        ];
        let group = RawPositionGroup { group_key: key, records };

        // div_ceil(4, 2) = 2 (paired-end estimate)
        assert_eq!(group.batch_weight(), 2);
    }

    #[test]
    fn test_raw_position_group_batch_weight_single() {
        let key = GroupKey::single(0, 100, 0, 0, 0, 12345);
        let records =
            vec![make_decoded(GroupKey::single(0, 100, 0, 0, 0, 11111), false, false, None)];
        let group = RawPositionGroup { group_key: key, records };

        // div_ceil(1, 2) = 1
        assert_eq!(group.batch_weight(), 1);
    }

    #[test]
    fn test_raw_position_group_memory_estimate() {
        let key = GroupKey::single(0, 100, 0, 0, 0, 12345);
        let records =
            vec![make_decoded(GroupKey::single(0, 100, 0, 0, 0, 11111), false, false, None)];
        let group = RawPositionGroup { group_key: key, records };

        // Should return a non-zero value
        assert!(group.estimate_heap_size() > 0);
    }
}

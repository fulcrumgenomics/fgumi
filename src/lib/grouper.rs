//! Grouper implementations for the 9-step pipeline.
//!
//! Each grouper is responsible for grouping decoded BAM records according
//! to command-specific rules and emitting complete groups for processing.
//!
//! Note: BAM record parsing is now handled by the Decode step in the pipeline.
//! Groupers receive pre-decoded `RecordBuf` vectors.

use std::io;

use noodles::sam::alignment::RecordBuf;

use crate::sort::bam_fields;
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
        records
            .into_iter()
            .map(|d| {
                d.into_record().ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        "SingleRecordGrouper requires parsed records, got raw bytes",
                    )
                })
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
    /// Records for current template (non-raw mode).
    current_records: Vec<RecordBuf>,
    /// Raw byte records for current template (raw-byte mode).
    current_raw_records: Vec<Vec<u8>>,
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
            current_raw_records: Vec::new(),
            pending_templates: VecDeque::new(),
        }
    }

    /// Flush current template to pending queue if non-empty.
    fn flush_current_template(&mut self) -> io::Result<()> {
        debug_assert!(
            self.current_raw_records.is_empty() || self.current_records.is_empty(),
            "mixed raw/parsed records in same template group"
        );
        if !self.current_raw_records.is_empty() {
            let template =
                Template::from_raw_records(std::mem::take(&mut self.current_raw_records))
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            self.pending_templates.push_back(template);
            self.current_name = None;
            self.current_name_hash = None;
        } else if !self.current_records.is_empty() {
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
        use crate::unified_pipeline::DecodedRecordData;

        // Group records by QNAME (using pre-computed name_hash for fast comparison)
        for decoded in records {
            let name_hash = decoded.key.name_hash;

            match decoded.data {
                DecodedRecordData::Raw(raw) => {
                    // Raw-byte mode: only extract name when starting a new template
                    match (self.current_name_hash, name_hash) {
                        (Some(current_hash), new_hash) if current_hash == new_hash => {
                            self.current_raw_records.push(raw);
                        }
                        _ => {
                            self.flush_current_template()?;
                            self.current_name = Some(bam_fields::read_name(&raw).to_vec());
                            self.current_name_hash = Some(name_hash);
                            self.current_raw_records.push(raw);
                        }
                    }
                }
                DecodedRecordData::Parsed(record) => {
                    let name = record.name().map(|n| Vec::from(<_ as AsRef<[u8]>>::as_ref(n)));
                    match (self.current_name_hash, name_hash) {
                        (Some(current_hash), new_hash) if current_hash == new_hash => {
                            self.current_records.push(record);
                        }
                        _ => {
                            self.flush_current_template()?;
                            self.current_name = name;
                            self.current_name_hash = Some(name_hash);
                            self.current_records.push(record);
                        }
                    }
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
        self.current_name.is_some()
            || !self.pending_templates.is_empty()
            || !self.current_raw_records.is_empty()
            || !self.current_records.is_empty()
    }
}

use noodles::sam::alignment::record::data::field::Tag;

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
    /// since unmapped reads have no CIGAR to report in an MC tag.
    fn validate_mc_tag(&mut self, decoded: &DecodedRecord) -> io::Result<()> {
        use crate::sort::bam_fields;

        if self.mc_validated {
            return Ok(());
        }

        use crate::unified_pipeline::DecodedRecordData;

        match &decoded.data {
            DecodedRecordData::Raw(raw) => {
                let flg = bam_fields::flags(raw);
                let is_paired = (flg & bam_fields::flags::PAIRED) != 0;
                let is_secondary = (flg & bam_fields::flags::SECONDARY) != 0;
                let is_supplementary = (flg & bam_fields::flags::SUPPLEMENTARY) != 0;
                let is_unmapped = (flg & bam_fields::flags::UNMAPPED) != 0;
                let is_mate_unmapped = (flg & bam_fields::flags::MATE_UNMAPPED) != 0;

                if is_paired
                    && !is_secondary
                    && !is_supplementary
                    && !is_unmapped
                    && !is_mate_unmapped
                {
                    if bam_fields::find_mc_tag_in_record(raw).is_none() {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "RecordPositionGrouper requires MC tags on paired-end reads. \
                             Run `fgumi zipper` to add MC tags before `fgumi group`.",
                        ));
                    }
                    self.mc_validated = true;
                }
            }
            DecodedRecordData::Parsed(record) => {
                let flags = record.flags();
                if flags.is_segmented()
                    && !flags.is_secondary()
                    && !flags.is_supplementary()
                    && !flags.is_unmapped()
                    && !flags.is_mate_unmapped()
                {
                    let mc_tag = Tag::from([b'M', b'C']);
                    if record.data().get(&mc_tag).is_none() {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "RecordPositionGrouper requires MC tags on paired-end reads. \
                             Run `fgumi zipper` to add MC tags before `fgumi group`.",
                        ));
                    }
                    self.mc_validated = true;
                }
            }
        }
        Ok(())
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
                if self
                    .current_records
                    .last()
                    .is_some_and(|last| last.key.name_hash == decoded.key.name_hash) =>
            {
                // Different position but same template (name_hash match with previous
                // record). This happens for paired reads with unmapped mates in
                // template-coordinate sorted input: R1 is mapped at some position while
                // R2 is unmapped (position -1:0), but they're adjacent by QNAME.
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
/// Groups records by `name_hash` (QNAME) and builds a [`Template`] from each group.
/// Records must be grouped by QNAME (guaranteed by template-coordinate sort order).
///
/// Designed to run in the parallel Process step after [`RecordPositionGrouper`]
/// emits [`RawPositionGroup`]s.
pub fn build_templates_from_records(records: Vec<DecodedRecord>) -> io::Result<Vec<Template>> {
    if records.is_empty() {
        return Ok(Vec::new());
    }

    // Check if this is raw-byte mode based on the first record
    use crate::unified_pipeline::DecodedRecordData;

    let raw_byte_mode = matches!(records[0].data, DecodedRecordData::Raw(_));

    if raw_byte_mode {
        return build_templates_from_raw_records(records);
    }

    let mut templates = Vec::new();
    let mut current_name_hash: Option<u64> = None;
    let mut current_records: Vec<RecordBuf> = Vec::new();

    for decoded in records {
        let name_hash = decoded.key.name_hash;
        let record = match decoded.data {
            DecodedRecordData::Parsed(r) => r,
            DecodedRecordData::Raw(_) => unreachable!("checked above"),
        };

        match current_name_hash {
            Some(h) if h == name_hash => {
                current_records.push(record);
            }
            Some(_) => {
                let template = Template::from_records(std::mem::take(&mut current_records))
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                templates.push(template);
                current_name_hash = Some(name_hash);
                current_records.push(record);
            }
            None => {
                current_name_hash = Some(name_hash);
                current_records.push(record);
            }
        }
    }

    // Flush last template
    if !current_records.is_empty() {
        let template = Template::from_records(current_records)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        templates.push(template);
    }

    Ok(templates)
}

/// Build templates from raw-byte mode decoded records.
fn build_templates_from_raw_records(records: Vec<DecodedRecord>) -> io::Result<Vec<Template>> {
    use crate::unified_pipeline::DecodedRecordData;

    let mut templates = Vec::new();
    let mut current_name_hash: Option<u64> = None;
    let mut current_raw: Vec<Vec<u8>> = Vec::new();

    for decoded in records {
        let name_hash = decoded.key.name_hash;
        let raw = match decoded.data {
            DecodedRecordData::Raw(v) => v,
            DecodedRecordData::Parsed(_) => unreachable!("raw-byte mode checked by caller"),
        };

        match current_name_hash {
            Some(h) if h == name_hash => {
                current_raw.push(raw);
            }
            Some(_) => {
                let template = Template::from_raw_records(std::mem::take(&mut current_raw))
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                templates.push(template);
                current_name_hash = Some(name_hash);
                current_raw.push(raw);
            }
            None => {
                current_name_hash = Some(name_hash);
                current_raw.push(raw);
            }
        }
    }

    if !current_raw.is_empty() {
        let template = Template::from_raw_records(current_raw)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        templates.push(template);
    }

    Ok(templates)
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
    // RecordPositionGrouper tests
    // ========================================================================

    use crate::sam::builder::RecordBuilder;

    /// Helper: create a `DecodedRecord` with the given `GroupKey` and flags/tags.
    fn make_decoded(
        key: GroupKey,
        paired: bool,
        first_segment: bool,
        mc: Option<&str>,
    ) -> DecodedRecord {
        let mut builder = RecordBuilder::new()
            .name("read1")
            .sequence("ACGT")
            .qualities(&[30, 30, 30, 30])
            .cigar("4M")
            .reference_sequence_id(0)
            .alignment_start(100)
            .paired(paired);
        if first_segment {
            builder = builder.first_segment(true);
        }
        if let Some(mc_val) = mc {
            builder = builder.tag("MC", mc_val);
        }
        DecodedRecord::new(builder.build(), key)
    }

    /// Helper: create a secondary/supplementary `DecodedRecord` with UNKNOWN key.
    fn make_secondary_decoded(name_hash: u64) -> DecodedRecord {
        let key = GroupKey { name_hash, ..GroupKey::default() };
        let record = RecordBuilder::new()
            .name("read1")
            .sequence("ACGT")
            .qualities(&[30, 30, 30, 30])
            .secondary(true)
            .build();
        DecodedRecord::new(record, key)
    }

    #[test]
    fn test_record_position_grouper_empty() {
        let mut grouper = RecordPositionGrouper::new();
        assert!(!grouper.has_pending());
        let result = grouper.finish().unwrap();
        assert!(result.is_none());
    }

    #[test]
    fn test_record_position_grouper_single_unpaired_record() {
        let mut grouper = RecordPositionGrouper::new();
        let key = GroupKey::single(0, 100, 0, 0, 0, 12345);
        let decoded = make_decoded(key, false, false, None);

        let groups = grouper.add_records(vec![decoded]).unwrap();
        assert!(groups.is_empty()); // Not emitted yet
        assert!(grouper.has_pending());

        let final_group = grouper.finish().unwrap().expect("should emit final group");
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

        let groups = grouper.add_records(records).unwrap();
        assert!(groups.is_empty()); // All same position — not emitted yet

        let final_group = grouper.finish().unwrap().expect("should emit final group");
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

        let groups = grouper.add_records(records).unwrap();
        // First two positions emitted when boundary detected
        assert_eq!(groups.len(), 2);
        assert_eq!(groups[0].group_key.pos1, 100);
        assert_eq!(groups[0].records.len(), 1);
        assert_eq!(groups[1].group_key.pos1, 200);
        assert_eq!(groups[1].records.len(), 1);

        // Third position still pending
        let final_group = grouper.finish().unwrap().expect("should emit final group");
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

        let groups = grouper.add_records(records).unwrap();
        assert!(groups.is_empty());

        let final_group = grouper.finish().unwrap().expect("should emit final group");
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

        let groups = grouper.add_records(records).unwrap();
        assert!(groups.is_empty()); // Same position — all in one group

        let final_group = grouper.finish().unwrap().expect("should emit final group");
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

        let groups = grouper.add_records(records).unwrap();
        assert_eq!(groups.len(), 2);
        assert_eq!(groups[0].records.len(), 3);
        assert_eq!(groups[1].records.len(), 2);

        let final_group = grouper.finish().unwrap().expect("should emit final group");
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
            let record = RecordBuilder::new()
                .name("read1")
                .sequence("ACGT")
                .qualities(&[30, 30, 30, 30])
                .cigar("4M")
                .reference_sequence_id(5)
                .alignment_start(100)
                .paired(true)
                .first_segment(true)
                .mate_unmapped(true)
                .build();
            DecodedRecord::new(record, r1_key)
        };

        // R2: unmapped, mate mapped — different position key
        let r2_key = GroupKey::single(-1, 0, 0, 0, 0, name_hash);
        let r2 = {
            let record = RecordBuilder::new()
                .name("read1")
                .sequence("TGCA")
                .qualities(&[30, 30, 30, 30])
                .paired(true)
                .first_segment(false)
                .unmapped(true)
                .build();
            DecodedRecord::new(record, r2_key)
        };

        // Verify position keys differ (this is the bug scenario)
        assert_ne!(r1_key.position_key(), r2_key.position_key());

        let groups = grouper.add_records(vec![r1, r2]).unwrap();
        // Both should be in the same group — no emission yet
        assert!(groups.is_empty());

        let final_group = grouper.finish().unwrap().expect("should emit final group");
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

        let groups = grouper.add_records(records).unwrap();
        // Position boundary with different name_hash → group emitted
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].records.len(), 1);
    }

    #[test]
    fn test_record_position_grouper_mc_validation_skips_unmapped_mate() {
        // Paired records with unmapped mates have no MC tag — validation should skip them.
        let mut grouper = RecordPositionGrouper::new();
        let key = GroupKey::single(0, 100, 0, 0, 0, 12345);
        let record = RecordBuilder::new()
            .name("read1")
            .sequence("ACGT")
            .qualities(&[30, 30, 30, 30])
            .cigar("4M")
            .reference_sequence_id(0)
            .alignment_start(100)
            .paired(true)
            .first_segment(true)
            .mate_unmapped(true)
            .build();
        let decoded = DecodedRecord::new(record, key);

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
        let result = build_templates_from_records(vec![]).unwrap();
        assert!(result.is_empty());
    }

    #[test]
    fn test_build_templates_single_record() {
        let key = GroupKey::single(0, 100, 0, 0, 0, 12345);
        let record = RecordBuilder::new()
            .name("read1")
            .sequence("ACGT")
            .qualities(&[30, 30, 30, 30])
            .cigar("4M")
            .reference_sequence_id(0)
            .alignment_start(100)
            .mapping_quality(60)
            .paired(true)
            .first_segment(true)
            .build();
        let decoded = DecodedRecord::new(record, key);

        let templates = build_templates_from_records(vec![decoded]).unwrap();
        assert_eq!(templates.len(), 1);
    }

    #[test]
    fn test_build_templates_paired_same_name_hash() {
        // R1 and R2 with same name_hash should produce one template
        let r1_key = GroupKey::paired(0, 100, 0, 0, 200, 1, 0, 0, 12345);
        let r2_key = GroupKey::paired(0, 100, 0, 0, 200, 1, 0, 0, 12345);

        let r1 = DecodedRecord::new(
            RecordBuilder::new()
                .name("read1")
                .sequence("ACGT")
                .qualities(&[30, 30, 30, 30])
                .cigar("4M")
                .reference_sequence_id(0)
                .alignment_start(100)
                .paired(true)
                .first_segment(true)
                .build(),
            r1_key,
        );
        let r2 = DecodedRecord::new(
            RecordBuilder::new()
                .name("read1")
                .sequence("TGCA")
                .qualities(&[30, 30, 30, 30])
                .cigar("4M")
                .reference_sequence_id(0)
                .alignment_start(200)
                .paired(true)
                .reverse_complement(true)
                .build(),
            r2_key,
        );

        let templates = build_templates_from_records(vec![r1, r2]).unwrap();
        assert_eq!(templates.len(), 1);
    }

    #[test]
    fn test_build_templates_multiple_qnames() {
        // Records with different name_hashes should produce separate templates
        let key1 = GroupKey::single(0, 100, 0, 0, 0, 11111);
        let key2 = GroupKey::single(0, 100, 0, 0, 0, 22222);
        let key3 = GroupKey::single(0, 100, 0, 0, 0, 33333);

        let records: Vec<DecodedRecord> = vec![
            DecodedRecord::new(
                RecordBuilder::new()
                    .name("readA")
                    .sequence("ACGT")
                    .qualities(&[30, 30, 30, 30])
                    .cigar("4M")
                    .reference_sequence_id(0)
                    .alignment_start(100)
                    .build(),
                key1,
            ),
            DecodedRecord::new(
                RecordBuilder::new()
                    .name("readB")
                    .sequence("ACGT")
                    .qualities(&[30, 30, 30, 30])
                    .cigar("4M")
                    .reference_sequence_id(0)
                    .alignment_start(100)
                    .build(),
                key2,
            ),
            DecodedRecord::new(
                RecordBuilder::new()
                    .name("readC")
                    .sequence("ACGT")
                    .qualities(&[30, 30, 30, 30])
                    .cigar("4M")
                    .reference_sequence_id(0)
                    .alignment_start(100)
                    .build(),
                key3,
            ),
        ];

        let templates = build_templates_from_records(records).unwrap();
        assert_eq!(templates.len(), 3);
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

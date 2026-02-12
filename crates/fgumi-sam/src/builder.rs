//! Builder for creating test SAM/BAM records and files.
//!
//! This module provides a fluent API for constructing SAM/BAM records for testing,
//! modeled after fgbio's `SamBuilder`. It supports both building individual records
//! and accumulating records into a collection that can be written to a file.
//!
//! ## Builders
//!
//! - [`SamBuilder`]: Accumulates records and manages headers for writing to files
//! - [`RecordBuilder`]: Creates individual records without header management (standalone)
//! - [`ConsensusTagsBuilder`]: Creates consensus-specific SAM tags
//!
//! ## Examples
//!
//! ```rust
//! use fgumi_sam::builder::{RecordBuilder, ConsensusTagsBuilder};
//!
//! // Create a simple record
//! let record = RecordBuilder::new()
//!     .name("read1")
//!     .sequence("ACGT")
//!     .build();
//!
//! // Create a consensus record with tags
//! let consensus = RecordBuilder::new()
//!     .name("consensus1")
//!     .sequence("ACGTACGT")
//!     .consensus_tags(ConsensusTagsBuilder::new().depth_max(10).depth_min(5))
//!     .build();
//! ```

use crate::to_smallest_signed_int;
use anyhow::Result;
use bstr::BString;
use noodles::core::Position;
use noodles::sam::Header;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record::cigar::Op;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value as BufValue;
use noodles::sam::alignment::record_buf::{QualityScores, RecordBuf, Sequence};
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::{Program, ReadGroup, ReferenceSequence};
use std::num::NonZeroUsize;
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};
use tempfile::{NamedTempFile, TempDir};

// Constants
pub const DEFAULT_READ_LENGTH: usize = 100;
pub const DEFAULT_BASE_QUALITY: u8 = 30;
pub const DEFAULT_MAPQ: u8 = 60;
pub const DEFAULT_SAMPLE_NAME: &str = "Sample";
pub const DEFAULT_READ_GROUP_ID: &str = "A";
pub const DEFAULT_REFERENCE_LENGTH: usize = 200_000_000;

/// Strand orientation for reads.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Plus,
    Minus,
}

impl Strand {
    #[must_use]
    pub fn is_negative(&self) -> bool {
        matches!(self, Strand::Minus)
    }
}

/// Builder for creating test SAM/BAM records and files.
///
/// This builder simplifies test record creation by providing methods to add
/// paired-end and single-end reads with sensible defaults. Records are accumulated
/// internally and can be written to a file or iterated over.
///
/// # Examples
///
/// ```rust
/// use fgumi_sam::builder::{SamBuilder, Strand};
///
/// // Create a builder with defaults
/// let mut builder = SamBuilder::new();
///
/// // Add a mapped pair
/// let pair = builder.add_pair()
///     .contig(0)
///     .start1(100)
///     .start2(200)
///     .build();
///
/// // Add a fragment
/// let frag = builder.add_frag()
///     .name("frag1")
///     .bases("ACGTACGT")
///     .contig(0)
///     .start(150)
///     .build();
///
/// // Write to a temp file
/// let path = builder.to_temp_file().unwrap();
/// ```
#[derive(Debug)]
pub struct SamBuilder {
    /// SAM header
    pub header: Header,
    /// Accumulated records
    records: Vec<RecordBuf>,
    /// Default read length
    read_length: usize,
    /// Default base quality
    base_quality: u8,
    /// Read group ID
    read_group_id: String,
    /// Counter for generating sequential names
    counter: AtomicU64,
}

impl Default for SamBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl SamBuilder {
    /// Creates a new builder with default settings.
    ///
    /// Default settings:
    /// - Read length: 100
    /// - Base quality: 30
    /// - Read group ID: "A"
    /// - Sample name: "Sample"
    /// - Reference sequences: chr1-chr22, chrX, chrY, chrM (200MB each)
    #[must_use]
    pub fn new() -> Self {
        Self::with_defaults(DEFAULT_READ_LENGTH, DEFAULT_BASE_QUALITY)
    }

    /// Creates a new builder with specified read length and base quality.
    ///
    /// # Panics
    ///
    /// Panics if `DEFAULT_REFERENCE_LENGTH` is zero (impossible in practice).
    #[must_use]
    pub fn with_defaults(read_length: usize, base_quality: u8) -> Self {
        let mut header = Header::builder();

        // Add reference sequences (chr1-22, X, Y, M)
        for i in 1..=22 {
            let name = format!("chr{i}");
            let map =
                Map::<ReferenceSequence>::new(NonZeroUsize::new(DEFAULT_REFERENCE_LENGTH).unwrap());
            header = header.add_reference_sequence(BString::from(name), map);
        }
        for chr in ["chrX", "chrY", "chrM"] {
            let map =
                Map::<ReferenceSequence>::new(NonZeroUsize::new(DEFAULT_REFERENCE_LENGTH).unwrap());
            header = header.add_reference_sequence(BString::from(chr), map);
        }

        // Add read group
        let rg_map = Map::<ReadGroup>::default();
        header = header.add_read_group(BString::from(DEFAULT_READ_GROUP_ID), rg_map);

        // Add program record
        let pg_map = Map::<Program>::default();
        header = header.add_program(BString::from("SamBuilder"), pg_map);

        let header = header.build();

        Self {
            header,
            records: Vec::new(),
            read_length,
            base_quality,
            read_group_id: DEFAULT_READ_GROUP_ID.to_string(),
            counter: AtomicU64::new(0),
        }
    }

    /// Creates a builder for unmapped reads (no reference sequences in header).
    #[must_use]
    pub fn new_unmapped() -> Self {
        let mut header = Header::builder();

        // Add read group
        let rg_map = Map::<ReadGroup>::default();
        header = header.add_read_group(BString::from(DEFAULT_READ_GROUP_ID), rg_map);

        // Add program record
        let pg_map = Map::<Program>::default();
        header = header.add_program(BString::from("SamBuilder"), pg_map);

        let header = header.build();

        Self {
            header,
            records: Vec::new(),
            read_length: DEFAULT_READ_LENGTH,
            base_quality: DEFAULT_BASE_QUALITY,
            read_group_id: DEFAULT_READ_GROUP_ID.to_string(),
            counter: AtomicU64::new(0),
        }
    }

    /// Creates a builder with a single reference sequence.
    ///
    /// # Panics
    ///
    /// Panics if `ref_length` is zero.
    #[must_use]
    pub fn with_single_ref(ref_name: &str, ref_length: usize) -> Self {
        let mut header = Header::builder();

        let map = Map::<ReferenceSequence>::new(NonZeroUsize::new(ref_length).unwrap());
        header = header.add_reference_sequence(BString::from(ref_name), map);

        // Add read group
        let rg_map = Map::<ReadGroup>::default();
        header = header.add_read_group(BString::from(DEFAULT_READ_GROUP_ID), rg_map);

        // Add program record
        let pg_map = Map::<Program>::default();
        header = header.add_program(BString::from("SamBuilder"), pg_map);

        let header = header.build();

        Self {
            header,
            records: Vec::new(),
            read_length: DEFAULT_READ_LENGTH,
            base_quality: DEFAULT_BASE_QUALITY,
            read_group_id: DEFAULT_READ_GROUP_ID.to_string(),
            counter: AtomicU64::new(0),
        }
    }

    /// Returns the next sequential name.
    fn next_name(&self) -> String {
        format!("{:04}", self.counter.fetch_add(1, Ordering::SeqCst))
    }

    /// Generates random bases of the configured read length.
    fn random_bases(&self) -> String {
        use std::collections::hash_map::DefaultHasher;
        use std::hash::{Hash, Hasher};

        let bases = [b'A', b'C', b'G', b'T'];
        let seed = self.counter.load(Ordering::SeqCst);

        let mut result = Vec::with_capacity(self.read_length);
        for i in 0..self.read_length {
            let mut hasher = DefaultHasher::new();
            (seed, i).hash(&mut hasher);
            let idx = usize::try_from(hasher.finish() % 4).expect("modulo 4 always fits in usize");
            result.push(bases[idx]);
        }
        String::from_utf8(result).unwrap()
    }

    /// Returns a reference to the accumulated records.
    #[must_use]
    pub fn records(&self) -> &[RecordBuf] {
        &self.records
    }

    /// Returns the number of accumulated records.
    #[must_use]
    pub fn len(&self) -> usize {
        self.records.len()
    }

    /// Returns true if no records have been added.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    /// Clears all accumulated records.
    pub fn clear(&mut self) {
        self.records.clear();
    }

    /// Returns an iterator over the accumulated records.
    pub fn iter(&self) -> impl Iterator<Item = &RecordBuf> {
        self.records.iter()
    }

    /// Pushes a raw record to the collection.
    ///
    /// This is useful when you need to add pre-built records to the collection.
    pub fn push_record(&mut self, record: RecordBuf) {
        self.records.push(record);
    }

    /// Starts building a paired-end read pair.
    ///
    /// Returns a `PairBuilder` that can be configured and then built.
    /// The resulting records are added to this builder and also returned.
    #[must_use]
    pub fn add_pair(&mut self) -> PairBuilder<'_> {
        PairBuilder::new(self)
    }

    /// Starts building a single-end (fragment) read.
    ///
    /// Returns a `FragBuilder` that can be configured and then built.
    /// The resulting record is added to this builder and also returned.
    #[must_use]
    pub fn add_frag(&mut self) -> FragBuilder<'_> {
        FragBuilder::new(self)
    }

    /// Writes accumulated records to a BAM file.
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be created or written.
    pub fn write_bam(&self, path: &Path) -> Result<()> {
        let file = std::fs::File::create(path)?;
        let mut writer = noodles::bam::io::Writer::new(file);
        writer.write_header(&self.header)?;

        for record in &self.records {
            writer.write_alignment_record(&self.header, record)?;
        }

        Ok(())
    }

    /// Writes accumulated records to a SAM file.
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be created or written.
    pub fn write_sam(&self, path: &Path) -> Result<()> {
        let file = std::fs::File::create(path)?;
        let mut writer = noodles::sam::io::Writer::new(file);
        writer.write_header(&self.header)?;

        for record in &self.records {
            writer.write_alignment_record(&self.header, record)?;
        }

        Ok(())
    }

    /// Writes to a temporary BAM file and returns the path.
    ///
    /// # Errors
    ///
    /// Returns an error if the temporary file cannot be created or written.
    pub fn to_temp_file(&self) -> Result<NamedTempFile> {
        let temp = NamedTempFile::new()?;
        self.write_bam(temp.path())?;
        Ok(temp)
    }
}

/// Builder for a paired-end read pair.
pub struct PairBuilder<'a> {
    parent: &'a mut SamBuilder,
    name: Option<String>,
    bases1: Option<String>,
    bases2: Option<String>,
    quals1: Option<Vec<u8>>,
    quals2: Option<Vec<u8>>,
    contig: usize,
    contig2: Option<usize>,
    start1: Option<usize>,
    start2: Option<usize>,
    cigar1: Option<String>,
    cigar2: Option<String>,
    mapq1: u8,
    mapq2: u8,
    strand1: Strand,
    strand2: Strand,
    attrs: Vec<(String, BufValue)>,
}

impl<'a> PairBuilder<'a> {
    fn new(parent: &'a mut SamBuilder) -> Self {
        Self {
            parent,
            name: None,
            bases1: None,
            bases2: None,
            quals1: None,
            quals2: None,
            contig: 0,
            contig2: None,
            start1: None,
            start2: None,
            cigar1: None,
            cigar2: None,
            mapq1: DEFAULT_MAPQ,
            mapq2: DEFAULT_MAPQ,
            strand1: Strand::Plus,
            strand2: Strand::Minus,
            attrs: Vec::new(),
        }
    }

    /// Sets the read name (same for both reads in pair).
    #[must_use]
    pub fn name(mut self, name: &str) -> Self {
        self.name = Some(name.to_string());
        self
    }

    /// Sets the bases for R1.
    #[must_use]
    pub fn bases1(mut self, bases: &str) -> Self {
        self.bases1 = Some(bases.to_string());
        self
    }

    /// Sets the bases for R2.
    #[must_use]
    pub fn bases2(mut self, bases: &str) -> Self {
        self.bases2 = Some(bases.to_string());
        self
    }

    /// Sets the quality scores for R1.
    #[must_use]
    pub fn quals1(mut self, quals: &[u8]) -> Self {
        self.quals1 = Some(quals.to_vec());
        self
    }

    /// Sets the quality scores for R2.
    #[must_use]
    pub fn quals2(mut self, quals: &[u8]) -> Self {
        self.quals2 = Some(quals.to_vec());
        self
    }

    /// Sets the reference sequence index for both reads.
    #[must_use]
    pub fn contig(mut self, contig: usize) -> Self {
        self.contig = contig;
        self
    }

    /// Sets a different reference sequence index for R2.
    #[must_use]
    pub fn contig2(mut self, contig: usize) -> Self {
        self.contig2 = Some(contig);
        self
    }

    /// Sets the alignment start for R1 (1-based). If not set, R1 is unmapped.
    #[must_use]
    pub fn start1(mut self, start: usize) -> Self {
        self.start1 = Some(start);
        self
    }

    /// Sets the alignment start for R2 (1-based). If not set, R2 is unmapped.
    #[must_use]
    pub fn start2(mut self, start: usize) -> Self {
        self.start2 = Some(start);
        self
    }

    /// Sets the CIGAR string for R1.
    #[must_use]
    pub fn cigar1(mut self, cigar: &str) -> Self {
        self.cigar1 = Some(cigar.to_string());
        self
    }

    /// Sets the CIGAR string for R2.
    #[must_use]
    pub fn cigar2(mut self, cigar: &str) -> Self {
        self.cigar2 = Some(cigar.to_string());
        self
    }

    /// Sets the mapping quality for R1.
    #[must_use]
    pub fn mapq1(mut self, mapq: u8) -> Self {
        self.mapq1 = mapq;
        self
    }

    /// Sets the mapping quality for R2.
    #[must_use]
    pub fn mapq2(mut self, mapq: u8) -> Self {
        self.mapq2 = mapq;
        self
    }

    /// Sets the strand for R1.
    #[must_use]
    pub fn strand1(mut self, strand: Strand) -> Self {
        self.strand1 = strand;
        self
    }

    /// Sets the strand for R2.
    #[must_use]
    pub fn strand2(mut self, strand: Strand) -> Self {
        self.strand2 = strand;
        self
    }

    /// Marks R1 as unmapped.
    #[must_use]
    pub fn unmapped1(mut self) -> Self {
        self.start1 = None;
        self
    }

    /// Marks R2 as unmapped.
    #[must_use]
    pub fn unmapped2(mut self) -> Self {
        self.start2 = None;
        self
    }

    /// Adds a tag to both reads.
    #[must_use]
    pub fn attr<V: Into<BufValue>>(mut self, tag: &str, value: V) -> Self {
        self.attrs.push((tag.to_string(), value.into()));
        self
    }

    /// Builds the pair, adds to parent, and returns references to both records.
    ///
    /// # Panics
    ///
    /// Panics if the alignment start positions are invalid.
    #[must_use]
    #[expect(clippy::too_many_lines, reason = "test builder with many configuration steps")]
    pub fn build(self) -> (RecordBuf, RecordBuf) {
        let name = self.name.unwrap_or_else(|| self.parent.next_name());
        let bases1 = self.bases1.unwrap_or_else(|| self.parent.random_bases());
        let bases2 = self.bases2.unwrap_or_else(|| self.parent.random_bases());
        let quals1 = self.quals1.unwrap_or_else(|| vec![self.parent.base_quality; bases1.len()]);
        let quals2 = self.quals2.unwrap_or_else(|| vec![self.parent.base_quality; bases2.len()]);
        let cigar1 = self.cigar1.unwrap_or_else(|| format!("{}M", bases1.len()));
        let cigar2 = self.cigar2.unwrap_or_else(|| format!("{}M", bases2.len()));
        let contig2 = self.contig2.unwrap_or(self.contig);

        let unmapped1 = self.start1.is_none();
        let unmapped2 = self.start2.is_none();

        // Build R1
        let mut first_read = RecordBuf::default();
        *first_read.name_mut() = Some(BString::from(name.as_bytes()));
        *first_read.sequence_mut() = Sequence::from(bases1.as_bytes().to_vec());
        *first_read.quality_scores_mut() = QualityScores::from(quals1);

        let mut flags1 = Flags::SEGMENTED | Flags::FIRST_SEGMENT;
        if unmapped1 {
            flags1 |= Flags::UNMAPPED;
        }
        if unmapped2 {
            flags1 |= Flags::MATE_UNMAPPED;
        }
        if self.strand1.is_negative() {
            flags1 |= Flags::REVERSE_COMPLEMENTED;
        }
        if self.strand2.is_negative() {
            flags1 |= Flags::MATE_REVERSE_COMPLEMENTED;
        }
        *first_read.flags_mut() = flags1;

        if !unmapped1 {
            *first_read.reference_sequence_id_mut() = Some(self.contig);
            *first_read.alignment_start_mut() =
                Some(Position::try_from(self.start1.unwrap()).unwrap());
            *first_read.cigar_mut() = parse_cigar(&cigar1).into_iter().collect();
            *first_read.mapping_quality_mut() = Some(
                noodles::sam::alignment::record::MappingQuality::try_from(self.mapq1).unwrap(),
            );
        }

        if !unmapped2 {
            *first_read.mate_reference_sequence_id_mut() = Some(contig2);
            *first_read.mate_alignment_start_mut() =
                Some(Position::try_from(self.start2.unwrap()).unwrap());
        }

        // Add read group tag
        let rg_tag = Tag::new(b'R', b'G');
        first_read
            .data_mut()
            .insert(rg_tag, BufValue::from(self.parent.read_group_id.clone()));

        // Add custom attributes
        for (tag_str, value) in &self.attrs {
            if tag_str.len() == 2 {
                let tag = Tag::new(tag_str.as_bytes()[0], tag_str.as_bytes()[1]);
                first_read.data_mut().insert(tag, value.clone());
            }
        }

        // Build R2
        let mut second_read = RecordBuf::default();
        *second_read.name_mut() = Some(BString::from(name.as_bytes()));
        *second_read.sequence_mut() = Sequence::from(bases2.as_bytes().to_vec());
        *second_read.quality_scores_mut() = QualityScores::from(quals2);

        let mut flags2 = Flags::SEGMENTED | Flags::LAST_SEGMENT;
        if unmapped2 {
            flags2 |= Flags::UNMAPPED;
        }
        if unmapped1 {
            flags2 |= Flags::MATE_UNMAPPED;
        }
        if self.strand2.is_negative() {
            flags2 |= Flags::REVERSE_COMPLEMENTED;
        }
        if self.strand1.is_negative() {
            flags2 |= Flags::MATE_REVERSE_COMPLEMENTED;
        }
        *second_read.flags_mut() = flags2;

        if !unmapped2 {
            *second_read.reference_sequence_id_mut() = Some(contig2);
            *second_read.alignment_start_mut() =
                Some(Position::try_from(self.start2.unwrap()).unwrap());
            *second_read.cigar_mut() = parse_cigar(&cigar2).into_iter().collect();
            *second_read.mapping_quality_mut() = Some(
                noodles::sam::alignment::record::MappingQuality::try_from(self.mapq2).unwrap(),
            );
        }

        if !unmapped1 {
            *second_read.mate_reference_sequence_id_mut() = Some(self.contig);
            *second_read.mate_alignment_start_mut() =
                Some(Position::try_from(self.start1.unwrap()).unwrap());
        }

        // Add read group tag
        second_read
            .data_mut()
            .insert(rg_tag, BufValue::from(self.parent.read_group_id.clone()));

        // Add custom attributes
        for (tag_str, value) in &self.attrs {
            if tag_str.len() == 2 {
                let tag = Tag::new(tag_str.as_bytes()[0], tag_str.as_bytes()[1]);
                second_read.data_mut().insert(tag, value.clone());
            }
        }

        // Calculate template length if both mapped to same contig
        if !unmapped1 && !unmapped2 && self.contig == contig2 {
            let pos1 = i32::try_from(self.start1.unwrap()).expect("start1 fits in i32");
            let pos2 = i32::try_from(self.start2.unwrap()).expect("start2 fits in i32");
            let end1 =
                pos1 + i32::try_from(bases1.len()).expect("bases1 length fits in i32") - 1;
            let end2 =
                pos2 + i32::try_from(bases2.len()).expect("bases2 length fits in i32") - 1;

            let (left, right) = if pos1 <= pos2 { (pos1, end2) } else { (pos2, end1) };
            let tlen = right - left + 1;

            *first_read.template_length_mut() = if pos1 <= pos2 { tlen } else { -tlen };
            *second_read.template_length_mut() = if pos2 <= pos1 { tlen } else { -tlen };
        }

        self.parent.records.push(first_read.clone());
        self.parent.records.push(second_read.clone());

        (first_read, second_read)
    }
}

/// Builder for a single-end (fragment) read.
pub struct FragBuilder<'a> {
    parent: &'a mut SamBuilder,
    name: Option<String>,
    bases: Option<String>,
    quals: Option<Vec<u8>>,
    contig: usize,
    start: Option<usize>,
    cigar: Option<String>,
    mapq: u8,
    strand: Strand,
    attrs: Vec<(String, BufValue)>,
}

impl<'a> FragBuilder<'a> {
    fn new(parent: &'a mut SamBuilder) -> Self {
        Self {
            parent,
            name: None,
            bases: None,
            quals: None,
            contig: 0,
            start: None,
            cigar: None,
            mapq: DEFAULT_MAPQ,
            strand: Strand::Plus,
            attrs: Vec::new(),
        }
    }

    /// Sets the read name.
    #[must_use]
    pub fn name(mut self, name: &str) -> Self {
        self.name = Some(name.to_string());
        self
    }

    /// Sets the bases.
    #[must_use]
    pub fn bases(mut self, bases: &str) -> Self {
        self.bases = Some(bases.to_string());
        self
    }

    /// Sets the quality scores.
    #[must_use]
    pub fn quals(mut self, quals: &[u8]) -> Self {
        self.quals = Some(quals.to_vec());
        self
    }

    /// Sets the reference sequence index.
    #[must_use]
    pub fn contig(mut self, contig: usize) -> Self {
        self.contig = contig;
        self
    }

    /// Sets the alignment start (1-based). If not set, read is unmapped.
    #[must_use]
    pub fn start(mut self, start: usize) -> Self {
        self.start = Some(start);
        self
    }

    /// Sets the CIGAR string.
    #[must_use]
    pub fn cigar(mut self, cigar: &str) -> Self {
        self.cigar = Some(cigar.to_string());
        self
    }

    /// Sets the mapping quality.
    #[must_use]
    pub fn mapq(mut self, mapq: u8) -> Self {
        self.mapq = mapq;
        self
    }

    /// Sets the strand.
    #[must_use]
    pub fn strand(mut self, strand: Strand) -> Self {
        self.strand = strand;
        self
    }

    /// Marks the read as unmapped.
    #[must_use]
    pub fn unmapped(mut self) -> Self {
        self.start = None;
        self
    }

    /// Adds a tag to the read.
    #[must_use]
    pub fn attr<V: Into<BufValue>>(mut self, tag: &str, value: V) -> Self {
        self.attrs.push((tag.to_string(), value.into()));
        self
    }

    /// Builds the record, adds to parent, and returns a clone of the record.
    ///
    /// # Panics
    ///
    /// Panics if the alignment start position is invalid.
    #[must_use]
    pub fn build(self) -> RecordBuf {
        let name = self.name.unwrap_or_else(|| self.parent.next_name());
        let bases = self.bases.unwrap_or_else(|| self.parent.random_bases());
        let quals = self.quals.unwrap_or_else(|| vec![self.parent.base_quality; bases.len()]);
        let cigar = self.cigar.unwrap_or_else(|| format!("{}M", bases.len()));

        let unmapped = self.start.is_none();

        let mut rec = RecordBuf::default();
        *rec.name_mut() = Some(BString::from(name.as_bytes()));
        *rec.sequence_mut() = Sequence::from(bases.as_bytes().to_vec());
        *rec.quality_scores_mut() = QualityScores::from(quals);

        let mut flags = Flags::empty();
        if unmapped {
            flags |= Flags::UNMAPPED;
        }
        if self.strand.is_negative() {
            flags |= Flags::REVERSE_COMPLEMENTED;
        }
        *rec.flags_mut() = flags;

        if !unmapped {
            *rec.reference_sequence_id_mut() = Some(self.contig);
            *rec.alignment_start_mut() = Some(Position::try_from(self.start.unwrap()).unwrap());
            *rec.cigar_mut() = parse_cigar(&cigar).into_iter().collect();
            *rec.mapping_quality_mut() =
                Some(noodles::sam::alignment::record::MappingQuality::try_from(self.mapq).unwrap());
        }

        // Add read group tag
        let rg_tag = Tag::new(b'R', b'G');
        rec.data_mut().insert(rg_tag, BufValue::from(self.parent.read_group_id.clone()));

        // Add custom attributes
        for (tag_str, value) in &self.attrs {
            if tag_str.len() == 2 {
                let tag = Tag::new(tag_str.as_bytes()[0], tag_str.as_bytes()[1]);
                rec.data_mut().insert(tag, value.clone());
            }
        }

        self.parent.records.push(rec.clone());
        rec
    }
}

/// Parses a CIGAR string into a vector of operations.
///
/// # Panics
///
/// Panics if the CIGAR string contains invalid characters or formatting.
#[must_use]
pub fn parse_cigar(cigar_str: &str) -> Vec<Op> {
    let mut ops = Vec::new();
    let mut num_str = String::new();

    for c in cigar_str.chars() {
        if c.is_ascii_digit() {
            num_str.push(c);
        } else {
            let len: usize = num_str.parse().expect("Invalid CIGAR: expected number");
            let kind = match c {
                'M' => Kind::Match,
                'I' => Kind::Insertion,
                'D' => Kind::Deletion,
                'N' => Kind::Skip,
                'S' => Kind::SoftClip,
                'H' => Kind::HardClip,
                'P' => Kind::Pad,
                '=' => Kind::SequenceMatch,
                'X' => Kind::SequenceMismatch,
                _ => panic!("Unknown CIGAR operation: {c}"),
            };
            ops.push(Op::new(kind, len));
            num_str.clear();
        }
    }

    ops
}

// ============================================================================
// Convenience Methods
// ============================================================================

/// Default reference length for test data.
pub const REFERENCE_LENGTH: usize = 1_000;

/// Program ID for mapped BAM files.
pub const MAPPED_PG_ID: &str = "mapped";

impl SamBuilder {
    /// Creates a new builder for mapped BAM files with a single reference.
    ///
    /// Creates a builder with "chr1" as the reference name and `REFERENCE_LENGTH` as length.
    #[must_use]
    pub fn new_mapped() -> Self {
        Self::with_single_ref("chr1", REFERENCE_LENGTH)
    }

    /// Adds a pair using positional arguments.
    ///
    /// A convenience method for adding pairs with common parameters.
    ///
    /// # Arguments
    /// * `name` - Read name
    /// * `start1` - Start position for R1 (None = unmapped)
    /// * `start2` - Start position for R2 (None = unmapped)
    /// * `strand1` - true = Plus, false = Minus for R1
    /// * `strand2` - true = Plus, false = Minus for R2
    /// * `attrs` - Attributes to add to both reads
    pub fn add_pair_with_attrs(
        &mut self,
        name: &str,
        start1: Option<usize>,
        start2: Option<usize>,
        strand1: bool,
        strand2: bool,
        attrs: &std::collections::HashMap<&str, BufValue>,
    ) {
        let mut pair = self
            .add_pair()
            .name(name)
            .strand1(if strand1 { Strand::Plus } else { Strand::Minus })
            .strand2(if strand2 { Strand::Plus } else { Strand::Minus });

        if let Some(s) = start1 {
            pair = pair.start1(s);
        }
        if let Some(s) = start2 {
            pair = pair.start2(s);
        }
        for (tag, value) in attrs {
            pair = pair.attr(tag, value.clone());
        }
        let _ = pair.build();
    }

    /// Adds a fragment using positional arguments.
    ///
    /// A convenience method for adding fragments with common parameters.
    ///
    /// # Arguments
    /// * `name` - Read name
    /// * `start` - Start position (None = unmapped)
    /// * `strand` - true = Plus, false = Minus
    /// * `attrs` - Attributes to add
    pub fn add_frag_with_attrs(
        &mut self,
        name: &str,
        start: Option<usize>,
        strand: bool,
        attrs: &std::collections::HashMap<&str, BufValue>,
    ) {
        let mut frag =
            self.add_frag().name(name).strand(if strand { Strand::Plus } else { Strand::Minus });

        if let Some(s) = start {
            frag = frag.start(s);
        }
        for (tag, value) in attrs {
            frag = frag.attr(tag, value.clone());
        }
        let _ = frag.build();
    }

    /// Writes accumulated records to a BAM file.
    ///
    /// Alias for `write_bam`.
    ///
    /// # Errors
    ///
    /// Returns an error if the BAM file cannot be created or written to.
    pub fn write(&self, path: &std::path::Path) -> Result<()> {
        self.write_bam(path)
    }
}

/// Helper to create a reference dictionary file.
///
/// # Panics
///
/// Panics if `ref_length` is zero, since a reference sequence length must be non-zero.
///
/// # Errors
///
/// Returns an error if the dictionary file cannot be created or written to.
pub fn create_ref_dict(
    dir: &TempDir,
    ref_name: &str,
    ref_length: usize,
) -> Result<std::path::PathBuf> {
    let dict_path = dir.path().join("ref.dict");
    let mut header = Header::builder();

    let map = Map::<ReferenceSequence>::new(NonZeroUsize::new(ref_length).unwrap());
    header = header.add_reference_sequence(BString::from(ref_name), map);

    let header = header.build();

    let file = std::fs::File::create(&dict_path)?;
    let mut writer = noodles::sam::io::Writer::new(file);
    writer.write_header(&header)?;

    Ok(dict_path)
}

// ============================================================================
// Standalone Record Builder
// ============================================================================

/// Builder for creating individual BAM/SAM records without a parent `SamBuilder`.
///
/// This builder simplifies test record creation by providing a chainable interface
/// for setting all record fields. All fields have sensible defaults. Unlike
/// `SamBuilder`, this creates records directly without accumulating them.
///
/// # Examples
///
/// ```rust
/// use fgumi_sam::builder::{RecordBuilder, ConsensusTagsBuilder};
/// use noodles::sam::alignment::record::Flags;
///
/// // Create a simple unmapped read
/// let record = RecordBuilder::new()
///     .name("read1")
///     .sequence("ACGT")
///     .qualities(&[30, 30, 30, 30])
///     .build();
///
/// // Create a mapped paired-end R1
/// let r1 = RecordBuilder::new()
///     .name("read1")
///     .sequence("ACGTACGT")
///     .qualities(&[30; 8])
///     .paired(true)
///     .first_segment(true)
///     .reference_sequence_id(0)
///     .alignment_start(100)
///     .cigar("8M")
///     .mapping_quality(60)
///     .build();
///
/// // Create a record with tags
/// let record = RecordBuilder::new()
///     .name("read1")
///     .sequence("ACGT")
///     .tag("RG", "A")
///     .tag("MI", "AAAA-CCCC/A")
///     .build();
/// ```
#[derive(Debug, Default)]
pub struct RecordBuilder {
    name: Option<Vec<u8>>,
    flags: Flags,
    reference_sequence_id: Option<usize>,
    alignment_start: Option<usize>,
    mapping_quality: Option<u8>,
    cigar: Option<String>,
    sequence: Vec<u8>,
    qualities: Vec<u8>,
    tags: Vec<(Tag, BufValue)>,
    mate_reference_sequence_id: Option<usize>,
    mate_alignment_start: Option<usize>,
    template_length: Option<i32>,
}

impl RecordBuilder {
    /// Creates a new builder with default values.
    #[must_use]
    pub fn new() -> Self {
        Self {
            name: None,
            flags: Flags::empty(),
            reference_sequence_id: None,
            alignment_start: None,
            mapping_quality: Some(60), // Default to high quality
            cigar: None,
            sequence: Vec::new(),
            qualities: Vec::new(),
            tags: Vec::new(),
            mate_reference_sequence_id: None,
            mate_alignment_start: None,
            template_length: None,
        }
    }

    /// Creates a new builder pre-configured for a typical mapped read.
    ///
    /// Sets common defaults:
    /// - `reference_sequence_id`: 0
    /// - `mapping_quality`: 60
    /// - CIGAR will be auto-generated as `{len}M` if not explicitly set
    ///
    /// # Examples
    ///
    /// ```rust
    /// use fgumi_sam::builder::RecordBuilder;
    /// use noodles::sam::alignment::record::Cigar;
    ///
    /// let record = RecordBuilder::mapped_read()
    ///     .name("read1")
    ///     .sequence("ACGTACGT")
    ///     .alignment_start(100)
    ///     .build();
    ///
    /// assert_eq!(record.reference_sequence_id(), Some(0));
    /// assert!(!record.cigar().is_empty()); // Auto-generated as "8M"
    /// ```
    #[must_use]
    pub fn mapped_read() -> Self {
        Self { reference_sequence_id: Some(0), ..Self::new() }
    }

    /// Sets the read name.
    #[must_use]
    pub fn name(mut self, name: &str) -> Self {
        self.name = Some(name.as_bytes().to_vec());
        self
    }

    /// Sets the sequence.
    #[must_use]
    pub fn sequence(mut self, seq: &str) -> Self {
        self.sequence = seq.as_bytes().to_vec();
        // Auto-generate qualities if not set
        if self.qualities.is_empty() {
            self.qualities = vec![DEFAULT_BASE_QUALITY; seq.len()];
        }
        self
    }

    /// Sets the quality scores (Phred+33 ASCII).
    #[must_use]
    pub fn qualities(mut self, quals: &[u8]) -> Self {
        self.qualities = quals.to_vec();
        self
    }

    /// Sets all flags at once.
    #[must_use]
    pub fn flags(mut self, flags: Flags) -> Self {
        self.flags = flags;
        self
    }

    /// Sets the paired flag.
    #[must_use]
    pub fn paired(mut self, paired: bool) -> Self {
        self.flags.set(Flags::SEGMENTED, paired);
        self
    }

    /// Sets the first segment (R1) flag. Implies paired.
    #[must_use]
    pub fn first_segment(mut self, is_first: bool) -> Self {
        self.flags.set(Flags::SEGMENTED, true);
        self.flags.set(Flags::FIRST_SEGMENT, is_first);
        if !is_first {
            self.flags.set(Flags::LAST_SEGMENT, true);
        }
        self
    }

    /// Sets the properly paired flag. Implies paired.
    #[must_use]
    pub fn properly_paired(mut self, properly_paired: bool) -> Self {
        if properly_paired {
            self.flags.set(Flags::SEGMENTED, true);
        }
        self.flags.set(Flags::PROPERLY_SEGMENTED, properly_paired);
        self
    }

    /// Sets the unmapped flag.
    #[must_use]
    pub fn unmapped(mut self, unmapped: bool) -> Self {
        self.flags.set(Flags::UNMAPPED, unmapped);
        self
    }

    /// Sets the reverse complement flag.
    #[must_use]
    pub fn reverse_complement(mut self, reverse: bool) -> Self {
        self.flags.set(Flags::REVERSE_COMPLEMENTED, reverse);
        self
    }

    /// Sets the secondary alignment flag.
    #[must_use]
    pub fn secondary(mut self, secondary: bool) -> Self {
        self.flags.set(Flags::SECONDARY, secondary);
        self
    }

    /// Sets the supplementary alignment flag.
    #[must_use]
    pub fn supplementary(mut self, supplementary: bool) -> Self {
        self.flags.set(Flags::SUPPLEMENTARY, supplementary);
        self
    }

    /// Sets the reference sequence ID (0-based).
    #[must_use]
    pub fn reference_sequence_id(mut self, id: usize) -> Self {
        self.reference_sequence_id = Some(id);
        self
    }

    /// Sets the alignment start position (1-based).
    #[must_use]
    pub fn alignment_start(mut self, pos: usize) -> Self {
        self.alignment_start = Some(pos);
        self
    }

    /// Sets the mapping quality.
    #[must_use]
    pub fn mapping_quality(mut self, mapq: u8) -> Self {
        self.mapping_quality = Some(mapq);
        self
    }

    /// Sets the CIGAR string.
    #[must_use]
    pub fn cigar(mut self, cigar: &str) -> Self {
        self.cigar = Some(cigar.to_string());
        self
    }

    /// Sets the mate reference sequence ID (0-based).
    #[must_use]
    pub fn mate_reference_sequence_id(mut self, id: usize) -> Self {
        self.mate_reference_sequence_id = Some(id);
        self
    }

    /// Sets the mate alignment start position (1-based).
    #[must_use]
    pub fn mate_alignment_start(mut self, pos: usize) -> Self {
        self.mate_alignment_start = Some(pos);
        self
    }

    /// Sets the template length (insert size).
    #[must_use]
    pub fn template_length(mut self, tlen: i32) -> Self {
        self.template_length = Some(tlen);
        self
    }

    /// Sets the mate reverse complement flag.
    #[must_use]
    pub fn mate_reverse_complement(mut self, reverse: bool) -> Self {
        self.flags.set(Flags::MATE_REVERSE_COMPLEMENTED, reverse);
        self
    }

    /// Sets the mate unmapped flag.
    #[must_use]
    pub fn mate_unmapped(mut self, unmapped: bool) -> Self {
        self.flags.set(Flags::MATE_UNMAPPED, unmapped);
        self
    }

    /// Adds a SAM tag.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use fgumi_sam::builder::RecordBuilder;
    /// let record = RecordBuilder::new()
    ///     .name("read1")
    ///     .tag("RG", "A")
    ///     .tag("MI", "AAAA-CCCC/A")
    ///     .build();
    /// ```
    #[must_use]
    pub fn tag<V: Into<BufValue>>(mut self, tag: &str, value: V) -> Self {
        let tag_bytes = tag.as_bytes();
        if tag_bytes.len() == 2 {
            let tag = Tag::from([tag_bytes[0], tag_bytes[1]]);
            self.tags.push((tag, value.into()));
        }
        self
    }

    /// Adds consensus tags using a `ConsensusTagsBuilder`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use fgumi_sam::builder::{RecordBuilder, ConsensusTagsBuilder};
    /// let record = RecordBuilder::new()
    ///     .sequence("ACGT")
    ///     .consensus_tags(
    ///         ConsensusTagsBuilder::new()
    ///             .depth_max(10)
    ///             .depth_min(5)
    ///             .error_rate(0.01)
    ///     )
    ///     .build();
    /// ```
    #[must_use]
    pub fn consensus_tags(mut self, builder: ConsensusTagsBuilder) -> Self {
        for (tag, value) in builder.build() {
            self.tags.push((tag, value));
        }
        self
    }

    /// Builds the `RecordBuf`.
    ///
    /// # Panics
    ///
    /// Panics if CIGAR string parsing fails (should only happen with invalid CIGAR).
    #[must_use]
    pub fn build(self) -> RecordBuf {
        let mut record = RecordBuf::default();

        // Set name
        if let Some(name) = self.name {
            *record.name_mut() = Some(name.into());
        }

        // Set flags
        *record.flags_mut() = self.flags;

        // Set reference and position
        if let Some(ref_id) = self.reference_sequence_id {
            *record.reference_sequence_id_mut() = Some(ref_id);
        }
        if let Some(pos) = self.alignment_start {
            *record.alignment_start_mut() =
                Some(Position::try_from(pos).expect("alignment_start must be >= 1"));
        }

        // Set mate reference and position
        if let Some(mate_ref_id) = self.mate_reference_sequence_id {
            *record.mate_reference_sequence_id_mut() = Some(mate_ref_id);
        }
        if let Some(mate_pos) = self.mate_alignment_start {
            *record.mate_alignment_start_mut() =
                Some(Position::try_from(mate_pos).expect("mate_alignment_start must be >= 1"));
        }
        if let Some(tlen) = self.template_length {
            *record.template_length_mut() = tlen;
        }

        // Set mapping quality
        if let Some(mapq) = self.mapping_quality {
            *record.mapping_quality_mut() = Some(
                noodles::sam::alignment::record::MappingQuality::try_from(mapq)
                    .expect("mapping_quality must be valid"),
            );
        }

        // Handle CIGAR and sequence auto-generation:
        // - If both set: use as-is
        // - If only sequence: auto-generate CIGAR as NM
        // - If only CIGAR: auto-generate sequence of appropriate length
        let (cigar_str, sequence) = match (self.cigar, self.sequence.is_empty()) {
            (Some(cigar), true) => {
                // CIGAR provided, no sequence: generate sequence from CIGAR
                let seq_len = cigar_seq_len(&cigar);
                let generated_seq: String =
                    (0..seq_len).map(|i| "ACGT".chars().nth(i % 4).unwrap()).collect();
                (cigar, generated_seq.into_bytes())
            }
            (Some(cigar), false) => (cigar, self.sequence),
            (None, false) => (format!("{}M", self.sequence.len()), self.sequence),
            (None, true) => (String::new(), Vec::new()),
        };

        if !cigar_str.is_empty() {
            let ops = parse_cigar(&cigar_str);
            *record.cigar_mut() = ops.into_iter().collect();
        }

        // Set sequence and qualities (auto-generate qualities if needed)
        let qualities = if self.qualities.is_empty() && !sequence.is_empty() {
            vec![DEFAULT_BASE_QUALITY; sequence.len()]
        } else {
            self.qualities
        };
        *record.sequence_mut() = Sequence::from(sequence);
        *record.quality_scores_mut() = QualityScores::from(qualities);

        // Add tags
        for (tag, value) in self.tags {
            record.data_mut().insert(tag, value);
        }

        record
    }
}

// ============================================================================
// Record Pair Builder (Standalone)
// ============================================================================

/// Builder for creating paired-end read pairs without a parent `SamBuilder`.
///
/// This builder simplifies test pair creation by providing a chainable interface
/// that builds both R1 and R2 with proper mate information and template length.
///
/// # Examples
///
/// ```rust
/// use fgumi_sam::builder::RecordPairBuilder;
///
/// let (r1, r2) = RecordPairBuilder::new()
///     .name("read1")
///     .r1_sequence("ACGTACGT")
///     .r1_start(100)
///     .r2_sequence("TGCATGCA")
///     .r2_start(200)
///     .tag("MI", "1/A")
///     .build();
///
/// assert!(r1.flags().is_first_segment());
/// assert!(r2.flags().is_last_segment());
/// ```
#[derive(Debug, Default)]
#[expect(
    clippy::struct_excessive_bools,
    reason = "test builder struct; bools represent independent SAM flags (r1_reverse, r2_reverse, secondary, supplementary)"
)]
pub struct RecordPairBuilder {
    name: Option<String>,
    r1_sequence: Option<String>,
    r2_sequence: Option<String>,
    r1_qualities: Option<Vec<u8>>,
    r2_qualities: Option<Vec<u8>>,
    r1_start: Option<usize>,
    r2_start: Option<usize>,
    r1_cigar: Option<String>,
    r2_cigar: Option<String>,
    reference_sequence_id: usize,
    r2_reference_sequence_id: Option<usize>,
    mapping_quality: u8,
    r1_reverse: bool,
    r2_reverse: bool,
    secondary: bool,
    supplementary: bool,
    tags: Vec<(String, BufValue)>,
    r1_tags: Vec<(String, BufValue)>,
    r2_tags: Vec<(String, BufValue)>,
}

impl RecordPairBuilder {
    /// Creates a new pair builder with defaults.
    ///
    /// Defaults:
    /// - `reference_sequence_id`: 0 (both reads)
    /// - `mapping_quality`: 60
    /// - R1: forward strand, R2: reverse strand (FR orientation)
    #[must_use]
    pub fn new() -> Self {
        Self { mapping_quality: 60, r2_reverse: true, ..Self::default() }
    }

    /// Sets the read name (shared by both reads).
    #[must_use]
    pub fn name(mut self, name: &str) -> Self {
        self.name = Some(name.to_string());
        self
    }

    /// Sets the R1 sequence.
    #[must_use]
    pub fn r1_sequence(mut self, seq: &str) -> Self {
        self.r1_sequence = Some(seq.to_string());
        self
    }

    /// Sets the R2 sequence.
    #[must_use]
    pub fn r2_sequence(mut self, seq: &str) -> Self {
        self.r2_sequence = Some(seq.to_string());
        self
    }

    /// Sets the R1 quality scores.
    #[must_use]
    pub fn r1_qualities(mut self, quals: &[u8]) -> Self {
        self.r1_qualities = Some(quals.to_vec());
        self
    }

    /// Sets the R2 quality scores.
    #[must_use]
    pub fn r2_qualities(mut self, quals: &[u8]) -> Self {
        self.r2_qualities = Some(quals.to_vec());
        self
    }

    /// Sets the R1 alignment start (1-based).
    #[must_use]
    pub fn r1_start(mut self, start: usize) -> Self {
        self.r1_start = Some(start);
        self
    }

    /// Sets the R2 alignment start (1-based).
    #[must_use]
    pub fn r2_start(mut self, start: usize) -> Self {
        self.r2_start = Some(start);
        self
    }

    /// Sets the R1 CIGAR string.
    #[must_use]
    pub fn r1_cigar(mut self, cigar: &str) -> Self {
        self.r1_cigar = Some(cigar.to_string());
        self
    }

    /// Sets the R2 CIGAR string.
    #[must_use]
    pub fn r2_cigar(mut self, cigar: &str) -> Self {
        self.r2_cigar = Some(cigar.to_string());
        self
    }

    /// Sets the reference sequence ID for both reads.
    #[must_use]
    pub fn reference_sequence_id(mut self, id: usize) -> Self {
        self.reference_sequence_id = id;
        self
    }

    /// Sets a different reference sequence ID for R2.
    #[must_use]
    pub fn r2_reference_sequence_id(mut self, id: usize) -> Self {
        self.r2_reference_sequence_id = Some(id);
        self
    }

    /// Sets the mapping quality for both reads.
    #[must_use]
    pub fn mapping_quality(mut self, mapq: u8) -> Self {
        self.mapping_quality = mapq;
        self
    }

    /// Sets R1 as reverse complemented.
    #[must_use]
    pub fn r1_reverse(mut self, reverse: bool) -> Self {
        self.r1_reverse = reverse;
        self
    }

    /// Sets R2 as reverse complemented.
    #[must_use]
    pub fn r2_reverse(mut self, reverse: bool) -> Self {
        self.r2_reverse = reverse;
        self
    }

    /// Marks both reads as secondary alignments.
    #[must_use]
    pub fn secondary(mut self, secondary: bool) -> Self {
        self.secondary = secondary;
        self
    }

    /// Marks both reads as supplementary alignments.
    #[must_use]
    pub fn supplementary(mut self, supplementary: bool) -> Self {
        self.supplementary = supplementary;
        self
    }

    /// Adds a tag to both reads.
    #[must_use]
    pub fn tag<V: Into<BufValue>>(mut self, tag: &str, value: V) -> Self {
        self.tags.push((tag.to_string(), value.into()));
        self
    }

    /// Adds a tag to R1 only.
    #[must_use]
    pub fn r1_tag<V: Into<BufValue>>(mut self, tag: &str, value: V) -> Self {
        self.r1_tags.push((tag.to_string(), value.into()));
        self
    }

    /// Adds a tag to R2 only.
    #[must_use]
    pub fn r2_tag<V: Into<BufValue>>(mut self, tag: &str, value: V) -> Self {
        self.r2_tags.push((tag.to_string(), value.into()));
        self
    }

    /// Builds the pair, returning (R1, R2).
    ///
    /// # Panics
    ///
    /// Panics if start positions or sequence lengths do not fit in `i32`, which is
    /// only possible with unreasonably large test values.
    #[must_use]
    pub fn build(self) -> (RecordBuf, RecordBuf) {
        let name = self.name.unwrap_or_else(|| "pair".to_string());
        let r1_seq = self.r1_sequence.unwrap_or_else(|| "ACGT".to_string());
        let r2_seq = self.r2_sequence.unwrap_or_else(|| "ACGT".to_string());
        let r1_quals =
            self.r1_qualities.unwrap_or_else(|| vec![DEFAULT_BASE_QUALITY; r1_seq.len()]);
        let r2_quals =
            self.r2_qualities.unwrap_or_else(|| vec![DEFAULT_BASE_QUALITY; r2_seq.len()]);
        let r1_cigar = self.r1_cigar.unwrap_or_else(|| format!("{}M", r1_seq.len()));
        let r2_cigar = self.r2_cigar.unwrap_or_else(|| format!("{}M", r2_seq.len()));
        let r2_ref_id = self.r2_reference_sequence_id.unwrap_or(self.reference_sequence_id);

        let r1_mapped = self.r1_start.is_some();
        let r2_mapped = self.r2_start.is_some();

        // Build R1
        let mut r1_builder = RecordBuilder::new()
            .name(&name)
            .sequence(&r1_seq)
            .qualities(&r1_quals)
            .paired(true)
            .first_segment(true)
            .reverse_complement(self.r1_reverse)
            .mate_reverse_complement(self.r2_reverse)
            .secondary(self.secondary)
            .supplementary(self.supplementary);

        if r1_mapped {
            r1_builder = r1_builder
                .reference_sequence_id(self.reference_sequence_id)
                .alignment_start(self.r1_start.unwrap())
                .cigar(&r1_cigar)
                .mapping_quality(self.mapping_quality);
        } else {
            r1_builder = r1_builder.unmapped(true);
        }

        if r2_mapped {
            r1_builder = r1_builder
                .mate_reference_sequence_id(r2_ref_id)
                .mate_alignment_start(self.r2_start.unwrap());
        } else {
            r1_builder = r1_builder.mate_unmapped(true);
        }

        for (tag, value) in &self.tags {
            r1_builder = r1_builder.tag(tag, value.clone());
        }
        for (tag, value) in &self.r1_tags {
            r1_builder = r1_builder.tag(tag, value.clone());
        }

        // Build R2
        let mut r2_builder = RecordBuilder::new()
            .name(&name)
            .sequence(&r2_seq)
            .qualities(&r2_quals)
            .paired(true)
            .first_segment(false) // Sets LAST_SEGMENT
            .reverse_complement(self.r2_reverse)
            .mate_reverse_complement(self.r1_reverse)
            .secondary(self.secondary)
            .supplementary(self.supplementary);

        if r2_mapped {
            r2_builder = r2_builder
                .reference_sequence_id(r2_ref_id)
                .alignment_start(self.r2_start.unwrap())
                .cigar(&r2_cigar)
                .mapping_quality(self.mapping_quality);
        } else {
            r2_builder = r2_builder.unmapped(true);
        }

        if r1_mapped {
            r2_builder = r2_builder
                .mate_reference_sequence_id(self.reference_sequence_id)
                .mate_alignment_start(self.r1_start.unwrap());
        } else {
            r2_builder = r2_builder.mate_unmapped(true);
        }

        for (tag, value) in &self.tags {
            r2_builder = r2_builder.tag(tag, value.clone());
        }
        for (tag, value) in &self.r2_tags {
            r2_builder = r2_builder.tag(tag, value.clone());
        }

        // Add MC (mate CIGAR) tags when both reads are mapped
        if r1_mapped && r2_mapped {
            r1_builder = r1_builder.tag("MC", r2_cigar.as_str());
            r2_builder = r2_builder.tag("MC", r1_cigar.as_str());
        }

        // Calculate template length if both mapped to same reference
        let mut first_read = r1_builder.build();
        let mut second_read = r2_builder.build();

        if r1_mapped && r2_mapped && self.reference_sequence_id == r2_ref_id {
            let pos1 = i32::try_from(self.r1_start.unwrap()).expect("r1_start fits in i32");
            let pos2 = i32::try_from(self.r2_start.unwrap()).expect("r2_start fits in i32");
            let end1 =
                pos1 + i32::try_from(r1_seq.len()).expect("r1_seq length fits in i32") - 1;
            let end2 =
                pos2 + i32::try_from(r2_seq.len()).expect("r2_seq length fits in i32") - 1;

            let (left, right) = if pos1 <= pos2 { (pos1, end2) } else { (pos2, end1) };
            let tlen = right - left + 1;

            *first_read.template_length_mut() = if pos1 <= pos2 { tlen } else { -tlen };
            *second_read.template_length_mut() = if pos2 <= pos1 { tlen } else { -tlen };
        }

        (first_read, second_read)
    }
}

// ============================================================================
// Consensus Tags Builder
// ============================================================================

/// Builder for consensus-specific SAM tags.
///
/// This builder creates the standard consensus tags used by fgumi consensus callers:
/// - `cD`: Maximum depth (coverage)
/// - `cM`: Minimum depth (coverage)
/// - `cE`: Error rate
/// - Per-base arrays for detailed metrics
///
/// For duplex consensus reads, additional strand-specific tags are supported:
/// - `aD`: A-strand (AB) depth
/// - `bD`: B-strand (BA) depth
/// - `aM`: A-strand minimum depth
/// - `bM`: B-strand minimum depth
/// - `aE`: A-strand error rate
/// - `bE`: B-strand error rate
#[derive(Debug, Default)]
pub struct ConsensusTagsBuilder {
    depth_max: Option<i32>,
    depth_min: Option<i32>,
    error_rate: Option<f32>,
    error_count: Option<i32>,
    per_base_depths: Option<Vec<u16>>,
    per_base_errors: Option<Vec<u16>>,
    // Duplex-specific fields
    ab_depth: Option<i32>,
    ba_depth: Option<i32>,
    ab_min_depth: Option<i32>,
    ba_min_depth: Option<i32>,
    ab_errors: Option<f32>,
    ba_errors: Option<f32>,
    ab_error_count: Option<i32>,
    ba_error_count: Option<i32>,
}

impl ConsensusTagsBuilder {
    /// Creates a new builder.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Sets the maximum depth (cD tag).
    #[must_use]
    pub fn depth_max(mut self, depth: i32) -> Self {
        self.depth_max = Some(depth);
        self
    }

    /// Sets the minimum depth (cM tag).
    #[must_use]
    pub fn depth_min(mut self, depth: i32) -> Self {
        self.depth_min = Some(depth);
        self
    }

    /// Sets the error rate (cE tag).
    #[must_use]
    pub fn error_rate(mut self, rate: f32) -> Self {
        self.error_rate = Some(rate);
        self
    }

    /// Sets per-base depths array.
    #[must_use]
    pub fn per_base_depths(mut self, depths: &[u16]) -> Self {
        self.per_base_depths = Some(depths.to_vec());
        self
    }

    /// Sets per-base errors array.
    #[must_use]
    pub fn per_base_errors(mut self, errors: &[u16]) -> Self {
        self.per_base_errors = Some(errors.to_vec());
        self
    }

    /// Sets the A-strand (AB) depth (aD tag) for duplex consensus.
    #[must_use]
    pub fn ab_depth(mut self, depth: i32) -> Self {
        self.ab_depth = Some(depth);
        self
    }

    /// Sets the B-strand (BA) depth (bD tag) for duplex consensus.
    #[must_use]
    pub fn ba_depth(mut self, depth: i32) -> Self {
        self.ba_depth = Some(depth);
        self
    }

    /// Sets the A-strand minimum depth (aM tag) for duplex consensus.
    #[must_use]
    pub fn ab_min_depth(mut self, depth: i32) -> Self {
        self.ab_min_depth = Some(depth);
        self
    }

    /// Sets the B-strand minimum depth (bM tag) for duplex consensus.
    #[must_use]
    pub fn ba_min_depth(mut self, depth: i32) -> Self {
        self.ba_min_depth = Some(depth);
        self
    }

    /// Sets the A-strand error rate (aE tag) for duplex consensus.
    #[must_use]
    pub fn ab_errors(mut self, rate: f32) -> Self {
        self.ab_errors = Some(rate);
        self
    }

    /// Sets the B-strand error rate (bE tag) for duplex consensus.
    #[must_use]
    pub fn ba_errors(mut self, rate: f32) -> Self {
        self.ba_errors = Some(rate);
        self
    }

    /// Sets the consensus error count (cE tag) as an integer.
    ///
    /// This is an alternative to `error_rate` for tools that store error counts rather than rates.
    #[must_use]
    pub fn tag_ce_as_int(mut self, count: i32) -> Self {
        self.error_count = Some(count);
        self
    }

    /// Sets the A-strand error count (aE tag) as an integer.
    #[must_use]
    pub fn ab_errors_int(mut self, count: i32) -> Self {
        self.ab_error_count = Some(count);
        self
    }

    /// Sets the B-strand error count (bE tag) as an integer.
    #[must_use]
    pub fn ba_errors_int(mut self, count: i32) -> Self {
        self.ba_error_count = Some(count);
        self
    }

    /// Builds the tags as a vector of (Tag, Value) pairs.
    #[must_use]
    pub fn build(self) -> Vec<(Tag, BufValue)> {
        let mut tags = Vec::new();

        // Standard consensus tags
        if let Some(depth) = self.depth_max {
            tags.push((Tag::from([b'c', b'D']), to_smallest_signed_int(depth)));
        }
        if let Some(depth) = self.depth_min {
            tags.push((Tag::from([b'c', b'M']), to_smallest_signed_int(depth)));
        }
        if let Some(rate) = self.error_rate {
            tags.push((Tag::from([b'c', b'E']), BufValue::from(rate)));
        }
        if let Some(count) = self.error_count {
            tags.push((Tag::from([b'c', b'E']), to_smallest_signed_int(count)));
        }
        if let Some(depths) = self.per_base_depths {
            tags.push((Tag::from([b'c', b'D']), BufValue::from(depths)));
        }
        if let Some(errors) = self.per_base_errors {
            tags.push((Tag::from([b'c', b'E']), BufValue::from(errors)));
        }

        // Duplex consensus tags
        if let Some(depth) = self.ab_depth {
            tags.push((Tag::from([b'a', b'D']), to_smallest_signed_int(depth)));
        }
        if let Some(depth) = self.ba_depth {
            tags.push((Tag::from([b'b', b'D']), to_smallest_signed_int(depth)));
        }
        if let Some(depth) = self.ab_min_depth {
            tags.push((Tag::from([b'a', b'M']), to_smallest_signed_int(depth)));
        }
        if let Some(depth) = self.ba_min_depth {
            tags.push((Tag::from([b'b', b'M']), to_smallest_signed_int(depth)));
        }
        if let Some(rate) = self.ab_errors {
            tags.push((Tag::from([b'a', b'E']), BufValue::from(rate)));
        }
        if let Some(rate) = self.ba_errors {
            tags.push((Tag::from([b'b', b'E']), BufValue::from(rate)));
        }
        if let Some(count) = self.ab_error_count {
            tags.push((Tag::from([b'a', b'E']), to_smallest_signed_int(count)));
        }
        if let Some(count) = self.ba_error_count {
            tags.push((Tag::from([b'b', b'E']), to_smallest_signed_int(count)));
        }

        tags
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::sam::alignment::record::Cigar;

    #[test]
    fn test_new_builder() {
        let builder = SamBuilder::new();
        assert!(builder.is_empty());
        assert!(!builder.header.reference_sequences().is_empty());
    }

    #[test]
    fn test_add_frag_unmapped() {
        let mut builder = SamBuilder::new_unmapped();
        let rec = builder.add_frag().name("read1").bases("ACGT").build();

        assert_eq!(rec.name().map(std::convert::AsRef::as_ref), Some(b"read1".as_ref()));
        assert!(rec.flags().is_unmapped());
        assert_eq!(builder.len(), 1);
    }

    #[test]
    fn test_add_frag_mapped() {
        let mut builder = SamBuilder::new();
        let rec = builder.add_frag().name("read1").bases("ACGT").contig(0).start(100).build();

        assert!(!rec.flags().is_unmapped());
        assert_eq!(rec.reference_sequence_id(), Some(0));
        assert_eq!(rec.alignment_start(), Some(Position::try_from(100).unwrap()));
    }

    #[test]
    fn test_add_pair_unmapped() {
        let mut builder = SamBuilder::new_unmapped();
        let (read_one, read_two) =
            builder.add_pair().name("pair1").bases1("AAAA").bases2("CCCC").build();

        assert!(read_one.flags().is_first_segment());
        assert!(read_two.flags().is_last_segment());
        assert!(read_one.flags().is_unmapped());
        assert!(read_two.flags().is_unmapped());
        assert_eq!(builder.len(), 2);
    }

    #[test]
    fn test_add_pair_mapped() {
        let mut builder = SamBuilder::new();
        let (read_one, read_two) = builder
            .add_pair()
            .name("pair1")
            .bases1("AAAA")
            .bases2("CCCC")
            .start1(100)
            .start2(200)
            .build();

        assert!(!read_one.flags().is_unmapped());
        assert!(!read_two.flags().is_unmapped());
        assert_eq!(read_one.alignment_start(), Some(Position::try_from(100).unwrap()));
        assert_eq!(read_two.alignment_start(), Some(Position::try_from(200).unwrap()));
    }

    #[test]
    fn test_add_pair_with_attrs() {
        let mut builder = SamBuilder::new_unmapped();
        let (read_one, read_two) =
            builder.add_pair().attr("MI", "test_umi").attr("RX", "ACGT").build();

        let mi_tag = Tag::new(b'M', b'I');
        assert!(read_one.data().get(&mi_tag).is_some());
        assert!(read_two.data().get(&mi_tag).is_some());
    }

    #[test]
    fn test_sequential_names() {
        let mut builder = SamBuilder::new_unmapped();
        let first = builder.add_frag().bases("ACGT").build();
        let second = builder.add_frag().bases("ACGT").build();
        let third = builder.add_frag().bases("ACGT").build();

        assert_eq!(first.name().map(std::convert::AsRef::as_ref), Some(b"0000".as_ref()));
        assert_eq!(second.name().map(std::convert::AsRef::as_ref), Some(b"0001".as_ref()));
        assert_eq!(third.name().map(std::convert::AsRef::as_ref), Some(b"0002".as_ref()));
    }

    #[test]
    fn test_write_to_temp_file() {
        let mut builder = SamBuilder::new();
        let _ = builder.add_frag().name("read1").bases("ACGT").contig(0).start(100).build();

        let temp = builder.to_temp_file().unwrap();
        assert!(temp.path().exists());
    }

    #[test]
    fn test_parse_cigar() {
        let ops = parse_cigar("10M2I5M");
        assert_eq!(ops.len(), 3);
        assert_eq!(ops[0].kind(), Kind::Match);
        assert_eq!(ops[0].len(), 10);
        assert_eq!(ops[1].kind(), Kind::Insertion);
        assert_eq!(ops[1].len(), 2);
        assert_eq!(ops[2].kind(), Kind::Match);
        assert_eq!(ops[2].len(), 5);
    }

    #[test]
    fn test_strand_setting() {
        let mut builder = SamBuilder::new();
        let rec =
            builder.add_frag().bases("ACGT").contig(0).start(100).strand(Strand::Minus).build();

        assert!(rec.flags().is_reverse_complemented());
    }

    #[test]
    fn test_pair_strands() {
        let mut builder = SamBuilder::new();
        let (read_one, read_two) = builder
            .add_pair()
            .start1(100)
            .start2(200)
            .strand1(Strand::Minus)
            .strand2(Strand::Plus)
            .build();

        assert!(read_one.flags().is_reverse_complemented());
        assert!(!read_two.flags().is_reverse_complemented());
        assert!(!read_one.flags().is_mate_reverse_complemented());
        assert!(read_two.flags().is_mate_reverse_complemented());
    }

    #[test]
    fn test_template_length_calculation() {
        let mut builder = SamBuilder::new();
        let (read_one, read_two) =
            builder.add_pair().bases1("AAAA").bases2("CCCC").start1(100).start2(110).build();

        // R1 at 100-103, R2 at 110-113
        // Template spans 100-113 = 14bp
        assert_eq!(read_one.template_length(), 14);
        assert_eq!(read_two.template_length(), -14);
    }

    // ========================================================================
    // RecordBuilder tests
    // ========================================================================

    #[test]
    fn test_record_builder_basic() {
        let record = RecordBuilder::new()
            .name("read1")
            .sequence("ACGT")
            .qualities(&[30, 30, 30, 30])
            .build();

        assert_eq!(record.name().map(std::convert::AsRef::as_ref), Some(b"read1".as_ref()));
        assert_eq!(record.sequence().as_ref(), b"ACGT");
        assert_eq!(record.quality_scores().as_ref(), &[30, 30, 30, 30]);
    }

    #[test]
    fn test_record_builder_paired() {
        let record =
            RecordBuilder::new().name("read1").sequence("ACGT").first_segment(true).build();

        assert!(record.flags().is_segmented());
        assert!(record.flags().is_first_segment());
    }

    #[test]
    fn test_record_builder_mapped() {
        let record = RecordBuilder::new()
            .name("read1")
            .sequence("ACGT")
            .reference_sequence_id(0)
            .alignment_start(100)
            .cigar("4M")
            .mapping_quality(60)
            .build();

        assert_eq!(record.reference_sequence_id(), Some(0));
        assert_eq!(record.alignment_start(), Some(Position::try_from(100).unwrap()));
    }

    #[test]
    fn test_record_builder_with_tags() {
        let record = RecordBuilder::new()
            .name("read1")
            .sequence("ACGT")
            .tag("RG", "A")
            .tag("MI", "AAAA-CCCC/A")
            .build();

        let rg_tag = Tag::from([b'R', b'G']);
        let mi_tag = Tag::from([b'M', b'I']);

        assert!(record.data().get(&rg_tag).is_some());
        assert!(record.data().get(&mi_tag).is_some());
    }

    #[test]
    fn test_record_builder_auto_qualities() {
        let record = RecordBuilder::new().name("read1").sequence("ACGTACGT").build();

        // Should auto-generate Q30 qualities
        assert_eq!(record.quality_scores().as_ref(), &[30; 8]);
    }

    #[test]
    fn test_record_builder_secondary_supplementary() {
        let secondary = RecordBuilder::new().sequence("ACGT").secondary(true).build();
        assert!(secondary.flags().is_secondary());

        let supplementary = RecordBuilder::new().sequence("ACGT").supplementary(true).build();
        assert!(supplementary.flags().is_supplementary());
    }

    // ========================================================================
    // ConsensusTagsBuilder tests
    // ========================================================================

    #[test]
    fn test_consensus_tags_builder() {
        let record = RecordBuilder::new()
            .sequence("ACGT")
            .consensus_tags(ConsensusTagsBuilder::new().depth_max(10).depth_min(5).error_rate(0.01))
            .build();

        let cd_tag = Tag::from([b'c', b'D']);
        let cm_tag = Tag::from([b'c', b'M']);
        let ce_tag = Tag::from([b'c', b'E']);

        assert!(record.data().get(&cd_tag).is_some());
        assert!(record.data().get(&cm_tag).is_some());
        assert!(record.data().get(&ce_tag).is_some());
    }

    #[test]
    fn test_consensus_tags_per_base() {
        let tags = ConsensusTagsBuilder::new()
            .per_base_depths(&[10, 20, 30])
            .per_base_errors(&[1, 2, 3])
            .build();

        assert_eq!(tags.len(), 2);
    }

    // ========================================================================
    // New builder tests (Phase 0)
    // ========================================================================

    #[test]
    fn test_mapped_read_convenience() {
        let record = RecordBuilder::mapped_read()
            .name("read1")
            .sequence("ACGTACGT")
            .alignment_start(100)
            .build();

        assert_eq!(record.reference_sequence_id(), Some(0));
        assert_eq!(record.mapping_quality().map(u8::from), Some(60));
        // CIGAR should be auto-generated as 8M
        assert!(!record.cigar().is_empty());
    }

    #[test]
    fn test_record_builder_auto_cigar() {
        let record = RecordBuilder::new().sequence("ACGTACGT").build();

        // CIGAR should be auto-generated as 8M when sequence is set
        let ops: Vec<_> = record.cigar().iter().map(|r| r.unwrap()).collect();
        assert_eq!(ops.len(), 1);
        assert_eq!(ops[0].kind(), Kind::Match);
        assert_eq!(ops[0].len(), 8);
    }

    #[test]
    fn test_record_builder_explicit_cigar_overrides_auto() {
        let record = RecordBuilder::new().sequence("ACGTACGT").cigar("4M2I2M").build();

        let ops: Vec<_> = record.cigar().iter().map(|r| r.unwrap()).collect();
        assert_eq!(ops.len(), 3);
        assert_eq!(ops[0].kind(), Kind::Match);
        assert_eq!(ops[0].len(), 4);
        assert_eq!(ops[1].kind(), Kind::Insertion);
        assert_eq!(ops[1].len(), 2);
    }

    #[test]
    fn test_record_pair_builder_basic() {
        let (read_one, read_two) = RecordPairBuilder::new()
            .name("pair1")
            .r1_sequence("AAAA")
            .r2_sequence("TTTT")
            .r1_start(100)
            .r2_start(200)
            .build();

        assert!(read_one.flags().is_first_segment());
        assert!(read_two.flags().is_last_segment());
        assert!(!read_one.flags().is_reverse_complemented());
        assert!(read_two.flags().is_reverse_complemented()); // FR default
    }

    #[test]
    fn test_record_pair_builder_fr_orientation() {
        let (read_one, read_two) = RecordPairBuilder::new()
            .r1_sequence("ACGT")
            .r2_sequence("TGCA")
            .r1_start(100)
            .r2_start(200)
            .build();

        // Default is FR: R1 forward, R2 reverse
        assert!(!read_one.flags().is_reverse_complemented());
        assert!(read_two.flags().is_reverse_complemented());

        // Mate flags should also be correct
        assert!(read_one.flags().is_mate_reverse_complemented());
        assert!(!read_two.flags().is_mate_reverse_complemented());
    }

    #[test]
    fn test_record_pair_builder_with_tags() {
        let (read_one, read_two) = RecordPairBuilder::new()
            .name("pair1")
            .r1_sequence("ACGT")
            .r2_sequence("TGCA")
            .r1_start(100)
            .r2_start(200)
            .tag("MI", "1/A")
            .tag("RX", "AAAA-TTTT")
            .build();

        let mi_tag = Tag::from([b'M', b'I']);
        let rx_tag = Tag::from([b'R', b'X']);
        assert!(read_one.data().get(&mi_tag).is_some());
        assert!(read_two.data().get(&mi_tag).is_some());
        assert!(read_one.data().get(&rx_tag).is_some());
        assert!(read_two.data().get(&rx_tag).is_some());
    }

    #[test]
    fn test_record_pair_builder_template_length() {
        let (read_one, read_two) = RecordPairBuilder::new()
            .r1_sequence("AAAA")
            .r2_sequence("CCCC")
            .r1_start(100)
            .r2_start(110)
            .build();

        // R1 at 100-103, R2 at 110-113
        // Template spans 100-113 = 14bp
        assert_eq!(read_one.template_length(), 14);
        assert_eq!(read_two.template_length(), -14);
    }

    #[test]
    fn test_record_pair_builder_unmapped() {
        let (read_one, read_two) = RecordPairBuilder::new()
            .name("unmapped_pair")
            .r1_sequence("ACGT")
            .r2_sequence("TGCA")
            .build();

        // When no start positions are set, reads are unmapped
        assert!(read_one.flags().is_unmapped());
        assert!(read_two.flags().is_unmapped());
        assert!(read_one.flags().is_mate_unmapped());
        assert!(read_two.flags().is_mate_unmapped());
    }

    #[test]
    fn test_record_pair_builder_per_read_tags() {
        let (read_one, read_two) = RecordPairBuilder::new()
            .name("pair1")
            .r1_sequence("ACGT")
            .r2_sequence("TGCA")
            .r1_start(100)
            .r2_start(200)
            .tag("MI", "UMI123") // Shared tag
            .r1_tag("MQ", 0i32) // R1-only: mate quality 0 (unmapped mate)
            .r2_tag("MQ", 60i32) // R2-only: mate quality 60
            .build();

        let mi_tag = Tag::from([b'M', b'I']);
        let mq_tag = Tag::from([b'M', b'Q']);

        // Both have MI tag
        assert!(read_one.data().get(&mi_tag).is_some());
        assert!(read_two.data().get(&mi_tag).is_some());

        // Each has different MQ value
        let read_one_mq = read_one.data().get(&mq_tag).expect("R1 should have MQ tag");
        let read_two_mq = read_two.data().get(&mq_tag).expect("R2 should have MQ tag");

        // Values can be stored as different int types, so check the value directly
        match read_one_mq {
            BufValue::Int8(v) => assert_eq!(*v, 0),
            BufValue::UInt8(v) => assert_eq!(*v, 0),
            BufValue::Int16(v) => assert_eq!(*v, 0),
            BufValue::UInt16(v) => assert_eq!(*v, 0),
            BufValue::Int32(v) => assert_eq!(*v, 0),
            BufValue::UInt32(v) => assert_eq!(*v, 0),
            _ => panic!("Expected integer type for MQ tag, got {read_one_mq:?}"),
        }
        match read_two_mq {
            BufValue::Int8(v) => assert_eq!(*v, 60),
            BufValue::UInt8(v) => assert_eq!(*v, 60),
            BufValue::Int16(v) => assert_eq!(*v, 60),
            BufValue::UInt16(v) => assert_eq!(*v, 60),
            BufValue::Int32(v) => assert_eq!(*v, 60),
            BufValue::UInt32(v) => assert_eq!(*v, 60),
            _ => panic!("Expected integer type for MQ tag, got {read_two_mq:?}"),
        }
    }

    #[test]
    fn test_record_pair_builder_secondary() {
        let (read_one, read_two) = RecordPairBuilder::new()
            .name("sec_pair")
            .r1_sequence("ACGT")
            .r2_sequence("TGCA")
            .r1_start(100)
            .r2_start(200)
            .secondary(true)
            .build();

        assert!(read_one.flags().is_secondary());
        assert!(read_two.flags().is_secondary());
        assert!(!read_one.flags().is_supplementary());
        assert!(!read_two.flags().is_supplementary());
    }

    #[test]
    fn test_record_pair_builder_supplementary() {
        let (read_one, read_two) = RecordPairBuilder::new()
            .name("sup_pair")
            .r1_sequence("ACGT")
            .r2_sequence("TGCA")
            .r1_start(100)
            .r2_start(200)
            .supplementary(true)
            .build();

        assert!(!read_one.flags().is_secondary());
        assert!(!read_two.flags().is_secondary());
        assert!(read_one.flags().is_supplementary());
        assert!(read_two.flags().is_supplementary());
    }

    #[test]
    fn test_duplex_consensus_tags() {
        let record = RecordBuilder::new()
            .sequence("ACGT")
            .consensus_tags(
                ConsensusTagsBuilder::new()
                    .depth_max(10)
                    .depth_min(5)
                    .error_rate(0.01)
                    .ab_depth(6)
                    .ba_depth(4)
                    .ab_min_depth(3)
                    .ba_min_depth(2)
                    .ab_errors(0.005)
                    .ba_errors(0.008),
            )
            .build();

        // Standard consensus tags
        let cd_tag = Tag::from([b'c', b'D']);
        let cm_tag = Tag::from([b'c', b'M']);
        let ce_tag = Tag::from([b'c', b'E']);
        assert!(record.data().get(&cd_tag).is_some());
        assert!(record.data().get(&cm_tag).is_some());
        assert!(record.data().get(&ce_tag).is_some());

        // Duplex-specific tags
        let ab_depth_tag = Tag::from([b'a', b'D']);
        let ba_depth_tag = Tag::from([b'b', b'D']);
        let ab_min_tag = Tag::from([b'a', b'M']);
        let ba_min_tag = Tag::from([b'b', b'M']);
        let ab_error_tag = Tag::from([b'a', b'E']);
        let ba_error_tag = Tag::from([b'b', b'E']);
        assert!(record.data().get(&ab_depth_tag).is_some());
        assert!(record.data().get(&ba_depth_tag).is_some());
        assert!(record.data().get(&ab_min_tag).is_some());
        assert!(record.data().get(&ba_min_tag).is_some());
        assert!(record.data().get(&ab_error_tag).is_some());
        assert!(record.data().get(&ba_error_tag).is_some());
    }
}

// ============================================================================
// Test Data Generators
// ============================================================================

/// Calculates the sequence length required for a CIGAR string.
///
/// This counts the number of bases that consume the query sequence:
/// M (match), I (insertion), S (soft clip), = (sequence match), X (mismatch).
///
/// # Examples
///
/// ```rust
/// use fgumi_sam::builder::cigar_seq_len;
/// assert_eq!(cigar_seq_len("50M"), 50);
/// assert_eq!(cigar_seq_len("10S40M"), 50);
/// assert_eq!(cigar_seq_len("10H50M5S"), 55);
/// assert_eq!(cigar_seq_len("10M5I30M5D10M"), 55); // Deletions don't consume query
/// ```
#[must_use]
pub fn cigar_seq_len(cigar: &str) -> usize {
    let ops = parse_cigar(cigar);
    ops.iter()
        .filter(|op| {
            matches!(
                op.kind(),
                noodles::sam::alignment::record::cigar::op::Kind::Match
                    | noodles::sam::alignment::record::cigar::op::Kind::Insertion
                    | noodles::sam::alignment::record::cigar::op::Kind::SoftClip
                    | noodles::sam::alignment::record::cigar::op::Kind::SequenceMatch
                    | noodles::sam::alignment::record::cigar::op::Kind::SequenceMismatch
            )
        })
        .map(|op| op.len())
        .sum()
}

/// Creates a vector with `n` copies of an item.
///
/// This is useful for quickly generating repeated test data.
///
/// # Examples
///
/// ```rust
/// use fgumi_sam::builder::repeat_n;
/// let umis = repeat_n("AAAA".to_string(), 5);
/// assert_eq!(umis.len(), 5);
/// assert!(umis.iter().all(|u| u == "AAAA"));
/// ```
#[must_use]
pub fn repeat_n<T: Clone>(item: T, n: usize) -> Vec<T> {
    vec![item; n]
}

/// Generates quality scores with a specified mean and length.
///
/// All quality scores will be set to the same value (uniform distribution).
///
/// # Examples
///
/// ```rust
/// use fgumi_sam::builder::uniform_qualities;
/// let quals = uniform_qualities(30, 100);
/// assert_eq!(quals.len(), 100);
/// assert!(quals.iter().all(|&q| q == 30));
/// ```
#[must_use]
pub fn uniform_qualities(quality: u8, length: usize) -> Vec<u8> {
    vec![quality; length]
}

/// Generates quality scores that degrade linearly from start to end.
///
/// Useful for simulating read quality degradation along the length of a read.
///
/// # Panics
///
/// Panics if an interpolated quality value does not fit in `u8`, which cannot happen
/// when `start` and `end` are both valid `u8` quality scores.
///
/// # Examples
///
/// ```rust
/// use fgumi_sam::builder::degrading_qualities;
/// let quals = degrading_qualities(40, 20, 10);
/// assert_eq!(quals.len(), 10);
/// assert_eq!(quals[0], 40);
/// assert_eq!(quals[9], 20);
/// ```
#[must_use]
#[expect(
    clippy::cast_possible_truncation,
    reason = "f64 -> i64 cast is safe: values are interpolated between two u8 inputs (0-255)"
)]
#[expect(
    clippy::cast_precision_loss,
    reason = "usize -> f64 cast is acceptable: length values in tests are small enough that precision loss is negligible"
)]
pub fn degrading_qualities(start: u8, end: u8, length: usize) -> Vec<u8> {
    if length == 0 {
        return Vec::new();
    }
    if length == 1 {
        return vec![start];
    }

    let start_f = f64::from(start);
    let end_f = f64::from(end);
    let step = (end_f - start_f) / (length - 1) as f64;

    (0..length)
        .map(|i| {
            let quality = (start_f + step * i as f64).round();
            u8::try_from(quality as i64).expect("quality score fits in u8")
        })
        .collect()
}

/// Creates a temporary FASTA file with the given sequences.
///
/// Each tuple in `sequences` contains (name, sequence).
///
/// # Examples
///
/// ```rust,no_run
/// use fgumi_sam::builder::create_test_fasta;
/// let fasta = create_test_fasta(&[
///     ("chr1", "ACGTACGTACGT"),
///     ("chr2", "GGGGCCCCAAAA"),
/// ]).unwrap();
/// assert!(fasta.path().exists());
/// ```
///
/// # Errors
///
/// Returns an error if the temporary file cannot be created or written to.
pub fn create_test_fasta(sequences: &[(&str, &str)]) -> std::io::Result<NamedTempFile> {
    use std::io::Write;
    let mut file = NamedTempFile::new()?;
    for (name, seq) in sequences {
        writeln!(file, ">{name}")?;
        writeln!(file, "{seq}")?;
    }
    file.flush()?;
    Ok(file)
}

/// Creates a simple two-chromosome test FASTA file.
///
/// Creates chr1 with "ACGTACGTACGT" and chr2 with "GGGGCCCCAAAA".
/// This matches the default test FASTA used in multiple test modules.
///
/// # Examples
///
/// ```rust,no_run
/// use fgumi_sam::builder::create_default_test_fasta;
/// let fasta = create_default_test_fasta().unwrap();
/// assert!(fasta.path().exists());
/// ```
///
/// # Errors
///
/// Returns an error if the temporary file cannot be created or written to.
pub fn create_default_test_fasta() -> std::io::Result<NamedTempFile> {
    create_test_fasta(&[("chr1", "ACGTACGTACGT"), ("chr2", "GGGGCCCCAAAA")])
}

#[cfg(test)]
mod generator_tests {
    use super::*;

    #[test]
    fn test_repeat_n() {
        let items = repeat_n(42, 5);
        assert_eq!(items, vec![42, 42, 42, 42, 42]);
    }

    #[test]
    fn test_repeat_n_strings() {
        let items = repeat_n("ACGT".to_string(), 3);
        assert_eq!(items.len(), 3);
        assert!(items.iter().all(|s| s == "ACGT"));
    }

    #[test]
    fn test_uniform_qualities() {
        let quals = uniform_qualities(30, 10);
        assert_eq!(quals, vec![30; 10]);
    }

    #[test]
    fn test_degrading_qualities() {
        let quals = degrading_qualities(40, 20, 5);
        assert_eq!(quals.len(), 5);
        assert_eq!(quals[0], 40);
        assert_eq!(quals[4], 20);
        // Check intermediate values are monotonically decreasing
        for i in 0..quals.len() - 1 {
            assert!(quals[i] >= quals[i + 1]);
        }
    }

    #[test]
    fn test_degrading_qualities_edge_cases() {
        // Length 0
        assert_eq!(degrading_qualities(40, 20, 0), Vec::<u8>::new());

        // Length 1
        assert_eq!(degrading_qualities(40, 20, 1), vec![40]);

        // Same start and end
        let quals = degrading_qualities(30, 30, 5);
        assert!(quals.iter().all(|&q| q == 30));
    }

    #[test]
    fn test_create_test_fasta() {
        let fasta = create_test_fasta(&[("chr1", "ACGT"), ("chr2", "GGGG")]).unwrap();
        assert!(fasta.path().exists());

        // Read the file and verify contents
        let contents = std::fs::read_to_string(fasta.path()).unwrap();
        assert!(contents.contains(">chr1"));
        assert!(contents.contains("ACGT"));
        assert!(contents.contains(">chr2"));
        assert!(contents.contains("GGGG"));
    }

    #[test]
    fn test_create_default_test_fasta() {
        let fasta = create_default_test_fasta().unwrap();
        assert!(fasta.path().exists());

        let contents = std::fs::read_to_string(fasta.path()).unwrap();
        assert!(contents.contains(">chr1"));
        assert!(contents.contains("ACGTACGTACGT"));
        assert!(contents.contains(">chr2"));
        assert!(contents.contains("GGGGCCCCAAAA"));
    }
}

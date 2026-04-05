//! Core UMI extraction logic for building unmapped BAM records from FASTQ data.
//!
//! This module contains the pure-function core of the `extract` command, extracted
//! from the binary crate so that the pipeline can call it directly without spawning
//! a subprocess.
//!
//! # Overview
//!
//! The main entry point is [`make_raw_records`], which takes a [`FastqSet`] (with
//! segments already split by read structure) and an [`ExtractParams`] configuration,
//! and produces raw BAM record bytes suitable for writing to a BAM file or parsing
//! into `Template` objects.
//!
//! # Quality Encoding
//!
//! The [`QualityEncoding`] enum handles automatic detection and conversion between
//! Phred+33 (standard/Sanger) and Phred+64 (Illumina 1.3-1.7) quality encodings.

use anyhow::{Result, bail, ensure};
use bstr::{BString, ByteSlice};
use fgumi_raw_bam::UnmappedBamRecordBuilder;
use fgumi_raw_bam::fields::flags;
use noodles::sam::Header;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value as BufValue;
use noodles::sam::alignment::record_buf::{QualityScores, RecordBuf, Sequence};
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::ReadGroup;
use noodles::sam::header::record::value::map::header::{group_order, sort_order};
use noodles::sam::header::record::value::map::read_group::tag as rg_tag;
use noodles::sam::header::record::value::{
    Map as HeaderRecordMap,
    map::{Header as HeaderRecord, Tag as HeaderTag},
};

use crate::fastq::{FastqSegment, FastqSet};
use crate::template::{self, Template};

/// Number of FASTQ records to sample for quality encoding detection.
pub const QUALITY_DETECTION_SAMPLE_SIZE: usize = 400;

/// Quality encoding type detected from FASTQ quality scores.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QualityEncoding {
    /// Phred+33 (Sanger/Illumina 1.8+) -- the modern standard.
    Standard,
    /// Phred+64 (Illumina 1.3-1.7) -- legacy encoding.
    Illumina,
}

impl QualityEncoding {
    /// Convert quality scores to standard numeric format (Phred+33 offset removed).
    ///
    /// For `Standard` encoding, subtracts 33 from each quality byte.
    /// For `Illumina` encoding, subtracts 64 from each quality byte.
    /// Uses saturating subtraction to avoid underflow.
    #[must_use]
    pub fn to_standard_numeric(self, quals: &[u8]) -> Vec<u8> {
        match self {
            QualityEncoding::Standard => quals.iter().map(|&q| q.saturating_sub(33)).collect(),
            QualityEncoding::Illumina => quals.iter().map(|&q| q.saturating_sub(64)).collect(),
        }
    }

    /// Detect quality encoding from a sample of quality score records.
    ///
    /// Uses heuristics based on the observed quality score range:
    /// - Scores below 59: definitely Phred+33
    /// - Scores >= 64 with a reasonable range (max >= 75): likely Phred+64
    /// - Ambiguous ranges: defaults to Phred+33 (more common in modern data)
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - No records are provided
    /// - Quality scores fall outside the printable ASCII range (33-126)
    pub fn detect(records: &[Vec<u8>]) -> Result<Self> {
        if records.is_empty() {
            bail!("Cannot detect quality encoding: no records provided");
        }

        let mut min_qual = u8::MAX;
        let mut max_qual = u8::MIN;
        let mut total_bases = 0;

        for qual in records {
            if qual.is_empty() {
                continue;
            }
            for &q in qual {
                min_qual = min_qual.min(q);
                max_qual = max_qual.max(q);
                total_bases += 1;
            }
        }

        // If all reads were empty, default to Standard
        if total_bases == 0 {
            return Ok(QualityEncoding::Standard);
        }

        if min_qual < 33 || max_qual > 126 {
            bail!(
                "Invalid quality scores detected: range [{min_qual}, {max_qual}]. \
                Quality scores must be in the printable ASCII range (33-126)"
            );
        }

        if min_qual < 59 {
            Ok(QualityEncoding::Standard)
        } else if min_qual >= 64 {
            if max_qual >= 75 {
                Ok(QualityEncoding::Illumina)
            } else {
                Ok(QualityEncoding::Standard)
            }
        } else {
            // Ambiguous range (59-63), default to Phred+33
            Ok(QualityEncoding::Standard)
        }
    }
}

/// Configuration for UMI extraction, capturing all user-configurable options
/// needed by [`make_raw_records`].
///
/// This is the library-side equivalent of the CLI `ExtractConfig`, renamed to
/// avoid confusion with CLI argument structs.
#[derive(Clone, Debug)]
pub struct ExtractParams {
    /// BAM read group ID to assign to all records.
    pub read_group_id: String,
    /// Two-byte BAM tag for storing UMI sequences (e.g., `RX`).
    pub umi_tag: [u8; 2],
    /// Two-byte BAM tag for storing cell barcode sequences (e.g., `CB`).
    pub cell_tag: [u8; 2],
    /// Optional two-byte BAM tag for UMI quality scores.
    pub umi_qual_tag: Option<[u8; 2]>,
    /// Optional two-byte BAM tag for cell barcode quality scores.
    pub cell_qual_tag: Option<[u8; 2]>,
    /// Optional two-byte BAM tag for storing all concatenated UMIs.
    pub single_tag: Option<[u8; 2]>,
    /// Whether to append UMI sequences to read names (e.g., `readname+ACGT`).
    pub annotate_read_names: bool,
    /// Whether to extract UMIs from the 8th colon-delimited field of read names.
    pub extract_umis_from_read_names: bool,
    /// Whether to store sample barcode quality scores in the `QT` tag.
    pub store_sample_barcode_qualities: bool,
}

/// Optional read group metadata fields beyond the required sample/library.
///
/// All fields default to `None` (omitted from the header). The `platform` field
/// defaults to `"illumina"` when `None`, matching SAM spec conventions.
#[derive(Clone, Debug, Default)]
pub struct ReadGroupMetadata {
    /// Sequencing platform (default: `"illumina"`).
    pub platform: Option<String>,
    /// Platform unit (e.g. `flowcell-barcode.lane.sample-barcode`).
    pub platform_unit: Option<String>,
    /// Platform model (e.g. `miseq`, `hiseq2500`).
    pub platform_model: Option<String>,
    /// Library or sample barcode sequence.
    pub barcode: Option<String>,
    /// Sequencing center name.
    pub sequencing_center: Option<String>,
    /// Predicted median insert size.
    pub predicted_insert_size: Option<u32>,
    /// Read group description.
    pub description: Option<String>,
    /// Date the run was produced.
    pub run_date: Option<String>,
}

/// Build a minimal unmapped BAM header with sort order, read group, and optional
/// metadata fields.
///
/// This is the shared header construction logic used by both the `extract` command
/// and the pipeline's in-process extract path. It produces a header with:
/// - `SO:unsorted GO:query`
/// - A read group with sample, library, platform, and any additional metadata
///
/// The caller is responsible for adding `@PG` records and any extra comments.
///
/// # Errors
///
/// Returns an error if the header builder fails.
pub fn build_unmapped_bam_header(
    read_group_id: &str,
    sample: &str,
    library: &str,
    metadata: &ReadGroupMetadata,
) -> Result<Header> {
    let mut header = Header::builder();

    let HeaderTag::Other(so_tag) = HeaderTag::from([b'S', b'O']) else { unreachable!() };
    let HeaderTag::Other(go_tag) = HeaderTag::from([b'G', b'O']) else { unreachable!() };
    let map = HeaderRecordMap::<HeaderRecord>::builder()
        .insert(so_tag, sort_order::UNSORTED)
        .insert(go_tag, group_order::QUERY)
        .build()?;
    header = header.set_header(map);

    let platform = metadata.platform.clone().unwrap_or_else(|| String::from("illumina"));
    let mut rg = Map::<ReadGroup>::builder()
        .insert(rg_tag::SAMPLE, sample.to_owned())
        .insert(rg_tag::LIBRARY, library.to_owned())
        .insert(rg_tag::PLATFORM, platform);

    if let Some(ref v) = metadata.platform_unit {
        rg = rg.insert(rg_tag::PLATFORM_UNIT, v.clone());
    }
    if let Some(ref v) = metadata.platform_model {
        rg = rg.insert(rg_tag::PLATFORM_MODEL, v.clone());
    }
    if let Some(ref v) = metadata.barcode {
        rg = rg.insert(rg_tag::BARCODE, v.clone());
    }
    if let Some(ref v) = metadata.sequencing_center {
        rg = rg.insert(rg_tag::SEQUENCING_CENTER, v.clone());
    }
    if let Some(v) = metadata.predicted_insert_size {
        rg = rg.insert(rg_tag::PREDICTED_MEDIAN_INSERT_SIZE, v.to_string());
    }
    if let Some(ref v) = metadata.description {
        rg = rg.insert(rg_tag::DESCRIPTION, v.clone());
    }
    if let Some(ref v) = metadata.run_date {
        rg = rg.insert(rg_tag::PRODUCED_AT, v.clone());
    }

    header = header.add_read_group(read_group_id.to_owned(), rg.build()?);

    Ok(header.build())
}

/// Output of [`make_raw_records`]: raw BAM record bytes and a record count.
///
/// The `data` field contains BAM records with `block_size` prefixes, ready for
/// writing to a BAM file via a raw writer.
pub struct RawExtractedRecords {
    /// BAM records with `block_size` prefixes, ready for the compress step.
    pub data: Vec<u8>,
    /// Number of BAM records in this batch.
    pub num_records: u64,
}

/// Joins byte slices with a separator, pre-allocating capacity.
///
/// Returns an empty `BString` if the iterator yields no items.
///
/// # Examples
///
/// ```
/// use fgumi_lib::extract::join_bytes_with_separator;
/// use bstr::BString;
///
/// let result = join_bytes_with_separator(
///     [b"ACG".as_slice(), b"TTA".as_slice()].into_iter(),
///     b'-',
/// );
/// assert_eq!(result, BString::from("ACG-TTA"));
/// ```
pub fn join_bytes_with_separator<'a>(
    segments: impl Iterator<Item = &'a [u8]>,
    separator: u8,
) -> BString {
    let segments: Vec<&[u8]> = segments.collect();
    if segments.is_empty() {
        return BString::default();
    }
    let total_len: usize = segments.iter().map(|s| s.len()).sum();
    let capacity = total_len + segments.len().saturating_sub(1);
    let mut result = Vec::with_capacity(capacity);
    for (i, seg) in segments.iter().enumerate() {
        if i > 0 {
            result.push(separator);
        }
        result.extend_from_slice(seg);
    }
    BString::from(result)
}

/// Extracts the read name and optionally a UMI from the 8th colon-delimited field.
///
/// The read name is cleaned by:
/// 1. Removing the leading `@` prefix (if present)
/// 2. Truncating at the first space (removing comment/description fields)
///
/// If `extract_umis` is true, the 8th colon-delimited field of the read name is
/// extracted as the UMI. Any `+` characters in the UMI are normalized to `-`.
///
/// # Returns
///
/// A tuple of `(read_name_bytes, optional_umi_bytes)`.
#[must_use]
pub fn extract_read_name_and_umi(header: &[u8], extract_umis: bool) -> (Vec<u8>, Option<Vec<u8>>) {
    // Remove @ prefix if present
    let name_bytes = if header.starts_with(b"@") { &header[1..] } else { header };

    // Split on space to get just the name part
    let name_part = name_bytes.find_byte(b' ').map_or(name_bytes, |pos| &name_bytes[..pos]);

    if !extract_umis {
        return (name_part.to_vec(), None);
    }

    // Count colons to find 8th field (UMI)
    let parts: Vec<&[u8]> = name_part.split(|&b| b == b':').collect();

    if parts.len() >= 8 {
        let mut umi = parts[7].to_vec();
        if !umi.is_empty() {
            // Normalize '+' to '-' in UMI
            for byte in &mut umi {
                if *byte == b'+' {
                    *byte = b'-';
                }
            }
            return (name_part.to_vec(), Some(umi));
        }
    }

    (name_part.to_vec(), None)
}

/// Collected barcode sequences, qualities, and final read name for a single template.
///
/// Produced by [`extract_barcodes`] and consumed by both [`make_raw_records`] and
/// [`extract_template_and_fastq`] to avoid duplicating the barcode collection and
/// UMI merge logic.
pub struct ExtractedBarcodes<'a> {
    /// Joined cell barcode sequences (hyphen-separated).
    pub cell_barcode: BString,
    /// Joined sample barcode sequences (hyphen-separated).
    pub sample_barcode: BString,
    /// Final UMI after merging read-name and read-structure UMIs.
    pub umi: BString,
    /// Joined UMI quality strings (space-separated).
    pub umi_qual: BString,
    /// Joined cell barcode quality strings (space-separated).
    pub cell_qual: BString,
    /// Joined sample barcode quality strings (space-separated).
    pub sample_qual: BString,
    /// Read name bytes (cleaned, without `@` prefix or comments).
    /// Borrows from `read_name_bytes` when annotation is off; owned when annotation is on.
    pub final_read_name: std::borrow::Cow<'a, [u8]>,
    /// UMI extracted from the read name, if any (used to suppress UMI qual output).
    pub umi_from_name: Option<Vec<u8>>,
}

/// Extract barcodes, UMI, and final read name from a [`FastqSet`].
///
/// Collects cell barcodes, sample barcodes, molecular barcodes, and the UMI from
/// the read name (if requested). Merges UMIs from the read structure and read name
/// according to fgumi conventions. Returns an [`ExtractedBarcodes`] that both
/// [`make_raw_records`] and [`extract_template_and_fastq`] consume.
#[must_use]
pub fn extract_barcodes<'a>(
    read_set: &FastqSet,
    params: &ExtractParams,
    read_name_bytes: &'a [u8],
) -> ExtractedBarcodes<'a> {
    let cell_barcode =
        join_bytes_with_separator(read_set.cell_barcode_segments().map(|s| s.seq.as_slice()), b'-');
    let cell_qual = join_bytes_with_separator(
        read_set.cell_barcode_segments().map(|s| s.quals.as_slice()),
        b' ',
    );
    let sample_barcode = join_bytes_with_separator(
        read_set.sample_barcode_segments().map(|s| s.seq.as_slice()),
        b'-',
    );
    let sample_qual = join_bytes_with_separator(
        read_set.sample_barcode_segments().map(|s| s.quals.as_slice()),
        b' ',
    );
    let umi_bs = join_bytes_with_separator(
        read_set.molecular_barcode_segments().map(|s| s.seq.as_slice()),
        b'-',
    );
    let umi_qual = join_bytes_with_separator(
        read_set.molecular_barcode_segments().map(|s| s.quals.as_slice()),
        b' ',
    );

    let (_, umi_from_name) =
        extract_read_name_and_umi(&read_set.header, params.extract_umis_from_read_names);

    let umi: BString = match (umi_bs.is_empty(), &umi_from_name) {
        (true, Some(from_name)) => BString::from(from_name.as_slice()),
        (true, None) => BString::default(),
        (false, Some(from_name)) => {
            let mut combined = Vec::with_capacity(from_name.len() + 1 + umi_bs.len());
            combined.extend_from_slice(from_name);
            combined.push(b'-');
            combined.extend_from_slice(umi_bs.as_bytes());
            BString::from(combined)
        }
        (false, None) => umi_bs,
    };

    let final_read_name = if params.annotate_read_names && !umi.is_empty() {
        let mut name = Vec::with_capacity(read_name_bytes.len() + 1 + umi.len());
        name.extend_from_slice(read_name_bytes);
        name.push(b'+');
        name.extend_from_slice(umi.as_bytes());
        std::borrow::Cow::Owned(name)
    } else {
        std::borrow::Cow::Borrowed(read_name_bytes)
    };

    ExtractedBarcodes {
        cell_barcode,
        sample_barcode,
        umi,
        umi_qual,
        cell_qual,
        sample_qual,
        final_read_name,
        umi_from_name,
    }
}

/// Build raw BAM records from a [`FastqSet`] with UMI extraction applied.
///
/// This is the core extraction function that converts a combined FASTQ read set
/// (with segments already split by read structure) into raw BAM record bytes.
/// Each template segment in the read set produces one BAM record.
///
/// For paired-end data (2 template segments), records are flagged as paired with
/// appropriate first/last segment flags.
///
/// # Arguments
///
/// * `read_set` - Combined FASTQ segments from all input files for one template
/// * `encoding` - Quality encoding detected from the input FASTQ files
/// * `params` - Extraction configuration (tags, read group, etc.)
///
/// # Returns
///
/// A [`RawExtractedRecords`] containing the raw BAM bytes and record count.
///
/// # Errors
///
/// Returns an error if the read set contains no template segments.
pub fn make_raw_records(
    read_set: &FastqSet,
    encoding: QualityEncoding,
    params: &ExtractParams,
) -> Result<RawExtractedRecords> {
    let templates: Vec<&FastqSegment> = read_set.template_segments().collect();

    let read_name = String::from_utf8_lossy(&read_set.header);
    ensure!(!templates.is_empty(), "No template segments found for read: {read_name}");

    let (read_name_bytes, _) =
        extract_read_name_and_umi(&read_set.header, params.extract_umis_from_read_names);
    let barcodes = extract_barcodes(read_set, params, &read_name_bytes);

    let num_templates = templates.len();
    let mut builder = UnmappedBamRecordBuilder::new();
    let mut data = Vec::new();

    for (index, template) in templates.iter().enumerate() {
        let mut flag = flags::UNMAPPED;
        if num_templates == 2 {
            flag |= flags::PAIRED | flags::MATE_UNMAPPED;
            if index == 0 {
                flag |= flags::FIRST_SEGMENT;
            } else {
                flag |= flags::LAST_SEGMENT;
            }
        }

        if template.seq.is_empty() {
            builder.build_record(&barcodes.final_read_name, flag, b"N", &[2u8]);
        } else {
            let numeric_quals = encoding.to_standard_numeric(&template.quals);
            builder.build_record(&barcodes.final_read_name, flag, &template.seq, &numeric_quals);
        }

        builder.append_string_tag(b"RG", params.read_group_id.as_bytes());

        if !barcodes.cell_barcode.is_empty() {
            builder.append_string_tag(&params.cell_tag, barcodes.cell_barcode.as_bytes());
        }

        if !barcodes.cell_qual.is_empty() {
            if let Some(ref qt) = params.cell_qual_tag {
                builder.append_string_tag(qt, barcodes.cell_qual.as_bytes());
            }
        }

        if !barcodes.sample_barcode.is_empty() {
            builder.append_string_tag(b"BC", barcodes.sample_barcode.as_bytes());
        }

        if params.store_sample_barcode_qualities && !barcodes.sample_qual.is_empty() {
            builder.append_string_tag(b"QT", barcodes.sample_qual.as_bytes());
        }

        if !barcodes.umi.is_empty() {
            builder.append_string_tag(&params.umi_tag, barcodes.umi.as_bytes());

            if let Some(ref st) = params.single_tag {
                builder.append_string_tag(st, barcodes.umi.as_bytes());
            }

            if barcodes.umi_from_name.is_none() && !barcodes.umi_qual.is_empty() {
                if let Some(ref qt) = params.umi_qual_tag {
                    builder.append_string_tag(qt, barcodes.umi_qual.as_bytes());
                }
            }
        }

        builder.write_with_block_size(&mut data);
        builder.clear();
    }

    Ok(RawExtractedRecords { data, num_records: num_templates as u64 })
}

/// Output of [`extract_template_and_fastq`]: an unmapped [`Template`] for zipper merge
/// and interleaved FASTQ text for the aligner.
pub struct ExtractedTemplateAndFastq {
    /// Unmapped template with UMI/barcode tags set, for zipper merge.
    pub template: Template,
    /// Interleaved FASTQ text (R1 then R2) for piping to the aligner stdin.
    /// Each record is formatted as `@name/N\nSEQ\n+\nQUAL\n`.
    pub fastq_bytes: Vec<u8>,
}

/// Build an unmapped [`Template`] and interleaved FASTQ text from a [`FastqSet`].
///
/// This is the in-process alternative to running `fgumi extract` followed by
/// `fgumi fastq`. It produces both outputs in a single pass:
///
/// 1. A `Template` containing unmapped `RecordBuf` objects with UMI/barcode tags
///    (identical to what `fgumi extract` would produce as a BAM file).
/// 2. Interleaved FASTQ text suitable for piping to an aligner's stdin.
///
/// The FASTQ text uses the same read name as the BAM records (including UMI
/// annotation if `params.annotate_read_names` is set) and appends `/1` or `/2`
/// suffixes for paired-end reads.
///
/// # Arguments
///
/// * `read_set` - Combined FASTQ segments from all input files for one template
/// * `encoding` - Quality encoding detected from the input FASTQ files
/// * `params` - Extraction configuration (tags, read group, etc.)
///
/// # Errors
///
/// Returns an error if the read set contains no template segments.
// single-pass extraction must handle all barcode/tag/FASTQ logic together
#[allow(clippy::too_many_lines)]
pub fn extract_template_and_fastq(
    read_set: &FastqSet,
    encoding: QualityEncoding,
    params: &ExtractParams,
) -> Result<ExtractedTemplateAndFastq> {
    let templates: Vec<&FastqSegment> = read_set.template_segments().collect();

    let read_name_display = String::from_utf8_lossy(&read_set.header);
    ensure!(!templates.is_empty(), "No template segments found for read: {read_name_display}");

    let (read_name_bytes, _) =
        extract_read_name_and_umi(&read_set.header, params.extract_umis_from_read_names);
    let barcodes = extract_barcodes(read_set, params, &read_name_bytes);

    let num_templates = templates.len();
    let mut template_builder = template::Builder::default();
    let mut fastq_bytes = Vec::new();

    for (index, tmpl_seg) in templates.iter().enumerate() {
        let mut sam_flags = Flags::UNMAPPED;
        if num_templates == 2 {
            sam_flags |= Flags::SEGMENTED | Flags::MATE_UNMAPPED;
            if index == 0 {
                sam_flags |= Flags::FIRST_SEGMENT;
            } else {
                sam_flags |= Flags::LAST_SEGMENT;
            }
        }

        let mut record = RecordBuf::default();
        *record.name_mut() = Some(BString::from(barcodes.final_read_name.as_ref()));
        *record.flags_mut() = sam_flags;

        // Build seq/qual, write FASTQ from them FIRST, then move into RecordBuf
        let (seq_bytes, qual_values): (Vec<u8>, Vec<u8>) = if tmpl_seg.seq.is_empty() {
            (b"N".to_vec(), vec![2u8])
        } else {
            (tmpl_seg.seq.clone(), encoding.to_standard_numeric(&tmpl_seg.quals))
        };

        // Build FASTQ text before moving seq/qual into the RecordBuf
        fastq_bytes.push(b'@');
        fastq_bytes.extend_from_slice(&barcodes.final_read_name);
        if num_templates == 2 {
            if index == 0 {
                fastq_bytes.extend_from_slice(b"/1");
            } else {
                fastq_bytes.extend_from_slice(b"/2");
            }
        }
        fastq_bytes.push(b'\n');

        fastq_bytes.extend_from_slice(&seq_bytes);
        fastq_bytes.extend_from_slice(b"\n+\n");

        if tmpl_seg.seq.is_empty() {
            fastq_bytes.push(b'#'); // Q2 in Phred+33
        } else {
            for &q in &qual_values {
                fastq_bytes.push(q + 33);
            }
        }
        fastq_bytes.push(b'\n');

        // Move seq/qual into RecordBuf (no clone needed)
        *record.sequence_mut() = Sequence::from(seq_bytes);
        *record.quality_scores_mut() = QualityScores::from(qual_values);

        let data = record.data_mut();

        let rg_tag = Tag::new(b'R', b'G');
        data.insert(rg_tag, BufValue::from(params.read_group_id.as_str()));

        if !barcodes.cell_barcode.is_empty() {
            let tag = Tag::new(params.cell_tag[0], params.cell_tag[1]);
            data.insert(tag, BufValue::from(barcodes.cell_barcode.to_str_lossy().as_ref()));
        }
        if !barcodes.cell_qual.is_empty() {
            if let Some(ref qt) = params.cell_qual_tag {
                let tag = Tag::new(qt[0], qt[1]);
                data.insert(tag, BufValue::from(barcodes.cell_qual.to_str_lossy().as_ref()));
            }
        }

        if !barcodes.sample_barcode.is_empty() {
            let tag = Tag::new(b'B', b'C');
            data.insert(tag, BufValue::from(barcodes.sample_barcode.to_str_lossy().as_ref()));
        }
        if params.store_sample_barcode_qualities && !barcodes.sample_qual.is_empty() {
            let tag = Tag::new(b'Q', b'T');
            data.insert(tag, BufValue::from(barcodes.sample_qual.to_str_lossy().as_ref()));
        }

        if !barcodes.umi.is_empty() {
            let tag = Tag::new(params.umi_tag[0], params.umi_tag[1]);
            data.insert(tag, BufValue::from(barcodes.umi.to_str_lossy().as_ref()));

            if let Some(ref st) = params.single_tag {
                let tag = Tag::new(st[0], st[1]);
                data.insert(tag, BufValue::from(barcodes.umi.to_str_lossy().as_ref()));
            }

            if barcodes.umi_from_name.is_none() && !barcodes.umi_qual.is_empty() {
                if let Some(ref qt) = params.umi_qual_tag {
                    let tag = Tag::new(qt[0], qt[1]);
                    data.insert(tag, BufValue::from(barcodes.umi_qual.to_str_lossy().as_ref()));
                }
            }
        }

        template_builder.push(record)?;
    }

    let template = template_builder.build()?;
    Ok(ExtractedTemplateAndFastq { template, fastq_bytes })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fastq::FastqSet;
    use read_structure::ReadStructure;
    use std::str::FromStr;

    /// Helper to build a simple `ExtractParams` for testing.
    fn test_params() -> ExtractParams {
        ExtractParams {
            read_group_id: "A".to_string(),
            umi_tag: *b"RX",
            cell_tag: *b"CB",
            umi_qual_tag: None,
            cell_qual_tag: None,
            single_tag: None,
            annotate_read_names: false,
            extract_umis_from_read_names: false,
            store_sample_barcode_qualities: false,
        }
    }

    #[test]
    fn test_quality_encoding_standard_detection() {
        // Typical Phred+33 qualities (ASCII 33-73 = Q0-Q40)
        let records = vec![b"IIIIIFFFF".to_vec(), b"HHHHH!!!!".to_vec()];
        let encoding = QualityEncoding::detect(&records).unwrap();
        assert_eq!(encoding, QualityEncoding::Standard);
    }

    #[test]
    fn test_quality_encoding_illumina_detection() {
        // Phred+64 qualities (ASCII 64-126)
        let records = vec![vec![100u8; 10], vec![110u8; 10]];
        let encoding = QualityEncoding::detect(&records).unwrap();
        assert_eq!(encoding, QualityEncoding::Illumina);
    }

    #[test]
    fn test_quality_encoding_empty_records() {
        let records: Vec<Vec<u8>> = vec![];
        assert!(QualityEncoding::detect(&records).is_err());
    }

    #[test]
    fn test_quality_encoding_all_empty_reads() {
        let records = vec![Vec::new(), Vec::new()];
        let encoding = QualityEncoding::detect(&records).unwrap();
        assert_eq!(encoding, QualityEncoding::Standard);
    }

    #[test]
    fn test_quality_encoding_to_standard_numeric() {
        let standard = QualityEncoding::Standard;
        assert_eq!(standard.to_standard_numeric(&[33, 53, 73]), vec![0, 20, 40]);

        let illumina = QualityEncoding::Illumina;
        assert_eq!(illumina.to_standard_numeric(&[64, 84, 104]), vec![0, 20, 40]);
    }

    #[test]
    fn test_join_bytes_with_separator() {
        let result =
            join_bytes_with_separator([b"ACG".as_slice(), b"TTA".as_slice()].into_iter(), b'-');
        assert_eq!(result, BString::from("ACG-TTA"));

        // Single segment
        let result = join_bytes_with_separator([b"ACG".as_slice()].into_iter(), b'-');
        assert_eq!(result, BString::from("ACG"));

        // Empty iterator
        let result = join_bytes_with_separator(std::iter::empty(), b'-');
        assert_eq!(result, BString::default());
    }

    #[test]
    fn test_extract_read_name_and_umi_no_extract() {
        let header = b"@INSTRUMENT:RUN:FLOWCELL:LANE:TILE:X:Y:UMI extra";
        let (name, umi) = extract_read_name_and_umi(header, false);
        assert_eq!(name, b"INSTRUMENT:RUN:FLOWCELL:LANE:TILE:X:Y:UMI");
        assert!(umi.is_none());
    }

    #[test]
    fn test_extract_read_name_and_umi_with_extract() {
        let header = b"@INSTRUMENT:RUN:FLOWCELL:LANE:TILE:X:Y:ACGT+TGCA";
        let (name, umi) = extract_read_name_and_umi(header, true);
        assert_eq!(name, b"INSTRUMENT:RUN:FLOWCELL:LANE:TILE:X:Y:ACGT+TGCA");
        // '+' should be normalized to '-'
        assert_eq!(umi, Some(b"ACGT-TGCA".to_vec()));
    }

    #[test]
    fn test_extract_read_name_and_umi_no_umi_field() {
        let header = b"@SHORT:NAME:ONLY";
        let (name, umi) = extract_read_name_and_umi(header, true);
        assert_eq!(name, b"SHORT:NAME:ONLY");
        assert!(umi.is_none());
    }

    #[test]
    fn test_make_raw_records_simple_paired() {
        let rs = ReadStructure::from_str("+T").unwrap();
        let r1 =
            FastqSet::from_record_with_structure(b"@read1", b"ACGT", b"IIII", &rs, &[]).unwrap();
        let r2 =
            FastqSet::from_record_with_structure(b"@read1", b"TGCA", b"FFFF", &rs, &[]).unwrap();
        let combined = FastqSet::combine_readsets(vec![r1, r2]);

        let result = make_raw_records(&combined, QualityEncoding::Standard, &test_params());
        assert!(result.is_ok());
        let batch = result.unwrap();
        assert_eq!(batch.num_records, 2);
        assert!(!batch.data.is_empty());
    }

    #[test]
    fn test_make_raw_records_with_umi() {
        let rs_umi = ReadStructure::from_str("4M+T").unwrap();
        let rs_template = ReadStructure::from_str("+T").unwrap();
        let r1 =
            FastqSet::from_record_with_structure(b"@read1", b"ACGTTTTT", b"IIIIFFFF", &rs_umi, &[])
                .unwrap();
        let r2 =
            FastqSet::from_record_with_structure(b"@read1", b"TGCA", b"FFFF", &rs_template, &[])
                .unwrap();
        let combined = FastqSet::combine_readsets(vec![r1, r2]);

        let params = test_params();
        let result = make_raw_records(&combined, QualityEncoding::Standard, &params);
        assert!(result.is_ok());
        let batch = result.unwrap();
        assert_eq!(batch.num_records, 2);
    }

    #[test]
    fn test_make_raw_records_no_templates_fails() {
        // A read set with only UMI segments and no template segments
        let rs = ReadStructure::from_str("+M").unwrap();
        let r1 =
            FastqSet::from_record_with_structure(b"@read1", b"ACGT", b"IIII", &rs, &[]).unwrap();

        let result = make_raw_records(&r1, QualityEncoding::Standard, &test_params());
        assert!(result.is_err());
    }
}

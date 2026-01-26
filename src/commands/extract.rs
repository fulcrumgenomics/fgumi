//! Extract UMIs from FASTQ files and create unmapped BAM.
//!
//! This module implements the `extract` command which reads FASTQ files, extracts UMI
//! sequences based on a read structure specification, and outputs an unmapped BAM file
//! with UMI sequences stored in the `RX` tag and template bases as the read sequence.
//!
//! # Read Structure
//!
//! The read structure uses fgbio-style notation (e.g., `8M12S+T` meaning 8bp molecular
//! barcode, 12bp to skip, then template bases to end of read).
//!
//! # Quality Encoding
//!
//! Automatically detects and handles both Phred+33 (standard) and Phred+64 (Illumina 1.3-1.7)
//! quality encodings.

use crate::commands::command::Command;
use crate::commands::common::{CompressionOptions, SchedulerOptions, ThreadingOptions};
use anyhow::{Result, bail, ensure};
use bstr::{BString, ByteSlice};
use clap::Parser;
use fgoxide::io::Io;
use fgumi_lib::bam_io::{BamWriter, create_bam_writer};
use fgumi_lib::fastq::FastqSegment;
use fgumi_lib::fastq::FastqSet;
use fgumi_lib::fastq::ReadSetIterator;
use fgumi_lib::grouper::FastqTemplate;
use fgumi_lib::logging::OperationTimer;
use fgumi_lib::progress::ProgressTracker;
use fgumi_lib::unified_pipeline::{
    FastqPipelineConfig, run_fastq_pipeline, serialize_bam_records_into,
};
use fgumi_lib::validation::validate_file_exists;
use log::{debug, info};
use noodles_bgzf::io::MultithreadedReader;

#[cfg(test)]
use fgumi_lib::bam_io::create_bam_reader;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record::data::field::Value;
use noodles::sam::alignment::record_buf::QualityScores;
use noodles::sam::header::Header;
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::ReadGroup;
use noodles::sam::header::record::value::map::builder::Builder;
use noodles::sam::header::record::value::map::header::group_order;
use noodles::sam::header::record::value::map::header::sort_order;
use noodles::sam::header::record::value::map::read_group::tag as rg_tag;
use noodles::sam::header::record::value::{
    Map as HeaderRecordMap,
    map::{Header as HeaderRecord, Tag as HeaderTag},
};
use read_structure::{ReadStructure, SegmentType};
use seq_io::fastq::Reader as FastqReader;
use seq_io::fastq::Record;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::str::FromStr;

const BUFFER_SIZE: usize = 1024 * 1024;
const QUALITY_DETECTION_SAMPLE_SIZE: usize = 400;

/// Compression format detected from file header
#[derive(Debug, Clone, Copy, PartialEq)]
enum CompressionFormat {
    /// BGZF format (blocked gzip, can be parallelized)
    Bgzf,
    /// Standard gzip format (cannot be parallelized)
    Gzip,
    /// Uncompressed file
    Plain,
}

/// Detect the compression format of a file by reading its header.
///
/// BGZF files are identified by:
/// - Gzip magic number (0x1f 0x8b)
/// - Deflate compression method (0x08)
/// - FEXTRA flag set (0x04)
/// - Extra field with SI1='B' (0x42), SI2='C' (0x43)
///
/// Regular gzip files have the magic number but don't have the BGZF extra field.
fn detect_compression_format(path: &Path) -> Result<CompressionFormat> {
    let mut file = File::open(path)?;
    let mut header = [0u8; 18];

    use std::io::Read;
    let bytes_read = file.read(&mut header)?;

    if bytes_read < 2 {
        return Ok(CompressionFormat::Plain);
    }

    // Check for gzip magic number
    if header[0] != 0x1f || header[1] != 0x8b {
        return Ok(CompressionFormat::Plain);
    }

    // It's gzip - check if it's BGZF
    if bytes_read >= 18 {
        // Check for deflate method (0x08) and FEXTRA flag (0x04)
        if header[2] == 0x08 && (header[3] & 0x04) != 0 {
            // Extra field starts at byte 10
            // XLEN is at bytes 10-11 (little-endian)
            let xlen = u16::from_le_bytes([header[10], header[11]]) as usize;
            if xlen >= 6 {
                // Check for BGZF signature: SI1='B', SI2='C'
                if header[12] == b'B' && header[13] == b'C' {
                    return Ok(CompressionFormat::Bgzf);
                }
            }
        }
    }

    Ok(CompressionFormat::Gzip)
}

/// Open a FASTQ file with automatic detection of compression format.
///
/// For BGZF-compressed files, uses noodles `MultithreadedReader` when threads > 1.
/// For regular gzip files, uses fgoxide (flate2/zlib-ng) single-threaded decompression.
/// For uncompressed files, opens directly.
///
/// # Arguments
/// * `path` - Path to the FASTQ file
/// * `threads` - Number of decompression threads (only used for BGZF)
///
/// # Returns
/// A boxed reader that implements `BufRead` + `Send`
fn open_fastq_reader(path: &Path, threads: usize) -> Result<Box<dyn BufRead + Send>> {
    let format = detect_compression_format(path)?;
    let fgio = Io::new(5, BUFFER_SIZE);

    match format {
        CompressionFormat::Bgzf if threads > 1 => {
            info!("Detected BGZF-compressed FASTQ, using {threads} decompression threads");
            let file = File::open(path)?;
            let worker_count = std::num::NonZero::new(threads).expect("threads > 1 checked above");
            let reader = MultithreadedReader::with_worker_count(worker_count, file);
            Ok(Box::new(BufReader::new(reader)))
        }
        CompressionFormat::Bgzf => {
            debug!("Detected BGZF-compressed FASTQ, using single-threaded decompression");
            Ok(fgio.new_reader(path)?)
        }
        CompressionFormat::Gzip => {
            debug!("Detected gzip-compressed FASTQ, using single-threaded decompression");
            Ok(fgio.new_reader(path)?)
        }
        CompressionFormat::Plain => {
            debug!("Detected uncompressed FASTQ");
            Ok(fgio.new_reader(path)?)
        }
    }
}

/// Quality encoding type
#[derive(Debug, Clone, Copy, PartialEq)]
enum QualityEncoding {
    Standard, // Phred+33 (Sanger)
    Illumina, // Phred+64 (Illumina 1.3-1.7)
}

impl QualityEncoding {
    /// Convert quality scores to standard numeric format (Phred+33)
    fn to_standard_numeric(self, quals: &[u8]) -> Vec<u8> {
        match self {
            QualityEncoding::Standard => quals.iter().map(|&q| q.saturating_sub(33)).collect(),
            QualityEncoding::Illumina => quals.iter().map(|&q| q.saturating_sub(64)).collect(),
        }
    }

    /// Detect encoding from sample records using robust heuristics
    ///
    /// This implements a more robust detection algorithm that:
    /// - Checks for quality scores in the Illumina 1.3-1.7 range (64-126)
    /// - Checks for quality scores in the Sanger/Illumina 1.8+ range (33-126)
    /// - Handles edge cases like empty reads or very short reads
    /// - Provides informative error messages for invalid encodings
    fn detect(records: &[Vec<u8>]) -> Result<Self> {
        if records.is_empty() {
            bail!("Cannot detect quality encoding: no records provided");
        }

        let mut min_qual = u8::MAX;
        let mut max_qual = u8::MIN;
        let mut total_bases = 0;

        // Collect statistics from all quality scores
        for qual in records {
            if qual.is_empty() {
                continue; // Skip empty reads
            }
            for &q in qual {
                min_qual = min_qual.min(q);
                max_qual = max_qual.max(q);
                total_bases += 1;
            }
        }

        // If all reads were empty, we can't detect encoding but we'll default to Standard
        if total_bases == 0 {
            return Ok(QualityEncoding::Standard);
        }

        // Quality scores should be printable ASCII (33-126)
        if min_qual < 33 || max_qual > 126 {
            bail!(
                "Invalid quality scores detected: range [{min_qual}, {max_qual}]. \
                Quality scores must be in the printable ASCII range (33-126)"
            );
        }

        // Detect encoding based on observed range
        // Phred+64 (Illumina 1.3-1.7) uses range 64-126
        // Phred+33 (Sanger/Illumina 1.8+) uses range 33-126
        // If we see scores below 59, it's definitely Phred+33
        // If we see scores in 59-63, it could be either (low quality Phred+33 or very low Phred+64)
        // If we only see scores >= 64, it's likely Phred+64 (but could be high-quality Phred+33)

        if min_qual < 59 {
            // Definitely Phred+33 (Sanger/Illumina 1.8+)
            Ok(QualityEncoding::Standard)
        } else if min_qual >= 64 {
            // Likely Phred+64 (Illumina 1.3-1.7), but warn if range suggests otherwise
            // Note: Modern data should be Phred+33, so this is increasingly rare
            if max_qual >= 75 {
                // Has a reasonable range for Phred+64
                Ok(QualityEncoding::Illumina)
            } else {
                // Narrow range, could be high-quality Phred+33
                // Default to Phred+33 for modern data
                Ok(QualityEncoding::Standard)
            }
        } else {
            // Ambiguous range (59-63). This is very unlikely in real data.
            // Low quality Phred+33 would be in this range (Q26-Q30)
            // Very low quality Phred+64 would also be here (Q-5 to Q-1)
            // Default to Phred+33 as it's more common and these are reasonable quality scores
            Ok(QualityEncoding::Standard)
        }
    }
}

/// Generates an unmapped BAM file from fastq files.  Takes in one or more fastq files (optionally
/// gzipped), each representing a different sequencing read (e.g. R1, R2, I1 or I2) and can use a set of read
/// structures to allocate bases in those reads to template reads, sample indices, unique molecular indices, cell
/// barcodes, or to designate bases to be skipped over.
///
/// Read structures are made up of `<number><operator>` pairs much like the CIGAR string in BAM files. Five kinds of
/// operators are recognized:
///
/// 1. `T` identifies a template read
/// 2. `B` identifies a sample barcode read
/// 3. `M` identifies a unique molecular index read
/// 4. `C` identifies a cell barcode read
/// 5. `S` identifies a set of bases that should be skipped or ignored
///
/// The last `<number><operator>` pair may be specified using a `+` sign instead of number to denote "all remaining
/// bases". This is useful if, e.g., FASTQs have been trimmed and contain reads of varying length.  For example
/// to convert a paired-end run with an index read and where the first 5 bases of R1 are a UMI and the second
/// five bases are monotemplate you might specify:
///
/// ```
/// --input r1.fq r2.fq i1.fq --read-structures 5M5S+T +T +B
/// ```
///
/// Alternative if you know your reads are of fixed length you could specify:
///
/// ```
/// --input r1.fq r2.fq i1.fq --read-structures 5M5S65T 75T 8B
/// ```
///
/// For more information on read structures see the
/// [Read Structure Wiki Page](https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures)
///
/// UMIs may be extracted from the read sequences, the read names, or both.  If `--extract-umis-from-read-names` is
/// specified, any UMIs present in the read names are extracted; read names are expected to be `:`-separated with
/// any UMIs present in the 8th field.  If this option is specified, the `--umi-qual-tag` option may not be used as
/// qualities are not available for UMIs in the read name. If UMI segments are present in the read structures those
/// will also be extracted.  If UMIs are present in both, the final UMIs are constructed by first taking the UMIs
/// from the read names, then adding a hyphen, then the UMIs extracted from the reads.
///
/// The same number of input files and read structures must be provided, with one exception: if supplying exactly
/// 1 or 2 fastq files, both of which are solely template reads, no read structures need be provided.
///
/// The output file will produce a BAM with reads in the same order as they appear in the fastq file.
///
/// # Design Note
///
/// This command combines the functionality of two separate fgbio tools:
/// 1. **`FastqToBam`** - Converts FASTQ files to unmapped BAM format
/// 2. **`ExtractUmisFromBam`** - Extracts UMIs from sequences/read names into BAM tags
///
/// By combining these operations, fgumi provides a more streamlined workflow that:
/// - Eliminates the need for intermediate BAM files
/// - Reduces I/O overhead
/// - Simplifies the command-line interface for common UMI-based workflows
///
/// The functionality is equivalent to running `fgbio FastqToBam | fgbio ExtractUmisFromBam`,
/// but performs both operations in a single pass over the input FASTQ files.
#[derive(Parser, Debug)]
#[command(
    name = "extract",
    author,
    version,
    about = "\x1b[38;5;30m[UMI EXTRACTION]\x1b[0m \x1b[36mExtract UMIs from FASTQ and create unmapped BAM\x1b[0m",
    long_about = r#"
Generates an unmapped BAM file from FASTQ files with UMI extraction.

Takes in one or more FASTQ files (optionally gzipped), each representing a different sequencing
read (e.g. R1, R2, I1 or I2) and can use a set of read structures to allocate bases in those
reads to template reads, sample indices, unique molecular indices, or to designate bases to be
skipped over.

Only template bases will be retained as read bases (stored in the `SEQ` field) as specified by
the read structure.

## Read Structures

Read structures are made up of `<number><operator>` pairs much like the CIGAR string in BAM files.
Five kinds of operators are recognized:

1. `T` identifies a template read
2. `B` identifies a sample barcode read
3. `M` identifies a unique molecular index read
4. `C` identifies a cell barcode read
5. `S` identifies a set of bases that should be skipped or ignored

The last `<number><operator>` pair may be specified using a `+` sign instead of number to denote
"all remaining bases". This is useful if, e.g., FASTQs have been trimmed and contain reads of
varying length.

For example, to convert a paired-end run with an index read and where the first 5 bases of R1 are
a UMI and the second five bases are monotemplate:

  fgumi extract --input r1.fq r2.fq i1.fq --read-structures 5M5S+T +T +B

Alternatively, if reads are fixed length:

  fgumi extract --input r1.fq r2.fq i1.fq --read-structures 5M5S65T 75T 8B

## UMI Extraction

A read structure should be provided for each read of a template. For paired end reads, two read
structures should be specified. The tags to store the molecular indices will be associated with
the molecular index segment(s) in the read structure based on the order specified. If only one
molecular index tag is given, then the molecular indices will be concatenated and stored in that
tag. In the resulting BAM file each end of a pair will contain the same molecular index tags and
values.

UMIs may be extracted from the read sequences, the read names, or both. If
`--extract-umis-from-read-names` is specified, any UMIs present in the read names are extracted;
read names are expected to be `:` -separated with any UMIs present in the 8th field. If UMI
segments are present in the read structures those will also be extracted. If UMIs are present in
both, the final UMIs are constructed by first taking the UMIs from the read names, then adding a
hyphen, then the UMIs extracted from the reads.
"#
)]
#[command(verbatim_doc_comment)]
#[allow(clippy::struct_excessive_bools)]
pub(crate) struct Extract {
    /// Input FASTQ files corresponding to each sequencing read (e.g. R1, I1, etc.)
    #[arg(long, short = 'i', required = true, num_args = 1..)]
    inputs: Vec<PathBuf>,

    /// Output BAM file to be written
    #[arg(long, short = 'o', required = true)]
    output: PathBuf,

    /// Read structures, one for each of the FASTQs (optional if 1-2 template-only FASTQs)
    #[arg(long, short = 'r', num_args = 0..)]
    read_structures: Vec<ReadStructure>,

    /// Tag in which to store molecular barcodes/UMIs
    #[arg(long, short = 'u', default_value = "RX")]
    umi_tag: String,

    /// Tag in which to store molecular barcode/UMI qualities
    #[arg(long, short = 'q')]
    umi_qual_tag: Option<String>,

    /// Tag in which to store the cellular barcodes
    #[arg(long, short = 'c', default_value = "CB")]
    cell_tag: String,

    /// Tag in which to store the cellular barcode qualities
    #[arg(long, short = 'C')]
    cell_qual_tag: Option<String>,

    /// Store the sample barcode qualities in the QT Tag
    #[arg(long, short = 'Q')]
    store_sample_barcode_qualities: bool,

    /// Extract UMI(s) from read names and prepend to UMIs from reads
    #[arg(long, short = 'n')]
    #[allow(clippy::struct_field_names)]
    extract_umis_from_read_names: bool,

    /// Annotate read names with UMIs (appends "+UMIs" to read names)
    #[arg(long, short = 'a')]
    annotate_read_names: bool,

    /// Single tag to store all concatenated UMIs (in addition to per-segment tags)
    #[arg(long, short = 's')]
    single_tag: Option<String>,

    /// Tag containing adapter clipping position to adjust (e.g. 'XT' from `MarkIlluminaAdapters`)
    #[arg(long)]
    clipping_attribute: Option<String>,

    /// Read group ID to use in the file header
    #[arg(long, default_value = "A")]
    read_group_id: String,

    /// The name of the sequenced sample
    #[arg(long, required = true)]
    sample: String,

    /// The name/ID of the sequenced library
    #[arg(long, required = true)]
    library: String,

    /// Library or Sample barcode sequence
    #[arg(long, short = 'b')]
    barcode: Option<String>,

    /// Sequencing Platform
    #[arg(long, default_value = "illumina")]
    platform: String,

    /// Platform unit (e.g. 'flowcell-barcode.lane.sample-barcode')
    #[arg(long)]
    platform_unit: Option<String>,

    /// Platform model to insert into the group header (ex. miseq, hiseq2500, hiseqX)
    #[arg(long)]
    platform_model: Option<String>,

    /// The sequencing center from which the data originated
    #[arg(long)]
    sequencing_center: Option<String>,

    /// Predicted median insert size, to insert into the read group header
    #[arg(long)]
    predicted_insert_size: Option<u32>,

    /// Description of the read group
    #[arg(long)]
    description: Option<String>,

    /// Comment(s) to include in the output file's header
    #[arg(long, num_args = 0..)]
    comment: Vec<String>,

    /// Date the run was produced, to insert into the read group header
    #[arg(long)]
    run_date: Option<String>,

    /// Threading options for parallel processing.
    #[command(flatten)]
    pub threading: ThreadingOptions,

    /// Compression options for output BAM.
    #[command(flatten)]
    pub compression: CompressionOptions,

    /// Scheduler and pipeline statistics options.
    #[command(flatten)]
    pub scheduler_opts: SchedulerOptions,
}

impl Extract {
    /// Get actual read structures (default to +T if none provided for 1-2 FASTQs)
    fn get_read_structures(&self) -> Result<Vec<ReadStructure>> {
        if self.read_structures.is_empty() && (1..=2).contains(&self.inputs.len()) {
            Ok(vec![ReadStructure::from_str("+T")?; self.inputs.len()])
        } else {
            Ok(self.read_structures.clone())
        }
    }

    /// Validate inputs
    fn validate(&self) -> Result<()> {
        let read_structures = self.get_read_structures()?;

        ensure!(
            self.inputs.len() == read_structures.len(),
            "input and read-structure must be supplied the same number of times."
        );

        let template_count: usize = read_structures
            .iter()
            .map(|rs| rs.segments_by_type(SegmentType::Template).count())
            .sum();
        ensure!(
            (1..=2).contains(&template_count),
            "read structures must contain 1-2 template reads total."
        );

        ensure!(
            !self.extract_umis_from_read_names || self.umi_qual_tag.is_none(),
            "Cannot extract UMI qualities when also extracting UMI from read names."
        );

        // Validate threads parameter (ThreadingOptions handles minimum of 1 internally)

        // Validate input files exist
        for input in &self.inputs {
            validate_file_exists(input, "Input FASTQ")?;
        }

        // Validate read structures are not empty
        for (i, rs) in read_structures.iter().enumerate() {
            ensure!(!rs.segments().is_empty(), "Read structure {} is empty", i + 1);
        }

        // Validate single_tag is 2 characters if specified
        if let Some(ref tag) = self.single_tag {
            ensure!(tag.len() == 2, "Single tag must be exactly 2 characters: {tag}");
            ensure!(tag != &self.umi_tag, "Single tag cannot be the same as umi-tag: {tag}");
        }

        // Validate clipping_attribute is 2 characters if specified
        if let Some(ref tag) = self.clipping_attribute {
            ensure!(tag.len() == 2, "Clipping attribute must be exactly 2 characters: {tag}");
        }

        // Note: clipping_attribute doesn't apply to FASTQ input (no existing clipping to adjust)
        // It's accepted for compatibility but has no effect in FASTQ→BAM mode

        Ok(())
    }

    /// Helper to conditionally add a tag/value pair to a read group
    ///
    /// If the value is Some, inserts the tag with the value into the read group builder.
    /// If the value is None, returns the builder unchanged.
    ///
    /// # Arguments
    /// * `rg` - The read group builder
    /// * `tag` - The tag to insert
    /// * `value` - Optional value to insert
    ///
    /// # Returns
    /// The read group builder, potentially with the tag added
    fn add_to_read_group(
        rg: Builder<ReadGroup>,
        tag: noodles::sam::header::record::value::map::tag::Other<rg_tag::Standard>,
        value: Option<&String>,
    ) -> Builder<ReadGroup> {
        if let Some(v) = value { rg.insert(tag, v.clone()) } else { rg }
    }

    /// Create SAM header
    fn create_header(&self, command_line: &str) -> Result<Header> {
        let mut header = Header::builder();

        // Sort and group order
        let HeaderTag::Other(so_tag) = HeaderTag::from([b'S', b'O']) else { unreachable!() };
        let HeaderTag::Other(go_tag) = HeaderTag::from([b'G', b'O']) else { unreachable!() };
        let map = HeaderRecordMap::<HeaderRecord>::builder()
            .insert(so_tag, sort_order::UNSORTED)
            .insert(go_tag, group_order::QUERY)
            .build()?;
        header = header.set_header(map);

        // Add comments
        for comment in &self.comment {
            header = header.add_comment(comment.clone());
        }

        // Create read group
        let mut rg = Map::<ReadGroup>::builder();
        rg = Self::add_to_read_group(rg, rg_tag::SAMPLE, Some(&self.sample.clone()));
        rg = Self::add_to_read_group(rg, rg_tag::LIBRARY, Some(&self.library.clone()));
        rg = Self::add_to_read_group(rg, rg_tag::BARCODE, self.barcode.as_ref());
        rg = Self::add_to_read_group(rg, rg_tag::PLATFORM, Some(&self.platform));
        rg = Self::add_to_read_group(rg, rg_tag::PLATFORM_UNIT, self.platform_unit.as_ref());
        rg = Self::add_to_read_group(rg, rg_tag::PLATFORM_MODEL, self.platform_model.as_ref());
        rg =
            Self::add_to_read_group(rg, rg_tag::SEQUENCING_CENTER, self.sequencing_center.as_ref());
        rg = Self::add_to_read_group(
            rg,
            rg_tag::PREDICTED_MEDIAN_INSERT_SIZE,
            self.predicted_insert_size.map(|i| i.to_string()).as_ref(),
        );
        rg = Self::add_to_read_group(rg, rg_tag::DESCRIPTION, self.description.as_ref());
        rg = Self::add_to_read_group(rg, rg_tag::PRODUCED_AT, self.run_date.as_ref());

        header = header.add_read_group(self.read_group_id.clone(), rg.build()?);

        // Add @PG record
        header = fgumi_lib::header::add_pg_to_builder(
            header,
            crate::version::VERSION.as_str(),
            command_line,
        )?;

        Ok(header.build())
    }

    /// Adds the given tag and value to the record's data.
    ///
    /// # Panics
    /// If the tag is not length two
    fn add_tag(
        data: &mut noodles::sam::alignment::record_buf::Data,
        tag: &[u8],
        value: &BString,
    ) -> Result<(), Box<dyn std::error::Error>> {
        assert!(tag.len() == 2, "Tag must have length two, found {}: {tag:?}", tag.len());
        // Create the tag
        let bytes = tag.as_bytes();
        let tag = Tag::from([bytes[0], bytes[1]]);
        // Create the value
        let bstr_value = value.as_bstr();
        // Insert it
        data.insert(tag, Value::String(bstr_value).try_into()?);
        Ok(())
    }

    /// Joins byte slices with a separator, pre-allocating capacity.
    /// Returns empty `BString` if iterator is empty.
    fn join_bytes_with_separator<'a>(
        segments: impl Iterator<Item = &'a [u8]>,
        separator: u8,
    ) -> BString {
        let segments: Vec<&[u8]> = segments.collect();
        if segments.is_empty() {
            return BString::default();
        }
        // Calculate total capacity needed
        let total_len: usize = segments.iter().map(|s| s.len()).sum();
        let capacity = total_len + segments.len().saturating_sub(1); // separators
        let mut result = Vec::with_capacity(capacity);
        for (i, seg) in segments.iter().enumerate() {
            if i > 0 {
                result.push(separator);
            }
            result.extend_from_slice(seg);
        }
        BString::from(result)
    }

    /// Extracts read name, optionally extracting UMI from the 8th colon-delimited field
    fn extract_read_name_and_umi(
        header: &Vec<u8>,
        extract_umis: bool,
    ) -> (Vec<u8>, Option<Vec<u8>>) {
        // Remove @ prefix if present
        let name_bytes = if header.starts_with(b"@") { &header[1..] } else { header.as_slice() };

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

    /// Validates that all read names match across the read sets
    fn validate_read_names_match(read_sets: &[FastqSet]) -> Result<()> {
        if read_sets.is_empty() {
            return Ok(());
        }

        // Extract the read name from the first header (removing @ prefix if present)
        let first_header = &read_sets[0].header;
        let first_name = if first_header.starts_with(b"@") {
            &first_header[1..]
        } else {
            first_header.as_slice()
        };

        // Strip space comments and /1, /2 suffixes for comparison
        let first_name_part = strip_read_suffix_extract(first_name);

        // Check that all other read sets have the same name
        for (i, read_set) in read_sets.iter().enumerate().skip(1) {
            let header = &read_set.header;
            let name = if header.starts_with(b"@") { &header[1..] } else { header.as_slice() };

            let name_part = strip_read_suffix_extract(name);

            if name_part != first_name_part {
                bail!(
                    "Read names do not match across FASTQs: '{}' vs '{}' (FASTQ index 0 vs {})",
                    String::from_utf8_lossy(first_name_part),
                    String::from_utf8_lossy(name_part),
                    i
                );
            }
        }

        Ok(())
    }

    /// Create SAM records from a read set
    #[allow(clippy::too_many_lines)]
    fn make_sam_records(
        &self,
        read_set: &FastqSet,
        encoding: QualityEncoding,
    ) -> Result<Vec<RecordBuf>> {
        let templates: Vec<&FastqSegment> = read_set.template_segments().collect();

        let read_name = String::from_utf8_lossy(&read_set.header);
        ensure!(!templates.is_empty(), "No template segments found for read: {read_name}");

        // Extract various barcode types as BString - use optimized join
        let cell_barcode_bs = Self::join_bytes_with_separator(
            read_set.cell_barcode_segments().map(|s| s.seq.as_slice()),
            b'-',
        );
        let cell_quals_bs = Self::join_bytes_with_separator(
            read_set.cell_barcode_segments().map(|s| s.quals.as_slice()),
            b' ',
        );
        let sample_barcode_bs = Self::join_bytes_with_separator(
            read_set.sample_barcode_segments().map(|s| s.seq.as_slice()),
            b'-',
        );
        let sample_quals_bs = Self::join_bytes_with_separator(
            read_set.sample_barcode_segments().map(|s| s.quals.as_slice()),
            b' ',
        );
        let umi_bs = Self::join_bytes_with_separator(
            read_set.molecular_barcode_segments().map(|s| s.seq.as_slice()),
            b'-',
        );
        let umi_qual_bs = Self::join_bytes_with_separator(
            read_set.molecular_barcode_segments().map(|s| s.quals.as_slice()),
            b' ',
        );

        // Extract UMI from read name if requested
        let (read_name, umi_from_name) =
            Self::extract_read_name_and_umi(&read_set.header, self.extract_umis_from_read_names);

        // Prepare final UMI as BString - avoid format! and unnecessary allocations
        let final_umi_bs: BString = match (umi_bs.is_empty(), &umi_from_name) {
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

        // Read group ID as BString
        let rgid_bs: BString = self.read_group_id.as_str().into();

        let mut records = Vec::new();

        for (index, template) in templates.iter().enumerate() {
            // Create record with proper API
            let mut rec = RecordBuf::default();

            // Set read name (optionally with UMI annotation)
            let final_read_name = if self.annotate_read_names && !final_umi_bs.is_empty() {
                format!("{}+{}", String::from_utf8_lossy(&read_name), final_umi_bs.as_bstr())
                    .into_bytes()
            } else {
                read_name.clone()
            };
            *rec.name_mut() = Some(final_read_name.into());

            // Set bases and qualities - if empty, substitute with single N @ Q2
            if template.seq.is_empty() {
                *rec.sequence_mut() = "N".as_bytes().to_vec().into();
                *rec.quality_scores_mut() = QualityScores::from(vec![2u8]);
            } else {
                *rec.sequence_mut() = template.seq.clone().into();
                let numeric_quals = encoding.to_standard_numeric(&template.quals);
                *rec.quality_scores_mut() = QualityScores::from(numeric_quals);
            }

            // Set flags for unmapped reads
            let mut flags = Flags::UNMAPPED;
            if templates.len() == 2 {
                flags |= Flags::SEGMENTED | Flags::MATE_UNMAPPED;
                if index == 0 {
                    flags |= Flags::FIRST_SEGMENT;
                } else {
                    flags |= Flags::LAST_SEGMENT;
                }
            }
            *rec.flags_mut() = flags;

            // Set mapping quality to 0 for unmapped reads (per SAM spec)
            *rec.mapping_quality_mut() = Some(
                noodles::sam::alignment::record::MappingQuality::new(0)
                    .expect("0 is a valid mapping quality"),
            );

            // Add tags using the mutable data accessor
            let data = rec.data_mut();

            // Read group
            data.insert(Tag::READ_GROUP, Value::String(rgid_bs.as_bstr()).try_into()?);

            // Cell barcode
            if !cell_barcode_bs.is_empty() {
                Self::add_tag(data, self.cell_tag.as_bytes(), &cell_barcode_bs).unwrap();
            }

            if !cell_quals_bs.is_empty() {
                if let Some(ref qual_tag) = self.cell_qual_tag {
                    Self::add_tag(data, qual_tag.as_bytes(), &cell_quals_bs).unwrap();
                }
            }

            // Sample barcode
            if !sample_barcode_bs.is_empty() {
                data.insert(
                    Tag::SAMPLE_BARCODE_SEQUENCE,
                    Value::String(sample_barcode_bs.as_bstr()).try_into()?,
                );
            }

            if self.store_sample_barcode_qualities && !sample_quals_bs.is_empty() {
                data.insert(
                    Tag::SAMPLE_BARCODE_QUALITY_SCORES,
                    Value::String(sample_quals_bs.as_bstr()).try_into()?,
                );
            }

            // UMI
            if !final_umi_bs.is_empty() {
                Self::add_tag(data, self.umi_tag.as_bytes(), &final_umi_bs).unwrap();

                // Single tag for all concatenated UMIs (if specified)
                if let Some(ref single_tag) = self.single_tag {
                    Self::add_tag(data, single_tag.as_bytes(), &final_umi_bs).unwrap();
                }

                // Only add UMI qualities if not extracted from read names
                if umi_from_name.is_none() && !umi_qual_bs.is_empty() {
                    if let Some(ref qual_tag) = self.umi_qual_tag {
                        Self::add_tag(data, qual_tag.as_bytes(), &umi_qual_bs).unwrap();
                    }
                }
            }

            records.push(rec);
        }

        Ok(records)
    }

    /// Process records in single-threaded mode
    ///
    /// Returns the number of records written.
    fn process_singlethreaded(
        &self,
        fq_iterators: &mut [ReadSetIterator],
        header: &Header,
        writer: &mut BamWriter,
        encoding: QualityEncoding,
    ) -> Result<u64> {
        let progress = ProgressTracker::new("Processed records").with_interval(1_000_000);
        let mut read_pair_count: usize = 0;

        loop {
            let mut next_read_sets = Vec::with_capacity(fq_iterators.len());
            for iter in fq_iterators.iter_mut() {
                if let Some(rec) = iter.next() {
                    next_read_sets.push(rec);
                } else {
                    break;
                }
            }

            if next_read_sets.is_empty() {
                break;
            }

            ensure!(next_read_sets.len() == fq_iterators.len(), "FASTQ sources out of sync");

            //  Validate read names match across all FASTQs
            Self::validate_read_names_match(&next_read_sets)?;

            let read_set = FastqSet::combine_readsets(next_read_sets);
            let sam_records = self.make_sam_records(&read_set, encoding)?;

            let num_records = sam_records.len();
            for record in &sam_records {
                writer.write_alignment_record(header, record)?;
            }

            read_pair_count += num_records;
            progress.log_if_needed(num_records as u64);
        }

        progress.log_final();
        Ok(read_pair_count as u64)
    }

    /// Process records using the 7-step pipeline for better thread scaling.
    ///
    /// This method uses the unified 7-step pipeline that provides better parallelism
    /// for FASTQ → BAM conversion by separating reading, decompression, grouping,
    /// processing, serialization, compression, and writing into distinct steps.
    /// Process records using the 7-step pipeline.
    ///
    /// Returns the number of records written.
    fn process_with_pipeline(
        &self,
        header: &Header,
        output: Box<dyn std::io::Write + Send>,
        encoding: QualityEncoding,
        read_structures: &[ReadStructure],
    ) -> Result<u64> {
        // Detect if all inputs are BGZF
        let all_bgzf = self.inputs.iter().all(|p| {
            detect_compression_format(p).map(|f| f == CompressionFormat::Bgzf).unwrap_or(false)
        });

        let num_threads = self.threading.threads.unwrap_or(1);
        let config =
            FastqPipelineConfig::new(num_threads, all_bgzf, self.compression.compression_level)
                .with_stats(self.scheduler_opts.collect_stats())
                .with_scheduler_strategy(self.scheduler_opts.strategy())
                .with_deadlock_timeout(self.scheduler_opts.deadlock_timeout_secs())
                .with_deadlock_recovery(self.scheduler_opts.deadlock_recover_enabled())
                .with_synchronized(true); // Extract always uses synchronized FASTQs

        info!(
            "Using unified 7-step pipeline: threads={}, bgzf={}, stats={}, scheduler={:?}",
            num_threads,
            all_bgzf,
            self.scheduler_opts.collect_stats(),
            self.scheduler_opts.strategy()
        );

        // For Gzip/Plain: open decompressed readers upfront
        // For BGZF: pass None, pipeline opens files directly
        let decomp_threads = self.threading.num_threads().max(1);
        let decompressed_readers: Option<Vec<Box<dyn BufRead + Send>>> = if all_bgzf {
            None
        } else {
            Some(
                self.inputs
                    .iter()
                    .map(|p| open_fastq_reader(p, decomp_threads))
                    .collect::<Result<Vec<_>>>()?,
            )
        };

        // Clone read_structures for the closure
        let read_structures = read_structures.to_vec();

        let records_written = run_fastq_pipeline(
            config,
            &self.inputs,
            decompressed_readers,
            header,
            output,
            // process_fn: FastqTemplate → Vec<RecordBuf>
            move |template: FastqTemplate| -> std::io::Result<Vec<RecordBuf>> {
                // Debug: Log template record count
                if template.records.len() != read_structures.len() {
                    log::warn!(
                        "Template has {} records but expected {} (read_structures.len())",
                        template.records.len(),
                        read_structures.len()
                    );
                }

                // Validate read names match (synchronized mode defers this from Group step)
                if template.records.len() >= 2 {
                    let base_name = strip_read_suffix_extract(&template.records[0].name);
                    for (i, record) in template.records.iter().enumerate().skip(1) {
                        let other_base = strip_read_suffix_extract(&record.name);
                        if base_name != other_base {
                            return Err(std::io::Error::new(
                                std::io::ErrorKind::InvalidData,
                                format!(
                                    "FASTQ files out of sync: R1 has '{}', R{} has '{}'",
                                    String::from_utf8_lossy(base_name),
                                    i + 1,
                                    String::from_utf8_lossy(other_base),
                                ),
                            ));
                        }
                    }
                }
                // Convert each FastqRecord to a FastqSet using its read structure
                let mut fastq_sets: Vec<FastqSet> = Vec::with_capacity(template.records.len());
                for (record, rs) in template.records.iter().zip(read_structures.iter()) {
                    let fastq_set = FastqSet::from_record_with_structure(
                        &record.name,
                        &record.sequence,
                        &record.quality,
                        rs,
                        &[], // No skip reasons
                    );
                    fastq_sets.push(fastq_set);
                }

                // Combine all FastqSets into one
                let combined = FastqSet::combine_readsets(fastq_sets);

                // Create SAM records
                let result = make_sam_records_static(&combined, encoding).map_err(|e| {
                    std::io::Error::new(std::io::ErrorKind::InvalidData, e.to_string())
                })?;

                // Debug: Check if we produce fewer records than template has
                if result.len() != template.records.len() {
                    log::warn!(
                        "Template with {} FASTQ records produced {} BAM records",
                        template.records.len(),
                        result.len()
                    );
                }

                Ok(result)
            },
            // serialize_fn: Vec<RecordBuf> → writes directly to buffer
            |records: Vec<RecordBuf>, header: &Header, buffer: &mut Vec<u8>| {
                serialize_bam_records_into(&records, header, buffer)
            },
        )?;

        info!("Wrote {records_written} records via 7-step pipeline");
        Ok(records_written)
    }
}

/// Static version of `make_sam_records` that doesn't require &self.
/// This is needed for the 7-step pipeline closure.
fn make_sam_records_static(
    read_set: &FastqSet,
    encoding: QualityEncoding,
) -> Result<Vec<RecordBuf>> {
    let templates: Vec<&FastqSegment> = read_set.template_segments().collect();

    let read_name = String::from_utf8_lossy(&read_set.header);
    ensure!(!templates.is_empty(), "No template segments found for read: {read_name}");

    // Extract various barcode types as BString
    let cell_barcode_bs = Extract::join_bytes_with_separator(
        read_set.cell_barcode_segments().map(|s| s.seq.as_slice()),
        b'-',
    );
    let sample_barcode_bs = Extract::join_bytes_with_separator(
        read_set.sample_barcode_segments().map(|s| s.seq.as_slice()),
        b'-',
    );
    let umi_bs = Extract::join_bytes_with_separator(
        read_set.molecular_barcode_segments().map(|s| s.seq.as_slice()),
        b'-',
    );

    // Use extract_read_name_and_umi to strip description after space (same as legacy path)
    let (read_name_bytes, _umi_from_name) =
        Extract::extract_read_name_and_umi(&read_set.header, false);

    // Default read group ID
    let rgid_bs: BString = "A".into();

    let mut records = Vec::new();

    for (index, template) in templates.iter().enumerate() {
        let mut rec = RecordBuf::default();

        // Set read name
        *rec.name_mut() = Some(read_name_bytes.clone().into());

        // Set bases and qualities - if empty, substitute with single N @ Q2
        if template.seq.is_empty() {
            *rec.sequence_mut() = "N".as_bytes().to_vec().into();
            *rec.quality_scores_mut() = QualityScores::from(vec![2u8]);
        } else {
            *rec.sequence_mut() = template.seq.clone().into();
            let numeric_quals = encoding.to_standard_numeric(&template.quals);
            *rec.quality_scores_mut() = QualityScores::from(numeric_quals);
        }

        // Set flags for unmapped reads
        let mut flags = Flags::UNMAPPED;
        if templates.len() == 2 {
            flags |= Flags::SEGMENTED | Flags::MATE_UNMAPPED;
            if index == 0 {
                flags |= Flags::FIRST_SEGMENT;
            } else {
                flags |= Flags::LAST_SEGMENT;
            }
        }
        *rec.flags_mut() = flags;

        // Set mapping quality to 0 for unmapped reads
        *rec.mapping_quality_mut() = Some(
            noodles::sam::alignment::record::MappingQuality::new(0)
                .expect("0 is a valid mapping quality"),
        );

        // Add tags using the mutable data accessor
        let data = rec.data_mut();

        // Read group
        data.insert(Tag::READ_GROUP, Value::String(rgid_bs.as_bstr()).try_into()?);

        // Sample barcode
        if !sample_barcode_bs.is_empty() {
            data.insert(
                Tag::SAMPLE_BARCODE_SEQUENCE,
                Value::String(sample_barcode_bs.as_bstr()).try_into()?,
            );
        }

        // UMI (using default RX tag)
        if !umi_bs.is_empty() {
            Extract::add_tag(data, b"RX", &umi_bs).unwrap();
        }

        // Cell barcode (using default CB tag)
        if !cell_barcode_bs.is_empty() {
            Extract::add_tag(data, b"CB", &cell_barcode_bs).unwrap();
        }

        records.push(rec);
    }

    Ok(records)
}

impl Command for Extract {
    fn execute(&self, command_line: &str) -> Result<()> {
        // Validate inputs
        self.validate()?;

        let timer = OperationTimer::new("Extracting UMIs");
        let read_structures = self.get_read_structures()?;

        // Detect quality encoding from first 400 records
        // Use a separate reader for sampling to avoid consuming records from the main reader
        let mut sample_quals = Vec::new();
        let temp_reader =
            FastqReader::with_capacity(open_fastq_reader(&self.inputs[0], 1)?, BUFFER_SIZE);
        for (i, result) in temp_reader.into_records().enumerate() {
            if i >= QUALITY_DETECTION_SAMPLE_SIZE {
                break;
            }
            if let Ok(rec) = result {
                sample_quals.push(rec.qual().to_vec());
            }
        }

        let encoding = QualityEncoding::detect(&sample_quals)?;

        // Create header with @PG record
        let header = self.create_header(command_line)?;

        let records_written = if self.threading.is_parallel() {
            // Unified 7-step pipeline mode (--threads N was specified)
            // Uses work-stealing for better CPU utilization with strict thread cap
            // The 7-step pipeline handles its own BGZF compression internally
            let output: Box<dyn std::io::Write + Send> =
                Box::new(std::io::BufWriter::new(File::create(&self.output)?));

            self.process_with_pipeline(&header, output, encoding, &read_structures)?
        } else {
            // Legacy mode: use noodles multi-threaded BAM writer
            // Determine thread count for parallel BGZF decompression
            let decomp_threads = self.threading.num_threads().max(1);

            // Open input FASTQs with automatic BGZF detection
            let fq_readers: Vec<Box<dyn BufRead + Send>> = self
                .inputs
                .iter()
                .map(|p| open_fastq_reader(p, decomp_threads))
                .collect::<Result<Vec<_>>>()?;

            let fq_sources: Vec<FastqReader<Box<dyn BufRead + Send>>> = fq_readers
                .into_iter()
                .map(|fq| FastqReader::with_capacity(fq, BUFFER_SIZE))
                .collect();

            // Create iterators
            let fq_iterators: Vec<ReadSetIterator> = fq_sources
                .into_iter()
                .zip(read_structures.iter())
                .map(|(source, rs)| ReadSetIterator::new(rs.clone(), source, Vec::new()))
                .collect();

            // Create output writer with multi-threaded BGZF
            let writer_threads = self.threading.num_threads();
            let mut writer = create_bam_writer(
                &self.output,
                &header,
                writer_threads,
                self.compression.compression_level,
            )?;

            // Single-threaded processing: original simple loop (iterators are borrowed)
            let mut fq_iterators = fq_iterators;
            let count =
                self.process_singlethreaded(&mut fq_iterators, &header, &mut writer, encoding)?;

            // Flush and finish the writer
            writer.into_inner().finish()?;

            count
        };

        timer.log_completion(records_written);
        Ok(())
    }
}

/// Strip /1, /2, or other read pair suffixes from a read name.
///
/// This is used to validate that reads at the same position in synchronized FASTQs
/// have matching names. Returns a slice of the name without the suffix.
fn strip_read_suffix_extract(name: &[u8]) -> &[u8] {
    // First strip any space-separated comment suffix (e.g., "read/1 extra" -> "read/1")
    let name = if let Some(space_pos) = name.iter().position(|&b| b == b' ') {
        &name[..space_pos]
    } else {
        name
    };

    // Then strip common pair suffixes: /1, /2, .1, .2, _1, _2, :1, :2
    if name.len() >= 2 {
        let last = name[name.len() - 1];
        let sep = name[name.len() - 2];
        if (last == b'1' || last == b'2')
            && (sep == b'/' || sep == b'.' || sep == b'_' || sep == b':')
        {
            return &name[..name.len() - 2];
        }
    }

    name
}

#[cfg(test)]
mod tests {
    use super::*;
    use bstr::BString;
    use noodles::sam::alignment::record_buf::data::field::Value as RecordBufValue;
    use rstest::rstest;
    use std::fs::File;
    use std::io::Write;
    use tempfile::TempDir;

    /// Create a FASTQ file for testing
    ///
    /// # Arguments
    /// * `dir` - Temporary directory to create the file in
    /// * `name` - Name of the FASTQ file
    /// * `records` - Array of tuples containing (name, sequence, quality) for each record
    ///
    /// # Returns
    /// Path to the created FASTQ file
    fn create_fastq(dir: &TempDir, name: &str, records: &[(&str, &str, &str)]) -> PathBuf {
        let path = dir.path().join(name);
        let mut file = File::create(&path).unwrap();
        for (name, seq, qual) in records {
            writeln!(file, "@{name}").unwrap();
            writeln!(file, "{seq}").unwrap();
            writeln!(file, "+").unwrap();
            writeln!(file, "{qual}").unwrap();
        }
        path
    }

    /// Read all records from a BAM file into a vector
    ///
    /// # Arguments
    /// * `path` - Path to the BAM file
    ///
    /// # Returns
    /// Vector of all BAM records in the file
    fn read_bam_records(path: &PathBuf) -> Vec<RecordBuf> {
        let (mut reader, header) = create_bam_reader(path, 1).unwrap();

        let mut records = Vec::new();
        for result in reader.record_bufs(&header) {
            records.push(result.unwrap());
        }
        records
    }

    /// Extract a string tag value from a BAM record
    ///
    /// # Arguments
    /// * `record` - The BAM record
    /// * `tag_name` - Two-character tag name (e.g., "RX", "MI")
    ///
    /// # Returns
    /// The string value of the tag, or None if not present or not a string
    fn get_tag_string(record: &RecordBuf, tag_name: &str) -> Option<String> {
        let tag_bytes = tag_name.as_bytes();
        let tag = Tag::from([tag_bytes[0], tag_bytes[1]]);
        record.data().get(&tag).and_then(|value| match value {
            RecordBufValue::String(s) => Some(String::from_utf8_lossy(s.as_ref()).to_string()),
            _ => None,
        })
    }

    #[test]
    fn test_single_fastq_no_read_structure() {
        let tmp = TempDir::new().unwrap();
        let r1 = create_fastq(
            &tmp,
            "r1.fq",
            &[("q1", "AAAAAAAAAA", "=========="), ("q2", "CCCCCCCCCC", "##########")],
        );
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1],
            output: output.clone(),
            read_structures: vec![],
            umi_tag: "RX".to_string(),
            umi_qual_tag: None,
            cell_tag: "CB".to_string(),
            cell_qual_tag: None,
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names: false,
            annotate_read_names: false,
            single_tag: None,
            clipping_attribute: None,
            read_group_id: "A".to_string(),
            sample: "foo".to_string(),
            library: "bar".to_string(),
            barcode: None,
            platform: "illumina".to_string(),
            platform_unit: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            comment: vec![],
            run_date: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        extract.execute("test").unwrap();

        let recs = read_bam_records(&output);
        assert_eq!(recs.len(), 2);

        let q1 = &recs[0];
        assert_eq!(q1.name().map(|n| n.as_bytes()), Some(b"q1".as_slice()));
        assert!(!q1.flags().is_segmented());
        assert_eq!(q1.sequence().as_ref(), b"AAAAAAAAAA");
        assert_eq!(q1.quality_scores().as_ref(), &[28; 10]);

        let q2 = &recs[1];
        assert_eq!(q2.name().map(|n| n.as_bytes()), Some(b"q2".as_slice()));
        assert!(!q2.flags().is_segmented());
        assert_eq!(q2.sequence().as_ref(), b"CCCCCCCCCC");
        assert_eq!(q2.quality_scores().as_ref(), &[2; 10]);
    }

    #[test]
    fn test_paired_end_no_read_structures() {
        let tmp = TempDir::new().unwrap();
        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "AAAAAAAAAA", "==========")]);
        let r2 = create_fastq(&tmp, "r2.fq", &[("q1", "CCCCCCCCCC", "##########")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1, r2],
            output: output.clone(),
            read_structures: vec![],
            umi_tag: "RX".to_string(),
            umi_qual_tag: None,
            cell_tag: "CB".to_string(),
            cell_qual_tag: None,
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names: false,
            annotate_read_names: false,
            single_tag: None,
            clipping_attribute: None,
            read_group_id: "A".to_string(),
            sample: "pip".to_string(),
            library: "pop".to_string(),
            barcode: None,
            platform: "illumina".to_string(),
            platform_unit: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            comment: vec![],
            run_date: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        extract.execute("test").unwrap();

        let recs = read_bam_records(&output);
        assert_eq!(recs.len(), 2);

        let r1 = &recs[0];
        assert_eq!(r1.name().map(|n| n.as_bytes()), Some(b"q1".as_slice()));
        assert!(r1.flags().is_segmented());
        assert!(r1.flags().is_first_segment());
        assert!(!r1.flags().is_last_segment());
        assert_eq!(r1.sequence().as_ref(), b"AAAAAAAAAA");
        assert_eq!(r1.quality_scores().as_ref(), &[28; 10]);

        let r2 = &recs[1];
        assert_eq!(r2.name().map(|n| n.as_bytes()), Some(b"q1".as_slice()));
        assert!(r2.flags().is_segmented());
        assert!(!r2.flags().is_first_segment());
        assert!(r2.flags().is_last_segment());
        assert_eq!(r2.sequence().as_ref(), b"CCCCCCCCCC");
        assert_eq!(r2.quality_scores().as_ref(), &[2; 10]);
    }

    #[test]
    fn test_paired_end_ignore_umi_qual_tag_when_no_umis() {
        let tmp = TempDir::new().unwrap();
        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "AAAAAAAAAA", "==========")]);
        let r2 = create_fastq(&tmp, "r2.fq", &[("q1", "CCCCCCCCCC", "##########")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1, r2],
            output: output.clone(),
            read_structures: vec![],
            umi_tag: "RX".to_string(),
            umi_qual_tag: Some("QX".to_string()),
            cell_tag: "CB".to_string(),
            cell_qual_tag: None,
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names: false,
            annotate_read_names: false,
            single_tag: None,
            clipping_attribute: None,
            read_group_id: "A".to_string(),
            sample: "pip".to_string(),
            library: "pop".to_string(),
            barcode: None,
            platform: "illumina".to_string(),
            platform_unit: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            comment: vec![],
            run_date: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        extract.execute("test").unwrap();

        let recs = read_bam_records(&output);
        assert_eq!(recs.len(), 2);
        // Should not have RX or QX tags since no UMIs
        assert!(get_tag_string(&recs[0], "RX").is_none());
        assert!(get_tag_string(&recs[0], "QX").is_none());
    }

    #[test]
    fn test_paired_end_with_inline_umi() {
        let tmp = TempDir::new().unwrap();
        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "ACGTAAAAAA", "IIII======")]);
        let r2 = create_fastq(&tmp, "r2.fq", &[("q1", "CCCCCCCCCC", "##########")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1, r2],
            output: output.clone(),
            read_structures: vec![
                ReadStructure::from_str("4M+T").unwrap(),
                ReadStructure::from_str("+T").unwrap(),
            ],
            umi_tag: "RX".to_string(),
            umi_qual_tag: None,
            cell_tag: "CB".to_string(),
            cell_qual_tag: None,
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names: false,
            annotate_read_names: false,
            single_tag: None,
            clipping_attribute: None,
            read_group_id: "A".to_string(),
            sample: "s".to_string(),
            library: "l".to_string(),
            barcode: None,
            platform: "illumina".to_string(),
            platform_unit: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            comment: vec![],
            run_date: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        extract.execute("test").unwrap();

        let recs = read_bam_records(&output);
        assert_eq!(recs.len(), 2);

        let r1 = &recs[0];
        assert_eq!(r1.name().map(|n| n.as_bytes()), Some(b"q1".as_slice()));
        assert_eq!(r1.sequence().as_ref(), b"AAAAAA");
        assert_eq!(r1.quality_scores().as_ref(), &[28; 6]);
        assert_eq!(get_tag_string(r1, "RX"), Some("ACGT".to_string()));

        let r2 = &recs[1];
        assert_eq!(r2.name().map(|n| n.as_bytes()), Some(b"q1".as_slice()));
        assert_eq!(r2.sequence().as_ref(), b"CCCCCCCCCC");
        assert_eq!(r2.quality_scores().as_ref(), &[2; 10]);
        assert_eq!(get_tag_string(r2, "RX"), Some("ACGT".to_string()));
    }

    #[test]
    fn test_paired_end_with_inline_umi_keep_qualities() {
        let tmp = TempDir::new().unwrap();
        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "ACGTAAAAAA", "IIII======")]);
        let r2 = create_fastq(&tmp, "r2.fq", &[("q1", "CCCCCCCCCC", "##########")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1, r2],
            output: output.clone(),
            read_structures: vec![
                ReadStructure::from_str("4M+T").unwrap(),
                ReadStructure::from_str("+T").unwrap(),
            ],
            umi_tag: "RX".to_string(),
            umi_qual_tag: Some("QX".to_string()),
            cell_tag: "CB".to_string(),
            cell_qual_tag: None,
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names: false,
            annotate_read_names: false,
            single_tag: None,
            clipping_attribute: None,
            read_group_id: "A".to_string(),
            sample: "s".to_string(),
            library: "l".to_string(),
            barcode: None,
            platform: "illumina".to_string(),
            platform_unit: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            comment: vec![],
            run_date: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        extract.execute("test").unwrap();

        let recs = read_bam_records(&output);
        assert_eq!(recs.len(), 2);

        assert_eq!(get_tag_string(&recs[0], "RX"), Some("ACGT".to_string()));
        assert_eq!(get_tag_string(&recs[0], "QX"), Some("IIII".to_string()));
        assert_eq!(get_tag_string(&recs[1], "RX"), Some("ACGT".to_string()));
        assert_eq!(get_tag_string(&recs[1], "QX"), Some("IIII".to_string()));
    }

    #[test]
    fn test_complex_read_structures_multiple_segments() {
        let tmp = TempDir::new().unwrap();
        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "AAACCCTTTAAAAA", "==============")]);
        let r2 = create_fastq(&tmp, "r2.fq", &[("q1", "GGGTTTAAACCCCC", "##############")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1, r2],
            output: output.clone(),
            read_structures: vec![
                ReadStructure::from_str("3B3M3B5T").unwrap(),
                ReadStructure::from_str("3B3M3B5T").unwrap(),
            ],
            umi_tag: "RX".to_string(),
            umi_qual_tag: None,
            cell_tag: "CB".to_string(),
            cell_qual_tag: None,
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names: false,
            annotate_read_names: false,
            single_tag: None,
            clipping_attribute: None,
            read_group_id: "A".to_string(),
            sample: "s".to_string(),
            library: "l".to_string(),
            barcode: None,
            platform: "illumina".to_string(),
            platform_unit: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            comment: vec![],
            run_date: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        extract.execute("test").unwrap();

        let recs = read_bam_records(&output);
        assert_eq!(recs.len(), 2);

        let r1 = &recs[0];
        assert_eq!(r1.sequence().as_ref(), b"AAAAA");
        assert_eq!(r1.quality_scores().as_ref(), &[28; 5]);
        assert_eq!(get_tag_string(r1, "RX"), Some("CCC-TTT".to_string()));
        assert_eq!(get_tag_string(r1, "BC"), Some("AAA-TTT-GGG-AAA".to_string()));

        let r2 = &recs[1];
        assert_eq!(r2.sequence().as_ref(), b"CCCCC");
        assert_eq!(r2.quality_scores().as_ref(), &[2; 5]);
        assert_eq!(get_tag_string(r2, "RX"), Some("CCC-TTT".to_string()));
        assert_eq!(get_tag_string(r2, "BC"), Some("AAA-TTT-GGG-AAA".to_string()));
    }

    #[test]
    fn test_complex_read_structures_with_umi_qualities() {
        let tmp = TempDir::new().unwrap();
        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "AAACCCTTTAAAAA", "===III========")]);
        let r2 = create_fastq(&tmp, "r2.fq", &[("q1", "GGGTTTAAACCCCC", "###JJJ########")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1, r2],
            output: output.clone(),
            read_structures: vec![
                ReadStructure::from_str("3B3M3B5T").unwrap(),
                ReadStructure::from_str("3B3M3B5T").unwrap(),
            ],
            umi_tag: "RX".to_string(),
            umi_qual_tag: Some("QX".to_string()),
            cell_tag: "CB".to_string(),
            cell_qual_tag: None,
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names: false,
            annotate_read_names: false,
            single_tag: None,
            clipping_attribute: None,
            read_group_id: "A".to_string(),
            sample: "s".to_string(),
            library: "l".to_string(),
            barcode: None,
            platform: "illumina".to_string(),
            platform_unit: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            comment: vec![],
            run_date: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        extract.execute("test").unwrap();

        let recs = read_bam_records(&output);
        assert_eq!(recs.len(), 2);

        assert_eq!(get_tag_string(&recs[0], "RX"), Some("CCC-TTT".to_string()));
        assert_eq!(get_tag_string(&recs[0], "QX"), Some("III JJJ".to_string()));
        assert_eq!(get_tag_string(&recs[1], "RX"), Some("CCC-TTT".to_string()));
        assert_eq!(get_tag_string(&recs[1], "QX"), Some("III JJJ".to_string()));
    }

    #[test]
    fn test_four_fastqs() {
        let tmp = TempDir::new().unwrap();
        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "AAAAAAAAAA", "==========")]);
        let r2 = create_fastq(&tmp, "r2.fq", &[("q1", "CCCCCCCCCC", "##########")]);
        let r3 = create_fastq(&tmp, "r3.fq", &[("q1", "ACGT", "????")]);
        let r4 = create_fastq(&tmp, "r4.fq", &[("q1", "GAGA", "????")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1, r2, r3, r4],
            output: output.clone(),
            read_structures: vec![
                ReadStructure::from_str("+T").unwrap(),
                ReadStructure::from_str("+T").unwrap(),
                ReadStructure::from_str("4B").unwrap(),
                ReadStructure::from_str("4M").unwrap(),
            ],
            umi_tag: "RX".to_string(),
            umi_qual_tag: None,
            cell_tag: "CB".to_string(),
            cell_qual_tag: None,
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names: false,
            annotate_read_names: false,
            single_tag: None,
            clipping_attribute: None,
            read_group_id: "A".to_string(),
            sample: "s".to_string(),
            library: "l".to_string(),
            barcode: None,
            platform: "illumina".to_string(),
            platform_unit: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            comment: vec![],
            run_date: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        extract.execute("test").unwrap();

        let recs = read_bam_records(&output);
        assert_eq!(recs.len(), 2);

        let r1 = &recs[0];
        assert_eq!(r1.sequence().as_ref(), b"AAAAAAAAAA");
        assert_eq!(r1.quality_scores().as_ref(), &[28; 10]);
        assert_eq!(get_tag_string(r1, "RX"), Some("GAGA".to_string()));
        assert_eq!(get_tag_string(r1, "BC"), Some("ACGT".to_string()));

        let r2 = &recs[1];
        assert_eq!(r2.sequence().as_ref(), b"CCCCCCCCCC");
        assert_eq!(r2.quality_scores().as_ref(), &[2; 10]);
        assert_eq!(get_tag_string(r2, "RX"), Some("GAGA".to_string()));
        assert_eq!(get_tag_string(r2, "BC"), Some("ACGT".to_string()));
    }

    #[test]
    fn test_four_fastqs_with_umi_qualities() {
        let tmp = TempDir::new().unwrap();
        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "AAAAAAAAAA", "==========")]);
        let r2 = create_fastq(&tmp, "r2.fq", &[("q1", "CCCCCCCCCC", "##########")]);
        let r3 = create_fastq(&tmp, "r3.fq", &[("q1", "ACGT", "????")]);
        let r4 = create_fastq(&tmp, "r4.fq", &[("q1", "GAGA", "????")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1, r2, r3, r4],
            output: output.clone(),
            read_structures: vec![
                ReadStructure::from_str("+T").unwrap(),
                ReadStructure::from_str("+T").unwrap(),
                ReadStructure::from_str("4B").unwrap(),
                ReadStructure::from_str("4M").unwrap(),
            ],
            umi_tag: "RX".to_string(),
            umi_qual_tag: Some("QX".to_string()),
            cell_tag: "CB".to_string(),
            cell_qual_tag: None,
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names: false,
            annotate_read_names: false,
            single_tag: None,
            clipping_attribute: None,
            read_group_id: "A".to_string(),
            sample: "s".to_string(),
            library: "l".to_string(),
            barcode: None,
            platform: "illumina".to_string(),
            platform_unit: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            comment: vec![],
            run_date: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        extract.execute("test").unwrap();

        let recs = read_bam_records(&output);
        assert_eq!(recs.len(), 2);

        assert_eq!(get_tag_string(&recs[0], "RX"), Some("GAGA".to_string()));
        assert_eq!(get_tag_string(&recs[0], "QX"), Some("????".to_string()));
        assert_eq!(get_tag_string(&recs[1], "RX"), Some("GAGA".to_string()));
        assert_eq!(get_tag_string(&recs[1], "QX"), Some("????".to_string()));
    }

    #[test]
    fn test_header_metadata() {
        let tmp = TempDir::new().unwrap();
        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "AAAAAAAAAA", "==========")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1],
            output: output.clone(),
            read_structures: vec![ReadStructure::from_str("10T").unwrap()],
            umi_tag: "RX".to_string(),
            umi_qual_tag: None,
            cell_tag: "CB".to_string(),
            cell_qual_tag: None,
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names: false,
            annotate_read_names: false,
            single_tag: None,
            clipping_attribute: None,
            read_group_id: "MyRG".to_string(),
            sample: "foo".to_string(),
            library: "bar".to_string(),
            barcode: Some("TATA-GAGA".to_string()),
            platform: "Illumina".to_string(),
            platform_unit: Some("pee-eww".to_string()),
            platform_model: Some("hiseq2500".to_string()),
            sequencing_center: Some("nowhere".to_string()),
            predicted_insert_size: Some(300),
            description: Some("Some reads!".to_string()),
            comment: vec!["hello world".to_string(), "comment two".to_string()],
            run_date: Some("2024-01-01".to_string()),
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        extract.execute("test").unwrap();

        let (_reader, header) = create_bam_reader(&output, 1).unwrap();

        // Check comments
        let comments: Vec<String> =
            header.comments().iter().map(std::string::ToString::to_string).collect();
        assert_eq!(comments, vec!["hello world", "comment two"]);

        // Check read group exists and has correct ID
        let rg_id = BString::from("MyRG");
        let rg = header.read_groups().get(&rg_id);
        assert!(rg.is_some(), "Read group MyRG should exist");

        // Verify the read group was written by checking the record has the RG tag
        let recs = read_bam_records(&output);
        assert_eq!(recs.len(), 1);
        assert_eq!(get_tag_string(&recs[0], "RG"), Some("MyRG".to_string()));
    }

    #[test]
    fn test_zero_length_reads() {
        // Test that zero-length reads are handled gracefully with variable-length read structures
        // (matching Scala/fgbio behavior where empty reads become "N" @ Q2)
        let tmp = TempDir::new().unwrap();
        let r1 = create_fastq(
            &tmp,
            "r1.fq",
            &[("q1", "AAAAAAAAAA", "=========="), ("q2", "", ""), ("q3", "GGGGGG", "IIIIII")],
        );
        let r2 = create_fastq(
            &tmp,
            "r2.fq",
            &[("q1", "TTTTTTTTTT", "~~~~~~~~~~"), ("q2", "", ""), ("q3", "", "")],
        );
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1, r2],
            output: output.clone(),
            read_structures: vec![
                ReadStructure::from_str("+T").unwrap(), // Variable length to allow zero-length
                ReadStructure::from_str("+T").unwrap(),
            ],
            umi_tag: "RX".to_string(),
            umi_qual_tag: None,
            cell_tag: "CB".to_string(),
            cell_qual_tag: None,
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names: false,
            annotate_read_names: false,
            single_tag: None,
            clipping_attribute: None,
            read_group_id: "A".to_string(),
            sample: "s".to_string(),
            library: "l".to_string(),
            barcode: None,
            platform: "illumina".to_string(),
            platform_unit: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            comment: vec![],
            run_date: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        extract.execute("test").unwrap();

        // Verify the output BAM was created and contains records
        let records = read_bam_records(&output);

        assert_eq!(records.len(), 6, "Should have 6 records (3 pairs)");

        // Check that q2 (the one with zero-length reads on both sides) has "N" @ Q2 for both reads
        let q2_r1 = &records[2]; // Third record (q2, R1)
        let q2_r2 = &records[3]; // Fourth record (q2, R2)

        assert_eq!(q2_r1.sequence().as_ref(), b"N", "Zero-length R1 should become 'N'");
        assert_eq!(q2_r1.quality_scores().as_ref(), &[2u8], "Zero-length R1 should have Q2");

        assert_eq!(q2_r2.sequence().as_ref(), b"N", "Zero-length R2 should become 'N'");
        assert_eq!(q2_r2.quality_scores().as_ref(), &[2u8], "Zero-length R2 should have Q2");

        // Check that q3 (zero-length on R2 only) has "N" @ Q2 for R2 but normal sequence for R1
        let q3_r1 = &records[4]; // Fifth record (q3, R1)
        let q3_r2 = &records[5]; // Sixth record (q3, R2)

        assert_eq!(q3_r1.sequence().as_ref(), b"GGGGGG", "Non-zero R1 should retain sequence");
        assert_eq!(q3_r2.sequence().as_ref(), b"N", "Zero-length R2 should become 'N'");
        assert_eq!(q3_r2.quality_scores().as_ref(), &[2u8], "Zero-length R2 should have Q2");
    }

    #[test]
    #[should_panic(expected = "had too few bases to demux")]
    fn test_zero_length_reads_with_fixed_structure_should_panic() {
        // Test that zero-length reads with fixed-length read structures still panic
        // (as they should - you can't have 0 bases when 10T are required)
        let tmp = TempDir::new().unwrap();
        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "", "")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1],
            output: output.clone(),
            read_structures: vec![
                ReadStructure::from_str("10T").unwrap(), // Fixed length - should panic
            ],
            umi_tag: "RX".to_string(),
            umi_qual_tag: None,
            cell_tag: "CB".to_string(),
            cell_qual_tag: None,
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names: false,
            annotate_read_names: false,
            single_tag: None,
            clipping_attribute: None,
            read_group_id: "A".to_string(),
            sample: "s".to_string(),
            library: "l".to_string(),
            barcode: None,
            platform: "illumina".to_string(),
            platform_unit: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            comment: vec![],
            run_date: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        extract.execute("test").unwrap();
    }

    #[test]
    fn test_extract_sample_barcode_qualities() {
        let tmp = TempDir::new().unwrap();
        let r1 = create_fastq(
            &tmp,
            "r1.fq",
            &[
                ("q1:2:3:4:5:6:7", "AAACCCAAAA", "ABCDEFGHIJ"),
                ("q2:2:3:4:5:6:7", "TAANNNAAAA", "BCDEFGHIJK"),
                ("q3:2:3:4:5:6:7", "GAACCCTCGA", "CDEFGHIJKL"),
            ],
        );
        let output = tmp.path().join("output.bam");

        // Test with store_sample_barcode_qualities = true
        let extract = Extract {
            inputs: vec![r1.clone()],
            output: output.clone(),
            read_structures: vec![ReadStructure::from_str("3B3M3B+T").unwrap()],
            umi_tag: "RX".to_string(),
            umi_qual_tag: None,
            cell_tag: "CB".to_string(),
            cell_qual_tag: None,
            store_sample_barcode_qualities: true,
            extract_umis_from_read_names: false,
            annotate_read_names: false,
            single_tag: None,
            clipping_attribute: None,
            read_group_id: "A".to_string(),
            sample: "s".to_string(),
            library: "l".to_string(),
            barcode: None,
            platform: "illumina".to_string(),
            platform_unit: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            comment: vec![],
            run_date: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        extract.execute("test").unwrap();

        let recs = read_bam_records(&output);
        assert_eq!(recs.len(), 3);
        assert_eq!(get_tag_string(&recs[0], "BC"), Some("AAA-AAA".to_string()));
        assert_eq!(get_tag_string(&recs[1], "BC"), Some("TAA-AAA".to_string()));
        assert_eq!(get_tag_string(&recs[2], "BC"), Some("GAA-TCG".to_string()));

        assert_eq!(get_tag_string(&recs[0], "QT"), Some("ABC GHI".to_string()));
        assert_eq!(get_tag_string(&recs[1], "QT"), Some("BCD HIJ".to_string()));
        assert_eq!(get_tag_string(&recs[2], "QT"), Some("CDE IJK".to_string()));

        // Test with store_sample_barcode_qualities = false
        let output2 = tmp.path().join("output2.bam");
        let extract2 = Extract {
            inputs: vec![r1],
            output: output2.clone(),
            read_structures: vec![ReadStructure::from_str("3B3M3B+T").unwrap()],
            umi_tag: "RX".to_string(),
            umi_qual_tag: None,
            cell_tag: "CB".to_string(),
            cell_qual_tag: None,
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names: false,
            annotate_read_names: false,
            single_tag: None,
            clipping_attribute: None,
            read_group_id: "A".to_string(),
            sample: "s".to_string(),
            library: "l".to_string(),
            barcode: None,
            platform: "illumina".to_string(),
            platform_unit: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            comment: vec![],
            run_date: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        extract2.execute("test").unwrap();

        let recs2 = read_bam_records(&output2);
        assert_eq!(recs2.len(), 3);
        assert!(get_tag_string(&recs2[0], "QT").is_none());
        assert!(get_tag_string(&recs2[1], "QT").is_none());
        assert!(get_tag_string(&recs2[2], "QT").is_none());
    }

    #[test]
    fn test_extract_umis_from_read_names() {
        let tmp = TempDir::new().unwrap();
        let r1 = create_fastq(
            &tmp,
            "r1.fq",
            &[
                ("q1:2:3:4:5:6:7:ACGT", "AAAAAAAAAA", "=========="),
                ("q2:2:3:4:5:6:7:TTGA", "TAAAAAAAAA", "=========="),
                ("q3:2:3:4:5:6:7", "TAAAAAAAAA", "=========="),
            ],
        );
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1],
            output: output.clone(),
            read_structures: vec![],
            umi_tag: "RX".to_string(),
            umi_qual_tag: None,
            cell_tag: "CB".to_string(),
            cell_qual_tag: None,
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names: true,
            annotate_read_names: false,
            single_tag: None,
            clipping_attribute: None,
            read_group_id: "A".to_string(),
            sample: "s".to_string(),
            library: "l".to_string(),
            barcode: None,
            platform: "illumina".to_string(),
            platform_unit: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            comment: vec![],
            run_date: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        extract.execute("test").unwrap();

        let recs = read_bam_records(&output);
        assert_eq!(recs.len(), 3);
        assert_eq!(get_tag_string(&recs[0], "RX"), Some("ACGT".to_string()));
        assert_eq!(get_tag_string(&recs[1], "RX"), Some("TTGA".to_string()));
        assert!(get_tag_string(&recs[2], "RX").is_none());
    }

    #[test]
    fn test_extract_umis_from_read_names_and_sequences() {
        let tmp = TempDir::new().unwrap();
        let r1 = create_fastq(
            &tmp,
            "r1.fq",
            &[
                ("q1:2:3:4:5:6:7:ACGT+CGTA", "GGNCCGAAAAAAA", "============="),
                ("q2:2:3:4:5:6:7:TTGA+TAAT", "TANAACAAAAAAA", "============="),
            ],
        );
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1],
            output: output.clone(),
            read_structures: vec![ReadStructure::from_str("2M1S2M+T").unwrap()],
            umi_tag: "RX".to_string(),
            umi_qual_tag: None,
            cell_tag: "CB".to_string(),
            cell_qual_tag: None,
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names: true,
            annotate_read_names: false,
            single_tag: None,
            clipping_attribute: None,
            read_group_id: "A".to_string(),
            sample: "s".to_string(),
            library: "l".to_string(),
            barcode: None,
            platform: "illumina".to_string(),
            platform_unit: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            comment: vec![],
            run_date: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        extract.execute("test").unwrap();

        let recs = read_bam_records(&output);
        assert_eq!(recs.len(), 2);
        assert_eq!(get_tag_string(&recs[0], "RX"), Some("ACGT-CGTA-GG-CC".to_string()));
        assert_eq!(get_tag_string(&recs[1], "RX"), Some("TTGA-TAAT-TA-AA".to_string()));
    }

    #[test]
    fn test_extract_cell_barcodes() {
        let tmp = TempDir::new().unwrap();
        let r1 = create_fastq(
            &tmp,
            "r1.fq",
            &[
                ("q1:2:3:4:5:6:7:ACGT+CGTA", "GGNCCGAAAAAAA", "============="),
                ("q2:2:3:4:5:6:7:TTGA+TAAT", "TANAACAAAAAAA", "============="),
            ],
        );
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1],
            output: output.clone(),
            read_structures: vec![ReadStructure::from_str("2C1S2C+T").unwrap()],
            umi_tag: "RX".to_string(),
            umi_qual_tag: None,
            cell_tag: "CB".to_string(),
            cell_qual_tag: Some("CY".to_string()),
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names: false,
            annotate_read_names: false,
            single_tag: None,
            clipping_attribute: None,
            read_group_id: "A".to_string(),
            sample: "s".to_string(),
            library: "l".to_string(),
            barcode: None,
            platform: "illumina".to_string(),
            platform_unit: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            comment: vec![],
            run_date: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        extract.execute("test").unwrap();

        let recs = read_bam_records(&output);
        assert_eq!(recs.len(), 2);
        assert_eq!(get_tag_string(&recs[0], "CB"), Some("GG-CC".to_string()));
        assert_eq!(get_tag_string(&recs[0], "CY"), Some("== ==".to_string()));
        assert_eq!(get_tag_string(&recs[1], "CB"), Some("TA-AA".to_string()));
        assert_eq!(get_tag_string(&recs[1], "CY"), Some("== ==".to_string()));
    }

    #[test]
    #[should_panic(expected = "Read names do not match")]
    fn test_fail_mismatched_read_names() {
        let tmp = TempDir::new().unwrap();
        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "AAAAAAAAAA", "==========")]);
        let r2 = create_fastq(&tmp, "r2.fq", &[("x1", "CCCCCCCCCC", "##########")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1, r2],
            output,
            read_structures: vec![],
            umi_tag: "RX".to_string(),
            umi_qual_tag: None,
            cell_tag: "CB".to_string(),
            cell_qual_tag: None,
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names: false,
            annotate_read_names: false,
            single_tag: None,
            clipping_attribute: None,
            read_group_id: "A".to_string(),
            sample: "s".to_string(),
            library: "l".to_string(),
            barcode: None,
            platform: "illumina".to_string(),
            platform_unit: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            comment: vec![],
            run_date: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        extract.execute("test").unwrap();
    }

    #[test]
    #[should_panic(expected = "out of sync")]
    fn test_fail_mismatched_read_counts() {
        let tmp = TempDir::new().unwrap();
        let r1 = create_fastq(
            &tmp,
            "r1.fq",
            &[("q1", "AAAAAAAAAA", "=========="), ("q2", "TTTTTTTTTT", "??????????")],
        );
        let r2 = create_fastq(&tmp, "r2.fq", &[("q1", "CCCCCCCCCC", "##########")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1, r2],
            output,
            read_structures: vec![],
            umi_tag: "RX".to_string(),
            umi_qual_tag: None,
            cell_tag: "CB".to_string(),
            cell_qual_tag: None,
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names: false,
            annotate_read_names: false,
            single_tag: None,
            clipping_attribute: None,
            read_group_id: "A".to_string(),
            sample: "s".to_string(),
            library: "l".to_string(),
            barcode: None,
            platform: "illumina".to_string(),
            platform_unit: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            comment: vec![],
            run_date: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        extract.execute("test").unwrap();
    }

    #[test]
    #[should_panic(expected = "must be supplied the same number of times")]
    fn test_fail_mismatched_inputs_and_read_structures() {
        let tmp = TempDir::new().unwrap();
        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "AAAAAAAAAA", "==========")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1],
            output,
            read_structures: vec![
                ReadStructure::from_str("+T").unwrap(),
                ReadStructure::from_str("+T").unwrap(),
            ],
            umi_tag: "RX".to_string(),
            umi_qual_tag: None,
            cell_tag: "CB".to_string(),
            cell_qual_tag: None,
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names: false,
            annotate_read_names: false,
            single_tag: None,
            clipping_attribute: None,
            read_group_id: "A".to_string(),
            sample: "s".to_string(),
            library: "l".to_string(),
            barcode: None,
            platform: "illumina".to_string(),
            platform_unit: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            comment: vec![],
            run_date: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        extract.execute("test").unwrap();
    }

    #[test]
    fn test_annotate_read_names_with_umis() {
        let tmp = TempDir::new().unwrap();
        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "ACGTAAAAAA", "IIII======")]);
        let r2 = create_fastq(&tmp, "r2.fq", &[("q1", "CCCCCCCCCC", "##########")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1, r2],
            output: output.clone(),
            read_structures: vec![
                ReadStructure::from_str("4M+T").unwrap(),
                ReadStructure::from_str("+T").unwrap(),
            ],
            umi_tag: "RX".to_string(),
            umi_qual_tag: None,
            cell_tag: "CB".to_string(),
            cell_qual_tag: None,
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names: false,
            annotate_read_names: true, // Enable read name annotation
            single_tag: None,
            clipping_attribute: None,
            read_group_id: "A".to_string(),
            sample: "s".to_string(),
            library: "l".to_string(),
            barcode: None,
            platform: "illumina".to_string(),
            platform_unit: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            comment: vec![],
            run_date: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        extract.execute("test").unwrap();

        let recs = read_bam_records(&output);
        assert_eq!(recs.len(), 2);

        // Check that read names have "+ACGT" appended
        let r1 = &recs[0];
        assert_eq!(r1.name().map(|n| n.as_bytes()), Some(b"q1+ACGT".as_slice()));
        assert_eq!(get_tag_string(r1, "RX"), Some("ACGT".to_string()));

        let r2 = &recs[1];
        assert_eq!(r2.name().map(|n| n.as_bytes()), Some(b"q1+ACGT".as_slice()));
        assert_eq!(get_tag_string(r2, "RX"), Some("ACGT".to_string()));
    }

    #[test]
    fn test_annotate_read_names_no_umis() {
        let tmp = TempDir::new().unwrap();
        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "AAAAAAAAAA", "==========")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1],
            output: output.clone(),
            read_structures: vec![],
            umi_tag: "RX".to_string(),
            umi_qual_tag: None,
            cell_tag: "CB".to_string(),
            cell_qual_tag: None,
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names: false,
            annotate_read_names: true, // Enable read name annotation
            single_tag: None,
            clipping_attribute: None,
            read_group_id: "A".to_string(),
            sample: "s".to_string(),
            library: "l".to_string(),
            barcode: None,
            platform: "illumina".to_string(),
            platform_unit: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            comment: vec![],
            run_date: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        extract.execute("test").unwrap();

        let recs = read_bam_records(&output);
        assert_eq!(recs.len(), 1);

        // Read name should NOT have "+UMI" appended since there are no UMIs
        let r1 = &recs[0];
        assert_eq!(r1.name().map(|n| n.as_bytes()), Some(b"q1".as_slice()));
    }

    #[test]
    fn test_single_tag() {
        let tmp = TempDir::new().unwrap();
        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "AAACCCTTTAAAAA", "==============")]);
        let r2 = create_fastq(&tmp, "r2.fq", &[("q1", "GGGTTTAAACCCCC", "##############")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1, r2],
            output: output.clone(),
            read_structures: vec![
                ReadStructure::from_str("3B3M3B5T").unwrap(),
                ReadStructure::from_str("3B3M3B5T").unwrap(),
            ],
            umi_tag: "RX".to_string(),
            umi_qual_tag: None,
            cell_tag: "CB".to_string(),
            cell_qual_tag: None,
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names: false,
            annotate_read_names: false,
            single_tag: Some("ZU".to_string()), // Use single tag
            clipping_attribute: None,
            read_group_id: "A".to_string(),
            sample: "s".to_string(),
            library: "l".to_string(),
            barcode: None,
            platform: "illumina".to_string(),
            platform_unit: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            comment: vec![],
            run_date: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        extract.execute("test").unwrap();

        let recs = read_bam_records(&output);
        assert_eq!(recs.len(), 2);

        // Check that both RX and ZU tags have the same UMI
        let r1 = &recs[0];
        assert_eq!(get_tag_string(r1, "RX"), Some("CCC-TTT".to_string()));
        assert_eq!(get_tag_string(r1, "ZU"), Some("CCC-TTT".to_string()));

        let r2 = &recs[1];
        assert_eq!(get_tag_string(r2, "RX"), Some("CCC-TTT".to_string()));
        assert_eq!(get_tag_string(r2, "ZU"), Some("CCC-TTT".to_string()));
    }

    #[test]
    #[should_panic(expected = "Single tag must be exactly 2 characters")]
    fn test_single_tag_invalid_length() {
        let tmp = TempDir::new().unwrap();
        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "ACGTAAAAAA", "IIII======")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1],
            output,
            read_structures: vec![ReadStructure::from_str("4M+T").unwrap()],
            umi_tag: "RX".to_string(),
            umi_qual_tag: None,
            cell_tag: "CB".to_string(),
            cell_qual_tag: None,
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names: false,
            annotate_read_names: false,
            single_tag: Some("INVALID".to_string()), // Too long
            clipping_attribute: None,
            read_group_id: "A".to_string(),
            sample: "s".to_string(),
            library: "l".to_string(),
            barcode: None,
            platform: "illumina".to_string(),
            platform_unit: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            comment: vec![],
            run_date: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        extract.execute("test").unwrap();
    }

    #[test]
    #[should_panic(expected = "Single tag cannot be the same as umi-tag")]
    fn test_single_tag_same_as_umi_tag() {
        let tmp = TempDir::new().unwrap();
        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "ACGTAAAAAA", "IIII======")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1],
            output,
            read_structures: vec![ReadStructure::from_str("4M+T").unwrap()],
            umi_tag: "RX".to_string(),
            umi_qual_tag: None,
            cell_tag: "CB".to_string(),
            cell_qual_tag: None,
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names: false,
            annotate_read_names: false,
            single_tag: Some("RX".to_string()), // Same as umi_tag
            clipping_attribute: None,
            read_group_id: "A".to_string(),
            sample: "s".to_string(),
            library: "l".to_string(),
            barcode: None,
            platform: "illumina".to_string(),
            platform_unit: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            comment: vec![],
            run_date: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        extract.execute("test").unwrap();
    }

    #[test]
    fn test_combined_annotate_read_names_and_single_tag() {
        let tmp = TempDir::new().unwrap();
        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "ACGTAAAAAA", "IIII======")]);
        let r2 = create_fastq(&tmp, "r2.fq", &[("q1", "CCCCCCCCCC", "##########")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1, r2],
            output: output.clone(),
            read_structures: vec![
                ReadStructure::from_str("4M+T").unwrap(),
                ReadStructure::from_str("+T").unwrap(),
            ],
            umi_tag: "RX".to_string(),
            umi_qual_tag: None,
            cell_tag: "CB".to_string(),
            cell_qual_tag: None,
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names: false,
            annotate_read_names: true, // Both features enabled
            single_tag: Some("ZU".to_string()),
            clipping_attribute: None,
            read_group_id: "A".to_string(),
            sample: "s".to_string(),
            library: "l".to_string(),
            barcode: None,
            platform: "illumina".to_string(),
            platform_unit: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            comment: vec![],
            run_date: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        extract.execute("test").unwrap();

        let recs = read_bam_records(&output);
        assert_eq!(recs.len(), 2);

        // Check read name annotation
        let r1 = &recs[0];
        assert_eq!(r1.name().map(|n| n.as_bytes()), Some(b"q1+ACGT".as_slice()));

        // Check both tags have UMI
        assert_eq!(get_tag_string(r1, "RX"), Some("ACGT".to_string()));
        assert_eq!(get_tag_string(r1, "ZU"), Some("ACGT".to_string()));

        let r2 = &recs[1];
        assert_eq!(r2.name().map(|n| n.as_bytes()), Some(b"q1+ACGT".as_slice()));
        assert_eq!(get_tag_string(r2, "RX"), Some("ACGT".to_string()));
        assert_eq!(get_tag_string(r2, "ZU"), Some("ACGT".to_string()));
    }

    #[test]
    #[should_panic(expected = "had too few bases to demux")]
    fn test_fail_read_too_short_for_structure() {
        // Test that we fail when a read is not long enough for the read structure
        // Read is 10 bases, but structure requires 8M + 8S = 16 bases minimum
        let tmp = TempDir::new().unwrap();
        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "AAAAAAAAAA", "==========")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1],
            output,
            read_structures: vec![ReadStructure::from_str("8M8S+T").unwrap()],
            umi_tag: "RX".to_string(),
            umi_qual_tag: None,
            cell_tag: "CB".to_string(),
            cell_qual_tag: None,
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names: false,
            annotate_read_names: false,
            single_tag: None,
            clipping_attribute: None,
            read_group_id: "A".to_string(),
            sample: "s".to_string(),
            library: "l".to_string(),
            barcode: None,
            platform: "illumina".to_string(),
            platform_unit: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            comment: vec![],
            run_date: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        // Should panic because read is too short
        extract.execute("test").unwrap();
    }

    #[test]
    fn test_variable_length_reads_with_plus() {
        // Test that the '+' operator in read structures handles variable-length reads
        let tmp = TempDir::new().unwrap();
        let r1 = create_fastq(
            &tmp,
            "r1.fq",
            &[
                ("q1", "AAAATTTTCCCCGGGGAAAAA", "====================="), // 21 bases
                ("q2", "AAAATTTTCCCCGGGGAAAAATTTT", "========================="), // 25 bases
                ("q3", "AAAATTTTCCCCGGGG", "================"), // 16 bases (exact minimum)
            ],
        );
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1],
            output: output.clone(),
            read_structures: vec![ReadStructure::from_str("4M4B4S+T").unwrap()],
            umi_tag: "RX".to_string(),
            umi_qual_tag: None,
            cell_tag: "CB".to_string(),
            cell_qual_tag: None,
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names: false,
            annotate_read_names: false,
            single_tag: None,
            clipping_attribute: None,
            read_group_id: "A".to_string(),
            sample: "s".to_string(),
            library: "l".to_string(),
            barcode: None,
            platform: "illumina".to_string(),
            platform_unit: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            comment: vec![],
            run_date: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        extract.execute("test").unwrap();

        let recs = read_bam_records(&output);
        assert_eq!(recs.len(), 3);

        // Check the template sequences have correct lengths
        // 4M + 4B + 4S = 12 bases consumed, rest is template
        assert_eq!(recs[0].sequence().len(), 9); // 21 - 12 = 9
        assert_eq!(recs[1].sequence().len(), 13); // 25 - 12 = 13
        assert_eq!(recs[2].sequence().len(), 4); // 16 - 12 = 4

        // Check UMI extraction
        assert_eq!(get_tag_string(&recs[0], "RX"), Some("AAAA".to_string()));
        assert_eq!(get_tag_string(&recs[1], "RX"), Some("AAAA".to_string()));
        assert_eq!(get_tag_string(&recs[2], "RX"), Some("AAAA".to_string()));
    }

    // ========================================================================
    // Quality Encoding Detection Tests
    // ========================================================================

    #[test]
    fn test_quality_encoding_detection_phred33() {
        // Test detection of standard Phred+33 encoding (Sanger/Illumina 1.8+)
        // Quality scores with ASCII values 33-126
        let records = vec![
            vec![33u8, 40, 50, 60, 70, 80], // Mix of low, medium, high quality
            vec![35u8, 45, 55, 65, 75, 85], // Another mix
            vec![40u8, 50, 60, 70, 80, 90], // Medium to high quality
        ];

        let encoding = QualityEncoding::detect(&records).unwrap();
        assert_eq!(encoding, QualityEncoding::Standard);
    }

    #[test]
    fn test_quality_encoding_detection_phred64() {
        // Test detection of Phred+64 encoding (Illumina 1.3-1.7)
        // Quality scores with ASCII values 64-126
        let records = vec![
            vec![64u8, 70, 80, 90, 100],  // Illumina 1.3-1.7 range
            vec![65u8, 75, 85, 95, 105],  // Another mix in that range
            vec![70u8, 80, 90, 100, 110], // Higher quality scores
        ];

        let encoding = QualityEncoding::detect(&records).unwrap();
        assert_eq!(encoding, QualityEncoding::Illumina);
    }

    #[test]
    fn test_quality_encoding_detection_high_quality_phred33() {
        // Test that high-quality Phred+33 data (all scores >= 64) is still detected as Phred+33
        // when the range is narrow (not spanning typical Phred+64 range)
        let records = vec![
            vec![64u8, 65, 66, 67, 68], // High quality, narrow range
            vec![64u8, 65, 66, 67, 69], // High quality, narrow range
        ];

        let encoding = QualityEncoding::detect(&records).unwrap();
        // Should be Standard because range is too narrow for typical Phred+64
        assert_eq!(encoding, QualityEncoding::Standard);
    }

    #[test]
    fn test_quality_encoding_detection_empty_reads() {
        // Test that all empty reads default to Standard encoding
        let records = vec![vec![], vec![], vec![]];

        let encoding = QualityEncoding::detect(&records).unwrap();
        assert_eq!(encoding, QualityEncoding::Standard);
    }

    #[test]
    fn test_quality_encoding_detection_mixed_empty_and_valid() {
        // Test that empty reads are skipped and detection works on valid ones
        let records = vec![
            vec![],             // Empty
            vec![40u8, 50, 60], // Valid Phred+33
            vec![],             // Empty
            vec![45u8, 55, 65], // Valid Phred+33
        ];

        let encoding = QualityEncoding::detect(&records).unwrap();
        assert_eq!(encoding, QualityEncoding::Standard);
    }

    #[test]
    fn test_quality_encoding_detection_empty_input() {
        // Test that empty input produces an error
        let records: Vec<Vec<u8>> = vec![];

        let result = QualityEncoding::detect(&records);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("no records provided"));
    }

    #[test]
    fn test_quality_encoding_detection_invalid_low() {
        // Test that quality scores below 33 produce an error
        let records = vec![vec![20u8, 30, 40, 50]]; // 20 is below valid range

        let result = QualityEncoding::detect(&records);
        assert!(result.is_err());
        let err_msg = result.unwrap_err().to_string();
        assert!(err_msg.contains("Invalid quality scores"));
        assert!(err_msg.contains("20"));
    }

    #[test]
    fn test_quality_encoding_detection_invalid_high() {
        // Test that quality scores above 126 produce an error
        let records = vec![vec![50u8, 60, 70, 127]]; // 127 is above valid range

        let result = QualityEncoding::detect(&records);
        assert!(result.is_err());
        let err_msg = result.unwrap_err().to_string();
        assert!(err_msg.contains("Invalid quality scores"));
        assert!(err_msg.contains("127"));
    }

    #[test]
    fn test_quality_encoding_detection_ambiguous_range() {
        // Test the ambiguous range (59-63) defaults to Standard
        // This is Q26-Q30 in Phred+33, which is reasonable quality
        let records = vec![
            vec![59u8, 60, 61, 62, 63], // Right in the ambiguous zone
            vec![60u8, 61, 62],         // Also ambiguous
        ];

        let encoding = QualityEncoding::detect(&records).unwrap();
        // Should default to Standard (Phred+33) as it's more common
        assert_eq!(encoding, QualityEncoding::Standard);
    }

    #[test]
    fn test_quality_encoding_to_standard_numeric_phred33() {
        // Test conversion from Phred+33 ASCII to numeric quality scores
        let encoding = QualityEncoding::Standard;
        let quals = vec![33u8, 43, 53, 63, 73]; // ASCII values
        let numeric = encoding.to_standard_numeric(&quals);
        assert_eq!(numeric, vec![0u8, 10, 20, 30, 40]); // Q scores
    }

    #[test]
    fn test_quality_encoding_to_standard_numeric_phred64() {
        // Test conversion from Phred+64 ASCII to numeric quality scores
        let encoding = QualityEncoding::Illumina;
        let quals = vec![64u8, 74, 84, 94, 104]; // ASCII values
        let numeric = encoding.to_standard_numeric(&quals);
        assert_eq!(numeric, vec![0u8, 10, 20, 30, 40]); // Q scores
    }

    #[test]
    fn test_quality_encoding_to_standard_numeric_empty() {
        // Test conversion of empty quality string
        let encoding = QualityEncoding::Standard;
        let quals: Vec<u8> = vec![];
        let numeric = encoding.to_standard_numeric(&quals);
        assert_eq!(numeric, Vec::<u8>::new());
    }

    #[test]
    fn test_zero_length_reads_integration_with_quality_detection() {
        // Integration test: verify zero-length reads work correctly with quality detection
        let tmp = TempDir::new().unwrap();
        let r1 = create_fastq(
            &tmp,
            "r1.fq",
            &[
                ("q1", "ACGTACGT", "IIIIIIII"), // Normal read, high quality (Phred+33)
                ("q2", "", ""),                 // Zero-length read
                ("q3", "TTTTTTTT", "########"), // Normal read, low quality (Phred+33)
            ],
        );
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1],
            output: output.clone(),
            read_structures: vec![ReadStructure::from_str("+T").unwrap()],
            umi_tag: "RX".to_string(),
            umi_qual_tag: None,
            cell_tag: "CB".to_string(),
            cell_qual_tag: None,
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names: false,
            annotate_read_names: false,
            single_tag: None,
            clipping_attribute: None,
            read_group_id: "A".to_string(),
            sample: "s".to_string(),
            library: "l".to_string(),
            barcode: None,
            platform: "illumina".to_string(),
            platform_unit: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            comment: vec![],
            run_date: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        // Should succeed without panicking
        extract.execute("test").unwrap();

        let records = read_bam_records(&output);
        assert_eq!(records.len(), 3);

        // Verify first read is normal
        assert_eq!(records[0].sequence().as_ref(), b"ACGTACGT");

        // Verify second read (zero-length) became "N" @ Q2
        assert_eq!(records[1].sequence().as_ref(), b"N");
        assert_eq!(records[1].quality_scores().as_ref(), &[2u8]);

        // Verify third read is normal
        assert_eq!(records[2].sequence().as_ref(), b"TTTTTTTT");
    }

    #[test]
    fn test_phred64_fastq_end_to_end() {
        // End-to-end test with actual Phred+64 encoded FASTQ
        let tmp = TempDir::new().unwrap();

        // Create FASTQ with Phred+64 quality scores
        // ASCII 64 = Q0, ASCII 70 = Q6, ASCII 80 = Q16, etc.
        let r1 = create_fastq(
            &tmp,
            "r1.fq",
            &[
                ("q1", "ACGTACGT", "@@DDHHLL"), // Phred+64: Q0,Q0,Q4,Q4,Q8,Q8,Q12,Q12
                ("q2", "GGGGGGGG", "PPPPPPPP"), // Phred+64: Q16 for all
            ],
        );
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1],
            output: output.clone(),
            read_structures: vec![ReadStructure::from_str("+T").unwrap()],
            umi_tag: "RX".to_string(),
            umi_qual_tag: None,
            cell_tag: "CB".to_string(),
            cell_qual_tag: None,
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names: false,
            annotate_read_names: false,
            single_tag: None,
            clipping_attribute: None,
            read_group_id: "A".to_string(),
            sample: "s".to_string(),
            library: "l".to_string(),
            barcode: None,
            platform: "illumina".to_string(),
            platform_unit: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            comment: vec![],
            run_date: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        extract.execute("test").unwrap();

        let records = read_bam_records(&output);
        assert_eq!(records.len(), 2);

        // Verify quality scores were correctly converted from Phred+64 to numeric
        // @@DDHHLL (ASCII 64,64,68,68,72,72,76,76) -> Q0,Q0,Q4,Q4,Q8,Q8,Q12,Q12
        assert_eq!(records[0].quality_scores().as_ref(), &[0u8, 0, 4, 4, 8, 8, 12, 12]);

        // PPPPPPPP (ASCII 80 repeated) -> Q16 repeated
        assert_eq!(records[1].quality_scores().as_ref(), &[16u8; 8]);
    }

    #[test]
    fn test_paired_end_with_different_read_structures() {
        // Test paired-end with different read structures for R1 and R2
        let tmp = TempDir::new().unwrap();

        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "AAAATTTTCCCCGGGG", "IIIIIIIIIIIIIIII")]);
        let r2 = create_fastq(&tmp, "r2.fq", &[("q1", "TTTTGGGG", "IIIIIIII")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1, r2],
            output: output.clone(),
            read_structures: vec![
                ReadStructure::from_str("4M4S+T").unwrap(), // R1: 4M UMI, 4S skip, rest template
                ReadStructure::from_str("4M+T").unwrap(),   // R2: 4M UMI, rest template
            ],
            umi_tag: "RX".to_string(),
            umi_qual_tag: None,
            cell_tag: "CB".to_string(),
            cell_qual_tag: None,
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names: false,
            annotate_read_names: false,
            single_tag: None,
            clipping_attribute: None,
            read_group_id: "A".to_string(),
            sample: "s".to_string(),
            library: "l".to_string(),
            barcode: None,
            platform: "illumina".to_string(),
            platform_unit: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            comment: vec![],
            run_date: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        extract.execute("test").unwrap();

        let records = read_bam_records(&output);
        assert_eq!(records.len(), 2); // R1 and R2

        // R1: Template should be "CCCCGGGG" (after skipping 4M and 4S)
        assert_eq!(records[0].sequence().as_ref(), b"CCCCGGGG");

        // R2: Template should be "GGGG" (after extracting 4M UMI)
        assert_eq!(records[1].sequence().as_ref(), b"GGGG");

        // Both should have same UMI: AAAA-TTTT
        let rx_tag = noodles::sam::alignment::record::data::field::Tag::from([b'R', b'X']);
        let r1_umi = records[0].data().get(&rx_tag).unwrap();
        let r2_umi = records[1].data().get(&rx_tag).unwrap();

        if let noodles::sam::alignment::record_buf::data::field::Value::String(s) = r1_umi {
            assert_eq!(<bstr::BString as AsRef<[u8]>>::as_ref(s), b"AAAA-TTTT");
        }
        if let noodles::sam::alignment::record_buf::data::field::Value::String(s) = r2_umi {
            assert_eq!(<bstr::BString as AsRef<[u8]>>::as_ref(s), b"AAAA-TTTT");
        }
    }

    #[test]
    fn test_multithreaded_extraction() {
        // Test that multi-threaded extraction works correctly
        let tmp = TempDir::new().unwrap();

        let r1 = create_fastq(
            &tmp,
            "r1.fq",
            &[
                ("q1", "AAAAACGTACGT", "IIIIIIIIIIII"),
                ("q2", "TTTTTTTTTTTT", "IIIIIIIIIIII"),
                ("q3", "CCCCGGGGAAAA", "IIIIIIIIIIII"),
            ],
        );
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1],
            output: output.clone(),
            read_structures: vec![ReadStructure::from_str("5M+T").unwrap()],
            umi_tag: "RX".to_string(),
            umi_qual_tag: None,
            cell_tag: "CB".to_string(),
            cell_qual_tag: None,
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names: false,
            annotate_read_names: false,
            single_tag: None,
            clipping_attribute: None,
            read_group_id: "A".to_string(),
            sample: "s".to_string(),
            library: "l".to_string(),
            barcode: None,
            platform: "illumina".to_string(),
            platform_unit: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            comment: vec![],
            run_date: None,
            threading: ThreadingOptions::new(4), // Use multiple threads
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        extract.execute("test").unwrap();

        let records = read_bam_records(&output);
        assert_eq!(records.len(), 3);

        // Verify all records were processed correctly (5M UMI extracted, rest is template)
        assert_eq!(records[0].sequence().as_ref(), b"CGTACGT");
        assert_eq!(records[1].sequence().as_ref(), b"TTTTTTT");
        assert_eq!(records[2].sequence().as_ref(), b"GGGAAAA");
    }

    #[test]
    fn test_multithreaded_extraction_preserves_order() {
        // Test that multi-threaded extraction preserves input order
        // This is critical for downstream tools that expect ordered output
        let tmp = TempDir::new().unwrap();

        // Create 100 reads with sequential names for easy order verification
        let reads: Vec<(&str, &str, &str)> = (0..100)
            .map(|i| {
                // Using a leak to get 'static lifetime - OK for tests
                let name: &'static str = Box::leak(format!("read_{i:03}").into_boxed_str());
                let seq: &'static str =
                    Box::leak(format!("AAAAA{}", "ACGT".repeat(10)).into_boxed_str());
                let qual: &'static str = Box::leak("I".repeat(45).into_boxed_str());
                (name, seq, qual)
            })
            .collect();

        let r1 = create_fastq(&tmp, "r1.fq", &reads);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1],
            output: output.clone(),
            read_structures: vec![ReadStructure::from_str("5M+T").unwrap()],
            umi_tag: "RX".to_string(),
            umi_qual_tag: None,
            cell_tag: "CB".to_string(),
            cell_qual_tag: None,
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names: false,
            annotate_read_names: false,
            single_tag: None,
            clipping_attribute: None,
            read_group_id: "A".to_string(),
            sample: "s".to_string(),
            library: "l".to_string(),
            barcode: None,
            platform: "illumina".to_string(),
            platform_unit: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            comment: vec![],
            run_date: None,
            threading: ThreadingOptions::new(8), // Use multiple threads
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        extract.execute("test").unwrap();

        let records = read_bam_records(&output);
        assert_eq!(records.len(), 100);

        // Verify output order matches input order by checking read names
        for (i, record) in records.iter().enumerate() {
            let expected_name = format!("read_{i:03}");
            let actual_name = record.name().map(|n| n.as_bytes()).unwrap_or_default();
            assert_eq!(
                actual_name,
                expected_name.as_bytes(),
                "Output order mismatch at position {i}: expected {expected_name}, got {}",
                String::from_utf8_lossy(actual_name)
            );
        }
    }

    #[test]
    fn test_sample_barcode_with_quality_tags_specified() {
        // Test that extract works with barcode and quality tag parameters
        let tmp = TempDir::new().unwrap();

        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "AAAAACGTACGT", "IIIIIIIIIIII")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1],
            output: output.clone(),
            read_structures: vec![ReadStructure::from_str("5B+T").unwrap()], // 5B for barcode
            umi_tag: "RX".to_string(),
            umi_qual_tag: Some("QX".to_string()),
            cell_tag: "CB".to_string(),
            cell_qual_tag: Some("CY".to_string()),
            store_sample_barcode_qualities: true,
            extract_umis_from_read_names: false,
            annotate_read_names: false,
            single_tag: None,
            clipping_attribute: None,
            read_group_id: "A".to_string(),
            sample: "s".to_string(),
            library: "l".to_string(),
            barcode: Some("AAAAA".to_string()),
            platform: "illumina".to_string(),
            platform_unit: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            comment: vec![],
            run_date: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        // Should succeed with all quality tag parameters specified
        extract.execute("test").unwrap();

        let records = read_bam_records(&output);
        assert_eq!(records.len(), 1);

        // Verify template sequence is correct (5B barcode extracted, rest is template)
        assert_eq!(records[0].sequence().as_ref(), b"CGTACGT");
    }

    #[test]
    fn test_mapping_quality_is_zero_for_unmapped_reads() -> Result<()> {
        let dir = TempDir::new()?;

        // Create simple paired-end FASTQs
        let fastq1 = create_fastq(&dir, "r1.fq", &[("read1", "ACGTACGTAC", "IIIIIIIIII")]);
        let fastq2 = create_fastq(&dir, "r2.fq", &[("read1", "TGCATGCATG", "IIIIIIIIII")]);
        let output = dir.path().join("output.bam");

        let extract = Extract {
            inputs: vec![fastq1, fastq2],
            output: output.clone(),
            read_structures: vec![],
            umi_tag: "RX".to_string(),
            umi_qual_tag: None,
            cell_tag: "CB".to_string(),
            cell_qual_tag: None,
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names: false,
            annotate_read_names: false,
            single_tag: None,
            clipping_attribute: None,
            read_group_id: "A".to_string(),
            sample: "sample".to_string(),
            library: "library".to_string(),
            barcode: None,
            platform: "illumina".to_string(),
            platform_unit: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            comment: vec![],
            run_date: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        extract.execute("test")?;

        let records = read_bam_records(&output);
        assert_eq!(records.len(), 2, "Should have 2 records for paired-end");

        // Verify MAPQ is 0 for all unmapped reads
        for record in &records {
            let mapq = record.mapping_quality();
            assert!(mapq.is_some(), "Mapping quality should be set (not None/255)");
            let mapq_value: u8 = mapq.unwrap().into();
            assert_eq!(mapq_value, 0, "Mapping quality should be 0 for unmapped reads");
        }

        Ok(())
    }

    #[test]
    fn test_compression_format_detection_plain() {
        // Create a plain FASTQ file (not compressed)
        let tmp = TempDir::new().unwrap();
        let plain_path = tmp.path().join("test.fq");
        let mut file = File::create(&plain_path).unwrap();
        use std::io::Write;
        writeln!(file, "@read1").unwrap();
        writeln!(file, "ACGT").unwrap();
        writeln!(file, "+").unwrap();
        writeln!(file, "IIII").unwrap();

        let format = detect_compression_format(&plain_path).unwrap();
        assert_eq!(format, CompressionFormat::Plain);
    }

    #[test]
    fn test_compression_format_detection_gzip() {
        use flate2::Compression;
        use flate2::write::GzEncoder;
        use std::io::Write;

        // Create a gzip-compressed FASTQ file
        let tmp = TempDir::new().unwrap();
        let gz_path = tmp.path().join("test.fq.gz");

        // Use flate2 to create a proper gzip file
        let file = File::create(&gz_path).unwrap();
        let mut encoder = GzEncoder::new(file, Compression::default());
        writeln!(encoder, "@read1").unwrap();
        writeln!(encoder, "ACGT").unwrap();
        writeln!(encoder, "+").unwrap();
        writeln!(encoder, "IIII").unwrap();
        encoder.finish().unwrap();

        let format = detect_compression_format(&gz_path).unwrap();
        assert_eq!(format, CompressionFormat::Gzip);
    }

    #[test]
    fn test_compression_format_detection_bgzf() {
        use noodles_bgzf::io::Writer as BgzfWriter;
        use std::io::Write;

        // Create a BGZF-compressed file using noodles
        let tmp = TempDir::new().unwrap();
        let bgzf_path = tmp.path().join("test.fq.bgz");

        let file = File::create(&bgzf_path).unwrap();
        let mut writer = BgzfWriter::new(file);
        writeln!(writer, "@read1").unwrap();
        writeln!(writer, "ACGT").unwrap();
        writeln!(writer, "+").unwrap();
        writeln!(writer, "IIII").unwrap();
        writer.finish().unwrap();

        let format = detect_compression_format(&bgzf_path).unwrap();
        assert_eq!(format, CompressionFormat::Bgzf);
    }

    /// Test that extraction works correctly across all threading modes.
    /// This parameterized test ensures both the single-threaded fast path (None)
    /// and the multi-threaded pipeline (Some(1), Some(2)) produce correct results.
    #[rstest]
    #[case::fast_path(ThreadingOptions::none())]
    #[case::pipeline_1(ThreadingOptions::new(1))]
    #[case::pipeline_2(ThreadingOptions::new(2))]
    fn test_threading_modes(#[case] threading: ThreadingOptions) -> Result<()> {
        let tmp = TempDir::new()?;

        let r1 = create_fastq(
            &tmp,
            "r1.fq",
            &[
                ("q1", "AAAAACGTACGT", "IIIIIIIIIIII"),
                ("q2", "TTTTTTTTTTTT", "IIIIIIIIIIII"),
                ("q3", "CCCCGGGGAAAA", "IIIIIIIIIIII"),
            ],
        );
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1],
            output: output.clone(),
            read_structures: vec![ReadStructure::from_str("5M+T").unwrap()],
            umi_tag: "RX".to_string(),
            umi_qual_tag: None,
            cell_tag: "CB".to_string(),
            cell_qual_tag: None,
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names: false,
            annotate_read_names: false,
            single_tag: None,
            clipping_attribute: None,
            read_group_id: "A".to_string(),
            sample: "s".to_string(),
            library: "l".to_string(),
            barcode: None,
            platform: "illumina".to_string(),
            platform_unit: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            comment: vec![],
            run_date: None,
            threading,
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        extract.execute("test")?;

        let records = read_bam_records(&output);
        assert_eq!(records.len(), 3);

        // Verify all records were processed correctly (5M UMI extracted, rest is template)
        assert_eq!(records[0].sequence().as_ref(), b"CGTACGT");
        assert_eq!(records[1].sequence().as_ref(), b"TTTTTTT");
        assert_eq!(records[2].sequence().as_ref(), b"GGGAAAA");

        Ok(())
    }
}

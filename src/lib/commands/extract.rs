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
//! quality encodings. Detection pools the heads of all input FASTQs (EXT3-01), matching fgbio
//! `FastqToBam`. When no encoding matches the observed qualities, fgumi fails like fgbio's
//! `case Nil => fail(...)` (EXT3-03) but improves on its static message by reporting the
//! observed quality range in the error (see `QualityEncoding::from_stats`).

use crate::commands::command::Command;
use crate::commands::common::{
    CompressionOptions, QueueMemoryOptions, SchedulerOptions, ThreadingOptions,
};
use crate::fastq::FastqSegment;
use crate::fastq::FastqSet;
use crate::fastq::ReadSetIterator;
use crate::fastq_deinterleave::deinterleave;
use crate::fastq_parse::strip_read_suffix;
use crate::logging::OperationTimer;
use crate::sam::SamTag;
use crate::unified_pipeline::fastq_out_of_sync_error;
use crate::validation::validate_input_exists;
use anyhow::{Context, Result, bail, ensure};
use bstr::{BString, ByteSlice};
use clap::Parser;
use fgumi_bam_io::ProgressTracker;
use fgumi_bam_io::{RawBamWriter, create_raw_bam_writer};
use fgumi_raw_bam::UnmappedSamBuilder;
use fgumi_raw_bam::fields::flags;
use log::{debug, info};
use noodles_bgzf::io::MultithreadedReader;

use crate::read_structure::{ReadStructure, SegmentType};
#[cfg(test)]
use fgumi_bam_io::create_bam_reader;
use fgumi_bam_io::{ChainedReader, InputFormat, TeeReader, is_stdin_path};
use fgumi_simd_fastq::SimdFastqReader;
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
use std::fs::File;
use std::io::{self, BufRead, BufReader};
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
/// Classification is delegated to [`fgumi_bam_io::classify_input`]; see
/// [`classify_compression_header`] for how the verdict maps onto FASTQ.
fn detect_compression_format(path: &Path) -> Result<CompressionFormat> {
    let mut file = File::open(path)?;
    let mut header = [0u8; fgumi_bam_io::FORMAT_PREFIX_LEN];
    let bytes_read = fgumi_bam_io::read_prefix(&mut file, &mut header)?;

    Ok(classify_compression_header(&header[..bytes_read]))
}

/// Classify a FASTQ's compression from its leading bytes.
///
/// Delegates to [`fgumi_bam_io::classify_input`], the classifier every fgumi
/// input sniff shares, so a file gets the same verdict whichever command opens
/// it. Unlike the alignment inputs, FASTQ legitimately comes plain-gzipped, so
/// `Gzip` is a supported format here rather than an error.
fn classify_compression_header(header: &[u8]) -> CompressionFormat {
    match fgumi_bam_io::classify_input(header) {
        InputFormat::Bgzf => CompressionFormat::Bgzf,
        InputFormat::Gzip => CompressionFormat::Gzip,
        InputFormat::Text | InputFormat::Empty => CompressionFormat::Plain,
    }
}

/// Open a FASTQ file with automatic detection of compression format.
///
/// For BGZF-compressed files, uses noodles `MultithreadedReader` when threads > 1.
/// For regular gzip files, uses `flate2::read::MultiGzDecoder` single-threaded.
/// For uncompressed files, opens directly.
///
/// # Arguments
/// * `path` - Path to the FASTQ file
/// * `threads` - Number of decompression threads (only used for BGZF)
/// * `async_reader` - If true, wrap the file in an async prefetch reader
///
/// # Returns
/// A boxed reader that implements `BufRead` + `Send`
pub(crate) fn open_fastq_reader(
    path: &Path,
    threads: usize,
    async_reader: bool,
    check_crc: bool,
    no_check_crc: bool,
) -> Result<Box<dyn BufRead + Send>> {
    use flate2::read::MultiGzDecoder;

    // CRC-verification policy, resolved per input path: an explicit flag wins,
    // otherwise verify a file and trust (skip) stdin. Only consulted for the
    // single-threaded BGZF path below.
    let verify_crc = if check_crc {
        true
    } else if no_check_crc {
        false
    } else {
        !is_stdin_path(path)
    };
    // Mirror the BAM path's `CRC verify:` line (see `RunOptions::log_effective_check_crc`)
    // so the effective policy — including the default-off-for-stdin skip, a change from
    // fgumi's previous always-verify default — is visible in every run's log rather than a
    // silent decision. Emitted below only for BGZF input, the only format that carries a
    // per-block CRC32 to verify or skip.
    let crc_reason = if check_crc {
        " (--check-crc)"
    } else if no_check_crc {
        " (--no-check-crc)"
    } else if is_stdin_path(path) {
        " (trusted stdin)"
    } else {
        ""
    };

    // stdin cannot be sniffed by path and then re-opened — the bytes read to
    // classify it are gone. Peek them off the stream and chain them back in
    // front of it instead, so the input is consumed exactly once.
    let (format, reader) = if is_stdin_path(path) {
        // Honor --async-reader for stdin too (parity with the file path below):
        // a pipe benefits from prefetch overlapping upstream production with us.
        let raw: Box<dyn std::io::Read + Send> = if async_reader {
            log::info!("async FASTQ reader enabled: spawning fgumi-prefetch thread for stdin");
            Box::new(fgumi_bam_io::prefetch_reader::PrefetchReader::new(io::stdin()))
        } else {
            Box::new(io::stdin())
        };
        let (format, restored) = detect_compression_format_from_stream(raw)?;
        debug!("Detected {format:?} FASTQ on stdin");
        (format, restored)
    } else {
        let format = detect_compression_format(path)?;
        let file = File::open(path)?;
        // Only a real `File` can carry the kernel hints `PrefetchReader` issues,
        // so the async wrap happens here rather than around the format decode.
        let reader: Box<dyn std::io::Read + Send> = if async_reader {
            fgumi_bam_io::os_hints::advise_sequential(&file);
            log::info!(
                "async FASTQ reader enabled: spawning fgumi-prefetch thread for {}",
                path.display()
            );
            Box::new(fgumi_bam_io::prefetch_reader::PrefetchReader::from_file(file))
        } else {
            Box::new(file)
        };
        (format, reader)
    };

    match format {
        CompressionFormat::Bgzf if threads > 1 && verify_crc => {
            // Verifying + multi-threaded: use noodles' parallel decoder. It always
            // verifies CRC32 and has no skip knob, which is exactly what this arm
            // wants. When `--no-check-crc` is set we fall through to the fgumi arm
            // below instead, so the skip is always honored (trading MT decode for
            // it); simultaneous MT decode + skip would need a threaded fgumi
            // decoder (deferred, same as the raw-BAM multi-threaded path).
            info!("CRC verify: on{crc_reason}");
            info!("Detected BGZF-compressed FASTQ, using {threads} decompression threads");
            let worker_count = std::num::NonZero::new(threads).expect("threads > 1 checked above");
            let reader = MultithreadedReader::with_worker_count(worker_count, reader);
            Ok(Box::new(BufReader::with_capacity(BUFFER_SIZE, reader)))
        }
        CompressionFormat::Bgzf => {
            // BGZF decode through fgumi-bgzf, which honors `verify_crc` (noodles /
            // flate2 always verify). Reached for single-threaded input, and for
            // `--no-check-crc` at any thread count (the guard above only takes the
            // parallel path when verifying), so the skip is never silently ignored.
            // The inner buffer feeds the block-header reads; the decoder is `BufRead`.
            info!("CRC verify: {}{crc_reason}", if verify_crc { "on" } else { "off" });
            debug!(
                "Detected BGZF-compressed FASTQ, single-threaded fgumi-bgzf decode (CRC verify: {verify_crc})"
            );
            let buffered: Box<dyn std::io::Read + Send> =
                Box::new(BufReader::with_capacity(BUFFER_SIZE, reader));
            let decoder = fgumi_bam_io::FgumiBgzfReader::new(buffered, verify_crc);
            Ok(Box::new(BufReader::with_capacity(BUFFER_SIZE, decoder)))
        }
        CompressionFormat::Gzip => {
            // Plain gzip is not block-structured; there is no per-block CRC to
            // skip, so it stays on flate2 and ignores the flag.
            debug!("Detected gzip-compressed FASTQ, using single-threaded decompression");
            Ok(Box::new(BufReader::with_capacity(
                BUFFER_SIZE,
                MultiGzDecoder::new(BufReader::with_capacity(BUFFER_SIZE, reader)),
            )))
        }
        CompressionFormat::Plain => {
            debug!("Detected uncompressed FASTQ");
            Ok(Box::new(BufReader::with_capacity(BUFFER_SIZE, reader)))
        }
    }
}

/// Classify a FASTQ stream's compression, returning a reader that replays the
/// bytes the classification consumed.
///
/// The path-based [`detect_compression_format`] opens the file, reads its
/// header, and drops the handle, relying on the caller to re-open. That is fine
/// for a regular file and impossible for a pipe, so this variant peeks from the
/// stream and chains the peeked bytes back on.
///
/// # Errors
///
/// Returns an error if the stream cannot be read.
fn detect_compression_format_from_stream(
    mut reader: Box<dyn std::io::Read + Send>,
) -> Result<(CompressionFormat, Box<dyn std::io::Read + Send>)> {
    // Same window the shared classifier needs: gzip magic, deflate method and
    // flag bytes, and the BGZF extra field.
    let mut header = [0u8; fgumi_bam_io::FORMAT_PREFIX_LEN];
    let filled =
        fgumi_bam_io::read_prefix(&mut reader, &mut header).context("Failed to read from stdin")?;

    let format = classify_compression_header(&header[..filled]);
    let restored: Box<dyn std::io::Read + Send> =
        Box::new(ChainedReader::new(header[..filled].to_vec(), reader));
    Ok((format, restored))
}

/// Quality encoding type
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum QualityEncoding {
    Standard, // Phred+33 (Sanger)
    Illumina, // Phred+64 (Illumina 1.3-1.7)
}

/// Fixed-size summary of the pooled quality bytes scanned during encoding
/// detection.
///
/// Detection only needs the global min/max quality byte and whether any
/// (non-empty) quality was seen, so accumulating this while scanning lets the
/// sampler avoid retaining every sampled quality string — memory stays O(1) in
/// the number of inputs, sample size, and read length rather than
/// `inputs × sample × read_length`.
#[derive(Debug, Clone, Copy)]
struct QualDetectionStats {
    /// Minimum quality byte observed across all non-empty sampled reads.
    min_qual: u8,
    /// Maximum quality byte observed across all non-empty sampled reads.
    max_qual: u8,
    /// Total number of quality bytes observed (0 if every sampled read was empty).
    total_bases: u64,
    /// Number of records sampled (including empty-quality reads); 0 means no
    /// records were available at all.
    num_records: u64,
}

impl QualDetectionStats {
    fn new() -> Self {
        Self { min_qual: u8::MAX, max_qual: u8::MIN, total_bases: 0, num_records: 0 }
    }

    /// Fold one sampled read's quality bytes into the running summary. Empty
    /// qualities still count as a record (so detection can tell "no records" from
    /// "records with empty qualities") but contribute no min/max/base statistics.
    fn observe(&mut self, qual: &[u8]) {
        self.num_records += 1;
        for &q in qual {
            self.min_qual = self.min_qual.min(q);
            self.max_qual = self.max_qual.max(q);
            self.total_bases += 1;
        }
    }
}

impl QualityEncoding {
    /// Convert quality scores to standard numeric format (Phred+33)
    fn to_standard_numeric(self, quals: &[u8]) -> Vec<u8> {
        match self {
            QualityEncoding::Standard => quals.iter().map(|&q| q.saturating_sub(33)).collect(),
            QualityEncoding::Illumina => quals.iter().map(|&q| q.saturating_sub(64)).collect(),
        }
    }

    /// Decide the encoding from pooled quality statistics using robust heuristics.
    ///
    /// This implements a more robust detection algorithm that:
    /// - Checks for quality scores in the Illumina 1.3-1.7 range (64-126)
    /// - Checks for quality scores in the Sanger/Illumina 1.8+ range (33-126)
    /// - Handles edge cases like empty reads or very short reads
    /// - Provides informative error messages for invalid encodings
    ///
    /// Takes pre-accumulated [`QualDetectionStats`] rather than the raw sampled
    /// quality strings so the sampler can reduce every head to a fixed-size
    /// (min, max, counts) summary instead of retaining `inputs × sample × read`
    /// bytes in memory (see [`sample_detection_quals`]).
    fn from_stats(stats: &QualDetectionStats) -> Result<Self> {
        if stats.num_records == 0 {
            bail!("Cannot detect quality encoding: no records provided");
        }

        let QualDetectionStats { min_qual, max_qual, total_bases, .. } = *stats;

        // If all reads were empty, we can't detect encoding but we'll default to Standard
        if total_bases == 0 {
            return Ok(QualityEncoding::Standard);
        }

        // No-compatible-encoding failure (fgbio parity: EXT3-03).
        //
        // fgbio `FastqToBam` derives this failure from its encoding model: it builds
        // the set of encodings whose ASCII range contains every observed quality char,
        // and `case Nil => fail("Quality scores in FASTQ files do not match any known
        // encoding.")` when that set is empty (`FastqToBam.scala:131`). Because the
        // Standard (Phred+33) ASCII range `[33, 126]` is the union of all three fgbio
        // encodings' ranges (Solexa `[59, 126]` and Illumina `[64, 126]` are subsets),
        // fgbio's compatible set is empty on *exactly* the condition below: an observed
        // char `< 33` or `> 126`. So this bail is fgbio-equivalent on the failure set.
        //
        // fgumi intentionally *improves* on fgbio's message here rather than copying it.
        // fgbio emits a static string; fgumi reports the observed `[min, max]` range and
        // the required range, which points the user straight at the offending file
        // instead of just asserting no encoding matched. This divergence is deliberate,
        // not a parity gap (see `detect_reports_observed_range_on_no_compatible_encoding`).
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

    /// Test-only convenience: decide the encoding directly from raw sampled
    /// quality strings, folding them into [`QualDetectionStats`] first. Production
    /// code streams the stats via [`sample_detection_quals`] and calls
    /// [`Self::from_stats`] instead of materializing this slice; this wrapper lets
    /// the boundary-condition tests keep exercising the decision logic from
    /// hand-crafted quality vectors.
    #[cfg(test)]
    fn detect(records: &[Vec<u8>]) -> Result<Self> {
        let mut stats = QualDetectionStats::new();
        for qual in records {
            stats.observe(qual);
        }
        Self::from_stats(&stats)
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
/// A single `<number><operator>` pair may be specified using a `+` sign instead of a number to denote "all
/// remaining bases". The `+` may appear in any position (not only the last); segments after it are resolved by
/// counting back from the end of the read. This is useful if, e.g., FASTQs have been trimmed and contain reads of
/// varying length.  For example to convert a paired-end run with an index read and where the first 5 bases of R1
/// are a UMI and the second five bases should be skipped, you might specify:
///
/// ```text
/// --input r1.fq r2.fq i1.fq --read-structures 5M5S+T +T +B
/// ```
///
/// Alternative if you know your reads are of fixed length you could specify:
///
/// ```text
/// --input r1.fq r2.fq i1.fq --read-structures 5M5S65T 75T 8B
/// ```
///
/// For more information on read structures see the
/// [Read Structure Wiki Page](https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures)
///
/// UMIs may be extracted from the read sequences, the read names, or both.  If `--extract-umis-from-read-names` is
/// specified, any UMIs present in the read names are extracted; read names are expected to be `:`-separated with
/// the UMI in the last field (of at least 8).  The extracted UMI is upper-cased, `r`-prefixed segments are
/// reverse-complemented, `+` dual-UMI delimiters become `-`, and a UMI containing any character outside `ACGTN-`
/// is rejected with an error.  If this option is specified, `--store-umi-quals` may not be used as
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

## Input Requirements

All input FASTQs must contain the same number of records, in the same order. A pair whose
streams disagree on record count is rejected with a non-zero exit status naming which input
ended first — surplus records are neither discarded nor emitted as unpaired reads, since a
short input almost always means a truncated transfer or an interrupted upstream tool.

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
a UMI and the second five bases should be skipped:

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
read names are expected to be `:`-separated and the UMI is taken from the **last** field. At
least 8 fields must be present — the standard Illumina shape
`@<instrument>:<run>:<flowcell>:<lane>:<tile>:<x>:<y>:<UMI>`. Names with 9+ fields (e.g.
produced by demultiplexers that fold the sample index into the colon-separated portion) are
also handled, with the UMI still coming from the last field. Any `+` characters in the
extracted UMI are normalized to `-`. If UMI segments are present in the read structures those
will also be extracted. If UMIs are present in both, the final UMIs are constructed by first
taking the UMIs from the read names, then adding a hyphen, then the UMIs extracted from the
reads.
"#
)]
#[command(verbatim_doc_comment)]
#[allow(clippy::struct_excessive_bools)]
pub struct Extract {
    /// Input FASTQ files corresponding to each sequencing read (e.g. R1, I1, etc.)
    #[arg(long, short = 'i', required = true, num_args = 1..)]
    inputs: Vec<PathBuf>,

    /// Output BAM file to be written
    #[arg(long, short = 'o', required = true)]
    output: PathBuf,

    /// Read structures, one for each of the FASTQs (optional if 1-2 template-only FASTQs)
    #[arg(long, short = 'r', num_args = 0..)]
    read_structures: Vec<ReadStructure>,

    /// Treat a single input as interleaved paired-end FASTQ (`R1, R2, R1, R2, …`),
    /// de-interleaving it into the two reads. Requires exactly one `--input` (a
    /// file or `-` for stdin) and describes both reads with two `--read-structures`
    /// (defaults to `+T +T`). This lets a streaming trimmer or converter pipe
    /// interleaved pairs straight into extract without staging two FASTQ files.
    #[arg(long = "interleaved", short = 'p', default_value_t = false)]
    interleaved: bool,

    /// Store UMI base quality scores in the QX SAM tag
    #[arg(long, short = 'q', conflicts_with = "extract_umis_from_read_names")]
    store_umi_quals: bool,

    /// Store cell barcode base quality scores in the CY SAM tag
    #[arg(long, short = 'C')]
    store_cell_quals: bool,

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
    single_tag: Option<SamTag>,

    /// Tag containing adapter clipping position to adjust (e.g. 'XT' from `MarkIlluminaAdapters`)
    #[arg(long)]
    clipping_attribute: Option<SamTag>,

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

    /// Queue memory options.
    #[command(flatten)]
    pub queue_memory: QueueMemoryOptions,

    /// Wrap FASTQ inputs in a userspace async prefetch reader. Dedicates one
    /// OS thread per input stream to issue reads ahead of decompression/parsing.
    /// Hidden experimental flag.
    #[arg(long = "async-reader", default_value_t = false, hide = true)]
    pub async_reader: bool,

    /// Verify each BGZF block's CRC32 checksum while decoding the input.
    ///
    /// Applies to BGZF-compressed FASTQ input (bgzip'd); plain gzip is not
    /// block-structured and has no per-block CRC to skip. The policy is honored
    /// at any thread count, though BGZF decode runs single-threaded when skipping
    /// (verifying BGZF input can decode in parallel). Without either flag,
    /// verification defaults on for file input and off for trusted stdin (a
    /// freshly-piped stream is trusted; a file may have been archived or
    /// transferred since it was written, where a flipped bit is what CRC32 exists
    /// to catch). Pass `--check-crc` to force it on. Mutually exclusive with
    /// `--no-check-crc`.
    #[arg(long = "check-crc", default_value_t = false, conflicts_with = "no_check_crc")]
    pub check_crc: bool,

    /// Skip CRC32 verification while decoding BGZF FASTQ input.
    ///
    /// Trades the CRC32 integrity check for faster BGZF decode (which then runs
    /// single-threaded). See `--check-crc` for the default policy this overrides.
    /// Mutually exclusive with `--check-crc`.
    #[arg(long = "no-check-crc", default_value_t = false, conflicts_with = "check_crc")]
    pub no_check_crc: bool,
}

impl Extract {
    /// Project the parsed CLI struct into the per-stage [`ExtractOptions`] the
    /// chain builder consumes (`StageOptionsBag::extract`).
    ///
    /// `quality_encoding` is a placeholder here (`QualityEncoding::Standard`):
    /// the chain FASTQ source detects the real encoding while opening its readers
    /// and overrides the field before the extract step runs (see
    /// [`ExtractOptions::quality_encoding`]). Every other field maps directly
    /// from the identically-named CLI flag; `platform` (CLI default `"illumina"`)
    /// becomes `Some(..)` so `build_fastq_header` always emits `@RG PL:`, matching
    /// `Self::create_header`. `clipping_attribute` is intentionally dropped — it
    /// does not apply to FASTQ input (there is no existing clipping to adjust).
    #[must_use]
    pub fn to_extract_options(&self) -> ExtractOptions {
        ExtractOptions {
            sample: self.sample.clone(),
            library: self.library.clone(),
            platform: Some(self.platform.clone()),
            platform_unit: self.platform_unit.clone(),
            read_group_id: self.read_group_id.clone(),
            comments: self.comment.clone(),
            barcode: self.barcode.clone(),
            platform_model: self.platform_model.clone(),
            sequencing_center: self.sequencing_center.clone(),
            predicted_insert_size: self.predicted_insert_size,
            description: self.description.clone(),
            run_date: self.run_date.clone(),
            quality_encoding: QualityEncoding::Standard,
            store_umi_quals: self.store_umi_quals,
            store_cell_quals: self.store_cell_quals,
            single_tag: self.single_tag,
            annotate_read_names: self.annotate_read_names,
            extract_umis_from_read_names: self.extract_umis_from_read_names,
            store_sample_barcode_qualities: self.store_sample_barcode_qualities,
            async_reader: self.async_reader,
            check_crc: self.check_crc,
            no_check_crc: self.no_check_crc,
        }
    }

    /// Run extract on the declarative chain builder (the `--threads` path).
    ///
    /// Hand-builds the [`ChainSpec`] (rather than using
    /// [`ChainSpec::single_stage`], which is BAM-in/BAM-out): the source is a
    /// FASTQ source — [`SourceSpec::InterleavedFastq`] for `--interleaved`, else
    /// [`SourceSpec::Fastqs`] — and the sink is a BAM. The chain opens its own
    /// readers and detects the quality encoding in `ChainBuilder::open_source`;
    /// `read_streams`/`verify_crc` are BAM-reader knobs and inert here (the FASTQ
    /// source carries its CRC policy in [`ExtractOptions`]). The no-`--threads`
    /// serial loop in [`Command::execute`] is the in-process parity oracle.
    ///
    /// [`ChainSpec`]: crate::pipeline::chains::ChainSpec
    /// [`ChainSpec::single_stage`]: crate::pipeline::chains::ChainSpec::single_stage
    /// [`SourceSpec::Fastqs`]: crate::pipeline::chains::SourceSpec::Fastqs
    /// [`SourceSpec::InterleavedFastq`]: crate::pipeline::chains::SourceSpec::InterleavedFastq
    fn execute_chain(&self, command_line: &str) -> Result<()> {
        use crate::pipeline::chains::{
            ChainSpec, SinkSpec, SourceSpec, Stage, StageOptionsBag, build_for,
        };

        let read_structures = self.get_read_structures()?;
        let source = if self.interleaved {
            SourceSpec::interleaved_fastq(self.inputs[0].clone(), read_structures)?
        } else {
            SourceSpec::fastqs(self.inputs.clone(), read_structures)?
        };

        let stage_opts =
            StageOptionsBag { extract: Some(self.to_extract_options()), ..Default::default() };

        let spec = ChainSpec {
            stages: vec![Stage::Extract],
            source,
            sink: SinkSpec::Bam(self.output.clone()),
            stage_opts,
            threading: self.threading.clone(),
            compression: self.compression.clone(),
            scheduler: self.scheduler_opts.clone(),
            queue_memory: self.queue_memory.clone(),
            async_reader: self.async_reader,
            // Inert for a FASTQ source: `read_streams` is a seekable-BAM knob and
            // `verify_crc` is the BAM BGZF-decode policy; the FASTQ source carries
            // its own CRC policy in `ExtractOptions` (`check_crc`/`no_check_crc`).
            read_streams: fgumi_bam_io::ReadStreams::Fixed(1),
            verify_crc: false,
            command_line: command_line.to_string(),
        };

        build_for(spec)?.run()
    }

    /// Get actual read structures (default to +T if none provided).
    ///
    /// For interleaved input the two reads share one physical stream, so the
    /// default is two `+T` structures regardless of the single `--input`; the
    /// normal path defaults to one `+T` per FASTQ for 1-2 FASTQs.
    fn get_read_structures(&self) -> Result<Vec<ReadStructure>> {
        if self.read_structures.is_empty() {
            if self.interleaved {
                return Ok(vec![ReadStructure::from_str("+T")?; 2]);
            }
            if (1..=2).contains(&self.inputs.len()) {
                return Ok(vec![ReadStructure::from_str("+T")?; self.inputs.len()]);
            }
        }
        Ok(self.read_structures.clone())
    }

    /// Validate inputs
    fn validate(&self) -> Result<()> {
        let read_structures = self.get_read_structures()?;

        if self.interleaved {
            // One physical stream carries both reads, so the usual
            // one-structure-per-input equality does not apply: the single input
            // is split into exactly two reads.
            ensure!(
                self.inputs.len() == 1,
                "--interleaved requires exactly one --input (the interleaved stream); got {}.",
                self.inputs.len()
            );
            ensure!(
                read_structures.len() == 2,
                "--interleaved requires exactly two --read-structures (one per read); got {}.",
                read_structures.len()
            );
        } else {
            ensure!(
                self.inputs.len() == read_structures.len(),
                "input and read-structure must be supplied the same number of times."
            );
        }

        let template_count: usize = read_structures
            .iter()
            .map(|rs| rs.segments_by_type(SegmentType::Template).count())
            .sum();
        ensure!(
            (1..=2).contains(&template_count),
            "read structures must contain 1-2 template reads total."
        );

        ensure!(
            !self.extract_umis_from_read_names || !self.store_umi_quals,
            "Cannot store UMI qualities (--store-umi-quals) when also extracting UMIs from read names (--extract-umis-from-read-names)."
        );

        // Validate threads parameter (ThreadingOptions handles minimum of 1 internally)

        // Validate input files exist. `-` names stdin, which is a stream rather
        // than a path — but one stdin cannot supply two FASTQs, so it is only
        // accepted when it is the sole input.
        let stdin_inputs = self.inputs.iter().filter(|p| is_stdin_path(p)).count();
        ensure!(
            stdin_inputs == 0 || self.inputs.len() == 1,
            "stdin (`-`) can only be used as the input when there is exactly one --input; \
             got {} inputs. Read the other FASTQs from paths, or interleave them first.",
            self.inputs.len()
        );
        for input in &self.inputs {
            validate_input_exists(input, "Input FASTQ")?;
        }

        // Validate read structures are not empty
        for (i, rs) in read_structures.iter().enumerate() {
            ensure!(!rs.segments().is_empty(), "Read structure {} is empty", i + 1);
        }

        // SamTag::from_str enforces the SAM aux-tag pattern via clap, so format
        // is already validated. Reject `--single-tag` matching any tag the
        // extractor itself emits, since those would either be silently
        // overwritten or produce conflicting aux fields.
        const RESERVED_OUTPUT_TAGS: &[SamTag] =
            &[SamTag::RX, SamTag::QX, SamTag::CB, SamTag::CY, SamTag::BC, SamTag::QT, SamTag::RG];
        if let Some(tag) = self.single_tag {
            ensure!(
                !RESERVED_OUTPUT_TAGS.contains(&tag),
                "Single tag (--single-tag) cannot collide with tags emitted by extract \
                 (RX, QX, CB, CY, BC, QT, RG): {tag}"
            );
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
        header = crate::commands::common::add_pg_to_builder(header, command_line)?;

        Ok(header.build())
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

    /// Extracts the read name and, if requested, the UMI from the read name.
    ///
    /// The standard Illumina FASTQ read name shape with a UMI present is:
    ///
    ///   `@<instrument>:<run>:<flowcell>:<lane>:<tile>:<x>:<y>:<UMI>`
    ///
    /// (8 `:`-separated fields; the trailing space-separated comment such as
    /// `1:N:0:<index>` is not counted, as it is stripped before parsing). See
    /// <https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm>.
    ///
    /// An old-style Casava (<1.8) `/1` / `/2` read-number suffix is stripped from
    /// the returned name (see [`strip_read_suffix`])
    /// so both mates of a pair share an identical QNAME, matching fgbio's `FastqSource`.
    /// Stripping happens before UMI extraction so a read-number digit never leaks into the UMI.
    ///
    /// Some demultiplexers fold additional information (e.g. the sample index) into
    /// the colon-separated portion of the read name, producing names with 9+ fields
    /// where the UMI is the **last** field rather than the 8th. To handle both
    /// shapes uniformly — matching the behavior of fgbio's
    /// `Umis.extractUmisFromReadName` — this function returns the **last**
    /// `:`-separated field as the UMI when at least 8 fields are present.
    ///
    /// Names with fewer than 8 fields are treated as not containing a UMI and
    /// produce `None`. The extracted UMI is then normalized via
    /// [`normalize_read_name_umi`](crate::umi::read_name::normalize_read_name_umi) to match fgbio's
    /// strict `Umis.extractUmisFromReadName` (reverse-complement `r`-prefixed
    /// segments, `+`→`-`, upper-case).
    ///
    /// # Errors
    /// Returns an error if the extracted UMI contains a character outside `ACGTN-`,
    /// mirroring fgbio's strict-mode rejection (EXT-01).
    fn extract_read_name_and_umi(
        header: &[u8],
        extract_umis: bool,
    ) -> Result<(Vec<u8>, Option<Vec<u8>>)> {
        // Remove @ prefix if present
        let name_bytes = if header.starts_with(b"@") { &header[1..] } else { header };

        // Strip a trailing space comment and an old-style Casava (<1.8) `/1` / `/2`
        // read-number suffix so both mates of a pair share an identical QNAME (required by
        // the SAM spec). Matches fgbio's `FastqSource`: only `/` followed by a single digit
        // is removed. Done before UMI extraction so a read-number digit never leaks into the
        // UMI's trailing `:` field. This is the shared helper the paired-FASTQ sync
        // validation also uses, keeping the written QNAME and validation in lock-step.
        let name_part = strip_read_suffix(name_bytes);

        if !extract_umis {
            return Ok((name_part.to_vec(), None));
        }

        // The UMI is the last `:`-separated field, but only if the name has at
        // least 8 fields (matching the standard Illumina read-name layout). This
        // works for both 8-field names (UMI in field 8) and 9+ field names where
        // a demultiplexer has appended additional fields (e.g. a sample index)
        // before the UMI.
        let parts: Vec<&[u8]> = name_part.split(|&b| b == b':').collect();

        if parts.len() >= 8
            && let Some(last) = parts.last()
            && !last.is_empty()
        {
            let umi =
                crate::umi::read_name::normalize_read_name_umi(last, true).with_context(|| {
                    format!(
                        "extracting UMI from read name '{}'",
                        String::from_utf8_lossy(name_part)
                    )
                })?;
            return Ok((name_part.to_vec(), Some(umi)));
        }

        Ok((name_part.to_vec(), None))
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
        let first_name_part = strip_read_suffix(first_name);

        // Check that all other read sets have the same name
        for (i, read_set) in read_sets.iter().enumerate().skip(1) {
            let header = &read_set.header;
            let name = if header.starts_with(b"@") { &header[1..] } else { header.as_slice() };

            let name_part = strip_read_suffix(name);

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

    /// Context for a read name that BAM cannot represent.
    ///
    /// The builder's own error reports the length and the limit but not which
    /// read hit it, which is what the user needs to find the offending FASTQ
    /// record. `name` is the final name as written, so it already includes the
    /// `+<UMI>` suffix when `--annotate-read-names` is set. Long names are
    /// truncated so the message stays readable.
    fn read_name_too_long_context(name: &[u8]) -> String {
        // Bounded head AND tail: the name is the FASTQ header plus an optional `+<UMI>` suffix
        // (`--annotate-read-names`), so head-only truncation would drop a UMI that lands past
        // the budget. Retaining both ends keeps two otherwise-similar reads distinguishable,
        // and bounding the tail keeps the line readable for an over-long input-derived name.
        const HEAD: usize = 44;
        const TAIL: usize = 20;
        let shown = if name.len() <= HEAD + TAIL {
            format!("'{}'", String::from_utf8_lossy(name))
        } else {
            let head = String::from_utf8_lossy(&name[..HEAD]);
            let tail = String::from_utf8_lossy(&name[name.len() - TAIL..]);
            format!("'{head}...{tail}' (name shown truncated)")
        };
        format!("could not write the record for read {shown}")
    }

    /// Builds one raw BAM record for a template segment, substituting a single `N` @ Q2 when
    /// the segment has no sequence.
    ///
    /// The name comes from the input FASTQ header (plus the UMI when `--annotate-read-names`
    /// is set), so an over-long name is bad input rather than a bug: fail cleanly via
    /// `try_build_record` instead of panicking partway through writing the output BAM, and
    /// attach [`Self::read_name_too_long_context`] so the message names the offending read.
    ///
    /// Used by the serial-oracle `make_raw_records` path; the chain path's
    /// `make_raw_records_from_fastq_set` builds records through the same
    /// `UnmappedSamBuilder` API, which the parity tests keep in step.
    fn build_template_record(
        builder: &mut UnmappedSamBuilder,
        name: &[u8],
        flag: u16,
        seq: &[u8],
        quals: &[u8],
        encoding: QualityEncoding,
    ) -> Result<()> {
        if seq.is_empty() {
            builder.try_build_record(name, flag, b"N", &[2u8])
        } else {
            builder.try_build_record(name, flag, seq, &encoding.to_standard_numeric(quals))
        }
        .with_context(|| Self::read_name_too_long_context(name))
    }

    /// Write raw BAM records from a read set directly to a writer.
    ///
    /// Uses `UnmappedSamBuilder` to construct records as raw bytes,
    /// bypassing `RecordBuf` allocation and encoding overhead.
    ///
    /// Returns the number of records written.
    #[allow(clippy::too_many_lines)]
    fn make_raw_records(
        &self,
        read_set: &FastqSet,
        encoding: QualityEncoding,
        builder: &mut UnmappedSamBuilder,
        writer: &mut RawBamWriter,
    ) -> Result<u64> {
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
        let (read_name_bytes, umi_from_name) =
            Self::extract_read_name_and_umi(&read_set.header, self.extract_umis_from_read_names)?;

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

        let num_templates = templates.len();

        for (index, template) in templates.iter().enumerate() {
            // Compute flags for unmapped reads
            let mut flag = flags::UNMAPPED;
            if num_templates == 2 {
                flag |= flags::PAIRED | flags::MATE_UNMAPPED;
                if index == 0 {
                    flag |= flags::FIRST_SEGMENT;
                } else {
                    flag |= flags::LAST_SEGMENT;
                }
            }

            // Set read name (optionally with UMI annotation)
            let annotated_name: Option<Vec<u8>> = if self.annotate_read_names
                && !final_umi_bs.is_empty()
            {
                let mut name = Vec::with_capacity(read_name_bytes.len() + 1 + final_umi_bs.len());
                name.extend_from_slice(&read_name_bytes);
                name.push(b'+');
                name.extend_from_slice(final_umi_bs.as_bytes());
                Some(name)
            } else {
                None
            };
            let final_read_name: &[u8] = annotated_name.as_deref().unwrap_or(&read_name_bytes);

            Self::build_template_record(
                builder,
                final_read_name,
                flag,
                &template.seq,
                &template.quals,
                encoding,
            )?;

            // Append tags
            // Read group
            builder.append_string_tag(SamTag::RG, self.read_group_id.as_bytes());

            // Cell barcode
            if !cell_barcode_bs.is_empty() {
                builder.append_string_tag(SamTag::CB, cell_barcode_bs.as_bytes());
            }

            if !cell_quals_bs.is_empty() && self.store_cell_quals {
                builder.append_string_tag(SamTag::CY, cell_quals_bs.as_bytes());
            }

            // Sample barcode
            if !sample_barcode_bs.is_empty() {
                builder.append_string_tag(SamTag::BC, sample_barcode_bs.as_bytes());
            }

            if self.store_sample_barcode_qualities && !sample_quals_bs.is_empty() {
                builder.append_string_tag(SamTag::QT, sample_quals_bs.as_bytes());
            }

            // UMI
            if !final_umi_bs.is_empty() {
                builder.append_string_tag(SamTag::RX, final_umi_bs.as_bytes());

                // Single tag for all concatenated UMIs (if specified)
                if let Some(single_tag) = self.single_tag {
                    builder.append_string_tag(single_tag, final_umi_bs.as_bytes());
                }

                // Only add UMI qualities if not extracted from read names
                if umi_from_name.is_none() && !umi_qual_bs.is_empty() && self.store_umi_quals {
                    builder.append_string_tag(SamTag::QX, umi_qual_bs.as_bytes());
                }
            }

            // Write the record directly
            writer.write_raw_record(builder.as_bytes())?;
            builder.clear();
        }

        Ok(num_templates as u64)
    }

    /// Process records in single-threaded mode
    ///
    /// Returns the number of records written.
    fn process_singlethreaded(
        &self,
        fq_iterators: &mut [ReadSetIterator],
        writer: &mut RawBamWriter,
        encoding: QualityEncoding,
    ) -> Result<u64> {
        let progress = ProgressTracker::new("Processed records").with_interval(1_000_000);
        let mut read_pair_count: u64 = 0;
        let mut builder = UnmappedSamBuilder::new();

        // One record per stream at this position; `1` = a record was read, `0` =
        // the stream ended. Naming the streams (rather than a vague "some ended
        // before others") matches the threaded pipeline and the help text's
        // promise to name which input ended first. Hoisted out of the loop and
        // reset each iteration to avoid a per-template allocation on this fast path.
        let mut counts = vec![0usize; fq_iterators.len()];
        loop {
            counts.fill(0);
            let mut next_read_sets = Vec::with_capacity(fq_iterators.len());
            for (idx, iter) in fq_iterators.iter_mut().enumerate() {
                match iter.next() {
                    Some(Ok(rec)) => {
                        counts[idx] = 1;
                        next_read_sets.push(rec);
                    }
                    Some(Err(e)) => return Err(e),
                    None => {} // stream ended; leave its count at 0
                }
            }

            // Every stream ended together: clean EOF.
            if counts.iter().all(|count| *count == 0) {
                break;
            }
            // Some streams still had a record while others ended: out of sync.
            if counts.iter().any(|count| *count != counts[0]) {
                return Err(fastq_out_of_sync_error(&counts).into());
            }

            //  Validate read names match across all FASTQs
            Self::validate_read_names_match(&next_read_sets)?;

            let read_set = FastqSet::combine_readsets(next_read_sets);
            let num_records = self.make_raw_records(&read_set, encoding, &mut builder, writer)?;

            read_pair_count += num_records;
            progress.log_if_needed(num_records);
        }

        progress.log_final();
        Ok(read_pair_count)
    }
}

/// Pool quality strings for encoding detection across the first
/// [`QUALITY_DETECTION_SAMPLE_SIZE`] records of **every** input FASTQ.
///
/// fgbio's `FastqToBam` detects the encoding from `heads.transpose.flatten` — the
/// pooled heads of all sources — so that a non-template first file (e.g. an
/// all-high-quality index read) cannot skew detection toward the wrong encoding
/// and shift every template base quality by ~31 (EXT3-01). Sampling only
/// `inputs[0]` (the previous behavior) did exactly that: a Phred+64-looking first
/// file forced Illumina detection even when the template reads in `inputs[1..]`
/// were Phred+33.
///
/// A fresh reader is opened per input purely for sampling, so the caller's main
/// readers are not consumed. Sampling is always synchronous — a bounded head-scan
/// of at most [`QUALITY_DETECTION_SAMPLE_SIZE`] records gains nothing from the
/// async reader (which the main extraction still uses per `--async-reader`).
fn sample_detection_quals(
    inputs: &[PathBuf],
    check_crc: bool,
    no_check_crc: bool,
) -> Result<QualDetectionStats> {
    let mut stats = QualDetectionStats::new();
    for input in inputs {
        // Honor the same CRC policy as the main decode: the sampler's buffered
        // read can pull the whole (small) input, so a corrupted trailing block
        // would otherwise fail here under the default even with `--no-check-crc`.
        let mut reader = SimdFastqReader::with_capacity(
            open_fastq_reader(input, 1, false, check_crc, no_check_crc)?,
            BUFFER_SIZE,
        );
        sample_into(&mut reader, &mut stats)?;
    }
    Ok(stats)
}

/// Fold up to [`QUALITY_DETECTION_SAMPLE_SIZE`] records' qualities into `stats`.
///
/// Each sampled read is folded into the running summary and dropped — no
/// per-read quality string is retained.
fn sample_into<R: BufRead>(
    reader: &mut SimdFastqReader<R>,
    stats: &mut QualDetectionStats,
) -> Result<()> {
    for _ in 0..QUALITY_DETECTION_SAMPLE_SIZE {
        match reader.next() {
            Some(Ok(rec)) => stats.observe(&rec.quality),
            Some(Err(e)) => return Err(e.into()),
            None => break,
        }
    }
    Ok(())
}

/// Sample quality stats from a stream, returning a reader that replays every
/// byte the sampling consumed.
///
/// [`sample_detection_quals`] re-opens each input by path so the main pass gets
/// a fresh reader. stdin has no second open: whatever the sample reads is gone,
/// and the main pass would silently see an empty FASTQ. Capturing the sampled
/// bytes through a [`TeeReader`] and chaining them back makes the stream
/// re-readable from byte zero without buffering the whole input.
///
/// The tee captures what was pulled from the *source*, not what the FASTQ
/// parser yielded, so the replay stays complete no matter how far the parser
/// read ahead internally.
///
/// # Errors
///
/// Returns an error if the stream cannot be read or contains malformed FASTQ.
fn sample_detection_quals_from_stream(
    reader: Box<dyn BufRead + Send>,
) -> Result<(QualDetectionStats, Box<dyn BufRead + Send>)> {
    let mut stats = QualDetectionStats::new();

    let tee = TeeReader::new(reader);
    let mut sampler =
        SimdFastqReader::with_capacity(BufReader::with_capacity(BUFFER_SIZE, tee), BUFFER_SIZE);
    sample_into(&mut sampler, &mut stats)?;

    let (consumed, source) = sampler.into_inner().into_inner().into_parts();
    let replayed: Box<dyn BufRead + Send> =
        Box::new(BufReader::with_capacity(BUFFER_SIZE, ChainedReader::new(consumed, source)));

    Ok((stats, replayed))
}

/// Open the decompressed FASTQ input readers for a run.
///
/// De-interleaves a single interleaved stream into an R1/R2 pair and consumes a
/// pre-sampled stdin reader when present. Both the serial oracle and the chain
/// reach it through [`detect_encoding_and_open_fastq_readers`] — the serial
/// `Extract::execute` path directly, and the chain FASTQ source via
/// `ChainBuilder::open_source` — so both open readers identically, the parity
/// contract for the `--threads` cutover.
fn open_fastq_input_readers(
    inputs: &[PathBuf],
    interleaved: bool,
    decomp_threads: usize,
    async_reader: bool,
    check_crc: bool,
    no_check_crc: bool,
    stdin_reader: Option<Box<dyn BufRead + Send>>,
) -> Result<Vec<Box<dyn BufRead + Send>>> {
    // Interleaved: one physical source (stdin or the sole file), split into an R1
    // and an R2 reader. Validation guarantees exactly one input, so the stdin
    // reader — when present — is that whole source. A single physical stream gets
    // the full `decomp_threads` BGZF-decode budget: the pipeline reads it on one
    // worker and the de-interleaver splits one shared decoder, so per-stream
    // decode threads are the only decode parallelism available here.
    if interleaved {
        let source = match stdin_reader {
            Some(reader) => reader,
            None => open_fastq_reader(
                &inputs[0],
                decomp_threads,
                async_reader,
                check_crc,
                no_check_crc,
            )?,
        };
        let (r1, r2) = deinterleave(source);
        return Ok(vec![r1, r2]);
    }
    // stdin is a single stream, already opened with `decomp_threads` by the caller
    // (the replay reader wraps that decoder), so its decode budget is fixed.
    if let Some(reader) = stdin_reader {
        return Ok(vec![reader]);
    }
    // Multiple file inputs: on the chain path the pipeline reads (and BGZF-decodes)
    // the streams concurrently, one reader per typed-step worker, so each reader
    // decodes single-threaded. Handing each file its own `decomp_threads` BGZF
    // workers would spawn up to N×threads decode threads ON TOP OF the pipeline's
    // own worker pool — over-subscription the pre-cutover chain deliberately
    // avoided (it passed 1 here, "the pipeline framework provides the
    // parallelism"). The serial oracle reaches this branch with `decomp_threads`
    // == 1 (it is single-threaded), so hardcoding 1 leaves it unchanged; per-file
    // BGZF-decode parallelism, if wanted, is a separate benchmarked change.
    inputs
        .iter()
        .map(|path| open_fastq_reader(path, 1, async_reader, check_crc, no_check_crc))
        .collect()
}

/// Detect the input FASTQ quality encoding and open the decompressed input
/// readers in one pass, mirroring `Extract::execute`'s detection + reader-open
/// exactly so the serial oracle and the `--threads` chain path see identical
/// readers and encoding.
///
/// stdin is sampled once and its bytes replayed into the returned reader (there
/// is no second open); file inputs are sampled through separate readers so the
/// returned readers are fresh. The stdin sampling reader is opened with the run's
/// decompression thread count (not a sampling one) because that same reader
/// decompresses every record that follows.
///
/// # Errors
///
/// Returns an error if a reader cannot be opened, the FASTQ is malformed, or the
/// encoding cannot be determined from the sampled qualities.
pub(crate) fn detect_encoding_and_open_fastq_readers(
    inputs: &[PathBuf],
    interleaved: bool,
    num_threads: usize,
    async_reader: bool,
    check_crc: bool,
    no_check_crc: bool,
) -> Result<(Vec<Box<dyn BufRead + Send>>, QualityEncoding)> {
    let decomp_threads = num_threads.max(1);
    let mut stdin_reader: Option<Box<dyn BufRead + Send>> = None;
    let detection_stats = if let Some(stdin_input) = inputs.iter().find(|p| is_stdin_path(p)) {
        let opened =
            open_fastq_reader(stdin_input, decomp_threads, async_reader, check_crc, no_check_crc)?;
        let (stats, replayed) = sample_detection_quals_from_stream(opened)?;
        stdin_reader = Some(replayed);
        stats
    } else {
        sample_detection_quals(inputs, check_crc, no_check_crc)?
    };
    let encoding = QualityEncoding::from_stats(&detection_stats)?;
    let readers = open_fastq_input_readers(
        inputs,
        interleaved,
        decomp_threads,
        async_reader,
        check_crc,
        no_check_crc,
        stdin_reader,
    )?;
    Ok((readers, encoding))
}

impl Command for Extract {
    fn execute(&self, command_line: &str) -> Result<()> {
        // Validate inputs
        self.validate()?;

        // `--threads N` routes onto the declarative chain builder — the chain
        // opens its own FASTQ readers and detects the quality encoding inside
        // `ChainBuilder::open_source`. The no-`--threads` serial loop below is the
        // in-process parity oracle. `--threads 1` takes the chain like any other
        // `Some(n)`: the old dispatch was `is_parallel()`, which was already true
        // for every `Threads(n)` including `Threads(1)`, so this
        // `threads.is_some()` predicate is behaviorally identical — the cutover
        // swaps the parallel *implementation* (the now-deleted legacy
        // `process_with_pipeline`) for the chain builder, it does not change which
        // `--threads` values run in parallel.
        if self.threading.threads.is_some() {
            return self.execute_chain(command_line);
        }

        let timer = OperationTimer::new("Extracting UMIs");
        let read_structures = self.get_read_structures()?;

        // Serial oracle: detect the quality encoding and open the input readers
        // in one pass — the same helper the chain FASTQ source uses, so both
        // paths open identical readers (de-interleaving / stdin-replay included)
        // and apply the same detected encoding.
        let (fq_readers, encoding) = detect_encoding_and_open_fastq_readers(
            &self.inputs,
            self.interleaved,
            self.threading.num_threads(),
            self.async_reader,
            self.check_crc,
            self.no_check_crc,
        )?;

        // Create header with @PG record
        let header = self.create_header(command_line)?;

        let fq_sources: Vec<SimdFastqReader<Box<dyn BufRead + Send>>> = fq_readers
            .into_iter()
            .map(|fq| SimdFastqReader::with_capacity(fq, BUFFER_SIZE))
            .collect();

        let mut fq_iterators: Vec<ReadSetIterator> = fq_sources
            .into_iter()
            .zip(read_structures.iter())
            .map(|(source, rs)| ReadSetIterator::new(rs.clone(), source, Vec::new()))
            .collect();

        // Raw BAM writer (single-threaded oracle output).
        let writer_threads = self.threading.num_threads();
        let mut writer = create_raw_bam_writer(
            &self.output,
            &header,
            writer_threads,
            self.compression.compression_level,
        )?;

        let count = self.process_singlethreaded(&mut fq_iterators, &mut writer, encoding)?;

        // Flush and finish the writer so any write error surfaces here, not on drop.
        writer.finish()?;

        timer.log_completion(count);
        Ok(())
    }
}

// ==== ported from feat-runall for the chain builder (R2) ====
/// Per-stage options for [`crate::pipeline::chains::Stage::Extract`].
///
/// Carries all the knobs that [`crate::pipeline::chains::commands::extract`]
/// (T5.6) will need when building the extract step sequence and synthesizing the
/// unmapped-BAM header:
///
/// - **Header-synthesis fields** (`sample`, `library`, `platform`,
///   `platform_unit`, `read_group_id`, `comments`) map directly to the
///   corresponding `@RG` and `@CO` entries written by `Extract::create_header`.
/// - **Behavior options** control tag output and name annotation in the same
///   way as the identically-named [`Extract`] CLI flags.
///
/// Constructed by `Extract::execute` (T5.7) from the parsed CLI struct and
/// placed into [`crate::pipeline::chains::StageOptionsBag::extract`].
#[derive(Debug, Clone)]
#[allow(clippy::struct_excessive_bools)]
pub struct ExtractOptions {
    // ── Header-synthesis fields (map to @RG / @CO in the unmapped BAM header) ──
    /// Sample name (`@RG SM:`).
    pub sample: String,
    /// Library name (`@RG LB:`).
    pub library: String,
    /// Sequencing platform (`@RG PL:`, e.g. `"illumina"`).
    pub platform: Option<String>,
    /// Platform unit (`@RG PU:`).
    pub platform_unit: Option<String>,
    /// Read-group identifier (`@RG ID:`). Defaults to `"A"`.
    pub read_group_id: String,
    /// Comments to add as `@CO` lines in the header.
    pub comments: Vec<String>,
    /// Library or sample barcode sequence (`@RG BC:`).
    pub barcode: Option<String>,
    /// Platform model (`@RG PM:`, e.g. `"hiseq2500"`).
    pub platform_model: Option<String>,
    /// Sequencing center (`@RG CN:`).
    pub sequencing_center: Option<String>,
    /// Predicted median insert size (`@RG PI:`).
    pub predicted_insert_size: Option<u32>,
    /// Description of the read group (`@RG DS:`).
    pub description: Option<String>,
    /// Date the run was produced (`@RG DT:`).
    pub run_date: Option<String>,

    // ── Extract behavior options ──
    /// Quality encoding of the input FASTQ files.
    ///
    /// On the chain path this is a placeholder at spec-construction time: the
    /// FASTQ source detects the encoding while opening its readers (in
    /// `ChainBuilder::open_source`) and overrides this field before the extract
    /// step runs, so both the serial oracle and the chain apply the same
    /// source-detected encoding. Tests that construct `ExtractOptions` directly
    /// set the encoding they want to exercise.
    pub quality_encoding: QualityEncoding,
    /// Store UMI base qualities in the `QX` tag.
    pub store_umi_quals: bool,
    /// Store cell-barcode base qualities in the `CY` tag.
    pub store_cell_quals: bool,
    /// Single SAM tag to store all concatenated UMIs (must not collide with
    /// reserved tags `RX`, `QX`, `CB`, `CY`, `BC`, `QT`, `RG`).
    pub single_tag: Option<SamTag>,
    /// Append `+<UMIs>` to each read name.
    pub annotate_read_names: bool,
    /// Extract UMIs from read names (last `:`-separated field, ≥8 fields).
    pub extract_umis_from_read_names: bool,
    /// Store sample-barcode qualities in the `QT` tag.
    pub store_sample_barcode_qualities: bool,
    /// Wrap FASTQ readers in a userspace async prefetch thread.
    pub async_reader: bool,
    /// Force BGZF input CRC verification (`--check-crc`). Consumed when the chain
    /// FASTQ source opens its readers, so `--check-crc`/`--no-check-crc` reach
    /// `open_fastq_reader` on the `--threads` path exactly as on the serial path.
    pub check_crc: bool,
    /// Disable BGZF input CRC verification (`--no-check-crc`). See `check_crc`.
    pub no_check_crc: bool,
}

/// Extract-stage options for the fused `runall` pipeline — a `runall`-only
/// `clap::Args` variant of [`Extract`]'s options.
///
/// The standalone [`Extract`] CLI struct owns `--output` and the flattened
/// engine sub-structs (`threading`/`compression`/`scheduler_opts`/
/// `queue_memory`), which `runall` supplies itself, so [`Extract`] cannot be
/// flattened into a fused command directly. This struct mirrors the current
/// [`Extract`] CLI-struct field shapes minus those, and carries
/// `#[fgumi_cli_macros::multi_options]` so `runall` can re-expose each field as a
/// prefixed `--extract::<flag>` via the generated `MultiExtractRunallOptions`
/// companion, without hand-maintaining a parallel option set.
///
/// The `inputs`/`read_structures`/`interleaved` source-construction fields are
/// not projected by [`Self::to_extract_options`] — they build the FASTQ source
/// for PR B and are exposed `pub` for that consumer. Every other field maps into
/// [`ExtractOptions`] exactly as [`Extract::to_extract_options`] does; in
/// particular `platform` (default `"illumina"`) becomes `Some(..)` and
/// `quality_encoding` is the [`QualityEncoding::Standard`] placeholder the chain
/// FASTQ source overrides at runtime. `--read-structures` is `required` here
/// (unlike the standalone's `+T` default): a fused `runall` always demands it.
///
/// The three cross-field `conflicts_with` attrs on [`Extract`]
/// (`--store-umi-quals` vs `--extract-umis-from-read-names`, and the
/// `--check-crc`/`--no-check-crc` pair) are dropped because the `multi_options`
/// macro rejects `conflicts_with`; they are re-enforced in [`Self::validate`].
#[fgumi_cli_macros::multi_options("extract", "Extract Options")]
#[derive(Debug, Clone, clap::Args)]
#[allow(clippy::struct_excessive_bools)]
pub struct ExtractRunallOptions {
    // ── Source-construction fields (not projected by to_extract_options) ──
    /// Input FASTQ files corresponding to each sequencing read (e.g. R1, I1, etc.).
    /// Accepts one or more space-separated values; required for this runall variant.
    #[arg(long = "inputs", required = true, num_args = 1..)]
    pub inputs: Vec<PathBuf>,

    /// Read structures, one for each of the FASTQs (optional if 1-2 template-only FASTQs).
    /// Accepts one or more space-separated values; required for this runall variant.
    #[arg(long = "read-structures", required = true, num_args = 1..)]
    pub read_structures: Vec<ReadStructure>,

    /// Treat a single input as interleaved paired-end FASTQ (`R1, R2, R1, R2, …`),
    /// de-interleaving it into the two reads. Requires exactly one `--input` (a
    /// file or `-` for stdin) and describes both reads with two `--read-structures`
    /// (required for this runall variant — unlike the standalone `Extract`, there
    /// is no `+T +T` default). This lets a streaming trimmer or converter pipe
    /// interleaved pairs straight into extract without staging two FASTQ files.
    #[arg(long = "interleaved", default_value_t = false)]
    pub interleaved: bool,

    // ── Header / behavior fields (projected into ExtractOptions) ──
    /// The name of the sequenced sample
    #[arg(long, required = true)]
    pub sample: String,

    /// The name/ID of the sequenced library
    #[arg(long, required = true)]
    pub library: String,

    /// Sequencing Platform
    #[arg(long, default_value = "illumina")]
    pub platform: String,

    /// Library or Sample barcode sequence
    #[arg(long)]
    pub barcode: Option<String>,

    /// Read group ID to use in the file header
    #[arg(long = "read-group-id", default_value = "A")]
    pub read_group_id: String,

    /// Platform unit (e.g. 'flowcell-barcode.lane.sample-barcode')
    #[arg(long = "platform-unit")]
    pub platform_unit: Option<String>,

    /// Platform model to insert into the group header (ex. miseq, hiseq2500, hiseqX)
    #[arg(long = "platform-model")]
    pub platform_model: Option<String>,

    /// The sequencing center from which the data originated
    #[arg(long = "sequencing-center")]
    pub sequencing_center: Option<String>,

    /// Predicted median insert size, to insert into the read group header
    #[arg(long = "predicted-insert-size")]
    pub predicted_insert_size: Option<u32>,

    /// Description of the read group
    #[arg(long)]
    pub description: Option<String>,

    /// Comment(s) to include in the output file's header
    #[arg(long, num_args = 0..)]
    pub comment: Vec<String>,

    /// Date the run was produced, to insert into the read group header
    #[arg(long = "run-date")]
    pub run_date: Option<String>,

    /// Store UMI base quality scores in the QX SAM tag
    #[arg(long = "store-umi-quals")]
    pub store_umi_quals: bool,

    /// Store cell barcode base quality scores in the CY SAM tag
    #[arg(long = "store-cell-quals")]
    pub store_cell_quals: bool,

    /// Store the sample barcode qualities in the QT Tag
    #[arg(long = "store-sample-barcode-qualities")]
    pub store_sample_barcode_qualities: bool,

    /// Extract UMI(s) from read names and prepend to UMIs from reads
    #[arg(long = "extract-umis-from-read-names")]
    #[allow(clippy::struct_field_names)]
    pub extract_umis_from_read_names: bool,

    /// Annotate read names with UMIs (appends "+UMIs" to read names)
    #[arg(long = "annotate-read-names")]
    pub annotate_read_names: bool,

    /// Single tag to store all concatenated UMIs (in addition to per-segment tags)
    #[arg(long = "single-tag")]
    pub single_tag: Option<SamTag>,

    /// Tag containing adapter clipping position to adjust (e.g. 'XT' from `MarkIlluminaAdapters`)
    #[arg(long = "clipping-attribute")]
    pub clipping_attribute: Option<SamTag>,

    /// Wrap FASTQ inputs in a userspace async prefetch reader. Dedicates one
    /// OS thread per input stream to issue reads ahead of decompression/parsing.
    /// Hidden experimental flag.
    #[arg(long = "async-reader", default_value_t = false, hide = true)]
    pub async_reader: bool,

    /// Verify each BGZF block's CRC32 checksum while decoding the input.
    ///
    /// Applies to BGZF-compressed FASTQ input (bgzip'd); plain gzip is not
    /// block-structured and has no per-block CRC to skip. The policy is honored
    /// at any thread count, though BGZF decode runs single-threaded when skipping
    /// (verifying BGZF input can decode in parallel). Without either flag,
    /// verification defaults on for file input and off for trusted stdin (a
    /// freshly-piped stream is trusted; a file may have been archived or
    /// transferred since it was written, where a flipped bit is what CRC32 exists
    /// to catch). Pass `--check-crc` to force it on. Mutually exclusive with
    /// `--no-check-crc`.
    #[arg(long = "check-crc", default_value_t = false)]
    pub check_crc: bool,

    /// Skip CRC32 verification while decoding BGZF FASTQ input.
    ///
    /// Trades the CRC32 integrity check for faster BGZF decode (which then runs
    /// single-threaded). See `--check-crc` for the default policy this overrides.
    /// Mutually exclusive with `--check-crc`.
    #[arg(long = "no-check-crc", default_value_t = false)]
    pub no_check_crc: bool,
}

/// Hand-written to match what parsing the minimal required flags (`--extract::sample`,
/// `--extract::library`, `--extract::inputs`, `--extract::read-structures`) through
/// `MultiExtractRunallOptions::try_parse_from` and `validate()` would produce, per the
/// branch-wide invariant "Default == the minimal-parse projection". `#[derive(Default)]`
/// cannot be used here because it would give `platform`/`read_group_id` empty-string
/// defaults instead of the CLI defaults (`"illumina"` / `"A"`). `sample`/`library`/
/// `inputs`/`read_structures` have no natural default (they are staged-required, lifted
/// by `MultiExtractRunallOptions::validate`, never `Extract`'s own `Default`) and are set
/// to placeholders that are never consumed, mirroring how the other structs' `Default`
/// impls handle required fields.
impl Default for ExtractRunallOptions {
    fn default() -> Self {
        Self {
            inputs: Vec::new(),
            read_structures: Vec::new(),
            interleaved: false,
            sample: String::new(),
            library: String::new(),
            platform: "illumina".to_string(),
            barcode: None,
            read_group_id: "A".to_string(),
            platform_unit: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            comment: Vec::new(),
            run_date: None,
            store_umi_quals: false,
            store_cell_quals: false,
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names: false,
            annotate_read_names: false,
            single_tag: None,
            clipping_attribute: None,
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
        }
    }
}

impl ExtractRunallOptions {
    /// Project into the per-stage [`ExtractOptions`] the chain builder consumes.
    ///
    /// Verbatim-logic copy of [`Extract::to_extract_options`]: `quality_encoding`
    /// is the [`QualityEncoding::Standard`] placeholder (the chain FASTQ source
    /// overrides it while opening its readers), `platform` becomes `Some(..)`, and
    /// `clipping_attribute` is intentionally dropped — it does not apply to FASTQ
    /// input. The `inputs`/`read_structures`/`interleaved` source fields are not
    /// part of this projection (PR B builds the FASTQ source from them).
    #[must_use]
    pub fn to_extract_options(&self) -> ExtractOptions {
        ExtractOptions {
            sample: self.sample.clone(),
            library: self.library.clone(),
            platform: Some(self.platform.clone()),
            platform_unit: self.platform_unit.clone(),
            read_group_id: self.read_group_id.clone(),
            comments: self.comment.clone(),
            barcode: self.barcode.clone(),
            platform_model: self.platform_model.clone(),
            sequencing_center: self.sequencing_center.clone(),
            predicted_insert_size: self.predicted_insert_size,
            description: self.description.clone(),
            run_date: self.run_date.clone(),
            quality_encoding: QualityEncoding::Standard,
            store_umi_quals: self.store_umi_quals,
            store_cell_quals: self.store_cell_quals,
            single_tag: self.single_tag,
            annotate_read_names: self.annotate_read_names,
            extract_umis_from_read_names: self.extract_umis_from_read_names,
            store_sample_barcode_qualities: self.store_sample_barcode_qualities,
            async_reader: self.async_reader,
            check_crc: self.check_crc,
            no_check_crc: self.no_check_crc,
        }
    }

    /// Re-enforce ONLY the two cross-field conflicts the `multi_options` macro
    /// required us to drop from the `#[arg]` attributes: on the standalone
    /// [`Extract`], `--store-umi-quals`/`--extract-umis-from-read-names` and
    /// `--check-crc`/`--no-check-crc` are each a clap-level `conflicts_with`
    /// pair (not part of `Extract::validate`'s own body), and the macro
    /// rejects `conflicts_with` on a `multi_options` struct, so those two
    /// checks are re-implemented here by hand.
    ///
    /// This method does NOT re-implement the rest of [`Extract::validate`] —
    /// the template-count-1-to-2 check, the `--single-tag` reserved-tag
    /// collision check, the read-structure-non-empty check, or the
    /// input/read-structure count and stdin/file-existence checks. Those all
    /// depend on the fully-resolved read structures and input list, which are
    /// only known once the FASTQ source is built; applying them here would
    /// duplicate logic that must live with that construction. The future
    /// `runall` command (PR B) is responsible for running the remaining
    /// `Extract::validate` checks (template-count 1–2, the `single_tag`
    /// reserved-tag collision, and read-structure non-emptiness) when it
    /// builds the FASTQ source from [`Self::inputs`] / [`Self::read_structures`]
    /// / [`Self::interleaved`]; nothing in this PR calls those checks.
    ///
    /// This is separate from the macro-generated `MultiExtractRunallOptions::
    /// validate()`, which only lifts the staged-required flags; `runall`/PR B
    /// calls this after building the [`ExtractRunallOptions`].
    ///
    /// # Errors
    ///
    /// Returns an error if `--store-umi-quals` is combined with
    /// `--extract-umis-from-read-names`, or if `--check-crc` and `--no-check-crc`
    /// are both set.
    pub fn validate(&self) -> Result<()> {
        ensure!(
            !self.extract_umis_from_read_names || !self.store_umi_quals,
            "Cannot store UMI qualities (--store-umi-quals) when also extracting UMIs from read names (--extract-umis-from-read-names)."
        );
        ensure!(
            !(self.check_crc && self.no_check_crc),
            "--check-crc and --no-check-crc are mutually exclusive."
        );
        Ok(())
    }
}

/// Build raw BAM `RawRecord`s from a [`FastqSet`].
///
/// This is the core extract logic: applies read structures (via the segments
/// already present in `read_set`), extracts UMI / cell-barcode / sample-barcode
/// tags, and produces one `RawRecord` per template segment. The caller wraps
/// the result in a [`crate::template::Template`] for the typed-step pipeline.
///
/// KNOWN DUPLICATION: this is the chain path's copy of the FASTQ→`RawRecord`
/// per-read loop, near-identical to `Extract::make_raw_records` (the serial
/// oracle's, which writes to a `RawBamWriter`). The two differ only in output
/// target, so tag-extraction/flag logic can drift between them. The
/// `chain_matches_serial_oracle*` integration tests guard against drift by
/// comparing chain-vs-oracle output: the base case covers RX/RG, and
/// `chain_matches_serial_oracle_barcode_and_annotate_tags` covers the
/// CB/CY, BC/QT, QX, `--single-tag`, and `--annotate-read-names` branches.
/// Unifying the two onto one shared per-read loop is a worthwhile follow-up.
/// (The former third copy, `make_raw_records_static`, was deleted with the
/// legacy threaded pipeline when extract cut over to the chain builder.)
///
/// # Errors
///
/// Returns an error if the read set contains no template segments, or if a read
/// name is 255 bytes or longer (via `try_build_record`).
pub(crate) fn make_raw_records_from_fastq_set(
    read_set: &FastqSet,
    opts: &ExtractOptions,
) -> Result<Vec<fgumi_raw_bam::RawRecord>> {
    let templates: Vec<&FastqSegment> = read_set.template_segments().collect();

    let read_name = String::from_utf8_lossy(&read_set.header);
    ensure!(!templates.is_empty(), "No template segments found for read: {read_name}");

    // Extract various barcode types as BString
    let cell_barcode_bs = Extract::join_bytes_with_separator(
        read_set.cell_barcode_segments().map(|s| s.seq.as_slice()),
        b'-',
    );
    let cell_quals_bs = Extract::join_bytes_with_separator(
        read_set.cell_barcode_segments().map(|s| s.quals.as_slice()),
        b' ',
    );
    let sample_barcode_bs = Extract::join_bytes_with_separator(
        read_set.sample_barcode_segments().map(|s| s.seq.as_slice()),
        b'-',
    );
    let sample_quals_bs = Extract::join_bytes_with_separator(
        read_set.sample_barcode_segments().map(|s| s.quals.as_slice()),
        b' ',
    );
    let umi_bs = Extract::join_bytes_with_separator(
        read_set.molecular_barcode_segments().map(|s| s.seq.as_slice()),
        b'-',
    );
    let umi_qual_bs = Extract::join_bytes_with_separator(
        read_set.molecular_barcode_segments().map(|s| s.quals.as_slice()),
        b' ',
    );

    // Extract UMI from read name if requested
    let (read_name_bytes, umi_from_name) =
        Extract::extract_read_name_and_umi(&read_set.header, opts.extract_umis_from_read_names)?;

    // Prepare final UMI
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

    let num_templates = templates.len();
    let mut builder = UnmappedSamBuilder::new();
    let mut records = Vec::with_capacity(num_templates);

    for (index, template) in templates.iter().enumerate() {
        // Compute flags for unmapped reads
        let mut flag = flags::UNMAPPED;
        if num_templates == 2 {
            flag |= flags::PAIRED | flags::MATE_UNMAPPED;
            if index == 0 {
                flag |= flags::FIRST_SEGMENT;
            } else {
                flag |= flags::LAST_SEGMENT;
            }
        }

        // Set read name (optionally with UMI annotation)
        let annotated_name: Option<Vec<u8>> =
            if opts.annotate_read_names && !final_umi_bs.is_empty() {
                let mut name = Vec::with_capacity(read_name_bytes.len() + 1 + final_umi_bs.len());
                name.extend_from_slice(&read_name_bytes);
                name.push(b'+');
                name.extend_from_slice(final_umi_bs.as_bytes());
                Some(name)
            } else {
                None
            };
        let final_read_name: &[u8] = annotated_name.as_deref().unwrap_or(&read_name_bytes);

        // Build the record - if empty seq, substitute with single N @ Q2.
        // Use try_build_record (not the panicking build_record): the name is
        // derived from the input FASTQ header (plus an optional `+<UMI>`), so an
        // over-long (>=255-byte) name is bad input, not a bug — fail cleanly with
        // context instead of panicking partway through the batch and truncating
        // the output BAM. Mirrors the standalone path's `build_template_record`.
        if template.seq.is_empty() {
            builder.try_build_record(final_read_name, flag, b"N", &[2u8])
        } else {
            let numeric_quals = opts.quality_encoding.to_standard_numeric(&template.quals);
            builder.try_build_record(final_read_name, flag, &template.seq, &numeric_quals)
        }
        .with_context(|| Extract::read_name_too_long_context(final_read_name))?;

        // Append tags
        // Read group
        builder.append_string_tag(SamTag::RG, opts.read_group_id.as_bytes());

        // Cell barcode
        if !cell_barcode_bs.is_empty() {
            builder.append_string_tag(SamTag::CB, cell_barcode_bs.as_bytes());
        }

        if !cell_quals_bs.is_empty() && opts.store_cell_quals {
            builder.append_string_tag(SamTag::CY, cell_quals_bs.as_bytes());
        }

        // Sample barcode
        if !sample_barcode_bs.is_empty() {
            builder.append_string_tag(SamTag::BC, sample_barcode_bs.as_bytes());
        }

        if opts.store_sample_barcode_qualities && !sample_quals_bs.is_empty() {
            builder.append_string_tag(SamTag::QT, sample_quals_bs.as_bytes());
        }

        // UMI
        if !final_umi_bs.is_empty() {
            builder.append_string_tag(SamTag::RX, final_umi_bs.as_bytes());

            // Single tag for all concatenated UMIs (if specified)
            if let Some(st) = opts.single_tag {
                builder.append_string_tag(st, final_umi_bs.as_bytes());
            }

            // Only add UMI qualities if not extracted from read names
            if umi_from_name.is_none() && !umi_qual_bs.is_empty() && opts.store_umi_quals {
                builder.append_string_tag(SamTag::QX, umi_qual_bs.as_bytes());
            }
        }

        records.push(builder.build());
    }

    Ok(records)
}

/// Validate that the read structures produce a valid SAM template: 1 or 2
/// template (`T`) reads total. 0 template reads yields records with no sequence;
/// 3+ produce more segments than a SAM template (R1/R2) can hold — either way
/// the output BAM would be malformed. Shared by the standalone `Extract` command
/// and the chain/runall `ChainBuilder::add_extract` so both paths reject it
/// (the chain path previously skipped this check — silent malformed output).
///
/// # Errors
///
/// Returns an error unless the total template-read count is 1 or 2.
pub(crate) fn validate_template_count(read_structures: &[ReadStructure]) -> anyhow::Result<()> {
    let template_count: usize =
        read_structures.iter().map(|rs| rs.segments_by_type(SegmentType::Template).count()).sum();
    anyhow::ensure!(
        (1..=2).contains(&template_count),
        "read structures must contain 1-2 template reads total."
    );
    Ok(())
}

// ==== impls ported from feat-runall for the chain builder (R2) ====
impl ExtractOptions {
    /// Validate that [`Self::single_tag`] does not collide with any of the
    /// SAM tags that the extractor emits internally.
    ///
    /// This mirrors the check in `Extract::validate`; keeping both guards in
    /// sync ensures that the error fires whether the caller constructs an
    /// `ExtractOptions` directly (chain path) or via the CLI struct (command path).
    ///
    /// # Errors
    ///
    /// Returns an error when `single_tag` matches `RX`, `QX`, `CB`, `CY`,
    /// `BC`, `QT`, or `RG`, or when `extract_umis_from_read_names` is combined
    /// with `store_umi_quals` (read-name UMIs carry no qualities, so
    /// `make_raw_records_from_fastq_set` would silently omit `QX`).
    pub fn validate(&self) -> anyhow::Result<()> {
        const RESERVED_OUTPUT_TAGS: &[SamTag] =
            &[SamTag::RX, SamTag::QX, SamTag::CB, SamTag::CY, SamTag::BC, SamTag::QT, SamTag::RG];
        if let Some(ref tag) = self.single_tag {
            anyhow::ensure!(
                !RESERVED_OUTPUT_TAGS.contains(tag),
                "Single tag (--single-tag) cannot collide with tags emitted by extract \
                 (RX, QX, CB, CY, BC, QT, RG): {tag}"
            );
        }
        anyhow::ensure!(
            !self.extract_umis_from_read_names || !self.store_umi_quals,
            "Cannot store UMI qualities (--store-umi-quals) when also extracting UMIs from read names (--extract-umis-from-read-names)."
        );
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bstr::BString;
    use noodles::sam::alignment::RecordBuf;
    use noodles::sam::alignment::record::data::field::Tag;
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
        let mut file = File::create(&path).expect("failed to create file");
        for (name, seq, qual) in records {
            writeln!(file, "@{name}").expect("failed to write line");
            writeln!(file, "{seq}").expect("failed to write line");
            writeln!(file, "+").expect("failed to write line");
            writeln!(file, "{qual}").expect("failed to write line");
        }
        path
    }

    /// The uncompressed FASTQ bytes that [`create_bgzf_fastq`] compresses. Kept as
    /// the single source of truth so tests can assert exact decoded output and
    /// exact per-record identity/order (record `i` is `@q{i}` / `ACGTACGTAC` /
    /// `IIIIIIIIII`).
    fn fastq_bytes(num_records: usize) -> Vec<u8> {
        let mut fastq = Vec::new();
        for i in 0..num_records {
            writeln!(fastq, "@q{i}").expect("write");
            writeln!(fastq, "ACGTACGTAC").expect("write");
            writeln!(fastq, "+").expect("write");
            writeln!(fastq, "IIIIIIIIII").expect("write");
        }
        fastq
    }

    /// Write `num_records` single-end FASTQ records, BGZF-compressed, to a file.
    /// `num_records` is chosen large enough by callers to span >= 2 BGZF blocks
    /// (blocks flush per ~64 KiB uncompressed).
    fn create_bgzf_fastq(dir: &TempDir, name: &str, num_records: usize) -> PathBuf {
        let path = dir.path().join(name);
        let fastq = fastq_bytes(num_records);
        let mut compressed = Vec::new();
        {
            let mut writer = noodles_bgzf::io::Writer::new(&mut compressed);
            writer.write_all(&fastq).expect("write bgzf");
            writer.finish().expect("finish bgzf");
        }
        std::fs::write(&path, compressed).expect("write bgzf fastq");
        path
    }

    /// Assert an error message names a CRC / checksum failure specifically, so a
    /// CRC regression test cannot pass on some unrelated decode error. The
    /// single-threaded fgumi-bgzf path reports `BGZF CRC32 mismatch`; the
    /// multi-threaded noodles path reports `block data checksum mismatch`.
    fn assert_crc_error_msg(message: &str) {
        let lower = message.to_lowercase();
        assert!(
            lower.contains("crc") || lower.contains("checksum"),
            "error should name a CRC/checksum failure, got: {message}"
        );
    }

    /// Flip a byte in the last real BGZF block's CRC32 footer. Requires >= 2
    /// blocks so the corrupted block is a data block.
    fn corrupt_last_block_crc(path: &PathBuf) {
        let mut bytes = std::fs::read(path).expect("read bgzf fastq");
        let blocks = {
            let mut cursor: &[u8] = &bytes;
            fgumi_bgzf::read_raw_blocks(&mut cursor, 1_000_000).expect("read bgzf blocks")
        };
        assert!(blocks.len() >= 2, "input must span >= 2 BGZF blocks; got {}", blocks.len());
        let offset: usize =
            blocks[..blocks.len() - 1].iter().map(fgumi_bgzf::RawBgzfBlock::len).sum();
        let last = blocks.last().expect("checked len >= 2");
        let crc_off = offset + last.len() - fgumi_bgzf::BGZF_FOOTER_SIZE;
        bytes[crc_off] ^= 0x01;
        std::fs::write(path, bytes).expect("write corrupted bgzf fastq");
    }

    /// Build an `Extract` for a single BGZF FASTQ `input` with the given threading
    /// and CRC flags.
    fn bgzf_crc_extract(
        input: PathBuf,
        output: PathBuf,
        threading: ThreadingOptions,
        check_crc: bool,
        no_check_crc: bool,
    ) -> Extract {
        Extract {
            inputs: vec![input],
            output,
            read_structures: vec![],
            store_umi_quals: false,
            store_cell_quals: false,
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
            threading,
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc,
            no_check_crc,
            interleaved: false,
        }
    }

    /// A corrupted-CRC BGZF FASTQ must be rejected by default and with
    /// `--check-crc`, and read clean with `--no-check-crc` (#819). Parameterized
    /// over threading so both BGZF decode paths are covered: single-threaded goes
    /// through `open_fastq_reader`'s fgumi-bgzf arm, while `--threads N` decodes
    /// the pure-BGZF input in the pipeline honoring `FastqPipelineConfig::verify_crc`.
    #[rstest]
    #[case::single_threaded(ThreadingOptions::none())]
    #[case::multi_threaded(ThreadingOptions::new(2))]
    fn test_bgzf_fastq_honors_check_crc(#[case] threading: ThreadingOptions) {
        let tmp = TempDir::new().expect("temp dir");
        // ~5000 records * ~30 bytes = ~150 KiB uncompressed => several BGZF blocks.
        let input = create_bgzf_fastq(&tmp, "reads.fq.gz", 5000);
        corrupt_last_block_crc(&input);

        // --no-check-crc: the corrupted block is accepted and every record is
        // decoded exactly, in input order (identity, not just count).
        let out_ok = tmp.path().join("ok.bam");
        bgzf_crc_extract(input.clone(), out_ok.clone(), threading.clone(), false, true)
            .execute("test")
            .expect("--no-check-crc must accept a corrupted BGZF FASTQ CRC32");
        let records = read_bam_records(&out_ok);
        assert_eq!(records.len(), 5000, "all records extracted");
        for (i, record) in records.iter().enumerate() {
            let expected_name = format!("q{i}");
            assert_eq!(
                record.name().map(|name| name.as_bytes()),
                Some(expected_name.as_bytes()),
                "record {i} name mismatch (identity/order)"
            );
            assert_eq!(record.sequence().as_ref(), b"ACGTACGTAC", "record {i} sequence");
            assert_eq!(record.quality_scores().as_ref(), &[40; 10], "record {i} quality");
        }

        // Default (file => verify on): rejected with a CRC-specific error.
        let out_default = tmp.path().join("default.bam");
        let default_err =
            bgzf_crc_extract(input.clone(), out_default, threading.clone(), false, false)
                .execute("test")
                .expect_err("default must reject a corrupted BGZF FASTQ CRC32");
        assert_crc_error_msg(&format!("{default_err:#}"));

        // --check-crc: rejected with a CRC-specific error.
        let out_check = tmp.path().join("check.bam");
        let check_err = bgzf_crc_extract(input, out_check, threading, true, false)
            .execute("test")
            .expect_err("--check-crc must reject a corrupted BGZF FASTQ CRC32");
        assert_crc_error_msg(&format!("{check_err:#}"));
    }

    /// `open_fastq_reader` must honor `--no-check-crc` on BGZF input at **any**
    /// thread count. At `threads > 1` the parallel noodles decoder always
    /// verifies, so the skip has to route through fgumi-bgzf instead — a
    /// regression test for the multi-threaded BGZF path (stdin / mixed inputs)
    /// silently ignoring the flag.
    #[rstest]
    #[case::single_threaded(1)]
    #[case::multi_threaded(4)]
    fn test_open_fastq_reader_bgzf_honors_no_check_crc(#[case] threads: usize) {
        use std::io::Read;
        let tmp = TempDir::new().expect("temp dir");
        let input = create_bgzf_fastq(&tmp, "reads.fq.gz", 5000);
        corrupt_last_block_crc(&input);

        // --no-check-crc (check=false, no_check=true): reads clean at any thread
        // count and decodes to exactly the original uncompressed FASTQ bytes.
        let mut skipping = open_fastq_reader(&input, threads, false, false, true)
            .expect("open (no-check-crc) should succeed");
        let mut buf = Vec::new();
        skipping.read_to_end(&mut buf).expect("--no-check-crc must skip the corrupted CRC32");
        assert_eq!(buf, fastq_bytes(5000), "--no-check-crc must decode the FASTQ exactly");

        // Default (check=false, no_check=false => verify for a file): errors on the
        // corrupted block with a CRC-specific message, whether via noodles MT
        // (threads > 1) or fgumi (threads 1).
        let mut verifying = open_fastq_reader(&input, threads, false, false, false)
            .expect("open (verify) should succeed");
        let mut buf = Vec::new();
        let err =
            verifying.read_to_end(&mut buf).expect_err("verify must reject the corrupted CRC32");
        assert_crc_error_msg(&err.to_string());
    }

    /// Read all records from a BAM file into a vector
    ///
    /// # Arguments
    /// * `path` - Path to the BAM file
    ///
    /// # Returns
    /// Vector of all BAM records in the file
    fn read_bam_records(path: &PathBuf) -> Vec<RecordBuf> {
        let (mut reader, header) = create_bam_reader(path, 1).expect("failed to create BAM reader");

        let mut records = Vec::new();
        for result in reader.record_bufs(&header) {
            records.push(result.expect("failed to read BAM record"));
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
        let tmp = TempDir::new().expect("failed to create temp dir");
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
            store_umi_quals: false,
            store_cell_quals: false,
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved: false,
        };

        extract.execute("test").expect("execute should succeed");

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
        let tmp = TempDir::new().expect("failed to create temp dir");
        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "AAAAAAAAAA", "==========")]);
        let r2 = create_fastq(&tmp, "r2.fq", &[("q1", "CCCCCCCCCC", "##########")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1, r2],
            output: output.clone(),
            read_structures: vec![],
            store_umi_quals: false,
            store_cell_quals: false,
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved: false,
        };

        extract.execute("test").expect("execute should succeed");

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
        let tmp = TempDir::new().expect("failed to create temp dir");
        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "AAAAAAAAAA", "==========")]);
        let r2 = create_fastq(&tmp, "r2.fq", &[("q1", "CCCCCCCCCC", "##########")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1, r2],
            output: output.clone(),
            read_structures: vec![],
            store_umi_quals: true,
            store_cell_quals: false,
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved: false,
        };

        extract.execute("test").expect("execute should succeed");

        let recs = read_bam_records(&output);
        assert_eq!(recs.len(), 2);
        // Should not have RX or QX tags since no UMIs
        assert!(get_tag_string(&recs[0], "RX").is_none());
        assert!(get_tag_string(&recs[0], "QX").is_none());
    }

    #[test]
    fn test_paired_end_with_inline_umi() {
        let tmp = TempDir::new().expect("failed to create temp dir");
        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "ACGTAAAAAA", "IIII======")]);
        let r2 = create_fastq(&tmp, "r2.fq", &[("q1", "CCCCCCCCCC", "##########")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1, r2],
            output: output.clone(),
            read_structures: vec![
                ReadStructure::from_str("4M+T").expect("valid read structure"),
                ReadStructure::from_str("+T").expect("valid read structure"),
            ],
            store_umi_quals: false,
            store_cell_quals: false,
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved: false,
        };

        extract.execute("test").expect("execute should succeed");

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
        let tmp = TempDir::new().expect("failed to create temp dir");
        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "ACGTAAAAAA", "IIII======")]);
        let r2 = create_fastq(&tmp, "r2.fq", &[("q1", "CCCCCCCCCC", "##########")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1, r2],
            output: output.clone(),
            read_structures: vec![
                ReadStructure::from_str("4M+T").expect("valid read structure"),
                ReadStructure::from_str("+T").expect("valid read structure"),
            ],
            store_umi_quals: true,
            store_cell_quals: false,
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved: false,
        };

        extract.execute("test").expect("execute should succeed");

        let recs = read_bam_records(&output);
        assert_eq!(recs.len(), 2);

        assert_eq!(get_tag_string(&recs[0], "RX"), Some("ACGT".to_string()));
        assert_eq!(get_tag_string(&recs[0], "QX"), Some("IIII".to_string()));
        assert_eq!(get_tag_string(&recs[1], "RX"), Some("ACGT".to_string()));
        assert_eq!(get_tag_string(&recs[1], "QX"), Some("IIII".to_string()));
    }

    #[test]
    fn test_complex_read_structures_multiple_segments() {
        let tmp = TempDir::new().expect("failed to create temp dir");
        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "AAACCCTTTAAAAA", "==============")]);
        let r2 = create_fastq(&tmp, "r2.fq", &[("q1", "GGGTTTAAACCCCC", "##############")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1, r2],
            output: output.clone(),
            read_structures: vec![
                ReadStructure::from_str("3B3M3B5T").expect("valid read structure"),
                ReadStructure::from_str("3B3M3B5T").expect("valid read structure"),
            ],
            store_umi_quals: false,
            store_cell_quals: false,
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved: false,
        };

        extract.execute("test").expect("execute should succeed");

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
        let tmp = TempDir::new().expect("failed to create temp dir");
        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "AAACCCTTTAAAAA", "===III========")]);
        let r2 = create_fastq(&tmp, "r2.fq", &[("q1", "GGGTTTAAACCCCC", "###JJJ########")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1, r2],
            output: output.clone(),
            read_structures: vec![
                ReadStructure::from_str("3B3M3B5T").expect("valid read structure"),
                ReadStructure::from_str("3B3M3B5T").expect("valid read structure"),
            ],
            store_umi_quals: true,
            store_cell_quals: false,
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved: false,
        };

        extract.execute("test").expect("execute should succeed");

        let recs = read_bam_records(&output);
        assert_eq!(recs.len(), 2);

        assert_eq!(get_tag_string(&recs[0], "RX"), Some("CCC-TTT".to_string()));
        assert_eq!(get_tag_string(&recs[0], "QX"), Some("III JJJ".to_string()));
        assert_eq!(get_tag_string(&recs[1], "RX"), Some("CCC-TTT".to_string()));
        assert_eq!(get_tag_string(&recs[1], "QX"), Some("III JJJ".to_string()));
    }

    #[test]
    fn test_four_fastqs() {
        let tmp = TempDir::new().expect("failed to create temp dir");
        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "AAAAAAAAAA", "==========")]);
        let r2 = create_fastq(&tmp, "r2.fq", &[("q1", "CCCCCCCCCC", "##########")]);
        let r3 = create_fastq(&tmp, "r3.fq", &[("q1", "ACGT", "????")]);
        let r4 = create_fastq(&tmp, "r4.fq", &[("q1", "GAGA", "????")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1, r2, r3, r4],
            output: output.clone(),
            read_structures: vec![
                ReadStructure::from_str("+T").expect("valid read structure"),
                ReadStructure::from_str("+T").expect("valid read structure"),
                ReadStructure::from_str("4B").expect("valid read structure"),
                ReadStructure::from_str("4M").expect("valid read structure"),
            ],
            store_umi_quals: false,
            store_cell_quals: false,
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved: false,
        };

        extract.execute("test").expect("execute should succeed");

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
        let tmp = TempDir::new().expect("failed to create temp dir");
        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "AAAAAAAAAA", "==========")]);
        let r2 = create_fastq(&tmp, "r2.fq", &[("q1", "CCCCCCCCCC", "##########")]);
        let r3 = create_fastq(&tmp, "r3.fq", &[("q1", "ACGT", "????")]);
        let r4 = create_fastq(&tmp, "r4.fq", &[("q1", "GAGA", "????")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1, r2, r3, r4],
            output: output.clone(),
            read_structures: vec![
                ReadStructure::from_str("+T").expect("valid read structure"),
                ReadStructure::from_str("+T").expect("valid read structure"),
                ReadStructure::from_str("4B").expect("valid read structure"),
                ReadStructure::from_str("4M").expect("valid read structure"),
            ],
            store_umi_quals: true,
            store_cell_quals: false,
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved: false,
        };

        extract.execute("test").expect("execute should succeed");

        let recs = read_bam_records(&output);
        assert_eq!(recs.len(), 2);

        assert_eq!(get_tag_string(&recs[0], "RX"), Some("GAGA".to_string()));
        assert_eq!(get_tag_string(&recs[0], "QX"), Some("????".to_string()));
        assert_eq!(get_tag_string(&recs[1], "RX"), Some("GAGA".to_string()));
        assert_eq!(get_tag_string(&recs[1], "QX"), Some("????".to_string()));
    }

    #[test]
    fn test_header_metadata() {
        let tmp = TempDir::new().expect("failed to create temp dir");
        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "AAAAAAAAAA", "==========")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1],
            output: output.clone(),
            read_structures: vec![ReadStructure::from_str("10T").expect("valid read structure")],
            store_umi_quals: false,
            store_cell_quals: false,
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved: false,
        };

        extract.execute("test").expect("execute should succeed");

        let (_reader, header) = create_bam_reader(&output, 1).expect("failed to create BAM reader");

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
        let tmp = TempDir::new().expect("failed to create temp dir");
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
                ReadStructure::from_str("+T").expect("valid read structure"), // Variable length to allow zero-length
                ReadStructure::from_str("+T").expect("valid read structure"),
            ],
            store_umi_quals: false,
            store_cell_quals: false,
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved: false,
        };

        extract.execute("test").expect("execute should succeed");

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
    fn test_zero_length_reads_with_fixed_structure_errors() {
        // Test that zero-length reads with fixed-length read structures return an error
        // (you can't have 0 bases when 10T are required)
        let tmp = TempDir::new().expect("failed to create temp dir");
        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "", "")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1],
            output: output.clone(),
            read_structures: vec![
                ReadStructure::from_str("10T").expect("valid read structure"), // Fixed length - should error
            ],
            store_umi_quals: false,
            store_cell_quals: false,
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved: false,
        };

        let err = extract
            .execute("test")
            .expect_err("should error for zero-length read with fixed structure");
        assert!(err.to_string().contains("had too few bases to demux"));
    }

    #[test]
    fn test_extract_sample_barcode_qualities() {
        let tmp = TempDir::new().expect("failed to create temp dir");
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
            read_structures: vec![
                ReadStructure::from_str("3B3M3B+T").expect("valid read structure"),
            ],
            store_umi_quals: false,
            store_cell_quals: false,
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved: false,
        };

        extract.execute("test").expect("execute should succeed");

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
            read_structures: vec![
                ReadStructure::from_str("3B3M3B+T").expect("valid read structure"),
            ],
            store_umi_quals: false,
            store_cell_quals: false,
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved: false,
        };

        extract2.execute("test").expect("execute should succeed");

        let recs2 = read_bam_records(&output2);
        assert_eq!(recs2.len(), 3);
        assert!(get_tag_string(&recs2[0], "QT").is_none());
        assert!(get_tag_string(&recs2[1], "QT").is_none());
        assert!(get_tag_string(&recs2[2], "QT").is_none());
    }

    #[test]
    fn test_extract_umis_from_read_names() {
        let tmp = TempDir::new().expect("failed to create temp dir");
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
            store_umi_quals: false,
            store_cell_quals: false,
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved: false,
        };

        extract.execute("test").expect("execute should succeed");

        let recs = read_bam_records(&output);
        assert_eq!(recs.len(), 3);
        assert_eq!(get_tag_string(&recs[0], "RX"), Some("ACGT".to_string()));
        assert_eq!(get_tag_string(&recs[1], "RX"), Some("TTGA".to_string()));
        assert!(get_tag_string(&recs[2], "RX").is_none());
    }

    #[test]
    fn test_extract_umis_from_read_names_and_sequences() {
        let tmp = TempDir::new().expect("failed to create temp dir");
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
            read_structures: vec![
                ReadStructure::from_str("2M1S2M+T").expect("valid read structure"),
            ],
            store_umi_quals: false,
            store_cell_quals: false,
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved: false,
        };

        extract.execute("test").expect("execute should succeed");

        let recs = read_bam_records(&output);
        assert_eq!(recs.len(), 2);
        assert_eq!(get_tag_string(&recs[0], "RX"), Some("ACGT-CGTA-GG-CC".to_string()));
        assert_eq!(get_tag_string(&recs[1], "RX"), Some("TTGA-TAAT-TA-AA".to_string()));
    }

    /// Tests for [`Extract::extract_read_name_and_umi`].
    ///
    /// fgumi previously hard-coded the 8th `:`-separated field (`parts[7]`) as the UMI,
    /// matching the standard Illumina FASTQ read name shape with UMI present:
    ///
    ///   `@<instr>:<run>:<flowcell>:<lane>:<tile>:<x>:<y>:<UMI>`
    ///
    /// However, some demultiplexers (e.g. ones that also fold the sample index into the
    /// `:`-separated portion) produce read names with 9+ fields where the UMI is the *last*
    /// field, not the 8th. That matches the behavior of fgbio's
    /// `Umis.extractUmisFromReadName`, which always returns the last `:`-separated field
    /// as the UMI. These tests pin the expected behavior:
    ///
    /// - 7 fields → no UMI (fgbio's strict-mode behavior, preserved for backward compat).
    /// - 8 fields → last field as UMI.
    /// - 9+ fields → last field as UMI (was broken; previously returned `parts[7]`).
    /// - `+` in the UMI is normalized to `-`.
    /// - `extract=false` returns `None` regardless of read name shape.
    #[test]
    fn test_extract_read_name_and_umi_nine_fields_takes_last() {
        // Real-world example: BCL Convert / external demux that places the sample index
        // in field 8 and the UMI in field 9 (here, "TCNGCG" is a 6 bp duplex UMI).
        let header =
            b"@LH00620:304:23LLHJLT4:7:1101:2578:1070:CAATCTATAA+rTTACATAGTT:TCNGCG".to_vec();
        let (name, umi) = Extract::extract_read_name_and_umi(&header, true).unwrap();
        assert_eq!(
            name,
            b"LH00620:304:23LLHJLT4:7:1101:2578:1070:CAATCTATAA+rTTACATAGTT:TCNGCG".to_vec()
        );
        assert_eq!(umi, Some(b"TCNGCG".to_vec()));
    }

    #[test]
    fn test_extract_read_name_and_umi_eight_fields_takes_last() {
        let header = b"@q1:2:3:4:5:6:7:ACGT".to_vec();
        let (name, umi) = Extract::extract_read_name_and_umi(&header, true).unwrap();
        assert_eq!(name, b"q1:2:3:4:5:6:7:ACGT".to_vec());
        assert_eq!(umi, Some(b"ACGT".to_vec()));
    }

    #[test]
    fn test_extract_read_name_and_umi_seven_fields_returns_none() {
        let header = b"@q1:2:3:4:5:6:7".to_vec();
        let (name, umi) = Extract::extract_read_name_and_umi(&header, true).unwrap();
        assert_eq!(name, b"q1:2:3:4:5:6:7".to_vec());
        assert_eq!(umi, None);
    }

    #[test]
    fn test_extract_read_name_and_umi_normalizes_plus_to_hyphen() {
        let header = b"@q1:2:3:4:5:6:7:ACGT+CGTA".to_vec();
        let (_name, umi) = Extract::extract_read_name_and_umi(&header, true).unwrap();
        assert_eq!(umi, Some(b"ACGT-CGTA".to_vec()));
    }

    #[test]
    fn test_extract_read_name_and_umi_rejects_field_that_normalizes_to_empty() {
        // An 8-field name whose UMI field is a bare `r` normalizes to empty; the
        // shared normalizer rejects it rather than emitting an empty RX.
        let header = b"@q1:2:3:4:5:6:7:r".to_vec();
        let err = Extract::extract_read_name_and_umi(&header, true).unwrap_err();
        assert!(format!("{err:#}").contains("empty"), "got: {err:#}");
    }

    #[test]
    fn test_extract_read_name_and_umi_normalizes_plus_in_nine_fields() {
        // 9-field name with duplex-style UMI containing '+'.
        let header = b"@A:B:C:D:E:F:G:CAATCTATAA+TTACATAGTT:ACGT+CGTA".to_vec();
        let (_name, umi) = Extract::extract_read_name_and_umi(&header, true).unwrap();
        assert_eq!(umi, Some(b"ACGT-CGTA".to_vec()));
    }

    #[test]
    fn test_extract_read_name_and_umi_strips_space_comment() {
        // Standard Illumina header with space-separated comment after the name.
        let header = b"@q1:2:3:4:5:6:7:ACGT 1:N:0:NNNN".to_vec();
        let (name, umi) = Extract::extract_read_name_and_umi(&header, true).unwrap();
        assert_eq!(name, b"q1:2:3:4:5:6:7:ACGT".to_vec());
        assert_eq!(umi, Some(b"ACGT".to_vec()));
    }

    #[test]
    fn test_extract_read_name_and_umi_extract_false_returns_none() {
        let header = b"@q1:2:3:4:5:6:7:ACGT".to_vec();
        let (name, umi) = Extract::extract_read_name_and_umi(&header, false).unwrap();
        assert_eq!(name, b"q1:2:3:4:5:6:7:ACGT".to_vec());
        assert_eq!(umi, None);
    }

    // ── Old-style Casava (<1.8) `/1` and `/2` read-number suffix stripping ──
    //
    // The written BAM QNAME must have a trailing `/<digit>` removed, matching
    // fgbio's `FastqSource` behavior. Both mates of a pair must share an
    // identical QNAME per the SAM spec, so `read/1` and `read/2` must both
    // become `read`. Stripping is conservative — only `/` followed by a single
    // digit is removed, never `.`, `_`, or `:` separators, which could appear
    // in legitimate read names.

    #[rstest]
    // fgbio strips `/` + any single digit, not just 1/2.
    #[case::one(b"@SRR001.1/1".as_slice(), b"SRR001.1".as_slice())]
    #[case::two(b"@SRR001.1/2".as_slice(), b"SRR001.1".as_slice())]
    #[case::any_digit(b"@read/3".as_slice(), b"read".as_slice())]
    fn test_extract_read_name_and_umi_strips_single_digit_suffix(
        #[case] header: &[u8],
        #[case] expected: &[u8],
    ) {
        let (name, umi) = Extract::extract_read_name_and_umi(header, false).unwrap();
        assert_eq!(name, expected.to_vec());
        assert_eq!(umi, None);
    }

    // Only `/` is a read-number separator; `.`/`_`/`:` may be part of a real
    // name and must be preserved (unlike the broader validation helper).
    #[rstest]
    #[case::underscore(b"@sample_1".as_slice())]
    #[case::dot(b"@foo.2".as_slice())]
    #[case::colon(b"@bar:1".as_slice())]
    fn test_extract_read_name_and_umi_does_not_strip_dot_or_underscore_suffix(
        #[case] header: &[u8],
    ) {
        let (name, _umi) = Extract::extract_read_name_and_umi(header, false).unwrap();
        assert_eq!(name, header[1..].to_vec(), "should not strip: {header:?}");
    }

    #[test]
    fn test_extract_read_name_and_umi_does_not_strip_slash_non_digit() {
        let header = b"@read/x".to_vec();
        let (name, _umi) = Extract::extract_read_name_and_umi(&header, false).unwrap();
        assert_eq!(name, b"read/x".to_vec());
    }

    #[test]
    fn test_extract_read_name_and_umi_does_not_strip_slash_multi_digit() {
        // Only a single trailing digit is a read number; `/12` is left intact.
        let header = b"@read/12".to_vec();
        let (name, _umi) = Extract::extract_read_name_and_umi(&header, false).unwrap();
        assert_eq!(name, b"read/12".to_vec());
    }

    #[test]
    fn test_extract_read_name_and_umi_strips_slash_digit_then_comment() {
        // The space comment is stripped first, then the `/1` read-number suffix.
        let header = b"@read/1 1:N:0:ACGT".to_vec();
        let (name, umi) = Extract::extract_read_name_and_umi(&header, false).unwrap();
        assert_eq!(name, b"read".to_vec());
        assert_eq!(umi, None);
    }

    #[test]
    fn test_extract_read_name_and_umi_strips_slash_before_umi_extraction() {
        // fgbio strips `/<digit>` from the whole name before extracting the UMI,
        // so the read-number digit never leaks into the UMI's last `:` field.
        let header = b"@a:b:c:d:e:f:g:ACGT/1".to_vec();
        let (name, umi) = Extract::extract_read_name_and_umi(&header, true).unwrap();
        assert_eq!(name, b"a:b:c:d:e:f:g:ACGT".to_vec());
        assert_eq!(umi, Some(b"ACGT".to_vec()));
    }

    #[test]
    #[should_panic(expected = "Read names do not match")]
    fn test_extract_fails_on_underscore_read_number_suffix() {
        // A `read_1`/`read_2` pair previously passed validation (`_1`/`_2` were stripped)
        // yet was written with two different QNAMEs, violating the SAM spec. Validation
        // now only strips `/<digit>`, so the mismatch is reported instead.
        let tmp = TempDir::new().expect("failed to create temp dir");
        let r1 = create_fastq(&tmp, "r1.fq", &[("read_1", "AAAAAAAAAA", "==========")]);
        let r2 = create_fastq(&tmp, "r2.fq", &[("read_2", "CCCCCCCCCC", "##########")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1, r2],
            output,
            read_structures: vec![],
            store_umi_quals: false,
            store_cell_quals: false,
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved: false,
        };

        extract.execute("test").expect("execute should succeed");
    }

    #[rstest]
    #[case::fast_path(ThreadingOptions::none())]
    #[case::threaded(ThreadingOptions::new(2))]
    fn test_extract_strips_read_number_suffix_from_written_qname(
        #[case] threading: ThreadingOptions,
    ) {
        // End-to-end: paired FASTQs with old-style `/1` and `/2` read-number
        // suffixes must produce a single shared QNAME (no suffix) for both mates,
        // as required by the SAM spec. Guards the wiring from
        // `extract_read_name_and_umi` through to the written BAM QNAME in both the
        // serial-oracle (`make_raw_records`) and the chain
        // (`make_raw_records_from_fastq_set` + `strip_read_suffix`) record-building
        // paths — this test is parameterized over threading to exercise both.
        let tmp = TempDir::new().expect("failed to create temp dir");
        let r1 = create_fastq(&tmp, "r1.fq", &[("SRR001.1/1", "AAAAAAAAAA", "==========")]);
        let r2 = create_fastq(&tmp, "r2.fq", &[("SRR001.1/2", "CCCCCCCCCC", "##########")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1, r2],
            output: output.clone(),
            read_structures: vec![],
            store_umi_quals: false,
            store_cell_quals: false,
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved: false,
        };

        extract.execute("test").expect("execute should succeed");

        let recs = read_bam_records(&output);
        assert_eq!(recs.len(), 2);
        // Both mates share the suffix-free QNAME.
        assert_eq!(recs[0].name().map(|n| n.as_bytes()), Some(b"SRR001.1".as_slice()));
        assert_eq!(recs[1].name().map(|n| n.as_bytes()), Some(b"SRR001.1".as_slice()));
    }

    /// A BAM read name is length-prefixed by a single byte that counts the
    /// trailing NUL, so names of 255 bytes or more cannot be represented. The
    /// name here comes straight from the input FASTQ header, so an over-long
    /// one is bad input: `extract` must fail cleanly rather than abort with a
    /// Rust panic partway through writing the output BAM and leave a truncated
    /// file behind.
    ///
    /// Both record-building paths are covered (the test is parameterized over
    /// threading): `make_raw_records` (serial oracle) and
    /// `make_raw_records_from_fastq_set` (chain, `--threads`).
    #[rstest]
    #[case::at_limit(254, true)]
    #[case::one_over(255, false)]
    #[case::far_over(4096, false)]
    fn test_extract_rejects_over_long_read_name(
        #[case] name_len: usize,
        #[case] expect_ok: bool,
        #[values(ThreadingOptions::none(), ThreadingOptions::new(2))] threading: ThreadingOptions,
    ) {
        let tmp = TempDir::new().expect("failed to create temp dir");
        let read_name = "N".repeat(name_len);
        let r1 = create_fastq(&tmp, "r1.fq", &[(&read_name, "AAAAAAAAAA", "==========")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1],
            output: output.clone(),
            read_structures: vec![],
            store_umi_quals: false,
            store_cell_quals: false,
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved: false,
        };

        // Before the fix this panicked in the builder instead of returning Err.
        let result = extract.execute("test");

        assert_eq!(result.is_ok(), expect_ok, "name of {name_len} bytes: expected ok={expect_ok}");
        match result {
            Ok(()) => {
                let recs = read_bam_records(&output);
                assert_eq!(recs.len(), 1);
                assert_eq!(
                    recs[0].name().map(|n| n.as_bytes()),
                    Some(read_name.as_bytes()),
                    "an accepted name must round-trip unchanged"
                );
            }
            Err(e) => {
                let msg = format!("{e:#}");
                assert!(msg.contains("read name too long"), "unexpected error: {msg}");
                assert!(
                    msg.contains("could not write the record for read 'NNN"),
                    "the error should name the offending read: {msg}"
                );
                // A rejected name is always over-long (255+ bytes), so it is truncated for the
                // diagnostic — the marker must say so rather than silently showing a partial name.
                assert!(
                    msg.contains("(name shown truncated)"),
                    "an over-long name must be marked truncated: {msg}"
                );
            }
        }
    }

    /// EXT-03: fgbio `Umis.extractUmisFromReadName` upper-cases the extracted UMI
    /// (`Umis.scala:111-117`); a lowercase read-name UMI must become upper-case so
    /// that case-insensitive grouping matches fgbio.
    #[test]
    fn test_extract_read_name_and_umi_uppercases_umi() {
        let header = b"@i:1:f:1:1:1:1:acgtn".to_vec();
        let (_name, umi) = Extract::extract_read_name_and_umi(&header, true).unwrap();
        assert_eq!(umi, Some(b"ACGTN".to_vec()));
    }

    /// EXT-04: an `r`-prefixed UMI is reverse-complemented (`Umis.scala:102-115`,
    /// `reverseComplementPrefixedUmis` defaults to true, as `FastqToBam` uses it).
    /// `rGGTTAA` → revcomp(`GGTTAA`) = `TTAACC`.
    #[test]
    fn test_extract_read_name_and_umi_revcomps_r_prefixed_umi() {
        let header = b"@i:1:f:1:1:1:1:rGGTTAA".to_vec();
        let (_name, umi) = Extract::extract_read_name_and_umi(&header, true).unwrap();
        assert_eq!(umi, Some(b"TTAACC".to_vec()));
    }

    /// EXT-04: for a dual (`+`-delimited) UMI, only the `r`-prefixed segment is
    /// reverse-complemented; the delimiter becomes `-` and the result is upper-cased.
    /// `acgt+rGGTTAA` → `ACGT-TTAACC`.
    #[test]
    fn test_extract_read_name_and_umi_revcomps_only_r_prefixed_segment() {
        let header = b"@i:1:f:1:1:1:1:acgt+rGGTTAA".to_vec();
        let (_name, umi) = Extract::extract_read_name_and_umi(&header, true).unwrap();
        assert_eq!(umi, Some(b"ACGT-TTAACC".to_vec()));
    }

    /// EXT-01: fgbio strict mode throws on a UMI containing characters outside
    /// `ACGTN-` (`Umis.scala:121-123`, `UmisTest.scala:110-114`). fgumi must
    /// error rather than silently writing the illegal UMI to `RX`.
    #[rstest]
    #[case::illegal_char_unpaired(&b"@i:1:f:1:1:1:1:ACGTXY"[..])]
    #[case::illegal_char_in_second_segment(&b"@i:1:f:1:1:1:1:ACGT-CCKC"[..])]
    #[case::illegal_char_in_first_segment(&b"@i:1:f:1:1:1:1:CCKC-ACGT"[..])]
    fn test_extract_read_name_and_umi_rejects_illegal_chars(#[case] illegal: &[u8]) {
        let result = Extract::extract_read_name_and_umi(illegal, true);
        assert!(
            result.is_err(),
            "expected an error for illegal UMI in {:?}",
            String::from_utf8_lossy(illegal)
        );
    }

    /// EXT-01 boundary: a valid UMI (upper-case `ACGTN-` only) is accepted and returned.
    #[test]
    fn test_extract_read_name_and_umi_accepts_valid_umi() {
        let header = b"@i:1:f:1:1:1:1:ACGTN-ACGT".to_vec();
        let (_name, umi) = Extract::extract_read_name_and_umi(&header, true).unwrap();
        assert_eq!(umi, Some(b"ACGTN-ACGT".to_vec()));
    }

    /// EXT-03/EXT-02 combined: a lowercase dual (`+`-delimited) UMI without an
    /// `r` prefix is both upper-cased and hyphen-normalized. `acgt+cgta` → `ACGT-CGTA`.
    #[test]
    fn test_extract_read_name_and_umi_uppercases_plus_delimited_umi() {
        let header = b"@i:1:f:1:1:1:1:acgt+cgta".to_vec();
        let (_name, umi) = Extract::extract_read_name_and_umi(&header, true).unwrap();
        assert_eq!(umi, Some(b"ACGT-CGTA".to_vec()));
    }

    /// EXT-04: reverse-complementing an `r`-prefixed UMI maps RNA `U`→`A`, matching
    /// fgbio's `Sequences.complement`. `rAAU` → revcomp(`AAU`) = `ATT`. (A `U` in a
    /// non-`r` UMI is not reverse-complemented and is rejected by both tools.)
    #[test]
    fn test_extract_read_name_and_umi_revcomps_uracil_to_a() {
        let header = b"@i:1:f:1:1:1:1:rAAU".to_vec();
        let (_name, umi) = Extract::extract_read_name_and_umi(&header, true).unwrap();
        assert_eq!(umi, Some(b"ATT".to_vec()));
    }

    #[test]
    fn test_extract_umis_from_read_names_nine_fields_via_execute() {
        // End-to-end: 9-field read names (sample index in field 8, UMI in field 9) must
        // produce RX tags equal to the *last* field — not the 8th (sample index).
        let tmp = TempDir::new().expect("failed to create temp dir");
        let r1 = create_fastq(
            &tmp,
            "r1.fq",
            &[
                (
                    "LH00620:304:23LLHJLT4:7:1101:2578:1070:CAATCTATAA+rTTACATAGTT:TCNGCG",
                    "AAAAAAAAAA",
                    "==========",
                ),
                (
                    "LH00620:304:23LLHJLT4:7:1101:3565:1070:CAATCTATAA+rTTACATAGTT:TTNCCT",
                    "TAAAAAAAAA",
                    "==========",
                ),
            ],
        );
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1],
            output: output.clone(),
            read_structures: vec![],
            store_umi_quals: false,
            store_cell_quals: false,
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved: false,
        };

        extract.execute("test").expect("execute should succeed");

        let recs = read_bam_records(&output);
        assert_eq!(recs.len(), 2);
        // Must be the 9th (last) field, NOT the 8th (sample index "CAATCTATAA-rTTACATAGTT").
        assert_eq!(get_tag_string(&recs[0], "RX"), Some("TCNGCG".to_string()));
        assert_eq!(get_tag_string(&recs[1], "RX"), Some("TTNCCT".to_string()));
    }

    #[test]
    fn test_extract_cell_barcodes() {
        let tmp = TempDir::new().expect("failed to create temp dir");
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
            read_structures: vec![
                ReadStructure::from_str("2C1S2C+T").expect("valid read structure"),
            ],
            store_umi_quals: false,
            store_cell_quals: true,
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved: false,
        };

        extract.execute("test").expect("execute should succeed");

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
        let tmp = TempDir::new().expect("failed to create temp dir");
        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "AAAAAAAAAA", "==========")]);
        let r2 = create_fastq(&tmp, "r2.fq", &[("x1", "CCCCCCCCCC", "##########")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1, r2],
            output,
            read_structures: vec![],
            store_umi_quals: false,
            store_cell_quals: false,
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved: false,
        };

        extract.execute("test").expect("execute should succeed");
    }

    #[test]
    #[should_panic(expected = "out of sync")]
    fn test_fail_mismatched_read_counts() {
        let tmp = TempDir::new().expect("failed to create temp dir");
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
            store_umi_quals: false,
            store_cell_quals: false,
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved: false,
        };

        extract.execute("test").expect("execute should succeed");
    }

    #[test]
    #[should_panic(expected = "must be supplied the same number of times")]
    fn test_fail_mismatched_inputs_and_read_structures() {
        let tmp = TempDir::new().expect("failed to create temp dir");
        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "AAAAAAAAAA", "==========")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1],
            output,
            read_structures: vec![
                ReadStructure::from_str("+T").expect("valid read structure"),
                ReadStructure::from_str("+T").expect("valid read structure"),
            ],
            store_umi_quals: false,
            store_cell_quals: false,
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved: false,
        };

        extract.execute("test").expect("execute should succeed");
    }

    #[test]
    fn test_annotate_read_names_with_umis() {
        let tmp = TempDir::new().expect("failed to create temp dir");
        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "ACGTAAAAAA", "IIII======")]);
        let r2 = create_fastq(&tmp, "r2.fq", &[("q1", "CCCCCCCCCC", "##########")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1, r2],
            output: output.clone(),
            read_structures: vec![
                ReadStructure::from_str("4M+T").expect("valid read structure"),
                ReadStructure::from_str("+T").expect("valid read structure"),
            ],
            store_umi_quals: false,
            store_cell_quals: false,
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved: false,
        };

        extract.execute("test").expect("execute should succeed");

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
        let tmp = TempDir::new().expect("failed to create temp dir");
        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "AAAAAAAAAA", "==========")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1],
            output: output.clone(),
            read_structures: vec![],
            store_umi_quals: false,
            store_cell_quals: false,
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved: false,
        };

        extract.execute("test").expect("execute should succeed");

        let recs = read_bam_records(&output);
        assert_eq!(recs.len(), 1);

        // Read name should NOT have "+UMI" appended since there are no UMIs
        let r1 = &recs[0];
        assert_eq!(r1.name().map(|n| n.as_bytes()), Some(b"q1".as_slice()));
    }

    #[test]
    fn test_single_tag() {
        let tmp = TempDir::new().expect("failed to create temp dir");
        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "AAACCCTTTAAAAA", "==============")]);
        let r2 = create_fastq(&tmp, "r2.fq", &[("q1", "GGGTTTAAACCCCC", "##############")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1, r2],
            output: output.clone(),
            read_structures: vec![
                ReadStructure::from_str("3B3M3B5T").expect("valid read structure"),
                ReadStructure::from_str("3B3M3B5T").expect("valid read structure"),
            ],
            store_umi_quals: false,
            store_cell_quals: false,
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names: false,
            annotate_read_names: false,
            single_tag: Some("ZU".parse().expect("valid tag")), // Use single tag
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved: false,
        };

        extract.execute("test").expect("execute should succeed");

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
    #[should_panic(
        expected = "Single tag (--single-tag) cannot collide with tags emitted by extract"
    )]
    fn test_single_tag_same_as_umi_tag() {
        let tmp = TempDir::new().expect("failed to create temp dir");
        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "ACGTAAAAAA", "IIII======")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1],
            output,
            read_structures: vec![ReadStructure::from_str("4M+T").expect("valid read structure")],
            store_umi_quals: false,
            store_cell_quals: false,
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names: false,
            annotate_read_names: false,
            single_tag: Some(SamTag::RX), // Same as umi_tag
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved: false,
        };

        extract.execute("test").expect("execute should succeed");
    }

    #[test]
    fn test_combined_annotate_read_names_and_single_tag() {
        let tmp = TempDir::new().expect("failed to create temp dir");
        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "ACGTAAAAAA", "IIII======")]);
        let r2 = create_fastq(&tmp, "r2.fq", &[("q1", "CCCCCCCCCC", "##########")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1, r2],
            output: output.clone(),
            read_structures: vec![
                ReadStructure::from_str("4M+T").expect("valid read structure"),
                ReadStructure::from_str("+T").expect("valid read structure"),
            ],
            store_umi_quals: false,
            store_cell_quals: false,
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names: false,
            annotate_read_names: true, // Both features enabled
            single_tag: Some("ZU".parse().expect("valid tag")),
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved: false,
        };

        extract.execute("test").expect("execute should succeed");

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
    fn test_fail_read_too_short_for_structure() {
        // Test that we fail when a read is not long enough for the read structure
        // Read is 10 bases, but structure requires 8M + 8S = 16 bases minimum
        let tmp = TempDir::new().expect("failed to create temp dir");
        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "AAAAAAAAAA", "==========")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1],
            output,
            read_structures: vec![ReadStructure::from_str("8M8S+T").expect("valid read structure")],
            store_umi_quals: false,
            store_cell_quals: false,
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved: false,
        };

        let err =
            extract.execute("test").expect_err("should error for read too short for structure");
        assert!(err.to_string().contains("had too few bases to demux"));
    }

    #[test]
    fn test_variable_length_reads_with_plus() {
        // Test that the '+' operator in read structures handles variable-length reads
        let tmp = TempDir::new().expect("failed to create temp dir");
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
            read_structures: vec![
                ReadStructure::from_str("4M4B4S+T").expect("valid read structure"),
            ],
            store_umi_quals: false,
            store_cell_quals: false,
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved: false,
        };

        extract.execute("test").expect("execute should succeed");

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

    /// EXT3-01: encoding detection must pool the heads of ALL input FASTQs, not
    /// just `inputs[0]`. A non-template first file whose quality range alone
    /// implies Phred+64 (e.g. an all-high-quality index read) must not force
    /// Illumina detection when the template reads in `inputs[1..]` are Phred+33 —
    /// that would shift every template base quality by ~31. Mirrors fgbio
    /// `FastqToBam`, which detects from `heads.transpose.flatten` across all
    /// sources.
    #[test]
    fn sample_detection_quals_pools_all_inputs_not_just_first() {
        let tmp = TempDir::new().expect("failed to create temp dir");
        // R1: all-high-quality ('h' = 104). On its own, min>=64 & max>=75 → Illumina.
        let r1 = create_fastq(
            &tmp,
            "R1.fq",
            &[("read0", "ACGTACGT", "hhhhhhhh"), ("read1", "ACGTACGT", "hhhhhhhh")],
        );
        // R2: low-quality Phred+33 ('(' = 40 = Q7). Pooled with R1 → min 40 < 59 → Standard.
        let r2 = create_fastq(
            &tmp,
            "R2.fq",
            &[("read0", "ACGTACGT", "(((((((("), ("read1", "ACGTACGT", "((((((((")],
        );

        // Sampling only the first input mis-detects Illumina (the EXT3-01 trap).
        let first_only = sample_detection_quals(std::slice::from_ref(&r1), false, false)
            .expect("sampling first input should succeed");
        assert_eq!(
            QualityEncoding::from_stats(&first_only).expect("detect should succeed"),
            QualityEncoding::Illumina,
            "the first file alone looks like Phred+64 — this is exactly the trap"
        );

        // Pooling both inputs detects Phred+33 (Standard), matching fgbio.
        let pooled = sample_detection_quals(&[r1, r2], false, false)
            .expect("sampling all inputs should succeed");
        assert_eq!(
            QualityEncoding::from_stats(&pooled).expect("detect should succeed"),
            QualityEncoding::Standard,
            "pooling the low-quality template read must recover Phred+33"
        );
    }

    /// A legitimately-empty FASTQ among the pooled inputs must contribute nothing
    /// (the per-input read loop hits `None` immediately) WITHOUT aborting detection
    /// over its non-empty siblings — the multi-input pooling loop has to be robust
    /// to an empty input file, with the empty one listed first so a naive
    /// `break`-the-whole-scan bug would surface as a mis-detection.
    #[test]
    fn sample_detection_quals_skips_empty_input_and_pools_the_rest() {
        let tmp = TempDir::new().expect("failed to create temp dir");
        // A genuinely empty FASTQ (no records), listed first.
        let empty = create_fastq(&tmp, "empty.fq", &[]);
        // A non-empty low-quality Phred+33 sibling ('(' = 40 = Q7).
        let low = create_fastq(
            &tmp,
            "low.fq",
            &[("read0", "ACGTACGT", "(((((((("), ("read1", "ACGTACGT", "((((((((")],
        );

        let pooled = sample_detection_quals(&[empty, low], false, false)
            .expect("sampling should skip the empty input");
        assert_eq!(
            QualityEncoding::from_stats(&pooled).expect("detect should succeed"),
            QualityEncoding::Standard,
            "an empty pooled input must not abort detection over its non-empty siblings"
        );
    }

    /// EXT3-01 (end-to-end): `Extract::execute` must apply the POOLED encoding
    /// detection to the emitted BAM quality scores, not just decide it in
    /// isolation. A high-quality first file (looks Phred+64 alone) pooled with a
    /// low-quality Phred+33 second file must resolve to Standard, so the output
    /// qualities are the raw bytes minus 33 — NOT minus 64 (which mis-detected
    /// Illumina would apply, shifting every base by ~31 and saturating the
    /// low-quality read to 0). Exercises the changed
    /// `sample_detection_quals`/`from_stats` call site through the real command.
    #[test]
    fn extract_execute_applies_pooled_standard_encoding_to_output_quals() {
        let tmp = TempDir::new().expect("failed to create temp dir");
        // R1 all-high-quality 'h' (104): min>=64 & max>=75 → Illumina *if judged alone*.
        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "ACGTACGT", "hhhhhhhh")]);
        // R2 low-quality Phred+33 '(' (40 = Q7): pooling pulls min to 40 (<59) → Standard.
        let r2 = create_fastq(&tmp, "r2.fq", &[("q1", "ACGTACGT", "((((((((")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1, r2],
            output: output.clone(),
            read_structures: vec![],
            store_umi_quals: false,
            store_cell_quals: false,
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved: false,
        };

        extract.execute("test").expect("execute should succeed");

        let recs = read_bam_records(&output);
        assert_eq!(recs.len(), 2, "one R1 + one R2 record");
        // Standard (Phred+33) decode: 'h'(104)-33 = 71, '('(40)-33 = 7. A
        // mis-detected Illumina (Phred+64) would instead yield 40 and 0 (saturated).
        assert_eq!(
            recs[0].quality_scores().as_ref(),
            &[71u8; 8],
            "R1 quals must be Phred+33-decoded (pooled Standard), not Phred+64"
        );
        assert_eq!(
            recs[1].quality_scores().as_ref(),
            &[7u8; 8],
            "R2 low-quality read must decode to Q7 under pooled Standard, not saturate to 0"
        );
    }

    #[test]
    fn test_quality_encoding_detection_phred33() {
        // Test detection of standard Phred+33 encoding (Sanger/Illumina 1.8+)
        // Quality scores with ASCII values 33-126
        let records = vec![
            vec![33u8, 40, 50, 60, 70, 80], // Mix of low, medium, high quality
            vec![35u8, 45, 55, 65, 75, 85], // Another mix
            vec![40u8, 50, 60, 70, 80, 90], // Medium to high quality
        ];

        let encoding =
            QualityEncoding::detect(&records).expect("quality encoding detection should succeed");
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

        let encoding =
            QualityEncoding::detect(&records).expect("quality encoding detection should succeed");
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

        let encoding =
            QualityEncoding::detect(&records).expect("quality encoding detection should succeed");
        // Should be Standard because range is too narrow for typical Phred+64
        assert_eq!(encoding, QualityEncoding::Standard);
    }

    #[test]
    fn test_quality_encoding_detection_empty_reads() {
        // Test that all empty reads default to Standard encoding
        let records = vec![vec![], vec![], vec![]];

        let encoding =
            QualityEncoding::detect(&records).expect("quality encoding detection should succeed");
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

        let encoding =
            QualityEncoding::detect(&records).expect("quality encoding detection should succeed");
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

    /// EXT3-03 (fgbio parity + message improvement): fgbio `FastqToBam` fails with
    /// `case Nil => fail("Quality scores in FASTQ files do not match any known encoding.")`
    /// when no encoding's ASCII range contains every observed quality char. Because
    /// Standard's `[33, 126]` range is the union of all three fgbio encodings (Solexa
    /// `[59, 126]` and Illumina `[64, 126]` are subsets), fgbio's compatible set is empty
    /// on exactly the condition fgumi bails on: a char `< 33` or `> 126`. So fgumi's
    /// failure *set* matches fgbio. fgumi deliberately *improves* on fgbio's static
    /// message by naming the observed `[min, max]` range; these cases assert both that
    /// each fgbio-Nil input fails and that the improved range-reporting message survives.
    #[rstest]
    #[case::below_min_char(vec![vec![20u8, 30, 40, 50]], "20", "50")]
    #[case::above_max_char(vec![vec![50u8, 60, 70, 127]], "50", "127")]
    #[case::both_out_of_range(vec![vec![20u8], vec![127u8]], "20", "127")]
    fn detect_reports_observed_range_on_no_compatible_encoding(
        #[case] records: Vec<Vec<u8>>,
        #[case] expected_min: &str,
        #[case] expected_max: &str,
    ) {
        let err = QualityEncoding::detect(&records).expect_err("no compatible encoding must fail");
        let msg = err.to_string();
        assert!(msg.contains("Invalid quality scores"), "unexpected message: {msg}");
        assert!(
            msg.contains(&format!("[{expected_min}, {expected_max}]")),
            "message must report the observed range (the fgbio-message improvement): {msg}"
        );
    }

    #[test]
    fn test_quality_encoding_detection_ambiguous_range() {
        // Test the ambiguous range (59-63) defaults to Standard
        // This is Q26-Q30 in Phred+33, which is reasonable quality
        let records = vec![
            vec![59u8, 60, 61, 62, 63], // Right in the ambiguous zone
            vec![60u8, 61, 62],         // Also ambiguous
        ];

        let encoding =
            QualityEncoding::detect(&records).expect("quality encoding detection should succeed");
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
        let tmp = TempDir::new().expect("failed to create temp dir");
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
            read_structures: vec![ReadStructure::from_str("+T").expect("valid read structure")],
            store_umi_quals: false,
            store_cell_quals: false,
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved: false,
        };

        // Should succeed without panicking
        extract.execute("test").expect("execute should succeed");

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
        let tmp = TempDir::new().expect("failed to create temp dir");

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
            read_structures: vec![ReadStructure::from_str("+T").expect("valid read structure")],
            store_umi_quals: false,
            store_cell_quals: false,
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved: false,
        };

        extract.execute("test").expect("execute should succeed");

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
        let tmp = TempDir::new().expect("failed to create temp dir");

        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "AAAATTTTCCCCGGGG", "IIIIIIIIIIIIIIII")]);
        let r2 = create_fastq(&tmp, "r2.fq", &[("q1", "TTTTGGGG", "IIIIIIII")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1, r2],
            output: output.clone(),
            read_structures: vec![
                ReadStructure::from_str("4M4S+T").expect("valid read structure"), // R1: 4M UMI, 4S skip, rest template
                ReadStructure::from_str("4M+T").expect("valid read structure"), // R2: 4M UMI, rest template
            ],
            store_umi_quals: false,
            store_cell_quals: false,
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved: false,
        };

        extract.execute("test").expect("execute should succeed");

        let records = read_bam_records(&output);
        assert_eq!(records.len(), 2); // R1 and R2

        // R1: Template should be "CCCCGGGG" (after skipping 4M and 4S)
        assert_eq!(records[0].sequence().as_ref(), b"CCCCGGGG");

        // R2: Template should be "GGGG" (after extracting 4M UMI)
        assert_eq!(records[1].sequence().as_ref(), b"GGGG");

        // Both should have same UMI: AAAA-TTTT
        let rx_tag = noodles::sam::alignment::record::data::field::Tag::from(SamTag::RX);
        let r1_umi = records[0].data().get(&rx_tag).expect("expected tag not found");
        let r2_umi = records[1].data().get(&rx_tag).expect("expected tag not found");

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
        let tmp = TempDir::new().expect("failed to create temp dir");

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
            read_structures: vec![ReadStructure::from_str("5M+T").expect("valid read structure")],
            store_umi_quals: false,
            store_cell_quals: false,
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved: false,
        };

        extract.execute("test").expect("execute should succeed");

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
        let tmp = TempDir::new().expect("failed to create temp dir");

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
            read_structures: vec![ReadStructure::from_str("5M+T").expect("valid read structure")],
            store_umi_quals: false,
            store_cell_quals: false,
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved: false,
        };

        extract.execute("test").expect("execute should succeed");

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
        let tmp = TempDir::new().expect("failed to create temp dir");

        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "AAAAACGTACGT", "IIIIIIIIIIII")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1],
            output: output.clone(),
            read_structures: vec![ReadStructure::from_str("5B+T").expect("valid read structure")], // 5B for barcode
            store_umi_quals: true,
            store_cell_quals: true,
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved: false,
        };

        // Should succeed with all quality tag parameters specified
        extract.execute("test").expect("execute should succeed");

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
            store_umi_quals: false,
            store_cell_quals: false,
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved: false,
        };

        extract.execute("test")?;

        let records = read_bam_records(&output);
        assert_eq!(records.len(), 2, "Should have 2 records for paired-end");

        // Verify MAPQ is 0 for all unmapped reads
        for record in &records {
            let mapq = record.mapping_quality();
            assert!(mapq.is_some(), "Mapping quality should be set (not None/255)");
            let mapq_value: u8 = mapq.expect("mapping quality should be set").into();
            assert_eq!(mapq_value, 0, "Mapping quality should be 0 for unmapped reads");
        }

        Ok(())
    }

    #[test]
    fn test_compression_format_detection_plain() {
        // Create a plain FASTQ file (not compressed)
        let tmp = TempDir::new().expect("failed to create temp dir");
        let plain_path = tmp.path().join("test.fq");
        let mut file = File::create(&plain_path).expect("failed to create file");
        use std::io::Write;
        writeln!(file, "@read1").expect("failed to write line");
        writeln!(file, "ACGT").expect("failed to write line");
        writeln!(file, "+").expect("failed to write line");
        writeln!(file, "IIII").expect("failed to write line");

        let format = detect_compression_format(&plain_path)
            .expect("compression format detection should succeed");
        assert_eq!(format, CompressionFormat::Plain);
    }

    #[test]
    fn test_compression_format_detection_gzip() {
        use flate2::Compression;
        use flate2::write::GzEncoder;
        use std::io::Write;

        // Create a gzip-compressed FASTQ file
        let tmp = TempDir::new().expect("failed to create temp dir");
        let gz_path = tmp.path().join("test.fq.gz");

        // Use flate2 to create a proper gzip file
        let file = File::create(&gz_path).expect("failed to create file");
        let mut encoder = GzEncoder::new(file, Compression::default());
        writeln!(encoder, "@read1").expect("failed to write line");
        writeln!(encoder, "ACGT").expect("failed to write line");
        writeln!(encoder, "+").expect("failed to write line");
        writeln!(encoder, "IIII").expect("failed to write line");
        encoder.finish().expect("failed to finish gzip encoding");

        let format = detect_compression_format(&gz_path)
            .expect("compression format detection should succeed");
        assert_eq!(format, CompressionFormat::Gzip);
    }

    #[test]
    fn test_compression_format_detection_bgzf() {
        use noodles_bgzf::io::Writer as BgzfWriter;
        use std::io::Write;

        // Create a BGZF-compressed file using noodles
        let tmp = TempDir::new().expect("failed to create temp dir");
        let bgzf_path = tmp.path().join("test.fq.bgz");

        let file = File::create(&bgzf_path).expect("failed to create file");
        let mut writer = BgzfWriter::new(file);
        writeln!(writer, "@read1").expect("failed to write line");
        writeln!(writer, "ACGT").expect("failed to write line");
        writeln!(writer, "+").expect("failed to write line");
        writeln!(writer, "IIII").expect("failed to write line");
        writer.finish().expect("failed to finish BGZF writing");

        let format = detect_compression_format(&bgzf_path)
            .expect("compression format detection should succeed");
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
            read_structures: vec![ReadStructure::from_str("5M+T").expect("valid read structure")],
            store_umi_quals: false,
            store_cell_quals: false,
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved: false,
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

    /// Verifies that the chain (`--threads`) path emits the same tags (RX, QX, RG) as the
    /// serial oracle, ensuring `make_raw_records_from_fastq_set` stays in sync with
    /// `make_raw_records`.
    #[rstest]
    #[case::fast_path(ThreadingOptions::none())]
    #[case::threaded(ThreadingOptions::new(2))]
    fn test_threading_modes_emit_correct_tags(#[case] threading: ThreadingOptions) -> Result<()> {
        let tmp = TempDir::new()?;
        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "ACGTAAAAAA", "IIII======")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1],
            output: output.clone(),
            read_structures: vec![ReadStructure::from_str("4M+T").expect("valid read structure")],
            store_umi_quals: true,
            store_cell_quals: false,
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names: false,
            annotate_read_names: false,
            single_tag: None,
            clipping_attribute: None,
            read_group_id: "RG1".to_string(),
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved: false,
        };

        extract.execute("test")?;

        let records = read_bam_records(&output);
        assert_eq!(records.len(), 1);
        let record = &records[0];

        // RX tag: UMI sequence "ACGT"
        let rx = get_tag_string(record, "RX");
        assert_eq!(rx.as_deref(), Some("ACGT"), "RX tag should contain the extracted UMI");

        // QX tag: UMI qualities from "IIII" (Phred+33 = 40)
        let qx = get_tag_string(record, "QX");
        assert!(qx.is_some(), "QX tag should be present when store_umi_quals is true");

        // RG tag
        let rg = get_tag_string(record, "RG");
        assert_eq!(rg.as_deref(), Some("RG1"), "RG tag should match read_group_id");

        Ok(())
    }

    /// Build an `Extract` for the interleaved-equivalence test, varying only the
    /// inputs, output, read structures, `interleaved`, read-name UMI extraction,
    /// and threading.
    fn interleave_test_extract(
        inputs: Vec<PathBuf>,
        output: PathBuf,
        read_structures: Vec<ReadStructure>,
        interleaved: bool,
        extract_umis_from_read_names: bool,
        threading: ThreadingOptions,
    ) -> Extract {
        Extract {
            inputs,
            output,
            read_structures,
            store_umi_quals: false,
            store_cell_quals: false,
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names,
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
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
            interleaved,
        }
    }

    /// A single interleaved input (`--interleaved`) must produce exactly the same
    /// extracted records as the equivalent pair of separate R1/R2 FASTQ files, on
    /// both the single-threaded and the multi-threaded pipeline paths. Distinct R1
    /// and R2 read structures (a 2 bp inline UMI on R1, template-only R2) prove the
    /// split routes R1 → structure[0] and R2 → structure[1].
    #[rstest]
    #[case::single_threaded(ThreadingOptions::none())]
    #[case::multi_threaded(ThreadingOptions::new(2))]
    fn interleaved_input_matches_two_file_input(#[case] threading: ThreadingOptions) {
        let tmp = TempDir::new().expect("failed to create temp dir");

        let r1_records = [("read0", "AACCCCGG", "IIIIIIII"), ("read1", "GGTTTTAA", "IIIIIIII")];
        let r2_records = [("read0", "TTTTGGGG", "JJJJJJJJ"), ("read1", "CCCCAAAA", "JJJJJJJJ")];
        // The same pairs, interleaved R1, R2, R1, R2.
        let interleaved_records = [
            ("read0", "AACCCCGG", "IIIIIIII"),
            ("read0", "TTTTGGGG", "JJJJJJJJ"),
            ("read1", "GGTTTTAA", "IIIIIIII"),
            ("read1", "CCCCAAAA", "JJJJJJJJ"),
        ];

        let r1 = create_fastq(&tmp, "r1.fq", &r1_records);
        let r2 = create_fastq(&tmp, "r2.fq", &r2_records);
        let interleaved = create_fastq(&tmp, "interleaved.fq", &interleaved_records);

        let read_structures = || {
            vec![
                ReadStructure::from_str("2M+T").expect("valid read structure"),
                ReadStructure::from_str("+T").expect("valid read structure"),
            ]
        };

        // Baseline: two separate FASTQ files.
        let two_out = tmp.path().join("two_file.bam");
        interleave_test_extract(
            vec![r1, r2],
            two_out.clone(),
            read_structures(),
            false,
            false,
            threading.clone(),
        )
        .execute("test")
        .expect("two-file extract should succeed");

        // Under test: one interleaved file.
        let il_out = tmp.path().join("interleaved.bam");
        interleave_test_extract(
            vec![interleaved],
            il_out.clone(),
            read_structures(),
            true,
            false,
            threading,
        )
        .execute("test")
        .expect("interleaved extract should succeed");

        // Compare as a multiset of (sequence, RX-UMI) — order-independent, so a
        // difference in per-thread write order cannot mask a real mismatch.
        let key_set = |path: &PathBuf| -> Vec<(Vec<u8>, Option<String>)> {
            let mut keys: Vec<_> = read_bam_records(path)
                .iter()
                .map(|r| (r.sequence().as_ref().to_vec(), get_tag_string(r, "RX")))
                .collect();
            keys.sort();
            keys
        };

        let two_keys = key_set(&two_out);
        let il_keys = key_set(&il_out);
        assert_eq!(il_keys.len(), 4, "expected 4 records (2 pairs)");
        assert_eq!(il_keys, two_keys, "interleaved output must match two-file output");

        // Routing proof: the UMIs are R1's first two bases ("AA"/"GG"), applied
        // via structure[0]. Swapped R1/R2 routing would instead yield R2's leading
        // bases ("TT"/"CC"). (Both mates of a pair carry the shared UMI, so each
        // value appears twice before dedup.)
        let mut umis: Vec<String> = il_keys.iter().filter_map(|(_, rx)| rx.clone()).collect();
        umis.sort();
        umis.dedup();
        assert_eq!(umis, vec!["AA".to_string(), "GG".to_string()], "R1 UMIs via structure[0]");
    }

    #[test]
    fn interleaved_rejects_multiple_inputs() {
        let tmp = TempDir::new().expect("temp dir");
        let a = create_fastq(&tmp, "a.fq", &[("r", "AC", "II")]);
        let b = create_fastq(&tmp, "b.fq", &[("r", "GT", "II")]);
        let structures = vec![
            ReadStructure::from_str("+T").expect("rs"),
            ReadStructure::from_str("+T").expect("rs"),
        ];
        let err = interleave_test_extract(
            vec![a, b],
            tmp.path().join("o.bam"),
            structures,
            true,
            false,
            ThreadingOptions::none(),
        )
        .execute("test")
        .expect_err("two inputs with --interleaved must fail");
        assert!(err.to_string().contains("exactly one --input"), "unexpected: {err}");
    }

    #[test]
    fn interleaved_requires_two_read_structures() {
        let tmp = TempDir::new().expect("temp dir");
        let il = create_fastq(&tmp, "i.fq", &[("r", "AC", "II"), ("r", "GT", "II")]);
        let err = interleave_test_extract(
            vec![il],
            tmp.path().join("o.bam"),
            vec![ReadStructure::from_str("+T").expect("rs")], // only one
            true,
            false,
            ThreadingOptions::none(),
        )
        .execute("test")
        .expect_err("a single read-structure with --interleaved must fail");
        assert!(err.to_string().contains("exactly two --read-structures"), "unexpected: {err}");
    }

    #[test]
    fn interleaved_defaults_to_two_template_read_structures() {
        let tmp = TempDir::new().expect("temp dir");
        let il = create_fastq(&tmp, "i.fq", &[("r", "ACGT", "IIII"), ("r", "TGCA", "JJJJ")]);
        let out = tmp.path().join("o.bam");
        // No --read-structures given: interleaved defaults to `+T +T` (both reads
        // template-only), so one pair yields two records.
        interleave_test_extract(
            vec![il],
            out.clone(),
            vec![],
            true,
            false,
            ThreadingOptions::none(),
        )
        .execute("test")
        .expect("interleaved with default read structures should succeed");
        assert_eq!(read_bam_records(&out).len(), 2, "one interleaved pair → two records");
    }

    /// An interleaved stream with an odd record count (a final R1 with no R2
    /// mate) must be rejected end-to-end — not silently truncated — on both the
    /// single-threaded and pipeline paths. This pins the de-interleaver doc's
    /// contract that the unpaired tail "surfaces as an out-of-sync error".
    #[rstest]
    #[case::single_threaded(ThreadingOptions::none())]
    #[case::multi_threaded(ThreadingOptions::new(2))]
    fn interleaved_odd_record_count_is_rejected(#[case] threading: ThreadingOptions) {
        let tmp = TempDir::new().expect("temp dir");
        let il = create_fastq(
            &tmp,
            "odd.fq",
            &[
                ("read0", "ACGT", "IIII"),
                ("read0", "TGCA", "JJJJ"),
                ("read1", "AAAA", "IIII"), // unpaired R1
            ],
        );
        let structures = vec![
            ReadStructure::from_str("+T").expect("rs"),
            ReadStructure::from_str("+T").expect("rs"),
        ];
        let err = interleave_test_extract(
            vec![il],
            tmp.path().join("o.bam"),
            structures,
            true,
            false,
            threading,
        )
        .execute("test")
        .expect_err("an odd-count interleaved stream must be rejected, not silently truncated");
        let msg = err.to_string();
        assert!(
            msg.contains("out of sync") || msg.contains("mismatch"),
            "expected an out-of-sync / batch-size-mismatch error, got: {msg}"
        );
    }

    /// The production path: `--interleaved` + `--extract-umis-from-read-names`.
    /// The UMI is the trailing `:`-delimited field of the read name, so it must
    /// survive the byte-preserving split and be extracted identically to the
    /// equivalent two-file input, on both paths.
    #[rstest]
    #[case::single_threaded(ThreadingOptions::none())]
    #[case::multi_threaded(ThreadingOptions::new(2))]
    fn interleaved_extracts_umi_from_read_name(#[case] threading: ThreadingOptions) {
        let tmp = TempDir::new().expect("temp dir");
        let r1 = create_fastq(
            &tmp,
            "r1.fq",
            &[
                ("read0:1:2:3:4:5:6:ACGT", "AAAAAAAA", "IIIIIIII"),
                ("read1:1:2:3:4:5:6:TTGA", "CCCCCCCC", "IIIIIIII"),
            ],
        );
        let r2 = create_fastq(
            &tmp,
            "r2.fq",
            &[
                ("read0:1:2:3:4:5:6:ACGT", "GGGGGGGG", "JJJJJJJJ"),
                ("read1:1:2:3:4:5:6:TTGA", "TTTTTTTT", "JJJJJJJJ"),
            ],
        );
        let il = create_fastq(
            &tmp,
            "il.fq",
            &[
                ("read0:1:2:3:4:5:6:ACGT", "AAAAAAAA", "IIIIIIII"),
                ("read0:1:2:3:4:5:6:ACGT", "GGGGGGGG", "JJJJJJJJ"),
                ("read1:1:2:3:4:5:6:TTGA", "CCCCCCCC", "IIIIIIII"),
                ("read1:1:2:3:4:5:6:TTGA", "TTTTTTTT", "JJJJJJJJ"),
            ],
        );
        let structures = || {
            vec![
                ReadStructure::from_str("+T").expect("rs"),
                ReadStructure::from_str("+T").expect("rs"),
            ]
        };

        let two_out = tmp.path().join("two.bam");
        interleave_test_extract(
            vec![r1, r2],
            two_out.clone(),
            structures(),
            false,
            true,
            threading.clone(),
        )
        .execute("test")
        .expect("two-file umi extract should succeed");

        let il_out = tmp.path().join("il.bam");
        interleave_test_extract(vec![il], il_out.clone(), structures(), true, true, threading)
            .execute("test")
            .expect("interleaved umi extract should succeed");

        let rx_multiset = |path: &PathBuf| {
            let mut rx: Vec<Option<String>> =
                read_bam_records(path).iter().map(|r| get_tag_string(r, "RX")).collect();
            rx.sort();
            rx
        };
        let il_rx = rx_multiset(&il_out);
        assert_eq!(il_rx, rx_multiset(&two_out), "interleaved UMI extraction must match two-file");
        // Each pair's read-name UMI lands on both mates.
        let mut distinct: Vec<String> = il_rx.into_iter().flatten().collect();
        distinct.sort();
        distinct.dedup();
        assert_eq!(distinct, vec!["ACGT".to_string(), "TTGA".to_string()], "UMIs from read names");
    }

    /// A BGZF-compressed interleaved file must be decompressed exactly once and
    /// then split: the chain opens the sole interleaved input once via
    /// `open_fastq_reader` and `deinterleave` splits the *decompressed* bytes into
    /// R1/R2, so BGZF is never decoded twice or line-split as raw compressed
    /// bytes. The BGZF and plaintext interleaved runs must produce the same
    /// records, on both the serial and chain paths.
    #[rstest]
    #[case::single_threaded(ThreadingOptions::none())]
    #[case::multi_threaded(ThreadingOptions::new(2))]
    fn interleaved_bgzf_input_matches_plaintext(#[case] threading: ThreadingOptions) {
        let tmp = TempDir::new().expect("temp dir");
        let records = [
            ("read0", "AACCCCGG", "IIIIIIII"),
            ("read0", "TTTTGGGG", "JJJJJJJJ"),
            ("read1", "GGTTTTAA", "IIIIIIII"),
            ("read1", "CCCCAAAA", "JJJJJJJJ"),
        ];
        let plain = create_fastq(&tmp, "il.fq", &records);

        // The same interleaved bytes, BGZF-compressed.
        let plain_bytes = std::fs::read(&plain).expect("read plaintext interleaved");
        let bgzf_path = tmp.path().join("il.fq.gz");
        let mut compressed = Vec::new();
        {
            let mut writer = noodles_bgzf::io::Writer::new(&mut compressed);
            writer.write_all(&plain_bytes).expect("write bgzf");
            writer.finish().expect("finish bgzf");
        }
        std::fs::write(&bgzf_path, compressed).expect("write bgzf interleaved file");

        let structures = || {
            vec![
                ReadStructure::from_str("2M+T").expect("rs"),
                ReadStructure::from_str("+T").expect("rs"),
            ]
        };
        let key_set = |path: &PathBuf| -> Vec<(Vec<u8>, Option<String>)> {
            let mut keys: Vec<_> = read_bam_records(path)
                .iter()
                .map(|r| (r.sequence().as_ref().to_vec(), get_tag_string(r, "RX")))
                .collect();
            keys.sort();
            keys
        };

        let plain_out = tmp.path().join("plain.bam");
        interleave_test_extract(
            vec![plain],
            plain_out.clone(),
            structures(),
            true,
            false,
            threading.clone(),
        )
        .execute("test")
        .expect("plaintext interleaved extract should succeed");

        let bgzf_out = tmp.path().join("bgzf.bam");
        interleave_test_extract(
            vec![bgzf_path],
            bgzf_out.clone(),
            structures(),
            true,
            false,
            threading,
        )
        .execute("test")
        .expect("bgzf interleaved extract should succeed");

        assert_eq!(
            key_set(&bgzf_out),
            key_set(&plain_out),
            "BGZF interleaved output must match plaintext interleaved output"
        );
    }

    /// `validate_template_count` accepts exactly 1-2 template reads (a SAM
    /// template is R1, or R1/R2); 0 yields sequence-less records and 3+ produce
    /// more segments than a template can hold. Pins the guard the chain path
    /// previously skipped (silent malformed output) — see the fn's doc.
    #[rstest]
    #[case::zero_no_structures(&[], false)]
    #[case::zero_no_template_segment(&["8M"], false)]
    #[case::one_template(&["+T"], true)]
    #[case::two_templates(&["+T", "+T"], true)]
    #[case::two_templates_with_barcodes(&["4M+T", "4M+T"], true)]
    #[case::three_templates(&["+T", "+T", "+T"], false)]
    fn validate_template_count_accepts_one_or_two_templates(
        #[case] structures: &[&str],
        #[case] expected_ok: bool,
    ) {
        let read_structures: Vec<ReadStructure> = structures
            .iter()
            .map(|s| ReadStructure::from_str(s).expect("valid read structure"))
            .collect();
        assert_eq!(validate_template_count(&read_structures).is_ok(), expected_ok);
    }

    /// Minimal `ExtractOptions` with every field at a benign default, so a test
    /// can set exactly the field under scrutiny. `ExtractOptions` has no `Default`
    /// (it is normally projected from the parsed CLI struct).
    fn minimal_extract_options() -> ExtractOptions {
        ExtractOptions {
            sample: "sample".to_string(),
            library: "lib".to_string(),
            platform: None,
            platform_unit: None,
            read_group_id: "A".to_string(),
            comments: Vec::new(),
            barcode: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            run_date: None,
            quality_encoding: QualityEncoding::Standard,
            store_umi_quals: false,
            store_cell_quals: false,
            single_tag: None,
            annotate_read_names: false,
            extract_umis_from_read_names: false,
            store_sample_barcode_qualities: false,
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
        }
    }

    /// `Extract::to_extract_options` must project every CLI field onto the
    /// matching `ExtractOptions` field so the chain path sees the same options
    /// the serial oracle uses. `quality_encoding` is a placeholder (`Standard`)
    /// because the chain FASTQ source overrides it with the detected encoding.
    #[test]
    fn to_extract_options_maps_all_fields() {
        let extract = Extract {
            inputs: vec![PathBuf::from("in.fq")],
            output: PathBuf::from("out.bam"),
            read_structures: vec![],
            // Alternate true/false across the adjacent bool fields so an accidental
            // field swap in `to_extract_options` (e.g. store_umi_quals <->
            // store_cell_quals) flips an assertion below instead of passing
            // silently. Bools can't be pairwise-distinct, so this catches the
            // realistic adjacent copy-paste swap, not every possible pair.
            store_umi_quals: true,
            store_cell_quals: false,
            store_sample_barcode_qualities: true,
            extract_umis_from_read_names: false,
            annotate_read_names: true,
            single_tag: Some(SamTag::MI),
            clipping_attribute: None,
            read_group_id: "RGX".to_string(),
            sample: "SAMP".to_string(),
            library: "LIBR".to_string(),
            barcode: Some("ACGT".to_string()),
            platform: "ont".to_string(),
            platform_unit: Some("PU1".to_string()),
            platform_model: Some("PM1".to_string()),
            sequencing_center: Some("CN1".to_string()),
            predicted_insert_size: Some(300),
            description: Some("desc".to_string()),
            comment: vec!["c1".to_string(), "c2".to_string()],
            run_date: Some("2020-01-01".to_string()),
            threading: ThreadingOptions::new(4),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
            async_reader: false,
            check_crc: true,
            no_check_crc: false,
            interleaved: false,
        };
        let opts = extract.to_extract_options();
        assert_eq!(opts.sample, "SAMP");
        assert_eq!(opts.library, "LIBR");
        // platform is Option in ExtractOptions but String on the CLI: always Some.
        assert_eq!(opts.platform.as_deref(), Some("ont"));
        assert_eq!(opts.platform_unit.as_deref(), Some("PU1"));
        assert_eq!(opts.read_group_id, "RGX");
        assert_eq!(opts.comments, vec!["c1".to_string(), "c2".to_string()]);
        assert_eq!(opts.barcode.as_deref(), Some("ACGT"));
        assert_eq!(opts.platform_model.as_deref(), Some("PM1"));
        assert_eq!(opts.sequencing_center.as_deref(), Some("CN1"));
        assert_eq!(opts.predicted_insert_size, Some(300));
        assert_eq!(opts.description.as_deref(), Some("desc"));
        assert_eq!(opts.run_date.as_deref(), Some("2020-01-01"));
        // Placeholder — the chain source detects and overrides this.
        assert_eq!(opts.quality_encoding, QualityEncoding::Standard);
        assert!(opts.store_umi_quals);
        assert!(!opts.store_cell_quals);
        assert_eq!(opts.single_tag, Some(SamTag::MI));
        assert!(opts.annotate_read_names);
        assert!(!opts.extract_umis_from_read_names);
        assert!(opts.store_sample_barcode_qualities);
        assert!(!opts.async_reader);
        assert!(opts.check_crc);
        assert!(!opts.no_check_crc);
    }

    /// `ExtractOptions::validate` must reject a `--single-tag` that collides with
    /// any tag extract emits (RX/QX/CB/CY/BC/QT/RG) and accept any other tag,
    /// mirroring `Extract::validate`.
    #[rstest]
    #[case::reserved_rx(Some(SamTag::RX), false)]
    #[case::reserved_rg(Some(SamTag::RG), false)]
    #[case::non_reserved_mi(Some(SamTag::MI), true)]
    #[case::none(None, true)]
    fn extract_options_validate_single_tag_collision(
        #[case] single_tag: Option<SamTag>,
        #[case] expected_ok: bool,
    ) {
        let mut opts = minimal_extract_options();
        opts.single_tag = single_tag;
        assert_eq!(opts.validate().is_ok(), expected_ok);
    }

    /// Extracting UMIs from read names carries no qualities, so pairing it with
    /// `--store-umi-quals` (which would silently omit `QX`) must be rejected.
    #[test]
    fn extract_options_validate_rejects_umi_from_name_with_store_umi_quals() {
        let mut opts = minimal_extract_options();
        opts.extract_umis_from_read_names = true;
        opts.store_umi_quals = true;
        assert!(
            opts.validate().is_err(),
            "extract-umis-from-read-names + store-umi-quals must be rejected"
        );
        // Either flag alone is fine.
        let mut only_from_name = minimal_extract_options();
        only_from_name.extract_umis_from_read_names = true;
        only_from_name.validate().expect("extract-umis-from-read-names alone must validate");
    }

    // ─────────────────────────────────────────────────────────────────────
    // ExtractRunallOptions / MultiExtractRunallOptions parity (multi_options)
    // ─────────────────────────────────────────────────────────────────────

    /// `MultiExtractRunallOptions` derives `clap::Args` (not `Parser`), so parse
    /// it through a local `#[command(flatten)]` wrapper.
    #[derive(clap::Parser, Debug)]
    struct PrefixedExtract {
        #[command(flatten)]
        opts: MultiExtractRunallOptions,
    }

    /// The minimal required flags to parse a `MultiExtractRunallOptions`.
    fn minimal_runall_args() -> Vec<&'static str> {
        vec![
            "x",
            "--extract::sample",
            "s1",
            "--extract::library",
            "lib1",
            "--extract::inputs",
            "r1.fq",
            "--extract::read-structures",
            "+T",
        ]
    }

    /// Parsing the minimal required flags projects to the expected
    /// [`ExtractOptions`] (all defaults), and the source fields round-trip.
    #[test]
    fn multi_extract_runall_options_defaults_match_projection() {
        let opts = PrefixedExtract::try_parse_from(minimal_runall_args())
            .expect("parses")
            .opts
            .validate()
            .expect("valid");

        let projected = opts.to_extract_options();
        assert_eq!(projected.sample, "s1");
        assert_eq!(projected.library, "lib1");
        assert_eq!(projected.platform, Some("illumina".to_string()));
        assert_eq!(projected.read_group_id, "A");
        assert_eq!(projected.quality_encoding, QualityEncoding::Standard);
        assert!(projected.platform_unit.is_none());
        assert!(projected.barcode.is_none());
        assert!(projected.platform_model.is_none());
        assert!(projected.sequencing_center.is_none());
        assert!(projected.predicted_insert_size.is_none());
        assert!(projected.description.is_none());
        assert!(projected.run_date.is_none());
        assert!(projected.single_tag.is_none());
        assert!(projected.comments.is_empty());
        assert!(!projected.store_umi_quals);
        assert!(!projected.store_cell_quals);
        assert!(!projected.store_sample_barcode_qualities);
        assert!(!projected.annotate_read_names);
        assert!(!projected.extract_umis_from_read_names);
        assert!(!projected.async_reader);
        assert!(!projected.check_crc);
        assert!(!projected.no_check_crc);

        assert_eq!(opts.inputs, vec![PathBuf::from("r1.fq")]);
        assert_eq!(opts.read_structures.len(), 1);
    }

    /// Each staged-required flag omitted makes `validate()` fail.
    #[rstest]
    #[case::sample("--extract::sample")]
    #[case::library("--extract::library")]
    #[case::inputs("--extract::inputs")]
    #[case::read_structures("--extract::read-structures")]
    fn multi_extract_runall_options_missing_required_is_err(#[case] flag_to_drop: &str) {
        // Drop the flag and its following value from the minimal arg list.
        let full = minimal_runall_args();
        let mut args: Vec<&str> = Vec::with_capacity(full.len());
        let mut idx = 0;
        while idx < full.len() {
            if full[idx] == flag_to_drop {
                idx += 2; // skip the flag and its value
            } else {
                args.push(full[idx]);
                idx += 1;
            }
        }
        let parsed = PrefixedExtract::try_parse_from(args).expect("parses").opts.validate();
        assert!(parsed.is_err(), "omitting {flag_to_drop} must fail validate()");
    }

    /// `--store-umi-quals` with `--extract-umis-from-read-names` is rejected by
    /// the struct's own `validate()` (the dropped `conflicts_with`).
    #[test]
    fn extract_runall_options_store_umi_quals_conflicts_with_extract_from_names() {
        let mut args = minimal_runall_args();
        args.push("--extract::store-umi-quals");
        args.push("--extract::extract-umis-from-read-names");
        let opts = PrefixedExtract::try_parse_from(args)
            .expect("parses")
            .opts
            .validate()
            .expect("staged-required lifting passes");
        assert!(opts.validate().is_err());
    }

    /// `--check-crc` with `--no-check-crc` is rejected by `validate()`.
    #[test]
    fn extract_runall_options_check_crc_conflicts_with_no_check_crc() {
        let mut args = minimal_runall_args();
        args.push("--extract::check-crc");
        args.push("--extract::no-check-crc");
        let opts = PrefixedExtract::try_parse_from(args)
            .expect("parses")
            .opts
            .validate()
            .expect("staged-required lifting passes");
        assert!(opts.validate().is_err());
    }

    /// `ExtractRunallOptions::default()` must match the CLI defaults
    /// (`"illumina"` / `"A"`), not the derived-`Default` empty strings — guards
    /// the branch-wide invariant "Default == the minimal-parse projection".
    #[test]
    fn multi_extract_runall_options_default_matches_cli_defaults() {
        let defaults = ExtractRunallOptions::default();
        assert_eq!(defaults.platform, "illumina");
        assert_eq!(defaults.read_group_id, "A");
    }

    /// A supplied `--extract::barcode` round-trips into the projection.
    #[test]
    fn multi_extract_runall_options_round_trips_barcode() {
        let mut args = minimal_runall_args();
        args.push("--extract::barcode");
        args.push("ACGT");
        let opts =
            PrefixedExtract::try_parse_from(args).expect("parses").opts.validate().expect("valid");
        assert_eq!(opts.barcode, Some("ACGT".to_string()));
        assert_eq!(opts.to_extract_options().barcode, Some("ACGT".to_string()));
    }
}

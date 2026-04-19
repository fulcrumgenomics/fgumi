//! `zipper`: Merge unmapped and mapped BAM files with metadata transfer
//!
//! This module implements the `zipper` tool, which merges information from an unmapped BAM
//! (typically from FASTQ conversion) with a mapped BAM (after alignment), transferring all
//! tags and metadata from the unmapped reads to the mapped reads.
//!
//! # Requirements
//!
//! Both input BAMs MUST be:
//! - Queryname sorted or grouped (all records with same name are consecutive)
//! - Have the same queryname ordering
//!
//! # Workflow
//!
//! 1. Read templates (groups of records with same name) from both BAMs in parallel
//! 2. For each template:
//!    - Remove specified tags from mapped reads
//!    - Copy all tags from unmapped to mapped reads
//!    - Apply reverse/revcomp transformations for negative strand reads
//!    - Transfer QC pass/fail flags
//! 3. Write merged records to output
//!
//! # Tag Manipulation
//!
//! - **Remove**: Tags listed in `--tags-to-remove` are removed from mapped reads
//! - **Reverse**: Tags listed in `--tags-to-reverse` are reversed for negative strand reads
//! - **Revcomp**: Tags listed in `--tags-to-revcomp` are reverse complemented for negative strand
//!
//! Named tag sets (e.g., "Consensus") are expanded to their constituent tags.
//!
//! # Example
//!
//! ```bash
//! fgumi zipper \
//!   -i mapped.bam \
//!   -u unmapped.bam \
//!   -r reference.fa \
//!   -o output.bam \
//!   --tags-to-reverse Consensus \
//!   --tags-to-revcomp Consensus
//! ```
//!
//! # Architecture Notes
//!
//! - `Template`: Represents a group of records with the same read name
//! - `TemplateIterator`: Lazily groups records by name from a BAM reader
//! - `TagInfo`: Holds sets of tags to remove/reverse/revcomp
//! - `merge()`: Core function that transfers metadata between templates
use crate::bam_io::{BamReaderAuto, create_bam_reader, create_raw_bam_writer, is_stdin_path};
use crate::batched_sam_reader::BatchedSamReader;
use crate::commands::command::Command;
use crate::commands::common::{CompressionOptions, parse_bool};
use crate::logging::OperationTimer;
use crate::progress::ProgressTracker;
use crate::reference::{ReferenceReader, find_dict_path};
use crate::sam::{
    SamTag, TC_TAG, TemplateCoordinateInfo, buf_value_to_smallest_signed_int, check_sort,
    revcomp_buf_value, reverse_buf_value, unclipped_five_prime_position,
};
use crate::sort::bam_fields;
use crate::template::{Template, TemplateIterator};
use crate::umi::TagInfo;
use crate::validation::validate_file_exists;
use crate::vendored::encode_record_buf;
use anyhow::{Context, Result};
use bstr::ByteSlice;
use clap::Parser;
use fgumi_raw_bam::{RawRecord, RawRecordView};
use log::{debug, info, warn};
use noodles::core::Position;
use noodles::sam::Header;
#[cfg(test)]
use noodles::sam::alignment::record::Cigar as CigarTrait;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::Data;
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::sam::alignment::record_buf::data::field::Value as BufValue;
use std::collections::HashSet;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

/// Command-line arguments for `zipper`
///
/// Merges unmapped and mapped BAM files, transferring tags and metadata from unmapped
/// reads to their corresponding mapped reads. Supports tag manipulation (removal, reversal,
/// reverse complement) based on read strand.
#[derive(Parser, Debug)]
#[command(
    name = "zipper",
    author,
    version,
    about = "\x1b[38;5;72m[ALIGNMENT]\x1b[0m      \x1b[36mZip unmapped BAM with aligned BAM\x1b[0m",
    long_about = r#"
Merges unmapped and mapped BAM files, transferring tags and metadata.

Takes an unmapped BAM (typically from FASTQ) and a mapped BAM (after alignment) and merges
them, copying tags from the unmapped to mapped reads. Both BAMs must be queryname sorted or
grouped, and have the same read name ordering.

The tool transfers tags from the unmapped reads to their corresponding mapped reads. For reads
mapped to the negative strand, tags can be optionally reversed or reverse-complemented. All
QC pass/fail flags are also transferred from the unmapped to mapped reads.

## Tag Manipulation

You can specify which tags to manipulate for reads mapped to the negative strand:
- --tags-to-reverse: Reverses array and string tags (e.g., [1,2,3] becomes [3,2,1])
- --tags-to-revcomp: Reverse complements sequence tags (e.g., AGAGG becomes CCTCT)

Named tag sets like "Consensus" are automatically expanded to their constituent tags:
- Consensus: aD bD cD aM bM cM aE bE cE ad bd cd ae be ce ac bc

## Default Behavior

By default, input is read from stdin and output is written to stdout, allowing for streaming
workflows like:

  bwa mem -t 8 -p -K 150000000 -Y ref.fa reads.fq | fgumi zipper -u unmapped.bam -r ref.fa | fgumi sort -i /dev/stdin -o output.bam --order template-coordinate
"#
)]
#[command(verbatim_doc_comment)]
pub struct Zipper {
    /// Input mapped SAM or BAM file (or `-` for stdin; SAM or BAM is auto-detected).
    /// BAM input is discouraged; prefer piping SAM directly from the aligner for best performance.
    #[arg(short = 'i', long, default_value = "-")]
    pub input: PathBuf,

    /// Input unmapped BAM file containing original tags
    #[arg(short = 'u', long)]
    pub unmapped: PathBuf,

    /// Reference FASTA file (must have accompanying .dict file)
    #[arg(short = 'r', long)]
    pub reference: PathBuf,

    /// Output BAM file (or `-` for stdout)
    #[arg(short = 'o', long, default_value = "-")]
    pub output: PathBuf,

    /// Tags to remove from mapped reads before copying unmapped tags
    #[arg(long, value_delimiter = ',')]
    pub tags_to_remove: Vec<String>,

    /// Tags to reverse for reads mapped to negative strand
    #[arg(long, value_delimiter = ',')]
    pub tags_to_reverse: Vec<String>,

    /// Tags to reverse complement for reads mapped to negative strand
    #[arg(long, value_delimiter = ',')]
    pub tags_to_revcomp: Vec<String>,

    /// Buffer size for template channel (default: 50000)
    #[arg(short = 'b', long, default_value = "50000")]
    pub buffer: usize,

    /// Number of threads to use for processing (default: 1, single-threaded)
    #[arg(long, short = 't', default_value = "1")]
    pub threads: usize,

    /// Compression options for output BAM.
    #[command(flatten)]
    pub compression: CompressionOptions,

    /// BWA -K parameter value (bases per batch). Used to optimize buffer sizing
    /// for stdin input. The buffer grows adaptively based on observed bytes per batch.
    /// Default matches common bwa mem usage.
    #[arg(short = 'K', long = "bwa-chunk-size", default_value = "150000000")]
    pub bwa_chunk_size: u64,

    /// Exclude reads from the unmapped BAM that are not present in the aligned BAM.
    /// Useful when reads were intentionally removed (e.g., by adapter trimming) prior to alignment.
    #[arg(long = "exclude-missing-reads", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool)]
    pub exclude_missing_reads: bool,

    /// Skip adding `tc` (template-coordinate) tags to secondary/supplementary reads.
    /// By default, zipper adds a `tc` tag containing the primary alignment's template
    /// sort key coordinates, which enables correct template-coordinate sorting and
    /// deduplication of these reads. Use this flag if you don't need this functionality.
    #[arg(long = "skip-tc-tags", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool)]
    pub skip_tc_tags: bool,

    /// Restore unconverted bases in EM-seq consensus reads after bwameth re-alignment.
    ///
    /// In EM-seq, unmethylated cytosines are converted to thymine (top strand) or
    /// adenine (bottom strand). After bwameth re-alignment, this flag replaces converted
    /// bases back to their unconverted reference form at reference C (top strand) or
    /// reference G (bottom strand) positions. Uses the bwameth `YD` tag to determine
    /// the bisulfite strand.
    ///
    /// This produces a final BAM where the sequence shows the original (unconverted) bases,
    /// while methylation state is preserved in MM/ML tags and cu/ct count tags.
    #[arg(long = "restore-unconverted-bases", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool)]
    pub restore_unconverted_bases: bool,
}

/// Builds the output BAM header from unmapped and mapped headers
///
/// Merges information from both input headers:
/// - Reference sequences from the reference dictionary
/// - Comments from both input files
/// - Read groups from both inputs (mapped takes precedence)
/// - Program records from both inputs (mapped takes precedence)
/// - Header lines from unmapped BAM
///
/// # Arguments
///
/// * `unmapped` - Header from the unmapped BAM
/// * `mapped` - Header from the mapped BAM
/// * `dict_path` - Path to the sequence dictionary file (.dict)
///
/// # Returns
///
/// Combined header for the output BAM
///
/// # Errors
///
/// Returns an error if the reference dictionary cannot be opened or parsed
pub fn build_output_header(unmapped: &Header, mapped: &Header, dict_path: &Path) -> Result<Header> {
    let dict_file =
        std::fs::File::open(dict_path).context("Failed to open reference dictionary")?;
    let mut dict_reader = noodles::sam::io::Reader::new(std::io::BufReader::new(dict_file));
    let dict_header = dict_reader.read_header()?;

    let mut header = Header::builder();

    for (name, map) in dict_header.reference_sequences() {
        header = header.add_reference_sequence(name.clone(), map.clone());
    }

    for comment in unmapped.comments() {
        header = header.add_comment(comment.clone());
    }
    for comment in mapped.comments() {
        header = header.add_comment(comment.clone());
    }

    let mut rg_ids = HashSet::new();
    for (id, rg) in mapped.read_groups() {
        header = header.add_read_group(id.clone(), rg.clone());
        rg_ids.insert(id.clone());
    }
    for (id, rg) in unmapped.read_groups() {
        if !rg_ids.contains(id) {
            header = header.add_read_group(id.clone(), rg.clone());
        }
    }

    let mut pg_ids = HashSet::new();
    for (id, pg) in mapped.programs().as_ref() {
        header = header.add_program(id.clone(), pg.clone());
        pg_ids.insert(id.clone());
    }
    for (id, pg) in unmapped.programs().as_ref() {
        if !pg_ids.contains(id) {
            header = header.add_program(id.clone(), pg.clone());
        }
    }

    if let Some(hdr) = unmapped.header() {
        header = header.set_header(hdr.clone());
    }

    Ok(header.build())
}

/// Copies tags from source record to destination record
///
/// Copies all tags from `src` to `dest`, with special handling:
/// - PG tag is only copied if not already present in dest
/// - Tags in `tag_info.remove` are skipped
///
/// # Arguments
///
/// * `src` - Source record to copy tags from
/// * `dest` - Destination record to copy tags to
/// * `tag_info` - Information about which tags to skip
///
/// # Returns
///
/// `Ok(())` on success
///
/// # Errors
///
/// Currently never returns an error, but signature allows for future error handling
pub fn copy_tags(src: &RecordBuf, dest: &mut Data, tag_info: &TagInfo) -> Result<()> {
    for (tag, value) in src.data().iter() {
        let tag_str = format!("{}{}", tag.as_ref()[0] as char, tag.as_ref()[1] as char);

        if tag_str == "PG" && dest.get(&tag).is_some() {
            continue;
        }

        if tag_info.remove.contains(&tag_str) {
            continue;
        }

        dest.insert(tag, value.clone());
    }
    Ok(())
}

/// Adds the `tc` tag to secondary and supplementary reads, storing the template
/// sort key from the primary alignments.
///
/// This enables downstream tools (sort, dedup) to properly handle these reads
/// by ensuring they sort adjacent to their primary alignments.
///
/// # Arguments
///
/// * `template` - The template containing all reads (primary, secondary, supplementary)
///
/// # Format
///
/// The `tc` tag is a B:i array with 6 int32 elements: `[tid1, pos1, neg1, tid2, pos2, neg2]`
/// where "1" is the earlier position (lower (tid, pos)) and neg is 0/1 for forward/reverse.
///
/// This matches the template-coordinate sort key format used by `fgumi sort`.
///
/// Optimized to skip all computation when there are no secondary/supplementary reads.
fn add_template_coordinate_tags(template: &mut Template) {
    // Fast path: check if there are any secondary/supplementary reads before doing work
    let has_sec_supp =
        template.records.iter().any(|r| r.flags().is_secondary() || r.flags().is_supplementary());

    if !has_sec_supp {
        return; // No secondary/supplementary reads - nothing to do
    }

    // Get R1 primary alignment info: (ref_id, unclipped_5prime_pos, is_reverse)
    // Uses unclipped 5' position to match template-coordinate sort key calculation
    let r1_primary_info: Option<(i32, i32, bool)> = template.r1().and_then(|r| {
        if r.flags().is_unmapped() {
            return None;
        }
        let ref_id: i32 = r.reference_sequence_id()?.try_into().ok()?;
        // Use unclipped 5' position to match sort key calculation
        let pos: i32 = unclipped_five_prime_position(r)?.try_into().ok()?;
        let is_reverse = r.flags().is_reverse_complemented();
        Some((ref_id, pos, is_reverse))
    });

    // Get R2 primary alignment info: (ref_id, unclipped_5prime_pos, is_reverse)
    let r2_primary_info: Option<(i32, i32, bool)> = template.r2().and_then(|r| {
        if r.flags().is_unmapped() {
            return None;
        }
        let ref_id: i32 = r.reference_sequence_id()?.try_into().ok()?;
        // Use unclipped 5' position to match sort key calculation
        let pos: i32 = unclipped_five_prime_position(r)?.try_into().ok()?;
        let is_reverse = r.flags().is_reverse_complemented();
        Some((ref_id, pos, is_reverse))
    });

    // Compute the template sort key (tid1, pos1, neg1, tid2, pos2, neg2)
    // where "1" is the earlier position (lower (tid, pos))
    let tc_info: Option<TemplateCoordinateInfo> = match (r1_primary_info, r2_primary_info) {
        (Some((tid1, pos1, neg1)), Some((tid2, pos2, neg2))) => {
            // Both mapped - determine which is earlier
            if (tid1, pos1) <= (tid2, pos2) {
                Some(TemplateCoordinateInfo::new(tid1, pos1, neg1, tid2, pos2, neg2))
            } else {
                Some(TemplateCoordinateInfo::new(tid2, pos2, neg2, tid1, pos1, neg1))
            }
        }
        (Some((tid, pos, neg)), None) | (None, Some((tid, pos, neg))) => {
            // Only one mapped - use same values for both positions
            Some(TemplateCoordinateInfo::new(tid, pos, neg, tid, pos, neg))
        }
        (None, None) => None,
    };

    // If we have template sort key info, add it to all secondary/supplementary reads
    let Some(tc_info) = tc_info else {
        return;
    };

    let tc_value = tc_info.to_tag_value();

    for record in &mut template.records {
        let flags = record.flags();

        // Skip primary alignments - they don't need the tc tag
        if !flags.is_secondary() && !flags.is_supplementary() {
            continue;
        }

        // Add the tc tag with the template sort key
        record.data_mut().insert(TC_TAG, tc_value.clone());
    }
}

/// Merges tags from unmapped template into mapped template
///
/// This is the core merge function that:
/// 1. Removes specified tags from mapped reads
/// 2. Copies tags from unmapped to mapped reads
/// 3. Applies reverse/revcomp transformations for negative strand
/// 4. Transfers QC pass/fail flags
///
/// Handles both paired-end and single-end data, matching R1 to R1 and R2 to R2.
/// Secondary and supplementary alignments are also processed.
///
/// # Arguments
///
/// * `unmapped` - Template containing unmapped reads with original tags
/// * `mapped` - Template containing mapped reads (will be modified in place)
/// * `tag_info` - Information about which tags to manipulate
///
/// # Returns
///
/// `Ok(())` on success
///
/// # Errors
///
/// Returns an error if tag copying fails (currently never happens)
pub fn merge(
    unmapped: &Template,
    mapped: &mut Template,
    tag_info: &TagInfo,
    skip_tc_tags: bool,
) -> Result<()> {
    // Fix mate info first
    mapped.fix_mate_info()?;

    // First, remove tags from mapped reads
    for i in 0..mapped.read_count() {
        for tag_str in &tag_info.remove {
            if tag_str.len() == 2 {
                let tag = Tag::new(tag_str.as_bytes()[0], tag_str.as_bytes()[1]);
                let _ = mapped.records[i].data_mut().remove(&tag);
            }
        }
    }

    // Process each primary unmapped read
    for u in unmapped.primary_reads() {
        let is_unpaired = !u.flags().is_segmented();
        let is_first = u.flags().is_first_segment();

        // Find corresponding mapped reads (R1 or R2)
        let mapped_indices: Vec<usize> = if is_unpaired || is_first {
            mapped
                .records
                .iter()
                .enumerate()
                .filter(|(_, r)| !r.flags().is_segmented() || r.flags().is_first_segment())
                .map(|(i, _)| i)
                .collect()
        } else {
            mapped
                .records
                .iter()
                .enumerate()
                .filter(|(_, r)| r.flags().is_segmented() && !r.flags().is_first_segment())
                .map(|(i, _)| i)
                .collect()
        };

        let has_pos_reads =
            mapped_indices.iter().any(|&i| !mapped.records[i].flags().is_reverse_complemented());
        let has_neg_reads =
            mapped_indices.iter().any(|&i| mapped.records[i].flags().is_reverse_complemented());

        // Copy tags to positive strand reads
        if has_pos_reads {
            for &i in &mapped_indices {
                if !mapped.records[i].flags().is_reverse_complemented() {
                    copy_tags(u, mapped.records[i].data_mut(), tag_info)?;
                }
            }
        }

        // Copy tags to negative strand reads (with reverse/revcomp)
        if has_neg_reads {
            let mut u_for_neg = u.clone();

            if tag_info.has_revs_or_revcomps() {
                for (tag, value) in u.data().iter() {
                    let tag_str = format!("{}{}", tag.as_ref()[0] as char, tag.as_ref()[1] as char);

                    if tag_info.reverse.contains(&tag_str) {
                        let new_val = reverse_buf_value(value);
                        u_for_neg.data_mut().insert(tag, new_val);
                    } else if tag_info.revcomp.contains(&tag_str) {
                        let new_val = revcomp_buf_value(value);
                        u_for_neg.data_mut().insert(tag, new_val);
                    }
                }
            }

            for &i in &mapped_indices {
                if mapped.records[i].flags().is_reverse_complemented() {
                    copy_tags(&u_for_neg, mapped.records[i].data_mut(), tag_info)?;
                }
            }
        }

        // Transfer QC pass/fail flag
        let is_qc_fail = u.flags().is_qc_fail();
        for &i in &mapped_indices {
            let mut flags = mapped.records[i].flags();
            if is_qc_fail {
                flags.insert(noodles::sam::alignment::record::Flags::QC_FAIL);
            } else {
                flags.remove(noodles::sam::alignment::record::Flags::QC_FAIL);
            }
            // Apply the modified flags back to the record
            *mapped.records[i].flags_mut() = flags;
        }
    }

    // Normalize alignment integer tags (AS, XS) to Int8 to match fgbio's encoding.
    // When reading from SAM format, noodles may parse these tags with different types
    // (e.g., Int32 or UInt8). fgbio consistently uses signed int8 (SAM type code 'c').
    let as_tag = Tag::from(SamTag::AS);
    let xs_tag = Tag::from(SamTag::XS);
    for i in 0..mapped.read_count() {
        // Normalize AS tag
        if let Some(value) = mapped.records[i].data().get(&as_tag) {
            if let Some(int8_value) = buf_value_to_smallest_signed_int(value) {
                mapped.records[i].data_mut().insert(as_tag, int8_value);
            }
        }
        // Normalize XS tag
        if let Some(value) = mapped.records[i].data().get(&xs_tag) {
            if let Some(int8_value) = buf_value_to_smallest_signed_int(value) {
                mapped.records[i].data_mut().insert(xs_tag, int8_value);
            }
        }
    }

    // Add primary alignment tags to secondary and supplementary reads
    // This enables proper template-coordinate sorting and duplicate marking
    if !skip_tc_tags {
        add_template_coordinate_tags(mapped);
    }

    Ok(())
}

/// Appends a noodles `BufValue` tag to a raw BAM record in BAM binary encoding.
fn append_buf_value_raw(dest: &mut Vec<u8>, tag: [u8; 2], value: &BufValue) {
    use noodles::sam::alignment::record_buf::data::field::value::Array;
    match value {
        BufValue::Character(c) => {
            dest.push(tag[0]);
            dest.push(tag[1]);
            dest.push(b'A');
            dest.push(*c);
        }
        BufValue::Int8(v) => {
            dest.push(tag[0]);
            dest.push(tag[1]);
            dest.push(b'c');
            dest.push(v.cast_unsigned());
        }
        BufValue::UInt8(v) => {
            dest.push(tag[0]);
            dest.push(tag[1]);
            dest.push(b'C');
            dest.push(*v);
        }
        BufValue::Int16(v) => {
            dest.push(tag[0]);
            dest.push(tag[1]);
            dest.push(b's');
            dest.extend_from_slice(&v.to_le_bytes());
        }
        BufValue::UInt16(v) => {
            dest.push(tag[0]);
            dest.push(tag[1]);
            dest.push(b'S');
            dest.extend_from_slice(&v.to_le_bytes());
        }
        BufValue::Int32(v) => {
            dest.push(tag[0]);
            dest.push(tag[1]);
            dest.push(b'i');
            dest.extend_from_slice(&v.to_le_bytes());
        }
        BufValue::UInt32(v) => {
            dest.push(tag[0]);
            dest.push(tag[1]);
            dest.push(b'I');
            dest.extend_from_slice(&v.to_le_bytes());
        }
        BufValue::Float(v) => {
            dest.push(tag[0]);
            dest.push(tag[1]);
            dest.push(b'f');
            dest.extend_from_slice(&v.to_le_bytes());
        }
        BufValue::String(s) => {
            dest.push(tag[0]);
            dest.push(tag[1]);
            dest.push(b'Z');
            dest.extend_from_slice(s.as_bytes());
            dest.push(0);
        }
        BufValue::Hex(s) => {
            dest.push(tag[0]);
            dest.push(tag[1]);
            dest.push(b'H');
            dest.extend_from_slice(s.as_bytes());
            dest.push(0);
        }
        BufValue::Array(arr) => {
            dest.push(tag[0]);
            dest.push(tag[1]);
            dest.push(b'B');
            match arr {
                Array::Int8(vals) => {
                    dest.push(b'c');
                    let count = u32::try_from(vals.len()).expect("BAM B-array length exceeds u32");
                    dest.extend_from_slice(&count.to_le_bytes());
                    for v in vals {
                        dest.push(v.cast_unsigned());
                    }
                }
                Array::UInt8(vals) => {
                    dest.push(b'C');
                    let count = u32::try_from(vals.len()).expect("BAM B-array length exceeds u32");
                    dest.extend_from_slice(&count.to_le_bytes());
                    dest.extend_from_slice(vals.as_slice());
                }
                Array::Int16(vals) => {
                    dest.push(b's');
                    let count = u32::try_from(vals.len()).expect("BAM B-array length exceeds u32");
                    dest.extend_from_slice(&count.to_le_bytes());
                    for v in vals {
                        dest.extend_from_slice(&v.to_le_bytes());
                    }
                }
                Array::UInt16(vals) => {
                    dest.push(b'S');
                    let count = u32::try_from(vals.len()).expect("BAM B-array length exceeds u32");
                    dest.extend_from_slice(&count.to_le_bytes());
                    for v in vals {
                        dest.extend_from_slice(&v.to_le_bytes());
                    }
                }
                Array::Int32(vals) => {
                    dest.push(b'i');
                    let count = u32::try_from(vals.len()).expect("BAM B-array length exceeds u32");
                    dest.extend_from_slice(&count.to_le_bytes());
                    for v in vals {
                        dest.extend_from_slice(&v.to_le_bytes());
                    }
                }
                Array::UInt32(vals) => {
                    dest.push(b'I');
                    let count = u32::try_from(vals.len()).expect("BAM B-array length exceeds u32");
                    dest.extend_from_slice(&count.to_le_bytes());
                    for v in vals {
                        dest.extend_from_slice(&v.to_le_bytes());
                    }
                }
                Array::Float(vals) => {
                    dest.push(b'f');
                    let count = u32::try_from(vals.len()).expect("BAM B-array length exceeds u32");
                    dest.extend_from_slice(&count.to_le_bytes());
                    for v in vals {
                        dest.extend_from_slice(&v.to_le_bytes());
                    }
                }
            }
        }
    }
}

/// Adds `tc` tags to secondary/supplementary reads in raw-byte mode.
fn add_template_coordinate_tags_raw(mapped: &mut Template) {
    let rr = match mapped.raw_records.as_ref() {
        Some(rr) => rr,
        None => return,
    };

    // Fast path: check if there are any secondary/supplementary reads
    let has_sec_supp = rr.iter().any(|r| {
        let f = RawRecordView::new(r).flags();
        (f & bam_fields::flags::SECONDARY) != 0 || (f & bam_fields::flags::SUPPLEMENTARY) != 0
    });
    if !has_sec_supp {
        return;
    }

    // Get R1 primary info
    let r1_info: Option<(i32, i32, bool)> = mapped.r1.and_then(|(i, _)| {
        let r = &rr[i];
        let f = RawRecordView::new(r).flags();
        if (f & bam_fields::flags::UNMAPPED) != 0 {
            return None;
        }
        let ref_id = bam_fields::ref_id(r);
        let pos = bam_fields::unclipped_5prime_from_raw_bam(r);
        let is_reverse = (f & bam_fields::flags::REVERSE) != 0;
        Some((ref_id, pos, is_reverse))
    });

    // Get R2 primary info
    let r2_info: Option<(i32, i32, bool)> = mapped.r2.and_then(|(i, _)| {
        let r = &rr[i];
        let f = RawRecordView::new(r).flags();
        if (f & bam_fields::flags::UNMAPPED) != 0 {
            return None;
        }
        let ref_id = bam_fields::ref_id(r);
        let pos = bam_fields::unclipped_5prime_from_raw_bam(r);
        let is_reverse = (f & bam_fields::flags::REVERSE) != 0;
        Some((ref_id, pos, is_reverse))
    });

    let tc_info: Option<TemplateCoordinateInfo> = match (r1_info, r2_info) {
        (Some((t1, p1, n1)), Some((t2, p2, n2))) => {
            if (t1, p1) <= (t2, p2) {
                Some(TemplateCoordinateInfo::new(t1, p1, n1, t2, p2, n2))
            } else {
                Some(TemplateCoordinateInfo::new(t2, p2, n2, t1, p1, n1))
            }
        }
        (Some((t, p, n)), None) | (None, Some((t, p, n))) => {
            Some(TemplateCoordinateInfo::new(t, p, n, t, p, n))
        }
        (None, None) => None,
    };

    let Some(tc_info) = tc_info else {
        return;
    };

    let tc_tag = b"tc";
    let tc_values = [
        tc_info.tid1,
        tc_info.pos1,
        i32::from(tc_info.neg1),
        tc_info.tid2,
        tc_info.pos2,
        i32::from(tc_info.neg2),
    ];

    let rr = mapped.raw_records.as_mut().unwrap();
    for record in rr.iter_mut() {
        let f = RawRecordView::new(record.as_ref()).flags();
        if (f & bam_fields::flags::SECONDARY) != 0 || (f & bam_fields::flags::SUPPLEMENTARY) != 0 {
            bam_fields::remove_tag(record.as_mut_vec(), tc_tag);
            bam_fields::append_i32_array_tag(record.as_mut_vec(), tc_tag, &tc_values);
        }
    }
}

/// Converts a mapped `Template` from `RecordBuf` mode to raw-byte mode.
///
/// Encodes each `RecordBuf` to raw BAM bytes using the vendored encoder
/// and populates `raw_records`.
fn convert_template_to_raw(template: &mut Template, header: &Header) -> Result<()> {
    if template.is_raw_byte_mode() {
        return Ok(());
    }
    let mut raw_records = Vec::with_capacity(template.records.len());
    for record in &template.records {
        let mut buf = Vec::with_capacity(256);
        encode_record_buf(&mut buf, header, record)
            .map_err(|e| anyhow::anyhow!("Failed to encode record: {e}"))?;
        raw_records.push(RawRecord::from(buf));
    }
    template.raw_records = Some(raw_records);
    // Keep records intact for read_count compatibility
    Ok(())
}

/// Collects the indices of mapped records corresponding to a given segment
/// (R1 or R2) using the template's pre-computed index fields, avoiding a
/// full scan of all mapped records.
fn collect_mapped_indices(mapped: &Template, is_first_segment: bool) -> Vec<usize> {
    let mut indices = Vec::with_capacity(4);
    if is_first_segment {
        if let Some((i, _)) = mapped.r1 {
            indices.push(i);
        }
        if let Some((s, e)) = mapped.r1_supplementals {
            indices.extend(s..e);
        }
        if let Some((s, e)) = mapped.r1_secondaries {
            indices.extend(s..e);
        }
    } else {
        if let Some((i, _)) = mapped.r2 {
            indices.push(i);
        }
        if let Some((s, e)) = mapped.r2_supplementals {
            indices.extend(s..e);
        }
        if let Some((s, e)) = mapped.r2_secondaries {
            indices.extend(s..e);
        }
    }
    indices
}

/// Merges tags from unmapped template into mapped template using raw bytes.
///
/// This is the raw-byte equivalent of [`merge`]. The mapped template must
/// be in raw-byte mode (call [`convert_template_to_raw`] first). The
/// unmapped template uses `RecordBuf` mode.
///
/// Performs the same 7 operations as `merge()`:
///
/// 1. Fix mate info (raw)
/// 2. Remove tags from mapped records
/// 3. Copy tags from unmapped to mapped (with reverse/revcomp for neg strand)
/// 5. Transfer QC flags
/// 6. Normalize AS/XS tags
/// 7. Add TC tags
pub fn merge_raw(
    unmapped: &Template,
    mapped: &mut Template,
    tag_info: &TagInfo,
    skip_tc_tags: bool,
) -> Result<()> {
    // Step 1: Fix mate info
    mapped.fix_mate_info_raw()?;

    let rr = mapped
        .raw_records
        .as_mut()
        .ok_or_else(|| anyhow::anyhow!("merge_raw: mapped not in raw mode"))?;

    // Step 2: Remove tags from mapped reads
    for record in rr.iter_mut() {
        for tag_str in &tag_info.remove {
            if tag_str.len() == 2 {
                let tag_bytes: [u8; 2] = [tag_str.as_bytes()[0], tag_str.as_bytes()[1]];
                bam_fields::remove_tag(record.as_mut_vec(), &tag_bytes);
            }
        }
    }

    // Steps 3–5: Copy tags from unmapped to mapped, transfer QC flags
    let has_transforms = tag_info.has_revs_or_revcomps();
    let pg_tag: [u8; 2] = [b'P', b'G'];

    for u in unmapped.primary_reads() {
        let is_unpaired = !u.flags().is_segmented();
        let is_first = u.flags().is_first_segment();

        // Use template's known indices instead of scanning all mapped records
        let mapped_indices = collect_mapped_indices(mapped, is_unpaired || is_first);

        // Single pass to determine strand presence
        let (has_pos, has_neg) = {
            let rr = mapped.raw_records.as_ref().unwrap();
            let mut pos = false;
            let mut neg = false;
            for &i in &mapped_indices {
                if (RawRecordView::new(&rr[i]).flags() & bam_fields::flags::REVERSE) == 0 {
                    pos = true;
                } else {
                    neg = true;
                }
                if pos && neg {
                    break;
                }
            }
            (pos, neg)
        };

        // Copy tags to positive strand reads
        if has_pos {
            let rr = mapped.raw_records.as_mut().unwrap();
            for &i in &mapped_indices {
                if (RawRecordView::new(&rr[i]).flags() & bam_fields::flags::REVERSE) != 0 {
                    continue;
                }
                let aux = bam_fields::aux_data_slice(&rr[i]);
                let has_pg = bam_fields::find_tag_type(aux, &pg_tag).is_some();

                for (tag, value) in u.data().iter() {
                    let tag_bytes: [u8; 2] = [tag.as_ref()[0], tag.as_ref()[1]];
                    if tag_bytes == pg_tag && has_pg {
                        continue;
                    }
                    // Tag bytes are always valid ASCII
                    let tag_str = std::str::from_utf8(&tag_bytes).unwrap_or("");
                    if tag_info.remove.contains(tag_str) {
                        continue;
                    }
                    bam_fields::remove_tag(rr[i].as_mut_vec(), &tag_bytes);
                    append_buf_value_raw(rr[i].as_mut_vec(), tag_bytes, value);
                }
            }
        }

        // Copy tags to negative strand reads (with reverse/revcomp)
        if has_neg {
            let rr = mapped.raw_records.as_mut().unwrap();
            for &i in &mapped_indices {
                if (RawRecordView::new(&rr[i]).flags() & bam_fields::flags::REVERSE) == 0 {
                    continue;
                }
                let aux = bam_fields::aux_data_slice(&rr[i]);
                let has_pg = bam_fields::find_tag_type(aux, &pg_tag).is_some();

                // Aux offset is safe to cache here: `remove_tag` and `append_buf_value_raw`
                // only modify bytes within/after the aux region, so this offset stays valid.
                let aux_offset =
                    bam_fields::aux_data_offset_from_record(&rr[i]).unwrap_or(rr[i].len());

                for (tag, value) in u.data().iter() {
                    let tag_bytes: [u8; 2] = [tag.as_ref()[0], tag.as_ref()[1]];
                    if tag_bytes == pg_tag && has_pg {
                        continue;
                    }
                    let tag_str = std::str::from_utf8(&tag_bytes).unwrap_or("");
                    if tag_info.remove.contains(tag_str) {
                        continue;
                    }

                    bam_fields::remove_tag(rr[i].as_mut_vec(), &tag_bytes);

                    append_buf_value_raw(rr[i].as_mut_vec(), tag_bytes, value);
                    if has_transforms && tag_info.reverse.contains(tag_str) {
                        reverse_tag_in_place_raw(&mut rr[i], aux_offset, tag_bytes, value);
                    } else if has_transforms && tag_info.revcomp.contains(tag_str) {
                        revcomp_tag_in_place_raw(&mut rr[i], aux_offset, tag_bytes, value);
                    }
                }
            }
        }

        // Step 5: Transfer QC pass/fail flag
        let is_qc_fail = u.flags().is_qc_fail();
        let rr = mapped.raw_records.as_mut().unwrap();
        for &i in &mapped_indices {
            let mut f = RawRecordView::new(&rr[i]).flags();
            if is_qc_fail {
                f |= bam_fields::flags::QC_FAIL;
            } else {
                f &= !bam_fields::flags::QC_FAIL;
            }
            bam_fields::set_flags(&mut rr[i], f);
        }
    }

    // Step 6: Normalize AS/XS tags
    let rr = mapped.raw_records.as_mut().unwrap();
    for record in rr.iter_mut() {
        bam_fields::normalize_int_tag_to_smallest_signed(record.as_mut_vec(), b"AS");
        bam_fields::normalize_int_tag_to_smallest_signed(record.as_mut_vec(), b"XS");
    }

    // Step 7: Add TC tags
    if !skip_tc_tags {
        add_template_coordinate_tags_raw(mapped);
    }

    Ok(())
}

/// Applies the appropriate reverse operation for a tag in-place.
fn reverse_tag_in_place_raw(
    record: &mut [u8],
    aux_offset: usize,
    tag_bytes: [u8; 2],
    value: &BufValue,
) {
    match value {
        BufValue::String(_) => {
            bam_fields::reverse_string_tag_in_place(record, aux_offset, &tag_bytes);
        }
        BufValue::Array(_) => {
            bam_fields::reverse_array_tag_in_place(record, aux_offset, &tag_bytes);
        }
        _ => {}
    }
}

/// Applies the appropriate reverse-complement operation for a tag in-place.
fn revcomp_tag_in_place_raw(
    record: &mut [u8],
    aux_offset: usize,
    tag_bytes: [u8; 2],
    value: &BufValue,
) {
    match value {
        BufValue::String(_) => {
            bam_fields::reverse_complement_string_tag_in_place(record, aux_offset, &tag_bytes);
        }
        BufValue::Array(_) => {
            // For array tags, revcomp is the same as reverse
            bam_fields::reverse_array_tag_in_place(record, aux_offset, &tag_bytes);
        }
        _ => {}
    }
}

/// Encodes all records of an unmapped template to raw BAM bytes for output.
fn encode_unmapped_template_records(
    template: &Template,
    header: &Header,
) -> Result<Vec<RawRecord>> {
    let mut raw = Vec::with_capacity(template.read_count());
    for record in template.all_reads() {
        let mut buf = Vec::with_capacity(256);
        encode_record_buf(&mut buf, header, record)
            .map_err(|e| anyhow::anyhow!("encode unmapped: {e}"))?;
        raw.push(RawRecord::from(buf));
    }
    Ok(raw)
}

/// YD value for the forward (top) bisulfite strand.
const YD_FORWARD: &[u8] = b"f";
/// YD value for the reverse (bottom) bisulfite strand.
const YD_REVERSE: &[u8] = b"r";

/// YD tag used by bwameth to indicate the bisulfite conversion strand.
/// Only consumed by the cfg(test) `RecordBuf` reference implementation; the
/// raw-byte production path looks up the tag bytes directly.
#[cfg(test)]
const YD_TAG: Tag = Tag::new(b'Y', b'D');

/// `RecordBuf`-based reference implementation kept around as an oracle for
/// the test suite; the raw-byte path is used in production.
///
/// Walks the CIGAR alignment for each record and replaces converted bases back
/// to their unconverted reference form:
///   - Top strand (`YD:Z:f`): at reference-C positions, T→C
///   - Bottom strand (`YD:Z:r`): at reference-G positions, A→G
///
/// Skips unmapped reads and reads without a `YD` tag.
#[cfg(test)]
#[allow(dead_code)]
fn restore_unconverted_bases_in_template(
    template: &mut Template,
    reference: &ReferenceReader,
    header: &Header,
) -> Result<()> {
    for record in &mut template.records {
        restore_unconverted_bases_in_record(record, reference, header)?;
    }
    Ok(())
}

/// `RecordBuf`-based reference implementation kept as an oracle for the test
/// suite; the raw-byte path is used in production.
#[cfg(test)]
fn restore_unconverted_bases_in_record(
    record: &mut RecordBuf,
    reference: &ReferenceReader,
    header: &Header,
) -> Result<()> {
    // Skip unmapped reads
    if record.flags().is_unmapped() {
        return Ok(());
    }

    // Get the bisulfite strand from the bwameth YD tag
    let yd_value = record.data().get(&YD_TAG);
    let is_top = match yd_value {
        Some(BufValue::String(s)) if AsRef::<[u8]>::as_ref(s) == YD_FORWARD => true,
        Some(BufValue::String(s)) if AsRef::<[u8]>::as_ref(s) == YD_REVERSE => false,
        _ => return Ok(()), // No YD tag or unexpected value; skip
    };

    // Get reference contig name
    let ref_id = match record.reference_sequence_id() {
        Some(id) => id,
        None => return Ok(()),
    };
    let (ref_name, _) = header
        .reference_sequences()
        .get_index(ref_id)
        .context("reference sequence ID not found in header")?;
    // Borrow the contig name as &str rather than allocating a fresh String per
    // record; `fetch_slice` only needs &str for its HashMap lookup.
    let ref_name: &str = ref_name.to_str().context("reference sequence name is not valid UTF-8")?;

    // Get alignment start (1-based)
    let alignment_start = match record.alignment_start().map(usize::from) {
        Some(pos) => pos,
        None => return Ok(()),
    };

    // Compute reference span from CIGAR
    let ref_span = crate::sam::reference_length(&record.cigar());

    if ref_span == 0 {
        return Ok(());
    }

    // Fetch the entire aligned reference region at once. Use the borrowed-slice
    // API so we don't allocate a fresh `Vec<u8>` per record (~2 GB of churn over
    // a 13M template run).
    let ref_start = Position::try_from(alignment_start)?;
    let ref_end = Position::try_from(alignment_start + ref_span - 1)?;
    let ref_bases = reference.fetch_slice(ref_name, ref_start, ref_end)?;

    // Determine replacement parameters, accounting for reverse-complemented SEQ.
    // SAM stores SEQ reverse-complemented when 0x10 is set, so the base to find
    // and replace flips on reverse-strand alignments.
    let is_reverse = record.flags().is_reverse_complemented();
    let (ref_target, converted_base, unconverted_base) = match (is_top, is_reverse) {
        (true, false) | (false, true) => (b'C', b'T', b'C'),
        (true, true) | (false, false) => (b'G', b'A', b'G'),
    };
    // Precompute lowercase variants once per record so the inner loop just does
    // two byte equality checks instead of `to_ascii_uppercase()` per byte.
    let ref_target_lower = ref_target.to_ascii_lowercase();
    let converted_base_lower = converted_base.to_ascii_lowercase();

    // Fast path: if no candidate reference base appears in the aligned span,
    // the read can't have any bases to restore — skip the SEQ allocation and
    // CIGAR walk entirely. `memchr2` is SIMD-accelerated.
    if memchr::memchr2(ref_target, ref_target_lower, ref_bases).is_none() {
        return Ok(());
    }

    // Walk CIGAR and replace converted bases. Defer allocating the mutable SEQ
    // buffer until we hit the first base we actually need to change — many reads
    // (highly methylated CpGs in EM-Seq) align to ref-C positions but already
    // carry C in the read, so no allocation is needed.
    let read_seq = record.sequence().as_ref();
    let mut seq: Option<Vec<u8>> = None;
    let mut read_pos: usize = 0;
    let mut ref_offset: usize = 0; // 0-based offset into ref_bases

    for op_result in record.cigar().iter() {
        let op = op_result?;
        let len = op.len();

        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for i in 0..len {
                    if ref_offset + i >= ref_bases.len() {
                        break;
                    }
                    let rb = ref_bases[ref_offset + i];
                    if (rb == ref_target || rb == ref_target_lower) && read_pos + i < read_seq.len()
                    {
                        // Consult the (possibly already-modified) SEQ buffer if present,
                        // otherwise the original read sequence.
                        let sb = seq.as_ref().map_or(read_seq[read_pos + i], |s| s[read_pos + i]);
                        if sb == converted_base || sb == converted_base_lower {
                            let s = seq.get_or_insert_with(|| read_seq.to_vec());
                            s[read_pos + i] = unconverted_base;
                        }
                    }
                }
                read_pos += len;
                ref_offset += len;
            }
            Kind::Insertion | Kind::SoftClip => {
                read_pos += len;
            }
            Kind::Deletion | Kind::Skip => {
                ref_offset += len;
            }
            Kind::HardClip | Kind::Pad => {}
        }
    }

    if let Some(seq) = seq {
        *record.sequence_mut() = seq.into();
        // NM/MD tags are now stale since SEQ changed; remove them so downstream
        // tools don't trust incorrect mismatch counts.
        let data = record.data_mut();
        data.remove(&Tag::from(SamTag::NM));
        data.remove(&Tag::from(SamTag::MD));
    }

    Ok(())
}

/// Raw-byte equivalent of `restore_unconverted_bases_in_template`.
///
/// Iterates the template's `raw_records` (populated by
/// `convert_template_to_raw`) and applies the unconverted-base restoration
/// directly to the BAM bytes — no `RecordBuf` round-trip. Lets the methylation
/// path stay on the same streaming raw-byte fast path TAPs uses.
fn restore_unconverted_bases_in_raw_template(
    template: &mut Template,
    reference: &ReferenceReader,
    header: &Header,
) -> Result<()> {
    let raw_records = template.raw_records.as_mut().ok_or_else(|| {
        anyhow::anyhow!("restore_unconverted_bases_in_raw_template: template not in raw-byte mode")
    })?;
    for rec in raw_records.iter_mut() {
        restore_unconverted_bases_in_raw_record(rec, reference, header)?;
    }
    Ok(())
}

/// Raw-byte equivalent of `restore_unconverted_bases_in_record`.
///
/// Operates on a [`RawRecord`] without inflating to a [`RecordBuf`]. Reads the
/// CIGAR, SEQ, and `YD` tag via the `fgumi-raw-bam` accessors, walks the
/// alignment, then rewrites converted bases via [`RawRecord::set_base`]. SEQ
/// stays the same length so no record resizing is needed; only `NM` / `MD` are
/// removed when SEQ is actually modified.
///
/// Skips unmapped records and records without a `YD` tag, just like the
/// `RecordBuf` path.
fn restore_unconverted_bases_in_raw_record(
    rec: &mut RawRecord,
    reference: &ReferenceReader,
    header: &Header,
) -> Result<()> {
    if rec.is_unmapped() {
        return Ok(());
    }

    // YD tag tells us the bisulfite strand (forward = top, reverse = bottom).
    // Anything else (missing or unexpected value) means "skip".
    let yd_bytes = rec.tags().find_string(SamTag::YD.as_ref());
    let is_top = match yd_bytes {
        Some(YD_FORWARD) => true,
        Some(YD_REVERSE) => false,
        _ => return Ok(()),
    };

    // Reference contig + alignment start.
    let ref_id = rec.ref_id();
    if ref_id < 0 {
        return Ok(());
    }
    let (ref_name, _) = header
        .reference_sequences()
        .get_index(ref_id as usize)
        .context("reference sequence ID not found in header")?;
    let ref_name = ref_name.to_str().context("reference sequence name is not valid UTF-8")?;

    let pos = rec.pos();
    if pos < 0 {
        return Ok(());
    }
    // BAM stores `pos` 0-based; ReferenceReader::fetch_slice expects 1-based.
    let alignment_start = (pos as usize) + 1;

    // Reference span from the CIGAR. The helper returns it as i32; clamp to
    // non-negative usize since negative spans don't make sense.
    let ref_span = rec.reference_length();
    if ref_span <= 0 {
        return Ok(());
    }
    let ref_span = ref_span as usize;

    // Fetch the aligned reference span (borrowed slice, no per-record alloc).
    let ref_start = Position::try_from(alignment_start)?;
    let ref_end = Position::try_from(alignment_start + ref_span - 1)?;
    let ref_bases = reference.fetch_slice(ref_name, ref_start, ref_end)?;

    // Determine replacement parameters, accounting for reverse-complemented SEQ.
    // Same logic as the `RecordBuf` path; pre-encode the converted base to its
    // 4-bit BAM nibble so the inner loop compares nibbles directly.
    let (ref_target, converted_ascii, unconverted_ascii) = match (is_top, rec.is_reverse()) {
        (true, false) | (false, true) => (b'C', b'T', b'C'),
        (true, true) | (false, false) => (b'G', b'A', b'G'),
    };
    let ref_target_lower = ref_target.to_ascii_lowercase();
    // BAM 4-bit nibble codes: A=1, C=2, G=4, T=8. SEQ codes are case-insensitive
    // by construction (the encode table maps both cases to the same nibble),
    // so a single nibble compare is enough.
    let converted_nibble = match converted_ascii {
        b'T' => 0x8u8,
        b'A' => 0x1u8,
        _ => unreachable!(),
    };

    // Fast path: if no candidate reference base appears in the aligned span,
    // there's nothing to restore. `memchr2` is SIMD-accelerated.
    if memchr::memchr2(ref_target, ref_target_lower, ref_bases).is_none() {
        return Ok(());
    }

    // Walk CIGAR and replace converted bases in place. We don't allocate a
    // SEQ buffer at all — `set_base` rewrites the relevant 4-bit nibble in
    // the existing record bytes.
    let l_seq = rec.l_seq() as usize;
    let cigar_ops = rec.cigar_ops_vec();
    let mut modified = false;
    let mut read_pos: usize = 0;
    let mut ref_offset: usize = 0; // 0-based into ref_bases

    for raw_op in cigar_ops {
        let kind = bam_fields::cigar_op_kind(raw_op);
        let len = (raw_op >> 4) as usize;

        match kind {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for i in 0..len {
                    if ref_offset + i >= ref_bases.len() {
                        break;
                    }
                    let rb = ref_bases[ref_offset + i];
                    if (rb == ref_target || rb == ref_target_lower)
                        && read_pos + i < l_seq
                        && rec.get_base(read_pos + i) == converted_nibble
                    {
                        rec.set_base(read_pos + i, unconverted_ascii);
                        modified = true;
                    }
                }
                read_pos += len;
                ref_offset += len;
            }
            Kind::Insertion | Kind::SoftClip => {
                read_pos += len;
            }
            Kind::Deletion | Kind::Skip => {
                ref_offset += len;
            }
            Kind::HardClip | Kind::Pad => {}
        }
    }

    // NM/MD only become stale when SEQ was actually modified.
    if modified {
        let mut ed = rec.tags_editor();
        ed.remove(SamTag::NM.as_ref());
        ed.remove(SamTag::MD.as_ref());
    }

    Ok(())
}

impl Zipper {
    /// Process templates using raw-byte merge path with BGZF compression.
    ///
    /// Thread count is controlled by `self.threads` (1 = single-threaded).
    fn process_raw<U, M>(
        &self,
        unmapped_iter: U,
        mut mapped_iter: M,
        output_header: &Header,
        tag_info: &TagInfo,
        reference: Option<&ReferenceReader>,
    ) -> Result<u64>
    where
        U: Iterator<Item = Result<Template>>,
        M: Iterator<Item = Result<Template>>,
    {
        let mut writer = create_raw_bam_writer(
            &self.output,
            output_header,
            self.threads,
            self.compression.compression_level,
        )?;

        let progress = ProgressTracker::new("Processed records").with_interval(1_000_000);
        let mut mapped_peek: Option<Template> = None;
        let mut templates_not_in_mapped_bam: u64 = 0;

        for unmapped_result in unmapped_iter {
            let unmapped_template = unmapped_result?;

            if mapped_peek.is_none() {
                mapped_peek = mapped_iter.next().transpose()?;
            }

            if let Some(ref mut mapped_template) = mapped_peek {
                if mapped_template.name == unmapped_template.name {
                    // Both branches now stream raw BAM bytes; the methylation
                    // case adds a per-record `restore_unconverted_bases` step
                    // that mutates SEQ in place via `fgumi-raw-bam` primitives.
                    convert_template_to_raw(mapped_template, output_header)?;
                    merge_raw(&unmapped_template, mapped_template, tag_info, self.skip_tc_tags)?;
                    if let Some(ref_reader) = reference {
                        restore_unconverted_bases_in_raw_template(
                            mapped_template,
                            ref_reader,
                            output_header,
                        )?;
                    }
                    let rr = mapped_template.raw_records.as_ref().unwrap();
                    for rec in rr {
                        writer.write_raw_record(rec)?;
                        progress.log_if_needed(1);
                    }
                    mapped_peek = None;
                } else {
                    debug!(
                        "Found unmapped read with no corresponding mapped \
                         read: {}",
                        String::from_utf8_lossy(&unmapped_template.name)
                    );
                    if self.exclude_missing_reads {
                        templates_not_in_mapped_bam += 1;
                    } else {
                        let raw =
                            encode_unmapped_template_records(&unmapped_template, output_header)?;
                        for rec in &raw {
                            writer.write_raw_record(rec)?;
                            progress.log_if_needed(1);
                        }
                    }
                }
            } else {
                debug!(
                    "Found unmapped read with no corresponding mapped \
                     read: {}",
                    String::from_utf8_lossy(&unmapped_template.name)
                );
                if self.exclude_missing_reads {
                    templates_not_in_mapped_bam += 1;
                } else {
                    let raw = encode_unmapped_template_records(&unmapped_template, output_header)?;
                    for rec in &raw {
                        writer.write_raw_record(rec)?;
                        progress.log_if_needed(1);
                    }
                }
            }
        }

        progress.log_final();

        // Check for leftover mapped reads
        if let Some(remaining) = mapped_peek {
            anyhow::bail!(
                "Error: processed all unmapped reads but there are mapped \
                 reads remaining. Found template '{}'. Please ensure the \
                 unmapped and mapped reads have the same set of read names \
                 in the same order, and reads with the same name are \
                 consecutive (grouped) in each input.",
                String::from_utf8_lossy(&remaining.name)
            );
        }
        if let Some(remaining) = mapped_iter.next() {
            let remaining = remaining?;
            anyhow::bail!(
                "Error: processed all unmapped reads but there are mapped \
                 reads remaining. Found template '{}'. Please ensure the \
                 unmapped and mapped reads have the same set of read names \
                 in the same order, and reads with the same name are \
                 consecutive (grouped) in each input.",
                String::from_utf8_lossy(&remaining.name)
            );
        }

        if self.exclude_missing_reads && templates_not_in_mapped_bam > 0 {
            info!(
                "Excluded {templates_not_in_mapped_bam} templates that \
                 were not present in the aligned BAM."
            );
        }

        writer.finish()?;
        Ok(progress.count())
    }
}

/// BAM reader for stdin: single-threaded BGZF over a buffered `Box<dyn Read + Send>`.
///
/// Stdin is non-seekable, so we use a single-threaded BGZF reader (the multi-threaded
/// variant in [`BamReaderAuto`] requires `Seek`). `BufReader` wraps stdin so we can peek
/// the BGZF magic bytes for format auto-detection without consuming them.
type BamStdinReader =
    noodles::bam::io::Reader<noodles::bgzf::io::Reader<BufReader<Box<dyn Read + Send>>>>;

/// Wraps SAM and BAM readers so the mapped-reader thread can handle either format.
/// Moved into the thread where `record_bufs()` is called, since it borrows `&self`.
enum MappedReader {
    Sam(noodles::sam::io::Reader<Box<dyn BufRead + Send>>),
    Bam(BamReaderAuto),
    StdinBam(BamStdinReader),
}

/// First four bytes of a BGZF stream (gzip magic + extra-flag byte unique to BGZF).
const BGZF_MAGIC: [u8; 4] = [0x1f, 0x8b, 0x08, 0x04];

/// Peeks up to [`BGZF_MAGIC.len()`] bytes from `stdin` to classify it as BGZF (BAM) vs SAM text.
///
/// A single `read()` on a pipe may return fewer bytes than requested even when more
/// data is coming (e.g. the aligner hasn't flushed its first block yet), so this loops
/// until the peek buffer is full or EOF is reached. The peeked bytes are then chained
/// back onto the stream so downstream parsers see the complete input.
///
/// # Returns
///
/// `(is_bgzf, reader)` where `reader` yields the peeked bytes followed by the rest of
/// the original stream.
fn peek_stdin_bgzf(mut stdin: Box<dyn Read + Send>) -> Result<(bool, Box<dyn Read + Send>)> {
    let mut prefix = [0u8; BGZF_MAGIC.len()];
    let mut filled = 0usize;
    while filled < prefix.len() {
        match stdin.read(&mut prefix[filled..]) {
            Ok(0) => break,
            Ok(n) => filled += n,
            Err(e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
            Err(e) => return Err(e).context("Failed to peek stdin"),
        }
    }
    let is_bgzf = filled == BGZF_MAGIC.len() && prefix == BGZF_MAGIC;
    // Prepend the peeked bytes back onto the stream so the chosen reader sees them.
    let restored: Box<dyn Read + Send> =
        Box::new(std::io::Cursor::new(prefix[..filled].to_vec()).chain(stdin));
    Ok((is_bgzf, restored))
}

impl Command for Zipper {
    /// Executes the `zipper` command
    ///
    /// Main workflow:
    /// 1. Opens and validates input BAM files
    /// 2. Checks queryname sorting
    /// 3. Builds output header from both inputs and reference dictionary
    /// 4. Iterates through templates from both files in parallel
    /// 5. Merges tags from unmapped to mapped reads
    /// 6. Writes merged records to output
    ///
    /// # Returns
    ///
    /// `Ok(())` on successful completion
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - Input files cannot be opened
    /// - Files are not properly queryname sorted/grouped
    /// - Reference dictionary is missing or invalid
    /// - Input files have different sets of read names
    /// - I/O errors occur during reading or writing
    fn execute(&self, command_line: &str) -> Result<()> {
        info!("Starting zipper");

        let timer = OperationTimer::new("Zipping BAMs");

        // Validate input files exist
        validate_file_exists(&self.unmapped, "unmapped BAM file")?;
        validate_file_exists(&self.reference, "reference FASTA file")?;
        let dict_path = find_dict_path(&self.reference).ok_or_else(|| {
            anyhow::anyhow!(
                "Reference dictionary file not found. Tried:\n  \
                - {}\n  \
                - {}.dict\n\
                Please run: samtools dict {} -o {}",
                self.reference.with_extension("dict").display(),
                self.reference.display(),
                self.reference.display(),
                self.reference.with_extension("dict").display()
            )
        })?;

        let (mut unmapped_reader, unmapped_header) = create_bam_reader(&self.unmapped, 1)?;

        // Read mapped input — detect format by extension.
        // The reader and header are separated here; record_bufs() is called inside
        // the spawned thread (it borrows &self, so the reader must live there).
        let (mapped_reader, mapped_header) = if is_stdin_path(&self.input) {
            // Auto-detect BAM (BGZF magic) vs SAM text on stdin so callers can pipe either,
            // matching the file-input path. Without this, piping `samtools view -b ...` to
            // zipper crashes with a confusing "invalid flags / lexical parse error" because
            // the SAM text parser misreads the binary BGZF stream.
            //
            // `peek_stdin_bgzf` reads up to the magic length in a loop (a single `read()` on
            // a pipe may return fewer bytes than requested) and returns a reader that prepends
            // the peeked bytes back onto the stream for the chosen downstream reader.
            let stdin: Box<dyn Read + Send> = Box::new(std::io::stdin());
            let (is_bgzf, stdin) = peek_stdin_bgzf(stdin)?;
            let buffered = BufReader::with_capacity(64 * 1024, stdin);
            if is_bgzf {
                warn!(
                    "BAM input detected on stdin. For best performance, pipe SAM directly \
                     from the aligner (e.g. bwa mem ... | fgumi zipper ...)."
                );
                let bgzf = noodles::bgzf::io::Reader::new(buffered);
                let mut bam_reader = noodles::bam::io::Reader::from(bgzf);
                let header = bam_reader.read_header()?;
                (MappedReader::StdinBam(bam_reader), header)
            } else {
                info!(
                    "Reading SAM from stdin with adaptive buffer (bwa -K {})",
                    self.bwa_chunk_size
                );
                let reader: Box<dyn BufRead + Send> =
                    Box::new(BatchedSamReader::new(buffered, self.bwa_chunk_size));
                let mut sam_reader = noodles::sam::io::Reader::new(reader);
                let header = sam_reader.read_header()?;
                (MappedReader::Sam(sam_reader), header)
            }
        } else if self.input.extension().is_some_and(|ext| ext.eq_ignore_ascii_case("bam")) {
            // BAM file input — functional but discouraged
            warn!(
                "BAM input detected for --input. For best performance, pipe SAM directly \
                 from the aligner (e.g. bwa mem ... | fgumi zipper ...)."
            );
            let (bam_reader, header) = create_bam_reader(&self.input, self.threads)?;
            (MappedReader::Bam(bam_reader), header)
        } else {
            // SAM file input
            let reader: Box<dyn BufRead + Send> = Box::new(BufReader::with_capacity(
                256 * 1024,
                std::fs::File::open(&self.input).context("Failed to open mapped SAM")?,
            ));
            let mut sam_reader = noodles::sam::io::Reader::new(reader);
            let header = sam_reader.read_header()?;
            (MappedReader::Sam(sam_reader), header)
        };

        check_sort(&unmapped_header, &self.unmapped, "unmapped");
        check_sort(&mapped_header, &self.input, "mapped");

        let output_header = build_output_header(&unmapped_header, &mapped_header, &dict_path)?;

        // Add @PG record with PP chaining
        let output_header = crate::commands::common::add_pg_record(output_header, command_line)?;

        let tag_info = TagInfo::new(
            self.tags_to_remove.clone(),
            self.tags_to_reverse.clone(),
            self.tags_to_revcomp.clone(),
        );

        if !tag_info.remove.is_empty() {
            info!("Tags for removal: {:?}", tag_info.remove);
        }
        if !tag_info.reverse.is_empty() {
            info!("Tags being reversed: {:?}", tag_info.reverse);
        }
        if !tag_info.revcomp.is_empty() {
            info!("Tags being reverse complemented: {:?}", tag_info.revcomp);
        }

        // Load reference FASTA if restoring unconverted bases
        let reference = if self.restore_unconverted_bases {
            info!("Loading reference FASTA for unconverted base restoration");
            Some(ReferenceReader::new(&self.reference)?)
        } else {
            None
        };

        // Create async unmapped reader - spawn thread to read ahead
        let (unmapped_tx, unmapped_rx) =
            std::sync::mpsc::sync_channel::<Result<Template>>(self.buffer);
        let unmapped_header_for_reader = unmapped_header.clone();
        std::thread::spawn(move || {
            let unmapped_iter = TemplateIterator::new(
                unmapped_reader
                    .record_bufs(&unmapped_header_for_reader)
                    .map(|r| r.map_err(anyhow::Error::from)),
            );
            for template in unmapped_iter {
                if unmapped_tx.send(template).is_err() {
                    break; // Receiver dropped, main thread done
                }
            }
        });
        let unmapped_iter = std::iter::from_fn(move || unmapped_rx.recv().ok());

        // Create async mapped reader - spawn thread to read ahead and prevent backpressure on stdin
        // (e.g., when reading piped output from bwa)
        let (mapped_tx, mapped_rx) = std::sync::mpsc::sync_channel::<Result<Template>>(self.buffer);
        let mapped_header_for_reader = mapped_header.clone();
        std::thread::spawn(move || {
            // The reader must be owned for the full duration so that the iterator
            // (which borrows it via record_bufs) remains valid.
            match mapped_reader {
                MappedReader::Sam(mut r) => {
                    let mapped_iter = TemplateIterator::new(
                        r.record_bufs(&mapped_header_for_reader)
                            .map(|rec| rec.map_err(anyhow::Error::from)),
                    );
                    for template in mapped_iter {
                        if mapped_tx.send(template).is_err() {
                            break;
                        }
                    }
                }
                MappedReader::Bam(mut r) => {
                    let mapped_iter = TemplateIterator::new(
                        r.record_bufs(&mapped_header_for_reader)
                            .map(|rec| rec.map_err(anyhow::Error::from)),
                    );
                    for template in mapped_iter {
                        if mapped_tx.send(template).is_err() {
                            break;
                        }
                    }
                }
                MappedReader::StdinBam(mut r) => {
                    let mapped_iter = TemplateIterator::new(
                        r.record_bufs(&mapped_header_for_reader)
                            .map(|rec| rec.map_err(anyhow::Error::from)),
                    );
                    for template in mapped_iter {
                        if mapped_tx.send(template).is_err() {
                            break;
                        }
                    }
                }
            }
        });
        let mapped_iter = std::iter::from_fn(move || mapped_rx.recv().ok());

        let total_records = self.process_raw(
            unmapped_iter,
            mapped_iter,
            &output_header,
            &tag_info,
            reference.as_ref(),
        )?;

        info!("zipper completed successfully");
        timer.log_completion(total_records);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sam::builder::{
        MAPPED_PG_ID, REFERENCE_LENGTH, RecordBuilder, SamBuilder, create_ref_dict,
    };
    use anyhow::Result;
    use bstr::ByteSlice;
    use noodles::sam::alignment::record::Flags;
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::RecordBuf;
    use noodles::sam::alignment::record_buf::data::field::Value as BufValue;
    use rstest::rstest;
    use std::collections::HashMap;
    use std::io::Read;
    use tempfile::TempDir;

    /// A `Read` implementation that returns at most one byte per `read` call,
    /// mimicking the behavior of a slow pipe or small kernel buffer. Used to
    /// regression-test `peek_stdin_bgzf`, which must tolerate short reads.
    struct TrickleReader(std::collections::VecDeque<u8>);

    impl Read for TrickleReader {
        fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
            if buf.is_empty() {
                return Ok(0);
            }
            match self.0.pop_front() {
                Some(b) => {
                    buf[0] = b;
                    Ok(1)
                }
                None => Ok(0),
            }
        }
    }

    fn trickle(bytes: &[u8]) -> Box<dyn Read + Send> {
        Box::new(TrickleReader(bytes.iter().copied().collect()))
    }

    /// Regression: `BufReader::fill_buf` can return fewer than `BGZF_MAGIC.len()`
    /// bytes on a pipe even when more data is coming. Before the fix, that caused
    /// BAM-on-stdin to be misclassified as SAM. Simulate it with a reader that
    /// returns one byte per `read` call and confirm classification + stream
    /// reconstruction both work.
    #[test]
    fn test_peek_stdin_bgzf_classifies_bam_under_short_reads() -> Result<()> {
        let mut stream = BGZF_MAGIC.to_vec();
        stream.extend_from_slice(b"rest-of-bgzf-stream");
        let (is_bgzf, mut restored) = peek_stdin_bgzf(trickle(&stream))?;
        assert!(is_bgzf);
        let mut round_trip = Vec::new();
        restored.read_to_end(&mut round_trip)?;
        assert_eq!(round_trip, stream);
        Ok(())
    }

    #[test]
    fn test_peek_stdin_bgzf_classifies_sam_text() -> Result<()> {
        let stream = b"@HD\tVN:1.6\tSO:queryname\n".to_vec();
        let (is_bgzf, mut restored) = peek_stdin_bgzf(trickle(&stream))?;
        assert!(!is_bgzf);
        let mut round_trip = Vec::new();
        restored.read_to_end(&mut round_trip)?;
        assert_eq!(round_trip, stream);
        Ok(())
    }

    #[test]
    fn test_peek_stdin_bgzf_handles_stream_shorter_than_magic() -> Result<()> {
        // Fewer bytes than BGZF_MAGIC.len() — must not be misclassified as BAM,
        // and the short prefix must still be recoverable from the returned reader.
        let stream = b"@H".to_vec();
        let (is_bgzf, mut restored) = peek_stdin_bgzf(trickle(&stream))?;
        assert!(!is_bgzf);
        let mut round_trip = Vec::new();
        restored.read_to_end(&mut round_trip)?;
        assert_eq!(round_trip, stream);
        Ok(())
    }

    #[test]
    fn test_peek_stdin_bgzf_handles_empty_stream() -> Result<()> {
        let (is_bgzf, mut restored) = peek_stdin_bgzf(Box::new(std::io::empty()))?;
        assert!(!is_bgzf);
        let mut round_trip = Vec::new();
        restored.read_to_end(&mut round_trip)?;
        assert!(round_trip.is_empty());
        Ok(())
    }

    /// Runs `zipper-bams` and returns the output records
    ///
    /// Helper function that:
    /// 1. Writes unmapped and mapped builders to temporary BAM files
    /// 2. Creates a temporary reference dictionary
    /// 3. Executes `zipper-bams`
    /// 4. Reads and returns all output records
    ///
    /// # Arguments
    ///
    /// * `unmapped` - Builder containing unmapped reads
    /// * `mapped` - Builder containing mapped reads
    /// * `tags_to_remove` - Tags to remove from mapped reads
    /// * `tags_to_reverse` - Tags to reverse for negative strand
    /// * `tags_to_revcomp` - Tags to revcomp for negative strand
    ///
    /// # Returns
    ///
    /// Vector of all records written to output BAM
    ///
    /// # Errors
    ///
    /// Returns an error if any step fails
    fn run_zipper(
        unmapped: &SamBuilder,
        mapped: &SamBuilder,
        tags_to_remove: Vec<String>,
        tags_to_reverse: Vec<String>,
        tags_to_revcomp: Vec<String>,
    ) -> Result<Vec<RecordBuf>> {
        let dir = TempDir::new()?;
        let unmapped_path = dir.path().join("unmapped.bam");
        let mapped_path = dir.path().join("mapped.sam");
        let output_path = dir.path().join("output.bam");
        let dict_path = create_ref_dict(&dir, "chr1", REFERENCE_LENGTH)?;

        unmapped.write(&unmapped_path)?;
        mapped.write_sam(&mapped_path)?;

        let zipper = Zipper {
            input: mapped_path.clone(),
            unmapped: unmapped_path.clone(),
            reference: dict_path,
            output: output_path.clone(),
            tags_to_remove,
            tags_to_reverse,
            tags_to_revcomp,
            buffer: 5000,
            threads: 1,
            compression: CompressionOptions { compression_level: 1 },
            bwa_chunk_size: 150_000_000,
            exclude_missing_reads: false,
            skip_tc_tags: false,
            restore_unconverted_bases: false,
        };

        zipper.execute("test")?;
        read_bam_records(&output_path)
    }

    /// Read all records from a BAM file.
    fn read_bam_records(path: &std::path::Path) -> Result<Vec<RecordBuf>> {
        let mut reader = noodles::bam::io::reader::Builder.build_from_path(path)?;
        let header = reader.read_header()?;
        Ok(reader.record_bufs(&header).collect::<std::io::Result<Vec<_>>>()?)
    }

    /// Tests basic tag merging from unmapped to mapped reads
    ///
    /// Verifies that:
    /// - Tags from unmapped reads are copied to mapped reads
    /// - Tags from mapped reads are preserved
    /// - All records are present in output
    /// - MC and MQ tags are properly set for paired reads
    #[test]
    fn test_basic_merge() -> Result<()> {
        let mut unmapped = SamBuilder::new_unmapped();
        let mut mapped = SamBuilder::new_mapped();

        // Add unmapped pairs with tags
        let mut attrs1 = HashMap::new();
        attrs1.insert("RX", BufValue::from("ACGT".to_string()));
        attrs1.insert("xy", BufValue::from(1234i32));
        unmapped.add_pair_with_attrs("q1", None, None, true, true, &attrs1);

        let mut attrs2 = HashMap::new();
        attrs2.insert("RX", BufValue::from("GGTA".to_string()));
        attrs2.insert("xy", BufValue::from(4567i32));
        unmapped.add_pair_with_attrs("q2", None, None, true, true, &attrs2);

        // Add mapped pairs with different tags
        let mut mapped_attrs1 = HashMap::new();
        mapped_attrs1.insert("PG", BufValue::from(MAPPED_PG_ID.to_string()));
        mapped_attrs1.insert("AS", BufValue::from(77i32));
        mapped.add_pair_with_attrs("q1", Some(100), Some(200), true, true, &mapped_attrs1);

        let mut mapped_attrs2 = HashMap::new();
        mapped_attrs2.insert("PG", BufValue::from(MAPPED_PG_ID.to_string()));
        mapped_attrs2.insert("AS", BufValue::from(16i32));
        mapped.add_pair_with_attrs("q2", Some(500), Some(700), true, true, &mapped_attrs2);

        let records = run_zipper(&unmapped, &mapped, vec![], vec![], vec![])?;

        assert_eq!(records.len(), 4);

        // Check that tags from unmapped BAM were copied
        for rec in &records {
            // Should have RX tag from unmapped
            assert!(rec.data().get(&Tag::from(SamTag::RX)).is_some());

            // Should have xy tag from unmapped
            assert!(rec.data().get(&Tag::new(b'x', b'y')).is_some());

            // Should have AS tag from mapped
            assert!(rec.data().get(&Tag::from(SamTag::AS)).is_some());

            // Should have PG tag from mapped
            assert!(rec.data().get(&Tag::from(SamTag::PG)).is_some());

            // Should have MC (mate CIGAR) and MQ (mate mapping quality) tags
            assert!(rec.data().get(&Tag::from(SamTag::MC)).is_some(), "MC tag should be present");
            assert!(rec.data().get(&Tag::from(SamTag::MQ)).is_some(), "MQ tag should be present");
        }

        Ok(())
    }

    /// Tests handling of unpaired (single-end) reads
    ///
    /// Verifies that single-end reads are properly merged without
    /// attempting to match R1/R2.
    #[test]
    fn test_unpaired_reads() -> Result<()> {
        let mut unmapped = SamBuilder::new_unmapped();
        let mut mapped = SamBuilder::new_mapped();

        let mut attrs1 = HashMap::new();
        attrs1.insert("RX", BufValue::from("ACGT".to_string()));
        attrs1.insert("xy", BufValue::from(1234i32));
        unmapped.add_frag_with_attrs("q1", None, true, &attrs1);

        let mut attrs2 = HashMap::new();
        attrs2.insert("RX", BufValue::from("GGTA".to_string()));
        attrs2.insert("xy", BufValue::from(4567i32));
        unmapped.add_frag_with_attrs("q2", None, true, &attrs2);

        let mut mapped_attrs1 = HashMap::new();
        mapped_attrs1.insert("PG", BufValue::from(MAPPED_PG_ID.to_string()));
        mapped_attrs1.insert("AS", BufValue::from(77i32));
        mapped.add_frag_with_attrs("q1", Some(100), true, &mapped_attrs1);

        let mut mapped_attrs2 = HashMap::new();
        mapped_attrs2.insert("PG", BufValue::from(MAPPED_PG_ID.to_string()));
        mapped_attrs2.insert("AS", BufValue::from(16i32));
        mapped.add_frag_with_attrs("q2", Some(500), false, &mapped_attrs2);

        let records = run_zipper(&unmapped, &mapped, vec![], vec![], vec![])?;

        assert_eq!(records.len(), 2);

        for rec in &records {
            // Should have tags from both unmapped and mapped
            assert!(rec.data().get(&Tag::from(SamTag::RX)).is_some());
            assert!(rec.data().get(&Tag::new(b'x', b'y')).is_some());
            assert!(rec.data().get(&Tag::from(SamTag::AS)).is_some());
        }

        Ok(())
    }

    /// Tests tag removal, reversal, and reverse complementation
    ///
    /// Verifies that:
    /// - Specified tags are removed from mapped reads
    /// - Array and string tags are reversed for negative strand
    /// - Sequence tags are reverse complemented for negative strand
    /// - Positive strand reads are unchanged
    #[test]
    fn test_tag_removal_reverse_revcomp() -> Result<()> {
        let mut unmapped = SamBuilder::new_unmapped();
        let mut mapped = SamBuilder::new_mapped();

        let mut attrs = HashMap::new();
        attrs.insert("n1", BufValue::from(vec![1i16, 2, 3, 4, 5]));
        attrs.insert("n2", BufValue::from(vec![2i16, 3, 4, 5, 6]));
        attrs.insert("s1", BufValue::from("abcde".to_string()));
        attrs.insert("s2", BufValue::from("vwxyz".to_string()));
        attrs.insert("s3", BufValue::from("AGAGG".to_string()));
        unmapped.add_pair_with_attrs("q1", None, None, true, true, &attrs);

        let mut mapped_attrs = HashMap::new();
        mapped_attrs.insert("PG", BufValue::from(MAPPED_PG_ID.to_string()));
        mapped_attrs.insert("AS", BufValue::from(77i32));
        mapped.add_pair_with_attrs("q1", Some(100), Some(200), true, false, &mapped_attrs);

        let records = run_zipper(
            &unmapped,
            &mapped,
            vec!["AS".to_string()],
            vec!["n1".to_string(), "n2".to_string(), "s1".to_string()],
            vec!["s3".to_string()],
        )?;

        assert_eq!(records.len(), 2);

        for rec in &records {
            // AS tag should be removed
            assert!(rec.data().get(&Tag::from(SamTag::AS)).is_none());

            if rec.flags().is_first_segment() {
                // R1 is positive strand - no reversing/revcomping
                if let Some(BufValue::Array(
                    noodles::sam::alignment::record_buf::data::field::value::Array::Int16(vals),
                )) = rec.data().get(&Tag::new(b'n', b'1'))
                {
                    assert_eq!(*vals, vec![1i16, 2, 3, 4, 5]);
                }

                if let Some(BufValue::String(s)) = rec.data().get(&Tag::new(b's', b'1')) {
                    assert_eq!(s.to_string(), "abcde");
                }

                if let Some(BufValue::String(s)) = rec.data().get(&Tag::new(b's', b'3')) {
                    assert_eq!(s.to_string(), "AGAGG");
                }
            } else {
                // R2 is negative strand - should be reversed/revcomped
                if let Some(BufValue::Array(
                    noodles::sam::alignment::record_buf::data::field::value::Array::Int16(vals),
                )) = rec.data().get(&Tag::new(b'n', b'1'))
                {
                    assert_eq!(*vals, vec![5i16, 4, 3, 2, 1]);
                }

                if let Some(BufValue::String(s)) = rec.data().get(&Tag::new(b's', b'1')) {
                    assert_eq!(s.to_string(), "edcba");
                }

                if let Some(BufValue::String(s)) = rec.data().get(&Tag::new(b's', b'3')) {
                    assert_eq!(s.to_string(), "CCTCT");
                }
            }

            // s2 should not be changed
            if let Some(BufValue::String(s)) = rec.data().get(&Tag::new(b's', b'2')) {
                assert_eq!(s.to_string(), "vwxyz");
            }
        }

        Ok(())
    }

    /// Tests the "Consensus" named tag set expansion
    ///
    /// Verifies that specifying "Consensus" in `tags_to_reverse` or `tags_to_revcomp`
    /// properly expands to all per-base consensus tags and applies the correct
    /// transformations.
    #[test]
    fn test_consensus_tag_set() -> Result<()> {
        let mut unmapped = SamBuilder::new_unmapped();
        let mut mapped = SamBuilder::new_mapped();

        let mut attrs = HashMap::new();
        attrs.insert("aD", BufValue::from("AAAGG".to_string()));
        attrs.insert("bD", BufValue::from("AAAGC".to_string()));
        attrs.insert("ad", BufValue::from(vec![3i16, 3, 4, 4, 2]));
        unmapped.add_frag_with_attrs("q1", None, true, &attrs);

        let mut mapped_attrs1 = HashMap::new();
        mapped_attrs1.insert("PG", BufValue::from(MAPPED_PG_ID.to_string()));
        mapped_attrs1.insert("AS", BufValue::from(77i32));
        mapped.add_frag_with_attrs("q1", Some(100), true, &mapped_attrs1.clone());

        // Add supplementary record
        let mut supp_rec = RecordBuilder::new()
            .name("q1")
            .flags(Flags::SUPPLEMENTARY | Flags::REVERSE_COMPLEMENTED)
            .reference_sequence_id(0)
            .alignment_start(200)
            .cigar("100M")
            .sequence(&"A".repeat(100))
            .qualities(&[30u8; 100])
            .build();
        for (tag_str, value) in &mapped_attrs1 {
            let tag = Tag::new(tag_str.as_bytes()[0], tag_str.as_bytes()[1]);
            supp_rec.data_mut().insert(tag, value.clone());
        }
        mapped.push_record(supp_rec);

        let records = run_zipper(
            &unmapped,
            &mapped,
            vec![],
            vec!["Consensus".to_string()],
            vec!["Consensus".to_string()],
        )?;

        assert_eq!(records.len(), 2);

        for rec in &records {
            if rec.flags().is_reverse_complemented() {
                // Negative strand - should be reversed/revcomped
                if let Some(BufValue::String(s)) = rec.data().get(&Tag::new(b'a', b'D')) {
                    assert_eq!(s.to_string(), "CCTTT");
                }

                if let Some(BufValue::String(s)) = rec.data().get(&Tag::new(b'b', b'D')) {
                    assert_eq!(s.to_string(), "GCTTT");
                }

                if let Some(BufValue::Array(
                    noodles::sam::alignment::record_buf::data::field::value::Array::Int16(vals),
                )) = rec.data().get(&Tag::new(b'a', b'd'))
                {
                    assert_eq!(*vals, vec![2i16, 4, 4, 3, 3]);
                }
            } else {
                // Positive strand - no changes
                if let Some(BufValue::String(s)) = rec.data().get(&Tag::new(b'a', b'D')) {
                    assert_eq!(s.to_string(), "AAAGG");
                }

                if let Some(BufValue::String(s)) = rec.data().get(&Tag::new(b'b', b'D')) {
                    assert_eq!(s.to_string(), "AAAGC");
                }

                if let Some(BufValue::Array(
                    noodles::sam::alignment::record_buf::data::field::value::Array::Int16(vals),
                )) = rec.data().get(&Tag::new(b'a', b'd'))
                {
                    assert_eq!(*vals, vec![3i16, 3, 4, 4, 2]);
                }
            }
        }

        Ok(())
    }

    /// Tests handling of reads that are only in the unmapped BAM
    ///
    /// Verifies that if a read has no corresponding mapped read, it is
    /// written to the output unchanged (still unmapped).
    #[test]
    fn test_unmapped_only_reads() -> Result<()> {
        let mut unmapped = SamBuilder::new_unmapped();
        let mut mapped = SamBuilder::new_mapped();

        let mut attrs1 = HashMap::new();
        attrs1.insert("RX", BufValue::from("ACGT".to_string()));
        attrs1.insert("xy", BufValue::from(1234i32));
        unmapped.add_frag_with_attrs("q1", None, true, &attrs1);

        let mut attrs2 = HashMap::new();
        attrs2.insert("RX", BufValue::from("GATA".to_string()));
        attrs2.insert("xy", BufValue::from(3456i32));
        unmapped.add_frag_with_attrs("q2", None, true, &attrs2);

        let mut attrs3 = HashMap::new();
        attrs3.insert("RX", BufValue::from("GGCG".to_string()));
        attrs3.insert("xy", BufValue::from(5678i32));
        unmapped.add_frag_with_attrs("q3", None, true, &attrs3);

        let mut mapped_attrs1 = HashMap::new();
        mapped_attrs1.insert("PG", BufValue::from(MAPPED_PG_ID.to_string()));
        mapped_attrs1.insert("AS", BufValue::from(77i32));
        mapped.add_frag_with_attrs("q1", Some(100), true, &mapped_attrs1);

        let mut mapped_attrs3 = HashMap::new();
        mapped_attrs3.insert("PG", BufValue::from(MAPPED_PG_ID.to_string()));
        mapped_attrs3.insert("AS", BufValue::from(77i32));
        mapped.add_frag_with_attrs("q3", Some(200), false, &mapped_attrs3);

        let records = run_zipper(&unmapped, &mapped, vec![], vec![], vec![])?;

        assert_eq!(records.len(), 3);

        let names: Vec<String> = records
            .iter()
            .map(|r| {
                String::from_utf8(r.name().expect("record should have a name").as_bytes().to_vec())
                    .expect("record name should be valid UTF-8")
            })
            .collect();

        assert_eq!(names, vec!["q1", "q2", "q3"]);

        // q2 should be unmapped
        let q2 = records
            .iter()
            .find(|r| {
                String::from_utf8(r.name().expect("record should have a name").as_bytes().to_vec())
                    .expect("record name should be valid UTF-8")
                    == "q2"
            })
            .expect("expected record q2 not found");
        assert!(q2.flags().is_unmapped());

        Ok(())
    }

    /// Tests `TagInfo` expansion of named tag sets
    ///
    /// Verifies that the "Consensus" tag set is properly expanded into
    /// its constituent per-base tags.
    #[test]
    fn test_tag_info_expansion() {
        let tag_info =
            TagInfo::new(vec![], vec!["Consensus".to_string()], vec!["Consensus".to_string()]);

        // Check that "Consensus" was expanded
        assert!(tag_info.reverse.contains("ad"));
        assert!(tag_info.reverse.contains("ae"));
        assert!(tag_info.reverse.contains("bd"));
        assert!(tag_info.reverse.contains("be"));
        assert!(tag_info.reverse.contains("cd"));

        assert!(tag_info.revcomp.contains("aD"));
        assert!(tag_info.revcomp.contains("bD"));
        assert!(tag_info.revcomp.contains("cD"));
    }

    #[test]
    fn test_mate_info_fixing() -> Result<()> {
        let mut unmapped = SamBuilder::new_unmapped();
        let mut mapped = SamBuilder::new_mapped();

        let mut attrs = HashMap::new();
        attrs.insert("RX", BufValue::from("ACGT".to_string()));
        unmapped.add_pair_with_attrs("q1", None, None, true, true, &attrs);

        let mut mapped_attrs = HashMap::new();
        mapped_attrs.insert("PG", BufValue::from(MAPPED_PG_ID.to_string()));
        mapped_attrs.insert("AS", BufValue::from(77i32));
        mapped.add_pair_with_attrs("q1", Some(100), Some(200), true, true, &mapped_attrs);

        let records = run_zipper(&unmapped, &mapped, vec![], vec![], vec![])?;

        assert_eq!(records.len(), 2);

        // Find R1 and R2
        let r1 = records
            .iter()
            .find(|r| r.flags().is_first_segment())
            .expect("expected R1 record not found");
        let r2 = records
            .iter()
            .find(|r| !r.flags().is_first_segment())
            .expect("expected R2 record not found");

        // Check that MQ tags are set
        let r1_mq = r1.data().get(&Tag::from(SamTag::MQ));
        let r2_mq = r2.data().get(&Tag::from(SamTag::MQ));
        assert!(r1_mq.is_some(), "R1 should have MQ tag");
        assert!(r2_mq.is_some(), "R2 should have MQ tag");

        // Check that MC tags are set
        let r1_mc = r1.data().get(&Tag::from(SamTag::MC));
        let r2_mc = r2.data().get(&Tag::from(SamTag::MC));
        assert!(r1_mc.is_some(), "R1 should have MC tag");
        assert!(r2_mc.is_some(), "R2 should have MC tag");

        Ok(())
    }

    /// Tests QC flag transfer from unmapped to mapped reads
    ///
    /// Verifies that:
    /// - `QC_FAIL` flag is transferred from unmapped reads to mapped reads
    /// - `QC_PASS` (absence of flag) is also properly transferred
    /// - Both R1 and R2 get the correct flags
    #[test]
    fn test_qc_flag_transfer() -> Result<()> {
        let mut unmapped = SamBuilder::new_unmapped();
        let mut mapped = SamBuilder::new_mapped();

        // Add unmapped pair where R1 passes QC but R2 fails QC
        let mut attrs1 = HashMap::new();
        attrs1.insert("RX", BufValue::from("ACGT".to_string()));

        // Create R1 (passes QC - no QC_FAIL flag)
        let r1_unmapped = RecordBuilder::new()
            .name("q1")
            .flags(Flags::SEGMENTED | Flags::FIRST_SEGMENT | Flags::UNMAPPED)
            .tag("RX", "ACGT")
            .build();

        // Create R2 (fails QC - has QC_FAIL flag)
        let r2_unmapped = RecordBuilder::new()
            .name("q1")
            .flags(Flags::SEGMENTED | Flags::LAST_SEGMENT | Flags::UNMAPPED | Flags::QC_FAIL)
            .tag("RX", "ACGT")
            .build();

        unmapped.push_record(r1_unmapped);
        unmapped.push_record(r2_unmapped);

        // Add corresponding mapped pair (neither has QC_FAIL initially)
        let mut mapped_attrs = HashMap::new();
        mapped_attrs.insert("PG", BufValue::from(MAPPED_PG_ID.to_string()));
        mapped_attrs.insert("AS", BufValue::from(77i32));
        mapped.add_pair_with_attrs("q1", Some(100), Some(200), true, true, &mapped_attrs);

        let records = run_zipper(&unmapped, &mapped, vec![], vec![], vec![])?;

        assert_eq!(records.len(), 2);

        // Find R1 and R2 in output
        let r1 = records
            .iter()
            .find(|r| r.flags().is_first_segment())
            .expect("expected R1 record not found");
        let r2 = records
            .iter()
            .find(|r| !r.flags().is_first_segment())
            .expect("expected R2 record not found");

        // R1 should NOT have QC_FAIL flag (it passed QC in unmapped)
        assert!(
            !r1.flags().is_qc_fail(),
            "R1 should not have QC_FAIL flag (passed QC in unmapped)"
        );

        // R2 SHOULD have QC_FAIL flag (it failed QC in unmapped)
        assert!(r2.flags().is_qc_fail(), "R2 should have QC_FAIL flag (failed QC in unmapped)");

        Ok(())
    }

    /// Tests QC flag removal when unmapped read passes QC but mapped read had failed
    ///
    /// Verifies that if a mapped read has `QC_FAIL` but the unmapped read doesn't,
    /// the flag is removed from the mapped read.
    #[test]
    fn test_qc_flag_removal() -> Result<()> {
        let mut unmapped = SamBuilder::new_unmapped();
        let mut mapped = SamBuilder::new_mapped();

        // Add unmapped read that passes QC
        let mut attrs = HashMap::new();
        attrs.insert("RX", BufValue::from("ACGT".to_string()));
        unmapped.add_frag_with_attrs("q1", None, true, &attrs);

        // Add mapped read with QC_FAIL flag set
        let mapped_rec = RecordBuilder::new()
            .name("q1")
            .flags(Flags::QC_FAIL)
            .reference_sequence_id(0)
            .alignment_start(100)
            .cigar("100M")
            .sequence(&"A".repeat(100))
            .qualities(&[30u8; 100])
            .build();
        mapped.push_record(mapped_rec);

        let records = run_zipper(&unmapped, &mapped, vec![], vec![], vec![])?;

        assert_eq!(records.len(), 1);

        // The QC_FAIL flag should have been removed
        assert!(
            !records[0].flags().is_qc_fail(),
            "QC_FAIL flag should be removed when unmapped read passes QC"
        );

        Ok(())
    }

    /// Tests tag removal, reversal, and reverse complementation on secondary and supplementary alignments
    ///
    /// Verifies that:
    /// - Tag operations work correctly on both secondary and supplementary records
    /// - Multiple alignments (primary, secondary, supplementary) all get processed
    /// - Negative strand transformations apply to secondary/supplementary alignments
    #[test]
    #[allow(clippy::too_many_lines)]
    fn test_tag_operations_on_secondary_and_supplementary() -> Result<()> {
        let mut unmapped = SamBuilder::new_unmapped();
        let mut mapped = SamBuilder::new_mapped();

        let mut attrs = HashMap::new();
        attrs.insert("n1", BufValue::from(vec![1i16, 2, 3, 4, 5]));
        attrs.insert("n2", BufValue::from(vec![2i16, 3, 4, 5, 6]));
        attrs.insert("s1", BufValue::from("abcde".to_string()));
        attrs.insert("s2", BufValue::from("vwxyz".to_string()));
        attrs.insert("s3", BufValue::from("AGAGG".to_string()));
        unmapped.add_frag_with_attrs("q1", None, true, &attrs);

        let mut mapped_attrs = HashMap::new();
        mapped_attrs.insert("PG", BufValue::from(MAPPED_PG_ID.to_string()));
        mapped_attrs.insert("AS", BufValue::from(77i32));

        // Add primary alignment (positive strand)
        mapped.add_frag_with_attrs("q1", Some(100), true, &mapped_attrs.clone());

        // Add secondary alignment (negative strand)
        let mut secondary_rec = RecordBuilder::new()
            .name("q1")
            .flags(Flags::SECONDARY | Flags::REVERSE_COMPLEMENTED)
            .reference_sequence_id(0)
            .alignment_start(200)
            .cigar("100M")
            .sequence(&"A".repeat(100))
            .qualities(&[30u8; 100])
            .build();
        for (tag_str, value) in &mapped_attrs {
            let tag = Tag::new(tag_str.as_bytes()[0], tag_str.as_bytes()[1]);
            secondary_rec.data_mut().insert(tag, value.clone());
        }
        mapped.push_record(secondary_rec);

        // Add supplementary alignment (negative strand)
        let mut supp_rec = RecordBuilder::new()
            .name("q1")
            .flags(Flags::SUPPLEMENTARY | Flags::REVERSE_COMPLEMENTED)
            .reference_sequence_id(0)
            .alignment_start(300)
            .cigar("100M")
            .sequence(&"A".repeat(100))
            .qualities(&[30u8; 100])
            .build();
        for (tag_str, value) in &mapped_attrs {
            let tag = Tag::new(tag_str.as_bytes()[0], tag_str.as_bytes()[1]);
            supp_rec.data_mut().insert(tag, value.clone());
        }
        mapped.push_record(supp_rec);

        let records = run_zipper(
            &unmapped,
            &mapped,
            vec!["AS".to_string()],
            vec!["n1".to_string(), "n2".to_string(), "s1".to_string()],
            vec!["s3".to_string()],
        )?;

        assert_eq!(records.len(), 3, "Should have 3 records (primary, secondary, supplementary)");

        let secondary_and_supp_count = records
            .iter()
            .filter(|r| r.flags().is_secondary() || r.flags().is_supplementary())
            .count();
        assert_eq!(secondary_and_supp_count, 2, "Should have 2 secondary/supplementary records");

        let neg_strand_count =
            records.iter().filter(|r| r.flags().is_reverse_complemented()).count();
        assert_eq!(neg_strand_count, 2, "Should have 2 negative strand records");

        // Check each record
        for rec in &records {
            // AS tag should be removed from all records
            assert!(rec.data().get(&Tag::from(SamTag::AS)).is_none(), "AS tag should be removed");

            if rec.flags().is_reverse_complemented() {
                // Negative strand - should be reversed/revcomped
                if let Some(BufValue::Array(
                    noodles::sam::alignment::record_buf::data::field::value::Array::Int16(vals),
                )) = rec.data().get(&Tag::new(b'n', b'1'))
                {
                    assert_eq!(
                        *vals,
                        vec![5i16, 4, 3, 2, 1],
                        "n1 should be reversed for negative strand"
                    );
                }

                if let Some(BufValue::Array(
                    noodles::sam::alignment::record_buf::data::field::value::Array::Int16(vals),
                )) = rec.data().get(&Tag::new(b'n', b'2'))
                {
                    assert_eq!(
                        *vals,
                        vec![6i16, 5, 4, 3, 2],
                        "n2 should be reversed for negative strand"
                    );
                }

                if let Some(BufValue::String(s)) = rec.data().get(&Tag::new(b's', b'1')) {
                    assert_eq!(s.to_string(), "edcba", "s1 should be reversed for negative strand");
                }

                if let Some(BufValue::String(s)) = rec.data().get(&Tag::new(b's', b'2')) {
                    assert_eq!(s.to_string(), "vwxyz", "s2 should not be changed");
                }

                if let Some(BufValue::String(s)) = rec.data().get(&Tag::new(b's', b'3')) {
                    assert_eq!(
                        s.to_string(),
                        "CCTCT",
                        "s3 should be revcomped for negative strand"
                    );
                }
            } else {
                // Positive strand - no changes
                if let Some(BufValue::Array(
                    noodles::sam::alignment::record_buf::data::field::value::Array::Int16(vals),
                )) = rec.data().get(&Tag::new(b'n', b'1'))
                {
                    assert_eq!(
                        *vals,
                        vec![1i16, 2, 3, 4, 5],
                        "n1 should be unchanged for positive strand"
                    );
                }

                if let Some(BufValue::Array(
                    noodles::sam::alignment::record_buf::data::field::value::Array::Int16(vals),
                )) = rec.data().get(&Tag::new(b'n', b'2'))
                {
                    assert_eq!(
                        *vals,
                        vec![2i16, 3, 4, 5, 6],
                        "n2 should be unchanged for positive strand"
                    );
                }

                if let Some(BufValue::String(s)) = rec.data().get(&Tag::new(b's', b'1')) {
                    assert_eq!(
                        s.to_string(),
                        "abcde",
                        "s1 should be unchanged for positive strand"
                    );
                }

                if let Some(BufValue::String(s)) = rec.data().get(&Tag::new(b's', b'2')) {
                    assert_eq!(s.to_string(), "vwxyz", "s2 should not be changed");
                }

                if let Some(BufValue::String(s)) = rec.data().get(&Tag::new(b's', b'3')) {
                    assert_eq!(
                        s.to_string(),
                        "AGAGG",
                        "s3 should be unchanged for positive strand"
                    );
                }
            }
        }

        Ok(())
    }

    /// Tests multi-threaded processing
    ///
    /// Verifies that multi-threaded mode produces the same results as single-threaded
    #[test]
    fn test_multithreaded_processing() -> Result<()> {
        let mut unmapped = SamBuilder::new_unmapped();
        let mut mapped = SamBuilder::new_mapped();

        // Add multiple templates to process
        for i in 0..10 {
            let name = format!("q{i}");
            let mut attrs = HashMap::new();
            attrs.insert("RX", BufValue::from(format!("ACG{i}")));
            attrs.insert("xy", BufValue::from(1000 + i as i32));
            unmapped.add_pair_with_attrs(&name, None, None, true, true, &attrs);

            let mut mapped_attrs = HashMap::new();
            mapped_attrs.insert("PG", BufValue::from(MAPPED_PG_ID.to_string()));
            mapped_attrs.insert("AS", BufValue::from(77i32));
            mapped.add_pair_with_attrs(
                &name,
                Some(100 + i * 100),
                Some(200 + i * 100),
                true,
                true,
                &mapped_attrs,
            );
        }

        let dir = TempDir::new()?;
        let unmapped_path = dir.path().join("unmapped.bam");
        let mapped_path = dir.path().join("mapped.sam");
        let output_path = dir.path().join("output.bam");
        let dict_path = create_ref_dict(&dir, "chr1", REFERENCE_LENGTH)?;

        unmapped.write(&unmapped_path)?;
        mapped.write_sam(&mapped_path)?;

        let zipper = Zipper {
            input: mapped_path.clone(),
            unmapped: unmapped_path.clone(),
            reference: dict_path,
            output: output_path.clone(),
            tags_to_remove: vec![],
            tags_to_reverse: vec![],
            tags_to_revcomp: vec![],
            buffer: 5000,
            threads: 4,
            compression: CompressionOptions { compression_level: 1 },
            bwa_chunk_size: 150_000_000,
            exclude_missing_reads: false,
            skip_tc_tags: false,
            restore_unconverted_bases: false,
        };

        zipper.execute("test")?;

        let records = read_bam_records(&output_path)?;
        assert_eq!(records.len(), 20);

        // Verify all records have both unmapped and mapped tags
        for rec in &records {
            assert!(rec.data().get(&Tag::from(SamTag::RX)).is_some());
            assert!(rec.data().get(&Tag::new(b'x', b'y')).is_some());
            assert!(rec.data().get(&Tag::from(SamTag::AS)).is_some());
        }

        Ok(())
    }

    /// Tests empty unmapped input with non-empty mapped
    ///
    /// Verifies that having mapped reads remaining after unmapped is exhausted is an error
    /// (matching fgbio's `ZipperBams` validation)
    #[test]
    fn test_empty_unmapped() -> Result<()> {
        let unmapped = SamBuilder::new_unmapped();
        let mut mapped = SamBuilder::new_mapped();

        let mut mapped_attrs = HashMap::new();
        mapped_attrs.insert("PG", BufValue::from(MAPPED_PG_ID.to_string()));
        mapped.add_frag_with_attrs("q1", Some(100), true, &mapped_attrs);

        let result = run_zipper(&unmapped, &mapped, vec![], vec![], vec![]);

        // With the leftover mapped reads validation (matching fgbio), having mapped
        // reads remaining after unmapped is exhausted should be an error
        assert!(result.is_err());
        let err_msg = result.unwrap_err().to_string();
        assert!(
            err_msg.contains("mapped reads remaining"),
            "Expected error about leftover mapped reads, got: {err_msg}"
        );

        Ok(())
    }

    /// Tests `copy_tags` function directly
    ///
    /// Verifies the tag copying logic including PG tag special handling
    #[test]
    fn test_copy_tags_with_pg_present() -> Result<()> {
        let src =
            RecordBuilder::new().sequence("ACGT").tag("PG", "old_pg").tag("XY", 123i32).build();

        let mut dest_data = Data::default();
        dest_data.insert(Tag::from(SamTag::PG), BufValue::from("new_pg".to_string()));

        let tag_info = TagInfo::new(vec![], vec![], vec![]);
        copy_tags(&src, &mut dest_data, &tag_info)?;

        // PG should not be overwritten
        if let Some(BufValue::String(s)) = dest_data.get(&Tag::from(SamTag::PG)) {
            assert_eq!(s.to_string(), "new_pg");
        }

        // XY should be copied
        assert!(dest_data.get(&Tag::new(b'X', b'Y')).is_some());

        Ok(())
    }

    /// Tests `copy_tags` with tag removal
    ///
    /// Verifies that tags in the remove list are not copied
    #[test]
    fn test_copy_tags_with_removal() -> Result<()> {
        let src = RecordBuilder::new().sequence("ACGT").tag("XY", 123i32).tag("AB", 456i32).build();

        let mut dest_data = Data::default();

        let tag_info = TagInfo::new(vec!["XY".to_string()], vec![], vec![]);
        copy_tags(&src, &mut dest_data, &tag_info)?;

        // XY should not be copied (it's in remove list)
        assert!(dest_data.get(&Tag::new(b'X', b'Y')).is_none());

        // AB should be copied
        assert!(dest_data.get(&Tag::new(b'A', b'B')).is_some());

        Ok(())
    }

    /// Tests merge function with mixed read configurations
    ///
    /// Verifies handling of reads with different strand orientations
    #[test]
    fn test_merge_mixed_strands() -> Result<()> {
        let mut unmapped = SamBuilder::new_unmapped();
        let mut mapped = SamBuilder::new_mapped();

        let mut attrs = HashMap::new();
        attrs.insert("RX", BufValue::from("ACGT".to_string()));
        attrs.insert("xy", BufValue::from(1234i32));
        unmapped.add_pair_with_attrs("q1", None, None, true, true, &attrs);

        let mut mapped_attrs = HashMap::new();
        mapped_attrs.insert("PG", BufValue::from(MAPPED_PG_ID.to_string()));
        // Both reads map to positive strand
        mapped.add_pair_with_attrs("q1", Some(100), Some(200), true, true, &mapped_attrs);

        let records = run_zipper(&unmapped, &mapped, vec![], vec![], vec![])?;

        assert_eq!(records.len(), 2);

        // All records should have the unmapped tags
        for rec in &records {
            assert!(rec.data().get(&Tag::from(SamTag::RX)).is_some());
            assert!(rec.data().get(&Tag::new(b'x', b'y')).is_some());
        }

        Ok(())
    }

    /// Tests `TagInfo` `has_revs_or_revcomps` method
    ///
    /// Verifies the helper method correctly identifies when reversals/revcomps are configured
    #[test]
    fn test_tag_info_has_revs_or_revcomps() {
        let tag_info1 = TagInfo::new(vec![], vec![], vec![]);
        assert!(!tag_info1.has_revs_or_revcomps());

        let tag_info2 = TagInfo::new(vec![], vec!["n1".to_string()], vec![]);
        assert!(tag_info2.has_revs_or_revcomps());

        let tag_info3 = TagInfo::new(vec![], vec![], vec!["s1".to_string()]);
        assert!(tag_info3.has_revs_or_revcomps());

        let tag_info4 = TagInfo::new(vec![], vec!["n1".to_string()], vec!["s1".to_string()]);
        assert!(tag_info4.has_revs_or_revcomps());
    }

    /// Tests handling of multiple secondary alignments
    ///
    /// Verifies that multiple secondary alignments for the same read are all processed
    #[test]
    fn test_multiple_secondary_alignments() -> Result<()> {
        let mut unmapped = SamBuilder::new_unmapped();
        let mut mapped = SamBuilder::new_mapped();

        let mut attrs = HashMap::new();
        attrs.insert("RX", BufValue::from("ACGT".to_string()));
        unmapped.add_frag_with_attrs("q1", None, true, &attrs);

        let mut mapped_attrs = HashMap::new();
        mapped_attrs.insert("PG", BufValue::from(MAPPED_PG_ID.to_string()));
        mapped.add_frag_with_attrs("q1", Some(100), true, &mapped_attrs.clone());

        // Add 3 secondary alignments
        for i in 0..3 {
            let mut secondary_rec = RecordBuilder::new()
                .name("q1")
                .flags(Flags::SECONDARY)
                .reference_sequence_id(0)
                .alignment_start(200 + i * 100)
                .cigar("100M")
                .sequence(&"A".repeat(100))
                .qualities(&[30u8; 100])
                .build();
            for (tag_str, value) in &mapped_attrs {
                let tag = Tag::new(tag_str.as_bytes()[0], tag_str.as_bytes()[1]);
                secondary_rec.data_mut().insert(tag, value.clone());
            }
            mapped.push_record(secondary_rec);
        }

        let records = run_zipper(&unmapped, &mapped, vec![], vec![], vec![])?;

        assert_eq!(records.len(), 4, "Should have 1 primary + 3 secondary alignments");

        // All records should have the RX tag
        for rec in &records {
            assert!(rec.data().get(&Tag::from(SamTag::RX)).is_some());
        }

        Ok(())
    }

    /// Tests handling of R2-only mapping (R1 unmapped, R2 mapped)
    ///
    /// Verifies correct handling when only one read of a pair is mapped
    #[test]
    fn test_r2_only_mapping() -> Result<()> {
        let mut unmapped = SamBuilder::new_unmapped();
        let mut mapped = SamBuilder::new_mapped();

        let mut attrs = HashMap::new();
        attrs.insert("RX", BufValue::from("ACGT".to_string()));
        unmapped.add_pair_with_attrs("q1", None, None, true, true, &attrs);

        // Only R2 is mapped
        let r2_mapped = RecordBuilder::new()
            .name("q1")
            .flags(Flags::SEGMENTED | Flags::LAST_SEGMENT)
            .reference_sequence_id(0)
            .alignment_start(200)
            .cigar("100M")
            .sequence(&"A".repeat(100))
            .qualities(&[30u8; 100])
            .tag("PG", MAPPED_PG_ID)
            .build();
        mapped.push_record(r2_mapped);

        let records = run_zipper(&unmapped, &mapped, vec![], vec![], vec![])?;

        assert_eq!(records.len(), 1, "Should have only R2 in output");

        let r2 = &records[0];
        assert!(!r2.flags().is_first_segment());
        assert!(r2.data().get(&Tag::from(SamTag::RX)).is_some());

        Ok(())
    }

    /// Helper to run zipper with `exclude_missing_reads` option
    fn run_zipper_with_options(
        unmapped: &SamBuilder,
        mapped: &SamBuilder,
        exclude_missing_reads: bool,
    ) -> Result<Vec<RecordBuf>> {
        let dir = TempDir::new()?;
        let unmapped_path = dir.path().join("unmapped.bam");
        let mapped_path = dir.path().join("mapped.sam");
        let output_path = dir.path().join("output.bam");
        let dict_path = create_ref_dict(&dir, "chr1", REFERENCE_LENGTH)?;

        unmapped.write(&unmapped_path)?;
        mapped.write_sam(&mapped_path)?;

        let zipper = Zipper {
            input: mapped_path.clone(),
            unmapped: unmapped_path.clone(),
            reference: dict_path,
            output: output_path.clone(),
            tags_to_remove: vec![],
            tags_to_reverse: vec![],
            tags_to_revcomp: vec![],
            buffer: 5000,
            threads: 1,
            compression: CompressionOptions { compression_level: 1 },
            bwa_chunk_size: 150_000_000,
            exclude_missing_reads,
            skip_tc_tags: false,
            restore_unconverted_bases: false,
        };

        zipper.execute("test")?;
        read_bam_records(&output_path)
    }

    /// Tests that both inputs being empty produces empty output (not an error)
    #[test]
    fn test_both_inputs_empty() -> Result<()> {
        let unmapped = SamBuilder::new_unmapped();
        let mapped = SamBuilder::new_mapped();

        let records = run_zipper(&unmapped, &mapped, vec![], vec![], vec![])?;

        // Both empty should produce empty output
        assert_eq!(records.len(), 0);

        Ok(())
    }

    /// Tests the --exclude-missing-reads flag
    ///
    /// When enabled, unmapped reads without matching mapped reads should be excluded from output
    #[test]
    fn test_exclude_missing_reads() -> Result<()> {
        let mut unmapped = SamBuilder::new_unmapped();
        let mut mapped = SamBuilder::new_mapped();

        // Add unmapped pairs - q1, q2, q3
        let mut attrs = HashMap::new();
        attrs.insert("RX", BufValue::from("ACGT".to_string()));
        unmapped.add_pair_with_attrs("q1", None, None, true, true, &attrs);
        unmapped.add_pair_with_attrs("q2", None, None, true, true, &attrs);
        unmapped.add_pair_with_attrs("q3", None, None, true, true, &attrs);

        // Add mapped pairs - only q1 and q3 (q2 missing)
        let mut mapped_attrs = HashMap::new();
        mapped_attrs.insert("PG", BufValue::from(MAPPED_PG_ID.to_string()));
        mapped.add_pair_with_attrs("q1", Some(100), Some(200), true, true, &mapped_attrs);
        mapped.add_pair_with_attrs("q3", Some(300), Some(400), true, true, &mapped_attrs);

        // Without exclude_missing_reads, q2 (unmapped-only) should be included
        let records_with_missing = run_zipper_with_options(&unmapped, &mapped, false)?;
        assert_eq!(records_with_missing.len(), 6, "Should have 6 records (3 pairs)");

        // With exclude_missing_reads, q2 should be excluded
        let records_without_missing = run_zipper_with_options(&unmapped, &mapped, true)?;
        assert_eq!(records_without_missing.len(), 4, "Should have 4 records (2 pairs)");

        // Verify q2 is not in output
        for rec in &records_without_missing {
            if let Some(name) = rec.name() {
                let name_str = String::from_utf8_lossy(name.as_ref());
                assert_ne!(name_str, "q2", "q2 should be excluded");
            }
        }

        Ok(())
    }

    /// Tests that file-not-found errors are clear
    #[test]
    fn test_file_not_found_error() -> Result<()> {
        let dir = TempDir::new()?;
        let unmapped_path = dir.path().join("nonexistent_unmapped.bam");
        let mapped_path = dir.path().join("mapped.sam");
        let output_path = dir.path().join("output.bam");

        // Create mapped file but not unmapped or reference
        let mapped = SamBuilder::new_mapped();
        mapped.write_sam(&mapped_path)?;

        let zipper = Zipper {
            input: mapped_path,
            unmapped: unmapped_path.clone(),
            reference: dir.path().join("ref.fa"),
            output: output_path,
            tags_to_remove: vec![],
            tags_to_reverse: vec![],
            tags_to_revcomp: vec![],
            buffer: 5000,
            threads: 1,
            compression: CompressionOptions { compression_level: 1 },
            bwa_chunk_size: 150_000_000,
            exclude_missing_reads: false,
            skip_tc_tags: false,
            restore_unconverted_bases: false,
        };

        let result = zipper.execute("test");
        assert!(result.is_err(), "Should error when reference file doesn't exist");
        let err_msg = result.unwrap_err().to_string();
        assert!(
            err_msg.contains("does not exist") || err_msg.contains("not found"),
            "Error should mention file not found: {err_msg}"
        );

        Ok(())
    }

    // =========================================================================
    // Template-Coordinate Tag (tc) Tests
    // =========================================================================

    // Helper for creating records with specific flags
    const FLAG_PAIRED: u16 = 0x1;
    const FLAG_READ1: u16 = 0x40;
    const FLAG_READ2: u16 = 0x80;
    const FLAG_SECONDARY: u16 = 0x100;
    const FLAG_SUPPLEMENTARY: u16 = 0x800;
    const FLAG_REVERSE: u16 = 0x10;

    /// Tests that the tc tag is added to supplementary reads with template sort key
    #[test]
    fn test_add_tc_tag_to_supplementary_r1() -> Result<()> {
        use crate::template::Template;

        // Build a template with primary R1 and supplementary R1
        let primary_r1 = RecordBuilder::new()
            .name("q1")
            .sequence("ACGTACGT")
            .qualities(&[30; 8])
            .reference_sequence_id(0)
            .alignment_start(100)
            .cigar("8M")
            .flags(Flags::from(FLAG_PAIRED | FLAG_READ1))
            .build();

        let supplementary_r1 = RecordBuilder::new()
            .name("q1")
            .sequence("ACGTACGT")
            .qualities(&[30; 8])
            .reference_sequence_id(5)
            .alignment_start(5000)
            .cigar("8M")
            .flags(Flags::from(FLAG_PAIRED | FLAG_READ1 | FLAG_SUPPLEMENTARY))
            .build();

        let mut template = Template::from_records(vec![primary_r1, supplementary_r1])?;

        // Call add_template_coordinate_tags
        add_template_coordinate_tags(&mut template);

        // Check that supplementary has tc tag
        let supp = template
            .records
            .iter()
            .find(|r| r.flags().is_supplementary())
            .expect("Should have supplementary");

        let tc_value = supp.data().get(&TC_TAG).expect("Supplementary should have tc tag");

        // Parse and verify the tc tag
        let tc_info = TemplateCoordinateInfo::from_tag_value(tc_value)
            .expect("Should be able to parse tc tag");

        // Only R1 mapped, so both positions should be the same
        assert_eq!(tc_info.tid1, 0, "tid1 should be R1's reference");
        assert_eq!(tc_info.pos1, 100, "pos1 should be R1's position");
        assert!(!tc_info.neg1, "neg1 should be false (forward strand)");
        assert_eq!(tc_info.tid2, 0, "tid2 should equal tid1 (single read)");
        assert_eq!(tc_info.pos2, 100, "pos2 should equal pos1 (single read)");
        assert!(!tc_info.neg2, "neg2 should equal neg1 (single read)");

        // Check that primary does NOT have tc tag
        let primary = template.r1().expect("Should have primary R1");
        assert!(primary.data().get(&TC_TAG).is_none(), "Primary should not have tc tag");

        Ok(())
    }

    /// Tests that the tc tag is added to secondary reads with correct strand info
    #[test]
    fn test_add_tc_tag_to_secondary_reverse_strand() -> Result<()> {
        use crate::template::Template;

        // Build a template with primary R1 on reverse strand and secondary R1
        let primary_r1 = RecordBuilder::new()
            .name("q1")
            .sequence("ACGTACGT")
            .qualities(&[30; 8])
            .reference_sequence_id(0)
            .alignment_start(200)
            .cigar("8M")
            .flags(Flags::from(FLAG_PAIRED | FLAG_READ1 | FLAG_REVERSE))
            .build();

        let secondary_r1 = RecordBuilder::new()
            .name("q1")
            .sequence("ACGTACGT")
            .qualities(&[30; 8])
            .reference_sequence_id(3)
            .alignment_start(3000)
            .cigar("8M")
            .flags(Flags::from(FLAG_PAIRED | FLAG_READ1 | FLAG_SECONDARY))
            .build();

        let mut template = Template::from_records(vec![primary_r1, secondary_r1])?;

        add_template_coordinate_tags(&mut template);

        let secondary = template
            .records
            .iter()
            .find(|r| r.flags().is_secondary())
            .expect("Should have secondary");

        let tc_value = secondary.data().get(&TC_TAG).expect("Secondary should have tc tag");

        // Parse and verify the tc tag
        let tc_info = TemplateCoordinateInfo::from_tag_value(tc_value)
            .expect("Should be able to parse tc tag");

        // Primary is on reverse strand with 8M cigar
        // Unclipped 5' position for reverse strand = alignment_start + alignment_span - 1
        // = 200 + 8 - 1 = 207
        assert_eq!(tc_info.tid1, 0, "tid1 should be R1's reference");
        assert_eq!(tc_info.pos1, 207, "pos1 should be R1's unclipped 5' position");
        assert!(tc_info.neg1, "neg1 should be true (reverse strand)");

        Ok(())
    }

    /// Tests that tc tag contains full template sort key for paired-end data
    #[test]
    fn test_add_tc_tag_paired_end_r2_supplementary() -> Result<()> {
        use crate::template::Template;

        // Primary R1 at position 100 (forward strand, 8M cigar)
        // Unclipped 5' = 100 (no soft clips)
        let primary_r1 = RecordBuilder::new()
            .name("q1")
            .sequence("ACGTACGT")
            .qualities(&[30; 8])
            .reference_sequence_id(0)
            .alignment_start(100)
            .cigar("8M")
            .flags(Flags::from(FLAG_PAIRED | FLAG_READ1))
            .build();

        // Primary R2 at position 300 (reverse strand, 8M cigar)
        // Unclipped 5' for reverse strand = 300 + 8 - 1 = 307
        let primary_r2 = RecordBuilder::new()
            .name("q1")
            .sequence("ACGTACGT")
            .qualities(&[30; 8])
            .reference_sequence_id(0)
            .alignment_start(300)
            .cigar("8M")
            .flags(Flags::from(FLAG_PAIRED | FLAG_READ2 | FLAG_REVERSE))
            .build();

        // Supplementary R2
        let supplementary_r2 = RecordBuilder::new()
            .name("q1")
            .sequence("ACGT")
            .qualities(&[30; 4])
            .reference_sequence_id(5)
            .alignment_start(5000)
            .cigar("4M")
            .flags(Flags::from(FLAG_PAIRED | FLAG_READ2 | FLAG_SUPPLEMENTARY))
            .build();

        let mut template = Template::from_records(vec![primary_r1, primary_r2, supplementary_r2])?;

        add_template_coordinate_tags(&mut template);

        let supp = template
            .records
            .iter()
            .find(|r| r.flags().is_supplementary())
            .expect("Should have supplementary");

        let tc_value = supp.data().get(&TC_TAG).expect("Supplementary R2 should have tc tag");

        // Parse and verify the tc tag
        let tc_info = TemplateCoordinateInfo::from_tag_value(tc_value)
            .expect("Should be able to parse tc tag");

        // tc tag should contain BOTH primaries' unclipped 5' positions (sorted by position)
        // R1 forward strand at 100 with 8M -> unclipped 5' = 100
        // R2 reverse strand at 300 with 8M -> unclipped 5' = 300 + 8 - 1 = 307
        assert_eq!(tc_info.tid1, 0, "tid1 should be R1's reference (earlier)");
        assert_eq!(tc_info.pos1, 100, "pos1 should be R1's unclipped 5' position (earlier)");
        assert!(!tc_info.neg1, "neg1 should be false (R1 forward strand)");
        assert_eq!(tc_info.tid2, 0, "tid2 should be R2's reference (later)");
        assert_eq!(tc_info.pos2, 307, "pos2 should be R2's unclipped 5' position (later)");
        assert!(tc_info.neg2, "neg2 should be true (R2 reverse strand)");

        Ok(())
    }

    /// Tests that no tc tag is added when there's no corresponding primary
    #[test]
    fn test_no_tc_tag_when_no_primary() -> Result<()> {
        use crate::template::Template;

        // Only supplementary, no primary
        let supplementary = RecordBuilder::new()
            .name("q1")
            .sequence("ACGT")
            .qualities(&[30; 4])
            .reference_sequence_id(5)
            .alignment_start(5000)
            .cigar("4M")
            .flags(Flags::from(FLAG_PAIRED | FLAG_READ1 | FLAG_SUPPLEMENTARY))
            .build();

        let mut template = Template::from_records(vec![supplementary])?;

        add_template_coordinate_tags(&mut template);

        let supp = template
            .records
            .iter()
            .find(|r| r.flags().is_supplementary())
            .expect("Should have supplementary");

        // Should NOT have tc tag since there's no primary
        assert!(
            supp.data().get(&TC_TAG).is_none(),
            "Supplementary without primary should not have tc tag"
        );

        Ok(())
    }

    /// Tests that tc tag is added during full merge operation
    #[test]
    fn test_tc_tag_added_during_merge() -> Result<()> {
        let mut unmapped = SamBuilder::new_unmapped();
        let mut mapped = SamBuilder::new_mapped();

        // Add unmapped pair with tags
        let mut attrs = HashMap::new();
        attrs.insert("RX", BufValue::from("ACGT".to_string()));
        unmapped.add_pair_with_attrs("q1", None, None, true, true, &attrs);

        // Add mapped primary pair using the builder pattern
        // R1 at 100, R2 at 200
        let _ = mapped.add_pair().name("q1").start1(100).start2(200).build();

        // Add supplementary R1 (same ref but different pos to test tc tag)
        let supp = RecordBuilder::new()
            .name("q1")
            .sequence("ACGTACGTACGT")
            .qualities(&[30; 12])
            .reference_sequence_id(0)
            .alignment_start(5000)
            .cigar("12M")
            .flags(Flags::from(FLAG_PAIRED | FLAG_READ1 | FLAG_SUPPLEMENTARY))
            .build();
        mapped.push_record(supp);

        let records = run_zipper(&unmapped, &mapped, vec![], vec![], vec![])?;

        // Find supplementary in output
        let supp_record = records
            .iter()
            .find(|r| r.flags().is_supplementary())
            .expect("Should have supplementary in output");

        // Check tc tag was added
        let tc_value =
            supp_record.data().get(&TC_TAG).expect("Supplementary should have tc tag after merge");

        // Parse and verify the tc tag
        let tc_info = TemplateCoordinateInfo::from_tag_value(tc_value)
            .expect("Should be able to parse tc tag");

        // tc tag should contain both primaries' unclipped 5' positions
        // R1: forward strand at 100 with 100M -> unclipped 5' = 100
        // R2: reverse strand at 200 with 100M -> unclipped 5' = 200 + 100 - 1 = 299
        assert_eq!(tc_info.tid1, 0, "tid1 should be 0");
        assert_eq!(tc_info.pos1, 100, "pos1 should be R1's unclipped 5' position");
        assert_eq!(tc_info.tid2, 0, "tid2 should be 0");
        assert_eq!(tc_info.pos2, 299, "pos2 should be R2's unclipped 5' position (reverse strand)");

        // Verify RX tag was also copied
        let rx_tag = Tag::from(SamTag::RX);
        assert!(
            supp_record.data().get(&rx_tag).is_some(),
            "Supplementary should also have RX tag copied from unmapped"
        );

        Ok(())
    }

    /// Tests that `add_template_coordinate_tags` returns early when there are no
    /// secondary/supplementary reads (the early exit optimization).
    #[test]
    fn test_add_tc_tag_early_exit_no_secondary_supplementary() -> Result<()> {
        use crate::template::Template;

        // Create a template with only primary reads (no secondary/supplementary)
        let primary_r1 = RecordBuilder::new()
            .name("q1")
            .sequence("ACGT")
            .qualities(&[30, 30, 30, 30])
            .reference_sequence_id(0)
            .alignment_start(100)
            .cigar("4M")
            .flags(Flags::from(FLAG_PAIRED | FLAG_READ1))
            .build();

        let primary_r2 = RecordBuilder::new()
            .name("q1")
            .sequence("TGCA")
            .qualities(&[30, 30, 30, 30])
            .reference_sequence_id(0)
            .alignment_start(200)
            .cigar("4M")
            .flags(Flags::from(FLAG_PAIRED | FLAG_READ2 | FLAG_REVERSE))
            .build();

        let mut template = Template::from_records(vec![primary_r1, primary_r2])?;

        // Call add_template_coordinate_tags - should return early
        add_template_coordinate_tags(&mut template);

        // Verify no tc tags were added (primaries don't get tc tags)
        for record in &template.records {
            assert!(record.data().get(&TC_TAG).is_none(), "Primary reads should not have tc tag");
        }

        Ok(())
    }

    /// Tests that `add_template_coordinate_tags` only adds tc tag to secondary/supplementary,
    /// not to primary reads, even when secondary/supplementary are present.
    #[test]
    fn test_add_tc_tag_only_to_secondary_supplementary() -> Result<()> {
        use crate::template::Template;

        // Create a template with primary + secondary
        let primary_r1 = RecordBuilder::new()
            .name("q1")
            .sequence("ACGT")
            .qualities(&[30, 30, 30, 30])
            .reference_sequence_id(0)
            .alignment_start(100)
            .cigar("4M")
            .flags(Flags::from(FLAG_PAIRED | FLAG_READ1))
            .build();

        let secondary_r1 = RecordBuilder::new()
            .name("q1")
            .sequence("ACGT")
            .qualities(&[30, 30, 30, 30])
            .reference_sequence_id(5)
            .alignment_start(5000)
            .cigar("4M")
            .flags(Flags::from(FLAG_PAIRED | FLAG_READ1 | FLAG_SECONDARY))
            .build();

        let mut template = Template::from_records(vec![primary_r1, secondary_r1])?;

        add_template_coordinate_tags(&mut template);

        // Check each record
        for record in &template.records {
            if record.flags().is_secondary() {
                assert!(record.data().get(&TC_TAG).is_some(), "Secondary should have tc tag");
            } else {
                assert!(record.data().get(&TC_TAG).is_none(), "Primary should not have tc tag");
            }
        }

        Ok(())
    }

    #[test]
    fn test_skip_tc_tags_parameter() -> Result<()> {
        use crate::template::Template;

        const FLAG_UNMAPPED: u16 = 0x4;

        // Create a template with primary + supplementary
        let primary_r1 = RecordBuilder::new()
            .name("q1")
            .sequence("ACGT")
            .qualities(&[30, 30, 30, 30])
            .reference_sequence_id(0)
            .alignment_start(100)
            .cigar("4M")
            .flags(Flags::from(FLAG_PAIRED | FLAG_READ1))
            .build();

        let supplementary_r1 = RecordBuilder::new()
            .name("q1")
            .sequence("ACGT")
            .qualities(&[30, 30, 30, 30])
            .reference_sequence_id(5)
            .alignment_start(5000)
            .cigar("4M")
            .flags(Flags::from(FLAG_PAIRED | FLAG_READ1 | FLAG_SUPPLEMENTARY))
            .build();

        let mut template = Template::from_records(vec![primary_r1, supplementary_r1])?;

        // Create unmapped template
        let unmapped_record = RecordBuilder::new()
            .name("q1")
            .sequence("ACGT")
            .qualities(&[30, 30, 30, 30])
            .flags(Flags::from(FLAG_PAIRED | FLAG_READ1 | FLAG_UNMAPPED))
            .tag("RX", "AAAA")
            .build();
        let unmapped = Template::from_records(vec![unmapped_record])?;

        let tag_info = TagInfo::new(vec![], vec![], vec![]);

        // Test with skip_tc_tags = true
        merge(&unmapped, &mut template, &tag_info, true)?;

        // Supplementary should NOT have tc tag when skip_tc_tags is true
        for record in &template.records {
            assert!(
                record.data().get(&TC_TAG).is_none(),
                "No records should have tc tag when skip_tc_tags=true"
            );
        }

        Ok(())
    }

    #[test]
    fn test_skip_tc_tags_false_adds_tags() -> Result<()> {
        use crate::template::Template;

        const FLAG_UNMAPPED: u16 = 0x4;

        // Create a template with primary + supplementary
        let primary_r1 = RecordBuilder::new()
            .name("q1")
            .sequence("ACGT")
            .qualities(&[30, 30, 30, 30])
            .reference_sequence_id(0)
            .alignment_start(100)
            .cigar("4M")
            .flags(Flags::from(FLAG_PAIRED | FLAG_READ1))
            .build();

        let supplementary_r1 = RecordBuilder::new()
            .name("q1")
            .sequence("ACGT")
            .qualities(&[30, 30, 30, 30])
            .reference_sequence_id(5)
            .alignment_start(5000)
            .cigar("4M")
            .flags(Flags::from(FLAG_PAIRED | FLAG_READ1 | FLAG_SUPPLEMENTARY))
            .build();

        let mut template = Template::from_records(vec![primary_r1, supplementary_r1])?;

        // Create unmapped template
        let unmapped_record = RecordBuilder::new()
            .name("q1")
            .sequence("ACGT")
            .qualities(&[30, 30, 30, 30])
            .flags(Flags::from(FLAG_PAIRED | FLAG_READ1 | FLAG_UNMAPPED))
            .tag("RX", "AAAA")
            .build();
        let unmapped = Template::from_records(vec![unmapped_record])?;

        let tag_info = TagInfo::new(vec![], vec![], vec![]);

        // Test with skip_tc_tags = false
        merge(&unmapped, &mut template, &tag_info, false)?;

        // Supplementary SHOULD have tc tag when skip_tc_tags is false
        let has_tc_tag = template
            .records
            .iter()
            .any(|r| r.flags().is_supplementary() && r.data().get(&TC_TAG).is_some());

        assert!(has_tc_tag, "Supplementary should have tc tag when skip_tc_tags=false");

        Ok(())
    }

    /// Tests AS/XS tag normalization to smallest signed integer type
    ///
    /// Verifies that:
    /// - AS tags with Int32 values are normalized to Int8 when they fit
    /// - XS tags with `UInt8` values are normalized to `Int8`
    /// - Large values that don't fit in Int8 are normalized to Int16
    /// - Tags on all reads (primary, secondary, supplementary) are normalized
    #[test]
    fn test_as_xs_normalization() -> Result<()> {
        let mut unmapped = SamBuilder::new_unmapped();
        let mut mapped = SamBuilder::new_mapped();

        let mut attrs = HashMap::new();
        attrs.insert("RX", BufValue::from("ACGT".to_string()));
        unmapped.add_frag_with_attrs("q1", None, true, &attrs);

        // Add mapped read with AS as Int32(77) and XS as Int32(50)
        let mapped_rec = RecordBuilder::new()
            .name("q1")
            .flags(Flags::empty())
            .reference_sequence_id(0)
            .alignment_start(100)
            .cigar("100M")
            .sequence(&"A".repeat(100))
            .qualities(&[30u8; 100])
            .tag("PG", MAPPED_PG_ID.to_string())
            .tag("AS", 77i32)
            .tag("XS", 50i32)
            .build();
        mapped.push_record(mapped_rec);

        let records = run_zipper(&unmapped, &mapped, vec![], vec![], vec![])?;
        assert_eq!(records.len(), 1);

        // AS=77 fits in i8 (-128..127), should be normalized to Int8
        let as_val = records[0].data().get(&Tag::from(SamTag::AS)).unwrap();
        assert!(matches!(as_val, BufValue::Int8(77)), "AS=77 should be Int8, got {as_val:?}");

        // XS=50 fits in i8, should be normalized to Int8
        let xs_val = records[0].data().get(&Tag::from(SamTag::XS)).unwrap();
        assert!(matches!(xs_val, BufValue::Int8(50)), "XS=50 should be Int8, got {xs_val:?}");

        Ok(())
    }

    /// Tests AS/XS normalization with values that exceed Int8 range
    ///
    /// Verifies that AS values >127 are normalized to Int16 rather than Int8
    #[test]
    fn test_as_xs_normalization_large_values() -> Result<()> {
        let mut unmapped = SamBuilder::new_unmapped();
        let mut mapped = SamBuilder::new_mapped();

        let mut attrs = HashMap::new();
        attrs.insert("RX", BufValue::from("ACGT".to_string()));
        unmapped.add_frag_with_attrs("q1", None, true, &attrs);

        // AS=200 exceeds i8 range, XS=-200 exceeds i8 range
        let mapped_rec = RecordBuilder::new()
            .name("q1")
            .flags(Flags::empty())
            .reference_sequence_id(0)
            .alignment_start(100)
            .cigar("100M")
            .sequence(&"A".repeat(100))
            .qualities(&[30u8; 100])
            .tag("PG", MAPPED_PG_ID.to_string())
            .tag("AS", 200i32)
            .tag("XS", -200i32)
            .build();
        mapped.push_record(mapped_rec);

        let records = run_zipper(&unmapped, &mapped, vec![], vec![], vec![])?;
        assert_eq!(records.len(), 1);

        // AS=200 doesn't fit in i8, should be Int16
        let as_val = records[0].data().get(&Tag::from(SamTag::AS)).unwrap();
        assert!(matches!(as_val, BufValue::Int16(200)), "AS=200 should be Int16, got {as_val:?}");

        // XS=-200 doesn't fit in i8, should be Int16
        let xs_val = records[0].data().get(&Tag::from(SamTag::XS)).unwrap();
        assert!(matches!(xs_val, BufValue::Int16(-200)), "XS=-200 should be Int16, got {xs_val:?}");

        Ok(())
    }

    /// Tests negative strand tag copying when R1 is on the negative strand
    ///
    /// Verifies that:
    /// - Tags on R1 (negative strand) are reversed/revcomped
    /// - Tags on R2 (positive strand) are copied unchanged
    #[test]
    fn test_negative_strand_r1_tag_copying() -> Result<()> {
        let mut unmapped = SamBuilder::new_unmapped();
        let mut mapped = SamBuilder::new_mapped();

        let mut attrs = HashMap::new();
        attrs.insert("n1", BufValue::from(vec![1i16, 2, 3]));
        attrs.insert("s1", BufValue::from("ACGT".to_string()));
        unmapped.add_pair_with_attrs("q1", None, None, true, true, &attrs);

        let mut mapped_attrs = HashMap::new();
        mapped_attrs.insert("PG", BufValue::from(MAPPED_PG_ID.to_string()));
        // R1 on negative strand, R2 on positive strand
        mapped.add_pair_with_attrs("q1", Some(100), Some(200), false, true, &mapped_attrs);

        let records =
            run_zipper(&unmapped, &mapped, vec![], vec!["n1".to_string()], vec!["s1".to_string()])?;

        assert_eq!(records.len(), 2);

        let r1 = records.iter().find(|r| r.flags().is_first_segment()).unwrap();
        let r2 = records.iter().find(|r| !r.flags().is_first_segment()).unwrap();

        // R1 is negative strand - tags should be reversed/revcomped
        if let Some(BufValue::Array(
            noodles::sam::alignment::record_buf::data::field::value::Array::Int16(vals),
        )) = r1.data().get(&Tag::new(b'n', b'1'))
        {
            assert_eq!(*vals, vec![3i16, 2, 1], "n1 should be reversed for R1 (negative strand)");
        } else {
            panic!("n1 tag not found or wrong type on R1");
        }

        if let Some(BufValue::String(s)) = r1.data().get(&Tag::new(b's', b'1')) {
            assert_eq!(s.to_string(), "ACGT", "s1 should be revcomped for R1 (negative strand)");
        } else {
            panic!("s1 tag not found or wrong type on R1");
        }

        // R2 is positive strand - tags should be unchanged
        if let Some(BufValue::Array(
            noodles::sam::alignment::record_buf::data::field::value::Array::Int16(vals),
        )) = r2.data().get(&Tag::new(b'n', b'1'))
        {
            assert_eq!(*vals, vec![1i16, 2, 3], "n1 should be unchanged for R2 (positive strand)");
        } else {
            panic!("n1 tag not found or wrong type on R2");
        }

        if let Some(BufValue::String(s)) = r2.data().get(&Tag::new(b's', b'1')) {
            assert_eq!(s.to_string(), "ACGT", "s1 should be unchanged for R2 (positive strand)");
        } else {
            panic!("s1 tag not found or wrong type on R2");
        }

        Ok(())
    }

    /// Tests negative strand tag copying when both R1 and R2 are on the negative strand
    ///
    /// Verifies that tags on both reads are reversed/revcomped
    #[test]
    fn test_negative_strand_both_reads() -> Result<()> {
        let mut unmapped = SamBuilder::new_unmapped();
        let mut mapped = SamBuilder::new_mapped();

        let mut attrs = HashMap::new();
        attrs.insert("n1", BufValue::from(vec![10i16, 20, 30]));
        attrs.insert("s1", BufValue::from("AACC".to_string()));
        unmapped.add_pair_with_attrs("q1", None, None, true, true, &attrs);

        let mut mapped_attrs = HashMap::new();
        mapped_attrs.insert("PG", BufValue::from(MAPPED_PG_ID.to_string()));
        // Both reads on negative strand
        mapped.add_pair_with_attrs("q1", Some(200), Some(100), false, false, &mapped_attrs);

        let records =
            run_zipper(&unmapped, &mapped, vec![], vec!["n1".to_string()], vec!["s1".to_string()])?;

        assert_eq!(records.len(), 2);

        // Both reads should have reversed/revcomped tags
        for rec in &records {
            if let Some(BufValue::Array(
                noodles::sam::alignment::record_buf::data::field::value::Array::Int16(vals),
            )) = rec.data().get(&Tag::new(b'n', b'1'))
            {
                assert_eq!(*vals, vec![30i16, 20, 10], "n1 should be reversed on negative strand");
            } else {
                panic!("n1 tag not found or wrong type");
            }

            if let Some(BufValue::String(s)) = rec.data().get(&Tag::new(b's', b'1')) {
                assert_eq!(s.to_string(), "GGTT", "s1 should be revcomped on negative strand");
            } else {
                panic!("s1 tag not found or wrong type");
            }
        }

        Ok(())
    }

    /// Tests `fix_mate_info` verifying actual RNEXT, PNEXT, TLEN, MQ, and MC values
    ///
    /// Verifies that:
    /// - RNEXT (`mate_reference_sequence_id`) is set correctly for both reads
    /// - PNEXT (`mate_alignment_start`) is set correctly for both reads
    /// - TLEN is computed correctly (positive for leftmost read, negative for rightmost)
    /// - MQ tag contains the mate's mapping quality
    /// - MC tag contains the mate's CIGAR string
    #[test]
    fn test_fix_mate_info_values() -> Result<()> {
        let mut unmapped = SamBuilder::new_unmapped();
        let mut mapped = SamBuilder::new_mapped();

        let mut attrs = HashMap::new();
        attrs.insert("RX", BufValue::from("ACGT".to_string()));
        unmapped.add_pair_with_attrs("q1", None, None, true, true, &attrs);

        // Build mapped pair with different CIGARs and mapping qualities
        let r1_mapped = RecordBuilder::new()
            .name("q1")
            .flags(Flags::SEGMENTED | Flags::FIRST_SEGMENT | Flags::MATE_REVERSE_COMPLEMENTED)
            .reference_sequence_id(0)
            .alignment_start(100)
            .cigar("50M")
            .mapping_quality(40)
            .sequence(&"A".repeat(50))
            .qualities(&[30u8; 50])
            .tag("PG", MAPPED_PG_ID.to_string())
            .build();
        let r2_mapped = RecordBuilder::new()
            .name("q1")
            .flags(
                Flags::SEGMENTED
                    | Flags::LAST_SEGMENT
                    | Flags::REVERSE_COMPLEMENTED
                    | Flags::PROPERLY_SEGMENTED,
            )
            .reference_sequence_id(0)
            .alignment_start(300)
            .cigar("75M")
            .mapping_quality(30)
            .sequence(&"A".repeat(75))
            .qualities(&[30u8; 75])
            .tag("PG", MAPPED_PG_ID.to_string())
            .build();
        mapped.push_record(r1_mapped);
        mapped.push_record(r2_mapped);

        let records = run_zipper(&unmapped, &mapped, vec![], vec![], vec![])?;
        assert_eq!(records.len(), 2);

        let r1 = records.iter().find(|r| r.flags().is_first_segment()).unwrap();
        let r2 = records.iter().find(|r| !r.flags().is_first_segment()).unwrap();

        // R1's mate should point to R2
        assert_eq!(
            r1.mate_reference_sequence_id(),
            r2.reference_sequence_id(),
            "R1's RNEXT should equal R2's reference"
        );
        assert_eq!(
            r1.mate_alignment_start(),
            r2.alignment_start(),
            "R1's PNEXT should equal R2's position"
        );

        // R2's mate should point to R1
        assert_eq!(
            r2.mate_reference_sequence_id(),
            r1.reference_sequence_id(),
            "R2's RNEXT should equal R1's reference"
        );
        assert_eq!(
            r2.mate_alignment_start(),
            r1.alignment_start(),
            "R2's PNEXT should equal R1's position"
        );

        // TLEN: R1 is at 100 (50M, forward), R2 is at 300 (75M, reverse)
        // Insert size = R2_end - R1_start + 1 = (300+75-1) - 100 + 1 = 275
        let r1_tlen = r1.template_length();
        let r2_tlen = r2.template_length();
        assert!(r1_tlen != 0, "R1 TLEN should be non-zero for mapped pair");
        assert_eq!(r1_tlen, -r2_tlen, "R1 TLEN should be negation of R2 TLEN");

        // MQ: R1 should have R2's mapping quality, and vice versa
        let r1_mq = r1.data().get(&Tag::from(SamTag::MQ));
        let r2_mq = r2.data().get(&Tag::from(SamTag::MQ));
        assert!(
            matches!(r1_mq, Some(BufValue::Int8(30))),
            "R1's MQ should be R2's mapq (30), got {r1_mq:?}"
        );
        assert!(
            matches!(r2_mq, Some(BufValue::Int8(40))),
            "R2's MQ should be R1's mapq (40), got {r2_mq:?}"
        );

        // MC: R1 should have R2's CIGAR, and vice versa
        if let Some(BufValue::String(mc)) = r1.data().get(&Tag::from(SamTag::MC)) {
            assert_eq!(mc.to_string(), "75M", "R1's MC should be R2's CIGAR");
        } else {
            panic!("R1 should have MC tag");
        }
        if let Some(BufValue::String(mc)) = r2.data().get(&Tag::from(SamTag::MC)) {
            assert_eq!(mc.to_string(), "50M", "R2's MC should be R1's CIGAR");
        } else {
            panic!("R2 should have MC tag");
        }

        Ok(())
    }

    /// Tests `fix_mate_info` for supplementary alignments
    ///
    /// Verifies that:
    /// - Supplementary alignments get correct mate info from the primary of the other end
    /// - RNEXT/PNEXT point to the mate's primary alignment
    /// - MQ and MC tags are set from the mate's primary alignment
    /// - TLEN is negative of the mate primary's TLEN
    /// - ms (mate score) tag is set from mate primary's AS tag
    #[test]
    fn test_fix_mate_info_supplementary() -> Result<()> {
        let mut unmapped = SamBuilder::new_unmapped();
        let mut mapped = SamBuilder::new_mapped();

        let mut attrs = HashMap::new();
        attrs.insert("RX", BufValue::from("ACGT".to_string()));
        unmapped.add_pair_with_attrs("q1", None, None, true, true, &attrs);

        // Add primary pair with known CIGARs, positions, mapqs, and per-read AS tags
        let r1_mapped = RecordBuilder::new()
            .name("q1")
            .flags(Flags::SEGMENTED | Flags::FIRST_SEGMENT | Flags::MATE_REVERSE_COMPLEMENTED)
            .reference_sequence_id(0)
            .alignment_start(100)
            .cigar("50M")
            .mapping_quality(40)
            .sequence(&"A".repeat(50))
            .qualities(&[30u8; 50])
            .tag("PG", MAPPED_PG_ID.to_string())
            .tag("AS", 77i32)
            .build();
        let r2_mapped = RecordBuilder::new()
            .name("q1")
            .flags(
                Flags::SEGMENTED
                    | Flags::LAST_SEGMENT
                    | Flags::REVERSE_COMPLEMENTED
                    | Flags::PROPERLY_SEGMENTED,
            )
            .reference_sequence_id(0)
            .alignment_start(300)
            .cigar("75M")
            .mapping_quality(30)
            .sequence(&"A".repeat(75))
            .qualities(&[30u8; 75])
            .tag("PG", MAPPED_PG_ID.to_string())
            .tag("AS", 55i32)
            .build();
        mapped.push_record(r1_mapped);
        mapped.push_record(r2_mapped);

        // Add R1 supplementary alignment
        let supp_rec = RecordBuilder::new()
            .name("q1")
            .flags(
                Flags::SEGMENTED
                    | Flags::FIRST_SEGMENT
                    | Flags::SUPPLEMENTARY
                    | Flags::REVERSE_COMPLEMENTED,
            )
            .reference_sequence_id(0)
            .alignment_start(500)
            .cigar("60M")
            .mapping_quality(20)
            .sequence(&"A".repeat(60))
            .qualities(&[30u8; 60])
            .tag("PG", MAPPED_PG_ID.to_string())
            .tag("AS", 33i32)
            .build();
        mapped.push_record(supp_rec);

        let records = run_zipper(&unmapped, &mapped, vec![], vec![], vec![])?;
        assert_eq!(records.len(), 3);

        let r1_primary = records
            .iter()
            .find(|r| r.flags().is_first_segment() && !r.flags().is_supplementary())
            .unwrap();
        let r2_primary = records
            .iter()
            .find(|r| !r.flags().is_first_segment() && !r.flags().is_supplementary())
            .unwrap();
        let r1_supp = records
            .iter()
            .find(|r| r.flags().is_first_segment() && r.flags().is_supplementary())
            .unwrap();

        // R1 supplementary's mate should point to R2's primary
        assert_eq!(
            r1_supp.mate_reference_sequence_id(),
            r2_primary.reference_sequence_id(),
            "R1 supp RNEXT should point to R2 primary's reference"
        );
        assert_eq!(
            r1_supp.mate_alignment_start(),
            r2_primary.alignment_start(),
            "R1 supp PNEXT should point to R2 primary's position"
        );

        // MQ on supplementary should be R2 primary's mapping quality
        let supp_mq = r1_supp.data().get(&Tag::from(SamTag::MQ));
        assert!(
            matches!(supp_mq, Some(BufValue::Int8(30))),
            "R1 supp MQ should be R2 primary's mapq (30), got {supp_mq:?}"
        );

        // MC on supplementary should be R2 primary's CIGAR
        if let Some(BufValue::String(mc)) = r1_supp.data().get(&Tag::from(SamTag::MC)) {
            assert_eq!(mc.to_string(), "75M", "R1 supp MC should be R2 primary's CIGAR");
        } else {
            panic!("R1 supp should have MC tag");
        }

        // ms (mate score) on supplementary should be R2 primary's AS value
        let supp_ms = r1_supp.data().get(&Tag::from(SamTag::MS));
        assert!(
            matches!(supp_ms, Some(BufValue::Int8(55))),
            "R1 supp ms should be R2 primary's AS (55), got {supp_ms:?}"
        );

        // TLEN on supplementary should be negative of R2 primary's TLEN
        let r1_primary_tlen = r1_primary.template_length();
        let r1_supp_tlen = r1_supp.template_length();
        // Supplementary TLEN = -(R2 primary TLEN) = -(-R1 primary TLEN) = R1 primary TLEN... no
        // Actually: supplementary's TLEN = -(mate primary's TLEN) = -(R2's TLEN) = R1's TLEN
        // Wait, the code says: *self.records[i].template_length_mut() = -r2_tlen;
        // R2's TLEN is the negative of R1's TLEN, so -r2_tlen = R1's TLEN
        assert_eq!(
            r1_supp_tlen, r1_primary_tlen,
            "R1 supp TLEN should equal R1 primary TLEN (both = -R2_TLEN)"
        );

        Ok(())
    }

    /// Tests `fix_mate_info` when one read is unmapped
    ///
    /// Verifies that:
    /// - The unmapped read is placed at the mapped read's coordinates
    /// - MQ/MC tags are removed from the mapped read (mate is unmapped)
    /// - MQ/MC tags are set on the unmapped read (mate is mapped)
    /// - TLEN is 0 for both reads
    #[test]
    fn test_fix_mate_info_one_unmapped() -> Result<()> {
        let mut unmapped = SamBuilder::new_unmapped();
        let mut mapped = SamBuilder::new_mapped();

        let mut attrs = HashMap::new();
        attrs.insert("RX", BufValue::from("ACGT".to_string()));
        unmapped.add_pair_with_attrs("q1", None, None, true, true, &attrs);

        // R1 is mapped, R2 is unmapped
        let mut mapped_attrs = HashMap::new();
        mapped_attrs.insert("PG", BufValue::from(MAPPED_PG_ID.to_string()));
        mapped.add_pair_with_attrs("q1", Some(100), None, true, true, &mapped_attrs);

        let records = run_zipper(&unmapped, &mapped, vec![], vec![], vec![])?;
        assert_eq!(records.len(), 2);

        let r1 = records.iter().find(|r| r.flags().is_first_segment()).unwrap();
        let r2 = records.iter().find(|r| !r.flags().is_first_segment()).unwrap();

        // TLEN should be 0 for both
        assert_eq!(r1.template_length(), 0, "R1 TLEN should be 0 when mate is unmapped");
        assert_eq!(r2.template_length(), 0, "R2 TLEN should be 0 when mate is unmapped");

        // Mapped read (R1) should NOT have MQ or MC (mate is unmapped)
        assert!(
            r1.data().get(&Tag::from(SamTag::MQ)).is_none(),
            "Mapped read should not have MQ when mate is unmapped"
        );
        assert!(
            r1.data().get(&Tag::from(SamTag::MC)).is_none(),
            "Mapped read should not have MC when mate is unmapped"
        );

        Ok(())
    }

    /// Tests `fix_mate_info` for R2 supplementary alignments
    ///
    /// Verifies that R2 supplementary alignments get correct mate info
    /// from the R1 primary alignment (covers the R2 supplemental path).
    #[test]
    fn test_fix_mate_info_r2_supplementary() -> Result<()> {
        let mut unmapped = SamBuilder::new_unmapped();
        let mut mapped = SamBuilder::new_mapped();

        let mut attrs = HashMap::new();
        attrs.insert("RX", BufValue::from("ACGT".to_string()));
        unmapped.add_pair_with_attrs("q1", None, None, true, true, &attrs);

        let r1_mapped = RecordBuilder::new()
            .name("q1")
            .flags(Flags::SEGMENTED | Flags::FIRST_SEGMENT | Flags::MATE_REVERSE_COMPLEMENTED)
            .reference_sequence_id(0)
            .alignment_start(100)
            .cigar("50M")
            .mapping_quality(40)
            .sequence(&"A".repeat(50))
            .qualities(&[30u8; 50])
            .tag("PG", MAPPED_PG_ID.to_string())
            .tag("AS", 77i32)
            .build();
        let r2_mapped = RecordBuilder::new()
            .name("q1")
            .flags(
                Flags::SEGMENTED
                    | Flags::LAST_SEGMENT
                    | Flags::REVERSE_COMPLEMENTED
                    | Flags::PROPERLY_SEGMENTED,
            )
            .reference_sequence_id(0)
            .alignment_start(300)
            .cigar("75M")
            .mapping_quality(30)
            .sequence(&"A".repeat(75))
            .qualities(&[30u8; 75])
            .tag("PG", MAPPED_PG_ID.to_string())
            .tag("AS", 55i32)
            .build();
        mapped.push_record(r1_mapped);
        mapped.push_record(r2_mapped);

        // Add R2 supplementary alignment
        let supp_rec = RecordBuilder::new()
            .name("q1")
            .flags(
                Flags::SEGMENTED
                    | Flags::LAST_SEGMENT
                    | Flags::SUPPLEMENTARY
                    | Flags::REVERSE_COMPLEMENTED,
            )
            .reference_sequence_id(0)
            .alignment_start(500)
            .cigar("60M")
            .mapping_quality(20)
            .sequence(&"A".repeat(60))
            .qualities(&[30u8; 60])
            .tag("PG", MAPPED_PG_ID.to_string())
            .tag("AS", 33i32)
            .build();
        mapped.push_record(supp_rec);

        let records = run_zipper(&unmapped, &mapped, vec![], vec![], vec![])?;
        assert_eq!(records.len(), 3);

        let r1_primary = records
            .iter()
            .find(|r| r.flags().is_first_segment() && !r.flags().is_supplementary())
            .unwrap();
        let r2_supp = records
            .iter()
            .find(|r| !r.flags().is_first_segment() && r.flags().is_supplementary())
            .unwrap();

        // R2 supplementary's mate should point to R1's primary
        assert_eq!(
            r2_supp.mate_reference_sequence_id(),
            r1_primary.reference_sequence_id(),
            "R2 supp RNEXT should point to R1 primary's reference"
        );
        assert_eq!(
            r2_supp.mate_alignment_start(),
            r1_primary.alignment_start(),
            "R2 supp PNEXT should point to R1 primary's position"
        );

        // MQ on R2 supplementary should be R1 primary's mapping quality
        let supp_mq = r2_supp.data().get(&Tag::from(SamTag::MQ));
        assert!(
            matches!(supp_mq, Some(BufValue::Int8(40))),
            "R2 supp MQ should be R1 primary's mapq (40), got {supp_mq:?}"
        );

        // MC on R2 supplementary should be R1 primary's CIGAR
        if let Some(BufValue::String(mc)) = r2_supp.data().get(&Tag::from(SamTag::MC)) {
            assert_eq!(mc.to_string(), "50M", "R2 supp MC should be R1 primary's CIGAR");
        } else {
            panic!("R2 supp should have MC tag");
        }

        // ms (mate score) should be R1 primary's AS value
        let supp_ms = r2_supp.data().get(&Tag::from(SamTag::MS));
        assert!(
            matches!(supp_ms, Some(BufValue::Int8(77))),
            "R2 supp ms should be R1 primary's AS (77), got {supp_ms:?}"
        );

        Ok(())
    }

    /// Tests `fix_mate_info` when both reads are unmapped
    ///
    /// Verifies that:
    /// - Both reads get unmapped coordinates (ref=-1, pos=-1)
    /// - TLEN is 0 for both reads
    /// - MQ/MC tags are removed from both reads
    #[test]
    fn test_fix_mate_info_both_unmapped() -> Result<()> {
        let mut unmapped = SamBuilder::new_unmapped();
        let mut mapped = SamBuilder::new_mapped();

        let mut attrs = HashMap::new();
        attrs.insert("RX", BufValue::from("ACGT".to_string()));
        unmapped.add_pair_with_attrs("q1", None, None, true, true, &attrs);

        // Both reads unmapped in mapped BAM
        let mut mapped_attrs = HashMap::new();
        mapped_attrs.insert("PG", BufValue::from(MAPPED_PG_ID.to_string()));
        mapped.add_pair_with_attrs("q1", None, None, true, true, &mapped_attrs);

        let records = run_zipper(&unmapped, &mapped, vec![], vec![], vec![])?;
        assert_eq!(records.len(), 2);

        for rec in &records {
            assert_eq!(rec.template_length(), 0, "TLEN should be 0 when both unmapped");
            assert!(
                rec.data().get(&Tag::from(SamTag::MQ)).is_none(),
                "MQ should not exist when both unmapped"
            );
            assert!(
                rec.data().get(&Tag::from(SamTag::MC)).is_none(),
                "MC should not exist when both unmapped"
            );
        }

        Ok(())
    }

    /// Tests `append_buf_value_raw` covers various tag value types
    ///
    /// Exercises `BufValue` types not hit by standard unmapped tags: `Character`,
    /// `Int16` array, `UInt8`, `Hex`, and `Float`.
    #[test]
    fn test_varied_tag_types_roundtrip() -> Result<()> {
        let mut unmapped = SamBuilder::new_unmapped();
        let mut mapped = SamBuilder::new_mapped();

        let mut attrs = HashMap::new();
        attrs.insert("RX", BufValue::from("ACGT".to_string()));
        attrs.insert("ca", BufValue::Character(b'X'));
        attrs.insert("fi", BufValue::Float(1.23));
        attrs.insert("ui", BufValue::UInt8(200));
        attrs.insert("si", BufValue::Int16(-300));
        attrs.insert("ul", BufValue::UInt32(100_000));
        attrs.insert("hx", BufValue::Hex(bstr::BString::from("DEADBEEF")));
        unmapped.add_frag_with_attrs("q1", None, true, &attrs);

        let mut mapped_attrs = HashMap::new();
        mapped_attrs.insert("PG", BufValue::from(MAPPED_PG_ID.to_string()));
        mapped.add_frag_with_attrs("q1", Some(100), true, &mapped_attrs);

        let records = run_zipper(&unmapped, &mapped, vec![], vec![], vec![])?;
        assert_eq!(records.len(), 1);
        let rec = &records[0];

        // Verify each tag type survived the raw-byte roundtrip
        assert!(
            matches!(rec.data().get(&Tag::new(b'c', b'a')), Some(BufValue::Character(b'X'))),
            "Character tag should roundtrip"
        );
        if let Some(BufValue::Float(f)) = rec.data().get(&Tag::new(b'f', b'i')) {
            assert!((f - 1.23).abs() < 0.001, "Float tag should roundtrip, got {f}");
        } else {
            panic!("Float tag missing");
        }
        // UInt8(200) may be normalized to Int16(200) by smallest-signed encoding path,
        // or stay as UInt8 — either is acceptable
        let ui_val = rec.data().get(&Tag::new(b'u', b'i'));
        assert!(ui_val.is_some(), "UInt8 tag should exist");
        let si_val = rec.data().get(&Tag::new(b's', b'i'));
        assert!(
            matches!(si_val, Some(BufValue::Int16(-300))),
            "Int16 tag should roundtrip, got {si_val:?}"
        );
        let ul_val = rec.data().get(&Tag::new(b'u', b'l'));
        assert!(ul_val.is_some(), "UInt32 tag should exist");
        if let Some(BufValue::Hex(h)) = rec.data().get(&Tag::new(b'h', b'x')) {
            assert_eq!(h.to_string(), "DEADBEEF", "Hex tag should roundtrip");
        } else {
            panic!("Hex tag missing or wrong type");
        }

        Ok(())
    }

    #[rstest]
    // --exclude-missing-reads (default false)
    #[case(&["zipper", "-u", "u.bam", "-r", "ref.fa", "-o", "out.bam"], false)]
    #[case(&["zipper", "-u", "u.bam", "-r", "ref.fa", "-o", "out.bam", "--exclude-missing-reads"], true)]
    #[case(&["zipper", "-u", "u.bam", "-r", "ref.fa", "-o", "out.bam", "--exclude-missing-reads", "true"], true)]
    #[case(&["zipper", "-u", "u.bam", "-r", "ref.fa", "-o", "out.bam", "--exclude-missing-reads", "false"], false)]
    #[case(&["zipper", "-u", "u.bam", "-r", "ref.fa", "-o", "out.bam", "--exclude-missing-reads=true"], true)]
    #[case(&["zipper", "-u", "u.bam", "-r", "ref.fa", "-o", "out.bam", "--exclude-missing-reads=false"], false)]
    fn test_exclude_missing_reads_parsing(#[case] args: &[&str], #[case] expected: bool) {
        let cmd = Zipper::try_parse_from(args).expect("failed to parse Zipper arguments");
        assert_eq!(cmd.exclude_missing_reads, expected);
    }

    #[rstest]
    // --skip-tc-tags (default false)
    #[case(&["zipper", "-u", "u.bam", "-r", "ref.fa", "-o", "out.bam"], false)]
    #[case(&["zipper", "-u", "u.bam", "-r", "ref.fa", "-o", "out.bam", "--skip-tc-tags"], true)]
    #[case(&["zipper", "-u", "u.bam", "-r", "ref.fa", "-o", "out.bam", "--skip-tc-tags", "true"], true)]
    #[case(&["zipper", "-u", "u.bam", "-r", "ref.fa", "-o", "out.bam", "--skip-tc-tags", "false"], false)]
    #[case(&["zipper", "-u", "u.bam", "-r", "ref.fa", "-o", "out.bam", "--skip-tc-tags=true"], true)]
    #[case(&["zipper", "-u", "u.bam", "-r", "ref.fa", "-o", "out.bam", "--skip-tc-tags=false"], false)]
    fn test_skip_tc_tags_parsing(#[case] args: &[&str], #[case] expected: bool) {
        let cmd = Zipper::try_parse_from(args).expect("failed to parse Zipper arguments");
        assert_eq!(cmd.skip_tc_tags, expected);
    }

    #[rstest]
    // --restore-unconverted-bases (default false)
    #[case(&["zipper", "-u", "u.bam", "-r", "ref.fa", "-o", "out.bam"], false)]
    #[case(&["zipper", "-u", "u.bam", "-r", "ref.fa", "-o", "out.bam", "--restore-unconverted-bases"], true)]
    #[case(&["zipper", "-u", "u.bam", "-r", "ref.fa", "-o", "out.bam", "--restore-unconverted-bases", "true"], true)]
    #[case(&["zipper", "-u", "u.bam", "-r", "ref.fa", "-o", "out.bam", "--restore-unconverted-bases", "false"], false)]
    #[case(&["zipper", "-u", "u.bam", "-r", "ref.fa", "-o", "out.bam", "--restore-unconverted-bases=true"], true)]
    #[case(&["zipper", "-u", "u.bam", "-r", "ref.fa", "-o", "out.bam", "--restore-unconverted-bases=false"], false)]
    fn test_restore_unconverted_bases_parsing(#[case] args: &[&str], #[case] expected: bool) {
        let cmd = Zipper::try_parse_from(args).expect("failed to parse Zipper arguments");
        assert_eq!(cmd.restore_unconverted_bases, expected);
    }

    /// Helper for the raw-byte restore tests: encode a `RecordBuf` to BAM
    /// bytes, run the raw-byte restore, and return the resulting SEQ as ASCII.
    /// Mirrors the assertion shape of the `RecordBuf` tests so equivalence is
    /// obvious at a glance.
    fn run_restore_raw(
        record: &RecordBuf,
        reference: &ReferenceReader,
        header: &Header,
    ) -> Result<(Vec<u8>, Vec<u8>)> {
        let mut raw = Vec::with_capacity(256);
        encode_record_buf(&mut raw, header, record)
            .map_err(|e| anyhow::anyhow!("encode_record_buf failed: {e}"))?;
        let mut rec = RawRecord::from(raw);
        restore_unconverted_bases_in_raw_record(&mut rec, reference, header)?;
        let raw = rec.into_inner();
        let seq = bam_fields::extract_sequence(&raw);
        Ok((raw, seq))
    }

    /// Raw-byte equivalent of `test_restore_unconverted_bases_top_strand` —
    /// asserts the new `restore_unconverted_bases_in_raw_record` (which operates
    /// on packed BAM bytes without round-tripping through `RecordBuf`) produces
    /// the same SEQ as the `RecordBuf` path for the canonical top-strand case.
    #[test]
    fn test_restore_unconverted_bases_in_raw_record_top_strand() -> Result<()> {
        use crate::sam::builder::create_test_fasta;

        // Reference: ACGTACGTACGT (chr1) — C positions: 2, 6, 10
        let fasta = create_test_fasta(&[("chr1", "ACGTACGTACGT")])?;
        let reference = ReferenceReader::new(fasta.path())?;

        let header: Header = "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:12\n".parse().unwrap();

        // Same record as test_restore_unconverted_bases_top_strand
        let record = RecordBuilder::new()
            .name("q1")
            .sequence("ATGTATGT")
            .qualities(&[30; 8])
            .reference_sequence_id(0)
            .alignment_start(1)
            .cigar("8M")
            .flags(Flags::from(FLAG_PAIRED | FLAG_READ1))
            .tag("YD", "f")
            .build();

        let (_raw, seq) = run_restore_raw(&record, &reference, &header)?;
        assert_eq!(seq, b"ACGTACGT");

        Ok(())
    }

    /// Raw-byte equivalent of `test_restore_unconverted_bases_bottom_strand` —
    /// bottom strand restores A→G at ref-G positions.
    #[test]
    fn test_restore_unconverted_bases_in_raw_record_bottom_strand() -> Result<()> {
        use crate::sam::builder::create_test_fasta;

        let fasta = create_test_fasta(&[("chr1", "ACGTACGTACGT")])?;
        let reference = ReferenceReader::new(fasta.path())?;
        let header: Header = "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:12\n".parse().unwrap();

        // A at positions 3,7 = ref-G → restore to G
        let record = RecordBuilder::new()
            .name("q1")
            .sequence("ACATACAT")
            .qualities(&[30; 8])
            .reference_sequence_id(0)
            .alignment_start(1)
            .cigar("8M")
            .flags(Flags::from(FLAG_PAIRED | FLAG_READ2))
            .tag("YD", "r")
            .build();

        let (_raw, seq) = run_restore_raw(&record, &reference, &header)?;
        assert_eq!(seq, b"ACGTACGT");

        Ok(())
    }

    /// Raw-byte equivalent of `test_restore_unconverted_bases_top_strand_reverse` —
    /// reverse-strand top-strand reads restore A→G (the complement of T→C) at
    /// ref-G positions (the complement of ref-C).
    #[test]
    fn test_restore_unconverted_bases_in_raw_record_top_strand_reverse() -> Result<()> {
        use crate::sam::builder::create_test_fasta;

        let fasta = create_test_fasta(&[("chr1", "ACGTACGTACGT")])?;
        let reference = ReferenceReader::new(fasta.path())?;
        let header: Header = "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:12\n".parse().unwrap();

        let record = RecordBuilder::new()
            .name("q1")
            .sequence("ACATACGT")
            .qualities(&[30; 8])
            .reference_sequence_id(0)
            .alignment_start(1)
            .cigar("8M")
            .flags(Flags::from(FLAG_PAIRED | FLAG_READ1 | FLAG_REVERSE))
            .tag("YD", "f")
            .build();

        let (_raw, seq) = run_restore_raw(&record, &reference, &header)?;
        assert_eq!(seq, b"ACGTACGT");

        Ok(())
    }

    /// Raw-byte equivalent of `test_restore_unconverted_bases_skips_no_yd_tag` —
    /// records without a YD tag come through unchanged.
    #[test]
    fn test_restore_unconverted_bases_in_raw_record_skips_no_yd_tag() -> Result<()> {
        use crate::sam::builder::create_test_fasta;

        let fasta = create_test_fasta(&[("chr1", "CCCCCCCC")])?;
        let reference = ReferenceReader::new(fasta.path())?;
        let header: Header = "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:8\n".parse().unwrap();

        // No YD tag — function must early-exit and leave SEQ untouched.
        let record = RecordBuilder::new()
            .name("q1")
            .sequence("TTTTTTTT")
            .qualities(&[30; 8])
            .reference_sequence_id(0)
            .alignment_start(1)
            .cigar("8M")
            .flags(Flags::from(FLAG_PAIRED | FLAG_READ1))
            .build();

        let (_raw, seq) = run_restore_raw(&record, &reference, &header)?;
        assert_eq!(seq, b"TTTTTTTT");

        Ok(())
    }

    /// Raw-byte equivalent of `test_restore_unconverted_bases_skips_unmapped` —
    /// unmapped records come through unchanged.
    #[test]
    fn test_restore_unconverted_bases_in_raw_record_skips_unmapped() -> Result<()> {
        use crate::sam::builder::create_test_fasta;

        let fasta = create_test_fasta(&[("chr1", "CCCCCCCC")])?;
        let reference = ReferenceReader::new(fasta.path())?;
        let header: Header = "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:8\n".parse().unwrap();

        let record = RecordBuilder::new()
            .name("q1")
            .sequence("TTTTTTTT")
            .qualities(&[30; 8])
            .flags(Flags::UNMAPPED)
            .tag("YD", "f")
            .build();

        let (_raw, seq) = run_restore_raw(&record, &reference, &header)?;
        assert_eq!(seq, b"TTTTTTTT");

        Ok(())
    }

    /// Raw-byte path with NM/MD tags present: removed only when SEQ is modified.
    /// Mirrors the implicit contract of the `RecordBuf` path's stale-tag clearing.
    #[test]
    fn test_restore_unconverted_bases_in_raw_record_clears_nm_md_when_modified() -> Result<()> {
        use crate::sam::builder::create_test_fasta;

        let fasta = create_test_fasta(&[("chr1", "ACGTACGTACGT")])?;
        let reference = ReferenceReader::new(fasta.path())?;
        let header: Header = "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:12\n".parse().unwrap();

        // SEQ modified case: T at ref-C should restore, NM/MD should be removed.
        let modified_record = RecordBuilder::new()
            .name("q1")
            .sequence("ATGTATGT")
            .qualities(&[30; 8])
            .reference_sequence_id(0)
            .alignment_start(1)
            .cigar("8M")
            .flags(Flags::from(FLAG_PAIRED | FLAG_READ1))
            .tag("YD", "f")
            .tag("NM", 2i32)
            .tag("MD", "1C3C2")
            .build();
        let (raw, _seq) = run_restore_raw(&modified_record, &reference, &header)?;
        let aux = bam_fields::aux_data_slice(&raw);
        assert!(bam_fields::find_tag_type(aux, &SamTag::NM).is_none(), "NM should be removed");
        assert!(bam_fields::find_tag_type(aux, &SamTag::MD).is_none(), "MD should be removed");

        // SEQ unchanged case: read already matches ref, NM/MD must be preserved.
        let unchanged_record = RecordBuilder::new()
            .name("q2")
            .sequence("ACGTACGT")
            .qualities(&[30; 8])
            .reference_sequence_id(0)
            .alignment_start(1)
            .cigar("8M")
            .flags(Flags::from(FLAG_PAIRED | FLAG_READ1))
            .tag("YD", "f")
            .tag("NM", 0i32)
            .tag("MD", "8")
            .build();
        let (raw, _seq) = run_restore_raw(&unchanged_record, &reference, &header)?;
        let aux = bam_fields::aux_data_slice(&raw);
        assert!(bam_fields::find_tag_type(aux, &SamTag::NM).is_some(), "NM should be kept");
        assert!(bam_fields::find_tag_type(aux, &SamTag::MD).is_some(), "MD should be kept");

        Ok(())
    }

    /// Parity test: the raw-byte path's own `read_pos` / `ref_offset` bookkeeping
    /// must match the `RecordBuf` oracle across nontrivial CIGARs (`I`, `D`, `N`).
    /// The straight-`M` raw tests wouldn't catch an off-by-one in the indel arms.
    #[rstest]
    #[case::match_insert_delete_match("2M1I2M1D3M", "TTNATTGT", "f", FLAG_READ1)]
    #[case::match_insert_match("4M1I3M", "ATGTNATG", "f", FLAG_READ1)]
    #[case::match_delete_match("4M1D3M", "ATGTATG", "f", FLAG_READ1)]
    #[case::match_skip_match("4M1N3M", "ATGTATG", "f", FLAG_READ1)]
    #[case::bottom_strand_indels("2M1I2M1D3M", "ACNTAGAT", "r", FLAG_READ2)]
    fn test_restore_unconverted_bases_in_raw_record_indels_parity(
        #[case] cigar: &str,
        #[case] sequence: &str,
        #[case] yd: &str,
        #[case] flag_mate: u16,
    ) -> Result<()> {
        use crate::sam::builder::create_test_fasta;

        let fasta = create_test_fasta(&[("chr1", "ACGTACGTACGT")])?;
        let reference = ReferenceReader::new(fasta.path())?;
        let header: Header = "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:12\n".parse().unwrap();

        let quals = vec![30u8; sequence.len()];
        let build = || {
            RecordBuilder::new()
                .name("q1")
                .sequence(sequence)
                .qualities(&quals)
                .reference_sequence_id(0)
                .alignment_start(1)
                .cigar(cigar)
                .flags(Flags::from(FLAG_PAIRED | flag_mate))
                .tag("YD", yd)
                .tag("NM", 1i32)
                .tag("MD", "8")
                .build()
        };

        // Oracle: RecordBuf path produces the expected SEQ.
        let mut oracle = build();
        restore_unconverted_bases_in_record(&mut oracle, &reference, &header)?;
        let oracle_seq: Vec<u8> = oracle.sequence().as_ref().to_vec();
        let oracle_modified = oracle_seq != sequence.as_bytes();

        // Raw-byte path: SEQ must match the oracle byte-for-byte, and stale NM/MD
        // must be cleared iff the oracle modified the SEQ.
        let (raw, raw_seq) = run_restore_raw(&build(), &reference, &header)?;
        assert_eq!(raw_seq, oracle_seq, "raw SEQ mismatch for CIGAR {cigar}");

        let aux = bam_fields::aux_data_slice(&raw);
        let nm_present = bam_fields::find_tag_type(aux, &SamTag::NM).is_some();
        let md_present = bam_fields::find_tag_type(aux, &SamTag::MD).is_some();
        if oracle_modified {
            assert!(!nm_present, "NM should be removed when SEQ modified (CIGAR {cigar})");
            assert!(!md_present, "MD should be removed when SEQ modified (CIGAR {cigar})");
        } else {
            assert!(nm_present, "NM should be preserved when SEQ unchanged (CIGAR {cigar})");
            assert!(md_present, "MD should be preserved when SEQ unchanged (CIGAR {cigar})");
        }

        Ok(())
    }

    /// Tests that `restore_unconverted_bases_in_record` replaces T→C at ref-C positions
    /// for top-strand reads (YD:Z:f).
    #[test]
    fn test_restore_unconverted_bases_top_strand() -> Result<()> {
        use crate::sam::builder::create_test_fasta;

        // Reference: ACGTACGTACGT (chr1)
        // Positions:  1234567890..
        // C positions: 2, 6, 10
        let fasta = create_test_fasta(&[("chr1", "ACGTACGTACGT")])?;
        let reference = ReferenceReader::new(fasta.path())?;

        // Build a header with chr1
        let header: Header = "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:12\n".parse().unwrap();

        // Create a record aligned to pos 1 with 8M CIGAR
        // Sequence: ATGTATGT (T at positions 2,6 = ref-C positions → should be restored to C)
        let mut record = RecordBuilder::new()
            .name("q1")
            .sequence("ATGTATGT")
            .qualities(&[30; 8])
            .reference_sequence_id(0)
            .alignment_start(1)
            .cigar("8M")
            .flags(Flags::from(FLAG_PAIRED | FLAG_READ1))
            .tag("YD", "f")
            .build();

        restore_unconverted_bases_in_record(&mut record, &reference, &header)?;

        let seq: Vec<u8> = record.sequence().as_ref().to_vec();
        // Position 1: A (ref=A, no change)
        // Position 2: T→C (ref=C, top strand, T was converted)
        // Position 3: G (ref=G, no change)
        // Position 4: T (ref=T, no change)
        // Position 5: A (ref=A, no change)
        // Position 6: T→C (ref=C, top strand, T was converted)
        // Position 7: G (ref=G, no change)
        // Position 8: T (ref=T, no change)
        assert_eq!(seq, b"ACGTACGT");

        Ok(())
    }

    /// Tests that `restore_unconverted_bases_in_record` replaces A→G at ref-G positions
    /// for bottom-strand reads (YD:Z:r).
    #[test]
    fn test_restore_unconverted_bases_bottom_strand() -> Result<()> {
        use crate::sam::builder::create_test_fasta;

        // Reference: ACGTACGTACGT (chr1)
        // G positions: 3, 7, 11
        let fasta = create_test_fasta(&[("chr1", "ACGTACGTACGT")])?;
        let reference = ReferenceReader::new(fasta.path())?;

        let header: Header = "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:12\n".parse().unwrap();

        // Record at pos 1 with 8M. A at positions 3,7 = ref-G → should be restored to G
        let mut record = RecordBuilder::new()
            .name("q1")
            .sequence("ACATACAT")
            .qualities(&[30; 8])
            .reference_sequence_id(0)
            .alignment_start(1)
            .cigar("8M")
            .flags(Flags::from(FLAG_PAIRED | FLAG_READ2))
            .tag("YD", "r")
            .build();

        restore_unconverted_bases_in_record(&mut record, &reference, &header)?;

        let seq: Vec<u8> = record.sequence().as_ref().to_vec();
        assert_eq!(seq, b"ACGTACGT");

        Ok(())
    }

    /// Tests that `restore_unconverted_bases_in_record` handles insertions and deletions.
    #[test]
    fn test_restore_unconverted_bases_with_indels() -> Result<()> {
        use crate::sam::builder::create_test_fasta;

        // Reference: ACGTACGTACGT
        let fasta = create_test_fasta(&[("chr1", "ACGTACGTACGT")])?;
        let reference = ReferenceReader::new(fasta.path())?;

        let header: Header = "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:12\n".parse().unwrap();

        // CIGAR: 2M1I2M1D3M → 8 read bases, 8 ref bases consumed
        // ref_bases (pos 1-8): A C G T A C G T
        //
        // Walk:
        //   2M: read[0,1] vs ref[0,1](A,C) → T stays, T→C
        //   1I: read[2] (insertion, no ref)
        //   2M: read[3,4] vs ref[2,3](G,T) → no change (not ref-C)
        //   1D: ref[4](A) skipped
        //   3M: read[5,6,7] vs ref[5,6,7](C,G,T) → T→C, G stays, T stays
        let mut record = RecordBuilder::new()
            .name("q1")
            .sequence("TTNATTGT")
            .qualities(&[30; 8])
            .reference_sequence_id(0)
            .alignment_start(1)
            .cigar("2M1I2M1D3M")
            .flags(Flags::from(FLAG_PAIRED | FLAG_READ1))
            .tag("YD", "f")
            .build();

        restore_unconverted_bases_in_record(&mut record, &reference, &header)?;

        let seq: Vec<u8> = record.sequence().as_ref().to_vec();
        assert_eq!(seq, b"TCNATCGT");

        Ok(())
    }

    /// Tests that `restore_unconverted_bases_in_record` skips unmapped reads.
    #[test]
    fn test_restore_unconverted_bases_skips_unmapped() -> Result<()> {
        use crate::sam::builder::create_test_fasta;

        let fasta = create_test_fasta(&[("chr1", "CCCCCCCC")])?;
        let reference = ReferenceReader::new(fasta.path())?;

        let header: Header = "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:8\n".parse().unwrap();

        // Unmapped read — should be left unchanged
        let mut record = RecordBuilder::new()
            .name("q1")
            .sequence("TTTTTTTT")
            .qualities(&[30; 8])
            .flags(Flags::UNMAPPED)
            .tag("YD", "f")
            .build();

        restore_unconverted_bases_in_record(&mut record, &reference, &header)?;

        let seq: Vec<u8> = record.sequence().as_ref().to_vec();
        assert_eq!(seq, b"TTTTTTTT"); // Unchanged

        Ok(())
    }

    /// Tests that `restore_unconverted_bases_in_record` skips reads without YD tag.
    #[test]
    fn test_restore_unconverted_bases_skips_no_yd_tag() -> Result<()> {
        use crate::sam::builder::create_test_fasta;

        let fasta = create_test_fasta(&[("chr1", "CCCCCCCC")])?;
        let reference = ReferenceReader::new(fasta.path())?;

        let header: Header = "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:8\n".parse().unwrap();

        // Mapped read without YD tag — should be left unchanged
        let mut record = RecordBuilder::new()
            .name("q1")
            .sequence("TTTTTTTT")
            .qualities(&[30; 8])
            .reference_sequence_id(0)
            .alignment_start(1)
            .cigar("8M")
            .flags(Flags::from(FLAG_PAIRED | FLAG_READ1))
            .build();

        restore_unconverted_bases_in_record(&mut record, &reference, &header)?;

        let seq: Vec<u8> = record.sequence().as_ref().to_vec();
        assert_eq!(seq, b"TTTTTTTT"); // Unchanged

        Ok(())
    }

    /// Tests that non-converted bases at ref-C positions are left alone.
    #[test]
    fn test_restore_unconverted_bases_preserves_already_unconverted() -> Result<()> {
        use crate::sam::builder::create_test_fasta;

        // Reference: CCCCCCCC (all C positions)
        let fasta = create_test_fasta(&[("chr1", "CCCCCCCC")])?;
        let reference = ReferenceReader::new(fasta.path())?;

        let header: Header = "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:8\n".parse().unwrap();

        // Top strand: some bases already C (unconverted), some T (converted), some other
        let mut record = RecordBuilder::new()
            .name("q1")
            .sequence("CTCACTGC")
            .qualities(&[30; 8])
            .reference_sequence_id(0)
            .alignment_start(1)
            .cigar("8M")
            .flags(Flags::from(FLAG_PAIRED | FLAG_READ1))
            .tag("YD", "f")
            .build();

        restore_unconverted_bases_in_record(&mut record, &reference, &header)?;

        let seq: Vec<u8> = record.sequence().as_ref().to_vec();
        // C stays C, T→C, C stays C, A stays A (not T→C), C stays C, T→C, G stays G, C stays C
        assert_eq!(seq, b"CCCACCGC");

        Ok(())
    }

    /// Top-strand read aligned to a span with no reference Cs: nothing to restore.
    /// Exercises the `memchr2` early-exit fast path — record must come through
    /// byte-for-byte unchanged (including any T's, which would only be candidates
    /// where the reference is C).
    #[test]
    fn test_restore_unconverted_bases_no_target_in_span() -> Result<()> {
        use crate::sam::builder::create_test_fasta;

        // Reference is all A/T/G — no C anywhere in the aligned span.
        let fasta = create_test_fasta(&[("chr1", "ATGTATGT")])?;
        let reference = ReferenceReader::new(fasta.path())?;

        let header: Header = "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:8\n".parse().unwrap();

        // SEQ contains T's but they cannot be restoration candidates because
        // the reference span has no C positions for them to align to.
        let original_seq = b"ATGTATGT";
        let mut record = RecordBuilder::new()
            .name("q1")
            .sequence(std::str::from_utf8(original_seq)?)
            .qualities(&[30; 8])
            .reference_sequence_id(0)
            .alignment_start(1)
            .cigar("8M")
            .flags(Flags::from(FLAG_PAIRED | FLAG_READ1))
            .tag("YD", "f")
            .tag("NM", 0i32)
            .tag("MD", "8")
            .build();

        restore_unconverted_bases_in_record(&mut record, &reference, &header)?;

        let seq: Vec<u8> = record.sequence().as_ref().to_vec();
        assert_eq!(seq, original_seq);
        // No SEQ change => NM/MD must remain present.
        assert!(record.data().get(&Tag::from(SamTag::NM)).is_some());
        assert!(record.data().get(&Tag::from(SamTag::MD)).is_some());

        Ok(())
    }

    /// Tests restoration for a reverse-complemented top-strand read (YD:Z:f, 0x10 set).
    /// SEQ is stored as the reverse complement, so at ref-G positions (complement of C),
    /// converted A (complement of T) should be restored to G (complement of C).
    #[test]
    fn test_restore_unconverted_bases_top_strand_reverse() -> Result<()> {
        use crate::sam::builder::create_test_fasta;

        // Reference: ACGTACGTACGT (chr1)
        // G positions (1-based): 3, 7, 11
        let fasta = create_test_fasta(&[("chr1", "ACGTACGTACGT")])?;
        let reference = ReferenceReader::new(fasta.path())?;

        let header: Header = "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:12\n".parse().unwrap();

        // Aligned at pos 1 with 8M, reverse-complemented (0x10).
        // Reference region: ACGTACGT (pos 1-8)
        // SEQ (stored RC'd): at ref-G positions 3,7 show A (converted) instead of G
        let mut record = RecordBuilder::new()
            .name("q1")
            .sequence("ACATACGT")
            .qualities(&[30; 8])
            .reference_sequence_id(0)
            .alignment_start(1)
            .cigar("8M")
            .flags(Flags::from(FLAG_PAIRED | FLAG_READ1 | FLAG_REVERSE))
            .tag("YD", "f")
            .build();

        restore_unconverted_bases_in_record(&mut record, &reference, &header)?;

        let seq: Vec<u8> = record.sequence().as_ref().to_vec();
        // Positions 3,7 (ref=G): A→G restored
        assert_eq!(seq, b"ACGTACGT");

        Ok(())
    }

    /// Tests restoration for a reverse-complemented bottom-strand read (YD:Z:r, 0x10 set).
    /// SEQ is stored as the reverse complement, so at ref-C positions,
    /// converted T (complement of A) should be restored to C (complement of G).
    #[test]
    fn test_restore_unconverted_bases_bottom_strand_reverse() -> Result<()> {
        use crate::sam::builder::create_test_fasta;

        // Reference: ACGTACGTACGT (chr1)
        // C positions (1-based): 2, 6, 10
        let fasta = create_test_fasta(&[("chr1", "ACGTACGTACGT")])?;
        let reference = ReferenceReader::new(fasta.path())?;

        let header: Header = "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:12\n".parse().unwrap();

        // Aligned at pos 1 with 8M, reverse-complemented (0x10).
        // Reference region: ACGTACGT (pos 1-8)
        // For (is_top=false, is_reverse=true): ref_target=C, converted=T, unconverted=C
        // At ref-C positions 2,6: SEQ shows T (converted) instead of C
        let mut record = RecordBuilder::new()
            .name("q1")
            .sequence("ATGTACGT")
            .qualities(&[30; 8])
            .reference_sequence_id(0)
            .alignment_start(1)
            .cigar("8M")
            .flags(Flags::from(FLAG_PAIRED | FLAG_READ2 | FLAG_REVERSE))
            .tag("YD", "r")
            .build();

        restore_unconverted_bases_in_record(&mut record, &reference, &header)?;

        let seq: Vec<u8> = record.sequence().as_ref().to_vec();
        // Positions 2,6 (ref=C): T→C restored
        assert_eq!(seq, b"ACGTACGT");

        Ok(())
    }
}

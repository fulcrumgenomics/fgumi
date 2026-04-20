//! Shared types and functions for metrics commands (duplex-metrics, simplex-metrics).
//!
//! This module contains the common infrastructure for reading UMI-grouped BAM files,
//! grouping templates by coordinate, filtering by genomic intervals, and performing
//! deterministic downsampling. Both `duplex_metrics` and `simplex_metrics` commands
//! build on these shared primitives.

use crate::bam_io::create_raw_bam_reader;
use crate::progress::ProgressTracker;
use crate::sam::SamTag;
use crate::template::TemplateIterator;
use anyhow::{Context, Result};
use fgumi_raw_bam::raw_record_to_record_buf;

use log::info;
use murmur3::murmur3_32;
use noodles::sam::alignment::record::Cigar;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::RecordBuf;
use std::path::Path;
use std::sync::OnceLock;

/// Standard downsampling fractions: 5%, 10%, 15%, ..., 100%.
pub const DOWNSAMPLING_FRACTIONS: [f64; 20] = [
    0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80,
    0.85, 0.90, 0.95, 1.00,
];

/// Cached R availability check (computed once per process).
static R_AVAILABLE: OnceLock<bool> = OnceLock::new();

/// Genomic interval for filtering, stored as 0-based half-open coordinates.
#[derive(Clone, Debug)]
pub struct Interval {
    /// Reference sequence name (e.g. "chr1").
    pub ref_name: String,
    /// 0-based start position (inclusive).
    pub start: i32,
    /// 0-based end position (exclusive).
    pub end: i32,
}

/// Read name and template information for downsampling.
#[derive(Clone)]
pub struct TemplateInfo {
    /// Molecular identifier tag value (e.g. "1/A").
    pub mi: String,
    /// Raw UMI tag value (e.g. "AAA-TTT").
    pub rx: String,
    /// Reference sequence name, if mapped.
    pub ref_name: Option<String>,
    /// Alignment start position (1-based), if mapped.
    pub position: Option<i32>,
    /// Alignment end position (1-based), if mapped.
    pub end_position: Option<i32>,
    /// Hash fraction for deterministic downsampling (computed once per template).
    pub hash_fraction: f64,
}

/// Grouping key matching fgbio's `ReadInfo` structure.
///
/// Fields are ordered so the earlier-mapping read comes first.
#[derive(Clone, PartialEq, Eq, Hash)]
pub struct ReadInfoKey {
    /// Reference sequence index for read 1.
    pub ref_index1: Option<usize>,
    /// Unclipped 5' position for read 1.
    pub start1: i32,
    /// `true` if read 1 is reverse-complemented.
    pub strand1: bool,
    /// Reference sequence index for read 2.
    pub ref_index2: Option<usize>,
    /// Unclipped 5' position for read 2.
    pub start2: i32,
    /// `true` if read 2 is reverse-complemented.
    pub strand2: bool,
}

/// Pre-computed metadata for a template within a coordinate group.
pub struct TemplateMetadata<'a> {
    /// Reference to the underlying template info.
    pub template: &'a TemplateInfo,
    /// MI tag value with strand suffix stripped (e.g. "1" from "1/A").
    pub base_umi: &'a str,
    /// `true` if this template belongs to the A strand.
    pub is_a_strand: bool,
    /// `true` if this template belongs to the B strand.
    pub is_b_strand: bool,
}

/// Computes the unclipped 5' position for a read, matching fgbio's `positionOf`.
///
/// For forward strand reads: `unclippedStart = alignmentStart - leading soft clips`
/// For reverse strand reads: `unclippedEnd = alignmentEnd + trailing soft clips`
pub fn unclipped_five_prime_position(record: &RecordBuf) -> Option<i32> {
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

/// Computes a hash value normalized to the [0, 1] range using Murmur3.
///
/// Hashes the read name using the Murmur3 32-bit algorithm with seed 42, then
/// normalizes the result to a floating-point value between 0 and 1. This is used
/// for deterministic downsampling where reads are assigned to fractions based on
/// their hash value. Matches the Scala implementation's hashing approach.
///
/// # Arguments
///
/// * `read_name` - The read name string to hash
///
/// # Returns
///
/// A floating-point value in the range [0, 1].
pub fn compute_hash_fraction(read_name: &str) -> f64 {
    // Use Murmur3 with seed 42 (matching Scala implementation)
    let hash = murmur3_32(&mut std::io::Cursor::new(read_name.as_bytes()), 42).unwrap_or(0);

    // Scala implementation uses Int (i32), takes absolute value, then normalizes
    // Convert u32 to i32 first to match Scala's behavior
    let positive_hash = (hash as i32).unsigned_abs() as f64;
    positive_hash / i32::MAX as f64
}

/// Parses an intervals file in BED or Picard interval list format.
///
/// Auto-detects the format: if any line starts with `@`, the file is treated as a
/// Picard interval list (1-based closed coordinates with a SAM header); otherwise
/// it is treated as BED (0-based half-open coordinates).
///
/// Intervals are stored internally using BED conventions (0-based half-open).
///
/// # Errors
///
/// Returns an error if the file cannot be read or lines cannot be parsed.
pub fn parse_intervals(path: &Path) -> Result<Vec<Interval>> {
    use std::fs::File;
    use std::io::{BufRead, BufReader};

    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut intervals = Vec::new();
    let mut is_interval_list = false;

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();

        // Skip empty lines and comments
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        // Skip SAM header lines (interval list format)
        if line.starts_with('@') {
            is_interval_list = true;
            continue;
        }

        let mut fields = line.splitn(4, '\t');
        let ref_name = fields.next().expect("splitn always yields at least one element");
        let start_str = fields.next();
        let end_str = fields.next();

        let (Some(start_str), Some(end_str)) = (start_str, end_str) else {
            let fmt = if is_interval_list { "interval list" } else { "BED" };
            anyhow::bail!("Invalid {fmt} line (needs at least 3 fields): {line}");
        };

        if is_interval_list {
            // Picard interval list: chr start end strand name (1-based, closed)
            let start: i32 = start_str
                .parse::<i32>()
                .map_err(|_| anyhow::anyhow!("Invalid start position: {start_str}"))?
                - 1; // Convert 1-based to 0-based
            let end: i32 =
                end_str.parse().map_err(|_| anyhow::anyhow!("Invalid end position: {end_str}"))?;
            // end stays the same: 1-based closed end == 0-based half-open end
            intervals.push(Interval { ref_name: ref_name.to_string(), start, end });
        } else {
            // BED format: chr start end [name] [score] [strand] (0-based, half-open)
            let start: i32 = start_str
                .parse()
                .map_err(|_| anyhow::anyhow!("Invalid start position: {start_str}"))?;
            let end: i32 =
                end_str.parse().map_err(|_| anyhow::anyhow!("Invalid end position: {end_str}"))?;
            intervals.push(Interval { ref_name: ref_name.to_string(), start, end });
        }
    }

    Ok(intervals)
}

/// Checks if a template's insert overlaps any provided interval.
///
/// Determines whether a template's genomic insert coordinates overlap with any
/// of the specified intervals. An insert overlaps an interval if any part of it
/// (from start to end position) overlaps the interval region. If no intervals are
/// provided, all templates are considered to overlap (no filtering).
///
/// # Arguments
///
/// * `template` - Template information including chromosome and positions
/// * `intervals` - Slice of intervals to check for overlap
///
/// # Returns
///
/// `true` if the template overlaps any interval or if no intervals are provided,
/// `false` if the template is unmapped or does not overlap any interval.
pub fn overlaps_intervals(template: &TemplateInfo, intervals: &[Interval]) -> bool {
    if intervals.is_empty() {
        return true; // No filtering if no intervals provided
    }

    if let (Some(ref_name), Some(start), Some(end)) =
        (&template.ref_name, template.position, template.end_position)
    {
        // Intervals are 0-based half-open; template positions are 1-based inclusive.
        // In 0-based half-open the template is [start-1, end), so the overlap
        // test is: (start-1) < interval.end && interval.start < end
        // which simplifies to: start <= interval.end && interval.start < end
        intervals.iter().any(|interval| {
            interval.ref_name == *ref_name && start <= interval.end && interval.start < end
        })
    } else {
        false // Unmapped reads or reads without proper coordinates don't overlap any interval
    }
}

/// Validates that a BAM file is not a consensus BAM.
///
/// Consensus BAMs (output from simplex/duplex callers) should not be used with
/// metrics tools. This checks the first valid R1 record and errors if it contains
/// consensus tags.
///
/// # Errors
///
/// Returns an error if the BAM file cannot be read or if it appears to be a consensus BAM.
pub fn validate_not_consensus_bam(input: &Path) -> Result<()> {
    use crate::consensus_tags::is_consensus;
    use crate::sort::bam_fields;
    use fgumi_raw_bam::{RawRecord, RawRecordView};

    let (mut reader, header) = create_raw_bam_reader(input, 1)?;

    // Look at the first valid R1 record
    let mut raw = RawRecord::new();
    loop {
        let n = reader.read_record(&mut raw).context("failed to read BAM record")?;
        if n == 0 {
            break; // EOF
        }

        // Only check R1 records that are paired and primary (raw flag checks).
        // Do not skip UNMAPPED: consensus BAMs are documented as unaligned, so
        // excluding unmapped records here would let consensus BAMs slip past this
        // guard unchecked.
        let flags = RawRecordView::new(&raw).flags();
        if (flags & bam_fields::flags::PAIRED) == 0
            || (flags & bam_fields::flags::FIRST_SEGMENT) == 0
            || (flags & bam_fields::flags::SECONDARY) != 0
            || (flags & bam_fields::flags::SUPPLEMENTARY) != 0
        {
            continue;
        }

        // Decode to RecordBuf for consensus-tag check and name extraction
        let record = raw_record_to_record_buf(&raw, &header)?;

        // Check if this is a consensus read
        if is_consensus(&record) {
            anyhow::bail!(
                "Input BAM file ({}) appears to contain consensus sequences. \
                This metrics tool cannot run on consensus BAMs, and instead requires \
                the UMI-grouped BAM generated by group which is run prior to consensus calling.\n\
                First R1 record '{}' has consensus SAM tags present.",
                input.display(),
                record.name().map_or_else(
                    || "<unnamed>".to_string(),
                    |n| String::from_utf8_lossy(n.as_ref()).to_string()
                )
            );
        }

        // Only need to check the first valid R1 record
        break;
    }

    Ok(())
}

/// Checks if R and required packages (ggplot2, scales) are available.
///
/// Result is cached for the lifetime of the process to avoid repeated subprocess spawns.
pub fn is_r_available() -> bool {
    use std::process::Command;

    *R_AVAILABLE.get_or_init(|| {
        Command::new("Rscript")
            .args(["-e", "stopifnot(require(ggplot2)); stopifnot(require(scales))"])
            .output()
            .map(|output| output.status.success())
            .unwrap_or(false)
    })
}

/// Executes an R script with the given arguments.
///
/// The R script content is written to a temporary file for execution. This ensures
/// the script is always available regardless of working directory or installation
/// location.
///
/// # Arguments
///
/// * `r_script_content` - The R script source code to execute
/// * `args` - Command-line arguments to pass to the R script
/// * `temp_file_name` - Base name for the temporary R script file
///
/// # Errors
///
/// Returns an error if the script cannot be written or R execution fails.
pub fn execute_r_script(r_script_content: &str, args: &[&str], temp_file_name: &str) -> Result<()> {
    use std::process::Command;

    // Write embedded R script to temp file
    let temp_dir = std::env::temp_dir();
    let r_script_path = temp_dir.join(temp_file_name);
    std::fs::write(&r_script_path, r_script_content)
        .context("Failed to write embedded R script to temp file")?;

    info!("Executing R script to generate PDF plots...");

    let output = Command::new("Rscript")
        .arg(&r_script_path)
        .args(args)
        .output()
        .context("Failed to execute Rscript command")?;

    // Clean up temp file (ignore errors)
    let _ = std::fs::remove_file(&r_script_path);

    if output.status.success() {
        Ok(())
    } else {
        let stderr = String::from_utf8_lossy(&output.stderr);
        anyhow::bail!(
            "R script execution failed with exit code {:?}. Error: {}",
            output.status.code(),
            stderr
        )
    }
}

/// Pre-computes metadata for each template in a coordinate group.
///
/// Parses the MI tag to determine strand assignment and extract the base UMI
/// (MI value without the `/A` or `/B` suffix).
pub fn compute_template_metadata(group: &[TemplateInfo]) -> Vec<TemplateMetadata<'_>> {
    group
        .iter()
        .map(|t| {
            let (base_umi, is_a, is_b) = if t.mi.ends_with("/A") {
                (&t.mi[..t.mi.len() - 2], true, false)
            } else if t.mi.ends_with("/B") {
                (&t.mi[..t.mi.len() - 2], false, true)
            } else {
                (t.mi.as_str(), false, false)
            };
            TemplateMetadata { template: t, base_umi, is_a_strand: is_a, is_b_strand: is_b }
        })
        .collect()
}

/// Reads a BAM file, groups templates by [`ReadInfoKey`], and calls a closure for each group.
///
/// This is the shared BAM processing loop used by both duplex-metrics and simplex-metrics.
/// Templates are streamed in coordinate order; when the [`ReadInfoKey`] changes, the
/// accumulated group is dispatched to the closure.
///
/// # Arguments
///
/// * `input` - Path to the input BAM file
/// * `intervals` - Intervals for filtering templates (empty = no filtering)
/// * `num_fractions` - Number of downsampling fractions (used to size the counts vector)
/// * `process_group` - Closure called for each coordinate group with `(group, fraction_counts)`
///
/// The SAM spec standard tags `MI` and `RX` are always used.
///
/// # Returns
///
/// A tuple of `(total_template_count, per_fraction_template_counts)`.
///
/// # Errors
///
/// Returns an error if the BAM file cannot be read, if required `MI`/`RX` tags
/// are missing on qualifying templates, or if tag values are invalid UTF-8.
pub fn process_templates_from_bam<F>(
    input: &Path,
    intervals: &[Interval],
    num_fractions: usize,
    mut process_group: F,
) -> Result<(usize, Vec<usize>)>
where
    F: FnMut(&[TemplateInfo], &mut Vec<usize>),
{
    let (reader, header) = create_raw_bam_reader(input, 1)?;

    let template_iter = TemplateIterator::new(reader);

    // Streaming approach: process groups as they arrive (assumes consecutive ReadInfo grouping)
    let mut current_group: Vec<TemplateInfo> = Vec::new();
    let mut current_key: Option<ReadInfoKey> = None;
    let mut template_count = 0;
    let progress = ProgressTracker::new("Processed records").with_interval(1_000_000);
    let mut fraction_template_counts: Vec<usize> = vec![0; num_fractions];
    let mi_tag = Tag::from(SamTag::MI);
    let umi_tag = Tag::from(SamTag::RX);

    for template in template_iter {
        let template = template?;
        if template.records().len() < 2 {
            continue;
        }

        // Decode raw records to RecordBuf for flag/tag/position access
        let record_bufs: Vec<RecordBuf> = template
            .records()
            .iter()
            .map(|r| raw_record_to_record_buf(r, &header))
            .collect::<anyhow::Result<Vec<_>>>()?;

        // Find R1 and R2 that pass fgbio's filtering criteria
        let r1 = record_bufs.iter().find(|r| {
            let f = r.flags();
            f.is_segmented()
                && !f.is_unmapped()
                && !f.is_mate_unmapped()
                && f.is_first_segment()
                && !f.is_secondary()
                && !f.is_supplementary()
        });
        let r2 = record_bufs.iter().find(|r| {
            let f = r.flags();
            f.is_segmented()
                && !f.is_unmapped()
                && !f.is_mate_unmapped()
                && f.is_last_segment()
                && !f.is_secondary()
                && !f.is_supplementary()
        });

        let (r1, r2) = match (r1, r2) {
            (Some(r1), Some(r2)) => (r1, r2),
            _ => continue,
        };

        // Get read name
        let read_name =
            r1.name().map(|n| String::from_utf8_lossy(n.as_ref()).to_string()).unwrap_or_default();

        // Get MI tag (molecular identifier)
        let mi = if let Some(noodles::sam::alignment::record_buf::data::field::Value::String(s)) =
            r1.data().get(&mi_tag)
        {
            String::from_utf8(s.iter().copied().collect::<Vec<u8>>())?
        } else {
            return Err(anyhow::anyhow!(
                "Read '{}' is missing the required MI tag. \
                 Metrics commands require standard MI/RX tags.",
                read_name
            ));
        };

        // Get UMI tag (raw UMI)
        let rx = if let Some(noodles::sam::alignment::record_buf::data::field::Value::String(s)) =
            r1.data().get(&umi_tag)
        {
            String::from_utf8(s.iter().copied().collect::<Vec<u8>>())?
        } else {
            return Err(anyhow::anyhow!(
                "Read '{}' is missing the required RX tag. \
                 Metrics commands require standard MI/RX tags.",
                read_name
            ));
        };

        // Get reference position and calculate insert coordinates
        let ref_name = if let Some(ref_id) = r1.reference_sequence_id() {
            header.reference_sequences().get_index(ref_id).map(|(name, _)| name.to_string())
        } else {
            None
        };

        // Calculate insert coordinates and ReadInfo key (matching fgbio)
        let (position, end_position, read_info_key) = {
            let r1_ref = r1.reference_sequence_id();
            let r2_ref = r2.reference_sequence_id();

            if r1_ref == r2_ref && r1_ref.is_some() {
                // Get unclipped 5' positions (matching fgbio's positionOf)
                let r1_5prime = unclipped_five_prime_position(r1);
                let r2_5prime = unclipped_five_prime_position(r2);
                let r1_strand = r1.flags().is_reverse_complemented();
                let r2_strand = r2.flags().is_reverse_complemented();

                if let (Some(s1), Some(s2)) = (r1_5prime, r2_5prime) {
                    // For insert coordinates, use alignment positions
                    let r1_start = r1.alignment_start().map(|p| usize::from(p) as i32);
                    let r2_start = r2.alignment_start().map(|p| usize::from(p) as i32);
                    let r1_end = r1.alignment_end().map(|p| usize::from(p) as i32);
                    let r2_end = r2.alignment_end().map(|p| usize::from(p) as i32);

                    let (pos, end) = match (r1_start, r2_start, r1_end, r2_end) {
                        (Some(rs1), Some(rs2), Some(re1), Some(re2)) => {
                            (rs1.min(rs2), re1.max(re2))
                        }
                        (Some(rs1), Some(rs2), _, _) => (rs1.min(rs2), rs1.max(rs2)),
                        _ => continue,
                    };

                    // Build ReadInfo key: order by (ref, 5' position) so earlier read is first
                    let key = if (r1_ref, s1) <= (r2_ref, s2) {
                        ReadInfoKey {
                            ref_index1: r1_ref,
                            start1: s1,
                            strand1: r1_strand,
                            ref_index2: r2_ref,
                            start2: s2,
                            strand2: r2_strand,
                        }
                    } else {
                        ReadInfoKey {
                            ref_index1: r2_ref,
                            start1: s2,
                            strand1: r2_strand,
                            ref_index2: r1_ref,
                            start2: s1,
                            strand2: r1_strand,
                        }
                    };

                    (pos, end, key)
                } else {
                    continue;
                }
            } else {
                continue;
            }
        };

        // Compute hash once for this template
        let hash_fraction = compute_hash_fraction(&read_name);

        let template_info = TemplateInfo {
            mi,
            rx,
            ref_name,
            position: Some(position),
            end_position: Some(end_position),
            hash_fraction,
        };

        // Check interval overlap
        if !overlaps_intervals(&template_info, intervals) {
            continue;
        }

        template_count += 1;
        progress.log_if_needed(2); // Each template has R1 and R2

        // Streaming: when ReadInfo key changes, process the accumulated group
        if current_key.as_ref() != Some(&read_info_key) && !current_group.is_empty() {
            process_group(&current_group, &mut fraction_template_counts);
            current_group.clear();
        }

        current_group.push(template_info);
        current_key = Some(read_info_key);
    }

    // Process the final group
    if !current_group.is_empty() {
        process_group(&current_group, &mut fraction_template_counts);
    }

    progress.log_final();
    Ok((template_count, fraction_template_counts))
}

//! Shared types and functions for metrics commands (duplex-metrics, simplex-metrics).
//!
//! This module contains the common infrastructure for reading UMI-grouped BAM files,
//! grouping templates by coordinate, filtering by genomic intervals, and performing
//! deterministic downsampling. Both `duplex_metrics` and `simplex_metrics` commands
//! build on these shared primitives.

use crate::bam_io::create_raw_bam_reader;
use crate::progress::ProgressTracker;
use crate::template::TemplateIterator;
use anyhow::{Context, Result};
use fgumi_raw_bam::{
    RawRecord, alignment_end_from_raw, aux_data_slice, find_string_tag_in_record, find_tag_type,
    flags as raw_flags, unclipped_5prime_from_raw_bam,
};

use log::info;
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
    pub ref_index1: usize,
    /// Unclipped 5' position for read 1.
    pub start1: i32,
    /// `true` if read 1 is reverse-complemented.
    pub strand1: bool,
    /// Reference sequence index for read 2.
    pub ref_index2: usize,
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
/// Delegates to [`unclipped_5prime_from_raw_bam`]. Includes both soft- and hard-clip
/// bases on the 5' side (matching htsjdk / fgbio semantics). Returns `None` for
/// unmapped records or records missing CIGAR ops.
fn unclipped_five_prime_position_raw(record: &RawRecord) -> Option<i32> {
    let flags = record.flags();
    if flags & raw_flags::UNMAPPED != 0 {
        return None;
    }
    if record.n_cigar_op() == 0 {
        return None;
    }
    Some(unclipped_5prime_from_raw_bam(record.as_ref()))
}

/// Computes an fgbio-compatible Murmur3 downsampling score.
///
/// Returns a value in `[0, 1]` for every hash except the Java `Int.MinValue`
/// overflow case, where fgbio/Scala's `math.abs` leaves `Int.MinValue`
/// unchanged and the quotient is slightly less than `-1`. Preserving this
/// quirk is required for byte-exact fgbio parity at every sampling fraction.
///
/// Mirrors fgbio's `CollectDuplexSeqMetrics` exactly:
///
/// ```scala
/// private val hasher = new htsjdk.samtools.util.Murmur3(42)
/// val intHash    = math.abs(hasher.hashUnencodedChars(rec.name))
/// val doubleHash = intHash / Int.MaxValue.toDouble
/// ```
///
/// The previous implementation used `murmur3::murmur3_32` over the UTF-8 bytes
/// of the read name. htsjdk's `hashUnencodedChars` walks the Java `char`
/// sequence (UTF-16 code units), so the two hashes diverge for every input and
/// produced a deterministic ~1% sampling bias vs. fgbio at every fraction.
///
/// For fgbio parity we port htsjdk's `Murmur3.hashUnencodedChars` byte-for-byte
/// and convert the read name to UTF-16 code units before hashing.
#[must_use]
pub fn compute_hash_fraction(read_name: &str) -> f64 {
    let chars: Vec<u16> = read_name.encode_utf16().collect();
    let hash = htsjdk_murmur3_hash_unencoded_chars(&chars, 42);
    // `wrapping_abs` mirrors Java `Math.abs` (which returns `Int.MinValue`
    // unchanged when the input is `Int.MinValue`) so the rare edge case
    // produces the same downsample decision as fgbio.
    f64::from(hash.wrapping_abs()) / f64::from(i32::MAX)
}

/// Port of htsjdk `Murmur3.hashUnencodedChars` (Apache-2.0; derived from
/// Guava's Apache-2.0 `Murmur3_32`; original `MurmurHash3` is public domain).
/// `chars` is the Java `CharSequence` / UTF-16 code units.
fn htsjdk_murmur3_hash_unencoded_chars(chars: &[u16], seed: i32) -> i32 {
    let mut h1: u32 = seed as u32;
    let length = chars.len();

    let mut i = 1;
    while i < length {
        let k1 = u32::from(chars[i - 1]) | (u32::from(chars[i]) << 16);
        h1 = murmur3_mix_h1(h1, murmur3_mix_k1(k1));
        i += 2;
    }

    if length & 1 == 1 {
        let k1 = murmur3_mix_k1(u32::from(chars[length - 1]));
        h1 ^= k1;
    }

    murmur3_fmix(h1, (2 * length) as u32) as i32
}

#[inline]
fn murmur3_mix_k1(mut k1: u32) -> u32 {
    k1 = k1.wrapping_mul(0xcc9e_2d51);
    k1 = k1.rotate_left(15);
    k1 = k1.wrapping_mul(0x1b87_3593);
    k1
}

#[inline]
fn murmur3_mix_h1(mut h1: u32, k1: u32) -> u32 {
    h1 ^= k1;
    h1 = h1.rotate_left(13);
    h1.wrapping_mul(5).wrapping_add(0xe654_6b64)
}

#[inline]
fn murmur3_fmix(mut h1: u32, length: u32) -> u32 {
    h1 ^= length;
    h1 ^= h1 >> 16;
    h1 = h1.wrapping_mul(0x85eb_ca6b);
    h1 ^= h1 >> 13;
    h1 = h1.wrapping_mul(0xc2b2_ae35);
    h1 ^= h1 >> 16;
    h1
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
    let (mut reader, _header) = create_raw_bam_reader(input, 1)?;

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
        let flags = raw.flags();
        if (flags & raw_flags::PAIRED) == 0
            || (flags & raw_flags::FIRST_SEGMENT) == 0
            || (flags & raw_flags::SECONDARY) != 0
            || (flags & raw_flags::SUPPLEMENTARY) != 0
        {
            continue;
        }

        // Consensus-tag check on raw aux bytes (mirrors
        // fgumi_consensus::tags::is_consensus: simplex = cD without aD+bD;
        // duplex = aD and bD). Avoids decoding the record to RecordBuf.
        let aux = aux_data_slice(raw.as_ref());
        let has_ad = find_tag_type(aux, b"aD").is_some();
        let has_bd = find_tag_type(aux, b"bD").is_some();
        let has_cd = find_tag_type(aux, b"cD").is_some();
        let is_duplex_consensus = has_ad && has_bd;
        let is_simplex_consensus = has_cd && !is_duplex_consensus;
        if is_simplex_consensus || is_duplex_consensus {
            let name = String::from_utf8_lossy(fgumi_raw_bam::read_name(raw.as_ref())).into_owned();
            anyhow::bail!(
                "Input BAM file ({}) appears to contain consensus sequences. \
                This metrics tool cannot run on consensus BAMs, and instead requires \
                the UMI-grouped BAM generated by group which is run prior to consensus calling.\n\
                First R1 record '{}' has consensus SAM tags present.",
                input.display(),
                name
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

    // fgbio R1/R2 filter: paired, both mapped, primary.
    let passes_filter = |r: &RawRecord, is_first: bool| -> bool {
        let f = r.flags();
        let seg_mask = if is_first { raw_flags::FIRST_SEGMENT } else { raw_flags::LAST_SEGMENT };
        (f & raw_flags::PAIRED) != 0
            && (f & raw_flags::UNMAPPED) == 0
            && (f & raw_flags::MATE_UNMAPPED) == 0
            && (f & seg_mask) != 0
            && (f & raw_flags::SECONDARY) == 0
            && (f & raw_flags::SUPPLEMENTARY) == 0
    };

    for template in template_iter {
        let template = template?;
        if template.records().len() < 2 {
            continue;
        }

        let r1 = template.records().iter().find(|r| passes_filter(r, true));
        let r2 = template.records().iter().find(|r| passes_filter(r, false));
        let (r1, r2) = match (r1, r2) {
            (Some(r1), Some(r2)) => (r1, r2),
            _ => continue,
        };

        let read_name = String::from_utf8_lossy(fgumi_raw_bam::read_name(r1.as_ref())).into_owned();
        let mi = required_z_tag(r1, *b"MI", &read_name)?;
        let rx = required_z_tag(r1, *b"RX", &read_name)?;

        // Filter already excluded unmapped reads, so tid >= 0 here; skip defensively.
        let r1_tid = r1.ref_id();
        let r2_tid = r2.ref_id();
        if r1_tid < 0 || r2_tid < 0 {
            continue;
        }
        let r1_ref = r1_tid as usize;
        let r2_ref = r2_tid as usize;
        let same_ref = r1_ref == r2_ref;

        // `ref_name` is always R1's reference. Interval overlap uses R1's own range
        // when R1 and R2 are on different chromosomes, matching fgbio
        // CollectDuplexSeqMetrics: `if (rec.refIndex == rec.mateRefIndex)
        // Bams.insertCoordinates(rec) else (rec.start, rec.end)`.
        let ref_name =
            header.reference_sequences().get_index(r1_ref).map(|(name, _)| name.to_string());

        // None here implies a malformed mapped record (no CIGAR); skip defensively.
        let (s1, s2) =
            match (unclipped_five_prime_position_raw(r1), unclipped_five_prime_position_raw(r2)) {
                (Some(s1), Some(s2)) => (s1, s2),
                _ => continue,
            };

        let r1_strand = (r1.flags() & raw_flags::REVERSE) != 0;
        let r2_strand = (r2.flags() & raw_flags::REVERSE) != 0;

        let r1_start = r1.pos() + 1;
        let r2_start = r2.pos() + 1;
        let r1_end = alignment_end_from_raw(r1.as_ref()).map(|e| e as i32);
        let r2_end = alignment_end_from_raw(r2.as_ref()).map(|e| e as i32);

        let (position, end_position) = if same_ref {
            match (r1_end, r2_end) {
                (Some(re1), Some(re2)) => (r1_start.min(r2_start), re1.max(re2)),
                _ => (r1_start.min(r2_start), r1_start.max(r2_start)),
            }
        } else {
            // No single insert interval spans both mates; use R1's own range so
            // interval filters still evaluate against R1's side of the pair.
            (r1_start, r1_end.unwrap_or(r1_start))
        };

        // ReadInfoKey fields are ordered so the earlier-mapping read comes first.
        let read_info_key = if (r1_ref, s1) <= (r2_ref, s2) {
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

        let hash_fraction = compute_hash_fraction(&read_name);

        let template_info = TemplateInfo {
            mi,
            rx,
            ref_name,
            position: Some(position),
            end_position: Some(end_position),
            hash_fraction,
        };

        if !overlaps_intervals(&template_info, intervals) {
            continue;
        }

        template_count += 1;
        progress.log_if_needed(2);

        // Flush the accumulated group when the ReadInfo key changes — input is
        // assumed to already be consecutively grouped by this key.
        if current_key.as_ref() != Some(&read_info_key) && !current_group.is_empty() {
            process_group(&current_group, &mut fraction_template_counts);
            current_group.clear();
        }

        current_group.push(template_info);
        current_key = Some(read_info_key);
    }

    if !current_group.is_empty() {
        process_group(&current_group, &mut fraction_template_counts);
    }

    progress.log_final();
    Ok((template_count, fraction_template_counts))
}

/// Extracts a required Z-typed aux tag from `record`, returning an error that
/// points at `read_name` when the tag is absent or not UTF-8.
fn required_z_tag(record: &RawRecord, tag: [u8; 2], read_name: &str) -> Result<String> {
    let tag_name = std::str::from_utf8(&tag).unwrap_or("??");
    let bytes = find_string_tag_in_record(record.as_ref(), &tag).ok_or_else(|| {
        anyhow::anyhow!(
            "Read '{read_name}' is missing the required {tag_name} tag. \
             Metrics commands require standard MI/RX tags."
        )
    })?;
    std::str::from_utf8(bytes)
        .map(str::to_string)
        .map_err(|e| anyhow::anyhow!("Read '{read_name}' {tag_name} tag is not UTF-8: {e}"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use fgumi_raw_bam::{SamBuilder as RawSamBuilder, flags as raw_flags, testutil::encode_op};
    use noodles::bam;
    use noodles::sam;
    use noodles::sam::alignment::io::Write as AlignmentWrite;
    use noodles::sam::alignment::record_buf::RecordBuf;
    use std::num::NonZeroUsize;
    use tempfile::NamedTempFile;

    fn test_header() -> sam::Header {
        use noodles::sam::header::record::value::Map;
        use noodles::sam::header::record::value::map::ReferenceSequence;
        sam::Header::builder()
            .add_reference_sequence(
                bstr::BString::from("chr1"),
                Map::<ReferenceSequence>::new(NonZeroUsize::new(248_956_422).expect("non-zero")),
            )
            .add_reference_sequence(
                bstr::BString::from("chr2"),
                Map::<ReferenceSequence>::new(NonZeroUsize::new(242_193_529).expect("non-zero")),
            )
            .build()
    }

    /// Build an R1/R2 pair with independent refs/positions for each mate.
    fn build_pair(
        name: &str,
        r1_ref: i32,
        r1_pos: i32,
        r2_ref: i32,
        r2_pos: i32,
        mi: &str,
    ) -> (RecordBuf, RecordBuf) {
        let seq = vec![b'A'; 100];
        let quals = vec![30u8; 100];
        let cigar = encode_op(0, 100); // 100M

        let mut b1 = RawSamBuilder::new();
        b1.read_name(name.as_bytes())
            .flags(raw_flags::PAIRED | raw_flags::FIRST_SEGMENT | raw_flags::MATE_REVERSE)
            .ref_id(r1_ref)
            .pos(r1_pos - 1)
            .mapq(60)
            .cigar_ops(&[cigar])
            .sequence(&seq)
            .qualities(&quals)
            .mate_ref_id(r2_ref)
            .mate_pos(r2_pos - 1);
        b1.add_string_tag(b"RX", b"ACGT-TGCA");
        b1.add_string_tag(b"MI", mi.as_bytes());
        let r1 = fgumi_raw_bam::raw_record_to_record_buf(&b1.build(), &sam::Header::default())
            .expect("decode r1");

        let mut b2 = RawSamBuilder::new();
        b2.read_name(name.as_bytes())
            .flags(raw_flags::PAIRED | raw_flags::LAST_SEGMENT | raw_flags::REVERSE)
            .ref_id(r2_ref)
            .pos(r2_pos - 1)
            .mapq(60)
            .cigar_ops(&[cigar])
            .sequence(&seq)
            .qualities(&quals)
            .mate_ref_id(r1_ref)
            .mate_pos(r1_pos - 1);
        b2.add_string_tag(b"RX", b"ACGT-TGCA");
        b2.add_string_tag(b"MI", mi.as_bytes());
        let r2 = fgumi_raw_bam::raw_record_to_record_buf(&b2.build(), &sam::Header::default())
            .expect("decode r2");

        (r1, r2)
    }

    fn write_test_bam(records: Vec<RecordBuf>) -> NamedTempFile {
        let file = NamedTempFile::new().expect("tempfile");
        let header = test_header();
        let mut writer =
            bam::io::writer::Builder.build_from_path(file.path()).expect("open writer");
        writer.write_header(&header).expect("write header");
        for r in &records {
            writer.write_alignment_record(&header, r).expect("write record");
        }
        drop(writer);
        file
    }

    use rstest::rstest;

    /// Reference values captured directly from htsjdk `Murmur3(42)
    /// .hashUnencodedChars(s)` against the 3.1.2 `Murmur3.class` on a set of
    /// read-name-shaped strings.  If this test fails, the Rust port has
    /// diverged from htsjdk; all fgbio-parity guarantees for downsampling are
    /// invalid until the port is corrected.
    #[rstest]
    #[case("", 142_593_372)]
    #[case("A", 309_601_938)]
    #[case("AB", 1_297_118_606)]
    #[case("ABC", 417_488_640)]
    #[case("read1", -958_943_510)]
    #[case("read2", 1_466_959_157)]
    #[case("read10", -87_319_652)]
    #[case("SRR099966.100", -1_840_920_289)]
    #[case("M00517:73:000000000-A5AEH:1:1101:15541:1541", 1_482_717_766)]
    #[case("NB500947:HT3JMBGX2:1:11101:19204:10048", -1_636_484_024)]
    fn test_murmur3_matches_htsjdk_reference_vectors(#[case] name: &str, #[case] expected: i32) {
        let chars: Vec<u16> = name.encode_utf16().collect();
        let got = htsjdk_murmur3_hash_unencoded_chars(&chars, 42);
        assert_eq!(got, expected, "Murmur3 mismatch on {name:?}");
    }

    /// The downsample fraction must be in `[0, 1]` for all non-`i32::MIN`
    /// hashes, matching fgbio's `math.abs(hash) / Int.MaxValue.toDouble`.
    #[rstest]
    #[case("")]
    #[case("A")]
    #[case("read1")]
    #[case("SRR099966.100")]
    #[case("a much longer read name here")]
    fn test_compute_hash_fraction_in_unit_range(#[case] name: &str) {
        let f = compute_hash_fraction(name);
        // Abs can produce up to Int.MaxValue, divided by itself == 1.0.
        assert!((0.0..=1.0).contains(&f), "compute_hash_fraction({name:?}) = {f}");
    }

    /// Regression test for fgbio parity: pairs whose mates map to different
    /// chromosomes must be kept (not silently dropped), matching fgbio's
    /// `CollectDuplexSeqMetrics` which retains inter-reference pairs and uses
    /// R1's own range for interval-overlap evaluation.
    #[test]
    fn test_inter_reference_pairs_are_retained() {
        // Two same-ref pairs (chr1:100 / chr1:100) + one inter-ref pair
        // (chr1:500 / chr2:500).  All three should be counted.
        let (s1r1, s1r2) = build_pair("same_1", 0, 100, 0, 300, "1");
        let (s2r1, s2r2) = build_pair("same_2", 0, 100, 0, 300, "2");
        let (ir1, ir2) = build_pair("inter_1", 0, 500, 1, 500, "3");
        let bam = write_test_bam(vec![s1r1, s1r2, s2r1, s2r2, ir1, ir2]);

        let mut groups: Vec<Vec<String>> = Vec::new();
        let (total, _) = process_templates_from_bam(bam.path(), &[], 1, |group, _| {
            groups.push(group.iter().map(|t| t.mi.clone()).collect());
        })
        .expect("process_templates_from_bam");

        assert_eq!(total, 3, "inter-reference pair must not be dropped");
        let mis: Vec<String> = groups.into_iter().flatten().collect();
        assert!(mis.contains(&"1".to_string()));
        assert!(mis.contains(&"2".to_string()));
        assert!(mis.contains(&"3".to_string()), "inter-ref pair's MI must be in output");
    }
}

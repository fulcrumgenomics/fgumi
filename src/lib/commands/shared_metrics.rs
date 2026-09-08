//! Shared types and functions for metrics commands (duplex-metrics, simplex-metrics).
//!
//! This module contains the common infrastructure for reading UMI-grouped BAM files,
//! grouping templates by coordinate, filtering by genomic intervals, and performing
//! deterministic downsampling. Both `duplex_metrics` and `simplex_metrics` commands
//! build on these shared primitives.

use crate::metrics::duplex::DuplexMetricsCollector;
use crate::metrics::simplex::SimplexMetricsCollector;
use crate::read_info::LibraryIndex;
use crate::sam::SamTag;
use crate::simple_umi_consensus::SimpleUmiConsensusCaller;
use crate::template::TemplateIterator;
use crate::umi::extract_mi_base;
use anyhow::{Context, Result};
use fgumi_bam_io::ProgressTracker;
use fgumi_bam_io::create_raw_bam_reader;
use fgumi_raw_bam::{
    AsTagBytes, RawRecord, alignment_end_from_raw, aux_data_slice, find_string_tag_in_record,
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
    /// `true` if read 1 (the first segment) is on the positive strand. Used to
    /// orient duplex UMIs to the F1R2 reading of the top strand (DXM-02).
    pub r1_positive: bool,
    /// Hash fraction for deterministic downsampling (computed once per template).
    pub hash_fraction: f64,
}

/// Grouping key matching fgbio's `ReadInfo` structure.
///
/// The two mate positions are ordered so the earlier-mapping read comes first;
/// `library` and `cell_barcode` are template-level (identical for both mates) and are
/// part of the key so that reads from different libraries or cells at the same
/// coordinate/strand form separate families, matching fgbio's `ReadInfo` (which carries
/// `library` and `cellBarcode`) and fgumi's own `group`/`dedup` grouping.
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
    /// Library index (from the read group's `LB`, via [`LibraryIndex`]); 0 = unknown.
    pub library: u16,
    /// Cell barcode (`CB` tag) value, or `None` when the tag is absent.
    pub cell_barcode: Option<Box<[u8]>>,
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

/// Computes the unclipped 5' position for a read.
///
/// Delegates to [`unclipped_5prime_from_raw_bam`], which subtracts both soft- **and**
/// hard-clip bases on the 5' side. This deliberately matches samtools
/// (`unclipped_start` in `bam.c`, used by `markdup` and template-coordinate sort) and
/// htsjdk `getUnclippedStart`, keeping this grouping key consistent with how fgumi's
/// own `group`/`dedup` compute coordinates. It intentionally diverges from fgbio's
/// `positionOf`/`unSoftClippedStart`, which counts soft clips only. Returns `None` for
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
    // SAM restricts read names to printable ASCII, and for those the UTF-16 code units htsjdk
    // hashes are exactly the bytes widened — so the byte entry point is identical here and
    // avoids a per-read `Vec<u16>`. A non-ASCII name violates the spec but must still hash the
    // way htsjdk would, so it takes the widening path.
    let hash = if read_name.is_ascii() {
        fgumi_raw_bam::hash::fgbio_read_name_rank(read_name.as_bytes())
    } else {
        let chars: Vec<u16> = read_name.encode_utf16().collect();
        fgumi_raw_bam::hash::htsjdk_murmur3_hash_unencoded_chars(&chars, 42)
    };
    // `wrapping_abs` mirrors Java `Math.abs` (which returns `Int.MinValue`
    // unchanged when the input is `Int.MinValue`) so the rare edge case
    // produces the same downsample decision as fgbio.
    f64::from(hash.wrapping_abs()) / f64::from(i32::MAX)
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

/// Whether `raw` is a primary record — neither secondary nor supplementary.
///
/// Deliberately does **not** exclude `UNMAPPED`: consensus BAMs are documented as
/// unaligned, so skipping unmapped records here would let exactly the input the
/// consensus guard exists to reject slip past unchecked. That makes this
/// predicate looser than [`process_templates_from_bam`]'s own R1/R2 filter, which
/// does require both mates mapped.
fn is_primary_record(raw: &RawRecord) -> bool {
    let flags = raw.flags();
    (flags & raw_flags::SECONDARY) == 0 && (flags & raw_flags::SUPPLEMENTARY) == 0
}

/// Whether `raw` is the record the consensus-BAM guard prefers to inspect.
///
/// A paired, primary R1 — the record fgbio's own check looks at, and the one the
/// metrics pass itself goes on to read tags from. See [`is_primary_record`] for
/// why `UNMAPPED` is not excluded.
pub(crate) fn is_consensus_guard_record(raw: &RawRecord) -> bool {
    let flags = raw.flags();
    (flags & raw_flags::PAIRED) != 0
        && (flags & raw_flags::FIRST_SEGMENT) != 0
        && is_primary_record(raw)
}

/// The record in `records` that the consensus-BAM guard inspects, if any.
///
/// A paired primary R1 when the template has one ([`is_consensus_guard_record`]);
/// otherwise the first primary record. The fallback is what makes the guard
/// unconditional: a fragment or single-end consensus BAM has no paired R1 at all,
/// so requiring one meant the check never ran, and the very input it exists to
/// reject fell through the `len() < 2` skip in
/// [`process_templates_from_bam`] to an empty metrics file rather than the
/// explicit error. Consensus tags sit on every consensus record, so any primary
/// record answers the question just as well.
///
/// Returns `None` only for a template holding nothing but secondary and
/// supplementary records, which leaves the caller to try the next template.
pub(crate) fn consensus_guard_record(records: &[RawRecord]) -> Option<&RawRecord> {
    records
        .iter()
        .find(|raw| is_consensus_guard_record(raw))
        .or_else(|| records.iter().find(|raw| is_primary_record(raw)))
}

/// Errors if `raw` carries consensus-calling tags.
///
/// Consensus BAMs (output from the simplex/duplex callers) must not be fed to
/// the metrics tools, which expect the UMI-grouped BAM produced by `group`.
///
/// The check reads raw aux bytes rather than decoding to a `RecordBuf`, via
/// [`fgumi_consensus::is_consensus`].
///
/// # Errors
///
/// Returns an error if `raw` appears to be a consensus read.
pub(crate) fn ensure_not_consensus_record(raw: &RawRecord, input: &Path) -> Result<()> {
    // The predicate lives in `fgumi-consensus` next to the record-level form, so
    // "what counts as a consensus read" has one definition rather than a copy
    // per consumer.
    if fgumi_consensus::is_consensus(aux_data_slice(raw.as_ref())) {
        let name = String::from_utf8_lossy(fgumi_raw_bam::read_name(raw.as_ref())).into_owned();
        anyhow::bail!(
            "Input BAM file ({}) appears to contain consensus sequences. \
            This metrics tool cannot run on consensus BAMs, and instead requires \
            the UMI-grouped BAM generated by group which is run prior to consensus calling.\n\
            Record '{}' has consensus SAM tags present.",
            input.display(),
            name
        );
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

/// Records one coordinate/strand group's family-size and UMI-consensus
/// contributions across every downsampling fraction. Shared by the
/// separate-pass `simplex-metrics` command and the inline consensus-metrics
/// accumulator — extracted verbatim from `SimplexMetrics::process_coordinate_group`.
///
/// `group` must already be interval-filtered by the caller.
pub(crate) fn record_simplex_coordinate_group(
    group: &[TemplateInfo],
    fractions: &[f64],
    collectors: &mut [SimplexMetricsCollector],
    umi_consensus_caller: &mut SimpleUmiConsensusCaller,
    fraction_template_counts: &mut [usize],
) -> Result<()> {
    use std::collections::HashMap;

    if group.is_empty() {
        return Ok(());
    }

    let metadata = compute_template_metadata(group);

    let mut base_umi_strands: HashMap<&str, (bool, bool)> = HashMap::new();
    for m in &metadata {
        let seen = base_umi_strands.entry(m.base_umi).or_default();
        seen.0 |= m.is_a_strand;
        seen.1 |= m.is_b_strand;
        if seen.0 && seen.1 {
            anyhow::bail!(
                "simplex-metrics received duplex-UMI data: base UMI '{}' has reads on \
                 both the /A and /B strands. Run duplex-metrics for duplex data.",
                m.base_umi
            );
        }
    }

    let last_fraction_idx = fractions.len() - 1;
    let mut ss_groups: HashMap<&str, usize> = HashMap::new();

    for (idx, &fraction) in fractions.iter().enumerate() {
        let downsampled: Vec<_> =
            metadata.iter().filter(|m| m.template.hash_fraction <= fraction).collect();

        if downsampled.is_empty() {
            continue;
        }

        fraction_template_counts[idx] += downsampled.len();
        collectors[idx].record_cs_family(downsampled.len());

        ss_groups.clear();
        for m in &downsampled {
            *ss_groups.entry(m.template.mi.as_str()).or_default() += 1;
        }
        for &ss_size in ss_groups.values() {
            collectors[idx].record_ss_family(ss_size);
        }

        if idx == last_fraction_idx {
            let mut umi_groups: HashMap<&str, Vec<&str>> = HashMap::new();
            for m in &downsampled {
                umi_groups.entry(m.base_umi).or_default().push(m.template.rx.as_str());
            }

            for rx_tags in umi_groups.values() {
                let split_rx: Vec<Vec<&str>> =
                    rx_tags.iter().map(|rx| rx.split('-').collect()).collect();
                // Use the maximum component count across the family, not the first
                // value's. `--no-umi` assigns one MI to every template while keeping
                // their raw `RX` tags, so a base-UMI family can hold ragged values
                // (e.g. "A-C" alongside "A-C-G"). Keying off the first value would
                // drop trailing components from longer values; `parts.get(pos)` below
                // still skips components absent from shorter values. Simplex
                // intentionally supports N-component UMIs, so no segment-count bail.
                let num_components = split_rx.iter().map(Vec::len).max().unwrap_or(0);

                for pos in 0..num_components {
                    let umis_at_pos: Vec<String> = split_rx
                        .iter()
                        .filter_map(|parts| parts.get(pos).map(|s| (*s).to_string()))
                        .collect();

                    if umis_at_pos.is_empty() {
                        continue;
                    }

                    let (consensus, _had_errors) = umi_consensus_caller.consensus(&umis_at_pos);
                    let raw_count = umis_at_pos.len();
                    let error_count = umis_at_pos.iter().filter(|u| **u != consensus).count();
                    collectors[idx].record_umi(&consensus, raw_count, error_count, true);
                }
            }
        }
    }
    Ok(())
}

/// Records one coordinate/strand group's family-size, duplex-family, and
/// UMI-consensus contributions across every downsampling fraction. Shared by
/// the separate-pass `duplex-metrics` command and the inline consensus-metrics
/// accumulator — extracted verbatim from `DuplexMetrics::process_coordinate_group`.
///
/// `group` must already be interval-filtered by the caller. `duplex_umi_counts`
/// replaces the `&DuplexMetrics` receiver the original method read
/// `self.duplex_umi_counts` from (the one field it needed) — the inline path
/// has no `DuplexMetrics` command-struct instance to supply.
pub(crate) fn record_duplex_coordinate_group(
    group: &[TemplateInfo],
    fractions: &[f64],
    collectors: &mut [DuplexMetricsCollector],
    umi_consensus_caller: &mut SimpleUmiConsensusCaller,
    fraction_template_counts: &mut [usize],
    duplex_umi_counts: bool,
) -> Result<()> {
    use std::collections::HashMap;

    if group.is_empty() {
        return Ok(());
    }

    // Pre-compute metadata once for the entire group
    let metadata = compute_template_metadata(group);

    // Hoist scratch buffers outside the 20-fraction loop. `HashMap::clear()`
    // preserves the outer bucket array across iterations — that's the
    // dominant allocator win on real cfDNA inputs with millions of
    // coordinate groups, and is the same pattern
    // `simplex_metrics::process_coordinate_group` already uses for its
    // `ss_groups` map.
    //
    // The inner `Vec<(&str, &str, bool)>` stored as `ds_groups` entry.2 is
    // still freed per entry on `.clear()`. Combined with the
    // `is_full_fraction` gate below, that inner Vec now allocates exactly
    // once per coordinate group (at the 100% fraction) rather than once
    // per fraction × per `or_default()` slot. The `bool` is `r1_positive`
    // (read 1 on the positive strand), used to orient the duplex UMI.
    let mut downsampled: Vec<&TemplateMetadata> = Vec::new();
    let mut ss_groups: HashMap<&str, usize> = HashMap::new();
    #[allow(clippy::type_complexity)]
    let mut ds_groups: HashMap<&str, (usize, usize, Vec<(&str, &str, bool)>)> = HashMap::new();

    // For each fraction: filter ONCE, then groupBy (like fgbio)
    for (idx, &fraction) in fractions.iter().enumerate() {
        // Filter once per fraction - equivalent to fgbio's downsampledGroup
        downsampled.clear();
        downsampled.extend(metadata.iter().filter(|m| m.template.hash_fraction <= fraction));

        if downsampled.is_empty() {
            continue;
        }

        // CS family size
        fraction_template_counts[idx] += downsampled.len();
        collectors[idx].record_cs_family(downsampled.len());

        let is_full_fraction = (fraction - 1.0_f64).abs() < 0.01;

        // Group by MI tag for SS families (like fgbio's groupBy)
        ss_groups.clear();
        for m in &downsampled {
            *ss_groups.entry(m.template.mi.as_str()).or_default() += 1;
        }
        for &ss_size in ss_groups.values() {
            collectors[idx].record_ss_family(ss_size);
        }

        // Group by base_umi for DS families with strand counts.
        // HashMap value: (a_count, b_count, mi_rx_pairs for UMI metrics).
        // The mi/rx pair vec is consumed only at the 100% fraction by
        // `record_duplex_umi_metrics`; skip populating it on downsampled fractions
        // to avoid O(group_size) wasted pushes per non-full fraction.
        ds_groups.clear();
        for m in &downsampled {
            let entry = ds_groups.entry(m.base_umi).or_default();
            if m.is_b_strand {
                entry.1 += 1;
            } else {
                // /A or an unsuffixed MI counts toward the AB strand. fgbio treats a
                // single unsuffixed MI as Pair(ab=n, ba=0) — a valid single-strand DS
                // family (CollectDuplexSeqMetrics). Without this, unsuffixed families
                // get ds_size=0 and are silently dropped by the family-size histogram
                // (which iterates 1..=max), undercounting ds_families (DXM3-03).
                entry.0 += 1;
            }
            if is_full_fraction {
                entry.2.push((
                    m.template.mi.as_str(),
                    m.template.rx.as_str(),
                    m.template.r1_positive,
                ));
            }
        }

        for (base_umi, (a_count, b_count, mi_rx_pairs)) in &ds_groups {
            let ds_size = a_count + b_count;
            collectors[idx].record_ds_family(ds_size);

            let (ab_count, ba_count) =
                if a_count >= b_count { (*a_count, *b_count) } else { (*b_count, *a_count) };

            collectors[idx].record_duplex_family(ab_count, ba_count);

            // Only collect UMI metrics for the 100% fraction. Pass the
            // (&str, &str) pairs directly — the underlying String storage
            // lives in `group`, which outlives this call, so the previous
            // `Vec<(String, String)>` clone was pure overhead.
            if is_full_fraction {
                record_duplex_umi_metrics(
                    &mut collectors[idx],
                    umi_consensus_caller,
                    mi_rx_pairs.as_slice(),
                    base_umi,
                    duplex_umi_counts,
                )?;
            }
        }
    }
    Ok(())
}

/// Updates UMI metrics for a duplex family
///
/// This method:
/// 1. Uses RX tags (raw UMI sequences) to extract individual UMI observations
/// 2. Separates by strand (/A and /B suffixes in MI tags), swapping UMI parts for B strand
/// 3. Calls consensus for each UMI position
/// 4. Records raw observations, errors, and unique observations for each individual UMI
/// 5. Records duplex UMI metrics if enabled
///
/// This matches the Scala implementation in CollectDuplexSeqMetrics.scala:407-431
///
/// Extracted from `DuplexMetrics::update_umi_metrics`; the `&self` receiver is
/// replaced by an explicit `duplex_umi_counts` parameter (the one
/// `self.duplex_umi_counts` field the method needed) so the inline path can call
/// this without a `DuplexMetrics` command-struct instance.
pub(crate) fn record_duplex_umi_metrics(
    collector: &mut DuplexMetricsCollector,
    umi_consensus_caller: &mut SimpleUmiConsensusCaller,
    group_pairs: &[(&str, &str, bool)],
    base_umi: &str,
    duplex_umi_counts: bool,
) -> Result<()> {
    // Collect the two UMI positions, each oriented to the F1R2 reading of the
    // top strand. umi1s holds the leading half, umi2s the trailing half.
    let mut umi1s = Vec::new();
    let mut umi2s = Vec::new();

    for &(mi, rx, r1_positive) in group_pairs {
        // Check if this MI tag belongs to the current base_umi family
        let mi_base = extract_mi_base(mi);

        if mi_base != base_umi {
            continue;
        }

        // Split the RX tag to get individual UMI parts. fgbio uses
        // `split("-", -1)`, which keeps a trailing empty field, and requires
        // exactly two parts (`case Array(u1, u2)`), throwing a `MatchError` on
        // anything else. Reject a malformed duplex UMI here with a clear error
        // rather than silently skipping it — matching fgbio's fail-fast and
        // fgumi's own `group`/`dedup`, which bail on non-2-segment paired UMIs
        // (DXM-03). Empty molecule-end halves (`-CCC`, `CCC-`) are still two
        // parts and are kept (DXM-01).
        let parts: Vec<&str> = rx.split('-').collect();
        if parts.len() != 2 {
            anyhow::bail!(
                "Duplex UMI did not contain 2 segments delimited by '-': '{rx}' (MI '{mi}')"
            );
        }

        // Do NOT skip empty molecule-end halves (e.g. `-CCC` or `CCC-`). fgbio
        // counts them, recording the empty half as an empty-string UMI, so
        // single-index / single-strand designs are not undercounted (DXM-01).

        // Orient by the actual R1 strand, not the MI `/A`,`/B` suffix (which is
        // not strand-reliable — e.g. an `/A` family can be on the negative
        // strand). If R1 is on the positive strand the molecule was read
        // top-strand-first, so the RX is already `u1-u2`; otherwise it was read
        // bottom-strand-first, so swap the halves. Because metrics collection is
        // R1-only, this reproduces fgbio's per-read F1R2 normalization
        // (CollectDuplexSeqMetrics.scala:407-408) and its duplex-UMI orientation
        // pick (:419-425), which then reduces to "lead with the positive-strand
        // read's half" (DXM-02).
        if r1_positive {
            umi1s.push(parts[0].to_string());
            umi2s.push(parts[1].to_string());
        } else {
            umi1s.push(parts[1].to_string());
            umi2s.push(parts[0].to_string());
        }
    }

    // Call consensus for each UMI position and record metrics
    let mut consensus_umis = Vec::new();

    if !umi1s.is_empty() {
        let (consensus, _had_errors) = umi_consensus_caller.consensus(&umi1s);
        let raw_count = umi1s.len();
        let error_count = umi1s.iter().filter(|u| **u != consensus).count();
        collector.record_umi(&consensus, raw_count, error_count, true);
        consensus_umis.push(consensus);
    }

    if !umi2s.is_empty() {
        let (consensus, _had_errors) = umi_consensus_caller.consensus(&umi2s);
        let raw_count = umi2s.len();
        let error_count = umi2s.iter().filter(|u| **u != consensus).count();
        collector.record_umi(&consensus, raw_count, error_count, true);
        consensus_umis.push(consensus);
    }

    // Record duplex UMI metrics if enabled
    if duplex_umi_counts && consensus_umis.len() == 2 {
        let duplex_umi = format!("{}-{}", consensus_umis[0], consensus_umis[1]);
        // Each read pair contributes one observation to the duplex UMI
        // (not two, even though we track each component separately)
        let total_raw = umi1s.len();

        // Count how many raw RX tags had errors (don't match either duplex orientation)
        let expected_duplex1 = format!("{}-{}", consensus_umis[0], consensus_umis[1]);
        let expected_duplex2 = format!("{}-{}", consensus_umis[1], consensus_umis[0]);
        let error_count = group_pairs
            .iter()
            .filter(|&&(mi, rx, _r1_positive)| {
                let mi_base = extract_mi_base(mi);
                mi_base == base_umi
                    && rx != expected_duplex1.as_str()
                    && rx != expected_duplex2.as_str()
            })
            .count();

        collector.record_duplex_umi(&duplex_umi, total_raw, error_count, true);
    }

    Ok(())
}

/// Builds [`TemplateInfo`] and [`ReadInfoKey`] from one already-selected R1/R2
/// raw-record pair. Extracted from `process_templates_from_bam`'s per-template
/// loop body so both the separate-pass metrics commands and the inline
/// consensus-metrics adapters (which source R1/R2 pairs differently — one
/// from a `Template`'s cached views, one from a re-paired `MiGroup`) share
/// one construction path and cannot numerically drift apart.
///
/// Returns `Ok(None)` for the same two defensive-skip cases the original
/// loop had: an unmapped reference id on either mate, or a missing CIGAR (no
/// unclipped 5' position). Returns `Err` only if a required `MI`/`RX` tag is
/// absent from `r1`.
pub(crate) fn build_template_info(
    r1: &RawRecord,
    r2: &RawRecord,
    header: &noodles::sam::Header,
    library_index: &LibraryIndex,
) -> Result<Option<(TemplateInfo, ReadInfoKey)>> {
    let read_name = String::from_utf8_lossy(fgumi_raw_bam::read_name(r1.as_ref())).into_owned();
    let mi = required_z_tag(r1, SamTag::MI, &read_name)?;
    let rx = required_z_tag(r1, SamTag::RX, &read_name)?;

    let r1_tid = r1.ref_id();
    let r2_tid = r2.ref_id();
    if r1_tid < 0 || r2_tid < 0 {
        return Ok(None);
    }
    let r1_ref = r1_tid as usize;
    let r2_ref = r2_tid as usize;
    let same_ref = r1_ref == r2_ref;

    let ref_name = header.reference_sequences().get_index(r1_ref).map(|(name, _)| name.to_string());

    let (s1, s2) =
        match (unclipped_five_prime_position_raw(r1), unclipped_five_prime_position_raw(r2)) {
            (Some(s1), Some(s2)) => (s1, s2),
            _ => return Ok(None),
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
        (r1_start, r1_end.unwrap_or(r1_start))
    };

    let library = find_string_tag_in_record(r1.as_ref(), SamTag::RG)
        .map_or(0u16, |rg| library_index.get(LibraryIndex::hash_rg(rg)));
    let cell_barcode: Option<Box<[u8]>> =
        find_string_tag_in_record(r1.as_ref(), SamTag::CB).map(Box::from);

    let (ref_index1, start1, strand1, ref_index2, start2, strand2) =
        if (r1_ref, s1, r1_strand) <= (r2_ref, s2, r2_strand) {
            (r1_ref, s1, r1_strand, r2_ref, s2, r2_strand)
        } else {
            (r2_ref, s2, r2_strand, r1_ref, s1, r1_strand)
        };
    let read_info_key = ReadInfoKey {
        ref_index1,
        start1,
        strand1,
        ref_index2,
        start2,
        strand2,
        library,
        cell_barcode,
    };

    let hash_fraction = compute_hash_fraction(&read_name);

    let template_info = TemplateInfo {
        mi,
        rx,
        ref_name,
        position: Some(position),
        end_position: Some(end_position),
        r1_positive: !r1_strand,
        hash_fraction,
    };

    Ok(Some((template_info, read_info_key)))
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
    F: FnMut(&[TemplateInfo], &mut Vec<usize>) -> Result<()>,
{
    let (reader, header) = create_raw_bam_reader(input, 1)?;

    // Library index (RG -> LB) for partitioning families by library, matching fgbio's
    // ReadInfo.library and fgumi's own group/dedup grouping.
    let library_index = LibraryIndex::from_header(&header);

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

    // Reject a consensus BAM against the first qualifying record of this pass,
    // rather than by re-opening the input beforehand. Folding it in is what lets
    // the metrics tools read from a pipe: they no longer need the input twice.
    let mut consensus_checked = false;

    for template in template_iter {
        let template = template?;

        if !consensus_checked && let Some(guard_record) = consensus_guard_record(template.records())
        {
            ensure_not_consensus_record(guard_record, input)?;
            consensus_checked = true;
        }

        if template.records().len() < 2 {
            continue;
        }

        let r1 = template.records().iter().find(|r| passes_filter(r, true));
        let r2 = template.records().iter().find(|r| passes_filter(r, false));
        let (r1, r2) = match (r1, r2) {
            (Some(r1), Some(r2)) => (r1, r2),
            _ => continue,
        };

        let Some((template_info, read_info_key)) =
            build_template_info(r1, r2, &header, &library_index)?
        else {
            continue;
        };

        if !overlaps_intervals(&template_info, intervals) {
            continue;
        }

        template_count += 1;
        progress.log_if_needed(2);

        // Flush the accumulated group when the ReadInfo key changes — input is
        // assumed to already be consecutively grouped by this key.
        if current_key.as_ref() != Some(&read_info_key) && !current_group.is_empty() {
            process_group(&current_group, &mut fraction_template_counts)?;
            current_group.clear();
        }

        current_group.push(template_info);
        current_key = Some(read_info_key);
    }

    if !current_group.is_empty() {
        process_group(&current_group, &mut fraction_template_counts)?;
    }

    progress.log_final();
    Ok((template_count, fraction_template_counts))
}

/// Extracts a required Z-typed aux tag from `record`, returning an error that
/// points at `read_name` when the tag is absent or not UTF-8.
fn required_z_tag(record: &RawRecord, tag: impl AsTagBytes, read_name: &str) -> Result<String> {
    let tag_bytes = *tag.as_tag_bytes();
    let tag_name = std::str::from_utf8(&tag_bytes).unwrap_or("??");
    let bytes = find_string_tag_in_record(record.as_ref(), tag).ok_or_else(|| {
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
pub(crate) mod tests {
    use super::*;
    use fgumi_raw_bam::{SamBuilder as RawSamBuilder, flags as raw_flags, testutil::encode_op};
    use noodles::bam;
    use noodles::sam;
    use noodles::sam::alignment::io::Write as AlignmentWrite;
    use noodles::sam::alignment::record_buf::RecordBuf;
    use std::num::NonZeroUsize;
    use tempfile::NamedTempFile;

    pub(crate) fn test_header() -> sam::Header {
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
    pub(crate) fn build_pair(
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
        b1.add_string_tag(SamTag::RX, b"ACGT-TGCA");
        b1.add_string_tag(SamTag::MI, mi.as_bytes());
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
        b2.add_string_tag(SamTag::RX, b"ACGT-TGCA");
        b2.add_string_tag(SamTag::MI, mi.as_bytes());
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
        let got = fgumi_raw_bam::hash::htsjdk_murmur3_hash_unencoded_chars(&chars, 42);
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
            Ok(())
        })
        .expect("process_templates_from_bam");

        assert_eq!(total, 3, "inter-reference pair must not be dropped");
        let mis: Vec<String> = groups.into_iter().flatten().collect();
        assert!(mis.contains(&"1".to_string()));
        assert!(mis.contains(&"2".to_string()));
        assert!(mis.contains(&"3".to_string()), "inter-ref pair's MI must be in output");
    }

    #[test]
    fn build_template_info_populates_all_template_info_and_read_info_key_fields() {
        // build_pair (this module's own helper, above) builds an R1/R2 pair via
        // fgumi_raw_bam::RawSamBuilder. Both mates: 100M, no soft-clips, no RG/CB
        // tags. R1 is FIRST_SEGMENT + MATE_REVERSE (forward strand) at 1-based pos
        // 100; R2 is LAST_SEGMENT + REVERSE (reverse strand) at 1-based pos 150.
        // Both are tagged MI="7", RX="ACGT-TGCA".
        //
        // Every expected value below is derived directly from those known inputs
        // (not by re-running build_template_info against process_templates_from_bam),
        // so a regression in the builder itself is caught rather than mirrored.
        let (r1_buf, r2_buf) = build_pair("read-1", 0, 100, 0, 150, "7");
        // `encode_record_buf_to_raw` needs a non-empty `@SQ` dictionary to round-trip a
        // mapped record's reference id (see its docstring); `sam::Header::default()` has
        // none, so use this module's `test_header()` (chr1 at index 0), matching build_pair's
        // ref_id 0 for both mates.
        let header = test_header();
        let r1 = fgumi_raw_bam::encode_record_buf_to_raw(&r1_buf, &header).expect("encode r1");
        let r2 = fgumi_raw_bam::encode_record_buf_to_raw(&r2_buf, &header).expect("encode r2");
        let library_index = LibraryIndex::from_header(&header);

        let (info, key) = build_template_info(&r1, &r2, &header, &library_index)
            .expect("build_template_info succeeds")
            .expect("both records map, so Some");

        // TemplateInfo: tags carried through verbatim; ref_name is always R1's ref
        // (chr1 at index 0). The two mates share a reference, so the template span
        // is min(start) .. max(end): min(100, 150) .. max(100+99, 150+99) = 100..249.
        // R1 is forward, so r1_positive is true. hash_fraction is keyed off the read
        // NAME "read-1" (independently pinned by compute_hash_fraction's own tests).
        assert_eq!(info.mi, "7");
        assert_eq!(info.rx, "ACGT-TGCA");
        assert_eq!(info.ref_name.as_deref(), Some("chr1"));
        assert_eq!(info.position, Some(100));
        assert_eq!(info.end_position, Some(249));
        assert!(info.r1_positive, "R1 is on the forward strand");
        // Bit-exact: both sides are the same deterministic murmur3 computation over
        // the read name, so this pins the hash to the read NAME without a tolerance
        // (and satisfies clippy::float_cmp).
        assert_eq!(info.hash_fraction.to_bits(), compute_hash_fraction("read-1").to_bits());

        // ReadInfoKey: mates are ordered so the earlier (ref, unclipped-5', strand)
        // sorts first. R1's forward 5' is its start (100); R2's reverse 5' is its
        // unclipped end (150 + 99 = 249). (0, 100, false) < (0, 249, true), so R1
        // is read 1 in the key. No RG -> library 0; no CB tag -> no cell barcode.
        assert_eq!(key.ref_index1, 0);
        assert_eq!(key.start1, 100);
        assert!(!key.strand1, "read 1 (R1) is forward");
        assert_eq!(key.ref_index2, 0);
        assert_eq!(key.start2, 249);
        assert!(key.strand2, "read 2 (R2) is reverse");
        assert_eq!(key.library, 0, "no @RG in the header, so unknown library");
        assert!(key.cell_barcode.is_none(), "no CB tag was set");
    }

    #[test]
    fn record_simplex_coordinate_group_rejects_mixed_strand_input() {
        // SIMM3-01: a base UMI with reads on BOTH /A and /B strands is
        // duplex data, rejected with a specific error pointing at
        // duplex-metrics.
        let group = vec![
            TemplateInfo {
                mi: "1/A".to_string(),
                rx: "AAA".to_string(),
                ref_name: Some("chr1".to_string()),
                position: Some(100),
                end_position: Some(150),
                r1_positive: true,
                hash_fraction: 0.1,
            },
            TemplateInfo {
                mi: "1/B".to_string(),
                rx: "TTT".to_string(),
                ref_name: Some("chr1".to_string()),
                position: Some(100),
                end_position: Some(150),
                r1_positive: false,
                hash_fraction: 0.2,
            },
        ];
        let fractions = [1.0];
        let mut collectors = vec![SimplexMetricsCollector::new()];
        let mut caller = SimpleUmiConsensusCaller::default();
        let mut counts = vec![0usize];

        let err = record_simplex_coordinate_group(
            &group,
            &fractions,
            &mut collectors,
            &mut caller,
            &mut counts,
        )
        .expect_err("mixed-strand input must be rejected");
        assert!(err.to_string().contains("duplex-UMI data"));
    }

    #[test]
    fn record_simplex_coordinate_group_uses_max_rx_component_count_for_ragged_families() {
        // A `--no-umi` base-UMI family can hold ragged `RX` values: one MI is
        // assigned to every template while their raw `RX` tags are preserved. Here
        // both templates share MI "1" (single strand, so no /A|/B mixing) but carry
        // a 2-component and a 3-component RX. The shorter value appears FIRST, so
        // keying `num_components` off the first entry would drop the trailing "GG"
        // component entirely. The max-length calculation must still record it.
        let group = vec![
            TemplateInfo {
                mi: "1".to_string(),
                rx: "AA-CC".to_string(),
                ref_name: Some("chr1".to_string()),
                position: Some(100),
                end_position: Some(150),
                r1_positive: true,
                hash_fraction: 0.1,
            },
            TemplateInfo {
                mi: "1".to_string(),
                rx: "AA-CC-GG".to_string(),
                ref_name: Some("chr1".to_string()),
                position: Some(100),
                end_position: Some(150),
                r1_positive: true,
                hash_fraction: 0.2,
            },
        ];
        let fractions = [1.0];
        let mut collectors = vec![SimplexMetricsCollector::new()];
        let mut caller = SimpleUmiConsensusCaller::default();
        let mut counts = vec![0usize];

        record_simplex_coordinate_group(
            &group,
            &fractions,
            &mut collectors,
            &mut caller,
            &mut counts,
        )
        .expect("records");

        let umis: Vec<String> = collectors[0].umi_metrics().into_iter().map(|m| m.umi).collect();
        assert!(umis.contains(&"AA".to_string()), "first component recorded: {umis:?}");
        assert!(umis.contains(&"CC".to_string()), "second component recorded: {umis:?}");
        assert!(
            umis.contains(&"GG".to_string()),
            "trailing component of the longer ragged RX must not be dropped: {umis:?}"
        );
    }

    #[test]
    fn record_simplex_coordinate_group_bucketing_is_deterministic_and_cumulative() {
        let group: Vec<TemplateInfo> = (0..20)
            .map(|i| TemplateInfo {
                mi: i.to_string(),
                rx: "AAA".to_string(),
                ref_name: Some("chr1".to_string()),
                position: Some(100),
                end_position: Some(150),
                r1_positive: true,
                hash_fraction: compute_hash_fraction(&format!("read-{i}")),
            })
            .collect();

        let run = || {
            let fractions = DOWNSAMPLING_FRACTIONS;
            let mut collectors: Vec<SimplexMetricsCollector> =
                fractions.iter().map(|_| SimplexMetricsCollector::new()).collect();
            let mut caller = SimpleUmiConsensusCaller::default();
            let mut counts = vec![0usize; fractions.len()];
            record_simplex_coordinate_group(
                &group,
                &fractions,
                &mut collectors,
                &mut caller,
                &mut counts,
            )
            .expect("records");
            counts
        };

        let counts_a = run();
        let counts_b = run();
        assert_eq!(counts_a, counts_b, "bucketing must be deterministic across repeated calls");

        for window in counts_a.windows(2) {
            assert!(
                window[1] >= window[0],
                "cumulative-superset violated: counts must be non-decreasing, got {counts_a:?}"
            );
        }
        assert_eq!(
            *counts_a.last().unwrap(),
            group.len(),
            "100% fraction must include every template"
        );
    }

    #[test]
    fn record_duplex_coordinate_group_rejects_malformed_rx_tag() {
        // DXM-03: an RX tag that doesn't split into exactly 2 '-'-delimited
        // segments is a malformed duplex UMI, rejected with a specific error
        // (matching fgbio's fail-fast MatchError on split("-", -1)).
        let group = vec![
            TemplateInfo {
                mi: "1/A".to_string(),
                rx: "AAA".to_string(), // malformed: no '-' delimiter
                ref_name: Some("chr1".to_string()),
                position: Some(100),
                end_position: Some(150),
                r1_positive: true,
                hash_fraction: 0.1,
            },
            TemplateInfo {
                mi: "1/B".to_string(),
                rx: "AAA".to_string(), // malformed: no '-' delimiter
                ref_name: Some("chr1".to_string()),
                position: Some(100),
                end_position: Some(150),
                r1_positive: false,
                hash_fraction: 0.2,
            },
        ];
        let fractions = [1.0];
        let mut collectors = vec![DuplexMetricsCollector::new(false)];
        let mut caller = SimpleUmiConsensusCaller::default();
        let mut counts = vec![0usize];

        let err = record_duplex_coordinate_group(
            &group,
            &fractions,
            &mut collectors,
            &mut caller,
            &mut counts,
            false,
        )
        .expect_err("malformed RX tag must be rejected");
        assert!(err.to_string().contains("did not contain 2 segments"));
    }
}

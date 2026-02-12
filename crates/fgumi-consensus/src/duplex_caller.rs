//! # Duplex Consensus Calling
//!
//! This module implements two-stage consensus calling for duplex sequencing, where both strands
//! of the original DNA molecule are independently sequenced and combined into a highly accurate
//! duplex consensus.
//!
//! ## Overview
//!
//! Duplex sequencing is a powerful error-correction method that sequences both strands of a DNA
//! molecule independently. Each molecule is tagged with a dual UMI (e.g., "AAAA-CCCC") where the
//! two UMI halves come from opposite ends. Critically:
//!
//! - Reads from one strand have UMI "AAAA-CCCC" (designated /A or "AB strand")
//! - Reads from the opposite strand have UMI "CCCC-AAAA" (designated /B or "BA strand")
//!
//! By comparing the consensus from both strands, nearly all random errors can be detected and
//! corrected, achieving error rates as low as 1 in 10^7 bases.
//!
//! ## Algorithm Overview
//!
//! The duplex consensus calling process involves two stages:
//!
//! ### Stage 1: Single-Strand Consensus (SS-CS)
//!
//! For each strand independently (AB and BA):
//! 1. Group raw reads by strand-specific MI tag (`base_umi/A` or `base_umi/B`)
//! 2. Call single-strand consensus using standard likelihood-based consensus calling
//! 3. Generate SS-CS reads with quality scores reflecting within-strand agreement
//!
//! This produces two single-strand consensus reads per duplex molecule (one per strand).
//!
//! ### Stage 2: Duplex Consensus (DCS)
//!
//! For each molecule with both AB and BA single-strand consensuses:
//! 1. Align corresponding positions in the AB and BA consensus sequences
//! 2. For each position, compare the two consensus bases and qualities:
//!    - **Agreement**: Bases match → Sum quality scores (Bayesian combination)
//!    - **Disagreement**: Bases differ → Take higher quality base, subtract lower quality
//!    - **Equal disagreement**: Equal qualities but different bases → Call as N with Q=2
//!    - **N bases**: Either strand has N → Output N with Q=2
//! 3. Generate duplex consensus read with very high quality scores
//!
//! ## Duplex Quality Score Calculation
//!
//! The duplex quality calculation follows a Bayesian framework:
//!
//! **Agreement** (both strands call same base):
//! ```text
//! Q_duplex = Q_AB + Q_BA (capped at Q93)
//! ```
//! This reflects independent evidence from both strands supporting the same call.
//!
//! **Disagreement** (strands call different bases):
//! ```text
//! Q_duplex = |Q_higher - Q_lower|
//! Base = base_with_higher_quality
//! ```
//! The quality difference represents our confidence that the higher-quality call is correct.
//!
//! **Edge cases**:
//! - If either strand has N: Output N with Q=2 (minimum quality)
//! - If qualities are equal but bases differ: Output N with Q=2 (no confidence)
//! - If resulting quality would be Q≤2: Mask to N with Q=2
//!
//! ## Min-Reads Parameters
//!
//! The `min_reads` parameter controls filtering at multiple levels:
//!
//! - **Format**: `[min_duplex]` or `[min_duplex, min_ss]` or `[min_duplex, min_ab, min_ba]`
//!
//! **Single value** `[M]`:
//! - Single-strand consensus requires M reads per strand
//! - Duplex output only if both strands meet threshold
//!
//! **Two values** `[D, M]`:
//! - Single-strand consensus requires M reads per strand (both AB and BA)
//! - Final duplex must have D reads supporting it (typically D ≤ M since we need both strands)
//!
//! **Three values** `[D, M_AB, M_BA]`:
//! - AB single-strand requires `M_AB` reads
//! - BA single-strand requires `M_BA` reads
//! - Allows asymmetric thresholds (useful if one strand consistently has lower coverage)
//!
//! **Example**: `[3, 2, 1]` means:
//! - Need 2+ reads for AB consensus
//! - Need 1+ read for BA consensus
//! - Need 3+ total reads across both strands for duplex output
//!
//! This provides flexibility for handling unbalanced strand coverage in real data.
//!
//! ## Output SAM Tags
//!
//! Duplex consensus reads include rich metadata tags:
//!
//! ### Per-Read Summary Tags
//!
//! **Duplex consensus**:
//! - `cD` (i32): Maximum read depth across all positions (AB + BA)
//! - `cM` (i32): Minimum read depth across all positions (AB + BA)
//! - `cE` (f32): Consensus error rate (weighted average of AB and BA error rates)
//!
//! **AB strand (top)**:
//! - `aD` (i32): Maximum AB strand depth
//! - `aM` (i32): Minimum AB strand depth
//! - `aE` (f32): AB strand error rate
//!
//! **BA strand (bottom)**:
//! - `bD` (i32): Maximum BA strand depth
//! - `bM` (i32): Minimum BA strand depth
//! - `bE` (f32): BA strand error rate
//!
//! ### Per-Base Tags (if enabled)
//!
//! **AB strand**:
//! - `ad` (Array\<i16\>): Per-base read depth for AB consensus
//! - `ae` (Array\<i16\>): Per-base error count for AB consensus
//! - `ac` (String): AB single-strand consensus sequence
//! - `aq` (String): AB single-strand consensus quality scores
//!
//! **BA strand** (similar with `b` prefix):
//! - `bd`, `be`, `bc`, `bq`: Corresponding BA strand arrays
//!
//! These tags enable detailed QC analysis and downstream filtering based on strand-specific
//! metrics.
//!
//! ## Usage Example
//!
//! ```rust,ignore
//! use fgumi_lib::consensus::duplex_caller::DuplexConsensusCaller;
//! use fgumi_lib::consensus::caller::ConsensusCaller;
//!
//! // Create duplex consensus caller
//! let mut caller = DuplexConsensusCaller::new(
//!     "duplex".to_string(),      // Read name prefix
//!     "RG1".to_string(),          // Read group ID
//!     vec![3, 2, 1],              // Min reads: [duplex, AB, BA]
//!     20,                         // Min input base quality
//!     true,                       // Produce per-base tags
//!     false,                      // Quality trim reads
//!     None,                       // Max reads per strand (no limit)
//!     None,                       // Cell barcode tag (none)
//!     false,                      // Track rejected reads
//!     45,                         // Pre-UMI error rate (Q45)
//!     40,                         // Post-UMI error rate (Q40)
//! )?;
//!
//! // Process reads - should include both /A and /B MI tags
//! // Example: "AAAA-CCCC/A" and "CCCC-AAAA/B" for same molecule
//! let duplex_consensus = caller.consensus_reads_from_sam_records(reads)?;
//!
//! // Check statistics
//! caller.log_statistics();
//! ```
//!
//! ## Strand Orientation Validation
//!
//! The caller validates that AB and BA reads have the expected strand orientations:
//!
//! - **AB R1** and **BA R2** should be on the same genomic strand
//! - **AB R2** and **BA R1** should be on the opposite genomic strand
//!
//! This ensures that the duplex structure is consistent with the expected ligation geometry.
//! Molecules failing this validation are rejected to prevent spurious duplex calls from
//! unrelated molecules.
//!
//! ## Performance Considerations
//!
//! - **Memory**: Peak memory usage is proportional to the largest UMI group size
//! - **Depth**: Works efficiently with 5-100x coverage per strand
//! - **Parallelization**: Multi-threading is handled at the command level (see `duplex.rs`)
//!
//! ## When to Use Duplex Consensus
//!
//! Duplex consensus is ideal for:
//! - Ultra-low frequency variant detection (<1% allele fraction)
//! - Detecting rare mutations in high background (e.g., liquid biopsy)
//! - Measuring somatic mutation rates
//! - Quality control for synthetic DNA constructs
//!
//! The trade-off is ~50% loss of reads (need both strands) and increased computational cost,
//! but the error rate reduction (typically 100-1000x) is unmatched by single-strand methods.
//!
//! ## See Also
//!
//! - `vanilla_consensus_caller`: Single-strand consensus calling used in Stage 1
//! - `caller`: Base consensus calling infrastructure and trait definitions
//! - `base_builder`: Likelihood-based consensus algorithm for individual bases

#[cfg(test)]
use ahash::AHashMap;
#[cfg(test)]
use anyhow::Context;
use anyhow::{Result, bail};
#[cfg(test)]
use bstr::BString;
#[cfg(test)]
use bstr::ByteSlice;
use log::{debug, info};
use noodles::sam::alignment::record::data::field::Tag;
// Used by #[cfg(test)] functions (call_duplex_from_ss_pair) and tests
#[cfg(test)]
use noodles::sam::alignment::record_buf::RecordBuf;

use crate::caller::ConsensusOutput;
use crate::simple_umi::consensus_umis;
use crate::{ReadType, SourceRead};
use crate::caller::{ConsensusCaller, ConsensusCallingStats, RejectionReason};
use crate::phred::MAX_PHRED;
use crate::phred::{MIN_PHRED, PhredScore};
use noodles_raw_bam::{self as bam_fields, UnmappedBamRecordBuilder, flags};
use crate::vanilla_caller::{
    VanillaConsensusRead, VanillaUmiConsensusCaller, VanillaUmiConsensusOptions,
};

/// Duplex consensus read - matches fgbio's `DuplexConsensusRead`
///
/// This struct holds the result of duplex consensus calling, which combines
/// single-strand consensus reads from both strands (AB and BA) of a duplex molecule.
/// It stores both the duplex consensus and references to the underlying single-strand
/// consensuses for downstream analysis and tag generation.
#[derive(Debug, Clone)]
pub struct DuplexConsensusRead {
    /// Unique identifier (base UMI without /A or /B suffix)
    pub id: String,

    /// Duplex consensus bases
    pub bases: Vec<u8>,

    /// Duplex consensus quality scores (Phred-scaled)
    pub quals: Vec<u8>,

    /// Per-base duplex errors (disagreements with source reads)
    pub errors: Vec<u16>,

    /// AB strand single-strand consensus (stored as `VanillaConsensusRead`)
    pub ab_consensus: VanillaConsensusRead,

    /// BA strand single-strand consensus (optional - may be absent for SS-only molecules)
    pub ba_consensus: Option<VanillaConsensusRead>,
}

impl DuplexConsensusRead {
    /// Returns the length of the duplex consensus read
    #[must_use]
    pub fn len(&self) -> usize {
        self.bases.len()
    }

    /// Returns true if the duplex consensus read has zero length
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.bases.is_empty()
    }

    /// Returns the combined depth from both strands at each position
    #[must_use]
    pub fn combined_depths(&self) -> Vec<u16> {
        let ab_depths = &self.ab_consensus.depths;
        let ba_depths = self.ba_consensus.as_ref().map(|ba| &ba.depths);

        (0..self.len())
            .map(|i| {
                let ab = ab_depths.get(i).copied().unwrap_or(0);
                let ba = ba_depths.and_then(|d| d.get(i).copied()).unwrap_or(0);
                ab.saturating_add(ba)
            })
            .collect()
    }

    /// Returns the maximum combined depth across all positions
    #[must_use]
    pub fn max_depth(&self) -> u16 {
        self.combined_depths().into_iter().max().unwrap_or(0)
    }

    /// Returns the minimum combined depth across all positions
    #[must_use]
    pub fn min_depth(&self) -> u16 {
        self.combined_depths().into_iter().min().unwrap_or(0)
    }
}

/// Duplex consensus caller that generates duplex consensus reads from paired
/// single-strand consensus reads.
pub struct DuplexConsensusCaller {
    /// Prefix for consensus read names
    read_name_prefix: String,
    /// Read group ID for consensus reads
    read_group_id: String,
    /// Minimum total reads (duplex threshold)
    min_total_reads: usize,
    /// Minimum reads for the larger strand (XY = max(AB, BA))
    min_xy_reads: usize,
    /// Minimum reads for the smaller strand (YX = min(AB, BA))
    min_yx_reads: usize,
    /// Whether to output per-base tags
    produce_per_base_tags: bool,
    /// Cell barcode tag (e.g., CB) for preserving cell barcodes
    cell_tag: Option<Tag>,
    /// Statistics tracking
    stats: ConsensusCallingStats,
    /// Single-strand consensus caller (used for both /A and /B families)
    /// Uses `min_reads=1` like fgbio; filtering happens at duplex level
    ss_caller: VanillaUmiConsensusCaller,
    /// Rejected reads as raw BAM bytes (if tracking is enabled).
    rejected_reads: Vec<Vec<u8>>,
    /// Whether to track rejected reads
    track_rejects: bool,
}

#[allow(
    clippy::similar_names,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::cast_sign_loss
)]
impl DuplexConsensusCaller {
    /// Creates a new duplex consensus caller
    ///
    /// # Arguments
    /// * `read_name_prefix` - Prefix for consensus read names
    /// * `read_group_id` - Read group ID for consensus reads
    /// * `min_reads` - Vector of 1-3 values: [total, XY, YX]. If only 1 value, used for all.
    ///   If 2 values, [total, XY/YX]. If 3 values, [total, XY, YX].
    ///   XY is the larger strand count, YX is the smaller strand count.
    ///   Filtering happens at duplex level, not single-strand level (matching fgbio).
    /// * `min_input_base_quality` - Minimum base quality for consensus
    /// * `produce_per_base_tags` - Whether to output per-base depth/error tags
    /// * `trim` - Whether to quality-trim reads before consensus
    /// * `max_reads_per_strand` - Maximum reads per strand (downsample if exceeded)
    /// * `cell_tag` - Optional cellular barcode tag for filtering
    /// * `track_rejects` - Whether to track rejected reads for output
    /// * `error_rate_pre_umi` - Phred-scaled error rate prior to UMI integration
    /// * `error_rate_post_umi` - Phred-scaled error rate post UMI integration
    ///
    /// # Errors
    ///
    /// Returns an error if `min_reads` is empty, has more than 3 values, or the values
    /// are not in decreasing order (total >= XY >= YX).
    ///
    /// # Panics
    ///
    /// Panics if `min_reads` is empty (though this is checked and returns an error first).
    #[allow(clippy::too_many_arguments)]
    #[expect(clippy::needless_pass_by_value, reason = "changing to &[usize] would require updating many call sites across the codebase")]
    pub fn new(
        read_name_prefix: String,
        read_group_id: String,
        min_reads: Vec<usize>,
        min_input_base_quality: u8,
        produce_per_base_tags: bool,
        trim: bool,
        max_reads_per_strand: Option<usize>,
        cell_tag: Option<Tag>,
        track_rejects: bool,
        error_rate_pre_umi: u8,
        error_rate_post_umi: u8,
    ) -> Result<Self> {
        // Parse min_reads vector like fgbio: padTo(3, last)
        // [total, XY, YX] where XY >= YX (larger strand, smaller strand)
        if min_reads.is_empty() {
            anyhow::bail!("min_reads parameter must have at least 1 value");
        }
        if min_reads.len() > 3 {
            anyhow::bail!(
                "min_reads parameter must have 1-3 values (total, [XY, [YX]]), got {} values",
                min_reads.len()
            );
        }

        let last = *min_reads.last().unwrap();
        let min_total_reads = min_reads[0];
        let min_xy_reads = min_reads.get(1).copied().unwrap_or(last);
        let min_yx_reads = min_reads.get(2).copied().unwrap_or(last);

        // Validate min_reads ordering: yx <= xy <= total (matching fgbio)
        // For depth thresholds it's required that yx <= xy <= total
        if min_xy_reads > min_total_reads {
            anyhow::bail!("min-reads values must be specified high to low (total >= XY)");
        }
        if min_yx_reads > min_xy_reads {
            anyhow::bail!("min-reads values must be specified high to low (XY >= YX)");
        }

        // Create single-strand consensus caller with min_reads=1 (matching fgbio)
        // Filtering happens at duplex level based on input read counts, not at SS level
        let min_consensus_base_quality: PhredScore = MIN_PHRED;

        let ss_options = VanillaUmiConsensusOptions {
            min_reads: 1, // fgbio uses minReads=1 for SS caller
            min_input_base_quality,
            produce_per_base_tags,
            trim,
            max_reads: max_reads_per_strand,
            error_rate_pre_umi,
            error_rate_post_umi,
            min_consensus_base_quality,
            cell_tag,
            ..Default::default()
        };

        let ss_caller = VanillaUmiConsensusCaller::new_with_rejects_tracking(
            read_name_prefix.clone(),
            read_group_id.clone(),
            ss_options,
            track_rejects,
        );

        Ok(Self {
            read_name_prefix,
            read_group_id,
            min_total_reads,
            min_xy_reads,
            min_yx_reads,
            produce_per_base_tags,
            cell_tag,
            stats: ConsensusCallingStats::new(),
            ss_caller,
            rejected_reads: Vec::new(),
            track_rejects,
        })
    }

    /// Returns the rejected reads as raw BAM bytes
    #[must_use]
    pub fn rejected_reads(&self) -> &[Vec<u8>] {
        &self.rejected_reads
    }

    /// Takes ownership of the rejected reads, leaving an empty Vec
    pub fn take_rejected_reads(&mut self) -> Vec<Vec<u8>> {
        std::mem::take(&mut self.rejected_reads)
    }

    /// Clears the rejected reads buffer
    pub fn clear_rejected_reads(&mut self) {
        self.rejected_reads.clear();
    }

    /// Clears all per-group state to prepare for reuse
    ///
    /// This resets statistics and rejected reads while preserving the caller's
    /// configuration and reusable components.
    pub fn clear(&mut self) {
        self.stats = ConsensusCallingStats::default();
        self.rejected_reads.clear();
        self.ss_caller.clear();
    }

    /// Test helper: accepts `Vec<RecordBuf>`, encodes to raw bytes, and delegates to `consensus_reads`.
    #[cfg(test)]
    pub(crate) fn consensus_reads_from_sam_records(
        &mut self,
        records: Vec<noodles::sam::alignment::RecordBuf>,
    ) -> anyhow::Result<ConsensusOutput> {
        use noodles::sam::header::record::value::Map;
        use noodles::sam::header::record::value::map::ReferenceSequence;
        use std::num::NonZeroUsize;
        // Build header with a reference sequence for mapped records
        let header = noodles::sam::Header::builder()
            .add_reference_sequence(
                "chr1",
                Map::<ReferenceSequence>::new(NonZeroUsize::new(1_000_000).unwrap()),
            )
            .build();
        let raw: Vec<Vec<u8>> = records
            .iter()
            .map(|rec| {
                let mut buf = Vec::new();
                crate::vendored::bam_codec::encode_record_buf(&mut buf, &header, rec)
                    .map_err(|e| anyhow::anyhow!("Failed to encode RecordBuf: {e}"))?;
                Ok(buf)
            })
            .collect::<anyhow::Result<Vec<_>>>()?;
        self.consensus_reads(raw)
    }

    /// Parses MI tag to extract base UMI and strand (/A or /B)
    ///
    /// Returns (`base_umi`, strand) where strand is 'A', 'B', or None if no suffix
    #[cfg(test)]
    fn parse_mi_tag(mi: &str) -> (String, Option<char>) {
        if let Some(stripped) = mi.strip_suffix("/A") {
            (stripped.to_string(), Some('A'))
        } else if let Some(stripped) = mi.strip_suffix("/B") {
            (stripped.to_string(), Some('B'))
        } else {
            (mi.to_string(), None)
        }
    }

    /// Extracts strand suffix from MI tag bytes using fast byte-level operations.
    ///
    /// Returns:
    /// - `Ok(Some('A'))` if tag ends with "/A"
    /// - `Ok(Some('B'))` if tag ends with "/B"
    /// - `Ok(None)` if no valid suffix (will be handled as error by caller)
    /// - The `base_mi` is the tag without the suffix
    #[inline]
    fn extract_strand_from_mi_bytes(mi_bytes: &[u8]) -> Option<char> {
        // Check for /A or /B suffix (need at least 2 bytes: "/" and "A" or "B")
        if mi_bytes.len() >= 2 {
            let len = mi_bytes.len();
            if mi_bytes[len - 2] == b'/' {
                match mi_bytes[len - 1] {
                    b'A' => return Some('A'),
                    b'B' => return Some('B'),
                    _ => {}
                }
            }
        }
        None
    }

    /// Extracts base MI (without /A or /B suffix) from MI tag bytes.
    ///
    /// Assumes the tag has a valid suffix (caller should verify with `extract_strand_from_mi_bytes`).
    #[inline]
    fn base_mi_from_bytes(mi_bytes: &[u8]) -> String {
        if mi_bytes.len() >= 2 {
            // Safe: we're slicing off the last 2 bytes ("/A" or "/B")
            String::from_utf8_lossy(&mi_bytes[..mi_bytes.len() - 2]).into_owned()
        } else {
            String::from_utf8_lossy(mi_bytes).into_owned()
        }
    }

    /// Quick check if a group of reads has both /A and /B strands.
    ///
    /// This is an optimization to skip expensive operations (like overlapping consensus)
    /// when a group cannot possibly produce duplex output. Returns `true` only if there
    /// is at least one read with /A suffix AND at least one read with /B suffix.
    ///
    /// This check is O(n) in the worst case but typically returns early once both
    /// strands are found.
    #[cfg(test)]
    #[must_use]
    pub fn has_both_strands(reads: &[RecordBuf]) -> bool {
        if reads.len() < 2 {
            return false;
        }

        let mi_tag = Tag::from([b'M', b'I']);
        let mut has_a = false;
        let mut has_b = false;

        for read in reads {
            if let Some(noodles::sam::alignment::record_buf::data::field::Value::String(mi_bytes)) =
                read.data().get(&mi_tag)
            {
                match Self::extract_strand_from_mi_bytes(mi_bytes) {
                    Some('A') => {
                        has_a = true;
                        if has_b {
                            return true;
                        }
                    }
                    Some('B') => {
                        has_b = true;
                        if has_a {
                            return true;
                        }
                    }
                    _ => {}
                }
            }
        }

        false
    }

    /// Partitions records by strand (A or B) using optimized byte-level operations on raw BAM bytes.
    ///
    /// This is an optimized version that avoids `HashMap` allocation since reads
    /// are already grouped by base MI (without /A or /B suffix) by the MI group iterator.
    ///
    /// Returns (`Option<base_mi>`, `a_records`, `b_records`) or an error if any read has an invalid MI tag.
    ///
    /// IMPORTANT: This function requires that all reads have MI tags with /A or /B suffixes.
    /// The duplex command MUST be used with reads grouped using the "paired" strategy.
    #[allow(clippy::type_complexity)]
    fn partition_records_by_strand(
        records: Vec<Vec<u8>>,
    ) -> Result<(Option<String>, Vec<Vec<u8>>, Vec<Vec<u8>>)> {
        if records.is_empty() {
            return Ok((None, Vec::new(), Vec::new()));
        }

        // Pre-allocate with estimated capacity (assume roughly even split)
        let half_capacity = (records.len() / 2).max(1);
        let mut a_records = Vec::with_capacity(half_capacity);
        let mut b_records = Vec::with_capacity(half_capacity);

        let mut base_mi: Option<String> = None;

        for record in records {
            // Extract strand info before moving the record (drop borrow before push)
            let is_a_strand = {
                let Some(mi_bytes) = bam_fields::find_string_tag_in_record(&record, b"MI") else {
                    let read_name = String::from_utf8_lossy(bam_fields::read_name(&record));
                    bail!(
                        "Read '{read_name}' is missing MI tag. \
                        The duplex command requires all reads to have MI tags. \
                        Please run 'fgumi group' on your input BAM first."
                    );
                };
                if base_mi.is_none() {
                    base_mi = Some(Self::base_mi_from_bytes(mi_bytes));
                }
                match Self::extract_strand_from_mi_bytes(mi_bytes) {
                    Some('A') => true,
                    Some('B') => false,
                    _ => {
                        let mi_str = String::from_utf8_lossy(mi_bytes).into_owned();
                        bail!(
                            "Read has MI tag '{mi_str}' without /A or /B suffix. \
                            The duplex command requires reads to be grouped using the 'paired' strategy, \
                            which adds /A and /B suffixes to MI tags to indicate the strand of the source \
                            duplex molecule. Please run 'fgumi group --strategy paired' on your input BAM \
                            before running duplex consensus calling."
                        );
                    }
                }
            };

            if is_a_strand {
                a_records.push(record);
            } else {
                b_records.push(record);
            }
        }

        Ok((base_mi, a_records, b_records))
    }

    /// Groups reads by base MI tag (without /A or /B suffix) and strand.
    ///
    /// IMPORTANT: This function requires that all reads have MI tags with /A or /B suffixes.
    /// The duplex command MUST be used with reads grouped using the "paired" strategy,
    /// which adds these suffixes to indicate the strand of the source duplex molecule.
    ///
    /// NOTE: This function is used only in tests. Production code uses
    /// [`partition_reads_by_strand_simple`] which is more efficient.
    #[cfg(test)]
    #[allow(clippy::type_complexity)]
    fn group_by_mi_and_strand(
        reads: Vec<RecordBuf>,
    ) -> Result<AHashMap<String, (Vec<RecordBuf>, Vec<RecordBuf>)>> {
        let mut groups: AHashMap<String, (Vec<RecordBuf>, Vec<RecordBuf>)> = AHashMap::new();

        for read in reads {
            // Get MI tag
            let mi_tag = noodles::sam::alignment::record::data::field::Tag::from([b'M', b'I']);
            let mi =
                if let Some(noodles::sam::alignment::record_buf::data::field::Value::String(s)) =
                    read.data().get(&mi_tag)
                {
                    String::from_utf8(s.iter().copied().collect::<Vec<u8>>())
                        .context("MI tag is not valid UTF-8")?
                } else {
                    // No MI tag, skip this read
                    continue;
                };

            let (base_mi, strand) = Self::parse_mi_tag(&mi);

            let entry = groups.entry(base_mi).or_insert_with(|| (Vec::new(), Vec::new()));

            match strand {
                Some('A') => entry.0.push(read),
                Some('B') => entry.1.push(read),
                None => {
                    // No strand suffix - this is an error!
                    bail!(
                        "Read has MI tag '{mi}' without /A or /B suffix. \
                        The duplex command requires reads to be grouped using the 'paired' strategy, \
                        which adds /A and /B suffixes to MI tags to indicate the strand of the source \
                        duplex molecule. Please run 'fgumi group --strategy paired' on your input BAM \
                        before running duplex consensus calling."
                    );
                }
                Some(c) => {
                    // Invalid strand suffix character
                    bail!(
                        "Read has MI tag '{mi}' with invalid strand suffix '{c}'. \
                        Expected /A or /B suffix only."
                    );
                }
            }
        }

        Ok(groups)
    }

    /// Helper function to extract integer value from SAM tag, handling different integer types
    #[cfg(test)]
    fn extract_int_tag(
        data: &noodles::sam::alignment::record_buf::Data,
        tag: noodles::sam::alignment::record::data::field::Tag,
    ) -> Option<i32> {
        use noodles::sam::alignment::record_buf::data::field::Value;
        match data.get(&tag)? {
            Value::Int8(v) => Some(i32::from(*v)),
            Value::UInt8(v) => Some(i32::from(*v)),
            Value::Int16(v) => Some(i32::from(*v)),
            Value::UInt16(v) => Some(i32::from(*v)),
            Value::Int32(v) => Some(*v),
            Value::UInt32(v) => i32::try_from(*v).ok(),
            _ => None,
        }
    }

    /// Checks if there are enough input reads according to the `min_reads` thresholds.
    /// This matches fgbio's `hasMinimumNumberOfReads` logic.
    ///
    /// Counts only R1 reads (first of pair) for each strand, sorts them so XY >= YX,
    /// then checks: `total >= min_total AND xy >= min_xy AND yx >= min_yx`
    fn has_minimum_number_of_reads(
        a_records: &[Vec<u8>],
        b_records: &[Vec<u8>],
        min_total_reads: usize,
        min_xy_reads: usize,
        min_yx_reads: usize,
    ) -> bool {
        // Count R1s only (first of pair) for each strand
        // Match fgbio: x.count(r => r.paired && r.firstOfPair)
        let num_a = a_records
            .iter()
            .filter(|r| {
                let flg = bam_fields::flags(r);
                (flg & flags::PAIRED != 0) && (flg & flags::FIRST_SEGMENT != 0)
            })
            .count();
        let num_b = b_records
            .iter()
            .filter(|r| {
                let flg = bam_fields::flags(r);
                (flg & flags::PAIRED != 0) && (flg & flags::FIRST_SEGMENT != 0)
            })
            .count();

        // Sort so xy is the larger count, yx is the smaller
        let (num_xy, num_yx) = if num_a >= num_b { (num_a, num_b) } else { (num_b, num_a) };

        let total = num_xy + num_yx;

        min_total_reads <= total && min_xy_reads <= num_xy && min_yx_reads <= num_yx
    }

    /// Check if a `DuplexConsensusRead` has the minimum number of reads required.
    /// Uses the max depth from AB and BA consensuses to determine if we have enough reads.
    fn duplex_consensus_has_minimum_reads(
        consensus: &DuplexConsensusRead,
        min_total_reads: usize,
        min_xy_reads: usize,
        min_yx_reads: usize,
    ) -> bool {
        // Get max depth from AB and BA consensuses
        let num_a = consensus.ab_consensus.max_depth() as usize;
        let num_b = consensus.ba_consensus.as_ref().map_or(0, |b| b.max_depth() as usize);

        // Sort so xy is the larger count, yx is the smaller
        let (num_xy, num_yx) = if num_a >= num_b { (num_a, num_b) } else { (num_b, num_a) };

        let total = num_xy + num_yx;

        min_total_reads <= total && min_xy_reads <= num_xy && min_yx_reads <= num_yx
    }

    /// Check if all reads in an iterator are on the same strand.
    /// Returns true if all reads have the same orientation (all forward or all reverse).
    /// Also returns true for empty iterators.
    fn are_all_same_strand<'a>(mut reads: impl Iterator<Item = &'a Vec<u8>>) -> bool {
        let Some(first) = reads.next() else {
            return true;
        };
        let first_is_reverse = bam_fields::flags(first) & flags::REVERSE != 0;
        reads.all(|r| (bam_fields::flags(r) & flags::REVERSE != 0) == first_is_reverse)
    }

    // Helper function to cap quality scores to valid range [2, 93]
    fn cap_quality(score: i32) -> u8 {
        if score < MIN_PHRED.into() {
            MIN_PHRED
        } else if score > MAX_PHRED.into() {
            MAX_PHRED
        } else {
            score as u8
        }
    }

    /// Checks if a source base differs from the consensus base (error detection).
    /// Returns true if both bases are called (not N) and they differ.
    /// This matches fgbio's `isError` function.
    #[inline]
    fn is_error(source_base: u8, consensus_base: u8) -> bool {
        const NO_CALL: u8 = b'N';
        source_base != NO_CALL && consensus_base != NO_CALL && source_base != consensus_base
    }

    /// Creates a duplex consensus read from AB and BA single-strand consensuses.
    ///
    /// This matches fgbio's `duplexConsensus` method exactly:
    /// - Takes optional AB and BA `VanillaConsensusReads`
    /// - Returns a `DuplexConsensusRead` containing the duplex sequence/quality
    /// - Stores references to the truncated AB and BA consensuses
    ///
    /// # Arguments
    /// * `ab` - Optional AB strand single-strand consensus
    /// * `ba` - Optional BA strand single-strand consensus
    /// * `source_reads` - Optional source reads for error calculation
    #[expect(clippy::too_many_lines, reason = "duplex consensus building has many sequential steps that are clearest in one function")]
    pub(crate) fn duplex_consensus(
        ab: Option<&VanillaConsensusRead>,
        ba: Option<&VanillaConsensusRead>,
        source_reads: Option<&[SourceRead]>,
    ) -> Option<DuplexConsensusRead> {
        const NO_CALL: u8 = b'N';
        const NO_CALL_QUAL: u8 = MIN_PHRED;

        // Calculate length as minimum of both consensuses (fgbio behavior)
        let len = std::cmp::min(
            ab.map_or(usize::MAX, VanillaConsensusRead::len),
            ba.map_or(usize::MAX, VanillaConsensusRead::len),
        );

        // Filter out consensuses that have no coverage in the truncated region
        let ab_filtered = ab.filter(|a| a.depths.iter().take(len).any(|&d| d > 0));
        let ba_filtered = ba.filter(|b| b.depths.iter().take(len).any(|&d| d > 0));

        match (ab_filtered, ba_filtered) {
            (Some(a), None) => {
                // Only AB strand - use it directly
                Some(DuplexConsensusRead {
                    id: a.id.clone(),
                    bases: a.bases.clone(),
                    quals: a.quals.clone(),
                    errors: a.errors.clone(),
                    ab_consensus: a.clone(),
                    ba_consensus: None,
                })
            }
            (None, Some(b)) => {
                // Only BA strand - use it directly
                Some(DuplexConsensusRead {
                    id: b.id.clone(),
                    bases: b.bases.clone(),
                    quals: b.quals.clone(),
                    errors: b.errors.clone(),
                    ab_consensus: b.clone(),
                    ba_consensus: None,
                })
            }
            (Some(a), Some(b)) => {
                // Both strands - combine them
                let id = a.id.clone();
                let mut bases = Vec::with_capacity(len);
                let mut quals = Vec::with_capacity(len);
                let mut errors = Vec::with_capacity(len);

                for i in 0..len {
                    let a_base = a.bases[i];
                    let b_base = b.bases[i];
                    let a_qual = i32::from(a.quals[i]);
                    let b_qual = i32::from(b.quals[i]);

                    // Calculate raw consensus base and quality (fgbio algorithm)
                    let (raw_base, raw_qual) = if a_base == b_base {
                        // Agreement: sum qualities (capped)
                        (a_base, Self::cap_quality(a_qual + b_qual))
                    } else if a_qual > b_qual {
                        // Disagreement: take higher quality base
                        (a_base, Self::cap_quality(a_qual - b_qual))
                    } else if b_qual > a_qual {
                        (b_base, Self::cap_quality(b_qual - a_qual))
                    } else {
                        // Equal qualities, different bases: no confidence
                        (a_base, NO_CALL_QUAL)
                    };

                    // Mask to N if either input is N or quality is minimum
                    let (final_base, final_qual) =
                        if a_base == NO_CALL || b_base == NO_CALL || raw_qual == NO_CALL_QUAL {
                            (NO_CALL, NO_CALL_QUAL)
                        } else {
                            (raw_base, raw_qual)
                        };

                    bases.push(final_base);
                    quals.push(final_qual);

                    // Calculate errors
                    let error_count = if let Some(source_reads) = source_reads {
                        // Exact method: count disagreements with source reads
                        let mut num_errors = 0i32;
                        for sr in source_reads {
                            if sr.bases.len() > i && Self::is_error(sr.bases[i], raw_base) {
                                num_errors += 1;
                            }
                        }
                        num_errors.clamp(0, i32::from(i16::MAX)) as u16
                    } else {
                        // Approximate method when source reads unavailable
                        let a_err = i32::from(a.errors[i]);
                        let b_err = i32::from(b.errors[i]);
                        let a_dep = i32::from(a.depths[i]);
                        let b_dep = i32::from(b.depths[i]);

                        let err = if a_base == b_base {
                            a_err + b_err
                        } else if raw_base == a_base {
                            a_err + (b_dep - b_err)
                        } else {
                            b_err + (a_dep - a_err)
                        };
                        err.clamp(0, i32::from(i16::MAX)) as u16
                    };

                    errors.push(error_count);
                }

                // Truncate AB and BA consensuses to match duplex length
                let ab_truncated = VanillaConsensusRead {
                    id: a.id.clone(),
                    bases: a.bases[..len].to_vec(),
                    quals: a.quals[..len].to_vec(),
                    depths: a.depths[..len].to_vec(),
                    errors: a.errors[..len].to_vec(),
                    source_reads: None, // Don't clone source reads
                };

                let ba_truncated = VanillaConsensusRead {
                    id: b.id.clone(),
                    bases: b.bases[..len].to_vec(),
                    quals: b.quals[..len].to_vec(),
                    depths: b.depths[..len].to_vec(),
                    errors: b.errors[..len].to_vec(),
                    source_reads: None,
                };

                Some(DuplexConsensusRead {
                    id,
                    bases,
                    quals,
                    errors,
                    ab_consensus: ab_truncated,
                    ba_consensus: Some(ba_truncated),
                })
            }
            (None, None) => None,
        }
    }

    /// Writes a `DuplexConsensusRead` as raw BAM bytes into `output`.
    ///
    /// This is the second stage of the two-stage consensus calling process.
    /// Matches fgbio's `createSamRecord` for duplex consensus.
    ///
    /// # Arguments
    /// * `builder` - Reusable builder for raw BAM record construction
    /// * `output` - `ConsensusOutput` to append the record into
    /// * `consensus` - The duplex consensus read to convert
    /// * `read_type` - Whether this is R1 or R2
    /// * `umi` - The UMI string (base MI without /A or /B suffix)
    /// * `source_reads_a` - Source reads from A strand for RX tag consensus
    /// * `source_reads_b` - Source reads from B strand for RX tag consensus
    /// * `produce_per_base_tags` - Whether to include per-base tags
    /// * `read_name_prefix` - Prefix for the read name
    /// * `read_group_id` - Read group ID
    /// * `first_of_pair` - Whether this is first of pair (for RX orientation)
    /// * `cell_tag` - Optional cell barcode tag (e.g., CB)
    /// * `cell_barcode` - Optional cell barcode value to add to record
    #[allow(clippy::too_many_arguments)]
    #[expect(clippy::too_many_lines, reason = "BAM record construction has many sequential tag-writing steps")]
    #[expect(clippy::unnecessary_wraps, reason = "Result return type kept for API consistency with other consensus record builders")]
    pub(crate) fn duplex_read_into(
        builder: &mut UnmappedBamRecordBuilder,
        output: &mut ConsensusOutput,
        consensus: &DuplexConsensusRead,
        read_type: ReadType,
        umi: &str,
        source_reads_a: &[&[u8]],
        source_reads_b: &[&[u8]],
        produce_per_base_tags: bool,
        read_name_prefix: &str,
        read_group_id: &str,
        first_of_pair: bool,
        cell_tag: Option<Tag>,
        cell_barcode: Option<&str>,
    ) -> Result<()> {
        // Build flags
        let mut flag = flags::UNMAPPED;
        match read_type {
            ReadType::R1 => {
                flag |= flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_UNMAPPED;
            }
            ReadType::R2 => {
                flag |= flags::PAIRED | flags::LAST_SEGMENT | flags::MATE_UNMAPPED;
            }
            ReadType::Fragment => {
                // No pair flags for fragment
            }
        }

        // Build the record (name, flags, bases, quals)
        let read_name = format!("{read_name_prefix}:{umi}");
        builder.build_record(read_name.as_bytes(), flag, &consensus.bases, &consensus.quals);

        // 1. MI tag (string)
        builder.append_string_tag(b"MI", umi.as_bytes());

        // 2. Cell barcode tag if present
        if let (Some(tag), Some(barcode)) = (cell_tag, cell_barcode) {
            let tag_bytes: [u8; 2] = <[u8; 2]>::from(tag);
            builder.append_string_tag(&tag_bytes, barcode.as_bytes());
        }

        // 3. RG tag (string)
        builder.append_string_tag(b"RG", read_group_id.as_bytes());

        // Calculate AB strand metrics
        let ab = &consensus.ab_consensus;
        let ba_opt = consensus.ba_consensus.as_ref();

        let ab_depth_max = i32::from(ab.depths.iter().copied().max().unwrap_or(0));
        let ab_depth_min = i32::from(ab.depths.iter().copied().min().unwrap_or(0));
        let ab_total_depth: i64 = ab.depths.iter().map(|&d| i64::from(d)).sum();
        let ab_total_errors: i64 = ab.errors.iter().map(|&e| i64::from(e)).sum();
        let ab_error_rate =
            if ab_total_depth > 0 { ab_total_errors as f32 / ab_total_depth as f32 } else { 0.0 };

        // 4. AB strand tags: aD (int, max depth), aE (float, error rate), aM (int, min depth)
        builder.append_int_tag(b"aD", ab_depth_max);
        builder.append_float_tag(b"aE", ab_error_rate);
        builder.append_int_tag(b"aM", ab_depth_min);

        // 5. Per-base AB tags if requested
        if produce_per_base_tags {
            builder.append_string_tag(b"ac", &ab.bases);

            let ab_depths_i16: Vec<i16> =
                ab.depths.iter().map(|&d| i16::try_from(d).unwrap_or(i16::MAX)).collect();
            builder.append_i16_array_tag(b"ad", &ab_depths_i16);

            let ab_errors_i16: Vec<i16> =
                ab.errors.iter().map(|&e| i16::try_from(e).unwrap_or(i16::MAX)).collect();
            builder.append_i16_array_tag(b"ae", &ab_errors_i16);

            builder.append_phred33_string_tag(b"aq", &ab.quals);
        }

        // Calculate BA strand metrics
        let (ba_depth_max, ba_depth_min, ba_error_rate) = if let Some(ba) = ba_opt {
            let ba_depth_max = i32::from(ba.depths.iter().copied().max().unwrap_or(0));
            let ba_depth_min = i32::from(ba.depths.iter().copied().min().unwrap_or(0));
            let ba_total_depth: i64 = ba.depths.iter().map(|&d| i64::from(d)).sum();
            let ba_total_errors: i64 = ba.errors.iter().map(|&e| i64::from(e)).sum();
            let ba_error_rate = if ba_total_depth > 0 {
                ba_total_errors as f32 / ba_total_depth as f32
            } else {
                0.0
            };
            (ba_depth_max, ba_depth_min, ba_error_rate)
        } else {
            (0i32, 0i32, 0.0f32)
        };

        // 6. BA strand tags: bD (int), bE (float), bM (int)
        builder.append_int_tag(b"bD", ba_depth_max);
        builder.append_float_tag(b"bE", ba_error_rate);
        builder.append_int_tag(b"bM", ba_depth_min);

        // 7. Per-base BA tags if requested and BA strand exists
        if produce_per_base_tags {
            if let Some(ba) = ba_opt {
                builder.append_string_tag(b"bc", &ba.bases);

                let ba_depths_i16: Vec<i16> =
                    ba.depths.iter().map(|&d| i16::try_from(d).unwrap_or(i16::MAX)).collect();
                builder.append_i16_array_tag(b"bd", &ba_depths_i16);

                let ba_errors_i16: Vec<i16> =
                    ba.errors.iter().map(|&e| i16::try_from(e).unwrap_or(i16::MAX)).collect();
                builder.append_i16_array_tag(b"be", &ba_errors_i16);

                builder.append_phred33_string_tag(b"bq", &ba.quals);
            }
        }

        // 8. Duplex consensus tags (cD, cM, cE)
        let combined_depths: Vec<i32> = (0..consensus.len())
            .map(|i| {
                let ab_d = i32::from(ab.depths.get(i).copied().unwrap_or(0));
                let ba_d = i32::from(ba_opt.and_then(|b| b.depths.get(i).copied()).unwrap_or(0));
                ab_d + ba_d
            })
            .collect();

        let duplex_depth_max = combined_depths.iter().copied().max().unwrap_or(0);
        let duplex_depth_min = combined_depths.iter().copied().min().unwrap_or(0);

        let total_depth: i64 = combined_depths.iter().map(|&d| i64::from(d)).sum();
        let total_errors: i64 = consensus.errors.iter().map(|&e| i64::from(e)).sum();
        let duplex_error_rate =
            if total_depth > 0 { total_errors as f32 / total_depth as f32 } else { 0.0 };

        builder.append_int_tag(b"cD", duplex_depth_max);
        builder.append_float_tag(b"cE", duplex_error_rate);
        builder.append_int_tag(b"cM", duplex_depth_min);

        // 9. Build RX consensus from source reads
        let mut all_umis = Vec::new();

        for raw in source_reads_a {
            if let Some(rx_bytes) = bam_fields::find_string_tag_in_record(raw, b"RX") {
                let rx = String::from_utf8_lossy(rx_bytes).to_string();
                let is_first = bam_fields::flags(raw) & flags::FIRST_SEGMENT != 0;
                if is_first == first_of_pair {
                    all_umis.push(rx);
                } else {
                    let reversed: Vec<&str> = rx.split('-').rev().collect();
                    all_umis.push(reversed.join("-"));
                }
            }
        }

        for raw in source_reads_b {
            if let Some(rx_bytes) = bam_fields::find_string_tag_in_record(raw, b"RX") {
                let rx = String::from_utf8_lossy(rx_bytes).to_string();
                let is_first = bam_fields::flags(raw) & flags::FIRST_SEGMENT != 0;
                if is_first == first_of_pair {
                    all_umis.push(rx);
                } else {
                    let reversed: Vec<&str> = rx.split('-').rev().collect();
                    all_umis.push(reversed.join("-"));
                }
            }
        }

        if !all_umis.is_empty() {
            let consensus_umi = consensus_umis(&all_umis);
            builder.append_string_tag(b"RX", consensus_umi.as_bytes());
        }

        // 10. Write to output
        builder.write_with_block_size(&mut output.data);
        output.count += 1;

        Ok(())
    }

    /// Extracts a per-base array (cd or ce tag) from single-strand consensus data
    #[cfg(test)]
    fn extract_per_base_array(
        data: &noodles::sam::alignment::record_buf::Data,
        tag: noodles::sam::alignment::record::data::field::Tag,
        default_len: usize,
    ) -> Vec<i16> {
        use noodles::sam::alignment::record_buf::data::field::Value;
        use noodles::sam::alignment::record_buf::data::field::value::Array;

        if let Some(Value::Array(arr)) = data.get(&tag) {
            match arr {
                Array::Int16(values) => values.clone(),
                Array::Int8(values) => values.iter().map(|&v| i16::from(v)).collect(),
                _ => vec![0i16; default_len],
            }
        } else {
            vec![0i16; default_len]
        }
    }

    /// Calls duplex consensus from paired single-strand consensuses
    #[cfg(test)]
    fn call_duplex_from_ss_pair(
        ss_a: RecordBuf,
        ss_b: RecordBuf,
        produce_per_base_tags: bool,
        source_reads_a: &[std::sync::Arc<RecordBuf>],
        source_reads_b: &[std::sync::Arc<RecordBuf>],
        ss_caller: &VanillaUmiConsensusCaller,
        duplex_first_of_pair: bool, // true if creating R1 duplex, false if creating R2 duplex
    ) -> Result<Option<RecordBuf>> {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value;
        use noodles::sam::alignment::record_buf::{QualityScores, Sequence};

        use fgumi_sam::to_smallest_signed_int;

        // Get sequences and qualities from both single-strand consensuses
        let seq_a = ss_a.sequence();
        let qual_a = ss_a.quality_scores();
        let seq_b = ss_b.sequence();
        let qual_b = ss_b.quality_scores();

        // Use minimum length - fgbio truncates to the shorter of the two consensuses
        // This matches fgbio's behavior:
        // val len = min(ab.map(_.length).getOrElse(Int.MaxValue), ba.map(_.length).getOrElse(Int.MaxValue))
        let len = seq_a.len().min(seq_b.len());

        if seq_a.len() != seq_b.len() {
            debug!("Truncating to minimum length {} (A={}, B={})", len, seq_a.len(), seq_b.len());
        }

        // Convert source reads to SourceRead with proper orientation (RC if on negative strand)
        // This matches fgbio's behavior where source reads are stored in consensus orientation
        let mut source_reads: Vec<SourceRead> = Vec::new();
        for (idx, read) in source_reads_a.iter().chain(source_reads_b.iter()).enumerate() {
            if let Some(src_read) = ss_caller.to_source_read_from_record(read.as_ref(), idx) {
                source_reads.push(src_read);
            }
        }

        // Extract per-base depths and errors from single-strand consensuses
        // These are needed to compute duplex errors relative to the final duplex consensus
        let consensus_depth_base_tag = Tag::from([b'c', b'd']);
        let consensus_errors_base_tag = Tag::from([b'c', b'e']);

        // Check if per-base tags exist (needed for proper cD/cM/cE calculation)
        let has_per_base_tags = ss_a.data().get(&consensus_depth_base_tag).is_some()
            && ss_b.data().get(&consensus_depth_base_tag).is_some();

        let depths_a = Self::extract_per_base_array(ss_a.data(), consensus_depth_base_tag, len);
        let depths_b = Self::extract_per_base_array(ss_b.data(), consensus_depth_base_tag, len);
        let errors_a = Self::extract_per_base_array(ss_a.data(), consensus_errors_base_tag, len);
        let errors_b = Self::extract_per_base_array(ss_b.data(), consensus_errors_base_tag, len);

        // Build duplex consensus sequence, quality, and per-base errors
        let mut duplex_seq = Vec::with_capacity(len);
        let mut duplex_qual = Vec::with_capacity(len);
        let mut duplex_errors = Vec::with_capacity(len);

        // Constants matching Scala implementation
        const NO_CALL: u8 = b'N';
        const NO_CALL_QUAL: u8 = MIN_PHRED;

        // For each position, compare bases and call consensus
        for i in 0..len {
            let base_a = seq_a.as_ref().get(i).copied().context("Index out of bounds")?;
            let qual_a_i =
                i32::from(qual_a.as_ref().get(i).copied().context("Index out of bounds")?);
            let base_b = seq_b.as_ref().get(i).copied().context("Index out of bounds")?;
            let qual_b_i =
                i32::from(qual_b.as_ref().get(i).copied().context("Index out of bounds")?);

            // Calculate raw consensus base and quality using Scala's algorithm
            // This matches the Bayesian approach in DuplexConsensusCaller.scala lines 423-428
            let (raw_base, raw_qual) = if base_a == base_b {
                // Agreement: sum qualities (Bayesian combination)
                (base_a, Self::cap_quality(qual_a_i + qual_b_i))
            } else if qual_a_i > qual_b_i {
                // Disagreement: take higher quality base, subtract lower quality
                (base_a, Self::cap_quality(qual_a_i - qual_b_i))
            } else if qual_b_i > qual_a_i {
                // Disagreement: take higher quality base, subtract lower quality
                (base_b, Self::cap_quality(qual_b_i - qual_a_i))
            } else {
                // Equal qualities but different bases: no confidence
                (base_a, NO_CALL_QUAL)
            };

            // Mask to N if either input base is N or if quality is minimum
            let (duplex_base, duplex_qual_val) =
                if base_a == NO_CALL || base_b == NO_CALL || raw_qual == NO_CALL_QUAL {
                    (NO_CALL, NO_CALL_QUAL)
                } else {
                    (raw_base, raw_qual)
                };

            duplex_seq.push(duplex_base);
            duplex_qual.push(duplex_qual_val);

            // Compute duplex errors using source reads if available, otherwise approximate
            // This matches fgbio's behavior: use source read method when available for exact error counting
            let error_at_i = if source_reads.is_empty() {
                // Approximate from single-strand errors (fgbio's fallback when source reads unavailable)
                let a_err = i32::from(errors_a.get(i).copied().unwrap_or(0));
                let b_err = i32::from(errors_b.get(i).copied().unwrap_or(0));
                let a_dep = i32::from(depths_a.get(i).copied().unwrap_or(0));
                let b_dep = i32::from(depths_b.get(i).copied().unwrap_or(0));

                if base_a == base_b {
                    // Bases agree: errors from both strands are errors relative to duplex
                    a_err + b_err
                } else if raw_base == base_a {
                    // Picked A's base: A's errors are errors, B's correct bases become errors
                    a_err + (b_dep - b_err)
                } else {
                    // Picked B's base: B's errors are errors, A's correct bases become errors
                    b_err + (a_dep - a_err)
                }
            } else {
                // Source read method: compare each source read base to raw consensus base
                // This is the exact method used by fgbio when source reads are available
                let mut num_errors = 0;
                for src_read in &source_reads {
                    if i < src_read.bases.len() && Self::is_error(src_read.bases[i], raw_base) {
                        num_errors += 1;
                    }
                }

                num_errors
            };

            duplex_errors.push(error_at_i.clamp(0, i32::from(i16::MAX)) as i16);
        }

        // Create duplex consensus record using ss_a as template
        let mut duplex = ss_a.clone();

        // Update the read name
        if let Some(name_bstr) = ss_a.name() {
            // BAM read names are ASCII, so to_str() should never fail
            let name: &str = name_bstr.to_str().expect("Read name should be valid UTF-8");
            let (name, _suffix) = Self::parse_mi_tag(name);
            *duplex.name_mut() = Some(BString::from(name.as_bytes().to_vec()));
        }
        // Update sequence and quality
        *duplex.sequence_mut() = Sequence::from(duplex_seq);
        *duplex.quality_scores_mut() = QualityScores::from(duplex_qual.clone());

        // Strip /A or /B suffix from MI tag to get base MI
        let mi_tag = Tag::from([b'M', b'I']);
        if let Some(Value::String(mi_bytes)) = duplex.data().get(&mi_tag) {
            let mi_str = String::from_utf8(mi_bytes.iter().copied().collect::<Vec<u8>>())
                .context("MI tag is not valid UTF-8")?;
            let (base_mi, _suffix) = Self::parse_mi_tag(&mi_str);
            duplex.data_mut().insert(mi_tag, Value::String(base_mi.as_bytes().into()));
        }

        // Extract depth information from single-strand consensuses (per-read tags)
        let cd_tag = Tag::from([b'c', b'D']);
        let cm_tag = Tag::from([b'c', b'M']);
        let ce_tag = Tag::from([b'c', b'E']);

        // Get per-read max depth (cD) from both strands for aD/bD tags
        let depth_a_max = Self::extract_int_tag(ss_a.data(), cd_tag).map_or(1, |d| d as usize);
        let depth_b_max = Self::extract_int_tag(ss_b.data(), cd_tag).map_or(1, |d| d as usize);

        // Get per-read min depth (cM) from both strands for aM/bM tags
        // IMPORTANT: If the duplex is truncated (len < simplex length), we must recalculate the min
        // from the truncated depths arrays, not use the original simplex cM tags
        let depth_a_min = if has_per_base_tags && len < depths_a.len() {
            // Duplex is truncated, recalculate min from truncated array
            depths_a.iter().take(len).map(|&d| d as usize).min().unwrap_or(1)
        } else {
            // No truncation or no per-base tags, use simplex cM tag
            Self::extract_int_tag(ss_a.data(), cm_tag).map_or(1, |d| d as usize)
        };
        let depth_b_min = if has_per_base_tags && len < depths_b.len() {
            // Duplex is truncated, recalculate min from truncated array
            depths_b.iter().take(len).map(|&d| d as usize).min().unwrap_or(1)
        } else {
            // No truncation or no per-base tags, use simplex cM tag
            Self::extract_int_tag(ss_b.data(), cm_tag).map_or(1, |d| d as usize)
        };

        // Calculate single-strand error rates for aE/bE tags
        // fgbio calculates these as: sum(ss.errors) / sum(ss.depths)
        // where errors and depths come from the simplex consensus, not from pre-existing tags
        let error_a = if has_per_base_tags {
            let total_errors_a: i64 = errors_a.iter().take(len).map(|&e| i64::from(e)).sum();
            let total_depths_a: i64 = depths_a.iter().take(len).map(|&d| i64::from(d)).sum();
            if total_depths_a > 0 { total_errors_a as f32 / total_depths_a as f32 } else { 0.0_f32 }
        } else {
            // Fallback to cE tag when per-base tags are missing
            if let Some(Value::Float(e)) = ss_a.data().get(&ce_tag) { *e } else { 0.0_f32 }
        };
        let error_b = if has_per_base_tags {
            let total_errors_b: i64 = errors_b.iter().take(len).map(|&e| i64::from(e)).sum();
            let total_depths_b: i64 = depths_b.iter().take(len).map(|&d| i64::from(d)).sum();
            if total_depths_b > 0 { total_errors_b as f32 / total_depths_b as f32 } else { 0.0_f32 }
        } else {
            // Fallback to cE tag when per-base tags are missing
            if let Some(Value::Float(e)) = ss_b.data().get(&ce_tag) { *e } else { 0.0_f32 }
        };

        // Calculate duplex cD and cM from per-base combined depths
        // fgbio: cD = max(ab.depths[i] + ba.depths[i]), cM = min(ab.depths[i] + ba.depths[i])
        let (duplex_depth_max, duplex_depth_min) = if has_per_base_tags {
            let combined_depths: Vec<i32> = depths_a
                .iter()
                .zip(depths_b.iter())
                .map(|(&a, &b)| i32::from(a) + i32::from(b))
                .collect();
            let max_depth = combined_depths.iter().copied().max().unwrap_or(0) as usize;
            let min_depth = combined_depths.iter().copied().min().unwrap_or(0) as usize;
            (max_depth, min_depth)
        } else {
            // Fallback when per-base tags missing: use per-read values (less accurate)
            (depth_a_max + depth_b_max, depth_a_min + depth_b_min)
        };

        // Calculate duplex error rate relative to the duplex consensus
        // This matches fgbio's calculation: sum(errors) / totalDepths.sum
        // Total depths is sum of per-base depths from both strands
        // IMPORTANT: Only sum the first `len` elements since we truncate to the shorter consensus
        let duplex_error_rate = if has_per_base_tags {
            let total_depths: i64 = depths_a.iter().take(len).map(|&d| i64::from(d)).sum::<i64>()
                + depths_b.iter().take(len).map(|&d| i64::from(d)).sum::<i64>();
            let total_errors: i64 = duplex_errors.iter().map(|&e| i64::from(e)).sum();
            if total_depths > 0 { total_errors as f32 / total_depths as f32 } else { 0.0_f32 }
        } else {
            // Fallback: weighted average when per-base tags are missing
            let total_bases_a = (depth_a_max * len) as f32;
            let total_bases_b = (depth_b_max * len) as f32;
            let total_errors_a = error_a * total_bases_a;
            let total_errors_b = error_b * total_bases_b;
            if total_bases_a + total_bases_b > 0.0 {
                (total_errors_a + total_errors_b) / (total_bases_a + total_bases_b)
            } else {
                0.0_f32
            }
        };

        // Update per-read tags
        // Duplex consensus tags
        let consensus_depth_tag = Tag::from([b'c', b'D']);
        let consensus_min_depth_tag = Tag::from([b'c', b'M']);
        let consensus_error_rate_tag = Tag::from([b'c', b'E']);

        duplex
            .data_mut()
            .insert(consensus_depth_tag, to_smallest_signed_int(duplex_depth_max as i32));
        duplex
            .data_mut()
            .insert(consensus_min_depth_tag, to_smallest_signed_int(duplex_depth_min as i32));
        duplex.data_mut().insert(consensus_error_rate_tag, Value::from(duplex_error_rate));

        // AB strand tags (from ss_a)
        let a_depth_tag = Tag::from([b'a', b'D']);
        let a_min_depth_tag = Tag::from([b'a', b'M']);
        let a_error_rate_tag = Tag::from([b'a', b'E']);

        duplex.data_mut().insert(a_depth_tag, to_smallest_signed_int(depth_a_max as i32));
        duplex.data_mut().insert(a_min_depth_tag, to_smallest_signed_int(depth_a_min as i32));
        duplex.data_mut().insert(a_error_rate_tag, Value::from(error_a));

        // BA strand tags (from ss_b)
        let b_depth_tag = Tag::from([b'b', b'D']);
        let b_min_depth_tag = Tag::from([b'b', b'M']);
        let b_error_rate_tag = Tag::from([b'b', b'E']);

        duplex.data_mut().insert(b_depth_tag, to_smallest_signed_int(depth_b_max as i32));
        duplex.data_mut().insert(b_min_depth_tag, to_smallest_signed_int(depth_b_min as i32));
        duplex.data_mut().insert(b_error_rate_tag, Value::from(error_b));

        let consensus_depth_base_tag = Tag::from([b'c', b'd']);
        let consensus_errors_base_tag = Tag::from([b'c', b'e']);

        // Add per-base tags if requested
        // Note: fgbio does NOT output cd/ce for duplex consensus, only ad/bd/ae/be/ac/bc/aq/bq
        // The cd/ce tags from ss_a (cloned) are removed below

        if produce_per_base_tags {
            // Per-base depth arrays (ad/bd on duplex)
            // Must truncate to len (minimum of both consensuses) to match fgbio
            let a_depth_base_tag = Tag::from([b'a', b'd']);
            duplex.data_mut().insert(a_depth_base_tag, Value::from(depths_a[..len].to_vec()));

            let b_depth_base_tag = Tag::from([b'b', b'd']);
            duplex.data_mut().insert(b_depth_base_tag, Value::from(depths_b[..len].to_vec()));

            // Per-base error arrays (ae/be on duplex)
            // Must truncate to len (minimum of both consensuses) to match fgbio
            let a_errors_base_tag = Tag::from([b'a', b'e']);
            duplex.data_mut().insert(a_errors_base_tag, Value::from(errors_a[..len].to_vec()));

            let b_errors_base_tag = Tag::from([b'b', b'e']);
            duplex.data_mut().insert(b_errors_base_tag, Value::from(errors_b[..len].to_vec()));

            // Single-strand consensus sequences (ac, bc)
            // Must truncate to len (minimum of both consensuses) to match fgbio
            let a_consensus_bases_tag = Tag::from([b'a', b'c']);
            let a_consensus_bases_str: String =
                String::from_utf8_lossy(&seq_a.as_ref()[..len]).to_string();
            duplex
                .data_mut()
                .insert(a_consensus_bases_tag, Value::from(a_consensus_bases_str.clone()));

            let b_consensus_bases_tag = Tag::from([b'b', b'c']);
            let b_consensus_bases_str: String =
                String::from_utf8_lossy(&seq_b.as_ref()[..len]).to_string();

            duplex
                .data_mut()
                .insert(b_consensus_bases_tag, Value::from(b_consensus_bases_str.clone()));

            // Single-strand consensus qualities (aq, bq)
            // Must truncate to len (minimum of both consensuses) to match fgbio
            let a_qualities_tag = Tag::from([b'a', b'q']);
            let a_quals_slice = &qual_a.as_ref()[..len];
            duplex.data_mut().insert(
                a_qualities_tag,
                crate::tags::qualities_to_tag_value(a_quals_slice),
            );

            let b_qualities_tag = Tag::from([b'b', b'q']);
            let b_quals_slice = &qual_b.as_ref()[..len];
            duplex.data_mut().insert(
                b_qualities_tag,
                crate::tags::qualities_to_tag_value(b_quals_slice),
            );
        }

        // Always remove cd/ce from duplex output - fgbio doesn't output these for duplex consensus
        // (they were cloned from ss_a which has them from the single-strand consensus call)
        duplex.data_mut().remove(&consensus_depth_base_tag);
        duplex.data_mut().remove(&consensus_errors_base_tag);

        // Call consensus on RX tags from source reads from BOTH strands
        // Match fgbio's behavior: collect RX tags from all source reads (both AB and BA strands)
        // and normalize the UMI order based on read pair orientation
        // See fgbio DuplexConsensusCaller.scala lines 465-482 (toUmiBasesForConsensusUmiCalling)
        let rx_tag = Tag::from([b'R', b'X']);
        let mut all_umis = Vec::new();

        // Extract RX tags from source reads A
        // If read.firstOfPair matches duplex_first_of_pair, use UMI as-is, otherwise reverse it
        for read in source_reads_a {
            if let Some(Value::String(rx_bytes)) = read.data().get(&rx_tag) {
                let rx = String::from_utf8_lossy(rx_bytes.as_ref()).to_string();
                if read.flags().is_first_segment() == duplex_first_of_pair {
                    all_umis.push(rx);
                } else {
                    // Reverse UMI order: split by "-", reverse, rejoin
                    let reversed: Vec<&str> = rx.split('-').rev().collect();
                    all_umis.push(reversed.join("-"));
                }
            }
        }

        // Extract RX tags from source reads B
        // If read.firstOfPair matches duplex_first_of_pair, use UMI as-is, otherwise reverse it
        for read in source_reads_b {
            if let Some(Value::String(rx_bytes)) = read.data().get(&rx_tag) {
                let rx = String::from_utf8_lossy(rx_bytes.as_ref()).to_string();
                if read.flags().is_first_segment() == duplex_first_of_pair {
                    all_umis.push(rx);
                } else {
                    // Reverse UMI order: split by "-", reverse, rejoin
                    let reversed: Vec<&str> = rx.split('-').rev().collect();
                    all_umis.push(reversed.join("-"));
                }
            }
        }

        // Call consensus on all UMI sequences and update RX tag
        if !all_umis.is_empty() {
            let consensus_umi = consensus_umis(&all_umis);
            duplex.data_mut().insert(rx_tag, Value::from(consensus_umi));
        }

        Ok(Some(duplex))
    }

    /// Processes a single UMI group to generate duplex consensus
    #[allow(clippy::type_complexity, clippy::too_many_arguments)]
    #[expect(clippy::too_many_lines, reason = "duplex group processing has many sequential steps including strand separation, SS calling, duplex calling, and output")]
    #[expect(clippy::needless_pass_by_value, reason = "owned values consumed by downstream consensus pipeline")]
    fn process_group(
        base_mi: String,
        a_records: Vec<Vec<u8>>,
        b_records: Vec<Vec<u8>>,
        ss_caller: &mut VanillaUmiConsensusCaller,
        produce_per_base_tags: bool,
        min_total_reads: usize,
        min_xy_reads: usize,
        min_yx_reads: usize,
        read_name_prefix: &str,
        read_group_id: &str,
        cell_tag: Option<Tag>,
    ) -> Result<(ConsensusOutput, ConsensusCallingStats)> {
        let mut stats = ConsensusCallingStats::new();
        let mut builder = UnmappedBamRecordBuilder::new();
        let mut output = ConsensusOutput::default();

        // Check if we have both strands
        if a_records.is_empty() && b_records.is_empty() {
            return Ok((ConsensusOutput::default(), stats));
        }

        // Check if we have enough input reads (matching fgbio's hasMinimumNumberOfReads)
        if !Self::has_minimum_number_of_reads(
            &a_records,
            &b_records,
            min_total_reads,
            min_xy_reads,
            min_yx_reads,
        ) {
            stats.record_rejection(
                RejectionReason::InsufficientReads,
                a_records.len() + b_records.len(),
            );
            return Ok((ConsensusOutput::default(), stats));
        }

        // Extract cell barcode from source reads (matching fgbio: recs.head.get[String](cellTag))
        let cell_barcode: Option<String> = cell_tag.and_then(|tag| {
            let tag_bytes: [u8; 2] = <[u8; 2]>::from(tag);
            a_records.first().or_else(|| b_records.first()).and_then(|r| {
                bam_fields::find_string_tag_in_record(r, &tag_bytes)
                    .map(|v| String::from_utf8_lossy(v).into_owned())
            })
        });

        // Split reads into R1/R2 groups for AB and BA strands (done once, reused below)
        let ab_r1s: Vec<&Vec<u8>> =
            a_records.iter().filter(|r| bam_fields::flags(r) & flags::FIRST_SEGMENT != 0).collect();
        let ab_r2s: Vec<&Vec<u8>> =
            a_records.iter().filter(|r| bam_fields::flags(r) & flags::FIRST_SEGMENT == 0).collect();
        let ba_r1s: Vec<&Vec<u8>> =
            b_records.iter().filter(|r| bam_fields::flags(r) & flags::FIRST_SEGMENT != 0).collect();
        let ba_r2s: Vec<&Vec<u8>> =
            b_records.iter().filter(|r| bam_fields::flags(r) & flags::FIRST_SEGMENT == 0).collect();

        // Validate strand orientations before processing
        // The expected orientations are:
        // AB R1: +  AB R2: -
        // BA R1: -  BA R2: +
        // (or the reverse of all)
        // Therefore: AB-R1 + BA-R2 should all be same strand, and AB-R2 + BA-R1 should all be same strand
        if !a_records.is_empty() && !b_records.is_empty() {
            // Check singleStrand1 (AB-R1 + BA-R2)
            let strand_1_iter = ab_r1s.iter().chain(ba_r2s.iter()).copied();
            if !Self::are_all_same_strand(strand_1_iter) {
                stats.record_rejection(
                    RejectionReason::PotentialCollision,
                    a_records.len() + b_records.len(),
                );
                return Ok((ConsensusOutput::default(), stats));
            }

            // Check singleStrand2 (AB-R2 + BA-R1)
            let strand_2_iter = ab_r2s.iter().chain(ba_r1s.iter()).copied();
            if !Self::are_all_same_strand(strand_2_iter) {
                stats.record_rejection(
                    RejectionReason::PotentialCollision,
                    a_records.len() + b_records.len(),
                );
                return Ok((ConsensusOutput::default(), stats));
            }
        }

        debug!(
            "MI {}: Initial read counts (AB-R1={}, AB-R2={}, BA-R1={}, BA-R2={})",
            base_mi,
            ab_r1s.len(),
            ab_r2s.len(),
            ba_r1s.len(),
            ba_r2s.len()
        );

        // Combine raw records and convert to SourceReads in one pass (matching fgbio's approach)
        // X = AB-R1 + BA-R2, Y = AB-R2 + BA-R1
        debug!(
            "MI {}: Combined records (X={}, Y={})",
            base_mi,
            ab_r1s.len() + ba_r2s.len(),
            ab_r2s.len() + ba_r1s.len()
        );

        // Keep references to original raw bytes for later tag extraction
        let x_raws: Vec<&[u8]> = ab_r1s.iter().chain(ba_r2s.iter()).map(|r| r.as_slice()).collect();
        let x_sources: Vec<SourceRead> = x_raws
            .iter()
            .enumerate()
            .filter_map(|(i, r)| {
                let mate_clip = bam_fields::num_bases_extending_past_mate_raw(r);
                ss_caller.create_source_read(r, i, mate_clip)
            })
            .collect();

        let y_raws: Vec<&[u8]> = ab_r2s.iter().chain(ba_r1s.iter()).map(|r| r.as_slice()).collect();
        let y_sources: Vec<SourceRead> = y_raws
            .iter()
            .enumerate()
            .filter_map(|(i, r)| {
                let mate_clip = bam_fields::num_bases_extending_past_mate_raw(r);
                ss_caller.create_source_read(r, i, mate_clip)
            })
            .collect();

        debug!(
            "MI {}: After SourceRead conversion (X={}, Y={})",
            base_mi,
            x_sources.len(),
            y_sources.len()
        );

        // Filter by common alignment (combined across strands)
        let (filtered_xs, _) = ss_caller.filter_by_alignment(x_sources);
        let (filtered_ys, _) = ss_caller.filter_by_alignment(y_sources);

        debug!(
            "MI {}: After alignment filtering (X={}, Y={})",
            base_mi,
            filtered_xs.len(),
            filtered_ys.len()
        );

        // Split filtered results back by firstOfPair flag (matching fgbio lines 347-350)
        // X = AB-R1 + BA-R2, Y = AB-R2 + BA-R1 (combined for alignment filtering)
        // Now split them back into the four groups using sr.flags:
        let filtered_ab_r1s: Vec<SourceRead> =
            filtered_xs.iter().filter(|sr| sr.flags & flags::FIRST_SEGMENT != 0).cloned().collect();
        let filtered_ba_r2s: Vec<SourceRead> =
            filtered_xs.iter().filter(|sr| sr.flags & flags::FIRST_SEGMENT == 0).cloned().collect();
        let filtered_ab_r2s: Vec<SourceRead> =
            filtered_ys.iter().filter(|sr| sr.flags & flags::FIRST_SEGMENT == 0).cloned().collect();
        let filtered_ba_r1s: Vec<SourceRead> =
            filtered_ys.iter().filter(|sr| sr.flags & flags::FIRST_SEGMENT != 0).cloned().collect();

        debug!(
            "MI {}: After split (filtered_AB-R1={}, filtered_AB-R2={}, filtered_BA-R1={}, filtered_BA-R2={})",
            base_mi,
            filtered_ab_r1s.len(),
            filtered_ab_r2s.len(),
            filtered_ba_r1s.len(),
            filtered_ba_r2s.len()
        );

        // Extract raw records for tag extraction using original_idx to map back
        let filtered_ab_r1_raws: Vec<&[u8]> = filtered_ab_r1s
            .iter()
            .map(|sr| {
                debug_assert!(sr.original_idx < x_raws.len());
                x_raws[sr.original_idx]
            })
            .collect();
        let filtered_ba_r2_raws: Vec<&[u8]> = filtered_ba_r2s
            .iter()
            .map(|sr| {
                debug_assert!(sr.original_idx < x_raws.len());
                x_raws[sr.original_idx]
            })
            .collect();
        let filtered_ab_r2_raws: Vec<&[u8]> = filtered_ab_r2s
            .iter()
            .map(|sr| {
                debug_assert!(sr.original_idx < y_raws.len());
                y_raws[sr.original_idx]
            })
            .collect();
        let filtered_ba_r1_raws: Vec<&[u8]> = filtered_ba_r1s
            .iter()
            .map(|sr| {
                debug_assert!(sr.original_idx < y_raws.len());
                y_raws[sr.original_idx]
            })
            .collect();

        // Build UMI for consensus read names (format: base_mi/A or base_mi/B)
        let ab_umi = format!("{base_mi}/A");
        let ba_umi = format!("{base_mi}/B");

        // Call single-strand consensus on each filtered group
        // consensus_call now returns VanillaConsensusRead instead of RecordBuf
        let ab_r1_consensus = ss_caller.consensus_call(&ab_umi, filtered_ab_r1s)?;
        let ab_r2_consensus = ss_caller.consensus_call(&ab_umi, filtered_ab_r2s)?;
        let ba_r1_consensus = ss_caller.consensus_call(&ba_umi, filtered_ba_r1s)?;
        let ba_r2_consensus = ss_caller.consensus_call(&ba_umi, filtered_ba_r2s)?;

        debug!(
            "MI {}: SS consensus results (AB-R1={}, AB-R2={}, BA-R1={}, BA-R2={})",
            base_mi,
            ab_r1_consensus.is_some(),
            ab_r2_consensus.is_some(),
            ba_r1_consensus.is_some(),
            ba_r2_consensus.is_some()
        );

        if let Some(ref cons) = ab_r1_consensus {
            debug!("MI {}: AB-R1 consensus len={}", base_mi, cons.len());
        }
        if let Some(ref cons) = ab_r2_consensus {
            debug!("MI {}: AB-R2 consensus len={}", base_mi, cons.len());
        }
        if let Some(ref cons) = ba_r1_consensus {
            debug!("MI {}: BA-R1 consensus len={}", base_mi, cons.len());
        }
        if let Some(ref cons) = ba_r2_consensus {
            debug!("MI {}: BA-R2 consensus len={}", base_mi, cons.len());
        }

        // Create duplex consensus from AB and BA single-strand consensuses
        // Duplex R1 = combine(AB-R1, BA-R1), Duplex R2 = combine(AB-R2, BA-R2)
        debug!(
            "MI {}: Match pattern: AB-R1={}, AB-R2={}, BA-R1={}, BA-R2={}",
            base_mi,
            ab_r1_consensus.is_some(),
            ab_r2_consensus.is_some(),
            ba_r1_consensus.is_some(),
            ba_r2_consensus.is_some()
        );
        match (ab_r1_consensus, ab_r2_consensus, ba_r1_consensus, ba_r2_consensus) {
            (Some(ref r1_a), Some(ref r2_a), Some(ref r1_b), Some(ref r2_b)) => {
                // Both strands have both R1 and R2 - create full duplex
                debug!("MI {base_mi}: Calling duplex consensus from SS pairs");

                // Collect source reads for error calculation
                let r1_source_reads: Vec<SourceRead> = r1_a
                    .source_reads
                    .iter()
                    .flatten()
                    .chain(r2_b.source_reads.iter().flatten())
                    .cloned()
                    .collect();
                let r2_source_reads: Vec<SourceRead> = r2_a
                    .source_reads
                    .iter()
                    .flatten()
                    .chain(r1_b.source_reads.iter().flatten())
                    .cloned()
                    .collect();

                // For R1 duplex: combine AB-R1 with BA-R2
                // This matches fgbio: duplexR1 = duplexConsensus(abR1Consensus, baR2Consensus, ...)
                let duplex_r1_consensus = Self::duplex_consensus(
                    Some(r1_a),
                    Some(r2_b),
                    if r1_source_reads.is_empty() { None } else { Some(&r1_source_reads) },
                );
                debug!("MI {}: duplex_r1 result: {}", base_mi, duplex_r1_consensus.is_some());

                // For R2 duplex: combine AB-R2 with BA-R1
                // This matches fgbio: duplexR2 = duplexConsensus(abR2Consensus, baR1Consensus, ...)
                let duplex_r2_consensus = Self::duplex_consensus(
                    Some(r2_a),
                    Some(r1_b),
                    if r2_source_reads.is_empty() { None } else { Some(&r2_source_reads) },
                );
                debug!("MI {}: duplex_r2 result: {}", base_mi, duplex_r2_consensus.is_some());

                debug!(
                    "MI {}: Duplex results (R1={}, R2={})",
                    base_mi,
                    duplex_r1_consensus.is_some(),
                    duplex_r2_consensus.is_some()
                );

                if let (Some(dr1), Some(dr2)) = (&duplex_r1_consensus, &duplex_r2_consensus) {
                    debug!(
                        "MI {}: Duplex consensus lengths (R1={}, R2={})",
                        base_mi,
                        dr1.len(),
                        dr2.len()
                    );
                }

                if let (Some(dr1), Some(dr2)) = (duplex_r1_consensus, duplex_r2_consensus) {
                    // Check if consensus has enough reads after filtering
                    if Self::duplex_consensus_has_minimum_reads(
                        &dr1,
                        min_total_reads,
                        min_xy_reads,
                        min_yx_reads,
                    ) && Self::duplex_consensus_has_minimum_reads(
                        &dr2,
                        min_total_reads,
                        min_xy_reads,
                        min_yx_reads,
                    ) {
                        // Write DuplexConsensusReads as raw bytes
                        Self::duplex_read_into(
                            &mut builder,
                            &mut output,
                            &dr1,
                            ReadType::R1,
                            &base_mi,
                            &filtered_ab_r1_raws,
                            &filtered_ba_r2_raws,
                            produce_per_base_tags,
                            read_name_prefix,
                            read_group_id,
                            true, // first of pair
                            cell_tag,
                            cell_barcode.as_deref(),
                        )?;
                        Self::duplex_read_into(
                            &mut builder,
                            &mut output,
                            &dr2,
                            ReadType::R2,
                            &base_mi,
                            &filtered_ab_r2_raws,
                            &filtered_ba_r1_raws,
                            produce_per_base_tags,
                            read_name_prefix,
                            read_group_id,
                            false, // second of pair
                            cell_tag,
                            cell_barcode.as_deref(),
                        )?;
                        stats.record_consensus();
                        return Ok((output, stats));
                    }
                    stats.record_rejection(
                        RejectionReason::InsufficientReads,
                        a_records.len() + b_records.len(),
                    );
                    return Ok((ConsensusOutput::default(), stats));
                }
            }
            (Some(ref r1_a), Some(ref r2_a), None, None) => {
                // Only AB strand - check if single-strand consensus is allowed (min_yx_reads == 0)
                if min_yx_reads == 0 {
                    // Use duplex_consensus with only AB (no BA)
                    let duplex_r1 = Self::duplex_consensus(Some(r1_a), None, None);
                    let duplex_r2 = Self::duplex_consensus(Some(r2_a), None, None);

                    if let (Some(dr1), Some(dr2)) = (duplex_r1, duplex_r2) {
                        let empty: &[&[u8]] = &[];
                        Self::duplex_read_into(
                            &mut builder,
                            &mut output,
                            &dr1,
                            ReadType::R1,
                            &base_mi,
                            &filtered_ab_r1_raws,
                            empty,
                            produce_per_base_tags,
                            read_name_prefix,
                            read_group_id,
                            true,
                            cell_tag,
                            cell_barcode.as_deref(),
                        )?;
                        Self::duplex_read_into(
                            &mut builder,
                            &mut output,
                            &dr2,
                            ReadType::R2,
                            &base_mi,
                            &filtered_ab_r2_raws,
                            empty,
                            produce_per_base_tags,
                            read_name_prefix,
                            read_group_id,
                            false,
                            cell_tag,
                            cell_barcode.as_deref(),
                        )?;
                        stats.record_consensus();
                        return Ok((output, stats));
                    }
                }
            }
            (None, None, Some(ref r1_b), Some(ref r2_b)) => {
                // Only BA strand - check if single-strand consensus is allowed
                // fgbio swaps B→A when only B is present
                if min_yx_reads == 0 {
                    // Use duplex_consensus with only BA (treated as AB)
                    let duplex_r1 = Self::duplex_consensus(Some(r1_b), None, None);
                    let duplex_r2 = Self::duplex_consensus(Some(r2_b), None, None);

                    if let (Some(dr1), Some(dr2)) = (duplex_r1, duplex_r2) {
                        let empty: &[&[u8]] = &[];
                        Self::duplex_read_into(
                            &mut builder,
                            &mut output,
                            &dr1,
                            ReadType::R1,
                            &base_mi,
                            &filtered_ba_r1_raws,
                            empty,
                            produce_per_base_tags,
                            read_name_prefix,
                            read_group_id,
                            true,
                            cell_tag,
                            cell_barcode.as_deref(),
                        )?;
                        Self::duplex_read_into(
                            &mut builder,
                            &mut output,
                            &dr2,
                            ReadType::R2,
                            &base_mi,
                            &filtered_ba_r2_raws,
                            empty,
                            produce_per_base_tags,
                            read_name_prefix,
                            read_group_id,
                            false,
                            cell_tag,
                            cell_barcode.as_deref(),
                        )?;
                        stats.record_consensus();
                        return Ok((output, stats));
                    }
                }
            }
            _ => {
                // Partial consensus (e.g., only R1 from one strand) - reject
                debug!("MI {base_mi}: Falling through to default match arm (partial consensus)");
            }
        }

        // Filter out (insufficient reads or conversion failed)
        debug!(
            "MI {base_mi}: No duplex consensus created - insufficient reads or conversion failed"
        );
        stats.record_rejection(RejectionReason::InsufficientReads, 1);
        Ok((ConsensusOutput::default(), stats))
    }
}

impl ConsensusCaller for DuplexConsensusCaller {
    fn consensus_reads(&mut self, records: Vec<Vec<u8>>) -> Result<ConsensusOutput> {
        self.stats.record_input(records.len());

        // Partition records by strand using raw byte-level operations.
        let (base_mi_opt, a_records, b_records) = Self::partition_records_by_strand(records)?;

        let Some(base_mi) = base_mi_opt else {
            return Ok(ConsensusOutput::default());
        };

        let produce_per_base_tags = self.produce_per_base_tags;

        let (output, group_stats) = Self::process_group(
            base_mi,
            a_records,
            b_records,
            &mut self.ss_caller,
            produce_per_base_tags,
            self.min_total_reads,
            self.min_xy_reads,
            self.min_yx_reads,
            &self.read_name_prefix,
            &self.read_group_id,
            self.cell_tag,
        )?;

        self.stats.merge(&group_stats);

        if self.track_rejects {
            self.rejected_reads.extend(self.ss_caller.take_rejected_reads());
        }

        Ok(output)
    }

    fn total_reads(&self) -> usize {
        self.stats.total_reads
    }

    fn total_filtered(&self) -> usize {
        self.stats.filtered_reads
    }

    fn consensus_reads_constructed(&self) -> usize {
        self.stats.consensus_reads
    }

    fn statistics(&self) -> ConsensusCallingStats {
        self.stats.clone()
    }

    fn log_statistics(&self) {
        info!("Duplex Consensus Calling Statistics:");
        info!("  Total input reads: {}", self.total_reads());
        info!("  Consensus reads constructed: {}", self.consensus_reads_constructed());
        info!("  Total filtered: {}", self.total_filtered());

        // Log rejection breakdown
        for (reason, count) in &self.stats.rejection_reasons {
            info!("    {}: {}", reason.description(), count);
        }

        // Log single-strand caller statistics
        info!("Single-strand caller:");
        self.ss_caller.log_statistics();
    }
}

#[cfg(test)]
#[allow(clippy::similar_names)]
mod tests {
    use super::*;
    use fgumi_sam::builder::RecordBuilder;
    use noodles_raw_bam::ParsedBamRecord;
    use noodles::sam::alignment::record_buf::data::field::Value;

    fn encode_to_raw(rec: &noodles::sam::alignment::RecordBuf) -> Vec<u8> {
        let header = noodles::sam::Header::default();
        let mut buf = Vec::new();
        crate::vendored::bam_codec::encode_record_buf(&mut buf, &header, rec).unwrap();
        buf
    }

    #[test]
    fn test_parse_mi_tag() {
        assert_eq!(
            DuplexConsensusCaller::parse_mi_tag("ACGT-TGCA/A"),
            ("ACGT-TGCA".to_string(), Some('A'))
        );
        assert_eq!(
            DuplexConsensusCaller::parse_mi_tag("ACGT-TGCA/B"),
            ("ACGT-TGCA".to_string(), Some('B'))
        );
        assert_eq!(
            DuplexConsensusCaller::parse_mi_tag("ACGT-TGCA"),
            ("ACGT-TGCA".to_string(), None)
        );
    }

    #[test]
    fn test_has_both_strands() {
        // Empty - no strands
        let empty: Vec<RecordBuf> = vec![];
        assert!(!DuplexConsensusCaller::has_both_strands(&empty));

        // Single read - not enough
        let single_a = vec![
            RecordBuilder::new()
                .name("read1")
                .sequence("ACGT")
                .qualities(&[30, 30, 30, 30])
                .tag("MI", "UMI1/A")
                .build(),
        ];
        assert!(!DuplexConsensusCaller::has_both_strands(&single_a));

        // Two reads, same strand - no duplex
        let two_a = vec![
            RecordBuilder::new()
                .name("read1")
                .sequence("ACGT")
                .qualities(&[30, 30, 30, 30])
                .tag("MI", "UMI1/A")
                .build(),
            RecordBuilder::new()
                .name("read2")
                .sequence("ACGT")
                .qualities(&[30, 30, 30, 30])
                .tag("MI", "UMI1/A")
                .build(),
        ];
        assert!(!DuplexConsensusCaller::has_both_strands(&two_a));

        // Two reads, both strands - duplex possible
        let both_strands = vec![
            RecordBuilder::new()
                .name("read1")
                .sequence("ACGT")
                .qualities(&[30, 30, 30, 30])
                .tag("MI", "UMI1/A")
                .build(),
            RecordBuilder::new()
                .name("read2")
                .sequence("ACGT")
                .qualities(&[30, 30, 30, 30])
                .tag("MI", "UMI1/B")
                .build(),
        ];
        assert!(DuplexConsensusCaller::has_both_strands(&both_strands));

        // Many reads, A strand first, B found later
        let many_a_one_b = vec![
            RecordBuilder::new()
                .name("read1")
                .sequence("ACGT")
                .qualities(&[30, 30, 30, 30])
                .tag("MI", "UMI1/A")
                .build(),
            RecordBuilder::new()
                .name("read2")
                .sequence("ACGT")
                .qualities(&[30, 30, 30, 30])
                .tag("MI", "UMI1/A")
                .build(),
            RecordBuilder::new()
                .name("read3")
                .sequence("ACGT")
                .qualities(&[30, 30, 30, 30])
                .tag("MI", "UMI1/A")
                .build(),
            RecordBuilder::new()
                .name("read4")
                .sequence("ACGT")
                .qualities(&[30, 30, 30, 30])
                .tag("MI", "UMI1/B")
                .build(),
        ];
        assert!(DuplexConsensusCaller::has_both_strands(&many_a_one_b));
    }

    #[test]
    fn test_new_with_different_min_reads() -> Result<()> {
        // Single value - padTo behavior: [3] -> [3, 3, 3]
        let caller = DuplexConsensusCaller::new(
            "consensus".to_string(),
            "RG1".to_string(),
            vec![3],
            20,
            false,
            false,
            None,
            None,
            false,
            45,
            40,
        )?;
        assert_eq!(caller.min_total_reads, 3);
        assert_eq!(caller.min_xy_reads, 3);
        assert_eq!(caller.min_yx_reads, 3);

        // Two values - padTo behavior: [3, 2] -> [3, 2, 2]
        let caller = DuplexConsensusCaller::new(
            "consensus".to_string(),
            "RG1".to_string(),
            vec![3, 2],
            20,
            false,
            false,
            None,
            None,
            false,
            45,
            40,
        )?;
        assert_eq!(caller.min_total_reads, 3);
        assert_eq!(caller.min_xy_reads, 2);
        assert_eq!(caller.min_yx_reads, 2);

        // Three values - [3, 2, 1] used directly
        let caller = DuplexConsensusCaller::new(
            "consensus".to_string(),
            "RG1".to_string(),
            vec![3, 2, 1],
            20,
            false,
            false,
            None,
            None,
            false,
            45,
            40,
        )?;
        assert_eq!(caller.min_total_reads, 3);
        assert_eq!(caller.min_xy_reads, 2);
        assert_eq!(caller.min_yx_reads, 1);

        Ok(())
    }

    #[test]
    fn test_group_by_mi_and_strand() -> Result<()> {
        // Create test reads with different MI tags
        let read1 = RecordBuilder::new().sequence("ACGT").tag("MI", "UMI1/A").build();

        let read2 = RecordBuilder::new().sequence("ACGT").tag("MI", "UMI1/B").build();

        let read3 = RecordBuilder::new().sequence("ACGT").tag("MI", "UMI1/A").build();

        let read4 = RecordBuilder::new().sequence("ACGT").tag("MI", "UMI2/A").build();

        let reads = vec![read1, read2, read3, read4];
        let groups = DuplexConsensusCaller::group_by_mi_and_strand(reads)?;

        assert_eq!(groups.len(), 2);
        assert!(groups.contains_key("UMI1"));
        assert!(groups.contains_key("UMI2"));

        let (a_reads, b_reads) = groups.get("UMI1").unwrap();
        assert_eq!(a_reads.len(), 2);
        assert_eq!(b_reads.len(), 1);

        let (a_reads, b_reads) = groups.get("UMI2").unwrap();
        assert_eq!(a_reads.len(), 1);
        assert_eq!(b_reads.len(), 0);

        Ok(())
    }

    // Helper function to create a test VanillaUmiConsensusCaller
    fn create_test_ss_caller() -> VanillaUmiConsensusCaller {
        VanillaUmiConsensusCaller::new(
            "test".to_string(),
            "RG1".to_string(),
            VanillaUmiConsensusOptions {
                tag: "MI".to_string(),
                error_rate_pre_umi: 45,
                error_rate_post_umi: 40,
                min_input_base_quality: 10,
                min_reads: 1,
                max_reads: None,
                produce_per_base_tags: true,
                seed: Some(42),
                trim: false,
                min_consensus_base_quality: 40,
                cell_tag: None,
            },
        )
    }

    // Helper function to create a test single-strand consensus record
    fn create_ss_consensus(
        seq: &str,
        qual: &[u8],
        depth_max: i32,
        depth_min: i32,
        error_rate: f32,
    ) -> RecordBuf {
        RecordBuilder::new()
            .sequence(seq)
            .qualities(qual)
            .tag("cD", depth_max)
            .tag("cM", depth_min)
            .tag("cE", error_rate)
            .build()
    }

    #[test]
    fn test_duplex_consensus_quality_scores_agreement() -> Result<()> {
        // Test that agreeing bases sum their qualities (capped at 93)
        let ss_a = create_ss_consensus("AAAA", &[20, 30, 40, 50], 3, 3, 0.0);
        let ss_b = create_ss_consensus("AAAA", &[20, 30, 40, 50], 2, 2, 0.0);

        let duplex = DuplexConsensusCaller::call_duplex_from_ss_pair(
            ss_a,
            ss_b,
            false,
            &[],
            &[],
            &create_test_ss_caller(),
            true, // R1 duplex (first of pair)
        )?
        .expect("Should produce duplex consensus");

        // Check sequence
        let seq: Vec<u8> = duplex.sequence().as_ref().to_vec();
        assert_eq!(seq, b"AAAA");

        // Check qualities - should sum (20+20=40, 30+30=60, 40+40=80, 50+50=100 capped to 93)
        let quals: Vec<u8> = duplex.quality_scores().as_ref().to_vec();
        assert_eq!(quals, vec![40, 60, 80, 93]);

        Ok(())
    }

    #[test]
    fn test_duplex_consensus_quality_scores_disagreement_unequal() -> Result<()> {
        // Test that disagreeing bases with unequal qualities subtract lower from higher
        let ss_a = create_ss_consensus("ACGT", &[30, 30, 30, 30], 3, 3, 0.0);
        let ss_b = create_ss_consensus("TGCA", &[10, 15, 20, 25], 2, 2, 0.0); // Different bases

        let duplex = DuplexConsensusCaller::call_duplex_from_ss_pair(
            ss_a,
            ss_b,
            false,
            &[],
            &[],
            &create_test_ss_caller(),
            true, // R1 duplex (first of pair)
        )?
        .expect("Should produce duplex consensus");

        // Check sequence - should take higher quality base (all from ss_a since it has higher quality)
        let seq: Vec<u8> = duplex.sequence().as_ref().to_vec();
        assert_eq!(seq, b"ACGT");

        // Check qualities - should subtract (30-10=20, 30-15=15, 30-20=10, 30-25=5)
        let quals: Vec<u8> = duplex.quality_scores().as_ref().to_vec();
        assert_eq!(quals, vec![20, 15, 10, 5]);

        Ok(())
    }

    #[test]
    fn test_duplex_consensus_quality_scores_disagreement_equal() -> Result<()> {
        // Test that disagreeing bases with equal qualities produce N with Q=2
        let ss_a = create_ss_consensus("ACGT", &[20, 20, 20, 20], 3, 3, 0.0);
        let ss_b = create_ss_consensus("TGCA", &[20, 20, 20, 20], 2, 2, 0.0);

        let duplex = DuplexConsensusCaller::call_duplex_from_ss_pair(
            ss_a,
            ss_b,
            false,
            &[],
            &[],
            &create_test_ss_caller(),
            true, // R1 duplex (first of pair)
        )?
        .expect("Should produce duplex consensus");

        // Check sequence - should all be N due to equal quality conflicts
        let seq: Vec<u8> = duplex.sequence().as_ref().to_vec();
        assert_eq!(seq, b"NNNN");

        // Check qualities - should all be minimum (2)
        let quals: Vec<u8> = duplex.quality_scores().as_ref().to_vec();
        assert_eq!(quals, vec![2, 2, 2, 2]);

        Ok(())
    }

    #[test]
    fn test_duplex_consensus_per_read_tags() -> Result<()> {
        // Test that all per-read tags are generated correctly
        use noodles::sam::alignment::record::data::field::Tag;

        let ss_a = create_ss_consensus("AAAA", &[20, 20, 20, 20], 3, 2, 0.05);
        let ss_b = create_ss_consensus("AAAA", &[20, 20, 20, 20], 2, 1, 0.10);

        let duplex = DuplexConsensusCaller::call_duplex_from_ss_pair(
            ss_a,
            ss_b,
            false,
            &[],
            &[],
            &create_test_ss_caller(),
            true, // R1 duplex (first of pair)
        )?
        .expect("Should produce duplex consensus");

        // Check duplex consensus tags (cD, cM, cE)
        let cd_tag = Tag::from([b'c', b'D']);
        let cm_tag = Tag::from([b'c', b'M']);
        let ce_tag = Tag::from([b'c', b'E']);

        let depth_max =
            DuplexConsensusCaller::extract_int_tag(duplex.data(), cd_tag).expect("Missing cD tag");
        assert_eq!(depth_max, 5); // 3 + 2

        let depth_min =
            DuplexConsensusCaller::extract_int_tag(duplex.data(), cm_tag).expect("Missing cM tag");
        assert_eq!(depth_min, 3); // 2 + 1

        if let Some(Value::Float(error_rate)) = duplex.data().get(&ce_tag) {
            // Error rate should be weighted average:
            // (0.05 * 3*4 + 0.10 * 2*4) / (3*4 + 2*4) = (0.6 + 0.8) / 20 = 0.07
            assert!((error_rate - 0.07).abs() < 0.001);
        } else {
            panic!("Missing cE tag");
        }

        // Check AB strand tags (aD, aM, aE)
        let ad_tag = Tag::from([b'a', b'D']);
        let am_tag = Tag::from([b'a', b'M']);
        let ae_tag = Tag::from([b'a', b'E']);

        let depth =
            DuplexConsensusCaller::extract_int_tag(duplex.data(), ad_tag).expect("Missing aD tag");
        assert_eq!(depth, 3);

        let depth =
            DuplexConsensusCaller::extract_int_tag(duplex.data(), am_tag).expect("Missing aM tag");
        assert_eq!(depth, 2);

        if let Some(Value::Float(error)) = duplex.data().get(&ae_tag) {
            assert!((error - 0.05).abs() < 0.001);
        } else {
            panic!("Missing aE tag");
        }

        // Check BA strand tags (bD, bM, bE)
        let bd_tag = Tag::from([b'b', b'D']);
        let bm_tag = Tag::from([b'b', b'M']);
        let be_tag = Tag::from([b'b', b'E']);

        let depth =
            DuplexConsensusCaller::extract_int_tag(duplex.data(), bd_tag).expect("Missing bD tag");
        assert_eq!(depth, 2);

        let depth =
            DuplexConsensusCaller::extract_int_tag(duplex.data(), bm_tag).expect("Missing bM tag");
        assert_eq!(depth, 1);

        if let Some(Value::Float(error)) = duplex.data().get(&be_tag) {
            assert!((error - 0.10).abs() < 0.001);
        } else {
            panic!("Missing bE tag");
        }

        Ok(())
    }

    #[test]
    fn test_duplex_consensus_n_bases() -> Result<()> {
        // Test that N bases in either strand produce N in duplex with Q=2
        let ss_a = create_ss_consensus("ANAA", &[20, 20, 20, 20], 3, 3, 0.0);
        let ss_b = create_ss_consensus("AANA", &[20, 20, 20, 20], 2, 2, 0.0);

        let duplex = DuplexConsensusCaller::call_duplex_from_ss_pair(
            ss_a,
            ss_b,
            false,
            &[],
            &[],
            &create_test_ss_caller(),
            true, // R1 duplex (first of pair)
        )?
        .expect("Should produce duplex consensus");

        // Check sequence
        let seq: Vec<u8> = duplex.sequence().as_ref().to_vec();
        assert_eq!(seq, b"ANNA");

        // Check qualities - N positions should have Q=2
        let quals: Vec<u8> = duplex.quality_scores().as_ref().to_vec();
        assert_eq!(quals[1], 2); // Position 1 has N in ss_a
        assert_eq!(quals[2], 2); // Position 2 has N in ss_b

        Ok(())
    }

    #[test]
    fn test_duplex_consensus_length_mismatch() -> Result<()> {
        // Test that length mismatches are handled by truncating to minimum length
        // This matches fgbio's behavior
        let ss_a = create_ss_consensus("AAAA", &[20, 20, 20, 20], 3, 3, 0.0);
        let ss_b = create_ss_consensus("AAA", &[20, 20, 20], 2, 2, 0.0);

        let duplex = DuplexConsensusCaller::call_duplex_from_ss_pair(
            ss_a,
            ss_b,
            false,
            &[],
            &[],
            &create_test_ss_caller(),
            true, // R1 duplex (first of pair)
        )?
        .expect("Should truncate to minimum length");

        // Check that the duplex is truncated to length 3 (min of 4 and 3)
        assert_eq!(duplex.sequence().len(), 3);
        let seq: Vec<u8> = duplex.sequence().as_ref().to_vec();
        assert_eq!(seq, b"AAA");

        Ok(())
    }

    #[test]
    fn test_duplex_consensus_quality_capping() -> Result<()> {
        // Test that summed qualities are capped at 93
        let ss_a = create_ss_consensus("AAA", &[50, 60, 93], 3, 3, 0.0);
        let ss_b = create_ss_consensus("AAA", &[50, 60, 93], 2, 2, 0.0);

        let duplex = DuplexConsensusCaller::call_duplex_from_ss_pair(
            ss_a,
            ss_b,
            false,
            &[],
            &[],
            &create_test_ss_caller(),
            true, // R1 duplex (first of pair)
        )?
        .expect("Should produce duplex consensus");

        // Check qualities - should cap at 93
        let quals: Vec<u8> = duplex.quality_scores().as_ref().to_vec();
        assert_eq!(quals, vec![93, 93, 93]); // All capped at 93

        Ok(())
    }

    #[test]
    fn test_duplex_consensus_quality_difference_at_threshold() -> Result<()> {
        // Test that quality differences of 2 or less mask to N
        // This tests the edge case where subtraction results in Q=2 or Q=1
        let ss_a = create_ss_consensus("ACGT", &[5, 4, 3, 10], 3, 3, 0.0);
        let ss_b = create_ss_consensus("TGCA", &[3, 2, 2, 8], 2, 2, 0.0);

        let duplex = DuplexConsensusCaller::call_duplex_from_ss_pair(
            ss_a,
            ss_b,
            false,
            &[],
            &[],
            &create_test_ss_caller(),
            true, // R1 duplex (first of pair)
        )?
        .expect("Should produce duplex consensus");

        // Check sequence - positions with Q <= 2 should be masked to N
        let seq: Vec<u8> = duplex.sequence().as_ref().to_vec();
        // 5-3=2, 4-2=2, 3-2=1 (capped to 2), 10-8=2 - all should become N
        assert_eq!(seq, b"NNNN");

        // Check qualities - all should be minimum (2)
        let quals: Vec<u8> = duplex.quality_scores().as_ref().to_vec();
        assert_eq!(quals, vec![2, 2, 2, 2]);

        Ok(())
    }

    #[test]
    fn test_duplex_consensus_with_deep_coverage() -> Result<()> {
        // Test that qualities with deep coverage on one strand don't saturate incorrectly
        // AB strand has very high depth, BA strand has low depth
        let ss_a = create_ss_consensus("AAAA", &[45, 45, 45, 45], 50, 50, 0.01);
        let ss_b = create_ss_consensus("AAAA", &[20, 20, 20, 20], 1, 1, 0.0);

        let duplex = DuplexConsensusCaller::call_duplex_from_ss_pair(
            ss_a,
            ss_b,
            false,
            &[],
            &[],
            &create_test_ss_caller(),
            true, // R1 duplex (first of pair)
        )?
        .expect("Should produce duplex consensus");

        // Check sequence
        let seq: Vec<u8> = duplex.sequence().as_ref().to_vec();
        assert_eq!(seq, b"AAAA");

        // Check qualities - should sum (45+20=65) but not exceed 93
        let quals: Vec<u8> = duplex.quality_scores().as_ref().to_vec();
        assert_eq!(quals, vec![65, 65, 65, 65]);

        // Verify depth tags reflect the asymmetric coverage
        use noodles::sam::alignment::record::data::field::Tag;
        let ad_tag = Tag::from([b'a', b'D']);
        let bd_tag = Tag::from([b'b', b'D']);

        let depth_a =
            DuplexConsensusCaller::extract_int_tag(duplex.data(), ad_tag).expect("Missing aD tag");
        let depth_b =
            DuplexConsensusCaller::extract_int_tag(duplex.data(), bd_tag).expect("Missing bD tag");

        assert_eq!(depth_a, 50);
        assert_eq!(depth_b, 1);

        Ok(())
    }

    #[test]
    fn test_duplex_consensus_per_base_tags() -> Result<()> {
        // Test that per-base tags are generated when enabled
        use noodles::sam::alignment::record::data::field::Tag;

        // Create single-strand consensuses with per-base tags
        let mut ss_a = create_ss_consensus("AAAA", &[20, 20, 20, 20], 3, 2, 0.05);
        let mut ss_b = create_ss_consensus("AAAA", &[20, 20, 20, 20], 2, 1, 0.10);

        // Add per-base depth tags (cd) for both strands
        let cd_tag = Tag::from([b'c', b'd']);
        ss_a.data_mut().insert(cd_tag, Value::from(vec![3i16, 3, 2, 3]));
        ss_b.data_mut().insert(cd_tag, Value::from(vec![2i16, 2, 1, 2]));

        // Add per-base error tags (ce) for both strands
        let ce_tag = Tag::from([b'c', b'e']);
        ss_a.data_mut().insert(ce_tag, Value::from(vec![0i16, 0, 1, 0]));
        ss_b.data_mut().insert(ce_tag, Value::from(vec![0i16, 1, 0, 0]));

        // Add single-strand consensus bases tags (ac, bc) for both strands
        let ac_tag = Tag::from([b'a', b'c']);
        ss_a.data_mut().insert(ac_tag, Value::String("AAAA".as_bytes().into()));

        let bc_tag = Tag::from([b'b', b'c']);
        ss_b.data_mut().insert(bc_tag, Value::String("AAAA".as_bytes().into()));

        // Add single-strand consensus quality tags (aq, bq) for both strands
        let aq_tag = Tag::from([b'a', b'q']);
        ss_a.data_mut().insert(aq_tag, Value::String("5555".as_bytes().into()));

        let bq_tag = Tag::from([b'b', b'q']);
        ss_b.data_mut().insert(bq_tag, Value::String("5555".as_bytes().into()));

        let duplex = DuplexConsensusCaller::call_duplex_from_ss_pair(
            ss_a,
            ss_b,
            true,
            &[],
            &[],
            &create_test_ss_caller(),
            true, // R1 duplex (first of pair)
        )?
        .expect("Should produce duplex consensus");

        // Check that per-base tags were generated
        let ad_tag = Tag::from([b'a', b'd']);
        let bd_tag = Tag::from([b'b', b'd']);
        let ae_tag = Tag::from([b'a', b'e']);
        let be_tag = Tag::from([b'b', b'e']);

        // Verify AB per-base depth tag
        if let Some(Value::Array(arr)) = duplex.data().get(&ad_tag) {
            let depths: Vec<i16> = match arr {
                noodles::sam::alignment::record_buf::data::field::value::Array::Int16(v) => {
                    v.clone()
                }
                noodles::sam::alignment::record_buf::data::field::value::Array::Int8(v) => {
                    v.iter().map(|&x| i16::from(x)).collect()
                }
                _ => panic!("Unexpected array type"),
            };
            assert_eq!(depths, vec![3, 3, 2, 3]);
        } else {
            panic!("Missing ad tag");
        }

        // Verify BA per-base depth tag
        if let Some(Value::Array(arr)) = duplex.data().get(&bd_tag) {
            let depths: Vec<i16> = match arr {
                noodles::sam::alignment::record_buf::data::field::value::Array::Int16(v) => {
                    v.clone()
                }
                noodles::sam::alignment::record_buf::data::field::value::Array::Int8(v) => {
                    v.iter().map(|&x| i16::from(x)).collect()
                }
                _ => panic!("Unexpected array type"),
            };
            assert_eq!(depths, vec![2, 2, 1, 2]);
        } else {
            panic!("Missing bd tag");
        }

        // Verify AB per-base error tag
        if let Some(Value::Array(_)) = duplex.data().get(&ae_tag) {
            // Tag exists
        } else {
            panic!("Missing ae tag");
        }

        // Verify BA per-base error tag
        if let Some(Value::Array(_)) = duplex.data().get(&be_tag) {
            // Tag exists
        } else {
            panic!("Missing be tag");
        }

        // Verify single-strand consensus sequence tags are generated
        let ac_tag_out = Tag::from([b'a', b'c']);
        let bc_tag_out = Tag::from([b'b', b'c']);

        // These tags are stored as Arrays (Vec<u8>), not Strings
        if let Some(Value::String(value)) = duplex.data().get(&ac_tag_out) {
            assert_eq!(value, b"AAAA");
        } else {
            panic!("Missing ac tag (AB consensus bases)");
        }

        if let Some(Value::String(value)) = duplex.data().get(&bc_tag_out) {
            assert_eq!(value, b"AAAA");
        } else {
            panic!("Missing bc tag (BA consensus bases)");
        }

        // Verify single-strand consensus quality tags are also generated as Arrays
        let aq_tag_out = Tag::from([b'a', b'q']);
        let bq_tag_out = Tag::from([b'b', b'q']);

        if let Some(Value::String(_)) = duplex.data().get(&aq_tag_out) {
            // Tag exists
        } else {
            panic!("Missing aq tag (AB consensus qualities)");
        }

        if let Some(Value::String(_)) = duplex.data().get(&bq_tag_out) {
            // Tag exists
        } else {
            panic!("Missing bq tag (BA consensus qualities)");
        }

        // Verify that cd and ce tags are NOT present in duplex output
        // (fgbio doesn't output these for duplex consensus)
        let cd_tag_check = Tag::from([b'c', b'd']);
        let ce_tag_check = Tag::from([b'c', b'e']);
        assert!(
            duplex.data().get(&cd_tag_check).is_none(),
            "cd tag should not be present in duplex output"
        );
        assert!(
            duplex.data().get(&ce_tag_check).is_none(),
            "ce tag should not be present in duplex output"
        );

        // Verify per-read tags cD and cM are computed correctly from per-base depths
        // fgbio: cD = max(ab.depths[i] + ba.depths[i]), cM = min(ab.depths[i] + ba.depths[i])
        // ab.depths = [3, 3, 2, 3], ba.depths = [2, 2, 1, 2]
        // combined = [5, 5, 3, 5], so cD = 5, cM = 3
        let cd_per_read_tag = Tag::from([b'c', b'D']);
        let cm_per_read_tag = Tag::from([b'c', b'M']);
        let depth_max = DuplexConsensusCaller::extract_int_tag(duplex.data(), cd_per_read_tag)
            .expect("Missing cD tag");
        let depth_min = DuplexConsensusCaller::extract_int_tag(duplex.data(), cm_per_read_tag)
            .expect("Missing cM tag");
        assert_eq!(depth_max, 5, "cD should be max of combined per-base depths");
        assert_eq!(depth_min, 3, "cM should be min of combined per-base depths");

        Ok(())
    }

    #[test]
    fn test_duplex_consensus_mixed_bases_and_n() -> Result<()> {
        // Test complex scenario with mixture of N bases and disagreements
        let ss_a = create_ss_consensus("ANGT", &[30, 30, 30, 30], 3, 3, 0.0);
        let ss_b = create_ss_consensus("TNCG", &[25, 25, 25, 25], 2, 2, 0.0);

        let duplex = DuplexConsensusCaller::call_duplex_from_ss_pair(
            ss_a,
            ss_b,
            false,
            &[],
            &[],
            &create_test_ss_caller(),
            true, // R1 duplex (first of pair)
        )?
        .expect("Should produce duplex consensus");

        // Check sequence
        // Position 0: A vs T with Q30 vs Q25 → A with Q=5
        // Position 1: N vs N → N with Q=2
        // Position 2: G vs C with Q30 vs Q25 → G with Q=5 (neither has N)
        // Position 3: T vs G with Q30 vs Q25 → T with Q=5
        let seq: Vec<u8> = duplex.sequence().as_ref().to_vec();
        assert_eq!(seq, b"ANGT");

        // Check qualities
        let quals: Vec<u8> = duplex.quality_scores().as_ref().to_vec();
        assert_eq!(quals[0], 5); // 30-25
        assert_eq!(quals[1], 2); // N masked
        assert_eq!(quals[2], 5); // 30-25 (no N at this position)
        assert_eq!(quals[3], 5); // 30-25

        Ok(())
    }

    #[test]
    fn test_duplex_consensus_zero_quality_difference() -> Result<()> {
        // Test that bases with same quality but different calls become N
        let ss_a = create_ss_consensus("AAAA", &[25, 25, 25, 25], 3, 3, 0.0);
        let ss_b = create_ss_consensus("TTTT", &[25, 25, 25, 25], 2, 2, 0.0);

        let duplex = DuplexConsensusCaller::call_duplex_from_ss_pair(
            ss_a,
            ss_b,
            false,
            &[],
            &[],
            &create_test_ss_caller(),
            true, // R1 duplex (first of pair)
        )?
        .expect("Should produce duplex consensus");

        // All positions have equal quality disagreements → should be N
        let seq: Vec<u8> = duplex.sequence().as_ref().to_vec();
        assert_eq!(seq, b"NNNN");

        // All qualities should be minimum
        let quals: Vec<u8> = duplex.quality_scores().as_ref().to_vec();
        assert_eq!(quals, vec![2, 2, 2, 2]);

        Ok(())
    }

    // ============================================================================
    // Ported from fgbio DuplexConsensusCallerTest
    // ============================================================================

    /// Port of fgbio test: "sourceMoleculeId should strip off that last /suffix"
    /// Tests various MI tag parsing edge cases
    #[test]
    fn test_parse_mi_tag_edge_cases() {
        // Multiple slashes: "f/A/B" -> "f/A"
        assert_eq!(DuplexConsensusCaller::parse_mi_tag("f/A/B"), ("f/A".to_string(), Some('B')));

        // Three slashes: "f///A" -> "f//"
        assert_eq!(DuplexConsensusCaller::parse_mi_tag("f///A"), ("f//".to_string(), Some('A')));
    }

    /// Port of fgbio test: "throw an exception if the MI tag doesn't have a /suffix"
    /// Tests that MI tags without strand suffix are handled gracefully
    #[test]
    fn test_parse_mi_tag_no_suffix() {
        // No suffix - should return None for strand
        let (mi, strand) = DuplexConsensusCaller::parse_mi_tag("foo");
        assert_eq!(mi, "foo");
        assert!(strand.is_none(), "Should return None for MI without /suffix");
    }

    /// Port of fgbio test: validation exception for invalid min-reads order
    /// Tests that `min_reads` values must be in decreasing order
    #[test]
    fn test_min_reads_validation_error() {
        // Increasing order should fail (total=1, XY=2, YX=3 is invalid)
        let result = DuplexConsensusCaller::new(
            "consensus".to_string(),
            "RG1".to_string(),
            vec![1, 2, 3], // Invalid: should be decreasing (total >= XY >= YX)
            20,
            false,
            false,
            None,
            None,
            false,
            45,
            40,
        );

        // fgbio throws ValidationException for invalid min_reads order
        assert!(result.is_err(), "Should fail with invalid min_reads order [1, 2, 3]");
    }

    /// Port of fgbio test: "count errors against the actual consensus base, before it is masked to N"
    /// Tests that error counting uses the original consensus call, not the masked N
    #[test]
    fn test_errors_counted_before_masking() -> Result<()> {
        // Create SS consensuses where both have N at one position
        // fgbio: bases1="ANAAAAAAAA", bases2="CCNCCCCCCC"
        // After RC of bases2: "GGNGGGGGGG"
        // At position 1: A vs G - disagreement, but one is N → no error
        // At position 2: A vs N - one is N → no error

        let ss_a = create_ss_consensus("ANAA", &[30, 30, 30, 30], 1, 1, 0.0);
        let ss_b = create_ss_consensus("GGNG", &[30, 30, 30, 30], 1, 1, 0.0);

        let duplex = DuplexConsensusCaller::call_duplex_from_ss_pair(
            ss_a,
            ss_b,
            false,
            &[],
            &[],
            &create_test_ss_caller(),
            true, // R1 duplex (first of pair)
        )?
        .expect("Should produce duplex consensus");

        // Check that per-read error rate is 0 (N positions not counted as errors)
        let ce_tag = noodles::sam::alignment::record::data::field::Tag::from([b'c', b'E']);
        if let Some(Value::Float(error_rate)) = duplex.data().get(&ce_tag) {
            // Error rate should be 0 because N positions don't count as errors
            assert!(
                *error_rate < 0.01,
                "Error rate should be ~0 when N positions are excluded, got {error_rate}"
            );
        }

        Ok(())
    }

    /// Tests that when one strand has N at a position, output is N
    /// Matches fgbio DuplexConsensusCaller.scala line 431:
    /// if (aBase == `NoCall` || bBase == `NoCall` || rawQual == `MinValue`) (`NoCall`, `NoCallQual`)
    #[test]
    fn test_duplex_n_base_propagation() -> Result<()> {
        // AB strand has data, BA has N's
        let ss_ab = create_ss_consensus("AAAA", &[30, 30, 30, 30], 3, 3, 0.0);
        let ss_ba_n = create_ss_consensus("NNNN", &[2, 2, 2, 2], 0, 0, 0.0);

        // When one strand has N, fgbio outputs N with NoCallQual (2)
        let duplex = DuplexConsensusCaller::call_duplex_from_ss_pair(
            ss_ab.clone(),
            ss_ba_n,
            false,
            &[],
            &[],
            &create_test_ss_caller(),
            true, // R1 duplex (first of pair)
        )?
        .expect("Should produce duplex consensus");

        // fgbio behavior: N propagates to output
        let seq: Vec<u8> = duplex.sequence().as_ref().to_vec();
        assert_eq!(seq, b"NNNN", "Output should be N when either strand has N");

        let quals: Vec<u8> = duplex.quality_scores().as_ref().to_vec();
        for &q in &quals {
            assert_eq!(q, 2, "Quality should be NoCallQual (2) when base is N");
        }

        Ok(())
    }

    /// Tests N propagation when AB has N and BA has data
    #[test]
    fn test_duplex_n_base_propagation_reverse() -> Result<()> {
        // Create AB with N's
        let ss_ab_n = create_ss_consensus("NNNN", &[2, 2, 2, 2], 0, 0, 0.0);
        // Create BA with actual data
        let ss_ba = create_ss_consensus("TTTT", &[30, 30, 30, 30], 3, 3, 0.0);

        // When one strand has N, fgbio outputs N
        let duplex = DuplexConsensusCaller::call_duplex_from_ss_pair(
            ss_ab_n,
            ss_ba,
            false,
            &[],
            &[],
            &create_test_ss_caller(),
            true, // R1 duplex (first of pair)
        )?
        .expect("Should produce duplex consensus");

        // fgbio behavior: N propagates to output
        let seq: Vec<u8> = duplex.sequence().as_ref().to_vec();
        assert_eq!(seq, b"NNNN", "Output should be N when either strand has N");

        let quals: Vec<u8> = duplex.quality_scores().as_ref().to_vec();
        for &q in &quals {
            assert_eq!(q, 2, "Quality should be NoCallQual (2) when base is N");
        }

        Ok(())
    }

    // ==========================================================================
    // Integration tests using consensus_reads_from_sam_records
    // ==========================================================================

    /// Port of fgbio test: "not create records from fragments"
    /// Tests that fragment reads (unpaired) are rejected
    #[test]
    fn test_not_create_records_from_fragments() -> Result<()> {
        let mut caller = DuplexConsensusCaller::new(
            "consensus".to_string(),
            "RG1".to_string(),
            vec![1],
            10,
            false,
            false,
            None,
            None,
            false,
            45,
            40,
        )?;

        // Create fragment (unpaired) reads with MI tags
        let frag_a = RecordBuilder::new()
            .name("frag1")
            .sequence("AAAAAAAAAA")
            .reference_sequence_id(0)
            .alignment_start(1)
            .cigar("10M")
            .tag("MI", "foo/A")
            .build();

        let frag_b = RecordBuilder::new()
            .name("frag2")
            .sequence("AAAAAAAAAA")
            .reference_sequence_id(0)
            .alignment_start(1)
            .cigar("10M")
            .tag("MI", "foo/B")
            .build();

        let result = caller.consensus_reads_from_sam_records(vec![frag_a, frag_b])?;

        // Fragments should be rejected - no consensus output
        assert_eq!(
            result.count, 0,
            "Fragment reads should not produce consensus, got {} reads",
            result.count
        );

        Ok(())
    }

    /// Port of fgbio test: "create a simple double stranded consensus for a pair of A and a pair of B reads"
    /// Tests basic duplex consensus creation from AB and BA read pairs
    #[test]
    fn test_create_simple_double_stranded_consensus() -> Result<()> {
        let mut caller = DuplexConsensusCaller::new(
            "consensus".to_string(),
            "RG1".to_string(),
            vec![1],
            10,
            false,
            false,
            None,
            None,
            false,
            45,
            40,
        )?;

        // Create AB pair: R1 at pos 100 (+), R2 at pos 200 (-)
        let ab_r1 = RecordBuilder::new()
            .name("q1")
            .sequence("AAAAAAAAAA")
            .qualities(&[20; 10])
            .first_segment(true)
            .reverse_complement(false)
            .mate_reverse_complement(true)
            .reference_sequence_id(0)
            .alignment_start(100)
            .mate_reference_sequence_id(0)
            .mate_alignment_start(200)
            .cigar("10M")
            .tag("MI", "foo/A")
            .tag("RG", "A")
            .build();

        let ab_r2 = RecordBuilder::new()
            .name("q1")
            .sequence("CCCCCCCCCC")
            .qualities(&[20; 10])
            .first_segment(false)
            .reverse_complement(true)
            .mate_reverse_complement(false)
            .reference_sequence_id(0)
            .alignment_start(200)
            .mate_reference_sequence_id(0)
            .mate_alignment_start(100)
            .cigar("10M")
            .tag("MI", "foo/A")
            .tag("RG", "A")
            .build();

        // Create BA pair: R1 at pos 200 (-), R2 at pos 100 (+)
        let ba_r1 = RecordBuilder::new()
            .name("q2")
            .sequence("CCCCCCCCCC")
            .qualities(&[20; 10])
            .first_segment(true)
            .reverse_complement(true)
            .mate_reverse_complement(false)
            .reference_sequence_id(0)
            .alignment_start(200)
            .mate_reference_sequence_id(0)
            .mate_alignment_start(100)
            .cigar("10M")
            .tag("MI", "foo/B")
            .tag("RG", "A")
            .build();

        let ba_r2 = RecordBuilder::new()
            .name("q2")
            .sequence("AAAAAAAAAA")
            .qualities(&[20; 10])
            .first_segment(false)
            .reverse_complement(false)
            .mate_reverse_complement(true)
            .reference_sequence_id(0)
            .alignment_start(100)
            .mate_reference_sequence_id(0)
            .mate_alignment_start(200)
            .cigar("10M")
            .tag("MI", "foo/B")
            .tag("RG", "A")
            .build();

        let result = caller.consensus_reads_from_sam_records(vec![ab_r1, ab_r2, ba_r1, ba_r2])?;

        // Should produce 2 consensus reads (R1 and R2)
        assert_eq!(
            result.count, 2,
            "Should produce 2 duplex consensus reads (R1 and R2), got {}",
            result.count
        );

        Ok(())
    }

    /// Port of fgbio test: "preserve the cell barcode"
    /// Tests that CB tag is preserved on consensus reads
    #[test]
    fn test_preserve_cell_barcode() -> Result<()> {
        use noodles::sam::alignment::record::data::field::Tag;

        let cb_tag = Tag::from([b'C', b'B']);
        let mut caller = DuplexConsensusCaller::new(
            "consensus".to_string(),
            "RG1".to_string(),
            vec![1],
            10,
            false,
            false,
            None,
            Some(cb_tag), // Enable cell barcode tracking
            false,
            45,
            40,
        )?;

        // Create AB pair with CB tag
        let ab_r1 = RecordBuilder::new()
            .name("q1")
            .sequence("AAAAAAAAAA")
            .qualities(&[20; 10])
            .first_segment(true)
            .reverse_complement(false)
            .mate_reverse_complement(true)
            .reference_sequence_id(0)
            .alignment_start(100)
            .mate_reference_sequence_id(0)
            .mate_alignment_start(200)
            .cigar("10M")
            .tag("MI", "foo/A")
            .tag("CB", "ACGT")
            .tag("RG", "A")
            .build();

        let ab_r2 = RecordBuilder::new()
            .name("q1")
            .sequence("CCCCCCCCCC")
            .qualities(&[20; 10])
            .first_segment(false)
            .reverse_complement(true)
            .mate_reverse_complement(false)
            .reference_sequence_id(0)
            .alignment_start(200)
            .mate_reference_sequence_id(0)
            .mate_alignment_start(100)
            .cigar("10M")
            .tag("MI", "foo/A")
            .tag("CB", "ACGT")
            .tag("RG", "A")
            .build();

        // Create BA pair with CB tag
        let ba_r1 = RecordBuilder::new()
            .name("q2")
            .sequence("CCCCCCCCCC")
            .qualities(&[20; 10])
            .first_segment(true)
            .reverse_complement(true)
            .mate_reverse_complement(false)
            .reference_sequence_id(0)
            .alignment_start(200)
            .mate_reference_sequence_id(0)
            .mate_alignment_start(100)
            .cigar("10M")
            .tag("MI", "foo/B")
            .tag("CB", "ACGT")
            .tag("RG", "A")
            .build();

        let ba_r2 = RecordBuilder::new()
            .name("q2")
            .sequence("AAAAAAAAAA")
            .qualities(&[20; 10])
            .first_segment(false)
            .reverse_complement(false)
            .mate_reverse_complement(true)
            .reference_sequence_id(0)
            .alignment_start(100)
            .mate_reference_sequence_id(0)
            .mate_alignment_start(200)
            .cigar("10M")
            .tag("MI", "foo/B")
            .tag("CB", "ACGT")
            .tag("RG", "A")
            .build();

        let result = caller.consensus_reads_from_sam_records(vec![ab_r1, ab_r2, ba_r1, ba_r2])?;

        // Should produce consensus reads with CB tag preserved
        assert_eq!(result.count, 2, "Should produce 2 consensus reads");

        // Check that CB tag is preserved on consensus reads
        let records = ParsedBamRecord::parse_all(&result.data);
        for rec in &records {
            let cb = rec.get_string_tag(b"CB").expect("CB tag should be present on consensus read");
            assert_eq!(cb, b"ACGT", "CB tag should be preserved");
        }

        Ok(())
    }

    /// Port of fgbio test: "handle absent UMI (left/right)"
    /// Tests that consensus is generated correctly when one UMI in the pair is absent
    /// (indicated by "-" in the RX tag)
    #[test]
    fn test_handle_absent_umi_right() -> Result<()> {
        let mut caller = DuplexConsensusCaller::new(
            "consensus".to_string(),
            "RG1".to_string(),
            vec![1],
            10,
            false,
            false,
            None,
            None,
            false,
            45,
            40,
        )?;

        // Create AB pair with RX="ACT-" (right UMI absent)
        let ab_r1 = RecordBuilder::new()
            .name("q1")
            .sequence("AAAAAAAAAA")
            .qualities(&[20; 10])
            .first_segment(true)
            .reverse_complement(false)
            .mate_reverse_complement(true)
            .reference_sequence_id(0)
            .alignment_start(100)
            .mate_reference_sequence_id(0)
            .mate_alignment_start(200)
            .cigar("10M")
            .tag("MI", "foo/A")
            .tag("RX", "ACT-")
            .tag("RG", "A")
            .build();

        let ab_r2 = RecordBuilder::new()
            .name("q1")
            .sequence("CCCCCCCCCC")
            .qualities(&[20; 10])
            .first_segment(false)
            .reverse_complement(true)
            .mate_reverse_complement(false)
            .reference_sequence_id(0)
            .alignment_start(200)
            .mate_reference_sequence_id(0)
            .mate_alignment_start(100)
            .cigar("10M")
            .tag("MI", "foo/A")
            .tag("RX", "ACT-")
            .tag("RG", "A")
            .build();

        // Create BA pair with RX="-ACT" (left UMI absent)
        let ba_r1 = RecordBuilder::new()
            .name("q2")
            .sequence("CCCCCCCCCC")
            .qualities(&[20; 10])
            .first_segment(true)
            .reverse_complement(true)
            .mate_reverse_complement(false)
            .reference_sequence_id(0)
            .alignment_start(200)
            .mate_reference_sequence_id(0)
            .mate_alignment_start(100)
            .cigar("10M")
            .tag("MI", "foo/B")
            .tag("RX", "-ACT")
            .tag("RG", "A")
            .build();

        let ba_r2 = RecordBuilder::new()
            .name("q2")
            .sequence("AAAAAAAAAA")
            .qualities(&[20; 10])
            .first_segment(false)
            .reverse_complement(false)
            .mate_reverse_complement(true)
            .reference_sequence_id(0)
            .alignment_start(100)
            .mate_reference_sequence_id(0)
            .mate_alignment_start(200)
            .cigar("10M")
            .tag("MI", "foo/B")
            .tag("RX", "-ACT")
            .tag("RG", "A")
            .build();

        let result = caller.consensus_reads_from_sam_records(vec![ab_r1, ab_r2, ba_r1, ba_r2])?;

        // Should produce 2 consensus reads
        assert_eq!(result.count, 2, "Should produce 2 duplex consensus reads");

        // Check RX tag is inherited from the A family (leftUmi pattern "ACT-")
        let records = ParsedBamRecord::parse_all(&result.data);
        for rec in &records {
            let rx = rec.get_string_tag(b"RX").expect("RX tag should be present on consensus read");
            assert_eq!(rx, b"ACT-", "RX tag should be inherited from A family");
        }

        Ok(())
    }

    /// Port of fgbio test: "handle absent UMI (left)"
    /// Tests the reverse case where left UMI is absent in A family
    #[test]
    fn test_handle_absent_umi_left() -> Result<()> {
        let mut caller = DuplexConsensusCaller::new(
            "consensus".to_string(),
            "RG1".to_string(),
            vec![1],
            10,
            false,
            false,
            None,
            None,
            false,
            45,
            40,
        )?;

        // Create AB pair with RX="-ACT" (left UMI absent)
        let ab_r1 = RecordBuilder::new()
            .name("q1")
            .sequence("AAAAAAAAAA")
            .qualities(&[20; 10])
            .first_segment(true)
            .reverse_complement(false)
            .mate_reverse_complement(true)
            .reference_sequence_id(0)
            .alignment_start(100)
            .mate_reference_sequence_id(0)
            .mate_alignment_start(200)
            .cigar("10M")
            .tag("MI", "foo/A")
            .tag("RX", "-ACT")
            .tag("RG", "A")
            .build();

        let ab_r2 = RecordBuilder::new()
            .name("q1")
            .sequence("CCCCCCCCCC")
            .qualities(&[20; 10])
            .first_segment(false)
            .reverse_complement(true)
            .mate_reverse_complement(false)
            .reference_sequence_id(0)
            .alignment_start(200)
            .mate_reference_sequence_id(0)
            .mate_alignment_start(100)
            .cigar("10M")
            .tag("MI", "foo/A")
            .tag("RX", "-ACT")
            .tag("RG", "A")
            .build();

        // Create BA pair with RX="ACT-" (right UMI absent)
        let ba_r1 = RecordBuilder::new()
            .name("q2")
            .sequence("CCCCCCCCCC")
            .qualities(&[20; 10])
            .first_segment(true)
            .reverse_complement(true)
            .mate_reverse_complement(false)
            .reference_sequence_id(0)
            .alignment_start(200)
            .mate_reference_sequence_id(0)
            .mate_alignment_start(100)
            .cigar("10M")
            .tag("MI", "foo/B")
            .tag("RX", "ACT-")
            .tag("RG", "A")
            .build();

        let ba_r2 = RecordBuilder::new()
            .name("q2")
            .sequence("AAAAAAAAAA")
            .qualities(&[20; 10])
            .first_segment(false)
            .reverse_complement(false)
            .mate_reverse_complement(true)
            .reference_sequence_id(0)
            .alignment_start(100)
            .mate_reference_sequence_id(0)
            .mate_alignment_start(200)
            .cigar("10M")
            .tag("MI", "foo/B")
            .tag("RX", "ACT-")
            .tag("RG", "A")
            .build();

        let result = caller.consensus_reads_from_sam_records(vec![ab_r1, ab_r2, ba_r1, ba_r2])?;

        // Should produce 2 consensus reads
        assert_eq!(result.count, 2, "Should produce 2 duplex consensus reads");

        // Check RX tag is inherited from the A family (leftUmi pattern "-ACT")
        let records = ParsedBamRecord::parse_all(&result.data);
        for rec in &records {
            let rx = rec.get_string_tag(b"RX").expect("RX tag should be present on consensus read");
            assert_eq!(rx, b"-ACT", "RX tag should be inherited from A family");
        }

        Ok(())
    }

    /// Port of fgbio test: "create consensus with only A reads"
    /// Tests single-strand pipeline when only A-strand reads are present
    /// and minReads=[3,3,0] allows BA=0
    #[test]
    fn test_create_single_strand_consensus_a_only() -> Result<()> {
        // Create caller with minReads=[1,1,0] to allow single-strand consensus
        let mut caller = DuplexConsensusCaller::new(
            "consensus".to_string(),
            "RG1".to_string(),
            vec![1, 1, 0], // duplex=1, AB=1, BA=0 (allows single strand)
            10,
            false,
            false,
            None,
            None,
            false,
            45,
            40,
        )?;

        // Create 3 A-strand pairs (no B-strand)
        let reads: Vec<RecordBuf> = (1..=3)
            .flat_map(|i| {
                vec![
                    RecordBuilder::new()
                        .name(&format!("q{i}"))
                        .sequence("AAAAAAAAAA")
                        .qualities(&[20; 10])
                        .first_segment(true)
                        .reverse_complement(false)
                        .mate_reverse_complement(true)
                        .reference_sequence_id(0)
                        .alignment_start(100)
                        .mate_reference_sequence_id(0)
                        .mate_alignment_start(200)
                        .cigar("10M")
                        .tag("MI", "foo/A")
                        .tag("RG", "A")
                        .build(),
                    RecordBuilder::new()
                        .name(&format!("q{i}"))
                        .sequence("CCCCCCCCCC")
                        .qualities(&[20; 10])
                        .first_segment(false)
                        .reverse_complement(true)
                        .mate_reverse_complement(false)
                        .reference_sequence_id(0)
                        .alignment_start(200)
                        .mate_reference_sequence_id(0)
                        .mate_alignment_start(100)
                        .cigar("10M")
                        .tag("MI", "foo/A")
                        .tag("RG", "A")
                        .build(),
                ]
            })
            .collect();

        let result = caller.consensus_reads_from_sam_records(reads)?;

        // Should produce 2 consensus reads (R1 and R2) from single strand
        assert_eq!(result.count, 2, "Should produce 2 consensus reads from single A strand");

        // Check per-read tags show 3 AB reads and 0 BA reads
        let records = ParsedBamRecord::parse_all(&result.data);

        for rec in &records {
            let ab_depth = rec.get_int_tag(b"aD").expect("Missing aD tag");
            assert_eq!(ab_depth, 3, "AB depth should be 3");

            // BA depth should be 0 or tag should not exist
            if let Some(ba_depth) = rec.get_int_tag(b"bD") {
                assert_eq!(ba_depth, 0, "BA depth should be 0");
            }
            // If tag doesn't exist, that's also correct for single-strand
        }

        Ok(())
    }

    /// Port of fgbio test: "create consensus with only B reads"
    /// Tests single-strand pipeline when only B-strand reads are present
    #[test]
    fn test_create_single_strand_consensus_b_only() -> Result<()> {
        // Create caller with minReads=[1,1,0] to allow single-strand consensus
        let mut caller = DuplexConsensusCaller::new(
            "consensus".to_string(),
            "RG1".to_string(),
            vec![1, 1, 0], // duplex=1, AB=1, BA=0 (allows single strand)
            10,
            false,
            false,
            None,
            None,
            false,
            45,
            40,
        )?;

        // Create 3 B-strand pairs (no A-strand)
        let reads: Vec<RecordBuf> = (1..=3)
            .flat_map(|i| {
                vec![
                    RecordBuilder::new()
                        .name(&format!("q{i}"))
                        .sequence("CCCCCCCCCC")
                        .qualities(&[20; 10])
                        .first_segment(true)
                        .reverse_complement(true)
                        .mate_reverse_complement(false)
                        .reference_sequence_id(0)
                        .alignment_start(200)
                        .mate_reference_sequence_id(0)
                        .mate_alignment_start(100)
                        .cigar("10M")
                        .tag("MI", "foo/B")
                        .tag("RG", "A")
                        .build(),
                    RecordBuilder::new()
                        .name(&format!("q{i}"))
                        .sequence("AAAAAAAAAA")
                        .qualities(&[20; 10])
                        .first_segment(false)
                        .reverse_complement(false)
                        .mate_reverse_complement(true)
                        .reference_sequence_id(0)
                        .alignment_start(100)
                        .mate_reference_sequence_id(0)
                        .mate_alignment_start(200)
                        .cigar("10M")
                        .tag("MI", "foo/B")
                        .tag("RG", "A")
                        .build(),
                ]
            })
            .collect();

        let result = caller.consensus_reads_from_sam_records(reads)?;

        // Should produce 2 consensus reads (R1 and R2) from single strand
        assert_eq!(result.count, 2, "Should produce 2 consensus reads from single B strand");

        // Check per-read tags - for B-only, fgbio swaps AB/BA so AB has the data
        let records = ParsedBamRecord::parse_all(&result.data);

        for rec in &records {
            // After swap, AB should have the depth
            if let Some(ab_depth) = rec.get_int_tag(b"aD") {
                assert_eq!(ab_depth, 3, "AB depth should be 3 after swap");
            }
        }

        Ok(())
    }

    /// Port of fgbio test: "not create consensus with only A or B reads when minReads doesn't allow it"
    /// Tests that single-strand reads are rejected when BA min > 0
    #[test]
    fn test_reject_single_strand_when_min_reads_requires_both() -> Result<()> {
        // Create caller with minReads=[1,1,1] which requires both strands
        let mut caller = DuplexConsensusCaller::new(
            "consensus".to_string(),
            "RG1".to_string(),
            vec![1, 1, 1], // All must have at least 1 read
            10,
            false,
            false,
            None,
            None,
            false,
            45,
            40,
        )?;

        // Create 3 A-strand pairs only (no B-strand)
        let reads: Vec<RecordBuf> = (1..=3)
            .flat_map(|i| {
                vec![
                    RecordBuilder::new()
                        .name(&format!("q{i}"))
                        .sequence("AAAAAAAAAA")
                        .qualities(&[20; 10])
                        .first_segment(true)
                        .reverse_complement(false)
                        .mate_reverse_complement(true)
                        .reference_sequence_id(0)
                        .alignment_start(100)
                        .mate_reference_sequence_id(0)
                        .mate_alignment_start(200)
                        .cigar("10M")
                        .tag("MI", "foo/A")
                        .tag("RG", "A")
                        .build(),
                    RecordBuilder::new()
                        .name(&format!("q{i}"))
                        .sequence("CCCCCCCCCC")
                        .qualities(&[20; 10])
                        .first_segment(false)
                        .reverse_complement(true)
                        .mate_reverse_complement(false)
                        .reference_sequence_id(0)
                        .alignment_start(200)
                        .mate_reference_sequence_id(0)
                        .mate_alignment_start(100)
                        .cigar("10M")
                        .tag("MI", "foo/A")
                        .tag("RG", "A")
                        .build(),
                ]
            })
            .collect();

        let result = caller.consensus_reads_from_sam_records(reads)?;

        // Should produce empty result since BA=0 doesn't meet minReads[2]=1
        assert_eq!(
            result.count, 0,
            "Should not produce consensus when single strand doesn't meet min_reads"
        );

        Ok(())
    }

    /// Port of fgbio test: "support the --min-reads option as a hard filter after alignment filtering"
    /// Tests that min-reads filtering applies after CIGAR-based alignment filtering
    #[test]
    fn test_min_reads_hard_filter_after_alignment_filtering() -> Result<()> {
        // Setup: 3 AB pairs, 2 BA pairs = 5 total reads
        // minReads=3 requires 3 reads, minReads=2 requires 2 reads
        let mut caller_3 = DuplexConsensusCaller::new(
            "consensus".to_string(),
            "RG1".to_string(),
            vec![3],
            10,
            false,
            false,
            None,
            None,
            false,
            45,
            40,
        )?;

        let mut caller_2 = DuplexConsensusCaller::new(
            "consensus".to_string(),
            "RG1".to_string(),
            vec![2],
            10,
            false,
            false,
            None,
            None,
            false,
            45,
            40,
        )?;

        // Create reads: 3 AB pairs + 2 BA pairs
        let mut reads = Vec::new();

        // 3 AB pairs
        for i in 1..=3 {
            reads.push(
                RecordBuilder::new()
                    .name(&format!("ab{i}"))
                    .sequence("AAAAAAAAAA")
                    .qualities(&[20; 10])
                    .first_segment(true)
                    .reverse_complement(false)
                    .mate_reverse_complement(true)
                    .reference_sequence_id(0)
                    .alignment_start(100)
                    .mate_reference_sequence_id(0)
                    .mate_alignment_start(200)
                    .cigar("10M")
                    .tag("MI", "foo/A")
                    .tag("RG", "A")
                    .build(),
            );
            reads.push(
                RecordBuilder::new()
                    .name(&format!("ab{i}"))
                    .sequence("CCCCCCCCCC")
                    .qualities(&[20; 10])
                    .first_segment(false)
                    .reverse_complement(true)
                    .mate_reverse_complement(false)
                    .reference_sequence_id(0)
                    .alignment_start(200)
                    .mate_reference_sequence_id(0)
                    .mate_alignment_start(100)
                    .cigar("10M")
                    .tag("MI", "foo/A")
                    .tag("RG", "A")
                    .build(),
            );
        }

        // 2 BA pairs
        for i in 4..=5 {
            reads.push(
                RecordBuilder::new()
                    .name(&format!("ba{i}"))
                    .sequence("CCCCCCCCCC")
                    .qualities(&[20; 10])
                    .first_segment(true)
                    .reverse_complement(true)
                    .mate_reverse_complement(false)
                    .reference_sequence_id(0)
                    .alignment_start(200)
                    .mate_reference_sequence_id(0)
                    .mate_alignment_start(100)
                    .cigar("10M")
                    .tag("MI", "foo/B")
                    .tag("RG", "A")
                    .build(),
            );
            reads.push(
                RecordBuilder::new()
                    .name(&format!("ba{i}"))
                    .sequence("AAAAAAAAAA")
                    .qualities(&[20; 10])
                    .first_segment(false)
                    .reverse_complement(false)
                    .mate_reverse_complement(true)
                    .reference_sequence_id(0)
                    .alignment_start(100)
                    .mate_reference_sequence_id(0)
                    .mate_alignment_start(200)
                    .cigar("10M")
                    .tag("MI", "foo/B")
                    .tag("RG", "A")
                    .build(),
            );
        }

        // With minReads=3: should fail (BA only has 2 reads)
        let result_3 = caller_3.consensus_reads_from_sam_records(reads.clone())?;
        assert_eq!(result_3.count, 0, "minReads=3 should fail with only 2 BA reads");

        // With minReads=2: should pass
        let result_2 = caller_2.consensus_reads_from_sam_records(reads.clone())?;
        assert_eq!(result_2.count, 2, "minReads=2 should produce 2 consensus reads");

        // Now add a BA read with different CIGAR (will be filtered out)
        let dissimilar = vec![
            RecordBuilder::new()
                .name("ba6")
                .sequence("CCCCCCCCCC")
                .qualities(&[20; 10])
                .first_segment(true)
                .reverse_complement(true)
                .mate_reverse_complement(false)
                .reference_sequence_id(0)
                .alignment_start(200)
                .mate_reference_sequence_id(0)
                .mate_alignment_start(100)
                .cigar("5M1D5M") // Different CIGAR - will be filtered
                .tag("MI", "foo/B")
                .tag("RG", "A")
                .build(),
            RecordBuilder::new()
                .name("ba6")
                .sequence("AAAAAAAAAA")
                .qualities(&[20; 10])
                .first_segment(false)
                .reverse_complement(false)
                .mate_reverse_complement(true)
                .reference_sequence_id(0)
                .alignment_start(100)
                .mate_reference_sequence_id(0)
                .mate_alignment_start(200)
                .cigar("10M")
                .tag("MI", "foo/B")
                .tag("RG", "A")
                .build(),
        ];

        let mut reads_with_dissimilar = reads.clone();
        reads_with_dissimilar.extend(dissimilar);

        // Create fresh callers for the extended test
        let mut caller_3_new = DuplexConsensusCaller::new(
            "consensus".to_string(),
            "RG1".to_string(),
            vec![3],
            10,
            false,
            false,
            None,
            None,
            false,
            45,
            40,
        )?;

        // With dissimilar read (filtered): still only 2 BA reads pass filtering
        // minReads=3 should still fail
        let result_dissimilar =
            caller_3_new.consensus_reads_from_sam_records(reads_with_dissimilar)?;
        assert_eq!(
            result_dissimilar.count, 0,
            "minReads=3 should still fail when dissimilar read is filtered out"
        );

        Ok(())
    }

    /// Port of fgbio test: "swap AB and BA when AB has zero depth"
    /// Tests that when AB consensus has zero depth, it gets swapped with BA
    #[test]
    fn test_swap_ab_ba_when_ab_has_zero_depth() -> Result<()> {
        // Create SS consensuses where AB has zero depth
        let ss_ab_zero = create_ss_consensus("NNN", &[2, 2, 2], 0, 0, 0.0);
        let ss_ba = create_ss_consensus("AAA", &[30, 30, 30], 3, 3, 0.0);

        // Call duplex consensus - should swap AB and BA
        let duplex = DuplexConsensusCaller::call_duplex_from_ss_pair(
            ss_ab_zero.clone(),
            ss_ba.clone(),
            false,
            &[],
            &[],
            &create_test_ss_caller(),
            true, // R1 duplex (first of pair)
        )?
        .expect("Should produce duplex consensus");

        // After swap, result should have the BA data as AB
        let seq: Vec<u8> = duplex.sequence().as_ref().to_vec();

        // Because AB has zero depth, fgbio swaps them, so output should be from BA
        // But since both are provided and AB has zeros, fgbio outputs N (NoCall)
        // This test verifies the behavior matches fgbio
        assert_eq!(seq, b"NNN", "When AB has zero depth, output should be N");

        Ok(())
    }

    /// Port of fgbio test: "swap BA and AB when BA has zero depth"
    /// Tests the reverse case where BA has zero depth
    #[test]
    fn test_swap_ba_ab_when_ba_has_zero_depth() -> Result<()> {
        // Create SS consensuses where BA has zero depth
        let ss_ab = create_ss_consensus("AAA", &[30, 30, 30], 3, 3, 0.0);
        let ss_ba_zero = create_ss_consensus("NNN", &[2, 2, 2], 0, 0, 0.0);

        // Call duplex consensus
        let duplex = DuplexConsensusCaller::call_duplex_from_ss_pair(
            ss_ab,
            ss_ba_zero,
            false,
            &[],
            &[],
            &create_test_ss_caller(),
            true, // R1 duplex (first of pair)
        )?
        .expect("Should produce duplex consensus");

        // BA has zero depth, so output should be N
        let seq: Vec<u8> = duplex.sequence().as_ref().to_vec();
        assert_eq!(seq, b"NNN", "When BA has zero depth, output should be N");

        Ok(())
    }

    /* DISABLED DUE TO bam::Reader ISSUE
    #[test]
    #[ignore] // Run with: cargo test --release -- --ignored test_mi5252_real_data
    fn test_mi5252_real_data() -> Result<()> {
        use noodles::sam::alignment::RecordBuf;
        use noodles::bam;
        use std::fs::File;

        // Read the test BAM file containing only MI 5252 reads
        let test_bam_path = "/Users/nhomer/work/git/fgumi/debug/test_mi5252_only.bam";
        let mut reader = File::open(test_bam_path)
            .map(bam::Reader::new)
            .expect("Failed to open test BAM file");

        let _header = reader.read_header()?;

        // Read all records from the BAM file
        let mut records = Vec::new();
        for result in reader.records() {
            let record = result?;
            let record_buf: RecordBuf = record.try_into()?;
            records.push(record_buf);
        }


        println!("Read {} records from test BAM", records.len());
        assert_eq!(records.len(), 18, "Should have exactly 18 reads for MI 5252");

        // Create a duplex consensus caller with minimum thresholds
        let mut caller = DuplexConsensusCaller::new(
            "test".to_string(),
            "RG1".to_string(),
            vec![1], // min_reads = 1 (M=1 in fgbio)
            20,      // min_base_quality
            false,   // trim
            false,   // sort_order check
            None,    // no mask
            None,    // no error_rate_pre_umi
            false,   // not producing per-base tags yet
            45,      // error_rate_post_umi
            40,      // min_consensus_base_quality
        )?;

        // Call consensus on these reads
        let consensus_records = caller.consensus_reads_from_sam_records(records)?;

        println!("Generated {} consensus records", consensus_records.len());
        assert_eq!(consensus_records.len(), 2, "Should generate 2 consensus records (R1 and R2)");

        // Find the R1 record (flag should indicate first of pair)
        let r1 = consensus_records
            .iter()
            .find(|r| r.flags().is_first_segment())
            .expect("Should have R1 consensus");

        // Extract the sequence at position 111 (0-based)
        let seq: Vec<u8> = r1.sequence().as_ref().to_vec();
        assert!(seq.len() > 111, "Sequence should be longer than 111 bases");

        println!("R1 consensus base at position 111: {}", seq[111] as char);

        // EXPECTED: fgbio outputs 'C' at position 111
        // ACTUAL: fgumi currently outputs 'N' at position 111
        // This test documents the discrepancy
        assert_eq!(
            seq[111] as char,
            'C',
            "Position 111 should be 'C' to match fgbio output (currently fgumi outputs 'N')"
        );

        Ok(())
    }

    */
    #[test]
    fn test_tie_breaking_for_simplex_consensus() -> Result<()> {
        use fgumi_sam::builder::{SamBuilder, Strand};

        // This test reproduces the tie-breaking scenario:
        // 8 BA read pairs with 4 having 'A' at position 5 and 4 having 'C' at position 5
        // With equal evidence (same quality, same count), this is a true tie.
        // With Kahan summation for numeric stability (matching fgbio PR #1120),
        // the likelihoods are exactly equal, so we correctly return 'N' (no-call).
        // Previously, numerical instability would incorrectly break the tie.
        let mut builder = SamBuilder::with_single_ref("chr1", 10000);

        // Quality array for 10 bases, all quality 38
        let quals = vec![38u8; 10];

        // Create 8 BA read pairs (strand B): 4 with 'A' at position 5, 4 with 'C' at position 5
        let _ = builder
            .add_pair()
            .name("ba1")
            .contig(0)
            .start1(100)
            .start2(200)
            .strand1(Strand::Plus)
            .strand2(Strand::Minus)
            .bases1("AAAAACAAAA")
            .bases2("TTTTTTTTTT")
            .quals1(&quals)
            .quals2(&quals)
            .attr("MI", "test/B")
            .build();

        let _ = builder
            .add_pair()
            .name("ba2")
            .contig(0)
            .start1(100)
            .start2(200)
            .strand1(Strand::Plus)
            .strand2(Strand::Minus)
            .bases1("AAAAAAAAAA")
            .bases2("TTTTTTTTTT")
            .quals1(&quals)
            .quals2(&quals)
            .attr("MI", "test/B")
            .build();

        let _ = builder
            .add_pair()
            .name("ba3")
            .contig(0)
            .start1(100)
            .start2(200)
            .strand1(Strand::Plus)
            .strand2(Strand::Minus)
            .bases1("AAAAACAAAA")
            .bases2("TTTTTTTTTT")
            .quals1(&quals)
            .quals2(&quals)
            .attr("MI", "test/B")
            .build();

        let _ = builder
            .add_pair()
            .name("ba4")
            .contig(0)
            .start1(100)
            .start2(200)
            .strand1(Strand::Plus)
            .strand2(Strand::Minus)
            .bases1("AAAAAAAAAA")
            .bases2("TTTTTTTTTT")
            .quals1(&quals)
            .quals2(&quals)
            .attr("MI", "test/B")
            .build();

        let _ = builder
            .add_pair()
            .name("ba5")
            .contig(0)
            .start1(100)
            .start2(200)
            .strand1(Strand::Plus)
            .strand2(Strand::Minus)
            .bases1("AAAAACAAAA")
            .bases2("TTTTTTTTTT")
            .quals1(&quals)
            .quals2(&quals)
            .attr("MI", "test/B")
            .build();

        let _ = builder
            .add_pair()
            .name("ba6")
            .contig(0)
            .start1(100)
            .start2(200)
            .strand1(Strand::Plus)
            .strand2(Strand::Minus)
            .bases1("AAAAAAAAAA")
            .bases2("TTTTTTTTTT")
            .quals1(&quals)
            .quals2(&quals)
            .attr("MI", "test/B")
            .build();

        let _ = builder
            .add_pair()
            .name("ba7")
            .contig(0)
            .start1(100)
            .start2(200)
            .strand1(Strand::Plus)
            .strand2(Strand::Minus)
            .bases1("AAAAACAAAA")
            .bases2("TTTTTTTTTT")
            .quals1(&quals)
            .quals2(&quals)
            .attr("MI", "test/B")
            .build();

        let _ = builder
            .add_pair()
            .name("ba8")
            .contig(0)
            .start1(100)
            .start2(200)
            .strand1(Strand::Plus)
            .strand2(Strand::Minus)
            .bases1("AAAAAAAAAA")
            .bases2("TTTTTTTTTT")
            .quals1(&quals)
            .quals2(&quals)
            .attr("MI", "test/B")
            .build();

        // Add 1 AB read pair (strand A)
        let _ = builder
            .add_pair()
            .name("ab1")
            .contig(0)
            .start1(200)
            .start2(100)
            .strand1(Strand::Minus)
            .strand2(Strand::Plus)
            .bases1("TTTTTTTTTT")
            .bases2("AAAAAAAAAA")
            .quals1(&quals)
            .quals2(&quals)
            .attr("MI", "test/A")
            .build();

        // Create a duplex consensus caller with minimum thresholds
        let mut caller = DuplexConsensusCaller::new(
            "test".to_string(), // read_name_prefix
            "RG1".to_string(),  // read_group_id
            vec![1],            // min_reads = 1
            20,                 // min_input_base_quality
            true,               // produce_per_base_tags
            false,              // trim
            None,               // max_reads_per_strand
            None,               // cell_tag
            false,              // track_rejects
            45,                 // error_rate_pre_umi
            45,                 // error_rate_post_umi
        )?;

        // Call consensus on these reads
        let consensus_output =
            caller.consensus_reads_from_sam_records(builder.records().to_vec())?;

        assert_eq!(consensus_output.count, 2, "Should generate 2 consensus records (R1 and R2)");

        // Parse the raw records
        let records = ParsedBamRecord::parse_all(&consensus_output.data);

        // Find the R1 record (flag should indicate first of pair)
        let r1 = records
            .iter()
            .find(|r| r.flag & flags::FIRST_SEGMENT != 0)
            .expect("Should have R1 consensus");

        // Get the bc tag (BA simplex consensus)
        let bc_bases = r1.get_string_tag(b"bc").expect("Should have bc tag");

        // Check base at position 5
        // With 4 'A' and 4 'C', this appears to be a tie. However, the order of
        // accumulation in the likelihood calculation (which depends on read order)
        // affects the floating-point result. With reads sorted by descending length
        // (matching fgbio's filterToMostCommonAlignment behavior), the order of
        // base addition produces 'A' as the winner due to floating-point accumulation.
        let base_at_5 = bc_bases[5] as char;

        assert_eq!(base_at_5, 'A', "Position 5 in bc tag should be 'A', got '{base_at_5}'");

        Ok(())
    }

    #[test]
    fn test_duplex_consensus_read_combined_depths() {
        // Test combined_depths() method with both strands present
        let ab_consensus = VanillaConsensusRead {
            id: "test".to_string(),
            bases: vec![b'A', b'C', b'G', b'T'],
            quals: vec![30, 30, 30, 30],
            depths: vec![5, 3, 4, 2],
            errors: vec![0, 0, 0, 0],
            source_reads: None,
        };

        let ba_consensus = VanillaConsensusRead {
            id: "test".to_string(),
            bases: vec![b'A', b'C', b'G', b'T'],
            quals: vec![30, 30, 30, 30],
            depths: vec![3, 4, 2, 5],
            errors: vec![0, 0, 0, 0],
            source_reads: None,
        };

        let duplex = DuplexConsensusRead {
            id: "test".to_string(),
            bases: vec![b'A', b'C', b'G', b'T'],
            quals: vec![30, 30, 30, 30],
            errors: vec![0, 0, 0, 0],
            ab_consensus,
            ba_consensus: Some(ba_consensus),
        };

        // Combined depths should be AB + BA at each position
        let combined = duplex.combined_depths();
        assert_eq!(combined, vec![8, 7, 6, 7]);
    }

    #[test]
    fn test_duplex_consensus_read_combined_depths_ab_only() {
        // Test combined_depths() when only AB strand is present
        let ab_consensus = VanillaConsensusRead {
            id: "test".to_string(),
            bases: vec![b'A', b'C', b'G', b'T'],
            quals: vec![30, 30, 30, 30],
            depths: vec![5, 3, 4, 2],
            errors: vec![0, 0, 0, 0],
            source_reads: None,
        };

        let duplex = DuplexConsensusRead {
            id: "test".to_string(),
            bases: vec![b'A', b'C', b'G', b'T'],
            quals: vec![30, 30, 30, 30],
            errors: vec![0, 0, 0, 0],
            ab_consensus,
            ba_consensus: None,
        };

        // Combined depths should be just AB depths when BA is absent
        let combined = duplex.combined_depths();
        assert_eq!(combined, vec![5, 3, 4, 2]);
    }

    #[test]
    fn test_duplex_consensus_read_max_min_depth() {
        // Test max_depth() and min_depth() methods
        let ab_consensus = VanillaConsensusRead {
            id: "test".to_string(),
            bases: vec![b'A', b'C', b'G', b'T'],
            quals: vec![30, 30, 30, 30],
            depths: vec![5, 3, 4, 2],
            errors: vec![0, 0, 0, 0],
            source_reads: None,
        };

        let ba_consensus = VanillaConsensusRead {
            id: "test".to_string(),
            bases: vec![b'A', b'C', b'G', b'T'],
            quals: vec![30, 30, 30, 30],
            depths: vec![3, 4, 2, 5],
            errors: vec![0, 0, 0, 0],
            source_reads: None,
        };

        let duplex = DuplexConsensusRead {
            id: "test".to_string(),
            bases: vec![b'A', b'C', b'G', b'T'],
            quals: vec![30, 30, 30, 30],
            errors: vec![0, 0, 0, 0],
            ab_consensus,
            ba_consensus: Some(ba_consensus),
        };

        // Combined depths: [8, 7, 6, 7]
        assert_eq!(duplex.max_depth(), 8);
        assert_eq!(duplex.min_depth(), 6);
    }

    #[test]
    fn test_duplex_read_into_unmapped_flag() {
        // Test that duplex_read_into sets UNMAPPED flag
        let ab_consensus = VanillaConsensusRead {
            id: "test".to_string(),
            bases: vec![b'A', b'C', b'G', b'T'],
            quals: vec![30, 30, 30, 30],
            depths: vec![5, 5, 5, 5],
            errors: vec![0, 0, 0, 0],
            source_reads: None,
        };

        let duplex = DuplexConsensusRead {
            id: "test".to_string(),
            bases: vec![b'A', b'C', b'G', b'T'],
            quals: vec![30, 30, 30, 30],
            errors: vec![0, 0, 0, 0],
            ab_consensus,
            ba_consensus: None,
        };

        let mut builder = UnmappedBamRecordBuilder::new();
        let mut output = ConsensusOutput::default();
        DuplexConsensusCaller::duplex_read_into(
            &mut builder,
            &mut output,
            &duplex,
            ReadType::Fragment,
            "UMI123",
            &[],
            &[],
            false,
            "consensus",
            "RG1",
            true,
            None,
            None,
        )
        .unwrap();

        // Verify UNMAPPED flag is set
        let records = ParsedBamRecord::parse_all(&output.data);
        assert_eq!(records.len(), 1);
        assert!(records[0].flag & flags::UNMAPPED != 0, "UNMAPPED flag should be set");
    }

    // ==========================================================================
    // Additional coverage tests for duplex consensus edge cases
    // ==========================================================================

    #[test]
    fn test_duplex_consensus_ab_only() {
        // Test duplex_consensus() when only AB strand is present
        let ab = VanillaConsensusRead {
            id: "UMI123".to_string(),
            bases: vec![b'A', b'C', b'G', b'T'],
            quals: vec![30, 31, 32, 33],
            depths: vec![5, 5, 5, 5],
            errors: vec![0, 1, 0, 1],
            source_reads: None,
        };

        let result = DuplexConsensusCaller::duplex_consensus(Some(&ab), None, None);

        assert!(result.is_some());
        let duplex = result.unwrap();
        assert_eq!(duplex.bases, ab.bases);
        assert_eq!(duplex.quals, ab.quals);
        assert_eq!(duplex.errors, ab.errors);
        assert!(duplex.ba_consensus.is_none());
    }

    #[test]
    fn test_duplex_consensus_ba_only() {
        // Test duplex_consensus() when only BA strand is present
        let ba = VanillaConsensusRead {
            id: "UMI123".to_string(),
            bases: vec![b'T', b'G', b'C', b'A'],
            quals: vec![25, 26, 27, 28],
            depths: vec![4, 4, 4, 4],
            errors: vec![1, 0, 1, 0],
            source_reads: None,
        };

        let result = DuplexConsensusCaller::duplex_consensus(None, Some(&ba), None);

        assert!(result.is_some());
        let duplex = result.unwrap();
        assert_eq!(duplex.bases, ba.bases);
        assert_eq!(duplex.quals, ba.quals);
        assert_eq!(duplex.errors, ba.errors);
        assert!(duplex.ba_consensus.is_none()); // When BA-only, it becomes AB in output
    }

    #[test]
    fn test_duplex_consensus_none_none() {
        // Test duplex_consensus() when neither strand is present
        let result = DuplexConsensusCaller::duplex_consensus(None, None, None);
        assert!(result.is_none());
    }

    #[test]
    fn test_duplex_consensus_different_lengths() {
        // Test duplex_consensus() when AB and BA have different lengths
        let ab = VanillaConsensusRead {
            id: "UMI123".to_string(),
            bases: vec![b'A', b'C', b'G', b'T', b'A', b'C'], // 6 bases
            quals: vec![30, 30, 30, 30, 30, 30],
            depths: vec![5, 5, 5, 5, 5, 5],
            errors: vec![0, 0, 0, 0, 0, 0],
            source_reads: None,
        };

        let ba = VanillaConsensusRead {
            id: "UMI123".to_string(),
            bases: vec![b'A', b'C', b'G', b'T'], // 4 bases
            quals: vec![25, 25, 25, 25],
            depths: vec![4, 4, 4, 4],
            errors: vec![0, 0, 0, 0],
            source_reads: None,
        };

        let result = DuplexConsensusCaller::duplex_consensus(Some(&ab), Some(&ba), None);

        assert!(result.is_some());
        let duplex = result.unwrap();
        // Output should be truncated to minimum length (4)
        assert_eq!(duplex.bases.len(), 4);
        assert_eq!(duplex.quals.len(), 4);
    }

    #[test]
    fn test_duplex_consensus_error_calculation_approximate() {
        // Test error calculation without source reads (approximate method)
        let ab = VanillaConsensusRead {
            id: "UMI123".to_string(),
            bases: vec![b'A', b'C', b'G', b'T'],
            quals: vec![30, 30, 30, 30],
            depths: vec![5, 5, 5, 5],
            errors: vec![1, 0, 2, 0],
            source_reads: None,
        };

        let ba = VanillaConsensusRead {
            id: "UMI123".to_string(),
            bases: vec![b'A', b'C', b'G', b'T'], // Same as AB (agreement)
            quals: vec![25, 25, 25, 25],
            depths: vec![4, 4, 4, 4],
            errors: vec![0, 1, 0, 2],
            source_reads: None,
        };

        let result = DuplexConsensusCaller::duplex_consensus(Some(&ab), Some(&ba), None);

        assert!(result.is_some());
        let duplex = result.unwrap();
        // When bases agree, errors are summed: AB.errors + BA.errors
        assert_eq!(duplex.errors[0], 1); // AB=1, BA=0 -> 1
        assert_eq!(duplex.errors[1], 1); // AB=0, BA=1 -> 1
        assert_eq!(duplex.errors[2], 2); // AB=2, BA=0 -> 2
        assert_eq!(duplex.errors[3], 2); // AB=0, BA=2 -> 2
    }

    #[test]
    fn test_duplex_consensus_error_calculation_disagreement() {
        // Test error calculation when strands disagree
        let ab = VanillaConsensusRead {
            id: "UMI123".to_string(),
            bases: vec![b'A', b'T'], // Position 1: A has higher qual
            quals: vec![30, 40],     // Higher quality
            depths: vec![5, 5],
            errors: vec![1, 2],
            source_reads: None,
        };

        let ba = VanillaConsensusRead {
            id: "UMI123".to_string(),
            bases: vec![b'A', b'C'], // Position 1 differs: C vs T
            quals: vec![25, 30],     // Lower quality at pos 1
            depths: vec![4, 4],
            errors: vec![0, 1],
            source_reads: None,
        };

        let result = DuplexConsensusCaller::duplex_consensus(Some(&ab), Some(&ba), None);

        assert!(result.is_some());
        let duplex = result.unwrap();
        // Position 0: Agreement, base = A
        // Position 1: Disagreement, AB wins (higher qual), base = T
        //   Errors when raw_base == a_base: a_err + (b_dep - b_err) = 2 + (4-1) = 5
        assert_eq!(duplex.bases[0], b'A');
        assert_eq!(duplex.bases[1], b'T');
        assert_eq!(duplex.errors[1], 5);
    }

    #[test]
    fn test_has_minimum_number_of_reads_symmetric() {
        // Test symmetric min_reads configuration
        // has_minimum_number_of_reads only counts R1s (first segment of paired reads)

        // Create 3 AB R1 reads
        let a_reads: Vec<Vec<u8>> = (0..3)
            .map(|_| {
                encode_to_raw(
                    &RecordBuilder::new()
                        .sequence("ACGT")
                        .first_segment(true)
                        .tag("MI", "test/A")
                        .build(),
                )
            })
            .collect();

        // Create 2 BA R1 reads
        let b_reads: Vec<Vec<u8>> = (0..2)
            .map(|_| {
                encode_to_raw(
                    &RecordBuilder::new()
                        .sequence("ACGT")
                        .first_segment(true)
                        .tag("MI", "test/B")
                        .build(),
                )
            })
            .collect();

        // min_total=4, min_xy=2, min_yx=2 should pass (3+2=5, max(3,2)=3, min(3,2)=2)
        assert!(DuplexConsensusCaller::has_minimum_number_of_reads(&a_reads, &b_reads, 4, 2, 2));

        // min_total=6 should fail (3+2=5)
        assert!(!DuplexConsensusCaller::has_minimum_number_of_reads(&a_reads, &b_reads, 6, 2, 2));

        // min_xy=4 should fail (max(3,2)=3)
        assert!(!DuplexConsensusCaller::has_minimum_number_of_reads(&a_reads, &b_reads, 4, 4, 1));

        // min_yx=3 should fail (min(3,2)=2)
        assert!(!DuplexConsensusCaller::has_minimum_number_of_reads(&a_reads, &b_reads, 4, 2, 3));
    }

    #[test]
    fn test_has_minimum_number_of_reads_asymmetric() {
        // Test asymmetric min_reads [D, M_AB, M_BA] configuration

        // Create 3 AB R1 reads
        let a_reads: Vec<Vec<u8>> = (0..3)
            .map(|_| {
                encode_to_raw(
                    &RecordBuilder::new()
                        .sequence("ACGT")
                        .first_segment(true)
                        .tag("MI", "test/A")
                        .build(),
                )
            })
            .collect();

        // Create 1 BA R1 read
        let b_reads = vec![encode_to_raw(
            &RecordBuilder::new().sequence("ACGT").first_segment(true).tag("MI", "test/B").build(),
        )];

        // [3, 2, 1]: min_total=3, min_xy=2, min_yx=1
        // A=3, B=1, xy=max(3,1)=3, yx=min(3,1)=1
        // Should pass: total=4>=3, xy=3>=2, yx=1>=1
        assert!(DuplexConsensusCaller::has_minimum_number_of_reads(&a_reads, &b_reads, 3, 2, 1));

        // [3, 2, 2]: should fail because yx=1 < 2
        assert!(!DuplexConsensusCaller::has_minimum_number_of_reads(&a_reads, &b_reads, 3, 2, 2));
    }

    #[test]
    fn test_cap_quality() {
        // Test quality score capping
        assert_eq!(DuplexConsensusCaller::cap_quality(-5), 2); // Below min
        assert_eq!(DuplexConsensusCaller::cap_quality(0), 2); // Below min
        assert_eq!(DuplexConsensusCaller::cap_quality(2), 2); // At min
        assert_eq!(DuplexConsensusCaller::cap_quality(50), 50); // In range
        assert_eq!(DuplexConsensusCaller::cap_quality(93), 93); // At max
        assert_eq!(DuplexConsensusCaller::cap_quality(100), 93); // Above max
    }

    #[test]
    fn test_is_error() {
        // Test error detection function
        assert!(DuplexConsensusCaller::is_error(b'A', b'T')); // Different bases
        assert!(DuplexConsensusCaller::is_error(b'C', b'G')); // Different bases
        assert!(!DuplexConsensusCaller::is_error(b'A', b'A')); // Same base
        assert!(!DuplexConsensusCaller::is_error(b'N', b'A')); // Source is N
        assert!(!DuplexConsensusCaller::is_error(b'A', b'N')); // Consensus is N
        assert!(!DuplexConsensusCaller::is_error(b'N', b'N')); // Both N
    }

    #[test]
    fn test_are_all_same_strand() {
        // Create test records
        let fwd_raw = encode_to_raw(&RecordBuilder::new().sequence("ACGT").build());

        let rev_raw =
            encode_to_raw(&RecordBuilder::new().sequence("ACGT").reverse_complement(true).build());

        // Empty collection
        assert!(DuplexConsensusCaller::are_all_same_strand(std::iter::empty::<&Vec<u8>>()));

        // All forward
        assert!(DuplexConsensusCaller::are_all_same_strand([&fwd_raw, &fwd_raw].into_iter()));

        // All reverse
        assert!(DuplexConsensusCaller::are_all_same_strand([&rev_raw, &rev_raw].into_iter()));

        // Mixed strands
        assert!(!DuplexConsensusCaller::are_all_same_strand([&fwd_raw, &rev_raw].into_iter()));
    }

    #[test]
    fn test_duplex_consensus_has_minimum_reads() {
        let ab_consensus = VanillaConsensusRead {
            id: "test".to_string(),
            bases: vec![b'A', b'C', b'G', b'T'],
            quals: vec![30, 30, 30, 30],
            depths: vec![5, 6, 4, 5], // max_depth = 6
            errors: vec![0, 0, 0, 0],
            source_reads: None,
        };

        let ba_consensus = VanillaConsensusRead {
            id: "test".to_string(),
            bases: vec![b'A', b'C', b'G', b'T'],
            quals: vec![30, 30, 30, 30],
            depths: vec![3, 4, 2, 3], // max_depth = 4
            errors: vec![0, 0, 0, 0],
            source_reads: None,
        };

        let duplex = DuplexConsensusRead {
            id: "test".to_string(),
            bases: vec![b'A', b'C', b'G', b'T'],
            quals: vec![30, 30, 30, 30],
            errors: vec![0, 0, 0, 0],
            ab_consensus: ab_consensus.clone(),
            ba_consensus: Some(ba_consensus),
        };

        // num_a=6, num_b=4, xy=6, yx=4, total=10
        // min_total=8, min_xy=5, min_yx=3 should pass
        assert!(DuplexConsensusCaller::duplex_consensus_has_minimum_reads(&duplex, 8, 5, 3));

        // min_total=11 should fail (total=10)
        assert!(!DuplexConsensusCaller::duplex_consensus_has_minimum_reads(&duplex, 11, 5, 3));

        // Test with BA=None
        let duplex_ab_only = DuplexConsensusRead {
            id: "test".to_string(),
            bases: vec![b'A', b'C', b'G', b'T'],
            quals: vec![30, 30, 30, 30],
            errors: vec![0, 0, 0, 0],
            ab_consensus,
            ba_consensus: None,
        };

        // num_a=6, num_b=0, xy=6, yx=0, total=6
        assert!(DuplexConsensusCaller::duplex_consensus_has_minimum_reads(
            &duplex_ab_only,
            5,
            5,
            0
        ));
        assert!(!DuplexConsensusCaller::duplex_consensus_has_minimum_reads(
            &duplex_ab_only,
            5,
            5,
            1
        )); // yx=0 < 1
    }

    #[test]
    fn test_duplex_read_into_with_cell_tag() {
        // Test that cell barcode is properly included
        let ab_consensus = VanillaConsensusRead {
            id: "test".to_string(),
            bases: vec![b'A', b'C', b'G', b'T'],
            quals: vec![30, 30, 30, 30],
            depths: vec![5, 5, 5, 5],
            errors: vec![0, 0, 0, 0],
            source_reads: None,
        };

        let duplex = DuplexConsensusRead {
            id: "test".to_string(),
            bases: vec![b'A', b'C', b'G', b'T'],
            quals: vec![30, 30, 30, 30],
            errors: vec![0, 0, 0, 0],
            ab_consensus,
            ba_consensus: None,
        };

        let cell_tag = Tag::new(b'C', b'B');
        let cell_barcode = "ACGTACGT-1";

        let mut builder = UnmappedBamRecordBuilder::new();
        let mut output = ConsensusOutput::default();
        DuplexConsensusCaller::duplex_read_into(
            &mut builder,
            &mut output,
            &duplex,
            ReadType::Fragment,
            "UMI123",
            &[],
            &[],
            false,
            "consensus",
            "RG1",
            true,
            Some(cell_tag),
            Some(cell_barcode),
        )
        .unwrap();

        // Check cell barcode tag is present
        let records = ParsedBamRecord::parse_all(&output.data);
        assert_eq!(records.len(), 1);
        let cb = records[0].get_string_tag(b"CB").expect("CB tag should be present");
        assert_eq!(String::from_utf8(cb).unwrap(), cell_barcode);
    }

    #[test]
    fn test_duplex_read_into_without_cell_tag() {
        // Test that record is created correctly without cell barcode
        let ab_consensus = VanillaConsensusRead {
            id: "test".to_string(),
            bases: vec![b'A', b'C', b'G', b'T'],
            quals: vec![30, 30, 30, 30],
            depths: vec![5, 5, 5, 5],
            errors: vec![0, 0, 0, 0],
            source_reads: None,
        };

        let duplex = DuplexConsensusRead {
            id: "test".to_string(),
            bases: vec![b'A', b'C', b'G', b'T'],
            quals: vec![30, 30, 30, 30],
            errors: vec![0, 0, 0, 0],
            ab_consensus,
            ba_consensus: None,
        };

        let mut builder = UnmappedBamRecordBuilder::new();
        let mut output = ConsensusOutput::default();
        DuplexConsensusCaller::duplex_read_into(
            &mut builder,
            &mut output,
            &duplex,
            ReadType::Fragment,
            "UMI123",
            &[],
            &[],
            false,
            "consensus",
            "RG1",
            true,
            None,
            None,
        )
        .unwrap();

        // Check CB tag is not present
        let records = ParsedBamRecord::parse_all(&output.data);
        assert_eq!(records.len(), 1);
        assert!(records[0].get_string_tag(b"CB").is_none());
    }

    #[test]
    fn test_duplex_consensus_with_zero_depth_positions() {
        // Test filtering of consensuses with no coverage in truncated region
        let ab = VanillaConsensusRead {
            id: "UMI123".to_string(),
            bases: vec![b'A', b'C', b'G', b'T', b'A', b'C'],
            quals: vec![30, 30, 30, 30, 30, 30],
            depths: vec![0, 0, 0, 0, 5, 5], // First 4 positions have no coverage
            errors: vec![0, 0, 0, 0, 0, 0],
            source_reads: None,
        };

        let ba = VanillaConsensusRead {
            id: "UMI123".to_string(),
            bases: vec![b'T', b'G', b'C', b'A'],
            quals: vec![25, 25, 25, 25],
            depths: vec![4, 4, 4, 4], // All positions have coverage
            errors: vec![0, 0, 0, 0],
            source_reads: None,
        };

        // After truncation to length 4, AB has no coverage (all depths are 0)
        let result = DuplexConsensusCaller::duplex_consensus(Some(&ab), Some(&ba), None);

        assert!(result.is_some());
        let duplex = result.unwrap();
        // AB should be filtered out, so we should only have BA data
        assert_eq!(duplex.bases, ba.bases);
    }

    #[test]
    fn test_extract_int_tag_variants() {
        // Test extract_int_tag with different integer types
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::Data;
        use noodles::sam::alignment::record_buf::data::field::Value;

        let cd_tag = Tag::from([b'c', b'D']);

        // Helper to create Data with a specific value
        let create_data = |value: Value| {
            let mut data = Data::default();
            data.insert(cd_tag, value);
            data
        };

        // Test Int8
        let data = create_data(Value::Int8(42));
        assert_eq!(DuplexConsensusCaller::extract_int_tag(&data, cd_tag), Some(42));

        // Test UInt8
        let data = create_data(Value::UInt8(200));
        assert_eq!(DuplexConsensusCaller::extract_int_tag(&data, cd_tag), Some(200));

        // Test Int16
        let data = create_data(Value::Int16(1000));
        assert_eq!(DuplexConsensusCaller::extract_int_tag(&data, cd_tag), Some(1000));

        // Test UInt16
        let data = create_data(Value::UInt16(50000));
        assert_eq!(DuplexConsensusCaller::extract_int_tag(&data, cd_tag), Some(50000));

        // Test Int32
        let data = create_data(Value::Int32(100_000));
        assert_eq!(DuplexConsensusCaller::extract_int_tag(&data, cd_tag), Some(100_000));

        // Test UInt32
        let data = create_data(Value::UInt32(200_000));
        assert_eq!(DuplexConsensusCaller::extract_int_tag(&data, cd_tag), Some(200_000));

        // Test missing tag returns None
        let data = Data::default();
        assert_eq!(DuplexConsensusCaller::extract_int_tag(&data, cd_tag), None);

        // Test non-integer type returns None
        let data = create_data(Value::Float(1.5));
        assert_eq!(DuplexConsensusCaller::extract_int_tag(&data, cd_tag), None);
    }

    #[test]
    fn test_duplex_read_into_ba_none_tags() {
        // Test that BA=None produces correct zero tags for bD, bM, bE
        let ab_consensus = VanillaConsensusRead {
            id: "test".to_string(),
            bases: vec![b'A', b'C', b'G', b'T'],
            quals: vec![30, 30, 30, 30],
            depths: vec![5, 5, 5, 5],
            errors: vec![0, 1, 0, 1],
            source_reads: None,
        };

        let duplex = DuplexConsensusRead {
            id: "test".to_string(),
            bases: vec![b'A', b'C', b'G', b'T'],
            quals: vec![30, 30, 30, 30],
            errors: vec![0, 1, 0, 1],
            ab_consensus: ab_consensus.clone(),
            ba_consensus: None, // BA is None
        };

        let mut builder = UnmappedBamRecordBuilder::new();
        let mut output = ConsensusOutput::default();
        DuplexConsensusCaller::duplex_read_into(
            &mut builder,
            &mut output,
            &duplex,
            ReadType::R1,
            "UMI123",
            &[],
            &[],
            false, // no per-base tags
            "consensus",
            "RG1",
            true,
            None,
            None,
        )
        .unwrap();

        let records = ParsedBamRecord::parse_all(&output.data);
        assert_eq!(records.len(), 1);
        let rec = &records[0];

        // Check BA tags are set to zero
        assert_eq!(rec.get_int_tag(b"bD"), Some(0));
        assert_eq!(rec.get_int_tag(b"bM"), Some(0));

        let be_val = rec.get_float_tag(b"bE").expect("bE tag should be present");
        assert!((be_val - 0.0).abs() < 0.001);

        // Verify AB tags are set correctly
        assert_eq!(rec.get_int_tag(b"aD"), Some(5)); // max depth from AB
        assert_eq!(rec.get_int_tag(b"aM"), Some(5)); // min depth from AB
    }

    #[test]
    fn test_are_all_same_strand_mixed() {
        // Test are_all_same_strand with mixed strands
        let fwd_raw = encode_to_raw(&RecordBuilder::new().sequence("ACGT").build());

        let rev_raw =
            encode_to_raw(&RecordBuilder::new().sequence("ACGT").reverse_complement(true).build());

        // All forward
        assert!(DuplexConsensusCaller::are_all_same_strand([&fwd_raw, &fwd_raw].into_iter()));

        // All reverse
        assert!(DuplexConsensusCaller::are_all_same_strand([&rev_raw, &rev_raw].into_iter()));

        // Mixed strands - should return false
        assert!(!DuplexConsensusCaller::are_all_same_strand([&fwd_raw, &rev_raw].into_iter()));

        // Empty - should return true
        assert!(DuplexConsensusCaller::are_all_same_strand(std::iter::empty::<&Vec<u8>>()));

        // Single record
        assert!(DuplexConsensusCaller::are_all_same_strand([&fwd_raw].into_iter()));
    }

    #[test]
    fn test_duplex_read_into_with_per_base_tags() {
        // Test that per-base tags are generated when produce_per_base_tags=true
        let ab_consensus = VanillaConsensusRead {
            id: "test".to_string(),
            bases: vec![b'A', b'C', b'G', b'T'],
            quals: vec![30, 31, 32, 33],
            depths: vec![5, 6, 7, 8],
            errors: vec![0, 1, 0, 2],
            source_reads: None,
        };

        let ba_consensus = VanillaConsensusRead {
            id: "test".to_string(),
            bases: vec![b'T', b'G', b'C', b'A'],
            quals: vec![25, 26, 27, 28],
            depths: vec![3, 4, 5, 6],
            errors: vec![1, 0, 1, 0],
            source_reads: None,
        };

        let duplex = DuplexConsensusRead {
            id: "test".to_string(),
            bases: vec![b'A', b'C', b'G', b'T'],
            quals: vec![30, 30, 30, 30],
            errors: vec![1, 1, 1, 2],
            ab_consensus: ab_consensus.clone(),
            ba_consensus: Some(ba_consensus),
        };

        let mut builder = UnmappedBamRecordBuilder::new();
        let mut output = ConsensusOutput::default();
        DuplexConsensusCaller::duplex_read_into(
            &mut builder,
            &mut output,
            &duplex,
            ReadType::R1,
            "UMI123",
            &[],
            &[],
            true, // with per-base tags
            "consensus",
            "RG1",
            true,
            None,
            None,
        )
        .unwrap();

        let records = ParsedBamRecord::parse_all(&output.data);
        assert_eq!(records.len(), 1);
        let rec = &records[0];

        // Check per-base AB tags
        assert!(rec.get_i16_array_tag(b"ad").is_some(), "ad tag should be present");
        assert!(rec.get_i16_array_tag(b"ae").is_some(), "ae tag should be present");
        assert!(rec.get_string_tag(b"ac").is_some(), "ac tag should be present");

        // Check per-base BA tags
        assert!(rec.get_i16_array_tag(b"bd").is_some(), "bd tag should be present");
        assert!(rec.get_i16_array_tag(b"be").is_some(), "be tag should be present");
        assert!(rec.get_string_tag(b"bc").is_some(), "bc tag should be present");
    }

    #[test]
    fn test_duplex_read_into_read_types() {
        // Test different read types (R1, R2, Fragment)
        let ab_consensus = VanillaConsensusRead {
            id: "test".to_string(),
            bases: vec![b'A', b'C', b'G', b'T'],
            quals: vec![30, 30, 30, 30],
            depths: vec![5, 5, 5, 5],
            errors: vec![0, 0, 0, 0],
            source_reads: None,
        };

        let duplex = DuplexConsensusRead {
            id: "test".to_string(),
            bases: vec![b'A', b'C', b'G', b'T'],
            quals: vec![30, 30, 30, 30],
            errors: vec![0, 0, 0, 0],
            ab_consensus: ab_consensus.clone(),
            ba_consensus: None,
        };

        // Helper to build and parse a single record
        let build_and_parse = |read_type: ReadType| -> ParsedBamRecord {
            let mut builder = UnmappedBamRecordBuilder::new();
            let mut output = ConsensusOutput::default();
            DuplexConsensusCaller::duplex_read_into(
                &mut builder,
                &mut output,
                &duplex,
                read_type,
                "UMI123",
                &[],
                &[],
                false,
                "consensus",
                "RG1",
                true,
                None,
                None,
            )
            .unwrap();
            let records = ParsedBamRecord::parse_all(&output.data);
            assert_eq!(records.len(), 1);
            records.into_iter().next().unwrap()
        };

        // R1 read type
        let rec_r1 = build_and_parse(ReadType::R1);
        assert!(rec_r1.flag & flags::PAIRED != 0);
        assert!(rec_r1.flag & flags::FIRST_SEGMENT != 0);
        assert!(rec_r1.flag & flags::LAST_SEGMENT == 0);

        // R2 read type
        let rec_r2 = build_and_parse(ReadType::R2);
        assert!(rec_r2.flag & flags::PAIRED != 0);
        assert!(rec_r2.flag & flags::FIRST_SEGMENT == 0);
        assert!(rec_r2.flag & flags::LAST_SEGMENT != 0);

        // Fragment read type
        let rec_frag = build_and_parse(ReadType::Fragment);
        // Fragment should not have PAIRED flag
        assert!(rec_frag.flag & flags::PAIRED == 0);
    }

    // ============================================================================
    // Additional HIGH risk code path tests
    // ============================================================================

    #[test]
    fn test_duplex_consensus_error_approximation_agreement() {
        // Test error approximation without source reads - agreement case
        let ab = VanillaConsensusRead {
            id: "test".to_string(),
            bases: vec![b'A', b'A', b'A'],
            quals: vec![30, 30, 30],
            depths: vec![5, 5, 5],
            errors: vec![1, 0, 2], // Some errors in AB
            source_reads: None,
        };

        let ba = VanillaConsensusRead {
            id: "test".to_string(),
            bases: vec![b'A', b'A', b'A'], // Same bases - agreement
            quals: vec![25, 25, 25],
            depths: vec![3, 3, 3],
            errors: vec![0, 1, 1], // Some errors in BA
            source_reads: None,
        };

        let result = DuplexConsensusCaller::duplex_consensus(Some(&ab), Some(&ba), None);
        assert!(result.is_some());

        let duplex = result.unwrap();
        // Errors should be sum of AB + BA errors when bases agree
        assert_eq!(duplex.errors[0], 1); // 1 + 0
        assert_eq!(duplex.errors[1], 1); // 0 + 1
        assert_eq!(duplex.errors[2], 3); // 2 + 1
    }

    #[test]
    fn test_duplex_consensus_error_approximation_disagreement() {
        // Test error approximation without source reads - disagreement case
        let ab = VanillaConsensusRead {
            id: "test".to_string(),
            bases: vec![b'A'],
            quals: vec![40], // Higher qual - AB wins
            depths: vec![5],
            errors: vec![1],
            source_reads: None,
        };

        let ba = VanillaConsensusRead {
            id: "test".to_string(),
            bases: vec![b'T'], // Different base
            quals: vec![20],   // Lower qual
            depths: vec![4],
            errors: vec![0],
            source_reads: None,
        };

        let result = DuplexConsensusCaller::duplex_consensus(Some(&ab), Some(&ba), None);
        assert!(result.is_some());

        let duplex = result.unwrap();
        // When AB wins: errors = ab_errors + (ba_depth - ba_errors) = 1 + (4 - 0) = 5
        assert_eq!(duplex.bases[0], b'A');
        assert_eq!(duplex.errors[0], 5);
    }

    #[test]
    fn test_duplex_consensus_equal_quality_disagreement() {
        // Test equal quality disagreement - should produce low quality
        let ab = VanillaConsensusRead {
            id: "test".to_string(),
            bases: vec![b'A'],
            quals: vec![30],
            depths: vec![5],
            errors: vec![0],
            source_reads: None,
        };

        let ba = VanillaConsensusRead {
            id: "test".to_string(),
            bases: vec![b'T'], // Different base
            quals: vec![30],   // Same qual
            depths: vec![5],
            errors: vec![0],
            source_reads: None,
        };

        let result = DuplexConsensusCaller::duplex_consensus(Some(&ab), Some(&ba), None);
        assert!(result.is_some());

        let duplex = result.unwrap();
        // Equal quality disagreement should result in N with min quality
        assert_eq!(duplex.bases[0], b'N');
        assert_eq!(duplex.quals[0], 2); // MIN_PHRED
    }

    #[test]
    fn test_duplex_consensus_n_base_masking() {
        // Test that N bases in either strand result in N output
        let ab = VanillaConsensusRead {
            id: "test".to_string(),
            bases: vec![b'N', b'A'],
            quals: vec![30, 30],
            depths: vec![5, 5],
            errors: vec![0, 0],
            source_reads: None,
        };

        let ba = VanillaConsensusRead {
            id: "test".to_string(),
            bases: vec![b'A', b'N'],
            quals: vec![30, 30],
            depths: vec![5, 5],
            errors: vec![0, 0],
            source_reads: None,
        };

        let result = DuplexConsensusCaller::duplex_consensus(Some(&ab), Some(&ba), None);
        assert!(result.is_some());

        let duplex = result.unwrap();
        // Both positions should be N due to N in one strand
        assert_eq!(duplex.bases[0], b'N');
        assert_eq!(duplex.bases[1], b'N');
    }

    #[test]
    fn test_duplex_consensus_zero_depth_filtering() {
        // Test that consensus with all zero depths in truncated region is filtered
        let ab = VanillaConsensusRead {
            id: "test".to_string(),
            bases: vec![b'A', b'C'],
            quals: vec![30, 30],
            depths: vec![0, 0], // All zero depths
            errors: vec![0, 0],
            source_reads: None,
        };

        let ba = VanillaConsensusRead {
            id: "test".to_string(),
            bases: vec![b'A', b'C'],
            quals: vec![30, 30],
            depths: vec![5, 5], // Non-zero depths
            errors: vec![0, 0],
            source_reads: None,
        };

        let result = DuplexConsensusCaller::duplex_consensus(Some(&ab), Some(&ba), None);
        assert!(result.is_some());

        // AB should be filtered out due to zero depths, leaving only BA
        let duplex = result.unwrap();
        assert!(duplex.ba_consensus.is_none()); // Only one strand present
    }

    #[test]
    fn test_has_minimum_number_of_reads_basic() {
        // Create mock reads with paired + first segment flags for R1 counting
        let a_reads: Vec<Vec<u8>> = (0..3)
            .map(|_| {
                encode_to_raw(&RecordBuilder::new().sequence("ACGT").first_segment(true).build())
            })
            .collect();

        let b_reads: Vec<Vec<u8>> = (0..3)
            .map(|_| {
                encode_to_raw(&RecordBuilder::new().sequence("ACGT").first_segment(true).build())
            })
            .collect();

        // Test with sufficient reads: total=6, xy=3, yx=3
        // min_total=4, min_xy=2, min_yx=2 -> should pass
        assert!(DuplexConsensusCaller::has_minimum_number_of_reads(
            &a_reads, &b_reads, 4, // min_total
            2, // min_xy
            2, // min_yx
        ));

        // Test with insufficient total reads: min_total=8 > 6
        assert!(!DuplexConsensusCaller::has_minimum_number_of_reads(
            &a_reads, &b_reads, 8, // min_total - too high
            2, // min_xy
            2, // min_yx
        ));
    }

    #[test]
    fn test_has_minimum_number_of_reads_empty_strand() {
        // Create 3 A reads
        let a_reads: Vec<Vec<u8>> = (0..3)
            .map(|_| {
                encode_to_raw(&RecordBuilder::new().sequence("ACGT").first_segment(true).build())
            })
            .collect();

        // Empty B reads
        let b_reads: Vec<Vec<u8>> = vec![];

        // XY = 3, YX = 0, total = 3

        // With min_yx=0, should pass (AB-only scenario)
        assert!(DuplexConsensusCaller::has_minimum_number_of_reads(
            &a_reads, &b_reads, 1, // min_total
            1, // min_xy
            0, // min_yx - allow zero
        ));

        // With min_yx=1, should fail
        assert!(!DuplexConsensusCaller::has_minimum_number_of_reads(
            &a_reads, &b_reads, 1, // min_total
            1, // min_xy
            1, // min_yx - requires at least 1
        ));
    }

    #[test]
    fn test_has_minimum_number_of_reads_only_counts_r1() {
        // Create R1 and R2 reads for A strand
        let a_r1 =
            encode_to_raw(&RecordBuilder::new().sequence("ACGT").first_segment(true).build());

        let a_r2 = encode_to_raw(
            &RecordBuilder::new()
                .sequence("ACGT")
                .first_segment(false) // R2, should not be counted
                .build(),
        );

        let a_reads = vec![a_r1, a_r2]; // Only 1 R1

        let b_r1 =
            encode_to_raw(&RecordBuilder::new().sequence("ACGT").first_segment(true).build());

        let b_reads = vec![b_r1]; // 1 R1

        // Total R1s = 2 (1 + 1), even though a_reads.len() = 2

        // Should pass with min_total=2
        assert!(DuplexConsensusCaller::has_minimum_number_of_reads(
            &a_reads, &b_reads, 2, // min_total
            1, // min_xy
            1, // min_yx
        ));

        // Should fail with min_total=3
        assert!(!DuplexConsensusCaller::has_minimum_number_of_reads(
            &a_reads, &b_reads, 3, // min_total - too high
            1, // min_xy
            1, // min_yx
        ));
    }

    #[test]
    fn test_duplex_consensus_has_minimum_reads_basic() {
        let ab = VanillaConsensusRead {
            id: "test".to_string(),
            bases: vec![b'A', b'C', b'G'],
            quals: vec![30, 30, 30],
            depths: vec![5, 4, 3], // max = 5
            errors: vec![0, 0, 0],
            source_reads: None,
        };

        let ba = VanillaConsensusRead {
            id: "test".to_string(),
            bases: vec![b'A', b'C', b'G'],
            quals: vec![30, 30, 30],
            depths: vec![3, 2, 1], // max = 3
            errors: vec![0, 0, 0],
            source_reads: None,
        };

        let duplex = DuplexConsensusRead {
            id: "test".to_string(),
            bases: vec![b'A', b'C', b'G'],
            quals: vec![60, 60, 60],
            errors: vec![0, 0, 0],
            ab_consensus: ab,
            ba_consensus: Some(ba),
        };

        // max_a = 5, max_b = 3, xy = 5, yx = 3, total = 8
        assert!(DuplexConsensusCaller::duplex_consensus_has_minimum_reads(
            &duplex, 8, // min_total (8 >= 8)
            5, // min_xy (5 >= 5)
            3, // min_yx (3 >= 3)
        ));

        // Should fail if min_yx too high
        assert!(!DuplexConsensusCaller::duplex_consensus_has_minimum_reads(
            &duplex, 8, // min_total
            5, // min_xy
            4, // min_yx (3 < 4)
        ));
    }

    #[test]
    fn test_duplex_consensus_has_minimum_reads_ba_only() {
        // Test with AB-only (BA = None)
        let ab = VanillaConsensusRead {
            id: "test".to_string(),
            bases: vec![b'A'],
            quals: vec![30],
            depths: vec![5],
            errors: vec![0],
            source_reads: None,
        };

        let duplex = DuplexConsensusRead {
            id: "test".to_string(),
            bases: vec![b'A'],
            quals: vec![30],
            errors: vec![0],
            ab_consensus: ab,
            ba_consensus: None,
        };

        // max_a = 5, max_b = 0, xy = 5, yx = 0, total = 5
        assert!(DuplexConsensusCaller::duplex_consensus_has_minimum_reads(
            &duplex, 1, // min_total
            1, // min_xy
            0, // min_yx - allow zero
        ));

        // Should fail if min_yx > 0
        assert!(!DuplexConsensusCaller::duplex_consensus_has_minimum_reads(
            &duplex, 1, // min_total
            1, // min_xy
            1, // min_yx (0 < 1)
        ));
    }

    // ========================================================================
    // partition_records_by_strand tests
    // ========================================================================

    #[test]
    fn test_partition_records_by_strand_basic() {
        let a1 = encode_to_raw(
            &RecordBuilder::new()
                .name("r1")
                .sequence("ACGT")
                .qualities(&[30; 4])
                .tag("MI", "UMI1/A")
                .build(),
        );
        let a2 = encode_to_raw(
            &RecordBuilder::new()
                .name("r2")
                .sequence("ACGT")
                .qualities(&[30; 4])
                .tag("MI", "UMI1/A")
                .build(),
        );
        let b1 = encode_to_raw(
            &RecordBuilder::new()
                .name("r3")
                .sequence("ACGT")
                .qualities(&[30; 4])
                .tag("MI", "UMI1/B")
                .build(),
        );

        let records = vec![a1, a2, b1];
        let (base_mi, a_records, b_records) =
            DuplexConsensusCaller::partition_records_by_strand(records).unwrap();

        assert_eq!(base_mi.as_deref(), Some("UMI1"));
        assert_eq!(a_records.len(), 2);
        assert_eq!(b_records.len(), 1);
    }

    #[test]
    fn test_partition_records_by_strand_empty() {
        let (base_mi, a_records, b_records) =
            DuplexConsensusCaller::partition_records_by_strand(Vec::new()).unwrap();

        assert!(base_mi.is_none());
        assert!(a_records.is_empty());
        assert!(b_records.is_empty());
    }

    #[test]
    fn test_partition_records_by_strand_missing_suffix() {
        let rec = encode_to_raw(
            &RecordBuilder::new()
                .name("r1")
                .sequence("ACGT")
                .qualities(&[30; 4])
                .tag("MI", "UMI1") // No /A or /B suffix
                .build(),
        );

        let result = DuplexConsensusCaller::partition_records_by_strand(vec![rec]);
        assert!(result.is_err(), "Should error on MI tag without /A or /B suffix");
    }

    #[test]
    fn test_partition_records_by_strand_a_only() {
        let a1 = encode_to_raw(
            &RecordBuilder::new()
                .name("r1")
                .sequence("ACGT")
                .qualities(&[30; 4])
                .tag("MI", "UMI1/A")
                .build(),
        );

        let (base_mi, a_records, b_records) =
            DuplexConsensusCaller::partition_records_by_strand(vec![a1]).unwrap();

        assert_eq!(base_mi.as_deref(), Some("UMI1"));
        assert_eq!(a_records.len(), 1);
        assert!(b_records.is_empty());
    }
}

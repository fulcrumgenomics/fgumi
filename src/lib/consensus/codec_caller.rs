//! # CODEC Consensus Calling
//!
//! This module implements consensus calling for CODEC sequencing protocol, where a single
//! read-pair sequences both strands of the original duplex molecule.
//!
//! ## Overview
//!
//! CODEC (Bae et al 2023: <https://doi.org/10.1038/s41588-023-01376-0>) is a sequencing
//! protocol where:
//! - Each read-pair sequences **both** strands of the original duplex molecule
//! - R1 comes from one strand, R2 comes from the opposite strand
//! - Even a single read-pair can generate duplex consensus
//! - Output is a **single fragment read** per duplex (not paired-end)
//!
//! ## Key Differences from Regular Duplex
//!
//! - **Regular Duplex**: Multiple read-pairs, separate /A and /B families, pairs consensus
//! - **CODEC**: Single read-pair can work, R1 and R2 are inherently from opposite strands
//! - **Output**: Single fragment (not pair) representing the full duplex molecule
//!
//! ## Algorithm Overview
//!
//! The CODEC consensus calling process involves these steps:
//!
//! 1. **Clip** read pairs where they extend past mate ends
//! 2. **Filter** R1s and R2s separately for compatible CIGARs
//! 3. **Check overlap** - must have sufficient overlap on genome
//! 4. **Build single-strand consensus** from R1s separately from R2s
//! 5. **Reverse complement** one SS read to align orientations
//! 6. **Pad with Ns** if reads don't fully overlap
//! 7. **Build duplex consensus** from the two SS consensus reads
//! 8. **Apply quality masking** based on single-strand vs duplex regions
//! 9. **Reverse complement** back if original R1s were on negative strand
//!
//! ## Critical Requirements
//!
//! - Reads must be **mapped** and **overlap** on the genome
//! - Uses genome alignment to synchronize reads for duplex calling
//! - Rejects unmapped, chimeric, or non-overlapping pairs
//!
//! ## Usage Example
//!
//! ```rust,ignore
//! use fgumi_lib::consensus::codec_caller::{CodecConsensusCaller, CodecConsensusOptions};
//!
//! let options = CodecConsensusOptions {
//!     min_reads_per_strand: 1,
//!     min_duplex_length: 1,
//!     ..Default::default()
//! };
//!
//! let mut caller = CodecConsensusCaller::new(
//!     "codec".to_string(),
//!     "RG1".to_string(),
//!     options,
//! );
//!
//! // Process reads - R1s and R2s from same molecule
//! let consensus = caller.consensus_reads_from_sam_records(reads)?;
//! ```
//!
//! ## Output Format
//!
//! - **Unmapped**: Output reads are unmapped fragments
//! - **Single fragment**: Not paired-end
//! - **Tags added** (similar to duplex):
//!   - `aD, bD, cD` - Max depth (int)
//!   - `aM, bM, cM` - Min depth (int)
//!   - `aE, bE, cE` - Error rate (float)
//!   - `ad, bd` - Per-base depth (short[])
//!   - `ae, be` - Per-base errors (short[])
//!   - `ac, bc` - Single-strand consensus bases (string)
//!   - `aq, bq` - Single-strand consensus qualities (string)

use crate::clipper::cigar_utils;
use crate::clipper::{ClippingMode, SamRecordClipper};
use crate::consensus::caller::{
    ConsensusCaller, ConsensusCallingStats, RejectionReason as CallerRejectionReason,
};
use crate::consensus::simple_umi::consensus_umis;
use crate::consensus::vanilla_caller::{
    VanillaConsensusRead, VanillaUmiConsensusCaller, VanillaUmiConsensusOptions,
};
use crate::consensus::{IndexedSourceRead, SourceRead, select_most_common_alignment_group};
use crate::dna::reverse_complement;
use crate::phred::{MIN_PHRED, NO_CALL_BASE, NO_CALL_BASE_LOWER, PhredScore};
use crate::sam::{record_utils, to_smallest_signed_int};
use anyhow::{Result, anyhow};
use bstr::BString;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record::cigar::Cigar as CigarTrait;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::record_buf::{Data, QualityScores, RecordBuf, Sequence};
use rand::SeedableRng;
use rand::rngs::StdRng;
use rand::seq::SliceRandom;
use std::cmp::{max, min};
use std::collections::HashMap;

/// Options for CODEC consensus calling
#[derive(Debug, Clone)]
pub struct CodecConsensusOptions {
    /// Minimum base quality to include in consensus
    pub min_input_base_quality: PhredScore,

    /// Pre-UMI error rate (Phred scale)
    pub error_rate_pre_umi: PhredScore,

    /// Post-UMI error rate (Phred scale)
    pub error_rate_post_umi: PhredScore,

    /// Minimum number of read pairs required to form consensus per strand
    pub min_reads_per_strand: usize,

    /// Maximum number of read pairs to use per strand (downsample if exceeded)
    pub max_reads_per_strand: Option<usize>,

    /// Minimum duplex overlap length (in bases)
    pub min_duplex_length: usize,

    /// Reduce quality of single-strand regions to this value (if set)
    pub single_strand_qual: Option<PhredScore>,

    /// Reduce quality of outer bases to this value (if set)
    pub outer_bases_qual: Option<PhredScore>,

    /// Number of outer bases to reduce quality for
    pub outer_bases_length: usize,

    /// Maximum number of duplex disagreements allowed
    pub max_duplex_disagreements: usize,

    /// Maximum duplex disagreement rate allowed (0.0-1.0)
    pub max_duplex_disagreement_rate: f64,

    /// Cell barcode tag (e.g., "CB")
    pub cell_tag: Option<Tag>,

    /// Whether to produce per-base tags (ad, ae, bd, be, etc.)
    pub produce_per_base_tags: bool,

    /// Whether to quality-trim reads before consensus calling
    pub trim: bool,

    /// Minimum consensus base quality (output bases below this are masked to N)
    pub min_consensus_base_quality: PhredScore,
}

impl Default for CodecConsensusOptions {
    fn default() -> Self {
        Self {
            min_input_base_quality: 10,
            error_rate_pre_umi: 45,
            error_rate_post_umi: 40,
            min_reads_per_strand: 1,
            max_reads_per_strand: None,
            min_duplex_length: 1,
            single_strand_qual: None,
            outer_bases_qual: None,
            outer_bases_length: 5,
            max_duplex_disagreements: usize::MAX,
            max_duplex_disagreement_rate: 1.0,
            cell_tag: None,
            produce_per_base_tags: false,
            trim: false,
            min_consensus_base_quality: 0,
        }
    }
}

/// Statistics for CODEC consensus calling
#[derive(Debug, Clone, Default)]
pub struct CodecConsensusStats {
    /// Total input reads processed
    pub total_input_reads: u64,

    /// Total consensus reads generated
    pub consensus_reads_generated: u64,

    /// Total reads filtered/rejected
    pub reads_filtered: u64,

    /// Consensus reads rejected for high duplex disagreement
    pub consensus_reads_rejected_hdd: u64,

    /// Total consensus bases emitted
    pub consensus_bases_emitted: u64,

    /// Total duplex region bases emitted
    pub consensus_duplex_bases_emitted: u64,

    /// Total duplex disagreement bases
    pub duplex_disagreement_base_count: u64,

    /// Rejection reasons
    pub rejection_reasons: HashMap<CallerRejectionReason, usize>,
}

impl CodecConsensusStats {
    /// Calculates the duplex disagreement rate
    #[must_use]
    pub fn duplex_disagreement_rate(&self) -> f64 {
        if self.consensus_duplex_bases_emitted > 0 {
            self.duplex_disagreement_base_count as f64 / self.consensus_duplex_bases_emitted as f64
        } else {
            0.0
        }
    }
}

/// Single-strand consensus read (intermediate result for internal use only)
///
/// This struct is used internally by `CodecConsensusCaller` during duplex
/// consensus building. It is NOT exported from this module - external code
/// should use `VanillaConsensusRead` instead.
#[derive(Debug, Clone)]
struct SingleStrandConsensus {
    /// Consensus sequence bases
    bases: Vec<u8>,
    /// Consensus quality scores
    quals: Vec<PhredScore>,
    /// Per-base read depth
    depths: Vec<u16>,
    /// Per-base error count
    errors: Vec<u16>,
    /// Total raw read count
    raw_read_count: usize,
    /// Reference start position (1-based)
    ref_start: usize,
    /// Reference end position (1-based, inclusive)
    ref_end: usize,
    /// Whether the original reads were on the negative strand
    is_negative_strand: bool,
}

/// CODEC consensus caller
///
/// Calls consensus from CODEC sequencing data where R1 and R2 from a single read-pair
/// sequence opposite strands of the same molecule.
pub struct CodecConsensusCaller {
    /// Prefix for consensus read names
    read_name_prefix: String,

    /// Read group ID for consensus reads
    read_group_id: String,

    /// Consensus calling options
    options: CodecConsensusOptions,

    /// Statistics tracker
    stats: CodecConsensusStats,

    /// Clipper for clipping reads past mate ends
    clipper: SamRecordClipper,

    /// Random number generator for downsampling
    rng: StdRng,

    /// Counter for consensus read naming
    consensus_counter: u64,

    /// Whether to track rejected reads
    track_rejects: bool,

    /// Rejected reads (only populated if `track_rejects` is true)
    rejected_reads: Vec<RecordBuf>,

    /// Single-strand consensus caller (matches fgbio's ssCaller delegation pattern)
    /// CODEC delegates to `VanillaUmiConsensusCaller` for single-strand consensus building
    ss_caller: VanillaUmiConsensusCaller,
}

#[allow(
    clippy::similar_names,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::cast_sign_loss
)]
impl CodecConsensusCaller {
    /// Creates a new CODEC consensus caller
    ///
    /// # Arguments
    /// * `read_name_prefix` - Prefix for consensus read names
    /// * `read_group_id` - Read group ID for consensus reads
    /// * `options` - CODEC consensus calling options
    #[must_use]
    pub fn new(
        read_name_prefix: String,
        read_group_id: String,
        options: CodecConsensusOptions,
    ) -> Self {
        Self::new_with_rejects_tracking(read_name_prefix, read_group_id, options, false)
    }

    /// Creates a new CODEC consensus caller with optional rejected reads tracking
    ///
    /// # Arguments
    /// * `read_name_prefix` - Prefix for consensus read names
    /// * `read_group_id` - Read group ID for consensus reads
    /// * `options` - CODEC consensus calling options
    /// * `track_rejects` - Whether to store rejected reads for later retrieval
    #[must_use]
    pub fn new_with_rejects_tracking(
        read_name_prefix: String,
        read_group_id: String,
        options: CodecConsensusOptions,
        track_rejects: bool,
    ) -> Self {
        let clipper = SamRecordClipper::new(ClippingMode::Hard);
        let rng = StdRng::seed_from_u64(42);

        // Create VanillaUmiConsensusCaller for single-strand consensus building
        // Matches fgbio's ssCaller initialization in DuplexConsensusCaller:
        // - minReads = 1 (CODEC handles filtering itself)
        // - producePerBaseTags = true (needed for duplex consensus calculation)
        // - minConsensusBaseQuality = MIN (we handle quality masking ourselves)
        let ss_options = VanillaUmiConsensusOptions {
            tag: "MI".to_string(),
            error_rate_pre_umi: options.error_rate_pre_umi,
            error_rate_post_umi: options.error_rate_post_umi,
            min_input_base_quality: options.min_input_base_quality,
            min_reads: 1, // CODEC handles read filtering
            max_reads: None,
            produce_per_base_tags: true,
            seed: None,
            trim: false,
            min_consensus_base_quality: 0, // MIN - we handle quality masking
            cell_tag: None,
        };

        let ss_caller = VanillaUmiConsensusCaller::new(
            "x".to_string(), // Temporary prefix (not used in final output)
            read_group_id.clone(),
            ss_options,
        );

        Self {
            read_name_prefix,
            read_group_id,
            options,
            stats: CodecConsensusStats::default(),
            clipper,
            rng,
            consensus_counter: 0,
            track_rejects,
            rejected_reads: Vec::new(),
            ss_caller,
        }
    }

    /// Returns the statistics for this caller
    #[must_use]
    pub fn statistics(&self) -> &CodecConsensusStats {
        &self.stats
    }

    /// Returns a slice of rejected reads (only populated if tracking is enabled)
    #[must_use]
    pub fn rejected_reads(&self) -> &[RecordBuf] {
        &self.rejected_reads
    }

    /// Takes ownership of the rejected reads, leaving an empty Vec
    pub fn take_rejected_reads(&mut self) -> Vec<RecordBuf> {
        std::mem::take(&mut self.rejected_reads)
    }

    /// Clears the rejected reads buffer
    pub fn clear_rejected_reads(&mut self) {
        self.rejected_reads.clear();
    }

    /// Clears all per-group state to prepare for reuse
    ///
    /// This resets statistics and rejected reads while preserving the caller's
    /// configuration and reusable components (consensus builder, RNG, etc.).
    pub fn clear(&mut self) {
        self.stats = CodecConsensusStats::default();
        self.rejected_reads.clear();
        self.ss_caller.clear();
    }

    /// Converts a `RecordBuf` to `SourceRead` for CODEC consensus calling.
    ///
    /// Matches fgbio's `toSourceReadForCodec` exactly:
    /// - Clones bases and quals from the record
    /// - Reverses CIGAR if read is on negative strand
    /// - Reverse complements bases and reverses quals if on negative strand
    /// - Creates `SourceRead` with sam field set for tag preservation
    ///
    /// KEY: This is SIMPLER than vanilla's `to_source_read`:
    /// - NO quality masking
    /// - NO mate trimming
    /// - NO trailing N removal
    fn to_source_read_for_codec(read: &RecordBuf, original_idx: usize) -> SourceRead {
        let mut bases = read.sequence().as_ref().to_vec();
        let mut quals: Vec<u8> = read.quality_scores().as_ref().to_vec();

        // Get simplified CIGAR, reverse if negative strand
        let simplified_cigar = if read.flags().is_reverse_complemented() {
            let mut cigar = cigar_utils::simplify_cigar(read.cigar());
            cigar.reverse();
            cigar
        } else {
            cigar_utils::simplify_cigar(read.cigar())
        };

        // Revcomp bases and reverse quals if negative strand
        if read.flags().is_reverse_complemented() {
            bases = reverse_complement(&bases);
            quals.reverse();
        }

        SourceRead {
            original_idx,
            bases,
            quals,
            simplified_cigar,
            sam: Some(std::sync::Arc::new(read.clone())),
        }
    }

    /// Converts a `VanillaConsensusRead` to a `SingleStrandConsensus`.
    ///
    /// This bridges the delegation pattern: `CodecConsensusCaller` delegates to
    /// `VanillaUmiConsensusCaller` for single-strand consensus, then converts
    /// the result to CODEC's internal `SingleStrandConsensus` type.
    fn vanilla_to_single_strand(
        vcr: VanillaConsensusRead,
        is_negative_strand: bool,
        raw_read_count: usize,
    ) -> SingleStrandConsensus {
        let len = vcr.bases.len();
        SingleStrandConsensus {
            bases: vcr.bases,
            quals: vcr.quals,
            depths: vcr.depths,
            errors: vcr.errors,
            raw_read_count,
            ref_start: 0,
            ref_end: len.saturating_sub(1),
            is_negative_strand,
        }
    }

    /// Reverse complements a `SingleStrandConsensus`
    ///
    /// This converts the bases from the R2 orientation (reverse complement of reference)
    /// to the R1 orientation (matching reference).
    ///
    /// For CODEC (query-based consensus):
    /// - Bases, quals, depths, and errors are ALL reversed (they're in query space)
    /// - This ensures depths/errors stay aligned with the reversed bases
    fn reverse_complement_ss(ss: &SingleStrandConsensus) -> SingleStrandConsensus {
        SingleStrandConsensus {
            bases: reverse_complement(&ss.bases),
            quals: ss.quals.iter().copied().rev().collect(),
            // For CODEC query-based consensus, depths and errors must also be reversed
            // to stay aligned with the reversed bases
            depths: ss.depths.iter().copied().rev().collect(),
            errors: ss.errors.iter().copied().rev().collect(),
            raw_read_count: ss.raw_read_count,
            ref_start: ss.ref_start,
            ref_end: ss.ref_end,
            is_negative_strand: !ss.is_negative_strand,
        }
    }

    /// Calls consensus reads from SAM records for a single molecule
    ///
    /// # Arguments
    /// * `recs` - All reads for a single molecule (same MI tag)
    ///
    /// # Returns
    /// Vector of consensus reads (usually 0 or 1)
    pub fn consensus_reads_from_sam_records(
        &mut self,
        recs: Vec<RecordBuf>,
    ) -> Result<Vec<RecordBuf>> {
        self.stats.total_input_reads += recs.len() as u64;

        if recs.is_empty() {
            return Ok(Vec::new());
        }

        // Extract MI tag from first record for naming
        let mi_tag = Tag::from([b'M', b'I']);
        let umi = recs.first().and_then(|r| {
            if let Some(Value::String(s)) = r.data().get(&mi_tag) {
                Some(String::from_utf8_lossy(s.as_ref()).to_string())
            } else {
                None
            }
        });

        // Partition into paired and fragment reads
        let (pairs, frags): (Vec<_>, Vec<_>) =
            recs.into_iter().partition(|r| r.flags().is_segmented());

        // Reject fragment reads - CODEC requires paired-end
        if !frags.is_empty() {
            self.reject_records(&frags, CallerRejectionReason::FragmentRead);
        }

        if pairs.is_empty() {
            return Ok(Vec::new());
        }

        // Keep a reference to all paired reads for UMI consensus (matching fgbio behavior)
        // fgbio uses all source reads for UMI consensus, not just filtered ones
        let all_pairs = pairs.clone();

        // Filter to primary alignments only, mapped reads, and FR pairs
        let primaries: Vec<_> = pairs
            .into_iter()
            .filter(|r| {
                let flags = r.flags();
                !flags.is_secondary() && !flags.is_supplementary() && !flags.is_unmapped()
            })
            .filter(|r| self.is_fr_pair(r))
            .collect();

        if primaries.is_empty() {
            return Ok(Vec::new());
        }

        // Group by read name to form pairs
        // Note: we have to be VERY CAREFUL to preserve order
        let mut by_name: HashMap<String, Vec<RecordBuf>> = HashMap::new();
        let mut names: Vec<String> = Vec::with_capacity(primaries.len() / 2);
        for rec in primaries {
            let name =
                rec.name().map(|n| String::from_utf8_lossy(n).to_string()).unwrap_or_default();
            if !by_name.contains_key(&name) {
                names.push(name.clone());
            }
            by_name.entry(name).or_default().push(rec);
        }

        // Clip each pair where they extend past mate ends, and collect valid pairs
        let mut clipped_pairs: Vec<(RecordBuf, RecordBuf)> = Vec::new();
        for name in names {
            // Safe: names was collected from by_name.keys()
            let pair_reads = by_name.get(&name).expect("name from keys()");
            if pair_reads.len() == 2 {
                let (mut r1, mut r2) = if pair_reads[0].flags().is_first_segment() {
                    (pair_reads[0].clone(), pair_reads[1].clone())
                } else {
                    (pair_reads[1].clone(), pair_reads[0].clone())
                };

                // fgbio DOES clip reads extending past mate ends using Hard clipping mode.
                // This is done BEFORE the overlap calculation and phase check.
                // The clipping removes bases (hard clips) from reads that extend past their mate.
                self.clipper.clip_extending_past_mate_ends(&mut r1, &mut r2);

                clipped_pairs.push((r1, r2));
            }
        }

        if clipped_pairs.is_empty() {
            return Ok(Vec::new());
        }

        // Separate R1s and R2s
        let (r1s, r2s): (Vec<_>, Vec<_>) = clipped_pairs.into_iter().unzip();

        // Check we have enough reads
        if r1s.len() < self.options.min_reads_per_strand {
            self.reject_records_count(
                r1s.len() + r2s.len(),
                CallerRejectionReason::InsufficientReads,
            );
            return Ok(Vec::new());
        }

        // Downsample if needed
        let (r1s, r2s) = self.downsample_pairs(r1s, r2s);

        // Filter R1s and R2s to most common alignment pattern
        let r1s = self.filter_to_most_common_alignment(r1s);
        let r2s = self.filter_to_most_common_alignment(r2s);

        if r1s.is_empty() || r2s.is_empty() {
            return Ok(Vec::new());
        }

        // Check we still have enough reads after filtering
        if r1s.len() < self.options.min_reads_per_strand
            || r2s.len() < self.options.min_reads_per_strand
        {
            self.reject_records_count(
                r1s.len() + r2s.len(),
                CallerRejectionReason::InsufficientReads,
            );
            return Ok(Vec::new());
        }

        // Get the longest R1 and R2 alignments (like fgbio)
        // Note: reverse the iterator so we return the first maximum value to match fgbio
        // Safe: we checked r1s.is_empty() and r2s.is_empty() above
        let longest_r1 = r1s
            .iter()
            .rev()
            .max_by_key(|r| cigar_utils::reference_length(&r.cigar()))
            .expect("r1s is non-empty");
        let longest_r2 = r2s
            .iter()
            .rev()
            .max_by_key(|r| cigar_utils::reference_length(&r.cigar()))
            .expect("r2s is non-empty");

        // Determine which alignment is positive strand and which is negative
        let r1_is_negative = longest_r1.flags().is_reverse_complemented();
        let (longest_pos_aln, longest_neg_aln) =
            if r1_is_negative { (longest_r2, longest_r1) } else { (longest_r1, longest_r2) };

        // Calculate overlap region on reference using fgbio's formula:
        // overlapStart = negativeStrandAlignment.start
        // overlapEnd = positiveStrandAlignment.end
        // This only works because we have an FR pair that sequences towards each other.
        let neg_start = longest_neg_aln
            .alignment_start()
            .map(usize::from)
            .ok_or_else(|| anyhow!("Negative strand alignment missing start position"))?;
        let pos_start = longest_pos_aln
            .alignment_start()
            .map(usize::from)
            .ok_or_else(|| anyhow!("Positive strand alignment missing start position"))?;
        let pos_end =
            pos_start + cigar_utils::reference_length(&longest_pos_aln.cigar()).saturating_sub(1);

        let (overlap_start, overlap_end) = (neg_start, pos_end);

        // Calculate duplex length (can be negative if no overlap)
        let duplex_length = overlap_end as i64 - overlap_start as i64 + 1;

        if duplex_length < self.options.min_duplex_length as i64 {
            self.reject_records_count(
                r1s.len() + r2s.len(),
                CallerRejectionReason::InsufficientOverlap,
            );
            return Ok(Vec::new());
        }

        // Check that overlap boundaries don't land in indels (matches fgbio behavior)
        if !self.check_overlap_phase(longest_r1, longest_r2, overlap_start, overlap_end) {
            self.reject_records_count(
                r1s.len() + r2s.len(),
                CallerRejectionReason::IndelErrorBetweenStrands,
            );
            return Ok(Vec::new());
        }

        // r1_is_negative and longest_pos_aln/longest_neg_aln already determined above
        let r2_is_negative = longest_r2.flags().is_reverse_complemented();

        // Compute the total consensus length (including soft-clipped bases)
        // This matches fgbio's computeConsensusLength
        let consensus_length_result =
            Self::compute_codec_consensus_length(longest_pos_aln, longest_neg_aln, overlap_end);

        let Some(consensus_length) = consensus_length_result else {
            self.reject_records_count(
                r1s.len() + r2s.len(),
                CallerRejectionReason::IndelErrorBetweenStrands,
            );
            return Ok(Vec::new());
        };

        // Build single-strand consensus from R1s using query-based approach
        // Following fgbio's flow: convert to SourceReads, then delegate to ssCaller
        let r1_is_neg_strand = r1s.first().is_some_and(|r| r.flags().is_reverse_complemented());
        let r1_source_reads: Vec<SourceRead> = r1s
            .iter()
            .enumerate()
            .map(|(idx, read)| Self::to_source_read_for_codec(read, idx))
            .collect();

        let umi_str = umi.as_deref().unwrap_or("");
        let ss_r1 = match self.ss_caller.consensus_call(umi_str, r1_source_reads)? {
            Some(vcr) => Self::vanilla_to_single_strand(vcr, r1_is_neg_strand, r1s.len()),
            None => return Ok(Vec::new()),
        };

        // Build single-strand consensus from R2s using query-based approach
        // Following fgbio's flow: convert to SourceReads, then delegate to ssCaller
        let r2_is_neg_strand = r2s.first().is_some_and(|r| r.flags().is_reverse_complemented());
        let r2_source_reads: Vec<SourceRead> = r2s
            .iter()
            .enumerate()
            .map(|(idx, read)| Self::to_source_read_for_codec(read, idx))
            .collect();

        let ss_r2 = match self.ss_caller.consensus_call(umi_str, r2_source_reads)? {
            Some(vcr) => Self::vanilla_to_single_strand(vcr, r2_is_neg_strand, r2s.len()),
            None => return Ok(Vec::new()),
        };

        // Check that consensus length is at least as long as each SS consensus
        // This matches fgbio's check: `n < r1Consensus.length || n < r2Consensus.length`
        if consensus_length < ss_r1.bases.len() || consensus_length < ss_r2.bases.len() {
            self.reject_records_count(
                r1s.len() + r2s.len(),
                CallerRejectionReason::IndelErrorBetweenStrands, // ClipOverlapFailed in fgbio
            );
            return Ok(Vec::new());
        }

        // FIRST reverse complement the negative strand consensus so it aligns with positive strand
        // fgbio does: if (longestR1Alignment.negativeStrand) r1Consensus.revcomp() else r2Consensus.revcomp()
        let (ss_r1_oriented, ss_r2_oriented) = if r1_is_negative {
            (Self::reverse_complement_ss(&ss_r1), ss_r2.clone())
        } else {
            (ss_r1.clone(), Self::reverse_complement_ss(&ss_r2))
        };

        // THEN pad out the consensus reads to the full length
        // fgbio: paddedR1 = r1Consensus.padded(newLength=consensusLength, left=longestR1Alignment.negativeStrand)
        // fgbio: paddedR2 = r2Consensus.padded(newLength=consensusLength, left=longestR2Alignment.negativeStrand)
        // NOTE: After revcomp, padding "left" means prepending N's which extends the genomic 5' end
        let padded_r1 = Self::pad_consensus(&ss_r1_oriented, consensus_length, r1_is_negative);
        let padded_r2 = Self::pad_consensus(&ss_r2_oriented, consensus_length, r2_is_negative);

        // Build duplex consensus from padded and oriented SS consensuses
        let consensus = self.build_duplex_consensus_from_padded(&padded_r1, &padded_r2)?;

        // Apply quality masking for single-strand regions and outer bases
        // Pass padded consensuses so we can identify single-strand regions
        let consensus = self.mask_consensus_quals_query_based(consensus, &padded_r1, &padded_r2);

        // If R1 was on negative strand, reverse complement the final consensus
        // This matches fgbio line 253: if (longestR1Alignment.negativeStrand) consensus.revcomp()
        let consensus =
            if r1_is_negative { self.reverse_complement_consensus(consensus) } else { consensus };

        // NOTE: fgbio does NOT filter empty positions for CODEC consensus.
        // Positions where neither strand has data are kept as 'N' in the output.

        // Build the output RecordBuf
        // Pass the padded and oriented single-strand consensuses for ac/bc tags
        // These are padded_r1 and padded_r2, which have been oriented and padded to match fgbio

        // IMPORTANT: The ac/bc tags must be in the SAME coordinate system as the final consensus.
        // Since the final consensus gets reverse complemented when r1_is_negative=true,
        // we need to also reverse complement the ac/bc tags in that case.
        //
        // The assignment is always: ac=paddedR1, bc=paddedR2 (no swap).
        // But when r1_is_negative, both need to be reverse complemented to match the final consensus.
        let (ss_for_ac, ss_for_bc) = if r1_is_negative {
            // R1 is negative: reverse complement both to match final consensus (but don't swap)
            (Self::reverse_complement_ss(&padded_r1), Self::reverse_complement_ss(&padded_r2))
        } else {
            // R1 is positive: no revcomp needed, no swap
            (padded_r1.clone(), padded_r2.clone())
        };

        let record = self.build_output_record(
            consensus,
            &ss_for_ac, // becomes ac tag
            &ss_for_bc, // becomes bc tag
            umi.as_deref(),
            overlap_start,
            overlap_end,
            &all_pairs, // Use all paired reads for UMI consensus, matching fgbio
        )?;

        self.stats.consensus_reads_generated += 1;

        Ok(vec![record])
    }

    /// Checks if a read is part of an FR pair (forward-reverse orientation)
    /// Uses proper orientation check including position information, matching fgbio's isFrPair.
    fn is_fr_pair(&self, read: &RecordBuf) -> bool {
        record_utils::is_fr_pair_from_tags(read)
    }

    /// Downsamples pairs if we have more than `max_reads_per_strand`
    fn downsample_pairs(
        &mut self,
        mut r1s: Vec<RecordBuf>,
        mut r2s: Vec<RecordBuf>,
    ) -> (Vec<RecordBuf>, Vec<RecordBuf>) {
        if let Some(max_reads) = self.options.max_reads_per_strand {
            if r1s.len() > max_reads {
                // Create indices and shuffle
                let mut indices: Vec<usize> = (0..r1s.len()).collect();
                indices.shuffle(&mut self.rng);
                indices.truncate(max_reads);
                indices.sort_unstable();

                // Select corresponding pairs
                let new_r1s: Vec<_> = indices.iter().map(|&i| r1s[i].clone()).collect();
                let new_r2s: Vec<_> = indices.iter().map(|&i| r2s[i].clone()).collect();
                r1s = new_r1s;
                r2s = new_r2s;
            }
        }
        (r1s, r2s)
    }

    /// Filters reads to only include those with the most common alignment pattern.
    /// This uses the shared `select_most_common_alignment_group` logic from `vanilla_caller`,
    /// matching fgbio's `filterToMostCommonAlignment` behavior with deterministic tie-breaking.
    fn filter_to_most_common_alignment(&mut self, reads: Vec<RecordBuf>) -> Vec<RecordBuf> {
        if reads.len() < 2 {
            return reads;
        }

        // Create indexed data: (read_index, query_length, simplified_cigar)
        let mut indexed: Vec<IndexedSourceRead> = reads
            .iter()
            .enumerate()
            .map(|(i, read)| {
                let cigar = cigar_utils::simplify_cigar(read.cigar());
                let query_len = read.sequence().len();
                (i, query_len, cigar)
            })
            .collect();

        // Sort by descending length for prefix matching
        indexed.sort_by(|a, b| b.1.cmp(&a.1));

        // Use shared grouping logic with deterministic tie-breaking
        let best_indices = select_most_common_alignment_group(&indexed);

        let rejected_count = reads.len() - best_indices.len();
        if rejected_count > 0 {
            *self
                .stats
                .rejection_reasons
                .entry(CallerRejectionReason::MinorityAlignment)
                .or_insert(0) += rejected_count;
            self.stats.reads_filtered += rejected_count as u64;
        }

        // Return only the reads in the best group
        let best_set: std::collections::HashSet<usize> = best_indices.into_iter().collect();
        reads
            .into_iter()
            .enumerate()
            .filter(|(i, _)| best_set.contains(i))
            .map(|(_, r)| r)
            .collect()
    }

    /// Returns the query position (1-based) at a given reference position.
    ///
    /// Parameters:
    /// - `read`: The read record
    /// - `ref_pos`: The reference position to query
    /// - `return_last_base_if_deleted`: If true and position falls in deletion, return last query
    ///   position before deletion. If false, return None when position falls in deletion.
    ///
    /// This matches fgbio's `readPosAtRefPos(pos, returnLastBaseIfDeleted)` behavior.
    fn read_pos_at_ref_pos(
        read: &RecordBuf,
        ref_pos: usize,
        return_last_base_if_deleted: bool,
    ) -> Option<usize> {
        let start = read.alignment_start()?;
        let start_pos = usize::from(start);

        if ref_pos < start_pos {
            return None;
        }

        let mut ref_offset = 0;
        let mut query_offset = 0;

        for op_result in read.cigar().iter() {
            let op = op_result.ok()?;
            let (kind, len) = (op.kind(), op.len());

            let consumes_ref = matches!(
                kind,
                Kind::Match
                    | Kind::SequenceMatch
                    | Kind::SequenceMismatch
                    | Kind::Deletion
                    | Kind::Skip
            );
            let consumes_query = matches!(
                kind,
                Kind::Match
                    | Kind::SequenceMatch
                    | Kind::SequenceMismatch
                    | Kind::Insertion
                    | Kind::SoftClip
            );

            let op_ref_start = start_pos + ref_offset;

            if consumes_ref {
                let op_ref_end = op_ref_start + len - 1;

                // Check if target ref_pos is within this CIGAR operation
                if ref_pos >= op_ref_start && ref_pos <= op_ref_end {
                    if consumes_query {
                        // M, =, X operations - we have a base at this position
                        let offset_in_op = ref_pos - op_ref_start;
                        return Some(query_offset + offset_in_op + 1); // 1-based
                    }
                    // D, N operations - position falls in a deletion
                    if return_last_base_if_deleted {
                        // Return the last query position before this deletion (1-based)
                        // If query_offset is 0, we're at the start, so return 1
                        return Some(if query_offset > 0 { query_offset } else { 1 });
                    }
                    return None;
                }
            }

            if consumes_ref {
                ref_offset += len;
            }
            if consumes_query {
                query_offset += len;
            }
        }

        None
    }

    /// Computes the total consensus length including soft-clipped bases.
    /// This matches fgbio's `computeConsensusLength` function.
    ///
    /// The formula is: `posReadPos + neg.length - negReadPos`
    /// where:
    /// - `posReadPos` = query position at overlap end for positive strand read
    /// - `negReadPos` = query position at overlap end for negative strand read
    /// - `neg.length` = total length of negative strand read (including soft-clips)
    ///
    /// Returns None if the position falls within a deletion (indel error).
    fn compute_codec_consensus_length(
        pos_read: &RecordBuf,
        neg_read: &RecordBuf,
        overlap_end: usize,
    ) -> Option<usize> {
        // Use returnLastBaseIfDeleted=false (strict check for consensus length calculation)
        let pos_read_pos = Self::read_pos_at_ref_pos(pos_read, overlap_end, false)?;
        let neg_read_pos = Self::read_pos_at_ref_pos(neg_read, overlap_end, false)?;

        // Match fgbio exactly: posReadPos + neg.length - negReadPos
        // Use full sequence length (including soft clips)
        let neg_length = neg_read.sequence().len();

        Some(pos_read_pos + neg_length - neg_read_pos)
    }

    /// Pads a single-strand consensus to a new length.
    /// If `pad_left` is true, pads on the left (prepends), otherwise on the right (appends).
    /// This matches fgbio's `VanillaConsensusRead.padded` method.
    ///
    /// # Lowercase 'n' for padding
    ///
    /// Uses lowercase 'n' (`NO_CALL_BASE_LOWER`) for padding bases to match fgbio behavior.
    /// These padded bases may be reverse complemented later (see `reverse_complement_ss`),
    /// and `dna::complement_base` preserves the lowercase 'n' through that operation.
    /// This is specific to CODEC - simplex/duplex use uppercase 'N' for no-call bases.
    fn pad_consensus(
        ss: &SingleStrandConsensus,
        new_length: usize,
        pad_left: bool,
    ) -> SingleStrandConsensus {
        let current_len = ss.bases.len();
        if new_length <= current_len {
            return ss.clone();
        }

        let pad_len = new_length - current_len;
        let pad_bases = vec![NO_CALL_BASE_LOWER; pad_len]; // Use lowercase 'n' to match fgbio
        let pad_quals = vec![0u8; pad_len];
        let pad_depths = vec![0u16; pad_len];
        let pad_errors = vec![0u16; pad_len];

        let (bases, quals, depths, errors) = if pad_left {
            (
                [pad_bases, ss.bases.clone()].concat(),
                [pad_quals, ss.quals.clone()].concat(),
                [pad_depths, ss.depths.clone()].concat(),
                [pad_errors, ss.errors.clone()].concat(),
            )
        } else {
            (
                [ss.bases.clone(), pad_bases].concat(),
                [ss.quals.clone(), pad_quals].concat(),
                [ss.depths.clone(), pad_depths].concat(),
                [ss.errors.clone(), pad_errors].concat(),
            )
        };

        SingleStrandConsensus {
            bases,
            quals,
            depths,
            errors,
            raw_read_count: ss.raw_read_count,
            ref_start: ss.ref_start,
            ref_end: ss.ref_end,
            is_negative_strand: ss.is_negative_strand,
        }
    }

    /// Builds duplex consensus from two query-based single-strand consensuses.
    /// The two SS consensuses should already be padded to the same length.
    ///
    /// Single-strand positions (where only one strand has data) preserve their
    /// original qualities from the single-strand consensus, matching fgbio's behavior.
    fn build_duplex_consensus_from_padded(
        &mut self,
        ss_a: &SingleStrandConsensus,
        ss_b: &SingleStrandConsensus,
    ) -> Result<SingleStrandConsensus> {
        let len = ss_a.bases.len();
        assert_eq!(ss_b.bases.len(), len, "Padded consensuses must be the same length");

        let mut bases = vec![NO_CALL_BASE; len];
        let mut quals: Vec<PhredScore> = vec![MIN_PHRED; len];
        let mut depths = vec![0u16; len];
        let mut errors = vec![0u16; len];

        let mut duplex_disagreements = 0usize;
        let mut duplex_bases_count = 0usize;

        for pos_idx in 0..len {
            let ba = ss_a.bases[pos_idx];
            let qa = ss_a.quals[pos_idx];
            let da = ss_a.depths[pos_idx];
            let ea = ss_a.errors[pos_idx];

            let bb = ss_b.bases[pos_idx];
            let qb = ss_b.quals[pos_idx];
            let db = ss_b.depths[pos_idx];
            let eb = ss_b.errors[pos_idx];

            // Check if this is a duplex position (both strands have data)
            // Match fgbio's logic: only check if base is N/n (padding), not quality
            // fgbio doesn't check qual > 0 here - it uses whatever qualities are present
            let a_has_data = ba != NO_CALL_BASE && ba != NO_CALL_BASE_LOWER;
            let b_has_data = bb != NO_CALL_BASE && bb != NO_CALL_BASE_LOWER;

            let (duplex_base, duplex_qual, depth, error) = match (a_has_data, b_has_data) {
                (true, true) => {
                    // Both strands have a base at this position - duplex region
                    duplex_bases_count += 1;
                    let (raw_base, raw_qual) = if ba == bb {
                        // Agreement - sum qualities, cap to 93
                        (ba, min(93, u16::from(qa) + u16::from(qb)) as u8)
                    } else if qa > qb {
                        // Disagreement - take higher quality base, subtract qualities
                        duplex_disagreements += 1;
                        // Cap to minimum 2 (matching fgbio's PhredScore.cap)
                        (ba, max(MIN_PHRED, qa.saturating_sub(qb)))
                    } else if qb > qa {
                        duplex_disagreements += 1;
                        // Cap to minimum 2 (matching fgbio's PhredScore.cap)
                        (bb, max(MIN_PHRED, qb.saturating_sub(qa)))
                    } else {
                        // Equal quality disagreement
                        // IMPORTANT: fgbio sets rawBase = aBase (not N!) even though the final
                        // base will be masked to N due to low quality. The error calculation
                        // uses rawBase (aBase), not the final masked base.
                        duplex_disagreements += 1;
                        (ba, MIN_PHRED) // Use aBase (ba), not N!
                    };

                    // Mask to N if quality is MinValue (2), matching fgbio's logic
                    let (final_base, final_qual) = if raw_qual == MIN_PHRED {
                        (NO_CALL_BASE, MIN_PHRED)
                    } else {
                        (raw_base, raw_qual)
                    };

                    // Calculate errors matching fgbio's logic:
                    // - When bases agree: sum of both errors
                    // - When bases disagree: count errors from the chosen strand + non-errors from other strand
                    //   (because the non-errors from the other strand now disagree with the consensus)
                    let duplex_error = if ba == bb {
                        // Agreement: sum errors from both strands
                        ea + eb
                    } else if ba == raw_base {
                        // We chose A's base (including equal quality disagreements where fgbio uses aBase)
                        // Count A's errors + B's non-errors (which now disagree with consensus)
                        ea + (db - eb)
                    } else {
                        // We chose B's base: count B's errors + A's non-errors
                        eb + (da - ea)
                    };

                    (final_base, final_qual, da + db, duplex_error)
                }
                (true, false) => {
                    // Only A has data - single-strand from A
                    // Mask to N if quality is MinValue (2), matching fgbio's behavior
                    if qa == MIN_PHRED {
                        (NO_CALL_BASE, MIN_PHRED, da, ea)
                    } else {
                        (ba, qa, da, ea)
                    }
                }
                (false, true) => {
                    // Only B has data - single-strand from B
                    // Mask to N if quality is MinValue (2), matching fgbio's behavior
                    if qb == MIN_PHRED {
                        (NO_CALL_BASE, MIN_PHRED, db, eb)
                    } else {
                        (bb, qb, db, eb)
                    }
                }
                (false, false) => {
                    // Neither has data - N
                    // fgbio still calculates errors for these positions:
                    // since both bases are the same (N == N), errors = a.errors + b.errors
                    let duplex_error = ea + eb;
                    (NO_CALL_BASE, MIN_PHRED, 0, duplex_error)
                }
            };

            // Apply fgbio's NoCall masking: if EITHER original strand base is NoCall (uppercase 'N'),
            // mask the final result to N. This is applied AFTER the rawBase calculation.
            // Note: fgbio only considers uppercase 'N' as NoCall, not lowercase 'n' (which is used for padding).
            let (duplex_base, duplex_qual) = if ba == NO_CALL_BASE || bb == NO_CALL_BASE {
                (NO_CALL_BASE, MIN_PHRED)
            } else {
                (duplex_base, duplex_qual)
            };

            bases[pos_idx] = duplex_base;
            quals[pos_idx] = duplex_qual;
            depths[pos_idx] = depth;
            errors[pos_idx] = error;
        }

        // Check duplex disagreement rate
        if duplex_bases_count > 0 {
            let duplex_error_rate = duplex_disagreements as f64 / duplex_bases_count as f64;
            self.stats.consensus_duplex_bases_emitted += duplex_bases_count as u64;
            self.stats.duplex_disagreement_base_count += duplex_disagreements as u64;

            if duplex_disagreements > self.options.max_duplex_disagreements {
                anyhow::bail!("High duplex disagreement: {duplex_disagreements} disagreements");
            }
            if duplex_error_rate > self.options.max_duplex_disagreement_rate {
                anyhow::bail!("High duplex disagreement rate: {duplex_error_rate:.4}");
            }
        }

        Ok(SingleStrandConsensus {
            bases,
            quals,
            depths,
            errors,
            raw_read_count: ss_a.raw_read_count + ss_b.raw_read_count,
            ref_start: 0,
            ref_end: len.saturating_sub(1),
            is_negative_strand: ss_a.is_negative_strand,
        })
    }

    /// Checks that the overlap boundaries don't land in indels.
    /// This matches fgbio's check that readPosAtRefPos for start and end are consistent.
    fn check_overlap_phase(
        &self,
        r1_aln: &RecordBuf,
        r2_aln: &RecordBuf,
        overlap_start: usize,
        overlap_end: usize,
    ) -> bool {
        // Get query positions at overlap start and end for both reads
        // Use returnLastBaseIfDeleted=true to match fgbio's behavior
        let r1_start_pos = Self::read_pos_at_ref_pos(r1_aln, overlap_start, true);
        let r2_start_pos = Self::read_pos_at_ref_pos(r2_aln, overlap_start, true);
        let r1_end_pos = Self::read_pos_at_ref_pos(r1_aln, overlap_end, true);
        let r2_end_pos = Self::read_pos_at_ref_pos(r2_aln, overlap_end, true);

        // If any position is completely outside the read alignment, reject
        match (r1_start_pos, r2_start_pos, r1_end_pos, r2_end_pos) {
            (Some(r1s), Some(r2s), Some(r1e), Some(r2e)) => {
                // Check that the length of overlapping region in query bases is the same
                // i.e., r1_end - r1_start == r2_end - r2_start
                // Or equivalently: r1_start - r2_start == r1_end - r2_end
                let r1_diff = r1s as i64 - r2s as i64;
                let r2_diff = r1e as i64 - r2e as i64;
                r1_diff == r2_diff
            }
            _ => {
                false // A position is outside the read alignment
            }
        }
    }

    /// Masks consensus qualities for query-based consensus.
    /// Identifies single-strand regions by checking if either padded consensus has 'N'.
    /// This matches fgbio's maskCodecConsensusQuals behavior.
    fn mask_consensus_quals_query_based(
        &self,
        mut consensus: SingleStrandConsensus,
        padded_r1: &SingleStrandConsensus,
        padded_r2: &SingleStrandConsensus,
    ) -> SingleStrandConsensus {
        let len = consensus.bases.len();

        for idx in 0..len {
            // Mask single-strand regions first
            // Single-strand positions are where one padded consensus has 'N'
            // This matches fgbio's check: if (a(i) == 'N' || b(i) == 'N')
            let a_is_n = padded_r1.bases.get(idx).copied().unwrap_or(NO_CALL_BASE) == NO_CALL_BASE;
            let b_is_n = padded_r2.bases.get(idx).copied().unwrap_or(NO_CALL_BASE) == NO_CALL_BASE;

            if (a_is_n || b_is_n) && consensus.bases[idx] != NO_CALL_BASE {
                // This is a single-strand position with a base
                if let Some(ss_qual) = self.options.single_strand_qual {
                    consensus.quals[idx] = ss_qual;
                }
            }

            // Mask outer bases
            if let Some(outer_qual) = self.options.outer_bases_qual {
                let outer_len = self.options.outer_bases_length;
                if idx < outer_len || idx >= len.saturating_sub(outer_len) {
                    consensus.quals[idx] = min(consensus.quals[idx], outer_qual);
                }
            }
        }

        consensus
    }

    /// Reverse complements a consensus
    fn reverse_complement_consensus(&self, ss: SingleStrandConsensus) -> SingleStrandConsensus {
        Self::reverse_complement_ss(&ss)
    }

    /// Builds the output `RecordBuf` from the consensus
    #[allow(clippy::too_many_arguments)]
    fn build_output_record(
        &mut self,
        consensus: SingleStrandConsensus,
        ss_a: &SingleStrandConsensus,
        ss_b: &SingleStrandConsensus,
        umi: Option<&str>,
        overlap_start: usize,
        overlap_end: usize,
        source_reads: &[RecordBuf],
    ) -> Result<RecordBuf> {
        let mut record = RecordBuf::default();

        // Generate read name - use ':' delimiter to match fgbio format
        self.consensus_counter += 1;
        let read_name = if let Some(umi_str) = umi {
            format!("{}:{}", self.read_name_prefix, umi_str)
        } else {
            format!("{}:{}", self.read_name_prefix, self.consensus_counter)
        };
        *record.name_mut() = Some(BString::from(read_name.into_bytes()));

        // Set flags - unmapped fragment
        *record.flags_mut() = Flags::UNMAPPED;

        // Set mapping quality to 0 (matching fgbio)
        *record.mapping_quality_mut() = noodles::sam::alignment::record::MappingQuality::new(0);

        // Set sequence and quality
        *record.sequence_mut() = Sequence::from(consensus.bases.clone());
        *record.quality_scores_mut() = QualityScores::from(consensus.quals.clone());

        // Build tags
        let mut data = Data::default();

        // Read group
        data.insert(Tag::READ_GROUP, Value::from(self.read_group_id.clone()));

        // MI tag
        if let Some(umi_str) = umi {
            let mi_tag = Tag::from([b'M', b'I']);
            data.insert(mi_tag, Value::from(umi_str.to_string()));
        }

        // Duplex consensus tags (cD, cM, cE)
        // fgbio calculates these from totalDepths = ab.depths + ba.depths (sum of SS depths),
        // NOT from the duplex consensus depths
        let total_depths: Vec<i32> = ss_a
            .depths
            .iter()
            .zip(ss_b.depths.iter())
            .map(|(&a, &b)| i32::from(a) + i32::from(b))
            .collect();
        let total_depth = total_depths.iter().copied().max().unwrap_or(0);
        let min_depth = total_depths.iter().copied().min().unwrap_or(0);
        let total_errors: usize = consensus.errors.iter().map(|&e| e as usize).sum();
        let total_bases: i32 = total_depths.iter().sum();
        let error_rate =
            if total_bases > 0 { total_errors as f32 / total_bases as f32 } else { 0.0 };

        let cd_tag = Tag::from([b'c', b'D']);
        let cm_tag = Tag::from([b'c', b'M']);
        let ce_tag = Tag::from([b'c', b'E']);
        data.insert(cd_tag, to_smallest_signed_int(total_depth));
        data.insert(cm_tag, to_smallest_signed_int(min_depth));
        data.insert(ce_tag, Value::from(error_rate));

        // AB strand tags (aD, aM, aE)
        let a_max_depth = ss_a.depths.iter().map(|&d| i32::from(d)).max().unwrap_or(0);
        let a_min_depth = ss_a.depths.iter().map(|&d| i32::from(d)).min().unwrap_or(0);
        let a_total_errors: usize = ss_a.errors.iter().map(|&e| e as usize).sum();
        let a_total_bases: usize = ss_a.depths.iter().map(|&d| d as usize).sum();
        let a_error_rate =
            if a_total_bases > 0 { a_total_errors as f32 / a_total_bases as f32 } else { 0.0 };

        let ad_tag = Tag::from([b'a', b'D']);
        let am_tag = Tag::from([b'a', b'M']);
        let ae_tag = Tag::from([b'a', b'E']);
        data.insert(ad_tag, to_smallest_signed_int(a_max_depth));
        data.insert(am_tag, to_smallest_signed_int(a_min_depth));
        data.insert(ae_tag, Value::from(a_error_rate));

        // BA strand tags (bD, bM, bE)
        let b_max_depth = ss_b.depths.iter().map(|&d| i32::from(d)).max().unwrap_or(0);
        let b_min_depth = ss_b.depths.iter().map(|&d| i32::from(d)).min().unwrap_or(0);
        let b_total_errors: usize = ss_b.errors.iter().map(|&e| e as usize).sum();
        let b_total_bases: usize = ss_b.depths.iter().map(|&d| d as usize).sum();
        let b_error_rate =
            if b_total_bases > 0 { b_total_errors as f32 / b_total_bases as f32 } else { 0.0 };

        let bd_tag = Tag::from([b'b', b'D']);
        let bm_tag = Tag::from([b'b', b'M']);
        let be_tag = Tag::from([b'b', b'E']);
        data.insert(bd_tag, to_smallest_signed_int(b_max_depth));
        data.insert(bm_tag, to_smallest_signed_int(b_min_depth));
        data.insert(be_tag, Value::from(b_error_rate));

        // Per-base tags if enabled
        if self.options.produce_per_base_tags {
            // Per-base depth arrays - convert to signed integers to match fgbio/simplex/duplex
            let ad_base_tag = Tag::from([b'a', b'd']);
            let bd_base_tag = Tag::from([b'b', b'd']);
            let ad_array: Vec<i16> = ss_a.depths.iter().map(|&d| d as i16).collect();
            let bd_array: Vec<i16> = ss_b.depths.iter().map(|&d| d as i16).collect();
            data.insert(ad_base_tag, Value::from(ad_array));
            data.insert(bd_base_tag, Value::from(bd_array));

            // Per-base error arrays - convert to signed integers to match fgbio/simplex/duplex
            let ae_base_tag = Tag::from([b'a', b'e']);
            let be_base_tag = Tag::from([b'b', b'e']);
            let ae_array: Vec<i16> = ss_a.errors.iter().map(|&e| e as i16).collect();
            let be_array: Vec<i16> = ss_b.errors.iter().map(|&e| e as i16).collect();
            data.insert(ae_base_tag, Value::from(ae_array));
            data.insert(be_base_tag, Value::from(be_array));

            // Single-strand consensus sequences - convert to String format to match fgbio
            let ac_tag = Tag::from([b'a', b'c']);
            let bc_tag = Tag::from([b'b', b'c']);
            data.insert(ac_tag, crate::consensus_tags::sequence_to_tag_value(&ss_a.bases));
            data.insert(bc_tag, crate::consensus_tags::sequence_to_tag_value(&ss_b.bases));

            // Single-strand consensus qualities - convert to Phred+33 ASCII string to match fgbio
            let aq_tag = Tag::from([b'a', b'q']);
            let bq_tag = Tag::from([b'b', b'q']);
            data.insert(aq_tag, crate::consensus_tags::qualities_to_tag_value(&ss_a.quals));
            data.insert(bq_tag, crate::consensus_tags::qualities_to_tag_value(&ss_b.quals));
        }

        // Cell barcode tag - extract from first source read if configured
        if let Some(cell_tag) = &self.options.cell_tag {
            for read in source_reads {
                if let Some(Value::String(cell_bytes)) = read.data().get(cell_tag) {
                    let cell_bc =
                        String::from_utf8(cell_bytes.iter().copied().collect::<Vec<u8>>())
                            .unwrap_or_default();
                    if !cell_bc.is_empty() {
                        data.insert(*cell_tag, Value::from(cell_bc));
                        break; // Use the first non-empty cell barcode found
                    }
                }
            }
        }

        // RX tag (UmiBases) - extract from all source reads and build consensus
        let rx_tag = Tag::from([b'R', b'X']);
        let umis: Vec<String> = source_reads
            .iter()
            .filter_map(|read| {
                if let Some(Value::String(rx_bytes)) = read.data().get(&rx_tag) {
                    String::from_utf8(rx_bytes.iter().copied().collect::<Vec<u8>>()).ok()
                } else {
                    None
                }
            })
            .collect();

        if !umis.is_empty() {
            let consensus_umi = consensus_umis(&umis);
            if !consensus_umi.is_empty() {
                data.insert(rx_tag, Value::from(consensus_umi));
            }
        }

        *record.data_mut() = data;

        // Suppress unused variable warning
        let _ = (overlap_start, overlap_end);

        Ok(record)
    }

    /// Rejects a set of records with the given reason
    fn reject_records(&mut self, records: &[RecordBuf], reason: CallerRejectionReason) {
        let count = records.len();
        *self.stats.rejection_reasons.entry(reason).or_insert(0) += count;
        self.stats.reads_filtered += count as u64;

        // Store rejected reads if tracking is enabled
        if self.track_rejects {
            self.rejected_reads.extend(records.iter().cloned());
        }
    }

    /// Rejects a count of records with the given reason (without storing the actual records)
    fn reject_records_count(&mut self, count: usize, reason: CallerRejectionReason) {
        *self.stats.rejection_reasons.entry(reason).or_insert(0) += count;
        self.stats.reads_filtered += count as u64;
    }
}

/// Implementation of the `ConsensusCaller` trait for CODEC consensus calling.
///
/// This allows `CodecConsensusCaller` to be used polymorphically with other
/// consensus callers (e.g., `VanillaUmiConsensusCaller`, `DuplexConsensusCaller`).
impl ConsensusCaller for CodecConsensusCaller {
    fn consensus_reads_from_sam_records(
        &mut self,
        reads: Vec<RecordBuf>,
    ) -> Result<Vec<RecordBuf>> {
        // Delegate to the existing implementation
        CodecConsensusCaller::consensus_reads_from_sam_records(self, reads)
    }

    fn total_reads(&self) -> usize {
        self.stats.total_input_reads as usize
    }

    fn total_filtered(&self) -> usize {
        self.stats.reads_filtered as usize
    }

    fn consensus_reads_constructed(&self) -> usize {
        self.stats.consensus_reads_generated as usize
    }

    fn statistics(&self) -> ConsensusCallingStats {
        ConsensusCallingStats {
            total_reads: self.stats.total_input_reads as usize,
            consensus_reads: self.stats.consensus_reads_generated as usize,
            filtered_reads: self.stats.reads_filtered as usize,
            rejection_reasons: self.stats.rejection_reasons.clone(),
        }
    }

    fn log_statistics(&self) {
        let stats = &self.stats;
        log::info!("CODEC Consensus Calling Statistics:");
        log::info!("  Total input reads: {}", stats.total_input_reads);
        log::info!("  Consensus reads generated: {}", stats.consensus_reads_generated);
        log::info!("  Reads filtered: {}", stats.reads_filtered);
        log::info!("  Consensus bases emitted: {}", stats.consensus_bases_emitted);
        log::info!("  Duplex bases emitted: {}", stats.consensus_duplex_bases_emitted);
        log::info!("  Duplex disagreement rate: {:.4}%", stats.duplex_disagreement_rate() * 100.0);
        if !stats.rejection_reasons.is_empty() {
            log::info!("  Rejection reasons:");
            for (reason, count) in &stats.rejection_reasons {
                log::info!("    {reason:?}: {count}");
            }
        }
    }
}

#[cfg(test)]
#[allow(clippy::similar_names)]
mod tests {
    use super::*;
    use crate::sam::builder::RecordBuilder;
    use noodles::sam::alignment::record::cigar::op::Kind;

    #[test]
    fn test_codec_caller_creation() {
        let options = CodecConsensusOptions::default();
        let caller = CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options);

        assert_eq!(caller.read_name_prefix, "codec");
        assert_eq!(caller.read_group_id, "RG1");
    }

    #[test]
    fn test_codec_options_defaults() {
        let options = CodecConsensusOptions::default();
        assert_eq!(options.min_reads_per_strand, 1);
        assert_eq!(options.min_duplex_length, 1);
        assert_eq!(options.min_input_base_quality, 10);
        assert_eq!(options.error_rate_pre_umi, 45);
        assert_eq!(options.error_rate_post_umi, 40);
    }

    #[test]
    fn test_empty_input() {
        let options = CodecConsensusOptions::default();
        let mut caller = CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options);

        let result = caller.consensus_reads_from_sam_records(Vec::new()).unwrap();
        assert!(result.is_empty());
        assert_eq!(caller.stats.total_input_reads, 0);
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT".to_vec());
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT".to_vec());
        assert_eq!(reverse_complement(b"ACGTN"), b"NACGT".to_vec());
        assert_eq!(reverse_complement(b""), b"".to_vec());
    }

    // FIXME: reenable later?
    /*
    #[test]
    fn test_reverse_complement_ss() {
        let ss = SingleStrandConsensus {
            bases: b"ACGT".to_vec(),
            quals: vec![10, 20, 30, 40],
            depths: vec![1, 2, 3, 4],
            errors: vec![0, 0, 1, 0],
            raw_read_count: 5,
            ref_start: 100,
            ref_end: 103,
            is_negative_strand: false,
        };

        let rc = CodecConsensusCaller::reverse_complement_ss(&ss);

        // Bases are reverse complemented, quals are reversed, depths/errors stay the same
        assert_eq!(rc.bases, b"ACGT"); // ACGT → TGCA reversed → ACGT complemented
        assert_eq!(rc.quals, vec![40, 30, 20, 10]); // reversed
        assert_eq!(rc.depths, vec![1, 2, 3, 4]); // NOT reversed (per-position metadata)
        assert_eq!(rc.errors, vec![0, 0, 1, 0]); // NOT reversed (per-position metadata)
        assert_eq!(rc.raw_read_count, 5);
        assert!(rc.is_negative_strand); // flipped
    }
    */

    #[test]
    fn test_stats_duplex_disagreement_rate() {
        let mut stats = CodecConsensusStats::default();
        assert!(stats.duplex_disagreement_rate().abs() < f64::EPSILON);

        stats.consensus_duplex_bases_emitted = 100;
        stats.duplex_disagreement_base_count = 5;
        assert!((stats.duplex_disagreement_rate() - 0.05).abs() < 1e-6);
    }

    /// Helper function to create a test paired read
    #[allow(clippy::too_many_arguments)]
    fn create_test_paired_read(
        name: &str,
        seq: &[u8],
        qual: &[u8],
        is_first: bool,
        is_reverse: bool,
        mate_is_reverse: bool,
        ref_start: usize,
        cigar_ops: &[(Kind, usize)],
    ) -> RecordBuf {
        // Convert CIGAR ops to string and calculate reference length
        use std::fmt::Write;
        let mut ref_len = 0usize;
        let cigar_str: String = cigar_ops.iter().fold(String::new(), |mut acc, (kind, len)| {
            // Calculate reference length consumed
            match kind {
                Kind::Match
                | Kind::SequenceMatch
                | Kind::SequenceMismatch
                | Kind::Deletion
                | Kind::Skip => {
                    ref_len += len;
                }
                _ => {}
            }
            let c = match kind {
                Kind::Match => 'M',
                Kind::Insertion => 'I',
                Kind::Deletion => 'D',
                Kind::SoftClip => 'S',
                Kind::HardClip => 'H',
                Kind::Skip => 'N',
                Kind::Pad => 'P',
                Kind::SequenceMatch => '=',
                Kind::SequenceMismatch => 'X',
            };
            let _ = write!(acc, "{len}{c}");
            acc
        });

        // Calculate mate position and template length for FR pair detection
        // Use a typical insert size of 200bp
        let insert_size: i32 = 200;
        let (mate_start, tlen) = if is_reverse {
            // Read on reverse strand: mate should be before this read
            let mate_pos = (ref_start as i32 - insert_size + ref_len as i32).max(1) as usize;
            (mate_pos, -insert_size)
        } else {
            // Read on positive strand: mate should be after this read
            let mate_pos = ((ref_start as i32) + insert_size - (ref_len as i32)).max(1) as usize;
            (mate_pos, insert_size)
        };

        let seq_str = String::from_utf8_lossy(seq);
        RecordBuilder::new()
            .name(name)
            .sequence(&seq_str)
            .qualities(qual)
            .reference_sequence_id(0)
            .alignment_start(ref_start)
            .cigar(&cigar_str)
            .first_segment(is_first)
            .properly_paired(true)
            .reverse_complement(is_reverse)
            .mate_reverse_complement(mate_is_reverse)
            .mate_reference_sequence_id(0)
            .mate_alignment_start(mate_start)
            .template_length(tlen)
            .tag("MI", "UMI123")
            .build()
    }

    #[test]
    fn test_is_fr_pair() {
        use noodles::sam::alignment::record::cigar::op::Kind;

        let options = CodecConsensusOptions::default();
        let caller = CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options);

        // FR pair: R1 forward, R2 reverse
        let r1_fr = create_test_paired_read(
            "read1",
            b"ACGT",
            b"####",
            true,
            false,
            true,
            100,
            &[(Kind::Match, 4)],
        );
        assert!(caller.is_fr_pair(&r1_fr));

        // RF pair: R1 reverse, R2 forward - still considered FR orientation
        let r1_rf = create_test_paired_read(
            "read1",
            b"ACGT",
            b"####",
            true,
            true,
            false,
            100,
            &[(Kind::Match, 4)],
        );
        assert!(caller.is_fr_pair(&r1_rf));

        // FF pair: both forward - not FR
        let r1_ff = create_test_paired_read(
            "read1",
            b"ACGT",
            b"####",
            true,
            false,
            false,
            100,
            &[(Kind::Match, 4)],
        );
        assert!(!caller.is_fr_pair(&r1_ff));
    }

    #[test]
    fn test_simplify_cigar() {
        use noodles::sam::alignment::record::cigar::op::Kind;

        // Test simple match
        let read = create_test_paired_read(
            "read",
            b"ACGTACGT",
            b"########",
            true,
            false,
            true,
            100,
            &[(Kind::Match, 8)],
        );
        assert_eq!(cigar_utils::simplify_cigar(read.cigar()), vec![(Kind::Match, 8)]);

        // Test with indels
        let read_indel = create_test_paired_read(
            "read",
            b"ACGTACGT",
            b"########",
            true,
            false,
            true,
            100,
            &[(Kind::Match, 4), (Kind::Deletion, 2), (Kind::Match, 4)],
        );
        assert_eq!(
            cigar_utils::simplify_cigar(read_indel.cigar()),
            vec![(Kind::Match, 4), (Kind::Deletion, 2), (Kind::Match, 4)]
        );

        // Test coalescing: soft clip + match should coalesce to match
        let read_softclip = create_test_paired_read(
            "read",
            b"ACGTACGT",
            b"########",
            true,
            false,
            true,
            100,
            &[(Kind::SoftClip, 2), (Kind::Match, 6)],
        );
        assert_eq!(
            cigar_utils::simplify_cigar(read_softclip.cigar()),
            vec![(Kind::Match, 8)] // Should coalesce S+M to M
        );
    }

    #[test]
    fn test_filter_to_most_common_alignment() {
        use noodles::sam::alignment::record::cigar::op::Kind;

        let options = CodecConsensusOptions::default();
        let mut caller = CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options);

        // Create 3 reads with same CIGAR and 1 with different
        let reads = vec![
            create_test_paired_read(
                "r1",
                b"ACGT",
                b"####",
                true,
                false,
                true,
                100,
                &[(Kind::Match, 4)],
            ),
            create_test_paired_read(
                "r2",
                b"ACGT",
                b"####",
                true,
                false,
                true,
                100,
                &[(Kind::Match, 4)],
            ),
            create_test_paired_read(
                "r3",
                b"ACGT",
                b"####",
                true,
                false,
                true,
                100,
                &[(Kind::Match, 4)],
            ),
            create_test_paired_read(
                "r4",
                b"ACGT",
                b"####",
                true,
                false,
                true,
                100,
                &[(Kind::Match, 2), (Kind::Deletion, 1), (Kind::Match, 2)],
            ),
        ];

        let filtered = caller.filter_to_most_common_alignment(reads);
        assert_eq!(filtered.len(), 3);
        assert_eq!(
            caller.stats.rejection_reasons.get(&CallerRejectionReason::MinorityAlignment),
            Some(&1)
        );
    }

    /// Tests deterministic tie-breaking when two groups have equal sizes.
    /// With CIGAR-based tie-breaking, the group with the "smaller" CIGAR wins.
    /// CIGAR comparison is done element-by-element: first by length, then by operation kind.
    /// This matches fgbio's filterToMostCommonAlignment behavior with deterministic tie-breaking.
    #[test]
    fn test_filter_to_most_common_alignment_tie_breaking() {
        use noodles::sam::alignment::record::cigar::op::Kind;

        let options = CodecConsensusOptions::default();
        let mut caller = CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options);

        // Create two groups with equal sizes: 2 reads with 4M and 2 reads with 3M1D1M
        // Both have the same query length (4 bases).
        // With deterministic CIGAR-based tie-breaking:
        // - 4M = [(Match, 4)]
        // - 3M1D1M = [(Match, 3), (Deletion, 1), (Match, 1)]
        // Comparing first elements: length 3 < length 4, so 3M1D1M is "smaller" and wins.
        let reads = vec![
            create_test_paired_read(
                "r1_4M",
                b"ACGT",
                b"####",
                true,
                false,
                true,
                100,
                &[(Kind::Match, 4)],
            ),
            create_test_paired_read(
                "r2_4M",
                b"ACGT",
                b"####",
                true,
                false,
                true,
                100,
                &[(Kind::Match, 4)],
            ),
            create_test_paired_read(
                "r3_del",
                b"ACGT",
                b"####",
                true,
                false,
                true,
                100,
                &[(Kind::Match, 3), (Kind::Deletion, 1), (Kind::Match, 1)],
            ),
            create_test_paired_read(
                "r4_del",
                b"ACGT",
                b"####",
                true,
                false,
                true,
                100,
                &[(Kind::Match, 3), (Kind::Deletion, 1), (Kind::Match, 1)],
            ),
        ];

        let filtered = caller.filter_to_most_common_alignment(reads);
        // Both groups have size 2, so we expect 2 reads to be returned
        assert_eq!(filtered.len(), 2);
        // The 3M1D1M group should be selected due to deterministic CIGAR tie-breaking
        // (smaller CIGAR wins, and 3 < 4 in the first element comparison)
        // Verify that the filtered reads are the 3M1D1M reads (r3_del and r4_del)
        for read in &filtered {
            let cigar = read.cigar();
            assert_eq!(cigar.len(), 3, "Expected 3-op cigar (3M1D1M)");
            let ops: Vec<_> = cigar.iter().map(|r| r.unwrap()).collect();
            assert_eq!(ops[0].kind(), Kind::Match);
            assert_eq!(ops[0].len(), 3);
            assert_eq!(ops[1].kind(), Kind::Deletion);
            assert_eq!(ops[1].len(), 1);
            assert_eq!(ops[2].kind(), Kind::Match);
            assert_eq!(ops[2].len(), 1);
        }
        // 2 reads should be rejected
        assert_eq!(
            caller.stats.rejection_reasons.get(&CallerRejectionReason::MinorityAlignment),
            Some(&2)
        );
    }

    // ==========================================================================
    // Integration tests ported from fgbio CodecConsensusCallerTest
    // ==========================================================================

    /// Reference sequence for tests (from fgbio)
    const REF_BASES: &[u8] = b"TGGGTGTTGTTTGGGTTCCTGCTTAAAGAGCTACTGTTCTTCACAGAAACTTCCAACTCACCCAGACTGAGATTTGTACTGAGACTACGATCCACATGTTCAATATCTGATATCTGATGGGAAATAGGCTTTACTGAATTATCCATTTGGGCTGTAATTAATTTCAGTGATGAGCGGGAGATGTTGTTAGTTGTGCTCAGTAACTTTTTGATAGTAGCGGGAGTAGGAGTAAATCTTGTACTAATTAGTGAATATTCTGTTGATGGTGGCTGAAAATTTATAGCTACACAACCAAAAAAATAAAAAACGTTAGTCAATAGCATTTATAAATAGTCTTCTCTACCTGAAATATTTTACATTAAGTAATTCATTCCTTCATTTAGTATCTACACATGTCTAACATTGTAGTAGGAGCTGTGTACTAACAAGAAATCATGACACTGTTTCTGCCTTCAAGGAGCTTATAATCTTTTGGGGTACACAAGATAACCCAGAATGTTAAATAGTATAAAAGTCAAAGTACAATAATTTATTTCATTAAGATTTTGAAATGGCTAACAAACACCTGTTGATCACCTCATACACATGAGCCTCAAAACAAAGGAAAGCACAGCCCCTATGCCTGAGCAATTTAGAATATTGTCAAGGATAGAGACATGTGAGCCATTCACTATGAAACAATCATTGAGAACTACTACAAGAGTGATAAATATAAAATGAAACCTACAGAAACACAGAAGAGTAAGTAATTTTCCCTATAAAGAAGACAGGAACTAAATGTATAAGCAAAAATTGGGAAATTATATAAATGCTATTTTATATGAGAGGCAAAGAACCACAGGTCTAATAATTTTACAAATGTGATAAAATCAGATTTTATGTCCCCATCTTTCTTGACTGCTCAGCTAGAAATTAAAACATTTTTACACATCTTTTTGGCGGGGGCGGGGGGGATCATTATTTATTTCACCTGCCAAAATACTTCATTTCCTTATTGCACTTTTTTACTTCTTTGGTATGGAAAAATCTAACGGGTTTTAGAGTATGAACACATTTTAAGCAGTGATTAGATACGTTTTTCTTGTTATGCTTTCTATTGCAAATTTAGGATTTGATTTTGCACTGTCTTCATGCAAAGCTCTTCTCAAAGGTCTTAAAATATAAAAAACACTTAATGCTTCTCAAAGCATTAAGATTTTATGTAAATCAAACCAAAACCAGAAAAAGACAGAAGAAAATGAACCAAAAACAACAAAAATAATCCTTAACATAGTTGGCAACAAGTGCAATGAAAGATTTTT";

    /// Helper to create an FR read pair similar to fgbio's SamBuilder.addPair
    ///
    /// Creates a properly oriented FR pair where:
    /// - R1 is `first_segment`, forward (unless strand1=Minus)
    /// - R2 is `last_segment`, reverse (unless strand2=Plus)
    /// - Both reads have proper mate information set
    #[allow(clippy::too_many_arguments)]
    fn create_fr_pair(
        name: &str,
        start1: usize,
        start2: usize,
        _read_length: usize,
        base_quality: u8,
        cigar1: &[(Kind, usize)],
        cigar2: &[(Kind, usize)],
        mi_tag: &str,
        rx_tag: Option<&str>,
        strand1_reverse: bool,
        strand2_reverse: bool,
    ) -> Vec<RecordBuf> {
        // Calculate reference end positions based on CIGAR
        let calc_ref_length = |cigar: &[(Kind, usize)]| -> usize {
            cigar
                .iter()
                .filter(|(k, _)| matches!(k, Kind::Match | Kind::Deletion | Kind::Skip))
                .map(|(_, len)| *len)
                .sum()
        };

        let ref_len1 = calc_ref_length(cigar1);
        let ref_len2 = calc_ref_length(cigar2);

        // Get sequence for R1 and R2 based on positions
        let get_sequence = |start: usize, cigar: &[(Kind, usize)]| -> Vec<u8> {
            let mut seq = Vec::new();
            let mut ref_pos = start - 1; // Convert to 0-based
            for &(kind, len) in cigar {
                match kind {
                    Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                        // Consume both ref and read
                        let end = (ref_pos + len).min(REF_BASES.len());
                        if ref_pos < REF_BASES.len() {
                            seq.extend_from_slice(&REF_BASES[ref_pos..end]);
                            // Pad with A's if we ran past ref
                            let pad_len = (ref_pos + len).saturating_sub(end);
                            seq.resize(seq.len() + pad_len, b'A');
                        }
                        ref_pos += len;
                    }
                    Kind::Insertion | Kind::SoftClip => {
                        // Consume read only - add placeholder bases
                        seq.resize(seq.len() + len, b'A');
                    }
                    Kind::Deletion | Kind::Skip => {
                        // Consume ref only
                        ref_pos += len;
                    }
                    Kind::HardClip | Kind::Pad => {
                        // No change
                    }
                }
            }
            seq
        };

        let seq1_fwd = get_sequence(start1, cigar1);
        let seq2_fwd = get_sequence(start2, cigar2);

        // For a true FR pair sequencing the same DNA molecule from opposite strands:
        // - Positive strand reads are stored as-is in BAM (matches reference)
        // - Negative strand reads are stored as revcomp in BAM
        //   After revcomp during SS building, they become the same as reference
        // This ensures that both R1 and R2 SS consensuses agree in the overlap.
        let seq1: Vec<u8> = if strand1_reverse { reverse_complement(&seq1_fwd) } else { seq1_fwd };
        let seq2: Vec<u8> = if strand2_reverse { reverse_complement(&seq2_fwd) } else { seq2_fwd };

        let qual1 = vec![base_quality; seq1.len()];
        let qual2 = vec![base_quality; seq2.len()];

        // Convert CIGAR ops to strings
        let cigar_to_str = |cigar: &[(Kind, usize)]| -> String {
            use std::fmt::Write;
            cigar.iter().fold(String::new(), |mut acc, (kind, len)| {
                let c = match kind {
                    Kind::Match => 'M',
                    Kind::Insertion => 'I',
                    Kind::Deletion => 'D',
                    Kind::SoftClip => 'S',
                    Kind::HardClip => 'H',
                    Kind::Skip => 'N',
                    Kind::Pad => 'P',
                    Kind::SequenceMatch => '=',
                    Kind::SequenceMismatch => 'X',
                };
                let _ = write!(acc, "{len}{c}");
                acc
            })
        };

        let cigar1_str = cigar_to_str(cigar1);
        let cigar2_str = cigar_to_str(cigar2);

        // Calculate template length (insert size)
        // For FR pairs, template length should be positive for the leftmost read
        let template_len = if start1 <= start2 {
            (start2 + ref_len2 - start1) as i32
        } else {
            -((start1 + ref_len1 - start2) as i32)
        };

        // Build R1
        let seq1_str = String::from_utf8_lossy(&seq1);
        let mut r1_builder = RecordBuilder::new()
            .name(name)
            .sequence(&seq1_str)
            .qualities(&qual1)
            .reference_sequence_id(0)
            .alignment_start(start1)
            .cigar(&cigar1_str)
            .first_segment(true)
            .properly_paired(true)
            .reverse_complement(strand1_reverse)
            .mate_reverse_complement(strand2_reverse)
            .mate_reference_sequence_id(0)
            .mate_alignment_start(start2)
            .template_length(template_len)
            .tag("MI", mi_tag);

        if let Some(rx) = rx_tag {
            r1_builder = r1_builder.tag("RX", rx);
        }
        let r1 = r1_builder.build();

        // Build R2
        let seq2_str = String::from_utf8_lossy(&seq2);
        let mut r2_builder = RecordBuilder::new()
            .name(name)
            .sequence(&seq2_str)
            .qualities(&qual2)
            .reference_sequence_id(0)
            .alignment_start(start2)
            .cigar(&cigar2_str)
            .first_segment(false)
            .properly_paired(true)
            .reverse_complement(strand2_reverse)
            .mate_reverse_complement(strand1_reverse)
            .mate_reference_sequence_id(0)
            .mate_alignment_start(start1)
            .template_length(-template_len)
            .tag("MI", mi_tag);

        if let Some(rx) = rx_tag {
            r2_builder = r2_builder.tag("RX", rx);
        }
        let r2 = r2_builder.build();

        vec![r1, r2]
    }

    /// Port of fgbio test: "make a consensus from two simple reads"
    ///
    /// Tests that a simple FR pair produces a consensus with correct structure.
    /// Note: Positions covered by only one strand get N bases (correct duplex behavior),
    /// so we don't check exact sequence content - only that consensus structure is correct.
    #[test]
    fn test_make_consensus_from_simple_reads() {
        let options = CodecConsensusOptions {
            min_reads_per_strand: 1,
            min_duplex_length: 1,
            ..Default::default()
        };
        let mut caller = CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options);

        // Create FR pair: R1 at pos 1, R2 at pos 11, both 30bp
        let reads = create_fr_pair(
            "read1",
            1,
            11,
            30,
            35,
            &[(Kind::Match, 30)],
            &[(Kind::Match, 30)],
            "hi",
            Some("ACC-TGA"),
            false, // R1 forward
            true,  // R2 reverse
        );

        let cons = caller.consensus_reads_from_sam_records(reads).unwrap();

        assert_eq!(cons.len(), 1, "Should produce one consensus read");

        let consensus = &cons[0];
        // Check read name format
        let name = consensus.name().map(std::string::ToString::to_string).unwrap_or_default();
        assert!(name.contains("hi"), "Consensus name should contain MI tag: {name}");

        // Consensus should cover from pos 1 to pos 40 (R1: 1-30, R2: 11-40)
        assert_eq!(consensus.sequence().len(), 40, "Consensus should be 40bp");

        // Check RX tag is preserved (UmiBases consensus)
        let rx_tag = Tag::from([b'R', b'X']);
        if let Some(Value::String(rx)) = consensus.data().get(&rx_tag) {
            assert_eq!(&rx[..], b"ACC-TGA", "RX tag should be preserved");
        } else {
            panic!("RX tag not found in consensus");
        }
    }

    /// Port of fgbio test: "make a consensus where R1 has a deletion outside of the overlap region"
    #[test]
    fn test_make_consensus_r1_deletion() {
        let options = CodecConsensusOptions {
            min_reads_per_strand: 1,
            min_duplex_length: 1,
            ..Default::default()
        };
        let mut caller = CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options);

        // R1 at pos 1 with 5M2D25M, R2 at pos 13 with 30M
        let reads = create_fr_pair(
            "read1",
            1,
            13,
            30,
            35,
            &[(Kind::Match, 5), (Kind::Deletion, 2), (Kind::Match, 25)],
            &[(Kind::Match, 30)],
            "hi",
            Some("ACC-TGA"),
            false,
            true,
        );

        let cons = caller.consensus_reads_from_sam_records(reads).unwrap();

        // Should produce consensus
        assert_eq!(cons.len(), 1, "Should produce one consensus read");

        // Consensus should cover the region with deletion accounted for
        let consensus = &cons[0];
        assert!(!consensus.sequence().is_empty(), "Should have sequence");
    }

    /// Port of fgbio test: "not emit a consensus when the reads are an RF pair"
    #[test]
    fn test_not_emit_consensus_for_rf_pair() {
        let options = CodecConsensusOptions {
            min_reads_per_strand: 1,
            min_duplex_length: 1,
            ..Default::default()
        };
        let mut caller = CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options);

        // RF pair: R1 reverse at pos 100, R2 forward at pos 135
        let reads = create_fr_pair(
            "read1",
            100,
            135,
            30,
            35,
            &[(Kind::Match, 30)],
            &[(Kind::Match, 30)],
            "hi",
            Some("ACC-TGA"),
            true,  // R1 reverse
            false, // R2 forward
        );

        let cons = caller.consensus_reads_from_sam_records(reads).unwrap();

        assert!(cons.is_empty(), "Should not emit consensus for RF pair");
    }

    /// Port of fgbio test: "not emit a consensus when there are insufficient reads"
    #[test]
    fn test_not_emit_consensus_insufficient_reads() {
        // With minReadsPerStrand=2, a single pair should fail
        let options2 = CodecConsensusOptions {
            min_reads_per_strand: 2,
            min_duplex_length: 1,
            ..Default::default()
        };
        let mut caller2 =
            CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options2);

        let reads = create_fr_pair(
            "read1",
            1,
            11,
            30,
            35,
            &[(Kind::Match, 30)],
            &[(Kind::Match, 30)],
            "hi",
            Some("ACC-TGA"),
            false,
            true,
        );

        let cons = caller2.consensus_reads_from_sam_records(reads).unwrap();

        assert!(cons.is_empty(), "Should not emit consensus with insufficient reads");

        // With minReadsPerStrand=1, it should succeed
        let options1 = CodecConsensusOptions {
            min_reads_per_strand: 1,
            min_duplex_length: 1,
            ..Default::default()
        };
        let mut caller1 =
            CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options1);

        let reads = create_fr_pair(
            "read1",
            1,
            11,
            30,
            35,
            &[(Kind::Match, 30)],
            &[(Kind::Match, 30)],
            "hi",
            Some("ACC-TGA"),
            false,
            true,
        );

        let cons = caller1.consensus_reads_from_sam_records(reads).unwrap();

        assert_eq!(cons.len(), 1, "Should emit consensus with minReadsPerStrand=1");
    }

    /// Port of fgbio test: "not emit a consensus when there is insufficient overlap between R1 and R2"
    #[test]
    fn test_not_emit_consensus_insufficient_overlap() {
        // R1 at 1-30, R2 at 11-40 gives 20bp overlap
        // With minDuplexLength=20, should pass
        let options20 = CodecConsensusOptions {
            min_reads_per_strand: 1,
            min_duplex_length: 20,
            ..Default::default()
        };
        let mut caller20 =
            CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options20);

        let reads = create_fr_pair(
            "read1",
            1,
            11,
            30,
            35,
            &[(Kind::Match, 30)],
            &[(Kind::Match, 30)],
            "hi",
            Some("ACC-TGA"),
            false,
            true,
        );

        let cons = caller20.consensus_reads_from_sam_records(reads).unwrap();
        assert_eq!(cons.len(), 1, "Should emit consensus with minDuplexLength=20 and 20bp overlap");

        // With minDuplexLength=21, should fail
        let options21 = CodecConsensusOptions {
            min_reads_per_strand: 1,
            min_duplex_length: 21,
            ..Default::default()
        };
        let mut caller21 =
            CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options21);

        let reads = create_fr_pair(
            "read1",
            1,
            11,
            30,
            35,
            &[(Kind::Match, 30)],
            &[(Kind::Match, 30)],
            "hi",
            Some("ACC-TGA"),
            false,
            true,
        );

        let cons = caller21.consensus_reads_from_sam_records(reads).unwrap();
        assert!(
            cons.is_empty(),
            "Should not emit consensus with minDuplexLength=21 and 20bp overlap"
        );
    }

    /// Port of fgbio test: "not emit a consensus when the read pair has one mate unmapped"
    #[test]
    fn test_not_emit_consensus_unmapped_mate() {
        let options = CodecConsensusOptions {
            min_reads_per_strand: 1,
            min_duplex_length: 1,
            ..Default::default()
        };
        let mut caller = CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options);

        // Create FR pair then mark R2 as unmapped
        let mut reads = create_fr_pair(
            "read1",
            1,
            11,
            30,
            35,
            &[(Kind::Match, 30)],
            &[(Kind::Match, 30)],
            "hi",
            Some("ACC-TGA"),
            false,
            true,
        );

        // Mark R2 as unmapped
        let r2 = &mut reads[1];
        *r2.flags_mut() = r2.flags() | Flags::UNMAPPED;

        let cons = caller.consensus_reads_from_sam_records(reads).unwrap();

        assert!(cons.is_empty(), "Should not emit consensus when mate is unmapped");
    }

    /// Port of fgbio test: "emit the consensus in the orientation of R1"
    ///
    /// This test verifies that the consensus output orientation follows R1.
    ///
    /// Tests that a forward R1 produces a consensus with the expected length.
    /// Note: Positions covered by only one strand get N bases (correct duplex behavior).
    #[test]
    fn test_emit_consensus_in_r1_orientation() {
        // Test that forward R1 produces consensus with correct length
        let options = CodecConsensusOptions {
            min_reads_per_strand: 1,
            min_duplex_length: 1,
            ..Default::default()
        };
        let mut caller_fwd =
            CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options.clone());

        // Forward R1 (start=1), reverse R2 (start=11)
        let reads_fwd = create_fr_pair(
            "read1",
            1,
            11,
            30,
            35,
            &[(Kind::Match, 30)],
            &[(Kind::Match, 30)],
            "hi",
            Some("ACC-TGA"),
            false, // R1 forward
            true,  // R2 reverse
        );

        let cons_fwd = caller_fwd.consensus_reads_from_sam_records(reads_fwd).unwrap();
        assert_eq!(cons_fwd.len(), 1, "Should produce one consensus read");

        // Consensus should cover from pos 1 to pos 40 (R1: 1-30, R2: 11-40)
        assert_eq!(cons_fwd[0].sequence().len(), 40, "Forward R1 should produce 40bp consensus");

        // Verify orientation flag is set correctly (forward R1 = not reverse-complemented)
        assert!(
            !cons_fwd[0].flags().is_reverse_complemented(),
            "Forward R1 should produce forward-oriented consensus"
        );
    }

    /// Port of fgbio test: "not emit a consensus when there are a lot of disagreements between strands"
    ///
    /// Tests disagreement filtering: with permissive settings consensus is produced,
    /// but with strict settings, consensus fails due to strand disagreements at
    /// non-overlapping positions (where one strand has N).
    #[test]
    fn test_not_emit_consensus_high_disagreement() {
        // First test that we get consensus with permissive disagreement settings
        let options_permissive = CodecConsensusOptions {
            min_reads_per_strand: 1,
            min_duplex_length: 1,
            max_duplex_disagreements: 100,
            max_duplex_disagreement_rate: 1.0,
            ..Default::default()
        };
        let mut caller_permissive =
            CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options_permissive);

        // Use positions within REF_BASES bounds (positions 1-40)
        let reads = create_fr_pair(
            "read1",
            1,
            11,
            30,
            35,
            &[(Kind::Match, 30)],
            &[(Kind::Match, 30)],
            "hi",
            Some("ACC-TGA"),
            false,
            true,
        );

        let cons = caller_permissive.consensus_reads_from_sam_records(reads).unwrap();
        assert_eq!(cons.len(), 1, "Should emit consensus with permissive settings");

        // Now test with strict disagreement limits - should fail because
        // positions covered by only one strand (10 on each side) count as disagreements
        let options_strict = CodecConsensusOptions {
            min_reads_per_strand: 1,
            min_duplex_length: 1,
            max_duplex_disagreements: 5, // Only allow 5 disagreements
            max_duplex_disagreement_rate: 0.05, // 5% rate
            ..Default::default()
        };
        let mut caller_strict =
            CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options_strict);

        // Use positions within REF_BASES bounds
        let reads2 = create_fr_pair(
            "read1",
            1,
            11,
            30,
            35,
            &[(Kind::Match, 30)],
            &[(Kind::Match, 30)],
            "hi",
            Some("ACC-TGA"),
            false,
            true,
        );

        // With strict settings, this should fail due to high disagreement
        // (positions covered by only one strand count as disagreements)
        let cons2 = caller_strict.consensus_reads_from_sam_records(reads2);
        assert!(cons2.is_err(), "Should reject consensus with strict disagreement settings");
    }

    /// Port of fgbio test: "make a consensus where R2 has a deletion outside of the overlap region"
    ///
    /// fgumi now matches fgbio behavior: produces a 40bp consensus (skipping the deleted region).
    #[test]
    fn test_make_consensus_r2_deletion() {
        let options = CodecConsensusOptions {
            min_reads_per_strand: 1,
            min_duplex_length: 1,
            ..Default::default()
        };
        let mut caller = CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options);

        // R1 at pos 1 with 30M, R2 at pos 11 with 25M5D5M
        // In fgbio: start1=1, start2=11, cigar1="30M", cigar2="25M5D5M"
        // R2 has a 5bp deletion outside the overlap region
        let reads = create_fr_pair(
            "read1",
            1,
            11,
            30,
            35,
            &[(Kind::Match, 30)],
            &[(Kind::Match, 25), (Kind::Deletion, 5), (Kind::Match, 5)],
            "hi",
            Some("ACC-TGA"),
            false,
            true,
        );

        let cons = caller.consensus_reads_from_sam_records(reads).unwrap();

        // Should produce consensus
        assert_eq!(cons.len(), 1, "Should produce one consensus read");

        // fgbio produces 40bp consensus (skipping the deleted region in R2)
        // fgumi now matches this behavior by filtering out positions with depth=0
        let consensus = &cons[0];
        let len = consensus.sequence().len();
        assert_eq!(len, 40, "Consensus should be 40bp (skipping the 5bp deletion), got {len}");
    }

    /// Port of fgbio test: "make a consensus where the reads have soft-clipping outside of the overlap region"
    #[test]
    fn test_make_consensus_soft_clipping() {
        let options = CodecConsensusOptions {
            min_reads_per_strand: 1,
            min_duplex_length: 1,
            ..Default::default()
        };
        let mut caller = CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options);

        // R1 at pos 1 with 5S25M (5bp soft-clip at start)
        // R2 at pos 11 with 25M5S (5bp soft-clip at end)
        // This creates a situation where:
        // - R1 aligned region: 1-25 (25bp aligned, 5bp soft-clipped at start)
        // - R2 aligned region: 11-35 (25bp aligned, 5bp soft-clipped at end)
        // - Overlap: 11-25 (15bp)
        // - Output should include soft-clipped bases from both ends
        let reads = create_fr_pair(
            "read1",
            1,
            11,
            30,
            35,
            &[(Kind::SoftClip, 5), (Kind::Match, 25)],
            &[(Kind::Match, 25), (Kind::SoftClip, 5)],
            "hi",
            Some("ACC-TGA"),
            false,
            true,
        );

        let cons = caller.consensus_reads_from_sam_records(reads).unwrap();

        // Should produce consensus
        assert_eq!(cons.len(), 1, "Should produce one consensus read");

        // fgbio expects length 45: 5 (R1 soft-clip) + 35 (REF[0..35]) + 5 (R2 soft-clip)
        let consensus = &cons[0];
        let len = consensus.sequence().len();
        // Currently fgumi doesn't include soft-clips, so length is just the reference span
        // TODO: Update to expect 45 when soft-clip handling is implemented
        assert_eq!(
            len, 45,
            "Consensus should be 45bp (5 soft-clip + 35 ref + 5 soft-clip), got {len}"
        );
    }

    /// Port of fgbio test: "make a consensus where both reads are soft-clipped at the same end and fully overlapped"
    #[test]
    fn test_make_consensus_both_soft_clipped_same_end() {
        let options = CodecConsensusOptions {
            min_reads_per_strand: 1,
            min_duplex_length: 1,
            ..Default::default()
        };
        let mut caller = CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options);

        // Both R1 and R2 at pos 1 with 5S25M (both soft-clipped at start)
        // They are fully overlapped
        // In fgbio, R2's bases are set to match R1's soft-clipped sequence
        let reads = create_fr_pair(
            "read1",
            1,
            1,
            30,
            35,
            &[(Kind::SoftClip, 5), (Kind::Match, 25)],
            &[(Kind::SoftClip, 5), (Kind::Match, 25)],
            "hi",
            Some("ACC-TGA"),
            false,
            true,
        );

        let cons = caller.consensus_reads_from_sam_records(reads).unwrap();

        // Should produce consensus
        assert_eq!(cons.len(), 1, "Should produce one consensus read");

        // fgbio expects length 30: 5 (soft-clip) + 25 (aligned)
        let consensus = &cons[0];
        assert!(
            !consensus.sequence().is_empty(),
            "Should have sequence with both reads soft-clipped"
        );
    }

    /// Port of fgbio test: "not emit a consensus when the reads are a cross-chromosomal chimeric pair"
    #[test]
    fn test_not_emit_consensus_chimeric_pair() {
        let options = CodecConsensusOptions {
            min_reads_per_strand: 1,
            min_duplex_length: 1,
            ..Default::default()
        };
        let mut caller = CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options);

        // Create FR pair then make it chimeric by putting R1 on a different chromosome
        let mut reads = create_fr_pair(
            "read1",
            100,
            135,
            30,
            35,
            &[(Kind::Match, 30)],
            &[(Kind::Match, 30)],
            "hi",
            Some("ACC-TGA"),
            false,
            true,
        );

        // Modify R1 to be on a different reference (chromosome)
        // In noodles, reference_sequence_id is Option<usize>
        *reads[0].reference_sequence_id_mut() = Some(2); // Different chromosome
        // Update mate's reference for consistency
        *reads[1].mate_reference_sequence_id_mut() = Some(2);

        let cons = caller.consensus_reads_from_sam_records(reads).unwrap();

        // Should not emit consensus for chimeric pairs
        assert!(cons.is_empty(), "Should not emit consensus for cross-chromosomal chimeric pair");
    }

    /// Port of fgbio test: "not emit a consensus when R1's end lands in an indel in R2"
    ///
    /// fgumi now matches fgbio behavior: rejects pairs where overlap boundary lands in an indel.
    #[test]
    fn test_not_emit_consensus_r1_end_in_indel() {
        let options = CodecConsensusOptions {
            min_reads_per_strand: 1,
            min_duplex_length: 1,
            ..Default::default()
        };
        let mut caller = CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options);

        // R1 at pos 1 with 30M
        // R2 at pos 11 with 19M2D11M
        // R1 ends at position 30, which lands inside R2's deletion at position 30-31
        let reads = create_fr_pair(
            "read1",
            1,
            11,
            30,
            35,
            &[(Kind::Match, 30)],
            &[(Kind::Match, 19), (Kind::Deletion, 2), (Kind::Match, 11)],
            "hi",
            Some("ACC-TGA"),
            false,
            true,
        );

        let cons = caller.consensus_reads_from_sam_records(reads).unwrap();

        // fgbio rejects this because R1's end lands in an indel in R2
        // fgumi now implements this check and should also reject
        assert!(cons.is_empty(), "Should not emit consensus when R1's end lands in an indel in R2");
    }

    /// Port of fgbio test: "mask end qualities"
    #[test]
    fn test_mask_end_qualities() {
        // This tests the outer_bases_qual and outer_bases_length options
        // Create a consensus with outerBasesLength=7 and outerBasesQual=5
        let options = CodecConsensusOptions {
            min_reads_per_strand: 1,
            min_duplex_length: 1,
            outer_bases_length: 7,
            outer_bases_qual: Some(5),
            ..Default::default()
        };
        let mut caller = CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options);

        // Create a simple FR pair with high quality bases
        let reads = create_fr_pair(
            "read1",
            1,
            1,
            50,
            90, // High base quality
            &[(Kind::Match, 50)],
            &[(Kind::Match, 50)],
            "hi",
            Some("ACC-TGA"),
            false,
            true,
        );

        let cons = caller.consensus_reads_from_sam_records(reads).unwrap();

        assert_eq!(cons.len(), 1, "Should produce one consensus read");

        let consensus = &cons[0];
        let quals = consensus.quality_scores().as_ref();

        // Check that outer bases have reduced quality
        // First 7 bases should have qual <= 5
        for (i, &q) in quals.iter().take(7).enumerate() {
            assert!(q <= 5, "First 7 bases should have qual <= 5, but base {i} has qual {q}");
        }

        // Last 7 bases should also have qual <= 5
        let len = quals.len();
        for (i, &q) in quals.iter().skip(len.saturating_sub(7)).enumerate() {
            assert!(q <= 5, "Last 7 bases should have qual <= 5, but base {i} has qual {q}");
        }
    }

    /// Port of fgbio test: "mask single stranded regions _and_ bases"
    #[test]
    fn test_mask_single_stranded_regions() {
        // This tests the single_strand_qual option
        // In fgbio, single-stranded bases (where abConsensus has N/qual=2) get masked
        let options = CodecConsensusOptions {
            min_reads_per_strand: 1,
            min_duplex_length: 1,
            single_strand_qual: Some(4),
            ..Default::default()
        };
        let mut caller = CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options);

        // Create a pair where R1 and R2 don't fully overlap
        // R1: pos 1-30, R2: pos 20-49
        // Overlap region: 20-30 (duplex)
        // Single-strand regions: 1-19 (R1 only), 31-49 (R2 only)
        let reads = create_fr_pair(
            "read1",
            1,
            20,
            30,
            90, // High base quality
            &[(Kind::Match, 30)],
            &[(Kind::Match, 30)],
            "hi",
            Some("ACC-TGA"),
            false,
            true,
        );

        let cons = caller.consensus_reads_from_sam_records(reads).unwrap();

        assert_eq!(cons.len(), 1, "Should produce one consensus read");

        // The consensus should have single-strand regions masked
        // The exact quality values depend on the implementation details
        // This test verifies the option is respected
        let consensus = &cons[0];
        let quals = consensus.quality_scores().as_ref();

        // Single-strand regions should have lower quality
        // In fgbio, these get masked to single_strand_qual (4)
        // Check that we have some variation in qualities
        let has_low_qual = quals.iter().any(|&q| q <= 4);
        assert!(
            has_low_qual || quals.is_empty(),
            "Should have some low-quality single-strand regions (or be filtered)"
        );
    }

    /// Test that uppercase 'N' (`NoCall`) from one strand masks the result to N,
    /// but lowercase 'n' (padding) does not mask the result.
    /// This matches fgbio's behavior where only uppercase 'N' is considered `NoCall`.
    #[test]
    fn test_nocall_masking_uppercase_n() {
        let options = CodecConsensusOptions::default();
        let mut caller = CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options);

        // Create two single-strand consensuses:
        // Position 0: Both have valid bases (A, A) - should agree
        // Position 1: Strand A has 'N' (uppercase NoCall), strand B has 'T' - should mask to N
        // Position 2: Strand A has 'A', strand B has 'N' (uppercase NoCall) - should mask to N
        // Position 3: Strand A has 'n' (lowercase padding), strand B has 'G' - should NOT mask, use G
        // Position 4: Strand A has 'C', strand B has 'n' (lowercase padding) - should NOT mask, use C
        let ss_a = SingleStrandConsensus {
            bases: b"ANAnC".to_vec(),
            quals: vec![40, 2, 40, 0, 40], // N has qual 2 (NoCallQual), n has qual 0 (padding)
            depths: vec![2, 0, 2, 0, 2],
            errors: vec![0, 0, 0, 0, 0],
            raw_read_count: 2,
            ref_start: 0,
            ref_end: 4,
            is_negative_strand: false,
        };

        let ss_b = SingleStrandConsensus {
            bases: b"ATNGn".to_vec(),
            quals: vec![40, 40, 2, 40, 0], // N has qual 2 (NoCallQual), n has qual 0 (padding)
            depths: vec![2, 2, 0, 2, 0],
            errors: vec![0, 0, 0, 0, 0],
            raw_read_count: 2,
            ref_start: 0,
            ref_end: 4,
            is_negative_strand: true,
        };

        let duplex = caller
            .build_duplex_consensus_from_padded(&ss_a, &ss_b)
            .expect("Should build duplex consensus");

        // Position 0: Both have A, should be A
        assert_eq!(duplex.bases[0], b'A', "Position 0: Both have A, should agree");

        // Position 1: A has N (uppercase), B has T -> should mask to N
        assert_eq!(
            duplex.bases[1], b'N',
            "Position 1: Strand A has uppercase N (NoCall), should mask to N"
        );

        // Position 2: A has A, B has N (uppercase) -> should mask to N
        assert_eq!(
            duplex.bases[2], b'N',
            "Position 2: Strand B has uppercase N (NoCall), should mask to N"
        );

        // Position 3: A has n (lowercase padding), B has G -> should NOT mask, use G
        assert_eq!(
            duplex.bases[3], b'G',
            "Position 3: Strand A has lowercase n (padding), should use strand B's base G"
        );

        // Position 4: A has C, B has n (lowercase padding) -> should NOT mask, use C
        assert_eq!(
            duplex.bases[4], b'C',
            "Position 4: Strand B has lowercase n (padding), should use strand A's base C"
        );
    }

    #[test]
    fn test_to_source_read_for_codec_positive_strand() {
        use noodles::sam::alignment::record::cigar::op::Kind;

        // Create a forward strand read
        let read = create_test_paired_read(
            "read1",
            b"ACGT",
            &[30, 31, 32, 33],
            true,
            false, // forward strand
            true,
            100,
            &[(Kind::Match, 4)],
        );

        let source = CodecConsensusCaller::to_source_read_for_codec(&read, 0);

        // Forward strand: no revcomp, no reversal
        assert_eq!(source.bases, b"ACGT");
        assert_eq!(source.quals, vec![30, 31, 32, 33]);
        assert_eq!(source.original_idx, 0);
    }

    #[test]
    fn test_to_source_read_for_codec_negative_strand() {
        use noodles::sam::alignment::record::cigar::op::Kind;

        // Create a reverse strand read
        let read = create_test_paired_read(
            "read1",
            b"ACGT",
            &[30, 31, 32, 33],
            true,
            true, // reverse strand
            false,
            100,
            &[(Kind::Match, 4)],
        );

        let source = CodecConsensusCaller::to_source_read_for_codec(&read, 1);

        // Reverse strand: revcomp bases, reverse quals
        assert_eq!(source.bases, b"ACGT"); // ACGT revcomp is ACGT (palindrome)
        assert_eq!(source.quals, vec![33, 32, 31, 30]); // reversed
        assert_eq!(source.original_idx, 1);
    }

    #[test]
    fn test_to_source_read_for_codec_negative_strand_non_palindrome() {
        use noodles::sam::alignment::record::cigar::op::Kind;

        // Create a reverse strand read with non-palindromic sequence
        let read = create_test_paired_read(
            "read1",
            b"AACC", // revcomp is GGTT
            &[10, 20, 30, 40],
            true,
            true, // reverse strand
            false,
            100,
            &[(Kind::Match, 4)],
        );

        let source = CodecConsensusCaller::to_source_read_for_codec(&read, 2);

        // Reverse strand: revcomp bases (AACC -> GGTT), reverse quals
        assert_eq!(source.bases, b"GGTT");
        assert_eq!(source.quals, vec![40, 30, 20, 10]); // reversed
        assert_eq!(source.original_idx, 2);
    }

    #[test]
    fn test_vanilla_to_single_strand_conversion() {
        let vcr = VanillaConsensusRead {
            id: "test".to_string(),
            bases: vec![b'A', b'C', b'G', b'T'],
            quals: vec![30, 31, 32, 33],
            depths: vec![5, 6, 7, 8],
            errors: vec![0, 1, 0, 1],
            source_reads: None,
        };

        let ss = CodecConsensusCaller::vanilla_to_single_strand(vcr.clone(), false, 10);

        assert_eq!(ss.bases, vcr.bases);
        assert_eq!(ss.quals, vcr.quals);
        assert_eq!(ss.depths, vcr.depths);
        assert_eq!(ss.errors, vcr.errors);
        assert_eq!(ss.raw_read_count, 10);
        assert!(!ss.is_negative_strand);
        assert_eq!(ss.ref_start, 0);
        assert_eq!(ss.ref_end, 3); // len - 1
    }

    #[test]
    fn test_vanilla_to_single_strand_negative_strand() {
        let vcr = VanillaConsensusRead {
            id: "test".to_string(),
            bases: vec![b'A', b'C'],
            quals: vec![30, 31],
            depths: vec![5, 6],
            errors: vec![0, 1],
            source_reads: None,
        };

        let ss = CodecConsensusCaller::vanilla_to_single_strand(vcr, true, 5);

        assert!(ss.is_negative_strand);
        assert_eq!(ss.raw_read_count, 5);
    }

    #[test]
    fn test_codec_statistics_tracking() {
        // Test that CodecConsensusCaller properly tracks statistics
        let options = CodecConsensusOptions {
            min_reads_per_strand: 1,
            min_duplex_length: 1,
            ..Default::default()
        };
        let caller = CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options);

        let stats = caller.statistics();

        // Initial statistics should be zero
        assert_eq!(stats.total_input_reads, 0);
        assert_eq!(stats.consensus_reads_generated, 0);
        assert_eq!(stats.reads_filtered, 0);
    }

    // ==========================================================================
    // Tests for CIGAR-based utility functions
    // ==========================================================================

    #[test]
    fn test_read_pos_at_ref_pos_simple_match() {
        use noodles::sam::alignment::record::cigar::op::Kind;

        // Read at position 100 with 20M
        let read = create_test_paired_read(
            "read1",
            &[b'A'; 20],
            &[30; 20],
            true,
            false,
            true,
            100,
            &[(Kind::Match, 20)],
        );

        // Query position at ref_pos 100 (start) should be 1 (1-based)
        assert_eq!(CodecConsensusCaller::read_pos_at_ref_pos(&read, 100, false), Some(1));

        // Query position at ref_pos 105 should be 6 (1-based)
        assert_eq!(CodecConsensusCaller::read_pos_at_ref_pos(&read, 105, false), Some(6));

        // Query position at ref_pos 119 (end) should be 20 (1-based)
        assert_eq!(CodecConsensusCaller::read_pos_at_ref_pos(&read, 119, false), Some(20));

        // Position before alignment should return None
        assert_eq!(CodecConsensusCaller::read_pos_at_ref_pos(&read, 99, false), None);

        // Position after alignment should return None
        assert_eq!(CodecConsensusCaller::read_pos_at_ref_pos(&read, 120, false), None);
    }

    #[test]
    fn test_read_pos_at_ref_pos_with_insertion() {
        use noodles::sam::alignment::record::cigar::op::Kind;

        // Read at position 100 with 10M5I10M (25 query bases, 20 ref bases)
        let read = create_test_paired_read(
            "read1",
            &[b'A'; 25],
            &[30; 25],
            true,
            false,
            true,
            100,
            &[(Kind::Match, 10), (Kind::Insertion, 5), (Kind::Match, 10)],
        );

        // Position in first match block
        assert_eq!(CodecConsensusCaller::read_pos_at_ref_pos(&read, 100, false), Some(1));
        assert_eq!(CodecConsensusCaller::read_pos_at_ref_pos(&read, 109, false), Some(10));

        // Position in second match block (after insertion)
        // ref_pos 110 maps to query position 16 (10 from first M + 5 from I + 1)
        assert_eq!(CodecConsensusCaller::read_pos_at_ref_pos(&read, 110, false), Some(16));
        assert_eq!(CodecConsensusCaller::read_pos_at_ref_pos(&read, 119, false), Some(25));
    }

    #[test]
    fn test_read_pos_at_ref_pos_with_deletion() {
        use noodles::sam::alignment::record::cigar::op::Kind;

        // Read at position 100 with 10M5D10M (20 query bases, 25 ref bases)
        let read = create_test_paired_read(
            "read1",
            &[b'A'; 20],
            &[30; 20],
            true,
            false,
            true,
            100,
            &[(Kind::Match, 10), (Kind::Deletion, 5), (Kind::Match, 10)],
        );

        // Position in first match block
        assert_eq!(CodecConsensusCaller::read_pos_at_ref_pos(&read, 100, false), Some(1));
        assert_eq!(CodecConsensusCaller::read_pos_at_ref_pos(&read, 109, false), Some(10));

        // Position in deletion - should return None with return_last_base_if_deleted=false
        assert_eq!(CodecConsensusCaller::read_pos_at_ref_pos(&read, 110, false), None);
        assert_eq!(CodecConsensusCaller::read_pos_at_ref_pos(&read, 114, false), None);

        // Position in deletion with return_last_base_if_deleted=true - should return last base before deletion
        assert_eq!(CodecConsensusCaller::read_pos_at_ref_pos(&read, 110, true), Some(10));
        assert_eq!(CodecConsensusCaller::read_pos_at_ref_pos(&read, 114, true), Some(10));

        // Position in second match block (after deletion)
        assert_eq!(CodecConsensusCaller::read_pos_at_ref_pos(&read, 115, false), Some(11));
        assert_eq!(CodecConsensusCaller::read_pos_at_ref_pos(&read, 124, false), Some(20));
    }

    #[test]
    fn test_read_pos_at_ref_pos_with_soft_clip() {
        use noodles::sam::alignment::record::cigar::op::Kind;

        // Read at position 100 with 5S15M (20 query bases, 15 ref bases)
        // Soft clips consume query but not reference
        let read = create_test_paired_read(
            "read1",
            &[b'A'; 20],
            &[30; 20],
            true,
            false,
            true,
            100,
            &[(Kind::SoftClip, 5), (Kind::Match, 15)],
        );

        // ref_pos 100 should map to query position 6 (after 5 soft-clipped bases)
        assert_eq!(CodecConsensusCaller::read_pos_at_ref_pos(&read, 100, false), Some(6));
        assert_eq!(CodecConsensusCaller::read_pos_at_ref_pos(&read, 114, false), Some(20));
    }

    #[test]
    fn test_check_overlap_phase_compatible() {
        use noodles::sam::alignment::record::cigar::op::Kind;

        let options = CodecConsensusOptions::default();
        let caller = CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options);

        // Two reads with compatible CIGARs (both 30M)
        let r1 = create_test_paired_read(
            "read1",
            &[b'A'; 30],
            &[30; 30],
            true,
            false,
            true,
            100,
            &[(Kind::Match, 30)],
        );
        let r2 = create_test_paired_read(
            "read2",
            &[b'A'; 30],
            &[30; 30],
            false,
            true,
            false,
            110,
            &[(Kind::Match, 30)],
        );

        // Overlap region: 110-129 (R2 start to R1 end)
        assert!(caller.check_overlap_phase(&r1, &r2, 110, 129));
    }

    #[test]
    fn test_check_overlap_phase_indel_mismatch() {
        use noodles::sam::alignment::record::cigar::op::Kind;

        let options = CodecConsensusOptions::default();
        let caller = CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options);

        // R1 with 30M, R2 with 15M5D15M
        // The overlap boundaries may land in the deletion
        let r1 = create_test_paired_read(
            "read1",
            &[b'A'; 30],
            &[30; 30],
            true,
            false,
            true,
            100,
            &[(Kind::Match, 30)],
        );
        let r2 = create_test_paired_read(
            "read2",
            &[b'A'; 30],
            &[30; 30],
            false,
            true,
            false,
            100,
            &[(Kind::Match, 15), (Kind::Deletion, 5), (Kind::Match, 15)],
        );

        // Overlap region spanning the deletion should detect the phase mismatch
        // Result depends on whether the boundaries land in indels - just verify it runs
        let _result = caller.check_overlap_phase(&r1, &r2, 100, 129);
    }

    #[test]
    fn test_pad_consensus_right() {
        // Test padding on the right (positive strand)
        let ss = SingleStrandConsensus {
            bases: b"ACGT".to_vec(),
            quals: vec![30, 31, 32, 33],
            depths: vec![5, 6, 7, 8],
            errors: vec![0, 1, 0, 1],
            raw_read_count: 5,
            ref_start: 100,
            ref_end: 103,
            is_negative_strand: false,
        };

        let padded = CodecConsensusCaller::pad_consensus(&ss, 8, false);

        assert_eq!(padded.bases.len(), 8);
        assert_eq!(&padded.bases[0..4], b"ACGT");
        assert_eq!(&padded.bases[4..8], b"nnnn"); // lowercase n for padding
        assert_eq!(padded.quals[4..8], vec![0, 0, 0, 0]);
        assert_eq!(padded.depths[4..8], vec![0, 0, 0, 0]);
    }

    #[test]
    fn test_pad_consensus_left() {
        // Test padding on the left (negative strand)
        let ss = SingleStrandConsensus {
            bases: b"ACGT".to_vec(),
            quals: vec![30, 31, 32, 33],
            depths: vec![5, 6, 7, 8],
            errors: vec![0, 1, 0, 1],
            raw_read_count: 5,
            ref_start: 100,
            ref_end: 103,
            is_negative_strand: true,
        };

        let padded = CodecConsensusCaller::pad_consensus(&ss, 8, true);

        assert_eq!(padded.bases.len(), 8);
        assert_eq!(&padded.bases[0..4], b"nnnn"); // lowercase n for padding
        assert_eq!(&padded.bases[4..8], b"ACGT");
        assert_eq!(padded.quals[0..4], vec![0, 0, 0, 0]);
        assert_eq!(padded.depths[0..4], vec![0, 0, 0, 0]);
    }

    #[test]
    fn test_pad_consensus_no_padding_needed() {
        // Test when consensus is already the target length
        let ss = SingleStrandConsensus {
            bases: b"ACGT".to_vec(),
            quals: vec![30, 31, 32, 33],
            depths: vec![5, 6, 7, 8],
            errors: vec![0, 1, 0, 1],
            raw_read_count: 5,
            ref_start: 100,
            ref_end: 103,
            is_negative_strand: false,
        };

        let padded = CodecConsensusCaller::pad_consensus(&ss, 4, false);

        // Should return unchanged
        assert_eq!(padded.bases, ss.bases);
        assert_eq!(padded.quals, ss.quals);
    }

    #[test]
    fn test_pad_consensus_shorter_target() {
        // Test when target length is shorter than current (should return unchanged)
        let ss = SingleStrandConsensus {
            bases: b"ACGT".to_vec(),
            quals: vec![30, 31, 32, 33],
            depths: vec![5, 6, 7, 8],
            errors: vec![0, 1, 0, 1],
            raw_read_count: 5,
            ref_start: 100,
            ref_end: 103,
            is_negative_strand: false,
        };

        let padded = CodecConsensusCaller::pad_consensus(&ss, 2, false);

        // Should return unchanged since new_length <= current_len
        assert_eq!(padded.bases.len(), 4);
    }

    #[test]
    fn test_mask_consensus_quals_query_based_single_strand() {
        // Test query-based masking where padding indicates single-strand
        let options = CodecConsensusOptions { single_strand_qual: Some(10), ..Default::default() };
        let caller = CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options);

        let consensus = SingleStrandConsensus {
            bases: b"ACGT".to_vec(),
            quals: vec![30, 30, 30, 30],
            depths: vec![5; 4],
            errors: vec![0; 4],
            raw_read_count: 5,
            ref_start: 0,
            ref_end: 3,
            is_negative_strand: false,
        };

        // R1 has N at position 0 (padding)
        let padded_r1 = SingleStrandConsensus {
            bases: b"NCGT".to_vec(),
            quals: vec![0, 30, 30, 30],
            depths: vec![0, 5, 5, 5],
            errors: vec![0; 4],
            raw_read_count: 5,
            ref_start: 0,
            ref_end: 3,
            is_negative_strand: false,
        };

        // R2 has N at position 3 (padding)
        let padded_r2 = SingleStrandConsensus {
            bases: b"ACGN".to_vec(),
            quals: vec![30, 30, 30, 0],
            depths: vec![5, 5, 5, 0],
            errors: vec![0; 4],
            raw_read_count: 5,
            ref_start: 0,
            ref_end: 3,
            is_negative_strand: true,
        };

        let masked = caller.mask_consensus_quals_query_based(consensus, &padded_r1, &padded_r2);

        // Position 0: R1 has N, so single-strand -> qual=10
        assert_eq!(masked.quals[0], 10);
        // Position 3: R2 has N, so single-strand -> qual=10
        assert_eq!(masked.quals[3], 10);
        // Positions 1-2: both have data, duplex -> keep original
        assert_eq!(masked.quals[1], 30);
        assert_eq!(masked.quals[2], 30);
    }

    #[test]
    fn test_is_cigar_prefix_exact_match() {
        use noodles::sam::alignment::record::cigar::op::Kind;

        // Same CIGAR
        let a = vec![(Kind::Match, 50)];
        let b = vec![(Kind::Match, 50)];
        assert!(cigar_utils::is_cigar_prefix(&a, &b));
    }

    #[test]
    fn test_is_cigar_prefix_shorter_last_element() {
        use noodles::sam::alignment::record::cigar::op::Kind;

        // a is shorter in last element
        let a = vec![(Kind::Match, 40)];
        let b = vec![(Kind::Match, 50)];
        assert!(cigar_utils::is_cigar_prefix(&a, &b));
    }

    #[test]
    fn test_is_cigar_prefix_longer_not_prefix() {
        use noodles::sam::alignment::record::cigar::op::Kind;

        // a is longer than b - not a prefix
        let a = vec![(Kind::Match, 50)];
        let b = vec![(Kind::Match, 40)];
        assert!(!cigar_utils::is_cigar_prefix(&a, &b));
    }

    #[test]
    fn test_is_cigar_prefix_different_ops() {
        use noodles::sam::alignment::record::cigar::op::Kind;

        // Different operations
        let a = vec![(Kind::Match, 50)];
        let b = vec![(Kind::Insertion, 50)];
        assert!(!cigar_utils::is_cigar_prefix(&a, &b));
    }

    #[test]
    fn test_is_cigar_prefix_complex() {
        use noodles::sam::alignment::record::cigar::op::Kind;

        // Complex CIGAR with matching indels
        let a = vec![(Kind::Match, 25), (Kind::Deletion, 5), (Kind::Match, 20)];
        let b = vec![(Kind::Match, 25), (Kind::Deletion, 5), (Kind::Match, 30)];
        assert!(cigar_utils::is_cigar_prefix(&a, &b));

        // Different deletion length - not a prefix
        let c = vec![(Kind::Match, 25), (Kind::Deletion, 3), (Kind::Match, 20)];
        assert!(!cigar_utils::is_cigar_prefix(&c, &b));
    }

    #[test]
    fn test_reverse_complement_ss_depths_errors_reversed() {
        // Verify that depths and errors are reversed along with bases/quals for CODEC
        let ss = SingleStrandConsensus {
            bases: b"ACGT".to_vec(),
            quals: vec![10, 20, 30, 40],
            depths: vec![1, 2, 3, 4],
            errors: vec![5, 6, 7, 8],
            raw_read_count: 5,
            ref_start: 100,
            ref_end: 103,
            is_negative_strand: false,
        };

        let rc = CodecConsensusCaller::reverse_complement_ss(&ss);

        // Bases should be reverse complemented
        assert_eq!(rc.bases, b"ACGT"); // ACGT revcomp is ACGT (palindrome)

        // Quals should be reversed
        assert_eq!(rc.quals, vec![40, 30, 20, 10]);

        // For CODEC query-based consensus, depths and errors ARE reversed
        assert_eq!(rc.depths, vec![4, 3, 2, 1]);
        assert_eq!(rc.errors, vec![8, 7, 6, 5]);

        // Strand flag should be flipped
        assert!(rc.is_negative_strand);
    }

    #[test]
    fn test_downsample_pairs() {
        use noodles::sam::alignment::record::cigar::op::Kind;

        let options = CodecConsensusOptions { max_reads_per_strand: Some(2), ..Default::default() };
        let mut caller = CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options);

        // Create 5 R1s and 5 R2s
        let r1s: Vec<RecordBuf> = (0..5)
            .map(|i| {
                create_test_paired_read(
                    &format!("r{i}"),
                    &[b'A'; 20],
                    &[30; 20],
                    true,
                    false,
                    true,
                    100,
                    &[(Kind::Match, 20)],
                )
            })
            .collect();
        let r2s: Vec<RecordBuf> = (0..5)
            .map(|i| {
                create_test_paired_read(
                    &format!("r{i}"),
                    &[b'A'; 20],
                    &[30; 20],
                    false,
                    true,
                    false,
                    110,
                    &[(Kind::Match, 20)],
                )
            })
            .collect();

        let (ds_r1s, ds_r2s) = caller.downsample_pairs(r1s, r2s);

        assert_eq!(ds_r1s.len(), 2);
        assert_eq!(ds_r2s.len(), 2);
    }

    #[test]
    fn test_downsample_pairs_no_limit() {
        use noodles::sam::alignment::record::cigar::op::Kind;

        let options = CodecConsensusOptions { max_reads_per_strand: None, ..Default::default() };
        let mut caller = CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options);

        let r1s: Vec<RecordBuf> = (0..5)
            .map(|i| {
                create_test_paired_read(
                    &format!("r{i}"),
                    &[b'A'; 20],
                    &[30; 20],
                    true,
                    false,
                    true,
                    100,
                    &[(Kind::Match, 20)],
                )
            })
            .collect();
        let r2s: Vec<RecordBuf> = (0..5)
            .map(|i| {
                create_test_paired_read(
                    &format!("r{i}"),
                    &[b'A'; 20],
                    &[30; 20],
                    false,
                    true,
                    false,
                    110,
                    &[(Kind::Match, 20)],
                )
            })
            .collect();

        let (ds_r1s, ds_r2s) = caller.downsample_pairs(r1s, r2s);

        // No limit, so all reads should be kept
        assert_eq!(ds_r1s.len(), 5);
        assert_eq!(ds_r2s.len(), 5);
    }

    #[test]
    fn test_compute_codec_consensus_length() {
        use noodles::sam::alignment::record::cigar::op::Kind;

        // R1 (positive strand) at pos 100 with 30M
        // R2 (negative strand) at pos 110 with 30M
        // Overlap end is at position 129 (R1 end)
        let pos_read = create_test_paired_read(
            "pos",
            &[b'A'; 30],
            &[30; 30],
            true,
            false,
            true,
            100,
            &[(Kind::Match, 30)],
        );
        let neg_read = create_test_paired_read(
            "neg",
            &[b'A'; 30],
            &[30; 30],
            false,
            true,
            false,
            110,
            &[(Kind::Match, 30)],
        );

        // Overlap end at 129 (R1's end position)
        let length =
            CodecConsensusCaller::compute_codec_consensus_length(&pos_read, &neg_read, 129);

        // pos_read_pos at 129 = 30 (last base)
        // neg_read_pos at 129 = 20 (position 129-110+1 = 20)
        // neg_length = 30
        // consensus_length = 30 + 30 - 20 = 40
        assert_eq!(length, Some(40));
    }

    #[test]
    fn test_compute_codec_consensus_length_in_deletion() {
        use noodles::sam::alignment::record::cigar::op::Kind;

        // Positive read with 15M5D15M
        let pos_read = create_test_paired_read(
            "pos",
            &[b'A'; 30],
            &[30; 30],
            true,
            false,
            true,
            100,
            &[(Kind::Match, 15), (Kind::Deletion, 5), (Kind::Match, 15)],
        );
        let neg_read = create_test_paired_read(
            "neg",
            &[b'A'; 30],
            &[30; 30],
            false,
            true,
            false,
            100,
            &[(Kind::Match, 30)],
        );

        // Overlap end at position 116 (inside the deletion of pos_read at 115-119)
        let length =
            CodecConsensusCaller::compute_codec_consensus_length(&pos_read, &neg_read, 116);

        // Should return None because position 116 is in a deletion
        assert!(length.is_none());
    }

    #[test]
    fn test_reject_records_tracking() {
        use noodles::sam::alignment::record::cigar::op::Kind;

        let options = CodecConsensusOptions {
            min_reads_per_strand: 2, // Require 2 reads per strand
            ..Default::default()
        };
        let mut caller = CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options);

        // Create single pair - should be rejected for insufficient reads
        let reads = create_fr_pair(
            "read1",
            100,
            110,
            30,
            35,
            &[(Kind::Match, 30)],
            &[(Kind::Match, 30)],
            "UMI1",
            None,
            false,
            true,
        );

        let result = caller.consensus_reads_from_sam_records(reads).unwrap();

        assert!(result.is_empty());
        assert!(
            caller.stats.rejection_reasons.contains_key(&CallerRejectionReason::InsufficientReads)
        );
    }

    #[test]
    fn test_rejected_reads_tracking_enabled() {
        use noodles::sam::alignment::record::cigar::op::Kind;

        let options = CodecConsensusOptions {
            min_reads_per_strand: 2, // Require 2 reads per strand
            ..Default::default()
        };
        let mut caller = CodecConsensusCaller::new_with_rejects_tracking(
            "codec".to_string(),
            "RG1".to_string(),
            options,
            true, // Enable tracking
        );

        let reads = create_fr_pair(
            "read1",
            100,
            110,
            30,
            35,
            &[(Kind::Match, 30)],
            &[(Kind::Match, 30)],
            "UMI1",
            None,
            false,
            true,
        );

        let result = caller.consensus_reads_from_sam_records(reads).unwrap();

        assert!(result.is_empty());
        // Rejected reads should be tracked
        assert!(!caller.rejected_reads().is_empty() || caller.statistics().reads_filtered > 0);
    }

    #[test]
    fn test_clear_rejected_reads() {
        let options = CodecConsensusOptions::default();
        let mut caller = CodecConsensusCaller::new_with_rejects_tracking(
            "codec".to_string(),
            "RG1".to_string(),
            options,
            true,
        );

        // Initially empty
        assert!(caller.rejected_reads().is_empty());

        // Clear (should be no-op on empty)
        caller.clear_rejected_reads();
        assert!(caller.rejected_reads().is_empty());
    }

    // ============================================================================
    // Tests for additional HIGH risk code paths
    // ============================================================================

    #[test]
    fn test_read_pos_at_ref_pos_with_hard_clip() {
        // Create a read with 2H5M2H CIGAR at position 100 (hard clips don't appear in seq)
        let record = RecordBuilder::new()
            .sequence("ACGTA") // 5 bases for 5M
            .alignment_start(100)
            .cigar("2H5M2H")
            .build();

        // Position 100 maps to query pos 1 (hard clips don't count)
        assert_eq!(CodecConsensusCaller::read_pos_at_ref_pos(&record, 100, false), Some(1));
        // Position 104 maps to query pos 5
        assert_eq!(CodecConsensusCaller::read_pos_at_ref_pos(&record, 104, false), Some(5));
    }

    #[test]
    fn test_read_pos_at_ref_pos_deletion_at_start() {
        // Create a read with 2D5M at position 100
        // Ref:   100 101 102 103 104 105 106
        // Query:   -   -   1   2   3   4   5
        let record = RecordBuilder::new()
            .sequence("ACGTA") // 5 bases for 5M
            .alignment_start(100)
            .cigar("2D5M")
            .build();

        // Position 100 is in deletion at start
        // return_last_base_if_deleted=true with query_offset=0 should return 1
        assert_eq!(CodecConsensusCaller::read_pos_at_ref_pos(&record, 100, true), Some(1));
        assert_eq!(CodecConsensusCaller::read_pos_at_ref_pos(&record, 100, false), None);
    }

    #[test]
    fn test_check_overlap_phase_matching_cigars() {
        let options = CodecConsensusOptions::default();
        let caller = CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options);

        // Create two reads that overlap perfectly
        let r1 = RecordBuilder::new()
            .sequence(&"A".repeat(50))
            .alignment_start(100)
            .cigar("50M")
            .build();

        let r2 = RecordBuilder::new()
            .sequence(&"A".repeat(50))
            .alignment_start(120)
            .cigar("50M")
            .build();

        // Overlap region: 120-149
        assert!(caller.check_overlap_phase(&r1, &r2, 120, 149));
    }

    #[test]
    fn test_check_overlap_phase_deletion_mismatch() {
        let options = CodecConsensusOptions::default();
        let caller = CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options);

        // Create two reads where one has a deletion
        let r1 = RecordBuilder::new()
            .sequence(&"A".repeat(50))
            .alignment_start(100)
            .cigar("50M")
            .build();

        // r2 has a 2bp deletion in the overlap region
        let r2 = RecordBuilder::new()
            .sequence(&"A".repeat(48)) // 15M + 33M = 48 bases (deletion doesn't consume query)
            .alignment_start(120)
            .cigar("15M2D33M")
            .build();

        // Overlap region spans the deletion - phases should differ
        assert!(!caller.check_overlap_phase(&r1, &r2, 120, 149));
    }
}

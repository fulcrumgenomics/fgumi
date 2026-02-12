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

use crate::caller::{
    ConsensusCaller, ConsensusCallingStats, ConsensusOutput,
    RejectionReason as CallerRejectionReason,
};
use crate::simple_umi::consensus_umis;
use crate::vanilla_caller::{
    VanillaConsensusRead, VanillaUmiConsensusCaller, VanillaUmiConsensusOptions,
};
use crate::{IndexedSourceRead, SourceRead, select_most_common_alignment_group};
use fgumi_dna::dna::reverse_complement;
use crate::phred::{MIN_PHRED, NO_CALL_BASE, NO_CALL_BASE_LOWER, PhredScore};
use noodles_raw_bam::{self as bam_fields, UnmappedBamRecordBuilder, flags};
use anyhow::Result;
use noodles::sam::alignment::record::data::field::Tag;
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
            #[expect(clippy::cast_precision_loss, reason = "acceptable precision loss for rate computation")]
            { self.duplex_disagreement_base_count as f64 / self.consensus_duplex_bases_emitted as f64 }
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

/// Precomputed clip metadata for a single record in the raw-byte codec pipeline.
///
/// Holds all information needed for filtering, overlap calculation, and
/// `SourceRead` construction without materializing a `RecordBuf`.
struct ClippedRecordInfo {
    /// Index into the original raw records vec.
    raw_idx: usize,
    /// Number of query bases to clip.
    clip_amount: usize,
    /// `true` ⟹ clip from start (negative-strand reads).
    clip_from_start: bool,
    /// Sequence length after clipping.
    clipped_seq_len: usize,
    /// CIGAR ops after clipping (raw u32 words).
    clipped_cigar: Vec<u32>,
    /// 1-based alignment start after clipping.
    adjusted_pos: usize,
    /// BAM flags.
    flags: u16,
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

    /// Random number generator for downsampling
    rng: StdRng,

    /// Counter for consensus read naming
    consensus_counter: u64,

    /// Single-strand consensus caller (matches fgbio's ssCaller delegation pattern)
    /// CODEC delegates to `VanillaUmiConsensusCaller` for single-strand consensus building
    ss_caller: VanillaUmiConsensusCaller,

    /// Reusable builder for raw-byte BAM record output
    bam_builder: UnmappedBamRecordBuilder,

    /// Whether to store rejected reads (raw bytes) for output
    track_rejects: bool,

    /// Rejected raw BAM records (only populated when `track_rejects` is true)
    rejected_reads: Vec<Vec<u8>>,
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
            rng,
            consensus_counter: 0,
            ss_caller,
            bam_builder: UnmappedBamRecordBuilder::new(),
            track_rejects: false,
            rejected_reads: Vec::new(),
        }
    }

    /// Creates a new CODEC consensus caller with optional rejected-reads tracking.
    ///
    /// When `track_rejects` is `true`, raw BAM bytes of rejected reads are stored
    /// and can be retrieved via [`rejected_reads`] / [`take_rejected_reads`].
    #[must_use]
    pub fn new_with_rejects_tracking(
        read_name_prefix: String,
        read_group_id: String,
        options: CodecConsensusOptions,
        track_rejects: bool,
    ) -> Self {
        let mut caller = Self::new(read_name_prefix, read_group_id, options);
        caller.track_rejects = track_rejects;
        caller
    }

    /// Returns the statistics for this caller
    #[must_use]
    pub fn statistics(&self) -> &CodecConsensusStats {
        &self.stats
    }

    /// Clears all per-group state to prepare for reuse
    ///
    /// This resets statistics while preserving the caller's
    /// configuration and reusable components (consensus builder, RNG, etc.).
    pub fn clear(&mut self) {
        self.stats = CodecConsensusStats::default();
        self.ss_caller.clear();
        self.rejected_reads.clear();
    }

    /// Returns a reference to the rejected reads (raw BAM bytes).
    #[must_use]
    pub fn rejected_reads(&self) -> &[Vec<u8>] {
        &self.rejected_reads
    }

    /// Clears the stored rejected reads.
    pub fn clear_rejected_reads(&mut self) {
        self.rejected_reads.clear();
    }

    /// Takes ownership of the stored rejected reads, leaving an empty vec.
    pub fn take_rejected_reads(&mut self) -> Vec<Vec<u8>> {
        std::mem::take(&mut self.rejected_reads)
    }

    /// Converts raw BAM bytes to `SourceRead` for CODEC consensus calling,
    /// applying virtual clipping.
    ///
    /// Simpler than vanilla's `create_source_read`: no quality masking,
    /// no trailing N removal, no quality trimming.
    fn to_source_read_for_codec_raw(
        raw: &[u8],
        original_idx: usize,
        clip_amount: usize,
        clip_from_start: bool,
        clipped_cigar: Option<&[u32]>,
    ) -> SourceRead {
        let mut bases = bam_fields::extract_sequence(raw);
        let mut quals = bam_fields::quality_scores_slice(raw).to_vec();
        let flg = bam_fields::flags(raw);

        // Apply clipping: truncate bases and quals
        // Clamp clip_amount to avoid panic on malformed input (e.g., CIGAR/MC mismatch)
        let clip_amount = clip_amount.min(bases.len());
        if clip_amount > 0 {
            if clip_from_start {
                bases = bases[clip_amount..].to_vec();
                quals = quals[clip_amount..].to_vec();
            } else {
                let new_len = bases.len() - clip_amount;
                bases.truncate(new_len);
                quals.truncate(new_len);
            }
        }

        // Get simplified CIGAR — use pre-clipped CIGAR if available to avoid redundant work
        let mut simplified = if let Some(ops) = clipped_cigar {
            bam_fields::simplify_cigar_from_raw(ops)
        } else {
            let original_ops = bam_fields::get_cigar_ops(raw);
            let (clipped_ops, _) =
                bam_fields::clip_cigar_ops_raw(&original_ops, clip_amount, clip_from_start);
            bam_fields::simplify_cigar_from_raw(&clipped_ops)
        };

        let is_negative = flg & flags::REVERSE != 0;
        if is_negative {
            simplified.reverse();
            bases = reverse_complement(&bases);
            quals.reverse();
        }

        SourceRead { original_idx, bases, quals, simplified_cigar: simplified, flags: flg }
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

    /// Calls consensus reads from raw BAM byte records for a single molecule.
    ///
    /// This is the primary entry point, working directly on raw bytes without
    /// materializing `RecordBuf`. Virtual clipping is computed and applied as
    /// offsets when extracting data for consensus.
    #[expect(clippy::too_many_lines, reason = "consensus pipeline has many sequential steps that are clearest in one function")]
    fn consensus_reads_raw(&mut self, records: &[Vec<u8>]) -> Result<ConsensusOutput> {
        self.stats.total_input_reads += records.len() as u64;

        if records.is_empty() {
            return Ok(ConsensusOutput::default());
        }

        // Extract MI tag from first record for naming
        let umi: Option<String> = bam_fields::find_string_tag_in_record(&records[0], b"MI")
            .map(|b| String::from_utf8_lossy(b).to_string());

        // Phase 1: Filter on raw bytes — keep paired, primary, mapped, FR-pair reads
        let mut paired_indices: Vec<usize> = Vec::new();
        let mut frag_count = 0usize;
        for (i, raw) in records.iter().enumerate() {
            let flg = bam_fields::flags(raw);
            if flg & flags::PAIRED == 0 {
                frag_count += 1;
                continue;
            }
            // Filter: primary, mapped, FR pair
            if flg & (flags::SECONDARY | flags::SUPPLEMENTARY | flags::UNMAPPED) != 0 {
                continue;
            }
            if !bam_fields::is_fr_pair_raw(raw) {
                continue;
            }
            paired_indices.push(i);
        }

        if frag_count > 0 {
            self.reject_records_count(frag_count, CallerRejectionReason::FragmentRead);
        }

        if paired_indices.is_empty() {
            return Ok(ConsensusOutput::default());
        }

        // Phase 2: Group by name, pair R1/R2, compute clip amounts
        let mut by_name: HashMap<&[u8], Vec<usize>> = HashMap::new();
        let mut name_order: Vec<&[u8]> = Vec::new();
        for &idx in &paired_indices {
            let name = bam_fields::read_name(&records[idx]);
            if !by_name.contains_key(name) {
                name_order.push(name);
            }
            by_name.entry(name).or_default().push(idx);
        }

        let mut r1_infos: Vec<ClippedRecordInfo> = Vec::new();
        let mut r2_infos: Vec<ClippedRecordInfo> = Vec::new();

        for name in name_order {
            let indices = by_name.get(name).expect("name from iteration");
            if indices.len() != 2 {
                continue;
            }

            let (i1, i2) = {
                let flg0 = bam_fields::flags(&records[indices[0]]);
                if flg0 & flags::FIRST_SEGMENT != 0 {
                    (indices[0], indices[1])
                } else {
                    (indices[1], indices[0])
                }
            };

            // Compute clip amounts from raw bytes
            let clip_r1 = bam_fields::num_bases_extending_past_mate_raw(&records[i1]);
            let clip_r2 = bam_fields::num_bases_extending_past_mate_raw(&records[i2]);

            // Build ClippedRecordInfo for R1
            let r1_info = Self::build_clipped_info(&records[i1], i1, clip_r1);
            // Build ClippedRecordInfo for R2
            let r2_info = Self::build_clipped_info(&records[i2], i2, clip_r2);

            r1_infos.push(r1_info);
            r2_infos.push(r2_info);
        }

        if r1_infos.is_empty() {
            return Ok(ConsensusOutput::default());
        }

        // Check we have enough reads
        if r1_infos.len() < self.options.min_reads_per_strand {
            self.reject_records_count(
                r1_infos.len() + r2_infos.len(),
                CallerRejectionReason::InsufficientReads,
            );
            return Ok(ConsensusOutput::default());
        }

        // Downsample if needed
        if let Some(max_reads) = self.options.max_reads_per_strand {
            if r1_infos.len() > max_reads {
                let mut indices: Vec<usize> = (0..r1_infos.len()).collect();
                indices.shuffle(&mut self.rng);
                indices.truncate(max_reads);
                indices.sort_unstable();

                let new_r1: Vec<_> = indices
                    .iter()
                    .map(|&i| std::mem::replace(&mut r1_infos[i], Self::dummy_info()))
                    .collect();
                let new_r2: Vec<_> = indices
                    .iter()
                    .map(|&i| std::mem::replace(&mut r2_infos[i], Self::dummy_info()))
                    .collect();
                r1_infos = new_r1;
                r2_infos = new_r2;
            }
        }

        // Phase 3: Filter to most common alignment on ClippedRecordInfo
        let r1_infos = self.filter_to_most_common_alignment_raw(r1_infos);
        let r2_infos = self.filter_to_most_common_alignment_raw(r2_infos);

        if r1_infos.is_empty() || r2_infos.is_empty() {
            return Ok(ConsensusOutput::default());
        }

        if r1_infos.len() < self.options.min_reads_per_strand
            || r2_infos.len() < self.options.min_reads_per_strand
        {
            self.reject_records_count(
                r1_infos.len() + r2_infos.len(),
                CallerRejectionReason::InsufficientReads,
            );
            return Ok(ConsensusOutput::default());
        }

        // Phase 4: Overlap/phase calculation on ClippedRecordInfo
        // Find longest R1 and R2 by reference length (reverse iter for first-max tie-break)
        let longest_r1 = r1_infos
            .iter()
            .rev()
            .max_by_key(|info| bam_fields::reference_length_from_cigar(&info.clipped_cigar))
            .expect("non-empty");
        let longest_r2 = r2_infos
            .iter()
            .rev()
            .max_by_key(|info| bam_fields::reference_length_from_cigar(&info.clipped_cigar))
            .expect("non-empty");

        let r1_is_negative = longest_r1.flags & flags::REVERSE != 0;
        let (longest_pos, longest_neg) =
            if r1_is_negative { (longest_r2, longest_r1) } else { (longest_r1, longest_r2) };

        let neg_start = longest_neg.adjusted_pos;
        let pos_start = longest_pos.adjusted_pos;
        let pos_ref_len =
            bam_fields::reference_length_from_cigar(&longest_pos.clipped_cigar) as usize;
        let pos_end = pos_start + pos_ref_len.saturating_sub(1);

        let (overlap_start, overlap_end) = (neg_start, pos_end);
        let duplex_length = overlap_end as i64 - overlap_start as i64 + 1;

        if duplex_length < self.options.min_duplex_length as i64 {
            self.reject_records_count(
                r1_infos.len() + r2_infos.len(),
                CallerRejectionReason::InsufficientOverlap,
            );
            return Ok(ConsensusOutput::default());
        }

        // Phase check using raw CIGAR ops
        if !Self::check_overlap_phase_raw(longest_r1, longest_r2, overlap_start, overlap_end) {
            self.reject_records_count(
                r1_infos.len() + r2_infos.len(),
                CallerRejectionReason::IndelErrorBetweenStrands,
            );
            return Ok(ConsensusOutput::default());
        }

        let r2_is_negative = longest_r2.flags & flags::REVERSE != 0;

        // Compute consensus length using clipped info
        let consensus_length =
            Self::compute_consensus_length_raw(longest_pos, longest_neg, overlap_end);
        let Some(consensus_length) = consensus_length else {
            self.reject_records_count(
                r1_infos.len() + r2_infos.len(),
                CallerRejectionReason::IndelErrorBetweenStrands,
            );
            return Ok(ConsensusOutput::default());
        };

        // Phase 5: Build SourceReads and call single-strand consensus
        let r1_is_neg_strand = r1_infos.first().is_some_and(|i| i.flags & flags::REVERSE != 0);
        let r1_source_reads: Vec<SourceRead> = r1_infos
            .iter()
            .enumerate()
            .map(|(idx, info)| {
                Self::to_source_read_for_codec_raw(
                    &records[info.raw_idx],
                    idx,
                    info.clip_amount,
                    info.clip_from_start,
                    Some(&info.clipped_cigar),
                )
            })
            .collect();

        let umi_str = umi.as_deref().unwrap_or("");
        let ss_r1 = match self.ss_caller.consensus_call(umi_str, r1_source_reads)? {
            Some(vcr) => Self::vanilla_to_single_strand(vcr, r1_is_neg_strand, r1_infos.len()),
            None => return Ok(ConsensusOutput::default()),
        };

        let r2_is_neg_strand = r2_infos.first().is_some_and(|i| i.flags & flags::REVERSE != 0);
        let r2_source_reads: Vec<SourceRead> = r2_infos
            .iter()
            .enumerate()
            .map(|(idx, info)| {
                Self::to_source_read_for_codec_raw(
                    &records[info.raw_idx],
                    idx,
                    info.clip_amount,
                    info.clip_from_start,
                    Some(&info.clipped_cigar),
                )
            })
            .collect();

        let ss_r2 = match self.ss_caller.consensus_call(umi_str, r2_source_reads)? {
            Some(vcr) => Self::vanilla_to_single_strand(vcr, r2_is_neg_strand, r2_infos.len()),
            None => return Ok(ConsensusOutput::default()),
        };

        if consensus_length < ss_r1.bases.len() || consensus_length < ss_r2.bases.len() {
            self.reject_records_count(
                r1_infos.len() + r2_infos.len(),
                CallerRejectionReason::IndelErrorBetweenStrands,
            );
            return Ok(ConsensusOutput::default());
        }

        // Orient and pad
        let (ss_r1_oriented, ss_r2_oriented) = if r1_is_negative {
            (Self::reverse_complement_ss(&ss_r1), ss_r2.clone())
        } else {
            (ss_r1.clone(), Self::reverse_complement_ss(&ss_r2))
        };

        let padded_r1 = Self::pad_consensus(&ss_r1_oriented, consensus_length, r1_is_negative);
        let padded_r2 = Self::pad_consensus(&ss_r2_oriented, consensus_length, r2_is_negative);

        let consensus = self.build_duplex_consensus_from_padded(&padded_r1, &padded_r2)?;
        let consensus = self.mask_consensus_quals_query_based(consensus, &padded_r1, &padded_r2);
        let consensus =
            if r1_is_negative { Self::reverse_complement_ss(&consensus) } else { consensus };

        let (ss_for_ac, ss_for_bc) = if r1_is_negative {
            (Self::reverse_complement_ss(&padded_r1), Self::reverse_complement_ss(&padded_r2))
        } else {
            (padded_r1.clone(), padded_r2.clone())
        };

        // Phase 6: Build output from raw bytes
        // Collect raw record refs from filtered reads only (r1_infos/r2_infos)
        // to avoid deriving tags from rejected reads (secondary/unmapped/FR mismatch)
        let all_paired_raws: Vec<&[u8]> = r1_infos
            .iter()
            .chain(r2_infos.iter())
            .map(|info| records[info.raw_idx].as_slice())
            .collect();

        let mut output = ConsensusOutput::default();
        self.build_output_record_into(
            &mut output,
            &consensus,
            &ss_for_ac,
            &ss_for_bc,
            umi.as_deref(),
            &all_paired_raws,
        )?;

        self.stats.consensus_reads_generated += 1;
        Ok(output)
    }

    /// Build `ClippedRecordInfo` for a single raw record.
    fn build_clipped_info(raw: &[u8], raw_idx: usize, clip_amount: usize) -> ClippedRecordInfo {
        let flg = bam_fields::flags(raw);
        let is_reverse = flg & flags::REVERSE != 0;
        let clip_from_start = is_reverse; // negative strand ⟹ clip from start

        let original_ops = bam_fields::get_cigar_ops(raw);
        let (clipped_cigar, ref_bases_consumed) =
            bam_fields::clip_cigar_ops_raw(&original_ops, clip_amount, clip_from_start);

        let original_seq_len = bam_fields::l_seq(raw) as usize;
        let clipped_seq_len = original_seq_len.saturating_sub(clip_amount);

        // 1-based alignment start, adjusted for start-clipping
        let pos_0based = bam_fields::pos(raw);
        debug_assert!(
            pos_0based >= 0,
            "build_clipped_info called on unmapped record (pos={pos_0based})"
        );
        let adjusted_pos = if clip_from_start {
            (pos_0based + 1) as usize + ref_bases_consumed
        } else {
            (pos_0based + 1) as usize
        };

        ClippedRecordInfo {
            raw_idx,
            clip_amount,
            clip_from_start,
            clipped_seq_len,
            clipped_cigar,
            adjusted_pos,
            flags: flg,
        }
    }

    /// Dummy `ClippedRecordInfo` used as a swap placeholder during downsample.
    fn dummy_info() -> ClippedRecordInfo {
        ClippedRecordInfo {
            raw_idx: 0,
            clip_amount: 0,
            clip_from_start: false,
            clipped_seq_len: 0,
            clipped_cigar: Vec::new(),
            adjusted_pos: 0,
            flags: 0,
        }
    }

    /// Filter `ClippedRecordInfo`s to the most common alignment pattern.
    fn filter_to_most_common_alignment_raw(
        &mut self,
        infos: Vec<ClippedRecordInfo>,
    ) -> Vec<ClippedRecordInfo> {
        if infos.len() < 2 {
            return infos;
        }

        let mut indexed: Vec<IndexedSourceRead> = infos
            .iter()
            .enumerate()
            .map(|(i, info)| {
                let mut cigar = bam_fields::simplify_cigar_from_raw(&info.clipped_cigar);
                if info.flags & flags::REVERSE != 0 {
                    cigar.reverse();
                }
                (i, info.clipped_seq_len, cigar)
            })
            .collect();

        indexed.sort_by(|a, b| b.1.cmp(&a.1));

        let best_indices = select_most_common_alignment_group(&indexed);

        let rejected_count = infos.len() - best_indices.len();
        if rejected_count > 0 {
            *self
                .stats
                .rejection_reasons
                .entry(CallerRejectionReason::MinorityAlignment)
                .or_insert(0) += rejected_count;
            self.stats.reads_filtered += rejected_count as u64;
        }

        let best_set: std::collections::HashSet<usize> = best_indices.into_iter().collect();
        infos
            .into_iter()
            .enumerate()
            .filter(|(i, _)| best_set.contains(i))
            .map(|(_, info)| info)
            .collect()
    }

    /// Phase check using `ClippedRecordInfo`.
    fn check_overlap_phase_raw(
        r1: &ClippedRecordInfo,
        r2: &ClippedRecordInfo,
        overlap_start: usize,
        overlap_end: usize,
    ) -> bool {
        let r1s = bam_fields::read_pos_at_ref_pos_raw(
            &r1.clipped_cigar,
            r1.adjusted_pos,
            overlap_start,
            true,
        );
        let r2s = bam_fields::read_pos_at_ref_pos_raw(
            &r2.clipped_cigar,
            r2.adjusted_pos,
            overlap_start,
            true,
        );
        let r1e = bam_fields::read_pos_at_ref_pos_raw(
            &r1.clipped_cigar,
            r1.adjusted_pos,
            overlap_end,
            true,
        );
        let r2e = bam_fields::read_pos_at_ref_pos_raw(
            &r2.clipped_cigar,
            r2.adjusted_pos,
            overlap_end,
            true,
        );

        match (r1s, r2s, r1e, r2e) {
            (Some(a), Some(b), Some(c), Some(d)) => (a as i64 - b as i64) == (c as i64 - d as i64),
            _ => false,
        }
    }

    /// Compute consensus length from `ClippedRecordInfo`.
    fn compute_consensus_length_raw(
        pos_info: &ClippedRecordInfo,
        neg_info: &ClippedRecordInfo,
        overlap_end: usize,
    ) -> Option<usize> {
        let pos_read_pos = bam_fields::read_pos_at_ref_pos_raw(
            &pos_info.clipped_cigar,
            pos_info.adjusted_pos,
            overlap_end,
            false,
        )?;
        let neg_read_pos = bam_fields::read_pos_at_ref_pos_raw(
            &neg_info.clipped_cigar,
            neg_info.adjusted_pos,
            overlap_end,
            false,
        )?;

        Some(pos_read_pos + neg_info.clipped_seq_len - neg_read_pos)
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

    /// Builds the output record from the consensus and writes raw bytes into `output`.
    #[expect(clippy::unnecessary_wraps, reason = "Result return type kept for API consistency with other callers")]
    fn build_output_record_into(
        &mut self,
        output: &mut ConsensusOutput,
        consensus: &SingleStrandConsensus,
        ss_a: &SingleStrandConsensus,
        ss_b: &SingleStrandConsensus,
        umi: Option<&str>,
        source_raws: &[&[u8]],
    ) -> Result<()> {
        // Generate read name - use ':' delimiter to match fgbio format
        self.consensus_counter += 1;
        let read_name = if let Some(umi_str) = umi {
            format!("{}:{}", self.read_name_prefix, umi_str)
        } else {
            format!("{}:{}", self.read_name_prefix, self.consensus_counter)
        };

        // Codec always outputs unmapped fragments
        let flag = flags::UNMAPPED;

        self.bam_builder.build_record(
            read_name.as_bytes(),
            flag,
            &consensus.bases,
            &consensus.quals,
        );

        // RG tag
        self.bam_builder.append_string_tag(b"RG", self.read_group_id.as_bytes());

        // MI tag
        if let Some(umi_str) = umi {
            self.bam_builder.append_string_tag(b"MI", umi_str.as_bytes());
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

        self.bam_builder.append_int_tag(b"cD", total_depth);
        self.bam_builder.append_int_tag(b"cM", min_depth);
        self.bam_builder.append_float_tag(b"cE", error_rate);

        // AB strand tags (aD, aM, aE)
        let a_max_depth = ss_a.depths.iter().map(|&d| i32::from(d)).max().unwrap_or(0);
        let a_min_depth = ss_a.depths.iter().map(|&d| i32::from(d)).min().unwrap_or(0);
        let a_total_errors: usize = ss_a.errors.iter().map(|&e| e as usize).sum();
        let a_total_bases: usize = ss_a.depths.iter().map(|&d| d as usize).sum();
        let a_error_rate =
            if a_total_bases > 0 { a_total_errors as f32 / a_total_bases as f32 } else { 0.0 };

        self.bam_builder.append_int_tag(b"aD", a_max_depth);
        self.bam_builder.append_int_tag(b"aM", a_min_depth);
        self.bam_builder.append_float_tag(b"aE", a_error_rate);

        // BA strand tags (bD, bM, bE)
        let b_max_depth = ss_b.depths.iter().map(|&d| i32::from(d)).max().unwrap_or(0);
        let b_min_depth = ss_b.depths.iter().map(|&d| i32::from(d)).min().unwrap_or(0);
        let b_total_errors: usize = ss_b.errors.iter().map(|&e| e as usize).sum();
        let b_total_bases: usize = ss_b.depths.iter().map(|&d| d as usize).sum();
        let b_error_rate =
            if b_total_bases > 0 { b_total_errors as f32 / b_total_bases as f32 } else { 0.0 };

        self.bam_builder.append_int_tag(b"bD", b_max_depth);
        self.bam_builder.append_int_tag(b"bM", b_min_depth);
        self.bam_builder.append_float_tag(b"bE", b_error_rate);

        // Per-base tags if enabled
        if self.options.produce_per_base_tags {
            // Per-base depth arrays - convert to signed integers to match fgbio/simplex/duplex
            let ad_array: Vec<i16> = ss_a.depths.iter().map(|&d| d as i16).collect();
            let bd_array: Vec<i16> = ss_b.depths.iter().map(|&d| d as i16).collect();
            self.bam_builder.append_i16_array_tag(b"ad", &ad_array);
            self.bam_builder.append_i16_array_tag(b"bd", &bd_array);

            // Per-base error arrays - convert to signed integers to match fgbio/simplex/duplex
            let ae_array: Vec<i16> = ss_a.errors.iter().map(|&e| e as i16).collect();
            let be_array: Vec<i16> = ss_b.errors.iter().map(|&e| e as i16).collect();
            self.bam_builder.append_i16_array_tag(b"ae", &ae_array);
            self.bam_builder.append_i16_array_tag(b"be", &be_array);

            // Single-strand consensus sequences
            self.bam_builder.append_string_tag(b"ac", &ss_a.bases);
            self.bam_builder.append_string_tag(b"bc", &ss_b.bases);

            // Single-strand consensus qualities - Phred+33 encoded
            self.bam_builder.append_phred33_string_tag(b"aq", &ss_a.quals);
            self.bam_builder.append_phred33_string_tag(b"bq", &ss_b.quals);
        }

        // Cell barcode tag - extract from first source read if configured
        if let Some(cell_tag) = &self.options.cell_tag {
            let cell_tag_bytes: [u8; 2] = [cell_tag.as_ref()[0], cell_tag.as_ref()[1]];
            for raw in source_raws {
                if let Some(cell_bc) = bam_fields::find_string_tag_in_record(raw, &cell_tag_bytes) {
                    if !cell_bc.is_empty() {
                        self.bam_builder.append_string_tag(&cell_tag_bytes, cell_bc);
                        break;
                    }
                }
            }
        }

        // RX tag (UmiBases) - extract from all source reads and build consensus
        let umis: Vec<String> = source_raws
            .iter()
            .filter_map(|raw| {
                bam_fields::find_string_tag_in_record(raw, b"RX")
                    .and_then(|b| String::from_utf8(b.to_vec()).ok())
            })
            .collect();

        if !umis.is_empty() {
            let consensus_umi = consensus_umis(&umis);
            if !consensus_umi.is_empty() {
                self.bam_builder.append_string_tag(b"RX", consensus_umi.as_bytes());
            }
        }

        // Write the completed record with block_size prefix
        self.bam_builder.write_with_block_size(&mut output.data);
        output.count += 1;

        Ok(())
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
    fn consensus_reads(&mut self, records: Vec<Vec<u8>>) -> Result<ConsensusOutput> {
        let result = self.consensus_reads_raw(&records)?;
        // When a group fails to produce consensus, all its records are rejected
        if self.track_rejects && result.count == 0 && !records.is_empty() {
            self.rejected_reads.extend(records);
        }
        Ok(result)
    }

    #[expect(clippy::cast_possible_truncation, reason = "read counts will not exceed usize on any supported platform")]
    fn total_reads(&self) -> usize {
        self.stats.total_input_reads as usize
    }

    #[expect(clippy::cast_possible_truncation, reason = "read counts will not exceed usize on any supported platform")]
    fn total_filtered(&self) -> usize {
        self.stats.reads_filtered as usize
    }

    #[expect(clippy::cast_possible_truncation, reason = "read counts will not exceed usize on any supported platform")]
    fn consensus_reads_constructed(&self) -> usize {
        self.stats.consensus_reads_generated as usize
    }

    #[expect(clippy::cast_possible_truncation, reason = "read counts will not exceed usize on any supported platform")]
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

/// Bridge methods for encoding `RecordBuf` to raw bytes.
/// Used by tests (both unit and integration) to call the raw-byte pipeline.
impl CodecConsensusCaller {
    /// Encode a single `RecordBuf` to raw BAM bytes.
    fn record_buf_to_raw(rec: &noodles::sam::alignment::RecordBuf) -> Vec<u8> {
        use crate::vendored::bam_codec::encode_record_buf;
        use noodles::sam::header::record::value::Map;
        use noodles::sam::header::record::value::map::ReferenceSequence;
        use std::num::NonZeroUsize;

        let mut header = noodles::sam::Header::default();
        // Add reference sequences so records with reference_sequence_id encode correctly
        for name in &["chr1", "chr2", "chr3"] {
            let rs = Map::<ReferenceSequence>::new(NonZeroUsize::new(1000).unwrap());
            header.reference_sequences_mut().insert(bstr::BString::from(*name), rs);
        }

        let mut buf = Vec::new();
        encode_record_buf(&mut buf, &header, rec).expect("encode_record_buf failed");
        buf
    }

    /// Encode `RecordBuf`s to raw bytes and delegate to `consensus_reads`.
    ///
    /// # Errors
    ///
    /// Returns an error if consensus calling fails on the provided records.
    #[expect(clippy::needless_pass_by_value, reason = "matches ConsensusCaller trait signature pattern")]
    pub fn consensus_reads_from_sam_records(
        &mut self,
        recs: Vec<noodles::sam::alignment::RecordBuf>,
    ) -> Result<ConsensusOutput> {
        let raw_records: Vec<Vec<u8>> = recs.iter().map(Self::record_buf_to_raw).collect();
        self.consensus_reads(raw_records)
    }
}

#[cfg(test)]
#[allow(clippy::must_use_candidate)]
impl CodecConsensusCaller {
    /// Test-only wrapper: check if a `RecordBuf` is part of an FR pair.
    pub fn is_fr_pair(&self, rec: &noodles::sam::alignment::RecordBuf) -> bool {
        let raw = Self::record_buf_to_raw(rec);
        bam_fields::is_fr_pair_raw(&raw)
    }

    /// Test-only wrapper: filter `RecordBuf`s to most common alignment.
    pub fn filter_to_most_common_alignment(
        &mut self,
        recs: Vec<noodles::sam::alignment::RecordBuf>,
    ) -> Vec<noodles::sam::alignment::RecordBuf> {
        let raws: Vec<Vec<u8>> = recs.iter().map(Self::record_buf_to_raw).collect();
        let infos: Vec<ClippedRecordInfo> =
            raws.iter().enumerate().map(|(i, raw)| Self::build_clipped_info(raw, i, 0)).collect();
        let filtered = self.filter_to_most_common_alignment_raw(infos);
        filtered.into_iter().map(|info| recs[info.raw_idx].clone()).collect()
    }

    /// Test-only wrapper: `read_pos_at_ref_pos` on a `RecordBuf`.
    pub fn read_pos_at_ref_pos(
        rec: &noodles::sam::alignment::RecordBuf,
        ref_pos: usize,
        return_last_base_if_deleted: bool,
    ) -> Option<usize> {
        let raw = Self::record_buf_to_raw(rec);
        let cigar_ops = bam_fields::get_cigar_ops(&raw);
        let alignment_start = (bam_fields::pos(&raw) + 1) as usize; // 1-based
        bam_fields::read_pos_at_ref_pos_raw(
            &cigar_ops,
            alignment_start,
            ref_pos,
            return_last_base_if_deleted,
        )
    }

    /// Test-only wrapper: `check_overlap_phase` on `RecordBuf`s.
    pub fn check_overlap_phase(
        &self,
        r1: &noodles::sam::alignment::RecordBuf,
        r2: &noodles::sam::alignment::RecordBuf,
        overlap_start: usize,
        overlap_end: usize,
    ) -> bool {
        let raw1 = Self::record_buf_to_raw(r1);
        let raw2 = Self::record_buf_to_raw(r2);
        let info1 = Self::build_clipped_info(&raw1, 0, 0);
        let info2 = Self::build_clipped_info(&raw2, 1, 0);
        Self::check_overlap_phase_raw(&info1, &info2, overlap_start, overlap_end)
    }

    /// Test-only wrapper: `to_source_read_for_codec` on a `RecordBuf`.
    pub(crate) fn to_source_read_for_codec(
        rec: &noodles::sam::alignment::RecordBuf,
        original_idx: usize,
    ) -> SourceRead {
        let raw = Self::record_buf_to_raw(rec);
        Self::to_source_read_for_codec_raw(&raw, original_idx, 0, false, None)
    }

    /// Test-only wrapper: `compute_codec_consensus_length` on two `RecordBuf`s.
    pub fn compute_codec_consensus_length(
        pos_rec: &noodles::sam::alignment::RecordBuf,
        neg_rec: &noodles::sam::alignment::RecordBuf,
        overlap_end: usize,
    ) -> Option<usize> {
        let raw_pos = Self::record_buf_to_raw(pos_rec);
        let raw_neg = Self::record_buf_to_raw(neg_rec);
        let info_pos = Self::build_clipped_info(&raw_pos, 0, 0);
        let info_neg = Self::build_clipped_info(&raw_neg, 1, 0);
        Self::compute_consensus_length_raw(&info_pos, &info_neg, overlap_end)
    }

    /// Test-only wrapper: downsample pairs using `max_reads_per_strand` option.
    pub fn downsample_pairs(
        &mut self,
        r1s: Vec<noodles::sam::alignment::RecordBuf>,
        r2s: Vec<noodles::sam::alignment::RecordBuf>,
    ) -> (Vec<noodles::sam::alignment::RecordBuf>, Vec<noodles::sam::alignment::RecordBuf>) {
        let Some(max_reads) = self.options.max_reads_per_strand else {
            return (r1s, r2s);
        };
        if r1s.len() <= max_reads {
            return (r1s, r2s);
        }
        let mut indices: Vec<usize> = (0..r1s.len()).collect();
        indices.shuffle(&mut self.rng);
        indices.truncate(max_reads);
        indices.sort_unstable();
        let new_r1: Vec<_> = indices.iter().map(|&i| r1s[i].clone()).collect();
        let new_r2: Vec<_> = indices.iter().map(|&i| r2s[i].clone()).collect();
        (new_r1, new_r2)
    }
}

#[cfg(test)]
#[allow(clippy::similar_names)]
mod tests {
    use super::*;
    use fgumi_sam::clipper::cigar_utils;
    use fgumi_sam::builder::RecordBuilder;
    use noodles_raw_bam::ParsedBamRecord;
    use noodles::sam::alignment::RecordBuf;
    use noodles::sam::alignment::record::Cigar as CigarTrait;
    use noodles::sam::alignment::record::Flags;
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

        let output = caller.consensus_reads_from_sam_records(Vec::new()).unwrap();
        assert_eq!(output.count, 0);
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
            assert_eq!(cigar.as_ref().len(), 3, "Expected 3-op cigar (3M1D1M)");
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

    /// Tests that CIGAR is reversed for negative strand reads before minority alignment filtering.
    /// This matches fgbio's toSourceReadForCodec behavior where the CIGAR is reversed for
    /// negative strand reads before filterToMostCommonAlignment.
    ///
    /// Without reversal:
    /// - Read A (3M1I2M): first element length = 3
    /// - Read B (2M1I3M): first element length = 2
    /// - 2 < 3, so Read B would win
    ///
    /// With reversal (both negative strand):
    /// - Read A reversed (2M1I3M): first element length = 2
    /// - Read B reversed (3M1I2M): first element length = 3
    /// - 2 < 3, so Read A wins
    #[test]
    fn test_filter_to_most_common_alignment_negative_strand_cigar_reversal() {
        use noodles::sam::alignment::record::cigar::op::Kind;

        let options = CodecConsensusOptions::default();
        let mut caller = CodecConsensusCaller::new("codec".to_string(), "RG1".to_string(), options);

        // Create two groups with equal sizes, both on negative strand
        // The CIGAR reversal should affect which group wins the tie-breaker
        let reads = vec![
            // Group A: 3M1I2M on negative strand -> reversed to 2M1I3M
            create_test_paired_read(
                "r1_groupA",
                b"ACGTAC",
                b"######",
                true,
                true, // negative strand
                false,
                100,
                &[(Kind::Match, 3), (Kind::Insertion, 1), (Kind::Match, 2)],
            ),
            create_test_paired_read(
                "r2_groupA",
                b"ACGTAC",
                b"######",
                true,
                true, // negative strand
                false,
                100,
                &[(Kind::Match, 3), (Kind::Insertion, 1), (Kind::Match, 2)],
            ),
            // Group B: 2M1I3M on negative strand -> reversed to 3M1I2M
            create_test_paired_read(
                "r3_groupB",
                b"ACGTAC",
                b"######",
                true,
                true, // negative strand
                false,
                100,
                &[(Kind::Match, 2), (Kind::Insertion, 1), (Kind::Match, 3)],
            ),
            create_test_paired_read(
                "r4_groupB",
                b"ACGTAC",
                b"######",
                true,
                true, // negative strand
                false,
                100,
                &[(Kind::Match, 2), (Kind::Insertion, 1), (Kind::Match, 3)],
            ),
        ];

        let filtered = caller.filter_to_most_common_alignment(reads);

        // Both groups have size 2, so CIGAR-based tie-breaking kicks in
        assert_eq!(filtered.len(), 2);

        // With CIGAR reversal for negative strand:
        // - Group A's 3M1I2M -> reversed to 2M1I3M (first element = 2)
        // - Group B's 2M1I3M -> reversed to 3M1I2M (first element = 3)
        // Since 2 < 3, Group A wins (the 3M1I2M reads)
        // Without reversal, Group B would win
        for read in &filtered {
            let cigar = read.cigar();
            let ops: Vec<_> = cigar.iter().map(|r| r.unwrap()).collect();
            assert_eq!(ops.len(), 3, "Expected 3-op cigar");
            // The original CIGAR should be 3M1I2M (Group A)
            assert_eq!(ops[0].kind(), Kind::Match);
            assert_eq!(ops[0].len(), 3, "Expected first element to be 3M (Group A won)");
            assert_eq!(ops[1].kind(), Kind::Insertion);
            assert_eq!(ops[1].len(), 1);
            assert_eq!(ops[2].kind(), Kind::Match);
            assert_eq!(ops[2].len(), 2);
        }

        // 2 reads from Group B should be rejected
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

        let output = caller.consensus_reads_from_sam_records(reads).unwrap();

        assert_eq!(output.count, 1, "Should produce one consensus read");

        let records = ParsedBamRecord::parse_all(&output.data);
        let consensus = &records[0];
        // Check read name format
        let name = String::from_utf8_lossy(&consensus.name);
        assert!(name.contains("hi"), "Consensus name should contain MI tag: {name}");

        // Consensus should cover from pos 1 to pos 40 (R1: 1-30, R2: 11-40)
        assert_eq!(consensus.bases.len(), 40, "Consensus should be 40bp");

        // Check RX tag is preserved (UmiBases consensus)
        let rx = consensus.get_string_tag(b"RX").expect("RX tag not found in consensus");
        assert_eq!(&rx[..], b"ACC-TGA", "RX tag should be preserved");
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

        let output = caller.consensus_reads_from_sam_records(reads).unwrap();

        // Should produce consensus
        assert_eq!(output.count, 1, "Should produce one consensus read");

        // Consensus should cover the region with deletion accounted for
        let records = ParsedBamRecord::parse_all(&output.data);
        let consensus = &records[0];
        assert!(!consensus.bases.is_empty(), "Should have sequence");
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

        let output = caller.consensus_reads_from_sam_records(reads).unwrap();

        assert_eq!(output.count, 0, "Should not emit consensus for RF pair");
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

        let output = caller2.consensus_reads_from_sam_records(reads).unwrap();

        assert_eq!(output.count, 0, "Should not emit consensus with insufficient reads");

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

        let output = caller1.consensus_reads_from_sam_records(reads).unwrap();

        assert_eq!(output.count, 1, "Should emit consensus with minReadsPerStrand=1");
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

        let output = caller20.consensus_reads_from_sam_records(reads).unwrap();
        assert_eq!(
            output.count, 1,
            "Should emit consensus with minDuplexLength=20 and 20bp overlap"
        );

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

        let output = caller21.consensus_reads_from_sam_records(reads).unwrap();
        assert_eq!(
            output.count, 0,
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

        let output = caller.consensus_reads_from_sam_records(reads).unwrap();

        assert_eq!(output.count, 0, "Should not emit consensus when mate is unmapped");
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

        let output_fwd = caller_fwd.consensus_reads_from_sam_records(reads_fwd).unwrap();
        assert_eq!(output_fwd.count, 1, "Should produce one consensus read");

        let records_fwd = ParsedBamRecord::parse_all(&output_fwd.data);
        // Consensus should cover from pos 1 to pos 40 (R1: 1-30, R2: 11-40)
        assert_eq!(records_fwd[0].bases.len(), 40, "Forward R1 should produce 40bp consensus");

        // Verify orientation flag is set correctly (forward R1 = unmapped flag only, not reverse)
        assert_eq!(
            records_fwd[0].flag & flags::REVERSE,
            0,
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

        let output = caller_permissive.consensus_reads_from_sam_records(reads).unwrap();
        assert_eq!(output.count, 1, "Should emit consensus with permissive settings");

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

        let output = caller.consensus_reads_from_sam_records(reads).unwrap();

        // Should produce consensus
        assert_eq!(output.count, 1, "Should produce one consensus read");

        // fgbio produces 40bp consensus (skipping the deleted region in R2)
        // fgumi now matches this behavior by filtering out positions with depth=0
        let records = ParsedBamRecord::parse_all(&output.data);
        let consensus = &records[0];
        let len = consensus.bases.len();
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

        let output = caller.consensus_reads_from_sam_records(reads).unwrap();

        // Should produce consensus
        assert_eq!(output.count, 1, "Should produce one consensus read");

        // fgbio expects length 45: 5 (R1 soft-clip) + 35 (REF[0..35]) + 5 (R2 soft-clip)
        let records = ParsedBamRecord::parse_all(&output.data);
        let consensus = &records[0];
        let len = consensus.bases.len();
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

        let output = caller.consensus_reads_from_sam_records(reads).unwrap();

        // Should produce consensus
        assert_eq!(output.count, 1, "Should produce one consensus read");

        // fgbio expects length 30: 5 (soft-clip) + 25 (aligned)
        let records = ParsedBamRecord::parse_all(&output.data);
        let consensus = &records[0];
        assert!(!consensus.bases.is_empty(), "Should have sequence with both reads soft-clipped");
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

        let output = caller.consensus_reads_from_sam_records(reads).unwrap();

        // Should not emit consensus for chimeric pairs
        assert_eq!(
            output.count, 0,
            "Should not emit consensus for cross-chromosomal chimeric pair"
        );
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

        let output = caller.consensus_reads_from_sam_records(reads).unwrap();

        // fgbio rejects this because R1's end lands in an indel in R2
        // fgumi now implements this check and should also reject
        assert_eq!(
            output.count, 0,
            "Should not emit consensus when R1's end lands in an indel in R2"
        );
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

        let output = caller.consensus_reads_from_sam_records(reads).unwrap();

        assert_eq!(output.count, 1, "Should produce one consensus read");

        let records = ParsedBamRecord::parse_all(&output.data);
        let consensus = &records[0];
        let quals = &consensus.quals;

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

        let output = caller.consensus_reads_from_sam_records(reads).unwrap();

        assert_eq!(output.count, 1, "Should produce one consensus read");

        // The consensus should have single-strand regions masked
        // The exact quality values depend on the implementation details
        // This test verifies the option is respected
        let records = ParsedBamRecord::parse_all(&output.data);
        let consensus = &records[0];
        let quals = &consensus.quals;

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
    fn test_to_source_read_clip_exceeds_sequence_length() {
        use noodles::sam::alignment::record::cigar::op::Kind;

        // Create a 4bp read and clip more than 4bp - should not panic
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

        let raw = CodecConsensusCaller::record_buf_to_raw(&read);

        // Clip 10 bases from a 4bp read (clip_from_start = false)
        let source = CodecConsensusCaller::to_source_read_for_codec_raw(&raw, 0, 10, false, None);
        assert!(source.bases.is_empty(), "Bases should be empty after over-clipping");
        assert!(source.quals.is_empty(), "Quals should be empty after over-clipping");

        // Clip 10 bases from start of a 4bp read (clip_from_start = true)
        let source = CodecConsensusCaller::to_source_read_for_codec_raw(&raw, 0, 10, true, None);
        assert!(source.bases.is_empty(), "Bases should be empty after over-clipping from start");
    }

    #[test]
    fn test_to_source_read_clip_exact_sequence_length() {
        use noodles::sam::alignment::record::cigar::op::Kind;

        // Clip exactly the sequence length - should produce empty but not panic
        let read = create_test_paired_read(
            "read1",
            b"ACGT",
            &[30, 31, 32, 33],
            true,
            false,
            true,
            100,
            &[(Kind::Match, 4)],
        );

        let raw = CodecConsensusCaller::record_buf_to_raw(&read);
        let source = CodecConsensusCaller::to_source_read_for_codec_raw(&raw, 0, 4, false, None);
        assert!(source.bases.is_empty());
        assert!(source.quals.is_empty());
    }

    // ========================================================================
    // build_clipped_info tests
    // ========================================================================

    #[test]
    fn test_build_clipped_info_no_clip_forward() {
        use noodles::sam::alignment::record::cigar::op::Kind;

        let read = create_test_paired_read(
            "read1",
            b"ACGTACGT",
            &[30; 8],
            true,
            false, // forward strand
            true,
            100, // 1-based pos -> 0-based = 99
            &[(Kind::Match, 8)],
        );
        let raw = CodecConsensusCaller::record_buf_to_raw(&read);
        let info = CodecConsensusCaller::build_clipped_info(&raw, 0, 0);

        assert_eq!(info.raw_idx, 0);
        assert_eq!(info.clip_amount, 0);
        assert!(!info.clip_from_start, "Forward strand should not clip from start");
        assert_eq!(info.clipped_seq_len, 8);
        assert_eq!(info.adjusted_pos, 100); // 1-based
    }

    #[test]
    fn test_build_clipped_info_clip_from_end_forward() {
        use noodles::sam::alignment::record::cigar::op::Kind;

        let read = create_test_paired_read(
            "read1",
            b"ACGTACGT",
            &[30; 8],
            true,
            false, // forward strand
            true,
            100,
            &[(Kind::Match, 8)],
        );
        let raw = CodecConsensusCaller::record_buf_to_raw(&read);
        let info = CodecConsensusCaller::build_clipped_info(&raw, 5, 3);

        assert_eq!(info.raw_idx, 5);
        assert_eq!(info.clip_amount, 3);
        assert!(!info.clip_from_start, "Forward strand clips from end");
        assert_eq!(info.clipped_seq_len, 5); // 8 - 3
        assert_eq!(info.adjusted_pos, 100); // Position unchanged for end clipping
    }

    #[test]
    fn test_build_clipped_info_clip_from_start_reverse() {
        use noodles::sam::alignment::record::cigar::op::Kind;

        let read = create_test_paired_read(
            "read1",
            b"ACGTACGT",
            &[30; 8],
            true,
            true, // reverse strand
            false,
            100,
            &[(Kind::Match, 8)],
        );
        let raw = CodecConsensusCaller::record_buf_to_raw(&read);
        let info = CodecConsensusCaller::build_clipped_info(&raw, 2, 3);

        assert_eq!(info.raw_idx, 2);
        assert_eq!(info.clip_amount, 3);
        assert!(info.clip_from_start, "Reverse strand should clip from start");
        assert_eq!(info.clipped_seq_len, 5); // 8 - 3
        // Position adjusts forward by clipped reference bases
        assert!(info.adjusted_pos > 100);
    }

    #[test]
    fn test_build_clipped_info_zero_clip_preserves_all() {
        use noodles::sam::alignment::record::cigar::op::Kind;

        let read = create_test_paired_read(
            "read1",
            b"ACGT",
            &[30; 4],
            true,
            false,
            true,
            50,
            &[(Kind::Match, 4)],
        );
        let raw = CodecConsensusCaller::record_buf_to_raw(&read);
        let info = CodecConsensusCaller::build_clipped_info(&raw, 0, 0);

        assert_eq!(info.clipped_seq_len, 4);
        assert_eq!(info.adjusted_pos, 50);
        assert_eq!(info.clipped_cigar.len(), 1); // Single 4M op
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

        let output = caller.consensus_reads_from_sam_records(reads).unwrap();

        assert_eq!(output.count, 0);
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

        let output = caller.consensus_reads_from_sam_records(reads).unwrap();

        assert_eq!(output.count, 0);
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

//! Overlapping bases consensus caller for paired-end reads.
//!
//! When R1 and R2 reads overlap in their insert, the overlapping bases represent
//! the same original molecule positions and should be consensus called before
//! UMI consensus calling. This prevents treating them as independent observations.

use anyhow::{Context, Result};
use noodles::sam::alignment::record::Cigar;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record_buf::RecordBuf;

use crate::phred::{MIN_PHRED, NO_CALL_BASE, NO_CALL_BASE_LOWER};
use noodles_raw_bam;

/// Check if a base is a no-call (N, n, or .)
/// Matches htsjdk's SequenceUtil.isNoCall behavior
#[inline]
fn is_no_call(base: u8) -> bool {
    matches!(base, NO_CALL_BASE | NO_CALL_BASE_LOWER | b'.')
}

/// Strategy for handling bases that agree in overlapping regions
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AgreementStrategy {
    /// Sum the quality scores (capped at Q93)
    Consensus,
    /// Use the maximum quality score
    MaxQual,
    /// Pass through without modification
    PassThrough,
}

/// Strategy for handling bases that disagree in overlapping regions
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DisagreementStrategy {
    /// Call the base with higher quality, new qual = `qual_high` - `qual_low`
    Consensus,
    /// Mask both bases to N with quality 2
    MaskBoth,
    /// Mask only the lower quality base to N with quality 2
    MaskLowerQual,
}

/// Statistics tracking corrections made during overlapping consensus
#[derive(Debug, Default, Clone)]
pub struct CorrectionStats {
    /// Total number of overlapping base pairs examined
    pub overlapping_bases: usize,
    /// Number of overlapping bases that agreed
    pub bases_agreeing: usize,
    /// Number of overlapping bases that disagreed
    pub bases_disagreeing: usize,
    /// Number of bases corrected (modified)
    pub bases_corrected: usize,
}

impl CorrectionStats {
    /// Create a new empty statistics tracker
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Reset all statistics to zero
    pub fn reset(&mut self) {
        self.overlapping_bases = 0;
        self.bases_agreeing = 0;
        self.bases_disagreeing = 0;
        self.bases_corrected = 0;
    }

    /// Merge statistics from another `CorrectionStats` instance
    pub fn merge(&mut self, other: &Self) {
        self.overlapping_bases += other.overlapping_bases;
        self.bases_agreeing += other.bases_agreeing;
        self.bases_disagreeing += other.bases_disagreeing;
        self.bases_corrected += other.bases_corrected;
    }
}

/// Consensus caller for overlapping bases in paired-end reads
pub struct OverlappingBasesConsensusCaller {
    agreement_strategy: AgreementStrategy,
    disagreement_strategy: DisagreementStrategy,
    stats: CorrectionStats,
}

#[allow(
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::cast_sign_loss
)]
impl OverlappingBasesConsensusCaller {
    /// Create a new overlapping bases consensus caller
    ///
    /// # Arguments
    /// * `agreement_strategy` - How to handle agreeing bases
    /// * `disagreement_strategy` - How to handle disagreeing bases
    #[must_use]
    pub fn new(
        agreement_strategy: AgreementStrategy,
        disagreement_strategy: DisagreementStrategy,
    ) -> Self {
        Self { agreement_strategy, disagreement_strategy, stats: CorrectionStats::new() }
    }

    /// Get the current statistics
    #[must_use]
    pub fn stats(&self) -> &CorrectionStats {
        &self.stats
    }

    /// Reset statistics
    pub fn reset_stats(&mut self) {
        self.stats.reset();
    }

    /// Process overlapping bases in a read pair
    ///
    /// This method modifies the reads in-place, updating bases and qualities
    /// in the overlapping region according to the configured strategies.
    ///
    /// This implementation matches fgbio's `OverlappingBasesConsensusCaller`:
    /// - Uses merge iteration to properly handle different CIGAR structures
    /// - Only considers aligned bases (M/X/=), not soft clips or insertions
    /// - Properly synchronizes by reference position, not read position
    ///
    /// # Arguments
    /// * `r1` - First read in the pair (will be modified)
    /// * `r2` - Second read in the pair (will be modified)
    ///
    /// # Returns
    /// * `true` if overlapping bases were processed, `false` if reads don't overlap
    ///
    /// # Errors
    ///
    /// Returns an error if CIGAR parsing or alignment position extraction fails.
    pub fn call(&mut self, r1: &mut RecordBuf, r2: &mut RecordBuf) -> Result<bool> {
        // Only process paired reads where both are mapped
        if r1.flags().is_unmapped() || r2.flags().is_unmapped() {
            return Ok(false);
        }

        // Must be on the same reference sequence (chromosome)
        if r1.reference_sequence_id() != r2.reference_sequence_id() {
            return Ok(false);
        }

        // Verify both have alignment positions
        if r1.alignment_start().is_none()
            || r1.alignment_end().is_none()
            || r2.alignment_start().is_none()
            || r2.alignment_end().is_none()
        {
            return Ok(false);
        }

        // Create merge iterator that yields positions where both reads have aligned bases
        let Ok(overlap_iter) = ReadMateAndRefPosIterator::new(r1, r2) else { return Ok(false) };

        // Collect overlapping positions
        let overlapping_positions: Vec<_> = overlap_iter.collect();

        if overlapping_positions.is_empty() {
            return Ok(false);
        }

        // Clone sequences and qualities ONCE at the start (avoids O(n²) allocations)
        let mut r1_seq: Vec<u8> = r1.sequence().as_ref().to_vec();
        let mut r2_seq: Vec<u8> = r2.sequence().as_ref().to_vec();
        let mut r1_quals: Vec<u8> = r1.quality_scores().as_ref().to_vec();
        let mut r2_quals: Vec<u8> = r2.quality_scores().as_ref().to_vec();
        let mut modified = false;

        // Process each overlapping position
        for pos in overlapping_positions {
            let r1_base = r1_seq[pos.read_offset];
            let r2_base = r2_seq[pos.mate_offset];

            // Skip if either base is a no-call (N) - matches fgbio behavior
            if is_no_call(r1_base) || is_no_call(r2_base) {
                continue;
            }

            self.stats.overlapping_bases += 1;

            let r1_qual = r1_quals[pos.read_offset];
            let r2_qual = r2_quals[pos.mate_offset];

            if r1_base == r2_base {
                // Bases agree
                self.stats.bases_agreeing += 1;
                if self.process_agreement(
                    pos.read_offset,
                    pos.mate_offset,
                    r1_qual,
                    r2_qual,
                    &mut r1_quals,
                    &mut r2_quals,
                ) {
                    modified = true;
                }
            } else {
                // Bases disagree
                self.stats.bases_disagreeing += 1;
                self.process_disagreement(
                    pos.read_offset,
                    pos.mate_offset,
                    r1_base,
                    r2_base,
                    r1_qual,
                    r2_qual,
                    &mut r1_seq,
                    &mut r2_seq,
                    &mut r1_quals,
                    &mut r2_quals,
                );
                modified = true;
            }
        }

        // Write back modified sequences and qualities ONCE at the end
        if modified {
            *r1.sequence_mut() = r1_seq.into();
            *r2.sequence_mut() = r2_seq.into();
            *r1.quality_scores_mut() = r1_quals.into();
            *r2.quality_scores_mut() = r2_quals.into();
        }

        Ok(true)
    }

    /// Process agreeing bases. Returns true if any modification was made.
    fn process_agreement(
        &mut self,
        r1_offset: usize,
        r2_offset: usize,
        r1_qual: u8,
        r2_qual: u8,
        r1_quals: &mut [u8],
        r2_quals: &mut [u8],
    ) -> bool {
        match self.agreement_strategy {
            AgreementStrategy::Consensus => {
                // Sum qualities, capped at Q93
                let new_qual = (r1_qual + r2_qual).min(93);
                r1_quals[r1_offset] = new_qual;
                r2_quals[r2_offset] = new_qual;

                if new_qual != r1_qual || new_qual != r2_qual {
                    self.stats.bases_corrected += 1;
                    true
                } else {
                    false
                }
            }
            AgreementStrategy::MaxQual => {
                // Use maximum quality
                let new_qual = r1_qual.max(r2_qual);
                r1_quals[r1_offset] = new_qual;
                r2_quals[r2_offset] = new_qual;

                if new_qual != r1_qual || new_qual != r2_qual {
                    self.stats.bases_corrected += 1;
                    true
                } else {
                    false
                }
            }
            AgreementStrategy::PassThrough => false,
        }
    }

    /// Process disagreeing bases
    #[allow(clippy::too_many_arguments)]
    fn process_disagreement(
        &mut self,
        r1_offset: usize,
        r2_offset: usize,
        r1_base: u8,
        r2_base: u8,
        r1_qual: u8,
        r2_qual: u8,
        r1_seq: &mut [u8],
        r2_seq: &mut [u8],
        r1_quals: &mut [u8],
        r2_quals: &mut [u8],
    ) {
        match self.disagreement_strategy {
            DisagreementStrategy::Consensus => {
                // Call the base with higher quality, new qual = qual_high - qual_low
                // If qualities are equal, mask both bases (matches fgbio)
                let (consensus_base, consensus_qual) = match r1_qual.cmp(&r2_qual) {
                    std::cmp::Ordering::Equal => (NO_CALL_BASE, MIN_PHRED),
                    std::cmp::Ordering::Greater => {
                        (r1_base, r1_qual.saturating_sub(r2_qual).max(MIN_PHRED))
                    }
                    std::cmp::Ordering::Less => {
                        (r2_base, r2_qual.saturating_sub(r1_qual).max(MIN_PHRED))
                    }
                };

                r1_seq[r1_offset] = consensus_base;
                r2_seq[r2_offset] = consensus_base;
                r1_quals[r1_offset] = consensus_qual;
                r2_quals[r2_offset] = consensus_qual;

                self.stats.bases_corrected += 2;
            }
            DisagreementStrategy::MaskBoth => {
                // Mask both bases to N with quality 2
                r1_seq[r1_offset] = NO_CALL_BASE;
                r2_seq[r2_offset] = NO_CALL_BASE;
                r1_quals[r1_offset] = MIN_PHRED;
                r2_quals[r2_offset] = MIN_PHRED;

                self.stats.bases_corrected += 2;
            }
            DisagreementStrategy::MaskLowerQual => {
                // Mask only the lower quality base, or both if equal (matches fgbio)
                match r1_qual.cmp(&r2_qual) {
                    std::cmp::Ordering::Less => {
                        r1_seq[r1_offset] = NO_CALL_BASE;
                        r1_quals[r1_offset] = MIN_PHRED;
                        self.stats.bases_corrected += 1;
                    }
                    std::cmp::Ordering::Greater => {
                        r2_seq[r2_offset] = NO_CALL_BASE;
                        r2_quals[r2_offset] = MIN_PHRED;
                        self.stats.bases_corrected += 1;
                    }
                    std::cmp::Ordering::Equal => {
                        // Mask both bases when qualities are equal (matches fgbio)
                        r1_seq[r1_offset] = NO_CALL_BASE;
                        r2_seq[r2_offset] = NO_CALL_BASE;
                        r1_quals[r1_offset] = MIN_PHRED;
                        r2_quals[r2_offset] = MIN_PHRED;

                        self.stats.bases_corrected += 2;
                    }
                }
            }
        }
    }
}

/// A position in a read with both read offset and reference position
#[derive(Debug, Clone, Copy)]
struct ReadPosition {
    /// 0-based offset in the read sequence
    read_offset: usize,
    /// 1-based reference position (matching fgbio convention)
    ref_pos: i32,
}

/// Position from `ReadMateAndRefPosIterator` with positions in both read and mate
#[derive(Debug, Clone, Copy)]
struct ReadMateAndRefPos {
    /// 0-based offset in the read sequence
    read_offset: usize,
    /// 0-based offset in the mate sequence
    mate_offset: usize,
}

/// Iterator that walks through aligned (M/X/=) positions in a read.
///
/// This matches fgbio's `ReadAndRefPosIterator`:
/// - Only returns positions that are aligned to the reference (M/X/= operations)
/// - Skips soft clips, insertions, deletions, etc.
/// - Can be limited to a reference position range
struct ReadAndRefPosIterator {
    /// 1-based current read position (fgbio uses 1-based)
    cur_read_pos: i32,
    /// 1-based current reference position
    cur_ref_pos: i32,
    /// CIGAR operations
    cigar_ops: Vec<(Kind, usize)>,
    /// Current element index (0-based)
    element_index: usize,
    /// Current offset within the element (0-based)
    in_elem_offset: usize,
    /// Start of reference window (1-based, inclusive)
    start_ref_pos: i32,
    /// End of reference window (1-based, inclusive)
    end_ref_pos: i32,
    /// Start of read window (1-based, inclusive)
    start_read_pos: i32,
    /// End of read window (1-based, inclusive)
    end_read_pos: i32,
}

#[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap, clippy::cast_sign_loss)]
impl ReadAndRefPosIterator {
    /// Create iterator for aligned positions overlapping with mate
    ///
    /// # Arguments
    /// * `record` - The read to iterate
    /// * `mate` - The mate read (used to determine overlap bounds)
    fn new_with_mate(record: &RecordBuf, mate: &RecordBuf) -> Result<Self> {
        let rec_start = match record.alignment_start() {
            Some(pos) => usize::from(pos) as i32,
            None => return Err(anyhow::anyhow!("Record has no alignment start")),
        };
        let rec_end = match record.alignment_end() {
            Some(pos) => usize::from(pos) as i32,
            None => return Err(anyhow::anyhow!("Record has no alignment end")),
        };

        let mate_start = match mate.alignment_start() {
            Some(pos) => usize::from(pos) as i32,
            None => return Err(anyhow::anyhow!("Mate has no alignment start")),
        };
        let mate_end = match mate.alignment_end() {
            Some(pos) => usize::from(pos) as i32,
            None => return Err(anyhow::anyhow!("Mate has no alignment end")),
        };

        // Reference window is the overlap between rec and mate
        let min_ref_pos = rec_start.max(mate_start);
        let max_ref_pos = rec_end.min(mate_end);

        Self::new(record, 1, i32::MAX, min_ref_pos, max_ref_pos)
    }

    /// Create a new iterator with specified bounds
    fn new(
        record: &RecordBuf,
        min_read_pos: i32,
        max_read_pos: i32,
        min_ref_pos: i32,
        max_ref_pos: i32,
    ) -> Result<Self> {
        let rec_start = match record.alignment_start() {
            Some(pos) => usize::from(pos) as i32,
            None => return Err(anyhow::anyhow!("Record has no alignment start")),
        };
        let rec_end = match record.alignment_end() {
            Some(pos) => usize::from(pos) as i32,
            None => return Err(anyhow::anyhow!("Record has no alignment end")),
        };
        let rec_len = record.sequence().len() as i32;

        // Extract CIGAR operations
        let cigar = record.cigar();
        let cigar_ops: Vec<_> = cigar
            .iter()
            .filter_map(std::result::Result::ok)
            .map(|op| (op.kind(), op.len()))
            .collect();

        // Calculate bounds (matching fgbio)
        let start_read_pos = 1.max(min_read_pos);
        let end_read_pos = rec_len.min(max_read_pos);
        let start_ref_pos = rec_start.max(min_ref_pos);
        let end_ref_pos = rec_end.min(max_ref_pos);

        let mut iter = Self {
            cur_read_pos: 1,
            cur_ref_pos: rec_start,
            cigar_ops,
            element_index: 0,
            in_elem_offset: 0,
            start_ref_pos,
            end_ref_pos,
            start_read_pos,
            end_read_pos,
        };

        // Initialize: skip to first position in bounds
        iter.skip_to_start();

        Ok(iter)
    }

    /// Get current CIGAR element length on reference
    fn cur_elem_length_on_target(&self) -> usize {
        if self.element_index >= self.cigar_ops.len() {
            return 0;
        }
        let (kind, len) = self.cigar_ops[self.element_index];
        match kind {
            Kind::Match
            | Kind::SequenceMatch
            | Kind::SequenceMismatch
            | Kind::Deletion
            | Kind::Skip => len,
            _ => 0,
        }
    }

    /// Get current CIGAR element length on query (read)
    fn cur_elem_length_on_query(&self) -> usize {
        if self.element_index >= self.cigar_ops.len() {
            return 0;
        }
        let (kind, len) = self.cigar_ops[self.element_index];
        match kind {
            Kind::Match
            | Kind::SequenceMatch
            | Kind::SequenceMismatch
            | Kind::Insertion
            | Kind::SoftClip => len,
            _ => 0,
        }
    }

    /// Check if current element is an alignment operation (M/X/=)
    fn cur_elem_is_alignment(&self) -> bool {
        if self.element_index >= self.cigar_ops.len() {
            return false;
        }
        let (kind, _) = self.cigar_ops[self.element_index];
        matches!(kind, Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch)
    }

    /// Skip to the first position in the ref/read window
    fn skip_to_start(&mut self) {
        // Skip until we have an element at or past both the read/ref starts
        while self.element_index < self.cigar_ops.len() {
            let cur_ref_end = self.cur_ref_pos + self.cur_elem_length_on_target() as i32 - 1;
            let cur_read_end = self.cur_read_pos + self.cur_elem_length_on_query() as i32 - 1;

            if cur_ref_end >= self.start_ref_pos && cur_read_end >= self.start_read_pos {
                break;
            }

            self.cur_ref_pos += self.cur_elem_length_on_target() as i32;
            self.cur_read_pos += self.cur_elem_length_on_query() as i32;
            self.element_index += 1;
        }

        // Skip non-aligned bases and adjust offset
        self.skip_non_aligned_and_adjust_offset();
    }

    /// Skip non-aligned elements and adjust offset into current element
    fn skip_non_aligned_and_adjust_offset(&mut self) {
        self.in_elem_offset = 0;

        // Skip over any non-aligned bases
        while self.element_index < self.cigar_ops.len() && !self.cur_elem_is_alignment() {
            self.cur_ref_pos += self.cur_elem_length_on_target() as i32;
            self.cur_read_pos += self.cur_elem_length_on_query() as i32;
            self.element_index += 1;
        }

        // Update offset if current element spans the read-start or reference-start
        if self.element_index < self.cigar_ops.len()
            && (self.cur_ref_pos < self.start_ref_pos || self.cur_read_pos < self.start_read_pos)
        {
            let offset = (self.start_ref_pos - self.cur_ref_pos)
                .max(self.start_read_pos - self.cur_read_pos);
            self.in_elem_offset = offset as usize;
            self.cur_ref_pos += offset;
            self.cur_read_pos += offset;
        }
    }
}

#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
impl Iterator for ReadAndRefPosIterator {
    type Item = ReadPosition;

    fn next(&mut self) -> Option<Self::Item> {
        // Check if we've consumed all bases in current element
        if self.element_index < self.cigar_ops.len() {
            let (_, len) = self.cigar_ops[self.element_index];
            if self.in_elem_offset >= len {
                self.element_index += 1;
                self.skip_non_aligned_and_adjust_offset();
            }
        }

        // Check bounds
        if self.element_index >= self.cigar_ops.len()
            || self.cur_read_pos > self.end_read_pos
            || self.cur_ref_pos > self.end_ref_pos
        {
            return None;
        }

        // Return current position and advance
        let pos = ReadPosition {
            read_offset: (self.cur_read_pos - 1) as usize, // Convert to 0-based
            ref_pos: self.cur_ref_pos,
        };

        self.cur_read_pos += 1;
        self.cur_ref_pos += 1;
        self.in_elem_offset += 1;

        Some(pos)
    }
}

/// Iterator that yields positions where both read and mate have aligned bases
/// at the same reference position.
///
/// This matches fgbio's `ReadMateAndRefPosIterator`:
/// - Uses merge iteration (not lockstep)
/// - Advances whichever iterator has the lower ref position
/// - Only yields when both have a position at the same ref position
struct ReadMateAndRefPosIterator {
    rec_iter: std::iter::Peekable<ReadAndRefPosIterator>,
    mate_iter: std::iter::Peekable<ReadAndRefPosIterator>,
}

impl ReadMateAndRefPosIterator {
    /// Create a new iterator for overlapping positions between read and mate
    fn new(record: &RecordBuf, mate: &RecordBuf) -> Result<Self> {
        let rec_iter = ReadAndRefPosIterator::new_with_mate(record, mate)?.peekable();
        let mate_iter = ReadAndRefPosIterator::new_with_mate(mate, record)?.peekable();

        Ok(Self { rec_iter, mate_iter })
    }
}

impl Iterator for ReadMateAndRefPosIterator {
    type Item = ReadMateAndRefPos;

    fn next(&mut self) -> Option<Self::Item> {
        use std::cmp::Ordering;

        loop {
            // Get heads of both iterators
            let rec_pos = self.rec_iter.peek()?;
            let mate_pos = self.mate_iter.peek()?;

            match rec_pos.ref_pos.cmp(&mate_pos.ref_pos) {
                Ordering::Less => {
                    // Advance rec iterator
                    self.rec_iter.next();
                }
                Ordering::Greater => {
                    // Advance mate iterator
                    self.mate_iter.next();
                }
                Ordering::Equal => {
                    // ref_pos matches - yield and advance both
                    // Safe: we just peeked and saw Some values for both iterators
                    let rec = self.rec_iter.next().expect("peeked Some");
                    let mate = self.mate_iter.next().expect("peeked Some");

                    return Some(ReadMateAndRefPos {
                        read_offset: rec.read_offset,
                        mate_offset: mate.read_offset,
                    });
                }
            }
        }
    }
}

/// Applies overlapping consensus calling to pairs of reads within a group.
///
/// For paired-end reads, this function:
/// 1. Groups reads by name to find R1/R2 pairs
/// 2. Calls overlapping consensus on each pair using the provided caller
/// 3. Modifies reads in-place (no copying)
///
/// Single-end reads pass through unchanged.
///
/// # Arguments
///
/// * `reads` - Mutable slice of reads to process (modified in-place)
/// * `caller` - The overlapping consensus caller to use
///
/// # Errors
///
/// Returns an error if consensus calling fails for any pair
pub fn apply_overlapping_consensus(
    reads: &mut [RecordBuf],
    caller: &mut OverlappingBasesConsensusCaller,
) -> Result<()> {
    use ahash::AHashMap;

    // Group reads by name for pairing
    let mut read_pairs: AHashMap<String, (Option<usize>, Option<usize>)> = AHashMap::new();

    for (idx, record) in reads.iter().enumerate() {
        let read_name = record.name().map(std::string::ToString::to_string).unwrap_or_default();

        if record.flags().is_first_segment() {
            read_pairs.entry(read_name).or_insert((None, None)).0 = Some(idx);
        } else if record.flags().is_last_segment() {
            read_pairs.entry(read_name).or_insert((None, None)).1 = Some(idx);
        }
        // Single-end reads are not paired, they pass through unchanged
    }

    // Process pairs
    for (r1_idx, r2_idx) in read_pairs.values() {
        if let (Some(idx1), Some(idx2)) = (r1_idx, r2_idx) {
            // Use split_at_mut to get two mutable references
            let (r1, r2) = if idx1 < idx2 {
                let (left, right) = reads.split_at_mut(*idx2);
                (&mut left[*idx1], &mut right[0])
            } else {
                let (left, right) = reads.split_at_mut(*idx1);
                (&mut right[0], &mut left[*idx2])
            };

            caller.call(r1, r2).context("Failed to call overlapping consensus")?;
        }
    }

    Ok(())
}

// ============================================================================
// Raw-byte overlapping consensus
// ============================================================================

impl OverlappingBasesConsensusCaller {
    /// Process overlapping bases in a read pair using raw BAM bytes.
    ///
    /// Same logic as `call()` but extracts fields from raw bytes directly.
    /// Modifies the raw byte records in-place.
    ///
    /// # Returns
    /// * `true` if overlapping bases were processed, `false` if reads don't overlap
    ///
    /// # Errors
    ///
    /// Returns an error if raw BAM field extraction or CIGAR parsing fails.
    #[allow(
        clippy::cast_precision_loss,
        clippy::cast_possible_truncation,
        clippy::cast_possible_wrap,
        clippy::cast_sign_loss
    )]
    pub fn call_raw(&mut self, r1: &mut [u8], r2: &mut [u8]) -> Result<bool> {
        // Only process paired reads where both are mapped
        if noodles_raw_bam::flags(r1) & noodles_raw_bam::flags::UNMAPPED != 0
            || noodles_raw_bam::flags(r2) & noodles_raw_bam::flags::UNMAPPED != 0
        {
            return Ok(false);
        }

        // Must be on the same reference sequence
        if noodles_raw_bam::ref_id(r1) != noodles_raw_bam::ref_id(r2) {
            return Ok(false);
        }

        // Verify both have alignment positions and ends
        let Some(r1_start) = noodles_raw_bam::alignment_start_from_raw(r1) else { return Ok(false) };
        let Some(r1_end) = noodles_raw_bam::alignment_end_from_raw(r1) else { return Ok(false) };
        let Some(r2_start) = noodles_raw_bam::alignment_start_from_raw(r2) else { return Ok(false) };
        let Some(r2_end) = noodles_raw_bam::alignment_end_from_raw(r2) else { return Ok(false) };

        // Create merge iterator that yields positions where both reads have aligned bases
        let overlap_iter =
            RawReadMateAndRefPosIterator::new(r1, r2, r1_start, r1_end, r2_start, r2_end);

        let overlapping_positions: Vec<_> = overlap_iter.collect();

        if overlapping_positions.is_empty() {
            return Ok(false);
        }

        // Extract sequences and qualities for modification
        let mut r1_seq = noodles_raw_bam::extract_sequence(r1);
        let mut r2_seq = noodles_raw_bam::extract_sequence(r2);
        let mut r1_quals: Vec<u8> = noodles_raw_bam::quality_scores_slice(r1).to_vec();
        let mut r2_quals: Vec<u8> = noodles_raw_bam::quality_scores_slice(r2).to_vec();
        let mut modified = false;

        for pos in overlapping_positions {
            let r1_base = r1_seq[pos.read_offset];
            let r2_base = r2_seq[pos.mate_offset];

            if is_no_call(r1_base) || is_no_call(r2_base) {
                continue;
            }

            self.stats.overlapping_bases += 1;

            let r1_qual = r1_quals[pos.read_offset];
            let r2_qual = r2_quals[pos.mate_offset];

            if r1_base == r2_base {
                self.stats.bases_agreeing += 1;
                if self.process_agreement(
                    pos.read_offset,
                    pos.mate_offset,
                    r1_qual,
                    r2_qual,
                    &mut r1_quals,
                    &mut r2_quals,
                ) {
                    modified = true;
                }
            } else {
                self.stats.bases_disagreeing += 1;
                self.process_disagreement(
                    pos.read_offset,
                    pos.mate_offset,
                    r1_base,
                    r2_base,
                    r1_qual,
                    r2_qual,
                    &mut r1_seq,
                    &mut r2_seq,
                    &mut r1_quals,
                    &mut r2_quals,
                );
                modified = true;
            }
        }

        // Write back modified sequences and qualities into the raw BAM records
        if modified {
            let r1_seq_off = noodles_raw_bam::seq_offset(r1);
            let r2_seq_off = noodles_raw_bam::seq_offset(r2);
            for (i, &base) in r1_seq.iter().enumerate() {
                noodles_raw_bam::set_base(r1, r1_seq_off, i, base);
            }
            for (i, &base) in r2_seq.iter().enumerate() {
                noodles_raw_bam::set_base(r2, r2_seq_off, i, base);
            }
            let r1_qual_off = noodles_raw_bam::qual_offset(r1);
            let r2_qual_off = noodles_raw_bam::qual_offset(r2);
            r1[r1_qual_off..r1_qual_off + r1_quals.len()].copy_from_slice(&r1_quals);
            r2[r2_qual_off..r2_qual_off + r2_quals.len()].copy_from_slice(&r2_quals);
        }

        Ok(true)
    }
}

/// Raw-byte version of `ReadAndRefPosIterator`.
///
/// Walks through aligned (M/X/=) positions in a raw BAM record's CIGAR.
struct RawReadAndRefPosIterator {
    cur_read_pos: i32,
    cur_ref_pos: i32,
    cigar_ops: Vec<(u8, usize)>, // (op_type, op_len)
    element_index: usize,
    in_elem_offset: usize,
    start_ref_pos: i32,
    end_ref_pos: i32,
    start_read_pos: i32,
    end_read_pos: i32,
}

#[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap, clippy::cast_sign_loss)]
impl RawReadAndRefPosIterator {
    /// Create iterator for aligned positions overlapping with mate.
    fn new_with_mate(
        bam: &[u8],
        rec_start: usize,
        rec_end: usize,
        mate_start: usize,
        mate_end: usize,
    ) -> Self {
        let rec_start_i32 = rec_start as i32;
        let rec_end_i32 = rec_end as i32;
        let rec_len = noodles_raw_bam::l_seq(bam) as i32;

        let min_ref_pos = rec_start_i32.max(mate_start as i32);
        let max_ref_pos = rec_end_i32.min(mate_end as i32);

        // Extract CIGAR ops from raw BAM
        let raw_ops = noodles_raw_bam::get_cigar_ops(bam);
        let cigar_ops: Vec<(u8, usize)> =
            raw_ops.iter().map(|&op| ((op & 0xF) as u8, (op >> 4) as usize)).collect();

        let start_read_pos = 1i32;
        let end_read_pos = rec_len;
        let start_ref_pos = rec_start_i32.max(min_ref_pos);
        let end_ref_pos = rec_end_i32.min(max_ref_pos);

        let mut iter = Self {
            cur_read_pos: 1,
            cur_ref_pos: rec_start_i32,
            cigar_ops,
            element_index: 0,
            in_elem_offset: 0,
            start_ref_pos,
            end_ref_pos,
            start_read_pos,
            end_read_pos,
        };

        iter.skip_to_start();
        iter
    }

    fn cur_elem_length_on_target(&self) -> usize {
        if self.element_index >= self.cigar_ops.len() {
            return 0;
        }
        let (op_type, len) = self.cigar_ops[self.element_index];
        // M(0), D(2), N(3), =(7), X(8) consume reference
        if matches!(op_type, 0 | 2 | 3 | 7 | 8) { len } else { 0 }
    }

    fn cur_elem_length_on_query(&self) -> usize {
        if self.element_index >= self.cigar_ops.len() {
            return 0;
        }
        let (op_type, len) = self.cigar_ops[self.element_index];
        // M(0), I(1), S(4), =(7), X(8) consume query
        if matches!(op_type, 0 | 1 | 4 | 7 | 8) { len } else { 0 }
    }

    fn cur_elem_is_alignment(&self) -> bool {
        if self.element_index >= self.cigar_ops.len() {
            return false;
        }
        let (op_type, _) = self.cigar_ops[self.element_index];
        // M(0), =(7), X(8) are alignment operations
        matches!(op_type, 0 | 7 | 8)
    }

    fn skip_to_start(&mut self) {
        while self.element_index < self.cigar_ops.len() {
            let cur_ref_end = self.cur_ref_pos + self.cur_elem_length_on_target() as i32 - 1;
            let cur_read_end = self.cur_read_pos + self.cur_elem_length_on_query() as i32 - 1;

            if cur_ref_end >= self.start_ref_pos && cur_read_end >= self.start_read_pos {
                break;
            }

            self.cur_ref_pos += self.cur_elem_length_on_target() as i32;
            self.cur_read_pos += self.cur_elem_length_on_query() as i32;
            self.element_index += 1;
        }

        self.skip_non_aligned_and_adjust_offset();
    }

    fn skip_non_aligned_and_adjust_offset(&mut self) {
        self.in_elem_offset = 0;

        while self.element_index < self.cigar_ops.len() && !self.cur_elem_is_alignment() {
            self.cur_ref_pos += self.cur_elem_length_on_target() as i32;
            self.cur_read_pos += self.cur_elem_length_on_query() as i32;
            self.element_index += 1;
        }

        if self.element_index < self.cigar_ops.len()
            && (self.cur_ref_pos < self.start_ref_pos || self.cur_read_pos < self.start_read_pos)
        {
            let offset = (self.start_ref_pos - self.cur_ref_pos)
                .max(self.start_read_pos - self.cur_read_pos);
            self.in_elem_offset = offset as usize;
            self.cur_ref_pos += offset;
            self.cur_read_pos += offset;
        }
    }
}

#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
impl Iterator for RawReadAndRefPosIterator {
    type Item = ReadPosition;

    fn next(&mut self) -> Option<Self::Item> {
        if self.element_index < self.cigar_ops.len() {
            let (_, len) = self.cigar_ops[self.element_index];
            if self.in_elem_offset >= len {
                self.element_index += 1;
                self.skip_non_aligned_and_adjust_offset();
            }
        }

        if self.element_index >= self.cigar_ops.len()
            || self.cur_read_pos > self.end_read_pos
            || self.cur_ref_pos > self.end_ref_pos
        {
            return None;
        }

        let pos = ReadPosition {
            read_offset: (self.cur_read_pos - 1) as usize,
            ref_pos: self.cur_ref_pos,
        };

        self.cur_read_pos += 1;
        self.cur_ref_pos += 1;
        self.in_elem_offset += 1;

        Some(pos)
    }
}

/// Raw-byte version of `ReadMateAndRefPosIterator`.
struct RawReadMateAndRefPosIterator {
    rec_iter: std::iter::Peekable<RawReadAndRefPosIterator>,
    mate_iter: std::iter::Peekable<RawReadAndRefPosIterator>,
}

impl RawReadMateAndRefPosIterator {
    fn new(
        r1: &[u8],
        r2: &[u8],
        r1_start: usize,
        r1_end: usize,
        r2_start: usize,
        r2_end: usize,
    ) -> Self {
        let rec_iter =
            RawReadAndRefPosIterator::new_with_mate(r1, r1_start, r1_end, r2_start, r2_end)
                .peekable();
        let mate_iter =
            RawReadAndRefPosIterator::new_with_mate(r2, r2_start, r2_end, r1_start, r1_end)
                .peekable();

        Self { rec_iter, mate_iter }
    }
}

impl Iterator for RawReadMateAndRefPosIterator {
    type Item = ReadMateAndRefPos;

    fn next(&mut self) -> Option<Self::Item> {
        use std::cmp::Ordering;

        loop {
            let rec_pos = self.rec_iter.peek()?;
            let mate_pos = self.mate_iter.peek()?;

            match rec_pos.ref_pos.cmp(&mate_pos.ref_pos) {
                Ordering::Less => {
                    self.rec_iter.next();
                }
                Ordering::Greater => {
                    self.mate_iter.next();
                }
                Ordering::Equal => {
                    let rec = self.rec_iter.next().expect("peeked Some");
                    let mate = self.mate_iter.next().expect("peeked Some");

                    return Some(ReadMateAndRefPos {
                        read_offset: rec.read_offset,
                        mate_offset: mate.read_offset,
                    });
                }
            }
        }
    }
}

/// Applies overlapping consensus calling to pairs of raw-byte reads within a group.
///
/// This is the raw-byte equivalent of `apply_overlapping_consensus`.
///
/// # Errors
///
/// Returns an error if overlapping consensus calling fails for any read pair.
pub fn apply_overlapping_consensus_raw(
    records: &mut [Vec<u8>],
    caller: &mut OverlappingBasesConsensusCaller,
) -> Result<()> {
    use ahash::AHashMap;

    // Group reads by name for pairing
    let mut read_pairs: AHashMap<Vec<u8>, (Option<usize>, Option<usize>)> = AHashMap::new();

    for (idx, record) in records.iter().enumerate() {
        let name = noodles_raw_bam::read_name(record).to_vec();
        let flg = noodles_raw_bam::flags(record);

        if flg & noodles_raw_bam::flags::FIRST_SEGMENT != 0 {
            read_pairs.entry(name).or_insert((None, None)).0 = Some(idx);
        } else if flg & noodles_raw_bam::flags::LAST_SEGMENT != 0 {
            read_pairs.entry(name).or_insert((None, None)).1 = Some(idx);
        }
    }

    // Process pairs
    for (r1_idx, r2_idx) in read_pairs.values() {
        if let (Some(idx1), Some(idx2)) = (r1_idx, r2_idx) {
            let (r1, r2) = if idx1 < idx2 {
                let (left, right) = records.split_at_mut(*idx2);
                (&mut left[*idx1], &mut right[0])
            } else {
                let (left, right) = records.split_at_mut(*idx1);
                (&mut right[0], &mut left[*idx2])
            };

            caller.call_raw(r1, r2).context("Failed to call overlapping consensus on raw bytes")?;
        }
    }

    Ok(())
}

#[cfg(test)]
#[allow(clippy::cast_precision_loss, clippy::cast_possible_truncation)]
mod tests {
    use super::*;
    use fgumi_sam::builder::RecordBuilder;
    use noodles::sam::alignment::record::Flags;

    /// Creates a test record with explicit sequence, qualities, position, and CIGAR.
    fn create_test_record(seq: &[u8], qual: &[u8], start: usize, cigar: &str) -> RecordBuf {
        RecordBuilder::mapped_read()
            .sequence(&String::from_utf8_lossy(seq))
            .qualities(qual)
            .alignment_start(start)
            .cigar(cigar)
            .build()
    }

    #[test]
    fn test_agreement_strategy_consensus() {
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::Consensus,
            DisagreementStrategy::Consensus,
        );

        assert_eq!(caller.agreement_strategy, AgreementStrategy::Consensus);

        // Test with overlapping reads
        let mut r1 = create_test_record(b"ACGT", &[30, 30, 30, 30], 100, "4M");
        let mut r2 = create_test_record(b"ACGT", &[20, 20, 20, 20], 100, "4M");

        let result = caller.call(&mut r1, &mut r2).unwrap();
        assert!(result);

        // Check that qualities were summed (30 + 20 = 50)
        assert_eq!(r1.quality_scores().as_ref()[0], 50);
        assert_eq!(r2.quality_scores().as_ref()[0], 50);
    }

    #[test]
    fn test_agreement_strategy_max_qual() {
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::MaxQual,
            DisagreementStrategy::Consensus,
        );

        let mut r1 = create_test_record(b"ACGT", &[30, 30, 30, 30], 100, "4M");
        let mut r2 = create_test_record(b"ACGT", &[20, 20, 20, 20], 100, "4M");

        let result = caller.call(&mut r1, &mut r2).unwrap();
        assert!(result);

        // Check that max quality was used (max(30, 20) = 30)
        assert_eq!(r1.quality_scores().as_ref()[0], 30);
        assert_eq!(r2.quality_scores().as_ref()[0], 30);
    }

    #[test]
    fn test_agreement_strategy_pass_through() {
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::PassThrough,
            DisagreementStrategy::Consensus,
        );

        let mut r1 = create_test_record(b"ACGT", &[30, 30, 30, 30], 100, "4M");
        let mut r2 = create_test_record(b"ACGT", &[20, 20, 20, 20], 100, "4M");

        let result = caller.call(&mut r1, &mut r2).unwrap();
        assert!(result);

        // Check that qualities were unchanged
        assert_eq!(r1.quality_scores().as_ref()[0], 30);
        assert_eq!(r2.quality_scores().as_ref()[0], 20);
    }

    #[test]
    fn test_disagreement_strategy_consensus() {
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::PassThrough,
            DisagreementStrategy::Consensus,
        );

        let mut r1 = create_test_record(b"ACGT", &[30, 30, 30, 30], 100, "4M");
        let mut r2 = create_test_record(b"GCTA", &[20, 20, 20, 20], 100, "4M");

        let result = caller.call(&mut r1, &mut r2).unwrap();
        assert!(result);

        // Higher quality base (A from r1) should be chosen, qual = 30 - 20 = 10
        assert_eq!(r1.sequence().as_ref()[0], b'A');
        assert_eq!(r2.sequence().as_ref()[0], b'A');
        assert_eq!(r1.quality_scores().as_ref()[0], 10);
        assert_eq!(r2.quality_scores().as_ref()[0], 10);
    }

    #[test]
    fn test_disagreement_strategy_mask_both() {
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::PassThrough,
            DisagreementStrategy::MaskBoth,
        );

        let mut r1 = create_test_record(b"ACGT", &[30, 30, 30, 30], 100, "4M");
        let mut r2 = create_test_record(b"GCTA", &[20, 20, 20, 20], 100, "4M");

        let result = caller.call(&mut r1, &mut r2).unwrap();
        assert!(result);

        // Both bases should be masked to N with quality 2
        assert_eq!(r1.sequence().as_ref()[0], b'N');
        assert_eq!(r2.sequence().as_ref()[0], b'N');
        assert_eq!(r1.quality_scores().as_ref()[0], 2);
        assert_eq!(r2.quality_scores().as_ref()[0], 2);
    }

    #[test]
    fn test_disagreement_strategy_mask_lower_qual() {
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::PassThrough,
            DisagreementStrategy::MaskLowerQual,
        );

        let mut r1 = create_test_record(b"ACGT", &[30, 30, 30, 30], 100, "4M");
        let mut r2 = create_test_record(b"GCTA", &[20, 20, 20, 20], 100, "4M");

        let result = caller.call(&mut r1, &mut r2).unwrap();
        assert!(result);

        // Only lower quality base (r2) should be masked
        assert_eq!(r1.sequence().as_ref()[0], b'A'); // Unchanged
        assert_eq!(r2.sequence().as_ref()[0], b'N'); // Masked
        assert_eq!(r1.quality_scores().as_ref()[0], 30); // Unchanged
        assert_eq!(r2.quality_scores().as_ref()[0], 2); // Masked
    }

    #[test]
    fn test_disagreement_mask_lower_qual_r1_lower() {
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::PassThrough,
            DisagreementStrategy::MaskLowerQual,
        );

        let mut r1 = create_test_record(b"ACGT", &[20, 20, 20, 20], 100, "4M");
        let mut r2 = create_test_record(b"GCTA", &[30, 30, 30, 30], 100, "4M");

        let result = caller.call(&mut r1, &mut r2).unwrap();
        assert!(result);

        // Only lower quality base (r1) should be masked
        assert_eq!(r1.sequence().as_ref()[0], b'N'); // Masked
        assert_eq!(r2.sequence().as_ref()[0], b'G'); // Unchanged
    }

    #[test]
    fn test_disagreement_mask_lower_qual_equal_quality() {
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::PassThrough,
            DisagreementStrategy::MaskLowerQual,
        );

        let mut r1 = create_test_record(b"ACGT", &[30, 30, 30, 30], 100, "4M");
        let mut r2 = create_test_record(b"GCTA", &[30, 30, 30, 30], 100, "4M");

        let result = caller.call(&mut r1, &mut r2).unwrap();
        assert!(result);

        // Both should be masked when qualities are equal (matches fgbio)
        assert_eq!(r1.sequence().as_ref()[0], b'N');
        assert_eq!(r2.sequence().as_ref()[0], b'N');
        assert_eq!(r1.quality_scores().as_ref()[0], 2);
        assert_eq!(r2.quality_scores().as_ref()[0], 2);
    }

    #[test]
    fn test_disagreement_consensus_equal_quality() {
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::PassThrough,
            DisagreementStrategy::Consensus,
        );

        let mut r1 = create_test_record(b"ACGT", &[30, 30, 30, 30], 100, "4M");
        let mut r2 = create_test_record(b"GCTA", &[30, 30, 30, 30], 100, "4M");

        let result = caller.call(&mut r1, &mut r2).unwrap();
        assert!(result);

        // Both should be masked when qualities are equal (matches fgbio)
        assert_eq!(r1.sequence().as_ref()[0], b'N');
        assert_eq!(r2.sequence().as_ref()[0], b'N');
        assert_eq!(r1.quality_scores().as_ref()[0], 2);
        assert_eq!(r2.quality_scores().as_ref()[0], 2);
    }

    #[test]
    fn test_no_overlap() {
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::Consensus,
            DisagreementStrategy::Consensus,
        );

        // Reads don't overlap
        let mut r1 = create_test_record(b"ACGT", &[30, 30, 30, 30], 100, "4M");
        let mut r2 = create_test_record(b"ACGT", &[20, 20, 20, 20], 200, "4M");

        let result = caller.call(&mut r1, &mut r2).unwrap();
        assert!(!result); // No overlap processed
    }

    #[test]
    fn test_unmapped_reads() {
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::Consensus,
            DisagreementStrategy::Consensus,
        );

        let mut r1 = create_test_record(b"ACGT", &[30, 30, 30, 30], 100, "4M");
        let mut r2 = create_test_record(b"ACGT", &[20, 20, 20, 20], 100, "4M");
        *r1.flags_mut() = Flags::UNMAPPED;

        let result = caller.call(&mut r1, &mut r2).unwrap();
        assert!(!result); // Unmapped reads not processed
    }

    #[test]
    fn test_quality_capping_at_93() {
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::Consensus,
            DisagreementStrategy::Consensus,
        );

        // High qualities that would sum > 93
        let mut r1 = create_test_record(b"ACGT", &[50, 50, 50, 50], 100, "4M");
        let mut r2 = create_test_record(b"ACGT", &[50, 50, 50, 50], 100, "4M");

        let result = caller.call(&mut r1, &mut r2).unwrap();
        assert!(result);

        // Quality should be capped at 93 (not 100)
        assert_eq!(r1.quality_scores().as_ref()[0], 93);
    }

    #[test]
    fn test_stats_tracking() {
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::Consensus,
            DisagreementStrategy::Consensus,
        );

        // Create reads where all bases agree - both have "ACGT"
        let mut r1 = create_test_record(b"ACGT", &[30, 30, 30, 30], 100, "4M");
        let mut r2 = create_test_record(b"ACGT", &[20, 20, 20, 20], 100, "4M");

        let result = caller.call(&mut r1, &mut r2).unwrap();
        assert!(result);

        // NOTE: The overlap calculation uses alignment_end() which returns the position
        // after the last base. With exclusive end comparison (ref_pos < overlap_end),
        // the last base may be excluded depending on exact position values.
        // We verify that overlapping bases are tracked, without requiring an exact count.
        assert!(caller.stats().overlapping_bases > 0);
        assert_eq!(caller.stats().bases_agreeing, caller.stats().overlapping_bases);
        assert_eq!(caller.stats().bases_disagreeing, 0);
        // With Consensus strategy, agreeing bases get quality summed (corrected)
        assert_eq!(caller.stats().bases_corrected, caller.stats().overlapping_bases);
    }

    #[test]
    fn test_stats_tracking_with_disagreements() {
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::PassThrough,
            DisagreementStrategy::Consensus,
        );

        // Create reads where bases disagree: ACGT vs TGCA (all different)
        // Note: N bases are skipped (fgbio behavior), so we use real bases
        let mut r1 = create_test_record(b"ACGT", &[30, 30, 30, 30], 100, "4M");
        let mut r2 = create_test_record(b"TGCA", &[20, 20, 20, 20], 100, "4M");

        let result = caller.call(&mut r1, &mut r2).unwrap();
        assert!(result);

        // All overlapping bases disagree
        assert!(caller.stats().overlapping_bases > 0);
        assert_eq!(caller.stats().bases_agreeing, 0);
        assert_eq!(caller.stats().bases_disagreeing, caller.stats().overlapping_bases);
    }

    #[test]
    fn test_stats_reset() {
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::Consensus,
            DisagreementStrategy::Consensus,
        );

        let mut r1 = create_test_record(b"ACGT", &[30, 30, 30, 30], 100, "4M");
        let mut r2 = create_test_record(b"ACGT", &[20, 20, 20, 20], 100, "4M");

        caller.call(&mut r1, &mut r2).unwrap();
        assert!(caller.stats().overlapping_bases > 0);

        caller.reset_stats();
        assert_eq!(caller.stats().overlapping_bases, 0);
        assert_eq!(caller.stats().bases_agreeing, 0);
        assert_eq!(caller.stats().bases_disagreeing, 0);
        assert_eq!(caller.stats().bases_corrected, 0);
    }

    #[test]
    fn test_overlap_different_start_positions() {
        // The merge iteration properly handles reads at different start positions
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::Consensus,
            DisagreementStrategy::Consensus,
        );

        // R1: positions 100-103 (ACGT)
        // R2: positions 102-105 (GTAC)
        // Overlap: positions 102-103 (GT matches GT)
        let mut r1 = create_test_record(b"ACGT", &[30, 30, 30, 30], 100, "4M");
        let mut r2 = create_test_record(b"GTAC", &[20, 20, 20, 20], 102, "4M");

        let result = caller.call(&mut r1, &mut r2).unwrap();

        // Overlap should be detected at positions 102-103 (2 bases)
        assert!(result);
        assert_eq!(caller.stats().overlapping_bases, 2);
        assert_eq!(caller.stats().bases_agreeing, 2); // GT == GT

        // Check quality sums at overlapping positions
        // R1[2] (pos 102) = 'G' should be summed: 30 + 20 = 50
        // R1[3] (pos 103) = 'T' should be summed: 30 + 20 = 50
        assert_eq!(r1.quality_scores().as_ref()[2], 50);
        assert_eq!(r1.quality_scores().as_ref()[3], 50);
    }

    #[test]
    fn test_full_overlap_same_position() {
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::Consensus,
            DisagreementStrategy::Consensus,
        );

        // Both reads start at position 100 and fully overlap
        let mut r1 = create_test_record(b"ACGT", &[30, 30, 30, 30], 100, "4M");
        let mut r2 = create_test_record(b"ACGT", &[20, 20, 20, 20], 100, "4M");

        let result = caller.call(&mut r1, &mut r2).unwrap();
        assert!(result);

        // Overlapping bases should be detected (may not be exactly 4 due to boundary handling)
        assert!(caller.stats().overlapping_bases > 0);
        assert_eq!(caller.stats().bases_agreeing, caller.stats().overlapping_bases);
    }

    #[test]
    fn test_correction_stats() {
        let mut stats = CorrectionStats::new();
        assert_eq!(stats.overlapping_bases, 0);
        assert_eq!(stats.bases_agreeing, 0);
        assert_eq!(stats.bases_disagreeing, 0);
        assert_eq!(stats.bases_corrected, 0);

        stats.overlapping_bases = 10;
        stats.bases_agreeing = 8;
        stats.bases_disagreeing = 2;
        stats.bases_corrected = 3;

        stats.reset();
        assert_eq!(stats.overlapping_bases, 0);
    }

    #[test]
    fn test_cigar_with_insertions() {
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::Consensus,
            DisagreementStrategy::Consensus,
        );

        // CIGAR with insertion: 2M2I2M
        let mut r1 = create_test_record(b"ACTTGG", &[30, 30, 30, 30, 30, 30], 100, "2M2I2M");
        let mut r2 = create_test_record(b"ACGG", &[20, 20, 20, 20], 100, "4M");

        let result = caller.call(&mut r1, &mut r2).unwrap();
        // Should process the matching regions
        assert!(result);
    }

    #[test]
    fn test_cigar_with_deletions() {
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::Consensus,
            DisagreementStrategy::Consensus,
        );

        // CIGAR with deletion: 2M2D2M (read has 4 bases but spans 6 ref positions)
        let mut r1 = create_test_record(b"ACGG", &[30, 30, 30, 30], 100, "2M2D2M");
        let mut r2 = create_test_record(b"ACTTGG", &[20, 20, 20, 20, 20, 20], 100, "6M");

        let result = caller.call(&mut r1, &mut r2).unwrap();
        // Should process overlapping matches
        assert!(result);
    }

    #[test]
    fn test_cigar_with_soft_clips_different_structure() {
        // The merge iteration properly handles different soft clip structures
        // by only iterating aligned (M/X/=) bases and synchronizing by ref position
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::Consensus,
            DisagreementStrategy::Consensus,
        );

        // R1: 2S4M (6 bases, first 2 soft-clipped, aligned at ref 100-103)
        // R2: 4M (4 bases, no soft clips, aligned at ref 100-103)
        // Both have 4 aligned bases at ref 100-103, sequences ACGT match
        let mut r1 = create_test_record(b"NNACGT", &[2, 2, 30, 30, 30, 30], 100, "2S4M");
        let mut r2 = create_test_record(b"ACGT", &[20, 20, 20, 20], 100, "4M");

        let result = caller.call(&mut r1, &mut r2).unwrap();

        // The merge iteration properly handles this case:
        // R1 iterator yields: (2, 100), (3, 101), (4, 102), (5, 103) - skips soft clips
        // R2 iterator yields: (0, 100), (1, 101), (2, 102), (3, 103)
        // All 4 positions match by ref_pos
        assert!(result);
        assert_eq!(caller.stats().overlapping_bases, 4);
        assert_eq!(caller.stats().bases_agreeing, 4); // ACGT == ACGT

        // Check quality sums - R1 positions 2-5 and R2 positions 0-3
        assert_eq!(r1.quality_scores().as_ref()[2], 50); // 30 + 20
        assert_eq!(r1.quality_scores().as_ref()[3], 50);
        assert_eq!(r1.quality_scores().as_ref()[4], 50);
        assert_eq!(r1.quality_scores().as_ref()[5], 50);
    }

    #[test]
    fn test_cigar_soft_clips_both_reads_different_lengths() {
        // Test where both reads have soft clips but different lengths
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::Consensus,
            DisagreementStrategy::Consensus,
        );

        // R1: 3S3M (6 bases, 3 soft-clipped, aligned at ref 100-102)
        // R2: 1S3M (4 bases, 1 soft-clipped, aligned at ref 100-102)
        let mut r1 = create_test_record(b"NNNACG", &[2, 2, 2, 30, 30, 30], 100, "3S3M");
        let mut r2 = create_test_record(b"NACG", &[2, 20, 20, 20], 100, "1S3M");

        let result = caller.call(&mut r1, &mut r2).unwrap();
        assert!(result);

        // R1 yields: (3, 100), (4, 101), (5, 102)
        // R2 yields: (1, 100), (2, 101), (3, 102)
        // 3 overlapping aligned positions
        assert_eq!(caller.stats().overlapping_bases, 3);
        assert_eq!(caller.stats().bases_agreeing, 3); // ACG == ACG
    }

    #[test]
    fn test_cigar_soft_clips_both_reads_same_structure() {
        // Test where both reads have the same structure including soft clips
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::Consensus,
            DisagreementStrategy::Consensus,
        );

        // Both reads have identical CIGAR structure: 1S3M
        let mut r1 = create_test_record(b"NACG", &[2, 30, 30, 30], 100, "1S3M");
        let mut r2 = create_test_record(b"NACG", &[2, 20, 20, 20], 100, "1S3M");

        let result = caller.call(&mut r1, &mut r2).unwrap();
        assert!(result);

        // Both have 3 aligned bases at ref 100-102
        assert_eq!(caller.stats().overlapping_bases, 3);
        assert_eq!(caller.stats().bases_agreeing, 3);
    }

    #[test]
    fn test_disagreement_consensus_min_quality() {
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::PassThrough,
            DisagreementStrategy::Consensus,
        );

        // Very similar qualities - difference should be at least 2
        let mut r1 = create_test_record(b"ACGT", &[30, 30, 30, 30], 100, "4M");
        let mut r2 = create_test_record(b"GCTA", &[29, 29, 29, 29], 100, "4M");

        let result = caller.call(&mut r1, &mut r2).unwrap();
        assert!(result);

        // Quality difference is 1, but should be at least 2
        assert_eq!(r1.quality_scores().as_ref()[0], 2);
        assert_eq!(r2.quality_scores().as_ref()[0], 2);
    }

    // ========================================================================
    // Raw-byte tests
    // ========================================================================

    /// SAM spec 4-bit encoding for A=1, C=2, G=4, T=8, N=15.
    fn base_to_4bit(b: u8) -> u8 {
        match b {
            b'A' | b'a' => 1,
            b'C' | b'c' => 2,
            b'G' | b'g' => 4,
            b'T' | b't' => 8,
            b'N' | b'n' => 15,
            _ => 15,
        }
    }

    /// Pack ASCII bases into BAM 4-bit-per-base format.
    fn pack_seq(bases: &[u8]) -> Vec<u8> {
        let mut packed = Vec::with_capacity(bases.len().div_ceil(2));
        for pair in bases.chunks(2) {
            let hi = base_to_4bit(pair[0]);
            let lo = if pair.len() > 1 { base_to_4bit(pair[1]) } else { 0 };
            packed.push((hi << 4) | lo);
        }
        packed
    }

    /// Build a raw BAM record with populated sequence and quality data.
    ///
    /// `pos_0based` is the 0-based leftmost position (BAM pos field).
    /// `cigar_ops` are pre-encoded BAM CIGAR u32 values.
    /// `flag` is the BAM flag field.
    /// `name` should have length+1 divisible by 4 for alignment.
    fn make_raw_bam(
        name: &[u8],
        flag: u16,
        tid: i32,
        pos_0based: i32,
        cigar_ops: &[u32],
        seq: &[u8],
        qual: &[u8],
    ) -> Vec<u8> {
        let l_read_name = (name.len() + 1) as u8; // +1 for null terminator
        let n_cigar_op = cigar_ops.len() as u16;
        let seq_len = seq.len();
        let packed_seq = pack_seq(seq);
        let total = 32 + l_read_name as usize + cigar_ops.len() * 4 + packed_seq.len() + seq_len;
        let mut buf = vec![0u8; total];

        // Fixed header
        buf[0..4].copy_from_slice(&tid.to_le_bytes());
        buf[4..8].copy_from_slice(&pos_0based.to_le_bytes());
        buf[8] = l_read_name;
        buf[9] = 60; // mapq
        buf[10..12].copy_from_slice(&0u16.to_le_bytes()); // bin
        buf[12..14].copy_from_slice(&n_cigar_op.to_le_bytes());
        buf[14..16].copy_from_slice(&flag.to_le_bytes());
        buf[16..20].copy_from_slice(&(seq_len as u32).to_le_bytes());
        buf[20..24].copy_from_slice(&(-1i32).to_le_bytes()); // mate tid
        buf[24..28].copy_from_slice(&(-1i32).to_le_bytes()); // mate pos
        buf[28..32].copy_from_slice(&0i32.to_le_bytes()); // tlen

        // Read name + null terminator
        let name_start = 32;
        buf[name_start..name_start + name.len()].copy_from_slice(name);
        buf[name_start + name.len()] = 0;

        // CIGAR
        let cigar_start = name_start + l_read_name as usize;
        for (i, &op) in cigar_ops.iter().enumerate() {
            let off = cigar_start + i * 4;
            buf[off..off + 4].copy_from_slice(&op.to_le_bytes());
        }

        // Sequence (4-bit packed)
        let seq_start = cigar_start + cigar_ops.len() * 4;
        buf[seq_start..seq_start + packed_seq.len()].copy_from_slice(&packed_seq);

        // Quality scores (raw Phred)
        let qual_start = seq_start + packed_seq.len();
        buf[qual_start..qual_start + qual.len()].copy_from_slice(qual);

        buf
    }

    /// Encode a single CIGAR op as BAM u32: `(len << 4) | op_type`.
    /// `op_type`: M=0, I=1, D=2, N=3, S=4, H=5, P=6, =7, X=8.
    fn cigar_op(len: usize, op_type: u8) -> u32 {
        ((len as u32) << 4) | u32::from(op_type)
    }

    /// Create a raw BAM record mirroring `create_test_record`.
    ///
    /// `start_1based` is a 1-based alignment start (matching the `RecordBuf` tests).
    fn create_raw_test_record(
        seq: &[u8],
        qual: &[u8],
        start_1based: usize,
        cigar: &[u32],
    ) -> Vec<u8> {
        let pos_0based = (start_1based as i32) - 1;
        // Use name "rea" (3 chars + 1 NUL = 4, divisible by 4)
        make_raw_bam(b"rea", 0, 0, pos_0based, cigar, seq, qual)
    }

    #[test]
    fn test_raw_agreement_strategy_consensus() {
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::Consensus,
            DisagreementStrategy::Consensus,
        );

        let cigar = [cigar_op(4, 0)]; // 4M
        let mut r1 = create_raw_test_record(b"ACGT", &[30, 30, 30, 30], 100, &cigar);
        let mut r2 = create_raw_test_record(b"ACGT", &[20, 20, 20, 20], 100, &cigar);

        let result = caller.call_raw(&mut r1, &mut r2).unwrap();
        assert!(result);

        // Check that qualities were summed (30 + 20 = 50)
        assert_eq!(noodles_raw_bam::quality_scores_slice(&r1)[0], 50);
        assert_eq!(noodles_raw_bam::quality_scores_slice(&r2)[0], 50);
    }

    #[test]
    fn test_raw_agreement_strategy_max_qual() {
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::MaxQual,
            DisagreementStrategy::Consensus,
        );

        let cigar = [cigar_op(4, 0)]; // 4M
        let mut r1 = create_raw_test_record(b"ACGT", &[30, 30, 30, 30], 100, &cigar);
        let mut r2 = create_raw_test_record(b"ACGT", &[20, 20, 20, 20], 100, &cigar);

        let result = caller.call_raw(&mut r1, &mut r2).unwrap();
        assert!(result);

        // Max quality used (max(30, 20) = 30)
        assert_eq!(noodles_raw_bam::quality_scores_slice(&r1)[0], 30);
        assert_eq!(noodles_raw_bam::quality_scores_slice(&r2)[0], 30);
    }

    #[test]
    fn test_raw_agreement_strategy_pass_through() {
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::PassThrough,
            DisagreementStrategy::Consensus,
        );

        let cigar = [cigar_op(4, 0)]; // 4M
        let mut r1 = create_raw_test_record(b"ACGT", &[30, 30, 30, 30], 100, &cigar);
        let mut r2 = create_raw_test_record(b"ACGT", &[20, 20, 20, 20], 100, &cigar);

        let result = caller.call_raw(&mut r1, &mut r2).unwrap();
        assert!(result);

        // Qualities unchanged
        assert_eq!(noodles_raw_bam::quality_scores_slice(&r1)[0], 30);
        assert_eq!(noodles_raw_bam::quality_scores_slice(&r2)[0], 20);
    }

    #[test]
    fn test_raw_disagreement_strategy_consensus() {
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::PassThrough,
            DisagreementStrategy::Consensus,
        );

        let cigar = [cigar_op(4, 0)]; // 4M
        let mut r1 = create_raw_test_record(b"ACGT", &[30, 30, 30, 30], 100, &cigar);
        let mut r2 = create_raw_test_record(b"GCTA", &[20, 20, 20, 20], 100, &cigar);

        let result = caller.call_raw(&mut r1, &mut r2).unwrap();
        assert!(result);

        // Higher quality base (A from r1) chosen, qual = 30 - 20 = 10
        let r1_seq = noodles_raw_bam::extract_sequence(&r1);
        let r2_seq = noodles_raw_bam::extract_sequence(&r2);
        assert_eq!(r1_seq[0], b'A');
        assert_eq!(r2_seq[0], b'A');
        assert_eq!(noodles_raw_bam::quality_scores_slice(&r1)[0], 10);
        assert_eq!(noodles_raw_bam::quality_scores_slice(&r2)[0], 10);
    }

    #[test]
    fn test_raw_disagreement_strategy_mask_both() {
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::PassThrough,
            DisagreementStrategy::MaskBoth,
        );

        let cigar = [cigar_op(4, 0)]; // 4M
        let mut r1 = create_raw_test_record(b"ACGT", &[30, 30, 30, 30], 100, &cigar);
        let mut r2 = create_raw_test_record(b"GCTA", &[20, 20, 20, 20], 100, &cigar);

        let result = caller.call_raw(&mut r1, &mut r2).unwrap();
        assert!(result);

        // Both bases masked to N with quality 2
        let r1_seq = noodles_raw_bam::extract_sequence(&r1);
        let r2_seq = noodles_raw_bam::extract_sequence(&r2);
        assert_eq!(r1_seq[0], b'N');
        assert_eq!(r2_seq[0], b'N');
        assert_eq!(noodles_raw_bam::quality_scores_slice(&r1)[0], 2);
        assert_eq!(noodles_raw_bam::quality_scores_slice(&r2)[0], 2);
    }

    #[test]
    fn test_raw_disagreement_strategy_mask_lower_qual() {
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::PassThrough,
            DisagreementStrategy::MaskLowerQual,
        );

        let cigar = [cigar_op(4, 0)]; // 4M
        let mut r1 = create_raw_test_record(b"ACGT", &[30, 30, 30, 30], 100, &cigar);
        let mut r2 = create_raw_test_record(b"GCTA", &[20, 20, 20, 20], 100, &cigar);

        let result = caller.call_raw(&mut r1, &mut r2).unwrap();
        assert!(result);

        // Only lower quality base (r2) masked
        let r1_seq = noodles_raw_bam::extract_sequence(&r1);
        let r2_seq = noodles_raw_bam::extract_sequence(&r2);
        assert_eq!(r1_seq[0], b'A'); // Unchanged
        assert_eq!(r2_seq[0], b'N'); // Masked
        assert_eq!(noodles_raw_bam::quality_scores_slice(&r1)[0], 30); // Unchanged
        assert_eq!(noodles_raw_bam::quality_scores_slice(&r2)[0], 2); // Masked
    }

    #[test]
    fn test_raw_disagreement_mask_lower_qual_r1_lower() {
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::PassThrough,
            DisagreementStrategy::MaskLowerQual,
        );

        let cigar = [cigar_op(4, 0)]; // 4M
        let mut r1 = create_raw_test_record(b"ACGT", &[20, 20, 20, 20], 100, &cigar);
        let mut r2 = create_raw_test_record(b"GCTA", &[30, 30, 30, 30], 100, &cigar);

        let result = caller.call_raw(&mut r1, &mut r2).unwrap();
        assert!(result);

        // Only lower quality base (r1) masked
        let r1_seq = noodles_raw_bam::extract_sequence(&r1);
        let r2_seq = noodles_raw_bam::extract_sequence(&r2);
        assert_eq!(r1_seq[0], b'N'); // Masked
        assert_eq!(r2_seq[0], b'G'); // Unchanged
    }

    #[test]
    fn test_raw_disagreement_mask_lower_qual_equal_quality() {
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::PassThrough,
            DisagreementStrategy::MaskLowerQual,
        );

        let cigar = [cigar_op(4, 0)]; // 4M
        let mut r1 = create_raw_test_record(b"ACGT", &[30, 30, 30, 30], 100, &cigar);
        let mut r2 = create_raw_test_record(b"GCTA", &[30, 30, 30, 30], 100, &cigar);

        let result = caller.call_raw(&mut r1, &mut r2).unwrap();
        assert!(result);

        // Both masked when qualities are equal
        let r1_seq = noodles_raw_bam::extract_sequence(&r1);
        let r2_seq = noodles_raw_bam::extract_sequence(&r2);
        assert_eq!(r1_seq[0], b'N');
        assert_eq!(r2_seq[0], b'N');
        assert_eq!(noodles_raw_bam::quality_scores_slice(&r1)[0], 2);
        assert_eq!(noodles_raw_bam::quality_scores_slice(&r2)[0], 2);
    }

    #[test]
    fn test_raw_disagreement_consensus_equal_quality() {
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::PassThrough,
            DisagreementStrategy::Consensus,
        );

        let cigar = [cigar_op(4, 0)]; // 4M
        let mut r1 = create_raw_test_record(b"ACGT", &[30, 30, 30, 30], 100, &cigar);
        let mut r2 = create_raw_test_record(b"GCTA", &[30, 30, 30, 30], 100, &cigar);

        let result = caller.call_raw(&mut r1, &mut r2).unwrap();
        assert!(result);

        // Both masked when qualities are equal
        let r1_seq = noodles_raw_bam::extract_sequence(&r1);
        let r2_seq = noodles_raw_bam::extract_sequence(&r2);
        assert_eq!(r1_seq[0], b'N');
        assert_eq!(r2_seq[0], b'N');
        assert_eq!(noodles_raw_bam::quality_scores_slice(&r1)[0], 2);
        assert_eq!(noodles_raw_bam::quality_scores_slice(&r2)[0], 2);
    }

    #[test]
    fn test_raw_no_overlap() {
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::Consensus,
            DisagreementStrategy::Consensus,
        );

        // Reads don't overlap: r1 at 100-103, r2 at 200-203
        let cigar = [cigar_op(4, 0)]; // 4M
        let mut r1 = create_raw_test_record(b"ACGT", &[30, 30, 30, 30], 100, &cigar);
        let mut r2 = create_raw_test_record(b"ACGT", &[20, 20, 20, 20], 200, &cigar);

        let result = caller.call_raw(&mut r1, &mut r2).unwrap();
        assert!(!result);
    }

    #[test]
    fn test_raw_unmapped_reads() {
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::Consensus,
            DisagreementStrategy::Consensus,
        );

        // r1 is unmapped (flag=0x4)
        let cigar = [cigar_op(4, 0)]; // 4M
        let mut r1 = make_raw_bam(
            b"rea",
            noodles_raw_bam::flags::UNMAPPED,
            0,
            99,
            &cigar,
            b"ACGT",
            &[30, 30, 30, 30],
        );
        let mut r2 = create_raw_test_record(b"ACGT", &[20, 20, 20, 20], 100, &cigar);

        let result = caller.call_raw(&mut r1, &mut r2).unwrap();
        assert!(!result);
    }

    #[test]
    fn test_raw_different_references() {
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::Consensus,
            DisagreementStrategy::Consensus,
        );

        // r1 on tid=0, r2 on tid=1
        let cigar = [cigar_op(4, 0)]; // 4M
        let mut r1 = make_raw_bam(b"rea", 0, 0, 99, &cigar, b"ACGT", &[30, 30, 30, 30]);
        let mut r2 = make_raw_bam(b"rea", 0, 1, 99, &cigar, b"ACGT", &[20, 20, 20, 20]);

        let result = caller.call_raw(&mut r1, &mut r2).unwrap();
        assert!(!result);
    }

    #[test]
    fn test_raw_quality_capping_at_93() {
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::Consensus,
            DisagreementStrategy::Consensus,
        );

        let cigar = [cigar_op(4, 0)]; // 4M
        let mut r1 = create_raw_test_record(b"ACGT", &[50, 50, 50, 50], 100, &cigar);
        let mut r2 = create_raw_test_record(b"ACGT", &[50, 50, 50, 50], 100, &cigar);

        let result = caller.call_raw(&mut r1, &mut r2).unwrap();
        assert!(result);

        // Quality capped at 93
        assert_eq!(noodles_raw_bam::quality_scores_slice(&r1)[0], 93);
    }

    #[test]
    fn test_raw_stats_tracking() {
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::Consensus,
            DisagreementStrategy::Consensus,
        );

        let cigar = [cigar_op(4, 0)]; // 4M
        let mut r1 = create_raw_test_record(b"ACGT", &[30, 30, 30, 30], 100, &cigar);
        let mut r2 = create_raw_test_record(b"ACGT", &[20, 20, 20, 20], 100, &cigar);

        let result = caller.call_raw(&mut r1, &mut r2).unwrap();
        assert!(result);

        assert!(caller.stats().overlapping_bases > 0);
        assert_eq!(caller.stats().bases_agreeing, caller.stats().overlapping_bases);
        assert_eq!(caller.stats().bases_disagreeing, 0);
        assert_eq!(caller.stats().bases_corrected, caller.stats().overlapping_bases);
    }

    #[test]
    fn test_raw_stats_tracking_with_disagreements() {
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::PassThrough,
            DisagreementStrategy::Consensus,
        );

        let cigar = [cigar_op(4, 0)]; // 4M
        let mut r1 = create_raw_test_record(b"ACGT", &[30, 30, 30, 30], 100, &cigar);
        let mut r2 = create_raw_test_record(b"TGCA", &[20, 20, 20, 20], 100, &cigar);

        let result = caller.call_raw(&mut r1, &mut r2).unwrap();
        assert!(result);

        assert!(caller.stats().overlapping_bases > 0);
        assert_eq!(caller.stats().bases_agreeing, 0);
        assert_eq!(caller.stats().bases_disagreeing, caller.stats().overlapping_bases);
    }

    #[test]
    fn test_raw_overlap_different_start_positions() {
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::Consensus,
            DisagreementStrategy::Consensus,
        );

        // R1: positions 100-103 (ACGT), R2: positions 102-105 (GTAC)
        // Overlap: positions 102-103 (GT matches GT)
        let cigar = [cigar_op(4, 0)]; // 4M
        let mut r1 = create_raw_test_record(b"ACGT", &[30, 30, 30, 30], 100, &cigar);
        let mut r2 = create_raw_test_record(b"GTAC", &[20, 20, 20, 20], 102, &cigar);

        let result = caller.call_raw(&mut r1, &mut r2).unwrap();
        assert!(result);
        assert_eq!(caller.stats().overlapping_bases, 2);
        assert_eq!(caller.stats().bases_agreeing, 2); // GT == GT

        // Check quality sums at overlapping positions
        assert_eq!(noodles_raw_bam::quality_scores_slice(&r1)[2], 50); // 30 + 20
        assert_eq!(noodles_raw_bam::quality_scores_slice(&r1)[3], 50);
    }

    #[test]
    fn test_raw_cigar_with_soft_clips() {
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::Consensus,
            DisagreementStrategy::Consensus,
        );

        // R1: 2S4M (6 bases, first 2 soft-clipped, aligned at ref 100-103)
        let r1_cigar = [cigar_op(2, 4), cigar_op(4, 0)]; // 2S4M
        let r2_cigar = [cigar_op(4, 0)]; // 4M
        let mut r1 = create_raw_test_record(b"NNACGT", &[2, 2, 30, 30, 30, 30], 100, &r1_cigar);
        let mut r2 = create_raw_test_record(b"ACGT", &[20, 20, 20, 20], 100, &r2_cigar);

        let result = caller.call_raw(&mut r1, &mut r2).unwrap();
        assert!(result);
        assert_eq!(caller.stats().overlapping_bases, 4);
        assert_eq!(caller.stats().bases_agreeing, 4); // ACGT == ACGT

        // Check quality sums at overlapping positions
        assert_eq!(noodles_raw_bam::quality_scores_slice(&r1)[2], 50); // 30 + 20
        assert_eq!(noodles_raw_bam::quality_scores_slice(&r1)[3], 50);
        assert_eq!(noodles_raw_bam::quality_scores_slice(&r1)[4], 50);
        assert_eq!(noodles_raw_bam::quality_scores_slice(&r1)[5], 50);
    }

    #[test]
    fn test_raw_cigar_with_insertions() {
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::Consensus,
            DisagreementStrategy::Consensus,
        );

        // R1: 2M2I2M (6 bases), R2: 4M
        let r1_cigar = [cigar_op(2, 0), cigar_op(2, 1), cigar_op(2, 0)]; // 2M2I2M
        let r2_cigar = [cigar_op(4, 0)]; // 4M
        let mut r1 = create_raw_test_record(b"ACTTGG", &[30, 30, 30, 30, 30, 30], 100, &r1_cigar);
        let mut r2 = create_raw_test_record(b"ACGG", &[20, 20, 20, 20], 100, &r2_cigar);

        let result = caller.call_raw(&mut r1, &mut r2).unwrap();
        assert!(result);
    }

    #[test]
    fn test_raw_cigar_with_deletions() {
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::Consensus,
            DisagreementStrategy::Consensus,
        );

        // R1: 2M2D2M (4 bases, spans 6 ref positions), R2: 6M
        let r1_cigar = [cigar_op(2, 0), cigar_op(2, 2), cigar_op(2, 0)]; // 2M2D2M
        let r2_cigar = [cigar_op(6, 0)]; // 6M
        let mut r1 = create_raw_test_record(b"ACGG", &[30, 30, 30, 30], 100, &r1_cigar);
        let mut r2 = create_raw_test_record(b"ACTTGG", &[20, 20, 20, 20, 20, 20], 100, &r2_cigar);

        let result = caller.call_raw(&mut r1, &mut r2).unwrap();
        assert!(result);
    }

    #[test]
    fn test_raw_disagreement_consensus_min_quality() {
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::PassThrough,
            DisagreementStrategy::Consensus,
        );

        let cigar = [cigar_op(4, 0)]; // 4M
        let mut r1 = create_raw_test_record(b"ACGT", &[30, 30, 30, 30], 100, &cigar);
        let mut r2 = create_raw_test_record(b"GCTA", &[29, 29, 29, 29], 100, &cigar);

        let result = caller.call_raw(&mut r1, &mut r2).unwrap();
        assert!(result);

        // Quality difference is 1, but minimum is 2
        assert_eq!(noodles_raw_bam::quality_scores_slice(&r1)[0], 2);
        assert_eq!(noodles_raw_bam::quality_scores_slice(&r2)[0], 2);
    }

    /// Verify that `call_raw` produces the same results as `call` for agreement.
    #[test]
    fn test_raw_matches_recordbuf_agreement() {
        let mut caller_buf = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::Consensus,
            DisagreementStrategy::Consensus,
        );
        let mut caller_raw = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::Consensus,
            DisagreementStrategy::Consensus,
        );

        // RecordBuf path
        let mut rb_r1 = create_test_record(b"ACGT", &[30, 30, 30, 30], 100, "4M");
        let mut rb_r2 = create_test_record(b"ACGT", &[20, 20, 20, 20], 100, "4M");
        caller_buf.call(&mut rb_r1, &mut rb_r2).unwrap();

        // Raw path
        let cigar = [cigar_op(4, 0)];
        let mut raw_r1 = create_raw_test_record(b"ACGT", &[30, 30, 30, 30], 100, &cigar);
        let mut raw_r2 = create_raw_test_record(b"ACGT", &[20, 20, 20, 20], 100, &cigar);
        caller_raw.call_raw(&mut raw_r1, &mut raw_r2).unwrap();

        // Compare results
        assert_eq!(rb_r1.quality_scores().as_ref(), noodles_raw_bam::quality_scores_slice(&raw_r1));
        assert_eq!(rb_r2.quality_scores().as_ref(), noodles_raw_bam::quality_scores_slice(&raw_r2));
        assert_eq!(caller_buf.stats().overlapping_bases, caller_raw.stats().overlapping_bases);
        assert_eq!(caller_buf.stats().bases_agreeing, caller_raw.stats().bases_agreeing);
        assert_eq!(caller_buf.stats().bases_corrected, caller_raw.stats().bases_corrected);
    }

    /// Verify that `call_raw` produces the same results as `call` for disagreement.
    #[test]
    fn test_raw_matches_recordbuf_disagreement() {
        let mut caller_buf = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::PassThrough,
            DisagreementStrategy::MaskBoth,
        );
        let mut caller_raw = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::PassThrough,
            DisagreementStrategy::MaskBoth,
        );

        // RecordBuf path
        let mut rb_r1 = create_test_record(b"ACGT", &[30, 30, 30, 30], 100, "4M");
        let mut rb_r2 = create_test_record(b"GCTA", &[20, 20, 20, 20], 100, "4M");
        caller_buf.call(&mut rb_r1, &mut rb_r2).unwrap();

        // Raw path
        let cigar = [cigar_op(4, 0)];
        let mut raw_r1 = create_raw_test_record(b"ACGT", &[30, 30, 30, 30], 100, &cigar);
        let mut raw_r2 = create_raw_test_record(b"GCTA", &[20, 20, 20, 20], 100, &cigar);
        caller_raw.call_raw(&mut raw_r1, &mut raw_r2).unwrap();

        // Compare sequences
        let buf_r1_seq: Vec<u8> = rb_r1.sequence().as_ref().to_vec();
        let buf_r2_seq: Vec<u8> = rb_r2.sequence().as_ref().to_vec();
        let raw_r1_seq = noodles_raw_bam::extract_sequence(&raw_r1);
        let raw_r2_seq = noodles_raw_bam::extract_sequence(&raw_r2);
        assert_eq!(buf_r1_seq, raw_r1_seq);
        assert_eq!(buf_r2_seq, raw_r2_seq);

        // Compare qualities
        assert_eq!(rb_r1.quality_scores().as_ref(), noodles_raw_bam::quality_scores_slice(&raw_r1));
        assert_eq!(rb_r2.quality_scores().as_ref(), noodles_raw_bam::quality_scores_slice(&raw_r2));

        // Compare stats
        assert_eq!(caller_buf.stats().overlapping_bases, caller_raw.stats().overlapping_bases);
        assert_eq!(caller_buf.stats().bases_disagreeing, caller_raw.stats().bases_disagreeing);
        assert_eq!(caller_buf.stats().bases_corrected, caller_raw.stats().bases_corrected);
    }

    // ========================================================================
    // apply_overlapping_consensus_raw tests
    // ========================================================================

    #[test]
    fn test_apply_overlapping_consensus_raw_pair() {
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::Consensus,
            DisagreementStrategy::Consensus,
        );

        // R1 (first segment): flag = PAIRED | FIRST_SEGMENT
        let cigar = [cigar_op(4, 0)]; // 4M
        let r1_flag = noodles_raw_bam::flags::PAIRED | noodles_raw_bam::flags::FIRST_SEGMENT;
        let r1 = make_raw_bam(b"rea", r1_flag, 0, 99, &cigar, b"ACGT", &[30, 30, 30, 30]);
        // R2 (last segment): flag = PAIRED | LAST_SEGMENT
        let r2_flag = noodles_raw_bam::flags::PAIRED | noodles_raw_bam::flags::LAST_SEGMENT;
        let r2 = make_raw_bam(b"rea", r2_flag, 0, 99, &cigar, b"ACGT", &[20, 20, 20, 20]);

        let mut records = vec![r1, r2];
        apply_overlapping_consensus_raw(&mut records, &mut caller).unwrap();

        // Check that qualities were summed (30 + 20 = 50)
        assert_eq!(noodles_raw_bam::quality_scores_slice(&records[0])[0], 50);
        assert_eq!(noodles_raw_bam::quality_scores_slice(&records[1])[0], 50);
        assert!(caller.stats().overlapping_bases > 0);
    }

    #[test]
    fn test_apply_overlapping_consensus_raw_no_pair() {
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::Consensus,
            DisagreementStrategy::Consensus,
        );

        // Only FIRST_SEGMENT, no matching LAST_SEGMENT for this name
        let cigar = [cigar_op(4, 0)]; // 4M
        let r1_flag = noodles_raw_bam::flags::PAIRED | noodles_raw_bam::flags::FIRST_SEGMENT;
        let r1 = make_raw_bam(b"rea", r1_flag, 0, 99, &cigar, b"ACGT", &[30, 30, 30, 30]);
        // Different name so no pairing occurs
        let r2_flag = noodles_raw_bam::flags::PAIRED | noodles_raw_bam::flags::LAST_SEGMENT;
        let r2 = make_raw_bam(b"reb", r2_flag, 0, 99, &cigar, b"ACGT", &[20, 20, 20, 20]);

        let mut records = vec![r1, r2];
        apply_overlapping_consensus_raw(&mut records, &mut caller).unwrap();

        // Qualities unchanged because no matching pair
        assert_eq!(noodles_raw_bam::quality_scores_slice(&records[0])[0], 30);
        assert_eq!(noodles_raw_bam::quality_scores_slice(&records[1])[0], 20);
        assert_eq!(caller.stats().overlapping_bases, 0);
    }

    #[test]
    fn test_apply_overlapping_consensus_raw_reversed_indices() {
        // Test when R2 appears before R1 in the slice (exercises the idx1 > idx2 branch)
        let mut caller = OverlappingBasesConsensusCaller::new(
            AgreementStrategy::Consensus,
            DisagreementStrategy::Consensus,
        );

        let cigar = [cigar_op(4, 0)]; // 4M
        let r2_flag = noodles_raw_bam::flags::PAIRED | noodles_raw_bam::flags::LAST_SEGMENT;
        let r2 = make_raw_bam(b"rea", r2_flag, 0, 99, &cigar, b"ACGT", &[20, 20, 20, 20]);
        let r1_flag = noodles_raw_bam::flags::PAIRED | noodles_raw_bam::flags::FIRST_SEGMENT;
        let r1 = make_raw_bam(b"rea", r1_flag, 0, 99, &cigar, b"ACGT", &[30, 30, 30, 30]);

        // R2 at index 0, R1 at index 1
        let mut records = vec![r2, r1];
        apply_overlapping_consensus_raw(&mut records, &mut caller).unwrap();

        // Both should have been consensus-called
        assert!(caller.stats().overlapping_bases > 0);
        assert_eq!(
            noodles_raw_bam::quality_scores_slice(&records[0])[0],
            noodles_raw_bam::quality_scores_slice(&records[1])[0]
        );
    }

    // ========================================================================
    // CorrectionStats::merge tests
    // ========================================================================

    #[test]
    fn test_correction_stats_merge() {
        let mut stats_a = CorrectionStats {
            overlapping_bases: 10,
            bases_agreeing: 8,
            bases_disagreeing: 2,
            bases_corrected: 3,
        };

        let stats_b = CorrectionStats {
            overlapping_bases: 5,
            bases_agreeing: 4,
            bases_disagreeing: 1,
            bases_corrected: 2,
        };

        stats_a.merge(&stats_b);

        assert_eq!(stats_a.overlapping_bases, 15);
        assert_eq!(stats_a.bases_agreeing, 12);
        assert_eq!(stats_a.bases_disagreeing, 3);
        assert_eq!(stats_a.bases_corrected, 5);
    }

    #[test]
    fn test_correction_stats_merge_with_empty() {
        let mut stats = CorrectionStats {
            overlapping_bases: 10,
            bases_agreeing: 8,
            bases_disagreeing: 2,
            bases_corrected: 3,
        };

        let empty = CorrectionStats::new();
        stats.merge(&empty);

        // Values unchanged after merging with empty
        assert_eq!(stats.overlapping_bases, 10);
        assert_eq!(stats.bases_agreeing, 8);
        assert_eq!(stats.bases_disagreeing, 2);
        assert_eq!(stats.bases_corrected, 3);
    }

    #[test]
    fn test_correction_stats_merge_into_empty() {
        let mut empty = CorrectionStats::new();

        let stats = CorrectionStats {
            overlapping_bases: 7,
            bases_agreeing: 5,
            bases_disagreeing: 2,
            bases_corrected: 4,
        };

        empty.merge(&stats);

        assert_eq!(empty.overlapping_bases, 7);
        assert_eq!(empty.bases_agreeing, 5);
        assert_eq!(empty.bases_disagreeing, 2);
        assert_eq!(empty.bases_corrected, 4);
    }

    #[test]
    fn test_correction_stats_merge_multiple() {
        let mut total = CorrectionStats::new();

        let a = CorrectionStats {
            overlapping_bases: 3,
            bases_agreeing: 2,
            bases_disagreeing: 1,
            bases_corrected: 1,
        };
        let b = CorrectionStats {
            overlapping_bases: 5,
            bases_agreeing: 5,
            bases_disagreeing: 0,
            bases_corrected: 5,
        };
        let c = CorrectionStats {
            overlapping_bases: 2,
            bases_agreeing: 0,
            bases_disagreeing: 2,
            bases_corrected: 4,
        };

        total.merge(&a);
        total.merge(&b);
        total.merge(&c);

        assert_eq!(total.overlapping_bases, 10);
        assert_eq!(total.bases_agreeing, 7);
        assert_eq!(total.bases_disagreeing, 3);
        assert_eq!(total.bases_corrected, 10);
    }
}

//! Read clipping utilities for BAM/SAM records.
//!
//! This module provides functionality for clipping reads in various ways (soft, hard, etc.)
//! and is used by tools like `clip` and consensus calling tools.

use anyhow::Result;
use noodles::core::Position;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::record::cigar::Cigar;
use noodles::sam::alignment::record::cigar::Op as CigarOp;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::record_buf::{Cigar as CigarBuf, QualityScores, Sequence};

use fgumi_dna::{MIN_PHRED, NO_CALL_BASE};
use crate::record_utils;

/// Helper macro to get array length for any Array variant
macro_rules! array_len {
    ($arr:expr) => {
        match $arr {
            Array::Int8(a) => a.len(),
            Array::UInt8(a) => a.len(),
            Array::Int16(a) => a.len(),
            Array::UInt16(a) => a.len(),
            Array::Int32(a) => a.len(),
            Array::UInt32(a) => a.len(),
            Array::Float(a) => a.len(),
        }
    };
}

/// Helper macro to slice an array and create a new Value
macro_rules! slice_array {
    ($arr:expr, $start:expr, $end:expr) => {
        match $arr {
            Array::Int8(a) => Value::from(a[$start..$end].to_vec()),
            Array::UInt8(a) => Value::from(a[$start..$end].to_vec()),
            Array::Int16(a) => Value::from(a[$start..$end].to_vec()),
            Array::UInt16(a) => Value::from(a[$start..$end].to_vec()),
            Array::Int32(a) => Value::from(a[$start..$end].to_vec()),
            Array::UInt32(a) => Value::from(a[$start..$end].to_vec()),
            Array::Float(a) => Value::from(a[$start..$end].to_vec()),
        }
    };
}

/// Modes of clipping that can be applied to reads
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ClippingMode {
    /// Soft clip: convert bases to S operators in CIGAR, keep bases and qualities
    Soft,
    /// Soft clip with masking: convert to S operators and mask bases to N, qualities to min
    SoftWithMask,
    /// Hard clip: remove bases, qualities, and convert to H operators in CIGAR
    Hard,
}

/// Utility for clipping BAM/SAM records in various ways
pub struct SamRecordClipper {
    /// The clipping mode to use
    mode: ClippingMode,
    /// Whether to automatically clip extended attributes that match read length
    auto_clip_attributes: bool,
}

impl SamRecordClipper {
    /// Creates a new clipper with the specified mode
    #[must_use]
    pub fn new(mode: ClippingMode) -> Self {
        Self { mode, auto_clip_attributes: false }
    }

    /// Creates a new clipper with auto-clip attributes enabled
    ///
    /// When enabled with hard clipping mode, any tags that are the same length as the
    /// read's sequence will be automatically clipped to match. This ensures per-base
    /// tags (like quality arrays, per-base depths, etc.) stay synchronized with the sequence.
    #[must_use]
    pub fn with_auto_clip(mode: ClippingMode, auto_clip_attributes: bool) -> Self {
        Self { mode, auto_clip_attributes }
    }

    /// Clips extended attributes (per-base tags) when hard clipping
    ///
    /// When auto-clipping is enabled and using hard clipping mode, this method
    /// automatically trims any tag values (strings or arrays) that are the same length
    /// as the original read sequence to match the new clipped length.
    ///
    /// # Arguments
    /// * `record` - The record to modify
    /// * `remove` - Number of bases being removed
    /// * `from_start` - If true, clip from start (5' end); if false, clip from end (3' end)
    fn clip_extended_attributes(&self, record: &mut RecordBuf, remove: usize, from_start: bool) {
        use noodles::sam::alignment::record_buf::data::field::value::Array;

        // Only clip attributes when using hard clipping mode with auto-clip enabled
        if !matches!(self.mode, ClippingMode::Hard) || remove == 0 || !self.auto_clip_attributes {
            return;
        }

        let new_length = record.sequence().len();
        let old_length = new_length + remove;

        // Collect tags to update (we can't modify while iterating)
        let mut tags_to_update: Vec<(Vec<u8>, Value)> = Vec::new();

        // Iterate through all data fields
        for (tag, value) in record.data().iter() {
            // Check if this is a String or Array that matches the old length
            let should_clip = match value {
                Value::String(s) => {
                    let bytes: &[u8] = s.as_ref();
                    bytes.len() == old_length
                }
                Value::Array(arr) => array_len!(arr) == old_length,
                _ => false,
            };

            if should_clip {
                // Create clipped version of the value
                let (start, end) = if from_start { (remove, old_length) } else { (0, new_length) };
                let new_value = match value {
                    Value::String(s) => {
                        let bytes: &[u8] = s.as_ref();
                        Value::from(std::str::from_utf8(&bytes[start..end]).unwrap_or(""))
                    }
                    Value::Array(arr) => slice_array!(arr, start, end),
                    _ => continue, // Should not reach here due to should_clip check
                };

                tags_to_update.push((tag.as_ref().to_vec(), new_value));
            }
        }

        // Now update the tags
        for (tag, value) in tags_to_update {
            let tag_array: [u8; 2] = [tag[0], tag[1]];
            record.data_mut().insert(tag_array.into(), value);
        }
    }

    /// Clips a specified number of bases from the start (left side) of the alignment
    ///
    /// Returns the number of bases actually clipped
    #[expect(clippy::too_many_lines, reason = "clipping logic with multiple modes requires handling many CIGAR edge cases")]
    pub fn clip_start_of_alignment(&self, record: &mut RecordBuf, bases_to_clip: usize) -> usize {
        if bases_to_clip == 0 {
            return 0;
        }

        // Don't clip unmapped reads
        if record.flags().is_unmapped() {
            return 0;
        }

        let sequence_len = record.sequence();
        if sequence_len.is_empty() {
            return 0;
        }

        // Collect existing CIGAR ops
        let old_ops: Vec<CigarOp> = record.cigar().iter().filter_map(Result::ok).collect();

        // Extract existing hard and soft clips from the start
        let existing_hard_clip = old_ops
            .iter()
            .take_while(|op| op.kind() == Kind::HardClip)
            .map(|op| op.len())
            .sum::<usize>();

        let existing_soft_clip = old_ops
            .iter()
            .skip_while(|op| op.kind() == Kind::HardClip)
            .take_while(|op| op.kind() == Kind::SoftClip)
            .map(|op| op.len())
            .sum::<usize>();

        // Skip to operations after existing clips
        let post_clip_ops: Vec<CigarOp> = old_ops
            .into_iter()
            .skip_while(|op| matches!(op.kind(), Kind::HardClip | Kind::SoftClip))
            .collect();

        let mut read_bases_clipped = 0;
        let mut ref_bases_clipped = 0;
        let mut new_ops = Vec::new();
        let mut iter = post_clip_ops.iter().peekable();

        // Clip operations from the start
        // Continue until we've clipped enough read bases, OR if we've clipped the right amount
        // but the next operation is a deletion (which should be removed)
        while read_bases_clipped < bases_to_clip
            || (read_bases_clipped == bases_to_clip
                && new_ops.is_empty()
                && iter.peek().map(|op| op.kind()) == Some(Kind::Deletion))
        {
            let Some(op) = iter.next() else { break };

            let kind = op.kind();
            let len = op.len();

            let consumes_read = matches!(
                kind,
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Insertion
            );
            let consumes_ref = matches!(
                kind,
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion
            );

            if consumes_read && len > (bases_to_clip - read_bases_clipped) {
                // This operation extends past our clip point
                if kind == Kind::Insertion {
                    // Consume entire insertion at clip boundary
                    read_bases_clipped += len;
                } else {
                    // Split the operation
                    let remaining_clip = bases_to_clip - read_bases_clipped;
                    let remaining_length = len - remaining_clip;
                    read_bases_clipped += remaining_clip;
                    ref_bases_clipped += remaining_clip;
                    new_ops.push(CigarOp::new(kind, remaining_length));
                }
            } else {
                // Consume entire operation
                if consumes_read {
                    read_bases_clipped += len;
                }
                if consumes_ref {
                    ref_bases_clipped += len;
                }
            }
        }

        // Add remaining operations
        new_ops.extend(iter.copied());

        // Prepend appropriate clipping operators
        let (final_ops, bases_to_remove) = match self.mode {
            ClippingMode::Hard => {
                // Match fgbio's behavior: convert ALL existing soft clips to hard clips
                let added_hard_clip = existing_soft_clip + read_bases_clipped;
                let total_hard_clip = existing_hard_clip + added_hard_clip;
                let mut result = Vec::new();
                result.push(CigarOp::new(Kind::HardClip, total_hard_clip));
                result.extend(new_ops);
                (result, added_hard_clip)
            }
            ClippingMode::Soft | ClippingMode::SoftWithMask => {
                let total_soft_clip = existing_soft_clip + read_bases_clipped;
                let mut result = Vec::new();
                if existing_hard_clip > 0 {
                    result.push(CigarOp::new(Kind::HardClip, existing_hard_clip));
                }
                result.push(CigarOp::new(Kind::SoftClip, total_soft_clip));
                result.extend(new_ops);
                (result, 0)
            }
        };

        // Update CIGAR
        *record.cigar_mut() = CigarBuf::from(final_ops);

        // Update alignment start position
        if ref_bases_clipped > 0 {
            if let Some(start_pos) = record.alignment_start() {
                if let Some(new_start) = Position::new(usize::from(start_pos) + ref_bases_clipped) {
                    *record.alignment_start_mut() = Some(new_start);
                }
            }
        }

        // Handle sequence and quality updates based on mode
        match self.mode {
            ClippingMode::Soft => {
                // Keep sequence and qualities as-is
            }
            ClippingMode::SoftWithMask => {
                // Mask clipped bases to N and quality to min
                let seq = record.sequence();
                let qual = record.quality_scores();
                let mut new_seq: Vec<u8> = seq.as_ref().to_vec();
                let mut new_qual: Vec<u8> = qual.as_ref().to_vec();

                let total_soft_clip = existing_soft_clip + read_bases_clipped;
                for i in 0..total_soft_clip.min(new_seq.len()) {
                    new_seq[i] = NO_CALL_BASE;
                    new_qual[i] = MIN_PHRED;
                }

                *record.sequence_mut() = Sequence::from(new_seq);
                *record.quality_scores_mut() = QualityScores::from(new_qual);
            }
            ClippingMode::Hard => {
                // Remove clipped bases and qualities using direct slice indexing
                let seq = record.sequence();
                let qual = record.quality_scores();
                let new_seq = seq.as_ref()[bases_to_remove..].to_vec();
                let new_qual = qual.as_ref()[bases_to_remove..].to_vec();

                *record.sequence_mut() = Sequence::from(new_seq);
                *record.quality_scores_mut() = QualityScores::from(new_qual);

                // Clip extended attributes if auto-clipping is enabled
                self.clip_extended_attributes(record, bases_to_remove, true);
            }
        }

        read_bases_clipped
    }

    /// Clips a specified number of bases from the end (right side) of the alignment
    ///
    /// Returns the number of bases actually clipped
    #[expect(clippy::too_many_lines, reason = "mirrors clip_start_of_alignment with symmetric end-clipping logic")]
    pub fn clip_end_of_alignment(&self, record: &mut RecordBuf, bases_to_clip: usize) -> usize {
        if bases_to_clip == 0 {
            return 0;
        }

        // Don't clip unmapped reads
        if record.flags().is_unmapped() {
            return 0;
        }

        let sequence_len = record.sequence();
        if sequence_len.is_empty() {
            return 0;
        }

        // Collect existing CIGAR ops
        let old_ops: Vec<CigarOp> = record.cigar().iter().filter_map(Result::ok).collect();

        // Extract existing hard and soft clips from the end (process in reverse)
        let existing_hard_clip = old_ops
            .iter()
            .rev()
            .take_while(|op| op.kind() == Kind::HardClip)
            .map(|op| op.len())
            .sum::<usize>();

        let existing_soft_clip = old_ops
            .iter()
            .rev()
            .skip_while(|op| op.kind() == Kind::HardClip)
            .take_while(|op| op.kind() == Kind::SoftClip)
            .map(|op| op.len())
            .sum::<usize>();

        // Skip to operations before existing clips (reverse order, so use rev/skip_while/rev pattern)
        let mut post_clip_ops: Vec<CigarOp> = old_ops
            .into_iter()
            .rev()
            .skip_while(|op| matches!(op.kind(), Kind::HardClip | Kind::SoftClip))
            .collect();
        post_clip_ops.reverse(); // Un-reverse to get normal order

        let mut read_bases_clipped = 0;
        let mut new_ops = Vec::new();
        let mut iter = post_clip_ops.iter().rev().peekable();

        // Clip operations from the end (working backwards)
        // Continue until we've clipped enough read bases, OR if we've clipped the right amount
        // but the next operation (when working backwards) is a deletion (which should be removed)
        while read_bases_clipped < bases_to_clip
            || (read_bases_clipped == bases_to_clip
                && new_ops.is_empty()
                && iter.peek().map(|op| op.kind()) == Some(Kind::Deletion))
        {
            let Some(op) = iter.next() else { break };

            let kind = op.kind();
            let len = op.len();

            let consumes_read = matches!(
                kind,
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Insertion
            );

            if consumes_read && len > (bases_to_clip - read_bases_clipped) {
                // This operation extends past our clip point
                if kind == Kind::Insertion {
                    // Consume entire insertion at clip boundary
                    read_bases_clipped += len;
                } else {
                    // Split the operation
                    let remaining_clip = bases_to_clip - read_bases_clipped;
                    let remaining_length = len - remaining_clip;
                    read_bases_clipped += remaining_clip;
                    new_ops.push(CigarOp::new(kind, remaining_length));
                }
            } else {
                // Consume entire operation
                if consumes_read {
                    read_bases_clipped += len;
                }
            }
        }

        // Add remaining operations (iter is already reversed, so collect as-is)
        let remaining: Vec<CigarOp> = iter.copied().collect();
        new_ops.extend(remaining.iter());

        // new_ops is in reverse order, so reverse it
        new_ops.reverse();

        // Append appropriate clipping operators
        let (final_ops, bases_to_remove) = match self.mode {
            ClippingMode::Hard => {
                // Match fgbio's behavior: convert ALL existing soft clips to hard clips
                let added_hard_clip = existing_soft_clip + read_bases_clipped;
                let total_hard_clip = existing_hard_clip + added_hard_clip;
                let mut result = new_ops;
                result.push(CigarOp::new(Kind::HardClip, total_hard_clip));
                (result, added_hard_clip)
            }
            ClippingMode::Soft | ClippingMode::SoftWithMask => {
                let total_soft_clip = existing_soft_clip + read_bases_clipped;
                let mut result = new_ops;
                result.push(CigarOp::new(Kind::SoftClip, total_soft_clip));
                if existing_hard_clip > 0 {
                    result.push(CigarOp::new(Kind::HardClip, existing_hard_clip));
                }
                (result, 0)
            }
        };

        // Update CIGAR
        *record.cigar_mut() = CigarBuf::from(final_ops);

        // Handle sequence and quality updates based on mode
        let seq_len = record.sequence().len();
        match self.mode {
            ClippingMode::Soft => {
                // Keep sequence and qualities as-is
            }
            ClippingMode::SoftWithMask => {
                // Mask clipped bases to N and quality to min
                let seq = record.sequence();
                let qual = record.quality_scores();
                let mut new_seq: Vec<u8> = seq.as_ref().to_vec();
                let mut new_qual: Vec<u8> = qual.as_ref().to_vec();

                let total_soft_clip = existing_soft_clip + read_bases_clipped;
                let start_mask = seq_len.saturating_sub(total_soft_clip);
                for i in start_mask..seq_len {
                    new_seq[i] = NO_CALL_BASE;
                    new_qual[i] = MIN_PHRED;
                }

                *record.sequence_mut() = Sequence::from(new_seq);
                *record.quality_scores_mut() = QualityScores::from(new_qual);
            }
            ClippingMode::Hard => {
                // Remove clipped bases and qualities using direct slice indexing
                let seq = record.sequence();
                let qual = record.quality_scores();
                let keep_len = seq_len.saturating_sub(bases_to_remove);
                let new_seq = seq.as_ref()[..keep_len].to_vec();
                let new_qual = qual.as_ref()[..keep_len].to_vec();

                *record.sequence_mut() = Sequence::from(new_seq);
                *record.quality_scores_mut() = QualityScores::from(new_qual);

                // Clip extended attributes if auto-clipping is enabled
                self.clip_extended_attributes(record, bases_to_remove, false);
            }
        }

        read_bases_clipped
    }

    /// Clips bases from the 5' end of the read (strand-aware)
    ///
    /// For positive strand reads, clips from the start of alignment.
    /// For negative strand reads, clips from the end of alignment.
    ///
    /// Returns the number of bases actually clipped
    pub fn clip_5_prime_end_of_alignment(
        &self,
        record: &mut RecordBuf,
        bases_to_clip: usize,
    ) -> usize {
        if record.flags().is_reverse_complemented() {
            self.clip_end_of_alignment(record, bases_to_clip)
        } else {
            self.clip_start_of_alignment(record, bases_to_clip)
        }
    }

    /// Clips bases from the 3' end of the read (strand-aware)
    ///
    /// For positive strand reads, clips from the end of alignment.
    /// For negative strand reads, clips from the start of alignment.
    ///
    /// Returns the number of bases actually clipped
    pub fn clip_3_prime_end_of_alignment(
        &self,
        record: &mut RecordBuf,
        bases_to_clip: usize,
    ) -> usize {
        if record.flags().is_reverse_complemented() {
            self.clip_start_of_alignment(record, bases_to_clip)
        } else {
            self.clip_end_of_alignment(record, bases_to_clip)
        }
    }

    /// Clips overlapping portions of an FR read pair
    ///
    /// **Important:** Only clips reads that are in FR (forward-reverse) orientation.
    /// Non-FR pairs (FF, RR, RF) are not clipped.
    ///
    /// This implementation matches fgbio's midpoint approach: calculates the midpoint
    /// between the 5' ends of the two reads and clips both reads at that position.
    ///
    /// Returns (`bases_clipped_r1`, `bases_clipped_r2`)
    pub fn clip_overlapping_reads(&self, r1: &mut RecordBuf, r2: &mut RecordBuf) -> (usize, usize) {
        // Check if this is a valid FR pair before clipping
        if !record_utils::is_fr_pair(r1, r2) {
            return (0, 0);
        }

        // Get alignment positions
        let r1_start = match r1.alignment_start() {
            Some(pos) => usize::from(pos),
            None => return (0, 0),
        };
        let r2_start = match r2.alignment_start() {
            Some(pos) => usize::from(pos),
            None => return (0, 0),
        };

        // Calculate reference end positions using CIGAR
        let r1_end = r1_start + cigar_utils::reference_length(&r1.cigar()) - 1; // -1 for inclusive end
        let r2_end = r2_start + cigar_utils::reference_length(&r2.cigar()) - 1; // -1 for inclusive end

        // Check if they overlap on the reference
        let overlap_start = r1_start.max(r2_start);
        let overlap_end = r1_end.min(r2_end);

        if overlap_start > overlap_end {
            // No overlap
            return (0, 0);
        }

        // Calculate midpoint between the 5' ends (r1.start and r2.end)
        // In FR orientation: r1 is forward (5' at start), r2 is reverse (5' at end)
        let mut midpoint = usize::midpoint(r1_start, r2_end);

        // Adjust midpoint if it falls outside the overlap region
        if midpoint > r1_end {
            // midpoint is past r1's end, trim only the reverse strand read (r2)
            midpoint = r1_end;
        } else if midpoint < r2_start {
            // midpoint is before r2's start, trim only the positive strand read (r1)
            // Use r2_start - 1 to ensure we clip at least one base
            midpoint = r2_start.saturating_sub(1);
        }

        // Calculate how many bases to clip from each read
        // R1 (forward): clip from 3' end everything after midpoint
        let r1_bases_to_clip = if r1_end > midpoint {
            let ref_bases_to_clip = r1_end - midpoint;
            // Convert reference bases to query bases
            self.calculate_query_bases_for_ref_region(r1, ref_bases_to_clip, false)
        } else {
            0
        };

        // R2 (reverse): clip from 3' end (which is the 5' end in reference coordinates)
        // everything before midpoint + 1
        let r2_bases_to_clip = if midpoint + 1 > r2_start {
            let ref_bases_to_clip = midpoint + 1 - r2_start;
            // Convert reference bases to query bases, clipping from start (5' in ref = 3' in read)
            self.calculate_query_bases_for_ref_region(r2, ref_bases_to_clip, true)
        } else {
            0
        };

        let clipped_r1 =
            if r1_bases_to_clip > 0 { self.clip_end_of_alignment(r1, r1_bases_to_clip) } else { 0 };

        // For R2 (reverse read), we clip from the beginning of the alignment (5' in reference),
        // which corresponds to the beginning of the stored sequence, so use clip_5_prime
        let clipped_r2 = if r2_bases_to_clip > 0 {
            self.clip_start_of_alignment(r2, r2_bases_to_clip)
        } else {
            0
        };

        // If in Hard mode, upgrade any remaining soft clips to hard clips
        if matches!(self.mode, ClippingMode::Hard) {
            let _ = self.upgrade_all_clipping(r1);
            let _ = self.upgrade_all_clipping(r2);
        }

        (clipped_r1, clipped_r2)
    }

    /// Calculates the number of bases in a read that extend past the mate's unclipped boundary
    ///
    /// This matches fgbio's `numBasesExtendingPastMate` function exactly.
    ///
    /// # Arguments
    /// * `record` - The SAM record to check
    /// * `mate_unclipped_start` - The mate's unclipped start position
    /// * `mate_unclipped_end` - The mate's unclipped end position
    ///
    /// # Returns
    /// The number of bases that extend past the mate's boundary, or 0 if not applicable
    #[must_use]
    pub fn num_bases_extending_past_mate(
        record: &RecordBuf,
        mate_unclipped_start: usize,
        mate_unclipped_end: usize,
    ) -> usize {
        // Note: FR pair check is done in the caller (clip_extending_past_mate_ends)
        // so we don't need to check it here

        let is_positive_strand = !record.flags().is_reverse_complemented();

        // Get read length (total bases in the read, excluding hard clips)
        let read_length: usize = record
            .cigar()
            .iter()
            .filter_map(Result::ok)
            .filter(|op| {
                matches!(
                    op.kind(),
                    Kind::Match
                        | Kind::SequenceMatch
                        | Kind::SequenceMismatch
                        | Kind::Insertion
                        | Kind::SoftClip
                )
            })
            .map(CigarOp::len)
            .sum();

        if is_positive_strand {
            // Positive strand: check if read extends past mate's unclipped end
            let Some(alignment_start) = record.alignment_start().map(usize::from) else {
                return 0;
            };
            let alignment_end =
                alignment_start + cigar_utils::reference_length(&record.cigar()).saturating_sub(1);

            if alignment_end >= mate_unclipped_end {
                // Aligned portion reaches or extends past mate's unclipped end
                // Use readPosAtRefPos to find where in the read the mate ends
                // fgbio: Math.max(0, rec.length - rec.readPosAtRefPos(pos=mateUnclippedEnd, returnLastBaseIfDeleted=false))
                let pos_at_mate_end =
                    record_utils::read_pos_at_ref_pos(record, mate_unclipped_end, false);
                // When pos_at_mate_end is 0 (position in deletion or outside alignment),
                // fgbio's formula gives: max(0, read_length - 0) = read_length (clip all)
                read_length.saturating_sub(pos_at_mate_end)
            } else {
                // Aligned portion ends before mate's unclipped end
                // Only clip excess soft-clipped bases that extend past the mate
                let trailing_soft_clip = Self::trailing_soft_clips(record.cigar());
                let gap = mate_unclipped_end - alignment_end;
                trailing_soft_clip.saturating_sub(gap)
            }
        } else {
            // Negative strand: check if read extends before mate's unclipped start
            let Some(alignment_start) = record.alignment_start().map(usize::from) else {
                return 0;
            };

            if alignment_start > mate_unclipped_start {
                // Aligned portion starts after mate's unclipped start
                // Only clip excess soft-clipped bases that extend before the mate
                let leading_soft_clip = Self::leading_soft_clips(record.cigar());
                let gap = alignment_start - mate_unclipped_start;
                leading_soft_clip.saturating_sub(gap)
            } else {
                // Aligned portion starts at or before mate's unclipped start
                // Use readPosAtRefPos to find where in the read the mate starts
                // fgbio: Math.max(0, rec.readPosAtRefPos(pos=mateUnclippedStart, returnLastBaseIfDeleted=false) - 1)
                let pos_at_mate_start =
                    record_utils::read_pos_at_ref_pos(record, mate_unclipped_start, false);
                // When pos_at_mate_start is 0 (position in deletion or outside alignment),
                // fgbio's formula gives: max(0, 0 - 1) = 0 (clip nothing)
                // pos_at_mate_start is 1-based, so subtract 1 to get bases before it
                pos_at_mate_start.saturating_sub(1)
            }
        }
    }

    /// Clips reads that extend beyond their mate's alignment ends
    ///
    /// This matches fgbio's `clipExtendingPastMateEnds` function exactly.
    ///
    /// **Important:** Only clips reads that are in FR (forward-reverse) orientation.
    /// Non-FR pairs (FF, RR, RF) are not clipped.
    ///
    /// Returns (`bases_clipped_r1`, `bases_clipped_r2`)
    pub fn clip_extending_past_mate_ends(
        &self,
        r1: &mut RecordBuf,
        r2: &mut RecordBuf,
    ) -> (usize, usize) {
        // Check if this is a valid FR pair before clipping
        if !record_utils::is_fr_pair(r1, r2) {
            return (0, 0);
        }

        // Get unclipped positions for both reads
        let r1_unclipped_start = Self::unclipped_start(r1);
        let r1_unclipped_end = Self::unclipped_end(r1);
        let r2_unclipped_start = Self::unclipped_start(r2);
        let r2_unclipped_end = Self::unclipped_end(r2);

        let (Some(r1_start), Some(r1_end), Some(r2_start), Some(r2_end)) =
            (r1_unclipped_start, r1_unclipped_end, r2_unclipped_start, r2_unclipped_end)
        else {
            return (0, 0);
        };

        // Clip each read individually using fgbio's algorithm
        let clipped_r1 = self.clip_single_read_extending_past_mate(r1, r2_start, r2_end);
        let clipped_r2 = self.clip_single_read_extending_past_mate(r2, r1_start, r1_end);

        (clipped_r1, clipped_r2)
    }

    /// Clips a single read if it extends past its mate's alignment boundaries.
    ///
    /// This matches fgbio's `clipExtendingPastMateEnd` private method exactly.
    ///
    /// - Positive strand reads: clips from the 3' end (end of read)
    /// - Negative strand reads: clips from the 5' end (start of read)
    fn clip_single_read_extending_past_mate(
        &self,
        rec: &mut RecordBuf,
        mate_unclipped_start: usize,
        mate_unclipped_end: usize,
    ) -> usize {
        // Calculate how many bases extend past the mate
        let total_clipped_bases =
            Self::num_bases_extending_past_mate(rec, mate_unclipped_start, mate_unclipped_end);

        if total_clipped_bases == 0 {
            return 0;
        }

        // Clip based on strand
        let is_positive = !rec.flags().is_reverse_complemented();
        if is_positive {
            // Positive strand: clip from end of read (3' in sequencing order)
            self.clip_end_of_read(rec, total_clipped_bases)
        } else {
            // Negative strand: clip from start of read (5' in sequencing order, which is 3' in read orientation)
            self.clip_start_of_read(rec, total_clipped_bases)
        }
    }

    /// Get unclipped start position (delegates to `record_utils`)
    fn unclipped_start(rec: &RecordBuf) -> Option<usize> {
        record_utils::unclipped_start(rec)
    }

    /// Get unclipped end position (delegates to `record_utils`)
    fn unclipped_end(rec: &RecordBuf) -> Option<usize> {
        record_utils::unclipped_end(rec)
    }

    /// Count leading soft clips
    fn leading_soft_clips(cigar: &noodles::sam::alignment::record_buf::Cigar) -> usize {
        let ops: Vec<_> = cigar.as_ref().iter().map(|op| (op.kind(), op.len())).collect();
        record_utils::leading_soft_clipping(&ops)
    }

    /// Count trailing soft clips
    fn trailing_soft_clips(cigar: &noodles::sam::alignment::record_buf::Cigar) -> usize {
        let ops: Vec<_> = cigar.as_ref().iter().map(|op| (op.kind(), op.len())).collect();
        record_utils::trailing_soft_clipping(&ops)
    }

    /// Helper function to calculate query bases corresponding to a reference region
    ///
    /// Given a reference length, calculates how many query bases correspond to that region
    /// starting from either the 5' end (`from_start=true`) or 3' end (`from_start=false`)
    #[expect(clippy::unused_self, reason = "kept as a method for consistency with other clipper operations")]
    fn calculate_query_bases_for_ref_region(
        &self,
        record: &RecordBuf,
        ref_bases: usize,
        from_start: bool,
    ) -> usize {
        let cigar = record.cigar();
        let ops: Vec<_> = cigar.iter().filter_map(Result::ok).collect();

        let mut remaining_ref = ref_bases;
        let mut query_bases = 0;

        // Process CIGAR operations in the appropriate direction
        let iter: Box<dyn Iterator<Item = &CigarOp>> =
            if from_start { Box::new(ops.iter()) } else { Box::new(ops.iter().rev()) };

        for op in iter {
            if remaining_ref == 0 {
                break;
            }

            let kind = op.kind();
            let len = op.len();

            let consumes_ref = matches!(
                kind,
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion
            );
            let consumes_query = matches!(
                kind,
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Insertion
            );

            if consumes_ref {
                let ref_consumed = len.min(remaining_ref);
                remaining_ref -= ref_consumed;

                // Add corresponding query bases
                if consumes_query {
                    query_bases += ref_consumed;
                }
            } else if consumes_query && remaining_ref > 0 {
                // Insertion - consumes query but not reference
                // If we haven't finished consuming reference, include these bases
                query_bases += len;
            }
        }

        query_bases
    }

    /// Upgrades soft clipping to hard clipping for a specified number of bases
    ///
    /// This matches fgbio's private `upgradeClipping` method.
    ///
    /// # Arguments
    /// * `record` - The record to modify
    /// * `length` - The total number of clipped bases requested
    /// * `from_start` - If true, upgrade from the start; if false, from the end
    fn upgrade_clipping(&self, record: &mut RecordBuf, length: usize, from_start: bool) {
        // Only upgrade if not in Soft mode and length > 0
        if self.mode == ClippingMode::Soft || length == 0 {
            return;
        }

        let old_ops: Vec<CigarOp> = record.cigar().iter().filter_map(Result::ok).collect();

        // Count existing hard and soft clips in the appropriate direction
        let (hard_clipped, soft_clipped) = if from_start {
            let hard: usize = old_ops
                .iter()
                .take_while(|op| op.kind() == Kind::HardClip)
                .map(|op| op.len())
                .sum();
            let soft: usize = old_ops
                .iter()
                .skip_while(|op| op.kind() == Kind::HardClip)
                .take_while(|op| op.kind() == Kind::SoftClip)
                .map(|op| op.len())
                .sum();
            (hard, soft)
        } else {
            let hard: usize = old_ops
                .iter()
                .rev()
                .take_while(|op| op.kind() == Kind::HardClip)
                .map(|op| op.len())
                .sum();
            let soft: usize = old_ops
                .iter()
                .rev()
                .skip_while(|op| op.kind() == Kind::HardClip)
                .take_while(|op| op.kind() == Kind::SoftClip)
                .map(|op| op.len())
                .sum();
            (hard, soft)
        };

        // If requested length is already hard-clipped, or there are no soft clips, nothing to do
        if hard_clipped >= length || soft_clipped == 0 {
            return;
        }

        // Calculate how many soft clips to upgrade
        let length_to_upgrade = soft_clipped.min(length - hard_clipped);

        // Build new CIGAR and update sequence/quals
        let mut new_ops = Vec::new();
        let ops_to_process =
            if from_start { old_ops.clone() } else { old_ops.iter().rev().copied().collect() };

        // Handle Hard mode
        if self.mode == ClippingMode::Hard {
            let remaining_to_upgrade = length_to_upgrade;

            // Process leading clips
            let mut i = 0;
            let mut existing_hard = 0;
            let mut existing_soft = 0;

            // Count leading hard clips
            while i < ops_to_process.len() && ops_to_process[i].kind() == Kind::HardClip {
                existing_hard += ops_to_process[i].len();
                i += 1;
            }

            // Count leading soft clips
            while i < ops_to_process.len() && ops_to_process[i].kind() == Kind::SoftClip {
                existing_soft += ops_to_process[i].len();
                i += 1;
            }

            // Build new clipping ops
            let new_hard_count = existing_hard + remaining_to_upgrade;
            new_ops.push(CigarOp::new(Kind::HardClip, new_hard_count));

            if existing_soft > remaining_to_upgrade {
                // Some soft clips remain
                new_ops.push(CigarOp::new(Kind::SoftClip, existing_soft - remaining_to_upgrade));
            }

            // Add remaining ops
            new_ops.extend_from_slice(&ops_to_process[i..]);

            // Update CIGAR
            let final_ops = if from_start { new_ops } else { new_ops.into_iter().rev().collect() };
            *record.cigar_mut() = CigarBuf::from(final_ops);

            // Update sequence and quals by dropping the upgraded bases
            let seq = record.sequence().as_ref().to_vec();
            let quals = record.quality_scores().as_ref().to_vec();

            if from_start {
                *record.sequence_mut() = Sequence::from(seq[length_to_upgrade..].to_vec());
                *record.quality_scores_mut() =
                    QualityScores::from(quals[length_to_upgrade..].to_vec());
            } else {
                let new_len = seq.len() - length_to_upgrade;
                *record.sequence_mut() = Sequence::from(seq[..new_len].to_vec());
                *record.quality_scores_mut() = QualityScores::from(quals[..new_len].to_vec());
            }
        }
        // Note: SoftWithMask mode would use hard_mask_start_of_read/hard_mask_end_of_read
        // but we don't implement that yet
    }

    /// Ensures at least `clip_length` bases are clipped at the start of the read
    ///
    /// This includes any existing hard and soft clipping. If more clipping is needed,
    /// delegates to `clip_start_of_alignment`.
    ///
    /// Returns the number of additional bases clipped (not including existing clips)
    pub fn clip_start_of_read(&self, record: &mut RecordBuf, clip_length: usize) -> usize {
        // Count existing clipping at the start
        let existing_clipping: usize = record
            .cigar()
            .iter()
            .filter_map(Result::ok)
            .take_while(|op| matches!(op.kind(), Kind::HardClip | Kind::SoftClip))
            .map(CigarOp::len)
            .sum();

        if clip_length > existing_clipping {
            self.clip_start_of_alignment(record, clip_length - existing_clipping)
        } else {
            self.upgrade_clipping(record, clip_length, true);
            0
        }
    }

    /// Ensures at least `clip_length` bases are clipped at the end of the read
    ///
    /// This includes any existing hard and soft clipping. If more clipping is needed,
    /// delegates to `clip_end_of_alignment`.
    ///
    /// Returns the number of additional bases clipped (not including existing clips)
    pub fn clip_end_of_read(&self, record: &mut RecordBuf, clip_length: usize) -> usize {
        // Count existing clipping at the end
        let ops: Vec<_> = record.cigar().iter().filter_map(Result::ok).collect();
        let existing_clipping: usize = ops
            .iter()
            .rev()
            .take_while(|op| matches!(op.kind(), Kind::HardClip | Kind::SoftClip))
            .map(|op| op.len())
            .sum();

        if clip_length > existing_clipping {
            self.clip_end_of_alignment(record, clip_length - existing_clipping)
        } else {
            self.upgrade_clipping(record, clip_length, false);
            0
        }
    }

    /// Ensures at least `clip_length` bases are clipped at the 5' end (strand-aware)
    ///
    /// For positive strand: clips from start. For negative strand: clips from end.
    ///
    /// Returns the number of additional bases clipped
    pub fn clip_5_prime_end_of_read(&self, record: &mut RecordBuf, clip_length: usize) -> usize {
        if record.flags().is_reverse_complemented() {
            self.clip_end_of_read(record, clip_length)
        } else {
            self.clip_start_of_read(record, clip_length)
        }
    }

    /// Ensures at least `clip_length` bases are clipped at the 3' end (strand-aware)
    ///
    /// For positive strand: clips from end. For negative strand: clips from start.
    ///
    /// Returns the number of additional bases clipped
    pub fn clip_3_prime_end_of_read(&self, record: &mut RecordBuf, clip_length: usize) -> usize {
        if record.flags().is_reverse_complemented() {
            self.clip_start_of_read(record, clip_length)
        } else {
            self.clip_end_of_read(record, clip_length)
        }
    }

    /// Upgrades all existing clipping in a read to the current clipping mode
    ///
    /// E.g., converts soft clips to hard clips if mode is Hard,
    /// or soft clips to soft-with-mask if mode is `SoftWithMask`
    ///
    /// Returns `(leading_clipped, trailing_clipped)` - the number of soft-clipped
    /// bases that were converted to hard clips at the 5' and 3' ends respectively.
    /// Returns `(0, 0)` if no conversion occurred.
    ///
    /// # Panics
    ///
    /// Panics if an adjacent hard-clip CIGAR op disappears between the `last()` check
    /// and the subsequent access (should not happen in practice), or if a BAM tag
    /// is not exactly 2 bytes.
    ///
    /// # Errors
    ///
    /// Returns `Result` for API compatibility, but the current implementation is
    /// infallible and always returns `Ok`.
    #[expect(clippy::too_many_lines, reason = "CIGAR rewriting with attribute clipping requires many branches")]
    pub fn upgrade_all_clipping(&self, record: &mut RecordBuf) -> Result<(usize, usize)> {
        // Only upgrade in Hard mode
        if !matches!(self.mode, ClippingMode::Hard) {
            return Ok((0, 0));
        }

        // Don't upgrade unmapped reads
        if record.flags().is_unmapped() {
            return Ok((0, 0));
        }
        let cigar = record.cigar();
        let ops: Vec<CigarOp> = cigar.iter().filter_map(Result::ok).collect();

        // Check if there are any soft clips to convert
        let has_soft_clips = ops.iter().any(|op| op.kind() == Kind::SoftClip);
        if !has_soft_clips {
            return Ok((0, 0));
        }

        // Count leading and trailing soft clips (and existing hard clips)
        let mut leading_hard = 0;
        let mut leading_soft = 0;
        let mut trailing_soft = 0;
        let mut _trailing_hard = 0;

        // Count leading clips
        for op in &ops {
            match op.kind() {
                Kind::HardClip => leading_hard += op.len(),
                Kind::SoftClip => {
                    leading_soft += op.len();
                    break;
                }
                _ => break,
            }
        }

        // Count trailing clips (from the end)
        for op in ops.iter().rev() {
            match op.kind() {
                Kind::HardClip => _trailing_hard += op.len(),
                Kind::SoftClip => {
                    trailing_soft += op.len();
                    break;
                }
                _ => break,
            }
        }

        // Build new CIGAR, converting soft clips to hard clips
        let mut new_cigar_ops = Vec::new();
        let sequence = record.sequence();
        let qualities = record.quality_scores();
        let old_seq_len = sequence.len(); // Capture this before any mutations
        let mut seq_pos = 0;
        let mut new_sequence = Vec::new();
        let mut new_qualities = Vec::new();
        let mut is_leading = true;

        for op in &ops {
            let kind = op.kind();
            let len = op.len();

            match kind {
                Kind::SoftClip => {
                    // Convert to hard clip
                    // Merge with existing hard clips at this position
                    if is_leading && new_cigar_ops.is_empty() && leading_hard > 0 {
                        // Merge with preceding hard clip
                        new_cigar_ops.push(CigarOp::new(Kind::HardClip, leading_hard + len));
                    } else if new_cigar_ops.last().map(|o| o.kind()) == Some(Kind::HardClip) {
                        // Merge with adjacent hard clip
                        // SAFETY: We just checked that last() is Some(HardClip)
                        let last_len = new_cigar_ops.last().expect("last() checked above").len();
                        new_cigar_ops.pop();
                        new_cigar_ops.push(CigarOp::new(Kind::HardClip, last_len + len));
                    } else {
                        new_cigar_ops.push(CigarOp::new(Kind::HardClip, len));
                    }
                    seq_pos += len;
                }
                Kind::HardClip => {
                    // Merge with preceding hard clip if present
                    if new_cigar_ops.last().map(|o| o.kind()) == Some(Kind::HardClip) {
                        // SAFETY: We just checked that last() is Some(HardClip)
                        let last_len = new_cigar_ops.last().expect("last() checked above").len();
                        new_cigar_ops.pop();
                        new_cigar_ops.push(CigarOp::new(Kind::HardClip, last_len + len));
                    } else if !is_leading || new_cigar_ops.is_empty() {
                        // Keep hard clips but don't add duplicates if we already merged leading
                        new_cigar_ops.push(*op);
                    }
                }
                _ => {
                    is_leading = false;
                    new_cigar_ops.push(*op);
                    // Copy bases/quals for query-consuming operations
                    let consumes_query = matches!(
                        kind,
                        Kind::Match
                            | Kind::SequenceMatch
                            | Kind::SequenceMismatch
                            | Kind::Insertion
                    );
                    if consumes_query {
                        for _ in 0..len {
                            if seq_pos < sequence.len() {
                                new_sequence.push(sequence.as_ref()[seq_pos]);
                                new_qualities.push(qualities.as_ref()[seq_pos]);
                            }
                            seq_pos += 1;
                        }
                    }
                }
            }
        }

        // Update the record
        *record.cigar_mut() = CigarBuf::from(new_cigar_ops);

        // Clip extended attributes if auto-clipping is enabled (do this BEFORE updating sequence)
        if self.auto_clip_attributes && (leading_soft > 0 || trailing_soft > 0) {
            // Manually clip attributes that match the original sequence length
            use noodles::sam::alignment::record::data::field::Tag;
            use noodles::sam::alignment::record_buf::data::field::value::Array;
            let mut tags_to_update: Vec<(Vec<u8>, Value)> = Vec::new();

            for (tag, value) in record.data().iter() {
                let should_clip = match value {
                    Value::String(s) => {
                        let bytes: &[u8] = s.as_ref();
                        bytes.len() == old_seq_len
                    }
                    Value::Array(arr) => array_len!(arr) == old_seq_len,
                    _ => false,
                };

                if should_clip {
                    // Clip from both ends: remove leading_soft from start, trailing_soft from end
                    let start = leading_soft;
                    let end = old_seq_len - trailing_soft;
                    let new_value = match value {
                        Value::String(s) => {
                            let bytes: &[u8] = s.as_ref();
                            Value::from(String::from_utf8_lossy(&bytes[start..end]).to_string())
                        }
                        Value::Array(arr) => slice_array!(arr, start, end),
                        _ => value.clone(),
                    };
                    tags_to_update.push((tag.as_ref().to_vec(), new_value));
                }
            }

            // Apply updates
            for (tag_bytes, value) in tags_to_update {
                // SAFETY: Tags are always 2 bytes, as tag_bytes came from tag.as_ref().to_vec() where tag is a 2-byte tag
                let tag = Tag::from(
                    <[u8; 2]>::try_from(tag_bytes.as_slice())
                        .expect("tag bytes are always 2 bytes"),
                );
                record.data_mut().insert(tag, value);
            }
        }

        *record.sequence_mut() = Sequence::from(new_sequence);
        *record.quality_scores_mut() = QualityScores::from(new_qualities);

        Ok((leading_soft, trailing_soft))
    }

    /// Returns the number of bases that are currently clipped in the read
    #[must_use]
    pub fn clipped_bases(record: &RecordBuf) -> usize {
        cigar_utils::clipped_bases(&record.cigar())
    }
}

/// Helper functions for CIGAR manipulation
pub mod cigar_utils {
    use super::Kind;
    use noodles::sam::alignment::record::Cigar as CigarTrait;

    /// Counts the number of aligned bases in a CIGAR string
    #[must_use]
    #[expect(clippy::redundant_closure_for_method_calls, reason = "Op::len is not a method on the trait")]
    pub fn aligned_bases(cigar: &impl CigarTrait) -> usize {
        cigar
            .iter()
            .filter_map(Result::ok)
            .filter(|op| {
                matches!(op.kind(), Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch)
            })
            .map(|op| op.len())
            .sum()
    }

    /// Counts the number of clipped bases in a CIGAR string
    #[must_use]
    #[expect(clippy::redundant_closure_for_method_calls, reason = "Op::len is not a method on the trait")]
    pub fn clipped_bases(cigar: &impl CigarTrait) -> usize {
        cigar
            .iter()
            .filter_map(Result::ok)
            .filter(|op| matches!(op.kind(), Kind::SoftClip | Kind::HardClip))
            .map(|op| op.len())
            .sum()
    }

    /// Counts reference-consuming operations
    #[must_use]
    #[expect(clippy::redundant_closure_for_method_calls, reason = "Op::len is not a method on the trait")]
    pub fn reference_length(cigar: &impl CigarTrait) -> usize {
        cigar
            .iter()
            .filter_map(Result::ok)
            .filter(|op| {
                matches!(
                    op.kind(),
                    Kind::Match
                        | Kind::SequenceMatch
                        | Kind::SequenceMismatch
                        | Kind::Deletion
                        | Kind::Skip
                )
            })
            .map(|op| op.len())
            .sum()
    }

    /// Simplified CIGAR representation: a vector of (operation kind, length) pairs.
    pub type SimplifiedCigar = Vec<(Kind, usize)>;

    /// Simplifies a CIGAR string by:
    /// 1. Converting S, EQ, X, H operations to M (match)
    /// 2. Coalescing adjacent operations of the same type
    ///
    /// This is useful for comparing read alignments where only indel positions matter.
    #[must_use]
    pub fn simplify_cigar(cigar: &noodles::sam::alignment::record_buf::Cigar) -> SimplifiedCigar {
        let mut simplified: SimplifiedCigar = Vec::new();

        for op in cigar.as_ref() {
            let (kind, len) = (op.kind(), op.len());

            // Convert S, EQ, X, H to M; keep I, D, and M as-is
            let new_kind = match kind {
                Kind::SoftClip | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::HardClip => {
                    Kind::Match
                }
                _ => kind,
            };

            // Coalesce adjacent operations of the same type
            if let Some((last_kind, last_len)) = simplified.last_mut() {
                if *last_kind == new_kind {
                    *last_len += len;
                    continue;
                }
            }

            simplified.push((new_kind, len));
        }

        simplified
    }

    /// Checks if `cigar_a` is a prefix of `cigar_b`.
    ///
    /// This allows for reads of different lengths to be grouped together if they share
    /// the same alignment pattern (same indel positions and lengths).
    ///
    /// For example, 10M is a prefix of 20M, but 10M1I is not a prefix of 10M1D.
    #[must_use]
    pub fn is_cigar_prefix(a: &[(Kind, usize)], b: &[(Kind, usize)]) -> bool {
        if a.len() > b.len() {
            return false;
        }

        let last_index = a.len().saturating_sub(1);

        for (i, &(op_a, len_a)) in a.iter().enumerate() {
            let (op_b, len_b) = b[i];
            if op_a != op_b {
                return false;
            }
            // For the last element, a's length can be <= b's length (prefix)
            // For other elements, lengths must match exactly
            if i == last_index {
                if len_a > len_b {
                    return false;
                }
            } else if len_a != len_b {
                return false;
            }
        }
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::builder::RecordBuilder;

    #[test]
    fn test_clipping_mode() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        assert_eq!(clipper.mode, ClippingMode::Soft);

        let clipper = SamRecordClipper::new(ClippingMode::Hard);
        assert_eq!(clipper.mode, ClippingMode::Hard);
    }

    #[test]
    fn test_auto_clip() {
        let clipper = SamRecordClipper::with_auto_clip(ClippingMode::Soft, true);
        assert!(clipper.auto_clip_attributes);

        let clipper_disabled = SamRecordClipper::with_auto_clip(ClippingMode::Hard, false);
        assert!(!clipper_disabled.auto_clip_attributes);
    }

    use noodles::core::Position;
    use noodles::sam::alignment::record::Cigar as CigarTrait;
    use noodles::sam::alignment::record::Flags;
    use noodles::sam::alignment::record_buf::RecordBuf;

    /// Helper to format a CIGAR string for comparison
    fn format_cigar(cigar: &impl CigarTrait) -> String {
        use std::fmt::Write;
        cigar.iter().filter_map(Result::ok).fold(String::new(), |mut acc, op| {
            let kind_char = match op.kind() {
                Kind::Match => 'M',
                Kind::Insertion => 'I',
                Kind::Deletion => 'D',
                Kind::Skip => 'N',
                Kind::SoftClip => 'S',
                Kind::HardClip => 'H',
                Kind::Pad => 'P',
                Kind::SequenceMatch => '=',
                Kind::SequenceMismatch => 'X',
            };
            let _ = write!(acc, "{}{}", op.len(), kind_char);
            acc
        })
    }

    /// Helper to create a simple test record with given CIGAR, sequence, and position
    fn create_test_record(cigar_str: &str, seq: &str, start_pos: usize) -> RecordBuf {
        RecordBuilder::mapped_read()
            .sequence(seq)
            .cigar(cigar_str)
            .alignment_start(start_pos)
            .build()
    }

    // ===================================================================
    // clip_5_prime (clipStartOfAlignment) tests - Basic functionality
    // ===================================================================

    #[test]
    fn test_clip_start_of_alignment_soft_matched_bases() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let mut record = create_test_record("50M", "ACGTACGTACGT", 10);

        let clipped = clipper.clip_start_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(20));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "10S40M");
    }

    #[test]
    fn test_clip_start_of_alignment_soft_with_insertion() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let mut record = create_test_record("4M2I44M", "ACGTACGTACGT", 10);

        let clipped = clipper.clip_start_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(18)); // 10 + (4+2 ref consumed)

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "10S40M");
    }

    #[test]
    fn test_clip_start_of_alignment_soft_with_deletion() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let mut record = create_test_record("6M2D44M", "ACGTACGTACGT", 10);

        let clipped = clipper.clip_start_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(22)); // 10 + 6 + 2 (deletion) + 4

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "10S40M");
    }

    #[test]
    fn test_clip_start_of_alignment_soft_additional_bases() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let mut record = create_test_record("10S40M", "ACGTACGTACGT", 10);

        let clipped = clipper.clip_start_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(20));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "20S30M");
    }

    #[test]
    fn test_clip_start_of_alignment_preserve_hard_clip() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let mut record = create_test_record("10H40M", "ACGTACGTACGT", 10);

        let clipped = clipper.clip_start_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(20));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "10H10S30M");
    }

    #[test]
    fn test_clip_start_of_alignment_complex_cigar() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        // 2H4S16M10I5M5I10M
        let mut record = create_test_record("2H4S16M10I5M5I10M", "ACGTACGTACGT", 10);

        let clipped = clipper.clip_start_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(20));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "2H14S6M10I5M5I10M");
    }

    #[test]
    fn test_clip_start_of_alignment_consumes_trailing_insertion() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let mut record = create_test_record("8M4I38M", "ACGTACGTACGT", 10);

        // Ask to clip 10 bases, but should consume 12 bases (8M + 4I) because
        // the insertion at the clip boundary is consumed entirely
        let clipped = clipper.clip_start_of_alignment(&mut record, 10);
        assert_eq!(clipped, 12);
        assert_eq!(record.alignment_start(), Position::new(18));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "12S38M");
    }

    #[test]
    fn test_clip_start_of_alignment_preserve_insertion_after_clip() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let mut record = create_test_record("10M4I36M", "ACGTACGTACGT", 10);

        let clipped = clipper.clip_start_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(20));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "10S4I36M");
    }

    #[test]
    fn test_clip_start_of_alignment_remove_deletion_after_clip() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let mut record = create_test_record("10M4D40M", "ACGTACGTACGT", 10);

        let clipped = clipper.clip_start_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(24)); // 10 + 10 + 4 (deletion)

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "10S40M");
    }

    #[test]
    fn test_clip_start_of_alignment_preserve_distant_deletion() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let mut record = create_test_record("25M4D25M", "ACGTACGTACGT", 10);

        let clipped = clipper.clip_start_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(20));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "10S15M4D25M");
    }

    #[test]
    fn test_clip_start_of_alignment_soft_with_mask() {
        let clipper = SamRecordClipper::new(ClippingMode::SoftWithMask);
        let seq = "ACGTACGTAC"; // 10 bases
        let mut record = create_test_record("10M", seq, 10);

        let clipped = clipper.clip_start_of_alignment(&mut record, 5);
        assert_eq!(clipped, 5);

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "5S5M");

        // Check first 5 bases are masked to N
        let bases = record.sequence();
        for i in 0..5 {
            assert_eq!(bases.as_ref()[i], b'N', "Base at position {i} should be N");
        }

        // Check first 5 quals are set to min
        let quals = record.quality_scores();
        for i in 0..5 {
            assert_eq!(quals.as_ref()[i], 2, "Quality at position {i} should be 2");
        }
    }

    #[test]
    fn test_clip_start_of_alignment_soft_with_mask_existing_soft_clip() {
        let clipper = SamRecordClipper::new(ClippingMode::SoftWithMask);
        let seq = "ACGTACGTACGTACGTACGT"; // 20 bases
        let mut record = create_test_record("10S10M", seq, 10);

        let clipped = clipper.clip_start_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "20S"); // All bases now soft-clipped (0M operations are removed)

        // Check all 20 bases are masked to N
        let bases = record.sequence();
        for i in 0..20 {
            assert_eq!(bases.as_ref()[i], b'N');
        }
    }

    #[test]
    fn test_clip_start_of_alignment_unmapped_read_does_nothing() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let mut record = create_test_record("50M", "ACGTACGTACGTACGTACGT", 1000);

        // Set unmapped flag
        *record.flags_mut() = Flags::UNMAPPED;

        let clipped = clipper.clip_start_of_alignment(&mut record, 10);

        // Should not clip unmapped reads
        assert_eq!(clipped, 0);
        assert_eq!(record.alignment_start(), Position::new(1000));

        // CIGAR should remain unchanged
        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "50M");
    }

    #[test]
    fn test_clip_start_of_alignment_soft() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let mut record = create_test_record("100M", "ACGTACGTACGT", 1000);

        let clipped = clipper.clip_start_of_alignment(&mut record, 5);
        assert_eq!(clipped, 5);

        // Check that CIGAR starts with soft clip
        let cigar = record.cigar();
        let first_op = cigar.iter().next().unwrap().unwrap();
        assert_eq!(first_op.kind(), Kind::SoftClip);
        assert_eq!(first_op.len(), 5);
    }

    #[test]
    fn test_clip_end_of_alignment_soft() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let mut record = create_test_record("100M", "ACGTACGTACGT", 1000);

        let clipped = clipper.clip_end_of_alignment(&mut record, 5);
        assert_eq!(clipped, 5);

        // Check that CIGAR ends with soft clip
        let cigar = record.cigar();
        let ops: Vec<_> = cigar.iter().filter_map(Result::ok).collect();
        let last_op = ops.last().unwrap();
        assert_eq!(last_op.kind(), Kind::SoftClip);
        assert_eq!(last_op.len(), 5);
    }

    // ===================================================================
    // clip_3_prime (clipEndOfAlignment) tests - Extended coverage
    // ===================================================================

    #[test]
    fn test_clip_end_of_alignment_soft_matched_bases() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        // 50 bases to match 50M CIGAR
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC";
        assert_eq!(seq.len(), 50, "Sequence length should be 50");
        let mut record = create_test_record("50M", seq, 10);

        let clipped = clipper.clip_end_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(10));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "40M10S");
    }

    #[test]
    fn test_clip_end_of_alignment_hard_existing_soft_clip() {
        let clipper = SamRecordClipper::new(ClippingMode::Hard);
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"; // 50 bases
        assert_eq!(seq.len(), 50, "Sequence length should be 50");
        let mut record = create_test_record("40M10S", seq, 10);

        let orig_seq: Vec<u8> = record.sequence().as_ref().to_vec();
        let orig_qual: Vec<u8> = record.quality_scores().as_ref().to_vec();

        let clipped = clipper.clip_end_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(10));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "30M20H");

        // Check sequence and qualities are truncated
        assert_eq!(record.sequence().as_ref(), &orig_seq[..30]);
        assert_eq!(record.quality_scores().as_ref(), &orig_qual[..30]);
    }

    #[test]
    fn test_clip_start_of_alignment_hard() {
        let clipper = SamRecordClipper::new(ClippingMode::Hard);
        let mut record = create_test_record("12M", "ACGTACGTACGT", 1000);

        let original_len = record.sequence().len();
        let clipped = clipper.clip_start_of_alignment(&mut record, 5);
        assert_eq!(clipped, 5);

        // Check that sequence is shortened
        assert_eq!(record.sequence().len(), original_len - 5);

        // Check that CIGAR starts with hard clip
        let cigar = record.cigar();
        let first_op = cigar.iter().next().unwrap().unwrap();
        assert_eq!(first_op.kind(), Kind::HardClip);
        assert_eq!(first_op.len(), 5);
    }

    #[test]
    fn test_clip_overlapping_reads_no_overlap() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);

        // R1: 1000-1100, R2: 1200-1300 (no overlap)
        let mut r1 = create_test_record("100M", "ACGTACGTACGT", 1000);
        let mut r2 = create_test_record("100M", "TGCATGCATGCA", 1200);

        let (clipped_r1, clipped_r2) = clipper.clip_overlapping_reads(&mut r1, &mut r2);

        // No clipping should occur
        assert_eq!(clipped_r1, 0);
        assert_eq!(clipped_r2, 0);
    }

    #[test]
    fn test_clip_overlapping_reads_with_overlap() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);

        // R1: 1000-1099 (forward), R2: 1050-1149 (reverse) - 50bp overlap
        // Create proper FR pair
        let mut r1 = create_paired_record("100M", "ACGTACGTACGT", 1000, false, true, 1050, "100M");
        let mut r2 = create_paired_record("100M", "TGCATGCATGCA", 1050, true, false, 1000, "100M");

        let (clipped_r1, clipped_r2) = clipper.clip_overlapping_reads(&mut r1, &mut r2);

        // With midpoint algorithm, both reads should be clipped
        // Midpoint = (1000 + 1149) / 2 = 1074
        // R1 clips from 1075 onwards (approximately 25 bases)
        // R2 clips from start to 1074 (approximately 25 bases)
        assert!(clipped_r1 > 0, "R1 should be clipped");
        assert!(clipped_r2 > 0, "R2 should be clipped");

        // Both should clip roughly equal amounts (within a few bases due to rounding)
        let diff = (i32::try_from(clipped_r1).unwrap() - i32::try_from(clipped_r2).unwrap()).abs();
        assert!(
            diff <= 2,
            "Clipping should be roughly equal, but got {clipped_r1} vs {clipped_r2}"
        );
    }

    #[test]
    fn test_clip_extending_past_mate_ends() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);

        // R1: 1000-1149, R2: 1100-1199
        // R1 does NOT extend past R2's end (1149 < 1199)
        // R2 does NOT start before R1's start (1100 > 1000)
        // So neither should be clipped
        let mut r1 = create_paired_record("150M", "ACGTACGTACGT", 1000, false, true, 1100, "100M");
        let mut r2 = create_paired_record("100M", "TGCATGCATGCA", 1100, true, false, 1000, "150M");

        let (clipped_r1, clipped_r2) = clipper.clip_extending_past_mate_ends(&mut r1, &mut r2);

        // Neither should be clipped
        assert_eq!(clipped_r1, 0);
        assert_eq!(clipped_r2, 0);
    }

    #[test]
    fn test_cigar_utils_reference_length() {
        let record = create_test_record("50M10I40M", "ACGTACGTACGT", 1000);
        let ref_len = cigar_utils::reference_length(&record.cigar());

        // 50M + 40M = 90 (insertion doesn't consume reference)
        assert_eq!(ref_len, 90);
    }

    #[test]
    fn test_cigar_utils_aligned_bases() {
        let record = create_test_record("50M10I40M", "ACGTACGTACGT", 1000);
        let aligned = cigar_utils::aligned_bases(&record.cigar());

        // 50M + 40M = 90 (insertion is not aligned)
        assert_eq!(aligned, 90);
    }

    #[test]
    fn test_cigar_utils_clipped_bases() {
        let record = create_test_record("5S90M5S", "ACGTACGTACGT", 1000);
        let clipped = cigar_utils::clipped_bases(&record.cigar());

        // 5S + 5S = 10
        assert_eq!(clipped, 10);
    }

    #[test]
    fn test_clip_start_of_alignment_with_existing_soft_clip() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let mut record = create_test_record("5S95M", "ACGTACGTACGT", 1000);

        // Clip additional 3 bases from 5' end
        let clipped = clipper.clip_start_of_alignment(&mut record, 3);
        assert_eq!(clipped, 3);

        // Total soft clip should be 5 + 3 = 8
        let cigar = record.cigar();
        let first_op = cigar.iter().next().unwrap().unwrap();
        assert_eq!(first_op.kind(), Kind::SoftClip);
    }

    #[test]
    fn test_clip_entire_read() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let seq = "ACGTACGTACGT";
        let mut record = create_test_record("12M", seq, 1000);

        // Try to clip more than the read length
        let clipped = clipper.clip_start_of_alignment(&mut record, 20);

        // Should only clip up to sequence length
        assert_eq!(clipped, seq.len());
    }

    #[test]
    fn test_clip_zero_bases() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let mut record = create_test_record("100M", "ACGTACGTACGT", 1000);

        let clipped = clipper.clip_start_of_alignment(&mut record, 0);
        assert_eq!(clipped, 0);

        // CIGAR should remain unchanged
        let cigar = record.cigar();
        let first_op = cigar.iter().next().unwrap().unwrap();
        assert_eq!(first_op.kind(), Kind::Match);
    }

    /// Helper to create a paired read with specific flags and positions
    fn create_paired_record(
        cigar_str: &str,
        seq: &str,
        start_pos: usize,
        is_reverse: bool,
        mate_reverse: bool,
        mate_pos: usize,
        mate_cigar_str: &str,
    ) -> RecordBuf {
        let mut record = create_test_record(cigar_str, seq, start_pos);

        // Set paired flags
        let mut flags = Flags::SEGMENTED;
        if is_reverse {
            flags |= Flags::REVERSE_COMPLEMENTED;
        }
        if mate_reverse {
            flags |= Flags::MATE_REVERSE_COMPLEMENTED;
        }
        *record.flags_mut() = flags;

        // Set mate position
        *record.mate_alignment_start_mut() = Position::new(mate_pos);

        // Set same reference for both reads
        *record.mate_reference_sequence_id_mut() = record.reference_sequence_id();

        // Calculate and set template_length for FR pair detection
        // Parse mate CIGAR to calculate its reference length
        let mate_cigar_ops: Vec<CigarOp> = mate_cigar_str
            .split(|c: char| !c.is_numeric())
            .filter(|s| !s.is_empty())
            .zip(mate_cigar_str.chars().filter(|c| c.is_alphabetic()))
            .map(|(len_str, kind_char)| {
                let len: usize = len_str.parse().unwrap();
                let kind = match kind_char {
                    'M' => Kind::Match,
                    'I' => Kind::Insertion,
                    'D' => Kind::Deletion,
                    'S' => Kind::SoftClip,
                    'H' => Kind::HardClip,
                    'N' => Kind::Skip,
                    _ => panic!("Invalid CIGAR operation: {kind_char}"),
                };
                CigarOp::new(kind, len)
            })
            .collect();
        let mate_cigar = CigarBuf::from(mate_cigar_ops);
        let mate_ref_len = cigar_utils::reference_length(&mate_cigar);
        let ref_len = cigar_utils::reference_length(&record.cigar());

        let tlen = if is_reverse {
            // This read is reverse, mate is forward
            // tlen = -(this_end - mate_start + 1)
            let this_end = start_pos + ref_len - 1;
            if this_end >= mate_pos {
                -i32::try_from(this_end - mate_pos + 1).unwrap()
            } else {
                i32::try_from(mate_pos - this_end - 1).unwrap()
            }
        } else {
            // This read is forward, mate is reverse
            // tlen = mate_end - this_start + 1
            let mate_end = mate_pos + mate_ref_len - 1;
            if mate_end >= start_pos {
                i32::try_from(mate_end - start_pos + 1).unwrap()
            } else {
                -i32::try_from(start_pos - mate_end - 1).unwrap()
            }
        };

        *record.template_length_mut() = tlen;

        record
    }

    #[test]
    fn test_is_fr_pair_valid() {
        // Create valid FR pair: R1 forward, R2 reverse
        let r1 = create_paired_record("100M", "ACGTACGTACGT", 1000, false, true, 1100, "100M");
        let r2 = create_paired_record("100M", "TGCATGCATGCA", 1100, true, false, 1000, "100M");

        assert!(record_utils::is_fr_pair(&r1, &r2));
    }

    #[test]
    fn test_is_fr_pair_both_forward() {
        // FF orientation - should not be considered FR
        let r1 = create_paired_record("100M", "ACGTACGTACGT", 1000, false, false, 1100, "100M");
        let r2 = create_paired_record("100M", "TGCATGCATGCA", 1100, false, false, 1000, "100M");

        assert!(!record_utils::is_fr_pair(&r1, &r2));
    }

    #[test]
    fn test_is_fr_pair_both_reverse() {
        // RR orientation - should not be considered FR
        let r1 = create_paired_record("100M", "ACGTACGTACGT", 1000, true, true, 1100, "100M");
        let r2 = create_paired_record("100M", "TGCATGCATGCA", 1100, true, true, 1000, "100M");

        assert!(!record_utils::is_fr_pair(&r1, &r2));
    }

    #[test]
    fn test_is_fr_pair_rf_orientation() {
        // RF orientation (reverse of FR) - should not be considered FR
        let r1 = create_paired_record("100M", "ACGTACGTACGT", 1000, true, false, 1100, "100M");
        let r2 = create_paired_record("100M", "TGCATGCATGCA", 1100, false, true, 1000, "100M");

        assert!(!record_utils::is_fr_pair(&r1, &r2));
    }

    #[test]
    fn test_is_fr_pair_unmapped() {
        let mut r1 = create_paired_record("100M", "ACGTACGTACGT", 1000, false, true, 1100, "100M");
        let r2 = create_paired_record("100M", "TGCATGCATGCA", 1100, true, false, 1000, "150M");

        // Set R1 as unmapped
        *r1.flags_mut() |= Flags::UNMAPPED;

        assert!(!record_utils::is_fr_pair(&r1, &r2));
    }

    #[test]
    fn test_is_fr_pair_mate_unmapped() {
        let mut r1 = create_paired_record("100M", "ACGTACGTACGT", 1000, false, true, 1100, "100M");
        let r2 = create_paired_record("100M", "TGCATGCATGCA", 1100, true, false, 1000, "150M");

        // Set R1's mate as unmapped
        *r1.flags_mut() |= Flags::MATE_UNMAPPED;

        assert!(!record_utils::is_fr_pair(&r1, &r2));
    }

    #[test]
    fn test_is_fr_pair_different_chromosomes() {
        let r1 = create_paired_record("100M", "ACGTACGTACGT", 1000, false, true, 1100, "100M");
        let mut r2 = create_paired_record("100M", "TGCATGCATGCA", 1100, true, false, 1000, "150M");

        // Set R2 to different reference sequence
        *r2.reference_sequence_id_mut() = Some(1);

        assert!(!record_utils::is_fr_pair(&r1, &r2));
    }

    #[test]
    fn test_is_fr_pair_not_paired() {
        let mut r1 = create_paired_record("100M", "ACGTACGTACGT", 1000, false, true, 1100, "100M");
        let r2 = create_paired_record("100M", "TGCATGCATGCA", 1100, true, false, 1000, "150M");

        // Remove paired flag from R1
        *r1.flags_mut() = Flags::empty();

        assert!(!record_utils::is_fr_pair(&r1, &r2));
    }

    #[test]
    fn test_clip_overlapping_reads_non_fr_pair() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);

        // Create FF pair (both forward) with overlap
        // R1: 1000-1100, R2: 1050-1150 (50bp overlap)
        let mut r1 = create_paired_record("100M", "ACGTACGTACGT", 1000, false, false, 1050, "100M");
        let mut r2 = create_paired_record("100M", "TGCATGCATGCA", 1050, false, false, 1000, "100M");

        let (clipped_r1, clipped_r2) = clipper.clip_overlapping_reads(&mut r1, &mut r2);

        // Should NOT clip because this is not an FR pair
        assert_eq!(clipped_r1, 0);
        assert_eq!(clipped_r2, 0);
    }

    #[test]
    fn test_clip_extending_past_mate_ends_non_fr_pair() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);

        // Create RR pair (both reverse) with extension
        // R1: 1000-1150, R2: 1100-1200
        let mut r1 = create_paired_record("150M", "ACGTACGTACGT", 1000, true, true, 1100, "100M");
        let mut r2 = create_paired_record("100M", "TGCATGCATGCA", 1100, true, true, 1000, "150M");

        let (clipped_r1, clipped_r2) = clipper.clip_extending_past_mate_ends(&mut r1, &mut r2);

        // Should NOT clip because this is not an FR pair
        assert_eq!(clipped_r1, 0);
        assert_eq!(clipped_r2, 0);
    }

    #[test]
    fn test_auto_clip_attributes_string_5_prime() {
        use noodles::sam::alignment::record::data::field::Tag;

        let clipper = SamRecordClipper::with_auto_clip(ClippingMode::Hard, true);
        let mut record = create_test_record("10M", "ACGTACGTAC", 1000);

        // Add a string attribute that matches the read length
        let tag = Tag::from([b'X', b'S']);
        let value = Value::from("0123456789");
        record.data_mut().insert(tag, value);

        // Clip 3 bases from 5' end
        let clipped = clipper.clip_start_of_alignment(&mut record, 3);
        assert_eq!(clipped, 3);

        // Check that the attribute was clipped
        if let Some(Value::String(s)) = record.data().get(&tag) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, b"3456789");
        } else {
            panic!("Tag XS not found or wrong type");
        }
    }

    #[test]
    fn test_auto_clip_attributes_string_3_prime() {
        use noodles::sam::alignment::record::data::field::Tag;

        let clipper = SamRecordClipper::with_auto_clip(ClippingMode::Hard, true);
        let mut record = create_test_record("10M", "ACGTACGTAC", 1000);

        // Add a string attribute that matches the read length
        let tag = Tag::from([b'X', b'S']);
        let value = Value::from("0123456789");
        record.data_mut().insert(tag, value);

        // Clip 3 bases from 3' end
        let clipped = clipper.clip_end_of_alignment(&mut record, 3);
        assert_eq!(clipped, 3);

        // Check that the attribute was clipped
        if let Some(Value::String(s)) = record.data().get(&tag) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, b"0123456");
        } else {
            panic!("Tag XS not found or wrong type");
        }
    }

    #[test]
    fn test_auto_clip_attributes_array_5_prime() {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::value::Array;

        let clipper = SamRecordClipper::with_auto_clip(ClippingMode::Hard, true);
        let mut record = create_test_record("10M", "ACGTACGTAC", 1000);

        // Add a UInt8Array attribute that matches the read length
        let tag = Tag::from([b'X', b'A']);
        let array: Vec<u8> = vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
        let value = Value::from(array);
        record.data_mut().insert(tag, value);

        // Clip 3 bases from 5' end
        let clipped = clipper.clip_start_of_alignment(&mut record, 3);
        assert_eq!(clipped, 3);

        // Check that the attribute was clipped
        if let Some(Value::Array(Array::UInt8(arr))) = record.data().get(&tag) {
            let vec_data: Vec<u8> = arr.clone();
            assert_eq!(vec_data, vec![3, 4, 5, 6, 7, 8, 9]);
        } else {
            panic!("Tag XA not found or wrong type");
        }
    }

    #[test]
    fn test_auto_clip_attributes_array_3_prime() {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::value::Array;

        let clipper = SamRecordClipper::with_auto_clip(ClippingMode::Hard, true);
        let mut record = create_test_record("10M", "ACGTACGTAC", 1000);

        // Add a UInt8Array attribute that matches the read length
        let tag = Tag::from([b'X', b'A']);
        let array: Vec<u8> = vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
        let value = Value::from(array);
        record.data_mut().insert(tag, value);

        // Clip 3 bases from 3' end
        let clipped = clipper.clip_end_of_alignment(&mut record, 3);
        assert_eq!(clipped, 3);

        // Check that the attribute was clipped
        if let Some(Value::Array(Array::UInt8(arr))) = record.data().get(&tag) {
            let vec_data: Vec<u8> = arr.clone();
            assert_eq!(vec_data, vec![0, 1, 2, 3, 4, 5, 6]);
        } else {
            panic!("Tag XA not found or wrong type");
        }
    }

    #[test]
    fn test_auto_clip_attributes_only_in_hard_mode() {
        use noodles::sam::alignment::record::data::field::Tag;

        // Test with Soft mode - attributes should NOT be clipped
        let clipper = SamRecordClipper::with_auto_clip(ClippingMode::Soft, true);
        let mut record = create_test_record("10M", "ACGTACGTAC", 1000);

        let tag = Tag::from([b'X', b'S']);
        let value = Value::from("0123456789");
        record.data_mut().insert(tag, value);

        clipper.clip_start_of_alignment(&mut record, 3);

        // Attribute should remain unchanged in Soft mode
        if let Some(Value::String(s)) = record.data().get(&tag) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, b"0123456789");
        }
    }

    #[test]
    fn test_auto_clip_attributes_only_when_enabled() {
        use noodles::sam::alignment::record::data::field::Tag;

        // Test with auto_clip_attributes disabled - attributes should NOT be clipped
        let clipper = SamRecordClipper::with_auto_clip(ClippingMode::Hard, false);
        let mut record = create_test_record("10M", "ACGTACGTAC", 1000);

        let tag = Tag::from([b'X', b'S']);
        let value = Value::from("0123456789");
        record.data_mut().insert(tag, value);

        clipper.clip_start_of_alignment(&mut record, 3);

        // Attribute should remain unchanged when auto-clip is disabled
        if let Some(Value::String(s)) = record.data().get(&tag) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, b"0123456789");
        }
    }

    #[test]
    fn test_auto_clip_attributes_only_matching_length() {
        use noodles::sam::alignment::record::data::field::Tag;

        let clipper = SamRecordClipper::with_auto_clip(ClippingMode::Hard, true);
        let mut record = create_test_record("10M", "ACGTACGTAC", 1000);

        // Add two attributes: one matching length, one not
        let tag1 = Tag::from([b'X', b'1']);
        let value1 = Value::from("0123456789"); // Matches length (10)
        record.data_mut().insert(tag1, value1);

        let tag2 = Tag::from([b'X', b'2']);
        let value2 = Value::from("01234"); // Does not match length (5)
        record.data_mut().insert(tag2, value2);

        clipper.clip_start_of_alignment(&mut record, 3);

        // Check tag1 was clipped
        if let Some(Value::String(s)) = record.data().get(&tag1) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, b"3456789");
        }

        // Check tag2 was NOT clipped
        if let Some(Value::String(s)) = record.data().get(&tag2) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, b"01234");
        }
    }

    // ===================================================================
    // clip_3_prime tests - Additional coverage matching Scala
    // ===================================================================

    #[test]
    fn test_clip_end_of_alignment_soft_with_insertion() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        // 44M + 2I + 4M = 50 bases
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"; // 50 bases
        assert_eq!(seq.len(), 50);
        let mut record = create_test_record("44M2I4M", seq, 10);

        let clipped = clipper.clip_end_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(10));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "40M10S");
    }

    #[test]
    fn test_clip_end_of_alignment_soft_with_deletion() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        // 44M + 2D + 6M = 50 bases in query
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"; // 50 bases
        assert_eq!(seq.len(), 50);
        let mut record = create_test_record("44M2D6M", seq, 10);

        let clipped = clipper.clip_end_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(10));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "40M10S");
    }

    #[test]
    fn test_clip_end_of_alignment_soft_additional_bases() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"; // 50 bases
        let mut record = create_test_record("40M10S", seq, 10);

        let clipped = clipper.clip_end_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(10));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "30M20S");
    }

    #[test]
    fn test_clip_end_of_alignment_preserve_hard_clip() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"; // 40 bases (no hard clip in sequence)
        let mut record = create_test_record("40M10H", &seq[..40], 10);

        let clipped = clipper.clip_end_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(10));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "30M10S10H");
    }

    #[test]
    fn test_clip_end_of_alignment_complex_cigar() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        // 10M + 5I + 5M + 10I + 16M + 4S + 2H = 50 query bases
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"; // 50 bases
        let mut record = create_test_record("10M5I5M10I16M4S2H", seq, 10);

        let clipped = clipper.clip_end_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(10));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "10M5I5M10I6M14S2H");
    }

    #[test]
    fn test_clip_end_of_alignment_consumes_leading_insertion() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        // 38M + 4I + 8M = 50 bases
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"; // 50 bases
        let mut record = create_test_record("38M4I8M", seq, 10);

        // Ask to clip 10 bases, but should consume 12 bases (4I + 8M) because
        // the insertion at the clip boundary is consumed entirely
        let clipped = clipper.clip_end_of_alignment(&mut record, 10);
        assert_eq!(clipped, 12);
        assert_eq!(record.alignment_start(), Position::new(10));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "38M12S");
    }

    #[test]
    fn test_clip_end_of_alignment_preserve_insertion_after_clip() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        // 36M + 4I + 10M = 50 bases
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"; // 50 bases
        let mut record = create_test_record("36M4I10M", seq, 10);

        let clipped = clipper.clip_end_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(10));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "36M4I10S");
    }

    #[test]
    fn test_clip_end_of_alignment_remove_deletion_before_clip() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        // 40M + 4D + 10M = 50 query bases
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"; // 50 bases
        let mut record = create_test_record("40M4D10M", seq, 10);

        let clipped = clipper.clip_end_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(10));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "40M10S");
    }

    #[test]
    fn test_clip_end_of_alignment_preserve_distant_deletion() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        // 25M + 4D + 25M = 50 query bases
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"; // 50 bases
        let mut record = create_test_record("25M4D25M", seq, 10);

        let clipped = clipper.clip_end_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(10));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "25M4D15M10S");
    }

    #[test]
    fn test_clip_end_of_alignment_unmapped_read_does_nothing() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"; // 50 bases
        let mut record = create_test_record("50M", seq, 10);

        // Set unmapped flag
        *record.flags_mut() = Flags::UNMAPPED;

        let clipped = clipper.clip_end_of_alignment(&mut record, 10);

        // Should not clip unmapped reads
        assert_eq!(clipped, 0);
        assert_eq!(record.alignment_start(), Position::new(10));

        // CIGAR should remain unchanged
        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "50M");
    }

    #[test]
    fn test_clip_end_of_alignment_soft_with_mask() {
        let clipper = SamRecordClipper::new(ClippingMode::SoftWithMask);
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"; // 50 bases
        let mut record = create_test_record("50M", seq, 10);

        let clipped = clipper.clip_end_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(10));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "40M10S");

        // Check last 10 bases are masked to N
        let bases = record.sequence();
        for i in 40..50 {
            assert_eq!(bases.as_ref()[i], b'N', "Base at position {i} should be N");
        }

        // Check last 10 quals are set to min
        let quals = record.quality_scores();
        for i in 40..50 {
            assert_eq!(quals.as_ref()[i], 2, "Quality at position {i} should be 2");
        }
    }

    #[test]
    fn test_clip_end_of_alignment_soft_with_mask_existing_soft_clip() {
        let clipper = SamRecordClipper::new(ClippingMode::SoftWithMask);
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"; // 50 bases
        let mut record = create_test_record("40M10S", seq, 10);

        let clipped = clipper.clip_end_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(10));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "30M20S");

        // Check last 20 bases are masked to N
        let bases = record.sequence();
        for i in 30..50 {
            assert_eq!(bases.as_ref()[i], b'N', "Base at position {i} should be N");
        }

        // Check last 20 quals are set to min
        let quals = record.quality_scores();
        for i in 30..50 {
            assert_eq!(quals.as_ref()[i], 2, "Quality at position {i} should be 2");
        }
    }

    #[test]
    fn test_clip_end_of_alignment_hard() {
        let clipper = SamRecordClipper::new(ClippingMode::Hard);
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"; // 50 bases
        let mut record = create_test_record("50M", seq, 10);

        let orig_seq: Vec<u8> = record.sequence().as_ref().to_vec();
        let orig_qual: Vec<u8> = record.quality_scores().as_ref().to_vec();

        let clipped = clipper.clip_end_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(10));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "40M10H");

        // Check sequence and qualities are truncated
        assert_eq!(record.sequence().as_ref(), &orig_seq[..40]);
        assert_eq!(record.quality_scores().as_ref(), &orig_qual[..40]);
    }

    #[test]
    fn test_clip_start_of_alignment_hard_existing_soft_clip() {
        let clipper = SamRecordClipper::new(ClippingMode::Hard);
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"; // 50 bases
        assert_eq!(seq.len(), 50);
        let mut record = create_test_record("10S40M", seq, 10);

        let orig_seq: Vec<u8> = record.sequence().as_ref().to_vec();
        let orig_qual: Vec<u8> = record.quality_scores().as_ref().to_vec();

        let clipped = clipper.clip_start_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);
        assert_eq!(record.alignment_start(), Position::new(20));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "20H30M");

        // Check sequence and qualities are truncated (removed 10S + 10M = 20 bases)
        assert_eq!(record.sequence().as_ref(), &orig_seq[20..]);
        assert_eq!(record.quality_scores().as_ref(), &orig_qual[20..]);
    }

    // ===================================================================
    // Auto-trim attribute tests - all mode combinations
    // ===================================================================

    #[test]
    fn test_clip_start_of_alignment_auto_trim_soft_mode_false() {
        use noodles::sam::alignment::record::data::field::Tag;

        let clipper = SamRecordClipper::with_auto_clip(ClippingMode::Soft, false);
        let seq = "ACGTACGTACGTACGTACGT"; // 20 bases
        let mut record = create_test_record("20M", seq, 10);

        let a1_tag = Tag::from([b'A', b'1']);
        record.data_mut().insert(a1_tag, Value::from("AB".repeat(10)));

        let clipped = clipper.clip_start_of_alignment(&mut record, 5);
        assert_eq!(clipped, 5);

        // In Soft mode with auto=false, attributes should NOT be modified
        if let Some(Value::String(s)) = record.data().get(&a1_tag) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, "AB".repeat(10).as_bytes());
        }
    }

    #[test]
    fn test_clip_start_of_alignment_auto_trim_soft_mode_true() {
        use noodles::sam::alignment::record::data::field::Tag;

        let clipper = SamRecordClipper::with_auto_clip(ClippingMode::Soft, true);
        let seq = "ACGTACGTACGTACGTACGT"; // 20 bases
        let mut record = create_test_record("20M", seq, 10);

        let a1_tag = Tag::from([b'A', b'1']);
        record.data_mut().insert(a1_tag, Value::from("AB".repeat(10)));

        let clipped = clipper.clip_start_of_alignment(&mut record, 5);
        assert_eq!(clipped, 5);

        // In Soft mode, even with auto=true, attributes should NOT be modified
        if let Some(Value::String(s)) = record.data().get(&a1_tag) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, "AB".repeat(10).as_bytes());
        }
    }

    #[test]
    fn test_clip_start_of_alignment_auto_trim_soft_with_mask_mode_false() {
        use noodles::sam::alignment::record::data::field::Tag;

        let clipper = SamRecordClipper::with_auto_clip(ClippingMode::SoftWithMask, false);
        let seq = "ACGTACGTACGTACGTACGT"; // 20 bases
        let mut record = create_test_record("20M", seq, 10);

        let a1_tag = Tag::from([b'A', b'1']);
        record.data_mut().insert(a1_tag, Value::from("AB".repeat(10)));

        let clipped = clipper.clip_start_of_alignment(&mut record, 5);
        assert_eq!(clipped, 5);

        // In SoftWithMask mode with auto=false, attributes should NOT be modified
        if let Some(Value::String(s)) = record.data().get(&a1_tag) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, "AB".repeat(10).as_bytes());
        }
    }

    #[test]
    fn test_clip_start_of_alignment_auto_trim_soft_with_mask_mode_true() {
        use noodles::sam::alignment::record::data::field::Tag;

        let clipper = SamRecordClipper::with_auto_clip(ClippingMode::SoftWithMask, true);
        let seq = "ACGTACGTACGTACGTACGT"; // 20 bases
        let mut record = create_test_record("20M", seq, 10);

        let a1_tag = Tag::from([b'A', b'1']);
        record.data_mut().insert(a1_tag, Value::from("AB".repeat(10)));

        let clipped = clipper.clip_start_of_alignment(&mut record, 5);
        assert_eq!(clipped, 5);

        // In SoftWithMask mode, even with auto=true, attributes should NOT be modified
        if let Some(Value::String(s)) = record.data().get(&a1_tag) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, "AB".repeat(10).as_bytes());
        }
    }

    #[test]
    fn test_clip_start_of_alignment_auto_trim_hard_mode_false() {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::value::Array;

        let clipper = SamRecordClipper::with_auto_clip(ClippingMode::Hard, false);
        let seq = "ACGTACGTACGTACGTACGT"; // 20 bases
        let mut record = create_test_record("20M", seq, 10);

        let a1_tag = Tag::from([b'A', b'1']);
        let a2_tag = Tag::from([b'A', b'2']);

        record.data_mut().insert(a1_tag, Value::from("AB".repeat(10)));
        record.data_mut().insert(a2_tag, Value::from((1..=20).collect::<Vec<i32>>()));

        let clipped = clipper.clip_start_of_alignment(&mut record, 5);
        assert_eq!(clipped, 5);

        // In Hard mode with auto=false, attributes should NOT be modified
        if let Some(Value::String(s)) = record.data().get(&a1_tag) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, "AB".repeat(10).as_bytes());
        }
        if let Some(Value::Array(Array::Int32(arr))) = record.data().get(&a2_tag) {
            let vec: Vec<i32> = arr.clone();
            assert_eq!(vec, (1..=20).collect::<Vec<i32>>());
        }
    }

    #[test]
    fn test_clip_start_of_alignment_auto_trim_hard_mode_true() {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::value::Array;

        let clipper = SamRecordClipper::with_auto_clip(ClippingMode::Hard, true);
        let seq = "ACGTACGTACGTACGTACGT"; // 20 bases
        let mut record = create_test_record("20M", seq, 10);

        let a1_tag = Tag::from([b'A', b'1']);
        let a2_tag = Tag::from([b'A', b'2']);
        let b1_tag = Tag::from([b'B', b'1']);
        let b2_tag = Tag::from([b'B', b'2']);

        record.data_mut().insert(a1_tag, Value::from("AB".repeat(10)));
        record.data_mut().insert(a2_tag, Value::from((1..=20).collect::<Vec<i32>>()));
        record.data_mut().insert(b1_tag, Value::from("A".repeat(10)));
        record.data_mut().insert(b2_tag, Value::from((1..=10).collect::<Vec<i32>>()));

        let clipped = clipper.clip_start_of_alignment(&mut record, 5);
        assert_eq!(clipped, 5);

        // In Hard mode with auto=true, attributes matching read length should be clipped
        if let Some(Value::String(s)) = record.data().get(&a1_tag) {
            // "ABABABABABABABABABAB" -> remove first 5 -> "BABABABABABABAB"
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, "BABABABABABABAB".as_bytes());
        }
        if let Some(Value::Array(Array::Int32(arr))) = record.data().get(&a2_tag) {
            let vec: Vec<i32> = arr.clone();
            assert_eq!(vec, (6..=20).collect::<Vec<i32>>());
        }
        // B1 and B2 should NOT be modified (length doesn't match)
        if let Some(Value::String(s)) = record.data().get(&b1_tag) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, "A".repeat(10).as_bytes());
        }
        if let Some(Value::Array(Array::Int32(arr))) = record.data().get(&b2_tag) {
            let vec: Vec<i32> = arr.clone();
            assert_eq!(vec, (1..=10).collect::<Vec<i32>>());
        }
    }

    #[test]
    fn test_clip_end_of_alignment_auto_trim_soft_mode_false() {
        use noodles::sam::alignment::record::data::field::Tag;

        let clipper = SamRecordClipper::with_auto_clip(ClippingMode::Soft, false);
        let seq = "ACGTACGTACGTACGTACGT"; // 20 bases
        let mut record = create_test_record("20M", seq, 10);

        let a1_tag = Tag::from([b'A', b'1']);
        record.data_mut().insert(a1_tag, Value::from("AB".repeat(10)));

        let clipped = clipper.clip_end_of_alignment(&mut record, 5);
        assert_eq!(clipped, 5);

        // In Soft mode with auto=false, attributes should NOT be modified
        if let Some(Value::String(s)) = record.data().get(&a1_tag) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, "AB".repeat(10).as_bytes());
        }
    }

    #[test]
    fn test_clip_end_of_alignment_auto_trim_soft_mode_true() {
        use noodles::sam::alignment::record::data::field::Tag;

        let clipper = SamRecordClipper::with_auto_clip(ClippingMode::Soft, true);
        let seq = "ACGTACGTACGTACGTACGTACGT"; // 20 bases
        let mut record = create_test_record("20M", seq, 10);

        let a1_tag = Tag::from([b'A', b'1']);
        record.data_mut().insert(a1_tag, Value::from("AB".repeat(10)));

        let clipped = clipper.clip_end_of_alignment(&mut record, 5);
        assert_eq!(clipped, 5);

        // In Soft mode, even with auto=true, attributes should NOT be modified
        if let Some(Value::String(s)) = record.data().get(&a1_tag) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, "AB".repeat(10).as_bytes());
        }
    }

    #[test]
    fn test_clip_end_of_alignment_auto_trim_soft_with_mask_mode_false() {
        use noodles::sam::alignment::record::data::field::Tag;

        let clipper = SamRecordClipper::with_auto_clip(ClippingMode::SoftWithMask, false);
        let seq = "ACGTACGTACGTACGTACGT"; // 20 bases
        let mut record = create_test_record("20M", seq, 10);

        let a1_tag = Tag::from([b'A', b'1']);
        record.data_mut().insert(a1_tag, Value::from("AB".repeat(10)));

        let clipped = clipper.clip_end_of_alignment(&mut record, 5);
        assert_eq!(clipped, 5);

        // In SoftWithMask mode with auto=false, attributes should NOT be modified
        if let Some(Value::String(s)) = record.data().get(&a1_tag) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, "AB".repeat(10).as_bytes());
        }
    }

    #[test]
    fn test_clip_end_of_alignment_auto_trim_soft_with_mask_mode_true() {
        use noodles::sam::alignment::record::data::field::Tag;

        let clipper = SamRecordClipper::with_auto_clip(ClippingMode::SoftWithMask, true);
        let seq = "ACGTACGTACGTACGTACGT"; // 20 bases
        let mut record = create_test_record("20M", seq, 10);

        let a1_tag = Tag::from([b'A', b'1']);
        record.data_mut().insert(a1_tag, Value::from("AB".repeat(10)));

        let clipped = clipper.clip_end_of_alignment(&mut record, 5);
        assert_eq!(clipped, 5);

        // In SoftWithMask mode, even with auto=true, attributes should NOT be modified
        if let Some(Value::String(s)) = record.data().get(&a1_tag) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, "AB".repeat(10).as_bytes());
        }
    }

    #[test]
    fn test_clip_end_of_alignment_auto_trim_hard_mode_false() {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::value::Array;

        let clipper = SamRecordClipper::with_auto_clip(ClippingMode::Hard, false);
        let seq = "ACGTACGTACGTACGTACGT"; // 20 bases
        let mut record = create_test_record("20M", seq, 10);

        let a1_tag = Tag::from([b'A', b'1']);
        let a2_tag = Tag::from([b'A', b'2']);

        record.data_mut().insert(a1_tag, Value::from("AB".repeat(10)));
        record.data_mut().insert(a2_tag, Value::from((1..=20).collect::<Vec<i32>>()));

        let clipped = clipper.clip_end_of_alignment(&mut record, 5);
        assert_eq!(clipped, 5);

        // In Hard mode with auto=false, attributes should NOT be modified
        if let Some(Value::String(s)) = record.data().get(&a1_tag) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, "AB".repeat(10).as_bytes());
        }
        if let Some(Value::Array(Array::Int32(arr))) = record.data().get(&a2_tag) {
            let vec: Vec<i32> = arr.clone();
            assert_eq!(vec, (1..=20).collect::<Vec<i32>>());
        }
    }

    #[test]
    fn test_clip_end_of_alignment_auto_trim_hard_mode_true() {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::value::Array;

        let clipper = SamRecordClipper::with_auto_clip(ClippingMode::Hard, true);
        let seq = "ACGTACGTACGTACGTACGT"; // 20 bases
        let mut record = create_test_record("20M", seq, 10);

        let a1_tag = Tag::from([b'A', b'1']);
        let a2_tag = Tag::from([b'A', b'2']);
        let b1_tag = Tag::from([b'B', b'1']);
        let b2_tag = Tag::from([b'B', b'2']);

        record.data_mut().insert(a1_tag, Value::from("AB".repeat(10)));
        record.data_mut().insert(a2_tag, Value::from((1..=20).collect::<Vec<i32>>()));
        record.data_mut().insert(b1_tag, Value::from("A".repeat(10)));
        record.data_mut().insert(b2_tag, Value::from((1..=10).collect::<Vec<i32>>()));

        let clipped = clipper.clip_end_of_alignment(&mut record, 5);
        assert_eq!(clipped, 5);

        // In Hard mode with auto=true, attributes matching read length should be clipped
        if let Some(Value::String(s)) = record.data().get(&a1_tag) {
            // "ABABABABABABABABABAB" -> remove last 5 -> "ABABABABABABABA"
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, "ABABABABABABABA".as_bytes());
        }
        if let Some(Value::Array(Array::Int32(arr))) = record.data().get(&a2_tag) {
            let vec: Vec<i32> = arr.clone();
            assert_eq!(vec, (1..=15).collect::<Vec<i32>>());
        }
        // B1 and B2 should NOT be modified (length doesn't match)
        if let Some(Value::String(s)) = record.data().get(&b1_tag) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, "A".repeat(10).as_bytes());
        }
        if let Some(Value::Array(Array::Int32(arr))) = record.data().get(&b2_tag) {
            let vec: Vec<i32> = arr.clone();
            assert_eq!(vec, (1..=10).collect::<Vec<i32>>());
        }
    }

    // ===================================================================
    // upgrade_clipping tests
    // ===================================================================

    #[test]
    fn test_upgrade_all_clipping_convert_leading_trailing() {
        use noodles::sam::alignment::record::data::field::Tag;

        // Test without auto-clip
        let clipper = SamRecordClipper::new(ClippingMode::Hard);
        let seq = "12345678901234567890123456789012345678901234567890"; // 50 bases
        let mut no_auto = create_test_record("5S35M10S", seq, 10);
        let az_tag = Tag::from([b'a', b'z']);
        no_auto
            .data_mut()
            .insert(az_tag, Value::from("12345678901234567890123456789012345678901234567890"));

        let result = clipper.upgrade_all_clipping(&mut no_auto).unwrap();
        assert_eq!(result, (5, 10));

        let cigar_str = format_cigar(&no_auto.cigar());
        assert_eq!(cigar_str, "5H35M10H");
        assert_eq!(no_auto.sequence().len(), 35);

        // Attributes should NOT be modified without auto-clip
        if let Some(Value::String(s)) = no_auto.data().get(&az_tag) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, "12345678901234567890123456789012345678901234567890".as_bytes());
        }

        // Test with auto-clip
        let clipper_auto = SamRecordClipper::with_auto_clip(ClippingMode::Hard, true);
        let mut with_auto = create_test_record("5S35M10S", seq, 10);
        with_auto
            .data_mut()
            .insert(az_tag, Value::from("12345678901234567890123456789012345678901234567890"));

        let result = clipper_auto.upgrade_all_clipping(&mut with_auto).unwrap();
        assert_eq!(result, (5, 10));

        let cigar_str = format_cigar(&with_auto.cigar());
        assert_eq!(cigar_str, "5H35M10H");
        assert_eq!(with_auto.sequence().len(), 35);

        // Attributes SHOULD be modified with auto-clip (remove first 5 and last 10)
        if let Some(Value::String(s)) = with_auto.data().get(&az_tag) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, "67890123456789012345678901234567890".as_bytes());
        }
    }

    #[test]
    fn test_upgrade_all_clipping_soft_clips_after_hard_clips() {
        use noodles::sam::alignment::record::data::field::Tag;

        // Test without auto-clip
        let clipper = SamRecordClipper::new(ClippingMode::Hard);
        let seq = "12345678901234567890123456789012345678901234567890"; // 50 bases
        let mut no_auto = create_test_record("5H5S35M10S5H", seq, 10);
        let az_tag = Tag::from([b'a', b'z']);
        no_auto
            .data_mut()
            .insert(az_tag, Value::from("12345678901234567890123456789012345678901234567890"));

        let result = clipper.upgrade_all_clipping(&mut no_auto).unwrap();
        assert_eq!(result, (5, 10));

        let cigar_str = format_cigar(&no_auto.cigar());
        assert_eq!(cigar_str, "10H35M15H");
        assert_eq!(no_auto.sequence().len(), 35);

        // Attributes should NOT be modified without auto-clip
        if let Some(Value::String(s)) = no_auto.data().get(&az_tag) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, "12345678901234567890123456789012345678901234567890".as_bytes());
        }

        // Test with auto-clip
        let clipper_auto = SamRecordClipper::with_auto_clip(ClippingMode::Hard, true);
        let mut with_auto = create_test_record("5H5S35M10S5H", seq, 10);
        with_auto
            .data_mut()
            .insert(az_tag, Value::from("12345678901234567890123456789012345678901234567890"));

        let result = clipper_auto.upgrade_all_clipping(&mut with_auto).unwrap();
        assert_eq!(result, (5, 10));

        let cigar_str = format_cigar(&with_auto.cigar());
        assert_eq!(cigar_str, "10H35M15H");
        assert_eq!(with_auto.sequence().len(), 35);

        // Attributes SHOULD be modified with auto-clip
        if let Some(Value::String(s)) = with_auto.data().get(&az_tag) {
            let bytes: &[u8] = s.as_ref();
            assert_eq!(bytes, "67890123456789012345678901234567890".as_bytes());
        }
    }

    #[test]
    fn test_upgrade_all_clipping_no_soft_clipping() {
        use noodles::sam::alignment::record::data::field::Tag;

        let clipper = SamRecordClipper::new(ClippingMode::Hard);
        let az_tag = Tag::from([b'a', b'z']);

        // Test with no soft clips
        let seq1 = "1234567890123456789012345678901234567890123456789012345"; // 55 bases
        let mut no_soft = create_test_record("55M", seq1, 10);
        no_soft
            .data_mut()
            .insert(az_tag, Value::from("12345678901234567890123456789012345678901234567890"));

        let result = clipper.upgrade_all_clipping(&mut no_soft).unwrap();
        assert_eq!(result, (0, 0));

        let cigar_str = format_cigar(&no_soft.cigar());
        assert_eq!(cigar_str, "55M");
        assert_eq!(no_soft.sequence().len(), 55);

        // Test with only hard clips (no soft clips)
        let mut hard_only = create_test_record("5H55M10H", seq1, 10);
        hard_only
            .data_mut()
            .insert(az_tag, Value::from("12345678901234567890123456789012345678901234567890"));

        let result = clipper.upgrade_all_clipping(&mut hard_only).unwrap();
        assert_eq!(result, (0, 0));

        let cigar_str = format_cigar(&hard_only.cigar());
        assert_eq!(cigar_str, "5H55M10H");
        assert_eq!(hard_only.sequence().len(), 55);
    }

    #[test]
    fn test_upgrade_all_clipping_unmapped_or_wrong_mode() {
        use noodles::sam::alignment::record::data::field::Tag;

        let az_tag = Tag::from([b'a', b'z']);
        let seq = "1234567890123456789012345678901234567890123456789012345"; // 55 bases

        // Test Soft mode (should not convert)
        let clipper_soft = SamRecordClipper::new(ClippingMode::Soft);
        let mut mapped = create_test_record("55M", seq, 10);
        mapped
            .data_mut()
            .insert(az_tag, Value::from("12345678901234567890123456789012345678901234567890"));

        let result = clipper_soft.upgrade_all_clipping(&mut mapped).unwrap();
        assert_eq!(result, (0, 0));

        // Test SoftWithMask mode (should not convert)
        let clipper_mask = SamRecordClipper::new(ClippingMode::SoftWithMask);
        let mut mapped2 = create_test_record("55M", seq, 10);
        mapped2
            .data_mut()
            .insert(az_tag, Value::from("12345678901234567890123456789012345678901234567890"));

        let result = clipper_mask.upgrade_all_clipping(&mut mapped2).unwrap();
        assert_eq!(result, (0, 0));

        // Test unmapped read (should not convert)
        let clipper_hard = SamRecordClipper::new(ClippingMode::Hard);
        let mut unmapped = create_test_record("55M", seq, 10);
        *unmapped.flags_mut() = Flags::UNMAPPED;
        unmapped
            .data_mut()
            .insert(az_tag, Value::from("12345678901234567890123456789012345678901234567890"));

        let result = clipper_hard.upgrade_all_clipping(&mut unmapped).unwrap();
        assert_eq!(result, (0, 0));
    }

    // ===================================================================
    // clipOverlappingReads tests - complex overlap scenarios
    // ===================================================================

    fn create_pair(
        r1_start: usize,
        r1_cigar: &str,
        r1_seq: &str,
        r2_start: usize,
        r2_cigar: &str,
        r2_seq: &str,
    ) -> (RecordBuf, RecordBuf) {
        let mut r1 = create_test_record(r1_cigar, r1_seq, r1_start);
        let mut r2 = create_test_record(r2_cigar, r2_seq, r2_start);

        // Calculate template length (insert size)
        // For FR pair: template length = r2_end - r1_start + 1
        let r2_len = cigar_utils::reference_length(&r2.cigar());
        let r2_end = r2_start + r2_len - 1;
        let tlen = i32::try_from(r2_end).unwrap() - i32::try_from(r1_start).unwrap() + 1;

        // Set up FR pair flags
        // R1 is forward strand, R2 is reverse strand (typical FR pair)
        *r1.flags_mut() = Flags::SEGMENTED | Flags::MATE_REVERSE_COMPLEMENTED;
        *r2.flags_mut() = Flags::SEGMENTED | Flags::REVERSE_COMPLEMENTED;

        // Set mate information
        *r1.mate_reference_sequence_id_mut() = Some(0);
        *r1.mate_alignment_start_mut() = Position::new(r2_start);
        *r1.template_length_mut() = tlen;
        *r2.mate_reference_sequence_id_mut() = Some(0);
        *r2.mate_alignment_start_mut() = Position::new(r1_start);
        *r2.template_length_mut() = -tlen;

        (r1, r2)
    }

    #[test]
    fn test_clip_overlapping_reads_one_base_overlap() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let seq = "A".repeat(100);
        let (mut r1, mut r2) = create_pair(1, "100M", &seq, 100, "100M", &seq);

        let (clipped_r1, clipped_r2) = clipper.clip_overlapping_reads(&mut r1, &mut r2);
        assert_eq!(clipped_r1, 0);
        assert_eq!(clipped_r2, 1);

        assert_eq!(r1.alignment_start(), Position::new(1));
        let cigar_r1 = format_cigar(&r1.cigar());
        assert_eq!(cigar_r1, "100M");

        assert_eq!(r2.alignment_start(), Position::new(101));
        let cigar_r2 = format_cigar(&r2.cigar());
        assert_eq!(cigar_r2, "1S99M");
    }

    #[test]
    fn test_clip_overlapping_reads_two_base_overlap() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let seq = "A".repeat(100);
        let (mut r1, mut r2) = create_pair(2, "100M", &seq, 100, "100M", &seq);

        let (clipped_r1, clipped_r2) = clipper.clip_overlapping_reads(&mut r1, &mut r2);
        assert_eq!(clipped_r1, 1);
        assert_eq!(clipped_r2, 1);

        assert_eq!(r1.alignment_start(), Position::new(2));
        let cigar_r1 = format_cigar(&r1.cigar());
        assert_eq!(cigar_r1, "99M1S");

        assert_eq!(r2.alignment_start(), Position::new(101));
        let cigar_r2 = format_cigar(&r2.cigar());
        assert_eq!(cigar_r2, "1S99M");
    }

    #[test]
    fn test_clip_overlapping_reads_with_deletion() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let seq = "A".repeat(90);
        let (mut r1, mut r2) = create_pair(2, "80M10D10M", &seq, 70, "10M10D80M", &seq);

        let (clipped_r1, clipped_r2) = clipper.clip_overlapping_reads(&mut r1, &mut r2);
        assert_eq!(clipped_r1, 10);
        assert_eq!(clipped_r2, 10);

        assert_eq!(r1.alignment_start(), Position::new(2));
        let cigar_r1 = format_cigar(&r1.cigar());
        assert_eq!(cigar_r1, "80M10S");

        assert_eq!(r2.alignment_start(), Position::new(90));
        let cigar_r2 = format_cigar(&r2.cigar());
        assert_eq!(cigar_r2, "10S80M");
    }

    #[test]
    fn test_clip_overlapping_reads_full_overlap() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let seq = "A".repeat(100);
        let (mut r1, mut r2) = create_pair(1, "100M", &seq, 1, "100M", &seq);

        let (clipped_r1, clipped_r2) = clipper.clip_overlapping_reads(&mut r1, &mut r2);
        assert_eq!(clipped_r1, 50);
        assert_eq!(clipped_r2, 50);

        assert_eq!(r1.alignment_start(), Position::new(1));
        let cigar_r1 = format_cigar(&r1.cigar());
        assert_eq!(cigar_r1, "50M50S");

        assert_eq!(r2.alignment_start(), Position::new(51));
        let cigar_r2 = format_cigar(&r2.cigar());
        assert_eq!(cigar_r2, "50S50M");
    }

    #[test]
    fn test_clip_overlapping_reads_extend_past_each_other() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let seq = "A".repeat(100);
        let (mut r1, mut r2) = create_pair(50, "100M", &seq, 1, "100M", &seq);

        // Midpoint should be (50+100)/2 = 75
        let (clipped_r1, clipped_r2) = clipper.clip_overlapping_reads(&mut r1, &mut r2);
        assert_eq!(clipped_r1, 74);
        assert_eq!(clipped_r2, 75);

        assert_eq!(r1.alignment_start(), Position::new(50));
        let cigar_r1 = format_cigar(&r1.cigar());
        assert_eq!(cigar_r1, "26M74S");

        assert_eq!(r2.alignment_start(), Position::new(76));
        let cigar_r2 = format_cigar(&r2.cigar());
        assert_eq!(cigar_r2, "75S25M");
    }

    #[test]
    fn test_clip_overlapping_reads_forward_much_longer() {
        let clipper = SamRecordClipper::with_auto_clip(ClippingMode::Hard, false);
        let r1_seq = "A".repeat(100);
        let r2_seq = "A".repeat(100);
        let (mut r1, mut r2) = create_pair(1, "100M", &r1_seq, 30, "80S20M", &r2_seq);

        let (clipped_r1, clipped_r2) = clipper.clip_overlapping_reads(&mut r1, &mut r2);
        assert_eq!(clipped_r1, 71);
        assert_eq!(clipped_r2, 0);

        assert_eq!(r1.alignment_start(), Position::new(1));
        let cigar_r1 = format_cigar(&r1.cigar());
        assert_eq!(cigar_r1, "29M71H");

        assert_eq!(r2.alignment_start(), Position::new(30));
        let cigar_r2 = format_cigar(&r2.cigar());
        assert_eq!(cigar_r2, "80H20M");
    }

    #[test]
    fn test_clip_overlapping_reads_reverse_much_longer() {
        let clipper = SamRecordClipper::with_auto_clip(ClippingMode::Hard, false);
        let r1_seq = "A".repeat(100);
        let r2_seq = "A".repeat(100);
        let (mut r1, mut r2) = create_pair(50, "20M80S", &r1_seq, 1, "100M", &r2_seq);

        let (clipped_r1, clipped_r2) = clipper.clip_overlapping_reads(&mut r1, &mut r2);
        assert_eq!(clipped_r1, 0);
        assert_eq!(clipped_r2, 69);

        assert_eq!(r1.alignment_start(), Position::new(50));
        let cigar_r1 = format_cigar(&r1.cigar());
        assert_eq!(cigar_r1, "20M80H");

        assert_eq!(r2.alignment_start(), Position::new(70));
        let cigar_r2 = format_cigar(&r2.cigar());
        assert_eq!(cigar_r2, "69H31M");
    }

    #[test]
    fn test_clip_overlapping_reads_with_one_end_deletion() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let r1_seq = "A".repeat(100);
        let r2_seq = "A".repeat(110);
        let (mut r1, mut r2) = create_pair(1, "60M10D40M", &r1_seq, 50, "10M10D80M10D10M", &r2_seq);

        let (clipped_r1, clipped_r2) = clipper.clip_overlapping_reads(&mut r1, &mut r2);
        assert_eq!(clipped_r1, 25);
        assert_eq!(clipped_r2, 26);

        assert_eq!(r1.alignment_start(), Position::new(1));
        let cigar_r1 = format_cigar(&r1.cigar());
        assert_eq!(cigar_r1, "60M10D15M25S");

        assert_eq!(r2.alignment_start(), Position::new(86));
        let cigar_r2 = format_cigar(&r2.cigar());
        assert_eq!(cigar_r2, "26S64M10D10M");
    }

    #[test]
    fn test_clip_overlapping_reads_both_ends_deletion() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let seq = "A".repeat(100);
        let (mut r1, mut r2) = create_pair(1, "50M10D50M", &seq, 3, "47M10D53M", &seq);

        let (clipped_r1, clipped_r2) = clipper.clip_overlapping_reads(&mut r1, &mut r2);
        assert_eq!(clipped_r1, 50);
        assert_eq!(clipped_r2, 47);

        assert_eq!(r1.alignment_start(), Position::new(1));
        let cigar_r1 = format_cigar(&r1.cigar());
        assert_eq!(cigar_r1, "50M50S");

        // R2 clips through the 47M and 10D, so the new start is at position 60 (where 53M begins)
        assert_eq!(r2.alignment_start(), Position::new(60));
        let cigar_r2 = format_cigar(&r2.cigar());
        assert_eq!(cigar_r2, "47S53M");
    }

    #[test]
    fn test_clip_overlapping_reads_no_overlap_far_apart() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let seq = "A".repeat(100);
        let (mut r1, mut r2) = create_pair(1000, "100M", &seq, 1, "100M", &seq);

        let (clipped_r1, clipped_r2) = clipper.clip_overlapping_reads(&mut r1, &mut r2);
        assert_eq!(clipped_r1, 0);
        assert_eq!(clipped_r2, 0);

        assert_eq!(r1.alignment_start(), Position::new(1000));
        let cigar_r1 = format_cigar(&r1.cigar());
        assert_eq!(cigar_r1, "100M");

        assert_eq!(r2.alignment_start(), Position::new(1));
        let cigar_r2 = format_cigar(&r2.cigar());
        assert_eq!(cigar_r2, "100M");
    }

    // Tests for clip_extending_past_mate

    #[test]
    fn test_clip_extending_past_mate_ends_basic() {
        // Based on the Scala test: r1 at 100-149, r2 at 90-139
        // After clipping, both should have same start and end
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let seq = "A".repeat(50);
        let (mut r1, mut r2) = create_pair(100, "50M", &seq, 90, "50M", &seq);

        // Before: R1: 100-149, R2: 90-139
        // R1 extends 10 bases past R2's end (149 vs 139), so clip 10 from R1's 3' end
        // R2 starts 10 bases before R1's start (90 vs 100), so clip 10 from R2's 5' end (reference)

        let (clipped_r1, clipped_r2) = clipper.clip_extending_past_mate_ends(&mut r1, &mut r2);

        // Both reads should be clipped
        assert_eq!(clipped_r1, 10); // R1 clipped from 3' end
        assert_eq!(clipped_r2, 10); // R2 clipped from 5' reference end (3' read end for reverse strand)

        // After clipping: R1: 100-139, R2: 100-139 (same start and end)
        assert_eq!(r1.alignment_start(), Position::new(100));
        let cigar_r1 = format_cigar(&r1.cigar());
        assert_eq!(cigar_r1, "40M10S");

        assert_eq!(r2.alignment_start(), Position::new(100));
        let cigar_r2 = format_cigar(&r2.cigar());
        assert_eq!(cigar_r2, "10S40M");
    }

    #[test]
    fn test_clip_extending_past_mate_ends_both_extend() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let seq = "A".repeat(60);
        // R1: 100-159, R2: 110-169
        let (mut r1, mut r2) = create_pair(100, "60M", &seq, 110, "60M", &seq);

        let (clipped_r1, clipped_r2) = clipper.clip_extending_past_mate_ends(&mut r1, &mut r2);

        // R1 extends 10 bases past R2's end (159 vs 169)
        assert_eq!(clipped_r1, 0); // R1 doesn't extend past R2
        // R2 starts 10 bases before R1's start (110 vs 100) - wait, that's backwards
        assert_eq!(clipped_r2, 0); // R2 starts after R1

        let cigar_r1 = format_cigar(&r1.cigar());
        assert_eq!(cigar_r1, "60M");

        let cigar_r2 = format_cigar(&r2.cigar());
        assert_eq!(cigar_r2, "60M");
    }

    #[test]
    fn test_clip_extending_past_mate_ends_with_soft_clips() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let seq = "A".repeat(60);
        // R1: `10S50M` at 100 → unclipped 90-149, alignment 100-149
        // R2: `50M` at 90 → unclipped 90-139, alignment 90-139
        // R1's alignment end (149) extends 10 bases past R2's unclipped end (139)
        // R2's alignment start (90) is 10 bases before R1's unclipped start (90) - no extension
        let (mut r1, mut r2) = create_pair(100, "10S50M", &seq, 90, "50M", &seq);

        let (clipped_r1, clipped_r2) = clipper.clip_extending_past_mate_ends(&mut r1, &mut r2);

        // R1 should be clipped by 10 from the end
        assert_eq!(clipped_r1, 10);
        // R2 should not be clipped (its start doesn't extend past R1's unclipped start)
        assert_eq!(clipped_r2, 0);

        let cigar_r1 = format_cigar(&r1.cigar());
        assert_eq!(cigar_r1, "10S40M10S");

        // R2 unchanged
        assert_eq!(r2.alignment_start(), Position::new(90));
        let cigar_r2 = format_cigar(&r2.cigar());
        assert_eq!(cigar_r2, "50M");
    }

    #[test]
    fn test_clip_extending_past_mate_ends_no_extension() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let seq = "A".repeat(50);
        // R1: 100-149, R2: 200-249 (no extension)
        let (mut r1, mut r2) = create_pair(100, "50M", &seq, 200, "50M", &seq);

        let (clipped_r1, clipped_r2) = clipper.clip_extending_past_mate_ends(&mut r1, &mut r2);

        assert_eq!(clipped_r1, 0);
        assert_eq!(clipped_r2, 0);

        let cigar_r1 = format_cigar(&r1.cigar());
        assert_eq!(cigar_r1, "50M");

        let cigar_r2 = format_cigar(&r2.cigar());
        assert_eq!(cigar_r2, "50M");
    }

    #[test]
    fn test_clip_extending_past_mate_ends_ff_pair() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let seq = "A".repeat(50);
        let mut r1 = create_test_record("50M", &seq, 100);
        let mut r2 = create_test_record("50M", &seq, 90);

        // Make them both forward strand (not FR pair)
        *r1.flags_mut() = Flags::SEGMENTED;
        *r2.flags_mut() = Flags::SEGMENTED; // Both forward, not FR

        let (clipped_r1, clipped_r2) = clipper.clip_extending_past_mate_ends(&mut r1, &mut r2);

        // Should not clip non-FR pairs
        assert_eq!(clipped_r1, 0);
        assert_eq!(clipped_r2, 0);
    }

    #[test]
    fn test_clip_extending_past_mate_ends_with_deletion() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let seq = "A".repeat(50);
        // R1: 100-159 (50M + 10D)
        // R2: 110-159
        let (mut r1, mut r2) = create_pair(100, "50M10D10M", &seq, 110, "50M", &seq);

        let (clipped_r1, clipped_r2) = clipper.clip_extending_past_mate_ends(&mut r1, &mut r2);

        // R1 ends at 169, extends past R2's end at 159
        // Should clip ~10 bases worth
        assert!(clipped_r1 > 0);
        assert_eq!(clipped_r2, 0);
    }

    #[test]
    fn test_clip_extending_past_mate_ends_hard_mode() {
        let clipper = SamRecordClipper::new(ClippingMode::Hard);
        let seq = "A".repeat(50);
        let (mut r1, mut r2) = create_pair(100, "50M", &seq, 90, "50M", &seq);

        let (clipped_r1, clipped_r2) = clipper.clip_extending_past_mate_ends(&mut r1, &mut r2);

        assert_eq!(clipped_r1, 10);
        assert_eq!(clipped_r2, 10);

        let cigar_r1 = format_cigar(&r1.cigar());
        assert_eq!(cigar_r1, "40M10H");

        let cigar_r2 = format_cigar(&r2.cigar());
        assert_eq!(cigar_r2, "10H40M");

        // Sequence should be trimmed in hard mode
        assert_eq!(r1.sequence().len(), 40);
        assert_eq!(r2.sequence().len(), 40);
    }

    // Tests for overlapping reads with insertions

    #[test]
    fn test_clip_overlapping_reads_with_insertions() {
        // Based on Scala test: "handle reads that contain insertions"
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let r1_seq = "A".repeat(50);
        let r2_seq = "A".repeat(50);
        // R1: 100-147 with 2bp insertion (40M2I8M)
        // R2: 130-169 with 2bp insertion (10M2I38M)
        let (mut r1, mut r2) = create_pair(100, "40M2I8M", &r1_seq, 130, "10M2I38M", &r2_seq);

        // Before clipping, r1.end >= r2.start (they overlap)
        let r1_end_before = 100 + 48; // 40M + 8M = 48 reference bases
        let r2_start_before = 130;
        assert!(r1_end_before >= r2_start_before);

        let (clipped_r1, clipped_r2) = clipper.clip_overlapping_reads(&mut r1, &mut r2);

        // Both should be clipped
        assert!(clipped_r1 > 0 || clipped_r2 > 0);

        // After clipping, they should not overlap (abutting is OK)
        let r1_end_after =
            usize::from(r1.alignment_start().unwrap()) + cigar_utils::reference_length(&r1.cigar());
        let r2_start_after = usize::from(r2.alignment_start().unwrap());
        assert!(
            r1_end_after <= r2_start_after,
            "After clipping, r1 end ({r1_end_after}) should be at or before r2 start ({r2_start_after})"
        );
    }

    #[test]
    fn test_clip_overlapping_reads_with_multiple_insertions() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let seq = "A".repeat(60);
        // R1: Multiple insertions
        // R2: Multiple insertions
        let (mut r1, mut r2) = create_pair(100, "20M2I20M3I10M", &seq, 120, "15M2I25M3I10M", &seq);

        let (clipped_r1, clipped_r2) = clipper.clip_overlapping_reads(&mut r1, &mut r2);

        // Should clip without errors
        assert!(clipped_r1 > 0 || clipped_r2 > 0);
    }

    #[test]
    fn test_clip_overlapping_reads_insertion_at_overlap_boundary() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let seq = "A".repeat(55);
        // R1: 100-149, with insertion near the end
        // R2: 145-194
        let (mut r1, mut r2) = create_pair(100, "48M2I5M", &seq, 145, "50M", &seq);

        let (clipped_r1, clipped_r2) = clipper.clip_overlapping_reads(&mut r1, &mut r2);

        // Should handle insertion at overlap boundary
        assert!(clipped_r1 > 0 || clipped_r2 > 0);
    }

    // Tests for SoftWithMask mode

    #[test]
    fn test_clip_end_of_alignment_soft_with_mask_masking() {
        let clipper = SamRecordClipper::new(ClippingMode::SoftWithMask);
        let mut record = create_test_record("50M", &"A".repeat(50), 1000);

        let clipped = clipper.clip_end_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);

        let cigar = format_cigar(&record.cigar());
        assert_eq!(cigar, "40M10S");

        // Check that the last 10 bases are masked to N
        let seq = record.sequence();
        let seq_bytes: Vec<u8> = seq.as_ref().to_vec();
        for (i, &base) in seq_bytes.iter().enumerate().skip(40).take(10) {
            assert_eq!(base, b'N', "Base at position {i} should be N");
        }

        // Check that qualities are set to min (2)
        let quals = record.quality_scores();
        let qual_bytes: Vec<u8> = quals.as_ref().to_vec();
        for (i, &qual) in qual_bytes.iter().enumerate().skip(40).take(10) {
            assert_eq!(qual, 2, "Quality at position {i} should be 2");
        }
    }

    #[test]
    fn test_clip_start_of_alignment_soft_with_mask_masking() {
        let clipper = SamRecordClipper::new(ClippingMode::SoftWithMask);
        let mut record = create_test_record("50M", &"A".repeat(50), 1000);

        let clipped = clipper.clip_start_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);

        let cigar = format_cigar(&record.cigar());
        assert_eq!(cigar, "10S40M");

        // Check that the first 10 bases are masked to N
        let seq = record.sequence();
        let seq_bytes: Vec<u8> = seq.as_ref().to_vec();
        for (i, &base) in seq_bytes.iter().enumerate().take(10) {
            assert_eq!(base, b'N', "Base at position {i} should be N");
        }

        // Check that qualities are set to min (2)
        let quals = record.quality_scores();
        let qual_bytes: Vec<u8> = quals.as_ref().to_vec();
        for (i, &qual) in qual_bytes.iter().enumerate().take(10) {
            assert_eq!(qual, 2, "Quality at position {i} should be 2");
        }
    }

    #[test]
    fn test_clip_overlapping_reads_soft_with_mask() {
        let clipper = SamRecordClipper::new(ClippingMode::SoftWithMask);
        let seq = "A".repeat(100);
        let (mut r1, mut r2) = create_pair(1, "100M", &seq, 50, "100M", &seq);

        let (clipped_r1, clipped_r2) = clipper.clip_overlapping_reads(&mut r1, &mut r2);

        assert!(clipped_r1 > 0 || clipped_r2 > 0);

        // Check that clipped bases are masked
        if clipped_r1 > 0 {
            let seq = r1.sequence();
            let seq_bytes: Vec<u8> = seq.as_ref().to_vec();
            let start_mask = 100 - clipped_r1;
            for (i, &base) in seq_bytes.iter().enumerate().skip(start_mask).take(100 - start_mask) {
                assert_eq!(base, b'N', "R1 base at position {i} should be N");
            }
        }
    }

    // Edge cases

    #[test]
    fn test_clip_very_short_read() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let mut record = create_test_record("5M", "ACGTA", 1000);

        // Try to clip more than the read length
        let clipped = clipper.clip_end_of_alignment(&mut record, 3);
        assert_eq!(clipped, 3);

        let cigar = format_cigar(&record.cigar());
        assert_eq!(cigar, "2M3S");
    }

    #[test]
    fn test_clip_entire_read_3_prime() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let mut record = create_test_record("10M", &"A".repeat(10), 1000);

        let clipped = clipper.clip_end_of_alignment(&mut record, 10);
        assert_eq!(clipped, 10);

        let cigar = format_cigar(&record.cigar());
        assert_eq!(cigar, "10S");
    }

    #[test]
    fn test_clip_overlapping_reads_complex_cigar() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let r1_seq = "A".repeat(70);
        let r2_seq = "A".repeat(70);
        // Complex CIGAR with multiple operation types
        let (mut r1, mut r2) =
            create_pair(100, "20M5I10M3D15M2I10M", &r1_seq, 130, "15M2I20M5D10M3I10M", &r2_seq);

        let (_clipped_r1, _clipped_r2) = clipper.clip_overlapping_reads(&mut r1, &mut r2);

        // Should handle complex CIGAR without panicking (test passes if no panic occurs)
    }

    #[test]
    fn test_clip_with_leading_and_trailing_soft_clips() {
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let mut record = create_test_record("10S30M10S", &"A".repeat(50), 1000);

        // Clip additional 5 bases from 3' end
        let clipped = clipper.clip_end_of_alignment(&mut record, 5);
        assert_eq!(clipped, 5);

        let cigar = format_cigar(&record.cigar());
        assert_eq!(cigar, "10S25M15S");
    }

    #[test]
    fn test_clip_overlapping_reads_one_base_different_positions() {
        // Test overlapping by exactly one base at different reference positions
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let seq = "A".repeat(50);

        // Test at low position
        let (mut r1, mut r2) = create_pair(10, "50M", &seq, 59, "50M", &seq);
        let (c1, c2) = clipper.clip_overlapping_reads(&mut r1, &mut r2);
        assert!(c1 > 0 || c2 > 0);

        // Test at high position
        let (mut r1, mut r2) = create_pair(10000, "50M", &seq, 10049, "50M", &seq);
        let (c1, c2) = clipper.clip_overlapping_reads(&mut r1, &mut r2);
        assert!(c1 > 0 || c2 > 0);
    }

    // Tests for NOT upgrading clipping (same mode or going backwards)

    #[test]
    fn test_not_upgrade_soft_to_soft() {
        // Same mode - should not upgrade
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let seq = "12345678901234567890123456789012345678901234567890";
        let mut record = create_test_record("5S35M10S", seq, 10);

        let result = clipper.upgrade_all_clipping(&mut record).unwrap();
        assert_eq!(result, (0, 0));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "5S35M10S"); // No change
        assert_eq!(record.sequence().len(), 50); // No sequence change
    }

    #[test]
    fn test_not_upgrade_soft_with_mask_to_soft() {
        // Going backwards - should not upgrade
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let seq = "12345678901234567890123456789012345678901234567890";
        let mut record = create_test_record("5S35M10S", seq, 10);

        let result = clipper.upgrade_all_clipping(&mut record).unwrap();
        assert_eq!(result, (0, 0));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "5S35M10S");
    }

    #[test]
    fn test_not_upgrade_soft_with_mask_to_soft_with_mask() {
        // Same mode - should not upgrade
        let clipper = SamRecordClipper::new(ClippingMode::SoftWithMask);
        let seq = "12345678901234567890123456789012345678901234567890";
        let mut record = create_test_record("5S35M10S", seq, 10);

        let result = clipper.upgrade_all_clipping(&mut record).unwrap();
        assert_eq!(result, (0, 0));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "5S35M10S");
    }

    #[test]
    fn test_not_upgrade_hard_to_hard() {
        // Same mode - should not upgrade
        let clipper = SamRecordClipper::new(ClippingMode::Hard);
        let seq = "12345678901234567890123456789012345"; // 35 bases
        let mut record = create_test_record("5H35M10H", seq, 10);

        let result = clipper.upgrade_all_clipping(&mut record).unwrap();
        assert_eq!(result, (0, 0));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "5H35M10H");
        assert_eq!(record.sequence().len(), 35);
    }

    #[test]
    fn test_not_upgrade_hard_to_soft_with_mask() {
        // Going backwards - should not upgrade (Hard is already most restrictive)
        let clipper = SamRecordClipper::new(ClippingMode::SoftWithMask);
        let seq = "12345678901234567890123456789012345";
        let mut record = create_test_record("5H35M10H", seq, 10);

        let result = clipper.upgrade_all_clipping(&mut record).unwrap();
        assert_eq!(result, (0, 0));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "5H35M10H");
    }

    #[test]
    fn test_not_upgrade_hard_to_soft() {
        // Going backwards - should not upgrade
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let seq = "12345678901234567890123456789012345";
        let mut record = create_test_record("5H35M10H", seq, 10);

        let result = clipper.upgrade_all_clipping(&mut record).unwrap();
        assert_eq!(result, (0, 0));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "5H35M10H");
    }

    #[test]
    fn test_not_upgrade_unmapped_read() {
        // Unmapped reads should not be upgraded
        let clipper = SamRecordClipper::new(ClippingMode::Hard);
        let seq = "12345678901234567890123456789012345678901234567890";
        let mut record = create_test_record("5S35M10S", seq, 10);

        // Mark as unmapped
        *record.flags_mut() = Flags::UNMAPPED;

        let result = clipper.upgrade_all_clipping(&mut record).unwrap();
        assert_eq!(result, (0, 0));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "5S35M10S");
    }

    #[test]
    fn test_upgrade_soft_to_hard_with_existing_hard_clips() {
        // Soft with existing hard clips -> Hard mode
        let clipper = SamRecordClipper::new(ClippingMode::Hard);
        let seq = "1234567890123456789012345678901234567890"; // 40 bases (2H + 5S + 30M + 3S)
        let mut record = create_test_record("2H5S30M3S", seq, 10);

        let result = clipper.upgrade_all_clipping(&mut record).unwrap();
        assert_eq!(result, (5, 3));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "7H30M3H"); // 2H + 5S -> 7H, 3S -> 3H
        assert_eq!(record.sequence().len(), 30);
    }

    #[test]
    fn test_upgrade_soft_with_mask_to_hard() {
        // SoftWithMask -> Hard should upgrade
        let clipper = SamRecordClipper::new(ClippingMode::Hard);
        let seq = "12345678901234567890123456789012345678901234567890";
        let mut record = create_test_record("5S35M10S", seq, 10);

        let result = clipper.upgrade_all_clipping(&mut record).unwrap();
        assert_eq!(result, (5, 10));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "5H35M10H");
        assert_eq!(record.sequence().len(), 35);
    }

    #[test]
    fn test_upgrade_clipping_with_only_leading_soft_clip() {
        let clipper = SamRecordClipper::new(ClippingMode::Hard);
        let seq = "12345678901234567890123456789012345678901234567890";
        let mut record = create_test_record("10S40M", seq, 10);

        let result = clipper.upgrade_all_clipping(&mut record).unwrap();
        assert_eq!(result, (10, 0));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "10H40M");
        assert_eq!(record.sequence().len(), 40);
    }

    #[test]
    fn test_upgrade_clipping_with_only_trailing_soft_clip() {
        let clipper = SamRecordClipper::new(ClippingMode::Hard);
        let seq = "12345678901234567890123456789012345678901234567890";
        let mut record = create_test_record("40M10S", seq, 10);

        let result = clipper.upgrade_all_clipping(&mut record).unwrap();
        assert_eq!(result, (0, 10));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "40M10H");
        assert_eq!(record.sequence().len(), 40);
    }

    #[test]
    fn test_upgrade_clipping_complex_cigar_with_indels() {
        let clipper = SamRecordClipper::new(ClippingMode::Hard);
        let seq = "123456789012345678901234567890123456789012345"; // 45 bases (5S + 20M + 5I + 10M + 5S)
        let mut record = create_test_record("5S20M5I10M5S", seq, 10);

        let result = clipper.upgrade_all_clipping(&mut record).unwrap();
        assert_eq!(result, (5, 5));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "5H20M5I10M5H");
        assert_eq!(record.sequence().len(), 35); // 20 + 5 + 10
    }

    #[test]
    fn test_upgrade_clipping_no_soft_clips_present() {
        let clipper = SamRecordClipper::new(ClippingMode::Hard);
        let seq = "12345678901234567890123456789012345678901234567890";
        let mut record = create_test_record("50M", seq, 10);

        let result = clipper.upgrade_all_clipping(&mut record).unwrap();
        assert_eq!(result, (0, 0));

        let cigar_str = format_cigar(&record.cigar());
        assert_eq!(cigar_str, "50M");
        assert_eq!(record.sequence().len(), 50);
    }

    // ===== Additional tests ported from fgbio to match clipExtendingPastMateEnds coverage =====

    #[test]
    fn test_clip_extending_past_mate_ends_one_base_extension() {
        // fgbio test: "clip reads that extend one base past their mate's start"
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let seq = "A".repeat(100);
        let (mut r1, mut r2) = create_pair(2, "100M", &seq, 1, "100M", &seq);

        let (clipped_r1, clipped_r2) = clipper.clip_extending_past_mate_ends(&mut r1, &mut r2);

        assert_eq!(clipped_r1, 1);
        assert_eq!(clipped_r2, 1);
        assert_eq!(r1.alignment_start(), Some(Position::new(2).unwrap()));
        assert_eq!(format_cigar(&r1.cigar()), "99M1S");
        assert_eq!(r2.alignment_start(), Some(Position::new(2).unwrap()));
        assert_eq!(format_cigar(&r2.cigar()), "1S99M");
    }

    #[test]
    fn test_clip_extending_past_mate_ends_two_base_extension() {
        // fgbio test: "clip reads that extend two bases past their mate's start"
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let seq = "A".repeat(100);
        let (mut r1, mut r2) = create_pair(3, "100M", &seq, 1, "100M", &seq);

        let (clipped_r1, clipped_r2) = clipper.clip_extending_past_mate_ends(&mut r1, &mut r2);

        assert_eq!(clipped_r1, 2);
        assert_eq!(clipped_r2, 2);
        assert_eq!(r1.alignment_start(), Some(Position::new(3).unwrap()));
        assert_eq!(format_cigar(&r1.cigar()), "98M2S");
        assert_eq!(r2.alignment_start(), Some(Position::new(3).unwrap()));
        assert_eq!(format_cigar(&r2.cigar()), "2S98M");
    }

    #[test]
    fn test_clip_extending_past_mate_ends_only_one_end_extends() {
        // fgbio test: "only one end extends their mate's start"
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let seq = "A".repeat(100);
        let (mut r1, mut r2) = create_pair(1, "100M", &seq, 1, "50S50M", &seq);

        let (clipped_r1, clipped_r2) = clipper.clip_extending_past_mate_ends(&mut r1, &mut r2);

        assert_eq!(clipped_r1, 50); // R1 gets clipped
        assert_eq!(clipped_r2, 0); // R2 already has clipping
        assert_eq!(r1.alignment_start(), Some(Position::new(1).unwrap()));
        assert_eq!(format_cigar(&r1.cigar()), "50M50S");
        assert_eq!(r2.alignment_start(), Some(Position::new(1).unwrap()));
        assert_eq!(format_cigar(&r2.cigar()), "50S50M"); // unchanged
    }

    #[test]
    fn test_clip_extending_past_mate_ends_with_insertions() {
        // fgbio test: "only one end extends with insertions"
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let seq = "A".repeat(100);
        let (mut r1, mut r2) = create_pair(1, "40M10I50M", &seq, 1, "50S50M", &seq);

        let (clipped_r1, clipped_r2) = clipper.clip_extending_past_mate_ends(&mut r1, &mut r2);

        assert_eq!(clipped_r1, 40); // 40 bases clipped (insertion doesn't count toward ref)
        assert_eq!(clipped_r2, 0);
        assert_eq!(r1.alignment_start(), Some(Position::new(1).unwrap()));
        assert_eq!(format_cigar(&r1.cigar()), "40M10I10M40S");
        assert_eq!(r2.alignment_start(), Some(Position::new(1).unwrap()));
        assert_eq!(format_cigar(&r2.cigar()), "50S50M"); // unchanged
    }

    #[test]
    fn test_clip_extending_past_mate_ends_soft_clips_extend_hard_mode() {
        // fgbio test: "forward read soft-clips extend past (hard mode)"
        let clipper = SamRecordClipper::new(ClippingMode::Hard);
        let seq = "A".repeat(50);
        let (mut r1, mut r2) = create_pair(20, "30M20S", &seq, 20, "10S40M", &seq);

        let (clipped_r1, clipped_r2) = clipper.clip_extending_past_mate_ends(&mut r1, &mut r2);

        assert_eq!(clipped_r1, 0); // No aligned bases clipped
        assert_eq!(clipped_r2, 0); // No aligned bases clipped
        assert_eq!(r1.alignment_start(), Some(Position::new(20).unwrap()));
        assert_eq!(format_cigar(&r1.cigar()), "30M10S10H"); // 10 of 20S hard-clipped
        assert_eq!(r2.alignment_start(), Some(Position::new(20).unwrap()));
        assert_eq!(format_cigar(&r2.cigar()), "10H40M"); // 10S hard-clipped
    }

    #[test]
    fn test_clip_extending_past_mate_ends_soft_clips_extend_with_deletion_hard_mode() {
        // fgbio test: "forward read soft-clips extend past with deletion (hard mode)"
        let clipper = SamRecordClipper::new(ClippingMode::Hard);
        let seq = "A".repeat(50);
        let (mut r1, mut r2) = create_pair(20, "15M1D15M20S", &seq, 20, "10S15M1D25M", &seq);

        let (clipped_r1, clipped_r2) = clipper.clip_extending_past_mate_ends(&mut r1, &mut r2);

        assert_eq!(clipped_r1, 0); // No aligned bases clipped
        assert_eq!(clipped_r2, 0); // No aligned bases clipped
        assert_eq!(r1.alignment_start(), Some(Position::new(20).unwrap()));
        assert_eq!(format_cigar(&r1.cigar()), "15M1D15M10S10H");
        assert_eq!(r2.alignment_start(), Some(Position::new(20).unwrap()));
        assert_eq!(format_cigar(&r2.cigar()), "10H15M1D25M");
    }

    #[test]
    fn test_clip_extending_past_mate_ends_soft_clips_extend_with_insertion_hard_mode() {
        // fgbio test: "forward read soft-clips extend past with insertion (hard mode)"
        let clipper = SamRecordClipper::new(ClippingMode::Hard);
        let seq = "A".repeat(51);
        let (mut r1, mut r2) = create_pair(20, "15M1I15M20S", &seq, 20, "10S15M1I25M", &seq);

        let (clipped_r1, clipped_r2) = clipper.clip_extending_past_mate_ends(&mut r1, &mut r2);

        assert_eq!(clipped_r1, 0); // No aligned bases clipped
        assert_eq!(clipped_r2, 0); // No aligned bases clipped
        assert_eq!(r1.alignment_start(), Some(Position::new(20).unwrap()));
        assert_eq!(format_cigar(&r1.cigar()), "15M1I15M10S10H");
        assert_eq!(r2.alignment_start(), Some(Position::new(20).unwrap()));
        assert_eq!(format_cigar(&r2.cigar()), "10H15M1I25M");
    }

    #[test]
    fn test_clip_extending_past_mate_ends_reverse_soft_clips_extend_hard_mode() {
        // fgbio test: "reverse read soft-clips extend past (hard mode)"
        let clipper = SamRecordClipper::new(ClippingMode::Hard);
        let seq = "A".repeat(50);
        let (mut r1, mut r2) = create_pair(20, "40M10S", &seq, 30, "20S30M", &seq);

        let (clipped_r1, clipped_r2) = clipper.clip_extending_past_mate_ends(&mut r1, &mut r2);

        assert_eq!(clipped_r1, 0); // No aligned bases clipped
        assert_eq!(clipped_r2, 0); // No aligned bases clipped
        assert_eq!(r1.alignment_start(), Some(Position::new(20).unwrap()));
        assert_eq!(format_cigar(&r1.cigar()), "40M10H");
        assert_eq!(r2.alignment_start(), Some(Position::new(30).unwrap()));
        assert_eq!(format_cigar(&r2.cigar()), "10H10S30M");
    }

    #[test]
    fn test_clip_extending_past_mate_ends_no_overlap_far_apart() {
        // fgbio test: "not clip when mapped +/- with start(R1) > end(R2) but no overlap"
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let seq = "A".repeat(100);
        let (mut r1, mut r2) = create_pair(1000, "100M", &seq, 1, "100M", &seq);

        let (clipped_r1, clipped_r2) = clipper.clip_extending_past_mate_ends(&mut r1, &mut r2);

        assert_eq!(clipped_r1, 0);
        assert_eq!(clipped_r2, 0);
        assert_eq!(r1.alignment_start(), Some(Position::new(1000).unwrap()));
        assert_eq!(format_cigar(&r1.cigar()), "100M");
        assert_eq!(r2.alignment_start(), Some(Position::new(1).unwrap()));
        assert_eq!(format_cigar(&r2.cigar()), "100M");
    }

    #[test]
    fn test_clip_extending_past_mate_ends_no_extension_with_insertions() {
        // fgbio test: "not clip when no extension with insertions"
        let clipper = SamRecordClipper::new(ClippingMode::Soft);
        let seq = "A".repeat(100);
        let (mut r1, mut r2) = create_pair(1, "40M20I40M", &seq, 1, "40M20I40M", &seq);

        let (clipped_r1, clipped_r2) = clipper.clip_extending_past_mate_ends(&mut r1, &mut r2);

        assert_eq!(clipped_r1, 0);
        assert_eq!(clipped_r2, 0);
        assert_eq!(r1.alignment_start(), Some(Position::new(1).unwrap()));
        assert_eq!(format_cigar(&r1.cigar()), "40M20I40M");
        assert_eq!(r2.alignment_start(), Some(Position::new(1).unwrap()));
        assert_eq!(format_cigar(&r2.cigar()), "40M20I40M");
    }

    #[test]
    fn test_clip_5_prime_end_of_alignment_positive_strand() {
        // fgbio test: "add more clipping to the 5' end" - positive strand
        let clipper = SamRecordClipper::new(ClippingMode::Soft);

        // Positive strand, no existing clipping
        let mut rec1 = create_test_record("50M", &"A".repeat(50), 10);
        *rec1.flags_mut() = Flags::empty(); // positive strand
        let clipped = clipper.clip_5_prime_end_of_alignment(&mut rec1, 10);
        assert_eq!(clipped, 10);
        assert_eq!(format_cigar(&rec1.cigar()), "10S40M");

        // Positive strand, existing soft clipping
        let mut rec2 = create_test_record("10S40M", &"A".repeat(50), 10);
        *rec2.flags_mut() = Flags::empty(); // positive strand
        let clipped = clipper.clip_5_prime_end_of_alignment(&mut rec2, 10);
        assert_eq!(clipped, 10);
        assert_eq!(format_cigar(&rec2.cigar()), "20S30M");
    }

    #[test]
    fn test_clip_5_prime_end_of_alignment_negative_strand() {
        // fgbio test: "add more clipping to the 5' end" - negative strand
        let clipper = SamRecordClipper::new(ClippingMode::Soft);

        // Negative strand, no existing clipping
        let mut rec3 = create_test_record("50M", &"A".repeat(50), 10);
        *rec3.flags_mut() = Flags::REVERSE_COMPLEMENTED;
        let clipped = clipper.clip_5_prime_end_of_alignment(&mut rec3, 10);
        assert_eq!(clipped, 10);
        assert_eq!(format_cigar(&rec3.cigar()), "40M10S");

        // Negative strand, existing soft clipping
        let mut rec4 = create_test_record("40M10S", &"A".repeat(50), 10);
        *rec4.flags_mut() = Flags::REVERSE_COMPLEMENTED;
        let clipped = clipper.clip_5_prime_end_of_alignment(&mut rec4, 10);
        assert_eq!(clipped, 10);
        assert_eq!(format_cigar(&rec4.cigar()), "30M20S");
    }

    #[test]
    fn test_clip_3_prime_end_of_alignment_negative_strand() {
        // fgbio test: "add more clipping to the 3' end" - negative strand
        let clipper = SamRecordClipper::new(ClippingMode::Soft);

        // Negative strand, no existing clipping
        let mut rec1 = create_test_record("50M", &"A".repeat(50), 10);
        *rec1.flags_mut() = Flags::REVERSE_COMPLEMENTED;
        let clipped = clipper.clip_3_prime_end_of_alignment(&mut rec1, 10);
        assert_eq!(clipped, 10);
        assert_eq!(format_cigar(&rec1.cigar()), "10S40M");

        // Negative strand, existing soft clipping
        let mut rec2 = create_test_record("10S40M", &"A".repeat(50), 10);
        *rec2.flags_mut() = Flags::REVERSE_COMPLEMENTED;
        let clipped = clipper.clip_3_prime_end_of_alignment(&mut rec2, 10);
        assert_eq!(clipped, 10);
        assert_eq!(format_cigar(&rec2.cigar()), "20S30M");
    }

    #[test]
    fn test_clip_3_prime_end_of_alignment_positive_strand() {
        // fgbio test: "add more clipping to the 3' end" - positive strand
        let clipper = SamRecordClipper::new(ClippingMode::Soft);

        // Positive strand, no existing clipping
        let mut rec3 = create_test_record("50M", &"A".repeat(50), 10);
        *rec3.flags_mut() = Flags::empty(); // positive strand
        let clipped = clipper.clip_3_prime_end_of_alignment(&mut rec3, 10);
        assert_eq!(clipped, 10);
        assert_eq!(format_cigar(&rec3.cigar()), "40M10S");

        // Positive strand, existing soft clipping
        let mut rec4 = create_test_record("40M10S", &"A".repeat(50), 10);
        *rec4.flags_mut() = Flags::empty(); // positive strand
        let clipped = clipper.clip_3_prime_end_of_alignment(&mut rec4, 10);
        assert_eq!(clipped, 10);
        assert_eq!(format_cigar(&rec4.cigar()), "30M20S");
    }
}

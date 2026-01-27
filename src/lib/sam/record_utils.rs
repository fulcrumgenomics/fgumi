//! Record-level utilities for SAM/BAM records.
//!
//! This module provides utilities for working with individual SAM records, including:
//! - Position mapping between read and reference coordinates
//! - FR pair detection using BAM tags
//! - Mate position calculations from MC tag
//! - CIGAR string parsing and analysis

use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::record::Cigar as CigarTrait;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value;

/// Pair orientation for paired-end reads.
///
/// This enum represents the orientation of read pairs based on
/// the relative positioning and strand of the two reads in a pair.
///
/// Note: This is NOT the orientation byte as used in Picard's `MarkDuplicates`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PairOrientation {
    /// Forward-Reverse orientation ("innie") - reads face each other:
    /// ```text
    /// 5' --F-->       <--R-- 5'
    /// ```
    /// The positive strand 5' position is less than the negative strand 5' position.
    FR,

    /// Reverse-Forward orientation ("outie") - reads face away from each other:
    /// ```text
    /// <--R-- 5'       5' --F-->
    /// ```
    /// The positive strand 5' position is greater than or equal to the negative strand 5' position.
    RF,

    /// Tandem orientation - both reads on the same strand:
    /// ```text
    /// 5' --F-->  5' --F-->   (both forward)
    /// <--R-- 5'  <--R-- 5'   (both reverse)
    /// ```
    Tandem,
}

/// Returns the 1-based read position corresponding to a reference position.
///
/// This matches fgbio's `readPosAtRefPos` function (via htsjdk's `getReadPositionAtReferencePosition`).
///
/// # Arguments
/// * `record` - The SAM record
/// * `ref_pos` - The 1-based reference position to query
/// * `return_last_base_if_deleted` - If true, returns the last aligned base position when
///   the `ref_pos` falls in a deletion; if false, returns 0
///
/// # Returns
/// * 1-based read position, or 0 if:
///   - The position is not found within the alignment
///   - The position falls in a deletion (when `return_last_base_if_deleted` is false)
///   - The record has no alignment start
#[must_use]
pub fn read_pos_at_ref_pos(
    record: &RecordBuf,
    ref_pos: usize,
    return_last_base_if_deleted: bool,
) -> usize {
    let Some(alignment_start) = record.alignment_start().map(usize::from) else {
        return 0;
    };

    // Walk the CIGAR tracking both read and reference positions
    let mut read_pos: usize = 0; // 0-based position in read
    let mut ref_cursor = alignment_start; // 1-based reference position
    let mut last_aligned_read_pos: usize = 0;

    for op_result in record.cigar().iter() {
        let Ok(op) = op_result else {
            continue;
        };
        let len = op.len();

        match op.kind() {
            // Consumes both read and reference
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                // Check if ref_pos falls within this aligned segment
                if ref_pos >= ref_cursor && ref_pos < ref_cursor + len {
                    // Found it! Return the 1-based read position
                    let offset = ref_pos - ref_cursor;
                    return read_pos + offset + 1; // +1 for 1-based
                }
                last_aligned_read_pos = read_pos + len; // Track last aligned position
                read_pos += len;
                ref_cursor += len;
            }
            // Consumes only read (insertion, soft clip)
            Kind::Insertion | Kind::SoftClip => {
                read_pos += len;
            }
            // Consumes only reference (deletion, skip)
            Kind::Deletion | Kind::Skip => {
                // Check if ref_pos falls within this deletion/skip
                if ref_pos >= ref_cursor && ref_pos < ref_cursor + len {
                    return if return_last_base_if_deleted && last_aligned_read_pos > 0 {
                        last_aligned_read_pos // 1-based (already 1-based because of how we track)
                    } else {
                        0
                    };
                }
                ref_cursor += len;
            }
            // Hard clip doesn't consume read bases (they're gone)
            Kind::HardClip | Kind::Pad => {}
        }
    }

    // Position not found within alignment
    0
}

/// Checks if a single read is part of an FR (forward-reverse) pair using BAM tags.
///
/// This function determines FR pair status from a single read using information
/// available in the BAM record (mate position, template length, strand flags).
/// This is useful when you don't have access to the mate record.
///
/// Returns `true` if:
/// - Read is paired
/// - Both read and mate are mapped
/// - Both are on the same reference
/// - The pair is in FR orientation (positive strand 5' < negative strand 5')
///
/// This matches fgbio's `isFrPair` algorithm using htsjdk's pair orientation logic.
#[must_use]
pub fn is_fr_pair_from_tags(read: &RecordBuf) -> bool {
    let flags = read.flags();

    // Must be paired
    if !flags.is_segmented() {
        return false;
    }

    // Both read and mate must be mapped
    if flags.is_unmapped() || flags.is_mate_unmapped() {
        return false;
    }

    // Must be on the same reference as mate
    let ref_id = read.reference_sequence_id();
    let mate_ref_id = read.mate_reference_sequence_id();
    if ref_id != mate_ref_id {
        return false;
    }

    // Check strand orientation - must be opposite strands for FR or RF
    let is_reverse = flags.is_reverse_complemented();
    let mate_reverse = flags.is_mate_reverse_complemented();
    if is_reverse == mate_reverse {
        // Same strand = TANDEM, not FR
        return false;
    }

    // Now determine if FR or RF using htsjdk's logic:
    // positiveStrandFivePrimePos = readIsOnReverseStrand ? mateStart : alignmentStart
    // negativeStrandFivePrimePos = readIsOnReverseStrand ? alignmentEnd : alignmentStart + insertSize
    let alignment_start = read.alignment_start().map_or(0, usize::from);
    let mate_start = read.mate_alignment_start().map_or(0, usize::from);
    let insert_size = read.template_length();

    let (positive_five_prime, negative_five_prime) = if is_reverse {
        // This read is on reverse strand, mate is on positive strand
        // Need alignment end for this read's 5' position
        let ref_len = reference_length(&read.cigar());
        let end = alignment_start + ref_len.saturating_sub(1);
        // positiveStrandFivePrimePos = mateStart
        // negativeStrandFivePrimePos = alignmentEnd (this read)
        (mate_start as i32, end as i32)
    } else {
        // This read is on positive strand, mate is on reverse strand
        // positiveStrandFivePrimePos = alignmentStart (this read)
        // negativeStrandFivePrimePos = alignmentStart + insertSize
        (alignment_start as i32, alignment_start as i32 + insert_size)
    };

    // FR if positive strand 5' < negative strand 5'
    positive_five_prime < negative_five_prime
}

/// Gets the mate's unclipped start position from the MC tag and mate position.
///
/// Calculates: `mateUnclippedStart = mateStart - leadingClipping`
///
/// Returns `None` if the MC tag is missing or invalid, or if mate position is missing.
#[must_use]
/// Returns the mate's unclipped start as a signed integer.
///
/// This uses signed arithmetic to properly handle cases where the unclipped start
/// would be before position 1 (i.e., negative in 1-based coordinates).
/// This matches fgbio's behavior which uses signed integers for these calculations.
pub fn mate_unclipped_start(read: &RecordBuf) -> Option<isize> {
    let mc_tag = Tag::from([b'M', b'C']);
    let mc_value = read.data().get(&mc_tag)?;

    let cigar_str = match mc_value {
        Value::String(s) => String::from_utf8_lossy(s.as_ref()).to_string(),
        _ => return None,
    };

    let mate_start = usize::from(read.mate_alignment_start()?) as isize;
    let ops = parse_cigar_string(&cigar_str);
    let leading_clip = leading_clipping(&ops) as isize;

    Some(mate_start - leading_clip)
}

/// Gets the mate's unclipped end position from the MC tag and mate position.
///
/// Calculates: `mateUnclippedEnd = mateStart + refLen - 1 + trailingClipping`
///
/// Returns `None` if the MC tag is missing or invalid, or if mate position is missing.
#[must_use]
pub fn mate_unclipped_end(read: &RecordBuf) -> Option<usize> {
    let mc_tag = Tag::from([b'M', b'C']);
    let mc_value = read.data().get(&mc_tag)?;

    let cigar_str = match mc_value {
        Value::String(s) => String::from_utf8_lossy(s.as_ref()).to_string(),
        _ => return None,
    };

    let mate_start = usize::from(read.mate_alignment_start()?);
    let ops = parse_cigar_string(&cigar_str);
    let ref_len = cigar_reference_length(&ops);
    let trailing_clip = trailing_clipping(&ops);

    // mateStart + refLen - 1 + trailingClipping
    Some(mate_start + ref_len - 1 + trailing_clip)
}

/// Calculates the number of bases in a read that extend past the mate's unclipped boundary.
///
/// This matches fgbio's `numBasesExtendingPastMate` function exactly.
/// Uses tag-based mate information (MC tag, mate position, template length).
///
/// Returns 0 if not an FR pair or if required tags are missing.
#[must_use]
pub fn num_bases_extending_past_mate(read: &RecordBuf) -> usize {
    // Only applies to FR pairs
    if !is_fr_pair_from_tags(read) {
        return 0;
    }

    let is_positive_strand = !read.flags().is_reverse_complemented();

    // Get the read's CIGAR ops for soft clip calculations
    let ops: Vec<(Kind, usize)> = read
        .cigar()
        .iter()
        .filter_map(std::result::Result::ok)
        .map(|op| (op.kind(), op.len()))
        .collect();

    // Get read length (total bases in the read, excluding hard clips)
    let read_length: usize = ops
        .iter()
        .filter(|(kind, _)| {
            matches!(
                kind,
                Kind::Match
                    | Kind::SequenceMatch
                    | Kind::SequenceMismatch
                    | Kind::Insertion
                    | Kind::SoftClip
            )
        })
        .map(|(_, len)| len)
        .sum();

    if is_positive_strand {
        // Positive strand: check if read extends past mate's unclipped end
        let Some(read_alignment_end) = alignment_end(read) else {
            return 0;
        };
        let mate_unclipped_end_pos = mate_unclipped_end(read).unwrap_or(usize::MAX);

        if read_alignment_end >= mate_unclipped_end_pos {
            // Aligned portion reaches or extends past mate's unclipped end
            // Use readPosAtRefPos to find where in the read the mate ends
            // fgbio: Math.max(0, rec.length - rec.readPosAtRefPos(pos=mateUnclippedEnd, returnLastBaseIfDeleted=false))
            let pos_at_mate_end = read_pos_at_ref_pos(read, mate_unclipped_end_pos, false);
            // When pos_at_mate_end is 0 (position in deletion or outside alignment),
            // fgbio's formula gives: max(0, read_length - 0) = read_length (clip all)
            read_length.saturating_sub(pos_at_mate_end)
        } else {
            // Aligned portion ends before mate's unclipped end
            // Only clip excess soft-clipped bases that extend past the mate
            let trailing_soft_clip = trailing_soft_clipping(&ops);
            let gap = mate_unclipped_end_pos - read_alignment_end;
            trailing_soft_clip.saturating_sub(gap)
        }
    } else {
        // Negative strand: check if read extends before mate's unclipped start
        let Some(read_alignment_start) = read.alignment_start().map(usize::from) else {
            return 0;
        };
        let read_alignment_start_signed = read_alignment_start as isize;
        let mate_unclipped_start_pos = mate_unclipped_start(read).unwrap_or(0);

        if read_alignment_start_signed <= mate_unclipped_start_pos {
            // Aligned portion starts at or before mate's unclipped start
            // Use readPosAtRefPos to find where in the read the mate starts
            // fgbio: Math.max(0, rec.readPosAtRefPos(pos=mateUnclippedStart, returnLastBaseIfDeleted=false) - 1)
            // Note: mate_unclipped_start_pos could be negative, but read_pos_at_ref_pos expects usize
            // If mate_unclipped_start_pos is negative, treat it as 0 for the lookup
            let lookup_pos =
                if mate_unclipped_start_pos < 0 { 0 } else { mate_unclipped_start_pos as usize };
            let pos_at_mate_start = read_pos_at_ref_pos(read, lookup_pos, false);
            // When pos_at_mate_start is 0 (position in deletion or outside alignment),
            // fgbio's formula gives: max(0, 0 - 1) = 0 (clip nothing)
            // pos_at_mate_start is 1-based, so subtract 1 to get bases before it
            pos_at_mate_start.saturating_sub(1)
        } else {
            // Aligned portion starts after mate's unclipped start
            // Only clip excess soft-clipped bases that extend before the mate
            // Use signed arithmetic to properly handle negative mate_unclipped_start
            let leading_soft_clip = leading_soft_clipping(&ops) as isize;
            let gap = read_alignment_start_signed - mate_unclipped_start_pos;
            // gap is always positive here (since read_alignment_start > mate_unclipped_start)
            // clip = leading_soft_clip - gap, clamped to 0
            if leading_soft_clip > gap { (leading_soft_clip - gap) as usize } else { 0 }
        }
    }
}

/// Gets the read's alignment end position (1-based, inclusive).
///
/// Calculated as: `alignment_start + reference_length - 1`
#[must_use]
pub fn alignment_end(read: &RecordBuf) -> Option<usize> {
    let start = usize::from(read.alignment_start()?);
    let ref_len = reference_length(&read.cigar());
    Some(start + ref_len - 1)
}

/// Parses a CIGAR string and returns it as a vector of (Kind, length) operations.
///
/// This is useful for parsing the MC (mate CIGAR) tag value.
#[must_use]
pub fn parse_cigar_string(cigar_str: &str) -> Vec<(Kind, usize)> {
    let mut ops = Vec::new();
    let mut num_str = String::new();

    for ch in cigar_str.chars() {
        if ch.is_ascii_digit() {
            num_str.push(ch);
        } else {
            let len: usize = num_str.parse().unwrap_or(0);
            num_str.clear();

            let kind = match ch {
                'M' => Kind::Match,
                'I' => Kind::Insertion,
                'D' => Kind::Deletion,
                'N' => Kind::Skip,
                'S' => Kind::SoftClip,
                'H' => Kind::HardClip,
                'P' => Kind::Pad,
                '=' => Kind::SequenceMatch,
                'X' => Kind::SequenceMismatch,
                _ => continue,
            };

            if len > 0 {
                ops.push((kind, len));
            }
        }
    }

    ops
}

/// Calculates the reference length consumed by CIGAR operations.
#[must_use]
pub fn cigar_reference_length(ops: &[(Kind, usize)]) -> usize {
    ops.iter()
        .filter_map(|(kind, len)| match kind {
            Kind::Match
            | Kind::Deletion
            | Kind::Skip
            | Kind::SequenceMatch
            | Kind::SequenceMismatch => Some(*len),
            _ => None,
        })
        .sum()
}

/// Calculates leading clipping (soft + hard) from CIGAR operations.
#[must_use]
pub fn leading_clipping(ops: &[(Kind, usize)]) -> usize {
    ops.iter()
        .take_while(|(kind, _)| matches!(kind, Kind::SoftClip | Kind::HardClip))
        .map(|(_, len)| *len)
        .sum()
}

/// Calculates trailing clipping (soft + hard) from CIGAR operations.
#[must_use]
pub fn trailing_clipping(ops: &[(Kind, usize)]) -> usize {
    ops.iter()
        .rev()
        .take_while(|(kind, _)| matches!(kind, Kind::SoftClip | Kind::HardClip))
        .map(|(_, len)| *len)
        .sum()
}

/// Calculates leading soft clipping only from CIGAR operations.
#[must_use]
pub fn leading_soft_clipping(ops: &[(Kind, usize)]) -> usize {
    ops.iter()
        .skip_while(|(kind, _)| *kind == Kind::HardClip)
        .take_while(|(kind, _)| *kind == Kind::SoftClip)
        .map(|(_, len)| *len)
        .sum()
}

/// Calculates trailing soft clipping only from CIGAR operations.
#[must_use]
pub fn trailing_soft_clipping(ops: &[(Kind, usize)]) -> usize {
    ops.iter()
        .rev()
        .skip_while(|(kind, _)| *kind == Kind::HardClip)
        .take_while(|(kind, _)| *kind == Kind::SoftClip)
        .map(|(_, len)| *len)
        .sum()
}

/// Collects CIGAR operations from a record into a Vec for use with clipping functions.
#[must_use]
fn cigar_to_ops(record: &RecordBuf) -> Vec<(Kind, usize)> {
    record.cigar().as_ref().iter().map(|op| (op.kind(), op.len())).collect()
}

/// Gets the unclipped start position of a read (alignment start minus leading clips).
///
/// This matches HTSJDK's `SAMRecord.getUnclippedStart()` behavior, which includes
/// both soft clips and hard clips.
///
/// Returns `None` for unmapped reads.
#[must_use]
pub fn unclipped_start(record: &RecordBuf) -> Option<usize> {
    if record.flags().is_unmapped() {
        return None;
    }
    let start = usize::from(record.alignment_start()?);
    let leading = leading_clipping(&cigar_to_ops(record));
    Some(start.saturating_sub(leading))
}

/// Gets the unclipped end position of a read (alignment end plus trailing clips).
///
/// This matches HTSJDK's `SAMRecord.getUnclippedEnd()` behavior, which includes
/// both soft clips and hard clips.
///
/// Returns `None` for unmapped reads.
#[must_use]
pub fn unclipped_end(record: &RecordBuf) -> Option<usize> {
    if record.flags().is_unmapped() {
        return None;
    }
    let start = usize::from(record.alignment_start()?);
    let ref_len = reference_length(&record.cigar());
    let trailing = trailing_clipping(&cigar_to_ops(record));
    // alignment_end = start + ref_len - 1
    // unclipped_end = alignment_end + trailing_clips
    Some(start + ref_len.saturating_sub(1) + trailing)
}

/// Gets the unclipped 5' position of a read.
///
/// For forward strand reads, returns the unclipped start position.
/// For reverse strand reads, returns the unclipped end position (the 5' end).
///
/// This matches fgbio's `positionOf` behavior in `GroupReadsByUmi`.
///
/// Returns `None` for unmapped reads.
#[must_use]
pub fn unclipped_five_prime_position(record: &RecordBuf) -> Option<usize> {
    if record.flags().is_unmapped() {
        return None;
    }
    if record.flags().is_reverse_complemented() {
        unclipped_end(record)
    } else {
        unclipped_start(record)
    }
}

/// Counts reference-consuming operations from a CIGAR.
#[allow(clippy::redundant_closure_for_method_calls)]
#[must_use]
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

/// Get pair orientation using htsjdk's algorithm.
/// This matches htsjdk's `SamPairUtil.getPairOrientation()` exactly.
#[must_use]
pub fn get_pair_orientation(record: &RecordBuf) -> PairOrientation {
    let is_reverse = record.flags().is_reverse_complemented();
    let mate_reverse = record.flags().is_mate_reverse_complemented();

    // Same strand = TANDEM
    if is_reverse == mate_reverse {
        return PairOrientation::Tandem;
    }

    // Now determine if FR or RF using htsjdk's logic
    let alignment_start = record.alignment_start().map_or(0, usize::from);
    let mate_start = record.mate_alignment_start().map_or(0, usize::from);
    let insert_size = record.template_length();

    let (positive_five_prime, negative_five_prime) = if is_reverse {
        // This read is on reverse strand, mate is on positive strand
        let ref_len = reference_length(&record.cigar());
        let end = alignment_start + ref_len.saturating_sub(1);
        (mate_start as i64, end as i64)
    } else {
        // This read is on positive strand, mate is on reverse strand
        (alignment_start as i64, alignment_start as i64 + i64::from(insert_size))
    };

    // FR if positive strand 5' < negative strand 5'
    if positive_five_prime < negative_five_prime {
        PairOrientation::FR
    } else {
        PairOrientation::RF
    }
}

/// Check if a read pair is in FR (forward-reverse) orientation.
///
/// FR orientation means one read is on the forward strand and one on the reverse strand,
/// with the positive strand 5' end before the negative strand 5' end.
///
/// This matches the behavior of fgbio's `isFrPair` check.
#[must_use]
pub fn is_fr_pair(r1: &RecordBuf, r2: &RecordBuf) -> bool {
    let r1_flags = r1.flags();
    let r2_flags = r2.flags();

    // Check if both reads are paired
    if !r1_flags.is_segmented() || !r2_flags.is_segmented() {
        return false;
    }

    // Check if both reads are mapped
    if r1_flags.is_unmapped() || r2_flags.is_unmapped() {
        return false;
    }

    // Check if both mates are mapped
    if r1_flags.is_mate_unmapped() || r2_flags.is_mate_unmapped() {
        return false;
    }

    // Check if both reads are on the same reference sequence
    let r1_ref = r1.reference_sequence_id();
    let r2_ref = r2.reference_sequence_id();
    if r1_ref != r2_ref {
        return false;
    }

    // Use htsjdk's pair orientation logic (matches fgbio)
    // FR means one read is forward, one is reverse, with positive strand 5' < negative strand 5'
    // This works regardless of which read is R1 vs R2
    let orientation = get_pair_orientation(r1);
    orientation == PairOrientation::FR
}

#[cfg(test)]
#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
mod tests {
    use super::*;
    use crate::sam::builder::{RecordBuilder, RecordPairBuilder};
    use noodles::sam::alignment::record::Flags;

    // Flag constants for test readability
    const FLAG_PAIRED: u16 = 0x1;
    const FLAG_READ1: u16 = 0x40;
    const FLAG_REVERSE: u16 = 0x10;
    const FLAG_MATE_REVERSE: u16 = 0x20;
    const FLAG_UNMAPPED: u16 = 0x4;
    const FLAG_MATE_UNMAPPED: u16 = 0x8;

    // =====================================================================
    // Tests for parse_cigar_string
    // =====================================================================

    #[test]
    fn test_parse_cigar_string() {
        let ops = parse_cigar_string("10M5I20M");
        assert_eq!(ops.len(), 3);
        assert_eq!(ops[0], (Kind::Match, 10));
        assert_eq!(ops[1], (Kind::Insertion, 5));
        assert_eq!(ops[2], (Kind::Match, 20));
    }

    #[test]
    fn test_parse_cigar_string_with_clips() {
        let ops = parse_cigar_string("5H10S50M10S5H");
        assert_eq!(ops.len(), 5);
        assert_eq!(ops[0], (Kind::HardClip, 5));
        assert_eq!(ops[1], (Kind::SoftClip, 10));
        assert_eq!(ops[2], (Kind::Match, 50));
        assert_eq!(ops[3], (Kind::SoftClip, 10));
        assert_eq!(ops[4], (Kind::HardClip, 5));
    }

    #[test]
    fn test_cigar_reference_length() {
        let ops = parse_cigar_string("10M5I20M5D10M");
        assert_eq!(cigar_reference_length(&ops), 45); // 10 + 20 + 5 + 10
    }

    #[test]
    fn test_leading_clipping() {
        let ops = parse_cigar_string("5H10S50M10S5H");
        assert_eq!(leading_clipping(&ops), 15); // 5H + 10S
    }

    #[test]
    fn test_trailing_clipping() {
        let ops = parse_cigar_string("5H10S50M10S5H");
        assert_eq!(trailing_clipping(&ops), 15); // 10S + 5H
    }

    #[test]
    fn test_leading_soft_clipping() {
        let ops = parse_cigar_string("5H10S50M10S5H");
        assert_eq!(leading_soft_clipping(&ops), 10); // Only 10S, not 5H
    }

    #[test]
    fn test_trailing_soft_clipping() {
        let ops = parse_cigar_string("5H10S50M10S5H");
        assert_eq!(trailing_soft_clipping(&ops), 10); // Only 10S, not 5H
    }

    // =====================================================================
    // Tests for is_fr_pair_from_tags
    // =====================================================================

    /// Helper to create a read for FR pair tests
    fn create_fr_test_read(
        name: &str,
        flags: u16,
        pos: usize,
        mate_pos: usize,
        tlen: i32,
        ref_id: usize,
        mate_ref_id: usize,
    ) -> RecordBuf {
        let is_reverse = (flags & FLAG_REVERSE) != 0;
        let is_mate_reverse = (flags & FLAG_MATE_REVERSE) != 0;
        let is_unmapped = (flags & FLAG_UNMAPPED) != 0;
        let is_mate_unmapped = (flags & FLAG_MATE_UNMAPPED) != 0;
        let is_paired = (flags & FLAG_PAIRED) != 0;
        let is_read1 = (flags & FLAG_READ1) != 0;

        let mut builder = RecordBuilder::new()
            .name(name)
            .sequence(&"A".repeat(100))
            .qualities(&[30u8; 100])
            .cigar("100M")
            .reference_sequence_id(ref_id)
            .alignment_start(pos)
            .mate_reference_sequence_id(mate_ref_id)
            .mate_alignment_start(mate_pos)
            .template_length(tlen)
            .reverse_complement(is_reverse)
            .mate_reverse_complement(is_mate_reverse)
            .unmapped(is_unmapped)
            .mate_unmapped(is_mate_unmapped);

        if is_paired {
            builder = builder.first_segment(is_read1);
        }

        builder.build()
    }

    #[test]
    fn test_is_fr_pair_true_for_fr_orientation() {
        // FR pair: positive strand read at pos 100, mate on negative strand
        // Insert size positive (300), so negative 5' = 100 + 300 = 400
        // positive 5' (100) < negative 5' (400) → FR
        let read = create_fr_test_read(
            "fr_pair",
            FLAG_PAIRED | FLAG_READ1 | FLAG_MATE_REVERSE,
            100,
            200,
            300,
            0,
            0,
        );
        assert!(is_fr_pair_from_tags(&read), "Should be FR pair");
    }

    #[test]
    fn test_is_fr_pair_false_for_rf_orientation() {
        // RF pair: negative insert size
        let read = create_fr_test_read(
            "rf_pair",
            FLAG_PAIRED | FLAG_READ1 | FLAG_MATE_REVERSE,
            63837,
            62870,
            -967,
            0,
            0,
        );
        assert!(!is_fr_pair_from_tags(&read), "Should NOT be FR pair: RF orientation");
    }

    #[test]
    fn test_is_fr_pair_false_for_tandem() {
        // Tandem: both reads on same strand
        let read = create_fr_test_read(
            "tandem_pair",
            FLAG_PAIRED | FLAG_READ1, // both positive strand
            100,
            200,
            300,
            0,
            0,
        );
        assert!(!is_fr_pair_from_tags(&read), "Should NOT be FR pair: tandem");
    }

    #[test]
    fn test_is_fr_pair_false_for_unpaired() {
        // Create unpaired read - no flags set means not SEGMENTED
        let record = RecordBuilder::new().name("unpaired").sequence("ACGT").build();
        assert!(!is_fr_pair_from_tags(&record), "Should NOT be FR pair: unpaired");
    }

    #[test]
    fn test_is_fr_pair_false_for_unmapped() {
        let read = create_fr_test_read(
            "unmapped",
            FLAG_PAIRED | FLAG_READ1 | FLAG_UNMAPPED | FLAG_MATE_REVERSE,
            100,
            200,
            300,
            0,
            0,
        );
        assert!(!is_fr_pair_from_tags(&read), "Should NOT be FR pair: unmapped");
    }

    #[test]
    fn test_is_fr_pair_false_for_mate_unmapped() {
        let read = create_fr_test_read(
            "mate_unmapped",
            FLAG_PAIRED | FLAG_READ1 | FLAG_MATE_UNMAPPED | FLAG_MATE_REVERSE,
            100,
            200,
            300,
            0,
            0,
        );
        assert!(!is_fr_pair_from_tags(&read), "Should NOT be FR pair: mate unmapped");
    }

    #[test]
    fn test_is_fr_pair_false_for_different_references() {
        let read = create_fr_test_read(
            "diff_ref",
            FLAG_PAIRED | FLAG_READ1 | FLAG_MATE_REVERSE,
            100,
            200,
            300,
            0,
            1, // different reference
        );
        assert!(!is_fr_pair_from_tags(&read), "Should NOT be FR pair: different refs");
    }

    #[test]
    fn test_is_fr_pair_false_for_rf_when_read_on_negative_strand() {
        // RF pair: read on negative strand, mate on positive strand
        // For negative strand read:
        //   positive 5' = mate_start (100)
        //   negative 5' = alignment_end (this read, pos 200 + 100bp = 300)
        // If mate_start (100) > alignment_end... wait, that would be FR
        // Let's make mate_start > alignment_end: mate at 400, read ends at 300
        // positive 5' (400) < negative 5' (300) is FALSE → RF
        let read = create_fr_test_read(
            "rf_negative_strand",
            FLAG_PAIRED | FLAG_READ1 | FLAG_REVERSE, // negative strand, mate positive
            200,                                     // alignment_start (alignment_end will be ~300)
            400,                                     // mate_alignment_start (positive strand 5')
            -200,                                    // template_length
            0,                                       // ref_id
            0,                                       // mate_ref_id
        );

        assert!(
            !is_fr_pair_from_tags(&read),
            "Should NOT be FR pair: RF orientation (mate 5' > read 5')"
        );
    }

    #[test]
    fn test_is_fr_pair_true_for_fr_when_read_on_negative_strand() {
        // FR pair: read on negative strand, mate on positive strand
        // For negative strand read:
        //   positive 5' = mate_start
        //   negative 5' = alignment_end (this read's end)
        // FR when positive 5' < negative 5'
        // mate at 100, read ends at 300 → 100 < 300 → FR
        let read = create_fr_test_read(
            "fr_negative_strand",
            FLAG_PAIRED | FLAG_READ1 | FLAG_REVERSE, // negative strand, mate positive
            200,                                     // alignment_start (alignment_end will be ~300)
            100,                                     // mate_alignment_start (positive strand 5')
            200,                                     // template_length
            0,                                       // ref_id
            0,                                       // mate_ref_id
        );

        assert!(is_fr_pair_from_tags(&read), "Should be FR pair: mate 5' (100) < read end (300)");
    }

    #[test]
    fn test_is_fr_pair_false_for_tandem_both_reverse() {
        // Tandem: both reads on negative strand
        let read = create_fr_test_read(
            "tandem_pair_reverse",
            FLAG_PAIRED | FLAG_READ1 | FLAG_REVERSE | FLAG_MATE_REVERSE, // both reverse
            100,
            200,
            300,
            0,
            0,
        );

        assert!(
            !is_fr_pair_from_tags(&read),
            "Should NOT be FR pair: tandem orientation (both reverse)"
        );
    }

    // =====================================================================
    // Tests for read_pos_at_ref_pos
    // =====================================================================

    /// Helper to create a read with specific CIGAR for `read_pos_at_ref_pos` tests.
    /// Sequence and qualities are auto-generated from CIGAR length.
    fn create_cigar_test_read(name: &str, pos: usize, cigar: &str) -> RecordBuf {
        RecordBuilder::new()
            .name(name)
            .cigar(cigar)
            .reference_sequence_id(0)
            .alignment_start(pos)
            .first_segment(true)
            .build()
    }

    #[test]
    fn test_read_pos_at_ref_pos_simple_match() {
        // Read at position 100 with 50M CIGAR
        let read = create_cigar_test_read("simple", 100, "50M");

        assert_eq!(read_pos_at_ref_pos(&read, 100, false), 1);
        assert_eq!(read_pos_at_ref_pos(&read, 110, false), 11);
        assert_eq!(read_pos_at_ref_pos(&read, 149, false), 50);
        assert_eq!(read_pos_at_ref_pos(&read, 150, false), 0); // outside
        assert_eq!(read_pos_at_ref_pos(&read, 99, false), 0); // before
    }

    #[test]
    fn test_read_pos_at_ref_pos_with_insertion() {
        // 10M5I10M at position 100
        let read = create_cigar_test_read("insertion", 100, "10M5I10M");

        assert_eq!(read_pos_at_ref_pos(&read, 105, false), 6);
        assert_eq!(read_pos_at_ref_pos(&read, 115, false), 21); // 10 + 5 insertion + 6
    }

    #[test]
    fn test_read_pos_at_ref_pos_with_deletion() {
        // 10M5D10M at position 100
        let read = create_cigar_test_read("deletion", 100, "10M5D10M");

        // Position within deletion
        assert_eq!(read_pos_at_ref_pos(&read, 112, false), 0);
        assert_eq!(read_pos_at_ref_pos(&read, 112, true), 10); // last aligned position
        assert_eq!(read_pos_at_ref_pos(&read, 115, false), 11); // after deletion
    }

    #[test]
    fn test_read_pos_at_ref_pos_with_soft_clip() {
        // 5S10M at position 100
        let read = create_cigar_test_read("softclip", 100, "5S10M");

        assert_eq!(read_pos_at_ref_pos(&read, 100, false), 6); // after 5bp soft clip
        assert_eq!(read_pos_at_ref_pos(&read, 105, false), 11);
    }

    // =====================================================================
    // Tests for mate_unclipped_start and mate_unclipped_end
    // =====================================================================

    /// Helper to create a read with MC tag
    fn create_mc_test_read(name: &str, mate_pos: usize, mc_cigar: &str) -> RecordBuf {
        RecordBuilder::new()
            .name(name)
            .sequence(&"A".repeat(50))
            .qualities(&[30u8; 50])
            .cigar("50M")
            .reference_sequence_id(0)
            .alignment_start(100)
            .mate_reference_sequence_id(0)
            .mate_alignment_start(mate_pos)
            .first_segment(true)
            .tag("MC", mc_cigar)
            .build()
    }

    #[test]
    fn test_mate_unclipped_start_no_clipping() {
        let read = create_mc_test_read("no_clip", 200, "50M");
        assert_eq!(mate_unclipped_start(&read), Some(200));
    }

    #[test]
    fn test_mate_unclipped_start_with_soft_clip() {
        // 10bp leading soft clip: unclipped start = 200 - 10 = 190
        let read = create_mc_test_read("soft_clip", 200, "10S50M");
        assert_eq!(mate_unclipped_start(&read), Some(190));
    }

    #[test]
    fn test_mate_unclipped_start_with_hard_clip() {
        // 5H + 10S = 15bp leading clip: unclipped start = 200 - 15 = 185
        let read = create_mc_test_read("hard_clip", 200, "5H10S50M");
        assert_eq!(mate_unclipped_start(&read), Some(185));
    }

    #[test]
    fn test_mate_unclipped_end_no_clipping() {
        // refLen = 50, unclipped end = 200 + 50 - 1 = 249
        let read = create_mc_test_read("no_clip", 200, "50M");
        assert_eq!(mate_unclipped_end(&read), Some(249));
    }

    #[test]
    fn test_mate_unclipped_end_with_soft_clip() {
        // refLen = 50, 10bp trailing soft clip: unclipped end = 200 + 50 - 1 + 10 = 259
        let read = create_mc_test_read("soft_clip", 200, "50M10S");
        assert_eq!(mate_unclipped_end(&read), Some(259));
    }

    #[test]
    fn test_mate_unclipped_end_with_deletion() {
        // refLen = 25 + 5 + 25 = 55: unclipped end = 200 + 55 - 1 = 254
        let read = create_mc_test_read("deletion", 200, "25M5D25M");
        assert_eq!(mate_unclipped_end(&read), Some(254));
    }

    #[test]
    fn test_mate_unclipped_end_with_trailing_hard_clip() {
        // refLen = 50, trailing = 5S + 3H = 8: unclipped end = 200 + 50 - 1 + 8 = 257
        let read = create_mc_test_read("hard_clip", 200, "50M5S3H");
        assert_eq!(mate_unclipped_end(&read), Some(257));
    }

    #[test]
    fn test_mate_unclipped_no_mc_tag() {
        // Read without MC tag
        let record = RecordBuilder::new()
            .name("no_mc")
            .sequence("ACGT")
            .reference_sequence_id(0)
            .alignment_start(100)
            .mate_reference_sequence_id(0)
            .mate_alignment_start(200)
            .first_segment(true)
            .build();

        assert_eq!(mate_unclipped_start(&record), None);
        assert_eq!(mate_unclipped_end(&record), None);
    }

    // =====================================================================
    // Tests for get_pair_orientation (ported from htsjdk SamPairUtilTest)
    // =====================================================================

    #[test]
    fn test_get_pair_orientation_fr_normal() {
        // FR (innie): positive strand at 100, negative strand 5' further right
        // Read on + strand at 100, mate on - strand, insert size 200
        // positive 5' = 100, negative 5' = 100 + 200 = 300 → 100 < 300 → FR
        let read = create_fr_test_read(
            "fr_normal",
            FLAG_PAIRED | FLAG_READ1 | FLAG_MATE_REVERSE,
            100,
            150,
            200,
            0,
            0,
        );
        assert_eq!(get_pair_orientation(&read), PairOrientation::FR);
    }

    #[test]
    fn test_get_pair_orientation_fr_overlapping() {
        // FR overlapping: reads face each other and overlap
        let read = create_fr_test_read(
            "fr_overlap",
            FLAG_PAIRED | FLAG_READ1 | FLAG_MATE_REVERSE,
            100,
            120,
            100,
            0,
            0,
        );
        assert_eq!(get_pair_orientation(&read), PairOrientation::FR);
    }

    #[test]
    fn test_get_pair_orientation_rf_outie() {
        // RF (outie): positive strand 5' >= negative strand 5'
        // Read on + strand at 200, mate on - strand at 100, negative insert size
        // positive 5' = 200, negative 5' = 200 + (-100) = 100 → 200 >= 100 → RF
        let read = create_fr_test_read(
            "rf_outie",
            FLAG_PAIRED | FLAG_READ1 | FLAG_MATE_REVERSE,
            200,
            100,
            -100,
            0,
            0,
        );
        assert_eq!(get_pair_orientation(&read), PairOrientation::RF);
    }

    #[test]
    fn test_get_pair_orientation_tandem_both_forward() {
        // Tandem: both reads on same strand (both forward)
        let read = create_fr_test_read(
            "tandem_fwd",
            FLAG_PAIRED | FLAG_READ1, // no FLAG_MATE_REVERSE, both forward
            100,
            200,
            200,
            0,
            0,
        );
        assert_eq!(get_pair_orientation(&read), PairOrientation::Tandem);
    }

    #[test]
    fn test_get_pair_orientation_tandem_both_reverse() {
        // Tandem: both reads on same strand (both reverse)
        let read = create_fr_test_read(
            "tandem_rev",
            FLAG_PAIRED | FLAG_READ1 | FLAG_REVERSE | FLAG_MATE_REVERSE,
            100,
            200,
            200,
            0,
            0,
        );
        assert_eq!(get_pair_orientation(&read), PairOrientation::Tandem);
    }

    #[test]
    fn test_get_pair_orientation_fr_from_reverse_strand() {
        // FR when read is on reverse strand: mate on + strand
        // For reverse strand read:
        //   positive 5' = mate_start (100)
        //   negative 5' = alignment_end (~200 + 100 - 1 = 299)
        // 100 < 299 → FR
        let read = create_fr_test_read(
            "fr_reverse",
            FLAG_PAIRED | FLAG_READ1 | FLAG_REVERSE, // read reverse, mate forward
            200,
            100,
            200,
            0,
            0,
        );
        assert_eq!(get_pair_orientation(&read), PairOrientation::FR);
    }

    #[test]
    fn test_get_pair_orientation_rf_from_reverse_strand() {
        // RF when read is on reverse strand
        // For reverse strand read:
        //   positive 5' = mate_start (400)
        //   negative 5' = alignment_end (~200 + 100 - 1 = 299)
        // 400 >= 299 → RF
        let read = create_fr_test_read(
            "rf_reverse",
            FLAG_PAIRED | FLAG_READ1 | FLAG_REVERSE, // read reverse, mate forward
            200,
            400,
            -200,
            0,
            0,
        );
        assert_eq!(get_pair_orientation(&read), PairOrientation::RF);
    }

    // =====================================================================
    // Tests for is_fr_pair (two-record version)
    // =====================================================================

    /// Helper to create a pair of reads for `is_fr_pair` tests
    fn create_read_pair(
        r1_pos: usize,
        r1_reverse: bool,
        r2_pos: usize,
        r2_reverse: bool,
        same_ref: bool,
    ) -> (RecordBuf, RecordBuf) {
        let mut builder = RecordPairBuilder::new()
            .name("read_pair")
            .r1_sequence(&"A".repeat(100))
            .r2_sequence(&"A".repeat(100))
            .r1_start(r1_pos)
            .r2_start(r2_pos)
            .r1_reverse(r1_reverse)
            .r2_reverse(r2_reverse);

        if !same_ref {
            builder = builder.r2_reference_sequence_id(1);
        }

        builder.build()
    }

    #[test]
    fn test_is_fr_pair_two_records_fr() {
        // FR pair: R1 forward at 100, R2 reverse at 200
        let (r1, r2) = create_read_pair(100, false, 200, true, true);
        assert!(is_fr_pair(&r1, &r2), "Should detect FR pair");
    }

    #[test]
    fn test_is_fr_pair_two_records_rf() {
        // RF pair: R1 reverse at 100, R2 forward at 50 (outie)
        // This creates an RF orientation where positive 5' > negative 5'
        let (r1, r2) = create_read_pair(200, false, 100, true, true);
        // With these positions, positive 5' = 200, and with negative tlen,
        // negative 5' = 200 + (-50) = 150, so 200 > 150 → RF
        let mut r1_rf = r1;
        *r1_rf.template_length_mut() = -50;
        assert!(!is_fr_pair(&r1_rf, &r2), "Should NOT detect FR pair for RF orientation");
    }

    #[test]
    fn test_is_fr_pair_two_records_tandem() {
        // Tandem: both forward
        let (r1, r2) = create_read_pair(100, false, 200, false, true);
        assert!(!is_fr_pair(&r1, &r2), "Should NOT detect FR pair for tandem");
    }

    #[test]
    fn test_is_fr_pair_two_records_different_refs() {
        // Different references
        let (r1, r2) = create_read_pair(100, false, 200, true, false);
        assert!(!is_fr_pair(&r1, &r2), "Should NOT detect FR pair for different refs");
    }

    #[test]
    fn test_is_fr_pair_two_records_unpaired() {
        // Create unpaired reads - no paired flags set
        let r1 = RecordBuilder::new()
            .name("unpaired")
            .sequence("ACGT")
            .reference_sequence_id(0)
            .alignment_start(100)
            .build();

        let r2 = RecordBuilder::new()
            .name("unpaired")
            .sequence("ACGT")
            .reference_sequence_id(0)
            .alignment_start(200)
            .build();

        assert!(!is_fr_pair(&r1, &r2), "Should NOT detect FR pair for unpaired reads");
    }

    #[test]
    fn test_is_fr_pair_two_records_unmapped() {
        // Create pair with unmapped read
        let (mut r1, r2) = create_read_pair(100, false, 200, true, true);
        *r1.flags_mut() = Flags::from(FLAG_PAIRED | FLAG_READ1 | FLAG_UNMAPPED | FLAG_MATE_REVERSE);
        assert!(!is_fr_pair(&r1, &r2), "Should NOT detect FR pair for unmapped read");
    }

    // =====================================================================
    // Tests for edge cases
    // =====================================================================

    #[test]
    fn test_parse_cigar_string_empty() {
        let ops = parse_cigar_string("");
        assert!(ops.is_empty());
    }

    #[test]
    fn test_parse_cigar_string_all_operations() {
        // Test all CIGAR operation types
        let ops = parse_cigar_string("10M5I3D2N7S4H1P6=8X");
        assert_eq!(ops.len(), 9);
        assert_eq!(ops[0], (Kind::Match, 10));
        assert_eq!(ops[1], (Kind::Insertion, 5));
        assert_eq!(ops[2], (Kind::Deletion, 3));
        assert_eq!(ops[3], (Kind::Skip, 2));
        assert_eq!(ops[4], (Kind::SoftClip, 7));
        assert_eq!(ops[5], (Kind::HardClip, 4));
        assert_eq!(ops[6], (Kind::Pad, 1));
        assert_eq!(ops[7], (Kind::SequenceMatch, 6));
        assert_eq!(ops[8], (Kind::SequenceMismatch, 8));
    }

    #[test]
    fn test_cigar_reference_length_with_all_ops() {
        // Reference-consuming ops: M, D, N, =, X
        // Non-consuming: I, S, H, P
        let ops = parse_cigar_string("10M5I3D2N7S4H1P6=8X");
        // Reference length = 10(M) + 3(D) + 2(N) + 6(=) + 8(X) = 29
        assert_eq!(cigar_reference_length(&ops), 29);
    }

    #[test]
    fn test_leading_clipping_only_hard_clip() {
        let ops = parse_cigar_string("10H50M");
        assert_eq!(leading_clipping(&ops), 10);
        assert_eq!(leading_soft_clipping(&ops), 0);
    }

    #[test]
    fn test_trailing_clipping_only_hard_clip() {
        let ops = parse_cigar_string("50M10H");
        assert_eq!(trailing_clipping(&ops), 10);
        assert_eq!(trailing_soft_clipping(&ops), 0);
    }

    #[test]
    fn test_alignment_end_simple() {
        let record = RecordBuilder::new()
            .sequence(&"A".repeat(50))
            .alignment_start(100)
            .cigar("50M")
            .build();
        // end = 100 + 50 - 1 = 149
        assert_eq!(alignment_end(&record), Some(149));
    }

    #[test]
    fn test_alignment_end_with_deletion() {
        let record = RecordBuilder::new()
            .sequence(&"A".repeat(50))
            .alignment_start(100)
            .cigar("25M5D25M")
            .build();
        // ref_len = 25 + 5 + 25 = 55, end = 100 + 55 - 1 = 154
        assert_eq!(alignment_end(&record), Some(154));
    }

    #[test]
    fn test_alignment_end_no_alignment_start() {
        // Record with no alignment start set
        let record = RecordBuilder::new().sequence("ACGT").build();
        assert_eq!(alignment_end(&record), None);
    }

    #[test]
    fn test_read_pos_at_ref_pos_complex_cigar() {
        // Test with skip (N) operation: 10M5N10M at position 100
        let read = create_cigar_test_read("skip", 100, "10M5N10M");

        // Position 105 is within first match block → read pos 6
        assert_eq!(read_pos_at_ref_pos(&read, 105, false), 6);
        // Position 112 is within skip → returns 0 (not in deletion)
        assert_eq!(read_pos_at_ref_pos(&read, 112, false), 0);
        // Position 115 is in second match block → read pos 11
        assert_eq!(read_pos_at_ref_pos(&read, 115, false), 11);
    }

    #[test]
    fn test_read_pos_at_ref_pos_with_sequence_match_mismatch() {
        // Test with =/X operations: 5=3X5= at position 100
        let read = create_cigar_test_read("eq_x", 100, "5=3X5=");

        assert_eq!(read_pos_at_ref_pos(&read, 100, false), 1);
        assert_eq!(read_pos_at_ref_pos(&read, 105, false), 6); // in mismatch block
        assert_eq!(read_pos_at_ref_pos(&read, 108, false), 9); // in second match block
    }

    // =========================================================================
    // Tests for unclipped_start, unclipped_end, unclipped_five_prime_position
    // These match HTSJDK/fgbio behavior including both soft AND hard clips
    // =========================================================================

    #[test]
    fn test_unclipped_start_no_clips() {
        // Simple 50M at position 100 → unclipped_start = 100
        let read = create_cigar_test_read("simple", 100, "50M");
        assert_eq!(unclipped_start(&read), Some(100));
    }

    #[test]
    fn test_unclipped_start_with_leading_soft_clip() {
        // 5S45M at position 100 → unclipped_start = 100 - 5 = 95
        let read = create_cigar_test_read("soft", 100, "5S45M");
        assert_eq!(unclipped_start(&read), Some(95));
    }

    #[test]
    fn test_unclipped_start_with_leading_hard_clip() {
        // 10H50M at position 100 → unclipped_start = 100 - 10 = 90
        let read = create_cigar_test_read("hard", 100, "10H50M");
        assert_eq!(unclipped_start(&read), Some(90));
    }

    #[test]
    fn test_unclipped_start_with_soft_and_hard_clips() {
        // Matches fgbio test: 5S45M10H at position 10 → unclipped_start = 10 - 5 = 5
        // (trailing hard clip doesn't affect unclipped_start)
        let read = create_cigar_test_read("both", 10, "5S45M10H");
        assert_eq!(unclipped_start(&read), Some(5));
    }

    #[test]
    fn test_unclipped_start_with_leading_hard_and_soft() {
        // 3H5S42M at position 100 → unclipped_start = 100 - 5 - 3 = 92
        let read = create_cigar_test_read("hard_soft", 100, "3H5S42M");
        assert_eq!(unclipped_start(&read), Some(92));
    }

    #[test]
    fn test_unclipped_end_no_clips() {
        // Simple 50M at position 100 → alignment_end = 149, unclipped_end = 149
        let read = create_cigar_test_read("simple", 100, "50M");
        assert_eq!(unclipped_end(&read), Some(149));
    }

    #[test]
    fn test_unclipped_end_with_trailing_soft_clip() {
        // 45M5S at position 100 → alignment_end = 144, unclipped_end = 144 + 5 = 149
        let read = create_cigar_test_read("soft", 100, "45M5S");
        assert_eq!(unclipped_end(&read), Some(149));
    }

    #[test]
    fn test_unclipped_end_with_trailing_hard_clip() {
        // 50M10H at position 100 → alignment_end = 149, unclipped_end = 149 + 10 = 159
        let read = create_cigar_test_read("hard", 100, "50M10H");
        assert_eq!(unclipped_end(&read), Some(159));
    }

    #[test]
    fn test_unclipped_end_with_soft_and_hard_clips() {
        // Matches fgbio test: 5S45M10H at position 10
        // alignment_end = 10 + 45 - 1 = 54
        // unclipped_end = 54 + 10 = 64
        let read = create_cigar_test_read("both", 10, "5S45M10H");
        assert_eq!(unclipped_end(&read), Some(64));
    }

    #[test]
    fn test_unclipped_end_with_trailing_soft_and_hard() {
        // 42M5S3H at position 100 → alignment_end = 141, unclipped_end = 141 + 5 + 3 = 149
        let read = create_cigar_test_read("soft_hard", 100, "42M5S3H");
        assert_eq!(unclipped_end(&read), Some(149));
    }

    #[test]
    fn test_unclipped_five_prime_forward_strand() {
        // Forward strand: 5' is at unclipped_start
        // 5S45M10H at position 10 → unclipped_start = 5
        let read = create_cigar_test_read("fwd", 10, "5S45M10H");
        assert_eq!(unclipped_five_prime_position(&read), Some(5));
    }

    #[test]
    fn test_unclipped_five_prime_reverse_strand() {
        // Reverse strand: 5' is at unclipped_end
        // 5S45M10H at position 10 → unclipped_end = 64
        let read = RecordBuilder::new()
            .name("rev")
            .sequence(&"A".repeat(60))
            .reference_sequence_id(0)
            .alignment_start(10)
            .cigar("5S45M10H")
            .flags(Flags::REVERSE_COMPLEMENTED)
            .build();
        assert_eq!(unclipped_five_prime_position(&read), Some(64));
    }

    #[test]
    fn test_unclipped_positions_unmapped_returns_none() {
        let read =
            RecordBuilder::new().name("unmapped").sequence("ACGT").flags(Flags::UNMAPPED).build();
        assert_eq!(unclipped_start(&read), None);
        assert_eq!(unclipped_end(&read), None);
        assert_eq!(unclipped_five_prime_position(&read), None);
    }

    #[test]
    fn test_unclipped_end_with_deletion() {
        // 25M5D25M at position 100
        // ref_len = 25 + 5 + 25 = 55
        // alignment_end = 100 + 55 - 1 = 154
        let read = create_cigar_test_read("del", 100, "25M5D25M");
        assert_eq!(unclipped_end(&read), Some(154));
    }

    #[test]
    fn test_unclipped_end_with_insertion() {
        // 25M5I25M at position 100
        // ref_len = 25 + 25 = 50 (insertion doesn't consume reference)
        // alignment_end = 100 + 50 - 1 = 149
        let read = create_cigar_test_read("ins", 100, "25M5I25M");
        assert_eq!(unclipped_end(&read), Some(149));
    }
}

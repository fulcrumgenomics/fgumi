//! Record-level utilities for SAM/BAM records.
//!
//! This module provides utilities for working with individual SAM records, including:
//! - Position mapping between read and reference coordinates
//! - FR pair detection using BAM tags
//! - Mate position calculations from MC tag
//! - CIGAR string parsing and analysis

use fgumi_raw_bam::{
    RawRecordView, alignment_start_from_raw, consumes_query, consumes_ref, find_mc_tag_in_record,
    flags, get_cigar_ops, mate_pos, parse_cigar_string_ops,
};
use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::record::Cigar as CigarTrait;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value;

/// Pair orientation for paired-end reads (defined in `fgumi-raw-bam`, re-exported here).
pub use fgumi_raw_bam::PairOrientation;

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
// RecordBuf kept: uses noodles typed CIGAR iterator (record.cigar().iter() -> Result<Op>); raw equivalent is read_pos_at_ref_pos_raw in fgumi-raw-bam.
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
/// - The pair is in FR orientation (positive strand 5' `<=` negative strand 5', per htsjdk 5.0.0)
///
/// This matches fgbio's per-record `isFrPair` using htsjdk's pair orientation logic (see
/// [`get_pair_orientation`]: on a forward-strand read the mate's end comes from `MC` when present,
/// else from `TLEN`). When both records are available, use [`is_fr_pair`].
// RecordBuf kept: accesses flags/ref-id/mate fields via noodles typed API; raw equivalent is is_fr_pair_raw in fgumi-raw-bam.
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
    if read.reference_sequence_id() != read.mate_reference_sequence_id() {
        return false;
    }

    // Delegate to get_pair_orientation for the FR/RF/Tandem determination
    get_pair_orientation(read) == PairOrientation::FR
}

/// Parses the MC (mate CIGAR) tag from a record, returning the parsed CIGAR ops
/// and the mate alignment start position.
///
/// Returns `None` if the MC tag is missing/invalid or mate position is missing.
// RecordBuf kept: reads MC tag via noodles typed data map (record.data().get()) returning Value::String.
fn parse_mate_cigar(read: &RecordBuf) -> Option<(Vec<(Kind, usize)>, usize)> {
    let mc_tag = Tag::from([b'M', b'C']);
    let mc_value = read.data().get(&mc_tag)?;

    let cigar_bytes = match mc_value {
        // Parse the raw bytes directly rather than gating on UTF-8. A valid
        // CIGAR prefix followed by non-UTF-8 trailing bytes is still usable, and
        // this keeps the typed path in parity with the raw-byte sibling
        // (`mate_unclipped_*_raw`), which parses the same bytes.
        Value::String(s) => s.as_ref(),
        _ => return None,
    };

    let mate_start = usize::from(read.mate_alignment_start()?);
    let ops = parse_cigar_bytes(cigar_bytes);
    if ops.is_empty() {
        return None;
    }
    Some((ops, mate_start))
}

/// Gets the mate's unclipped start position from the MC tag and mate position.
///
/// Calculates: `mateUnclippedStart = mateStart - leadingClipping`
///
/// Returns the mate's unclipped start as a signed integer.
///
/// This uses signed arithmetic to properly handle cases where the unclipped start
/// would be before position 1 (i.e., negative in 1-based coordinates).
/// This matches fgbio's behavior which uses signed integers for these calculations.
///
/// Returns `None` if the MC tag is missing or invalid, or if mate position is missing.
// RecordBuf kept: delegates to parse_mate_cigar which uses noodles typed data API.
#[must_use]
pub fn mate_unclipped_start(read: &RecordBuf) -> Option<isize> {
    let (ops, mate_start) = parse_mate_cigar(read)?;

    #[expect(
        clippy::cast_possible_wrap,
        reason = "genomic positions fit in isize on all supported platforms (64-bit)"
    )]
    let mate_start_signed = mate_start as isize;
    #[expect(
        clippy::cast_possible_wrap,
        reason = "CIGAR clipping lengths are small relative to isize::MAX"
    )]
    let leading_clip = leading_clipping(&ops) as isize;

    Some(mate_start_signed - leading_clip)
}

/// Gets the mate's unclipped end position from the MC tag and mate position.
///
/// Calculates: `mateUnclippedEnd = mateStart + refLen - 1 + trailingClipping`
///
/// Returns `None` if the MC tag is missing or invalid, or if mate position is missing.
// RecordBuf kept: delegates to parse_mate_cigar which uses noodles typed data API.
#[must_use]
pub fn mate_unclipped_end(read: &RecordBuf) -> Option<usize> {
    let (ops, mate_start) = parse_mate_cigar(read)?;
    let ref_len = cigar_reference_length(&ops);
    let trailing_clip = trailing_clipping(&ops);

    // mateStart + refLen - 1 + trailingClipping
    Some(mate_start + ref_len - 1 + trailing_clip)
}

/// Gets the read's alignment end position (1-based, inclusive).
///
/// Calculated as: `alignment_start + reference_length - 1`
// RecordBuf kept: reads alignment_start via noodles typed Position and uses typed CIGAR; raw equivalent is alignment_end_from_raw in fgumi-raw-bam.
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
    parse_cigar_bytes(cigar_str.as_bytes())
}

/// Parses a CIGAR from raw bytes, as [`parse_cigar_string`] does from `&str`.
///
/// `MC` is extracted as bytes (`find_mc_tag_in_record`) so that the sort's
/// per-record ingest path does not pay a UTF-8 validation it cannot use, and
/// this is the entry point that consumes it. Operators are ASCII by
/// construction, so a non-ASCII byte falls into the same "unknown operator"
/// arm a non-ASCII `char` did: the two functions agree on every input.
#[must_use]
pub fn parse_cigar_bytes(cigar: &[u8]) -> Vec<(Kind, usize)> {
    let mut ops = Vec::new();
    let mut num_str = String::new();

    for &b in cigar {
        if b.is_ascii_digit() {
            num_str.push(char::from(b));
        } else {
            let len: usize = num_str.parse().unwrap_or(0);
            num_str.clear();

            let kind = match b {
                b'M' => Kind::Match,
                b'I' => Kind::Insertion,
                b'D' => Kind::Deletion,
                b'N' => Kind::Skip,
                b'S' => Kind::SoftClip,
                b'H' => Kind::HardClip,
                b'P' => Kind::Pad,
                b'=' => Kind::SequenceMatch,
                b'X' => Kind::SequenceMismatch,
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
// RecordBuf kept: accesses record.cigar() via noodles typed CigarBuf to extract ops.
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
// RecordBuf kept: checks record.flags() and alignment_start() via noodles typed API; raw equivalent is unclipped_5prime_from_raw_bam in fgumi-raw-bam.
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
// RecordBuf kept: checks record.flags() and alignment_start() via noodles typed API.
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

/// Gets the un-*soft*-clipped start position of a read (alignment start minus leading
/// **soft** clips only; hard clips are ignored because the hard-clipped bases are
/// physically absent from the read).
///
/// This matches fgbio's `SamRecord.unSoftClippedStart`, used when bounding a mate's
/// window for `clipExtendingPastMate` (see fgbio `SamRecordClipper`). Contrast
/// [`unclipped_start`], which includes hard clips (htsjdk `getUnclippedStart`).
///
/// Returns `None` for unmapped reads.
#[must_use]
pub fn unsoftclipped_start(record: &RecordBuf) -> Option<usize> {
    if record.flags().is_unmapped() {
        return None;
    }
    let start = usize::from(record.alignment_start()?);
    let leading = leading_soft_clipping(&cigar_to_ops(record));
    Some(start.saturating_sub(leading))
}

/// Gets the un-*soft*-clipped end position of a read (alignment end plus trailing
/// **soft** clips only; hard clips are ignored because the hard-clipped bases are
/// physically absent from the read).
///
/// This matches fgbio's `SamRecord.unSoftClippedEnd`. Contrast [`unclipped_end`],
/// which includes hard clips (htsjdk `getUnclippedEnd`).
///
/// Returns `None` for unmapped reads or records with no CIGAR ops (matching the raw
/// sibling [`unsoftclipped_end_raw`]).
#[must_use]
pub fn unsoftclipped_end(record: &RecordBuf) -> Option<usize> {
    if record.flags().is_unmapped() {
        return None;
    }
    let start = usize::from(record.alignment_start()?);
    let ops = cigar_to_ops(record);
    if ops.is_empty() {
        return None;
    }
    let ref_len = reference_length(&record.cigar());
    let trailing = trailing_soft_clipping(&ops);
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
// RecordBuf kept: delegates to unclipped_start/unclipped_end which use noodles typed API; raw equivalent is unclipped_5prime_from_raw_bam in fgumi-raw-bam.
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

/// Sums the lengths of the reference-consuming operations of a CIGAR, saturating at
/// `usize::MAX` rather than overflowing on an adversarial CIGAR.
// RecordBuf kept: parameterized over &impl CigarTrait (noodles trait) to accept both CigarBuf and inline Cigar; raw equivalent is reference_length_from_raw_bam in fgumi-raw-bam.
#[expect(
    clippy::redundant_closure_for_method_calls,
    reason = "Op::len is not a method on the trait"
)]
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
        .fold(0usize, usize::saturating_add)
}

/// htsjdk's `getAlignmentEnd()` for a 1-based `start` and reference span `ref_len`:
/// `start + ref_len - 1`, unclamped (a zero span ends at `start - 1`), in `i64`. A span or end
/// too large for `i64` saturates at `i64::MAX` rather than wrapping, so it stays past any real
/// position.
fn alignment_end_i64(start: usize, ref_len: usize) -> i64 {
    let start = i64::try_from(start).unwrap_or(i64::MAX);
    let ref_len = i64::try_from(ref_len).unwrap_or(i64::MAX);
    start.saturating_add(ref_len).saturating_sub(1)
}

/// The mate's inclusive alignment end from the `MC` tag and mate position, htsjdk's
/// `CoordMath.getEnd(mateStart, mateCigar.getReferenceLength())` = `mateStart + refLen - 1`,
/// unclamped (a zero span ends at `mateStart - 1`).
///
/// Returns `None` when `MC` or the mate position is missing, and when `MC` is not a whole, valid
/// CIGAR for a mapped record per `fgumi_raw_bam::parse_cigar_string_ops` (a dangling run length
/// like `100M5`, an unknown operator, a zero or oversized run length, a misplaced clip, or no
/// reference-consuming operation). Unlike the tolerant [`parse_mate_cigar`], which keeps a valid
/// prefix for the mate-coordinate helpers, a malformed tag here is treated as unusable so it
/// cannot override `TLEN`.
// RecordBuf kept: reads the MC tag via the noodles data map.
fn mate_alignment_end_from_mc(record: &RecordBuf) -> Option<i64> {
    let Value::String(mc) = record.data().get(&Tag::from([b'M', b'C']))? else {
        return None;
    };
    let mate_start = usize::from(record.mate_alignment_start()?);
    let mate_ops = parse_cigar_string_ops(std::str::from_utf8(mc.as_ref()).ok()?)?;
    // Each run length is at most 2^28 - 1, so the checked sum fails only on an absurd op count.
    let mate_ref_len = mate_ops
        .iter()
        .filter(|&&op| consumes_ref(op & 0xF))
        .try_fold(0i64, |acc, &op| acc.checked_add(i64::from(op >> 4)))?;
    i64::try_from(mate_start).ok()?.checked_add(mate_ref_len)?.checked_sub(1)
}

/// Get pair orientation using htsjdk's algorithm, from this record alone.
///
/// This matches htsjdk 5.0.0's `SamPairUtil.getPairOrientation()` (htsjdk#1771): the pair is FR
/// iff the positive-strand 5' `<=` the negative-strand 5' (`SamPairUtil.java:136`), so a 5' tie
/// is FR. On the reverse-strand record the negative-strand 5' is its own alignment end, htsjdk's
/// `getAlignmentEnd()` (`start + referenceLength - 1`, so a zero reference span ends at
/// `start - 1`). On the forward-strand record it is the mate's alignment end: from the `MC` tag
/// when present (`CoordMath.getEnd(mateStart, mateReferenceLength)`), else from `TLEN`
/// (`CoordMath.getEnd(start, TLEN)`). An `MC` that is not a whole, valid CIGAR (see
/// `mate_alignment_end_from_mc`) is treated as absent.
///
/// Being per-record, this can still disagree with the mate's actual alignment when `MC` is
/// absent or the mate fields are stale. When both records are available, use [`is_fr_pair`].
// RecordBuf kept: reads flags/positions via noodles typed Flags and Position API; serves noodles-typed callers.
#[must_use]
pub fn get_pair_orientation(record: &RecordBuf) -> PairOrientation {
    let is_reverse = record.flags().is_reverse_complemented();
    let mate_reverse = record.flags().is_mate_reverse_complemented();

    // Same strand = TANDEM
    if is_reverse == mate_reverse {
        return PairOrientation::Tandem;
    }

    let alignment_start = record.alignment_start().map_or(0, usize::from);
    let mate_start = record.mate_alignment_start().map_or(0, usize::from);
    let insert_size = record.template_length();

    // Genomic positions always fit in i64
    #[expect(
        clippy::cast_possible_wrap,
        reason = "genomic positions are non-negative and fit in i64 on all supported platforms"
    )]
    let (positive_five_prime, negative_five_prime) = if is_reverse {
        // This read is on reverse strand, mate is on positive strand. htsjdk's
        // `getAlignmentEnd()` is `start + referenceLength - 1`, deliberately unclamped: a read
        // with no reference-consuming ops ends at `start - 1`, so at its mate's start it is RF.
        let ref_len = reference_length(&record.cigar());
        (mate_start as i64, alignment_end_i64(alignment_start, ref_len))
    } else if let Some(mate_end) = mate_alignment_end_from_mc(record) {
        // Forward read with a usable `MC` tag: the mate's alignment end, htsjdk's
        // `CoordMath.getEnd(mateStart, mateCigar.getReferenceLength())`.
        (alignment_start as i64, mate_end)
    } else {
        // Forward read without a usable `MC` (absent or malformed): htsjdk
        // `CoordMath.getEnd(alignmentStart, insertSize)` = start + insertSize - 1. With `<=` below this classifies exactly as the former
        // `start + insertSize` with `<` did (`a <= b - 1` iff `a < b`).
        (alignment_start as i64, alignment_start as i64 + i64::from(insert_size) - 1)
    };

    // FR if positive strand 5' <= negative strand 5' (htsjdk 5.0.0): a 5' tie is FR.
    if positive_five_prime <= negative_five_prime {
        PairOrientation::FR
    } else {
        PairOrientation::RF
    }
}

/// Check if a read pair is in FR (forward-reverse) orientation.
///
/// Pair-level and order-independent: both records must be paired, mapped with mapped mates, on
/// the same reference and on opposite strands, and the pair is FR iff the positive-strand
/// record's alignment start `<=` the negative-strand record's alignment end (htsjdk 5.0.0
/// `SamPairUtil.getPairOrientation`, so a 5' tie is FR). Both positions come from the records'
/// own alignments, never from `TLEN`, the `MC` tag or the mate-position fields, so the answer is
/// the same whichever record is passed first. This is the `RecordBuf` counterpart of
/// `fgumi_raw_bam::is_primary_fr_pair_raw` and matches fgbio's `SamRecord.isFrPair(r1, r2)`.
// RecordBuf kept: reads flags/ref-ids via noodles typed Flags API; serves noodles-typed callers alongside raw is_primary_fr_pair_raw in fgumi-raw-bam.
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
    if r1.reference_sequence_id() != r2.reference_sequence_id() {
        return false;
    }

    // Opposite strands
    let r1_reverse = r1_flags.is_reverse_complemented();
    if r1_reverse == r2_flags.is_reverse_complemented() {
        return false;
    }

    // FR iff positive-strand start <= negative-strand alignment end (unclamped, as htsjdk's
    // `getAlignmentEnd()`), both from the records themselves, in i64 (the end saturates).
    // A missing start is not coordinate 0 (SAM's POS 0 means "unavailable"), so no FR pair.
    let (positive, negative) = if r1_reverse { (r2, r1) } else { (r1, r2) };
    let (Some(positive_start), Some(negative_start)) =
        (positive.alignment_start(), negative.alignment_start())
    else {
        return false;
    };
    let positive_start = i64::try_from(usize::from(positive_start)).unwrap_or(i64::MAX);
    let negative_end =
        alignment_end_i64(usize::from(negative_start), reference_length(&negative.cigar()));
    positive_start <= negative_end
}

// ============================================================================
// Raw-byte siblings — take `&[u32]` CIGAR ops or a full `&[u8]` BAM record
// ============================================================================

/// Calculates leading clipping (soft + hard) from raw u32 CIGAR operations.
///
/// This is the raw-byte sibling of [`leading_clipping`] that operates directly
/// on the packed BAM CIGAR op array returned by [`get_cigar_ops`], avoiding
/// allocation of a typed `(Kind, usize)` vector.
#[must_use]
pub fn leading_clipping_raw(cigar_ops: &[u32]) -> usize {
    cigar_ops
        .iter()
        .take_while(|&&op| {
            let t = op & 0xF;
            t == 4 || t == 5 // SoftClip or HardClip
        })
        .map(|&op| (op >> 4) as usize)
        .sum()
}

/// Calculates trailing clipping (soft + hard) from raw u32 CIGAR operations.
///
/// This is the raw-byte sibling of [`trailing_clipping`].
#[must_use]
pub fn trailing_clipping_raw(cigar_ops: &[u32]) -> usize {
    cigar_ops
        .iter()
        .rev()
        .take_while(|&&op| {
            let t = op & 0xF;
            t == 4 || t == 5 // SoftClip or HardClip
        })
        .map(|&op| (op >> 4) as usize)
        .sum()
}

/// Calculates leading soft clipping only from raw u32 CIGAR operations.
///
/// This is the raw-byte sibling of [`leading_soft_clipping`]. Hard clips at
/// the start are skipped; only the first soft-clip run is counted.
#[must_use]
pub fn leading_soft_clipping_raw(cigar_ops: &[u32]) -> usize {
    cigar_ops
        .iter()
        .skip_while(|&&op| (op & 0xF) == 5) // skip leading HardClip
        .take_while(|&&op| (op & 0xF) == 4) // count SoftClip run
        .map(|&op| (op >> 4) as usize)
        .sum()
}

/// Calculates trailing soft clipping only from raw u32 CIGAR operations.
///
/// This is the raw-byte sibling of [`trailing_soft_clipping`]. Hard clips at
/// the end are skipped; only the last soft-clip run is counted.
#[must_use]
pub fn trailing_soft_clipping_raw(cigar_ops: &[u32]) -> usize {
    cigar_ops
        .iter()
        .rev()
        .skip_while(|&&op| (op & 0xF) == 5) // skip trailing HardClip
        .take_while(|&&op| (op & 0xF) == 4) // count SoftClip run
        .map(|&op| (op >> 4) as usize)
        .sum()
}

/// Gets the mate's unclipped start position from a raw BAM record's MC tag and mate position.
///
/// This is the raw-byte sibling of [`mate_unclipped_start`].
///
/// Returns `None` if the MC tag is missing/invalid or mate position is missing.
#[must_use]
pub fn mate_unclipped_start_raw(bam: &[u8]) -> Option<isize> {
    let mc_cigar = find_mc_tag_in_record(bam)?;
    let mate_start_0based = mate_pos(bam);
    if mate_start_0based < 0 {
        return None;
    }
    let ops = parse_cigar_bytes(mc_cigar);
    if ops.is_empty() {
        return None;
    }
    // Convert 0-based BAM pos to 1-based, then subtract leading clips
    let mate_start_1based = mate_start_0based + 1;
    let leading = leading_clipping(&ops);
    // Use isize to handle the rare case where unclipped start would be before pos 1
    #[expect(
        clippy::cast_possible_wrap,
        reason = "genomic positions fit in isize on all supported platforms (64-bit)"
    )]
    Some(mate_start_1based as isize - leading as isize)
}

/// Gets the mate's unclipped end position from a raw BAM record's MC tag and mate position.
///
/// This is the raw-byte sibling of [`mate_unclipped_end`].
///
/// Returns `None` if the MC tag is missing/invalid or mate position is missing.
#[must_use]
pub fn mate_unclipped_end_raw(bam: &[u8]) -> Option<usize> {
    let mc_cigar = find_mc_tag_in_record(bam)?;
    let mate_start_0based = mate_pos(bam);
    if mate_start_0based < 0 {
        return None;
    }
    // mate_start_0based >= 0 is verified above, so the cast is safe.
    #[expect(
        clippy::cast_sign_loss,
        reason = "mate_start_0based is verified non-negative by the guard above"
    )]
    let mate_start_1based = (mate_start_0based + 1) as usize;
    let ops = parse_cigar_bytes(mc_cigar);
    if ops.is_empty() {
        return None;
    }
    let ref_len = cigar_reference_length(&ops);
    let trailing = trailing_clipping(&ops);
    Some(mate_start_1based + ref_len - 1 + trailing)
}

/// Gets the unclipped start position of a read from raw BAM bytes.
///
/// This is the raw-byte sibling of [`unclipped_start`].
/// Returns `None` for unmapped reads (no alignment start).
#[must_use]
pub fn unclipped_start_raw(bam: &[u8]) -> Option<usize> {
    if RawRecordView::new(bam).flags() & flags::UNMAPPED != 0 {
        return None;
    }
    let start = alignment_start_from_raw(bam)?;
    let cigar_ops = get_cigar_ops(bam);
    let leading = leading_clipping_raw(&cigar_ops);
    Some(start.saturating_sub(leading))
}

/// Gets the unclipped end position of a read from raw BAM bytes.
///
/// This is the raw-byte sibling of [`unclipped_end`].
/// Returns `None` for unmapped reads or records with no CIGAR ops.
#[must_use]
pub fn unclipped_end_raw(bam: &[u8]) -> Option<usize> {
    if RawRecordView::new(bam).flags() & flags::UNMAPPED != 0 {
        return None;
    }
    let start = alignment_start_from_raw(bam)?;
    let cigar_ops = get_cigar_ops(bam);
    if cigar_ops.is_empty() {
        return None;
    }
    let ref_len = cigar_reference_length_raw(&cigar_ops);
    let trailing = trailing_clipping_raw(&cigar_ops);
    // alignment_end = start + ref_len - 1; unclipped_end = alignment_end + trailing
    Some(start + ref_len.saturating_sub(1) + trailing)
}

/// Raw-byte sibling of [`unsoftclipped_start`]: alignment start minus leading **soft**
/// clips only (hard clips ignored). Matches fgbio `SamRecord.unSoftClippedStart`.
///
/// Returns `None` for unmapped reads.
#[must_use]
pub fn unsoftclipped_start_raw(bam: &[u8]) -> Option<usize> {
    if RawRecordView::new(bam).flags() & flags::UNMAPPED != 0 {
        return None;
    }
    let start = alignment_start_from_raw(bam)?;
    let cigar_ops = get_cigar_ops(bam);
    let leading = leading_soft_clipping_raw(&cigar_ops);
    Some(start.saturating_sub(leading))
}

/// Raw-byte sibling of [`unsoftclipped_end`]: alignment end plus trailing **soft**
/// clips only (hard clips ignored). Matches fgbio `SamRecord.unSoftClippedEnd`.
///
/// Returns `None` for unmapped reads or records with no CIGAR ops.
#[must_use]
pub fn unsoftclipped_end_raw(bam: &[u8]) -> Option<usize> {
    if RawRecordView::new(bam).flags() & flags::UNMAPPED != 0 {
        return None;
    }
    let start = alignment_start_from_raw(bam)?;
    let cigar_ops = get_cigar_ops(bam);
    if cigar_ops.is_empty() {
        return None;
    }
    let ref_len = cigar_reference_length_raw(&cigar_ops);
    let trailing = trailing_soft_clipping_raw(&cigar_ops);
    Some(start + ref_len.saturating_sub(1) + trailing)
}

/// Returns the reference length consumed by raw u32 CIGAR operations.
///
/// This is the raw-byte sibling of [`cigar_reference_length`].
/// Uses the packed BAM CIGAR op format (low 4 bits = op type, high 28 bits = length).
#[must_use]
pub fn cigar_reference_length_raw(cigar_ops: &[u32]) -> usize {
    cigar_ops
        .iter()
        .map(|&op| {
            let t = op & 0xF;
            let len = (op >> 4) as usize;
            if consumes_ref(t) { len } else { 0 }
        })
        .sum()
}

/// Returns the query length consumed by raw u32 CIGAR operations (includes soft clips).
///
/// Uses the packed BAM CIGAR op format. This counts M/I/S/=/X operations.
#[must_use]
pub fn cigar_query_length_raw(cigar_ops: &[u32]) -> usize {
    cigar_ops
        .iter()
        .map(|&op| {
            let t = op & 0xF;
            let len = (op >> 4) as usize;
            if consumes_query(t) { len } else { 0 }
        })
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::builder::{RecordBuilder, RecordPairBuilder};
    use noodles::sam::alignment::record::Flags;
    use rstest::rstest;

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
        // Insert size positive (300), so negative 5' = 100 + 300 - 1 = 399
        // positive 5' (100) <= negative 5' (399) -> FR
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
        // positive 5' (400) < negative 5' (300) is FALSE -> RF
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
        // FR when positive 5' <= negative 5'
        // mate at 100, read ends at 300 -> 100 <= 300 -> FR
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

    #[test]
    fn test_mate_unclipped_invalid_utf8_mc_tag() {
        use noodles::sam::alignment::record_buf::data::field::value::Value;
        let mc_tag = Tag::from([b'M', b'C']);
        // Create a read then replace the MC tag with invalid UTF-8 bytes
        let mut record = create_mc_test_read("invalid_utf8", 200, "50M");
        record.data_mut().insert(mc_tag, Value::String(vec![0xFF, 0xFE].into()));

        assert_eq!(mate_unclipped_start(&record), None);
        assert_eq!(mate_unclipped_end(&record), None);
    }

    /// A valid CIGAR prefix followed by non-UTF-8 trailing bytes must yield the
    /// same mate coordinates on the typed `RecordBuf` path as on the raw-byte
    /// sibling (which parses the bytes directly). The byte parser stops at the
    /// prefix; the typed path used to reject the whole value on a UTF-8 gate,
    /// disagreeing with the raw path on identical input.
    #[rstest]
    #[case::leading_soft_clip(b"10S50M", 190, 249)]
    #[case::trailing_soft_clip(b"50M10S", 200, 259)]
    fn test_mate_unclipped_valid_prefix_non_utf8_mc_matches_raw(
        #[case] mc_prefix: &[u8],
        #[case] expected_start: isize,
        #[case] expected_end: usize,
    ) {
        let mc_tag = Tag::from([b'M', b'C']);
        // The BAM encoder rejects non-UTF-8 strings, so encode with a valid
        // placeholder trailing byte ('x'), then overwrite it with 0xFF in both
        // the raw buffer and the typed record. 'x' and 0xFF are both non-CIGAR
        // bytes the byte parser skips, so the coordinate is identical either
        // way — the point is a genuinely non-UTF-8 MC value on both paths.
        let mut placeholder = mc_prefix.to_vec();
        placeholder.push(b'x');

        let mut record = create_mc_test_read("valid_prefix_non_utf8", 200, "50M");
        record.data_mut().insert(mc_tag, Value::String(placeholder.clone().into()));
        let mut raw = to_raw(&record);

        // Inject the non-UTF-8 byte into the raw record's MC value.
        let pos = raw
            .windows(placeholder.len())
            .position(|w| w == placeholder.as_slice())
            .expect("placeholder MC value present in raw record");
        raw[pos + placeholder.len() - 1] = 0xFF;

        // And into the typed record.
        let mut non_utf8 = mc_prefix.to_vec();
        non_utf8.push(0xFF);
        record.data_mut().insert(mc_tag, Value::String(non_utf8.into()));

        assert_eq!(mate_unclipped_start(&record), Some(expected_start));
        assert_eq!(mate_unclipped_end(&record), Some(expected_end));
        // Parity with the raw-byte sibling on the identical non-UTF-8 input.
        assert_eq!(mate_unclipped_start_raw(&raw), mate_unclipped_start(&record));
        assert_eq!(mate_unclipped_end_raw(&raw), mate_unclipped_end(&record));
    }

    #[test]
    fn test_mate_unclipped_empty_cigar_mc_tag() {
        let record = create_mc_test_read("empty_cigar", 200, "");
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
        // positive 5' = 100, negative 5' = 100 + 200 - 1 = 299 -> 100 <= 299 -> FR
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
        // RF (outie): positive strand 5' > negative strand 5'
        // Read on + strand at 200, mate on - strand at 100, negative insert size
        // positive 5' = 200, negative 5' = 200 + (-100) - 1 = 99 -> 200 > 99 -> RF
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
        // 100 < 299 -> FR
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
        // 400 >= 299 -> RF
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

    /// A paired, mapped read on reference 0 with an explicit CIGAR, for the 5'-tie tests.
    fn create_tie_test_read(
        flags: u16,
        pos: usize,
        cigar: &str,
        mate_pos: usize,
        tlen: i32,
    ) -> RecordBuf {
        let query_len: usize = cigar
            .split(|c: char| c.is_ascii_alphabetic())
            .zip(cigar.chars().filter(char::is_ascii_alphabetic))
            .filter(|(_, op)| matches!(op, 'M' | 'I' | 'S' | '=' | 'X'))
            .map(|(len, _)| len.parse::<usize>().expect("CIGAR op length"))
            .sum();
        RecordBuilder::new()
            .name("tie")
            .sequence(&"A".repeat(query_len))
            .qualities(&vec![30u8; query_len])
            .cigar(cigar)
            .reference_sequence_id(0)
            .alignment_start(pos)
            .mate_reference_sequence_id(0)
            .mate_alignment_start(mate_pos)
            .template_length(tlen)
            .reverse_complement(flags & FLAG_REVERSE != 0)
            .mate_reverse_complement(flags & FLAG_MATE_REVERSE != 0)
            .first_segment(flags & FLAG_READ1 != 0)
            .build()
    }

    /// htsjdk 5.0.0 `SamPairUtilTest` "dovetail 5' tie" (`testGetPairOrientationDataProvider`):
    /// read 1 forward `100M` @ 100, read 2 reverse `100M` @ 1, so both 5' ends are at 100.
    /// `SamPairUtil.setMateInfo` writes TLEN +1 on read 1 and -1 on read 2 (`computeInsertSize`'s
    /// `+1` adjustment when the 5' ends are equal). htsjdk asserts FR from both ends.
    ///
    /// - forward arm: `100 <= CoordMath.getEnd(100, 1) = 100`, FR;
    /// - reverse arm: `mateStart 100 <= alignmentEnd 100`, FR.
    ///
    /// Each arm fails under a strict `<`, so this pins the inclusive comparison in
    /// `get_pair_orientation` and, through it, `is_fr_pair_from_tags` and `is_fr_pair` in both
    /// argument orders.
    #[test]
    fn test_get_pair_orientation_htsjdk_dovetail_five_prime_tie_is_fr() {
        let fwd =
            create_tie_test_read(FLAG_PAIRED | FLAG_READ1 | FLAG_MATE_REVERSE, 100, "100M", 1, 1);
        let rev = create_tie_test_read(FLAG_PAIRED | FLAG_REVERSE, 1, "100M", 100, -1);

        assert_eq!(get_pair_orientation(&fwd), PairOrientation::FR, "forward arm");
        assert_eq!(get_pair_orientation(&rev), PairOrientation::FR, "reverse arm");
        assert!(is_fr_pair_from_tags(&fwd));
        assert!(is_fr_pair_from_tags(&rev));
        assert!(is_fr_pair(&fwd, &rev));
        assert!(is_fr_pair(&rev, &fwd));
    }

    /// One base past the tie stays RF: reverse `100M` @ 1 ends at 100 and its forward mate
    /// starts at 101 (htsjdk's "nojump outie"; TLEN -2 / +2 from `computeInsertSize`).
    #[test]
    fn test_get_pair_orientation_one_past_five_prime_tie_is_rf() {
        let fwd = create_tie_test_read(FLAG_PAIRED | FLAG_MATE_REVERSE, 101, "100M", 1, -2);
        let rev = create_tie_test_read(FLAG_PAIRED | FLAG_READ1 | FLAG_REVERSE, 1, "100M", 101, 2);

        assert_eq!(get_pair_orientation(&fwd), PairOrientation::RF, "forward arm");
        assert_eq!(get_pair_orientation(&rev), PairOrientation::RF, "reverse arm");
    }

    /// A reverse read with no reference-consuming CIGAR ops (`100S`) at its forward mate's start
    /// is RF. htsjdk `SAMRecord.getAlignmentEnd` is `start + referenceLength - 1`, i.e.
    /// `start - 1` for a zero span, so the reverse arm compares `100 <= 99`. Clamping the end to
    /// `start` would make the inclusive test `100 <= 100` pass and call it FR.
    #[test]
    fn test_get_pair_orientation_zero_reference_span_reverse_read_is_rf() {
        let rev = create_tie_test_read(FLAG_PAIRED | FLAG_REVERSE, 100, "100S", 100, 0);
        assert_eq!(get_pair_orientation(&rev), PairOrientation::RF);
        assert!(!is_fr_pair_from_tags(&rev));
    }

    /// The fgbio#1190 dovetail: forward `100M` at 50 (50-149) and reverse `50S50M` at 1 (aligned
    /// 1-50), so both 5' ends are at 50, with bwa's `TLEN` 0 on both records. Returns
    /// `(forward, reverse)`; `mc` adds each record's `MC` tag.
    fn coincident_dovetail(mc: bool) -> (RecordBuf, RecordBuf) {
        let build = |flags: u16, pos: usize, cigar: &str, mate_pos: usize, mate_cigar: &str| {
            let mut b = RecordBuilder::new()
                .name("dovetail")
                .sequence(&"A".repeat(100))
                .qualities(&[30u8; 100])
                .cigar(cigar)
                .reference_sequence_id(0)
                .alignment_start(pos)
                .mate_reference_sequence_id(0)
                .mate_alignment_start(mate_pos)
                .template_length(0)
                .reverse_complement(flags & FLAG_REVERSE != 0)
                .mate_reverse_complement(flags & FLAG_MATE_REVERSE != 0)
                .first_segment(flags & FLAG_READ1 != 0)
                .paired(flags & FLAG_PAIRED != 0);
            if mc {
                b = b.tag("MC", mate_cigar);
            }
            b.build()
        };
        (
            build(FLAG_PAIRED | FLAG_READ1 | FLAG_MATE_REVERSE, 50, "100M", 1, "50S50M"),
            build(FLAG_PAIRED | FLAG_REVERSE, 1, "50S50M", 50, "100M"),
        )
    }

    /// Per record, htsjdk 5.0.0 takes the forward record's mate end from `MC` when present, so
    /// the coincident dovetail is FR from both records; without `MC` the forward record falls
    /// back to `TLEN` 0 (`50 <= 50 + 0 - 1` is false) and reports RF, while the reverse record's
    /// own alignment end still gives FR.
    #[rstest]
    #[case::with_mc(true, PairOrientation::FR)]
    #[case::without_mc(false, PairOrientation::RF)]
    fn test_get_pair_orientation_forward_record_uses_mc(
        #[case] mc: bool,
        #[case] expected_forward: PairOrientation,
    ) {
        let (fwd, rev) = coincident_dovetail(mc);
        assert_eq!(get_pair_orientation(&fwd), expected_forward, "forward record");
        assert_eq!(is_fr_pair_from_tags(&fwd), expected_forward == PairOrientation::FR);
        assert_eq!(get_pair_orientation(&rev), PairOrientation::FR, "reverse record");
        assert!(is_fr_pair_from_tags(&rev));
    }

    /// Per record with `MC`, the forward arm is the mate's `MC`-derived end, not `TLEN`: a `TLEN`
    /// that alone would say RF does not override an `MC` that says FR.
    #[test]
    fn test_get_pair_orientation_forward_record_prefers_mc_over_tlen() {
        let fwd = create_tie_test_read(FLAG_PAIRED | FLAG_MATE_REVERSE, 200, "100M", 150, -5);
        assert_eq!(get_pair_orientation(&fwd), PairOrientation::RF, "TLEN arm: 200 > 194");
        let mut with_mc = fwd;
        with_mc.data_mut().insert(Tag::from([b'M', b'C']), Value::from("100M"));
        assert_eq!(get_pair_orientation(&with_mc), PairOrientation::FR, "MC arm: 200 <= 249");
    }

    /// An `MC` with a run length no BAM CIGAR op can hold (past 2^28 - 1) is treated as unusable,
    /// so the forward arm falls back to `TLEN` instead of panicking on overflow or wrapping to a
    /// bogus span. Forward
    /// `100M` @ 200 with mate @ 300 and `TLEN` 0: the `TLEN` arm is `200 <= 199`, RF, whereas a
    /// wrapped span would put the mate's end at or just below 300 and call it FR.
    #[rstest]
    #[case::sum_overflows_usize("18446744073709551615M1M")]
    #[case::single_op_exceeds_i64("18446744073709551615M")]
    #[case::span_exceeds_i64_with_start("9223372036854775807M")]
    #[case::run_length_one_past_bam_limit("268435456M")]
    fn test_get_pair_orientation_overflowing_mc_falls_back_to_tlen(#[case] mc: &str) {
        let mut fwd = create_tie_test_read(FLAG_PAIRED | FLAG_MATE_REVERSE, 200, "100M", 300, 0);
        fwd.data_mut().insert(Tag::from([b'M', b'C']), Value::from(mc));
        assert_eq!(get_pair_orientation(&fwd), PairOrientation::RF);
        assert!(!is_fr_pair_from_tags(&fwd));
    }

    /// A malformed `MC` is unusable, so the forward arm falls back to `TLEN` rather than reading
    /// a valid-looking prefix: forward `100M` @ 200, mate @ 150, `TLEN` 0. A valid `100M` puts the
    /// mate's end at 249 (`200 <= 249`, FR); `100M5` (a dangling run length) and other malformed
    /// tags fall back to the `TLEN` arm (`200 <= 199`, RF).
    #[rstest]
    #[case::valid("100M", PairOrientation::FR)]
    #[case::dangling_run_length("100M5", PairOrientation::RF)]
    #[case::unknown_operator("100Mfoo", PairOrientation::RF)]
    #[case::zero_length_op("0M100M", PairOrientation::RF)]
    #[case::no_reference_consuming_op("100S", PairOrientation::RF)]
    fn test_get_pair_orientation_malformed_mc_falls_back_to_tlen(
        #[case] mc: &str,
        #[case] expected: PairOrientation,
    ) {
        let mut fwd = create_tie_test_read(FLAG_PAIRED | FLAG_MATE_REVERSE, 200, "100M", 150, 0);
        fwd.data_mut().insert(Tag::from([b'M', b'C']), Value::from(mc));
        assert_eq!(get_pair_orientation(&fwd), expected);
        assert_eq!(is_fr_pair_from_tags(&fwd), expected == PairOrientation::FR);
    }

    /// `is_fr_pair` is pair-level: the MC-less coincident dovetail is FR whichever record is
    /// passed first, although the forward record's own `TLEN` says RF.
    #[rstest]
    fn test_is_fr_pair_mc_less_coincident_dovetail(#[values(false, true)] reverse_first: bool) {
        let (fwd, rev) = coincident_dovetail(false);
        let fr = if reverse_first { is_fr_pair(&rev, &fwd) } else { is_fr_pair(&fwd, &rev) };
        assert!(fr);
    }

    /// A reverse record whose CIGAR reference span does not fit in `i64` saturates rather than
    /// wrapping: reverse `<usize::MAX>M` @ 1 ends far past its forward mate's start of 100, so the
    /// pair is FR. Casting the span with `as i64` would wrap it to -1, end the reverse read at -1
    /// and call the pair RF.
    ///
    /// A second reference-consuming op whose length pushes the sum past `usize::MAX` saturates
    /// too, rather than panicking on overflow or wrapping the span to zero (which would end the
    /// reverse read at 0 and call the pair RF).
    #[rstest]
    #[case::single_op(&[usize::MAX])]
    #[case::sum_overflows_usize(&[usize::MAX, 1])]
    fn test_pair_orientation_reverse_span_exceeding_i64_saturates(#[case] lens: &[usize]) {
        use noodles::sam::alignment::record::cigar::Op;
        let fwd =
            create_tie_test_read(FLAG_PAIRED | FLAG_READ1 | FLAG_MATE_REVERSE, 100, "100M", 1, 0);
        let mut rev = create_tie_test_read(FLAG_PAIRED | FLAG_REVERSE, 1, "100M", 100, 0);
        *rev.cigar_mut() = lens.iter().map(|&len| Op::new(Kind::Match, len)).collect();

        assert_eq!(get_pair_orientation(&rev), PairOrientation::FR, "reverse arm");
        assert!(is_fr_pair(&fwd, &rev), "(fwd, rev)");
        assert!(is_fr_pair(&rev, &fwd), "(rev, fwd)");
    }

    /// `is_fr_pair` needs both alignment starts: a record with no start is not an FR pair, rather
    /// than being placed at coordinate 0 (`0 <= 199` would call it FR in both orders).
    #[rstest]
    #[case::forward_missing_start(true, false)]
    #[case::reverse_missing_start(false, true)]
    fn test_is_fr_pair_missing_alignment_start(
        #[case] clear_forward: bool,
        #[case] clear_reverse: bool,
    ) {
        let mut fwd =
            create_tie_test_read(FLAG_PAIRED | FLAG_READ1 | FLAG_MATE_REVERSE, 100, "100M", 100, 0);
        let mut rev = create_tie_test_read(FLAG_PAIRED | FLAG_REVERSE, 100, "100M", 100, 0);
        if clear_forward {
            *fwd.alignment_start_mut() = None;
        }
        if clear_reverse {
            *rev.alignment_start_mut() = None;
        }
        assert!(!is_fr_pair(&fwd, &rev), "(fwd, rev)");
        assert!(!is_fr_pair(&rev, &fwd), "(rev, fwd)");
    }

    /// `is_fr_pair` reads the two alignments, not the mate fields: here the forward record's
    /// `MC`, `TLEN` and mate position describe an FR mate, but the reverse record actually ends at
    /// 199, before the forward start of 300, so the pair is RF.
    #[test]
    fn test_is_fr_pair_ignores_stale_mate_fields() {
        let mut fwd = create_tie_test_read(FLAG_MATE_REVERSE | FLAG_READ1, 300, "100M", 500, 300);
        fwd.data_mut().insert(Tag::from([b'M', b'C']), Value::from("100M"));
        let rev = create_tie_test_read(FLAG_REVERSE, 100, "100M", 300, -300);
        let (mut fwd, mut rev) = (fwd, rev);
        *fwd.flags_mut() |= Flags::SEGMENTED;
        *rev.flags_mut() |= Flags::SEGMENTED;
        assert!(is_fr_pair_from_tags(&fwd), "the forward record's own fields say FR");
        assert!(!is_fr_pair(&fwd, &rev));
        assert!(!is_fr_pair(&rev, &fwd));
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
        // negative 5' = 200 + (-50) - 1 = 149, so 200 > 149 -> RF
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

        // Position 105 is within first match block -> read pos 6
        assert_eq!(read_pos_at_ref_pos(&read, 105, false), 6);
        // Position 112 is within skip -> returns 0 (not in deletion)
        assert_eq!(read_pos_at_ref_pos(&read, 112, false), 0);
        // Position 115 is in second match block -> read pos 11
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
        // Simple 50M at position 100 -> unclipped_start = 100
        let read = create_cigar_test_read("simple", 100, "50M");
        assert_eq!(unclipped_start(&read), Some(100));
    }

    #[test]
    fn test_unclipped_start_with_leading_soft_clip() {
        // 5S45M at position 100 -> unclipped_start = 100 - 5 = 95
        let read = create_cigar_test_read("soft", 100, "5S45M");
        assert_eq!(unclipped_start(&read), Some(95));
    }

    #[test]
    fn test_unclipped_start_with_leading_hard_clip() {
        // 10H50M at position 100 -> unclipped_start = 100 - 10 = 90
        let read = create_cigar_test_read("hard", 100, "10H50M");
        assert_eq!(unclipped_start(&read), Some(90));
    }

    #[test]
    fn test_unclipped_start_with_soft_and_hard_clips() {
        // Matches fgbio test: 5S45M10H at position 10 -> unclipped_start = 10 - 5 = 5
        // (trailing hard clip doesn't affect unclipped_start)
        let read = create_cigar_test_read("both", 10, "5S45M10H");
        assert_eq!(unclipped_start(&read), Some(5));
    }

    #[test]
    fn test_unclipped_start_with_leading_hard_and_soft() {
        // 3H5S42M at position 100 -> unclipped_start = 100 - 5 - 3 = 92
        let read = create_cigar_test_read("hard_soft", 100, "3H5S42M");
        assert_eq!(unclipped_start(&read), Some(92));
    }

    #[test]
    fn test_unclipped_end_no_clips() {
        // Simple 50M at position 100 -> alignment_end = 149, unclipped_end = 149
        let read = create_cigar_test_read("simple", 100, "50M");
        assert_eq!(unclipped_end(&read), Some(149));
    }

    #[test]
    fn test_unclipped_end_with_trailing_soft_clip() {
        // 45M5S at position 100 -> alignment_end = 144, unclipped_end = 144 + 5 = 149
        let read = create_cigar_test_read("soft", 100, "45M5S");
        assert_eq!(unclipped_end(&read), Some(149));
    }

    #[test]
    fn test_unclipped_end_with_trailing_hard_clip() {
        // 50M10H at position 100 -> alignment_end = 149, unclipped_end = 149 + 10 = 159
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
        // 42M5S3H at position 100 -> alignment_end = 141, unclipped_end = 141 + 5 + 3 = 149
        let read = create_cigar_test_read("soft_hard", 100, "42M5S3H");
        assert_eq!(unclipped_end(&read), Some(149));
    }

    // R2-CLIP-03: un-*soft*-clipped positions ignore hard clips (fgbio `unSoftClipped*`),
    // unlike `unclipped_*` (htsjdk `getUnclipped*`, which counts hard clips too).
    #[test]
    fn test_unsoftclipped_start_ignores_leading_hard_clip() {
        // 5H10S35M at position 100: unclipped_start = 100 - 5 - 10 = 85 (counts hard);
        // unsoftclipped_start = 100 - 10 = 90 (soft only).
        let read = create_cigar_test_read("hs", 100, "5H10S35M");
        assert_eq!(unclipped_start(&read), Some(85));
        assert_eq!(unsoftclipped_start(&read), Some(90));
    }

    #[test]
    fn test_unsoftclipped_end_ignores_trailing_hard_clip() {
        // 35M10S5H at position 100: alignment_end = 134;
        // unclipped_end = 134 + 10 + 5 = 149 (counts hard);
        // unsoftclipped_end = 134 + 10 = 144 (soft only).
        let read = create_cigar_test_read("sh", 100, "35M10S5H");
        assert_eq!(unclipped_end(&read), Some(149));
        assert_eq!(unsoftclipped_end(&read), Some(144));
    }

    #[test]
    fn test_unsoftclipped_no_hard_matches_unclipped() {
        // With no hard clips, soft-only equals the full unclipped positions.
        let read = create_cigar_test_read("soft", 100, "10S30M10S");
        assert_eq!(unsoftclipped_start(&read), unclipped_start(&read));
        assert_eq!(unsoftclipped_end(&read), unclipped_end(&read));
    }

    #[test]
    fn test_unsoftclipped_unmapped_returns_none() {
        let read = RecordBuilder::new().name("u").sequence("ACGT").flags(Flags::UNMAPPED).build();
        assert_eq!(unsoftclipped_start(&read), None);
        assert_eq!(unsoftclipped_end(&read), None);
    }

    #[test]
    fn test_unsoftclipped_end_empty_cigar_returns_none() {
        // A mapped record with no CIGAR ops has no end position: the typed function must
        // return None, matching the raw sibling `unsoftclipped_end_raw` (not `Some(start)`).
        let read = create_cigar_test_read("empty", 100, "");
        assert!(read.cigar().as_ref().is_empty(), "test setup: CIGAR must be empty");
        assert_eq!(unsoftclipped_end(&read), None);
    }

    #[test]
    fn test_unclipped_five_prime_forward_strand() {
        // Forward strand: 5' is at unclipped_start
        // 5S45M10H at position 10 -> unclipped_start = 5
        let read = create_cigar_test_read("fwd", 10, "5S45M10H");
        assert_eq!(unclipped_five_prime_position(&read), Some(5));
    }

    #[test]
    fn test_unclipped_five_prime_reverse_strand() {
        // Reverse strand: 5' is at unclipped_end
        // 5S45M10H at position 10 -> unclipped_end = 64
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

    // =========================================================================
    // Helpers for raw-sibling tests
    // =========================================================================

    /// Encode a `RecordBuf` to raw BAM bytes for testing raw-byte helpers.
    ///
    /// Builds a header with a single dummy reference sequence (length 100 000) so
    /// that records using `reference_sequence_id(0)` encode without errors.
    fn to_raw(record: &RecordBuf) -> Vec<u8> {
        use noodles::sam::header::record::value::Map;
        use noodles::sam::header::record::value::map::ReferenceSequence;
        use std::num::NonZeroUsize;

        let ref_seq = Map::<ReferenceSequence>::new(
            NonZeroUsize::new(100_000).expect("ref length must be nonzero"),
        );
        let header =
            noodles::sam::Header::builder().add_reference_sequence(b"chr1", ref_seq).build();
        fgumi_raw_bam::encode_record_buf_to_raw(record, &header)
            .expect("encode_record_buf_to_raw should succeed")
            .into_inner()
    }

    /// Convert a typed `(Kind, usize)` CIGAR op pair to the raw u32 BAM packing.
    ///
    /// Encoding: high 28 bits = length, low 4 bits = op type (0=M, 1=I, …, 8=X).
    fn kind_to_raw_op(kind: Kind, len: usize) -> u32 {
        let op_type: u32 = match kind {
            Kind::Match => 0,
            Kind::Insertion => 1,
            Kind::Deletion => 2,
            Kind::Skip => 3,
            Kind::SoftClip => 4,
            Kind::HardClip => 5,
            Kind::Pad => 6,
            Kind::SequenceMatch => 7,
            Kind::SequenceMismatch => 8,
        };
        // CIGAR op lengths are small (genomic context); the cast from usize to u32 is safe.
        let len_u32: u32 = len.try_into().expect("CIGAR op length exceeds u32 range");
        (len_u32 << 4) | op_type
    }

    /// Build a raw `Vec<u32>` CIGAR from a CIGAR string.
    fn raw_ops(cigar_str: &str) -> Vec<u32> {
        parse_cigar_string(cigar_str).into_iter().map(|(k, l)| kind_to_raw_op(k, l)).collect()
    }

    // =========================================================================
    // Tests for leading_clipping_raw / trailing_clipping_raw
    // =========================================================================

    #[test]
    fn test_leading_clipping_raw_matches_typed() {
        // 5H10S50M10S5H — same fixture as test_leading_clipping
        let cigar_str = "5H10S50M10S5H";
        let typed_ops = parse_cigar_string(cigar_str);
        assert_eq!(leading_clipping_raw(&raw_ops(cigar_str)), leading_clipping(&typed_ops));
        assert_eq!(leading_clipping_raw(&raw_ops(cigar_str)), 15); // 5H + 10S
    }

    #[test]
    fn test_trailing_clipping_raw_matches_typed() {
        let cigar_str = "5H10S50M10S5H";
        let typed_ops = parse_cigar_string(cigar_str);
        assert_eq!(trailing_clipping_raw(&raw_ops(cigar_str)), trailing_clipping(&typed_ops));
        assert_eq!(trailing_clipping_raw(&raw_ops(cigar_str)), 15); // 10S + 5H
    }

    #[test]
    fn test_leading_soft_clipping_raw_matches_typed() {
        let cigar_str = "5H10S50M10S5H";
        let typed_ops = parse_cigar_string(cigar_str);
        assert_eq!(
            leading_soft_clipping_raw(&raw_ops(cigar_str)),
            leading_soft_clipping(&typed_ops)
        );
        assert_eq!(leading_soft_clipping_raw(&raw_ops(cigar_str)), 10); // only 10S
    }

    #[test]
    fn test_trailing_soft_clipping_raw_matches_typed() {
        let cigar_str = "5H10S50M10S5H";
        let typed_ops = parse_cigar_string(cigar_str);
        assert_eq!(
            trailing_soft_clipping_raw(&raw_ops(cigar_str)),
            trailing_soft_clipping(&typed_ops)
        );
        assert_eq!(trailing_soft_clipping_raw(&raw_ops(cigar_str)), 10); // only 10S
    }

    #[test]
    fn test_clipping_raw_no_clips() {
        // 50M — no clipping at all
        assert_eq!(leading_clipping_raw(&raw_ops("50M")), 0);
        assert_eq!(trailing_clipping_raw(&raw_ops("50M")), 0);
        assert_eq!(leading_soft_clipping_raw(&raw_ops("50M")), 0);
        assert_eq!(trailing_soft_clipping_raw(&raw_ops("50M")), 0);
    }

    #[test]
    fn test_clipping_raw_empty() {
        let empty: Vec<u32> = vec![];
        assert_eq!(leading_clipping_raw(&empty), 0);
        assert_eq!(trailing_clipping_raw(&empty), 0);
        assert_eq!(leading_soft_clipping_raw(&empty), 0);
        assert_eq!(trailing_soft_clipping_raw(&empty), 0);
    }

    // =========================================================================
    // Tests for cigar_reference_length_raw / cigar_query_length_raw
    // =========================================================================

    #[test]
    fn test_cigar_reference_length_raw() {
        // 10M5I20M5D10M — ref: 10+20+5+10=45
        let cigar_str = "10M5I20M5D10M";
        let typed_ops = parse_cigar_string(cigar_str);
        assert_eq!(
            cigar_reference_length_raw(&raw_ops(cigar_str)),
            cigar_reference_length(&typed_ops)
        );
        assert_eq!(cigar_reference_length_raw(&raw_ops(cigar_str)), 45);
    }

    #[test]
    fn test_cigar_query_length_raw() {
        // 10M5I20M5D10M — query-consuming (M+I+M+M): 10+5+20+10 = 45
        let cigar_str = "10M5I20M5D10M";
        assert_eq!(cigar_query_length_raw(&raw_ops(cigar_str)), 45);
    }

    // =========================================================================
    // Tests for unclipped_start_raw / unclipped_end_raw
    // =========================================================================

    #[test]
    fn test_unclipped_start_raw_no_clips() {
        let record = create_cigar_test_read("nosc", 100, "50M");
        let raw = to_raw(&record);
        assert_eq!(unclipped_start_raw(&raw), unclipped_start(&record));
        assert_eq!(unclipped_start_raw(&raw), Some(100));
    }

    #[test]
    fn test_unclipped_start_raw_with_soft_clip() {
        let record = create_cigar_test_read("sc", 100, "10S40M");
        let raw = to_raw(&record);
        assert_eq!(unclipped_start_raw(&raw), unclipped_start(&record));
        assert_eq!(unclipped_start_raw(&raw), Some(90)); // 100 - 10
    }

    #[test]
    fn test_unclipped_start_raw_with_hard_and_soft() {
        let record = create_cigar_test_read("hsc", 100, "5H10S35M");
        let raw = to_raw(&record);
        assert_eq!(unclipped_start_raw(&raw), unclipped_start(&record));
        assert_eq!(unclipped_start_raw(&raw), Some(85)); // 100 - 5H - 10S
    }

    #[test]
    fn test_unclipped_start_raw_unmapped_returns_none() {
        let record = RecordBuilder::new().name("u").sequence("ACGT").flags(Flags::UNMAPPED).build();
        let raw = to_raw(&record);
        assert_eq!(unclipped_start_raw(&raw), None);
    }

    #[test]
    fn test_unclipped_end_raw_no_clips() {
        let record = create_cigar_test_read("nosc", 100, "50M");
        let raw = to_raw(&record);
        assert_eq!(unclipped_end_raw(&raw), unclipped_end(&record));
        assert_eq!(unclipped_end_raw(&raw), Some(149)); // 100 + 50 - 1
    }

    #[test]
    fn test_unclipped_end_raw_with_soft_clip() {
        let record = create_cigar_test_read("sc", 100, "40M10S");
        let raw = to_raw(&record);
        assert_eq!(unclipped_end_raw(&raw), unclipped_end(&record));
        assert_eq!(unclipped_end_raw(&raw), Some(149)); // 100+40-1+10 = 149
    }

    #[test]
    fn test_unclipped_end_raw_with_hard_clip() {
        let record = create_cigar_test_read("hc", 100, "45M5H");
        let raw = to_raw(&record);
        assert_eq!(unclipped_end_raw(&raw), unclipped_end(&record));
        assert_eq!(unclipped_end_raw(&raw), Some(149)); // 100+45-1+5
    }

    #[test]
    fn test_unclipped_end_raw_unmapped_returns_none() {
        let record = RecordBuilder::new().name("u").sequence("ACGT").flags(Flags::UNMAPPED).build();
        let raw = to_raw(&record);
        assert_eq!(unclipped_end_raw(&raw), None);
    }

    // R2-CLIP-03 (raw): soft-only unclipped positions ignore hard clips.
    #[test]
    fn test_unsoftclipped_start_raw_ignores_hard_clip() {
        let record = create_cigar_test_read("hs", 100, "5H10S35M");
        let raw = to_raw(&record);
        assert_eq!(unsoftclipped_start_raw(&raw), unsoftclipped_start(&record));
        assert_eq!(unsoftclipped_start_raw(&raw), Some(90)); // 100 - 10S (hard ignored)
        assert_eq!(unclipped_start_raw(&raw), Some(85)); // 100 - 5H - 10S
    }

    #[test]
    fn test_unsoftclipped_end_raw_ignores_hard_clip() {
        let record = create_cigar_test_read("sh", 100, "35M10S5H");
        let raw = to_raw(&record);
        assert_eq!(unsoftclipped_end_raw(&raw), unsoftclipped_end(&record));
        assert_eq!(unsoftclipped_end_raw(&raw), Some(144)); // 134 + 10S (hard ignored)
        assert_eq!(unclipped_end_raw(&raw), Some(149)); // 134 + 10S + 5H
    }

    #[test]
    fn test_unsoftclipped_raw_unmapped_returns_none() {
        let record = RecordBuilder::new().name("u").sequence("ACGT").flags(Flags::UNMAPPED).build();
        let raw = to_raw(&record);
        assert_eq!(unsoftclipped_start_raw(&raw), None);
        assert_eq!(unsoftclipped_end_raw(&raw), None);
    }

    // =========================================================================
    // Tests for mate_unclipped_start_raw / mate_unclipped_end_raw
    // =========================================================================

    #[test]
    fn test_mate_unclipped_start_raw_no_clip() {
        let record = create_mc_test_read("mc", 200, "50M");
        let raw = to_raw(&record);
        assert_eq!(mate_unclipped_start_raw(&raw), mate_unclipped_start(&record));
        assert_eq!(mate_unclipped_start_raw(&raw), Some(200));
    }

    #[test]
    fn test_mate_unclipped_start_raw_with_soft_clip() {
        // 10S50M: leading clip = 10, unclipped start = 200 - 10 = 190
        let record = create_mc_test_read("mc_sc", 200, "10S50M");
        let raw = to_raw(&record);
        assert_eq!(mate_unclipped_start_raw(&raw), mate_unclipped_start(&record));
        assert_eq!(mate_unclipped_start_raw(&raw), Some(190));
    }

    #[test]
    fn test_mate_unclipped_end_raw_no_clip() {
        // 50M: ref_len=50, unclipped_end = 200 + 50 - 1 = 249
        let record = create_mc_test_read("mc", 200, "50M");
        let raw = to_raw(&record);
        assert_eq!(mate_unclipped_end_raw(&raw), mate_unclipped_end(&record));
        assert_eq!(mate_unclipped_end_raw(&raw), Some(249));
    }

    #[test]
    fn test_mate_unclipped_end_raw_with_soft_clip() {
        // 50M10S: ref_len=50, trailing=10, unclipped_end = 200 + 50 - 1 + 10 = 259
        let record = create_mc_test_read("mc_sc", 200, "50M10S");
        let raw = to_raw(&record);
        assert_eq!(mate_unclipped_end_raw(&raw), mate_unclipped_end(&record));
        assert_eq!(mate_unclipped_end_raw(&raw), Some(259));
    }

    #[test]
    fn test_mate_unclipped_start_raw_no_mc_tag_returns_none() {
        let record = RecordBuilder::new()
            .name("no_mc")
            .sequence(&"A".repeat(50))
            .qualities(&[30u8; 50])
            .cigar("50M")
            .reference_sequence_id(0)
            .alignment_start(100)
            .build();
        let raw = to_raw(&record);
        assert_eq!(mate_unclipped_start_raw(&raw), None);
        assert_eq!(mate_unclipped_end_raw(&raw), None);
    }

    /// Raw paths must reject empty MC CIGARs to match the typed path's `None` contract.
    #[test]
    fn test_mate_unclipped_raw_empty_mc_cigar_returns_none() {
        let record = create_mc_test_read("empty_mc", 200, "");
        let raw = to_raw(&record);
        assert_eq!(mate_unclipped_start_raw(&raw), None);
        assert_eq!(mate_unclipped_end_raw(&raw), None);
        // typed path behaves the same way
        assert_eq!(mate_unclipped_start(&record), None);
        assert_eq!(mate_unclipped_end(&record), None);
    }
}

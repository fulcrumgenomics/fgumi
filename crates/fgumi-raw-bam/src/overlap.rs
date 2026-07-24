use crate::SamTag;
use crate::cigar::{
    consumes_query, get_cigar_ops, reference_length_from_cigar, reference_length_from_raw_bam,
};
use crate::fields::{RawRecordView, aux_data_slice, flags, mate_pos, mate_ref_id, template_length};
use crate::tags::RawTagsView;

/// Check if a single read is part of an FR (forward-reverse) pair using raw BAM bytes.
///
/// This is the raw-byte equivalent of `record_utils::is_fr_pair_from_tags`.
/// Returns `true` if the read is paired, both read and mate are mapped,
/// on the same reference, and in FR orientation (positive strand 5' < negative strand 5').
#[must_use]
pub fn is_fr_pair_raw(bam: &[u8]) -> bool {
    let v = RawRecordView::new(bam);
    let flg = v.flags();

    // Must be paired
    if flg & flags::PAIRED == 0 {
        return false;
    }

    // Both read and mate must be mapped
    if flg & flags::UNMAPPED != 0 || flg & flags::MATE_UNMAPPED != 0 {
        return false;
    }

    // Must be on the same reference
    let this_ref_id = v.ref_id();
    let m_ref_id = mate_ref_id(bam);
    if this_ref_id != m_ref_id {
        return false;
    }

    // Must be on opposite strands for FR or RF
    let is_reverse = flg & flags::REVERSE != 0;
    let mate_is_reverse = flg & flags::MATE_REVERSE != 0;
    if is_reverse == mate_is_reverse {
        return false;
    }

    // Determine if FR or RF using htsjdk's logic:
    // positiveStrandFivePrimePos = readIsOnReverseStrand ? mateStart : alignmentStart
    // negativeStrandFivePrimePos = readIsOnReverseStrand ? alignmentEnd : alignmentStart + insertSize
    let alignment_start = v.pos() + 1; // 1-based
    let m_start = mate_pos(bam) + 1; // 1-based
    let insert_size = template_length(bam);

    let (positive_five_prime, negative_five_prime) = if is_reverse {
        // This read is on reverse strand, mate is on positive strand
        let ref_len = reference_length_from_raw_bam(bam);
        let end = alignment_start + (ref_len - 1).max(0);
        (m_start, end)
    } else {
        // This read is on positive strand, mate is on reverse strand
        (alignment_start, alignment_start + insert_size)
    };

    // FR if positive strand 5' < negative strand 5'
    positive_five_prime < negative_five_prime
}

/// Symmetric per-pair FR classification for a template's two primary reads.
///
/// This is the raw-byte port of fgbio `CodecConsensusCaller.isPrimaryFrPair(a, b)`
/// (`CodecConsensusCaller.scala:424-433`). Unlike [`is_fr_pair_raw`], which is a *per-record*
/// test, this evaluates a *pair* and is order-independent: after confirming both reads are
/// mapped with mapped mates, on the same reference, and on opposite strands, it derives the FR
/// orientation from the **reverse-strand record only**. That branch of the orientation test is
/// CIGAR-derived (`is_fr_pair_raw`'s reverse arm), so it is independent of TLEN and gives the
/// same answer regardless of argument order — avoiding the htsjdk per-record asymmetry
/// (samtools/htsjdk#1771) that mis-drops dovetail pairs whose aligned ends coincide.
///
/// Returns `true` iff `(a, b)` is a single primary FR pair.
#[must_use]
pub fn is_primary_fr_pair_raw(a: &[u8], b: &[u8]) -> bool {
    let (va, vb) = (RawRecordView::new(a), RawRecordView::new(b));
    let (fa, fb) = (va.flags(), vb.flags());

    // Both reads mapped.
    if fa & flags::UNMAPPED != 0 || fb & flags::UNMAPPED != 0 {
        return false;
    }
    // Both mates mapped.
    if fa & flags::MATE_UNMAPPED != 0 || fb & flags::MATE_UNMAPPED != 0 {
        return false;
    }
    // Same reference.
    if va.ref_id() != vb.ref_id() {
        return false;
    }
    // Opposite strands.
    let a_reverse = fa & flags::REVERSE != 0;
    let b_reverse = fb & flags::REVERSE != 0;
    if a_reverse == b_reverse {
        return false;
    }
    // Evaluate FR orientation on the reverse-strand record only (the symmetric branch).
    let reverse_record = if a_reverse { a } else { b };
    is_fr_pair_raw(reverse_record)
}

/// Number of bases a read extends past its mate for FR pairs, deriving the mate's boundary
/// from the read's own **MC tag** (the mate CIGAR). Returns `0` for non-FR pairs or when the
/// MC tag is absent/invalid.
///
/// The mate boundary uses **soft-only** unclipped coordinates (fgbio `mateUnSoftClippedStart`/
/// `mateUnSoftClippedEnd`) — hard clips do NOT extend it. This deliberately differs from the
/// samtools-style soft+hard unclipped 5' position used for template-coordinate sorting/grouping
/// (see [`crate::cigar::mate_unclipped_5prime`]): overlap clipping counts only bases physically
/// present in the mate's sequence, whereas the sort/group key follows samtools' soft+hard
/// convention. Used by the simplex (vanilla) and duplex consensus callers, mirroring fgbio's
/// base `UmiConsensusCaller.numBasesExtendingPastMate`. The codec caller instead uses
/// [`num_bases_extending_past_mate_vs_mate_raw`], which reads the mate boundary from the mate
/// record in hand rather than its MC tag (see CODEC3-04).
#[must_use]
pub fn num_bases_extending_past_mate_raw(bam: &[u8]) -> usize {
    // Only applies to FR pairs (matching the RecordBuf-based `num_bases_extending_past_mate`)
    if !is_fr_pair_raw(bam) {
        return 0;
    }

    // Need MC tag for mate CIGAR information
    let aux = aux_data_slice(bam);
    let Some(mc_bytes) = RawTagsView::new(aux).find_string(SamTag::MC) else {
        return 0;
    };
    let Ok(mc_cigar) = std::str::from_utf8(mc_bytes) else {
        return 0;
    };

    // A malformed MC CIGAR fails closed to 0 (no clip), honoring the doc contract above
    // ("Returns 0 ... when the MC tag is absent/invalid").
    let Some((mate_unclipped_start, mate_unclipped_end)) =
        mate_soft_unclipped_from_mc(mate_pos(bam) + 1, mc_cigar)
    else {
        return 0;
    };
    bases_extending_past_mate(bam, mate_unclipped_start, mate_unclipped_end)
}

/// Number of bases a read extends past its mate for FR pairs, deriving the mate's boundary
/// from the **mate record in hand** rather than the read's MC tag.
///
/// Codec (and any caller holding both primary reads) uses this so overlap clipping does not
/// silently no-op when the MC tag is missing — mirroring fgbio's `updateMateCigars`, which
/// backfills the mate CIGAR from the in-group mate before clipping. Uses the same soft-only
/// mate boundary as [`num_bases_extending_past_mate_raw`].
///
/// Because both records are in hand, FR classification uses the symmetric per-pair
/// [`is_primary_fr_pair_raw`] rather than the per-record [`is_fr_pair_raw`]. The per-record test is
/// asymmetric for dovetail pairs (htsjdk/samtools#1771): a valid dovetail-forward pair can be
/// mis-classified as non-FR from the forward read's perspective, which would zero its clip even
/// though the caller (e.g. the codec caller's `is_primary_fr_pair_raw` gate) accepted the pair.
#[must_use]
pub fn num_bases_extending_past_mate_vs_mate_raw(rec: &[u8], mate: &[u8]) -> usize {
    if !is_primary_fr_pair_raw(rec, mate) {
        return 0;
    }
    let mate_pos_1based = RawRecordView::new(mate).pos() + 1;
    let mate_ops = get_cigar_ops(mate);
    let (mate_unclipped_start, mate_unclipped_end) =
        mate_soft_unclipped_from_ops(mate_pos_1based, &mate_ops);
    bases_extending_past_mate(rec, mate_unclipped_start, mate_unclipped_end)
}

/// Core of the overlap-clip computation shared by the MC-based and mate-record-based entry
/// points: given the mate's soft-only unclipped start/end (1-based), returns how many bases at
/// the read's 3' (sequencing) end extend past the mate. Mirrors fgbio
/// `SamRecordClipper.numBasesExtendingPastMate`.
#[must_use]
fn bases_extending_past_mate(
    bam: &[u8],
    mate_unclipped_start: i32,
    mate_unclipped_end: i32,
) -> usize {
    let v = RawRecordView::new(bam);
    let is_reverse = v.flags() & flags::REVERSE != 0;
    let this_pos = v.pos() + 1; // 1-based

    let cigar_ops = get_cigar_ops(bam);
    let read_length: usize = cigar_ops
        .iter()
        .map(|&op| {
            let op_type = op & 0xF;
            let op_len = (op >> 4) as usize;
            if consumes_query(op_type) { op_len } else { 0 }
        })
        .sum();

    if is_reverse {
        // Negative strand: check if read extends before mate's unclipped start
        if this_pos <= mate_unclipped_start {
            compute_bases_before_ref_pos(&cigar_ops, this_pos, mate_unclipped_start)
        } else {
            // Only clip excess soft-clipped bases at the start.
            // saturating: an adversarial MC-tag/mate CIGAR can push mate_unclipped_start far
            // negative (leading_soft saturates to i32::MAX upstream), so a plain subtraction
            // would overflow i32 here — panic in debug, wrap to garbage in release.
            let leading_sc = leading_soft_clip_from_ops(&cigar_ops);
            let gap = this_pos.saturating_sub(mate_unclipped_start).cast_unsigned() as usize;
            leading_sc.saturating_sub(gap)
        }
    } else {
        // Positive strand: check if read extends past mate's unclipped end.
        // saturating: an oversized decoded CIGAR can drive `ref_len` up to `i32::MAX`, so a plain
        // `this_pos + ref_len - 1` would overflow i32 before clipping can fail closed — panic in
        // debug, wrap to garbage in release.
        let ref_len = reference_length_from_cigar(&cigar_ops);
        // `-1` first, matching the mate-boundary helpers: subtracting it after the
        // saturating add clamps a genuinely-at-the-ceiling end inward to `i32::MAX - 1`.
        let alignment_end = this_pos.saturating_sub(1).saturating_add(ref_len);

        if alignment_end >= mate_unclipped_end {
            // Compute read position at mate's unclipped end
            let bases_past = compute_bases_past_ref_pos(&cigar_ops, this_pos, mate_unclipped_end);
            read_length.saturating_sub(bases_past)
        } else {
            // Only clip excess soft-clipped bases.
            // saturating for symmetry with the reverse-strand branch: mate_unclipped_end only
            // saturates upward (bounded by i32::MAX) so this is not currently reachable-overflow,
            // but saturating the subtraction closes the analogous gap defensively and cheaply.
            let trailing_sc = trailing_soft_clip_from_ops(&cigar_ops);
            let gap = mate_unclipped_end.saturating_sub(alignment_end).cast_unsigned() as usize;
            trailing_sc.saturating_sub(gap)
        }
    }
}

/// Mate's soft-only unclipped start/end (1-based) from an MC-tag CIGAR **string** and the
/// mate's 1-based leftmost position. Only soft clips count — hard clips are excluded — matching
/// fgbio `mateUnSoftClippedStart`/`mateUnSoftClippedEnd`.
fn mate_soft_unclipped_from_mc(mate_pos_1based: i32, mc_cigar: &str) -> Option<(i32, i32)> {
    let (leading_soft, ref_len, trailing_soft) = parse_soft_clips_and_ref_len(mc_cigar)?;
    // saturating: a saturated ref_len/soft-clip (from an adversarial MC CIGAR) must
    // clamp the boundary, not overflow-panic (debug) / wrap to a negative (release).
    // The inclusive-end `-1` is applied *first*, for the reason given on
    // `cigar::unclipped_other_end`: applied last it would saturate the sum to
    // `i32::MAX` and then walk it back to `i32::MAX - 1`, clamping inward rather
    // than outward -- the opposite of what the fail-closed contract wants. As a
    // plain `- 1` it also underflow-panicked in debug when `mate_pos_1based` was
    // `i32::MIN` and the spans were zero.
    Some((
        mate_pos_1based.saturating_sub(leading_soft),
        mate_pos_1based.saturating_sub(1).saturating_add(ref_len).saturating_add(trailing_soft),
    ))
}

/// Mate's soft-only unclipped start/end (1-based) from the mate record's own CIGAR **ops**.
fn mate_soft_unclipped_from_ops(mate_pos_1based: i32, mate_ops: &[u32]) -> (i32, i32) {
    // Both soft-clip totals clamp to i32::MAX (not 0) when they overflow i32: an oversized
    // trailing clip that collapsed to 0 would shrink mate_unclipped_end, making the read look
    // like it extends further past the mate and over-clipping the consensus. Fail closed by
    // clamping the boundary outward, matching leading_soft.
    let leading_soft = i32::try_from(leading_soft_clip_from_ops(mate_ops)).unwrap_or(i32::MAX);
    let trailing_soft = i32::try_from(trailing_soft_clip_from_ops(mate_ops)).unwrap_or(i32::MAX);
    let ref_len = reference_length_from_cigar(mate_ops);
    // The inclusive-end `-1` is applied first; see `mate_soft_unclipped_from_mc`.
    (
        mate_pos_1based.saturating_sub(leading_soft),
        mate_pos_1based.saturating_sub(1).saturating_add(ref_len).saturating_add(trailing_soft),
    )
}

/// Parses `(leading_soft_clips, reference_length, trailing_soft_clips)` from a CIGAR string,
/// counting only soft (`S`) clips — hard clips (`H`) are ignored, giving the soft-only unclipped
/// extent fgbio uses for overlap clipping.
///
/// Returns `None` for any CIGAR that is not a structurally valid SAM CIGAR for a mapped record,
/// so malformed MC-tag input fails closed rather than yielding a partial boundary. Rejected:
/// - an unknown operator byte (e.g. the `f`/`o` in `10Mfoo5S`, or a non-ASCII byte);
/// - an operator with no preceding run-length (a bare `MS`) or a zero run-length (`0M`, `10M0S`);
/// - a trailing run-length with no operator (`10M5`), or an empty CIGAR;
/// - a soft (`S`) or hard (`H`) clip that is not at the end of the CIGAR — `S` may only sit at the
///   ends (inside any `H`), and `H` only as the first/last operation (e.g. `10M5S10M` is invalid);
/// - a CIGAR with no reference-consuming operation (`10S`), which a mapped mate cannot have.
fn parse_soft_clips_and_ref_len(cigar: &str) -> Option<(i32, i32, i32)> {
    // Phase 1: tokenize into (run_length, operator) pairs. Accumulate the run-length from digit
    // bytes directly rather than slicing the string by byte offset — slicing on a non-ASCII byte
    // (a malformed CIGAR) would panic on a char boundary. Reject lexical errors here.
    let mut tokens: Vec<(i32, u8)> = Vec::new();
    let mut num = 0i32;
    let mut have_digits = false;
    for &c in cigar.as_bytes() {
        if c.is_ascii_digit() {
            // saturating: an adversarial MC CIGAR can sum a run-length past i32::MAX; clamp
            // rather than overflow-panic (debug) / wrap (release).
            num = num.saturating_mul(10).saturating_add(i32::from(c - b'0'));
            have_digits = true;
            continue;
        }
        // Every operator needs a positive run-length, and must be a known CIGAR operator.
        if !have_digits || num == 0 {
            return None;
        }
        if !matches!(c, b'M' | b'I' | b'D' | b'N' | b'S' | b'H' | b'P' | b'=' | b'X') {
            return None;
        }
        tokens.push((num, c));
        num = 0;
        have_digits = false;
    }
    // A dangling run-length with no operator (`10M5`), or an empty CIGAR, is malformed.
    if have_digits || tokens.is_empty() {
        return None;
    }

    // Phase 2: validate structure and accumulate the soft-only boundary in one pass. Soft clips
    // must sit at the ends (inside any hard clip); hard clips only as the first/last operation.
    let last = tokens.len() - 1;
    let mut leading_soft = 0i32;
    let mut trailing_soft = 0i32;
    let mut ref_len = 0i32;
    let mut saw_ref_op = false;
    for (i, &(len, op)) in tokens.iter().enumerate() {
        match op {
            b'M' | b'D' | b'N' | b'=' | b'X' => {
                ref_len = ref_len.saturating_add(len);
                saw_ref_op = true;
            }
            b'I' | b'P' => {} // interior op: no boundary contribution, no placement constraint
            b'S' => {
                let leading = tokens[..i].iter().all(|&(_, o)| o == b'H');
                let trailing = tokens[i + 1..].iter().all(|&(_, o)| o == b'H');
                if !leading && !trailing {
                    return None; // internal soft clip: invalid placement
                }
                if saw_ref_op {
                    trailing_soft = trailing_soft.saturating_add(len);
                } else {
                    leading_soft = leading_soft.saturating_add(len);
                }
            }
            b'H' if i == 0 || i == last => {} // hard clip only permitted at the ends
            // An internal hard clip is invalid placement; the tokenizer already rejected unknown
            // operators, so every other operator byte here also fails closed.
            _ => return None,
        }
    }
    // A mapped mate's CIGAR must consume reference; a fully-clipped CIGAR (`10S`) cannot.
    if !saw_ref_op {
        return None;
    }

    Some((leading_soft, ref_len, trailing_soft))
}

/// Whether to return the read position at or past the target reference position,
/// or the read position before it.
#[derive(Clone, Copy)]
enum RefPosMode {
    /// Return the 1-based read position at the target (for positive strand clipping).
    AtOrPast,
    /// Return the 1-based read position before the target (for negative strand clipping).
    Before,
}

/// Compute the read position relative to a reference position by walking the CIGAR.
///
/// Returns the 1-based read position at (or just before) the given reference position,
/// or 0 if the position falls in a deletion/skip or outside the alignment.
fn compute_read_pos_at_ref(
    cigar_ops: &[u32],
    alignment_start_1based: i32,
    target_ref_pos: i32,
    mode: RefPosMode,
) -> usize {
    let mut ref_pos = alignment_start_1based;
    let mut read_pos: usize = 0;

    for &op in cigar_ops {
        let op_type = op & 0xF;
        let op_len = (op >> 4) as usize;

        match op_type {
            0 | 7 | 8 => {
                // M, =, X: consume both query and reference
                for _ in 0..op_len {
                    read_pos += 1;
                    if ref_pos == target_ref_pos {
                        return match mode {
                            RefPosMode::AtOrPast => read_pos,
                            RefPosMode::Before => read_pos.saturating_sub(1),
                        };
                    }
                    ref_pos += 1;
                }
            }
            1 => {
                // I: consume query only
                read_pos += op_len;
            }
            4 => {
                // S: consume query only
                read_pos += op_len;
            }
            2 | 3 => {
                // D, N: consume reference only
                for _ in 0..op_len {
                    if ref_pos == target_ref_pos {
                        return 0; // Position in deletion
                    }
                    ref_pos += 1;
                }
            }
            _ => {}
        }
    }

    0
}

/// Compute number of read bases at or past a reference position (for positive strand).
///
/// Returns the 1-based read position at the given reference position,
/// or 0 if the position falls in a deletion or outside the alignment.
fn compute_bases_past_ref_pos(
    cigar_ops: &[u32],
    alignment_start_1based: i32,
    target_ref_pos: i32,
) -> usize {
    compute_read_pos_at_ref(cigar_ops, alignment_start_1based, target_ref_pos, RefPosMode::AtOrPast)
}

/// Compute number of read bases before a reference position (for negative strand).
fn compute_bases_before_ref_pos(
    cigar_ops: &[u32],
    alignment_start_1based: i32,
    target_ref_pos: i32,
) -> usize {
    compute_read_pos_at_ref(cigar_ops, alignment_start_1based, target_ref_pos, RefPosMode::Before)
}

/// Count trailing soft clips from CIGAR ops.
fn trailing_soft_clip_from_ops(cigar_ops: &[u32]) -> usize {
    let mut trailing = 0usize;
    for &op in cigar_ops.iter().rev() {
        let op_type = op & 0xF;
        let op_len = (op >> 4) as usize;
        match op_type {
            4 => trailing += op_len, // S
            5 => {}                  // H - skip
            _ => break,
        }
    }
    trailing
}

/// Count leading soft clips from CIGAR ops.
fn leading_soft_clip_from_ops(cigar_ops: &[u32]) -> usize {
    let mut leading = 0usize;
    for &op in cigar_ops {
        let op_type = op & 0xF;
        let op_len = (op >> 4) as usize;
        match op_type {
            4 => leading += op_len, // S
            5 => {}                 // H - skip
            _ => break,
        }
    }
    leading
}

#[cfg(test)]
#[allow(clippy::identity_op)]
mod tests {
    use super::*;
    use crate::testutil::*;
    use rstest::rstest;

    // ========================================================================
    // parse_soft_clips_and_ref_len tests
    // ========================================================================

    /// `(leading_soft, ref_len, trailing_soft)`, counting only soft clips (H excluded).
    /// Any CIGAR that is not a structurally valid SAM CIGAR for a mapped record fails closed
    /// to `None` rather than yielding a partial boundary.
    #[rstest]
    // --- valid CIGARs ---
    #[case::simple_match("40M", Some((0, 40, 0)))]
    #[case::leading_soft("5S10M", Some((5, 10, 0)))]
    #[case::trailing_soft("10M3S", Some((0, 10, 3)))]
    #[case::both_soft("5S10M3S", Some((5, 10, 3)))]
    #[case::hard_excluded("5H10M2H", Some((0, 10, 0)))]
    #[case::soft_and_hard("2H5S10M3S2H", Some((5, 10, 3)))]
    #[case::ref_ops_summed("5M2D3N4M", Some((0, 14, 0)))]
    #[case::insertion_ignored("5M2I3M", Some((0, 8, 0)))]
    // A malformed/adversarial MC CIGAR whose ref ops sum past i32::MAX must
    // saturate rather than overflow the `ref_len += num` accumulation.
    #[case::oversized_ref_ops_saturate("2000000000M2000000000M", Some((0, i32::MAX, 0)))]
    // --- lexically malformed: unknown byte, missing/zero run-length, dangling length, empty ---
    #[case::non_ascii_rejected("10M\u{20ac}5S", None)]
    #[case::embedded_garbage_rejected("10Mfoo5S", None)]
    #[case::bare_operator_rejected("MS", None)]
    #[case::dangling_length_rejected("10M5", None)]
    #[case::empty_rejected("", None)]
    #[case::zero_length_ref_rejected("0M", None)]
    #[case::zero_length_soft_rejected("10M0S", None)]
    // --- structurally malformed: misplaced clips, or no reference-consuming operation ---
    #[case::internal_soft_rejected("10M5S10M", None)]
    #[case::internal_hard_rejected("10M5H10M", None)]
    #[case::leading_hard_after_soft_rejected("5S2H10M", None)]
    #[case::no_ref_op_soft_only_rejected("36S", None)]
    #[case::no_ref_op_insertion_only_rejected("5S5I", None)]
    fn test_parse_soft_clips_and_ref_len(
        #[case] cigar: &str,
        #[case] expected: Option<(i32, i32, i32)>,
    ) {
        assert_eq!(parse_soft_clips_and_ref_len(cigar), expected);
    }

    /// A crafted (valid-UTF-8) MC-tag CIGAR with an op length that saturates
    /// `ref_len` must not overflow the boundary arithmetic
    /// (`mate_pos_1based - 1 + ref_len + trailing_soft`). Pre-fix this panicked in
    /// debug/test and wrapped to a bogus negative boundary in release, corrupting
    /// the overlap-clip amount. A well-formed BAM cannot reach this (real ref
    /// lengths are chromosome-bounded), but the MC Z-string is never length-checked.
    #[rstest]
    #[case::from_mc(mate_soft_unclipped_from_mc(100, "9999999999M"))]
    fn mate_soft_unclipped_saturates_on_oversized_cigar(#[case] bounds: Option<(i32, i32)>) {
        let (start, end) = bounds.expect("well-formed (if oversized) CIGAR parses");
        assert_eq!(start, 100, "no leading soft clip -> start stays at mate_pos");
        assert!(end >= start, "end boundary must stay >= start, not wrap negative");
        assert_eq!(end, i32::MAX, "saturated ref_len clamps the end boundary");
    }

    /// The inclusive-end `-1` must be applied before the saturating additions, in
    /// both mate-boundary helpers.
    ///
    /// Applied last it saturates the sum to `i32::MAX` and then walks it back, so a
    /// boundary whose exact value *is* `i32::MAX` was reported as `i32::MAX - 1` --
    /// clamping the fail-closed boundary inward by one base rather than outward.
    /// This is the `overlap.rs` sibling of `cigar::unclipped_other_end`.
    #[rstest]
    // Exact ceiling, no saturation warranted: MAX + 1M - 1 == MAX.
    #[case::from_mc_exactly_max(mate_soft_unclipped_from_mc(i32::MAX, "1M"), i32::MAX)]
    // Genuinely past the ceiling: clamping outward to MAX is correct.
    #[case::from_mc_past_max(mate_soft_unclipped_from_mc(i32::MAX, "10M"), i32::MAX)]
    // The ordinary case must be untouched by the reordering: 100 + 10 - 1 == 109.
    #[case::from_mc_ordinary(mate_soft_unclipped_from_mc(100, "10M"), 109)]
    // A trailing soft clip still extends the boundary: 100 + 10 + 5 - 1 == 114.
    #[case::from_mc_trailing_soft(mate_soft_unclipped_from_mc(100, "10M5S"), 114)]
    fn mate_soft_unclipped_end_applies_the_inclusive_minus_one_first(
        #[case] bounds: Option<(i32, i32)>,
        #[case] expected_end: i32,
    ) {
        let (_, end) = bounds.expect("well-formed CIGAR parses");
        assert_eq!(end, expected_end);
    }

    /// Same contract for the mate-record-based helper, which shares the formula.
    #[rstest]
    #[case::exactly_max(i32::MAX, &[encode_op(0, 1)], i32::MAX)]
    #[case::past_max(i32::MAX, &[encode_op(0, 10)], i32::MAX)]
    #[case::ordinary(100, &[encode_op(0, 10)], 109)]
    #[case::trailing_soft(100, &[encode_op(0, 10), encode_op(4, 5)], 114)]
    fn mate_soft_unclipped_from_ops_applies_the_inclusive_minus_one_first(
        #[case] mate_pos_1based: i32,
        #[case] ops: &[u32],
        #[case] expected_end: i32,
    ) {
        let (_, end) = mate_soft_unclipped_from_ops(mate_pos_1based, ops);
        assert_eq!(end, expected_end);
    }

    /// A malformed MC CIGAR fails closed through the whole MC-tag boundary path:
    /// `mate_soft_unclipped_from_mc` returns `None` rather than a partial boundary.
    #[rstest]
    #[case::embedded_garbage("10Mfoo5S")]
    #[case::bare_operator("MS")]
    #[case::non_ascii("10M\u{20ac}5S")]
    #[case::internal_soft("10M5S10M")]
    #[case::no_ref_op("36S")]
    #[case::empty("")]
    fn mate_soft_unclipped_from_mc_rejects_malformed(#[case] cigar: &str) {
        assert_eq!(mate_soft_unclipped_from_mc(100, cigar), None);
    }

    /// A mate record whose trailing soft clips sum past `i32::MAX` must clamp the mate
    /// boundary outward (to `i32::MAX`), not collapse to `0`. Collapsing to `0` shrank
    /// `mate_unclipped_end`, making the read look like it extends further past the mate
    /// and over-clipping the consensus.
    #[test]
    fn mate_soft_unclipped_from_ops_clamps_oversized_trailing_soft() {
        // A CIGAR op length is stored in 28 bits, so a single op maxes at 2^28 - 1; nine of
        // them sum past i32::MAX. The leading 10M anchors the ref op so the S ops count as
        // *trailing* soft clips.
        let max_op_len = (1usize << 28) - 1;
        let mut ops = vec![encode_op(0, 10)]; // 10M
        ops.extend(vec![encode_op(4, max_op_len); 9]); // 9 oversized trailing S
        let (start, end) = mate_soft_unclipped_from_ops(100, &ops);
        assert_eq!(start, 100, "no leading soft clip -> start stays at mate_pos");
        assert_eq!(
            end,
            i32::MAX,
            "oversized trailing soft clip clamps the boundary outward, not to a shrunken 109"
        );
    }

    // ========================================================================
    // is_fr_pair_raw tests
    // ========================================================================

    #[test]
    fn test_is_fr_pair_raw_not_paired() {
        // Not paired => false
        let rec = make_bam_bytes(0, 100, 0, b"rea", &[encode_op(0, 10)], 10, 0, 200, &[]);
        assert!(!is_fr_pair_raw(&rec));
    }

    #[test]
    fn test_is_fr_pair_raw_unmapped() {
        // Paired but unmapped => false
        let rec = make_bam_bytes(
            0,
            100,
            flags::PAIRED | flags::UNMAPPED,
            b"rea",
            &[encode_op(0, 10)],
            10,
            0,
            200,
            &[],
        );
        assert!(!is_fr_pair_raw(&rec));
    }

    #[test]
    fn test_is_fr_pair_raw_mate_unmapped() {
        // Paired, mapped but mate unmapped => false
        let rec = make_bam_bytes(
            0,
            100,
            flags::PAIRED | flags::MATE_UNMAPPED,
            b"rea",
            &[encode_op(0, 10)],
            10,
            -1,
            -1,
            &[],
        );
        assert!(!is_fr_pair_raw(&rec));
    }

    #[test]
    fn test_is_fr_pair_raw_different_references() {
        // Paired, both mapped, but different references => false
        let rec = make_bam_bytes(
            0,
            100,
            flags::PAIRED | flags::MATE_REVERSE,
            b"rea",
            &[encode_op(0, 10)],
            10,
            1,
            200,
            &[],
        );
        assert!(!is_fr_pair_raw(&rec));
    }

    #[test]
    fn test_is_fr_pair_raw_same_strand_ff() {
        // Paired, same reference, but both forward (FF) => false
        let rec =
            make_bam_bytes(0, 100, flags::PAIRED, b"rea", &[encode_op(0, 10)], 10, 0, 200, &[]);
        assert!(!is_fr_pair_raw(&rec));
    }

    #[test]
    fn test_is_fr_pair_raw_same_strand_rr() {
        // Paired, same reference, both reverse (RR) => false
        let rec = make_bam_bytes(
            0,
            100,
            flags::PAIRED | flags::REVERSE | flags::MATE_REVERSE,
            b"rea",
            &[encode_op(0, 10)],
            10,
            0,
            200,
            &[],
        );
        assert!(!is_fr_pair_raw(&rec));
    }

    #[test]
    fn test_is_fr_pair_raw_fr_positive_strand_read() {
        // FR pair: this read is forward, mate is reverse, on same reference
        // positive_five_prime = alignment_start = 101
        // negative_five_prime = alignment_start + insert_size = 101 + 200 = 301
        // 101 < 301 => FR => true
        let rec = make_bam_bytes_with_tlen(
            0,
            100,
            flags::PAIRED | flags::MATE_REVERSE,
            b"rea",
            &[encode_op(0, 10)],
            10,
            0,
            200,
            200,
            &[],
        );
        assert!(is_fr_pair_raw(&rec));
    }

    #[test]
    fn test_is_fr_pair_raw_fr_negative_strand_read() {
        // FR pair: this read is reverse, mate is forward
        // positive_five_prime = mate_start = 101
        // negative_five_prime = alignment_end = 101 + 10 - 1 = 110
        // Since mate at 101 < end at 110, this is FR => true
        let rec = make_bam_bytes_with_tlen(
            0,
            100,
            flags::PAIRED | flags::REVERSE,
            b"rea",
            &[encode_op(0, 10)],
            10,
            0,
            100,
            -10,
            &[],
        );
        assert!(is_fr_pair_raw(&rec));
    }

    #[test]
    fn test_is_fr_pair_raw_rf_orientation() {
        // RF pair: this read is forward, mate is reverse, but mate is upstream
        // Read at pos 200, mate at pos 100
        // positive_five_prime = alignment_start = 201
        // negative_five_prime = alignment_start + insert_size = 201 + (-100) = 101
        // 201 > 101 => NOT FR (it's RF) => false
        let rec = make_bam_bytes_with_tlen(
            0,
            200,
            flags::PAIRED | flags::MATE_REVERSE,
            b"rea",
            &[encode_op(0, 10)],
            10,
            0,
            100,
            -100,
            &[],
        );
        assert!(!is_fr_pair_raw(&rec));
    }

    // ========================================================================
    // is_primary_fr_pair_raw tests (symmetric per-pair FR classification)
    // ========================================================================

    #[test]
    fn test_is_primary_fr_pair_raw_symmetric_on_dovetail() {
        // A dovetail FR pair on which the *per-record* check disagrees:
        //  - forward read has TLEN = -90, so is_fr_pair_raw(fwd) uses the TLEN branch
        //    (0 < -90 => false) and wrongly reports NOT FR (CODEC-01);
        //  - reverse read (100M @ 61..160, mate 5' at 101) uses the CIGAR branch
        //    (101 < 160 => true) and reports FR.
        // is_primary_fr_pair_raw evaluates the reverse record only, so it is symmetric
        // and returns true regardless of argument order.
        let fwd = make_bam_bytes_with_tlen(
            0,
            100,
            flags::PAIRED | flags::MATE_REVERSE | flags::FIRST_SEGMENT,
            b"dtl",
            &[encode_op(0, 50)],
            50,
            0,
            60,
            -90,
            &[],
        );
        let rev = make_bam_bytes_with_tlen(
            0,
            60,
            flags::PAIRED | flags::REVERSE | flags::LAST_SEGMENT,
            b"dtl",
            &[encode_op(0, 100)],
            100,
            0,
            100,
            90,
            &[],
        );

        // Per-record is asymmetric — this is the bug the per-pair check avoids.
        assert!(!is_fr_pair_raw(&fwd), "forward per-record check mis-reports NOT FR");
        assert!(is_fr_pair_raw(&rev), "reverse per-record check reports FR");

        // Per-pair is symmetric and correct in both argument orders.
        assert!(is_primary_fr_pair_raw(&fwd, &rev));
        assert!(is_primary_fr_pair_raw(&rev, &fwd));
    }

    #[test]
    fn test_is_primary_fr_pair_raw_negative_cases() {
        let fr_fwd = make_bam_bytes_with_tlen(
            0,
            100,
            flags::PAIRED | flags::MATE_REVERSE,
            b"p",
            &[encode_op(0, 50)],
            50,
            0,
            150,
            100,
            &[],
        );
        // RF pair: reverse read upstream of its forward mate -> not FR.
        let rf_rev = make_bam_bytes_with_tlen(
            0,
            100,
            flags::PAIRED | flags::REVERSE,
            b"p",
            &[encode_op(0, 100)],
            100,
            0,
            200,
            100,
            &[],
        );
        let rf_fwd = make_bam_bytes_with_tlen(
            0,
            200,
            flags::PAIRED | flags::MATE_REVERSE,
            b"p",
            &[encode_op(0, 50)],
            50,
            0,
            100,
            -100,
            &[],
        );
        assert!(!is_primary_fr_pair_raw(&rf_rev, &rf_fwd), "RF pair is not a primary FR pair");
        assert!(!is_primary_fr_pair_raw(&rf_fwd, &rf_rev), "RF pair is not FR either order");

        // Same strand (FF): opposite-strand precondition fails.
        let ff_a = make_bam_bytes_with_tlen(
            0,
            10,
            flags::PAIRED,
            b"p",
            &[encode_op(0, 30)],
            30,
            0,
            50,
            70,
            &[],
        );
        let ff_b = make_bam_bytes_with_tlen(
            0,
            50,
            flags::PAIRED,
            b"p",
            &[encode_op(0, 30)],
            30,
            0,
            10,
            -70,
            &[],
        );
        assert!(!is_primary_fr_pair_raw(&ff_a, &ff_b), "same-strand pair is not FR");

        // Different references.
        let xchrom = make_bam_bytes_with_tlen(
            1,
            100,
            flags::PAIRED | flags::REVERSE,
            b"p",
            &[encode_op(0, 50)],
            50,
            0,
            100,
            0,
            &[],
        );
        assert!(!is_primary_fr_pair_raw(&fr_fwd, &xchrom), "cross-chromosomal is not FR");
    }

    // ========================================================================
    // compute_bases_past_ref_pos tests
    // ========================================================================

    #[test]
    fn test_compute_bases_past_ref_pos_simple_match() {
        // 10M starting at ref pos 100 (1-based)
        // target_ref_pos = 105: should find read_pos at offset 5
        let cigar = &[encode_op(0, 10)]; // 10M
        let result = compute_bases_past_ref_pos(cigar, 100, 105);
        assert_eq!(result, 6); // 1-based: read pos 6 at ref pos 105
    }

    #[test]
    fn test_compute_bases_past_ref_pos_at_start() {
        // 10M starting at ref pos 100
        // target_ref_pos = 100: first position
        let cigar = &[encode_op(0, 10)];
        let result = compute_bases_past_ref_pos(cigar, 100, 100);
        assert_eq!(result, 1);
    }

    #[test]
    fn test_compute_bases_past_ref_pos_at_end() {
        // 10M starting at ref pos 100
        // target_ref_pos = 109: last position
        let cigar = &[encode_op(0, 10)];
        let result = compute_bases_past_ref_pos(cigar, 100, 109);
        assert_eq!(result, 10);
    }

    #[test]
    fn test_compute_bases_past_ref_pos_past_alignment() {
        // 10M starting at ref pos 100
        // target_ref_pos = 110: beyond alignment
        let cigar = &[encode_op(0, 10)];
        let result = compute_bases_past_ref_pos(cigar, 100, 110);
        assert_eq!(result, 0); // outside alignment
    }

    #[test]
    fn test_compute_bases_past_ref_pos_with_insertion() {
        // 5M3I5M: insertion adds 3 query bases without consuming reference
        // At ref 100: 5M covers ref 100-104, 3I adds 3 query bases,
        // then 5M covers ref 105-109
        // target=107: in second 5M, offset 2 from ref 105
        // query pos = 5 (from first M) + 3 (from I) + 3 (offset in second M) = 11
        let cigar = &[encode_op(0, 5), encode_op(1, 3), encode_op(0, 5)]; // 5M3I5M
        let result = compute_bases_past_ref_pos(cigar, 100, 107);
        assert_eq!(result, 11);
    }

    #[test]
    fn test_compute_bases_past_ref_pos_in_deletion() {
        // 5M3D5M: deletion spans ref 105-107 without consuming query
        // target=106: falls in the deletion
        let cigar = &[encode_op(0, 5), encode_op(2, 3), encode_op(0, 5)]; // 5M3D5M
        let result = compute_bases_past_ref_pos(cigar, 100, 106);
        assert_eq!(result, 0); // position is in a deletion
    }

    #[test]
    fn test_compute_bases_past_ref_pos_with_soft_clip() {
        // 3S10M: soft clip consumes 3 query bases but no reference
        // Alignment starts at ref 100, so 10M covers ref 100-109
        // target=102: offset 2 in the M, query pos = 3 (from S) + 3 = 6
        let cigar = &[encode_op(4, 3), encode_op(0, 10)]; // 3S10M
        let result = compute_bases_past_ref_pos(cigar, 100, 102);
        assert_eq!(result, 6);
    }

    // ========================================================================
    // compute_bases_before_ref_pos tests
    // ========================================================================

    #[test]
    fn test_compute_bases_before_ref_pos_simple_match() {
        // 10M starting at ref pos 100
        // target_ref_pos = 105: read_pos increments to 6, but returns 6-1=5
        let cigar = &[encode_op(0, 10)];
        let result = compute_bases_before_ref_pos(cigar, 100, 105);
        assert_eq!(result, 5);
    }

    #[test]
    fn test_compute_bases_before_ref_pos_at_start() {
        // 10M starting at ref pos 100
        // target_ref_pos = 100: first position, read_pos=1, returns 1-1=0
        let cigar = &[encode_op(0, 10)];
        let result = compute_bases_before_ref_pos(cigar, 100, 100);
        assert_eq!(result, 0);
    }

    #[test]
    fn test_compute_bases_before_ref_pos_past_alignment() {
        // 10M starting at ref pos 100
        // target_ref_pos = 110: beyond alignment
        let cigar = &[encode_op(0, 10)];
        let result = compute_bases_before_ref_pos(cigar, 100, 110);
        assert_eq!(result, 0);
    }

    #[test]
    fn test_compute_bases_before_ref_pos_with_insertion() {
        // 5M3I5M: at ref 107 (in second M block)
        // query consumed: 5(M) + 3(I) + 3(into second M) = 11, returns 11-1=10
        let cigar = &[encode_op(0, 5), encode_op(1, 3), encode_op(0, 5)]; // 5M3I5M
        let result = compute_bases_before_ref_pos(cigar, 100, 107);
        assert_eq!(result, 10);
    }

    #[test]
    fn test_compute_bases_before_ref_pos_in_deletion() {
        // 5M3D5M: deletion at ref 105-107
        // target=106: falls in deletion
        let cigar = &[encode_op(0, 5), encode_op(2, 3), encode_op(0, 5)]; // 5M3D5M
        let result = compute_bases_before_ref_pos(cigar, 100, 106);
        assert_eq!(result, 0);
    }

    #[test]
    fn test_compute_bases_before_ref_pos_with_soft_clip() {
        // 3S10M: soft clip consumes 3 query bases
        // At ref 102: query = 3(S) + 3(into M) = 6, returns 6-1=5
        let cigar = &[encode_op(4, 3), encode_op(0, 10)]; // 3S10M
        let result = compute_bases_before_ref_pos(cigar, 100, 102);
        assert_eq!(result, 5);
    }

    // ========================================================================
    // num_bases_extending_past_mate_raw tests
    // ========================================================================

    #[test]
    fn test_num_bases_extending_past_mate_raw_not_paired() {
        // Not paired => 0
        let rec = make_bam_bytes(0, 100, 0, b"rea", &[encode_op(0, 10)], 10, 0, 200, &[]);
        assert_eq!(num_bases_extending_past_mate_raw(&rec), 0);
    }

    #[test]
    fn test_num_bases_extending_past_mate_raw_unmapped() {
        // Paired but unmapped => 0
        let rec = make_bam_bytes(
            0,
            100,
            flags::PAIRED | flags::UNMAPPED | flags::MATE_REVERSE,
            b"rea",
            &[encode_op(0, 10)],
            10,
            0,
            200,
            &[],
        );
        assert_eq!(num_bases_extending_past_mate_raw(&rec), 0);
    }

    #[test]
    fn test_num_bases_extending_past_mate_raw_mate_unmapped() {
        // Paired, mapped but mate unmapped => 0
        let rec = make_bam_bytes(
            0,
            100,
            flags::PAIRED | flags::MATE_UNMAPPED | flags::MATE_REVERSE,
            b"rea",
            &[encode_op(0, 10)],
            10,
            -1,
            -1,
            &[],
        );
        assert_eq!(num_bases_extending_past_mate_raw(&rec), 0);
    }

    #[test]
    fn test_num_bases_extending_past_mate_raw_same_strand() {
        // Both same strand => 0
        let mut aux = Vec::new();
        aux.extend_from_slice(b"MCZ10M\x00");
        let rec = make_bam_bytes(
            0,
            100,
            flags::PAIRED, // both forward
            b"rea",
            &[encode_op(0, 10)],
            10,
            0,
            200,
            &aux,
        );
        assert_eq!(num_bases_extending_past_mate_raw(&rec), 0);
    }

    #[test]
    fn test_num_bases_extending_past_mate_raw_different_references() {
        // Different references => 0
        let mut aux = Vec::new();
        aux.extend_from_slice(b"MCZ10M\x00");
        let rec = make_bam_bytes(
            0,
            100,
            flags::PAIRED | flags::MATE_REVERSE,
            b"rea",
            &[encode_op(0, 10)],
            10,
            1,
            200,
            &aux,
        );
        assert_eq!(num_bases_extending_past_mate_raw(&rec), 0);
    }

    #[test]
    fn test_num_bases_extending_past_mate_raw_no_mc_tag() {
        // Paired FR but no MC tag => 0
        let rec = make_bam_bytes_with_tlen(
            0,
            100,
            flags::PAIRED | flags::MATE_REVERSE,
            b"rea",
            &[encode_op(0, 10)],
            10,
            0,
            200,
            110, // TLEN: valid FR pair
            &[],
        );
        assert_eq!(num_bases_extending_past_mate_raw(&rec), 0);
    }

    #[test]
    fn test_num_bases_extending_past_mate_raw_positive_strand_overlap() {
        // Positive strand read extending past reverse mate's unclipped end
        // Read: forward at pos 100, 20M (alignment_end = 101 + 20 - 1 = 120, 1-based)
        // Mate: reverse at pos 105, MC=10M (mate unclipped_end = 106 + 10 - 1 = 115, 1-based)
        // alignment_end (120) >= mate_unclipped_end (115), so need to compute bases past ref 115
        // At ref 115 (1-based), starting from alignment_start=101:
        // offset in 20M = 115 - 101 + 1 = 15 -> read_pos 15
        // read_length = 20
        // result = 20 - 15 = 5
        let mut aux = Vec::new();
        aux.extend_from_slice(b"MCZ10M\x00");
        let rec = make_bam_bytes_with_tlen(
            0,
            100, // 0-based pos
            flags::PAIRED | flags::MATE_REVERSE,
            b"rea",
            &[encode_op(0, 20)],
            20,
            0,
            105, // 0-based mate pos
            20,  // TLEN: FR pair spanning 20 bases
            &aux,
        );
        assert_eq!(num_bases_extending_past_mate_raw(&rec), 5);
    }

    #[test]
    fn test_num_bases_extending_past_mate_raw_positive_strand_no_overlap() {
        // Positive strand read NOT extending past reverse mate's unclipped end
        // Read: forward at pos 100, 10M (alignment_end = 101 + 10 - 1 = 110)
        // Mate: reverse at pos 200, MC=10M (mate unclipped_end = 201 + 10 - 1 = 210)
        // alignment_end (110) < mate_unclipped_end (210), check trailing soft clips
        // No trailing soft clips => trailing_sc.saturating_sub(gap) = 0
        let mut aux = Vec::new();
        aux.extend_from_slice(b"MCZ10M\x00");
        let rec = make_bam_bytes_with_tlen(
            0,
            100,
            flags::PAIRED | flags::MATE_REVERSE,
            b"rea",
            &[encode_op(0, 10)],
            10,
            0,
            200,
            110, // TLEN: FR pair spanning pos 100 to mate end 209
            &aux,
        );
        assert_eq!(num_bases_extending_past_mate_raw(&rec), 0);
    }

    /// An oversized reference-consuming read CIGAR must not overflow the positive-strand
    /// `alignment_end = this_pos + ref_len - 1` computation. Eight max-length `M` ops sum to
    /// `8 * (2^28 - 1) = 2_147_483_640` (just under `i32::MAX`, so `reference_length_from_cigar`
    /// itself does not overflow), and `this_pos (101) + ref_len` then exceeds `i32::MAX` — a plain
    /// subtraction would panic in debug / wrap in release. The saturating computation instead
    /// clamps `alignment_end` to `i32::MAX`, which still exceeds the mate's unclipped end (210),
    /// so the read is measured to extend past the mate by `read_length - 110` bases.
    #[test]
    fn test_num_bases_extending_past_mate_raw_oversized_ref_cigar_saturates() {
        let max_op_len = (1usize << 28) - 1; // CIGAR op length is 28 bits
        let cigar: Vec<u32> = vec![encode_op(0, max_op_len); 8]; // 8 oversized M ops
        let read_length = 8 * max_op_len;
        let mut aux = Vec::new();
        aux.extend_from_slice(b"MCZ10M\x00");
        let rec = make_bam_bytes_with_tlen(
            0,
            100, // 0-based pos -> this_pos = 101
            flags::PAIRED | flags::MATE_REVERSE,
            b"rea",
            &cigar,
            10, // l_seq is unused by the overlap math; keep the record small
            0,
            200, // mate reverse at 200 -> mate unclipped_end = 210
            100, // TLEN > 0 so the forward read classifies as a valid FR pair
            &aux,
        );
        // read_pos at ref 210 (alignment_start 101) = 110 -> read_length - 110.
        assert_eq!(num_bases_extending_past_mate_raw(&rec), read_length - 110);
    }

    #[test]
    fn test_num_bases_extending_past_mate_raw_negative_strand_overlap() {
        // Negative strand read extending before forward mate's unclipped start
        // Read: reverse at pos 100, 20M
        // Mate: forward at pos 105, MC=10M (mate unclipped_start = 106)
        // this_pos (101) <= mate_unclipped_start (106), so compute bases before ref pos 106
        // 20M from pos 101: at ref 106, read_pos = 6, returns 6-1 = 5
        let mut aux = Vec::new();
        aux.extend_from_slice(b"MCZ10M\x00");
        let rec = make_bam_bytes(
            0,
            100, // 0-based pos
            flags::PAIRED | flags::REVERSE,
            b"rea",
            &[encode_op(0, 20)],
            20,
            0,
            105, // 0-based mate pos
            &aux,
        );
        assert_eq!(num_bases_extending_past_mate_raw(&rec), 5);
    }

    #[test]
    fn test_num_bases_extending_past_mate_raw_negative_strand_no_overlap() {
        // Negative strand read NOT extending before forward mate
        // Read: reverse at pos 200, 10M
        // Mate: forward at pos 100, MC=10M (mate unclipped_start = 101)
        // this_pos (201) > mate_unclipped_start (101), check leading soft clips
        // No leading soft clips, gap = 201 - 101 = 100, 0.saturating_sub(100) = 0
        let mut aux = Vec::new();
        aux.extend_from_slice(b"MCZ10M\x00");
        let rec = make_bam_bytes(
            0,
            200,
            flags::PAIRED | flags::REVERSE,
            b"rea",
            &[encode_op(0, 10)],
            10,
            0,
            100,
            &aux,
        );
        assert_eq!(num_bases_extending_past_mate_raw(&rec), 0);
    }

    #[test]
    fn test_num_bases_extending_past_mate_raw_negative_strand_gap_with_soft_clip() {
        // Negative strand read with soft clip, this_pos > mate_unclipped_start
        // Read: reverse at pos 110 (0-based), 3S10M (query_len=13)
        // Mate: forward at pos 105 (0-based), MC=10M (mate unclipped_start = 106)
        // this_pos = 111, mate_unclipped_start = 106
        // this_pos > mate_unclipped_start, so gap = 111 - 106 = 5
        // leading_soft_clip = 3, 3.saturating_sub(5) = 0
        let mut aux = Vec::new();
        aux.extend_from_slice(b"MCZ10M\x00");
        let rec = make_bam_bytes(
            0,
            110,
            flags::PAIRED | flags::REVERSE,
            b"rea",
            &[encode_op(4, 3), encode_op(0, 10)],
            13,
            0,
            105,
            &aux,
        );
        assert_eq!(num_bases_extending_past_mate_raw(&rec), 0);
    }

    #[test]
    fn test_num_bases_extending_past_mate_raw_reverse_oversized_mc_leading_soft_no_overflow() {
        // Regression (CodeRabbit): an adversarial MC-tag CIGAR with an oversized leading
        // soft clip saturates `leading_soft` to i32::MAX, pushing the mate's soft-only
        // unclipped start far negative. In the reverse-strand `this_pos > mate_unclipped_start`
        // branch, the gap subtraction `this_pos - mate_unclipped_start` then exceeds i32::MAX.
        // Pre-fix this panicked in debug/test (overflow) and wrapped to garbage in release;
        // it must saturate instead and clip nothing.
        //
        // Read: reverse at pos 100 (0-based) => this_pos = 101, 20M (no leading soft clip)
        // Mate: forward at pos 0 (0-based) => mate_pos_1based = 1
        //   MC = "9999999999S10M" => leading_soft saturates to i32::MAX
        //   => mate_unclipped_start = 1 - i32::MAX (well below this_pos)
        // gap = 101 - (1 - i32::MAX) overflows i32; saturating_sub clamps it, and with no
        // leading soft clip the result is leading_sc(0).saturating_sub(gap) = 0.
        let mut aux = Vec::new();
        aux.extend_from_slice(b"MCZ9999999999S10M\x00");
        let rec = make_bam_bytes(
            0,
            100,
            flags::PAIRED | flags::REVERSE,
            b"rea",
            &[encode_op(0, 20)],
            20,
            0,
            0,
            &aux,
        );
        assert_eq!(num_bases_extending_past_mate_raw(&rec), 0);
    }

    #[test]
    fn test_num_bases_extending_past_mate_raw_positive_strand_gap_with_soft_clip() {
        // Positive strand read with trailing soft clip, alignment_end < mate_unclipped_end
        // Read: forward at pos 100 (0-based), 10M3S (query_len=13)
        // Mate: reverse at pos 200 (0-based), MC=10M (mate unclipped_end = 201+10-1=210)
        // alignment_end = 101+10-1 = 110, mate_unclipped_end = 210
        // 110 < 210, gap = 210 - 110 = 100
        // trailing_soft_clip = 3, 3.saturating_sub(100) = 0
        let mut aux = Vec::new();
        aux.extend_from_slice(b"MCZ10M\x00");
        let rec = make_bam_bytes_with_tlen(
            0,
            100,
            flags::PAIRED | flags::MATE_REVERSE,
            b"rea",
            &[encode_op(0, 10), encode_op(4, 3)],
            13,
            0,
            200,
            110, // TLEN: valid FR pair
            &aux,
        );
        assert_eq!(num_bases_extending_past_mate_raw(&rec), 0);
    }

    // ========================================================================
    // trailing_soft_clip_from_ops / leading_soft_clip_from_ops tests
    // ========================================================================

    #[test]
    fn test_trailing_soft_clip_from_ops_none() {
        let cigar = &[encode_op(0, 10)]; // 10M
        assert_eq!(trailing_soft_clip_from_ops(cigar), 0);
    }

    #[test]
    fn test_trailing_soft_clip_from_ops_with_soft() {
        let cigar = &[encode_op(0, 10), encode_op(4, 5)]; // 10M5S
        assert_eq!(trailing_soft_clip_from_ops(cigar), 5);
    }

    #[test]
    fn test_trailing_soft_clip_from_ops_with_hard_after_soft() {
        let cigar = &[encode_op(0, 10), encode_op(4, 5), encode_op(5, 3)]; // 10M5S3H
        assert_eq!(trailing_soft_clip_from_ops(cigar), 5);
    }

    #[test]
    fn test_leading_soft_clip_from_ops_none() {
        let cigar = &[encode_op(0, 10)]; // 10M
        assert_eq!(leading_soft_clip_from_ops(cigar), 0);
    }

    #[test]
    fn test_leading_soft_clip_from_ops_with_soft() {
        let cigar = &[encode_op(4, 5), encode_op(0, 10)]; // 5S10M
        assert_eq!(leading_soft_clip_from_ops(cigar), 5);
    }

    #[test]
    fn test_leading_soft_clip_from_ops_with_hard_before_soft() {
        let cigar = &[encode_op(5, 3), encode_op(4, 5), encode_op(0, 10)]; // 3H5S10M
        assert_eq!(leading_soft_clip_from_ops(cigar), 5);
    }

    // ========================================================================
    // Regression: chimeric/split reads with non-FR orientation
    // ========================================================================

    #[test]
    fn test_num_bases_extending_past_mate_raw_non_fr_chimeric_forward() {
        // Reproduces the MI=807 bug from SRR6109273 simplex equivalency failure.
        // R1: flag=97 (paired, mate_reverse, first_in_pair), pos=11576620,
        //     CIGAR=145M124S, mate_pos=11576412, TLEN=-28, MC=87S182M
        // The reads are on opposite strands and same ref, but NOT in FR orientation
        // (TLEN is negative => RF orientation). The old RecordBuf code checked
        // is_fr_pair_from_tags and returned 0. The raw-byte code was missing this
        // check and incorrectly returned 269 (entire read length), causing the
        // consensus read to be trimmed to zero bases and dropped.
        let mut aux = Vec::new();
        aux.extend_from_slice(b"MCZ87S182M\x00");
        let rec = make_bam_bytes_with_tlen(
            0,
            11_576_620,                                 // 0-based pos
            flags::PAIRED | flags::MATE_REVERSE | 0x40, // flag=97 (0x40=first_in_pair)
            b"rea",
            &[encode_op(0, 145), encode_op(4, 124)], // 145M124S
            269,                                     // seq_len = 145 + 124
            0,                                       // same ref
            11_576_412,                              // 0-based mate pos
            -28,                                     // TLEN
            &aux,
        );
        // Not FR pair, so should return 0 (no clipping)
        assert_eq!(num_bases_extending_past_mate_raw(&rec), 0);
    }

    #[test]
    fn test_num_bases_extending_past_mate_raw_soft_only_mate_end() {
        // CODEC3-03: the mate boundary must count SOFT clips only (fgbio mateUnSoftClippedEnd),
        // not hard clips. Positive read 40M @100 (alignment end ref 140). Reverse mate MC=30M5H
        // at ref 101:
        //   soft-only unclipped end = 101 + 30 - 1 = 130  -> clip 40 - readPos(ref130)=30 = 10
        //   (old soft+hard end 135 would have clipped only 5).
        let mut aux = Vec::new();
        aux.extend_from_slice(b"MCZ30M5H\x00");
        let rec = make_bam_bytes_with_tlen(
            0,
            100,
            flags::PAIRED | flags::MATE_REVERSE,
            b"rea",
            &[encode_op(0, 40)], // 40M
            40,
            0,
            100,
            40, // positive TLEN -> FR
            &aux,
        );
        assert_eq!(num_bases_extending_past_mate_raw(&rec), 10);
    }

    #[test]
    fn test_num_bases_extending_past_mate_vs_mate_raw_no_mc_uses_mate() {
        // CODEC3-04: with the MC tag absent, the MC-based path under-clips (returns 0), but the
        // mate-record path derives the (soft-only) boundary from the mate in hand and still clips.
        let rec = make_bam_bytes_with_tlen(
            0,
            100,
            flags::PAIRED | flags::MATE_REVERSE,
            b"rea",
            &[encode_op(0, 40)], // 40M
            40,
            0,
            100,
            40,
            &[], // no MC tag
        );
        let mate = make_bam_bytes_with_tlen(
            0,
            100,
            flags::PAIRED | flags::REVERSE,
            b"rea",
            &[encode_op(0, 30), encode_op(5, 5)], // 30M5H
            30,
            0,
            100,
            -40,
            &[],
        );
        // MC-based path: no MC tag -> silent no-clip.
        assert_eq!(num_bases_extending_past_mate_raw(&rec), 0);
        // mate-record path: clips 10 using the mate's real soft-only boundary (ref end 130).
        assert_eq!(num_bases_extending_past_mate_vs_mate_raw(&rec, &mate), 10);
    }

    #[test]
    fn test_num_bases_extending_past_mate_vs_mate_raw_symmetric_on_dovetail() {
        // A dovetail-forward pair whose forward read has a negative TLEN, so the per-record
        // is_fr_pair_raw(fwd) wrongly reports NOT FR (htsjdk/samtools#1771). The pair IS a primary
        // FR pair, and the forward read extends 50 bases past the mate's (soft-only) end. Using the
        // symmetric per-pair guard, the clip is computed (50); the old per-record guard returned 0.
        let fwd = make_bam_bytes_with_tlen(
            0,
            100,
            flags::PAIRED | flags::MATE_REVERSE | 0x40,
            b"dtl",
            &[encode_op(0, 60)], // 60M -> ref 101..160
            60,
            0,
            60,
            -50, // negative TLEN -> is_fr_pair_raw(fwd) == false
            &[],
        );
        let mate = make_bam_bytes_with_tlen(
            0,
            60,
            flags::PAIRED | flags::REVERSE | 0x80,
            b"dtl",
            &[encode_op(0, 50)], // 50M -> ref 61..110 (soft-only end 110)
            50,
            0,
            100,
            50,
            &[],
        );
        // The per-record check is asymmetric and mis-reports the forward read as non-FR ...
        assert!(!is_fr_pair_raw(&fwd));
        // ... but the pair is a valid primary FR pair, so the clip is still computed.
        assert!(is_primary_fr_pair_raw(&fwd, &mate));
        assert_eq!(num_bases_extending_past_mate_vs_mate_raw(&fwd, &mate), 50);
    }

    #[test]
    fn test_num_bases_extending_past_mate_vs_mate_matches_mc_without_hard_clips() {
        // Sanity: when the mate has no hard clips, the mate-record path and the MC-tag path agree.
        let mut aux = Vec::new();
        aux.extend_from_slice(b"MCZ30M\x00");
        let rec = make_bam_bytes_with_tlen(
            0,
            100,
            flags::PAIRED | flags::MATE_REVERSE,
            b"rea",
            &[encode_op(0, 40)],
            40,
            0,
            100,
            40,
            &aux,
        );
        let mate = make_bam_bytes_with_tlen(
            0,
            100,
            flags::PAIRED | flags::REVERSE,
            b"rea",
            &[encode_op(0, 30)],
            30,
            0,
            100,
            -40,
            &[],
        );
        let via_mc = num_bases_extending_past_mate_raw(&rec);
        let via_mate = num_bases_extending_past_mate_vs_mate_raw(&rec, &mate);
        assert_eq!(via_mc, via_mate);
        assert_eq!(via_mc, 10); // 40M end ref 140, mate soft end ref 130 -> clip 40 - readPos(130)=30 = 10
    }

    #[test]
    fn test_num_bases_extending_past_mate_raw_non_fr_chimeric_reverse() {
        // The mate of the above read:
        // R2: flag=145 (paired, reverse, second_in_pair), pos=11576412,
        //     CIGAR=87S182M, mate_pos=11576620, TLEN=28, MC=145M124S
        let mut aux = Vec::new();
        aux.extend_from_slice(b"MCZ145M124S\x00");
        let rec = make_bam_bytes_with_tlen(
            0,
            11_576_412,                            // 0-based pos
            flags::PAIRED | flags::REVERSE | 0x80, // flag=145 (0x80=second_in_pair)
            b"rea",
            &[encode_op(4, 87), encode_op(0, 182)], // 87S182M
            269,                                    // seq_len = 87 + 182
            0,                                      // same ref
            11_576_620,                             // 0-based mate pos
            28,                                     // TLEN
            &aux,
        );
        // Not FR pair, so should return 0 (no clipping)
        assert_eq!(num_bases_extending_past_mate_raw(&rec), 0);
    }
}

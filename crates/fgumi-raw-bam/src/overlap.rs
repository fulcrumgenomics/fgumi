use crate::SamTag;
use crate::cigar::{
    consumes_ref, get_cigar_ops, query_length_from_cigar, reference_length_from_raw_bam,
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

/// FR classification for the MC-tag consensus path, evaluated per-*pair* rather than via the
/// per-record TLEN arm of [`is_fr_pair_raw`].
///
/// [`is_fr_pair_raw`]'s forward-strand arm derives the negative-strand 5' end from `TLEN`, which
/// misclassifies dovetail FR pairs whose aligned ends coincide (htsjdk/samtools#1771) and zeroes
/// the read-through clip for the forward read. The simplex/duplex callers hold the read plus its
/// `MC` tag (mate CIGAR), so the reverse-strand record's CIGAR-derived orientation — the branch
/// [`is_primary_fr_pair_raw`] evaluates — can be reconstructed for either strand:
/// - a reverse-strand read *is* the reverse record, so its own CIGAR-based arm (exactly
///   [`is_fr_pair_raw`]) is already correct and TLEN-independent;
/// - a forward-strand read's mate is the reverse record: its leftmost (5') is `mate_pos` and its
///   inclusive alignment end is `mate_pos + mate_ref_len - 1` (from the `MC` CIGAR), and the pair
///   is FR iff the forward read's own 5' (its leftmost) is `<` that end — the same
///   `positive_five_prime < negative_five_prime` test [`is_fr_pair_raw`]'s reverse arm applies.
///
/// `mate_ref_len` is the reference span of the mate's CIGAR (from the `MC` tag), clamped by
/// [`saturating_reference_length`]. This matches `is_primary_fr_pair_raw(read, mate)` on dovetail
/// pairs without needing the mate record in hand.
#[must_use]
fn is_fr_pair_with_mate_cigar_raw(bam: &[u8], mate_ref_len: i32) -> bool {
    let v = RawRecordView::new(bam);
    let flg = v.flags();

    // Paired, both read and mate mapped, same reference, opposite strands — identical to the
    // non-orientation guards in `is_fr_pair_raw`.
    if flg & flags::PAIRED == 0 {
        return false;
    }
    if flg & flags::UNMAPPED != 0 || flg & flags::MATE_UNMAPPED != 0 {
        return false;
    }
    if v.ref_id() != mate_ref_id(bam) {
        return false;
    }
    let is_reverse = flg & flags::REVERSE != 0;
    let mate_is_reverse = flg & flags::MATE_REVERSE != 0;
    if is_reverse == mate_is_reverse {
        return false;
    }

    if is_reverse {
        // This read is the reverse record; its own arm is CIGAR-derived and TLEN-independent.
        return is_fr_pair_raw(bam);
    }

    // Forward read: evaluate the reverse mate's CIGAR-derived arm. positive strand 5' (this
    // read's leftmost) < negative strand 5' (the mate's inclusive alignment end).
    let this_start = v.pos() + 1; // 1-based
    let mate_start = mate_pos(bam) + 1; // 1-based
    let mate_end = mate_start.saturating_add((mate_ref_len - 1).max(0));
    this_start < mate_end
}

/// Number of bases a read extends past its mate for FR pairs, taking the mate's alignment from
/// the read's own **MC tag** (the mate CIGAR). Returns `0` for non-FR pairs or when the MC tag
/// is absent/invalid.
///
/// The distance is measured in **query** bases: the read's query bases past the last reference
/// position it shares with its mate, less the mate's own query bases past that position (the
/// first shared position, for a negative-strand read, whose 3' end is its leftmost). Hard clips
/// never contribute — the mate's overhang counts only bases physically present in its sequence,
/// matching fgbio's soft-only `mateUnSoftClippedStart`/
/// `mateUnSoftClippedEnd`. This deliberately differs from the samtools-style soft+hard unclipped
/// 5' position used for template-coordinate sorting/grouping (see
/// [`crate::cigar::mate_unclipped_5prime`]), which follows samtools' soft+hard convention.
///
/// Used by the simplex (vanilla) and duplex consensus callers, in the place of fgbio's base
/// `UmiConsensusCaller.numBasesExtendingPastMate`. The codec caller instead uses
/// [`num_bases_extending_past_mate_vs_mate_raw`], which reads the mate's alignment from the mate
/// record in hand rather than its MC tag (see CODEC3-04).
#[must_use]
pub fn num_bases_extending_past_mate_raw(bam: &[u8]) -> usize {
    // Need the MC tag for mate CIGAR information. Parsed before the FR gate because the gate
    // now classifies orientation from the mate CIGAR (see below); with no usable MC there is
    // nothing to compute the overhang from anyway, so absent/invalid MC fails closed to 0.
    let aux = aux_data_slice(bam);
    let Some(mc_bytes) = RawTagsView::new(aux).find_string(SamTag::MC) else {
        return 0;
    };
    let Ok(mc_cigar) = std::str::from_utf8(mc_bytes) else {
        return 0;
    };

    // A malformed MC CIGAR fails closed to 0 (no clip), honoring the doc contract above
    // ("Returns 0 ... when the MC tag is absent/invalid").
    let Some(mate_ops) = parse_mc_cigar_ops(mc_cigar) else {
        return 0;
    };

    // Only applies to FR pairs. Classify per-pair from the read plus its mate CIGAR rather than
    // the per-record `is_fr_pair_raw`, whose forward-strand TLEN arm misclassifies dovetail FR
    // pairs (#839 / htsjdk/samtools#1771) and would zero the read-through clip for the forward
    // read, leaking adapter into the simplex/duplex consensus.
    if !is_fr_pair_with_mate_cigar_raw(bam, saturating_reference_length(&mate_ops)) {
        return 0;
    }
    bases_extending_past_mate(bam, mate_pos(bam) + 1, &mate_ops)
}

/// Number of bases a read extends past its mate for FR pairs, taking the mate's alignment from
/// the **mate record in hand** rather than the read's MC tag.
///
/// Codec (and any caller holding both primary reads) uses this so overlap clipping does not
/// silently no-op when the MC tag is missing — mirroring fgbio's `updateMateCigars`, which
/// backfills the mate CIGAR from the in-group mate before clipping. Measures the same query
/// distance as [`num_bases_extending_past_mate_raw`].
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
    bases_extending_past_mate(rec, mate_pos_1based, &mate_ops)
}

/// Core of the overlap-clip computation shared by the MC-based and mate-record-based entry
/// points: given the mate's 1-based leftmost position and its CIGAR ops, returns how many bases
/// at the read's 3' (sequencing) end extend past the mate.
///
/// Pure geometry: takes strand, positions, and BAM-encoded CIGAR ops (`(len << 4) | op_code`,
/// same encoding as [`crate::cigar::get_cigar_ops`]) rather than a raw record, so the count can
/// be reached from decoded CIGAR ops without a raw record in hand.
///
/// # Distances are measured in query space
///
/// Both distances are counted in **query** bases from a reference position the two reads
/// actually share, per fgbio's `SamRecordClipper` rework (fulcrumgenomics/fgbio#1090, mirrored
/// as fulcrumgenomics/fgumi#752):
///
/// 1. take the last reference position both alignments cover (the first, for a negative-strand
///    read, whose 3' end is its leftmost);
/// 2. count the query bases each read has past it;
/// 3. clip the read down to the mate's count.
///
/// The superseded formula built a *reference* coordinate — the mate's soft-only unclipped end —
/// by adding the mate's trailing soft clip, a *query* distance, to its alignment end, and then
/// looked up the read position there. Those two spaces only agree across an ungapped alignment.
/// With an indel near the read's 3' end the boundary lands somewhere the read never sequenced:
/// inside a deletion it has no read position at all, and the lookup's "no such position" was
/// then subtracted from the read length, clipping the **entire read**. An insertion moved the
/// boundary the other way and under-clipped by the insertion's length. Reading the deleted base
/// as the last one before it (htsjdk's `returnLastBaseIfDeleted`) does not fix this — it only
/// replaces one reference-space answer with another, and still carries the mate's soft clip in
/// the wrong space.
///
/// # When the read stops short of the mate's alignment
///
/// A read whose alignment does not reach its mate's has no shared reference position to anchor
/// on, and nothing but its soft-clipped tail can reach past the mate. The estimate is then the
/// historical one: extrapolate the tail one query base per reference base against the mate's
/// soft-only unclipped boundary — what fgbio's `SamRecordClipper.numBasesExtendingPastMate` has
/// done since fulcrumgenomics/fgbio#842, which added it precisely because answering 0 here made
/// `FilterConsensusReads` pass reads it should have failed.
///
/// Read-through is entirely possible in this geometry, so answering 0 here would leave adapter
/// on the read. Both alignments lie within the insert, so a short insert read from both ends
/// leaves them disjoint whenever each alignment stops short of the other's — which soft clipping
/// at a variant or a mis-called base is enough to produce. The extrapolation is exact whenever
/// the region between the two alignments is ungapped, and it is continuous with the query-space
/// count above, which reduces to exactly this expression as the shared span shrinks to a single
/// reference position.
///
/// The conflation this function otherwise avoids cannot bite here: neither alignment places a
/// base between the two, so neither can place an indel between them either. An indel in the
/// sample that falls in that gap is invisible to both reads and shifts the estimate by its own
/// length, in either direction; the result is capped by the read's soft clipping regardless, so
/// no aligned base is ever removed on the strength of it.
///
/// # When the read lies entirely on the far side of the mate
///
/// A read whose 3' end points *away* from a mate it does not overlap describes an outward-facing
/// template rather than read-through, and nothing is clipped. Measuring it against the anchor
/// like any other read would count every one of its query bases as past the mate and take its
/// aligned bases — the destructive answer fulcrumgenomics/fgumi#752 exists to remove — on a
/// geometry this function cannot verify.
///
/// It is reachable two ways, one per strand. [`is_fr_pair_raw`] mirrors htsjdk's *per-record*
/// orientation test, whose forward-strand arm reads `TLEN`, so a forward record can be accepted
/// as one end of an FR pair while the mate CIGAR it carries places the mate to its left.
/// (htsjdk 5 prefers the mate CIGAR when the MC tag is present, which is why fgbio does not admit
/// that one.) The reverse-strand arm is CIGAR-derived, but it reads the mate position recorded
/// *on this record*, while [`num_bases_extending_past_mate_vs_mate_raw`] measures against the
/// mate record in hand — so a stale mate-position field admits the mirror image.
#[must_use]
pub fn bases_extending_past_mate_ops(
    is_reverse: bool,
    this_pos_1based: i32,
    this_ops: &[u32],
    mate_pos_1based: i32,
    mate_ops: &[u32],
) -> usize {
    let read_end = alignment_end_1based(this_pos_1based, this_ops);
    let mate_end = alignment_end_1based(mate_pos_1based, mate_ops);

    if is_reverse {
        // Negative strand: the read's 3' end is its leftmost, so it is measured against the
        // mate's alignment *start*.
        if this_pos_1based > mate_end {
            // The read begins after the mate's alignment ends; see "stops short" above. The gap
            // is non-negative by construction (mate_unclipped_start <= mate_pos <= mate_end <
            // this_pos), and saturating only because an adversarial mate CIGAR can push
            // mate_unclipped_start to i32::MIN, where a plain subtraction would overflow i32 —
            // a panic in debug, garbage in release.
            let (mate_unclipped_start, _) = mate_soft_unclipped_from_ops(mate_pos_1based, mate_ops);
            let gap = this_pos_1based.saturating_sub(mate_unclipped_start).cast_unsigned() as usize;
            return leading_soft_clip_from_ops(this_ops).saturating_sub(gap);
        }
        if read_end < mate_pos_1based {
            // The mate lies entirely to the right of a read whose 3' end faces left; see "far
            // side" above.
            return 0;
        }
        // Equalize the query bases the two reads have *before* the first reference position they
        // share.
        let first_shared = this_pos_1based.max(mate_pos_1based);
        let read_before = query_bases_before_ref_pos(this_ops, this_pos_1based, first_shared);
        let mate_before = query_bases_before_ref_pos(mate_ops, mate_pos_1based, first_shared);
        read_before.saturating_sub(mate_before)
    } else {
        // Positive strand: the read's 3' end is its rightmost, so it is measured against the
        // mate's alignment *end*.
        if read_end < mate_pos_1based {
            // The read ends before the mate's alignment begins; see "stops short" above. The gap
            // is non-negative by construction (read_end < mate_pos <= mate_end <=
            // mate_unclipped_end); saturating for symmetry with the reverse branch.
            let (_, mate_unclipped_end) = mate_soft_unclipped_from_ops(mate_pos_1based, mate_ops);
            let gap = mate_unclipped_end.saturating_sub(read_end).cast_unsigned() as usize;
            return trailing_soft_clip_from_ops(this_ops).saturating_sub(gap);
        }
        if mate_end < this_pos_1based {
            // The mate lies entirely to the left of a read whose 3' end faces right; see "far
            // side" above.
            return 0;
        }
        // Equalize the query bases past the last reference position they share.
        let last_shared = read_end.min(mate_end);
        let read_past = query_bases_past_ref_pos(this_ops, this_pos_1based, last_shared);
        let mate_past = query_bases_past_ref_pos(mate_ops, mate_pos_1based, last_shared);
        read_past.saturating_sub(mate_past)
    }
}

/// Thin adapter over [`bases_extending_past_mate_ops`]: reads strand, position, and CIGAR from
/// the raw record and delegates the geometry to the pure core.
fn bases_extending_past_mate(bam: &[u8], mate_pos_1based: i32, mate_ops: &[u32]) -> usize {
    let v = RawRecordView::new(bam);
    let is_reverse = v.flags() & flags::REVERSE != 0;
    let this_pos = v.pos() + 1; // 1-based
    let cigar_ops = get_cigar_ops(bam);
    bases_extending_past_mate_ops(is_reverse, this_pos, &cigar_ops, mate_pos_1based, mate_ops)
}

/// 1-based inclusive alignment end for an alignment starting at `pos_1based`.
///
/// saturating: an oversized decoded CIGAR can drive the reference span up to `i32::MAX`, so a
/// plain `pos + ref_len - 1` would overflow i32 before clipping can fail closed — panic in
/// debug, wrap to garbage in release. The `-1` is applied *first*, matching the mate-boundary
/// helper: subtracting it after the saturating add would clamp a genuinely-at-the-ceiling end
/// inward to `i32::MAX - 1`.
fn alignment_end_1based(pos_1based: i32, cigar_ops: &[u32]) -> i32 {
    pos_1based.saturating_sub(1).saturating_add(saturating_reference_length(cigar_ops))
}

/// Reference-consuming span of `cigar_ops`, clamped to `i32::MAX` instead of overflowing.
///
/// [`crate::cigar::reference_length_from_cigar`] sums with a plain `+=`, which is fine for a
/// CIGAR decoded from a BAM record but not for one decoded from an MC tag: the tag is a
/// free-form `Z` string whose operations are never length-checked against a real reference, so
/// enough of them can overflow the sum. Every boundary this module derives is fail-closed, so
/// clamping outward is the right answer here.
fn saturating_reference_length(cigar_ops: &[u32]) -> i32 {
    let mut ref_len = 0i32;
    for &op in cigar_ops {
        if consumes_ref(op & 0xF) {
            ref_len = ref_len.saturating_add(i32::try_from(op >> 4).unwrap_or(i32::MAX));
        }
    }
    ref_len
}

/// Mate's soft-only unclipped start/end (1-based) from the mate's CIGAR **ops**. Only soft clips
/// count — hard clips are excluded — matching fgbio `mateUnSoftClippedStart`/`mateUnSoftClippedEnd`.
///
/// Used only by [`bases_extending_past_mate_ops`]'s stops-short branch: these boundaries mix a query
/// distance (the soft clip) into a reference coordinate, which is sound only when no aligned base
/// lies between the two reads.
fn mate_soft_unclipped_from_ops(mate_pos_1based: i32, mate_ops: &[u32]) -> (i32, i32) {
    // Both soft-clip totals clamp to i32::MAX (not 0) when they overflow i32: an oversized
    // trailing clip that collapsed to 0 would shrink mate_unclipped_end, making the read look
    // like it extends further past the mate and over-clipping the consensus. Fail closed by
    // clamping the boundary outward, matching leading_soft.
    let leading_soft = i32::try_from(leading_soft_clip_from_ops(mate_ops)).unwrap_or(i32::MAX);
    let trailing_soft = i32::try_from(trailing_soft_clip_from_ops(mate_ops)).unwrap_or(i32::MAX);
    let ref_len = saturating_reference_length(mate_ops);
    // The inclusive-end `-1` is applied first; see `alignment_end_1based`.
    (
        mate_pos_1based.saturating_sub(leading_soft),
        mate_pos_1based.saturating_sub(1).saturating_add(ref_len).saturating_add(trailing_soft),
    )
}

/// The largest run length a BAM CIGAR operation can carry: the length occupies the upper 28 bits
/// of the packed `u32`. An MC-tag run length beyond it cannot describe a real alignment.
const MAX_CIGAR_OP_LEN: u32 = (1 << 28) - 1;

/// Parses an MC-tag CIGAR **string** into BAM-packed CIGAR ops (`(len << 4) | op_code`), the same
/// encoding [`crate::cigar::get_cigar_ops`] yields for a record's own CIGAR, so the MC-tag and
/// mate-record paths share one set of CIGAR walkers.
///
/// Returns `None` for any CIGAR that is not a structurally valid SAM CIGAR for a mapped record,
/// so malformed MC-tag input fails closed rather than yielding a partial answer. Rejected:
/// - an unknown operator byte (e.g. the `f`/`o` in `10Mfoo5S`, or a non-ASCII byte);
/// - an operator with no preceding run-length (a bare `MS`) or a zero run-length (`0M`, `10M0S`);
/// - a run-length past [`MAX_CIGAR_OP_LEN`], which no BAM CIGAR operation can hold;
/// - a trailing run-length with no operator (`10M5`), or an empty CIGAR;
/// - a soft (`S`) or hard (`H`) clip that is not at the end of the CIGAR — `S` may only sit at the
///   ends (inside any `H`), and `H` only as the first/last operation (e.g. `10M5S10M` is invalid);
/// - a CIGAR with no reference-consuming operation (`10S`), which a mapped mate cannot have.
fn parse_mc_cigar_ops(cigar: &str) -> Option<Vec<u32>> {
    // Phase 1: tokenize into (run_length, operator) pairs. Accumulate the run-length from digit
    // bytes directly rather than slicing the string by byte offset — slicing on a non-ASCII byte
    // (a malformed CIGAR) would panic on a char boundary. Reject lexical errors here.
    let mut tokens: Vec<(u32, u8)> = Vec::new();
    let mut num = 0u32;
    let mut have_digits = false;
    for &c in cigar.as_bytes() {
        if c.is_ascii_digit() {
            // saturating, then range-checked: an adversarial MC CIGAR can run the digits past
            // any bound, and clamping rather than wrapping keeps the check below decisive.
            num = num.saturating_mul(10).saturating_add(u32::from(c - b'0'));
            if num > MAX_CIGAR_OP_LEN {
                return None;
            }
            have_digits = true;
            continue;
        }
        // Every operator needs a positive run-length, and must be a known CIGAR operator.
        if !have_digits || num == 0 || cigar_op_code(c).is_none() {
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

    // Phase 2: validate structure and pack the ops in one pass. Soft clips must sit at the ends
    // (inside any hard clip); hard clips only as the first/last operation.
    let last = tokens.len() - 1;
    let mut ops = Vec::with_capacity(tokens.len());
    let mut saw_ref_op = false;
    for (i, &(len, op)) in tokens.iter().enumerate() {
        // Phase 1 accepted every operator, so this lookup cannot fail.
        let op_code = cigar_op_code(op)?;
        match op {
            b'M' | b'D' | b'N' | b'=' | b'X' => saw_ref_op = true,
            b'I' | b'P' => {} // interior op: no placement constraint
            b'S' => {
                let leading = tokens[..i].iter().all(|&(_, o)| o == b'H');
                let trailing = tokens[i + 1..].iter().all(|&(_, o)| o == b'H');
                if !leading && !trailing {
                    return None; // internal soft clip: invalid placement
                }
            }
            b'H' if i == 0 || i == last => {} // hard clip only permitted at the ends
            // An internal hard clip is invalid placement; the tokenizer already rejected unknown
            // operators, so every other operator byte here also fails closed.
            _ => return None,
        }
        // `len <= MAX_CIGAR_OP_LEN` (28 bits), so the shift cannot drop bits.
        ops.push((len << 4) | op_code);
    }
    // A mapped mate's CIGAR must consume reference; a fully-clipped CIGAR (`10S`) cannot.
    if !saw_ref_op {
        return None;
    }

    Some(ops)
}

/// BAM operation code for a CIGAR operator character, in `BAM_CIGAR_STR` (`MIDNSHP=X`) order, or
/// `None` for a byte that is not a CIGAR operator. Doubles as the tokenizer's operator check, so
/// the set of accepted operators is defined in exactly one place.
fn cigar_op_code(op: u8) -> Option<u32> {
    Some(match op {
        b'M' => 0,
        b'I' => 1,
        b'D' => 2,
        b'N' => 3,
        b'S' => 4,
        b'H' => 5,
        b'P' => 6,
        b'=' => 7,
        b'X' => 8,
        _ => return None,
    })
}

/// Which side of `target_ref_pos` the base sitting *on* it belongs to.
#[derive(Clone, Copy)]
enum RefPosBoundary {
    /// Count the query base aligned to `target_ref_pos` itself.
    Inclusive,
    /// Stop one reference base short of `target_ref_pos`.
    Exclusive,
}

/// Number of query bases the alignment places at or before `target_ref_pos` (`Inclusive`), or
/// strictly before it (`Exclusive`).
///
/// A query-only operation (`I`/`S`) occupies the gap between the reference position last consumed
/// and the next one, so it belongs to whichever side of `target_ref_pos` that gap falls on — a
/// leading soft clip precedes the alignment start and always counts, a trailing soft clip follows
/// the alignment end and never does, and an interior insertion counts only when the alignment has
/// not yet reached past the target. Both modes agree on those; they differ only in whether the
/// base aligned to `target_ref_pos` itself is counted.
///
/// A `target_ref_pos` inside a deletion is *not* a special case: the deletion consumes reference
/// without query, so the walk simply carries past the target and every query base after the
/// deletion falls on its far side. That is the whole point — this is what the superseded
/// read-position lookup could not express, and where it returned "no such position" instead.
fn query_bases_up_to_ref_pos(
    cigar_ops: &[u32],
    alignment_start_1based: i32,
    target_ref_pos: i32,
    boundary: RefPosBoundary,
) -> usize {
    // Widened to i64 so an alignment start far below the target (an adversarial mate position)
    // cannot overflow the span subtraction below.
    let target = i64::from(target_ref_pos);
    // Whether the reference base at `target` is on this side of the boundary.
    let inclusive = match boundary {
        RefPosBoundary::Inclusive => 1i64,
        RefPosBoundary::Exclusive => 0i64,
    };

    let mut ref_pos = alignment_start_1based;
    let mut query_bases = 0usize;
    for &op in cigar_ops {
        // Everything from here on sits past the target, on the reference and in the query.
        if i64::from(ref_pos) > target {
            break;
        }
        let op_type = op & 0xF;
        let op_len = (op >> 4) as usize;
        match op_type {
            0 | 7 | 8 => {
                // M, =, X: consume both query and reference. Take only the bases on this side of
                // the boundary; if any are left over, the rest of the CIGAR is past it.
                let span = target - i64::from(ref_pos) + inclusive;
                let take = op_len.min(usize::try_from(span.max(0)).unwrap_or(usize::MAX));
                query_bases += take;
                ref_pos = ref_pos.saturating_add(i32::try_from(op >> 4).unwrap_or(i32::MAX));
                if take < op_len {
                    break;
                }
            }
            // I, S: consume query only, in the gap the walk has reached — which is on this side
            // of the boundary, or the loop would have broken above.
            1 | 4 => query_bases += op_len,
            // D, N: consume reference only.
            2 | 3 => ref_pos = ref_pos.saturating_add(i32::try_from(op >> 4).unwrap_or(i32::MAX)),
            _ => {} // H, P: neither
        }
    }

    query_bases
}

/// Number of query bases the alignment places strictly past `target_ref_pos` — the 3' overhang of
/// a positive-strand read measured from a shared reference position.
fn query_bases_past_ref_pos(
    cigar_ops: &[u32],
    alignment_start_1based: i32,
    target_ref_pos: i32,
) -> usize {
    let at_or_before = query_bases_up_to_ref_pos(
        cigar_ops,
        alignment_start_1based,
        target_ref_pos,
        RefPosBoundary::Inclusive,
    );
    query_length_from_cigar(cigar_ops).saturating_sub(at_or_before)
}

/// Number of query bases the alignment places strictly before `target_ref_pos` — the 3' overhang
/// of a negative-strand read, whose 3' end is its leftmost.
fn query_bases_before_ref_pos(
    cigar_ops: &[u32],
    alignment_start_1based: i32,
    target_ref_pos: i32,
) -> usize {
    query_bases_up_to_ref_pos(
        cigar_ops,
        alignment_start_1based,
        target_ref_pos,
        RefPosBoundary::Exclusive,
    )
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
    // bases_extending_past_mate_ops tests
    // ========================================================================

    #[test]
    fn ops_core_matches_record_entry_for_indel_pair() {
        // 100bp read, positive strand, 70M10I20M then read-through; mate 80S20M at 1019.
        // Encode via testutil::encode_op (M=0,I=1,S=4).
        use crate::testutil::encode_op;
        let this_ops = vec![encode_op(0, 70), encode_op(1, 10), encode_op(0, 20)];
        let mate_ops = vec![encode_op(4, 80), encode_op(0, 20)];
        // Positive strand read at 1000, mate at 1019.
        let via_ops = bases_extending_past_mate_ops(false, 1000, &this_ops, 1019, &mate_ops);
        // The last shared reference position is the mate's alignment end: mate 80S20M@1019 ends
        // at 1019 + 20 - 1 = 1038, and the read 70M10I20M@1000 spans ref 1000-1089. Through ref
        // 1038 the read has aligned 39 query bases (ref 1000-1038); the remaining 61 (31 aligned
        // bases over ref 1039-1069, the 10I insertion, and the final 20M) extend past the mate,
        // which has none, so the clip is 100 - 39 = 61.
        assert_eq!(via_ops, 61, "exact query-space clip past the mate");

        // The private adapter must agree with the pure ops-core when given a raw record
        // carrying the same strand/pos/cigar: this proves it extracts is_reverse/pos/cigar
        // correctly and delegates to the core unchanged.
        let rec = make_bam_bytes_with_tlen(
            0,
            999, // 0-based pos (1-based 1000)
            flags::PAIRED,
            b"rea",
            &this_ops,
            100,
            0,
            1018, // 0-based mate pos (1-based 1019)
            0,
            &[],
        );
        let via_adapter = bases_extending_past_mate(&rec, 1019, &mate_ops);
        assert_eq!(via_ops, via_adapter, "adapter must agree with the pure ops-core");
    }

    // ========================================================================
    // parse_mc_cigar_ops tests
    // ========================================================================

    /// Summarizes parsed MC ops as `(leading_soft, ref_len, trailing_soft)` — counting only soft
    /// clips (H excluded) — so the case table below reads in the terms the boundary math uses.
    fn mc_summary(cigar: &str) -> Option<(i32, i32, i32)> {
        parse_mc_cigar_ops(cigar).map(|ops| {
            (
                i32::try_from(leading_soft_clip_from_ops(&ops)).unwrap_or(i32::MAX),
                saturating_reference_length(&ops),
                i32::try_from(trailing_soft_clip_from_ops(&ops)).unwrap_or(i32::MAX),
            )
        })
    }

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
    // A CIGAR whose ref ops sum past i32::MAX must saturate rather than overflow the
    // accumulation. Each op is at the 28-bit CIGAR maximum, so this is the widest span any
    // structurally valid CIGAR can describe; nine of them exceed i32::MAX.
    #[case::many_max_length_ops_saturate(&"268435455M".repeat(9), Some((0, i32::MAX, 0)))]
    // --- lexically malformed: unknown byte, missing/zero/oversized run-length, dangling, empty ---
    // A run length past the 28-bit CIGAR field cannot describe a real alignment.
    #[case::out_of_range_run_length_rejected("2000000000M", None)]
    #[case::out_of_range_run_length_absurd_rejected("9999999999M", None)]
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
    fn test_parse_mc_cigar_ops(#[case] cigar: &str, #[case] expected: Option<(i32, i32, i32)>) {
        assert_eq!(mc_summary(cigar), expected);
    }

    /// A crafted (valid-UTF-8) MC-tag CIGAR whose reference span saturates `ref_len` must not
    /// overflow the boundary arithmetic (`mate_pos_1based - 1 + ref_len + trailing_soft`).
    /// Pre-fix this panicked in debug/test and wrapped to a bogus negative boundary in release,
    /// corrupting the overlap-clip amount. A well-formed BAM cannot reach this (real ref lengths
    /// are chromosome-bounded), but the MC Z-string is never length-checked.
    #[test]
    fn mate_soft_unclipped_saturates_on_oversized_cigar() {
        let ops = parse_mc_cigar_ops(&"268435455M".repeat(9))
            .expect("well-formed (if oversized) CIGAR parses");
        let (start, end) = mate_soft_unclipped_from_ops(100, &ops);
        assert_eq!(start, 100, "no leading soft clip -> start stays at mate_pos");
        assert!(end >= start, "end boundary must stay >= start, not wrap negative");
        assert_eq!(end, i32::MAX, "saturated ref_len clamps the end boundary");
    }

    /// The inclusive-end `-1` must be applied before the saturating additions in the
    /// mate-boundary helper.
    ///
    /// Applied last it saturates the sum to `i32::MAX` and then walks it back, so a
    /// boundary whose exact value *is* `i32::MAX` was reported as `i32::MAX - 1` --
    /// clamping the fail-closed boundary inward by one base rather than outward.
    /// This is the `overlap.rs` sibling of `cigar::unclipped_other_end`.
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

    /// A malformed MC CIGAR fails closed through the whole MC-tag path: `parse_mc_cigar_ops`
    /// returns `None` rather than a partial CIGAR, and the entry point then clips nothing.
    #[rstest]
    #[case::embedded_garbage("10Mfoo5S")]
    #[case::bare_operator("MS")]
    #[case::non_ascii("10M\u{20ac}5S")]
    #[case::internal_soft("10M5S10M")]
    #[case::no_ref_op("36S")]
    #[case::empty("")]
    #[case::out_of_range_run_length("9999999999M")]
    fn mc_cigar_path_rejects_malformed(#[case] cigar: &str) {
        assert_eq!(parse_mc_cigar_ops(cigar), None);

        // A read whose mate CIGAR cannot be trusted is not clipped at all.
        let mut aux = Vec::new();
        aux.extend_from_slice(b"MCZ");
        aux.extend_from_slice(cigar.as_bytes());
        aux.push(0);
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
        assert_eq!(num_bases_extending_past_mate_raw(&rec), 0);
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
    // query_bases_past_ref_pos tests
    // ========================================================================

    /// Query bases lying strictly past a reference position. Query-only operations belong to the
    /// side of the boundary their gap falls on: a leading soft clip never counts, a trailing soft
    /// clip always does, and an interior insertion counts once the walk is past the target.
    #[rstest]
    // 10M at 100-109; ref 105 is the 6th base, leaving 4 past it.
    #[case::simple_match(&[encode_op(0, 10)], 100, 105, 4)]
    // The first aligned base leaves the other 9 past it.
    #[case::at_start(&[encode_op(0, 10)], 100, 100, 9)]
    // The last aligned base leaves nothing past it.
    #[case::at_end(&[encode_op(0, 10)], 100, 109, 0)]
    // A target past the whole alignment leaves nothing past it.
    #[case::past_alignment(&[encode_op(0, 10)], 100, 110, 0)]
    // 5M3I5M at 100-109 (13 query bases): ref 107 is the 3rd base of the second 5M, so 11 query
    // bases are at or before it (5M + 3I + 3M) and 2 are past.
    #[case::with_insertion(&[encode_op(0, 5), encode_op(1, 3), encode_op(0, 5)], 100, 107, 2)]
    // 5M3D5M at 100-112 with the deletion spanning 105-107: a target *inside* the deletion has no
    // aligned query base, and every base after the deletion is past it. This is the fgumi#752
    // case -- the superseded read-position lookup answered "no such position" (0) here, which the
    // caller then read as "the whole read is past the mate".
    #[case::in_deletion(&[encode_op(0, 5), encode_op(2, 3), encode_op(0, 5)], 100, 106, 5)]
    // 3S10M at 100-109 (13 query bases): the leading soft clip is never past the target, and ref
    // 102 is the 6th query base, leaving 7 past it.
    #[case::with_soft_clip(&[encode_op(4, 3), encode_op(0, 10)], 100, 102, 7)]
    fn test_query_bases_past_ref_pos(
        #[case] cigar: &[u32],
        #[case] alignment_start: i32,
        #[case] target: i32,
        #[case] expected: usize,
    ) {
        assert_eq!(query_bases_past_ref_pos(cigar, alignment_start, target), expected);
    }

    // ========================================================================
    // query_bases_before_ref_pos tests
    // ========================================================================

    /// Query bases lying strictly before a reference position — the mirror image, used for the
    /// negative strand, whose 3' end is its leftmost.
    #[rstest]
    // 10M at 100-109; 5 bases precede ref 105.
    #[case::simple_match(&[encode_op(0, 10)], 100, 105, 5)]
    // Nothing precedes the first aligned base.
    #[case::at_start(&[encode_op(0, 10)], 100, 100, 0)]
    // A target past the whole alignment has all 10 bases before it.
    #[case::past_alignment(&[encode_op(0, 10)], 100, 110, 10)]
    // 5M3I5M: 5M + 3I + 2 more aligned bases precede ref 107.
    #[case::with_insertion(&[encode_op(0, 5), encode_op(1, 3), encode_op(0, 5)], 100, 107, 10)]
    // 5M3D5M with the deletion spanning 105-107: only the 5 bases before the deletion precede a
    // target inside it.
    #[case::in_deletion(&[encode_op(0, 5), encode_op(2, 3), encode_op(0, 5)], 100, 106, 5)]
    // 3S10M: the leading soft clip precedes the alignment, so it always counts.
    #[case::with_soft_clip(&[encode_op(4, 3), encode_op(0, 10)], 100, 102, 5)]
    // A leading soft clip at the alignment start still precedes it -- the case that makes this
    // the mirror of `past`, not simply `past` with the target shifted by one.
    #[case::leading_soft_clip_at_start(&[encode_op(4, 3), encode_op(0, 10)], 100, 100, 3)]
    fn test_query_bases_before_ref_pos(
        #[case] cigar: &[u32],
        #[case] alignment_start: i32,
        #[case] target: i32,
        #[case] expected: usize,
    ) {
        assert_eq!(query_bases_before_ref_pos(cigar, alignment_start, target), expected);
    }

    /// Every query base falls on exactly one side of a boundary, whichever mode is asked for.
    #[rstest]
    #[case::simple_match(&[encode_op(0, 10)], 100, 105)]
    #[case::with_insertion(&[encode_op(0, 5), encode_op(1, 3), encode_op(0, 5)], 100, 107)]
    #[case::in_deletion(&[encode_op(0, 5), encode_op(2, 3), encode_op(0, 5)], 100, 106)]
    #[case::both_soft_clips(&[encode_op(4, 3), encode_op(0, 10), encode_op(4, 4)], 100, 104)]
    #[case::hard_clips_excluded(&[encode_op(5, 6), encode_op(0, 10), encode_op(5, 2)], 100, 104)]
    fn query_bases_partition_the_read(
        #[case] cigar: &[u32],
        #[case] alignment_start: i32,
        #[case] target: i32,
    ) {
        let query_length = query_length_from_cigar(cigar);
        let at_or_before =
            query_bases_up_to_ref_pos(cigar, alignment_start, target, RefPosBoundary::Inclusive);
        assert_eq!(
            at_or_before + query_bases_past_ref_pos(cigar, alignment_start, target),
            query_length
        );
        // The base aligned to the target itself is the only difference between the two modes.
        let before = query_bases_before_ref_pos(cigar, alignment_start, target);
        assert!(before <= at_or_before, "`before` counts a subset of `at or before`");
        assert!(
            at_or_before - before <= 1,
            "at most the single base aligned to the target separates the two modes"
        );
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
        // Regression (CodeRabbit): an adversarial MC-tag CIGAR with an enormous leading soft
        // clip pushed the mate's soft-only unclipped start far negative, and the reverse-strand
        // gap subtraction `this_pos - mate_unclipped_start` then exceeded i32::MAX -- a debug
        // overflow panic, garbage in release. A run length that large is no longer accepted at
        // all (it cannot fit a BAM CIGAR's 28-bit length field), so the MC-tag path now fails
        // closed one step earlier; either way the read must be left unclipped.
        //
        // Read: reverse at pos 100 (0-based) => this_pos = 101, 20M (no leading soft clip)
        // Mate: forward at pos 0 (0-based) => mate_pos_1based = 1, MC = "9999999999S10M"
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

    /// The saturation the MC-tag path can no longer reach is still reachable from a mate
    /// *record*, whose CIGAR can hold arbitrarily many maximum-length soft clips: nine of them
    /// sum past `i32::MAX`, pushing the mate's soft-only unclipped start far negative, and the
    /// reverse-strand gap subtraction `this_pos - mate_unclipped_start` must saturate rather
    /// than overflow (a debug panic, garbage in release).
    #[test]
    fn test_num_bases_extending_past_mate_vs_mate_raw_reverse_oversized_leading_soft_no_overflow() {
        let max_op_len = (1usize << 28) - 1; // a CIGAR op length is 28 bits
        let mut mate_ops = vec![encode_op(4, max_op_len); 9]; // 9 maximal leading soft clips
        mate_ops.push(encode_op(0, 10)); // 10M
        let mate = make_bam_bytes_with_tlen(
            0,
            0,
            flags::PAIRED | flags::MATE_REVERSE,
            b"rea",
            &mate_ops,
            10, // l_seq is unused by the overlap math; keep the record small
            0,
            100,
            110,
            &[],
        );
        let rec = make_bam_bytes_with_tlen(
            0,
            100, // this_pos = 101, past the mate's alignment (1-10): no shared reference position
            flags::PAIRED | flags::REVERSE,
            b"rea",
            &[encode_op(0, 20)], // no leading soft clip -> nothing to clip
            20,
            0,
            0,
            -110,
            &[],
        );
        assert_eq!(num_bases_extending_past_mate_vs_mate_raw(&rec, &mate), 0);
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

    /// #839 regression: the MC-tag path (`num_bases_extending_past_mate_raw`, used by the
    /// simplex/duplex consensus callers) must trim read-through on a dovetail FR pair whose
    /// FORWARD read has a negative TLEN. The per-record `is_fr_pair_raw(fwd)` mis-reports the
    /// forward read as non-FR (htsjdk/samtools#1771), which used to zero the clip and leak the
    /// read-through into the consensus. The pair is a valid primary FR pair, so the 40 bases of
    /// read-through past the mate must still be trimmed — matching the symmetric mate-record
    /// path `num_bases_extending_past_mate_vs_mate_raw`.
    #[test]
    fn test_num_bases_extending_past_mate_raw_symmetric_on_dovetail_forward() {
        // forward 50M50S @ 1-based 101, TLEN = -90, MC = 100M; reverse 100M @ 1-based 61.
        // 50M covers ref 101..150; the 50S read-through covers ref 151..200 in read space.
        // The mate (100M) ends at ref 160, so 40 soft bases (ref 161..200) extend past it.
        let mut aux = Vec::new();
        aux.extend_from_slice(b"MCZ100M\x00");
        let fwd = make_bam_bytes_with_tlen(
            0,
            100, // 1-based 101
            flags::PAIRED | flags::MATE_REVERSE | 0x40,
            b"dtl",
            &[encode_op(0, 50), encode_op(4, 50)], // 50M50S
            100,
            0,
            60,  // mate 1-based 61
            -90, // negative TLEN -> is_fr_pair_raw(fwd) == false
            &aux,
        );
        // The reverse mate (100M @ 1-based 61), so the MC-tag path can be pinned against the
        // symmetric mate-record path on the exact same geometry.
        let mate = make_bam_bytes_with_tlen(
            0,
            60, // 1-based 61
            flags::PAIRED | flags::REVERSE | 0x80,
            b"dtl",
            &[encode_op(0, 100)], // 100M -> ref 61..160
            100,
            0,
            100, // mate (forward) 1-based 101
            90,  // positive TLEN
            &[],
        );

        // The per-record check mis-reports the forward read as non-FR ...
        assert!(!is_fr_pair_raw(&fwd));
        // ... but the pair is a valid primary FR pair, so both the MC-tag path (simplex/duplex)
        // and the symmetric mate-record path (CODEC / clip) must trim the same 40 bases of
        // read-through. Pinning them together makes the parity claim load-bearing.
        assert!(is_primary_fr_pair_raw(&fwd, &mate));
        assert_eq!(num_bases_extending_past_mate_raw(&fwd), 40);
        assert_eq!(num_bases_extending_past_mate_vs_mate_raw(&fwd, &mate), 40);
    }

    // ========================================================================
    // Reference-vs-query distance conflation (fgumi#752 / fgbio#1090)
    // ========================================================================

    /// An FR pair that reads through into the adapter with a 1-base deletion three bases
    /// from the positive read's 3' end — the shape from fulcrumgenomics/fgumi#752.
    ///
    /// - positive: `2S124M1D3M` at 1-based 83,585,781 (reference span 83,585,781-83,585,908)
    /// - negative: `3S124M2S`   at 1-based 83,585,780 (reference span 83,585,780-83,585,903)
    ///
    /// Both reads are 129 bases. The negative read's soft-only unclipped *reference* end is
    /// 83,585,903 + 2 = 83,585,905 — a reference coordinate built by adding a **query**
    /// distance (its 2-base trailing soft clip). That coordinate lands squarely inside the
    /// positive read's `1D`, which is exactly the collision fgbio#1090 describes.
    fn deletion_at_boundary_pair() -> (Vec<u8>, Vec<u8>) {
        let mut fwd_aux = Vec::new();
        fwd_aux.extend_from_slice(b"MCZ3S124M2S\x00");
        let fwd = make_bam_bytes_with_tlen(
            0,
            83_585_780, // 0-based -> 1-based 83,585,781
            flags::PAIRED | flags::MATE_REVERSE | flags::FIRST_SEGMENT,
            b"del",
            &[encode_op(4, 2), encode_op(0, 124), encode_op(2, 1), encode_op(0, 3)],
            129,
            0,
            83_585_779,
            124,
            &fwd_aux,
        );

        let mut rev_aux = Vec::new();
        rev_aux.extend_from_slice(b"MCZ2S124M1D3M\x00");
        let rev = make_bam_bytes_with_tlen(
            0,
            83_585_779, // 0-based -> 1-based 83,585,780
            flags::PAIRED | flags::REVERSE | flags::LAST_SEGMENT,
            b"del",
            &[encode_op(4, 3), encode_op(0, 124), encode_op(4, 2)],
            129,
            0,
            83_585_780,
            -124,
            &rev_aux,
        );

        (fwd, rev)
    }

    /// The clip point must be a **query** distance measured from a reference position the
    /// two reads share, not a query distance added to a reference coordinate.
    ///
    /// The two reads share reference positions up to 83,585,903 (the negative read's
    /// alignment end). Past it the positive read has 4 query bases (`129 - 125`, since its
    /// `2S` plus 123 aligned bases carry it to query position 125) and the negative read
    /// has 2 (its trailing soft clip), so the positive read must give up `4 - 2 = 2` bases.
    ///
    /// Before the fix the boundary was the negative read's soft-only unclipped reference end
    /// (83,585,905), which falls inside the positive read's `1D`; the read-position lookup
    /// reported "no such position" and the clip became the read's entire 129 bases.
    #[test]
    fn test_deletion_at_mate_boundary_clips_the_query_distance() {
        let (fwd, rev) = deletion_at_boundary_pair();
        assert_eq!(num_bases_extending_past_mate_raw(&fwd), 2, "MC-tag path (simplex/duplex)");
        assert_eq!(
            num_bases_extending_past_mate_vs_mate_raw(&fwd, &rev),
            2,
            "mate-record path (codec)"
        );
    }

    /// The negative read of the same pair is clipped by the mirror-image rule: the two reads
    /// first share reference position 83,585,781 (the positive read's alignment start), before
    /// which the negative read has 4 query bases and the positive read has 2, so the negative
    /// read gives up 2. It was already correct — no deletion sits at *its* boundary — and must
    /// stay that way.
    #[test]
    fn test_deletion_at_mate_boundary_leaves_the_unaffected_mate_alone() {
        let (fwd, rev) = deletion_at_boundary_pair();
        assert_eq!(num_bases_extending_past_mate_raw(&rev), 2);
        assert_eq!(num_bases_extending_past_mate_vs_mate_raw(&rev, &fwd), 2);
    }

    /// The whole-read over-clip must be gone at the source, not merely clamped downstream.
    ///
    /// Consumers clamp the clip amount to the read length, so the pre-fix answer (129 of 129
    /// bases) did not panic or truncate a buffer — it silently deleted the read from its
    /// family. Pinning "strictly less than the read length" fails loudly if the boundary math
    /// ever regresses to answering with the whole read again.
    #[test]
    fn test_deletion_at_mate_boundary_does_not_consume_the_whole_read() {
        let (fwd, rev) = deletion_at_boundary_pair();
        let read_length = 129;
        assert!(
            num_bases_extending_past_mate_vs_mate_raw(&fwd, &rev) < read_length,
            "an over-clip that reaches the read length is the fgumi#752 failure"
        );
    }

    /// fgbio#1090's own worked example, which is the insertion half of the same defect:
    ///
    /// - positive: `70M10I23M47S` at 1-based 1,234,500 (reference span 1,234,500-1,234,592)
    /// - negative: `50S70M30S`    at 1-based 1,234,500 (reference span 1,234,500-1,234,569)
    ///
    /// The reads last share reference position 1,234,569; past it the positive read has 80
    /// query bases (`10I` + `23M` + `47S`) and the negative read has 30 (its trailing soft
    /// clip), so the positive read must give up 50 — fgbio#1090's stated answer, leaving
    /// `70M10I20M50S`.
    ///
    /// Before the fix the negative read's soft-only unclipped reference end (1,234,599) sat
    /// *past* the positive read's alignment end, so the clip was estimated from the trailing
    /// soft clip alone (`47 - 7 = 40`) and under-clipped by exactly the insertion's 10 bases.
    #[test]
    fn test_insertion_before_mate_boundary_is_not_under_clipped() {
        let mut fwd_aux = Vec::new();
        fwd_aux.extend_from_slice(b"MCZ50S70M30S\x00");
        let fwd = make_bam_bytes_with_tlen(
            0,
            1_234_499,
            flags::PAIRED | flags::MATE_REVERSE | flags::FIRST_SEGMENT,
            b"ins",
            &[encode_op(0, 70), encode_op(1, 10), encode_op(0, 23), encode_op(4, 47)],
            150,
            0,
            1_234_499,
            93,
            &fwd_aux,
        );
        let rev = make_bam_bytes_with_tlen(
            0,
            1_234_499,
            flags::PAIRED | flags::REVERSE | flags::LAST_SEGMENT,
            b"ins",
            &[encode_op(4, 50), encode_op(0, 70), encode_op(4, 30)],
            150,
            0,
            1_234_499,
            -93,
            &[],
        );
        assert_eq!(num_bases_extending_past_mate_raw(&fwd), 50);
        assert_eq!(num_bases_extending_past_mate_vs_mate_raw(&fwd, &rev), 50);
    }

    /// Control: with no indel anywhere in either read, the clip is unchanged. The equalization
    /// is exact for an ungapped overlap — every query base past the shared reference position
    /// is one reference base past it — so an indel-free pair must produce the same answer it
    /// always did, from both entry points and on both strands.
    ///
    /// `40M` at 101-140 against `10S30M10S` at 101-130: the reads last share 130, past which
    /// the positive read has 10 query bases and the negative read has its 10-base trailing
    /// soft clip, so nothing is clipped. They first share 101, before which the positive read
    /// has 0 and the negative read has 10, so the negative read gives up 10.
    #[test]
    fn test_ungapped_overlap_is_unchanged() {
        let mut fwd_aux = Vec::new();
        fwd_aux.extend_from_slice(b"MCZ10S30M10S\x00");
        let fwd = make_bam_bytes_with_tlen(
            0,
            100,
            flags::PAIRED | flags::MATE_REVERSE | flags::FIRST_SEGMENT,
            b"pln",
            &[encode_op(0, 40)],
            40,
            0,
            100,
            40,
            &fwd_aux,
        );
        let mut rev_aux = Vec::new();
        rev_aux.extend_from_slice(b"MCZ40M\x00");
        let rev = make_bam_bytes_with_tlen(
            0,
            100,
            flags::PAIRED | flags::REVERSE | flags::LAST_SEGMENT,
            b"pln",
            &[encode_op(4, 10), encode_op(0, 30), encode_op(4, 10)],
            50,
            0,
            100,
            -40,
            &rev_aux,
        );
        assert_eq!(num_bases_extending_past_mate_raw(&fwd), 0, "positive strand, MC path");
        assert_eq!(
            num_bases_extending_past_mate_vs_mate_raw(&fwd, &rev),
            0,
            "positive strand, mate-record path"
        );
        assert_eq!(num_bases_extending_past_mate_raw(&rev), 10, "negative strand, MC path");
        assert_eq!(
            num_bases_extending_past_mate_vs_mate_raw(&rev, &fwd),
            10,
            "negative strand, mate-record path"
        );
    }

    /// A read lying entirely on the far side of its mate is left alone: its 3' end points away
    /// from the mate, which is an outward-facing template rather than read-through.
    ///
    /// The per-record FR test reads forward-strand orientation off `TLEN`, which is never
    /// cross-checked against the MC tag, so a record can claim a forward-oriented pair while its
    /// mate is actually mapped to its left. Forward `20M` at 301-320 against a reverse mate
    /// `10M5S` at 101-110 is such a record. Measuring it against the last shared position like
    /// any other read would put the anchor at 110, count all 20 of its bases as past it against
    /// the mate's 5 trailing soft-clipped bases, and take 15 of its 20 *aligned* bases — the same
    /// destructive answer the superseded formula gave (all 20, from the "no such position"
    /// sentinel at a coordinate before the alignment began), on a geometry that cannot be
    /// verified. Nothing may be clipped here.
    #[test]
    fn test_num_bases_extending_past_mate_raw_read_entirely_past_mate() {
        let mut aux = Vec::new();
        aux.extend_from_slice(b"MCZ10M5S\x00");
        let rec = make_bam_bytes_with_tlen(
            0,
            300, // 0-based -> 1-based 301
            flags::PAIRED | flags::MATE_REVERSE,
            b"rea",
            &[encode_op(0, 20)],
            20,
            0,
            100, // 0-based -> 1-based 101, to the *left* of this read
            100, // TLEN > 0, so the per-record FR test accepts the pair
            &aux,
        );
        assert_eq!(num_bases_extending_past_mate_raw(&rec), 0);
    }

    /// The negative-strand mirror of the case above, which only the mate-record entry point can
    /// reach: [`is_fr_pair_raw`]'s reverse arm compares the mate position recorded *on the read*
    /// against the read's own end, while [`num_bases_extending_past_mate_vs_mate_raw`] measures
    /// against the mate record in hand, and the two disagree when that field is stale.
    ///
    /// Reverse `20M` at 101-120 whose mate-position field says 101, against a forward mate `10M`
    /// actually at 301-310. The pair passes the FR test on the recorded field, and the read's 3'
    /// end (its leftmost) then points away from the mate. Anchoring at 301 would count all 20 of
    /// the read's bases as before it and clip the read away entirely.
    #[test]
    fn test_num_bases_extending_past_mate_vs_mate_raw_read_entirely_before_mate() {
        let rec = make_bam_bytes_with_tlen(
            0,
            100, // 0-based -> 1-based 101
            flags::PAIRED | flags::REVERSE,
            b"rea",
            &[encode_op(0, 20)],
            20,
            0,
            100, // stale mate position: 1-based 101, left of this read's end
            -100,
            &[],
        );
        let mate = make_bam_bytes_with_tlen(
            0,
            300, // 0-based -> 1-based 301, entirely to the right of the read
            flags::PAIRED | flags::MATE_REVERSE,
            b"rea",
            &[encode_op(0, 10)],
            10,
            0,
            100,
            100,
            &[],
        );
        assert_eq!(num_bases_extending_past_mate_vs_mate_raw(&rec, &mate), 0);
    }

    /// Builds one FR pair at 1-based `s1`/`s2` with the given CIGARs and returns the clip each
    /// read is given by the MC-tag entry point, as `(forward, reverse)`.
    fn clips_for_fr_pair(s1: i32, c1: &str, s2: i32, c2: &str) -> (usize, usize) {
        let fwd_ops = parse_mc_cigar_ops(c1).expect("valid forward CIGAR");
        let rev_ops = parse_mc_cigar_ops(c2).expect("valid reverse CIGAR");
        // TLEN spans the forward read's leftmost base to the reverse read's rightmost.
        let tlen = alignment_end_1based(s2, &rev_ops) - s1 + 1;

        let mut fwd_aux = Vec::new();
        fwd_aux.extend_from_slice(format!("MCZ{c2}\0").as_bytes());
        let fwd = make_bam_bytes_with_tlen(
            0,
            s1 - 1,
            flags::PAIRED | flags::MATE_REVERSE | flags::FIRST_SEGMENT,
            b"pair",
            &fwd_ops,
            query_length_from_cigar(&fwd_ops),
            0,
            s2 - 1,
            tlen,
            &fwd_aux,
        );

        let mut rev_aux = Vec::new();
        rev_aux.extend_from_slice(format!("MCZ{c1}\0").as_bytes());
        let rev = make_bam_bytes_with_tlen(
            0,
            s2 - 1,
            flags::PAIRED | flags::REVERSE | flags::LAST_SEGMENT,
            b"pair",
            &rev_ops,
            query_length_from_cigar(&rev_ops),
            0,
            s1 - 1,
            -tlen,
            &rev_aux,
        );

        (num_bases_extending_past_mate_raw(&fwd), num_bases_extending_past_mate_raw(&rev))
    }

    /// Both alignments lie within the insert, so a pair that read through into the adapter can
    /// still leave them disjoint. Each 100-base read here aligns only its first 20 bases and
    /// soft-clips the rest, over an insert of 39 reference bases.
    ///
    /// With the reverse read at 1019 the two alignments share exactly one reference position,
    /// 1019, past which the forward read has all 80 of its soft-clipped bases and the reverse
    /// read has 19 — so 61 bases of read-through come off each. Moving the reverse read one base
    /// to the right removes that shared position, and the answer must not jump: the insert is one
    /// base longer, so one fewer base of each read is read-through. Extrapolating the read's
    /// soft-clipped tail against the mate's soft-only unclipped boundary is what carries the
    /// count across that boundary; answering 0 there would leave 60 bases of adapter on a read
    /// whose neighbour in this table gives up 61.
    #[test]
    fn test_disjoint_alignments_still_clip_their_read_through() {
        assert_eq!(clips_for_fr_pair(1000, "20M80S", 1019, "80S20M"), (61, 61), "one shared base");
        assert_eq!(clips_for_fr_pair(1000, "20M80S", 1020, "80S20M"), (60, 60), "no shared base");
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

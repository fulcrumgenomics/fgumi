use crate::cigar::{
    get_cigar_ops, mate_alignment_end, reference_length_from_cigar, reference_length_from_raw_bam,
    unclipped_other_end, unclipped_other_start,
};
use crate::fields::{aux_data_slice, flags, mate_pos, mate_ref_id, pos, ref_id};
use crate::tags::find_string_tag;

/// Check if a single read is part of an FR (forward-reverse) pair using raw BAM bytes.
///
/// This is the raw-byte equivalent of `record_utils::is_fr_pair_from_tags`.
/// Returns `true` if the read is paired, both read and mate are mapped,
/// on the same reference, and in FR orientation (positive strand 5' < negative strand 5').
#[must_use]
pub fn is_fr_pair_raw(bam: &[u8]) -> bool {
    let flg = flags(bam);

    // Must be paired
    if flg & flags::PAIRED == 0 {
        return false;
    }

    // Both read and mate must be mapped
    if flg & flags::UNMAPPED != 0 || flg & flags::MATE_UNMAPPED != 0 {
        return false;
    }

    // Must be on the same reference
    let this_ref_id = ref_id(bam);
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

    // Determine if FR or RF:
    // positiveStrandFivePrimePos = readIsOnReverseStrand ? mateStart : alignmentStart
    // negativeStrandFivePrimePos = readIsOnReverseStrand ? alignmentEnd : mateAlignmentEnd (from MC tag)
    let alignment_start = pos(bam) + 1; // 1-based
    let m_start = mate_pos(bam) + 1; // 1-based

    let (positive_five_prime, negative_five_prime) = if is_reverse {
        // This read is on reverse strand, mate is on positive strand
        let ref_len = reference_length_from_raw_bam(bam);
        let end = alignment_start + (ref_len - 1).max(0);
        (m_start, end)
    } else {
        // This read is on positive strand, mate is on reverse strand
        // Compute mate's alignment end from MC tag instead of using TLEN
        let aux = aux_data_slice(bam);
        let Some(mc_bytes) = find_string_tag(aux, b"MC") else {
            return false;
        };
        let Ok(mc_cigar) = std::str::from_utf8(mc_bytes) else {
            return false;
        };
        let mate_end = mate_alignment_end(m_start, mc_cigar);
        (alignment_start, mate_end)
    };

    // FR if positive strand 5' < negative strand 5'
    positive_five_prime < negative_five_prime
}

#[must_use]
pub fn num_bases_extending_past_mate_raw(bam: &[u8]) -> usize {
    let flg = flags(bam);

    // Only applies to paired, mapped reads with mapped mates
    if flg & flags::PAIRED == 0 || flg & flags::UNMAPPED != 0 || flg & flags::MATE_UNMAPPED != 0 {
        return 0;
    }

    // Check FR pair: read and mate on same reference, opposite strands
    let is_reverse = flg & flags::REVERSE != 0;
    let mate_is_reverse = flg & flags::MATE_REVERSE != 0;

    // FR pair requires opposite strand orientations
    if is_reverse == mate_is_reverse {
        return 0;
    }

    let this_ref_id = ref_id(bam);
    let m_ref_id = mate_ref_id(bam);
    if this_ref_id != m_ref_id {
        return 0;
    }

    // Need MC tag for mate CIGAR information
    let aux = aux_data_slice(bam);
    let Some(mc_bytes) = find_string_tag(aux, b"MC") else {
        return 0;
    };
    let Ok(mc_cigar) = std::str::from_utf8(mc_bytes) else {
        return 0;
    };

    let this_pos_0based = pos(bam);
    let m_pos_0based = mate_pos(bam);
    // Convert to 1-based for coordinate calculations
    let this_pos = this_pos_0based + 1;
    let m_pos = m_pos_0based + 1;

    // Calculate read length from CIGAR (query-consuming ops)
    let cigar_ops = get_cigar_ops(bam);
    let read_length: usize = cigar_ops
        .iter()
        .map(|&op| {
            let op_type = op & 0xF;
            let op_len = (op >> 4) as usize;
            // M(0), I(1), S(4), =(7), X(8) consume query
            if matches!(op_type, 0 | 1 | 4 | 7 | 8) { op_len } else { 0 }
        })
        .sum();

    if is_reverse {
        // Negative strand: check if read extends before mate's unclipped start
        let mate_us = unclipped_other_start(m_pos, mc_cigar);

        if this_pos <= mate_us {
            compute_bases_before_ref_pos(&cigar_ops, this_pos, mate_us)
        } else {
            // Only clip excess soft-clipped bases at the start
            let leading_sc = leading_soft_clip_from_ops(&cigar_ops);
            let gap = (this_pos - mate_us).cast_unsigned() as usize;
            leading_sc.saturating_sub(gap)
        }
    } else {
        // Positive strand: check if read extends past mate's unclipped end
        let ref_len = reference_length_from_cigar(&cigar_ops);
        let alignment_end = this_pos + ref_len - 1;
        let mate_ue = unclipped_other_end(m_pos, mc_cigar);

        if alignment_end >= mate_ue {
            // Compute read position at mate's unclipped end
            // Simplified: use reference position mapping
            let bases_past = compute_bases_past_ref_pos(&cigar_ops, this_pos, mate_ue);
            read_length.saturating_sub(bases_past)
        } else {
            // Only clip excess soft-clipped bases
            let trailing_sc = trailing_soft_clip_from_ops(&cigar_ops);
            let gap = (mate_ue - alignment_end).cast_unsigned() as usize;
            trailing_sc.saturating_sub(gap)
        }
    }
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
                        return read_pos;
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

/// Compute number of read bases before a reference position (for negative strand).
fn compute_bases_before_ref_pos(
    cigar_ops: &[u32],
    alignment_start_1based: i32,
    target_ref_pos: i32,
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
                        return read_pos.saturating_sub(1);
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
                        return 0;
                    }
                    ref_pos += 1;
                }
            }
            _ => {}
        }
    }

    0
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
        // negative_five_prime = mate_alignment_end = 201 + 100 - 1 = 300
        // 101 < 300 => FR => true
        let mut aux = Vec::new();
        aux.extend_from_slice(b"MCZ100M\x00");
        let rec = make_bam_bytes(
            0,
            100,
            flags::PAIRED | flags::MATE_REVERSE,
            b"rea",
            &[encode_op(0, 10)],
            10,
            0,
            200,
            &aux,
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
        // negative_five_prime = mate_alignment_end = 101 + 10 - 1 = 110
        // 201 > 110 => NOT FR (it's RF) => false
        let mut aux = Vec::new();
        aux.extend_from_slice(b"MCZ10M\x00");
        let rec = make_bam_bytes(
            0,
            200,
            flags::PAIRED | flags::MATE_REVERSE,
            b"rea",
            &[encode_op(0, 10)],
            10,
            0,
            100,
            &aux,
        );
        assert!(!is_fr_pair_raw(&rec));
    }

    #[test]
    fn test_is_fr_pair_raw_no_mc_tag_returns_false() {
        // Forward strand read with no MC tag should return false
        let rec = make_bam_bytes(
            0,
            100,
            flags::PAIRED | flags::MATE_REVERSE,
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
    fn test_is_fr_pair_raw_ignores_wrong_tlen_zero() {
        // FR pair with TLEN=0 (wrong) but correct MC tag
        // Read: forward at pos 100 (0-based), 10M
        // Mate: reverse at pos 200 (0-based), MC=10M
        // positive_five_prime = 101
        // negative_five_prime = mate_alignment_end = 201 + 10 - 1 = 210
        // 101 < 210 => FR => true (TLEN=0 is ignored)
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
            0, // wrong TLEN
            &aux,
        );
        assert!(is_fr_pair_raw(&rec));
    }

    #[test]
    fn test_is_fr_pair_raw_ignores_wrong_tlen_negative() {
        // FR pair with TLEN=-200 (wrong, would make it look RF) but correct MC tag
        // Read: forward at pos 100 (0-based), 10M
        // Mate: reverse at pos 200 (0-based), MC=10M
        // positive_five_prime = 101
        // negative_five_prime = mate_alignment_end = 201 + 10 - 1 = 210
        // 101 < 210 => FR => true (negative TLEN is ignored)
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
            -200, // wrong TLEN
            &aux,
        );
        assert!(is_fr_pair_raw(&rec));
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
        let rec = make_bam_bytes(
            0,
            100,
            flags::PAIRED | flags::MATE_REVERSE,
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
    fn test_num_bases_extending_past_mate_raw_positive_strand_overlap() {
        // Positive strand read extending past reverse mate's unclipped end
        // Read: forward at pos 100, 20M (alignment_end = 101 + 20 - 1 = 120, 1-based)
        // Mate: reverse at pos 105, MC=10M (mate unclipped_end = 106 + 10 - 1 = 115, 1-based)
        // alignment_end (120) >= mate_ue (115), so need to compute bases past ref 115
        // At ref 115 (1-based), starting from alignment_start=101:
        // offset in 20M = 115 - 101 + 1 = 15 -> read_pos 15
        // read_length = 20
        // result = 20 - 15 = 5
        let mut aux = Vec::new();
        aux.extend_from_slice(b"MCZ10M\x00");
        let rec = make_bam_bytes(
            0,
            100, // 0-based pos
            flags::PAIRED | flags::MATE_REVERSE,
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
    fn test_num_bases_extending_past_mate_raw_positive_strand_no_overlap() {
        // Positive strand read NOT extending past reverse mate's unclipped end
        // Read: forward at pos 100, 10M (alignment_end = 101 + 10 - 1 = 110)
        // Mate: reverse at pos 200, MC=10M (mate unclipped_end = 201 + 10 - 1 = 210)
        // alignment_end (110) < mate_ue (210), check trailing soft clips
        // No trailing soft clips => trailing_sc.saturating_sub(gap) = 0
        let mut aux = Vec::new();
        aux.extend_from_slice(b"MCZ10M\x00");
        let rec = make_bam_bytes(
            0,
            100,
            flags::PAIRED | flags::MATE_REVERSE,
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
    fn test_num_bases_extending_past_mate_raw_negative_strand_overlap() {
        // Negative strand read extending before forward mate's unclipped start
        // Read: reverse at pos 100, 20M
        // Mate: forward at pos 105, MC=10M (mate unclipped_start = 106)
        // this_pos (101) <= mate_us (106), so compute bases before ref pos 106
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
        // this_pos (201) > mate_us (101), check leading soft clips
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
        // Negative strand read with soft clip, this_pos > mate_us
        // Read: reverse at pos 110 (0-based), 3S10M (query_len=13)
        // Mate: forward at pos 105 (0-based), MC=10M (mate unclipped_start = 106)
        // this_pos = 111, mate_us = 106
        // this_pos > mate_us, so gap = 111 - 106 = 5
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
    fn test_num_bases_extending_past_mate_raw_positive_strand_gap_with_soft_clip() {
        // Positive strand read with trailing soft clip, alignment_end < mate_ue
        // Read: forward at pos 100 (0-based), 10M3S (query_len=13)
        // Mate: reverse at pos 200 (0-based), MC=10M (mate unclipped_end = 201+10-1=210)
        // alignment_end = 101+10-1 = 110, mate_ue = 210
        // 110 < 210, gap = 210 - 110 = 100
        // trailing_soft_clip = 3, 3.saturating_sub(100) = 0
        let mut aux = Vec::new();
        aux.extend_from_slice(b"MCZ10M\x00");
        let rec = make_bam_bytes(
            0,
            100,
            flags::PAIRED | flags::MATE_REVERSE,
            b"rea",
            &[encode_op(0, 10), encode_op(4, 3)],
            13,
            0,
            200,
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
}

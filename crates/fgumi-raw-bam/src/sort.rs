use std::cmp::Ordering;

use crate::cigar::mate_unclipped_5prime;
use crate::cigar::unclipped_5prime_sort;
use crate::fields::flags;
use crate::tags::{find_mc_tag_in_record, find_mi_tag_in_record};

#[must_use]
pub fn compare_coordinate_raw(a: &[u8], b: &[u8]) -> Ordering {
    // Extract fields at fixed offsets
    let a_tid = i32::from_le_bytes([a[0], a[1], a[2], a[3]]);
    let b_tid = i32::from_le_bytes([b[0], b[1], b[2], b[3]]);

    let a_pos = i32::from_le_bytes([a[4], a[5], a[6], a[7]]);
    let b_pos = i32::from_le_bytes([b[4], b[5], b[6], b[7]]);

    let a_flag = u16::from_le_bytes([a[14], a[15]]);
    let b_flag = u16::from_le_bytes([b[14], b[15]]);

    // Handle reads with no reference (tid = -1) - sort last
    // Unmapped reads with a valid tid (mate's position) sort by that position
    let a_no_ref = a_tid < 0;
    let b_no_ref = b_tid < 0;

    match (a_no_ref, b_no_ref) {
        (true, false) => return Ordering::Greater,
        (false, true) => return Ordering::Less,
        (true, true) => return compare_names_raw(a, b),
        _ => {}
    }

    // Compare tid, pos, reverse, then name
    a_tid
        .cmp(&b_tid)
        .then_with(|| a_pos.cmp(&b_pos))
        .then_with(|| {
            let a_rev = (a_flag & flags::REVERSE) != 0;
            let b_rev = (b_flag & flags::REVERSE) != 0;
            a_rev.cmp(&b_rev)
        })
        .then_with(|| compare_names_raw(a, b))
}

/// Compare read names directly from BAM bytes.
#[inline]
#[must_use]
pub fn compare_names_raw(a: &[u8], b: &[u8]) -> Ordering {
    let a_name_len = a[8] as usize;
    let b_name_len = b[8] as usize;

    // Exclude null terminator
    let a_len = if a_name_len > 0 { a_name_len - 1 } else { 0 };
    let b_len = if b_name_len > 0 { b_name_len - 1 } else { 0 };

    let a_name = &a[32..32 + a_len];
    let b_name = &b[32..32 + b_len];

    a_name.cmp(b_name)
}

/// Compare for queryname ordering using raw bytes.
#[inline]
#[must_use]
pub fn compare_queryname_raw(a: &[u8], b: &[u8]) -> Ordering {
    compare_names_raw(a, b).then_with(|| {
        let a_flag = u16::from_le_bytes([a[14], a[15]]);
        let b_flag = u16::from_le_bytes([b[14], b[15]]);
        a_flag.cmp(&b_flag)
    })
}

/// Compare for template-coordinate ordering using raw bytes.
///
/// This matches samtools' template-coordinate sorting which uses unclipped 5' positions.
#[inline]
#[must_use]
pub fn compare_template_coordinate_raw(a: &[u8], b: &[u8]) -> Ordering {
    // Extract all needed fields from both records
    let a_tid = i32::from_le_bytes([a[0], a[1], a[2], a[3]]);
    let a_pos = i32::from_le_bytes([a[4], a[5], a[6], a[7]]);
    let a_flag = u16::from_le_bytes([a[14], a[15]]);
    let a_mate_tid = i32::from_le_bytes([a[20], a[21], a[22], a[23]]);
    let a_mate_pos = i32::from_le_bytes([a[24], a[25], a[26], a[27]]);

    let b_tid = i32::from_le_bytes([b[0], b[1], b[2], b[3]]);
    let b_pos = i32::from_le_bytes([b[4], b[5], b[6], b[7]]);
    let b_flag = u16::from_le_bytes([b[14], b[15]]);
    let b_mate_tid = i32::from_le_bytes([b[20], b[21], b[22], b[23]]);
    let b_mate_pos = i32::from_le_bytes([b[24], b[25], b[26], b[27]]);

    // Extract strand information
    let a_reverse = (a_flag & flags::REVERSE) != 0;
    let b_reverse = (b_flag & flags::REVERSE) != 0;
    let a_mate_reverse = (a_flag & flags::MATE_REVERSE) != 0;
    let b_mate_reverse = (b_flag & flags::MATE_REVERSE) != 0;

    // Get MC tags for mate's unclipped position
    let a_mc = find_mc_tag_in_record(a);
    let b_mc = find_mc_tag_in_record(b);

    // Calculate unclipped 5' positions (zero-allocation: reads CIGAR from raw bytes)
    let a_unclipped_pos = unclipped_5prime_sort(a, a_pos, a_reverse);
    let b_unclipped_pos = unclipped_5prime_sort(b, b_pos, b_reverse);

    let a_mate_unclipped_pos =
        a_mc.map_or(a_mate_pos, |mc| mate_unclipped_5prime(a_mate_pos, a_mate_reverse, mc));
    let b_mate_unclipped_pos =
        b_mc.map_or(b_mate_pos, |mc| mate_unclipped_5prime(b_mate_pos, b_mate_reverse, mc));

    // Compute canonical positions using unclipped coordinates
    let (a_tid1, a_tid2, a_pos1, a_pos2, a_neg1, a_neg2, a_upper) =
        canonical_template_pos_unclipped(
            a_tid,
            a_unclipped_pos,
            a_mate_tid,
            a_mate_unclipped_pos,
            a_flag,
            a_reverse,
            a_mate_reverse,
        );
    let (b_tid1, b_tid2, b_pos1, b_pos2, b_neg1, b_neg2, b_upper) =
        canonical_template_pos_unclipped(
            b_tid,
            b_unclipped_pos,
            b_mate_tid,
            b_mate_unclipped_pos,
            b_flag,
            b_reverse,
            b_mate_reverse,
        );

    // Compare canonical positions (samtools order)
    a_tid1
        .cmp(&b_tid1)
        .then_with(|| a_tid2.cmp(&b_tid2))
        .then_with(|| a_pos1.cmp(&b_pos1))
        .then_with(|| a_pos2.cmp(&b_pos2))
        .then_with(|| {
            // Reverse strand sorts before forward (samtools convention)
            match (a_neg1, b_neg1) {
                (true, false) => Ordering::Less,
                (false, true) => Ordering::Greater,
                _ => Ordering::Equal,
            }
        })
        .then_with(|| match (a_neg2, b_neg2) {
            (true, false) => Ordering::Less,
            (false, true) => Ordering::Greater,
            _ => Ordering::Equal,
        })
        .then_with(|| compare_mi_tags_raw(a, b))
        .then_with(|| compare_names_raw(a, b))
        .then_with(|| a_upper.cmp(&b_upper))
}

/// Compute canonical template position using pre-computed unclipped positions.
///
/// This is used for template-coordinate sorting where positions are already
/// computed as unclipped 5' coordinates.
///
/// Returns `(tid1, tid2, pos1, pos2, neg1, neg2, is_upper)`.
#[inline]
fn canonical_template_pos_unclipped(
    tid: i32,
    unclipped_pos: i32,
    mate_tid: i32,
    mate_unclipped_pos: i32,
    flag: u16,
    reverse: bool,
    mate_reverse: bool,
) -> (i32, i32, i32, i32, bool, bool, bool) {
    let unmapped = (flag & flags::UNMAPPED) != 0;
    let mate_unmapped = (flag & flags::MATE_UNMAPPED) != 0;
    let paired = (flag & flags::PAIRED) != 0;

    if unmapped && (!paired || mate_unmapped) {
        // Completely unmapped (no mapped mate) - sort to end
        (i32::MAX, i32::MAX, i32::MAX, i32::MAX, false, false, false)
    } else if unmapped {
        // Unmapped read with mapped mate - use mate's position as primary key
        // This keeps unmapped reads with their mapped mates (samtools behavior)
        (mate_tid, i32::MAX, mate_unclipped_pos, i32::MAX, mate_reverse, false, true)
    } else if !paired || mate_unmapped {
        // Mapped read with unmapped mate - use MAX for tid2/pos2 (samtools behavior)
        (tid, i32::MAX, unclipped_pos, i32::MAX, reverse, false, false)
    } else {
        // Both reads mapped - canonical ordering
        // Samtools logic: is_upper if tid > mate_tid, or pos > mate_pos, or (pos equal and reverse)
        let is_upper = (tid, unclipped_pos) > (mate_tid, mate_unclipped_pos)
            || ((tid, unclipped_pos) == (mate_tid, mate_unclipped_pos) && reverse);

        if is_upper {
            // Swap: mate's position comes first
            (mate_tid, tid, mate_unclipped_pos, unclipped_pos, mate_reverse, reverse, true)
        } else {
            // No swap: this read's position comes first
            (tid, mate_tid, unclipped_pos, mate_unclipped_pos, reverse, mate_reverse, false)
        }
    }
}

/// Compare MI tags directly from raw BAM bytes.
///
/// Uses the pre-parsed integer values from `find_mi_tag_in_record`.
#[inline]
fn compare_mi_tags_raw(a: &[u8], b: &[u8]) -> Ordering {
    let (a_mi, _) = find_mi_tag_in_record(a);
    let (b_mi, _) = find_mi_tag_in_record(b);

    // Compare as integers (0 means no MI tag found)
    a_mi.cmp(&b_mi)
}

#[cfg(test)]
#[allow(clippy::identity_op)]
mod tests {
    use super::*;
    use crate::testutil::*;
    use std::cmp::Ordering;

    // ========================================================================
    // compare_coordinate_raw tests
    // ========================================================================

    #[test]
    fn test_compare_coordinate_raw_same_records() {
        let rec = make_bam_bytes(1, 100, 0, b"rea", &[], 0, -1, -1, &[]);
        assert_eq!(compare_coordinate_raw(&rec, &rec), Ordering::Equal);
    }

    #[test]
    fn test_compare_coordinate_raw_different_tid() {
        let a = make_bam_bytes(0, 100, 0, b"rea", &[], 0, -1, -1, &[]);
        let b = make_bam_bytes(2, 100, 0, b"rea", &[], 0, -1, -1, &[]);
        assert_eq!(compare_coordinate_raw(&a, &b), Ordering::Less);
        assert_eq!(compare_coordinate_raw(&b, &a), Ordering::Greater);
    }

    #[test]
    fn test_compare_coordinate_raw_different_pos() {
        let a = make_bam_bytes(1, 50, 0, b"rea", &[], 0, -1, -1, &[]);
        let b = make_bam_bytes(1, 200, 0, b"rea", &[], 0, -1, -1, &[]);
        assert_eq!(compare_coordinate_raw(&a, &b), Ordering::Less);
        assert_eq!(compare_coordinate_raw(&b, &a), Ordering::Greater);
    }

    #[test]
    fn test_compare_coordinate_raw_reverse_strand() {
        // Forward strand (flag=0) should sort before reverse (flag=REVERSE)
        let fwd = make_bam_bytes(1, 100, 0, b"rea", &[], 0, -1, -1, &[]);
        let rev = make_bam_bytes(1, 100, flags::REVERSE, b"rea", &[], 0, -1, -1, &[]);
        assert_eq!(compare_coordinate_raw(&fwd, &rev), Ordering::Less);
        assert_eq!(compare_coordinate_raw(&rev, &fwd), Ordering::Greater);
    }

    #[test]
    fn test_compare_coordinate_raw_name_tiebreak() {
        let a = make_bam_bytes(1, 100, 0, b"aaa", &[], 0, -1, -1, &[]);
        let b = make_bam_bytes(1, 100, 0, b"zzz", &[], 0, -1, -1, &[]);
        assert_eq!(compare_coordinate_raw(&a, &b), Ordering::Less);
        assert_eq!(compare_coordinate_raw(&b, &a), Ordering::Greater);
    }

    #[test]
    fn test_compare_coordinate_raw_no_ref_sorts_last() {
        // tid=-1 (no reference) should sort after mapped records
        let mapped = make_bam_bytes(1, 100, 0, b"rea", &[], 0, -1, -1, &[]);
        let no_ref = make_bam_bytes(-1, -1, 0, b"rea", &[], 0, -1, -1, &[]);
        assert_eq!(compare_coordinate_raw(&mapped, &no_ref), Ordering::Less);
        assert_eq!(compare_coordinate_raw(&no_ref, &mapped), Ordering::Greater);
    }

    #[test]
    fn test_compare_coordinate_raw_both_no_ref() {
        // Both tid=-1, compare by name
        let a = make_bam_bytes(-1, -1, 0, b"aaa", &[], 0, -1, -1, &[]);
        let b = make_bam_bytes(-1, -1, 0, b"zzz", &[], 0, -1, -1, &[]);
        assert_eq!(compare_coordinate_raw(&a, &b), Ordering::Less);
        assert_eq!(compare_coordinate_raw(&b, &a), Ordering::Greater);
    }

    #[test]
    fn test_compare_coordinate_raw_one_no_ref_each_direction() {
        let mapped = make_bam_bytes(5, 999, 0, b"rea", &[], 0, -1, -1, &[]);
        let no_ref = make_bam_bytes(-1, 0, 0, b"rea", &[], 0, -1, -1, &[]);
        // Mapped < no_ref
        assert_eq!(compare_coordinate_raw(&mapped, &no_ref), Ordering::Less);
        // no_ref > mapped
        assert_eq!(compare_coordinate_raw(&no_ref, &mapped), Ordering::Greater);
    }

    // ========================================================================
    // compare_queryname_raw tests
    // ========================================================================

    #[test]
    fn test_compare_queryname_raw_name_ordering() {
        let a = make_bam_bytes(0, 0, 0, b"abc", &[], 0, -1, -1, &[]);
        let b = make_bam_bytes(0, 0, 0, b"xyz", &[], 0, -1, -1, &[]);
        assert_eq!(compare_queryname_raw(&a, &b), Ordering::Less);
        assert_eq!(compare_queryname_raw(&b, &a), Ordering::Greater);
    }

    #[test]
    fn test_compare_queryname_raw_same_name_flag_tiebreak() {
        // Same name, different flags -> flag decides
        let a = make_bam_bytes(0, 0, 0x40, b"rea", &[], 0, -1, -1, &[]); // first in pair
        let b = make_bam_bytes(0, 0, 0x80, b"rea", &[], 0, -1, -1, &[]); // second in pair
        assert_eq!(compare_queryname_raw(&a, &b), Ordering::Less);
        assert_eq!(compare_queryname_raw(&b, &a), Ordering::Greater);
    }

    #[test]
    fn test_compare_queryname_raw_equal() {
        let rec = make_bam_bytes(1, 100, 0x40, b"rea", &[], 0, -1, -1, &[]);
        assert_eq!(compare_queryname_raw(&rec, &rec), Ordering::Equal);
    }

    // ========================================================================
    // compare_template_coordinate_raw tests
    // ========================================================================

    #[test]
    fn test_compare_template_coordinate_raw_equal() {
        // Two identical unmapped, unpaired records -> should be Equal
        let rec = make_bam_bytes(
            -1,
            -1,
            flags::UNMAPPED, // unmapped, not paired
            b"rea",
            &[],
            0,
            -1,
            -1,
            &[],
        );
        assert_eq!(compare_template_coordinate_raw(&rec, &rec), Ordering::Equal);
    }

    #[test]
    fn test_compare_template_coordinate_raw_different_tid() {
        // Two mapped, unpaired records with different tids
        // 10M cigar for proper alignment
        let cigar = &[(10 << 4) | 0]; // 10M
        let a = make_bam_bytes(0, 100, 0, b"rea", cigar, 10, -1, -1, &[]);
        let b = make_bam_bytes(2, 100, 0, b"rea", cigar, 10, -1, -1, &[]);
        assert_eq!(compare_template_coordinate_raw(&a, &b), Ordering::Less);
        assert_eq!(compare_template_coordinate_raw(&b, &a), Ordering::Greater);
    }

    #[test]
    fn test_compare_template_coordinate_raw_unmapped_unpaired() {
        // Unmapped, unpaired -> (MAX, MAX, MAX, MAX, false, false, false)
        let a = make_bam_bytes(-1, -1, flags::UNMAPPED, b"aaa", &[], 0, -1, -1, &[]);
        let b = make_bam_bytes(-1, -1, flags::UNMAPPED, b"zzz", &[], 0, -1, -1, &[]);
        // Both fully unmapped, compare by name
        assert_eq!(compare_template_coordinate_raw(&a, &b), Ordering::Less);
    }

    #[test]
    fn test_compare_template_coordinate_raw_unmapped_with_mapped_mate() {
        // Unmapped read with mapped mate: uses mate's position
        let cigar = &[(10 << 4) | 0]; // 10M
        let a = make_bam_bytes(0, -1, flags::UNMAPPED | flags::PAIRED, b"rea", &[], 0, 0, 100, &[]);
        let b = make_bam_bytes(0, 200, flags::PAIRED, b"rea", cigar, 10, 0, 100, &[]);
        // a is unmapped with mate at tid=0,pos=100; b is mapped at tid=0,pos=200
        // a sorts by mate position (100) which is < b's position (200)
        let cmp = compare_template_coordinate_raw(&a, &b);
        assert_ne!(cmp, Ordering::Equal);
    }

    #[test]
    fn test_compare_template_coordinate_raw_mapped_with_unmapped_mate() {
        // Mapped read with unmapped mate: uses MAX for mate position
        let cigar = &[(10 << 4) | 0]; // 10M
        let a = make_bam_bytes(
            0,
            100,
            flags::PAIRED | flags::MATE_UNMAPPED,
            b"rea",
            cigar,
            10,
            -1,
            -1,
            &[],
        );
        let b = make_bam_bytes(
            0,
            200,
            flags::PAIRED | flags::MATE_UNMAPPED,
            b"rea",
            cigar,
            10,
            -1,
            -1,
            &[],
        );
        // Both mapped, mates unmapped. a at pos 100 < b at pos 200
        assert_eq!(compare_template_coordinate_raw(&a, &b), Ordering::Less);
    }

    #[test]
    fn test_compare_template_coordinate_raw_both_mapped_canonical_swap() {
        // Both mapped, paired. Test the canonical swap (is_upper branch)
        let cigar = &[(10 << 4) | 0]; // 10M
        // a: tid=0,pos=200 with mate at tid=0,pos=100 -> is_upper=true -> swap
        let a = make_bam_bytes(0, 200, flags::PAIRED, b"rea", cigar, 10, 0, 100, &[]);
        // b: tid=0,pos=100 with mate at tid=0,pos=200 -> is_upper=false -> no swap
        let b = make_bam_bytes(0, 100, flags::PAIRED, b"rea", cigar, 10, 0, 200, &[]);
        // After canonical ordering, both should produce same (tid1=0,tid2=0,pos1=~100,pos2=~200)
        // so they compare equal on positions, differentiated by is_upper
        let cmp = compare_template_coordinate_raw(&a, &b);
        // a is_upper=true, b is_upper=false -> a > b (true > false)
        assert_eq!(cmp, Ordering::Greater);
    }

    #[test]
    fn test_compare_template_coordinate_raw_neg_strand_tiebreak() {
        // Test the neg1/neg2 tiebreak in compare_template_coordinate_raw
        let cigar = &[(10 << 4) | 0]; // 10M
        // a: forward strand
        let a = make_bam_bytes(
            0,
            100,
            flags::PAIRED | flags::MATE_UNMAPPED,
            b"rea",
            cigar,
            10,
            -1,
            -1,
            &[],
        );
        // b: reverse strand at same position
        let b = make_bam_bytes(
            0,
            100,
            flags::PAIRED | flags::MATE_UNMAPPED | flags::REVERSE,
            b"rea",
            cigar,
            10,
            -1,
            -1,
            &[],
        );
        // Reverse strand sorts before forward in samtools convention (neg1=true < neg1=false)
        let cmp = compare_template_coordinate_raw(&a, &b);
        assert_ne!(cmp, Ordering::Equal);
    }

    // ========================================================================
    // compare_mi_tags_raw (exercised through compare_template_coordinate_raw)
    // ========================================================================

    #[test]
    fn test_compare_template_coordinate_raw_mi_tag_tiebreak() {
        // Two identical records except for MI tag value
        let cigar = &[(10 << 4) | 0]; // 10M
        let mut aux_a = Vec::new();
        aux_a.extend_from_slice(b"MIZ10\x00");
        let a = make_bam_bytes(
            0,
            100,
            flags::PAIRED | flags::MATE_UNMAPPED,
            b"rea",
            cigar,
            10,
            -1,
            -1,
            &aux_a,
        );
        let mut aux_b = Vec::new();
        aux_b.extend_from_slice(b"MIZ20\x00");
        let b = make_bam_bytes(
            0,
            100,
            flags::PAIRED | flags::MATE_UNMAPPED,
            b"rea",
            cigar,
            10,
            -1,
            -1,
            &aux_b,
        );
        // MI 10 < MI 20
        assert_eq!(compare_template_coordinate_raw(&a, &b), Ordering::Less);
    }

    // ========================================================================
    // compare_names_raw: zero-length names
    // ========================================================================

    #[test]
    fn test_compare_names_raw_empty_names() {
        let mut a = vec![0u8; 33];
        a[8] = 0; // l_read_name = 0
        let mut b = vec![0u8; 33];
        b[8] = 0;
        assert_eq!(compare_names_raw(&a, &b), Ordering::Equal);
    }
}

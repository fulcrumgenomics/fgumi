use std::cmp::Ordering;

use crate::fields::{RawRecordView, flags};

#[must_use]
pub fn compare_coordinate_raw(a: &[u8], b: &[u8]) -> Ordering {
    let av = RawRecordView::new(a);
    let bv = RawRecordView::new(b);

    let a_tid = av.ref_id();
    let b_tid = bv.ref_id();

    let a_pos = av.pos();
    let b_pos = bv.pos();

    let a_flag = av.flags();
    let b_flag = bv.flags();

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
    RawRecordView::new(a).read_name().cmp(RawRecordView::new(b).read_name())
}

/// Compare for queryname ordering using raw bytes.
#[inline]
#[must_use]
pub fn compare_queryname_raw(a: &[u8], b: &[u8]) -> Ordering {
    compare_names_raw(a, b)
        .then_with(|| RawRecordView::new(a).flags().cmp(&RawRecordView::new(b).flags()))
}

// ============================================================================
// Natural string comparators
// ============================================================================

/// Natural string comparison that handles numeric runs.
///
/// Compares strings such that "read1" < "read2" < "read10".
/// Uses `get_unchecked` to eliminate per-byte bounds checks in the hot loop,
/// matching the performance characteristics of samtools' `strnum_cmp`. Compares
/// numeric runs by digit count and first differing digit, avoiding integer
/// parsing entirely.
///
/// **Note:** Unlike [`natural_compare_nul`] (which matches samtools exactly),
/// this function uses `a.len().cmp(&b.len())` as the final tiebreaker, so
/// strings that differ only in leading zeros (e.g. `"01"` vs `"1"`) are *not*
/// equal — the longer string sorts after the shorter one.
///
/// **Not used in production.** The sort pipeline uses [`natural_compare_nul`]
/// (via `RawQuerynameKey`) instead. This function is exported for benchmarks
/// that compare sorting strategies.
#[must_use]
#[inline]
#[allow(unsafe_code)]
pub fn natural_compare(a: &[u8], b: &[u8]) -> Ordering {
    // Helper: read byte at index without bounds check.
    // SAFETY contract: caller must ensure `idx < slice.len()`.
    #[inline]
    unsafe fn at(s: &[u8], idx: usize) -> u8 {
        debug_assert!(idx < s.len());
        // SAFETY: caller guarantees idx < s.len().
        unsafe { *s.get_unchecked(idx) }
    }

    let alen = a.len();
    let blen = b.len();
    let mut pa: usize = 0;
    let mut pb: usize = 0;

    while pa < alen && pb < blen {
        // SAFETY: pa < alen, pb < blen (loop invariant).
        let ca = unsafe { at(a, pa) };
        let cb = unsafe { at(b, pb) };

        if !ca.is_ascii_digit() || !cb.is_ascii_digit() {
            // At least one non-digit: plain ASCII comparison (matches samtools).
            if ca != cb {
                return ca.cmp(&cb);
            }
            pa += 1;
            pb += 1;
        } else {
            // Both digits: compare numeric runs without integer parsing.
            // Skip leading zeros.
            // SAFETY: pa < alen on entry, incremented only while byte is b'0' (a digit),
            // so we stop before or at alen.
            while pa < alen && unsafe { at(a, pa) } == b'0' {
                pa += 1;
            }
            while pb < blen && unsafe { at(b, pb) } == b'0' {
                pb += 1;
            }
            // Skip matching digits.
            // SAFETY: pa < alen && pb < blen checked before each access.
            while pa < alen
                && pb < blen
                && unsafe { at(a, pa) } == unsafe { at(b, pb) }
                && unsafe { at(a, pa) }.is_ascii_digit()
            {
                pa += 1;
                pb += 1;
            }
            // Save first differing digit, then count remaining digits in both
            // runs simultaneously (interleaved, matching samtools).
            let a_is_digit = pa < alen && unsafe { at(a, pa) }.is_ascii_digit();
            let b_is_digit = pb < blen && unsafe { at(b, pb) }.is_ascii_digit();
            let diff = if a_is_digit && b_is_digit {
                // SAFETY: just confirmed pa < alen and pb < blen.
                let d = i16::from(unsafe { at(a, pa) }) - i16::from(unsafe { at(b, pb) });
                pa += 1;
                pb += 1;
                // Count remaining digits in both simultaneously
                while pa < alen
                    && pb < blen
                    && unsafe { at(a, pa) }.is_ascii_digit()
                    && unsafe { at(b, pb) }.is_ascii_digit()
                {
                    pa += 1;
                    pb += 1;
                }
                d
            } else {
                0
            };
            // Check which run is longer
            let a_more = pa < alen && unsafe { at(a, pa) }.is_ascii_digit();
            let b_more = pb < blen && unsafe { at(b, pb) }.is_ascii_digit();
            if a_more {
                return Ordering::Greater;
            } else if b_more {
                return Ordering::Less;
            } else if diff != 0 {
                return if diff > 0 { Ordering::Greater } else { Ordering::Less };
            }
        }
    }

    alen.cmp(&blen)
}

/// Natural string comparison on null-terminated byte strings using raw pointer
/// walking — no per-byte bounds checks.
///
/// This is a faithful port of samtools' `strnum_cmp` from `bam_sort.c`. It
/// compares strings such that `"read1" < "read2" < "read10"`, treating
/// consecutive digit runs as numbers (leading zeros are skipped, then digit
/// count determines magnitude, with first-differing-digit as tiebreaker).
///
/// # Safety
///
/// Both `a` and `b` must point to valid, null-terminated byte sequences.
#[must_use]
#[inline]
#[allow(unsafe_code)]
pub unsafe fn natural_compare_nul(a: *const u8, b: *const u8) -> Ordering {
    let mut pa = a;
    let mut pb = b;

    // SAFETY: caller guarantees both pointers are null-terminated, so every
    // dereference below will hit a null byte before going out of bounds.
    unsafe {
        while *pa != 0 && *pb != 0 {
            if !(*pa).is_ascii_digit() || !(*pb).is_ascii_digit() {
                if *pa != *pb {
                    return (*pa).cmp(&(*pb));
                }
                pa = pa.add(1);
                pb = pb.add(1);
            } else {
                // Skip leading zeros.
                while *pa == b'0' {
                    pa = pa.add(1);
                }
                while *pb == b'0' {
                    pb = pb.add(1);
                }
                // Skip matching significant digits.
                while (*pa).is_ascii_digit() && *pa == *pb {
                    pa = pa.add(1);
                    pb = pb.add(1);
                }
                // Record first differing digit, then count remaining digits.
                let a_digit = (*pa).is_ascii_digit();
                let b_digit = (*pb).is_ascii_digit();
                let diff = i16::from(*pa) - i16::from(*pb);
                if a_digit && b_digit {
                    pa = pa.add(1);
                    pb = pb.add(1);
                }
                while (*pa).is_ascii_digit() && (*pb).is_ascii_digit() {
                    pa = pa.add(1);
                    pb = pb.add(1);
                }
                if (*pa).is_ascii_digit() {
                    return Ordering::Greater;
                } else if (*pb).is_ascii_digit() {
                    return Ordering::Less;
                } else if diff != 0 {
                    return if diff > 0 { Ordering::Greater } else { Ordering::Less };
                }
            }
        }

        if *pa != 0 {
            Ordering::Greater
        } else if *pb != 0 {
            Ordering::Less
        } else {
            Ordering::Equal
        }
    }
}

#[cfg(test)]
#[allow(clippy::identity_op)]
mod tests {
    use rstest::rstest;

    use super::*;
    use crate::testutil::*;
    use std::cmp::Ordering;

    // ========================================================================
    // natural_compare tests
    // ========================================================================

    #[rstest]
    #[case(b"read1", b"read2", Ordering::Less)]
    #[case(b"read2", b"read10", Ordering::Less)]
    #[case(b"read10", b"read10", Ordering::Equal)]
    #[case(b"01", b"1", Ordering::Greater)] // leading zeros: longer string wins
    #[case(b"001", b"01", Ordering::Greater)]
    #[case(b"a", b"b", Ordering::Less)]
    #[case(b"", b"", Ordering::Equal)]
    #[case(b"", b"a", Ordering::Less)]
    #[case(b"9", b"a", Ordering::Less)] // digits before letters in ASCII
    #[case(b"a0b", b"a00b", Ordering::Less)] // embedded leading zeros
    #[case(b"SRR1.9", b"SRR1.10", Ordering::Less)]
    fn test_natural_compare(#[case] a: &[u8], #[case] b: &[u8], #[case] expected: Ordering) {
        assert_eq!(natural_compare(a, b), expected);
    }

    // ========================================================================
    // natural_compare_nul tests
    // ========================================================================

    #[allow(unsafe_code)]
    fn compare_nul(a: &[u8], b: &[u8]) -> Ordering {
        let mut a_nul = a.to_vec();
        a_nul.push(0);
        let mut b_nul = b.to_vec();
        b_nul.push(0);
        unsafe { natural_compare_nul(a_nul.as_ptr(), b_nul.as_ptr()) }
    }

    #[rstest]
    #[case(b"read1", b"read2", Ordering::Less)]
    #[case(b"read2", b"read10", Ordering::Less)]
    #[case(b"01", b"1", Ordering::Equal)] // nul variant: leading zeros are equal
    #[case(b"a", b"b", Ordering::Less)]
    #[case(b"", b"", Ordering::Equal)]
    #[case(b"SRR1.9", b"SRR1.10", Ordering::Less)]
    fn test_natural_compare_nul(#[case] a: &[u8], #[case] b: &[u8], #[case] expected: Ordering) {
        assert_eq!(compare_nul(a, b), expected);
    }

    // ========================================================================
    // proptest: natural_compare agrees with natural_compare_nul on NUL-free inputs
    // (except for leading-zero tiebreaker behavior)
    // ========================================================================

    proptest::proptest! {
        #[test]
        #[allow(unsafe_code)]
        fn test_natural_compare_agrees_with_nul(
            a in proptest::collection::vec(1..=255u8, 0..30),
            b in proptest::collection::vec(1..=255u8, 0..30),
        ) {
            let slice_cmp = natural_compare(&a, &b);
            let nul_cmp = compare_nul(&a, &b);
            // Both must agree on Less/Greater; they may disagree on Equal vs
            // non-Equal for leading-zero cases (natural_compare uses length
            // tiebreaker, natural_compare_nul treats them as equal).
            if slice_cmp != Ordering::Equal && nul_cmp != Ordering::Equal {
                proptest::prop_assert_eq!(slice_cmp, nul_cmp);
            }
        }
    }

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

    // ========================================================================
    // Candidate-comparator equivalence with the production oracle (S3-024)
    // ========================================================================
    //
    // `benches/core_functions.rs` holds candidate queryname comparators
    // (`compare_illumina_structured`, the samtools `strnum_cmp` port, and the
    // `natural_compare_nul` wrapper) used to benchmark sort strategies. A
    // "fast" candidate must not be promoted on speed alone without agreeing
    // with the production ordering, so this module ports the candidate
    // comparators as test-only helpers and asserts every pairwise comparison
    // over `generate_name_pairs`-style data matches `natural_compare`. (The
    // bench is `harness = false`, so a `#[cfg(test)]` module *inside* the bench
    // never runs — the equivalence oracle belongs in this crate's test module.)
    //
    // `compare_u64_hash` is intentionally excluded: it is documented as a
    // degenerate O(1) upper bound that byte-reverses ordering and collapses
    // shared-prefix names to `Equal`, so it is deliberately NOT equivalent.

    /// Faithful port of samtools' `strnum_cmp` over NUL-terminated byte slices,
    /// mirroring the bench candidate. Used here only as an equivalence reference
    /// for `natural_compare`.
    ///
    /// This is a *safe*, index-based port (no raw-pointer walk): the equivalence
    /// oracle runs only in tests, which have no hot-path perf budget, so there is
    /// no reason to introduce `unsafe` here. The caller NUL-terminates both
    /// slices, so the `!= 0` guards stop the walk at the terminator and every
    /// index stays in bounds.
    fn strnum_cmp_samtools(a: &[u8], b: &[u8]) -> i32 {
        let (mut pa, mut pb) = (0usize, 0usize);
        while a[pa] != 0 && b[pb] != 0 {
            if !a[pa].is_ascii_digit() || !b[pb].is_ascii_digit() {
                if a[pa] != b[pb] {
                    return i32::from(a[pa]) - i32::from(b[pb]);
                }
                pa += 1;
                pb += 1;
            } else {
                while a[pa] == b'0' {
                    pa += 1;
                }
                while b[pb] == b'0' {
                    pb += 1;
                }
                while a[pa].is_ascii_digit() && a[pa] == b[pb] {
                    pa += 1;
                    pb += 1;
                }
                let diff = i32::from(a[pa]) - i32::from(b[pb]);
                while a[pa].is_ascii_digit() && b[pb].is_ascii_digit() {
                    pa += 1;
                    pb += 1;
                }
                if a[pa].is_ascii_digit() {
                    return 1;
                } else if b[pb].is_ascii_digit() {
                    return -1;
                } else if diff != 0 {
                    return diff;
                }
            }
        }
        if a[pa] != 0 {
            1
        } else if b[pb] != 0 {
            -1
        } else {
            0
        }
    }

    fn compare_strnum_samtools(a: &[u8], b: &[u8]) -> Ordering {
        let mut a_nul = a.to_vec();
        a_nul.push(0);
        let mut b_nul = b.to_vec();
        b_nul.push(0);
        strnum_cmp_samtools(&a_nul, &b_nul).cmp(&0)
    }

    fn parse_int_fast(bytes: &[u8]) -> u64 {
        let mut val: u64 = 0;
        for &b in bytes {
            if b.is_ascii_digit() {
                val = val * 10 + u64::from(b - b'0');
            } else {
                break;
            }
        }
        val
    }

    /// Port of the bench's Illumina-structured comparator: compare the first
    /// four colon-delimited fields as bytes, then tile/x/y as integers.
    fn compare_illumina_structured(a: &[u8], b: &[u8]) -> Ordering {
        let mut a_fields = a.splitn(7, |&c| c == b':');
        let mut b_fields = b.splitn(7, |&c| c == b':');
        for _ in 0..4 {
            let af = a_fields.next().unwrap_or(b"");
            let bf = b_fields.next().unwrap_or(b"");
            match af.cmp(bf) {
                Ordering::Equal => {}
                ord => return ord,
            }
        }
        for _ in 0..3 {
            let af = a_fields.next().unwrap_or(b"");
            let bf = b_fields.next().unwrap_or(b"");
            match parse_int_fast(af).cmp(&parse_int_fast(bf)) {
                Ordering::Equal => {}
                ord => return ord,
            }
        }
        Ordering::Equal
    }

    /// Generate realistic sorted name pairs, mirroring the bench's
    /// `generate_name_pairs`. Numbers are non-zero-padded so the samtools-strnum
    /// and natural comparators agree (no leading-zero ambiguity).
    fn generate_name_pairs(style: &str, count: usize) -> Vec<(Vec<u8>, Vec<u8>)> {
        let mut pairs = Vec::with_capacity(count);
        match style {
            "illumina" => {
                let prefixes =
                    ["A00132:53:HFHJKDSXX", "A00132:54:HFH2JDSXX", "A00132:55:HFHKKDSXX"];
                for i in 0..count {
                    let p1 = prefixes[i % 3];
                    let p2 = prefixes[(i + 1) % 3];
                    let a = format!(
                        "{p1}:{}:{}:{}:{}",
                        (i % 4) + 1,
                        1100 + (i % 500),
                        1000 + (i * 17 % 30000),
                        1000 + (i * 31 % 50000)
                    );
                    let b = format!(
                        "{p2}:{}:{}:{}:{}",
                        ((i + 1) % 4) + 1,
                        1100 + ((i + 7) % 500),
                        1000 + ((i + 3) * 17 % 30000),
                        1000 + ((i + 5) * 31 % 50000)
                    );
                    pairs.push((a.into_bytes(), b.into_bytes()));
                }
            }
            "illumina_same_prefix" => {
                for i in 0..count {
                    let a = format!(
                        "A00132:53:HFHJKDSXX:{}:{}:{}:{}",
                        (i % 4) + 1,
                        1100 + (i % 500),
                        1000 + (i * 17 % 30000),
                        1000 + (i * 31 % 50000)
                    );
                    let b = format!(
                        "A00132:53:HFHJKDSXX:{}:{}:{}:{}",
                        ((i + 1) % 4) + 1,
                        1100 + ((i + 7) % 500),
                        1000 + ((i + 3) * 17 % 30000),
                        1000 + ((i + 5) * 31 % 50000)
                    );
                    pairs.push((a.into_bytes(), b.into_bytes()));
                }
            }
            "srr" => {
                for i in 0..count {
                    let a = format!("SRR6109273.{}", 100_000 + i * 7);
                    let b = format!("SRR6109273.{}", 100_000 + (i + 1) * 7);
                    pairs.push((a.into_bytes(), b.into_bytes()));
                }
            }
            other => panic!("unknown name style: {other}"),
        }
        pairs
    }

    fn assert_candidate_matches_oracle(style: &str, candidate: impl Fn(&[u8], &[u8]) -> Ordering) {
        for (a, b) in generate_name_pairs(style, 5_000) {
            // Forward agreement with the oracle.
            assert_eq!(
                candidate(&a, &b),
                natural_compare(&a, &b),
                "candidate disagrees with natural_compare (style={style}): {:?} vs {:?}",
                String::from_utf8_lossy(&a),
                String::from_utf8_lossy(&b),
            );
            // Reverse agreement — pins antisymmetry, so a candidate cannot match
            // the oracle on `(a, b)` while disagreeing on `(b, a)`.
            assert_eq!(
                candidate(&b, &a),
                natural_compare(&b, &a),
                "candidate disagrees with natural_compare on reversed pair (style={style}): \
                 {:?} vs {:?}",
                String::from_utf8_lossy(&b),
                String::from_utf8_lossy(&a),
            );
            // Reflexivity — every name must compare equal to itself.
            assert_eq!(
                candidate(&a, &a),
                Ordering::Equal,
                "candidate is not reflexive (style={style}): {:?}",
                String::from_utf8_lossy(&a),
            );
            assert_eq!(
                candidate(&b, &b),
                Ordering::Equal,
                "candidate is not reflexive (style={style}): {:?}",
                String::from_utf8_lossy(&b),
            );
        }
    }

    #[test]
    fn strnum_samtools_candidate_matches_natural_compare() {
        // samtools strnum_cmp agrees with natural_compare on non-zero-padded
        // numeric runs, which every generated style produces.
        for style in ["illumina", "illumina_same_prefix", "srr"] {
            assert_candidate_matches_oracle(style, compare_strnum_samtools);
        }
    }

    #[test]
    #[allow(unsafe_code)]
    fn natural_nul_candidate_matches_natural_compare() {
        // The production `natural_compare_nul` (via the `compare_nul` wrapper)
        // agrees with the slice-based `natural_compare` when names carry no
        // leading-zero ambiguity (none of the styles do).
        for style in ["illumina", "illumina_same_prefix", "srr"] {
            assert_candidate_matches_oracle(style, compare_nul);
        }
    }

    #[test]
    fn illumina_structured_candidate_matches_natural_compare() {
        // The field-parsing Illumina comparator is only a faithful model for
        // Illumina-format names, so it is checked against the oracle on the
        // Illumina styles only.
        for style in ["illumina", "illumina_same_prefix"] {
            assert_candidate_matches_oracle(style, compare_illumina_structured);
        }
    }
}

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
    // SAFETY: caller guarantees both pointers are null-terminated, and the walk
    // reads index `i` of a side only after reading every byte before it as
    // non-NUL (it advances a side solely past a byte it has just read as a
    // digit or as a non-NUL byte equal to the other side's), so every read is
    // at or before that side's terminator.
    strnum_walk(|i| unsafe { *a.add(i) }, |i| unsafe { *b.add(i) }).0
}

/// [`natural_compare_nul`] over slices: the terminator of each side is its first
/// NUL byte or, when it has none, the slice end. Returns the ordering and the
/// index each side had reached when the walk decided it — the index of the last
/// byte read, where `a.len()` (or `b.len()`) means the walk read that side's
/// virtual terminator. The indices never exceed the slice lengths.
///
/// A caller comparing *prefixes* of longer strings uses the indices to tell a
/// final answer from one that depends on bytes past a prefix: if the walk never
/// reached a truncated side's end, every byte it read is a real byte of the full
/// string and the full walk decides the same way.
#[must_use]
pub fn natural_cmp_bounded(a: &[u8], b: &[u8]) -> (Ordering, usize, usize) {
    // The byte at `i`, or the terminator past the end.
    strnum_walk(|i| a.get(i).copied().unwrap_or(0), |i| b.get(i).copied().unwrap_or(0))
}

/// samtools' `strnum_cmp` walk, shared by [`natural_compare_nul`] and
/// [`natural_cmp_bounded`] so the run sort, the merge and the owned keys apply
/// one state machine. `at_a(i)` / `at_b(i)` are each side's byte `i`, 0 at its
/// terminator; the walk never asks for an index past a side's terminator.
/// Returns the ordering and the index each side had reached.
#[allow(clippy::inline_always)] // each comparator gets its own copy of the hot walk
#[inline(always)]
fn strnum_walk(at_a: impl Fn(usize) -> u8, at_b: impl Fn(usize) -> u8) -> (Ordering, usize, usize) {
    let mut pa = 0;
    let mut pb = 0;
    while at_a(pa) != 0 && at_b(pb) != 0 {
        let (ca, cb) = (at_a(pa), at_b(pb));
        if !ca.is_ascii_digit() || !cb.is_ascii_digit() {
            if ca != cb {
                return (ca.cmp(&cb), pa, pb);
            }
            pa += 1;
            pb += 1;
        } else {
            // Skip leading zeros.
            while at_a(pa) == b'0' {
                pa += 1;
            }
            while at_b(pb) == b'0' {
                pb += 1;
            }
            // Skip matching significant digits.
            while at_a(pa).is_ascii_digit() && at_a(pa) == at_b(pb) {
                pa += 1;
                pb += 1;
            }
            // Record the first differing digit, then count the remaining digits.
            let a_digit = at_a(pa).is_ascii_digit();
            let b_digit = at_b(pb).is_ascii_digit();
            let diff = i16::from(at_a(pa)) - i16::from(at_b(pb));
            if a_digit && b_digit {
                pa += 1;
                pb += 1;
            }
            while at_a(pa).is_ascii_digit() && at_b(pb).is_ascii_digit() {
                pa += 1;
                pb += 1;
            }
            if at_a(pa).is_ascii_digit() {
                return (Ordering::Greater, pa, pb);
            } else if at_b(pb).is_ascii_digit() {
                return (Ordering::Less, pa, pb);
            } else if diff != 0 {
                return (if diff > 0 { Ordering::Greater } else { Ordering::Less }, pa, pb);
            }
        }
    }
    let ord = if at_a(pa) != 0 {
        Ordering::Greater
    } else if at_b(pb) != 0 {
        Ordering::Less
    } else {
        Ordering::Equal
    };
    (ord, pa, pb)
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
    // natural_cmp_bounded tests
    // ========================================================================

    #[rstest]
    #[case(b"read1", b"read2", Ordering::Less)]
    #[case(b"read2", b"read10", Ordering::Less)]
    #[case(b"01", b"1", Ordering::Equal)]
    #[case(b"001", b"01", Ordering::Equal)]
    #[case(b"a", b"b", Ordering::Less)]
    #[case(b"", b"", Ordering::Equal)]
    #[case(b"", b"a", Ordering::Less)]
    #[case(b"9", b"a", Ordering::Less)]
    #[case(b"a0b", b"a00b", Ordering::Equal)]
    #[case(b"SRR1.9", b"SRR1.10", Ordering::Less)]
    #[case(b"ab\0zz", b"ab\0aa", Ordering::Equal)] // the first NUL terminates
    #[case(b"pab\0", b"pab", Ordering::Equal)]
    fn test_natural_cmp_bounded(#[case] a: &[u8], #[case] b: &[u8], #[case] expected: Ordering) {
        assert_eq!(natural_cmp_bounded(a, b).0, expected);
        assert_eq!(natural_cmp_bounded(b, a).0, expected.reverse());
    }

    /// The indices name the last byte read: a digit run decided at the window
    /// edge reads the virtual terminator of the shorter side.
    #[rstest]
    #[case(b"x:12", b"x:13", (Ordering::Less, 4, 4))]
    #[case(b"x:3333", b"x:3334", (Ordering::Less, 6, 6))]
    #[case(b"pab", b"pac", (Ordering::Less, 2, 2))]
    #[case(b"pab", b"pab", (Ordering::Equal, 3, 3))]
    #[case(b"pab", b"pabc", (Ordering::Less, 3, 3))]
    fn natural_cmp_bounded_reports_the_indices_it_read(
        #[case] a: &[u8],
        #[case] b: &[u8],
        #[case] expected: (Ordering, usize, usize),
    ) {
        assert_eq!(natural_cmp_bounded(a, b), expected);
    }

    proptest::proptest! {
        /// Equal to `natural_compare_nul` over NUL-appended copies, for arbitrary
        /// bytes including embedded NULs; the indices stay within the slices.
        #[test]
        fn natural_cmp_bounded_agrees_with_natural_compare_nul(
            a in proptest::collection::vec(proptest::prop_oneof![
                8 => proptest::sample::select(b"0123456789".to_vec()),
                4 => proptest::sample::select(b"ab:._".to_vec()),
                1 => proptest::prelude::any::<u8>(),
            ], 0..30),
            b in proptest::collection::vec(proptest::prop_oneof![
                8 => proptest::sample::select(b"0123456789".to_vec()),
                4 => proptest::sample::select(b"ab:._".to_vec()),
                1 => proptest::prelude::any::<u8>(),
            ], 0..30),
        ) {
            let (ord, ia, ib) = natural_cmp_bounded(&a, &b);
            proptest::prop_assert_eq!(ord, compare_nul(&a, &b));
            proptest::prop_assert!(ia <= a.len() && ib <= b.len());
        }
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
}

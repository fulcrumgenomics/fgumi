//! Template-coordinate sorting for simulated BAM files.
//!
//! This module provides the sort key structures matching the template-coordinate
//! sort order used by `samtools sort --template-coordinate`.
//!
//! Reference implementation:
//! - samtools: `bam_sort.c` - `template_coordinate_key_t` and `bam1_cmp_template_coordinate`
//!
//! Note: fgbio's `SamOrder.scala` `TemplateCoordinateKey` has a different strand comparison
//! (forward before reverse) vs samtools (reverse before forward). We match samtools since
//! that's what `fgumi group` validates against.

use std::cmp::Ordering;

/// Template-coordinate sort key matching the sort order used by `fgumi group` and
/// `samtools sort --template-coordinate`.
///
/// Sort order (matching samtools `bam_sort.c` `bam1_cmp_template_coordinate`):
/// 1. tid1 - reference ID of earlier read
/// 2. tid2 - reference ID of later read
/// 3. pos1 - unclipped 5' position of earlier read
/// 4. pos2 - unclipped 5' position of later read
/// 5. neg1 - strand of earlier read (reverse before forward in samtools)
/// 6. neg2 - strand of later read
/// 7. library (omitted - constant for simulated data)
/// 8. mid - MI tag without /A /B suffix
/// 9. name - read name (or molecule name prefix for sorting)
/// 10. `is_upper_of_pair`
///
/// For unclipped 5' positions:
/// - Forward strand: `unclipped_start` (leftmost position)
/// - Reverse strand: `unclipped_end` (rightmost position, which is the 5' end)
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct TemplateCoordKey {
    /// Reference ID of earlier read (Int32.MAX for unmapped)
    pub tid1: i32,
    /// Reference ID of later read (Int32.MAX for unmapped)
    pub tid2: i32,
    /// Unclipped 5' position of earlier read
    pub pos1: i64,
    /// Unclipped 5' position of later read
    pub pos2: i64,
    /// Strand of earlier read (true = reverse/negative)
    pub neg1: bool,
    /// Strand of later read (true = reverse/negative)
    pub neg2: bool,
    /// Molecular identifier (MI tag value, without /A /B suffix)
    pub mid: String,
    /// Read name (or molecule name prefix for molecule-level sorting)
    pub name: String,
    /// True if this read is the "upper" (later) of the pair
    pub is_upper_of_pair: bool,
}

impl TemplateCoordKey {
    /// Create a sort key for an F1R2 read pair (R1 forward, R2 reverse).
    ///
    /// This is the standard orientation for Illumina paired-end sequencing.
    ///
    /// # Arguments
    /// * `tid` - Reference sequence ID (chromosome index)
    /// * `pos1` - R1 alignment start position (0-based)
    /// * `insert_size` - Insert size (template length)
    /// * `mid` - Molecular identifier (empty string if none)
    /// * `name` - Read/molecule name for tie-breaking
    pub fn for_f1r2_pair(
        tid: i32,
        pos1: usize,
        insert_size: usize,
        mid: String,
        name: String,
    ) -> Self {
        // For F1R2 pairs:
        // - R1 (forward) at pos1: unclipped 5' = pos1 (unclipped_start)
        // - R2 (reverse): unclipped 5' = unclipped_end = pos1 + insert_size - 1
        let pos2 = (pos1 + insert_size).saturating_sub(1);

        Self {
            tid1: tid,
            tid2: tid,
            pos1: pos1 as i64,
            pos2: pos2 as i64,
            neg1: false, // R1 is forward
            neg2: true,  // R2 is reverse
            mid,
            name,
            is_upper_of_pair: false, // R1 is always "earlier" for F1R2
        }
    }
}

/// Compare MI tags using the samtools/fgbio algorithm.
///
/// The comparison:
/// 1. Strips trailing "/X" suffix (e.g., "/A", "/B") if present
/// 2. Compares by length first (shorter sorts before longer)
/// 3. Then compares lexicographically
fn compare_mid(mid1: &str, mid2: &str) -> Ordering {
    // Strip trailing "/X" suffix if present
    let mid1 = if mid1.len() >= 2 && mid1.as_bytes()[mid1.len() - 2] == b'/' {
        &mid1[..mid1.len() - 2]
    } else {
        mid1
    };
    let mid2 = if mid2.len() >= 2 && mid2.as_bytes()[mid2.len() - 2] == b'/' {
        &mid2[..mid2.len() - 2]
    } else {
        mid2
    };

    // Compare by length first, then lexicographically
    mid1.len().cmp(&mid2.len()).then_with(|| mid1.cmp(mid2))
}

impl Ord for TemplateCoordKey {
    fn cmp(&self, other: &Self) -> Ordering {
        // Match samtools bam1_cmp_template_coordinate order exactly
        self.tid1
            .cmp(&other.tid1)
            .then_with(|| self.tid2.cmp(&other.tid2))
            .then_with(|| self.pos1.cmp(&other.pos1))
            .then_with(|| self.pos2.cmp(&other.pos2))
            // samtools: neg=true sorts before neg=false (reverse before forward)
            .then_with(|| match (self.neg1, other.neg1) {
                (true, false) => Ordering::Less,
                (false, true) => Ordering::Greater,
                _ => Ordering::Equal,
            })
            .then_with(|| match (self.neg2, other.neg2) {
                (true, false) => Ordering::Less,
                (false, true) => Ordering::Greater,
                _ => Ordering::Equal,
            })
            // library comparison omitted - constant for simulated data
            .then_with(|| compare_mid(&self.mid, &other.mid))
            .then_with(|| self.name.cmp(&other.name))
            .then_with(|| self.is_upper_of_pair.cmp(&other.is_upper_of_pair))
    }
}

impl PartialOrd for TemplateCoordKey {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper to create a key with all fields specified.
    #[allow(clippy::too_many_arguments)]
    fn key(
        tid1: i32,
        tid2: i32,
        pos1: i64,
        pos2: i64,
        neg1: bool,
        neg2: bool,
        mid: &str,
        name: &str,
        is_upper: bool,
    ) -> TemplateCoordKey {
        TemplateCoordKey {
            tid1,
            tid2,
            pos1,
            pos2,
            neg1,
            neg2,
            mid: mid.to_string(),
            name: name.to_string(),
            is_upper_of_pair: is_upper,
        }
    }

    // =========================================================================
    // Basic for_f1r2_pair tests
    // =========================================================================

    #[test]
    fn test_f1r2_pair_key() {
        let key =
            TemplateCoordKey::for_f1r2_pair(0, 1000, 300, String::new(), "mol001".to_string());
        assert_eq!(key.tid1, 0);
        assert_eq!(key.tid2, 0);
        assert_eq!(key.pos1, 1000);
        assert_eq!(key.pos2, 1299); // 1000 + 300 - 1
        assert!(!key.neg1);
        assert!(key.neg2);
        assert!(!key.is_upper_of_pair);
    }

    #[test]
    fn test_f1r2_pair_key_with_mi() {
        let key =
            TemplateCoordKey::for_f1r2_pair(1, 500, 250, "42/A".to_string(), "read1".to_string());
        assert_eq!(key.tid1, 1);
        assert_eq!(key.tid2, 1);
        assert_eq!(key.pos1, 500);
        assert_eq!(key.pos2, 749); // 500 + 250 - 1
        assert_eq!(key.mid, "42/A");
        assert_eq!(key.name, "read1");
    }

    #[test]
    fn test_f1r2_pair_key_zero_insert() {
        // Edge case: insert_size = 1 (minimum valid)
        let key = TemplateCoordKey::for_f1r2_pair(0, 100, 1, String::new(), "x".to_string());
        assert_eq!(key.pos1, 100);
        assert_eq!(key.pos2, 100); // 100 + 1 - 1 = 100
    }

    // =========================================================================
    // Sort order tests - tid (reference ID)
    // =========================================================================

    #[test]
    fn test_sort_by_tid1() {
        let key1 = key(0, 0, 1000, 1100, false, true, "", "a", false);
        let key2 = key(1, 0, 500, 600, false, true, "", "a", false);
        // tid1=0 sorts before tid1=1, regardless of position
        assert!(key1 < key2);
    }

    #[test]
    fn test_sort_by_tid2_when_tid1_equal() {
        let key1 = key(0, 0, 1000, 1100, false, true, "", "a", false);
        let key2 = key(0, 1, 1000, 1100, false, true, "", "a", false);
        // tid2=0 sorts before tid2=1
        assert!(key1 < key2);
    }

    // =========================================================================
    // Sort order tests - position
    // =========================================================================

    #[test]
    fn test_sort_by_pos1() {
        let key1 = TemplateCoordKey::for_f1r2_pair(0, 1000, 300, String::new(), "a".to_string());
        let key2 = TemplateCoordKey::for_f1r2_pair(0, 2000, 300, String::new(), "b".to_string());
        assert!(key1 < key2);
    }

    #[test]
    fn test_sort_by_pos2_when_pos1_equal() {
        let key1 = TemplateCoordKey::for_f1r2_pair(0, 1000, 200, String::new(), "a".to_string());
        let key2 = TemplateCoordKey::for_f1r2_pair(0, 1000, 300, String::new(), "a".to_string());
        // key1.pos2 = 1199, key2.pos2 = 1299
        assert!(key1 < key2);
    }

    #[test]
    fn test_sort_by_pos1_negative_positions() {
        // Some coordinate systems use negative positions (though rare)
        let key1 = key(0, 0, -100, 0, false, true, "", "a", false);
        let key2 = key(0, 0, 0, 100, false, true, "", "a", false);
        assert!(key1 < key2);
    }

    // =========================================================================
    // Sort order tests - strand (neg1, neg2)
    // samtools: reverse (neg=true) sorts BEFORE forward (neg=false)
    // =========================================================================

    #[test]
    fn test_sort_by_neg1_reverse_before_forward() {
        // samtools convention: reverse strand sorts before forward
        let key_rev = key(0, 0, 100, 200, true, false, "", "a", false);
        let key_fwd = key(0, 0, 100, 200, false, false, "", "a", false);
        assert!(
            key_rev < key_fwd,
            "reverse strand (neg1=true) should sort before forward (neg1=false)"
        );
    }

    #[test]
    fn test_sort_by_neg2_reverse_before_forward() {
        let key_rev = key(0, 0, 100, 200, false, true, "", "a", false);
        let key_fwd = key(0, 0, 100, 200, false, false, "", "a", false);
        assert!(
            key_rev < key_fwd,
            "reverse strand (neg2=true) should sort before forward (neg2=false)"
        );
    }

    #[test]
    fn test_sort_neg1_takes_precedence_over_neg2() {
        // When neg1 differs, it takes precedence over neg2
        let key1 = key(0, 0, 100, 200, true, false, "", "a", false); // neg1=true
        let key2 = key(0, 0, 100, 200, false, true, "", "a", false); // neg1=false
        assert!(key1 < key2, "neg1=true should sort before neg1=false regardless of neg2");
    }

    // =========================================================================
    // Sort order tests - MI tag (compare_mid)
    // =========================================================================

    #[test]
    fn test_compare_mid_strips_suffix() {
        assert_eq!(compare_mid("123/A", "123/B"), Ordering::Equal);
        assert_eq!(compare_mid("123/A", "123"), Ordering::Equal);
        assert_eq!(compare_mid("123/B", "123"), Ordering::Equal);
    }

    #[test]
    fn test_compare_mid_strips_various_suffixes() {
        // Any single character after / should be stripped
        assert_eq!(compare_mid("42/A", "42/B"), Ordering::Equal);
        assert_eq!(compare_mid("42/X", "42/Y"), Ordering::Equal);
        assert_eq!(compare_mid("42/1", "42/2"), Ordering::Equal);
    }

    #[test]
    fn test_compare_mid_length_first() {
        // Shorter MI sorts before longer
        assert_eq!(compare_mid("1", "12"), Ordering::Less);
        assert_eq!(compare_mid("99", "100"), Ordering::Less);
        assert_eq!(compare_mid("9", "10"), Ordering::Less);
    }

    #[test]
    fn test_compare_mid_lexicographic() {
        assert_eq!(compare_mid("abc", "abd"), Ordering::Less);
        assert_eq!(compare_mid("123", "124"), Ordering::Less);
        assert_eq!(compare_mid("aaa", "aab"), Ordering::Less);
    }

    #[test]
    fn test_compare_mid_empty_strings() {
        assert_eq!(compare_mid("", ""), Ordering::Equal);
        assert_eq!(compare_mid("", "1"), Ordering::Less);
        assert_eq!(compare_mid("1", ""), Ordering::Greater);
    }

    #[test]
    fn test_compare_mid_single_char() {
        // Single char MIs should compare lexicographically
        assert_eq!(compare_mid("1", "2"), Ordering::Less);
        assert_eq!(compare_mid("a", "b"), Ordering::Less);
        assert_eq!(compare_mid("0", "9"), Ordering::Less);
    }

    #[test]
    fn test_compare_mid_with_suffix_vs_without() {
        // "1/A" stripped = "1", "1" = "1" → Equal
        assert_eq!(compare_mid("1/A", "1"), Ordering::Equal);
        // "12/B" stripped = "12", vs "1" → "12" > "1" by length
        assert_eq!(compare_mid("12/B", "1"), Ordering::Greater);
    }

    #[test]
    fn test_compare_mid_numeric_ordering() {
        // Numeric MIs: "1" < "2" < "9" < "10" < "11" < "99" < "100"
        // This is length-first, then lexicographic
        let mids = ["1", "2", "9", "10", "11", "99", "100", "101"];
        for i in 0..mids.len() - 1 {
            assert_eq!(
                compare_mid(mids[i], mids[i + 1]),
                Ordering::Less,
                "{} should sort before {}",
                mids[i],
                mids[i + 1]
            );
        }
    }

    #[test]
    fn test_sort_by_mid() {
        let mut key1 =
            TemplateCoordKey::for_f1r2_pair(0, 1000, 300, "1".to_string(), "a".to_string());
        let mut key2 =
            TemplateCoordKey::for_f1r2_pair(0, 1000, 300, "12".to_string(), "a".to_string());
        // Shorter MI sorts first
        assert!(key1 < key2);

        // Same length, lexicographic
        key1.mid = "12".to_string();
        key2.mid = "13".to_string();
        assert!(key1 < key2);
    }

    // =========================================================================
    // Sort order tests - name
    // =========================================================================

    #[test]
    fn test_sort_by_name_when_positions_equal() {
        let key1 =
            TemplateCoordKey::for_f1r2_pair(0, 1000, 300, String::new(), "mol001".to_string());
        let key2 =
            TemplateCoordKey::for_f1r2_pair(0, 1000, 300, String::new(), "mol002".to_string());
        assert!(key1 < key2);
    }

    #[test]
    fn test_sort_by_name_lexicographic() {
        let key1 = key(0, 0, 100, 200, false, true, "1", "aaa", false);
        let key2 = key(0, 0, 100, 200, false, true, "1", "aab", false);
        assert!(key1 < key2);
    }

    #[test]
    fn test_sort_by_name_empty_vs_nonempty() {
        let key1 = key(0, 0, 100, 200, false, true, "", "", false);
        let key2 = key(0, 0, 100, 200, false, true, "", "a", false);
        assert!(key1 < key2);
    }

    // =========================================================================
    // Sort order tests - is_upper_of_pair
    // =========================================================================

    #[test]
    fn test_sort_by_is_upper_of_pair() {
        let key_lower = key(0, 0, 100, 200, false, true, "1", "read", false);
        let key_upper = key(0, 0, 100, 200, false, true, "1", "read", true);
        assert!(key_lower < key_upper, "is_upper_of_pair=false should sort before true");
    }

    // =========================================================================
    // Sort order tests - comprehensive ordering (based on fgbio tests)
    // =========================================================================

    /// Test based on fgbio's "sort by molecular identifier then name" test.
    /// Adapted for samtools strand ordering (reverse before forward).
    #[test]
    fn test_sort_by_mi_then_name_fgbio_style() {
        // Create keys similar to fgbio test but with samtools strand ordering
        // fgbio test has pairs at pos 100 and 200, various MIs
        let mut keys = [
            // MI="0/A" at pos 200 (A strand: F1R2)
            key(0, 0, 200, 209, false, true, "0/A", "ab0", false),
            // MI="1/A" at pos 100
            key(0, 0, 100, 109, false, true, "1/A", "ab1", false),
            key(0, 0, 100, 109, false, true, "1/A", "ab2", false),
            // MI="2/A" at pos 100
            key(0, 0, 100, 109, false, true, "2/A", "ab3", false),
            // MI="0/B" at pos 200 (B strand: R1F2, so neg1=true, neg2=false)
            key(0, 0, 200, 209, true, false, "0/B", "ba0", false),
            // MI="1/B" at pos 100
            key(0, 0, 100, 109, true, false, "1/B", "ba1", false),
            key(0, 0, 100, 109, true, false, "1/B", "ba2", false),
            // MI="2/B" at pos 100
            key(0, 0, 100, 109, true, false, "2/B", "ba3", false),
        ];

        keys.sort();

        // Expected order (pos 100 before 200, then by strand, then by MI, then by name):
        // pos=100, neg1=true (B strand) before neg1=false (A strand)
        // Within same strand: MI sorted by length then lex ("1" < "2")
        // Within same MI: name sorted lex

        let names: Vec<&str> = keys.iter().map(|k| k.name.as_str()).collect();

        // At pos 100: B strand first (neg1=true), then A strand (neg1=false)
        // B strand at pos 100: MI "1/B" (ba1, ba2), then MI "2/B" (ba3)
        // A strand at pos 100: MI "1/A" (ab1, ab2), then MI "2/A" (ab3)
        // Then pos 200: B strand (ba0), A strand (ab0)
        let expected = vec!["ba1", "ba2", "ba3", "ab1", "ab2", "ab3", "ba0", "ab0"];

        assert_eq!(names, expected);
    }

    /// Test that position takes precedence over everything else (except tid).
    #[test]
    fn test_position_takes_precedence() {
        let key_pos100 = key(0, 0, 100, 200, false, true, "999", "zzz", false);
        let key_pos200 = key(0, 0, 200, 300, true, false, "1", "aaa", false);
        assert!(
            key_pos100 < key_pos200,
            "lower position should sort first regardless of other fields"
        );
    }

    /// Test that tid takes precedence over position.
    #[test]
    fn test_tid_takes_precedence_over_position() {
        let key_tid0 = key(0, 0, 99999, 100_000, false, true, "", "a", false);
        let key_tid1 = key(1, 1, 1, 2, false, true, "", "a", false);
        assert!(key_tid0 < key_tid1, "tid=0 should sort before tid=1 regardless of position");
    }

    /// Test complete sort order hierarchy.
    #[test]
    fn test_complete_sort_order_hierarchy() {
        // Create keys that differ only in the field being tested
        let baseline = key(0, 0, 100, 200, false, true, "1", "name", false);

        // tid1 difference
        let diff_tid1 = key(1, 0, 100, 200, false, true, "1", "name", false);
        assert!(baseline < diff_tid1);

        // tid2 difference (tid1 same)
        let diff_tid2 = key(0, 1, 100, 200, false, true, "1", "name", false);
        assert!(baseline < diff_tid2);

        // pos1 difference
        let diff_pos1 = key(0, 0, 200, 200, false, true, "1", "name", false);
        assert!(baseline < diff_pos1);

        // pos2 difference (pos1 same)
        let diff_pos2 = key(0, 0, 100, 300, false, true, "1", "name", false);
        assert!(baseline < diff_pos2);

        // neg1 difference (reverse before forward)
        let baseline_rev = key(0, 0, 100, 200, true, true, "1", "name", false);
        assert!(baseline_rev < baseline);

        // neg2 difference
        let diff_neg2 = key(0, 0, 100, 200, false, false, "1", "name", false);
        assert!(baseline < diff_neg2);

        // mid difference
        let diff_mid = key(0, 0, 100, 200, false, true, "2", "name", false);
        assert!(baseline < diff_mid);

        // name difference
        let diff_name = key(0, 0, 100, 200, false, true, "1", "zzzz", false);
        assert!(baseline < diff_name);

        // is_upper difference
        let diff_upper = key(0, 0, 100, 200, false, true, "1", "name", true);
        assert!(baseline < diff_upper);
    }

    // =========================================================================
    // Sorting a collection of keys
    // =========================================================================

    #[test]
    fn test_sort_collection() {
        let mut keys = [
            key(1, 1, 100, 200, false, true, "2", "b", false),
            key(0, 0, 100, 200, false, true, "1", "a", false),
            key(0, 0, 50, 150, false, true, "1", "a", false),
            key(0, 0, 100, 200, true, false, "1", "a", false), // reverse strand
            key(0, 0, 100, 200, false, true, "1", "b", false),
        ];

        keys.sort();

        // Expected order:
        // 1. tid=0, pos1=50
        // 2. tid=0, pos1=100, neg1=true (reverse)
        // 3. tid=0, pos1=100, neg1=false, name="a"
        // 4. tid=0, pos1=100, neg1=false, name="b"
        // 5. tid=1
        assert_eq!(keys[0].pos1, 50);
        assert!(keys[1].neg1);
        assert_eq!(keys[2].name, "a");
        assert_eq!(keys[3].name, "b");
        assert_eq!(keys[4].tid1, 1);
    }

    // =========================================================================
    // Edge cases
    // =========================================================================

    #[test]
    fn test_equal_keys() {
        let key1 = key(0, 0, 100, 200, false, true, "1", "name", false);
        let key2 = key(0, 0, 100, 200, false, true, "1", "name", false);
        assert_eq!(key1.cmp(&key2), Ordering::Equal);
        assert_eq!(key1, key2);
    }

    #[test]
    fn test_unmapped_reads_use_max_tid() {
        // Convention: unmapped reads use i32::MAX for tid
        let mapped = key(0, 0, 100, 200, false, true, "", "a", false);
        let unmapped = key(i32::MAX, i32::MAX, 0, 0, false, false, "", "a", false);
        assert!(mapped < unmapped, "mapped reads should sort before unmapped");
    }

    #[test]
    fn test_compare_mid_does_not_strip_short_strings() {
        // Strings shorter than 2 chars should not be modified
        assert_eq!(compare_mid("", ""), Ordering::Equal);
        assert_eq!(compare_mid("a", "a"), Ordering::Equal);
        assert_eq!(compare_mid("/", "/"), Ordering::Equal);
    }

    #[test]
    fn test_compare_mid_requires_slash_at_correct_position() {
        // "/A" should strip to "" (the slash is at position len-2)
        // "A/" should NOT strip (slash not at len-2)
        assert_eq!(compare_mid("/A", ""), Ordering::Equal);
        assert_eq!(compare_mid("A/", "A/"), Ordering::Equal); // No stripping
    }
}

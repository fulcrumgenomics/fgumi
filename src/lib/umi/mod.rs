//! UMI (Unique Molecular Identifier) utilities
//!
//! This module provides functionality for working with UMIs including:
//! - UMI assignment strategies (identity, edit-distance, adjacency, paired)
//! - Tag set management for consensus calling
//! - UMI grouping and family analysis

pub mod assigner;

// Re-export commonly used items
pub use assigner::{
    AdjacencyUmiAssigner, IdentityUmiAssigner, PairedUmiAssigner, SimpleErrorUmiAssigner, Strategy,
    Umi, UmiAssigner,
};
// MoleculeId is now defined in template module
pub use crate::template::MoleculeId;

use std::collections::HashSet;

/// Named sets of tags for common use cases
///
/// Provides predefined tag sets that can be referenced by name (e.g., "Consensus")
/// instead of listing individual tags.
pub struct TagSets;

impl TagSets {
    /// Returns the list of Consensus tags that should be reversed
    ///
    /// These are per-base tags from fgbio consensus callers that need to be
    /// reversed when the read is mapped to the negative strand.
    ///
    /// # Returns
    ///
    /// Vector of tag names: `["ad", "ae", "bd", "be", "cd"]`
    #[must_use]
    pub fn consensus_reverse() -> Vec<String> {
        vec![
            "ad".to_string(),
            "ae".to_string(),
            "bd".to_string(),
            "be".to_string(),
            "cd".to_string(),
        ]
    }

    /// Returns the list of Consensus tags that should be reverse complemented
    ///
    /// These are per-base sequence tags from fgbio consensus callers that need
    /// to be reverse complemented when the read is mapped to the negative strand.
    ///
    /// # Returns
    ///
    /// Vector of tag names: `["aD", "bD", "cD"]`
    #[must_use]
    pub fn consensus_revcomp() -> Vec<String> {
        vec!["aD".to_string(), "bD".to_string(), "cD".to_string()]
    }
}

/// Information about which tags to remove, reverse, or reverse complement
///
/// Stores sets of tag names that should be manipulated during the merge operation.
/// Named tag sets (like "Consensus") are expanded into their constituent tags.
#[derive(Debug, Clone)]
pub struct TagInfo {
    /// Tags to remove from mapped reads
    pub remove: HashSet<String>,
    /// Tags to reverse for negative strand reads
    pub reverse: HashSet<String>,
    /// Tags to reverse complement for negative strand reads
    pub revcomp: HashSet<String>,
}

impl TagInfo {
    /// Creates a new `TagInfo` from lists of tag names
    ///
    /// Expands named tag sets (e.g., "Consensus") into their constituent tags.
    /// Individual tag names are added as-is.
    ///
    /// # Arguments
    ///
    /// * `remove` - Tags to remove from mapped reads
    /// * `reverse` - Tags to reverse for negative strand (or "Consensus")
    /// * `revcomp` - Tags to reverse complement for negative strand (or "Consensus")
    ///
    /// # Examples
    ///
    /// ```
    /// use fgumi_lib::umi::TagInfo;
    /// let tag_info = TagInfo::new(
    ///     vec!["AS".to_string()],
    ///     vec!["Consensus".to_string()],
    ///     vec!["Consensus".to_string()],
    /// );
    /// ```
    #[must_use]
    pub fn new(remove: Vec<String>, reverse: Vec<String>, revcomp: Vec<String>) -> Self {
        let mut reverse_set = HashSet::new();
        let mut revcomp_set = HashSet::new();

        for tag in reverse {
            if tag == "Consensus" {
                reverse_set.extend(TagSets::consensus_reverse());
            } else {
                reverse_set.insert(tag);
            }
        }

        for tag in revcomp {
            if tag == "Consensus" {
                revcomp_set.extend(TagSets::consensus_revcomp());
            } else {
                revcomp_set.insert(tag);
            }
        }

        TagInfo { remove: remove.into_iter().collect(), reverse: reverse_set, revcomp: revcomp_set }
    }

    /// Checks if any tags need reversal or reverse complementation
    ///
    /// # Returns
    ///
    /// `true` if there are any tags to reverse or revcomp, `false` otherwise
    #[must_use]
    pub fn has_revs_or_revcomps(&self) -> bool {
        !self.reverse.is_empty() || !self.revcomp.is_empty()
    }
}

/// Result of validating a UMI string.
///
/// Used by both `group` and `dedup` commands to consistently filter UMIs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UmiValidation {
    /// UMI is valid with the given number of bases (ACGT, case-insensitive)
    Valid(usize),
    /// UMI contains 'N' (uppercase only, matching fgbio behavior)
    ContainsN,
}

/// Validates a UMI string, matching fgbio's behavior.
///
/// This function:
/// - Counts valid DNA bases (A, C, G, T - case insensitive)
/// - Rejects UMIs containing uppercase 'N' (ambiguous base)
/// - Skips dashes (paired UMI separator) and other characters including lowercase 'n'
///
/// This matches fgbio's `GroupReadsByUmi` which:
/// 1. Only rejects uppercase 'N': `.filter(r => !umi.contains('N'))`
/// 2. Counts bases case-insensitively: `umi.toUpperCase.count(c => isUpperACGTN(c))`
///
/// # Arguments
///
/// * `umi` - The UMI string to validate
///
/// # Returns
///
/// * `UmiValidation::Valid(count)` - UMI is valid with `count` DNA bases
/// * `UmiValidation::ContainsN` - UMI contains uppercase 'N'
///
/// # Examples
///
/// ```
/// use fgumi_lib::umi::{validate_umi, UmiValidation};
///
/// assert_eq!(validate_umi(b"ACGT"), UmiValidation::Valid(4));
/// assert_eq!(validate_umi(b"acgt"), UmiValidation::Valid(4));
/// assert_eq!(validate_umi(b"ACGT-TGCA"), UmiValidation::Valid(8));
/// assert_eq!(validate_umi(b"ACNT"), UmiValidation::ContainsN);
/// assert_eq!(validate_umi(b"acnt"), UmiValidation::Valid(3)); // lowercase 'n' skipped
/// ```
#[must_use]
pub fn validate_umi(umi: &[u8]) -> UmiValidation {
    let mut base_count = 0usize;
    for &b in umi {
        match b {
            b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't' => base_count += 1,
            b'N' => return UmiValidation::ContainsN,
            _ => {} // Skip dash, lowercase 'n', other chars (matches fgbio)
        }
    }
    UmiValidation::Valid(base_count)
}

/// Extracts the molecular identifier base, removing any trailing /A, /B suffix.
///
/// MI tags in duplex sequencing have suffixes `/A` or `/B` to indicate
/// which strand the read came from. This function strips those specific suffixes
/// to get the base molecular identifier for grouping purposes.
///
/// Only `/A` and `/B` suffixes are stripped; other suffixes are preserved.
///
/// # Examples
///
/// ```
/// use fgumi_lib::umi::extract_mi_base;
///
/// assert_eq!(extract_mi_base("123"), "123");
/// assert_eq!(extract_mi_base("123/A"), "123");
/// assert_eq!(extract_mi_base("123/B"), "123");
/// assert_eq!(extract_mi_base("abc/456/A"), "abc/456");
/// assert_eq!(extract_mi_base("123/C"), "123/C");  // not stripped
/// ```
#[must_use]
pub fn extract_mi_base(mi: &str) -> &str {
    if let Some(stripped) = mi.strip_suffix("/A") {
        stripped
    } else if let Some(stripped) = mi.strip_suffix("/B") {
        stripped
    } else {
        mi
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_consensus_reverse_returns_correct_tags() {
        let tags = TagSets::consensus_reverse();
        assert_eq!(tags.len(), 5);
        assert_eq!(tags, vec!["ad", "ae", "bd", "be", "cd"]);
    }

    #[test]
    fn test_consensus_revcomp_returns_correct_tags() {
        let tags = TagSets::consensus_revcomp();
        assert_eq!(tags.len(), 3);
        assert_eq!(tags, vec!["aD", "bD", "cD"]);
    }

    #[test]
    fn test_taginfo_new_with_empty_lists() {
        let tag_info = TagInfo::new(vec![], vec![], vec![]);
        assert!(tag_info.remove.is_empty());
        assert!(tag_info.reverse.is_empty());
        assert!(tag_info.revcomp.is_empty());
    }

    #[test]
    fn test_taginfo_new_with_individual_tags() {
        let tag_info = TagInfo::new(
            vec!["AS".to_string(), "NM".to_string()],
            vec!["BQ".to_string()],
            vec!["E2".to_string(), "U2".to_string()],
        );
        assert_eq!(tag_info.remove.len(), 2);
        assert!(tag_info.remove.contains("AS"));
        assert!(tag_info.remove.contains("NM"));
        assert_eq!(tag_info.reverse.len(), 1);
        assert!(tag_info.reverse.contains("BQ"));
        assert_eq!(tag_info.revcomp.len(), 2);
        assert!(tag_info.revcomp.contains("E2"));
        assert!(tag_info.revcomp.contains("U2"));
    }

    #[test]
    fn test_taginfo_new_with_consensus_reverse() {
        let tag_info = TagInfo::new(vec![], vec!["Consensus".to_string()], vec![]);
        assert_eq!(tag_info.reverse.len(), 5);
        assert!(tag_info.reverse.contains("ad"));
        assert!(tag_info.reverse.contains("ae"));
        assert!(tag_info.reverse.contains("bd"));
        assert!(tag_info.reverse.contains("be"));
        assert!(tag_info.reverse.contains("cd"));
    }

    #[test]
    fn test_taginfo_new_with_consensus_revcomp() {
        let tag_info = TagInfo::new(vec![], vec![], vec!["Consensus".to_string()]);
        assert_eq!(tag_info.revcomp.len(), 3);
        assert!(tag_info.revcomp.contains("aD"));
        assert!(tag_info.revcomp.contains("bD"));
        assert!(tag_info.revcomp.contains("cD"));
    }

    #[test]
    fn test_taginfo_new_with_consensus_and_individual_tags() {
        let tag_info = TagInfo::new(
            vec!["AS".to_string()],
            vec!["Consensus".to_string(), "BQ".to_string()],
            vec!["Consensus".to_string(), "E2".to_string()],
        );
        // Remove set
        assert_eq!(tag_info.remove.len(), 1);
        assert!(tag_info.remove.contains("AS"));
        // Reverse set should have Consensus tags + BQ
        assert_eq!(tag_info.reverse.len(), 6);
        assert!(tag_info.reverse.contains("ad"));
        assert!(tag_info.reverse.contains("ae"));
        assert!(tag_info.reverse.contains("bd"));
        assert!(tag_info.reverse.contains("be"));
        assert!(tag_info.reverse.contains("cd"));
        assert!(tag_info.reverse.contains("BQ"));
        // Revcomp set should have Consensus tags + E2
        assert_eq!(tag_info.revcomp.len(), 4);
        assert!(tag_info.revcomp.contains("aD"));
        assert!(tag_info.revcomp.contains("bD"));
        assert!(tag_info.revcomp.contains("cD"));
        assert!(tag_info.revcomp.contains("E2"));
    }

    #[test]
    fn test_taginfo_new_with_duplicate_tags() {
        let tag_info = TagInfo::new(
            vec!["AS".to_string(), "AS".to_string()],
            vec!["BQ".to_string(), "BQ".to_string()],
            vec!["E2".to_string(), "E2".to_string()],
        );
        // HashSet should deduplicate
        assert_eq!(tag_info.remove.len(), 1);
        assert_eq!(tag_info.reverse.len(), 1);
        assert_eq!(tag_info.revcomp.len(), 1);
    }

    #[test]
    fn test_taginfo_new_with_multiple_consensus_references() {
        let tag_info = TagInfo::new(
            vec![],
            vec!["Consensus".to_string(), "Consensus".to_string()],
            vec!["Consensus".to_string(), "Consensus".to_string()],
        );
        // Should not duplicate consensus tags
        assert_eq!(tag_info.reverse.len(), 5);
        assert_eq!(tag_info.revcomp.len(), 3);
    }

    #[test]
    fn test_has_revs_or_revcomps_with_both_empty() {
        let tag_info = TagInfo::new(vec!["AS".to_string()], vec![], vec![]);
        assert!(!tag_info.has_revs_or_revcomps());
    }

    #[test]
    fn test_has_revs_or_revcomps_with_reverse_only() {
        let tag_info = TagInfo::new(vec![], vec!["BQ".to_string()], vec![]);
        assert!(tag_info.has_revs_or_revcomps());
    }

    #[test]
    fn test_has_revs_or_revcomps_with_revcomp_only() {
        let tag_info = TagInfo::new(vec![], vec![], vec!["E2".to_string()]);
        assert!(tag_info.has_revs_or_revcomps());
    }

    #[test]
    fn test_has_revs_or_revcomps_with_both() {
        let tag_info = TagInfo::new(vec![], vec!["BQ".to_string()], vec!["E2".to_string()]);
        assert!(tag_info.has_revs_or_revcomps());
    }

    #[test]
    fn test_has_revs_or_revcomps_with_consensus() {
        let tag_info = TagInfo::new(vec![], vec!["Consensus".to_string()], vec![]);
        assert!(tag_info.has_revs_or_revcomps());
    }

    #[test]
    fn test_taginfo_clone() {
        let tag_info =
            TagInfo::new(vec!["AS".to_string()], vec!["BQ".to_string()], vec!["E2".to_string()]);
        let cloned = tag_info.clone();
        assert_eq!(cloned.remove.len(), tag_info.remove.len());
        assert_eq!(cloned.reverse.len(), tag_info.reverse.len());
        assert_eq!(cloned.revcomp.len(), tag_info.revcomp.len());
    }

    #[test]
    fn test_extract_mi_base_simple() {
        assert_eq!(extract_mi_base("123"), "123");
        assert_eq!(extract_mi_base("A"), "A");
    }

    #[test]
    fn test_extract_mi_base_with_duplex_suffix() {
        // Only /A and /B suffixes are stripped
        assert_eq!(extract_mi_base("123/A"), "123");
        assert_eq!(extract_mi_base("123/B"), "123");
        assert_eq!(extract_mi_base("123/456/A"), "123/456");
        // Other suffixes are preserved
        assert_eq!(extract_mi_base("123/C"), "123/C");
        assert_eq!(extract_mi_base("123/ReallyLongSuffix"), "123/ReallyLongSuffix");
    }

    // =========================================================================
    // validate_umi tests - matching fgbio behavior
    // =========================================================================

    #[test]
    fn test_validate_umi_uppercase_acgt() {
        assert_eq!(validate_umi(b"ACGT"), UmiValidation::Valid(4));
        assert_eq!(validate_umi(b"AAAAAAAA"), UmiValidation::Valid(8));
        assert_eq!(validate_umi(b"TTTTTTTT"), UmiValidation::Valid(8));
    }

    #[test]
    fn test_validate_umi_lowercase_acgt() {
        // Lowercase bases should be counted (matches fgbio's toUpperCase behavior)
        assert_eq!(validate_umi(b"acgt"), UmiValidation::Valid(4));
        assert_eq!(validate_umi(b"AcGt"), UmiValidation::Valid(4));
    }

    #[test]
    fn test_validate_umi_uppercase_n_rejected() {
        // Uppercase N should be rejected (matches fgbio's .contains('N'))
        assert_eq!(validate_umi(b"ACNT"), UmiValidation::ContainsN);
        assert_eq!(validate_umi(b"NACGT"), UmiValidation::ContainsN);
        assert_eq!(validate_umi(b"ACGTN"), UmiValidation::ContainsN);
        assert_eq!(validate_umi(b"NNNN"), UmiValidation::ContainsN);
    }

    #[test]
    fn test_validate_umi_lowercase_n_skipped() {
        // Lowercase 'n' should be skipped, not rejected (matches fgbio)
        // fgbio's .contains('N') is case-sensitive
        assert_eq!(validate_umi(b"ACnT"), UmiValidation::Valid(3));
        assert_eq!(validate_umi(b"acnt"), UmiValidation::Valid(3));
        assert_eq!(validate_umi(b"nnnn"), UmiValidation::Valid(0));
    }

    #[test]
    fn test_validate_umi_dash_skipped() {
        // Dash is skipped for paired UMIs
        assert_eq!(validate_umi(b"ACGT-TGCA"), UmiValidation::Valid(8));
        assert_eq!(validate_umi(b"----"), UmiValidation::Valid(0));
    }

    #[test]
    fn test_validate_umi_other_chars_skipped() {
        // Other characters are silently skipped
        assert_eq!(validate_umi(b"ACGT+TGCA"), UmiValidation::Valid(8));
        assert_eq!(validate_umi(b"AC GT"), UmiValidation::Valid(4));
    }

    #[test]
    fn test_validate_umi_empty() {
        assert_eq!(validate_umi(b""), UmiValidation::Valid(0));
    }

    #[test]
    fn test_validate_umi_mixed_case_with_uppercase_n() {
        // Mixed case with uppercase N should still be rejected
        assert_eq!(validate_umi(b"acNt"), UmiValidation::ContainsN);
        assert_eq!(validate_umi(b"AcNt"), UmiValidation::ContainsN);
    }
}

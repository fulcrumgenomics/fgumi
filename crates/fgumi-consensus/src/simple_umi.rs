//! Simple UMI consensus calling for metrics collection.
//!
//! This module provides position-by-position consensus calling for UMI sequences,
//! matching fgbio's `SimpleConsensusCaller` behavior. It uses likelihood-based
//! consensus calling at each position.

use super::base_builder::ConsensusBaseBuilder;

/// DNA bases recognized for consensus calling
const DNA_BASES: [u8; 5] = [b'A', b'C', b'G', b'T', b'N'];

/// Default error rate pre-labeling (Q90 - effectively no error)
const DEFAULT_ERROR_RATE_PRE: u8 = 90;

/// Default error rate post-labeling (Q90 - effectively no error)
const DEFAULT_ERROR_RATE_POST: u8 = 90;

/// Default quality score for observations (Q20)
const DEFAULT_Q_ERROR: u8 = 20;

/// Simple consensus caller for UMI sequences.
///
/// Uses position-by-position likelihood-based consensus calling, matching
/// fgbio's `SimpleConsensusCaller` behavior. All observations are treated
/// with equal quality (Q20 by default).
///
/// Non-DNA characters (like '-') are preserved from the first sequence
/// if all sequences have the same character at that position.
pub struct SimpleConsensusCaller {
    /// Consensus base builder for likelihood-based calling
    builder: ConsensusBaseBuilder,
    /// Quality score assigned to each observation
    q_error: u8,
}

impl SimpleConsensusCaller {
    /// Creates a new consensus caller with specified error rates.
    ///
    /// # Arguments
    /// * `error_rate_pre_labeling` - Phred-scaled error rate prior to UMI labeling
    /// * `error_rate_post_labeling` - Phred-scaled error rate after UMI labeling
    /// * `q_error` - Quality score to assign to each observation
    #[must_use]
    pub fn new(error_rate_pre_labeling: u8, error_rate_post_labeling: u8, q_error: u8) -> Self {
        Self {
            builder: ConsensusBaseBuilder::new(error_rate_pre_labeling, error_rate_post_labeling),
            q_error,
        }
    }

    /// Calls consensus on a list of sequences of the same length.
    ///
    /// # Arguments
    /// * `sequences` - Vector of sequences to build consensus from
    ///
    /// # Returns
    /// The consensus sequence string.
    ///
    /// # Panics
    /// * If `sequences` is empty
    /// * If sequences have different lengths
    /// * If sequences have different non-DNA characters at the same position
    /// * If a position has mixed DNA and non-DNA characters
    #[must_use]
    pub fn call_consensus(&mut self, sequences: &[String]) -> String {
        assert!(!sequences.is_empty(), "Can't call consensus on an empty set of sequences!");

        if sequences.len() == 1 {
            return sequences[0].clone();
        }

        let first = &sequences[0];
        let seq_len = first.len();

        // Verify all sequences have the same length
        assert!(
            sequences.iter().all(|s| s.len() == seq_len),
            "Sequences must all have the same length. Found {sequences:?}"
        );

        let mut result = String::with_capacity(seq_len);
        let sequences_len = sequences.len();

        for i in 0..seq_len {
            self.builder.reset();
            let mut non_dna_count = 0;
            let first_char = first.as_bytes()[i];

            for seq in sequences {
                let c = seq.as_bytes()[i];
                if Self::is_dna_base(c) {
                    self.builder.add(c, self.q_error);
                } else {
                    non_dna_count += 1;
                    // Verify all non-DNA bases are the same character
                    assert!(
                        first_char == c,
                        "Sequences must have character '{}' at position {}, found '{}'",
                        first_char as char,
                        i,
                        c as char
                    );
                }
            }

            if non_dna_count == 0 {
                // All DNA bases - call consensus
                let (base, _qual) = self.builder.call();
                result.push(base as char);
            } else if non_dna_count == sequences_len {
                // All non-DNA - preserve from first sequence
                result.push(first_char as char);
            } else {
                panic!(
                    "Sequences contained a mix of DNA and non-DNA characters at offset {i}: {sequences:?}"
                );
            }
        }

        result
    }

    /// Checks if a byte is a DNA base (A, C, G, T, N - case insensitive)
    fn is_dna_base(b: u8) -> bool {
        let upper = b.to_ascii_uppercase();
        DNA_BASES.contains(&upper)
    }
}

impl Default for SimpleConsensusCaller {
    fn default() -> Self {
        Self::new(DEFAULT_ERROR_RATE_PRE, DEFAULT_ERROR_RATE_POST, DEFAULT_Q_ERROR)
    }
}

// Keep the old API for backward compatibility with duplex_metrics
/// Simple consensus caller for UMIs within a tag family.
///
/// This is a wrapper around `SimpleConsensusCaller` that provides the old API
/// for backward compatibility.
pub struct SimpleUmiConsensusCaller {
    caller: SimpleConsensusCaller,
}

impl SimpleUmiConsensusCaller {
    /// Creates a new simple UMI consensus caller
    #[must_use]
    pub fn new() -> Self {
        Self { caller: SimpleConsensusCaller::default() }
    }

    /// Calls consensus on a list of UMIs, returning the consensus UMI and whether errors were corrected.
    ///
    /// Uses position-by-position likelihood-based consensus calling, matching
    /// fgbio's `SimpleConsensusCaller` behavior.
    ///
    /// # Arguments
    /// * `umis` - Vector of UMI strings observed in a family
    ///
    /// # Returns
    /// Tuple of (`consensus_umi`, `had_errors`) where `had_errors` indicates if any UMIs differed from consensus
    #[must_use]
    pub fn consensus(&mut self, umis: &[String]) -> (String, bool) {
        if umis.is_empty() {
            return (String::new(), false);
        }

        if umis.len() == 1 {
            return (umis[0].clone(), false);
        }

        // Check all UMIs have the same length - if not, return first (graceful degradation)
        let first_len = umis[0].len();
        if !umis.iter().all(|u| u.len() == first_len) {
            return (umis[0].clone(), true);
        }

        let consensus = self.caller.call_consensus(umis);
        let had_errors = umis.iter().any(|umi| *umi != consensus);

        (consensus, had_errors)
    }

    /// Calls consensus on paired UMIs (duplex UMIs with /A and /B suffixes)
    ///
    /// Returns (`consensus_base_umi`, `had_errors_in_a`, `had_errors_in_b`)
    #[must_use]
    pub fn consensus_duplex(&mut self, umis: &[String]) -> (String, bool, bool) {
        // Split UMIs into /A and /B groups
        let mut a_umis = Vec::new();
        let mut b_umis = Vec::new();

        for umi in umis {
            if umi.ends_with("/A") {
                a_umis.push(umi[..umi.len() - 2].to_string());
            } else if umi.ends_with("/B") {
                b_umis.push(umi[..umi.len() - 2].to_string());
            }
        }

        if a_umis.is_empty() && b_umis.is_empty() {
            return (String::new(), false, false);
        }

        // Get consensus for each strand
        let (a_consensus, a_errors) =
            if a_umis.is_empty() { (String::new(), false) } else { self.consensus(&a_umis) };

        let (b_consensus, b_errors) =
            if b_umis.is_empty() { (String::new(), false) } else { self.consensus(&b_umis) };

        // Return the base UMI (without /A or /B suffix)
        let base_umi = if a_consensus.is_empty() { b_consensus } else { a_consensus };

        (base_umi, a_errors, b_errors)
    }
}

impl Default for SimpleUmiConsensusCaller {
    fn default() -> Self {
        Self::new()
    }
}

/// Calls consensus on a set of UMI sequences using likelihood-based consensus.
///
/// This is the single shared function for UMI consensus calling across all
/// consensus callers (Vanilla, Duplex, CODEC), matching fgbio's behavior.
///
/// # Arguments
/// * `umis` - Slice of UMI strings observed in a family
///
/// # Returns
/// The consensus UMI string, or empty string if input is empty
#[must_use]
pub fn consensus_umis(umis: &[String]) -> String {
    if umis.is_empty() {
        return String::new();
    }
    if umis.len() == 1 {
        return umis[0].clone();
    }
    let mut caller = SimpleConsensusCaller::default();
    caller.call_consensus(umis)
}

#[cfg(test)]
mod tests {
    use super::*;

    // ============================================
    // Tests ported from fgbio SimpleConsensusCallerTest
    // ============================================

    #[test]
    #[should_panic(expected = "empty")]
    fn test_fail_if_no_sequences_given() {
        let mut caller = SimpleConsensusCaller::default();
        let _ = caller.call_consensus(&[]);
    }

    #[test]
    #[should_panic(expected = "same length")]
    fn test_fail_if_sequences_have_different_lengths() {
        let mut caller = SimpleConsensusCaller::default();
        let _ = caller.call_consensus(&["A".to_string(), "AC".to_string()]);
    }

    #[test]
    #[should_panic(expected = "mix of DNA and non-DNA")]
    fn test_fail_if_mixed_dna_and_non_dna() {
        let mut caller = SimpleConsensusCaller::default();
        // Third sequence has 'A' where others have '-'
        let _ = caller.call_consensus(&[
            "GATT-ACA".to_string(),
            "GATT-ACA".to_string(),
            "GATTAACA".to_string(),
        ]);
    }

    #[test]
    #[should_panic(expected = "must have character")]
    fn test_fail_if_non_dna_chars_differ() {
        let mut caller = SimpleConsensusCaller::default();
        // Different non-DNA characters at the same position
        let _ = caller.call_consensus(&["GATT-ACA".to_string(), "GATT+ACA".to_string()]);
    }

    #[test]
    fn test_consensus_from_sequences_that_agree() {
        let mut caller = SimpleConsensusCaller::default();
        assert_eq!(caller.call_consensus(&["A".to_string(), "A".to_string()]), "A");
        assert_eq!(
            caller.call_consensus(&["GATTACA".to_string(), "GATTACA".to_string()]),
            "GATTACA"
        );
    }

    #[test]
    fn test_consensus_from_sequences_that_differ() {
        let mut caller = SimpleConsensusCaller::default();

        // Four equal bases should result in N (tie)
        assert_eq!(
            caller.call_consensus(&[
                "A".to_string(),
                "C".to_string(),
                "G".to_string(),
                "T".to_string()
            ]),
            "N"
        );

        // Majority C should win
        assert_eq!(
            caller.call_consensus(&[
                "A".to_string(),
                "C".to_string(),
                "C".to_string(),
                "C".to_string()
            ]),
            "C"
        );

        // Majority C should win (different order)
        assert_eq!(
            caller.call_consensus(&[
                "C".to_string(),
                "C".to_string(),
                "C".to_string(),
                "A".to_string()
            ]),
            "C"
        );

        // With N bases - they don't contribute to consensus
        assert_eq!(
            caller.call_consensus(&[
                "GATTACA".to_string(),
                "GATTACA".to_string(),
                "GATTACA".to_string(),
                "NNNNNNN".to_string()
            ]),
            "GATTACA"
        );
    }

    #[test]
    fn test_gracefully_handle_non_acgtn_bases() {
        let mut caller = SimpleConsensusCaller::default();

        // Non-DNA character in middle
        assert_eq!(
            caller.call_consensus(&[
                "GATT-ACA".to_string(),
                "GATT-ACA".to_string(),
                "GATT-ACA".to_string()
            ]),
            "GATT-ACA"
        );

        // Non-DNA character at start
        assert_eq!(caller.call_consensus(&["XGAT".to_string(), "XGAT".to_string()]), "XGAT");

        // Non-DNA character at end
        assert_eq!(caller.call_consensus(&["GATY".to_string(), "GATY".to_string()]), "GATY");
    }

    // ============================================
    // Tests for SimpleUmiConsensusCaller wrapper
    // ============================================

    #[test]
    fn test_single_umi() {
        let mut caller = SimpleUmiConsensusCaller::default();
        let umis = vec!["ACGT".to_string()];
        let (consensus, had_errors) = caller.consensus(&umis);
        assert_eq!(consensus, "ACGT");
        assert!(!had_errors);
    }

    #[test]
    fn test_all_identical() {
        let mut caller = SimpleUmiConsensusCaller::default();
        let umis = vec!["ACGT".to_string(), "ACGT".to_string(), "ACGT".to_string()];
        let (consensus, had_errors) = caller.consensus(&umis);
        assert_eq!(consensus, "ACGT");
        assert!(!had_errors);
    }

    #[test]
    fn test_majority_vote() {
        let mut caller = SimpleUmiConsensusCaller::default();
        // Position-by-position: position 3 has G, G, G, G -> G wins
        let umis =
            vec!["ACGT".to_string(), "ACGT".to_string(), "ACGT".to_string(), "ACGG".to_string()];
        let (consensus, had_errors) = caller.consensus(&umis);
        assert_eq!(consensus, "ACGT");
        assert!(had_errors);
    }

    #[test]
    fn test_duplex_consensus() {
        let mut caller = SimpleUmiConsensusCaller::default();
        let umis = vec![
            "ACGT/A".to_string(),
            "ACGT/A".to_string(),
            "ACGG/A".to_string(),
            "TGCA/B".to_string(),
            "TGCA/B".to_string(),
        ];
        let (base_umi, a_errors, b_errors) = caller.consensus_duplex(&umis);
        assert_eq!(base_umi, "ACGT");
        assert!(a_errors); // ACGG differs from ACGT
        assert!(!b_errors); // All /B are TGCA
    }

    #[test]
    fn test_empty_umis() {
        let mut caller = SimpleUmiConsensusCaller::default();
        let (consensus, had_errors) = caller.consensus(&[]);
        assert_eq!(consensus, "");
        assert!(!had_errors);
    }

    #[test]
    fn test_different_length_umis_returns_first() {
        let mut caller = SimpleUmiConsensusCaller::default();
        let umis = vec!["ACGT".to_string(), "AC".to_string()];
        let (consensus, had_errors) = caller.consensus(&umis);
        assert_eq!(consensus, "ACGT"); // Returns first UMI
        assert!(had_errors); // Marks as had errors
    }

    #[test]
    fn test_position_by_position_consensus() {
        let mut caller = SimpleUmiConsensusCaller::default();
        // This tests the key difference from whole-string majority vote
        // Input: ["AACC", "CCAA"]
        // Position 0: A, C -> tie -> N
        // Position 1: A, C -> tie -> N
        // Position 2: C, A -> tie -> N
        // Position 3: C, A -> tie -> N
        let umis = vec!["AACC".to_string(), "CCAA".to_string()];
        let (consensus, had_errors) = caller.consensus(&umis);
        assert_eq!(consensus, "NNNN"); // All positions are ties
        assert!(had_errors);
    }

    #[test]
    fn test_position_consensus_with_majority() {
        let mut caller = SimpleUmiConsensusCaller::default();
        // Position 0: A, A, C -> A wins
        // Position 1: C, C, A -> C wins
        // Position 2: G, G, G -> G wins
        // Position 3: T, T, T -> T wins
        let umis = vec!["ACGT".to_string(), "ACGT".to_string(), "CAGT".to_string()];
        let (consensus, had_errors) = caller.consensus(&umis);
        assert_eq!(consensus, "ACGT");
        assert!(had_errors);
    }
}

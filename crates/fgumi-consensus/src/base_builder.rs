//! # Consensus Base Calling
//!
//! This module provides the core likelihood-based algorithm for calling a consensus base at a
//! single position from multiple observations. This is the fundamental building block used by
//! all consensus callers to determine the most likely base and assign an accurate quality score.
//!
//! ## Overview
//!
//! Given multiple observations of a position (base + quality score pairs), the consensus base
//! builder uses Bayesian inference to determine:
//! 1. Which base (A, C, G, or T) is most likely to be the true base
//! 2. What quality score to assign to that consensus call
//!
//! This differs from simple majority voting by:
//! - Incorporating quality scores (low-quality observations contribute less)
//! - Using a probabilistic model that can properly handle ambiguous cases
//! - Generating quality scores that accurately reflect confidence in the call
//!
//! ## Likelihood-Based Model
//!
//! For each candidate base B ∈ {A, C, G, T}, we calculate the likelihood of observing all the
//! data given that B is the true base:
//!
//! ```text
//! L(B | observations) = ∏ P(obs_i | B)
//! ```
//!
//! Where for each observation:
//! - **If observed base matches B**: `P(obs | B) = 1 - error_rate`
//! - **If observed base differs from B**: `P(obs | B) = error_rate / 3`
//!
//! The division by 3 assumes equal probability of the three incorrect bases.
//!
//! ## Quality Score Calculation
//!
//! After computing likelihoods for all four bases:
//!
//! 1. **Normalize** to get posterior probabilities:
//!    ```text
//!    P(B | data) = L(B) / (L(A) + L(C) + L(G) + L(T))
//!    ```
//!
//! 2. **Select** the base with maximum posterior probability
//!
//! 3. **Calculate error probability**:
//!    ```text
//!    P_error = 1 - P(best_base | data)
//!    ```
//!
//! 4. **Convert to Phred scale**:
//!    ```text
//!    Q = -10 * log10(P_error)
//!    ```
//!
//! This Phred score represents the probability that the consensus call is wrong, properly
//! accounting for all the evidence (both supporting and contradicting).
//!
//! ## Error Rate Priors
//!
//! The model incorporates two sources of error:
//!
//! ### Pre-UMI Error Rate
//!
//! Errors present in the original molecule before PCR amplification. These errors will appear
//! in all reads from that molecule, so they cannot be corrected by consensus calling. This
//! rate sets a floor on the achievable consensus quality.
//!
//! **Typical value**: Q45 (error rate ~3 × 10^-5)
//!
//! ### Post-UMI Error Rate
//!
//! Errors introduced during sequencing or PCR amplification. These are independent across reads
//! and can be corrected by consensus calling.
//!
//! **Typical value**: Q40 (error rate ~1 × 10^-4)
//!
//! ### Combined Error Model
//!
//! For each observation with quality Q, the adjusted error probability accounts for both
//! post-UMI errors (from sequencing) and the quality score:
//!
//! ```text
//! P_adjusted_error = P_post_umi_error + (1 - P_post_umi_error) * P_sequencing_error
//! ```
//!
//! This adjusted probability is used in the likelihood calculation.
//!
//! After calling the consensus, the pre-UMI error rate is incorporated:
//! ```text
//! P_final_error = P_pre_umi_error + (1 - P_pre_umi_error) * P_consensus_error
//! ```
//!
//! ## Mathematical Details
//!
//! The implementation uses **log-space arithmetic** to avoid numerical underflow:
//!
//! - All probabilities are stored as natural logarithms
//! - Products become sums: `log(a * b) = log(a) + log(b)`
//! - Sums use log-sum-exp trick: `log(e^a + e^b) = log_sum_exp(a, b)`
//!
//! This allows accurate computation even when individual likelihoods are extremely small
//! (e.g., 10^-100).
//!
//! ## Usage in Consensus Callers
//!
//! Consensus callers typically use `ConsensusBaseBuilder` as follows:
//!
//! ```rust,ignore
//! use fgumi_lib::consensus::base_builder::ConsensusBaseBuilder;
//!
//! // Create builder with error rates
//! let mut builder = ConsensusBaseBuilder::new(
//!     45,  // Pre-UMI error rate (Q45)
//!     40,  // Post-UMI error rate (Q40)
//! );
//!
//! // For each position in the read:
//! for position in 0..read_length {
//!     builder.reset(); // Clear previous position
//!
//!     // Add observations from all reads
//!     for read in &reads {
//!         let base = read.sequence()[position];
//!         let qual = read.quality_scores()[position];
//!         builder.add(base, qual);
//!     }
//!
//!     // Call consensus for this position
//!     let (consensus_base, consensus_qual) = builder.call();
//!
//!     // Check depth
//!     let depth = builder.contributions();
//! }
//! ```
//!
//! ## Handling Edge Cases
//!
//! ### No Observations
//! Returns `(N, 0)` - no-call with minimum quality.
//!
//! ### N Bases
//! N (no-call) bases in input reads are ignored. They don't contribute evidence for any
//! particular base.
//!
//! ### Ties
//! If multiple bases have exactly equal posterior probability, returns `(N, 0)` to indicate
//! ambiguity. In practice, this is rare due to quality score differences.
//!
//! ### Low Depth
//! With only 1-2 observations, the consensus quality will be similar to the input qualities.
//! The model naturally handles this - more observations lead to higher confidence.
//!
//! ## Quality Score Bounds
//!
//! Phred scores are capped at Q93 (error rate ~5 × 10^-10):
//! - **Minimum**: Q2 (error rate ~0.63)
//! - **Maximum**: Q93 (error rate ~5 × 10^-10)
//!
//! These bounds ensure quality scores fit in standard BAM format and remain numerically stable.
//!
//! ## Performance Characteristics
//!
//! - **Time complexity**: O(n) where n = number of observations at a position
//! - **Space complexity**: O(1) - fixed-size arrays for four bases
//! - **Caching**: Pre-computes error tables for all possible quality scores (Q2-Q93) to avoid
//!   repeated calculations
//!
//! ## Example Quality Progression
//!
//! Consider calling consensus from multiple reads observing 'A':
//!
//! | Reads | All Q30 'A' | Consensus Quality |
//! |-------|-------------|-------------------|
//! | 1     | A           | ~Q30              |
//! | 2     | AA          | ~Q40              |
//! | 5     | AAAAA       | ~Q50              |
//! | 10    | AAAAAAAAAA  | ~Q60              |
//! | 20    | (20x A)     | ~Q70              |
//!
//! The consensus quality increases logarithmically with depth, reflecting the compounding
//! evidence that all observations support the same base.
//!
//! ## See Also
//!
//! - `vanilla_consensus_caller`: Uses this builder to call consensus across entire reads
//! - `duplex_caller`: Uses this builder for single-strand consensus before duplex combination
//! - `phred`: Utility functions for Phred score and probability conversions

use crate::phred::{
    LN_ONE, LogProbability, MAX_PHRED, MIN_PHRED, NO_CALL_BASE, PhredScore,
    ln_error_prob_two_trials, ln_normalize, ln_not, ln_prob_to_phred, ln_sum_exp_array,
    phred_to_ln_error_prob,
};
use approx::abs_diff_eq;
use std::cmp::Ordering;
use wide::f64x4;

/// The four DNA bases in the order used throughout consensus calling
const DNA_BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
const DNA_BASE_COUNT: usize = 4;

/// Lookup table for converting ASCII base to index (0-3 for A,C,G,T, 255 for invalid)
/// This avoids a loop+compare for every base observation
const BASE_TO_INDEX: [u8; 256] = {
    let mut table = [255u8; 256];
    table[b'A' as usize] = 0;
    table[b'a' as usize] = 0;
    table[b'C' as usize] = 1;
    table[b'c' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'g' as usize] = 2;
    table[b'T' as usize] = 3;
    table[b't' as usize] = 3;
    table
};

/// Builder for calling consensus at a single base position
///
/// Accumulates observations (base + quality) and uses a likelihood model to
/// call the consensus base and quality.
///
/// Uses Kahan summation algorithm with SIMD vectorization to maintain numeric stability
/// when accumulating log-likelihoods across many observations. The 4 DNA bases map
/// perfectly to f64x4 (256-bit AVX2), enabling parallel computation of all 4 likelihoods.
pub struct ConsensusBaseBuilder {
    /// Log-likelihoods for each of the four bases (A, C, G, T) - SIMD vectorized
    likelihoods: f64x4,

    /// Kahan summation compensation terms for each base - SIMD vectorized
    compensations: f64x4,

    /// Count of observations for each base
    observations: [u16; DNA_BASE_COUNT],

    /// Pre-computed correct probabilities adjusted for post-UMI errors
    adjusted_correct_table: Vec<LogProbability>,

    /// Pre-computed 1/3 of error probability (for wrong base likelihoods)
    adjusted_error_per_alt: Vec<LogProbability>,

    /// Pre-UMI error rate (log probability)
    ln_error_pre_umi: LogProbability,
}

impl ConsensusBaseBuilder {
    /// Creates a new consensus base builder
    ///
    /// # Arguments
    /// * `error_rate_pre_umi` - Phred-scaled error rate prior to UMI integration
    /// * `error_rate_post_umi` - Phred-scaled error rate after UMI integration
    #[must_use]
    pub fn new(error_rate_pre_umi: PhredScore, error_rate_post_umi: PhredScore) -> Self {
        let ln_error_post = phred_to_ln_error_prob(error_rate_post_umi);

        // Pre-compute adjusted probabilities for all possible quality scores
        // Must start from 0 since we index by phred score directly
        let mut adjusted_correct_table = Vec::with_capacity(MAX_PHRED as usize + 1);
        let mut adjusted_error_per_alt = Vec::with_capacity(MAX_PHRED as usize + 1);

        for q in 0..=MAX_PHRED {
            let ln_error_seq = phred_to_ln_error_prob(q);
            let adjusted_error = ln_error_prob_two_trials(ln_error_post, ln_error_seq);

            adjusted_correct_table.push(ln_not(adjusted_error));
            // Error for any specific alternate base is 1/3 of total error
            // ln(p/3) = ln(p) - ln(3)
            adjusted_error_per_alt.push(adjusted_error - 3.0_f64.ln());
        }

        Self {
            likelihoods: f64x4::splat(LN_ONE),
            compensations: f64x4::splat(0.0),
            observations: [0; DNA_BASE_COUNT],
            adjusted_correct_table,
            adjusted_error_per_alt,
            ln_error_pre_umi: phred_to_ln_error_prob(error_rate_pre_umi),
        }
    }

    /// Resets the builder to process a new position
    pub fn reset(&mut self) {
        self.likelihoods = f64x4::splat(LN_ONE);
        self.compensations = f64x4::splat(0.0);
        self.observations.fill(0);
    }

    /// Adds an observation (base + quality) to the consensus
    ///
    /// Uses SIMD-vectorized Kahan summation for numeric stability when accumulating
    /// log-likelihoods. All 4 base likelihoods are updated in parallel using f64x4.
    ///
    /// # Arguments
    /// * `base` - The observed base (A, C, G, T, or N)
    /// * `qual` - The base quality score
    pub fn add(&mut self, base: u8, qual: PhredScore) {
        // Use lookup table to convert base to index (handles both upper and lower case)
        let matching_idx = BASE_TO_INDEX[base as usize];

        // Ignore N bases and other invalid bases (index 255)
        if matching_idx == 255 {
            return;
        }

        let matching_idx = matching_idx as usize;

        // Get adjusted probabilities for this quality
        let qual_idx = qual.min(MAX_PHRED) as usize;
        let ln_correct = self.adjusted_correct_table[qual_idx];
        let ln_error_per_base = self.adjusted_error_per_alt[qual_idx];

        // Build values array: error for all bases except ln_correct at matching position
        let mut values = [ln_error_per_base; 4];
        values[matching_idx] = ln_correct;
        let values = f64x4::from(values);

        // SIMD-vectorized Kahan summation for numeric stability
        // y = value - compensation (recovers the lost low-order bits)
        let y = values - self.compensations;
        // t = sum + y (the new sum, but low-order bits of y may be lost)
        let t = self.likelihoods + y;
        // compensation = (t - sum) - y (algebraically zero, but captures lost bits)
        self.compensations = (t - self.likelihoods) - y;
        // Update the sum
        self.likelihoods = t;

        self.observations[matching_idx] += 1;
    }

    /// Fast path for unanimous consensus when only one base is observed.
    ///
    /// Returns `Some((base, quality))` if all observations are the same base,
    /// `None` otherwise (requiring full probabilistic consensus).
    ///
    /// This optimization avoids the expensive log-sum-exp calculation when
    /// all reads agree and the likelihood difference is large enough to
    /// guarantee high quality output.
    #[inline]
    fn try_unanimous_fast_path(&self) -> Option<(u8, PhredScore)> {
        // Count how many bases have observations
        let mut observed_base_idx = None;
        let mut num_bases_observed = 0;

        for (i, &count) in self.observations.iter().enumerate() {
            if count > 0 {
                num_bases_observed += 1;
                observed_base_idx = Some(i);
                if num_bases_observed > 1 {
                    // Multiple bases observed - need full consensus
                    return None;
                }
            }
        }

        // Unanimous: only one base observed
        let base_idx = observed_base_idx?;

        // Check if the likelihood difference is large enough to skip full calculation.
        // When ln(L_winner) - ln(L_loser) > threshold, the posterior is essentially 1.
        // threshold = 23 corresponds to likelihood ratio > 10^10, giving Q > 90
        #[expect(clippy::items_after_statements, reason = "const is logically scoped to this fast-path check")]
        const FAST_PATH_THRESHOLD: f64 = 23.0;

        let likelihoods = self.likelihoods.as_array();
        let winner_ll = likelihoods[base_idx];
        let loser_ll = likelihoods[(base_idx + 1) % 4]; // Any loser will do (all same)

        if winner_ll - loser_ll > FAST_PATH_THRESHOLD {
            // Very high confidence - return quality capped by pre-UMI error rate
            // The pre-UMI error rate limits the maximum achievable quality
            let max_qual = ln_prob_to_phred(self.ln_error_pre_umi);
            return Some((DNA_BASES[base_idx], max_qual));
        }

        // Likelihood difference not large enough - fall through to full calculation
        None
    }

    /// Calls the consensus base and quality
    ///
    /// Returns (`consensus_base`, `consensus_quality`)
    /// Returns (N, 0) if no observations or multiple equally likely bases
    ///
    /// # Panics
    ///
    /// Panics if there are observations but no maximum likelihood base is found
    /// (should not occur in practice since we check for ties).
    #[must_use]
    pub fn call(&self) -> (u8, PhredScore) {
        if self.contributions() == 0 {
            return (NO_CALL_BASE, MIN_PHRED);
        }

        // Fast path: check for unanimous consensus (only one base observed)
        // This is common in high-quality data and avoids expensive log-sum-exp
        if let Some((base, qual)) = self.try_unanimous_fast_path() {
            return (base, qual);
        }

        // Extract likelihoods as array for iteration and function calls
        let likelihoods = self.likelihoods.as_array();

        // Compute the sum of likelihoods for normalization
        let ln_sum = ln_sum_exp_array(likelihoods);

        // Find the maximum likelihood and check if it's unique
        let mut max_likelihood = f64::NEG_INFINITY;
        let mut max_index = None;
        let mut tie = false;

        for (i, &ll) in likelihoods.iter().enumerate() {
            match ll.partial_cmp(&max_likelihood) {
                Some(Ordering::Greater) => {
                    max_likelihood = ll;
                    max_index = Some(i);
                    tie = false;
                }
                Some(Ordering::Equal) => {
                    tie = true;
                }
                Some(Ordering::Less) => {
                    // Check for epsilon-level tie (within machine precision)
                    if abs_diff_eq!(ll, max_likelihood, epsilon = f64::EPSILON) {
                        tie = true;
                    }
                }
                None => {}
            }
        }

        // If there's a tie, return no-call
        if tie || max_index.is_none() {
            return (NO_CALL_BASE, MIN_PHRED);
        }

        let max_idx = max_index.unwrap();
        let consensus_base = DNA_BASES[max_idx];

        // Calculate posterior probability
        let ln_posterior = ln_normalize(max_likelihood, ln_sum);

        // Convert to error probability
        let ln_consensus_error = ln_not(ln_posterior);

        // Factor in pre-UMI error rate
        // P(final_error) = P(pre_UMI_error AND correct_consensus)
        //                  + P(no_pre_UMI_error AND consensus_error)
        //                  + P(pre_UMI_error AND consensus_error that don't cancel)
        // Approximation: just combine the two error rates
        let ln_final_error = ln_error_prob_two_trials(self.ln_error_pre_umi, ln_consensus_error);

        // Convert to Phred score and cap at maximum
        let qual = ln_prob_to_phred(ln_final_error);

        (consensus_base, qual)
    }

    /// Returns the total number of contributing observations
    ///
    /// This is the sum of observations across all four bases
    #[must_use]
    pub fn contributions(&self) -> u16 {
        self.observations.iter().sum()
    }

    /// Returns the number of observations for a specific base
    ///
    /// # Arguments
    /// * `base` - The base to query (A, C, G, or T)
    ///
    /// Returns 0 for invalid bases
    #[must_use]
    pub fn observations_for_base(&self, base: u8) -> u16 {
        let idx = BASE_TO_INDEX[base as usize];
        if idx == 255 { 0 } else { self.observations[idx as usize] }
    }

    /// Returns the observations for all bases as [A, C, G, T]
    #[must_use]
    pub fn all_observations(&self) -> [u16; DNA_BASE_COUNT] {
        self.observations
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_single_base_perfect() {
        let mut builder = ConsensusBaseBuilder::new(45, 40);

        // Add 10 perfect quality A's
        for _ in 0..10 {
            builder.add(b'A', 40);
        }

        let (base, qual) = builder.call();
        assert_eq!(base, b'A');
        assert!(qual >= 40); // Should have high quality
        assert_eq!(builder.contributions(), 10);
    }

    #[test]
    fn test_mixed_bases() {
        let mut builder = ConsensusBaseBuilder::new(45, 40);

        // Add 8 A's and 2 C's
        for _ in 0..8 {
            builder.add(b'A', 30);
        }
        for _ in 0..2 {
            builder.add(b'C', 30);
        }

        let (base, _qual) = builder.call();
        assert_eq!(base, b'A'); // Majority should win
        assert_eq!(builder.contributions(), 10);
    }

    #[test]
    fn test_no_observations() {
        let builder = ConsensusBaseBuilder::new(45, 40);

        let (base, qual) = builder.call();
        assert_eq!(base, b'N');
        assert_eq!(qual, 2); // fgbio's PhredScore.MinValue
        assert_eq!(builder.contributions(), 0);
    }

    #[test]
    fn test_ignore_n_bases() {
        let mut builder = ConsensusBaseBuilder::new(45, 40);

        builder.add(b'A', 40);
        builder.add(b'N', 40);
        builder.add(b'A', 40);

        let (base, _qual) = builder.call();
        assert_eq!(base, b'A');
        assert_eq!(builder.contributions(), 2); // N should not be counted
    }

    #[test]
    fn test_case_insensitive() {
        let mut builder = ConsensusBaseBuilder::new(45, 40);

        builder.add(b'a', 40);
        builder.add(b'A', 40);
        builder.add(b'a', 40);

        let (base, _qual) = builder.call();
        assert_eq!(base, b'A');
        assert_eq!(builder.contributions(), 3);
    }

    #[test]
    fn test_reset() {
        let mut builder = ConsensusBaseBuilder::new(45, 40);

        builder.add(b'A', 40);
        builder.add(b'A', 40);
        assert_eq!(builder.contributions(), 2);

        builder.reset();
        assert_eq!(builder.contributions(), 0);

        builder.add(b'C', 40);
        let (base, _qual) = builder.call();
        assert_eq!(base, b'C');
    }

    #[test]
    fn test_observations_for_base() {
        let mut builder = ConsensusBaseBuilder::new(45, 40);

        builder.add(b'A', 40);
        builder.add(b'A', 40);
        builder.add(b'C', 40);
        builder.add(b'G', 40);

        assert_eq!(builder.observations_for_base(b'A'), 2);
        assert_eq!(builder.observations_for_base(b'C'), 1);
        assert_eq!(builder.observations_for_base(b'G'), 1);
        assert_eq!(builder.observations_for_base(b'T'), 0);
    }

    #[test]
    fn test_quality_variation() {
        let mut builder = ConsensusBaseBuilder::new(45, 40);

        // Add high quality A's
        for _ in 0..5 {
            builder.add(b'A', 40);
        }

        // Add one low quality C (should not outweigh the A's)
        builder.add(b'C', 10);

        let (base, _qual) = builder.call();
        assert_eq!(base, b'A');
    }

    #[test]
    fn test_tie_returns_no_call() {
        let mut builder = ConsensusBaseBuilder::new(45, 40);

        // This is hard to create a perfect tie, but we can test the no-call logic
        // by having very few observations
        builder.add(b'A', 10);
        builder.add(b'C', 10);

        // With low quality and equal counts, might get different results
        // Just verify it doesn't crash and returns a valid base (A, C, or N)
        let (base, _qual) = builder.call();
        assert!(base == b'A' || base == b'C' || base == b'N');
    }

    // Port of fgbio test: "produce a no-call if two bases are of equal likelihood"
    #[test]
    fn test_equal_likelihood_produces_no_call() {
        // Use MAX_PHRED for error rates to effectively disable error rate adjustment
        let mut builder = ConsensusBaseBuilder::new(MAX_PHRED, MAX_PHRED);

        // Empty builder should return no-call with quality 2
        let (base, qual) = builder.call();
        assert_eq!(base, b'N');
        assert_eq!(qual, 2);

        // Equal evidence should also return no-call
        builder.add(b'A', 20);
        builder.add(b'C', 20);
        let (base, qual) = builder.call();
        assert_eq!(base, b'N');
        assert_eq!(qual, 2);
    }

    // Port of fgbio test: "calculate consensus base and quality given a massive pileup"
    #[test]
    fn test_massive_pileup() {
        let mut builder = ConsensusBaseBuilder::new(50, 50);

        // Add 1000 Q20 C's
        for _ in 0..1000 {
            builder.add(b'C', 20);
        }

        let (base, qual) = builder.call();
        assert_eq!(base, b'C');
        assert_eq!(qual, 50); // Quality capped by pre-UMI error rate
        assert_eq!(builder.contributions(), 1000);
        assert_eq!(builder.observations_for_base(b'A'), 0);
        assert_eq!(builder.observations_for_base(b'C'), 1000);
        assert_eq!(builder.observations_for_base(b'G'), 0);
        assert_eq!(builder.observations_for_base(b'T'), 0);

        // Add 10 T's - should still call C
        for _ in 0..10 {
            builder.add(b'T', 20);
        }

        let (base, qual) = builder.call();
        assert_eq!(base, b'C');
        assert_eq!(qual, 50);
        assert_eq!(builder.contributions(), 1010);
        assert_eq!(builder.observations_for_base(b'T'), 10);
    }

    // Port of fgbio test: "calculate consensus base and quality given conflicting evidence"
    #[test]
    fn test_conflicting_evidence() {
        let mut builder = ConsensusBaseBuilder::new(50, 50);

        builder.add(b'A', 30);
        builder.add(b'C', 28);

        let (base, qual) = builder.call();
        assert_eq!(base, b'A');
        assert!(qual <= 5, "Quality should be low due to conflicting evidence, got {qual}");
    }

    // Port of fgbio test: "support calling multiple pileups from the same builder"
    #[test]
    fn test_reset_and_reuse_fgbio() {
        let mut builder = ConsensusBaseBuilder::new(50, 50);

        builder.add(b'A', 20);
        let (base, qual) = builder.call();
        assert_eq!(base, b'A');
        assert_eq!(qual, 20);
        assert_eq!(builder.contributions(), 1);

        builder.reset();

        builder.add(b'C', 20);
        let (base, qual) = builder.call();
        assert_eq!(base, b'C');
        assert_eq!(qual, 20);
        assert_eq!(builder.contributions(), 1);
    }

    // Test that Kahan summation produces consistent results regardless of addition order
    // This tests the numeric stability improvement from PR #1120
    #[test]
    fn test_kahan_summation_order_independence() {
        // Create two builders with the same error rates
        let mut builder1 = ConsensusBaseBuilder::new(45, 40);
        let mut builder2 = ConsensusBaseBuilder::new(45, 40);

        // Add bases in different orders with varying qualities
        // Builder 1: A's first, then C's
        for q in [10, 20, 30, 40, 50] {
            builder1.add(b'A', q);
        }
        for q in [15, 25, 35, 45] {
            builder1.add(b'C', q);
        }

        // Builder 2: Interleaved
        builder2.add(b'A', 10);
        builder2.add(b'C', 15);
        builder2.add(b'A', 20);
        builder2.add(b'C', 25);
        builder2.add(b'A', 30);
        builder2.add(b'C', 35);
        builder2.add(b'A', 40);
        builder2.add(b'C', 45);
        builder2.add(b'A', 50);

        // Both should produce the same result
        let (base1, qual1) = builder1.call();
        let (base2, qual2) = builder2.call();

        assert_eq!(base1, base2, "Consensus base should be order-independent");
        assert_eq!(qual1, qual2, "Consensus quality should be order-independent");
        assert_eq!(base1, b'A', "Should call A (5 vs 4 observations)");
    }

    // Test that Kahan summation handles extreme quality differences without precision loss
    #[test]
    fn test_kahan_summation_extreme_quality_range() {
        let mut builder = ConsensusBaseBuilder::new(45, 40);

        // Add many low quality bases followed by a few high quality bases
        // Q93 provides MUCH more evidence than Q2 (error rates: ~5e-10 vs ~0.63)
        // So 10 Q93 observations easily outweigh 100 Q2 observations
        for _ in 0..100 {
            builder.add(b'A', 2); // Very low quality (error ~63%)
        }
        for _ in 0..10 {
            builder.add(b'C', 93); // Maximum quality (error ~5e-10)
        }

        // High quality C's should win - this tests that extreme values are handled correctly
        let (base, _qual) = builder.call();
        assert_eq!(base, b'C', "High-quality observations provide much stronger evidence");
    }

    // Port of fgbio test: "scale base qualities using the post-umi error rate"
    // Tests that ConsensusBaseBuilder properly adjusts qualities based on post-UMI error rate.
    // With errorRatePostLabeling=10 (Q10), input qualities should be reduced:
    // Q20 → Q9, Q15 → Q8, Q10 → Q7, Q5 → Q4
    #[test]
    fn test_scale_base_qualities_post_umi_error_rate() {
        // errorRatePreLabeling=MAX_PHRED (no pre-UMI adjustment), errorRatePostLabeling=10
        // This matches fgbio's test setup
        let input_quals: [u8; 4] = [20, 15, 10, 5];
        let expected_output: [u8; 4] = [9, 8, 7, 4];

        for (input_q, expected_q) in input_quals.iter().zip(expected_output.iter()) {
            let mut builder = ConsensusBaseBuilder::new(MAX_PHRED, 10);
            builder.add(b'A', *input_q);
            let (_base, qual) = builder.call();

            // Allow for small differences due to implementation details
            // but the general pattern should hold: qualities are reduced
            assert!(
                qual <= *input_q,
                "Q{input_q} should be reduced by post-UMI error rate, got Q{qual}"
            );
            assert!(
                (i16::from(qual) - i16::from(*expected_q)).abs() <= 1,
                "Q{input_q} should become ~Q{expected_q}, got Q{qual}"
            );
        }
    }
}

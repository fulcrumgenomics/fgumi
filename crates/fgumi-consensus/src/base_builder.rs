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
//! ```text
//! use fgumi_consensus::base_builder::ConsensusBaseBuilder;
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
//! - `vanilla_caller`: Uses this builder to call consensus across entire reads
//! - `duplex_caller`: Uses this builder for single-strand consensus before duplex combination
//! - `phred`: Utility functions for Phred score and probability conversions

use crate::phred::{
    LN_ONE, LogProbability, MAX_PHRED, MIN_PHRED, NO_CALL_BASE, PhredScore,
    ln_error_prob_two_trials, ln_normalize, ln_not, ln_prob_to_phred, ln_sum_exp_array,
    phred_to_ln_error_prob,
};
use approx::abs_diff_eq;
use std::cmp::Ordering;
use std::sync::OnceLock;
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

    /// Count of observations for each base.
    ///
    /// `u32` (not `u16`) so very deep UMI families — a single position can exceed
    /// `65_535` reads in high-duplication libraries — accumulate without wrapping.
    /// The per-base depth is later clamped to `i16::MAX` only at tag emission, to
    /// match fgbio's `Short` consensus-depth representation.
    observations: [u32; DNA_BASE_COUNT],

    /// Pre-computed correct probabilities adjusted for post-UMI errors
    adjusted_correct_table: Vec<LogProbability>,

    /// Pre-computed 1/3 of error probability (for wrong base likelihoods)
    adjusted_error_per_alt: Vec<LogProbability>,

    /// Pre-UMI error rate (log probability)
    ln_error_pre_umi: LogProbability,

    /// Smallest winner-minus-loser log-likelihood gap at which the unanimous fast
    /// path is provably exact for this `ln_error_pre_umi`.
    ///
    /// Derived once per pre-UMI error rate by [`cached_min_exact_unanimous_gap`];
    /// see [`min_exact_unanimous_gap`] for why a single constant cannot serve every
    /// pre-UMI error rate.
    min_unanimous_gap: f64,
}

/// Smallest winner-minus-loser log-likelihood gap at which a unanimous call's
/// quality is guaranteed to round to the pre-UMI cap.
///
/// The fast path returns `ln_prob_to_phred(ln_error_pre_umi)` — the cap — without
/// running the full calculation. That is only correct once the consensus error has
/// fallen far enough below the pre-UMI error that combining them no longer changes
/// the rounded Phred. How far "far enough" is depends on the pre-UMI rate: the
/// larger it is, the smaller the consensus error must be relative to it, so the
/// required gap grows as the pre-UMI quality rises. A single hardcoded threshold is
/// therefore too small above some pre-UMI quality (inflating qualities) and
/// needlessly large below it (skipping fewer calls than it safely could).
///
/// For a unanimous position the four likelihoods are one winner and three equal
/// losers `gap` below it, so the posterior is `1 / (1 + 3e^-gap)` independent of the
/// winner's absolute value. That makes the required gap a pure function of
/// `ln_error_pre_umi`, resolved here by bisection over the real calculation rather
/// than a closed form — the quality chain floors, clamps and has a short-circuit
/// branch in `ln_error_prob_two_trials`, and reproducing all of that in algebra is
/// how the threshold drifts out of step with the code again. Quality is monotone
/// non-decreasing in the gap and bounded above by the cap, so the predicate
/// "rounds to the cap" flips exactly once and bisection is exact.
///
/// Returns the known-good (upper) end of the final bracket, so the answer is never
/// optimistic.
fn min_exact_unanimous_gap(ln_error_pre_umi: LogProbability) -> f64 {
    /// Gap wide enough that the consensus error underflows to nothing; the true
    /// answer is far below this for every representable pre-UMI rate.
    const MAX_GAP: f64 = 256.0;
    /// Enough halvings of `MAX_GAP` to converge past f64 resolution. The bracket
    /// collapses below an ULP of the returned gap by ~56 halvings across the whole
    /// representable pre-UMI range, so 64 is comfortably converged with headroom;
    /// further halvings are no-ops that leave the result bit-identical.
    const ITERATIONS: usize = 64;

    let cap = ln_prob_to_phred(ln_error_pre_umi);

    // Mirrors `call_full` for a unanimous position, using the same log-sum-exp
    // helper so the two agree bit for bit. The winner is pinned at 0 and the losers
    // at -gap; `ln_sum_exp_array` factors out the max, so the absolute offset the
    // real likelihoods carry does not change the result.
    let rounds_to_cap = |gap: f64| {
        let ln_sum = ln_sum_exp_array(&[0.0, -gap, -gap, -gap]);
        let ln_posterior = ln_normalize(0.0, ln_sum);
        let ln_consensus_error = ln_not(ln_posterior);
        let ln_final_error = ln_error_prob_two_trials(ln_error_pre_umi, ln_consensus_error);
        ln_prob_to_phred(ln_final_error) == cap
    };

    if rounds_to_cap(0.0) {
        return 0.0;
    }

    let (mut too_small, mut wide_enough) = (0.0_f64, MAX_GAP);
    for _ in 0..ITERATIONS {
        let midpoint = 0.5 * (too_small + wide_enough);
        if rounds_to_cap(midpoint) {
            wide_enough = midpoint;
        } else {
            too_small = midpoint;
        }
    }
    wide_enough
}

/// Memoized [`min_exact_unanimous_gap`], one slot per possible pre-UMI Phred.
///
/// The bisection above costs ~200 iterations of the full quality chain — several
/// thousand `exp`/`log1p` evaluations — which is far too expensive to repeat per
/// builder: `consensus_umis` constructs a fresh `SimpleConsensusCaller`, and hence a
/// fresh `ConsensusBaseBuilder`, for every UMI family. Since `PhredScore` is a `u8`,
/// 256 slots cover the input domain exactly, so every `u8` is memoized without
/// clamping. `--error-rate-pre-umi` accepts any `u8` (nothing validates the range),
/// and the cache intentionally supports all 256 of them.
static UNANIMOUS_GAP_CACHE: [OnceLock<f64>; 256] = [const { OnceLock::new() }; 256];

/// Returns the memoized minimum exact unanimous gap for a pre-UMI error rate.
///
/// The gap is a pure function of `error_rate_pre_umi`: it depends only on
/// `phred_to_ln_error_prob` of that one `u8`, so memoizing on the `u8` returns the
/// bit-identical `f64` the direct computation would have produced. The rate is taken
/// as a Phred score rather than as an already-converted `LogProbability` so that the
/// cache key and the value's provenance cannot drift apart.
fn cached_min_exact_unanimous_gap(error_rate_pre_umi: PhredScore) -> f64 {
    *UNANIMOUS_GAP_CACHE[error_rate_pre_umi as usize]
        .get_or_init(|| min_exact_unanimous_gap(phred_to_ln_error_prob(error_rate_pre_umi)))
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

        let ln_error_pre_umi = phred_to_ln_error_prob(error_rate_pre_umi);

        Self {
            likelihoods: f64x4::splat(LN_ONE),
            compensations: f64x4::splat(0.0),
            observations: [0; DNA_BASE_COUNT],
            adjusted_correct_table,
            adjusted_error_per_alt,
            ln_error_pre_umi,
            min_unanimous_gap: cached_min_exact_unanimous_gap(error_rate_pre_umi),
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
        // Beyond `min_unanimous_gap` the combined error rounds to the pre-UMI cap, so
        // the full calculation cannot return anything else. The bound is derived per
        // pre-UMI error rate in `min_exact_unanimous_gap`; it is not a constant,
        // because how wide the gap must be depends on that rate.
        let likelihoods = self.likelihoods.as_array();
        let winner_ll = likelihoods[base_idx];
        let loser_ll = likelihoods[(base_idx + 1) % 4]; // Any loser will do (all same)

        if winner_ll - loser_ll > self.min_unanimous_gap {
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

        self.call_full()
    }

    /// The full log-sum-exp consensus calculation, with no fast path.
    ///
    /// [`Self::call`] short-circuits to [`Self::try_unanimous_fast_path`] when that
    /// is provably exact; this is the answer it must agree with. Kept separate so
    /// the agreement is directly testable — see
    /// `test_unanimous_fast_path_matches_full_calculation`.
    fn call_full(&self) -> (u8, PhredScore) {
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

        let max_idx = max_index.expect("max_index is Some after is_none() check above");
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
    pub fn contributions(&self) -> u32 {
        self.observations.iter().sum()
    }

    /// Returns the number of observations for a specific base
    ///
    /// # Arguments
    /// * `base` - The base to query (A, C, G, or T)
    ///
    /// Returns 0 for invalid bases
    #[must_use]
    pub fn observations_for_base(&self, base: u8) -> u32 {
        let idx = BASE_TO_INDEX[base as usize];
        if idx == 255 { 0 } else { self.observations[idx as usize] }
    }

    /// Returns the observations for all bases as [A, C, G, T]
    #[must_use]
    pub fn all_observations(&self) -> [u32; DNA_BASE_COUNT] {
        self.observations
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    /// The per-base error count in the consensus callers is
    /// `contributions() - observations_for_base(base)`. Both operands must stay
    /// unsaturated `u32` counters: the callers clamp only the *stored* cD/cE
    /// values to fgbio's `Short` ceiling, and clamping the operands instead
    /// would collapse the difference toward zero on a deep pileup and understate
    /// cE. This pins that on a pileup deep enough that a `u16`/`i16` counter
    /// would have saturated.
    #[test]
    fn test_deep_pileup_keeps_error_count_exact_past_the_short_ceiling() {
        let short_ceiling = u32::from(i16::MAX.unsigned_abs()); // 32_767
        let matching = short_ceiling + 5_000;
        let mismatching = short_ceiling + 1_000;

        let mut builder = ConsensusBaseBuilder::new(45, 40);
        for _ in 0..matching {
            builder.add(b'A', 30);
        }
        for _ in 0..mismatching {
            builder.add(b'C', 30);
        }

        let depth = builder.contributions();
        assert_eq!(depth, matching + mismatching, "depth must not saturate");
        assert!(depth > short_ceiling, "precondition: the pileup exceeds the Short ceiling");

        let (base, _qual) = builder.call();
        assert_eq!(base, b'A', "the majority base should win");

        // The subtraction is exact because neither operand is clamped. Were both
        // saturated at the ceiling first, this difference would collapse to 0.
        let error_count = depth - builder.observations_for_base(base);
        assert_eq!(error_count, mismatching, "error count must be exact, not a clamped difference");
        assert_ne!(error_count, 0, "a saturated-operand subtraction would collapse to zero here");
    }

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
    fn test_contributions_do_not_overflow_at_u16_boundary() {
        // A UMI family deeper than u16::MAX (e.g. high-duplication cfDNA) must not
        // wrap the per-base observation counter. With a u16 counter, 70_000 single-base
        // observations wrapped to 70_000 % 65_536 = 4_464; the depth must stay 70_000.
        let mut builder = ConsensusBaseBuilder::new(45, 40);
        let n: u32 = 70_000;
        for _ in 0..n {
            builder.add(b'A', 40);
        }

        assert_eq!(builder.contributions(), n, "consensus depth must not wrap at the u16 boundary");
        assert_eq!(
            builder.observations_for_base(b'A'),
            n,
            "per-base observation count must not wrap at the u16 boundary"
        );
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

    /// SIMPLEX3-01: a Q0 observation of a *minority* base drives that base's
    /// likelihood lane to `−∞`. The consensus base must stay the majority (`A`) and
    /// the quality must NOT exceed the quality of the same pileup *without* the
    /// contradicting read — a conflicting observation cannot increase confidence.
    /// The pre-fix `ln_sum_exp_array` returned `−∞` for the normalizer whenever any
    /// lane was `−∞`, driving the posterior to `+∞` and inflating the quality to the
    /// pre-UMI cap (observed as fgumi Q45 vs fgbio Q44).
    #[test]
    fn test_neg_inf_lane_does_not_inflate_quality() {
        let clean = {
            let mut b = ConsensusBaseBuilder::new(45, 40);
            b.add(b'A', 30);
            b.add(b'A', 30);
            b.call()
        };
        let with_q0_conflict = {
            let mut b = ConsensusBaseBuilder::new(45, 40);
            b.add(b'A', 30);
            b.add(b'A', 30);
            b.add(b'C', 0); // Q0 -> lane[C] = -inf
            b.call()
        };
        assert_eq!(with_q0_conflict.0, b'A', "majority base must remain A");
        // Pin the exact fgbio-compatible qualities: the clean 2×A pileup caps at Q44,
        // and adding the Q0 (probability-0) `C` must leave it at Q44 — the pre-fix bug
        // produced Q45 by driving the normalizer to +∞. Asserting the exact value (not
        // just `conflict <= clean`) pins the fgbio baseline against future drift.
        assert_eq!(clean.1, 44, "clean 2xA pileup must be fgbio's Q44");
        assert_eq!(
            with_q0_conflict.1, 44,
            "a Q0 conflicting base must leave the quality at fgbio's Q44, not inflate it"
        );
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

    // ========================================================================
    // Unanimous fast-path exactness
    // ========================================================================

    /// The unanimous fast path must return exactly what the full calculation
    /// returns, at every pre-UMI error rate the CLI accepts.
    ///
    /// `-1/--error-rate-pre-umi` takes any `u8` and nothing validates it, but the
    /// fast path's gap threshold was a single hardcoded constant reasoned about for
    /// one quality regime. The gap actually required for the capped answer to be
    /// provably correct grows with the pre-UMI rate, so above some pre-UMI quality
    /// the shortcut fired on calls whose true answer is below the cap and inflated
    /// them.
    ///
    /// Sweeping the pre-UMI rate is the part that matters, and it is the axis the
    /// rest of the suite never varies: everything else runs at the Q45 default,
    /// where the old constant happened to be safe, so this bug was invisible to it.
    /// Pre-UMI and post-UMI are the `#[case]` axes so a regression names the rate it
    /// broke at; depth and observation quality are swept inside the body — together
    /// they only walk the gap across the threshold (it grows linearly with depth),
    /// and giving them their own cases multiplies into thousands of near-duplicate
    /// tests. The assertion still names the exact combination that failed.
    #[rstest]
    fn test_unanimous_fast_path_matches_full_calculation(
        #[values(20, 30, 40, 45, 50, 60, 69, 70, 80, 90, 93)] error_rate_pre_umi: PhredScore,
        #[values(10, 30, 40)] error_rate_post_umi: PhredScore,
    ) {
        for depth in [1usize, 2, 3, 4, 5, 8, 12, 20, 50] {
            for obs_qual in [10u8, 20, 30, 40, 60, 93] {
                let mut builder =
                    ConsensusBaseBuilder::new(error_rate_pre_umi, error_rate_post_umi);
                for _ in 0..depth {
                    builder.add(b'A', obs_qual);
                }

                let full = builder.call_full();
                let fast = builder.call();

                assert_eq!(
                    fast, full,
                    "unanimous fast path disagreed with the full calculation at \
                     pre-UMI Q{error_rate_pre_umi}, post-UMI Q{error_rate_post_umi}, \
                     depth {depth}, obs Q{obs_qual}: fast={fast:?} full={full:?}"
                );
            }
        }
    }

    /// The derived bound must both preserve the old threshold's safety at typical
    /// pre-UMI rates and be tighter than it there — the fix is not a slowdown at the
    /// default. Above Q69 the required gap exceeds the old hardcoded 23.0, which is
    /// exactly the range where the constant was silently inflating qualities.
    #[rstest]
    #[case::q40(40, false)]
    #[case::q45_default(45, false)]
    #[case::q69_last_safe(69, false)]
    #[case::q70_first_unsafe(70, true)]
    #[case::q93_max(93, true)]
    fn test_derived_gap_brackets_the_old_hardcoded_threshold(
        #[case] error_rate_pre_umi: PhredScore,
        #[case] expected_above_old_threshold: bool,
    ) {
        /// The constant this bound replaces.
        const OLD_THRESHOLD: f64 = 23.0;

        let gap = min_exact_unanimous_gap(phred_to_ln_error_prob(error_rate_pre_umi));
        assert_eq!(
            gap > OLD_THRESHOLD,
            expected_above_old_threshold,
            "pre-UMI Q{error_rate_pre_umi} needs a gap of {gap:.2}; expected it to be \
             {} 23.0",
            if expected_above_old_threshold { "above" } else { "at or below" }
        );
    }

    /// The required gap rises monotonically with pre-UMI quality — the property that
    /// makes a single constant unable to serve the whole range.
    #[test]
    fn test_required_gap_increases_with_pre_umi_quality() {
        let gaps: Vec<f64> =
            (0..=MAX_PHRED).map(|q| min_exact_unanimous_gap(phred_to_ln_error_prob(q))).collect();

        for (q, pair) in gaps.windows(2).enumerate() {
            assert!(
                pair[1] >= pair[0],
                "gap must not shrink as pre-UMI quality rises, but Q{} needs {:.4} and \
                 Q{} needs {:.4}",
                q,
                pair[0],
                q + 1,
                pair[1]
            );
        }
    }

    /// The memoized gap must equal the directly computed one for every representable
    /// pre-UMI Phred, bit for bit. Sweeps the whole `u8` domain — including the
    /// values above `MAX_PHRED`, which the CLI accepts unvalidated — because the cache
    /// is indexed by the raw `u8` and must be exact across every index it can ever be given.
    #[test]
    fn test_cached_gap_matches_uncached_across_the_whole_u8_domain() {
        for q in 0..=u8::MAX {
            let cached = cached_min_exact_unanimous_gap(q);
            let direct = min_exact_unanimous_gap(phred_to_ln_error_prob(q));
            assert_eq!(
                cached.to_bits(),
                direct.to_bits(),
                "memoized gap for Q{q} ({cached}) differs from the direct computation ({direct})"
            );
        }
    }

    /// The gap a builder stores must be derived from the **pre**-UMI rate, not the
    /// post-UMI one. The asymmetric cases are what discriminate a transposed
    /// argument; `umi_consensus_defaults` is symmetric and cannot, but is kept
    /// because Q90/Q90 is the pair `SimpleConsensusCaller` actually uses.
    #[rstest]
    #[case::default_cli_rates(45, 40)]
    #[case::pre_below_post(40, 45)]
    #[case::umi_consensus_defaults(90, 90)]
    #[case::max_pre_low_post(93, 2)]
    #[case::low_pre_max_post(2, 93)]
    fn test_builder_gap_is_keyed_on_the_pre_umi_rate(
        #[case] error_rate_pre_umi: PhredScore,
        #[case] error_rate_post_umi: PhredScore,
    ) {
        let builder = ConsensusBaseBuilder::new(error_rate_pre_umi, error_rate_post_umi);
        let expected = min_exact_unanimous_gap(phred_to_ln_error_prob(error_rate_pre_umi));

        assert_eq!(
            builder.min_unanimous_gap.to_bits(),
            expected.to_bits(),
            "builder built with pre={error_rate_pre_umi}, post={error_rate_post_umi} stored gap \
             {} but the pre-UMI rate requires {expected}",
            builder.min_unanimous_gap
        );
    }
}

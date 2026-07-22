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
//! Returns `(N, MIN_PHRED)` - a no-call at the minimum quality, Q2.
//!
//! ### N Bases
//! N (no-call) bases in input reads are ignored. They don't contribute evidence for any
//! particular base.
//!
//! ### Ties
//! If the greatest likelihood is not unique, returns `(N, MIN_PHRED)` to indicate ambiguity.
//! "Not unique" means within a few ULPs, not bit-for-bit equal: the likelihoods are running
//! sums, so a separation narrower than summation noise is not evidence. The tolerance is
//! relative to the magnitude, so it holds at every family size — see `TIE_TOLERANCE_ULPS`.
//! In practice ties are rare, since quality scores differ between observations.
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
use approx::ulps_eq;
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

/// The post-UMI-adjusted per-quality probability tables used by `ConsensusBaseBuilder`.
///
/// Both tables are indexed by the observed base quality (`0..=MAX_PHRED`) and depend *only*
/// on the post-UMI error rate they were computed for. Because that error rate and the base
/// quality are both integer Phred scores drawn from a tiny finite domain, the whole table is
/// a pure function of a single `u8` — so it can be computed once and shared, rather than
/// recomputed per `ConsensusBaseBuilder`.
struct AdjustedProbabilityTables {
    /// Log probability that an observation of quality `q` is correct, adjusted for post-UMI error.
    correct: Box<[LogProbability]>,

    /// Log probability of any one *specific* alternate base (total error / 3) for quality `q`.
    error_per_alt: Box<[LogProbability]>,
}

impl AdjustedProbabilityTables {
    /// Computes both tables for a given post-UMI error rate.
    ///
    /// This performs exactly the arithmetic the per-builder loop used to perform inline, in the
    /// same order, so the resulting `f64` values are bit-identical to the un-memoized version.
    fn compute(error_rate_post_umi: PhredScore) -> Self {
        let ln_error_post = phred_to_ln_error_prob(error_rate_post_umi);
        let ln_three = 3.0_f64.ln();

        let mut correct = Vec::with_capacity(MAX_PHRED as usize + 1);
        let mut error_per_alt = Vec::with_capacity(MAX_PHRED as usize + 1);

        for q in 0..=MAX_PHRED {
            let ln_error_seq = phred_to_ln_error_prob(q);
            let adjusted_error = ln_error_prob_two_trials(ln_error_post, ln_error_seq);

            correct.push(ln_not(adjusted_error));
            // Error for any specific alternate base is 1/3 of total error
            // ln(p/3) = ln(p) - ln(3)
            error_per_alt.push(adjusted_error - ln_three);
        }

        Self {
            correct: correct.into_boxed_slice(),
            error_per_alt: error_per_alt.into_boxed_slice(),
        }
    }
}

/// Lazily-populated cache of `AdjustedProbabilityTables`, one slot per possible post-UMI Phred.
///
/// The post-UMI error rate is a `u8`, so 256 slots cover the entire input domain exactly — no
/// clamping, and therefore no behaviour change for out-of-range Phred values. Each slot is
/// initialized at most once; after that, building a `ConsensusBaseBuilder` is a pair of atomic
/// loads instead of 94 `exp`/`log1p` evaluations plus two heap allocations. This matters because
/// `consensus_umis()` constructs a fresh caller (and thus a fresh builder) for every UMI family.
static ADJUSTED_TABLE_CACHE: [OnceLock<AdjustedProbabilityTables>; 256] =
    [const { OnceLock::new() }; 256];

/// Returns the shared, memoized adjusted-probability tables for a post-UMI error rate.
fn adjusted_tables(error_rate_post_umi: PhredScore) -> &'static AdjustedProbabilityTables {
    ADJUSTED_TABLE_CACHE[error_rate_post_umi as usize]
        .get_or_init(|| AdjustedProbabilityTables::compute(error_rate_post_umi))
}

/// Tolerance, in ULPs, within which two log-likelihoods are considered tied.
///
/// The likelihoods are running sums of log-probabilities, so their magnitude grows with
/// the family size: values near `-500` are routine for a large family. Any tie tolerance
/// must therefore be *relative* to that magnitude. An absolute epsilon such as
/// [`f64::EPSILON`] (the ULP of `1.0`, ~2.2e-16) is meaningless at that scale, because a
/// single ULP near `-500` is already ~5.7e-14 — roughly 250x larger — so no two distinct
/// values could ever compare equal. The failure inverts below unit scale, where the same
/// absolute epsilon is *far too wide*: a shallow family of high-quality observations sits
/// near `-1e-4`, and one `f64::EPSILON` spans hundreds of ULPs there, so clearly-separated
/// bases would read as tied. ULP comparison scales automatically in both directions, so this
/// one constant behaves consistently for a 2-read family and a 2000-read family alike —
/// provided the absolute epsilon is explicitly disabled at the comparison site (see
/// [`unique_max_index`]).
const TIE_TOLERANCE_ULPS: u32 = 4;

/// Returns the index of the strictly-greatest likelihood, or `None` when that maximum is
/// tied with another base and the call is therefore ambiguous (a no-call).
///
/// Two properties matter here, and both were violated by the previous inline implementation:
///
/// * **Order independence.** Ties are measured against the *final* maximum, not a running
///   one. Comparing against a running maximum meant a near-tie was only ever detected when
///   the tied lane happened to appear *after* the maximum; for likelihoods arriving in
///   ascending order every comparison took the `Greater` branch and no tie was recorded at
///   all. That made the outcome depend on the fixed `A,C,G,T` lane order rather than on the
///   data.
/// * **Scale independence.** See [`TIE_TOLERANCE_ULPS`].
///
/// `NaN` lanes are ignored. An array with no finite entry (no evidence for any base) yields
/// `None`, since there is nothing to call.
fn unique_max_index(likelihoods: &[LogProbability; DNA_BASE_COUNT]) -> Option<usize> {
    let mut max = f64::NEG_INFINITY;
    let mut max_index = None;
    for (index, &likelihood) in likelihoods.iter().enumerate() {
        if likelihood > max {
            max = likelihood;
            max_index = Some(index);
        }
    }

    // No finite evidence for any base (every lane `-inf` or `NaN`) is not a callable position.
    let max_index = max_index.filter(|_| max.is_finite())?;

    // The maximum always matches itself, so a unique maximum yields exactly one match.
    //
    // `epsilon = 0.0` is required, not cosmetic: `ulps_eq!` applies an *absolute* epsilon
    // before its ULP comparison, and that epsilon defaults to `f64::EPSILON` — the ULP of
    // `1.0`. Likelihoods below unit scale are routine (a single high-quality observation
    // contributes only ~`-1e-4`), and there one default epsilon spans hundreds of ULPs, so
    // leaving it at the default would silently widen the window far past
    // `TIE_TOLERANCE_ULPS` and re-introduce the scale dependence this function exists to
    // remove. Zeroing it makes the rule purely relative; exact equality still ties, since
    // `|a - b| <= 0.0` holds when `a == b`.
    let tied_lanes = likelihoods
        .iter()
        .filter(|&&ll| ulps_eq!(ll, max, epsilon = 0.0, max_ulps = TIE_TOLERANCE_ULPS))
        .count();

    (tied_lanes == 1).then_some(max_index)
}

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

    /// Pre-computed correct probabilities adjusted for post-UMI errors.
    ///
    /// Borrowed from the process-wide `ADJUSTED_TABLE_CACHE`; see `adjusted_tables`.
    adjusted_correct_table: &'static [LogProbability],

    /// Pre-computed 1/3 of error probability (for wrong base likelihoods).
    ///
    /// Borrowed from the process-wide `ADJUSTED_TABLE_CACHE`; see `adjusted_tables`.
    adjusted_error_per_alt: &'static [LogProbability],

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
        // The adjusted probability tables depend only on the post-UMI error rate, and both it
        // and the base quality are integer Phred scores, so the tables are memoized process-wide
        // and shared rather than recomputed (94 `exp`/`log1p` evaluations) per builder.
        let tables = adjusted_tables(error_rate_post_umi);

        Self {
            likelihoods: f64x4::splat(LN_ONE),
            compensations: f64x4::splat(0.0),
            observations: [0; DNA_BASE_COUNT],
            adjusted_correct_table: &tables.correct,
            adjusted_error_per_alt: &tables.error_per_alt,
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
        #[expect(
            clippy::items_after_statements,
            reason = "const is logically scoped to this fast-path check"
        )]
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
    /// Returns (N, [`MIN_PHRED`]) if there are no observations, or if the greatest likelihood
    /// is tied with another base — "tied" meaning within a few ULPs, not bit-for-bit equal
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

        // Find the maximum likelihood, which must be unique; a tie is a no-call.
        let Some(max_idx) = unique_max_index(likelihoods) else {
            return (NO_CALL_BASE, MIN_PHRED);
        };
        let max_likelihood = likelihoods[max_idx];
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

    /// Recomputes one adjusted-probability entry with the formula the builder used inline
    /// before the tables were memoized, in the same order, so a drift in `compute` shows up
    /// as a bit-level difference rather than an approximate one.
    fn reference_adjusted_entry(
        error_rate_post_umi: PhredScore,
        qual: PhredScore,
    ) -> (LogProbability, LogProbability) {
        let ln_error_post = phred_to_ln_error_prob(error_rate_post_umi);
        let ln_error_seq = phred_to_ln_error_prob(qual);
        let adjusted_error = ln_error_prob_two_trials(ln_error_post, ln_error_seq);
        (ln_not(adjusted_error), adjusted_error - 3.0_f64.ln())
    }

    /// The memoized tables must reproduce the pre-memoization inline math exactly, across the
    /// full `0..=MAX_PHRED` quality domain, for every post-UMI error rate the callers use.
    #[rstest]
    #[case::post_umi_q0(0)]
    #[case::post_umi_q10(10)]
    #[case::post_umi_default_q40(40)]
    #[case::post_umi_q45(45)]
    #[case::post_umi_max_phred(MAX_PHRED)]
    #[case::post_umi_above_max_phred(u8::MAX)]
    fn test_adjusted_tables_match_the_inline_formula_bit_for_bit(
        #[case] error_rate_post_umi: PhredScore,
    ) {
        let tables = adjusted_tables(error_rate_post_umi);

        assert_eq!(tables.correct.len(), MAX_PHRED as usize + 1, "table must cover every Phred");
        assert_eq!(tables.error_per_alt.len(), MAX_PHRED as usize + 1);

        for qual in 0..=MAX_PHRED {
            let (correct, error_per_alt) = reference_adjusted_entry(error_rate_post_umi, qual);
            let actual_correct = tables.correct[qual as usize];
            let actual_error_per_alt = tables.error_per_alt[qual as usize];

            // Compared as bit patterns, not approximately: the memoized computation is meant to
            // be the identical arithmetic in the identical order, so any drift is a regression.
            assert_eq!(
                actual_correct.to_bits(),
                correct.to_bits(),
                "correct probability drifted at post-UMI Q{error_rate_post_umi}, base Q{qual}: \
                 {actual_correct} != {correct}"
            );
            assert_eq!(
                actual_error_per_alt.to_bits(),
                error_per_alt.to_bits(),
                "alt-error probability drifted at post-UMI Q{error_rate_post_umi}, base Q{qual}: \
                 {actual_error_per_alt} != {error_per_alt}"
            );
        }
    }

    /// The cache must hand out the *same* tables on every call, otherwise the memoization is
    /// not actually happening and every builder pays the transcendental math again.
    #[rstest]
    #[case::post_umi_q0(0)]
    #[case::post_umi_default_q40(40)]
    #[case::post_umi_above_max_phred(u8::MAX)]
    fn test_adjusted_tables_are_memoized(#[case] error_rate_post_umi: PhredScore) {
        let first = adjusted_tables(error_rate_post_umi);
        let second = adjusted_tables(error_rate_post_umi);

        assert!(std::ptr::eq(first, second), "the cache slot must be initialized at most once");
        assert!(std::ptr::eq(first.correct.as_ptr(), second.correct.as_ptr()));
        assert!(std::ptr::eq(first.error_per_alt.as_ptr(), second.error_per_alt.as_ptr()));
    }

    /// Distinct post-UMI error rates must land in distinct cache slots — a shared slot would
    /// silently apply one family's error model to another's.
    #[test]
    fn test_adjusted_tables_are_keyed_by_post_umi_error_rate() {
        let q40 = adjusted_tables(40);
        let q45 = adjusted_tables(45);

        assert!(!std::ptr::eq(q40, q45), "different error rates must not share a slot");
        assert_ne!(q40.correct[30].to_bits(), q45.correct[30].to_bits());
    }

    /// The builder must borrow the shared tables rather than computing its own copy.
    #[test]
    fn test_builder_borrows_the_shared_tables() {
        let first = ConsensusBaseBuilder::new(45, 40);
        let second = ConsensusBaseBuilder::new(10, 40);

        assert!(std::ptr::eq(
            first.adjusted_correct_table.as_ptr(),
            second.adjusted_correct_table.as_ptr()
        ));
        assert!(std::ptr::eq(
            first.adjusted_error_per_alt.as_ptr(),
            second.adjusted_error_per_alt.as_ptr()
        ));
    }

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

    /// One ULP at the given magnitude, i.e. the gap to the next representable `f64`.
    pub(super) fn one_ulp(at: f64) -> f64 {
        let next = f64::from_bits(at.abs().to_bits() + 1);
        next - at.abs()
    }

    /// The value exactly `ulps` steps further from zero than `at` (more negative when `at` is).
    pub(super) fn step_away_from_zero(at: f64, ulps: u32) -> f64 {
        let magnitude = f64::from_bits(at.abs().to_bits() + u64::from(ulps));
        if at.is_sign_negative() { -magnitude } else { magnitude }
    }

    /// Exact ties must be no-calls regardless of lane order.
    ///
    /// These cases pass on the pre-fix implementation too (an exact tie always hit its
    /// `Equal` branch); they are kept as regression guards, not as evidence of the bug.
    #[rstest]
    #[case::tied_pair_leading([-500.0, -500.0, -600.0, -700.0], None)]
    #[case::tied_pair_trailing([-700.0, -600.0, -500.0, -500.0], None)]
    #[case::tied_pair_split([-500.0, -600.0, -500.0, -700.0], None)]
    #[case::unique_max_first([-100.0, -600.0, -500.0, -700.0], Some(0))]
    #[case::unique_max_last([-600.0, -500.0, -700.0, -100.0], Some(3))]
    #[case::all_four_tied([-500.0, -500.0, -500.0, -500.0], None)]
    fn exact_ties_are_no_calls(
        #[case] likelihoods: [LogProbability; DNA_BASE_COUNT],
        #[case] expected: Option<usize>,
    ) {
        assert_eq!(unique_max_index(&likelihoods), expected);
    }

    /// A near-tie must resolve the same way whichever lane order it arrives in.
    ///
    /// The pre-fix implementation compared each lane against a *running* maximum, so a
    /// near-tie was only detected when the tied lane appeared after the maximum. At unit
    /// scale (where the old absolute `f64::EPSILON` was comparable to one ULP) the same two
    /// values returned `None` in descending order but `Some(1)` in ascending order — the
    /// answer depended on the fixed `A,C,G,T` lane order rather than on the data.
    #[rstest]
    #[case::descending([-1.0, -1.0 - f64::EPSILON, -1.0e9, -1.0e9])]
    #[case::ascending([-1.0 - f64::EPSILON, -1.0, -1.0e9, -1.0e9])]
    fn near_ties_are_order_independent(#[case] likelihoods: [LogProbability; DNA_BASE_COUNT]) {
        assert_eq!(unique_max_index(&likelihoods), None);
    }

    /// The tie tolerance must scale with the magnitude of the likelihoods.
    ///
    /// A sum of log-probabilities for a large family sits far from `1.0`, where the old
    /// absolute `f64::EPSILON` tolerance was orders of magnitude smaller than a single ULP
    /// and so could never fire. A one-ULP separation is noise from summation order, not
    /// real evidence, and must read as a tie at every scale.
    #[rstest]
    #[case::single_observation(-1.0e-4)]
    #[case::sub_unit(-0.01)]
    #[case::near_zero(-1.0)]
    #[case::small(-10.0)]
    #[case::typical_family(-500.0)]
    #[case::large_family(-5000.0)]
    fn one_ulp_separation_is_a_tie_at_any_scale(#[case] magnitude: f64) {
        let ulp = one_ulp(magnitude);
        let likelihoods = [magnitude, magnitude - ulp, -1.0e9, -1.0e9];
        assert_eq!(
            unique_max_index(&likelihoods),
            None,
            "lanes separated by 1 ULP at {magnitude} must tie (ulp = {ulp:e})"
        );
    }

    /// A separation far wider than the tolerance must still resolve to a definite call,
    /// at the same magnitudes — otherwise the fix would simply no-call everything.
    #[rstest]
    #[case::single_observation(-1.0e-4)]
    #[case::sub_unit(-0.01)]
    #[case::near_zero(-1.0)]
    #[case::small(-10.0)]
    #[case::typical_family(-500.0)]
    #[case::large_family(-5000.0)]
    fn clear_separation_is_not_a_tie_at_any_scale(#[case] magnitude: f64) {
        let likelihoods = [magnitude, magnitude - 1.0, -1.0e9, -1.0e9];
        assert_eq!(unique_max_index(&likelihoods), Some(0));
    }

    /// Below unit scale the tie window must stay exactly [`TIE_TOLERANCE_ULPS`] wide.
    ///
    /// `ulps_eq!` checks an *absolute* epsilon before comparing ULPs, and that epsilon defaults
    /// to `f64::EPSILON` — the ULP of `1.0`. Likelihoods below `1.0` are routine (one
    /// high-quality observation contributes only ~`-1e-4`), and there the default epsilon spans
    /// hundreds of ULPs. Leaving it at the default would no-call lanes separated far beyond the
    /// tolerance, which is the same scale dependence — just inverted — that the fix removes at
    /// the large-magnitude end. This pins the explicit `epsilon = 0.0`.
    #[rstest]
    #[case::tenth(-0.1)]
    #[case::hundredth(-0.01)]
    #[case::single_observation(-1.0e-4)]
    fn tie_window_stays_ulp_relative_below_unit_scale(#[case] magnitude: f64) {
        let separated = step_away_from_zero(magnitude, TIE_TOLERANCE_ULPS + 1);
        let gap = (magnitude - separated).abs();

        assert!(
            gap < f64::EPSILON,
            "precondition: a {}-ULP gap at {magnitude} ({gap:e}) must fall inside the default \
             absolute epsilon ({:e}), or this case would pass without the fix",
            TIE_TOLERANCE_ULPS + 1,
            f64::EPSILON
        );

        let likelihoods = [magnitude, separated, -1.0e9, -1.0e9];
        assert_eq!(
            unique_max_index(&likelihoods),
            Some(0),
            "a {}-ULP gap at {magnitude} must resolve to a call, not a tie",
            TIE_TOLERANCE_ULPS + 1
        );
    }

    #[rstest]
    #[case::all_neg_inf([f64::NEG_INFINITY; DNA_BASE_COUNT], None)]
    #[case::single_finite_lane([f64::NEG_INFINITY, -3.0, f64::NEG_INFINITY, f64::NEG_INFINITY], Some(1))]
    #[case::nan_lanes_ignored([f64::NAN, -3.0, f64::NAN, -9.0], Some(1))]
    #[case::all_nan([f64::NAN; DNA_BASE_COUNT], None)]
    fn tie_detection_handles_non_finite_lanes(
        #[case] likelihoods: [LogProbability; DNA_BASE_COUNT],
        #[case] expected: Option<usize>,
    ) {
        assert_eq!(unique_max_index(&likelihoods), expected);
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
}

/// Parity — and *documented, intentional divergence* — against fgbio's tie rule.
///
/// fgbio decides "is the maximum unique?" in `MathUtil.maxWithIndex(requireUniqueMaximum=true)`
/// (`src/main/scala/com/fulcrumgenomics/util/MathUtil.scala`), which
/// `ConsensusCaller.call()` uses to emit a no-call. Rather than pin a captured fixture, these
/// tests port that method exactly and run both rules over the same inputs, so the oracle is
/// generated programmatically and the boundary between "we match fgbio" and "we deliberately
/// do not" is checked rather than asserted in prose.
///
/// fgbio's rule has the two defects this module's [`unique_max_index`] was written to remove:
///
/// 1. **Order dependence.** The tie test lives in the `else` branch of `v > max`, so it only
///    fires for a lane that arrives *after* the maximum. The same two values therefore
///    no-call in one lane order and call in the other.
/// 2. **Scale dependence.** The tolerance is the absolute constant `2^-52` (numerically
///    `f64::EPSILON`). At the magnitudes a real family reaches it is orders of magnitude
///    below one ULP, so it can never fire; below unit scale it spans hundreds of ULPs and
///    fires far too readily.
///
/// The divergence is confined to near-ties: on exact ties and on clearly-separated lanes —
/// which is what real pileups almost always produce, and what the ported fgbio behaviour
/// tests in the module above exercise end to end — the two rules agree.
#[cfg(test)]
mod fgbio_tie_oracle_tests {
    use super::tests::{one_ulp, step_away_from_zero};
    use super::*;
    use rstest::rstest;

    /// fgbio's `MathUtil.epsilon`, transcribed as `1.0 / Math.pow(2, 52)` (`2^52` spelled out
    /// so the constant stays exact and independent of `f64::EPSILON`, which
    /// `fgbio_epsilon_is_f64_epsilon` then proves it equals).
    const FGBIO_EPSILON: f64 = 1.0 / 4_503_599_627_370_496.0;

    /// A faithful port of fgbio's `MathUtil.maxWithIndex(ds, requireUniqueMaximum = true)`,
    /// returning `None` where fgbio returns a `maxIndex` of `-1`.
    ///
    /// Kept structurally identical to the Scala original — running maximum, tie test only in
    /// the `else` branch, `NaN` lanes skipped — so the divergences below are the real ones and
    /// not artifacts of a paraphrase. fgbio throws when every lane is `NaN`; that arm is
    /// mapped to `None` because a panic is not a comparable outcome.
    fn fgbio_unique_max_index(likelihoods: &[LogProbability; DNA_BASE_COUNT]) -> Option<usize> {
        let mut max = f64::MIN;
        let mut max_index: i32 = -1;
        let mut assigned = false;

        for (index, &value) in likelihoods.iter().enumerate() {
            if value.is_nan() {
                continue;
            }
            if !assigned || value > max {
                max = value;
                max_index = i32::try_from(index).expect("four lanes fit in an i32");
                assigned = true;
            } else if (value - max).abs() <= FGBIO_EPSILON {
                max_index = -1;
            }
        }

        if !assigned || max_index < 0 {
            return None;
        }
        Some(usize::try_from(max_index).expect("checked non-negative above"))
    }

    /// The constant fgbio compares against is precisely `f64::EPSILON`, so the "absolute
    /// epsilon" characterisation of fgbio's rule above is exact, not approximate.
    #[test]
    fn fgbio_epsilon_is_f64_epsilon() {
        assert_eq!(FGBIO_EPSILON.to_bits(), f64::EPSILON.to_bits());
    }

    /// Where a pileup is unambiguous — an exact tie, or lanes separated far beyond either
    /// tolerance — fgumi and fgbio must agree. This is the parity half of the contract.
    #[rstest]
    #[case::exact_tie_leading([-500.0, -500.0, -600.0, -700.0], None)]
    #[case::exact_tie_trailing([-700.0, -600.0, -500.0, -500.0], None)]
    #[case::all_four_tied([-500.0; DNA_BASE_COUNT], None)]
    #[case::clear_max_first([-100.0, -600.0, -500.0, -700.0], Some(0))]
    #[case::clear_max_last([-600.0, -500.0, -700.0, -100.0], Some(3))]
    #[case::clear_max_sub_unit([-0.01, -1.0, -2.0, -3.0], Some(0))]
    #[case::nan_lanes_ignored([f64::NAN, -3.0, f64::NAN, -9.0], Some(1))]
    fn unambiguous_pileups_match_fgbio(
        #[case] likelihoods: [LogProbability; DNA_BASE_COUNT],
        #[case] expected: Option<usize>,
    ) {
        assert_eq!(unique_max_index(&likelihoods), expected, "fgumi");
        assert_eq!(fgbio_unique_max_index(&likelihoods), expected, "fgbio oracle");
    }

    /// Divergence 1 — order dependence. fgbio answers differently for the same two values
    /// depending only on which lane they land in; fgumi no-calls either way.
    #[test]
    fn fgumi_is_order_independent_where_fgbio_is_not() {
        let descending = [-1.0, -1.0 - f64::EPSILON, -1.0e9, -1.0e9];
        let ascending = [-1.0 - f64::EPSILON, -1.0, -1.0e9, -1.0e9];

        assert_eq!(unique_max_index(&descending), None, "fgumi ties regardless of lane order");
        assert_eq!(unique_max_index(&ascending), None, "fgumi ties regardless of lane order");

        // The documented fgbio defect: the tie is seen only when the tied lane trails the max.
        assert_eq!(fgbio_unique_max_index(&descending), None);
        assert_eq!(
            fgbio_unique_max_index(&ascending),
            Some(1),
            "fgbio's running-maximum tie test cannot fire on an ascending lane order"
        );
    }

    /// Divergence 2a — at family scale fgbio's absolute epsilon is far below one ULP, so a
    /// one-ULP separation (summation-order noise, not evidence) reads as a definite call.
    #[rstest]
    #[case::typical_family(-500.0)]
    #[case::large_family(-5000.0)]
    fn fgumi_ties_summation_noise_that_fgbio_calls(#[case] magnitude: f64) {
        let ulp = one_ulp(magnitude);
        assert!(
            ulp > FGBIO_EPSILON,
            "precondition: one ULP at {magnitude} exceeds fgbio's epsilon"
        );

        let likelihoods = [magnitude, magnitude - ulp, -1.0e9, -1.0e9];

        assert_eq!(unique_max_index(&likelihoods), None, "fgumi treats 1 ULP as a tie");
        assert_eq!(
            fgbio_unique_max_index(&likelihoods),
            Some(0),
            "fgbio's absolute epsilon cannot reach one ULP at {magnitude}"
        );
    }

    /// Divergence 2b — the inverse, below unit scale: fgbio's epsilon spans hundreds of ULPs
    /// and no-calls lanes that are separated well past the tolerance. This is the case the
    /// explicit `epsilon = 0.0` in [`unique_max_index`] exists to avoid inheriting.
    #[rstest]
    #[case::tenth(-0.1)]
    #[case::hundredth(-0.01)]
    #[case::single_observation(-1.0e-4)]
    fn fgumi_calls_separations_that_fgbio_no_calls_below_unit_scale(#[case] magnitude: f64) {
        let separated = step_away_from_zero(magnitude, TIE_TOLERANCE_ULPS + 1);
        let gap = (magnitude - separated).abs();
        assert!(
            gap < FGBIO_EPSILON,
            "precondition: the gap at {magnitude} is inside fgbio's epsilon"
        );

        // Descending lane order, so fgbio's running-maximum tie test does fire — isolating the
        // scale defect from the order defect.
        let likelihoods = [magnitude, separated, -1.0e9, -1.0e9];

        assert_eq!(
            unique_max_index(&likelihoods),
            Some(0),
            "{} ULPs at {magnitude} is a real separation for fgumi",
            TIE_TOLERANCE_ULPS + 1
        );
        assert_eq!(
            fgbio_unique_max_index(&likelihoods),
            None,
            "fgbio's absolute epsilon swallows {} ULPs at {magnitude}",
            TIE_TOLERANCE_ULPS + 1
        );
    }

    /// Divergence 3 — a position with no evidence at all. fgbio assigns the first `-inf` lane
    /// as the maximum and calls it; fgumi requires a finite maximum and no-calls.
    #[test]
    fn fgumi_no_calls_an_all_neg_inf_position_that_fgbio_calls() {
        let likelihoods = [f64::NEG_INFINITY; DNA_BASE_COUNT];

        assert_eq!(unique_max_index(&likelihoods), None, "no finite evidence is not callable");
        assert_eq!(
            fgbio_unique_max_index(&likelihoods),
            Some(0),
            "fgbio assigns on the first non-NaN lane, even at -inf"
        );
    }
}

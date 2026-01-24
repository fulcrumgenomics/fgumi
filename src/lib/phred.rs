//! Phred score utilities and probability calculations.
//!
//! This module provides functions for working with Phred quality scores and converting
//! between Phred scores and probabilities. All probability calculations are done in
//! log-space (natural log) for numerical stability.
//!
//! The implementations in this module are designed to exactly match fgbio's
//! NumericTypes.scala for consistent consensus calling results.
//!
//! Key references:
//! - Equation (7) and (10) from <https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf>

use std::f64::consts::LN_10;

/// Natural log of 2, used as threshold in log1mexp (Equation 7)
const LN_TWO: f64 = std::f64::consts::LN_2; // ln(2) ≈ 0.693

/// Natural log of 4/3, used in probabilityOfErrorTwoTrials
const LN_FOUR_THIRDS: f64 = 0.287_682_072_451_780_9; // ln(4/3)

/// Natural log of 1.0 (which is 0)
pub const LN_ONE: f64 = 0.0;

/// Minimum Phred score (Q2, matching fgbio's PhredScore.MinValue)
pub const MIN_PHRED: u8 = 2;

/// Maximum Phred score we handle (Q93, matching `SAMUtils.MAX_PHRED_SCORE`)
pub const MAX_PHRED: u8 = 93;

/// No-call base character (matches fgbio's `NoCallBase`)
pub const NO_CALL_BASE: u8 = b'N';

/// Lowercase no-call base character
pub const NO_CALL_BASE_LOWER: u8 = b'n';

/// Precision constant used in Phred score conversion (matching fgbio)
const PHRED_PRECISION: f64 = 0.001;

/// Phred score type
pub type PhredScore = u8;

/// Log probability type (natural log)
pub type LogProbability = f64;

/// Converts a Phred score to a log probability of error
///
/// Phred score Q relates to error probability P by: Q = -10 * log10(P)
/// So: P = 10^(-Q/10)
/// And: ln(P) = ln(10^(-Q/10)) = -Q * ln(10) / 10
///
/// # Examples
/// ```
/// use fgumi_lib::phred::phred_to_ln_error_prob;
///
/// // Q10 corresponds to 10% error rate
/// let ln_error = phred_to_ln_error_prob(10);
/// assert!((ln_error - 0.1_f64.ln()).abs() < 1e-10);
///
/// // Q20 corresponds to 1% error rate
/// let ln_error = phred_to_ln_error_prob(20);
/// assert!((ln_error - 0.01_f64.ln()).abs() < 1e-10);
///
/// // Q30 corresponds to 0.1% error rate
/// let ln_error = phred_to_ln_error_prob(30);
/// assert!((ln_error - 0.001_f64.ln()).abs() < 1e-10);
/// ```
#[inline]
#[must_use]
pub fn phred_to_ln_error_prob(phred: PhredScore) -> LogProbability {
    -f64::from(phred) * LN_10 / 10.0
}

/// Converts a Phred score to a log probability of being correct (1 - error)
///
/// For numerical stability with high quality scores, we use:
/// ln(1 - P) ≈ ln(1 - e^ln(P))
///
/// # Examples
/// ```
/// use fgumi_lib::phred::phred_to_ln_correct_prob;
///
/// // Q30 = 0.1% error, so 99.9% correct
/// let ln_correct = phred_to_ln_correct_prob(30);
/// assert!((ln_correct - 0.999_f64.ln()).abs() < 1e-6);
///
/// // Q20 = 1% error, so 99% correct
/// let ln_correct = phred_to_ln_correct_prob(20);
/// assert!((ln_correct - 0.99_f64.ln()).abs() < 1e-6);
/// ```
#[inline]
#[must_use]
pub fn phred_to_ln_correct_prob(phred: PhredScore) -> LogProbability {
    let ln_error = phred_to_ln_error_prob(phred);
    ln_one_minus_exp(ln_error)
}

/// Converts a log probability to a Phred score
///
/// Q = -10 * log10(P) = -10 * ln(P) / ln(10)
///
/// Uses floor(value + precision) matching fgbio's PhredScore.fromLogProbability.
/// The result is clamped to the valid Phred score range [2, 93].
///
/// # Examples
/// ```
/// use fgumi_lib::phred::ln_prob_to_phred;
///
/// // 1% error (0.01) corresponds to Q20
/// let phred = ln_prob_to_phred(0.01_f64.ln());
/// assert_eq!(phred, 20);
///
/// // 0.1% error (0.001) corresponds to Q30
/// let phred = ln_prob_to_phred(0.001_f64.ln());
/// assert_eq!(phred, 30);
///
/// // Very low error probability is clamped to Q93
/// let phred = ln_prob_to_phred(1e-20_f64.ln());
/// assert_eq!(phred, 93);
/// ```
#[inline]
#[must_use]
pub fn ln_prob_to_phred(ln_prob: LogProbability) -> PhredScore {
    // Match fgbio's PhredScore.fromLogProbability:
    // if (lnProbError < MaxValueAsLogDouble) MaxValue
    // else Math.floor(-10.0 * (lnProbError/ Ln10) + Precision).toByte
    let max_value_as_log = phred_to_ln_error_prob(MAX_PHRED);
    if ln_prob < max_value_as_log {
        return MAX_PHRED;
    }
    let phred = (-10.0 * ln_prob / LN_10 + PHRED_PRECISION).floor();
    phred.clamp(f64::from(MIN_PHRED), f64::from(MAX_PHRED)) as PhredScore
}

/// Precise computation of log(1 + exp(x)).
///
/// Implements Equation (10) from <https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf>
/// Matches fgbio's `log1pexp` function exactly.
///
/// Thresholds from the paper and fgbio:
/// - x <= -37:   exp(x) is so small that log(1 + exp(x)) ≈ exp(x)
/// - x <= 18:    use log1p(exp(x)) for precision
/// - x <= 33.3:  use x + exp(-x) approximation
/// - x > 33.3:   exp(-x) is negligible, so log(1 + exp(x)) ≈ x
#[inline]
fn log1pexp(x: f64) -> f64 {
    if x <= -37.0 {
        x.exp()
    } else if x <= 18.0 {
        x.exp().ln_1p()
    } else if x <= 33.3 {
        x + (-x).exp()
    } else {
        x
    }
}

/// Computes ln(1 - e^x) for x < 0 in a numerically stable way.
///
/// This is equivalent to fgbio's `log1mexp(a)` where `a = -x` (positive).
/// Implements Equation (7) from <https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf>
///
/// For x >= -ln(2): use ln(-expm1(x)) to avoid catastrophic cancellation when e^x is close to 1
/// For x < -ln(2):  use ln1p(-exp(x)) which is stable when e^x is small
#[inline]
fn ln_one_minus_exp(x: f64) -> f64 {
    if x >= 0.0 {
        f64::NEG_INFINITY
    } else if x >= -LN_TWO {
        // When x is close to 0 (x >= -ln(2)), e^x is close to 1, so 1 - e^x is close to 0
        // Use expm1 to avoid catastrophic cancellation
        // ln(1 - e^x) = ln(-(e^x - 1)) = ln(-expm1(x))
        (-x.exp_m1()).ln()
    } else {
        // When x < -ln(2), e^x < 0.5, so 1 - e^x > 0.5, safe to use ln1p
        // ln(1 - e^x) = ln(1 + (-e^x)) = ln1p(-e^x)
        (-x.exp()).ln_1p()
    }
}

/// Computes log(exp(a) - exp(b)) where a > b.
///
/// Matches fgbio's `aOrNotB` function.
/// Uses: a + log1mexp(a - b) for numerical stability.
#[inline]
fn ln_a_minus_b(a: f64, b: f64) -> f64 {
    if b.is_infinite() && b < 0.0 {
        a
    } else if (a - b).abs() < f64::EPSILON {
        f64::NEG_INFINITY
    } else {
        // a + log(1 - exp(-(a-b))) = a + log1mexp(a-b)
        // In our convention, ln_one_minus_exp takes negative values, but here b < a
        a + ln_one_minus_exp(b - a)
    }
}

/// Computes the probability of observing an error across two independent trials.
///
/// Matches fgbio's `probabilityOfErrorTwoTrials` exactly.
///
/// Given two error probabilities p1 and p2 (in log space), this models the probability
/// of either or both reads being in error. For DNA sequencing with 4 bases:
///
/// f(X, Y) = X + Y - 4/3 * X * Y
///
/// where the 4/3 factor accounts for the probability of two errors resulting in
/// the same incorrect base (1/3 chance).
///
/// This is equivalent to:
/// P(error) = p1*(1-p2) + (1-p1)*p2 + p1*p2*(2/3)
///
/// The quick approximation is used when one probability is much smaller than the other
/// (difference of 6 or more in log space corresponds to ~400x difference).
///
/// # Examples
/// ```
/// use fgumi_lib::phred::ln_error_prob_two_trials;
///
/// // If both trials have 10% error, combined error using 4/3 formula:
/// // f(0.1, 0.1) = 0.1 + 0.1 - 4/3*0.1*0.1 = 0.2 - 0.0133... ≈ 0.1867
/// let ln_p1 = 0.1_f64.ln();
/// let ln_p2 = 0.1_f64.ln();
/// let result = ln_error_prob_two_trials(ln_p1, ln_p2);
/// let expected = 0.1 + 0.1 - (4.0/3.0) * 0.1 * 0.1;
/// assert!((result.exp() - expected).abs() < 1e-10);
/// ```
#[must_use]
pub fn ln_error_prob_two_trials(ln_p1: LogProbability, ln_p2: LogProbability) -> LogProbability {
    // Match fgbio's probabilityOfErrorTwoTrials exactly
    // Ensure ln_p1 >= ln_p2 (larger error probability first)
    let (ln_p1, ln_p2) = if ln_p1 < ln_p2 { (ln_p2, ln_p1) } else { (ln_p1, ln_p2) };

    // Quick approximation when one error probability is much smaller
    // A difference of 6 in log space is ~400x in linear space
    if ln_p1 - ln_p2 >= 6.0 {
        return ln_p1;
    }

    // Full formula: f(X, Y) = X + Y - 4/3 * X * Y
    // In log space:
    // term1 = log(X + Y) = ln_sum_exp(ln_p1, ln_p2)
    // term2 = log(4/3 * X * Y) = LN_FOUR_THIRDS + ln_p1 + ln_p2
    // result = log(term1 - term2) = ln_a_minus_b(term1, term2)
    let term1 = ln_sum_exp(ln_p1, ln_p2); // log(X + Y)
    let term2 = LN_FOUR_THIRDS + ln_p1 + ln_p2; // log(4/3 * X * Y)
    ln_a_minus_b(term1, term2)
}

/// Computes log(a + b) given log(a) and log(b)
///
/// Matches fgbio's `or` function using log1pexp for numerical stability:
/// log(a + b) = log(a) + log(1 + exp(log(b) - log(a)))
///            = log(a) + log1pexp(log(b) - log(a))
///
/// This is essential for adding probabilities while working in log-space, avoiding
/// numerical underflow that would occur when converting back to linear space.
///
/// # Examples
/// ```
/// use fgumi_lib::phred::ln_sum_exp;
///
/// // ln(0.1) + ln(0.2) should give ln(0.3)
/// let result = ln_sum_exp(0.1_f64.ln(), 0.2_f64.ln());
/// assert!((result - 0.3_f64.ln()).abs() < 1e-10);
///
/// // Works with very small probabilities
/// let result = ln_sum_exp(1e-100_f64.ln(), 2e-100_f64.ln());
/// assert!((result - 3e-100_f64.ln()).abs() < 1e-10);
/// ```
#[must_use]
pub fn ln_sum_exp(ln_a: LogProbability, ln_b: LogProbability) -> LogProbability {
    if ln_a.is_infinite() && ln_a < 0.0 {
        return ln_b;
    }
    if ln_b.is_infinite() && ln_b < 0.0 {
        return ln_a;
    }

    // Match fgbio's `or` function: ensure ln_a <= ln_b, then use log1pexp
    let (ln_a, ln_b) = if ln_b < ln_a { (ln_b, ln_a) } else { (ln_a, ln_b) };
    ln_a + log1pexp(ln_b - ln_a)
}

/// Computes log(sum(exp(values))) for an array of log probabilities
///
/// This is the normalization constant for converting log-likelihoods to
/// log-probabilities (posteriors). Uses the log-sum-exp trick for numerical
/// stability across all values.
///
/// # Examples
/// ```
/// use fgumi_lib::phred::ln_sum_exp_array;
///
/// // ln(0.1) + ln(0.2) + ln(0.3) should give ln(0.6)
/// let values = vec![0.1_f64.ln(), 0.2_f64.ln(), 0.3_f64.ln()];
/// let result = ln_sum_exp_array(&values);
/// assert!((result - 0.6_f64.ln()).abs() < 1e-10);
///
/// // Empty array returns negative infinity
/// let result = ln_sum_exp_array(&[]);
/// assert!(result.is_infinite() && result < 0.0);
/// ```
#[must_use]
pub fn ln_sum_exp_array(values: &[LogProbability]) -> LogProbability {
    if values.is_empty() {
        return f64::NEG_INFINITY;
    }

    let mut min_value = f64::INFINITY;
    let mut min_index = 0;
    for (i, value) in values.iter().enumerate() {
        if *value < min_value {
            min_index = i;
            min_value = *value;
        }
    }
    if min_value.is_infinite() {
        return min_value;
    }
    let mut sum = min_value;
    for (i, value) in values.iter().enumerate() {
        if i != min_index {
            sum = ln_sum_exp(sum, *value);
        }
    }
    sum
}

/// Normalizes a log probability by a normalization constant
///
/// Returns `ln(exp(ln_value)` / `exp(ln_norm)`) = `ln_value` - `ln_norm`
#[inline]
#[must_use]
pub fn ln_normalize(ln_value: LogProbability, ln_norm: LogProbability) -> LogProbability {
    ln_value - ln_norm
}

/// Computes ln(1 - exp(x)) = ln(NOT(exp(x)))
#[inline]
#[must_use]
pub fn ln_not(x: LogProbability) -> LogProbability {
    ln_one_minus_exp(x)
}

#[cfg(test)]
#[allow(clippy::similar_names)]
mod tests {
    use super::*;

    #[test]
    fn test_phred_to_ln_error() {
        // Q10 = 10% error = 0.1
        let q10_error = phred_to_ln_error_prob(10);
        assert!((q10_error - 0.1_f64.ln()).abs() < 1e-10);

        // Q20 = 1% error = 0.01
        let q20_error = phred_to_ln_error_prob(20);
        assert!((q20_error - 0.01_f64.ln()).abs() < 1e-10);

        // Q30 = 0.1% error = 0.001
        let q30_error = phred_to_ln_error_prob(30);
        assert!((q30_error - 0.001_f64.ln()).abs() < 1e-10);
    }

    #[test]
    fn test_phred_round_trip() {
        // Note: Q0 and Q1 will be clamped to MIN_PHRED (2) on the round trip
        for q in [MIN_PHRED, 10, 20, 30, 40, 50, 60] {
            let ln_p = phred_to_ln_error_prob(q);
            let q_back = ln_prob_to_phred(ln_p);
            assert_eq!(q, q_back);
        }
    }

    #[test]
    fn test_ln_sum_exp() {
        // ln(0.1) + ln(0.2) should give ln(0.3)
        let result = ln_sum_exp(0.1_f64.ln(), 0.2_f64.ln());
        assert!((result - 0.3_f64.ln()).abs() < 1e-10);
    }

    #[test]
    fn test_ln_sum_exp_array() {
        // ln(0.1) + ln(0.2) + ln(0.3) should give ln(0.6)
        let values = vec![0.1_f64.ln(), 0.2_f64.ln(), 0.3_f64.ln()];
        let result = ln_sum_exp_array(&values);
        assert!((result - 0.6_f64.ln()).abs() < 1e-10);
    }

    #[test]
    fn test_error_two_trials() {
        // If both trials have 10% error, combined error using 4/3 formula:
        // f(0.1, 0.1) = 0.1 + 0.1 - 4/3*0.1*0.1 = 0.2 - 0.01333... ≈ 0.18667
        let ln_p1 = 0.1_f64.ln();
        let ln_p2 = 0.1_f64.ln();
        let result = ln_error_prob_two_trials(ln_p1, ln_p2);
        let expected = 0.1 + 0.1 - (4.0 / 3.0) * 0.1 * 0.1;
        assert!((result.exp() - expected).abs() < 1e-10);
    }

    /// Ported from fgbio's NumericTypesTest.scala
    /// Tests the comprehensive probabilityOfErrorTwoTrials formula
    #[test]
    fn test_error_two_trials_comprehensive_fgbio() {
        // From fgbio: for (i <- 1 to 100; j <- 1 to 100)
        // expected = (p1*(1-p2)) + ((1-p1)*p2) + (p1 * p2 * 2/3)
        // This simplifies to: p1 + p2 - 4/3*p1*p2
        for i in 1..=100 {
            for j in 1..=100 {
                let p1 = 1.0 / f64::from(i);
                let p2 = 1.0 / f64::from(j);
                let expected = (p1 * (1.0 - p2)) + ((1.0 - p1) * p2) + (p1 * p2 * (2.0 / 3.0));
                let ln_p1 = p1.ln();
                let ln_p2 = p2.ln();
                let actual = ln_error_prob_two_trials(ln_p1, ln_p2);
                assert!(
                    (actual.exp() - expected).abs() < 0.0001,
                    "Failed for i={}, j={}: actual={}, expected={}",
                    i,
                    j,
                    actual.exp(),
                    expected
                );
            }
        }
    }

    /// Ported from fgbio's `NumericTypesTest`: "or in log space"
    #[test]
    fn test_ln_sum_exp_fgbio() {
        // exp(or(10, 10)) shouldBe exp(10)*2
        let result = ln_sum_exp(10.0, 10.0);
        assert!((result.exp() - 10.0_f64.exp() * 2.0).abs() < 0.00001);

        // exp(or(10, 20)) shouldBe exp(10)+exp(20)
        let result = ln_sum_exp(10.0, 20.0);
        assert!((result.exp() - (10.0_f64.exp() + 20.0_f64.exp())).abs() < 0.00001);

        // exp(or(20, 10)) shouldBe exp(20)+exp(10)
        let result = ln_sum_exp(20.0, 10.0);
        assert!((result.exp() - (20.0_f64.exp() + 10.0_f64.exp())).abs() < 0.00001);

        // or(10, LnZero) shouldBe exp(10)
        let result = ln_sum_exp(10.0, f64::NEG_INFINITY);
        assert!((result.exp() - 10.0_f64.exp()).abs() < 0.00001);

        // or(LnZero, 10) shouldBe exp(10)
        let result = ln_sum_exp(f64::NEG_INFINITY, 10.0);
        assert!((result.exp() - 10.0_f64.exp()).abs() < 0.00001);

        // or(-718.39..., -8.40...) shouldBe -8.40... (very small + small = small)
        let result = ln_sum_exp(-718.394_775_628_242_3, -8.404_216_861_178_751);
        assert!((result - (-8.404_216_861_178_751)).abs() < 0.00001);
    }

    /// Ported from fgbio's `NumericTypesTest`: "aOrNotB in log space"
    #[test]
    fn test_ln_a_minus_b_fgbio() {
        let q10 = phred_to_ln_error_prob(10); // ln(0.1)
        let q20 = phred_to_ln_error_prob(20); // ln(0.01)

        // aOrNotB(10, 10) shouldBe LnZero
        assert!(ln_a_minus_b(10.0, 10.0).is_infinite() && ln_a_minus_b(10.0, 10.0) < 0.0);

        // aOrNotB(q10, q10) shouldBe LnZero
        assert!(ln_a_minus_b(q10, q10).is_infinite() && ln_a_minus_b(q10, q10) < 0.0);

        // exp(aOrNotB(q10, q20)) shouldBe (0.1-0.01)
        let result = ln_a_minus_b(q10, q20);
        assert!((result.exp() - 0.09).abs() < 0.00001, "Expected 0.09, got {}", result.exp());

        // aOrNotB(LnTen, LnZero) shouldBe LnTen
        let ln_ten = 10.0_f64.ln();
        let result = ln_a_minus_b(ln_ten, f64::NEG_INFINITY);
        assert!((result - ln_ten).abs() < 0.00001);
    }

    /// Ported from fgbio's `NumericTypesTest`: "1 - probability in log space"
    #[test]
    fn test_ln_one_minus_exp_fgbio() {
        let q10 = phred_to_ln_error_prob(10); // ln(0.1)
        let q20 = phred_to_ln_error_prob(20); // ln(0.01)

        // exp(not(q10)) shouldBe 0.9
        let result = ln_one_minus_exp(q10);
        assert!((result.exp() - 0.9).abs() < 0.00001);

        // exp(not(q20)) shouldBe 0.99
        let result = ln_one_minus_exp(q20);
        assert!((result.exp() - 0.99).abs() < 0.00001);

        // exp(not(ln(0.90))) shouldBe 0.1
        let result = ln_one_minus_exp(0.90_f64.ln());
        assert!((result.exp() - 0.1).abs() < 0.00001);

        // exp(not(ln(0.99))) shouldBe 0.01
        let result = ln_one_minus_exp(0.99_f64.ln());
        assert!((result.exp() - 0.01).abs() < 0.00001);

        // exp(not(LnZero)) shouldBe 1.0
        let result = ln_one_minus_exp(f64::NEG_INFINITY);
        assert!((result.exp() - 1.0).abs() < 0.00001);
    }

    /// Ported from fgbio's `NumericTypesTest`: phred score conversions
    #[test]
    fn test_phred_conversions_fgbio() {
        // fromLogProbability(toLogProbability(0.0)) shouldBe MaxValue
        let ln_zero = 0.0_f64.ln(); // NEG_INFINITY
        assert_eq!(ln_prob_to_phred(ln_zero), MAX_PHRED);

        // fromLogProbability(toLogProbability(0.1)) shouldBe 10
        let ln_01 = 0.1_f64.ln();
        assert_eq!(ln_prob_to_phred(ln_01), 10);

        // fromLogProbability(toLogProbability(0.5)) shouldBe 3
        let ln_05 = 0.5_f64.ln();
        assert_eq!(ln_prob_to_phred(ln_05), 3);

        // fromLogProbability(toLogProbability(1.0)) shouldBe 0
        let ln_1 = 1.0_f64.ln(); // 0.0
        // Note: fgbio returns 0, but our MIN_PHRED is 2, so we clamp to 2
        assert_eq!(ln_prob_to_phred(ln_1), MIN_PHRED);
    }

    /// Test log1pexp implementation (Equation 10)
    #[test]
    fn test_log1pexp() {
        // Test threshold regions
        // x <= -37: returns exp(x)
        assert!((log1pexp(-50.0) - (-50.0_f64).exp()).abs() < 1e-10);
        assert!((log1pexp(-37.0) - (-37.0_f64).exp().ln_1p()).abs() < 1e-10);

        // x <= 18: returns log1p(exp(x))
        assert!((log1pexp(0.0) - (1.0_f64 + 1.0).ln()).abs() < 1e-10);
        assert!((log1pexp(10.0) - (1.0 + 10.0_f64.exp()).ln()).abs() < 1e-10);

        // x <= 33.3: returns x + exp(-x)
        assert!((log1pexp(25.0) - (25.0 + (-25.0_f64).exp())).abs() < 1e-10);

        // x > 33.3: returns x
        assert!((log1pexp(40.0) - 40.0).abs() < 1e-10);
        assert!((log1pexp(100.0) - 100.0).abs() < 1e-10);
    }

    #[test]
    fn test_ln_normalize() {
        // Normalizing ln(0.2) by ln(0.5) should give ln(0.4)
        let result = ln_normalize(0.2_f64.ln(), 0.5_f64.ln());
        assert!((result - 0.4_f64.ln()).abs() < 1e-10);
    }

    // =====================================================================
    // Additional edge case tests
    // =====================================================================

    #[test]
    fn test_phred_boundary_values() {
        // Test boundary Phred scores
        // Q0 = 100% error (probability = 1.0)
        let q0_error = phred_to_ln_error_prob(0);
        assert!((q0_error - 0.0).abs() < 1e-10); // ln(1.0) = 0

        // Q2 = MIN_PHRED
        let q2_error = phred_to_ln_error_prob(MIN_PHRED);
        let expected_q2 = 10.0_f64.powf(-0.2); // 10^(-2/10) ≈ 0.631
        assert!((q2_error.exp() - expected_q2).abs() < 1e-6);

        // Q93 = MAX_PHRED
        let q93_error = phred_to_ln_error_prob(MAX_PHRED);
        let expected_q93 = 10.0_f64.powf(-9.3); // 10^(-93/10) ≈ 5e-10
        assert!((q93_error.exp() - expected_q93).abs() < 1e-15);
    }

    #[test]
    fn test_phred_to_ln_correct_boundary() {
        // Q0: error = 1.0, correct = 0.0 → ln(0) = -∞
        let q0_correct = phred_to_ln_correct_prob(0);
        assert!(q0_correct.is_infinite() && q0_correct < 0.0);

        // Q93: error ≈ 5e-10, correct ≈ 1.0 → ln(correct) ≈ 0
        let q93_correct = phred_to_ln_correct_prob(MAX_PHRED);
        assert!(q93_correct < 0.0 && q93_correct.abs() < 1e-9);
    }

    #[test]
    fn test_ln_prob_to_phred_clamping() {
        // Very high error (low quality) should clamp to MIN_PHRED
        let very_high_error = 0.9_f64.ln(); // 90% error
        assert_eq!(ln_prob_to_phred(very_high_error), MIN_PHRED);

        // 100% error should also clamp to MIN_PHRED
        // ln(1.0) = 0, which gives Q = 0, clamped to MIN_PHRED
        assert_eq!(ln_prob_to_phred(0.0), MIN_PHRED);

        // Very low error should clamp to MAX_PHRED
        let very_low_error = 1e-15_f64.ln();
        assert_eq!(ln_prob_to_phred(very_low_error), MAX_PHRED);
    }

    #[test]
    fn test_ln_sum_exp_array_edge_cases() {
        // Empty array
        let result = ln_sum_exp_array(&[]);
        assert!(result.is_infinite() && result < 0.0);

        // Single element
        let result = ln_sum_exp_array(&[0.5_f64.ln()]);
        assert!((result - 0.5_f64.ln()).abs() < 1e-10);

        // Two elements (should match ln_sum_exp)
        let a = 0.2_f64.ln();
        let b = 0.3_f64.ln();
        let array_result = ln_sum_exp_array(&[a, b]);
        let pair_result = ln_sum_exp(a, b);
        assert!((array_result - pair_result).abs() < 1e-10);

        // All negative infinity
        let result = ln_sum_exp_array(&[f64::NEG_INFINITY, f64::NEG_INFINITY]);
        assert!(result.is_infinite() && result < 0.0);
    }

    #[test]
    fn test_ln_sum_exp_very_small_values() {
        // Very small probabilities (should not underflow)
        let a = 1e-300_f64.ln();
        let b = 2e-300_f64.ln();
        let result = ln_sum_exp(a, b);
        let expected = 3e-300_f64.ln();
        assert!((result - expected).abs() < 1.0); // Large tolerance for extreme values
    }

    #[test]
    fn test_ln_not_edge_cases() {
        // ln_not(ln(0)) = ln(1 - 0) = ln(1) = 0
        let result = ln_not(f64::NEG_INFINITY);
        assert!((result - 0.0).abs() < 1e-10);

        // ln_not(ln(1)) = ln(1 - 1) = ln(0) = -∞
        let result = ln_not(0.0);
        assert!(result.is_infinite() && result < 0.0);

        // ln_not(ln(0.5)) = ln(0.5)
        let result = ln_not(0.5_f64.ln());
        assert!((result - 0.5_f64.ln()).abs() < 1e-10);

        // Positive value (invalid probability) should return -∞
        let result = ln_not(1.0);
        assert!(result.is_infinite() && result < 0.0);
    }

    #[test]
    fn test_error_two_trials_quick_approximation() {
        // When one probability is much smaller (>6 in log space), use quick approximation
        let large_error = 0.5_f64.ln(); // ln(0.5) ≈ -0.693
        let small_error = 1e-6_f64.ln(); // ln(1e-6) ≈ -13.8

        let result = ln_error_prob_two_trials(large_error, small_error);
        // Should approximately equal the larger error
        assert!((result - large_error).abs() < 0.01);

        // Also test with arguments swapped
        let result2 = ln_error_prob_two_trials(small_error, large_error);
        assert!((result2 - large_error).abs() < 0.01);
    }

    #[test]
    fn test_log1pexp_boundary_values() {
        // At threshold -37
        let at_threshold = log1pexp(-37.0);
        assert!(at_threshold > 0.0);
        assert!(at_threshold < 1e-15);

        // At threshold 18
        let at_18 = log1pexp(18.0);
        let expected = (1.0 + 18.0_f64.exp()).ln();
        assert!((at_18 - expected).abs() < 1e-10);

        // At threshold 33.3
        let at_33 = log1pexp(33.3);
        let expected = 33.3 + (-33.3_f64).exp();
        assert!((at_33 - expected).abs() < 1e-10);
    }

    #[test]
    fn test_constants() {
        // Verify constants are correct
        assert!((LN_TWO - std::f64::consts::LN_2).abs() < 1e-15);
        assert!((LN_FOUR_THIRDS - (4.0_f64 / 3.0).ln()).abs() < 1e-15);
        assert!(LN_ONE.abs() < 1e-15); // LN_ONE should be 0.0
        assert_eq!(MIN_PHRED, 2);
        assert_eq!(MAX_PHRED, 93);
        assert_eq!(NO_CALL_BASE, b'N');
        assert_eq!(NO_CALL_BASE_LOWER, b'n');
    }
}

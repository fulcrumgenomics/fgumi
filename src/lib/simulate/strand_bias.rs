//! Strand bias models for duplex sequencing simulation.
//!
//! Models the distribution of reads between A and B strands in duplex
//! sequencing, accounting for PCR amplification bias where one strand
//! may start amplifying in earlier rounds.

use rand::Rng;
use rand_distr::{Beta, Distribution};

/// Beta-distributed A/B strand ratio model for duplex sequencing.
///
/// In duplex sequencing, molecules are sequenced from both strands (A and B).
/// Due to PCR amplification dynamics, the ratio of A to B strand reads is
/// rarely exactly 50/50. This model uses a Beta distribution to simulate
/// realistic strand ratios.
///
/// # Parameters
///
/// - `alpha` and `beta` control the shape of the distribution
/// - When alpha == beta, the distribution is symmetric around 0.5
/// - Higher values concentrate the distribution closer to 0.5
/// - Lower values allow more extreme ratios
///
/// # Examples
///
/// ```
/// use fgumi_lib::simulate::StrandBiasModel;
/// use fgumi_lib::simulate::create_rng;
///
/// // Symmetric distribution centered at 0.5
/// let model = StrandBiasModel::default();
/// let mut rng = create_rng(Some(42));
///
/// // Split 10 reads between A and B strands
/// let (a_count, b_count) = model.split_reads(10, &mut rng);
/// assert_eq!(a_count + b_count, 10);
/// ```
#[derive(Debug, Clone)]
pub struct StrandBiasModel {
    /// Beta distribution alpha parameter
    pub alpha: f64,
    /// Beta distribution beta parameter
    pub beta: f64,
}

impl Default for StrandBiasModel {
    fn default() -> Self {
        // alpha=5, beta=5 gives a symmetric distribution
        // concentrated around 0.5 but with realistic variation
        Self { alpha: 5.0, beta: 5.0 }
    }
}

#[allow(clippy::cast_precision_loss, clippy::cast_possible_truncation)]
impl StrandBiasModel {
    /// Create a new strand bias model with custom parameters.
    ///
    /// # Arguments
    ///
    /// * `alpha` - Beta distribution alpha parameter (must be > 0)
    /// * `beta` - Beta distribution beta parameter (must be > 0)
    #[must_use]
    pub fn new(alpha: f64, beta: f64) -> Self {
        Self { alpha, beta }
    }

    /// Create a model with no strand bias (always 50/50 split).
    ///
    /// Uses very high alpha/beta values to concentrate near 0.5.
    #[must_use]
    pub fn no_bias() -> Self {
        Self { alpha: 1000.0, beta: 1000.0 }
    }

    /// Sample an A-strand fraction from the distribution.
    ///
    /// # Arguments
    ///
    /// * `rng` - Random number generator
    ///
    /// # Returns
    ///
    /// A value in [0, 1] representing the fraction of reads on the A strand.
    pub fn sample_a_fraction(&self, rng: &mut impl Rng) -> f64 {
        let dist = Beta::new(self.alpha, self.beta).expect("Invalid beta parameters");
        dist.sample(rng)
    }

    /// Split a total read count between A and B strands.
    ///
    /// # Arguments
    ///
    /// * `total` - Total number of reads to split
    /// * `rng` - Random number generator
    ///
    /// # Returns
    ///
    /// A tuple of (`a_count`, `b_count`) where `a_count` + `b_count` == total.
    pub fn split_reads(&self, total: usize, rng: &mut impl Rng) -> (usize, usize) {
        if total == 0 {
            return (0, 0);
        }

        let a_fraction = self.sample_a_fraction(rng);
        let a_count = (total as f64 * a_fraction).round() as usize;
        let b_count = total.saturating_sub(a_count);

        (a_count, b_count)
    }

    /// Split a total read count, ensuring both strands have at least `min_per_strand` reads.
    ///
    /// If total < 2 * `min_per_strand`, the constraint cannot be satisfied and
    /// the split is done without the constraint.
    ///
    /// # Arguments
    ///
    /// * `total` - Total number of reads to split
    /// * `min_per_strand` - Minimum reads required per strand
    /// * `rng` - Random number generator
    ///
    /// # Returns
    ///
    /// A tuple of (`a_count`, `b_count`) where `a_count` + `b_count` == total.
    pub fn split_reads_with_minimum(
        &self,
        total: usize,
        min_per_strand: usize,
        rng: &mut impl Rng,
    ) -> (usize, usize) {
        if total < 2 * min_per_strand {
            // Cannot satisfy constraint, fall back to unconstrained split
            return self.split_reads(total, rng);
        }

        // Reserve minimum for each strand, then distribute remainder
        let remainder = total - 2 * min_per_strand;
        let a_fraction = self.sample_a_fraction(rng);
        let a_extra = (remainder as f64 * a_fraction).round() as usize;

        let a_count = min_per_strand + a_extra;
        let b_count = total - a_count;

        (a_count, b_count)
    }
}

#[cfg(test)]
#[allow(clippy::cast_precision_loss)]
mod tests {
    use super::*;
    use crate::simulate::create_rng;

    #[test]
    fn test_default_parameters() {
        let model = StrandBiasModel::default();
        assert!((model.alpha - 5.0).abs() < f64::EPSILON);
        assert!((model.beta - 5.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_a_fraction_in_range() {
        let model = StrandBiasModel::default();
        let mut rng = create_rng(Some(42));

        for _ in 0..100 {
            let frac = model.sample_a_fraction(&mut rng);
            assert!((0.0..=1.0).contains(&frac));
        }
    }

    #[test]
    fn test_symmetric_distribution() {
        let model = StrandBiasModel::new(5.0, 5.0);
        let mut rng = create_rng(Some(42));

        let fractions: Vec<f64> = (0..1000).map(|_| model.sample_a_fraction(&mut rng)).collect();
        let mean: f64 = fractions.iter().sum::<f64>() / fractions.len() as f64;

        // With symmetric alpha=beta, mean should be close to 0.5
        assert!((mean - 0.5).abs() < 0.05);
    }

    #[test]
    fn test_split_reads_sum() {
        let model = StrandBiasModel::default();
        let mut rng = create_rng(Some(42));

        for total in [0, 1, 5, 10, 100] {
            let (a, b) = model.split_reads(total, &mut rng);
            assert_eq!(a + b, total);
        }
    }

    #[test]
    fn test_split_reads_zero() {
        let model = StrandBiasModel::default();
        let mut rng = create_rng(Some(42));

        let (a, b) = model.split_reads(0, &mut rng);
        assert_eq!(a, 0);
        assert_eq!(b, 0);
    }

    #[test]
    fn test_no_bias() {
        let model = StrandBiasModel::no_bias();
        let mut rng = create_rng(Some(42));

        // With very high alpha/beta, fraction should be very close to 0.5
        let fractions: Vec<f64> = (0..100).map(|_| model.sample_a_fraction(&mut rng)).collect();
        let mean: f64 = fractions.iter().sum::<f64>() / fractions.len() as f64;

        assert!((mean - 0.5).abs() < 0.01);
    }

    #[test]
    fn test_split_with_minimum() {
        let model = StrandBiasModel::default();
        let mut rng = create_rng(Some(42));

        // With 10 reads and min 2, both strands should have at least 2
        for _ in 0..100 {
            let (a, b) = model.split_reads_with_minimum(10, 2, &mut rng);
            assert!(a >= 2);
            assert!(b >= 2);
            assert_eq!(a + b, 10);
        }
    }

    #[test]
    fn test_split_with_minimum_unsatisfiable() {
        let model = StrandBiasModel::default();
        let mut rng = create_rng(Some(42));

        // With 3 reads and min 2, constraint cannot be satisfied
        // Should fall back to unconstrained split
        let (a, b) = model.split_reads_with_minimum(3, 2, &mut rng);
        assert_eq!(a + b, 3);
    }

    #[test]
    fn test_reproducibility() {
        let model = StrandBiasModel::default();

        let mut rng1 = create_rng(Some(42));
        let mut rng2 = create_rng(Some(42));

        let fractions1: Vec<f64> = (0..10).map(|_| model.sample_a_fraction(&mut rng1)).collect();
        let fractions2: Vec<f64> = (0..10).map(|_| model.sample_a_fraction(&mut rng2)).collect();

        assert_eq!(fractions1, fractions2);
    }

    #[test]
    fn test_split_reads_one() {
        let model = StrandBiasModel::default();
        let mut rng = create_rng(Some(42));

        // With 1 read, should go to either A or B
        for _ in 0..100 {
            let (a, b) = model.split_reads(1, &mut rng);
            assert_eq!(a + b, 1);
            assert!(a == 0 || a == 1);
        }
    }

    #[test]
    fn test_split_reads_large() {
        let model = StrandBiasModel::new(5.0, 5.0);
        let mut rng = create_rng(Some(42));

        let (a, b) = model.split_reads(10000, &mut rng);
        assert_eq!(a + b, 10000);

        // With symmetric distribution, should be roughly equal
        let frac_a = a as f64 / 10000.0;
        assert!(frac_a > 0.3 && frac_a < 0.7, "Expected roughly balanced split, got {frac_a}");
    }

    #[test]
    fn test_asymmetric_a_heavy() {
        let model = StrandBiasModel::new(10.0, 2.0);
        let mut rng = create_rng(Some(42));

        let splits: Vec<(usize, usize)> =
            (0..100).map(|_| model.split_reads(100, &mut rng)).collect();

        let avg_a: f64 = splits.iter().map(|(a, _)| *a as f64).sum::<f64>() / splits.len() as f64;
        let avg_b: f64 = splits.iter().map(|(_, b)| *b as f64).sum::<f64>() / splits.len() as f64;

        // A-biased, so avg_a should be higher
        assert!(avg_a > avg_b, "Expected A > B, got A={avg_a}, B={avg_b}");
        assert!(avg_a > 70.0, "Expected strong A bias");
    }

    #[test]
    fn test_asymmetric_b_heavy() {
        let model = StrandBiasModel::new(2.0, 10.0);
        let mut rng = create_rng(Some(42));

        let splits: Vec<(usize, usize)> =
            (0..100).map(|_| model.split_reads(100, &mut rng)).collect();

        let avg_a: f64 = splits.iter().map(|(a, _)| *a as f64).sum::<f64>() / splits.len() as f64;
        let avg_b: f64 = splits.iter().map(|(_, b)| *b as f64).sum::<f64>() / splits.len() as f64;

        // B-biased, so avg_b should be higher
        assert!(avg_b > avg_a, "Expected B > A, got A={avg_a}, B={avg_b}");
        assert!(avg_b > 70.0, "Expected strong B bias");
    }

    #[test]
    fn test_split_with_minimum_zero() {
        let model = StrandBiasModel::default();
        let mut rng = create_rng(Some(42));

        // min_per_strand = 0 is same as unconstrained
        for _ in 0..100 {
            let (a, b) = model.split_reads_with_minimum(10, 0, &mut rng);
            assert_eq!(a + b, 10);
        }
    }

    #[test]
    fn test_split_with_minimum_half() {
        let model = StrandBiasModel::default();
        let mut rng = create_rng(Some(42));

        // min_per_strand = 5 with total = 10 means exactly 5/5
        for _ in 0..100 {
            let (a, b) = model.split_reads_with_minimum(10, 5, &mut rng);
            assert_eq!(a, 5);
            assert_eq!(b, 5);
        }
    }

    #[test]
    fn test_split_with_minimum_respects_bounds() {
        let model = StrandBiasModel::new(8.0, 2.0); // Heavily A-biased
        let mut rng = create_rng(Some(42));

        // Even with heavy A bias, should respect minimum of 3 per strand
        for _ in 0..100 {
            let (a, b) = model.split_reads_with_minimum(20, 3, &mut rng);
            assert!(a >= 3, "A={a} should be >= 3");
            assert!(b >= 3, "B={b} should be >= 3");
            assert_eq!(a + b, 20);
        }
    }

    #[test]
    fn test_high_concentration() {
        // Very high alpha/beta should concentrate around 0.5
        let model = StrandBiasModel::new(100.0, 100.0);
        let mut rng = create_rng(Some(42));

        let fractions: Vec<f64> = (0..100).map(|_| model.sample_a_fraction(&mut rng)).collect();

        // Mean should be very close to 0.5
        let mean: f64 = fractions.iter().sum::<f64>() / fractions.len() as f64;
        assert!((mean - 0.5).abs() < 0.05, "Mean {mean} too far from 0.5 with high concentration");

        // Most fractions should be within 0.15 of 0.5
        let close_count = fractions.iter().filter(|&&f| (f - 0.5).abs() < 0.15).count();
        assert!(close_count > 90, "Only {close_count}/100 fractions within 0.15 of 0.5");
    }

    #[test]
    fn test_low_concentration() {
        // Very low alpha/beta should give extreme values
        let model = StrandBiasModel::new(0.5, 0.5);
        let mut rng = create_rng(Some(42));

        let fractions: Vec<f64> = (0..1000).map(|_| model.sample_a_fraction(&mut rng)).collect();

        // Should have some extreme values
        let extreme_low = fractions.iter().filter(|&&f| f < 0.2).count();
        let extreme_high = fractions.iter().filter(|&&f| f > 0.8).count();

        assert!(extreme_low > 100, "Expected more extreme low values, got {extreme_low}");
        assert!(extreme_high > 100, "Expected more extreme high values, got {extreme_high}");
    }

    #[test]
    fn test_split_reads_consistency() {
        let model = StrandBiasModel::default();
        let mut rng = create_rng(Some(42));

        // Verify that splits always sum to total
        for total in 0..50 {
            let (a, b) = model.split_reads(total, &mut rng);
            assert_eq!(a + b, total, "Split ({a}, {b}) doesn't sum to {total}");
        }
    }

    #[test]
    fn test_a_fraction_distribution_shape() {
        let model = StrandBiasModel::new(2.0, 5.0);
        let mut rng = create_rng(Some(42));

        let fractions: Vec<f64> = (0..10000).map(|_| model.sample_a_fraction(&mut rng)).collect();

        // Expected mean of Beta(2, 5) is 2/(2+5) = 0.286
        let mean: f64 = fractions.iter().sum::<f64>() / fractions.len() as f64;
        let expected_mean = 2.0 / 7.0;

        assert!(
            (mean - expected_mean).abs() < 0.02,
            "Mean {mean} not close to expected {expected_mean}"
        );
    }
}

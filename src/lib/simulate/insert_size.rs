//! Insert size distribution models for paired-end simulation.
//!
//! Provides log-normal insert size distribution to match real library prep
//! characteristics with realistic asymmetric distribution.

use rand::Rng;
use rand_distr::{Distribution, LogNormal};

/// Log-normal insert size distribution model.
///
/// Real paired-end sequencing libraries have insert size distributions that are
/// approximately log-normal with:
/// - Sharp left tail (minimum viable fragment size)
/// - Extended right tail (larger fragments)
///
/// # Examples
///
/// ```
/// use fgumi_lib::simulate::InsertSizeModel;
/// use fgumi_lib::simulate::create_rng;
///
/// let model = InsertSizeModel::default();
/// let mut rng = create_rng(Some(42));
///
/// // Sample insert sizes
/// let sizes: Vec<usize> = (0..100).map(|_| model.sample(&mut rng)).collect();
///
/// // Check that sizes are within bounds
/// for size in &sizes {
///     assert!(*size >= model.min);
///     assert!(*size <= model.max);
/// }
/// ```
#[derive(Debug, Clone)]
pub struct InsertSizeModel {
    /// Mean insert size
    pub mean: f64,
    /// Standard deviation of insert size
    pub stddev: f64,
    /// Minimum insert size (hard floor)
    pub min: usize,
    /// Maximum insert size (hard ceiling)
    pub max: usize,
}

impl Default for InsertSizeModel {
    fn default() -> Self {
        Self { mean: 300.0, stddev: 50.0, min: 50, max: 800 }
    }
}

#[allow(clippy::cast_precision_loss, clippy::cast_possible_truncation, clippy::cast_sign_loss)]
impl InsertSizeModel {
    /// Create a new insert size model with custom parameters.
    ///
    /// # Arguments
    ///
    /// * `mean` - Mean insert size
    /// * `stddev` - Standard deviation
    /// * `min` - Minimum insert size
    /// * `max` - Maximum insert size
    #[must_use]
    pub fn new(mean: f64, stddev: f64, min: usize, max: usize) -> Self {
        Self { mean, stddev, min, max }
    }

    /// Sample an insert size from the distribution.
    ///
    /// Uses log-normal distribution with parameters derived from the
    /// specified mean and standard deviation, then clamps to [min, max].
    ///
    /// # Arguments
    ///
    /// * `rng` - Random number generator
    ///
    /// # Returns
    ///
    /// An insert size clamped to [min, max].
    pub fn sample(&self, rng: &mut impl Rng) -> usize {
        // Convert mean/stddev to log-normal parameters
        // For log-normal: E[X] = exp(mu + sigma^2/2), Var[X] = (exp(sigma^2) - 1) * exp(2*mu + sigma^2)
        // Solving for mu and sigma:
        let variance = self.stddev.powi(2);
        let mean_sq = self.mean.powi(2);
        let sigma_sq = (1.0 + variance / mean_sq).ln();
        let sigma = sigma_sq.sqrt();
        let mu = self.mean.ln() - sigma_sq / 2.0;

        let dist = LogNormal::new(mu, sigma).expect("Invalid log-normal parameters");
        let sample = dist.sample(rng);

        // Clamp to valid range
        (sample.round() as usize).clamp(self.min, self.max)
    }

    /// Sample multiple insert sizes.
    ///
    /// # Arguments
    ///
    /// * `count` - Number of samples to generate
    /// * `rng` - Random number generator
    ///
    /// # Returns
    ///
    /// Vector of insert sizes.
    pub fn sample_n(&self, count: usize, rng: &mut impl Rng) -> Vec<usize> {
        (0..count).map(|_| self.sample(rng)).collect()
    }
}

#[cfg(test)]
#[allow(clippy::cast_precision_loss)]
mod tests {
    use super::*;
    use crate::simulate::create_rng;

    #[test]
    fn test_default_values() {
        let model = InsertSizeModel::default();
        assert!((model.mean - 300.0).abs() < f64::EPSILON);
        assert!((model.stddev - 50.0).abs() < f64::EPSILON);
        assert_eq!(model.min, 50);
        assert_eq!(model.max, 800);
    }

    #[test]
    fn test_samples_within_bounds() {
        let model = InsertSizeModel::default();
        let mut rng = create_rng(Some(42));

        for _ in 0..1000 {
            let size = model.sample(&mut rng);
            assert!(size >= model.min, "Size {size} below min {}", model.min);
            assert!(size <= model.max, "Size {size} above max {}", model.max);
        }
    }

    #[test]
    fn test_distribution_centered_around_mean() {
        let model = InsertSizeModel::new(300.0, 50.0, 50, 800);
        let mut rng = create_rng(Some(42));

        let samples: Vec<usize> = (0..10000).map(|_| model.sample(&mut rng)).collect();
        let mean: f64 = samples.iter().map(|&x| x as f64).sum::<f64>() / samples.len() as f64;

        // Mean should be reasonably close to target (within 10%)
        assert!((mean - 300.0).abs() < 30.0, "Sample mean {mean} too far from expected 300.0");
    }

    #[test]
    fn test_custom_model() {
        let model = InsertSizeModel::new(200.0, 30.0, 100, 400);
        let mut rng = create_rng(Some(42));

        let size = model.sample(&mut rng);
        assert!(size >= 100);
        assert!(size <= 400);
    }

    #[test]
    fn test_sample_n() {
        let model = InsertSizeModel::default();
        let mut rng = create_rng(Some(42));

        let samples = model.sample_n(100, &mut rng);
        assert_eq!(samples.len(), 100);
    }

    #[test]
    fn test_reproducibility() {
        let model = InsertSizeModel::default();

        let mut rng1 = create_rng(Some(42));
        let mut rng2 = create_rng(Some(42));

        let samples1: Vec<usize> = (0..10).map(|_| model.sample(&mut rng1)).collect();
        let samples2: Vec<usize> = (0..10).map(|_| model.sample(&mut rng2)).collect();

        assert_eq!(samples1, samples2);
    }

    #[test]
    fn test_very_narrow_range() {
        // Model where min == max
        let model = InsertSizeModel::new(100.0, 10.0, 100, 100);
        let mut rng = create_rng(Some(42));

        for _ in 0..100 {
            let size = model.sample(&mut rng);
            assert_eq!(size, 100);
        }
    }

    #[test]
    fn test_high_variance() {
        let model = InsertSizeModel::new(300.0, 150.0, 50, 1000);
        let mut rng = create_rng(Some(42));

        let samples: Vec<usize> = (0..1000).map(|_| model.sample(&mut rng)).collect();

        // Check that we have variation
        let min = *samples.iter().min().unwrap();
        let max = *samples.iter().max().unwrap();
        assert!(max - min > 200, "Expected high variance but got narrow range");
    }

    #[test]
    fn test_low_variance() {
        let model = InsertSizeModel::new(300.0, 5.0, 250, 350);
        let mut rng = create_rng(Some(42));

        let samples: Vec<usize> = (0..1000).map(|_| model.sample(&mut rng)).collect();

        // With low stddev, most samples should be close to mean
        let count_near_mean = samples.iter().filter(|&&s| (290..=310).contains(&s)).count();
        assert!(count_near_mean > 800, "Expected most samples near mean with low variance");
    }

    #[test]
    fn test_asymmetric_clamping() {
        // Mean is near min, so clamping should be more frequent on left
        let model = InsertSizeModel::new(100.0, 50.0, 75, 200);
        let mut rng = create_rng(Some(42));

        let samples: Vec<usize> = (0..1000).map(|_| model.sample(&mut rng)).collect();
        let at_min = samples.iter().filter(|&&s| s == 75).count();

        // Should have some samples clamped to min
        assert!(at_min > 0, "Expected some samples clamped to min");
    }

    #[test]
    fn test_sample_n_empty() {
        let model = InsertSizeModel::default();
        let mut rng = create_rng(Some(42));

        let samples = model.sample_n(0, &mut rng);
        assert!(samples.is_empty());
    }

    #[test]
    fn test_sample_n_large() {
        let model = InsertSizeModel::default();
        let mut rng = create_rng(Some(42));

        let samples = model.sample_n(10000, &mut rng);
        assert_eq!(samples.len(), 10000);

        // All should be within bounds
        for &size in &samples {
            assert!(size >= model.min && size <= model.max);
        }
    }

    #[test]
    fn test_min_greater_than_mean() {
        // Edge case: min > mean (unusual but should clamp)
        let model = InsertSizeModel::new(100.0, 20.0, 150, 200);
        let mut rng = create_rng(Some(42));

        for _ in 0..100 {
            let size = model.sample(&mut rng);
            assert!(size >= 150);
        }
    }

    #[test]
    fn test_statistics() {
        let model = InsertSizeModel::new(300.0, 50.0, 50, 800);
        let mut rng = create_rng(Some(42));

        let samples: Vec<usize> = (0..50000).map(|_| model.sample(&mut rng)).collect();
        let mean: f64 = samples.iter().map(|&x| x as f64).sum::<f64>() / samples.len() as f64;
        let variance: f64 =
            samples.iter().map(|&x| (x as f64 - mean).powi(2)).sum::<f64>() / samples.len() as f64;
        let stddev = variance.sqrt();

        // Check mean is close to target (within 5%)
        assert!((mean - 300.0).abs() < 15.0, "Mean {mean} not close enough to 300");

        // Check stddev is roughly in range (log-normal will have slightly different stddev)
        assert!(stddev > 30.0 && stddev < 80.0, "Stddev {stddev} out of expected range");
    }
}

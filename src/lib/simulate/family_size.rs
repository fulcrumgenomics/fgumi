//! Family size distribution models for UMI simulation.
//!
//! Provides multiple distribution types for modeling the number of reads
//! per UMI family (molecule):
//! - Log-normal (default)
//! - Negative binomial (alternative PCR model)
//! - Empirical (loaded from `fgumi group -f` histogram output)

use anyhow::{Context, Result};
use rand::{Rng, RngExt};
use rand_distr::{Distribution, Gamma, LogNormal, Poisson};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Family size distribution for UMI simulation.
///
/// Determines how many reads are generated per molecule (UMI family).
///
/// # Variants
///
/// - `LogNormal` - Natural for PCR amplification, parameterized by mean and stddev
/// - `NegativeBinomial` - Alternative PCR model, parameterized by r and p
/// - `Empirical` - Loaded from histogram file (e.g., `fgumi group -f` output)
///
/// # Examples
///
/// ```
/// use fgumi_lib::simulate::FamilySizeDistribution;
/// use fgumi_lib::simulate::create_rng;
///
/// let dist = FamilySizeDistribution::log_normal(3.0, 2.0);
/// let mut rng = create_rng(Some(42));
///
/// let family_size = dist.sample(&mut rng, 1);
/// assert!(family_size >= 1);
/// ```
#[derive(Debug, Clone)]
pub enum FamilySizeDistribution {
    /// Log-normal distribution (natural for PCR)
    LogNormal {
        /// Mean family size
        mean: f64,
        /// Standard deviation
        stddev: f64,
    },
    /// Negative binomial distribution
    NegativeBinomial {
        /// Number of successes (shape parameter)
        r: f64,
        /// Probability of success
        p: f64,
    },
    /// Empirical distribution from histogram
    Empirical {
        /// Cumulative distribution: (`family_size`, `cumulative_probability`)
        cdf: Vec<(usize, f64)>,
    },
}

impl Default for FamilySizeDistribution {
    fn default() -> Self {
        Self::LogNormal { mean: 3.0, stddev: 2.0 }
    }
}

#[allow(clippy::cast_precision_loss, clippy::cast_possible_truncation, clippy::cast_sign_loss)]
impl FamilySizeDistribution {
    /// Create a log-normal family size distribution.
    ///
    /// # Arguments
    ///
    /// * `mean` - Mean family size
    /// * `stddev` - Standard deviation
    #[must_use]
    pub fn log_normal(mean: f64, stddev: f64) -> Self {
        Self::LogNormal { mean, stddev }
    }

    /// Create a negative binomial family size distribution.
    ///
    /// # Arguments
    ///
    /// * `r` - Number of successes (shape parameter)
    /// * `p` - Probability of success
    #[must_use]
    pub fn negative_binomial(r: f64, p: f64) -> Self {
        Self::NegativeBinomial { r, p }
    }

    /// Load an empirical distribution from a histogram file.
    ///
    /// The file should be TSV format with columns:
    /// - `family_size` (required)
    /// - `count` (required)
    /// - Other columns are ignored
    ///
    /// This is compatible with `fgumi group -f` output.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the histogram TSV file
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be read or parsed.
    pub fn from_histogram<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        let file = File::open(path)
            .with_context(|| format!("Failed to open histogram file: {}", path.display()))?;
        let reader = BufReader::new(file);

        let mut counts: Vec<(usize, u64)> = Vec::new();
        let mut header_seen = false;
        let mut family_size_col = 0;
        let mut count_col = 1;

        for line_result in reader.lines() {
            let line = line_result.context("Failed to read line")?;
            let fields: Vec<&str> = line.split('\t').collect();

            if !header_seen {
                // Parse header to find column indices
                for (i, field) in fields.iter().enumerate() {
                    if *field == "family_size" {
                        family_size_col = i;
                    } else if *field == "count" {
                        count_col = i;
                    }
                }
                header_seen = true;
                continue;
            }

            if fields.len() <= family_size_col.max(count_col) {
                continue; // Skip malformed lines
            }

            let family_size: usize = fields[family_size_col]
                .parse()
                .with_context(|| format!("Invalid family_size: {}", fields[family_size_col]))?;
            let count: u64 = fields[count_col]
                .parse()
                .with_context(|| format!("Invalid count: {}", fields[count_col]))?;

            counts.push((family_size, count));
        }

        if counts.is_empty() {
            anyhow::bail!("No data found in histogram file");
        }

        // Sort by family size
        counts.sort_by_key(|(size, _)| *size);

        // Convert to CDF
        let total: u64 = counts.iter().map(|(_, c)| c).sum();
        let mut cumulative = 0u64;
        let cdf: Vec<(usize, f64)> = counts
            .into_iter()
            .map(|(size, count)| {
                cumulative += count;
                (size, cumulative as f64 / total as f64)
            })
            .collect();

        Ok(Self::Empirical { cdf })
    }

    /// Sample a family size from the distribution.
    ///
    /// # Arguments
    ///
    /// * `rng` - Random number generator
    /// * `min_size` - Minimum family size (samples below this are clamped)
    ///
    /// # Returns
    ///
    /// A family size >= `min_size`.
    pub fn sample(&self, rng: &mut impl Rng, min_size: usize) -> usize {
        let size = match self {
            Self::LogNormal { mean, stddev } => {
                // Convert mean/stddev to log-normal parameters
                let variance = stddev.powi(2);
                let mean_sq = mean.powi(2);
                let sigma_sq = (1.0 + variance / mean_sq).ln();
                let sigma = sigma_sq.sqrt();
                let mu = mean.ln() - sigma_sq / 2.0;

                let dist = LogNormal::new(mu, sigma).expect("Invalid log-normal parameters");
                dist.sample(rng).round() as usize
            }
            Self::NegativeBinomial { r, p } => {
                // Negative binomial via Gamma-Poisson mixture:
                // If X ~ Gamma(r, p/(1-p)) and Y|X ~ Poisson(X), then Y ~ NegBinomial(r, p)
                let scale = *p / (1.0 - *p);
                let gamma = Gamma::new(*r, scale).expect("Invalid gamma parameters");
                let lambda = gamma.sample(rng);
                if lambda <= 0.0 {
                    0
                } else {
                    let poisson = Poisson::new(lambda).expect("Invalid poisson parameter");
                    poisson.sample(rng) as usize
                }
            }
            Self::Empirical { cdf } => {
                let u: f64 = rng.random();
                cdf.iter().find(|(_, cum_prob)| *cum_prob >= u).map_or(1, |(size, _)| *size)
            }
        };

        size.max(min_size)
    }

    /// Sample multiple family sizes.
    ///
    /// # Arguments
    ///
    /// * `count` - Number of samples
    /// * `rng` - Random number generator
    /// * `min_size` - Minimum family size
    ///
    /// # Returns
    ///
    /// Vector of family sizes.
    pub fn sample_n(&self, count: usize, rng: &mut impl Rng, min_size: usize) -> Vec<usize> {
        (0..count).map(|_| self.sample(rng, min_size)).collect()
    }
}

#[cfg(test)]
#[allow(clippy::cast_precision_loss)]
mod tests {
    use super::*;
    use crate::simulate::create_rng;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_log_normal_default() {
        let dist = FamilySizeDistribution::default();
        let mut rng = create_rng(Some(42));

        for _ in 0..100 {
            let size = dist.sample(&mut rng, 1);
            assert!(size >= 1);
        }
    }

    #[test]
    fn test_log_normal_min_size() {
        let dist = FamilySizeDistribution::log_normal(1.5, 0.5);
        let mut rng = create_rng(Some(42));

        for _ in 0..100 {
            let size = dist.sample(&mut rng, 3);
            assert!(size >= 3);
        }
    }

    #[test]
    fn test_negative_binomial() {
        let dist = FamilySizeDistribution::negative_binomial(2.0, 0.5);
        let mut rng = create_rng(Some(42));

        for _ in 0..100 {
            let size = dist.sample(&mut rng, 1);
            assert!(size >= 1);
        }
    }

    #[test]
    fn test_empirical_from_histogram() -> Result<()> {
        let mut temp = NamedTempFile::new()?;
        writeln!(temp, "family_size\tcount\tfraction")?;
        writeln!(temp, "1\t100\t0.5")?;
        writeln!(temp, "2\t50\t0.25")?;
        writeln!(temp, "3\t30\t0.15")?;
        writeln!(temp, "4\t20\t0.10")?;

        let dist = FamilySizeDistribution::from_histogram(temp.path())?;
        let mut rng = create_rng(Some(42));

        // Sample many times and check distribution
        let samples: Vec<usize> = (0..1000).map(|_| dist.sample(&mut rng, 1)).collect();

        // All samples should be in range 1-4
        for &size in &samples {
            assert!((1..=4).contains(&size));
        }

        // Size 1 should be most common (50% of weight)
        let count_1 = samples.iter().filter(|&&s| s == 1).count();
        assert!(count_1 > 400 && count_1 < 600); // Roughly 50% ± 10%

        Ok(())
    }

    #[test]
    fn test_sample_n() {
        let dist = FamilySizeDistribution::default();
        let mut rng = create_rng(Some(42));

        let samples = dist.sample_n(100, &mut rng, 1);
        assert_eq!(samples.len(), 100);
    }

    #[test]
    fn test_reproducibility() {
        let dist = FamilySizeDistribution::log_normal(3.0, 2.0);

        let mut rng1 = create_rng(Some(42));
        let mut rng2 = create_rng(Some(42));

        let samples1: Vec<usize> = (0..10).map(|_| dist.sample(&mut rng1, 1)).collect();
        let samples2: Vec<usize> = (0..10).map(|_| dist.sample(&mut rng2, 1)).collect();

        assert_eq!(samples1, samples2);
    }

    #[test]
    fn test_log_normal_high_variance() {
        let dist = FamilySizeDistribution::log_normal(5.0, 10.0);
        let mut rng = create_rng(Some(42));

        let samples: Vec<usize> = (0..1000).map(|_| dist.sample(&mut rng, 1)).collect();

        // With high variance, we should see a range of sizes
        let max = *samples.iter().max().unwrap();
        assert!(max > 10, "Expected some large family sizes with high variance");
    }

    #[test]
    fn test_log_normal_low_variance() {
        let dist = FamilySizeDistribution::log_normal(3.0, 0.5);
        let mut rng = create_rng(Some(42));

        let samples: Vec<usize> = (0..1000).map(|_| dist.sample(&mut rng, 1)).collect();

        // With low variance, most should be close to mean
        let count_near_mean = samples.iter().filter(|&&s| (2..=5).contains(&s)).count();
        assert!(count_near_mean > 700, "Expected most samples near mean");
    }

    #[test]
    fn test_negative_binomial_different_params() {
        let dist1 = FamilySizeDistribution::negative_binomial(1.0, 0.5);
        let dist2 = FamilySizeDistribution::negative_binomial(5.0, 0.5);
        let mut rng1 = create_rng(Some(42));
        let mut rng2 = create_rng(Some(42));

        let samples1: Vec<usize> = (0..1000).map(|_| dist1.sample(&mut rng1, 1)).collect();
        let samples2: Vec<usize> = (0..1000).map(|_| dist2.sample(&mut rng2, 1)).collect();

        let mean1: f64 = samples1.iter().map(|&x| x as f64).sum::<f64>() / samples1.len() as f64;
        let mean2: f64 = samples2.iter().map(|&x| x as f64).sum::<f64>() / samples2.len() as f64;

        // Higher r should give higher mean
        assert!(mean2 > mean1);
    }

    #[test]
    fn test_negative_binomial_high_p() {
        // In the Gamma-Poisson mixture, scale = p/(1-p)
        // With high p, scale is large, so lambda is large, so family sizes are large
        let dist = FamilySizeDistribution::negative_binomial(2.0, 0.9);
        let mut rng = create_rng(Some(42));

        let samples: Vec<usize> = (0..1000).map(|_| dist.sample(&mut rng, 1)).collect();
        let mean: f64 = samples.iter().map(|&x| x as f64).sum::<f64>() / samples.len() as f64;

        // Expected mean: r * p / (1-p) = 2 * 0.9 / 0.1 = 18
        // Allow some tolerance for sampling variance
        assert!(mean > 10.0 && mean < 30.0, "Mean {mean} not close to expected 18 for p=0.9");
    }

    #[test]
    fn test_empirical_single_value() -> Result<()> {
        let mut temp = NamedTempFile::new()?;
        writeln!(temp, "family_size\tcount")?;
        writeln!(temp, "5\t100")?;

        let dist = FamilySizeDistribution::from_histogram(temp.path())?;
        let mut rng = create_rng(Some(42));

        // All samples should be 5
        for _ in 0..100 {
            let size = dist.sample(&mut rng, 1);
            assert_eq!(size, 5);
        }

        Ok(())
    }

    #[test]
    fn test_empirical_with_extra_columns() -> Result<()> {
        let mut temp = NamedTempFile::new()?;
        writeln!(temp, "extra_col\tfamily_size\tcount\tanother")?;
        writeln!(temp, "foo\t1\t50\tbar")?;
        writeln!(temp, "foo\t2\t50\tbar")?;

        let dist = FamilySizeDistribution::from_histogram(temp.path())?;
        let mut rng = create_rng(Some(42));

        // Should still work with extra columns
        for _ in 0..100 {
            let size = dist.sample(&mut rng, 1);
            assert!(size == 1 || size == 2);
        }

        Ok(())
    }

    #[test]
    fn test_empirical_unsorted_input() -> Result<()> {
        let mut temp = NamedTempFile::new()?;
        writeln!(temp, "family_size\tcount")?;
        writeln!(temp, "5\t20")?;
        writeln!(temp, "1\t50")?; // Out of order
        writeln!(temp, "3\t30")?;

        let dist = FamilySizeDistribution::from_histogram(temp.path())?;
        let mut rng = create_rng(Some(42));

        // Should still sample correctly
        let samples: Vec<usize> = (0..1000).map(|_| dist.sample(&mut rng, 1)).collect();

        // All should be in range 1, 3, or 5
        for &size in &samples {
            assert!(size == 1 || size == 3 || size == 5);
        }

        Ok(())
    }

    #[test]
    fn test_empirical_empty_file() {
        let temp = NamedTempFile::new().unwrap();

        let result = FamilySizeDistribution::from_histogram(temp.path());
        assert!(result.is_err());
    }

    #[test]
    fn test_empirical_header_only() {
        let mut temp = NamedTempFile::new().unwrap();
        writeln!(temp, "family_size\tcount").unwrap();

        let result = FamilySizeDistribution::from_histogram(temp.path());
        assert!(result.is_err());
    }

    #[test]
    fn test_min_size_larger_than_samples() {
        let dist = FamilySizeDistribution::log_normal(1.5, 0.5);
        let mut rng = create_rng(Some(42));

        // min_size of 10 should always return at least 10
        for _ in 0..100 {
            let size = dist.sample(&mut rng, 10);
            assert!(size >= 10);
        }
    }

    #[test]
    fn test_sample_n_empty() {
        let dist = FamilySizeDistribution::default();
        let mut rng = create_rng(Some(42));

        let samples = dist.sample_n(0, &mut rng, 1);
        assert!(samples.is_empty());
    }

    #[test]
    fn test_sample_n_large() {
        let dist = FamilySizeDistribution::log_normal(3.0, 2.0);
        let mut rng = create_rng(Some(42));

        let samples = dist.sample_n(10000, &mut rng, 1);
        assert_eq!(samples.len(), 10000);

        // All should respect min_size
        for &size in &samples {
            assert!(size >= 1);
        }
    }

    #[test]
    fn test_log_normal_statistics() {
        let dist = FamilySizeDistribution::log_normal(3.0, 2.0);
        let mut rng = create_rng(Some(42));

        let samples: Vec<usize> = (0..50000).map(|_| dist.sample(&mut rng, 1)).collect();
        let mean: f64 = samples.iter().map(|&x| x as f64).sum::<f64>() / samples.len() as f64;

        // Mean should be close to target (log-normal can have different behavior due to rounding)
        assert!(
            mean > 2.0 && mean < 4.5,
            "Mean {mean} not in expected range for log_normal(3.0, 2.0)"
        );
    }

    #[test]
    fn test_empirical_cdf_order() -> Result<()> {
        let mut temp = NamedTempFile::new()?;
        writeln!(temp, "family_size\tcount")?;
        writeln!(temp, "1\t100")?; // 50%
        writeln!(temp, "2\t80")?; // 40%
        writeln!(temp, "10\t20")?; // 10%

        let dist = FamilySizeDistribution::from_histogram(temp.path())?;
        let mut rng = create_rng(Some(42));

        let samples: Vec<usize> = (0..10000).map(|_| dist.sample(&mut rng, 1)).collect();

        let count_1 = samples.iter().filter(|&&s| s == 1).count();
        let count_2 = samples.iter().filter(|&&s| s == 2).count();
        let count_10 = samples.iter().filter(|&&s| s == 10).count();

        // Verify proportions roughly match (with tolerance)
        let frac_1 = count_1 as f64 / samples.len() as f64;
        let frac_2 = count_2 as f64 / samples.len() as f64;
        let frac_10 = count_10 as f64 / samples.len() as f64;

        assert!((frac_1 - 0.5).abs() < 0.05, "frac_1 = {frac_1}");
        assert!((frac_2 - 0.4).abs() < 0.05, "frac_2 = {frac_2}");
        assert!((frac_10 - 0.1).abs() < 0.03, "frac_10 = {frac_10}");

        Ok(())
    }
}

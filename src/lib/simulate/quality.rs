//! Position-dependent quality score models for simulation.
//!
//! Models realistic Illumina-style quality score patterns:
//! - Warmup phase: Quality ramps up from initial value to peak
//! - Peak phase: Quality stays at maximum
//! - Decay phase: Quality degrades toward end of read

use rand::Rng;
use rand_distr::{Distribution, Normal};

/// Position-dependent quality model with warmup, peak, and decay phases.
///
/// Models realistic Illumina sequencing quality patterns where:
/// - Early bases have lower quality (phasing warmup)
/// - Middle bases have peak quality
/// - Late bases show quality degradation (phasing/prephasing accumulation)
///
/// # Examples
///
/// ```
/// use fgumi_lib::simulate::PositionQualityModel;
/// use fgumi_lib::simulate::create_rng;
///
/// let model = PositionQualityModel::default();
/// let mut rng = create_rng(Some(42));
///
/// // Generate quality for a 150bp read
/// let qualities: Vec<u8> = (0..150)
///     .map(|pos| model.quality_at(pos, &mut rng))
///     .collect();
///
/// // Early positions should be lower than peak
/// assert!(qualities[0] < qualities[50]);
/// // Late positions should be lower than peak
/// assert!(qualities[149] < qualities[50]);
/// ```
#[derive(Debug, Clone)]
pub struct PositionQualityModel {
    /// Number of bases in warmup phase before peak quality
    pub warmup_bases: usize,
    /// Starting quality score during warmup phase
    pub warmup_start: u8,
    /// Peak quality score (achieved after warmup)
    pub peak_quality: u8,
    /// Position where quality decay begins
    pub decay_start: usize,
    /// Quality drop per base after decay starts
    pub decay_rate: f64,
    /// Minimum quality floor
    pub min_quality: u8,
    /// Standard deviation of Gaussian noise added to quality
    pub noise_stddev: f64,
}

impl Default for PositionQualityModel {
    fn default() -> Self {
        Self {
            warmup_bases: 10,
            warmup_start: 25,
            peak_quality: 37,
            decay_start: 100,
            decay_rate: 0.08,
            min_quality: 2,
            noise_stddev: 2.0,
        }
    }
}

#[allow(clippy::cast_precision_loss, clippy::cast_possible_truncation, clippy::cast_sign_loss)]
impl PositionQualityModel {
    /// Create a new quality model with custom parameters.
    #[must_use]
    pub fn new(
        warmup_bases: usize,
        warmup_start: u8,
        peak_quality: u8,
        decay_start: usize,
        decay_rate: f64,
        min_quality: u8,
        noise_stddev: f64,
    ) -> Self {
        Self {
            warmup_bases,
            warmup_start,
            peak_quality,
            decay_start,
            decay_rate,
            min_quality,
            noise_stddev,
        }
    }

    /// Generate a quality score for the given read position.
    ///
    /// # Arguments
    ///
    /// * `position` - 0-based position within the read
    /// * `rng` - Random number generator for noise
    ///
    /// # Returns
    ///
    /// A quality score clamped to [`min_quality`, 41].
    ///
    /// # Panics
    ///
    /// Panics if `noise_stddev` produces an invalid normal distribution.
    pub fn quality_at(&self, position: usize, rng: &mut impl Rng) -> u8 {
        let base_q = if position < self.warmup_bases {
            // Linear ramp from warmup_start to peak
            let frac = position as f64 / self.warmup_bases as f64;
            f64::from(self.warmup_start) + frac * f64::from(self.peak_quality - self.warmup_start)
        } else if position < self.decay_start {
            // Peak quality region
            f64::from(self.peak_quality)
        } else {
            // Linear decay
            let decay = (position - self.decay_start) as f64 * self.decay_rate;
            (f64::from(self.peak_quality) - decay).max(f64::from(self.min_quality))
        };

        // Add Gaussian noise
        let noise_dist = Normal::new(0.0, self.noise_stddev).unwrap();
        let noisy = base_q + noise_dist.sample(rng);

        // Clamp to valid range [min_quality, 41]
        noisy.round().clamp(f64::from(self.min_quality), 41.0) as u8
    }

    /// Generate quality scores for an entire read.
    ///
    /// # Arguments
    ///
    /// * `read_length` - Length of the read
    /// * `rng` - Random number generator
    ///
    /// # Returns
    ///
    /// Vector of quality scores, one per base.
    pub fn generate_qualities(&self, read_length: usize, rng: &mut impl Rng) -> Vec<u8> {
        (0..read_length).map(|pos| self.quality_at(pos, rng)).collect()
    }
}

/// R1/R2 quality bias model.
///
/// Models the empirically observed phenomenon that R2 reads tend to have
/// slightly lower quality scores than R1 reads.
///
/// # Examples
///
/// ```
/// use fgumi_lib::simulate::ReadPairQualityBias;
///
/// let bias = ReadPairQualityBias::default();
///
/// // Apply to R1 (no change)
/// assert_eq!(bias.apply(30, false), 30);
///
/// // Apply to R2 (offset applied)
/// assert_eq!(bias.apply(30, true), 28);
/// ```
#[derive(Debug, Clone)]
pub struct ReadPairQualityBias {
    /// Quality offset for R2 reads (typically negative)
    pub r2_offset: i8,
}

impl Default for ReadPairQualityBias {
    fn default() -> Self {
        Self { r2_offset: -2 }
    }
}

impl ReadPairQualityBias {
    /// Create a new bias model with custom R2 offset.
    #[must_use]
    pub fn new(r2_offset: i8) -> Self {
        Self { r2_offset }
    }

    /// Apply the bias to a quality score.
    ///
    /// # Arguments
    ///
    /// * `quality` - Original quality score
    /// * `is_r2` - Whether this is an R2 read
    ///
    /// # Returns
    ///
    /// Adjusted quality score, clamped to [2, 41].
    #[must_use]
    #[allow(clippy::cast_sign_loss)]
    pub fn apply(&self, quality: u8, is_r2: bool) -> u8 {
        if is_r2 {
            let adjusted = i16::from(quality) + i16::from(self.r2_offset);
            adjusted.clamp(2, 41) as u8
        } else {
            quality
        }
    }

    /// Apply the bias to a vector of quality scores.
    ///
    /// # Arguments
    ///
    /// * `qualities` - Original quality scores
    /// * `is_r2` - Whether this is an R2 read
    ///
    /// # Returns
    ///
    /// Adjusted quality scores.
    #[must_use]
    pub fn apply_to_vec(&self, qualities: &[u8], is_r2: bool) -> Vec<u8> {
        qualities.iter().map(|&q| self.apply(q, is_r2)).collect()
    }
}

#[cfg(test)]
#[allow(clippy::cast_precision_loss)]
mod tests {
    use super::*;
    use crate::simulate::create_rng;

    #[test]
    fn test_quality_model_default() {
        let model = PositionQualityModel::default();
        assert_eq!(model.warmup_bases, 10);
        assert_eq!(model.peak_quality, 37);
    }

    #[test]
    fn test_quality_warmup_phase() {
        let model = PositionQualityModel {
            noise_stddev: 0.0, // No noise for deterministic test
            ..Default::default()
        };
        let mut rng = create_rng(Some(42));

        // Position 0 should be at warmup_start
        let q0 = model.quality_at(0, &mut rng);
        assert_eq!(q0, model.warmup_start);

        // Position at end of warmup should be at peak
        let q_end_warmup = model.quality_at(model.warmup_bases, &mut rng);
        assert_eq!(q_end_warmup, model.peak_quality);
    }

    #[test]
    fn test_quality_peak_phase() {
        let model = PositionQualityModel { noise_stddev: 0.0, ..Default::default() };
        let mut rng = create_rng(Some(42));

        // Mid-read should be at peak
        let q_mid = model.quality_at(50, &mut rng);
        assert_eq!(q_mid, model.peak_quality);
    }

    #[test]
    fn test_quality_decay_phase() {
        let model = PositionQualityModel {
            noise_stddev: 0.0,
            decay_start: 100,
            decay_rate: 0.1,
            ..Default::default()
        };
        let mut rng = create_rng(Some(42));

        // At decay_start, should still be peak
        let q_start = model.quality_at(100, &mut rng);
        assert_eq!(q_start, model.peak_quality);

        // 10 positions into decay, should drop by 1
        let q_decay = model.quality_at(110, &mut rng);
        assert_eq!(q_decay, model.peak_quality - 1);
    }

    #[test]
    fn test_quality_clamped_to_min() {
        let model = PositionQualityModel {
            noise_stddev: 0.0,
            decay_start: 50,
            decay_rate: 1.0, // Aggressive decay
            min_quality: 10,
            ..Default::default()
        };
        let mut rng = create_rng(Some(42));

        // Far into decay should be clamped to min
        let q_late = model.quality_at(200, &mut rng);
        assert_eq!(q_late, model.min_quality);
    }

    #[test]
    fn test_generate_qualities_length() {
        let model = PositionQualityModel::default();
        let mut rng = create_rng(Some(42));

        let quals = model.generate_qualities(150, &mut rng);
        assert_eq!(quals.len(), 150);
    }

    #[test]
    fn test_r2_bias_default() {
        let bias = ReadPairQualityBias::default();
        assert_eq!(bias.r2_offset, -2);
    }

    #[test]
    fn test_r2_bias_apply() {
        let bias = ReadPairQualityBias::new(-3);

        // R1 unchanged
        assert_eq!(bias.apply(30, false), 30);

        // R2 reduced
        assert_eq!(bias.apply(30, true), 27);
    }

    #[test]
    fn test_r2_bias_clamped() {
        let bias = ReadPairQualityBias::new(-10);

        // Should clamp to minimum of 2
        assert_eq!(bias.apply(5, true), 2);
    }

    #[test]
    fn test_r2_bias_apply_to_vec() {
        let bias = ReadPairQualityBias::new(-2);
        let quals = vec![30, 35, 40];

        let r1_quals = bias.apply_to_vec(&quals, false);
        assert_eq!(r1_quals, vec![30, 35, 40]);

        let r2_quals = bias.apply_to_vec(&quals, true);
        assert_eq!(r2_quals, vec![28, 33, 38]);
    }

    #[test]
    fn test_quality_model_custom_parameters() {
        let model = PositionQualityModel::new(5, 20, 40, 50, 0.2, 5, 0.0);
        let mut rng = create_rng(Some(42));

        // Position 0 should be at warmup_start
        assert_eq!(model.quality_at(0, &mut rng), 20);

        // Position at end of warmup should be at peak
        assert_eq!(model.quality_at(5, &mut rng), 40);

        // Position after decay should be lower
        let q_decay = model.quality_at(55, &mut rng);
        // 55 - 50 = 5 positions into decay, 5 * 0.2 = 1.0 drop
        assert_eq!(q_decay, 39);
    }

    #[test]
    fn test_quality_model_zero_warmup() {
        let model = PositionQualityModel::new(0, 25, 37, 100, 0.08, 2, 0.0);
        let mut rng = create_rng(Some(42));

        // With zero warmup, position 0 should be at peak
        assert_eq!(model.quality_at(0, &mut rng), 37);
    }

    #[test]
    fn test_quality_model_immediate_decay() {
        let model = PositionQualityModel::new(0, 25, 37, 0, 0.5, 2, 0.0);
        let mut rng = create_rng(Some(42));

        // Position 0 at peak (decay_start = 0 means decay starts at position 0)
        assert_eq!(model.quality_at(0, &mut rng), 37);

        // Position 10 should be decayed
        let q10 = model.quality_at(10, &mut rng);
        // 10 * 0.5 = 5 drop
        assert_eq!(q10, 32);
    }

    #[test]
    fn test_quality_model_high_noise() {
        let model = PositionQualityModel::new(10, 25, 37, 100, 0.08, 2, 5.0);
        let mut rng = create_rng(Some(42));

        // Generate many samples at peak position and verify variance
        let samples: Vec<u8> = (0..1000).map(|_| model.quality_at(50, &mut rng)).collect();

        // With high noise, we should see variation
        let min = *samples.iter().min().unwrap();
        let max = *samples.iter().max().unwrap();
        assert!(max - min > 5, "Expected variation with high noise");
    }

    #[test]
    fn test_quality_model_quality_clamped_to_max() {
        let model = PositionQualityModel::new(10, 25, 50, 100, 0.08, 2, 0.0);
        let mut rng = create_rng(Some(42));

        // Even with peak_quality > 41, should clamp to 41
        let q = model.quality_at(50, &mut rng);
        assert!(q <= 41);
    }

    #[test]
    fn test_generate_qualities_shape() {
        let model = PositionQualityModel {
            noise_stddev: 0.0,
            warmup_bases: 10,
            decay_start: 50,
            ..Default::default()
        };
        let mut rng = create_rng(Some(42));

        let quals = model.generate_qualities(100, &mut rng);

        // Check shape: warmup -> peak -> decay
        assert!(quals[0] < quals[15]); // warmup is lower than peak
        assert!(quals[90] < quals[30]); // decay is lower than peak
    }

    #[test]
    fn test_r2_bias_zero_offset() {
        let bias = ReadPairQualityBias::new(0);

        assert_eq!(bias.apply(30, false), 30);
        assert_eq!(bias.apply(30, true), 30); // No change with zero offset
    }

    #[test]
    fn test_r2_bias_upper_clamp() {
        let bias = ReadPairQualityBias::new(10);

        // With positive offset, should clamp to max 41
        assert_eq!(bias.apply(40, true), 41);
        assert_eq!(bias.apply(35, true), 41);
    }

    #[test]
    fn test_r2_bias_empty_vec() {
        let bias = ReadPairQualityBias::default();
        let empty: Vec<u8> = vec![];

        assert_eq!(bias.apply_to_vec(&empty, false), Vec::<u8>::new());
        assert_eq!(bias.apply_to_vec(&empty, true), Vec::<u8>::new());
    }

    #[test]
    fn test_quality_model_very_long_read() {
        let model = PositionQualityModel {
            noise_stddev: 0.0,
            decay_start: 100,
            decay_rate: 0.1,
            min_quality: 5,
            ..Default::default()
        };
        let mut rng = create_rng(Some(42));

        // Very long read should decay to min_quality
        let quals = model.generate_qualities(500, &mut rng);
        assert_eq!(quals[450], 5); // Should be at floor
    }

    #[test]
    fn test_quality_at_reproducibility() {
        let model = PositionQualityModel::default();

        let mut rng1 = create_rng(Some(42));
        let mut rng2 = create_rng(Some(42));

        let q1 = model.quality_at(50, &mut rng1);
        let q2 = model.quality_at(50, &mut rng2);

        assert_eq!(q1, q2);
    }
}

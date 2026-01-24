//! Shared CLI arguments for simulation commands.

use clap::Args;
use std::path::PathBuf;

/// Common simulation options shared across all simulate subcommands.
#[derive(Args, Debug, Clone)]
pub struct SimulationCommon {
    /// Random seed for reproducibility
    #[arg(long = "seed")]
    pub seed: Option<u64>,

    /// Number of molecules to simulate
    #[arg(short = 'n', long = "num-molecules", default_value = "1000")]
    pub num_molecules: usize,

    /// Read length in bases
    #[arg(short = 'l', long = "read-length", default_value = "150")]
    pub read_length: usize,

    /// UMI length in bases
    #[arg(short = 'u', long = "umi-length", default_value = "8")]
    pub umi_length: usize,
}

/// Quality model options.
#[derive(Args, Debug, Clone)]
pub struct QualityArgs {
    /// Number of bases before peak quality is reached
    #[arg(long = "warmup-bases", default_value = "10")]
    pub warmup_bases: usize,

    /// Starting quality score during warmup phase
    #[arg(long = "warmup-quality", default_value = "25")]
    pub warmup_quality: u8,

    /// Peak quality score (Phred)
    #[arg(long = "peak-quality", default_value = "37")]
    pub peak_quality: u8,

    /// Position where quality decay begins
    #[arg(long = "decay-start", default_value = "100")]
    pub decay_start: usize,

    /// Quality drop per base after decay starts
    #[arg(long = "decay-rate", default_value = "0.08")]
    pub decay_rate: f64,

    /// Standard deviation of quality noise
    #[arg(long = "quality-noise", default_value = "2.0")]
    pub quality_noise: f64,

    /// Quality offset for R2 reads (typically negative)
    #[arg(long = "r2-quality-offset", default_value = "-2", allow_hyphen_values = true)]
    pub r2_quality_offset: i8,
}

impl QualityArgs {
    /// Convert to a [`PositionQualityModel`].
    pub fn to_quality_model(&self) -> fgumi_lib::simulate::PositionQualityModel {
        fgumi_lib::simulate::PositionQualityModel::new(
            self.warmup_bases,
            self.warmup_quality,
            self.peak_quality,
            self.decay_start,
            self.decay_rate,
            2, // min_quality
            self.quality_noise,
        )
    }

    /// Convert to a [`ReadPairQualityBias`].
    pub fn to_quality_bias(&self) -> fgumi_lib::simulate::ReadPairQualityBias {
        fgumi_lib::simulate::ReadPairQualityBias::new(self.r2_quality_offset)
    }
}

/// Insert size distribution options.
#[derive(Args, Debug, Clone)]
pub struct InsertSizeArgs {
    /// Mean insert size
    #[arg(long = "insert-size-mean", default_value = "300.0")]
    pub insert_size_mean: f64,

    /// Insert size standard deviation
    #[arg(long = "insert-size-stddev", default_value = "50.0")]
    pub insert_size_stddev: f64,

    /// Minimum insert size
    #[arg(long = "insert-size-min", default_value = "50")]
    pub insert_size_min: usize,

    /// Maximum insert size
    #[arg(long = "insert-size-max", default_value = "800")]
    pub insert_size_max: usize,
}

impl InsertSizeArgs {
    /// Convert to an [`InsertSizeModel`].
    pub fn to_insert_size_model(&self) -> fgumi_lib::simulate::InsertSizeModel {
        fgumi_lib::simulate::InsertSizeModel::new(
            self.insert_size_mean,
            self.insert_size_stddev,
            self.insert_size_min,
            self.insert_size_max,
        )
    }
}

/// Family size distribution options.
#[derive(Args, Debug, Clone)]
pub struct FamilySizeArgs {
    /// Family size distribution: "lognormal", "negbin", or path to histogram file
    #[arg(long = "family-size-dist", default_value = "lognormal")]
    pub family_size_dist: String,

    /// Mean family size (for lognormal)
    #[arg(long = "family-size-mean", default_value = "3.0")]
    pub family_size_mean: f64,

    /// Family size standard deviation (for lognormal)
    #[arg(long = "family-size-stddev", default_value = "2.0")]
    pub family_size_stddev: f64,

    /// r parameter for negative binomial
    #[arg(long = "family-size-r", default_value = "2.0")]
    pub family_size_r: f64,

    /// p parameter for negative binomial
    #[arg(long = "family-size-p", default_value = "0.5")]
    pub family_size_p: f64,

    /// Minimum reads per family
    #[arg(long = "min-family-size", default_value = "1")]
    pub min_family_size: usize,
}

impl FamilySizeArgs {
    /// Convert to a [`FamilySizeDistribution`].
    pub fn to_family_size_distribution(
        &self,
    ) -> anyhow::Result<fgumi_lib::simulate::FamilySizeDistribution> {
        match self.family_size_dist.as_str() {
            "lognormal" => Ok(fgumi_lib::simulate::FamilySizeDistribution::log_normal(
                self.family_size_mean,
                self.family_size_stddev,
            )),
            "negbin" => Ok(fgumi_lib::simulate::FamilySizeDistribution::negative_binomial(
                self.family_size_r,
                self.family_size_p,
            )),
            path => {
                // Treat as a path to a histogram file
                fgumi_lib::simulate::FamilySizeDistribution::from_histogram(path)
            }
        }
    }
}

/// Strand bias options for duplex simulation.
#[derive(Args, Debug, Clone)]
pub struct StrandBiasArgs {
    /// Beta distribution alpha for A/B strand ratio
    #[arg(long = "strand-alpha", default_value = "5.0")]
    pub strand_alpha: f64,

    /// Beta distribution beta for A/B strand ratio
    #[arg(long = "strand-beta", default_value = "5.0")]
    pub strand_beta: f64,
}

impl StrandBiasArgs {
    /// Convert to a [`StrandBiasModel`].
    pub fn to_strand_bias_model(&self) -> fgumi_lib::simulate::StrandBiasModel {
        fgumi_lib::simulate::StrandBiasModel::new(self.strand_alpha, self.strand_beta)
    }
}

/// Reference options for mapped reads.
#[derive(Args, Debug, Clone)]
pub struct ReferenceArgs {
    /// Reference FASTA file (sequences sampled from here)
    #[arg(short = 'r', long = "reference")]
    pub reference: Option<PathBuf>,

    /// Synthetic reference name (only used if no --reference)
    #[arg(long = "ref-name", default_value = "chr1")]
    pub ref_name: String,

    /// Synthetic reference length (only used if no --reference).
    /// A larger value prevents position collisions with many molecules.
    #[arg(long = "ref-length", default_value = "250000000")]
    pub ref_length: usize,
}

/// Position distribution options for mapped reads.
#[derive(Args, Debug, Clone)]
pub struct PositionDistArgs {
    /// Number of genomic positions to use (default: same as num-molecules)
    #[arg(long = "num-positions")]
    pub num_positions: Option<usize>,

    /// Number of unique UMIs per position
    #[arg(long = "umis-per-position", default_value = "1")]
    pub umis_per_position: usize,
}

#[cfg(test)]
mod tests {
    use super::*;
    use fgumi_lib::simulate::create_rng;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_quality_args_default_model() {
        let args = QualityArgs {
            warmup_bases: 10,
            warmup_quality: 25,
            peak_quality: 37,
            decay_start: 100,
            decay_rate: 0.08,
            quality_noise: 2.0,
            r2_quality_offset: -2,
        };

        let model = args.to_quality_model();
        assert_eq!(model.warmup_bases, 10);
        assert_eq!(model.warmup_start, 25);
        assert_eq!(model.peak_quality, 37);
        assert_eq!(model.decay_start, 100);
        assert!((model.decay_rate - 0.08).abs() < f64::EPSILON);
    }

    #[test]
    fn test_quality_args_custom_values() {
        let args = QualityArgs {
            warmup_bases: 5,
            warmup_quality: 20,
            peak_quality: 40,
            decay_start: 80,
            decay_rate: 0.1,
            quality_noise: 1.5,
            r2_quality_offset: -3,
        };

        let model = args.to_quality_model();
        assert_eq!(model.warmup_bases, 5);
        assert_eq!(model.peak_quality, 40);

        let bias = args.to_quality_bias();
        assert_eq!(bias.r2_offset, -3);
    }

    #[test]
    fn test_quality_bias_positive_offset() {
        let args = QualityArgs {
            warmup_bases: 10,
            warmup_quality: 25,
            peak_quality: 37,
            decay_start: 100,
            decay_rate: 0.08,
            quality_noise: 2.0,
            r2_quality_offset: 3, // Positive offset (unusual but valid)
        };

        let bias = args.to_quality_bias();
        assert_eq!(bias.apply(30, true), 33);
    }

    #[test]
    fn test_insert_size_args_default() {
        let args = InsertSizeArgs {
            insert_size_mean: 300.0,
            insert_size_stddev: 50.0,
            insert_size_min: 50,
            insert_size_max: 800,
        };

        let model = args.to_insert_size_model();
        assert!((model.mean - 300.0).abs() < f64::EPSILON);
        assert_eq!(model.min, 50);
        assert_eq!(model.max, 800);
    }

    #[test]
    fn test_insert_size_args_narrow_range() {
        let args = InsertSizeArgs {
            insert_size_mean: 200.0,
            insert_size_stddev: 10.0,
            insert_size_min: 180,
            insert_size_max: 220,
        };

        let model = args.to_insert_size_model();
        let mut rng = create_rng(Some(42));

        // All samples should be within narrow range
        for _ in 0..100 {
            let size = model.sample(&mut rng);
            assert!((180..=220).contains(&size));
        }
    }

    #[test]
    fn test_family_size_args_lognormal() {
        let args = FamilySizeArgs {
            family_size_dist: "lognormal".to_string(),
            family_size_mean: 3.0,
            family_size_stddev: 2.0,
            family_size_r: 2.0,
            family_size_p: 0.5,
            min_family_size: 1,
        };

        let dist = args.to_family_size_distribution().unwrap();
        let mut rng = create_rng(Some(42));
        let size = dist.sample(&mut rng, 1);
        assert!(size >= 1);
    }

    #[test]
    fn test_family_size_args_negbin() {
        let args = FamilySizeArgs {
            family_size_dist: "negbin".to_string(),
            family_size_mean: 3.0,
            family_size_stddev: 2.0,
            family_size_r: 2.0,
            family_size_p: 0.5,
            min_family_size: 1,
        };

        let dist = args.to_family_size_distribution().unwrap();
        let mut rng = create_rng(Some(42));
        let size = dist.sample(&mut rng, 1);
        assert!(size >= 1);
    }

    #[test]
    fn test_family_size_args_from_histogram() -> anyhow::Result<()> {
        let mut temp = NamedTempFile::new()?;
        writeln!(temp, "family_size\tcount")?;
        writeln!(temp, "1\t50")?;
        writeln!(temp, "2\t30")?;
        writeln!(temp, "3\t20")?;
        temp.flush()?;

        let args = FamilySizeArgs {
            family_size_dist: temp.path().to_string_lossy().to_string(),
            family_size_mean: 3.0,
            family_size_stddev: 2.0,
            family_size_r: 2.0,
            family_size_p: 0.5,
            min_family_size: 1,
        };

        let dist = args.to_family_size_distribution()?;
        let mut rng = create_rng(Some(42));

        // Sample and verify all sizes are from histogram
        for _ in 0..100 {
            let size = dist.sample(&mut rng, 1);
            assert!((1..=3).contains(&size));
        }

        Ok(())
    }

    #[test]
    fn test_family_size_args_invalid_histogram() {
        let args = FamilySizeArgs {
            family_size_dist: "/nonexistent/path/histogram.tsv".to_string(),
            family_size_mean: 3.0,
            family_size_stddev: 2.0,
            family_size_r: 2.0,
            family_size_p: 0.5,
            min_family_size: 1,
        };

        let result = args.to_family_size_distribution();
        assert!(result.is_err());
    }

    #[test]
    fn test_strand_bias_args_symmetric() {
        let args = StrandBiasArgs { strand_alpha: 5.0, strand_beta: 5.0 };

        let model = args.to_strand_bias_model();
        assert!((model.alpha - 5.0).abs() < f64::EPSILON);
        assert!((model.beta - 5.0).abs() < f64::EPSILON);

        let mut rng = create_rng(Some(42));
        let fractions: Vec<f64> = (0..1000).map(|_| model.sample_a_fraction(&mut rng)).collect();
        let mean: f64 = fractions.iter().sum::<f64>() / fractions.len() as f64;

        // Mean should be close to 0.5 for symmetric distribution
        assert!((mean - 0.5).abs() < 0.05);
    }

    #[test]
    fn test_strand_bias_args_a_biased() {
        let args = StrandBiasArgs { strand_alpha: 8.0, strand_beta: 2.0 };

        let model = args.to_strand_bias_model();
        let mut rng = create_rng(Some(42));

        let fractions: Vec<f64> = (0..1000).map(|_| model.sample_a_fraction(&mut rng)).collect();
        let mean: f64 = fractions.iter().sum::<f64>() / fractions.len() as f64;

        // Mean should be biased toward A (> 0.5)
        assert!(mean > 0.7);
    }

    #[test]
    fn test_strand_bias_args_b_biased() {
        let args = StrandBiasArgs { strand_alpha: 2.0, strand_beta: 8.0 };

        let model = args.to_strand_bias_model();
        let mut rng = create_rng(Some(42));

        let fractions: Vec<f64> = (0..1000).map(|_| model.sample_a_fraction(&mut rng)).collect();
        let mean: f64 = fractions.iter().sum::<f64>() / fractions.len() as f64;

        // Mean should be biased toward B (< 0.5)
        assert!(mean < 0.3);
    }

    #[test]
    fn test_quality_args_zero_warmup() {
        let args = QualityArgs {
            warmup_bases: 0,
            warmup_quality: 25,
            peak_quality: 37,
            decay_start: 100,
            decay_rate: 0.08,
            quality_noise: 2.0,
            r2_quality_offset: -2,
        };

        let model = args.to_quality_model();
        assert_eq!(model.warmup_bases, 0);
    }

    #[test]
    fn test_quality_args_high_peak_quality() {
        let args = QualityArgs {
            warmup_bases: 10,
            warmup_quality: 30,
            peak_quality: 41, // Max Phred quality
            decay_start: 100,
            decay_rate: 0.08,
            quality_noise: 0.0, // No noise for predictable testing
            r2_quality_offset: 0,
        };

        let model = args.to_quality_model();
        let mut rng = create_rng(Some(42));

        // Generate qualities and check peak region
        let quals = model.generate_qualities(50, &mut rng);
        assert!(!quals.is_empty());
    }

    #[test]
    fn test_insert_size_args_min_equals_max() {
        let args = InsertSizeArgs {
            insert_size_mean: 200.0,
            insert_size_stddev: 50.0,
            insert_size_min: 200,
            insert_size_max: 200,
        };

        let model = args.to_insert_size_model();
        let mut rng = create_rng(Some(42));

        // All samples should be exactly 200 when min == max
        for _ in 0..10 {
            let size = model.sample(&mut rng);
            assert_eq!(size, 200);
        }
    }

    #[test]
    fn test_family_size_args_high_min() {
        let args = FamilySizeArgs {
            family_size_dist: "lognormal".to_string(),
            family_size_mean: 3.0,
            family_size_stddev: 2.0,
            family_size_r: 2.0,
            family_size_p: 0.5,
            min_family_size: 5, // High minimum
        };

        let dist = args.to_family_size_distribution().unwrap();
        let mut rng = create_rng(Some(42));

        for _ in 0..100 {
            let size = dist.sample(&mut rng, 5);
            assert!(size >= 5, "Size {size} should be >= min 5");
        }
    }

    #[test]
    fn test_strand_bias_split_reads() {
        let args = StrandBiasArgs { strand_alpha: 5.0, strand_beta: 5.0 };
        let model = args.to_strand_bias_model();
        let mut rng = create_rng(Some(42));

        for total in [0, 1, 2, 5, 10, 100] {
            let (a, b) = model.split_reads(total, &mut rng);
            assert_eq!(a + b, total, "A ({a}) + B ({b}) should equal total ({total})");
        }
    }

    #[test]
    fn test_strand_bias_split_zero_total() {
        let args = StrandBiasArgs { strand_alpha: 5.0, strand_beta: 5.0 };
        let model = args.to_strand_bias_model();
        let mut rng = create_rng(Some(42));

        let (a, b) = model.split_reads(0, &mut rng);
        assert_eq!(a, 0);
        assert_eq!(b, 0);
    }

    #[test]
    fn test_strand_bias_split_one_read() {
        let args = StrandBiasArgs { strand_alpha: 5.0, strand_beta: 5.0 };
        let model = args.to_strand_bias_model();
        let mut rng = create_rng(Some(42));

        // With 1 read, should go to either A or B
        let (a, b) = model.split_reads(1, &mut rng);
        assert_eq!(a + b, 1);
        assert!(a <= 1 && b <= 1);
    }

    #[test]
    fn test_quality_args_zero_noise() {
        let args = QualityArgs {
            warmup_bases: 10,
            warmup_quality: 25,
            peak_quality: 37,
            decay_start: 100,
            decay_rate: 0.08,
            quality_noise: 0.0, // No noise
            r2_quality_offset: 0,
        };

        let model = args.to_quality_model();
        assert!((model.noise_stddev - 0.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_insert_size_sample_distribution() {
        let args = InsertSizeArgs {
            insert_size_mean: 300.0,
            insert_size_stddev: 50.0,
            insert_size_min: 50,
            insert_size_max: 800,
        };

        let model = args.to_insert_size_model();
        let mut rng = create_rng(Some(42));

        let samples: Vec<usize> = (0..1000).map(|_| model.sample(&mut rng)).collect();
        let mean: f64 = samples.iter().map(|&s| s as f64).sum::<f64>() / samples.len() as f64;

        // Mean should be close to 300
        assert!(mean > 280.0 && mean < 320.0, "Mean {mean} not close to expected 300");
    }

    #[test]
    fn test_family_size_distribution_type_case_insensitive() {
        // Test that "LOGNORMAL" and "lognormal" both work
        let args_lower = FamilySizeArgs {
            family_size_dist: "lognormal".to_string(),
            family_size_mean: 3.0,
            family_size_stddev: 2.0,
            family_size_r: 2.0,
            family_size_p: 0.5,
            min_family_size: 1,
        };

        // Should not panic
        let _ = args_lower.to_family_size_distribution().unwrap();
    }
}

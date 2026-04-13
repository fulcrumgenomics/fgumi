//! Shared CLI arguments and utilities for simulation commands.

use super::sort::TemplateCoordKey;
use crate::commands::common::MethylationModeArg;
use anyhow::{Context, Result, bail};
use clap::Args;
use fgumi_consensus::MethylationMode;
use fgumi_consensus::methylation::is_cpg_context;
use log::info;
use noodles::fasta;
use rand::{Rng, RngExt};
use std::fs::File;
use std::io::BufReader;
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

    /// Standard deviation of quality noise (must be finite and >= 0)
    #[arg(long = "quality-noise", default_value = "2.0", value_parser = parse_noise_stddev)]
    pub quality_noise: f64,

    /// Quality offset for R2 reads (typically negative)
    #[arg(long = "r2-quality-offset", default_value = "-2", allow_hyphen_values = true)]
    pub r2_quality_offset: i8,
}

/// Parse and validate a noise standard deviation value for quality score simulation.
fn parse_noise_stddev(s: &str) -> Result<f64, String> {
    let val: f64 = s.parse().map_err(|e| format!("invalid float: {e}"))?;
    if !val.is_finite() || val < 0.0 {
        return Err(format!("quality-noise must be finite and >= 0.0, got {val}"));
    }
    Ok(val)
}

impl QualityArgs {
    /// Convert to a [`PositionQualityModel`].
    pub fn to_quality_model(&self) -> crate::simulate::PositionQualityModel {
        crate::simulate::PositionQualityModel::new(
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
    pub fn to_quality_bias(&self) -> crate::simulate::ReadPairQualityBias {
        crate::simulate::ReadPairQualityBias::new(self.r2_quality_offset)
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
    pub fn to_insert_size_model(&self) -> crate::simulate::InsertSizeModel {
        crate::simulate::InsertSizeModel::new(
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
    ) -> anyhow::Result<crate::simulate::FamilySizeDistribution> {
        match self.family_size_dist.as_str() {
            "lognormal" => Ok(crate::simulate::FamilySizeDistribution::log_normal(
                self.family_size_mean,
                self.family_size_stddev,
            )),
            "negbin" => Ok(crate::simulate::FamilySizeDistribution::negative_binomial(
                self.family_size_r,
                self.family_size_p,
            )),
            path => {
                // Treat as a path to a histogram file
                crate::simulate::FamilySizeDistribution::from_histogram(path)
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
    pub fn to_strand_bias_model(&self) -> crate::simulate::StrandBiasModel {
        crate::simulate::StrandBiasModel::new(self.strand_alpha, self.strand_beta)
    }
}

/// Methylation simulation options shared across simulate subcommands.
#[derive(Args, Debug, Clone)]
pub struct MethylationArgs {
    /// Methylation chemistry mode. When set, enables methylation-aware base
    /// conversion in simulated reads. Requires --reference for commands that
    /// generate read bases (fastq-reads, mapped-reads, grouped-reads).
    #[arg(long = "methylation-mode", value_enum)]
    pub methylation_mode: Option<MethylationModeArg>,

    /// Fraction of `CpG` cytosines that are methylated [0.0-1.0].
    /// Methylated `CpG`s are protected from conversion in EM-Seq and are
    /// targets for conversion in TAPs.
    #[arg(long = "cpg-methylation-rate", default_value = "0.75")]
    pub cpg_methylation_rate: f64,

    /// Enzymatic conversion efficiency for target cytosines [0.0-1.0].
    /// In EM-Seq, this is the probability that an unmethylated C is converted to T.
    /// In TAPs, this is the probability that a methylated C is converted to T.
    #[arg(long = "conversion-rate", default_value = "0.98")]
    pub conversion_rate: f64,
}

impl MethylationArgs {
    /// Resolves the optional CLI arg to a [`MethylationConfig`].
    pub fn resolve(&self) -> MethylationConfig {
        MethylationConfig {
            mode: crate::commands::common::resolve_methylation_mode(self.methylation_mode),
            cpg_methylation_rate: self.cpg_methylation_rate,
            conversion_rate: self.conversion_rate,
        }
    }

    /// Validates that rate parameters are in [0.0, 1.0] and finite.
    pub fn validate(&self) -> anyhow::Result<()> {
        validate_rate(self.cpg_methylation_rate, "cpg-methylation-rate")?;
        validate_rate(self.conversion_rate, "conversion-rate")?;
        Ok(())
    }
}

/// Resolved methylation simulation parameters.
#[derive(Debug, Clone, Copy)]
pub struct MethylationConfig {
    /// Methylation chemistry mode.
    pub mode: MethylationMode,
    /// Fraction of `CpG` cytosines that are methylated.
    pub cpg_methylation_rate: f64,
    /// Enzymatic conversion efficiency.
    pub conversion_rate: f64,
}

/// Validates that a rate is a finite value in [0.0, 1.0].
fn validate_rate(value: f64, name: &str) -> anyhow::Result<()> {
    if !value.is_finite() || !(0.0..=1.0).contains(&value) {
        anyhow::bail!("--{name} must be a finite value between 0.0 and 1.0, got {value}");
    }
    Ok(())
}

/// Reference options for simulate commands.
#[derive(Args, Debug, Clone)]
pub struct ReferenceArgs {
    /// Reference FASTA file for sampling template sequences and building BAM headers.
    #[arg(short = 'r', long = "reference", required = true)]
    pub reference: PathBuf,
}

/// Minimum contig length (bp) considered usable for sampling. Contigs shorter than this
/// are skipped when loading the reference so that the 1000-bp N-check window and typical
/// insert sizes have room to fit.
const MIN_CONTIG_LENGTH: usize = 1000;

/// Loaded reference genome for sampling template sequences.
pub(super) struct ReferenceGenome {
    names: Vec<String>,
    sequences: Vec<Vec<u8>>,
    cumulative_lengths: Vec<usize>,
    total_length: usize,
}

impl ReferenceGenome {
    /// Load a reference genome from a FASTA file.
    pub fn load<P: AsRef<std::path::Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        info!("Loading reference from {}", path.display());

        let file = File::open(path)
            .with_context(|| format!("Failed to open reference: {}", path.display()))?;
        let reader = BufReader::new(file);
        let mut fasta_reader = fasta::io::Reader::new(reader);

        let mut names = Vec::new();
        let mut sequences = Vec::new();
        let mut cumulative_lengths = Vec::new();
        let mut total_length = 0usize;

        for result in fasta_reader.records() {
            let record = result.with_context(|| "Failed to read FASTA record")?;
            let name = std::str::from_utf8(record.name())
                .with_context(|| "Invalid chromosome name")?
                .to_string();
            let seq: Vec<u8> =
                record.sequence().as_ref().iter().map(|&b| b.to_ascii_uppercase()).collect();

            if seq.len() >= MIN_CONTIG_LENGTH {
                // Only include chromosomes with sufficient length
                total_length += seq.len();
                cumulative_lengths.push(total_length);
                names.push(name);
                sequences.push(seq);
            }
        }

        if sequences.is_empty() {
            bail!("No valid sequences found in reference FASTA");
        }

        info!("Loaded {} chromosomes, total {} bp", sequences.len(), total_length);

        Ok(Self { names, sequences, cumulative_lengths, total_length })
    }

    /// Returns the chromosome name at the given index.
    pub fn name(&self, chrom_idx: usize) -> &str {
        &self.names[chrom_idx]
    }

    /// Sample a random position and return (`chrom_idx`, position, sequence).
    /// Returns None if the position contains N bases or if `length` is zero or
    /// exceeds `total_length`.
    pub fn sample_sequence(
        &self,
        length: usize,
        rng: &mut impl Rng,
    ) -> Option<(usize, usize, Vec<u8>)> {
        if length == 0 || length > self.total_length {
            return None;
        }
        // Exclusive upper bound keeps `genome_pos < total_length`, which guarantees
        // `partition_point` returns a valid index into `cumulative_lengths`.
        let start_bound = self.total_length - length + 1;
        // Try up to 10 times to find a valid position without N bases
        for _ in 0..10 {
            // Pick a random position in the genome
            let genome_pos = rng.random_range(0..start_bound);

            // Find which chromosome this falls in (binary search on sorted cumulative lengths)
            let chrom_idx = self.cumulative_lengths.partition_point(|&cum| cum <= genome_pos);

            let chrom_start =
                if chrom_idx == 0 { 0 } else { self.cumulative_lengths[chrom_idx - 1] };
            let local_pos = genome_pos - chrom_start;

            let seq = &self.sequences[chrom_idx];
            if local_pos + length > seq.len() {
                continue;
            }

            let template = &seq[local_pos..local_pos + length];

            // Check for N bases
            if template.iter().any(|&b| b == b'N' || b == b'n') {
                continue;
            }

            return Some((chrom_idx, local_pos, template.to_vec()));
        }
        None
    }

    /// Returns the total genome length (sum of all loaded chromosome sequences).
    pub fn total_length(&self) -> usize {
        self.total_length
    }

    /// Return the subsequence at a genome-wide position, mapping to the correct
    /// chromosome automatically. The position wraps around the total genome length.
    /// Returns None if out of bounds or if the sequence contains N bases.
    pub fn sequence_at_genome_pos(&self, genome_pos: usize, length: usize) -> Option<Vec<u8>> {
        if self.total_length == 0 {
            return None;
        }
        let genome_pos = genome_pos % self.total_length;
        let chrom_idx = self.cumulative_lengths.partition_point(|&cum| cum <= genome_pos);
        let chrom_start = if chrom_idx == 0 { 0 } else { self.cumulative_lengths[chrom_idx - 1] };
        let local_pos = genome_pos - chrom_start;
        self.sequence_at(chrom_idx, local_pos, length)
    }

    /// Return the subsequence at a specific chromosome and position.
    /// Returns None if out of bounds or if the sequence contains N bases.
    pub fn sequence_at(&self, chrom_idx: usize, pos: usize, length: usize) -> Option<Vec<u8>> {
        if chrom_idx >= self.sequences.len() {
            return None;
        }
        let seq = &self.sequences[chrom_idx];
        if pos + length > seq.len() {
            return None;
        }
        let subseq = &seq[pos..pos + length];
        if subseq.iter().any(|&b| b == b'N' || b == b'n') {
            return None;
        }
        Some(subseq.to_vec())
    }

    /// Build a BAM [`Header`] with `@SQ` lines for every loaded contig.
    pub(super) fn build_bam_header(&self) -> noodles::sam::header::Header {
        use bstr::BString;
        use noodles::sam::header::Header;
        use noodles::sam::header::record::value::Map;
        use noodles::sam::header::record::value::map::ReferenceSequence;
        use std::num::NonZeroUsize;

        let mut builder = Header::builder();
        for (name, seq) in self.names.iter().zip(self.sequences.iter()) {
            let length = NonZeroUsize::try_from(seq.len()).expect("chromosome length must be > 0");
            let ref_seq = Map::<ReferenceSequence>::new(length);
            builder = builder.add_reference_sequence(BString::from(name.as_str()), ref_seq);
        }
        builder.build()
    }

    /// Returns the length of the longest loaded chromosome.
    pub(super) fn max_contig_length(&self) -> usize {
        self.sequences.iter().map(|s| s.len()).max().unwrap_or(0)
    }

    /// Returns the number of loaded chromosomes.
    #[allow(dead_code)] // used by grouped-reads and consensus-reads (upcoming tasks)
    pub(super) fn num_chromosomes(&self) -> usize {
        self.sequences.len()
    }

    /// Returns the length of the chromosome at the given index.
    #[allow(dead_code)] // used by grouped-reads and consensus-reads (upcoming tasks)
    pub(super) fn chromosome_length(&self, chrom_idx: usize) -> usize {
        self.sequences[chrom_idx].len()
    }

    /// Pre-sample `num_positions` random loci as `(chrom_idx, local_pos)` tuples.
    ///
    /// Positions are drawn uniformly across the genome and are checked against
    /// a `MIN_CONTIG_LENGTH` bp window for N bases. Sampling retries internally up to
    /// `num_positions * 100` attempts before panicking.
    pub(super) fn sample_positions(
        &self,
        num_positions: usize,
        rng: &mut impl Rng,
    ) -> Vec<(usize, usize)> {
        const WINDOW: usize = MIN_CONTIG_LENGTH;
        let max_attempts = num_positions.saturating_mul(100).max(1);
        let mut positions = Vec::with_capacity(num_positions);
        let mut attempts = 0usize;

        while positions.len() < num_positions {
            assert!(
                attempts < max_attempts,
                "sample_positions: exhausted {max_attempts} attempts to find \
                 {num_positions} N-free positions in the reference"
            );
            attempts += 1;

            // Pick a genome-wide position and map to (chrom_idx, local_pos)
            let genome_pos = rng.random_range(0..self.total_length);
            let chrom_idx = self.cumulative_lengths.partition_point(|&cum| cum <= genome_pos);
            let chrom_start =
                if chrom_idx == 0 { 0 } else { self.cumulative_lengths[chrom_idx - 1] };
            let local_pos = genome_pos - chrom_start;

            // Check a 1000bp window for N bases
            let seq = &self.sequences[chrom_idx];
            let window_end = (local_pos + WINDOW).min(seq.len());
            let window_start = local_pos.min(window_end);
            let window = &seq[window_start..window_end];
            if window.iter().any(|&b| b == b'N' || b == b'n') {
                continue;
            }

            positions.push((chrom_idx, local_pos));
        }
        positions
    }
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

/// Generate a random DNA sequence of the given length.
pub(super) fn generate_random_sequence(len: usize, rng: &mut impl Rng) -> Vec<u8> {
    const BASES: &[u8] = b"ACGT";
    let mut seq = Vec::with_capacity(len);
    for _ in 0..len {
        seq.push(BASES[rng.random_range(0..4)]);
    }
    seq
}

/// Pad a sequence to the target length with random bases, or truncate if too long.
pub(super) fn pad_sequence(mut seq: Vec<u8>, target_len: usize, rng: &mut impl Rng) -> Vec<u8> {
    while seq.len() < target_len {
        seq.push(b"ACGT"[rng.random_range(0..4)]);
    }
    seq.truncate(target_len);
    seq
}

/// Apply methylation conversion to a read sequence in-place.
///
/// For each position in `read_seq` that aligns to a reference C (top strand) or G
/// (bottom strand), determines the `CpG` context and applies stochastic conversion
/// based on the methylation mode and rates.
pub(super) fn apply_methylation_conversion(
    read_seq: &mut [u8],
    ref_seq: &[u8],
    ref_offset: usize,
    is_top_strand: bool,
    config: &MethylationConfig,
    rng: &mut impl Rng,
) {
    if !config.mode.is_enabled() {
        return;
    }

    for (i, base) in read_seq.iter_mut().enumerate() {
        let ref_pos = ref_offset + i;
        if ref_pos >= ref_seq.len() {
            break;
        }

        let ref_base = ref_seq[ref_pos].to_ascii_uppercase();

        if is_top_strand && ref_base == b'C' {
            let cpg = is_cpg_context(ref_seq, ref_pos, true);
            if is_conversion_target(config.mode, cpg, config.cpg_methylation_rate, rng)
                && rng.random::<f64>() < config.conversion_rate
            {
                *base = b'T';
            }
        } else if !is_top_strand && ref_base == b'G' {
            let cpg = is_cpg_context(ref_seq, ref_pos, false);
            if is_conversion_target(config.mode, cpg, config.cpg_methylation_rate, rng)
                && rng.random::<f64>() < config.conversion_rate
            {
                *base = b'A';
            }
        }
    }
}

/// Determines whether a cytosine at this position is a conversion target.
///
/// Returns true if the base should be converted (subject to `conversion_rate`).
///
/// - EM-Seq converts **unmethylated** C: non-`CpG` always, `CpG` when not methylated
/// - TAPs converts **methylated** C: `CpG` when methylated, non-`CpG` never
fn is_conversion_target(
    mode: MethylationMode,
    is_cpg: bool,
    cpg_methylation_rate: f64,
    rng: &mut impl Rng,
) -> bool {
    if is_cpg {
        let methylated = rng.random::<f64>() < cpg_methylation_rate;
        match mode {
            MethylationMode::EmSeq => !methylated, // convert unmethylated
            MethylationMode::Taps => methylated,   // convert methylated
            MethylationMode::Disabled => false,
        }
    } else {
        // Non-CpG cytosines are unmethylated
        match mode {
            MethylationMode::EmSeq => true, // unmethylated = target
            MethylationMode::Taps => false, // unmethylated = not a target
            MethylationMode::Disabled => false,
        }
    }
}

/// Lightweight molecule info for position-first sorting.
#[derive(Debug)]
pub(super) struct MoleculeInfo {
    pub mol_id: usize,
    pub seed: u64,
    pub sort_key: TemplateCoordKey,
}

impl Ord for MoleculeInfo {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.sort_key.cmp(&other.sort_key)
    }
}

impl PartialOrd for MoleculeInfo {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for MoleculeInfo {
    fn eq(&self, other: &Self) -> bool {
        self.sort_key == other.sort_key
    }
}

impl Eq for MoleculeInfo {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simulate::create_rng;
    use rstest::rstest;
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

        let dist =
            args.to_family_size_distribution().expect("lognormal distribution should be created");
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

        let dist =
            args.to_family_size_distribution().expect("negbin distribution should be created");
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

        let dist = args
            .to_family_size_distribution()
            .expect("lognormal distribution with high min should be created");
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

    #[rstest]
    #[case("0.0", 0.0)]
    #[case("2.5", 2.5)]
    #[case("100.0", 100.0)]
    fn test_parse_noise_stddev_accepts_valid(#[case] input: &str, #[case] expected: f64) {
        let parsed = parse_noise_stddev(input).expect("valid noise stddev should parse");
        assert!((parsed - expected).abs() < f64::EPSILON);
    }

    #[rstest]
    #[case("-0.1")]
    #[case("NaN")]
    #[case("inf")]
    #[case("-inf")]
    #[case("abc")]
    fn test_parse_noise_stddev_rejects_invalid(#[case] input: &str) {
        assert!(parse_noise_stddev(input).is_err(), "input should be rejected: {input}");
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
        let _ = args_lower
            .to_family_size_distribution()
            .expect("lowercase distribution name should be accepted");
    }

    fn make_molecule_info(mol_id: usize, tid1: i32, pos1: i64) -> MoleculeInfo {
        MoleculeInfo {
            mol_id,
            seed: 0,
            sort_key: TemplateCoordKey {
                tid1,
                tid2: 0,
                pos1,
                pos2: 0,
                neg1: false,
                neg2: false,
                mid: String::new(),
                name: String::new(),
                is_upper_of_pair: false,
            },
        }
    }

    #[test]
    fn test_molecule_info_ordering() {
        let a = make_molecule_info(0, 1, 100);
        let b = make_molecule_info(1, 1, 200);
        let c = make_molecule_info(2, 2, 50);

        // Same tid, earlier position comes first
        assert!(a < b);
        assert!(b > a);
        // Higher tid comes later regardless of position
        assert!(b < c);

        // PartialOrd consistent with Ord
        assert_eq!(a.partial_cmp(&b), Some(std::cmp::Ordering::Less));

        // Equal sort keys are equal
        let a2 = make_molecule_info(99, 1, 100);
        assert_eq!(a, a2);
        assert_eq!(a.cmp(&a2), std::cmp::Ordering::Equal);
    }

    // ========================================================================
    // ReferenceGenome tests
    // ========================================================================

    /// Create a temp FASTA file with a single chromosome of the given sequence.
    fn write_test_fasta(seq: &[u8]) -> NamedTempFile {
        let mut f = NamedTempFile::new().unwrap();
        writeln!(f, ">chr1").unwrap();
        f.write_all(seq).unwrap();
        writeln!(f).unwrap();
        f.flush().unwrap();
        f
    }

    #[test]
    fn test_reference_genome_load_and_sample() {
        // 1500 bases to exceed the 1000-bp minimum
        let seq = b"ACGT".repeat(375);
        let fasta = write_test_fasta(&seq);
        let genome = ReferenceGenome::load(fasta.path()).unwrap();
        assert_eq!(genome.name(0), "chr1");
        assert!(genome.sequence_at(0, 0, 1500).is_some());
    }

    #[test]
    fn test_reference_genome_sequence_at_valid() {
        let seq = b"ACGT".repeat(375);
        let fasta = write_test_fasta(&seq);
        let genome = ReferenceGenome::load(fasta.path()).unwrap();
        let subseq = genome.sequence_at(0, 4, 8).unwrap();
        assert_eq!(subseq, b"ACGTACGT");
    }

    #[test]
    fn test_reference_genome_sequence_at_out_of_bounds() {
        let seq = b"ACGT".repeat(375);
        let fasta = write_test_fasta(&seq);
        let genome = ReferenceGenome::load(fasta.path()).unwrap();
        assert!(genome.sequence_at(0, 1495, 10).is_none());
    }

    #[test]
    fn test_reference_genome_sequence_at_invalid_chrom() {
        let seq = b"ACGT".repeat(375);
        let fasta = write_test_fasta(&seq);
        let genome = ReferenceGenome::load(fasta.path()).unwrap();
        assert!(genome.sequence_at(99, 0, 10).is_none());
    }

    #[test]
    fn test_reference_genome_sequence_at_n_bases() {
        let mut seq = b"ACGT".repeat(375);
        seq[10] = b'N';
        let fasta = write_test_fasta(&seq);
        let genome = ReferenceGenome::load(fasta.path()).unwrap();
        // Region containing the N should return None
        assert!(genome.sequence_at(0, 8, 4).is_none());
        // Region before the N should succeed
        assert!(genome.sequence_at(0, 0, 4).is_some());
    }

    #[test]
    fn test_reference_genome_skips_short_sequences() {
        let mut f = NamedTempFile::new().unwrap();
        writeln!(f, ">short").unwrap();
        // 500 bases - below 1000 minimum
        let short = b"ACGT".repeat(125);
        f.write_all(&short).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">long").unwrap();
        let long = b"ACGT".repeat(375);
        f.write_all(&long).unwrap();
        writeln!(f).unwrap();
        f.flush().unwrap();

        let genome = ReferenceGenome::load(f.path()).unwrap();
        assert_eq!(genome.name(0), "long");
        assert!(genome.sequence_at(1, 0, 1).is_none()); // only one chromosome loaded
    }

    #[test]
    fn test_reference_genome_total_length() {
        let seq = b"ACGT".repeat(375);
        let fasta = write_test_fasta(&seq);
        let genome = ReferenceGenome::load(fasta.path()).unwrap();
        assert_eq!(genome.total_length(), 1500);
    }

    #[test]
    fn test_reference_genome_sequence_at_genome_pos_single_chrom() {
        let seq = b"ACGT".repeat(375);
        let fasta = write_test_fasta(&seq);
        let genome = ReferenceGenome::load(fasta.path()).unwrap();
        // Position 4 should map to chr1:4
        let subseq = genome.sequence_at_genome_pos(4, 8).unwrap();
        assert_eq!(subseq, b"ACGTACGT");
    }

    #[test]
    fn test_reference_genome_sequence_at_genome_pos_wraps_around() {
        let seq = b"ACGT".repeat(375);
        let fasta = write_test_fasta(&seq);
        let genome = ReferenceGenome::load(fasta.path()).unwrap();
        // Position beyond total_length should wrap
        let direct = genome.sequence_at_genome_pos(4, 8).unwrap();
        let wrapped = genome.sequence_at_genome_pos(4 + genome.total_length(), 8).unwrap();
        assert_eq!(direct, wrapped);
    }

    #[test]
    fn test_reference_genome_sequence_at_genome_pos_multi_chrom() {
        let mut f = NamedTempFile::new().unwrap();
        writeln!(f, ">chr1").unwrap();
        let chr1 = b"AAAA".repeat(375); // 1500 bp of A's
        f.write_all(&chr1).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">chr2").unwrap();
        let chr2 = b"CCCC".repeat(375); // 1500 bp of C's
        f.write_all(&chr2).unwrap();
        writeln!(f).unwrap();
        f.flush().unwrap();

        let genome = ReferenceGenome::load(f.path()).unwrap();
        assert_eq!(genome.total_length(), 3000);

        // Position 0 -> chr1, should be A's
        let from_chr1 = genome.sequence_at_genome_pos(0, 4).unwrap();
        assert_eq!(from_chr1, b"AAAA");

        // Position 1500 -> chr2, should be C's
        let from_chr2 = genome.sequence_at_genome_pos(1500, 4).unwrap();
        assert_eq!(from_chr2, b"CCCC");
    }

    #[test]
    fn test_reference_genome_sequence_at_genome_pos_boundary() {
        let seq = b"ACGT".repeat(375);
        let fasta = write_test_fasta(&seq);
        let genome = ReferenceGenome::load(fasta.path()).unwrap();
        // Request that spans past end of chromosome should fail
        assert!(genome.sequence_at_genome_pos(1490, 20).is_none());
    }

    #[test]
    fn test_reference_genome_sample_sequence_returns_valid() {
        let seq = b"ACGT".repeat(375);
        let fasta = write_test_fasta(&seq);
        let genome = ReferenceGenome::load(fasta.path()).unwrap();
        let mut rng = create_rng(Some(42));
        let result = genome.sample_sequence(100, &mut rng);
        assert!(result.is_some());
        let (chrom_idx, _pos, subseq) = result.unwrap();
        assert_eq!(chrom_idx, 0);
        assert_eq!(subseq.len(), 100);
    }

    #[test]
    fn test_reference_genome_sample_sequence_exact_fit() {
        // length == total_length should succeed (not panic from empty range)
        let seq = b"ACGT".repeat(375); // 1500 bp
        let fasta = write_test_fasta(&seq);
        let genome = ReferenceGenome::load(fasta.path()).unwrap();
        let mut rng = create_rng(Some(42));
        let result = genome.sample_sequence(genome.total_length(), &mut rng);
        assert!(result.is_some());
        let (_chrom_idx, pos, subseq) = result.unwrap();
        assert_eq!(pos, 0);
        assert_eq!(subseq.len(), genome.total_length());
    }

    #[test]
    fn test_reference_genome_sample_sequence_too_large() {
        // length > total_length should return None (not panic)
        let seq = b"ACGT".repeat(375); // 1500 bp
        let fasta = write_test_fasta(&seq);
        let genome = ReferenceGenome::load(fasta.path()).unwrap();
        let mut rng = create_rng(Some(42));
        let result = genome.sample_sequence(genome.total_length() + 1, &mut rng);
        assert!(result.is_none());
    }

    #[test]
    fn test_reference_genome_sample_sequence_zero_length() {
        // length == 0 must not panic and must return None (no valid 0-length sample)
        let seq = b"ACGT".repeat(375);
        let fasta = write_test_fasta(&seq);
        let genome = ReferenceGenome::load(fasta.path()).unwrap();
        let mut rng = create_rng(Some(42));
        let result = genome.sample_sequence(0, &mut rng);
        assert!(result.is_none());
    }

    // ========================================================================
    // MethylationArgs tests
    // ========================================================================

    #[test]
    fn test_methylation_args_defaults() {
        let args = MethylationArgs {
            methylation_mode: None,
            cpg_methylation_rate: 0.75,
            conversion_rate: 0.98,
        };
        assert!((args.cpg_methylation_rate - 0.75).abs() < f64::EPSILON);
        assert!((args.conversion_rate - 0.98).abs() < f64::EPSILON);
    }

    #[test]
    fn test_methylation_args_resolve_disabled() {
        let args = MethylationArgs {
            methylation_mode: None,
            cpg_methylation_rate: 0.75,
            conversion_rate: 0.98,
        };
        assert_eq!(args.resolve().mode, MethylationMode::Disabled);
    }

    #[test]
    fn test_methylation_args_resolve_emseq() {
        let args = MethylationArgs {
            methylation_mode: Some(MethylationModeArg::EmSeq),
            cpg_methylation_rate: 0.75,
            conversion_rate: 0.98,
        };
        let config = args.resolve();
        assert_eq!(config.mode, MethylationMode::EmSeq);
        assert!((config.cpg_methylation_rate - 0.75).abs() < f64::EPSILON);
        assert!((config.conversion_rate - 0.98).abs() < f64::EPSILON);
    }

    #[test]
    fn test_methylation_args_resolve_taps() {
        let args = MethylationArgs {
            methylation_mode: Some(MethylationModeArg::Taps),
            cpg_methylation_rate: 0.75,
            conversion_rate: 0.98,
        };
        assert_eq!(args.resolve().mode, MethylationMode::Taps);
    }

    #[test]
    fn test_methylation_args_validate_valid_rates() {
        for rate in [0.0, 0.5, 1.0] {
            let args = MethylationArgs {
                methylation_mode: None,
                cpg_methylation_rate: rate,
                conversion_rate: rate,
            };
            assert!(args.validate().is_ok(), "rate {rate} should be valid");
        }
    }

    #[rstest]
    #[case(-0.1)]
    #[case(1.1)]
    #[case(f64::NAN)]
    #[case(f64::INFINITY)]
    #[case(f64::NEG_INFINITY)]
    fn test_methylation_args_validate_invalid_cpg_rate(#[case] rate: f64) {
        let args = MethylationArgs {
            methylation_mode: None,
            cpg_methylation_rate: rate,
            conversion_rate: 0.98,
        };
        assert!(args.validate().is_err(), "cpg rate {rate} should be invalid");
    }

    #[rstest]
    #[case(-0.1)]
    #[case(1.1)]
    #[case(f64::NAN)]
    #[case(f64::INFINITY)]
    fn test_methylation_args_validate_invalid_conversion_rate(#[case] rate: f64) {
        let args = MethylationArgs {
            methylation_mode: None,
            cpg_methylation_rate: 0.75,
            conversion_rate: rate,
        };
        assert!(args.validate().is_err(), "conversion rate {rate} should be invalid");
    }

    // ========================================================================
    // apply_methylation_conversion tests
    // ========================================================================

    /// Helper to apply conversion with deterministic rates.
    #[allow(clippy::too_many_arguments)]
    fn convert(
        seq: &mut [u8],
        ref_seq: &[u8],
        ref_offset: usize,
        is_top: bool,
        mode: MethylationMode,
        cpg_rate: f64,
        conv_rate: f64,
        seed: u64,
    ) {
        let config =
            MethylationConfig { mode, cpg_methylation_rate: cpg_rate, conversion_rate: conv_rate };
        let mut rng = create_rng(Some(seed));
        apply_methylation_conversion(seq, ref_seq, ref_offset, is_top, &config, &mut rng);
    }

    #[test]
    fn test_emseq_cpg_all_methylated_no_conversion() {
        // EM-Seq: methylated CpG = protected, should NOT convert
        // cpg_methylation_rate=1.0 means all CpGs are methylated
        let ref_seq = b"ACGTACGT";
        let mut read = ref_seq.to_vec();
        convert(&mut read, ref_seq, 0, true, MethylationMode::EmSeq, 1.0, 1.0, 42);
        // CpG Cs (at positions 1 and 5) should stay as C (methylated = protected in EM-Seq)
        assert_eq!(read[1], b'C', "CpG C should be protected when methylated");
        assert_eq!(read[5], b'C', "CpG C should be protected when methylated");
    }

    #[test]
    fn test_emseq_cpg_all_unmethylated_full_conversion() {
        // EM-Seq: unmethylated CpG = target, should convert C->T
        // cpg_methylation_rate=0.0 means all CpGs are unmethylated
        let ref_seq = b"ACGTACGT";
        let mut read = ref_seq.to_vec();
        convert(&mut read, ref_seq, 0, true, MethylationMode::EmSeq, 0.0, 1.0, 42);
        // CpG Cs should convert to T
        assert_eq!(read[1], b'T', "unmethylated CpG C should convert to T");
        assert_eq!(read[5], b'T', "unmethylated CpG C should convert to T");
    }

    #[test]
    fn test_emseq_non_cpg_c_always_converts() {
        // Non-CpG Cs are unmethylated, always targets in EM-Seq
        // ref = "ACCTA" -> C at pos 1 (non-CpG, followed by C), C at pos 2 (non-CpG, followed by T)
        let ref_seq = b"ACCTA";
        let mut read = ref_seq.to_vec();
        convert(&mut read, ref_seq, 0, true, MethylationMode::EmSeq, 0.75, 1.0, 42);
        assert_eq!(read[1], b'T', "non-CpG C should convert to T in EM-Seq");
        assert_eq!(read[2], b'T', "non-CpG C should convert to T in EM-Seq");
    }

    #[test]
    fn test_taps_cpg_all_methylated_full_conversion() {
        // TAPs: methylated CpG = target, should convert C->T
        let ref_seq = b"ACGTACGT";
        let mut read = ref_seq.to_vec();
        convert(&mut read, ref_seq, 0, true, MethylationMode::Taps, 1.0, 1.0, 42);
        assert_eq!(read[1], b'T', "methylated CpG C should convert in TAPs");
        assert_eq!(read[5], b'T', "methylated CpG C should convert in TAPs");
    }

    #[test]
    fn test_taps_cpg_all_unmethylated_no_conversion() {
        // TAPs: unmethylated CpG = not a target, should stay as C
        let ref_seq = b"ACGTACGT";
        let mut read = ref_seq.to_vec();
        convert(&mut read, ref_seq, 0, true, MethylationMode::Taps, 0.0, 1.0, 42);
        assert_eq!(read[1], b'C', "unmethylated CpG C should not convert in TAPs");
        assert_eq!(read[5], b'C', "unmethylated CpG C should not convert in TAPs");
    }

    #[test]
    fn test_taps_non_cpg_c_never_converts() {
        // Non-CpG Cs are unmethylated, never targets in TAPs
        let ref_seq = b"ACCTA";
        let mut read = ref_seq.to_vec();
        convert(&mut read, ref_seq, 0, true, MethylationMode::Taps, 0.75, 1.0, 42);
        assert_eq!(read[1], b'C', "non-CpG C should not convert in TAPs");
        assert_eq!(read[2], b'C', "non-CpG C should not convert in TAPs");
    }

    #[test]
    fn test_bottom_strand_emseq_converts_g_to_a() {
        // Bottom strand: G at CpG context = unmethylated target in EM-Seq
        // ref = "ACGT" -> G at pos 2, preceded by C -> CpG context
        let ref_seq = b"ACGT";
        let mut read = ref_seq.to_vec();
        convert(&mut read, ref_seq, 0, false, MethylationMode::EmSeq, 0.0, 1.0, 42);
        assert_eq!(read[2], b'A', "bottom strand unmethylated CpG G should convert to A");
    }

    #[test]
    fn test_bottom_strand_taps_converts_g_to_a() {
        // Bottom strand: G at CpG context = methylated target in TAPs
        let ref_seq = b"ACGT";
        let mut read = ref_seq.to_vec();
        convert(&mut read, ref_seq, 0, false, MethylationMode::Taps, 1.0, 1.0, 42);
        assert_eq!(read[2], b'A', "bottom strand methylated CpG G should convert to A in TAPs");
    }

    #[test]
    fn test_non_target_bases_unchanged_top_strand() {
        let ref_seq = b"AGTAGT";
        let mut read = ref_seq.to_vec();
        convert(&mut read, ref_seq, 0, true, MethylationMode::EmSeq, 0.0, 1.0, 42);
        // No Cs in this sequence, nothing should change
        assert_eq!(read, b"AGTAGT");
    }

    #[test]
    fn test_non_target_bases_unchanged_bottom_strand() {
        let ref_seq = b"ACTACT";
        let mut read = ref_seq.to_vec();
        convert(&mut read, ref_seq, 0, false, MethylationMode::EmSeq, 0.0, 1.0, 42);
        // No Gs in this sequence, nothing should change on bottom strand
        assert_eq!(read, b"ACTACT");
    }

    #[test]
    fn test_disabled_mode_no_conversion() {
        let ref_seq = b"ACGTACGT";
        let mut read = ref_seq.to_vec();
        convert(&mut read, ref_seq, 0, true, MethylationMode::Disabled, 0.0, 1.0, 42);
        assert_eq!(read, ref_seq, "Disabled mode should not modify any bases");
    }

    #[test]
    fn test_empty_sequence() {
        let ref_seq = b"";
        let mut read: Vec<u8> = vec![];
        convert(&mut read, ref_seq, 0, true, MethylationMode::EmSeq, 0.75, 0.98, 42);
        assert!(read.is_empty());
    }

    #[test]
    fn test_ref_offset_nonzero() {
        // Read starts at offset 2 in the reference
        let ref_seq = b"AACGTAA";
        //                  ^ offset 2 = C, CpG context
        let mut read = b"CGT".to_vec();
        convert(&mut read, ref_seq, 2, true, MethylationMode::EmSeq, 0.0, 1.0, 42);
        assert_eq!(read[0], b'T', "C at ref_offset=2 (CpG) should convert");
    }

    #[test]
    fn test_conversion_rate_zero_no_conversion() {
        // Even with unmethylated non-CpG C, conversion_rate=0 means no conversion
        let ref_seq = b"ACCTA";
        let mut read = ref_seq.to_vec();
        convert(&mut read, ref_seq, 0, true, MethylationMode::EmSeq, 0.0, 0.0, 42);
        assert_eq!(read[1], b'C', "conversion_rate=0 should prevent conversion");
        assert_eq!(read[2], b'C', "conversion_rate=0 should prevent conversion");
    }

    #[test]
    fn test_probabilistic_emseq_cpg_partial_methylation() {
        // With cpg_methylation_rate=0.5, roughly half of CpG Cs should convert
        let ref_seq = b"CG"; // single CpG
        let mut converted_count = 0;
        let trials = 10_000;
        for seed in 0..trials {
            let mut read = ref_seq.to_vec();
            convert(&mut read, ref_seq, 0, true, MethylationMode::EmSeq, 0.5, 1.0, seed);
            if read[0] == b'T' {
                converted_count += 1;
            }
        }
        // Expected: ~50% convert (unmethylated) with conversion_rate=1.0
        let fraction = converted_count as f64 / trials as f64;
        assert!(
            (fraction - 0.5).abs() < 0.05,
            "Expected ~50% conversion at CpG with methylation_rate=0.5, got {fraction:.3}"
        );
    }

    #[test]
    fn test_probabilistic_emseq_non_cpg_partial_conversion_rate() {
        // Non-CpG C with conversion_rate=0.5 should convert ~50% of the time
        let ref_seq = b"ACT"; // C at pos 1, non-CpG (followed by T)
        let mut converted_count = 0;
        let trials = 10_000;
        for seed in 0..trials {
            let mut read = ref_seq.to_vec();
            convert(&mut read, ref_seq, 0, true, MethylationMode::EmSeq, 0.75, 0.5, seed);
            if read[1] == b'T' {
                converted_count += 1;
            }
        }
        let fraction = converted_count as f64 / trials as f64;
        assert!(
            (fraction - 0.5).abs() < 0.05,
            "Expected ~50% conversion with conversion_rate=0.5, got {fraction:.3}"
        );
    }

    #[test]
    fn test_conversion_rate_zero_leaves_bases_unchanged() {
        let ref_seq = b"CACACACACACACACAC"; // non-CpG Cs
        let mut read = ref_seq.to_vec();
        convert(&mut read, ref_seq, 0, true, MethylationMode::EmSeq, 0.0, 0.0, 42);
        assert_eq!(read, ref_seq, "conversion_rate=0 should leave all bases unchanged");
    }

    #[test]
    fn test_conversion_rate_one_converts_all_targets() {
        // EM-Seq, cpg_methylation_rate=0 means all CpGs unmethylated -> all Cs are targets
        let ref_seq = b"CACACACACACACACAC"; // all non-CpG Cs
        let mut read = ref_seq.to_vec();
        convert(&mut read, ref_seq, 0, true, MethylationMode::EmSeq, 0.0, 1.0, 42);
        for (i, &b) in read.iter().enumerate() {
            if ref_seq[i] == b'C' {
                assert_eq!(b, b'T', "position {i}: C should be converted with rate=1.0");
            } else {
                assert_eq!(b, ref_seq[i], "position {i}: non-C should be unchanged");
            }
        }
    }

    #[test]
    fn test_disabled_mode_never_converts() {
        let ref_seq = b"CACGTCACGTCACGT";
        let mut read = ref_seq.to_vec();
        convert(&mut read, ref_seq, 0, true, MethylationMode::Disabled, 0.75, 1.0, 42);
        assert_eq!(read, ref_seq, "Disabled mode should never modify bases");
    }

    // ========================================================================
    // ReferenceGenome::build_bam_header, num_chromosomes, chromosome_length,
    // and sample_positions tests
    // ========================================================================

    /// Create a temp FASTA with four contigs: chr1 (2000bp), chr2 (1500bp),
    /// chr3 (1800bp), and `short_contig` (500bp, below the 1000bp minimum).
    fn write_multi_contig_fasta() -> NamedTempFile {
        let mut f = NamedTempFile::new().unwrap();
        writeln!(f, ">chr1").unwrap();
        f.write_all(&b"ACGT".repeat(500)).unwrap(); // 2000bp
        writeln!(f).unwrap();
        writeln!(f, ">chr2").unwrap();
        f.write_all(&b"CCGG".repeat(375)).unwrap(); // 1500bp
        writeln!(f).unwrap();
        writeln!(f, ">chr3").unwrap();
        f.write_all(&b"AATT".repeat(450)).unwrap(); // 1800bp
        writeln!(f).unwrap();
        writeln!(f, ">short_contig").unwrap();
        f.write_all(&b"ACGT".repeat(125)).unwrap(); // 500bp, should be skipped
        writeln!(f).unwrap();
        f.flush().unwrap();
        f
    }

    #[test]
    fn test_build_bam_header_has_all_contigs() {
        let fasta = write_multi_contig_fasta();
        let genome = ReferenceGenome::load(fasta.path()).unwrap();
        let header = genome.build_bam_header();

        let ref_seqs: Vec<_> = header.reference_sequences().keys().collect();
        assert_eq!(ref_seqs.len(), 3, "short_contig should be excluded");
        let names: Vec<&str> =
            ref_seqs.iter().map(|k| std::str::from_utf8(k.as_ref()).unwrap()).collect();
        assert!(names.contains(&"chr1"));
        assert!(names.contains(&"chr2"));
        assert!(names.contains(&"chr3"));
        assert!(!names.contains(&"short_contig"));
    }

    #[test]
    fn test_build_bam_header_contig_lengths() {
        let fasta = write_multi_contig_fasta();
        let genome = ReferenceGenome::load(fasta.path()).unwrap();
        let header = genome.build_bam_header();

        let ref_seqs = header.reference_sequences();
        let chr1_len: usize = ref_seqs.get(&bstr::BString::from("chr1")).unwrap().length().get();
        let chr2_len: usize = ref_seqs.get(&bstr::BString::from("chr2")).unwrap().length().get();
        let chr3_len: usize = ref_seqs.get(&bstr::BString::from("chr3")).unwrap().length().get();
        assert_eq!(chr1_len, 2000);
        assert_eq!(chr2_len, 1500);
        assert_eq!(chr3_len, 1800);
    }

    #[test]
    fn test_num_chromosomes() {
        let fasta = write_multi_contig_fasta();
        let genome = ReferenceGenome::load(fasta.path()).unwrap();
        assert_eq!(genome.num_chromosomes(), 3);
    }

    #[test]
    fn test_chromosome_length() {
        let fasta = write_multi_contig_fasta();
        let genome = ReferenceGenome::load(fasta.path()).unwrap();
        assert_eq!(genome.chromosome_length(0), 2000);
        assert_eq!(genome.chromosome_length(1), 1500);
        assert_eq!(genome.chromosome_length(2), 1800);
    }

    #[test]
    fn test_sample_positions_count_and_bounds() {
        let fasta = write_multi_contig_fasta();
        let genome = ReferenceGenome::load(fasta.path()).unwrap();
        let mut rng = create_rng(Some(42));
        let positions = genome.sample_positions(20, &mut rng);
        assert_eq!(positions.len(), 20);
        for (chrom_idx, local_pos) in &positions {
            assert!(*chrom_idx < genome.num_chromosomes());
            assert!(*local_pos < genome.chromosome_length(*chrom_idx));
        }
    }

    #[test]
    fn test_sample_positions_deterministic() {
        let fasta = write_multi_contig_fasta();
        let genome = ReferenceGenome::load(fasta.path()).unwrap();
        let mut rng1 = create_rng(Some(99));
        let mut rng2 = create_rng(Some(99));
        let pos1 = genome.sample_positions(10, &mut rng1);
        let pos2 = genome.sample_positions(10, &mut rng2);
        assert_eq!(pos1, pos2);
    }

    #[test]
    fn test_sample_positions_spans_chromosomes() {
        let fasta = write_multi_contig_fasta();
        let genome = ReferenceGenome::load(fasta.path()).unwrap();
        let mut rng = create_rng(Some(7));
        let positions = genome.sample_positions(100, &mut rng);
        let unique_chroms: std::collections::HashSet<usize> =
            positions.iter().map(|(c, _)| *c).collect();
        assert!(
            unique_chroms.len() >= 2,
            "100 positions should span at least 2 chromosomes, got {}",
            unique_chroms.len()
        );
    }
}

//! Shared CLI arguments and utilities for simulation commands.

use crate::commands::common::{
    MemoryLimit, MemoryReserve, MethylationModeArg, parse_memory, parse_memory_reserve,
    resolve_memory_budget,
};
use crate::commands::sort::{TMP_DIRS_ENV, resolve_tmp_dirs};
use anyhow::{Context, Result, anyhow, bail};
use bytesize::ByteSize;
use clap::Args;
use crossbeam_channel::{Receiver, Sender, bounded};
use fgumi_consensus::MethylationMode;
use fgumi_consensus::methylation::is_cpg_context;
use fgumi_raw_bam::RawRecord;
use fgumi_sort::{RawExternalSorter, SortOrder};
use log::info;
use noodles::fasta;
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::header::{self as HeaderRecord, Tag as HeaderTag};
use rand::{Rng, RngExt};
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;
use std::thread::JoinHandle;

/// Build the `@HD` header map carrying the `SO`/`GO`/`SS` sort-order tags for `order`.
///
/// Shared by the `simulate` subcommands so they emit the same `SO`/`GO`/`SS` values —
/// including the `SS = <sort-order>:<sub-sort>` prefix — that `fgumi sort` writes (see
/// [`SortOrder::header_ss_tag`]). `GO` and `SS` are inserted only when the order defines
/// them, so this is safe for any [`SortOrder`], not just template-coordinate.
#[must_use]
pub fn build_sort_order_header_map(order: SortOrder) -> Map<HeaderRecord::Header> {
    let HeaderTag::Other(so_tag) = HeaderTag::from([b'S', b'O']) else { unreachable!() };
    let HeaderTag::Other(go_tag) = HeaderTag::from([b'G', b'O']) else { unreachable!() };
    let HeaderTag::Other(ss_tag) = HeaderTag::from([b'S', b'S']) else { unreachable!() };

    let mut builder = Map::<HeaderRecord::Header>::builder().insert(so_tag, order.header_so_tag());
    if let Some(go) = order.header_go_tag() {
        builder = builder.insert(go_tag, go);
    }
    if let Some(ss) = order.header_ss_tag() {
        builder = builder.insert(ss_tag, ss);
    }
    builder.build().expect("header map with valid SO/GO/SS tags")
}

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
    /// Convert to a [`PositionQualityModel`](crate::simulate::PositionQualityModel).
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

    /// Convert to a [`ReadPairQualityBias`](crate::simulate::ReadPairQualityBias).
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
    /// Convert to an [`InsertSizeModel`](crate::simulate::InsertSizeModel).
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
    /// Convert to a [`FamilySizeDistribution`](crate::simulate::FamilySizeDistribution).
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
    /// Convert to a [`StrandBiasModel`](crate::simulate::StrandBiasModel).
    pub fn to_strand_bias_model(&self) -> crate::simulate::StrandBiasModel {
        crate::simulate::StrandBiasModel::new(self.strand_alpha, self.strand_beta)
    }
}

/// Methylation simulation options shared across simulate subcommands.
///
/// Models a directional EM-seq/TAPs library as `holodeck simulate`/`methylate` do: each
/// `CpG` has a methylation state fixed for the run per strand (non-`CpG` cytosines are
/// unmethylated), and chemistry converts each original strand of a molecule once, so all
/// reads of a family and both mates of a pair share it.
#[derive(Args, Debug, Clone)]
pub struct MethylationArgs {
    /// Methylation chemistry mode. When set, enables methylation-aware base
    /// conversion in simulated reads. Requires --reference for commands that
    /// generate read bases (fastq-reads, mapped-reads, grouped-reads).
    #[arg(long = "methylation-mode", value_enum)]
    pub methylation_mode: Option<MethylationModeArg>,

    /// Fraction of `CpG`s that are methylated [0.0-1.0], drawn once per `CpG` for the run
    /// (from --seed). Methylated `CpG`s are protected from conversion in EM-Seq and are
    /// targets for conversion in TAPs.
    #[arg(long = "cpg-methylation-rate", default_value = "0.75")]
    pub cpg_methylation_rate: f64,

    /// Probability that a methylated `CpG` is methylated on one strand only [0.0-1.0].
    #[arg(long = "hemimethylation-rate", default_value = "0.01")]
    pub hemimethylation_rate: f64,

    /// Enzymatic conversion efficiency for target cytosines [0.0-1.0].
    /// In EM-Seq, this is the probability that an unmethylated C is converted to T.
    /// In TAPs, this is the probability that a methylated C is converted to T.
    #[arg(
        long = "methylation-conversion-rate",
        alias = "conversion-rate",
        default_value = "0.999"
    )]
    pub conversion_rate: f64,

    /// Fraction of molecule strands whose conversion fails as a whole [0.0-1.0]; a failed
    /// strand converts at 1 - conversion rate, keeping almost all its cytosines.
    #[arg(long = "methylation-failure-rate", default_value = "0.01")]
    pub failure_rate: f64,
}

impl MethylationArgs {
    /// Resolves the optional CLI arg to a [`MethylationConfig`]. The fixed per-`CpG`
    /// methylation state is seeded from `seed` (the run's `--seed`), or randomly when unset.
    pub fn resolve(&self, seed: Option<u64>) -> MethylationConfig {
        MethylationConfig {
            mode: crate::commands::common::resolve_methylation_mode(self.methylation_mode),
            cpg_methylation_rate: self.cpg_methylation_rate,
            conversion_rate: self.conversion_rate,
            hemimethylation_rate: self.hemimethylation_rate,
            failure_rate: self.failure_rate,
            table_seed: seed.unwrap_or_else(rand::random),
        }
    }

    /// Validates that rate parameters are in [0.0, 1.0] and finite.
    pub fn validate(&self) -> anyhow::Result<()> {
        validate_rate(self.cpg_methylation_rate, "cpg-methylation-rate")?;
        validate_rate(self.hemimethylation_rate, "hemimethylation-rate")?;
        validate_rate(self.conversion_rate, "methylation-conversion-rate")?;
        validate_rate(self.failure_rate, "methylation-failure-rate")?;
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
    /// Probability that a methylated `CpG` is methylated on one strand only.
    pub hemimethylation_rate: f64,
    /// Fraction of molecule strands whose conversion fails as a whole.
    pub failure_rate: f64,
    /// Seed of the run's fixed per-`CpG` methylation state.
    pub table_seed: u64,
}

/// A noise-free test baseline: no hemimethylation, no conversion failures and a fixed table
/// seed, unlike the CLI defaults (see [`MethylationArgs`]).
#[cfg(test)]
impl Default for MethylationConfig {
    fn default() -> Self {
        Self {
            mode: MethylationMode::Disabled,
            cpg_methylation_rate: 0.75,
            conversion_rate: 0.999,
            hemimethylation_rate: 0.0,
            failure_rate: 0.0,
            table_seed: 0,
        }
    }
}

/// Validates that a rate is a finite value in [0.0, 1.0].
pub(super) fn validate_rate(value: f64, name: &str) -> anyhow::Result<()> {
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

    /// Returns the whole (uppercased) sequence of the chromosome at the given index.
    pub fn contig(&self, chrom_idx: usize) -> &[u8] {
        &self.sequences[chrom_idx]
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
    #[allow(dead_code)] // used only in tests
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

    /// Build a BAM `Header` with `@SQ` lines for every loaded contig.
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
    #[allow(dead_code)] // used only in tests
    pub(super) fn num_chromosomes(&self) -> usize {
        self.sequences.len()
    }

    /// Returns the length of the chromosome at the given index.
    #[allow(dead_code)] // used only in tests
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

/// Sort-engine resource options for the simulation commands that stream their
/// records through `fgumi-sort`.
///
/// These mirror the equivalent `fgumi sort` flags (same names, same defaults)
/// so a simulation that has to spill behaves — and is tuned — exactly like a
/// standalone sort.
#[derive(Args, Debug, Clone)]
pub struct SortResourceArgs {
    /// Maximum memory for in-memory sorting of the generated records.
    ///
    /// Default is "768M" per thread (matching `fgumi sort` and samtools). Pass
    /// "auto" to detect system memory and subtract --memory-reserve. When the
    /// limit is reached, sorted chunks spill to temporary files.
    #[arg(short = 'm', long = "max-memory", default_value = "768M", value_parser = parse_memory)]
    pub max_memory: MemoryLimit,

    /// Memory to reserve for other processes when --max-memory=auto.
    ///
    /// Ignored when --max-memory is set to an explicit value.
    #[arg(long = "memory-reserve", default_value = "auto", value_parser = parse_memory_reserve)]
    pub memory_reserve: MemoryReserve,

    /// Scale the memory limit by thread count (samtools behavior).
    ///
    /// When enabled (default), --max-memory specifies memory per thread.
    #[arg(long = "memory-per-thread", value_name = "true|false", default_value = "true", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = clap::builder::BoolishValueParser::new(), hide_possible_values = true)]
    pub memory_per_thread: bool,

    /// Temporary directory for sort spill files. Repeatable.
    ///
    /// If no flags are given and `FGUMI_TMP_DIRS` is set, its value is parsed as
    /// a `PATH`-style list and used instead. If neither is provided, the system
    /// default temp directory is used — note that on many modern Linux distros
    /// `/tmp` is a RAM-backed tmpfs, so point this at real disk for large runs.
    #[arg(short = 'T', long = "tmp-dir", action = clap::ArgAction::Append)]
    pub tmp_dirs: Vec<PathBuf>,
}

impl SortResourceArgs {
    /// Apply the resolved memory budget and temp directories to `sorter`.
    ///
    /// `threads` is the command's thread count, which scales the per-thread
    /// memory budget exactly as it does in `fgumi sort`.
    ///
    /// # Errors
    ///
    /// Returns an error if the memory budget cannot be resolved (e.g. zero
    /// threads, or an overflowing per-thread limit).
    pub fn apply(&self, sorter: RawExternalSorter, threads: usize) -> Result<RawExternalSorter> {
        let effective_memory = resolve_memory_budget(
            self.max_memory,
            self.memory_reserve,
            threads,
            self.memory_per_thread,
        )?;
        info!("Sort memory: {}", ByteSize(effective_memory as u64));

        let mut sorter = sorter.memory_limit(effective_memory);

        // For auto mode, cap the initial buffer pre-allocation at 768 MiB/thread
        // so an `auto` budget on a large host doesn't allocate the whole budget
        // upfront for what may be a tiny simulation. Mirrors `fgumi sort`.
        if matches!(self.max_memory, MemoryLimit::Auto) {
            let init = 768_usize
                .checked_mul(1024 * 1024)
                .and_then(|b| b.checked_mul(threads))
                .ok_or_else(|| anyhow::anyhow!("initial auto buffer size overflowed"))?;
            sorter = sorter.initial_capacity(effective_memory.min(init));
        }

        let env_value = std::env::var(TMP_DIRS_ENV).ok();
        let tmp_dirs = resolve_tmp_dirs(&self.tmp_dirs, env_value.as_deref());
        if !tmp_dirs.is_empty() {
            let joined =
                tmp_dirs.iter().map(|p| p.display().to_string()).collect::<Vec<_>>().join(", ");
            info!("Sort temp directories: {joined}");
            sorter = sorter.temp_dirs(tmp_dirs);
        }

        Ok(sorter)
    }
}

/// Number of record *batches* in flight between the generator thread and the
/// sort's accumulate loop.
///
/// With `RECORD_BATCH_SIZE` this buffers ~16k records (a few MiB). Deep enough
/// that the generator is not stalled by short pauses in the sort (e.g. a
/// spill), small enough that the queue itself is not a meaningful memory
/// consumer next to the sort buffer.
const RECORD_CHANNEL_CAPACITY: usize = 64;

/// Records per batch handed across the channel.
///
/// Records move in batches, not one at a time: a per-record send costs more
/// than the whole intermediate-BAM round-trip this streaming design removes
/// (measured: ~7% slower end-to-end when unbatched). 256 matches the batch size
/// the sort's own read-ahead reader uses.
const RECORD_BATCH_SIZE: usize = 256;

/// A batch of generated records, or the generator's terminal error.
type RecordBatch = Result<Vec<RawRecord>>;

/// Sending half of the generator → sorter record stream.
///
/// The simulation commands generate records on a background thread and hand
/// them straight to [`RawExternalSorter::sort_records`], so no unsorted
/// intermediate BAM is ever written: runs that fit in memory touch disk only
/// for the final output, and runs that don't spill exactly like `fgumi sort`.
pub(super) struct RecordSink {
    sender: Sender<RecordBatch>,
    batch: Vec<RawRecord>,
}

impl RecordSink {
    /// Create a sink/stream pair to connect a generator thread to the sorter.
    ///
    /// Pass the receiver through [`into_record_stream`] to get the per-record
    /// iterator [`RawExternalSorter::sort_records`] takes; the sink goes to the
    /// generator.
    ///
    /// The generator must **own** the sink and drop it when it finishes:
    /// dropping the sender is what signals end-of-stream. A sink kept alive
    /// elsewhere (e.g. left in the caller's scope) leaves the sort blocked
    /// waiting for records that will never arrive. Call [`Self::finish`] to
    /// deliver the last partial batch; a sink dropped with records still
    /// buffered fails the sort rather than silently truncating it (see the
    /// [`Drop`] impl).
    pub(super) fn new() -> (Self, Receiver<RecordBatch>) {
        let (sender, receiver) = bounded(RECORD_CHANNEL_CAPACITY);
        (Self { sender, batch: Vec::with_capacity(RECORD_BATCH_SIZE) }, receiver)
    }

    /// Hand one record to the sorter, sending it once a batch has accumulated.
    ///
    /// Returns `false` once the sorter has hung up, which only happens when the
    /// sort itself failed. Generators must stop on `false` — continuing would
    /// buffer forever against a receiver that is gone — and let the sort's own
    /// error surface as the reported failure.
    pub(super) fn send(&mut self, record: RawRecord) -> bool {
        self.batch.push(record);
        if self.batch.len() < RECORD_BATCH_SIZE {
            return true;
        }
        self.flush()
    }

    /// Send any buffered records, ending the stream cleanly.
    ///
    /// Returns `false` if the sorter had already hung up, exactly as
    /// [`Self::send`] does.
    pub(super) fn finish(mut self) -> bool {
        self.flush()
    }

    /// Report a generator failure to the sorter, which aborts the sort and
    /// propagates the error instead of writing a truncated output.
    ///
    /// Buffered records are dropped: the stream is being abandoned, and sending
    /// them would only delay the abort.
    pub(super) fn fail(mut self, error: anyhow::Error) {
        // Deliberate: this abort is the error the caller should see, so clear
        // the buffer before dropping rather than let `Drop` report a second one.
        self.batch.clear();
        // If the sorter already hung up it is failing on its own error, which
        // is the one the user should see; dropping ours is correct.
        let _ = self.sender.send(Err(error));
    }

    /// Send the current batch if it is non-empty, leaving the buffer empty.
    fn flush(&mut self) -> bool {
        if self.batch.is_empty() {
            return true;
        }
        let batch = std::mem::replace(&mut self.batch, Vec::with_capacity(RECORD_BATCH_SIZE));
        self.sender.send(Ok(batch)).is_ok()
    }
}

impl Drop for RecordSink {
    /// Report records still buffered when the sink is abandoned.
    ///
    /// Dropping the sender is how end-of-stream is signalled, so a sink dropped
    /// with a partial batch — a `?` on the generator's path, or a panic —
    /// looked exactly like a clean end of input: the sort finished and wrote a
    /// complete, short output. Sending the loss as the stream's terminal error
    /// makes the sort abort instead, so no truncated output is written.
    ///
    /// [`Self::finish`] and [`Self::fail`] both leave the buffer empty, so this
    /// is a no-op after either of them.
    fn drop(&mut self) {
        if self.batch.is_empty() {
            return;
        }
        let unflushed = self.batch.len();
        // A closed receiver means the sort already failed on its own error,
        // which is the one to report; ours would only mask it.
        let _ = self.sender.send(Err(anyhow!(
            "The record generator was abandoned with {unflushed} generated record(s) still \
             buffered; they were never handed to the sort"
        )));
    }
}

/// Adapt a [`RecordSink`]'s batched receiver into the per-record stream
/// [`RawExternalSorter::sort_records`] consumes.
///
/// Batching is a transport detail of the generator → sorter handoff, so it is
/// unwrapped here rather than pushed into the sort's API.
pub(super) fn into_record_stream(
    receiver: Receiver<RecordBatch>,
) -> impl Iterator<Item = Result<RawRecord>> + Send + 'static {
    receiver.into_iter().flat_map(|batch| -> Box<dyn Iterator<Item = Result<RawRecord>> + Send> {
        // One `Box` per batch, not per record.
        match batch {
            Ok(records) => Box::new(records.into_iter().map(Ok)),
            Err(e) => Box::new(std::iter::once(Err(e))),
        }
    })
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

/// Introduce random substitution errors into `seq` in place.
///
/// Uses a direct alternative-base lookup table to avoid retry loops. When `error_rate <= 0.0` this
/// returns early without drawing any RNG, so a zero (or non-positive) rate leaves the RNG stream
/// untouched and callers may invoke it unconditionally. For `error_rate > 0.0` it draws exactly one
/// `rng.random::<f64>()` uniform per base and substitutes that base when the draw falls below the
/// rate.
pub(super) fn introduce_errors_inplace(seq: &mut [u8], error_rate: f64, rng: &mut impl Rng) {
    // Fast-path/footgun guard: at a non-positive rate there is nothing to substitute, and
    // returning early keeps the RNG stream untouched (no draws) so a zero rate is deterministic
    // and safe to call unconditionally.
    if error_rate <= 0.0 {
        return;
    }

    // Lookup table: for each base, provides 3 alternatives.
    // A -> C,G,T; C -> A,G,T; G -> A,C,T; T -> A,C,G
    const ALTERNATIVES: [&[u8; 3]; 256] = {
        let mut table: [&[u8; 3]; 256] = [b"ACG"; 256]; // Default (shouldn't be used)
        table[b'A' as usize] = b"CGT";
        table[b'C' as usize] = b"AGT";
        table[b'G' as usize] = b"ACT";
        table[b'T' as usize] = b"ACG";
        table
    };

    for base in seq.iter_mut() {
        if rng.random::<f64>() < error_rate {
            let alts = ALTERNATIVES[*base as usize];
            *base = alts[rng.random_range(0..3)];
        }
    }
}

/// Derive a deterministic RNG dedicated to body-error injection for one mate of one read.
///
/// Body-error injection ([`introduce_errors_inplace`]) must not draw from the
/// molecule-generation RNG. If it did, enabling `--error-rate > 0` would advance the shared
/// stream and thereby perturb read qualities, mate/strand assignments, and every later
/// stochastic draw — so an error-injected run would differ from the error-free run in far more
/// than the read bodies. Seeding a dedicated stream keeps the molecule RNG byte-identical whether
/// or not errors are injected, so the error-injected output differs from the error-free output
/// *only* in the read-body bases.
///
/// The seed is an FNV-1a hash of the globally-unique `read_name` folded with `is_first_mate` (which
/// simply distinguishes a pair's two mates so each gets an independent error stream). It consumes
/// no draws from the caller's RNG, so the caller's stream is unaffected by the presence of errors.
pub(super) fn body_error_rng(read_name: &str, is_first_mate: bool) -> rand::rngs::StdRng {
    const FNV_OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
    const FNV_PRIME: u64 = 0x0000_0100_0000_01b3;
    let mut hash = FNV_OFFSET;
    for &byte in read_name.as_bytes() {
        hash = (hash ^ u64::from(byte)).wrapping_mul(FNV_PRIME);
    }
    hash = (hash ^ u64::from(is_first_mate)).wrapping_mul(FNV_PRIME);
    crate::simulate::create_rng(Some(hash))
}

/// Pad a sequence to the target length with random bases, or truncate if too long.
pub(super) fn pad_sequence(mut seq: Vec<u8>, target_len: usize, rng: &mut impl Rng) -> Vec<u8> {
    while seq.len() < target_len {
        seq.push(b"ACGT"[rng.random_range(0..4)]);
    }
    seq.truncate(target_len);
    seq
}

/// Where a template sits in the reference, for looking up `CpG` context and methylation state.
#[derive(Debug, Clone, Copy)]
pub(super) struct TemplateLocus<'a> {
    /// Index of the contig, or `None` for a template not drawn from the reference.
    pub chrom_idx: Option<usize>,
    /// The whole contig (or the template itself when not drawn from the reference).
    pub contig: &'a [u8],
    /// 0-based offset of the template's first base in `contig`.
    pub start: usize,
}

impl<'a> TemplateLocus<'a> {
    /// A template not drawn from the reference: its own bases are the context.
    pub fn standalone(template: &'a [u8]) -> Self {
        Self { chrom_idx: None, contig: template, start: 0 }
    }
}

/// Converts one strand of a molecule, returning the converted template in genomic orientation.
///
/// Models a directional EM-seq/bisulfite/TAPs library as `holodeck` does: chemistry acts once
/// on the original strand the molecule's reads derive from, so every read of a family (PCR
/// copies of that strand) shares the result, and both mates of a pair are cut from it. The
/// top strand's convertible cytosines are reference `C` (converted to `T`); the bottom
/// strand's are reference `G` (converted to `A` in genomic orientation).
///
/// Whether a cytosine converts depends on its methylation state, fixed for the run per `CpG`
/// per strand ([`MethylationConfig::is_methylated`]; non-`CpG` cytosines are unmethylated), and
/// on the chemistry: EM-seq converts unmethylated cytosines, TAPs methylated ones, each with
/// probability `conversion_rate`. A molecule strand drawn as a conversion failure (probability
/// `failure_rate`) converts at `1 - conversion_rate` instead. Returns the template unchanged,
/// without drawing from `rng`, when methylation is disabled.
pub(super) fn convert_molecule_strand(
    template: &[u8],
    locus: TemplateLocus<'_>,
    is_top_strand: bool,
    config: &MethylationConfig,
    rng: &mut impl Rng,
) -> Vec<u8> {
    let mut converted = template.to_vec();
    if !config.mode.is_enabled() {
        return converted;
    }
    let failed = config.failure_rate > 0.0 && rng.random::<f64>() < config.failure_rate;
    let rate = if failed { 1.0 - config.conversion_rate } else { config.conversion_rate };
    let (target, converted_base) = if is_top_strand { (b'C', b'T') } else { (b'G', b'A') };
    for (i, base) in converted.iter_mut().enumerate() {
        if *base != target {
            continue;
        }
        let methylated = config.is_methylated(locus, locus.start + i, is_top_strand);
        let should_convert = match config.mode {
            MethylationMode::EmSeq => !methylated,
            MethylationMode::Taps => methylated,
            MethylationMode::Disabled => false,
        };
        if should_convert && rng.random::<f64>() < rate {
            *base = converted_base;
        }
    }
    converted
}

impl MethylationConfig {
    /// Logs the methylation model's settings.
    pub(super) fn log_settings(&self) {
        log::info!("  Methylation mode: {:?}", self.mode);
        log::info!("  CpG methylation rate: {}", self.cpg_methylation_rate);
        log::info!("  Hemimethylation rate: {}", self.hemimethylation_rate);
        log::info!("  Conversion rate: {}", self.conversion_rate);
        log::info!("  Conversion failure rate: {}", self.failure_rate);
        log::info!("  Methylation table seed: {}", self.table_seed);
    }

    /// Whether the cytosine at `pos` of `locus.contig` on the given strand is methylated: the
    /// `C` of a `CpG` on the top strand, the `G` on the bottom strand. Non-`CpG` cytosines are
    /// unmethylated. A `CpG` is methylated with probability `cpg_methylation_rate`, on both
    /// strands unless hemimethylated (probability `hemimethylation_rate`, one strand chosen at
    /// random). The draws are a pure function of `table_seed`, the contig and the `CpG`, so the
    /// state is fixed for the run.
    pub(super) fn is_methylated(
        &self,
        locus: TemplateLocus<'_>,
        pos: usize,
        is_top_strand: bool,
    ) -> bool {
        if !is_cpg_context(locus.contig, pos, is_top_strand) {
            return false;
        }
        let cpg = if is_top_strand { pos } else { pos - 1 };
        // A template not drawn from the reference is keyed by its own bases, so two random
        // templates do not share one methylation pattern; the top bit keeps the key apart from
        // contig indices.
        let chrom = locus.chrom_idx.map_or_else(
            || locus.contig.iter().fold(0, |h, &b| splitmix64(h ^ u64::from(b))) | 1 << 63,
            |c| c as u64,
        );
        let draw = |salt: u64| {
            // Hash the position before mixing in the salt: `cpg ^ salt` would alias the draw for
            // one salt at one CpG onto another salt's draw at a nearby CpG.
            let h = splitmix64(
                self.table_seed ^ splitmix64(chrom ^ splitmix64(splitmix64(cpg as u64) ^ salt)),
            );
            #[expect(clippy::cast_precision_loss, reason = "53-bit uniform in [0, 1)")]
            let u = (h >> 11) as f64 / (1u64 << 53) as f64;
            u
        };
        if draw(0) >= self.cpg_methylation_rate {
            return false;
        }
        if draw(1) < self.hemimethylation_rate {
            // Hemimethylated: methylated on one strand only.
            let methylated_top = draw(2) < 0.5;
            return methylated_top == is_top_strand;
        }
        true
    }
}

/// `SplitMix64` finalizer: a fast, well-mixed 64-bit hash.
fn splitmix64(x: u64) -> u64 {
    let mut z = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Join the writer thread, then decide which of the two failures to report.
///
/// Every `simulate` subcommand that streams records to a writer thread ends the
/// same way: the generation pass returns a channel `SendError` or `Ok(())`, and
/// the writer thread returns the counts it accumulated. Joining is unconditional
/// because the writer owns the output file — returning while it is still running
/// leaves the file being written behind the error.
///
/// A send failure is a *symptom*: the channel only closes because the writer
/// already died, so the writer's error is the cause and is reported in preference
/// to it. The send failure is reported only when the writer itself succeeded,
/// which means the receiver was dropped for some other reason.
///
/// # Errors
///
/// Returns an error if the writer thread panicked, if the writer thread returned
/// an error, or — only when the writer succeeded — if generation failed to send.
pub fn join_writer_result<T, E: std::fmt::Display>(
    writer_handle: JoinHandle<Result<T>>,
    generation_result: std::result::Result<(), E>,
) -> Result<T> {
    let writer_result = writer_handle.join().map_err(|_| anyhow!("Writer thread panicked"))?;

    if let Err(e) = generation_result {
        writer_result?;
        return Err(anyhow!("Failed to send record to writer: {e}"));
    }

    writer_result
}

/// Converts a read's bases and qualities from read (sequencing) orientation to the reference
/// orientation BAM stores: reverse-complemented and reversed for a reverse-strand record,
/// unchanged (borrowed) otherwise.
pub(super) fn to_reference_orientation<'a>(
    seq: &'a [u8],
    quals: &'a [u8],
    is_reverse: bool,
) -> (std::borrow::Cow<'a, [u8]>, std::borrow::Cow<'a, [u8]>) {
    use std::borrow::Cow;
    if is_reverse {
        let quals: Vec<u8> = quals.iter().rev().copied().collect();
        (Cow::Owned(crate::dna::reverse_complement(seq)), Cow::Owned(quals))
    } else {
        (Cow::Borrowed(seq), Cow::Borrowed(quals))
    }
}

/// BAM-encoded CIGAR of a simulated mapped read of `read_len` bases whose first `aligned_len`
/// bases (in read orientation) come from the template and the rest is read-through padding,
/// as when the insert is shorter than the read. The padding sits at the read's 3' end, which
/// is the right end of a forward record and the left end of a reverse one, and is soft-clipped
/// so the aligned bases match the reference at the record's position.
pub(super) fn read_through_cigar(
    read_len: usize,
    aligned_len: usize,
    is_reverse: bool,
) -> Vec<u32> {
    // BAM CIGAR encoding: (length << 4) | op_code; 0 = M (alignment match), 4 = S (soft clip).
    let op =
        |len: usize, code: u32| (u32::try_from(len).expect("read length fits u32") << 4) | code;
    let aligned_len = aligned_len.min(read_len);
    let clipped = read_len - aligned_len;
    match (clipped, is_reverse) {
        (0, _) => vec![op(aligned_len, 0)],
        (_, false) => vec![op(aligned_len, 0), op(clipped, 4)],
        (_, true) => vec![op(clipped, 4), op(aligned_len, 0)],
    }
}

/// Renders BAM-encoded CIGAR ops (as produced by [`read_through_cigar`]) as a SAM CIGAR string.
pub(super) fn cigar_string(ops: &[u32]) -> String {
    use std::fmt::Write;
    ops.iter().fold(String::new(), |mut cigar, op| {
        let code = if op & 0xF == 4 { 'S' } else { 'M' };
        write!(cigar, "{}{code}", op >> 4).expect("write to String is infallible");
        cigar
    })
}

/// Test fixtures shared by the `simulate` generators.
#[cfg(test)]
pub(super) mod test_support {
    use std::io::Write;

    /// A 2 kb single-contig (`chr1`) reference whose reverse complement differs from itself,
    /// so a reverse-strand record carrying reverse-complemented SEQ cannot match it by accident
    /// (unlike `ACGT` repeats, which are their own reverse complement). Returns the FASTA and
    /// the contig bases.
    pub(in crate::commands::simulate) fn asymmetric_reference() -> (tempfile::NamedTempFile, Vec<u8>)
    {
        let mut state: u32 = 0x1234_5678;
        let bases: Vec<u8> = (0..2000)
            .map(|_| {
                state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                b"ACGT"[(state >> 30) as usize]
            })
            .collect();
        let mut fasta = tempfile::NamedTempFile::new().unwrap();
        writeln!(fasta, ">chr1").unwrap();
        fasta.write_all(&bases).unwrap();
        writeln!(fasta).unwrap();
        fasta.flush().unwrap();
        (fasta, bases)
    }

    /// Asserts each mate's `MC` tag names the other mate's actual CIGAR.
    pub(in crate::commands::simulate) fn assert_mate_cigars_agree(
        r1: &noodles::sam::alignment::RecordBuf,
        r2: &noodles::sam::alignment::RecordBuf,
    ) {
        use noodles::sam::alignment::record::cigar::op::Kind;
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value;
        let cigar = |r: &noodles::sam::alignment::RecordBuf| {
            r.cigar().as_ref().iter().fold(String::new(), |mut cigar, op| {
                let code = match op.kind() {
                    Kind::Match => 'M',
                    Kind::SoftClip => 'S',
                    kind => panic!("unexpected CIGAR op {kind:?} in a simulated record"),
                };
                cigar.push_str(&op.len().to_string());
                cigar.push(code);
                cigar
            })
        };
        let mc = |r: &noodles::sam::alignment::RecordBuf| match r.data().get(&Tag::MATE_CIGAR) {
            Some(Value::String(s)) => String::from_utf8_lossy(s.as_ref()).into_owned(),
            other => panic!("MC must be a string, got {other:?}"),
        };
        assert_eq!(mc(r1), cigar(r2), "R1's MC must be R2's CIGAR");
        assert_eq!(mc(r2), cigar(r1), "R2's MC must be R1's CIGAR");
    }

    /// Asserts that a mapped simulated record's SEQ equals the reference over its aligned
    /// (`M`) bases, whatever its strand; soft-clipped bases (read-through padding) are skipped.
    /// Returns `(is_reverse, soft_clipped_bases)`.
    pub(in crate::commands::simulate) fn assert_seq_matches_reference(
        record: &noodles::sam::alignment::RecordBuf,
        reference: &[u8],
    ) -> (bool, usize) {
        use noodles::sam::alignment::record::cigar::op::Kind;
        let seq: Vec<u8> = record.sequence().as_ref().to_vec();
        let mut ref_pos = usize::from(record.alignment_start().expect("mapped record")) - 1;
        let mut read_pos = 0;
        let mut clipped = 0;
        let (mut aligned_seq, mut aligned_ref) = (Vec::new(), Vec::new());
        for op in record.cigar().as_ref() {
            match op.kind() {
                Kind::Match => {
                    aligned_seq.extend_from_slice(&seq[read_pos..read_pos + op.len()]);
                    aligned_ref.extend_from_slice(&reference[ref_pos..ref_pos + op.len()]);
                    read_pos += op.len();
                    ref_pos += op.len();
                }
                Kind::SoftClip => {
                    read_pos += op.len();
                    clipped += op.len();
                }
                kind => panic!("unexpected CIGAR op {kind:?} in a simulated record"),
            }
        }
        assert_eq!(read_pos, seq.len(), "CIGAR must cover SEQ");
        assert_eq!(
            String::from_utf8_lossy(&aligned_seq),
            String::from_utf8_lossy(&aligned_ref),
            "aligned SEQ must be in reference orientation (reverse = {})",
            record.flags().is_reverse_complemented()
        );
        (record.flags().is_reverse_complemented(), clipped)
    }
}

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

    /// Parses `MethylationArgs` alone from `args`.
    fn parse_methylation_args(args: &[&str]) -> MethylationArgs {
        #[derive(clap::Parser)]
        struct Cli {
            #[command(flatten)]
            methylation: MethylationArgs,
        }
        let mut argv = vec!["cli"];
        argv.extend_from_slice(args);
        <Cli as clap::Parser>::try_parse_from(argv).expect("valid args").methylation
    }

    /// Defaults match `holodeck simulate`/`methylate`: conversion 0.999, failure 0.01,
    /// hemimethylation 0.01.
    #[test]
    fn test_methylation_args_holodeck_defaults() {
        let config = parse_methylation_args(&["--methylation-mode", "em-seq"]).resolve(Some(7));
        assert!((config.conversion_rate - 0.999).abs() < f64::EPSILON);
        assert!((config.failure_rate - 0.01).abs() < f64::EPSILON);
        assert!((config.hemimethylation_rate - 0.01).abs() < f64::EPSILON);
        assert!((config.cpg_methylation_rate - 0.75).abs() < f64::EPSILON);
        assert_eq!(config.table_seed, 7, "the methylation state follows --seed");
    }

    /// `--methylation-conversion-rate` is the holodeck name; `--conversion-rate` stays an alias.
    #[rstest]
    #[case::holodeck_name("--methylation-conversion-rate")]
    #[case::legacy_alias("--conversion-rate")]
    fn test_methylation_args_conversion_rate_names(#[case] flag: &str) {
        let args = parse_methylation_args(&[flag, "0.9"]);
        assert!((args.resolve(None).conversion_rate - 0.9).abs() < f64::EPSILON);
    }

    #[rstest]
    #[case::failure("--methylation-failure-rate")]
    #[case::hemimethylation("--hemimethylation-rate")]
    fn test_methylation_args_validate_new_rates(#[case] flag: &str) {
        assert!(parse_methylation_args(&[flag, "0.5"]).validate().is_ok());
        assert!(parse_methylation_args(&[flag, "1.5"]).validate().is_err());
    }

    #[test]
    fn test_methylation_args_cli_defaults() {
        #[derive(clap::Parser)]
        struct Cli {
            #[command(flatten)]
            methylation: MethylationArgs,
        }
        let args = <Cli as clap::Parser>::try_parse_from(["simulate"]).unwrap().methylation;
        assert!(args.methylation_mode.is_none());
        assert!((args.cpg_methylation_rate - 0.75).abs() < f64::EPSILON);
        assert!((args.hemimethylation_rate - 0.01).abs() < f64::EPSILON);
        assert!((args.conversion_rate - 0.999).abs() < f64::EPSILON);
        assert!((args.failure_rate - 0.01).abs() < f64::EPSILON);
    }

    #[test]
    fn test_methylation_args_resolve_disabled() {
        let args = MethylationArgs {
            methylation_mode: None,
            cpg_methylation_rate: 0.75,
            conversion_rate: 0.98,
            hemimethylation_rate: 0.0,
            failure_rate: 0.0,
        };
        assert_eq!(args.resolve(None).mode, MethylationMode::Disabled);
    }

    #[test]
    fn test_methylation_args_resolve_emseq() {
        let args = MethylationArgs {
            methylation_mode: Some(MethylationModeArg::EmSeq),
            cpg_methylation_rate: 0.75,
            conversion_rate: 0.98,
            hemimethylation_rate: 0.0,
            failure_rate: 0.0,
        };
        let config = args.resolve(None);
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
            hemimethylation_rate: 0.0,
            failure_rate: 0.0,
        };
        assert_eq!(args.resolve(None).mode, MethylationMode::Taps);
    }

    #[test]
    fn test_methylation_args_validate_valid_rates() {
        for rate in [0.0, 0.5, 1.0] {
            let args = MethylationArgs {
                methylation_mode: None,
                cpg_methylation_rate: rate,
                conversion_rate: rate,
                hemimethylation_rate: 0.0,
                failure_rate: 0.0,
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
            hemimethylation_rate: 0.0,
            failure_rate: 0.0,
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
            hemimethylation_rate: 0.0,
            failure_rate: 0.0,
        };
        assert!(args.validate().is_err(), "conversion rate {rate} should be invalid");
    }

    // ========================================================================
    // convert_molecule_strand tests
    // ========================================================================

    /// A template of `n` repeats of `ACGT` (`CpG`s at 1-2, 5-6, ...), its own locus.
    fn cpg_template(n: usize) -> Vec<u8> {
        b"ACGT".repeat(n)
    }

    fn emseq(cpg_rate: f64, hemi_rate: f64, table_seed: u64) -> MethylationConfig {
        MethylationConfig {
            mode: MethylationMode::EmSeq,
            cpg_methylation_rate: cpg_rate,
            conversion_rate: 1.0,
            hemimethylation_rate: hemi_rate,
            failure_rate: 0.0,
            table_seed,
        }
    }

    /// Methylation is a property of the genome, not of a read: converting the same locus
    /// twice, with different chemistry draws, gives the same calls when conversion is perfect.
    #[test]
    fn test_convert_methylation_state_is_fixed_per_cpg() {
        let template = cpg_template(64);
        let config = emseq(0.5, 0.0, 7);
        let locus = TemplateLocus::standalone(&template);
        let first =
            convert_molecule_strand(&template, locus, true, &config, &mut create_rng(Some(1)));
        let second =
            convert_molecule_strand(&template, locus, true, &config, &mut create_rng(Some(2)));
        assert_eq!(first, second);
        assert_ne!(first, template, "some CpGs should be unmethylated (converted)");
        assert!(
            (0..64).any(|k| first[4 * k + 1] == b'C'),
            "some CpGs should be methylated (protected)"
        );
    }

    /// Without hemimethylation a `CpG` is methylated on both strands or neither: the top
    /// strand keeps its C exactly where the bottom strand keeps the G of the same `CpG`.
    #[test]
    fn test_convert_methylation_is_symmetric_without_hemimethylation() {
        let template = cpg_template(64);
        let config = emseq(0.5, 0.0, 11);
        let locus = TemplateLocus::standalone(&template);
        let top =
            convert_molecule_strand(&template, locus, true, &config, &mut create_rng(Some(3)));
        let bottom =
            convert_molecule_strand(&template, locus, false, &config, &mut create_rng(Some(4)));
        for k in 0..64 {
            assert_eq!(top[4 * k + 1] == b'C', bottom[4 * k + 2] == b'G', "CpG {k}");
        }
    }

    /// With every methylated `CpG` hemimethylated, exactly one strand of each keeps its base.
    #[test]
    fn test_convert_hemimethylation_drops_one_strand() {
        let template = cpg_template(64);
        let config = emseq(1.0, 1.0, 13);
        let locus = TemplateLocus::standalone(&template);
        let top =
            convert_molecule_strand(&template, locus, true, &config, &mut create_rng(Some(5)));
        let bottom =
            convert_molecule_strand(&template, locus, false, &config, &mut create_rng(Some(6)));
        for k in 0..64 {
            assert_ne!(top[4 * k + 1] == b'C', bottom[4 * k + 2] == b'G', "CpG {k}");
        }
    }

    /// Two templates not drawn from the reference get their own methylation states, even with
    /// `CpG`s at the same offsets.
    #[test]
    fn test_convert_standalone_templates_do_not_share_methylation() {
        let a_template = cpg_template(64);
        let mut b_template = a_template.clone();
        b_template[0] = b'T';
        let config = emseq(0.5, 0.0, 3);
        let convert = |template: &[u8]| {
            convert_molecule_strand(
                template,
                TemplateLocus::standalone(template),
                true,
                &config,
                &mut create_rng(Some(9)),
            )
        };
        let (a, b) = (convert(&a_template), convert(&b_template));
        // The `CpG` cytosines sit at offsets 1, 5, 9, ... in both templates.
        assert!((1..a.len()).step_by(4).any(|i| a[i] != b[i]));
    }

    /// The fixed state comes from the run's table seed.
    #[test]
    fn test_convert_methylation_state_depends_on_table_seed() {
        let template = cpg_template(64);
        let locus = TemplateLocus::standalone(&template);
        let a = convert_molecule_strand(
            &template,
            locus,
            true,
            &emseq(0.5, 0.0, 1),
            &mut create_rng(Some(9)),
        );
        let b = convert_molecule_strand(
            &template,
            locus,
            true,
            &emseq(0.5, 0.0, 2),
            &mut create_rng(Some(9)),
        );
        assert_ne!(a, b);
    }

    /// A failed molecule strand converts at `1 - conversion_rate`: with perfect chemistry,
    /// not at all.
    #[test]
    fn test_convert_failed_molecule_keeps_its_cytosines() {
        let template = b"ACCTACCTACCT".to_vec();
        let config = MethylationConfig { failure_rate: 1.0, ..emseq(0.0, 0.0, 1) };
        let converted = convert_molecule_strand(
            &template,
            TemplateLocus::standalone(&template),
            true,
            &config,
            &mut create_rng(Some(1)),
        );
        assert_eq!(converted, template);
    }

    /// Draws for different purposes at nearby `CpG`s are independent: with every methylated
    /// `CpG` hemimethylated, which strand carries the mark must not follow the neighbouring
    /// `CpG`'s methylation state (it would if the salt aliased one position onto another).
    #[test]
    fn test_hemimethylated_strand_is_independent_of_neighbour() {
        let contig = b"CG".repeat(4000);
        let locus = TemplateLocus { chrom_idx: Some(0), contig: &contig, start: 0 };
        let config = MethylationConfig {
            mode: fgumi_consensus::MethylationMode::EmSeq,
            cpg_methylation_rate: 0.5,
            conversion_rate: 1.0,
            hemimethylation_rate: 1.0,
            failure_rate: 0.0,
            table_seed: 7,
        };
        let methylated_any = |p: usize| {
            config.is_methylated(locus, p, true) || config.is_methylated(locus, p + 1, false)
        };
        let (mut same, mut total) = (0usize, 0usize);
        for p in (0..contig.len() - 4).step_by(4) {
            if methylated_any(p) {
                total += 1;
                if config.is_methylated(locus, p, true) == methylated_any(p + 2) {
                    same += 1;
                }
            }
        }
        assert!(total > 500, "enough methylated CpGs: {total}");
        #[expect(clippy::cast_precision_loss, reason = "small test counts")]
        let fraction = same as f64 / total as f64;
        assert!(fraction < 0.65, "strand choice follows the neighbour: {same}/{total}");
    }

    /// `CpG` context comes from the contig, not the template: a template ending in the C of
    /// a `CpG` whose G lies just past it still treats that C as `CpG` (methylated, protected).
    #[test]
    fn test_convert_cpg_context_uses_contig() {
        let contig = b"AACGAA";
        let template = &contig[..3]; // "AAC"
        let locus = TemplateLocus { chrom_idx: Some(0), contig, start: 0 };
        let converted = convert_molecule_strand(
            template,
            locus,
            true,
            &emseq(1.0, 0.0, 1),
            &mut create_rng(Some(1)),
        );
        assert_eq!(converted, b"AAC");
    }

    /// Disabled methylation returns the template and draws nothing from the RNG, so runs
    /// without methylation are unchanged.
    #[test]
    fn test_convert_disabled_is_identity_without_rng_draws() {
        let template = cpg_template(8);
        let mut rng = create_rng(Some(1));
        let converted = convert_molecule_strand(
            &template,
            TemplateLocus::standalone(&template),
            true,
            &MethylationConfig::default(),
            &mut rng,
        );
        assert_eq!(converted, template);
        assert_eq!(rng.random::<u64>(), create_rng(Some(1)).random::<u64>());
    }

    // ========================================================================
    // convert_molecule_strand chemistry tests
    // ========================================================================

    /// Helper: `template` converted as one molecule strand (its own locus), no failures or
    /// hemimethylation.
    fn convert(
        template: &[u8],
        is_top: bool,
        mode: MethylationMode,
        cpg_rate: f64,
        conv_rate: f64,
        seed: u64,
    ) -> Vec<u8> {
        // The seed drives both the fixed methylation state and the chemistry draws.
        let config = MethylationConfig {
            mode,
            cpg_methylation_rate: cpg_rate,
            conversion_rate: conv_rate,
            table_seed: seed,
            ..MethylationConfig::default()
        };
        let mut rng = create_rng(Some(seed));
        convert_molecule_strand(
            template,
            TemplateLocus::standalone(template),
            is_top,
            &config,
            &mut rng,
        )
    }

    #[test]
    fn test_emseq_cpg_all_methylated_no_conversion() {
        // EM-Seq: methylated CpG = protected, should NOT convert
        // cpg_methylation_rate=1.0 means all CpGs are methylated
        let ref_seq = b"ACGTACGT";
        let read = convert(ref_seq, true, MethylationMode::EmSeq, 1.0, 1.0, 42);
        // CpG Cs (at positions 1 and 5) should stay as C (methylated = protected in EM-Seq)
        assert_eq!(read[1], b'C', "CpG C should be protected when methylated");
        assert_eq!(read[5], b'C', "CpG C should be protected when methylated");
    }

    #[test]
    fn test_emseq_cpg_all_unmethylated_full_conversion() {
        // EM-Seq: unmethylated CpG = target, should convert C->T
        // cpg_methylation_rate=0.0 means all CpGs are unmethylated
        let ref_seq = b"ACGTACGT";
        let read = convert(ref_seq, true, MethylationMode::EmSeq, 0.0, 1.0, 42);
        // CpG Cs should convert to T
        assert_eq!(read[1], b'T', "unmethylated CpG C should convert to T");
        assert_eq!(read[5], b'T', "unmethylated CpG C should convert to T");
    }

    #[test]
    fn test_emseq_non_cpg_c_always_converts() {
        // Non-CpG Cs are unmethylated, always targets in EM-Seq
        // ref = "ACCTA" -> C at pos 1 (non-CpG, followed by C), C at pos 2 (non-CpG, followed by T)
        let ref_seq = b"ACCTA";
        let read = convert(ref_seq, true, MethylationMode::EmSeq, 0.75, 1.0, 42);
        assert_eq!(read[1], b'T', "non-CpG C should convert to T in EM-Seq");
        assert_eq!(read[2], b'T', "non-CpG C should convert to T in EM-Seq");
    }

    #[test]
    fn test_taps_cpg_all_methylated_full_conversion() {
        // TAPs: methylated CpG = target, should convert C->T
        let ref_seq = b"ACGTACGT";
        let read = convert(ref_seq, true, MethylationMode::Taps, 1.0, 1.0, 42);
        assert_eq!(read[1], b'T', "methylated CpG C should convert in TAPs");
        assert_eq!(read[5], b'T', "methylated CpG C should convert in TAPs");
    }

    #[test]
    fn test_taps_cpg_all_unmethylated_no_conversion() {
        // TAPs: unmethylated CpG = not a target, should stay as C
        let ref_seq = b"ACGTACGT";
        let read = convert(ref_seq, true, MethylationMode::Taps, 0.0, 1.0, 42);
        assert_eq!(read[1], b'C', "unmethylated CpG C should not convert in TAPs");
        assert_eq!(read[5], b'C', "unmethylated CpG C should not convert in TAPs");
    }

    #[test]
    fn test_taps_non_cpg_c_never_converts() {
        // Non-CpG Cs are unmethylated, never targets in TAPs
        let ref_seq = b"ACCTA";
        let read = convert(ref_seq, true, MethylationMode::Taps, 0.75, 1.0, 42);
        assert_eq!(read[1], b'C', "non-CpG C should not convert in TAPs");
        assert_eq!(read[2], b'C', "non-CpG C should not convert in TAPs");
    }

    #[test]
    fn test_bottom_strand_emseq_converts_g_to_a() {
        // Bottom strand: G at CpG context = unmethylated target in EM-Seq
        // ref = "ACGT" -> G at pos 2, preceded by C -> CpG context
        let ref_seq = b"ACGT";
        let read = convert(ref_seq, false, MethylationMode::EmSeq, 0.0, 1.0, 42);
        assert_eq!(read[2], b'A', "bottom strand unmethylated CpG G should convert to A");
    }

    #[test]
    fn test_bottom_strand_taps_converts_g_to_a() {
        // Bottom strand: G at CpG context = methylated target in TAPs
        let ref_seq = b"ACGT";
        let read = convert(ref_seq, false, MethylationMode::Taps, 1.0, 1.0, 42);
        assert_eq!(read[2], b'A', "bottom strand methylated CpG G should convert to A in TAPs");
    }

    #[test]
    fn test_non_target_bases_unchanged_top_strand() {
        let ref_seq = b"AGTAGT";
        let read = convert(ref_seq, true, MethylationMode::EmSeq, 0.0, 1.0, 42);
        // No Cs in this sequence, nothing should change
        assert_eq!(read, b"AGTAGT");
    }

    #[test]
    fn test_non_target_bases_unchanged_bottom_strand() {
        let ref_seq = b"ACTACT";
        let read = convert(ref_seq, false, MethylationMode::EmSeq, 0.0, 1.0, 42);
        // No Gs in this sequence, nothing should change on bottom strand
        assert_eq!(read, b"ACTACT");
    }

    #[test]
    fn test_disabled_mode_no_conversion() {
        let ref_seq = b"ACGTACGT";
        let read = convert(ref_seq, true, MethylationMode::Disabled, 0.0, 1.0, 42);
        assert_eq!(read, ref_seq, "Disabled mode should not modify any bases");
    }

    #[test]
    fn test_empty_sequence() {
        let ref_seq = b"";
        let read = convert(ref_seq, true, MethylationMode::EmSeq, 0.75, 0.98, 42);
        assert!(read.is_empty());
    }

    #[test]
    fn test_ref_offset_nonzero() {
        // The template starts at offset 2 of the contig; its C at contig offset 2 is a CpG C.
        let contig = b"AACGTAA";
        let config = MethylationConfig {
            mode: MethylationMode::EmSeq,
            cpg_methylation_rate: 0.0,
            conversion_rate: 1.0,
            ..MethylationConfig::default()
        };
        let locus = TemplateLocus { chrom_idx: Some(0), contig, start: 2 };
        let read =
            convert_molecule_strand(&contig[2..5], locus, true, &config, &mut create_rng(Some(42)));
        assert_eq!(read, b"TGT", "unmethylated CpG C at contig offset 2 should convert");
    }

    #[test]
    fn test_conversion_rate_zero_no_conversion() {
        // Even with unmethylated non-CpG C, conversion_rate=0 means no conversion
        let ref_seq = b"ACCTA";
        let read = convert(ref_seq, true, MethylationMode::EmSeq, 0.0, 0.0, 42);
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
            let read = convert(ref_seq, true, MethylationMode::EmSeq, 0.5, 1.0, seed);
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
            let read = convert(ref_seq, true, MethylationMode::EmSeq, 0.75, 0.5, seed);
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
        let read = convert(ref_seq, true, MethylationMode::EmSeq, 0.0, 0.0, 42);
        assert_eq!(read, ref_seq, "conversion_rate=0 should leave all bases unchanged");
    }

    #[test]
    fn test_conversion_rate_one_converts_all_targets() {
        // EM-Seq, cpg_methylation_rate=0 means all CpGs unmethylated -> all Cs are targets
        let ref_seq = b"CACACACACACACACAC"; // all non-CpG Cs
        let read = convert(ref_seq, true, MethylationMode::EmSeq, 0.0, 1.0, 42);
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
        let read = convert(ref_seq, true, MethylationMode::Disabled, 0.75, 1.0, 42);
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

    // ========================================================================
    // RecordSink tests
    // ========================================================================

    /// Build a minimal named record for the sink tests.
    fn sink_record(name: &str) -> RawRecord {
        let mut builder = fgumi_raw_bam::SamBuilder::new();
        builder
            .read_name(name.as_bytes())
            .flags(fgumi_raw_bam::flags::UNMAPPED)
            .sequence(b"ACGT")
            .qualities(&[30; 4]);
        builder.build()
    }

    /// The names the sink actually delivered, in order, plus any error it sent.
    fn drain(receiver: &Receiver<RecordBatch>) -> (Vec<String>, Vec<String>) {
        let mut names = Vec::new();
        let mut errors = Vec::new();
        for batch in receiver {
            match batch {
                Ok(records) => names.extend(
                    records
                        .iter()
                        .map(|r| String::from_utf8_lossy(fgumi_raw_bam::read_name(r)).into_owned()),
                ),
                Err(e) => errors.push(format!("{e:#}")),
            }
        }
        (names, errors)
    }

    /// `finish` must deliver the last partial batch, by record identity.
    #[test]
    fn test_record_sink_finish_delivers_the_partial_batch() {
        let (mut sink, receiver) = RecordSink::new();
        assert!(sink.send(sink_record("readA")), "the receiver is still open");
        assert!(sink.send(sink_record("readB")), "the receiver is still open");
        assert!(sink.finish(), "finish must succeed while the receiver is open");

        let (names, errors) = drain(&receiver);
        assert_eq!(names, vec!["readA".to_string(), "readB".to_string()]);
        assert!(errors.is_empty(), "a clean finish must not report an error: {errors:?}");
    }

    /// Dropping a sink with records still buffered must abort the sort.
    ///
    /// Dropping the sender is what signals end-of-stream, so a sink abandoned
    /// with a partial batch — a `?` on the generator's path, or a panic —
    /// used to look exactly like a clean end of input: the sort finished and
    /// wrote a complete, short BAM. Reporting the loss as the stream's terminal
    /// error is what makes the sort fail instead.
    #[test]
    fn test_record_sink_drop_reports_the_unflushed_batch() {
        let (mut sink, receiver) = RecordSink::new();
        assert!(sink.send(sink_record("readA")), "the receiver is still open");
        drop(sink);

        let (names, errors) = drain(&receiver);
        assert!(
            names.is_empty(),
            "an abandoned partial batch must not be delivered as if it were complete: {names:?}"
        );
        assert_eq!(errors.len(), 1, "the stream must end in exactly one error: {errors:?}");
        assert!(
            errors[0].contains("1 generated record"),
            "the error must name how many records were left unflushed, got: {}",
            errors[0]
        );
    }

    /// A finished sink has nothing left, so its drop must stay silent.
    #[test]
    fn test_record_sink_drop_after_finish_reports_nothing() {
        let (mut sink, receiver) = RecordSink::new();
        assert!(sink.send(sink_record("readA")), "the receiver is still open");
        assert!(sink.finish(), "finish must succeed while the receiver is open");

        let (names, errors) = drain(&receiver);
        assert_eq!(names, vec!["readA".to_string()]);
        assert!(errors.is_empty(), "a finished sink must not also report a loss: {errors:?}");
    }

    /// `fail` abandons its buffer deliberately, and must report only its own error.
    #[test]
    fn test_record_sink_fail_reports_only_the_generator_error() {
        let (mut sink, receiver) = RecordSink::new();
        assert!(sink.send(sink_record("readA")), "the receiver is still open");
        sink.fail(anyhow!("generator exploded"));

        let (names, errors) = drain(&receiver);
        assert!(names.is_empty(), "an aborted stream delivers no records: {names:?}");
        assert_eq!(errors.len(), 1, "the abort must not be followed by a drop error: {errors:?}");
        assert!(errors[0].contains("generator exploded"), "got: {}", errors[0]);
    }
}

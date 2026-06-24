//! Generate consensus BAM with tags for filter.

use crate::commands::command::Command;
use crate::commands::common::{CompressionOptions, parse_bool};
use crate::commands::simulate::common::{MethylationArgs, ReferenceGenome, StrandBiasArgs};
use crate::commands::simulate::region_to_bin;
use crate::dna::reverse_complement;
use crate::sam::SamTag;
use crate::simulate::{StrandBiasModel, create_rng};
use anyhow::{Context, Result};
use clap::Parser;
use crossbeam_channel::bounded;
use fgumi_bam_io::ProgressTracker;
use fgumi_bam_io::create_raw_bam_writer;
use fgumi_consensus::MethylationMode;
use fgumi_consensus::methylation::{
    MethylationAnnotation, MethylationEvidence, build_mm_ml_tags, is_cpg_context,
};
use fgumi_raw_bam::{RawRecord, SamBuilder, flags as raw_flags};
use log::info;
use noodles::sam::header::Header;
use rand::{Rng, RngExt};
use rand_distr::{Distribution, LogNormal, Normal};
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::sync::Arc;
use std::thread;

/// Generate mapped BAM with consensus tags for `fgumi filter`.
#[derive(Parser, Debug)]
#[command(
    name = "consensus-reads",
    about = "Generate consensus BAM with tags for filter",
    long_about = r#"
Generate synthetic consensus reads with proper consensus tags.

The output is a mapped BAM suitable for input to `fgumi filter`.
Reads contain consensus tags (cD, cM, cE, cd, ce) for filtering.
"#
)]
pub struct ConsensusReads {
    /// Output BAM file (mapped)
    #[arg(short = 'o', long = "output", required = true)]
    pub output: PathBuf,

    /// Output truth TSV file for validation
    #[arg(long = "truth")]
    pub truth_output: Option<PathBuf>,

    /// Number of consensus read pairs to generate
    #[arg(short = 'n', long = "num-reads", default_value = "1000")]
    pub num_reads: usize,

    /// Read length in bases
    #[arg(short = 'l', long = "read-length", default_value = "150")]
    pub read_length: usize,

    /// Random seed for reproducibility
    #[arg(long = "seed")]
    pub seed: Option<u64>,

    /// Number of writer threads
    #[arg(short = 't', long = "threads", default_value = "1")]
    pub threads: usize,

    /// Compression options for output BAM.
    #[command(flatten)]
    pub compression: CompressionOptions,

    /// Minimum consensus depth (cM tag). Must satisfy 0 <= min <= max <= 32767
    /// (per-base depth arrays are emitted as i16).
    #[arg(long = "min-depth", default_value = "1")]
    pub min_depth: i32,

    /// Maximum consensus depth (cD tag). Must satisfy 0 <= min <= max <= 32767
    /// (per-base depth arrays are emitted as i16).
    #[arg(long = "max-depth", default_value = "10")]
    pub max_depth: i32,

    /// Mean depth for sampling
    #[arg(long = "depth-mean", default_value = "5.0")]
    pub depth_mean: f64,

    /// Depth standard deviation
    #[arg(long = "depth-stddev", default_value = "2.0")]
    pub depth_stddev: f64,

    /// Mean per-base error rate used to derive the simulated error *count*
    /// (`ce = round(read_length * error_rate)`). NOTE: this is the per-base
    /// error / read-length input rate, not the value of the emitted read-level
    /// `cE` tag. The `cE`/`aE`/`bE` tags are written as the real-caller rate
    /// `sum(per-base errors) / sum(per-base depths)` (matching
    /// `fgumi simplex`/`duplex` output), so the emitted `cE` is roughly this
    /// requested rate divided by the mean simulated depth.
    #[arg(long = "error-rate-mean", default_value = "0.01")]
    pub error_rate_mean: f64,

    /// Error rate standard deviation
    #[arg(long = "error-rate-stddev", default_value = "0.005")]
    pub error_rate_stddev: f64,

    /// Generate duplex consensus tags (aD, bD, aM, bM, aE, bE)
    #[arg(long = "duplex", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool)]
    pub duplex: bool,

    /// Base quality for consensus reads
    #[arg(long = "consensus-quality", default_value = "40")]
    pub consensus_quality: u8,

    #[command(flatten)]
    pub strand_bias: StrandBiasArgs,

    #[command(flatten)]
    pub methylation: MethylationArgs,

    /// Mean depth for methylation count sampling (cu + ct per position).
    /// Stddev is set to half the mean.
    #[arg(long = "methylation-depth-mean", default_value = "5.0")]
    pub methylation_depth_mean: f64,

    /// Reference FASTA file for sampling template sequences and building BAM headers.
    #[arg(short = 'r', long = "reference", required = true)]
    pub reference: PathBuf,
}

/// Header line for the truth TSV emitted by `--truth`.
///
/// The depth columns (`cD`/`cM`/`aD`/`bD`/`aM`/`bM`) carry the same integer
/// counts as the matching BAM tags. The error columns are suffixed `_count`
/// (`cE_count`/`aE_count`/`bE_count`) because the BAM's read-level `cE`/`aE`/`bE`
/// tags are FLOAT error *rates* (`sum(per-base errors) / sum(per-base depths)`,
/// the real-caller formula — see `build_consensus_record`), whereas the truth
/// records the simulated per-read error *count*. Naming them `_count` keeps the
/// TSV from silently disagreeing with the BAM it validates.
const TRUTH_TSV_HEADER: &str =
    "read_name\tchrom\tpos\tstrand\tcD\tcM\tcE_count\taD\tbD\taM\tbM\taE_count\tbE_count";

/// A generated consensus read pair ready for output.
struct ConsensusReadPair {
    read_name: String,
    r1_record: RawRecord,
    r2_record: RawRecord,
    /// Truth data: per-read integer counts
    /// `(cD, cM, cE_count, aD, bD, aM, bM, aE_count, bE_count)`. The error
    /// entries are counts, not the FLOAT rates emitted as the BAM `cE`/`aE`/`bE`
    /// tags — see [`TRUTH_TSV_HEADER`].
    truth: (i32, i32, i32, i32, i32, i32, i32, i32, i32),
    /// Chromosome name for truth output.
    chrom_name: String,
    /// 0-based local position within the chromosome.
    local_pos: usize,
    /// Whether the molecule was sampled from the top strand.
    is_top_strand: bool,
}

/// Parameters needed for parallel consensus generation.
struct GenerationParams {
    read_length: usize,
    min_depth: i32,
    max_depth: i32,
    depth_mean: f64,
    depth_stddev: f64,
    error_rate_mean: f64,
    error_rate_stddev: f64,
    duplex: bool,
    consensus_quality: u8,
    /// Methylation mode (Disabled if not set).
    methylation_mode: MethylationMode,
    /// Distribution for methylation count sampling.
    methylation_depth_dist: LogNormal<f64>,
    /// `CpG` methylation rate.
    cpg_methylation_rate: f64,
    /// Enzymatic conversion rate for target cytosines.
    conversion_rate: f64,
    /// Loaded reference genome for sampling template sequences.
    ref_genome: Arc<ReferenceGenome>,
}

/// Channel capacity for buffering read pairs between producer and writer threads.
const CHANNEL_CAPACITY: usize = 1_000;

impl Command for ConsensusReads {
    fn execute(&self, command_line: &str) -> Result<()> {
        let methylation = self.methylation.resolve();
        self.methylation.validate()?;
        if methylation.mode.is_enabled()
            && (!self.methylation_depth_mean.is_finite() || self.methylation_depth_mean <= 0.0)
        {
            anyhow::bail!(
                "--methylation-depth-mean must be finite and positive when methylation is enabled, got {}",
                self.methylation_depth_mean
            );
        }
        validate_depth_bounds(self.min_depth, self.max_depth)?;

        info!("Generating consensus reads");
        info!("  Output: {}", self.output.display());
        info!("  Num reads: {}", self.num_reads);
        info!("  Read length: {}", self.read_length);
        info!("  Duplex: {}", self.duplex);
        info!("  Depth range: {}-{}", self.min_depth, self.max_depth);
        info!("  Threads: {}", self.threads);
        if methylation.mode.is_enabled() {
            info!("  Methylation mode: {:?}", methylation.mode);
            info!("  Methylation depth mean: {}", self.methylation_depth_mean);
            info!("  CpG methylation rate: {}", methylation.cpg_methylation_rate);
        }

        // Load reference genome
        let ref_genome = Arc::new(ReferenceGenome::load(&self.reference)?);

        // Validate that the reference has at least one contig >= read_length
        if ref_genome.max_contig_length() < self.read_length {
            anyhow::bail!(
                "No reference contig is >= read length ({} bp). \
                 The longest contig is {} bp. Use a larger reference or shorter --read-length.",
                self.read_length,
                ref_genome.max_contig_length(),
            );
        }

        // Build header from reference contigs
        let ref_header = ref_genome.build_bam_header();
        let mut header_builder = Header::builder();
        for (name, map) in ref_header.reference_sequences() {
            header_builder = header_builder.add_reference_sequence(name.clone(), map.clone());
        }
        header_builder = crate::commands::common::add_pg_to_builder(header_builder, command_line)?;
        let header = header_builder.build();

        // Set up shared parameters
        let params = Arc::new(GenerationParams {
            read_length: self.read_length,
            min_depth: self.min_depth,
            max_depth: self.max_depth,
            depth_mean: self.depth_mean,
            depth_stddev: self.depth_stddev,
            error_rate_mean: self.error_rate_mean,
            error_rate_stddev: self.error_rate_stddev,
            duplex: self.duplex,
            consensus_quality: self.consensus_quality,
            methylation_mode: methylation.mode,
            methylation_depth_dist: if methylation.mode.is_enabled() {
                create_depth_distribution(
                    self.methylation_depth_mean,
                    self.methylation_depth_mean / 2.0,
                )
            } else {
                // Unused when methylation is disabled; use safe defaults
                create_depth_distribution(5.0, 2.5)
            },
            cpg_methylation_rate: methylation.cpg_methylation_rate,
            conversion_rate: methylation.conversion_rate,
            ref_genome: Arc::clone(&ref_genome),
        });

        let strand_bias_model = Arc::new(self.strand_bias.to_strand_bias_model());

        // Generate seeds for reproducibility
        let mut seed_rng = create_rng(self.seed);
        let read_seeds: Vec<u64> = (0..self.num_reads).map(|_| seed_rng.random()).collect();

        // Create bounded channel for streaming read pairs to writer
        let (sender, receiver) = bounded::<ConsensusReadPair>(CHANNEL_CAPACITY);

        // Clone paths for writer thread
        let output_path = self.output.clone();
        let truth_path = self.truth_output.clone();
        let compression_level = self.compression.compression_level;
        let writer_threads = self.threads;
        let header_clone = header.clone();

        // Spawn writer thread with multi-threaded BGZF compression
        let writer_handle = thread::spawn(move || -> Result<u64> {
            let mut writer = create_raw_bam_writer(
                &output_path,
                &header_clone,
                writer_threads,
                compression_level,
            )?;

            // Create truth file if requested
            let mut truth_writer = if let Some(ref truth_path) = truth_path {
                let truth_file = File::create(truth_path)
                    .with_context(|| format!("Failed to create {}", truth_path.display()))?;
                let mut w = BufWriter::new(truth_file);
                writeln!(w, "{TRUTH_TSV_HEADER}")?;
                Some(w)
            } else {
                None
            };

            let mut read_count = 0u64;
            let progress = ProgressTracker::new("Generated consensus pairs").with_interval(100_000);

            // Receive and write read pairs as they arrive
            for pair in receiver {
                read_count += 1;
                progress.log_if_needed(1);

                writer.write_raw_record(pair.r1_record.as_ref())?;
                writer.write_raw_record(pair.r2_record.as_ref())?;

                // Write truth
                if let Some(ref mut tw) = truth_writer {
                    let (cd, cm, ce, ad, bd, am, bm, ae, be) = pair.truth;
                    let strand_char = if pair.is_top_strand { '+' } else { '-' };
                    writeln!(
                        tw,
                        "{}\t{}\t{}\t{strand_char}\t{cd}\t{cm}\t{ce}\t{ad}\t{bd}\t{am}\t{bm}\t{ae}\t{be}",
                        pair.read_name, pair.chrom_name, pair.local_pos
                    )?;
                }
            }

            progress.log_final();

            if let Some(ref mut tw) = truth_writer {
                tw.flush()?;
            }

            writer.finish()?;

            Ok(read_count)
        });

        // Configure thread pool for generation
        let gen_threads = if self.threads <= 1 { 1 } else { self.threads.max(2) };
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(gen_threads)
            .build()
            .with_context(|| "Failed to create thread pool")?;

        // Generate reads in parallel and stream to writer
        let generation_result: Result<(), crossbeam_channel::SendError<ConsensusReadPair>> = pool
            .install(|| {
                read_seeds.into_par_iter().enumerate().try_for_each(|(read_idx, seed)| {
                    let pair = generate_consensus_pair(read_idx, seed, &params, &strand_bias_model);
                    sender.send(pair)
                })
            });

        // Drop sender to signal writer thread that we're done
        drop(sender);

        // Check generation result
        if let Err(e) = generation_result {
            return Err(anyhow::anyhow!("Failed to send record to writer: {e}"));
        }

        // Wait for writer thread to finish
        let read_count =
            writer_handle.join().map_err(|_| anyhow::anyhow!("Writer thread panicked"))??;

        info!("Generated {read_count} consensus read pairs");
        info!("Done");

        Ok(())
    }
}

/// Generate a single consensus read pair.
fn generate_consensus_pair(
    read_idx: usize,
    seed: u64,
    params: &GenerationParams,
    strand_bias_model: &StrandBiasModel,
) -> ConsensusReadPair {
    let mut rng = create_rng(Some(seed));

    let read_name = format!("consensus_{read_idx:08}");

    // Create distributions for this read
    let depth_dist = create_depth_distribution(params.depth_mean, params.depth_stddev);
    let error_dist = Normal::new(params.error_rate_mean, params.error_rate_stddev)
        .expect("Invalid error distribution parameters");

    // Sample sequence from reference
    let (chrom_idx, local_pos, seq) = params
        .ref_genome
        .sample_sequence(params.read_length, &mut rng)
        .expect("Failed to sample sequence from reference");
    let quals = vec![params.consensus_quality; params.read_length];

    // Generate consensus depth (cD)
    let cd = sample_depth(&depth_dist, params.min_depth, params.max_depth, &mut rng);

    // Generate minimum depth (cM) - at most cD
    let cm = sample_depth(&depth_dist, params.min_depth, cd, &mut rng).min(cd);

    // Generate the error *count* from the sampled per-base error rate. NOTE:
    // `ce` is `read_length * error_rate` (a count), which is then spread
    // one-per-base in the per-base arrays. The emitted read-level cE/aE/bE tag
    // is NOT this rate but `sum(per-base errors) / sum(per-base depths)` — the
    // real `vanilla_caller` formula (see `build_consensus_record` and the
    // `--error-rate-mean` doc).
    let error_rate = error_dist.sample(&mut rng).clamp(0.0, 1.0);
    let ce = (params.read_length as f64 * error_rate).round() as i32;

    // For duplex mode, generate strand-specific depths
    let (ad, bd, am, bm, ae, be) = if params.duplex {
        // Split total depth between A and B strands
        let a_frac = strand_bias_model.sample_a_fraction(&mut rng);

        let ad = ((cd as f64) * a_frac).round() as i32;
        let bd = cd - ad;

        // Min depths for each strand. The real duplex/codec callers (and fgbio)
        // derive cM = min(aD_i + bD_i); because both per-strand depth arrays ramp
        // monotonically to their last base, that minimum equals aM + bM. Split cM
        // by the same strand fraction used for cD so aM + bM == cM exactly, while
        // keeping aM <= aD and bM <= bD. The clamp bounds are well-ordered because
        // cM <= cD == aD + bD, so (cM - bD) <= aD and (cM - bD) <= cM.
        let am = (((cm as f64) * a_frac).round() as i32).clamp((cm - bd).max(0), ad.min(cm));
        let bm = cm - am;

        // Errors distributed proportionally
        let ae = ((ce as f64) * a_frac).round() as i32;
        let be = ce - ae;

        (ad, bd, am, bm, ae, be)
    } else {
        (0, 0, 0, 0, 0, 0)
    };

    // A single-base read cannot carry distinct max/min depths: its only base is
    // emitted at the minimum depth (see `fill_per_base`'s `read_len == 1` arm),
    // so collapse each summary maximum (cD/aD/bD) down to its minimum
    // (cM/aM/bM) here — upstream of both the record tags and the truth tuple —
    // so they stay self-consistent and cE/aE/bE are computed over the right depth.
    let (cd, ad, bd) = if params.read_length == 1 { (cm, am, bm) } else { (cd, ad, bd) };

    // Generate methylation annotation if enabled
    let methylation = if params.methylation_mode.is_enabled() {
        Some(generate_methylation_annotation(
            &seq,
            params.methylation_mode,
            params.cpg_methylation_rate,
            params.conversion_rate,
            &params.methylation_depth_dist,
            params.duplex,
            &mut rng,
        ))
    } else {
        None
    };

    // Coin flip for strand orientation
    let is_top_strand: bool = rng.random();
    let r1_is_reverse = !is_top_strand;

    // Compute methylation per read: the forward read gets the original annotation,
    // the reverse read gets the reversed annotation. Which read is forward depends
    // on the strand coin flip.
    let reverse_methylation = methylation.as_ref().map(|m| m.reverse());
    let (r1_methylation, r2_methylation) = if r1_is_reverse {
        // R1F2: R1 is reverse, R2 is forward
        (reverse_methylation.as_ref(), methylation.as_ref())
    } else {
        // F1R2: R1 is forward, R2 is reverse
        (methylation.as_ref(), reverse_methylation.as_ref())
    };

    // Build R1 record
    let r1_seq = if is_top_strand { seq.clone() } else { reverse_complement(&seq) };
    let r1_record = build_consensus_record(
        &format!("{read_name}/1"),
        &r1_seq,
        &quals,
        true, // is_first
        r1_is_reverse,
        chrom_idx,
        local_pos,
        cd,
        cm,
        ce,
        if params.duplex { Some((ad, bd, am, bm, ae, be)) } else { None },
        r1_methylation,
        params.methylation_mode,
    );

    // Build R2 record (opposite strand orientation)
    let r2_seq = if is_top_strand { reverse_complement(&seq) } else { seq.clone() };
    let r2_record = build_consensus_record(
        &format!("{read_name}/2"),
        &r2_seq,
        &quals,
        false, // is_first
        !r1_is_reverse,
        chrom_idx,
        local_pos,
        cd,
        cm,
        ce,
        if params.duplex { Some((ad, bd, am, bm, ae, be)) } else { None },
        r2_methylation,
        params.methylation_mode,
    );

    let chrom_name = params.ref_genome.name(chrom_idx).to_string();

    ConsensusReadPair {
        read_name,
        r1_record,
        r2_record,
        truth: (cd, cm, ce, ad, bd, am, bm, ae, be),
        chrom_name,
        local_pos,
        is_top_strand,
    }
}

/// Generated methylation data for a consensus read.
struct MethylationData {
    /// Combined annotation (both strands merged for simplex, or total for duplex).
    annotation: MethylationAnnotation,
    /// AB strand annotation (duplex only).
    ab_annotation: Option<MethylationAnnotation>,
    /// BA strand annotation (duplex only).
    ba_annotation: Option<MethylationAnnotation>,
}

impl MethylationData {
    /// Returns a copy with all annotations reversed.
    ///
    /// Used for R2 records where the sequence is reverse-complemented:
    /// the per-position methylation evidence must be reversed to match.
    fn reverse(&self) -> Self {
        Self {
            annotation: self.annotation.reverse(),
            ab_annotation: self.ab_annotation.as_ref().map(MethylationAnnotation::reverse),
            ba_annotation: self.ba_annotation.as_ref().map(MethylationAnnotation::reverse),
        }
    }
}

/// Generate methylation annotation for a consensus sequence.
fn generate_methylation_annotation(
    seq: &[u8],
    mode: MethylationMode,
    cpg_methylation_rate: f64,
    conversion_rate: f64,
    methylation_depth_dist: &LogNormal<f64>,
    duplex: bool,
    rng: &mut impl Rng,
) -> MethylationData {
    let mut evidence = Vec::with_capacity(seq.len());
    let mut ab_evidence = if duplex { Some(Vec::with_capacity(seq.len())) } else { None };
    let mut ba_evidence = if duplex { Some(Vec::with_capacity(seq.len())) } else { None };

    for i in 0..seq.len() {
        let base = seq[i].to_ascii_uppercase();
        let is_ref_c = base == b'C';

        if !is_ref_c {
            let zero_ev =
                MethylationEvidence { is_ref_c: false, unconverted_count: 0, converted_count: 0 };
            evidence.push(zero_ev.clone());
            if let Some(ref mut ab) = ab_evidence {
                ab.push(zero_ev.clone());
            }
            if let Some(ref mut ba) = ba_evidence {
                ba.push(zero_ev);
            }
            continue;
        }

        // This is a C position - determine if CpG
        let cpg = is_cpg_context(seq, i, true);

        // Determine methylation probability
        let methylated_prob = if cpg { cpg_methylation_rate } else { 0.0 };

        // Sample total depth for this position
        let total_depth = sample_depth(methylation_depth_dist, 1, 100, rng) as u32;

        // Sample how many reads show the "unconverted" (methylated in EM-Seq) signal
        let methylated = rng.random::<f64>() < methylated_prob;

        // In EM-Seq: unconverted = methylated (protected), converted = unmethylated
        // In TAPs: unconverted = unmethylated, converted = methylated
        //
        // The conversion_rate controls what fraction of targeted bases are converted.
        // For the "mostly converted" branch: converted ~= total * conversion_rate
        // For the "mostly unconverted" branch: unconverted ~= total * conversion_rate
        // (i.e., conversion_rate reflects the efficiency of the chemistry)
        let (unconverted, converted) = match mode {
            MethylationMode::EmSeq => {
                if cpg && methylated {
                    // Methylated CpG in EM-Seq: mostly unconverted (protected)
                    let converted_count =
                        (total_depth as f64 * (1.0 - conversion_rate)).round() as u32;
                    (total_depth - converted_count, converted_count)
                } else {
                    // Unmethylated (non-CpG or unmethylated CpG): mostly converted
                    let unconverted_count =
                        (total_depth as f64 * (1.0 - conversion_rate)).round() as u32;
                    (unconverted_count, total_depth - unconverted_count)
                }
            }
            MethylationMode::Taps => {
                if cpg && methylated {
                    // Methylated CpG in TAPs: mostly converted (target)
                    let unconverted_count =
                        (total_depth as f64 * (1.0 - conversion_rate)).round() as u32;
                    (unconverted_count, total_depth - unconverted_count)
                } else {
                    // Unmethylated: mostly unconverted (not a target). The minority
                    // leakage (at most `1 - conversion_rate` fraction) is the number
                    // of spuriously converted reads.
                    let leaked_converted =
                        (total_depth as f64 * (1.0 - conversion_rate)).round() as u32;
                    (total_depth - leaked_converted, leaked_converted)
                }
            }
            MethylationMode::Disabled => (0, 0),
        };

        let ev = MethylationEvidence {
            is_ref_c: true,
            unconverted_count: unconverted,
            converted_count: converted,
        };

        if duplex {
            // Split counts between AB and BA strands
            let ab_frac = rng.random_range(30..=70) as f64 / 100.0;
            let ab_unconverted = (unconverted as f64 * ab_frac).round() as u32;
            let ab_converted = (converted as f64 * ab_frac).round() as u32;
            let ba_unconverted = unconverted - ab_unconverted;
            let ba_converted = converted - ab_converted;

            if let Some(ref mut ab) = ab_evidence {
                ab.push(MethylationEvidence {
                    is_ref_c: true,
                    unconverted_count: ab_unconverted,
                    converted_count: ab_converted,
                });
            }
            if let Some(ref mut ba) = ba_evidence {
                ba.push(MethylationEvidence {
                    is_ref_c: true,
                    unconverted_count: ba_unconverted,
                    converted_count: ba_converted,
                });
            }
        }

        evidence.push(ev);
    }

    MethylationData {
        annotation: MethylationAnnotation { evidence },
        ab_annotation: ab_evidence.map(|e| MethylationAnnotation { evidence: e }),
        ba_annotation: ba_evidence.map(|e| MethylationAnnotation { evidence: e }),
    }
}

fn create_depth_distribution(mean: f64, stddev: f64) -> LogNormal<f64> {
    // Convert mean/stddev to log-normal parameters
    let variance = stddev.powi(2);
    let mean_sq = mean.powi(2);
    let sigma_sq = (1.0 + variance / mean_sq).ln();
    let sigma = sigma_sq.sqrt();
    let mu = mean.ln() - sigma_sq / 2.0;

    LogNormal::new(mu, sigma).expect("Invalid log-normal parameters")
}

fn sample_depth(dist: &LogNormal<f64>, min: i32, max: i32, rng: &mut impl Rng) -> i32 {
    let sample = dist.sample(rng).round() as i32;
    sample.clamp(min, max)
}

/// Validates the consensus depth bounds before generation. The per-base `cd`/`ce`
/// (`CD_BASES`/`CE_BASES`) and strand (`AD_BASES`/`BD_BASES`) arrays are emitted
/// as `i16`; without this guard `fill_per_base`'s saturating cast would silently
/// make those arrays diverge from the `cD`/`cM` summary tags (and the truth TSV)
/// once a depth exceeds `i16::MAX`, and a negative `--min-depth` could emit
/// negative depths. Reject those inputs up front instead.
fn validate_depth_bounds(min_depth: i32, max_depth: i32) -> Result<()> {
    if min_depth < 0 || max_depth < min_depth || max_depth > i32::from(i16::MAX) {
        anyhow::bail!(
            "--min-depth/--max-depth must satisfy 0 <= min <= max <= {}, got {}..={}",
            i16::MAX,
            min_depth,
            max_depth,
        );
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn build_consensus_record(
    name: &str,
    seq: &[u8],
    quals: &[u8],
    is_first: bool,
    is_reverse: bool,
    chrom_idx: usize,
    local_pos: usize,
    cd: i32,
    cm: i32,
    ce: i32,
    duplex_tags: Option<(i32, i32, i32, i32, i32, i32)>,
    methylation: Option<&MethylationData>,
    methylation_mode: MethylationMode,
) -> RawRecord {
    let is_top_strand = !is_reverse;

    // Build flags: PAIRED + PROPER_PAIR + (FIRST_SEGMENT or LAST_SEGMENT) + reverse-strand bits.
    // Consensus records are simulated proper paired alignments; downstream tools commonly
    // filter on the 0x2 bit, so set it explicitly.
    let segment_flag = if is_first { raw_flags::FIRST_SEGMENT } else { raw_flags::LAST_SEGMENT };
    let reverse_flag = if is_reverse { raw_flags::REVERSE } else { 0 };
    let mate_reverse_flag = if is_reverse { 0 } else { raw_flags::MATE_REVERSE };
    let flags = raw_flags::PAIRED
        | raw_flags::PROPER_PAIR
        | segment_flag
        | reverse_flag
        | mate_reverse_flag;

    // Single CIGAR op: {seq.len()}M (match). Consensus records are gap-free.
    let n = u32::try_from(seq.len()).expect("sequence length fits u32");
    // BAM CIGAR encoding: (length << 4) | op_code. op_code 0 = M (alignment match).
    let cigar_ops: Vec<u32> = if n > 0 { vec![n << 4] } else { Vec::new() };

    // Pre-compute bin for the alignment range (use unmapped bin for empty seq).
    let bin = if n > 0 {
        let alignment_start_1based =
            u32::try_from(local_pos + 1).expect("alignment start fits u32");
        let alignment_end_1based = alignment_start_1based + n - 1;
        region_to_bin(Some(alignment_start_1based), Some(alignment_end_1based))
    } else {
        region_to_bin(None, None)
    };

    let chrom_idx_i32 = i32::try_from(chrom_idx).expect("chrom_idx fits i32");
    let local_pos_i32 = i32::try_from(local_pos).expect("local_pos fits i32");

    let mut b = SamBuilder::new();
    b.read_name(name.as_bytes())
        .flags(flags)
        .ref_id(chrom_idx_i32)
        .pos(local_pos_i32)
        .mapq(60)
        .bin(bin)
        .mate_ref_id(chrom_idx_i32)
        .mate_pos(local_pos_i32) // R1 and R2 at same position
        .template_length(0)
        .cigar_ops(&cigar_ops)
        .sequence(seq)
        .qualities(quals);

    // `fgumi filter` masks bases using the PER-BASE depth/error arrays and reads
    // the read-level error tags (cE / aE / bE) as FLOAT rates — exactly what
    // `fgumi simplex`/`duplex` emit (see `vanilla_caller`/`duplex_caller`). Simplex
    // filtering reads the combined `CD_BASES`/`CE_BASES`; duplex filtering reads the
    // per-strand `AD_BASES`/`BD_BASES`. Without the per-base arrays every base reads
    // as zero depth and the read is rejected outright; with the error tags written
    // as integer counts, filter's `find_float_tag` skips them. Emit them faithfully
    // so the output honors its documented "suitable for input to `fgumi filter`"
    // contract.
    //
    // Each per-base depth array is a linear ramp from its max at the first base
    // down to its min at the last base, with the error count spread one-per-base.
    // The ramp — rather than "max everywhere with min only at the last base" — lets
    // a `fgumi filter --min-reads` threshold mask a realistic contiguous low-depth
    // *run* of interior bases instead of at most the single final base.
    //
    // Summary tags follow the real callers:
    //   * Simplex: `cD`/`cM` are the ramp max/min; `cE = sum(ce)/sum(cd)`.
    //   * Duplex: only the per-strand `aD`/`aM` (ramp max/min) and `bD`/`bM` are
    //     ramp endpoints. The combined `cD`/`cM`/`cE` are derived from the
    //     element-wise SUM of the two strand arrays — `cD = max(aD_i + bD_i)`,
    //     `cM = min(aD_i + bD_i)`, `cE = sum(ae + be) / sum(aD + bD)` — matching
    //     `duplex_caller`/`codec_caller`/fgbio, which emit no combined per-base
    //     array for duplex (see the `build_consensus_record` emission block below).
    //
    // The read-level error rate is the resulting `sum(errors)/sum(depths)` — the
    // formula the real callers emit (errors over total depth), NOT the simulator's
    // per-base `--error-rate-mean` input. See the `--error-rate-mean` arg doc for
    // the relationship to the requested rate.
    let read_len = seq.len();

    // Per-base depth/error arrays are written into reusable thread-local scratch
    // buffers rather than freshly-allocated `Vec`s. `build_consensus_record` runs
    // once per read (millions for a benchmark fixture) under `into_par_iter`, and
    // each record previously allocated two `Vec<i16>` for simplex plus four more
    // for duplex. Since `add_array_i16` copies the slice into the record's aux
    // bytes, the arrays are transient and the same buffers can be refilled. There
    // are two depth/error pairs — `(a_depths, a_errors, b_depths, b_errors)`:
    // simplex uses only the first pair, while duplex needs both strands live at
    // once to derive the combined `cD`/`cM`/`cE` summary from their element-wise
    // sum. The thread-local keeps the reuse safe across the parallel generation
    // loop (each rayon worker gets its own buffers).
    type PerBaseScratch = (Vec<i16>, Vec<i16>, Vec<i16>, Vec<i16>);
    thread_local! {
        static PER_BASE_SCRATCH: std::cell::RefCell<PerBaseScratch> =
            const {
                std::cell::RefCell::new((Vec::new(), Vec::new(), Vec::new(), Vec::new()))
            };
    }
    // Fills `depths`/`errors` (cleared first) with the per-base arrays for one
    // strand: depth is a linear ramp from `max_depth` (first base) to
    // `min_depth` (last base) so interior bases span the full [min, max] range;
    // error is one-per-base from the front for `error_count` bases. A
    // single-base read carries `min_depth` so it still honors the summary
    // minimum (it cannot carry both endpoints).
    let fill_per_base = |depths: &mut Vec<i16>,
                         errors: &mut Vec<i16>,
                         max_depth: i32,
                         min_depth: i32,
                         error_count: i32| {
        let clamp = |v: i32| i16::try_from(v).unwrap_or(i16::MAX);
        depths.clear();
        depths.reserve(read_len);
        match read_len {
            0 => {}
            1 => depths.push(clamp(min_depth)),
            _ => {
                // depth[i] = round(max - (max - min) * i / (read_len - 1)),
                // giving depth[0] == max and depth[read_len-1] == min.
                let span = i64::from(max_depth) - i64::from(min_depth);
                let last = (read_len - 1) as i64;
                for i in 0..read_len {
                    let v = i64::from(max_depth) - (span * i as i64 + last / 2) / last;
                    depths.push(clamp(i32::try_from(v).unwrap_or(i32::MAX)));
                }
            }
        }
        errors.clear();
        errors.resize(read_len, 0i16);
        let n = usize::try_from(error_count).unwrap_or(0).min(read_len);
        for slot in errors.iter_mut().take(n) {
            *slot = 1;
        }
    };
    let error_rate = |errors: &[i16], depths: &[i16]| -> f32 {
        let total_errors: i64 = errors.iter().map(|&e| i64::from(e)).sum();
        let total_depth: i64 = depths.iter().map(|&d| i64::from(d)).sum();
        if total_depth > 0 { (total_errors as f64 / total_depth as f64) as f32 } else { 0.0 }
    };

    PER_BASE_SCRATCH.with_borrow_mut(|(a_depths, a_errors, b_depths, b_errors)| {
        match duplex_tags {
            // Simplex: `cD`/`cM`/`cE` and the per-base `cd`/`ce` arrays all derive
            // from the single consensus depth/error ramp (matches `vanilla_caller`
            // and fgbio's `VanillaUmiConsensusCaller`, which DO emit per-base
            // `cd`/`ce` for single-strand consensus).
            None => {
                fill_per_base(a_depths, a_errors, cd, cm, ce);
                b.add_int_tag(SamTag::CD, cd)
                    .add_int_tag(SamTag::CM, cm)
                    .add_float_tag(SamTag::CE, error_rate(a_errors, a_depths));
                b.add_array_i16(SamTag::CD_BASES, a_depths)
                    .add_array_i16(SamTag::CE_BASES, a_errors);
            }
            // Duplex: `cD`/`cM`/`cE` are derived from the element-wise SUM of the
            // two strand depth/error arrays — `cD = max(aD_i + bD_i)`, `cM = min`,
            // `cE = total errors / total combined depth` — and ONLY the per-strand
            // `cd`/`ce` arrays (AD_BASES/AE_BASES, BD_BASES/BE_BASES) are emitted.
            // No combined per-base CD_BASES/CE_BASES is written. This mirrors the
            // real `duplex_caller`/`codec_caller` and fgbio's
            // `DuplexConsensusCaller`, which omit the per-base combined arrays for
            // duplex (`filter`'s duplex path reads the strand arrays). The passed
            // `cd`/`cm`/`ce` are intentionally unused here.
            Some((ad, bd, am, bm, ae, be)) => {
                fill_per_base(a_depths, a_errors, ad, am, ae);
                fill_per_base(b_depths, b_errors, bd, bm, be);

                // Combined per-base depth/error summary across both strands.
                let mut cd_max = 0i32;
                let mut cd_min = i32::MAX;
                let mut total_depth = 0i64;
                for (&a_d, &b_d) in a_depths.iter().zip(b_depths.iter()) {
                    let total = i32::from(a_d) + i32::from(b_d);
                    cd_max = cd_max.max(total);
                    cd_min = cd_min.min(total);
                    total_depth += i64::from(total);
                }
                let total_errors: i64 =
                    a_errors.iter().chain(b_errors.iter()).map(|&e| i64::from(e)).sum();
                // Empty read (read_len == 0): no per-base entries → zero summary.
                let cd_min = if a_depths.is_empty() { 0 } else { cd_min };
                let combined_rate = if total_depth > 0 {
                    (total_errors as f64 / total_depth as f64) as f32
                } else {
                    0.0
                };
                b.add_int_tag(SamTag::CD, cd_max)
                    .add_int_tag(SamTag::CM, cd_min)
                    .add_float_tag(SamTag::CE, combined_rate);

                // Per-strand summary (aD/bD/aM/bM int, aE/bE float rate) plus the
                // per-base strand depth/error arrays `filter` needs.
                b.add_int_tag(SamTag::AD, ad)
                    .add_int_tag(SamTag::BD, bd)
                    .add_int_tag(SamTag::AM, am)
                    .add_int_tag(SamTag::BM, bm);
                b.add_float_tag(SamTag::AE, error_rate(a_errors, a_depths));
                b.add_array_i16(SamTag::AD_BASES, a_depths)
                    .add_array_i16(SamTag::AE_BASES, a_errors);
                b.add_float_tag(SamTag::BE, error_rate(b_errors, b_depths));
                b.add_array_i16(SamTag::BD_BASES, b_depths)
                    .add_array_i16(SamTag::BE_BASES, b_errors);
            }
        }
    });

    // Add methylation tags if enabled.
    if let Some(meth) = methylation {
        // cu/ct tags (unconverted/converted counts per position).
        let cu: Vec<i16> = meth.annotation.unconverted_counts();
        let ct: Vec<i16> = meth.annotation.converted_counts();
        b.add_array_i16(SamTag::CU, &cu).add_array_i16(SamTag::CT, &ct);

        // MM/ML tags (SAM spec methylation tags).
        if let Some((mm, ml)) =
            build_mm_ml_tags(seq, &meth.annotation, is_top_strand, methylation_mode)
        {
            b.add_string_tag(SamTag::MM, mm.as_bytes()).add_array_u8(SamTag::ML, &ml);
        }

        // Duplex per-strand tags.
        if let (Some(ab), Some(ba)) = (&meth.ab_annotation, &meth.ba_annotation) {
            let au: Vec<i16> = ab.unconverted_counts();
            let at: Vec<i16> = ab.converted_counts();
            let bu: Vec<i16> = ba.unconverted_counts();
            let bt: Vec<i16> = ba.converted_counts();
            b.add_array_i16(SamTag::AU, &au)
                .add_array_i16(SamTag::AT, &at)
                .add_array_i16(SamTag::BU, &bu)
                .add_array_i16(SamTag::BT, &bt);

            // Per-strand MM tags (am/bm).
            if let Some(am) = fgumi_consensus::methylation::build_mm_tag_no_ml(
                seq,
                ab,
                is_top_strand,
                methylation_mode,
            ) {
                b.add_string_tag(SamTag::AM_BASES, am.as_bytes());
            }
            if let Some(bm) = fgumi_consensus::methylation::build_mm_tag_no_ml(
                seq,
                ba,
                is_top_strand,
                methylation_mode,
            ) {
                b.add_string_tag(SamTag::BM_BASES, bm.as_bytes());
            }
        }
    }

    b.build()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::commands::simulate::common::generate_random_sequence;
    use crate::simulate::create_rng;

    #[test]
    fn truth_tsv_header_names_error_columns_as_counts() {
        // The BAM emits read-level cE/aE/bE as FLOAT error rates, but the truth
        // TSV stores the simulated integer error *counts* (`ConsensusReadPair.
        // truth`). The error columns must be suffixed `_count` so the TSV never
        // silently disagrees with the BAM it is meant to validate; the depth
        // columns keep their bare tag names because they match the integer BAM
        // tags exactly.
        let cols: Vec<&str> = TRUTH_TSV_HEADER.split('\t').collect();
        assert_eq!(
            cols,
            vec![
                "read_name",
                "chrom",
                "pos",
                "strand",
                "cD",
                "cM",
                "cE_count",
                "aD",
                "bD",
                "aM",
                "bM",
                "aE_count",
                "bE_count",
            ],
        );
        // The nine data columns line up 1:1 with the 9-element `truth` tuple.
        assert_eq!(
            cols.len() - 4,
            9,
            "data columns (after read_name/chrom/pos/strand) must match the truth tuple arity"
        );
        // No bare rate-valued tag name leaks through as a count column.
        for bare in ["cE", "aE", "bE"] {
            assert!(
                !cols.contains(&bare),
                "bare `{bare}` would collide with the BAM's rate-valued tag; use `{bare}_count`"
            );
        }
    }

    /// Decode a raw BAM record into a noodles `RecordBuf` for higher-level test
    /// assertions.
    fn to_record_buf(raw: &RawRecord) -> noodles::sam::alignment::RecordBuf {
        fgumi_raw_bam::raw_record_to_record_buf(raw, &noodles::sam::Header::default())
            .expect("raw_record_to_record_buf failed in test")
    }

    #[test]
    fn test_generate_random_sequence_length() {
        let mut rng = create_rng(Some(42));
        for len in [0, 1, 8, 100, 300] {
            let seq = generate_random_sequence(len, &mut rng);
            assert_eq!(seq.len(), len);
        }
    }

    #[test]
    fn test_generate_random_sequence_valid_bases() {
        let mut rng = create_rng(Some(42));
        let seq = generate_random_sequence(1000, &mut rng);
        for &base in &seq {
            assert!(
                base == b'A' || base == b'C' || base == b'G' || base == b'T',
                "Invalid base: {}",
                base as char
            );
        }
    }

    #[test]
    fn test_reverse_complement_basic() {
        assert_eq!(reverse_complement(b"A"), b"T");
        assert_eq!(reverse_complement(b"T"), b"A");
        assert_eq!(reverse_complement(b"C"), b"G");
        assert_eq!(reverse_complement(b"G"), b"C");
    }

    #[test]
    fn test_reverse_complement_sequence() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
    }

    #[test]
    fn test_create_depth_distribution() {
        let dist = create_depth_distribution(5.0, 2.0);
        let mut rng = create_rng(Some(42));

        // Sample many times and check mean is reasonable
        let samples: Vec<f64> = (0..1000).map(|_| dist.sample(&mut rng)).collect();
        let mean: f64 = samples.iter().sum::<f64>() / samples.len() as f64;

        // Mean should be close to 5.0
        assert!(mean > 3.0 && mean < 7.0, "Mean {mean} not close to expected 5.0");
    }

    #[test]
    fn test_sample_depth_clamping() {
        let dist = create_depth_distribution(5.0, 2.0);
        let mut rng = create_rng(Some(42));

        for _ in 0..100 {
            let depth = sample_depth(&dist, 2, 8, &mut rng);
            assert!((2..=8).contains(&depth), "Depth {depth} out of range [2, 8]");
        }
    }

    #[test]
    fn test_build_consensus_record_simplex() {
        let seq = b"ACGTACGT";
        let quals = vec![40; 8];

        let raw = build_consensus_record(
            "test_read/1",
            seq,
            &quals,
            true,  // is_first
            false, // is_reverse
            0,     // chrom_idx
            100,   // local_pos
            5,     // cD
            3,     // cM
            1,     // cE
            None,  // no duplex tags
            None,
            MethylationMode::Disabled,
        );
        let record = to_record_buf(&raw);

        assert!(record.name().is_some());
        let flags = record.flags();
        assert!(!flags.is_unmapped());
        assert!(flags.is_first_segment());
        assert_eq!(record.reference_sequence_id(), Some(0));
    }

    /// Consensus reads must carry the per-base `cd`/`ce` arrays (`CD_BASES` /
    /// `CE_BASES`) in addition to the read-level `cD`/`cM`/`cE` — `fgumi filter`
    /// masks bases using them, and without them it rejects every read. Regression
    /// for the simulator omitting these arrays (its output was unusable by filter).
    #[test]
    fn test_build_consensus_record_emits_per_base_cd_ce_arrays() {
        let seq = b"ACGTACGTAC"; // length 10
        let quals = vec![40u8; seq.len()];
        let (cd, cm, ce) = (9_i32, 4_i32, 2_i32);

        let raw = build_consensus_record(
            "per_base/1",
            seq,
            &quals,
            true,
            false,
            0,
            100,
            cd,
            cm,
            ce,
            None,
            None,
            MethylationMode::Disabled,
        );

        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value as BufValue;
        use noodles::sam::alignment::record_buf::data::field::value::Array;
        let record = to_record_buf(&raw);
        let read_i16_array = |tag: SamTag| -> Vec<i16> {
            match record.data().get(&Tag::from(tag)) {
                Some(BufValue::Array(Array::Int16(vals))) => vals.clone(),
                other => panic!("expected an i16 array tag, got {other:?}"),
            }
        };
        let cd_bases = read_i16_array(SamTag::CD_BASES);
        let ce_bases = read_i16_array(SamTag::CE_BASES);

        // One entry per base.
        assert_eq!(cd_bases.len(), seq.len(), "per-base depth length == read length");
        assert_eq!(ce_bases.len(), seq.len(), "per-base error length == read length");

        // Per-base depth is consistent with the read-level summary tags.
        assert_eq!(i32::from(*cd_bases.iter().max().unwrap()), cd, "max per-base depth == cD");
        assert_eq!(i32::from(*cd_bases.iter().min().unwrap()), cm, "min per-base depth == cM");

        // Per-base errors sum to the read-level error count cE.
        let total_errors: i32 = ce_bases.iter().map(|&e| i32::from(e)).sum();
        let total_depth: i32 = cd_bases.iter().map(|&d| i32::from(d)).sum();
        assert_eq!(total_errors, ce, "per-base errors sum to cE");

        // The read-level cE tag is a FLOAT error rate (what `fgumi filter` reads
        // via `find_float_tag`), not an integer count, and equals sum(ce)/sum(cd).
        match record.data().get(&Tag::from(SamTag::CE)) {
            Some(BufValue::Float(rate)) => {
                let expected = total_errors as f64 / total_depth as f64;
                assert!(
                    (f64::from(*rate) - expected).abs() < 1e-6,
                    "cE should equal sum(ce)/sum(cd): got {rate}, expected {expected}"
                );
            }
            other => panic!("cE should be a float error rate, got {other:?}"),
        }
    }

    /// Duplex consensus reads must also carry the per-base strand depth/error
    /// arrays (aD/aE/bD/bE) and float aE/bE rates — `fgumi filter`'s duplex path
    /// masks bases from these arrays, so their absence makes it reject everything.
    #[test]
    fn test_build_consensus_record_duplex_emits_per_base_strand_arrays() {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value as BufValue;
        use noodles::sam::alignment::record_buf::data::field::value::Array;

        let seq = b"ACGTACGTAC"; // length 10
        let quals = vec![40u8; seq.len()];
        // cD/cM/cE then duplex (aD, bD, aM, bM, aE, bE).
        let raw = build_consensus_record(
            "dx/1",
            seq,
            &quals,
            true,
            false,
            0,
            100,
            10,
            6,
            3,
            Some((6, 4, 3, 2, 2, 1)),
            None,
            MethylationMode::Disabled,
        );
        let record = to_record_buf(&raw);
        let read_i16 = |tag: SamTag| -> Vec<i16> {
            match record.data().get(&Tag::from(tag)) {
                Some(BufValue::Array(Array::Int16(vals))) => vals.clone(),
                other => panic!("expected an i16 array tag, got {other:?}"),
            }
        };

        for tag in [SamTag::AD_BASES, SamTag::AE_BASES, SamTag::BD_BASES, SamTag::BE_BASES] {
            assert_eq!(read_i16(tag).len(), seq.len(), "per-base strand array spans the read");
        }
        // Strand depth arrays honor the per-strand max/min summary tags.
        let ad_bases = read_i16(SamTag::AD_BASES);
        assert_eq!(i32::from(*ad_bases.iter().max().unwrap()), 6, "max per-base aD == aD");
        assert_eq!(i32::from(*ad_bases.iter().min().unwrap()), 3, "min per-base aD == aM");
        let bd_bases = read_i16(SamTag::BD_BASES);
        assert_eq!(i32::from(*bd_bases.iter().max().unwrap()), 4, "max per-base bD == bD");
        assert_eq!(i32::from(*bd_bases.iter().min().unwrap()), 2, "min per-base bD == bM");
        // Per-strand errors sum to the strand error counts.
        assert_eq!(read_i16(SamTag::AE_BASES).iter().map(|&e| i32::from(e)).sum::<i32>(), 2);
        assert_eq!(read_i16(SamTag::BE_BASES).iter().map(|&e| i32::from(e)).sum::<i32>(), 1);
        // aE/bE are float rates derived from their per-strand arrays.
        for (rate_tag, error_tag, depth_tag) in [
            (SamTag::AE, SamTag::AE_BASES, SamTag::AD_BASES),
            (SamTag::BE, SamTag::BE_BASES, SamTag::BD_BASES),
        ] {
            let errors: i32 = read_i16(error_tag).iter().map(|&e| i32::from(e)).sum();
            let depths: i32 = read_i16(depth_tag).iter().map(|&d| i32::from(d)).sum();
            let expected = errors as f64 / depths as f64;
            match record.data().get(&Tag::from(rate_tag)) {
                Some(BufValue::Float(rate)) => assert!(
                    (f64::from(*rate) - expected).abs() < 1e-6,
                    "{rate_tag:?} should equal sum(errors)/sum(depths): got {rate}, expected {expected}"
                ),
                other => panic!("{rate_tag:?} should be a float error rate, got {other:?}"),
            }
        }
    }

    /// Duplex reads must derive `cD`/`cM`/`cE` from the element-wise SUM of the
    /// two strand depth arrays (`cD = max(aD_i + bD_i)`, `cM = min`, `cE = total
    /// errors / total depth`) and must NOT emit a combined per-base
    /// `CD_BASES`/`CE_BASES` array — exactly matching `duplex_caller`,
    /// `codec_caller`, and fgbio's `DuplexConsensusCaller`, which omit the
    /// per-base combined arrays for duplex (`filter`'s duplex path reads the
    /// strand arrays).
    #[test]
    fn test_build_consensus_record_duplex_cd_is_strand_sum() {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value as BufValue;
        use noodles::sam::alignment::record_buf::data::field::value::Array;

        let seq = b"ACGTACGTAC"; // length 10
        let quals = vec![40u8; seq.len()];
        // Intentionally inconsistent passed cd/cm/ce: the duplex path must IGNORE
        // them and recompute the summary from the strand arrays.
        let raw = build_consensus_record(
            "dx/1",
            seq,
            &quals,
            true,
            false,
            0,
            100,
            999,                      // bogus cd — must be ignored in duplex
            999,                      // bogus cm — must be ignored in duplex
            999,                      // bogus ce — must be ignored in duplex
            Some((6, 4, 3, 2, 2, 1)), // (ad, bd, am, bm, ae, be)
            None,
            MethylationMode::Disabled,
        );
        let record = to_record_buf(&raw);
        let int_tag = |tag: SamTag| -> i32 {
            match record.data().get(&Tag::from(*tag)) {
                Some(BufValue::Int32(v)) => *v,
                Some(BufValue::Int8(v)) => i32::from(*v),
                Some(BufValue::Int16(v)) => i32::from(*v),
                Some(BufValue::UInt8(v)) => i32::from(*v),
                Some(BufValue::UInt16(v)) => i32::from(*v),
                other => panic!("expected int tag for {:?}, got {other:?}", *tag),
            }
        };
        let float_tag = |tag: SamTag| -> f32 {
            match record.data().get(&Tag::from(*tag)) {
                Some(BufValue::Float(v)) => *v,
                other => panic!("expected float tag for {:?}, got {other:?}", *tag),
            }
        };
        let i16_array = |tag: SamTag| -> Vec<i16> {
            match record.data().get(&Tag::from(*tag)) {
                Some(BufValue::Array(Array::Int16(v))) => v.clone(),
                other => panic!("expected i16 array for {:?}, got {other:?}", *tag),
            }
        };

        // No combined per-base arrays in duplex (matches fgbio/duplex/codec callers).
        assert!(
            record.data().get(&Tag::from(*SamTag::CD_BASES)).is_none(),
            "CD_BASES must NOT be emitted in duplex mode"
        );
        assert!(
            record.data().get(&Tag::from(*SamTag::CE_BASES)).is_none(),
            "CE_BASES must NOT be emitted in duplex mode"
        );

        let ad = i16_array(SamTag::AD_BASES);
        let bd = i16_array(SamTag::BD_BASES);
        let ae = i16_array(SamTag::AE_BASES);
        let be = i16_array(SamTag::BE_BASES);
        let combined_depth: Vec<i32> =
            ad.iter().zip(&bd).map(|(&a, &b)| i32::from(a) + i32::from(b)).collect();

        // cD/cM are the max/min of the per-base strand-sum.
        assert_eq!(
            int_tag(SamTag::CD),
            *combined_depth.iter().max().unwrap(),
            "cD == max(aD_i + bD_i)"
        );
        assert_eq!(
            int_tag(SamTag::CM),
            *combined_depth.iter().min().unwrap(),
            "cM == min(aD_i + bD_i)"
        );
        // cE is total errors over total combined depth.
        let total_err: i32 = ae.iter().chain(&be).map(|&e| i32::from(e)).sum();
        let total_dep: i32 = combined_depth.iter().sum();
        let expected_ce = total_err as f64 / total_dep as f64;
        assert!(
            (f64::from(float_tag(SamTag::CE)) - expected_ce).abs() < 1e-6,
            "cE == total errors / total combined depth"
        );
    }

    #[test]
    fn test_validate_depth_bounds() {
        // Valid ranges pass.
        assert!(validate_depth_bounds(0, 0).is_ok());
        assert!(validate_depth_bounds(1, 10).is_ok());
        assert!(validate_depth_bounds(5, 5).is_ok());
        assert!(validate_depth_bounds(0, i32::from(i16::MAX)).is_ok());
        // Negative minimum is rejected.
        assert!(validate_depth_bounds(-1, 10).is_err());
        // max < min is rejected.
        assert!(validate_depth_bounds(10, 5).is_err());
        // max above i16::MAX is rejected (per-base arrays are i16; would saturate).
        assert!(validate_depth_bounds(1, i32::from(i16::MAX) + 1).is_err());
    }

    /// In duplex mode the per-strand minimum depths must sum to the total minimum
    /// (`aM + bM == cM`) and each must not exceed its strand max — so the emitted
    /// `cM` (derived from the strand-sum) matches the truth tuple. Exercised
    /// across seeds and a wide depth range so the sampled cM is usually < cD.
    #[test]
    fn test_generate_consensus_pair_duplex_strand_mins_sum_to_total() {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value as BufValue;
        use tempfile::NamedTempFile;

        let mut fasta = NamedTempFile::new().unwrap();
        writeln!(fasta, ">chr1").unwrap();
        fasta.write_all(&b"ACGT".repeat(500)).unwrap();
        writeln!(fasta).unwrap();
        fasta.flush().unwrap();
        let ref_genome = Arc::new(ReferenceGenome::load(fasta.path()).unwrap());
        let params = Arc::new(GenerationParams {
            read_length: 10,
            min_depth: 1,
            max_depth: 20,
            depth_mean: 10.0,
            depth_stddev: 5.0,
            error_rate_mean: 0.1,
            error_rate_stddev: 0.05,
            duplex: true,
            consensus_quality: 40,
            methylation_mode: MethylationMode::Disabled,
            methylation_depth_dist: create_depth_distribution(5.0, 2.5),
            cpg_methylation_rate: 0.75,
            conversion_rate: 0.98,
            ref_genome: Arc::clone(&ref_genome),
        });
        let strand_bias = StrandBiasModel::new(5.0, 5.0);

        let int_tag = |rec: &noodles::sam::alignment::RecordBuf, tag: SamTag| -> i32 {
            match rec.data().get(&Tag::from(*tag)) {
                Some(BufValue::Int32(v)) => *v,
                Some(BufValue::Int8(v)) => i32::from(*v),
                Some(BufValue::Int16(v)) => i32::from(*v),
                Some(BufValue::UInt8(v)) => i32::from(*v),
                Some(BufValue::UInt16(v)) => i32::from(*v),
                other => panic!("expected int tag for {:?}, got {other:?}", *tag),
            }
        };

        for seed in 0..64u64 {
            let pair = generate_consensus_pair(0, seed, &params, &strand_bias);
            let (cd, cm, _ce, ad, bd, am, bm, _ae, _be) = pair.truth;
            assert_eq!(ad + bd, cd, "truth aD + bD == cD (seed {seed})");
            assert_eq!(am + bm, cm, "truth aM + bM == cM (seed {seed})");
            assert!(am <= ad && bm <= bd, "each strand min <= its max (seed {seed})");

            // Emitted cD/cM (derived from the strand-sum) match the truth tuple.
            let r1 = to_record_buf(&pair.r1_record);
            assert_eq!(int_tag(&r1, SamTag::CD), cd, "BAM cD == truth cD (seed {seed})");
            assert_eq!(int_tag(&r1, SamTag::CM), cm, "BAM cM == truth cM (seed {seed})");
        }
    }

    #[test]
    fn test_build_consensus_record_duplex() {
        let seq = b"ACGTACGT";
        let quals = vec![40; 8];

        let raw = build_consensus_record(
            "test_read/1",
            seq,
            &quals,
            true,                     // is_first
            false,                    // is_reverse
            0,                        // chrom_idx
            100,                      // local_pos
            10,                       // cD
            5,                        // cM
            2,                        // cE
            Some((6, 4, 3, 2, 1, 1)), // aD, bD, aM, bM, aE, bE
            None,
            MethylationMode::Disabled,
        );
        let record = to_record_buf(&raw);

        assert!(record.name().is_some());
        let flags = record.flags();
        assert!(!flags.is_unmapped());
    }

    #[test]
    fn test_build_consensus_record_r2_flags() {
        let seq = b"ACGT";
        let quals = vec![40; 4];

        let raw = build_consensus_record(
            "test_read/2",
            seq,
            &quals,
            false, // is_first = false means R2
            true,  // is_reverse
            0,     // chrom_idx
            100,   // local_pos
            5,
            3,
            1,
            None,
            None,
            MethylationMode::Disabled,
        );
        let record = to_record_buf(&raw);

        let flags = record.flags();
        assert!(flags.is_last_segment());
        assert!(!flags.is_first_segment());
        assert!(!flags.is_unmapped());
    }

    #[test]
    fn test_depth_min_less_than_max() {
        let dist = create_depth_distribution(5.0, 2.0);
        let mut rng = create_rng(Some(42));

        for _ in 0..100 {
            let cd = sample_depth(&dist, 1, 10, &mut rng);
            let cm = sample_depth(&dist, 1, cd, &mut rng).min(cd);
            assert!(cm <= cd, "cM ({cm}) should be <= cD ({cd})");
        }
    }

    #[test]
    fn test_sample_depth_respects_min() {
        let dist = create_depth_distribution(5.0, 2.0);
        let mut rng = create_rng(Some(42));

        for _ in 0..100 {
            let depth = sample_depth(&dist, 3, 10, &mut rng);
            assert!(depth >= 3, "Depth {depth} should be >= min 3");
        }
    }

    #[test]
    fn test_sample_depth_respects_max() {
        let dist = create_depth_distribution(50.0, 10.0);
        let mut rng = create_rng(Some(42));

        for _ in 0..100 {
            let depth = sample_depth(&dist, 1, 5, &mut rng);
            assert!(depth <= 5, "Depth {depth} should be <= max 5");
        }
    }

    #[test]
    fn test_create_depth_distribution_high_mean() {
        let dist = create_depth_distribution(100.0, 20.0);
        let mut rng = create_rng(Some(42));

        let samples: Vec<f64> = (0..1000).map(|_| dist.sample(&mut rng)).collect();
        let mean: f64 = samples.iter().sum::<f64>() / samples.len() as f64;

        // Mean should be close to 100
        assert!(mean > 80.0 && mean < 120.0, "Mean {mean} not close to expected 100");
    }

    #[test]
    fn test_build_consensus_record_mapped_flags() {
        let seq = b"ACGT";
        let quals = vec![40; 4];

        let raw = build_consensus_record(
            "test/1",
            seq,
            &quals,
            true,  // is_first
            false, // is_reverse
            0,     // chrom_idx
            100,   // local_pos
            5,
            3,
            1,
            None,
            None,
            MethylationMode::Disabled,
        );
        let record = to_record_buf(&raw);

        let flags = record.flags();
        assert!(!flags.is_mate_unmapped());
        assert!(!flags.is_unmapped());
        assert_eq!(record.reference_sequence_id(), Some(0));
        assert!(record.alignment_start().is_some());
    }

    #[test]
    fn test_build_consensus_record_segmented_flag() {
        let seq = b"ACGT";
        let quals = vec![40; 4];

        let raw = build_consensus_record(
            "test/1",
            seq,
            &quals,
            true,  // is_first
            false, // is_reverse
            0,     // chrom_idx
            100,   // local_pos
            5,
            3,
            1,
            None,
            None,
            MethylationMode::Disabled,
        );
        let record = to_record_buf(&raw);

        let flags = record.flags();
        assert!(flags.is_segmented());
    }

    #[test]
    fn test_build_consensus_record_sequence_stored() {
        let seq = b"ACGTACGTACGT";
        let quals = vec![40; 12];

        let raw = build_consensus_record(
            "test/1",
            seq,
            &quals,
            true,  // is_first
            false, // is_reverse
            0,     // chrom_idx
            100,   // local_pos
            5,
            3,
            1,
            None,
            None,
            MethylationMode::Disabled,
        );
        let record = to_record_buf(&raw);

        // Verify sequence length matches
        assert_eq!(record.sequence().len(), 12);
    }

    #[test]
    fn test_build_consensus_record_qualities_stored() {
        let seq = b"ACGT";
        let quals = vec![10, 20, 30, 40];

        let raw = build_consensus_record(
            "test/1",
            seq,
            &quals,
            true,  // is_first
            false, // is_reverse
            0,     // chrom_idx
            100,   // local_pos
            5,
            3,
            1,
            None,
            None,
            MethylationMode::Disabled,
        );
        let record = to_record_buf(&raw);

        let record_quals: Vec<u8> = record.quality_scores().iter().collect();
        assert_eq!(record_quals, quals);
    }

    #[test]
    fn test_build_consensus_record_empty_sequence() {
        let seq: &[u8] = b"";
        let quals: Vec<u8> = vec![];

        let raw = build_consensus_record(
            "test/1",
            seq,
            &quals,
            true,  // is_first
            false, // is_reverse
            0,     // chrom_idx
            100,   // local_pos
            5,
            3,
            1,
            None,
            None,
            MethylationMode::Disabled,
        );
        let record = to_record_buf(&raw);

        assert!(record.name().is_some());
    }

    #[test]
    fn test_build_consensus_record_long_sequence() {
        let seq = vec![b'A'; 500];
        let quals = vec![40; 500];

        let raw = build_consensus_record(
            "test/1",
            &seq,
            &quals,
            true,  // is_first
            false, // is_reverse
            0,     // chrom_idx
            100,   // local_pos
            5,
            3,
            1,
            None,
            None,
            MethylationMode::Disabled,
        );
        let record = to_record_buf(&raw);

        assert_eq!(record.sequence().len(), 500);
    }

    #[test]
    fn test_build_consensus_record_zero_depth() {
        let seq = b"ACGT";
        let quals = vec![40; 4];

        // Edge case: zero depth (shouldn't happen normally but test it)
        let raw = build_consensus_record(
            "test/1",
            seq,
            &quals,
            true,  // is_first
            false, // is_reverse
            0,     // chrom_idx
            100,   // local_pos
            0,
            0,
            0,
            None,
            None,
            MethylationMode::Disabled,
        );
        let record = to_record_buf(&raw);

        assert!(record.name().is_some());
    }

    #[test]
    fn test_build_consensus_record_high_error_count() {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value as BufValue;
        use noodles::sam::alignment::record_buf::data::field::value::Array;

        let seq = b"ACGT";
        let quals = vec![40; 4];

        // High error count (cE=100, far more errors than the 4 bases): the
        // per-base error array saturates at one error per base.
        let (cd, cm, ce) = (10_i32, 5_i32, 100_i32);
        let raw = build_consensus_record(
            "test/1",
            seq,
            &quals,
            true,  // is_first
            false, // is_reverse
            0,     // chrom_idx
            100,   // local_pos
            cd,
            cm,
            ce,
            None,
            None,
            MethylationMode::Disabled,
        );
        let record = to_record_buf(&raw);
        assert!(record.name().is_some());

        let read_i16_array = |tag: SamTag| -> Vec<i16> {
            match record.data().get(&Tag::from(tag)) {
                Some(BufValue::Array(Array::Int16(vals))) => vals.clone(),
                other => panic!("expected an i16 array tag, got {other:?}"),
            }
        };
        let cd_bases = read_i16_array(SamTag::CD_BASES);
        let ce_bases = read_i16_array(SamTag::CE_BASES);

        // Errors saturate at one-per-base: every base carries an error.
        assert_eq!(ce_bases, vec![1, 1, 1, 1], "per-base errors saturate at one per base");
        let total_errors: i32 = ce_bases.iter().map(|&e| i32::from(e)).sum();
        let total_depth: i32 = cd_bases.iter().map(|&d| i32::from(d)).sum();
        assert_eq!(total_errors, seq.len() as i32, "saturated error count == read length");

        // The read-level cE tag is the saturated rate sum(ce)/sum(cd), matching
        // the real-caller formula (S8-001), pinned here so the saturation path is
        // tested rather than silently divergent.
        match record.data().get(&Tag::from(SamTag::CE)) {
            Some(BufValue::Float(rate)) => {
                let expected = f64::from(total_errors) / f64::from(total_depth);
                assert!(
                    (f64::from(*rate) - expected).abs() < 1e-6,
                    "cE should equal sum(ce)/sum(cd): got {rate}, expected {expected}"
                );
            }
            other => panic!("cE should be a float error rate, got {other:?}"),
        }
    }

    /// S8-002: the per-base depth array is a ramp spanning [cM, cD], so interior
    /// bases — not just the final one — fall below a `--min-reads` threshold
    /// between cM and cD. Endpoints honor the summary tags (max == cD at the
    /// first base, min == cM at the last).
    #[test]
    fn test_build_consensus_record_per_base_depth_is_a_ramp() {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value as BufValue;
        use noodles::sam::alignment::record_buf::data::field::value::Array;

        let seq = b"ACGTACGTAC"; // length 10
        let quals = vec![40u8; seq.len()];
        let (cd, cm, ce) = (20_i32, 4_i32, 1_i32);
        let raw = build_consensus_record(
            "ramp/1",
            seq,
            &quals,
            true,
            false,
            0,
            100,
            cd,
            cm,
            ce,
            None,
            None,
            MethylationMode::Disabled,
        );
        let record = to_record_buf(&raw);
        let cd_bases = match record.data().get(&Tag::from(SamTag::CD_BASES)) {
            Some(BufValue::Array(Array::Int16(vals))) => vals.clone(),
            other => panic!("expected CD_BASES i16 array, got {other:?}"),
        };

        assert_eq!(cd_bases.len(), seq.len());
        assert_eq!(i32::from(cd_bases[0]), cd, "first base depth == cD (max)");
        assert_eq!(i32::from(cd_bases[seq.len() - 1]), cm, "last base depth == cM (min)");
        // Monotonically non-increasing ramp (no flat max-everywhere plateau).
        for w in cd_bases.windows(2) {
            assert!(w[0] >= w[1], "depth ramp must be non-increasing: {w:?}");
        }
        // A threshold strictly between cM and cD masks MORE than one base — the
        // whole gap S8-002 was about (the old layout masked at most one).
        let threshold = cd.midpoint(cm);
        let masked = cd_bases.iter().filter(|&&d| i32::from(d) < threshold).count();
        assert!(masked > 1, "a mid threshold should mask a run of bases, masked={masked}");
    }

    /// A single-base read carries the summary minimum (it cannot carry both
    /// endpoints); the ramp must not panic on `read_len == 1`.
    #[test]
    fn test_build_consensus_record_single_base_depth_honors_min() {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value as BufValue;
        use noodles::sam::alignment::record_buf::data::field::value::Array;

        let seq = b"A";
        let quals = vec![40u8; 1];
        let (cd, cm, ce) = (9_i32, 3_i32, 0_i32);
        let raw = build_consensus_record(
            "single/1",
            seq,
            &quals,
            true,
            false,
            0,
            100,
            cd,
            cm,
            ce,
            None,
            None,
            MethylationMode::Disabled,
        );
        let record = to_record_buf(&raw);
        let cd_bases = match record.data().get(&Tag::from(SamTag::CD_BASES)) {
            Some(BufValue::Array(Array::Int16(vals))) => vals.clone(),
            other => panic!("expected CD_BASES i16 array, got {other:?}"),
        };
        assert_eq!(cd_bases, vec![3], "single base honors the summary minimum cM");
    }

    #[test]
    fn test_generate_random_sequence_reproducibility() {
        let mut rng1 = create_rng(Some(42));
        let mut rng2 = create_rng(Some(42));

        let seq1 = generate_random_sequence(100, &mut rng1);
        let seq2 = generate_random_sequence(100, &mut rng2);

        assert_eq!(seq1, seq2);
    }

    #[test]
    fn test_reverse_complement_double() {
        let seq = b"ACGTACGT";
        let rc = reverse_complement(seq);
        let rc_rc = reverse_complement(&rc);
        assert_eq!(rc_rc, seq.to_vec());
    }

    #[test]
    fn test_reverse_complement_empty() {
        let empty: &[u8] = b"";
        assert_eq!(reverse_complement(empty), Vec::<u8>::new());
    }

    #[test]
    fn test_reverse_complement_unknown_bases() {
        assert_eq!(reverse_complement(b"N"), b"N");
        assert_eq!(reverse_complement(b"X"), b"X");
        assert_eq!(reverse_complement(b"ANCG"), b"CGNT");
    }

    #[test]
    fn test_build_consensus_record_duplex_strand_depths_sum() {
        let seq = b"ACGT";
        let quals = vec![40; 4];

        // aD + bD should equal cD (approximately - due to rounding this might differ slightly)
        let cd = 10;
        let ad = 6;
        let bd = 4;

        let raw = build_consensus_record(
            "test/1",
            seq,
            &quals,
            true,  // is_first
            false, // is_reverse
            0,     // chrom_idx
            100,   // local_pos
            cd,
            5,
            2,
            Some((ad, bd, 3, 2, 1, 1)),
            None,
            MethylationMode::Disabled,
        );
        let record = to_record_buf(&raw);

        assert!(record.name().is_some());
        // The record should be valid with these tag values
    }

    #[test]
    fn test_build_consensus_record_various_names() {
        let seq = b"ACGT";
        let quals = vec![40; 4];

        for name in ["read1/1", "consensus_00000001/1", "test-read_123/2", "a/1"] {
            let raw = build_consensus_record(
                name,
                seq,
                &quals,
                true,
                false,
                0,
                100,
                5,
                3,
                1,
                None,
                None,
                MethylationMode::Disabled,
            );
            let record = to_record_buf(&raw);
            assert!(record.name().is_some());
        }
    }

    fn test_methylation_depth_dist() -> LogNormal<f64> {
        create_depth_distribution(5.0, 2.5)
    }

    #[test]
    fn test_generate_methylation_annotation_emseq() {
        // Sequence with CpG sites: positions 0(C),1(G) is CpG; position 4(C),5(A) is non-CpG
        let seq = b"CGAACAAT";
        let mut rng = create_rng(Some(42));
        let dist = test_methylation_depth_dist();

        let data = generate_methylation_annotation(
            seq,
            MethylationMode::EmSeq,
            1.0, // all CpG methylated
            0.98,
            &dist,
            false,
            &mut rng,
        );

        // Position 0 is C in CpG context -> methylated in EM-Seq -> mostly unconverted
        assert!(data.annotation.evidence[0].is_ref_c);
        assert!(data.annotation.evidence[0].unconverted_count > 0);

        // Position 4 is C but not CpG -> unmethylated in EM-Seq -> mostly converted
        assert!(data.annotation.evidence[4].is_ref_c);
        assert!(
            data.annotation.evidence[4].converted_count
                >= data.annotation.evidence[4].unconverted_count
        );

        // Position 1 is G -> not a ref C
        assert!(!data.annotation.evidence[1].is_ref_c);

        // cu/ct should have correct lengths
        let cu = data.annotation.unconverted_counts();
        let ct = data.annotation.converted_counts();
        assert_eq!(cu.len(), seq.len());
        assert_eq!(ct.len(), seq.len());
    }

    #[test]
    fn test_generate_methylation_annotation_taps() {
        let seq = b"CGAACAAT";
        let mut rng = create_rng(Some(42));
        let dist = test_methylation_depth_dist();

        let data = generate_methylation_annotation(
            seq,
            MethylationMode::Taps,
            1.0, // all CpG methylated
            0.98,
            &dist,
            false,
            &mut rng,
        );

        // Position 0 is C in CpG -> methylated in TAPs -> mostly converted (target)
        assert!(data.annotation.evidence[0].is_ref_c);
        assert!(
            data.annotation.evidence[0].converted_count
                >= data.annotation.evidence[0].unconverted_count
        );

        // Position 4 is C not CpG -> unmethylated in TAPs -> mostly unconverted (not a target)
        assert!(data.annotation.evidence[4].is_ref_c);
        assert!(
            data.annotation.evidence[4].unconverted_count
                >= data.annotation.evidence[4].converted_count
        );
    }

    #[test]
    fn test_build_consensus_record_with_methylation() {
        let seq = b"CGAACGAT";
        let quals = vec![40; 8];
        let mut rng = create_rng(Some(42));
        let dist = test_methylation_depth_dist();

        let meth = generate_methylation_annotation(
            seq,
            MethylationMode::EmSeq,
            0.75,
            0.98,
            &dist,
            false,
            &mut rng,
        );

        let raw = build_consensus_record(
            "test/1",
            seq,
            &quals,
            true,  // is_first
            false, // is_reverse
            0,     // chrom_idx
            100,   // local_pos
            5,
            3,
            1,
            None,
            Some(&meth),
            MethylationMode::EmSeq,
        );
        let record = to_record_buf(&raw);

        // Verify methylation tags are present
        assert!(record.name().is_some());
        let cu: Vec<i16> = meth.annotation.unconverted_counts();
        let ct: Vec<i16> = meth.annotation.converted_counts();
        assert_eq!(cu.len(), seq.len());
        assert_eq!(ct.len(), seq.len());

        // Verify at C positions cu+ct > 0
        for (i, &b) in seq.iter().enumerate() {
            if b == b'C' {
                assert!(cu[i] + ct[i] > 0, "position {i}: C should have non-zero cu+ct");
            } else {
                assert_eq!(cu[i], 0, "position {i}: non-C should have cu=0");
                assert_eq!(ct[i], 0, "position {i}: non-C should have ct=0");
            }
        }
    }

    #[test]
    fn test_build_consensus_record_duplex_with_methylation() {
        let seq = b"CGAACGAT";
        let quals = vec![40; 8];
        let mut rng = create_rng(Some(42));
        let dist = test_methylation_depth_dist();

        let meth = generate_methylation_annotation(
            seq,
            MethylationMode::EmSeq,
            0.75,
            0.98,
            &dist,
            true, // duplex
            &mut rng,
        );

        // Verify duplex methylation data was generated
        let ab = meth.ab_annotation.as_ref().expect("AB annotation should be present");
        let ba = meth.ba_annotation.as_ref().expect("BA annotation should be present");

        // Verify per-strand counts have correct lengths
        assert_eq!(ab.unconverted_counts().len(), seq.len());
        assert_eq!(ba.converted_counts().len(), seq.len());

        let raw = build_consensus_record(
            "test/1",
            seq,
            &quals,
            true,  // is_first
            false, // is_reverse
            0,     // chrom_idx
            100,   // local_pos
            10,
            5,
            2,
            Some((6, 4, 3, 2, 1, 1)),
            Some(&meth),
            MethylationMode::EmSeq,
        );
        let record = to_record_buf(&raw);

        assert!(record.name().is_some());
    }

    #[test]
    fn test_no_methylation_tags_when_disabled() {
        let seq = b"CGAACGAT";
        let quals = vec![40; 8];

        // When no methylation data is passed, record should still build fine
        let raw = build_consensus_record(
            "test/1",
            seq,
            &quals,
            true,  // is_first
            false, // is_reverse
            0,     // chrom_idx
            100,   // local_pos
            5,
            3,
            1,
            None,
            None,
            MethylationMode::Disabled,
        );
        let record = to_record_buf(&raw);

        assert!(record.name().is_some());
    }

    #[test]
    fn test_methylation_depth_mean_validation_rejects_zero() {
        let cmd = ConsensusReads {
            output: PathBuf::from("/dev/null"),
            truth_output: None,
            num_reads: 1,
            read_length: 10,
            seed: Some(42),
            threads: 1,
            compression: CompressionOptions { compression_level: 1 },
            min_depth: 1,
            max_depth: 10,
            depth_mean: 5.0,
            depth_stddev: 2.0,
            error_rate_mean: 0.01,
            error_rate_stddev: 0.005,
            duplex: false,
            consensus_quality: 40,
            strand_bias: StrandBiasArgs { strand_alpha: 5.0, strand_beta: 5.0 },
            methylation: MethylationArgs {
                methylation_mode: Some(crate::commands::common::MethylationModeArg::EmSeq),
                cpg_methylation_rate: 0.75,
                conversion_rate: 0.98,
            },
            methylation_depth_mean: 0.0,
            reference: PathBuf::from("dummy.fa"),
        };
        let result = cmd.execute("test");
        assert!(result.is_err());
        let msg = result.unwrap_err().to_string();
        assert!(
            msg.contains("methylation-depth-mean"),
            "error should mention methylation-depth-mean, got: {msg}"
        );
    }

    #[test]
    fn test_methylation_depth_mean_validation_rejects_negative() {
        let cmd = ConsensusReads {
            output: PathBuf::from("/dev/null"),
            truth_output: None,
            num_reads: 1,
            read_length: 10,
            seed: Some(42),
            threads: 1,
            compression: CompressionOptions { compression_level: 1 },
            min_depth: 1,
            max_depth: 10,
            depth_mean: 5.0,
            depth_stddev: 2.0,
            error_rate_mean: 0.01,
            error_rate_stddev: 0.005,
            duplex: false,
            consensus_quality: 40,
            strand_bias: StrandBiasArgs { strand_alpha: 5.0, strand_beta: 5.0 },
            methylation: MethylationArgs {
                methylation_mode: Some(crate::commands::common::MethylationModeArg::EmSeq),
                cpg_methylation_rate: 0.75,
                conversion_rate: 0.98,
            },
            methylation_depth_mean: -1.0,
            reference: PathBuf::from("dummy.fa"),
        };
        let result = cmd.execute("test");
        assert!(result.is_err());
    }

    #[test]
    fn test_methylation_depth_mean_validation_rejects_nan() {
        let cmd = ConsensusReads {
            output: PathBuf::from("/dev/null"),
            truth_output: None,
            num_reads: 1,
            read_length: 10,
            seed: Some(42),
            threads: 1,
            compression: CompressionOptions { compression_level: 1 },
            min_depth: 1,
            max_depth: 10,
            depth_mean: 5.0,
            depth_stddev: 2.0,
            error_rate_mean: 0.01,
            error_rate_stddev: 0.005,
            duplex: false,
            consensus_quality: 40,
            strand_bias: StrandBiasArgs { strand_alpha: 5.0, strand_beta: 5.0 },
            methylation: MethylationArgs {
                methylation_mode: Some(crate::commands::common::MethylationModeArg::EmSeq),
                cpg_methylation_rate: 0.75,
                conversion_rate: 0.98,
            },
            methylation_depth_mean: f64::NAN,
            reference: PathBuf::from("dummy.fa"),
        };
        let result = cmd.execute("test");
        assert!(result.is_err());
    }

    #[test]
    fn test_methylation_depth_mean_validation_accepts_when_disabled() {
        use std::io::Write as IoWrite;
        use tempfile::NamedTempFile;

        let mut fasta = NamedTempFile::new().unwrap();
        writeln!(fasta, ">chr1").unwrap();
        fasta.write_all(&b"ACGT".repeat(500)).unwrap();
        writeln!(fasta).unwrap();
        fasta.flush().unwrap();

        // When methylation is disabled, invalid methylation_depth_mean should not error
        let cmd = ConsensusReads {
            output: PathBuf::from("/dev/null"),
            truth_output: None,
            num_reads: 0,
            read_length: 10,
            seed: Some(42),
            threads: 1,
            compression: CompressionOptions { compression_level: 1 },
            min_depth: 1,
            max_depth: 10,
            depth_mean: 5.0,
            depth_stddev: 2.0,
            error_rate_mean: 0.01,
            error_rate_stddev: 0.005,
            duplex: false,
            consensus_quality: 40,
            strand_bias: StrandBiasArgs { strand_alpha: 5.0, strand_beta: 5.0 },
            methylation: MethylationArgs {
                methylation_mode: None,
                cpg_methylation_rate: 0.75,
                conversion_rate: 0.98,
            },
            methylation_depth_mean: 0.0,
            reference: fasta.path().to_path_buf(),
        };
        // Should succeed since methylation is disabled
        let result = cmd.execute("test");
        assert!(result.is_ok(), "disabled methylation should not validate depth mean");
    }

    #[test]
    fn test_methylation_annotation_reverse() {
        let annotation = MethylationAnnotation {
            evidence: vec![
                MethylationEvidence { is_ref_c: true, unconverted_count: 10, converted_count: 2 },
                MethylationEvidence { is_ref_c: false, unconverted_count: 0, converted_count: 0 },
                MethylationEvidence { is_ref_c: true, unconverted_count: 3, converted_count: 7 },
            ],
        };
        let reversed = annotation.reverse();
        assert_eq!(reversed.evidence.len(), 3);
        assert!(reversed.evidence[0].is_ref_c);
        assert_eq!(reversed.evidence[0].unconverted_count, 3);
        assert_eq!(reversed.evidence[0].converted_count, 7);
        assert!(!reversed.evidence[1].is_ref_c);
        assert!(reversed.evidence[2].is_ref_c);
        assert_eq!(reversed.evidence[2].unconverted_count, 10);
        assert_eq!(reversed.evidence[2].converted_count, 2);
    }

    #[test]
    fn test_methylation_data_reverse() {
        let annotation = MethylationAnnotation {
            evidence: vec![
                MethylationEvidence { is_ref_c: true, unconverted_count: 10, converted_count: 2 },
                MethylationEvidence { is_ref_c: false, unconverted_count: 0, converted_count: 0 },
            ],
        };
        let ab = MethylationAnnotation {
            evidence: vec![
                MethylationEvidence { is_ref_c: true, unconverted_count: 6, converted_count: 1 },
                MethylationEvidence { is_ref_c: false, unconverted_count: 0, converted_count: 0 },
            ],
        };
        let ba = MethylationAnnotation {
            evidence: vec![
                MethylationEvidence { is_ref_c: true, unconverted_count: 4, converted_count: 1 },
                MethylationEvidence { is_ref_c: false, unconverted_count: 0, converted_count: 0 },
            ],
        };
        let data = MethylationData { annotation, ab_annotation: Some(ab), ba_annotation: Some(ba) };
        let reversed = data.reverse();

        // Main annotation reversed
        assert!(!reversed.annotation.evidence[0].is_ref_c);
        assert!(reversed.annotation.evidence[1].is_ref_c);
        assert_eq!(reversed.annotation.evidence[1].unconverted_count, 10);

        // AB reversed
        let ab_rev = reversed.ab_annotation.unwrap();
        assert!(!ab_rev.evidence[0].is_ref_c);
        assert!(ab_rev.evidence[1].is_ref_c);
        assert_eq!(ab_rev.evidence[1].unconverted_count, 6);

        // BA reversed
        let ba_rev = reversed.ba_annotation.unwrap();
        assert!(!ba_rev.evidence[0].is_ref_c);
        assert!(ba_rev.evidence[1].is_ref_c);
        assert_eq!(ba_rev.evidence[1].unconverted_count, 4);
    }

    #[test]
    fn test_r2_methylation_cu_ct_tags_are_reversed() {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value;
        use noodles::sam::alignment::record_buf::data::field::value::Array;

        // Create a sequence with C positions at known offsets:
        // "CGAT" has C at pos 0 (CpG context) and nothing else
        let seq = b"CGAT";
        let quals = vec![40; 4];

        // Build methylation annotation with distinct counts at each position
        let annotation = MethylationAnnotation {
            evidence: vec![
                MethylationEvidence { is_ref_c: true, unconverted_count: 10, converted_count: 2 },
                MethylationEvidence { is_ref_c: false, unconverted_count: 0, converted_count: 0 },
                MethylationEvidence { is_ref_c: false, unconverted_count: 0, converted_count: 0 },
                MethylationEvidence { is_ref_c: false, unconverted_count: 0, converted_count: 0 },
            ],
        };
        let meth = MethylationData { annotation, ab_annotation: None, ba_annotation: None };

        // Build R1 record (uses original methylation)
        let r1_raw = build_consensus_record(
            "test/1",
            seq,
            &quals,
            true,  // is_first
            false, // is_reverse
            0,     // chrom_idx
            100,   // local_pos
            5,
            3,
            1,
            None,
            Some(&meth),
            MethylationMode::EmSeq,
        );
        let r1 = to_record_buf(&r1_raw);

        // Build R2 record with reversed methylation (as generate_consensus_pair does)
        let r2_seq = reverse_complement(seq);
        let r2_meth = meth.reverse();
        let r2_raw = build_consensus_record(
            "test/2",
            &r2_seq,
            &quals,
            false, // is_first
            true,  // is_reverse
            0,     // chrom_idx
            100,   // local_pos
            5,
            3,
            1,
            None,
            Some(&r2_meth),
            MethylationMode::EmSeq,
        );
        let r2 = to_record_buf(&r2_raw);

        // Extract cu tags from both records
        let cu_tag = Tag::from([b'c', b'u']);
        let ct_tag = Tag::from([b'c', b't']);

        let r1_cu = match r1.data().get(&cu_tag) {
            Some(Value::Array(Array::Int16(arr))) => arr.clone(),
            other => panic!("Expected Int16 array for cu tag on R1, got {other:?}"),
        };
        let r1_ct = match r1.data().get(&ct_tag) {
            Some(Value::Array(Array::Int16(arr))) => arr.clone(),
            other => panic!("Expected Int16 array for ct tag on R1, got {other:?}"),
        };

        let r2_cu = match r2.data().get(&cu_tag) {
            Some(Value::Array(Array::Int16(arr))) => arr.clone(),
            other => panic!("Expected Int16 array for cu tag on R2, got {other:?}"),
        };
        let r2_ct = match r2.data().get(&ct_tag) {
            Some(Value::Array(Array::Int16(arr))) => arr.clone(),
            other => panic!("Expected Int16 array for ct tag on R2, got {other:?}"),
        };

        // R2's cu/ct should be the reverse of R1's cu/ct
        let mut r1_cu_rev = r1_cu.clone();
        r1_cu_rev.reverse();
        assert_eq!(
            r2_cu, r1_cu_rev,
            "R2 cu tag should be reversed relative to R1: R1={r1_cu:?}, R2={r2_cu:?}"
        );

        let mut r1_ct_rev = r1_ct.clone();
        r1_ct_rev.reverse();
        assert_eq!(
            r2_ct, r1_ct_rev,
            "R2 ct tag should be reversed relative to R1: R1={r1_ct:?}, R2={r2_ct:?}"
        );
    }

    #[test]
    fn test_r2_per_base_depth_error_arrays_are_read_oriented_not_reversed() {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value;
        use noodles::sam::alignment::record_buf::data::field::value::Array;

        // Unlike the methylation cu/ct arrays (which ARE reversed for R2 because
        // they annotate reference-oriented evidence), the per-base depth/error
        // arrays cd/ce (and the duplex ad/ae/bd/be) are built fresh from the
        // scalar summaries and are *read-oriented*: each read's "low-depth tail"
        // is anchored at its own last base. `fgumi simplex`/`duplex` emit cd/ce
        // in read order, so the simulator deliberately does NOT reverse them for
        // R2. This test pins that intended orientation so a future change that
        // started reversing them (or stopped anchoring the min at the last
        // index) is caught.
        let seq = b"ACGTA";
        let quals = vec![40; seq.len()];

        // Distinct max/min so the anchored min is identifiable, and a non-zero
        // error count so the per-base error array is non-trivial. Duplex tags are
        // present so ad/ae/bd/be are emitted too.
        let cd = 8; // max depth (array fill)
        let cm = 3; // min depth (anchored at last base)
        let ce = 2; // error count (spread one-per-base from the front)
        let duplex_tags = (7, 6, 2, 1, 1, 1); // (ad, bd, am, bm, ae, be)
        let duplex = Some(duplex_tags);

        let r1_raw = build_consensus_record(
            "test/1",
            seq,
            &quals,
            true,
            false,
            0,
            100,
            cd,
            cm,
            ce,
            duplex,
            None,
            MethylationMode::Disabled,
        );
        let r1 = to_record_buf(&r1_raw);

        // R2 is the reverse complement, matching generate_consensus_pair.
        let r2_seq = reverse_complement(seq);
        let r2_raw = build_consensus_record(
            "test/2",
            &r2_seq,
            &quals,
            false,
            true,
            0,
            100,
            cd,
            cm,
            ce,
            duplex,
            None,
            MethylationMode::Disabled,
        );
        let r2 = to_record_buf(&r2_raw);

        let array_tag = |rec: &noodles::sam::alignment::RecordBuf, tag: SamTag| -> Vec<i16> {
            match rec.data().get(&Tag::from(*tag)) {
                Some(Value::Array(Array::Int16(arr))) => arr.clone(),
                other => panic!("expected Int16 array for tag {:?}, got {other:?}", *tag),
            }
        };

        // Every per-base array must be IDENTICAL between R1 and R2 (read-oriented,
        // not reversed), and the depth array's min must sit at the LAST index on
        // both reads. Duplex mode emits only the per-strand arrays (no combined
        // CD_BASES/CE_BASES — see `build_consensus_record`).
        for tag in [SamTag::AD_BASES, SamTag::AE_BASES, SamTag::BD_BASES, SamTag::BE_BASES] {
            let r1_arr = array_tag(&r1, tag);
            let r2_arr = array_tag(&r2, tag);
            assert_eq!(
                r1_arr, r2_arr,
                "per-base array {:?} must be read-oriented (same for R1 and R2), \
                 not reversed: R1={r1_arr:?}, R2={r2_arr:?}",
                *tag
            );
            let mut reversed = r1_arr.clone();
            reversed.reverse();
            // Sanity: with our asymmetric inputs the array is not palindromic, so
            // "same as R1" is a strictly stronger claim than "reverse of R1".
            assert_ne!(
                r1_arr, reversed,
                "test input for {:?} must be non-palindromic to be meaningful",
                *tag
            );
        }

        // The strand depth arrays anchor their min at the last base on BOTH reads
        // (aM/bM are the per-strand minima), so the expected last-base values are
        // `am` and `bm` from the duplex tuple `(ad, bd, am, bm, ae, be)`.
        let (_ad, _bd, am, bm, _ae, _be) = duplex_tags;
        let r1_ad = array_tag(&r1, SamTag::AD_BASES);
        let r2_ad = array_tag(&r2, SamTag::AD_BASES);
        assert_eq!(*r1_ad.last().unwrap(), am as i16, "R1 aD min anchored at last base");
        assert_eq!(*r2_ad.last().unwrap(), am as i16, "R2 aD min anchored at last base");
        let r1_bd = array_tag(&r1, SamTag::BD_BASES);
        let r2_bd = array_tag(&r2, SamTag::BD_BASES);
        assert_eq!(*r1_bd.last().unwrap(), bm as i16, "R1 bD min anchored at last base");
        assert_eq!(*r2_bd.last().unwrap(), bm as i16, "R2 bD min anchored at last base");
    }

    #[test]
    fn test_consensus_reads_produces_mapped_records() {
        use std::io::Write as IoWrite;
        use tempfile::NamedTempFile;

        let mut fasta = NamedTempFile::new().unwrap();
        writeln!(fasta, ">chr1").unwrap();
        fasta.write_all(&b"ACGT".repeat(500)).unwrap();
        writeln!(fasta).unwrap();
        fasta.flush().unwrap();

        let ref_genome = Arc::new(ReferenceGenome::load(fasta.path()).unwrap());
        let params = Arc::new(GenerationParams {
            read_length: 50,
            min_depth: 1,
            max_depth: 10,
            depth_mean: 5.0,
            depth_stddev: 2.0,
            error_rate_mean: 0.01,
            error_rate_stddev: 0.005,
            duplex: false,
            consensus_quality: 40,
            methylation_mode: MethylationMode::Disabled,
            methylation_depth_dist: create_depth_distribution(5.0, 2.5),
            cpg_methylation_rate: 0.75,
            conversion_rate: 0.98,
            ref_genome: Arc::clone(&ref_genome),
        });
        let strand_bias = StrandBiasModel::new(5.0, 5.0);

        let pair = generate_consensus_pair(0, 42, &params, &strand_bias);
        let r1 = to_record_buf(&pair.r1_record);
        let r2 = to_record_buf(&pair.r2_record);

        // Records should be mapped
        assert!(!r1.flags().is_unmapped());
        assert!(!r2.flags().is_unmapped());

        // Records should have reference sequence ID
        assert!(r1.reference_sequence_id().is_some());
        assert!(r2.reference_sequence_id().is_some());

        // Records should have alignment start
        assert!(r1.alignment_start().is_some());
        assert!(r2.alignment_start().is_some());

        // Chrom name should be set
        assert_eq!(pair.chrom_name, "chr1");
    }

    #[test]
    fn single_base_read_collapses_summary_depths() {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value as BufValue;
        use noodles::sam::alignment::record_buf::data::field::value::Array;
        use std::io::Write as IoWrite;
        use tempfile::NamedTempFile;

        let mut fasta = NamedTempFile::new().unwrap();
        writeln!(fasta, ">chr1").unwrap();
        fasta.write_all(&b"ACGT".repeat(500)).unwrap();
        writeln!(fasta).unwrap();
        fasta.flush().unwrap();
        let ref_genome = Arc::new(ReferenceGenome::load(fasta.path()).unwrap());

        // A wide depth range so the sampled cD would normally exceed cM.
        let params = Arc::new(GenerationParams {
            read_length: 1,
            min_depth: 1,
            max_depth: 20,
            depth_mean: 10.0,
            depth_stddev: 5.0,
            error_rate_mean: 0.2,
            error_rate_stddev: 0.05,
            duplex: true,
            consensus_quality: 40,
            methylation_mode: MethylationMode::Disabled,
            methylation_depth_dist: create_depth_distribution(5.0, 2.5),
            cpg_methylation_rate: 0.75,
            conversion_rate: 0.98,
            ref_genome: Arc::clone(&ref_genome),
        });
        let strand_bias = StrandBiasModel::new(5.0, 5.0);

        let int_eq = |rec: &noodles::sam::alignment::RecordBuf, a: SamTag, b: SamTag| {
            rec.data().get(&Tag::from(*a)) == rec.data().get(&Tag::from(*b))
        };
        let i16_array = |rec: &noodles::sam::alignment::RecordBuf, tag: SamTag| -> Vec<i16> {
            match rec.data().get(&Tag::from(*tag)) {
                Some(BufValue::Array(Array::Int16(v))) => v.clone(),
                other => panic!("expected i16 array for {:?}, got {other:?}", *tag),
            }
        };

        // A one-base read cannot carry distinct max/min depths. Across seeds, each
        // summary max must equal its min — in the truth tuple and the BAM tags —
        // and the lone per-base depth must match the collapsed summary.
        for seed in 0..64u64 {
            let pair = generate_consensus_pair(0, seed, &params, &strand_bias);
            let (cd, cm, _ce, ad, bd, am, bm, _ae, _be) = pair.truth;
            assert_eq!(cd, cm, "truth cD must equal cM for a 1-base read (seed {seed})");
            assert_eq!(ad, am, "truth aD must equal aM for a 1-base read (seed {seed})");
            assert_eq!(bd, bm, "truth bD must equal bM for a 1-base read (seed {seed})");

            let r1 = to_record_buf(&pair.r1_record);
            assert!(int_eq(&r1, SamTag::CD, SamTag::CM), "BAM cD==cM (seed {seed})");
            assert!(int_eq(&r1, SamTag::AD, SamTag::AM), "BAM aD==aM (seed {seed})");
            assert!(int_eq(&r1, SamTag::BD, SamTag::BM), "BAM bD==bM (seed {seed})");
            // Duplex omits the combined CD_BASES; the lone per-base entry lives on
            // each strand array, and aM + bM == cM == the combined single depth.
            assert!(
                r1.data().get(&Tag::from(*SamTag::CD_BASES)).is_none(),
                "duplex must not emit CD_BASES (seed {seed})"
            );
            assert_eq!(
                i16_array(&r1, SamTag::AD_BASES),
                vec![am as i16],
                "single per-base aD must match the collapsed strand min (seed {seed})"
            );
            assert_eq!(
                i16_array(&r1, SamTag::BD_BASES),
                vec![bm as i16],
                "single per-base bD must match the collapsed strand min (seed {seed})"
            );
            assert_eq!(am + bm, cm, "aM + bM == cM for a 1-base duplex read (seed {seed})");
        }
    }
}

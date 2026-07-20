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

    /// Minimum consensus depth (cM tag)
    #[arg(long = "min-depth", default_value = "1")]
    pub min_depth: i32,

    /// Maximum consensus depth (cD tag)
    #[arg(long = "max-depth", default_value = "10")]
    pub max_depth: i32,

    /// Mean depth for sampling
    #[arg(long = "depth-mean", default_value = "5.0")]
    pub depth_mean: f64,

    /// Depth standard deviation
    #[arg(long = "depth-stddev", default_value = "2.0")]
    pub depth_stddev: f64,

    /// Mean error rate (cE tag)
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

/// A generated consensus read pair ready for output.
struct ConsensusReadPair {
    read_name: String,
    r1_record: RawRecord,
    r2_record: RawRecord,
    /// Truth data: (cD, cM, cE, aD, bD, aM, bM, aE, bE)
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
                writeln!(w, "read_name\tchrom\tpos\tstrand\tcD\tcM\tcE\taD\tbD\taM\tbM\taE\tbE")?;
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

    // Generate error count based on error rate
    let error_rate = error_dist.sample(&mut rng).clamp(0.0, 1.0);
    let ce = (params.read_length as f64 * error_rate).round() as i32;

    // For duplex mode, generate strand-specific depths
    let (ad, bd, am, bm, ae, be) = if params.duplex {
        // Split total depth between A and B strands
        let a_frac = strand_bias_model.sample_a_fraction(&mut rng);

        let ad = ((cd as f64) * a_frac).round() as i32;
        let bd = cd - ad;

        // Min depths for each strand. These must satisfy aM + bM == cM, because a
        // duplex consensus position's combined depth is the sum of its two strand
        // depths — the real duplex caller derives cM as the per-base minimum of
        // (ab_i + ba_i). Splitting cM by the same strand fraction as cD keeps the
        // truth file and the emitted tags mutually consistent.
        //
        // The split is feasible because cM <= cD == aD + bD; clamping the A share to
        // [cM - bD, min(aD, cM)] guarantees both aM <= aD and bM = cM - aM <= bD.
        let am = (((cm as f64) * a_frac).round() as i32).clamp((cm - bd).max(0), ad.min(cm));
        let bm = cm - am;

        // Errors distributed proportionally
        let ae = ((ce as f64) * a_frac).round() as i32;
        let be = ce - ae;

        (ad, bd, am, bm, ae, be)
    } else {
        (0, 0, 0, 0, 0, 0)
    };

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
    // `fgumi simplex`/`duplex` emit (see `vanilla_caller`). Without the per-base
    // arrays every base reads as zero depth and the read is rejected outright;
    // with the error tags written as integer counts, filter's `find_float_tag`
    // skips them. Emit both faithfully so the output honors its documented
    // "suitable for input to `fgumi filter`" contract.
    //
    // Build per-base arrays consistent with the read-level summary tags: depth
    // spans [min, max] (min == cM/aM/bM, max == cD/aD/bD) and the total error
    // count is spread one-per-base. The read-level error tag is the resulting
    // rate (sum of per-base errors / sum of per-base depths).
    let read_len = seq.len();
    let per_base_arrays =
        |max_depth: i32, min_depth: i32, error_count: i32| -> (Vec<i16>, Vec<i16>) {
            let clamp = |v: i32| i16::try_from(v).unwrap_or(i16::MAX);
            let mut depths = vec![clamp(max_depth); read_len];
            if read_len >= 2 {
                // Anchor the per-base minimum at the summary min (max stays the fill).
                depths[read_len - 1] = clamp(min_depth);
            }
            let mut errors = vec![0i16; read_len];
            let n = usize::try_from(error_count).unwrap_or(0).min(read_len);
            for slot in errors.iter_mut().take(n) {
                *slot = 1;
            }
            (depths, errors)
        };
    let error_rate = |errors: &[i16], depths: &[i16]| -> f32 {
        let total_errors: i64 = errors.iter().map(|&e| i64::from(e)).sum();
        let total_depth: i64 = depths.iter().map(|&d| i64::from(d)).sum();
        if total_depth > 0 { (total_errors as f64 / total_depth as f64) as f32 } else { 0.0 }
    };

    if let Some((ad, bd, am, bm, ae, be)) = duplex_tags {
        // Duplex: emit per-strand summary (aD/bD/aM/bM int, aE/bE float rate) plus
        // the per-base strand arrays `fgumi filter` masks with (`mask_duplex_bases`
        // reads AD/AE/BD/BE only).
        //
        // The combined cD/cM/cE scalars are derived from the strand sums rather than
        // the independently sampled values, matching `duplex_caller`, which computes
        // them over the per-base combined depth (ab_i + ba_i).
        //
        // No CD_BASES/CE_BASES here: the real duplex caller never emits the combined
        // arrays, so writing them would make simulated duplex BAMs diverge from
        // anything production produces.
        let (ad_bases, ae_bases) = per_base_arrays(ad, am, ae);
        let (bd_bases, be_bases) = per_base_arrays(bd, bm, be);

        let combined_depths: Vec<i16> = ad_bases
            .iter()
            .zip(&bd_bases)
            .map(|(&a, &b)| i16::try_from(i32::from(a) + i32::from(b)).unwrap_or(i16::MAX))
            .collect();
        let combined_errors: Vec<i16> = ae_bases
            .iter()
            .zip(&be_bases)
            .map(|(&a, &b)| i16::try_from(i32::from(a) + i32::from(b)).unwrap_or(i16::MAX))
            .collect();

        // Sum the scalars rather than taking max/min over `combined_depths`: the two
        // agree whenever the arrays carry the min anchor, but `per_base_arrays` only
        // anchors it when `read_len >= 2`, so a 1-base read would report cM == cD (and
        // a 0-base read cD == cM == 0) and drift from the truth TSV again.
        b.add_int_tag(SamTag::CD, ad + bd)
            .add_int_tag(SamTag::CM, am + bm)
            .add_float_tag(SamTag::CE, error_rate(&combined_errors, &combined_depths));

        b.add_int_tag(SamTag::AD, ad)
            .add_int_tag(SamTag::BD, bd)
            .add_int_tag(SamTag::AM, am)
            .add_int_tag(SamTag::BM, bm)
            .add_float_tag(SamTag::AE, error_rate(&ae_bases, &ad_bases))
            .add_float_tag(SamTag::BE, error_rate(&be_bases, &bd_bases));
        b.add_array_i16(SamTag::AD_BASES, &ad_bases)
            .add_array_i16(SamTag::AE_BASES, &ae_bases)
            .add_array_i16(SamTag::BD_BASES, &bd_bases)
            .add_array_i16(SamTag::BE_BASES, &be_bases);
    } else {
        // Simplex: combined summary + the per-base arrays `mask_bases` reads.
        let (cd_bases, ce_bases) = per_base_arrays(cd, cm, ce);
        b.add_int_tag(SamTag::CD, cd)
            .add_int_tag(SamTag::CM, cm)
            .add_float_tag(SamTag::CE, error_rate(&ce_bases, &cd_bases));
        b.add_array_i16(SamTag::CD_BASES, &cd_bases).add_array_i16(SamTag::CE_BASES, &ce_bases);
    }

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
    use rstest::rstest;

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
        let seq = b"ACGT";
        let quals = vec![40; 4];

        // High error count (more errors than bases)
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
            100,
            None,
            None,
            MethylationMode::Disabled,
        );
        let record = to_record_buf(&raw);

        assert!(record.name().is_some());
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

    /// Build `GenerationParams` over a small synthetic reference, duplex-configurable.
    fn strand_test_params(fasta: &tempfile::NamedTempFile, duplex: bool) -> Arc<GenerationParams> {
        strand_test_params_with_read_length(fasta, duplex, 50)
    }

    /// As [`strand_test_params`], with an explicit read length so tests can reach the
    /// short-read boundaries (`--read-length 1` skips the per-base min anchor).
    fn strand_test_params_with_read_length(
        fasta: &tempfile::NamedTempFile,
        duplex: bool,
        read_length: usize,
    ) -> Arc<GenerationParams> {
        let ref_genome = Arc::new(ReferenceGenome::load(fasta.path()).unwrap());
        Arc::new(GenerationParams {
            read_length,
            min_depth: 1,
            max_depth: 30,
            depth_mean: 12.0,
            depth_stddev: 6.0,
            error_rate_mean: 0.02,
            error_rate_stddev: 0.01,
            duplex,
            consensus_quality: 40,
            methylation_mode: MethylationMode::Disabled,
            methylation_depth_dist: create_depth_distribution(5.0, 2.5),
            cpg_methylation_rate: 0.75,
            conversion_rate: 0.98,
            ref_genome,
        })
    }

    fn small_test_fasta() -> tempfile::NamedTempFile {
        use std::io::Write as IoWrite;
        let mut fasta = tempfile::NamedTempFile::new().unwrap();
        writeln!(fasta, ">chr1").unwrap();
        fasta.write_all(&b"ACGT".repeat(500)).unwrap();
        writeln!(fasta).unwrap();
        fasta.flush().unwrap();
        fasta
    }

    /// Duplex strand minimum depths must sum to the combined minimum, and each
    /// strand's minimum must not exceed its own maximum.
    ///
    /// Previously `aM` and `bM` were each computed as `strand_depth.min(cM)`, so both
    /// could equal `cM` and `aM + bM` routinely exceeded it — making the truth file
    /// and the emitted tags mutually inconsistent. The strand fraction is sampled per
    /// read, so the invariants are checked as a property over generated seeds rather
    /// than at a handful of fixed points; the reference and params are built once and
    /// shared across cases because loading them is the expensive part.
    #[test]
    fn test_duplex_strand_minimums_sum_to_combined_minimum() {
        use proptest::prelude::*;

        let fasta = small_test_fasta();
        let params = strand_test_params(&fasta, true);
        let strand_bias = StrandBiasModel::new(5.0, 5.0);

        proptest!(|(seed in any::<u64>())| {
            let pair = generate_consensus_pair(0, seed, &params, &strand_bias);
            let aux = fgumi_raw_bam::aux_data_slice(&pair.r1_record);

            let get = |tag: SamTag| -> i64 {
                fgumi_raw_bam::find_int_tag(aux, tag)
                    .unwrap_or_else(|| panic!("tag {tag:?} missing for seed {seed}"))
            };
            let (cd, cm, ad, bd, am, bm) = (
                get(SamTag::CD),
                get(SamTag::CM),
                get(SamTag::AD),
                get(SamTag::BD),
                get(SamTag::AM),
                get(SamTag::BM),
            );

            // The truth TSV and the emitted tags must not drift apart: the truth
            // tuple carries the sampled values while duplex CD/CM are derived from
            // the strand sums, and only the aM+bM==cM invariant keeps them equal.
            let (truth_cd, truth_cm) = (pair.truth.0, pair.truth.1);
            prop_assert_eq!(
                i64::from(truth_cd),
                cd,
                "truth cD ({}) disagrees with emitted CD ({}) for seed {}",
                truth_cd,
                cd,
                seed
            );
            prop_assert_eq!(
                i64::from(truth_cm),
                cm,
                "truth cM ({}) disagrees with emitted CM ({}) for seed {}",
                truth_cm,
                cm,
                seed
            );

            prop_assert_eq!(am + bm, cm, "aM + bM != cM (aM={} bM={} cM={})", am, bm, cm);
            prop_assert_eq!(ad + bd, cd, "aD + bD != cD (aD={} bD={} cD={})", ad, bd, cd);
            prop_assert!(am <= ad, "aM ({}) > aD ({})", am, ad);
            prop_assert!(bm <= bd, "bM ({}) > bD ({})", bm, bd);
            prop_assert!(am >= 0 && bm >= 0, "negative strand minimum (aM={} bM={})", am, bm);
        });
    }

    /// The duplex cD/cM tags must match the truth TSV at short read lengths too.
    ///
    /// `per_base_arrays` only anchors the per-base minimum when `read_len >= 2`, so
    /// deriving the scalars as max/min over the combined per-base depths would report
    /// `cM == cD` for a 1-base read — the same truth-vs-tag drift this change fixes.
    #[rstest]
    #[case::single_base(1)]
    #[case::two_bases(2)]
    #[case::typical(50)]
    fn test_duplex_combined_scalars_match_truth_at_short_read_lengths(#[case] read_length: usize) {
        let fasta = small_test_fasta();
        let params = strand_test_params_with_read_length(&fasta, true, read_length);
        let strand_bias = StrandBiasModel::new(5.0, 5.0);

        for seed in 0..8u64 {
            let pair = generate_consensus_pair(0, seed, &params, &strand_bias);
            let aux = fgumi_raw_bam::aux_data_slice(&pair.r1_record);
            let cd = fgumi_raw_bam::find_int_tag(aux, SamTag::CD).expect("CD missing");
            let cm = fgumi_raw_bam::find_int_tag(aux, SamTag::CM).expect("CM missing");
            let (truth_cd, truth_cm) = (pair.truth.0, pair.truth.1);
            assert_eq!(
                i64::from(truth_cd),
                cd,
                "truth cD ({truth_cd}) != emitted CD ({cd}) at read_length {read_length}, seed {seed}"
            );
            assert_eq!(
                i64::from(truth_cm),
                cm,
                "truth cM ({truth_cm}) != emitted CM ({cm}) at read_length {read_length}, seed {seed}"
            );
        }
    }

    /// Duplex records must not carry the combined per-base arrays: the real
    /// `duplex_caller` emits only the per-strand AD/AE/BD/BE arrays, so emitting
    /// `CD_BASES`/`CE_BASES` would make simulated duplex BAMs unfaithful fixtures.
    /// Simplex records must still carry them (`fgumi filter`'s `mask_bases` reads them).
    #[test]
    fn test_duplex_omits_combined_per_base_arrays_simplex_keeps_them() {
        let fasta = small_test_fasta();
        let strand_bias = StrandBiasModel::new(5.0, 5.0);

        let duplex_params = strand_test_params(&fasta, true);
        let duplex_pair = generate_consensus_pair(0, 7, &duplex_params, &strand_bias);
        let duplex_aux = fgumi_raw_bam::aux_data_slice(&duplex_pair.r1_record);
        assert!(
            fgumi_raw_bam::find_array_tag(duplex_aux, SamTag::CD_BASES).is_none(),
            "duplex record must not carry CD_BASES (no real caller emits it)"
        );
        assert!(
            fgumi_raw_bam::find_array_tag(duplex_aux, SamTag::CE_BASES).is_none(),
            "duplex record must not carry CE_BASES (no real caller emits it)"
        );
        // The per-strand arrays filter masks with must still be present.
        for tag in [SamTag::AD_BASES, SamTag::AE_BASES, SamTag::BD_BASES, SamTag::BE_BASES] {
            assert!(
                fgumi_raw_bam::find_array_tag(duplex_aux, tag).is_some(),
                "duplex record missing {tag:?}"
            );
        }

        let simplex_params = strand_test_params(&fasta, false);
        let simplex_pair = generate_consensus_pair(0, 7, &simplex_params, &strand_bias);
        let simplex_aux = fgumi_raw_bam::aux_data_slice(&simplex_pair.r1_record);
        assert!(
            fgumi_raw_bam::find_array_tag(simplex_aux, SamTag::CD_BASES).is_some(),
            "simplex record must keep CD_BASES for filter's mask_bases"
        );
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
}

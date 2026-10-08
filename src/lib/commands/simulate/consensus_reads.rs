//! Generate consensus BAM with tags for filter.

use crate::commands::command::Command;
use crate::commands::common::CompressionOptions;
use crate::commands::simulate::common::{
    MethylationArgs, MethylationConfig, ReferenceGenome, StrandBiasArgs, TemplateLocus,
    convert_molecule_strand, join_writer_result,
};
use crate::commands::simulate::region_to_bin;
use crate::dna::reverse_complement;
use crate::sam::SamTag;
use crate::simulate::{StrandBiasModel, close_output, create_rng};
use anyhow::{Context, Result};
use clap::Parser;
use crossbeam_channel::bounded;
use fgumi_bam_io::OutputFile;
use fgumi_bam_io::ProgressTracker;
use fgumi_bam_io::create_raw_bam_writer;
use fgumi_consensus::MethylationMode;
use fgumi_consensus::methylation::{ConversionPattern, MethylationAnnotation, MethylationEvidence};
use fgumi_raw_bam::{RawRecord, SamBuilder, flags as raw_flags};
use log::info;
use noodles::sam::header::Header;
use rand::{Rng, RngExt};
use rand_distr::{Distribution, LogNormal, Normal};
use rayon::prelude::*;
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
    #[arg(long = "duplex", value_name = "true|false", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = clap::builder::BoolishValueParser::new(), hide_possible_values = true)]
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
    /// Methylation model (mode `Disabled` if not set).
    methylation: MethylationConfig,
    /// Distribution for methylation count sampling.
    methylation_depth_dist: LogNormal<f64>,
    /// Loaded reference genome for sampling template sequences.
    ref_genome: Arc<ReferenceGenome>,
}

/// Channel capacity for buffering read pairs between producer and writer threads.
const CHANNEL_CAPACITY: usize = 1_000;

impl Command for ConsensusReads {
    fn execute(&self, command_line: &str) -> Result<()> {
        let methylation = self.methylation.resolve(self.seed);
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
            methylation.log_settings();
            info!("  Methylation depth mean: {}", self.methylation_depth_mean);
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
            methylation,
            methylation_depth_dist: if methylation.mode.is_enabled() {
                create_depth_distribution(
                    self.methylation_depth_mean,
                    self.methylation_depth_mean / 2.0,
                )
            } else {
                // Unused when methylation is disabled; use safe defaults
                create_depth_distribution(5.0, 2.5)
            },
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
                let truth_file = OutputFile::create(truth_path)
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

            if let (Some(tw), Some(truth_path)) = (truth_writer, truth_path.as_ref()) {
                close_output(tw, truth_path)?;
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

        let read_count = join_writer_result(writer_handle, generation_result)?;

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

    // Coin flip for strand orientation: R1 reads the molecule's top strand forward, or its
    // bottom strand reverse.
    let is_top_strand: bool = rng.random();
    let r1_is_reverse = !is_top_strand;

    // Methylation: one conversion per original strand of the molecule, as the consensus
    // callers see it (see `simulate_molecule_methylation`).
    let molecule = params.methylation.mode.is_enabled().then(|| {
        let locus = TemplateLocus {
            chrom_idx: Some(chrom_idx),
            contig: params.ref_genome.contig(chrom_idx),
            start: local_pos,
        };
        simulate_molecule_methylation(
            &seq,
            locus,
            is_top_strand,
            params.duplex,
            &params.methylation,
            &params.methylation_depth_dist,
            &mut rng,
        )
    });
    let record_bases = |is_r2: bool, is_reverse: bool| match &molecule {
        Some(molecule) => {
            let (bases, methylation) = molecule.record(is_r2, is_reverse);
            (bases, Some(methylation))
        }
        None if is_reverse => (reverse_complement(&seq), None),
        None => (seq.clone(), None),
    };

    // Build R1 record
    let (r1_seq, r1_methylation) = record_bases(false, r1_is_reverse);
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
        r1_methylation.as_ref(),
        params.methylation.mode,
    );

    // Build R2 record (opposite strand orientation)
    let (r2_seq, r2_methylation) = record_bases(true, !r1_is_reverse);
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
        r2_methylation.as_ref(),
        params.methylation.mode,
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

/// Methylation of one consensus record, in its read orientation.
struct MethylationData {
    /// Evidence of the strand family of the record's own read type (the AB strand of a
    /// duplex), informative at `C` for R1 and `G` for R2.
    ab: MethylationAnnotation,
    /// Duplex only: evidence of the other strand (BA), at the complementary base.
    ba: Option<MethylationAnnotation>,
}

/// Methylation of one consensus molecule over the read span, in genomic orientation.
struct MoleculeMethylation {
    /// The span as the records show it: each strand's cytosines after its conversion.
    converted: Vec<u8>,
    /// Evidence of the strand R1 reads (AB), at its cytosines.
    ab: Vec<MethylationEvidence>,
    /// Duplex only: evidence of the other strand (BA), at its cytosines.
    ba: Option<Vec<MethylationEvidence>>,
}

/// Simulates the methylation evidence a consensus caller would record for one molecule.
///
/// Each original strand is converted once (`convert_molecule_strand`): every read of its
/// family is a PCR copy of it, so at each of that strand's cytosines all `depth` reads show
/// the same base and the counts are either all unconverted or all converted. The top strand's
/// cytosines are reference `C`, the bottom strand's reference `G` (genomic orientation). A
/// simplex consensus has the strand R1 reads (`r1_strand_is_top`); a duplex consensus has both.
fn simulate_molecule_methylation(
    span: &[u8],
    locus: TemplateLocus<'_>,
    r1_strand_is_top: bool,
    duplex: bool,
    config: &MethylationConfig,
    depth_dist: &LogNormal<f64>,
    rng: &mut impl Rng,
) -> MoleculeMethylation {
    let mut converted = span.to_vec();
    let mut strand_evidence = |is_top: bool, rng: &mut _| -> Vec<MethylationEvidence> {
        let strand_bases = convert_molecule_strand(span, locus, is_top, config, rng);
        let cytosine = if is_top { b'C' } else { b'G' };
        span.iter()
            .zip(&strand_bases)
            .enumerate()
            .map(|(i, (&reference, &base))| {
                if reference != cytosine {
                    return MethylationEvidence::default();
                }
                converted[i] = base;
                let depth = sample_depth(depth_dist, 1, 100, rng) as u32;
                let (unconverted_count, converted_count) =
                    if base == reference { (depth, 0) } else { (0, depth) };
                MethylationEvidence { informative: true, unconverted_count, converted_count }
            })
            .collect()
    };
    let ab = strand_evidence(r1_strand_is_top, rng);
    let ba = duplex.then(|| strand_evidence(!r1_strand_is_top, rng));
    MoleculeMethylation { converted, ab, ba }
}

impl MoleculeMethylation {
    /// The converted bases and methylation of one record of this molecule, in its read
    /// orientation.
    ///
    /// The AB strand's cytosines appear as `C` in R1's read orientation and `G` in R2's (R2
    /// is the copy complementary to the original strand); a duplex record's BA strand the
    /// other way round. A simplex record keeps these bases; a duplex record's SEQ is restored
    /// to the molecule's sequence when it is written (as the duplex caller writes it).
    fn record(&self, is_r2: bool, is_reverse: bool) -> (Vec<u8>, MethylationData) {
        let (own, other) = if is_r2 {
            (ConversionPattern::GToA, ConversionPattern::CToT)
        } else {
            (ConversionPattern::CToT, ConversionPattern::GToA)
        };
        let orient = |evidence: &[MethylationEvidence], pattern| {
            let mut evidence = evidence.to_vec();
            if is_reverse {
                evidence.reverse();
            }
            MethylationAnnotation { evidence, pattern }
        };
        let ab = orient(&self.ab, own);
        let ba = self.ba.as_deref().map(|ba| orient(ba, other));
        let bases =
            if is_reverse { reverse_complement(&self.converted) } else { self.converted.clone() };
        (bases, MethylationData { ab, ba })
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

    // A duplex methylation record is written as the duplex caller writes it, by the same code:
    // SEQ restored to the molecule's sequence, with MM/ML/MN, am/bm and the counts.
    let duplex_methylation = methylation.and_then(|meth| {
        meth.ba.as_ref().map(|ba| {
            fgumi_consensus::methylation::duplex_methylation_tags(
                seq,
                &[&meth.ab, ba],
                methylation_mode,
            )
        })
    });
    let seq: &[u8] = duplex_methylation.as_ref().map_or(seq, |tags| &tags.seq);

    // `seq`/`quals` are in read orientation; BAM stores SEQ/QUAL of a reverse-strand record in
    // reference orientation.
    let (stored_seq, stored_quals) =
        super::common::to_reference_orientation(seq, quals, is_reverse);
    // Per-base arrays are stored co-oriented with SEQ (reversed on reverse-strand records), the
    // orientation `fgumi filter` reads by default: a simulated BAM is filtered without
    // `--reverse-per-base-tags`. MM/ML index the sequenced read, so they stay in read orientation.
    let oriented = |values: &[i16]| -> Vec<i16> {
        if is_reverse { values.iter().rev().copied().collect() } else { values.to_vec() }
    };

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
        .sequence(&stored_seq)
        .qualities(&stored_quals);

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
        b.add_array_i16(SamTag::AD_BASES, &oriented(&ad_bases))
            .add_array_i16(SamTag::AE_BASES, &oriented(&ae_bases))
            .add_array_i16(SamTag::BD_BASES, &oriented(&bd_bases))
            .add_array_i16(SamTag::BE_BASES, &oriented(&be_bases));
    } else {
        // Simplex: combined summary + the per-base arrays `mask_bases` reads.
        let (cd_bases, ce_bases) = per_base_arrays(cd, cm, ce);
        b.add_int_tag(SamTag::CD, cd)
            .add_int_tag(SamTag::CM, cm)
            .add_float_tag(SamTag::CE, error_rate(&ce_bases, &cd_bases));
        b.add_array_i16(SamTag::CD_BASES, &oriented(&cd_bases))
            .add_array_i16(SamTag::CE_BASES, &oriented(&ce_bases));
    }

    // Methylation tags, as the consensus callers write them. Simplex: the counts only. Duplex:
    // the per-strand calls and counts, the combined MM/ML/MN and counts (see above).
    if let Some(tags) = &duplex_methylation {
        for ((mm_tag, u_tag, t_tag), strand) in
            [(SamTag::AM_BASES, SamTag::AU, SamTag::AT), (SamTag::BM_BASES, SamTag::BU, SamTag::BT)]
                .into_iter()
                .zip(&tags.strands)
        {
            if let Some(mm) = &strand.mm {
                b.add_string_tag(mm_tag, mm.as_bytes());
            }
            b.add_array_i16(u_tag, &oriented(&strand.unconverted))
                .add_array_i16(t_tag, &oriented(&strand.converted));
        }
        if let Some((mm, ml)) = &tags.mm_ml {
            b.add_string_tag(SamTag::MM, mm.as_bytes())
                .add_array_u8(SamTag::ML, ml)
                .add_int_tag(SamTag::MN, i32::try_from(seq.len()).unwrap_or(i32::MAX));
        }
        b.add_array_i16(SamTag::CU, &oriented(&tags.unconverted))
            .add_array_i16(SamTag::CT, &oriented(&tags.converted));
    } else if let Some(meth) = methylation {
        b.add_array_i16(SamTag::CU, &oriented(&meth.ab.unconverted_counts()))
            .add_array_i16(SamTag::CT, &oriented(&meth.ab.converted_counts()));
    }

    b.build()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::commands::simulate::common::generate_random_sequence;
    use crate::simulate::create_rng;
    use fgumi_raw_bam::RawRecordView;
    use rstest::rstest;

    /// Decode a raw BAM record into a noodles `RecordBuf` for higher-level test
    /// assertions.
    fn to_record_buf(raw: &RawRecord) -> noodles::sam::alignment::RecordBuf {
        fgumi_raw_bam::raw_record_to_record_buf(raw, &noodles::sam::Header::default())
            .expect("raw_record_to_record_buf failed in test")
    }

    /// A record's methylation-relevant content, in read orientation.
    struct ReadView {
        is_r2: bool,
        seq: Vec<u8>,
        reference: Vec<u8>,
        counts: [Option<Vec<i16>>; 6], // cu, ct, au, at, bu, bt
        mm: Option<String>,
        ml: Option<Vec<u8>>,
        mn: Option<i64>,
        am: Option<String>,
        bm: Option<String>,
    }

    /// Decodes `raw` against `contig` into read orientation. A reverse record stores SEQ and
    /// the per-base count arrays in reference orientation, so both are flipped back here.
    fn read_view(raw: &RawRecord, contig: &[u8]) -> ReadView {
        let view = RawRecordView::new(raw.as_ref());
        let stored = view.sequence_vec();
        let start = usize::try_from(view.pos()).unwrap();
        let genomic_reference = &contig[start..start + stored.len()];
        let is_reverse = view.is_reverse();
        let (seq, reference) = if is_reverse {
            (reverse_complement(&stored), reverse_complement(genomic_reference))
        } else {
            (stored, genomic_reference.to_vec())
        };
        let aux = fgumi_raw_bam::aux_data_slice(raw.as_ref());
        let array = |tag: SamTag| {
            fgumi_raw_bam::find_array_tag(aux, tag).map(|a| {
                let mut values: Vec<i16> =
                    fgumi_raw_bam::array_tag_to_vec_u16(&a).into_iter().map(|v| v as i16).collect();
                if is_reverse {
                    values.reverse();
                }
                values
            })
        };
        let string = |tag: SamTag| {
            fgumi_raw_bam::find_string_tag(aux, tag).map(|v| String::from_utf8(v.to_vec()).unwrap())
        };
        ReadView {
            is_r2: view.flags() & raw_flags::LAST_SEGMENT != 0,
            seq,
            reference,
            counts: [
                array(SamTag::CU),
                array(SamTag::CT),
                array(SamTag::AU),
                array(SamTag::AT),
                array(SamTag::BU),
                array(SamTag::BT),
            ],
            mm: string(SamTag::MM),
            ml: fgumi_raw_bam::find_array_tag(aux, SamTag::ML).map(|a| a.data.to_vec()),
            mn: fgumi_raw_bam::find_int_tag(aux, SamTag::MN),
            am: string(SamTag::AM_BASES),
            bm: string(SamTag::BM_BASES),
        }
    }

    /// EM-seq parameters over a random 4 kb contig; returns the params and the contig.
    fn methylation_params(
        mode: MethylationMode,
        duplex: bool,
    ) -> (Arc<GenerationParams>, Vec<u8>, tempfile::NamedTempFile) {
        use std::io::Write as IoWrite;
        let mut contig_rng = create_rng(Some(99));
        let contig: Vec<u8> = (0..4000).map(|_| b"ACGT"[contig_rng.random_range(0..4)]).collect();
        let mut fasta = tempfile::NamedTempFile::new().unwrap();
        writeln!(fasta, ">chr1").unwrap();
        fasta.write_all(&contig).unwrap();
        writeln!(fasta).unwrap();
        fasta.flush().unwrap();
        let params = Arc::new(GenerationParams {
            read_length: 60,
            min_depth: 1,
            max_depth: 10,
            depth_mean: 5.0,
            depth_stddev: 2.0,
            error_rate_mean: 0.01,
            error_rate_stddev: 0.005,
            duplex,
            consensus_quality: 40,
            methylation: MethylationConfig {
                mode,
                cpg_methylation_rate: 0.5,
                conversion_rate: 1.0,
                hemimethylation_rate: 0.0,
                failure_rate: 0.0,
                table_seed: 3,
            },
            methylation_depth_dist: create_depth_distribution(5.0, 2.5),
            ref_genome: Arc::new(ReferenceGenome::load(fasta.path()).unwrap()),
        });
        (params, contig, fasta)
    }

    /// Checks one record (read orientation) against what the consensus callers emit.
    #[allow(clippy::too_many_lines)]
    fn assert_matches_caller_output(read: &ReadView, duplex: bool, mode: MethylationMode) {
        let len = read.seq.len();
        let [cu, ct, au, at, bu, bt] = &read.counts;
        let (cu, ct) = (cu.as_ref().expect("cu"), ct.as_ref().expect("ct"));
        // The record's own read type is informative at C (R1) / G (R2) in read orientation;
        // a duplex record also carries the other strand at the complementary base.
        let (own, other) = if read.is_r2 { (b'G', b'C') } else { (b'C', b'G') };
        let converted_base = |base: u8| if base == b'C' { b'T' } else { b'A' };
        for i in 0..len {
            let reference = read.reference[i];
            let informative = reference == own || (duplex && reference == other);
            let depth = cu[i] + ct[i];
            assert_eq!(depth > 0, informative, "position {i}: evidence only at informative bases");
            if informative {
                assert!(cu[i] == 0 || ct[i] == 0, "position {i}: a family is one converted strand");
            }
            // Simplex keeps the converted bases; duplex is the molecule's sequence.
            let expected_base =
                if !duplex && ct[i] > 0 { converted_base(reference) } else { reference };
            assert_eq!(read.seq[i], expected_base, "position {i}: SEQ");
        }
        if duplex {
            let (au, at, bu, bt) = (
                au.as_ref().unwrap(),
                at.as_ref().unwrap(),
                bu.as_ref().unwrap(),
                bt.as_ref().unwrap(),
            );
            for i in 0..len {
                assert_eq!(
                    au[i] + at[i] > 0,
                    read.reference[i] == own,
                    "position {i}: AB at own class"
                );
                assert_eq!(
                    bu[i] + bt[i] > 0,
                    read.reference[i] == other,
                    "position {i}: BA at other class"
                );
                assert_eq!(
                    (cu[i], ct[i]),
                    (au[i] + bu[i], at[i] + bt[i]),
                    "position {i}: cu/ct sums"
                );
            }
        } else {
            assert!(au.is_none() && bu.is_none(), "simplex records carry no per-strand counts");
        }
        if duplex {
            let mm = read.mm.as_deref().expect("duplex MM");
            assert_eq!(read.mn, Some(i64::try_from(len).unwrap()), "MN");
            assert!(mm.starts_with("C+m?") && mm.contains(";G-m?"), "duplex MM: {mm}");
            assert_mm_ml_match_counts(read, mm, mode);
            // am/bm: one group per strand, by the read type of that strand's reads here (AB is
            // the record's own type), listing exactly the strand's sites with evidence.
            let [_, _, au, at, bu, bt] = &read.counts;
            let strand_sites = |u: &Option<Vec<i16>>, t: &Option<Vec<i16>>| -> Vec<usize> {
                let (u, t) = (u.as_ref().unwrap(), t.as_ref().unwrap());
                (0..len).filter(|&i| u[i] + t[i] > 0).collect()
            };
            let (am_header, bm_header) =
                if read.is_r2 { ("G-m?", "C+m?") } else { ("C+m?", "G-m?") };
            for (tag, value, header, sites) in [
                ("am", &read.am, am_header, strand_sites(au, at)),
                ("bm", &read.bm, bm_header, strand_sites(bu, bt)),
            ] {
                let value = value.as_deref().unwrap_or_else(|| panic!("duplex {tag}"));
                assert_eq!(
                    decode_single_group(value, &read.seq),
                    (header.to_string(), sites),
                    "{tag}: {value}"
                );
            }
        } else {
            assert!(read.mm.is_none() && read.mn.is_none(), "no MM/MN on simplex");
            assert!(read.am.is_none() && read.bm.is_none(), "no am/bm on simplex");
        }
    }

    /// Decodes a one-group MM-format string against SEQ into its header and listed positions.
    fn decode_single_group(mm: &str, seq: &[u8]) -> (String, Vec<usize>) {
        let group = mm.strip_suffix(';').expect("terminated group");
        assert!(!group.contains(';'), "one group: {mm}");
        let mut fields = group.split(',');
        let header = fields.next().unwrap().to_string();
        let base = header.as_bytes()[0];
        let tracked: Vec<usize> = (0..seq.len()).filter(|&i| seq[i] == base).collect();
        let mut ordinal = 0;
        let positions = fields
            .map(|skip| {
                ordinal += skip.parse::<usize>().unwrap();
                let pos = tracked[ordinal];
                ordinal += 1;
                pos
            })
            .collect();
        (header, positions)
    }

    /// Decodes MM against SEQ and checks it lists exactly the tracked bases with evidence, each
    /// with the probability the counts imply: a family is one converted strand, so a site is
    /// all unconverted (EM-seq 255, TAPs 0) or all converted (EM-seq 0, TAPs 255).
    fn assert_mm_ml_match_counts(read: &ReadView, mm: &str, mode: MethylationMode) {
        let [cu, ct, ..] = &read.counts;
        let (cu, ct) = (cu.as_ref().unwrap(), ct.as_ref().unwrap());
        let ml = read.ml.as_ref().expect("duplex ML");
        let mut ml_values = ml.iter();
        for group in mm.trim_end_matches(';').split(';') {
            let mut fields = group.split(',');
            let base = fields.next().unwrap().as_bytes()[0];
            let tracked: Vec<usize> =
                (0..read.seq.len()).filter(|&i| read.seq[i] == base).collect();
            let mut listed = Vec::new();
            let mut ordinal = 0usize;
            for skip in fields {
                ordinal += skip.parse::<usize>().unwrap();
                listed.push(tracked[ordinal]);
                ordinal += 1;
            }
            let with_evidence: Vec<usize> =
                tracked.iter().copied().filter(|&i| cu[i] + ct[i] > 0).collect();
            assert_eq!(listed, with_evidence, "MM lists exactly the {} sites with evidence", base);
            for &i in &listed {
                let methylated = match mode {
                    MethylationMode::Taps => ct[i] > 0,
                    _ => cu[i] > 0,
                };
                let expected = if methylated { 255 } else { 0 };
                assert_eq!(ml_values.next(), Some(&expected), "ML at {i}");
            }
        }
        assert!(ml_values.next().is_none(), "ML has no extra values");
    }

    /// Simulated consensus records carry what `simplex`/`duplex` emit after this model: simplex
    /// SEQ with converted bases and no MM/ML, duplex SEQ as the molecule's sequence with
    /// MM/ML/MN, evidence only at the read type's informative base, family-coherent counts,
    /// and, for duplex, AB and BA counts at disjoint, complementary bases.
    #[rstest]
    #[case::simplex(MethylationMode::EmSeq, false)]
    #[case::duplex(MethylationMode::EmSeq, true)]
    #[case::taps_simplex(MethylationMode::Taps, false)]
    #[case::taps_duplex(MethylationMode::Taps, true)]
    fn test_methylation_records_match_caller_output(
        #[case] mode: MethylationMode,
        #[case] duplex: bool,
    ) {
        let (params, contig, _fasta) = methylation_params(mode, duplex);
        let strand_bias = StrandBiasModel::new(5.0, 5.0);
        let mut orientations = std::collections::BTreeSet::new();
        for seed in 0..20u64 {
            let pair = generate_consensus_pair(0, seed, &params, &strand_bias);
            let r1 = read_view(&pair.r1_record, &contig);
            let r2 = read_view(&pair.r2_record, &contig);
            assert_matches_caller_output(&r1, duplex, mode);
            assert_matches_caller_output(&r2, duplex, mode);
            let mut r1_cu = r1.counts[0].clone().unwrap();
            r1_cu.reverse();
            assert_eq!(Some(r1_cu), r2.counts[0], "both mates read one molecule: mirrored counts");
            orientations.insert(pair.is_top_strand);
        }
        assert_eq!(orientations.len(), 2, "both strands exercised");
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
                hemimethylation_rate: 0.0,
                failure_rate: 0.0,
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
                hemimethylation_rate: 0.0,
                failure_rate: 0.0,
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
                hemimethylation_rate: 0.0,
                failure_rate: 0.0,
            },
            methylation_depth_mean: f64::NAN,
            reference: PathBuf::from("dummy.fa"),
        };
        let result = cmd.execute("test");
        assert!(result.is_err());
    }

    /// End to end, MM/ML follow the caller being simulated: none on simplex records (converted
    /// SEQ), present on duplex records (the molecule's sequence).
    #[rstest]
    #[case::simplex(false)]
    #[case::duplex(true)]
    fn test_methylation_mm_follows_the_caller(#[case] duplex: bool) {
        use std::io::Write as IoWrite;
        let mut fasta = tempfile::NamedTempFile::new().unwrap();
        writeln!(fasta, ">chr1").unwrap();
        let mut contig_rng = create_rng(Some(7));
        let contig: Vec<u8> = (0..2000).map(|_| b"ACGT"[contig_rng.random_range(0..4)]).collect();
        fasta.write_all(&contig).unwrap();
        writeln!(fasta).unwrap();
        fasta.flush().unwrap();
        let dir = tempfile::TempDir::new().unwrap();
        let output = dir.path().join("out.bam");

        let args = vec![
            "consensus-reads".to_string(),
            "-o".to_string(),
            output.display().to_string(),
            "-r".to_string(),
            fasta.path().display().to_string(),
            "-n".to_string(),
            "5".to_string(),
            "-l".to_string(),
            "40".to_string(),
            "--seed".to_string(),
            "1".to_string(),
            "--methylation-mode".to_string(),
            "em-seq".to_string(),
            format!("--duplex={duplex}"),
        ];
        let cmd = ConsensusReads::try_parse_from(&args).expect("valid args");
        cmd.execute("test").expect("simulation succeeds");
        let mut reader = noodles::bam::io::reader::Builder.build_from_path(&output).unwrap();
        let header = reader.read_header().unwrap();
        let mm_tag = SamTag::MM.to_noodles_tag();
        let mut records = 0;
        for record in reader.records() {
            let record = record.unwrap();
            let buf =
                noodles::sam::alignment::RecordBuf::try_from_alignment_record(&header, &record)
                    .unwrap();
            assert_eq!(buf.data().get(&mm_tag).is_some(), duplex);
            records += 1;
        }
        assert_eq!(records, 10, "five pairs written");
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
                hemimethylation_rate: 0.0,
                failure_rate: 0.0,
            },
            methylation_depth_mean: 0.0,
            reference: fasta.path().to_path_buf(),
        };
        // Should succeed since methylation is disabled
        let result = cmd.execute("test");
        assert!(result.is_ok(), "disabled methylation should not validate depth mean");
    }

    /// Per-base arrays of a reverse-strand record are stored in reference orientation, co-
    /// oriented with SEQ (as `fgumi filter` expects by default, without
    /// `--reverse-per-base-tags`): each array is the reverse of the same record built forward.
    #[rstest::rstest]
    #[case::simplex(None)]
    #[case::duplex(Some((6, 4, 3, 2, 2, 1)))]
    fn test_reverse_record_per_base_arrays_follow_seq(
        #[case] duplex_tags: Option<(i32, i32, i32, i32, i32, i32)>,
    ) {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value;
        use noodles::sam::alignment::record_buf::data::field::value::Array;

        let seq = b"CGATCA";
        let quals = vec![40; seq.len()];
        let evidence = |u: u32, t: u32| MethylationEvidence {
            informative: u + t > 0,
            unconverted_count: u,
            converted_count: t,
        };
        let annotation = MethylationAnnotation {
            evidence: vec![
                evidence(9, 1),
                evidence(0, 0),
                evidence(0, 0),
                evidence(0, 0),
                evidence(2, 7),
                evidence(0, 0),
            ],
            pattern: fgumi_consensus::methylation::ConversionPattern::CToT,
        };
        let meth = MethylationData {
            ab: annotation.clone(),
            ba: duplex_tags.map(|_| MethylationAnnotation {
                pattern: fgumi_consensus::methylation::ConversionPattern::GToA,
                ..annotation.clone()
            }),
        };
        let build = |is_reverse: bool| {
            to_record_buf(&build_consensus_record(
                "test/1",
                seq,
                &quals,
                true,
                is_reverse,
                0,
                100,
                5,
                3,
                2,
                duplex_tags,
                Some(&meth),
                MethylationMode::EmSeq,
            ))
        };
        let (forward, reverse) = (build(false), build(true));

        let per_base = |r: &noodles::sam::alignment::RecordBuf, tag: SamTag| match r
            .data()
            .get(&Tag::from(*tag))
        {
            Some(Value::Array(Array::Int16(arr))) => Some(arr.clone()),
            None => None,
            other => panic!("expected an i16 array for {tag:?}, got {other:?}"),
        };
        let tags: &[SamTag] = if duplex_tags.is_some() {
            &[
                SamTag::AD_BASES,
                SamTag::AE_BASES,
                SamTag::BD_BASES,
                SamTag::BE_BASES,
                SamTag::CU,
                SamTag::CT,
                SamTag::AU,
                SamTag::AT,
                SamTag::BU,
                SamTag::BT,
            ]
        } else {
            &[SamTag::CD_BASES, SamTag::CE_BASES, SamTag::CU, SamTag::CT]
        };
        for &tag in tags {
            let mut expected = per_base(&forward, tag).unwrap_or_else(|| panic!("{tag:?} set"));
            expected.reverse();
            assert_eq!(per_base(&reverse, tag), Some(expected), "{tag:?}");
        }
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
            methylation: MethylationConfig::default(),
            methylation_depth_dist: create_depth_distribution(5.0, 2.5),
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

    /// Every consensus record's SEQ is in reference orientation, so it matches the reference
    /// at its position on both strands, for simplex and duplex consensus.
    #[rstest::rstest]
    #[case::simplex(false)]
    #[case::duplex(true)]
    fn test_records_seq_matches_reference_on_both_strands(#[case] duplex: bool) {
        use crate::commands::simulate::common::test_support::{
            assert_seq_matches_reference, asymmetric_reference,
        };
        let (fasta, reference) = asymmetric_reference();
        let params = strand_test_params(&fasta, duplex);
        let strand_bias = StrandBiasModel::new(5.0, 5.0);

        let mut strands_seen = [false; 2];
        for seed in 0..16u64 {
            let pair = generate_consensus_pair(0, seed, &params, &strand_bias);
            for raw in [&pair.r1_record, &pair.r2_record] {
                let (reverse, _) = assert_seq_matches_reference(&to_record_buf(raw), &reference);
                strands_seen[usize::from(reverse)] = true;
            }
        }
        assert_eq!(strands_seen, [true, true], "both strands must be exercised");
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
            methylation: MethylationConfig::default(),
            methylation_depth_dist: create_depth_distribution(5.0, 2.5),
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

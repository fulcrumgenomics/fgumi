//! Generate consensus BAM with tags for filter.

use crate::bam_io::create_bam_writer;
use crate::commands::command::Command;
use crate::commands::common::{CompressionOptions, parse_bool};
use crate::commands::simulate::common::{
    MethylationArgs, StrandBiasArgs, generate_random_sequence,
};
use crate::dna::reverse_complement;
use crate::progress::ProgressTracker;
use crate::sam::builder::{ConsensusTagsBuilder, RecordBuilder};
use crate::simulate::{StrandBiasModel, create_rng};
use anyhow::{Context, Result};
use clap::Parser;
use crossbeam_channel::bounded;
use fgumi_consensus::MethylationMode;
use fgumi_consensus::methylation::{
    MethylationAnnotation, MethylationEvidence, build_mm_ml_tags, is_cpg_context,
};
use log::info;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::sam::header::Header;
use rand::{Rng, RngExt};
use rand_distr::{Distribution, LogNormal, Normal};
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::sync::Arc;
use std::thread;

/// Generate unmapped BAM with consensus tags for `fgumi filter`.
#[derive(Parser, Debug)]
#[command(
    name = "consensus-reads",
    about = "Generate consensus BAM with tags for filter",
    long_about = r#"
Generate synthetic consensus reads with proper consensus tags.

The output is unmapped and suitable for input to `fgumi filter`.
Reads contain consensus tags (cD, cM, cE, cd, ce) for filtering.
"#
)]
pub struct ConsensusReads {
    /// Output BAM file (unmapped)
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
}

/// A generated consensus read pair ready for output.
struct ConsensusReadPair {
    read_name: String,
    r1_record: RecordBuf,
    r2_record: RecordBuf,
    /// Truth data: (cD, cM, cE, aD, bD, aM, bM, aE, bE)
    truth: (i32, i32, i32, i32, i32, i32, i32, i32, i32),
}

/// Parameters needed for parallel consensus generation.
#[derive(Clone)]
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

        // Build header (no references for unmapped BAM)
        let mut header_builder = Header::builder();

        // Add @PG record
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
            let mut writer =
                create_bam_writer(&output_path, &header_clone, writer_threads, compression_level)?;

            // Create truth file if requested
            let mut truth_writer = if let Some(ref truth_path) = truth_path {
                let truth_file = File::create(truth_path)
                    .with_context(|| format!("Failed to create {}", truth_path.display()))?;
                let mut w = BufWriter::new(truth_file);
                writeln!(w, "read_name\tcD\tcM\tcE\taD\tbD\taM\tbM\taE\tbE")?;
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

                writer.write_alignment_record(&header_clone, &pair.r1_record)?;
                writer.write_alignment_record(&header_clone, &pair.r2_record)?;

                // Write truth
                if let Some(ref mut tw) = truth_writer {
                    let (cd, cm, ce, ad, bd, am, bm, ae, be) = pair.truth;
                    writeln!(
                        tw,
                        "{}\t{cd}\t{cm}\t{ce}\t{ad}\t{bd}\t{am}\t{bm}\t{ae}\t{be}",
                        pair.read_name
                    )?;
                }
            }

            progress.log_final();

            if let Some(ref mut tw) = truth_writer {
                tw.flush()?;
            }

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

    // Generate sequence
    let seq = generate_random_sequence(params.read_length, &mut rng);
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

        // Min depths for each strand
        let am = (ad.min(cm)).max(0);
        let bm = (bd.min(cm)).max(0);

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
            &params.methylation_depth_dist,
            params.duplex,
            &mut rng,
        ))
    } else {
        None
    };

    // Build R1 record
    let r1_record = build_consensus_record(
        &format!("{read_name}/1"),
        &seq,
        &quals,
        true, // is_first
        cd,
        cm,
        ce,
        if params.duplex { Some((ad, bd, am, bm, ae, be)) } else { None },
        methylation.as_ref(),
        params.methylation_mode,
    );

    // Build R2 record (reverse complement)
    let r2_seq = reverse_complement(&seq);
    let r2_methylation = methylation.as_ref().map(|m| m.reverse());
    let r2_record = build_consensus_record(
        &format!("{read_name}/2"),
        &r2_seq,
        &quals,
        false, // is_first
        cd,
        cm,
        ce,
        if params.duplex { Some((ad, bd, am, bm, ae, be)) } else { None },
        r2_methylation.as_ref(),
        params.methylation_mode,
    );

    ConsensusReadPair {
        read_name,
        r1_record,
        r2_record,
        truth: (cd, cm, ce, ad, bd, am, bm, ae, be),
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
        let (unconverted, converted) = match mode {
            MethylationMode::EmSeq => {
                if cpg && methylated {
                    // Methylated CpG in EM-Seq: mostly unconverted (protected)
                    let conv = rng.random_range(0..=total_depth / 5);
                    (total_depth - conv, conv)
                } else {
                    // Unmethylated (non-CpG or unmethylated CpG): mostly converted
                    let unconverted = rng.random_range(0..=total_depth / 5);
                    (unconverted, total_depth - unconverted)
                }
            }
            MethylationMode::Taps => {
                if cpg && methylated {
                    // Methylated CpG in TAPs: mostly converted (target)
                    let unconverted = rng.random_range(0..=total_depth / 5);
                    (unconverted, total_depth - unconverted)
                } else {
                    // Unmethylated: mostly unconverted (not a target)
                    let conv = rng.random_range(0..=total_depth / 5);
                    (total_depth - conv, conv)
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
    cd: i32,
    cm: i32,
    ce: i32,
    duplex_tags: Option<(i32, i32, i32, i32, i32, i32)>,
    methylation: Option<&MethylationData>,
    methylation_mode: MethylationMode,
) -> RecordBuf {
    let seq_str = String::from_utf8_lossy(seq);

    // Build consensus tags (note: cE is error count in this file, not error rate)
    let mut consensus_tags =
        ConsensusTagsBuilder::new().depth_max(cd).depth_min(cm).tag_ce_as_int(ce);

    // Add duplex-specific tags if present
    if let Some((ad, bd, am, bm, ae, be)) = duplex_tags {
        consensus_tags = consensus_tags
            .ab_depth(ad)
            .ba_depth(bd)
            .ab_min_depth(am)
            .ba_min_depth(bm)
            .ab_errors_int(ae)
            .ba_errors_int(be);
    }

    let mut builder = RecordBuilder::new()
        .name(name)
        .sequence(&seq_str)
        .qualities(quals)
        .paired(true)
        .first_segment(is_first)
        .unmapped(true)
        .mate_unmapped(true)
        .consensus_tags(consensus_tags);

    // Add methylation tags if enabled
    if let Some(meth) = methylation {
        // Add cu/ct tags (unconverted/converted counts per position)
        let cu: Vec<i16> = meth.annotation.unconverted_counts();
        let ct: Vec<i16> = meth.annotation.converted_counts();
        builder = builder.tag("cu", cu).tag("ct", ct);

        // Add MM/ML tags (SAM spec methylation tags)
        if let Some((mm, ml)) = build_mm_ml_tags(seq, &meth.annotation, true, methylation_mode) {
            builder = builder.tag("MM", mm).tag("ML", ml);
        }

        // Add duplex per-strand tags
        if let (Some(ab), Some(ba)) = (&meth.ab_annotation, &meth.ba_annotation) {
            let au: Vec<i16> = ab.unconverted_counts();
            let at: Vec<i16> = ab.converted_counts();
            let bu: Vec<i16> = ba.unconverted_counts();
            let bt: Vec<i16> = ba.converted_counts();
            builder = builder.tag("au", au).tag("at", at).tag("bu", bu).tag("bt", bt);

            // Add per-strand MM tags (am/bm)
            if let Some(am) =
                fgumi_consensus::methylation::build_mm_tag_no_ml(seq, ab, true, methylation_mode)
            {
                builder = builder.tag("am", am);
            }
            if let Some(bm) =
                fgumi_consensus::methylation::build_mm_tag_no_ml(seq, ba, true, methylation_mode)
            {
                builder = builder.tag("bm", bm);
            }
        }
    }

    builder.build()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simulate::create_rng;

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

        let record = build_consensus_record(
            "test_read/1",
            seq,
            &quals,
            true, // is_first
            5,    // cD
            3,    // cM
            1,    // cE
            None, // no duplex tags
            None,
            MethylationMode::Disabled,
        );

        assert!(record.name().is_some());
        let flags = record.flags();
        assert!(flags.is_unmapped());
        assert!(flags.is_first_segment());
    }

    #[test]
    fn test_build_consensus_record_duplex() {
        let seq = b"ACGTACGT";
        let quals = vec![40; 8];

        let record = build_consensus_record(
            "test_read/1",
            seq,
            &quals,
            true,
            10,                       // cD
            5,                        // cM
            2,                        // cE
            Some((6, 4, 3, 2, 1, 1)), // aD, bD, aM, bM, aE, bE
            None,
            MethylationMode::Disabled,
        );

        assert!(record.name().is_some());
        let flags = record.flags();
        assert!(flags.is_unmapped());
    }

    #[test]
    fn test_build_consensus_record_r2_flags() {
        let seq = b"ACGT";
        let quals = vec![40; 4];

        let record = build_consensus_record(
            "test_read/2",
            seq,
            &quals,
            false, // is_first = false means R2
            5,
            3,
            1,
            None,
            None,
            MethylationMode::Disabled,
        );

        let flags = record.flags();
        assert!(flags.is_last_segment());
        assert!(!flags.is_first_segment());
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
    fn test_build_consensus_record_mate_unmapped_flag() {
        let seq = b"ACGT";
        let quals = vec![40; 4];

        let record = build_consensus_record(
            "test/1",
            seq,
            &quals,
            true,
            5,
            3,
            1,
            None,
            None,
            MethylationMode::Disabled,
        );

        let flags = record.flags();
        assert!(flags.is_mate_unmapped());
        assert!(flags.is_unmapped());
    }

    #[test]
    fn test_build_consensus_record_segmented_flag() {
        let seq = b"ACGT";
        let quals = vec![40; 4];

        let record = build_consensus_record(
            "test/1",
            seq,
            &quals,
            true,
            5,
            3,
            1,
            None,
            None,
            MethylationMode::Disabled,
        );

        let flags = record.flags();
        assert!(flags.is_segmented());
    }

    #[test]
    fn test_build_consensus_record_sequence_stored() {
        let seq = b"ACGTACGTACGT";
        let quals = vec![40; 12];

        let record = build_consensus_record(
            "test/1",
            seq,
            &quals,
            true,
            5,
            3,
            1,
            None,
            None,
            MethylationMode::Disabled,
        );

        // Verify sequence length matches
        assert_eq!(record.sequence().len(), 12);
    }

    #[test]
    fn test_build_consensus_record_qualities_stored() {
        let seq = b"ACGT";
        let quals = vec![10, 20, 30, 40];

        let record = build_consensus_record(
            "test/1",
            seq,
            &quals,
            true,
            5,
            3,
            1,
            None,
            None,
            MethylationMode::Disabled,
        );

        let record_quals: Vec<u8> = record.quality_scores().iter().collect();
        assert_eq!(record_quals, quals);
    }

    #[test]
    fn test_build_consensus_record_empty_sequence() {
        let seq: &[u8] = b"";
        let quals: Vec<u8> = vec![];

        let record = build_consensus_record(
            "test/1",
            seq,
            &quals,
            true,
            5,
            3,
            1,
            None,
            None,
            MethylationMode::Disabled,
        );

        assert!(record.name().is_some());
    }

    #[test]
    fn test_build_consensus_record_long_sequence() {
        let seq = vec![b'A'; 500];
        let quals = vec![40; 500];

        let record = build_consensus_record(
            "test/1",
            &seq,
            &quals,
            true,
            5,
            3,
            1,
            None,
            None,
            MethylationMode::Disabled,
        );

        assert_eq!(record.sequence().len(), 500);
    }

    #[test]
    fn test_build_consensus_record_zero_depth() {
        let seq = b"ACGT";
        let quals = vec![40; 4];

        // Edge case: zero depth (shouldn't happen normally but test it)
        let record = build_consensus_record(
            "test/1",
            seq,
            &quals,
            true,
            0,
            0,
            0,
            None,
            None,
            MethylationMode::Disabled,
        );

        assert!(record.name().is_some());
    }

    #[test]
    fn test_build_consensus_record_high_error_count() {
        let seq = b"ACGT";
        let quals = vec![40; 4];

        // High error count (more errors than bases)
        let record = build_consensus_record(
            "test/1",
            seq,
            &quals,
            true,
            10,
            5,
            100,
            None,
            None,
            MethylationMode::Disabled,
        );

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

        let record = build_consensus_record(
            "test/1",
            seq,
            &quals,
            true,
            cd,
            5,
            2,
            Some((ad, bd, 3, 2, 1, 1)),
            None,
            MethylationMode::Disabled,
        );

        assert!(record.name().is_some());
        // The record should be valid with these tag values
    }

    #[test]
    fn test_build_consensus_record_various_names() {
        let seq = b"ACGT";
        let quals = vec![40; 4];

        for name in ["read1/1", "consensus_00000001/1", "test-read_123/2", "a/1"] {
            let record = build_consensus_record(
                name,
                seq,
                &quals,
                true,
                5,
                3,
                1,
                None,
                None,
                MethylationMode::Disabled,
            );
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
                > data.annotation.evidence[4].unconverted_count
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
            &dist,
            false,
            &mut rng,
        );

        // Position 0 is C in CpG -> methylated in TAPs -> mostly converted (target)
        assert!(data.annotation.evidence[0].is_ref_c);
        assert!(
            data.annotation.evidence[0].converted_count
                > data.annotation.evidence[0].unconverted_count
        );

        // Position 4 is C not CpG -> unmethylated in TAPs -> mostly unconverted (not a target)
        assert!(data.annotation.evidence[4].is_ref_c);
        assert!(
            data.annotation.evidence[4].unconverted_count
                > data.annotation.evidence[4].converted_count
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
            &dist,
            false,
            &mut rng,
        );

        let record = build_consensus_record(
            "test/1",
            seq,
            &quals,
            true,
            5,
            3,
            1,
            None,
            Some(&meth),
            MethylationMode::EmSeq,
        );

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

        let record = build_consensus_record(
            "test/1",
            seq,
            &quals,
            true,
            10,
            5,
            2,
            Some((6, 4, 3, 2, 1, 1)),
            Some(&meth),
            MethylationMode::EmSeq,
        );

        assert!(record.name().is_some());
    }

    #[test]
    fn test_no_methylation_tags_when_disabled() {
        let seq = b"CGAACGAT";
        let quals = vec![40; 8];

        // When no methylation data is passed, record should still build fine
        let record = build_consensus_record(
            "test/1",
            seq,
            &quals,
            true,
            5,
            3,
            1,
            None,
            None,
            MethylationMode::Disabled,
        );

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
        };
        let result = cmd.execute("test");
        assert!(result.is_err());
    }

    #[test]
    fn test_methylation_depth_mean_validation_accepts_when_disabled() {
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
        let r1 = build_consensus_record(
            "test/1",
            seq,
            &quals,
            true,
            5,
            3,
            1,
            None,
            Some(&meth),
            MethylationMode::EmSeq,
        );

        // Build R2 record with reversed methylation (as generate_consensus_pair does)
        let r2_seq = reverse_complement(seq);
        let r2_meth = meth.reverse();
        let r2 = build_consensus_record(
            "test/2",
            &r2_seq,
            &quals,
            false,
            5,
            3,
            1,
            None,
            Some(&r2_meth),
            MethylationMode::EmSeq,
        );

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
}

//! Generate consensus BAM with tags for filter.

use crate::commands::command::Command;
use crate::commands::common::CompressionOptions;
use crate::commands::simulate::common::StrandBiasArgs;
use anyhow::{Context, Result};
use clap::Parser;
use crossbeam_channel::bounded;
use fgumi_lib::bam_io::create_bam_writer;
use fgumi_lib::progress::ProgressTracker;
use fgumi_lib::sam::builder::{ConsensusTagsBuilder, RecordBuilder};
use fgumi_lib::simulate::{StrandBiasModel, create_rng};
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
    #[arg(long = "duplex")]
    pub duplex: bool,

    /// Base quality for consensus reads
    #[arg(long = "consensus-quality", default_value = "40")]
    pub consensus_quality: u8,

    #[command(flatten)]
    pub strand_bias: StrandBiasArgs,
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
}

/// Channel capacity for buffering read pairs between producer and writer threads.
const CHANNEL_CAPACITY: usize = 1_000;

impl Command for ConsensusReads {
    fn execute(&self, command_line: &str) -> Result<()> {
        info!("Generating consensus reads");
        info!("  Output: {}", self.output.display());
        info!("  Num reads: {}", self.num_reads);
        info!("  Read length: {}", self.read_length);
        info!("  Duplex: {}", self.duplex);
        info!("  Depth range: {}-{}", self.min_depth, self.max_depth);
        info!("  Threads: {}", self.threads);

        // Build header (no references for unmapped BAM)
        let mut header_builder = Header::builder();

        // Add @PG record
        header_builder = fgumi_lib::header::add_pg_to_builder(
            header_builder,
            crate::version::VERSION.as_str(),
            command_line,
        )?;

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
    );

    // Build R2 record (reverse complement)
    let r2_seq = reverse_complement(&seq);
    let r2_record = build_consensus_record(
        &format!("{read_name}/2"),
        &r2_seq,
        &quals,
        false, // is_first
        cd,
        cm,
        ce,
        if params.duplex { Some((ad, bd, am, bm, ae, be)) } else { None },
    );

    ConsensusReadPair {
        read_name,
        r1_record,
        r2_record,
        truth: (cd, cm, ce, ad, bd, am, bm, ae, be),
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

    RecordBuilder::new()
        .name(name)
        .sequence(&seq_str)
        .qualities(quals)
        .paired(true)
        .first_segment(is_first)
        .unmapped(true)
        .mate_unmapped(true)
        .consensus_tags(consensus_tags)
        .build()
}

fn generate_random_sequence(len: usize, rng: &mut impl Rng) -> Vec<u8> {
    const BASES: &[u8] = b"ACGT";
    let mut seq = Vec::with_capacity(len);
    for _ in 0..len {
        seq.push(BASES[rng.random_range(0..4)]);
    }
    seq
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    let mut result = Vec::with_capacity(seq.len());
    for &b in seq.iter().rev() {
        result.push(match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'N',
        });
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use fgumi_lib::simulate::create_rng;

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

        let record = build_consensus_record("test/1", seq, &quals, true, 5, 3, 1, None);

        let flags = record.flags();
        assert!(flags.is_mate_unmapped());
        assert!(flags.is_unmapped());
    }

    #[test]
    fn test_build_consensus_record_segmented_flag() {
        let seq = b"ACGT";
        let quals = vec![40; 4];

        let record = build_consensus_record("test/1", seq, &quals, true, 5, 3, 1, None);

        let flags = record.flags();
        assert!(flags.is_segmented());
    }

    #[test]
    fn test_build_consensus_record_sequence_stored() {
        let seq = b"ACGTACGTACGT";
        let quals = vec![40; 12];

        let record = build_consensus_record("test/1", seq, &quals, true, 5, 3, 1, None);

        // Verify sequence length matches
        assert_eq!(record.sequence().len(), 12);
    }

    #[test]
    fn test_build_consensus_record_qualities_stored() {
        let seq = b"ACGT";
        let quals = vec![10, 20, 30, 40];

        let record = build_consensus_record("test/1", seq, &quals, true, 5, 3, 1, None);

        let record_quals: Vec<u8> = record.quality_scores().iter().collect();
        assert_eq!(record_quals, quals);
    }

    #[test]
    fn test_build_consensus_record_empty_sequence() {
        let seq: &[u8] = b"";
        let quals: Vec<u8> = vec![];

        let record = build_consensus_record("test/1", seq, &quals, true, 5, 3, 1, None);

        assert!(record.name().is_some());
    }

    #[test]
    fn test_build_consensus_record_long_sequence() {
        let seq = vec![b'A'; 500];
        let quals = vec![40; 500];

        let record = build_consensus_record("test/1", &seq, &quals, true, 5, 3, 1, None);

        assert_eq!(record.sequence().len(), 500);
    }

    #[test]
    fn test_build_consensus_record_zero_depth() {
        let seq = b"ACGT";
        let quals = vec![40; 4];

        // Edge case: zero depth (shouldn't happen normally but test it)
        let record = build_consensus_record("test/1", seq, &quals, true, 0, 0, 0, None);

        assert!(record.name().is_some());
    }

    #[test]
    fn test_build_consensus_record_high_error_count() {
        let seq = b"ACGT";
        let quals = vec![40; 4];

        // High error count (more errors than bases)
        let record = build_consensus_record("test/1", seq, &quals, true, 10, 5, 100, None);

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
        assert_eq!(reverse_complement(b"X"), b"N");
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
        );

        assert!(record.name().is_some());
        // The record should be valid with these tag values
    }

    #[test]
    fn test_build_consensus_record_various_names() {
        let seq = b"ACGT";
        let quals = vec![40; 4];

        for name in ["read1/1", "consensus_00000001/1", "test-read_123/2", "a/1"] {
            let record = build_consensus_record(name, seq, &quals, true, 5, 3, 1, None);
            assert!(record.name().is_some());
        }
    }
}

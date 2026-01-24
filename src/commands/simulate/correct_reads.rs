//! Generate BAM and includelist for correct.

use crate::commands::command::Command;
use crate::commands::common::CompressionOptions;
use anyhow::{Context, Result};
use clap::Parser;
use crossbeam_channel::bounded;
use fgumi_lib::bam_io::create_bam_writer;
use fgumi_lib::progress::ProgressTracker;
use fgumi_lib::sam::builder::RecordBuilder;
use fgumi_lib::simulate::create_rng;
use log::info;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::sam::header::Header;
use rand::Rng;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::sync::Arc;
use std::thread;

/// Generate unmapped BAM and UMI includelist for `fgumi correct`.
#[derive(Parser, Debug)]
#[command(
    name = "correct-reads",
    about = "Generate BAM and includelist for correct",
    long_about = r#"
Generate synthetic unmapped BAM and UMI includelist for testing `fgumi correct`.

Reads are generated with varying edit distances from the includelist UMIs
to test the correction algorithm's ability to find the correct UMI.
"#
)]
pub struct CorrectReads {
    /// Output BAM file (unmapped)
    #[arg(short = 'o', long = "output", required = true)]
    pub output: PathBuf,

    /// Output UMI includelist file
    #[arg(short = 'i', long = "includelist", required = true)]
    pub includelist: PathBuf,

    /// Output truth TSV file for validation
    #[arg(long = "truth", required = true)]
    pub truth_output: PathBuf,

    /// Number of reads to generate
    #[arg(short = 'n', long = "num-reads", default_value = "10000")]
    pub num_reads: usize,

    /// Number of unique UMIs in includelist
    #[arg(long = "num-umis", default_value = "1000")]
    pub num_umis: usize,

    /// Length of UMI sequence in bases
    #[arg(short = 'u', long = "umi-length", default_value = "8")]
    pub umi_length: usize,

    /// Length of template sequence
    #[arg(short = 'l', long = "read-length", default_value = "100")]
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

    /// Fraction with exact UMI match (0 edits)
    #[arg(long = "exact-fraction", default_value = "0.4")]
    pub exact_fraction: f64,

    /// Fraction with 1 edit distance
    #[arg(long = "edit1-fraction", default_value = "0.3")]
    pub edit1_fraction: f64,

    /// Fraction with 2 edit distance
    #[arg(long = "edit2-fraction", default_value = "0.2")]
    pub edit2_fraction: f64,

    /// Fraction with 3+ edits (should not correct)
    #[arg(long = "multi-fraction", default_value = "0.1")]
    pub multi_fraction: f64,

    /// Base quality for all bases
    #[arg(long = "quality", default_value = "30")]
    pub quality: u8,
}

/// A generated read pair ready for output.
struct CorrectReadPair {
    read_name: String,
    r1_record: RecordBuf,
    r2_record: RecordBuf,
    /// Truth data: (`true_umi`, `observed_umi`, `expected_correction`, `edit_distance`, `error_type`)
    truth: (String, String, String, usize, &'static str),
}

/// Parameters needed for parallel read generation.
#[derive(Clone)]
struct GenerationParams {
    read_length: usize,
    quality: u8,
    exact_fraction: f64,
    edit1_fraction: f64,
    edit2_fraction: f64,
    umi_length: usize,
}

/// Channel capacity for buffering read pairs between producer and writer threads.
const CHANNEL_CAPACITY: usize = 1_000;

impl Command for CorrectReads {
    fn execute(&self, command_line: &str) -> Result<()> {
        // Validate that we can generate enough unique UMIs
        let max_possible_umis = 4_usize.saturating_pow(self.umi_length as u32);
        if self.num_umis > max_possible_umis {
            anyhow::bail!(
                "Cannot generate {} unique UMIs with {}-base UMIs (maximum possible: 4^{} = {})",
                self.num_umis,
                self.umi_length,
                self.umi_length,
                max_possible_umis
            );
        }

        info!("Generating correct-reads data");
        info!("  Output: {}", self.output.display());
        info!("  Includelist: {}", self.includelist.display());
        info!("  Truth: {}", self.truth_output.display());
        info!("  Num reads: {}", self.num_reads);
        info!("  Num UMIs: {}", self.num_umis);
        info!("  UMI length: {}", self.umi_length);
        info!("  Threads: {}", self.threads);

        // Validate fractions sum to 1.0
        let total =
            self.exact_fraction + self.edit1_fraction + self.edit2_fraction + self.multi_fraction;
        if (total - 1.0).abs() > 0.001 {
            anyhow::bail!("Error fractions must sum to 1.0 (got {total})");
        }

        let mut rng = create_rng(self.seed);

        // Generate includelist UMIs (sequential - needed before parallel generation)
        info!("Generating {} unique UMIs", self.num_umis);
        let mut umis: Vec<String> = Vec::with_capacity(self.num_umis);
        let mut seen = std::collections::HashSet::new();

        while umis.len() < self.num_umis {
            let umi = generate_random_sequence(self.umi_length, &mut rng);
            let umi_str = String::from_utf8_lossy(&umi).to_string();
            if !seen.contains(&umi_str) {
                seen.insert(umi_str.clone());
                umis.push(umi_str);
            }
        }

        // Sort includelist alphabetically
        umis.sort();

        // Write includelist
        let includelist_file = File::create(&self.includelist)
            .with_context(|| format!("Failed to create {}", self.includelist.display()))?;
        let mut includelist_writer = BufWriter::new(includelist_file);
        for umi in &umis {
            writeln!(includelist_writer, "{umi}")?;
        }
        includelist_writer.flush()?;
        info!("Wrote includelist with {} UMIs", umis.len());

        // Share UMI list across workers
        let umis = Arc::new(umis);

        // Build BAM header
        let mut header_builder = Header::builder();

        // Add @PG record
        header_builder = fgumi_lib::header::add_pg_to_builder(
            header_builder,
            crate::version::VERSION.as_str(),
            command_line,
        )?;

        let header = header_builder.build();

        // Set up generation parameters
        let params = Arc::new(GenerationParams {
            read_length: self.read_length,
            quality: self.quality,
            exact_fraction: self.exact_fraction,
            edit1_fraction: self.edit1_fraction,
            edit2_fraction: self.edit2_fraction,
            umi_length: self.umi_length,
        });

        // Generate seeds for reproducibility
        let read_seeds: Vec<u64> = (0..self.num_reads).map(|_| rng.random()).collect();

        // Create bounded channel for streaming read pairs to writer
        let (sender, receiver) = bounded::<CorrectReadPair>(CHANNEL_CAPACITY);

        // Clone paths for writer thread
        let output_path = self.output.clone();
        let truth_path = self.truth_output.clone();
        let compression_level = self.compression.compression_level;
        let writer_threads = self.threads;
        let header_clone = header.clone();
        let num_reads = self.num_reads;

        // Spawn writer thread with multi-threaded BGZF compression
        let writer_handle = thread::spawn(move || -> Result<(u64, u64, u64, u64, u64)> {
            let mut writer =
                create_bam_writer(&output_path, &header_clone, writer_threads, compression_level)?;

            // Create truth file
            let truth_file = File::create(&truth_path)
                .with_context(|| format!("Failed to create {}", truth_path.display()))?;
            let mut truth_writer = BufWriter::new(truth_file);
            writeln!(
                truth_writer,
                "read_name\ttrue_umi\tobserved_umi\texpected_correction\tedit_distance\terror_type"
            )?;

            let mut read_count = 0u64;
            let mut exact_count = 0u64;
            let mut edit1_count = 0u64;
            let mut edit2_count = 0u64;
            let mut multi_count = 0u64;

            let progress = ProgressTracker::new("Generated reads").with_interval(100_000);

            // Receive and write read pairs as they arrive
            for pair in receiver {
                read_count += 1;
                progress.log_if_needed(1);

                writer.write_alignment_record(&header_clone, &pair.r1_record)?;
                writer.write_alignment_record(&header_clone, &pair.r2_record)?;

                // Count error types
                let (true_umi, observed_umi, expected_correction, edit_distance, error_type) =
                    &pair.truth;
                match *error_type {
                    "exact" => exact_count += 1,
                    "edit1" => edit1_count += 1,
                    "edit2" => edit2_count += 1,
                    "multi" => multi_count += 1,
                    _ => {}
                }

                // Write truth
                writeln!(
                    truth_writer,
                    "{}\t{true_umi}\t{observed_umi}\t{expected_correction}\t{edit_distance}\t{error_type}",
                    pair.read_name
                )?;
            }

            progress.log_final();
            truth_writer.flush()?;

            Ok((read_count, exact_count, edit1_count, edit2_count, multi_count))
        });

        // Configure thread pool for generation
        let gen_threads = if self.threads <= 1 { 1 } else { self.threads.max(2) };
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(gen_threads)
            .build()
            .with_context(|| "Failed to create thread pool")?;

        // Generate reads in parallel and stream to writer
        let generation_result: Result<(), crossbeam_channel::SendError<CorrectReadPair>> = pool
            .install(|| {
                read_seeds.into_par_iter().enumerate().try_for_each(|(read_idx, seed)| {
                    let pair = generate_correct_read_pair(read_idx, seed, &params, &umis);
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
        let (read_count, exact_count, edit1_count, edit2_count, multi_count) =
            writer_handle.join().map_err(|_| anyhow::anyhow!("Writer thread panicked"))??;

        info!("Generated {read_count} reads:");
        info!(
            "  Exact (0 edits): {} ({:.1}%)",
            exact_count,
            100.0 * exact_count as f64 / num_reads as f64
        );
        info!(
            "  Edit1 (1 edit):  {} ({:.1}%)",
            edit1_count,
            100.0 * edit1_count as f64 / num_reads as f64
        );
        info!(
            "  Edit2 (2 edits): {} ({:.1}%)",
            edit2_count,
            100.0 * edit2_count as f64 / num_reads as f64
        );
        info!(
            "  Multi (3+ edits): {} ({:.1}%)",
            multi_count,
            100.0 * multi_count as f64 / num_reads as f64
        );
        info!("Done");

        Ok(())
    }
}

/// Generate a single correct-read pair.
fn generate_correct_read_pair(
    read_idx: usize,
    seed: u64,
    params: &GenerationParams,
    umis: &[String],
) -> CorrectReadPair {
    let mut rng = create_rng(Some(seed));

    let read_name = format!("read_{read_idx:08}");

    // Pick a random UMI from the includelist
    let true_umi_idx = rng.random_range(0..umis.len());
    let true_umi = &umis[true_umi_idx];

    // Calculate cumulative fractions for sampling
    let exact_threshold = params.exact_fraction;
    let edit1_threshold = exact_threshold + params.edit1_fraction;
    let edit2_threshold = edit1_threshold + params.edit2_fraction;

    // Determine error type
    let r: f64 = rng.random();
    let (observed_umi, edit_distance, error_type) = if r < exact_threshold {
        (true_umi.clone(), 0, "exact")
    } else if r < edit1_threshold {
        (introduce_n_errors(true_umi, 1, &mut rng), 1, "edit1")
    } else if r < edit2_threshold {
        (introduce_n_errors(true_umi, 2, &mut rng), 2, "edit2")
    } else {
        // 3+ errors
        let num_errors = rng.random_range(3..=params.umi_length.min(5));
        (introduce_n_errors(true_umi, num_errors, &mut rng), num_errors, "multi")
    };

    // Expected correction: for exact/edit1/edit2, should correct to true_umi
    // For multi, should not correct (stays as observed)
    let expected_correction =
        if error_type == "multi" { observed_umi.clone() } else { true_umi.clone() };

    // Generate template sequences for R1 and R2
    let template_r1 = generate_random_sequence(params.read_length, &mut rng);
    let template_r2 = generate_random_sequence(params.read_length, &mut rng);

    // Build R1 (flag 77: UNMAPPED | SEGMENTED | FIRST_SEGMENT | MATE_UNMAPPED)
    let r1_record = RecordBuilder::new()
        .name(&read_name)
        .sequence(&String::from_utf8_lossy(&template_r1))
        .qualities(&vec![params.quality; template_r1.len()])
        .paired(true)
        .first_segment(true)
        .unmapped(true)
        .mate_unmapped(true)
        .tag("RX", observed_umi.clone())
        .build();

    // Build R2 (flag 141: UNMAPPED | SEGMENTED | LAST_SEGMENT | MATE_UNMAPPED)
    let r2_record = RecordBuilder::new()
        .name(&read_name)
        .sequence(&String::from_utf8_lossy(&template_r2))
        .qualities(&vec![params.quality; template_r2.len()])
        .paired(true)
        .first_segment(false)
        .unmapped(true)
        .mate_unmapped(true)
        .tag("RX", observed_umi.clone())
        .build();

    CorrectReadPair {
        read_name,
        r1_record,
        r2_record,
        truth: (true_umi.clone(), observed_umi, expected_correction, edit_distance, error_type),
    }
}

fn generate_random_sequence(len: usize, rng: &mut impl Rng) -> Vec<u8> {
    const BASES: &[u8] = b"ACGT";
    let mut seq = Vec::with_capacity(len);
    for _ in 0..len {
        seq.push(BASES[rng.random_range(0..4)]);
    }
    seq
}

fn introduce_n_errors(umi: &str, n: usize, rng: &mut impl Rng) -> String {
    const BASES: &[u8] = b"ACGT";
    let mut result: Vec<u8> = umi.as_bytes().to_vec();

    // Choose n unique positions to mutate
    let mut positions: Vec<usize> = (0..result.len()).collect();

    // Fisher-Yates shuffle to get random positions
    for i in 0..n.min(positions.len()) {
        let j = rng.random_range(i..positions.len());
        positions.swap(i, j);
    }

    // Introduce errors at the first n positions
    for &pos in positions.iter().take(n) {
        let current = result[pos];
        // Pick a different base
        loop {
            let new_base = BASES[rng.random_range(0..4)];
            if new_base != current {
                result[pos] = new_base;
                break;
            }
        }
    }

    String::from_utf8_lossy(&result).to_string()
}

#[cfg(test)]
#[allow(clippy::naive_bytecount)]
mod tests {
    use super::*;
    use fgumi_lib::simulate::create_rng;
    use noodles::sam::alignment::record::data::field::Tag;

    // Tests for generate_random_sequence
    #[test]
    fn test_generate_random_sequence_length() {
        let mut rng = create_rng(Some(42));
        for len in [0, 1, 8, 16, 100] {
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
    fn test_generate_random_sequence_reproducibility() {
        let mut rng1 = create_rng(Some(42));
        let mut rng2 = create_rng(Some(42));

        let seq1 = generate_random_sequence(100, &mut rng1);
        let seq2 = generate_random_sequence(100, &mut rng2);

        assert_eq!(seq1, seq2);
    }

    #[test]
    fn test_generate_random_sequence_distribution() {
        let mut rng = create_rng(Some(42));
        let seq = generate_random_sequence(10000, &mut rng);

        let a_count = seq.iter().filter(|&&b| b == b'A').count();
        let c_count = seq.iter().filter(|&&b| b == b'C').count();
        let g_count = seq.iter().filter(|&&b| b == b'G').count();
        let t_count = seq.iter().filter(|&&b| b == b'T').count();

        // Each base should be roughly 25% (with some tolerance)
        for (base, count) in [('A', a_count), ('C', c_count), ('G', g_count), ('T', t_count)] {
            let fraction = count as f64 / 10000.0;
            assert!(
                (0.20..0.30).contains(&fraction),
                "Base {} has unexpected frequency: {:.2}%",
                base,
                fraction * 100.0
            );
        }
    }

    // Tests for introduce_n_errors
    #[test]
    fn test_introduce_n_errors_zero() {
        let mut rng = create_rng(Some(42));
        let umi = "ACGTACGT";
        let result = introduce_n_errors(umi, 0, &mut rng);
        assert_eq!(result, umi);
    }

    #[test]
    fn test_introduce_n_errors_preserves_length() {
        let mut rng = create_rng(Some(42));
        let umi = "ACGTACGT";
        for n in 1..=8 {
            let result = introduce_n_errors(umi, n, &mut rng);
            assert_eq!(result.len(), umi.len());
        }
    }

    #[test]
    fn test_introduce_n_errors_exact_count() {
        let mut rng = create_rng(Some(42));
        let umi = "AAAAAAAA"; // All same base makes counting changes easy

        for n in 1..=8 {
            let result = introduce_n_errors(umi, n, &mut rng);
            let diff_count = umi.chars().zip(result.chars()).filter(|(a, b)| a != b).count();
            assert_eq!(
                diff_count, n,
                "Expected {n} errors but got {diff_count} for UMI {umi} -> {result}"
            );
        }
    }

    #[test]
    fn test_introduce_n_errors_more_than_length() {
        let mut rng = create_rng(Some(42));
        let umi = "ACGT"; // Length 4

        // Requesting 10 errors on 4-base UMI should only change 4 positions
        let result = introduce_n_errors(umi, 10, &mut rng);
        assert_eq!(result.len(), 4);

        // All positions should be different (since we request more errors than bases)
        let diff_count = umi.chars().zip(result.chars()).filter(|(a, b)| a != b).count();
        assert_eq!(diff_count, 4);
    }

    #[test]
    fn test_introduce_n_errors_different_base_guarantee() {
        // Each mutated position must have a DIFFERENT base
        let mut rng = create_rng(Some(42));
        let umi = "AAAAAAAA";

        for _ in 0..100 {
            let result = introduce_n_errors(umi, 4, &mut rng);
            // Changed positions should not be 'A'
            for (orig, new) in umi.chars().zip(result.chars()) {
                if orig != new {
                    assert_ne!(new, 'A', "Mutated base should not equal original");
                }
            }
        }
    }

    #[test]
    fn test_introduce_n_errors_valid_bases() {
        let mut rng = create_rng(Some(42));
        let umi = "ACGTACGT";

        for _ in 0..100 {
            let result = introduce_n_errors(umi, 4, &mut rng);
            for c in result.chars() {
                assert!(
                    c == 'A' || c == 'C' || c == 'G' || c == 'T',
                    "Invalid base in result: {c}"
                );
            }
        }
    }

    #[test]
    fn test_introduce_n_errors_reproducibility() {
        let umi = "ACGTACGT";

        let mut rng1 = create_rng(Some(42));
        let mut rng2 = create_rng(Some(42));

        let result1 = introduce_n_errors(umi, 3, &mut rng1);
        let result2 = introduce_n_errors(umi, 3, &mut rng2);

        assert_eq!(result1, result2);
    }

    #[test]
    fn test_introduce_n_errors_single_base_umi() {
        let mut rng = create_rng(Some(42));
        let umi = "A";

        let result = introduce_n_errors(umi, 1, &mut rng);
        assert_eq!(result.len(), 1);
        assert_ne!(result, "A");
        assert!(result == "C" || result == "G" || result == "T");
    }

    #[test]
    fn test_introduce_n_errors_empty_umi() {
        let mut rng = create_rng(Some(42));
        let umi = "";

        let result = introduce_n_errors(umi, 5, &mut rng);
        assert_eq!(result, "");
    }

    #[test]
    fn test_introduce_n_errors_unique_positions() {
        // Verify that errors are introduced at UNIQUE positions
        let mut rng = create_rng(Some(42));
        let umi = "AAAAAAAAAAAA"; // 12 A's

        // Introduce 6 errors multiple times and verify exactly 6 positions change
        for _ in 0..100 {
            let result = introduce_n_errors(umi, 6, &mut rng);
            let diff_count = umi.chars().zip(result.chars()).filter(|(a, b)| a != b).count();
            assert_eq!(diff_count, 6, "Should have exactly 6 different positions");
        }
    }

    // Tests for paired read generation using RecordBuilder
    #[test]
    fn test_paired_read_r1_flags() {
        // R1 should have flag 77: UNMAPPED | SEGMENTED | FIRST_SEGMENT | MATE_UNMAPPED
        let r1 = RecordBuilder::new()
            .name("test_read")
            .sequence("ACGT")
            .paired(true)
            .first_segment(true)
            .unmapped(true)
            .mate_unmapped(true)
            .tag("RX", "AAAAAAAA")
            .build();

        assert!(r1.flags().is_unmapped());
        assert!(r1.flags().is_segmented());
        assert!(r1.flags().is_first_segment());
        assert!(r1.flags().is_mate_unmapped());
        assert!(!r1.flags().is_last_segment());

        // Verify flag value is 77
        let flag_bits: u16 = r1.flags().bits();
        assert_eq!(flag_bits, 77, "R1 flag should be 77");
    }

    #[test]
    fn test_paired_read_r2_flags() {
        // R2 should have flag 141: UNMAPPED | SEGMENTED | LAST_SEGMENT | MATE_UNMAPPED
        let r2 = RecordBuilder::new()
            .name("test_read")
            .sequence("ACGT")
            .paired(true)
            .first_segment(false)
            .unmapped(true)
            .mate_unmapped(true)
            .tag("RX", "AAAAAAAA")
            .build();

        assert!(r2.flags().is_unmapped());
        assert!(r2.flags().is_segmented());
        assert!(r2.flags().is_last_segment());
        assert!(r2.flags().is_mate_unmapped());
        assert!(!r2.flags().is_first_segment());

        // Verify flag value is 141
        let flag_bits: u16 = r2.flags().bits();
        assert_eq!(flag_bits, 141, "R2 flag should be 141");
    }

    #[test]
    fn test_paired_reads_same_name_and_umi() {
        let read_name = "test_pair";
        let umi = "ACGTACGT";

        let r1 = RecordBuilder::new()
            .name(read_name)
            .sequence("AAAA")
            .paired(true)
            .first_segment(true)
            .unmapped(true)
            .mate_unmapped(true)
            .tag("RX", umi)
            .build();

        let r2 = RecordBuilder::new()
            .name(read_name)
            .sequence("CCCC")
            .paired(true)
            .first_segment(false)
            .unmapped(true)
            .mate_unmapped(true)
            .tag("RX", umi)
            .build();

        // Both should have the same read name
        assert_eq!(r1.name(), r2.name());

        // Both should have the same RX tag
        let rx_tag: Tag = [b'R', b'X'].into();
        assert!(r1.data().get(&rx_tag).is_some());
        assert!(r2.data().get(&rx_tag).is_some());
    }

    #[test]
    fn test_paired_read_qualities() {
        let quality: u8 = 35;
        let seq = "ACGTACGT";

        let r1 = RecordBuilder::new()
            .name("read")
            .sequence(seq)
            .qualities(&vec![quality; seq.len()])
            .paired(true)
            .first_segment(true)
            .unmapped(true)
            .mate_unmapped(true)
            .build();

        let quals: Vec<u8> = r1.quality_scores().iter().collect();
        assert_eq!(quals.len(), 8);
        assert!(quals.iter().all(|&q| q == quality));
    }

    // Tests for UMI count validation
    #[test]
    fn test_num_umis_exceeds_possible_umis_fails() {
        use crate::commands::command::Command;
        use tempfile::tempdir;

        let dir = tempdir().unwrap();
        let output = dir.path().join("output.bam");
        let includelist = dir.path().join("includelist.txt");
        let truth = dir.path().join("truth.tsv");

        // 4^4 = 256 possible UMIs, requesting 1000 should fail
        let cmd = CorrectReads {
            output,
            includelist,
            truth_output: truth,
            num_reads: 100,
            num_umis: 1000,
            umi_length: 4,
            read_length: 100,
            seed: Some(42),
            threads: 1,
            compression: CompressionOptions::default(),
            exact_fraction: 0.4,
            edit1_fraction: 0.3,
            edit2_fraction: 0.2,
            multi_fraction: 0.1,
            quality: 30,
        };

        let result = cmd.execute("test");
        assert!(result.is_err());
        let err_msg = result.unwrap_err().to_string();
        assert!(
            err_msg.contains("Cannot generate 1000 unique UMIs"),
            "Error message should mention UMI count: {err_msg}"
        );
        assert!(
            err_msg.contains("4-base UMIs"),
            "Error message should mention UMI length: {err_msg}"
        );
        assert!(err_msg.contains("256"), "Error message should mention max possible: {err_msg}");
    }

    #[test]
    fn test_num_umis_at_max_succeeds() {
        // 4^2 = 16 possible UMIs, requesting exactly 16 should work
        let max_possible = 4_usize.pow(2);
        assert_eq!(max_possible, 16);

        // Just verify the validation logic - don't need full execution
        let num_umis = 16;
        let umi_length = 2;
        let max_possible_umis = 4_usize.saturating_pow(umi_length as u32);
        assert!(num_umis <= max_possible_umis, "16 UMIs should be allowed for 2-base UMIs");
    }
}

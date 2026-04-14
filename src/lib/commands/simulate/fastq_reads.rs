//! Generate paired-end FASTQ files with UMI sequences.

use crate::commands::command::Command;
use crate::commands::common::parse_bool;
use crate::commands::simulate::common::{
    FamilySizeArgs, InsertSizeArgs, MethylationArgs, MethylationConfig, QualityArgs,
    ReferenceGenome, SimulationCommon, apply_methylation_conversion, generate_random_sequence,
};
use crate::progress::ProgressTracker;
use crate::simulate::{FastqWriter, create_rng};
use anyhow::{Context, Result, bail};
use clap::Parser;
use crossbeam_channel::bounded;
use fgumi_dna::dna::complement_base;
use log::{info, warn};
use rand::{Rng, RngExt};
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::thread;

/// Generate paired-end FASTQ files with UMI sequences for input to `fgumi extract`.
#[derive(Parser, Debug)]
#[command(
    name = "fastq-reads",
    about = "Generate R1/R2 FASTQ.gz files for extract",
    long_about = r#"
Generate synthetic paired-end FASTQ files with UMI sequences.

The output is suitable for input to `fgumi extract`. Each read pair contains:
- R1: UMI sequence + template (according to --read-structure-r1)
- R2: Template only (according to --read-structure-r2)

Quality scores follow a realistic model with warmup, peak, and decay phases.
"#
)]
pub struct FastqReads {
    /// Output R1 FASTQ.gz file
    #[arg(short = '1', long = "r1", required = true)]
    pub r1_output: PathBuf,

    /// Output R2 FASTQ.gz file
    #[arg(short = '2', long = "r2", required = true)]
    pub r2_output: PathBuf,

    /// Output truth TSV file for validation
    #[arg(long = "truth", required = true)]
    pub truth_output: PathBuf,

    /// Read structure for R1 (fgbio notation, e.g., "8M+T")
    #[arg(long = "read-structure-r1", default_value = "8M+T")]
    pub read_structure_r1: String,

    /// Read structure for R2 (fgbio notation, e.g., "+T")
    #[arg(long = "read-structure-r2", default_value = "+T")]
    pub read_structure_r2: String,

    /// Generate duplex-style reads (A/B strand pairs)
    #[arg(long = "duplex", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool)]
    pub duplex: bool,

    /// Reference FASTA file for sampling template sequences.
    #[arg(short = 'r', long = "reference", required = true)]
    pub reference: PathBuf,

    /// UMI includelist file (one UMI per line).
    /// When provided, UMIs are sampled from this list instead of generated randomly.
    /// The includelist also determines UMI length (overriding --umi-length).
    #[arg(short = 'i', long = "includelist")]
    pub includelist: Option<PathBuf>,

    /// Number of threads for parallel molecule generation
    #[arg(short = 't', long = "threads", default_value = "1")]
    pub threads: usize,

    #[command(flatten)]
    pub common: SimulationCommon,

    #[command(flatten)]
    pub quality: QualityArgs,

    #[command(flatten)]
    pub family_size: FamilySizeArgs,

    #[command(flatten)]
    pub insert_size: InsertSizeArgs,

    #[command(flatten)]
    pub methylation: MethylationArgs,
}

/// A single read record ready for output.
struct ReadRecord {
    read_name: String,
    r1_seq: Vec<u8>,
    r1_quals: Vec<u8>,
    r2_seq: Vec<u8>,
    r2_quals: Vec<u8>,
    umi: Vec<u8>,
    mol_id: usize,
    family_id: usize,
    strand: &'static str,
    /// Chromosome name.
    chrom: String,
    /// 0-based start position.
    pos: usize,
}

/// Parameters needed for parallel molecule generation.
#[derive(Clone)]
struct GenerationParams {
    umi_length: usize,
    read_length: usize,
    min_family_size: usize,
    r2_quality_offset: i8,
    /// Generate duplex-style reads (A/B strand pairs).
    duplex: bool,
    /// When set, UMIs are sampled from this list instead of generated randomly.
    includelist: Option<Vec<Vec<u8>>>,
    methylation: MethylationConfig,
}

/// Load UMI sequences from an includelist file (one UMI per line).
fn load_includelist(path: &Path) -> Result<Vec<Vec<u8>>> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open includelist: {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut umis: Vec<(usize, Vec<u8>)> = Vec::new(); // (1-based line number, UMI)

    for (line_num, line) in reader.lines().enumerate() {
        let line = line.with_context(|| format!("Failed to read line {}", line_num + 1))?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let umi = trimmed.as_bytes().to_ascii_uppercase();
        if !umi.iter().all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
            bail!("Invalid UMI at line {}: '{}' (only A, C, G, T allowed)", line_num + 1, trimmed,);
        }
        umis.push((line_num + 1, umi));
    }

    if umis.is_empty() {
        bail!("Includelist is empty: {}", path.display());
    }

    // Validate all UMIs have the same length
    let expected_len = umis[0].1.len();
    for (file_line, umi) in &umis {
        if umi.len() != expected_len {
            bail!(
                "UMI length mismatch at line {}: length {} but expected {} (from first UMI)",
                file_line,
                umi.len(),
                expected_len,
            );
        }
    }

    Ok(umis.into_iter().map(|(_, umi)| umi).collect())
}

/// Channel capacity for buffering molecule batches between producer and writer threads.
/// Each batch contains all reads for one molecule (typically 1-10 reads).
const CHANNEL_CAPACITY: usize = 1_000;

impl Command for FastqReads {
    fn execute(&self, _command_line: &str) -> Result<()> {
        let includelist = if let Some(ref path) = self.includelist {
            let umis = load_includelist(path)?;
            info!("Loaded {} UMIs from includelist (length={})", umis.len(), umis[0].len());
            Some(umis)
        } else {
            None
        };

        // Includelist determines UMI length, overriding --umi-length
        let umi_length = if let Some(ref list) = includelist {
            let list_len = list[0].len();
            if list_len > self.common.read_length {
                bail!(
                    "Includelist UMI length {} exceeds read length {}",
                    list_len,
                    self.common.read_length,
                );
            }
            if self.common.umi_length != list_len {
                warn!(
                    "--umi-length {} overridden by includelist UMI length {}",
                    self.common.umi_length, list_len,
                );
            }
            list_len
        } else {
            if self.common.umi_length > self.common.read_length {
                bail!(
                    "UMI length {} exceeds read length {}",
                    self.common.umi_length,
                    self.common.read_length,
                );
            }
            self.common.umi_length
        };

        // Validate methylation args
        let methylation = self.methylation.resolve();
        self.methylation.validate()?;

        info!("Generating FASTQ reads");
        info!("  Output R1: {}", self.r1_output.display());
        info!("  Output R2: {}", self.r2_output.display());
        info!("  Truth: {}", self.truth_output.display());
        info!("  Num molecules: {}", self.common.num_molecules);
        info!("  Read length: {}", self.common.read_length);
        info!("  UMI length: {}", umi_length);
        if let Some(ref path) = self.includelist {
            info!("  Includelist: {}", path.display());
        }
        info!("  Duplex: {}", self.duplex);
        info!("  Threads: {}", self.threads);
        info!("  Reference: {}", self.reference.display());
        if methylation.mode.is_enabled() {
            info!("  Methylation mode: {:?}", methylation.mode);
            info!("  CpG methylation rate: {}", methylation.cpg_methylation_rate);
            info!("  Conversion rate: {}", methylation.conversion_rate);
        }

        let reference = Arc::new(ReferenceGenome::load(&self.reference)?);

        // Validate that the reference has at least one contig >= read_length
        if reference.max_contig_length() < self.common.read_length {
            bail!(
                "No reference contig is >= read length ({} bp). \
                 The longest contig is {} bp. Use a larger reference or shorter --read-length.",
                self.common.read_length,
                reference.max_contig_length(),
            );
        }

        // Use at least 2 threads for generation if multi-threaded
        // Reserve some threads for gzip compression
        let gen_threads = if self.threads <= 1 { 1 } else { self.threads.max(2) };
        let compress_threads = self.threads;

        let quality_model = Arc::new(self.quality.to_quality_model());
        let quality_bias = Arc::new(self.quality.to_quality_bias());
        let family_dist = Arc::new(self.family_size.to_family_size_distribution()?);
        let insert_model = Arc::new(self.insert_size.to_insert_size_model());

        let params = Arc::new(GenerationParams {
            umi_length,
            read_length: self.common.read_length,
            min_family_size: self.family_size.min_family_size,
            r2_quality_offset: self.quality.r2_quality_offset,
            duplex: self.duplex,
            includelist,
            methylation,
        });

        // Generate molecule IDs with seeds for reproducibility
        let mut seed_rng = create_rng(self.common.seed);
        let mol_seeds: Vec<u64> =
            (0..self.common.num_molecules).map(|_| seed_rng.random()).collect();

        // Create bounded channel for streaming record batches to writer
        // Each batch contains all reads for one molecule
        let (sender, receiver) = bounded::<Vec<ReadRecord>>(CHANNEL_CAPACITY);

        // Clone paths for writer thread
        let r1_path = self.r1_output.clone();
        let r2_path = self.r2_output.clone();
        let truth_path = self.truth_output.clone();
        // Spawn writer thread with multi-threaded gzip compression
        let writer_handle = thread::spawn(move || -> Result<u64> {
            let mut r1_writer = FastqWriter::with_threads(&r1_path, compress_threads)?;
            let mut r2_writer = FastqWriter::with_threads(&r2_path, compress_threads)?;
            let truth_file = File::create(&truth_path)
                .with_context(|| format!("Failed to create {}", truth_path.display()))?;
            let mut truth_writer = BufWriter::new(truth_file);

            writeln!(
                truth_writer,
                "read_name\ttrue_umi\tmolecule_id\tfamily_id\tstrand\tchrom\tpos"
            )?;

            let mut read_count = 0u64;
            let progress = ProgressTracker::new("Generated read pairs").with_interval(100_000);

            // Receive and write record batches as they arrive
            for batch in receiver {
                for record in batch {
                    read_count += 1;
                    progress.log_if_needed(1);
                    r1_writer.write_record(&record.read_name, &record.r1_seq, &record.r1_quals)?;
                    r2_writer.write_record(&record.read_name, &record.r2_seq, &record.r2_quals)?;

                    // Write truth
                    writeln!(
                        truth_writer,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        record.read_name,
                        String::from_utf8_lossy(&record.umi),
                        record.mol_id,
                        record.family_id,
                        record.strand,
                        record.chrom,
                        record.pos
                    )?;
                }
            }

            progress.log_final();
            r1_writer.finish()?;
            r2_writer.finish()?;
            truth_writer.flush()?;

            Ok(read_count)
        });

        // Configure thread pool for generation
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(gen_threads)
            .build()
            .with_context(|| "Failed to create thread pool")?;

        // Generate reads in parallel and stream to writer
        let generation_result: Result<(), crossbeam_channel::SendError<Vec<ReadRecord>>> = pool
            .install(|| {
                mol_seeds.into_par_iter().enumerate().try_for_each(|(mol_id, seed)| {
                    let batch = generate_molecule_reads(
                        mol_id,
                        seed,
                        &params,
                        &quality_model,
                        &quality_bias,
                        &family_dist,
                        &insert_model,
                        &reference,
                    );
                    sender.send(batch)
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

        info!("Generated {read_count} read pairs");
        info!("Done");

        Ok(())
    }
}

/// Generate all reads for a single molecule.
#[allow(clippy::too_many_arguments)]
fn generate_molecule_reads(
    mol_id: usize,
    seed: u64,
    params: &GenerationParams,
    quality_model: &crate::simulate::PositionQualityModel,
    quality_bias: &crate::simulate::ReadPairQualityBias,
    family_dist: &crate::simulate::FamilySizeDistribution,
    insert_model: &crate::simulate::InsertSizeModel,
    reference: &ReferenceGenome,
) -> Vec<ReadRecord> {
    let mut rng = create_rng(Some(seed));
    let mut records = Vec::new();

    // Generate UMIs for this molecule (one for R1, one for R2)
    // For duplex sequencing with paired grouping strategy, B strand needs swapped UMIs
    let (umi_r1, umi_r2) = if let Some(ref includelist) = params.includelist {
        // Sample from the includelist
        let idx_r1 = rng.random_range(0..includelist.len());
        let idx_r2 = rng.random_range(0..includelist.len());
        (includelist[idx_r1].clone(), includelist[idx_r2].clone())
    } else {
        (
            generate_random_sequence(params.umi_length, &mut rng),
            generate_random_sequence(params.umi_length, &mut rng),
        )
    };

    // Generate family size
    let family_size = family_dist.sample(&mut rng, params.min_family_size);

    // Generate insert size for this molecule
    let insert_size = insert_model.sample(&mut rng);

    // Generate template sequence (shared by all reads in family)
    let template_len = insert_size.saturating_sub(params.umi_length);
    let (chrom_idx, pos, template) = reference
        .sample_sequence(template_len, &mut rng)
        .expect("Failed to sample sequence from reference");
    let chrom = reference.name(chrom_idx).to_string();

    // Pre-calculate template slicing bounds (both reads have UMI prefix)
    let template_len = params.read_length.saturating_sub(params.umi_length);
    let r2_start = template.len().saturating_sub(template_len);
    let r2_end = template.len();

    // Reusable buffers for sequence building
    let mut r1_seq_buf = Vec::with_capacity(params.read_length);
    let mut r2_seq_buf = Vec::with_capacity(params.read_length);

    let mut read_counter = 0usize;

    // For duplex mode, each read pair is independently assigned to A or B strand (coin flip)
    // For non-duplex mode, all reads are A strand (no UMI swapping or orientation change)
    for read_idx in 0..family_size {
        let (strand, family_idx): (&'static str, usize) = if params.duplex {
            if rng.random_bool(0.5) { ("A", 0) } else { ("B", 1) }
        } else {
            ("A", 0)
        };

        // For duplex paired strategy: A strand uses umi_r1-umi_r2, B strand uses umi_r2-umi_r1
        let (r1_umi, r2_umi) = if strand == "B" {
            (&umi_r2, &umi_r1) // Swapped for B strand
        } else {
            (&umi_r1, &umi_r2) // Normal for A strand
        };

        // Build combined UMI for truth output (format: R1_UMI-R2_UMI)
        let combined_umi: Vec<u8> = {
            let mut umi = r1_umi.clone();
            umi.push(b'-');
            umi.extend_from_slice(r2_umi);
            umi
        };

        read_counter += 1;
        let read_name = format!("mol{mol_id:06}_read{read_counter:04}");

        // Build R1 and R2 sequences differently based on strand
        // For A strand: R1 maps forward, R2 maps reverse (standard FR orientation)
        // For B strand: R1 maps reverse, R2 maps forward (RF orientation at same position)
        // This allows fgbio's paired strategy to find duplex pairs at the same position
        r1_seq_buf.clear();
        r2_seq_buf.clear();

        if strand == "A" {
            // A strand: standard forward orientation
            // R1 reads the top strand (forward) — apply C→T conversion
            r1_seq_buf.extend_from_slice(r1_umi);
            let r1_template_end = template_len.min(template.len());
            let mut r1_template = template[..r1_template_end].to_vec();
            apply_methylation_conversion(
                &mut r1_template,
                &template,
                0,
                true, // top strand
                &params.methylation,
                &mut rng,
            );
            r1_seq_buf.extend_from_slice(&r1_template);

            // R2 reads the bottom strand (reverse) — apply G→A conversion, then RC
            r2_seq_buf.extend_from_slice(r2_umi);
            let mut r2_template = template[r2_start..r2_end].to_vec();
            apply_methylation_conversion(
                &mut r2_template,
                &template,
                r2_start,
                false, // bottom strand
                &params.methylation,
                &mut rng,
            );
            reverse_complement_into(&r2_template, &mut r2_seq_buf);
        } else {
            // B strand: comes from the complementary DNA strand of the same molecule
            // For duplex sequencing, A and B strand reads should align to the SAME positions
            // but with opposite orientations (A=FR, B=RF).
            //
            // A strand: R1 at template START (forward), R2 at template END (reverse)
            // B strand: R1 at template END (reverse), R2 at template START (forward)
            //
            // This means:
            // - B R1 covers the same region as A R2 (template end), sequenced from revcomp
            // - B R2 covers the same region as A R1 (template start), sequenced forward

            // B R1 reads the bottom strand (template end, reverse) — apply G→A, then RC
            r1_seq_buf.extend_from_slice(r1_umi);
            let mut r1_template = template[r2_start..r2_end].to_vec();
            apply_methylation_conversion(
                &mut r1_template,
                &template,
                r2_start,
                false, // bottom strand
                &params.methylation,
                &mut rng,
            );
            reverse_complement_into(&r1_template, &mut r1_seq_buf);

            // B R2 reads the top strand (template start, forward) — apply C→T
            r2_seq_buf.extend_from_slice(r2_umi);
            let r1_template_end = template_len.min(template.len());
            let mut r2_template = template[..r1_template_end].to_vec();
            apply_methylation_conversion(
                &mut r2_template,
                &template,
                0,
                true, // top strand
                &params.methylation,
                &mut rng,
            );
            r2_seq_buf.extend_from_slice(&r2_template);
        }

        // Introduce UMI errors directly into buffer
        introduce_errors_inplace(&mut r1_seq_buf[..params.umi_length], 0.01, &mut rng);

        // Pad R1 to read length if needed
        pad_sequence_inplace(&mut r1_seq_buf, params.read_length, &mut rng);

        // Introduce UMI errors in R2 as well
        introduce_errors_inplace(&mut r2_seq_buf[..params.umi_length], 0.01, &mut rng);

        // Pad R2 to read length if needed
        pad_sequence_inplace(&mut r2_seq_buf, params.read_length, &mut rng);

        // Generate quality scores
        let r1_quals = quality_model.generate_qualities(r1_seq_buf.len(), &mut rng);
        let r2_quals_raw = quality_model.generate_qualities(r2_seq_buf.len(), &mut rng);

        // Apply R2 quality bias only if offset is non-zero
        let r2_quals = if params.r2_quality_offset != 0 {
            quality_bias.apply_to_vec(&r2_quals_raw, true)
        } else {
            r2_quals_raw
        };

        records.push(ReadRecord {
            read_name,
            r1_seq: r1_seq_buf.clone(),
            r1_quals,
            r2_seq: r2_seq_buf.clone(),
            r2_quals,
            umi: combined_umi.clone(),
            mol_id,
            family_id: family_idx * 1000 + read_idx,
            strand,
            chrom: chrom.clone(),
            pos,
        });
    }

    records
}

/// Reverse complement a DNA sequence into an existing buffer.
fn reverse_complement_into(seq: &[u8], output: &mut Vec<u8>) {
    output.reserve(seq.len());
    for &b in seq.iter().rev() {
        output.push(complement_base(b));
    }
}

/// Introduce random substitution errors in-place.
/// Uses direct lookup tables to avoid retry loops.
fn introduce_errors_inplace(seq: &mut [u8], error_rate: f64, rng: &mut impl Rng) {
    // Lookup table: for each base, provides 3 alternatives
    // A(65) -> C, G, T; C(67) -> A, G, T; G(71) -> A, C, T; T(84) -> A, C, G
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

/// Pad sequence in-place to target length with random bases, or truncate if too long.
fn pad_sequence_inplace(seq: &mut Vec<u8>, target_len: usize, rng: &mut impl Rng) {
    const BASES: &[u8] = b"ACGT";
    seq.reserve(target_len.saturating_sub(seq.len()));
    while seq.len() < target_len {
        seq.push(BASES[rng.random_range(0..4)]);
    }
    seq.truncate(target_len);
}

#[cfg(test)]
#[allow(clippy::naive_bytecount)]
mod tests {
    use super::*;
    use crate::simulate::create_rng;

    /// Reverse complement a DNA sequence (test wrapper).
    fn reverse_complement(seq: &[u8]) -> Vec<u8> {
        let mut result = Vec::with_capacity(seq.len());
        reverse_complement_into(seq, &mut result);
        result
    }

    /// Introduce random substitution errors (test wrapper).
    fn introduce_errors(seq: &[u8], error_rate: f64, rng: &mut impl Rng) -> Vec<u8> {
        let mut result = seq.to_vec();
        introduce_errors_inplace(&mut result, error_rate, rng);
        result
    }

    /// Pad sequence to target length with random bases (test wrapper).
    fn pad_sequence(mut seq: Vec<u8>, target_len: usize, rng: &mut impl Rng) -> Vec<u8> {
        pad_sequence_inplace(&mut seq, target_len, rng);
        seq
    }

    #[test]
    fn test_generate_random_sequence_length() {
        let mut rng = create_rng(Some(42));
        for len in [0, 1, 10, 100, 500] {
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
                "Invalid base: {base}"
            );
        }
    }

    #[test]
    fn test_generate_random_sequence_distribution() {
        let mut rng = create_rng(Some(42));
        let seq = generate_random_sequence(10000, &mut rng);

        let count_a = seq.iter().filter(|&&b| b == b'A').count();
        let count_c = seq.iter().filter(|&&b| b == b'C').count();
        let count_g = seq.iter().filter(|&&b| b == b'G').count();
        let count_t = seq.iter().filter(|&&b| b == b'T').count();

        // Each base should be roughly 25% (within 5%)
        for (name, count) in [("A", count_a), ("C", count_c), ("G", count_g), ("T", count_t)] {
            let frac = count as f64 / 10000.0;
            assert!((0.20..0.30).contains(&frac), "{name} fraction {frac} out of expected range");
        }
    }

    #[test]
    fn test_generate_random_sequence_reproducible() {
        let mut rng1 = create_rng(Some(42));
        let mut rng2 = create_rng(Some(42));

        let seq1 = generate_random_sequence(100, &mut rng1);
        let seq2 = generate_random_sequence(100, &mut rng2);

        assert_eq!(seq1, seq2);
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
        // ACGT -> complement TGCA -> reverse ACGT
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");

        // AAAT -> complement TTTA -> reverse ATTT
        assert_eq!(reverse_complement(b"AAAT"), b"ATTT");

        // GGGC -> complement CCCG -> reverse GCCC
        assert_eq!(reverse_complement(b"GGGC"), b"GCCC");
    }

    #[test]
    fn test_reverse_complement_empty() {
        assert_eq!(reverse_complement(b""), Vec::<u8>::new());
    }

    #[test]
    fn test_reverse_complement_single() {
        for (base, expected) in [(b'A', b'T'), (b'T', b'A'), (b'C', b'G'), (b'G', b'C')] {
            assert_eq!(reverse_complement(&[base]), vec![expected]);
        }
    }

    #[test]
    fn test_reverse_complement_unknown_bases() {
        // Unknown bases are returned unchanged by the library complement_base
        assert_eq!(reverse_complement(b"N"), b"N");
        assert_eq!(reverse_complement(b"X"), b"X");
        assert_eq!(reverse_complement(b"ANT"), b"ANT"); // A->T, N->N, T->A, reversed = ANT
    }

    #[test]
    fn test_reverse_complement_double() {
        let mut rng = create_rng(Some(42));
        let seq = generate_random_sequence(100, &mut rng);

        // Double reverse complement should give back original
        let rc = reverse_complement(&seq);
        let rc_rc = reverse_complement(&rc);

        assert_eq!(seq, rc_rc);
    }

    #[test]
    fn test_introduce_errors_zero_rate() {
        let mut rng = create_rng(Some(42));
        let seq = b"ACGTACGTACGT";

        let result = introduce_errors(seq, 0.0, &mut rng);
        assert_eq!(result, seq.to_vec());
    }

    #[test]
    fn test_introduce_errors_full_rate() {
        let mut rng = create_rng(Some(42));
        let seq = b"AAAAAAAAAA"; // 10 A's

        let result = introduce_errors(seq, 1.0, &mut rng);

        // All bases should be different (not A)
        for &base in &result {
            assert_ne!(base, b'A', "Expected all bases to be mutated");
        }
    }

    #[test]
    fn test_introduce_errors_valid_bases() {
        let mut rng = create_rng(Some(42));
        let seq = generate_random_sequence(1000, &mut rng);

        let result = introduce_errors(&seq, 0.5, &mut rng);

        // All bases should still be valid
        for &base in &result {
            assert!(
                base == b'A' || base == b'C' || base == b'G' || base == b'T',
                "Invalid base after mutation: {base}"
            );
        }
    }

    #[test]
    fn test_introduce_errors_length_preserved() {
        let mut rng = create_rng(Some(42));
        let seq = generate_random_sequence(100, &mut rng);

        let result = introduce_errors(&seq, 0.1, &mut rng);
        assert_eq!(result.len(), seq.len());
    }

    #[test]
    fn test_introduce_errors_approximate_rate() {
        let mut rng = create_rng(Some(42));
        let seq = vec![b'A'; 10000]; // All A's so we can count errors easily

        let result = introduce_errors(&seq, 0.1, &mut rng);
        let error_count = result.iter().filter(|&&b| b != b'A').count();

        // Should be roughly 10% errors (within 2%)
        let error_rate = error_count as f64 / 10000.0;
        assert!(
            (error_rate - 0.1).abs() < 0.02,
            "Error rate {error_rate} not close to expected 0.1"
        );
    }

    #[test]
    fn test_pad_sequence_already_correct_length() {
        let mut rng = create_rng(Some(42));
        let seq = b"ACGTACGT".to_vec();

        let result = pad_sequence(seq.clone(), 8, &mut rng);
        assert_eq!(result, seq);
    }

    #[test]
    fn test_pad_sequence_needs_padding() {
        let mut rng = create_rng(Some(42));
        let seq = b"ACGT".to_vec();

        let result = pad_sequence(seq, 10, &mut rng);
        assert_eq!(result.len(), 10);

        // First 4 bases should be original
        assert_eq!(&result[..4], b"ACGT");

        // Remaining bases should be valid DNA
        for &base in &result[4..] {
            assert!(base == b'A' || base == b'C' || base == b'G' || base == b'T');
        }
    }

    #[test]
    fn test_pad_sequence_needs_truncation() {
        let mut rng = create_rng(Some(42));
        let seq = b"ACGTACGTACGT".to_vec();

        let result = pad_sequence(seq, 5, &mut rng);
        assert_eq!(result.len(), 5);
        assert_eq!(result, b"ACGTA".to_vec());
    }

    #[test]
    fn test_pad_sequence_empty_to_length() {
        let mut rng = create_rng(Some(42));
        let seq = Vec::new();

        let result = pad_sequence(seq, 10, &mut rng);
        assert_eq!(result.len(), 10);

        // All bases should be valid
        for &base in &result {
            assert!(base == b'A' || base == b'C' || base == b'G' || base == b'T');
        }
    }

    #[test]
    fn test_pad_sequence_to_zero() {
        let mut rng = create_rng(Some(42));
        let seq = b"ACGT".to_vec();

        let result = pad_sequence(seq, 0, &mut rng);
        assert!(result.is_empty());
    }

    #[test]
    fn test_load_includelist_valid() {
        let dir = tempfile::tempdir().expect("failed to create temp dir");
        let path = dir.path().join("umis.txt");
        std::fs::write(&path, "AACAC\nAAGGA\nAATGC\n").expect("failed to write file");

        let umis = load_includelist(&path).expect("failed to load includelist");
        assert_eq!(umis.len(), 3);
        assert_eq!(umis[0], b"AACAC");
        assert_eq!(umis[1], b"AAGGA");
        assert_eq!(umis[2], b"AATGC");
    }

    #[test]
    fn test_load_includelist_lowercase_uppercased() {
        let dir = tempfile::tempdir().expect("failed to create temp dir");
        let path = dir.path().join("umis.txt");
        std::fs::write(&path, "aacac\naagga\n").expect("failed to write file");

        let umis = load_includelist(&path).expect("failed to load includelist");
        assert_eq!(umis[0], b"AACAC");
        assert_eq!(umis[1], b"AAGGA");
    }

    #[test]
    fn test_load_includelist_skips_blank_lines() {
        let dir = tempfile::tempdir().expect("failed to create temp dir");
        let path = dir.path().join("umis.txt");
        std::fs::write(&path, "AACAC\n\nAAGGA\n  \nAATGC\n").expect("failed to write file");

        let umis = load_includelist(&path).expect("failed to load includelist");
        assert_eq!(umis.len(), 3);
    }

    #[test]
    fn test_load_includelist_rejects_invalid_bases() {
        let dir = tempfile::tempdir().expect("failed to create temp dir");
        let path = dir.path().join("umis.txt");
        std::fs::write(&path, "AACAC\nAANGA\n").expect("failed to write file");

        let result = load_includelist(&path);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("Invalid UMI"));
    }

    #[test]
    fn test_load_includelist_rejects_mismatched_lengths() {
        let dir = tempfile::tempdir().expect("failed to create temp dir");
        let path = dir.path().join("umis.txt");
        std::fs::write(&path, "AACAC\nAAGG\n").expect("failed to write file");

        let result = load_includelist(&path);
        assert!(result.is_err());
        let msg = result.unwrap_err().to_string();
        assert!(msg.contains("length mismatch"));
        // Should report the actual file line number (line 2), not a filtered index
        assert!(msg.contains("line 2"), "Expected file line number in error: {msg}");
    }

    #[test]
    fn test_load_includelist_mismatched_length_reports_correct_line_with_blanks() {
        let dir = tempfile::tempdir().expect("failed to create temp dir");
        let path = dir.path().join("umis.txt");
        // Blank line between valid UMIs; the short UMI is on file line 4
        std::fs::write(&path, "AACAC\n\nAAGGA\nAAGG\n").expect("failed to write file");

        let result = load_includelist(&path);
        assert!(result.is_err());
        let msg = result.unwrap_err().to_string();
        assert!(msg.contains("line 4"), "Expected file line 4 in error: {msg}");
    }

    #[test]
    fn test_load_includelist_rejects_empty() {
        let dir = tempfile::tempdir().expect("failed to create temp dir");
        let path = dir.path().join("umis.txt");
        std::fs::write(&path, "\n\n").expect("failed to write file");

        let result = load_includelist(&path);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("empty"));
    }

    #[test]
    fn test_execute_rejects_includelist_umi_exceeding_read_length() {
        let dir = tempfile::tempdir().expect("failed to create temp dir");
        let includelist_path = dir.path().join("umis.txt");
        std::fs::write(&includelist_path, "AAAAACCCCC\n").expect("failed to write file"); // length 10

        let r1 = dir.path().join("r1.fq.gz");
        let r2 = dir.path().join("r2.fq.gz");
        let truth = dir.path().join("truth.tsv");

        let cmd = FastqReads::try_parse_from([
            "fastq-reads",
            "-1",
            r1.to_str().expect("path should be valid UTF-8"),
            "-2",
            r2.to_str().expect("path should be valid UTF-8"),
            "--truth",
            truth.to_str().expect("path should be valid UTF-8"),
            "-r",
            "dummy.fa",
            "-i",
            includelist_path.to_str().expect("path should be valid UTF-8"),
            "--read-length",
            "5", // shorter than UMI
        ])
        .expect("valid CLI arguments");

        let result = cmd.execute("");
        assert!(result.is_err());
        let msg = result.unwrap_err().to_string();
        assert!(
            msg.contains("exceeds read length"),
            "Expected 'exceeds read length' in error: {msg}",
        );
    }

    #[test]
    fn test_execute_rejects_umi_length_exceeding_read_length() {
        let dir = tempfile::tempdir().expect("failed to create temp dir");
        let r1 = dir.path().join("r1.fq.gz");
        let r2 = dir.path().join("r2.fq.gz");
        let truth = dir.path().join("truth.tsv");

        let cmd = FastqReads::try_parse_from([
            "fastq-reads",
            "-1",
            r1.to_str().expect("path should be valid UTF-8"),
            "-2",
            r2.to_str().expect("path should be valid UTF-8"),
            "--truth",
            truth.to_str().expect("path should be valid UTF-8"),
            "-r",
            "dummy.fa",
            "--umi-length",
            "100",
            "--read-length",
            "50",
        ])
        .expect("valid CLI arguments");

        let result = cmd.execute("");
        assert!(result.is_err());
        let msg = result.unwrap_err().to_string();
        assert!(
            msg.contains("exceeds read length"),
            "Expected 'exceeds read length' in error: {msg}",
        );
    }

    #[test]
    fn test_generate_molecule_reads_with_includelist() {
        use std::io::Write as IoWrite;
        use tempfile::NamedTempFile;

        let mut fasta = NamedTempFile::new().unwrap();
        writeln!(fasta, ">chr1").unwrap();
        fasta.write_all(&b"ACGT".repeat(500)).unwrap();
        writeln!(fasta).unwrap();
        fasta.flush().unwrap();
        let ref_genome = ReferenceGenome::load(fasta.path()).unwrap();

        let umis = vec![b"AAAAA".to_vec(), b"CCCCC".to_vec(), b"GGGGG".to_vec()];
        let params = GenerationParams {
            umi_length: 5,
            read_length: 50,
            min_family_size: 1,
            r2_quality_offset: 0,
            duplex: false,
            includelist: Some(umis.clone()),
            methylation: MethylationConfig {
                mode: fgumi_consensus::MethylationMode::Disabled,
                cpg_methylation_rate: 0.75,
                conversion_rate: 0.98,
            },
        };

        let quality_model =
            crate::simulate::PositionQualityModel::new(10, 25, 37, 100, 0.08, 2, 2.0);
        let quality_bias = crate::simulate::ReadPairQualityBias::new(0);
        let family_dist = crate::simulate::FamilySizeDistribution::log_normal(3.0, 1.0);
        let insert_model = crate::simulate::InsertSizeModel::new(150.0, 30.0, 50, 500);

        // Generate molecules and verify UMIs come from includelist
        for seed in 0..20u64 {
            let records = generate_molecule_reads(
                0,
                seed,
                &params,
                &quality_model,
                &quality_bias,
                &family_dist,
                &insert_model,
                &ref_genome,
            );
            for record in &records {
                // The truth UMI (before errors) should be from the includelist
                // UMI format is "R1-R2", split and check each half
                let umi_str = String::from_utf8_lossy(&record.umi);
                let parts: Vec<&str> = umi_str.split('-').collect();
                assert_eq!(parts.len(), 2);
                // Note: the actual read sequences have errors introduced, but the
                // truth UMI in the record should be from the includelist
                assert!(
                    umis.iter().any(|u| u == parts[0].as_bytes()),
                    "R1 UMI '{}' not in includelist",
                    parts[0],
                );
                assert!(
                    umis.iter().any(|u| u == parts[1].as_bytes()),
                    "R2 UMI '{}' not in includelist",
                    parts[1],
                );
            }
        }
    }

    #[test]
    fn test_emseq_fastq_reads_converts_non_cpg_c() {
        // With EM-Seq, non-CpG C should be converted to T (on top strand R1)
        // Use a reference that has a known pattern: all C's with no adjacent G
        // Template: "CACACACACA..." (no CpG dinucleotides)
        use crate::commands::simulate::common::apply_methylation_conversion;

        let template = b"CACACACACACACACACACAC".to_vec();
        let mut r1 = template.clone();
        let mut rng = create_rng(Some(42));

        let config = MethylationConfig {
            mode: fgumi_consensus::MethylationMode::EmSeq,
            cpg_methylation_rate: 0.75,
            conversion_rate: 1.0,
        };
        apply_methylation_conversion(&mut r1, &template, 0, true, &config, &mut rng);

        // All C's should be converted to T (non-CpG in EM-Seq)
        for (i, &b) in r1.iter().enumerate() {
            if template[i] == b'C' {
                assert_eq!(b, b'T', "position {i}: non-CpG C should be converted to T in EM-Seq");
            } else {
                assert_eq!(b, template[i], "position {i}: A should remain unchanged");
            }
        }
    }

    #[test]
    fn test_taps_fastq_reads_preserves_non_cpg_c() {
        // With TAPs, non-CpG C should NOT be converted
        use crate::commands::simulate::common::apply_methylation_conversion;

        let template = b"CACACACACACACACACACAC".to_vec();
        let mut r1 = template.clone();
        let mut rng = create_rng(Some(42));

        let config = MethylationConfig {
            mode: fgumi_consensus::MethylationMode::Taps,
            cpg_methylation_rate: 0.75,
            conversion_rate: 1.0,
        };
        apply_methylation_conversion(&mut r1, &template, 0, true, &config, &mut rng);

        // No changes — non-CpG C is not a target in TAPs
        assert_eq!(r1, template);
    }

    #[test]
    fn test_reads_in_same_family_differ_with_methylation() {
        // Each read independently samples conversion, so two reads from
        // the same molecule should (usually) differ when rates are partial
        use crate::commands::simulate::common::apply_methylation_conversion;

        // Template with CpG sites
        let template = b"ACGTCGATCGACGTCGATCG".to_vec();
        let mut different_count = 0;

        let config = MethylationConfig {
            mode: fgumi_consensus::MethylationMode::EmSeq,
            cpg_methylation_rate: 0.5,
            conversion_rate: 0.5,
        };

        for seed in 0..50u64 {
            let mut r1_a = template.clone();
            let mut r1_b = template.clone();
            let mut rng_a = create_rng(Some(seed * 2));
            let mut rng_b = create_rng(Some(seed * 2 + 1));

            apply_methylation_conversion(&mut r1_a, &template, 0, true, &config, &mut rng_a);
            apply_methylation_conversion(&mut r1_b, &template, 0, true, &config, &mut rng_b);

            if r1_a != r1_b {
                different_count += 1;
            }
        }

        // With 50% CpG methylation rate, most pairs should differ
        assert!(different_count > 10, "Expected most read pairs to differ, got {different_count}");
    }
}

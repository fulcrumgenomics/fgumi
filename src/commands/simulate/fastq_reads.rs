//! Generate paired-end FASTQ files with UMI sequences.

use crate::commands::command::Command;
use crate::commands::simulate::common::{
    FamilySizeArgs, InsertSizeArgs, QualityArgs, SimulationCommon,
};
use anyhow::{Context, Result, bail};
use clap::Parser;
use crossbeam_channel::bounded;
use fgumi_lib::progress::ProgressTracker;
use fgumi_lib::simulate::{FastqWriter, create_rng};
use log::info;
use noodles::fasta;
use rand::{Rng, RngExt};
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::PathBuf;
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
    #[arg(long = "duplex")]
    pub duplex: bool,

    /// Reference FASTA file for sampling template sequences.
    /// If provided, templates are sampled from this reference instead of random bases.
    /// This produces reads that will map back to the reference.
    #[arg(short = 'r', long = "reference")]
    pub reference: Option<PathBuf>,

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
    /// Chromosome name (if reference-based)
    chrom: Option<String>,
    /// 1-based start position (if reference-based)
    pos: Option<usize>,
}

/// Loaded reference genome for sampling template sequences.
struct ReferenceGenome {
    /// Chromosome names
    names: Vec<String>,
    /// Chromosome sequences (uppercase)
    sequences: Vec<Vec<u8>>,
    /// Cumulative lengths for weighted random selection
    cumulative_lengths: Vec<usize>,
    /// Total genome length
    total_length: usize,
}

impl ReferenceGenome {
    /// Load a reference genome from a FASTA file.
    fn load(path: &PathBuf) -> Result<Self> {
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

            if seq.len() >= 1000 {
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

    /// Sample a random position and return (`chrom_idx`, position, sequence).
    /// Returns None if the position contains N bases.
    fn sample_sequence(
        &self,
        length: usize,
        rng: &mut impl Rng,
    ) -> Option<(usize, usize, Vec<u8>)> {
        // Try up to 10 times to find a valid position without N bases
        for _ in 0..10 {
            // Pick a random position in the genome
            let genome_pos = rng.random_range(0..self.total_length.saturating_sub(length));

            // Find which chromosome this falls in
            let chrom_idx =
                self.cumulative_lengths.iter().position(|&cum| cum > genome_pos).unwrap_or(0);

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
}

/// Parameters needed for parallel molecule generation.
#[derive(Clone)]
struct GenerationParams {
    umi_length: usize,
    read_length: usize,
    duplex: bool,
    min_family_size: usize,
    r2_quality_offset: i8,
}

/// Channel capacity for buffering molecule batches between producer and writer threads.
/// Each batch contains all reads for one molecule (typically 1-10 reads).
const CHANNEL_CAPACITY: usize = 1_000;

impl Command for FastqReads {
    fn execute(&self, _command_line: &str) -> Result<()> {
        info!("Generating FASTQ reads");
        info!("  Output R1: {}", self.r1_output.display());
        info!("  Output R2: {}", self.r2_output.display());
        info!("  Truth: {}", self.truth_output.display());
        info!("  Num molecules: {}", self.common.num_molecules);
        info!("  Read length: {}", self.common.read_length);
        info!("  UMI length: {}", self.common.umi_length);
        info!("  Duplex: {}", self.duplex);
        info!("  Threads: {}", self.threads);
        if let Some(ref path) = self.reference {
            info!("  Reference: {}", path.display());
        }

        // Load reference genome if provided
        let reference = if let Some(ref path) = self.reference {
            Some(Arc::new(ReferenceGenome::load(path)?))
        } else {
            None
        };

        // Use at least 2 threads for generation if multi-threaded
        // Reserve some threads for gzip compression
        let gen_threads = if self.threads <= 1 { 1 } else { self.threads.max(2) };
        let compress_threads = self.threads;

        let quality_model = Arc::new(self.quality.to_quality_model());
        let quality_bias = Arc::new(self.quality.to_quality_bias());
        let family_dist = Arc::new(self.family_size.to_family_size_distribution()?);
        let insert_model = Arc::new(self.insert_size.to_insert_size_model());

        let params = Arc::new(GenerationParams {
            umi_length: self.common.umi_length,
            read_length: self.common.read_length,
            duplex: self.duplex,
            min_family_size: self.family_size.min_family_size,
            r2_quality_offset: self.quality.r2_quality_offset,
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
        let has_reference = reference.is_some();

        // Spawn writer thread with multi-threaded gzip compression
        let writer_handle = thread::spawn(move || -> Result<u64> {
            let mut r1_writer = FastqWriter::with_threads(&r1_path, compress_threads)?;
            let mut r2_writer = FastqWriter::with_threads(&r2_path, compress_threads)?;
            let truth_file = File::create(&truth_path)
                .with_context(|| format!("Failed to create {}", truth_path.display()))?;
            let mut truth_writer = BufWriter::new(truth_file);

            // Write truth header (include chrom/pos if reference-based)
            if has_reference {
                writeln!(
                    truth_writer,
                    "read_name\ttrue_umi\tmolecule_id\tfamily_id\tstrand\tchrom\tpos"
                )?;
            } else {
                writeln!(truth_writer, "read_name\ttrue_umi\tmolecule_id\tfamily_id\tstrand")?;
            }

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
                    if has_reference {
                        writeln!(
                            truth_writer,
                            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                            record.read_name,
                            String::from_utf8_lossy(&record.umi),
                            record.mol_id,
                            record.family_id,
                            record.strand,
                            record.chrom.as_deref().unwrap_or("."),
                            record.pos.map_or_else(|| ".".to_string(), |p| (p + 1).to_string())
                        )?;
                    } else {
                        writeln!(
                            truth_writer,
                            "{}\t{}\t{}\t{}\t{}",
                            record.read_name,
                            String::from_utf8_lossy(&record.umi),
                            record.mol_id,
                            record.family_id,
                            record.strand
                        )?;
                    }
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
                        reference.as_deref(),
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
    quality_model: &fgumi_lib::simulate::PositionQualityModel,
    quality_bias: &fgumi_lib::simulate::ReadPairQualityBias,
    family_dist: &fgumi_lib::simulate::FamilySizeDistribution,
    insert_model: &fgumi_lib::simulate::InsertSizeModel,
    reference: Option<&ReferenceGenome>,
) -> Vec<ReadRecord> {
    let mut rng = create_rng(Some(seed));
    let mut records = Vec::new();

    // Generate UMIs for this molecule (one for R1, one for R2)
    // For duplex sequencing with paired grouping strategy, B strand needs swapped UMIs
    let umi_r1 = generate_random_sequence(params.umi_length, &mut rng);
    let umi_r2 = generate_random_sequence(params.umi_length, &mut rng);

    // Generate family size
    let family_size = family_dist.sample(&mut rng, params.min_family_size);

    // Generate insert size for this molecule
    let insert_size = insert_model.sample(&mut rng);

    // Generate template sequence (shared by all reads in family)
    // Either sample from reference or generate random
    let template_len = insert_size.saturating_sub(params.umi_length);
    let (template, chrom, pos) = if let Some(ref_genome) = reference {
        // Sample from reference
        if let Some((chrom_idx, local_pos, seq)) =
            ref_genome.sample_sequence(template_len, &mut rng)
        {
            let chrom_name = ref_genome.names[chrom_idx].clone();
            (seq, Some(chrom_name), Some(local_pos))
        } else {
            // Fallback to random if sampling failed (e.g., all N regions)
            (generate_random_sequence(template_len, &mut rng), None, None)
        }
    } else {
        // Generate random template
        (generate_random_sequence(template_len, &mut rng), None, None)
    };

    // Pre-calculate template slicing bounds (both reads have UMI prefix)
    let r1_template_len = params.read_length.saturating_sub(params.umi_length);
    let r2_template_len = params.read_length.saturating_sub(params.umi_length);
    let r2_start = template.len().saturating_sub(r2_template_len);
    let r2_end = template.len();

    // Reusable buffers for sequence building
    let mut r1_seq_buf = Vec::with_capacity(params.read_length);
    let mut r2_seq_buf = Vec::with_capacity(params.read_length);

    let mut read_counter = 0usize;

    // For duplex mode, each read pair is independently assigned to A or B strand (coin flip)
    // For non-duplex mode, all reads are A strand
    for read_idx in 0..family_size {
        // Determine strand for this read pair
        let (strand, family_idx): (&'static str, usize) = if params.duplex {
            // 50/50 coin flip for A or B strand
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
            // R1 = UMI + template start (maps forward)
            r1_seq_buf.extend_from_slice(r1_umi);
            let r1_template_end = r1_template_len.min(template.len());
            r1_seq_buf.extend_from_slice(&template[..r1_template_end]);

            // R2 = UMI + revcomp(template end) (maps reverse)
            r2_seq_buf.extend_from_slice(r2_umi);
            reverse_complement_into(&template[r2_start..r2_end], &mut r2_seq_buf);
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

            // B R1 = UMI + revcomp(template end) - same region as A R2
            r1_seq_buf.extend_from_slice(r1_umi);
            reverse_complement_into(&template[r2_start..r2_end], &mut r1_seq_buf);

            // B R2 = UMI + template start - same region as A R1
            r2_seq_buf.extend_from_slice(r2_umi);
            let r1_template_end = r1_template_len.min(template.len());
            r2_seq_buf.extend_from_slice(&template[..r1_template_end]);
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

/// Generate a random DNA sequence.
fn generate_random_sequence(len: usize, rng: &mut impl Rng) -> Vec<u8> {
    const BASES: &[u8] = b"ACGT";
    let mut seq = Vec::with_capacity(len);
    for _ in 0..len {
        seq.push(BASES[rng.random_range(0..4)]);
    }
    seq
}

/// Reverse complement a DNA sequence into an existing buffer.
fn reverse_complement_into(seq: &[u8], output: &mut Vec<u8>) {
    output.reserve(seq.len());
    for &b in seq.iter().rev() {
        output.push(match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'N',
        });
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
    use fgumi_lib::simulate::create_rng;

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
        // Unknown bases should become N
        assert_eq!(reverse_complement(b"N"), b"N");
        assert_eq!(reverse_complement(b"X"), b"N");
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
}

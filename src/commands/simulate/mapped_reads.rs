//! Generate template-coordinate sorted BAM for group.
//!
//! Uses a memory-efficient approach: sort molecule IDs by their position first,
//! then generate records in sorted order (streaming). This avoids storing all
//! records in memory before sorting.

use super::sort::TemplateCoordKey;
use crate::commands::command::Command;
use crate::commands::common::CompressionOptions;
use crate::commands::simulate::common::{
    FamilySizeArgs, InsertSizeArgs, PositionDistArgs, QualityArgs, ReferenceArgs, SimulationCommon,
};
use anyhow::{Context, Result};
use bstr::BString;
use clap::Parser;
use fgumi_lib::bam_io::create_bam_writer;
use fgumi_lib::progress::ProgressTracker;
use fgumi_lib::sam::builder::RecordBuilder;
use fgumi_lib::simulate::{
    FamilySizeDistribution, InsertSizeModel, PositionQualityModel, ReadPairQualityBias, create_rng,
};
use log::info;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::sam::header::Header;
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::ReferenceSequence;
use noodles::sam::header::record::value::map::header::{self as HeaderRecord, Tag as HeaderTag};
use rand::{Rng, RngExt};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::num::NonZeroUsize;
use std::path::PathBuf;

/// Generate template-coordinate sorted BAM with paired alignments for `fgumi group`.
#[derive(Parser, Debug)]
#[command(
    name = "mapped-reads",
    about = "Generate template-coord sorted BAM for group",
    long_about = r#"
Generate synthetic paired-end BAM files with proper alignments.

The output is template-coordinate sorted and suitable for input to `fgumi group`.
Reads contain RX tags with UMI sequences.
"#
)]
pub struct MappedReads {
    /// Output BAM file (template-coordinate sorted)
    #[arg(short = 'o', long = "output", required = true)]
    pub output: PathBuf,

    /// Output truth TSV file for validation
    #[arg(long = "truth", required = true)]
    pub truth_output: PathBuf,

    /// Mapping quality for aligned reads
    #[arg(long = "mapq", default_value = "60")]
    pub mapq: u8,

    /// Fraction of reads to leave unmapped
    #[arg(long = "unmapped-fraction", default_value = "0.0")]
    pub unmapped_fraction: f64,

    /// Number of writer threads
    #[arg(short = 't', long = "threads", default_value = "1")]
    pub threads: usize,

    /// Compression options for output BAM.
    #[command(flatten)]
    pub compression: CompressionOptions,

    #[command(flatten)]
    pub common: SimulationCommon,

    #[command(flatten)]
    pub quality: QualityArgs,

    #[command(flatten)]
    pub family_size: FamilySizeArgs,

    #[command(flatten)]
    pub insert_size: InsertSizeArgs,

    #[command(flatten)]
    pub reference: ReferenceArgs,

    #[command(flatten)]
    pub position_dist: PositionDistArgs,
}

/// Lightweight molecule info for position-first sorting.
struct MoleculeInfo {
    mol_id: usize,
    seed: u64,
    sort_key: TemplateCoordKey,
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

/// Parameters needed for molecule generation.
struct GenerationParams {
    read_length: usize,
    umi_length: usize,
    mapq: u8,
    min_family_size: usize,
    num_positions: usize,
    ref_length: usize,
    quality_model: PositionQualityModel,
    quality_bias: ReadPairQualityBias,
    family_dist: FamilySizeDistribution,
    insert_model: InsertSizeModel,
}

impl Command for MappedReads {
    fn execute(&self, command_line: &str) -> Result<()> {
        info!("Generating mapped reads");
        info!("  Output: {}", self.output.display());
        info!("  Truth: {}", self.truth_output.display());
        info!("  Num molecules: {}", self.common.num_molecules);
        info!("  Read length: {}", self.common.read_length);
        info!("  UMI length: {}", self.common.umi_length);
        info!("  Threads: {}", self.threads);

        // Determine positions
        let num_positions = self.position_dist.num_positions.unwrap_or(self.common.num_molecules);

        // Build header with template-coordinate sort order
        let ref_name = self.reference.ref_name.clone();
        let ref_length = self.reference.ref_length;
        let mut header = Header::builder();

        // Add sort order tags: SO:unsorted, GO:query, SS:template-coordinate
        let HeaderTag::Other(so_tag) = HeaderTag::from([b'S', b'O']) else { unreachable!() };
        let HeaderTag::Other(go_tag) = HeaderTag::from([b'G', b'O']) else { unreachable!() };
        let HeaderTag::Other(ss_tag) = HeaderTag::from([b'S', b'S']) else { unreachable!() };

        let header_map = Map::<HeaderRecord::Header>::builder()
            .insert(so_tag, "unsorted")
            .insert(go_tag, "query")
            .insert(ss_tag, "template-coordinate")
            .build()
            .unwrap();
        header = header.set_header(header_map);

        let length = NonZeroUsize::try_from(ref_length).expect("Reference length must be > 0");
        let ref_seq: Map<ReferenceSequence> = Map::<ReferenceSequence>::new(length);
        header = header.add_reference_sequence(BString::from(&*ref_name), ref_seq);

        // Add @PG record
        header = fgumi_lib::header::add_pg_to_builder(
            header,
            crate::version::VERSION.as_str(),
            command_line,
        )?;

        let header = header.build();

        // Set up generation parameters
        let params = GenerationParams {
            read_length: self.common.read_length,
            umi_length: self.common.umi_length,
            mapq: self.mapq,
            min_family_size: self.family_size.min_family_size,
            num_positions,
            ref_length,
            quality_model: self.quality.to_quality_model(),
            quality_bias: self.quality.to_quality_bias(),
            family_dist: self.family_size.to_family_size_distribution()?,
            insert_model: self.insert_size.to_insert_size_model(),
        };

        // Generate seeds for reproducibility
        let mut seed_rng = create_rng(self.common.seed);

        // MEMORY-EFFICIENT APPROACH: Sort molecule IDs by template-coordinate key first,
        // then generate records in sorted order (streaming).
        //
        // We need to pre-compute insert_size to get pos2 for the sort key, since
        // template-coordinate sorting considers both ends of the template.
        info!("Computing molecule positions and sort keys...");
        let mut molecules: Vec<MoleculeInfo> = (0..self.common.num_molecules)
            .map(|mol_id| {
                let seed: u64 = seed_rng.random();
                let pos1 = compute_position(mol_id, num_positions, ref_length);

                // Pre-compute insert_size using the molecule's seed (same RNG sequence as generation)
                let mut mol_rng = create_rng(Some(seed));
                // Skip UMI generation RNG calls
                for _ in 0..params.umi_length {
                    let _: usize = mol_rng.random_range(0..4);
                }
                // Skip family_size RNG call
                let _ = params.family_dist.sample(&mut mol_rng, params.min_family_size);
                // Get insert_size
                let insert_size = params.insert_model.sample(&mut mol_rng);

                // Build template-coordinate sort key using the shared sort module.
                // for_f1r2_pair handles pos2 calculation (unclipped 5' of R2 = pos1 + insert - 1)
                let sort_key = TemplateCoordKey::for_f1r2_pair(
                    0, // tid - all simulated reads on reference 0
                    pos1,
                    insert_size,
                    String::new(), // empty mid for mapped-reads (no MI tag yet)
                    format!("mol{mol_id:08}"),
                );
                MoleculeInfo { mol_id, seed, sort_key }
            })
            .collect();

        // Sort molecules by template-coordinate key (matching samtools bam1_cmp_template_coordinate)
        info!("Sorting {} molecules by template-coordinate...", molecules.len());
        molecules.sort_unstable();

        // Set up writer with multi-threaded BGZF compression
        let mut writer = create_bam_writer(
            &self.output,
            &header,
            self.threads,
            self.compression.compression_level,
        )?;

        // Create truth file
        let truth_file = File::create(&self.truth_output)
            .with_context(|| format!("Failed to create {}", self.truth_output.display()))?;
        let mut truth_writer = BufWriter::new(truth_file);
        writeln!(truth_writer, "read_name\ttrue_umi\tmolecule_id\tchrom\tposition\tstrand")?;

        // Generate and write records in position-sorted order (streaming)
        info!("Generating and writing records in sorted order...");
        let progress = ProgressTracker::new("Processed molecules").with_interval(100_000);
        let mut total_pairs = 0usize;

        for mol_info in molecules {
            progress.log_if_needed(1);

            // Generate all read pairs for this molecule
            let pairs = generate_molecule_reads(mol_info.mol_id, mol_info.seed, &params);

            for (r1, r2, read_name, umi_str) in pairs {
                writer.write_alignment_record(&header, &r1)?;
                writer.write_alignment_record(&header, &r2)?;
                writeln!(
                    truth_writer,
                    "{}\t{}\t{}\t{}\t{}\tA",
                    read_name, umi_str, mol_info.mol_id, ref_name, mol_info.sort_key.pos1
                )?;
                total_pairs += 1;
            }
        }

        progress.log_final();
        truth_writer.flush()?;

        info!("Generated {total_pairs} read pairs");
        info!("Done");

        Ok(())
    }
}

/// Compute the genomic position for a molecule based on its ID.
#[inline]
fn compute_position(mol_id: usize, num_positions: usize, ref_length: usize) -> usize {
    let position_idx = mol_id % num_positions;
    ((position_idx as f64 / num_positions as f64) * (ref_length - 1000) as f64) as usize + 100
}

/// Generate all read pairs for a single molecule.
/// Returns Vec of (`r1_record`, `r2_record`, `read_name`, `umi_str`) tuples.
fn generate_molecule_reads(
    mol_id: usize,
    seed: u64,
    params: &GenerationParams,
) -> Vec<(RecordBuf, RecordBuf, String, String)> {
    let mut rng = create_rng(Some(seed));

    // Note: position is computed externally and passed via mol_id for reproducibility
    // We need to recompute it here for the record generation
    let position = compute_position(mol_id, params.num_positions, params.ref_length);

    // Generate UMI
    let umi = generate_random_sequence(params.umi_length, &mut rng);
    let umi_str = String::from_utf8_lossy(&umi).to_string();

    // Generate family size
    let family_size = params.family_dist.sample(&mut rng, params.min_family_size);

    // Generate insert size
    let insert_size = params.insert_model.sample(&mut rng);

    let mut pairs = Vec::with_capacity(family_size);

    for read_idx in 0..family_size {
        let read_name = format!("mol{mol_id:08}_read{read_idx:04}");

        // Generate template sequence
        let template = generate_random_sequence(insert_size, &mut rng);

        // R1: forward strand at position
        let r1_seq: Vec<u8> = template.iter().take(params.read_length).copied().collect();
        let r1_seq = pad_sequence(r1_seq, params.read_length, &mut rng);

        // R2: reverse strand at position + insert_size - read_length
        let r2_start = insert_size.saturating_sub(params.read_length);
        let r2_template: Vec<u8> =
            template.iter().skip(r2_start).take(params.read_length).copied().collect();
        let r2_seq = reverse_complement(&r2_template);
        let r2_seq = pad_sequence(r2_seq, params.read_length, &mut rng);

        // Quality scores
        let r1_quals = params.quality_model.generate_qualities(params.read_length, &mut rng);
        let r2_quals_raw = params.quality_model.generate_qualities(params.read_length, &mut rng);
        let r2_quals = params.quality_bias.apply_to_vec(&r2_quals_raw, true);

        // Mate cigar is based on the mate's read length (same for both reads)
        let mate_cigar = format!("{}M", params.read_length);

        // Build R1 record
        let r1_record = build_record(
            &read_name,
            &r1_seq,
            &r1_quals,
            0, // ref_id
            position,
            params.mapq,
            true,  // is_first
            false, // is_reverse
            position + insert_size - params.read_length,
            insert_size as i32,
            &umi_str,
            &mate_cigar,
            params.mapq,
        );

        // Build R2 record
        let r2_record = build_record(
            &read_name,
            &r2_seq,
            &r2_quals,
            0,
            position + insert_size - params.read_length,
            params.mapq,
            false, // is_first
            true,  // is_reverse
            position,
            -(insert_size as i32),
            &umi_str,
            &mate_cigar,
            params.mapq,
        );

        pairs.push((r1_record, r2_record, read_name, umi_str.clone()));
    }

    pairs
}

#[allow(clippy::too_many_arguments)]
fn build_record(
    name: &str,
    seq: &[u8],
    quals: &[u8],
    ref_id: usize,
    pos: usize,
    mapq: u8,
    is_first: bool,
    is_reverse: bool,
    mate_pos: usize,
    tlen: i32,
    umi: &str,
    mate_cigar: &str,
    mate_mapq: u8,
) -> RecordBuf {
    let seq_str = String::from_utf8_lossy(seq);

    RecordBuilder::new()
        .name(name)
        .sequence(&seq_str)
        .qualities(quals)
        .reference_sequence_id(ref_id)
        .alignment_start(pos + 1) // Convert 0-based to 1-based
        .mapping_quality(mapq)
        .paired(true)
        .first_segment(is_first)
        .reverse_complement(is_reverse)
        .mate_reverse_complement(!is_reverse)
        .mate_reference_sequence_id(ref_id)
        .mate_alignment_start(mate_pos + 1) // Convert 0-based to 1-based
        .template_length(tlen)
        .tag("RX", umi)
        .tag("MC", mate_cigar)
        .tag("MQ", i32::from(mate_mapq))
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

fn pad_sequence(mut seq: Vec<u8>, target_len: usize, rng: &mut impl Rng) -> Vec<u8> {
    while seq.len() < target_len {
        seq.push(b"ACGT"[rng.random_range(0..4)]);
    }
    seq.truncate(target_len);
    seq
}

#[cfg(test)]
mod tests {
    use super::*;
    use fgumi_lib::simulate::create_rng;

    // Tests for generate_random_sequence
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
    fn test_generate_random_sequence_reproducibility() {
        let mut rng1 = create_rng(Some(42));
        let mut rng2 = create_rng(Some(42));

        let seq1 = generate_random_sequence(100, &mut rng1);
        let seq2 = generate_random_sequence(100, &mut rng2);

        assert_eq!(seq1, seq2);
    }

    // Tests for reverse_complement
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
        assert_eq!(reverse_complement(b"AACG"), b"CGTT");
    }

    #[test]
    fn test_reverse_complement_empty() {
        let empty: &[u8] = b"";
        assert_eq!(reverse_complement(empty), Vec::<u8>::new());
    }

    #[test]
    fn test_reverse_complement_unknown_base() {
        assert_eq!(reverse_complement(b"N"), b"N");
        assert_eq!(reverse_complement(b"X"), b"N");
        assert_eq!(reverse_complement(b"ANCG"), b"CGNT");
    }

    #[test]
    fn test_reverse_complement_double_rc() {
        let seq = b"ACGTACGTACGT";
        let rc = reverse_complement(seq);
        let rc_rc = reverse_complement(&rc);
        assert_eq!(rc_rc, seq.to_vec());
    }

    // Tests for pad_sequence
    #[test]
    fn test_pad_sequence_already_correct_length() {
        let mut rng = create_rng(Some(42));
        let seq = vec![b'A', b'C', b'G', b'T'];
        let padded = pad_sequence(seq.clone(), 4, &mut rng);
        assert_eq!(padded, seq);
    }

    #[test]
    fn test_pad_sequence_needs_padding() {
        let mut rng = create_rng(Some(42));
        let seq = vec![b'A', b'C'];
        let padded = pad_sequence(seq, 6, &mut rng);
        assert_eq!(padded.len(), 6);
        assert_eq!(&padded[0..2], b"AC");
        // Rest should be valid bases
        for &base in &padded[2..] {
            assert!(base == b'A' || base == b'C' || base == b'G' || base == b'T');
        }
    }

    #[test]
    fn test_pad_sequence_truncates() {
        let mut rng = create_rng(Some(42));
        let seq = vec![b'A', b'C', b'G', b'T', b'A', b'C'];
        let padded = pad_sequence(seq, 3, &mut rng);
        assert_eq!(padded, vec![b'A', b'C', b'G']);
    }

    #[test]
    fn test_pad_sequence_empty_to_length() {
        let mut rng = create_rng(Some(42));
        let seq: Vec<u8> = vec![];
        let padded = pad_sequence(seq, 5, &mut rng);
        assert_eq!(padded.len(), 5);
        for &base in &padded {
            assert!(base == b'A' || base == b'C' || base == b'G' || base == b'T');
        }
    }

    #[test]
    fn test_pad_sequence_to_zero() {
        let mut rng = create_rng(Some(42));
        let seq = vec![b'A', b'C', b'G', b'T'];
        let padded = pad_sequence(seq, 0, &mut rng);
        assert!(padded.is_empty());
    }

    // Tests for build_record
    #[test]
    fn test_build_record_r1_basic() {
        let seq = b"ACGT";
        let quals = vec![30, 30, 30, 30];
        let record = build_record(
            "test_read",
            seq,
            &quals,
            0,
            100,
            60,
            true,
            false,
            200,
            150,
            "AAAAAAAA",
            "4M",
            60,
        );

        // Check name
        assert!(record.name().is_some());

        // Check flags for R1
        let flags = record.flags();
        assert!(flags.is_segmented());
        assert!(flags.is_first_segment());
        assert!(!flags.is_last_segment());
        assert!(!flags.is_reverse_complemented());
        assert!(flags.is_mate_reverse_complemented());
    }

    #[test]
    fn test_build_record_r2_basic() {
        let seq = b"ACGT";
        let quals = vec![30, 30, 30, 30];
        let record = build_record(
            "test_read",
            seq,
            &quals,
            0,
            200,
            60,
            false,
            true,
            100,
            -150,
            "AAAAAAAA",
            "4M",
            60,
        );

        // Check flags for R2
        let flags = record.flags();
        assert!(flags.is_segmented());
        assert!(!flags.is_first_segment());
        assert!(flags.is_last_segment());
        assert!(flags.is_reverse_complemented());
        assert!(!flags.is_mate_reverse_complemented());
    }

    #[test]
    fn test_build_record_positions() {
        let seq = b"ACGTACGT";
        let quals = vec![30; 8];
        let record =
            build_record("read", seq, &quals, 0, 1000, 60, true, false, 1100, 200, "AAA", "8M", 60);

        // Position should be 1-based in BAM
        assert_eq!(record.alignment_start().unwrap().get(), 1001);
        assert_eq!(record.mate_alignment_start().unwrap().get(), 1101);
        assert_eq!(record.template_length(), 200);
    }

    #[test]
    fn test_build_record_mapping_quality() {
        let seq = b"ACGT";
        let quals = vec![30; 4];

        // Note: mapq 255 is treated as "unavailable" by noodles, so skip it
        for mapq in [0, 30, 60] {
            let record = build_record(
                "read", seq, &quals, 0, 100, mapq, true, false, 200, 100, "AAA", "4M", mapq,
            );
            assert_eq!(record.mapping_quality().unwrap().get(), mapq);
        }
    }

    #[test]
    fn test_build_record_quality_scores() {
        let seq = b"ACGTACGT";
        let quals = vec![10, 20, 30, 40, 30, 20, 10, 5];
        let record =
            build_record("read", seq, &quals, 0, 100, 60, true, false, 200, 100, "AAA", "8M", 60);

        let record_quals: Vec<u8> = record.quality_scores().iter().collect();
        assert_eq!(record_quals, quals);
    }

    #[test]
    fn test_build_record_reference_id() {
        let seq = b"ACGT";
        let quals = vec![30; 4];

        for ref_id in [0, 1, 5] {
            let record = build_record(
                "read", seq, &quals, ref_id, 100, 60, true, false, 200, 100, "AAA", "4M", 60,
            );
            assert_eq!(record.reference_sequence_id(), Some(ref_id));
            assert_eq!(record.mate_reference_sequence_id(), Some(ref_id));
        }
    }

    #[test]
    fn test_build_record_negative_tlen() {
        let seq = b"ACGT";
        let quals = vec![30; 4];
        let record =
            build_record("read", seq, &quals, 0, 200, 60, false, true, 100, -150, "AAA", "4M", 60);

        assert_eq!(record.template_length(), -150);
    }
}

//! Generate grouped BAM with MI tags for consensus calling.
//!
//! Uses a memory-efficient approach: sort molecule IDs by their position first,
//! then generate records in sorted order (streaming). This avoids storing all
//! records in memory before sorting.

use super::sort::TemplateCoordKey;
use crate::commands::command::Command;
use crate::commands::common::CompressionOptions;
use crate::commands::simulate::common::{
    FamilySizeArgs, InsertSizeArgs, PositionDistArgs, QualityArgs, ReferenceArgs, SimulationCommon,
    StrandBiasArgs,
};
use anyhow::{Context, Result};
use bstr::BString;
use clap::Parser;
use fgumi_lib::bam_io::create_bam_writer;
use fgumi_lib::progress::ProgressTracker;
use fgumi_lib::sam::builder::RecordBuilder;
use fgumi_lib::simulate::{
    FamilySizeDistribution, InsertSizeModel, PositionQualityModel, ReadPairQualityBias,
    StrandBiasModel, create_rng,
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

/// Generate template-coordinate sorted BAM with MI tags for consensus callers.
#[derive(Parser, Debug)]
#[command(
    name = "grouped-reads",
    about = "Generate grouped BAM with MI tags for consensus calling",
    long_about = r#"
Generate synthetic BAM files with MI (molecule ID) tags already assigned.

The output is template-coordinate sorted and suitable for input to
`fgumi simplex`, `fgumi duplex`, or `fgumi codec`.
"#
)]
pub struct GroupedReads {
    /// Output BAM file (template-coordinate sorted)
    #[arg(short = 'o', long = "output", required = true)]
    pub output: PathBuf,

    /// Output truth TSV file for validation
    #[arg(long = "truth", required = true)]
    pub truth_output: PathBuf,

    /// Generate duplex-style MI tags (e.g., "1/A", "1/B")
    #[arg(long = "duplex")]
    pub duplex: bool,

    /// Mapping quality for aligned reads
    #[arg(long = "mapq", default_value = "60")]
    pub mapq: u8,

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

    #[command(flatten)]
    pub strand_bias: StrandBiasArgs,
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
    duplex: bool,
    min_family_size: usize,
    num_positions: usize,
    ref_length: usize,
    quality_model: PositionQualityModel,
    quality_bias: ReadPairQualityBias,
    family_dist: FamilySizeDistribution,
    insert_model: InsertSizeModel,
    strand_bias_model: StrandBiasModel,
}

impl Command for GroupedReads {
    fn execute(&self, command_line: &str) -> Result<()> {
        info!("Generating grouped reads");
        info!("  Output: {}", self.output.display());
        info!("  Truth: {}", self.truth_output.display());
        info!("  Duplex: {}", self.duplex);
        info!("  Num molecules: {}", self.common.num_molecules);
        info!("  Read length: {}", self.common.read_length);
        info!("  UMI length: {}", self.common.umi_length);
        info!("  Threads: {}", self.threads);

        // Determine positions
        let num_positions = self.position_dist.num_positions.unwrap_or(self.common.num_molecules);
        let ref_length = self.reference.ref_length;

        // Check for position collisions that would cause UMI conflicts
        let usable_bases = ref_length.saturating_sub(1000);
        let bases_per_position = usable_bases as f64 / num_positions as f64;
        if bases_per_position < 1.0 {
            let suggested_ref_length = num_positions * 2;
            anyhow::bail!(
                "Position collision: {num_positions} positions cannot fit in {ref_length} bp reference ({bases_per_position:.2} bp/position). \
                 Increase --ref-length to at least {suggested_ref_length} or reduce --num-molecules."
            );
        } else if bases_per_position < 10.0 {
            log::warn!(
                "Low position spacing ({bases_per_position:.1} bp/position) may cause UMI collisions. \
                 Consider increasing --ref-length."
            );
        }

        // Build header with template-coordinate sort order
        let ref_name = self.reference.ref_name.clone();
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
        header = crate::commands::common::add_pg_to_builder(header, command_line)?;

        let header = header.build();

        // Set up generation parameters
        let params = GenerationParams {
            read_length: self.common.read_length,
            umi_length: self.common.umi_length,
            mapq: self.mapq,
            duplex: self.duplex,
            min_family_size: self.family_size.min_family_size,
            num_positions,
            ref_length,
            quality_model: self.quality.to_quality_model(),
            quality_bias: self.quality.to_quality_bias(),
            family_dist: self.family_size.to_family_size_distribution()?,
            insert_model: self.insert_size.to_insert_size_model(),
            strand_bias_model: self.strand_bias.to_strand_bias_model(),
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
                // MI tag is mol_id.to_string() (or "{mol_id}/A" for duplex, but /A /B suffix is stripped)
                let sort_key = TemplateCoordKey::for_f1r2_pair(
                    0, // tid - all simulated reads on reference 0
                    pos1,
                    insert_size,
                    mol_id.to_string(), // MI tag (stripped suffix matches this)
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
        writeln!(
            truth_writer,
            "read_name\ttrue_umi\tmolecule_id\tmi_tag\tchrom\tposition\tstrand"
        )?;

        // Generate and write records in position-sorted order (streaming)
        info!("Generating and writing records in sorted order...");
        let progress = ProgressTracker::new("Processed molecules").with_interval(100_000);
        let mut total_pairs = 0usize;

        for mol_info in molecules {
            progress.log_if_needed(1);

            // Generate all read pairs for this molecule
            let pairs = generate_molecule_reads(mol_info.mol_id, mol_info.seed, &params);

            for (r1, r2, read_name, umi_str, mi_tag, strand, insert_size) in pairs {
                // Write reads in template-coordinate order (earlier 5' position first).
                // For A strand (F1R2): R1 forward at pos, R2 reverse at pos+insert-read_len
                //   R1's 5' = pos, R2's 5' = pos+insert-1 → R1 first (always)
                // For B strand (R1F2): R1 reverse at pos, R2 forward at pos+insert-read_len
                //   R1's 5' = pos+read_len-1, R2's 5' = pos+insert-read_len
                //   R2 first if insert <= 2*read_len-1, else R1 first
                //   (When equal, samtools puts forward strand first)
                let r2_first = strand == 'B' && insert_size < 2 * params.read_length;
                if r2_first {
                    writer.write_alignment_record(&header, &r2)?;
                    writer.write_alignment_record(&header, &r1)?;
                } else {
                    writer.write_alignment_record(&header, &r1)?;
                    writer.write_alignment_record(&header, &r2)?;
                }
                writeln!(
                    truth_writer,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    read_name,
                    umi_str,
                    mol_info.mol_id,
                    mi_tag,
                    ref_name,
                    mol_info.sort_key.pos1,
                    strand
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
/// Returns Vec of (`r1_record`, `r2_record`, `read_name`, `umi_str`, `mi_tag`, strand, `insert_size`) tuples.
fn generate_molecule_reads(
    mol_id: usize,
    seed: u64,
    params: &GenerationParams,
) -> Vec<(RecordBuf, RecordBuf, String, String, String, char, usize)> {
    let mut rng = create_rng(Some(seed));

    // Compute position for this molecule
    let position = compute_position(mol_id, params.num_positions, params.ref_length);

    // Generate UMI
    let umi = generate_random_sequence(params.umi_length, &mut rng);
    let umi_str = String::from_utf8_lossy(&umi).to_string();

    // Generate family size
    let family_size = params.family_dist.sample(&mut rng, params.min_family_size);

    // Generate insert size
    let insert_size = params.insert_model.sample(&mut rng);

    let mut pairs = Vec::new();

    if params.duplex {
        // Split reads between A and B strands
        let (a_count, b_count) = params.strand_bias_model.split_reads(family_size, &mut rng);

        // Generate A strand reads
        for read_idx in 0..a_count {
            let read_name = format!("mol{mol_id:08}_readA{read_idx:04}");
            let mi_tag = format!("{mol_id}/A");

            let (r1_record, r2_record) = generate_read_pair_records(
                &read_name,
                &umi_str,
                &mi_tag,
                position,
                insert_size,
                params.read_length,
                params.mapq,
                'A',
                &params.quality_model,
                &params.quality_bias,
                &mut rng,
            );

            pairs.push((
                r1_record,
                r2_record,
                read_name,
                umi_str.clone(),
                mi_tag,
                'A',
                insert_size,
            ));
        }

        // Generate B strand reads
        for read_idx in 0..b_count {
            let read_name = format!("mol{mol_id:08}_readB{read_idx:04}");
            let mi_tag = format!("{mol_id}/B");

            let (r1_record, r2_record) = generate_read_pair_records(
                &read_name,
                &umi_str,
                &mi_tag,
                position,
                insert_size,
                params.read_length,
                params.mapq,
                'B',
                &params.quality_model,
                &params.quality_bias,
                &mut rng,
            );

            pairs.push((
                r1_record,
                r2_record,
                read_name,
                umi_str.clone(),
                mi_tag,
                'B',
                insert_size,
            ));
        }
    } else {
        // Simplex mode - all reads get same MI tag
        let mi_tag = mol_id.to_string();

        for read_idx in 0..family_size {
            let read_name = format!("mol{mol_id:08}_read{read_idx:04}");

            let (r1_record, r2_record) = generate_read_pair_records(
                &read_name,
                &umi_str,
                &mi_tag,
                position,
                insert_size,
                params.read_length,
                params.mapq,
                '+',
                &params.quality_model,
                &params.quality_bias,
                &mut rng,
            );

            pairs.push((
                r1_record,
                r2_record,
                read_name,
                umi_str.clone(),
                mi_tag.clone(),
                '+',
                insert_size,
            ));
        }
    }

    pairs
}

#[allow(clippy::too_many_arguments)]
fn generate_read_pair_records(
    read_name: &str,
    umi_str: &str,
    mi_tag: &str,
    position: usize,
    insert_size: usize,
    read_length: usize,
    mapq: u8,
    strand: char,
    quality_model: &PositionQualityModel,
    quality_bias: &ReadPairQualityBias,
    rng: &mut impl Rng,
) -> (RecordBuf, RecordBuf) {
    // Generate template sequence
    let template = generate_random_sequence(insert_size, rng);

    // R1: forward strand at position
    let r1_seq: Vec<u8> = template.iter().take(read_length).copied().collect();
    let r1_seq = pad_sequence(r1_seq, read_length, rng);

    // R2: reverse strand at position + insert_size - read_length
    let r2_start = insert_size.saturating_sub(read_length);
    let r2_template: Vec<u8> = template.iter().skip(r2_start).take(read_length).copied().collect();
    let r2_seq = reverse_complement(&r2_template);
    let r2_seq = pad_sequence(r2_seq, read_length, rng);

    // Quality scores
    let r1_quals = quality_model.generate_qualities(read_length, rng);
    let r2_quals_raw = quality_model.generate_qualities(read_length, rng);
    let r2_quals = quality_bias.apply_to_vec(&r2_quals_raw, true);

    // Mate cigar is based on the mate's read length (same for both reads)
    let mate_cigar = format!("{read_length}M");

    // Strand orientation: A strand is normal (R1 forward, R2 reverse),
    // B strand is opposite (R1 reverse, R2 forward) to represent the
    // complementary DNA strand in duplex sequencing.
    let (r1_is_reverse, r2_is_reverse) = match strand {
        'B' => (true, false), // B strand: R1 reverse, R2 forward
        _ => (false, true),   // A strand/simplex: R1 forward, R2 reverse
    };

    // Build R1 record
    let r1_record = build_record(
        read_name,
        &r1_seq,
        &r1_quals,
        0, // ref_id
        position,
        mapq,
        true, // is_first
        r1_is_reverse,
        position + insert_size - read_length,
        insert_size as i32,
        umi_str,
        mi_tag,
        &mate_cigar,
        mapq,
    );

    // Build R2 record
    let r2_record = build_record(
        read_name,
        &r2_seq,
        &r2_quals,
        0,
        position + insert_size - read_length,
        mapq,
        false, // is_first
        r2_is_reverse,
        position,
        -(insert_size as i32),
        umi_str,
        mi_tag,
        &mate_cigar,
        mapq,
    );

    (r1_record, r2_record)
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
    mi_tag: &str,
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
        .tag("MI", mi_tag)
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
        assert_eq!(reverse_complement(b"AACG"), b"CGTT");
    }

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
    }

    #[test]
    fn test_build_record_with_mi_tag() {
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
            "42/A",
            "4M",
            60,
        );

        // Check name exists
        assert!(record.name().is_some());

        // Check flags for R1
        let flags = record.flags();
        assert!(flags.is_segmented());
        assert!(flags.is_first_segment());
    }

    #[test]
    fn test_build_record_simplex_mi_tag() {
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
            "42",
            "4M",
            60,
        );

        assert!(record.name().is_some());
    }

    #[test]
    fn test_generate_read_pair_records() {
        let mut rng = create_rng(Some(42));
        let quality_model = PositionQualityModel::default();
        let quality_bias = ReadPairQualityBias::default();

        let (r1, r2) = generate_read_pair_records(
            "read_001",
            "ACGTACGT",
            "0/A",
            1000,
            300,
            150,
            60,
            'A',
            &quality_model,
            &quality_bias,
            &mut rng,
        );

        // Check R1 flags
        assert!(r1.flags().is_first_segment());
        assert!(!r1.flags().is_reverse_complemented());

        // Check R2 flags
        assert!(r2.flags().is_last_segment());
        assert!(r2.flags().is_reverse_complemented());
    }

    #[test]
    fn test_generate_read_pair_records_b_strand() {
        let mut rng = create_rng(Some(42));
        let quality_model = PositionQualityModel::default();
        let quality_bias = ReadPairQualityBias::default();

        let (r1, r2) = generate_read_pair_records(
            "read_001",
            "ACGTACGT",
            "5/B",
            1000,
            300,
            150,
            60,
            'B',
            &quality_model,
            &quality_bias,
            &mut rng,
        );

        // B strand has opposite orientation: R1 reverse, R2 forward
        assert!(r1.flags().is_first_segment());
        assert!(r1.flags().is_reverse_complemented());
        assert!(r2.flags().is_last_segment());
        assert!(!r2.flags().is_reverse_complemented());
    }

    #[test]
    fn test_generate_read_pair_records_simplex() {
        let mut rng = create_rng(Some(42));
        let quality_model = PositionQualityModel::default();
        let quality_bias = ReadPairQualityBias::default();

        let (r1, r2) = generate_read_pair_records(
            "read_001",
            "ACGTACGT",
            "10",
            1000,
            300,
            150,
            60,
            '+',
            &quality_model,
            &quality_bias,
            &mut rng,
        );

        // Simplex has same orientation as A strand
        assert!(r1.flags().is_first_segment());
        assert!(!r1.flags().is_reverse_complemented());
        assert!(r2.flags().is_last_segment());
        assert!(r2.flags().is_reverse_complemented());
    }

    #[test]
    fn test_generate_read_pair_records_small_insert() {
        let mut rng = create_rng(Some(42));
        let quality_model = PositionQualityModel::default();
        let quality_bias = ReadPairQualityBias::default();

        // Insert size smaller than read length
        let (r1, r2) = generate_read_pair_records(
            "read_001",
            "ACGTACGT",
            "0/A",
            1000,
            50, // insert_size < read_length
            150,
            60,
            'A',
            &quality_model,
            &quality_bias,
            &mut rng,
        );

        // Both records should still be valid
        assert!(r1.flags().is_first_segment());
        assert!(r2.flags().is_last_segment());
    }

    #[test]
    fn test_generate_read_pair_records_reproducibility() {
        let quality_model = PositionQualityModel::default();
        let quality_bias = ReadPairQualityBias::default();

        let mut rng1 = create_rng(Some(42));
        let mut rng2 = create_rng(Some(42));

        let (r1_a, r2_a) = generate_read_pair_records(
            "read_001",
            "ACGTACGT",
            "0/A",
            1000,
            300,
            150,
            60,
            'A',
            &quality_model,
            &quality_bias,
            &mut rng1,
        );

        let (r1_b, r2_b) = generate_read_pair_records(
            "read_001",
            "ACGTACGT",
            "0/A",
            1000,
            300,
            150,
            60,
            'A',
            &quality_model,
            &quality_bias,
            &mut rng2,
        );

        // Both should produce same flags
        assert_eq!(r1_a.flags(), r1_b.flags());
        assert_eq!(r2_a.flags(), r2_b.flags());
    }

    #[test]
    fn test_build_record_r2_flags() {
        let seq = b"ACGT";
        let quals = vec![30, 30, 30, 30];
        let record = build_record(
            "test_read",
            seq,
            &quals,
            0,
            200,
            60,
            false, // R2
            true,  // reverse
            100,
            -150,
            "AAAAAAAA",
            "42/B",
            "4M",
            60,
        );

        let flags = record.flags();
        assert!(flags.is_last_segment());
        assert!(!flags.is_first_segment());
        assert!(flags.is_reverse_complemented());
        assert!(!flags.is_mate_reverse_complemented());
    }

    #[test]
    fn test_build_record_positions() {
        let seq = b"ACGTACGT";
        let quals = vec![30; 8];
        let record = build_record(
            "read", seq, &quals, 0, 1000, 60, true, false, 1100, 200, "AAA", "0/A", "8M", 60,
        );

        // Position should be 1-based in BAM
        assert_eq!(record.alignment_start().unwrap().get(), 1001);
        assert_eq!(record.mate_alignment_start().unwrap().get(), 1101);
        assert_eq!(record.template_length(), 200);
    }

    #[test]
    fn test_build_record_mapping_quality() {
        let seq = b"ACGT";
        let quals = vec![30; 4];

        for mapq in [0, 30, 60] {
            let record = build_record(
                "read", seq, &quals, 0, 100, mapq, true, false, 200, 100, "AAA", "0", "4M", mapq,
            );
            assert_eq!(record.mapping_quality().unwrap().get(), mapq);
        }
    }

    #[test]
    fn test_build_record_quality_scores() {
        let seq = b"ACGTACGT";
        let quals = vec![10, 20, 30, 40, 30, 20, 10, 5];
        let record = build_record(
            "read", seq, &quals, 0, 100, 60, true, false, 200, 100, "AAA", "0", "8M", 60,
        );

        let record_quals: Vec<u8> = record.quality_scores().iter().collect();
        assert_eq!(record_quals, quals);
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
    fn test_reverse_complement_empty() {
        let empty: &[u8] = b"";
        assert_eq!(reverse_complement(empty), Vec::<u8>::new());
    }

    #[test]
    fn test_reverse_complement_unknown_base() {
        assert_eq!(reverse_complement(b"N"), b"N");
        assert_eq!(reverse_complement(b"ANCG"), b"CGNT");
    }

    /// Helper to test collision detection logic
    fn check_collision(num_positions: usize, ref_length: usize) -> Result<(), String> {
        let usable_bases = ref_length.saturating_sub(1000);
        let bases_per_position = usable_bases as f64 / num_positions as f64;
        if bases_per_position < 1.0 {
            Err(format!(
                "Position collision: {num_positions} positions cannot fit in {ref_length} bp reference ({bases_per_position:.2} bp/position)"
            ))
        } else {
            Ok(())
        }
    }

    #[test]
    fn test_collision_detection_error() {
        // 5M molecules in 10MB reference = ~2 bp/position, but modulo causes collisions
        // With 5M positions and 10M - 1000 usable bases, we get ~1.8 bp/position
        let result = check_collision(5_000_000, 10_000_000);
        assert!(result.is_ok(), "Should pass with ~1.8 bp/position");

        // 10M molecules in 10MB reference = <1 bp/position, should fail
        let result = check_collision(10_000_000, 10_000_000);
        assert!(result.is_err(), "Should fail with <1 bp/position");

        // 20M molecules in 10MB reference = definitely fails
        let result = check_collision(20_000_000, 10_000_000);
        assert!(result.is_err());
    }

    #[test]
    fn test_collision_detection_success() {
        // 5M molecules in 250MB reference = ~50 bp/position, plenty of room
        let result = check_collision(5_000_000, 250_000_000);
        assert!(result.is_ok());

        // 1M molecules in 250MB reference = ~250 bp/position
        let result = check_collision(1_000_000, 250_000_000);
        assert!(result.is_ok());
    }

    #[test]
    fn test_a_strand_orientation() {
        let mut rng = create_rng(Some(42));
        let quality_model = PositionQualityModel::default();
        let quality_bias = ReadPairQualityBias::default();

        let (r1, r2) = generate_read_pair_records(
            "read_001",
            "ACGTACGT",
            "1/A",
            1000,
            300,
            150,
            60,
            'A',
            &quality_model,
            &quality_bias,
            &mut rng,
        );

        // A strand: R1 forward, R2 reverse
        assert!(r1.flags().is_first_segment());
        assert!(!r1.flags().is_reverse_complemented());
        assert!(r2.flags().is_last_segment());
        assert!(r2.flags().is_reverse_complemented());
    }
}

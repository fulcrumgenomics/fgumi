//! Generate grouped BAM with MI tags for consensus calling.
//!
//! Uses a memory-efficient approach: sort molecule IDs by their position first,
//! then generate records in sorted order (streaming). This avoids storing all
//! records in memory before sorting.

use super::sort::TemplateCoordKey;
use crate::bam_io::create_bam_writer;
use crate::commands::command::Command;
use crate::commands::common::{CompressionOptions, parse_bool};
use crate::commands::simulate::common::{
    FamilySizeArgs, InsertSizeArgs, MethylationArgs, MethylationConfig, MoleculeInfo,
    PositionDistArgs, QualityArgs, ReferenceArgs, ReferenceGenome, SimulationCommon,
    StrandBiasArgs, apply_methylation_conversion, generate_random_sequence, pad_sequence,
};
use crate::dna::reverse_complement;
use crate::progress::ProgressTracker;
use crate::sam::builder::RecordBuilder;
use crate::simulate::{
    FamilySizeDistribution, InsertSizeModel, PositionQualityModel, ReadPairQualityBias,
    StrandBiasModel, create_rng,
};
use anyhow::{Context, Result};
use clap::Parser;
use log::info;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::sam::header::Header;
use noodles::sam::header::record::value::map::header::{self as HeaderRecord, Tag as HeaderTag};
use rand::{Rng, RngExt};
use std::fs::File;
use std::io::{BufWriter, Write};
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
    #[arg(long = "duplex", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool)]
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

    #[command(flatten)]
    pub methylation: MethylationArgs,
}

/// Parameters needed for molecule generation.
struct GenerationParams {
    read_length: usize,
    umi_length: usize,
    mapq: u8,
    duplex: bool,
    min_family_size: usize,
    quality_model: PositionQualityModel,
    quality_bias: ReadPairQualityBias,
    family_dist: FamilySizeDistribution,
    insert_model: InsertSizeModel,
    strand_bias_model: StrandBiasModel,
    methylation: MethylationConfig,
}

impl Command for GroupedReads {
    fn execute(&self, command_line: &str) -> Result<()> {
        // Validate methylation args
        let methylation = self.methylation.resolve();
        self.methylation.validate()?;

        info!("Generating grouped reads");
        info!("  Output: {}", self.output.display());
        info!("  Truth: {}", self.truth_output.display());
        info!("  Duplex: {}", self.duplex);
        info!("  Num molecules: {}", self.common.num_molecules);
        info!("  Read length: {}", self.common.read_length);
        info!("  UMI length: {}", self.common.umi_length);
        info!("  Threads: {}", self.threads);
        if methylation.mode.is_enabled() {
            info!("  Methylation mode: {:?}", methylation.mode);
            info!("  CpG methylation rate: {}", methylation.cpg_methylation_rate);
            info!("  Conversion rate: {}", methylation.conversion_rate);
        }

        // Load reference genome
        let ref_genome = ReferenceGenome::load(&self.reference.reference)?;

        // Validate that the reference has at least one contig >= read_length
        if ref_genome.max_contig_length() < self.common.read_length {
            anyhow::bail!(
                "No reference contig is >= read length ({} bp). \
                 The longest contig is {} bp. Use a larger reference or shorter --read-length.",
                self.common.read_length,
                ref_genome.max_contig_length(),
            );
        }

        // Build header from reference contigs + sort order tags + @PG
        let ref_header = ref_genome.build_bam_header();
        let mut header_builder = Header::builder();

        // Add sort order tags: SO:unsorted, GO:query, SS:template-coordinate
        let HeaderTag::Other(so_tag) = HeaderTag::from([b'S', b'O']) else { unreachable!() };
        let HeaderTag::Other(go_tag) = HeaderTag::from([b'G', b'O']) else { unreachable!() };
        let HeaderTag::Other(ss_tag) = HeaderTag::from([b'S', b'S']) else { unreachable!() };

        let header_map =
            noodles::sam::header::record::value::Map::<HeaderRecord::Header>::builder()
                .insert(so_tag, "unsorted")
                .insert(go_tag, "query")
                .insert(ss_tag, "template-coordinate")
                .build()
                .expect("header map with valid SO/GO/SS tags");
        header_builder = header_builder.set_header(header_map);

        for (name, map) in ref_header.reference_sequences() {
            header_builder = header_builder.add_reference_sequence(name.clone(), map.clone());
        }

        header_builder = crate::commands::common::add_pg_to_builder(header_builder, command_line)?;

        let header = header_builder.build();

        // Set up generation parameters
        let params = GenerationParams {
            read_length: self.common.read_length,
            umi_length: self.common.umi_length,
            mapq: self.mapq,
            duplex: self.duplex,
            min_family_size: self.family_size.min_family_size,
            quality_model: self.quality.to_quality_model(),
            quality_bias: self.quality.to_quality_bias(),
            family_dist: self.family_size.to_family_size_distribution()?,
            insert_model: self.insert_size.to_insert_size_model(),
            strand_bias_model: self.strand_bias.to_strand_bias_model(),
            methylation,
        };

        // Pre-sample positions from reference (use a dedicated RNG so molecule seeds
        // stay deterministic and uncorrelated with position sampling)
        let num_positions = self.position_dist.num_positions.unwrap_or(self.common.num_molecules);
        if num_positions == 0 {
            anyhow::bail!("--num-positions must be greater than 0");
        }

        // Collision check using total_length
        let usable_bases = ref_genome.total_length().saturating_sub(1000);
        let bases_per_position = usable_bases as f64 / num_positions as f64;
        if bases_per_position < 1.0 {
            anyhow::bail!(
                "Too many positions ({num_positions}) for reference of size {} bp. \
                 Reduce --num-positions or use a larger reference.",
                ref_genome.total_length()
            );
        } else if bases_per_position < 10.0 {
            log::warn!(
                "Low position spacing ({bases_per_position:.1} bp/position) may cause UMI collisions."
            );
        }

        let mut pos_rng = create_rng(self.common.seed);
        let position_table = ref_genome.sample_positions(num_positions, &mut pos_rng);

        // Use a different seed for molecule generation to avoid correlation with positions
        let mut seed_rng = create_rng(self.common.seed.map(|s| s.wrapping_add(1)));

        // MEMORY-EFFICIENT APPROACH: Sort molecule IDs by template-coordinate key first,
        // then generate records in sorted order (streaming).
        //
        // We need to pre-compute insert_size to get pos2 for the sort key, since
        // template-coordinate sorting considers both ends of the template.
        info!("Computing molecule positions and sort keys...");
        let mut molecules: Vec<MoleculeInfo> = (0..self.common.num_molecules)
            .map(|mol_id| {
                let seed: u64 = seed_rng.random();
                let pos_idx = mol_id % num_positions;
                let (chrom_idx, local_pos) = position_table[pos_idx];

                // Pre-compute insert_size using the molecule's seed (same RNG sequence
                // as generation)
                let mut mol_rng = create_rng(Some(seed));
                // Skip UMI generation RNG calls
                for _ in 0..params.umi_length {
                    let _: usize = mol_rng.random_range(0..4);
                }
                // Skip family_size RNG call
                let _ = params.family_dist.sample(&mut mol_rng, params.min_family_size);
                // Get insert_size
                let insert_size = params.insert_model.sample(&mut mol_rng);

                // Sort key — for_f1r2_pair gives the canonical template-coordinate sort
                // key. This works for both F1R2 and R1F2 orientations because the
                // template covers the same genomic positions regardless of strand.
                let sort_key = TemplateCoordKey::for_f1r2_pair(
                    chrom_idx as i32, // real tid
                    local_pos,
                    insert_size,
                    mol_id.to_string(), // MI tag (stripped suffix matches this)
                    format!("mol{mol_id:08}"),
                );
                MoleculeInfo { mol_id, seed, sort_key, is_unmapped: false }
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

            let pos_idx = mol_info.mol_id % num_positions;
            let (pos_chrom_idx, pos_local_pos) = position_table[pos_idx];

            // Generate all read pairs for this molecule. If the pre-sampled locus
            // could not provide a valid template, `generate_molecule_reads` falls
            // back to sampling a new locus and returns the effective coordinates
            // so BAM records and truth rows stay in sync with the emitted sequence.
            let (pairs, chrom_idx, local_pos) = generate_molecule_reads(
                mol_info.mol_id,
                mol_info.seed,
                pos_chrom_idx,
                pos_local_pos,
                &params,
                &ref_genome,
            );

            for (r1, r2, read_name, umi_str, mi_tag, is_top_strand, insert_size) in pairs {
                // Write reads in template-coordinate order (earlier 5' position first).
                // For F1R2 (top strand): R1 forward at pos, R2 reverse at pos+insert-read_len
                //   R1's 5' = pos, R2's 5' = pos+insert-1 → R1 first (always)
                // For R1F2 (!top strand): R1 reverse at end, R2 forward at start
                //   R2 first if insert < 2*read_len, else R1 first
                let r2_first = !is_top_strand && insert_size < 2 * params.read_length;
                if r2_first {
                    writer.write_alignment_record(&header, &r2)?;
                    writer.write_alignment_record(&header, &r1)?;
                } else {
                    writer.write_alignment_record(&header, &r1)?;
                    writer.write_alignment_record(&header, &r2)?;
                }
                let strand_char = if is_top_strand { '+' } else { '-' };
                writeln!(
                    truth_writer,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    read_name,
                    umi_str,
                    mol_info.mol_id,
                    mi_tag,
                    ref_genome.name(chrom_idx),
                    local_pos,
                    strand_char,
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

/// Generate all read pairs for a single molecule.
///
/// A single R1/R2 pair plus identifiers used for truth/metadata:
/// `(r1, r2, read_name, umi_str, mi, is_top_strand, family_size)`.
type MoleculeReadPair = (RecordBuf, RecordBuf, String, String, String, bool, usize);

/// Result of generating all read pairs for a molecule: the pairs plus the
/// effective `(chrom_idx, local_pos)` at which the template actually lives.
type MoleculeReadsResult = (Vec<MoleculeReadPair>, usize, usize);

/// Returns the pair vector together with the effective `(chrom_idx, local_pos)` —
/// the fallback may re-sample a different locus when the pre-sampled one doesn't
/// yield a valid template, and callers must propagate those coordinates to BAM
/// records and the truth TSV so reported positions match the emitted sequence.
#[allow(clippy::too_many_arguments)]
fn generate_molecule_reads(
    mol_id: usize,
    seed: u64,
    chrom_idx: usize,
    local_pos: usize,
    params: &GenerationParams,
    ref_genome: &ReferenceGenome,
) -> MoleculeReadsResult {
    let mut rng = create_rng(Some(seed));

    // Generate UMI
    let umi = generate_random_sequence(params.umi_length, &mut rng);
    let umi_str = String::from_utf8_lossy(&umi).to_string();

    // Generate family size
    let family_size = params.family_dist.sample(&mut rng, params.min_family_size);

    // Generate insert size
    let insert_size = params.insert_model.sample(&mut rng);

    // 50/50 strand coin flip: determines genomic strand of origin
    let is_top_strand: bool = rng.random();

    // Get template from reference at the pre-sampled position, falling back to a
    // random position if the exact location doesn't yield a valid sequence. When
    // falling back, adopt the new locus so records and truth stay aligned with the
    // sequence actually emitted.
    let (eff_chrom_idx, eff_local_pos, shared_template) =
        match ref_genome.sequence_at(chrom_idx, local_pos, insert_size) {
            Some(seq) => (chrom_idx, local_pos, Some(seq)),
            None => match ref_genome.sample_sequence(insert_size, &mut rng) {
                Some((c, p, seq)) => (c, p, Some(seq)),
                None => (chrom_idx, local_pos, None),
            },
        };
    let chrom_idx = eff_chrom_idx;
    let local_pos = eff_local_pos;

    let mut pairs = Vec::new();

    if params.duplex {
        // Split reads between A and B strands
        let (a_count, b_count) = params.strand_bias_model.split_reads(family_size, &mut rng);

        // A reads: orientation follows the coin flip
        let a_is_top = is_top_strand;
        for read_idx in 0..a_count {
            let read_name = format!("mol{mol_id:08}_readA{read_idx:04}");
            let mi_tag = format!("{mol_id}/A");

            let (r1_record, r2_record) = generate_read_pair_records(
                &read_name,
                &umi_str,
                &mi_tag,
                chrom_idx,
                local_pos,
                insert_size,
                params.read_length,
                params.mapq,
                a_is_top,
                &params.quality_model,
                &params.quality_bias,
                &params.methylation,
                shared_template.as_deref(),
                &mut rng,
            );

            pairs.push((
                r1_record,
                r2_record,
                read_name,
                umi_str.clone(),
                mi_tag,
                a_is_top,
                insert_size,
            ));
        }

        // B reads: opposite orientation of A
        let b_is_top = !is_top_strand;
        for read_idx in 0..b_count {
            let read_name = format!("mol{mol_id:08}_readB{read_idx:04}");
            let mi_tag = format!("{mol_id}/B");

            let (r1_record, r2_record) = generate_read_pair_records(
                &read_name,
                &umi_str,
                &mi_tag,
                chrom_idx,
                local_pos,
                insert_size,
                params.read_length,
                params.mapq,
                b_is_top,
                &params.quality_model,
                &params.quality_bias,
                &params.methylation,
                shared_template.as_deref(),
                &mut rng,
            );

            pairs.push((
                r1_record,
                r2_record,
                read_name,
                umi_str.clone(),
                mi_tag,
                b_is_top,
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
                chrom_idx,
                local_pos,
                insert_size,
                params.read_length,
                params.mapq,
                is_top_strand,
                &params.quality_model,
                &params.quality_bias,
                &params.methylation,
                shared_template.as_deref(),
                &mut rng,
            );

            pairs.push((
                r1_record,
                r2_record,
                read_name,
                umi_str.clone(),
                mi_tag.clone(),
                is_top_strand,
                insert_size,
            ));
        }
    }

    (pairs, chrom_idx, local_pos)
}

#[allow(clippy::too_many_arguments)]
fn generate_read_pair_records(
    read_name: &str,
    umi_str: &str,
    mi_tag: &str,
    chrom_idx: usize,
    local_pos: usize,
    insert_size: usize,
    read_length: usize,
    mapq: u8,
    is_top_strand: bool,
    quality_model: &PositionQualityModel,
    quality_bias: &ReadPairQualityBias,
    methylation: &MethylationConfig,
    shared_template: Option<&[u8]>,
    rng: &mut impl Rng,
) -> (RecordBuf, RecordBuf) {
    // Use shared reference template or generate random per-read.
    // Avoid copying the shared template — only borrow it.
    let random_template;
    let template: &[u8] = if let Some(shared) = shared_template {
        shared
    } else {
        random_template = generate_random_sequence(insert_size, rng);
        &random_template
    };

    // Compute forward-read sequence (from start of template, top strand)
    let fwd_end = read_length.min(template.len());
    let mut fwd_seq: Vec<u8> = template[..fwd_end].to_vec();
    apply_methylation_conversion(
        &mut fwd_seq,
        template,
        0,
        true, // top strand
        methylation,
        rng,
    );
    let fwd_seq = pad_sequence(fwd_seq, read_length, rng);

    // Compute reverse-read sequence (from end of template, bottom strand, revcomped)
    let rev_start = insert_size.saturating_sub(read_length);
    let rev_end = read_length.min(template.len().saturating_sub(rev_start));
    let mut rev_template: Vec<u8> = template[rev_start..rev_start + rev_end].to_vec();
    apply_methylation_conversion(
        &mut rev_template,
        template,
        rev_start,
        false, // bottom strand
        methylation,
        rng,
    );
    let rev_seq = reverse_complement(&rev_template);
    let rev_seq = pad_sequence(rev_seq, read_length, rng);

    // Quality scores
    let r1_quals = quality_model.generate_qualities(read_length, rng);
    let r2_quals_raw = quality_model.generate_qualities(read_length, rng);
    let r2_quals = quality_bias.apply_to_vec(&r2_quals_raw, true);

    let mate_cigar = format!("{read_length}M");

    // Assign to R1/R2 based on strand coin flip.
    // tlen sign convention: positive for the leftmost read, negative for the rightmost.
    // When insert_size < read_length the reverse read starts at local_pos (mirroring
    // rev_start's `saturating_sub`), so use `saturating_sub` to avoid a usize underflow.
    let rev_pos = local_pos + insert_size.saturating_sub(read_length);
    let (r1_seq, r2_seq, r1_is_reverse, r1_pos, r2_pos, r1_tlen) = if is_top_strand {
        // F1R2: R1=forward at start, R2=reverse at end
        (fwd_seq, rev_seq, false, local_pos, rev_pos, insert_size as i32)
    } else {
        // R1F2: R1=reverse at end, R2=forward at start
        (rev_seq, fwd_seq, true, rev_pos, local_pos, -(insert_size as i32))
    };

    let r1_record = build_record(
        read_name,
        &r1_seq,
        &r1_quals,
        chrom_idx,
        r1_pos,
        mapq,
        true,          // is_first
        r1_is_reverse, // is_reverse
        r2_pos,
        r1_tlen,
        umi_str,
        mi_tag,
        &mate_cigar,
        mapq,
    );

    let r2_record = build_record(
        read_name,
        &r2_seq,
        &r2_quals,
        chrom_idx,
        r2_pos,
        mapq,
        false,          // is_first
        !r1_is_reverse, // is_reverse (opposite of R1)
        r1_pos,
        -r1_tlen,
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simulate::create_rng;

    const DISABLED_METHYLATION: MethylationConfig = MethylationConfig {
        mode: fgumi_consensus::MethylationMode::Disabled,
        cpg_methylation_rate: 0.75,
        conversion_rate: 0.98,
    };

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
            0,    // chrom_idx
            1000, // local_pos
            300,
            150,
            60,
            true, // is_top_strand
            &quality_model,
            &quality_bias,
            &DISABLED_METHYLATION,
            None,
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
            0,    // chrom_idx
            1000, // local_pos
            300,
            150,
            60,
            false, // is_top_strand (B strand = bottom)
            &quality_model,
            &quality_bias,
            &DISABLED_METHYLATION,
            None,
            &mut rng,
        );

        // Bottom strand has opposite orientation: R1 reverse, R2 forward
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
            0,    // chrom_idx
            1000, // local_pos
            300,
            150,
            60,
            true, // is_top_strand (simplex top)
            &quality_model,
            &quality_bias,
            &DISABLED_METHYLATION,
            None,
            &mut rng,
        );

        // Simplex top strand has same orientation as F1R2
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
            0,    // chrom_idx
            1000, // local_pos
            50,   // insert_size < read_length
            150,
            60,
            true, // is_top_strand
            &quality_model,
            &quality_bias,
            &DISABLED_METHYLATION,
            None,
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
            0,    // chrom_idx
            1000, // local_pos
            300,
            150,
            60,
            true, // is_top_strand
            &quality_model,
            &quality_bias,
            &DISABLED_METHYLATION,
            None,
            &mut rng1,
        );

        let (r1_b, r2_b) = generate_read_pair_records(
            "read_001",
            "ACGTACGT",
            "0/A",
            0,    // chrom_idx
            1000, // local_pos
            300,
            150,
            60,
            true, // is_top_strand
            &quality_model,
            &quality_bias,
            &DISABLED_METHYLATION,
            None,
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
        assert_eq!(
            record.alignment_start().expect("record should have alignment start").get(),
            1001
        );
        assert_eq!(
            record.mate_alignment_start().expect("record should have mate alignment start").get(),
            1101
        );
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
            assert_eq!(
                record.mapping_quality().expect("record should have mapping quality").get(),
                mapq
            );
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
            0,    // chrom_idx
            1000, // local_pos
            300,
            150,
            60,
            true, // is_top_strand
            &quality_model,
            &quality_bias,
            &DISABLED_METHYLATION,
            None,
            &mut rng,
        );

        // Top strand: R1 forward, R2 reverse
        assert!(r1.flags().is_first_segment());
        assert!(!r1.flags().is_reverse_complemented());
        assert!(r2.flags().is_last_segment());
        assert!(r2.flags().is_reverse_complemented());
    }

    // ========================================================================
    // Methylation tests
    // ========================================================================

    #[test]
    fn test_emseq_grouped_reads_strand_conversion() {
        // Top strand: R1=top (C->T), R2=bottom (G->A then RC)
        // Template: non-CpG Cs at known positions
        let template = b"CACACACACACACACAC"; // no CpG
        let config = MethylationConfig {
            mode: fgumi_consensus::MethylationMode::EmSeq,
            cpg_methylation_rate: 0.75,
            conversion_rate: 1.0,
        };

        let quality_model = PositionQualityModel::default();
        let quality_bias = ReadPairQualityBias::default();
        let mut rng = create_rng(Some(42));

        let (r1, _r2) = generate_read_pair_records(
            "test_meth",
            "AAAAAAAA",
            "1/A",
            0, // chrom_idx
            0, // local_pos
            template.len(),
            template.len(),
            60,
            true, // is_top_strand
            &quality_model,
            &quality_bias,
            &config,
            Some(template),
            &mut rng,
        );

        // R1 (top strand, forward) should have C->T conversions at non-CpG C positions
        let r1_seq: Vec<u8> = r1.sequence().as_ref().to_vec();
        for (i, &b) in template.iter().enumerate() {
            if i >= r1_seq.len() {
                break;
            }
            if b == b'C' {
                assert_eq!(r1_seq[i], b'T', "R1 position {i}: non-CpG C should convert to T");
            }
        }
    }

    #[test]
    fn test_bottom_strand_flips_conversion_orientation() {
        // Bottom strand: R1 gets reverse-complemented sequence from end of template
        // Top strand: R1 gets forward sequence from start of template
        // Template needs both C and G (non-CpG) so both strands have convertible bases
        let template = b"CATAGATAGATAGATA"; // has C and G, no CpG dinucleotides
        let emseq_config = MethylationConfig {
            mode: fgumi_consensus::MethylationMode::EmSeq,
            cpg_methylation_rate: 0.75,
            conversion_rate: 1.0,
        };
        let disabled_config = DISABLED_METHYLATION;

        let quality_model = PositionQualityModel::default();
        let quality_bias = ReadPairQualityBias::default();

        // Generate top strand pair (F1R2)
        let mut rng_a = create_rng(Some(42));
        let (r1_a, _) = generate_read_pair_records(
            "test",
            "AAAAAAAA",
            "1/A",
            0, // chrom_idx
            0, // local_pos
            template.len(),
            template.len(),
            60,
            true, // is_top_strand
            &quality_model,
            &quality_bias,
            &emseq_config,
            Some(template),
            &mut rng_a,
        );

        // Generate bottom strand pair with disabled methylation to check it differs
        let mut rng_b = create_rng(Some(42));
        let (r1_b_no_meth, _) = generate_read_pair_records(
            "test",
            "AAAAAAAA",
            "1/B",
            0, // chrom_idx
            0, // local_pos
            template.len(),
            template.len(),
            60,
            false, // is_top_strand (bottom strand)
            &quality_model,
            &quality_bias,
            &disabled_config,
            Some(template),
            &mut rng_b,
        );

        // Generate bottom strand pair with EM-Seq
        let mut rng_b2 = create_rng(Some(42));
        let (r1_b_meth, _) = generate_read_pair_records(
            "test",
            "AAAAAAAA",
            "1/B",
            0, // chrom_idx
            0, // local_pos
            template.len(),
            template.len(),
            60,
            false, // is_top_strand (bottom strand)
            &quality_model,
            &quality_bias,
            &emseq_config,
            Some(template),
            &mut rng_b2,
        );

        // Top strand R1 (forward, top) should differ from bottom strand R1 (reverse, revcomped)
        let r1_a_seq: Vec<u8> = r1_a.sequence().as_ref().to_vec();
        let r1_b_seq: Vec<u8> = r1_b_meth.sequence().as_ref().to_vec();
        // Bottom strand R1 is reverse-complemented from end of template with bottom-strand
        // methylation, so conversion pattern should differ from top strand R1
        assert_ne!(r1_a_seq, r1_b_seq, "Top and bottom strand R1 should differ");

        // Bottom strand with methylation should differ from bottom strand without
        let r1_b_no_meth_seq: Vec<u8> = r1_b_no_meth.sequence().as_ref().to_vec();
        assert_ne!(
            r1_b_seq, r1_b_no_meth_seq,
            "Bottom strand with/without methylation should differ"
        );
    }

    // ========================================================================
    // Real reference coordinate tests
    // ========================================================================

    #[test]
    fn test_grouped_reads_uses_real_ref_coordinates() {
        use std::io::Write as IoWrite;
        use tempfile::NamedTempFile;

        let mut fasta = NamedTempFile::new().unwrap();
        writeln!(fasta, ">chr1").unwrap();
        fasta.write_all(&b"ACGT".repeat(500)).unwrap();
        writeln!(fasta).unwrap();
        writeln!(fasta, ">chr2").unwrap();
        fasta.write_all(&b"CCGG".repeat(375)).unwrap();
        writeln!(fasta).unwrap();
        fasta.flush().unwrap();

        let ref_genome = ReferenceGenome::load(fasta.path()).unwrap();
        let mut rng = create_rng(Some(42));
        let positions = ref_genome.sample_positions(10, &mut rng);

        let params = GenerationParams {
            read_length: 50,
            umi_length: 8,
            mapq: 60,
            duplex: false,
            min_family_size: 1,
            quality_model: PositionQualityModel::default(),
            quality_bias: ReadPairQualityBias::default(),
            family_dist: FamilySizeDistribution::log_normal(1.0, 0.1),
            insert_model: InsertSizeModel::new(100.0, 5.0, 80, 120),
            strand_bias_model: StrandBiasModel::no_bias(),
            methylation: DISABLED_METHYLATION,
        };

        let (chrom_idx, local_pos) = positions[0];
        let (pairs, eff_chrom, _eff_pos) =
            generate_molecule_reads(0, 42, chrom_idx, local_pos, &params, &ref_genome);
        assert!(!pairs.is_empty());

        for (r1, r2, _, _, _, _, _) in &pairs {
            assert_eq!(r1.reference_sequence_id(), Some(eff_chrom));
            assert_eq!(r2.reference_sequence_id(), Some(eff_chrom));
        }
    }

    #[test]
    fn test_duplex_strand_coin_flip_produces_both_orientations_for_a_reads() {
        use std::io::Write as IoWrite;
        use tempfile::NamedTempFile;

        let mut fasta = NamedTempFile::new().unwrap();
        writeln!(fasta, ">chr1").unwrap();
        fasta.write_all(&b"ACGT".repeat(500)).unwrap();
        writeln!(fasta).unwrap();
        fasta.flush().unwrap();

        let ref_genome = ReferenceGenome::load(fasta.path()).unwrap();
        let params = GenerationParams {
            read_length: 50,
            umi_length: 8,
            mapq: 60,
            duplex: true,
            min_family_size: 2,
            quality_model: PositionQualityModel::default(),
            quality_bias: ReadPairQualityBias::default(),
            family_dist: FamilySizeDistribution::log_normal(2.0, 0.5),
            insert_model: InsertSizeModel::new(100.0, 5.0, 80, 120),
            strand_bias_model: StrandBiasModel::no_bias(),
            methylation: DISABLED_METHYLATION,
        };

        let mut a_saw_top = false;
        let mut a_saw_bottom = false;

        for seed in 0u64..100 {
            let (pairs, _chrom, _pos) =
                generate_molecule_reads(seed as usize, seed, 0, 500, &params, &ref_genome);

            for (_, _, _, _, mi_tag, is_top, _) in &pairs {
                if mi_tag.ends_with("/A") {
                    if *is_top {
                        a_saw_top = true;
                    } else {
                        a_saw_bottom = true;
                    }
                }
            }
        }

        assert!(a_saw_top, "A reads should sometimes be F1R2 (top strand)");
        assert!(a_saw_bottom, "A reads should sometimes be R1F2 (bottom strand)");
    }
}

//! Generate template-coordinate sorted BAM for group.
//!
//! Uses a memory-efficient approach: sort molecule IDs by their position first,
//! then generate records in sorted order (streaming). This avoids storing all
//! records in memory before sorting.

use super::sort::TemplateCoordKey;
use crate::bam_io::create_bam_writer;
use crate::commands::command::Command;
use crate::commands::common::CompressionOptions;
use crate::commands::simulate::common::{
    FamilySizeArgs, InsertSizeArgs, MethylationArgs, MethylationConfig, MoleculeInfo,
    PositionDistArgs, QualityArgs, ReferenceArgs, ReferenceGenome, SimulationCommon,
    apply_methylation_conversion, generate_random_sequence, pad_sequence, validate_rate,
};
use crate::dna::reverse_complement;
use crate::progress::ProgressTracker;
use crate::sam::builder::RecordBuilder;
use crate::simulate::{
    FamilySizeDistribution, InsertSizeModel, PositionQualityModel, ReadPairQualityBias, create_rng,
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

    /// Fraction of molecules (templates) to emit as entirely unmapped, in `[0.0, 1.0]`.
    ///
    /// For each molecule, a uniform random `f64` in `[0, 1)` is drawn from a
    /// dedicated RNG; when below this fraction, all read pairs of the molecule
    /// are emitted with `UNMAPPED` and `MATE_UNMAPPED` set on both mates (no
    /// reference, position, CIGAR, MAPQ, or MC/MQ tags), preserving pair
    /// consistency. Sequence, qualities, read names, and the RX tag are still
    /// populated so downstream tools see well-formed unmapped pairs.
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

    #[command(flatten)]
    pub methylation: MethylationArgs,
}

/// Parameters needed for molecule generation.
struct GenerationParams {
    read_length: usize,
    umi_length: usize,
    mapq: u8,
    min_family_size: usize,
    quality_model: PositionQualityModel,
    quality_bias: ReadPairQualityBias,
    family_dist: FamilySizeDistribution,
    insert_model: InsertSizeModel,
    methylation: MethylationConfig,
}

impl Command for MappedReads {
    fn execute(&self, command_line: &str) -> Result<()> {
        // Validate methylation args
        let methylation = self.methylation.resolve();
        self.methylation.validate()?;

        validate_rate(self.unmapped_fraction, "unmapped-fraction")?;

        info!("Generating mapped reads");
        info!("  Output: {}", self.output.display());
        info!("  Truth: {}", self.truth_output.display());
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

        // Re-add reference sequences from ref_genome header
        for (name, map) in ref_header.reference_sequences() {
            header_builder = header_builder.add_reference_sequence(name.clone(), map.clone());
        }

        // Add @PG record
        header_builder = crate::commands::common::add_pg_to_builder(header_builder, command_line)?;

        let header = header_builder.build();

        // Set up generation parameters
        let params = GenerationParams {
            read_length: self.common.read_length,
            umi_length: self.common.umi_length,
            mapq: self.mapq,
            min_family_size: self.family_size.min_family_size,
            quality_model: self.quality.to_quality_model(),
            quality_bias: self.quality.to_quality_bias(),
            family_dist: self.family_size.to_family_size_distribution()?,
            insert_model: self.insert_size.to_insert_size_model(),
            methylation,
        };

        // Pre-sample positions from reference (use a dedicated RNG so molecule seeds
        // stay deterministic and uncorrelated with position sampling)
        let num_positions = self.position_dist.num_positions.unwrap_or(self.common.num_molecules);
        if num_positions == 0 {
            anyhow::bail!("--num-positions must be greater than 0");
        }

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

        // Dedicated RNG for the per-molecule unmapped decision. Kept separate
        // from `seed_rng` and `mol_rng` so that toggling `--unmapped-fraction`
        // does not perturb UMI/sequence content for molecules that remain
        // mapped.
        let mut unmapped_rng = create_rng(self.common.seed.map(|s| s.wrapping_add(2)));

        // MEMORY-EFFICIENT APPROACH: Sort molecule IDs by template-coordinate key first,
        // then generate records in sorted order (streaming). This avoids storing all
        // records in memory before sorting.
        //
        // We need to pre-compute insert_size to get pos2 for the sort key, since
        // template-coordinate sorting considers both ends of the template.
        info!("Computing molecule positions and sort keys...");
        let mut molecules: Vec<MoleculeInfo> = (0..self.common.num_molecules)
            .map(|mol_id| {
                let seed: u64 = seed_rng.random();
                // Draw per-molecule unmapped decision from the dedicated RNG so it
                // stays deterministic across builds and uncorrelated with the other
                // streams. Always draw (no short-circuit on `unmapped_fraction`) so
                // the stream advances identically across molecules regardless of the
                // fraction value.
                let is_unmapped = unmapped_rng.random::<f64>() < self.unmapped_fraction;

                let sort_key = if is_unmapped {
                    TemplateCoordKey::for_unmapped_pair(format!("mol{mol_id:08}"))
                } else {
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
                    TemplateCoordKey::for_f1r2_pair(
                        chrom_idx as i32, // real tid
                        local_pos,
                        insert_size,
                        String::new(), // empty mid for mapped-reads (no MI tag yet)
                        format!("mol{mol_id:08}"),
                    )
                };

                MoleculeInfo { mol_id, seed, sort_key, is_unmapped }
            })
            .collect();

        // Sort molecules by template-coordinate key
        // (matching samtools bam1_cmp_template_coordinate)
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
                mol_info.is_unmapped,
                &params,
                &ref_genome,
            );

            for (r1, r2, read_name, umi_str, is_top_strand) in pairs {
                if mol_info.is_unmapped {
                    // Unmapped pairs: emit R1 then R2 in the conventional order.
                    // Truth rows record unmapped-sentinel values for chrom/pos/strand
                    // so downstream validators can distinguish unmapped from mapped
                    // molecules.
                    writer.write_alignment_record(&header, &r1)?;
                    writer.write_alignment_record(&header, &r2)?;
                    writeln!(
                        truth_writer,
                        "{}\t{}\t{}\t*\t*\t*",
                        read_name, umi_str, mol_info.mol_id
                    )?;
                } else {
                    // Template-coordinate order: write the record whose reference start
                    // is lower first so pairs are emitted in coordinate order. R1F2
                    // pairs have R1 (reverse read) at the higher coordinate, so emit R2
                    // first when R1's alignment_start is strictly greater.
                    let r1_pos = r1.alignment_start().map(|p| p.get()).unwrap_or(0);
                    let r2_pos = r2.alignment_start().map(|p| p.get()).unwrap_or(0);
                    let r2_first = r1_pos > r2_pos;
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
                        "{}\t{}\t{}\t{}\t{}\t{}",
                        read_name,
                        umi_str,
                        mol_info.mol_id,
                        ref_genome.name(chrom_idx),
                        local_pos,
                        strand_char,
                    )?;
                }
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

/// A single R1/R2 pair plus identifiers used for truth/metadata:
/// `(r1, r2, read_name, umi_str, is_top_strand)`.
type MoleculeReadPair = (RecordBuf, RecordBuf, String, String, bool);

/// Result of generating all read pairs for a molecule: the pairs plus the
/// effective `(chrom_idx, local_pos)` at which the template actually lives.
type MoleculeReadsResult = (Vec<MoleculeReadPair>, usize, usize);

/// Generate all read pairs for a single molecule.
///
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
    is_unmapped: bool,
    params: &GenerationParams,
    ref_genome: &ReferenceGenome,
) -> MoleculeReadsResult {
    let mut rng = create_rng(Some(seed));

    // Generate UMI
    let umi = generate_random_sequence(params.umi_length, &mut rng);
    let umi_str = String::from_utf8_lossy(&umi).to_string();

    // Generate family size
    let family_size = params.family_dist.sample(&mut rng, params.min_family_size);

    // Generate insert size (drawn even for unmapped molecules so the RNG
    // sequence matches the pre-pass insert_size pre-computation)
    let insert_size = params.insert_model.sample(&mut rng);

    // 50/50 strand coin flip: determines genomic strand of origin (unused for
    // unmapped pairs, but drawn to keep the RNG sequence identical to the
    // mapped case for a given seed)
    let is_top_strand: bool = rng.random();

    if is_unmapped {
        return generate_unmapped_molecule_reads(
            mol_id,
            chrom_idx,
            local_pos,
            &umi_str,
            family_size,
            params,
            &mut rng,
        );
    }

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

    let mut pairs = Vec::with_capacity(family_size);

    for read_idx in 0..family_size {
        let read_name = format!("mol{mol_id:08}_read{read_idx:04}");

        // Use shared reference template or generate random per-read.
        // Avoid cloning the shared template — only borrow it.
        let random_template;
        let template: &[u8] = if let Some(ref shared) = shared_template {
            shared
        } else {
            random_template = generate_random_sequence(insert_size, &mut rng);
            &random_template
        };

        // Compute forward-read sequence (from start of template, top strand)
        let fwd_end = params.read_length.min(template.len());
        let mut fwd_seq: Vec<u8> = template[..fwd_end].to_vec();
        apply_methylation_conversion(
            &mut fwd_seq,
            template,
            0,
            true, // top strand
            &params.methylation,
            &mut rng,
        );
        let fwd_seq = pad_sequence(fwd_seq, params.read_length, &mut rng);

        // Compute reverse-read sequence (from end of template, bottom strand,
        // reverse complemented)
        let rev_start = insert_size.saturating_sub(params.read_length);
        let rev_end = params.read_length.min(template.len().saturating_sub(rev_start));
        let mut rev_template: Vec<u8> = template[rev_start..rev_start + rev_end].to_vec();
        apply_methylation_conversion(
            &mut rev_template,
            template,
            rev_start,
            false, // bottom strand
            &params.methylation,
            &mut rng,
        );
        let rev_seq = reverse_complement(&rev_template);
        let rev_seq = pad_sequence(rev_seq, params.read_length, &mut rng);

        // Quality scores
        let r1_quals = params.quality_model.generate_qualities(params.read_length, &mut rng);
        let r2_quals_raw = params.quality_model.generate_qualities(params.read_length, &mut rng);
        let r2_quals = params.quality_bias.apply_to_vec(&r2_quals_raw, true);

        // Mate cigar is based on the mate's read length (same for both reads)
        let mate_cigar = format!("{}M", params.read_length);

        // Assign to R1/R2 based on strand coin flip.
        // tlen sign convention: positive for the leftmost read, negative for the
        // rightmost. When insert_size < read_length the reverse read starts at
        // local_pos (mirroring rev_start's `saturating_sub`), so use
        // `saturating_sub` here to avoid a usize underflow.
        let rev_pos = local_pos + insert_size.saturating_sub(params.read_length);
        let (r1_seq, r2_seq, r1_is_reverse, r1_pos, r2_pos, r1_tlen) = if is_top_strand {
            // F1R2: R1=forward at start, R2=reverse at end
            (fwd_seq, rev_seq, false, local_pos, rev_pos, insert_size as i32)
        } else {
            // R1F2: R1=reverse at end, R2=forward at start
            (rev_seq, fwd_seq, true, rev_pos, local_pos, -(insert_size as i32))
        };

        // Build R1 record
        let r1_record = build_record(
            &read_name,
            &r1_seq,
            &r1_quals,
            chrom_idx,
            r1_pos,
            params.mapq,
            true,          // is_first
            r1_is_reverse, // is_reverse
            r2_pos,
            r1_tlen,
            &umi_str,
            &mate_cigar,
            params.mapq,
        );

        // Build R2 record
        let r2_record = build_record(
            &read_name,
            &r2_seq,
            &r2_quals,
            chrom_idx,
            r2_pos,
            params.mapq,
            false,          // is_first
            !r1_is_reverse, // is_reverse (opposite of R1)
            r1_pos,
            -r1_tlen,
            &umi_str,
            &mate_cigar,
            params.mapq,
        );

        pairs.push((r1_record, r2_record, read_name, umi_str.clone(), is_top_strand));
    }

    (pairs, chrom_idx, local_pos)
}

/// Generate all read pairs for a molecule that has been flagged as unmapped.
///
/// Sequence content and qualities are still produced per-read (so BAMs contain
/// well-formed unmapped records) but no reference lookups are performed. Every
/// emitted pair has `UNMAPPED` and `MATE_UNMAPPED` set and carries only the RX
/// tag — no MC/MQ, no alignment position, no CIGAR.
fn generate_unmapped_molecule_reads(
    mol_id: usize,
    chrom_idx: usize,
    local_pos: usize,
    umi_str: &str,
    family_size: usize,
    params: &GenerationParams,
    rng: &mut impl Rng,
) -> MoleculeReadsResult {
    let mut pairs = Vec::with_capacity(family_size);

    for read_idx in 0..family_size {
        let read_name = format!("mol{mol_id:08}_read{read_idx:04}");

        let r1_seq_raw = generate_random_sequence(params.read_length, rng);
        let r1_seq = pad_sequence(r1_seq_raw, params.read_length, rng);
        let r2_seq_raw = generate_random_sequence(params.read_length, rng);
        let r2_seq = pad_sequence(r2_seq_raw, params.read_length, rng);

        let r1_quals = params.quality_model.generate_qualities(params.read_length, rng);
        let r2_quals_raw = params.quality_model.generate_qualities(params.read_length, rng);
        let r2_quals = params.quality_bias.apply_to_vec(&r2_quals_raw, true);

        let r1_record = build_unmapped_record(&read_name, &r1_seq, &r1_quals, true, umi_str);
        let r2_record = build_unmapped_record(&read_name, &r2_seq, &r2_quals, false, umi_str);

        // `is_top_strand=false` for unmapped pairs: there is no meaningful
        // strand of origin, so we report a consistent sentinel.
        pairs.push((r1_record, r2_record, read_name, umi_str.to_string(), false));
    }

    // Unmapped molecules have no meaningful genomic coordinates; return the
    // pre-sampled locus only so the tuple shape matches the mapped case. The
    // caller does not use these fields for unmapped pairs (truth rows emit
    // `*` sentinels instead).
    (pairs, chrom_idx, local_pos)
}

fn build_unmapped_record(
    name: &str,
    seq: &[u8],
    quals: &[u8],
    is_first: bool,
    umi: &str,
) -> RecordBuf {
    let seq_str = String::from_utf8_lossy(seq);

    // Omit cigar/reference/position/mate-ref/mate-pos/MAPQ/TLEN/MC/MQ; only
    // sequence, qualities, PAIRED/FIRST|LAST/UNMAPPED/MATE_UNMAPPED flags,
    // and the RX tag are valid on an unmapped pair. `.cigar("")` forces an
    // empty CIGAR — without it, `RecordBuilder` auto-generates `{len}M` from
    // the sequence.
    // `RecordBuilder::new()` defaults MAPQ to 60, which is wrong for an unmapped
    // read (SAM spec §1.4.5: MAPQ on unmapped reads should be 0). Override.
    RecordBuilder::new()
        .name(name)
        .sequence(&seq_str)
        .qualities(quals)
        .cigar("")
        .mapping_quality(0)
        .paired(true)
        .first_segment(is_first)
        .unmapped(true)
        .mate_unmapped(true)
        .tag("RX", umi)
        .build()
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simulate::create_rng;

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
        assert_eq!(reverse_complement(b"X"), b"X");
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

        // Note: mapq 255 is treated as "unavailable" by noodles, so skip it
        for mapq in [0, 30, 60] {
            let record = build_record(
                "read", seq, &quals, 0, 100, mapq, true, false, 200, 100, "AAA", "4M", mapq,
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

    // ========================================================================
    // Methylation tests
    // ========================================================================

    #[test]
    fn test_emseq_mapped_reads_conversion() {
        // Template with non-CpG Cs that should convert C->T in EM-Seq
        let template = b"CACACACACACACACAC".to_vec(); // no CpG dinucleotides
        let config = MethylationConfig {
            mode: fgumi_consensus::MethylationMode::EmSeq,
            cpg_methylation_rate: 0.75,
            conversion_rate: 1.0,
        };

        let mut r1 = template[..8].to_vec();
        let mut rng = create_rng(Some(42));
        apply_methylation_conversion(&mut r1, &template, 0, true, &config, &mut rng);

        // All non-CpG Cs should convert to T
        for (i, &b) in r1.iter().enumerate() {
            if template[i] == b'C' {
                assert_eq!(b, b'T', "position {i}: non-CpG C should convert");
            }
        }
    }

    #[test]
    fn test_reads_in_family_share_template_with_reference() {
        // With a reference, all reads in a family should start from the same template
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut fasta = NamedTempFile::new().unwrap();
        writeln!(fasta, ">chr1").unwrap();
        let seq = b"ACGT".repeat(500); // 2000 bp
        fasta.write_all(&seq).unwrap();
        writeln!(fasta).unwrap();
        fasta.flush().unwrap();

        let ref_genome = ReferenceGenome::load(fasta.path()).unwrap();
        let params = GenerationParams {
            read_length: 50,
            umi_length: 8,
            mapq: 60,
            min_family_size: 3,
            quality_model: crate::simulate::PositionQualityModel::new(
                10, 25, 37, 100, 0.08, 2, 0.0,
            ),
            quality_bias: crate::simulate::ReadPairQualityBias::new(0),
            family_dist: crate::simulate::FamilySizeDistribution::log_normal(5.0, 1.0),
            insert_model: crate::simulate::InsertSizeModel::new(100.0, 10.0, 80, 120),
            methylation: MethylationConfig {
                mode: fgumi_consensus::MethylationMode::Disabled,
                cpg_methylation_rate: 0.75,
                conversion_rate: 0.98,
            },
        };

        let (pairs, _chrom, _pos) =
            generate_molecule_reads(0, 42, 0, 500, false, &params, &ref_genome);
        // Multiple reads should be generated (family_size >= 3)
        assert!(pairs.len() >= 3, "Expected at least 3 reads, got {}", pairs.len());
    }

    // ========================================================================
    // Real-coordinate and strand tests
    // ========================================================================

    #[test]
    fn test_mapped_reads_uses_real_ref_coordinates() {
        use std::io::Write as IoWrite;
        use tempfile::NamedTempFile;

        // Build multi-contig test reference
        let mut fasta = NamedTempFile::new().unwrap();
        writeln!(fasta, ">chr1").unwrap();
        fasta.write_all(&b"ACGT".repeat(500)).unwrap(); // 2000bp
        writeln!(fasta).unwrap();
        writeln!(fasta, ">chr2").unwrap();
        fasta.write_all(&b"CCGG".repeat(375)).unwrap(); // 1500bp
        writeln!(fasta).unwrap();
        writeln!(fasta, ">chr3").unwrap();
        fasta.write_all(&b"AATT".repeat(450)).unwrap(); // 1800bp
        writeln!(fasta).unwrap();
        fasta.flush().unwrap();

        let ref_genome = ReferenceGenome::load(fasta.path()).unwrap();
        let mut rng = create_rng(Some(42));
        let positions = ref_genome.sample_positions(10, &mut rng);

        let params = GenerationParams {
            read_length: 50,
            umi_length: 8,
            mapq: 60,
            min_family_size: 1,
            quality_model: crate::simulate::PositionQualityModel::new(
                10, 25, 37, 100, 0.08, 2, 0.0,
            ),
            quality_bias: crate::simulate::ReadPairQualityBias::new(0),
            family_dist: crate::simulate::FamilySizeDistribution::log_normal(5.0, 1.0),
            insert_model: crate::simulate::InsertSizeModel::new(100.0, 10.0, 80, 120),
            methylation: MethylationConfig {
                mode: fgumi_consensus::MethylationMode::Disabled,
                cpg_methylation_rate: 0.75,
                conversion_rate: 0.98,
            },
        };

        // Test that generate_molecule_reads produces records with correct ref_id
        let (chrom_idx, local_pos) = positions[0];
        let (pairs, eff_chrom, eff_pos) =
            generate_molecule_reads(0, 42, chrom_idx, local_pos, false, &params, &ref_genome);
        assert!(!pairs.is_empty());

        for (r1, r2, _, _, _) in &pairs {
            // Both reads should reference the effective chromosome (primary or fallback)
            assert_eq!(r1.reference_sequence_id(), Some(eff_chrom));
            assert_eq!(r2.reference_sequence_id(), Some(eff_chrom));

            // R1 or R2 should be at eff_pos (depending on strand)
            let r1_pos = r1.alignment_start().map(|p| p.get()).unwrap_or(0);
            let r2_pos = r2.alignment_start().map(|p| p.get()).unwrap_or(0);
            // One of the reads should be at eff_pos + 1 (1-based)
            assert!(
                r1_pos == eff_pos + 1 || r2_pos == eff_pos + 1,
                "Expected one read at position {} (1-based), got r1={}, r2={}",
                eff_pos + 1,
                r1_pos,
                r2_pos
            );
        }
    }

    #[test]
    fn test_strand_coin_flip_produces_both_orientations() {
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
            min_family_size: 1,
            quality_model: crate::simulate::PositionQualityModel::new(
                10, 25, 37, 100, 0.08, 2, 0.0,
            ),
            quality_bias: crate::simulate::ReadPairQualityBias::new(0),
            family_dist: crate::simulate::FamilySizeDistribution::log_normal(1.0, 0.1),
            insert_model: crate::simulate::InsertSizeModel::new(100.0, 5.0, 80, 120),
            methylation: MethylationConfig {
                mode: fgumi_consensus::MethylationMode::Disabled,
                cpg_methylation_rate: 0.75,
                conversion_rate: 0.98,
            },
        };

        let mut saw_f1r2 = false;
        let mut saw_r1f2 = false;

        // Generate many molecules with different seeds to verify both orientations
        for seed in 0u64..100 {
            let (pairs, _chrom, _pos) =
                generate_molecule_reads(seed as usize, seed, 0, 500, false, &params, &ref_genome);
            for (r1, _, _, _, is_top) in &pairs {
                if *is_top {
                    saw_f1r2 = true;
                    assert!(!r1.flags().is_reverse_complemented(), "F1R2: R1 should be forward");
                } else {
                    saw_r1f2 = true;
                    assert!(r1.flags().is_reverse_complemented(), "R1F2: R1 should be reverse");
                }
            }
        }

        assert!(saw_f1r2, "Should see at least one F1R2 molecule");
        assert!(saw_r1f2, "Should see at least one R1F2 molecule");
    }

    // ========================================================================
    // --unmapped-fraction tests (issue #302)
    // ========================================================================

    fn basic_unmapped_test_params(min_family_size: usize) -> GenerationParams {
        GenerationParams {
            read_length: 50,
            umi_length: 8,
            mapq: 60,
            min_family_size,
            quality_model: crate::simulate::PositionQualityModel::new(
                10, 25, 37, 100, 0.08, 2, 0.0,
            ),
            quality_bias: crate::simulate::ReadPairQualityBias::new(0),
            family_dist: crate::simulate::FamilySizeDistribution::log_normal(5.0, 1.0),
            insert_model: crate::simulate::InsertSizeModel::new(100.0, 10.0, 80, 120),
            methylation: MethylationConfig {
                mode: fgumi_consensus::MethylationMode::Disabled,
                cpg_methylation_rate: 0.75,
                conversion_rate: 0.98,
            },
        }
    }

    /// When a molecule is flagged as unmapped, every emitted R1/R2 pair must
    /// have `UNMAPPED` and `MATE_UNMAPPED` set, no reference/position/CIGAR/MAPQ,
    /// and carry the same read-name, sequence length, and RX tag as a mapped
    /// pair would.
    #[test]
    fn test_unmapped_molecule_produces_unmapped_records() {
        use std::io::Write as IoWrite;
        use tempfile::NamedTempFile;

        let mut fasta = NamedTempFile::new().unwrap();
        writeln!(fasta, ">chr1").unwrap();
        fasta.write_all(&b"ACGT".repeat(500)).unwrap();
        writeln!(fasta).unwrap();
        fasta.flush().unwrap();

        let ref_genome = ReferenceGenome::load(fasta.path()).unwrap();
        let params = basic_unmapped_test_params(3);

        let (pairs, _chrom, _pos) =
            generate_molecule_reads(0, 42, 0, 500, true, &params, &ref_genome);
        assert!(!pairs.is_empty(), "molecule should still emit pairs when unmapped");

        for (r1, r2, _, umi, _) in &pairs {
            for r in [r1, r2] {
                let flags = r.flags();
                assert!(flags.is_segmented(), "PAIRED flag should be set");
                assert!(flags.is_unmapped(), "UNMAPPED flag should be set");
                assert!(flags.is_mate_unmapped(), "MATE_UNMAPPED flag should be set");
                assert!(r.reference_sequence_id().is_none(), "ref_id should be unset");
                assert!(r.alignment_start().is_none(), "pos should be unset");
                assert!(r.mate_reference_sequence_id().is_none(), "mate_ref_id should be unset");
                assert!(r.mate_alignment_start().is_none(), "mate_pos should be unset");
                // SAM spec §1.4.5: MAPQ on unmapped reads should be 0.
                assert_eq!(
                    r.mapping_quality().map(|m| m.get()),
                    Some(0),
                    "MAPQ should be 0 on unmapped reads"
                );
                assert_eq!(r.template_length(), 0, "TLEN should be 0");
                assert_eq!(r.cigar().as_ref().len(), 0, "CIGAR should be empty");
                assert_eq!(r.sequence().len(), params.read_length, "sequence length preserved");
                assert_eq!(
                    r.quality_scores().as_ref().len(),
                    params.read_length,
                    "qualities length preserved"
                );
            }
            assert!(r1.flags().is_first_segment(), "R1 should be FIRST_SEGMENT");
            assert!(r2.flags().is_last_segment(), "R2 should be LAST_SEGMENT");
            assert_eq!(umi.len(), params.umi_length, "UMI length preserved");
        }
    }

    /// Unmapped molecules must not carry a MC/MQ tag; there is no mate
    /// alignment to describe.
    #[test]
    fn test_unmapped_molecule_has_no_mate_alignment_tags() {
        use noodles::sam::alignment::record::data::field::Tag;
        use std::io::Write as IoWrite;
        use tempfile::NamedTempFile;

        let mut fasta = NamedTempFile::new().unwrap();
        writeln!(fasta, ">chr1").unwrap();
        fasta.write_all(&b"ACGT".repeat(500)).unwrap();
        writeln!(fasta).unwrap();
        fasta.flush().unwrap();

        let ref_genome = ReferenceGenome::load(fasta.path()).unwrap();
        let params = basic_unmapped_test_params(2);

        let (pairs, _chrom, _pos) =
            generate_molecule_reads(1, 99, 0, 500, true, &params, &ref_genome);
        assert!(!pairs.is_empty());

        let mc = Tag::from([b'M', b'C']);
        let mq = Tag::from([b'M', b'Q']);
        let rx = Tag::from([b'R', b'X']);
        for (r1, r2, _, _, _) in &pairs {
            for r in [r1, r2] {
                assert!(r.data().get(&mc).is_none(), "MC tag must not be present");
                assert!(r.data().get(&mq).is_none(), "MQ tag must not be present");
                assert!(r.data().get(&rx).is_some(), "RX tag must still be present");
            }
        }
    }

    /// Passing `is_unmapped=false` must produce records identical to the
    /// pre-fix (mapped) behavior: reference, position, CIGAR, and MAPQ are
    /// all populated.
    #[test]
    fn test_mapped_path_unchanged_when_is_unmapped_false() {
        use std::io::Write as IoWrite;
        use tempfile::NamedTempFile;

        let mut fasta = NamedTempFile::new().unwrap();
        writeln!(fasta, ">chr1").unwrap();
        fasta.write_all(&b"ACGT".repeat(500)).unwrap();
        writeln!(fasta).unwrap();
        fasta.flush().unwrap();

        let ref_genome = ReferenceGenome::load(fasta.path()).unwrap();
        let params = basic_unmapped_test_params(1);

        let (pairs, _chrom, _pos) =
            generate_molecule_reads(0, 42, 0, 500, false, &params, &ref_genome);
        assert!(!pairs.is_empty());

        for (r1, r2, _, _, _) in &pairs {
            for r in [r1, r2] {
                let flags = r.flags();
                assert!(!flags.is_unmapped(), "mapped path should not set UNMAPPED");
                assert!(!flags.is_mate_unmapped(), "mapped path should not set MATE_UNMAPPED");
                assert!(r.reference_sequence_id().is_some());
                assert!(r.alignment_start().is_some());
                assert!(r.mapping_quality().is_some());
                assert_ne!(r.cigar().as_ref().len(), 0, "mapped path should populate CIGAR");
            }
        }
    }

    /// Helper to build a minimal reference FASTA on disk for end-to-end tests.
    fn write_minimal_fasta(dir: &std::path::Path) -> std::path::PathBuf {
        use std::io::Write as IoWrite;
        let fasta = dir.join("ref.fa");
        let mut f = std::fs::File::create(&fasta).unwrap();
        writeln!(f, ">chr1").unwrap();
        f.write_all(&b"ACGT".repeat(1500)).unwrap(); // 6000 bp
        writeln!(f).unwrap();
        fasta
    }

    /// End-to-end: running `mapped-reads --unmapped-fraction 1.0` must produce
    /// a BAM in which every record has the UNMAPPED flag set.
    #[test]
    fn test_execute_unmapped_fraction_one_produces_all_unmapped() {
        use clap::Parser;
        use tempfile::tempdir;

        let dir = tempdir().unwrap();
        let fasta = write_minimal_fasta(dir.path());
        let out_bam = dir.path().join("out.bam");
        let truth = dir.path().join("truth.tsv");

        let cmd = MappedReads::try_parse_from([
            "mapped-reads",
            "-o",
            out_bam.to_str().unwrap(),
            "--truth",
            truth.to_str().unwrap(),
            "-r",
            fasta.to_str().unwrap(),
            "--num-molecules",
            "20",
            "--unmapped-fraction",
            "1.0",
            "--seed",
            "42",
        ])
        .expect("CLI parse");

        cmd.execute("test").expect("execute() succeeded");

        let mut reader = noodles::bam::io::reader::Builder.build_from_path(&out_bam).unwrap();
        let header = reader.read_header().unwrap();
        let mut n_records = 0;
        for result in reader.records() {
            let record = result.unwrap();
            assert!(
                record.flags().is_unmapped(),
                "every record must be UNMAPPED when --unmapped-fraction=1.0"
            );
            assert!(
                record.flags().is_mate_unmapped(),
                "every record must be MATE_UNMAPPED when --unmapped-fraction=1.0"
            );
            n_records += 1;
            // suppress unused-header warning
            let _ = &header;
        }
        assert!(n_records > 0, "expected some records, got 0");
    }

    /// End-to-end: running `mapped-reads --unmapped-fraction 0.0` (the default)
    /// must produce a BAM with zero unmapped records — the original mapped
    /// behavior is preserved.
    #[test]
    fn test_execute_unmapped_fraction_zero_produces_no_unmapped() {
        use clap::Parser;
        use tempfile::tempdir;

        let dir = tempdir().unwrap();
        let fasta = write_minimal_fasta(dir.path());
        let out_bam = dir.path().join("out.bam");
        let truth = dir.path().join("truth.tsv");

        let cmd = MappedReads::try_parse_from([
            "mapped-reads",
            "-o",
            out_bam.to_str().unwrap(),
            "--truth",
            truth.to_str().unwrap(),
            "-r",
            fasta.to_str().unwrap(),
            "--num-molecules",
            "20",
            "--unmapped-fraction",
            "0.0",
            "--seed",
            "42",
        ])
        .expect("CLI parse");

        cmd.execute("test").expect("execute() succeeded");

        let mut reader = noodles::bam::io::reader::Builder.build_from_path(&out_bam).unwrap();
        let _ = reader.read_header().unwrap();
        for result in reader.records() {
            let record = result.unwrap();
            assert!(
                !record.flags().is_unmapped(),
                "no record should be UNMAPPED when --unmapped-fraction=0.0"
            );
        }
    }

    /// End-to-end: with a fractional `--unmapped-fraction`, the BAM must
    /// contain BOTH mapped and unmapped records, and all unmapped records
    /// must sort after all mapped records (matching the `for_unmapped_pair`
    /// sort-key convention).
    #[test]
    fn test_execute_unmapped_fraction_mixed_sorts_unmapped_last() {
        use clap::Parser;
        use tempfile::tempdir;

        let dir = tempdir().unwrap();
        let fasta = write_minimal_fasta(dir.path());
        let out_bam = dir.path().join("out.bam");
        let truth = dir.path().join("truth.tsv");

        let cmd = MappedReads::try_parse_from([
            "mapped-reads",
            "-o",
            out_bam.to_str().unwrap(),
            "--truth",
            truth.to_str().unwrap(),
            "-r",
            fasta.to_str().unwrap(),
            "--num-molecules",
            "50",
            "--unmapped-fraction",
            "0.5",
            "--seed",
            "42",
        ])
        .expect("CLI parse");

        cmd.execute("test").expect("execute() succeeded");

        let mut reader = noodles::bam::io::reader::Builder.build_from_path(&out_bam).unwrap();
        let _ = reader.read_header().unwrap();
        let mut mapped_count = 0;
        let mut unmapped_count = 0;
        let mut first_unmapped_idx: Option<usize> = None;
        let mut last_mapped_idx: Option<usize> = None;
        for (idx, result) in reader.records().enumerate() {
            let record = result.unwrap();
            if record.flags().is_unmapped() {
                unmapped_count += 1;
                first_unmapped_idx.get_or_insert(idx);
            } else {
                mapped_count += 1;
                last_mapped_idx = Some(idx);
            }
        }
        assert!(mapped_count > 0, "expected some mapped records, got {mapped_count}");
        assert!(unmapped_count > 0, "expected some unmapped records, got {unmapped_count}");
        if let (Some(first_unmapped), Some(last_mapped)) = (first_unmapped_idx, last_mapped_idx) {
            assert!(
                last_mapped < first_unmapped,
                "all mapped records must sort before any unmapped record \
                 (last mapped idx {last_mapped}, first unmapped idx {first_unmapped})"
            );
        }
    }

    /// `--unmapped-fraction` outside `[0.0, 1.0]` must be rejected with a
    /// descriptive error, not silently ignored.
    #[rstest::rstest]
    #[case("-0.1")]
    #[case("1.5")]
    #[case("NaN")]
    fn test_execute_rejects_out_of_range_unmapped_fraction(#[case] bad: &str) {
        use clap::Parser;
        use tempfile::tempdir;

        let dir = tempdir().unwrap();
        let fasta = write_minimal_fasta(dir.path());
        let out_bam = dir.path().join("out.bam");
        let truth = dir.path().join("truth.tsv");

        // Use `--flag=value` so clap doesn't mis-parse negative values as
        // separate flags.
        let flag_arg = format!("--unmapped-fraction={bad}");
        let cmd = MappedReads::try_parse_from([
            "mapped-reads",
            "-o",
            out_bam.to_str().unwrap(),
            "--truth",
            truth.to_str().unwrap(),
            "-r",
            fasta.to_str().unwrap(),
            "--num-molecules",
            "5",
            &flag_arg,
            "--seed",
            "1",
        ])
        .expect("CLI parse");

        let result = cmd.execute("test");
        assert!(result.is_err(), "value {bad} should be rejected");
        let msg = result.unwrap_err().to_string();
        assert!(
            msg.contains("--unmapped-fraction"),
            "error message should mention the flag name, got: {msg}"
        );
    }
}

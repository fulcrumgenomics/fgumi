//! Review consensus variant calls by extracting supporting reads.
//!
//! This tool extracts consensus reads containing variants and their supporting
//! raw reads to facilitate manual review of variant calls. It creates filtered
//! BAM files and a detailed TSV report.

use crate::commands::common::parse_bool;
use crate::logging::OperationTimer;
use crate::reference::find_dict_path;
use crate::sam::SamTag;
use crate::umi::extract_mi_base;
use crate::validation::validate_file_exists;
use crate::variant_review::{
    BaseCounts, ConsensusVariantReviewInfo, Variant, format_insert_string, read_number_suffix,
};
use anyhow::{Result, bail};
use clap::Parser;
use fgumi_bam_io::ProgressTracker;
use fgumi_raw_bam::{
    BAM_BASE_TO_ASCII, CigarKind, IndexedRawBamReader, RawBamReader, RawRecord,
    find_string_tag_in_record, write_raw_record,
};
use log::info;
use std::collections::HashSet;
use std::path::{Path, PathBuf};

use super::command::Command;

/// Extracts the source molecule id from a record's `MI` tag, truncated at the
/// last `/` (see [`extract_mi_base`]).
///
/// Mirrors fgbio `ReviewConsensusVariants.toMi`: a record that reaches this
/// point is expected to carry an `MI` tag, and it is an error if it does not
/// (fgbio throws `IllegalStateException`; `ReviewConsensusVariantsTest.scala:166`).
fn to_mi(record: &RawRecord) -> Result<String> {
    match find_string_tag_in_record(record, SamTag::MI) {
        Some(mi_bytes) => Ok(extract_mi_base(std::str::from_utf8(mi_bytes)?).to_string()),
        None => {
            bail!("{} did not have a value for tag MI", String::from_utf8_lossy(record.read_name()))
        }
    }
}

/// Formats a genotype the way htsjdk's `Genotype.getGenotypeString` does, so the
/// `genotype` column matches fgbio's (`variant.genotype.getOrElse("NA")`).
///
/// Each allele position is mapped to its base string (`0` => reference,
/// `i` => `i`-th alternate, missing => `.`) and the base strings are joined in
/// genotype order. htsjdk preserves allele order (it does not sort, so `1/0`
/// renders as `T/A`, not `A/T`); only the separator differs — `|` for a fully
/// phased genotype, otherwise `/`.
fn format_genotype_string(
    alleles: &[noodles::vcf::variant::record_buf::samples::sample::value::genotype::Allele],
    reference: &str,
    alternates: &[String],
) -> String {
    use noodles::vcf::variant::record::samples::series::value::genotype::Phasing;

    let strings: Vec<String> = alleles
        .iter()
        .map(|allele| match allele.position() {
            None => ".".to_string(),
            Some(0) => reference.to_string(),
            Some(i) => alternates.get(i - 1).cloned().unwrap_or_else(|| ".".to_string()),
        })
        .collect();

    // htsjdk treats a genotype as phased only when every separator is `|`, which
    // noodles records as `Phased` on each allele after the first.
    let phased = alleles.len() > 1 && alleles[1..].iter().all(|a| a.phasing() == Phasing::Phased);
    let separator = if phased { "|" } else { "/" };

    strings.join(separator)
}

/// Reviews consensus variant calls by extracting relevant reads.
///
/// Creates filtered BAM files containing consensus reads with variants
/// and their supporting raw reads, plus a detailed review TSV file.
#[derive(Parser, Debug)]
#[command(
    name = "review",
    author,
    version,
    about = "\x1b[38;5;173m[POST-CONSENSUS]\x1b[0m \x1b[36mExtract data to review variant calls from consensus reads\x1b[0m",
    long_about = r#"
Extracts data to make reviewing of variant calls from consensus reads easier.

Creates a list of variant sites from the input VCF (SNPs only) or IntervalList then extracts all
the consensus reads that do not contain a reference allele at the variant sites, and all raw reads
that contributed to those consensus reads. This will include consensus reads that carry the
alternate allele, a third allele, a no-call or a spanning deletion at the variant site.

Reads are correlated between consensus and grouped BAMs using a molecule ID stored in an optional
attribute, `MI` by default. In order to support paired molecule IDs where two or more molecule IDs
are related (e.g. see the Paired assignment strategy in `group`) the molecule ID is truncated at
the last `/` if present (e.g. `1/A => 1` and `2 => 2`).

Both input BAMs must be coordinate sorted and indexed.

## Output Files

A pair of output BAMs are created:
- **<output>.consensus.bam**: Contains the relevant consensus reads from the consensus BAM
- **<output>.grouped.bam**: Contains the relevant raw reads from the grouped BAM

A review file **<output>.txt** is also created. The review file contains details on each variant
position along with detailed information on each consensus read that supports the variant. If the
`--sample` argument is supplied and the input is VCF, genotype information for that sample will be
retrieved. If the sample name isn't supplied and the VCF contains only a single sample then those
genotypes will be used.

The `--maf` parameter controls which variants get detailed per-read information in the output file.
Only variants with a minor allele frequency at or below this threshold will have detailed information
written.
"#
)]
pub struct Review {
    /// Input VCF or `IntervalList` of variant locations
    #[arg(short = 'i', long = "input")]
    pub input: PathBuf,

    /// BAM file of consensus reads used to call variants
    #[arg(short = 'c', long = "consensus-bam")]
    pub consensus_bam: PathBuf,

    /// BAM file of grouped raw reads used to build consensuses
    #[arg(short = 'g', long = "grouped-bam")]
    pub grouped_bam: PathBuf,

    /// Reference FASTA file
    #[arg(short = 'r', long = "ref")]
    pub reference: PathBuf,

    /// Output prefix for generated files
    #[arg(short = 'o', long = "output")]
    pub output: PathBuf,

    /// Name of sample being reviewed (for VCF genotype extraction)
    #[arg(short = 's', long = "sample")]
    pub sample: Option<String>,

    /// Ignore N bases in consensus reads
    #[arg(short = 'N', long = "ignore-ns", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool)]
    pub ignore_ns: bool,

    /// Only output detailed information for variants at or below this MAF
    #[arg(short = 'm', long = "maf", default_value = "0.05")]
    pub maf: f64,
}

impl Command for Review {
    fn execute(&self, command_line: &str) -> Result<()> {
        info!("Review");
        info!("  Input: {}", self.input.display());
        info!("  Consensus BAM: {}", self.consensus_bam.display());
        info!("  Grouped BAM: {}", self.grouped_bam.display());
        info!("  Reference: {}", self.reference.display());
        info!("  Output prefix: {}", self.output.display());
        info!("  Ignore Ns: {}", self.ignore_ns);
        info!("  MAF threshold: {}", self.maf);

        let timer = OperationTimer::new("Reviewing variants");

        // Validate inputs
        validate_file_exists(&self.input, "input file")?;
        validate_file_exists(&self.consensus_bam, "consensus BAM")?;
        validate_file_exists(&self.grouped_bam, "grouped BAM")?;
        validate_file_exists(&self.reference, "reference file")?;

        // Validate FASTA index and dictionary
        self.validate_reference_files()?;

        // Determine file type and load variants
        info!("Loading variants...");
        let variants = if self.is_vcf_file(&self.input) {
            self.load_variants_from_vcf()?
        } else {
            self.load_variants_from_intervals()?
        };

        info!("Loaded {} variant positions", variants.len());

        // Order variants by the reference sequence dictionary so the review file's rows
        // are emitted in coordinate order, matching fgbio (issue #497). The consensus BAM
        // header carries the same sequences, in the same order, as the reference `.dict`.
        let variants = {
            use noodles::bam;
            let mut reader = bam::io::reader::Builder.build_from_path(&self.consensus_bam)?;
            let header = reader.read_header()?;
            Self::order_variants_by_dictionary(variants, &header)?
        };

        // Extract consensus reads with variants
        info!("Extracting consensus reads with variants...");
        let mi_set = self.extract_consensus_reads(&variants, command_line)?;
        info!("Found {} unique molecular identifiers", mi_set.len());

        // Extract grouped reads matching those MIs
        info!("Extracting grouped reads...");
        self.extract_grouped_reads(&mi_set, command_line)?;

        // Create BAM indexes for output files
        info!("Creating BAM indexes...");
        let consensus_out_path = self.output.with_extension("consensus.bam");
        let grouped_out_path = self.output.with_extension("grouped.bam");

        use noodles::bam;
        use noodles::bam::bai;
        let consensus_index = bam::fs::index(&consensus_out_path)?;
        let mut consensus_index_writer = bai::io::Writer::new(std::fs::File::create(
            consensus_out_path.with_extension("bam.bai"),
        )?);
        consensus_index_writer.write_index(&consensus_index)?;

        let grouped_index = bam::fs::index(&grouped_out_path)?;
        let mut grouped_index_writer = bai::io::Writer::new(std::fs::File::create(
            grouped_out_path.with_extension("bam.bai"),
        )?);
        grouped_index_writer.write_index(&grouped_index)?;

        // Generate detailed review file
        info!("Generating detailed review file...");
        self.generate_review_file(&variants)?;

        info!("Done!");
        timer.log_completion(variants.len() as u64);
        Ok(())
    }
}

impl Review {
    /// Load a BAM index, preferring the `<prefix>.bam.bai` sidecar and falling back
    /// to `<prefix>.bai` for callers that used samtools' default naming.
    fn read_bam_index(path: &Path) -> Result<noodles::bam::bai::Index> {
        use noodles::bam;

        let bam_bai = path.with_extension("bam.bai");
        if bam_bai.exists() {
            return Ok(bam::bai::fs::read(&bam_bai)?);
        }

        let bai = path.with_extension("bai");
        if bai.exists() {
            return Ok(bam::bai::fs::read(&bai)?);
        }

        bail!(
            "Missing BAM index for {}. Tried {} and {}",
            path.display(),
            bam_bai.display(),
            bai.display()
        );
    }

    /// Orders variants by the reference sequence dictionary (coordinate order).
    ///
    /// fgbio processes variants by zipping `variants.iterator` with a
    /// coordinate-ordered `SamLocusIterator` and `require`s the two stay in sync
    /// (`ReviewConsensusVariants.scala:241`), so its review rows are emitted in
    /// coordinate order and it errors on out-of-order input. fgumi previously used
    /// input-file order, diverging for non-coordinate-sorted variants (issue #497).
    /// Sorting here reproduces fgbio's row order regardless of the input order.
    ///
    /// The dictionary order is taken from `header`'s reference sequences, which match
    /// the reference `.dict` fgbio reads.
    ///
    /// # Errors
    ///
    /// Returns an error if a variant's contig is absent from the sequence dictionary
    /// (fgbio likewise requires each variant interval's contig to exist in the dict).
    fn order_variants_by_dictionary(
        variants: Vec<Variant>,
        header: &noodles::sam::Header,
    ) -> Result<Vec<Variant>> {
        use std::collections::HashMap;

        let order: HashMap<String, usize> = header
            .reference_sequences()
            .keys()
            .enumerate()
            .map(|(idx, name)| (String::from_utf8_lossy(name).to_string(), idx))
            .collect();

        for variant in &variants {
            if !order.contains_key(&variant.chrom) {
                bail!(
                    "Variant contig '{}' is not present in the consensus BAM sequence dictionary",
                    variant.chrom
                );
            }
        }

        // Stable sort by (dictionary index, position); ties keep input order, matching
        // fgbio's stable `sortBy` for equal keys.
        let mut sorted = variants;
        sorted.sort_by(|a, b| order[&a.chrom].cmp(&order[&b.chrom]).then(a.pos.cmp(&b.pos)));
        Ok(sorted)
    }

    /// Checks if a file is a VCF file based on its extension.
    ///
    /// Determines whether a file path represents a VCF file by examining its extension.
    /// Recognizes both `.vcf` and `.vcf.gz` extensions (case-insensitive).
    ///
    /// # Arguments
    ///
    /// * `path` - The file path to check
    ///
    /// # Returns
    ///
    /// `true` if the file has a VCF extension, `false` otherwise.
    fn is_vcf_file(&self, path: &Path) -> bool {
        if let Some(ext) = path.extension() {
            let ext_str = ext.to_string_lossy().to_lowercase();
            ext_str == "vcf" || ext_str == "gz" && path.to_string_lossy().ends_with(".vcf.gz")
        } else {
            false
        }
    }

    /// Validates that the reference FASTA has required index and dictionary files.
    ///
    /// Checks for the presence of:
    /// - `.fai` index file (required for random access)
    /// - `.dict` sequence dictionary file (required for coordinate validation)
    ///
    /// # Returns
    ///
    /// `Ok(())` if both files exist and are readable.
    ///
    /// # Errors
    ///
    /// Returns an error with helpful message if either file is missing, including
    /// the command to generate the missing file.
    fn validate_reference_files(&self) -> Result<()> {
        // Check for FASTA index (.fai)
        let fai_path = self.reference.with_extension("fa.fai");
        let fai_path_alt = format!("{}.fai", self.reference.display());
        let has_fai = fai_path.exists() || std::path::Path::new(&fai_path_alt).exists();

        if !has_fai {
            bail!(
                "The reference file has not been indexed. Please run:\n  \
                samtools faidx {}",
                self.reference.display()
            );
        }

        // Check for sequence dictionary (.dict)
        // Supports both fgbio/HTSJDK convention (ref.dict) and GATK convention (ref.fa.dict)
        if find_dict_path(&self.reference).is_none() {
            bail!(
                "The reference file has no sequence dictionary. Tried:\n  \
                - {}\n  \
                - {}.dict\n\
                Please run: samtools dict {} -o {}",
                self.reference.with_extension("dict").display(),
                self.reference.display(),
                self.reference.display(),
                self.reference.with_extension("dict").display()
            );
        }

        Ok(())
    }

    /// Loads variant positions from a VCF file.
    ///
    /// Parses a VCF file to extract SNP variant positions. Filters variants based on MAF
    /// (minor allele frequency) if a sample is specified. Extracts genotype and filter
    /// information for each variant.
    ///
    /// # Returns
    ///
    /// A vector of `Variant` objects containing chromosome, position, reference base,
    /// and optional genotype/filter information.
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - VCF file cannot be opened or parsed
    /// - Reference FASTA cannot be read
    /// - Variant positions are invalid
    fn load_variants_from_vcf(&self) -> Result<Vec<Variant>> {
        use noodles::vcf;

        let mut reader = vcf::io::reader::Builder::default().build_from_path(&self.input)?;
        let header = reader.read_header()?;
        let mut variants = Vec::new();

        // Load reference for getting ref bases
        let mut fasta_reader = noodles::fasta::io::indexed_reader::Builder::default()
            .build_from_path(&self.reference)?;

        // Determine which sample to use for genotype extraction
        let sample_name = if let Some(ref name) = self.sample {
            Some(name.clone())
        } else if header.sample_names().len() == 1 {
            Some(header.sample_names()[0].clone())
        } else {
            None
        };

        for result in reader.record_bufs(&header) {
            let record = result?;

            // Get chromosome name
            let chrom = record.reference_sequence_name().to_string();

            // Get position
            let pos = usize::from(
                record
                    .variant_start()
                    .ok_or_else(|| anyhow::anyhow!("Missing variant position"))?,
            ) as i32;

            // Get reference base from FASTA
            let ref_base = self.get_reference_base(&mut fasta_reader, &chrom, pos)?;

            // Only process SNPs (single base variants)
            if record.reference_bases().len() != 1 {
                continue;
            }

            // Apply MAF filtering if sample specified
            if let Some(ref sname) = sample_name {
                if let Some(maf) = Self::calculate_maf(&record, &header, sname) {
                    if maf > self.maf {
                        continue; // Skip variants above MAF threshold
                    }
                }
            }

            // Get genotype if sample specified
            let genotype = if let Some(ref sname) = sample_name {
                Self::extract_genotype_for_sample(&record, &header, sname)
            } else {
                None
            };

            // Get filters
            let filters = Self::extract_filters(&record, &header);

            let mut variant = Variant::new(chrom, pos, ref_base);
            if let Some(gt) = genotype {
                variant = variant.with_genotype(gt);
            }
            variant = variant.with_filters(filters);

            variants.push(variant);
        }

        Ok(variants)
    }

    /// Loads variant positions from an interval list file.
    ///
    /// Parses an interval list file (BED-like format) to extract variant positions.
    /// Creates a variant for each position within each interval by reading the reference
    /// base from the FASTA file.
    ///
    /// # Returns
    ///
    /// A vector of `Variant` objects, one for each position in all intervals.
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - Interval file cannot be opened or parsed
    /// - Reference FASTA cannot be read
    /// - Interval coordinates are invalid
    fn load_variants_from_intervals(&self) -> Result<Vec<Variant>> {
        use std::fs::File;
        use std::io::{BufRead, BufReader};

        let file = File::open(&self.input)?;
        let reader = BufReader::new(file);
        let mut variants = Vec::new();

        // Load reference
        let mut fasta_reader = noodles::fasta::io::indexed_reader::Builder::default()
            .build_from_path(&self.reference)?;

        for line in reader.lines() {
            let line = line?;
            let line = line.trim();

            // Skip headers and comments
            if line.is_empty() || line.starts_with('@') || line.starts_with('#') {
                continue;
            }

            // Parse interval format: chr\tstart\tend
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 3 {
                continue;
            }

            let chrom = fields[0].to_string();
            let start: i32 = fields[1].parse()?;
            let end: i32 = fields[2].parse()?;

            // Create a variant for each position in the interval
            for pos in start..=end {
                let ref_base = self.get_reference_base(&mut fasta_reader, &chrom, pos)?;
                variants.push(Variant::new(chrom.clone(), pos, ref_base));
            }
        }

        Ok(variants)
    }

    /// Retrieves the reference base at a specific genomic position.
    ///
    /// Queries the reference FASTA file to extract the base at the given chromosome
    /// and position. Returns 'N' if the position is invalid or cannot be read.
    ///
    /// # Arguments
    ///
    /// * `reader` - Indexed FASTA reader for accessing the reference
    /// * `chrom` - Chromosome/contig name
    /// * `pos` - 1-based genomic position
    ///
    /// # Returns
    ///
    /// The uppercase reference base character, or 'N' if unavailable.
    ///
    /// # Errors
    ///
    /// Returns an error if the region query fails or position conversion fails.
    fn get_reference_base(
        &self,
        reader: &mut noodles::fasta::io::IndexedReader<
            noodles::fasta::io::BufReader<std::fs::File>,
        >,
        chrom: &str,
        pos: i32,
    ) -> Result<char> {
        use noodles::core::Position;

        let region = noodles::core::Region::new(
            chrom,
            Position::try_from(pos as usize)?..=Position::try_from(pos as usize)?,
        );

        let record = reader.query(&region)?;
        let sequence = record.sequence();

        if sequence.is_empty() {
            Ok('N')
        } else {
            // Get first base from sequence
            let base = sequence.as_ref().first().copied().unwrap_or(b'N');
            Ok((base as char).to_ascii_uppercase())
        }
    }

    /// Extract genotype for a specific sample.
    ///
    /// Returns the genotype rendered like htsjdk's `Genotype.getGenotypeString`
    /// (e.g. `A/T`, `A|T`) — the format fgbio writes to the `genotype` column —
    /// or `None` if the sample has no `GT` value.
    fn extract_genotype_for_sample(
        record: &noodles::vcf::variant::RecordBuf,
        header: &noodles::vcf::Header,
        sample_name: &str,
    ) -> Option<String> {
        use noodles::vcf::variant::record_buf::samples::sample::Value;

        // Find the sample index
        let sample_idx = header.sample_names().iter().position(|s| s == sample_name)?;

        let samples = record.samples();
        let gt_idx = samples.keys().as_ref().iter().position(|k| k == "GT")?;
        let sample = samples.get_index(sample_idx)?;

        match sample.values().get(gt_idx) {
            Some(Some(Value::Genotype(genotype))) => Some(format_genotype_string(
                genotype.as_ref(),
                record.reference_bases(),
                record.alternate_bases().as_ref(),
            )),
            _ => None,
        }
    }

    /// Extract filters from VCF record
    fn extract_filters(
        record: &noodles::vcf::variant::RecordBuf,
        header: &noodles::vcf::Header,
    ) -> String {
        use noodles::vcf::variant::record::Filters;

        let filters = record.filters();

        // Check if PASS (no filters)
        if filters.is_pass() {
            return "PASS".to_string();
        }

        // Collect all filter IDs and sort them
        let mut filter_names: Vec<String> = filters
            .iter(header)
            .filter_map(std::result::Result::ok)
            .map(std::string::ToString::to_string)
            .collect();

        if filter_names.is_empty() {
            "PASS".to_string()
        } else {
            filter_names.sort();
            filter_names.join(",")
        }
    }

    /// Calculates minor allele frequency from genotype info.
    ///
    /// Follows Scala logic exactly:
    /// 1. Try to extract AF (allele frequency) attribute
    /// 2. Fall back to calculating from AD (allelic depth) if AF not present
    /// 3. Return None if neither is available
    ///
    /// # Arguments
    ///
    /// * `record` - VCF record
    /// * `header` - VCF header
    /// * `sample_name` - Name of sample to extract MAF for
    ///
    /// # Returns
    ///
    /// MAF as f64 between 0 and 1, or None if cannot be calculated
    fn calculate_maf(
        record: &noodles::vcf::variant::RecordBuf,
        header: &noodles::vcf::Header,
        sample_name: &str,
    ) -> Option<f64> {
        use noodles::vcf::variant::record::samples::Sample;

        // Find the sample index
        let sample_idx = header.sample_names().iter().position(|s| s == sample_name)?;

        let samples = record.samples();
        let keys = samples.keys();

        // Try to get AF (allele frequency) field first - matches Scala's getExtendedAttribute("AF")
        if let Some(af_idx) = keys.as_ref().iter().position(|k| k == "AF") {
            if let Some(sample_values) = samples.get_index(sample_idx) {
                if let Some(Ok(Some(value))) = sample_values.get_index(header, af_idx) {
                    // Convert to string and parse - handles String, Float, or Array types
                    let value_str = format!("{value:?}");
                    // Handle array format like "[0.01]" or simple format like "0.01"
                    let cleaned = value_str.trim_matches(|c| c == '[' || c == ']' || c == '"');
                    // Take first value if comma-separated
                    if let Some(first_val) = cleaned.split(',').next() {
                        if let Ok(af) = first_val.trim().parse::<f64>() {
                            return Some(af);
                        }
                    }
                }
            }
        }

        // Fall back to calculating from AD (allele depth) - matches Scala's getAD
        if let Some(ad_idx) = keys.as_ref().iter().position(|k| k == "AD") {
            if let Some(sample_values) = samples.get_index(sample_idx) {
                if let Some(Ok(Some(value))) = sample_values.get_index(header, ad_idx) {
                    // AD should be an array like "[100,2]" or "100,2"
                    let value_str = format!("{value:?}");
                    let cleaned = value_str.trim_matches(|c| c == '[' || c == ']' || c == '"');
                    let depths: Vec<i32> =
                        cleaned.split(',').filter_map(|s| s.trim().parse::<i32>().ok()).collect();

                    if !depths.is_empty() {
                        let total: i32 = depths.iter().sum();
                        if total > 0 {
                            let ref_depth = depths[0];
                            // Scala: 1 - ad(0) / ad.sum.toDouble
                            return Some(1.0 - (ref_depth as f64 / total as f64));
                        }
                    }
                }
            }
        }

        None
    }

    /// Extracts consensus reads containing non-reference bases at variant positions.
    ///
    /// Queries the consensus BAM for reads overlapping each variant position. Filters reads
    /// to keep only those with non-reference bases at the variant site, extracting their
    /// molecular identifiers (MI tags) for downstream processing.
    ///
    /// # Arguments
    ///
    /// * `variants` - Slice of variant positions to check
    ///
    /// # Returns
    ///
    /// A `HashSet` of molecular identifier base strings (MI tags with strand suffix removed).
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - Consensus BAM cannot be opened or queried
    /// - Output BAM cannot be created
    /// - Record processing fails
    fn extract_consensus_reads(
        &self,
        variants: &[Variant],
        command_line: &str,
    ) -> Result<HashSet<String>> {
        use noodles::bam;

        let index = Self::read_bam_index(&self.consensus_bam)?;
        let mut reader = IndexedRawBamReader::from_path(&self.consensus_bam, index)?;
        let header = reader.read_header()?;

        // Add @PG record with PP chaining
        let header = crate::commands::common::add_pg_record(header, command_line)?;

        // Create output BAM
        let consensus_out_path = self.output.with_extension("consensus.bam");
        let mut writer = bam::io::writer::Builder.build_from_path(&consensus_out_path)?;
        writer.write_header(&header)?;

        let mut mi_set = HashSet::new();

        // Query each variant region
        for variant in variants {
            let start = noodles::core::Position::try_from(variant.pos as usize)?;
            let region = noodles::core::Region::new(variant.chrom.as_str(), start..=start);

            for result in reader.query(&header, &region)? {
                let record = result?;

                // Check if this read has a non-reference base at the variant position
                if self.has_non_reference_base(&record, variant)? {
                    // fgbio records `toMi(rec)` in its source set for every
                    // non-reference consensus read, erroring if the MI tag is
                    // absent (ReviewConsensusVariantsTest.scala:166).
                    let mi_base = to_mi(&record)?;
                    mi_set.insert(mi_base);

                    // Write raw record to output BAM.
                    write_raw_record(writer.get_mut(), &record)?;
                }
            }
        }

        Ok(mi_set)
    }

    /// Extracts grouped raw reads matching molecular identifiers from consensus reads.
    ///
    /// Reads through the grouped BAM file sequentially, extracting all reads whose MI tags
    /// match the provided set of molecular identifiers. These are the raw reads that
    /// contributed to the consensus reads with variants.
    ///
    /// # Arguments
    ///
    /// * `mi_set` - Set of molecular identifiers to match
    ///
    /// # Returns
    ///
    /// `Ok(())` on success.
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - Grouped BAM cannot be opened or read
    /// - Output BAM cannot be created or written
    fn extract_grouped_reads(&self, mi_set: &HashSet<String>, command_line: &str) -> Result<()> {
        use noodles::bam;

        // Use a plain (non-indexed) reader for the sequential scan; raw-byte mode avoids
        // the noodles decode/encode round-trip on every record.
        let mut bam_reader = bam::io::reader::Builder.build_from_path(&self.grouped_bam)?;
        let header = bam_reader.read_header()?;

        // Add @PG record with PP chaining (must happen before we consume the reader below)
        let header = crate::commands::common::add_pg_record(header, command_line)?;

        let grouped_out_path = self.output.with_extension("grouped.bam");
        let mut writer = bam::io::writer::Builder.build_from_path(&grouped_out_path)?;
        writer.write_header(&header)?;

        // Extract the inner BGZF reader and wrap in RawBamReader for zero-copy record iteration.
        let mut raw_reader = RawBamReader::new(bam_reader.into_inner());
        let mut raw_rec = RawRecord::new();

        // Read through entire grouped BAM and extract matching reads
        let progress = ProgressTracker::new("Processed grouped reads").with_interval(1_000_000);
        loop {
            let bytes_read = raw_reader.read_record(&mut raw_rec)?;
            if bytes_read == 0 {
                break; // EOF
            }
            progress.log_if_needed(1);

            // Extract MI tag directly from raw bytes; no RecordBuf decode needed.
            // Unlike fgbio's `query(variants)` (which only touches variant-overlapping
            // reads and errors via `toMi` on a missing MI), this is a full sequential
            // scan of the grouped BAM, so a missing MI is skipped rather than an error:
            // an MI-less read cannot be in `mi_set` and erroring here would be stricter
            // than fgbio for reads that don't overlap any variant. Variant-overlapping
            // grouped reads are still MI-guarded in `generate_review_file` via `to_mi`.
            if let Some(mi_bytes) = find_string_tag_in_record(&raw_rec, SamTag::MI) {
                let mi = std::str::from_utf8(mi_bytes)?;
                let mi_base = extract_mi_base(mi);

                if mi_set.contains(mi_base) {
                    write_raw_record(writer.get_mut(), &raw_rec)?;
                }
            }
        }

        progress.log_final();
        Ok(())
    }

    /// Generates the detailed review TSV file with per-variant and per-read information.
    ///
    /// Creates a comprehensive review file containing:
    /// - Variant-level information (position, genotype, filters)
    /// - Consensus-level base counts across all reads
    /// - Per-consensus-read information (base call, quality, insert string)
    /// - Raw read base counts for each consensus read
    ///
    /// # Arguments
    ///
    /// * `variants` - Slice of variants to include in the review
    ///
    /// # Returns
    ///
    /// `Ok(())` on success.
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - BAM files cannot be queried
    /// - Review TSV file cannot be created or written
    /// - Record processing fails
    fn generate_review_file(&self, variants: &[Variant]) -> Result<()> {
        let review_path = self.output.with_extension("txt");

        // Open both BAMs with indexed readers that yield RawRecord directly.
        let con_index = Self::read_bam_index(&self.consensus_bam)?;
        let mut consensus_reader = IndexedRawBamReader::from_path(&self.consensus_bam, con_index)?;
        let consensus_header = consensus_reader.read_header()?;

        let grp_index = Self::read_bam_index(&self.grouped_bam)?;
        let mut grouped_reader = IndexedRawBamReader::from_path(&self.grouped_bam, grp_index)?;
        let grouped_header = grouped_reader.read_header()?;

        let mut all_metrics = Vec::new();

        // Process each variant
        for variant in variants {
            // Query consensus BAM for this position
            let start = noodles::core::Position::try_from(variant.pos as usize)?;
            let region = noodles::core::Region::new(variant.chrom.as_str(), start..=start);

            // Collect all consensus reads at this position
            let mut consensus_reads: Vec<RawRecord> = Vec::new();
            for result in consensus_reader.query(&consensus_header, &region)? {
                consensus_reads.push(result?);
            }

            // Build consensus-level base counts
            let mut consensus_counts = BaseCounts::default();
            for record in &consensus_reads {
                if let Some(base) = self.get_base_at_position(record, variant.pos)? {
                    consensus_counts
                        .add_base(Self::normalize_base_for_variant(base, variant.ref_base));
                }
            }

            // Collect this variant's rows so they can be sorted by MI + read-number
            // before writing, matching fgbio (ReviewConsensusVariants.scala:258).
            let mut variant_rows: Vec<(String, ConsensusVariantReviewInfo)> = Vec::new();

            // Process each consensus read
            for record in consensus_reads {
                // Get the base at the variant position
                // If None, check if it's a spanning deletion
                let read_base = match self.get_base_at_position(&record, variant.pos)? {
                    Some(b) => Self::normalize_base_for_variant(b, variant.ref_base) as char,
                    // fgbio iterates `SamLocusIterator.getRecordAndOffsets`, which
                    // excludes reads with a deletion (or no base) at the locus, so
                    // spanning-deletion reads produce no detail row
                    // (ReviewConsensusVariantsTest.scala:250-251). Such reads are
                    // still extracted to the output BAMs by `has_non_reference_base`.
                    None => continue,
                };

                // Skip if it's the reference base (unless we're not ignoring Ns)
                if read_base == variant.ref_base {
                    continue;
                }

                // Skip N bases if requested
                if self.ignore_ns && read_base == 'N' {
                    continue;
                }

                // Get the quality at this position
                let read_qual = self.get_quality_at_position(&record, variant.pos)?.unwrap_or(0);

                // Source molecule id for this consensus read (fgbio calls
                // `toMi(rec)`, which errors when the MI tag is absent).
                let mi_base = to_mi(&record)?;

                // Format consensus read name with /1 or /2 suffix
                let read_name = String::from_utf8_lossy(record.read_name()).to_string();
                let is_first = record.is_first_segment();
                let consensus_read_name = format!("{}{}", read_name, read_number_suffix(is_first));

                // fgbio sorts each variant's rows by the string `toMi(rec) +
                // readNumberSuffix(rec)` (ReviewConsensusVariants.scala:258), so e.g.
                // "10/1" sorts before "2/1". This key mirrors that (MI base, not the
                // full read name, then the /1 or /2 suffix).
                let sort_key = format!("{}{}", mi_base, read_number_suffix(is_first));

                // Generate insert string
                let insert_string = format_insert_string(&record, &consensus_header);

                // Query grouped BAM for raw reads with matching MI
                let raw_counts = {
                    use std::collections::HashSet;

                    let start = noodles::core::Position::try_from(variant.pos as usize)?;
                    let region = noodles::core::Region::new(variant.chrom.as_str(), start..=start);

                    let mut counts = BaseCounts::default();
                    let mut seen_reads = HashSet::new();
                    let expected_read_num =
                        if consensus_read_name.ends_with("/1") { "/1" } else { "/2" };

                    for result in grouped_reader.query(&grouped_header, &region)? {
                        let rec = result?;

                        // Source molecule id for this grouped read. fgbio buckets
                        // the locus's grouped reads with `toMi(rec)` (erroring on a
                        // missing MI tag) before matching them to the consensus MI.
                        let read_mi_base = to_mi(&rec)?;
                        if read_mi_base != mi_base {
                            continue;
                        }

                        // Check read number matches
                        let is_first = rec.is_first_segment();
                        let read_num = read_number_suffix(is_first);
                        if read_num != expected_read_num {
                            continue;
                        }

                        // Ensure we only count each read once
                        let rn = String::from_utf8_lossy(rec.read_name()).to_string();
                        if !seen_reads.insert(rn) {
                            continue;
                        }

                        // Get base at position
                        if let Some(base) = self.get_base_at_position(&rec, variant.pos)? {
                            counts
                                .add_base(Self::normalize_base_for_variant(base, variant.ref_base));
                        }
                    }

                    counts
                };

                // Create the review info record
                let info = ConsensusVariantReviewInfo {
                    chrom: variant.chrom.clone(),
                    pos: variant.pos,
                    ref_allele: variant.ref_base.to_string(),
                    genotype: variant.genotype.clone().unwrap_or_else(|| "NA".to_string()),
                    filters: variant.filters.clone().unwrap_or_else(|| "NA".to_string()),
                    consensus_a: consensus_counts.a,
                    consensus_c: consensus_counts.c,
                    consensus_g: consensus_counts.g,
                    consensus_t: consensus_counts.t,
                    consensus_n: consensus_counts.n,
                    consensus_read: consensus_read_name,
                    consensus_insert: insert_string,
                    consensus_call: read_base,
                    consensus_qual: read_qual,
                    a: raw_counts.a,
                    c: raw_counts.c,
                    g: raw_counts.g,
                    t: raw_counts.t,
                    n: raw_counts.n,
                };

                variant_rows.push((sort_key, info));
            }

            // Emit this variant's rows in MI + read-number order; a stable sort keeps
            // BAM-query order for ties, matching fgbio's stable `sortBy`.
            variant_rows.sort_by(|a, b| a.0.cmp(&b.0));
            all_metrics.extend(variant_rows.into_iter().map(|(_, info)| info));
        }

        // Write TSV
        let mut writer = csv::WriterBuilder::new().delimiter(b'\t').from_path(&review_path)?;

        for metric in all_metrics {
            writer.serialize(&metric)?;
        }

        writer.flush()?;
        info!("Wrote review metrics to {}", review_path.display());
        Ok(())
    }

    /// Get the ASCII base at a specific reference position in a read.
    ///
    /// Returns the ASCII byte (e.g. b'A', b'C', b'N') at `ref_pos`, or `None` if the
    /// position is not covered by this read.
    fn get_base_at_position(&self, record: &RawRecord, ref_pos: i32) -> Result<Option<u8>> {
        if let Some(offset) = self.calculate_read_offset(record, ref_pos)? {
            let l_seq = record.l_seq() as usize;
            if offset < l_seq {
                // record.get_base returns the 4-bit BAM code (0–15); convert to ASCII.
                let code = record.get_base(offset) as usize;
                return Ok(Some(BAM_BASE_TO_ASCII[code]));
            }
        }
        Ok(None)
    }

    /// Map BAM's "same as reference" `=` sentinel to the actual reference base and
    /// uppercase the result so downstream comparisons against `variant.ref_base` and
    /// `BaseCounts::add_base` behave correctly.
    fn normalize_base_for_variant(base: u8, ref_base: char) -> u8 {
        if base == b'=' { ref_base.to_ascii_uppercase() as u8 } else { base.to_ascii_uppercase() }
    }

    /// Get the quality score at a specific reference position in a read.
    fn get_quality_at_position(&self, record: &RawRecord, ref_pos: i32) -> Result<Option<u8>> {
        if let Some(offset) = self.calculate_read_offset(record, ref_pos)? {
            let qual = record.quality_scores();
            if offset < qual.len() {
                return Ok(Some(qual[offset]));
            }
        }
        Ok(None)
    }

    /// Check if a record has a non-reference base at the variant position.
    fn has_non_reference_base(&self, record: &RawRecord, variant: &Variant) -> Result<bool> {
        // Skip unmapped reads
        if record.is_unmapped() {
            return Ok(false);
        }

        // Get alignment start / end (1-based)
        let start = match record.alignment_start_1based() {
            Some(pos) => pos as i32,
            None => return Ok(false),
        };
        let end = match record.alignment_end_1based() {
            Some(pos) => pos as i32,
            None => return Ok(false),
        };

        // Check if variant position overlaps this read
        if variant.pos < start || variant.pos > end {
            return Ok(false);
        }

        // Get the base at the variant position using CIGAR-aware mapping
        let read_offset = self.calculate_read_offset(record, variant.pos)?;

        if let Some(offset) = read_offset {
            let l_seq = record.l_seq() as usize;
            if offset < l_seq {
                // record.get_base returns the 4-bit BAM code; convert to ASCII then to char.
                let code = record.get_base(offset) as usize;
                let base_char =
                    Self::normalize_base_for_variant(BAM_BASE_TO_ASCII[code], variant.ref_base)
                        as char;

                // Check if base differs from reference
                if base_char != variant.ref_base {
                    // Also check if we should ignore Ns
                    if self.ignore_ns && base_char == 'N' {
                        return Ok(false);
                    }
                    return Ok(true);
                }
            }
        } else {
            // read_offset is None — check if it's a spanning deletion
            if self.is_spanning_deletion(record, variant.pos)? {
                return Ok(true); // Spanning deletions are considered non-reference
            }
        }

        Ok(false)
    }

    /// Check if a reference position falls within a deletion in the read.
    fn is_spanning_deletion(&self, record: &RawRecord, ref_pos: i32) -> Result<bool> {
        let start = match record.alignment_start_1based() {
            Some(pos) => pos as i32,
            None => return Ok(false),
        };
        let end = match record.alignment_end_1based() {
            Some(pos) => pos as i32,
            None => return Ok(false),
        };

        if ref_pos < start || ref_pos > end {
            return Ok(false); // Position not covered by read at all
        }

        let mut current_ref_pos = start;

        for op in record.cigar_ops_typed() {
            let len = op.len() as i32;
            match op.kind() {
                CigarKind::Match | CigarKind::SequenceMatch | CigarKind::SequenceMismatch => {
                    current_ref_pos += len;
                }
                CigarKind::Deletion => {
                    // Check if target position is within this deletion
                    if ref_pos >= current_ref_pos && ref_pos < current_ref_pos + len {
                        return Ok(true);
                    }
                    current_ref_pos += len;
                }
                CigarKind::Skip => {
                    // N (skipped region) consumes reference but is not a deletion allele.
                    current_ref_pos += len;
                }
                CigarKind::Insertion
                | CigarKind::SoftClip
                | CigarKind::HardClip
                | CigarKind::Pad => {
                    // These don't affect reference positions
                }
            }
        }

        Ok(false)
    }

    /// Calculate the read (query) offset for a reference position using CIGAR.
    ///
    /// Returns `None` if the reference position falls within a deletion or
    /// is not covered by the read.
    fn calculate_read_offset(&self, record: &RawRecord, ref_pos: i32) -> Result<Option<usize>> {
        let start = match record.alignment_start_1based() {
            Some(pos) => pos as i32,
            None => return Ok(None),
        };

        let mut current_ref_pos = start;
        let mut current_read_pos: usize = 0;

        for op in record.cigar_ops_typed() {
            let len = op.len() as usize;
            match op.kind() {
                CigarKind::Match | CigarKind::SequenceMatch | CigarKind::SequenceMismatch => {
                    if current_ref_pos + len as i32 > ref_pos {
                        // Target position is in this op
                        let offset_in_op = (ref_pos - current_ref_pos) as usize;
                        return Ok(Some(current_read_pos + offset_in_op));
                    }
                    current_ref_pos += len as i32;
                    current_read_pos += len;
                }
                CigarKind::Insertion | CigarKind::SoftClip => {
                    current_read_pos += len;
                }
                CigarKind::Deletion | CigarKind::Skip => {
                    if current_ref_pos + len as i32 > ref_pos {
                        // Position is deleted in this read
                        return Ok(None);
                    }
                    current_ref_pos += len as i32;
                }
                CigarKind::HardClip | CigarKind::Pad => {
                    // These don't affect positions
                }
            }
        }

        Ok(None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
    use std::path::PathBuf;
    use tempfile::TempDir;

    // Test utilities for creating synthetic test data
    mod test_utils {
        use super::*;
        use fgumi_raw_bam::{
            SamBuilder as RawSamBuilder, flags, raw_record_to_record_buf, testutil::encode_op,
        };
        use noodles::sam::Header;
        use noodles::sam::alignment::RecordBuf;
        use noodles::sam::header::record::value::{Map, map::ReferenceSequence};
        use std::num::NonZeroUsize;

        /// Create a simple reference FASTA with two chromosomes
        pub fn create_test_reference(dir: &TempDir) -> PathBuf {
            let ref_path = dir.path().join("ref.fa");
            let fai_path = dir.path().join("ref.fa.fai");
            let dict_path = dir.path().join("ref.dict");

            // Create FASTA with 100bp of A's on chr1, 100bp of C's on chr2
            let fasta_content = b">chr1\n\
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n\
>chr2\n\
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n";

            std::fs::write(&ref_path, fasta_content).expect("failed to write file");

            // Create FAI
            let fai_content = b"chr1\t100\t6\t100\t101\nchr2\t100\t114\t100\t101\n";
            std::fs::write(&fai_path, fai_content).expect("failed to write file");

            // Create dictionary
            let dict_content = b"@HD\tVN:1.5\tSO:unsorted\n\
@SQ\tSN:chr1\tLN:100\n\
@SQ\tSN:chr2\tLN:100\n";
            std::fs::write(&dict_path, dict_content).expect("failed to write file");

            ref_path
        }

        /// Create a SAM header with sequence dictionary
        pub fn create_test_header() -> Header {
            use bstr::BString;

            // Create header with coordinate sort order
            let mut header = Header::builder()
                .add_reference_sequence(
                    BString::from("chr1"),
                    Map::<ReferenceSequence>::new(
                        NonZeroUsize::try_from(100).expect("non-zero value required"),
                    ),
                )
                .add_reference_sequence(
                    BString::from("chr2"),
                    Map::<ReferenceSequence>::new(
                        NonZeroUsize::try_from(100).expect("non-zero value required"),
                    ),
                )
                .build();

            // Parse and set the sort order in the header line
            let header_str = "@HD\tVN:1.6\tSO:coordinate\n";
            let parsed_header: Header = header_str.parse().expect("valid parse input");
            if let Some(hd) = parsed_header.header() {
                header = Header::builder()
                    .set_header(hd.clone())
                    .add_reference_sequence(
                        BString::from("chr1"),
                        Map::<ReferenceSequence>::new(
                            NonZeroUsize::try_from(100).expect("non-zero value required"),
                        ),
                    )
                    .add_reference_sequence(
                        BString::from("chr2"),
                        Map::<ReferenceSequence>::new(
                            NonZeroUsize::try_from(100).expect("non-zero value required"),
                        ),
                    )
                    .build();
            }

            header
        }

        fn to_record_buf(raw: fgumi_raw_bam::RawRecord) -> RecordBuf {
            raw_record_to_record_buf(&raw, &Header::default())
                .expect("raw_record_to_record_buf failed in test")
        }

        fn parse_cigar_to_ops(cigar_str: &str) -> Vec<u32> {
            let mut ops = Vec::new();
            let mut len_buf = String::new();
            for ch in cigar_str.chars() {
                if ch.is_ascii_digit() {
                    len_buf.push(ch);
                } else {
                    let len: usize = len_buf.parse().expect("valid cigar length");
                    len_buf.clear();
                    let op_code: u32 = match ch {
                        'M' => 0,
                        'I' => 1,
                        'D' => 2,
                        'N' => 3,
                        'S' => 4,
                        'H' => 5,
                        'P' => 6,
                        '=' => 7,
                        'X' => 8,
                        other => panic!("unknown CIGAR op: {other}"),
                    };
                    ops.push(encode_op(op_code, len));
                }
            }
            ops
        }

        /// Create a simple read pair
        #[allow(clippy::too_many_arguments, clippy::cast_possible_truncation)]
        pub fn create_read_pair(
            name: &str,
            chrom_idx: usize,
            start1: i32,
            start2: i32,
            bases1: &[u8],
            bases2: &[u8],
            mi: &str,
            cigar1: Option<&str>,
        ) -> (RecordBuf, RecordBuf) {
            let cigar1_str = cigar1.map_or_else(|| format!("{}M", bases1.len()), String::from);
            let cigar1_ops = parse_cigar_to_ops(&cigar1_str);
            let cigar2_ops = [encode_op(0, bases2.len())];

            let mut b1 = RawSamBuilder::new();
            b1.read_name(name.as_bytes())
                .flags(flags::PAIRED | flags::FIRST_SEGMENT)
                .ref_id(chrom_idx as i32)
                .pos(start1 - 1)
                .mapq(60)
                .cigar_ops(&cigar1_ops)
                .sequence(bases1)
                .qualities(&vec![45u8; bases1.len()])
                .mate_ref_id(chrom_idx as i32)
                .mate_pos(start2 - 1);
            b1.add_string_tag(SamTag::MI, mi.as_bytes());
            let r1 = to_record_buf(b1.build());

            let mut b2 = RawSamBuilder::new();
            b2.read_name(name.as_bytes())
                .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE)
                .ref_id(chrom_idx as i32)
                .pos(start2 - 1)
                .mapq(60)
                .cigar_ops(&cigar2_ops)
                .sequence(bases2)
                .qualities(&vec![45u8; bases2.len()])
                .mate_ref_id(chrom_idx as i32)
                .mate_pos(start1 - 1);
            b2.add_string_tag(SamTag::MI, mi.as_bytes());
            let r2 = to_record_buf(b2.build());

            (r1, r2)
        }

        /// Create a test BAM file with synthetic reads
        pub fn create_test_bams(dir: &TempDir) -> (PathBuf, PathBuf) {
            use noodles::bam;
            use noodles::sam::alignment::io::Write;

            let raw_path = dir.path().join("raw.bam");
            let consensus_path = dir.path().join("consensus.bam");

            let header = create_test_header();

            // Create raw BAM
            let mut raw_writer = bam::io::Writer::new(
                std::fs::File::create(&raw_path).expect("failed to create file"),
            );
            raw_writer.write_header(&header).expect("failed to write BAM header");

            // Create consensus BAM
            let mut con_writer = bam::io::Writer::new(
                std::fs::File::create(&consensus_path).expect("failed to create file"),
            );
            con_writer.write_header(&header).expect("failed to write BAM header");

            // Add reads with variant at chr1:10 (A->T)
            let (r1, r2) =
                create_read_pair("A1", 0, 6, 50, b"AAAATAAAAA", b"AAAAAAAAAA", "A", None);
            raw_writer.write_alignment_record(&header, &r1).expect("failed to write BAM record");
            raw_writer.write_alignment_record(&header, &r2).expect("failed to write BAM record");

            let (r1, r2) =
                create_read_pair("A2", 0, 6, 50, b"AAAATAAAAG", b"AAAAAAAAAA", "A", None);
            raw_writer.write_alignment_record(&header, &r1).expect("failed to write BAM record");
            raw_writer.write_alignment_record(&header, &r2).expect("failed to write BAM record");

            let (r1, r2) = create_read_pair("A", 0, 6, 50, b"AAAATAAAAN", b"AAAAAAAAAA", "A", None);
            con_writer.write_alignment_record(&header, &r1).expect("failed to write BAM record");
            con_writer.write_alignment_record(&header, &r2).expect("failed to write BAM record");

            // Add reads at chr1:20 (A->C)
            let (r1, r2) =
                create_read_pair("B1", 0, 16, 50, b"AAAACAAAAA", b"AAAAAAAAAA", "B", None);
            raw_writer.write_alignment_record(&header, &r1).expect("failed to write BAM record");
            raw_writer.write_alignment_record(&header, &r2).expect("failed to write BAM record");

            let (r1, r2) =
                create_read_pair("B", 0, 16, 50, b"AAAACAAAAA", b"AAAAAAAAAA", "B", None);
            con_writer.write_alignment_record(&header, &r1).expect("failed to write BAM record");
            con_writer.write_alignment_record(&header, &r2).expect("failed to write BAM record");

            // Reference read at chr1:20
            let (r1, r2) =
                create_read_pair("C1", 0, 17, 50, b"AAAAAAAAAA", b"AAAAAAAAAA", "C", None);
            raw_writer.write_alignment_record(&header, &r1).expect("failed to write BAM record");
            raw_writer.write_alignment_record(&header, &r2).expect("failed to write BAM record");

            let (r1, r2) =
                create_read_pair("C", 0, 17, 50, b"AAAAAAAAAA", b"AAAAAAAAAA", "C", None);
            con_writer.write_alignment_record(&header, &r1).expect("failed to write BAM record");
            con_writer.write_alignment_record(&header, &r2).expect("failed to write BAM record");

            // Reads at chr1:30 with various edge cases
            // D: spanning deletion
            let (r1, r2) = create_read_pair(
                "D1",
                0,
                25,
                60,
                b"AAAAAAAAAA",
                b"AAAAAAAAAA",
                "D",
                Some("4M4D6M"),
            );
            raw_writer.write_alignment_record(&header, &r1).expect("failed to write BAM record");
            raw_writer.write_alignment_record(&header, &r2).expect("failed to write BAM record");

            let (r1, r2) =
                create_read_pair("D", 0, 25, 60, b"AAAAAAAAAA", b"AAAAAAAAAA", "D", Some("4M4D6M"));
            con_writer.write_alignment_record(&header, &r1).expect("failed to write BAM record");
            con_writer.write_alignment_record(&header, &r2).expect("failed to write BAM record");

            // E: no-call (N)
            let (r1, r2) =
                create_read_pair("E1", 0, 26, 60, b"AAAANAAAAA", b"AAAAAAAAAA", "E", None);
            raw_writer.write_alignment_record(&header, &r1).expect("failed to write BAM record");
            raw_writer.write_alignment_record(&header, &r2).expect("failed to write BAM record");

            let (r1, r2) =
                create_read_pair("E", 0, 26, 60, b"AAAANAAAAA", b"AAAAAAAAAA", "E", None);
            con_writer.write_alignment_record(&header, &r1).expect("failed to write BAM record");
            con_writer.write_alignment_record(&header, &r2).expect("failed to write BAM record");

            // F: variant allele (G)
            let (r1, r2) =
                create_read_pair("F1", 0, 27, 60, b"AAAGAAAAAA", b"AAAAAAAAAA", "F", None);
            raw_writer.write_alignment_record(&header, &r1).expect("failed to write BAM record");
            raw_writer.write_alignment_record(&header, &r2).expect("failed to write BAM record");

            let (r1, r2) =
                create_read_pair("F", 0, 27, 60, b"AAAGAAAAAA", b"AAAAAAAAAA", "F", None);
            con_writer.write_alignment_record(&header, &r1).expect("failed to write BAM record");
            con_writer.write_alignment_record(&header, &r2).expect("failed to write BAM record");

            // Reads at chr2:20 where both ends overlap variant
            let (r1, r2) =
                create_read_pair("H1", 1, 15, 19, b"CCCCCTCCCC", b"CTCCCCCCCC", "H", None);
            raw_writer.write_alignment_record(&header, &r1).expect("failed to write BAM record");
            raw_writer.write_alignment_record(&header, &r2).expect("failed to write BAM record");

            let (r1, r2) =
                create_read_pair("H", 1, 15, 19, b"CCCCCTCCCC", b"CTCCCCCCCC", "H", None);
            con_writer.write_alignment_record(&header, &r1).expect("failed to write BAM record");
            con_writer.write_alignment_record(&header, &r2).expect("failed to write BAM record");

            // Close writers to flush data
            drop(raw_writer);
            drop(con_writer);

            // Create BAM index files using noodles bam::fs::index
            use noodles::bam::bai;

            // Index raw BAM
            let raw_index_path = raw_path.with_extension("bam.bai");
            let raw_index = bam::fs::index(&raw_path).expect("failed to index BAM file");
            let mut raw_index_writer = bai::io::Writer::new(
                std::fs::File::create(&raw_index_path).expect("failed to create file"),
            );
            raw_index_writer.write_index(&raw_index).expect("failed to write BAM index");

            // Index consensus BAM
            let consensus_index_path = consensus_path.with_extension("bam.bai");
            let con_index = bam::fs::index(&consensus_path).expect("failed to index BAM file");
            let mut con_index_writer = bai::io::Writer::new(
                std::fs::File::create(&consensus_index_path).expect("failed to create file"),
            );
            con_index_writer.write_index(&con_index).expect("failed to write BAM index");

            (raw_path, consensus_path)
        }

        /// Create a test VCF with variants
        pub fn create_test_vcf(dir: &TempDir) -> PathBuf {
            let vcf_path = dir.path().join("variants.vcf");
            let vcf_content = b"##fileformat=VCFv4.2\n\
##contig=<ID=chr1,length=100>\n\
##contig=<ID=chr2,length=100>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">\n\
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele Depth\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttumor\n\
chr1\t10\t.\tA\tT\t.\tPASS\t.\tGT:AF\t0/1:0.01\n\
chr1\t20\t.\tA\tC\t.\tPASS\t.\tGT:AF\t0/1:0.01\n\
chr1\t30\t.\tA\tG\t.\tPASS\t.\tGT:AF\t0/1:0.01\n\
chr2\t20\t.\tC\tT\t.\tPASS\t.\tGT:AD\t0/1:100,2\n";
            std::fs::write(&vcf_path, vcf_content).expect("failed to write file");
            vcf_path
        }

        /// Create an empty VCF
        pub fn create_empty_vcf(dir: &TempDir) -> PathBuf {
            let vcf_path = dir.path().join("empty.vcf");
            let vcf_content = b"##fileformat=VCFv4.2\n\
##contig=<ID=chr1,length=100>\n\
##contig=<ID=chr2,length=100>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttumor\n";
            std::fs::write(&vcf_path, vcf_content).expect("failed to write file");
            vcf_path
        }

        /// Create an empty interval list
        pub fn create_empty_interval_list(dir: &TempDir) -> PathBuf {
            let interval_path = dir.path().join("empty.interval_list");
            let interval_content = b"@HD\tVN:1.5\tSO:unsorted\n\
@SQ\tSN:chr1\tLN:100\n\
@SQ\tSN:chr2\tLN:100\n";
            std::fs::write(&interval_path, interval_content).expect("failed to write file");
            interval_path
        }

        /// One read-pair spec for [`create_ordering_bams`]:
        /// `(name, chrom_idx, r1_start, r1_variant_base, mi, r2_override)`, where
        /// `r2_override` is `Some((r2_start, r2_base))` to place R2 over a queried locus.
        type OrderingRead<'a> = (&'a str, usize, i32, u8, &'a str, Option<(i32, u8)>);

        /// Build consensus + grouped BAMs for the output-row-ordering test (issue #497).
        ///
        /// Each tuple is `(name, chrom_idx, start, variant_base, mi, r2_override)`. The
        /// read pair's R1 is forward at `start` and carries `variant_base` at offset 4
        /// (reference position `start + 4`). By default R2 is reverse at position 50 with
        /// all-reference bases, off every variant, so it yields no review row. When
        /// `r2_override` is `Some((r2_start, r2_base))`, R2 is instead placed at `r2_start`
        /// carrying `r2_base` at offset 4, so it overlaps the queried locus and (when the
        /// base is non-reference) emits a `/2` row — letting a single MI produce both a
        /// `/1` and a `/2` row to exercise the read-number tie-break. Read pairs are
        /// written in the given slice order, so the natural BAM/query order can differ from
        /// fgbio's MI+read-number sorted order. The same pairs are written to the grouped
        /// BAM so the raw-count lookups find a matching MI.
        ///
        /// Returns `(grouped_bam, consensus_bam)`.
        pub fn create_ordering_bams(dir: &TempDir, reads: &[OrderingRead]) -> (PathBuf, PathBuf) {
            use noodles::bam;
            use noodles::bam::bai;
            use noodles::sam::alignment::io::Write;

            let grouped_path = dir.path().join("ordering_grouped.bam");
            let consensus_path = dir.path().join("ordering_consensus.bam");
            let header = create_test_header();

            let mut grp_writer = bam::io::Writer::new(
                std::fs::File::create(&grouped_path).expect("failed to create grouped BAM"),
            );
            grp_writer.write_header(&header).expect("failed to write grouped header");
            let mut con_writer = bam::io::Writer::new(
                std::fs::File::create(&consensus_path).expect("failed to create consensus BAM"),
            );
            con_writer.write_header(&header).expect("failed to write consensus header");

            for (name, chrom_idx, start, variant_base, mi, r2_override) in reads {
                let mut bases1 = vec![b'A'; 10];
                bases1[4] = *variant_base;
                // R2 defaults to position 50 with all-reference bases (off every queried
                // variant); when overridden it moves to `r2_start` and carries `r2_base` at
                // offset 4 so it overlaps the locus and can emit a `/2` row.
                let mut bases2 = vec![b'A'; 10];
                let start2 = match r2_override {
                    Some((r2_start, r2_base)) => {
                        bases2[4] = *r2_base;
                        *r2_start
                    }
                    None => 50,
                };
                let (r1, r2) =
                    create_read_pair(name, *chrom_idx, *start, start2, &bases1, &bases2, mi, None);
                for rec in [&r1, &r2] {
                    con_writer
                        .write_alignment_record(&header, rec)
                        .expect("failed to write consensus record");
                    grp_writer
                        .write_alignment_record(&header, rec)
                        .expect("failed to write grouped record");
                }
            }

            drop(grp_writer);
            drop(con_writer);

            for path in [&consensus_path, &grouped_path] {
                let index = bam::fs::index(path).expect("failed to index BAM");
                let mut writer = bai::io::Writer::new(
                    std::fs::File::create(path.with_extension("bam.bai"))
                        .expect("failed to create BAM index"),
                );
                writer.write_index(&index).expect("failed to write BAM index");
            }

            (grouped_path, consensus_path)
        }

        /// VCF for the ordering test: the chr1:10 variant is listed *before* chr1:5, so
        /// the on-disk order differs from coordinate order (issue #497). No sample column,
        /// so callers pass `sample: None`.
        pub fn create_ordering_vcf(dir: &TempDir) -> PathBuf {
            let vcf_path = dir.path().join("ordering.vcf");
            let vcf_content = b"##fileformat=VCFv4.2\n\
##contig=<ID=chr1,length=100>\n\
##contig=<ID=chr2,length=100>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t10\t.\tA\tT\t.\tPASS\t.\n\
chr1\t5\t.\tA\tC\t.\tPASS\t.\n";
            std::fs::write(&vcf_path, vcf_content).expect("failed to write file");
            vcf_path
        }
    }

    #[test]
    fn test_default_review_parameters() {
        let review = Review {
            input: PathBuf::from("variants.vcf"),
            consensus_bam: PathBuf::from("consensus.bam"),
            grouped_bam: PathBuf::from("grouped.bam"),
            reference: PathBuf::from("ref.fa"),
            output: PathBuf::from("output"),
            sample: None,
            ignore_ns: false,
            maf: 0.05,
        };

        assert!((review.maf - 0.05).abs() < f64::EPSILON);
        assert!(!review.ignore_ns);
        assert!(review.sample.is_none());
    }

    #[test]
    fn test_review_with_sample_name() {
        let review = Review {
            input: PathBuf::from("variants.vcf"),
            consensus_bam: PathBuf::from("consensus.bam"),
            grouped_bam: PathBuf::from("grouped.bam"),
            reference: PathBuf::from("ref.fa"),
            output: PathBuf::from("output"),
            sample: Some("SAMPLE1".to_string()),
            ignore_ns: false,
            maf: 0.05,
        };

        assert_eq!(review.sample, Some("SAMPLE1".to_string()));
    }

    #[test]
    fn test_review_ignore_ns_enabled() {
        let review = Review {
            input: PathBuf::from("variants.vcf"),
            consensus_bam: PathBuf::from("consensus.bam"),
            grouped_bam: PathBuf::from("grouped.bam"),
            reference: PathBuf::from("ref.fa"),
            output: PathBuf::from("output"),
            sample: None,
            ignore_ns: true,
            maf: 0.05,
        };

        assert!(review.ignore_ns);
    }

    #[test]
    fn test_review_custom_maf_threshold() {
        let review = Review {
            input: PathBuf::from("variants.vcf"),
            consensus_bam: PathBuf::from("consensus.bam"),
            grouped_bam: PathBuf::from("grouped.bam"),
            reference: PathBuf::from("ref.fa"),
            output: PathBuf::from("output"),
            sample: None,
            ignore_ns: false,
            maf: 0.01,
        };

        assert!((review.maf - 0.01).abs() < f64::EPSILON);
    }

    #[test]
    fn test_is_vcf_file_detection() {
        let review = Review {
            input: PathBuf::from("test.vcf"),
            consensus_bam: PathBuf::from("consensus.bam"),
            grouped_bam: PathBuf::from("grouped.bam"),
            reference: PathBuf::from("ref.fa"),
            output: PathBuf::from("output"),
            sample: None,
            ignore_ns: false,
            maf: 0.05,
        };

        assert!(review.is_vcf_file(&PathBuf::from("variants.vcf")));
        assert!(review.is_vcf_file(&PathBuf::from("variants.vcf.gz")));
        assert!(!review.is_vcf_file(&PathBuf::from("intervals.bed")));
        assert!(!review.is_vcf_file(&PathBuf::from("file.txt")));
    }

    #[test]
    fn test_review_output_path_extensions() {
        let review = Review {
            input: PathBuf::from("variants.vcf"),
            consensus_bam: PathBuf::from("consensus.bam"),
            grouped_bam: PathBuf::from("grouped.bam"),
            reference: PathBuf::from("ref.fa"),
            output: PathBuf::from("/path/to/output"),
            sample: None,
            ignore_ns: false,
            maf: 0.05,
        };

        let consensus_path = review.output.with_extension("consensus.bam");
        let grouped_path = review.output.with_extension("grouped.bam");
        let review_path = review.output.with_extension("txt");

        assert_eq!(consensus_path, PathBuf::from("/path/to/output.consensus.bam"));
        assert_eq!(grouped_path, PathBuf::from("/path/to/output.grouped.bam"));
        assert_eq!(review_path, PathBuf::from("/path/to/output.txt"));
    }

    // ------------------------------------------------------------------
    // Unit tests for the fgbio-parity helpers (REV-01 genotype, REV-05 MI guard)
    // ------------------------------------------------------------------

    /// Each case pairs a genotype — a list of `(allele position, is_phased)` tuples — with the
    /// reference base, the ALT bases, and the htsjdk-formatted string `format_genotype_string`
    /// must produce. Using an `#[rstest]` case table reports each parity scenario independently
    /// on failure.
    #[rstest]
    // Unphased het 0/1 (REF=A, ALT=[T]) => base strings joined with `/`.
    #[case::unphased_het(&[(Some(0), false), (Some(1), false)], "A", &["T"], "A/T")]
    // htsjdk preserves genotype order (it does NOT sort), so 1/0 => "T/A".
    // Verified against fgbio ReviewConsensusVariants on a `GT=1/0` record.
    #[case::unphased_order_preserved(&[(Some(1), false), (Some(0), false)], "A", &["T"], "T/A")]
    // Phased 0|1 preserves order and joins with `|`.
    #[case::phased(&[(Some(0), true), (Some(1), true)], "A", &["T"], "A|T")]
    // Multi-allelic phased using the second alternate (index 2 => ALT[1]).
    #[case::phased_multiallelic(&[(Some(1), true), (Some(2), true)], "A", &["C", "G"], "C|G")]
    // A missing allele renders as `.`.
    #[case::missing_allele(&[(None, false), (Some(1), false)], "A", &["T"], "./T")]
    fn test_format_genotype_string_matches_htsjdk(
        #[case] alleles: &[(Option<usize>, bool)],
        #[case] reference: &str,
        #[case] alternates: &[&str],
        #[case] expected: &str,
    ) {
        use noodles::vcf::variant::record::samples::series::value::genotype::Phasing;
        use noodles::vcf::variant::record_buf::samples::sample::value::genotype::Allele;

        let genotype: Vec<Allele> = alleles
            .iter()
            .map(|&(position, phased)| {
                Allele::new(position, if phased { Phasing::Phased } else { Phasing::Unphased })
            })
            .collect();
        let alternates: Vec<String> = alternates.iter().map(|s| (*s).to_string()).collect();

        assert_eq!(format_genotype_string(&genotype, reference, &alternates), expected);
    }

    #[test]
    fn test_to_mi_truncates_and_errors_when_absent() {
        use fgumi_raw_bam::SamBuilder as RawSamBuilder;

        // Present MI with a non-/A,/B suffix is truncated at the last `/` (REV-02).
        let mut b = RawSamBuilder::new();
        b.read_name(b"withmi").sequence(b"A").qualities(&[30]);
        b.add_string_tag(SamTag::MI, b"5/456");
        assert_eq!(to_mi(&b.build()).expect("MI present"), "5");

        // A leading `/` is preserved (fgbio's `slash > 0` guard): only the *last* `/`
        // truncates, so `/abc/456` keeps the leading-slash prefix `/abc`.
        let mut b_lead = RawSamBuilder::new();
        b_lead.read_name(b"leadmi").sequence(b"A").qualities(&[30]);
        b_lead.add_string_tag(SamTag::MI, b"/abc/456");
        assert_eq!(to_mi(&b_lead.build()).expect("MI present"), "/abc");

        // Absent MI is an error, mirroring fgbio `toMi`
        // (ReviewConsensusVariantsTest.scala:166).
        let mut b2 = RawSamBuilder::new();
        b2.read_name(b"nomi").sequence(b"A").qualities(&[30]);
        let err = to_mi(&b2.build()).expect_err("missing MI must error");
        assert!(err.to_string().contains("nomi"), "error should name the read: {err}");
    }

    // Integration tests based on Scala test suite

    #[test]
    fn test_empty_vcf_produces_empty_outputs() {
        use noodles::bam;

        let temp_dir = TempDir::new().expect("failed to create temp dir");
        let ref_path = test_utils::create_test_reference(&temp_dir);
        let (raw_path, consensus_path) = test_utils::create_test_bams(&temp_dir);
        let vcf_path = test_utils::create_empty_vcf(&temp_dir);
        let output_path = temp_dir.path().join("output");

        let review = Review {
            input: vcf_path,
            consensus_bam: consensus_path,
            grouped_bam: raw_path,
            reference: ref_path,
            output: output_path.clone(),
            sample: None,
            ignore_ns: false,
            maf: 0.05,
        };

        review.execute("test").expect("execute should succeed");

        // Verify output files exist
        let con_out = output_path.with_extension("consensus.bam");
        let raw_out = output_path.with_extension("grouped.bam");
        let txt_out = output_path.with_extension("txt");

        assert!(con_out.exists());
        assert!(raw_out.exists());
        assert!(txt_out.exists());

        // Verify BAMs are empty
        let mut con_reader = bam::io::indexed_reader::Builder::default()
            .build_from_path(&con_out)
            .expect("failed to open indexed BAM");
        let con_header = con_reader.read_header().expect("failed to read BAM header");
        let mut con_record = noodles::sam::alignment::RecordBuf::default();
        assert_eq!(
            con_reader
                .read_record_buf(&con_header, &mut con_record)
                .expect("failed to read BAM record"),
            0
        );

        let mut raw_reader = bam::io::indexed_reader::Builder::default()
            .build_from_path(&raw_out)
            .expect("failed to open indexed BAM");
        let raw_header = raw_reader.read_header().expect("failed to read BAM header");
        let mut raw_record = noodles::sam::alignment::RecordBuf::default();
        assert_eq!(
            raw_reader
                .read_record_buf(&raw_header, &mut raw_record)
                .expect("failed to read BAM record"),
            0
        );
    }

    #[test]
    fn test_empty_interval_list_produces_empty_outputs() {
        use noodles::bam;

        let temp_dir = TempDir::new().expect("failed to create temp dir");
        let ref_path = test_utils::create_test_reference(&temp_dir);
        let (raw_path, consensus_path) = test_utils::create_test_bams(&temp_dir);
        let interval_path = test_utils::create_empty_interval_list(&temp_dir);
        let output_path = temp_dir.path().join("output");

        let review = Review {
            input: interval_path,
            consensus_bam: consensus_path,
            grouped_bam: raw_path,
            reference: ref_path,
            output: output_path.clone(),
            sample: None,
            ignore_ns: false,
            maf: 0.05,
        };

        review.execute("test").expect("execute should succeed");

        // Verify output files exist and are empty
        let con_out = output_path.with_extension("consensus.bam");
        let raw_out = output_path.with_extension("grouped.bam");

        assert!(con_out.exists());
        assert!(raw_out.exists());

        // Verify BAMs are empty
        let mut con_reader = bam::io::indexed_reader::Builder::default()
            .build_from_path(&con_out)
            .expect("failed to open indexed BAM");
        let con_header = con_reader.read_header().expect("failed to read BAM header");
        let mut con_record = noodles::sam::alignment::RecordBuf::default();
        assert_eq!(
            con_reader
                .read_record_buf(&con_header, &mut con_record)
                .expect("failed to read BAM record"),
            0
        );
    }

    #[test]
    fn test_extracts_correct_reads_for_variants() {
        use noodles::bam;

        let temp_dir = TempDir::new().expect("failed to create temp dir");
        let ref_path = test_utils::create_test_reference(&temp_dir);
        let (raw_path, consensus_path) = test_utils::create_test_bams(&temp_dir);
        let vcf_path = test_utils::create_test_vcf(&temp_dir);
        let output_path = temp_dir.path().join("output");

        let review = Review {
            input: vcf_path,
            consensus_bam: consensus_path,
            grouped_bam: raw_path,
            reference: ref_path,
            output: output_path.clone(),
            sample: Some("tumor".to_string()),
            ignore_ns: false,
            maf: 0.05,
        };

        review.execute("test").expect("execute should succeed");

        // Verify output files exist
        let con_out = output_path.with_extension("consensus.bam");
        let raw_out = output_path.with_extension("grouped.bam");
        let txt_out = output_path.with_extension("txt");

        assert!(con_out.exists());
        assert!(raw_out.exists());
        assert!(txt_out.exists());

        // Read consensus BAM and verify read names
        let mut con_reader = bam::io::indexed_reader::Builder::default()
            .build_from_path(&con_out)
            .expect("failed to open indexed BAM");
        let con_header = con_reader.read_header().expect("failed to read BAM header");

        let mut consensus_reads = Vec::new();
        let mut con_record = noodles::sam::alignment::RecordBuf::default();
        while con_reader
            .read_record_buf(&con_header, &mut con_record)
            .expect("failed to read BAM record")
            > 0
        {
            let name = String::from_utf8_lossy(
                con_record.name().expect("record should have name").as_ref(),
            )
            .to_string();
            consensus_reads.push(name);
        }

        // Should contain reads A, B, D, E, F, H (first and second of pair)
        // Note: Based on Scala test, we expect specific reads
        assert!(consensus_reads.contains(&"A".to_string()));
        assert!(consensus_reads.contains(&"B".to_string()));
        assert!(consensus_reads.contains(&"E".to_string()));
        assert!(consensus_reads.contains(&"F".to_string()));
        assert!(consensus_reads.contains(&"H".to_string()));

        // Read raw BAM and verify read names
        let mut raw_reader = bam::io::indexed_reader::Builder::default()
            .build_from_path(&raw_out)
            .expect("failed to open indexed BAM");
        let raw_header = raw_reader.read_header().expect("failed to read BAM header");

        let mut raw_reads = Vec::new();
        let mut raw_record = noodles::sam::alignment::RecordBuf::default();
        while raw_reader
            .read_record_buf(&raw_header, &mut raw_record)
            .expect("failed to read BAM record")
            > 0
        {
            let name = String::from_utf8_lossy(
                raw_record.name().expect("record should have name").as_ref(),
            )
            .to_string();
            raw_reads.push(name);
        }

        // Should contain A1, A2, B1, D1, E1, F1, H1
        assert!(raw_reads.contains(&"A1".to_string()));
        assert!(raw_reads.contains(&"A2".to_string()));
        assert!(raw_reads.contains(&"B1".to_string()));
        assert!(raw_reads.contains(&"E1".to_string()));
        assert!(raw_reads.contains(&"F1".to_string()));
        assert!(raw_reads.contains(&"H1".to_string()));
    }

    #[test]
    fn test_review_tsv_contains_correct_information() {
        use std::io::BufRead;

        let temp_dir = TempDir::new().expect("failed to create temp dir");
        let ref_path = test_utils::create_test_reference(&temp_dir);
        let (raw_path, consensus_path) = test_utils::create_test_bams(&temp_dir);
        let vcf_path = test_utils::create_test_vcf(&temp_dir);
        let output_path = temp_dir.path().join("output");

        let review = Review {
            input: vcf_path,
            consensus_bam: consensus_path,
            grouped_bam: raw_path,
            reference: ref_path,
            output: output_path.clone(),
            sample: Some("tumor".to_string()),
            ignore_ns: false,
            maf: 0.05,
        };

        review.execute("test").expect("execute should succeed");

        let txt_out = output_path.with_extension("txt");
        assert!(txt_out.exists());

        // Read TSV and verify it has content
        let file = std::fs::File::open(&txt_out).expect("failed to open file");
        let reader = std::io::BufReader::new(file);
        let lines: Vec<String> = reader.lines().map(|l| l.expect("failed to read line")).collect();

        // Should have header plus data rows
        assert!(lines.len() > 1, "TSV should have header and data rows");

        // Check header contains expected columns. The reference-allele column is
        // named `ref` to match fgbio's `ConsensusVariantReviewInfo` (REV-04), not
        // `ref_allele`.
        let header = &lines[0];
        let columns: Vec<&str> = header.split('\t').collect();
        assert!(columns.contains(&"chrom"));
        assert!(columns.contains(&"pos"));
        assert!(columns.contains(&"ref"), "reference column must be `ref`, got: {header}");
        assert!(!columns.contains(&"ref_allele"), "column must not be `ref_allele`: {header}");
        assert!(columns.contains(&"consensus_read"));
        assert!(columns.contains(&"consensus_call"));
    }

    #[test]
    fn test_order_variants_by_dictionary() {
        // create_test_header defines chr1 then chr2.
        let header = test_utils::create_test_header();

        // Variants in a scrambled (non-coordinate) order.
        let variants = vec![
            Variant::new("chr2".to_string(), 20, 'C'),
            Variant::new("chr1".to_string(), 30, 'A'),
            Variant::new("chr1".to_string(), 10, 'A'),
            Variant::new("chr1".to_string(), 20, 'A'),
        ];

        let sorted = Review::order_variants_by_dictionary(variants, &header)
            .expect("ordering should succeed");
        let got: Vec<(&str, i32)> = sorted.iter().map(|v| (v.chrom.as_str(), v.pos)).collect();
        assert_eq!(
            got,
            vec![("chr1", 10), ("chr1", 20), ("chr1", 30), ("chr2", 20)],
            "variants must be ordered by dictionary index then position"
        );

        // A contig absent from the dictionary is an error, matching fgbio's requirement
        // that every variant interval exists in the sequence dictionary.
        let bad = vec![Variant::new("chrX".to_string(), 1, 'A')];
        assert!(
            Review::order_variants_by_dictionary(bad, &header).is_err(),
            "an unknown contig must be rejected"
        );
    }

    #[test]
    fn test_review_output_row_order_matches_fgbio() {
        use std::io::BufRead;

        let temp_dir = TempDir::new().expect("failed to create temp dir");
        let ref_path = test_utils::create_test_reference(&temp_dir);

        // chr1:5 has one supporting read; chr1:10 has three with MIs "2", "3", "10"
        // written in that (non-sorted) order. MI "10" additionally carries a non-reference
        // R2 over chr1:10, so it emits both a `/1` and a `/2` row. fgbio emits variants in
        // coordinate order (chr1:5 before chr1:10) and sorts each variant's rows by
        // MI+read-number as strings, so the chr1:10 rows come out
        // "10/1" < "10/2" < "2/1" < "3/1" — exercising the /1-before-/2 tie-break.
        let reads = [
            ("r_p5", 0usize, 1i32, b'C', "1", None),
            ("r_mi2", 0usize, 6i32, b'T', "2", None),
            ("r_mi3", 0usize, 6i32, b'T', "3", None),
            ("r_mi10", 0usize, 6i32, b'T', "10", Some((6i32, b'T'))),
        ];
        let (grouped_path, consensus_path) = test_utils::create_ordering_bams(&temp_dir, &reads);
        let vcf_path = test_utils::create_ordering_vcf(&temp_dir);
        let output_path = temp_dir.path().join("output");

        let review = Review {
            input: vcf_path,
            consensus_bam: consensus_path,
            grouped_bam: grouped_path,
            reference: ref_path,
            output: output_path.clone(),
            sample: None,
            ignore_ns: false,
            maf: 0.05,
        };
        review.execute("test").expect("execute should succeed");

        let txt_out = output_path.with_extension("txt");
        let file = std::fs::File::open(&txt_out).expect("failed to open review txt");
        let lines: Vec<String> = std::io::BufReader::new(file)
            .lines()
            .map(|l| l.expect("failed to read line"))
            .collect();

        let cols: Vec<&str> = lines[0].split('\t').collect();
        let pos_idx = cols.iter().position(|c| *c == "pos").expect("pos column");
        let read_idx =
            cols.iter().position(|c| *c == "consensus_read").expect("consensus_read column");

        let rows: Vec<(String, String)> = lines[1..]
            .iter()
            .map(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                (f[pos_idx].to_string(), f[read_idx].to_string())
            })
            .collect();

        let expected = vec![
            ("5".to_string(), "r_p5/1".to_string()),
            ("10".to_string(), "r_mi10/1".to_string()),
            ("10".to_string(), "r_mi10/2".to_string()),
            ("10".to_string(), "r_mi2/1".to_string()),
            ("10".to_string(), "r_mi3/1".to_string()),
        ];
        assert_eq!(
            rows, expected,
            "review rows must match fgbio coordinate + MI/read-number ordering"
        );
        // The same MI ("10") appears with both read numbers, so /1 must sort before /2.
        let mi10: Vec<&String> = rows
            .iter()
            .filter(|(pos, read)| pos == "10" && read.starts_with("r_mi10/"))
            .map(|(_, read)| read)
            .collect();
        assert_eq!(
            mi10,
            vec!["r_mi10/1", "r_mi10/2"],
            "for the same MI, the /1 row must sort before the /2 row"
        );
    }

    #[test]
    fn test_spanning_deletions_handled_correctly() {
        use noodles::bam;

        let temp_dir = TempDir::new().expect("failed to create temp dir");
        let ref_path = test_utils::create_test_reference(&temp_dir);
        let (raw_path, consensus_path) = test_utils::create_test_bams(&temp_dir);
        let vcf_path = test_utils::create_test_vcf(&temp_dir);
        let output_path = temp_dir.path().join("output");

        let review = Review {
            input: vcf_path,
            consensus_bam: consensus_path,
            grouped_bam: raw_path,
            reference: ref_path,
            output: output_path.clone(),
            sample: Some("tumor".to_string()),
            ignore_ns: false,
            maf: 0.05,
        };

        review.execute("test").expect("execute should succeed");

        // Read consensus BAM - should include D (spanning deletion)
        let con_out = output_path.with_extension("consensus.bam");
        let mut con_reader = bam::io::indexed_reader::Builder::default()
            .build_from_path(&con_out)
            .expect("failed to open indexed BAM");
        let con_header = con_reader.read_header().expect("failed to read BAM header");

        let mut consensus_reads = Vec::new();
        let mut con_record = noodles::sam::alignment::RecordBuf::default();
        while con_reader
            .read_record_buf(&con_header, &mut con_record)
            .expect("failed to read BAM record")
            > 0
        {
            let name = String::from_utf8_lossy(
                con_record.name().expect("record should have name").as_ref(),
            )
            .to_string();
            consensus_reads.push(name);
        }

        // D should be extracted (spanning deletion at chr1:30)
        assert!(consensus_reads.contains(&"D".to_string()));

        // REV-03: but the spanning-deletion read must NOT get a row in the review
        // .txt (fgbio's SamLocusIterator.getRecordAndOffsets excludes deletions;
        // ReviewConsensusVariantsTest.scala:250-251).
        let txt_out = output_path.with_extension("txt");
        let content = std::fs::read_to_string(&txt_out).expect("failed to read review txt");
        assert!(
            !content.contains("\tD/1\t"),
            "spanning-deletion read D should not have a review row:\n{content}"
        );
        assert!(
            !content.contains("\t*\t"),
            "no row should carry a spanning-deletion `*` call:\n{content}"
        );

        // REV-03: the grouped/raw spanning-deletion read (D1, MI=D) is still extracted to
        // the grouped output BAM — only the review .txt row is suppressed, not the read
        // itself (fgbio `readBamRecs(rawOut)` contains `D1/1`;
        // ReviewConsensusVariantsTest.scala:249).
        let grouped_out = output_path.with_extension("grouped.bam");
        let mut grp_reader = bam::io::indexed_reader::Builder::default()
            .build_from_path(&grouped_out)
            .expect("failed to open grouped BAM");
        let grp_header = grp_reader.read_header().expect("failed to read grouped BAM header");
        let mut grouped_reads = Vec::new();
        let mut grp_record = noodles::sam::alignment::RecordBuf::default();
        while grp_reader
            .read_record_buf(&grp_header, &mut grp_record)
            .expect("failed to read grouped BAM record")
            > 0
        {
            grouped_reads.push(
                String::from_utf8_lossy(
                    grp_record.name().expect("record should have name").as_ref(),
                )
                .to_string(),
            );
        }
        assert!(
            grouped_reads.contains(&"D1".to_string()),
            "spanning-deletion grouped read D1 must still be extracted: {grouped_reads:?}"
        );
    }

    #[test]
    fn test_no_calls_handled_correctly() {
        let temp_dir = TempDir::new().expect("failed to create temp dir");
        let ref_path = test_utils::create_test_reference(&temp_dir);
        let (raw_path, consensus_path) = test_utils::create_test_bams(&temp_dir);
        let vcf_path = test_utils::create_test_vcf(&temp_dir);
        let output_path = temp_dir.path().join("output");

        let review = Review {
            input: vcf_path.clone(),
            consensus_bam: consensus_path.clone(),
            grouped_bam: raw_path.clone(),
            reference: ref_path.clone(),
            output: output_path.clone(),
            sample: Some("tumor".to_string()),
            ignore_ns: false, // Include N bases
            maf: 0.05,
        };

        review.execute("test").expect("execute should succeed");

        // Read TSV and verify N bases are included
        let txt_out = output_path.with_extension("txt");
        let content = std::fs::read_to_string(&txt_out).expect("failed to read file");

        // Should find E consensus read with N base
        assert!(content.contains("E/"), "Should contain consensus read E");

        // Now test with ignore_ns enabled
        let output_path2 = temp_dir.path().join("output2");
        let review2 = Review {
            input: vcf_path,
            consensus_bam: consensus_path,
            grouped_bam: raw_path,
            reference: ref_path,
            output: output_path2.clone(),
            sample: Some("tumor".to_string()),
            ignore_ns: true, // Ignore N bases
            maf: 0.05,
        };

        review2.execute("test").expect("execute should succeed");

        let txt_out2 = output_path2.with_extension("txt");
        let content2 = std::fs::read_to_string(&txt_out2).expect("failed to read file");

        // E should not be in the TSV (N bases ignored)
        let lines: Vec<&str> = content2.lines().collect();
        let e_count = lines.iter().filter(|l| l.contains("E/")).count();
        assert_eq!(e_count, 0, "Should not contain consensus read E when ignoring Ns");
    }

    #[test]
    fn test_both_ends_overlapping_variant() {
        use noodles::bam;

        let temp_dir = TempDir::new().expect("failed to create temp dir");
        let ref_path = test_utils::create_test_reference(&temp_dir);
        let (raw_path, consensus_path) = test_utils::create_test_bams(&temp_dir);
        let vcf_path = test_utils::create_test_vcf(&temp_dir);
        let output_path = temp_dir.path().join("output");

        let review = Review {
            input: vcf_path,
            consensus_bam: consensus_path,
            grouped_bam: raw_path,
            reference: ref_path,
            output: output_path.clone(),
            sample: Some("tumor".to_string()),
            ignore_ns: false,
            maf: 0.05,
        };

        review.execute("test").expect("execute should succeed");

        // Read consensus BAM - should include both H/1 and H/2
        let con_out = output_path.with_extension("consensus.bam");
        let mut con_reader = bam::io::indexed_reader::Builder::default()
            .build_from_path(&con_out)
            .expect("failed to open indexed BAM");
        let con_header = con_reader.read_header().expect("failed to read BAM header");

        let mut h_reads = 0;
        let mut con_record = noodles::sam::alignment::RecordBuf::default();
        while con_reader
            .read_record_buf(&con_header, &mut con_record)
            .expect("failed to read BAM record")
            > 0
        {
            let name = String::from_utf8_lossy(
                con_record.name().expect("record should have name").as_ref(),
            )
            .to_string();
            if name == "H" {
                h_reads += 1;
            }
        }

        // Should have both /1 and /2 reads for H
        assert_eq!(h_reads, 2, "Should extract both ends of pair H overlapping variant at chr2:20");
    }

    #[test]
    fn test_missing_fasta_index_fails() {
        let temp_dir = TempDir::new().expect("failed to create temp dir");
        let ref_path = temp_dir.path().join("ref.fa");

        // Create FASTA without .fai
        let fasta_content = b">chr1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n";
        std::fs::write(&ref_path, fasta_content).expect("failed to write file");

        let (raw_path, consensus_path) = test_utils::create_test_bams(&temp_dir);
        let vcf_path = test_utils::create_empty_vcf(&temp_dir);
        let output_path = temp_dir.path().join("output");

        let review = Review {
            input: vcf_path,
            consensus_bam: consensus_path,
            grouped_bam: raw_path,
            reference: ref_path,
            output: output_path,
            sample: None,
            ignore_ns: false,
            maf: 0.05,
        };

        // Should fail due to missing .fai
        let result = review.execute("test");
        assert!(result.is_err(), "Should fail when FASTA index is missing");
    }

    #[test]
    fn test_missing_fasta_dict_fails() {
        let temp_dir = TempDir::new().expect("failed to create temp dir");
        let ref_path = temp_dir.path().join("ref.fa");
        let fai_path = temp_dir.path().join("ref.fa.fai");

        // Create FASTA with .fai but without .dict
        let fasta_content = b">chr1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n";
        std::fs::write(&ref_path, fasta_content).expect("failed to write file");

        let fai_content = b"chr1\t100\t6\t100\t101\n";
        std::fs::write(&fai_path, fai_content).expect("failed to write file");

        let (raw_path, consensus_path) = test_utils::create_test_bams(&temp_dir);
        let vcf_path = test_utils::create_empty_vcf(&temp_dir);
        let output_path = temp_dir.path().join("output");

        let review = Review {
            input: vcf_path,
            consensus_bam: consensus_path,
            grouped_bam: raw_path,
            reference: ref_path,
            output: output_path,
            sample: None,
            ignore_ns: false,
            maf: 0.05,
        };

        // Should fail due to missing .dict (now that we validate it)
        let result = review.execute("test");
        assert!(result.is_err(), "Should fail when FASTA dictionary is missing");
    }
}

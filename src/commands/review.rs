//! Review consensus variant calls by extracting supporting reads.
//!
//! This tool extracts consensus reads containing variants and their supporting
//! raw reads to facilitate manual review of variant calls. It creates filtered
//! BAM files and a detailed TSV report.

use anyhow::{Result, bail};
use clap::Parser;
use fgumi_lib::logging::OperationTimer;
use fgumi_lib::progress::ProgressTracker;
use fgumi_lib::umi::extract_mi_base;
use fgumi_lib::validation::validate_file_exists;
use fgumi_lib::variant_review::{
    BaseCounts, ConsensusVariantReviewInfo, Variant, format_insert_string, read_number_suffix,
};
use log::info;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record::Cigar;
use noodles::sam::alignment::record_buf::RecordBuf;
use std::collections::HashSet;
use std::path::{Path, PathBuf};

use super::command::Command;

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
    #[arg(short = 'N', long = "ignore-ns", default_value = "false")]
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
        let dict_path = self.reference.with_extension("dict");
        if !dict_path.exists() {
            bail!(
                "The reference file has no sequence dictionary. Please run:\n  \
                samtools dict {} -o {}",
                self.reference.display(),
                dict_path.display()
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

    /// Extract genotype for a specific sample
    ///
    /// Returns the genotype string (e.g., "0/1", "1|1") for the specified sample.
    /// Implementation note: This is simplified due to noodles API complexity.
    /// For production use with complex VCF files, consider enhancing with proper genotype parsing.
    fn extract_genotype_for_sample(
        record: &noodles::vcf::variant::RecordBuf,
        header: &noodles::vcf::Header,
        sample_name: &str,
    ) -> Option<String> {
        use noodles::vcf::variant::record::samples::Sample;

        // Find the sample index
        let sample_idx = header.sample_names().iter().position(|s| s == sample_name)?;

        // Get the samples
        let samples = record.samples();
        let keys = samples.keys();

        // Try to get the GT field for this sample
        let gt_key = "GT";

        // Iterate through sample keys to find GT
        if let Some(gt_idx) = keys.as_ref().iter().position(|k| k == gt_key) {
            // Access the sample value at this index
            if let Some(sample_values) = samples.get_index(sample_idx) {
                if let Some(Ok(Some(gt_value))) = sample_values.get_index(header, gt_idx) {
                    // Convert to string representation
                    return Some(format!("{gt_value:?}"));
                }
            }
        }

        None
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

        let mut reader =
            bam::io::indexed_reader::Builder::default().build_from_path(&self.consensus_bam)?;
        let header = reader.read_header()?;

        // Add @PG record with PP chaining
        let header = fgumi_lib::header::add_pg_record(
            header,
            crate::version::VERSION.as_str(),
            command_line,
        )?;

        // Create output BAM
        let consensus_out_path = self.output.with_extension("consensus.bam");
        let mut writer = bam::io::writer::Builder.build_from_path(&consensus_out_path)?;
        writer.write_header(&header)?;

        let mut mi_set = HashSet::new();

        // Query each variant region
        for variant in variants {
            let start = noodles::core::Position::try_from(variant.pos as usize)?;
            let region = noodles::core::Region::new(variant.chrom.as_str(), start..=start);

            let query = reader.query(&header, &region)?;

            for result in query.records() {
                let bam_record = result?;
                // Convert to RecordBuf for processing
                let record = RecordBuf::try_from_alignment_record(&header, &bam_record)?;

                // Check if this read has a non-reference base at the variant position
                if self.has_non_reference_base(&record, variant)? {
                    // Extract MI tag
                    let mi_tag =
                        noodles::sam::alignment::record::data::field::Tag::from([b'M', b'I']);
                    if let Some(noodles::sam::alignment::record_buf::data::field::Value::String(
                        mi_bytes,
                    )) = record.data().get(&mi_tag)
                    {
                        let mi = String::from_utf8(mi_bytes.iter().copied().collect())?;
                        let mi_base = extract_mi_base(&mi).to_string();
                        mi_set.insert(mi_base);

                        // Write to output BAM
                        writer.write_alignment_record(&header, &record)?;
                    }
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

        let mut reader =
            bam::io::indexed_reader::Builder::default().build_from_path(&self.grouped_bam)?;
        let header = reader.read_header()?;

        // Add @PG record with PP chaining
        let header = fgumi_lib::header::add_pg_record(
            header,
            crate::version::VERSION.as_str(),
            command_line,
        )?;

        let grouped_out_path = self.output.with_extension("grouped.bam");
        let mut writer = bam::io::writer::Builder.build_from_path(&grouped_out_path)?;
        writer.write_header(&header)?;

        // Read through entire grouped BAM and extract matching reads
        let progress = ProgressTracker::new("Processed grouped reads").with_interval(1_000_000);
        let mut record = RecordBuf::default();
        while reader.read_record_buf(&header, &mut record)? != 0 {
            progress.log_if_needed(1);

            // Extract MI tag
            let mi_tag = noodles::sam::alignment::record::data::field::Tag::from([b'M', b'I']);
            if let Some(noodles::sam::alignment::record_buf::data::field::Value::String(mi_bytes)) =
                record.data().get(&mi_tag)
            {
                let mi = String::from_utf8(mi_bytes.iter().copied().collect())?;
                let mi_base = extract_mi_base(&mi);

                if mi_set.contains(mi_base) {
                    writer.write_alignment_record(&header, &record)?;
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
        use noodles::bam;

        let review_path = self.output.with_extension("txt");

        // Open both BAMs
        let mut consensus_reader =
            bam::io::indexed_reader::Builder::default().build_from_path(&self.consensus_bam)?;
        let consensus_header = consensus_reader.read_header()?;

        let mut grouped_reader =
            bam::io::indexed_reader::Builder::default().build_from_path(&self.grouped_bam)?;
        let grouped_header = grouped_reader.read_header()?;

        let mut all_metrics = Vec::new();

        // Process each variant
        for variant in variants {
            // Query consensus BAM for this position
            let start = noodles::core::Position::try_from(variant.pos as usize)?;
            let region = noodles::core::Region::new(variant.chrom.as_str(), start..=start);

            let query = consensus_reader.query(&consensus_header, &region)?;

            // Collect all consensus reads at this position
            let mut consensus_reads = Vec::new();
            for result in query.records() {
                let bam_record = result?;
                let record = RecordBuf::try_from_alignment_record(&consensus_header, &bam_record)?;
                consensus_reads.push(record);
            }

            // Build consensus-level base counts
            let mut consensus_counts = BaseCounts::default();
            for record in &consensus_reads {
                if let Some(base) = self.get_base_at_position(record, variant.pos)? {
                    consensus_counts.add_base(base);
                }
            }

            // Process each consensus read
            for record in consensus_reads {
                // Get the base at the variant position
                // If None, check if it's a spanning deletion
                let read_base = match self.get_base_at_position(&record, variant.pos)? {
                    Some(b) => (b as char).to_ascii_uppercase(),
                    None => {
                        // Check if this is a spanning deletion (position is within a deletion)
                        if self.is_spanning_deletion(&record, variant.pos)? {
                            '*' // Use '*' to represent spanning deletion
                        } else {
                            continue; // Position not covered by read
                        }
                    }
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

                // Extract MI tag
                let mi_tag = noodles::sam::alignment::record::data::field::Tag::from([b'M', b'I']);
                let mi_str = match record.data().get(&mi_tag) {
                    Some(noodles::sam::alignment::record_buf::data::field::Value::String(
                        mi_bytes,
                    )) => String::from_utf8(mi_bytes.iter().copied().collect())?,
                    _ => continue,
                };

                let mi_base = extract_mi_base(&mi_str);

                // Format consensus read name with /1 or /2 suffix
                let read_name =
                    String::from_utf8_lossy(record.name().map_or(b"unnamed", |n| n.as_ref()))
                        .to_string();
                let is_first = record.flags().is_first_segment();
                let consensus_read_name = format!("{}{}", read_name, read_number_suffix(is_first));

                // Generate insert string
                let insert_string = format_insert_string(&record, &consensus_header);

                // Query grouped BAM for raw reads with matching MI
                let raw_counts = {
                    use std::collections::HashSet;

                    let start = noodles::core::Position::try_from(variant.pos as usize)?;
                    let region = noodles::core::Region::new(variant.chrom.as_str(), start..=start);
                    let query = grouped_reader.query(&grouped_header, &region)?;

                    let mut counts = BaseCounts::default();
                    let mut seen_reads = HashSet::new();
                    let expected_read_num =
                        if consensus_read_name.ends_with("/1") { "/1" } else { "/2" };

                    for result in query.records() {
                        let bam_record = result?;
                        let record =
                            RecordBuf::try_from_alignment_record(&grouped_header, &bam_record)?;

                        // Extract MI tag
                        let mi_tag =
                            noodles::sam::alignment::record::data::field::Tag::from([b'M', b'I']);
                        let mi_str = match record.data().get(&mi_tag) {
                            Some(
                                noodles::sam::alignment::record_buf::data::field::Value::String(
                                    mi_bytes,
                                ),
                            ) => String::from_utf8(mi_bytes.iter().copied().collect())?,
                            _ => continue,
                        };

                        // Check if MI matches
                        let read_mi_base = extract_mi_base(&mi_str);
                        if read_mi_base != mi_base {
                            continue;
                        }

                        // Check read number matches
                        let is_first = record.flags().is_first_segment();
                        let read_num = read_number_suffix(is_first);
                        if read_num != expected_read_num {
                            continue;
                        }

                        // Ensure we only count each read once
                        let read_name = String::from_utf8_lossy(
                            record.name().map_or(b"unnamed", |n| n.as_ref()),
                        )
                        .to_string();
                        if !seen_reads.insert(read_name) {
                            continue;
                        }

                        // Get base at position
                        if let Some(base) = self.get_base_at_position(&record, variant.pos)? {
                            counts.add_base(base);
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

                all_metrics.push(info);
            }
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

    /// Get the base at a specific reference position in a read
    fn get_base_at_position(&self, record: &RecordBuf, ref_pos: i32) -> Result<Option<u8>> {
        if let Some(offset) = self.calculate_read_offset(record, ref_pos)? {
            let sequence = record.sequence();
            if offset < sequence.len() {
                return Ok(Some(sequence.get(offset).unwrap()));
            }
        }
        Ok(None)
    }

    /// Get the quality score at a specific reference position in a read
    fn get_quality_at_position(&self, record: &RecordBuf, ref_pos: i32) -> Result<Option<u8>> {
        if let Some(offset) = self.calculate_read_offset(record, ref_pos)? {
            let quality_scores = record.quality_scores();
            if offset < quality_scores.len() {
                return Ok(Some(quality_scores.as_ref()[offset]));
            }
        }
        Ok(None)
    }

    /// Check if a record has a non-reference base at the variant position
    fn has_non_reference_base(&self, record: &RecordBuf, variant: &Variant) -> Result<bool> {
        // Skip unmapped reads
        if record.flags().is_unmapped() {
            return Ok(false);
        }

        // Get alignment start
        let start = match record.alignment_start() {
            Some(pos) => usize::from(pos) as i32,
            None => return Ok(false),
        };

        // Get alignment end
        let end = match record.alignment_end() {
            Some(pos) => usize::from(pos) as i32,
            None => return Ok(false),
        };

        // Check if variant position overlaps this read
        if variant.pos < start || variant.pos > end {
            return Ok(false);
        }

        // Get the base at the variant position using CIGAR-aware mapping
        let read_offset = self.calculate_read_offset(record, variant.pos)?;

        if let Some(offset) = read_offset {
            let sequence = record.sequence();
            if offset < sequence.len() {
                let base = sequence.get(offset).unwrap();
                let base_char = (base as char).to_ascii_uppercase();

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
            // read_offset is None - check if it's a spanning deletion
            if self.is_spanning_deletion(record, variant.pos)? {
                return Ok(true); // Spanning deletions are considered non-reference
            }
        }

        Ok(false)
    }

    /// Calculate read offset for a reference position using CIGAR
    /// Check if a reference position falls within a deletion in the read
    fn is_spanning_deletion(&self, record: &RecordBuf, ref_pos: i32) -> Result<bool> {
        use noodles::sam::alignment::record::cigar::op::Kind;

        let start = match record.alignment_start() {
            Some(pos) => usize::from(pos) as i32,
            None => return Ok(false),
        };

        // Check if position is within the read's reference span
        let end = match record.alignment_end() {
            Some(pos) => usize::from(pos) as i32,
            None => return Ok(false),
        };

        if ref_pos < start || ref_pos > end {
            return Ok(false); // Position not covered by read at all
        }

        let mut current_ref_pos = start;

        let cigar = record.cigar();
        let ops: Vec<_> = cigar.iter().filter_map(std::result::Result::ok).collect();

        for op in ops {
            let len = op.len();
            let kind = op.kind();

            match kind {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    current_ref_pos += len as i32;
                }
                Kind::Deletion | Kind::Skip => {
                    // Check if target position is within this deletion
                    // Deletions span from current_ref_pos to current_ref_pos + len - 1
                    if ref_pos >= current_ref_pos && ref_pos < current_ref_pos + len as i32 {
                        return Ok(true);
                    }
                    current_ref_pos += len as i32;
                }
                Kind::Insertion | Kind::SoftClip | Kind::HardClip | Kind::Pad => {
                    // These don't affect reference positions
                }
            }
        }

        Ok(false)
    }

    fn calculate_read_offset(&self, record: &RecordBuf, ref_pos: i32) -> Result<Option<usize>> {
        use noodles::sam::alignment::record::cigar::op::Kind;

        let start = match record.alignment_start() {
            Some(pos) => usize::from(pos) as i32,
            None => return Ok(None),
        };

        let mut current_ref_pos = start;
        let mut current_read_pos = 0;

        let cigar = record.cigar();
        let ops: Vec<_> = cigar.iter().filter_map(std::result::Result::ok).collect();

        for op in ops {
            let len = op.len();
            let kind = op.kind();

            match kind {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    if current_ref_pos + len as i32 > ref_pos {
                        // Target position is in this op
                        let offset_in_op = (ref_pos - current_ref_pos) as usize;
                        return Ok(Some(current_read_pos + offset_in_op));
                    }
                    current_ref_pos += len as i32;
                    current_read_pos += len;
                }
                Kind::Insertion | Kind::SoftClip => {
                    current_read_pos += len;
                }
                Kind::Deletion | Kind::Skip => {
                    if current_ref_pos + len as i32 > ref_pos {
                        // Position is deleted in this read
                        return Ok(None);
                    }
                    current_ref_pos += len as i32;
                }
                Kind::HardClip | Kind::Pad => {
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
    use std::path::PathBuf;
    use tempfile::TempDir;

    // Test utilities for creating synthetic test data
    mod test_utils {
        use super::*;
        use fgumi_lib::sam::builder::RecordBuilder;
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

            std::fs::write(&ref_path, fasta_content).unwrap();

            // Create FAI
            let fai_content = b"chr1\t100\t6\t100\t101\nchr2\t100\t114\t100\t101\n";
            std::fs::write(&fai_path, fai_content).unwrap();

            // Create dictionary
            let dict_content = b"@HD\tVN:1.5\tSO:unsorted\n\
@SQ\tSN:chr1\tLN:100\n\
@SQ\tSN:chr2\tLN:100\n";
            std::fs::write(&dict_path, dict_content).unwrap();

            ref_path
        }

        /// Create a SAM header with sequence dictionary
        pub fn create_test_header() -> Header {
            use bstr::BString;

            // Create header with coordinate sort order
            let mut header = Header::builder()
                .add_reference_sequence(
                    BString::from("chr1"),
                    Map::<ReferenceSequence>::new(NonZeroUsize::try_from(100).unwrap()),
                )
                .add_reference_sequence(
                    BString::from("chr2"),
                    Map::<ReferenceSequence>::new(NonZeroUsize::try_from(100).unwrap()),
                )
                .build();

            // Parse and set the sort order in the header line
            let header_str = "@HD\tVN:1.6\tSO:coordinate\n";
            let parsed_header: Header = header_str.parse().unwrap();
            if let Some(hd) = parsed_header.header() {
                header = Header::builder()
                    .set_header(hd.clone())
                    .add_reference_sequence(
                        BString::from("chr1"),
                        Map::<ReferenceSequence>::new(NonZeroUsize::try_from(100).unwrap()),
                    )
                    .add_reference_sequence(
                        BString::from("chr2"),
                        Map::<ReferenceSequence>::new(NonZeroUsize::try_from(100).unwrap()),
                    )
                    .build();
            }

            header
        }

        /// Create a simple read pair
        #[allow(clippy::too_many_arguments)]
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
            use noodles::sam::alignment::record::Flags;

            let cigar1_str = cigar1.map_or_else(|| format!("{}M", bases1.len()), String::from);

            let r1 = RecordBuilder::new()
                .name(name)
                .flags(Flags::SEGMENTED | Flags::FIRST_SEGMENT)
                .reference_sequence_id(chrom_idx)
                .alignment_start(start1 as usize)
                .mapping_quality(60)
                .cigar(&cigar1_str)
                .sequence(&String::from_utf8_lossy(bases1))
                .qualities(&vec![45u8; bases1.len()])
                .mate_reference_sequence_id(chrom_idx)
                .mate_alignment_start(start2 as usize)
                .tag("MI", mi)
                .build();

            let r2 = RecordBuilder::new()
                .name(name)
                .flags(Flags::SEGMENTED | Flags::LAST_SEGMENT | Flags::REVERSE_COMPLEMENTED)
                .reference_sequence_id(chrom_idx)
                .alignment_start(start2 as usize)
                .mapping_quality(60)
                .cigar(&format!("{}M", bases2.len()))
                .sequence(&String::from_utf8_lossy(bases2))
                .qualities(&vec![45u8; bases2.len()])
                .mate_reference_sequence_id(chrom_idx)
                .mate_alignment_start(start1 as usize)
                .tag("MI", mi)
                .build();

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
            let mut raw_writer = bam::io::Writer::new(std::fs::File::create(&raw_path).unwrap());
            raw_writer.write_header(&header).unwrap();

            // Create consensus BAM
            let mut con_writer =
                bam::io::Writer::new(std::fs::File::create(&consensus_path).unwrap());
            con_writer.write_header(&header).unwrap();

            // Add reads with variant at chr1:10 (A->T)
            let (r1, r2) =
                create_read_pair("A1", 0, 6, 50, b"AAAATAAAAA", b"AAAAAAAAAA", "A", None);
            raw_writer.write_alignment_record(&header, &r1).unwrap();
            raw_writer.write_alignment_record(&header, &r2).unwrap();

            let (r1, r2) =
                create_read_pair("A2", 0, 6, 50, b"AAAATAAAAG", b"AAAAAAAAAA", "A", None);
            raw_writer.write_alignment_record(&header, &r1).unwrap();
            raw_writer.write_alignment_record(&header, &r2).unwrap();

            let (r1, r2) = create_read_pair("A", 0, 6, 50, b"AAAATAAAAN", b"AAAAAAAAAA", "A", None);
            con_writer.write_alignment_record(&header, &r1).unwrap();
            con_writer.write_alignment_record(&header, &r2).unwrap();

            // Add reads at chr1:20 (A->C)
            let (r1, r2) =
                create_read_pair("B1", 0, 16, 50, b"AAAACAAAAA", b"AAAAAAAAAA", "B", None);
            raw_writer.write_alignment_record(&header, &r1).unwrap();
            raw_writer.write_alignment_record(&header, &r2).unwrap();

            let (r1, r2) =
                create_read_pair("B", 0, 16, 50, b"AAAACAAAAA", b"AAAAAAAAAA", "B", None);
            con_writer.write_alignment_record(&header, &r1).unwrap();
            con_writer.write_alignment_record(&header, &r2).unwrap();

            // Reference read at chr1:20
            let (r1, r2) =
                create_read_pair("C1", 0, 17, 50, b"AAAAAAAAAA", b"AAAAAAAAAA", "C", None);
            raw_writer.write_alignment_record(&header, &r1).unwrap();
            raw_writer.write_alignment_record(&header, &r2).unwrap();

            let (r1, r2) =
                create_read_pair("C", 0, 17, 50, b"AAAAAAAAAA", b"AAAAAAAAAA", "C", None);
            con_writer.write_alignment_record(&header, &r1).unwrap();
            con_writer.write_alignment_record(&header, &r2).unwrap();

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
            raw_writer.write_alignment_record(&header, &r1).unwrap();
            raw_writer.write_alignment_record(&header, &r2).unwrap();

            let (r1, r2) =
                create_read_pair("D", 0, 25, 60, b"AAAAAAAAAA", b"AAAAAAAAAA", "D", Some("4M4D6M"));
            con_writer.write_alignment_record(&header, &r1).unwrap();
            con_writer.write_alignment_record(&header, &r2).unwrap();

            // E: no-call (N)
            let (r1, r2) =
                create_read_pair("E1", 0, 26, 60, b"AAAANAAAAA", b"AAAAAAAAAA", "E", None);
            raw_writer.write_alignment_record(&header, &r1).unwrap();
            raw_writer.write_alignment_record(&header, &r2).unwrap();

            let (r1, r2) =
                create_read_pair("E", 0, 26, 60, b"AAAANAAAAA", b"AAAAAAAAAA", "E", None);
            con_writer.write_alignment_record(&header, &r1).unwrap();
            con_writer.write_alignment_record(&header, &r2).unwrap();

            // F: variant allele (G)
            let (r1, r2) =
                create_read_pair("F1", 0, 27, 60, b"AAAGAAAAAA", b"AAAAAAAAAA", "F", None);
            raw_writer.write_alignment_record(&header, &r1).unwrap();
            raw_writer.write_alignment_record(&header, &r2).unwrap();

            let (r1, r2) =
                create_read_pair("F", 0, 27, 60, b"AAAGAAAAAA", b"AAAAAAAAAA", "F", None);
            con_writer.write_alignment_record(&header, &r1).unwrap();
            con_writer.write_alignment_record(&header, &r2).unwrap();

            // Reads at chr2:20 where both ends overlap variant
            let (r1, r2) =
                create_read_pair("H1", 1, 15, 19, b"CCCCCTCCCC", b"CTCCCCCCCC", "H", None);
            raw_writer.write_alignment_record(&header, &r1).unwrap();
            raw_writer.write_alignment_record(&header, &r2).unwrap();

            let (r1, r2) =
                create_read_pair("H", 1, 15, 19, b"CCCCCTCCCC", b"CTCCCCCCCC", "H", None);
            con_writer.write_alignment_record(&header, &r1).unwrap();
            con_writer.write_alignment_record(&header, &r2).unwrap();

            // Close writers to flush data
            drop(raw_writer);
            drop(con_writer);

            // Create BAM index files using noodles bam::fs::index
            use noodles::bam::bai;

            // Index raw BAM
            let raw_index_path = raw_path.with_extension("bam.bai");
            let raw_index = bam::fs::index(&raw_path).unwrap();
            let mut raw_index_writer =
                bai::io::Writer::new(std::fs::File::create(&raw_index_path).unwrap());
            raw_index_writer.write_index(&raw_index).unwrap();

            // Index consensus BAM
            let consensus_index_path = consensus_path.with_extension("bam.bai");
            let con_index = bam::fs::index(&consensus_path).unwrap();
            let mut con_index_writer =
                bai::io::Writer::new(std::fs::File::create(&consensus_index_path).unwrap());
            con_index_writer.write_index(&con_index).unwrap();

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
            std::fs::write(&vcf_path, vcf_content).unwrap();
            vcf_path
        }

        /// Create an empty VCF
        pub fn create_empty_vcf(dir: &TempDir) -> PathBuf {
            let vcf_path = dir.path().join("empty.vcf");
            let vcf_content = b"##fileformat=VCFv4.2\n\
##contig=<ID=chr1,length=100>\n\
##contig=<ID=chr2,length=100>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttumor\n";
            std::fs::write(&vcf_path, vcf_content).unwrap();
            vcf_path
        }

        /// Create an empty interval list
        pub fn create_empty_interval_list(dir: &TempDir) -> PathBuf {
            let interval_path = dir.path().join("empty.interval_list");
            let interval_content = b"@HD\tVN:1.5\tSO:unsorted\n\
@SQ\tSN:chr1\tLN:100\n\
@SQ\tSN:chr2\tLN:100\n";
            std::fs::write(&interval_path, interval_content).unwrap();
            interval_path
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

    // Integration tests based on Scala test suite

    #[test]
    fn test_empty_vcf_produces_empty_outputs() {
        use noodles::bam;

        let temp_dir = TempDir::new().unwrap();
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

        review.execute("test").unwrap();

        // Verify output files exist
        let con_out = output_path.with_extension("consensus.bam");
        let raw_out = output_path.with_extension("grouped.bam");
        let txt_out = output_path.with_extension("txt");

        assert!(con_out.exists());
        assert!(raw_out.exists());
        assert!(txt_out.exists());

        // Verify BAMs are empty
        let mut con_reader =
            bam::io::indexed_reader::Builder::default().build_from_path(&con_out).unwrap();
        let con_header = con_reader.read_header().unwrap();
        let mut con_record = noodles::sam::alignment::RecordBuf::default();
        assert_eq!(con_reader.read_record_buf(&con_header, &mut con_record).unwrap(), 0);

        let mut raw_reader =
            bam::io::indexed_reader::Builder::default().build_from_path(&raw_out).unwrap();
        let raw_header = raw_reader.read_header().unwrap();
        let mut raw_record = noodles::sam::alignment::RecordBuf::default();
        assert_eq!(raw_reader.read_record_buf(&raw_header, &mut raw_record).unwrap(), 0);
    }

    #[test]
    fn test_empty_interval_list_produces_empty_outputs() {
        use noodles::bam;

        let temp_dir = TempDir::new().unwrap();
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

        review.execute("test").unwrap();

        // Verify output files exist and are empty
        let con_out = output_path.with_extension("consensus.bam");
        let raw_out = output_path.with_extension("grouped.bam");

        assert!(con_out.exists());
        assert!(raw_out.exists());

        // Verify BAMs are empty
        let mut con_reader =
            bam::io::indexed_reader::Builder::default().build_from_path(&con_out).unwrap();
        let con_header = con_reader.read_header().unwrap();
        let mut con_record = noodles::sam::alignment::RecordBuf::default();
        assert_eq!(con_reader.read_record_buf(&con_header, &mut con_record).unwrap(), 0);
    }

    #[test]
    fn test_extracts_correct_reads_for_variants() {
        use noodles::bam;

        let temp_dir = TempDir::new().unwrap();
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

        review.execute("test").unwrap();

        // Verify output files exist
        let con_out = output_path.with_extension("consensus.bam");
        let raw_out = output_path.with_extension("grouped.bam");
        let txt_out = output_path.with_extension("txt");

        assert!(con_out.exists());
        assert!(raw_out.exists());
        assert!(txt_out.exists());

        // Read consensus BAM and verify read names
        let mut con_reader =
            bam::io::indexed_reader::Builder::default().build_from_path(&con_out).unwrap();
        let con_header = con_reader.read_header().unwrap();

        let mut consensus_reads = Vec::new();
        let mut con_record = noodles::sam::alignment::RecordBuf::default();
        while con_reader.read_record_buf(&con_header, &mut con_record).unwrap() > 0 {
            let name = String::from_utf8_lossy(con_record.name().unwrap().as_ref()).to_string();
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
        let mut raw_reader =
            bam::io::indexed_reader::Builder::default().build_from_path(&raw_out).unwrap();
        let raw_header = raw_reader.read_header().unwrap();

        let mut raw_reads = Vec::new();
        let mut raw_record = noodles::sam::alignment::RecordBuf::default();
        while raw_reader.read_record_buf(&raw_header, &mut raw_record).unwrap() > 0 {
            let name = String::from_utf8_lossy(raw_record.name().unwrap().as_ref()).to_string();
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

        let temp_dir = TempDir::new().unwrap();
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

        review.execute("test").unwrap();

        let txt_out = output_path.with_extension("txt");
        assert!(txt_out.exists());

        // Read TSV and verify it has content
        let file = std::fs::File::open(&txt_out).unwrap();
        let reader = std::io::BufReader::new(file);
        let lines: Vec<String> = reader.lines().map(|l| l.unwrap()).collect();

        // Should have header plus data rows
        assert!(lines.len() > 1, "TSV should have header and data rows");

        // Check header contains expected columns
        let header = &lines[0];
        assert!(header.contains("chrom"));
        assert!(header.contains("pos"));
        assert!(header.contains("ref_allele"));
        assert!(header.contains("consensus_read"));
        assert!(header.contains("consensus_call"));
    }

    #[test]
    fn test_spanning_deletions_handled_correctly() {
        use noodles::bam;

        let temp_dir = TempDir::new().unwrap();
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

        review.execute("test").unwrap();

        // Read consensus BAM - should include D (spanning deletion)
        let con_out = output_path.with_extension("consensus.bam");
        let mut con_reader =
            bam::io::indexed_reader::Builder::default().build_from_path(&con_out).unwrap();
        let con_header = con_reader.read_header().unwrap();

        let mut consensus_reads = Vec::new();
        let mut con_record = noodles::sam::alignment::RecordBuf::default();
        while con_reader.read_record_buf(&con_header, &mut con_record).unwrap() > 0 {
            let name = String::from_utf8_lossy(con_record.name().unwrap().as_ref()).to_string();
            consensus_reads.push(name);
        }

        // D should be extracted (spanning deletion at chr1:30)
        assert!(consensus_reads.contains(&"D".to_string()));
    }

    #[test]
    fn test_no_calls_handled_correctly() {
        let temp_dir = TempDir::new().unwrap();
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

        review.execute("test").unwrap();

        // Read TSV and verify N bases are included
        let txt_out = output_path.with_extension("txt");
        let content = std::fs::read_to_string(&txt_out).unwrap();

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

        review2.execute("test").unwrap();

        let txt_out2 = output_path2.with_extension("txt");
        let content2 = std::fs::read_to_string(&txt_out2).unwrap();

        // E should not be in the TSV (N bases ignored)
        let lines: Vec<&str> = content2.lines().collect();
        let e_count = lines.iter().filter(|l| l.contains("E/")).count();
        assert_eq!(e_count, 0, "Should not contain consensus read E when ignoring Ns");
    }

    #[test]
    fn test_both_ends_overlapping_variant() {
        use noodles::bam;

        let temp_dir = TempDir::new().unwrap();
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

        review.execute("test").unwrap();

        // Read consensus BAM - should include both H/1 and H/2
        let con_out = output_path.with_extension("consensus.bam");
        let mut con_reader =
            bam::io::indexed_reader::Builder::default().build_from_path(&con_out).unwrap();
        let con_header = con_reader.read_header().unwrap();

        let mut h_reads = 0;
        let mut con_record = noodles::sam::alignment::RecordBuf::default();
        while con_reader.read_record_buf(&con_header, &mut con_record).unwrap() > 0 {
            let name = String::from_utf8_lossy(con_record.name().unwrap().as_ref()).to_string();
            if name == "H" {
                h_reads += 1;
            }
        }

        // Should have both /1 and /2 reads for H
        assert_eq!(h_reads, 2, "Should extract both ends of pair H overlapping variant at chr2:20");
    }

    #[test]
    fn test_missing_fasta_index_fails() {
        let temp_dir = TempDir::new().unwrap();
        let ref_path = temp_dir.path().join("ref.fa");

        // Create FASTA without .fai
        let fasta_content = b">chr1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n";
        std::fs::write(&ref_path, fasta_content).unwrap();

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
        let temp_dir = TempDir::new().unwrap();
        let ref_path = temp_dir.path().join("ref.fa");
        let fai_path = temp_dir.path().join("ref.fa.fai");

        // Create FASTA with .fai but without .dict
        let fasta_content = b">chr1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n";
        std::fs::write(&ref_path, fasta_content).unwrap();

        let fai_content = b"chr1\t100\t6\t100\t101\n";
        std::fs::write(&fai_path, fai_content).unwrap();

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

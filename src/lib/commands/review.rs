//! Review consensus variant calls by extracting supporting reads.
//!
//! This tool extracts consensus reads containing variants and their supporting
//! raw reads to facilitate manual review of variant calls. It creates filtered
//! BAM files and a detailed TSV report.

use crate::logging::OperationTimer;
use crate::metrics::write_metrics;
use crate::reference::find_dict_path;
use crate::sam::SamTag;
use crate::umi::extract_mi_base;
use crate::validation::validate_file_exists;
use crate::variant_review::{
    BaseCounts, ConsensusVariantReviewInfo, Variant, format_insert_string, read_number_suffix,
};
use anyhow::{Context, Result, bail};
use clap::Parser;
use fgumi_raw_bam::{
    BAM_BASE_TO_ASCII, CigarKind, IndexedRawBamReader, RawRecord, find_string_tag_in_record,
    write_raw_record,
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
/// Borrows from the record so the per-read MI checks allocate nothing.
fn to_mi(record: &RawRecord) -> Result<&str> {
    match find_string_tag_in_record(record, SamTag::MI) {
        Some(mi_bytes) => Ok(extract_mi_base(std::str::from_utf8(mi_bytes)?)),
        None => {
            bail!(
                "{} did not have a value for tag MI; review needs the MI tag that fgumi group \
                 and the consensus callers write",
                String::from_utf8_lossy(record.read_name())
            )
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
the consensus reads that do not contain a reference allele at the variant sites, and the raw reads
overlapping the variant sites that contributed to those consensus reads. This will include consensus
reads that carry the alternate allele, a third allele, a no-call or a spanning deletion at the
variant site.

Reads are correlated between consensus and grouped BAMs using a molecule ID stored in an optional
attribute, `MI` by default. In order to support paired molecule IDs where two or more molecule IDs
are related (e.g. see the Paired assignment strategy in `group`) the molecule ID is truncated at
the last `/` if present (e.g. `1/A => 1` and `2 => 2`). Every consensus read that is extracted, and
every raw read that overlaps a variant site, must carry the molecule ID; a missing one is an error.

Both input BAMs must be coordinate sorted and indexed.

## Output Files

A pair of output BAMs are created:
- **<output>.consensus.bam**: Contains the relevant consensus reads from the consensus BAM
- **<output>.grouped.bam**: Contains the raw reads from the grouped BAM that overlap a variant site
  and share a molecule ID with an extracted consensus read (mates that overlap no variant site are
  not included)

Each suffix is appended to the whole `--output` prefix, including any dots it contains, so
`--output out.v1` writes `out.v1.consensus.bam`, `out.v1.grouped.bam` and `out.v1.txt`. A trailing
path separator is dropped (`out/` writes `out.consensus.bam` beside `out`), and the prefix must name
a file, not `.`, `..` or `/`. No output may be the same file as an input.

A review file **<output>.txt** is also created. The review file contains details on each variant
position along with detailed information on each consensus read that supports the variant. If the
`--sample` argument is supplied and the input is VCF, genotype information for that sample will be
retrieved. If the sample name isn't supplied and the VCF contains only a single sample then those
genotypes will be used.

The `--maf` parameter controls which variants get detailed per-read information in the output file.
Only variants with a minor allele frequency at or below this threshold will have detailed information
written.

## Known divergences from fgbio `ReviewConsensusVariants`

These are intentional, documented behavior differences; the parity-relevant ones
(SNP selection, base-count dedup, VCF REF, filters default, read-number suffix,
MAF gating) are matched exactly.

- **Overlapping interval lists are not merged (REV3-06):** fgbio applies
  `IntervalList.uniqued(false)` (sorts and merges overlapping/abutting intervals);
  fgumi emits one variant per position per interval, so overlapping intervals in
  the input produce duplicate rows. Provide a pre-merged interval list to match.
- **CIGAR `N` (ref-skip) sites (REV3-07):** an RNA read that skips over a variant
  site with an `N` operator is extracted by fgbio (its `readPosAtRefPos==0` fires
  for both `D` and `N` gaps) but not by fgumi. RNA-only; irrelevant to DNA
  consensus review.
- **Vendor-QC-fail reads (REV3-09):** fgbio's `SamLocusIterator` pileup runs with an
  emptied filter list (`setSamFilters(Collections.emptyList())`) and a mapping-quality
  cutoff of 0, so its only hardcoded exclusions are unmapped and vendor-QC-fail reads —
  secondary/supplementary reads are *not* skipped. fgumi's region queries additionally
  include vendor-QC-fail reads in counts and rows. Consensus reads essentially never set
  the QC-fail flag, so this is a no-op in practice.
- **Multi-valued `AF` (REV3-11):** fgbio parses a genotype `AF` with `toDouble`
  (which throws on a comma-separated list); fgumi takes the first value. Divergent
  only for multiallelic `AF` FORMAT fields.
- **CLI / output surface (REV3-12):** the N-ignoring flag is `--ignore-ns` (`-N`)
  rather than fgbio's `--ignore-ns-in-consensus-reads`; fgumi always writes
  `.bam.bai` sidecars for the output BAMs.
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
    ///
    /// Each suffix is appended to the whole prefix, dots included (`out.v1` →
    /// `out.v1.<suffix>`). A trailing path separator is dropped, so `out/` writes
    /// `out.<suffix>` beside `out`, not inside it. The prefix must name a file, not
    /// `.`, `..` or `/`.
    #[arg(short = 'o', long = "output")]
    pub output: PathBuf,

    /// Name of sample being reviewed (for VCF genotype extraction)
    #[arg(short = 's', long = "sample")]
    pub sample: Option<String>,

    /// Ignore N bases in consensus reads
    #[arg(short = 'N', long = "ignore-ns", value_name = "true|false", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = clap::builder::BoolishValueParser::new(), hide_possible_values = true)]
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
        validate_file_exists(&self.input, "input")?;
        validate_file_exists(&self.consensus_bam, "consensus BAM")?;
        validate_file_exists(&self.grouped_bam, "grouped BAM")?;
        validate_file_exists(&self.reference, "reference")?;

        // Refuse an `--output` prefix that names no file, or whose outputs would overwrite an
        // input, before any writer opens (`File::create` truncates the input mid-read).
        self.reject_outputs_aliasing_inputs()?;

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

        let consensus_out_path = self.output_path("consensus.bam")?;
        let grouped_out_path = self.output_path("grouped.bam")?;

        // Extract consensus reads with variants
        info!("Extracting consensus reads with variants...");
        let (mi_set, consensus_staged) =
            self.extract_consensus_reads(&variants, &consensus_out_path, command_line)?;
        info!("Found {} unique molecular identifiers", mi_set.len());

        // Extract grouped reads matching those MIs
        info!("Extracting grouped reads...");
        let grouped_staged =
            self.extract_grouped_reads(&variants, &mi_set, &grouped_out_path, command_line)?;

        // Both BAMs are complete: move them into place. An error in either extraction
        // returns above, dropping (deleting) the staged files, so a failed run never
        // leaves a truncated BAM at an output path. The renames and index writes below
        // are not atomic as a set: a failure part-way can leave a new BAM beside an
        // older run's grouped BAM or index.
        consensus_staged.persist(&consensus_out_path)?;
        grouped_staged.persist(&grouped_out_path)?;

        // Create BAM indexes for output files
        info!("Creating BAM indexes...");

        // Both paths are built with a `.bam` suffix above, so the sidecar naming
        // is unchanged here; going through the shared helper keeps every writer
        // on one convention rather than re-deriving it per call site.
        fgumi_bam_io::write_bai_sidecar(&consensus_out_path)?;
        fgumi_bam_io::write_bai_sidecar(&grouped_out_path)?;

        // Generate detailed review file
        info!("Generating detailed review file...");
        self.generate_review_file(&variants)?;

        info!("Done!");
        timer.log_completion(variants.len() as u64);
        Ok(())
    }
}

impl Review {
    /// Build an output path by appending `.{suffix}` to the full `--output` prefix.
    ///
    /// The prefix is kept intact even when it contains dots (e.g. `out.v1` →
    /// `out.v1.consensus.bam`), matching fgbio; `Path::with_extension` would instead
    /// replace the prefix's last dotted component. Trailing path separators are dropped
    /// by [`append_suffix`](crate::commands::common::append_suffix), so `out/` writes
    /// `out.consensus.bam`.
    ///
    /// # Errors
    ///
    /// Returns an error if the `--output` prefix does not name a file.
    fn output_path(&self, suffix: &str) -> Result<PathBuf> {
        crate::commands::common::append_suffix(&self.output, suffix)
    }

    /// Reject a run whose outputs (both BAMs, their `.bai` sidecars and the review `.txt`)
    /// would overwrite one of its inputs or an input BAM's index, e.g. `-o s -c
    /// s.consensus.bam`. Identity is by file (dev+inode on unix), so symlinks and hard links
    /// are caught too.
    ///
    /// # Errors
    ///
    /// Returns an error if the `--output` prefix does not name a file, or naming both paths
    /// if an output is the same file as an input.
    fn reject_outputs_aliasing_inputs(&self) -> Result<()> {
        let consensus_out = self.output_path("consensus.bam")?;
        let grouped_out = self.output_path("grouped.bam")?;
        let review_out = self.output_path("txt")?;
        let consensus_out_bai = fgumi_bam_io::bai_sidecar_path(&consensus_out);
        let grouped_out_bai = fgumi_bam_io::bai_sidecar_path(&grouped_out);
        // `read_bam_index` reads either index form, and re-reads them after the output
        // sidecars are written, so an index an output lands on is an input too.
        let consensus_bai = fgumi_bam_io::bai_sidecar_path(&self.consensus_bam);
        let consensus_alt_bai = self.consensus_bam.with_extension("bai");
        let grouped_bai = fgumi_bam_io::bai_sidecar_path(&self.grouped_bam);
        let grouped_alt_bai = self.grouped_bam.with_extension("bai");
        let inputs: [(&Path, &str); 8] = [
            (&self.input, "--input"),
            (&self.consensus_bam, "--consensus-bam"),
            (&consensus_bai, "--consensus-bam index"),
            (&consensus_alt_bai, "--consensus-bam index"),
            (&self.grouped_bam, "--grouped-bam"),
            (&grouped_bai, "--grouped-bam index"),
            (&grouped_alt_bai, "--grouped-bam index"),
            (&self.reference, "--ref"),
        ];
        let outputs: [(&Path, &str); 5] = [
            (&consensus_out, "--output consensus BAM"),
            (&consensus_out_bai, "--output consensus BAM index"),
            (&grouped_out, "--output grouped BAM"),
            (&grouped_out_bai, "--output grouped BAM index"),
            (&review_out, "--output review file"),
        ];
        crate::commands::common::reject_writes_aliasing_inputs(&inputs, &outputs)
    }

    /// Load a BAM index, preferring the samtools sidecar (`.bai` appended to the
    /// full BAM path, e.g. `foo.bam` → `foo.bam.bai`) and falling back to the
    /// extension-replaced form (`foo.bam` → `foo.bai`) that some tools emit.
    ///
    /// The lookup must derive the primary path the same way the writers do, or a
    /// BAM whose name does not end in `.bam` would be indexed at one path and
    /// searched for at another.
    fn read_bam_index(path: &Path) -> Result<noodles::bam::bai::Index> {
        use noodles::bam;

        let bam_bai = fgumi_bam_io::bai_sidecar_path(path);
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

            // Only process SNPs, matching fgbio's `.filter(_.isSNP)`: the REF is a
            // single base AND there is at least one ALT and every ALT is a single
            // base. This drops insertions/deletions (`A -> ATG`), mixed sites
            // (`A -> C,ATG`), and monomorphic records (no ALT) that fgbio excludes.
            let alternate_bases = record.alternate_bases();
            let is_snp = record.reference_bases().len() == 1
                && !alternate_bases.as_ref().is_empty()
                && alternate_bases.as_ref().iter().all(|alt| alt.len() == 1);
            if !is_snp {
                continue;
            }

            // Get chromosome name
            let chrom = record.reference_sequence_name().to_string();

            // Get position
            let pos = usize::from(
                record
                    .variant_start()
                    .ok_or_else(|| anyhow::anyhow!("Missing variant position"))?,
            ) as i32;

            // Reference base comes from the VCF REF allele (fgbio uses
            // `v.getReference.getBases()(0)`), not the FASTA — so a VCF whose REF
            // disagrees with this FASTA still classifies reads against the VCF's
            // REF. Safe to index [0] here: `is_snp` guaranteed REF length 1.
            let ref_base = record
                .reference_bases()
                .chars()
                .next()
                .ok_or_else(|| anyhow::anyhow!("VCF record has an empty REF allele"))?
                .to_ascii_uppercase();

            // Apply MAF filtering if sample specified. fgbio keeps a variant when
            // `mafFromGenotype(...).forall(_ <= maf)`, i.e. only when the MAF is
            // absent or `<= maf`. `!(maf <= self.maf)` reproduces that: it drops
            // both `maf > self.maf` and a NaN MAF (e.g. AD summing to zero, where
            // fgbio computes `1 - 0/0 = NaN` and `NaN <= maf` is false).
            if let Some(ref sname) = sample_name
                && let Some(maf) = Self::calculate_maf(&record, &header, sname)
            {
                // Keep only when `maf <= threshold`; a MAF above the threshold
                // *and* a NaN MAF (partial comparison is false) both fall through
                // to `continue`, matching fgbio's `forall(_ <= maf)`.
                let within_threshold = maf <= self.maf;
                if !within_threshold {
                    continue;
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
    /// Mirrors fgbio's `mafFromGenotype`:
    /// 1. Try to extract the `AF` (allele frequency) attribute.
    /// 2. Fall back to `1 - ad(0)/ad.sum` from `AD` (allelic depth) if `AF` is absent.
    /// 3. Return `None` if neither is available.
    ///
    /// # Arguments
    ///
    /// * `record` - VCF record
    /// * `header` - VCF header
    /// * `sample_name` - Name of sample to extract MAF for
    ///
    /// # Returns
    ///
    /// A finite MAF, or `None` when neither `AF` nor `AD` is present. Returns
    /// `NaN` when `AD` is present but sums to zero (fgbio computes `1 - 0/0`);
    /// the caller treats a `NaN` MAF as "exclude this variant".
    fn calculate_maf(
        record: &noodles::vcf::variant::RecordBuf,
        header: &noodles::vcf::Header,
        sample_name: &str,
    ) -> Option<f64> {
        use noodles::vcf::variant::record_buf::samples::sample::Value;
        use noodles::vcf::variant::record_buf::samples::sample::value::Array;

        // Find the sample index
        let sample_idx = header.sample_names().iter().position(|s| s == sample_name)?;

        let samples = record.samples();
        let keys = samples.keys();
        let sample = samples.get_index(sample_idx)?;

        // Try AF (allele frequency) first, matching fgbio's `getExtendedAttribute("AF")`.
        // Read the typed sample value directly — a `format!("{value:?}")` Debug parse
        // never matches noodles' `Float(..)` / `Array(Float(..))` rendering, so it
        // silently fails and disables MAF filtering entirely.
        if let Some(af_idx) = keys.as_ref().iter().position(|k| k == "AF")
            && let Some(Some(value)) = sample.values().get(af_idx)
            && let Some(af) = Self::value_as_f64(value)
        {
            return Some(af);
        }

        // Fall back to AD (allele depth), matching fgbio's `gt.getAD` → `1 - ad(0)/ad.sum`.
        if let Some(ad_idx) = keys.as_ref().iter().position(|k| k == "AD")
            && let Some(Some(Value::Array(Array::Integer(depths)))) = sample.values().get(ad_idx)
            && !depths.is_empty()
        {
            // fgbio: `1 - ad(0) / ad.sum.toDouble`. When AD sums to zero this
            // is `1 - 0/0 = NaN`; the caller's `!(maf <= threshold)` check
            // then excludes the variant, matching fgbio's `forall(_ <= maf)`.
            // Missing (`.`) AD entries are treated as 0 (htsjdk's AD is a
            // dense int array).
            let total: i32 = depths.iter().map(|d| d.unwrap_or(0)).sum();
            let ref_depth = depths[0].unwrap_or(0);
            return Some(1.0 - (f64::from(ref_depth) / f64::from(total)));
        }

        None
    }

    /// Extracts the leading numeric value from a typed VCF sample `Value` (used to
    /// read `AF`). Scalars convert directly; arrays take the first present element;
    /// strings parse their first comma-separated token (fgbio does `AF.toDouble`).
    fn value_as_f64(
        value: &noodles::vcf::variant::record_buf::samples::sample::Value,
    ) -> Option<f64> {
        use noodles::vcf::variant::record_buf::samples::sample::Value;
        use noodles::vcf::variant::record_buf::samples::sample::value::Array;

        match value {
            Value::Float(f) => Some(f64::from(*f)),
            Value::Integer(i) => Some(f64::from(*i)),
            Value::String(s) => s.split(',').next()?.trim().parse::<f64>().ok(),
            Value::Array(Array::Float(vs)) => vs.iter().flatten().next().map(|f| f64::from(*f)),
            Value::Array(Array::Integer(vs)) => vs.iter().flatten().next().map(|i| f64::from(*i)),
            Value::Array(Array::String(vs)) => {
                vs.iter().flatten().next()?.split(',').next()?.trim().parse::<f64>().ok()
            }
            _ => None,
        }
    }

    /// Builds the single-base query region for `variant` (fgbio's `Variant` is a
    /// one-base `Locatable`).
    fn variant_region(variant: &Variant) -> Result<noodles::core::Region> {
        let pos = noodles::core::Position::try_from(variant.pos as usize)?;
        Ok(noodles::core::Region::new(variant.chrom.as_str(), pos..=pos))
    }

    /// Writes the reads of `input` that overlap any variant and that `keep` accepts to a
    /// staged BAM beside `output`, the shared body of fgbio's
    /// `consensusIn.query(variants).filter(..)` and `groupedIn.query(variants).filter(..)`.
    ///
    /// A single multi-interval query (one point region per variant) merges the regions'
    /// BAI chunks and scans them once, so each overlapping read is visited exactly once
    /// and in coordinate order, and the output is coordinate-sorted by construction.
    /// `make_keep` receives the output header (input header plus `@HD`/`@PG`) and returns
    /// the per-read filter.
    ///
    /// The BAM is written to a temporary file in `output`'s directory and returned
    /// finished; the caller renames it to `output` with `persist`. Dropping it instead
    /// (on any error) deletes it, so no partial BAM is left at `output`.
    fn extract_variant_reads<K>(
        input: &Path,
        output: &Path,
        variants: &[Variant],
        command_line: &str,
        make_keep: impl FnOnce(&noodles::sam::Header) -> K,
    ) -> Result<tempfile::NamedTempFile>
    where
        K: FnMut(&RawRecord) -> Result<bool>,
    {
        use noodles::bam;

        let index = Self::read_bam_index(input)?;
        let mut reader = IndexedRawBamReader::from_path(input, index)?;
        let header = reader.read_header()?;

        // Synthesize @HD VN:1.6 SO:unsorted when the input lacks one (match fgbio).
        let header = crate::commands::common::ensure_hd_record(header)?;

        // Add @PG record with PP chaining
        let header = crate::commands::common::add_pg_record(header, command_line)?;

        // Stage beside the output so the final rename stays on one filesystem. The file
        // is created with the default `0o666 & !umask` mode, not `NamedTempFile`'s
        // owner-only `0o600`, so the renamed BAM has the mode a direct write would give.
        let dir = output.parent().filter(|d| !d.as_os_str().is_empty()).unwrap_or(Path::new("."));
        let staged = tempfile::Builder::new()
            .prefix(".fgumi-review-")
            .suffix(".bam.tmp")
            .make_in(dir, |path| {
                std::fs::OpenOptions::new().write(true).create_new(true).open(path)
            })
            .with_context(|| format!("creating a staging file for {}", output.display()))?;
        let mut writer = bam::io::Writer::new(staged.as_file().try_clone()?);
        writer.write_header(&header)?;

        let regions = variants.iter().map(Self::variant_region).collect::<Result<Vec<_>>>()?;
        let mut keep = make_keep(&header);
        for result in reader.query_intervals(&header, &regions)? {
            let record = result?;
            if keep(&record).with_context(|| format!("reviewing reads from {}", input.display()))? {
                write_raw_record(writer.get_mut(), &record)?;
            }
        }

        // Finish by value: `try_finish` would leave the BGZF writer live, and its `Drop`
        // would append a second EOF block. A failed final flush is an error here rather
        // than a truncated BAM.
        let file = writer.into_inner().finish()?;
        file.sync_all()?;
        Ok(staged)
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
    /// * `output` - Final `.consensus.bam` path; the BAM is staged beside it
    ///
    /// # Returns
    ///
    /// A `HashSet` of molecular identifier base strings (MI tags with strand suffix
    /// removed), and the finished BAM, staged by [`extract_variant_reads`] for the caller
    /// to `persist` to `output`.
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - Consensus BAM cannot be opened or queried
    /// - Output BAM cannot be created or written
    /// - Record processing fails, including a non-reference read with no `MI` tag
    ///
    /// [`extract_variant_reads`]: Self::extract_variant_reads
    fn extract_consensus_reads(
        &self,
        variants: &[Variant],
        output: &Path,
        command_line: &str,
    ) -> Result<(HashSet<String>, tempfile::NamedTempFile)> {
        use std::collections::HashMap;

        let mut mi_set = HashSet::new();
        let staged = Self::extract_variant_reads(
            &self.consensus_bam,
            output,
            variants,
            command_line,
            |header| {
                // Group variants by reference ID so each read is only tested against
                // variants on its own contig — `has_non_reference_base` matches on
                // position, not contig.
                let ref_id_by_name: HashMap<String, i32> = header
                    .reference_sequences()
                    .keys()
                    .enumerate()
                    .map(|(idx, name)| (String::from_utf8_lossy(name).into_owned(), idx as i32))
                    .collect();
                let mut variants_by_ref_id: HashMap<i32, Vec<&Variant>> = HashMap::new();
                for variant in variants {
                    if let Some(&ref_id) = ref_id_by_name.get(&variant.chrom) {
                        variants_by_ref_id.entry(ref_id).or_default().push(variant);
                    }
                }
                // Sort each reference's variants by position so a read only has to test
                // the variants inside its aligned span, found by binary search below,
                // rather than every variant on the contig.
                for ref_variants in variants_by_ref_id.values_mut() {
                    ref_variants.sort_by_key(|variant| variant.pos);
                }

                let mi_set = &mut mi_set;
                move |record: &RawRecord| {
                    let Some(overlapping) = variants_by_ref_id.get(&record.ref_id()) else {
                        return Ok(false); // read on a reference with no reviewed variant
                    };

                    // Only variants within the read's aligned span can be non-reference
                    // in it (`has_non_reference_base` returns false outside
                    // `[start, end]`). Restrict the scan to that position range via
                    // binary search on the sorted variants.
                    let (Some(start), Some(end)) =
                        (record.alignment_start_1based(), record.alignment_end_1based())
                    else {
                        return Ok(false); // no aligned span: not non-reference at a variant
                    };
                    #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
                    let (start, end) = (start as i32, end as i32);
                    let lo = overlapping.partition_point(|variant| variant.pos < start);
                    let hi = overlapping.partition_point(|variant| variant.pos <= end);

                    // fgbio's `nonReferenceAtAnyVariant`: keep the read if it is
                    // non-reference at any variant it overlaps.
                    let mut non_reference = false;
                    for &variant in &overlapping[lo..hi] {
                        if self.has_non_reference_base(record, variant)? {
                            non_reference = true;
                            break;
                        }
                    }
                    if non_reference {
                        // fgbio records `toMi(rec)` for every non-reference consensus
                        // read, erroring if the MI tag is absent
                        // (ReviewConsensusVariantsTest.scala:166).
                        mi_set.insert(to_mi(record)?.to_owned());
                    }
                    Ok(non_reference)
                }
            },
        )?;
        Ok((mi_set, staged))
    }

    /// Extracts the grouped raw reads that overlap a variant and came from a reviewed
    /// consensus molecule.
    ///
    /// Mirrors fgbio's `groupedIn.query(variants).filter(rec => sources.contains(toMi(rec)))`:
    /// each grouped read overlapping a variant is visited once and kept when its MI base
    /// is in `mi_set`. Mates and other reads of the molecule that overlap no variant are
    /// not written.
    ///
    /// # Arguments
    ///
    /// * `variants` - Slice of variant positions to query
    /// * `mi_set` - Set of molecular identifiers to match
    /// * `output` - Final `.grouped.bam` path; the BAM is staged beside it
    ///
    /// # Returns
    ///
    /// The finished BAM, staged by [`extract_variant_reads`] for the caller to `persist` to
    /// `output`.
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - Grouped BAM or its index cannot be opened or queried
    /// - Output BAM cannot be created or written
    /// - A grouped read overlapping a variant has no `MI` tag (fgbio's `toMi` throws)
    ///
    /// [`extract_variant_reads`]: Self::extract_variant_reads
    fn extract_grouped_reads(
        &self,
        variants: &[Variant],
        mi_set: &HashSet<String>,
        output: &Path,
        command_line: &str,
    ) -> Result<tempfile::NamedTempFile> {
        Self::extract_variant_reads(&self.grouped_bam, output, variants, command_line, |_| {
            // fgbio calls `toMi` on every variant-overlapping grouped read, so a missing
            // MI is an error here, not a skip.
            |record: &RawRecord| Ok(mi_set.contains(to_mi(record)?))
        })
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
        let review_path = self.output_path("txt")?;

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
            let region = Self::variant_region(variant)?;

            // Collect the consensus reads at this position. fgbio builds the counts and
            // rows from a `SamLocusIterator`, which skips every read flagged unmapped
            // (htsjdk `AbstractLocusIterator`), including one placed at the locus that
            // still carries a CIGAR; skip those here too.
            let mut consensus_reads: Vec<RawRecord> = Vec::new();
            for result in consensus_reader.query(&consensus_header, &region)? {
                let record = result?;
                if !record.is_unmapped() {
                    consensus_reads.push(record);
                }
            }

            // Build consensus-level base counts, deduplicated by distinct read name
            // per base to match fgbio's `BaseCounts` (BaseCounts.scala:45-46:
            // `groupBy(base).map((ch, rs) => rs.map(_.getReadName).distinct.size)`).
            // Without the dedup, a short-insert FR pair whose R1 and R2 both overlap
            // the locus double-counts the site A/C/G/T/N (which is repeated on every
            // row for the variant).
            let mut consensus_counts = BaseCounts::default();
            let mut counted_bases: HashSet<(u8, Vec<u8>)> = HashSet::new();
            for record in &consensus_reads {
                if let Some(base) = self.get_base_at_position(record, variant.pos)? {
                    let normalized = Self::normalize_base_for_variant(base, variant.ref_base);
                    if counted_bases.insert((normalized, record.read_name().to_vec())) {
                        consensus_counts.add_base(normalized);
                    }
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
                let mi_base = to_mi(&record)?.to_owned();

                // Format consensus read name with /1 or /2 suffix. fgbio's
                // `readNumberSuffix` uses `/2` only for paired second-of-pair reads;
                // unpaired reads get `/1`, so key off `paired && second`, not
                // `!first` (which mislabels unpaired reads).
                let read_name = String::from_utf8_lossy(record.read_name()).to_string();
                let is_second_of_pair = record.is_paired() && record.is_last_segment();
                let consensus_read_name =
                    format!("{}{}", read_name, read_number_suffix(is_second_of_pair));

                // fgbio sorts each variant's rows by the string `toMi(rec) +
                // readNumberSuffix(rec)` (ReviewConsensusVariants.scala:258), so e.g.
                // "10/1" sorts before "2/1". This key mirrors that (MI base, not the
                // full read name, then the /1 or /2 suffix).
                let sort_key = format!("{}{}", mi_base, read_number_suffix(is_second_of_pair));

                // Generate insert string
                let insert_string = format_insert_string(&record, &consensus_header);

                // Query grouped BAM for raw reads with matching MI
                let raw_counts = {
                    use std::collections::HashSet;

                    let mut counts = BaseCounts::default();
                    let mut seen_reads = HashSet::new();
                    // Reuse the suffix appended to `consensus_read_name` above rather
                    // than re-parsing it back out of the formatted string.
                    let expected_read_num = read_number_suffix(is_second_of_pair);

                    for result in grouped_reader.query(&grouped_header, &region)? {
                        let rec = result?;

                        // fgbio's raw counts come from a `SamLocusIterator` too, which
                        // skips unmapped reads, so a placed unmapped mate is not counted.
                        if rec.is_unmapped() {
                            continue;
                        }

                        // Source molecule id for this grouped read (fgbio's `toMi`).
                        // `extract_grouped_reads` has already rejected any MI-less read
                        // overlapping a variant, so this cannot fail on a valid run.
                        if to_mi(&rec)? != mi_base {
                            continue;
                        }

                        // Check read number matches (same `/2 iff paired && second`
                        // rule as the consensus read above, for a consistent key).
                        let is_second_of_pair = rec.is_paired() && rec.is_last_segment();
                        let read_num = read_number_suffix(is_second_of_pair);
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
                    // fgbio renders `filters.getOrElse("PASS")` — a variant with no
                    // recorded filters (every interval-list variant, and any VCF
                    // record that is PASS) shows `PASS`, not `NA`.
                    filters: variant.filters.clone().unwrap_or_else(|| "PASS".to_string()),
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

        // The shared metrics writer replaces the file atomically and writes a header-only
        // file when every variant was filtered out, as fgbio does.
        write_metrics(&review_path, &all_metrics, "review")?;
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

        /// Create a single unpaired mapped read (no PAIRED flag) with an `MI` tag.
        #[allow(clippy::cast_possible_truncation)]
        pub fn create_single_read(
            name: &str,
            chrom_idx: usize,
            start: i32,
            bases: &[u8],
            cigar: Option<&str>,
            mi: &str,
        ) -> RecordBuf {
            let cigar_str = cigar.map_or_else(|| format!("{}M", bases.len()), String::from);
            let cigar_ops = parse_cigar_to_ops(&cigar_str);
            let mut b = RawSamBuilder::new();
            b.read_name(name.as_bytes())
                .flags(0) // unpaired, mapped
                .ref_id(chrom_idx as i32)
                .pos(start - 1)
                .mapq(60)
                .cigar_ops(&cigar_ops)
                .sequence(bases)
                .qualities(&vec![45u8; bases.len()]);
            b.add_string_tag(SamTag::MI, mi.as_bytes());
            to_record_buf(b.build())
        }

        /// Create a single unpaired mapped read that carries NO `MI` tag.
        #[allow(clippy::cast_possible_truncation)]
        pub fn create_single_read_without_mi(
            name: &str,
            chrom_idx: usize,
            start: i32,
            bases: &[u8],
        ) -> RecordBuf {
            let cigar_ops = [encode_op(0, bases.len())];
            let mut b = RawSamBuilder::new();
            b.read_name(name.as_bytes())
                .flags(0) // unpaired, mapped
                .ref_id(chrom_idx as i32)
                .pos(start - 1)
                .mapq(60)
                .cigar_ops(&cigar_ops)
                .sequence(bases)
                .qualities(&vec![45u8; bases.len()]);
            to_record_buf(b.build())
        }

        /// Create a single unpaired read placed at `start` with `flags` and `cigar`,
        /// carrying `mi` when given. Covers records the other helpers cannot build: an
        /// unmapped read placed at a position that keeps a CIGAR, and a mapped read with
        /// a zero-length reference span.
        #[allow(clippy::cast_possible_truncation)]
        pub fn create_placed_read(
            name: &str,
            chrom_idx: usize,
            start: i32,
            flags: u16,
            cigar: &str,
            bases: &[u8],
            mi: Option<&str>,
        ) -> RecordBuf {
            let mut b = RawSamBuilder::new();
            b.read_name(name.as_bytes())
                .flags(flags)
                .ref_id(chrom_idx as i32)
                .pos(start - 1)
                .cigar_ops(&parse_cigar_to_ops(cigar))
                .sequence(bases)
                .qualities(&vec![45u8; bases.len()]);
            if let Some(mi) = mi {
                b.add_string_tag(SamTag::MI, mi.as_bytes());
            }
            to_record_buf(b.build())
        }

        /// Create a read pair whose R1 is mapped at `start` (all-`M` CIGAR) and whose R2
        /// is unmapped but placed at R1's position, both carrying `mi`. R2 carries the
        /// same all-`M` CIGAR as R1 when `mate_cigar` is set, and no CIGAR otherwise.
        #[allow(clippy::cast_possible_truncation)]
        pub fn create_pair_with_placed_unmapped_mate(
            name: &str,
            chrom_idx: usize,
            start: i32,
            bases: &[u8],
            mi: &str,
            mate_cigar: bool,
        ) -> (RecordBuf, RecordBuf) {
            let cigar_ops = [encode_op(0, bases.len())];
            let mut b1 = RawSamBuilder::new();
            b1.read_name(name.as_bytes())
                .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_UNMAPPED)
                .ref_id(chrom_idx as i32)
                .pos(start - 1)
                .mapq(60)
                .cigar_ops(&cigar_ops)
                .sequence(bases)
                .qualities(&vec![45u8; bases.len()])
                .mate_ref_id(chrom_idx as i32)
                .mate_pos(start - 1);
            b1.add_string_tag(SamTag::MI, mi.as_bytes());

            let mut b2 = RawSamBuilder::new();
            b2.read_name(name.as_bytes())
                .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::UNMAPPED)
                .ref_id(chrom_idx as i32)
                .pos(start - 1)
                .cigar_ops(if mate_cigar { &cigar_ops } else { &[] })
                .sequence(bases)
                .qualities(&vec![45u8; bases.len()])
                .mate_ref_id(chrom_idx as i32)
                .mate_pos(start - 1);
            b2.add_string_tag(SamTag::MI, mi.as_bytes());

            (to_record_buf(b1.build()), to_record_buf(b2.build()))
        }

        /// Create a fully unmapped read pair (no reference or position) carrying `mi`.
        pub fn create_unmapped_pair(name: &str, mi: &str) -> (RecordBuf, RecordBuf) {
            let make = |segment: u16, strand: u16| {
                let mut b = RawSamBuilder::new();
                b.read_name(name.as_bytes())
                    .flags(
                        flags::PAIRED | flags::UNMAPPED | flags::MATE_UNMAPPED | segment | strand,
                    )
                    .sequence(b"ACGTACGTAC")
                    .qualities(&[45u8; 10]);
                b.add_string_tag(SamTag::MI, mi.as_bytes());
                to_record_buf(b.build())
            };
            (
                make(flags::FIRST_SEGMENT, flags::MATE_REVERSE),
                make(flags::LAST_SEGMENT, flags::REVERSE),
            )
        }

        /// Build the `(grouped, consensus)` BAM pair from fgbio's
        /// `ReviewConsensusVariantsTest` fixture, verbatim: reads A-H at the four variant
        /// loci (chr1:10/20/30, chr2:20), including the `G` pair that overlaps chr1:30
        /// but is reference there, and the fully unmapped `X1`/`X` pairs at the end of the
        /// grouped BAM.
        pub fn create_fgbio_fixture_bams(dir: &TempDir) -> (PathBuf, PathBuf) {
            let a10 = b"AAAAAAAAAA".as_slice();
            // (name, chrom, start1, start2, cigar1, mi, bases1, bases2)
            type PairSpec = (
                &'static str,
                usize,
                i32,
                i32,
                Option<&'static str>,
                &'static str,
                &'static [u8],
                &'static [u8],
            );
            let spec: [PairSpec; 8] = [
                ("A", 0, 6, 50, None, "A", b"AAAATAAAAN", a10),
                ("B", 0, 16, 50, None, "B", b"AAAACAAAAA", a10),
                ("C", 0, 17, 50, None, "C", a10, a10),
                ("D", 0, 25, 60, Some("4M4D6M"), "D", a10, a10),
                ("E", 0, 26, 60, None, "E", b"AAAANAAAAA", a10),
                ("F", 0, 27, 60, None, "F", b"AAAGAAAAAA", a10),
                ("G", 0, 28, 60, None, "G", b"AAAAAAAGAA", a10),
                ("H", 1, 15, 19, None, "H", b"CCCCCTCCCC", b"CTCCCCCCCC"),
            ];
            let mut consensus = Vec::new();
            let mut grouped = Vec::new();
            for (name, chrom, start1, start2, cigar1, mi, bases1, bases2) in spec {
                let (r1, r2) =
                    create_read_pair(name, chrom, start1, start2, bases1, bases2, mi, cigar1);
                consensus.extend([r1, r2]);
                // Raw reads are named `<consensus>1`; `A` also has a second raw read pair
                // `A2` whose R1 differs from the consensus at the last base.
                let raw_bases1: &[u8] = if name == "A" { b"AAAATAAAAA" } else { bases1 };
                let (r1, r2) = create_read_pair(
                    &format!("{name}1"),
                    chrom,
                    start1,
                    start2,
                    raw_bases1,
                    bases2,
                    mi,
                    cigar1,
                );
                grouped.extend([r1, r2]);
            }
            let (r1, r2) = create_read_pair("A2", 0, 6, 50, b"AAAATAAAAG", a10, "A", None);
            grouped.extend([r1, r2]);
            for name in ["X1", "X"] {
                let (r1, r2) = create_unmapped_pair(name, "X");
                grouped.extend([r1, r2]);
            }

            let grouped_path = dir.path().join("raw.bam");
            let consensus_path = dir.path().join("consensus.bam");
            // `X` exists only in the grouped BAM (fgbio adds both unmapped pairs to `raw`).
            write_indexed_bam(&grouped_path, &grouped);
            write_indexed_bam(&consensus_path, &consensus);
            (grouped_path, consensus_path)
        }

        /// Collect fgbio-style read ids (`name`, plus `/1` or `/2` for paired reads),
        /// sorted, from any BAM.
        pub fn bam_read_ids(path: &std::path::Path) -> Vec<String> {
            use noodles::bam;
            let mut reader =
                bam::io::reader::Builder.build_from_path(path).expect("open bam for read ids");
            let header = reader.read_header().expect("read header");
            let mut ids: Vec<String> = reader
                .record_bufs(&header)
                .map(|result| {
                    let record = result.expect("read record");
                    let name =
                        String::from_utf8_lossy(record.name().expect("name").as_ref()).to_string();
                    let flags = record.flags();
                    let suffix = match (flags.is_segmented(), flags.is_first_segment()) {
                        (false, _) => "",
                        (true, true) => "/1",
                        (true, false) => "/2",
                    };
                    format!("{name}{suffix}")
                })
                .collect();
            ids.sort();
            ids
        }

        /// Write `records` to a coordinate-sorted, indexed BAM at `path`.
        pub fn write_indexed_bam(path: &std::path::Path, records: &[RecordBuf]) {
            use noodles::bam;
            use noodles::bam::bai;
            use noodles::sam::alignment::io::Write;

            let header = create_test_header();
            let mut writer = bam::io::Writer::new(std::fs::File::create(path).expect("create bam"));
            writer.write_header(&header).expect("write header");
            // Sort by (ref_id, pos) so the BAM is coordinate-ordered for indexing.
            let mut sorted: Vec<&RecordBuf> = records.iter().collect();
            sorted.sort_by_key(|r| {
                (
                    r.reference_sequence_id().unwrap_or(usize::MAX),
                    r.alignment_start().map(usize::from).unwrap_or(0),
                )
            });
            for rec in sorted {
                writer.write_alignment_record(&header, rec).expect("write record");
            }
            drop(writer);

            let index_path = fgumi_bam_io::bai_sidecar_path(path);
            let index = bam::fs::index(path).expect("index bam");
            let mut index_writer =
                bai::io::Writer::new(std::fs::File::create(&index_path).expect("create bai"));
            index_writer.write_index(&index).expect("write bai");
        }

        /// Parse a review `.txt` into rows of string fields (row 0 is the header).
        pub fn read_review_rows(txt_path: &std::path::Path) -> Vec<Vec<String>> {
            std::fs::read_to_string(txt_path)
                .expect("read review txt")
                .lines()
                .map(|line| line.split('\t').map(str::to_string).collect())
                .collect()
        }

        /// Expected `review` output path: `.{suffix}` appended to the full `--output` prefix.
        ///
        /// An independent oracle for `Review::output_path` — deliberately not
        /// `Path::with_extension`, which replaces the prefix's last dotted component.
        pub fn suffixed(prefix: &std::path::Path, suffix: &str) -> PathBuf {
            let mut s = prefix.as_os_str().to_owned();
            s.push(".");
            s.push(suffix);
            PathBuf::from(s)
        }

        /// Collect the read names present in an indexed output BAM (`.consensus.bam` or
        /// `.grouped.bam`).
        pub fn consensus_bam_read_names(path: &std::path::Path) -> Vec<String> {
            use noodles::bam;
            let mut reader = bam::io::indexed_reader::Builder::default()
                .build_from_path(path)
                .expect("open consensus bam");
            let header = reader.read_header().expect("read header");
            let mut names = Vec::new();
            let mut record = RecordBuf::default();
            while reader.read_record_buf(&header, &mut record).expect("read record") > 0 {
                names.push(
                    String::from_utf8_lossy(record.name().expect("name").as_ref()).to_string(),
                );
            }
            names
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
            let raw_index_path = fgumi_bam_io::bai_sidecar_path(&raw_path);
            let raw_index = bam::fs::index(&raw_path).expect("failed to index BAM file");
            let mut raw_index_writer = bai::io::Writer::new(
                std::fs::File::create(&raw_index_path).expect("failed to create file"),
            );
            raw_index_writer.write_index(&raw_index).expect("failed to write BAM index");

            // Index consensus BAM
            let consensus_index_path = fgumi_bam_io::bai_sidecar_path(&consensus_path);
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
                    std::fs::File::create(fgumi_bam_io::bai_sidecar_path(path))
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

    /// Output paths append each suffix to the full `--output` prefix, so a dotted prefix
    /// keeps every component, and a trailing path separator is dropped rather than
    /// producing hidden files inside the named directory.
    #[rstest]
    #[case::plain_prefix(
        "/path/to/output",
        "/path/to/output.consensus.bam",
        "/path/to/output.grouped.bam",
        "/path/to/output.txt"
    )]
    #[case::dotted_prefix(
        "/path/to/output.v1",
        "/path/to/output.v1.consensus.bam",
        "/path/to/output.v1.grouped.bam",
        "/path/to/output.v1.txt"
    )]
    #[case::trailing_separator(
        "/path/to/output/",
        "/path/to/output.consensus.bam",
        "/path/to/output.grouped.bam",
        "/path/to/output.txt"
    )]
    #[case::relative_trailing_separators(
        "out.v1//",
        "out.v1.consensus.bam",
        "out.v1.grouped.bam",
        "out.v1.txt"
    )]
    fn test_review_output_path_extensions(
        #[case] prefix: &str,
        #[case] expected_consensus: &str,
        #[case] expected_grouped: &str,
        #[case] expected_txt: &str,
    ) {
        let review = Review {
            input: PathBuf::from("variants.vcf"),
            consensus_bam: PathBuf::from("consensus.bam"),
            grouped_bam: PathBuf::from("grouped.bam"),
            reference: PathBuf::from("ref.fa"),
            output: PathBuf::from(prefix),
            sample: None,
            ignore_ns: false,
            maf: 0.05,
        };

        assert_eq!(review.output_path("consensus.bam").unwrap(), PathBuf::from(expected_consensus));
        assert_eq!(review.output_path("grouped.bam").unwrap(), PathBuf::from(expected_grouped));
        assert_eq!(review.output_path("txt").unwrap(), PathBuf::from(expected_txt));
    }

    /// Which `review` input a guard test places at a chosen file name.
    #[derive(Clone, Copy, PartialEq, Eq)]
    enum ReviewInput {
        Input,
        ConsensusBam,
        GroupedBam,
        Reference,
    }

    /// A `review` whose four inputs exist as empty files in `dir` (the alias guard runs before
    /// any input is parsed), with `role` named `role_name` and `prefix` as `--output`.
    fn review_with_inputs(
        dir: &TempDir,
        role: ReviewInput,
        role_name: &str,
        prefix: &str,
    ) -> Review {
        let path_for = |r: ReviewInput, default: &str| {
            let path = dir.path().join(if r == role { role_name } else { default });
            std::fs::write(&path, b"").expect("create input");
            path
        };
        Review {
            input: path_for(ReviewInput::Input, "variants.interval_list"),
            consensus_bam: path_for(ReviewInput::ConsensusBam, "in.consensus.bam"),
            grouped_bam: path_for(ReviewInput::GroupedBam, "in.grouped.bam"),
            reference: path_for(ReviewInput::Reference, "ref.fa"),
            output: dir.path().join(prefix),
            sample: None,
            ignore_ns: false,
            maf: 0.05,
        }
    }

    /// An `--output` prefix whose derived outputs (either BAM, its `.bai` sidecar, or the
    /// review `.txt`) are the same file as an input, or as an input BAM's index, is rejected
    /// before any writer opens and names both paths. The dotted-prefix cases collide only
    /// because the whole prefix is now kept (`sample.v1` → `sample.v1.consensus.bam`); the
    /// old extension-replacing naming wrote `sample.consensus.bam` instead.
    #[rstest]
    #[case::consensus_bam(
        "sample",
        ReviewInput::ConsensusBam,
        "sample.consensus.bam",
        "--output consensus BAM",
        "sample.consensus.bam",
        "--consensus-bam",
        "sample.consensus.bam"
    )]
    #[case::dotted_prefix_consensus_bam(
        "sample.v1",
        ReviewInput::ConsensusBam,
        "sample.v1.consensus.bam",
        "--output consensus BAM",
        "sample.v1.consensus.bam",
        "--consensus-bam",
        "sample.v1.consensus.bam"
    )]
    #[case::trailing_separator_consensus_bam(
        "sample/",
        ReviewInput::ConsensusBam,
        "sample.consensus.bam",
        "--output consensus BAM",
        "sample.consensus.bam",
        "--consensus-bam",
        "sample.consensus.bam"
    )]
    #[case::grouped_bam(
        "sample",
        ReviewInput::GroupedBam,
        "sample.grouped.bam",
        "--output grouped BAM",
        "sample.grouped.bam",
        "--grouped-bam",
        "sample.grouped.bam"
    )]
    #[case::dotted_prefix_grouped_bam(
        "sample.v1",
        ReviewInput::GroupedBam,
        "sample.v1.grouped.bam",
        "--output grouped BAM",
        "sample.v1.grouped.bam",
        "--grouped-bam",
        "sample.v1.grouped.bam"
    )]
    #[case::review_txt(
        "sample",
        ReviewInput::Input,
        "sample.txt",
        "--output review file",
        "sample.txt",
        "--input",
        "sample.txt"
    )]
    #[case::dotted_prefix_review_txt(
        "sample.v1",
        ReviewInput::Reference,
        "sample.v1.txt",
        "--output review file",
        "sample.v1.txt",
        "--ref",
        "sample.v1.txt"
    )]
    #[case::consensus_bam_sidecar(
        "sample",
        ReviewInput::Input,
        "sample.consensus.bam.bai",
        "--output consensus BAM index",
        "sample.consensus.bam.bai",
        "--input",
        "sample.consensus.bam.bai"
    )]
    #[case::grouped_bam_sidecar(
        "sample",
        ReviewInput::Input,
        "sample.grouped.bam.bai",
        "--output grouped BAM index",
        "sample.grouped.bam.bai",
        "--input",
        "sample.grouped.bam.bai"
    )]
    #[case::input_bam_index(
        "sample",
        ReviewInput::ConsensusBam,
        "sample.consensus.bam.bam",
        "--output consensus BAM index",
        "sample.consensus.bam.bai",
        "--consensus-bam index",
        "sample.consensus.bam.bai"
    )]
    fn test_review_rejects_an_output_that_aliases_an_input(
        #[case] prefix: &str,
        #[case] role: ReviewInput,
        #[case] role_name: &str,
        #[case] output_flag: &str,
        #[case] output_name: &str,
        #[case] input_flag: &str,
        #[case] input_name: &str,
    ) {
        let dir = TempDir::new().unwrap();
        let review = review_with_inputs(&dir, role, role_name, prefix);
        // The aliased file must exist for the guard to compare identities (an index, here).
        std::fs::write(dir.path().join(input_name), b"").unwrap();
        let files_before = std::fs::read_dir(dir.path()).unwrap().count();

        let err = review.execute("test").unwrap_err();
        assert_eq!(
            err.to_string(),
            format!(
                "{output_flag} '{}' is the same file as {input_flag} '{}'; choose a different path",
                dir.path().join(output_name).display(),
                dir.path().join(input_name).display(),
            )
        );
        assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), files_before, "nothing written");
    }

    /// A symlink onto an existing output is the same file, so it is caught by identity rather
    /// than by name.
    #[cfg(unix)]
    #[test]
    fn test_review_rejects_an_input_symlinked_to_an_output() {
        let dir = TempDir::new().unwrap();
        let review = review_with_inputs(&dir, ReviewInput::ConsensusBam, "link.bam", "sample");
        let output = dir.path().join("sample.consensus.bam");
        std::fs::write(&output, b"").unwrap();
        std::fs::remove_file(&review.consensus_bam).unwrap();
        std::os::unix::fs::symlink(&output, &review.consensus_bam).unwrap();

        let err = review.execute("test").unwrap_err();
        assert_eq!(
            err.to_string(),
            format!(
                "--output consensus BAM '{}' is the same file as --consensus-bam '{}'; choose a \
                 different path",
                output.display(),
                review.consensus_bam.display(),
            )
        );
    }

    /// An `--output` prefix that names no file is rejected before any output is written,
    /// rather than writing dot-named files such as `/.txt` or `out/...consensus.bam`.
    #[rstest]
    #[case::cur_dir(".")]
    #[case::cur_dir_trailing_separator("./")]
    #[case::parent_dir("..")]
    #[case::trailing_parent_dir("out/..")]
    #[case::root("/")]
    fn test_review_rejects_a_prefix_that_names_no_file(#[case] prefix: &str) {
        let dir = TempDir::new().unwrap();
        let mut review =
            review_with_inputs(&dir, ReviewInput::Input, "variants.interval_list", "x");
        review.output = PathBuf::from(prefix);
        let files_before = std::fs::read_dir(dir.path()).unwrap().count();

        let err = review.execute("test").unwrap_err();
        assert_eq!(
            err.to_string(),
            format!(
                "output prefix '{prefix}' does not name a file: it must end in a file name \
                 (e.g. `out` or `dir/out`), not `.`, `..` or a root"
            )
        );
        assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), files_before, "nothing written");
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
        let con_out = test_utils::suffixed(&output_path, "consensus.bam");
        let raw_out = test_utils::suffixed(&output_path, "grouped.bam");
        let txt_out = test_utils::suffixed(&output_path, "txt");

        assert!(con_out.exists());
        assert!(raw_out.exists());
        // Header-only, like every other metrics output (and fgbio), not a 0-byte file.
        let txt = std::fs::read_to_string(&txt_out).expect("failed to read review txt");
        assert_eq!(txt, format!("{}\n", ConsensusVariantReviewInfo::tsv_header()));

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
        let con_out = test_utils::suffixed(&output_path, "consensus.bam");
        let raw_out = test_utils::suffixed(&output_path, "grouped.bam");
        let txt_out = test_utils::suffixed(&output_path, "txt");

        assert!(con_out.exists());
        assert!(raw_out.exists());
        let txt = std::fs::read_to_string(&txt_out).expect("failed to read review txt");
        assert_eq!(txt, format!("{}\n", ConsensusVariantReviewInfo::tsv_header()));

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

    /// Output file names are the full `--output` prefix plus a literal suffix, so a prefix
    /// containing dots must be kept intact rather than having its last component replaced.
    ///
    /// Ported from fgbio `ReviewConsensusVariantsTest.scala:188` ("create empty BAMs when given
    /// an empty interval list as input"), which names its output prefix
    /// `makeTempFile("review_consensus.", ".out")` (i.e. `review_consensus.<random>.out`) and
    /// expects `<outBase>.consensus.bam`, `<outBase>.grouped.bam` and `<outBase>.txt`
    /// (`ReviewConsensusVariants.scala:186-187,212`). The exact set of files written next to the
    /// prefix is asserted: both BAMs, their `.bam.bai` sidecars, and the review text file. As in
    /// fgbio, both BAMs hold no records and the review file holds no rows (fgumi writes the
    /// header line, as for every other metrics output).
    ///
    /// A prefix with a trailing path separator (`out/`) names the files `out.*` next to it, not
    /// hidden `.consensus.bam` etc. inside an `out` directory, matching fgbio's `Path`
    /// normalization; `expected_stem` is the file-name stem every output must carry.
    #[rstest]
    #[case::fgbio_temp_file_prefix("review_consensus.123.out", "review_consensus.123.out")]
    #[case::single_version_component("out.v1", "out.v1")]
    #[case::many_components("a.b.c.d", "a.b.c.d")]
    #[case::trailing_separator("out/", "out")]
    fn review_output_prefix_with_dots_is_preserved(
        #[case] prefix: &str,
        #[case] expected_stem: &str,
    ) {
        let temp_dir = TempDir::new().expect("failed to create temp dir");
        let ref_path = test_utils::create_test_reference(&temp_dir);
        let (raw_path, consensus_path) = test_utils::create_test_bams(&temp_dir);
        let interval_path = test_utils::create_empty_interval_list(&temp_dir);
        let out_dir = temp_dir.path().join("outputs");
        std::fs::create_dir(&out_dir).expect("failed to create output dir");

        let review = Review {
            input: interval_path,
            consensus_bam: consensus_path,
            grouped_bam: raw_path,
            reference: ref_path,
            output: out_dir.join(prefix),
            sample: None,
            ignore_ns: false,
            maf: 0.05,
        };

        review.execute("test").expect("execute should succeed");

        let mut produced: Vec<String> = std::fs::read_dir(&out_dir)
            .expect("failed to list output dir")
            .map(|e| {
                e.expect("failed to read dir entry").file_name().to_string_lossy().into_owned()
            })
            .collect();
        produced.sort();

        let mut expected: Vec<String> =
            [".consensus.bam", ".consensus.bam.bai", ".grouped.bam", ".grouped.bam.bai", ".txt"]
                .iter()
                .map(|suffix| format!("{expected_stem}{suffix}"))
                .collect();
        expected.sort();

        assert_eq!(produced, expected);

        let consensus_out = out_dir.join(format!("{expected_stem}.consensus.bam"));
        let grouped_out = out_dir.join(format!("{expected_stem}.grouped.bam"));
        let txt_out = out_dir.join(format!("{expected_stem}.txt"));
        assert_eq!(test_utils::consensus_bam_read_names(&consensus_out), Vec::<String>::new());
        assert_eq!(test_utils::consensus_bam_read_names(&grouped_out), Vec::<String>::new());
        assert_eq!(
            std::fs::read_to_string(&txt_out).expect("failed to read review txt"),
            format!("{}\n", ConsensusVariantReviewInfo::tsv_header())
        );
    }

    /// Ports `ReviewConsensusVariantsTest.scala:227` ("extract the right reads given a set
    /// of loci"), using fgbio's fixture verbatim (including the `G` pair, which overlaps
    /// chr1:30 but is reference there, and the unmapped `X1`/`X` pairs). fgbio writes only
    /// grouped reads that *overlap a variant locus* and whose MI is in the review set, so
    /// the mates that sit off every locus (`A1/2`, `A2/2`, `B1/2`, `D1/2`, `E1/2`, `F1/2`)
    /// must not appear in `.grouped.bam`.
    #[test]
    fn test_review_grouped_bam_contains_only_reads_overlapping_variants() {
        let temp_dir = TempDir::new().expect("failed to create temp dir");
        let ref_path = test_utils::create_test_reference(&temp_dir);
        let (raw_path, consensus_path) = test_utils::create_fgbio_fixture_bams(&temp_dir);
        let vcf_path = test_utils::create_test_vcf(&temp_dir);
        let output_path = temp_dir.path().join("review_consensus_out");

        Review {
            input: vcf_path,
            consensus_bam: consensus_path,
            grouped_bam: raw_path,
            reference: ref_path,
            output: output_path.clone(),
            sample: Some("tumor".to_string()),
            ignore_ns: false,
            maf: 0.05,
        }
        .execute("test")
        .expect("execute should succeed");

        // `bam_read_ids` sorts, matching fgbio's `contain theSameElementsAs`.
        assert_eq!(
            test_utils::bam_read_ids(&test_utils::suffixed(&output_path, "grouped.bam")),
            vec!["A1/1", "A2/1", "B1/1", "D1/1", "E1/1", "F1/1", "H1/1", "H1/2"]
        );
        assert_eq!(
            test_utils::bam_read_ids(&test_utils::suffixed(&output_path, "consensus.bam")),
            vec!["A/1", "B/1", "D/1", "E/1", "F/1", "H/1", "H/2"]
        );

        // The metrics don't contain reads with spanning deletions, so `D/1` is present in
        // the consensus BAM above and absent here (fgbio test comment, `:250`).
        let rows = test_utils::read_review_rows(&test_utils::suffixed(&output_path, "txt"));
        let mut metric_reads: Vec<&str> =
            rows[1..].iter().map(|r| r[COL_CONSENSUS_READ].as_str()).collect();
        metric_reads.sort_unstable();
        assert_eq!(metric_reads, vec!["A/1", "B/1", "E/1", "F/1", "H/1", "H/2"]);
    }

    /// Builds a `Review` for a single chr1:10 `A>T` variant over a one-read consensus BAM
    /// (`CON`, MI `A`, `T` at chr1:10, spanning chr1:6-15) and a grouped BAM holding
    /// `grouped`. The output prefix is `<temp_dir>/output`.
    fn review_with_one_variant(
        temp_dir: &TempDir,
        grouped: &[noodles::sam::alignment::RecordBuf],
    ) -> Review {
        let consensus = test_utils::create_single_read("CON", 0, 6, b"AAAATAAAAA", None, "A");
        review_with_one_variant_and_consensus(temp_dir, &[consensus], grouped)
    }

    /// As [`review_with_one_variant`], with the consensus BAM holding `consensus`.
    fn review_with_one_variant_and_consensus(
        temp_dir: &TempDir,
        consensus: &[noodles::sam::alignment::RecordBuf],
        grouped: &[noodles::sam::alignment::RecordBuf],
    ) -> Review {
        let ref_path = test_utils::create_test_reference(temp_dir);
        let consensus_path = temp_dir.path().join("consensus.bam");
        let grouped_path = temp_dir.path().join("grouped.bam");
        test_utils::write_indexed_bam(&consensus_path, consensus);
        test_utils::write_indexed_bam(&grouped_path, grouped);

        let vcf = write_file(
            temp_dir,
            "one.vcf",
            "##fileformat=VCFv4.2\n\
##contig=<ID=chr1,length=100>\n\
##contig=<ID=chr2,length=100>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t10\t.\tA\tT\t.\tPASS\t.\n",
        );
        Review {
            input: vcf,
            consensus_bam: consensus_path,
            grouped_bam: grouped_path,
            reference: ref_path,
            output: temp_dir.path().join("output"),
            sample: None,
            ignore_ns: false,
            maf: 0.05,
        }
    }

    /// The grouped reads for the MI-less tests: `RAW1` (MI `A`, chr1:6-15, `T` at the
    /// chr1:10 variant) plus `NOMI`, which has no `MI` tag and starts at `nomi_start`.
    fn grouped_reads_with_mi_less_read(nomi_start: i32) -> Vec<noodles::sam::alignment::RecordBuf> {
        let bases = b"AAAATAAAAA";
        vec![
            test_utils::create_single_read("RAW1", 0, 6, bases, None, "A"),
            test_utils::create_single_read_without_mi("NOMI", 0, nomi_start, bases),
        ]
    }

    /// A grouped read that overlaps a variant locus but has no `MI` tag is an error, not a
    /// silent skip. fgbio applies `toMi` (`ReviewConsensusVariants.scala:88`) to every read
    /// from `groupedIn.query(variants)` (`ReviewConsensusVariants.scala:206`), throwing
    /// `IllegalStateException("<name> did not have a value for tag MI")`. fgbio has no test
    /// for this path; the message is copied from `toMi`.
    ///
    /// The grouped-read extraction is called directly as well as through `execute`: the
    /// review-file raw-count loop raises the same message for the same read, so an
    /// `execute`-only assertion would stay green if the extraction went back to skipping.
    #[test]
    fn test_review_errors_on_variant_overlapping_grouped_read_without_mi() {
        let temp_dir = TempDir::new().expect("failed to create temp dir");
        // `NOMI` spans chr1:8-17, which covers the variant at chr1:10.
        let review = review_with_one_variant(&temp_dir, &grouped_reads_with_mi_less_read(8));

        let expected = format!(
            "reviewing reads from {}: NOMI did not have a value for tag MI; review needs the MI \
             tag that fgumi group and the consensus callers write",
            review.grouped_bam.display()
        );
        let variants = [Variant::new("chr1".to_string(), 10, 'A')];
        let mi_set = HashSet::from(["A".to_string()]);
        let err = review
            .extract_grouped_reads(
                &variants,
                &mi_set,
                &test_utils::suffixed(&review.output, "grouped.bam"),
                "test",
            )
            .expect_err("grouped-read extraction must reject the MI-less read");
        assert_eq!(format!("{err:#}"), expected);

        let err = review.execute("test").expect_err("review must reject the MI-less read");
        assert_eq!(format!("{err:#}"), expected);
    }

    /// A failed run leaves no BAM at either output path, and no staging file behind: the
    /// grouped extraction fails on an MI-less read after `RAW1` has been written and after
    /// the consensus BAM is complete, so writing in place would leave a truncated
    /// `.grouped.bam` that still ends in a BGZF EOF block, beside a `.consensus.bam`.
    #[test]
    fn test_review_failure_leaves_no_output_bams() {
        let temp_dir = TempDir::new().expect("failed to create temp dir");
        // `NOMI` sorts after `RAW1` (chr1:6), so `RAW1` is written before the error.
        let review = review_with_one_variant(&temp_dir, &grouped_reads_with_mi_less_read(8));
        let mut before: Vec<_> = std::fs::read_dir(temp_dir.path())
            .expect("read dir")
            .map(|entry| entry.expect("entry").file_name())
            .collect();
        before.sort();

        review.execute("test").expect_err("review must reject the MI-less read");

        let mut after: Vec<_> = std::fs::read_dir(temp_dir.path())
            .expect("read dir")
            .map(|entry| entry.expect("entry").file_name())
            .collect();
        after.sort();
        assert_eq!(after, before, "a failed run must not add files to the output directory");
    }

    /// A mapped grouped read with a zero-length reference span (here CIGAR `10S`) placed at
    /// the variant is not returned by fgbio's `groupedIn.query(variants)`: htsjdk's
    /// `compareIntervalToRecord` gives it alignment end `start - 1`, which misses a point
    /// query at `start`. fgbio therefore never calls `toMi` on it, so an MI-less one is not
    /// an error and is not written.
    #[test]
    fn test_review_skips_mapped_zero_span_grouped_read_at_a_variant() {
        let temp_dir = TempDir::new().expect("failed to create temp dir");
        let grouped = [
            test_utils::create_single_read("RAW1", 0, 6, b"AAAATAAAAA", None, "A"),
            test_utils::create_placed_read("ZERO", 0, 10, 0, "10S", b"TAAAAAAAAA", None),
        ];
        let review = review_with_one_variant(&temp_dir, &grouped);
        review.execute("test").expect("a zero-span read is never inspected");
        assert_eq!(
            test_utils::bam_read_ids(&test_utils::suffixed(&review.output, "grouped.bam")),
            vec!["RAW1"]
        );
    }

    /// Each output BAM ends with exactly one 28-byte BGZF EOF block: finishing the BGZF
    /// writer with `try_finish` and then dropping it would append a second.
    #[test]
    fn test_review_output_bams_end_with_one_eof_block() {
        let temp_dir = TempDir::new().expect("failed to create temp dir");
        let grouped = [test_utils::create_single_read("RAW1", 0, 6, b"AAAATAAAAA", None, "A")];
        let review = review_with_one_variant(&temp_dir, &grouped);
        review.execute("test").expect("execute should succeed");

        let eof = fgumi_bgzf::BGZF_EOF.as_slice();
        for suffix in ["consensus.bam", "grouped.bam"] {
            let bytes =
                std::fs::read(test_utils::suffixed(&review.output, suffix)).expect("read bam");
            let (body, tail) = bytes.split_at(bytes.len() - eof.len());
            assert_eq!(tail, eof, "{suffix} must end with a BGZF EOF block");
            assert!(!body.ends_with(eof), "{suffix} must end with exactly one BGZF EOF block");
        }
    }

    /// Counterpart to the test above: fgbio only calls `toMi` on grouped reads returned by
    /// `query(variants)` (`ReviewConsensusVariants.scala:206`), so an MI-less grouped read
    /// that overlaps no variant is never inspected. It must neither error nor be written to
    /// `.grouped.bam`.
    #[test]
    fn test_review_ignores_grouped_read_without_mi_that_overlaps_no_variant() {
        let temp_dir = TempDir::new().expect("failed to create temp dir");
        // `NOMI` spans chr1:60-69, which overlaps no variant.
        let review = review_with_one_variant(&temp_dir, &grouped_reads_with_mi_less_read(60));
        review
            .execute("test")
            .expect("an MI-less grouped read off every variant locus must not error");
        assert_eq!(
            test_utils::bam_read_ids(&test_utils::suffixed(&review.output, "grouped.bam")),
            vec!["RAW1"]
        );
    }

    /// An unmapped mate placed at its mapped mate's position is a one-base read at that
    /// position for the grouped query, as in htsjdk, whose
    /// `BAMQueryMultipleIntervalsIteratorFilter.compareIntervalToRecord` sets such a
    /// read's alignment end to its alignment start. fgbio writes
    /// `groupedIn.query(variants)` reads with a reviewed MI
    /// (`ReviewConsensusVariants.scala:206`), so `RAWP/2`, placed on the chr1:10 variant,
    /// is written, and `OFFP/2`, placed at chr1:6, is not, even though its 10M CIGAR
    /// would reach chr1:10 if it were honored. fgbio has no test for this
    /// path.
    #[test]
    fn test_review_grouped_bam_includes_placed_unmapped_mate_only_at_a_variant() {
        let temp_dir = TempDir::new().expect("failed to create temp dir");
        let (on_r1, on_r2) = test_utils::create_pair_with_placed_unmapped_mate(
            "RAWP",
            0,
            10,
            b"TAAAAAAAAA",
            "A",
            false,
        );
        // `OFFP/2` keeps a 10M CIGAR, which would reach chr1:10 if it were honored.
        let (off_r1, off_r2) = test_utils::create_pair_with_placed_unmapped_mate(
            "OFFP",
            0,
            6,
            b"AAAATAAAAA",
            "A",
            true,
        );
        let review = review_with_one_variant(&temp_dir, &[off_r1, off_r2, on_r1, on_r2]);
        review.execute("test").expect("execute should succeed");
        assert_eq!(
            test_utils::bam_read_ids(&test_utils::suffixed(&review.output, "grouped.bam")),
            vec!["OFFP/1", "RAWP/1", "RAWP/2"]
        );
    }

    /// Unmapped reads placed at the variant that still carry a CIGAR are not counted and get
    /// no row. fgbio builds the review file from `SamLocusIterator`, and htsjdk's
    /// `AbstractLocusIterator` skips every read flagged unmapped. `UNMCON` (consensus, MI
    /// `B`) and `UNMRAW` (grouped, MI `A`) are placed at chr1:10 with a 10M CIGAR and a `T`
    /// there, so walking their CIGARs would add a spurious `UNMCON/1` row and count `T`
    /// twice in both the consensus and the raw columns. `UNMRAW` is still written to
    /// `.grouped.bam`, as fgbio's `groupedIn.query(variants)` returns it.
    #[test]
    fn test_review_file_skips_placed_unmapped_reads_with_a_cigar() {
        let temp_dir = TempDir::new().expect("failed to create temp dir");
        let bases = b"TAAAAAAAAA";
        let consensus = [
            test_utils::create_single_read("CON", 0, 6, b"AAAATAAAAA", None, "A"),
            test_utils::create_placed_read(
                "UNMCON",
                0,
                10,
                fgumi_raw_bam::flags::UNMAPPED,
                "10M",
                bases,
                Some("B"),
            ),
        ];
        let grouped = [
            test_utils::create_single_read("RAW1", 0, 6, b"AAAATAAAAA", None, "A"),
            test_utils::create_placed_read(
                "UNMRAW",
                0,
                10,
                fgumi_raw_bam::flags::UNMAPPED,
                "10M",
                bases,
                Some("A"),
            ),
        ];
        let review = review_with_one_variant_and_consensus(&temp_dir, &consensus, &grouped);
        review.execute("test").expect("execute should succeed");

        assert_eq!(
            test_utils::bam_read_ids(&test_utils::suffixed(&review.output, "consensus.bam")),
            vec!["CON"]
        );
        assert_eq!(
            test_utils::bam_read_ids(&test_utils::suffixed(&review.output, "grouped.bam")),
            vec!["RAW1", "UNMRAW"]
        );
        let rows = test_utils::read_review_rows(&test_utils::suffixed(&review.output, "txt"));
        let expected: Vec<String> = ["chr1", "10", "A", "NA", "PASS", "0", "0", "0", "1", "0"]
            .into_iter()
            .chain(["CON/1", "NA", "T", "45", "0", "0", "0", "1", "0"])
            .map(str::to_string)
            .collect();
        assert_eq!(rows[1..], [expected]);
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

        let txt_out = test_utils::suffixed(&output_path, "txt");
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

        let txt_out = test_utils::suffixed(&output_path, "txt");
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
        let con_out = test_utils::suffixed(&output_path, "consensus.bam");
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
        let txt_out = test_utils::suffixed(&output_path, "txt");
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
        let grouped_out = test_utils::suffixed(&output_path, "grouped.bam");
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
        let txt_out = test_utils::suffixed(&output_path, "txt");
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

        let txt_out2 = test_utils::suffixed(&output_path2, "txt");
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
        let con_out = test_utils::suffixed(&output_path, "consensus.bam");
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

    /// Helper: run `review` over the given VCF/interval input against the shared
    /// `create_test_bams` fixture and return the parsed review rows (header at [0]).
    fn run_review_over_test_bams(
        temp_dir: &TempDir,
        input: PathBuf,
        sample: Option<String>,
        maf: f64,
    ) -> (PathBuf, Vec<Vec<String>>) {
        let ref_path = test_utils::create_test_reference(temp_dir);
        let (raw_path, consensus_path) = test_utils::create_test_bams(temp_dir);
        let output_path = temp_dir.path().join("output");
        let review = Review {
            input,
            consensus_bam: consensus_path,
            grouped_bam: raw_path,
            reference: ref_path,
            output: output_path.clone(),
            sample,
            ignore_ns: false,
            maf,
        };
        review.execute("test").expect("execute should succeed");
        let rows = test_utils::read_review_rows(&test_utils::suffixed(&output_path, "txt"));
        (output_path, rows)
    }

    fn write_file(dir: &TempDir, name: &str, content: &str) -> PathBuf {
        let path = dir.path().join(name);
        std::fs::write(&path, content).expect("write fixture");
        path
    }

    // Column indices in the review TSV.
    const COL_CHROM: usize = 0;
    const COL_POS: usize = 1;
    const COL_REF: usize = 2;
    const COL_FILTERS: usize = 4;
    const COL_T: usize = 8;
    const COL_CONSENSUS_READ: usize = 10;

    /// REV3-02: site base counts dedup overlapping mates by read name. Read `H` has
    /// both ends overlapping chr2:20 calling `T`; fgbio's `BaseCounts` counts the
    /// template once (`T`=1), not once per mate (`T`=2).
    #[test]
    fn test_rev3_02_overlapping_mates_not_double_counted() {
        let temp_dir = TempDir::new().expect("temp dir");
        let vcf = test_utils::create_test_vcf(&temp_dir);
        let (_out, rows) =
            run_review_over_test_bams(&temp_dir, vcf, Some("tumor".to_string()), 0.05);

        let chr2_20_t: Vec<&str> = rows[1..]
            .iter()
            .filter(|r| r[COL_CHROM] == "chr2" && r[COL_POS] == "20")
            .map(|r| r[COL_T].as_str())
            .collect();
        assert!(!chr2_20_t.is_empty(), "expected rows at chr2:20");
        for t in chr2_20_t {
            assert_eq!(t, "1", "overlapping mates of H must count once, not twice");
        }
    }

    /// REV3-14 (subsumes REV3-08): AF- and AD-based MAF gating actually filters.
    /// The previous `format!(\"{{value:?}}\")` parse never matched noodles' typed
    /// values, so every variant was kept regardless of AF/AD.
    #[test]
    fn test_rev3_14_maf_filtering_af_and_ad() {
        let temp_dir = TempDir::new().expect("temp dir");
        let vcf = write_file(
            &temp_dir,
            "maf.vcf",
            "##fileformat=VCFv4.2\n\
##contig=<ID=chr1,length=100>\n\
##contig=<ID=chr2,length=100>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">\n\
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele Depth\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttumor\n\
chr1\t10\t.\tA\tT\t.\tPASS\t.\tGT:AF\t0/1:0.01\n\
chr1\t20\t.\tA\tC\t.\tPASS\t.\tGT:AF\t0/1:0.9\n\
chr1\t30\t.\tA\tG\t.\tPASS\t.\tGT:AD\t0/1:1,99\n\
chr2\t20\t.\tC\tT\t.\tPASS\t.\tGT:AD\t0/1:99,1\n",
        );
        let (_out, rows) =
            run_review_over_test_bams(&temp_dir, vcf, Some("tumor".to_string()), 0.05);

        let present: std::collections::HashSet<(String, String)> =
            rows[1..].iter().map(|r| (r[COL_CHROM].clone(), r[COL_POS].clone())).collect();
        assert!(present.contains(&("chr1".into(), "10".into())), "AF=0.01 must be kept");
        assert!(present.contains(&("chr2".into(), "20".into())), "AD 99,1 (maf 0.01) must be kept");
        assert!(!present.contains(&("chr1".into(), "20".into())), "AF=0.9 must be filtered out");
        assert!(
            !present.contains(&("chr1".into(), "30".into())),
            "AD 1,99 (maf 0.99) must be filtered out"
        );
    }

    /// Builds a single-sample VCF record + header carrying the given typed `AF`/`AD`
    /// sample values and runs [`Review::calculate_maf`] against it, mirroring the
    /// production lookup path (prefer `AF`, fall back to `AD`). A `None` argument omits
    /// that key entirely so the "attribute absent" branch is exercised too.
    fn calculate_maf_for(
        af: Option<noodles::vcf::variant::record_buf::samples::sample::Value>,
        ad: Option<noodles::vcf::variant::record_buf::samples::sample::Value>,
    ) -> Option<f64> {
        use noodles::vcf::Header as VcfHeader;
        use noodles::vcf::variant::RecordBuf;
        use noodles::vcf::variant::record_buf::Samples;
        use noodles::vcf::variant::record_buf::samples::Keys;

        let mut key_names: Vec<String> = Vec::new();
        let mut values = Vec::new();
        if let Some(value) = af {
            key_names.push("AF".to_string());
            values.push(Some(value));
        }
        if let Some(value) = ad {
            key_names.push("AD".to_string());
            values.push(Some(value));
        }
        let keys: Keys = key_names.into_iter().collect();
        let samples = Samples::new(keys, vec![values]);
        let record = RecordBuf::builder().set_samples(samples).build();
        let header = VcfHeader::builder().add_sample_name("s").build();
        Review::calculate_maf(&record, &header, "s")
    }

    /// Neither `AF` nor `AD` present yields `None` (no MAF gate), matching fgbio's
    /// `mafFromGenotype` returning `None`. A fixed case, so a plain test rather than a
    /// property.
    #[test]
    fn no_af_or_ad_yields_none() {
        assert_eq!(calculate_maf_for(None, None), None);
    }

    proptest::proptest! {
        /// A scalar float `AF` is read back verbatim (bit-for-bit) and takes precedence
        /// over any `AD` present, since `calculate_maf` tries `AF` first. Exercises the
        /// `Value::Float` branch of `value_as_f64` over the realistic finite AF domain.
        #[test]
        fn af_scalar_float_takes_precedence(
            af in proptest::num::f32::NORMAL | proptest::num::f32::ZERO,
            ref_depth in 0i32..10_000,
            alt_depth in 0i32..10_000,
        ) {
            use noodles::vcf::variant::record_buf::samples::sample::Value;
            use noodles::vcf::variant::record_buf::samples::sample::value::Array;
            let ad = Value::Array(Array::Integer(vec![Some(ref_depth), Some(alt_depth)]));
            let got = calculate_maf_for(Some(Value::Float(af)), Some(ad));
            proptest::prop_assert_eq!(got.map(f64::to_bits), Some(f64::from(af).to_bits()));
        }

        /// `value_as_f64` widens a scalar float / integer to exactly the same `f64`.
        #[test]
        fn value_as_f64_scalar_widens_exactly(
            f in proptest::num::f32::NORMAL | proptest::num::f32::ZERO,
            i in proptest::num::i32::ANY,
        ) {
            use noodles::vcf::variant::record_buf::samples::sample::Value;
            // Compare bit patterns: `clippy::pedantic`'s `float_cmp` forbids `==` on floats,
            // and both sides apply the identical widening so the bits match exactly.
            proptest::prop_assert_eq!(
                Review::value_as_f64(&Value::Float(f)).map(f64::to_bits),
                Some(f64::from(f).to_bits())
            );
            proptest::prop_assert_eq!(
                Review::value_as_f64(&Value::Integer(i)).map(f64::to_bits),
                Some(f64::from(i).to_bits())
            );
        }

        /// `value_as_f64` on a float array returns the first present (`Some`) element, or
        /// `None` when every element is missing or the array is empty.
        #[test]
        fn value_as_f64_array_takes_first_present(
            elements in proptest::collection::vec(
                proptest::option::of(proptest::num::f32::NORMAL | proptest::num::f32::ZERO),
                0..8,
            ),
        ) {
            use noodles::vcf::variant::record_buf::samples::sample::Value;
            use noodles::vcf::variant::record_buf::samples::sample::value::Array;
            let expected = elements.iter().flatten().next().map(|f| f64::from(*f));
            let got = Review::value_as_f64(&Value::Array(Array::Float(elements)));
            proptest::prop_assert_eq!(got.map(f64::to_bits), expected.map(f64::to_bits));
        }

        /// `value_as_f64` parses only the first comma-separated token of a string, matching
        /// fgbio's `AF.toDouble` on the leading value (the noodles `Value::String` branch).
        #[test]
        fn value_as_f64_string_parses_first_token(
            head in proptest::num::f64::NORMAL | proptest::num::f64::ZERO,
            tail in "[a-zA-Z, ]{0,8}",
        ) {
            use noodles::vcf::variant::record_buf::samples::sample::Value;
            let text = format!("{head},{tail}");
            proptest::prop_assert_eq!(
                Review::value_as_f64(&Value::String(text)).map(f64::to_bits),
                Some(head.to_bits())
            );
        }

        /// With no `AF`, MAF is `1 - ad[0] / sum(ad)`; when the depths sum to zero the
        /// result is `NaN` (fgbio's `1 - 0/0`), which the caller treats as "exclude".
        #[test]
        fn ad_only_matches_one_minus_ref_over_sum(
            depths in proptest::collection::vec(0i32..10_000, 1..6),
        ) {
            use noodles::vcf::variant::record_buf::samples::sample::Value;
            use noodles::vcf::variant::record_buf::samples::sample::value::Array;
            let total: i32 = depths.iter().sum();
            let ad = Value::Array(Array::Integer(depths.iter().map(|d| Some(*d)).collect()));
            let got = calculate_maf_for(None, Some(ad)).expect("AD present yields Some");
            if total == 0 {
                proptest::prop_assert!(got.is_nan(), "AD summing to zero must be NaN, got {}", got);
            } else {
                let expected = 1.0 - f64::from(depths[0]) / f64::from(total);
                proptest::prop_assert_eq!(got.to_bits(), expected.to_bits());
            }
        }

        /// An `AF` string that does not parse as a float falls through to the `AD` fallback
        /// rather than silently disabling MAF filtering (the pre-REV3-14 bug).
        #[test]
        fn unparseable_af_string_falls_back_to_ad(
            garbage in "[a-zA-Z]{1,8}",
            ref_depth in 1i32..10_000,
            alt_depth in 1i32..10_000,
        ) {
            use noodles::vcf::variant::record_buf::samples::sample::Value;
            use noodles::vcf::variant::record_buf::samples::sample::value::Array;
            // `[a-zA-Z]{1,8}` can spell the case-insensitive `inf`/`infinity`/`nan`
            // forms that `f64::from_str` parses as genuine floats; those take the `AF`
            // path, not the `AD` fallback this test asserts. Only truly unparseable
            // garbage exercises the fallback, so skip any sample that parses.
            proptest::prop_assume!(garbage.parse::<f64>().is_err());
            let af = Value::String(garbage);
            let ad = Value::Array(Array::Integer(vec![Some(ref_depth), Some(alt_depth)]));
            let got = calculate_maf_for(Some(af), Some(ad)).expect("AD fallback yields Some");
            let total = ref_depth + alt_depth;
            let expected = 1.0 - f64::from(ref_depth) / f64::from(total);
            proptest::prop_assert_eq!(got.to_bits(), expected.to_bits());
        }
    }

    /// REV3-01: only SNPs are reviewed. An insertion (`A -> AT`) at a covered site
    /// is dropped, so no read is extracted or reviewed for it.
    #[test]
    fn test_rev3_01_insertion_variant_dropped() {
        let temp_dir = TempDir::new().expect("temp dir");
        let vcf = write_file(
            &temp_dir,
            "ins.vcf",
            "##fileformat=VCFv4.2\n\
##contig=<ID=chr1,length=100>\n\
##contig=<ID=chr2,length=100>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t10\t.\tA\tAT\t.\tPASS\t.\n",
        );
        let (out, rows) = run_review_over_test_bams(&temp_dir, vcf, None, 0.05);
        // The insertion is dropped, so the review file is header-only (no data rows).
        assert_eq!(rows.len(), 1, "insertion-only input must yield only the header: {rows:?}");
        let names =
            test_utils::consensus_bam_read_names(&test_utils::suffixed(&out, "consensus.bam"));
        assert!(
            names.is_empty(),
            "no reads should be extracted for a dropped insertion: {names:?}"
        );
    }

    /// REV3-03: on the VCF path the reference allele comes from the VCF `REF`, not
    /// the FASTA. Here the VCF `REF` (`G`) deliberately disagrees with the FASTA
    /// (`A`) at chr1:10; the `ref` column must show `G`.
    #[test]
    fn test_rev3_03_ref_base_from_vcf_not_fasta() {
        let temp_dir = TempDir::new().expect("temp dir");
        let vcf = write_file(
            &temp_dir,
            "ref.vcf",
            "##fileformat=VCFv4.2\n\
##contig=<ID=chr1,length=100>\n\
##contig=<ID=chr2,length=100>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t10\t.\tG\tA\t.\tPASS\t.\n",
        );
        let (_out, rows) = run_review_over_test_bams(&temp_dir, vcf, None, 0.05);
        let chr1_10: Vec<&Vec<String>> =
            rows[1..].iter().filter(|r| r[COL_CHROM] == "chr1" && r[COL_POS] == "10").collect();
        assert!(!chr1_10.is_empty(), "expected rows at chr1:10");
        for r in chr1_10 {
            assert_eq!(
                r[COL_REF], "G",
                "ref column must be the VCF REF (G), not the FASTA base (A)"
            );
        }
    }

    /// REV3-04: interval-list variants have no filters, and fgbio renders the
    /// default as `PASS` (not `NA`).
    #[test]
    fn test_rev3_04_interval_list_filters_default_pass() {
        let temp_dir = TempDir::new().expect("temp dir");
        let interval = write_file(
            &temp_dir,
            "sites.interval_list",
            "@HD\tVN:1.5\tSO:unsorted\n\
@SQ\tSN:chr1\tLN:100\n\
@SQ\tSN:chr2\tLN:100\n\
chr1\t10\t10\t+\tsnp10\n",
        );
        let (_out, rows) = run_review_over_test_bams(&temp_dir, interval, None, 0.05);
        let chr1_10: Vec<&Vec<String>> =
            rows[1..].iter().filter(|r| r[COL_CHROM] == "chr1" && r[COL_POS] == "10").collect();
        assert!(!chr1_10.is_empty(), "expected rows at chr1:10");
        for r in chr1_10 {
            assert_eq!(r[COL_FILTERS], "PASS", "interval-list filters default must be PASS");
        }
    }

    /// REV3-05: a consensus read that is non-reference at two or more variant loci
    /// is written to `.consensus.bam` exactly once (not once per locus).
    #[test]
    fn test_rev3_05_multi_locus_read_written_once() {
        let temp_dir = TempDir::new().expect("temp dir");
        let ref_path = test_utils::create_test_reference(&temp_dir);

        // A single read spanning chr1:10..30 (21bp), non-ref (T) at 10 and (G) at 30.
        let mut bases = vec![b'A'; 21];
        bases[0] = b'T'; // chr1:10
        bases[20] = b'G'; // chr1:30
        let read = test_utils::create_single_read("SPAN", 0, 10, &bases, None, "S");
        let consensus_path = temp_dir.path().join("consensus.bam");
        let grouped_path = temp_dir.path().join("grouped.bam");
        test_utils::write_indexed_bam(&consensus_path, std::slice::from_ref(&read));
        test_utils::write_indexed_bam(&grouped_path, std::slice::from_ref(&read));

        let vcf = write_file(
            &temp_dir,
            "two.vcf",
            "##fileformat=VCFv4.2\n\
##contig=<ID=chr1,length=100>\n\
##contig=<ID=chr2,length=100>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t10\t.\tA\tT\t.\tPASS\t.\n\
chr1\t30\t.\tA\tG\t.\tPASS\t.\n",
        );
        let output_path = temp_dir.path().join("output");
        Review {
            input: vcf,
            consensus_bam: consensus_path,
            grouped_bam: grouped_path,
            reference: ref_path,
            output: output_path.clone(),
            sample: None,
            ignore_ns: false,
            maf: 0.05,
        }
        .execute("test")
        .expect("execute should succeed");

        let names = test_utils::consensus_bam_read_names(&test_utils::suffixed(
            &output_path,
            "consensus.bam",
        ));
        let span_count = names.iter().filter(|n| *n == "SPAN").count();
        assert_eq!(
            span_count, 1,
            "multi-locus read must be written once, got {span_count}: {names:?}"
        );
        // The grouped extraction's single multi-interval query also writes it once.
        assert_eq!(
            test_utils::bam_read_ids(&test_utils::suffixed(&output_path, "grouped.bam")),
            ["SPAN"]
        );
    }

    /// REV3-05b: `.consensus.bam` is written in coordinate order regardless of which
    /// variant locus first selected each read. fgbio's single `consensusIn.query(variants)`
    /// yields matching reads coordinate-sorted; fgumi queries per variant, so a read that
    /// is non-reference only at a *later* variant but starts *before* a read emitted at an
    /// earlier variant would land out of order without the explicit sort — leaving the
    /// indexed output BAM non-coordinate-sorted.
    #[test]
    fn test_rev3_05b_consensus_output_coordinate_sorted() {
        let temp_dir = TempDir::new().expect("temp dir");
        let ref_path = test_utils::create_test_reference(&temp_dir);

        // LATE_20: 5bp at chr1:20, non-ref (T) at chr1:20 only.
        let mut late_bases = vec![b'A'; 5];
        late_bases[0] = b'T'; // chr1:20
        let late = test_utils::create_single_read("LATE_20", 0, 20, &late_bases, None, "L");

        // EARLY_10: 21bp spanning chr1:10..30, reference (A) at chr1:20 but non-ref (G) at
        // chr1:30. It is only selected at the chr1:30 locus, yet starts before LATE_20.
        let mut early_bases = vec![b'A'; 21];
        early_bases[20] = b'G'; // chr1:30
        let early = test_utils::create_single_read("EARLY_10", 0, 10, &early_bases, None, "E");

        let consensus_path = temp_dir.path().join("consensus.bam");
        let grouped_path = temp_dir.path().join("grouped.bam");
        // `write_indexed_bam` coordinate-sorts its input, so the *input* BAM is valid; the
        // ordering under test is purely the *output* the extraction writes.
        let reads = [late.clone(), early.clone()];
        test_utils::write_indexed_bam(&consensus_path, &reads);
        test_utils::write_indexed_bam(&grouped_path, &reads);

        // Variants listed chr1:20 before chr1:30 (already coordinate order).
        let vcf = write_file(
            &temp_dir,
            "sortcheck.vcf",
            "##fileformat=VCFv4.2\n\
##contig=<ID=chr1,length=100>\n\
##contig=<ID=chr2,length=100>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t20\t.\tA\tT\t.\tPASS\t.\n\
chr1\t30\t.\tA\tG\t.\tPASS\t.\n",
        );
        let output_path = temp_dir.path().join("output");
        Review {
            input: vcf,
            consensus_bam: consensus_path,
            grouped_bam: grouped_path,
            reference: ref_path,
            output: output_path.clone(),
            sample: None,
            ignore_ns: false,
            maf: 0.05,
        }
        .execute("test")
        .expect("execute should succeed");

        // The output must be coordinate-sorted: EARLY_10 (chr1:10) before LATE_20 (chr1:20),
        // even though LATE_20 was selected first (at the earlier chr1:20 variant locus).
        let names = test_utils::consensus_bam_read_names(&test_utils::suffixed(
            &output_path,
            "consensus.bam",
        ));
        assert_eq!(
            names,
            vec!["EARLY_10".to_string(), "LATE_20".to_string()],
            "consensus output must be coordinate-sorted, got {names:?}"
        );
    }

    /// REV3-13: an unpaired consensus read gets a `/1` suffix (fgbio uses `/2` only
    /// for paired second-of-pair reads).
    #[test]
    fn test_rev3_13_unpaired_read_suffix_is_slash_one() {
        let temp_dir = TempDir::new().expect("temp dir");
        let ref_path = test_utils::create_test_reference(&temp_dir);

        let read = test_utils::create_single_read("FRAG", 0, 6, b"AAAATAAAAA", None, "S");
        let consensus_path = temp_dir.path().join("consensus.bam");
        let grouped_path = temp_dir.path().join("grouped.bam");
        test_utils::write_indexed_bam(&consensus_path, std::slice::from_ref(&read));
        test_utils::write_indexed_bam(&grouped_path, std::slice::from_ref(&read));

        let vcf = write_file(
            &temp_dir,
            "snp.vcf",
            "##fileformat=VCFv4.2\n\
##contig=<ID=chr1,length=100>\n\
##contig=<ID=chr2,length=100>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t10\t.\tA\tT\t.\tPASS\t.\n",
        );
        let output_path = temp_dir.path().join("output");
        Review {
            input: vcf,
            consensus_bam: consensus_path,
            grouped_bam: grouped_path,
            reference: ref_path,
            output: output_path.clone(),
            sample: None,
            ignore_ns: false,
            maf: 0.05,
        }
        .execute("test")
        .expect("execute should succeed");

        let rows = test_utils::read_review_rows(&test_utils::suffixed(&output_path, "txt"));
        let frag: Vec<&Vec<String>> =
            rows[1..].iter().filter(|r| r[COL_CONSENSUS_READ].starts_with("FRAG")).collect();
        assert!(!frag.is_empty(), "expected a row for FRAG");
        for r in frag {
            assert_eq!(r[COL_CONSENSUS_READ], "FRAG/1", "unpaired read must use the /1 suffix");
        }
    }
}
